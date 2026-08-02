#!/usr/bin/env Rscript

# Single-cell RNA-seq preprocessing for BMO and egress-blood samples
#
# This script performs the same ordered operations for every sample:
#   1. import the raw 10x Genomics count matrix;
#   2. calculate and apply cell-level quality-control thresholds;
#   3. normalize, scale, reduce dimensions, and cluster the data;
#   4. identify and remove predicted doublets with DoubletFinder;
#   5. estimate ambient RNA contamination with DecontX;
#   6. remove highly contaminated cells and repeat dimensional reduction;
#   7. save an auditable object after each major checkpoint.
#
# Run this file from the repository root. To process only one dataset, use:
#   Rscript analysis/01_preprocessing/01_preprocess_scrna_samples.R BMO
#   Rscript analysis/01_preprocessing/01_preprocess_scrna_samples.R EGRESS_BLOOD


suppressPackageStartupMessages({
  library(celda)
  library(DoubletFinder)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(Seurat)
})

# -----------------------------------------------------------------------------
# Step 0: define project paths and select samples
# -----------------------------------------------------------------------------

project_root <- normalizePath(
  Sys.getenv("BMO_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)

manifest_path <- file.path(project_root, "config", "sample_manifest.csv")
sample_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

command_line_args <- commandArgs(trailingOnly = TRUE)
selected_dataset <- if (length(command_line_args) == 0) "ALL" else command_line_args[[1]]

if (selected_dataset != "ALL") {
  sample_manifest <- sample_manifest[
    sample_manifest$dataset == selected_dataset,
    ,
    drop = FALSE
  ]
}

if (nrow(sample_manifest) == 0) {
  stop("No samples matched the requested dataset: ", selected_dataset)
}

processed_root <- file.path(project_root, "data", "processed", "preprocessing")
figure_root <- file.path(project_root, "results", "figures", "preprocessing")
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Step 1: define the two algorithm-specific operations used below
# -----------------------------------------------------------------------------

detect_doublets <- function(
    seurat_object,
    dimensions,
    expected_doublet_rate,
    use_sctransform = FALSE,
    adjust_for_homotypic_doublets = FALSE,
    annotation_column = NULL,
    workers = 1) {
  # Test candidate pK values and retain the value with the largest BC metric.
  parameter_sweep <- paramSweep(
    seurat_object,
    PCs = dimensions,
    sct = use_sctransform,
    num.cores = workers
  )
  sweep_summary <- summarizeSweep(parameter_sweep, GT = FALSE)
  pk_summary <- find.pK(sweep_summary)
  selected_pk <- as.numeric(as.character(
    pk_summary$pK[which.max(pk_summary$BCmetric)]
  ))

  # Convert the expected rate into an expected number of doublets.
  expected_doublets <- round(expected_doublet_rate * ncol(seurat_object))
  if (adjust_for_homotypic_doublets) {
    if (is.null(annotation_column)) {
      stop("annotation_column is required for homotypic doublet adjustment.")
    }
    homotypic_fraction <- modelHomotypic(
      seurat_object@meta.data[[annotation_column]]
    )
    expected_doublets <- round(expected_doublets * (1 - homotypic_fraction))
  }

  message(
    "DoubletFinder parameters: pK = ", selected_pk,
    "; expected doublets = ", expected_doublets
  )

  # Run DoubletFinder and rename its generated metadata columns predictably.
  seurat_object <- doubletFinder(
    seurat_object,
    PCs = dimensions,
    pN = 0.25,
    pK = selected_pk,
    nExp = expected_doublets,
    sct = use_sctransform
  )

  pann_column <- grep("^pANN", names(seurat_object@meta.data), value = TRUE)
  class_column <- grep(
    "^DF.classifications",
    names(seurat_object@meta.data),
    value = TRUE
  )

  if (length(pann_column) != 1 || length(class_column) != 1) {
    stop("DoubletFinder did not produce one unambiguous result column.")
  }

  seurat_object$pANN <- seurat_object@meta.data[[pann_column]]
  seurat_object$DF.classify <- seurat_object@meta.data[[class_column]]
  seurat_object@meta.data[[pann_column]] <- NULL
  seurat_object@meta.data[[class_column]] <- NULL
  seurat_object
}

remove_ambient_rna <- function(seurat_object, cluster_column, random_seed = 1) {
  # Extract raw counts with the Seurat v5 layer interface.
  raw_counts <- GetAssayData(
    object = seurat_object,
    assay = "RNA",
    layer = "counts"
  )
  cluster_labels <- seurat_object@meta.data[[cluster_column]]

  # DecontX excludes all-zero genes, so they are removed before estimation.
  expressed_counts <- raw_counts[rowSums(raw_counts) > 0, , drop = FALSE]
  decontx_result <- decontX(
    expressed_counts,
    z = cluster_labels,
    verbose = TRUE,
    seed = random_seed
  )

  corrected_counts <- decontx_result$decontXcounts

  # Restore any excluded genes as zero rows and recover the original gene order.
  missing_genes <- setdiff(rownames(raw_counts), rownames(corrected_counts))
  if (length(missing_genes) > 0) {
    zero_rows <- Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(raw_counts),
      dimnames = list(missing_genes, colnames(raw_counts)),
      sparse = TRUE
    )
    corrected_counts <- rbind(corrected_counts, zero_rows)
  }
  corrected_counts <- corrected_counts[rownames(raw_counts), , drop = FALSE]

  # Replace the RNA count layer and retain the estimated contamination fraction.
  seurat_object[["RNA"]] <- SetAssayData(
    seurat_object[["RNA"]],
    layer = "counts",
    new.data = as(round(corrected_counts), "dgCMatrix")
  )
  seurat_object$decontx_contamination <- decontx_result$contamination
  seurat_object
}

# -----------------------------------------------------------------------------
# Step 2: process each sample in the manifest
# -----------------------------------------------------------------------------

for (sample_index in seq_len(nrow(sample_manifest))) {
  sample_parameters <- sample_manifest[sample_index, , drop = FALSE]
  sample_id <- sample_parameters$sample_id
  dataset_id <- sample_parameters$dataset

  message("\n============================================================")
  message("Processing sample: ", sample_id, " [", dataset_id, "]")
  message("============================================================")

  raw_dataset_directory <- if (dataset_id == "BMO") "bmo" else "egress_blood"
  input_directory <- file.path(
    project_root,
    "data",
    "raw",
    raw_dataset_directory,
    sample_parameters$input_directory
  )
  sample_output_directory <- file.path(processed_root, dataset_id, sample_id)
  sample_figure_directory <- file.path(figure_root, dataset_id, sample_id)
  dir.create(sample_output_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(sample_figure_directory, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(input_directory)) {
    stop("Missing 10x Genomics input directory: ", input_directory)
  }

  dimensions <- seq_len(sample_parameters$pca_dimensions)

  # Step 2.1: import the raw matrix and construct the Seurat object.
  raw_counts <- Read10X(data.dir = input_directory)
  sample_object <- CreateSeuratObject(
    counts = raw_counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  sample_object[["percent.mt"]] <- PercentageFeatureSet(
    sample_object,
    pattern = "^MT-"
  )

  # Step 2.2: save pre-filter QC distributions and apply registered thresholds.
  qc_before_filtering <- VlnPlot(
    sample_object,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
  )
  ggsave(
    filename = file.path(sample_figure_directory, "01_qc_before_filtering.pdf"),
    plot = qc_before_filtering,
    width = 10,
    height = 4
  )

  sample_object <- subset(
    sample_object,
    subset = nFeature_RNA > sample_parameters$min_features &
      nFeature_RNA < sample_parameters$max_features &
      percent.mt < sample_parameters$max_mito_percent
  )

  # Step 2.3: run the initial Seurat dimensional-reduction workflow.
  sample_object <- NormalizeData(sample_object, verbose = FALSE)
  sample_object <- FindVariableFeatures(
    sample_object,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  sample_object <- ScaleData(
    sample_object,
    features = rownames(sample_object),
    verbose = FALSE
  )
  sample_object <- RunPCA(
    sample_object,
    features = VariableFeatures(sample_object),
    verbose = FALSE
  )
  sample_object <- FindNeighbors(
    sample_object,
    dims = dimensions,
    verbose = FALSE
  )
  sample_object <- FindClusters(
    sample_object,
    resolution = sample_parameters$clustering_resolution,
    verbose = FALSE
  )
  sample_object <- RunUMAP(sample_object, dims = dimensions, verbose = FALSE)
  sample_object <- RunTSNE(sample_object, dims = dimensions, verbose = FALSE)

  saveRDS(
    sample_object,
    file.path(sample_output_directory, "01_post_qc.rds")
  )

  # Step 2.4: detect doublets and retain predicted singlets.
  sample_object <- detect_doublets(
    seurat_object = sample_object,
    dimensions = dimensions,
    expected_doublet_rate = sample_parameters$expected_doublet_rate
  )

  doublet_plot <- DimPlot(
    sample_object,
    reduction = "umap",
    group.by = "DF.classify",
    cols = c("Doublet" = "#252525", "Singlet" = "#D7301F")
  ) + ggtitle(paste(sample_id, "doublet classification"))
  ggsave(
    filename = file.path(sample_figure_directory, "02_doublet_classification.pdf"),
    plot = doublet_plot,
    width = 6,
    height = 5
  )

  sample_object <- subset(
    sample_object,
    subset = DF.classify == "Singlet"
  )

  # Recalculate the representation because removing cells changes the manifold.
  sample_object <- NormalizeData(sample_object, verbose = FALSE)
  sample_object <- FindVariableFeatures(
    sample_object,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  sample_object <- ScaleData(
    sample_object,
    features = rownames(sample_object),
    verbose = FALSE
  )
  sample_object <- RunPCA(
    sample_object,
    features = VariableFeatures(sample_object),
    verbose = FALSE
  )
  sample_object <- FindNeighbors(
    sample_object,
    dims = dimensions,
    verbose = FALSE
  )
  sample_object <- FindClusters(
    sample_object,
    resolution = sample_parameters$clustering_resolution,
    verbose = FALSE
  )
  sample_object <- RunUMAP(sample_object, dims = dimensions, verbose = FALSE)
  sample_object <- RunTSNE(sample_object, dims = dimensions, verbose = FALSE)

  saveRDS(
    sample_object,
    file.path(sample_output_directory, "02_post_doublet_removal.rds")
  )

  # Step 2.5: estimate ambient RNA contamination and filter high-scoring cells.
  sample_object <- remove_ambient_rna(
    seurat_object = sample_object,
    cluster_column = "seurat_clusters",
    random_seed = 20260730
  )

  contamination_plot <- VlnPlot(
    sample_object,
    features = "decontx_contamination",
    group.by = "seurat_clusters",
    pt.size = 0
  ) +
    geom_hline(
      yintercept = sample_parameters$decontx_max_contamination,
      color = "#D7301F",
      linetype = "dashed"
    ) +
    ggtitle(paste(sample_id, "estimated ambient RNA contamination"))
  ggsave(
    filename = file.path(sample_figure_directory, "03_decontx_contamination.pdf"),
    plot = contamination_plot,
    width = 8,
    height = 5
  )

  keep_cells <- sample_object$decontx_contamination <
    sample_parameters$decontx_max_contamination
  sample_object <- sample_object[, keep_cells]

  # Step 2.6: recalculate the final low-dimensional representation.
  sample_object <- NormalizeData(sample_object, verbose = FALSE)
  sample_object <- FindVariableFeatures(
    sample_object,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  sample_object <- ScaleData(
    sample_object,
    features = rownames(sample_object),
    verbose = FALSE
  )
  sample_object <- RunPCA(
    sample_object,
    features = VariableFeatures(sample_object),
    verbose = FALSE
  )
  sample_object <- FindNeighbors(
    sample_object,
    dims = dimensions,
    verbose = FALSE
  )
  sample_object <- FindClusters(
    sample_object,
    resolution = sample_parameters$clustering_resolution,
    verbose = FALSE
  )
  sample_object <- RunUMAP(sample_object, dims = dimensions, verbose = FALSE)
  sample_object <- RunTSNE(sample_object, dims = dimensions, verbose = FALSE)

  final_embedding_plot <-
    DimPlot(sample_object, reduction = "umap", label = TRUE) +
    DimPlot(sample_object, reduction = "tsne", label = TRUE)
  ggsave(
    filename = file.path(sample_figure_directory, "04_final_embeddings.pdf"),
    plot = final_embedding_plot,
    width = 11,
    height = 5
  )

  saveRDS(
    sample_object,
    file.path(sample_output_directory, "03_post_decontx.rds")
  )

  # Save the software versions beside every final object.
  writeLines(
    capture.output(sessionInfo()),
    file.path(sample_output_directory, "session_info.txt")
  )

  rm(raw_counts, sample_object)
  invisible(gc())
}

message("\nPreprocessing completed for ", nrow(sample_manifest), " sample(s).")

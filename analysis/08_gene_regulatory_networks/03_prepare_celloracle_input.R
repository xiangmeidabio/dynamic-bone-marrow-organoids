#!/usr/bin/env Rscript

# Prepare erythroid data for CellOracle
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. recluster the integrated BMO erythroid object;
#   2. map clusters to three marker-supported maturation stages;
#   3. inspect stage and condition composition; and
#   4. export counts, normalized data, metadata, and embeddings to H5AD.

project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

library(Seurat)
library(scop)
library(scCustomize)

BMO_Erythroid <- readRDS(file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))

BMO_Erythroid <- FindNeighbors(BMO_Erythroid, reduction = "integrated.cca", dims = 1:10)
BMO_Erythroid <- FindClusters(BMO_Erythroid, resolution = 3)

BMO_Erythroid <- RunUMAP(BMO_Erythroid, reduction = "integrated.cca", dims = 1:10)
BMO_Erythroid <- RunTSNE(BMO_Erythroid, reduction = "integrated.cca", dims = 1:10)

BMO_Erythroid <- RenameIdents(
  BMO_Erythroid,
  `5`  = "Erythroid Stage3",
  `20` = "Erythroid Stage3",
  `4`  = "Erythroid Stage3",
  `2`  = "Erythroid Stage3",
  `28` = "Erythroid Stage3",
  `29` = "Erythroid Stage3",
  `24` = "Erythroid Stage3",
  `1`  = "Erythroid Stage3",
  `13` = "Erythroid Stage3",

  `7`  = "Erythroid Stage2",
  `0`  = "Erythroid Stage2",
  `32` = "Erythroid Stage2",
  `17` = "Erythroid Stage2",
  `18` = "Erythroid Stage2",
  `11` = "Erythroid Stage2",
  `26` = "Erythroid Stage2",
  `12` = "Erythroid Stage2",
  `21` = "Erythroid Stage2",

  `30` = "Erythroid Stage1",
  `22` = "Erythroid Stage1",
  `25` = "Erythroid Stage1",
  `8` = "Erythroid Stage1",
  `16` = "Erythroid Stage1",
  `6` = "Erythroid Stage1",
  `15` = "Erythroid Stage1",
  `9`  = "Erythroid Stage1",
  `3`  = "Erythroid Stage1",
  `31` = "Erythroid Stage1",
  `14` = "Erythroid Stage1",
  `19` = "Erythroid Stage1",
  `27` = "Erythroid Stage1",
  `10` = "Erythroid Stage1",
  `23` = "Erythroid Stage1"
)

DimPlot(BMO_Erythroid, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
DimPlot(BMO_Erythroid, reduction = "umap", label = TRUE) + ggtitle("UMAP")
DimPlot(BMO_Erythroid, reduction = "tsne", label = TRUE, group.by = "seurat_clusters") + ggtitle("t-SNE")
DimPlot(BMO_Erythroid, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + ggtitle("UMAP")


BMO_Erythroid@meta.data$celltype <- ""
BMO_Erythroid@meta.data$celltype[BMO_Erythroid@active.ident %in% c("Erythroid Stage1")] <- "Erythroid Stage1"
BMO_Erythroid@meta.data$celltype[BMO_Erythroid@active.ident %in% c("Erythroid Stage2")] <- "Erythroid Stage2"
BMO_Erythroid@meta.data$celltype[BMO_Erythroid@active.ident %in% c("Erythroid Stage3")] <- "Erythroid Stage3"

BMO_Erythroid$celltype <- factor(
  BMO_Erythroid$celltype,
  levels = c("Erythroid Stage1", "Erythroid Stage2", "Erythroid Stage3")
)

Idents(BMO_Erythroid) <- "celltype"

DimPlot(BMO_Erythroid, reduction = "umap", label = TRUE, group.by = "celltype") + ggtitle("UMAP")
DimPlot(BMO_Erythroid, reduction = "tsne", label = TRUE, group.by = "celltype") + ggtitle("TSNE")

saveRDS(BMO_Erythroid, file = file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))
BMO_Erythroid <- readRDS(file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))

markers_to_check <- c(
  "KLF1", "EPOR", "TFRC", "GATA1",
  "ALAS2", "FECH", "HMBS", "UROS", "PPOX", "ABCB6", "SLC25A37", "HBA1", "HBA2", "HBG1", "HBG2","HBB",
  "GYPA", "GYPB", "GYPC", "GYPE", "DMTN", "MPP1",
  "AHSP", "BNIP3L", "BLVRB", "FOXO3", "NFE2"
)

DotPlot(BMO_Erythroid, features = markers_to_check, group.by = "celltype") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

FeaturePlot(BMO_Erythroid, features = c("HBB"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "tsne")

FeaturePlot(BMO_Erythroid, features = c("HBB"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

FeaturePlot(BMO_Erythroid, features = c("PTPRC"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

DotPlot(BMO_Erythroid, features = markers_to_check, group.by = "celltype") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

table(BMO_Erythroid@meta.data$group)
table(BMO_Erythroid@meta.data$celltype)

library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

custom_cols <- c(
  "Erythroid Stage1" = "#FDD49E",
  "Erythroid Stage2" = "#FC8D59",
  "Erythroid Stage3" = "#B30000"
)

BMO_Erythroid$celltype <- factor(BMO_Erythroid$celltype, levels = names(custom_cols))

DimPlot(BMO_Erythroid,
        reduction = "tsne",
        group.by = "celltype",
        label = TRUE,
        cols = custom_cols,
        pt.size = 0.8) +
  theme(aspect.ratio = 1)

DimPlot(BMO_Erythroid,
        reduction = "tsne",
        split.by = "group",
        label = TRUE,
        cols = custom_cols,
        pt.size = 0.8) +
  theme(aspect.ratio = 1)

custom_cols <- c(
  "Erythroid Stage1" = "#FDD49E",
  "Erythroid Stage2" = "#FC8D59",
  "Erythroid Stage3" = "#B30000"
)

base_data <- BMO_Erythroid@meta.data %>%
  dplyr::filter(group != "") %>%
  dplyr::count(group, celltype) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(perc = n / sum(n)) %>%
  dplyr::ungroup()

base_data$celltype <- factor(base_data$celltype, levels = names(custom_cols))

data_fig1 <- base_data %>%
  dplyr::filter(group %in% c("Static.25d", "Dynamic.25d"))

data_fig1$group <- factor(data_fig1$group, levels = c("Static.25d", "Dynamic.25d"))

p1 <- ggplot(data_fig1, aes(x = group, y = perc, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = custom_cols) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(y = "Percentage of Cells", x = "Group", fill = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 11),
    legend.title = element_blank()
  )
p1

saveRDS(BMO_Erythroid, file = file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))

#seurat to h5ad-----------------
BMO_Erythroid <- readRDS(file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))
dim(BMO_Erythroid)

#BMO_Erythroid <- FindVariableFeatures(BMO_Erythroid, nfeatures = 3000)

#VariableFeatures(BMO_Erythroid)[1:3000]

#gene <- VariableFeatures(BMO_Erythroid)[1:3000]

#BMO_Erythroid <- BMO_Erythroid[gene,]
#dim(BMO_Erythroid)

DefaultAssay(BMO_Erythroid) <- 'RNA'
dim(BMO_Erythroid)

library(reticulate)
py_config()
celloracle_output_directory <- file.path(project_root, "results", "objects", "grn")
dir.create(celloracle_output_directory, recursive = TRUE, showWarnings = FALSE)
as.anndata(
  x = BMO_Erythroid,
  file_path = celloracle_output_directory,
  file_name = "BMO_Erythroid_allgene.h5ad",
  assay = "RNA",
  main_layer = "counts",
  other_layers = c("data"),
  transfer_dimreduc = TRUE,
  verbose = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_prepare_celloracle.txt")
)

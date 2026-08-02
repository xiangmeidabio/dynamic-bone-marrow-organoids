#!/usr/bin/env Rscript

# Erythroid subclustering and lineage analysis
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. extract early and late erythroid cells from egress blood;
#   2. integrate and recluster the egress-blood erythroid compartment;
#   3. assign broad haematopoietic lineages in the BMO compartment;
#   4. extract, integrate, and recluster BMO erythroid cells; and
#   5. save both erythroid objects for downstream comparisons.

project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

library(tradeSeq)
library(Seurat)
library(cowplot)
library(ggplot2)
library(Matrix)
library(dplyr)
library(RColorBrewer)
library(GenomeInfoDbData)
library(GenomeInfoDb)
library(slingshot)

Haematopoietic <- readRDS(file.path(project_root, "data", "processed", "Haematopoietic.rds"))

Haematopoietic@meta.data$group <- ""
table(Haematopoietic@meta.data$orig.ident)
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.25d.1")] <- "Dynamic.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.25d.2")] <- "Dynamic.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.25d.3")] <- "Dynamic.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.31d")] <- "Dynamic.31d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Static.25d.1")] <- "Static.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Static.25d.2")] <- "Static.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Static.25d.3")] <- "Static.25d"

table(Haematopoietic@meta.data$group)

Haematopoietic@meta.data$celltype <- ""
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("17")] <- "HSC"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("0","1","4","5","6","7","10","12","13")] <- "Erythroid"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("14")] <- "Mast"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("19")] <- "DC"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("3", "9")] <- "Macrophage"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("2", "11", "8")] <- "Monocyte"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("15")] <- "Basophil"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("16")] <- "Eosinophil"
Haematopoietic@meta.data$celltype[Haematopoietic@meta.data$seurat_clusters %in% c("18")] <- "Neutrophil"

mk_barcodes <- WhichCells(Haematopoietic, expression = PPBP > 1 | PF4 > 1)
valid_mk <- intersect(mk_barcodes, colnames(Haematopoietic))
Haematopoietic@meta.data[valid_mk, "celltype"] <- "Megakaryocyte"

table(Haematopoietic$celltype)

custom_cols <- c(
  "HSC"              = "#555555",
  "Megakaryocyte"    = "#FFD92F",
  "Erythroid"        = "#E41A1C",
  "Neutrophil"       = "#1F78B4",
  "Eosinophil"       = "#A6CEE3",
  "Basophil"         = "#33A02C",
  "Mast"             = "#B2DF8A",
  "Monocyte"         = "#FDBF6F",
  "Macrophage"       = "#B15928",
  "DC"               = "#CAB2D6"
)

Haematopoietic$celltype <- factor(Haematopoietic$celltype, levels = names(custom_cols))

DimPlot(Haematopoietic,
        reduction = "tsne",
        group.by = "celltype",
        label = TRUE,
        cols = custom_cols,
        pt.size = 0.8) +
  theme(aspect.ratio = 1)

custom_cols <- c(
  "HSC"              = "#555555",
  "Megakaryocyte"    = "#FFD92F",
  "Erythroid"        = "#E41A1C",
  "Neutrophil"       = "#1F78B4",
  "Eosinophil"       = "#A6CEE3",
  "Basophil"         = "#33A02C",
  "Mast"             = "#B2DF8A",
  "Monocyte"         = "#FDBF6F",
  "Macrophage"       = "#B15928",
  "DC"               = "#CAB2D6"
)

plot_data <- Haematopoietic@meta.data %>%
  dplyr::filter(group != "") %>%
  dplyr::count(group, celltype) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(perc = n / sum(n)) %>%
  dplyr::ungroup()

group_levels <- c("Static.25d", "Dynamic.25d", "Dynamic.31d")
plot_data$group <- factor(plot_data$group, levels = group_levels)

plot_data$celltype <- factor(plot_data$celltype, levels = names(custom_cols))

BMO_Erythroid <- subset(Haematopoietic, subset = celltype %in% c("Erythroid"))

DimPlot(BMO_Erythroid, reduction = "tsne", group.by = "celltype", label = TRUE)

BMO_Erythroid[["RNA"]] <- split(BMO_Erythroid[["RNA"]], f = BMO_Erythroid$orig.ident)

BMO_Erythroid <- NormalizeData(BMO_Erythroid)
BMO_Erythroid <- FindVariableFeatures(BMO_Erythroid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(BMO_Erythroid)
BMO_Erythroid <- ScaleData(BMO_Erythroid, features = all.genes)

BMO_Erythroid <- RunPCA(BMO_Erythroid, features = VariableFeatures(object = BMO_Erythroid))

BMO_Erythroid <- IntegrateLayers(
  object = BMO_Erythroid,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE,
  k.weight = 40)


BMO_Erythroid[["RNA"]] <- JoinLayers(BMO_Erythroid[["RNA"]])

#DimHeatmap(BMO_Erythroid, dims = 1:30, cells = 500, balanced = TRUE)
ElbowPlot(BMO_Erythroid, ndims = 30)

BMO_Erythroid <- FindNeighbors(BMO_Erythroid, reduction = "integrated.cca", dims = 1:10)
BMO_Erythroid <- FindClusters(BMO_Erythroid, resolution = 3)

BMO_Erythroid <- RunUMAP(BMO_Erythroid, reduction = "integrated.cca", dims = 1:10)
BMO_Erythroid <- RunTSNE(BMO_Erythroid, reduction = "integrated.cca", dims = 1:10)

p1 <- DimPlot(BMO_Erythroid, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
p2 <- DimPlot(BMO_Erythroid, reduction = "umap", label = TRUE) + ggtitle("UMAP")

p1+p2

DimPlot(BMO_Erythroid, reduction = "umap", label = TRUE, group.by = "orig.ident") + ggtitle("UMAP")
DimPlot(BMO_Erythroid, reduction = "tsne", label = TRUE, group.by = "orig.ident") + ggtitle("TSNE")


# Erythroid
FeaturePlot(BMO_Erythroid, features = c("TFRC","GYPA","SLC4A1","KLF1"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

FeaturePlot(BMO_Erythroid, features = c("TFRC","GYPA","SLC4A1","KLF1"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "tsne")

FeaturePlot(BMO_Erythroid, features = c("PTPRC"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

FeaturePlot(BMO_Erythroid, features = c("PTPRC"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "tsne")

library(Seurat)

DotPlot(BMO_Erythroid, features = "PTPRC") +
  scale_color_gradient(low = "grey", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(BMO_Erythroid, file = file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))
BMO_Erythroid <- readRDS(file.path(project_root, "data", "processed", "BMO_Erythroid.rds"))

FeaturePlot(BMO_Erythroid, features = c("TFRC","GYPA"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "tsne")

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_erythroid_subclustering.txt")
)

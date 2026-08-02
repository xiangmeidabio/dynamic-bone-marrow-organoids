#!/usr/bin/env Rscript

# BMO integration and cell-type annotation
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. merge the seven post-DecontX BMO sample objects;
#   2. integrate samples in a shared CCA representation;
#   3. separate stromal and haematopoietic compartments;
#   4. recluster each compartment and assign marker-supported cell labels;
#   5. add experimental groups and calculate cell-type proportions; and
#   6. save compartment-specific and combined annotated objects.

set.seed(20260730)
project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
sample_files <- c(
  "Dynamic.25d.1" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Dynamic.24d", "03_post_decontx.rds"),
  "Dynamic.25d.2" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Dynamic.25d.1", "03_post_decontx.rds"),
  "Dynamic.25d.3" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Dynamic.25d.2", "03_post_decontx.rds"),
  "Dynamic.31d" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Dynamic.31d", "03_post_decontx.rds"),
  "Static.25d.1" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Static.24d", "03_post_decontx.rds"),
  "Static.25d.2" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Static.25d.1", "03_post_decontx.rds"),
  "Static.25d.3" = file.path(project_root, "data", "processed", "preprocessing", "BMO", "Static.25d.2", "03_post_decontx.rds")
)

missing_sample_files <- sample_files[!file.exists(sample_files)]
if (length(missing_sample_files) > 0) {
  stop("Missing preprocessed BMO objects: ", paste(missing_sample_files, collapse = ", "))
}

list_scRNA <- lapply(sample_files, readRDS)

# Step 1: merge samples and calculate the initial PCA representation.
scRNA_all <- merge(list_scRNA[[1]], y = list_scRNA[-1], add.cell.ids = names(sample_files))

rm(list_scRNA)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#scRNA_all[["RNA"]] <- split(scRNA_all[["RNA"]], f = scRNA_all$orig.ident)

scRNA_all <- NormalizeData(scRNA_all)
scRNA_all <- FindVariableFeatures(scRNA_all, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(scRNA_all)
scRNA_all <- ScaleData(scRNA_all, features = all.genes)

scRNA_all <- RunPCA(scRNA_all, features = VariableFeatures(object = scRNA_all))

saveRDS(scRNA_all, file = file.path(project_root, "data", "processed", "1_scRNA_all_PCA.rds"))
scRNA_all <- readRDS(file.path(project_root, "data", "processed", "1_scRNA_all_PCA.rds"))

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Step 2: integrate samples with Seurat v5 CCA integration.
scRNA_all <- IntegrateLayers(
  object = scRNA_all,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE
)

scRNA_all[["RNA"]] <- JoinLayers(scRNA_all[["RNA"]])

saveRDS(scRNA_all, file = file.path(project_root, "data", "processed", "2_scRNA_all_CCA.rds"))
scRNA_all <- readRDS(file.path(project_root, "data", "processed", "2_scRNA_all_CCA.rds"))

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
DimHeatmap(scRNA_all, dims = 1:30, cells = 500, balanced = TRUE)
ElbowPlot(scRNA_all, ndims = 50)

scRNA_all <- FindNeighbors(scRNA_all, reduction = "integrated.cca", dims = 1:12)
scRNA_all <- FindClusters(scRNA_all, resolution = 0.8)

#scRNA_all <- RunUMAP(scRNA_all, reduction = "integrated.cca", dims = 1:12)
scRNA_all <- RunTSNE(scRNA_all, reduction = "integrated.cca", dims = 1:12)

#p1 <- DimPlot(scRNA_all, reduction = "umap", label = TRUE) + ggtitle("UMAP")
p2 <- DimPlot(scRNA_all, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")


p2
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
p1 <- DimPlot(scRNA_all, reduction = "umap", group.by = "orig.ident") + ggtitle("Samples (CCA Integrated)")
p2 <- DimPlot(scRNA_all, reduction = "umap", label = TRUE) + ggtitle("Clusters")

p1 + p2

#saveRDS(scRNA_all, file = file.path(project_root, "data", "processed", "3_Final_Integrated_BMO_CCA.rds"))

DimPlot(scRNA_all, reduction = "tsne",label = TRUE, group.by = "orig.ident")
DimPlot(scRNA_all, reduction = "tsne",label = TRUE, split.by = "orig.ident")

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Stromal(14,20,3,12,21,1,6,0,4,5,2,7,19)
FeaturePlot(scRNA_all, features = c("DCN","PLVAP","PDGFRA",
                                    "COL1A1","PDGFA"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "tsne")

# Haematopoietic(18,22,17,10,11,13,8,15,9,16)
FeaturePlot(scRNA_all, features = c("CD14","ALAS2","GP9",
                                    "TPSB2","CD7", "CD1C",
                                    "KLF1", "PTPRC"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "tsne")

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Step 3: isolate and recluster the stromal compartment.
Stromal <- scRNA_all[, scRNA_all@meta.data$seurat_clusters %in% c(14, 20, 3, 12, 21, 1, 6, 0, 4, 5, 2, 7, 19)]
DimPlot(Stromal, reduction = "tsne")

Stromal[["RNA"]] <- split(Stromal[["RNA"]], f = Stromal$orig.ident)

Stromal <- NormalizeData(Stromal)
Stromal <- FindVariableFeatures(Stromal, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Stromal)
Stromal <- ScaleData(Stromal, features = all.genes)

Stromal <- RunPCA(Stromal, features = VariableFeatures(object = Stromal))

Stromal <- IntegrateLayers(
  object = Stromal,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE
)

Stromal[["RNA"]] <- JoinLayers(Stromal[["RNA"]])

ElbowPlot(Stromal, ndims = 50)

Stromal <- FindNeighbors(Stromal, reduction = "integrated.cca", dims = 1:15)
Stromal <- FindClusters(Stromal, resolution = 0.8)

Stromal <- RunTSNE(Stromal, reduction = "integrated.cca", dims = 1:15)
Stromal <- RunUMAP(Stromal, reduction = "integrated.cca", dims = 1:15)

DimPlot(Stromal, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
DimPlot(Stromal, reduction = "umap", label = TRUE) + ggtitle("UMAP")


saveRDS(Stromal, file = file.path(project_root, "data", "processed", "Stromal.rds"))
Stromal <- readRDS(file.path(project_root, "data", "processed", "Stromal.rds"))

DimPlot(Stromal, reduction = "umap", label = TRUE) + ggtitle("UMAP")
DimPlot(Stromal, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")

DimPlot(Stromal, reduction = "tsne", label = TRUE,group.by = "seurat_clusters") + ggtitle("t-SNE")

all_stromal_markers <- c(
  "CDH5", "CD34", "PECAM1", "ICAM2", "EMCN", "FLT1", "KITLG",   # Endothelial
  
  "PDGFRA", "PDGFRB", "NGFR", "ENG", "COL3A1", "NES",
  "CXCL12", "LEPR", "VCAM1", "THY1",                         # MSC
  
  "AEBP1", "THY1", "RUNX2",                                    # Osteochondral progenitor
  
  "CSPG4", "HAPLN1", "SOX9", "COL2A1",                         # Chondrocyte
  
  "NDNF", "LIMCH1", "TNC",                                     # Osteoblast progenitor
  
  "DCN", "CPE", "ALPL", "SPP1",                                # Osteoblast
  
  "LPL", "APOE", "FABP4",                                      # Adipocyte progenitor
  
  "ANGPT1", "RGS5", "MCAM", "CSPG4", "PDGFRB",
  "ACTA2", "TAGLN", "DES"                                     # Pericyte
)

#keep only unique ones
all_stromal_markers <- unique(all_stromal_markers)


#plot custom dot plot
stroma_ipsc_dot <- DotPlot(
  Stromal,
  assay = "RNA",
  features = all_stromal_markers,
  cols = c("white", "#16a085", "#079992"),
  col.min = 0,
  #col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  #idents = NULL,
  #group.by = NULL,
  #split.by = NULL,
  #cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)  +
  #theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size = 10, angle = 90, vjust = 0.5, hjust=1, face = "italic"),
        #axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2)
  ) +
  coord_fixed()

stroma_ipsc_dot

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
Haematopoietic = scRNA_all[,scRNA_all@meta.data$seurat_clusters %in% c(18,22,17,10,11,13,8,15,9,16)]
DimPlot(Haematopoietic, reduction = "tsne")

Haematopoietic[["RNA"]] <- split(Haematopoietic[["RNA"]], f = Haematopoietic$orig.ident)

Haematopoietic <- NormalizeData(Haematopoietic)
Haematopoietic <- FindVariableFeatures(Haematopoietic, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Haematopoietic)
Haematopoietic <- ScaleData(Haematopoietic, features = all.genes)

Haematopoietic <- RunPCA(Haematopoietic, features = VariableFeatures(object = Haematopoietic))

Haematopoietic <- IntegrateLayers(
    object = Haematopoietic,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = TRUE
)

Haematopoietic[["RNA"]] <- JoinLayers(Haematopoietic[["RNA"]])

ElbowPlot(Haematopoietic, ndims = 50)

Haematopoietic <- FindNeighbors(Haematopoietic, reduction = "integrated.cca", dims = 1:10)
Haematopoietic <- FindClusters(Haematopoietic, resolution = 0.8)

Haematopoietic <- RunTSNE(Haematopoietic, reduction = "integrated.cca", dims = 1:10)
Haematopoietic <- RunUMAP(Haematopoietic, reduction = "integrated.cca", dims = 1:10)

DimPlot(Haematopoietic, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
DimPlot(Haematopoietic, reduction = "umap", label = TRUE) + ggtitle("UMAP")

#rm(scRNA_all)
#gc()

hematopoietic_markers <- c(
  "CD34", "ITGA6", "ATXN1",          # HSC
  "PPBP", "PF4", "GP9",              # Megakaryocyte
  "TFRC", "GYPA", "SLC4A1", "KLF1", "ALAS2",   # Erythroid
  "MPO", "ELANE", "CSF3R",          # Neutrophil
  "EPX",                              # Eosinophil
  "FCER1A",                           # Basophil
  "TPSAB1", "KIT", "TPSB2",          # Mast
  "S100A8", "S100A9", "CSF3R",       # Monocyte
  "CD14", "CD68", "CD163"            # Macrophage
)


haem_ipsc_dot <- DotPlot(
  Haematopoietic,
  #assay = "RNA",
  features = all_haem_markers,
  cols = c("white", "#eb0616", "#eb2f06"),
  col.min = 0,
  #col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  #idents = NULL,
  #group.by = NULL,
  #split.by = NULL,
  #cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)  +
  #theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size = 10, angle = 90, vjust = 0.5, hjust=1, face = "italic"),
        #axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2)
  ) +
  coord_fixed()

haem_ipsc_dot

Haematopoietic <- RenameIdents(Haematopoietic,
                    `17`="HSC",

                    `13`="Erythroid", `12`="Erythroid", `10`="Erythroid",
                    `7`="Erythroid", `6`="Erythroid", `5`="Erythroid", `4`="Erythroid",
                    `1`="Erythroid", `0`="Erythroid",

                    `15`="Megakaryocyte",

                    `16`="Basophil/Eosinophil",

                    `18` ="Neutrophil",

                    `14` = "Mast",

                    `12` ="Monocyte", `8` ="Monocyte", `2` ="Monocyte", `11` ="Monocyte",

                    `3` = "Macrophage", `9` = "Macrophage",

                    `19` = "DC")

Haematopoietic$celltype <- Idents(Haematopoietic)
head(Haematopoietic@meta.data)

DimPlot(Haematopoietic, label = TRUE)

saveRDS(Haematopoietic, file = file.path(project_root, "data", "processed", "Haematopoietic.rds"))
Haematopoietic <- readRDS(file.path(project_root, "data", "processed", "Haematopoietic.rds"))

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)


man_haem_bright <- c(
  "HSC"                = "#818181",
  "Erythroid"          = "#E41A1C",
  "Megakaryocyte"      = "#FF7F00",
  "Mast"               = "#FFFF33",
  "Basophil/Eosinophil"= "#33A02C",
  "Neutrophil"         = "#F781BF",
  "Monocyte"           = "#FDBF6F",
  "Macrophage"         = "#A6CEE3",
  "DC"                 = "#6A3D9A"
)

centers <- plot_data %>%
  group_by(cell_type) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

haem_perfect_final <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(fill = cell_type),
    shape = 21,
    color = "black",
    size = 1.2,
    stroke = 0.05,
    alpha = 0.9
  ) +
  scale_fill_manual(values = man_haem_bright) +
  geom_text_repel(
    data = centers,
    aes(label = cell_type),
    size = 3.5,
    fontface = "bold",
    box.padding = 0.05,
    segment.color = "grey50"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(arrow = arrow(length = unit(0.1, "cm"), type = "closed")),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(hjust = 0, face = "italic"),
    legend.position = "none",
    aspect.ratio = 1
  ) +
  labs(x = "UMAP1", y = "UMAP2")

print(haem_perfect_final)
# ggsave("Final_Bright_Haem_UMAP.pdf", haem_perfect_final, width = 6, height = 6)

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
Haematopoietic@meta.data$group <- ""
table(Haematopoietic@meta.data$orig.ident)
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.25d.1")] <- "Dynamic.24d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.25d.2")] <- "Dynamic.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.25d.3")] <- "Dynamic.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Dynamic.31d")] <- "Dynamic.31d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Static.25d.1")] <- "Static.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Static.25d.2")] <- "Static.25d"
Haematopoietic@meta.data$group[Haematopoietic@meta.data$orig.ident %in% c("Static.25d.3")] <- "Static.25d"

table(Haematopoietic@meta.data$group)
# Figure 1C - Bar chart showing the counts of each lineage assayed

## Figure 1C (left)
cal1_cols <- c("#E0B0FF", "#A7C7E7", "#AFE1AF", "#BDB5D5", "#FFB6C1", "#F28C28", "#DD3F4E")
cal1_cols <- c("#818181","#E41A1C","#6A3D9A","#fb8d62","#A6CEE3",
               "#9E0142","#B2DF8A","#FB9A99",
               "#CAB2D6","#33A02C","#FDBF6F")

cal1_cols <- c("#818181","#E41A1C","#FF7F00","#FFFF33","#F781BF",
               "#33A02C","#A6CEE3","#6A3D9A",
               "#FDBF6F")

cal1_cols <- c("#818181","#E41A1C","#FF7F00","#FFFF33","#F781BF",
               "#33A02C","#A6CEE3","#6A3D9A",
               "#FDBF6F")

p2 <- Haematopoietic@meta.data %>%
  ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(celltype)), fill = celltype)) +
  geom_bar(stat = 'count') +
  labs(x = 'Cell count', y = NULL) +
  scale_fill_manual(name = "Cell lineage",
                    values = cal1_cols,
                    #labels = c("HSPC", "Myeloid", "Lymphoid", "Meg/E", "Mesenchymal", "Endothelial", "Muscle")
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = 14, color = 'black'))

cell_counts <- as.data.frame(table(Haematopoietic@meta.data[["celltype"]],
                                   Haematopoietic@meta.data[["group"]]))

## Figure 1C (right)
p3 <- ggplot(data = cell_counts, aes(x = forcats::fct_rev(Var2),y = Freq, fill = Var1)) +
  geom_bar(position="fill",stat="identity") +
  coord_flip() +
  labs(x = NULL, y = 'Cell Lineage Frequency') +
  #scale_x_discrete(labels = c("H41", "H39", "H38", "H36", "H35", "H34", "H33", "H32", "H24", "H23", "H21", "H14")) +
  scale_fill_manual(name = "Cell lineage",
                    values = cal1_cols) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(size = 14, color = 'black'))

library(tidyverse)
library(patchwork)
library(RColorBrewer)

fig1c <- p2 + p3 + plot_layout(guides = "collect") &
  plot_annotation(
    title = "Atlas Composition",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))
  )

fig1c

ggsave(
  file.path(project_root, "results", "figures", "bmo_lineage_counts_and_proportions.pdf"),
  fig1c,
  width = 9,
  height = 4
)

library(dplyr)
library(ggplot2)
library(scales)

cell_props <- Haematopoietic@meta.data %>%
  group_by(group, celltype) %>%
  summarise(cnt = n()) %>%
  mutate(perc = cnt / sum(cnt)) %>%
  ungroup()

cell_props_matrix <- cell_props %>%
  pivot_wider(names_from = celltype, values_from = perc, values_fill = 0) %>%
  column_to_rownames("group") %>%
  as.matrix()

dim(cell_props_matrix)
class(cell_props_matrix)  # "matrix"


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

ggplot(plot_data, aes(x = group, y = perc, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = custom_cols) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(y = "Percentage of Cells",
       x = "Group",
       fill = NULL)
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

library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

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

base_data <- Haematopoietic@meta.data %>%
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

data_fig2 <- base_data %>%
  dplyr::filter(group %in% c("Dynamic.25d", "Dynamic.31d"))

data_fig2$group <- factor(data_fig2$group, levels = c("Dynamic.25d", "Dynamic.31d"))

p2 <- ggplot(data_fig2, aes(x = group, y = perc, fill = celltype)) +
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

print(p1)
print(p2)

library(dplyr)
library(ggplot2)
library(scales)

cell_props <- Haematopoietic@meta.data %>%
  dplyr::group_by(group, celltype) %>%
  dplyr::summarise(cnt = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(perc = cnt / sum(cnt)) %>%
  dplyr::ungroup()

library(tidyr)
library(tibble)

cell_props_matrix <- cell_props %>%
  dplyr::select(group, celltype, perc) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = perc, values_fill = 0) %>%
  tibble::column_to_rownames("group") %>%
  as.matrix()

print(dim(cell_props_matrix))
head(cell_props_matrix)

library(writexl)
library(dplyr)
library(tibble)

df_to_save <- as.data.frame(cell_props_matrix) %>%
  rownames_to_column("Group")

write_xlsx(df_to_save, path = file.path(project_root, "results", "tables", "BMO_Haematopoietic_cell_type_proportions.xlsx"))


Stromal <- readRDS(file.path(project_root, "data", "processed", "Stromal.rds"))

Stromal@meta.data$group <- ""
table(Stromal@meta.data$orig.ident)
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Dynamic.25d.1")] <- "Dynamic.25d"
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Dynamic.25d.2")] <- "Dynamic.25d"
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Dynamic.25d.3")] <- "Dynamic.25d"
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Dynamic.31d")] <- "Dynamic.31d"
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Static.25d.1")] <- "Static.25d"
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Static.25d.2")] <- "Static.25d"
Stromal@meta.data$group[Stromal@meta.data$orig.ident %in% c("Static.25d.3")] <- "Static.25d"
table(Stromal@meta.data$group)

Stromal@meta.data$celltype <- ""
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("41", "28", "2", "12", "13", "40")] <- "Osteoblasts"
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("15", "16", "22", "23")] <- "Osteoblasts pre."
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("31", "33")] <- "Cycling Osteoblasts pre."
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("7",  "24", "34")] <- "Osteochondral pre."
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("32", "38")] <- "Chondrocytes"
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("1", "3", "5", "6", "8", "26", "29", "35", "11", "17", "19", "0", "4", "10", "14", "20", "30")] <- "MSCs"
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("9", "37")] <- "Pericytes"
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("21")] <- "Adipocytes pre."
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("18", "27", "39")] <- "Proliferating MSCs"
Stromal@meta.data$celltype[Stromal@meta.data$seurat_clusters %in% c("25", "36")] <- "Endothelial Cells"
table(Stromal@meta.data$celltype)

library(ggplot2)
library(dplyr)
library(reshape2)

lineage_order <- c(
  "MSCs", "Proliferating MSCs",
  "Osteochondral pre.", "Chondrocytes",
  "Osteoblasts pre.", "Cycling Osteoblasts pre.", "Osteoblasts",
  "Adipocytes pre.",
  "Pericytes", "Endothelial Cells"
)

Stromal$celltype <- factor(Stromal$celltype, levels = lineage_order)

custom_cols <- c(
  "MSCs" = "#D9E3E4",
  "Proliferating MSCs" = "#E39A94",
  "Osteochondral pre." = "#B3CDE3",
  "Chondrocytes" = "#80B1D3",
  "Osteoblasts pre." = "#CCEBC5",
  "Cycling Osteoblasts pre." = "#B2DF8A",
  "Osteoblasts" = "#33A02C",
  "Adipocytes pre." = "#FDB462",
  "Pericytes" = "#F2D7D0",
  "Endothelial Cells" = "#BC80BD"
)

DimPlot(Stromal, reduction = "tsne", group.by = "celltype",
        label = TRUE, repel = TRUE, cols = custom_cols) +
  theme(aspect.ratio = 1)

base_data <- Stromal@meta.data %>%
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

data_fig2 <- base_data %>%
  dplyr::filter(group %in% c("Dynamic.25d", "Dynamic.31d"))

data_fig2$group <- factor(data_fig2$group, levels = c("Dynamic.25d", "Dynamic.31d"))

p2 <- ggplot(data_fig2, aes(x = group, y = perc, fill = celltype)) +
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

print(p1)
print(p2)

library(dplyr)
library(ggplot2)
library(scales)

cell_props <- Stromal@meta.data %>%
  dplyr::group_by(group, celltype) %>%
  dplyr::summarise(cnt = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(perc = cnt / sum(cnt)) %>%
  dplyr::ungroup()

library(tidyr)
library(tibble)

cell_props_matrix <- cell_props %>%
  dplyr::select(group, celltype, perc) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = perc, values_fill = 0) %>%
  tibble::column_to_rownames("group") %>%
  as.matrix()

print(dim(cell_props_matrix))
head(cell_props_matrix)

library(writexl)
library(dplyr)
library(tibble)

df_to_save <- as.data.frame(cell_props_matrix) %>%
  rownames_to_column("Group")

write_xlsx(df_to_save, path = file.path(project_root, "results", "tables", "BMO_Stromal_cell_type_proportions.xlsx"))

# Save the final annotated stromal object and a combined object used by CellChat.
saveRDS(
  Stromal,
  file.path(project_root, "data", "processed", "Stromal_clean_2_celltype.rds")
)
BMO_combined <- merge(Haematopoietic, y = Stromal)
saveRDS(
  BMO_combined,
  file.path(project_root, "data", "processed", "BMO_combined.rds")
)


plot_data <- FetchData(Stromal, vars = c("PIEZO1", "group", "celltype"))

head(plot_data)

library(tidyr)
plot_data <- plot_data %>%
  separate(group, into = c("Condition", "Day"), sep = "\\.", remove = FALSE)

plot_data_sub <- plot_data %>%
  filter(Day %in% c("25d"))

p <- ggplot(plot_data_sub, aes(x = Condition, y = PIEZO1, fill = Condition)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  facet_grid(Day ~ celltype) +
  stat_compare_means(comparisons = list(c("Dynamic", "Static")),
                     label = "p.format",
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8))

p

p <- ggplot(plot_data_sub, aes(x = Condition, y = PIEZO1, fill = Condition)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +

  stat_summary(fun = "mean", geom = "point",
               shape = 21, size = 2.5, fill = "black", color = "black") +

  stat_summary(fun = "mean", geom = "text",
               aes(label = after_stat(round(y, 3))),
               vjust = -1.5, size = 3, color = "black") +

  facet_grid(Day ~ celltype) +
  stat_compare_means(comparisons = list(c("Dynamic", "Static")),
                     label = "p.format",
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8))

p

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_annotate_bmo.txt")
)

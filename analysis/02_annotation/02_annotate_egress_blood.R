#!/usr/bin/env Rscript

# Egress-blood integration and cell-type annotation
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. merge the four post-doublet egress-blood objects;
#   2. integrate samples with CCA and recluster the combined object;
#   3. inspect canonical lineage markers and assign cluster identities;
#   4. add dynamic/static group labels; and
#   5. export the annotated object, marker plot, and composition table.

set.seed(20260730)
project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Seurat)
  library(tidyverse)
  library(writexl)
})

# The original analysis intentionally used the post-doublet objects for this
# integration. Ambient-RNA-corrected objects remain available as an alternative.
sample_files <- c(
  "blood.Dynamic.31d" = file.path(project_root, "data", "processed", "preprocessing", "EGRESS_BLOOD", "blood.Dynamic.31d", "02_post_doublet_removal.rds"),
  "blood.Dynamic2.31d" = file.path(project_root, "data", "processed", "preprocessing", "EGRESS_BLOOD", "blood.Dynamic2.31d", "02_post_doublet_removal.rds"),
  "blood.Static.31d" = file.path(project_root, "data", "processed", "preprocessing", "EGRESS_BLOOD", "blood.Static.31d", "02_post_doublet_removal.rds"),
  "blood.Static.33d" = file.path(project_root, "data", "processed", "preprocessing", "EGRESS_BLOOD", "blood.Static.33d", "02_post_doublet_removal.rds")
)

missing_sample_files <- sample_files[!file.exists(sample_files)]
if (length(missing_sample_files) > 0) {
  stop("Missing preprocessed egress-blood objects: ", paste(missing_sample_files, collapse = ", "))
}

list_scRNA <- lapply(sample_files, readRDS)

scRNA_all <- merge(list_scRNA[[1]], y = list_scRNA[-1], add.cell.ids = names(sample_files))

rm(list_scRNA)

#scRNA_all[["RNA"]] <- split(scRNA_all[["RNA"]], f = scRNA_all$orig.ident)

scRNA_all <- NormalizeData(scRNA_all)
scRNA_all <- FindVariableFeatures(scRNA_all, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(scRNA_all)
scRNA_all <- ScaleData(scRNA_all, features = all.genes)

scRNA_all <- RunPCA(scRNA_all, features = VariableFeatures(object = scRNA_all))

#saveRDS(scRNA_all, file = file.path(project_root, "data", "processed", "1_scRNA_all_PCA.rds"))
#scRNA_all <- readRDS(file.path(project_root, "data", "processed", "1_scRNA_all_PCA.rds"))

scRNA_all <- IntegrateLayers(
  object = scRNA_all,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE
)

scRNA_all[["RNA"]] <- JoinLayers(scRNA_all[["RNA"]])

saveRDS(
  scRNA_all,
  file.path(project_root, "data", "processed", "2_blood_scRNA_all_CCA_before_annotation.rds")
)

DimHeatmap(scRNA_all, dims = 1:30, cells = 500, balanced = TRUE)
ElbowPlot(scRNA_all, ndims = 50)

scRNA_all <- FindNeighbors(scRNA_all, reduction = "integrated.cca", dims = 1:20) #20
scRNA_all <- FindClusters(scRNA_all, resolution = 2)

scRNA_all <- RunUMAP(scRNA_all, reduction = "integrated.cca", dims = 1:20)
scRNA_all <- RunTSNE(scRNA_all, reduction = "integrated.cca", dims = 1:20)

p1 <- DimPlot(scRNA_all, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
p2 <- DimPlot(scRNA_all, reduction = "umap", label = TRUE) + ggtitle("UMAP")

p1+p2

p1 <- DimPlot(scRNA_all, reduction = "umap", group.by = "orig.ident") + ggtitle("Samples (CCA Integrated)")
p2 <- DimPlot(scRNA_all, reduction = "umap", label = TRUE) + ggtitle("Clusters")

p1 + p2

saveRDS(scRNA_all, file = file.path(project_root, "data", "processed", "4_Final_Integrated_Blood_CCA.rds"))

DimPlot(scRNA_all, reduction = "umap",label = TRUE, group.by = "orig.ident")
DimPlot(scRNA_all, reduction = "umap",label = TRUE, split.by = "orig.ident")

all_haem_markers <- c("CD34", "ITGA6", "ATXN1",   #HSC
                      "PPBP","PF4", #Megakaryocyte
                      "TFRC","GYPA","SLC4A1","KLF1", #Erythroid
                      "MPO", "ELANE", #Neutrophil
                      "EPX",  #Eosinophil
                      "FCER1A",  # Basophil
                      "TPSAB1", "KIT", # Mast
                      "S100A8", "S100A9", "CSF3R", # Monocyte
                      "CD14", "CD68", "CD163",  #Macrophage
                      "HLA-DPA1","CD1C","CLEC10A") # DCs


# HSC/MPP
FeaturePlot(scRNA_all, features = c("CD34", "THY1", "MECOM"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

FeaturePlot(scRNA_all, features = c("CD34", "ITGA6", "ATXN1"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")


# Erythroid
FeaturePlot(scRNA_all, features = c("TFRC","GYPA","SLC4A1","KLF1"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

FeaturePlot(scRNA_all, features = c("PTPRC"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# MK
FeaturePlot(scRNA_all, features = c("PPBP","PF4"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# Monocyte
FeaturePlot(scRNA_all, features = c("S100A8", "S100A9", "CSF3R"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# Macrophage
FeaturePlot(scRNA_all, features = c("CD14","CD68","CD163"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# DC
FeaturePlot(scRNA_all, features = c("HLA-DPA1","CD1C","CLEC10A"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# Neutrophil
FeaturePlot(scRNA_all, features = c("MPO", "ELANE"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# Basophil
FeaturePlot(scRNA_all, features = c("FCER1A"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# Eosinophil
FeaturePlot(scRNA_all, features = c("EPX"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# Mast
FeaturePlot(scRNA_all, features = c("TPSAB1"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

# NK
FeaturePlot(scRNA_all, features = c("NCAM1","FCGR3A","KIR3DL1","KLRC1"), max.cutoff = 3, min.cutoff=0,
            cols = c("grey","red"), reduction = "umap")

scRNA_all <- RenameIdents(scRNA_all,
                    `33`="HSC",

                    `0`="Late Erythroid", `1`="Late Erythroid", `2`="Late Erythroid", `28`="Late Erythroid",
                    `31`="Late Erythroid", `23`="Late Erythroid", `8`="Late Erythroid", `7`="Late Erythroid",
                    `6`="Late Erythroid", `15`="Late Erythroid", `10`="Late Erythroid", `16`="Late Erythroid",
                    `20`="Late Erythroid", `30`="Late Erythroid", `13`="Late Erythroid", `4`="Late Erythroid",
                    `22`="Late Erythroid", `25`="Late Erythroid", `12`="Late Erythroid", `11`="Late Erythroid",
                    `26`="Late Erythroid", `32`="Late Erythroid", `19`="Late Erythroid", `9`="Late Erythroid",
                    `14`="Late Erythroid",

                    `17`="Early Erythroid",

                    `34`="Megakaryocyte", `18`="Megakaryocyte", `38`="Megakaryocyte",

                    `5`="Basophil", `21`="Basophil", `35`="Basophil",

                    `27`="Eosinophil",

                    `36` ="Neutrophil",

                    `29` = "Mast",

                    `24` ="Monocyte",

                    `3` = "Macrophage",

                    `37` = "DC")

scRNA_all$celltype <- Idents(scRNA_all)
head(scRNA_all@meta.data)

levels_order <- c(
  "HSC",
  "Early Erythroid",
  "Late Erythroid",
  "Megakaryocyte",
  "Monocyte",
  "Macrophage",
  "DC",
  "Neutrophil",
  "Basophil",
  "Eosinophil",
  "Mast"
)

scRNA_all$celltype <- factor(scRNA_all$celltype, levels = levels_order)
Idents(scRNA_all) <- "celltype"

DimPlot(scRNA_all, label = TRUE)

saveRDS(scRNA_all, file = file.path(project_root, "data", "processed", "BMO_egress_Blood.rds"))
BMO_egress_Blood <- readRDS(file.path(project_root, "data", "processed", "BMO_egress_Blood.rds"))

# Plot lineage markers only after the cell-type metadata has been assigned.
haem_ipsc_dot <- DotPlot(
  BMO_egress_Blood,
  features = all_haem_markers,
  cols = c("white", "#EB0616", "#EB2F06"),
  col.min = 0,
  dot.min = 0,
  dot.scale = 6,
  group.by = "celltype",
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  coord_fixed()

ggsave(
  file.path(project_root, "results", "figures", "egress_blood_lineage_markers.pdf"),
  haem_ipsc_dot,
  width = 8,
  height = 6
)

###QC
VlnPlot(scRNA_all, features = c("nFeature_RNA", "percent.mt"), group.by = "seurat_clusters", pt.size = 0)
Stromal_clean <- subset(scRNA_all, idents = c("5", "7"), invert = TRUE)

custom_cols <- c(
  "HSC"              = "#555555",

  "Megakaryocyte"    = "#FFD92F",

  "Early Erythroid"  = "#FB9A99",
  "Late Erythroid"   = "#E41A1C",

  "Neutrophil"       = "#1F78B4",
  "Eosinophil"       = "#A6CEE3",

  "Basophil"         = "#33A02C",
  "Mast"             = "#B2DF8A",

  "Monocyte"         = "#FDBF6F",
  "Macrophage"       = "#B15928",
  "DC"               = "#CAB2D6"
)

levels_order <- c("HSC", "Megakaryocyte", "Early Erythroid", "Late Erythroid",
                  "Neutrophil", "Eosinophil", "Basophil", "Mast",
                  "Monocyte", "Macrophage", "DC")

scRNA_all$celltype <- factor(scRNA_all$celltype, levels = levels_order)

DimPlot(scRNA_all, reduction = "umap", group.by = "celltype",
        label = TRUE, repel = TRUE, cols = custom_cols) +
  theme(aspect.ratio = 1)


BMO_egress_Blood@meta.data$group <- ""
table(BMO_egress_Blood@meta.data$orig.ident)
BMO_egress_Blood@meta.data$group[BMO_egress_Blood@meta.data$orig.ident %in% c("blood.Dynamic.31d")] <- "Dynamic"
BMO_egress_Blood@meta.data$group[BMO_egress_Blood@meta.data$orig.ident %in% c("blood.Dynamic2.31d")] <- "Dynamic"
BMO_egress_Blood@meta.data$group[BMO_egress_Blood@meta.data$orig.ident %in% c("blood.Static.31d")] <- "Static"
BMO_egress_Blood@meta.data$group[BMO_egress_Blood@meta.data$orig.ident %in% c("blood.Static.33d")] <- "Static"

table(BMO_egress_Blood@meta.data$group)

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

BMO_egress_Blood$celltype <- factor(BMO_egress_Blood$celltype, levels = names(custom_cols))

DimPlot(BMO_egress_Blood,
        reduction = "tsne",
        group.by = "celltype",
        label = TRUE,
        cols = custom_cols,
        pt.size = 0.8) +
  theme(aspect.ratio = 1)

lineage_order <- c("HSC", "Megakaryocyte", "Early Erythroid", "Late Erythroid",
                  "Neutrophil", "Eosinophil", "Basophil", "Mast",
                  "Monocyte", "Macrophage", "DC")

BMO_egress_Blood$celltype <- factor(BMO_egress_Blood$celltype, levels = lineage_order)

custom_cols <- c(
  "HSC"              = "#555555",
  "Megakaryocyte"    = "#FFD92F",
  "Early Erythroid"  = "#FB9A99",
  "Late Erythroid"   = "#E41A1C",
  "Neutrophil"       = "#1F78B4",
  "Eosinophil"       = "#A6CEE3",
  "Basophil"         = "#33A02C",
  "Mast"             = "#B2DF8A",
  "Monocyte"         = "#FDBF6F",
  "Macrophage"       = "#B15928",
  "DC"               = "#CAB2D6"
)

DimPlot(BMO_egress_Blood, reduction = "umap", group.by = "celltype",
        label = TRUE, repel = TRUE, cols = custom_cols) +
  theme(aspect.ratio = 1)

base_data <- BMO_egress_Blood@meta.data %>%
  dplyr::filter(group != "") %>%
  dplyr::count(group, celltype) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(perc = n / sum(n)) %>%
  dplyr::ungroup()

base_data$celltype <- factor(base_data$celltype, levels = names(custom_cols))

base_data$group <- factor(base_data$group, levels = c("Static", "Dynamic"))

p1 <- ggplot(base_data, aes(x = group, y = perc, fill = celltype)) +
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

library(dplyr)
library(ggplot2)
library(scales)

cell_props <- BMO_egress_Blood@meta.data %>%
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

write_xlsx(df_to_save, path = file.path(project_root, "results", "tables", "BMO_egress_Blood_cell_type_proportions.xlsx"))

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_annotate_egress_blood.txt")
)

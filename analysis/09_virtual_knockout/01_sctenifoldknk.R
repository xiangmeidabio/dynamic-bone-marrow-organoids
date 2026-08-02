#!/usr/bin/env Rscript

# scTenifoldKnk virtual-knockout analysis
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. select highly variable genes plus the perturbation target;
#   2. infer scTenifoldKnk networks and simulate target-gene knockout;
#   3. repeat the perturbation within dynamic and static erythroid subsets; and
#   4. export significant genes and publication figures.

set.seed(20260730)
project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

library(remotes)
library(scTenifoldKnk)
library(slingshot)
library(Seurat)
library(devtools)
library(Seurat)
library(cowplot)
library(ggplot2)
  library(Matrix)
  library(openxlsx)
library(dplyr)
library(tradeSeq)
library(RColorBrewer)

scRNA <- readRDS(file.path(project_root, "data", "processed", "Blood_Erythroid.rds"))

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

levels_order <- c("HSC", "Megakaryocyte", "Erythroid",
                  "Neutrophil", "Eosinophil", "Basophil", "Mast",
                  "Monocyte", "Macrophage", "DC")

scRNA$celltype <- factor(scRNA$celltype, levels = levels_order)

DimPlot(scRNA, reduction = "umap", group.by = "celltype",
        label = TRUE, repel = TRUE, cols = custom_cols) +
  theme(aspect.ratio = 1)

VlnPlot(scRNA, features = "BSG", split.by = "condition")

library(ggpubr)
library(ggplot2)

p <- VlnPlot(scRNA, features = "BSG", split.by = "condition", group.by = "condition", y.max = 5)

condition_levels <- unique(na.omit(as.character(scRNA$condition)))
if (length(condition_levels) != 2) {
  stop("The condition column must contain exactly two groups for this comparison.")
}
my_comparisons <- list(condition_levels)

p + stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",
                       label.y = 4.5) +
  stat_summary(fun = "mean", geom = "point", color = "red", size = 2,
               position = position_dodge(0.9))

countMat <- GetAssayData(scRNA, layer = "counts")

scRNA <- FindVariableFeatures(
  object = scRNA,
  selection.method = "vst",
  nfeatures = 2000
)

hvgs <- VariableFeatures(scRNA)
target_gene <- "BSG"

genes_use <- unique(c(target_gene, hvgs))
genes_use <- intersect(genes_use, rownames(countMat))

# 2000 HVGs + target gene, all cells
data <- countMat[genes_use, , drop = FALSE]

dim(data) #2001 21608

result <- scTenifoldKnk(
  countMatrix = data,
  gKO = target_gene,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = 500
)

#saveRDS(result, file = file.path(project_root, "data", "processed", "scTenifoldKnk_Blood_Erythroid_result_2000genes_ncells5000_nNet10.rds"))
saveRDS(result, file = file.path(project_root, "data", "processed", "scTenifoldKnk_Blood_Erythroid_result_2000genes_ncells500_nNet10.rds"))

library(igraph)
plotKO(result,
       gKO ="BSG",
       annotate = FALSE,
       q = 0.99)

library(scTenifoldKnk)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(remotes)
library(forcats)

df=result$diffRegulation
df_KO <- result[["diffRegulation"]]
df_KO$log2FC <- log2(df_KO$FC)

genes_focus <- c("HBB", "HBA1", "HBA2", "HBD")

df_focus <- df_KO %>%
  filter(gene %in% genes_focus)

top20 <- df_KO %>%
  filter(p.value < 0.05 & log2FC > 0) %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 50)

p <- df_focus %>%
  mutate(gene = fct_reorder(gene, log2FC)) %>%
  ggplot(aes(x = gene, y = log2FC, fill = log2FC)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  geom_text(
    aes(label = round(log2FC, 1)),
    hjust = -0.15,
    size = 3.6,
    color = "grey20"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
  labs(
    title = "Differentially Regulated Genes",
    subtitle = paste0("Ranked by log2FC (n = ", nrow(df_KO), ")"),
    x = NULL,
    y = "log2FC"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40")
  );
p

significant_threshold <- 0.05
differential_genes <- result$diffRegulation %>%
  mutate(log_fold_change = log2(FC)) %>%
  filter(gene != target_gene)
  #%>%
  #filter(p.value < significant_threshold)

numeric_columns <- 2:7
differential_genes[, numeric_columns] <- sapply(differential_genes[, numeric_columns], as.numeric)

output_filename <- paste0('1_', target_gene, "_virtual_knockout_results.xlsx")
write.xlsx(
  differential_genes,
  file.path(project_root, "results", "tables", output_filename)
)
head(differential_genes)

genes_focus <- c("HBB", "HBA1", "HBA2", "HBD")

df_focus <- differential_genes %>%
  filter(gene %in% genes_focus)

top_perturbed <- differential_genes %>%
  arrange(desc(log_fold_change)) %>%
  slice_head(n = 20)

ggplot(top_perturbed,
       aes(x = reorder(gene, log_fold_change),
           y = log_fold_change,
           size = -log10(p.value),
           color = log_fold_change)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#3498DB", high = "#E74C3C",
                       name = "Perturbation\nMagnitude") +
  coord_flip() +
  labs(title = "Top Perturbed Genes",
       subtitle = paste("KO:", target_gene),
       x = "Gene",
       y = "Log2(Perturbation Magnitude)",
       size = "-log10(p-value)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  scale_size_continuous(range = c(3, 8))

ggsave(file.path(project_root, "results", "figures", "5_perturbation_dot_plot.pdf"), width = 5.6, height = 6.5)


table(scRNA@meta.data[["group"]])

scRNA@meta.data$condition <- ""
scRNA@meta.data$condition[scRNA@meta.data$group %in% c("blood.Dynamic.31d")] <- "Dynamic"
scRNA@meta.data$condition[scRNA@meta.data$group %in% c("blood.Static.31d", "blood.Static.33d")] <- "Static"
table(scRNA@meta.data[["condition"]])

Dynamic <- subset(scRNA, subset = condition %in% c("Dynamic"))
Static <- subset(scRNA, subset = condition %in% c("Static"))

Dynamic_countMat <- GetAssayData(Dynamic, layer = "counts")
Static_countMat <- GetAssayData(Static, layer = "counts")

Dynamic_scRNA <- FindVariableFeatures(
  object = Dynamic,
  selection.method = "vst",
  nfeatures = 1000
)
Static_scRNA <- FindVariableFeatures(
  object = Static,
  selection.method = "vst",
  nfeatures = 1000
)

Dynamic_hvgs <- VariableFeatures(Dynamic_scRNA)
Static_hvgs <- VariableFeatures(Static_scRNA)

target_gene <- "BSG"

Dynamic_genes_use <- unique(c(target_gene, Dynamic_hvgs))
Static_genes_use <- unique(c(target_gene, Static_hvgs))

Dynamic_genes_use <- intersect(Dynamic_genes_use, rownames(Dynamic_countMat))
Dynamic_data <- Dynamic_countMat[Dynamic_genes_use, , drop = FALSE]

Static_genes_use <- intersect(Static_genes_use, rownames(Static_countMat))
Static_data <- Static_countMat[Static_genes_use, , drop = FALSE]

Dynamic_cell_mean <- Matrix::colMeans(Dynamic_data)
Static_cell_mean <- Matrix::colMeans(Static_data)

Dynamic_top_cell_idx <- order(Dynamic_cell_mean, decreasing = TRUE)[1:min(1000, ncol(Dynamic_data))]
Dynamic_data <- Dynamic_data[, Dynamic_top_cell_idx, drop = FALSE]

Static_top_cell_idx <- order(Static_cell_mean, decreasing = TRUE)[1:min(1000, ncol(Static_data))]
Static_data <- Static_data[, Static_top_cell_idx, drop = FALSE]

dim(Dynamic_data)
dim(Static_data)

Dynamic_result <- scTenifoldKnk(
  countMatrix = Dynamic_data,
  gKO = target_gene,
  qc = FALSE,
  qc_mtThreshold = 0.1,
  qc_minLSize = 1000,
  nc_nNet = 10,
  nc_nCells = 500
)

Static_result <- scTenifoldKnk(
  countMatrix = Static_data,
  gKO = target_gene,
  qc = FALSE,
  qc_mtThreshold = 0.1,
  qc_minLSize = 1000,
  nc_nNet = 10,
  nc_nCells = 500
)


library(igraph)
plotKO(Dynamic_result,
       gKO ="BSG",
       annotate = FALSE,
       q = 0.99)

plotKO(Static_result,
       gKO ="BSG",
       annotate = FALSE,
       q = 0.99)

library(scTenifoldKnk)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(remotes)
library(forcats)

Dynamic_df_KO <- Dynamic_result[["diffRegulation"]]
Dynamic_df_KO$log2FC <- log2(Dynamic_df_KO$FC)

Static_df_KO <- Static_result[["diffRegulation"]]
Static_df_KO$log2FC <- log2(Static_df_KO$FC)

genes_focus <- c("HBB", "HBA1", "HBA2", "HBD")

Dynamic_df_focus <- Dynamic_df_KO %>%
  filter(gene %in% genes_focus)

Static_df_focus <- Static_df_KO %>%
  filter(gene %in% genes_focus)

#Dynamic
Dynamic_p <- Dynamic_df_focus %>%
  mutate(gene = fct_reorder(gene, log2FC)) %>%
  ggplot(aes(x = gene, y = log2FC, fill = log2FC)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  geom_text(
    aes(label = round(log2FC, 1)),
    hjust = -0.15,
    size = 3.6,
    color = "grey20"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
  labs(
    title = "Differentially Regulated Genes",
    subtitle = paste0("Ranked by log2FC (n = ", nrow(df_KO), ")"),
    x = NULL,
    y = "log2FC"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40")
  );

Dynamic_p

#Static
Static_p <- Static_df_focus %>%
  mutate(gene = fct_reorder(gene, log2FC)) %>%
  ggplot(aes(x = gene, y = log2FC, fill = log2FC)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  geom_text(
    aes(label = round(log2FC, 1)),
    hjust = -0.15,
    size = 3.6,
    color = "grey20"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
  labs(
    title = "Differentially Regulated Genes",
    subtitle = paste0("Ranked by log2FC (n = ", nrow(df_KO), ")"),
    x = NULL,
    y = "log2FC"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40")
  );

Static_p

saveRDS(Dynamic_result, file = file.path(project_root, "data", "processed", "scTenifoldKnk_Dynamic_Blood_Erythroid_result.rds"))
saveRDS(Static_result, file = file.path(project_root, "data", "processed", "scTenifoldKnk_Static_Blood_Erythroid_result.rds"))


target_order <- c("Static.25d", "Dynamic.25d", "Dynamic.31d")

Haematopoietic@meta.data$group <- factor(
  Haematopoietic@meta.data$group,
  levels = target_order
)

p <- VlnPlot(Haematopoietic, features = "BSG", split.by = "group", group.by = "group", y.max = 5)
p <- VlnPlot(Haematopoietic, features = "GDF15", split.by = "group", group.by = "group", y.max = 7.1)

table(Haematopoietic@meta.data$group)

my_comparisons <- list(c("Dynamic.25d", "Static.25d"))

p + stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "wilcox.test",
                       label.y = 6.5) +
  stat_summary(fun = "mean", geom = "point", color = "red", size = 2,
               position = position_dodge(0.9)) +
  theme(legend.position = "none")


Ery_dyn25 <- subset(Haematopoietic, subset = group == "Dynamic.25d" & celltype == "Erythroid")
dim(Ery_dyn25) #[1] 32104  6056

sum(GetAssayData(Ery_dyn25, layer = "data")["BSG", ] > 0) #[1] 4965

Ery_dyn25 <- NormalizeData(Ery_dyn25)

Ery_dyn25 <- FindVariableFeatures(Ery_dyn25, selection.method = "vst", nfeatures = 3000)

hvgs <- VariableFeatures(Ery_dyn25)
target_gene <- "BSG"
genes_use <- unique(c(target_gene, hvgs))

data_knk <- GetAssayData(Ery_dyn25, layer = "data")[genes_use, ]

result <- scTenifoldKnk(
  countMatrix = data_knk,
  gKO = target_gene,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = 500
)

saveRDS(result, file = file.path(project_root, "data", "processed", "scTenifoldKnk_bmo_Dynamic_Blood_Erythroid_result.rds"))
Dynamic_result <- readRDS(file.path(project_root, "data", "processed", "scTenifoldKnk_bmo_Dynamic_Blood_Erythroid_result.rds"))

significant_threshold <- 0.05
differential_genes <- Dynamic_result$diffRegulation %>%
  mutate(log_fold_change = log2(FC)) %>%
  filter(gene != target_gene) %>%
  filter(p.value < significant_threshold)

numeric_columns <- 2:7
differential_genes[, numeric_columns] <- sapply(differential_genes[, numeric_columns], as.numeric)

library(openxlsx)

output_filename <- paste0('1_Dynamic_', target_gene, "_significant_virtual_knockout_results.xlsx")
write.xlsx(
  differential_genes,
  file.path(project_root, "results", "tables", output_filename)
)
head(differential_genes)

genes_focus <- c("HBB", "HBA1", "HBA2", "HBD")

df_focus <- differential_genes %>%
  filter(gene %in% genes_focus)

top_perturbed <- differential_genes %>%
  arrange(desc(log_fold_change)) %>%
  slice_head(n = 20)

ggplot(top_perturbed,
       aes(x = reorder(gene, log_fold_change),
           y = log_fold_change,
           size = -log10(p.value),
           color = log_fold_change)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#3498DB", high = "#E74C3C",
                       name = "Perturbation\nMagnitude") +
  coord_flip() +
  labs(title = "Top Perturbed Genes",
       subtitle = paste("KO:", target_gene),
       x = "Gene",
       y = "Log2(Perturbation Magnitude)",
       size = "-log10(p-value)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  scale_size_continuous(range = c(3, 8))

ggsave(file.path(project_root, "results", "figures", "2_Dynamic_perturbation_dot_plot.pdf"), width = 5.6, height = 6.5)
ggsave(file.path(project_root, "results", "figures", "2_Dynamic_pvalue_perturbation_dot_plot.pdf"), width = 5.6, height = 4)

table(Haematopoietic@meta.data$group)

Ery_sta25 <- subset(Haematopoietic, subset = group == "Static.25d" & celltype == "Erythroid")
dim(Ery_sta25) #[1] 32104  1314

sum(GetAssayData(Ery_sta25, layer = "data")["BSG", ] > 0) #[1] 1063

Ery_sta25 <- NormalizeData(Ery_sta25)

Ery_sta25 <- FindVariableFeatures(Ery_sta25, selection.method = "vst", nfeatures = 3000)

hvgs <- VariableFeatures(Ery_sta25)
target_gene <- "BSG"
genes_use <- unique(c(target_gene, hvgs))

data_knk <- GetAssayData(Ery_sta25, layer = "data")[genes_use, ]
dim(data_knk) #[1] 3001 1314

result <- scTenifoldKnk(
  countMatrix = data_knk,
  gKO = target_gene,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = 500
)

saveRDS(result, file = file.path(project_root, "data", "processed", "scTenifoldKnk_bmo_Static_Blood_Erythroid_result.rds"))
Static_result <- readRDS(file.path(project_root, "data", "processed", "scTenifoldKnk_bmo_Static_Blood_Erythroid_result.rds"))

significant_threshold <- 0.05
differential_genes <- Static_result$diffRegulation %>%
  mutate(log_fold_change = log2(FC)) %>%
  filter(gene != target_gene) %>%
  filter(p.value < significant_threshold)

numeric_columns <- 2:7
differential_genes[, numeric_columns] <- sapply(differential_genes[, numeric_columns], as.numeric)

output_filename <- paste0('1_Static_result', target_gene, "_significant_virtual_knockout_results.xlsx")
write.xlsx(
  differential_genes,
  file.path(project_root, "results", "tables", output_filename)
)
head(differential_genes)

genes_focus <- c("HBB", "HBA1", "HBA2", "HBD")
genes_focus <- c("GDF15")
genes_focus <- c("HBM","FTL")

df_focus <- differential_genes %>%
  filter(gene %in% genes_focus)

top_perturbed <- differential_genes %>%
  arrange(desc(log_fold_change)) %>%
  slice_head(n = 20)

top_perturbed$p.value[top_perturbed$p.value == 0] <- 1e-300

ggplot(top_perturbed,
       aes(x = reorder(gene, log_fold_change),
           y = log_fold_change,
           size = -log10(p.value),
           color = log_fold_change)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "#3498DB", high = "#E74C3C",
                       name = "Perturbation\nMagnitude") +
  coord_flip() +
  labs(title = "Top Perturbed Genes",
       subtitle = paste("KO:", target_gene),
       x = "Gene",
       y = "Log2(Perturbation Magnitude)",
       size = "-log10(p-value)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  scale_size_continuous(range = c(3, 8))

ggsave(file.path(project_root, "results", "figures", "2_Static_perturbation_dot_plot.pdf"), width = 5.6, height = 6.5)
ggsave(file.path(project_root, "results", "figures", "2_Static_pvalue_perturbation_dot_plot.pdf"), width = 5.6, height = 4)

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_sctenifoldknk.txt")
)

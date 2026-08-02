#!/usr/bin/env Rscript

# pySCENIC downstream analysis and visualization
#
# Run from the repository root. Statements remain in analysis order so that
# every biological decision is visible and auditable.
#
# Workflow:
#   1. read the pySCENIC loom file and extract regulons and AUC scores;
#   2. compare regulon activity between experimental groups;
#   3. identify erythroid-maturation targets of selected regulons; and
#   4. generate tabular exports and regulatory-network visualizations.

project_root <- normalizePath(Sys.getenv("BMO_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dir.create(file.path(project_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

#############################################################################


#############################################################################
#############################################################################

library(AUCell)
library(SCopeLoomR)
library(SCENIC)

library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)

#############################################################################
#############################################################################

sce_SCENIC <- open_loom(
  file.path(project_root, "results", "objects", "grn", "BMO_Erythroid_SCENIC.loom")
)


#############################################################################
#############################################################################

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name = "Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)

class(regulons)

head(regulons)

#############################################################################
#############################################################################

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name = "RegulonsAUC")

regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

dim(getAUC(regulonAUC))

getAUC(regulonAUC)[1:5, 1:5]

#############################################################################
#############################################################################

auc_matrix <- getAUC(regulonAUC)

auc_df <- as.data.frame(t(auc_matrix))

dim(auc_df)
head(auc_df[, 1:5])

#############################################################################
#############################################################################

regulon_var <- apply(auc_df, 2, var)
top_regulons <- names(sort(regulon_var, decreasing = TRUE))[1:20]

heatmap_mat <- as.matrix(auc_df[, top_regulons])

heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

pheatmap(
  t(heatmap_mat_scaled),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  main = "Top 20 Regulons AUC Heatmap"
)

#############################################################################
#############################################################################

regulon_name <- colnames(auc_df)[1]
regulon_name <- "CEBPB(+)"
regulon_name <- "MAFG(+)"

plot_df <- data.frame(
  Cell = rownames(auc_df),
  AUC = auc_df[, regulon_name]
)

ggplot(plot_df, aes(x = AUC)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30) +
  theme_bw() +
  labs(
    title = paste("AUC Distribution of", regulon_name),
    x = "AUC score",
    y = "Cell count"
  )

#############################################################################
#############################################################################

plot_df <- plot_df %>%
  arrange(desc(AUC)) %>%
  mutate(Cell_order = 1:n())

ggplot(plot_df, aes(x = Cell_order, y = AUC)) +
  geom_line(color = "tomato") +
  theme_bw() +
  labs(
    title = paste("AUC Ranking of", regulon_name),
    x = "Cells ranked by AUC",
    y = "AUC score"
  )

#############################################################################
#############################################################################

regulons[[1]]

regulons[["MAFG(+)"]]
regulons[["CEBPB(+)"]]


#############################################################################
#############################################################################

close_loom(sce_SCENIC)

#############################################################################
#############################################################################


#############################################################################
#############################################################################

target_gene <- "GDF15"

gdf15_regulons <- names(regulons)[sapply(regulons, function(x) target_gene %in% x)]

cat("Regulons containing", target_gene, ":\n")
print(gdf15_regulons)  #[1] "CEBPB(+)" "MAFG(+)"

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# auc_mat: rownames = regulons, colnames = cells
auc_df <- as.data.frame(t(auc_mat))
auc_df$Sample <- rownames(auc_df)

auc_df <- auc_df %>%
  mutate(
    Prefix = str_extract(Sample, "^[^_]+"),
    Group  = str_extract(Prefix, "^[^.]+"),
    Day    = str_extract(Prefix, "(?<=\\.)[^.]+")
  )

auc_long <- auc_df %>%
  pivot_longer(
    cols = -c(Sample, Prefix, Group, Day),
    names_to = "Regulon",
    values_to = "AUC"
  )

head(auc_long)
table(auc_long$Group)
table(auc_long$Day)
table(auc_long$Prefix)

library(dplyr)

auc_long <- auc_long %>%
  mutate(
    condition = case_when(
      Prefix %in% c("Dynamic.24d", "Dynamic.25d.1", "Dynamic.25d.2") ~ "Dynamic.25d",
      Prefix %in% c("Static.24d", "Static.25d.1", "Static.25d.2") ~ "Static.25d",
      TRUE ~ Prefix
    )
  )
table(auc_long$condition)

#key_regs <- c("CEBPB(+)", "MAFG(+)", "KLF1(+)", "GATA1(+)", "TAL1(+)", "NFE2(+)")
key_regs <- c("CEBPB(+)")
key_regs <- c("MAFG(+)")

colnames(plot_df)

plot_df <- auc_long %>%
  filter(Regulon %in% key_regs)

ggplot(plot_df, aes(x = condition, y = AUC, fill = condition)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ Regulon, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggpubr)
library(dplyr)
library(ggplot2)
library(ggpubr)

plot_df <- auc_long %>%
  filter(Regulon %in% key_regs) %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Static.25d", "Dynamic.25d", "Dynamic.31d")
    )
  )

my_comparisons <- list(
  c("Static.25d", "Dynamic.25d")
)

p <- ggplot(plot_df, aes(x = condition, y = AUC, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, width = 0.65) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif"
  ) +
  facet_wrap(~ Regulon, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Regulon AUC"
  )

p

library(dplyr)
library(ggplot2)
library(ggpubr)

plot_df <- plot_df %>%
  mutate(
    condition = factor(
      condition,
      levels = c("Static.25d", "Dynamic.25d", "Dynamic.31d")
    )
  )

my_comparisons <- list(
  c("Static.25d", "Dynamic.25d")
)

cols_use <- c(
  "Static.25d"  = "#D98B7E",
  "Dynamic.25d" = "#4DB36A",
  "Dynamic.31d" = "#F2B84B"
)

p <- ggplot(plot_df, aes(x = condition, y = AUC, fill = condition)) +
  geom_violin(
    width = 0.85,
    trim = FALSE,
    scale = "width",
    color = "black",
    linewidth = 0.35,
    alpha = 0.85
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.35,
    alpha = 0.85
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.35,
    size = 4
  ) +
  facet_wrap(~ Regulon, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = cols_use) +
  scale_color_manual(values = cols_use) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      size = 10
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10
    ),
    axis.title.y = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      color = "black",
      face = "bold",
      size = 11
    ),
    legend.position = "none",
    panel.spacing = unit(0.9, "lines")
  ) +
  labs(
    x = NULL,
    y = "Regulon AUC"
  )
p
ggsave(
  file.path(project_root, "results", "figures", "F6E.pdf"),
  plot = p,
  width = 4,
  height = 4,
  units = "in"
)

#sce_SCENIC <- open_loom("BMO_Erythroid_SCENIC.loom")

regulons <- get_regulons(sce_SCENIC, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons)

cebpb_targets <- regulons[["CEBPB(+)"]]
mafg_targets  <- regulons[["MAFG(+)"]]

ery_maturation_genes <- c(
   "GATA1","KLF1","TAL1","NFE2","ALAS2","AHSP","SLC4A1","EPB42","GYPA",
   "HBA1","HBA2","HBB","HBD","HBM","FECH","TFRC","BPGM","ANK1","SPTA1","SPTB"
 )

cebpb_hits <- intersect(cebpb_targets, ery_maturation_genes)
mafg_hits  <- intersect(mafg_targets, ery_maturation_genes)

cat("CEBPB regulon name:", cebpb_name, "\n")
cat("MAFG regulon name:", mafg_name, "\n\n")

cat("Number of CEBPB targets:", length(cebpb_targets), "\n")
cat("Number of MAFG targets:", length(mafg_targets), "\n\n")

cat("CEBPB targets overlapping erythroid maturation genes:\n")
print(cebpb_hits)

cat("\nMAFG targets overlapping erythroid maturation genes:\n")
print(mafg_hits)

targets_df <- bind_rows(
  data.frame(TF = "CEBPB", Target = cebpb_targets),
  data.frame(TF = "MAFG",  Target = mafg_targets)
) %>%
  mutate(Erythroid_maturation_related = Target %in% ery_maturation_genes)

write.csv(
  targets_df,
  file.path(project_root, "results", "tables", "CEBPB_MAFG_all_targets_annotated.csv"),
  row.names = FALSE
)

result_df <- bind_rows(
  data.frame(TF = "CEBPB", Gene = cebpb_hits),
  data.frame(TF = "MAFG",  Gene = mafg_hits)
)

write.csv(
  result_df,
  file.path(project_root, "results", "tables", "CEBPB_MAFG_erythroid_maturation_targets.csv"),
  row.names = FALSE
)


target_genes <- unique(c(
  "GYPA", "NFE2", "HBM", "HBA2", "HBA1", "AHSP", "FECH", "KLF1", "ALAS2",
  "TAL1", "BPGM"
))

target_genes

library(dplyr)
library(stringr)

exprMat <- get_dgem(sce_SCENIC)
exprMat_log <- log2(exprMat + 1)

cell_meta <- data.frame(Cell = colnames(exprMat), stringsAsFactors = FALSE) %>%
  mutate(
    Condition = str_extract(Cell, "^[^_]+")
  )

table(cell_meta$Condition)


cell_meta$Condition <- factor(
  cell_meta$Condition,
  levels = c("Static.25d","Dynamic.25d", "Dynamic.31d")
)

genes_use <- intersect(target_genes, rownames(exprMat))
genes_use
setdiff(target_genes, genes_use)

library(tidyr)
library(ggplot2)

plot_df <- exprMat[genes_use, , drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Cell") %>%
  left_join(cell_meta, by = "Cell") %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "Gene",
    values_to = "Expression"
  )

head(plot_df)


library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

edges <- data.frame(
  from = "CEBPB",
  to = c("HBA1", "HBA2", "HBM", "NFE2", "GYPA", "GDF15", "ALAS2", "KLF1", "FECH", "AHSP")
)

nodes <- data.frame(
  name = c("CEBPB", "HBA1", "HBA2", "HBM", "NFE2", "GYPA", "GDF15", "ALAS2", "KLF1", "FECH", "AHSP"),
  x = c(0, -3.2, -2.6, -1.0, 1.0, 2.6, 3.2, 2.6, 1.0, -1.0, -2.6),
  y = c(0,  0.0,  1.8,  3.2, 3.2, 1.8, 0.0,-1.8,-3.2,-3.2,-1.8)
)

scale_factor <- 1.45

nodes <- nodes %>%
  mutate(
    x = ifelse(name == "CEBPB", x, x * scale_factor),
    y = ifelse(name == "CEBPB", y, y * scale_factor),
    fill_color = case_when(
      name == "CEBPB" ~ "black",
      name == "GDF15" ~ "blue",
      TRUE ~ "red"
    ),
    label_x = ifelse(name == "CEBPB", x, x * 1.16),
    label_y = ifelse(name == "CEBPB", y, y * 1.16)
  )

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +
  geom_edge_arc(
    strength = 0.22,
    color = "grey70",
    width = 0.9
  ) +
  geom_node_point(
    aes(x = x, y = y, fill = I(fill_color)),
    shape = 21,
    size = 7,
    color = "grey35",
    stroke = 0.9
  ) +
  geom_node_text(
    aes(x = label_x, y = label_y, label = name),
    color = "black",
    size = 10/.pt,
    fontface = "italic"
  ) +
  theme_void() +
  coord_equal() +
  expand_limits(x = c(-6, 6), y = c(-6, 6)) +
  theme(legend.position = "none")


library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

target_genes <- c(
  "KLF1", "NFE2", "FOXO3",
  "HBA1", "HBA2", "HBG1", "HBG2", "AHSP",
  "ALAS2", "FECH", "HMBS", "UROS", "PPOX", "ABCB6", "SLC25A37", "BLVRB",
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN", "MPP1",
  "EPOR", "BNIP3L"
)

edges <- data.frame(
  from = "CEBPB",
  to   = target_genes
)

n      <- length(target_genes)
radius <- 5.5
angles <- pi/2 - (0:(n - 1)) * 2 * pi / n

nodes <- data.frame(
  name = c("CEBPB", target_genes),
  x    = c(0, radius * cos(angles)),
  y    = c(0, radius * sin(angles))
)

nodes <- nodes %>%
  mutate(
    fill_color = case_when(
      name == "CEBPB" ~ "#FDCC8A",
      name == "GDF15" ~ "#4393C3",
      TRUE            ~ "#D9D9D9"
    ),
    border_color = case_when(
      name == "CEBPB" ~ "#E08214",
      name == "GDF15" ~ "#2166AC",
      TRUE            ~ "#969696"
    ),
    angle   = atan2(y, x),
    label_x = ifelse(name == "CEBPB", 0, x + 1.3 * cos(angle)),
    label_y = ifelse(name == "CEBPB", 0, y + 1.3 * sin(angle)),
    hjust = ifelse(name == "CEBPB", 0.5,
                   ifelse(cos(angle) < -0.1, 1,
                          ifelse(cos(angle) > 0.1, 0, 0.5)))
  )

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +
  geom_edge_arc(
    strength = 0.15,
    color    = "grey75",
    width    = 0.5
  ) +
  geom_node_point(
    aes(x = x, y = y,
        fill  = I(fill_color),
        color = I(border_color)),
    shape  = 21,
    size   = ifelse(nodes$name == "CEBPB", 6, 3.5),
    stroke = 1.2
  ) +
  geom_node_text(
    aes(x = label_x, y = label_y, label = name, hjust = hjust),
    color    = "black",
    size     = 8 / .pt,
    fontface = "italic"
  ) +
  theme_void() +
  coord_equal() +
  expand_limits(x = c(-9, 9), y = c(-9, 9)) +
  theme(legend.position = "none")


library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

target_genes <- c(
  # Erythroid differentiation
  "KLF1", "EPOR", "NFE2", "FOXO3",
  # Heme / hemoglobin synthesis
  "AHSP", "ALAS2", "FECH", "HMBS", "UROS", "PPOX",
  "ABCB6", "SLC25A37", "HBA1", "HBA2", "HBG1", "HBG2",
  # Membrane & cytoskeleton
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN", "MPP1",
  # Terminal maturation & homeostasis
  "BNIP3L", "BLVRB",
  # GDF15
  "GDF15"
)

erythroid_diff <- c("KLF1", "EPOR", "NFE2", "FOXO3")
heme_hb        <- c("AHSP", "ALAS2", "FECH", "HMBS", "UROS", "PPOX",
                    "ABCB6", "SLC25A37", "HBA1", "HBA2", "HBG1", "HBG2")
membrane       <- c("GYPA", "GYPB", "GYPE", "GYPC", "DMTN", "MPP1")
terminal_mat   <- c("BNIP3L", "BLVRB")

edges <- data.frame(
  from = "CEBPB",
  to   = target_genes
)

n      <- length(target_genes)
radius <- 5.5
angles <- pi/2 - (0:(n - 1)) * 2 * pi / n

nodes <- data.frame(
  name = c("CEBPB", target_genes),
  x    = c(0, radius * cos(angles)),
  y    = c(0, radius * sin(angles))
)

nodes <- nodes %>%
  mutate(
    category = case_when(
      name == "CEBPB"          ~ "CEBPB",
      name == "GDF15"          ~ "GDF15",
      name %in% erythroid_diff ~ "Erythroid differentiation",
      name %in% heme_hb        ~ "Heme/hemoglobin synthesis",
      name %in% membrane       ~ "Membrane & cytoskeleton",
      name %in% terminal_mat   ~ "Terminal maturation & homeostasis",
      TRUE                     ~ "Other"
    ),
    angle   = atan2(y, x),
    label_x = ifelse(name == "CEBPB", 0, x + 1.3 * cos(angle)),
    label_y = ifelse(name == "CEBPB", 0, y + 1.3 * sin(angle)),
    hjust   = ifelse(name == "CEBPB", 0.5,
                     ifelse(cos(angle) < -0.1, 1,
                            ifelse(cos(angle) > 0.1, 0, 0.5)))
  )

cat_levels <- c(
  "CEBPB",
  "Erythroid differentiation",
  "Heme/hemoglobin synthesis",
  "Membrane & cytoskeleton",
  "Terminal maturation & homeostasis",
  "GDF15"
)
nodes$category <- factor(nodes$category, levels = cat_levels)

cat_colors <- c(
  "CEBPB"                     = "#FDCC8A",
  "Erythroid differentiation"           = "#E41A1C",
  "Heme/hemoglobin synthesis"           = "#66C2A5",
  "Membrane & cytoskeleton"             = "#8DA0CB",
  "Terminal maturation & homeostasis"    = "#FC8D62",
  "GDF15"                                = "#4393C3"
)

cat_borders <- c(
  "CEBPB"                     = "#E08214",
  "Erythroid differentiation"           = "#B2182B",
  "Heme/hemoglobin synthesis"           = "#2CA25F",
  "Membrane & cytoskeleton"             = "#5E81AC",
  "Terminal maturation & homeostasis"    = "#D95F02",
  "GDF15"                               = "#2166AC"
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +
  geom_edge_arc(
    strength = 0.15,
    color    = "grey78",
    width    = 0.4
  ) +
  geom_node_point(
    aes(fill = category, color = category),
    shape  = 21,
    size   = ifelse(V(g)$name == "CEBPB", 7, 4),
    stroke = 1.2
  ) +
  scale_fill_manual(values = cat_colors, name = "Functional category") +
  scale_color_manual(values = cat_borders, name = "Functional category") +
  geom_node_text(
    aes(x = label_x, y = label_y, label = name, hjust = hjust),
    color    = "black",
    size     = 8 / .pt,
    fontface = "italic"
  ) +
  theme_void() +
  coord_equal() +
  expand_limits(x = c(-9, 9), y = c(-9, 9)) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 8.5),
    legend.key.size  = unit(0.5, "cm"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(
    fill  = guide_legend(override.aes = list(size = 4)),
    color = guide_legend(override.aes = list(size = 4))
  )


library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(grid)

# ===============================
# Selected CEBPB regulon targets
# pySCENIC-inferred, erythroid-related subset
# ===============================

target_genes <- c(
  # Erythroid differentiation
  "KLF1", "EPOR", "NFE2", "FOXO3",

  # Heme / hemoglobin synthesis
  "ABCB6", "ALAS2", "HMBS", "PPOX",
  "HBA1", "HBA2", "HBG1", "HBG2",

  # Membrane & cytoskeleton
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN",

  # Terminal maturation / stress response
  "BNIP3L", "BLVRB", "GDF15"
)

erythroid_diff <- c("KLF1", "EPOR", "NFE2", "FOXO3")

terminal_stress <- c("BNIP3L", "BLVRB", "GDF15")

heme_hb <- c(
  "ABCB6", "ALAS2", "HMBS", "PPOX",
  "HBA1", "HBA2", "HBG1", "HBG2"
)

membrane <- c(
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN"
)


edges <- data.frame(
  from = "CEBPB",
  to   = target_genes
)

n <- length(target_genes)
radius <- 5.5
angles <- pi / 2 - (0:(n - 1)) * 2 * pi / n

nodes <- data.frame(
  name = c("CEBPB", target_genes),
  x    = c(0, radius * cos(angles)),
  y    = c(0, radius * sin(angles))
)

nodes <- nodes %>%
  mutate(
    annotation = case_when(
      name == "CEBPB" ~ "CEBPB",
      name %in% erythroid_diff ~ "Erythroid differentiation",
      name %in% heme_hb ~ "Heme/hemoglobin synthesis",
      name %in% membrane ~ "Membrane & cytoskeleton",
      name %in% terminal_stress ~ "Terminal maturation & stress response",
      TRUE ~ "Other"
    ),
    is_tf = ifelse(name == "CEBPB", "TF", "Target"),
    angle = atan2(y, x),
    label_x = ifelse(name == "CEBPB", 0, x + 1.15 * cos(angle)),
    label_y = ifelse(name == "CEBPB", 0, y + 1.15 * sin(angle)),
    hjust = ifelse(
      name == "CEBPB", 0.5,
      ifelse(cos(angle) < -0.1, 1,
             ifelse(cos(angle) > 0.1, 0, 0.5))
    ),
    label_face = ifelse(name == "CEBPB", "bold.italic", "italic")
  )

annotation_levels <- c(
  "CEBPB",
  "Erythroid differentiation",
  "Heme/hemoglobin synthesis",
  "Membrane & cytoskeleton",
  "Terminal maturation & stress response"
)

nodes$annotation <- factor(nodes$annotation, levels = annotation_levels)

anno_fill <- c(
  "CEBPB"                                  = "#FDCC8A",
  "Erythroid differentiation"             = "#E41A1C",
  "Heme/hemoglobin synthesis"             = "#66C2A5",
  "Membrane & cytoskeleton"               = "#8DA0CB",
  "Terminal maturation & stress response" = "#FC8D62"
)

anno_border <- c(
  "CEBPB"                                  = "#E08214",
  "Erythroid differentiation"             = "#B2182B",
  "Heme/hemoglobin synthesis"             = "#2CA25F",
  "Membrane & cytoskeleton"               = "#5E81AC",
  "Terminal maturation & stress response" = "#D95F02"
)

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

p <- ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +

  geom_edge_arc(
    strength = 0.12,
    color = "grey78",
    width = 0.35,
    alpha = 0.8
  ) +

  geom_node_point(
    aes(fill = annotation, color = annotation, size = is_tf),
    shape = 21,
    stroke = 1.1
  ) +

  scale_size_manual(
    values = c("TF" = 7.2, "Target" = 4.2),
    guide = "none"
  ) +

  scale_fill_manual(
    values = anno_fill,
    name = "Functional annotation"
  ) +

  scale_color_manual(
    values = anno_border,
    name = "Functional annotation"
  ) +

  geom_node_text(
    aes(
      x = label_x,
      y = label_y,
      label = name,
      hjust = hjust,
      fontface = label_face
    ),
    color = "black",
    size = 8 / .pt
  ) +

  labs(
    title = "Selected CEBPB Regulon Targets",
    caption = "Selected erythroid-related target genes from the pySCENIC-inferred CEBPB regulon; this subset does not represent the complete CEBPB regulon."
  ) +

  theme_void() +
  coord_equal() +
  expand_limits(x = c(-9, 9), y = c(-8.5, 8.5)) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12,
      color = "black",
      margin = margin(b = 5)
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 7.5,
      color = "grey30",
      margin = margin(t = 5)
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 10,
      color = "black"
    ),
    legend.text = element_text(
      size = 8.5,
      color = "black"
    ),
    legend.key.size = unit(0.45, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +

  guides(
    fill = guide_legend(
      override.aes = list(size = 4.5, stroke = 1)
    ),
    color = guide_legend(
      override.aes = list(size = 4.5, stroke = 1)
    )
  )

p

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(grid)

# ===============================
# Selected CEBPB regulon targets
# ===============================

target_genes <- c(
  "KLF1", "EPOR", "NFE2", "FOXO3",
  "ABCB6", "ALAS2", "HMBS", "PPOX",
  "HBA1", "HBA2", "HBG1", "HBG2",
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN",
  "BNIP3L", "BLVRB", "GDF15"
)

edges <- data.frame(
  from = "CEBPB",
  to   = target_genes
)

n <- length(target_genes)
radius <- 5.5
angles <- pi / 2 - (0:(n - 1)) * 2 * pi / n

nodes <- data.frame(
  name = c("CEBPB", target_genes),
  x    = c(0, radius * cos(angles)),
  y    = c(0, radius * sin(angles))
)

nodes <- nodes %>%
  mutate(
    is_tf = ifelse(name == "CEBPB", "TF", "Target"),
    angle = atan2(y, x),
    label_x = ifelse(name == "CEBPB", -0.25, x + 0.65 * cos(angle)),
    label_y = ifelse(name == "CEBPB", -0.45, y + 0.65 * sin(angle)),
    hjust = ifelse(
      name == "CEBPB", 0.5,
      ifelse(cos(angle) < -0.1, 1,
             ifelse(cos(angle) > 0.1, 0, 0.5))
    )
  )

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

p <- ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +

  geom_edge_arc(
    strength = 0.18,
    color = "grey65",
    width = 0.6,
    alpha = 0.75
  ) +

  geom_node_point(
    data = function(x) dplyr::filter(x, is_tf == "Target"),
    shape = 21,
    fill = "#D9DEE3",
    color = "#6F777D",
    size = 4.0,
    stroke = 0.8
  ) +

  geom_node_point(
    data = function(x) dplyr::filter(x, is_tf == "TF"),
    shape = 21,
    fill = "#5B5CE2",
    color = "#3B3E8C",
    size = 6.2,
    stroke = 1.0
  ) +

  geom_node_text(
    data = function(x) dplyr::filter(x, is_tf == "TF"),
    aes(
      x = label_x,
      y = label_y,
      label = name
    ),
    color = "blue",
    fontface = "bold.italic",
    size = 7
  ) +

  labs(
    title = expression(italic(CEBPB)~"gene module")
  ) +

  theme_void() +
  coord_equal() +
  expand_limits(x = c(-7, 7), y = c(-7, 7)) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "plain",
      size = 16,
      color = "black",
      margin = margin(b = 2)
    ),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

p


library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(grid)

# ===============================
# Selected CEBPB regulon targets
# pySCENIC-inferred, erythroid-related subset
# ===============================

target_genes <- c(
  # Erythroid differentiation
  "KLF1", "EPOR", "NFE2", "FOXO3",

  # Heme / hemoglobin synthesis
  "ABCB6", "ALAS2", "HMBS", "PPOX",
  "HBA1", "HBA2", "HBG1", "HBG2",

  # Membrane & cytoskeleton
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN",

  # Terminal maturation / stress response
  "BNIP3L", "BLVRB", "GDF15"
)

erythroid_diff <- c("KLF1", "EPOR", "NFE2", "FOXO3")

terminal_stress <- c("BNIP3L", "BLVRB", "GDF15")

heme_hb <- c(
  "ABCB6", "ALAS2", "HMBS", "PPOX",
  "HBA1", "HBA2", "HBG1", "HBG2"
)

membrane <- c(
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN"
)

edges <- data.frame(
  from = "CEBPB",
  to   = target_genes
)

n <- length(target_genes)
radius <- 5.5
angles <- pi / 2 - (0:(n - 1)) * 2 * pi / n

nodes <- data.frame(
  name = c("CEBPB", target_genes),
  x    = c(0, radius * cos(angles)),
  y    = c(0, radius * sin(angles))
)

nodes <- nodes %>%
  mutate(
    annotation = case_when(
      name == "CEBPB" ~ "CEBPB",
      name %in% erythroid_diff ~ "Erythroid differentiation",
      name %in% heme_hb ~ "Heme/hemoglobin synthesis",
      name %in% membrane ~ "Membrane & cytoskeleton",
      name %in% terminal_stress ~ "Terminal maturation & stress response",
      TRUE ~ "Other"
    ),
    is_tf = ifelse(name == "CEBPB", "TF", "Target"),
    angle = atan2(y, x),
    label_x = ifelse(name == "CEBPB", 0, x + 1.15 * cos(angle)),
    label_y = ifelse(name == "CEBPB", 0, y + 1.15 * sin(angle)),
    hjust = ifelse(
      name == "CEBPB", 0.5,
      ifelse(cos(angle) < -0.1, 1,
             ifelse(cos(angle) > 0.1, 0, 0.5))
    ),
    label_face = ifelse(name == "CEBPB", "bold.italic", "italic")
  )

annotation_levels <- c(
  "CEBPB",
  "Erythroid differentiation",
  "Heme/hemoglobin synthesis",
  "Membrane & cytoskeleton",
  "Terminal maturation & stress response"
)

nodes$annotation <- factor(nodes$annotation, levels = annotation_levels)

# ============================================================
# ============================================================

anno_fill <- c(
  "CEBPB"                                  = "#F3C677",
  "Erythroid differentiation"             = "#D95F5F",
  "Heme/hemoglobin synthesis"             = "#73BFA3",
  "Membrane & cytoskeleton"               = "#8FA8C8",
  "Terminal maturation & stress response" = "#E89A6A"
)

anno_border <- c(
  "CEBPB"                                  = "#C67C1A",
  "Erythroid differentiation"             = "#B73535",
  "Heme/hemoglobin synthesis"             = "#2F8F73",
  "Membrane & cytoskeleton"               = "#5C789E",
  "Terminal maturation & stress response" = "#C66A36"
)

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

p <- ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +

  geom_edge_arc(
    strength = 0.12,
    color = "grey70",
    width = 0.75,
    alpha = 0.9
  ) +

  geom_node_point(
    aes(fill = annotation, color = annotation, size = is_tf),
    shape = 21,
    stroke = 1.2
  ) +

  scale_size_manual(
    values = c("TF" = 7.4, "Target" = 4.4),
    guide = "none"
  ) +

  scale_fill_manual(
    values = anno_fill,
    guide = "none"
  ) +

  scale_color_manual(
    values = anno_border,
    guide = "none"
  ) +

  geom_node_text(
    aes(
      x = label_x,
      y = label_y,
      label = name,
      hjust = hjust,
      fontface = label_face
    ),
    color = "black",
    size = 8 / .pt
  ) +

  labs(
    title = "Selected CEBPB Regulon Targets",
    caption = "Selected erythroid-related target genes from the pySCENIC-inferred CEBPB regulon."
  ) +

  theme_void() +
  coord_equal() +
  expand_limits(x = c(-9, 9), y = c(-8.5, 8.5)) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12,
      color = "black",
      margin = margin(b = 5)
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 7.5,
      color = "grey30",
      margin = margin(t = 5)
    ),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

p

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(grid)

# ===============================
# Selected CEBPB regulon targets
# pySCENIC-inferred, erythroid-related subset
# ===============================

target_genes <- c(
  "KLF1", "EPOR", "NFE2", "FOXO3",
  "ABCB6", "ALAS2", "HMBS", "PPOX",
  "HBA1", "HBA2", "HBG1", "HBG2",
  "GYPA", "GYPB", "GYPE", "GYPC", "DMTN",
  "BNIP3L", "BLVRB", "GDF15"
)

edges <- data.frame(
  from = "CEBPB",
  to   = target_genes
)

n <- length(target_genes)
radius <- 5.5
angles <- pi / 2 - (0:(n - 1)) * 2 * pi / n

nodes <- data.frame(
  name = c("CEBPB", target_genes),
  x    = c(0, radius * cos(angles)),
  y    = c(0, radius * sin(angles))
)

nodes <- nodes %>%
  mutate(
    node_type = ifelse(name == "CEBPB", "TF", "Target"),
    angle = atan2(y, x),

    label_x = ifelse(name == "CEBPB", -0.15, x + 1.15 * cos(angle)),
    label_y = ifelse(name == "CEBPB", -0.35, y + 1.15 * sin(angle)),

    hjust = ifelse(
      name == "CEBPB", 0.5,
      ifelse(cos(angle) < -0.1, 1,
             ifelse(cos(angle) > 0.1, 0, 0.5))
    ),

    label_face = ifelse(name == "CEBPB", "bold.italic", "italic"),
    label_color = ifelse(name == "CEBPB", "#1F2AFF", "black")
  )

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = TRUE
)

# ===============================
# ===============================

p <- ggraph(g, layout = "manual", x = nodes$x, y = nodes$y) +

  geom_edge_arc(
    strength = 0.13,
    color = "grey65",
    width = 0.85,
    alpha = 0.85
  ) +

  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "Target"),
    shape = 21,
    fill = "#D9DEE3",
    color = "#6F777D",
    size = 4.3,
    stroke = 1.1
  ) +

  geom_node_point(
    data = function(x) dplyr::filter(x, node_type == "TF"),
    shape = 21,
    fill = "#5A5CE2",
    color = "#3D3F9B",
    size = 7.2,
    stroke = 1.1
  ) +

  geom_node_text(
    data = function(x) dplyr::filter(x, node_type == "Target"),
    aes(
      x = label_x,
      y = label_y,
      label = name,
      hjust = hjust
    ),
    color = "black",
    fontface = "italic",
    size = 8 / .pt,
    show.legend = FALSE
  ) +

  geom_node_text(
    data = function(x) dplyr::filter(x, node_type == "TF"),
    aes(
      x = label_x,
      y = label_y,
      label = name
    ),
    color = "#1F2AFF",
    fontface = "bold.italic",
    size = 8 / .pt,
    show.legend = FALSE
  ) +

  labs(
    title = expression(italic(CEBPB)~"gene module"),
    caption = "Selected erythroid-related target genes from the pySCENIC-inferred CEBPB regulon."
  ) +

  theme_void() +
  coord_equal() +
  expand_limits(x = c(-9, 9), y = c(-8.5, 8.5)) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "plain",
      size = 16,
      color = "black",
      margin = margin(b = 4)
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 7.5,
      color = "grey35",
      margin = margin(t = 5)
    ),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

p

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "data", "processed", "session_info_pyscenic_downstream.txt")
)

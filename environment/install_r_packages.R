#!/usr/bin/env Rscript

# Install packages required by the R scripts and notebooks.
#
# The original notebooks record R 4.2.3 but do not contain package versions.
# Run this installer inside a new renv project, execute the full workflow, and
# commit the resulting renv.lock only after validating the generated figures.

cran_packages <- c(
  "Cairo",
  "data.table",
  "devtools",
  "dplyr",
  "forcats",
  "ggraph",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggsci",
  "igraph",
  "IRkernel",
  "openxlsx",
  "patchwork",
  "pheatmap",
  "plotly",
  "RColorBrewer",
  "remotes",
  "renv",
  "reshape2",
  "reticulate",
  "rstatix",
  "scales",
  "Seurat",
  "SeuratObject",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse",
  "writexl"
)

missing_cran_packages <- setdiff(
  cran_packages,
  rownames(installed.packages())
)
if (length(missing_cran_packages) > 0) {
  install.packages(missing_cran_packages, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioconductor_packages <- c(
  "AUCell",
  "BiocParallel",
  "celda",
  "clusterProfiler",
  "ComplexHeatmap",
  "GenomeInfoDb",
  "GenomeInfoDbData",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "SCENIC",
  "SingleCellExperiment",
  "slingshot",
  "tradeSeq"
)

missing_bioconductor_packages <- setdiff(
  bioconductor_packages,
  rownames(installed.packages())
)
if (length(missing_bioconductor_packages) > 0) {
  BiocManager::install(missing_bioconductor_packages, ask = FALSE, update = FALSE)
}

github_packages <- c(
  "aertslab/SCopeLoomR",
  "chris-mcginnis-ucsf/DoubletFinder",
  "JinmiaoChenLab/CellChat",
  "cailab-tamu/scTenifoldKnk",
  "mengxu98/scop",
  "mengxu98/thisplot",
  "samuel-marsh/scCustomize@develop"
)

for (github_package in github_packages) {
  remotes::install_github(github_package, upgrade = "never")
}

message("Package installation completed.")
message("Next: validate the workflow, then run renv::snapshot().")

# Bone Marrow Organoid Single-Cell RNA-seq Analysis

## Overview

This repository contains the analysis code associated with a single-cell RNA-sequencing study of human bone marrow organoids (BMOs) cultured under static and dynamic conditions.

The workflow includes sample-level preprocessing, quality control, doublet detection, ambient RNA correction, Seurat-based integration and cell-type annotation, reference mapping, cell-type composition analysis, erythroid maturation analysis, cell–cell communication analysis, gene regulatory network inference, and in silico knockout analysis.

The scripts are organized according to the analysis workflow. Key analytical parameters, including quality-control thresholds, principal components, clustering resolutions, cell-type annotation markers, cluster-to-cell-type mappings, and selected target genes, are recorded directly in the corresponding scripts to facilitate review and reproducibility.

## Repository structure

```text
.
├── analysis/
│   ├── 01_preprocessing/             # Quality control, DoubletFinder, and DecontX
│   ├── 02_annotation/                # BMO and egressed blood cell annotation
│   ├── 03_reference_integration/     # Integration with bone marrow and organoid references
│   ├── 04_projection_mapping/        # Reference projection and label transfer
│   ├── 05_composition/               # Cell-type composition analysis
│   ├── 06_erythroid/                 # Erythroid differential analyses
│   ├── 07_cell_communication/        # CellChat analysis
│   ├── 08_gene_regulatory_networks/  # pySCENIC and CellOracle analyses
│   └── 09_virtual_knockout/          # scTenifoldKnk analysis
└── environment/                      # R and Python environment information

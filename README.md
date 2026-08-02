# Bone marrow organoid single-cell RNA-seq analysis

## Overview

This repository contains the analysis code associated with a single-cell RNA-sequencing study of human bone marrow organoids (BMOs) cultured under static and dynamic conditions.

The workflow includes sample-level quality control, doublet detection, ambient RNA correction, Seurat-based data integration and cell-type annotation, reference integration and mapping, cell-type composition analysis, erythroid-focused analysis, cell–cell communication analysis, gene regulatory network inference, and in silico knockout analysis.

The analysis is organized into numbered directories following the main computational workflow. Key parameters, including quality-control thresholds, principal components, clustering resolutions, cell-type markers, cluster-to-cell-type mappings, comparison groups, and selected target genes, are recorded directly in the corresponding scripts.

## Repository structure

```text
.
├── analysis/
│   ├── 01_preprocessing/               # Quality control, DoubletFinder, and DecontX
│   ├── 02_annotation/                  # BMO and egressed hematopoietic cell annotation
│   ├── 03_reference_integration/       # Integration with bone marrow and organoid references
│   ├── 04_projection_mapping/          # Reference projection and label transfer
│   ├── 05_composition/                 # Cell-type composition analysis
│   ├── 06_erythroid/                   # Erythroid subclustering, differential expression, and ABM mapping
│   ├── 07_cell_communication/          # CellChat analysis
│   ├── 08_gene_regulatory_networks/    # pySCENIC and CellOracle analyses
│   └── 09_virtual_knockout/            # scTenifoldKnk analysis
├── environment/
│   ├── README.md                       # Environment setup instructions
│   ├── install_r_packages.R            # R package installation script
│   └── requirements.txt                # Python package requirements
└── README.md
```

A detailed description of the purpose and execution order of each analysis stage is provided in:

```text
analysis/README.md
```

## Requirements

The analyses were performed using:

- R 4.2.3
- Python 3.10.18

Core software and packages include:

- Seurat
- DoubletFinder
- DecontX and celda
- Slingshot
- CellChat
- pySCENIC
- CellOracle
- scTenifoldKnk

Instructions for setting up the R and Python environments are provided in:

```text
environment/README.md
```

R packages can be installed using:

```bash
Rscript environment/install_r_packages.R
```

Python packages can be installed using:

```bash
pip install -r environment/requirements.txt
```

Some analyses, particularly multi-dataset integration, cell–cell communication analysis, and gene regulatory network inference, may require substantial memory and runtime.

## Data requirements

Raw sequencing data, processed analysis objects, and large external reference datasets are not stored in this repository.

The workflow requires:

- raw 10x Genomics count matrices for the BMO samples;
- raw 10x Genomics count matrices for the egressed hematopoietic cell samples;
- adult bone marrow reference data;
- fetal bone marrow reference data;
- published bone marrow organoid reference data; and
- motif annotation and ranking databases required for pySCENIC.

The required input objects and local file paths are described in the corresponding analysis scripts.

Users should update the input and output paths according to their local computing environment before running the analyses.

## Analysis workflow

The numbered directories under `analysis/` follow the main order of the computational workflow.

### 1. Preprocessing

Scripts in:

```text
analysis/01_preprocessing/
```

perform sample-level quality control, doublet detection using DoubletFinder, and ambient RNA correction using DecontX.

BMO samples and egressed hematopoietic cell samples are processed separately before downstream integration and annotation.

### 2. Cell-type annotation

Scripts in:

```text
analysis/02_annotation/
```

integrate the processed samples and perform dimensionality reduction, clustering, and cell-type annotation.

Cell identities were assigned manually at the cluster level according to canonical marker expression. The marker sets and cluster-to-cell-type mappings used in the study are explicitly recorded in the corresponding scripts.

### 3. Reference integration

Scripts and notebooks in:

```text
analysis/03_reference_integration/
```

integrate BMO populations with adult bone marrow, fetal bone marrow, and published bone marrow organoid reference datasets.

These analyses are used to compare the transcriptional features of BMO-derived stromal and hematopoietic populations with native human bone marrow and previously reported organoid models.

### 4. Projection and mapping

Scripts and notebooks in:

```text
analysis/04_projection_mapping/
```

perform reference projection, label transfer, and mapping of BMO-derived cells to corresponding native bone marrow populations.

### 5. Cell-type composition

Scripts and notebooks in:

```text
analysis/05_composition/
```

compare the relative abundance of stromal and hematopoietic populations across experimental groups and reference datasets.

### 6. Erythroid analysis

The `06_erythroid/` directory contains three erythroid-focused analyses:

```text
analysis/06_erythroid/
├── 01_erythroid_subclustering.R
├── 02_ABM_vs_DBMO_differential_expression.ipynb
└── 03_map_DBMO_erythroid_to_ABM.ipynb
```

These files should be executed in numerical order.

#### 6.1 Erythroid subclustering

```text
01_erythroid_subclustering.R
```

extracts erythroid cells, performs erythroid-specific dimensionality reduction and clustering, and prepares the erythroid populations for downstream comparisons.

#### 6.2 Differential expression between ABM and DBMO erythroid cells

```text
02_ABM_vs_DBMO_differential_expression.ipynb
```

compares erythroid cells derived from adult bone marrow (ABM) and dynamic bone marrow organoids (DBMOs) and identifies differentially expressed genes between the two groups.

#### 6.3 Mapping DBMO erythroid cells to ABM

```text
03_map_DBMO_erythroid_to_ABM.ipynb
```

maps DBMO erythroid cells to erythroid populations in the ABM reference dataset to evaluate their transcriptional similarity and maturation-state distribution.

### 7. Cell–cell communication

Scripts in:

```text
analysis/07_cell_communication/
```

use CellChat to compare predicted signaling interactions between static BMOs and DBMOs, with particular emphasis on communication between endothelial and erythroid populations.

### 8. Gene regulatory network analysis

Scripts and notebooks in:

```text
analysis/08_gene_regulatory_networks/
```

perform gene regulatory network inference and downstream analysis using pySCENIC and CellOracle.

The pySCENIC workflow is used to identify transcription factor regulons and their predicted target genes. CellOracle is used for regulatory network analysis and computational perturbation of selected transcription factors.

Python notebooks should be executed using the validated Python environment after the required input files have been prepared.

### 9. Virtual knockout analysis

Scripts in:

```text
analysis/09_virtual_knockout/
```

perform in silico gene perturbation and virtual knockout analysis using scTenifoldKnk.

## Reproducibility

To support review and reproducibility:

- sample-specific quality-control parameters are retained in the preprocessing scripts;
- dimensionality reduction and clustering parameters are explicitly recorded;
- cell-type marker sets and manual cluster mappings are documented;
- comparison groups and selected genes are defined directly in the relevant scripts;
- fixed random seeds are used where applicable;
- absolute local paths and usernames have been removed where possible; and
- analysis files are organized according to their intended execution order.

Complete reproduction requires access to the study data and the external reference datasets used in the original analyses.

Software packages and reference resources may change over time. Users should therefore verify package compatibility and resource versions before rerunning the complete workflow.

## Data availability

The single-cell RNA-sequencing data generated in this study have been deposited in the Gene Expression Omnibus and will be made publicly available upon publication of the associated article.

External reference datasets should be obtained from their original repositories and used in accordance with the corresponding data-access and usage requirements.

Large intermediate analysis objects and generated results are not hosted in this GitHub repository.

## Citation

Please cite the associated manuscript when using or adapting the analysis code provided in this repository.

The complete manuscript citation and DOI will be added upon publication.

# Bone Marrow Organoid Single-Cell RNA-seq Analysis

## Overview

This repository contains the analysis code associated with a single-cell RNA-sequencing study of human bone marrow organoids (BMOs) cultured under static and dynamic conditions.

The workflow includes sample-level quality control, doublet detection, ambient RNA correction, Seurat-based data integration and cell-type annotation, reference integration and mapping, cell-type composition analysis, erythroid maturation analysis, cell–cell communication analysis, gene regulatory network inference, and in silico knockout analysis.

The analysis is organized into numbered directories following the main computational workflow. Key parameters, including quality-control thresholds, principal components, clustering resolutions, cell-type markers, cluster annotations, and selected target genes, are recorded directly in the corresponding scripts.

## Repository structure

```text
.
├── analysis/
│   ├── 01_preprocessing/               # Quality control, DoubletFinder, and DecontX
│   ├── 02_annotation/                  # BMO and egressed hematopoietic cell annotation
│   ├── 03_reference_integration/       # Integration with bone marrow and organoid references
│   ├── 04_projection_mapping/          # Reference projection and label transfer
│   ├── 05_composition/                 # Cell-type composition analysis
│   ├── 06_erythroid/                   # Erythroid maturation and trajectory analyses
│   ├── 07_cell_communication/          # CellChat analysis
│   ├── 08_gene_regulatory_networks/    # pySCENIC and CellOracle analyses
│   ├── 09_virtual_knockout/            # scTenifoldKnk analysis
│   └── README.md                       # Detailed description of the analysis workflow
├── environment/
│   ├── README.md                       # Environment setup instructions
│   ├── install_r_packages.R            # R package installation script
│   └── requirements.txt                # Python package requirements
└── README.md
```

A detailed description of the purpose, input, output, and execution order of each analysis stage is provided in:

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

Raw sequencing data, processed Seurat objects, and large external reference datasets are not stored in this repository.

The workflow requires:

- raw 10x Genomics count matrices for the BMO and egressed hematopoietic cell samples;
- adult bone marrow reference data;
- fetal bone marrow reference data;
- published bone marrow organoid reference data; and
- motif annotation and ranking databases required for pySCENIC.

The required input objects and file paths are described in the corresponding analysis scripts and in `analysis/README.md`.

Users should update local data paths where necessary before running the analysis.

## Running the analysis

Run the analyses from the repository root and follow the numbered directories in order.

### 1. Preprocessing

Scripts in:

```text
analysis/01_preprocessing/
```

perform sample-level quality control, doublet detection using DoubletFinder, and ambient RNA correction using DecontX.

Separate preprocessing is performed for BMO samples and egressed hematopoietic cell samples.

### 2. Cell-type annotation

Scripts in:

```text
analysis/02_annotation/
```

integrate the processed samples and perform clustering and cell-type annotation.

Cell identities were assigned manually at the cluster level according to canonical marker expression. The marker sets and cluster-to-cell-type mappings used in the study are recorded explicitly in the annotation scripts.

### 3. Reference integration

Notebooks and scripts in:

```text
analysis/03_reference_integration/
```

compare BMO populations with adult bone marrow, fetal bone marrow, and published bone marrow organoid reference datasets.

### 4. Projection and mapping

Analyses in:

```text
analysis/04_projection_mapping/
```

perform reference projection, label transfer, and mapping of BMO-derived cells to corresponding native bone marrow populations.

### 5. Cell-type composition

Analyses in:

```text
analysis/05_composition/
```

compare the relative abundance of stromal and hematopoietic populations across experimental groups and reference datasets.

### 6. Erythroid analysis

Analyses in:

```text
analysis/06_erythroid/
```

include erythroid subclustering, differential expression analysis, maturation-state comparison, cell-cycle analysis, trajectory inference, and comparison with native bone marrow erythroid populations.

### 7. Cell–cell communication

Scripts in:

```text
analysis/07_cell_communication/
```

use CellChat to compare predicted signaling interactions between static BMOs and dynamic BMOs, with particular emphasis on endothelial–erythroid communication.

### 8. Gene regulatory network analysis

Scripts and notebooks in:

```text
analysis/08_gene_regulatory_networks/
```

perform gene regulatory network inference and downstream analysis using pySCENIC and CellOracle.

The Python notebooks should be executed using the validated Python environment after the required input files have been prepared.

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
- selected genes and comparison groups are defined directly in the relevant scripts;
- fixed random seeds are used where applicable; and
- the analyses are organized according to their intended execution order.

Complete reproduction requires access to the study data and the external reference datasets used in the original analyses.

Because some software packages and reference resources may change over time, users should verify package compatibility before rerunning the complete workflow.

## Data availability

The single-cell RNA-sequencing data associated with this study have been deposited in the Gene Expression Omnibus under accession:

```text
GSE286476
```

External reference datasets should be obtained from their original repositories and used in accordance with the corresponding data-access and usage requirements.

Large intermediate objects and generated analysis results are not distributed through GitHub.

## Citation

Please cite the associated manuscript when using or adapting the analysis code provided in this repository.

The complete manuscript citation and DOI will be added upon publication.

## License

Licensing information will be provided upon public release.

# Analysis workflow

The repository follows the order in which data are transformed. Every numbered directory represents one analysis stage; generated objects are passed forward through explicit files rather than an interactive R session.

```text
Raw 10x matrices
  -> quality control, DoubletFinder, DecontX
  -> BMO and egress-blood integration
  -> cell-type annotation
  -> reference integration and projection
  -> composition and erythroid analyses
  -> CellChat, pySCENIC, CellOracle, and scTenifoldKnk
```

## Stage 1: preprocessing

`analysis/01_preprocessing/01_preprocess_scrna_samples.R` reads `config/sample_manifest.csv` and applies the same seven visible operations to each sample. It writes post-QC, post-doublet, and post-DecontX checkpoint objects plus diagnostic PDFs and `sessionInfo()` output.

## Stage 2: integration and annotation

- `analysis/02_annotation/01_annotate_bmo.R` integrates the seven BMO libraries, separates stromal and haematopoietic compartments, assigns cluster labels, and writes annotated objects.
- `analysis/02_annotation/02_annotate_egress_blood.R` integrates the four egress-blood libraries and assigns haematopoietic labels.

Cluster-to-cell-type mappings are kept directly in the scripts because they are biological decisions that reviewers must be able to inspect. They should be cross-checked against the final manuscript figures before release.

## Stages 3–5: external references and composition

The R notebooks integrate adult bone marrow, fetal bone marrow, BMO, and published organoid references; calculate correlations; transfer labels; and compare cell-type composition. Notebook output cells are cleared intentionally. Execute the notebooks in numerical order and save final figures to `results/figures/`.

## Stage 6: erythroid analyses

The erythroid directory contains subclustering, adult-bone-marrow comparison, label transfer, and Slingshot trajectory analysis. Inputs are the annotated objects from Stage 2 and the integrated reference objects from Stage 3.

## Stage 7: cell-cell communication

`analysis/07_cell_communication/01_cellchat_comparison.R` compares static and dynamic day-25 BMO samples using the human CellChat database, with a minimum of 10 cells per communicating group.

## Stage 8: regulatory-network analyses

1. `01_run_pyscenic.sh` runs GRN inference, motif pruning, and AUCell scoring.
2. `02_pyscenic_downstream.R` reads the loom result and performs regulon analyses.
3. `03_prepare_celloracle_input.R` exports the erythroid Seurat object to H5AD.
4. `04_celloracle_analysis.ipynb` constructs and perturbs the CellOracle network.

## Stage 9: virtual knockout

`analysis/09_virtual_knockout/01_sctenifoldknk.R` contains the scTenifoldKnk perturbation analyses. The target gene, selected variable genes, cell subsampling, network count, and significance threshold remain visible in the script.

## Reproducibility boundary

The source repository contains analysis logic but not raw or processed biological data. Exact reproduction additionally requires the data accessions, package lock files, reference-database releases, random seeds, and manuscript-to-code mapping listed in `docs/reproducibility_checklist.md`.

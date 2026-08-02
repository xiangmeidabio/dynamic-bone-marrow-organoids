# Bone marrow organoid single-cell RNA-seq analysis

This repository contains the analysis code for a bone marrow organoid (BMO) single-cell RNA-sequencing study. The workflow covers sample-level quality control, doublet and ambient-RNA removal, Seurat integration and annotation, reference mapping, cell-type composition, erythroid maturation, cell-cell communication, gene-regulatory networks, and virtual knockout analysis.

The repository is organized for transparent review and computational reproduction. Analysis decisions such as thresholds, principal-component ranges, clustering resolutions, cluster-to-cell-type mappings, and target genes remain visible in the numbered scripts. Small helper functions are used only where an algorithm must be applied identically to multiple samples.

## Repository structure

```text
.
├── analysis/
│   ├── 01_preprocessing/          # QC, DoubletFinder, and DecontX
│   ├── 02_annotation/             # BMO and egress-blood integration/annotation
│   ├── 03_reference_integration/  # ABM, FBM, BMO, and organoid references
│   ├── 04_projection_mapping/     # Reference projection and label transfer
│   ├── 05_composition/            # Cell-type proportion comparisons
│   ├── 06_erythroid/              # Erythroid DE and trajectory analyses
│   ├── 07_cell_communication/     # CellChat
│   ├── 08_gene_regulatory_networks/ # pySCENIC and CellOracle
│   └── 09_virtual_knockout/       # scTenifoldKnk
├── config/                        # Registered sample-specific parameters
├── data/                          # Local data layout; large files are Git-ignored
├── docs/                          # Workflow and release documentation
├── environment/                   # R/Python environment reconstruction
└── results/                       # Generated figures, tables, and objects
```

See `docs/analysis_workflow.md` for the scientific purpose and dependencies of every stage.

## Requirements

The original notebook metadata records R 4.2.3 and Python 3.10.18. Exact package versions were not present in the source files and therefore have not been invented. Instructions for reconstructing, validating, and locking both environments are provided in `environment/README.md`.

Core software includes Seurat v5, DoubletFinder, DecontX/celda, Slingshot, CellChat, pySCENIC, CellOracle, and scTenifoldKnk. Some analyses require substantial memory and runtime.

## Data setup

Raw and reference data are not committed to Git. Arrange the inputs as documented in `data/README.md`. At minimum:

1. place the 11 raw 10x Genomics matrices under `data/raw/`;
2. place adult bone marrow, fetal bone marrow, and published organoid reference objects under `data/reference/`; and
3. place pySCENIC ranking, motif, and transcription-factor resources under `data/reference/pyscenic/`.

The complete sample-to-directory mapping and all sample-specific preprocessing parameters are in `config/sample_manifest.csv`.

## Running the analysis

Run commands from the repository root. All scripts also accept the repository location through the `BMO_PROJECT_ROOT` environment variable.

### 1. Preprocess the BMO samples

```bash
Rscript analysis/01_preprocessing/01_preprocess_scrna_samples.R BMO
```

For each sample, this produces a post-QC object, a post-doublet object, a post-DecontX object, diagnostic figures, and an R session report.

### 2. Preprocess the egress-blood samples

```bash
Rscript analysis/01_preprocessing/01_preprocess_scrna_samples.R EGRESS_BLOOD
```

### 3. Integrate and annotate cells

```bash
Rscript analysis/02_annotation/01_annotate_bmo.R
Rscript analysis/02_annotation/02_annotate_egress_blood.R
```

Manual cluster mappings are explicit in these scripts. Review them against marker plots and the final manuscript before accepting the generated labels.

### 4. Run reference and downstream analyses

Open the notebooks in directories `03_reference_integration` through `06_erythroid` and execute them in numerical order. Notebook outputs are intentionally cleared from Git; figures and tables are written to `results/`.

Run the script-based analyses as needed:

```bash
Rscript analysis/06_erythroid/01_erythroid_subclustering.R
Rscript analysis/07_cell_communication/01_cellchat_comparison.R
bash analysis/08_gene_regulatory_networks/01_run_pyscenic.sh
Rscript analysis/08_gene_regulatory_networks/02_pyscenic_downstream.R
Rscript analysis/08_gene_regulatory_networks/03_prepare_celloracle_input.R
Rscript analysis/09_virtual_knockout/01_sctenifoldknk.R
```

Execute `analysis/08_gene_regulatory_networks/04_celloracle_analysis.ipynb` with the validated Python kernel after preparing the H5AD input.

## Outputs and provenance

- `data/processed/` contains intermediate and final analysis objects.
- `results/figures/` contains publication figures.
- `results/tables/` contains exported statistics and gene tables.
- `results/objects/` contains pySCENIC and CellOracle network objects.

These directories are excluded from Git by default. Deposit large reproducibility objects in an appropriate data repository and publish checksums and stable links.

## Reproducibility notes

- A fixed random seed is set in publication-facing scripts.
- Absolute local paths and usernames have been removed.
- Notebook outputs and execution counters are cleared.
- The preprocessing manifest preserves the thresholds from the original analysis.
- R and Python package locks must still be generated from the final validated environment.
- Full numerical verification requires the study data and reference objects, which are not included here.

The mandatory pre-release checks are listed in `docs/reproducibility_checklist.md`.

## Citation

Complete `CITATION.cff.template`, rename it to `CITATION.cff`, and update this section with the paper and archived software DOI. Cite an immutable release rather than the moving default branch.

## Data and code availability

A manuscript-ready template is provided in `docs/code_availability_statement.md`. Replace every bracketed field before submission.

## License

No software license has been selected because that choice requires agreement from the authors and institution. Add the approved license before public release; without a license, others may view the repository but do not automatically receive permission to reuse the code.

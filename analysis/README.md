# Analysis index

Run the numbered directories in order. Within a directory, run files in numerical order unless the table marks an analysis as optional.

| Stage | File | Main input | Main output |
|---|---|---|---|
| 01 | `01_preprocessing/01_preprocess_scrna_samples.R` | Raw 10x matrices and sample manifest | Per-sample QC, doublet-filtered, and DecontX RDS objects |
| 02 | `02_annotation/01_annotate_bmo.R` | BMO post-DecontX objects | Annotated stromal, haematopoietic, and combined BMO objects |
| 02 | `02_annotation/02_annotate_egress_blood.R` | Egress-blood post-doublet objects | Annotated egress-blood object |
| 03 | `03_reference_integration/*.ipynb` | Annotated BMO and external references | Integrated reference objects and correlation figures |
| 04 | `04_projection_mapping/*.ipynb` | Integrated references | Projection figures and scores |
| 05 | `05_composition/*.ipynb` | Integrated stromal/haematopoietic objects | Composition tables and figures |
| 06 | `06_erythroid/01_erythroid_subclustering.R` | Annotated BMO and egress-blood objects | Erythroid subset objects |
| 06 | Remaining erythroid notebooks | Erythroid objects and ABM reference | DE, label-transfer, and trajectory results |
| 07 | `07_cell_communication/01_cellchat_comparison.R` | Combined annotated BMO object | CellChat objects and comparison figures |
| 08 | `08_gene_regulatory_networks/*` | Erythroid object plus GRN databases | pySCENIC and CellOracle networks/figures |
| 09 | `09_virtual_knockout/01_sctenifoldknk.R` | Erythroid objects | Virtual-knockout tables, objects, and figures |

Every file resolves paths from the repository root. R notebooks use the `ir` kernel; the CellOracle notebook uses Python 3.10. Execute notebooks from top to bottom because later cells use objects created by earlier cells.

# Data organization

Large data files are intentionally excluded from Git. Deposit public data in an archival repository and record accession numbers in this file and in the manuscript Data Availability statement.

## Raw 10x Genomics matrices

Place each Cell Ranger filtered-feature matrix directory at the path shown below. Each directory must contain `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` (or the corresponding uncompressed files).

```text
data/raw/
├── bmo/
│   ├── Dynamic-24d/
│   ├── Dynamic-25d-1/
│   ├── Dynamic-25d-2/
│   ├── Dynamic-31d/
│   ├── Static-24d/
│   ├── Static-25d-1/
│   └── Static-25d-2/
└── egress_blood/
    ├── blood-yao-31/
    ├── blood-yao2-31/
    ├── blood-ctr-31/
    └── blood-ctr-33/
```

The mapping between directory names, biological sample IDs, and preprocessing parameters is recorded in `config/sample_manifest.csv`.

## Reference data

Place external reference objects under `data/reference/`:

```text
data/reference/
├── Adult_Bone_Marrow.RData
├── FBM_obj_2.RData
├── organoids23/
│   └── Organoids23.RData
└── pyscenic/
    ├── allTFs_hg38.txt
    ├── hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
    └── motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

For every external dataset, add its accession/DOI, license, download date, checksum, genome build, and any preprocessing performed before use.

## Processed data

Analysis scripts write derived objects below `data/processed/`. These files are ignored by Git because RDS, H5AD, and loom objects can be large and may contain controlled or unpublished data. If processed objects are required for reproduction, deposit them in Zenodo, Figshare, Dryad, GEO, or another appropriate repository and provide a manifest with SHA-256 checksums.

## Data manifest to complete before submission

| Dataset | Repository/accession | Genome build | Public release date | Notes |
|---|---|---|---|---|
| BMO scRNA-seq | TODO | GRCh38 (verify) | TODO | Seven organoid samples |
| Egress-blood scRNA-seq | TODO | GRCh38 (verify) | TODO | Four samples |
| Adult bone marrow reference | TODO | TODO | Public | Add source publication |
| Fetal bone marrow reference | TODO | TODO | Public | Add source publication |
| 2023 organoid reference | TODO | TODO | Public | Add source publication |

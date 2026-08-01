# Preprocessing

`01_preprocess_scrna_samples.R` replaces the two duplicated preprocessing scripts in the original analysis. The biological parameters were not averaged or guessed; each original threshold is recorded in `config/sample_manifest.csv`.

The script writes three checkpoint objects per sample:

1. `01_post_qc.rds`
2. `02_post_doublet_removal.rds`
3. `03_post_decontx.rds`

This checkpoint design makes every change to the expression matrix traceable. Diagnostic PDFs and a `session_info.txt` file are saved with the final objects.

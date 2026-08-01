# Configuration files

`sample_manifest.csv` records every sample-specific parameter used during preprocessing. Values were transcribed from the original scripts rather than silently harmonized. This preserves the analysis decisions made for each library and makes them auditable.

Input directories are interpreted relative to:

- `data/raw/bmo/` for `BMO` samples; and
- `data/raw/egress_blood/` for `EGRESS_BLOOD` samples.

Before rerunning the analysis, verify the QC thresholds and expected doublet rates against the final sequencing metrics reported in the manuscript.

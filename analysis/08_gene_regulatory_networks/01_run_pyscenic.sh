#!/usr/bin/env bash

# Run the three pySCENIC command-line stages in a reproducible order.
# Execute this script from the repository root after installing pySCENIC.

set -euo pipefail

PROJECT_ROOT="${BMO_PROJECT_ROOT:-$(pwd)}"
INPUT_LOOM="${PROJECT_ROOT}/data/processed/BMO_Erythroid_sce.loom"
TF_LIST="${PROJECT_ROOT}/data/reference/pyscenic/allTFs_hg38.txt"
MOTIF_ANNOTATIONS="${PROJECT_ROOT}/data/reference/pyscenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
RANKINGS_DATABASE="${PROJECT_ROOT}/data/reference/pyscenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
OUTPUT_DIRECTORY="${PROJECT_ROOT}/results/objects/grn"
WORKERS="${PYSCENIC_WORKERS:-8}"

mkdir -p "${OUTPUT_DIRECTORY}"

# Step 1: infer transcription-factor-to-target co-expression modules.
pyscenic grn \
  --num_workers "${WORKERS}" \
  --sparse \
  --method grnboost2 \
  --output "${OUTPUT_DIRECTORY}/grn.csv" \
  "${INPUT_LOOM}" \
  "${TF_LIST}"

# Step 2: prune modules using motif enrichment to define regulons.
pyscenic ctx \
  --num_workers "${WORKERS}" \
  --output "${OUTPUT_DIRECTORY}/regulons.csv" \
  --expression_mtx_fname "${INPUT_LOOM}" \
  --mode custom_multiprocessing \
  --annotations_fname "${MOTIF_ANNOTATIONS}" \
  "${OUTPUT_DIRECTORY}/grn.csv" \
  "${RANKINGS_DATABASE}"

# Step 3: calculate regulon activity in every cell.
pyscenic aucell \
  --num_workers "${WORKERS}" \
  --output "${OUTPUT_DIRECTORY}/BMO_Erythroid_SCENIC.loom" \
  "${INPUT_LOOM}" \
  "${OUTPUT_DIRECTORY}/regulons.csv"

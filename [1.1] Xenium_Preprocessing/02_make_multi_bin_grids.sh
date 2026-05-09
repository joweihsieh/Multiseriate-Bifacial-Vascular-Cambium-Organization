#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash 02_make_multi_bin_grids.sh <transcripts_xyz.csv> [neighbor_k]" >&2
  exit 1
fi

INPUT_CSV="$1"
NEIGHBOR_K="${2:-2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID_SCRIPT="${SCRIPT_DIR}/xenium_grid_um_from_csv.R"

#for BIN_UM in 1 2 3 4 5 6 7 8 9 10; do
for BIN_UM in 5; do
  OUTDIR="grid${BIN_UM}um_out"
  echo "Running ${BIN_UM} um -> ${OUTDIR}"
  Rscript "${GRID_SCRIPT}" \
    --input "${INPUT_CSV}" \
    --outdir "${OUTDIR}" \
    --bin_um "${BIN_UM}" \
    --neighbor_k "${NEIGHBOR_K}"
done

#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <seed> <loop> <base_beta_file> <out_csv> [J] [Time] [start_sd] [maxit]" >&2
  exit 2
fi

SEED="$1"
LOOP="$2"
BASE_BETA_FILE="$3"
OUT_CSV="$4"
J="${5:-100}"
TIME="${6:-3}"
START_SD="${7:-0.1}"
MAXIT="${8:-5000}"
START_ID="${LSB_JOBINDEX:-1}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${SCRIPT_DIR}/gmm_start_worker.R" \
  --seed "${SEED}" \
  --loop "${LOOP}" \
  --start-id "${START_ID}" \
  --base-beta-file "${BASE_BETA_FILE}" \
  --out-csv "${OUT_CSV}" \
  --J "${J}" \
  --Time "${TIME}" \
  --start-sd "${START_SD}" \
  --maxit "${MAXIT}"

#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <seed> <out_dir> [loops] [J] [Time] [start_sd] [maxit]" >&2
  exit 2
fi

SEED="$1"
OUT_DIR="$2"
LOOPS="${3:-6}"
J="${4:-100}"
TIME="${5:-3}"
START_SD="${6:-0.1}"
MAXIT="${7:-8000}"
START_ID="${LSB_JOBINDEX:-1}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORKER_ARGS=(
  --seed "${SEED}"
  --start-id "${START_ID}"
  --loops "${LOOPS}"
  --out-dir "${OUT_DIR}"
  --J "${J}"
  --Time "${TIME}"
  --start-sd "${START_SD}"
  --maxit "${MAXIT}"
)

if command -v "${RSCRIPT:-Rscript}" >/dev/null 2>&1; then
  exec "${RSCRIPT:-Rscript}" "${SCRIPT_DIR}/gmm_start_worker.R" "${WORKER_ARGS[@]}"
elif command -v "${R_BIN:-R}" >/dev/null 2>&1; then
  exec "${R_BIN:-R}" --vanilla --slave --file="${SCRIPT_DIR}/gmm_start_worker.R" --args "${WORKER_ARGS[@]}"
else
  echo "Neither Rscript nor R is available. Load an R module before submitting." >&2
  exit 127
fi

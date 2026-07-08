#!/usr/bin/env bash
set -euo pipefail

SEED="${SEED:-1}"
LOOPS="${LOOPS:-20}"
STARTS="${STARTS:-120}"
TARGET="${TARGET:-100}"
TOL="${TOL:-1e-7}"
OUT_DIR="${OUT_DIR:-hpc_runs/seed${SEED}}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if command -v Rscript >/dev/null 2>&1; then
  exec Rscript gmm_hpc_driver.R \
    --seed "${SEED}" \
    --loops "${LOOPS}" \
    --starts "${STARTS}" \
    --target "${TARGET}" \
    --tol "${TOL}" \
    --out-dir "${OUT_DIR}" \
    "$@"
fi

if command -v R >/dev/null 2>&1; then
  exec R --vanilla --slave --file=gmm_hpc_driver.R --args \
    --seed "${SEED}" \
    --loops "${LOOPS}" \
    --starts "${STARTS}" \
    --target "${TARGET}" \
    --tol "${TOL}" \
    --out-dir "${OUT_DIR}" \
    "$@"
fi

echo "Neither Rscript nor R is available. Load an R module first." >&2
exit 127

#!/usr/bin/env bash
set -euo pipefail

SEED="${SEED:-1}"
LOOPS="${LOOPS:-20}"
STARTS="${STARTS:-120}"
TARGET="${TARGET:-100}"
TOL="${TOL:-1e-7}"
OUT_DIR="${OUT_DIR:-hpc_runs/seed${SEED}}"
QUEUE="${QUEUE:-${QUEUE_NAME:-serial}}"
WALL_TIME="${WALL_TIME:-14:00}"
MEMORY_GB="${MEMORY_GB:-8}"
R_MODULE="${R_MODULE:-R/4.4.0}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

DRIVER_ARGS=(
  --seed "${SEED}"
  --loops "${LOOPS}"
  --starts "${STARTS}"
  --target "${TARGET}"
  --tol "${TOL}"
  --out-dir "${OUT_DIR}"
  --queue "${QUEUE}"
  --wall-time "${WALL_TIME}"
  --memory-gb "${MEMORY_GB}"
  --r-module "${R_MODULE}"
  --r-libs-user "${R_LIBS_USER_PATH}"
)

if command -v Rscript >/dev/null 2>&1; then
  exec Rscript gmm_hpc_driver.R "${DRIVER_ARGS[@]}" "$@"
fi

if command -v R >/dev/null 2>&1; then
  exec R --vanilla --slave --file=gmm_hpc_driver.R --args "${DRIVER_ARGS[@]}" "$@"
fi

echo "Neither Rscript nor R is available. Load an R module first." >&2
exit 127

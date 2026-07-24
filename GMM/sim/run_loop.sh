#!/usr/bin/env bash
set -euo pipefail

SEED="${SEED:-1}"
J="${J:-200}"
LOOPS="${LOOPS:-20}"
START_LOOP="${START_LOOP:-1}"
STARTS="${STARTS:-10}"
TARGET="${TARGET:-1}"
TARGET_TYPE="${TARGET_TYPE:-converged}"
TOL="${TOL:-1e-7}"
OUT_DIR="${OUT_DIR:-hpc_runs/seed${SEED}}"
QUEUE="${QUEUE:-${QUEUE_NAME:-serial}}"
WALL_TIME="${WALL_TIME:-14:00}"
MEMORY_GB="${MEMORY_GB:-8}"
R_MODULE="${R_MODULE:-R/4.4.0}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
REQUIRE_CONVERGED_BEST="${REQUIRE_CONVERGED_BEST:-1}"
RUN_INFERENCE="${RUN_INFERENCE:-1}"
INFERENCE_PROB="${INFERENCE_PROB:-0.95}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

DRIVER_ARGS=(
  --seed "${SEED}"
  --J "${J}"
  --loops "${LOOPS}"
  --start-loop "${START_LOOP}"
  --starts "${STARTS}"
  --target "${TARGET}"
  --target-type "${TARGET_TYPE}"
  --tol "${TOL}"
  --out-dir "${OUT_DIR}"
  --queue "${QUEUE}"
  --wall-time "${WALL_TIME}"
  --memory-gb "${MEMORY_GB}"
  --r-module "${R_MODULE}"
  --r-libs-user "${R_LIBS_USER_PATH}"
  --require-converged-best "${REQUIRE_CONVERGED_BEST}"
  --run-inference "${RUN_INFERENCE}"
  --inference-prob "${INFERENCE_PROB}"
)

if command -v Rscript >/dev/null 2>&1; then
  exec Rscript gmm_hpc_driver.R "${DRIVER_ARGS[@]}" "$@"
fi

if command -v R >/dev/null 2>&1; then
  exec R --vanilla --slave --file=gmm_hpc_driver.R --args "${DRIVER_ARGS[@]}" "$@"
fi

echo "Neither Rscript nor R is available. Load an R module first." >&2
exit 127

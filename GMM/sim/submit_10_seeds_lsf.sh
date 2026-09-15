#!/usr/bin/env bash
set -euo pipefail

SEED_FIRST="${SEED_FIRST:-1}"
SEED_LAST="${SEED_LAST:-10}"
J="${J:-200}"
LOOPS="${LOOPS:-6}"
STARTS="${STARTS:-20}"
START_SD="${START_SD:-0.1}"
MAXIT="${MAXIT:-8000}"
OUTER_TOL="${OUTER_TOL:-1e-5}"
RUN_INFERENCE="${RUN_INFERENCE:-1}"
INFERENCE_PROB="${INFERENCE_PROB:-0.95}"
RUN_ROOT="${RUN_ROOT:-hpc_runs/J${J}_10seeds}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

for seed in $(seq "${SEED_FIRST}" "${SEED_LAST}"); do
  echo "Submitting seed ${seed}, J=${J}, starts=${STARTS}, loops=${LOOPS}"
  SEED="${seed}" \
  J="${J}" \
  LOOPS="${LOOPS}" \
  STARTS="${STARTS}" \
  START_SD="${START_SD}" \
  MAXIT="${MAXIT}" \
  OUTER_TOL="${OUTER_TOL}" \
  RUN_INFERENCE="${RUN_INFERENCE}" \
  INFERENCE_PROB="${INFERENCE_PROB}" \
  OUT_DIR="${RUN_ROOT}/seed${seed}" \
    "${SCRIPT_DIR}/submit_controller_lsf.sh"
done

echo "Submitted seeds ${SEED_FIRST}-${SEED_LAST}."
echo "Output root: ${RUN_ROOT}"

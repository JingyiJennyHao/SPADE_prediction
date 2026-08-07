#!/usr/bin/env bash
set -euo pipefail

SEED_FIRST="${SEED_FIRST:-1}"
SEED_LAST="${SEED_LAST:-10}"
J="${J:-200}"
LOOPS="${LOOPS:-20}"
STARTS="${STARTS:-50}"
TARGET="${TARGET:-20}"
TARGET_TYPE="${TARGET_TYPE:-converged}"
START_SD="${START_SD:-0.1}"
MAXIT="${MAXIT:-8000}"
RUN_ROOT="${RUN_ROOT:-hpc_runs/J${J}_10seeds}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

for seed in $(seq "${SEED_FIRST}" "${SEED_LAST}"); do
  echo "Submitting seed ${seed}, J=${J}, starts=${STARTS}, target=${TARGET} ${TARGET_TYPE}"
  SEED="${seed}" \
  J="${J}" \
  LOOPS="${LOOPS}" \
  STARTS="${STARTS}" \
  TARGET="${TARGET}" \
  TARGET_TYPE="${TARGET_TYPE}" \
  START_SD="${START_SD}" \
  MAXIT="${MAXIT}" \
  REQUIRE_CONVERGED_BEST=1 \
  RUN_INFERENCE=1 \
  OUT_DIR="${RUN_ROOT}/seed${seed}" \
    "${SCRIPT_DIR}/submit_controller_lsf.sh"
done

echo "Submitted seeds ${SEED_FIRST}-${SEED_LAST}."
echo "Output root: ${RUN_ROOT}"

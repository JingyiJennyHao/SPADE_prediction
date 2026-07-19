#!/usr/bin/env bash
set -euo pipefail

SEED="${SEED:-1}"
LOOPS="${LOOPS:-20}"
START_LOOP="${START_LOOP:-1}"
STARTS="${STARTS:-120}"
TARGET="${TARGET:-70}"
TARGET_TYPE="${TARGET_TYPE:-completed}"
TOL="${TOL:-1e-7}"
OUT_DIR="${OUT_DIR:-hpc_runs/seed${SEED}}"
QUEUE="${QUEUE:-${QUEUE_NAME:-serial}}"
WORKER_QUEUE="${WORKER_QUEUE:-serial}"
WALL_TIME="${WALL_TIME:-14:00}"
MEMORY_GB="${MEMORY_GB:-8}"
CONTROLLER_WALL_TIME="${CONTROLLER_WALL_TIME:-72:00}"
CONTROLLER_MEMORY_GB="${CONTROLLER_MEMORY_GB:-2}"
R_MODULE="${R_MODULE:-R/4.4.0}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
POLL_SECONDS="${POLL_SECONDS:-300}"
MAX_WAIT_HOURS="${MAX_WAIT_HOURS:-72}"
REQUIRE_CONVERGED_BEST="${REQUIRE_CONVERGED_BEST:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p "${OUT_DIR}/logs"

bsub <<EOF
#!/bin/bash
#BSUB -J gmm_controller_s${SEED}
#BSUB -q ${QUEUE}
#BSUB -W ${CONTROLLER_WALL_TIME}
#BSUB -n 1
#BSUB -R "rusage[mem=${CONTROLLER_MEMORY_GB}GB]"
#BSUB -o ${OUT_DIR}/logs/controller_%J.out
#BSUB -e ${OUT_DIR}/logs/controller_%J.err
set -euo pipefail
cd "${SCRIPT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"

if command -v Rscript >/dev/null 2>&1; then
  Rscript gmm_hpc_driver.R \\
    --seed "${SEED}" \\
    --loops "${LOOPS}" \\
    --start-loop "${START_LOOP}" \\
    --starts "${STARTS}" \\
    --target "${TARGET}" \\
    --target-type "${TARGET_TYPE}" \\
    --tol "${TOL}" \\
    --out-dir "${OUT_DIR}" \\
    --queue "${WORKER_QUEUE}" \\
    --wall-time "${WALL_TIME}" \\
    --memory-gb "${MEMORY_GB}" \\
    --r-module "${R_MODULE}" \\
    --r-libs-user "${R_LIBS_USER_PATH}" \\
    --poll-seconds "${POLL_SECONDS}" \\
    --max-wait-hours "${MAX_WAIT_HOURS}" \\
    --require-converged-best "${REQUIRE_CONVERGED_BEST}"
else
  R --vanilla --slave --file=gmm_hpc_driver.R --args \\
    --seed "${SEED}" \\
    --loops "${LOOPS}" \\
    --start-loop "${START_LOOP}" \\
    --starts "${STARTS}" \\
    --target "${TARGET}" \\
    --target-type "${TARGET_TYPE}" \\
    --tol "${TOL}" \\
    --out-dir "${OUT_DIR}" \\
    --queue "${WORKER_QUEUE}" \\
    --wall-time "${WALL_TIME}" \\
    --memory-gb "${MEMORY_GB}" \\
    --r-module "${R_MODULE}" \\
    --r-libs-user "${R_LIBS_USER_PATH}" \\
    --poll-seconds "${POLL_SECONDS}" \\
    --max-wait-hours "${MAX_WAIT_HOURS}" \\
    --require-converged-best "${REQUIRE_CONVERGED_BEST}"
fi
EOF

echo "Submitted GMM controller."
echo "Driver log: ${OUT_DIR}/logs/driver_seed${SEED}.log"

#!/usr/bin/env bash
set -euo pipefail

# Submit one independent start-to-finish task for every start point in a seed.
# The historical filename is kept so existing launch commands continue to
# work; this script no longer submits a controller job or waits between loops.
SEED="${SEED:-1}"
J="${J:-200}"
LOOPS="${LOOPS:-6}"
STARTS="${STARTS:-20}"
START_SD="${START_SD:-0.1}"
MAXIT="${MAXIT:-8000}"
OUTER_TOL="${OUTER_TOL:-1e-5}"
OUT_DIR="${OUT_DIR:-hpc_runs/seed${SEED}}"
QUEUE="${QUEUE:-${QUEUE_NAME:-serial}}"
WALL_TIME="${WALL_TIME:-72:00}"
MEMORY_GB="${MEMORY_GB:-8}"
R_MODULE="${R_MODULE:-R/4.4.0}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
RUN_INFERENCE="${RUN_INFERENCE:-1}"
INFERENCE_PROB="${INFERENCE_PROB:-0.95}"

if [[ ! "${SEED}" =~ ^[0-9]+$ || ! "${LOOPS}" =~ ^[1-9][0-9]*$ ||
      ! "${STARTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "SEED, LOOPS, and STARTS must be positive integers." >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p "${OUT_DIR}/logs"

array_submit_output="$(bsub <<EOF
#!/bin/bash
#BSUB -J gmm_s${SEED}[1-${STARTS}]
#BSUB -q ${QUEUE}
#BSUB -W ${WALL_TIME}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB}GB]"
#BSUB -o ${OUT_DIR}/logs/start_%I.out
#BSUB -e ${OUT_DIR}/logs/start_%I.err
set -euo pipefail
cd "${SCRIPT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"

if command -v Rscript >/dev/null 2>&1; then
  Rscript gmm_start_worker.R \\
    --seed "${SEED}" \\
    --start-id "\${LSB_JOBINDEX}" \\
    --loops "${LOOPS}" \\
    --out-dir "${OUT_DIR}" \\
    --J "${J}" \\
    --Time 3 \\
    --start-sd "${START_SD}" \\
    --maxit "${MAXIT}" \\
    --outer-tol "${OUTER_TOL}"
else
  R --vanilla --slave --file=gmm_start_worker.R --args \\
    --seed "${SEED}" \\
    --start-id "\${LSB_JOBINDEX}" \\
    --loops "${LOOPS}" \\
    --out-dir "${OUT_DIR}" \\
    --J "${J}" \\
    --Time 3 \\
    --start-sd "${START_SD}" \\
    --maxit "${MAXIT}" \\
    --outer-tol "${OUTER_TOL}"
fi
EOF
)"

array_job_id="$(printf '%s\n' "${array_submit_output}" |
  sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n 1)"
if [ -z "${array_job_id}" ]; then
  echo "Could not parse the start-point array job id from bsub output:" >&2
  echo "${array_submit_output}" >&2
  exit 1
fi

finalize_submit_output="$(bsub -w "ended(${array_job_id})" <<EOF
#!/bin/bash
#BSUB -J gmm_finalize_s${SEED}
#BSUB -q ${QUEUE}
#BSUB -W ${WALL_TIME}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB}GB]"
#BSUB -o ${OUT_DIR}/logs/finalize.out
#BSUB -e ${OUT_DIR}/logs/finalize.err
set -euo pipefail
cd "${SCRIPT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"

if command -v Rscript >/dev/null 2>&1; then
  Rscript gmm_finalize_seed.R \\
    --seed "${SEED}" \\
    --loops "${LOOPS}" \\
    --starts "${STARTS}" \\
    --out-dir "${OUT_DIR}" \\
    --J "${J}" \\
    --Time 3 \\
    --run-inference "${RUN_INFERENCE}" \\
    --inference-prob "${INFERENCE_PROB}"
else
  R --vanilla --slave --file=gmm_finalize_seed.R --args \\
    --seed "${SEED}" \\
    --loops "${LOOPS}" \\
    --starts "${STARTS}" \\
    --out-dir "${OUT_DIR}" \\
    --J "${J}" \\
    --Time 3 \\
    --run-inference "${RUN_INFERENCE}" \\
    --inference-prob "${INFERENCE_PROB}"
fi
EOF
 )"

echo "Submitted ${STARTS} start-to-finish jobs for seed ${SEED}."
echo "Submitted finalizer for seed ${SEED} after array job ${array_job_id}."
echo "${array_submit_output}"
echo "${finalize_submit_output}"
echo "Results: ${OUT_DIR}"

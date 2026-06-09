#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$HOME/SPADE/hpc_long_loop}"
DATA_PATH="${DATA_PATH:-$HOME/SPADE/dat_spade_yr_0212.RDS}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_DIR/results/manual_stage1}"
GROUP_COUNT="${GROUP_COUNT:-210}"
TOTAL_TASKS="${TOTAL_TASKS:-130}"
RANDOM_TASKS="${RANDOM_TASKS:-30}"
JOB_NAME="${JOB_NAME:-manual_s1}"
WALL_TIME="${WALL_TIME:-14:00}"
MEMORY_GB="${MEMORY_GB:-8}"
QUEUE_NAME="${QUEUE_NAME:-}"
R_MODULE="${R_MODULE:-R/4.4}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
WEIGHTS_FILE="${WEIGHTS_FILE:-}"
BETA_REF_FILE="${BETA_REF_FILE:-}"

if (( RANDOM_TASKS < 0 || RANDOM_TASKS > TOTAL_TASKS )); then
  echo "RANDOM_TASKS must be between 0 and TOTAL_TASKS." >&2
  exit 1
fi

mkdir -p "${RESULTS_DIR}/runs" "${RESULTS_DIR}/summary" "${RESULTS_DIR}/logs"

QUEUE_DIRECTIVE=""
if [[ -n "${QUEUE_NAME}" ]]; then
  QUEUE_DIRECTIVE="#BSUB -q ${QUEUE_NAME}"
fi

bsub <<EOF
#!/bin/bash
#BSUB -J ${JOB_NAME}[1-${TOTAL_TASKS}]
#BSUB -W ${WALL_TIME}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB}GB]"
#BSUB -o ${RESULTS_DIR}/logs/${JOB_NAME}.out
#BSUB -e ${RESULTS_DIR}/logs/${JOB_NAME}.err
${QUEUE_DIRECTIVE}
set -euo pipefail
cd "${PROJECT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"
TASK_ID=\${LSB_JOBINDEX}
if (( TASK_ID <= ${RANDOM_TASKS} )); then
  START_FAMILY="random"
elif [[ -n "${BETA_REF_FILE}" ]]; then
  START_FAMILY="beta_ref"
else
  START_FAMILY="fixed"
fi
SEED=\$((500000 + TASK_ID))
CMD=(Rscript stage1_beta_runner.R
  --groups "${GROUP_COUNT}"
  --task-id "\${TASK_ID}"
  --start-family "\${START_FAMILY}"
  --seed "\${SEED}"
  --data "${DATA_PATH}"
  --out "${RESULTS_DIR}/runs/task_\$(printf "%03d" "\${TASK_ID}").rds")
if [[ -n "${WEIGHTS_FILE}" ]]; then
  CMD+=(--weights-file "${WEIGHTS_FILE}")
fi
if [[ -n "${BETA_REF_FILE}" ]]; then
  CMD+=(--beta-ref-file "${BETA_REF_FILE}")
fi
"\${CMD[@]}"
EOF

echo "Submitted stage 1 batch."
echo "Results directory: ${RESULTS_DIR}"
echo "Shared log: ${RESULTS_DIR}/logs/${JOB_NAME}.out"
echo "Shared err: ${RESULTS_DIR}/logs/${JOB_NAME}.err"
echo "After all tasks finish, run:"
echo "cd ${PROJECT_DIR}"
echo "module load ${R_MODULE}"
echo "export R_LIBS_USER=\"${R_LIBS_USER_PATH}\""
echo "Rscript collect_stage1_results.R --input-dir ${RESULTS_DIR}/runs --out-csv ${RESULTS_DIR}/summary/${JOB_NAME}.csv"

#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$HOME/SPADE_prediction/hpc_long_loop}"
DATA_PATH="${DATA_PATH:-$HOME/SPADE/dat_spade_yr_0212.RDS}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_DIR/results/manual_stage2}"
GROUP_COUNT="${GROUP_COUNT:-210}"
TOTAL_TASKS="${TOTAL_TASKS:-130}"
JOB_NAME="${JOB_NAME:-manual_s2}"
WALL_TIME="${WALL_TIME:-14:00}"
MEMORY_GB="${MEMORY_GB:-8}"
QUEUE_NAME="${QUEUE_NAME:-}"
R_MODULE="${R_MODULE:-R/4.4}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
BETA_FILE="${BETA_FILE:-}"

if [[ -z "${BETA_FILE}" ]]; then
  echo "Set BETA_FILE to the beta_hat.rds you want to use." >&2
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
SEED=\$((600000 + TASK_ID))
Rscript stage2_weight_runner.R \
  --groups "${GROUP_COUNT}" \
  --task-id "\${TASK_ID}" \
  --seed "\${SEED}" \
  --beta-file "${BETA_FILE}" \
  --data "${DATA_PATH}" \
  --out "${RESULTS_DIR}/runs/task_\$(printf "%03d" "\${TASK_ID}").rds"
EOF

echo "Submitted stage 2 batch."
echo "Results directory: ${RESULTS_DIR}"
echo "Shared log: ${RESULTS_DIR}/logs/${JOB_NAME}.out"
echo "Shared err: ${RESULTS_DIR}/logs/${JOB_NAME}.err"
echo "After all tasks finish, run:"
echo "cd ${PROJECT_DIR}"
echo "module load ${R_MODULE}"
echo "export R_LIBS_USER=\"${R_LIBS_USER_PATH}\""
echo "Rscript collect_stage2_results.R --input-dir ${RESULTS_DIR}/runs --out-csv ${RESULTS_DIR}/summary/${JOB_NAME}.csv"

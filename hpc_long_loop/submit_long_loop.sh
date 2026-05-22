#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$HOME/SPADE/hpc_long_loop}"
RESULTS_ROOT="${RESULTS_ROOT:-$PROJECT_DIR/results}"
DATA_PATH="${DATA_PATH:-$HOME/SPADE/dat_spade_yr_0212.RDS}"
GROUP_COUNT="${GROUP_COUNT:-210}"
ITERATIONS="${ITERATIONS:-3}"
INITIAL_TOTAL_TASKS="${INITIAL_TOTAL_TASKS:-200}"
INITIAL_RANDOM_TASKS="${INITIAL_RANDOM_TASKS:-100}"
REFINED_TOTAL_TASKS="${REFINED_TOTAL_TASKS:-200}"
REFINED_BETA_REF_TASKS="${REFINED_BETA_REF_TASKS:-100}"
JOB_NAME_PREFIX="${JOB_NAME_PREFIX:-clloop}"
QUEUE_NAME="${QUEUE_NAME:-}"
WALL_TIME_STAGE1="${WALL_TIME_STAGE1:-10:00}"
WALL_TIME_STAGE2="${WALL_TIME_STAGE2:-10:00}"
MEMORY_GB_STAGE1="${MEMORY_GB_STAGE1:-8}"
MEMORY_GB_STAGE2="${MEMORY_GB_STAGE2:-8}"
R_MODULE="${R_MODULE:-R/4.4}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

mkdir -p "${RESULTS_ROOT}"

QUEUE_ARGS=()
if [[ -n "${QUEUE_NAME}" ]]; then
  QUEUE_ARGS=(-q "${QUEUE_NAME}")
fi

submit_stage1() {
  local iteration="$1"
  local total_tasks="$2"
  local leading_family_tasks="$3"
  local leading_family="$4"
  local trailing_family="$5"
  local weights_file="$6"
  local beta_ref_file="$7"
  local depends_on="$8"
  local iter_dir="${RESULTS_ROOT}/iter_$(printf "%02d" "${iteration}")"
  local log_dir="${iter_dir}/logs"
  mkdir -p "${iter_dir}/stage1" "${iter_dir}/stage1_summary" "${log_dir}"

  local dependency_args=()
  if [[ -n "${depends_on}" ]]; then
    dependency_args=(-w "done(${depends_on})")
  fi

  local output
  output=$(bsub "${dependency_args[@]}" "${QUEUE_ARGS[@]}" <<EOF
#!/bin/bash
#BSUB -J ${JOB_NAME_PREFIX}_i$(printf "%02d" "${iteration}")_s1[1-${total_tasks}]
#BSUB -W ${WALL_TIME_STAGE1}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB_STAGE1}GB]"
#BSUB -o ${log_dir}/stage1_%J_%I.out
#BSUB -e ${log_dir}/stage1_%J_%I.err
set -euo pipefail
cd "${PROJECT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"
TASK_ID=\${LSB_JOBINDEX}
if (( TASK_ID <= ${leading_family_tasks} )); then
  START_FAMILY="${leading_family}"
else
  START_FAMILY="${trailing_family}"
fi
SEED=\$((100000 * ${iteration} + TASK_ID))
OUT_FILE="${iter_dir}/stage1/task_\$(printf "%03d" "\${TASK_ID}").rds"
CMD=(Rscript stage1_beta_runner.R
  --groups "${GROUP_COUNT}"
  --task-id "\${TASK_ID}"
  --start-family "\${START_FAMILY}"
  --seed "\${SEED}"
  --data "${DATA_PATH}"
  --out "\${OUT_FILE}")
if [[ -n "${weights_file}" ]]; then
  CMD+=(--weights-file "${weights_file}")
fi
if [[ -n "${beta_ref_file}" ]]; then
  CMD+=(--beta-ref-file "${beta_ref_file}")
fi
"\${CMD[@]}"
EOF
)
  echo "${output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p'
}

submit_summary() {
  local iteration="$1"
  local stage_name="$2"
  local depends_on="$3"
  local iter_dir="${RESULTS_ROOT}/iter_$(printf "%02d" "${iteration}")"
  local log_dir="${iter_dir}/logs"
  local input_dir="${iter_dir}/${stage_name}"
  local output_path="${iter_dir}/${stage_name}_summary/$(if [[ "${stage_name}" == "stage1" ]]; then echo beta_hat; else echo weights_hat; fi).rds"
  mkdir -p "${iter_dir}/${stage_name}_summary" "${log_dir}"

  local script_name
  if [[ "${stage_name}" == "stage1" ]]; then
    script_name="stage1_summarize.R"
  else
    script_name="stage2_summarize.R"
  fi

  local output
  output=$(bsub -w "done(${depends_on})" "${QUEUE_ARGS[@]}" <<EOF
#!/bin/bash
#BSUB -J ${JOB_NAME_PREFIX}_i$(printf "%02d" "${iteration}")_${stage_name}_summary
#BSUB -W 01:00
#BSUB -n 1
#BSUB -R "rusage[mem=2GB]"
#BSUB -o ${log_dir}/${stage_name}_summary_%J.out
#BSUB -e ${log_dir}/${stage_name}_summary_%J.err
set -euo pipefail
cd "${PROJECT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"
Rscript "${script_name}" --input-dir "${input_dir}" --out "${output_path}"
EOF
)
  echo "${output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p'
}

submit_stage2() {
  local iteration="$1"
  local beta_summary_job="$2"
  local iter_dir="${RESULTS_ROOT}/iter_$(printf "%02d" "${iteration}")"
  local log_dir="${iter_dir}/logs"
  local beta_file="${iter_dir}/stage1_summary/beta_hat.rds"
  mkdir -p "${iter_dir}/stage2" "${iter_dir}/stage2_summary" "${log_dir}"

  local output
  output=$(bsub -w "done(${beta_summary_job})" "${QUEUE_ARGS[@]}" <<EOF
#!/bin/bash
#BSUB -J ${JOB_NAME_PREFIX}_i$(printf "%02d" "${iteration}")_s2
#BSUB -W ${WALL_TIME_STAGE2}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB_STAGE2}GB]"
#BSUB -o ${log_dir}/stage2_%J.out
#BSUB -e ${log_dir}/stage2_%J.err
set -euo pipefail
cd "${PROJECT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"
Rscript stage2_weight_runner.R \
  --groups "${GROUP_COUNT}" \
  --task-id 1 \
  --seed $((200000 + ${iteration})) \
  --beta-file "${beta_file}" \
  --data "${DATA_PATH}" \
  --out "${iter_dir}/stage2/task_001.rds"
EOF
)
  echo "${output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p'
}

previous_weights_file=""
previous_beta_file=""
last_job=""

for ((iter = 1; iter <= ITERATIONS; iter++)); do
  if (( iter == 1 )); then
    stage1_job=$(submit_stage1 "${iter}" "${INITIAL_TOTAL_TASKS}" "${INITIAL_RANDOM_TASKS}" "random" "fixed" "" "" "")
  else
    stage1_job=$(submit_stage1 "${iter}" "${REFINED_TOTAL_TASKS}" "${REFINED_BETA_REF_TASKS}" "beta_ref" "fixed" "${previous_weights_file}" "${previous_beta_file}" "${last_job}")
  fi

  stage1_summary_job=$(submit_summary "${iter}" "stage1" "${stage1_job}")
  stage2_job=$(submit_stage2 "${iter}" "${stage1_summary_job}")
  stage2_summary_job=$(submit_summary "${iter}" "stage2" "${stage2_job}")

  previous_beta_file="${RESULTS_ROOT}/iter_$(printf "%02d" "${iter}")/stage1_summary/beta_hat.rds"
  previous_weights_file="${RESULTS_ROOT}/iter_$(printf "%02d" "${iter}")/stage2_summary/weights_hat.rds"
  last_job="${stage2_summary_job}"

  echo "Iteration ${iter}: stage1 job ${stage1_job}, stage1 summary ${stage1_summary_job}, stage2 job ${stage2_job}, stage2 summary ${stage2_summary_job}"
done

echo "All iterations submitted. Final dependency chain ends at job ${last_job}."
echo "Results root: ${RESULTS_ROOT}"

#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$HOME/SPADE/hpc_long_loop}"
RESULTS_ROOT="${RESULTS_ROOT:-$PROJECT_DIR/results/loop}"
GROUP_COUNT="${GROUP_COUNT:-210}"

STAGE1_TOTAL_TASKS="${STAGE1_TOTAL_TASKS:-130}"
STAGE1_RANDOM_TASKS="${STAGE1_RANDOM_TASKS:-30}"
STAGE2_TOTAL_TASKS="${STAGE2_TOTAL_TASKS:-15}"

MIN_STAGE1_RESULTS="${MIN_STAGE1_RESULTS:-100}"
MIN_STAGE2_RESULTS="${MIN_STAGE2_RESULTS:-10}"

MAX_ROUNDS="${MAX_ROUNDS:-50}"
START_ROUND="${START_ROUND:-1}"
WEIGHT_TOL="${WEIGHT_TOL:-1e-7}"
POLL_SECONDS="${POLL_SECONDS:-300}"
POST_KILL_WAIT_SECONDS="${POST_KILL_WAIT_SECONDS:-20}"
KEEP_RUN_RDS="${KEEP_RUN_RDS:-0}"

R_MODULE="${R_MODULE:-R/4.4}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
DATA_PATH="${DATA_PATH:-$HOME/SPADE/dat_spade_yr_0212.RDS}"

CONTROLLER_LOG="${CONTROLLER_LOG:-$RESULTS_ROOT/run_loop_controller.log}"

cd "${PROJECT_DIR}"
mkdir -p "${RESULTS_ROOT}"
touch "${CONTROLLER_LOG}"
exec >>"${CONTROLLER_LOG}" 2>&1

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log() {
  echo "[$(timestamp)] $*"
}

die() {
  log "ERROR: $*"
  exit 1
}

on_error() {
  local exit_code="$1"
  local line_no="$2"
  log "Loop stopped at line ${line_no} with exit code ${exit_code}."
  exit "${exit_code}"
}

trap 'on_error $? $LINENO' ERR

ensure_module_cmd() {
  if command -v module >/dev/null 2>&1; then
    return
  fi

  if [[ -f /etc/profile.d/modules.sh ]]; then
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    source /usr/share/Modules/init/bash
  elif [[ -f /etc/profile.d/z00_lmod.sh ]]; then
    source /etc/profile.d/z00_lmod.sh
  elif [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    source /usr/share/lmod/lmod/init/bash
  fi

  command -v module >/dev/null 2>&1 || die "Could not initialize the 'module' command in this shell."
}

load_r_env() {
  ensure_module_cmd
  module load "${R_MODULE}"
  export R_LIBS_USER="${R_LIBS_USER_PATH}"
}

count_rds() {
  local target_dir="$1"
  if [[ ! -d "${target_dir}" ]]; then
    echo 0
    return
  fi
  find "${target_dir}" -maxdepth 1 -type f -name '*.rds' | wc -l | awk '{print $1}'
}

require_min_results() {
  local target_dir="$1"
  local min_count="$2"
  local label="$3"

  local actual_count
  actual_count=$(count_rds "${target_dir}")

  if (( actual_count < min_count )); then
    die "${label}: only ${actual_count} result files found, need at least ${min_count}."
  fi

  log "${label}: found ${actual_count} result files, meeting minimum ${min_count}."
}

cleanup_run_rds() {
  local target_dir="$1"
  local label="$2"

  if [[ "${KEEP_RUN_RDS}" == "1" ]]; then
    log "${label}: KEEP_RUN_RDS=1, keeping task-level RDS files in ${target_dir}."
    return
  fi

  if [[ ! -d "${target_dir}" ]]; then
    log "${label}: run directory ${target_dir} does not exist, nothing to clean."
    return
  fi

  local rds_count
  rds_count=$(count_rds "${target_dir}")
  if (( rds_count == 0 )); then
    log "${label}: no task-level RDS files to clean in ${target_dir}."
    return
  fi

  find "${target_dir}" -maxdepth 1 -type f -name '*.rds' -delete
  log "${label}: deleted ${rds_count} task-level RDS files from ${target_dir}."
}

job_active() {
  local job_id="$1"
  if bjobs -noheader "${job_id}" >/dev/null 2>&1; then
    return 0
  fi
  return 1
}

cancel_remaining_job_elements() {
  local job_id="$1"
  log "Cancelling remaining elements of job ${job_id}."
  bkill "${job_id}" >/dev/null 2>&1 || true
}

wait_until_ready_or_finished() {
  local job_id="$1"
  local target_dir="$2"
  local expected_count="$3"
  local min_count="$4"
  local label="$5"

  log "${label}: waiting for at least ${min_count}/${expected_count} result files from job ${job_id}."

  while true; do
    local current_count
    current_count=$(count_rds "${target_dir}")

    if (( current_count >= min_count )); then
      log "${label}: reached minimum result threshold ${current_count}/${expected_count} (minimum ${min_count})."
      if job_active "${job_id}"; then
        cancel_remaining_job_elements "${job_id}"
        sleep "${POST_KILL_WAIT_SECONDS}"
      fi
      break
    fi

    if job_active "${job_id}"; then
      log "${label}: job ${job_id} still active, ${current_count}/${expected_count} result files present."
      sleep "${POLL_SECONDS}"
    else
      log "${label}: job ${job_id} is no longer active with ${current_count}/${expected_count} result files present."
      break
    fi
  done

  local final_count
  final_count=$(count_rds "${target_dir}")
  log "${label}: proceeding with ${final_count}/${expected_count} result files present."
}

ensure_round_dirs() {
  local round="$1"
  mkdir -p \
    "${RESULTS_ROOT}/round$(printf "%02d" "${round}")/stage1/runs" \
    "${RESULTS_ROOT}/round$(printf "%02d" "${round}")/stage1/summary" \
    "${RESULTS_ROOT}/round$(printf "%02d" "${round}")/stage2/runs" \
    "${RESULTS_ROOT}/round$(printf "%02d" "${round}")/stage2/summary"
}

stage1_round_dir() {
  local round="$1"
  echo "${RESULTS_ROOT}/round$(printf "%02d" "${round}")/stage1"
}

stage2_round_dir() {
  local round="$1"
  echo "${RESULTS_ROOT}/round$(printf "%02d" "${round}")/stage2"
}

stage1_summary_file() {
  local round="$1"
  echo "$(stage1_round_dir "${round}")/summary/beta_hat$(printf "%02d" "${round}").rds"
}

stage2_summary_file() {
  local round="$1"
  echo "$(stage2_round_dir "${round}")/summary/weights_hat$(printf "%02d" "${round}").rds"
}

submit_stage1_round() {
  local round="$1"
  local weights_file="$2"
  local beta_ref_file="$3"
  local stage1_dir
  stage1_dir="$(stage1_round_dir "${round}")"

  RESULTS_DIR="${stage1_dir}" \
  GROUP_COUNT="${GROUP_COUNT}" \
  TOTAL_TASKS="${STAGE1_TOTAL_TASKS}" \
  RANDOM_TASKS="${STAGE1_RANDOM_TASKS}" \
  JOB_NAME="s1_round$(printf "%02d" "${round}")" \
  DATA_PATH="${DATA_PATH}" \
  WEIGHTS_FILE="${weights_file}" \
  BETA_REF_FILE="${beta_ref_file}" \
  ./submit_stage1_manual.sh
}

submit_stage2_round() {
  local round="$1"
  local beta_file="$2"
  local stage2_dir
  stage2_dir="$(stage2_round_dir "${round}")"

  RESULTS_DIR="${stage2_dir}" \
  GROUP_COUNT="${GROUP_COUNT}" \
  TOTAL_TASKS="${STAGE2_TOTAL_TASKS}" \
  JOB_NAME="s2_round$(printf "%02d" "${round}")" \
  DATA_PATH="${DATA_PATH}" \
  BETA_FILE="${beta_file}" \
  ./submit_stage2_manual.sh
}

summarize_stage1_round() {
  local round="$1"
  local stage1_dir
  stage1_dir="$(stage1_round_dir "${round}")"

  log "Round ${round}: starting stage 1 summary."
  load_r_env

  Rscript collect_stage1_results.R \
    --input-dir "${stage1_dir}/runs" \
    --out-csv "${stage1_dir}/summary/s1_round$(printf "%02d" "${round}").csv"

  Rscript stage1_summarize.R \
    --input-dir "${stage1_dir}/runs" \
    --out "${stage1_dir}/summary/beta_hat$(printf "%02d" "${round}").rds"

  log "Round ${round}: completed stage 1 summary."
}

summarize_stage2_round() {
  local round="$1"
  local stage2_dir
  stage2_dir="$(stage2_round_dir "${round}")"

  log "Round ${round}: starting stage 2 summary."
  load_r_env

  Rscript collect_stage2_results.R \
    --input-dir "${stage2_dir}/runs" \
    --out-csv "${stage2_dir}/summary/s2_round$(printf "%02d" "${round}").csv"

  Rscript stage2_summarize.R \
    --input-dir "${stage2_dir}/runs" \
    --out "${stage2_dir}/summary/weights_hat$(printf "%02d" "${round}").rds"

  log "Round ${round}: completed stage 2 summary."
}

weight_diff() {
  local prev_file="$1"
  local curr_file="$2"
  load_r_env
  Rscript -e '
    prev <- readRDS(commandArgs(trailingOnly = TRUE)[1])$weights_hat
    curr <- readRDS(commandArgs(trailingOnly = TRUE)[2])$weights_hat
    cat(max(abs(curr - prev)), "\n")
  ' "${prev_file}" "${curr_file}"
}

log "Project directory: ${PROJECT_DIR}"
log "Results root: ${RESULTS_ROOT}"
log "Group count: ${GROUP_COUNT}"
log "Stage 1 tasks per round: ${STAGE1_TOTAL_TASKS}"
log "Stage 1 random tasks per round: ${STAGE1_RANDOM_TASKS}"
log "Stage 1 minimum results: ${MIN_STAGE1_RESULTS}"
log "Stage 2 tasks per round: ${STAGE2_TOTAL_TASKS}"
log "Stage 2 minimum results: ${MIN_STAGE2_RESULTS}"
log "Keep task-level RDS files: ${KEEP_RUN_RDS}"
log "Weight convergence tolerance: ${WEIGHT_TOL}"
log "Poll seconds: ${POLL_SECONDS}"
log "Post-kill wait seconds: ${POST_KILL_WAIT_SECONDS}"
log "Controller log: ${CONTROLLER_LOG}"

prev_weight_file=""
prev_beta_file=""

if (( START_ROUND > 1 )); then
  prev_round=$((START_ROUND - 1))
  prev_weight_file="$(stage2_summary_file "${prev_round}")"
  prev_beta_file="$(stage1_summary_file "${prev_round}")"

  if [[ ! -f "${prev_weight_file}" || ! -f "${prev_beta_file}" ]]; then
    die "START_ROUND=${START_ROUND} requires ${prev_beta_file} and ${prev_weight_file}."
  fi
fi

for ((round = START_ROUND; round <= MAX_ROUNDS; round++)); do
  ensure_round_dirs "${round}"
  curr_beta_file="$(stage1_summary_file "${round}")"
  curr_weight_file="$(stage2_summary_file "${round}")"

  log "========== Round ${round}: Stage 1 =========="

  if [[ -f "${curr_beta_file}" ]]; then
    log "Round ${round}: existing beta summary found at ${curr_beta_file}; skipping stage 1 submission."
  else
    if (( round == 1 )); then
      stage1_output="$(submit_stage1_round "${round}" "" "")"
    else
      stage1_output="$(submit_stage1_round "${round}" "${prev_weight_file}" "${prev_beta_file}")"
    fi

    printf '%s\n' "${stage1_output}"
    stage1_job_id="$(printf '%s\n' "${stage1_output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n 1)"
    [[ -n "${stage1_job_id}" ]] || die "Could not parse stage 1 job ID for round ${round}."

    log "Round ${round}: submitted stage 1 as job ${stage1_job_id}."
    wait_until_ready_or_finished "${stage1_job_id}" "$(stage1_round_dir "${round}")/runs" "${STAGE1_TOTAL_TASKS}" "${MIN_STAGE1_RESULTS}" "Round ${round} stage 1"
    require_min_results "$(stage1_round_dir "${round}")/runs" "${MIN_STAGE1_RESULTS}" "Round ${round} stage 1"
  fi

  if [[ ! -f "${curr_beta_file}" ]]; then
    summarize_stage1_round "${round}"
  fi

  [[ -f "${curr_beta_file}" ]] || die "Round ${round}: stage 1 summary did not create ${curr_beta_file}."
  log "Round ${round}: saved beta summary to ${curr_beta_file}"
  cleanup_run_rds "$(stage1_round_dir "${round}")/runs" "Round ${round} stage 1"

  log "========== Round ${round}: Stage 2 =========="

  if [[ -f "${curr_weight_file}" ]]; then
    log "Round ${round}: existing weight summary found at ${curr_weight_file}; skipping stage 2 submission."
  else
    stage2_output="$(submit_stage2_round "${round}" "${curr_beta_file}")"
    printf '%s\n' "${stage2_output}"
    stage2_job_id="$(printf '%s\n' "${stage2_output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n 1)"
    [[ -n "${stage2_job_id}" ]] || die "Could not parse stage 2 job ID for round ${round}."

    log "Round ${round}: submitted stage 2 as job ${stage2_job_id}."
    wait_until_ready_or_finished "${stage2_job_id}" "$(stage2_round_dir "${round}")/runs" "${STAGE2_TOTAL_TASKS}" "${MIN_STAGE2_RESULTS}" "Round ${round} stage 2"
    require_min_results "$(stage2_round_dir "${round}")/runs" "${MIN_STAGE2_RESULTS}" "Round ${round} stage 2"
  fi

  if [[ ! -f "${curr_weight_file}" ]]; then
    summarize_stage2_round "${round}"
  fi

  [[ -f "${curr_weight_file}" ]] || die "Round ${round}: stage 2 summary did not create ${curr_weight_file}."
  log "Round ${round}: saved weight summary to ${curr_weight_file}"
  cleanup_run_rds "$(stage2_round_dir "${round}")/runs" "Round ${round} stage 2"

  if [[ -n "${prev_weight_file}" ]]; then
    diff_val="$(weight_diff "${prev_weight_file}" "${curr_weight_file}")"
    log "Round ${round}: max abs weight change = ${diff_val}"

    if awk "BEGIN {exit !(${diff_val} < ${WEIGHT_TOL})}"; then
      log "Weight convergence reached at round ${round} with tolerance ${WEIGHT_TOL}."
      exit 0
    fi
  fi

  prev_beta_file="${curr_beta_file}"
  prev_weight_file="${curr_weight_file}"
  log "Round ${round}: moving to next round."
done

log "Reached MAX_ROUNDS=${MAX_ROUNDS} without satisfying the weight tolerance ${WEIGHT_TOL}."

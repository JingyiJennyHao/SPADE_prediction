#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$HOME/SPADE_prediction/GMM}"
RESULTS_ROOT="${RESULTS_ROOT:-$PROJECT_DIR/results/iterations}"
GROUP_COUNT="${GROUP_COUNT:-210}"
TOTAL_TASKS="${TOTAL_TASKS:-120}"
MIN_RESULTS="${MIN_RESULTS:-100}"
MAX_ITERATIONS="${MAX_ITERATIONS:-20}"
START_ITERATION="${START_ITERATION:-1}"
POLL_SECONDS="${POLL_SECONDS:-300}"
POST_KILL_WAIT_SECONDS="${POST_KILL_WAIT_SECONDS:-20}"
KEEP_RUN_RDS="${KEEP_RUN_RDS:-0}"

R_MODULE="${R_MODULE:-R/4.4}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
DATA_PATH="${DATA_PATH:-$HOME/SPADE/dat_spade_yr_0212.RDS}"
INITIAL_BETA_FILE="${INITIAL_BETA_FILE:?INITIAL_BETA_FILE is required for iteration 1.}"
START_SCALE="${START_SCALE:-0.1}"
RIDGE="${RIDGE:-1e-6}"
MAXIT="${MAXIT:-8000}"
UNCENTERED_V="${UNCENTERED_V:-0}"
BETA_TOL="${BETA_TOL:-0}"
CONTROLLER_LOG="${CONTROLLER_LOG:-$RESULTS_ROOT/gmm_iteration_controller.log}"

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
  log "GMM iteration loop stopped at line ${line_no} with exit code ${exit_code}."
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
      log "${label}: reached minimum result threshold ${current_count}/${expected_count}."
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
      log "${label}: job ${job_id} finished with ${current_count}/${expected_count} result files present."
      break
    fi
  done
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

iteration_dir() {
  local iter="$1"
  echo "${RESULTS_ROOT}/iter$(printf "%02d" "${iter}")"
}

summary_file() {
  local iter="$1"
  echo "$(iteration_dir "${iter}")/summary/beta_hat$(printf "%02d" "${iter}").rds"
}

ensure_iteration_dirs() {
  local iter="$1"
  mkdir -p "$(iteration_dir "${iter}")/runs" "$(iteration_dir "${iter}")/summary" "$(iteration_dir "${iter}")/logs"
}

submit_iteration() {
  local iter="$1"
  local beta0_file="$2"
  local iter_dir
  iter_dir="$(iteration_dir "${iter}")"

  RESULTS_DIR="${iter_dir}" \
  PROJECT_DIR="${PROJECT_DIR}" \
  GROUP_COUNT="${GROUP_COUNT}" \
  TOTAL_TASKS="${TOTAL_TASKS}" \
  JOB_NAME="gmm_iter$(printf "%02d" "${iter}")" \
  DATA_PATH="${DATA_PATH}" \
  BETA0_FILE="${beta0_file}" \
  START_SCALE="${START_SCALE}" \
  RIDGE="${RIDGE}" \
  MAXIT="${MAXIT}" \
  UNCENTERED_V="${UNCENTERED_V}" \
  R_MODULE="${R_MODULE}" \
  R_LIBS_USER_PATH="${R_LIBS_USER_PATH}" \
  ./submit_gmm_iteration.sh
}

summarize_iteration() {
  local iter="$1"
  local iter_dir
  iter_dir="$(iteration_dir "${iter}")"

  log "Iteration ${iter}: starting GMM summary."
  load_r_env

  Rscript gmm_summarize.R \
    --input-dir "${iter_dir}/runs" \
    --out "${iter_dir}/summary/beta_hat$(printf "%02d" "${iter}").rds"

  log "Iteration ${iter}: completed GMM summary."
}

beta_diff() {
  local prev_file="$1"
  local curr_file="$2"
  load_r_env
  Rscript -e '
    read_beta <- function(path) {
      item <- readRDS(path)
      if (!is.null(item$beta_hat)) return(as.numeric(item$beta_hat))
      if (is.numeric(item)) return(as.numeric(item))
      stop(sprintf("Could not extract beta from %s", path), call. = FALSE)
    }
    prev <- read_beta(commandArgs(trailingOnly = TRUE)[1])
    curr <- read_beta(commandArgs(trailingOnly = TRUE)[2])
    cat(max(abs(curr - prev)), "\n")
  ' "${prev_file}" "${curr_file}"
}

if (( TOTAL_TASKS < MIN_RESULTS )); then
  die "TOTAL_TASKS=${TOTAL_TASKS} must be at least MIN_RESULTS=${MIN_RESULTS}."
fi

log "Project directory: ${PROJECT_DIR}"
log "Results root: ${RESULTS_ROOT}"
log "Group count: ${GROUP_COUNT}"
log "Tasks per iteration: ${TOTAL_TASKS}"
log "Minimum results per iteration: ${MIN_RESULTS}"
log "Initial beta file: ${INITIAL_BETA_FILE}"
log "Start scale: ${START_SCALE}"
log "Ridge: ${RIDGE}"
log "Beta tolerance: ${BETA_TOL}"
log "Keep task-level RDS files: ${KEEP_RUN_RDS}"
log "Controller log: ${CONTROLLER_LOG}"

prev_beta_file="${INITIAL_BETA_FILE}"
if (( START_ITERATION > 1 )); then
  prev_iter=$((START_ITERATION - 1))
  prev_beta_file="$(summary_file "${prev_iter}")"
  [[ -f "${prev_beta_file}" ]] || die "START_ITERATION=${START_ITERATION} requires ${prev_beta_file}."
fi

for ((iter = START_ITERATION; iter <= MAX_ITERATIONS; iter++)); do
  ensure_iteration_dirs "${iter}"
  curr_beta_file="$(summary_file "${iter}")"

  log "========== GMM Iteration ${iter} =========="

  if [[ -f "${curr_beta_file}" ]]; then
    log "Iteration ${iter}: existing beta summary found at ${curr_beta_file}; skipping submission."
  else
    submit_output="$(submit_iteration "${iter}" "${prev_beta_file}")"
    printf '%s\n' "${submit_output}"
    job_id="$(printf '%s\n' "${submit_output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n 1)"
    [[ -n "${job_id}" ]] || die "Could not parse GMM job ID for iteration ${iter}."

    log "Iteration ${iter}: submitted GMM starts as job ${job_id}."
    wait_until_ready_or_finished "${job_id}" "$(iteration_dir "${iter}")/runs" "${TOTAL_TASKS}" "${MIN_RESULTS}" "Iteration ${iter}"
    require_min_results "$(iteration_dir "${iter}")/runs" "${MIN_RESULTS}" "Iteration ${iter}"
  fi

  if [[ ! -f "${curr_beta_file}" ]]; then
    summarize_iteration "${iter}"
  fi

  [[ -f "${curr_beta_file}" ]] || die "Iteration ${iter}: summary did not create ${curr_beta_file}."
  log "Iteration ${iter}: saved best beta summary to ${curr_beta_file}."

  if [[ "${BETA_TOL}" != "0" ]]; then
    diff_val="$(beta_diff "${prev_beta_file}" "${curr_beta_file}")"
    log "Iteration ${iter}: max abs beta change = ${diff_val}"
    if awk "BEGIN {exit !(${diff_val} < ${BETA_TOL})}"; then
      log "Beta convergence reached at iteration ${iter} with tolerance ${BETA_TOL}."
      cleanup_run_rds "$(iteration_dir "${iter}")/runs" "Iteration ${iter}"
      exit 0
    fi
  fi

  cleanup_run_rds "$(iteration_dir "${iter}")/runs" "Iteration ${iter}"
  prev_beta_file="${curr_beta_file}"
  log "Iteration ${iter}: moving to next iteration with beta0=${prev_beta_file}."
done

log "Reached MAX_ITERATIONS=${MAX_ITERATIONS}."

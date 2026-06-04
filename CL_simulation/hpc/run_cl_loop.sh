#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RESULTS_ROOT="${RESULTS_ROOT:-${PROJECT_DIR}/results}"
ROUNDS="${ROUNDS:-20}"
BETA_TASKS="${BETA_TASKS:-100}"
WEIGHT_TASKS="${WEIGHT_TASKS:-100}"
OPTIM_MAXIT="${OPTIM_MAXIT:-8000}"
GROUP_COUNT="${GROUP_COUNT:-200}"
SIM_SEED="${SIM_SEED:-12345}"
R_MODULE="${R_MODULE:-R/4.4}"
R_LIBS_USER_PATH="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
WALL_TIME_BETA="${WALL_TIME_BETA:-14:00}"
WALL_TIME_WEIGHT="${WALL_TIME_WEIGHT:-14:00}"
MEMORY_GB_BETA="${MEMORY_GB_BETA:-8}"
MEMORY_GB_WEIGHT="${MEMORY_GB_WEIGHT:-8}"
QUEUE_NAME="${QUEUE_NAME:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
JOB_PREFIX="${JOB_PREFIX:-clsim}"

cd "${PROJECT_DIR}"
mkdir -p "${RESULTS_ROOT}"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log() {
  echo "[$(timestamp)] $*"
}

die() {
  echo "[$(timestamp)] ERROR: $*" >&2
  exit 1
}

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
}

load_r_env() {
  ensure_module_cmd
  module load "${R_MODULE}"
  export R_LIBS_USER="${R_LIBS_USER_PATH}"
}

job_active() {
  local job_id="$1"
  bjobs -noheader "${job_id}" >/dev/null 2>&1
}

wait_for_job() {
  local job_id="$1"
  local label="$2"

  while job_active "${job_id}"; do
    log "${label}: job ${job_id} still active."
    sleep "${POLL_SECONDS}"
  done
  log "${label}: job ${job_id} finished or left the queue."
}

round_dir() {
  printf "%s/round%02d" "${RESULTS_ROOT}" "$1"
}

beta_csv() {
  printf "%s/beta_results_round%02d.csv" "$(round_dir "$1")" "$1"
}

weight_csv() {
  printf "%s/weight_results_round%02d.csv" "$(round_dir "$1")" "$1"
}

beta_file() {
  printf "%s/best_beta_round%02d.rds" "$(round_dir "$1")" "$1"
}

weight_file() {
  printf "%s/best_weight_round%02d.rds" "$(round_dir "$1")" "$1"
}

queue_directive() {
  if [[ -n "${QUEUE_NAME}" ]]; then
    printf '#BSUB -q %s\n' "${QUEUE_NAME}"
  fi
}

submit_beta_stage() {
  local round="$1"
  local weights_file="$2"
  local beta_ref_file="$3"
  local dir
  dir="$(round_dir "${round}")"
  mkdir -p "${dir}/logs"

  local extra_weights=""
  local extra_beta=""
  if [[ -n "${weights_file}" ]]; then
    extra_weights="--weights-file \"${weights_file}\""
  fi
  if [[ -n "${beta_ref_file}" ]]; then
    extra_beta="--beta-ref-file \"${beta_ref_file}\""
  fi

  local output
  output=$(bsub <<EOF
#!/bin/bash
#BSUB -J ${JOB_PREFIX}_b_r$(printf "%02d" "${round}")[1-${BETA_TASKS}]
#BSUB -W ${WALL_TIME_BETA}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB_BETA}GB]"
#BSUB -o ${dir}/logs/beta_%J_%I.out
#BSUB -e ${dir}/logs/beta_%J_%I.err
$(queue_directive)
set -euo pipefail
cd "${PROJECT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"
TASK_ID=\${LSB_JOBINDEX}
SEED=\$((1000000 + ${round} * 10000 + TASK_ID))
Rscript beta_job.R \
  --round "${round}" \
  --task-id "\${TASK_ID}" \
  --seed "\${SEED}" \
  --sim-seed "${SIM_SEED}" \
  --groups "${GROUP_COUNT}" \
  --out-csv "$(beta_csv "${round}")" \
  --maxit "${OPTIM_MAXIT}" \
  ${extra_weights} \
  ${extra_beta}
EOF
)
  printf '%s\n' "${output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n 1
}

submit_weight_stage() {
  local round="$1"
  local beta_summary="$2"
  local alpha_ref_file="$3"
  local dir
  dir="$(round_dir "${round}")"
  mkdir -p "${dir}/logs"

  local extra_alpha=""
  if [[ -n "${alpha_ref_file}" ]]; then
    extra_alpha="--alpha-ref-file \"${alpha_ref_file}\""
  fi

  local output
  output=$(bsub <<EOF
#!/bin/bash
#BSUB -J ${JOB_PREFIX}_w_r$(printf "%02d" "${round}")[1-${WEIGHT_TASKS}]
#BSUB -W ${WALL_TIME_WEIGHT}
#BSUB -n 1
#BSUB -R "rusage[mem=${MEMORY_GB_WEIGHT}GB]"
#BSUB -o ${dir}/logs/weight_%J_%I.out
#BSUB -e ${dir}/logs/weight_%J_%I.err
$(queue_directive)
set -euo pipefail
cd "${PROJECT_DIR}"
module load "${R_MODULE}"
export R_LIBS_USER="${R_LIBS_USER_PATH}"
TASK_ID=\${LSB_JOBINDEX}
SEED=\$((2000000 + ${round} * 10000 + TASK_ID))
Rscript weight_job.R \
  --round "${round}" \
  --task-id "\${TASK_ID}" \
  --seed "\${SEED}" \
  --sim-seed "${SIM_SEED}" \
  --groups "${GROUP_COUNT}" \
  --beta-file "${beta_summary}" \
  --out-csv "$(weight_csv "${round}")" \
  --maxit "${OPTIM_MAXIT}" \
  ${extra_alpha}
EOF
)
  printf '%s\n' "${output}" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n 1
}

log "Project directory: ${PROJECT_DIR}"
log "Results root: ${RESULTS_ROOT}"
log "Rounds: ${ROUNDS}"
log "Beta tasks per round: ${BETA_TASKS}"
log "Weight tasks per round: ${WEIGHT_TASKS}"
log "optim maxit: ${OPTIM_MAXIT}"
log "Groups: ${GROUP_COUNT}"
log "Simulation seed: ${SIM_SEED}"

previous_beta_file=""
previous_weight_file=""

for ((round = 1; round <= ROUNDS; round++)); do
  dir="$(round_dir "${round}")"
  mkdir -p "${dir}/logs"

  log "========== Round ${round}: beta-step =========="
  if [[ ! -f "$(beta_file "${round}")" ]]; then
    if [[ -f "$(beta_csv "${round}")" ]]; then
      log "Removing existing beta CSV for this round: $(beta_csv "${round}")"
      rm -f "$(beta_csv "${round}")" "$(beta_csv "${round}").lock"
    fi

    beta_job_id="$(submit_beta_stage "${round}" "${previous_weight_file}" "${previous_beta_file}")"
    [[ -n "${beta_job_id}" ]] || die "Could not parse beta job id for round ${round}."
    log "Round ${round}: submitted beta array job ${beta_job_id}."
    wait_for_job "${beta_job_id}" "Round ${round} beta-step"

    load_r_env
    Rscript summarize_beta.R \
      --input-csv "$(beta_csv "${round}")" \
      --out "$(beta_file "${round}")"
  else
    log "Round ${round}: beta summary already exists, skipping beta-step."
  fi

  log "========== Round ${round}: weight-step =========="
  if [[ ! -f "$(weight_file "${round}")" ]]; then
    if [[ -f "$(weight_csv "${round}")" ]]; then
      log "Removing existing weight CSV for this round: $(weight_csv "${round}")"
      rm -f "$(weight_csv "${round}")" "$(weight_csv "${round}").lock"
    fi

    weight_job_id="$(submit_weight_stage "${round}" "$(beta_file "${round}")" "${previous_weight_file}")"
    [[ -n "${weight_job_id}" ]] || die "Could not parse weight job id for round ${round}."
    log "Round ${round}: submitted weight array job ${weight_job_id}."
    wait_for_job "${weight_job_id}" "Round ${round} weight-step"

    load_r_env
    Rscript summarize_weight.R \
      --input-csv "$(weight_csv "${round}")" \
      --out "$(weight_file "${round}")"
  else
    log "Round ${round}: weight summary already exists, skipping weight-step."
  fi

  previous_beta_file="$(beta_file "${round}")"
  previous_weight_file="$(weight_file "${round}")"
  log "Round ${round}: complete."
done

log "Completed ${ROUNDS} rounds."

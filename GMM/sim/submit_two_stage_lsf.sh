#!/usr/bin/env bash
# One controller per seed: independent Stage 1 paths, then Stage 2 refinement array.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: bash submit_two_stage_lsf.sh config_single_replication.R" >&2
  exit 2
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
export GMM_TWO_STAGE_SCRIPT_DIR="$SCRIPT_DIR"
export GMM_TWO_STAGE_CONFIG="$CONFIG_PATH"
export GMM_TWO_STAGE_R_MODULE="${R_MODULE:-R/4.4.0}"
export GMM_TWO_STAGE_R_LIBS="${R_LIBS_USER_PATH:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
export GMM_STAGE1_BACKEND=lsf
export GMM_TWO_STAGE_QUEUE="${QUEUE:-serial}"
export GMM_TWO_STAGE_WALL="${WALL_TIME:-72:00}"
export GMM_TWO_STAGE_MEMORY="${MEMORY_GB:-8}"
mkdir -p logs
bsub -J gmm_two_stage -q "${QUEUE:-serial}" -W "${WALL_TIME:-72:00}" \
  -n 1 -R "rusage[mem=${MEMORY_GB:-8}GB]" \
  -o 'logs/two_stage_%J.out' -e 'logs/two_stage_%J.err' <<'JOB'
#!/usr/bin/env bash
set -euo pipefail
module load "$GMM_TWO_STAGE_R_MODULE"
export R_LIBS_USER="$GMM_TWO_STAGE_R_LIBS"
Rscript "$GMM_TWO_STAGE_SCRIPT_DIR/run_two_stage_sim.R" "$GMM_TWO_STAGE_CONFIG"
JOB

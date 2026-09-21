#!/usr/bin/env bash
set -euo pipefail

# Submit one LSF array for a seed. Each array task owns one start point and
# runs all requested loops sequentially.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/submit_controller_lsf.sh" "$@"

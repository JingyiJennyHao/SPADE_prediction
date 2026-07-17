#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if command -v "${RSCRIPT:-Rscript}" >/dev/null 2>&1; then
  exec "${RSCRIPT:-Rscript}" "${SCRIPT_DIR}/gmm_hpc_driver.R" "$@"
fi

if command -v "${R_BIN:-R}" >/dev/null 2>&1; then
  exec "${R_BIN:-R}" --vanilla --slave --file="${SCRIPT_DIR}/gmm_hpc_driver.R" --args "$@"
fi

echo "Neither Rscript nor R is available. Load an R module first, for example: module avail R" >&2
exit 127

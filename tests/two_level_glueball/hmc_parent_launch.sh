#!/usr/bin/env bash
set -euo pipefail

# Basic local run of gauge-only HMC
# Usage: ./hmc_parent_launch.sh [hmc_parent_ini]

REPO_DIR="${REPO_DIR:-$PWD/chroma}"
CHROMA_BUILD="${CHROMA_BUILD:-/private/tmp/chroma-build/chroma-localbinarydb-0pppatch}"
HMC_BIN="${HMC_BIN:-$CHROMA_BUILD/mainprogs/main/hmc}"
hmc_parent_ini="${1:-hmc_parent.ini.xml}"
hmc_parent_out="${2:-hmc_parent.out.xml}"
hmc_parent_log="${3:-hmc_parent.log.xml}"

"$HMC_BIN" \
  -i "$hmc_parent_ini" \
  -o "$hmc_parent_out" \
  -l "$hmc_parent_log"

echo "Done."
echo "Input XML:  $hmc_parent_ini"
echo "Output XML: $hmc_parent_out"
echo "Log XML:    $hmc_parent_log"

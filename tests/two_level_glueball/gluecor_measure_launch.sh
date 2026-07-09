#!/usr/bin/env bash
set -euo pipefail

# Simple gluecor CSV launcher for every config in one directory.
# Usage: ./gluecor_measure_launch.sh <config_dir>

CHROMA_BUILD="${CHROMA_BUILD:-/private/tmp/chroma-build/chroma-localbinarydb-0pppatch}"
MEASURE_BIN="${MEASURE_BIN:-$CHROMA_BUILD/mainprogs/main/gluecor_measure}"

NROW="${NROW:-8 8 8 16}"
DECAY_DIR="${DECAY_DIR:-3}"
BL_LEVEL="${BL_LEVEL:-1}"
BLK_ACCU="${BLK_ACCU:-1.0e-5}"
BLK_MAX="${BLK_MAX:-50}"

CFG_TYPE="${CFG_TYPE:-auto}"
CFG_TYPE_LIME="${CFG_TYPE_LIME:-SZINQIO}"
CFG_TYPE_SCIDAC="${CFG_TYPE_SCIDAC:-SCIDAC}"

usage() {
  echo "Usage: $0 <config_dir>" >&2
  echo "Set NROW=\"n0 n1 n2 n3\" to match the configs in that directory." >&2
}

infer_cfg_type() {
  local cfg_file="$1"

  if [[ "$CFG_TYPE" != "auto" ]]; then
    printf '%s\n' "$CFG_TYPE"
    return
  fi

  case "$cfg_file" in
    *.scidac)
      printf '%s\n' "$CFG_TYPE_SCIDAC"
      ;;
    *.lime)
      printf '%s\n' "$CFG_TYPE_LIME"
      ;;
    *)
      echo "ERROR: could not infer cfg type for ${cfg_file}" >&2
      echo "Set CFG_TYPE, CFG_TYPE_LIME, or CFG_TYPE_SCIDAC explicitly." >&2
      exit 1
      ;;
  esac
}

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

if [[ ! -d "$1" ]]; then
  echo "ERROR: config directory not found: $1" >&2
  exit 1
fi

if [[ ! -x "$MEASURE_BIN" ]]; then
  echo "ERROR: measurement binary not found or not executable: $MEASURE_BIN" >&2
  exit 1
fi

read -r -a nrow_args <<< "$NROW"
if [[ "${#nrow_args[@]}" -ne 4 ]]; then
  echo "ERROR: NROW must contain four integers, got: $NROW" >&2
  exit 1
fi

INPUT_DIR="$(cd "$1" && pwd)"
OUTPUT_ROOT="${OUTPUT_ROOT:-$INPUT_DIR/gluecor_output}"

mkdir -p "$OUTPUT_ROOT"

count=0

while IFS= read -r -d '' cfg_file; do
  cfg_name="$(basename "$cfg_file")"
  cfg_stem="${cfg_name%.*}"
  cfg_type="$(infer_cfg_type "$cfg_file")"
  csv_file="$OUTPUT_ROOT/${cfg_stem}.gluecor.csv"

  echo "Measuring ${cfg_name}"
  "$MEASURE_BIN" \
    --nrow "${nrow_args[@]}" \
    --cfg-type "$cfg_type" \
    --decay-dir "$DECAY_DIR" \
    --bl-level "$BL_LEVEL" \
    --blk-accu "$BLK_ACCU" \
    --blk-max "$BLK_MAX" \
    "$cfg_file" \
    "$csv_file"
  echo "  csv: ${csv_file}"

  count=$((count + 1))
done < <(find "$INPUT_DIR" -maxdepth 1 \( -type f -o -type l \) \( -name '*.lime' -o -name '*.scidac' \) -print0)

if [[ "$count" -eq 0 ]]; then
  echo "ERROR: no .lime or .scidac config files found in ${INPUT_DIR}" >&2
  exit 1
fi

echo "Done."
echo "Processed configs: ${count}"
echo "Output directory: ${OUTPUT_ROOT}"

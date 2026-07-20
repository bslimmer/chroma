#!/usr/bin/env bash
set -euo pipefail

# Simple gluecor CSV launcher for every config in one directory.
# Usage: ./gluecor_measure_launch.sh <config_dir>

CHROMA_BUILD="${CHROMA_BUILD:-/private/tmp/chroma-build-codex-gluecor/chroma-localbinarydb}"
MEASURE_BIN="${MEASURE_BIN:-}"

#Below is for child matrices
NROW="${NROW:-8 8 8 9}"
#Below is for parent matrices
#NROW="${NROW:-8 8 8 16}"
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
  echo "Optional overrides: CHROMA_BUILD=/path/to/chroma-build or MEASURE_BIN=/path/to/gluecor_measure." >&2
}

resolve_measure_bin() {
  local candidate

  if [[ -n "$MEASURE_BIN" ]]; then
    printf '%s\n' "$MEASURE_BIN"
    return
  fi

  for candidate in \
    "$CHROMA_BUILD/mainprogs/main/gluecor_measure" \
    "/private/tmp/chroma-build-codex-gluecor/chroma-localbinarydb/mainprogs/main/gluecor_measure" \
    "/private/tmp/chroma-build-codex-gluecor/chroma-localbinarydb-0pppatch/mainprogs/main/gluecor_measure" \
    "/private/tmp/chroma-build/chroma-localbinarydb-0pppatch/mainprogs/main/gluecor_measure"
  do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
  done

  printf '%s\n' "$CHROMA_BUILD/mainprogs/main/gluecor_measure"
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

MEASURE_BIN="$(resolve_measure_bin)"

if [[ ! -x "$MEASURE_BIN" ]]; then
  echo "ERROR: measurement binary not found or not executable: $MEASURE_BIN" >&2
  echo "Set MEASURE_BIN to the executable path or CHROMA_BUILD to the matching build tree." >&2
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

total_count=0
measured_count=0
skipped_count=0

while IFS= read -r -d '' cfg_file; do
  total_count=$((total_count + 1))
  cfg_name="$(basename "$cfg_file")"
  cfg_stem="${cfg_name%.*}"
  cfg_type="$(infer_cfg_type "$cfg_file")"
  csv_file="$OUTPUT_ROOT/${cfg_stem}.gluecor.csv"

  if [[ -e "$csv_file" ]]; then
    echo "Skipping ${cfg_name}: found existing csv ${csv_file}"
    skipped_count=$((skipped_count + 1))
    continue
  fi

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

  measured_count=$((measured_count + 1))
done < <(find "$INPUT_DIR" -maxdepth 1 \( -type f -o -type l \) \( -name '*.lime' -o -name '*.scidac' \) -print0)

if [[ "$total_count" -eq 0 ]]; then
  echo "ERROR: no .lime or .scidac config files found in ${INPUT_DIR}" >&2
  exit 1
fi

echo "Done."
echo "Found configs: ${total_count}"
echo "Measured configs: ${measured_count}"
echo "Skipped existing: ${skipped_count}"
echo "Output directory: ${OUTPUT_ROOT}"

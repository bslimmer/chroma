#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: generate_two_level_0pp_xml_bundle.sh [options] [output_root [first_update [last_update [outer_step [child_updates [discard_updates [child_save_interval [bl_level_selected]]]]]]]]

Defaults match the first concrete 8^4, 40 x 10 target:
  output_root         cfgs/two_level_0pp_40x10
  first_update        100
  last_update         4000
  outer_step          100
  child_updates       10
  discard_updates     0
  child_save_interval 1
  bl_level_selected   1

Options:
  --parent-seed N
      Parent HMC RNG seed. Default: 37.
  --reuse-parent-root PATH
      Reuse the aligned retained parent configs already saved under PATH.
      The generator will point split/measurement XML for the reused prefix at
      PATH/parent_outer_cfg_<update>.lime.
  --reuse-parent-last-update N
      Last retained parent update to reuse from --reuse-parent-root. The
      reused prefix is FIRST_UPDATE..N in steps of OUTER_STEP. If N < LAST_UPDATE,
      the generated parent HMC XML starts a fresh HMC leg from the saved gauge
      config PATH/parent_outer_cfg_N.lime and only produces the missing later
      retained parent samples under OUTPUT_ROOT.

The script writes concrete run XML under:
  <output_root>/xml/
EOF
}

POSITIONAL_ARGS=()
PARENT_SEED=37
REUSE_PARENT_ROOT=""
REUSE_PARENT_LAST_UPDATE=""

while (( $# > 0 )); do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --parent-seed)
      if [[ -z "${2:-}" ]]; then
        echo "--parent-seed requires an integer argument" >&2
        exit 1
      fi
      PARENT_SEED="$2"
      shift 2
      ;;
    --reuse-parent-root)
      if [[ -z "${2:-}" ]]; then
        echo "--reuse-parent-root requires a path argument" >&2
        exit 1
      fi
      REUSE_PARENT_ROOT="$2"
      shift 2
      ;;
    --reuse-parent-last-update)
      if [[ -z "${2:-}" ]]; then
        echo "--reuse-parent-last-update requires an integer argument" >&2
        exit 1
      fi
      REUSE_PARENT_LAST_UPDATE="$2"
      shift 2
      ;;
    --)
      shift
      while (( $# > 0 )); do
        POSITIONAL_ARGS+=("$1")
        shift
      done
      ;;
    -*)
      echo "unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}"

OUTPUT_ROOT="${1:-cfgs/two_level_0pp_40x10}"
FIRST_UPDATE="${2:-100}"
LAST_UPDATE="${3:-4000}"
OUTER_STEP="${4:-100}"
CHILD_UPDATES="${5:-10}"
DISCARD_UPDATES="${6:-0}"
CHILD_SAVE_INTERVAL="${7:-1}"
BL_LEVEL_SELECTED="${8:-1}"

CHILD0_SEED_BASE=1000
CHILD1_SEED_BASE=2000
PARENT_TOTAL_UPDATES="${LAST_UPDATE}"
XML_ROOT="${OUTPUT_ROOT}/xml"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if (( OUTER_STEP <= 0 )); then
  echo "outer_step must be positive" >&2
  exit 1
fi
if (( LAST_UPDATE < FIRST_UPDATE )); then
  echo "last_update must be >= first_update" >&2
  exit 1
fi
if (( (LAST_UPDATE - FIRST_UPDATE) % OUTER_STEP != 0 )); then
  echo "outer updates must form an arithmetic progression with the given step" >&2
  exit 1
fi
if (( CHILD_UPDATES <= DISCARD_UPDATES )); then
  echo "child_updates must exceed discard_updates" >&2
  exit 1
fi
if (( CHILD_SAVE_INTERVAL <= 0 )); then
  echo "child_save_interval must be positive" >&2
  exit 1
fi
if (( BL_LEVEL_SELECTED < 0 )); then
  echo "bl_level_selected must be non-negative" >&2
  exit 1
fi
if [[ -n "${REUSE_PARENT_ROOT}" && -z "${REUSE_PARENT_LAST_UPDATE}" ]]; then
  echo "--reuse-parent-last-update is required when --reuse-parent-root is set" >&2
  exit 1
fi
if [[ -n "${REUSE_PARENT_LAST_UPDATE}" && -z "${REUSE_PARENT_ROOT}" ]]; then
  echo "--reuse-parent-root is required when --reuse-parent-last-update is set" >&2
  exit 1
fi

retained_span=$(( CHILD_UPDATES - DISCARD_UPDATES ))
if (( retained_span % CHILD_SAVE_INTERVAL != 0 )); then
  echo "child_updates - discard_updates must be divisible by child_save_interval" >&2
  exit 1
fi

EXPECTED_MEASUREMENTS=$(( retained_span / CHILD_SAVE_INTERVAL ))
OUTER_COUNT=$(( (LAST_UPDATE - FIRST_UPDATE) / OUTER_STEP + 1 ))

PARENT_START_CFG_TYPE="UNIT"
PARENT_START_CFG_FILE="DUMMY"
PARENT_START_UPDATE=0
PARENT_UPDATES_THIS_RUN="${LAST_UPDATE}"
REUSED_OUTER_COUNT=0

if [[ -n "${REUSE_PARENT_ROOT}" ]]; then
  if [[ ! -d "${REUSE_PARENT_ROOT}" ]]; then
    echo "reuse-parent-root does not exist: ${REUSE_PARENT_ROOT}" >&2
    exit 1
  fi
  if (( REUSE_PARENT_LAST_UPDATE < FIRST_UPDATE )); then
    echo "reuse-parent-last-update must be >= first_update" >&2
    exit 1
  fi
  if (( REUSE_PARENT_LAST_UPDATE > LAST_UPDATE )); then
    echo "reuse-parent-last-update must be <= last_update" >&2
    exit 1
  fi
  if (( (REUSE_PARENT_LAST_UPDATE - FIRST_UPDATE) % OUTER_STEP != 0 )); then
    echo "reuse-parent-last-update must lie on the requested outer-sample progression" >&2
    exit 1
  fi

  for (( update_no = FIRST_UPDATE; update_no <= REUSE_PARENT_LAST_UPDATE; update_no += OUTER_STEP )); do
    reuse_cfg_file="${REUSE_PARENT_ROOT}/parent_outer_cfg_${update_no}.lime"
    if [[ ! -f "${reuse_cfg_file}" ]]; then
      echo "missing reused parent config: ${reuse_cfg_file}" >&2
      exit 1
    fi
  done

  PARENT_START_CFG_TYPE="SCIDAC"
  PARENT_START_CFG_FILE="${REUSE_PARENT_ROOT}/parent_outer_cfg_${REUSE_PARENT_LAST_UPDATE}.lime"
  PARENT_START_UPDATE="${REUSE_PARENT_LAST_UPDATE}"
  PARENT_UPDATES_THIS_RUN=$(( LAST_UPDATE - REUSE_PARENT_LAST_UPDATE ))
  REUSED_OUTER_COUNT=$(( (REUSE_PARENT_LAST_UPDATE - FIRST_UPDATE) / OUTER_STEP + 1 ))
fi

mkdir -p "${XML_ROOT}"

parent_cfg_file_for_update() {
  local update_no="$1"

  if [[ -n "${REUSE_PARENT_ROOT}" ]] && (( update_no <= REUSE_PARENT_LAST_UPDATE )); then
    printf '%s\n' "${REUSE_PARENT_ROOT}/parent_outer_cfg_${update_no}.lime"
  else
    printf '%s\n' "${OUTPUT_ROOT}/parent_outer_cfg_${update_no}.lime"
  fi
}

render_template() {
  local template_path="$1"
  local output_path="$2"
  local outer_update="$3"
  local outer_sample_id="$4"
  local child_id="$5"
  local child_seed="$6"
  local meas_update="$7"
  local stream_id="$8"
  local parent_cfg_file="$9"

  sed \
    -e "s|@OUTPUT_ROOT@|${OUTPUT_ROOT}|g" \
    -e "s|@OUTER_SAMPLE_ID@|${outer_sample_id}|g" \
    -e "s|@OUTER_UPDATE@|${outer_update}|g" \
    -e "s|@CHILD_ID@|${child_id}|g" \
    -e "s|@CHILD_SEED@|${child_seed}|g" \
    -e "s|@PARENT_TOTAL_UPDATES@|${PARENT_TOTAL_UPDATES}|g" \
    -e "s|@OUTER_SAVE_INTERVAL@|${OUTER_STEP}|g" \
    -e "s|@PARENT_SEED@|${PARENT_SEED}|g" \
    -e "s|@PARENT_START_CFG_TYPE@|${PARENT_START_CFG_TYPE}|g" \
    -e "s|@PARENT_START_CFG_FILE@|${PARENT_START_CFG_FILE}|g" \
    -e "s|@PARENT_START_UPDATE@|${PARENT_START_UPDATE}|g" \
    -e "s|@PARENT_UPDATES_THIS_RUN@|${PARENT_UPDATES_THIS_RUN}|g" \
    -e "s|@PARENT_CFG_FILE@|${parent_cfg_file}|g" \
    -e "s|@CHILD_UPDATES@|${CHILD_UPDATES}|g" \
    -e "s|@DISCARD_UPDATES@|${DISCARD_UPDATES}|g" \
    -e "s|@CHILD_SAVE_INTERVAL@|${CHILD_SAVE_INTERVAL}|g" \
    -e "s|@EXPECTED_MEASUREMENTS@|${EXPECTED_MEASUREMENTS}|g" \
    -e "s|@BL_LEVEL_SELECTED@|${BL_LEVEL_SELECTED}|g" \
    -e "s|@MEAS_UPDATE@|${meas_update}|g" \
    -e "s|@STREAM_ID@|${stream_id}|g" \
    -e "s|@UPDATE_NO@|${outer_update}|g" \
    -e "s|@CFG_FILE@|${OUTPUT_ROOT}/parent_outer_cfg_${outer_update}.lime|g" \
    "${template_path}" > "${output_path}"
}

render_block_template() {
  local template_path="$1"
  local output_path="$2"
  local block_placeholder="$3"
  local block_value="$4"
  local outer_update="$5"
  local outer_sample_id="$6"
  local child_id="$7"

  sed \
    -e "s|@OUTPUT_ROOT@|${OUTPUT_ROOT}|g" \
    -e "s|@OUTER_SAMPLE_ID@|${outer_sample_id}|g" \
    -e "s|@OUTER_UPDATE@|${outer_update}|g" \
    -e "s|@CHILD_ID@|${child_id}|g" \
    -e "s|@DISCARD_UPDATES@|${DISCARD_UPDATES}|g" \
    -e "s|@EXPECTED_MEASUREMENTS@|${EXPECTED_MEASUREMENTS}|g" \
    -e "s|@BL_LEVEL_SELECTED@|${BL_LEVEL_SELECTED}|g" \
    "${template_path}" | \
    BLOCK_PLACEHOLDER="${block_placeholder}" \
    BLOCK_VALUE="${block_value}" \
    perl -0pe 's/\Q$ENV{BLOCK_PLACEHOLDER}\E/$ENV{BLOCK_VALUE}/g' > "${output_path}"
}

render_template \
  "${SCRIPT_DIR}/hmc_parent.0pp_two_level_outer.template.ini.xml" \
  "${XML_ROOT}/hmc_parent.0pp_two_level_outer.ini.xml" \
  "${FIRST_UPDATE}" \
  "outer_${FIRST_UPDATE}" \
  "0" \
  "${CHILD0_SEED_BASE}" \
  "0" \
  "0" \
  "${PARENT_START_CFG_FILE}"

outer_sample_summary_lines=""
sample_index=0
for (( update_no = FIRST_UPDATE; update_no <= LAST_UPDATE; update_no += OUTER_STEP )); do
  sample_index=$(( sample_index + 1 ))
  outer_sample_id="outer_${update_no}"
  child0_seed=$(( CHILD0_SEED_BASE + sample_index ))
  child1_seed=$(( CHILD1_SEED_BASE + sample_index ))
  parent_cfg_file="$(parent_cfg_file_for_update "${update_no}")"

  mkdir -p "${OUTPUT_ROOT}/${outer_sample_id}"

  render_template \
    "${SCRIPT_DIR}/gauge_subdomain_split.0pp_two_level.template.ini.xml" \
    "${XML_ROOT}/gauge_subdomain_split.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}" \
    "0" \
    "0" \
    "${parent_cfg_file}"

  render_template \
    "${SCRIPT_DIR}/hmc_child.temporal_zone_glueball_0pp_2lvl.template.ini.xml" \
    "${XML_ROOT}/hmc_child0.temporal_zone_glueball_0pp_2lvl.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}" \
    "0" \
    "0" \
    "${parent_cfg_file}"

  render_template \
    "${SCRIPT_DIR}/hmc_child.temporal_zone_glueball_0pp_2lvl.template.ini.xml" \
    "${XML_ROOT}/hmc_child1.temporal_zone_glueball_0pp_2lvl.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "1" \
    "${child1_seed}" \
    "0" \
    "0" \
    "${parent_cfg_file}"

  render_template \
    "${SCRIPT_DIR}/measure_glueball_0pp_parent.template.ini.xml" \
    "${XML_ROOT}/measure_glueball_0pp_parent.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "-1" \
    "0" \
    "0" \
    "0" \
    "${parent_cfg_file}"

  child0_measurement_lines=""
  child1_measurement_lines=""
  for (( meas_update = CHILD_SAVE_INTERVAL; meas_update <= CHILD_UPDATES; meas_update += CHILD_SAVE_INTERVAL )); do
    if (( meas_update <= DISCARD_UPDATES )); then
      continue
    fi

    render_template \
      "${SCRIPT_DIR}/measure_glueball_0pp_child.template.ini.xml" \
      "${XML_ROOT}/measure_glueball_0pp_child0.meas_${meas_update}.${outer_sample_id}.ini.xml" \
      "${update_no}" \
      "${outer_sample_id}" \
      "0" \
      "${child0_seed}" \
      "${meas_update}" \
      "0" \
      "${parent_cfg_file}"

    render_template \
      "${SCRIPT_DIR}/measure_glueball_0pp_child.template.ini.xml" \
      "${XML_ROOT}/measure_glueball_0pp_child1.meas_${meas_update}.${outer_sample_id}.ini.xml" \
      "${update_no}" \
      "${outer_sample_id}" \
      "1" \
      "${child1_seed}" \
      "${meas_update}" \
      "0" \
      "${parent_cfg_file}"

    child0_measurement_lines="${child0_measurement_lines}      <elem>${OUTPUT_ROOT}/${outer_sample_id}/glueball_0pp.child0.stream0.meas_${meas_update}.summary.xml</elem>
"
    child1_measurement_lines="${child1_measurement_lines}      <elem>${OUTPUT_ROOT}/${outer_sample_id}/glueball_0pp.child1.stream0.meas_${meas_update}.summary.xml</elem>
"
  done

  render_template \
    "${SCRIPT_DIR}/glueball_0pp_parent_window.template.check.ini.xml" \
    "${XML_ROOT}/glueball_0pp_parent_window.${outer_sample_id}.check.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "-1" \
    "0" \
    "0" \
    "0" \
    "${parent_cfg_file}"

  render_block_template \
    "${SCRIPT_DIR}/glueball_0pp_child.template.check.ini.xml" \
    "${XML_ROOT}/glueball_0pp_child0.${outer_sample_id}.check.ini.xml" \
    "@MEASUREMENT_FILES@" \
    "${child0_measurement_lines}" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0"

  render_block_template \
    "${SCRIPT_DIR}/glueball_0pp_child.template.check.ini.xml" \
    "${XML_ROOT}/glueball_0pp_child1.${outer_sample_id}.check.ini.xml" \
    "@MEASUREMENT_FILES@" \
    "${child1_measurement_lines}" \
    "${update_no}" \
    "${outer_sample_id}" \
    "1"

  render_template \
    "${SCRIPT_DIR}/glueball_0pp_outer_sample.template.check.ini.xml" \
    "${XML_ROOT}/glueball_0pp_outer_sample.${outer_sample_id}.check.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}" \
    "0" \
    "0" \
    "${parent_cfg_file}"

  outer_sample_summary_lines="${outer_sample_summary_lines}      <elem>${OUTPUT_ROOT}/${outer_sample_id}/glueball_0pp.two_level.summary.xml</elem>
"
done

cat > "${XML_ROOT}/glueball_0pp_two_level.check.ini.xml" <<EOF
<?xml version="1.0"?>
<glueball_0pp_corr_check>
  <Input>
    <outer_sample_summaries>
${outer_sample_summary_lines}    </outer_sample_summaries>
  </Input>

  <Output>
    <summary_file>${OUTPUT_ROOT}/glueball_0pp.two_level.outer_ensemble.summary.xml</summary_file>
    <csv_file>${OUTPUT_ROOT}/glueball_0pp.two_level.outer_ensemble.csv</csv_file>
  </Output>

  <Reducer>
    <mode>TWO_LEVEL_OUTER_ENSEMBLE</mode>
    <compare_delta_t_parent_set>2 3 4 5 6</compare_delta_t_parent_set>
  </Reducer>

  <Checks>
    <min_outer_samples>${OUTER_COUNT}</min_outer_samples>
  </Checks>
</glueball_0pp_corr_check>
EOF

echo "Generated two-level 0++ XML bundle under ${XML_ROOT}"
echo "  parent_hmc: ${XML_ROOT}/hmc_parent.0pp_two_level_outer.ini.xml"
echo "  outer_samples: ${OUTER_COUNT} (${FIRST_UPDATE}..${LAST_UPDATE} step ${OUTER_STEP})"
if [[ -n "${REUSE_PARENT_ROOT}" ]]; then
  echo "  reused_parent_prefix: ${REUSED_OUTER_COUNT} samples from ${REUSE_PARENT_ROOT} (${FIRST_UPDATE}..${REUSE_PARENT_LAST_UPDATE} step ${OUTER_STEP})"
  if (( PARENT_UPDATES_THIS_RUN > 0 )); then
    echo "  parent_extension_start: ${PARENT_START_CFG_FILE} (fresh HMC leg, StartUpdateNum=${PARENT_START_UPDATE})"
  else
    echo "  parent_extension_start: none (requested outer samples fully covered by reused parent prefix)"
  fi
else
  echo "  reused_parent_prefix: none"
fi
echo "  retained_child_measurements_per_stream: ${EXPECTED_MEASUREMENTS}"
echo "  outer_ensemble_check: ${XML_ROOT}/glueball_0pp_two_level.check.ini.xml"

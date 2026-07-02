#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: generate_two_level_0mp_xml_bundle.sh [output_root [first_update [last_update [outer_step [child_updates [discard_updates [child_qactden_frequency]]]]]]]

Defaults match the current "one-hour" preset:
  output_root               cfgs/two_level_0mp_1h
  first_update              100
  last_update               1200
  outer_step                100
  child_updates             180
  discard_updates           20
  child_qactden_frequency   5

The script writes concrete run XML under:
  <output_root>/xml/
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

OUTPUT_ROOT="${1:-cfgs/two_level_0mp_1h}"
FIRST_UPDATE="${2:-100}"
LAST_UPDATE="${3:-1200}"
OUTER_STEP="${4:-100}"
CHILD_UPDATES="${5:-180}"
DISCARD_UPDATES="${6:-20}"
CHILD_QACTDEN_FREQUENCY="${7:-5}"

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
if (( CHILD_QACTDEN_FREQUENCY <= 0 )); then
  echo "child_qactden_frequency must be positive" >&2
  exit 1
fi

retained_span=$(( CHILD_UPDATES - DISCARD_UPDATES ))
if (( retained_span % CHILD_QACTDEN_FREQUENCY != 0 )); then
  echo "child_updates - discard_updates must be divisible by child_qactden_frequency" >&2
  exit 1
fi

EXPECTED_MEASUREMENTS=$(( retained_span / CHILD_QACTDEN_FREQUENCY ))
OUTER_COUNT=$(( (LAST_UPDATE - FIRST_UPDATE) / OUTER_STEP + 1 ))

mkdir -p "${XML_ROOT}"

render_template() {
  local template_path="$1"
  local output_path="$2"
  local outer_update="$3"
  local outer_sample_id="$4"
  local child_id="$5"
  local child_seed="$6"

  sed \
    -e "s|@OUTPUT_ROOT@|${OUTPUT_ROOT}|g" \
    -e "s|@OUTER_SAMPLE_ID@|${outer_sample_id}|g" \
    -e "s|@OUTER_UPDATE@|${outer_update}|g" \
    -e "s|@CHILD_ID@|${child_id}|g" \
    -e "s|@CHILD_SEED@|${child_seed}|g" \
    -e "s|@PARENT_TOTAL_UPDATES@|${PARENT_TOTAL_UPDATES}|g" \
    -e "s|@OUTER_SAVE_INTERVAL@|${OUTER_STEP}|g" \
    -e "s|@CHILD_UPDATES@|${CHILD_UPDATES}|g" \
    -e "s|@DISCARD_UPDATES@|${DISCARD_UPDATES}|g" \
    -e "s|@CHILD_QACTDEN_FREQUENCY@|${CHILD_QACTDEN_FREQUENCY}|g" \
    -e "s|@EXPECTED_MEASUREMENTS@|${EXPECTED_MEASUREMENTS}|g" \
    "${template_path}" > "${output_path}"
}

render_template \
  "${SCRIPT_DIR}/hmc_parent.0mp_two_level_outer.template.ini.xml" \
  "${XML_ROOT}/hmc_parent.0mp_two_level_outer.ini.xml" \
  "${FIRST_UPDATE}" \
  "outer_${FIRST_UPDATE}" \
  "0" \
  "${CHILD0_SEED_BASE}"

outer_sample_summary_lines=""
sample_index=0
for (( update_no = FIRST_UPDATE; update_no <= LAST_UPDATE; update_no += OUTER_STEP )); do
  sample_index=$(( sample_index + 1 ))
  outer_sample_id="outer_${update_no}"
  child0_seed=$(( CHILD0_SEED_BASE + sample_index ))
  child1_seed=$(( CHILD1_SEED_BASE + sample_index ))

  mkdir -p "${OUTPUT_ROOT}/${outer_sample_id}"

  render_template \
    "${SCRIPT_DIR}/gauge_subdomain_split.0mp_two_level.template.ini.xml" \
    "${XML_ROOT}/gauge_subdomain_split.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}"

  render_template \
    "${SCRIPT_DIR}/hmc_child.temporal_zone_qactden_0mp_2lvl.template.ini.xml" \
    "${XML_ROOT}/hmc_child0.temporal_zone_qactden_0mp_2lvl.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}"

  render_template \
    "${SCRIPT_DIR}/hmc_child.temporal_zone_qactden_0mp_2lvl.template.ini.xml" \
    "${XML_ROOT}/hmc_child1.temporal_zone_qactden_0mp_2lvl.${outer_sample_id}.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "1" \
    "${child1_seed}"

  render_template \
    "${SCRIPT_DIR}/qactden_0mp_parent_window.template.check.ini.xml" \
    "${XML_ROOT}/qactden_0mp_parent_window.${outer_sample_id}.check.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}"

  render_template \
    "${SCRIPT_DIR}/qactden_0mp_child.template.check.ini.xml" \
    "${XML_ROOT}/qactden_0mp_child0.${outer_sample_id}.check.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}"

  render_template \
    "${SCRIPT_DIR}/qactden_0mp_child.template.check.ini.xml" \
    "${XML_ROOT}/qactden_0mp_child1.${outer_sample_id}.check.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "1" \
    "${child1_seed}"

  render_template \
    "${SCRIPT_DIR}/qactden_0mp_outer_sample.template.check.ini.xml" \
    "${XML_ROOT}/qactden_0mp_outer_sample.${outer_sample_id}.check.ini.xml" \
    "${update_no}" \
    "${outer_sample_id}" \
    "0" \
    "${child0_seed}"

  outer_sample_summary_lines="${outer_sample_summary_lines}      <elem>${OUTPUT_ROOT}/${outer_sample_id}/two_level.summary.xml</elem>
"
done

cat > "${XML_ROOT}/qactden_0mp_two_level.check.ini.xml" <<EOF
<?xml version="1.0"?>
<qactden_0mp_corr_check>
  <Input>
    <outer_sample_summaries>
${outer_sample_summary_lines}    </outer_sample_summaries>
  </Input>

  <Output>
    <summary_file>${OUTPUT_ROOT}/two_level.outer_ensemble.summary.xml</summary_file>
    <csv_file>${OUTPUT_ROOT}/two_level.outer_ensemble.csv</csv_file>
  </Output>

  <Reducer>
    <mode>TWO_LEVEL_OUTER_ENSEMBLE</mode>
    <compare_delta_t_parent_set>5 6 7</compare_delta_t_parent_set>
  </Reducer>

  <Checks>
    <min_outer_samples>${OUTER_COUNT}</min_outer_samples>
  </Checks>
</qactden_0mp_corr_check>
EOF

echo "Generated two-level XML bundle under ${XML_ROOT}"
echo "  parent_hmc: ${XML_ROOT}/hmc_parent.0mp_two_level_outer.ini.xml"
echo "  outer_samples: ${OUTER_COUNT} (${FIRST_UPDATE}..${LAST_UPDATE} step ${OUTER_STEP})"
echo "  retained_child_measurements_per_stream: ${EXPECTED_MEASUREMENTS}"
echo "  outer_ensemble_check: ${XML_ROOT}/qactden_0mp_two_level.check.ini.xml"

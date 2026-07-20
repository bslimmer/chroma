#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------

CHROMA_BUILD="${CHROMA_BUILD:-/private/tmp/chroma-build/chroma-localbinarydb-0pppatch}"

SPLIT_BIN="${SPLIT_BIN:-$CHROMA_BUILD/mainprogs/main/gauge_subdomain_split}"
HMC_BIN="${HMC_BIN:-$CHROMA_BUILD/mainprogs/main/hmc}"

# Parent lattice.
PARENT_NROW="${PARENT_NROW:-8 8 8 16}"

# Split direction and cuts.
# For 8^3 x 16, this splits into:
#   child0: global t = [0,1,2,3,4,5,6,7,8]
#   child1: global t = [8,9,10,11,12,13,14,15,0]
TDIR="${TDIR:-3}"
CUT0="${CUT0:-0}"
CUT1="${CUT1:-8}"
FROZEN_WIDTH="${FROZEN_WIDTH:-1}"

# Derived for the symmetric 8^3 x 16, cut 0/8, width 1 case.
# If you change CUT0/CUT1/FROZEN_WIDTH, update CHILD_NROW and frozen intervals.
CHILD_NROW="${CHILD_NROW:-8 8 8 9}"
FROZEN_A_START="${FROZEN_A_START:-0}"
FROZEN_A_END="${FROZEN_A_END:-0}"
FROZEN_B_START="${FROZEN_B_START:-8}"
FROZEN_B_END="${FROZEN_B_END:-8}"

BETA="${BETA:-5.7}"
TAU0="${TAU0:-1.0}"
NSTEPS="${NSTEPS:-20}"

NPROD="${NPROD:-10000}"
NTHISRUN="${NTHISRUN:-1000}"
SAVE_INTERVAL="${SAVE_INTERVAL:-10}"

# Parent configs to process.
# Override like:
#   PARENT_CONFIGS="cfg_10.lime cfg_20.lime" ./split_and_run_child_hmc.sh
PARENT_CONFIGS="${PARENT_CONFIGS:-hmc_parent_cfg_10.lime}"

OUT_ROOT="${OUT_ROOT:-child_hmc_runs}"
# Set to 0 to force the original sequential child-stream behavior.
CHILD_HMC_PARALLEL="${CHILD_HMC_PARALLEL:-1}"

mkdir -p "$OUT_ROOT"

# --------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------

write_split_xml() {
  local parent_cfg="$1"
  local child0_cfg="$2"
  local child1_cfg="$3"
  local sidecar="$4"
  local xml_file="$5"

  cat > "$xml_file" <<XML
<?xml version="1.0"?>
<gauge_subdomain_split>
  <Parent>
    <nrow>$PARENT_NROW</nrow>
    <Cfg>
      <cfg_type>SZINQIO</cfg_type>
      <cfg_file>$parent_cfg</cfg_file>
      <parallel_io>false</parallel_io>
    </Cfg>
  </Parent>

  <Param>
    <t_dir>$TDIR</t_dir>
    <cut0>$CUT0</cut0>
    <cut1>$CUT1</cut1>
    <frozen_width>$FROZEN_WIDTH</frozen_width>
  </Param>

  <Child0>
    <cfg_file>$child0_cfg</cfg_file>
    <volfmt>SINGLEFILE</volfmt>
  </Child0>

  <Child1>
    <cfg_file>$child1_cfg</cfg_file>
    <volfmt>SINGLEFILE</volfmt>
  </Child1>

  <Sidecar>
    <file>$sidecar</file>
  </Sidecar>
</gauge_subdomain_split>
XML
}

write_child_hmc_xml() {
  local child_cfg="$1"
  local save_prefix="$2"
  local seed0="$3"
  local xml_file="$4"

  cat > "$xml_file" <<XML
<?xml version="1.0"?>
<Params>
  <MCControl>
    <Cfg>
      <cfg_type>SCIDAC</cfg_type>
      <cfg_file>$child_cfg</cfg_file>
      <parallel_io>false</parallel_io>
    </Cfg>

    <RNG>
      <Seed>
        <elem>$seed0</elem>
        <elem>0</elem>
        <elem>0</elem>
        <elem>0</elem>
      </Seed>
    </RNG>

    <StartUpdateNum>0</StartUpdateNum>
    <NWarmUpUpdates>0</NWarmUpUpdates>
    <NProductionUpdates>$NPROD</NProductionUpdates>
    <NUpdatesThisRun>$NTHISRUN</NUpdatesThisRun>

    <SaveInterval>$SAVE_INTERVAL</SaveInterval>
    <SavePrefix>$save_prefix</SavePrefix>
    <SaveVolfmt>SINGLEFILE</SaveVolfmt>

    <ReproCheckP>false</ReproCheckP>
    <ReverseCheckP>false</ReverseCheckP>
    <MonitorForces>false</MonitorForces>
  </MCControl>

  <HMCTrj>
    <Monomials>
      <elem>
        <Name>GAUGE_MONOMIAL</Name>

        <GaugeAction>
          <Name>WILSON_GAUGEACT</Name>
          <beta>$BETA</beta>

          <GaugeBC>
            <Name>TEMPORAL_ZONE_GAUGEBC</Name>
            <t_dir>$TDIR</t_dir>
            <zero_intervals>
              <elem>
                <t_start>$FROZEN_A_START</t_start>
                <t_end>$FROZEN_A_END</t_end>
              </elem>
              <elem>
                <t_start>$FROZEN_B_START</t_start>
                <t_end>$FROZEN_B_END</t_end>
              </elem>
            </zero_intervals>
          </GaugeBC>
        </GaugeAction>

        <NamedObject>
          <monomial_id>gauge</monomial_id>
        </NamedObject>
      </elem>
    </Monomials>

    <Hamiltonian>
      <monomial_ids>
        <elem>gauge</elem>
      </monomial_ids>
    </Hamiltonian>

    <MDIntegrator>
      <tau0>$TAU0</tau0>

      <Integrator>
        <Name>LCM_STS_LEAPFROG</Name>
        <n_steps>$NSTEPS</n_steps>
        <monomial_ids>
          <elem>gauge</elem>
        </monomial_ids>
      </Integrator>
    </MDIntegrator>

    <nrow>$CHILD_NROW</nrow>
  </HMCTrj>
</Params>
XML
}

run_child_hmc() {
  local child_label="$1"
  local child_hmc_xml="$2"
  local child_out_xml="$3"
  local child_log_xml="$4"

  echo "Running ${child_label} HMC with frozen boundary slices..."
  "$HMC_BIN" \
    -i "$child_hmc_xml" \
    -o "$child_out_xml" \
    -l "$child_log_xml"
}

# --------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------

for parent_cfg in $PARENT_CONFIGS; do
  if [ ! -f "$parent_cfg" ]; then
    echo "ERROR: parent config not found: $parent_cfg"
    exit 1
  fi

  tag="$(basename "$parent_cfg")"
  tag="${tag%.lime}"
  tag="${tag%.scidac}"

  run_dir="$OUT_ROOT/$tag"
  mkdir -p "$run_dir"

  parent_abs="$(cd "$(dirname "$parent_cfg")" && pwd)/$(basename "$parent_cfg")"

  child0_cfg="$run_dir/child0.${tag}.scidac"
  child1_cfg="$run_dir/child1.${tag}.scidac"
  sidecar="$run_dir/${tag}.sidecar.xml"

  split_xml="$run_dir/split.${tag}.ini.xml"
  split_out="$run_dir/split.${tag}.out.xml"
  split_log="$run_dir/split.${tag}.log.xml"

  child0_hmc_xml="$run_dir/hmc_child0.${tag}.ini.xml"
  child1_hmc_xml="$run_dir/hmc_child1.${tag}.ini.xml"
  child0_out_xml="$run_dir/hmc_child0.${tag}.out.xml"
  child1_out_xml="$run_dir/hmc_child1.${tag}.out.xml"
  child0_log_xml="$run_dir/hmc_child0.${tag}.log.xml"
  child1_log_xml="$run_dir/hmc_child1.${tag}.log.xml"

  child0_prefix="$run_dir/child0_hmc.${tag}"
  child1_prefix="$run_dir/child1_hmc.${tag}"

  echo "============================================================"
  echo "Parent config: $parent_cfg"
  echo "Run dir:       $run_dir"
  echo "============================================================"

  write_split_xml "$parent_abs" "$child0_cfg" "$child1_cfg" "$sidecar" "$split_xml"

  echo "Splitting parent into child configs..."
  "$SPLIT_BIN" \
    -i "$split_xml" \
    -o "$split_out" \
    -l "$split_log"

  echo "Writing child HMC XML..."
  write_child_hmc_xml "$child0_cfg" "$child0_prefix" 23 "$child0_hmc_xml"
  write_child_hmc_xml "$child1_cfg" "$child1_prefix" 23 "$child1_hmc_xml"

  if [ "$CHILD_HMC_PARALLEL" = "1" ]; then
    echo "Launching child HMC streams in parallel..."
    run_child_hmc "child0" "$child0_hmc_xml" "$child0_out_xml" "$child0_log_xml" &
    child0_pid=$!
    run_child_hmc "child1" "$child1_hmc_xml" "$child1_out_xml" "$child1_log_xml" &
    child1_pid=$!

    child0_status=0
    child1_status=0
    wait "$child0_pid" || child0_status=$?
    wait "$child1_pid" || child1_status=$?

    if [ "$child0_status" -ne 0 ] || [ "$child1_status" -ne 0 ]; then
      echo "ERROR: child HMC failed (child0=${child0_status}, child1=${child1_status})"
      exit 1
    fi
  else
    run_child_hmc "child0" "$child0_hmc_xml" "$child0_out_xml" "$child0_log_xml"
    run_child_hmc "child1" "$child1_hmc_xml" "$child1_out_xml" "$child1_log_xml"
  fi

  echo "Done with $parent_cfg"
  echo
done

echo "All parent configs processed."

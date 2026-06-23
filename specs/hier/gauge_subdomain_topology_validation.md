# Gauge Subdomain Topology Validation Spec

## Goal

Define a standalone end-to-end validation workflow for the split child-lattice
gauge-only HMC use case using raw topological charge density measurements on
the child lattices.

The validation target is:

1. start from one cold parent lattice
2. run short whole-lattice gauge-only HMC
3. split the warmed parent into two equal child lattices
4. run independent child gauge-only HMC with frozen duplicated temporal
   boundary blocks
5. measure a child-local timeseries of raw topological charge restricted to the
   child interior
6. compare the two child timeseries for statistical consistency

This is a workflow-validation spec. It is not yet a production topology study
and it is not yet a `0-+` glueball spectroscopy spec.

## Scope

This first version intentionally stays narrow:

- gauge-only HMC only
- no fermion monomials
- no SMD path
- raw topology only
- no Wilson-flow smoothing in the first version
- one whole-lattice warmup run before splitting
- two independent child runs with different RNG seeds
- one child-local topology timeseries comparison
- one explicit frozen-link invariance check
- small-lattice first instance, with parameters kept configurable for later
  large-lattice runs

This first version does not require:

- stitched-parent physics comparison
- file-backed topology reducers
- Wilson-flowed topology
- long-run autocorrelation analysis
- a dedicated analysis helper in the repository

## Relationship To Existing Specs

This spec is standalone, but it builds on the following existing contracts:

- split/stitch geometry and sidecar metadata:
  `specs/hier/gauge_subdomain_split.md`
- child-run force freezing:
  `specs/hier/temporal_zone_gaugebc.md`
- child-run refreshed-momentum freezing:
  `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`

This spec is intentionally separate from the plaquette-based validation in:

- `specs/hier/gauge_subdomain_gauge_hmc_validation.md`

The intent is to make this document reusable later as a template for other
child-local observable studies.

## Current Repo Facts

- `mainprogs/main/hmc.cc` accepts `InlineMeasurements` under `MCControl` and
  emits their outputs in `InlineObservables`.
- `QTOP_NAIVE` already exists as an inline measurement of one total gluonic
  topological charge:
  - XML/factory name: `QTOP_NAIVE`
  - implementation: `lib/meas/inline/glue/inline_qnaive.cc`
- `QACTDEN` already exists as an inline measurement of the local raw action
  density and local raw topological charge density:
  - XML/factory name: `QACTDEN`
  - implementation: `lib/meas/inline/glue/inline_qactden.cc`
  - core density routine: `lib/meas/glue/qactden.cc`
- `QTOP_NAIVE` is not sufficient for this study because it only returns one
  whole-child total `Q`, whereas this validation wants a quantity restricted to
  the non-frozen child interior.
- `QACTDEN` writes a full `LatticeReal` field named `naiveTopCharge` into XML.
  That is acceptable for the intended small-lattice first version, but it is
  not an ideal transport format for large-lattice production runs.
- `QACTDEN` is not ultralocal. The implementation in `lib/meas/glue/qactden.cc`
  uses forward and backward nearest-neighbor shifts in every plane. The choice
  of which local time slices count as the "interior" therefore matters.
- `TEMPORAL_ZONE_GAUGEBC` plus the HMC gauge-monomial momentum BC autodiscovery
  together provide the intended child-run freezing semantics on the duplicated
  boundary intervals.

## Decisions Fixed For This Draft

The following choices are fixed for this draft based on current discussion:

- use raw topology in the first version
- do not use Wilson flow in the first version
- compare a topology timeseries, not only one final-state value
- use different RNG seeds in the two child runs
- keep trajectory counts configurable
- use `100` whole-lattice warmup trajectories in the suggested first instance
- use `100` child trajectories in the suggested first instance
- use measurement cadence `5` trajectories
- discard the first `20` child trajectories in the first-instance comparison
- use postprocessing of emitted `QACTDEN` outputs rather than requiring a
  dedicated reducer/helper inside the repository

## Confirmed First-Version Conventions

The following conventions are fixed for the first version:

1. Cold-start convention:
   - the parent run starts from a cold gauge field
   - state this explicitly as `cfg_type = UNIT`

2. Child-interior definition for `QACTDEN`:
   - use `RootInterior`
   - include all root sites whose child-local time coordinate is not in a
     frozen slice
   - this deliberately favors a broader first-version signal on the small child
     lattice, even though the `QACTDEN` stencil reaches one link beyond the
     root site

For the existing symmetric `4 4 4 5` child geometry with frozen local slices
`[0,0]` and `[4,4]`, this gives:

- `active_topology_slices = {1, 2, 3}`

## Chosen Flexible Geometry Template

Keep the validation geometry parameterized by:

- `parent_nrow`
- `t_dir`
- `cut0`
- `cut1`
- `frozen_width`
- parent HMC trajectory count
- child HMC trajectory count
- child measurement cadence

The first-instance geometry should remain aligned with the existing symmetric
split workflow so the comparison is easy to reason about:

- `parent_nrow = 4 4 4 8`
- `t_dir = 3`
- `cut0 = 0`
- `cut1 = 4`
- `frozen_width = 1`

Derived child geometry:

- `child0_nrow = 4 4 4 5`
- `child1_nrow = 4 4 4 5`
- `child0_frozen_local_intervals = [0,0] and [4,4]`
- `child1_frozen_local_intervals = [0,0] and [4,4]`

This spec should keep these values as configurable inputs so the same workflow
can later be reused on larger lattices.

For this first-instance geometry:

- `active_topology_slices = {1, 2, 3}`

## Planned Test Artifacts

Add a small topology-validation bundle under `tests/gauge_subdomain_split/`.

Suggested first-version inputs:

- `hmc_parent.topology_warmup.ini.xml`
  - whole-lattice gauge-only HMC input
  - cold parent start with `cfg_type = UNIT`
  - short warmup run
  - saves one warmed parent config for splitting
- `gauge_subdomain_split.topology_validation.ini.xml`
  - split input
  - reads the warmed parent config
  - writes `child0` and `child1` configs plus sidecar XML
- `hmc_child0.temporal_zone_qactden.ini.xml`
  - gauge-only child HMC input for `child0`
  - includes inline `QACTDEN`
- `hmc_child1.temporal_zone_qactden.ini.xml`
  - same as `child0` except for filenames and RNG seed

Generated outputs from these inputs are workspace artifacts and should not be
committed. For local runs in the source checkout, create a dedicated `cfgs/`
subdirectory and route generated configs, restart XML, split sidecars, and
command-line XML outputs there rather than the repository root.

The first version does not require separate one-shot `chroma` measurement
inputs because `QACTDEN` can run inline during child HMC.

## Parent Warmup Requirements

The parent warmup step should use `mainprogs/main/hmc.cc` with:

- one `GAUGE_MONOMIAL`
- a simple gauge action such as `WILSON_GAUGEACT`
- ordinary periodic gauge BC
- `cfg_type = UNIT`
- no fermion monomials
- one saved final warmed parent config

For the first-instance run:

- use `100` trajectories
- keep the integrator and gauge-action parameters configurable
- use the same gauge action family intended for the child runs so the split
  study does not change more ingredients than necessary

This step is only intended to move away from the exact cold start before the
split.

## Child HMC Requirements

The child evolution step should use `mainprogs/main/hmc.cc` with:

- one `GAUGE_MONOMIAL`
- a simple gauge action such as `WILSON_GAUGEACT`
- `TEMPORAL_ZONE_GAUGEBC` under the gauge action
- `zero_intervals` taken exactly from the child-local frozen intervals in the
  split sidecar
- child-local `nrow`
- a saved final child config

For the first-instance run:

- use `100` trajectories per child
- use different RNG seeds for `child0` and `child1`
- keep `beta`, `tau0`, integrator type, and `n_steps` the same between the two
  children
- use inline `QACTDEN` with frequency `5`
- keep reverse-check noise disabled unless it is specifically being debugged

Required implementation assumption:

- the HMC trajectory implements
  `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`

## Inline Measurement Plan

The first-version child HMC inputs should include:

```xml
<InlineMeasurements>
  <elem>
    <Name>QACTDEN</Name>
    <Frequency>5</Frequency>
    <Param>
      <version>1</version>
    </Param>
    <NamedObject>
      <gauge_id>default_gauge_field</gauge_id>
    </NamedObject>
  </elem>
</InlineMeasurements>
```

This produces a `QActDen` block in the child HMC XML output at the chosen
measurement cadence.

The first version intentionally accepts the cost of XML field output because
the target lattice is small. A larger-lattice follow-up may want a different
transport path for the density field.

## Observable Definition

For each child measurement time `n`, let:

- `q_n(x)` be the local field written as `naiveTopCharge` by `QACTDEN`
- `active_topology_slices` be the set of child-local root time slices that are
  not in any frozen interval

Define the interior topological charge timeseries point:

```text
Q_int^(child)(n) = sum_{x : t_local(x) in active_topology_slices} q_n(x)
```

Properties:

- this is a child-local partial topological charge, not a full-lattice
  topological charge
- it is intentionally restricted away from the duplicated frozen boundary
  slices
- in the raw first version it is not expected to be quantized

## Timeseries Definition

For the first-instance run:

- child trajectory count: `100`
- measurement cadence: every `5` trajectories
- raw measurement indices:
  `{5, 10, 15, ..., 100}`
- discard the first four measurements
  `{5, 10, 15, 20}`
- retained measurements:
  `{25, 30, 35, ..., 100}`
- retained count per child: `M = 16`

This first version treats the retained sequence as a short comparison
timeseries, not as a full production ensemble.

## Statistical Comparison Rule

For each child `c in {0,1}`, compute over the retained `M` measurements:

```text
mu_c = (1 / M) * sum_{m=1..M} Q_int^(c)(m)
```

Use the ordinary sample standard deviation:

```text
s_c^2 = (1 / (M - 1)) * sum_{m=1..M} (Q_int^(c)(m) - mu_c)^2
```

and define the first-version standard error estimate:

```text
SE_c = s_c / sqrt(M)
```

Then define:

```text
Delta_mu = |mu_0 - mu_1|
SE_Delta = sqrt(SE_0^2 + SE_1^2)
```

First-version acceptance criterion:

```text
Delta_mu <= 3 * SE_Delta
```

If `SE_Delta = 0`, interpret the comparison as:

- pass if `mu_0 == mu_1`
- fail otherwise

This is intentionally a pragmatic short-run validation criterion. A later
stronger study may replace the unblocked `SE` estimate with a blocked or
autocorrelation-aware treatment.

## Required Checks

The first validation should require all of the following:

1. The parent warmup HMC run completes successfully and writes one warmed
   parent config.
2. The split tool completes successfully and writes the expected sidecar.
3. The sidecar reports the expected child-local frozen intervals.
4. Both child HMC runs complete successfully.
5. Both child HMC runs emit the expected number of `QACTDEN` measurements.
6. Frozen-link invariance:
   for each child, all link directions rooted on the frozen local slices in the
   final evolved child config must match the initial split child config exactly.
7. Postprocessing successfully reduces each emitted `naiveTopCharge` field to
   one `Q_int` value per retained child measurement.
8. The child means satisfy:
   `|mu_0 - mu_1| <= 3 * sqrt(SE_0^2 + SE_1^2)`.

## Interpretation Notes

- A successful result here is a consistency check, not yet a physics claim.
- For a short raw-topology small-lattice study, both child means may remain
  close to zero. That is acceptable for this first validation.
- If the retained timeseries is numerically trivial and therefore not very
  informative, the natural next steps are:
  - run longer
  - increase the measurement window
  - move to a larger lattice
  - add Wilson flow in a follow-up spec

## Workflow

The intended first validation workflow is:

1. run whole-lattice gauge-only HMC from a cold parent start
2. save one warmed parent config
3. split the warmed parent into `child0`, `child1`, and one sidecar XML file
4. preserve the written child configs as the initial child reference states
5. run child gauge-only HMC on `child0` with inline `QACTDEN`
6. run child gauge-only HMC on `child1` with inline `QACTDEN`
7. postprocess the two child HMC XML logs to extract each emitted
   `naiveTopCharge` field
8. reduce each field to one `Q_int` value using `active_topology_slices`
9. discard the first `20` trajectories worth of measurements
10. compute `mu_0`, `mu_1`, `SE_0`, `SE_1`
11. compare the child means with the `3 * SE_Delta` criterion
12. separately compare the final child configs against the initial child
    configs on the frozen local slices

## Non-Goals For The First Version

- no Wilson-flowed topology in the first version
- no `0-+` correlator or spectroscopy extraction
- no stitched-parent comparison
- no production-quality uncertainty analysis
- no dedicated reducer executable
- no requirement that the first implementation run efficiently on large
  lattices

## Scaling Notes

- This draft is intentionally suitable for the small symmetric child geometry
  already used in local validation work.
- The geometry, trajectory counts, and measurement cadence should remain
  configurable so the same workflow can later be moved to larger lattices.
- Because `QACTDEN` currently writes the full density field to XML, a future
  supercomputer-scale version may want either:
  - a file-backed output path for the density field, or
  - a lightweight reducer that emits only the interior-integrated `Q_int`
    timeseries
- A natural follow-up document is a Wilson-flowed version of this same
  topology-validation workflow.

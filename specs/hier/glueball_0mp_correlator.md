# 0-+ Glueball Correlator Spec

## Goal

Define a first implementation and validation plan for a zero-momentum `0-+`
gluonic correlator built by reusing the existing `QACTDEN` inline measurement.

The immediate target is not a believable mass determination on local toy
lattices. The immediate target is:

1. emit the local pseudoscalar gluonic density with existing inline machinery
2. reduce that field offline into a per-configuration Euclidean-time correlator
3. validate the reduction on a full periodic lattice first
4. run the same reduction on split child lattices with an interior-only
   measurement window
5. compare the two child correlator ensembles for statistical consistency

This is therefore an implementation-and-validation spec. It is not yet a
production spectroscopy campaign spec.

## Scope

This first version intentionally stays narrow:

- gauge-only HMC only
- no fermion monomials
- no SMD path
- zero momentum only
- raw `QACTDEN`-based topology only
- no new file under `lib/meas/inline/glue/` in the first version
- no Wilson-flow smoothing in the first version
- no variational basis of multiple `0-+` operators
- no believable mass fit requirement on local test volumes
- full-lattice validation first, split-child comparison second

This first version does not require:

- a production-quality `0-+` mass estimate
- a stitched-parent spectroscopy comparison
- a flowed or cooled topology measurement
- a generalized eigenvalue analysis
- long-run autocorrelation analysis
- large-lattice optimization

## Relationship To Existing Specs

This spec is standalone, but it builds on the following existing work:

- split/stitch geometry and sidecar metadata:
  `specs/hier/gauge_subdomain_split.md`
- child-run force freezing:
  `specs/hier/temporal_zone_gaugebc.md`
- child-run refreshed-momentum freezing:
  `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`
- child topology-validation workflow using raw `QACTDEN`:
  `specs/hier/gauge_subdomain_topology_validation.md`

This spec is intentionally separate from:

- `specs/hier/gauge_subdomain_gauge_hmc_validation.md`
  - simple plaquette-based child validation

## Current Repo Facts

- `mainprogs/main/hmc.cc` accepts `InlineMeasurements` under `MCControl` and
  emits their outputs in `InlineObservables`.
- `QACTDEN` already exists and writes:
  - a local gluonic action-density field
  - a local raw topological-charge-density field named `naiveTopCharge`
  - implementation: `lib/meas/inline/glue/inline_qactden.cc`
  - core density routine: `lib/meas/glue/qactden.cc`
- The local `naiveTopCharge` field from `QACTDEN` is the current repository
  object most directly related to the `0-+` channel.
- `QACTDEN` is not ultralocal. Its core routine uses forward and backward
  nearest-neighbor shifts in every plane, so the topological-density value at
  one rooted timeslice depends on adjacent timeslices.
- `PLAQ_DENSITY` writes plaquette action densities to a separate XML file. It
  is useful for plaquette-based validation, but it does not contain the
  pseudoscalar `F \tilde F` information needed for a `0-+` correlator.
- `QTOP_NAIVE` computes one total topological charge, not a timeslice operator
  or correlator:
  - implementation: `lib/meas/inline/glue/inline_qnaive.cc`
- The old file `lib/meas/glue/qtopcor.cc` is unfinished and computes a radial
  four-dimensional correlation function, not the timeslice correlator needed
  here.
- Existing glueball code in `lib/meas/glue/gluecor.cc` constructs `0++`,
  `1+-`, and `2++` correlators from fuzzy plaquette operators. It does not
  implement the pseudoscalar `0-+` channel.
- `GLUEBALL_OPS` already shows repository precedent for:
  - zero-momentum projection via `SftMom`
  - optional gauge-link smearing before gluonic-operator construction
  - field-strength-based gluonic operator construction
  - inline measurement output written through HMC XML plus sidecar DB files
- `GLUEBALL_OPS` does not currently construct the `QACTDEN` topological-density
  operator or a `0-+` correlator from it.
- In the split workflow, `TEMPORAL_ZONE_GAUGEBC` freezes forces and refreshed
  momenta on selected child-local temporal intervals but does not modify the
  stored gauge links themselves.

## Physics Target And Observable Definition

### Local Operator

For this draft, define the local pseudoscalar gluonic density by reusing the
`QACTDEN` local field:

- `q(x) = naiveTopCharge(x)`

This is the lattice field that will be used as the first raw `0-+` operator.

### Zero-Momentum Timeslice Operator

On one gauge configuration, define the zero-momentum timeslice operator:

- `O_q(t) = sum_{x spatial at timeslice t} q(x, t)`

For the first version, the primary stored operator should be the spatial sum,
not the spatial average. The reducer should also record the spatial volume so
any later normalization change remains recoverable offline.

### Full-Lattice Periodic Correlator

For a periodic full lattice with temporal extent `T`, define the
per-configuration correlator:

- `C_cfg(dt) = (1 / T) * sum_{t0=0}^{T-1} O_q(t0 + dt mod T) * O_q(t0)`

for `dt = 0, 1, ..., dt_max`, with `dt_max <= floor(T/2)`.

### Child-Lattice Interior No-Wrap Correlator

For a split child lattice, the correlator must not assume periodic source/sink
wrapping across the duplicated frozen boundary zones.

Define:

- `frozen_intervals`
  - child-local temporal intervals held fixed by the split-child workflow
- `support_guard`
  - extra number of child-local timeslices excluded on each side of each frozen
    interval because the chosen operator reaches beyond its rooted timeslice

For the first `QACTDEN`-based operator:

- `support_guard = 1`

because the density construction uses forward and backward nearest-neighbor
shifts.

Define the child-local safe timeslice set:

- `safe_slices`
  - all child-local timeslices not in a frozen interval and not within
    `support_guard` of a frozen interval

Then define, for one child configuration and one separation `dt`,

- `source_slices(dt) = { t0 in safe_slices | t0 + dt is also in safe_slices }`

with no wrap-around.

The per-configuration child correlator is:

- `C_cfg(dt) = (1 / N_src(dt)) * sum_{t0 in source_slices(dt)} O_q(t0 + dt) * O_q(t0)`

where `N_src(dt) = |source_slices(dt)|`.

Only separations with `N_src(dt) > 0` are defined.

## Current-Geometry Implication

The current small split-child geometry used for topology validation is:

- parent `4 4 4 8`
- child `4 4 4 5`
- frozen child-local slices `{0, 4}`

For the first `QACTDEN`-based `support_guard = 1` operator, this leaves only
one safe child-local timeslice:

- `safe_slices = {2}`

Therefore:

- the current `4 4 4 5` child geometry is acceptable only for `dt = 0`
  plumbing checks
- it is not acceptable for a nonzero-separation child-correlator validation

The split-child comparison phase in this spec therefore needs a slightly larger
temporal geometry than the existing topology-timeseries test.

## Implementation Strategy

### First-Version Reuse Plan

The first implementation should reuse the existing `QACTDEN` inline
measurement directly.

The measurement path is:

1. run `QACTDEN` inline during HMC
2. read `QActDen/naiveTopCharge` from the emitted HMC XML
3. reduce that local field offline into:
   - `O_q(t)`
   - `C_cfg(dt)`
   - the metadata needed to interpret the measurement window

The first version should not add:

- a new inline measurement in `lib/meas/inline/glue/`
- a new library measurement helper in `lib/meas/glue/`

unless later scaling pressure justifies it.

### Preferred First Reducer Form

The preferred first reducer is a focused checker executable under
`mainprogs/tests/`, because it keeps the workflow reproducible inside the
repository without creating new inline-measurement surface area.

Suggested first-version checker:

- `mainprogs/tests/t_qactden_0mp_corr.cc`

This checker should:

1. read one HMC XML file
2. locate every `QActDen` block
3. reduce each `naiveTopCharge` field into the timeslice operator `O_q(t)`
4. construct either:
   - the periodic full-lattice correlator, or
   - the child interior no-wrap correlator
5. write a concise XML or text summary for regression checking
6. optionally emit CSV for human inspection

This checker may later be refactored into shared library code if the same
reduction logic is needed in multiple places. The first version does not
require that refactor.

### Reducer Inputs

The reducer/checker should accept:

- HMC XML input file
- `decay_dir`
- `max_dt`
- `mode`
  - `PERIODIC_ALL_T`
  - `INTERIOR_NOWRAP`
- `support_guard`
  - first-version default `1`
- `blocked_intervals`
  - required only for `INTERIOR_NOWRAP`

For ensemble comparison mode, it should also accept:

- discard count
- selected `dt` set to compare

### Reducer Outputs

For each `QActDen` block, the reducer should produce:

- `update_no`
- `timeslice_operator`
  - one `O_q(t)` value for each child-local or full-lattice timeslice
- `correlator`
  - one `C_cfg(dt)` value for each emitted separation
- `num_sources(dt)`
- `source_slices(dt)`
- `safe_slices`
  - for child no-wrap mode

The first version does not need to write these back into a Chroma inline XML
stream. Standalone checker output is sufficient.

## Validation Ladder

Validation should happen in two phases.

### Phase 1: Full Periodic Lattice Correctness

This phase validates that the `QACTDEN` reduction is mathematically correct on
ordinary periodic lattices before any split-child boundary complications are
introduced.

#### Suggested Artifacts

- `tests/glueball_0mp/hmc_full_lattice.qactden_0mp.ini.xml`
  - periodic gauge-only HMC
  - includes inline `QACTDEN`
- `tests/glueball_0mp/qactden_0mp_corr.full_lattice.check.ini.xml`
  - checker input
- `mainprogs/tests/t_qactden_0mp_corr.cc`
  - focused checker executable

Optional trivial smoke:

- `tests/glueball_0mp/measure_unit_qactden_0mp.ini.xml`
  - one-shot unit-field measurement input

#### Suggested First-Instance Geometry

- `nrow = 4 4 4 8`
- `t_dir = 3`
- periodic gauge BC
- short gauge-only HMC run
- measurement cadence `5`

The geometry can stay small because this phase checks construction correctness,
not spectroscopy quality.

#### Required Checks

1. `QACTDEN` registers successfully and runs under `hmc`.
2. The number of emitted `QActDen` blocks matches the requested frequency.
3. For every emitted update, the reducer emits exactly one `O_q(t)` value per
   timeslice.
4. For every emitted update, the reducer emits `num_sources(dt) = T` for every
   periodic separation `dt`.
5. For every emitted update, the reducer constructs the periodic wrap-around
   correlator using all `t0 mod T`.
6. For a trivial unit-field smoke input, every `O_q(t)` and every `C_cfg(dt)`
   should vanish within floating-point tolerance.

Optional development cross-check:

- compare the checker output against a lightweight independent script during
  early implementation, but this is not a required repo artifact.

This phase should pass before any split-child comparison work is considered
complete.

### Phase 2: Split Child-Lattice Comparison

This phase validates that the same `QACTDEN` reduction runs correctly on split
child lattices when the correlator is restricted to the non-frozen interior
window.

#### Recommended First Child-Correlator Geometry

Use a slightly larger temporal parent than the current topology-timeseries
validation:

- `parent_nrow = 4 4 4 12`
- `t_dir = 3`
- `cut0 = 0`
- `cut1 = 6`
- `frozen_width = 1`

Derived child geometry:

- `child_nrow = 4 4 4 7`
- `child_frozen_local_intervals = [0,0] and [6,6]`
- `child_active_root_slices = {1, 2, 3, 4, 5}`
- `support_guard = 1`
- `safe_slices = {2, 3, 4}`

This gives:

- `num_sources(0) = 3`
- `num_sources(1) = 2`
- `num_sources(2) = 1`

For the first comparison pass, define:

- `compare_dt_set = { dt | num_sources(dt) >= 2 }`

which gives:

- `compare_dt_set = {0, 1}`

`dt = 2` should still be emitted for plumbing visibility, but it should not be
used in the first pass/fail child-child comparison criterion.

#### Suggested Artifacts

Under `tests/gauge_subdomain_split/`:

- `hmc_parent.0mp_corr_warmup.ini.xml`
  - periodic parent warmup
- `gauge_subdomain_split.0mp_corr_validation.ini.xml`
  - split input
- `hmc_child0.temporal_zone_qactden_0mp.ini.xml`
  - child0 gauge-only HMC
  - includes inline `QACTDEN`
- `hmc_child1.temporal_zone_qactden_0mp.ini.xml`
  - child1 gauge-only HMC
  - includes inline `QACTDEN`
- `gauge_subdomain_qactden_0mp_corr.check.ini.xml`
  - checker input for child comparison

The same `t_qactden_0mp_corr` checker may validate both the full-lattice and
child-lattice workflows if its input format is kept generic enough.

#### Suggested First-Instance Run Controls

Reuse the first-pass control philosophy from the topology-timeseries test:

- `100` parent warmup trajectories
- `100` child trajectories per child
- measurement cadence `5`
- discard the first `20` child trajectories in the ensemble comparison
- use different RNG seeds in the two child runs

These values should remain configurable.

#### Child-Phase Required Checks

1. The full-lattice correctness phase has already passed.
2. Parent warmup, split, and both child HMC runs complete successfully.
3. Each child HMC run emits the expected number of `QActDen` blocks.
4. The reducer derives child `safe_slices` and `num_sources(dt)` values that
   match the geometry implied by:
   - the child-local frozen intervals
   - `support_guard = 1`
   - no wrap-around
5. For every retained child measurement, the reducer emits:
   - one `O_q(t)` value per child-local timeslice
   - one `C_cfg(dt)` value for every defined separation
6. For every `dt` in `compare_dt_set`, the retained child correlator means are
   statistically consistent between `child0` and `child1`.

For the first-pass ensemble comparison, define:

- `mu_c(dt)`
  - retained-measurement mean of child `c` at separation `dt`
- `SE_c(dt)`
  - sample standard deviation divided by `sqrt(M)` for the retained series
- `Delta_mu(dt) = abs(mu_0(dt) - mu_1(dt))`
- `SE_Delta(dt) = sqrt(SE_0(dt)^2 + SE_1(dt)^2)`

Pass the child-child comparison if, for every `dt` in `compare_dt_set`,

- `Delta_mu(dt) <= 3 * SE_Delta(dt)`

This is intentionally the same first-pass style used in the raw-topology
timeseries validation. A later production workflow may replace it with a
blocking or autocorrelation-aware treatment.

## Non-Goals For The First Version

- no claim of a believable `0-+` mass on the local test lattices
- no new inline measurement under `lib/meas/inline/glue/`
- no Wilson-flowed or cooled operator in the first pass
- no variational operator basis
- no GEVP analysis
- no stitched-parent mass comparison
- no attempt to retrofit the unfinished radial `qtopcor` path
- no requirement that the first implementation be efficient enough for
  supercomputer-scale `QACTDEN` field-output volumes

## Follow-Up Path To Actual Mass Extraction

Once the first implementation and validation pass works correctly, the natural
follow-up path is:

1. add a flowed or otherwise smoothed `0-+` operator variant
2. if needed, replace full-field XML postprocessing with a lighter-weight
   reduced-output path
3. save larger periodic and split-child ensembles
4. perform ensemble averaging and uncertainty analysis offline
5. extract effective masses and fit windows from the resulting correlators
6. only then treat the split-child observable as a candidate spectroscopy
   measurement rather than a code-path validation

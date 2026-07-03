# Gauge-Subdomain Two-Level 0++ Workflow Spec

## Goal

Define a first two-level measurement workflow that combines:

1. parent-lattice gauge evolution
2. split child-lattice gauge-only HMC with frozen duplicated temporal boundary
   blocks
3. subdomain-localized `0++` timeslice operators built from existing spatial
   Wilson-loop / plaquette glueball code

into one nested estimator for a cross-subdomain scalar glueball correlator.

The immediate target is not a production spectroscopy campaign. The immediate
target is:

1. generate a sequence of level-0 parent boundary samples
2. split each parent sample into two child lattices plus sidecar metadata
3. generate independent level-1 child ensembles at fixed boundary data
4. reduce each child ensemble to conditional means of the localized scalar
   timeslice operator
5. form a factorized cross-subdomain `0++` correlator from those conditional
   means
6. validate that estimator against a direct parent-window measurement on a
   small lattice
7. stage the workflow toward a first even-split `8^4` campaign with `40`
   retained parent samples and `10` retained child measurements per child
   domain per retained parent sample

This is a workflow-and-estimator spec. It is not yet a performance-tuning
document and it is not yet a mass-fit campaign spec.

## Scope

This first version intentionally stays narrow:

- gauge-only HMC only
- two temporal child domains only
- existing spatial-plaquette / Wilson-loop `0++` operator only
- one explicit blocking level per run
- zero momentum only
- one two-level estimator built from operators in opposite child interiors
- one small-lattice validation path against a direct parent-window observable
- current split sidecar geometry contract only

This first version does not require:

- fermion monomials
- SMD support
- Wilson flow or cooling in the first pass
- a variational basis of multiple `0++` operators
- stitched-parent postprocessing
- full reconstruction of the translationally averaged periodic correlator
- production-scale cost optimization
- an automatic workflow launcher
- more than two child subdomains

## Relationship To Existing Specs

This spec sits on top of the following existing work:

- split/stitch geometry and persisted sidecar metadata:
  `specs/hier/gauge_subdomain_split.md`
- child-run force freezing:
  `specs/hier/temporal_zone_gaugebc.md`
- child-run refreshed-momentum freezing:
  `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`
- structural template for the level-0 / level-1 nesting and validation ladder:
  `specs/hier/gauge_subdomain_two_level_0mp_workflow.md`

This spec also relies on existing `0++` glueball measurement code already in
the repository:

- scalar-glueball correlator construction from spatial plaquettes:
  `lib/meas/glue/gluecor.cc`
- driver for unblocked and blocked glueball measurements:
  `lib/meas/glue/fuzglue.cc`
- existing inline Wilson-loop and glueball-operator surfaces:
  `lib/meas/inline/glue/{inline_wilslp.cc,inline_fuzwilp.cc,inline_glueball_ops.cc}`

The conceptual pattern is the same standard two-level one used in the
multilevel lattice literature: for each level-0 boundary sample, average
subdomain-local observables independently at level 1, then average the
product of those level-1 means over the level-0 ensemble.

This document intentionally follows the same shape as the `0-+` workflow spec,
but it differs in one important operator-level detail:

- the first `0++` operator is local in Euclidean time, so the first-pass
  temporal support guard is `0`
- unlike `QACTDEN`, the repository does not yet expose that `0++` timeslice
  operator in a turnkey per-update HMC output surface, so the first
  measurement surface should be a thin wrapper around the existing
  scalar-glueball code

## Current Repo Facts

- `lib/meas/glue/gluecor.cc` already constructs `0++`, `1+-`, and `2++`
  correlators from purely spatial plaquettes on each timeslice.
- Its internal scalar timeslice operator is `op0[t]`, which is exactly the
  object most closely aligned with the first two-level `0++` estimator.
- The `0++` construction in `gluecor.cc` uses only spatial planes orthogonal
  to the decay direction.
- Therefore the rooted operator on timeslice `t` does not reach onto adjacent
  timeslices, even when a blocked spatial plaquette level is used.
- `lib/meas/glue/fuzglue.cc` already drives that same scalar-operator family
  across unblocked and blocked spatial levels.
- Within that existing infrastructure, the first blocked level
  `bl_level_selected = 1` corresponds to the first nontrivial square blocked
  plaquette, i.e. a `2 x 2` square on the underlying fine lattice.
- `GLUEBALL_OPS` already shows repository precedent for:
  - zero-momentum projection
  - optional link smearing
  - DB-backed gluonic operator storage
- `GLUEBALL_OPS` does not directly emit the first-pass scalar timeslice
  operator contract needed here.
- `WILSLP` and `FUZZED_WILSON_LOOP` already measure Wilson loops inline, but
  they do not directly expose the sidecar-aware `0++` timeslice operator or
  the no-wrap reductions needed for the two-level estimator.
- In the split workflow, `TEMPORAL_ZONE_GAUGEBC` freezes forces and refreshed
  momenta on selected child-local temporal intervals but does not modify the
  stored gauge links themselves.
- There is not yet:
  - a `QACTDEN`-like measurement surface that emits raw `0++` timeslice
    operators at retained HMC updates
  - a child no-wrap reducer mode for scalar timeslice operators
  - a sidecar-aware cross-subdomain pair builder for `0++`
  - a two-level factorized estimator
  - a workflow spec for the level-0 / level-1 nesting in the scalar channel

## Why This Is The Right Next Step

The existing pieces already cover most of the mechanics we need:

- the parent lattice can be split into ordinary child configs
- the child configs can be evolved independently while holding the duplicated
  temporal boundary intervals fixed
- the scalar `0++` operator can already be built from spatial Wilson loops /
  plaquettes without any additional temporal support buffer

What is still missing is the measurement contract that says:

1. what counts as one level-0 sample
2. what counts as one level-1 child sample conditioned on that level-0 sample
3. how child-local reduced observables are combined into one factorized
   estimator
4. what direct comparison should validate the result before any production use
5. what thin measurement surface should expose the scalar timeslice operator
   without inventing a large new operator stack

That is the purpose of this spec.

## Level Structure

### Level 0: Parent Boundary Ensemble

For this spec, one level-0 sample is:

- one parent gauge configuration `U_k`
- selected from an ordinary whole-lattice gauge-only HMC history
- together with the split sidecar derived from `U_k`
- together with one matching retained `0++` measurement product for that saved
  parent configuration

The level-0 index is `k = 1, ..., N0`.

The sidecar derived from `U_k` determines the fixed boundary data seen by the
level-1 child runs:

- the duplicated frozen temporal boundary blocks
- the child-local frozen intervals
- the child-local to parent-global time maps

Important contract:

- the two child ensembles used in one factorized measurement must come from the
  same level-0 sample `U_k`
- mixing child measurements from different level-0 samples is invalid

### Outer-Sample Identifier Convention

The first version should use one simple stable identifier for each level-0
sample:

- `outer_sample_id = outer_<update_no>`

where:

- `<update_no>` is the decimal parent HMC update number of the retained saved
  parent configuration

Required first-pass convention:

- this `outer_sample_id` is the primary cross-file identifier for one level-0
  sample
- the numeric `update_no` should also be retained as machine-readable metadata
- filenames may use additional formatting, but the logical sample identifier
  should remain `outer_<update_no>`

### Level 1: Fixed-Boundary Child Ensembles

For one fixed level-0 sample `k` and one child `c in {0,1}`, a level-1 sample
is one child gauge configuration produced by gauge-only HMC on that child
lattice while:

- using the split child config from `U_k` as the starting state
- using the child-local `TEMPORAL_ZONE_GAUGEBC` intervals from the sidecar
- using the corresponding refreshed-momentum mask
- keeping the duplicated frozen boundary blocks fixed

The first version should organize level-1 sampling as one or more explicit
child HMC streams for each fixed level-0 sample.

Required first-pass convention:

- the workflow must allow `N1_streams >= 1`
- if `N1_streams > 1`, each stream must:
  - start from the same split child config derived from `U_k`
  - use a distinct RNG seed
- the first validation instance should use `N1_streams = 1`
- a later higher-statistics study may raise `N1_streams` above `1` without
  changing the estimator definition

If multiple retained `0++` measurements are collected from those streams, they
are all treated as level-1 observations conditioned on the same `U_k`.

Required first-pass stream-ordering convention:

- assign each child stream a deterministic `stream_id`
- order `stream_id` values monotonically within each `(outer_sample_id, child)`
  group

## Observable Building Blocks

### Selected Operator Family And Blocking Level

Reuse the existing scalar-glueball operator family already implemented in
`lib/meas/glue/gluecor.cc`.

For the first version, one run should select one explicit blocking level:

- `bl_level_selected`

Recommended initial default:

- `bl_level_selected = 1`

This keeps the first measurement surface inside the already-existing blocked
plaquette infrastructure while making the first operator an explicit `2 x 2`
square on the fine lattice.

Recommended first validation policy:

- use `bl_level_selected = 1` for the first plumbing and correctness pass
- interpret that first operator as the existing-infrastructure `2 x 2` square
  blocked plaquette
- postpone larger blocked / fuzzed operator levels until the `2 x 2` two-level
  workflow is validated end to end

### Local Spatial Wilson-Loop Building Block

For a chosen `bl_level_selected`, let:

- `P^{(bl)}_{mu,nu}(x,t)`

denote the real trace of the blocked purely spatial plaquette rooted at
spatial position `x`, timeslice `t`, and spatial plane `(mu,nu)`.

Only spatial planes orthogonal to the decay direction participate in the first
`0++` operator.

For the first implementation:

- `bl_level_selected = 1`
- the first operator is therefore a square `2 x 2` spatial plaquette on the
  fine lattice, built using the existing blocked-link glueball path

### Child Timeslice Operator

For level-0 sample `k`, child `c`, retained level-1 measurement `m`, and
child-local time `t`, define:

- `O_{k,c,m}(t) = sum_x sum_{(mu,nu) spatial planes} P^{(bl)}_{mu,nu}(x,t)`

where:

- the spatial sum is over that child timeslice
- the plane sum is over the three purely spatial plaquette orientations

This is the first-pass `0++` timeslice operator.

Repository-alignment note:

- this is the object most directly corresponding to `op0[t]` in
  `lib/meas/glue/gluecor.cc`

### Periodic Full-Lattice Reference Correlator

For a periodic full lattice with temporal extent `T`, define the raw
timeslice-sum correlator:

- `C_cfg^{periodic}(dt) = (1 / T) * sum_{t0=0}^{T-1} O_cfg(t0 + dt mod T) * O_cfg(t0)`

for `dt = 0, 1, ..., dt_max`, with `dt_max <= floor(T/2)`.

Important normalization note:

- the first stored object in this workflow should remain the correlator of raw
  timeslice sums
- the existing `gluecor.cc` output `glue0(dt)` uses a four-volume
  normalization
- therefore the reducer should record the spatial volume so an exact
  `gluecor`-normalized comparison remains recoverable offline

### Safe Child Timeslices

The usable child-local timeslices are determined by the split sidecar alone:

- begin from the child-local frozen intervals in the sidecar
- use `support_guard = 0` for the first scalar `0++` operator
- define `safe_slices_c` as the child-local timeslices not in a frozen
  interval

The first two-level estimator should only use timeslices in `safe_slices_c`.

This is the key operator-level difference from the first `0-+` workflow:

- no additional measurement buffer zone is needed next to a frozen interval
- the first scalar operator is local in Euclidean time

### Conditional Child Mean

For one fixed level-0 sample `k` and child `c`, define the conditional child
mean timeslice operator:

- `bar_O_{k,c}(t) = (1 / M_{k,c}) * sum_m O_{k,c,m}(t)`

where:

- `m` runs over all retained level-1 measurements for that `(k,c)` pair
- `M_{k,c}` is the retained measurement count for that `(k,c)` pair

Important interpretation note:

- `bar_O_{k,c}(t)` is a conditional expectation estimate at fixed boundary
  sample `U_k`
- only the outer average over `k` produces the full two-level estimator

### Secondary Child-Local Correlator Observable

In addition to the timeslice operator needed for the factorized cross-domain
estimator, the first implementation should preserve child-local no-wrap
correlator information as a sibling reporting product.

For one retained level-1 child measurement, define:

- `C_{k,c,m}^{local}(delta_t_child)`

using:

- child-local safe slices derived from the sidecar
- `support_guard = 0`
- no wrap-around

More explicitly, define:

- `source_slices_c(delta_t_child) = { t0 in safe_slices_c | t0 + delta_t_child is also in safe_slices_c }`

and then:

- `C_{k,c,m}^{local}(delta_t_child) = (1 / N_src_c(delta_t_child)) * sum_{t0 in source_slices_c(delta_t_child)} O_{k,c,m}(t0 + delta_t_child) * O_{k,c,m}(t0)`

where `N_src_c(delta_t_child) = |source_slices_c(delta_t_child)|`.

This child-local correlator is not part of the first factorized estimator
definition, but it should be emitted and summarized so the workflow can report:

- the subdomain-local correlator statistics
- the cross-subdomain factorized correlator statistics

side by side.

### Inner-Sample Weighting Convention

The first version should avoid ambiguous weighting across child streams.

Required first-pass convention:

- each child stream for a fixed `(k,c)` should use the same:
  - warmup discard
  - total updates
  - measurement stride
  - retained measurement count
- every retained `0++` measurement within a fixed `(k,c)` pair carries equal
  weight in `bar_O_{k,c}(t)`
- every retained post-warmup `0++` measurement should contribute to the
  retained child-local correlator sample as well

Therefore, if `N1_streams` streams are used and each retains `M_stream`
measurements, then:

- `M_{k,c} = N1_streams * M_stream`

and `bar_O_{k,c}(t)` is the ordinary pooled mean over those retained
measurements.

### Child-Summary Data Product

The first reducer stage should emit, for each fixed `(k,c)` pair:

- sidecar-derived metadata:
  - `outer_sample_id`
  - `child_id`
  - `child_nrow`
  - `child_local_to_global_t`
  - `frozen_local_intervals`
  - `support_guard`
  - `safe_slices`
- measurement metadata:
  - operator family identifier
  - `bl_level_selected`
  - spatial volume
  - any chosen spatial blocking / smearing controls
- reduction metadata:
  - `stream_count`
  - `stream_id` list
  - `retained_measurements_per_stream`
  - `retained_measurement_count`
  - retained `update_no` list
- per-measurement observables:
  - `O_{k,c,m}(t)` for every retained measurement and child-local timeslice
  - `C_{k,c,m}^{local}(delta_t_child)` for every retained measurement and
    every defined child-local no-wrap separation
- conditional mean observable:
  - `bar_O_{k,c}(t)` for every child-local timeslice
  - `bar_C_{k,c}^{local}(delta_t_child)` for every defined child-local
    no-wrap separation
- first-pass uncertainty summaries over the retained level-1 sample:
  - sample mean of `C_{k,c,m}^{local}(delta_t_child)`
  - sample standard deviation of `C_{k,c,m}^{local}(delta_t_child)`
  - sample standard error of `C_{k,c,m}^{local}(delta_t_child)`

This child-summary product is the natural input to the later two-level
cross-domain assembly, and it also preserves the child-local correlator
statistics as a sibling reporting product.

Required first-pass ordering convention:

- retained child measurements should be ordered by `(stream_id, update_no)`
- if two emitted products contain the same retained child measurements, they
  should use the same ordering

## First Estimator To Target

### Cross-Subdomain Localized Correlator

The first estimator should target source and sink operators that lie in the
interiors of opposite child domains.

This keeps the factorization clean:

- the source operator depends on active links in `child0`
- the sink operator depends on active links in `child1`
- conditioned on the fixed level-0 boundary sample, those two measurements can
  be averaged independently

This first version deliberately does not attempt to reconstruct the full
periodic correlator averaged over every parent timeslice.

### Global-Time Pair Construction

Let:

- `g_0(t)` be the child0 local-to-global time map from the sidecar
- `g_1(t)` be the child1 local-to-global time map from the sidecar
- `W_0 subset safe_slices_0`
- `W_1 subset safe_slices_1`

For the first version, choose:

- `W_0 = safe_slices_0`
- `W_1 = safe_slices_1`

For a positive parent-global timeslice separation `delta_t_parent`, define:

- `Pairs(delta_t_parent) = { (t0, t1) in W_0 x W_1 | g_1(t1) - g_0(t0) = delta_t_parent }`

with:

- no wrap-around in the first version
- only the forward ordering from `child0` to `child1`

Only separations with
`N_pairs(delta_t_parent) = |Pairs(delta_t_parent)| > 0` are defined.

### Per-Outer-Sample Two-Level Estimator

For one level-0 sample `k`, define:

- `C_k^{2lvl}(delta_t_parent) = (1 / N_pairs(delta_t_parent)) * sum_{(t0,t1) in Pairs(delta_t_parent)} bar_O_{k,1}(t1) * bar_O_{k,0}(t0)`

This is the first target quantity emitted by the two-level workflow.

### Outer Average

The ensemble estimator is then:

- `C^{2lvl}(delta_t_parent) = (1 / N0) * sum_k C_k^{2lvl}(delta_t_parent)`

The first implementation should report:

- the per-outer-sample `C_k^{2lvl}(delta_t_parent)`
- the outer-sample mean
- the outer-sample standard error
- the pair counts `N_pairs(delta_t_parent)`

### Pair Orientation And Normalization Convention

The first version should make the pairing and normalization fully explicit.

Required convention:

- only positive parent-global timeslice separations are used
- only forward cross-domain pairs from `child0` to `child1` are used
- each admissible pair in `Pairs(delta_t_parent)` carries equal weight
- no additional division by child temporal extent or spatial volume is applied
  inside `C_k^{2lvl}(delta_t_parent)`

Interpretation:

- the first stored object is a correlator of timeslice sums, not timeslice
  averages
- any later volume normalization can be applied offline because the geometry is
  recorded separately
- if exact comparison to the existing `gluecor` four-volume normalization is
  desired, the reducer should emit that as a derived reporting view rather than
  baking it into the first two-level estimator definition

The first version should also emit the explicit pair list for each
`delta_t_parent`, not only the pair count, so later auditing can verify the
exact parent-global timeslices that contributed.

Required first-pass pair-ordering convention:

- emit the explicit pair list sorted by
  `(delta_t_parent, child0_local_t, child1_local_t)`

## Direct Validation Observable

Before treating the two-level estimator as production-ready, validate it
against a direct parent-window measurement on the same geometry.

For the saved parent configuration `U_k`, define the parent timeslice operator:

- `O_k^{parent}(t) = O_{0++}^{parent}(t)`

using the same scalar operator family and the same `bl_level_selected` as the
child measurements.

First-pass measurement convention:

- the parent comparator should come from a thin post-HMC measurement pass on
  the retained saved parent configurations
- unlike the `0-+` path, the first version should not assume the existing HMC
  XML already contains the raw timeslice operator
- the first implementation should therefore use saved parent and child configs
  as the measurement surface for correctness bring-up
- that measurement pass must apply the same:
  - operator family
  - `bl_level_selected`
  - blocking / smearing metadata
  - normalization metadata

Then define the localized direct comparator:

- `C_k^{parent-window}(delta_t_parent) = (1 / N_pairs(delta_t_parent)) * sum_{(t0,t1) in Pairs(delta_t_parent)} O_k^{parent}(g_1(t1)) * O_k^{parent}(g_0(t0))`

This comparator uses:

- the same parent-global time pairs as the two-level estimator
- no factorization
- no periodic wrap

It is not the production measurement path. It is the first correctness check.

### Expected Relation Between The Two Estimators

The direct parent-window comparator and the factorized two-level estimator are
not expected to match configuration by configuration.

The intended relation is:

- `C_k^{parent-window}(delta_t_parent)` is one direct product measured on the
  saved parent field `U_k`
- `C_k^{2lvl}(delta_t_parent)` is an estimate of the same localized quantity
  after replacing each child factor by its level-1 conditional mean at fixed
  boundary data

Therefore the first validation compares:

- outer-ensemble means across `k`

and not:

- equality on each individual `k`

This distinction should remain explicit in any future checker or regression
test wording.

## First-Version Sampling Plan

### Outer Parent Chain

The level-0 parent workflow should be:

1. run ordinary whole-lattice gauge-only HMC on the parent lattice
2. save parent configurations at a configurable stride
3. retain one matching `0++` measurement product on the same saved update
   numbers
4. treat each saved parent configuration plus its matching `0++` measurement
   product as one level-0 sample
5. split each saved parent configuration into:
   - `child0`
   - `child1`
   - sidecar XML

Suggested configurable controls:

- `N0_warmup`
- `N0_samples`
- `N0_stride`
- `N0_measure_stride`
- parent gauge action and integrator parameters

Required first-pass convention:

- choose the retained scalar-measurement cadence so it matches the retained
  level-0 sample cadence
- in the simplest first version, use `N0_measure_stride = N0_stride`
- define one level-0 sample only at updates where:
  - a parent config is saved
  - a retained parent `0++` measurement product exists for that same update
    number
- preserve the same one-to-one contract through `outer_sample_id` and
  `update_no` when running the thin post-HMC measurement pass

Required first-pass control-matching convention:

- use identical values across the outer parent run, `child0`, and `child1` for
  every shared gauge-HMC control field that is meaningful on all three runs
- this includes, at minimum:
  - gauge action family and parameters
  - integrator type and integrator parameters
  - trajectory length / MD step controls
  - retained measurement cadence conventions
- only geometry- or role-specific fields may differ, such as:
  - `nrow`
  - gauge boundary-condition details
  - split sidecar and config I/O paths
  - RNG seeds

### Inner Child Sampling

For each level-0 sample `k` and each child `c`:

1. start from the split child config derived from `U_k`
2. launch `N1_streams` child HMC streams
   - if `N1_streams > 1`, start every stream from that same split child config
   - if `N1_streams > 1`, use distinct RNG seeds
3. optionally discard `N1_warmup` updates per stream
4. measure the selected `0++` operator every `N1_measure_stride` updates
5. retain every remaining post-warmup measurement for:
   - the conditional child timeslice mean
   - the child-local no-wrap correlator sample

Suggested configurable controls:

- `N1_streams`
- `N1_warmup`
- `N1_updates`
- `N1_measure_stride`
- `N1_retained_measurements_target`
- child gauge action and integrator parameters
- scalar-operator controls:
  - `bl_level_selected`
  - any chosen spatial blocking / smearing controls

Important first-pass convention:

- the spec remains generic in `N1_streams`
- the first validation instance should use `N1_streams = 1`
- later studies may increase `N1_streams` without changing the estimator
- define `N1_warmup` as a number of discarded child HMC updates at the start
  of each child stream
- do not interpret `N1_warmup` as a number of discarded retained measurements
- define the headline child-sampling target by retained measurements, not raw
  update count
- in particular, `N1_retained_measurements_target = 10` means ten retained
  post-warmup, post-thinning `0++` measurements for each fixed
  `(outer_sample_id, child_id)` group in the first concrete campaign target
- it is sufficient that every retained child measurement be identified
  unambiguously within the child stream data products
- the preferred identifier for one retained child measurement is:
  - `(outer_sample_id, child_id, stream_id, update_no)`

Because the first measurement surface is a thin post-HMC pass, the workflow
must also:

- save or otherwise retain the child configurations needed to support every
  retained measurement point
- keep the mapping from retained child config to
  `(outer_sample_id, child_id, stream_id, update_no)` explicit

## Recommended First Geometry

Use the first even-split lattice already targeted by the intended workflow:

- `parent_nrow = 8 8 8 8`
- `t_dir = 3`
- `cut0 = 0`
- `cut1 = 4`
- `frozen_width = 1`

Derived child geometry:

- `child0_nrow = 8 8 8 5`
- `child1_nrow = 8 8 8 5`
- frozen child-local intervals:
  - `[0,0]`
  - `[4,4]`
- `support_guard = 0`
- `safe_slices_0 = {1, 2, 3}`
- `safe_slices_1 = {1, 2, 3}`

Global safe timeslice maps:

- `child0 global safe slices = {1, 2, 3}`
- `child1 global safe slices = {5, 6, 7}`

For this geometry, the cross-subdomain pair counts are:

- `delta_t_parent = 2` -> `N_pairs = 1`
- `delta_t_parent = 3` -> `N_pairs = 2`
- `delta_t_parent = 4` -> `N_pairs = 3`
- `delta_t_parent = 5` -> `N_pairs = 2`
- `delta_t_parent = 6` -> `N_pairs = 1`

Recommended first comparison set:

- use all available cross-subdomain parent separations
- for the current first geometry, this evaluates to `{2, 3, 4, 5, 6}`

The edge separations `delta_t_parent = 2` and `delta_t_parent = 6` should
be kept in the first pass/fail criterion as ordinary measured points rather
than emitted only for visibility.

For the same geometry, the child-local no-wrap source counts are:

- `delta_t_child = 0` -> `N_src = 3`
- `delta_t_child = 1` -> `N_src = 2`
- `delta_t_child = 2` -> `N_src = 1`

Again, the edge case `delta_t_child = 2` should be emitted but should not
define the first pass/fail criterion.

## First Concrete Campaign Target

Once the basic plumbing is in place, the first concrete target campaign should
be organized around the same even-split `8^4` geometry with:

- `N0_samples = 40`
- `N1_streams = 1`
- `N1_retained_measurements_target = 10` per child per retained outer sample

Recommended provisional defaults for that first campaign target:

- `bl_level_selected = 1`
- `N1_measure_stride = 1`
- `N1_warmup = 0`

With those provisional defaults, this implies:

- `40` retained parent samples
- `2 * 40 * 10 = 800` retained child measurements across both child domains

If warmup discard or measurement thinning is later desired, the workflow should
keep the retained update lists explicit and should continue to interpret the
headline `40 x 10` description as:

- `40` retained outer samples
- `10` retained post-discard, post-thinning child measurements per child
  domain per retained outer sample

## Planned Repository Artifacts

The first implementation will likely need a small XML bundle under
`tests/gauge_subdomain_split/` plus a small full-lattice periodic check under
`tests/glueball_0pp/`.

Suggested checked-in artifacts:

- `tests/glueball_0pp/measure_full_lattice.glueball_0pp.template.ini.xml`
  - full-lattice `0++` measurement template at one selected blocking level
- `tests/glueball_0pp/glueball_0pp_corr.full_lattice.check.ini.xml`
  - periodic full-lattice reducer/checker template
- `tests/gauge_subdomain_split/hmc_parent.0pp_two_level_outer.template.ini.xml`
  - whole-lattice parent HMC template for level-0 samples
- `tests/gauge_subdomain_split/gauge_subdomain_split.0pp_two_level.template.ini.xml`
  - split input template for one saved parent sample
- `tests/gauge_subdomain_split/hmc_child.temporal_zone_glueball_0pp_2lvl.template.ini.xml`
  - level-1 child HMC template parameterized by `child_id`
- `tests/gauge_subdomain_split/measure_glueball_0pp_parent.template.ini.xml`
  - retained parent-config measurement template
- `tests/gauge_subdomain_split/measure_glueball_0pp_child.template.ini.xml`
  - retained child-config measurement template
- `tests/gauge_subdomain_split/glueball_0pp_parent_window.template.check.ini.xml`
  - reducer/checker template for the direct parent-window comparator
- `tests/gauge_subdomain_split/glueball_0pp_child.template.check.ini.xml`
  - reducer/checker template for one child summary
- `tests/gauge_subdomain_split/glueball_0pp_outer_sample.template.check.ini.xml`
  - reducer/checker template for one factorized outer sample
- one small generator script:
  - `tests/gauge_subdomain_split/generate_two_level_0pp_xml_bundle.sh`
  - expands those templates into concrete run XML for a chosen outer-sample set
  - writes the concrete XML under a checkout-local `cfgs/.../xml/` directory
  - also writes the outer-ensemble reducer input with the explicit list of
    generated outer-sample summaries

Generated run inputs and outputs should live under a dedicated checkout-local
hierarchy, for example:

- `cfgs/two_level_0pp/xml/...`
- `cfgs/two_level_0pp/outer_100/...`
- `cfgs/two_level_0pp/outer_200/...`

or, for the first concrete `8^4`, `40 x 10` campaign target:

- `cfgs/two_level_0pp_40x10/xml/...`
- `cfgs/two_level_0pp_40x10/outer_100/...`

This tree should hold:

- parent configs
- split sidecars
- child starting configs
- retained measurement inputs and outputs
- child HMC XML logs
- reducer summaries
- CSV outputs

## Proposed Checker Surfaces

The preferred first implementation path is still to keep the reduction logic in
a focused checker, but unlike the `0-+` path the raw measurement input should
come from a thin dedicated post-HMC `0++` measurement step rather than
directly from HMC XML.

### Stage 0: Parent Window Summary Input

Suggested first-pass input blocks:

- `Input/measurement_file`
  - one retained parent `0++` measurement product
- `Geometry/nrow`
  - parent lattice size
- `Geometry/decay_dir`
  - temporal direction
- `Reducer/mode = PARENT_WINDOW`
- `Reducer/sidecar_file`
- `Reducer/outer_sample_id`
- `Reducer/retained_update_no`
- `Reducer/bl_level_selected`
- `Reducer/compare_delta_t_parent_set`
- `Output/summary_file`
- `Output/csv_file`
- `Checks/expected_pairs`

Suggested parent-window outputs:

- copied input metadata
- `Pairs(delta_t_parent)` as explicit
  `(child0_local_t, child1_local_t, parent_t0, parent_t1)` tuples
- parent `timeslice_operator`
- optional periodic reference correlator in both:
  - raw timeslice-sum normalization
  - derived `gluecor` normalization
- `C_k^{parent-window}(delta_t_parent)`

### Stage 1: Child Summary Input

Suggested first-pass input blocks:

- `Input/measurement_files`
  - one `elem` per retained child `0++` measurement for a fixed
    `(outer_sample_id, child_id)` group
- `Geometry/nrow`
  - child lattice size
- `Geometry/decay_dir`
  - temporal direction
- `Reducer/mode = INTERIOR_NOWRAP_CHILD`
- `Reducer/support_guard = 0`
- `Reducer/sidecar_file`
- `Reducer/child_id`
- `Reducer/bl_level_selected`
- `Reducer/discard_updates`
  - or equivalent retained-update selection rule
- `Output/summary_file`
- `Output/csv_file`
- `Checks/expected_measurements`
- `Checks/expected_safe_slices`
- `Checks/expected_num_sources`

Suggested child-summary outputs:

- one XML summary file with:
  - copied input metadata
  - `safe_slices`
  - retained update numbers
  - per-measurement `timeslice_operator`
  - `conditional_mean_timeslice_operator`
  - per-measurement `child_local_correlator`
  - `conditional_mean_child_local_correlator`
  - `child_local_correlator_stddev`
  - `child_local_correlator_stderr`
- optional CSV with rows keyed by:
  - `update_no`
  - `series`
  - `index`
  - `value`

### Stage 2: One-Outer-Sample Assembly Input

The second stage should consume the reduced child summaries plus one parent
measurement summary for the same outer sample.

Suggested first-pass input blocks:

- `OuterSample/id`
- `OuterSample/sidecar_file`
- `Child0/summary_file`
- `Child1/summary_file`
- `Parent/summary_file`
- `Reducer/mode = TWO_LEVEL_CROSS_DOMAIN`
- `Reducer/compare_delta_t_parent_set`
- `Output/summary_file`
- `Output/csv_file`
- `Checks/expected_pairs`

Suggested one-outer-sample outputs:

- `Pairs(delta_t_parent)` as explicit
  `(child0_local_t, child1_local_t, parent_t0, parent_t1)` tuples
- `C_k^{2lvl}(delta_t_parent)`
- `C_k^{parent-window}(delta_t_parent)`
- `Delta_k(delta_t_parent)`
- copied child-local correlator summaries for:
  - `child0`
  - `child1`

### Stage 3: Outer-Ensemble Assembly Input

The final stage can remain a checker mode or a separate lightweight reducer,
but the contract should be explicit in the spec.

Suggested first-pass input blocks:

- `Input/outer_sample_summaries`
  - one `elem` per level-0 sample
- `Reducer/mode = TWO_LEVEL_OUTER_ENSEMBLE`
- `Reducer/compare_delta_t_parent_set`
- `Output/summary_file`
- `Output/csv_file`
- `Checks/min_outer_samples`

Suggested outer-ensemble outputs:

- `N0`
- `compare_delta_t_parent_set`
- outer-sample means of:
  - `C^{2lvl}(delta_t_parent)`
  - `C^{parent-window}(delta_t_parent)`
  - `Delta(delta_t_parent)`
- child-local no-wrap correlator reporting for each child:
  - mean of `bar_C_{k,c}^{local}(delta_t_child)` over outer samples
  - standard error of `bar_C_{k,c}^{local}(delta_t_child)` over outer samples
- outer-sample standard errors of the same quantities
- pass/fail status for the first `3 * SE` acceptance rule applied only to the
  parent-window comparison observable

Required first-pass outer-summary ordering convention:

- outer-sample summaries should be ordered by increasing retained parent
  `update_no`
- equivalently, by increasing `outer_sample_id = outer_<update_no>`

## Preferred First Measurement And Reducer Shape

Unlike the `0-+` path, the repository does not already have a turnkey inline
surface that exposes the raw scalar timeslice operator at retained HMC updates.

The preferred first implementation should therefore be:

1. add a thin measurement helper that reuses the existing scalar-glueball
   construction from `gluecor.cc`
2. expose one selected blocking level per run
3. emit the raw scalar timeslice operator `O_{0++}(t)` plus the metadata
   needed for exact offline interpretation
4. run that helper on retained saved parent and child configs for correctness
   bring-up
5. pair that helper with one focused reducer/checker under `mainprogs/tests/`

Suggested first-version checker:

- `mainprogs/tests/t_glueball_0pp_corr.cc`

This checker should support modes analogous to the `0-+` path:

1. `PERIODIC_ALL_T`
2. `INTERIOR_NOWRAP_CHILD`
3. `PARENT_WINDOW`
4. `TWO_LEVEL_CROSS_DOMAIN`
5. `TWO_LEVEL_OUTER_ENSEMBLE`

Important design preference:

- do not route the first version through the full `GLUEBALL_OPS` SDB stack
  unless there is a concrete later need for that richer operator surface
- do not generalize `WILSLP` or `FUZZED_WILSON_LOOP` first, because they do
  not directly expose the scalar timeslice operator contract needed by the
  factorized estimator
- do not add a dedicated inline HMC surface first, because the thin post-HMC
  pass on retained configs is sufficient for the initial correctness target

If a later inline HMC surface becomes desirable, it should emit exactly the
same operator contract so the reducer inputs do not need to change shape.

## Validation Ladder

Validation should happen in four phases.

### Phase 1: Existing Periodic Full-Lattice Scalar Reduction

Required checks:

1. the selected `0++` measurement surface emits the expected scalar timeslice
   operator `O_{0++}(t)`
2. the offline periodic correlator built from that operator reproduces the
   same chosen normalization convention every time
3. the derived `gluecor`-normalized periodic correlator matches the existing
   `Glueball_0pp/glue0` output from the reused scalar-glueball code at the
   same `bl_level_selected`

Required status before continuing:

- periodic scalar reduction works
- `O_{0++}(t)` construction is correct
- periodic correlator construction is correct

### Phase 2: Child No-Wrap Reduction

Required checks:

1. the reducer derives the expected child `safe_slices`
2. the reducer emits one `O_{0++}(t)` value per child-local timeslice
3. only the safe slices contribute to retained child-local summaries
4. the child pair counts implied by the sidecar and `support_guard = 0` match
   the expected geometry

### Phase 3: Single-Outer-Sample Two-Level Plumbing

Required checks:

1. one parent sample is split successfully
2. both child fixed-boundary runs complete successfully
3. frozen-link invariance still holds on both children
4. the reducer emits `bar_O_{k,0}(t)` and `bar_O_{k,1}(t)`
5. the reducer emits `C_k^{2lvl}(delta_t_parent)` for every defined
   `delta_t_parent`
6. the emitted `N_pairs(delta_t_parent)` values match the sidecar geometry

This phase is a plumbing check, not yet a statistical correctness check.

### Phase 4: Small-Ensemble Correctness Against Direct Parent Windows

This is the first real acceptance criterion.

The child-local no-wrap correlators remain part of the reported output in this
phase, but they do not define a separate first-pass pass/fail gate.

For each `delta_t_parent` in `compare_delta_t_parent_set`, compute over the
`N0` level-0 samples:

- `Delta_k(delta_t_parent) = C_k^{2lvl}(delta_t_parent) - C_k^{parent-window}(delta_t_parent)`
- `mu_Delta(delta_t_parent) = (1 / N0) * sum_k Delta_k(delta_t_parent)`

Use the ordinary sample standard deviation:

- `s_Delta(delta_t_parent)^2 = (1 / (N0 - 1)) * sum_k (Delta_k(delta_t_parent) - mu_Delta(delta_t_parent))^2`

and the standard error:

- `SE_Delta(delta_t_parent) = s_Delta(delta_t_parent) / sqrt(N0)`

First-pass acceptance criterion:

- `|mu_Delta(delta_t_parent)| <= 3 * SE_Delta(delta_t_parent)` for every
  `delta_t_parent` in `compare_delta_t_parent_set`

If `SE_Delta(delta_t_parent) = 0`, interpret the check as:

- pass if `mu_Delta(delta_t_parent) = 0`
- fail otherwise

This paired outer-sample comparison is preferred because both estimators are
built from the same saved parent samples.

No additional first-pass requirement is imposed that `child0` and `child1`
child-local correlator summaries agree with each other. Those summaries are
reported for inspection, not used as an independent acceptance criterion in
the first validation workflow.

The intended first concrete target for this phase is:

- `parent_nrow = 8 8 8 8`
- `N0 = 40`
- `N1_retained_measurements_target = 10`

## Important Invariants

The workflow must preserve the following invariants:

1. The child-local frozen intervals come only from the split sidecar.
2. The level-1 child runs for one outer sample must all use those same frozen
   intervals.
3. A factorized child pair is valid only if both children come from the same
   outer sample.
4. The first version uses `support_guard = 0`, so it excludes only frozen
   slices and not their immediate temporal neighbors.
5. The first two-level estimator uses only safe child-local timeslices.
6. The first version uses no wrap-around in the cross-subdomain pair builder.
7. Validation compares the new estimator against a direct parent-window
   observable before any production claims are made.

## Non-Goals For The First Version

- no claim of improved efficiency on day one
- no full translationally averaged parent correlator reconstruction
- no same-child contribution assembly into a larger estimator
- no stitched-parent measurement path
- no production `0++` mass extraction
- no Wilson-flowed or cooled operator
- no variational basis over multiple blocking levels
- no automated campaign orchestration
- no more-than-two-domain generalization

## Follow-Up Path

Once this first workflow is specified and validated, the natural follow-up path
is:

1. decide whether the first estimator should stay cross-subdomain-localized or
   be extended into a fuller correlator assembly
2. decide whether the first measurement path should remain a thin post-HMC
   measurement step or later be replaced by a dedicated inline surface
3. evaluate whether one selected blocking level is sufficient or whether a
   small blocked/fuzzed operator basis should be added
4. add flowed, cooled, or otherwise improved scalar operators
5. measure larger ensembles and perform proper uncertainty analysis
6. only then use the two-level workflow as a real spectroscopy tool

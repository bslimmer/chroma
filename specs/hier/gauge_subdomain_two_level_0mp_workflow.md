# Gauge-Subdomain Two-Level 0-+ Workflow Spec

## Goal

Define a first two-level measurement workflow that combines:

1. parent-lattice gauge evolution
2. split child-lattice gauge-only HMC with frozen duplicated temporal boundary
   blocks
3. subdomain-localized `QACTDEN`-based `0-+` timeslice operators

into one nested estimator for a cross-subdomain glueball correlator.

The immediate target is not a production spectroscopy campaign. The immediate
target is:

1. generate a sequence of level-0 parent boundary samples
2. split each parent sample into two child lattices plus sidecar metadata
3. generate independent level-1 child ensembles at fixed boundary data
4. reduce each child ensemble to conditional means of the localized
   pseudoscalar timeslice operator
5. form a factorized cross-subdomain `0-+` correlator from those conditional
   means
6. validate that estimator against a direct parent-window measurement on a
   small lattice

This is a workflow-and-estimator spec. It is not yet a performance-tuning
document and it is not yet a mass-fit campaign spec.

## Scope

This first version intentionally stays narrow:

- gauge-only HMC only
- two temporal child domains only
- `QACTDEN` / `naiveTopCharge` only
- zero momentum only
- offline reduction only
- one two-level estimator built from operators in opposite child interiors
- one small-lattice validation path against a direct parent-window observable
- current split sidecar geometry contract only

This first version does not require:

- fermion monomials
- SMD support
- Wilson flow or any smoothing in the first pass
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
- child-local topology validation:
  `specs/hier/gauge_subdomain_topology_validation.md`
- first `QACTDEN`-based `0-+` operator and correlator reduction:
  `specs/hier/glueball_0mp_correlator.md`

The conceptual pattern is the standard two-level one used in the multilevel
lattice literature: for each level-0 boundary sample, average subdomain-local
observables independently at level 1, then average the product of those
level-1 means over the level-0 ensemble.

## Current Repo Facts

- `gauge_subdomain_split` already writes:
  - `child0` config
  - `child1` config
  - sidecar XML with:
    - parent and child lattice sizes
    - child-local to parent-global time maps
    - child-local frozen intervals
- `TEMPORAL_ZONE_GAUGEBC` plus the HMC gauge-monomial momentum-mask
  autodiscovery already define the fixed-boundary child evolution semantics.
- `QACTDEN` already emits the local field needed for the first `0-+` operator:
  - `naiveTopCharge`
- `mainprogs/tests/t_qactden_0mp_corr.cc` already parses `QACTDEN` output and
  reduces the periodic full-lattice correlator, but it currently implements
  only `PERIODIC_ALL_T`.
- There is not yet:
  - a child no-wrap reducer mode
  - a sidecar-aware cross-subdomain pair builder
  - a two-level factorized estimator
  - a workflow spec for the level-0 / level-1 nesting

## Why This Is The Right Next Step

The existing pieces already cover the mechanics we need:

- the parent lattice can be split into ordinary child configs
- the child configs can be evolved independently while holding the duplicated
  temporal boundary intervals fixed
- the child-local `0-+` operator can be obtained from inline `QACTDEN`

What is still missing is the measurement contract that says:

1. what counts as one level-0 sample
2. what counts as one level-1 child sample conditioned on that level-0 sample
3. how child-local reduced observables are combined into one factorized
   estimator
4. what direct comparison should validate the result before any production use

That is the purpose of this spec.

## Level Structure

### Level 0: Parent Boundary Ensemble

For this spec, one level-0 sample is:

- one parent gauge configuration `U_k`
- selected from an ordinary whole-lattice gauge-only HMC history
- together with the split sidecar derived from `U_k`

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

If multiple retained `QACTDEN` measurements are collected from those streams,
they are all treated as level-1 observations conditioned on the same `U_k`.

Required first-pass stream-ordering convention:

- assign each child stream a deterministic `stream_id`
- order `stream_id` values monotonically within each `(outer_sample_id, child)`
  group

## Observable Building Blocks

### Local Density

Reuse the first-version local `0-+` operator from
`specs/hier/glueball_0mp_correlator.md`:

- `q(x) = naiveTopCharge(x)`

### Child Timeslice Operator

For level-0 sample `k`, child `c`, retained level-1 measurement `m`, and
child-local time `t`, define:

- `O_{k,c,m}(t) = sum_x q_{k,c,m}(x,t)`

where the spatial sum is over that child timeslice.

### Safe Child Timeslices

The usable child-local timeslices are determined exactly as in
`specs/hier/glueball_0mp_correlator.md`:

- begin from the child-local frozen intervals in the sidecar
- use `support_guard = 1` for the first `QACTDEN` operator
- define `safe_slices_c` as the child-local timeslices not in a frozen
  interval and not within `support_guard` of one

The first two-level estimator should only use timeslices in `safe_slices_c`.

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
- because `q(x)` is pseudoscalar, `bar_O_{k,c}(t)` need not vanish for one
  fixed `U_k`
- only the outer average over `k` enforces the full-ensemble symmetry

### Secondary Child-Local Correlator Observable

In addition to the timeslice operator needed for the factorized cross-domain
estimator, the first implementation should preserve the child-local no-wrap
correlator information from the earlier `0-+` workflow.

For one retained level-1 child measurement, define:

- `C_{k,c,m}^{local}(delta_t_child)`

using the same no-wrap child-interior construction from
`specs/hier/glueball_0mp_correlator.md`, with:

- child-local safe slices derived from the sidecar
- `support_guard = 1`
- no wrap-around

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
- every retained `QACTDEN` measurement within a fixed `(k,c)` pair carries
  equal weight in `bar_O_{k,c}(t)`
- every retained post-warmup `QACTDEN` measurement should contribute
  to the retained child-local correlator sample as well

Therefore, if `N1_streams` streams are used and each retains `M_stream`
measurements, then:

- `M_{k,c} = N1_streams * M_stream`

and `bar_O_{k,c}(t)` is the ordinary pooled mean over those retained
measurements.

This deliberately avoids, in the first version, any extra policy question
about whether longer streams should carry more weight than shorter ones.

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

- `O_k^{parent}(t) = sum_x q_k^{parent}(x,t)`

using the same `QACTDEN`-based field on the full parent lattice.

First-pass measurement convention:

- the parent comparator should come from inline `QACTDEN` emitted during the
  outer parent HMC run
- it should not require a separate one-shot parent measurement pass in the
  first version
- the retained parent `QACTDEN` measurements used for the comparator must
  correspond to the same retained update numbers as the saved level-0 parent
  configurations

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

1. run ordinary whole-lattice gauge-only HMC on the parent lattice with inline
   `QACTDEN`
2. save parent configurations at a configurable stride
3. retain parent `QACTDEN` measurements on the same saved update numbers
4. treat each saved parent configuration plus its matching parent `QACTDEN`
   record as one level-0 sample
5. split each saved parent configuration into:
   - `child0`
   - `child1`
   - sidecar XML

Suggested configurable controls:

- `N0_warmup`
- `N0_samples`
- `N0_stride`
- `N0_qactden_stride`
- parent gauge action and integrator parameters

Required first-pass convention:

- choose the retained parent `QACTDEN` cadence so it matches the retained
  level-0 sample cadence
- in the simplest first version, use `N0_qactden_stride = N0_stride`
- define one level-0 sample only at updates where:
  - a parent config is saved
  - a retained parent `QACTDEN` record exists for that same update number
- do not allow extra retained parent `QACTDEN` measurements without matching
  saved parent configs in the first version

Required first-pass control-matching convention:

- use identical values across the outer parent run, `child0`, and `child1` for
  every shared gauge-HMC control field that is meaningful on all three runs
- this includes, at minimum:
  - gauge action family and parameters
  - integrator type and integrator parameters
  - trajectory length / MD step controls
  - retained `QACTDEN` cadence conventions
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
4. measure `QACTDEN` every `N1_measure_stride` updates
5. retain every remaining post-warmup `QACTDEN` measurement for:
   - the conditional child timeslice mean
   - the child-local no-wrap correlator sample

Suggested configurable controls:

- `N1_streams`
- `N1_warmup`
- `N1_updates`
- `N1_measure_stride`
- child gauge action and integrator parameters

Important first-pass convention:

- the spec remains generic in `N1_streams`
- the first validation instance should use `N1_streams = 1`
- later studies may increase `N1_streams` without changing the estimator
- define `N1_warmup` as a number of discarded child HMC updates at the start
  of each child stream
- do not interpret `N1_warmup` as a number of discarded retained measurements
- unlike the parent level-0 sampling rule, no strict requirement is imposed
  that child save cadence match child `QACTDEN` measurement cadence
- it is sufficient that every retained child measurement be identified
  unambiguously within the child stream data products
- the preferred identifier for one retained child measurement is:
  - `(outer_sample_id, child_id, stream_id, update_no)`

## Recommended First Geometry

Reuse the larger child-correlator geometry already identified in
`specs/hier/glueball_0mp_correlator.md`:

- `parent_nrow = 4 4 4 12`
- `t_dir = 3`
- `cut0 = 0`
- `cut1 = 6`
- `frozen_width = 1`

Derived child geometry:

- `child0_nrow = 4 4 4 7`
- `child1_nrow = 4 4 4 7`
- frozen child-local intervals:
  - `[0,0]`
  - `[6,6]`
- `support_guard = 1`
- `safe_slices_0 = {2, 3, 4}`
- `safe_slices_1 = {2, 3, 4}`

Global safe timeslice maps:

- `child0 global safe slices = {2, 3, 4}`
- `child1 global safe slices = {8, 9, 10}`

For this geometry, the cross-subdomain pair counts are:

- `delta_t_parent = 4` -> `N_pairs = 1`
- `delta_t_parent = 5` -> `N_pairs = 2`
- `delta_t_parent = 6` -> `N_pairs = 3`
- `delta_t_parent = 7` -> `N_pairs = 2`
- `delta_t_parent = 8` -> `N_pairs = 1`

Recommended first comparison set:

- derive `compare_delta_t_parent_set` from the rule
  `{ delta_t_parent | N_pairs(delta_t_parent) >= 2 }`
- for the current first geometry, this evaluates to `{5, 6, 7}`

The edge separations `delta_t_parent = 4` and `delta_t_parent = 8` should
still be emitted for visibility, but they should not define the first
pass/fail criterion.

## Planned Repository Artifacts

The first implementation will likely need a small XML bundle under
`tests/gauge_subdomain_split/`.

Suggested checked-in artifacts:

- `hmc_parent.0mp_two_level_outer.template.ini.xml`
  - whole-lattice parent HMC template for level-0 samples
  - includes inline `QACTDEN` at the retained outer-sample cadence
- `gauge_subdomain_split.0mp_two_level.template.ini.xml`
  - split input template for one saved parent sample
- `hmc_child.temporal_zone_qactden_0mp_2lvl.template.ini.xml`
  - level-1 child HMC template parameterized by `child_id`
- `qactden_0mp_parent_window.template.check.ini.xml`
  - reducer/checker template for the direct parent-window comparator
- `qactden_0mp_child.template.check.ini.xml`
  - reducer/checker template for one child summary
- `qactden_0mp_outer_sample.template.check.ini.xml`
  - reducer/checker template for one factorized outer sample
- one small generator script
  - expands those templates into concrete run XML for a chosen outer-sample set
  - writes the concrete XML under a checkout-local `cfgs/.../xml/` directory
  - also writes the outer-ensemble reducer input with the explicit list of
    generated outer-sample summaries

Generated run inputs and outputs should live under a dedicated checkout-local
hierarchy, for example:

- `cfgs/two_level_0mp/xml/...`
- `cfgs/two_level_0mp/outer_100/...`
- `cfgs/two_level_0mp/outer_200/...`

This tree should hold:

- parent configs
- split sidecars
- child starting configs
- child HMC XML logs
- reducer summaries
- CSV outputs

## Proposed Checker Surfaces

The preferred first implementation path is still to extend
`mainprogs/tests/t_qactden_0mp_corr.cc`, but the XML surfaces should be
planned now so the later code does not grow ad hoc.

### Stage 0: Parent Window Summary Input

The outer parent HMC run should be reduced directly from its HMC XML log
rather than through a separate parent measurement job.

Suggested first-pass input blocks:

- `Input/file`
  - outer parent HMC XML log
- `Geometry/nrow`
  - parent lattice size
- `Geometry/decay_dir`
  - temporal direction
- `Reducer/mode = PARENT_WINDOW`
- `Reducer/sidecar_file`
- `Reducer/outer_sample_id`
- `Reducer/retained_update_no`
- `Reducer/compare_delta_t_parent_set`
- `Output/summary_file`
- `Output/csv_file`
- `Checks/expected_pairs`

Suggested parent-window outputs:

- copied input metadata
- `Pairs(delta_t_parent)` as explicit `(child0_local_t, child1_local_t, parent_t0, parent_t1)` tuples
- parent `timeslice_operator`
- `C_k^{parent-window}(delta_t_parent)`

### Stage 1: Child Summary Input

The existing checker already accepts one HMC XML plus geometry/reducer
parameters. The child-summary mode should extend that shape rather than replace
it.

Suggested first-pass input blocks:

- `Input/file`
  - child HMC XML log
- `Geometry/nrow`
  - child lattice size
- `Geometry/decay_dir`
  - temporal direction
- `Reducer/mode = INTERIOR_NOWRAP_CHILD`
- `Reducer/support_guard`
- `Reducer/sidecar_file`
- `Reducer/child_id`
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
  - direct parent-window source reduced from the outer parent HMC XML
- `Reducer/mode = TWO_LEVEL_CROSS_DOMAIN`
- `Reducer/compare_delta_t_parent_set`
- `Output/summary_file`
- `Output/csv_file`
- `Checks/expected_pairs`

Suggested one-outer-sample outputs:

- `Pairs(delta_t_parent)` as explicit `(child0_local_t, child1_local_t, parent_t0, parent_t1)` tuples
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

## Preferred First Reducer Shape

The preferred first implementation is to extend
`mainprogs/tests/t_qactden_0mp_corr.cc` rather than invent a separate
measurement stack immediately.

The natural staged extension is:

1. add `INTERIOR_NOWRAP` child reduction mode
2. add sidecar-aware safe-slice and pair construction
3. add a `TWO_LEVEL_CROSS_DOMAIN` mode that consumes paired child summaries for
   one level-0 sample
4. add a direct parent-window mode for validation
5. add an outer-ensemble mode that consumes one summary per level-0 sample

If that checker becomes unwieldy, a later refactor into shared library code or
a second focused executable is acceptable, but it should not be the first
design move.

## Validation Ladder

Validation should happen in four phases.

### Phase 1: Existing Periodic Full-Lattice Reduction

This phase is already defined in `specs/hier/glueball_0mp_correlator.md`.

Required status before continuing:

- periodic `QACTDEN` reduction works
- `O_q(t)` construction is correct
- periodic correlator construction is correct

### Phase 2: Child No-Wrap Reduction

Required checks:

1. the reducer derives the expected child `safe_slices`
2. the reducer emits one `O_q(t)` value per child-local timeslice
3. only the safe slices contribute to retained child-local summaries
4. the child pair counts implied by the sidecar and `support_guard = 1` match
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

## Important Invariants

The workflow must preserve the following invariants:

1. The child-local frozen intervals come only from the split sidecar.
2. The level-1 child runs for one outer sample must all use those same frozen
   intervals.
3. A factorized child pair is valid only if both children come from the same
   outer sample.
4. The first two-level estimator uses only safe child-local timeslices.
5. The first version uses no wrap-around in the cross-subdomain pair builder.
6. Validation compares the new estimator against a direct parent-window
   observable before any production claims are made.

## Non-Goals For The First Version

- no claim of improved efficiency on day one
- no full translationally averaged parent correlator reconstruction
- no same-child contribution assembly into a larger estimator
- no stitched-parent measurement path
- no production `0-+` mass extraction
- no Wilson-flowed operator
- no automated campaign orchestration
- no more-than-two-domain generalization

## Follow-Up Path

Once this first workflow is specified and validated, the natural follow-up path
is:

1. decide whether the first estimator should stay cross-subdomain-localized or
   be extended into a fuller correlator assembly
2. evaluate whether `QACTDEN` output volume requires a lighter-weight transport
   format
3. add flowed or otherwise smoothed `0-+` operators
4. measure larger ensembles and perform proper uncertainty analysis
5. only then use the two-level workflow as a real spectroscopy tool

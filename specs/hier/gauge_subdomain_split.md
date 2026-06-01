# Gauge Subdomain Split Spec

## Goal

Add a gauge-only split/stitch utility that decomposes one parent lattice into
two child lattices along the temporal direction, duplicates two shared frozen
time-slice boundary regions into both children, and later stitches the two
children back into a single parent gauge field.

Chosen public name:

- Feature/doc name: `gauge_subdomain_split`
- C++ types:
  - `GaugeSubdomainSplitParams`
  - `GaugeSubdomainSplitChild`
  - `GaugeSubdomainSplitResult`
- C++ APIs:
  - `splitGaugeSubdomains(...)`
  - `stitchGaugeSubdomains(...)`
- Proposed files:
  - `lib/util/gauge/gauge_subdomain_split.h`
  - `lib/util/gauge/gauge_subdomain_split.cc`
- Focused test:
  - `mainprogs/tests/t_gauge_subdomain_split.cc`

This first version is an in-memory C++ utility only. It is not an XML-facing
feature and it does not itself run HMC.

This first version has no halo feature. The shared frozen boundary is exactly
the duplicated time-slice blocks.

Child lattices are non-periodic in their local temporal direction.

## Current Repo Facts

- Gauge fields in the target scope use:
  ```c++
  multi1d<LatticeColorMatrix>
  ```
- There is currently no split/stitch utility in the repository for gauge
  fields.
- Time-slice subset support already exists in:
  - `lib/util/ft/time_slice_set.h`
  - `lib/util/ft/time_slice_set.cc`
- Gauge-configuration truncation support already exists in:
  - `lib/io/writeszin.cc` via `writeSzinTrunc(...)`
  That helper writes one truncated time window, but it does not construct two
  child lattices, duplicate shared frozen boundaries, or stitch child lattices
  back together.
- The current `TEMPORAL_ZONE_GAUGEBC` spec in
  `specs/hier/temporal_zone_gaugebc.md` freezes only the selected time-slice
  intervals. That is the intended frozen-boundary semantics for this split/stitch utility too.

## API Surface

Use a small parameter type plus result metadata:

```c++
struct GaugeSubdomainSplitParams {
  int t_dir;
  int cut0;
  int cut1;
  int frozen_width;
};

struct GaugeSubdomainSplitInterval {
  int t_start;
  int t_end;
};

struct GaugeSubdomainSplitChild {
  multi1d<LatticeColorMatrix> u;
  multi1d<int> nrow;
  multi1d<int> local_to_global_t;
  multi1d<GaugeSubdomainSplitInterval> frozen_local_intervals;
};

struct GaugeSubdomainSplitResult {
  GaugeSubdomainSplitParams param;
  multi1d<int> parent_nrow;
  GaugeSubdomainSplitChild child0;
  GaugeSubdomainSplitChild child1;
};
```

Suggested entry points:

```c++
GaugeSubdomainSplitResult
splitGaugeSubdomains(const multi1d<LatticeColorMatrix>& parent_u,
                     const GaugeSubdomainSplitParams& param);

multi1d<LatticeColorMatrix>
stitchGaugeSubdomains(const GaugeSubdomainSplitResult& split);
```

The first version may let `stitchGaugeSubdomains(...)` consume the current
state of `split.child0.u` and `split.child1.u`, so downstream code can update
those child gauge fields in place before stitching.

## Parameter Semantics

- `t_dir` is required and identifies which lattice axis is the temporal
  direction for this lattice.
- This utility may split only along `t_dir`. It is not a generic spatial-domain
  splitter.
- `cut0` and `cut1` are the global parent-lattice time coordinates where the
  two frozen shared boundary blocks begin.
- `frozen_width` is the shared frozen width in time slices and must satisfy
  `frozen_width > 0`.

Interpretation:

- Let `Lt = Layout::lattSize()[t_dir]`.
- Normalize the two cuts to increasing global time order:
  - `cut_left = min(cut0, cut1)`
  - `cut_right = max(cut0, cut1)`
- Define a forward-wrap helper:
  ```c++
  advance(t, n) = (t + n) % Lt
  ```
- Define the two frozen boundary blocks as ordered lists of global times:
  ```c++
  F0[k] = advance(cut_left,  k), 0 <= k < frozen_width
  F1[k] = advance(cut_right, k), 0 <= k < frozen_width
  ```

Wrapping frozen blocks are allowed. For example, if `Lt = 16`, `cut_right = 14`,
and `frozen_width = 3`, then `F1 = [14, 15, 0]`.

## Validity Rules

The split must reject invalid parameter sets with a clear error.

Required checks:

- `0 <= t_dir < Nd`
- `0 <= cut0 < Lt`
- `0 <= cut1 < Lt`
- `cut0 != cut1`
- `frozen_width > 0`
- The two frozen boundary blocks must not overlap.
- The two frozen boundary blocks must not touch.
- The chosen cuts and width must leave two nonempty independent subdomains
  between the shared frozen blocks.

Equivalent operational rule:

- If `A` is the forward interval after `F0` and before `F1`, and `B` is the
  forward interval after `F1` and before `F0`, then both `A` and `B` must have
  strictly positive length.

That is, the first version requires two nonempty active temporal subdomains plus
two shared frozen boundary blocks.

## Split Construction

Let:

```c++
lenA = (cut_right - (cut_left + frozen_width) + Lt) % Lt;
lenB = (cut_left - (cut_right + frozen_width) + Lt) % Lt;
```

After validation, `lenA > 0` and `lenB > 0`.

Define the two active parent-lattice temporal subdomains as ordered lists:

```c++
A[k] = advance(cut_left  + frozen_width, k), 0 <= k < lenA
B[k] = advance(cut_right + frozen_width, k), 0 <= k < lenB
```

Then build the child-local temporal orderings:

- `child0` local time order: `F0 + A + F1`
- `child1` local time order: `F1 + B + F0`

In other words:

- `child0` owns active subdomain `A` and shares both frozen boundary blocks.
- `child1` owns active subdomain `B` and shares both frozen boundary blocks.

The child time coordinates are renumbered locally from `0`.

Child extents:

- All non-temporal extents match the parent lattice exactly.
- Only `nrow[t_dir]` changes:
  - `child0.nrow[t_dir] = 2 * frozen_width + lenA`
  - `child1.nrow[t_dir] = 2 * frozen_width + lenB`

Required metadata:

- `child0.local_to_global_t` stores the global parent time coordinate for each
  child-local time slice in `child0`.
- `child1.local_to_global_t` stores the same for `child1`.
- `frozen_local_intervals` must contain exactly two inclusive intervals for each
  child:
  - first frozen block: `[0, frozen_width - 1]`
  - second frozen block:
    `[child_nrow_t - frozen_width, child_nrow_t - 1]`

This metadata is intended to be sufficient for later stitch logic and for
downstream temporal boundary-condition wiring on the child lattices.

## Gauge-Field Copy Semantics

For each child-local site:

1. Map its child-local time coordinate through `local_to_global_t`.
2. Keep all non-temporal coordinates unchanged.
3. Copy all `Nd` rooted gauge links from the matching parent site into the
   child site.

That includes the `mu = t_dir` links rooted on the local last time slice, even
though their forward endpoint lies outside the child lattice when the child is
treated as non-periodic.

The split utility should copy those link values unchanged. The downstream
non-periodic gauge-evolution code is then responsible for not constructing
plaquettes or other paths that wrap across the local temporal boundary.

No extra data is duplicated beyond the child-local time slices listed above.
There is no halo slice and no separately duplicated outgoing-link region in the
first version.

## Frozen Boundary Semantics

The shared frozen boundary blocks are exactly the duplicated time-slice regions
`F0` and `F1`.

Consequences:

- All link directions rooted on those frozen time slices belong to the shared
  frozen boundary.
- There is no additional incoming-link or outgoing-link special case in this
  first version.
- There is no halo option in this first version.

This matches the updated `TEMPORAL_ZONE_GAUGEBC` spec, which freezes only the
selected time-slice intervals themselves.

Downstream child-lattice HMC evolution is expected to wire those two child-local
frozen intervals into `TEMPORAL_ZONE_GAUGEBC` or equivalent future non-periodic
boundary logic. That downstream evolution work is out of scope for this spec.

## Stitch Semantics

`stitchGaugeSubdomains(...)` reconstructs a parent gauge field from the two
child gauge fields plus the split metadata.

Use the same normalized cut ordering and `local_to_global_t` metadata that were
used to construct the children.

Ownership rules:

- Active subdomain `A` comes from `child0`.
- Active subdomain `B` comes from `child1`.
- Shared frozen blocks `F0` and `F1` must match exactly between `child0` and
  `child1`.

Stitch procedure:

1. Validate that the child metadata is consistent with the stored parent shape
   and split parameters.
2. Allocate a parent-sized gauge field.
3. Copy all rooted links on active global times in `A` from `child0`.
4. Copy all rooted links on active global times in `B` from `child1`.
5. For each rooted link on global times in `F0` and `F1`:
   - compare the `child0` and `child1` values exactly
   - abort with a clear error if any duplicated frozen-boundary link differs
   - otherwise copy one checked-equal version into the parent

The first version requires exact equality, not a tolerance-based comparison.
Since the frozen boundary links are duplicated in memory and are not supposed to
be altered, any mismatch should be treated as a bug or a downstream contract
violation.

Because there is no halo in this first version:

- the first active slice after a frozen block belongs only to its active child
  subdomain
- rooted links on that active slice are stitched from that owning child only
- only rooted links on the frozen slices themselves are exact-checked as shared
  data

## Non-Periodic Child-Lattice Contract

The child lattices are non-periodic in their local temporal direction.

This spec intentionally stops at the split/stitch utility and its metadata
contract. It does not specify the downstream gauge-action, integrator, or HMC
changes required to evolve those child lattices without computing plaquettes
across the local temporal edges.

That downstream non-periodic child-evolution support is required for the full
workflow, but it is a separate feature from this split/stitch utility.

## No XML Or Factory Wiring

This first version is not a factory-registered or XML-facing feature.

So the first implementation should not add:

- a new factory name
- a `GaugeBC` registration hook
- a new `tests/*.xml` smoke input just for this utility

If a later workflow wants file-backed or XML-driven split/stitch operations,
that should be specified separately.

## Build System Updates

When the utility is implemented, add the new utility files and focused test to
both build systems:

- `lib/Makefile.am`
  - add `util/gauge/gauge_subdomain_split.h`
  - add `util/gauge/gauge_subdomain_split.cc`
- `lib/CMakeLists.txt`
  - add the same header and source
- `mainprogs/tests/Makefile.am`
  - add `t_gauge_subdomain_split`
- `mainprogs/tests/CMakeLists.txt`
  - add the same test executable

## Tests

Add at least one focused code test executable for the split/stitch utility.

Suggested focused tests:

1. Non-wrapping split:
   - use a small parent lattice such as `Lt = 8`
   - choose `cut0 = 1`, `cut1 = 5`, `frozen_width = 1`
   - verify:
     - `child0.local_to_global_t = [1, 2, 3, 4, 5]`
     - `child1.local_to_global_t = [5, 6, 7, 0, 1]`
     - each child reports frozen local intervals `[0,0]` and `[4,4]`
2. Wrapping frozen block:
   - use a case such as `Lt = 16`, `cut0 = 4`, `cut1 = 14`, `frozen_width = 3`
   - verify the wrapped block is represented as `[14, 15, 0]`
   - verify child-local time orderings and extents
3. Exact round trip:
   - split a nontrivial parent gauge field
   - stitch immediately without modifying either child
   - verify the reconstructed parent matches exactly
4. Interior-only evolution surrogate:
   - split a parent gauge field
   - change only `child0` active subdomain `A`
   - change only `child1` active subdomain `B`
   - leave both duplicated frozen blocks untouched
   - stitch and verify the reconstructed parent contains:
     - updated `A` from `child0`
     - updated `B` from `child1`
     - exact frozen blocks from the shared boundary
5. Frozen-block mismatch detection:
   - change any rooted link on a duplicated frozen slice in only one child
   - verify `stitchGaugeSubdomains(...)` aborts with a clear error
6. Invalid parameter rejection:
   - invalid `t_dir`
   - invalid cuts
   - `cut0 == cut1`
   - `frozen_width <= 0`
   - overlapping frozen blocks
   - touching frozen blocks
   - any case leaving `A` or `B` empty

No XML smoke input is required in the first version because this utility is a
pure C++ helper, not an XML-facing feature.

## Validation Checklist

1. Build `chroma-hier`.
2. Run `t_gauge_subdomain_split` successfully.
3. Verify `splitGaugeSubdomains(...)` constructs the expected child time
   orderings in both a non-wrapping and a wrapping case.
4. Verify the child-local frozen intervals are exactly the first and last
   `frozen_width` time slices.
5. Verify split followed immediately by stitch reproduces the original parent
   gauge field exactly.
6. Verify interior-only modifications stitch back into the correct parent
   subdomains.
7. Verify any mismatch in duplicated frozen-boundary links aborts with a clear
   error.
8. Verify invalid parameter sets are rejected with clear errors.

## Scope Notes

- Gauge only. Fermion-field split/stitch support is out of scope.
- No halo feature in the first version.
- No incoming-link or outgoing-link special cases beyond the links already
  rooted on the duplicated frozen time slices.
- No XML or factory interface in the first version.
- No file-based split/stitch workflow in the first version.
- No downstream non-periodic HMC or gauge-action integration in the first
  version.
- Exact frozen-boundary equality is required at stitch time in the first
  version; tolerance-based comparison is out of scope.

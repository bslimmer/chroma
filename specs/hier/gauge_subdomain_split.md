# Gauge Subdomain Split Spec

## Goal

Define a gauge-only split/stitch workflow that decomposes one parent lattice
into two child lattices along the temporal direction, duplicates two shared
frozen boundary blocks into both children, evolves the children in separate
runs as ordinary gauge configurations, and later stitches the two evolved
children back into one parent gauge field.

This spec intentionally avoids introducing wrapper datatypes for parent and
child gauge fields. Parent and child gauge configurations remain ordinary:

```c++
multi1d<LatticeColorMatrix>
```

The split/stitch contract is carried by a small geometry plan plus a persisted
sidecar XML file.

## Design Principles

- Parent and child gauge fields remain ordinary gauge configurations.
- Split geometry, field extraction, and stitching are separate concerns.
- The child configurations should be usable by ordinary child-lattice HMC runs.
- The shared frozen boundary is represented explicitly and identically in both
  children.
- Split metadata lives outside the child gauge configs in a separate sidecar
  XML file.

## Current Repo Facts

- Gauge fields in scope are:
  ```c++
  multi1d<LatticeColorMatrix>
  ```
- The repository already has ordinary gauge-config read/write entry points:
  - `lib/util/gauge/gauge_startup.cc` via `gaugeStartup(...)`
  - `lib/io/gauge_io.cc` via `writeGauge(...)`
- The repository already has time-slice subset support in:
  - `lib/util/ft/time_slice_set.h`
  - `lib/util/ft/time_slice_set.cc`
- The repository already has truncated gauge-write support in:
  - `lib/io/writeszin.cc` via `writeSzinTrunc(...)`
- Existing gauge-config I/O understands full configurations for the active
  layout, but does not know the split maps, duplicated frozen boundary blocks,
  or stitch ownership rules defined here.
- `TEMPORAL_ZONE_GAUGEBC` already provides the intended child-run freezing
  semantics: freeze the selected time slices themselves and leave stored gauge
  links unchanged.

## Scope

This first version specifies:

1. parameter validation
2. split geometry construction
3. child-local to parent-global time maps
4. child-local frozen intervals
5. split-aware gauge-field extraction
6. exact-match stitch rules on duplicated frozen slices
7. persisted sidecar metadata
8. user-facing split/stitch tools

This spec does not define the downstream gauge-action or integrator changes
required to run non-periodic child-lattice HMC. It defines the geometry and
data contract that those runs will use.

## Public Types

Use a small parameter type:

```c++
struct GaugeSubdomainSplitParams {
  int t_dir;
  int cut0;
  int cut1;
  int frozen_width;
};
```

Use a small interval type for frozen local time windows:

```c++
struct GaugeSubdomainSplitInterval {
  int t_start;
  int t_end;
};
```

Use one split-plan type containing geometry and mapping metadata only:

```c++
struct GaugeSubdomainSplitPlan {
  GaugeSubdomainSplitParams param;
  multi1d<int> parent_nrow;
  multi1d<int> child0_nrow;
  multi1d<int> child1_nrow;
  multi1d<int> child0_local_to_global_t;
  multi1d<int> child1_local_to_global_t;
  multi1d<GaugeSubdomainSplitInterval> child0_frozen_local_intervals;
  multi1d<GaugeSubdomainSplitInterval> child1_frozen_local_intervals;
};
```

Important non-goals:

- no `GaugeSubdomainSplitChild`
- no `GaugeSubdomainSplitResult`
- no new gauge-field container type for parent or child configurations

## Library API

The library layer should expose a split-plan builder plus split/stitch helpers.

Important implementation note:

- QDP lattice objects are bound to the active `Layout`.
- Parent and child lattices therefore cannot all remain live as
  `multi1d<LatticeColorMatrix>` objects across different lattice sizes in one
  process at the same time.
- The shared split/stitch core should therefore operate on layout-independent
  host-side link buffers and provide small helpers to snapshot from and
  materialize back into ordinary gauge fields.

Required plan builder:

```c++
GaugeSubdomainSplitPlan
makeGaugeSubdomainSplitPlan(const multi1d<int>& parent_nrow,
                            const GaugeSubdomainSplitParams& param);
```

Required snapshot helper:

```c++
void
snapshotGaugeField(const multi1d<LatticeColorMatrix>& u,
                   multi1d<ColorMatrix>& site_links);
```

Required materialization helper:

```c++
void
materializeGaugeField(const multi1d<int>& nrow,
                      const multi1d<ColorMatrix>& site_links,
                      multi1d<LatticeColorMatrix>& u);
```

Required extraction helper:

```c++
void
extractGaugeSubdomain(const multi1d<ColorMatrix>& parent_links,
                      const GaugeSubdomainSplitPlan& plan,
                      int child_id,
                      multi1d<ColorMatrix>& child_links);
```

Contract:

- `child_id` is `0` or `1`
- `parent_links` contains one `ColorMatrix` per rooted link on the parent
  lattice, stored in a deterministic site-major ordering
- `child_links` receives the extracted child-lattice rooted links in the same
  ordering convention
- the helper fills `child_links` by copying from `parent_links` according to
  the plan

Required stitch helper:

```c++
void
stitchGaugeSubdomains(const GaugeSubdomainSplitPlan& plan,
                      const multi1d<ColorMatrix>& child0_links,
                      const multi1d<ColorMatrix>& child1_links,
                      multi1d<ColorMatrix>& parent_links);
```

Contract:

- `child0_links` and `child1_links` contain child-lattice rooted links in the
  same deterministic ordering used by `snapshotGaugeField(...)`
- the helper validates the duplicated frozen slices exactly
- the helper reconstructs the parent rooted-link buffer in `parent_links`

Ordinary gauge-field workflows use these helpers in sequence:

1. create the active parent or child layout
2. snapshot an ordinary `multi1d<LatticeColorMatrix>` gauge field into a host
   buffer
3. recreate `Layout` for the next lattice size as needed
4. extract or stitch in the host-buffer representation
5. materialize the result back into an ordinary
   `multi1d<LatticeColorMatrix>` gauge field for I/O or evolution

The primary user workflow is through the user-facing tools below. These library
helpers exist so the tools and focused tests can share one split/stitch core.

## User-Facing Tools

Provide two user-facing executables:

- `gauge_subdomain_split`
- `gauge_subdomain_stitch`

Suggested source locations:

- `mainprogs/main/gauge_subdomain_split.cc`
- `mainprogs/main/gauge_subdomain_stitch.cc`

These tools should sit on top of the library split-plan/extract/stitch logic
and operate on ordinary gauge configs plus one explicit sidecar XML file.

### `gauge_subdomain_split`

Responsibilities:

1. read one parent config using the ordinary gauge startup path
2. build the split plan
3. materialize `child0` and `child1` as ordinary child-layout gauge configs
4. write the two child configs as QIO configs
5. write one sidecar XML file containing the persisted split metadata

### `gauge_subdomain_stitch`

Responsibilities:

1. read the sidecar XML file
2. read the evolved `child0` and `child1` configs
3. reconstruct the split plan from the sidecar metadata
4. validate exact agreement on duplicated frozen slices
5. stitch the child configs back into one parent config
6. write the reconstructed parent config as a QIO config

### Tool I/O Rules

- Input configs should use the existing `Cfg_t` / `gaugeStartup(...)`
  conventions.
- Input parent and child configs may use any format already supported by
  `gaugeStartup(...)`.
- Output configs are ordinary QIO gauge configs written through
  `writeGauge(...)`.
- Split metadata must not be embedded into the child gauge configs.
- The sidecar XML is the only persisted split metadata artifact required by
  this spec.

## Parameter Semantics

- `t_dir` identifies the temporal direction for the lattice.
- This feature splits only along `t_dir`.
- `cut0` and `cut1` are the global parent-lattice time coordinates where the
  two frozen shared boundary blocks begin.
- `frozen_width` is the shared frozen width in time slices and must satisfy
  `frozen_width > 0`.

Interpretation:

```c++
Lt = parent_nrow[t_dir]
cut_left  = min(cut0, cut1)
cut_right = max(cut0, cut1)
advance(t, n) = (t + n) % Lt
```

Frozen boundary blocks:

```c++
F0[k] = advance(cut_left,  k), 0 <= k < frozen_width
F1[k] = advance(cut_right, k), 0 <= k < frozen_width
```

Wrapping frozen blocks are allowed.

## Validity Rules

Reject invalid parameter sets with a clear error.

Required checks:

- `0 <= t_dir < Nd`
- `0 <= cut0 < Lt`
- `0 <= cut1 < Lt`
- `cut0 != cut1`
- `frozen_width > 0`
- frozen blocks `F0` and `F1` must not overlap
- frozen blocks `F0` and `F1` must not touch
- both active subdomains between the frozen blocks must be nonempty

Operational rule:

- let `A` be the forward interval after `F0` and before `F1`
- let `B` be the forward interval after `F1` and before `F0`
- both `A` and `B` must have strictly positive length

## Split Geometry

Define:

```c++
lenA = (cut_right - (cut_left + frozen_width) + Lt) % Lt;
lenB = (cut_left - (cut_right + frozen_width) + Lt) % Lt;
```

Then:

```c++
A[k] = advance(cut_left  + frozen_width, k), 0 <= k < lenA
B[k] = advance(cut_right + frozen_width, k), 0 <= k < lenB
```

Child-local time orderings:

- `child0`: `F0 + A + F1`
- `child1`: `F1 + B + F0`

Child extents:

- all non-temporal extents match the parent
- only `nrow[t_dir]` changes
- `child0_nrow[t_dir] = 2 * frozen_width + lenA`
- `child1_nrow[t_dir] = 2 * frozen_width + lenB`

Required metadata:

- `child0_local_to_global_t`
- `child1_local_to_global_t`
- `child0_frozen_local_intervals`
- `child1_frozen_local_intervals`

Frozen local intervals:

- each child has exactly two intervals
- first interval: `[0, frozen_width - 1]`
- second interval:
  `[child_nrow[t_dir] - frozen_width, child_nrow[t_dir] - 1]`

This split plan is the canonical geometry description used by the extraction
helpers, stitch helper, user-facing tools, and sidecar XML.

## Gauge-Field Copy Semantics

When a child gauge field is materialized, each child-local time slice maps back
to one parent-global time slice through the relevant `local_to_global_t` array.

Copy rule:

1. map child-local time to parent-global time
2. keep non-temporal coordinates unchanged
3. copy all rooted links `mu = 0 ... Nd - 1` from the matching parent site

This includes the `mu = t_dir` links rooted on the local last child time slice,
even though the forward endpoint lies outside the child lattice once the child
is treated as non-periodic.

The split utility should copy those link values unchanged. Downstream
child-lattice evolution code is then responsible for not constructing paths
that wrap across the local temporal boundary.

## Frozen Boundary Semantics

The shared frozen boundary is exactly the duplicated time-slice blocks `F0` and
`F1`.

Consequences:

- all link directions rooted on those frozen slices are shared/frozen
- there is no halo slice in this version
- there are no extra incoming-link or outgoing-link regions
- each child exposes exactly the two local frozen intervals stored in the plan

Those child-local frozen intervals are the values that downstream child-lattice
HMC should feed into `TEMPORAL_ZONE_GAUGEBC`.

## Stitch Semantics

Stitching reconstructs a parent gauge field using:

- the split plan
- a `child0` gauge field on `child0_nrow`
- a `child1` gauge field on `child1_nrow`

Ownership rules:

- active subdomain `A` comes from `child0`
- active subdomain `B` comes from `child1`
- duplicated frozen blocks `F0` and `F1` must match exactly between children

Required stitch behavior:

1. validate the supplied child geometries against the stored plan
2. allocate or reconstruct a parent-sized gauge field
3. copy active region `A` from `child0`
4. copy active region `B` from `child1`
5. compare duplicated frozen slices exactly
6. abort with a clear error on any frozen-slice mismatch

Exact equality is required in the first version. Since the frozen boundary data
is duplicated and should remain unchanged during child evolution, any mismatch
is a contract violation.

## Separate-Run Workflow

The intended workflow is:

1. read a parent config on the parent layout
2. build the split plan
3. materialize and write `child0` and `child1` as ordinary QIO gauge configs
4. run ordinary child-lattice HMC in separate runs
5. read the evolved child configs later
6. stitch them back into a reconstructed parent config

This workflow assumes sequential layout recreation or separate executables/runs,
not simultaneous live parent/child lattice objects.

## Persisted Metadata

The split/stitch workflow must preserve:

- the split parameters
- the parent lattice shape
- the child lattice shapes
- the child-local to parent-global time maps, or enough equivalent information
  to reconstruct them exactly
- the child-local frozen intervals needed to configure
  `TEMPORAL_ZONE_GAUGEBC` for each child run

The storage location is fixed by this spec:

- use one explicit sidecar XML file
- do not embed the split metadata into the child gauge configs themselves
- do not rely on remembered filenames or in-memory state across runs

The sidecar is path-agnostic. It describes the split geometry and frozen
intervals, not a registry of specific child-config filenames.

## Sidecar XML Contract

The sidecar XML is the canonical persisted description of the split.

Suggested top-level shape:

```xml
<GaugeSubdomainSplitInfo>
  <Param>
    <t_dir>3</t_dir>
    <cut0>1</cut0>
    <cut1>5</cut1>
    <frozen_width>1</frozen_width>
  </Param>
  <Parent>
    <nrow>4 4 4 8</nrow>
  </Parent>
  <Child0>
    <nrow>4 4 4 5</nrow>
    <local_to_global_t>1 2 3 4 5</local_to_global_t>
    <frozen_local_intervals>
      <elem><t_start>0</t_start><t_end>0</t_end></elem>
      <elem><t_start>4</t_start><t_end>4</t_end></elem>
    </frozen_local_intervals>
  </Child0>
  <Child1>
    <nrow>4 4 4 5</nrow>
    <local_to_global_t>5 6 7 0 1</local_to_global_t>
    <frozen_local_intervals>
      <elem><t_start>0</t_start><t_end>0</t_end></elem>
      <elem><t_start>4</t_start><t_end>4</t_end></elem>
    </frozen_local_intervals>
  </Child1>
</GaugeSubdomainSplitInfo>
```

Required properties of the sidecar:

- it is sufficient to reconstruct the split plan exactly in a later stitch run
- it is sufficient to configure `TEMPORAL_ZONE_GAUGEBC` on each child lattice
- it does not change the fact that the child gauge configs themselves remain
  ordinary configs
- it may be passed independently of child-config filenames so evolved child
  configs can be stitched later without rewriting the sidecar

## Tool XML Shape

The split tool should consume a Chroma-style XML input containing:

- one parent `Cfg`
- one `GaugeSubdomainSplitParams`
- two QIO child output targets
- one sidecar XML output path

The stitch tool should consume a Chroma-style XML input containing:

- one sidecar XML input path
- one `child0` input `Cfg`
- one `child1` input `Cfg`
- one QIO parent output target

Because the child configs are meant to behave like ordinary configs, child HMC
inputs should continue to reference the child configs directly and should pull
their frozen-interval information from the sidecar when constructing the
`TEMPORAL_ZONE_GAUGEBC` block.

## Build System Updates

When this spec is implemented, add the new library utility and user-facing
tools to both build systems.

Expected library files:

- `lib/util/gauge/gauge_subdomain_split.h`
- `lib/util/gauge/gauge_subdomain_split.cc`

Expected executable sources:

- `mainprogs/main/gauge_subdomain_split.cc`
- `mainprogs/main/gauge_subdomain_stitch.cc`

Build list updates:

- `lib/Makefile.am`
- `lib/CMakeLists.txt`
- `mainprogs/main/Makefile.am`
- `mainprogs/main/CMakeLists.txt`

## Tests

Add at least one focused code test executable for the core split-plan and
stitch logic.

Suggested focused tests:

1. Non-wrapping split geometry:
   - use a small parent lattice such as `Lt = 8`
   - choose `cut0 = 1`, `cut1 = 5`, `frozen_width = 1`
   - verify:
     - `child0_local_to_global_t = [1, 2, 3, 4, 5]`
     - `child1_local_to_global_t = [5, 6, 7, 0, 1]`
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
   - stitch and verify the reconstructed parent contains updated `A`, updated
     `B`, and exact frozen blocks
5. Frozen-block mismatch detection:
   - change any rooted link on a duplicated frozen slice in only one child
   - verify `stitchGaugeSubdomains(...)` throws or reports a clear error
6. Invalid parameter rejection:
   - invalid `t_dir`
   - invalid cuts
   - `cut0 == cut1`
   - `frozen_width <= 0`
   - overlapping frozen blocks
   - touching frozen blocks
   - any case leaving `A` or `B` empty

Suggested end-to-end tool tests:

1. Split a parent config into two QIO child configs plus one sidecar XML file.
2. Reload the sidecar XML and both child configs.
3. Stitch back into a parent QIO config.
4. Verify exact frozen-boundary agreement checks.

## Validation Checklist

1. Build the new library utility and user-facing tools successfully.
2. Run the focused split/stitch test executable successfully.
3. Verify both a non-wrapping and wrapping case produce the expected child time
   orderings.
4. Verify the child-local frozen intervals are exactly the first and last
   `frozen_width` child time slices.
5. Verify split followed immediately by stitch reproduces the original parent
   gauge field exactly.
6. Verify interior-only child modifications stitch back into the correct parent
   subdomains.
7. Verify any duplicated frozen-boundary mismatch aborts with a clear error.
8. Verify the sidecar is sufficient to configure `TEMPORAL_ZONE_GAUGEBC` on the
   child lattices.

## Scope Notes

- Gauge only. Fermion-field split/stitch support is out of scope.
- The first user-facing write path is QIO-only.
- Output-format parity with the broader set of supported input config types is
  out of scope for the first version.
- No halo feature in the first version.
- No extra incoming-link or outgoing-link special cases beyond links already
  rooted on the duplicated frozen slices.

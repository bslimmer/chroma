# Gauge Subdomain HMC Spec

## Status

Implementation draft for the first version. The core workflow, naming, restart
contract, and child-edge semantics below are intended as the working contract
for implementation.

## Goal

Define the next feature layer above:

- `TEMPORAL_ZONE_GAUGEBC`
- `splitGaugeSubdomains(...)`
- `stitchGaugeSubdomains(...)`

so that two child lattices can be evolved independently with quenched HMC and
stitched back into one parent gauge field afterward.

This draft is specifically about the child-evolution workflow. The split/stitch
utility itself is already specified separately in
`specs/hier/gauge_subdomain_split.md`.

The `subdomain_hmc` workflow is intended to reuse that existing split/stitch
implementation directly. It does not define a second split algorithm or a
second stitch algorithm in XML. The XML surface only decides when to call the
existing utility and how to drive child evolution around it.

## Constraints Already Agreed

- Minimize the amount of new code.
- Do not alter the original Chroma HMC implementation in place.
- Prefer an additive feature built on top of existing HMC outputs and existing
  split/stitch functionality.
- First target scope is gauge-only child evolution.
- A wrapper class, wrapper driver, or equivalent additive layer is acceptable if
  it is genuinely needed.
- The first user-facing version should be XML-facing and runnable through a
  normal input workflow.
- The XML-facing entry point should be a dedicated new main program that mirrors
  the structure of `mainprogs/main/hmc.cc`, rather than a new mode inside the
  existing `hmc` executable.
- The dedicated executable name should be `subdomain_hmc`.
- The top-level XML block should be `SubdomainHMC`.
- The feature should be flexible enough to accommodate improved gauge actions,
  not only the simplest Wilson-plaquette case.
- The design must not rely on "the duplicated frozen boundaries happen to make
  periodic evolution acceptable". Child non-periodicity must be explicit.
- A child-specific fixed-boundary helper is acceptable, especially if it reuses
  existing Schrödinger-style masking and the temporal-zone interval semantics
  already introduced in this workstream.
- In the first version, that child-specific fixed-boundary helper may be wired
  only through the new child-HMC driver rather than exposed immediately as a
  general standalone XML/factory BC type.
- The workflow must support multiple child trajectories before stitching back to
  the parent lattice.
- The workflow must support child-lattice measurements on the independently
  evolved child trajectory chains before stitching.
- After the split, `child0` and `child1` should be treated as independent
  chains. They may run different numbers of trajectories and different
  measurement schedules before stitch.
- Stitching should be optional. A run may stop after advancing and saving the
  independent child chains without stitching them back to a parent in that same
  invocation.
- The XML surface should provide separate measurement blocks for `child0` and
  `child1`, plus an optional post-stitch parent measurement block.
- The workflow should support saving and restoring the child chains separately so
  the independent chain state is preserved across restart.
- The required per-child restart payload should include at least:
  - the child gauge field
  - the child momenta
  - the child RNG state
  - the child update counter
- Restart artifacts should use one manifest plus separate per-child payload
  files.
- When a stitched parent output is requested, the workflow should keep both the
  stitched parent output and the child restart artifacts by default.
- The XML surface should use one shared HMC definition for both children, with
  optional per-child overrides rather than two fully separate full HMC blocks.
- The child-specific fixed-boundary helper name should be
  `SUBDOMAIN_FIXED_GAUGEBC`.
- Any final workflow must preserve the duplicated frozen boundary slices so
  stitch-time exact equality remains valid.

## Current Repo Facts

- `TEMPORAL_ZONE_GAUGEBC` currently zeroes force-like fields through
  `GaugeBC::zero(P&)`.
- `TEMPORAL_ZONE_GAUGEBC` intentionally leaves `modify(Q&)` as a no-op.
- The current split/stitch utility duplicates two frozen temporal blocks into
  both children and records child-local frozen intervals.
- The current split/stitch spec explicitly stops before downstream child HMC.
- Existing HMC refreshes momenta internally and updates all gauge links on the
  active `Layout`.
- Existing gauge actions and link updates are written against the ordinary
  periodic lattice behavior of the active `Layout`.
- The repository already contains Schrödinger-functional gauge BC support and
  leapfrog tests that combine those BCs with multiple gauge-action families,
  including Wilson-like, anisotropic, and rectangle-based actions.

## Existing Repo Support That Looks Relevant

There is already meaningful support in the repository for fixed gauge
boundaries:

- Schrödinger-functional gauge BC classes can:
  - modify gauge links seen by gauge actions
  - zero gauge-like force fields on masked boundary links
- The existing gauge-state path applies `GaugeBC::modify(...)` to the copied
  gauge field stored in the state object.
- The leapfrog test suite already exercises gauge monomials with
  `SCHROEDINGER_NONPERT_GAUGEBC` for several gauge-action families.

This means the repository is not starting from zero on fixed-boundary gauge
evolution.

## Why Existing Repo Support Is Not Yet Sufficient

The current repository support is still not the full child-sublattice workflow
needed here.

### 1. Existing HMC Still Refreshes And Updates Everything

The current HMC trajectory code refreshes momenta internally for every link and
the current MD leap updates every link on the active lattice.

So even if a gauge action later zeroes boundary forces, the current trajectory
path does not by itself guarantee that child boundary links remain exactly
frozen.

### 2. Existing Gauge BC Support Is Boundary-Condition Logic, Not A Generic Non-Periodic-Lattice Mode

Current gauge actions and link-path construction still operate on the ordinary
periodic `Layout`.

So the repository does not appear to offer a generic switch saying "treat this
whole lattice as non-periodic in the temporal direction". Instead, current
non-periodic behavior is encoded through BC-specific masking and
action-specific logic.

### 3. Existing Schrödinger BCs Are Not A Direct Match For Split Children

The existing Schrödinger gauge BC classes are designed for outer decay
boundaries with specific boundary-field construction rules.

They are not yet a direct expression of:

- "freeze these child-local boundary slices"
- "use the already copied child boundary links as the fixed values"
- "support the exact path exclusions required by split/stitch child lattices"

So existing Schrödinger support is strong precedent, but not yet a drop-in
solution for this subdomain workflow.

## Design Requirements

Any accepted design will need to address both of these:

1. Frozen child-boundary slices must not acquire nonzero molecular-dynamics
   momentum.
2. The child evolution must explicitly forbid the unwanted local temporal
   wrap-around paths near the child edges.

More concretely, this must honor the split/stitch contract that:

- the child gauge field still contains the copied `mu = t_dir` links rooted on
  the local last time slice
- those links may be frozen and carried through stitch
- but no plaquette, rectangle, or longer improved-action loop may be
  constructed if it would step forward off the child lattice and wrap back
  through the periodic `Layout`

## Current Best-Fit Direction

Given the decisions recorded above, the current best-fit direction is:

- keep stock HMC, Hamiltonian, integrator, and gauge-action code unchanged
- add a child-specific fixed-boundary helper that:
  - reinstalls the duplicated child boundary links as fixed links in
    `modify(Q&)`
  - zeroes gauge-force fields on the corresponding masked links in `zero(P&)`
  - reuses the Schrödinger-style `loop_extent` idea so improved actions can
    declare how far temporal loops reach near the boundary
- add an additive child-HMC workflow layer around the existing machinery
- add a dedicated XML-facing main program that performs:
  - split
  - repeated and independent `child0` and `child1` trajectory evolution
  - separate child-local measurement schedules for `child0` and `child1`
  - optional stitch on explicit workflow completion
  - optional parent measurements after stitch
  - separate child-chain save/restart through one manifest plus child payloads
  - retention of child restart artifacts even when a stitched parent is also
    written

This is the chosen first-version direction for the implementation draft.

However, the fixed-boundary helper and momentum masking are not by themselves
enough to satisfy the split/stitch edge-link requirement. Ordinary periodic
gauge-action code still builds loops with `shift(...)` across the active
`Layout`. So the final design must also include child-specific path exclusion
semantics for any loop that would step off a child temporal edge.

## Chosen Public Names

The first version should use these public names:

- executable: `subdomain_hmc`
- top-level XML block: `SubdomainHMC`
- child fixed-boundary helper: `SUBDOMAIN_FIXED_GAUGEBC`

## Is Momentum Zeroing The Only Additional Functionality Needed?

Not quite.

If the child-specific fixed-boundary helper described above exists, then
momentum masking is the main missing trajectory-level feature. But the complete
workflow still needs a few distinct responsibilities:

1. Child fixed-boundary helper:
   - installs the exact duplicated child boundary links seen by the gauge state
   - zeroes forces on the masked boundary region
   - carries the boundary-width or loop-extent rule needed by improved actions
2. Momentum suppression:
   - ensures refreshed MD momenta vanish on the child frozen intervals so the
     link update does not move those links
3. Non-periodic path exclusion:
   - ensures copied edge-rooted `mu = t_dir` child links are not used to form
     plaquettes or longer loops that would wrap across the local temporal edge
   - scales from plaquettes to improved-action loops via an explicit temporal
     reach or `loop_extent` rule
4. XML-facing orchestration:
   - exposes a runnable workflow that splits a parent field
   - supports repeated child-HMC updates before stitch
   - supports separate `child0` and `child1` measurement blocks before stitch
   - supports separate save/restart of the independent child chains
   - writes restart artifacts as one manifest plus separate per-child payloads
   - stitches the result only when the workflow requests it
   - can run optional parent measurements after stitch
   - keeps both stitched-parent and child-chain outputs by default when both
     are requested

So if we count only the new HMC-trajectory behavior, then "force the momentum to
be zero in those regions" is one main extra piece. If we count the full
feature, it is one key piece among a small set of additive components, and it
is not sufficient without explicit path exclusion at the child temporal edges.

## Independent Child Chains Requirement

The child lattices are not just temporary one-trajectory work buffers.

The intended workflow is to create independently evolved child trajectory chains
that may:

- run for many HMC trajectories before any stitch happens
- perform measurements directly on the child lattices while those chains are
  being generated
- use different trajectory counts or different measurement schedules on
  `child0` and `child1`
- stitch back to a parent lattice only after the desired child evolution window
  has completed

The workflow should preserve the child-chain state directly, rather than
requiring an intermediate stitched-parent checkpoint as the only restart path.
At minimum, that means the design should allow separate persistence of the child
gauge fields and the core child-chain state needed to resume the independent
updates faithfully.

Required per-child restart contents for the first version:

- child gauge field
- child momenta
- child RNG state
- child update counter

This means the driver design must account for more than a one-shot
split-update-stitch helper. It needs chain-oriented orchestration semantics.

## Chosen Architectural Direction

The first implementation should follow the additive-wrapper-plus-explicit
non-periodic-support direction.

Chosen shape:

- Leave the stock HMC, Hamiltonian, integrator, and existing gauge-action
  implementations untouched in place.
- Add a dedicated XML-facing main program `subdomain_hmc` that mirrors the
  high-level structure of `mainprogs/main/hmc.cc`.
- Add a child-HMC orchestration layer that owns:
  - parent start vs child-restart startup
  - split and stitch invocation
  - per-child `Layout` switching
  - per-child trajectory loops
  - per-child measurement dispatch
  - per-child checkpoint/restart output
- Add a driver-owned child fixed-boundary helper named
  `SUBDOMAIN_FIXED_GAUGEBC`.
- Add additive momentum-masking behavior for child frozen intervals.
- Add child-specific non-periodic path-exclusion support so plaquettes and
  improved-action loops never wrap across the child-local temporal edge.

This is the required first-version architecture because a wrapper-only approach
without explicit path exclusion would still allow ordinary periodic loop
construction at the copied child edge.

## Chosen HMC XML Contract Shape

The XML surface should use one shared HMC definition for both children, with
optional per-child overrides.

Meaning:

- one common XML block defines the default monomials, Hamiltonian, and
  integrator
- both children inherit that common definition
- a child may optionally replace whole top-level HMC sub-blocks if it needs to
  diverge

The first version should support replacement overrides at the block level, not
fine-grained XML patch/merge semantics. In practice, a child override may
replace:

- `Monomials`
- `Hamiltonian`
- `MDIntegrator`

The child `nrow` is never provided in XML for the HMC definition. It is derived
from the split metadata or restart manifest.

## Concrete XML Surface

The top-level input should be:

```xml
<SubdomainHMC>
  <Start>...</Start>
  <SharedHMC>...</SharedHMC>
  <Child0>...</Child0>
  <Child1>...</Child1>
  <Outputs>...</Outputs>
</SubdomainHMC>
```

### `Start`

This block selects one of two startup modes.

Mode A: start from a parent gauge configuration

```xml
<Start>
  <Mode>FROM_PARENT_CFG</Mode>
  <Cfg>
    ...
  </Cfg>
  <Split>
    <t_dir>...</t_dir>
    <cut0>...</cut0>
    <cut1>...</cut1>
    <frozen_width>...</frozen_width>
  </Split>
</Start>
```

Mode B: restart from saved child-chain state

```xml
<Start>
  <Mode>FROM_CHILD_RESTART</Mode>
  <ChildRestartManifest>...</ChildRestartManifest>
</Start>
```

Required semantics:

- `FROM_PARENT_CFG` reads one parent gauge field and immediately calls
  `splitGaugeSubdomains(...)`.
- `FROM_CHILD_RESTART` reconstructs the child chains directly from the saved
  manifest and child payload files.
- When starting from child restart, the manifest is authoritative for split
  metadata, child extents, frozen intervals, RNG state, and update counters.
- The `<Split>` block is just the XML-facing transport for the existing
  `GaugeSubdomainSplitParams` inputs. It should map directly onto the existing
  split utility rather than introducing alternate split semantics here.

### `SharedHMC`

This block contains the default HMC definition used by both children.

```xml
<SharedHMC>
  <Boundary>
    <GaugeBC>SUBDOMAIN_FIXED_GAUGEBC</GaugeBC>
    <loop_extent>...</loop_extent>
  </Boundary>
  <Monomials>
    ...
  </Monomials>
  <Hamiltonian>
    ...
  </Hamiltonian>
  <MDIntegrator>
    ...
  </MDIntegrator>
</SharedHMC>
```

Required semantics:

- `GaugeBC` is fixed to `SUBDOMAIN_FIXED_GAUGEBC` for this workflow.
- that XML spelling is part of the `subdomain_hmc` contract and does not by
  itself imply first-version standalone factory exposure outside this driver
- `loop_extent` is the maximum temporal reach that the child non-periodic logic
  must honor when deciding which loops are legal near the child edge.
- The driver injects child-local frozen intervals, child-local copied boundary
  links, and the active temporal direction into the helper at runtime.

### `Child0` and `Child1`

Each child block should contain:

```xml
<Child0>
  <ChainControl>
    ...
  </ChainControl>
  <HMCOverrides>
    ...
  </HMCOverrides>
  <InlineMeasurements>
    ...
  </InlineMeasurements>
</Child0>
```

`Child1` has the same shape.

`ChainControl` should intentionally mirror the trajectory-control parts of
`MCControl` from `hmc.cc`, but omit:

- `Cfg`
- embedded inline-measurement XML

The required `ChainControl` fields for the first version are:

- `RNG`
- `StartUpdateNum`
- `NWarmUpUpdates`
- `NProductionUpdates`
- `NUpdatesThisRun`
- `SaveInterval`
- `SavePrefix`
- `SaveVolfmt`
- `ParallelIO`
- `ReproCheckP`
- `ReproCheckFrequency` when enabled
- `ReverseCheckP`
- `ReverseCheckFrequency` when enabled
- `MonitorForces`

Child-specific behavior:

- `child0` and `child1` may have different update counts
- `child0` and `child1` may have different RNG states
- `child0` and `child1` may have different save prefixes and save cadence
- `child0` and `child1` may have different inline-measurement schedules

Override semantics:

- `HMCOverrides` is optional
- if present, it may replace whole `Monomials`, `Hamiltonian`, or
  `MDIntegrator` blocks from `SharedHMC`
- if absent, the child uses the shared definition unchanged

### `Outputs`

The output block should contain stitch and final-output policy.

```xml
<Outputs>
  <StitchAtEnd>true</StitchAtEnd>
  <WriteStitchedParent>true</WriteStitchedParent>
  <StitchedParentCfg>
    ...
  </StitchedParentCfg>
  <KeepChildRestartArtifacts>true</KeepChildRestartArtifacts>
  <PostStitchInlineMeasurements>
    ...
  </PostStitchInlineMeasurements>
</Outputs>
```

Required semantics:

- `StitchAtEnd` is optional and may be `false`
- if `StitchAtEnd` is `false`, the run may end after advancing and checkpointing
  the child chains
- if `StitchAtEnd` is `true`, stitch is attempted only after both child chains
  complete their requested updates
- `WriteStitchedParent` may be `true` only when `StitchAtEnd` is `true`
- `PostStitchInlineMeasurements` is optional and runs on the stitched parent
  field after stitch succeeds
- `KeepChildRestartArtifacts` should default to `true`

## Program Control Flow

The new main program should mirror `hmc.cc` structurally, but replace the
single-chain call to `doHMC(...)` with child-chain orchestration.

Required first-version flow:

1. Parse `SubdomainHMC` input.
2. Register the same broad factory families that `hmc.cc` needs, plus the new
   child-specific pieces required by this workflow.
3. Resolve startup mode:
   - `FROM_PARENT_CFG`: read parent field and split it
   - `FROM_CHILD_RESTART`: load manifest, child gauge fields, child momenta,
     child RNG state, child update counters, and split metadata
4. For each child independently:
   - switch to the child `nrow`
   - build the effective HMC definition from `SharedHMC` plus that child's
     optional overrides
   - construct runtime `SUBDOMAIN_FIXED_GAUGEBC` state using the child's copied
     frozen boundary links and frozen local intervals
   - run the requested warm-up/production updates for that child
   - apply child-local inline measurements on that child schedule
   - write child restart artifacts on that child save cadence
5. After both child chains complete, optionally stitch:
   - restore parent `Layout`
   - call `stitchGaugeSubdomains(...)`
   - abort with a clear error if duplicated frozen-boundary links no longer
     match exactly
   - run optional post-stitch parent measurements
   - write optional stitched parent output
6. On normal completion, preserve the latest child restart state even if stitch
   also ran.

The split and stitch steps above are required to call the already implemented
`splitGaugeSubdomains(...)` and `stitchGaugeSubdomains(...)` APIs, carrying
forward their validation and metadata contract from
`specs/hier/gauge_subdomain_split.md`.

The two child chains are logically independent. The implementation may execute
them sequentially in one process, but the semantics must be those of separate
chains with separate counters, RNG state, measurement schedules, and restart
payloads.

## Exact Responsibility Split

### `SUBDOMAIN_FIXED_GAUGEBC`

This helper is responsible for:

- storing the exact copied child boundary links that must remain fixed
- identifying the two child-local frozen intervals
- reinstalling the fixed boundary links in `modify(Q&)`
- zeroing gauge-like force fields on masked frozen links in `zero(P&)`
- exposing the `loop_extent` value needed by the non-periodic loop logic

It is not sufficient by itself to provide strict child non-periodicity.

### Momentum Masking

The additive momentum-masking layer is responsible for:

- ensuring MD momenta vanish on all links rooted on frozen child slices
- preserving exact frozen-boundary equality across child trajectories
- cooperating with stock HMC machinery without editing stock HMC source in
  place

### Non-Periodic Path Exclusion

The non-periodic child evolution layer is responsible for:

- allowing copied edge-rooted temporal links to exist in the child gauge field
- forbidding any plaquette, rectangle, or longer loop that would step off the
  child lattice and wrap through the periodic `Layout`
- honoring `loop_extent` when determining which temporal-near-edge loops are
  illegal

Operational rule:

- a loop is illegal if evaluating it on the periodic child `Layout` would rely
  on local temporal wrap-around that does not exist in the intended child
  geometry

For `loop_extent = 1`, this excludes offending plaquette construction.
For larger `loop_extent`, the same rule extends to improved-action loops with
longer temporal reach.

### Orchestration Layer

The orchestration layer is responsible for:

- choosing startup mode
- split and stitch calls
- per-child `Layout` switching
- per-child HMC object construction
- per-child measurement execution
- per-child save/restart handling
- optional parent-after-stitch measurements and output

## Restart Artifact Layout

Restart output should use one manifest plus separate per-child payload files.

Manifest contents must include:

- split parameters
- parent `nrow`
- `child0.nrow` and `child1.nrow`
- `child0.local_to_global_t` and `child1.local_to_global_t`
- `child0.frozen_local_intervals` and `child1.frozen_local_intervals`
- child update counters
- child RNG state
- file references for child gauge-field payloads
- file references for child momenta payloads
- enough versioning information to reject incompatible restart layouts cleanly

Separate payload files must store at least:

- `child0` gauge field
- `child1` gauge field
- `child0` momenta
- `child1` momenta

Required behavior:

- the manifest is the authoritative restart entry point
- restart must reconstruct the independent child chains directly, without
  requiring an intermediate stitched parent
- child restart artifacts remain valid whether or not the same run also wrote a
  stitched parent field

## Validation Targets

The first implementation should add focused validation for both the workflow and
the child-edge semantics.

Required checks:

1. Parent-start split/evolve/stitch smoke:
   - start from a parent gauge field
   - evolve both children
   - stitch successfully
   - confirm no frozen-boundary mismatch is reported
2. Child restart round-trip:
   - start from parent
   - evolve children
   - save child restart artifacts
   - restart from manifest
   - continue both child chains without re-splitting from parent
3. Frozen-boundary invariance:
   - verify all links rooted on frozen child slices remain exact across child
     trajectories
4. Independent child scheduling:
   - verify `child0` and `child1` can run different trajectory counts and
     different measurement schedules in one invocation
5. Edge-link path exclusion for plaquettes:
   - verify copied `mu = t_dir` links on the local last child slice do not
     participate in wrap-around plaquette construction
6. Edge-link path exclusion for improved actions:
   - verify `loop_extent > 1` excludes offending longer temporal loops near the
     child edge
7. Stitch failure on duplicated-boundary mismatch:
   - verify stitch aborts clearly if duplicated frozen links are no longer equal

Recommended artifact additions:

- new dedicated main program source under `mainprogs/main/`
- at least one focused test executable for child-edge helper/path semantics
- at least one XML smoke input for `subdomain_hmc`

## Remaining Questions

No additional user-level physics or workflow questions are required to begin the
first implementation from this draft. Any remaining work is implementation
detail, not missing spec intent.

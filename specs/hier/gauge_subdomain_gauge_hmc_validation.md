# Gauge Subdomain Gauge-HMC Validation Spec

## Goal

Define a focused end-to-end validation workflow for the split child-lattice
gauge-only HMC use case:

1. generate one parent gauge field
2. split it into two child lattices
3. run independent gauge-only HMC on the children
4. freeze both force and refreshed momentum on the duplicated temporal
   boundary blocks
5. measure one simple observable on the evolved child configs
6. compare the child results

This is a validation spec for the first child-run experiment. It is not a full
production workflow spec.

## Scope

This first validation intentionally stays narrow:

- gauge-only HMC only
- no fermion monomials
- no SMD path
- deterministic parent start
- symmetric split geometry
- one final-state observable comparison
- one explicit frozen-link invariance check

This spec does not require a stitched-parent physics comparison in the first
version, though that is a natural later extension.

## Current Repo Facts

- The split/stitch geometry and sidecar contract are defined in
  `specs/hier/gauge_subdomain_split.md`.
- `TEMPORAL_ZONE_GAUGEBC` zeroes gauge-like force/update fields on selected
  local time intervals and leaves stored gauge links unchanged.
- The planned modern HMC momentum-masking design is specified separately in
  `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`.
- `mainprogs/main/hmc.cc` already supports gauge-only HMC runs and already
  emits full-lattice plaquette measurements.
- Full-lattice plaquette measurements are not sufficient for an
  updating-region-only comparison.
- `PLAQ_DENSITY` already exists and can write per-root-site plaquette action
  densities for a saved gauge config.

## Chosen Validation Geometry

Use a symmetric split so the two child runs have identical local geometry.

Parent setup:

- `nrow = 4 4 4 8`
- `cfg_type = UNIT`
- `t_dir = 3`
- `cut0 = 0`
- `cut1 = 4`
- `frozen_width = 1`

Derived split plan:

- `F0 = {0}`
- `F1 = {4}`
- `A = {1, 2, 3}`
- `B = {5, 6, 7}`
- `child0_local_to_global_t = [0, 1, 2, 3, 4]`
- `child1_local_to_global_t = [4, 5, 6, 7, 0]`
- `child0_nrow = 4 4 4 5`
- `child1_nrow = 4 4 4 5`
- `child0_frozen_local_intervals = [0,0] and [4,4]`
- `child1_frozen_local_intervals = [0,0] and [4,4]`
- each child updating local interval is `[1,3]`

Rationale:

- Both children have the same local lattice size.
- A `UNIT` parent makes the initial child configs exactly identical.
- Using the same HMC seed on both children makes child-child comparisons very
  sharp.
- Child-child equality by itself is not enough to prove freezing is correct,
  because both children could drift in the same incorrect way. The validation
  must therefore also check that frozen links remain unchanged relative to the
  initial child configs.

## Planned Test Artifacts

Add a small test bundle under `tests/gauge_subdomain_split/`:

- `gauge_subdomain_split.gauge_hmc_validation.ini.xml`
  - split input
  - parent start is `UNIT`
  - writes `child0.validation.scidac`, `child1.validation.scidac`, and one
    sidecar XML file
- `hmc_child0.temporal_zone.ini.xml`
  - gauge-only child HMC input for `child0`
- `hmc_child1.temporal_zone.ini.xml`
  - same as `child0` except for filenames
- `measure_child0_plaq_density.ini.xml`
  - one-shot `chroma` input that reads the evolved `child0` config and writes
    a `PLAQ_DENSITY` file
- `measure_child1_plaq_density.ini.xml`
  - same as `child0` except for filenames

Generated outputs from these inputs are workspace artifacts and should not be
committed.

Optional later helper:

- a small focused reducer or comparison executable may be added if reducing the
  `PLAQ_DENSITY` files manually becomes too awkward

## Child HMC Requirements

The child evolution step should use `mainprogs/main/hmc.cc` with:

- one `GAUGE_MONOMIAL`
- a simple gauge action such as `WILSON_GAUGEACT`
- `TEMPORAL_ZONE_GAUGEBC` under the gauge action
- `zero_intervals` taken exactly from the child-local frozen intervals stored
  in the split sidecar
- `Cfg` pointing to the split child `SCIDAC` file
- `nrow = 4 4 4 5`

For the first validation:

- use the same RNG seed in both child runs
- use the same `beta`, `tau0`, integrator type, and `n_steps`
- keep the run short: a few production updates and one saved final config
- disable reverse-check noise unless it is specifically being debugged

Required implementation assumption:

- the HMC trajectory implements
  `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`

## Observable Definition

The first comparison observable should avoid ambiguity from plaquettes that
straddle the child-local temporal boundary.

Use the final evolved child configs and define the child-local
`active_root_spatial_plaquette`:

1. Run `PLAQ_DENSITY` on the final child config.
2. Interpret its written field as rooted plaquette action densities
   ```c++
   s_{mu,nu}(x) = 1 - (1 / Nc) * ReTr P_{mu,nu}(x)
   ```
3. Restrict root sites `x` to those whose local time coordinate lies in the
   child updating interval.
4. Restrict plaquette planes to purely spatial planes:
   `mu < nu`, `mu != t_dir`, `nu != t_dir`.
5. Average over those sites and planes, then convert back to a plaquette:
   ```c++
   P_active_spatial =
     1 - mean_{x in active roots, spatial planes}(s_{mu,nu}(x))
   ```

For this geometry with `t_dir = 3`, the included spatial planes are:

- `(0,1)`
- `(0,2)`
- `(1,2)`

This observable is intentionally rooted on active slices only. It is simple,
local, and does not mix in temporal plaquettes that necessarily touch the
frozen boundary slice.

## Required Checks

The first validation should require all of the following:

1. The split tool completes successfully and writes the expected sidecar.
2. The sidecar reports child-local frozen intervals `[0,0]` and `[4,4]` for
   both children.
3. Both child HMC runs complete successfully.
4. Frozen-link invariance:
   for each child, all link directions rooted on the frozen local slices in the
   final evolved child config must match the initial split child config exactly.
5. Nontrivial evolution:
   `active_root_spatial_plaquette` on each final child config must differ from
   the unit-field value `1.0` by a detectable nonzero amount.
6. Child-child consistency:
   `active_root_spatial_plaquette(child0)` and
   `active_root_spatial_plaquette(child1)` must agree within a tight tolerance
   for this symmetric same-seed setup.

Optional later stronger checks:

- exact final child-config equality between `child0` and `child1`
- stitch the evolved children back to a parent and compare a parent observable

## Workflow

The intended first validation workflow is:

1. run `gauge_subdomain_split` with the `UNIT` parent input
2. preserve the written child configs as the initial-reference child states
3. run child gauge-only HMC on `child0`
4. run child gauge-only HMC on `child1`
5. run one-shot `PLAQ_DENSITY` measurements on the final `child0` and `child1`
   configs
6. reduce those density outputs over the active local time slices and spatial
   planes
7. compare the two reduced observables
8. separately compare the final child configs against the initial child configs
   on the frozen local slices

The frozen-link invariance check is the direct test of the freezing idea.
The active-region plaquette comparison is the simple child-chain observable
comparison that accompanies it.

## Non-Goals For The First Version

- no long-run statistical study
- no fermion-force validation
- no SMD validation
- no per-update active-region timeseries
- no generalized asymmetric comparison geometry
- no requirement that the first implementation produce a single one-command
  harness

## Open Implementation Dependencies

- Implement `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`.
- Decide whether active-region reduction is done by:
  - a small dedicated helper, or
  - documented postprocessing of `PLAQ_DENSITY` outputs.
- If this validation later moves beyond `UNIT` parents, add a follow-up check
  that compares stitched-parent observables as well as child-local observables.

# HMC Gauge-Monomial Momentum BC Autodiscovery Spec

## Goal

Add a trajectory-level momentum-masking hook to the modern HMC path by
autodiscovering a nontrivial gauge boundary-condition source from the gauge
monomial or monomials actually referenced by the Hamiltonian.

Chosen behavior:

- do not add any new user XML for a separate momentum BC
- use the gauge monomial as the source of truth for momentum masking
- after momentum refresh and `taproj(...)`, apply the same `zero(P&)` behavior
  to the refreshed momentum field before MD integration begins

This keeps force masking and momentum masking coupled to the same gauge-action
boundary-condition configuration.

## Scope

First-version scope:

- `mainprogs/main/hmc.cc`
- `mainprogs/main/const_hmc.cc`
- `lib/update/molecdyn/hmc/lcm_hmc.h`
- `lib/update/molecdyn/hmc/const_lcm_hmc.h`
- named gauge monomials of type
  `Monomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>`
- autodiscovery from `Hamiltonian/monomial_ids`
- no new XML knobs

Explicitly out of scope in the first version:

- `mainprogs/main/smd.cc`
- `lib/update/molecdyn/smd/lcm_smd.h`
- fermion-side BC autodiscovery
- generic frozen-link guarantees for mixed gauge-plus-fermion HMC
- any change to `GaugeBC::modify(Q&)`

Important limitation:

- This hook is necessary for frozen-link behavior, but it is not sufficient by
  itself unless every force term that acts on the same links also honors the
  same masked region.
- For the immediate child-lattice validation use case, the intended target is
  gauge-only HMC.

## Current Repo Facts

- `GaugeBC<P,Q>` already exposes:
  - `zero(P& ds_u) const`
  - `nontrivialP() const`
- `GaugeMonomial` wraps a `GaugeAction`, and the gauge action already owns the
  gauge BC used for force construction.
- `ConstGaugeMonomial` derives from `GaugeMonomial`.
- `readNamedMonomialArray(...)` loads named monomials into `TheNamedObjMap`.
- `ExactHamiltonianParams` reads the list of monomial ids actually used by the
  trajectory from `Hamiltonian/monomial_ids`.
- `LatColMatHMCTrj::refreshP(...)` and `ConstLatColMatHMCTrj::refreshP(...)`
  currently refresh momenta and project with `taproj(...)`, but do not apply a
  BC mask afterward.
- The legacy HMC path did explicitly zero refreshed momenta with the gauge BC.

## Design Principles

- No duplicate user configuration for a separate momentum BC.
- Discover the masking source only from monomials actually referenced by the
  Hamiltonian, not from every named monomial defined in the input.
- Keep trajectory code simple: the trajectory applies a discovered mask source,
  but it does not perform XML or named-object discovery itself.
- Reuse the gauge monomial itself as the masking source instead of extending
  `GaugeAction`, `CreateGaugeState`, or `ExactHamiltonian` handle plumbing just
  to recover a `GaugeBC` handle.
- Be conservative when multiple nontrivial gauge monomials are present:
  accept only if their `zero(P&)` behavior is compatible for the current
  lattice.

## No XML Changes

The first implementation intentionally changes no user XML surface:

- no new `<MomentumBC>` block
- no new `HMCTrj` parameter
- no new `Hamiltonian` parameter
- no new `Monomials` parameter

The gauge monomial configuration remains the only source of truth.

## Public API Additions

Add small public query/apply helpers to `GaugeMonomial` in
`lib/update/molecdyn/monomial/gauge_monomial.h`:

```c++
const GaugeBC<P,Q>& getGaugeBC() const;
bool hasNontrivialGaugeBC() const;
void zeroGaugeLikeField(P& field) const;
```

Semantics:

- `getGaugeBC()` forwards to the wrapped gauge action
- `hasNontrivialGaugeBC()` forwards to `getGaugeBC().nontrivialP()`
- `zeroGaugeLikeField(...)` forwards to `getGaugeBC().zero(...)`

These are read-only helpers. They do not change monomial force construction.

`ConstGaugeMonomial` inherits these helpers through `GaugeMonomial`.

## Discovery Helper

Add one small helper under the HMC code path, for example:

- `lib/update/molecdyn/hmc/gauge_monomial_momentum_bc.h`
- `lib/update/molecdyn/hmc/gauge_monomial_momentum_bc.cc`

Suggested signature:

```c++
Handle<GaugeMonomial>
discoverMomentumMaskingGaugeMonomial(const multi1d<std::string>& monomial_ids);
```

Contract:

- inspect only the monomial ids supplied in `Hamiltonian/monomial_ids`
- return an empty handle if no nontrivial gauge BC source is found
- return one compatible `GaugeMonomial` handle if one or more compatible
  nontrivial gauge BC sources are found
- abort with a clear message if incompatible nontrivial gauge BC sources are
  found

## Discovery Algorithm

Use the following algorithm:

1. Iterate over `monomial_ids` in Hamiltonian order.
2. For each id, look up the named monomial handle in `TheNamedObjMap`.
3. Try to cast the handle to `Handle<GaugeMonomial>`.
4. If the cast fails, ignore that monomial and continue.
5. If the cast succeeds but `hasNontrivialGaugeBC()` is false, ignore it and
   continue.
6. If this is the first nontrivial gauge monomial found, store it as the
   candidate masking source.
7. For each later nontrivial gauge monomial, compare its `zero(P&)` behavior
   against the candidate masking source on a deterministic probe field.
8. If the resulting masked probe fields differ, abort with a clear error that
   names both monomial ids.
9. If the masked probe fields match, treat the two monomials as compatible and
   keep the original candidate handle.
10. If no nontrivial gauge monomial was found, return an empty handle.

### Compatibility Test

The compatibility test should compare actual `zero(P&)` behavior, not XML text
or pointer identity.

Suggested probe field:

```c++
multi1d<LatticeColorMatrix> probe(Nd);
for (int mu = 0; mu < Nd; ++mu) {
  probe[mu] = Real(mu + 1);
}
```

Then:

```c++
probe_a = probe;
probe_b = probe;

candidate->zeroGaugeLikeField(probe_a);
other->zeroGaugeLikeField(probe_b);
```

Compatibility requirement:

- `probe_a` and `probe_b` must match exactly on all link directions

Why compare `zero(P&)` behavior directly:

- this is the exact operation the trajectory will apply
- it is generic across gauge BC classes
- it avoids introducing BC-handle plumbing just to compare configuration
- it avoids fragile pointer-identity checks for semantically equivalent BCs

## Trajectory Changes

Extend the modern HMC trajectory classes to accept one optional masking source.

### `LatColMatHMCTrj`

Update the constructor in `lib/update/molecdyn/hmc/lcm_hmc.h` to accept:

```c++
Handle<GaugeMonomial> momentum_mask_source
```

Store it as a private member.

After the existing refresh loop and `taproj(...)`, add:

```c++
if (!momentum_mask_source.null() &&
    momentum_mask_source->hasNontrivialGaugeBC()) {
  momentum_mask_source->zeroGaugeLikeField(s.getP());
}
```

Apply the mask once to the full `multi1d<LatticeColorMatrix>` momentum field,
after all directions have been refreshed and projected.

### `ConstLatColMatHMCTrj`

Apply the same constructor and post-refresh masking pattern in
`lib/update/molecdyn/hmc/const_lcm_hmc.h`.

## Driver Wiring

### `mainprogs/main/hmc.cc`

After:

1. `readNamedMonomialArray(...)`
2. `ExactHamiltonianParams ham_params(...)`

call:

```c++
Handle<GaugeMonomial> momentum_mask_source =
  discoverMomentumMaskingGaugeMonomial(ham_params.monomial_ids);
```

Then construct:

```c++
LatColMatHMCTrj theHMCTrj(H_MC, Integrator, momentum_mask_source);
```

### `mainprogs/main/const_hmc.cc`

Use the same discovery helper and pass the result into
`ConstLatColMatHMCTrj`.

## Error Behavior

Abort with a clear message if:

- a monomial id listed in the Hamiltonian cannot be found in `TheNamedObjMap`
- a named monomial handle lookup fails unexpectedly
- two or more nontrivial gauge monomials referenced by the Hamiltonian produce
  different `zero(P&)` results on the compatibility probe field

Do not abort when:

- there are no gauge monomials in the Hamiltonian
- all discovered gauge monomials have trivial gauge BCs
- nongauge monomials appear in the Hamiltonian

In those cases, return an empty handle and preserve existing behavior.

## Why Not Extend `ExactHamiltonian`?

This first implementation intentionally avoids making `ExactHamiltonian`
responsible for BC discovery.

Reasons:

- `ExactHamiltonian` currently binds exact monomials for energy and internal
  field refresh, but it does not expose gauge-specific structure
- the needed discovery inputs are already available in the driver:
  - named monomials have already been instantiated
  - Hamiltonian monomial ids have already been parsed
- the trajectory only needs a masking source, not full Hamiltonian-level BC
  semantics

This keeps the first patch smaller and more local.

## Tests

Add one focused code test and one trajectory smoke test.

### Focused Discovery Test

Add a small executable, for example:

- `mainprogs/tests/t_hmc_momentum_bc_autodiscovery.cc`

Suggested checks:

1. One gauge monomial with `PERIODIC_GAUGEBC`
   - discovery returns empty
2. One gauge monomial with `TEMPORAL_ZONE_GAUGEBC`
   - discovery returns non-empty
   - applying `zeroGaugeLikeField(...)` matches the expected masked region
3. Two compatible nontrivial gauge monomials
   - discovery succeeds
4. Two incompatible nontrivial gauge monomials
   - discovery aborts with a clear error

The focused test may construct monomials through existing XML/factory paths and
use the named object map exactly as `hmc.cc` does.

### HMC Smoke Test

Add a small gauge-only HMC input under `tests/hmc/`, for example:

- `tests/hmc/hmc.wilson.temporal_zone_gaugebc.ini.xml`

Suggested setup:

- `cfg_type = UNIT`
- one `GAUGE_MONOMIAL`
- `WILSON_GAUGEACT`
- `TEMPORAL_ZONE_GAUGEBC`
- a short trajectory
- save the final gauge config

Suggested validation:

1. links rooted on frozen slices remain exactly unchanged from the initial unit
   config
2. at least one rooted link on an active slice changes

The higher-level split-child workflow in
`specs/hier/gauge_subdomain_gauge_hmc_validation.md` then supplies the
end-to-end child-lattice integration coverage.

## Build System Updates

If a new helper source/header pair is added under `lib/`, update both:

- `lib/Makefile.am`
- `lib/CMakeLists.txt`

If a new focused test executable is added, update both:

- `mainprogs/tests/Makefile.am`
- `mainprogs/tests/CMakeLists.txt`

## Validation Checklist

1. Build the updated library and HMC drivers.
2. Verify existing periodic gauge-only HMC behavior is unchanged when no
   nontrivial gauge BC is discovered.
3. Verify autodiscovery finds the temporal-zone gauge monomial for the HMC
   smoke input.
4. Verify refreshed momenta are zeroed on the frozen slices before the first
   gauge update.
5. Verify frozen rooted links remain unchanged after the trajectory.
6. Verify active rooted links do evolve.
7. Verify the child-lattice validation workflow can depend on this hook without
   adding a separate momentum-BC XML block.

## Follow-On Work

- Apply the same discovery helper and post-refresh masking pattern to SMD.
- If mixed gauge-plus-fermion child HMC becomes a target, define the fermion
  force-side compatibility story explicitly.
- If future use cases need broader BC introspection, reconsider whether
  `ExactHamiltonian` should expose discovered gauge-masking metadata directly.

# AGENTS.md

This file is the working handoff for humans and coding agents contributing to this checkout. Update it whenever the build recipe, active feature work, or test strategy changes.

## Project Snapshot

- Chroma is the Jefferson Lab lattice QCD codebase in this repository.
- This checkout supports both legacy Autotools and CMake. In current local work, the CMake path is the most reliable documented route.
- Chroma depends on QDP++/QDPXX. For this checkout, the tested local recipe uses a QDPXX worktree based on `origin/eloy/localbinarydb`.

## Repo Map

- `lib/`: core library code. Gauge boundary-condition work lives under `lib/actions/gauge/gaugebcs/`. Gauge split/stitch utility work lives under `lib/util/gauge/`.
- `mainprogs/main/`: production executables.
- `mainprogs/tests/`: focused test executables and small integration checks.
- `tests/`: XML inputs, reference outputs, and regression fixtures.
- `docs/`: user and developer documentation. The current local build notes live in `docs/build_local_qdpxx.md`.
- `specs/hier/`: design notes and implementation specs for in-flight features.
- `scripts/`: build and bootstrap helpers. The current tested local bootstrap is `scripts/bootstrap_local_qdpxx_build.sh`.
- `other_libs/`: bundled or submodule-backed dependencies used by Chroma builds.
- `build/`: local build and dependency scratch space in this checkout. Treat it as workspace state, not primary source, unless you are intentionally changing the bootstrap flow.

## Current Workstream

The current feature work is centered on two related gauge-only temporal-boundary efforts:

- force-suppression gauge boundary condition:
  - public XML/factory name: `TEMPORAL_ZONE_GAUGEBC`
  - implementation: `lib/actions/gauge/gaugebcs/temporal_zone_gaugebc.{h,cc}`
  - factory wiring: `lib/actions/gauge/gaugebcs/gaugebc_aggregate.cc` and `lib/actions/gauge/gaugebcs/gaugebcs.h`
  - focused test executable: `mainprogs/tests/t_temporal_zone_gaugebc.cc`
  - smoke input: `tests/t_leapfrog/t_leapfrog.temporal_zone_gaugebc.ini.xml`
  - spec: `specs/hier/temporal_zone_gaugebc.md`
- temporal gauge split/stitch utility:
  - public utility names: `splitGaugeSubdomains(...)` and `stitchGaugeSubdomains(...)`
  - implementation: `lib/util/gauge/gauge_subdomain_split.{h,cc}`
  - focused test executable: `mainprogs/tests/t_gauge_subdomain_split.cc`
  - spec: `specs/hier/gauge_subdomain_split.md`
- subdomain quenched HMC driver:
  - public executable name: `subdomain_hmc`
  - main implementation: `mainprogs/main/subdomain_hmc.cc`
  - child BC helper: `SUBDOMAIN_FIXED_GAUGEBC`
  - implementation: `lib/actions/gauge/gaugebcs/subdomain_fixed_gaugebc.{h,cc}`
  - focused helper test: `mainprogs/tests/t_subdomain_fixed_gaugebc.cc`
  - smoke input: `tests/subdomain_hmc/subdomain_hmc.wilson.ini.xml`
  - spec: `specs/hier/gauge_subdomain_hmc.md`

Important behavior note:

- `TEMPORAL_ZONE_GAUGEBC` currently uses `GaugeBC::zero(P&)` to suppress gauge-like force/update fields on the selected time intervals only.
- It does not modify the stored gauge links in `modify(Q&)`.
- `gauge_subdomain_split` currently targets an in-memory gauge-only split/stitch workflow and depends on layout switching. Child gauge fields are stored as `multi1d<LatticeColorMatrix>` with metadata describing their child-local temporal ordering.

## Build And Test

Preferred local recipe for this checkout:

```bash
./scripts/bootstrap_local_qdpxx_build.sh all
```

Useful modes:

```bash
./scripts/bootstrap_local_qdpxx_build.sh bootstrap
./scripts/bootstrap_local_qdpxx_build.sh build-tests
./scripts/bootstrap_local_qdpxx_build.sh run-tests
```

Notes about the current recipe:

- The source path contains spaces, so the bootstrap script builds through a no-space alias at `/private/tmp/chroma-ws`.
- The tested QDPXX branch for this checkout is `origin/eloy/localbinarydb`.
- The helper script may patch the temporary QDPXX worktree for Apple clang compatibility by adding `<array>` to `include/qdp_map_obj_disk.h`.
- Detailed manual commands and expected outputs are documented in `docs/build_local_qdpxx.md`.

Current verification targets:

- `t_temporal_zone_gaugebc` should exit successfully.
- `t_gauge_subdomain_split` should exit successfully.
- `t_leapfrog` should accept `TEMPORAL_ZONE_GAUGEBC` and complete using `tests/t_leapfrog/t_leapfrog.temporal_zone_gaugebc.ini.xml`.
- `t_subdomain_fixed_gaugebc` should exit successfully.
- `subdomain_hmc` should complete the parent-start split/evolve/stitch/save smoke using `tests/subdomain_hmc/subdomain_hmc.wilson.ini.xml`.

Current known gaps in the new subdomain-HMC work:

- post-stitch parent inline measurements are still unstable; the committed smoke input leaves `PostStitchInlineMeasurements` empty for now.
- `FROM_CHILD_RESTART` is not yet passing end-to-end; resume currently fails while reading the child restart manifest payloads.

Legacy build context still matters:

- `README` documents the traditional out-of-tree Autotools flow with `../configure --with-qdp=...` and `make`.
- If you touch source lists or install-facing build behavior, keep both the Autotools and CMake paths in sync.

## Editing Rules For This Repo

- When adding or removing library sources, update both `lib/Makefile.am` and `lib/CMakeLists.txt`.
- When adding or removing test executables, update both `mainprogs/tests/Makefile.am` and `mainprogs/tests/CMakeLists.txt`.
- Factory-registered features should be wired into the relevant aggregate registration point and umbrella include header, not only added as standalone files.
- XML-facing features should usually ship with both a focused executable test and a smoke or regression XML input under `tests/`.
- Utilities that switch `Layout` internally should document which lattice extent must be active before their returned `Lattice*` objects are used.
- If the intended behavior is still being reasoned through, capture it in `specs/hier/*.md` rather than leaving it only in chat history.
- Prefer changing real source under `lib/`, `mainprogs/`, `tests/`, `docs/`, and `specs/`. Avoid editing `build/deps/src/qdpxx` unless you are intentionally changing the local dependency workflow.

## Working Tree Hygiene

- Check `git status` before staging. This checkout can accumulate local artifacts that are not fully ignored.
- Common local-only artifacts seen during current work include `build/`, `XMLDAT`, `.DS_Store`, and `.vscode/`.
- Keep generated outputs and temporary workspace files out of commits unless the change is explicitly about the bootstrap or test workflow itself.

## How To Update This File

When work continues, keep these sections current:

- `Current Workstream`: what feature or subsystem is active right now
- `Build And Test`: the latest known-good commands and fragile assumptions
- `Editing Rules For This Repo`: any newly discovered wiring requirements
- `Working Tree Hygiene`: new generated artifacts or ignore-list gaps worth remembering

If a future effort replaces the local QDPXX recipe, update this file and `docs/build_local_qdpxx.md` together.

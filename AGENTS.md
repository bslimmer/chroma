# AGENTS.md

This file is the working handoff for humans and coding agents contributing to this checkout. Update it whenever the build recipe, active feature work, or test strategy changes.

## Project Snapshot

- Chroma is the Jefferson Lab lattice QCD codebase in this repository.
- This checkout supports both legacy Autotools and CMake. In current local work, the CMake path is the most reliable documented route.
- Chroma depends on QDP++/QDPXX. For this checkout, the tested local recipe uses a QDPXX worktree based on `origin/eloy/localbinarydb`.

## Repo Map

- `lib/`: core library code. Gauge boundary-condition work lives under `lib/actions/gauge/gaugebcs/`.
- `mainprogs/main/`: production executables.
- `mainprogs/tests/`: focused test executables and small integration checks.
- `tests/`: XML inputs, reference outputs, and regression fixtures.
- `docs/`: user and developer documentation. The current local build notes live in `docs/build_local_qdpxx.md`.
- `specs/hier/`: design notes and implementation specs for in-flight features.
- `scripts/`: build and bootstrap helpers. The current tested local bootstrap is `scripts/bootstrap_local_qdpxx_build.sh`.
- `other_libs/`: bundled or submodule-backed dependencies used by Chroma builds.
- `build/`: local build and dependency scratch space in this checkout. Treat it as workspace state, not primary source, unless you are intentionally changing the bootstrap flow.

## Current Workstream

The current feature work has four closely related threads:

1. A new force-suppression gauge boundary condition:
   - public XML/factory name: `TEMPORAL_ZONE_GAUGEBC`
   - implementation: `lib/actions/gauge/gaugebcs/temporal_zone_gaugebc.{h,cc}`
   - factory wiring: `lib/actions/gauge/gaugebcs/gaugebc_aggregate.cc` and `lib/actions/gauge/gaugebcs/gaugebcs.h`
   - focused test executable: `mainprogs/tests/t_temporal_zone_gaugebc.cc`
   - smoke input: `tests/t_leapfrog/t_leapfrog.temporal_zone_gaugebc.ini.xml`
   - spec: `specs/hier/temporal_zone_gaugebc.md`
   - HMC momentum-hook spec: `specs/hier/hmc_gauge_momentum_bc_autodiscovery.md`
   - HMC momentum-hook helper: `lib/update/molecdyn/hmc/gauge_monomial_momentum_bc.{h,cc}`
   - HMC/const-HMC wiring: `mainprogs/main/{hmc,const_hmc}.cc` and `lib/update/molecdyn/hmc/{lcm_hmc,const_lcm_hmc}.h`
   - focused test executable: `mainprogs/tests/t_hmc_momentum_bc_autodiscovery.cc`

2. A gauge subdomain split/stitch workflow for separate child-lattice runs:
   - spec: `specs/hier/gauge_subdomain_split.md`
   - child gauge-HMC validation draft: `specs/hier/gauge_subdomain_gauge_hmc_validation.md`
   - child topology validation draft: `specs/hier/gauge_subdomain_topology_validation.md`
   - `0-+` correlator draft: `specs/hier/glueball_0mp_correlator.md`
     - first-version plan reuses inline `QACTDEN` plus offline reduction, with
       no new inline glue measurement required
   - two-level `0-+` workflow spec:
     `specs/hier/gauge_subdomain_two_level_0mp_workflow.md`
     - first target is a level-0 parent-boundary ensemble plus level-1 child
       conditional averages, before any production implementation or tuning
   - two-level `0++` workflow draft:
     `specs/hier/gauge_subdomain_two_level_0pp_workflow.md`
     - first-version plan reuses the existing scalar glueball / spatial
       Wilson-loop operator path with `support_guard = 0`
     - first selected operator is the existing-infrastructure `2 x 2` square
       blocked plaquette (`bl_level_selected = 1`)
     - first implementation assumes a thin post-HMC measurement pass on saved
       parent and child configs rather than a new inline HMC emitter
     - first concrete target is an even split `8^4` geometry with a nominal
       `40` retained parent samples and `10` retained post-discard,
       post-thinning child measurements per child domain per retained parent
       sample
   - periodic / two-level `0++` correlator checker:
     `mainprogs/tests/t_glueball_0pp_corr.cc`
     - measures blocked spatial `0++` operators from saved configs via the
       existing `gluecor` path and then reduces them offline across the
       periodic, parent-window, child-summary, cross-domain, and
       outer-ensemble stages
   - standalone thin gluecor CSV tool:
     `mainprogs/main/gluecor_measure.cc`
     - reads one saved gauge config, applies the usual blocked-link gluecor
       path, and writes the underlying `op0[t]` array as a simple CSV
     - defaults to `decay_dir = 3`, `bl_level = 1`, `BlkAccu = 1.0e-5`, and
       `BlkMax = 50`, with command-line overrides for general use
   - local directory launcher:
     `tests/two_level_glueball/gluecor_measure_launch.sh`
     - walks one directory of `.lime` or `.scidac` configs and launches the
       standalone `gluecor_measure` executable once per file
     - expects one lattice geometry per directory and takes that geometry from
       `NROW="n0 n1 n2 n3"` in the environment
     - skips any child config whose matching
       `gluecor_output/<cfg_stem>.gluecor.csv` already exists
     - now emits both the existing raw correlators and connected
       `V_s^{-2}`-normalized correlators in its measurement summaries and CSV
       outputs
     - the tracked `8^4` workflow now compares all available cross-domain
       separations for `support_guard = 0`, namely `delta_t_parent = 2..6`
     - `TWO_LEVEL_OUTER_ENSEMBLE` also accepts
       `Reducer/require_pass = false` for low-stat thin validation passes where
       writing the summary matters more than enforcing the final `3 sigma`
       consistency gate
   - first periodic `0++` XML bundle:
     `tests/glueball_0pp/measure_unit_glueball_0pp.ini.xml`,
     `tests/glueball_0pp/measure_full_lattice.glueball_0pp.template.ini.xml`,
     `tests/glueball_0pp/glueball_0pp_corr.unit.check.ini.xml`, and
     `tests/glueball_0pp/glueball_0pp_corr.full_lattice.check.ini.xml`
   - first two-level `0++` template bundle:
     `tests/gauge_subdomain_split/hmc_parent.0pp_two_level_outer.template.ini.xml`,
     `tests/gauge_subdomain_split/gauge_subdomain_split.0pp_two_level.template.ini.xml`,
     `tests/gauge_subdomain_split/hmc_child.temporal_zone_glueball_0pp_2lvl.template.ini.xml`,
     `tests/gauge_subdomain_split/measure_glueball_0pp_parent.template.ini.xml`,
     `tests/gauge_subdomain_split/measure_glueball_0pp_child.template.ini.xml`,
     `tests/gauge_subdomain_split/glueball_0pp_parent_window.template.check.ini.xml`,
     `tests/gauge_subdomain_split/glueball_0pp_child.template.check.ini.xml`,
     `tests/gauge_subdomain_split/glueball_0pp_outer_sample.template.check.ini.xml`,
     and generator
     `tests/gauge_subdomain_split/generate_two_level_0pp_xml_bundle.sh`
     - first higher-statistics preset targets the first `8^4` evenly split
       production plan with `40` retained outer samples and `10` retained
       post-thinning child measurements per child stream
     - the generator also supports reusing an aligned retained-parent prefix
       from an existing `cfgs/...` directory via `--reuse-parent-root` and
       `--reuse-parent-last-update`, then generating only the missing later
       parent samples from the last reused saved gauge config
     - current provisional `8^4` HMC tuning uses `tau0 = 0.1`,
       `n_steps = 15` for the parent `8^4` stream and `tau0 = 0.1`,
       `n_steps = 6` for the child `8^3 x 5` fixed-boundary streams,
       targeting about `80%` acceptance
   - periodic `QACTDEN` correlator checker: `mainprogs/tests/t_qactden_0mp_corr.cc`
   - first periodic XML bundle: `tests/glueball_0mp/hmc_full_lattice.qactden_0mp.ini.xml`, `tests/glueball_0mp/qactden_0mp_corr.full_lattice.check.ini.xml`, and `tests/glueball_0mp/measure_unit_qactden_0mp.ini.xml`
   - first two-level XML bundle: `tests/gauge_subdomain_split/hmc_parent.0mp_two_level_outer.ini.xml`, `tests/gauge_subdomain_split/gauge_subdomain_split.0mp_two_level.outer_{100,200}.ini.xml`, `tests/gauge_subdomain_split/hmc_child{0,1}.temporal_zone_qactden_0mp_2lvl.outer_{100,200}.ini.xml`, `tests/gauge_subdomain_split/qactden_0mp_parent_window.outer_{100,200}.check.ini.xml`, `tests/gauge_subdomain_split/qactden_0mp_child{0,1}.outer_{100,200}.check.ini.xml`, `tests/gauge_subdomain_split/qactden_0mp_outer_sample.outer_{100,200}.check.ini.xml`, and `tests/gauge_subdomain_split/qactden_0mp_two_level.check.ini.xml`
   - higher-statistics two-level template bundle: `tests/gauge_subdomain_split/hmc_parent.0mp_two_level_outer.template.ini.xml`, `tests/gauge_subdomain_split/gauge_subdomain_split.0mp_two_level.template.ini.xml`, `tests/gauge_subdomain_split/hmc_child.temporal_zone_qactden_0mp_2lvl.template.ini.xml`, `tests/gauge_subdomain_split/qactden_0mp_parent_window.template.check.ini.xml`, `tests/gauge_subdomain_split/qactden_0mp_child.template.check.ini.xml`, `tests/gauge_subdomain_split/qactden_0mp_outer_sample.template.check.ini.xml`, and generator `tests/gauge_subdomain_split/generate_two_level_0mp_xml_bundle.sh`
     - first higher-statistics preset targets the former one-hour plan with `12` outer samples and `32` retained child measurements per child stream
   - library implementation: `lib/util/gauge/gauge_subdomain_split.{h,cc}`
   - user-facing tools: `mainprogs/main/gauge_subdomain_split.cc` and `mainprogs/main/gauge_subdomain_stitch.cc`
   - focused test executable: `mainprogs/tests/t_gauge_subdomain_split.cc`
   - example tool inputs: `tests/gauge_subdomain_split/gauge_subdomain_split.ini.xml` and `tests/gauge_subdomain_split/gauge_subdomain_stitch.ini.xml`
   - standalone child-sized HMC smoke input: `tests/gauge_subdomain_split/hmc_temporal_zone_child_smoke.ini.xml`
   - child gauge-HMC validation bundle: `tests/gauge_subdomain_split/gauge_subdomain_split.gauge_hmc_validation.ini.xml`, `tests/gauge_subdomain_split/hmc_child{0,1}.temporal_zone.ini.xml`, `tests/gauge_subdomain_split/measure_child{0,1}_plaq_density.ini.xml`, and `tests/gauge_subdomain_split/gauge_subdomain_gauge_hmc_validation.check.ini.xml`
   - child topology validation bundle: `tests/gauge_subdomain_split/hmc_parent.topology_warmup.ini.xml`, `tests/gauge_subdomain_split/gauge_subdomain_split.topology_validation.ini.xml`, and `tests/gauge_subdomain_split/hmc_child{0,1}.temporal_zone_qactden.ini.xml`
   - focused validation checker: `mainprogs/tests/t_gauge_subdomain_gauge_hmc_validation.cc`
   - intended workflow: split one parent config into two ordinary child configs, evolve the children in separate runs with frozen temporal boundary intervals, then stitch them back into a parent config
   - first-version assumption: user-facing split/stitch tools write QIO outputs and persist split metadata in a sidecar XML file
   - current validation assumption: child gauge-only HMC must freeze both force and refreshed momentum on the duplicated boundary intervals

3. A factorized distillation propagator-writing workflow built on split child
   HMC infrastructure:
   - draft spec outline: `specs/hier/factorized_prop.md`
   - intended first measurement surface:
     `lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc`
   - current draft goal: define the frozen-region / unfrozen-region block
     decomposition, the exact identity-inserted cross-subdomain propagator
     writing with left/right tensors and the assembled cross-subdomain
     propagator, the frozen-boundary XML input path, and the child-local
     output objects before any meson two-point or two-level estimator spec is
     added
   - current first-version decisions:
     `Gamma*` is the union of both frozen temporal boundary blocks, frozen
     boundary XML enters through `NamedObject`, and the assembled
     `D_cross^{-1}` output is written separately from the `L` / `R`
     intermediates
   - overlapping-domain (Ce/Giusti/Schaefer) variant of the parent-lattice
     test path: `specs/hier/factorized_prop_overlapping.md`, enabled with
     `NamedObject/overlapping_factorization = true`; projector on the first
     internal frozen slice from the source side (`c1+1`, `fw-2`), two cut-field
     solvers (`u_src`, `u_snk`); `num_intermediate_projectors <= 0` means no
     projector in this mode

4. Child-local factorized baryon tensors (LEFT sink bracket / RIGHT source
   bracket of the factorized baryon two-point function):
   - spec: `specs/hier/factorized_baryon_tensors.md`
   - measurement `FACTORIZED_BARYON_TENSOR_DISTILLATION_SUPERB`:
     `lib/meas/inline/hadron/inline_factorized_baryon_tensor_distillation_superb_w.{h,cc}`
     (registered in `inline_hadron_aggregate.cc` under `BUILD_SB`)
   - shared helpers (geometry, cut links, parent temporal phase, parent-file
     projector reader), also used by the prototype:
     `lib/util/ferm/factorized_distillation.{h,cc}`
   - projectors live on frozen slices adjacent to the source child's active
     region and are read from one parent-lattice colorvec file
     (`projector_colorvec_files`, MOD or fingerprint s3t); both sides are
     solved from the projector sources; all XML times are parent-global
   - the fermion state must not smear in time (`orthog_dir` = decay dir)
   - checker: `mainprogs/tests/t_factorized_baryon_tensor.cc`
   - end-to-end validation: `tests/factorized_baryon_tensor/run_validation.sh`
     (outputs under `cfgs/factorized_baryon_tensor/`)

Important behavior note:

- `TEMPORAL_ZONE_GAUGEBC` currently uses `GaugeBC::zero(P&)` to suppress gauge-like force/update fields on the selected time intervals only.
- It does not modify the stored gauge links in `modify(Q&)`.
- Modern `hmc` and `const_hmc` now autodiscover a compatible gauge-monomial BC source from `Hamiltonian/monomial_ids` and apply its `zero(P&)` mask to refreshed momenta after `taproj(...)`, with no separate momentum-BC XML.

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

Higher-statistics two-level XML generation:

```bash
./tests/gauge_subdomain_split/generate_two_level_0mp_xml_bundle.sh
./tests/gauge_subdomain_split/generate_two_level_0pp_xml_bundle.sh
./tests/gauge_subdomain_split/generate_two_level_0pp_xml_bundle.sh \
  --reuse-parent-root cfgs/two_level_0pp_10x10_baseline_jul03 \
  --reuse-parent-last-update 140 \
  cfgs/two_level_0pp_extend 50 1040 10 10 0 1 1
```

Notes about the current recipe:

- The source path contains spaces, so the bootstrap script builds through a no-space alias at `/private/tmp/chroma-ws`.
- The tested QDPXX branch for this checkout is `origin/eloy/localbinarydb`.
- Chroma's CMake build requires `cmake 3.17+`; older cluster modules such as
  `3.13.x` should be treated as unsupported for this branch.
- Older Linux cluster GCC defaults may compile QIO in pre-C99 mode; the local
  QDPXX helper now defaults `QDPXX_CFLAGS` to `-O2 -std=gnu99` so
  `other_libs/qio/lib/dml/DML_utils.c` builds cleanly.
- `./scripts/bootstrap_local_qdpxx_build.sh build-tests` now also builds
  `mainprogs/main/gluecor_measure` in addition to the focused test and smoke
  targets.
- If that ref already exists locally and network fetches are unavailable, set `QDPXX_SKIP_FETCH=1`; the bootstrap script will reuse the local ref, repopulate QDPXX `other_libs/` from the local source checkout when needed, and reuse an already-populated `other_libs/qdp-lapack` tree without mutating `.git/modules`.
- If a fork-based clone resolves `other_libs/qdp-lapack` against the wrong
  GitHub namespace, either add an `upstream` remote for the canonical Chroma
  repo or set `QDP_LAPACK_URL` explicitly before running the bootstrap helper.
- The bootstrap helper now falls back to `cd ... && git ...`, a plain local
  QDPXX clone, and `git checkout -B` when a cluster Git is too old for
  `git -C`, `git worktree`, or `git switch`.
- The helper script may patch the temporary QDPXX worktree for Apple clang compatibility by adding `<array>` to `include/qdp_map_obj_disk.h`.
- The current Chroma tree also includes the `Serializable::serialID()` return-type compatibility fix in `lib/util/ferm/key_val_db.h` required by the tested QDPXX/filedb combination.
- Detailed manual commands and expected outputs are documented in `docs/build_local_qdpxx.md`.

Current verification targets:

- `t_temporal_zone_gaugebc` should exit successfully.
- `t_gauge_subdomain_split` should exit successfully.
- `t_gauge_subdomain_gauge_hmc_validation` should complete successfully after running the split plus paired child-HMC validation workflow.
- `t_qactden_0mp_corr` should build successfully and validate the periodic `QACTDEN` timeslice reduction using `tests/glueball_0mp/qactden_0mp_corr.full_lattice.check.ini.xml` once the corresponding HMC XML has been produced.
- `t_qactden_0mp_corr` should also support the two-level parent-window, child-summary, cross-domain, and outer-ensemble reducer modes using the `tests/gauge_subdomain_split/qactden_0mp_*.check.ini.xml` bundle once the corresponding parent and child HMC XML logs have been produced under `cfgs/two_level_0mp/` or from a generated template bundle such as `cfgs/two_level_0mp_1h/`.
- `t_glueball_0pp_corr` should build successfully, pass the unit-gauge
  measurement smoke `tests/glueball_0pp/measure_unit_glueball_0pp.ini.xml`,
  and validate the periodic blocked-plaquette reduction using
  `tests/glueball_0pp/glueball_0pp_corr.unit.check.ini.xml`.
- `t_glueball_0pp_corr` should also support the full-lattice periodic reducer
  using `tests/glueball_0pp/glueball_0pp_corr.full_lattice.check.ini.xml`
  after a saved-config measurement pass has produced the expected summary file.
- `t_glueball_0pp_corr` should also support the two-level parent-window,
  child-summary, cross-domain, and outer-ensemble reducer modes using the
  `tests/gauge_subdomain_split/glueball_0pp_*.check.ini.xml` bundle once the
  corresponding parent and child measurement summaries have been produced under
  a generated template bundle such as `cfgs/two_level_0pp_40x10/`.
- `t_glueball_0pp_corr` should also support a geometry-matched thin validation
  pass against existing saved split configs, such as the local
  `cfgs/two_level_0mp/` run, where a low-stat outer-ensemble summary may use
  `Reducer/require_pass = false`.
- `gluecor_measure` should build successfully and accept
  `--nrow n0 n1 n2 n3 <cfg_file> <csv_file>` for thin saved-config `op0[t]`
  dumps without routing through `t_glueball_0pp_corr`.
- `t_hmc_momentum_bc_autodiscovery` should exit successfully.
- `t_leapfrog` should accept `TEMPORAL_ZONE_GAUGEBC` and complete using `tests/t_leapfrog/t_leapfrog.temporal_zone_gaugebc.ini.xml`.
- `hmc` should complete a one-update child-sized gauge-only smoke run using `tests/gauge_subdomain_split/hmc_temporal_zone_child_smoke.ini.xml`.
- `tests/factorized_baryon_tensor/run_validation.sh` (needs a superbblas +
  PRIMME build with `chroma`, `gauge_subdomain_split` and
  `t_factorized_baryon_tensor`; set `CHROMA_BUILD`) should end with "All
  factorized baryon tensor validation steps passed", for both
  `FROZEN_WIDTH=1` and `FROZEN_WIDTH=2`.
- Known QDPXX issue: with the tested `localbinarydb` QDPXX branch, writing
  FILEDB (MOD) colorvec files (`CREATE_COLORVECS_SUPERB` with
  `write_fingerprint = false`) aborts with SIGFPE in
  `MapObjectDisk::convertToPrivate` where `std::streamoff` is 64 bit; the
  validation script therefore only writes s3t colorvec files.

Legacy build context still matters:

- `README` documents the traditional out-of-tree Autotools flow with `../configure --with-qdp=...` and `make`.
- If you touch source lists or install-facing build behavior, keep both the Autotools and CMake paths in sync.

## Editing Rules For This Repo

- When adding or removing library sources, update both `lib/Makefile.am` and `lib/CMakeLists.txt`.
- When adding or removing test executables, update both `mainprogs/tests/Makefile.am` and `mainprogs/tests/CMakeLists.txt`.
- Factory-registered features should be wired into the relevant aggregate registration point and umbrella include header, not only added as standalone files.
- XML-facing features should usually ship with both a focused executable test and a smoke or regression XML input under `tests/`.
- For new XML-driven local test workflows, create a checkout-local `cfgs/` subdirectory before running and route generated configs, restart XML, split sidecars, command-line `-o` XML outputs, and similar run artifacts there instead of the repository root so those outputs stay easy to ignore.
- If the intended behavior is still being reasoned through, capture it in `specs/hier/*.md` rather than leaving it only in chat history.
- Prefer changing real source under `lib/`, `mainprogs/`, `tests/`, `docs/`, and `specs/`. Avoid editing `build/deps/src/qdpxx` unless you are intentionally changing the local dependency workflow.

## Working Tree Hygiene

- Check `git status` before staging. This checkout can accumulate local artifacts that are not fully ignored.
- Common local-only artifacts seen during current work include `build/`, `XMLDAT`, `.DS_Store`, and `.vscode/`.
- New local XML-driven tests should prefer a checkout-local `cfgs/` directory for generated configs, restart files, and command-line XML outputs rather than writing them into the repository root. Create that directory before launching runs that expect it.
- For the two-level `0-+` workflow, prefer nested output directories such as
  `cfgs/two_level_0mp/outer_*/` or `cfgs/two_level_0mp_1h/outer_*/` so parent
  saves, split sidecars, child HMC logs, reducer summaries, and any generated
  workflow XML stay grouped by level-0 sample or run preset.
- For the two-level `0++` workflow, prefer nested output directories such as
  `cfgs/two_level_0pp_40x10/outer_*/` so parent saves, split sidecars, child
  HMC logs, post-HMC measurement summaries, reducer summaries, and generated
  workflow XML stay grouped by level-0 sample.
- `tests/two_level_glueball/gluecor_measure_launch.sh` writes CSV outputs under
  `gluecor_output/` inside the input directory unless `OUTPUT_ROOT` is
  overridden.
- `tests/two_level_glueball/gluecor_measure_launch.sh` now defaults to the
  `/private/tmp/chroma-build-codex-gluecor/chroma-localbinarydb` build tree,
  falls back across common local build variants, and can still be overridden
  with `CHROMA_BUILD` or `MEASURE_BIN`.
- Additional local smoke artifacts now seen in the repo root include `child_temporal_zone_smoke_cfg_*.lime` and `child_temporal_zone_smoke_restart_*.xml` when running the child HMC smoke input from the checkout root.
- Keep generated outputs and temporary workspace files out of commits unless the change is explicitly about the bootstrap or test workflow itself.

## How To Update This File

When work continues, keep these sections current:

- `Current Workstream`: what feature or subsystem is active right now
- `Build And Test`: the latest known-good commands and fragile assumptions
- `Editing Rules For This Repo`: any newly discovered wiring requirements
- `Working Tree Hygiene`: new generated artifacts or ignore-list gaps worth remembering

If a future effort replaces the local QDPXX recipe, update this file and `docs/build_local_qdpxx.md` together.

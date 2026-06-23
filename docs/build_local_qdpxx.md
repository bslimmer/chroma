# Local QDPXX Build Recipe

This is the tested local recipe for building the current Chroma tree with CMake, wiring it to a compatible QDPXX install, and running the temporal-zone gauge boundary-condition tests.

## Why this recipe exists

- Chroma `3.43+` expects `QDP++ 1.44+`, but this tree also expects the older filedb API shape used by `origin/eloy/localbinarydb`.
- The matching QDPXX branch provides `LocalBinaryBufferWriter` and the `Serializable::serialID() const` signature that this Chroma tree expects.
- QDPXX recursive configure is fragile when the source path contains spaces, so this recipe builds through a no-space symlink at `/private/tmp/chroma-ws`.
- Current Apple clang needs `<array>` included in QDPXX `include/qdp_map_obj_disk.h`; the helper script patches that into the temporary QDPXX worktree before building.
- This repo already includes the current Chroma-side compatibility fixes needed by this recipe, including [lib/actions/ferm/linop/eoprec_clover_orbifold_linop_w.cc](/Users/benslimmer/Documents/Physics%20Research/Frontier%20Work/DD/chroma-dd/chroma/lib/actions/ferm/linop/eoprec_clover_orbifold_linop_w.cc:171) and [lib/util/ferm/key_val_db.h](/Users/benslimmer/Documents/Physics%20Research/Frontier%20Work/DD/chroma-dd/chroma/lib/util/ferm/key_val_db.h:36), so the remaining setup is in the build recipe.

## Prerequisites

- `cmake`
- `pkg-config` or Homebrew `pkgconf`
- `git`
- `make`
- `autoreconf` from `autoconf`/`automake`/`libtool`

On Homebrew that usually means:

```bash
brew install cmake pkgconf autoconf automake libtool
```

## One command

From the repo root:

```bash
./scripts/bootstrap_local_qdpxx_build.sh all
```

That will:

1. reuse `build/deps/src/qdpxx` if it already exists, otherwise clone QDPXX
2. fetch `origin/eloy/localbinarydb`
3. create or reuse a dedicated QDPXX worktree at `/private/tmp/qdpxx-localbinarydb`
4. build and install QDPXX under `/private/tmp/chroma-prefix/qdpxx-localbinarydb`
5. write a small `QDPXXConfig.cmake` wrapper so Chroma can consume that install with `find_package(QDPXX)`
6. configure Chroma in `/private/tmp/chroma-build/chroma-localbinarydb`
7. build `t_temporal_zone_gaugebc`, `t_leapfrog`, `t_hmc_momentum_bc_autodiscovery`, `t_gauge_subdomain_split`, `t_gauge_subdomain_gauge_hmc_validation`, plus `hmc` and `gauge_subdomain_split`
8. run the focused executable tests

## Script modes

```bash
./scripts/bootstrap_local_qdpxx_build.sh bootstrap
./scripts/bootstrap_local_qdpxx_build.sh build-tests
./scripts/bootstrap_local_qdpxx_build.sh run-tests
```

- `bootstrap` stops after QDPXX install and Chroma configure.
- `build-tests` builds the focused boundary-condition, subdomain-validation, and smoke-run targets.
- `run-tests` runs `t_temporal_zone_gaugebc`, `t_hmc_momentum_bc_autodiscovery`, and the `t_leapfrog` temporal-zone smoke input after building them.

## Useful overrides

The script is parameterized with environment variables if you want different locations:

```bash
JOBS=8 \
BUILD_ROOT=/private/tmp/chroma-build \
PREFIX_ROOT=/private/tmp/chroma-prefix \
CHROMA_ALIAS=/private/tmp/chroma-ws \
QDPXX_GIT_DIR=/path/to/qdpxx-clone \
QDPXX_SKIP_FETCH=1 \
QDPXX_WORKTREE=/private/tmp/qdpxx-localbinarydb \
QDPXX_BUILD=/private/tmp/chroma-build/qdpxx-localbinarydb \
QDPXX_PREFIX=/private/tmp/chroma-prefix/qdpxx-localbinarydb \
CHROMA_BUILD=/private/tmp/chroma-build/chroma-localbinarydb \
./scripts/bootstrap_local_qdpxx_build.sh run-tests
```

If `QDPXX_GIT_DIR` does not exist, the script clones from `QDPXX_REMOTE_URL`, which defaults to `https://github.com/usqcd-software/qdpxx.git`.

If you already have `origin/eloy/localbinarydb` in the local QDPXX clone and do not want to depend on a fresh network fetch, set:

```bash
QDPXX_SKIP_FETCH=1
```

The script will then reuse the locally available ref and still verify that it exists before building.

If recursive QDPXX submodule fetches are also unavailable, the script will fall back to copying the already-populated `other_libs/{filedb,libintrin,qio,xpath_reader}` trees from the local QDPXX source checkout into the temporary worktree.

## Manual commands

If you would rather run the recipe by hand, these are the important steps:

```bash
git -C build/deps/src/qdpxx fetch origin eloy/localbinarydb
git -C build/deps/src/qdpxx worktree add -B localbinarydb /private/tmp/qdpxx-localbinarydb origin/eloy/localbinarydb
git -C /private/tmp/qdpxx-localbinarydb submodule update --init --recursive
autoreconf -fi
```

If submodule fetches are unavailable, copy the populated `other_libs/` trees from `build/deps/src/qdpxx` into `/private/tmp/qdpxx-localbinarydb` before running `autoreconf -fi`.

Patch the QDPXX worktree once for Apple clang by adding `#include <array>` to `include/qdp_map_obj_disk.h`, then build:

```bash
mkdir -p /private/tmp/chroma-build/qdpxx-localbinarydb
cd /private/tmp/chroma-build/qdpxx-localbinarydb
CXXFLAGS='-std=c++11' /private/tmp/qdpxx-localbinarydb/configure \
  --prefix=/private/tmp/chroma-prefix/qdpxx-localbinarydb \
  --enable-parallel-arch=scalar
make -j8
make install
```

Write the wrapper package under `/private/tmp/chroma-prefix/qdpxx-localbinarydb/lib/cmake/QDPXX/`, then configure Chroma:

```bash
ln -sfn "$PWD" /private/tmp/chroma-ws
cmake -S /private/tmp/chroma-ws -B /private/tmp/chroma-build/chroma-localbinarydb \
  -DQDPXX_DIR=/private/tmp/chroma-prefix/qdpxx-localbinarydb/lib/cmake/QDPXX \
  -DCMAKE_CXX_STANDARD=11 \
  -DCMAKE_CXX_EXTENSIONS=OFF
```

If `other_libs/qdp-lapack` is already populated in the checkout, no extra submodule command is needed before configuring Chroma.

Build the validation-related executables:

```bash
cmake --build /private/tmp/chroma-build/chroma-localbinarydb \
  --target \
    t_temporal_zone_gaugebc \
    t_leapfrog \
    t_hmc_momentum_bc_autodiscovery \
    t_gauge_subdomain_split \
    t_gauge_subdomain_gauge_hmc_validation \
    hmc \
    gauge_subdomain_split \
  -j8
```

Run the focused executable tests:

```bash
/private/tmp/chroma-build/chroma-localbinarydb/mainprogs/tests/t_temporal_zone_gaugebc

/private/tmp/chroma-build/chroma-localbinarydb/mainprogs/tests/t_hmc_momentum_bc_autodiscovery

/private/tmp/chroma-build/chroma-localbinarydb/mainprogs/tests/t_leapfrog \
  -i "$PWD/tests/t_leapfrog/t_leapfrog.temporal_zone_gaugebc.ini.xml" \
  -o /private/tmp/chroma-build/chroma-localbinarydb/mainprogs/tests/t_leapfrog.temporal_zone_gaugebc.out.xml \
  -l /private/tmp/chroma-build/chroma-localbinarydb/mainprogs/tests/t_leapfrog.temporal_zone_gaugebc.log.xml
```

## Expected results

- `t_temporal_zone_gaugebc` should exit `0`.
- `t_hmc_momentum_bc_autodiscovery` should exit `0`.
- `t_gauge_subdomain_split` should exit `0`.
- `t_gauge_subdomain_gauge_hmc_validation` should build successfully and pass after running the split plus paired child-HMC workflow.
- `t_leapfrog` should parse `TEMPORAL_ZONE_GAUGEBC`, complete the trajectory, and write XML outputs under the Chroma build tree.
- `hmc` should complete `tests/gauge_subdomain_split/hmc_temporal_zone_child_smoke.ini.xml`.
- `gauge_subdomain_split` should complete the validation split input and produce child configs plus a sidecar XML.

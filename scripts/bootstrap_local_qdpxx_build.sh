#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_ROOT="${BUILD_ROOT:-/private/tmp/chroma-build}"
PREFIX_ROOT="${PREFIX_ROOT:-/private/tmp/chroma-prefix}"
CHROMA_ALIAS="${CHROMA_ALIAS:-/private/tmp/chroma-ws}"
QDPXX_REMOTE_URL="${QDPXX_REMOTE_URL:-https://github.com/usqcd-software/qdpxx.git}"
QDPXX_REF="${QDPXX_REF:-origin/eloy/localbinarydb}"
QDPXX_BRANCH_NAME="${QDPXX_BRANCH_NAME:-localbinarydb}"
QDPXX_SKIP_FETCH="${QDPXX_SKIP_FETCH:-0}"
QDPXX_WORKTREE="${QDPXX_WORKTREE:-/private/tmp/qdpxx-localbinarydb}"
QDPXX_BUILD="${QDPXX_BUILD:-${BUILD_ROOT}/qdpxx-localbinarydb}"
QDPXX_PREFIX="${QDPXX_PREFIX:-${PREFIX_ROOT}/qdpxx-localbinarydb}"
CHROMA_BUILD="${CHROMA_BUILD:-${BUILD_ROOT}/chroma-localbinarydb}"
QDPXX_CONFIG_DIR="${QDPXX_PREFIX}/lib/cmake/QDPXX"
LEAPFROG_INPUT="${REPO_ROOT}/tests/t_leapfrog/t_leapfrog.temporal_zone_gaugebc.ini.xml"
LEAPFROG_OUT="${CHROMA_BUILD}/mainprogs/tests/t_leapfrog.temporal_zone_gaugebc.out.xml"
LEAPFROG_LOG="${CHROMA_BUILD}/mainprogs/tests/t_leapfrog.temporal_zone_gaugebc.log.xml"

if [ -d "${REPO_ROOT}/build/deps/src/qdpxx/.git" ]; then
  QDPXX_GIT_DIR="${QDPXX_GIT_DIR:-${REPO_ROOT}/build/deps/src/qdpxx}"
else
  QDPXX_GIT_DIR="${QDPXX_GIT_DIR:-${BUILD_ROOT}/qdpxx-origin}"
fi

if command -v getconf >/dev/null 2>&1; then
  DEFAULT_JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
elif command -v sysctl >/dev/null 2>&1; then
  DEFAULT_JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
else
  DEFAULT_JOBS=4
fi
JOBS="${JOBS:-${DEFAULT_JOBS}}"

say() {
  printf '\n==> %s\n' "$1"
}

die() {
  printf 'error: %s\n' "$1" >&2
  exit 1
}

require_tool() {
  command -v "$1" >/dev/null 2>&1 || die "missing required tool: $1"
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [bootstrap|build-tests|run-tests|all]

bootstrap   Build and install QDPXX, write the CMake wrapper package, and
            configure Chroma in ${CHROMA_BUILD}.
build-tests bootstrap, then build the temporal-zone and subdomain validation targets.
run-tests   build-tests, then run the focused executable tests.
all         same as run-tests.

Environment overrides:
  JOBS
  BUILD_ROOT
  PREFIX_ROOT
  CHROMA_ALIAS
  QDPXX_GIT_DIR
  QDPXX_REMOTE_URL
  QDPXX_REF
  QDPXX_BRANCH_NAME
  QDPXX_SKIP_FETCH
  QDPXX_WORKTREE
  QDPXX_BUILD
  QDPXX_PREFIX
  CHROMA_BUILD
EOF
}

ensure_prereqs() {
  require_tool git
  require_tool cmake
  require_tool autoreconf
  require_tool make
  require_tool perl
  require_tool pkg-config
}

ensure_qdpxx_source() {
  if [ ! -d "${QDPXX_GIT_DIR}/.git" ]; then
    say "Cloning QDPXX source into ${QDPXX_GIT_DIR}"
    mkdir -p "$(dirname "${QDPXX_GIT_DIR}")"
    git clone "${QDPXX_REMOTE_URL}" "${QDPXX_GIT_DIR}"
  fi

  if [ "${QDPXX_SKIP_FETCH}" = "1" ]; then
    say "Skipping fetch for ${QDPXX_REF}; using locally available refs"
  else
    say "Fetching ${QDPXX_REF}"
    if ! git -C "${QDPXX_GIT_DIR}" fetch origin eloy/localbinarydb; then
      if git -C "${QDPXX_GIT_DIR}" rev-parse --verify --quiet "${QDPXX_REF}" >/dev/null; then
        say "Fetch failed, but ${QDPXX_REF} exists locally; continuing"
      else
        die "failed to fetch ${QDPXX_REF} and no local copy is available"
      fi
    fi
  fi

  git -C "${QDPXX_GIT_DIR}" rev-parse --verify --quiet "${QDPXX_REF}" >/dev/null || \
    die "QDPXX ref not found locally: ${QDPXX_REF}"
}

populate_qdpxx_submodule_from_local_source() {
  local submodule
  local src
  local dst

  submodule="$1"
  src="${QDPXX_GIT_DIR}/other_libs/${submodule}"
  dst="${QDPXX_WORKTREE}/other_libs/${submodule}"

  [ -d "${src}" ] || die "local QDPXX source is missing submodule tree: ${src}"

  say "Populating ${submodule} from local source checkout"
  mkdir -p "${dst}"
  (
    cd "${src}"
    tar --exclude='.git' -cf - .
  ) | (
    cd "${dst}"
    tar -xf -
  )

  [ -f "${dst}/configure.ac" ] || [ -f "${dst}/Makefile.am" ] || \
    die "failed to populate ${submodule} into ${dst}"
}

sync_qdpxx_submodules() {
  say "Syncing QDPXX submodules"
  if git -C "${QDPXX_WORKTREE}" submodule update --init --recursive; then
    return
  fi

  if [ "${QDPXX_SKIP_FETCH}" != "1" ]; then
    die "failed to sync QDPXX submodules"
  fi

  say "Submodule update could not reach remotes; falling back to local source trees"
  populate_qdpxx_submodule_from_local_source filedb
  populate_qdpxx_submodule_from_local_source libintrin
  populate_qdpxx_submodule_from_local_source qio
  populate_qdpxx_submodule_from_local_source xpath_reader
}

qdpxx_worktree_has_only_bootstrap_changes() {
  local status
  local path

  while IFS= read -r status; do
    [ -n "${status}" ] || continue
    path="${status#?? }"
    case "${path}" in
      INSTALL|config/depcomp|config/install-sh|config/missing|config/mkinstalldirs|include/qdp_map_obj_disk.h|configure~)
        ;;
      *)
        return 1
        ;;
    esac
  done < <(git -C "${QDPXX_WORKTREE}" status --porcelain)

  return 0
}

ensure_qdpxx_worktree() {
  local qdpxx_status

  if [ -e "${QDPXX_WORKTREE}" ] && ! git -C "${QDPXX_WORKTREE}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    die "QDPXX_WORKTREE exists but is not a git worktree: ${QDPXX_WORKTREE}"
  fi

  if ! git -C "${QDPXX_WORKTREE}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    say "Creating QDPXX worktree ${QDPXX_WORKTREE}"
    mkdir -p "$(dirname "${QDPXX_WORKTREE}")"
    git -C "${QDPXX_GIT_DIR}" worktree add -B "${QDPXX_BRANCH_NAME}" "${QDPXX_WORKTREE}" "${QDPXX_REF}"
  else
    qdpxx_status="$(git -C "${QDPXX_WORKTREE}" status --porcelain)"
    if [ -n "${qdpxx_status}" ]; then
      if qdpxx_worktree_has_only_bootstrap_changes && \
         [ "$(git -C "${QDPXX_WORKTREE}" rev-parse HEAD)" = "$(git -C "${QDPXX_WORKTREE}" rev-parse "${QDPXX_REF}")" ]; then
        say "Reusing dirty QDPXX worktree with only local bootstrap changes"
      else
        die "refusing to reuse a dirty QDPXX worktree: ${QDPXX_WORKTREE}"
      fi
    else
      say "Updating existing QDPXX worktree ${QDPXX_WORKTREE}"
      if [ "${QDPXX_SKIP_FETCH}" != "1" ]; then
        if ! git -C "${QDPXX_WORKTREE}" fetch origin eloy/localbinarydb; then
          say "Worktree fetch failed; continuing with locally available refs"
        fi
      fi
      git -C "${QDPXX_WORKTREE}" switch "${QDPXX_BRANCH_NAME}" >/dev/null 2>&1 || \
        git -C "${QDPXX_WORKTREE}" switch -c "${QDPXX_BRANCH_NAME}" "${QDPXX_REF}"
      git -C "${QDPXX_WORKTREE}" merge --ff-only "${QDPXX_REF}"
    fi
  fi

  sync_qdpxx_submodules
}

patch_qdpxx_for_clang() {
  local file
  file="${QDPXX_WORKTREE}/include/qdp_map_obj_disk.h"

  if ! grep -Fq '#include <array>' "${file}"; then
    say "Patching ${file} for Apple clang"
    if grep -Fq '#include <vector>' "${file}"; then
      perl -0pi -e 's/#include <vector>\n/#include <vector>\n#include <array>\n/' "${file}"
    elif grep -Fq '#include <limits>' "${file}"; then
      perl -0pi -e 's/#include <limits>\n/#include <array>\n#include <limits>\n/' "${file}"
    else
      die "could not find insertion point in ${file}"
    fi
  fi
}

build_qdpxx() {
  say "Bootstrapping QDPXX"
  (
    cd "${QDPXX_WORKTREE}"
    autoreconf -fi
  )

  say "Configuring QDPXX in ${QDPXX_BUILD}"
  mkdir -p "${QDPXX_BUILD}" "${QDPXX_PREFIX}"
  (
    cd "${QDPXX_BUILD}"
    CXXFLAGS='-std=c++11' "${QDPXX_WORKTREE}/configure" \
      --prefix="${QDPXX_PREFIX}" \
      --enable-parallel-arch=scalar
    make -j"${JOBS}"
    make install
  )
}

write_qdpxx_cmake_wrapper() {
  say "Writing QDPXX CMake package wrapper into ${QDPXX_CONFIG_DIR}"
  mkdir -p "${QDPXX_CONFIG_DIR}"

  cat > "${QDPXX_CONFIG_DIR}/QDPXXConfig.cmake" <<'EOF'
if(TARGET QDPXX::qdp)
  return()
endif()

if(NOT TARGET qdp)
  include("${CMAKE_CURRENT_LIST_DIR}/../../../share/FindQDPXX.cmake")
endif()

if(TARGET qdp AND NOT TARGET QDPXX::qdp)
  add_library(QDPXX::qdp INTERFACE IMPORTED)
  target_link_libraries(QDPXX::qdp INTERFACE qdp)
endif()

if(NOT TARGET QDPXX::qdp)
  message(FATAL_ERROR "QDPXXConfig.cmake could not create QDPXX::qdp")
endif()
EOF

  cat > "${QDPXX_CONFIG_DIR}/QDPXXConfigVersion.cmake" <<'EOF'
set(PACKAGE_VERSION "1")
set(PACKAGE_VERSION_COMPATIBLE TRUE)
set(PACKAGE_VERSION_EXACT TRUE)
EOF
}

configure_chroma() {
  say "Creating no-space source alias ${CHROMA_ALIAS}"
  mkdir -p "$(dirname "${CHROMA_ALIAS}")"
  ln -sfn "${REPO_ROOT}" "${CHROMA_ALIAS}"

  if [ -f "${REPO_ROOT}/other_libs/qdp-lapack/CMakeLists.txt" ] && \
     [ -d "${REPO_ROOT}/other_libs/qdp-lapack/include" ] && \
     [ -d "${REPO_ROOT}/other_libs/qdp-lapack/lib" ]; then
    say "Using existing Chroma qdp-lapack checkout"
  else
    say "Initializing required Chroma submodule"
    git -C "${REPO_ROOT}" submodule update --init other_libs/qdp-lapack
  fi

  say "Configuring Chroma in ${CHROMA_BUILD}"
  cmake -S "${CHROMA_ALIAS}" -B "${CHROMA_BUILD}" \
    -DQDPXX_DIR="${QDPXX_CONFIG_DIR}" \
    -DCMAKE_CXX_STANDARD=11 \
    -DCMAKE_CXX_EXTENSIONS=OFF
}

bootstrap() {
  ensure_prereqs
  ensure_qdpxx_source
  ensure_qdpxx_worktree
  patch_qdpxx_for_clang
  build_qdpxx
  write_qdpxx_cmake_wrapper
  configure_chroma
}

build_tests() {
  bootstrap

  say "Building test targets"
  cmake --build "${CHROMA_BUILD}" \
    --target \
      t_temporal_zone_gaugebc \
      t_leapfrog \
      t_hmc_momentum_bc_autodiscovery \
      t_gauge_subdomain_split \
      t_gauge_subdomain_gauge_hmc_validation \
      gluecor_measure \
      hmc \
      gauge_subdomain_split \
    -j"${JOBS}"
}

run_tests() {
  build_tests

  say "Running t_temporal_zone_gaugebc"
  "${CHROMA_BUILD}/mainprogs/tests/t_temporal_zone_gaugebc"

  say "Running t_hmc_momentum_bc_autodiscovery"
  "${CHROMA_BUILD}/mainprogs/tests/t_hmc_momentum_bc_autodiscovery"

  say "Running t_leapfrog temporal-zone smoke test"
  "${CHROMA_BUILD}/mainprogs/tests/t_leapfrog" \
    -i "${LEAPFROG_INPUT}" \
    -o "${LEAPFROG_OUT}" \
    -l "${LEAPFROG_LOG}"

  cat <<EOF

Outputs:
  t_temporal_zone_gaugebc: ${CHROMA_BUILD}/mainprogs/tests/t_temporal_zone_gaugebc
  t_hmc_momentum_bc_autodiscovery: ${CHROMA_BUILD}/mainprogs/tests/t_hmc_momentum_bc_autodiscovery
  t_leapfrog log xml:      ${LEAPFROG_LOG}
  t_leapfrog output xml:   ${LEAPFROG_OUT}
  hmc executable:          ${CHROMA_BUILD}/mainprogs/main/hmc
  gluecor_measure:         ${CHROMA_BUILD}/mainprogs/main/gluecor_measure
  gauge_subdomain_split:   ${CHROMA_BUILD}/mainprogs/main/gauge_subdomain_split
  t_gauge_subdomain_gauge_hmc_validation:
                           ${CHROMA_BUILD}/mainprogs/tests/t_gauge_subdomain_gauge_hmc_validation
EOF
}

MODE="${1:-all}"

case "${MODE}" in
  bootstrap)
    bootstrap
    ;;
  build-tests)
    build_tests
    ;;
  run-tests|all)
    run_tests
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage
    exit 1
    ;;
esac

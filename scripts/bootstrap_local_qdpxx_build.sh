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
build-tests bootstrap, then build t_temporal_zone_gaugebc,
            t_gauge_subdomain_split, and t_leapfrog.
run-tests   build-tests, then run all three tests.
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

  say "Fetching ${QDPXX_REF}"
  git -C "${QDPXX_GIT_DIR}" fetch origin eloy/localbinarydb
}

ensure_qdpxx_worktree() {
  if [ -e "${QDPXX_WORKTREE}" ] && [ ! -d "${QDPXX_WORKTREE}/.git" ]; then
    die "QDPXX_WORKTREE exists but is not a git worktree: ${QDPXX_WORKTREE}"
  fi

  if [ ! -d "${QDPXX_WORKTREE}/.git" ]; then
    say "Creating QDPXX worktree ${QDPXX_WORKTREE}"
    mkdir -p "$(dirname "${QDPXX_WORKTREE}")"
    git -C "${QDPXX_GIT_DIR}" worktree add -B "${QDPXX_BRANCH_NAME}" "${QDPXX_WORKTREE}" "${QDPXX_REF}"
  else
    if [ -n "$(git -C "${QDPXX_WORKTREE}" status --porcelain)" ]; then
      die "refusing to reuse a dirty QDPXX worktree: ${QDPXX_WORKTREE}"
    fi

    say "Updating existing QDPXX worktree ${QDPXX_WORKTREE}"
    git -C "${QDPXX_WORKTREE}" fetch origin eloy/localbinarydb
    git -C "${QDPXX_WORKTREE}" switch "${QDPXX_BRANCH_NAME}" >/dev/null 2>&1 || \
      git -C "${QDPXX_WORKTREE}" switch -c "${QDPXX_BRANCH_NAME}" "${QDPXX_REF}"
    git -C "${QDPXX_WORKTREE}" merge --ff-only "${QDPXX_REF}"
  fi

  say "Syncing QDPXX submodules"
  git -C "${QDPXX_WORKTREE}" submodule update --init --recursive
}

patch_qdpxx_for_clang() {
  local file
  file="${QDPXX_WORKTREE}/include/qdp_map_obj_disk.h"

  if ! grep -Fq '#include <array>' "${file}"; then
    grep -Fq '#include <vector>' "${file}" || die "could not find insertion point in ${file}"
    say "Patching ${file} for Apple clang"
    perl -0pi -e 's/#include <vector>\n/#include <vector>\n#include <array>\n/' "${file}"
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
include("${CMAKE_CURRENT_LIST_DIR}/../../../share/FindQDPXX.cmake")

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

  say "Initializing required Chroma submodule"
  git -C "${REPO_ROOT}" submodule update --init other_libs/qdp-lapack

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
    --target t_temporal_zone_gaugebc t_gauge_subdomain_split t_leapfrog \
    -j"${JOBS}"
}

run_tests() {
  build_tests

  say "Running t_temporal_zone_gaugebc"
  "${CHROMA_BUILD}/mainprogs/tests/t_temporal_zone_gaugebc"

  say "Running t_gauge_subdomain_split"
  "${CHROMA_BUILD}/mainprogs/tests/t_gauge_subdomain_split"

  say "Running t_leapfrog temporal-zone smoke test"
  "${CHROMA_BUILD}/mainprogs/tests/t_leapfrog" \
    -i "${LEAPFROG_INPUT}" \
    -o "${LEAPFROG_OUT}" \
    -l "${LEAPFROG_LOG}"

  cat <<EOF

Outputs:
  t_temporal_zone_gaugebc: ${CHROMA_BUILD}/mainprogs/tests/t_temporal_zone_gaugebc
  t_gauge_subdomain_split: ${CHROMA_BUILD}/mainprogs/tests/t_gauge_subdomain_split
  t_leapfrog log xml:      ${LEAPFROG_LOG}
  t_leapfrog output xml:   ${LEAPFROG_OUT}
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

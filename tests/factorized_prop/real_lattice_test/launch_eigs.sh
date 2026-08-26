#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------

CHROMA_BUILD="/qcd/work/JLabLQCD/slimmer/chromaform-rocm/build/chroma-quda-qdp-jit-double-nd4-cmake-superbblas-hip"

CHROMA_BIN="${CHROMA_BIN:-$CHROMA_BUILD/mainprogs/main/chroma}"

#QPhix settings
QPHIX_BY="${QPHIX_BY:-4}"
QPHIX_BZ="${QPHIX_BZ:-2}"
QPHIX_PXY="${QPHIX_PXY:-1}"
QPHIX_PXYZ="${QPHIX_PXYZ:-1}"
QPHIX_MINCT="${QPHIX_MINCT:-1}"
QPHIX_NCORES="${QPHIX_NCORES:-2}"
QPHIX_SY="${QPHIX_SY:-1}"
QPHIX_SZ="${QPHIX_SZ:-1}"

mpirun -np 2 \
	"$CHROMA_BIN" -i eigs.ini.xml -o eigs.out.xml \
	-geom 1 1 1 2 \
    	-by "$QPHIX_BY" \
    	-bz "$QPHIX_BZ" \
    	-pxy "$QPHIX_PXY" \
    	-pxyz "$QPHIX_PXYZ" \
    	-minct "$QPHIX_MINCT" \
    	-c "$QPHIX_NCORES" \
    	-sy "$QPHIX_SY" \
    	-sz "$QPHIX_SZ"    

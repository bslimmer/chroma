#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------

CHROMA_BUILD="${CHROMA_BUILD:-/u/home/slimmer/mychroma_test/chromaform/build/chroma-mgproto-qphix-qdpxx-double-nd4-avx2-superbblas-cpu}"

CHROMA_BIN="${CHROMA_BIN:-$CHROMA_BUILD/mainprogs/main/chroma}"

#QPhix settings
QPHIX_BY="${QPHIX_BY:-1}"
QPHIX_BZ="${QPHIX_BZ:-1}"
QPHIX_PXY="${QPHIX_PXY:-1}"
QPHIX_PXYZ="${QPHIX_PXYZ:-1}"
QPHIX_MINCT="${QPHIX_MINCT:-1}"
QPHIX_NCORES="${QPHIX_NCORES:-1}"
QPHIX_SY="${QPHIX_SY:-1}"
QPHIX_SZ="${QPHIX_SZ:-1}"

mpirun -np 1 \
	-x OMP_NUM_THREADS=1 \
	"$CHROMA_BIN" -i test_prop_and_matelem_distillation_superb.ini.xml -o test_factorized_prop.out.xml \
    	-by "$QPHIX_BY" \
    	-bz "$QPHIX_BZ" \
    	-pxy "$QPHIX_PXY" \
    	-pxyz "$QPHIX_PXYZ" \
    	-minct "$QPHIX_MINCT" \
    	-c "$QPHIX_NCORES" \
    	-sy "$QPHIX_SY" \
    	-sz "$QPHIX_SZ"    

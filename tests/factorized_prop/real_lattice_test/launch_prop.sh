#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# User settings
# --------------------------------------------------------------------

CHROMA_BUILD="/qcd/work/JLabLQCD/slimmer/chromaform-rocm/build/chroma-quda-qdp-jit-double-nd4-cmake-superbblas-hip"
#CHROMA_BUILD="${CHROMA_BUILD:-/qcd/work/JLabLQCD/slimmer/chromaform/build/chroma-mgproto-qphix-qdpxx-double-nd4-avx512-superbblas-cpu}"

CHROMA_BIN="${CHROMA_BIN:-$CHROMA_BUILD/mainprogs/main/chroma}"

#Trying to get multiple GPUs to work
#export SB_MPI_GPU=1
#export SB_CACHEGB_GPU=8   
#export QUDA_ENABLE_P2P=0
#export QUDA_ENABLE_GDR=0
#export QUDA_ENABLE_NVSHMEM=0
#export QUDA_ENABLE_MPS=0
#export OMP_NUM_THREADS=1

#export ROCR_VISIBLE_DEVICES=0

#QPhix settings
QPHIX_BY="${QPHIX_BY:-4}"
QPHIX_BZ="${QPHIX_BZ:-2}"
QPHIX_PXY="${QPHIX_PXY:-1}"
QPHIX_PXYZ="${QPHIX_PXYZ:-1}"
QPHIX_MINCT="${QPHIX_MINCT:-1}"
QPHIX_NCORES="${QPHIX_NCORES:-4}"
QPHIX_SY="${QPHIX_SY:-1}"
QPHIX_SZ="${QPHIX_SZ:-1}"

mpirun -np 4 \
	gdb --batch -ex 'b QDP_abort' -ex r -ex bt --args "$CHROMA_BIN" -i test_prop_and_matelem_distillation_superb_exact.ini.xml -o exact_prop.out.xml \
	-geom 1 1 2 2 \
    	-by "$QPHIX_BY" \
    	-bz "$QPHIX_BZ" \
    	-pxy "$QPHIX_PXY" \
    	-pxyz "$QPHIX_PXYZ" \
    	-minct "$QPHIX_MINCT" \
    	-c "$QPHIX_NCORES" \
    	-sy "$QPHIX_SY" \
    	-sz "$QPHIX_SZ" 	

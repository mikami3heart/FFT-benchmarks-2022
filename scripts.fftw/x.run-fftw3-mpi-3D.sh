#! /bin/bash
#PJM -N FFTW-3D-MPI-MEASURE
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=04:00:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=48"
#PJM -j
#PJM -S

set -x
hostname
date
WD=${HOME}/tmp/check_fftw_3d_mpi_measure
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*
rm -rf ${WD}/output.*

SRC=${HOME}/fftw/src_fftw3

FFTW_DIR=${HOME}/fftw/github-fujitsu/local

CFLAGS="-Kfast,openmp -I${FFTW_DIR}/include "
FFLAGS="-Kfast,openmp -I${FFTW_DIR}/include "
LDFLAGS="-Kfast,openmp --linkfortran "
LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3_mpi -lfftw3 "
cp $SRC/fftw3_complex_3D.mpi.c main.c
mpifcc  -o fftw.ex $CFLAGS  main.c  $LDFLAGS

export LD_LIBRARY_PATH=${FFTW_DIR}/lib:${LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1

#	for i in 4 8 12 24 48
for i in 1 4 16 48
do
NPROCS=${i}
mpiexec -n ${NPROCS} ./fftw.ex
#	mpiexec -n ${NPROCS} --std stdout.1node-${NPROCS}procs.txt ./fftw.ex
#	cat stdout.fft3d_mpi.1node-${NPROCS}procs.txt
done

more output.*/0/*/*

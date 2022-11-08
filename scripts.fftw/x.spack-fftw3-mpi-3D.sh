#! /bin/bash
#PJM -N SPACK-FFTW3-3D-MPI-MEASURE
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=02:00:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=48"
#PJM -j
#PJM -S

module list
source /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load fujitsu-fftw/ytok4j4

export LANG=C
echo job manager reserved the following resources
echo PJM_NODE         : ${PJM_NODE}
echo PJM_MPI_PROC     : ${PJM_MPI_PROC}
echo PJM_PROC_BY_NODE : ${PJM_PROC_BY_NODE}

set -x
hostname
date
WD=${HOME}/tmp/spack_fftw3_3d_mpi
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

SRC=${HOME}/fftw/src
#	FFTW_DIR=${HOME}/fftw/fftw-3.3.10/root
FFTW_DIR=/vol0004/apps/oss/spack-v0.16.2/opt/spack/linux-rhel8-a64fx/fj-4.6.1/fujitsu-fftw-master-ytok4j4l2qa355nopq524zx6ecktzjg2

CFLAGS="-Kfast -I${FFTW_DIR}/include "
FFLAGS="-Kfast -I${FFTW_DIR}/include -Nlst=t "
LDFLAGS="-Kfast --linkfortran "
LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3_mpi -lfftw3 "
cp $SRC/fftw3_complex_3D.mpi.c main.c
mpifcc  -o fftw.ex $CFLAGS  main.c  $LDFLAGS


for i in 4 8 12 24 48
do
export OMP_NUM_THREADS=1
NPROCS=${i}
mpiexec -n ${NPROCS} --std stdout.fft3d_mpi.${NPROCS}procs.txt ./fftw.ex
cat stdout.fft3d_mpi.${NPROCS}procs.txt
done


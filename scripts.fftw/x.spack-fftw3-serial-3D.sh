#! /bin/bash
#PJM -N SPACK-FFTW3-3D-SERIAL
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:30:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -S

module list
source /vol0004/apps/oss/spack/share/spack/setup-env.sh
#   spack find -lx | grep fftw
#   spack load /ytok4j4
spack load fujitsu-fftw/ytok4j4

set -x
hostname
date
WD=${HOME}/tmp/spack_fftw3_3d_serial
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

SRC=${HOME}/fftw/src
#	FFTW_DIR=${HOME}/fftw/fftw-3.3.10/root
FFTW_DIR=/vol0004/apps/oss/spack-v0.16.2/opt/spack/linux-rhel8-a64fx/fj-4.6.1/fujitsu-fftw-master-ytok4j4l2qa355nopq524zx6ecktzjg2

CFLAGS="-Kfast -I${FFTW_DIR}/include "
FFLAGS="-Kfast -I${FFTW_DIR}/include -Nlst=t "
LDFLAGS="-Kfast "
LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 "
cp $SRC/main_fftw3_complex_3d.serial.F main.f90
mpifrt  -o fftw.ex $FFLAGS  main.f90  $LDFLAGS

#	for i in 4 8 12 24 48
for i in 1
do
export OMP_NUM_THREADS=1
NPROCS=${i}
time mpiexec -n ${NPROCS} --std stdout.fft3d_mpi.1node-${NPROCS}procs.txt ./fftw.ex
cat stdout.fft3d_mpi.1node-${NPROCS}procs.txt
done


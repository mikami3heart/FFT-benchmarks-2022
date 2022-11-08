#! /bin/bash
#PJM -N SPACK-FFTW3-1D-SERIAL
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:30:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -S

module list
#	ENV_COPY=${C_INCLUDE_PATH}

source /vol0004/apps/oss/spack/share/spack/setup-env.sh
#   spack find -lx | grep fftw
#   spack find -lv fftw
#   spack load /ytok4j4
#	spack load fujitsu-fftw/ytok4j4
#	FFTW_DIR=/vol0004/apps/oss/spack-v0.16.2/opt/spack/linux-rhel8-a64fx/fj-4.6.1/fujitsu-fftw-master-ytok4j4l2qa355nopq524zx6ecktzjg2
#	FFTW_DIR=/vol0004/apps/oss/spack-v0.16.2/opt/spack/linux-rhel8-a64fx/fj-4.6.1/fujitsu-fftw-master-ytok4j4l2qa355nopq524zx6ecktzjg2

spack load fujitsu-fftw%fj@4.8.0

set -x
hostname
date
WD=${HOME}/tmp/spack_fftw3_1d_serial
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

SRC=${HOME}/fftw/src_fftw3

FFTW_INCLUDE=`echo "${C_INCLUDE_PATH}"| cut -f1 -d':'`
FFTW_DIR=`dirname "${FFTW_INCLUDE}"`

CFLAGS="-Kfast -I${FFTW_DIR}/include "
FFLAGS="-Kfast -I${FFTW_DIR}/include -Nlst=t "
LDFLAGS="-L${FFTW_DIR}/lib -lfftw3 "
cp $SRC/main_fftw3_complex_1d.serial.F main.F90
#	mpifrtpx  -o fftw.ex $FFLAGS  main.F90  $LDFLAGS
#	exit
mpifrt  -o fftw.ex $FFLAGS  main.F90  $LDFLAGS
file fftw.ex
ldd fftw.ex

export OMP_NUM_THREADS=1
./fftw.ex

echo Also run using mpiexec
rm stdout.*
export OMP_NUM_THREADS=1
mpiexec -n 1 ./fftw.ex
cat stdout.*
exit

#	for i in 4 8 12 24 48
for i in 1
do
export OMP_NUM_THREADS=1
NPROCS=${i}
time mpiexec -n ${NPROCS} --std stdout.fft3d_mpi.1node-${NPROCS}procs.txt ./fftw.ex
cat stdout.fft3d_mpi.1node-${NPROCS}procs.txt
done


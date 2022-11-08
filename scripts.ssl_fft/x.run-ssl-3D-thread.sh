#! /bin/bash
#PJM -N SSL2FFT-3D-THREAD
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:30:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -S

module list

set -x
hostname
date
WD=${HOME}/tmp/check_ssl2_dm_v3dcft2
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

SRC=${HOME}/ssl_fft/src_ssl2fft

CFLAGS="-Kfast -I${FFTW_DIR}/include "
FFLAGS="-Kfast,openmp -I${FFTW_DIR}/include -Nlst=t -w "
LDFLAGS="-Kfast,openmp "
LDFLAGS+="-SSL2BLAMP --linkstl=libfjc++ "
cp $SRC/main_ssl2fft_complex_3D.thread.F main.F90
#	frtpx  -o ssl2fft.ex $FFLAGS  main.F90  $LDFLAGS
#	exit
frt  -o ssl2fft.ex $FFLAGS  main.F90  $LDFLAGS

#	for i in 1
for i in 1 2 4 6 8 10 12
do
export OMP_NUM_THREADS=${i}
./ssl2fft.ex
done
exit


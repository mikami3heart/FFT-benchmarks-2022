#! /bin/bash
#	#PJM -N FFTW3-3D-many1D
#PJM -N FFTW3-3D-many1D-THREADS
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:60:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -S

module list

set -x
hostname
date
WD=${HOME}/tmp/check_fftw3_3d_many1d_threads
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

#	FFTW_DIR=${HOME}/fftw/fftw-3.3.10/root
FFTW_DIR=${HOME}/fftw/github-fujitsu/local

CFLAGS="-Kfast,openmp -I${FFTW_DIR}/include "
FFLAGS="-Kfast,openmp -I${FFTW_DIR}/include -Nlst=t "
LDFLAGS="-Kfast "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 --linkstl=libfjc++ "
LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3_omp -lfftw3 --linkstl=libfjc++ "

SRC=${HOME}/fftw/src_fftw3
cp $SRC/main_fftw3_complex_3D_many1d.F main.F90
#	mpifrtpx  -o fftw.ex $FFLAGS  main.F90  $LDFLAGS
#	exit
mpifrt  -o fftw.ex $FFLAGS  main.F90  $LDFLAGS

export LD_LIBRARY_PATH=${FFTW_DIR}/lib:${LD_LIBRARY_PATH}
export OMP_STACKSIZE=32M

#	for plan in FFTW_ESTIMATE FFTW_MEASURE FFTW_PATIENT FFTW_EXHAUSTIVE
for plan in FFTW_ESTIMATE FFTW_MEASURE
do
for nthread in 1 2 4 6 8 10 12
do
export MY_FFTW_PLAN=${plan}
export OMP_NUM_THREADS=${nthread}
./fftw.ex
done
done



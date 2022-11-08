#! /bin/bash
#PJM -N FFTW3-1D-C-SERIAL
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:05:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -S

module list

set -x
hostname
date
WD=${HOME}/tmp/check_fftw3_1D_serial
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

SRC=${HOME}/fftw/src_fftw3

#	FFTW_DIR=${HOME}/fftw/fftw-3.3.10/root
FFTW_DIR=${HOME}/fftw/github-fujitsu/local

CFLAGS="-Kfast -I${FFTW_DIR}/include "
FFLAGS="-Kfast -I${FFTW_DIR}/include -Nlst=t "
LDFLAGS="-Kfast "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3_omp  -lfftw3_mpi -lfftw3 -lfjprofmpif "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3_omp  -lfftw3_mpi -lfftw3 -lfjprofmpi "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 "
#	LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 --linkstl=libfjc++ "
LDFLAGS+="-L${FFTW_DIR}/lib -lfftw3 --linkfortran "

cp $SRC/fftw3_complex_1D.serial.c main.c
#	fccpx  -o fftw.ex $CFLAGS  main.c  $LDFLAGS	# NG. needing libmpi
mpifccpx  -o fftw.ex $CFLAGS  main.c  $LDFLAGS

export LD_LIBRARY_PATH=${FFTW_DIR}/lib:${LD_LIBRARY_PATH}

./fftw.ex

#	export OMP_NUM_THREADS=1
#	export OMP_STACKSIZE=32M
#	for i in FFTW_ESTIMATE FFTW_MEASURE FFTW_PATIENT FFTW_EXHAUSTIVE
#	do
#	export MY_FFTW_PLAN=${i}
#	./fftw.ex
#	done

exit

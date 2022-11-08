#! /bin/bash
#PJM -N MY-FFTW3-1D-PMLIB
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:30:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -S

module list
export LANG=C
echo job manager reserved the following resources
echo PJM_NODE         : ${PJM_NODE}
echo PJM_MPI_PROC     : ${PJM_MPI_PROC}
echo PJM_PROC_BY_NODE : ${PJM_PROC_BY_NODE}

set -x
hostname
date
WD=${HOME}/tmp/check_fftw3_1d_pmlib
mkdir -p ${WD}
cd ${WD}; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
rm ${WD}/*

#	SRCDIR=${PMLIB_DIR}/doc/src_tutorial
#	SRCDIR=${HOME}/pmlib/src_tests/src_test_threads
SRCDIR=${HOME}/fftw/src
FFTW_DIR=${HOME}/fftw/fftw-3.3.10/root

################################################################
# on Login node
#	PAPI_DIR=/opt/FJSVxos/devkit/aarch64/rfs/usr
# on compute node
PAPI_DIR=/usr
PAPI_LDFLAGS="-L${PAPI_DIR}/lib64 -lpapi -lpfm "
PAPI_INCLUDES="-I${PAPI_DIR}/include "

PMLIB_DIR=${HOME}/pmlib/usr_local_pmlib/fugaku-precise
PMLIB_INCLUDES="-I${PMLIB_DIR}/include "	# needed for C++ and C programs
# Choose MPI (1) or serial (2)
#	PMLIB_LDFLAGS="-L${PMLIB_DIR}/lib -lPMmpi -lpapi_ext -lpower_ext "	# (1) MPI
PMLIB_LDFLAGS="-L${PMLIB_DIR}/lib -lPM -lpapi_ext -lpower_ext "	# (2) serial

POWER_DIR="/opt/FJSVtcs/pwrm/aarch64"
POWER_LDFLAGS="-L${POWER_DIR}/lib64 -lpwr "
POWER_INCLUDES="-I${POWER_DIR}/include "

INCLUDES="${PAPI_INCLUDES} ${POWER_INCLUDES} ${PMLIB_INCLUDES} "
#	MPI_LDFLAGS="-lmpi "
LDFLAGS="${PAPI_LDFLAGS} ${POWER_LDFLAGS} ${PMLIB_LDFLAGS} ${MPI_LDFLAGS} "
#	LDFLAGS+="--linkfortran "	# for C++ linker
#	LDFLAGS+="-lstdc++ "		# for fortran linker
LDFLAGS+="-lstdc++ "		# for fortran linker
################################################################

OMPFLAGS="-Kopenmp "
FFLAGS="-Kfast ${OMPFLAGS} -Cpp ${DEBUG} -I${FFTW_DIR}/include "
CXXFLAGS="-Kfast ${OMPFLAGS} -Nclang --std=c++11 -w ${DEBUG} -I${FFTW_DIR}/include "
CFLAGS="-Kfast ${OMPFLAGS} -Nclang --std=c11 -w ${DEBUG} -I${FFTW_DIR}/include "
LDFLAGS+="-Kfast ${OMPFLAGS}  -L${FFTW_DIR}/lib -lfftw3 -w "

cp $SRCDIR/main_fftw3_complex_1d.pmlib.F main.F90
frt  -o fftw.ex ${FFLAGS} main.F90 ${LDFLAGS}

xospastop
export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE
export XOS_MMM_L_PAGING_POLICY=prepage:demand:demand

export OMP_NUM_THREADS=1
export HWPC_CHOOSER=FLOPS
export PMLIB_REPORT=FULL
export POWER_CHOOSER=NODE

./fftw.ex


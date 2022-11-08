#! /bin/bash
#PJM -N BUILD-FFTW-GITHUB
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=06:00:00"
#PJM --rsc-list "node=1"
#PJM --mpi "max-proc-per-node=4"
#   #PJM --mpi "proc=4"
#PJM -j
#   #PJM -s
#   #PJM -S

module list

set -x
#	SRC_DIR=$HOME/fftw/fftw-3.3.10
SRC_DIR=$HOME/fftw/github-fujitsu/fftw3
INSTALL_PATH=$HOME/fftw/github-fujitsu/local


cd $SRC_DIR; if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi

export CFLAGS="-Kfast,openmp -Nlst=t"
#	export CFLAGS="-Kfast,openmp -Nclang -ffj-list=t"
export FFLAGS="-Kfast,openmp -Nlst=t"
export LDFLAGS="-Kfast,openmp --linkfortran"

# Cross Compile

touch ChangeLog
autoreconf --verbose --install --symlink --force

./configure                        \
	CC=mpifccpx CXX=mpiFCCpx FC=mpifrtpx F77=mpifrtpx \
    --host=aarch64-unknown-linux-gnu \
    --build=x86_64-cross-linux-gnu   \
    --enable-armv8-cntvct-el0        \
    --enable-sve                     \
    --enable-fma                     \
    --enable-fortran                 \
    --enable-openmp                  \
    --enable-mpi                  \
    --enable-shared                  \
    --prefix="$INSTALL_PATH"         \
    ac_cv_prog_f77_v='-###'          \
    OPENMP_CFLAGS='-Kopenmp'

make
make install
exit

# Native Compile
./configure --enable-openmp --enable-threads \
	--enable-mpi --enable-fortran \
	--prefix=${SRC_DIR}/root \
	CC=mpifcc CXX=mpiFCC FC=mpifrt F77=mpifrt

make
make install


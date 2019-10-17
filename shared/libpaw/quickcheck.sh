#!/bin/bash
#
## Copyright (C) 2017-2019 Yann Pouillon <devops@materialsevolution.es>

# Note: this script is for maintainers and testers working with GCC

# Decide whether to check configure options
do_check_config="no"

# Stop at first error and echo commands
set -ev

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "src/m_paw_atom.F90"; then
  echo "This is not a LibPAW source tree - aborting now"
  exit 1
fi

# Make sure the directory tree is writable
chmod -R u+w .

# Check that we are in the correct directory
test -s "configure.ac" -a -s "src/m_paw_atom.F90" || exit 0

# Init build parameters
export DBGFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"

# Prepare source tree
./wipeout.sh
./autogen.sh

# Check configure script
if test "${do_check_config}" = "yes"; then
  mkdir tmp-config
  cd tmp-config
  echo ""
  echo "### BASIC ###"
  ../configure
  sleep 1
  echo ""
  echo "### SERIAL ###"
  ../configure \
    CC="gcc" CFLAGS="${DBGFLAGS}" FC="gfortran" FCFLAGS="${DBGFLAGS}"
  sleep 1
  echo ""
  echo "### MPI(env) ###"
  ../configure \
    CC="mpicc" CFLAGS="${DBGFLAGS}" FC="mpif90" FCFLAGS="${DBGFLAGS}"
  sleep 1
  echo ""
  echo "### MPI(dir) ###"
  ../configure \
    --with-mpi=/usr CFLAGS="${DBGFLAGS}" FCFLAGS="${DBGFLAGS}"
  sleep 1
  echo ""
  echo "### MPI(yon) ###"
  ../configure \
    --with-mpi CFLAGS="${DBGFLAGS}" FCFLAGS="${DBGFLAGS}"
  sleep 1
  echo ""
  echo "### MPI(yon) + MULTICONFIG ###"
  ../configure \
    --enable-multiconfig --with-mpi CFLAGS="${DBGFLAGS}" FCFLAGS="${DBGFLAGS}"
  sleep 1
  cd ..
fi

# Check export
mkdir tmp-export
cd tmp-export
../configure
make export
tar xvzf libpaw-[0-9].[0-9].[0-9].tar.gz
cd libpaw-[0-9].[0-9].[0-9]


# Check default build (requires ABINIT to build due to symbol names)
abinit_common_libs="../common/src/libabinit_common.a"
if test -s "${abinit_common_libs}"; then
  mkdir tmp-minimal
  cd tmp-minimal
  ../configure \
    CC="gcc" CFLAGS="${DBGFLAGS}" FC="gfortran" FCFLAGS="${DBGFLAGS}"
  sleep 1
  make dist
  make
  make check
  make clean
  make -j4
  mkdir install-minimal
  make install DESTDIR="${PWD}/install-minimal"
  ls -lR install-minimal | tee install-minimal.log
  cd ..
fi

# Check parallel build
mkdir tmp-mpi
cd tmp-mpi
../configure \
  CC="mpicc" CFLAGS="${DBGFLAGS}" FC="mpif90" FCFLAGS="${DBGFLAGS}"
sleep 1
make
make clean && make -j4
make check
cd ..

# Make distcheck
mkdir tmp-distcheck
cd tmp-distcheck
../configure \
  CC="gcc" CFLAGS="${DBGFLAGS}" FC="gfortran" FCFLAGS="${DBGFLAGS}"
sleep 1
make distcheck -j4
make distcleancheck
cd ..

# Go back to the top
cd ..
cd ..

# Clean-up the mess
rm -rf tmp-config tmp-export

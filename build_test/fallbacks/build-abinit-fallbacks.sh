#!/bin/bash

# Init
fallbacks_prefix="/home/vtrinquet/Softwares_Packages/abinit/build_test/fallbacks/install_fb/gnu/11.3"

# Find and Unpack tarball
tarfile=`basename $(ls /home/vtrinquet/Softwares_Packages/abinit/build_test/../fallbacks/*.tar.gz)`
source=${tarfile%.tar.gz}

mkdir -p $source && tar -xzf /home/vtrinquet/Softwares_Packages/abinit/build_test/../fallbacks/$tarfile -C $source --strip-components=1
cd $source

# Configure
./configure \
  --prefix="${fallbacks_prefix}" \
  --with-tardir="${HOME}/.abinit/tarballs" \
  --with-linalg-incs="" \
  --with-linalg-libs="-lopenblas" \
  --with-fc-vendor="gnu" \
  --with-fc-version="11.3" \
  --disable-bigdft \
  --disable-atompaw \
  --disable-wannier90 \
  LIBS_NETCDF4_FORTRAN="-ldl -lm -lz" \
  CC="mpicc" \
  CXX="mpic++" \
  FC="mpifort"

make -j 4 install
rc=`echo $?`

if test "$rc" = "0"; then
  printf "$(tput bold)----------------------------------------------------------------------$(tput sgr0)\n\n"
  echo "The fallbacks are now ready to use."; \
  echo "You can link these fallbacks with Abinit by copying the following options to your ac9 file.";

  list_of_fbks=( libxc hdf5 netcdf4 netcdf4_fortran linalg xmlf90 libpsml wannier90 )
  for i in "${list_of_fbks[@]}"; do
    if test "`${fallbacks_prefix}/bin/abinit-fallbacks-config --enabled ${i}`" = "yes"; then
      Prefix=`${fallbacks_prefix}/bin/abinit-fallbacks-config --libs ${i}`
      printf "\n$(tput bold)"
      echo "with_${i}=${Prefix}" | sed '-e s/-L//;  s/\/lib //; s/netcdf4/netcdf/; s/-l.*$//'
      printf "$(tput sgr0)"
    fi
  done
  printf "\n"
else
  printf "We have detected a problem while generating fallbacks : contact Abinit's team\n"
fi

exit

#!/bin/bash

# Init
Fallbacks="bigdft atompaw wannier90 libpsml xmlf90 libxc hdf5 netcdf4 netcdf4_fortran"

Disabled="--disable-bigdft --disable-atompaw --disable-wannier90 --disable-libpsml --disable-xmlf90 --disable-libxc --disable-hdf5 --disable-netcdf4 --disable-netcdf4_fortran"

Optional_fbks="--disable-bigdft --disable-atompaw --disable-wannier90 --disable-libpsml --disable-xmlf90"

if [[ -z $1 ]]; then
   Disabled=$Optional_fbks
   if [ "dir" != "fb" ]; then
	   Disabled+=' --disable-libxc ' 
   fi 	   
   if [ "dir" != "fb" ]  && [ "dir" != "fb" ]  && [  "dir" != "fb" ]; then
	   Disabled+=' --disable-hdf5  --disable-netcdf4 --disable-netcdf4_fortran' 
   fi 	   
fi
while [[ "$#" -gt 0 ]]; do
   if [ ! "`echo -n $Fallbacks | grep -Fo ${1}`" ]; then
       echo "Unknown parameter: $1";
       exit 1
   fi
   Disabled=`echo $Disabled | sed "s/--disable-${1}//"`
   shift
done

# Init
fallbacks_prefix="/home/buildbot/ABINIT3/eos_gnu_13.2_mpich/lygatsika_spectrum-slicing-optim/fallbacks/install_fb/gnu/13.2"

# Find and Unpack tarball
tarfile=`basename $(ls /home/buildbot/ABINIT3/eos_gnu_13.2_mpich/lygatsika_spectrum-slicing-optim/fallbacks/*.tar.gz)`
source=${tarfile%.tar.gz}

mkdir -p $source && tar -xzf /home/buildbot/ABINIT3/eos_gnu_13.2_mpich/lygatsika_spectrum-slicing-optim/fallbacks/$tarfile -C $source --strip-components=1
cd $source

# Configure
./configure \
  --prefix="${fallbacks_prefix}" \
  --with-tardir="${HOME}/.abinit/tarballs" \
  --with-linalg-incs="-m64 -I/opt/intel/oneapi/mkl/2023.1.0/include" \
  --with-linalg-libs="-m64 -L/opt/intel/oneapi/mkl/2023.1.0/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -lcurl" \
  --with-fc-vendor="gnu" \
  --with-fc-version="13.2" \
  `echo $Disabled` \
  LIBS_NETCDF4_FORTRAN="-ldl -lm -lz" \
  CC="/usr/local/mpich-4.2.2_gnu-13.2/bin/mpicc" \
  CXX="/usr/local/mpich-4.2.2_gnu-13.2/bin/mpic++" \
  FC="/usr/local/mpich-4.2.2_gnu-13.2/bin/mpif90"

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

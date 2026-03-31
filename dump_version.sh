#!/bin/bash

####################
FBV_FILE="config/specs/fbversion.conf"
if test ! -f "$FBV_FILE"; then
   FBV_FILE="../$FBV_FILE"
fi
fbv=`egrep "\[.*\]" $FBV_FILE | head -n 1 | sed -e 's/\[/_/;s/\]//'`
FALLBACKS_PATH=${FALLBACKS_HOME}${fbv}/${FB}/
###################

echo "========== PYTHON ===================================================="
echo
def_python_version=$(python --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1)
python2_version=$(python2 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1)
python3_version=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1)
echo "which python : `which python`"
echo
python -V
echo
if [ "${def_python_version}" == "2" ]; then
     echo "cmd : pip2 freeze"
     echo
     pip2 freeze
     echo
     echo "which python3 : `which python3`"
else
     echo "cmd : pip3 freeze"
     echo
     pip3 freeze
fi
echo

echo "========== PERL ======================================================"
perl -v | head -n 4

echo "========== NETCDF ===================================================="
echo
which nc-config >& /dev/null
if [ "`echo $?`" == "0" ]; then
  nc-config --all
else
  if test "`${FALLBACKS_PATH}/netcdf4/default/bin/nc-config`"; then
    ${FALLBACKS_PATH}/netcdf4/default/bin/nc-config --all
  else
    echo "Can't check netcdf..."
  fi
fi
echo

echo "========== NETCDF Fortran ============================================"
echo
which nf-config >& /dev/null
if [ "`echo $?`" == "0" ]; then
  nc-config --all
else
  if test "`${FALLBACKS_PATH}/netcdf4_fortran/default/bin/nf-config`"; then
    ${FALLBACKS_PATH}/netcdf4_fortran/default/bin/nf-config --all
  else
    echo "Can't check netcdf fortran..."
  fi
fi
echo

echo "========== ABINIT links =============================================="
echo
if test "`which ldd`"; then
   echo "cmd = ldd src/98_main/abinit"
   echo
   ldd src/98_main/abinit
else
   echo "cmd = otool -L src/98_main/abinit"
   otool -L src/98_main/abinit
   echo
fi
echo

echo "========== MPI ======================================================="
echo

cmd_build_pref=""
which mpif90 >& /dev/null
if [ "`echo $?`" == "0" ]; then
   cmd_build_pref="mpiexec -n 1 "
   echo -e "cmd = which mpif90 :\n\n     `which mpif90`"
   echo -e "----------------------------------------------------------------------\n"
   echo -e "cmd = mpif90 --version :\n\n     `mpif90 --version`"
   echo -e "----------------------------------------------------------------------\n"
   echo -e "cmd = mpif90 -show :\n\n     `mpif90 -show`"
   echo -e "----------------------------------------------------------------------\n"
else
   echo "   SERIAL version..."
fi
echo

echo "========== ABINIT BUILD =============================================="
echo
echo "cmd = src/98_main/abinit -b"
echo
${cmd_build_pref}./src/98_main/abinit -b
echo

echo "========== FALLBACKS VERSIONS ========================================"
echo
echo "cmd = egrep -m 1  -B 1 -A 11 '[' $FBV_FILE"
egrep -m 1  -B 1 -A 11 "\[" $FBV_FILE
echo

echo "========== SHELL ENV ================================================="
echo
echo "cmd = PrintEnv.sh"
echo
PrintEnv.sh
echo

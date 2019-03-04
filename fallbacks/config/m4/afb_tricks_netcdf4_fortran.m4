# -*- Autoconf -*-
#
# Copyright (C) 2006-2017 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Tricks for external packages
#



# AFB_TRICKS_NETCDF4_FORTRAN(FC_VENDOR,FC_VERSION)
# ----------------------------------------
#
# Applies tricks and workarounds to have the NetCDF4 Fortran
# C libraries correctly linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_NETCDF4_FORTRAN],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], [], [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_netcdf4_fortran_tricks="no"
  afb_netcdf4_fortran_tricky_vars=""
  tmp_netcdf4_fortran_num_tricks=1
  tmp_netcdf4_fortran_cnt_tricks=0

  dnl Configure tricks
  if test "${afb_netcdf4_fortran_cfgflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying NetCDF4 Fortran tricks (vendor: $1, version: $2, flags: config)])

    dnl Internal NetCDF4 Fortran parameters
    CFGFLAGS_NETCDF4_FORTRAN="${CFGFLAGS_NETCDF4_FORTRAN} --enable-large-file-tests --disable-shared"
    if test "${afb_hdf5_ok}" = "yes"; then
      CFGFLAGS_NETCDF4_FORTRAN="${CFGFLAGS_NETCDF4_FORTRAN} --enable-parallel-tests"
    fi

    dnl Finish
    tmp_netcdf4_fortran_cnt_tricks=`expr ${tmp_netcdf4_fortran_cnt_tricks} \+ 1`
    afb_netcdf4_fortran_tricky_vars="${afb_netcdf4_fortran_tricky_vars} CFGFLAGS"
  else
    AC_MSG_NOTICE([CFGFLAGS_NETCDF4_FORTRAN set => skipping NetCDF4 Fortran config tricks])
  fi

  dnl CPP tricks
  if test "${afb_netcdf4_fortran_cppflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying NetCDF4 Fortran tricks (vendor: $1, version: $2, flags: C preprocessing)])

    CPPFLAGS_NETCDF4_FORTRAN="${CPPFLAGS_NETCDF4_FORTRAN} \$(afb_netcdf4_incs) "

    dnl Finish
    tmp_netcdf4_fortran_cnt_tricks=`expr ${tmp_netcdf4_fortran_cnt_tricks} \+ 1`
    afb_netcdf4_fortran_tricky_vars="${afb_netcdf4_fortran_tricky_vars} CFGFLAGS"
  else
    AC_MSG_NOTICE([CPPFLAGS_NETCDF4_FORTRAN set => skipping NetCDF4 Fortran  preprocessing tricks])
  fi

  dnl Count applied tricks
  case "${tmp_netcdf4_fortran_cnt_tricks}" in
    0)
      afb_netcdf4_fortran_tricks="no"
      ;;
    ${tmp_netcdf4_fortran_num_tricks})
      afb_netcdf4_fortran_tricks="yes"
      ;;
    *)
      afb_netcdf4_fortran_tricks="partial"
      ;;
  esac
  unset tmp_netcdf4_fortran_cnt_tricks
  unset tmp_netcdf4_fortran_num_tricks
]) # AFB_TRICKS_NETCDF4_FORTRAN

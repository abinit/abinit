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



# AFB_TRICKS_NETCDF4(FC_VENDOR,FC_VERSION)
# ----------------------------------------
#
# Applies tricks and workarounds to have the NetCDF4
# C libraries correctly linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_NETCDF4],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], [], [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_netcdf4_tricks="no"
  afb_netcdf4_tricky_vars=""
  tmp_netcdf4_num_tricks=1
  tmp_netcdf4_cnt_tricks=0

  dnl Configure tricks
  if test "${afb_netcdf4_cfgflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying NetCDF4 tricks (vendor: $1, version: $2, flags: config)])

    dnl Internal NetCDF4 parameters
    CFGFLAGS_NETCDF4="${CFGFLAGS_NETCDF4} --disable-dap --disable-examples --disable-hdf4 --disable-v2 --disable-shared"
    if test "${afb_hdf5_ok}" = "yes"; then
      CFGFLAGS_NETCDF4="${CFGFLAGS_NETCDF4} --enable-parallel-tests"
    else
      CFGFLAGS_NETCDF4="${CFGFLAGS_NETCDF4} --disable-netcdf-4"
    fi

    dnl Finish
    tmp_netcdf4_cnt_tricks=`expr ${tmp_netcdf4_cnt_tricks} \+ 1`
    afb_netcdf4_tricky_vars="${afb_netcdf4_tricky_vars} CFGFLAGS"
  else
    AC_MSG_NOTICE([CFGFLAGS_NETCDF4 set => skipping NetCDF4 config tricks])
  fi

  dnl C tricks
  if test "${afb_netcdf4_cflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying LibXC tricks (vendor: $1, version: $2, flags: C)])
    CFLAGS_NETCDF4="${CFLAGS_NETCDF4} -fPIC"

    dnl Finish
    tmp_netcdf4_cnt_tricks=`expr ${tmp_netcdf4_cnt_tricks} \+ 1`
    afb_netcdf4_tricky_vars="${afb_netcdf4_tricky_vars} CFLAGS"
  else
    AC_MSG_NOTICE([CFLAGS_NETCDF4 set => skipping LibXC C tricks])
  fi

  dnl Count applied tricks
  case "${tmp_netcdf4_cnt_tricks}" in
    0)
      afb_netcdf4_tricks="no"
      ;;
    ${tmp_netcdf4_num_tricks})
      afb_netcdf4_tricks="yes"
      ;;
    *)
      afb_netcdf4_tricks="partial"
      ;;
  esac
  unset tmp_netcdf4_cnt_tricks
  unset tmp_netcdf4_num_tricks
]) # AFB_TRICKS_NETCDF4

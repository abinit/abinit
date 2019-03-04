# -*- Autoconf -*-
#
# Copyright (C) 2014-2017 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_NETCDF4()
# -------------------
#
# Check whether the specified NetCDF4 library is working.
#
AC_DEFUN([AFB_CHECK_NETCDF4],[
  dnl Init
  afb_netcdf4_default_libs="-lnetcdf"
  afb_netcdf4_ok="unknown"
  afb_netcdf4_par_ok="unknown"

  dnl Prepare environment
  tmp_saved_with_netcdf4="${with_netcdf4}"
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${afb_netcdf4_incs}"
  AC_MSG_CHECKING([for NetCDF4 libraries to try])
  if test "${afb_netcdf4_libs}" = ""; then
    LIBS="${afb_netcdf4_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_netcdf4_default_libs}])
  else
    LIBS="${afb_netcdf4_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_netcdf4_libs}])
  fi

  dnl Use existing NetCDF4 macro from Autoconf Archive,
  dnl first parallel, then serial if it fails
  dnl We have to do this because we prefer parallel, while the macro
  dnl prefers serial
  AX_LIB_NETCDF4([parallel])
  if test "${with_netcdf4}" = "no"; then
    AC_MSG_NOTICE([no parallel NetCDF4 found, looking for a serial one])
    with_netcdf4="${tmp_saved_with_netcdf4}"
    AX_LIB_NETCDF4([serial])
  fi

  dnl Set internal variables according to the results of AX_LIB_NETCDF4
  dnl See ax_lib_netcdf4.m4 for more information
  if test "${with_netcdf4}" = "yes"; then
    AC_MSG_CHECKING([which NetCDF4 version we have])
    AC_MSG_RESULT([${NETCDF4_VERSION}])
    afb_netcdf4_ok="yes"
    if test "${afb_netcdf4_incs}" = ""; then
      afb_netcdf4_incs="${NETCDF4_CPPFLAGS}"
    fi
    if test "${afb_netcdf4_libs}" = ""; then
      afb_netcdf4_libs="${NETCDF4_LIBS}"
    fi
  else
    afb_netcdf4_ok="no"
  fi
  if test "${with_netcdf4_parallel}" = "yes"; then
    afb_netcdf4_par_ok="yes"
  else
    afb_netcdf4_par_ok="no"
  fi

  dnl Restore environment
  with_netcdf4="${tmp_saved_with_netcdf4}"
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_NETCDF4

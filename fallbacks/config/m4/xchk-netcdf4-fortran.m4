# -*- Autoconf -*-
#
# Copyright (C) 2014-2017 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_NETCDF4_FORTRAN()
# ---------------------------
#
# Check whether the specified NetCDF4-Fortran library is working.
#
AC_DEFUN([AFB_CHECK_NETCDF4_FORTRAN],[
  dnl Init
  afb_netcdf4_fortran_default_libs="-lnetcdff"
  afb_netcdf4_fortran_ok="unknown"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${afb_netcdf4_fortran_incs}"
  AC_MSG_CHECKING([for NetCDF4-Fortran libraries to try])
  if test "${afb_netcdf4_fortran_libs}" = ""; then
    LIBS="${afb_netcdf4_fortran_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_netcdf4_fortran_default_libs}])
  else
    LIBS="${afb_netcdf4_fortran_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_netcdf4_fortran_libs}])
  fi

  dnl Set internal variables according to the results of AX_LIB_NETCDF4
  dnl if available (see ax_lib_netcdf4.m4 for more information)
  if test "${with_netcdf4_fortran}" = "yes"; then
    AC_MSG_CHECKING([which NetCDF4-Fortran version we have])
    AC_MSG_RESULT([${NETCDF4_VERSION}])
    afb_netcdf4_fortran_ok="yes"
    if test "${afb_netcdf4_fortran_incs}" = ""; then
      afb_netcdf4_fortran_incs="${NETCDF4_FFLAGS}"
    fi
    if test "${afb_netcdf4_fortran_libs}" = ""; then
      afb_netcdf4_fortran_libs="${NETCDF4_FLIBS}"
    fi
  fi

  dnl Explicitly check Fortran support if not confirmed yet
  if test "${afb_netcdf4_fortran_ok}" != "yes"; then
    AC_MSG_CHECKING([whether NetCDF4-Fortran works])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use netcdf
        character(len=*), parameter :: path = "dummy"
        integer :: mode, ncerr, ncid
        ncerr = nf90_open(path,mode,ncid)
      ]])], [afb_netcdf4_fortran_ok="yes"], [afb_netcdf4_fortran_ok="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_netcdf4_fortran_ok}])
  fi

  dnl Set internal variables if successful
  if test "${afb_netcdf4_fortran_ok}" = "yes"; then
    if test "${afb_netcdf4_fortran_libs}" = ""; then
      afb_netcdf4_fortran_libs="${afb_netcdf4_fortran_default_libs}"
    fi
  fi

  dnl Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_NETCDF4_FORTRAN

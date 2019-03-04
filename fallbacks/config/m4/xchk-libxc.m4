# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_LIBXC(API_MAJOR, API_MINOR, API_MAX)
# --------------------------------------------------------
#
# Check whether the specified LibXC library is working.
#
AC_DEFUN([AFB_CHECK_LIBXC],[
  dnl Init
  afb_libxc_default_libs="-lxc"
  afb_libxc_has_hdrs="unknown"
  afb_libxc_has_mods="unknown"
  afb_libxc_has_incs="unknown"
  afb_libxc_has_libs="unknown"
  afb_libxc_version="unknown"
  afb_libxc_ext_ok="yes"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${afb_libxc_incs}"
  FCFLAGS="${FCFLAGS} ${afb_libxc_incs}"
  if test "${afb_libxc_libs}" = ""; then
    AC_MSG_CHECKING([for LibXC libraries to try])
    LIBS="${afb_libxc_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_libxc_default_libs}])
  else
    LIBS="${afb_libxc_libs} ${LIBS}"
  fi

  dnl Look for C includes
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([xc.h xc_funcs.h],[afb_libxc_has_hdrs="yes"],[afb_libxc_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for Fortran includes
  AC_MSG_CHECKING([for LibXC Fortran modules])
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
    [[
      use xc_f90_lib_m
    ]])], [afb_libxc_has_mods="yes"], [afb_libxc_has_mods="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${afb_libxc_has_mods}])

  dnl Check status of include files
  if test "${afb_libxc_has_hdrs}" = "yes" -a \
          "${afb_libxc_has_mods}" = "yes"; then
    afb_libxc_has_incs="yes"
  else
    afb_libxc_has_incs="no"
  fi

  dnl Check whether the Fortran wrappers work
  if test "${afb_libxc_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether LibXC has Fortran support])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use xc_f90_lib_m
        integer :: i
        type(xc_f90_pointer_t) :: info
        i = xc_f90_info_number(info)
      ]])], [afb_libxc_has_libs="yes"], [afb_libxc_has_libs="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_libxc_has_libs}])
  fi

  dnl Check that we have the correct LibXC version
  if test "${afb_libxc_has_incs}" = "yes" -a \
          "${afb_libxc_has_libs}" = "yes"; then
    AC_MSG_CHECKING([whether this is LibXC version $1.$2])
    AC_LANG_PUSH([C])
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
      [[
#include <stdio.h>
#include "xc.h"
      ]],
      [[
        int major = -1, minor = -1, micro = -1;
        FILE *tmp;
        xc_version(&major, &minor, &micro);
        tmp = fopen("conftest.tmp", "w");
        fprintf(tmp, "%d.%d\n", major, minor);
        fclose(tmp);
        if ( (major != $1) || (minor < $2) || (minor > $3) ) {
          return 1; }
      ]])], [afb_libxc_version="yes"], [afb_libxc_version="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${afb_libxc_version}])
    if test "${afb_libxc_version}" = "no" -a -s "conftest.tmp"; then
      tmp_libxc_version_found=`cat conftest.tmp`
      AC_MSG_WARN([found LibXC API version ${tmp_libxc_version_found}])
    fi
    rm -f conftest.tmp
  fi

  dnl Final adjustments
  if test "${afb_libxc_has_incs}" = "yes" -a \
          "${afb_libxc_has_libs}" = "yes" -a \
          "${afb_libxc_version}" = "yes"; then
    afb_libxc_ext_ok="yes"
    if test "${afb_libxc_libs}" = ""; then
      afb_libxc_libs="${afb_libxc_default_libs}"
    fi
  else
    afb_libxc_ext_ok="no"
  fi

  dnl Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_LIBXC

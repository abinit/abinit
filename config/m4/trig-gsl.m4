# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for the GSL library
#



# ABI_TRIGGER_GSL()
# -----------------
#
# Check whether the GSL library is working.
#
AC_DEFUN([ABI_TRIGGER_GSL],[
  dnl Init
  abi_gsl_has_incs="no"
  abi_gsl_has_libs="no"
  abi_gsl_serial="no"
  abi_gsl_mpi="no"
  abi_gsl_fcflags=""
  abi_gsl_ldflags=""
  abi_gsl_incs="${with_gsl_incs}"
  abi_gsl_libs="${with_gsl_libs}"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_gsl_incs} ${CPPFLAGS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${with_gsl_libs} ${LIBS}"
  AC_LANG_PUSH([C])

  dnl Look for a configurator
  AC_CHECK_PROGS([GSL_CONFIG],[gsl-config])
  if test "${GSL_CONFIG}" != ""; then
    if test "${abi_gsl_incs}" = ""; then
      AC_MSG_CHECKING([for GSL include flags])
      abi_gsl_incs=`${GSL_CONFIG} --cflags`
      CPPFLAGS="${abi_gsl_incs} ${CPPFLAGS}"
      AC_MSG_RESULT([${abi_gsl_incs}])
    fi
    if test "${abi_gsl_libs}" = ""; then
      AC_MSG_CHECKING([for GSL link flags])
      abi_gsl_libs=`${GSL_CONFIG} --libs`
      LIBS="${abi_gsl_libs} ${LIBS}"
      AC_MSG_RESULT([${abi_gsl_libs}])
    fi
  fi

  dnl Look for includes
  AC_CHECK_HEADERS([gsl/gsl_sf_gamma.h],[abi_gsl_has_incs="yes"],[abi_gsl_has_incs="no"])

  dnl Look for libraries and routines
  if test "${abi_gsl_libs}" = ""; then
    LIBS="-lgslcblas -lm ${LIBS}"
    AC_CHECK_LIB([gsl],[gsl_sf_gamma],[abi_gsl_has_libs="yes"],[abi_gsl_has_libs="no"])
    if test "${abi_gsl_has_libs}" = "yes"; then
      abi_gsl_libs="-lgsl -lgslcblas -lm"
    fi
  else
    AC_MSG_CHECKING([whether the specified GSL library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#include "gsl/gsl_sf_gamma.h"
      ]],
      [[
        double x,y;
        x = 1.0;
        y = gsl_sf_gamma(x);
      ]])], [abi_gsl_has_libs="yes"], [abi_gsl_has_libs="no"])
    AC_MSG_RESULT([${abi_gsl_has_libs}])
  fi

  dnl Take final decision
  if test "${abi_gsl_has_incs}" = "yes" -a \
          "${abi_gsl_has_libs}" = "yes"; then
    abi_gsl_serial="yes"
    AC_DEFINE([HAVE_GSL],1,[Define to 1 if you have the GNU Scientific Library.])
    abi_gsl_fcflags="${abi_gsl_fcflags}"
    abi_gsl_ldflags="${abi_gsl_ldflags}"
    abi_gsl_incs="${abi_gsl_incs}"
    abi_gsl_libs="${abi_gsl_libs}"
  else
    ABI_MSG_NOTICE([connectors-failure],[GSL detection failure])
    if test "${with_gsl_libs}" = ""; then
      AC_MSG_ERROR([GSL is not available])
    else
      AC_MSG_ERROR([the specified GSL link parameters do not work])
    fi
  fi

  dnl Restore build environment
  AC_LANG_POP([C])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the library
  AC_SUBST(abi_gsl_fcflags)
  AC_SUBST(abi_gsl_ldflags)
  AC_SUBST(abi_gsl_incs)
  AC_SUBST(abi_gsl_libs)
]) # ABI_TRIGGER_GSL

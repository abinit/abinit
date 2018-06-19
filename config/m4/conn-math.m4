# -*- Autoconf -*-
#
# Copyright (C) 2005-2018 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for optimized math libraries
#



# _ABI_MATH_CHECK_GSL()
# ---------------------
#
# Check whether the GSL library is working.
#
AC_DEFUN([_ABI_MATH_CHECK_GSL],[
  dnl Init
  abi_math_gsl_has_incs="no"
  abi_math_gsl_has_libs="no"
  abi_math_gsl_serial="no"
  abi_math_gsl_mpi="no"
  abi_math_gsl_fcflags=""
  abi_math_gsl_ldflags=""
  abi_math_gsl_incs="${with_math_incs}"
  abi_math_gsl_libs="${with_math_libs}"

  dnl Need to switch to C
  AC_LANG_PUSH([C])

  dnl Look for a configurator
  AC_CHECK_PROGS([GSL_CONFIG],[gsl-config])
  if test "${GSL_CONFIG}" != ""; then
    if test "${abi_math_gsl_incs}" = ""; then
      AC_MSG_CHECKING([for GSL include flags])
      abi_math_gsl_incs=`${GSL_CONFIG} --cflags`
      CPPFLAGS="${abi_math_gsl_incs} ${CPPFLAGS}"
      AC_MSG_RESULT([${abi_math_gsl_incs}])
    fi
    if test "${abi_math_gsl_libs}" = ""; then
      AC_MSG_CHECKING([for GSL link flags])
      abi_math_gsl_libs=`${GSL_CONFIG} --libs`
      LIBS="${abi_math_gsl_libs} ${LIBS}"
      AC_MSG_RESULT([${abi_math_gsl_libs}])
    fi
  fi

  dnl Look for includes
  AC_CHECK_HEADERS([gsl/gsl_sf_gamma.h],[abi_math_gsl_has_incs="yes"],[abi_math_gsl_has_incs="no"])

  dnl Look for libraries and routines
  if test "${abi_math_gsl_libs}" = ""; then
    LIBS="-lgslcblas -lm ${LIBS}"
    AC_CHECK_LIB([gsl],[gsl_sf_gamma],[abi_math_gsl_has_libs="yes"],[abi_math_gsl_has_libs="no"])
    if test "${abi_math_gsl_has_libs}" = "yes"; then
      abi_math_gsl_libs="-lgsl -lgslcblas -lm"
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
      ]])], [abi_math_gsl_has_libs="yes"], [abi_math_gsl_has_libs="no"])
    AC_MSG_RESULT([${abi_math_gsl_has_libs}])
  fi

  dnl Take final decision
  if test "${abi_math_gsl_has_incs}" = "yes" -a \
          "${abi_math_gsl_has_libs}" = "yes"; then
    abi_math_gsl_serial="yes"
  fi

  dnl Restore previous language
  AC_LANG_POP([C])
]) # _ABI_MATH_CHECK_GSL



# _ABI_MATH_CHECK_MLIB()
# ----------------------
#
# Check whether the MLIB library is working.
#
AC_DEFUN([_ABI_MATH_CHECK_MLIB],[
  dnl Init
  abi_math_mlib_serial="no"
  abi_math_mlib_mpi="no"
  abi_math_mlib_fcflags=""
  abi_math_mlib_ldflags=""
  abi_math_mlib_incs=""
  abi_math_mlib_libs=""

  dnl Look for libraries and routines
  if test "${with_math_libs}" = ""; then
    AC_CHECK_LIB([veclib],[vrpbrmrb])
    if test "${ac_cv_lib_veclib_vrpbrmrb}" = "yes"; then
      abi_math_mlib_serial="yes"
      abi_math_mlib_libs="-lveclib"
    fi
    if test "${enable_mpi}" = "yes" -a \
            "${abi_math_mlib_serial}" = "yes"; then
      abi_math_mlib_mpi="no"
    fi
  else
    dnl FIXME: implement something
    dnl _ABI_MATH_CHECK_USER
    dnl abi_math_mlib_serial="${abi_math_user_serial}"
    dnl abi_math_mlib_mpi="${abi_math_user_mpi}"
    AC_MSG_WARN([library check not implemented])
    abi_math_mlib_serial="yes"
    if test "${enable_mpi}" = "yes"; then
      abi_math_mlib_mpi="yes"
    fi
    if test "${abi_math_mlib_serial}" = "yes"; then
      abi_math_mlib_incs="${with_math_incs}"
      abi_math_mlib_libs="${with_math_libs}"
    fi
  fi
]) # _ABI_MATH_CHECK_MLIB



# ABI_CONNECT_MATH()
# -------------------
#
# Sets all variables needed to handle the optimized math libraries.
#
AC_DEFUN([ABI_CONNECT_MATH],[
  dnl Initial setup
  lib_math_flavor="${with_math_flavor}"
  lib_math_fcflags=""
  lib_math_ldflags=""
  lib_math_incs=""
  lib_math_libs=""
  abi_math_serial="no"
  abi_math_mpi="no"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_math_incs} ${CPPFLAGS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${with_math_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested math support])
  AC_MSG_RESULT([${with_math_flavor}])

  dnl Look for external math libraries
  if test "${with_math_flavor}" != "none"; then

    case "${with_math_flavor}" in

      custom)
        if test "${with_math_libs}" == ""; then
          AC_MSG_ERROR([you must specify custom math libraries (--with-math-libs)])
        fi
        abi_math_serial="yes"
        abi_math_mpi="yes"
        lib_math_incs="${with_math_incs}"
        lib_math_libs="${with_math_libs}"
        ;;

      gsl)
        _ABI_MATH_CHECK_GSL
        abi_math_serial="${abi_math_gsl_serial}"
        abi_math_mpi="${abi_math_gsl_mpi}"
        if test "${abi_math_serial}" = "yes"; then
          AC_DEFINE([HAVE_GSL],1,[Define to 1 if you have the GNU Scientific Library.])
          lib_math_fcflags="${abi_math_gsl_fcflags}"
          lib_math_ldflags="${abi_math_gsl_ldflags}"
          lib_math_incs="${abi_math_gsl_incs}"
          lib_math_libs="${abi_math_gsl_libs}"
        fi
        ;;

      mlib)
        _ABI_MATH_CHECK_MLIB
        abi_math_serial="${abi_math_mlib_serial}"
        abi_math_mpi="${abi_math_mlib_mpi}"
        if test "${abi_math_serial}" = "yes"; then
          AC_DEFINE([HAVE_LINALG_MLIB],1,[Define to 1 if you have the HP MLib Library.])
          lib_math_fcflags="${abi_math_mlib_fcflags}"
          lib_math_ldflags="${abi_math_mlib_ldflags}"
          lib_math_incs="${abi_math_mlib_incs}"
          lib_math_libs="${abi_math_mlib_libs}"
        fi
        ;;

      *)
        AC_MSG_ERROR([unknown math flavor '${with_math_flavor}'])
        ;;

    esac

  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual math support])
  AC_MSG_RESULT([${lib_math_flavor}])
  if test "${lib_math_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    if test "${with_math_libs}" = ""; then
      AC_MSG_ERROR([the requested ${with_math_flavor} math flavor is not available])
    else
      AC_MSG_ERROR([the specified math libraries do not work])
    fi
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_math_flavor)
  AC_SUBST(lib_math_fcflags)
  AC_SUBST(lib_math_ldflags)
  AC_SUBST(lib_math_incs)
  AC_SUBST(lib_math_libs)
]) # ABI_CONNECT_MATH

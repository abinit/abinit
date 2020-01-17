# -*- Autoconf -*-
#
# Copyright (C) 2011-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for optimized algorithmic libraries
#



# _ABI_ALGO_CHECK_LEVMAR()
# ------------------------
#
# Check whether the LEVMAR library is working.
#
AC_DEFUN([_ABI_ALGO_CHECK_LEVMAR],[
  dnl Init
  abi_algo_levmar_has_incs="no"
  abi_algo_levmar_has_libs="no"
  abi_algo_levmar_serial="no"
  abi_algo_levmar_mpi="no"
  abi_algo_levmar_fcflags=""
  abi_algo_levmar_ldflags=""
  abi_algo_levmar_incs="${with_algo_incs}"
  abi_algo_levmar_libs="${with_algo_libs}"

  dnl Need to switch to C
  AC_LANG_PUSH([C])

  dnl Look for includes
  AC_CHECK_HEADERS([levmar.h],[abi_algo_levmar_has_incs="yes"],[abi_algo_levmar_has_incs="no"])

  dnl Look for libraries and routines
  if test "${abi_algo_levmar_libs}" = ""; then
    AC_CHECK_LIB([levmar],[dlevmar_dif],[abi_algo_levmar_has_libs="yes"],[abi_algo_levmar_has_libs="no"])
    if test "${abi_algo_levmar_has_libs}" = "yes"; then
      abi_algo_levmar_libs="-llevmar"
    fi
  else
    AC_MSG_CHECKING([whether the specified LEVMAR library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#include <stdlib.h>
#include <levmar.h>

        void dfit_function(double *p, double *y, int m, int n, void *adata)
        {
          p = 0;
        }
      ]],
      [[
        int ret;
        int c_npoles = 1;
        int c_nvals  = 1;
        int nparam   = 1;

        double adata[1];
        double p[1];
        double yvals[1];

        double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

        ret=dlevmar_dif(dfit_function, p, yvals, nparam, c_nvals, 5000, \
          opts, info, NULL, NULL, (void *)&adata);
      ]])], [abi_algo_levmar_has_libs="yes"], [abi_algo_levmar_has_libs="no"])
    AC_MSG_RESULT([${abi_algo_levmar_has_libs}])
  fi

  dnl Take final decision
  if test "${abi_algo_levmar_has_incs}" = "yes" -a \
          "${abi_algo_levmar_has_libs}" = "yes"; then
    abi_algo_levmar_serial="yes"
  fi

  dnl Restore previous language
  AC_LANG_POP([C])
]) # _ABI_ALGO_CHECK_LEVMAR



# ABI_CONNECT_ALGO()
# ------------------
#
# Sets all variables needed to handle the optimized algorithmic libraries.
#
AC_DEFUN([ABI_CONNECT_ALGO],[
  dnl Initial setup
  lib_algo_flavor="${with_algo_flavor}"
  lib_algo_fcflags=""
  lib_algo_ldflags=""
  lib_algo_incs=""
  lib_algo_libs=""
  abi_algo_serial="no"
  abi_algo_mpi="no"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_algo_incs} ${CPPFLAGS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${with_algo_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested algorithmic support])
  AC_MSG_RESULT([${with_algo_flavor}])

  dnl Look for external algorithmic libraries
  if test "${with_algo_flavor}" != "none"; then

    case "${with_algo_flavor}" in

      custom)
        if test "${with_algo_incs}" == ""; then
          AC_MSG_ERROR([you must specify custom algorithmic includes (--with-algo-incs)])
        fi
        if test "${with_algo_libs}" == ""; then
          AC_MSG_ERROR([you must specify custom algorithmic libraries (--with-algo-libs)])
        fi
        abi_algo_serial="yes"
        abi_algo_mpi="yes"
        lib_algo_incs="${with_algo_incs}"
        lib_algo_libs="${with_algo_libs}"
        ;;

      levmar)
        if test "${abi_linalg_serial}" == "no"; then
          AC_MSG_ERROR([levmar support only works with external linear algebra libraries])
        fi
        LIBS="${lib_linalg_libs} ${LIBS}"
        _ABI_ALGO_CHECK_LEVMAR
        abi_algo_serial="${abi_algo_levmar_serial}"
        abi_algo_mpi="${abi_algo_levmar_mpi}"
        if test "${abi_algo_serial}" = "yes"; then
          AC_DEFINE([HAVE_LEVMAR],1,[Define to 1 if you have the Levenberg-Marquardt algorithmic library.])
          lib_algo_fcflags="${abi_algo_levmar_fcflags}"
          lib_algo_ldflags="${abi_algo_levmar_ldflags}"
          lib_algo_incs="${abi_algo_levmar_incs}"
          lib_algo_libs="${abi_algo_levmar_libs}"
        fi
        ;;

      *)
        AC_MSG_ERROR([unknown algorithmic flavor '${with_algo_flavor}'])
        ;;

    esac

  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual algorithmic support])
  AC_MSG_RESULT([${lib_algo_flavor}])
  if test "${lib_algo_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    if test "${with_algo_libs}" = ""; then
      AC_MSG_ERROR([the requested ${with_algo_flavor} algorithmic flavor is not available])
    else
      AC_MSG_ERROR([the specified algorithmic libraries do not work])
    fi
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_algo_flavor)
  AC_SUBST(lib_algo_fcflags)
  AC_SUBST(lib_algo_ldflags)
  AC_SUBST(lib_algo_incs)
  AC_SUBST(lib_algo_libs)
]) # ABI_CONNECT_ALGO

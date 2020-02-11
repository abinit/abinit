# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for Levmar
#



# ABI_TRIGGER_LEVMAR()
# --------------------
#
# Check whether the Levmar library is working.
#
AC_DEFUN([ABI_TRIGGER_LEVMAR],[
  dnl Initial setup
  abi_levmar_fcflags=""
  abi_levmar_ldflags=""
  abi_levmar_incs=""
  abi_levmar_libs=""
  abi_levmar_has_incs="no"
  abi_levmar_has_libs="no"
  abi_levmar_serial="no"
  abi_levmar_mpi="no"
  abi_levmar_incs="${with_levmar_incs}"
  abi_levmar_libs="${with_levmar_libs}"

  dnl Check for prerequisites
  if test "${abi_linalg_serial}" = "no"; then
    AC_MSG_ERROR([Levmar support only works with external linear algebra libraries])
  fi

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_linalg_incs} ${with_levmar_incs} ${CPPFLAGS}"
  LDFLAGS="${CC_LDFLAGS}"
  LIBS="${abi_levmar_libs} ${abi_linalg_libs} ${LIBS}"
  AC_LANG_PUSH([C])

  dnl Look for includes
  AC_CHECK_HEADERS([levmar.h],[abi_levmar_has_incs="yes"],[abi_levmar_has_incs="no"])

  dnl Look for libraries and routines
  if test "${abi_levmar_libs}" = ""; then
    AC_CHECK_LIB([levmar],[dlevmar_dif],[abi_levmar_has_libs="yes"],[abi_levmar_has_libs="no"])
    if test "${abi_levmar_has_libs}" = "yes"; then
      abi_levmar_libs="-llevmar"
    fi
  fi

  AC_MSG_CHECKING([whether the specified Levmar library works])
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
    ]])], [abi_levmar_has_libs="yes"], [abi_levmar_has_libs="no"])
  AC_MSG_RESULT([${abi_levmar_has_libs}])

  dnl Take final decision
  if test "${abi_levmar_has_incs}" = "yes" -a \
          "${abi_levmar_has_libs}" = "yes"; then
    abi_levmar_serial="yes"
    AC_DEFINE([HAVE_LEVMAR],1,[Define to 1 if you have the Levenberg-Marquardt algorithmic library.])
    abi_levmar_fcflags="${abi_levmar_fcflags}"
    abi_levmar_ldflags="${abi_levmar_ldflags}"
    abi_levmar_incs="${abi_levmar_incs}"
    abi_levmar_libs="${abi_levmar_libs}"
  fi

  dnl Restore build environment
  AC_LANG_POP([C])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the library
  AC_SUBST(abi_levmar_fcflags)
  AC_SUBST(abi_levmar_ldflags)
  AC_SUBST(abi_levmar_incs)
  AC_SUBST(abi_levmar_libs)
]) # ABI_TRIGGER_LEVMAR

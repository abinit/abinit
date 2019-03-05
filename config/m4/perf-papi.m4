# -*- Autoconf -*-
#
# Copyright (C) 2014-2015 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for the PAPI timer library
#



# ABI_PERF_CHECK_PAPI()
# ------------------
#
# Check whether the PAPI library is working.
#
AC_DEFUN([ABI_PERF_CHECK_PAPI],[
  dnl Init
  abi_papi_has_incs="no"
  abi_papi_has_libs="no"
  abi_papi_serial="no"
  abi_papi_mpi="no"
  abi_papi_fcflags=""
  abi_papi_ldflags=""
  abi_papi_incs="${with_papi_incs}"
  abi_papi_libs="${with_papi_libs}"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_papi_incs} ${CPPFLAGS}"
  LIBS="${with_papi_libs} ${LIBS}"

  dnl Add rt support if available on the machine
  if test "${with_papi_libs}" = ""; then
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADERS([time.h])
    AC_CHECK_LIB([rt], [clock_gettime], [abi_papi_rt_libs="-lrt"], [abi_papi_rt_libs=""])
    AC_CHECK_FUNCS([clock_gettime])
    AC_LANG_POP([C])
    abi_papi_libs="${abi_papi_rt_libs} ${abi_papi_libs}"
  fi

  dnl Look for C headers
  LDFLAGS="${CC_LDFLAGS}"
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([papi.h f90papi.h],[abi_papi_has_incs="yes"],[abi_papi_has_incs="no"])
  if test "${abi_papi_has_incs}" = "yes"; then
    AC_CHECK_HEADERS([papi/papi_sf_gamma.h])
  fi
  AC_LANG_POP([C])

  dnl Look for libraries and routines
  AC_LANG_PUSH([Fortran])
  if test "${abi_papi_libs}" = ""; then
    AC_CHECK_LIB([papi],[PAPIf_library_init],[abi_papi_has_libs="yes"],[abi_papi_has_libs="no"])
    if test "${abi_papi_has_libs}" = "yes"; then
      abi_papi_libs="-lpapi"
    fi
  else
    AC_MSG_CHECKING([whether the specified PAPI library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
#if defined HAVE_F90PAPI_H
#include "papi/papi_sf_gamma.h"
#endif
        call PAPIf_library_init
      ]])], [abi_papi_has_libs="yes"], [abi_papi_has_libs="no"])
    AC_MSG_RESULT([${abi_papi_has_libs}])
  fi
  AC_LANG_POP([Fortran])

  dnl Take final decision
  if test "${abi_papi_has_incs}" = "yes" -a \
          "${abi_papi_has_libs}" = "yes"; then
    abi_papi_serial="yes"
    AC_DEFINE([HAVE_PAPI],1,[Define to 1 if you have the PAPI library.])
  else
    ABI_MSG_NOTICE([connectors-failure],[PAPI detection failure])
    if test "${with_papi_libs}" = ""; then
      AC_MSG_ERROR([PAPI support is not available])
    else
      AC_MSG_ERROR([the specified PAPI libraries do not work])
    fi
  fi

  dnl Restore build environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the library
  AC_SUBST(abi_papi_fcflags)
  AC_SUBST(abi_papi_ldflags)
  AC_SUBST(abi_papi_incs)
  AC_SUBST(abi_papi_libs)
]) # ABI_PERF_CHECK_PAPI

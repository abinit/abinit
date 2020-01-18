# -*- Autoconf -*-
#
# Copyright (C) 2005-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for optimized timer libraries
#



# _ABI_TIMER_CHECK_PAPI()
# -----------------------
#
# Check whether the PAPI library is working.
#
AC_DEFUN([_ABI_TIMER_CHECK_PAPI],[
  dnl Init
  abi_timer_papi_has_incs="no"
  abi_timer_papi_has_libs="no"
  abi_timer_papi_serial="no"
  abi_timer_papi_mpi="no"
  abi_timer_papi_fcflags=""
  abi_timer_papi_ldflags=""
  abi_timer_papi_incs="${with_timer_incs}"
  abi_timer_papi_libs="${with_timer_libs}"

  dnl Look for includes
  AC_LANG_PUSH([C])
  dnl AC_CHECK_HEADERS([f90papi.h],[abi_timer_papi_has_incs="yes"],[abi_timer_papi_has_incs="no"])
  AC_CHECK_HEADERS([papi.h],[abi_timer_papi_has_incs="yes"],[abi_timer_papi_has_incs="no"])
  AC_LANG_POP([C])

  dnl Look for libraries and routines
  if test "${abi_timer_papi_libs}" = ""; then
    AC_CHECK_LIB([papi],[PAPIf_library_init],[abi_timer_papi_has_libs="yes"],[abi_timer_papi_has_libs="no"])
    if test "${abi_timer_papi_has_libs}" = "yes"; then
      abi_timer_papi_libs="-lpapi"
    fi
  else
    AC_MSG_CHECKING([whether the specified PAPI library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
#if defined HAVE_F90PAPI_H
#include "papi/papi_sf_gamma.h"
#endif
        call PAPIf_library_init
      ]])], [abi_timer_papi_has_libs="yes"], [abi_timer_papi_has_libs="no"])
    AC_MSG_RESULT([${abi_timer_papi_has_libs}])
  fi

  dnl Take final decision
  if test "${abi_timer_papi_has_incs}" = "yes" -a \
          "${abi_timer_papi_has_libs}" = "yes"; then
    abi_timer_papi_serial="yes"
  fi
]) # _ABI_TIMER_CHECK_PAPI



# ABI_CONNECT_TIMER()
# -------------------
#
# Sets all variables needed to handle the optimized timer libraries.
#
AC_DEFUN([ABI_CONNECT_TIMER],[
  dnl Initial setup
  lib_timer_flavor="${with_timer_flavor}"
  lib_timer_fcflags=""
  lib_timer_ldflags=""
  lib_timer_incs=""
  lib_timer_libs=""
  abi_timer_serial="no"
  abi_timer_mpi="no"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_timer_incs} ${CPPFLAGS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${with_timer_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested timer support])
  AC_MSG_RESULT([${with_timer_flavor}])

  dnl Look for external timer libraries
  if test "${with_timer_flavor}" != "none"; then

    case "${with_timer_flavor}" in

      abinit)
        AC_DEFINE([HAVE_TIMER_ABINIT],1,[Define to 1 if you want to use the Abinit timer.])
        abi_timer_serial="yes"
        if test "${enable_mpi}" = "yes"; then
          abi_timer_mpi="yes"
        fi
        ;;

      custom)
        if test "${with_timer_libs}" == ""; then
          AC_MSG_ERROR([you must specify custom timer libraries (--with-timer-libs)])
        fi
        abi_timer_serial="yes"
        abi_timer_mpi="yes"
        lib_timer_incs="${with_timer_incs}"
        lib_timer_libs="${with_timer_libs}"
        ;;

      gptl)
        AC_MSG_ERROR([not implemented])
        ;;

      papi)
        _ABI_TIMER_CHECK_PAPI
        abi_timer_serial="${abi_timer_papi_serial}"
        abi_timer_mpi="${abi_timer_papi_mpi}"
        if test "${abi_timer_serial}" = "yes"; then
          AC_DEFINE([HAVE_PAPI],1,[Define to 1 if you have the PAPI library.])
          lib_timer_fcflags="${abi_timer_papi_fcflags}"
          lib_timer_ldflags="${abi_timer_papi_ldflags}"
          lib_timer_incs="${abi_timer_papi_incs}"
          lib_timer_libs="${abi_timer_papi_libs}"
        fi
        ;;

        *)
          AC_MSG_ERROR([unknown timer flavor '${with_timer_flavor}'])
          ;;

    esac

  fi

  dnl Add rt support if available on the machine.
  AC_LANG_PUSH(C)
  AC_CHECK_HEADERS([time.h])
  AC_CHECK_LIB([rt], [clock_gettime], [abi_timer_rt_libs="-lrt"], [abi_timer_rt_libs=""])
  lib_timer_libs=$lib_timer_libs" "$abi_timer_rt_libs
  AC_CHECK_FUNCS([clock_gettime])
  AC_LANG_POP(C)

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual timer support])
  AC_MSG_RESULT([${lib_timer_flavor}])
  if test "${lib_timer_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    if test "${with_timer_libs}" = ""; then
      AC_MSG_ERROR([the requested ${with_timer_flavor} timer flavor is not available])
    else
      AC_MSG_ERROR([the specified timer libraries do not work])
    fi
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_timer_flavor)
  AC_SUBST(lib_timer_fcflags)
  AC_SUBST(lib_timer_ldflags)
  AC_SUBST(lib_timer_incs)
  AC_SUBST(lib_timer_libs)
]) # ABI_CONNECT_TIMER

# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Comprehensive exchange-correlation functionals library - LibXC
#



# ABI_LIBXC_DETECT(API_MAJOR, API_MINOR_LOW)
# ------------------------------------------
#
# Check whether the LibXC library is working.
#
AC_DEFUN([ABI_LIBXC_DETECT],[
  dnl Check parameters
  dnl m4_if([$1],[],AC_MSG_ERROR([$0 requires API major version number]))
  dnl m4_if([$2],[],AC_MSG_ERROR([$0 requires API minimum minor version number]))

  dnl Init
  abi_libxc_has_hdrs="unknown"
  abi_libxc_has_libs="unknown"
  abi_libxc_has_mods="unknown"
  abi_libxc_has_fort="unknown"
  abi_libxc_version="unknown"
  abi_libxc_ok="unknown"
  abi_libxc_fcflags=""
  abi_libxc_ldflags=""
  abi_libxc_incs="${with_libxc_incs}"
  abi_libxc_libs="${with_libxc_libs}"

  dnl Display input parameters
  AC_MSG_CHECKING([whether LibXC includes have been specified])
  if test "${with_libxc_incs}" = ""; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi
  AC_MSG_CHECKING([whether LibXC libraries have been specified])
  if test "${with_libxc_libs}" = ""; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_libxc_incs}"
  FCFLAGS="${FCFLAGS} ${abi_libxc_incs}"
  LIBS="${abi_libxc_libs} ${LIBS}"

  dnl Look for C includes
  LDFLAGS="${CC_LDFLAGS}"
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([xc.h xc_funcs.h],
    [abi_libxc_has_hdrs="yes"],[abi_libxc_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for C libraries and routines
  if test "${abi_libxc_has_hdrs}" = "yes"; then
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    if test "${abi_libxc_libs}" = ""; then
      AC_SEARCH_LIBS([xc_func_init],[xc],
        [abi_libxc_has_libs="c"], [abi_libxc_has_libs="no"])
      if test "${ac_cv_search_xc_func_init}" != "no"; then
        if test "${ac_cv_search_xc_func_init}" != "none required"; then
          abi_libxc_libs="${ac_cv_search_xc_func_init}"
        fi
      fi
    else
      AC_MSG_CHECKING([whether specified LibXC C libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [[
#include <xc.h>
        ]],
        [[
          struct xc_func_type *p;
          xc_func_init(p, 0, 1);
          return 0;
        ]])], [abi_libxc_has_libs="c"], [abi_libxc_has_libs="no"])
      tmp_libxc_has_libs=`echo "${abi_libxc_has_libs}" | sed -e 's/^c$/yes/'`
      AC_MSG_RESULT([${tmp_libxc_has_libs}])
      unset tmp_libxc_has_libs
    fi
    AC_LANG_POP([C])
  fi

  dnl Check that we have the correct LibXC version
  if test "${abi_libxc_has_libs}" = "c"; then
    AC_MSG_CHECKING([whether this is LibXC >= $1.$2])
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
      [[
#include "xc.h"
      ]],
      [[
        int major = -1, minor = -1;
        xc_version(&major, &minor);
        if ( (major < $1) || ((major == $1) && (minor < $2)) ) {
          return 1; }
      ]])], [abi_libxc_version="yes"], [abi_libxc_version="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${abi_libxc_version}])
  fi

  dnl Look for Fortran libraries and routines
  if test "${abi_libxc_has_libs}" = "c" -a \
          "${abi_libxc_version}" = "yes"; then
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    if test "${with_libxc_libs}" = ""; then
      AC_SEARCH_LIBS([xc_f90_info_number], [xcf90],
        [abi_libxc_has_libs="yes"], [abi_libxc_has_libs="no"])
      if test "${ac_cv_search_xc_f90_info_number}" != "no"; then
        if test "${ac_cv_search_xc_f90_info_number}" != "none required"; then
          abi_libxc_libs="${ac_cv_search_xc_f90_info_number} ${abi_libxc_libs}"
        fi
      fi
    else
      AC_MSG_CHECKING([whether specified LibXC Fortran libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          call xc_f90_info_number
        ]])], [abi_libxc_has_libs="yes"], [abi_libxc_has_libs="no"])
      AC_MSG_RESULT([${abi_libxc_has_libs}])
    fi
    AC_LANG_POP([Fortran])
  fi

  dnl Look for Fortran includes
  dnl Note: must be done after the libraries have been discovered
  if test "${abi_libxc_has_libs}" = "yes"; then
    LDFLAGS="${FC_LDFLAGS}"
    ABI_FC_MOD_INCS([xc_f90_lib_m])
    FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
    if test "${abi_fc_mod_incs_ok}" = "yes"; then
      abi_libxc_has_mods="yes"
      if test "${abi_libxc_incs}" = ""; then
        abi_libxc_incs="${fc_mod_incs}"
      fi
    fi
  fi

  dnl Check whether the Fortran wrappers work
  if test "${abi_libxc_has_mods}" = "yes"; then
    AC_MSG_CHECKING([whether LibXC Fortran wrappers work])
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use xc_f90_lib_m
        integer :: i
        type(xc_f90_pointer_t) :: info
        i = xc_f90_info_number(info)
      ]])], [abi_libxc_has_fort="yes"], [abi_libxc_has_fort="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_libxc_has_fort}])
  fi

  dnl Set status of minimum LibXC support
  if test "${abi_libxc_has_hdrs}" = "yes" -a \
          "${abi_libxc_has_libs}" = "yes" -a \
          "${abi_libxc_has_mods}" = "yes" -a \
          "${abi_libxc_has_fort}" = "yes" -a \
          "${abi_libxc_version}" = "yes"; then
    abi_libxc_ok="yes"
  else
    abi_libxc_ok="no"
  fi

  dnl Define preprocessing options
  if test "${abi_libxc_ok}" = "yes" -o \
          "${enable_fallbacks}" = "yes"; then
    AC_DEFINE([HAVE_LIBXC], 1,
      [Define to 1 if you have the LibXC library.])
  fi

  dnl Propagate information to fallbacks
  if test "${abi_libxc_ok}" = "yes"; then
    abi_libxc_fallback="no"
  else
    if test "${with_libxc_incs}" = "" -a "${with_libxc_libs}" = ""; then
      AC_MSG_WARN([falling back to limited developer LibXC version])
      abi_libxc_fallback="yes"
      abi_fallbacks="${abi_fallbacks} libxc"
      abi_libxc_fcflags=""
      abi_libxc_ldflags=""
      abi_libxc_incs=""
      abi_libxc_libs=""
    else
      ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
      AC_MSG_ERROR([external LibXC support does not work])
    fi
  fi

  dnl Restore build environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Final report
  AC_MSG_CHECKING([for final LibXC Fortran flags])
  AC_MSG_RESULT(['${abi_libxc_fcflags}'])
  AC_MSG_CHECKING([for final LibXC link flags])
  AC_MSG_RESULT(['${abi_libxc_ldflags}'])
  AC_MSG_CHECKING([for final LibXC include flags])
  AC_MSG_RESULT(['${abi_libxc_incs}'])
  AC_MSG_CHECKING([for final LibXC library flags])
  AC_MSG_RESULT(['${abi_libxc_libs}'])
  AC_MSG_CHECKING([whether we need a LibXC fallback])
  AC_MSG_RESULT([${abi_libxc_fallback}])
  if test "${abi_libxc_fallback}" = "yes"; then
    AC_MSG_WARN([fallbacks have been requested
                  you should NOT run production calculations])
  fi

  dnl Substitute variables needed for the use of the libraries
  AC_SUBST(abi_libxc_fcflags)
  AC_SUBST(abi_libxc_ldflags)
  AC_SUBST(abi_libxc_incs)
  AC_SUBST(abi_libxc_libs)
]) # ABI_LIBXC_DETECT

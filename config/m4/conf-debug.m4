# -*- Autoconf -*-
#
# Copyright (C) 2009-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Debugging
#



# ABI_DEBUG_INIT(DBG_MODE)
# ------------------------
#
# Sets debugging parameters according to the requested mode.
#
AC_DEFUN([ABI_DEBUG_INIT],[
  dnl Init
  abi_debug_flavor="$1"
  abi_debug_verbose="no"

  dnl Display debugging status
  AC_MSG_CHECKING([debugging status])
  case "${abi_debug_flavor}" in
    none)
      AC_MSG_RESULT([disabled])
      CPPFLAGS_DEBUG=""
      CFLAGS_DEBUG=""
      CC_LDFLAGS_DEBUG=""
      CC_LIBS_DEBUG=""
      CXXFLAGS_DEBUG=""
      CXX_LDFLAGS_DEBUG=""
      CXX_LIBS_DEBUG=""
      FCFLAGS_DEBUG=""
      FC_LDFLAGS_DEBUG=""
      FC_LIBS_DEBUG=""
      ARFLAGS_DEBUG=""
      ;;
    custom)
      AC_MSG_RESULT([enabled (using user-specified debug flags)])
      ;;
    *)
      AC_MSG_RESULT([enabled (profile: ${abi_debug_flavor})])
      CPPFLAGS_DEBUG=""
      CFLAGS_DEBUG=""
      CC_LDFLAGS_DEBUG=""
      CC_LIBS_DEBUG=""
      CXXFLAGS_DEBUG=""
      CXX_LDFLAGS_DEBUG=""
      CXX_LIBS_DEBUG=""
      FCFLAGS_DEBUG=""
      FC_LDFLAGS_DEBUG=""
      FC_LIBS_DEBUG=""
      ARFLAGS_DEBUG=""
      ;;
  esac

  dnl Get debug flags from database for all profiles
  if test "${abi_debug_flavor}" != "none" -a \
          "${abi_debug_flavor}" != "custom"; then

    dnl Set debug flags for the C compiler
    if test "${CFLAGS}" = ""; then
      if test "${ac_cv_prog_cc_g}" = "yes"; then
        CFLAGS_DEBUG="-g"
        AC_MSG_NOTICE([setting C debug flags to '-g'])
      fi
      ABI_CC_DBGFLAGS
    else
      AC_MSG_NOTICE([debugging profile overriden by user-defined CFLAGS])
    fi

    dnl Set debug flags for the C++ compiler
    if test "${CXXFLAGS}" = ""; then
      if test "${ac_cv_prog_cxx_g}" = "yes"; then
        CXXFLAGS_DEBUG="-g"
        AC_MSG_NOTICE([setting C++ debug flags to '-g'])
      fi
      ABI_CXX_DBGFLAGS
    else
      AC_MSG_NOTICE([debugging profile overriden by user-defined CXXFLAGS])
    fi

    dnl Set debug flags for the Fortran compiler
    if test "${FCFLAGS}" = ""; then
      if test "${ac_cv_prog_fc_g}" = "yes"; then
        FCFLAGS_DEBUG="-g"
        AC_MSG_NOTICE([setting Fortran debug flags to '-g'])
      fi
      ABI_FC_DBGFLAGS
    else
      AC_MSG_NOTICE([debugging profile overriden by user-defined FCFLAGS])
    fi

  fi

  dnl Enable source and/or verbose debugging for selected profiles
  case "${abi_debug_flavor}" in
    verbose|enhanced)
      abi_debug_verbose="yes"
      ;;
    paranoid|naughty)
      abi_source_debug_enable="yes"
      abi_debug_verbose="yes"
      ;;
  esac

  dnl Define DEBUG_MODE preprocessing option
  AC_MSG_CHECKING([whether to activate debug mode in source files])
  if test "${abi_source_debug_enable}" = "yes"; then
    AC_DEFINE([DEBUG_MODE],1,
      [Define to 1 to build debugging instructions in the source code.])
  fi
  AC_MSG_RESULT([${abi_source_debug_enable}])

  dnl Define DEBUG_VERBOSE preprocessing option
  AC_MSG_CHECKING([whether to activate verbose debug messages in source files])
  if test "${abi_debug_verbose}" = "yes"; then
    AC_DEFINE([DEBUG_VERBOSE],1,
      [Define to 1 to turn on verbose debug messages in the source code.])
  fi
  AC_MSG_RESULT([${abi_debug_verbose}])
]) # ABI_DEBUG_INIT

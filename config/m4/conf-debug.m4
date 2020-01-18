# -*- Autoconf -*-
#
# Copyright (C) 2009-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Debugging
#



# ABI_DEBUG_INIT(DBG_MODE,OPT_MODE)
# ---------------------------------
#
# Sets debugging parameters according to the requested mode.
#
AC_DEFUN([ABI_DEBUG_INIT],[
  dnl Init
  abi_debug_mode="$1"
  abi_optim_mode="$2"

  dnl Display debugging status
  AC_MSG_CHECKING([debugging status])
  case "${abi_debug_mode}" in
    no)
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
    yes)
      AC_MSG_RESULT([enabled (using user-specified debug flags)])
      ;;
    *)
      AC_MSG_RESULT([enabled (profile mode: ${abi_debug_mode})])
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

  dnl Init debug flags
  if test "${abi_debug_mode}" != "no" -a "${abi_debug_mode}" != "yes"; then

    dnl Init debug flags for the C compiler
    if test "${CFLAGS}" = ""; then
      if test "${ac_cv_prog_cc_g}" = "yes"; then
        CFLAGS_DEBUG="-g"
        AC_MSG_NOTICE([setting C debug flags to '-g'])
      fi
    fi

    dnl Init debug flags for the C++ compiler
    if test "${CXXFLAGS}" = ""; then
      if test "${ac_cv_prog_cxx_g}" = "yes"; then
        CXXFLAGS_DEBUG="-g"
        AC_MSG_NOTICE([setting C++ debug flags to '-g'])
      fi
    fi

    dnl Init debug flags for the Fortran compiler
    if test "${FCFLAGS}" = ""; then
      if test "${ac_cv_prog_fc_g}" = "yes"; then
        FCFLAGS_DEBUG="-g"
        AC_MSG_NOTICE([setting Fortran debug flags to '-g'])
      fi
    fi
  fi

  dnl Define DEBUG_MODE preprocessing option
  AC_MSG_CHECKING([whether to activate debug mode in source files])
  src_debug_mode="no"
  if test \( "${abi_debug_mode}" != "no" -a "${abi_debug_mode}" != "yes" -a \
             "${abi_debug_mode}" != "basic" -a \
             "${abi_debug_mode}" != "verbose" \) \
       -o \( "${abi_debug_mode}" = "yes" -a "${abi_optim_mode}" = "no" \); then
    AC_DEFINE([DEBUG_MODE],1,[Define to 1 to turn on debug mode.])
    src_debug_mode="yes"
  fi
  AC_MSG_RESULT([${src_debug_mode}])

  dnl Define DEBUG_VERBOSE preprocessing option
  AC_MSG_CHECKING([whether to activate debug verbosity in source files])
  src_debug_verbose="no"
  if test "${src_debug_mode}" = "yes" -o \
          "${abi_debug_mode}" = "verbose"; then
    AC_DEFINE([DEBUG_VERBOSE],1,[Define to 1 to turn on verbose debug messages.])
    src_debug_verbose="yes"
  fi
  AC_MSG_RESULT([${src_debug_verbose}])
]) # ABI_DEBUG_INIT

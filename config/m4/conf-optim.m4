# -*- Autoconf -*-
#
# Copyright (C) 2009-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Optimizations
#



# ABI_OPTIM_INIT(OPT_MODE)
# ------------------------
#
# Sets optimization parameters according to the requested mode.
#
AC_DEFUN([ABI_OPTIM_INIT],[
  dnl Display optimization status
  AC_MSG_CHECKING([optimization status])
  case "${abi_optim_flavor}" in
    none)
      AC_MSG_RESULT([disabled])
      CPPFLAGS_OPTIM=""
      CFLAGS_OPTIM="-O0"
      CC_LDFLAGS_OPTIM=""
      CC_LIBS_OPTIM=""
      CXXFLAGS_OPTIM="-O0"
      CXX_LDFLAGS_OPTIM=""
      CXX_LIBS_OPTIM=""
      FCFLAGS_OPTIM="-O0"
      FC_LDFLAGS_OPTIM=""
      FC_LIBS_OPTIM=""
      ARFLAGS_OPTIM=""
      ;;
    custom)
      AC_MSG_RESULT([enabled (using user-specified flags)])
      ;;
    *)
      AC_MSG_RESULT([enabled (profile: ${abi_optim_flavor})])
      CPPFLAGS_OPTIM=""
      CFLAGS_OPTIM=""
      CC_LDFLAGS_OPTIM=""
      CC_LIBS_OPTIM=""
      CXXFLAGS_OPTIM=""
      CXX_LDFLAGS_OPTIM=""
      CXX_LIBS_OPTIM=""
      FCFLAGS_OPTIM=""
      FC_LDFLAGS_OPTIM=""
      FC_LIBS_OPTIM=""
      ARFLAGS_OPTIM=""
      ;;
  esac

  dnl Get optimization flags from database for all profiles
  if test "${with_optim_flavor}" != "none" -a \
          "${with_optim_flavor}" != "custom"; then
    if test "${CFLAGS}" = ""; then
      ABI_CC_OPTFLAGS
    else
      AC_MSG_NOTICE([optimization profile overriden by user-defined CFLAGS])
    fi
    if test "${CXXFLAGS}" = ""; then
      ABI_CXX_OPTFLAGS
    else
      AC_MSG_NOTICE([optimization profile overriden by user-defined CXXFLAGS])
    fi
    if test "${FCFLAGS}" = ""; then
      ABI_FC_OPTFLAGS
    else
      AC_MSG_NOTICE([optimization profile overriden by user-defined FCFLAGS])
    fi
  fi
]) # ABI_OPTIM_INIT

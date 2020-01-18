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
  dnl Init
  abi_optim_mode="$1"

  dnl Display optimization status
  AC_MSG_CHECKING([optimization status])
  case "${abi_optim_mode}" in
    no)
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
    yes)
      AC_MSG_RESULT([enabled (using user-specified flags)])
      ;;
    *)
      AC_MSG_RESULT([enabled (profile mode: ${abi_optim_mode})])
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
]) # ABI_OPTIM_INIT

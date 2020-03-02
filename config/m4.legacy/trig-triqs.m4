# -*- Autoconf -*-
#
# Copyright (C) 2005-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for TRIQS libraries
#



# ABI_TRIGGER_TRIQS()
# -------------------
#
# Check whether the TRIQS library is working.
#
AC_DEFUN([ABI_TRIGGER_TRIQS], [
  dnl Init
  abi_triqs_fcflags=""
  abi_triqs_ldflags=""
  abi_triqs_has_incs="no"
  abi_triqs_has_libs="no"
  abi_triqs_incs="${with_triqs_incs}"
  abi_triqs_libs="${with_triqs_libs}"

  if test "${enable_triqs_v2_0}" = "yes" || test "${enable_triqs_v1_4}" = "yes"; then

    AC_MSG_WARN([TRIQS detection is not implemented - FIXME])

  fi

  dnl Inform Automake
  AM_CONDITIONAL([DO_BUILD_67_TRIQS_EXT], [test "${enable_triqs_v2_0}" = "yes" || test "${enable_triqs_v1_4}" = "yes"])

  dnl Substitute variables
  AC_SUBST(abi_triqs_fcflags)
  AC_SUBST(abi_triqs_ldflags)
  AC_SUBST(abi_triqs_incs)
  AC_SUBST(abi_triqs_libs)
]) # ABI_TRIGGER_TRIQS

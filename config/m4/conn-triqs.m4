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



# ABI_CONNECT_TRIQS()
# -------------------
#
# Check whether the TRIQS library is working.
#
AC_DEFUN([ABI_CONNECT_TRIQS], [
  dnl Init
  lib_triqs_incs=""
  lib_triqs_libs=""

  if test "${enable_triqs_v2_0}" = "yes" || test "${enable_triqs_v1_4}" = "yes"; then

    lib_triqs_incs="${with_triqs_incs}"
    lib_triqs_libs="${with_triqs_libs}"

    dnl FIXME: M4 code to detect external libraries should go here

  fi

  dnl Inform Automake
  AM_CONDITIONAL([DO_BUILD_67_TRIQS_EXT], [test "${enable_triqs_v2_0}" = "yes" || test "${enable_triqs_v1_4}" = "yes"])

  dnl Substitute variables
  AC_SUBST(lib_triqs_incs)
  AC_SUBST(lib_triqs_libs)
  AC_SUBST(lib_triqs_fcflags)
  AC_SUBST(lib_triqs_ldflags)

]) # ABI_CONNECT_TRIQS

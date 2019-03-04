# -*- Autoconf -*-
#
# Copyright (C) 2006-2017 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Tricks for external packages
#



# AFB_TRICKS_LINALG(FC_VENDOR,FC_VERSION)
# ---------------------------------------
#
# Applies tricks and workarounds to have the optimized linear algebra
# libraries correctly linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_LINALG],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], [], [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_linalg_tricks="no"
  afb_linalg_tricky_vars=""

  dnl Fortran tricks
  if test "${afb_linalg_fcflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying linear algebra tricks (vendor: $1, version: $2, flags: Fortran)])

    case "$1" in
      ibm)
        FCFLAGS_LINALG="${FCFLAGS_LINALG} ${FCFLAGS_FIXEDFORM}"
        ;;
    esac

    dnl Finish
    afb_linalg_tricks="yes"
    afb_linalg_tricky_vars="${afb_linalg_tricky_vars} FCFLAGS"
  else
    AC_MSG_NOTICE([FCFLAGS_LINALG set => skipping linear algebra Fortran tricks])
  fi
]) # AFB_TRICKS_LINALG

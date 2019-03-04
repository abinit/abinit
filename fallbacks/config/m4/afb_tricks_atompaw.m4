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



# AFB_TRICKS_ATOMPAW(FC_VENDOR,FC_VERSION)
# ----------------------------------------
#
# Applies tricks and workarounds to have the AtomPAW library correctly
# linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_ATOMPAW],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], [], [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_atompaw_tricks="no"
  afb_atompaw_tricky_vars=""

  dnl Configure tricks
  if test "${afb_atompaw_cfgflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying AtomPAW tricks (vendor: $1, version: $2, flags: config)])

    dnl Linear algebra
    tmpflags_atompaw='--with-linalg-libs="$(afb_linalg_libs)"'
    CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} ${tmpflags_atompaw}"

    dnl LibXC
    CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} --enable-libxc"
    tmpflags_atompaw='--with-libxc-incs="$(afb_libxc_incs)"'
    CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} ${tmpflags_atompaw}"
    tmpflags_atompaw='--with-libxc-libs="$(afb_libxc_libs)"'
    CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} ${tmpflags_atompaw}"

    dnl Force static build (shared libraries fail to build)
    CFGFLAGS_ATOMPAW="${CFGFLAGS_ATOMPAW} --enable-static --disable-shared"

    dnl Finish
    test "${afb_atompaw_tricks}" = "no" && afb_atompaw_tricks="yes"
    afb_atompaw_tricky_vars="${afb_atompaw_tricky_vars} CFGFLAGS"
    unset tmpflags_atompaw
  else
    AC_MSG_NOTICE([CFGFLAGS_ATOMPAW set => skipping AtomPAW config tricks])
    test "${afb_atompaw_tricks}" = "yes" && afb_atompaw_tricks="partial"
  fi
]) # AFB_TRICKS_ATOMPAW

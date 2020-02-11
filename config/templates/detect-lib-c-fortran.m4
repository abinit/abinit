# -*- Autoconf -*-
#
# Copyright (C) 2015 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# M4 template to detect C + Fortran libraries
#



# _ABI_CHECK_LIBS_@@NAME@@()
# -----@@dashes@@-------------
#
# Check whether the @@Name@@ library works.
#
AC_DEFUN([_ABI_CHECK_LIBS_@@NAME@@],[
  dnl Init
  abi_@@name@@_ok="unknown"
  abi_@@name@@_has_hdrs="unknown"
  abi_@@name@@_has_libs="unknown"
  abi_@@name@@_has_mods="unknown"
  abi_@@name@@_has_fort="unknown"
  abi_@@name@@_has_mpi="unknown"
  abi_@@name@@_has_omp="unknown"
  abi_@@name@@_has_thr="unknown"
  abi_@@name@@_fcflags="@@fcflags@@"
  abi_@@name@@_ldflags="@@ldflags@@"
  abi_@@name@@_incs="@@includes@@"
  abi_@@name@@_libs="@@libraries@@"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_@@name@@_incs}"
  FCFLAGS="${FCFLAGS} ${abi_@@name@@_incs}"
  LIBS="${abi_@@name@@_libs} ${LIBS}"

  dnl Look for C includes
  LDFLAGS="${CC_LDFLAGS}"
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([@@headers@@],
    [abi_@@name@@_has_hdrs="yes"], [abi_@@name@@_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for C libraries and routines
  if test "${abi_@@name@@_has_hdrs}" = "yes"; then
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    if test "@@libraries@@" = ""; then
      @@searchclibs@@
    else
      @@checkclibs@@
    fi
    AC_LANG_POP([C])
  fi

  dnl Look for Fortran libraries and routines
  if test "${abi_@@name@@_has_libs}" = "c"; then
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    if test "${with_@@name@@_libs}" = ""; then
      @@searchflibs@@
    else
      @@checkflibs@@
    fi
    AC_LANG_POP([Fortran])
  fi

  dnl Look for Fortran includes
  dnl Note: must be done after the libraries have been discovered
  @@checkfincs@@

  dnl Look for Fortran modules
  dnl Note: must be done after the libraries have been discovered
  @@checkfmods@@

  dnl Check Fortran support
  if test "${abi_@@name@@_has_mods}" = "yes"; then
    AC_MSG_CHECKING([whether FFTW3 Fortran wrappers work])
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        @@checkfprog@@
      ]])], [abi_@@name@@_has_fort="yes"], [abi_@@name@@_has_fort="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_@@name@@_has_fort}])
  fi

  dnl Set status of minimum FFTW3 support
  @@setstatus@@

  dnl Check for MPI support
  @@checkmpi@@

  dnl Check for OpenMP support
  @@checkomp@@

  dnl Check for POSIX threads support
  @@checkpth@@

  dnl Restore environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE
]) # _ABI_CHECK_LIBS_@@NAME@@

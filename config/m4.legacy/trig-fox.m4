# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for the FoX library
#



# ABI_TRIGGER_FOX()
# -----------------
#
# Check whether the FoX library is working.
#
AC_DEFUN([ABI_TRIGGER_FOX],[
  dnl Init
  abi_fox_default_libs="-lFoX_sax -lFoX_utils -lFoX_fsys -lFoX_common"
  abi_fox_has_incs="no"
  abi_fox_has_libs="no"
  abi_fox_serial="no"
  abi_fox_mpi="no"
  abi_fox_fcflags=""
  abi_fox_ldflags=""
  abi_fox_incs="${with_fox_incs}"
  abi_fox_libs="${with_fox_libs}"
  abi_test_fox="no"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  LDFLAGS="${FC_LDFLAGS}"
  AC_LANG_PUSH([Fortran])
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  FCFLAGS="${FCFLAGS} ${abi_fox_incs}"
  if test "${abi_fox_libs}" = ""; then
    AC_MSG_CHECKING([for FoX libraries to try])
    LIBS="${abi_fox_default_libs} ${LIBS}"
    AC_MSG_RESULT([${abi_fox_default_libs}])
  else
    LIBS="${abi_fox_libs} ${LIBS}"
  fi

  dnl Look for includes
  ABI_FC_MOD_INCS([fox_sax])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" != "unknown"; then
    abi_fox_has_incs="yes"
  fi

  dnl Look for libraries and routines
  if test "${abi_fox_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether FoX libraries work])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use fox_sax
        type(xml_t) :: xt
        call open_xml_file(xt,"conftest.xml")
      ]])], [abi_fox_has_libs="yes"], [abi_fox_has_libs="no"])
    AC_MSG_RESULT([${abi_fox_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_fox_has_incs}" = "yes" -a \
          "${abi_fox_has_libs}" = "yes"; then
    abi_fox_serial="yes"
    if test "${with_fox_libs}" = ""; then
      abi_fox_libs="${abi_fox_default_libs}"
    fi
  fi
  if test "${abi_fox_serial}" = "yes" -o \
          "${enable_fallbacks}" = "yes"; then
    AC_DEFINE([HAVE_FOX],1,
      [Define to 1 if you have the FoX library.])
    abi_test_fox="yes"
  fi
  if test "${abi_fox_serial}" = "yes"; then
    abi_fox_incs="${abi_fox_incs}"
    abi_fox_libs="${abi_fox_libs}"
  else
    if test "${with_fox_libs}" = ""; then
      AC_MSG_WARN([falling back to internal FoX version])
      abi_fallbacks="${abi_fallbacks} fox"
    else
      ABI_MSG_NOTICE([connectors-failure],[FoX detection failure])
      AC_MSG_ERROR([external FoX support does not work])
    fi
  fi

  dnl Check for MPI support
  if test "${enable_mpi_io}" = "yes" -a \
          "${abi_fox_serial}" = "yes"; then
    abi_fox_mpi="yes"
  fi

  dnl Restore environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the library
  AC_SUBST(abi_fox_fcflags)
  AC_SUBST(abi_fox_incs)
  AC_SUBST(abi_fox_ldflags)
  AC_SUBST(abi_fox_libs)
]) # ABI_TRIGGER_FOX

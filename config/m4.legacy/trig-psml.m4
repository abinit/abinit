# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for Density-Functional Theory libraries
#



# ABI_TRIGGER_PSML()
# -----------------------
#
# Check whether the PSML library is working.
#
AC_DEFUN([ABI_TRIGGER_PSML],[
  dnl Init
  abi_psml_has_incs="no"
  abi_psml_has_libs="no"
  abi_psml_serial="no"
  abi_psml_mpi="no"
  abi_psml_fcflags=""
  abi_psml_ldflags=""
  abi_psml_incs="${with_psml_incs}"
  abi_psml_libs="${with_psml_libs}"

  dnl PSML requires linear algebra support
  if test "${abi_linalg_serial}" = "yes"; then
    abi_psml_prereqs="yes"
  else
    AC_MSG_WARN([psml requires missing linear algebra support])
    if test "${abi_dft_linalg_fallback}" != "yes"; then
      abi_psml_prereqs="no"
    fi
  fi

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${abi_psml_libs} ${abi_linalg_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Look for libraries and routines
  if test "${abi_psml_libs}" = ""; then
    AC_SEARCH_LIBS([wannier_run],[wannier psml],
      [abi_psml_has_incs="yes"; abi_psml_has_libs="yes"],
      [abi_psml_has_libs="no"])
    if test "${abi_psml_has_libs}" = "yes"; then
      if test "${ac_cv_search_wannier_run}" != "none required"; then
        abi_psml_libs="${ac_cv_search_wannier_run}"
      fi
    fi
  else
    AC_MSG_CHECKING([whether the specified PSML library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call wannier_run
      ]])],
      [abi_psml_has_incs="yes"; abi_psml_has_libs="yes"],
      [abi_psml_has_libs="no"])
    AC_MSG_RESULT([${abi_psml_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_psml_has_incs}" = "yes" -a \
          "${abi_psml_has_libs}" = "yes"; then
    abi_psml_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes"; then
    if test "${abi_wannier_serial}" = "yes"; then
      abi_wannier_mpi="yes"
    fi
  fi


  if test "${abi_psml_serial}" = "yes" -o \
          "${enable_fallbacks}" = "yes"; then
    AC_DEFINE([HAVE_PSML],1,
      [Define to 1 if you have the PSML library.])
  fi
  if test "${abi_psml_serial}" = "yes"; then
    abi_psml_incs="${abi_psml_incs}"
    abi_psml_libs="${abi_psml_libs}"
  fi

  dnl Rebuild actual flavor
  if test "${abi_psml_fallback}" = "yes"; then
    abi_fallbacks="${abi_fallbacks} ${abi_psml_base_flavor}"
    if test "${abi_psml_prereqs}" = "no"; then
      ABI_MSG_NOTICE([connectors-failure],[PSML detection failure])
      AC_MSG_ERROR([prerequisites for PSML not found])
    fi
  else
    if test "${abi_psml_serial}" = "no"; then
      if test "${abi_psml_libs}" = "" -a \
              "${abi_psml_prereqs}" != "no"; then
        AC_MSG_WARN([falling back to internal ${abi_dft_flavor} version])
        abi_fallbacks="${abi_fallbacks} ${abi_dft_flavor}"
      else
        ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
        AC_MSG_ERROR([external PSML support does not work])
      fi
    fi
  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the libraries
  AC_SUBST(abi_psml_fcflags)
  AC_SUBST(abi_psml_incs)
  AC_SUBST(abi_psml_ldflags)
  AC_SUBST(abi_psml_libs)
]) # ABI_TRIGGER_PSML

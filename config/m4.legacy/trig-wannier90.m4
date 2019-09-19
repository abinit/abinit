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



# ABI_TRIGGER_WANNIER90()
# -----------------------
#
# Check whether the Wannier90 library is working.
#
AC_DEFUN([ABI_TRIGGER_WANNIER90],[
  dnl Init
  abi_wannier90_has_incs="no"
  abi_wannier90_has_libs="no"
  abi_wannier90_serial="no"
  abi_wannier90_mpi="no"
  abi_wannier90_fcflags=""
  abi_wannier90_ldflags=""
  abi_wannier90_incs="${with_wannier90_incs}"
  abi_wannier90_libs="${with_wannier90_libs}"

  dnl Wannier90 requires linear algebra support
  if test "${abi_linalg_serial}" = "yes"; then
    abi_wannier90_prereqs="yes"
  else
    AC_MSG_WARN([wannier90 requires missing linear algebra support])
    if test "${abi_dft_linalg_fallback}" != "yes"; then
      abi_wannier90_prereqs="no"
    fi
  fi

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${abi_wannier90_libs} ${abi_linalg_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Look for libraries and routines
  if test "${abi_wannier90_libs}" = ""; then
    AC_SEARCH_LIBS([wannier_run],[wannier wannier90],
      [abi_wannier90_has_incs="yes"; abi_wannier90_has_libs="yes"],
      [abi_wannier90_has_libs="no"])
    if test "${abi_wannier90_has_libs}" = "yes"; then
      if test "${ac_cv_search_wannier_run}" != "none required"; then
        abi_wannier90_libs="${ac_cv_search_wannier_run}"
      fi
    fi
  else
    AC_MSG_CHECKING([whether the specified Wannier90 library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call wannier_run
      ]])],
      [abi_wannier90_has_incs="yes"; abi_wannier90_has_libs="yes"],
      [abi_wannier90_has_libs="no"])
    AC_MSG_RESULT([${abi_wannier90_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_wannier90_has_incs}" = "yes" -a \
          "${abi_wannier90_has_libs}" = "yes"; then
    abi_wannier90_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes"; then
    if test "${abi_wannier_serial}" = "yes"; then
      abi_wannier_mpi="yes"
    fi
  fi


  if test "${abi_wannier90_serial}" = "yes" -o \
          "${enable_fallbacks}" = "yes"; then
    AC_DEFINE([HAVE_WANNIER90],1,
      [Define to 1 if you have the Wannier90 library.])
  fi
  if test "${abi_wannier90_serial}" = "yes"; then
    abi_wannier90_incs="${abi_wannier90_incs}"
    abi_wannier90_libs="${abi_wannier90_libs}"
  fi

  dnl Rebuild actual flavor
  if test "${abi_wannier90_fallback}" = "yes"; then
    abi_fallbacks="${abi_fallbacks} ${abi_wannier90_base_flavor}"
    if test "${abi_wannier90_prereqs}" = "no"; then
      ABI_MSG_NOTICE([connectors-failure],[Wannier90 detection failure])
      AC_MSG_ERROR([prerequisites for Wannier90 not found])
    fi
  else
    if test "${abi_wannier90_serial}" = "no"; then
      if test "${abi_wannier90_libs}" = "" -a \
              "${abi_wannier90_prereqs}" != "no"; then
        AC_MSG_WARN([falling back to internal ${abi_dft_flavor} version])
        abi_fallbacks="${abi_fallbacks} ${abi_dft_flavor}"
      else
        ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
        AC_MSG_ERROR([external Wannier90 support does not work])
      fi
    fi
  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the libraries
  AC_SUBST(abi_wannier90_fcflags)
  AC_SUBST(abi_wannier90_incs)
  AC_SUBST(abi_wannier90_ldflags)
  AC_SUBST(abi_wannier90_libs)
]) # ABI_TRIGGER_WANNIER90

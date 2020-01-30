# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for the BigDFT library
#



# ABI_BIGDFT_DETECT()
# --------------------
#
# Check whether the BigDFT library is working.
#
AC_DEFUN([ABI_BIGDFT_DETECT],[
  dnl Init
  abi_bigdft_has_mods="unknown"
  abi_bigdft_has_libs="unknown"
  abi_bigdft_ok="unknown"
  abi_bigdft_fcflags=""
  abi_bigdft_ldflags=""
  abi_bigdft_incs="${with_bigdft_incs}"
  abi_bigdft_libs="${with_bigdft_libs}"

  dnl Check for prerequisites
  if test "${abi_libxc_ok}" = "yes" -a \
          "${abi_linalg_ok}" = "yes"; then

    dnl Display input parameters
    AC_MSG_CHECKING([whether BigDFT includes have been specified])
    if test "${with_bigdft_incs}" = ""; then
      AC_MSG_RESULT([no])
    else
      AC_MSG_RESULT([yes])
    fi
    AC_MSG_CHECKING([whether BigDFT libraries have been specified])
    if test "${with_bigdft_libs}" = ""; then
      AC_MSG_RESULT([no])
    else
      AC_MSG_RESULT([yes])
    fi

    dnl Prepare environment
    ABI_ENV_BACKUP
    abi_saved_LIBS="${LIBS}"
    LDFLAGS="${FC_LDFLAGS}"
    FCFLAGS="${FCFLAGS} ${abi_bigdft_incs}"
    LIBS="${abi_bigdft_libs} ${abi_libxc_libs} ${abi_linalg_libs} ${LIBS}"

    dnl Look for libraries
    if test "${with_bigdft_libs}" = ""; then
      AC_LANG_PUSH([Fortran])
      AC_SEARCH_LIBS([bigdft1], [abinit],
        [abi_bigdft_has_libs="abinit"], [abi_bigdft_has_libs="no"])
      if test "${ac_cv_search_bigdft1}" != "no"; then
        if test "${ac_cv_search_bigdft1}" != "none required"; then
          abi_bigdft_libs="${ac_cv_search_bigdft1} ${abi_bigdft_libs}"
        fi
      fi
      if test "${abi_bigdft_has_libs}" = "abinit"; then
        AC_SEARCH_LIBS([bigdft2], [poissonsolver],
          [abi_bigdft_has_libs="poisson"], [abi_bigdft_has_libs="no"])
        if test "${ac_cv_search_bigdft2}" != "no"; then
          if test "${ac_cv_search_bigdft2}" != "none required"; then
            abi_bigdft_libs="${ac_cv_search_bigdft2} ${abi_bigdft_libs}"
          fi
        fi
      fi
      if test "${abi_bigdft_has_libs}" = "poisson"; then
        AC_SEARCH_LIBS([bigdft3], [bigdft],
          [abi_bigdft_has_libs="yes"], [abi_bigdft_has_libs="no"])
        if test "${ac_cv_search_bigdft3}" != "no"; then
          if test "${ac_cv_search_bigdft3}" != "none required"; then
            abi_bigdft_libs="${ac_cv_search_bigdft3} ${abi_bigdft_libs}"
          fi
        fi
      fi
      AC_LANG_POP([Fortran])
    else
      AC_LANG_PUSH([Fortran])
      AC_MSG_CHECKING([whether specified BigDFT Fortran libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          call bigdft1
          call bigdft2
          call bigdft3
        ]])], [abi_bigdft_has_libs="yes"], [abi_bigdft_has_libs="no"])
      AC_MSG_RESULT([${abi_bigdft_has_fort}])
      AC_LANG_POP([Fortran])
    fi

    dnl Look for Fortran modules
    ABI_FC_MOD_INCS([bigdft_api])
    FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
    if test "${abi_fc_mod_incs_ok}" != "unknown"; then
      abi_bigdft_has_mods="yes"
    fi

    dnl Look for libraries and routines
    AC_MSG_CHECKING([whether BigDFT libraries work])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use bigdft_api
        implicit none
        integer iproc
        type(input_variables) :: inputs
        type(atoms_data) :: atoms
        type(restart_objects) :: rst
        character(len=*),parameter :: routine = "conftest"
        call init_restart_objects(iproc,inputs,atoms,rst,routine)
      ]])], [abi_bigdft_has_fort="yes"], [abi_bigdft_has_fort="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_bigdft_has_fort}])

    dnl Take final decision for the serial case
    if test "${abi_bigdft_has_libs}" = "yes" -a \
            "${abi_bigdft_has_mods}" = "yes" -a \
            "${abi_bigdft_has_fort}" = "yes"; then
      abi_bigdft_ok="yes"
    fi

    dnl Check for MPI support
    if test "${enable_mpi}" = "yes"; then
      if test "${abi_bigdft_ok}" = "yes"; then
        abi_bigdft_mpi="yes"
      fi
    fi

    dnl Propagate information
    if test "${abi_bigdft_ok}" = "yes" -o \
            "${enable_fallbacks}" = "yes"; then
      AC_DEFINE([HAVE_BIGDFT],1,
        [Define to 1 if you have the BigDFT library.])
    fi

    if test "${abi_bigdft_fallback}" = "yes"; then
      if test "${abi_bigdft_prereqs}" = "no"; then
        ABI_MSG_NOTICE([connectors-failure],[BigDFT detection failure])
        AC_MSG_ERROR([prerequisites for BigDFT not found])
      fi
      abi_fallbacks="${abi_fallbacks} ${abi_bigdft_base_flavor}"
    else
      if test "${abi_bigdft_ok}" = "no"; then
        if test "${abi_bigdft_libs}" = "" -a "${abi_bigdft_prereqs}" != "no"; then
          AC_MSG_WARN([falling back to internal BigDFT version])
          abi_fallbacks="${abi_fallbacks} ${abi_dft_flavor}"
        else
          ABI_MSG_NOTICE([connectors-failure],[BigDFT detection failure])
          AC_MSG_ERROR([external BigDFT support does not work])
        fi
      fi
    fi

    dnl Restore build environment
    LIBS="${abi_saved_LIBS}"
    ABI_ENV_RESTORE

  else

    abi_bigdft_ok="no"

  fi

  dnl Substitute variables needed for the use of the libraries
  AC_SUBST(abi_bigdft_fcflags)
  AC_SUBST(abi_bigdft_ldflags)
  AC_SUBST(abi_bigdft_incs)
  AC_SUBST(abi_bigdft_libs)
]) # ABI_BIGDFT_DETECT

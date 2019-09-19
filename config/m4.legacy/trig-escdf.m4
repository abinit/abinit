# -*- Autoconf -*-
#
# Copyright (C) 2015 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for ESCDF I/O libraries
#



# ABI_TRIGGER_ESCDF()
# -------------------
#
# Check whether the ESCDF libraries are working.
#
AC_DEFUN([ABI_TRIGGER_ESCDF],[
  dnl Init
  abi_escdf_ok="unknown"
  abi_escdf_has_hdrs="unknown"
  abi_escdf_has_libs="unknown"
  abi_escdf_has_mods="unknown"
  abi_escdf_has_fort="unknown"
  abi_escdf_has_mpi="unknown"
  abi_escdf_fcflags=""
  abi_escdf_ldflags=""
  abi_escdf_incs="${with_hpio_incs}"
  abi_escdf_libs="${with_hpio_libs}"

  dnl No option set-up because ESCDF is a flavor of HPIO

  dnl Check for NetCDF presence
  if test "${abi_netcdf_ok}" = "yes"; then

    dnl Prepare environment
    ABI_ENV_BACKUP
    abi_saved_LIBS="${LIBS}"
    CPPFLAGS="${CPPFLAGS} ${abi_netcdf_incs} ${abi_escdf_incs}"
    FCFLAGS="${FCFLAGS} ${abi_netcdf_incs} ${abi_escdf_incs}"
    LIBS="${abi_escdf_libs} ${abi_netcdf_libs} ${LIBS}"

    dnl Look for C includes
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADERS([escdf.h],
      [abi_escdf_has_hdrs="yes"], [abi_escdf_has_hdrs="no"])
    AC_LANG_POP([C])

    dnl Look for C libraries and routines
    if test "${abi_escdf_has_hdrs}" = "yes"; then
      LDFLAGS="${CC_LDFLAGS}"
      AC_LANG_PUSH([C])
      if test "${with_escdf_libs}" = ""; then
        AC_SEARCH_LIBS([escdf_open], [escdf],
          [abi_escdf_has_libs="c"], [abi_escdf_has_libs="no"])
        if test "${ac_cv_search_dftc_open}" != "no"; then
          if test "${ac_cv_search_dftc_open}" != "none required"; then
            abi_escdf_libs="${ac_cv_search_dftc_open}"
          fi
        fi
      else
        AC_MSG_CHECKING([whether specified ESCDF C libraries work])
        AC_LINK_IFELSE([AC_LANG_PROGRAM(
          [[
#include <escdf.h>
          ]],
          [[
            int dcid;
            return escdf_open('conftest.dft', dcid);
          ]])], [abi_escdf_has_libs="c"], [abi_escdf_has_libs="no"])
        tmp_escdf_has_libs=`echo "${abi_escdf_has_libs}" | sed -e 's/^c$/yes/'`
        AC_MSG_RESULT([${tmp_escdf_has_libs}])
        unset tmp_escdf_has_libs
      fi
      AC_LANG_POP([C])
    fi

    dnl Look for Fortran libraries and routines
    if test "${abi_escdf_has_libs}" = "c"; then
      LDFLAGS="${FC_LDFLAGS}"
      AC_LANG_PUSH([Fortran])
      if test "${with_escdf_libs}" = ""; then
        AC_SEARCH_LIBS([dftf_open], [escdff],
          [abi_escdf_has_libs="yes"], [abi_escdf_has_libs="no"])
        if test "${ac_cv_search_dftf_open}" != "no"; then
          if test "${ac_cv_search_dftf_open}" != "none required"; then
            abi_escdf_libs="${ac_cv_search_dftf_open} ${abi_escdf_libs}"
          fi
        fi
      else
        AC_MSG_CHECKING([whether specified ESCDF Fortran libraries work])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([],
          [[
            call escdff_open
          ]])], [abi_escdf_has_libs="yes"], [abi_escdf_has_libs="no"])
        AC_MSG_RESULT([${abi_escdf_has_libs}])
      fi
      AC_LANG_POP([Fortran])
    fi

    dnl Look for Fortran includes
    dnl Note: must be done after the libraries have been discovered
    if test "${abi_escdf_has_libs}" = "yes"; then
      LDFLAGS="${FC_LDFLAGS}"
      ABI_FC_MOD_INCS([escdf])
      FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
      if test "${abi_fc_mod_incs_ok}" = "yes"; then
        abi_escdf_has_mods="yes"
        if test "${abi_escdf_incs}" = ""; then
          abi_escdf_incs="${fc_mod_incs}"
        fi
      fi
    fi

    dnl Check Fortran support
    if test "${abi_escdf_has_mods}" = "yes"; then
      AC_MSG_CHECKING([whether ESCDF Fortran wrappers work])
      LDFLAGS="${FC_LDFLAGS}"
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use escdf
          character(len=*), parameter :: path = "dummy"
          integer :: mode, ncerr, ncid
          ncerr = escdff_open(path,mode,ncid)
        ]])], [abi_escdf_has_fort="yes"], [abi_escdf_has_fort="no"])
      AC_LANG_POP([Fortran])
      AC_MSG_RESULT([${abi_escdf_has_fort}])
    fi

    dnl Set status of minimum ESCDF support
    if test "${abi_escdf_has_hdrs}" = "yes" -a \
            "${abi_escdf_has_libs}" = "yes" -a \
            "${abi_escdf_has_mods}" = "yes" -a \
            "${abi_escdf_has_fort}" = "yes"; then
      abi_escdf_ok="yes"
    else
      abi_escdf_ok="no"
    fi

    dnl Check for MPI support
    if test "${enable_mpi_io}" = "yes" -a \
            "${abi_escdf_ok}" = "yes"; then
      AC_MSG_CHECKING([whether ESCDF supports MPI I/O])
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use escdf
          character(len=*), parameter :: path = "dummy"
          integer :: cmode, comm, info, ncerr, ncid
          ncerr = escdff_open_par(path, cmode, comm, info, ncid)
        ]])], [abi_escdf_has_mpi="yes"], [abi_escdf_has_mpi="no"])
      AC_LANG_POP([Fortran])
      AC_MSG_RESULT([${abi_escdf_has_mpi}])
    fi

    dnl Restore environment
    LIBS="${abi_saved_LIBS}"
    ABI_ENV_RESTORE

  else

    abi_escdf_ok="no"

  fi
]) # ABI_TRIGGER_ESCDF

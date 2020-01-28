# -*- Autoconf -*-
#
# Copyright (C) 2015 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for NetCDF I/O libraries
#



# ABI_TRIGGER_NETCDF()
# --------------------
#
# Check whether the NetCDF libraries are working.
#
AC_DEFUN([ABI_TRIGGER_NETCDF],[
  dnl Init
  abi_netcdf_ok="unknown"
  abi_netcdf_has_hdrs="unknown"
  abi_netcdf_has_libs="unknown"
  abi_netcdf_has_mods="unknown"
  abi_netcdf_has_fort="unknown"
  abi_netcdf_has_mpi="unknown"
  abi_netcdf_fcflags=""
  abi_netcdf_ldflags=""
  abi_netcdf_incs="${with_netcdf_incs}"
  abi_netcdf_libs="${with_netcdf_libs}"

  dnl No option set-up because NetCDF is a flavor of HPIO

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_netcdf_incs}"
  FCFLAGS="${FCFLAGS} ${abi_netcdf_incs}"
  LIBS="${abi_netcdf_libs} ${LIBS}"

  dnl Look for C includes
  LDFLAGS="${CC_LDFLAGS}"
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([netcdf.h],
    [abi_netcdf_has_hdrs="yes"], [abi_netcdf_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for C libraries and routines
  if test "${abi_netcdf_has_hdrs}" = "yes"; then
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    if test "${with_netcdf_libs}" = ""; then
      AC_SEARCH_LIBS([nc_open], [netcdf],
        [abi_netcdf_has_libs="c"], [abi_netcdf_has_libs="no"])
      if test "${ac_cv_search_nc_open}" != "no"; then
        if test "${ac_cv_search_nc_open}" != "none required"; then
          abi_netcdf_libs="${ac_cv_search_nc_open}"
        fi
      fi
    else
      AC_MSG_CHECKING([whether specified NetCDF C libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [[
#include <netcdf.h>
        ]],
        [[
          int ncid;
          return nc_open('conftest.nc', NC_WRITE, ncid);
        ]])], [abi_netcdf_has_libs="c"], [abi_netcdf_has_libs="no"])
      tmp_netcdf_has_libs=`echo "${abi_netcdf_has_libs}" | sed -e 's/^c$/yes/'`
      AC_MSG_RESULT([${tmp_netcdf_has_libs}])
      unset tmp_netcdf_has_libs
    fi
    AC_LANG_POP([C])
  fi

  dnl Look for Fortran libraries and routines
  if test "${abi_netcdf_has_libs}" = "c"; then
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    if test "${with_netcdf_libs}" = ""; then
      AC_SEARCH_LIBS([nf_open], [netcdff],
        [abi_netcdf_has_libs="yes"], [abi_netcdf_has_libs="no"])
      if test "${ac_cv_search_nf_open}" != "no"; then
        if test "${ac_cv_search_nf_open}" != "none required"; then
          abi_netcdf_libs="${ac_cv_search_nf_open} ${abi_netcdf_libs}"
        fi
      fi
    else
      AC_MSG_CHECKING([whether specified NetCDF Fortran libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          call nf_open
        ]])], [abi_netcdf_has_libs="yes"], [abi_netcdf_has_libs="no"])
      AC_MSG_RESULT([${abi_netcdf_has_libs}])
    fi
    AC_LANG_POP([Fortran])
  fi

  dnl Look for Fortran includes
  dnl Note: must be done after the libraries have been discovered
  if test "${abi_netcdf_has_libs}" = "yes"; then
    LDFLAGS="${FC_LDFLAGS}"
    ABI_FC_MOD_INCS([netcdf])
    FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
    if test "${abi_fc_mod_incs_ok}" = "yes"; then
      abi_netcdf_has_mods="yes"
      if test "${abi_netcdf_incs}" = ""; then
        abi_netcdf_incs="${fc_mod_incs}"
      fi
    fi
  fi

  dnl Check Fortran support
  if test "${abi_netcdf_has_mods}" = "yes"; then
    AC_MSG_CHECKING([whether NetCDF Fortran wrappers work])
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use netcdf
        character(len=*), parameter :: path = "dummy"
        integer :: mode, ncerr, ncid
        ncerr = nf90_open(path,mode,ncid)
      ]])], [abi_netcdf_has_fort="yes"], [abi_netcdf_has_fort="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_netcdf_has_fort}])
  fi

  dnl Set status of minimum NetCDF support
  if test "${abi_netcdf_has_hdrs}" = "yes" -a \
          "${abi_netcdf_has_libs}" = "yes" -a \
          "${abi_netcdf_has_mods}" = "yes" -a \
          "${abi_netcdf_has_fort}" = "yes"; then
    abi_netcdf_ok="yes"
  else
    abi_netcdf_ok="no"
  fi

  dnl Check for MPI support
  if test "${enable_mpi_io}" = "yes" -a \
          "${abi_netcdf_ok}" = "yes"; then
    AC_MSG_CHECKING([whether NetCDF supports MPI I/O])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use netcdf
        character(len=*), parameter :: path = "dummy"
        integer :: cmode, comm, info, ncerr, ncid
        ncerr = nf90_open_par(path, cmode, comm, info, ncid)
      ]])], [abi_netcdf_has_mpi="yes"], [abi_netcdf_has_mpi="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_netcdf_has_mpi}])
  fi

  dnl Restore environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables
  AC_SUBST(abi_netcdf_fcflags)
  AC_SUBST(abi_netcdf_ldflags)
  AC_SUBST(abi_netcdf_incs)
  AC_SUBST(abi_netcdf_libs)
]) # ABI_TRIGGER_NETCDF

# -*- Autoconf -*-
#
# Copyright (C) 2015 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for ETSF I/O libraries
#



# ABI_TRIGGER_ETSF_IO()
# ---------------------
#
# Check whether the ETSF_IO libraries are working.
#
AC_DEFUN([ABI_TRIGGER_ETSF_IO],[
  dnl Init
  abi_etsf_io_default_libs="-letsf_io_utils -letsf_io -letsf_io_low_level"
  abi_etsf_io_ok="unknown"
  abi_etsf_io_has_libs="unknown"
  abi_etsf_io_has_mods="unknown"
  abi_etsf_io_has_mpi="no"
  abi_etsf_io_fcflags=""
  abi_etsf_io_ldflags=""
  abi_etsf_io_incs="${with_hpio_incs}"
  abi_etsf_io_libs="${with_hpio_libs}"

  dnl No option set-up because ETSF_IO is a flavor of HPIO

  dnl Check for NetCDF presence
  if test "${abi_netcdf_ok}" = "yes"; then

    dnl Prepare environment
    ABI_ENV_BACKUP
    abi_saved_LIBS="${abi_netcdf_libs} ${LIBS}"
    LDFLAGS="${FC_LDFLAGS}"
    FCFLAGS="${FCFLAGS} ${abi_netcdf_incs} ${abi_etsf_io_incs}"
    if test "${abi_etsf_io_libs}" = ""; then
      AC_MSG_CHECKING([for ETSF_IO libraries to try])
      abi_etsf_io_libs="${abi_etsf_io_default_libs}"
      AC_MSG_RESULT([${abi_etsf_io_libs}])
    fi
    LIBS="${abi_etsf_io_libs} ${abi_netcdf_libs} ${LIBS}"

    dnl Look for includes
    ABI_FC_MOD_INCS([etsf_io])
    FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
    if test "${abi_fc_mod_incs_ok}" != "unknown"; then
      abi_etsf_io_has_mods="yes"
      if test "${abi_etsf_io_incs}" = ""; then
        abi_etsf_io_incs="${fc_mod_incs}"
      fi
    fi

    dnl Look for libraries and routines
    if test "${abi_etsf_io_has_mods}" = "yes"; then
      AC_MSG_CHECKING([whether ETSF_IO libraries work])
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use etsf_io_low_level
          use etsf_io
          use etsf_io_tools
          character(len=etsf_charlen),allocatable :: atoms(:)
          integer :: ncid
          logical :: lstat
          type(etsf_io_low_error) :: err
          call etsf_io_tools_get_atom_names(ncid,atoms,lstat,err)
        ]])], [abi_etsf_io_has_libs="yes"], [abi_etsf_io_has_libs="no"])
      AC_LANG_POP([Fortran])
      AC_MSG_RESULT([${abi_etsf_io_has_libs}])
    fi

    dnl Set status of minimum ETSF_IO support
    if test "${abi_etsf_io_has_libs}" = "yes" -a \
            "${abi_etsf_io_has_mods}" = "yes"; then
      abi_etsf_io_ok="yes"
    fi

    dnl Restore environment
    LIBS="${abi_saved_LIBS}"
    ABI_ENV_RESTORE

  else

    abi_etsf_io_ok="no"

  fi

  dnl Substitute variables
  AC_SUBST(abi_etsf_io_fcflags)
  AC_SUBST(abi_etsf_io_ldflags)
  AC_SUBST(abi_etsf_io_incs)
  AC_SUBST(abi_etsf_io_libs)
]) # ABI_TRIGGER_ETSF_IO

# -*- Autoconf -*-
#
# Copyright (C) 2015 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for high-performance I/O libraries
#



# ABI_HPIO_DETECT()
# ---------------
#
# Sets all variables needed to take benefit from high-performance I/O libraries.
#
AC_DEFUN([ABI_HPIO_DETECT],[
  dnl Initial setup
  abi_hpio_has_escdf="no"
  abi_hpio_has_netcdf="no"
  abi_hpio_flavor="${with_hpio_flavor}"
  abi_hpio_fcflags=""
  abi_hpio_ldflags=""
  abi_hpio_incs="${HPIO_CPPFLAGS}"
  abi_hpio_libs="${HPIO_LIBS}"
  abi_netcdf_fallback="unknown"

  dnl Display input parameters
  AC_MSG_CHECKING([for requested HPIO flavor])
  AC_MSG_RESULT([${with_hpio_flavor}])
  AC_MSG_CHECKING([whether NetCDF includes have been specified])
  if test "${with_netcdf_incs}" = ""; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi
  AC_MSG_CHECKING([whether NetCDF libraries have been specified])
  if test "${with_netcdf_libs}" = ""; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi
  AC_MSG_CHECKING([whether HPIO includes have been specified])
  if test "${with_hpio_incs}" = ""; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi
  AC_MSG_CHECKING([whether HPIO libraries have been specified])
  if test "${with_hpio_libs}" = ""; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi

  dnl Check for NetCDF first
  ABI_TRIGGER_NETCDF

  dnl Propagate NetCDF information to source code
  if test "${abi_netcdf_ok}" = "yes" -o \
          "${enable_fallbacks}" = "yes"; then
    AC_DEFINE([HAVE_NETCDF], 1,
      [Define to 1 if you have the NetCDF library.])
    abi_hpio_incs="${abi_netcdf_incs}"
    abi_hpio_libs="${abi_netcdf_libs}"
    abi_netcdf_fallback="no"
  fi

  dnl Inform user about NetCDF status
  if test "${abi_netcdf_ok}" = "no"; then
    if test "${with_netcdf_incs}" = "" -a "${with_netcdf_libs}" = ""; then
      AC_MSG_WARN([falling back to limited developer NetCDF version])
      abi_netcdf_fallback="yes"
      abi_fallbacks="${abi_fallbacks} netcdf"
      abi_netcdf_fcflags=""
      abi_netcdf_ldflags=""
      abi_netcdf_incs=""
      abi_netcdf_libs=""
    else
      ABI_MSG_NOTICE([connectors-failure],[NetCDF detection failure])
      AC_MSG_ERROR([external NetCDF support does not work])
    fi
  fi
  AC_MSG_CHECKING([for final NetCDF Fortran flags])
  AC_MSG_RESULT(['${abi_netcdf_fcflags}'])
  AC_MSG_CHECKING([for final NetCDF link flags])
  AC_MSG_RESULT(['${abi_netcdf_ldflags}'])
  AC_MSG_CHECKING([for final NetCDF include flags])
  AC_MSG_RESULT(['${abi_netcdf_incs}'])
  AC_MSG_CHECKING([for final NetCDF library flags])
  AC_MSG_RESULT(['${abi_netcdf_libs}'])

  dnl Look for high-performance I/O libraries
  case "${abi_hpio_flavor}" in

    auto)
      if test "${abi_netcdf_fallback}" = "yes"; then
        abi_hpio_flavor="netcdf"
      else
        abi_hpio_flavor="escdf"
        ABI_TRIGGER_ESCDF
        if test "${abi_escdf_ok}" != "yes"; then
          abi_hpio_flavor="none"
        fi
      fi
      ;;

    escdf)
      if test "${abi_netcdf_fallback}" = "yes"; then
        AC_MSG_ERROR([the NetCDF fallback is too limited for ESCDF
                  please install a full recent version of NetCDF if you wish
                  to use ESCDF])
      else
        ABI_TRIGGER_ESCDF
      fi
      ;;

    none)
      AC_MSG_WARN([Fortran I/O is deprecated and will be disabled in the future])
      ;;

    *)
      AC_MSG_ERROR([unknown high-performance I/O flavor '${abi_hpio_flavor}'])
      ;;

  esac

  dnl Detection completed => propagate information
  case "${abi_hpio_flavor}" in

    escdf)
      if test "${abi_escdf_ok}" = "yes"; then
        AC_DEFINE([HAVE_ESCDF], 1,
          [Define to 1 if you have the ESCDF library.])
        abi_hpio_fcflags="${abi_hpio_fcflags} ${abi_escdf_fcflags}"
        abi_hpio_ldflags="${abi_hpio_ldflags} ${abi_escdf_ldflags}"
        abi_hpio_incs="${abi_hpio_incs} ${abi_escdf_incs}"
        abi_hpio_libs="${abi_escdf_libs} ${abi_hpio_libs}"
      else
        if test "${with_escdf_incs}" = "" -a "${with_escdf_libs}" = ""; then
          AC_MSG_WARN([falling back to deprecated Fortran I/O])
          abi_hpio_flavor="none"
          abi_escdf_fcflags=""
          abi_escdf_ldflags=""
          abi_escdf_incs=""
          abi_escdf_libs=""
        else
          ABI_MSG_NOTICE([connectors-failure],[ESCDF detection failure])
          AC_MSG_ERROR([external ESCDF support does not work])
        fi
      fi
      ;;

  esac

  dnl Propagate downgrade from ESCDF to Fortran I/O
  if test "${abi_hpio_flavor}" = "none"; then
    AC_DEFINE([HAVE_FORTRAN_IO],1,
      [Define to 1 if you have Fortran I/O only (deprecated).])
    abi_hpio_fcflags=""
    abi_hpio_ldflags=""
    abi_hpio_incs=""
    abi_hpio_libs=""
    abi_netcdf_fallback="no"
  fi

  dnl Final report
  AC_MSG_CHECKING([for final HPIO flavor])
  AC_MSG_RESULT(['${abi_hpio_flavor}'])
  AC_MSG_CHECKING([for final HPIO Fortran flags])
  AC_MSG_RESULT(['${abi_hpio_fcflags}'])
  AC_MSG_CHECKING([for final HPIO link flags])
  AC_MSG_RESULT(['${abi_hpio_ldflags}'])
  AC_MSG_CHECKING([for final HPIO include flags])
  AC_MSG_RESULT(['${abi_hpio_incs}'])
  AC_MSG_CHECKING([for final HPIO library flags])
  AC_MSG_RESULT(['${abi_hpio_libs}'])
  AC_MSG_CHECKING([whether we need a NetCDF fallback])
  AC_MSG_RESULT([${abi_netcdf_fallback}])
  if test "${abi_netcdf_fallback}" = "yes"; then
    AC_MSG_WARN([fallbacks have been requested
                  you should NOT run production calculations])
  fi

  dnl Substitute variables needed for the use of the libraries
  AC_SUBST(abi_hpio_flavor)
  AC_SUBST(abi_hpio_fcflags)
  AC_SUBST(abi_hpio_ldflags)
  AC_SUBST(abi_hpio_incs)
  AC_SUBST(abi_hpio_libs)
]) # ABI_HPIO_DETECT

# -*- Autoconf -*-
#
# Copyright (C) 2005-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for transferable file I/O libraries
#



# _ABI_TRIO_CHECK_ETSF_IO()
# -------------------------
#
# Check whether the ETSF_IO library is working.
#
AC_DEFUN([_ABI_TRIO_CHECK_ETSF_IO],[
  dnl Init
  abi_trio_etsf_io_default_libs="-letsf_io_utils -letsf_io"
  abi_trio_etsf_io_has_incs="no"
  abi_trio_etsf_io_has_libs="no"
  abi_trio_etsf_io_serial="no"
  abi_trio_etsf_io_mpi="no"
  abi_trio_etsf_io_fcflags=""
  abi_trio_etsf_io_ldflags=""
  abi_trio_etsf_io_incs="${with_etsf_io_incs}"
  abi_trio_etsf_io_libs="${with_etsf_io_libs}"

  dnl Prepare environment
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  FCFLAGS="${FCFLAGS} ${abi_trio_etsf_io_incs}"
  if test "${with_etsf_io_libs}" = ""; then
    AC_MSG_CHECKING([for ETSF_IO libraries to try])
    LIBS="${abi_trio_etsf_io_default_libs} ${LIBS}"
    AC_MSG_RESULT([${abi_trio_etsf_io_default_libs}])
  else
    LIBS="${abi_trio_etsf_io_libs} ${LIBS}"
  fi

  dnl Look for includes
  ABI_FC_MOD_INCS([etsf_io])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" != "unknown"; then
    abi_trio_etsf_io_has_incs="yes"
  fi

  dnl Look for libraries and routines
  if test "${abi_trio_etsf_io_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether ETSF_IO libraries work])
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
      ]])], [abi_trio_etsf_io_has_libs="yes"], [abi_trio_etsf_io_has_libs="no"])
    AC_MSG_RESULT([${abi_trio_etsf_io_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_trio_etsf_io_has_incs}" = "yes" -a \
          "${abi_trio_etsf_io_has_libs}" = "yes"; then
    abi_trio_etsf_io_serial="yes"
    if test "${with_etsf_io_libs}" = ""; then
      abi_trio_etsf_io_libs="${abi_trio_etsf_io_default_libs}"
    fi
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes" -a \
          "${abi_trio_etsf_io_serial}" = "yes"; then
    abi_trio_etsf_io_mpi="yes"
  fi

  dnl Restore environment
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_TRIO_CHECK_ETSF_IO



# _ABI_TRIO_CHECK_PSML()
# ----------------------
#
# Check whether the PSML library is working.
#
AC_DEFUN([_ABI_TRIO_CHECK_PSML],[
  dnl Init
  abi_trio_psml_default_libs="-lpsml"
  abi_trio_psml_has_incs="no"
  abi_trio_psml_has_libs="no"
  abi_trio_psml_serial="no"
  abi_trio_psml_mpi="no"
  abi_trio_psml_fcflags=""
  abi_trio_psml_ldflags=""
  abi_trio_psml_incs="${with_psml_incs}"
  abi_trio_psml_libs="${with_psml_libs}"

  dnl Prepare environment
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  FCFLAGS="${FCFLAGS} ${abi_trio_psml_incs}"
  if test "${abi_trio_psml_libs}" = ""; then
    AC_MSG_CHECKING([for PSML libraries to try])
    LIBS="${abi_trio_psml_default_libs} ${LIBS}"
    AC_MSG_RESULT([${abi_trio_psml_default_libs}])
  else
    LIBS="${abi_trio_psml_libs} ${LIBS}"
  fi

  dnl Look for includes
  ABI_FC_MOD_INCS([m_psml])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" != "unknown"; then
    abi_trio_psml_has_incs="yes"
  fi

  dnl Look for libraries and routines
  if test "${abi_trio_psml_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether PSML libraries work])

    AC_LINK_IFELSE( AC_LANG_SOURCE([[
        subroutine psml_die(str)
          character, intent(in) :: str 
          STOP 
        end subroutine psml_die 

        program test
        use m_psml
        type(ps_t) :: psxml
        call ps_destroy(psxml)
        end program
      ]]), [abi_trio_psml_has_libs="yes"], [abi_trio_psml_has_libs="no"])
    AC_MSG_RESULT([${abi_trio_psml_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_trio_psml_has_incs}" = "yes" -a \
          "${abi_trio_psml_has_libs}" = "yes"; then
    abi_trio_psml_serial="yes"
    if test "${with_psml_libs}" = ""; then
      abi_trio_psml_libs="${abi_trio_psml_default_libs}"
    fi
  fi

  dnl Check for MPI support
  if test "${enable_mpi_io}" = "yes" -a \
          "${abi_trio_psml_serial}" = "yes"; then
    abi_trio_psml_mpi="yes"
  fi

  dnl Restore environment
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_TRIO_CHECK_PSML



# _ABI_TRIO_CHECK_NETCDF()
# ------------------------
#
# Check whether the NetCDF library is working.
#
AC_DEFUN([_ABI_TRIO_CHECK_NETCDF],[
  dnl Init
  abi_trio_netcdf_has_incs="no"
  abi_trio_netcdf_has_libs="no"
  abi_trio_netcdf_serial="no"
  abi_trio_netcdf_mpi="no"
  abi_trio_netcdf_fcflags=""
  abi_trio_netcdf_ldflags=""
  abi_trio_netcdf_incs="${with_netcdf_incs}"
  abi_trio_netcdf_libs="${with_netcdf_libs}"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_trio_netcdf_incs}"
  FCFLAGS="${FCFLAGS} ${abi_trio_netcdf_incs}"

  dnl Look for C includes
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([netcdf.h],[abi_trio_netcdf_has_incs="yes"],[abi_trio_netcdf_has_incs="no"])
  AC_LANG_POP([C])

  dnl Look for libraries and routines
  if test "${abi_trio_netcdf_libs}" = ""; then
    AC_LANG_PUSH([C])
    AC_SEARCH_LIBS([nc_open],[netcdf])
    AC_LANG_POP([C])
    if test "${ac_cv_search_nc_open}" != "no"; then
      if test "${ac_cv_search_nc_open}" != "none required"; then
        abi_trio_netcdf_libs="${ac_cv_search_nc_open}"
        LIBS="${abi_trio_netcdf_libs} ${LIBS}"
      fi
      AC_SEARCH_LIBS([nf_open],[netcdf netcdff],
        [abi_trio_netcdf_has_libs="yes"],[abi_trio_netcdf_has_libs="no"])
      if test "${abi_trio_netcdf_has_libs}" = "yes"; then
        if test "${ac_cv_search_nf_open}" != "none required"; then
          abi_trio_netcdf_libs="${ac_cv_search_nf_open} ${abi_trio_netcdf_libs}"
        fi
      fi
    fi
  else
    LIBS="${abi_trio_netcdf_libs} ${LIBS}"
  fi

  dnl Look for Fortran includes
  dnl Note: must be done after the libraries have been discovered
  ABI_FC_MOD_INCS([netcdf])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" = "unknown"; then
    abi_trio_netcdf_has_incs="no"
  fi

  dnl Check Fortran support
  if test "${abi_trio_netcdf_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether NetCDF Fortran wrappers work])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use netcdf
        character(len=*), parameter :: path = "dummy"
        integer :: mode, ncerr, ncid
        ncerr = nf90_open(path,mode,ncid)
      ]])], [abi_trio_netcdf_has_libs="yes"], [abi_trio_netcdf_has_libs="no"])
    AC_MSG_RESULT([${abi_trio_netcdf_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_trio_netcdf_has_incs}" = "yes" -a \
          "${abi_trio_netcdf_has_libs}" = "yes"; then
    abi_trio_netcdf_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi_io}" = "yes" -a \
          "${abi_trio_netcdf_serial}" = "yes"; then
    AC_MSG_CHECKING([whether NetCDF supports MPI I/O])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use netcdf
        character(len=*), parameter :: path = "dummy"
        integer :: cmode, comm, info, ncerr, ncid
        ncerr = nf90_open_par(path, cmode, comm, info, ncid)
      ]])], [abi_trio_netcdf_mpi="yes"], [abi_trio_netcdf_mpi="no"])
    AC_MSG_RESULT([${abi_trio_netcdf_mpi}])
  fi

  dnl Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"

  dnl Make sure LIBS are properly set for ETSF_IO
  if test "${abi_trio_netcdf_serial}" = "yes"; then
    CPPFLAGS="${CPPFLAGS} ${abi_trio_netcdf_incs}"
    FCFLAGS="${FCFLAGS} ${abi_trio_netcdf_incs}"
    LIBS="${abi_trio_netcdf_libs} ${LIBS}"
  fi
]) # _ABI_TRIO_CHECK_NETCDF



# _ABI_TRIO_CHECK_YAML()
# ----------------------
#
# Check whether the YAML library is working.
#
AC_DEFUN([_ABI_TRIO_CHECK_YAML],[
  dnl Init
  abi_trio_yaml_default_libs="-lfyaml"
  abi_trio_yaml_has_incs="no"
  abi_trio_yaml_has_libs="no"
  abi_trio_yaml_serial="no"
  abi_trio_yaml_mpi="no"
  abi_trio_yaml_fcflags=""
  abi_trio_yaml_ldflags=""
  abi_trio_yaml_incs="${with_yaml_incs}"
  abi_trio_yaml_libs="${with_yaml_libs}"

  dnl Prepare environment
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  FCFLAGS="${FCFLAGS} ${abi_trio_yaml_incs}"
  if test "${with_yaml_libs}" = ""; then
    AC_MSG_CHECKING([for YAML libraries to try])
    LIBS="${abi_trio_yaml_default_libs} ${LIBS}"
    AC_MSG_RESULT([${abi_trio_yaml_default_libs}])
  else
    LIBS="${abi_trio_yaml_libs} ${LIBS}"
  fi

  dnl Look for includes
  ABI_FC_MOD_INCS([yaml_output])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" != "unknown"; then
    abi_trio_yaml_has_incs="yes"
  fi

  dnl Look for libraries and routines
  if test "${abi_trio_yaml_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether the YAML library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use yaml_output
        use dictionaries
        type(dictionary), pointer :: dict
        call yaml_new_document()
        call dict_init(dict)
        call set(dict//'tot_ncpus', 4)
      ]])], [abi_trio_yaml_has_libs="yes"], [abi_trio_yaml_has_libs="no"])
    AC_MSG_RESULT([${abi_trio_yaml_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_trio_yaml_has_incs}" = "yes" -a \
          "${abi_trio_yaml_has_libs}" = "yes"; then
    abi_trio_yaml_serial="yes"
    if test "${with_yaml_libs}" = ""; then
      abi_trio_yaml_libs="${abi_trio_yaml_default_libs}"
    fi
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes" -a \
          "${abi_trio_yaml_serial}" = "yes"; then
    abi_trio_yaml_mpi="yes"
  fi

  dnl Restore environment
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_TRIO_CHECK_YAML



# ABI_CONNECT_TRIO()
# ------------------
#
# Sets all variables needed to handle the transferable I/O libraries.
#
AC_DEFUN([ABI_CONNECT_TRIO],[
  dnl Initial setup
  abi_test_etsf_io="no"
  abi_test_netcdf="no"
  lib_trio_flavor="${with_trio_flavor}"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  LDFLAGS="${FC_LDFLAGS}"
  AC_LANG_PUSH([Fortran])

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested transferable I/O support])
  AC_MSG_RESULT([${with_trio_flavor}])

  dnl Look for external I/O libraries
  if test "${with_trio_flavor}" != "none"; then

    dnl Make sure NetCDF is looked for before ETSF_IO
    abi_trio_iter=`echo "${lib_trio_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s",[$]1}'`
    abi_trio_tmp="${abi_trio_iter}"
    for abi_trio_flavor in ${abi_trio_iter}; do
      if test "${abi_trio_flavor}" = "netcdf"; then
        abi_trio_tmp=`echo "${abi_trio_iter}" | sed -e 's/netcdf//'`
        abi_trio_tmp="netcdf ${abi_trio_tmp}"
      fi
      if test "${abi_trio_flavor}" = "netcdf-fallback"; then
        abi_trio_tmp=`echo "${abi_trio_iter}" | sed -e 's/netcdf-fallback//'`
        abi_trio_tmp="netcdf-fallback ${abi_trio_tmp}"
      fi
    done
    abi_trio_iter="${abi_trio_tmp}"

    for abi_trio_flavor in ${abi_trio_iter}; do

      dnl Check if the user has requested a fallback
      tmp_trio_base_flavor=`echo "${abi_trio_flavor}" | cut -d- -f1`
      AC_MSG_CHECKING([whether to select a fallback for ${tmp_trio_base_flavor}])
      tmp_trio_fallback=`echo "${abi_trio_flavor}" | cut -s -d- -f2`
      if test "${tmp_trio_fallback}" = "fallback"; then
        tmp_trio_fallback="yes"
      else
        tmp_trio_fallback="no"
      fi
      AC_MSG_RESULT([${tmp_trio_fallback}])
      if test "${tmp_trio_fallback}" = "yes" -a \
              "${enable_fallbacks}" = "no"; then
        AC_MSG_ERROR([fallback requested while fallbacks have been globally disabled])
      fi

      dnl Look for TRIO libraries
      case "${abi_trio_flavor}" in

        etsf_io)
          if test "${abi_trio_netcdf_serial}" = "yes"; then
            abi_trio_etsf_io_prereqs="yes"
            _ABI_TRIO_CHECK_ETSF_IO
          else
            AC_MSG_WARN([ETSF_IO requires missing NetCDF support])
            if test "${abi_trio_netcdf_fallback}" != "yes"; then
              abi_trio_etsf_io_prereqs="no"
            fi
            abi_trio_etsf_io_serial="no"
            abi_trio_etsf_io_mpi="no"
          fi
          if test "${abi_trio_etsf_io_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_ETSF_IO],1,
              [Define to 1 if you have the ETSF_IO library.])
            abi_test_etsf_io="yes"
          fi
          if test "${abi_trio_etsf_io_serial}" = "yes"; then
            lib_etsf_io_incs="${abi_trio_etsf_io_incs}"
            lib_etsf_io_libs="${abi_trio_etsf_io_libs}"
          fi
          ;;

        etsf_io-fallback)
          if test "${abi_trio_netcdf_serial}" != "yes" -a \
                  "${abi_trio_netcdf_fallback}" != "yes"; then
            AC_MSG_WARN([ETSF_IO requires missing NetCDF support])
            abi_trio_etsf_io_prereqs="no"
            abi_trio_etsf_io_serial="no"
            abi_trio_etsf_io_mpi="no"
          fi
          ;;

        psml)
          _ABI_TRIO_CHECK_PSML
          if test "${abi_trio_psml_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_PSML],1,
              [Define to 1 if you have the PSML library.])
            abi_test_psml="yes"
          fi
          if test "${abi_trio_psml_serial}" = "yes"; then
            lib_psml_incs="${abi_trio_psml_incs}"
            lib_psml_libs="${abi_trio_psml_libs}"
          else
            ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
            AC_ERROR([no working version of psml has been found])
          fi
          ;;

        netcdf)
          _ABI_TRIO_CHECK_NETCDF
          if test "${abi_trio_netcdf_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_NETCDF],1,
              [Define to 1 if you have the NetCDF library.])
            abi_test_netcdf="yes"
          fi
          if test "${abi_trio_netcdf_serial}" = "yes"; then
            lib_netcdf_incs="${abi_trio_netcdf_incs}"
            lib_netcdf_libs="${abi_trio_netcdf_libs}"
          elif test "${enable_fallbacks}" = "yes"; then
            abi_trio_netcdf_fallback="yes"
          fi

          if test "${abi_trio_netcdf_mpi}" = "yes"; then
            AC_DEFINE([HAVE_NETCDF_MPI],1,
              [Define to 1 if you have MPI-IO support in the NetCDF library.])
          fi
          ;;

        yaml)
          _ABI_TRIO_CHECK_YAML
          if test "${abi_trio_yaml_serial}" = "yes"; then
            AC_DEFINE([HAVE_YAML],1,
              [Define to 1 if you have the YAML library.])
            abi_test_yaml="yes"
            lib_yaml_incs="${abi_trio_yaml_incs}"
            lib_yaml_libs="${abi_trio_yaml_libs}"
          else
            ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
            AC_MSG_ERROR([the YAML library is absent or unusable])
          fi
          ;;

        *)
          if test "${tmp_trio_fallback}" = "no"; then
            AC_MSG_ERROR([unknown transferable I/O flavor '${abi_trio_flavor}'])
          fi
          ;;

      esac

      dnl Rebuild actual flavor
      if test "${tmp_trio_fallback}" = "yes"; then
        eval abi_trio_${tmp_trio_base_flavor}_fallback="yes"
        abi_fallbacks="${abi_fallbacks} ${tmp_trio_base_flavor}"
        tmp_trio_flavor="${tmp_trio_flavor}+${abi_trio_flavor}"
        tmp_trio_prereqs=`eval echo \$\{abi_trio_${tmp_trio_base_flavor}_prereqs\}`
        if test "${tmp_trio_prereqs}" = "no"; then
          ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
          AC_MSG_ERROR([prerequisites for ${abi_trio_flavor} not found])
        fi
      else
        tmp_trio_prereqs=`eval echo \$\{abi_trio_${abi_trio_flavor}_prereqs\}`
        tmp_trio_serial=`eval echo \$\{abi_trio_${abi_trio_flavor}_serial\}`
        tmp_trio_libs=`eval echo \$\{with_${abi_trio_flavor}_libs\}`
        if test "${tmp_trio_serial}" = "no"; then
          if test "${tmp_trio_libs}" = "" -a "${tmp_trio_prereqs}" != "no"; then
            AC_MSG_WARN([falling back to internal ${abi_trio_flavor} version])
            abi_fallbacks="${abi_fallbacks} ${abi_trio_flavor}"
            tmp_trio_flavor="${tmp_trio_flavor}+${abi_trio_flavor}-fallback"
          else
            ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
            AC_MSG_ERROR([external ${abi_trio_flavor} support does not work])
          fi
        else
          tmp_trio_flavor="${tmp_trio_flavor}+${abi_trio_flavor}"
        fi
      fi

    done

  fi

  dnl Restore environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  lib_trio_flavor=`echo "${tmp_trio_flavor}" | sed -e 's/^\+//;s/\+[$]//'`
  if test "${lib_trio_flavor}" = ""; then
    lib_trio_flavor="none"
  fi
  AC_MSG_CHECKING([for the actual transferable I/O support])
  AC_MSG_RESULT([${lib_trio_flavor}])

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_trio_flavor)
  AC_SUBST(lib_etsf_io_fcflags)
  AC_SUBST(lib_etsf_io_incs)
  AC_SUBST(lib_etsf_io_ldflags)
  AC_SUBST(lib_etsf_io_libs)
  AC_SUBST(lib_psml_fcflags)
  AC_SUBST(lib_psml_incs)
  AC_SUBST(lib_psml_ldflags)
  AC_SUBST(lib_psml_libs)
  AC_SUBST(lib_netcdf_fcflags)
  AC_SUBST(lib_netcdf_incs)
  AC_SUBST(lib_netcdf_ldflags)
  AC_SUBST(lib_netcdf_libs)
  AC_SUBST(lib_yaml_fcflags)
  AC_SUBST(lib_yaml_incs)
  AC_SUBST(lib_yaml_ldflags)
  AC_SUBST(lib_yaml_libs)
]) # ABI_CONNECT_TRIO

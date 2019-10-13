## Copyright (C) 2019 Yann Pouillon

#
# NetCDF I/O library - Fortran interface
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_NETCDF_FORTRAN_INIT], [
  # Init
  sd_netcdf_fortran_cppflags=""
  sd_netcdf_fortran_cflags=""
  sd_netcdf_fortran_cxxflags=""
  sd_netcdf_fortran_fcflags=""
  sd_netcdf_fortran_ldflags=""
  sd_netcdf_fortran_libs=""
  sd_netcdf_fortran_prefix=""
  sd_netcdf_fortran_enable=""
  sd_netcdf_fortran_init="unknown"
  sd_netcdf_fortran_mpi_ok="unknown"
  sd_netcdf_fortran_ok="unknown"

  # Set adjustable parameters
  sd_netcdf_fortran_options="$1"
  sd_netcdf_fortran_libs_def="$2"
  sd_netcdf_fortran_cppflags_def="$3"
  sd_netcdf_fortran_cflags_def="$4"
  sd_netcdf_fortran_cxxflags_def="$5"
  sd_netcdf_fortran_fcflags_def="$6"
  sd_netcdf_fortran_ldflags_def="$7"

  # Process options
  for kwd in ${sd_netcdf_fortran_options}; do
    case "${kwd}" in
      auto)
        sd_netcdf_fortran_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_netcdf_fortran_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_netcdf_fortran_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg NetCDF Fortran interface option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_netcdf_fortran_status}" && sd_netcdf_fortran_status="optional"
  test -z "${sd_netcdf_fortran_policy}" && sd_netcdf_fortran_policy="fail"
  test -z "${sd_netcdf_fortran_enable_def}" && sd_netcdf_fortran_enable_def="no"
  test -z "${sd_netcdf_fortran_libs_def}" && sd_netcdf_fortran_libs_def="-lnetcdff"
  case "${sd_netcdf_fortran_status}" in
    implicit|required)
      sd_netcdf_fortran_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([netcdf-fortran],
    [AS_HELP_STRING([--with-netcdf-fortran],
      [Install prefix of the NetCDF Fortran interface library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_netcdf_fortran_enable="${withval}"
        sd_netcdf_fortran_init="yon"
      else
        sd_netcdf_fortran_enable="yes"
        sd_netcdf_fortran_init="dir"
        sd_netcdf_fortran_prefix="${withval}"
      fi],
    [ sd_netcdf_fortran_enable="${sd_netcdf_fortran_enable_def}"; sd_netcdf_fortran_init="def"])

  # Declare environment variables
  AC_ARG_VAR([NETCDF_FORTRAN_CPPFLAGS], [C preprocessing flags for NetCDF Fortran interface.])
  AC_ARG_VAR([NETCDF_FORTRAN_CFLAGS], [C flags for NetCDF Fortran interface.])
  AC_ARG_VAR([NETCDF_FORTRAN_CXXFLAGS], [C++ flags for NetCDF Fortran interface.])
  AC_ARG_VAR([NETCDF_FORTRAN_FCFLAGS], [Fortran flags for NetCDF Fortran interface.])
  AC_ARG_VAR([NETCDF_FORTRAN_LDFLAGS], [Linker flags for NetCDF Fortran interface.])
  AC_ARG_VAR([NETCDF_FORTRAN_LIBS], [Library flags for NetCDF Fortran interface.])

  # Detect use of environment variables
  if test "${sd_netcdf_fortran_enable}" = "yes" -o "${sd_netcdf_fortran_enable}" = "auto"; then
    tmp_netcdf_fortran_vars="${NETCDF_FORTRAN_CPPFLAGS}${NETCDF_FORTRAN_CFLAGS}${NETCDF_FORTRAN_CXXFLAGS}${NETCDF_FORTRAN_FCFLAGS}${NETCDF_FORTRAN_LDFLAGS}${NETCDF_FORTRAN_LIBS}"
    if test "${sd_netcdf_fortran_init}" = "def" -a \
            ! -z "${tmp_netcdf_fortran_vars}"; then
      sd_netcdf_fortran_enable="yes"
      sd_netcdf_fortran_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_NETCDF_FORTRAN_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_netcdf_fortran_enable}" = "yes" -o "${sd_netcdf_fortran_enable}" = "auto"; then

    # Set NetCDF Fortran interface-specific flags
    case "${sd_netcdf_fortran_init}" in

      def|yon)
        sd_netcdf_fortran_cppflags="${sd_netcdf_fortran_cppflags_def}"
        sd_netcdf_fortran_cflags="${sd_netcdf_fortran_cflags_def}"
        sd_netcdf_fortran_cxxflags="${sd_netcdf_fortran_cxxflags_def}"
        sd_netcdf_fortran_fcflags="${sd_netcdf_fortran_fcflags_def}"
        sd_netcdf_fortran_ldflags="${sd_netcdf_fortran_ldflags_def}"
        sd_netcdf_fortran_libs="${sd_netcdf_fortran_libs_def}"
        ;;

      dir)
        sd_netcdf_fortran_cppflags="${sd_netcdf_fortran_cppflags_def} -I${sd_etcdf_fortran_prefix}/include"
        sd_netcdf_fortran_cflags="${sd_netcdf_fortran_cflags_def}"
        sd_netcdf_fortran_cxxflags="${sd_netcdf_fortran_cxxflags_def}"
        sd_netcdf_fortran_fcflags="${sd_netcdf_fortran_fcflags_def} -I${sd_netcdf_fortran_prefix}/include"
        sd_netcdf_fortran_ldflags="${sd_netcdf_fortran_ldflags_def}"
        sd_netcdf_fortran_libs="-L${sd_netcdf_fortran_prefix}/lib ${sd_netcdf_fortran_libs_def} ${sd_netcdf_fortran_libs}"
        ;;

      env)
        sd_netcdf_fortran_cppflags="${sd_netcdf_fortran_cppflags_def}"
        sd_netcdf_fortran_cflags="${sd_netcdf_fortran_cflags_def}"
        sd_netcdf_fortran_cxxflags="${sd_netcdf_fortran_cxxflags_def}"
        sd_netcdf_fortran_fcflags="${sd_netcdf_fortran_fcflags_def}"
        sd_netcdf_fortran_ldflags="${sd_netcdf_fortran_ldflags_def}"
        sd_netcdf_fortran_libs="${sd_netcdf_fortran_libs_def}"
        test ! -z "${NETCDF_FORTRAN_CPPFLAGS}" && sd_netcdf_fortran_cppflags="${NETCDF_FORTRAN_CPPFLAGS}"
        test ! -z "${NETCDF_FORTRAN_CFLAGS}" && sd_netcdf_fortran_cflags="${NETCDF_FORTRAN_CFLAGS}"
        test ! -z "${NETCDF_FORTRAN_CXXFLAGS}" && sd_netcdf_fortran_cxxflags="${NETCDF_FORTRAN_CXXFLAGS}"
        test ! -z "${NETCDF_FORTRAN_FCFLAGS}" && sd_netcdf_fortran_fcflags="${NETCDF_FORTRAN_FCFLAGS}"
        test ! -z "${NETCDF_FORTRAN_LDFLAGS}" && sd_netcdf_fortran_ldflags="${NETCDF_FORTRAN_LDFLAGS}"
        test ! -z "${NETCDF_FORTRAN_LIBS}" && sd_netcdf_fortran_libs="${NETCDF_FORTRAN_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for the NetCDF Fortran interface: '${sd_netcdf_fortran_init}'])
        ;;

    esac

  fi

  # Display configuration
  _SD_NETCDF_FORTRAN_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_netcdf_fortran_options)
  AC_SUBST(sd_netcdf_fortran_enable_def)
  AC_SUBST(sd_netcdf_fortran_policy)
  AC_SUBST(sd_netcdf_fortran_status)
  AC_SUBST(sd_netcdf_fortran_enable)
  AC_SUBST(sd_netcdf_fortran_init)
  AC_SUBST(sd_netcdf_fortran_ok)
  AC_SUBST(sd_netcdf_fortran_cppflags)
  AC_SUBST(sd_netcdf_fortran_cflags)
  AC_SUBST(sd_netcdf_fortran_fcflags)
  AC_SUBST(sd_netcdf_fortran_ldflags)
  AC_SUBST(sd_netcdf_fortran_libs)
  AC_SUBST(with_netcdf)

  # Clean-up
  unset tmp_netcdf_fortran_vars
]) # SD_NETCDF_FORTRAN_INIT


AC_DEFUN([SD_NETCDF_FORTRAN_DETECT], [
  # Display configuration
  _SD_NETCDF_FORTRAN_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_netcdf_fortran_enable}" = "auto" -o "${sd_netcdf_fortran_enable}" = "yes"; then
    _SD_NETCDF_FORTRAN_CHECK_USE

    if test "${sd_netcdf_fortran_ok}" = "yes"; then
      CPPFLAGS="${CPPFLAGS} ${sd_netcdf_fortran_cppflags}"
      CFLAGS="${CFLAGS} ${sd_netcdf_fortran_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_netcdf_fortran_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_netcdf_fortran_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_netcdf_fortran_ldflags}"
      LIBS="${sd_netcdf_fortran_libs} ${LIBS}"

      AC_DEFINE([HAVE_NETCDF_FORTRAN], 1,
        [Define to 1 if you have the NetCDF Fortran interface library.])

      if test "${sd_mpi_ok}" = "yes" -a "${sd_netcdf_fortran_mpi_ok}" = "yes"; then
        AC_DEFINE([HAVE_NETCDF_FORTRAN_MPI], 1,
          [Define to 1 if you have a parallel NetCDF Fortran interface library.])
      fi
    else
      if test "${sd_netcdf_fortran_status}" = "optional" -a \
              "${sd_netcdf_fortran_init}" = "def"; then
        sd_netcdf_fortran_enable="no"
        sd_netcdf_fortran_cppflags=""
        sd_netcdf_fortran_cflags=""
        sd_netcdf_fortran_cxxflags=""
        sd_netcdf_fortran_fcflags=""
        sd_netcdf_fortran_ldflags=""
        sd_netcdf_fortran_libs=""
      else
        AC_MSG_FAILURE([invalid NetCDF Fortran interface configuration])
      fi
    fi
  else
    sd_netcdf_fortran_enable="no"
    sd_netcdf_fortran_cppflags=""
    sd_netcdf_fortran_cflags=""
    sd_netcdf_fortran_cxxflags=""
    sd_netcdf_fortran_fcflags=""
    sd_netcdf_fortran_ldflags=""
    sd_netcdf_fortran_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_NETCDF_FORTRAN_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_netcdf_fortran_cppflags}"
  CFLAGS="${CFLAGS} ${sd_netcdf_fortran_cflags}"
  CXXFLAGS="${CXXFLAGS} ${sd_netcdf_fortran_cxxflags}"
  FCFLAGS="${FCFLAGS} ${sd_netcdf_fortran_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_netcdf_fortran_ldflags}"
  LIBS="${sd_netcdf_fortran_libs} ${LIBS}"

  # Check NetCDF Fortran API
  AC_MSG_CHECKING([whether the NetCDF Fortran interface works])
  if test "${sd_netcdf_ok}" = "yes"; then
    for tmp_incs in "" "-I/usr/include"; do
      FCFLAGS="${FCFLAGS} ${tmp_incs}"
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use netcdf
          character(len=*), parameter :: path = "dummy"
          integer :: mode, ncerr, ncid
          ncerr = nf90_open(path,mode,ncid)
        ]])], [sd_netcdf_fortran_ok="yes"], [sd_netcdf_fortran_ok="no"])
      AC_LANG_POP([Fortran])
      if test "${sd_netcdf_fortran_ok}" = "yes"; then
        test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
        break
      fi
    done
  else
    sd_netcdf_fortran_ok="no"
  fi
  AC_MSG_RESULT([${sd_netcdf_fortran_ok}])
  unset tmp_incs

  # Check if we can do parallel I/O
  if test "${sd_netcdf_fortran_ok}" = "yes" -a "${sd_mpi_ok}" = "yes"; then
    if test "${sd_hdf5_mpi_ok}" = "yes"; then
      AC_MSG_CHECKING([whether the NetCDF Fortran interface has parallel I/O])
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use mpi
          use netcdf
          integer :: ierr, ncid
          ierr = nf90_create("conftest.nc", ior(NF90_NETCDF4, NF90_MPIPOSIX), &
            ncid, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
        ]])], [sd_netcdf_fortran_mpi_ok="yes"], [sd_netcdf_fortran_mpi_ok="no"])
      AC_LANG_POP([Fortran])
      AC_MSG_RESULT([${sd_netcdf_fortran_mpi_ok}])
    else
      sd_netcdf_fortran_mpi_ok="no"
    fi
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_NETCDF_FORTRAN_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_NETCDF_FORTRAN_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_netcdf_fortran_invalid="no"
  if test "${sd_netcdf_fortran_enable_def}" != "auto" -a \
          "${sd_netcdf_fortran_enable_def}" != "no" -a \
          "${sd_netcdf_fortran_enable_def}" != "yes"; then
    case "${sd_netcdf_fortran_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_netcdf_fortran_enable_def = '${sd_netcdf_fortran_enable_def}'])
        ;;
      skip)
        tmp_netcdf_fortran_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_netcdf_fortran_enable_def = '${sd_netcdf_fortran_enable_def}'])
        tmp_netcdf_fortran_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_netcdf_fortran_invalid}" = "yes"; then
    if test "${sd_netcdf_fortran_status}" = "required"; then
      sd_netcdf_fortran_enable_def="yes"
    else
      sd_netcdf_fortran_enable_def="no"
    fi
    tmp_netcdf_fortran_invalid="no"
    AC_MSG_NOTICE([setting sd_netcdf_fortran_enable_def to '${sd_netcdf_fortran_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_netcdf_fortran_invalid="no"
  if test "${sd_netcdf_fortran_status}" = "implicit" -o \
          "${sd_netcdf_fortran_status}" = "required"; then
    if test "${sd_netcdf_fortran_enable}" = "no"; then
      case "${sd_netcdf_fortran_policy}" in
        fail)
          AC_MSG_ERROR([The NetCDF Fortran interface package is required and cannot be disabled
                  See https://www.unidata.ucar.edu/software/netcdf/ for
                  details on how to install it.])
          ;;
        skip)
          tmp_netcdf_fortran_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The NetCDF Fortran interface package is required and cannot be disabled])
          tmp_netcdf_fortran_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_netcdf_fortran_enable}" = "auto"; then
      AC_MSG_NOTICE([setting NetCDF Fortran interface trigger to yes])
      sd_netcdf_fortran_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_netcdf_fortran_invalid}" = "yes"; then
    case "${sd_netcdf_fortran_status}" in
      implicit|required)
        sd_netcdf_fortran_enable="yes"
        ;;
      optional)
        sd_netcdf_fortran_enable="no"
        ;;
    esac
    tmp_netcdf_fortran_invalid="no"
    AC_MSG_NOTICE([setting sd_netcdf_fortran_enable to '${sd_netcdf_fortran_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_netcdf_fortran_vars="${NETCDF_FORTRAN_CPPFLAGS}${NETCDF_FORTRAN_CFLAGS}${NETCDF_FORTRAN_FCFLAGS}${NETCDF_FORTRAN_LDFLAGS}${NETCDF_FORTRAN_LIBS}"
  tmp_netcdf_fortran_invalid="no"
  if test ! -z "${tmp_netcdf_fortran_vars}" -a ! -z "${sd_netcdf_fortran_prefix}"; then
    case "${sd_netcdf_fortran_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for NetCDF Fortran interface
                  Please use NETCDF_FORTRAN_FCFLAGS + NETCDF_FORTRAN_LIBS or --with-netcdf,
                  not both.])
        ;;
      skip)
        tmp_netcdf_fortran_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for NetCDF Fortran interface])
        tmp_netcdf_fortran_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_netcdf_fortran_vars}"; then
    sd_netcdf_fortran_enable="yes"
    sd_netcdf_fortran_init="env"
    if test "${tmp_netcdf_fortran_invalid}" = "yes"; then
      tmp_netcdf_fortran_invalid="no"
      AC_MSG_NOTICE([overriding --with-netcdf with NETCDF_FORTRAN_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_netcdf_fortran_status}" = "implicit"; then
    if test "${sd_netcdf_fortran_fcflags}" != ""; then
      sd_netcdf_fortran_fcflags=""
      AC_MSG_NOTICE([resetting NetCDF Fortran interface Fortran flags (implicit package)])
    fi
    if test "${sd_netcdf_fortran_ldflags}" != ""; then
      sd_netcdf_fortran_ldflags=""
      AC_MSG_NOTICE([resetting NetCDF Fortran interface linker flags (implicit package)])
    fi
    if test "${sd_netcdf_fortran_libs}" != ""; then
      sd_netcdf_fortran_libs=""
      AC_MSG_NOTICE([resetting NetCDF Fortran interface library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_netcdf_fortran_enable}" = "implicit"; then
    sd_netcdf_fortran_fcflags=""
    sd_netcdf_fortran_ldflags=""
    sd_netcdf_fortran_libs=""
  fi

  # Clean-up
  unset tmp_netcdf_fortran_invalid
  unset tmp_netcdf_fortran_vars
]) # _SD_NETCDF_FORTRAN_CHECK_CONFIG


AC_DEFUN([_SD_NETCDF_FORTRAN_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable the NetCDF Fortran interface])
  AC_MSG_RESULT([${sd_netcdf_fortran_enable}])
  if test "${sd_netcdf_fortran_enable}" != "no"; then
    AC_MSG_CHECKING([how the NetCDF Fortran interface parameters have been set])
    AC_MSG_RESULT([${sd_netcdf_fortran_init}])
    AC_MSG_CHECKING([for the NetCDF Fortran interface C preprocessing flags])
    if test "${sd_netcdf_fortran_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_fortran_cppflags}])
    fi
    AC_MSG_CHECKING([for the NetCDF Fortran interface C flags])
    if test "${sd_netcdf_fortran_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_fortran_cflags}])
    fi
    if test "${sd_netcdf_fortran_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([for the NetCDF Fortran interface C++ flags])
      if test "${sd_netcdf_fortran_cxxflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_netcdf_fortran_cxxflags}])
      fi
    fi
    if test "${sd_netcdf_fortran_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([for the NetCDF Fortran interface Fortran flags])
      if test "${sd_netcdf_fortran_fcflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_netcdf_fortran_fcflags}])
      fi
    fi
    AC_MSG_CHECKING([for the NetCDF Fortran interface linker flags])
    if test "${sd_netcdf_fortran_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_fortran_ldflags}])
    fi
    AC_MSG_CHECKING([for the NetCDF Fortran interface library flags])
    if test "${sd_netcdf_fortran_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_fortran_libs}])
    fi
  fi
]) # _SD_NETCDF_FORTRAN_DUMP_CONFIG

## Copyright (C) 2019 Yann Pouillon

#
# NetCDF I/O library - C++ interface
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_NETCDF_CXX_INIT], [
  # Init
  sd_netcdf_cxx_cppflags=""
  sd_netcdf_cxx_cflags=""
  sd_netcdf_cxx_cxxflags=""
  sd_netcdf_cxx_fcflags=""
  sd_netcdf_cxx_ldflags=""
  sd_netcdf_cxx_libs=""
  sd_netcdf_cxx_prefix=""
  sd_netcdf_cxx_enable=""
  sd_netcdf_cxx_init="unknown"
  sd_netcdf_cxx_mpi_ok="unknown"
  sd_netcdf_cxx_ok="unknown"

  # Set adjustable parameters
  sd_netcdf_cxx_options="$1"
  sd_netcdf_cxx_libs_def="$2"
  sd_netcdf_cxx_cppflags_def="$3"
  sd_netcdf_cxx_cflags_def="$4"
  sd_netcdf_cxx_cxxflags_def="$5"
  sd_netcdf_cxx_fcflags_def="$6"
  sd_netcdf_cxx_ldflags_def="$7"

  # Process options
  for kwd in ${sd_netcdf_cxx_options}; do
    case "${kwd}" in
      auto)
        sd_netcdf_cxx_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_netcdf_cxx_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_netcdf_cxx_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg NetCDF C++ interface option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_netcdf_cxx_status}" && sd_netcdf_cxx_status="optional"
  test -z "${sd_netcdf_cxx_policy}" && sd_netcdf_cxx_policy="fail"
  test -z "${sd_netcdf_cxx_enable_def}" && sd_netcdf_cxx_enable_def="no"
  test -z "${sd_netcdf_cxx_libs_def}" && sd_netcdf_cxx_libs_def="-lnetcdf_c++"
  case "${sd_netcdf_cxx_status}" in
    implicit|required)
      sd_netcdf_cxx_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([netcdf-cxx],
    [AS_HELP_STRING([--with-netcdf-cxx],
      [Install prefix of the NetCDF C++ interface library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_netcdf_cxx_enable="${withval}"
        sd_netcdf_cxx_init="yon"
      else
        sd_netcdf_cxx_enable="yes"
        sd_netcdf_cxx_init="dir"
        sd_netcdf_cxx_prefix="${withval}"
      fi],
    [ sd_netcdf_cxx_enable="${sd_netcdf_cxx_enable_def}"; sd_netcdf_cxx_init="def"])

  # Declare environment variables
  AC_ARG_VAR([NETCDF_CXX_CPPFLAGS], [C preprocessing flags for NetCDF C++ interface.])
  AC_ARG_VAR([NETCDF_CXX_CFLAGS], [C flags for NetCDF C++ interface.])
  AC_ARG_VAR([NETCDF_CXX_CXXFLAGS], [C++ flags for NetCDF C++ interface.])
  AC_ARG_VAR([NETCDF_CXX_FCFLAGS], [Fortran flags for NetCDF C++ interface.])
  AC_ARG_VAR([NETCDF_CXX_LDFLAGS], [Linker flags for NetCDF C++ interface.])
  AC_ARG_VAR([NETCDF_CXX_LIBS], [Library flags for NetCDF C++ interface.])

  # Detect use of environment variables
  if test "${sd_netcdf_cxx_enable}" = "yes" -o "${sd_netcdf_cxx_enable}" = "auto"; then
    tmp_netcdf_cxx_vars="${NETCDF_CXX_CPPFLAGS}${NETCDF_CXX_CFLAGS}${NETCDF_CXX_CXXFLAGS}${NETCDF_CXX_FCFLAGS}${NETCDF_CXX_LDFLAGS}${NETCDF_CXX_LIBS}"
    if test "${sd_netcdf_cxx_init}" = "def" -a \
            ! -z "${tmp_netcdf_cxx_vars}"; then
      sd_netcdf_cxx_enable="yes"
      sd_netcdf_cxx_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_NETCDF_CXX_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_netcdf_cxx_enable}" = "yes" -o "${sd_netcdf_cxx_enable}" = "auto"; then

    # Set NetCDF C++ interface-specific flags
    case "${sd_netcdf_cxx_init}" in

      def|yon)
        sd_netcdf_cxx_cppflags="${sd_netcdf_cxx_cppflags_def}"
        sd_netcdf_cxx_cflags="${sd_netcdf_cxx_cflags_def}"
        sd_netcdf_cxx_cxxflags="${sd_netcdf_cxx_cxxflags_def}"
        sd_netcdf_cxx_fcflags="${sd_netcdf_cxx_fcflags_def}"
        sd_netcdf_cxx_ldflags="${sd_netcdf_cxx_ldflags_def}"
        sd_netcdf_cxx_libs="${sd_netcdf_cxx_libs_def}"
        ;;

      dir)
        sd_netcdf_cxx_cppflags="${sd_netcdf_cxx_cppflags_def} -I${sd_netcdf_cxx_prefix}/include"
        sd_netcdf_cxx_cflags="${sd_netcdf_cxx_cflags_def}"
        sd_netcdf_cxx_cxxflags="${sd_netcdf_cxx_cxxflags_def}"
        sd_netcdf_cxx_fcflags="${sd_netcdf_cxx_fcflags_def} -I${sd_netcdf_cxx_prefix}/include"
        sd_netcdf_cxx_ldflags="${sd_netcdf_cxx_ldflags_def}"
        sd_netcdf_cxx_libs="-L${sd_netcdf_cxx_prefix}/lib ${sd_netcdf_cxx_libs_def} ${sd_netcdf_cxx_libs}"
        ;;

      env)
        sd_netcdf_cxx_cppflags="${sd_netcdf_cxx_cppflags_def}"
        sd_netcdf_cxx_cflags="${sd_netcdf_cxx_cflags_def}"
        sd_netcdf_cxx_cxxflags="${sd_netcdf_cxx_cxxflags_def}"
        sd_netcdf_cxx_fcflags="${sd_netcdf_cxx_fcflags_def}"
        sd_netcdf_cxx_ldflags="${sd_netcdf_cxx_ldflags_def}"
        sd_netcdf_cxx_libs="${sd_netcdf_cxx_libs_def}"
        test ! -z "${NETCDF_CXX_CPPFLAGS}" && sd_netcdf_cxx_cppflags="${NETCDF_CXX_CPPFLAGS}"
        test ! -z "${NETCDF_CXX_CFLAGS}" && sd_netcdf_cxx_cflags="${NETCDF_CXX_CFLAGS}"
        test ! -z "${NETCDF_CXX_CXXFLAGS}" && sd_netcdf_cxx_cxxflags="${NETCDF_CXX_CXXFLAGS}"
        test ! -z "${NETCDF_CXX_FCFLAGS}" && sd_netcdf_cxx_fcflags="${NETCDF_CXX_FCFLAGS}"
        test ! -z "${NETCDF_CXX_LDFLAGS}" && sd_netcdf_cxx_ldflags="${NETCDF_CXX_LDFLAGS}"
        test ! -z "${NETCDF_CXX_LIBS}" && sd_netcdf_cxx_libs="${NETCDF_CXX_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for the NetCDF C++ interface: '${sd_netcdf_cxx_init}'])
        ;;

    esac

  fi

  # Display configuration
  _SD_NETCDF_CXX_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_netcdf_cxx_options)
  AC_SUBST(sd_netcdf_cxx_enable_def)
  AC_SUBST(sd_netcdf_cxx_policy)
  AC_SUBST(sd_netcdf_cxx_status)
  AC_SUBST(sd_netcdf_cxx_enable)
  AC_SUBST(sd_netcdf_cxx_init)
  AC_SUBST(sd_netcdf_cxx_ok)
  AC_SUBST(sd_netcdf_cxx_cppflags)
  AC_SUBST(sd_netcdf_cxx_cflags)
  AC_SUBST(sd_netcdf_cxx_fcflags)
  AC_SUBST(sd_netcdf_cxx_ldflags)
  AC_SUBST(sd_netcdf_cxx_libs)
  AC_SUBST(with_netcdf)

  # Clean-up
  unset tmp_netcdf_cxx_vars
]) # SD_NETCDF_CXX_INIT


AC_DEFUN([SD_NETCDF_CXX_DETECT], [
  # Display configuration
  _SD_NETCDF_CXX_DUMP_CONFIG

  # Check if NetCDF has been enabled
  if test "${sd_netcdf_ok}" != "yes"; then
    if test "${sd_netcdf_cxx_init}" = "def"; then
      sd_netcdf_cxx_enable="no"
    else
      AC_MSG_ERROR([NetCDF is not available])
    fi
  fi

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_netcdf_cxx_enable}" = "auto" -o "${sd_netcdf_cxx_enable}" = "yes"; then
    _SD_NETCDF_CXX_CHECK_USE

    if test "${sd_netcdf_cxx_ok}" = "yes"; then
      CPPFLAGS="${CPPFLAGS} ${sd_netcdf_cxx_cppflags}"
      CFLAGS="${CFLAGS} ${sd_netcdf_cxx_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_netcdf_cxx_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_netcdf_cxx_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_netcdf_cxx_ldflags}"
      LIBS="${sd_netcdf_cxx_libs} ${LIBS}"

      AC_DEFINE([HAVE_NETCDF_CXX], 1,
        [Define to 1 if you have the NetCDF C++ interface library.])

      if test "${sd_mpi_ok}" = "yes" -a "${sd_netcdf_cxx_mpi_ok}" = "yes"; then
        AC_DEFINE([HAVE_NETCDF_CXX_MPI], 1,
          [Define to 1 if you have a parallel NetCDF C++ interface library.])
      fi
    else
      if test "${sd_netcdf_cxx_status}" = "optional" -a \
              "${sd_netcdf_cxx_init}" = "def"; then
        sd_netcdf_cxx_enable="no"
        sd_netcdf_cxx_cppflags=""
        sd_netcdf_cxx_cflags=""
        sd_netcdf_cxx_cxxflags=""
        sd_netcdf_cxx_fcflags=""
        sd_netcdf_cxx_ldflags=""
        sd_netcdf_cxx_libs=""
      else
        AC_MSG_FAILURE([invalid NetCDF C++ interface configuration])
      fi
    fi
  else
    sd_netcdf_cxx_enable="no"
    sd_netcdf_cxx_cppflags=""
    sd_netcdf_cxx_cflags=""
    sd_netcdf_cxx_cxxflags=""
    sd_netcdf_cxx_fcflags=""
    sd_netcdf_cxx_ldflags=""
    sd_netcdf_cxx_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_NETCDF_CXX_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_netcdf_cxx_cppflags}"
  CFLAGS="${CFLAGS} ${sd_netcdf_cxx_cflags}"
  CXXFLAGS="${CXXFLAGS} ${sd_netcdf_cxx_cxxflags}"
  FCFLAGS="${FCFLAGS} ${sd_netcdf_cxx_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_netcdf_cxx_ldflags}"
  LIBS="${sd_netcdf_cxx_libs} ${LIBS}"

  # Check NetCDF C++ API
  AC_MSG_CHECKING([whether the NetCDF C++ interface works])
  if test "${sd_netcdf_ok}" = "yes"; then
    AC_LANG_PUSH([C++])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#       include <netcdf4>
        using namespace netcdf;
      ]],
      [[
        NcFile dataFile("conftest.nc", NcFile::replace);
        Create(dataFile);
      ]])], [sd_netcdf_cxx_ok="yes"], [sd_netcdf_cxx_ok="no"])
    AC_LANG_POP([C++])
  else
    sd_netcdf_cxx_ok="no"
  fi
  AC_MSG_RESULT([${sd_netcdf_cxx_ok}])

  # Check if we can do parallel I/O
  if test "${sd_netcdf_cxx_ok}" = "yes" -a "${sd_mpi_ok}" = "yes"; then
    if test "${sd_hdf5_mpi_ok}" = "yes"; then
      AC_MSG_CHECKING([whether the NetCDF C++ interface has parallel I/O in Fortran])
      AC_LANG_PUSH([C++])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
#         include <mpi.h>
#         include <netcdf4>
          using namespace::netcdf;
        ]],
        [[
          MPI::Comm & mpiComm = MPI::COMM_WORLD;
          mpiComm.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
          MPI::Info mpiInfo = MPI::INFO_NULL;
          netcdf::NcFile ncFile(mpiComm, mpiInfo, "conftest.nc", netcdf::NcFile::Replace);
        ]])], [sd_netcdf_cxx_mpi_ok="yes"], [sd_netcdf_cxx_mpi_ok="no"])
      AC_LANG_POP([C++])
      AC_MSG_RESULT([${sd_netcdf_cxx_mpi_ok}])
    else
      sd_netcdf_cxx_mpi_ok="no"
    fi
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_NETCDF_CXX_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_NETCDF_CXX_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_netcdf_cxx_invalid="no"
  if test "${sd_netcdf_cxx_enable_def}" != "auto" -a \
          "${sd_netcdf_cxx_enable_def}" != "no" -a \
          "${sd_netcdf_cxx_enable_def}" != "yes"; then
    case "${sd_netcdf_cxx_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_netcdf_cxx_enable_def = '${sd_netcdf_cxx_enable_def}'])
        ;;
      skip)
        tmp_netcdf_cxx_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_netcdf_cxx_enable_def = '${sd_netcdf_cxx_enable_def}'])
        tmp_netcdf_cxx_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_netcdf_cxx_invalid}" = "yes"; then
    if test "${sd_netcdf_cxx_status}" = "required"; then
      sd_netcdf_cxx_enable_def="yes"
    else
      sd_netcdf_cxx_enable_def="no"
    fi
    tmp_netcdf_cxx_invalid="no"
    AC_MSG_NOTICE([setting sd_netcdf_cxx_enable_def to '${sd_netcdf_cxx_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_netcdf_cxx_invalid="no"
  if test "${sd_netcdf_cxx_status}" = "implicit" -o \
          "${sd_netcdf_cxx_status}" = "required"; then
    if test "${sd_netcdf_cxx_enable}" = "no"; then
      case "${sd_netcdf_cxx_policy}" in
        fail)
          AC_MSG_ERROR([The NetCDF C++ interface package is required and cannot be disabled
                  See https://www.unidata.ucar.edu/software/netcdf/ for
                  details on how to install it.])
          ;;
        skip)
          tmp_netcdf_cxx_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The NetCDF C++ interface package is required and cannot be disabled])
          tmp_netcdf_cxx_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_netcdf_cxx_enable}" = "auto"; then
      AC_MSG_NOTICE([setting NetCDF C++ interface trigger to yes])
      sd_netcdf_cxx_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_netcdf_cxx_invalid}" = "yes"; then
    case "${sd_netcdf_cxx_status}" in
      implicit|required)
        sd_netcdf_cxx_enable="yes"
        ;;
      optional)
        sd_netcdf_cxx_enable="no"
        ;;
    esac
    tmp_netcdf_cxx_invalid="no"
    AC_MSG_NOTICE([setting sd_netcdf_cxx_enable to '${sd_netcdf_cxx_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_netcdf_cxx_vars="${NETCDF_CXX_CPPFLAGS}${NETCDF_CXX_CFLAGS}${NETCDF_CXX_FCFLAGS}${NETCDF_CXX_LDFLAGS}${NETCDF_CXX_LIBS}"
  tmp_netcdf_cxx_invalid="no"
  if test ! -z "${tmp_netcdf_cxx_vars}" -a ! -z "${sd_netcdf_cxx_prefix}"; then
    case "${sd_netcdf_cxx_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for NetCDF C++ interface
                  Please use NETCDF_CXX_FCFLAGS + NETCDF_CXX_LIBS or --with-netcdf,
                  not both.])
        ;;
      skip)
        tmp_netcdf_cxx_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for NetCDF C++ interface])
        tmp_netcdf_cxx_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_netcdf_cxx_vars}"; then
    sd_netcdf_cxx_enable="yes"
    sd_netcdf_cxx_init="env"
    if test "${tmp_netcdf_cxx_invalid}" = "yes"; then
      tmp_netcdf_cxx_invalid="no"
      AC_MSG_NOTICE([overriding --with-netcdf with NETCDF_CXX_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_netcdf_cxx_status}" = "implicit"; then
    if test "${sd_netcdf_cxx_fcflags}" != ""; then
      sd_netcdf_cxx_fcflags=""
      AC_MSG_NOTICE([resetting NetCDF C++ interface Fortran flags (implicit package)])
    fi
    if test "${sd_netcdf_cxx_ldflags}" != ""; then
      sd_netcdf_cxx_ldflags=""
      AC_MSG_NOTICE([resetting NetCDF C++ interface linker flags (implicit package)])
    fi
    if test "${sd_netcdf_cxx_libs}" != ""; then
      sd_netcdf_cxx_libs=""
      AC_MSG_NOTICE([resetting NetCDF C++ interface library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_netcdf_cxx_enable}" = "implicit"; then
    sd_netcdf_cxx_fcflags=""
    sd_netcdf_cxx_ldflags=""
    sd_netcdf_cxx_libs=""
  fi

  # Clean-up
  unset tmp_netcdf_cxx_invalid
  unset tmp_netcdf_cxx_vars
]) # _SD_NETCDF_CXX_CHECK_CONFIG


AC_DEFUN([_SD_NETCDF_CXX_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable the NetCDF C++ interface])
  AC_MSG_RESULT([${sd_netcdf_cxx_enable}])
  if test "${sd_netcdf_cxx_enable}" != "no"; then
    AC_MSG_CHECKING([how the NetCDF C++ interface parameters have been set])
    AC_MSG_RESULT([${sd_netcdf_cxx_init}])
    AC_MSG_CHECKING([for the NetCDF C++ interface C preprocessing flags])
    if test "${sd_netcdf_cxx_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_cxx_cppflags}])
    fi
    AC_MSG_CHECKING([for the NetCDF C++ interface C flags])
    if test "${sd_netcdf_cxx_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_cxx_cflags}])
    fi
    if test "${sd_netcdf_cxx_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([for the NetCDF C++ interface C++ flags])
      if test "${sd_netcdf_cxx_cxxflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_netcdf_cxx_cxxflags}])
      fi
    fi
    if test "${sd_netcdf_cxx_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([for the NetCDF C++ interface Fortran flags])
      if test "${sd_netcdf_cxx_fcflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_netcdf_cxx_fcflags}])
      fi
    fi
    AC_MSG_CHECKING([for the NetCDF C++ interface linker flags])
    if test "${sd_netcdf_cxx_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_cxx_ldflags}])
    fi
    AC_MSG_CHECKING([for the NetCDF C++ interface library flags])
    if test "${sd_netcdf_cxx_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_cxx_libs}])
    fi
  fi
]) # _SD_NETCDF_CXX_DUMP_CONFIG

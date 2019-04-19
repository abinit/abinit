## Copyright (C) 2019 Yann Pouillon

#
# Exchange-Correlation functionals library (NetCDF)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_NETCDF_INIT], [
  # Init
  sd_netcdf_cppflags=""
  sd_netcdf_cflags=""
  sd_netcdf_cxxflags=""
  sd_netcdf_fcflags=""
  sd_netcdf_ldflags=""
  sd_netcdf_libs=""
  sd_netcdf_enable=""
  sd_netcdf_init="unknown"
  sd_netcdf_c_ok="unknown"
  sd_netcdf_fortran_ok="unknown"
  sd_netcdf_ok="unknown"

  # Set adjustable parameters
  sd_netcdf_options="$1"
  sd_netcdf_libs_def="$2"
  sd_netcdf_cppflags_def="$3"
  sd_netcdf_cflags_def="$4"
  sd_netcdf_cxxflags_def="$5"
  sd_netcdf_fcflags_def="$6"
  sd_netcdf_ldflags_def="$7"

  # Process options
  sd_netcdf_enable_cxx="yes"
  sd_netcdf_enable_def="auto"
  sd_netcdf_enable_fc="yes"
  sd_netcdf_policy="fail"
  sd_netcdf_status="optional"
  for kwd in ${sd_netcdf_options}; do
    case "${kwd}" in
      auto|no|yes)
        sd_netcdf_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_netcdf_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_netcdf_policy="${kwd}"
        ;;
      no-cxx)
        sd_netcdf_enable_cxx="no"
        ;;
      no-fortran)
        sd_netcdf_enable_fc="no"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg NetCDF option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_netcdf_enable_def}" && sd_netcdf_enable_def="auto"
  test -z "${sd_netcdf_status}" && sd_netcdf_status="optional"
  test -z "${sd_netcdf_policy}" && sd_netcdf_policy="fail"
  if test "${sd_netcdf_enable_fc}" = "yes"; then
    test -z "${sd_netcdf_libs_def}" && sd_netcdf_libs_def="-lnetcdff -lnetcdf"
  else
    test -z "${sd_netcdf_libs_def}" && sd_netcdf_libs_def="-lnetcdf"
  fi

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([netcdf],
    [AS_HELP_STRING([--with-netcdf],
      [Install prefix of the NetCDF library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_netcdf_enable="${withval}"
        sd_netcdf_init="yon"
      else
        sd_netcdf_enable="yes"
        sd_netcdf_init="dir"
      fi],
    [ sd_netcdf_enable="${sd_netcdf_enable_def}"; sd_netcdf_init="def"])

  # Declare environment variables
  AC_ARG_VAR([NETCDF_CPPFLAGS], [C preprocessing flags for NetCDF.])
  AC_ARG_VAR([NETCDF_CFLAGS], [C flags for NetCDF.])
  AC_ARG_VAR([NETCDF_CXXFLAGS], [C++ flags for NetCDF.])
  AC_ARG_VAR([NETCDF_FCFLAGS], [Fortran flags for NetCDF.])
  AC_ARG_VAR([NETCDF_LDFLAGS], [Linker flags for NetCDF.])
  AC_ARG_VAR([NETCDF_LIBS], [Library flags for NetCDF.])

  # Detect use of environment variables
  if test "${sd_netcdf_enable}" = "yes" -o "${sd_netcdf_enable}" = "auto"; then
    tmp_netcdf_vars="${NETCDF_CPPFLAGS}${NETCDF_CFLAGS}${NETCDF_LDFLAGS}${NETCDF_LIBS}"
    if test "${sd_netcdf_enable_cxx}" = "yes"; then
      tmp_netcdf_vars="${tmp_netcdf_vars}${NETCDF_CXXFLAGS}"
    fi
    if test "${sd_netcdf_enable_fc}" = "yes"; then
      tmp_netcdf_vars="${tmp_netcdf_vars}${NETCDF_FCFLAGS}"
    fi
    if test "${sd_netcdf_init}" = "def" -a ! -z "${tmp_netcdf_vars}"; then
      sd_netcdf_enable="yes"
      sd_netcdf_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_NETCDF_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_netcdf_enable}" = "yes" -o "${sd_netcdf_enable}" = "auto"; then

    # Set NetCDF-specific flags
    case "${sd_netcdf_init}" in

      def|yon)
        sd_netcdf_cppflags="${sd_netcdf_cppflags_def}"
        sd_netcdf_cflags="${sd_netcdf_cflags_def}"
        test "${sd_netcdf_enable_cxx}" = "yes" && \
          sd_netcdf_cxxflags="${sd_netcdf_cxxflags_def}"
        test "${sd_netcdf_enable_fc}" = "yes" && \
          sd_netcdf_fcflags="${sd_netcdf_fcflags_def}"
        sd_netcdf_ldflags="${sd_netcdf_ldflags_def}"
        sd_netcdf_libs="${sd_netcdf_libs_def}"
        ;;

      dir)
        sd_netcdf_cppflags="-I${with_netcdf}/include"
        sd_netcdf_cflags="${sd_netcdf_cflags_def}"
        test "${sd_netcdf_enable_cxx}" = "yes" && \
          sd_netcdf_cxxflags="${sd_netcdf_cxxflags_def}"
        test "${sd_netcdf_enable_fc}" = "yes" && \
          sd_netcdf_fcflags="${sd_netcdf_fcflags_def} -I${with_netcdf}/include"
        sd_netcdf_ldflags="${sd_netcdf_ldflags_def}"
        sd_netcdf_libs="-L${with_netcdf}/lib ${sd_netcdf_libs_def}"
        ;;

      env)
        sd_netcdf_cppflags="${sd_netcdf_cppflags_def}"
        sd_netcdf_cflags="${sd_netcdf_cflags_def}"
        test "${sd_netcdf_enable_cxx}" = "yes" && \
          sd_netcdf_cxxflags="${sd_netcdf_cxxflags_def}"
        test "${sd_netcdf_enable_fc}" = "yes" && \
          sd_netcdf_fcflags="${sd_netcdf_fcflags_def}"
        sd_netcdf_ldflags="${sd_netcdf_ldflags_def}"
        sd_netcdf_libs="${sd_netcdf_libs_def}"
        test ! -z "${NETCDF_CPPFLAGS}" && sd_netcdf_cppflags="${NETCDF_CPPFLAGS}"
        test ! -z "${NETCDF_CFLAGS}" && sd_netcdf_cflags="${NETCDF_CFLAGS}"
        if test "${sd_netcdf_enable_cxx}" = "yes"; then
          test ! -z "${NETCDF_CXXFLAGS}" && sd_netcdf_cxxflags="${NETCDF_CXXFLAGS}"
        fi
        if test "${sd_netcdf_enable_fc}" = "yes"; then
          test ! -z "${NETCDF_FCFLAGS}" && sd_netcdf_fcflags="${NETCDF_FCFLAGS}"
        fi
        test ! -z "${NETCDF_LDFLAGS}" && sd_netcdf_ldflags="${NETCDF_LDFLAGS}"
        test ! -z "${NETCDF_LIBS}" && sd_netcdf_libs="${NETCDF_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for NetCDF: '${sd_netcdf_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_netcdf_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_netcdf_init="esl"
    sd_netcdf_cppflags=""
    sd_netcdf_cflags=""
    sd_netcdf_cxxflags=""
    sd_netcdf_fcflags=""
    sd_netcdf_ldflags=""
    sd_netcdf_libs=""
  fi

  # Display configuration
  _SD_NETCDF_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_netcdf_options)
  AC_SUBST(sd_netcdf_enable_def)
  AC_SUBST(sd_netcdf_enable_cxx)
  AC_SUBST(sd_netcdf_enable_fc)
  AC_SUBST(sd_netcdf_policy)
  AC_SUBST(sd_netcdf_status)
  AC_SUBST(sd_netcdf_enable)
  AC_SUBST(sd_netcdf_init)
  AC_SUBST(sd_netcdf_ok)
  AC_SUBST(sd_netcdf_cppflags)
  AC_SUBST(sd_netcdf_cflags)
  AC_SUBST(sd_netcdf_fcflags)
  AC_SUBST(sd_netcdf_ldflags)
  AC_SUBST(sd_netcdf_libs)
  AC_SUBST(with_netcdf)

  # Clean-up
  unset tmp_netcdf_vars
]) # SD_NETCDF_INIT


AC_DEFUN([SD_NETCDF_DETECT], [
  # Display configuration
  _SD_NETCDF_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_netcdf_enable}" = "auto" -o "${sd_netcdf_enable}" = "yes"; then
    _SD_NETCDF_CHECK_USE

    if test "${sd_netcdf_ok}" = "yes"; then
      if test "${sd_netcdf_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_netcdf_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_netcdf_fcflags}"
        LIBS="${sd_netcdf_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_netcdf_ldflags}"

      AC_DEFINE([HAVE_NETCDF], 1,
        [Define to 1 if you have the NetCDF library.])
    else
      if test "${sd_netcdf_status}" = "optional" -a \
              "${sd_netcdf_init}" = "def"; then
        sd_netcdf_enable="no"
        sd_netcdf_cppflags=""
        sd_netcdf_cflags=""
        sd_netcdf_cxxflags=""
        sd_netcdf_fcflags=""
        sd_netcdf_ldflags=""
        sd_netcdf_libs=""
      else
        AC_MSG_FAILURE([invalid NetCDF configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_NETCDF_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_netcdf_init}" = "esl"; then
    AC_MSG_NOTICE([will look for NetCDF in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_netcdf_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_netcdf_cppflags}"
    CFLAGS="${CFLAGS} ${sd_netcdf_cflags}"
    FCFLAGS="${FCFLAGS} ${sd_netcdf_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_netcdf_ldflags}"
    LIBS="${sd_netcdf_libs} ${LIBS}"
  fi

  # Check NetCDF C API
  AC_MSG_CHECKING([whether the NetCDF library works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#     include <netcdf.h>
    ]],
    [[
      int ncid;
      return nc_open('conftest.nc', NC_WRITE, ncid);
    ]])], [sd_netcdf_c_ok="yes"], [sd_netcdf_c_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_netcdf_c_ok}])

  # Check NetCDF Fortran API
  if test "${sd_netcdf_c_ok}" = "yes" -a "${sd_netcdf_enable_fc}" = "yes"; then
    AC_MSG_CHECKING([whether the NetCDF Fortran interface works])
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
    AC_MSG_RESULT([${sd_netcdf_fortran_ok}])
  fi
  unset tmp_incs

  # Combine the available results
  sd_netcdf_ok="no"
  if test "${sd_netcdf_enable_fc}" = "yes"; then
    if test "${sd_netcdf_c_ok}" = "yes" -a \
            "${sd_netcdf_fortran_ok}" = "yes"; then
      sd_netcdf_ok="yes"
    fi
  else
    if test "${sd_netcdf_c_ok}" = "yes"; then
      sd_netcdf_ok="yes"
    fi
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_NETCDF_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_NETCDF_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_netcdf_invalid="no"
  if test "${sd_netcdf_enable_def}" != "auto" -a \
          "${sd_netcdf_enable_def}" != "no" -a \
          "${sd_netcdf_enable_def}" != "yes"; then
    case "${sd_netcdf_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_netcdf_enable_def = '${sd_netcdf_enable_def}'])
        ;;
      skip)
        tmp_netcdf_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_netcdf_enable_def = '${sd_netcdf_enable_def}'])
        tmp_netcdf_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_netcdf_invalid}" = "yes"; then
    if test "${sd_netcdf_status}" = "required"; then
      sd_netcdf_enable_def="yes"
    else
      sd_netcdf_enable_def="no"
    fi
    tmp_netcdf_invalid="no"
    AC_MSG_NOTICE([setting sd_netcdf_enable_def to '${sd_netcdf_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_netcdf_invalid="no"
  if test "${sd_netcdf_status}" = "implicit" -o \
          "${sd_netcdf_status}" = "required"; then
    if test "${sd_netcdf_enable}" = "no"; then
      case "${sd_netcdf_policy}" in
        fail)
          AC_MSG_ERROR([The NetCDF package is required and cannot be disabled
                  See https://www.unidata.ucar.edu/software/netcdf/ for
                  details on how to install it.])
          ;;
        skip)
          tmp_netcdf_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The NetCDF package is required and cannot be disabled])
          tmp_netcdf_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_netcdf_enable}" = "auto"; then
      AC_MSG_NOTICE([setting NetCDF trigger to yes])
      sd_netcdf_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_netcdf_invalid}" = "yes"; then
    case "${sd_netcdf_status}" in
      implicit|required)
        sd_netcdf_enable="yes"
        ;;
      optional)
        sd_netcdf_enable="no"
        ;;
    esac
    tmp_netcdf_invalid="no"
    AC_MSG_NOTICE([setting sd_netcdf_enable to '${sd_netcdf_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_netcdf_vars="${NETCDF_FCFLAGS}${NETCDF_LDFLAGS}${NETCDF_LIBS}"
  tmp_netcdf_invalid="no"
  if test ! -z "${tmp_netcdf_vars}" -a ! -z "${with_netcdf}"; then
    case "${sd_netcdf_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for NetCDF
                  Please use NETCDF_FCFLAGS + NETCDF_LIBS or --with-netcdf,
                  not both.])
        ;;
      skip)
        tmp_netcdf_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for NetCDF])
        tmp_netcdf_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_netcdf_vars}"; then
    sd_netcdf_enable="yes"
    sd_netcdf_init="env"
    if test "${tmp_netcdf_invalid}" = "yes"; then
      tmp_netcdf_invalid="no"
      AC_MSG_NOTICE([overriding --with-netcdf with NETCDF_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_netcdf_status}" = "implicit"; then
    if test "${sd_netcdf_fcflags}" != ""; then
      sd_netcdf_fcflags=""
      AC_MSG_NOTICE([resetting NetCDF Fortran flags (implicit package)])
    fi
    if test "${sd_netcdf_ldflags}" != ""; then
      sd_netcdf_ldflags=""
      AC_MSG_NOTICE([resetting NetCDF linker flags (implicit package)])
    fi
    if test "${sd_netcdf_libs}" != ""; then
      sd_netcdf_libs=""
      AC_MSG_NOTICE([resetting NetCDF library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_netcdf_enable}" = "implicit"; then
    sd_netcdf_fcflags=""
    sd_netcdf_ldflags=""
    sd_netcdf_libs=""
  fi

  # Clean-up
  unset tmp_netcdf_invalid
  unset tmp_netcdf_vars
]) # _SD_NETCDF_CHECK_CONFIG


AC_DEFUN([_SD_NETCDF_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable NetCDF])
  AC_MSG_RESULT([${sd_netcdf_enable}])
  if test "${sd_netcdf_enable}" != "no"; then
    AC_MSG_CHECKING([how NetCDF parameters have been set])
    AC_MSG_RESULT([${sd_netcdf_init}])
    AC_MSG_CHECKING([for NetCDF C preprocessing flags])
    if test "${sd_netcdf_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_cppflags}])
    fi
    AC_MSG_CHECKING([for NetCDF C flags])
    if test "${sd_netcdf_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_cflags}])
    fi
    if test "${sd_netcdf_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([for NetCDF C++ flags])
      if test "${sd_netcdf_cxxflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_netcdf_cxxflags}])
      fi
    fi
    if test "${sd_netcdf_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([for NetCDF Fortran flags])
      if test "${sd_netcdf_fcflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_netcdf_fcflags}])
      fi
    fi
    AC_MSG_CHECKING([for NetCDF linker flags])
    if test "${sd_netcdf_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_ldflags}])
    fi
    AC_MSG_CHECKING([for NetCDF library flags])
    if test "${sd_netcdf_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_netcdf_libs}])
    fi
  fi
]) # _SD_NETCDF_DUMP_CONFIG

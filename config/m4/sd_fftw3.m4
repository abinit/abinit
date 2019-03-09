## Copyright (C) 2019 Yann Pouillon

#
# Fastest Fourier Transform in the West library (FFTW3)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_FFTW3_DETECT], [
  # Display configuration
  _SD_FFTW3_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  _SD_FFTW3_CHECK_USE

  # Update build flags
  if test "${sd_fftw3_ok}" = "yes"; then
    FCFLAGS="${FCFLAGS} ${sd_fftw3_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_fftw3_ldflags}"
    LIBS="${sd_fftw3_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([invalid FFTW3 configuration])
  fi
])


AC_DEFUN([SD_FFTW3_INIT], [
  # Init
  sd_fftw3_cflags=""
  sd_fftw3_cppflags=""
  sd_fftw3_enable="unknown"
  sd_fftw3_enable_def="$1"
  sd_fftw3_fcflags=""
  sd_fftw3_fcflags_def="$2"
  sd_fftw3_init="unknown"
  sd_fftw3_ldflags=""
  sd_fftw3_ldflags_def="$3"
  sd_fftw3_libs=""
  sd_fftw3_libs_def="$4"
  sd_fftw3_mpi_ok="unknown"
  sd_fftw3_ok="unknown"
  sd_fftw3_policy="$6"
  sd_fftw3_status="$5"

  # Make sure default values are correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_FFTW3_CHECK_DEFAULTS
  fi

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([fftw3],
    [AS_HELP_STRING([--with-fftw3],
      [Install prefix of the Flexible Data Format I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_fftw3_enable="${withval}"
        sd_fftw3_init="yon"
      else
        sd_fftw3_enable="yes"
        sd_fftw3_init="dir"
      fi],
    [ sd_fftw3_enable="${sd_fftw3_enable_def}"; sd_fftw3_init="def"])

  # Declare environment variables
  AC_ARG_VAR([FFTW3_FCFLAGS],
    [Fortran flags for FFTW3 (conflicts with the --with-fftw3 option).])
  AC_ARG_VAR([FFTW3_LDFLAGS],
    [Linker flags for FFTW3 to prepend at link-time (conflicts with the --with-fftw3 option).])
  AC_ARG_VAR([FFTW3_LIBS],
    [Library flags for FFTW3 to append at link-time (conflicts with the --with-fftw3 option).])

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_FFTW3_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  case "${sd_fftw3_init}" in
    def|yon)
      sd_fftw3_fcflags="${sd_fftw3_fcflags_def}"
      sd_fftw3_ldflags="${sd_fftw3_ldflags_def}"
      sd_fftw3_libs="${sd_fftw3_libs_def}"
      ;;
    dir)
      sd_fftw3_fcflags="-I${with_fftw3}/include"
      sd_fftw3_ldflags="${sd_fftw3_ldflags_def}"
      sd_fftw3_libs="-L${with_fftw3}/lib ${sd_fftw3_libs_def}"
      ;;
    env)
      sd_fftw3_fcflags="${FFTW3_FCFLAGS}"
      sd_fftw3_ldflags="${FFTW3_LDFLAGS}"
      sd_fftw3_libs="${FFTW3_LIBS}"
      ;;
    *)
      AC_MSG_ERROR([invalid init type for FFTW3: '${sd_fftw3_init}'])
      ;;
  esac

  # Display configuration
  _SD_FFTW3_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_fftw3_enable)
  AC_SUBST(sd_fftw3_enable_def)
  AC_SUBST(sd_fftw3_fcflags)
  AC_SUBST(sd_fftw3_init)
  AC_SUBST(sd_fftw3_ldflags)
  AC_SUBST(sd_fftw3_libs)
  AC_SUBST(sd_fftw3_ok)
  AC_SUBST(with_fftw3)
]) # SD_FFTW3_INIT


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_FFTW3_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_fftw3_invalid="no"
  if test "${sd_fftw3_enable_def}" != "auto" -a \
          "${sd_fftw3_enable_def}" != "no" -a \
          "${sd_fftw3_enable_def}" != "yes"; then
    case "${sd_fftw3_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_fftw3_enable_def = '${sd_fftw3_enable_def}'])
        ;;
      skip)
        tmp_fftw3_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_fftw3_enable_def = '${sd_fftw3_enable_def}'])
        tmp_fftw3_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_fftw3_invalid}" = "yes"; then
    if test "${sd_fftw3_status}" = "required"; then
      sd_fftw3_enable_def="yes"
    else
      sd_fftw3_enable_def="no"
    fi
    tmp_fftw3_invalid="no"
    AC_MSG_NOTICE([setting sd_fftw3_enable_def to '${sd_fftw3_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_fftw3_invalid="no"
  if test "${sd_fftw3_status}" = "implicit" -o \
          "${sd_fftw3_status}" = "required"; then
    if test "${sd_fftw3_enable}" = "no"; then
      case "${sd_fftw3_policy}" in
        fail)
          AC_MSG_ERROR([The FFTW3 package is required and cannot be disabled
                  See https://launchpad.net/fftw3 for details on how to
                  install it.])
          ;;
        skip)
          tmp_fftw3_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The FFTW3 package is required and cannot be disabled])
          tmp_fftw3_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_fftw3_enable}" = "auto"; then
      AC_MSG_NOTICE([setting FFTW3 trigger to yes])
      sd_fftw3_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_fftw3_invalid}" = "yes"; then
    case "${sd_fftw3_status}" in
      implicit|required)
        sd_fftw3_enable="yes"
        ;;
      optional)
        sd_fftw3_enable="no"
        ;;
    esac
    tmp_fftw3_invalid="no"
    AC_MSG_NOTICE([setting sd_fftw3_enable to '${sd_fftw3_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_fftw3_vars="${FFTW3_FCFLAGS}${FFTW3_LDFLAGS}${FFTW3_LIBS}"
  tmp_fftw3_invalid="no"
  if test ! -z "${tmp_fftw3_vars}" -a ! -z "${with_fftw3}"; then
    case "${sd_fftw3_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for FFTW3
                  Please use FFTW3_FCFLAGS + FFTW3_LIBS or --with-fftw3,
                  not both.])
        ;;
      skip)
        tmp_fftw3_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for FFTW3])
        tmp_fftw3_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_fftw3_vars}"; then
    sd_fftw3_enable="yes"
    sd_fftw3_init="env"
    if test "${tmp_fftw3_invalid}" = "yes"; then
      tmp_fftw3_invalid="no"
      AC_MSG_NOTICE([overriding --with-fftw3 with FFTW3_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_fftw3_status}" = "implicit"; then
    if test "${sd_fftw3_fcflags}" != ""; then
      sd_fftw3_fcflags=""
      AC_MSG_NOTICE([resetting FFTW3 Fortran flags (implicit package)])
    fi
    if test "${sd_fftw3_ldflags}" != ""; then
      sd_fftw3_ldflags=""
      AC_MSG_NOTICE([resetting FFTW3 linker flags (implicit package)])
    fi
    if test "${sd_fftw3_libs}" != ""; then
      sd_fftw3_libs=""
      AC_MSG_NOTICE([resetting FFTW3 library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_fftw3_enable}" = "implicit"; then
    sd_fftw3_fcflags=""
    sd_fftw3_ldflags=""
    sd_fftw3_libs=""
  fi

  # Clean-up
  unset tmp_fftw3_invalid
  unset tmp_fftw3_vars
]) # _SD_FFTW3_CHECK_CONFIG


AC_DEFUN([_SD_FFTW3_CHECK_DEFAULTS], [
  # Policy and status must be defined before initialisation
  # Note: this is a developer error, hence we abort unconditionally
  if test "${sd_fftw3_policy}" != "fail" -a \
          "${sd_fftw3_policy}" != "skip" -a \
          "${sd_fftw3_policy}" != "warn"; then
    AC_MSG_ERROR([invalid policy sd_fftw3_policy='${sd_fftw3_policy}'
                  Valid policies for broken configurations are:
                      'fail', 'skip', or 'warn'])
  fi
  if test "${sd_fftw3_status}" != "implicit" -a \
          "${sd_fftw3_status}" != "optional" -a \
          "${sd_fftw3_status}" != "required"; then
    AC_MSG_ERROR([invalid policy sd_fftw3_status='${sd_fftw3_status}'
                  Valid dependency statuses are:
                      'implicit', 'optional', or 'required'])
  fi

  # Set default trigger if undefined
  if test -z "${sd_fftw3_enable_def}"; then
    case "${sd_fftw3_status}" in
      optional)
        sd_fftw3_enable_def="auto"
        ;;
      required)
        sd_fftw3_enable_def="yes"
        ;;
    esac
  fi
]) # _SD_FFTW3_CHECK_DEFAULTS


AC_DEFUN([_SD_FFTW3_CHECK_USE], [
  # Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_CFLAGS="${CFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LDFLAGS="${LDFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${sd_fftw3_cppflags}"
  CFLAGS="${CFLAGS} ${sd_fftw3_cflags}"
  FCFLAGS="${FCFLAGS} ${sd_fftw3_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_fftw3_ldflags}"
  LIBS="${sd_fftw3_libs} ${LIBS}"

  # Check FFTW3 C API
  AC_MSG_CHECKING([whether the FFTW3 library works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <fftw3.h>
    ]],
    [[
      fftw_plan *plan;
      fftw_complex *a1, *a2;
      fftw_execute_dft(plan, a1, a2);
    ]])], [sd_fftw3_ok="yes"], [sd_fftw3_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_fftw3_ok}])

  # Check FFTW3 MPI C API
  if test "${sd_fftw3_ok}" = "yes" -a "${sd_mpi_enable}" = "yes"; then
    AC_MSG_CHECKING([whether the FFTW3 MPI library works])
    AC_LANG_PUSH([C])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#include <fftw3-mpi.h>
      ]],
      [[
        fftw_mpi_init();
      ]])], [sd_fftw3_mpi_ok="yes"], [sd_fftw3_mpi_ok="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${sd_fftw3_ok}])
  else
    sd_fftw3_mpi_ok="no"
  fi

  # Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  CFLAGS="${tmp_saved_CFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LDFLAGS="${tmp_saved_LDFLAGS}"
  LIBS="${tmp_saved_LIBS}"
  unset tmp_saved_FCFLAGS
  unset tmp_saved_LDFLAGS
  unset tmp_saved_LIBS
]) # _SD_FFTW3_CHECK_USE


AC_DEFUN([_SD_FFTW3_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether FFTW3 is enabled])
  AC_MSG_RESULT([${sd_fftw3_enable}])
  AC_MSG_CHECKING([how FFTW3 parameters have been set])
  AC_MSG_RESULT([${sd_fftw3_init}])
  AC_MSG_CHECKING([for FFTW3 C preprocessing flags])
  AC_MSG_RESULT([${sd_fftw3_cppflags}])
  AC_MSG_CHECKING([for FFTW3 C flags])
  AC_MSG_RESULT([${sd_fftw3_cflags}])
  AC_MSG_CHECKING([for FFTW3 Fortran flags])
  AC_MSG_RESULT([${sd_fftw3_fcflags}])
  AC_MSG_CHECKING([for FFTW3 linker flags])
  AC_MSG_RESULT([${sd_fftw3_ldflags}])
  AC_MSG_CHECKING([for FFTW3 library flags])
  AC_MSG_RESULT([${sd_fftw3_libs}])
]) # _SD_FFTW3_DUMP_CONFIG

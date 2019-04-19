## Copyright (C) 2019 Yann Pouillon

#
# Fastest Fourier Transform in the West library (PFFT)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_PFFT_DETECT], [
  # Display configuration
  _SD_PFFT_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  _SD_PFFT_CHECK_USE

  # Update build flags
  if test "${sd_pfft_ok}" = "yes"; then
    FCFLAGS="${FCFLAGS} ${sd_pfft_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_pfft_ldflags}"
    LIBS="${sd_pfft_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([invalid PFFT configuration])
  fi
])


AC_DEFUN([SD_PFFT_INIT], [
  # Init
  sd_pfft_cflags=""
  sd_pfft_cppflags=""
  sd_pfft_enable="unknown"
  sd_pfft_enable_def="$1"
  sd_pfft_fcflags=""
  sd_pfft_fcflags_def="$2"
  sd_pfft_init="unknown"
  sd_pfft_ldflags=""
  sd_pfft_ldflags_def="$3"
  sd_pfft_libs=""
  sd_pfft_libs_def="$4"
  sd_pfft_ok="unknown"
  sd_pfft_policy="$6"
  sd_pfft_status="$5"

  # Make sure default values are correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_PFFT_CHECK_DEFAULTS
  fi

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([pfft],
    [AS_HELP_STRING([--with-pfft],
      [Install prefix of the Flexible Data Format I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_pfft_enable="${withval}"
        sd_pfft_init="yon"
      else
        sd_pfft_enable="yes"
        sd_pfft_init="dir"
      fi],
    [ sd_pfft_enable="${sd_pfft_enable_def}"; sd_pfft_init="def"])

  # Declare environment variables
  AC_ARG_VAR([PFFT_FCFLAGS],
    [Fortran flags for PFFT (conflicts with the --with-pfft option).])
  AC_ARG_VAR([PFFT_LDFLAGS],
    [Linker flags for PFFT to prepend at link-time (conflicts with the --with-pfft option).])
  AC_ARG_VAR([PFFT_LIBS],
    [Library flags for PFFT to append at link-time (conflicts with the --with-pfft option).])

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_PFFT_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  case "${sd_pfft_init}" in
    def|yon)
      sd_pfft_fcflags="${sd_pfft_fcflags_def}"
      sd_pfft_ldflags="${sd_pfft_ldflags_def}"
      sd_pfft_libs="${sd_pfft_libs_def}"
      ;;
    dir)
      sd_pfft_fcflags="-I${with_pfft}/include"
      sd_pfft_ldflags="${sd_pfft_ldflags_def}"
      sd_pfft_libs="-L${with_pfft}/lib ${sd_pfft_libs_def}"
      ;;
    env)
      sd_pfft_fcflags="${PFFT_FCFLAGS}"
      sd_pfft_ldflags="${PFFT_LDFLAGS}"
      sd_pfft_libs="${PFFT_LIBS}"
      ;;
    *)
      AC_MSG_ERROR([invalid init type for PFFT: '${sd_pfft_init}'])
      ;;
  esac

  # Display configuration
  _SD_PFFT_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_pfft_enable)
  AC_SUBST(sd_pfft_enable_def)
  AC_SUBST(sd_pfft_fcflags)
  AC_SUBST(sd_pfft_init)
  AC_SUBST(sd_pfft_ldflags)
  AC_SUBST(sd_pfft_libs)
  AC_SUBST(sd_pfft_ok)
  AC_SUBST(with_pfft)
]) # SD_PFFT_INIT


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_PFFT_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_pfft_invalid="no"
  if test "${sd_pfft_enable_def}" != "auto" -a \
          "${sd_pfft_enable_def}" != "no" -a \
          "${sd_pfft_enable_def}" != "yes"; then
    case "${sd_pfft_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_pfft_enable_def = '${sd_pfft_enable_def}'])
        ;;
      skip)
        tmp_pfft_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_pfft_enable_def = '${sd_pfft_enable_def}'])
        tmp_pfft_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_pfft_invalid}" = "yes"; then
    if test "${sd_pfft_status}" = "required"; then
      sd_pfft_enable_def="yes"
    else
      sd_pfft_enable_def="no"
    fi
    tmp_pfft_invalid="no"
    AC_MSG_NOTICE([setting sd_pfft_enable_def to '${sd_pfft_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_pfft_invalid="no"
  if test "${sd_pfft_status}" = "implicit" -o \
          "${sd_pfft_status}" = "required"; then
    if test "${sd_pfft_enable}" = "no"; then
      case "${sd_pfft_policy}" in
        fail)
          AC_MSG_ERROR([The PFFT package is required and cannot be disabled
                  See https://github.com/mpip/pfft for details on how to
                  install it.])
          ;;
        skip)
          tmp_pfft_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The PFFT package is required and cannot be disabled])
          tmp_pfft_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_pfft_enable}" = "auto"; then
      AC_MSG_NOTICE([setting PFFT trigger to yes])
      sd_pfft_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_pfft_invalid}" = "yes"; then
    case "${sd_pfft_status}" in
      implicit|required)
        sd_pfft_enable="yes"
        ;;
      optional)
        sd_pfft_enable="no"
        ;;
    esac
    tmp_pfft_invalid="no"
    AC_MSG_NOTICE([setting sd_pfft_enable to '${sd_pfft_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_pfft_vars="${PFFT_FCFLAGS}${PFFT_LDFLAGS}${PFFT_LIBS}"
  tmp_pfft_invalid="no"
  if test ! -z "${tmp_pfft_vars}" -a ! -z "${with_pfft}"; then
    case "${sd_pfft_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for PFFT
                  Please use PFFT_FCFLAGS + PFFT_LIBS or --with-pfft,
                  not both.])
        ;;
      skip)
        tmp_pfft_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for PFFT])
        tmp_pfft_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_pfft_vars}"; then
    sd_pfft_enable="yes"
    sd_pfft_init="env"
    if test "${tmp_pfft_invalid}" = "yes"; then
      tmp_pfft_invalid="no"
      AC_MSG_NOTICE([overriding --with-pfft with PFFT_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_pfft_status}" = "implicit"; then
    if test "${sd_pfft_fcflags}" != ""; then
      sd_pfft_fcflags=""
      AC_MSG_NOTICE([resetting PFFT Fortran flags (implicit package)])
    fi
    if test "${sd_pfft_ldflags}" != ""; then
      sd_pfft_ldflags=""
      AC_MSG_NOTICE([resetting PFFT linker flags (implicit package)])
    fi
    if test "${sd_pfft_libs}" != ""; then
      sd_pfft_libs=""
      AC_MSG_NOTICE([resetting PFFT library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_pfft_enable}" = "implicit"; then
    sd_pfft_fcflags=""
    sd_pfft_ldflags=""
    sd_pfft_libs=""
  fi

  # Clean-up
  unset tmp_pfft_invalid
  unset tmp_pfft_vars
]) # _SD_PFFT_CHECK_CONFIG


AC_DEFUN([_SD_PFFT_CHECK_DEFAULTS], [
  # Policy and status must be defined before initialisation
  # Note: this is a developer error, hence we abort unconditionally
  if test "${sd_pfft_policy}" != "fail" -a \
          "${sd_pfft_policy}" != "skip" -a \
          "${sd_pfft_policy}" != "warn"; then
    AC_MSG_ERROR([invalid policy sd_pfft_policy='${sd_pfft_policy}'
                  Valid policies for broken configurations are:
                      'fail', 'skip', or 'warn'])
  fi
  if test "${sd_pfft_status}" != "implicit" -a \
          "${sd_pfft_status}" != "optional" -a \
          "${sd_pfft_status}" != "required"; then
    AC_MSG_ERROR([invalid policy sd_pfft_status='${sd_pfft_status}'
                  Valid dependency statuses are:
                      'implicit', 'optional', or 'required'])
  fi

  # Set default trigger if undefined
  if test -z "${sd_pfft_enable_def}"; then
    case "${sd_pfft_status}" in
      optional)
        sd_pfft_enable_def="auto"
        ;;
      required)
        sd_pfft_enable_def="yes"
        ;;
    esac
  fi
]) # _SD_PFFT_CHECK_DEFAULTS


AC_DEFUN([_SD_PFFT_CHECK_USE], [
  # Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_CFLAGS="${CFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LDFLAGS="${LDFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${sd_pfft_cppflags}"
  CFLAGS="${CFLAGS} ${sd_pfft_cflags}"
  FCFLAGS="${FCFLAGS} ${sd_pfft_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_pfft_ldflags}"
  LIBS="${sd_pfft_libs} ${LIBS}"

  # Check PFFT C API
  AC_MSG_CHECKING([whether the PFFT library works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <pfft.h>
    ]],
    [[
      pfft_init();
    ]])], [sd_pfft_ok="yes"], [sd_pfft_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_pfft_ok}])

  # Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  CFLAGS="${tmp_saved_CFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LDFLAGS="${tmp_saved_LDFLAGS}"
  LIBS="${tmp_saved_LIBS}"
  unset tmp_saved_FCFLAGS
  unset tmp_saved_LDFLAGS
  unset tmp_saved_LIBS
]) # _SD_PFFT_CHECK_USE


AC_DEFUN([_SD_PFFT_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether PFFT is enabled])
  AC_MSG_RESULT([${sd_pfft_enable}])
  AC_MSG_CHECKING([how PFFT parameters have been set])
  AC_MSG_RESULT([${sd_pfft_init}])
  AC_MSG_CHECKING([for PFFT C preprocessing flags])
  AC_MSG_RESULT([${sd_pfft_cppflags}])
  AC_MSG_CHECKING([for PFFT C flags])
  AC_MSG_RESULT([${sd_pfft_cflags}])
  AC_MSG_CHECKING([for PFFT Fortran flags])
  AC_MSG_RESULT([${sd_pfft_fcflags}])
  AC_MSG_CHECKING([for PFFT linker flags])
  AC_MSG_RESULT([${sd_pfft_ldflags}])
  AC_MSG_CHECKING([for PFFT library flags])
  AC_MSG_RESULT([${sd_pfft_libs}])
]) # _SD_PFFT_DUMP_CONFIG

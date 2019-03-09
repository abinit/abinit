## Copyright (C) 2019 Yann Pouillon

#
# Flexible Data Format I/O library (LibXC)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_LIBXC_DETECT], [
  # Display configuration
  _SD_LIBXC_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_libxc_enable}" = "yes"; then
    _SD_LIBXC_CHECK_USE

    if test "${sd_libxc_ok}" = "yes"; then
      if test "${sd_libxc_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_libxc_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_libxc_fcflags}"
        LIBS="${sd_libxc_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_libxc_ldflags}"
    else
      AC_MSG_FAILURE([invalid LibXC configuration])
    fi
  fi
])


AC_DEFUN([SD_LIBXC_INIT], [
  # Init
  sd_libxc_cflags=""
  sd_libxc_cppflags=""
  sd_libxc_enable="unknown"
  sd_libxc_enable_def="$1"
  sd_libxc_fcflags=""
  sd_libxc_fcflags_def="$2"
  sd_libxc_fortran_ok="unknown"
  sd_libxc_init="unknown"
  sd_libxc_ldflags=""
  sd_libxc_ldflags_def="$3"
  sd_libxc_libs=""
  sd_libxc_libs_def="$4"
  sd_libxc_ok="unknown"
  sd_libxc_policy="$6"
  sd_libxc_status="$5"

  # Make sure default values are correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_LIBXC_CHECK_DEFAULTS
  fi

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([libxc],
    [AS_HELP_STRING([--with-libxc],
      [Install prefix of the Flexible Data Format I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_libxc_enable="${withval}"
        sd_libxc_init="yon"
      else
        sd_libxc_enable="yes"
        sd_libxc_init="dir"
      fi],
    [ sd_libxc_enable="${sd_libxc_enable_def}"; sd_libxc_init="def"])

  # Declare environment variables
  AC_ARG_VAR([LIBXC_FCFLAGS],
    [Fortran flags for LibXC (conflicts with the --with-libxc option).])
  AC_ARG_VAR([LIBXC_LDFLAGS],
    [Linker flags for LibXC to prepend at link-time (conflicts with the --with-libxc option).])
  AC_ARG_VAR([LIBXC_LIBS],
    [Library flags for LibXC to append at link-time (conflicts with the --with-libxc option).])

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_LIBXC_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  case "${sd_libxc_init}" in
    def|yon)
      sd_libxc_fcflags="${sd_libxc_fcflags_def}"
      sd_libxc_ldflags="${sd_libxc_ldflags_def}"
      sd_libxc_libs="${sd_libxc_libs_def}"
      ;;
    dir)
      sd_libxc_fcflags="-I${with_libxc}/include"
      sd_libxc_ldflags="${sd_libxc_ldflags_def}"
      sd_libxc_libs="-L${with_libxc}/lib ${sd_libxc_libs_def}"
      ;;
    env)
      sd_libxc_fcflags="${LIBXC_FCFLAGS}"
      sd_libxc_ldflags="${LIBXC_LDFLAGS}"
      sd_libxc_libs="${LIBXC_LIBS}"
      ;;
    *)
      AC_MSG_ERROR([invalid init type for LibXC: '${sd_libxc_init}'])
      ;;
  esac

  # Override the default configuration if the ESL Bundle is available
  if test "${sd_libxc_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_libxc_init="esl"
    sd_libxc_fcflags=""
    sd_libxc_ldflags=""
    sd_libxc_libs=""
  fi

  # Display configuration
  _SD_LIBXC_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_libxc_enable)
  AC_SUBST(sd_libxc_enable_def)
  AC_SUBST(sd_libxc_fcflags)
  AC_SUBST(sd_libxc_init)
  AC_SUBST(sd_libxc_ldflags)
  AC_SUBST(sd_libxc_libs)
  AC_SUBST(sd_libxc_ok)
  AC_SUBST(with_libxc)
]) # SD_LIBXC_INIT


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_LIBXC_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_libxc_invalid="no"
  if test "${sd_libxc_enable_def}" != "auto" -a \
          "${sd_libxc_enable_def}" != "no" -a \
          "${sd_libxc_enable_def}" != "yes"; then
    case "${sd_libxc_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_libxc_enable_def = '${sd_libxc_enable_def}'])
        ;;
      skip)
        tmp_libxc_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_libxc_enable_def = '${sd_libxc_enable_def}'])
        tmp_libxc_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_libxc_invalid}" = "yes"; then
    if test "${sd_libxc_status}" = "required"; then
      sd_libxc_enable_def="yes"
    else
      sd_libxc_enable_def="no"
    fi
    tmp_libxc_invalid="no"
    AC_MSG_NOTICE([setting sd_libxc_enable_def to '${sd_libxc_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_libxc_invalid="no"
  if test "${sd_libxc_status}" = "implicit" -o \
          "${sd_libxc_status}" = "required"; then
    if test "${sd_libxc_enable}" = "no"; then
      case "${sd_libxc_policy}" in
        fail)
          AC_MSG_ERROR([The LibXC package is required and cannot be disabled
                  See https://launchpad.net/libxc for details on how to
                  install it.])
          ;;
        skip)
          tmp_libxc_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The LibXC package is required and cannot be disabled])
          tmp_libxc_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_libxc_enable}" = "auto"; then
      AC_MSG_NOTICE([setting LibXC trigger to yes])
      sd_libxc_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_libxc_invalid}" = "yes"; then
    case "${sd_libxc_status}" in
      implicit|required)
        sd_libxc_enable="yes"
        ;;
      optional)
        sd_libxc_enable="no"
        ;;
    esac
    tmp_libxc_invalid="no"
    AC_MSG_NOTICE([setting sd_libxc_enable to '${sd_libxc_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_libxc_vars="${LIBXC_FCFLAGS}${LIBXC_LDFLAGS}${LIBXC_LIBS}"
  tmp_libxc_invalid="no"
  if test ! -z "${tmp_libxc_vars}" -a ! -z "${with_libxc}"; then
    case "${sd_libxc_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for LibXC
                  Please use LIBXC_FCFLAGS + LIBXC_LIBS or --with-libxc,
                  not both.])
        ;;
      skip)
        tmp_libxc_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for LibXC])
        tmp_libxc_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_libxc_vars}"; then
    sd_libxc_enable="yes"
    sd_libxc_init="env"
    if test "${tmp_libxc_invalid}" = "yes"; then
      tmp_libxc_invalid="no"
      AC_MSG_NOTICE([overriding --with-libxc with LIBXC_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_libxc_status}" = "implicit"; then
    if test "${sd_libxc_fcflags}" != ""; then
      sd_libxc_fcflags=""
      AC_MSG_NOTICE([resetting LibXC Fortran flags (implicit package)])
    fi
    if test "${sd_libxc_ldflags}" != ""; then
      sd_libxc_ldflags=""
      AC_MSG_NOTICE([resetting LibXC linker flags (implicit package)])
    fi
    if test "${sd_libxc_libs}" != ""; then
      sd_libxc_libs=""
      AC_MSG_NOTICE([resetting LibXC library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_libxc_enable}" = "implicit"; then
    sd_libxc_fcflags=""
    sd_libxc_ldflags=""
    sd_libxc_libs=""
  fi

  # Clean-up
  unset tmp_libxc_invalid
  unset tmp_libxc_vars
]) # _SD_LIBXC_CHECK_CONFIG


AC_DEFUN([_SD_LIBXC_CHECK_DEFAULTS], [
  # Policy and status must be defined before initialisation
  # Note: this is a developer error, hence we abort unconditionally
  if test "${sd_libxc_policy}" != "fail" -a \
          "${sd_libxc_policy}" != "skip" -a \
          "${sd_libxc_policy}" != "warn"; then
    AC_MSG_ERROR([invalid policy sd_libxc_policy='${sd_libxc_policy}'
                  Valid policies for broken configurations are:
                      'fail', 'skip', or 'warn'])
  fi
  if test "${sd_libxc_status}" != "implicit" -a \
          "${sd_libxc_status}" != "optional" -a \
          "${sd_libxc_status}" != "required"; then
    AC_MSG_ERROR([invalid policy sd_libxc_status='${sd_libxc_status}'
                  Valid dependency statuses are:
                      'implicit', 'optional', or 'required'])
  fi

  # Set default trigger if undefined
  if test -z "${sd_libxc_enable_def}"; then
    case "${sd_libxc_status}" in
      optional)
        sd_libxc_enable_def="auto"
        ;;
      required)
        sd_libxc_enable_def="yes"
        ;;
    esac
  fi
]) # _SD_LIBXC_CHECK_DEFAULTS


AC_DEFUN([_SD_LIBXC_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_libxc_init}" = "esl"; then
    AC_MSG_NOTICE([will look for LibXC in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_libxc_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_libxc_cppflags}"
    CFLAGS="${CFLAGS} ${sd_libxc_cflags}"
    FCFLAGS="${FCFLAGS} ${sd_libxc_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_libxc_ldflags}"
    LIBS="${sd_libxc_libs} ${LIBS}"
  fi

  # Check LibXC C API
  AC_MSG_CHECKING([whether the LibXC library works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <xc.h>
    ]],
    [[
      int ma, mi, mu;
      xc_version(&ma, &mi, &mu);
    ]])], [sd_libxc_ok="yes"], [sd_libxc_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_libxc_ok}])

  # Check LibXC Fortran API
  AC_MSG_CHECKING([whether the LibXC Fortran interface works])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use libxc_m
      implicit none
      integer :: ma, mi, mu
      call xc_f90_version(ma, mi, mu)
    ]])], [sd_libxc_fortran_ok="yes"], [sd_libxc_fortran_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_libxc_fortran_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_LIBXC_CHECK_USE


AC_DEFUN([_SD_LIBXC_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether LibXC is enabled])
  AC_MSG_RESULT([${sd_libxc_enable}])
  AC_MSG_CHECKING([how LibXC parameters have been set])
  AC_MSG_RESULT([${sd_libxc_init}])
  AC_MSG_CHECKING([for LibXC C preprocessing flags])
  AC_MSG_RESULT([${sd_libxc_cppflags}])
  AC_MSG_CHECKING([for LibXC C flags])
  AC_MSG_RESULT([${sd_libxc_cflags}])
  AC_MSG_CHECKING([for LibXC Fortran flags])
  AC_MSG_RESULT([${sd_libxc_fcflags}])
  AC_MSG_CHECKING([for LibXC linker flags])
  AC_MSG_RESULT([${sd_libxc_ldflags}])
  AC_MSG_CHECKING([for LibXC library flags])
  AC_MSG_RESULT([${sd_libxc_libs}])
]) # _SD_LIBXC_DUMP_CONFIG

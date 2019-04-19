## Copyright (C) 2019 Yann Pouillon

#
# Performance Application Programming Interface (PAPI)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_PAPI_INIT], [
  # Init
  sd_papi_cppflags=""
  sd_papi_cflags=""
  sd_papi_ldflags=""
  sd_papi_libs=""
  sd_papi_enable=""
  sd_papi_init="unknown"
  sd_papi_ok="unknown"

  # Set adjustable parameters
  sd_papi_options="$1"
  sd_papi_libs_def="$2"
  sd_papi_cppflags_def="$3"
  sd_papi_cflags_def="$4"
  sd_papi_cxxflags_def="$5"
  sd_papi_fcflags_def="$6"
  sd_papi_ldflags_def="$7"
  sd_papi_enable_def=""
  sd_papi_policy=""
  sd_papi_status=""

  # Process options
  for kwd in ${sd_papi_options}; do
    case "${kwd}" in
      auto|no|yes)
        sd_papi_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_papi_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_papi_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg PAPI option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_papi_enable_def}" && sd_papi_enable_def="auto"
  test -z "${sd_papi_policy}" && sd_papi_policy="fail"
  test -z "${sd_papi_status}" && sd_papi_status="optional"
  test -z "${sd_papi_libs_def}" && sd_papi_libs_def="-lpapi"

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([papi],
    [AS_HELP_STRING([--with-papi],
      [Install prefix of the Flexible Data Format I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_papi_enable="${withval}"
        sd_papi_init="yon"
      else
        sd_papi_enable="yes"
        sd_papi_init="dir"
      fi],
    [ sd_papi_enable="${sd_papi_enable_def}"; sd_papi_init="def"])

  # Declare environment variables
  AC_ARG_VAR([PAPI_CPPFLAGS], [C preprocessing flags for PAPI.])
  AC_ARG_VAR([PAPI_CFLAGS], [C flags for PAPI.])
  AC_ARG_VAR([PAPI_LDFLAGS], [Linker flags for PAPI.])
  AC_ARG_VAR([PAPI_LIBS], [Library flags for PAPI.])

  # Detect use of environment variables
  if test "${sd_papi_enable}" = "yes" -o "${sd_papi_enable}" = "auto"; then
    tmp_papi_vars="${PAPI_CPPFLAGS}${PAPI_CFLAGS}${PAPI_LDFLAGS}${PAPI_LIBS}"
    if test "${sd_papi_init}" = "def" -a ! -z "${tmp_papi_vars}"; then
      sd_papi_enable="yes"
      sd_papi_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_PAPI_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_papi_enable}" = "yes" -o "${sd_papi_enable}" = "auto"; then

    # Set PAPI-specific flags
    case "${sd_papi_init}" in

      def|yon)
        sd_papi_cppflags="${sd_papi_cppflags_def}"
        sd_papi_cflags="${sd_papi_cflags_def}"
        sd_papi_ldflags="${sd_papi_ldflags_def}"
        sd_papi_libs="${sd_papi_libs_def}"
        ;;

      dir)
        sd_papi_cppflags="-I${with_papi}/include"
        sd_papi_cflags="${sd_papi_cflags_def}"
        sd_papi_ldflags="${sd_papi_ldflags_def}"
        sd_papi_libs="-L${with_papi}/lib ${sd_papi_libs_def}"
        ;;

      env)
        sd_papi_cppflags="${sd_papi_cppflags_def}"
        sd_papi_cflags="${sd_papi_cflags_def}"
        sd_papi_ldflags="${sd_papi_ldflags_def}"
        sd_papi_libs="${sd_papi_libs_def}"
        test ! -z "${PAPI_CPPFLAGS}" && sd_papi_cppflags="${PAPI_CPPFLAGS}"
        test ! -z "${PAPI_CFLAGS}" && sd_papi_cflags="${PAPI_CFLAGS}"
        test ! -z "${PAPI_LDFLAGS}" && sd_papi_ldflags="${PAPI_LDFLAGS}"
        test ! -z "${PAPI_LIBS}" && sd_papi_libs="${PAPI_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for PAPI: '${sd_papi_init}'])
        ;;

    esac

  fi

  # Display configuration
  _SD_PAPI_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_papi_options)
  AC_SUBST(sd_papi_enable_def)
  AC_SUBST(sd_papi_policy)
  AC_SUBST(sd_papi_status)
  AC_SUBST(sd_papi_enable)
  AC_SUBST(sd_papi_init)
  AC_SUBST(sd_papi_ok)
  AC_SUBST(sd_papi_cppflags)
  AC_SUBST(sd_papi_cflags)
  AC_SUBST(sd_papi_ldflags)
  AC_SUBST(sd_papi_libs)
  AC_SUBST(with_papi)

  # Clean-up
  unset tmp_papi_vars
]) # SD_PAPI_INIT


AC_DEFUN([SD_PAPI_DETECT], [
  # Display configuration
  _SD_PAPI_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_papi_enable}" = "auto" -o "${sd_papi_enable}" = "yes"; then
    _SD_PAPI_CHECK_USE

    if test "${sd_papi_ok}" = "yes"; then
      if test "${sd_papi_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_papi_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_papi_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_papi_ldflags}"
    else
      if test "${sd_papi_status}" = "optional" -a \
              "${sd_papi_init}" = "def"; then
        sd_papi_enable="no"
        sd_papi_cppflags=""
        sd_papi_cflags=""
        sd_papi_ldflags=""
        sd_papi_libs=""
      else
        AC_MSG_FAILURE([invalid PAPI configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_PAPI_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_papi_init}" = "esl"; then
    AC_MSG_NOTICE([will look for PAPI in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_papi_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_papi_cppflags}"
    CFLAGS="${CFLAGS} ${sd_papi_cflags}"
    LDFLAGS="${LDFLAGS} ${sd_papi_ldflags}"
    LIBS="${sd_papi_libs} ${LIBS}"
  fi

  # Add rt support if available on the machine
  if test "${PAPI_LIBS}" = ""; then
    AC_LANG_PUSH([C])
    AC_CHECK_HEADERS([time.h])
    AC_CHECK_LIB([rt], [clock_gettime], [sd_papi_rt_libs="-lrt"], [sd_papi_rt_libs=""])
    AC_CHECK_FUNCS([clock_gettime])
    AC_LANG_POP([C])
    sd_papi_libs="${sd_papi_rt_libs} ${sd_papi_libs}"
  fi

  # Check C headers
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([papi.h f90papi.h],
    [sd_papi_has_hdrs="yes"], [sd_papi_has_hdrs="no"])
  AC_LANG_POP([C])

  # Check Fortran support
  if test "${sd_papi_has_hdrs}" = "yes"; then
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([whether the specified PAPI library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
#       if defined HAVE_F90PAPI_H
#       include "f90papi.h"
#       endif
        call PAPIf_library_init
      ]])], [sd_papi_ok="yes"], [sd_papi_ok="no"])
    AC_MSG_RESULT([${sd_papi_ok}])
    AC_LANG_POP([Fortran])
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_PAPI_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_PAPI_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_papi_invalid="no"
  if test "${sd_papi_enable_def}" != "auto" -a \
          "${sd_papi_enable_def}" != "no" -a \
          "${sd_papi_enable_def}" != "yes"; then
    case "${sd_papi_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_papi_enable_def = '${sd_papi_enable_def}'])
        ;;
      skip)
        tmp_papi_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_papi_enable_def = '${sd_papi_enable_def}'])
        tmp_papi_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_papi_invalid}" = "yes"; then
    if test "${sd_papi_status}" = "required"; then
      sd_papi_enable_def="yes"
    else
      sd_papi_enable_def="no"
    fi
    tmp_papi_invalid="no"
    AC_MSG_NOTICE([setting sd_papi_enable_def to '${sd_papi_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_papi_invalid="no"
  if test "${sd_papi_status}" = "implicit" -o \
          "${sd_papi_status}" = "required"; then
    if test "${sd_papi_enable}" = "no"; then
      case "${sd_papi_policy}" in
        fail)
          AC_MSG_ERROR([The PAPI package is required and cannot be disabled
                  See https://icl.utk.edu/papi/index.html for details on how
                  to install it.])
          ;;
        skip)
          tmp_papi_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The PAPI package is required and cannot be disabled])
          tmp_papi_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_papi_enable}" = "auto"; then
      AC_MSG_NOTICE([setting PAPI trigger to yes])
      sd_papi_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_papi_invalid}" = "yes"; then
    case "${sd_papi_status}" in
      implicit|required)
        sd_papi_enable="yes"
        ;;
      optional)
        sd_papi_enable="no"
        ;;
    esac
    tmp_papi_invalid="no"
    AC_MSG_NOTICE([setting sd_papi_enable to '${sd_papi_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_papi_vars="${PAPI_CFLAGS}${PAPI_LDFLAGS}${PAPI_LIBS}"
  tmp_papi_invalid="no"
  if test ! -z "${tmp_papi_vars}" -a ! -z "${with_papi}"; then
    case "${sd_papi_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for PAPI
                  Please use PAPI_CFLAGS + PAPI_LIBS or --with-papi,
                  not both.])
        ;;
      skip)
        tmp_papi_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for PAPI])
        tmp_papi_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_papi_vars}"; then
    sd_papi_enable="yes"
    sd_papi_init="env"
    if test "${tmp_papi_invalid}" = "yes"; then
      tmp_papi_invalid="no"
      AC_MSG_NOTICE([overriding --with-papi with PAPI_{CPPFLAGS,CFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_papi_status}" = "implicit"; then
    if test "${sd_papi_ldflags}" != ""; then
      sd_papi_ldflags=""
      AC_MSG_NOTICE([resetting PAPI linker flags (implicit package)])
    fi
    if test "${sd_papi_libs}" != ""; then
      sd_papi_libs=""
      AC_MSG_NOTICE([resetting PAPI library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_papi_enable}" = "implicit"; then
    sd_papi_ldflags=""
    sd_papi_libs=""
  fi

  # Clean-up
  unset tmp_papi_invalid
  unset tmp_papi_vars
]) # _SD_PAPI_CHECK_CONFIG


AC_DEFUN([_SD_PAPI_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable PAPI])
  AC_MSG_RESULT([${sd_papi_enable}])
  if test "${sd_papi_enable}" != "no"; then
    AC_MSG_CHECKING([how PAPI parameters have been set])
    AC_MSG_RESULT([${sd_papi_init}])
    AC_MSG_CHECKING([for PAPI C preprocessing flags])
    if test "${sd_papi_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_papi_cppflags}])
    fi
    AC_MSG_CHECKING([for PAPI C flags])
    if test "${sd_papi_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_papi_cflags}])
    fi
    AC_MSG_CHECKING([for PAPI linker flags])
    if test "${sd_papi_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_papi_ldflags}])
    fi
    AC_MSG_CHECKING([for PAPI library flags])
    if test "${sd_papi_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_papi_libs}])
    fi
  fi
]) # _SD_PAPI_DUMP_CONFIG

## Copyright (C) 2019 Yann Pouillon

#
# Massively Parallel FFT library (PFFT)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_PFFT_INIT], [
  # Init
  sd_pfft_cppflags=""
  sd_pfft_cflags=""
  sd_pfft_ldflags=""
  sd_pfft_libs=""
  sd_pfft_enable=""
  sd_pfft_init="unknown"
  sd_pfft_ok="unknown"

  # Set adjustable parameters
  sd_pfft_options="$1"
  sd_pfft_libs_def="$2"
  sd_pfft_cppflags_def="$3"
  sd_pfft_cflags_def="$4"
  sd_pfft_cxxflags_def="$5"
  sd_pfft_fcflags_def="$6"
  sd_pfft_ldflags_def="$7"
  sd_pfft_enable_def=""
  sd_pfft_policy=""
  sd_pfft_status=""

  # Process options
  for kwd in ${sd_pfft_options}; do
    case "${kwd}" in
      auto|no|yes)
        sd_pfft_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_pfft_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_pfft_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg PFFT option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_pfft_enable_def}" && sd_pfft_enable_def="auto"
  test -z "${sd_pfft_policy}" && sd_pfft_policy="fail"
  test -z "${sd_pfft_status}" && sd_pfft_status="optional"
  test -z "${sd_pfft_libs_def}" && sd_pfft_libs_def="-lpfft"

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([pfft],
    [AS_HELP_STRING([--with-pfft],
      [Install prefix of the PFFT library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_pfft_enable="${withval}"
        sd_pfft_init="yon"
      else
        sd_pfft_enable="yes"
        sd_pfft_init="dir"
      fi],
    [ sd_pfft_enable="${sd_pfft_enable_def}"; sd_pfft_init="def"])

  # Declare environment variables
  AC_ARG_VAR([PFFT_CPPFLAGS], [C preprocessing flags for PFFT.])
  AC_ARG_VAR([PFFT_CFLAGS], [C flags for PFFT.])
  AC_ARG_VAR([PFFT_LDFLAGS], [Linker flags for PFFT.])
  AC_ARG_VAR([PFFT_LIBS], [Library flags for PFFT.])

  # Detect use of environment variables
  if test "${sd_pfft_enable}" = "yes" -o "${sd_pfft_enable}" = "auto"; then
    tmp_pfft_vars="${PFFT_CPPFLAGS}${PFFT_CFLAGS}${PFFT_LDFLAGS}${PFFT_LIBS}"
    if test "${sd_pfft_init}" = "def" -a ! -z "${tmp_pfft_vars}"; then
      sd_pfft_enable="yes"
      sd_pfft_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_PFFT_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_pfft_enable}" = "yes" -o "${sd_pfft_enable}" = "auto"; then

    # Set PFFT-specific flags
    case "${sd_pfft_init}" in

      def|yon)
        sd_pfft_cppflags="${sd_pfft_cppflags_def}"
        sd_pfft_cflags="${sd_pfft_cflags_def}"
        sd_pfft_ldflags="${sd_pfft_ldflags_def}"
        sd_pfft_libs="${sd_pfft_libs_def}"
        ;;

      dir)
        sd_pfft_cppflags="-I${with_pfft}/include"
        sd_pfft_cflags="${sd_pfft_cflags_def}"
        sd_pfft_ldflags="${sd_pfft_ldflags_def}"
        sd_pfft_libs="-L${with_pfft}/lib ${sd_pfft_libs_def}"
        ;;

      env)
        sd_pfft_cppflags="${sd_pfft_cppflags_def}"
        sd_pfft_cflags="${sd_pfft_cflags_def}"
        sd_pfft_ldflags="${sd_pfft_ldflags_def}"
        sd_pfft_libs="${sd_pfft_libs_def}"
        test ! -z "${PFFT_CPPFLAGS}" && sd_pfft_cppflags="${PFFT_CPPFLAGS}"
        test ! -z "${PFFT_CFLAGS}" && sd_pfft_cflags="${PFFT_CFLAGS}"
        test ! -z "${PFFT_LDFLAGS}" && sd_pfft_ldflags="${PFFT_LDFLAGS}"
        test ! -z "${PFFT_LIBS}" && sd_pfft_libs="${PFFT_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for PFFT: '${sd_pfft_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_pfft_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_pfft_init="esl"
    sd_pfft_cppflags=""
    sd_pfft_cflags=""
    sd_pfft_ldflags=""
    sd_pfft_libs=""
  fi

  # Display configuration
  _SD_PFFT_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_pfft_options)
  AC_SUBST(sd_pfft_enable_def)
  AC_SUBST(sd_pfft_policy)
  AC_SUBST(sd_pfft_status)
  AC_SUBST(sd_pfft_enable)
  AC_SUBST(sd_pfft_init)
  AC_SUBST(sd_pfft_ok)
  AC_SUBST(sd_pfft_cppflags)
  AC_SUBST(sd_pfft_cflags)
  AC_SUBST(sd_pfft_ldflags)
  AC_SUBST(sd_pfft_libs)
  AC_SUBST(with_pfft)

  # Clean-up
  unset tmp_pfft_vars
]) # SD_PFFT_INIT


AC_DEFUN([SD_PFFT_DETECT], [
  # Display configuration
  _SD_PFFT_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_pfft_enable}" = "auto" -o "${sd_pfft_enable}" = "yes"; then
    _SD_PFFT_CHECK_USE

    if test "${sd_pfft_ok}" = "yes"; then
      if test "${sd_pfft_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_pfft_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_pfft_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_pfft_ldflags}"
    else
      if test "${sd_pfft_status}" = "optional" -a \
              "${sd_pfft_init}" = "def"; then
        sd_pfft_enable="no"
        sd_pfft_cppflags=""
        sd_pfft_cflags=""
        sd_pfft_ldflags=""
        sd_pfft_libs=""
      else
        AC_MSG_FAILURE([invalid PFFT configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_PFFT_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_pfft_init}" = "esl"; then
    AC_MSG_NOTICE([will look for PFFT in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_pfft_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_pfft_cppflags}"
    CFLAGS="${CFLAGS} ${sd_pfft_cflags}"
    LDFLAGS="${LDFLAGS} ${sd_pfft_ldflags}"
    LIBS="${sd_pfft_libs} ${LIBS}"
  fi

  # Check PFFT C API
  AC_MSG_CHECKING([whether the PFFT library works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#     include <pfft.h>
    ]],
    [[
      pfft_init();
    ]])], [sd_pfft_ok="yes"], [sd_pfft_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_pfft_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_PFFT_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
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
                  See https://github.com/mpip/pfft for details on how
                  to install it.])
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
  tmp_pfft_vars="${PFFT_CFLAGS}${PFFT_LDFLAGS}${PFFT_LIBS}"
  tmp_pfft_invalid="no"
  if test ! -z "${tmp_pfft_vars}" -a ! -z "${with_pfft}"; then
    case "${sd_pfft_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for PFFT
                  Please use PFFT_CFLAGS + PFFT_LIBS or --with-pfft,
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
      AC_MSG_NOTICE([overriding --with-pfft with PFFT_{CPPFLAGS,CFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_pfft_status}" = "implicit"; then
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
    sd_pfft_cppflags=""
    sd_pfft_cflags=""
    sd_pfft_ldflags=""
    sd_pfft_libs=""
  fi

  # Clean-up
  unset tmp_pfft_invalid
  unset tmp_pfft_vars
]) # _SD_PFFT_CHECK_CONFIG


AC_DEFUN([_SD_PFFT_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable PFFT])
  AC_MSG_RESULT([${sd_pfft_enable}])
  if test "${sd_pfft_enable}" != "no"; then
    AC_MSG_CHECKING([how PFFT parameters have been set])
    AC_MSG_RESULT([${sd_pfft_init}])
    AC_MSG_CHECKING([for PFFT C preprocessing flags])
    if test "${sd_pfft_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_pfft_cppflags}])
    fi
    AC_MSG_CHECKING([for PFFT C flags])
    if test "${sd_pfft_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_pfft_cflags}])
    fi
    AC_MSG_CHECKING([for PFFT linker flags])
    if test "${sd_pfft_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_pfft_ldflags}])
    fi
    AC_MSG_CHECKING([for PFFT library flags])
    if test "${sd_pfft_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_pfft_libs}])
    fi
  fi
]) # _SD_PFFT_DUMP_CONFIG

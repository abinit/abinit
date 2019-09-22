## Copyright (C) 2019 Yann Pouillon

#
# ABINIT Common low-level library
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_ABINIT_COMMON_INIT], [
  # Init
  sd_abinit_common_cppflags=""
  sd_abinit_common_cflags=""
  sd_abinit_common_cxxflags=""
  sd_abinit_common_fcflags=""
  sd_abinit_common_ldflags=""
  sd_abinit_common_libs=""
  sd_abinit_common_enable=""
  sd_abinit_common_init="unknown"
  sd_abinit_common_ok="unknown"

  # Set adjustable parameters
  sd_abinit_common_options="$1"
  sd_abinit_common_libs_def="$2"
  sd_abinit_common_cppflags_def="$3"
  sd_abinit_common_cflags_def="$4"
  sd_abinit_common_cxxflags_def="$5"
  sd_abinit_common_fcflags_def="$6"
  sd_abinit_common_ldflags_def="$7"
  sd_abinit_common_enable_def=""
  sd_abinit_common_policy=""
  sd_abinit_common_status=""

  # Process options
  for kwd in ${sd_abinit_common_options}; do
    case "${kwd}" in
      auto)
        sd_abinit_common_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_abinit_common_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_abinit_common_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg ABINIT Common option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_abinit_common_libs_def}" && sd_abinit_common_libs_def="-labinit_common"
  test -z "${sd_abinit_common_policy}" && sd_abinit_common_policy="fail"
  test -z "${sd_abinit_common_status}" && sd_abinit_common_status="optional"
  test -z "${sd_abinit_common_enable_def}" && sd_abinit_common_enable_def="no"
  case "${sd_abinit_common_status}" in
    implicit|required)
      sd_abinit_common_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([abinit_common],
    [AS_HELP_STRING([--with-abinit_common],
      [Install prefix of the ABINIT Common library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_abinit_common_enable="${withval}"
        sd_abinit_common_init="yon"
      else
        sd_abinit_common_enable="yes"
        sd_abinit_common_init="dir"
      fi],
    [ sd_abinit_common_enable="${sd_abinit_common_enable_def}"; sd_abinit_common_init="def"])

  # Declare environment variables
  AC_ARG_VAR([ABINIT_COMMON_CPPFLAGS], [C preprocessing flags for ABINIT Common.])
  AC_ARG_VAR([ABINIT_COMMON_FCFLAGS], [Fortran flags for ABINIT Common.])
  AC_ARG_VAR([ABINIT_COMMON_LDFLAGS], [Linker flags for ABINIT Common.])
  AC_ARG_VAR([ABINIT_COMMON_LIBS], [Library flags for ABINIT Common.])

  # Detect use of environment variables
  if test "${sd_abinit_common_enable}" = "yes" -o "${sd_abinit_common_enable}" = "auto"; then
    tmp_abinit_common_vars="${ABINIT_COMMON_CPPFLAGS}${ABINIT_COMMON_FCFLAGS}${ABINIT_COMMON_LDFLAGS}${ABINIT_COMMON_LIBS}"
    if test "${sd_abinit_common_init}" = "def" -a ! -z "${tmp_abinit_common_vars}"; then
      sd_abinit_common_enable="yes"
      sd_abinit_common_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_ABINIT_COMMON_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_abinit_common_enable}" = "yes" -o "${sd_abinit_common_enable}" = "auto"; then

    # Set ABINIT Common-specific flags
    case "${sd_abinit_common_init}" in

      def|yon)
        sd_abinit_common_cppflags="${sd_abinit_common_cppflags_def}"
        sd_abinit_common_fcflags="${sd_abinit_common_fcflags_def}"
        sd_abinit_common_ldflags="${sd_abinit_common_ldflags_def}"
        sd_abinit_common_libs="${sd_abinit_common_libs_def}"
        ;;

      dir)
        sd_abinit_common_cppflags="-I${with_abinit_common}/include"
        sd_abinit_common_fcflags="${sd_abinit_common_fcflags_def} -I${with_abinit_common}/include"
        sd_abinit_common_ldflags="${sd_abinit_common_ldflags_def}"
        sd_abinit_common_libs="-L${with_abinit_common}/lib ${sd_abinit_common_libs_def}"
        ;;

      env)
        sd_abinit_common_cppflags="${sd_abinit_common_cppflags_def}"
        sd_abinit_common_fcflags="${sd_abinit_common_fcflags_def}"
        sd_abinit_common_ldflags="${sd_abinit_common_ldflags_def}"
        sd_abinit_common_libs="${sd_abinit_common_libs_def}"
        test ! -z "${ABINIT_COMMON_CPPFLAGS}" && sd_abinit_common_cppflags="${ABINIT_COMMON_CPPFLAGS}"
        test ! -z "${ABINIT_COMMON_FCFLAGS}" && sd_abinit_common_fcflags="${ABINIT_COMMON_FCFLAGS}"
        test ! -z "${ABINIT_COMMON_LDFLAGS}" && sd_abinit_common_ldflags="${ABINIT_COMMON_LDFLAGS}"
        test ! -z "${ABINIT_COMMON_LIBS}" && sd_abinit_common_libs="${ABINIT_COMMON_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for ABINIT Common: '${sd_abinit_common_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_abinit_common_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_abinit_common_init="esl"
    sd_abinit_common_cppflags=""
    sd_abinit_common_fcflags=""
    sd_abinit_common_ldflags=""
    sd_abinit_common_libs=""
  fi

  # Display configuration
  _SD_ABINIT_COMMON_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_abinit_common_options)
  AC_SUBST(sd_abinit_common_enable_def)
  AC_SUBST(sd_abinit_common_policy)
  AC_SUBST(sd_abinit_common_status)
  AC_SUBST(sd_abinit_common_enable)
  AC_SUBST(sd_abinit_common_init)
  AC_SUBST(sd_abinit_common_ok)
  AC_SUBST(sd_abinit_common_cppflags)
  AC_SUBST(sd_abinit_common_fcflags)
  AC_SUBST(sd_abinit_common_ldflags)
  AC_SUBST(sd_abinit_common_libs)
  AC_SUBST(with_abinit_common)

  # Clean-up
  unset tmp_abinit_common_vars
]) # SD_ABINIT_COMMON_INIT


AC_DEFUN([SD_ABINIT_COMMON_DETECT], [
  # Display configuration
  _SD_ABINIT_COMMON_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_abinit_common_enable}" = "auto" -o "${sd_abinit_common_enable}" = "yes"; then
    _SD_ABINIT_COMMON_CHECK_USE

    if test "${sd_abinit_common_ok}" = "yes"; then
      if test "${sd_abinit_common_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_abinit_common_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_abinit_common_fcflags}"
        LIBS="${sd_abinit_common_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_abinit_common_ldflags}"

      AC_DEFINE([HAVE_ABINIT_COMMON], 1,
        [Define to 1 if you have the ABINIT Common library.])
    else
      if test "${sd_abinit_common_status}" = "optional" -a \
              "${sd_abinit_common_init}" = "def"; then
        sd_abinit_common_enable="no"
        sd_abinit_common_cppflags=""
        sd_abinit_common_fcflags=""
        sd_abinit_common_ldflags=""
        sd_abinit_common_libs=""
      else
        AC_MSG_FAILURE([invalid ABINIT Common configuration])
      fi
    fi
  else
    sd_abinit_common_enable="no"
    sd_abinit_common_cppflags=""
    sd_abinit_common_fcflags=""
    sd_abinit_common_ldflags=""
    sd_abinit_common_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_ABINIT_COMMON_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_abinit_common_init}" = "esl"; then
    AC_MSG_NOTICE([will look for ABINIT Common in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_abinit_common_libs_def}])
  else
    # FIXME: dirty hack to get the detection work in ABINIT (YP)
    sd_linalg_libs="${abi_linalg_libs}"
    CPPFLAGS="${CPPFLAGS} ${sd_abinit_common_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_abinit_common_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_abinit_common_ldflags}"
    LIBS="${sd_abinit_common_libs} ${sd_linalg_libs} ${LIBS}"
  fi

  # Check ABINIT Common API
  AC_MSG_CHECKING([whether the ABINIT Common Fortran interface works])
  for tmp_incs in "" "-I/usr/include"; do
    FCFLAGS="${FCFLAGS} ${tmp_incs}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use defs_basis
      ]])], [sd_abinit_common_ok="yes"], [sd_abinit_common_ok="no"])
    AC_LANG_POP([Fortran])
    if test "${sd_abinit_common_ok}" = "yes"; then
      test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
      break
    fi
  done
  AC_MSG_RESULT([${sd_abinit_common_ok}])
  unset tmp_incs

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_ABINIT_COMMON_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_ABINIT_COMMON_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_abinit_common_invalid="no"
  if test "${sd_abinit_common_enable_def}" != "auto" -a \
          "${sd_abinit_common_enable_def}" != "no" -a \
          "${sd_abinit_common_enable_def}" != "yes"; then
    case "${sd_abinit_common_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_abinit_common_enable_def = '${sd_abinit_common_enable_def}'])
        ;;
      skip)
        tmp_abinit_common_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_abinit_common_enable_def = '${sd_abinit_common_enable_def}'])
        tmp_abinit_common_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_abinit_common_invalid}" = "yes"; then
    if test "${sd_abinit_common_status}" = "required"; then
      sd_abinit_common_enable_def="yes"
    else
      sd_abinit_common_enable_def="no"
    fi
    tmp_abinit_common_invalid="no"
    AC_MSG_NOTICE([setting sd_abinit_common_enable_def to '${sd_abinit_common_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_abinit_common_invalid="no"
  if test "${sd_abinit_common_status}" = "implicit" -o \
          "${sd_abinit_common_status}" = "required"; then
    if test "${sd_abinit_common_enable}" = "no"; then
      case "${sd_abinit_common_policy}" in
        fail)
          AC_MSG_ERROR([The ABINIT Common package is required and cannot be disabled
                  See https://launchpad.net/abinit_common for details on how to
                  install it.])
          ;;
        skip)
          tmp_abinit_common_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The ABINIT Common package is required and cannot be disabled])
          tmp_abinit_common_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_abinit_common_enable}" = "auto"; then
      AC_MSG_NOTICE([setting ABINIT Common trigger to yes])
      sd_abinit_common_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_abinit_common_invalid}" = "yes"; then
    case "${sd_abinit_common_status}" in
      implicit|required)
        sd_abinit_common_enable="yes"
        ;;
      optional)
        sd_abinit_common_enable="no"
        ;;
    esac
    tmp_abinit_common_invalid="no"
    AC_MSG_NOTICE([setting sd_abinit_common_enable to '${sd_abinit_common_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_abinit_common_vars="${ABINIT_COMMON_CPPFLAGS}${ABINIT_COMMON_FCFLAGS}${ABINIT_COMMON_LDFLAGS}${ABINIT_COMMON_LIBS}"
  tmp_abinit_common_invalid="no"
  if test ! -z "${tmp_abinit_common_vars}" -a ! -z "${with_abinit_common}"; then
    case "${sd_abinit_common_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for ABINIT Common
                  Please use ABINIT_COMMON_FCFLAGS + ABINIT_COMMON_LIBS or --with-abinit_common,
                  not both.])
        ;;
      skip)
        tmp_abinit_common_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for ABINIT Common])
        tmp_abinit_common_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_abinit_common_vars}"; then
    sd_abinit_common_enable="yes"
    sd_abinit_common_init="env"
    if test "${tmp_abinit_common_invalid}" = "yes"; then
      tmp_abinit_common_invalid="no"
      AC_MSG_NOTICE([overriding --with-abinit_common with ABINIT_COMMON_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_abinit_common_status}" = "implicit"; then
    if test "${sd_abinit_common_fcflags}" != ""; then
      sd_abinit_common_fcflags=""
      AC_MSG_NOTICE([resetting ABINIT Common Fortran flags (implicit package)])
    fi
    if test "${sd_abinit_common_ldflags}" != ""; then
      sd_abinit_common_ldflags=""
      AC_MSG_NOTICE([resetting ABINIT Common linker flags (implicit package)])
    fi
    if test "${sd_abinit_common_libs}" != ""; then
      sd_abinit_common_libs=""
      AC_MSG_NOTICE([resetting ABINIT Common library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_abinit_common_enable}" = "implicit"; then
    sd_abinit_common_fcflags=""
    sd_abinit_common_ldflags=""
    sd_abinit_common_libs=""
  fi

  # Clean-up
  unset tmp_abinit_common_invalid
  unset tmp_abinit_common_vars
]) # _SD_ABINIT_COMMON_CHECK_CONFIG


AC_DEFUN([_SD_ABINIT_COMMON_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable ABINIT Common])
  AC_MSG_RESULT([${sd_abinit_common_enable}])
  if test "${sd_abinit_common_enable}" != "no"; then
    AC_MSG_CHECKING([how ABINIT Common parameters have been set])
    AC_MSG_RESULT([${sd_abinit_common_init}])
    AC_MSG_CHECKING([for ABINIT Common C preprocessing flags])
    if test "${sd_abinit_common_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_abinit_common_cppflags}])
    fi
    AC_MSG_CHECKING([for ABINIT Common Fortran flags])
    if test "${sd_abinit_common_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_abinit_common_fcflags}])
    fi
    AC_MSG_CHECKING([for ABINIT Common linker flags])
    if test "${sd_abinit_common_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_abinit_common_ldflags}])
    fi
    AC_MSG_CHECKING([for ABINIT Common library flags])
    if test "${sd_abinit_common_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_abinit_common_libs}])
    fi
  fi
]) # _SD_ABINIT_COMMON_DUMP_CONFIG

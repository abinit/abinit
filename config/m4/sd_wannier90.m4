## Copyright (C) 2019 Yann Pouillon

#
# Maximally-Localized Wannier Functions library (Wannier90)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_WANNIER90_INIT], [
  # Init
  sd_wannier90_cppflags=""
  sd_wannier90_fcflags=""
  sd_wannier90_ldflags=""
  sd_wannier90_libs=""
  sd_wannier90_enable=""
  sd_wannier90_init="unknown"
  sd_wannier90_ok="unknown"

  # Set adjustable parameters
  sd_wannier90_options="$1"
  sd_wannier90_libs_def="$2"
  sd_wannier90_cppflags_def="$3"
  sd_wannier90_cflags_def="$4"
  sd_wannier90_cxxflags_def="$5"
  sd_wannier90_fcflags_def="$6"
  sd_wannier90_ldflags_def="$7"
  sd_wannier90_enable_def=""
  sd_wannier90_policy=""
  sd_wannier90_status=""

  # Process options
  for kwd in ${sd_wannier90_options}; do
    case "${kwd}" in
      auto|no|yes)
        sd_wannier90_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_wannier90_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_wannier90_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg Wannier90 option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_wannier90_enable_def}" && sd_wannier90_enable_def="auto"
  test -z "${sd_wannier90_policy}" && sd_wannier90_policy="fail"
  test -z "${sd_wannier90_status}" && sd_wannier90_status="optional"
  test -z "${sd_wannier90_libs_def}" && sd_wannier90_libs_def="-lwannier"

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([wannier90],
    [AS_HELP_STRING([--with-wannier90],
      [Install prefix of the PSML I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_wannier90_enable="${withval}"
        sd_wannier90_init="yon"
      else
        sd_wannier90_enable="yes"
        sd_wannier90_init="dir"
      fi],
    [ sd_wannier90_enable="${sd_wannier90_enable_def}"; sd_wannier90_init="def"])

  # Declare environment variables
  AC_ARG_VAR([WANNIER90_CPPFLAGS], [C preprocessing flags for Wannier90.])
  AC_ARG_VAR([WANNIER90_FCFLAGS], [Fortran flags for Wannier90.])
  AC_ARG_VAR([WANNIER90_LDFLAGS], [Linker flags for Wannier90.])
  AC_ARG_VAR([WANNIER90_LIBS], [Library flags for Wannier90.])

  # Detect use of environment variables
  if test "${sd_wannier90_enable}" = "yes" -o "${sd_wannier90_enable}" = "auto"; then
    tmp_wannier90_vars="${WANNIER90_CPPFLAGS}${WANNIER90_FCFLAGS}${WANNIER90_LDFLAGS}${WANNIER90_LIBS}"
    if test "${sd_wannier90_init}" = "def" -a ! -z "${tmp_wannier90_vars}"; then
      sd_wannier90_enable="yes"
      sd_wannier90_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_WANNIER90_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_wannier90_enable}" = "yes" -o "${sd_wannier90_enable}" = "auto"; then

    # Set Wannier90-specific flags
    case "${sd_wannier90_init}" in

      def|yon)
        sd_wannier90_cppflags="${sd_wannier90_cppflags_def}"
        sd_wannier90_fcflags="${sd_wannier90_fcflags_def}"
        sd_wannier90_ldflags="${sd_wannier90_ldflags_def}"
        sd_wannier90_libs="${sd_wannier90_libs_def}"
        ;;

      dir)
        sd_wannier90_cppflags="-I${with_wannier90}/include"
        sd_wannier90_fcflags="${sd_wannier90_fcflags_def} -I${with_wannier90}/include"
        sd_wannier90_ldflags="${sd_wannier90_ldflags_def}"
        sd_wannier90_libs="-L${with_wannier90}/lib ${sd_wannier90_libs_def}"
        ;;

      env)
        sd_wannier90_cppflags="${sd_wannier90_cppflags_def}"
        sd_wannier90_fcflags="${sd_wannier90_fcflags_def}"
        sd_wannier90_ldflags="${sd_wannier90_ldflags_def}"
        sd_wannier90_libs="${sd_wannier90_libs_def}"
        test ! -z "${WANNIER90_CPPFLAGS}" && sd_wannier90_cppflags="${WANNIER90_CPPFLAGS}"
        test ! -z "${WANNIER90_FCFLAGS}" && sd_wannier90_fcflags="${WANNIER90_FCFLAGS}"
        test ! -z "${WANNIER90_LDFLAGS}" && sd_wannier90_ldflags="${WANNIER90_LDFLAGS}"
        test ! -z "${WANNIER90_LIBS}" && sd_wannier90_libs="${WANNIER90_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for Wannier90: '${sd_wannier90_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_wannier90_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_wannier90_init="esl"
    sd_wannier90_cppflags=""
    sd_wannier90_fcflags=""
    sd_wannier90_ldflags=""
    sd_wannier90_libs=""
  fi

  # Display configuration
  _SD_WANNIER90_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_wannier90_options)
  AC_SUBST(sd_wannier90_enable_def)
  AC_SUBST(sd_wannier90_policy)
  AC_SUBST(sd_wannier90_status)
  AC_SUBST(sd_wannier90_enable)
  AC_SUBST(sd_wannier90_init)
  AC_SUBST(sd_wannier90_ok)
  AC_SUBST(sd_wannier90_cppflags)
  AC_SUBST(sd_wannier90_fcflags)
  AC_SUBST(sd_wannier90_ldflags)
  AC_SUBST(sd_wannier90_libs)
  AC_SUBST(with_wannier90)

  # Clean-up
  unset tmp_wannier90_vars
]) # SD_WANNIER90_INIT


AC_DEFUN([SD_WANNIER90_DETECT], [
  # Display configuration
  _SD_WANNIER90_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_wannier90_enable}" = "auto" -o "${sd_wannier90_enable}" = "yes"; then
    _SD_WANNIER90_CHECK_USE

    if test "${sd_wannier90_ok}" = "yes"; then
      if test "${sd_wannier90_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_wannier90_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_wannier90_fcflags}"
        LIBS="${sd_wannier90_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_wannier90_ldflags}"

      AC_DEFINE([HAVE_WANNIER90], 1,
        [Define to 1 if you have the Wannier90 library.])
    else
      if test "${sd_wannier90_status}" = "optional" -a \
              "${sd_wannier90_init}" = "def"; then
        sd_wannier90_enable="no"
        sd_wannier90_cppflags=""
        sd_wannier90_fcflags=""
        sd_wannier90_ldflags=""
        sd_wannier90_libs=""
      else
        AC_MSG_FAILURE([invalid Wannier90 configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_WANNIER90_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_wannier90_init}" = "esl"; then
    AC_MSG_NOTICE([will look for Wannier90 in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_wannier90_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_wannier90_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_wannier90_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_wannier90_ldflags}"
    LIBS="${sd_wannier90_libs} ${LIBS}"
  fi

  # Check Wannier90 API
  AC_MSG_CHECKING([whether the Wannier90 Fortran interface works])
  for tmp_incs in "" "-I/usr/include"; do
    FCFLAGS="${FCFLAGS} ${tmp_incs}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call wannier_run
      ]])], [sd_wannier90_ok="yes"], [sd_wannier90_ok="no"])
    AC_LANG_POP([Fortran])
    if test "${sd_wannier90_ok}" = "yes"; then
      test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
      break
    fi
  done
  AC_MSG_RESULT([${sd_wannier90_ok}])
  unset tmp_incs

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_WANNIER90_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_WANNIER90_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_wannier90_invalid="no"
  if test "${sd_wannier90_enable_def}" != "auto" -a \
          "${sd_wannier90_enable_def}" != "no" -a \
          "${sd_wannier90_enable_def}" != "yes"; then
    case "${sd_wannier90_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_wannier90_enable_def = '${sd_wannier90_enable_def}'])
        ;;
      skip)
        tmp_wannier90_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_wannier90_enable_def = '${sd_wannier90_enable_def}'])
        tmp_wannier90_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_wannier90_invalid}" = "yes"; then
    if test "${sd_wannier90_status}" = "required"; then
      sd_wannier90_enable_def="yes"
    else
      sd_wannier90_enable_def="no"
    fi
    tmp_wannier90_invalid="no"
    AC_MSG_NOTICE([setting sd_wannier90_enable_def to '${sd_wannier90_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_wannier90_invalid="no"
  if test "${sd_wannier90_status}" = "implicit" -o \
          "${sd_wannier90_status}" = "required"; then
    if test "${sd_wannier90_enable}" = "no"; then
      case "${sd_wannier90_policy}" in
        fail)
          AC_MSG_ERROR([The Wannier90 package is required and cannot be disabled
                  See https://launchpad.net/wannier90 for details on how to
                  install it.])
          ;;
        skip)
          tmp_wannier90_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The Wannier90 package is required and cannot be disabled])
          tmp_wannier90_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_wannier90_enable}" = "auto"; then
      AC_MSG_NOTICE([setting Wannier90 trigger to yes])
      sd_wannier90_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_wannier90_invalid}" = "yes"; then
    case "${sd_wannier90_status}" in
      implicit|required)
        sd_wannier90_enable="yes"
        ;;
      optional)
        sd_wannier90_enable="no"
        ;;
    esac
    tmp_wannier90_invalid="no"
    AC_MSG_NOTICE([setting sd_wannier90_enable to '${sd_wannier90_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_wannier90_vars="${WANNIER90_FCFLAGS}${WANNIER90_LDFLAGS}${WANNIER90_LIBS}"
  tmp_wannier90_invalid="no"
  if test ! -z "${tmp_wannier90_vars}" -a ! -z "${with_wannier90}"; then
    case "${sd_wannier90_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for Wannier90
                  Please use WANNIER90_FCFLAGS + WANNIER90_LIBS or --with-wannier90,
                  not both.])
        ;;
      skip)
        tmp_wannier90_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for Wannier90])
        tmp_wannier90_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_wannier90_vars}"; then
    sd_wannier90_enable="yes"
    sd_wannier90_init="env"
    if test "${tmp_wannier90_invalid}" = "yes"; then
      tmp_wannier90_invalid="no"
      AC_MSG_NOTICE([overriding --with-wannier90 with WANNIER90_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_wannier90_status}" = "implicit"; then
    if test "${sd_wannier90_fcflags}" != ""; then
      sd_wannier90_fcflags=""
      AC_MSG_NOTICE([resetting Wannier90 Fortran flags (implicit package)])
    fi
    if test "${sd_wannier90_ldflags}" != ""; then
      sd_wannier90_ldflags=""
      AC_MSG_NOTICE([resetting Wannier90 linker flags (implicit package)])
    fi
    if test "${sd_wannier90_libs}" != ""; then
      sd_wannier90_libs=""
      AC_MSG_NOTICE([resetting Wannier90 library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_wannier90_enable}" = "implicit"; then
    sd_wannier90_fcflags=""
    sd_wannier90_ldflags=""
    sd_wannier90_libs=""
  fi

  # Clean-up
  unset tmp_wannier90_invalid
  unset tmp_wannier90_vars
]) # _SD_WANNIER90_CHECK_CONFIG


AC_DEFUN([_SD_WANNIER90_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable Wannier90])
  AC_MSG_RESULT([${sd_wannier90_enable}])
  if test "${sd_wannier90_enable}" != "no"; then
    AC_MSG_CHECKING([how Wannier90 parameters have been set])
    AC_MSG_RESULT([${sd_wannier90_init}])
    AC_MSG_CHECKING([for Wannier90 C preprocessing flags])
    if test "${sd_wannier90_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_wannier90_cppflags}])
    fi
    AC_MSG_CHECKING([for Wannier90 Fortran flags])
    if test "${sd_wannier90_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_wannier90_fcflags}])
    fi
    AC_MSG_CHECKING([for Wannier90 linker flags])
    if test "${sd_wannier90_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_wannier90_ldflags}])
    fi
    AC_MSG_CHECKING([for Wannier90 library flags])
    if test "${sd_wannier90_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_wannier90_libs}])
    fi
  fi
]) # _SD_WANNIER90_DUMP_CONFIG

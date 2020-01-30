## Copyright (C) 2019 Yann Pouillon

#
# LibPAW Projector-Augmented Waves library
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_LIBPAW_INIT], [
  # Init
  sd_libpaw_cppflags=""
  sd_libpaw_cflags=""
  sd_libpaw_cxxflags=""
  sd_libpaw_fcflags=""
  sd_libpaw_ldflags=""
  sd_libpaw_libs=""
  sd_libpaw_enable=""
  sd_libpaw_init="unknown"
  sd_libpaw_ok="unknown"

  # Set adjustable parameters
  sd_libpaw_options="$1"
  sd_libpaw_libs_def="$2"
  sd_libpaw_cppflags_def="$3"
  sd_libpaw_cflags_def="$4"
  sd_libpaw_cxxflags_def="$5"
  sd_libpaw_fcflags_def="$6"
  sd_libpaw_ldflags_def="$7"
  sd_libpaw_enable_def=""
  sd_libpaw_policy=""
  sd_libpaw_status=""

  # Process options
  for kwd in ${sd_libpaw_options}; do
    case "${kwd}" in
      auto)
        sd_libpaw_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_libpaw_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_libpaw_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg LibPAW option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_libpaw_libs_def}" && sd_libpaw_libs_def="-lpaw"
  test -z "${sd_libpaw_policy}" && sd_libpaw_policy="fail"
  test -z "${sd_libpaw_status}" && sd_libpaw_status="optional"
  test -z "${sd_libpaw_enable_def}" && sd_libpaw_enable_def="no"
  case "${sd_libpaw_status}" in
    implicit|required)
      sd_libpaw_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([libpaw],
    [AS_HELP_STRING([--with-libpaw],
      [Install prefix of the LibPAW library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_libpaw_enable="${withval}"
        sd_libpaw_init="yon"
      else
        sd_libpaw_enable="yes"
        sd_libpaw_init="dir"
      fi],
    [ sd_libpaw_enable="${sd_libpaw_enable_def}"; sd_libpaw_init="def"])

  # Declare environment variables
  AC_ARG_VAR([LIBPAW_CPPFLAGS], [C preprocessing flags for LibPAW.])
  AC_ARG_VAR([LIBPAW_FCFLAGS], [Fortran flags for LibPAW.])
  AC_ARG_VAR([LIBPAW_LDFLAGS], [Linker flags for LibPAW.])
  AC_ARG_VAR([LIBPAW_LIBS], [Library flags for LibPAW.])

  # Detect use of environment variables
  if test "${sd_libpaw_enable}" = "yes" -o "${sd_libpaw_enable}" = "auto"; then
    tmp_libpaw_vars="${LIBPAW_CPPFLAGS}${LIBPAW_FCFLAGS}${LIBPAW_LDFLAGS}${LIBPAW_LIBS}"
    if test "${sd_libpaw_init}" = "def" -a ! -z "${tmp_libpaw_vars}"; then
      sd_libpaw_enable="yes"
      sd_libpaw_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_LIBPAW_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_libpaw_enable}" = "yes" -o "${sd_libpaw_enable}" = "auto"; then

    # Set LibPAW-specific flags
    case "${sd_libpaw_init}" in

      def|yon)
        sd_libpaw_cppflags="${sd_libpaw_cppflags_def}"
        sd_libpaw_fcflags="${sd_libpaw_fcflags_def}"
        sd_libpaw_ldflags="${sd_libpaw_ldflags_def}"
        sd_libpaw_libs="${sd_libpaw_libs_def}"
        ;;

      dir)
        sd_libpaw_cppflags="-I${with_libpaw}/include"
        sd_libpaw_fcflags="${sd_libpaw_fcflags_def} -I${with_libpaw}/include"
        sd_libpaw_ldflags="${sd_libpaw_ldflags_def}"
        sd_libpaw_libs="-L${with_libpaw}/lib ${sd_libpaw_libs_def}"
        ;;

      env)
        sd_libpaw_cppflags="${sd_libpaw_cppflags_def}"
        sd_libpaw_fcflags="${sd_libpaw_fcflags_def}"
        sd_libpaw_ldflags="${sd_libpaw_ldflags_def}"
        sd_libpaw_libs="${sd_libpaw_libs_def}"
        test ! -z "${LIBPAW_CPPFLAGS}" && sd_libpaw_cppflags="${LIBPAW_CPPFLAGS}"
        test ! -z "${LIBPAW_FCFLAGS}" && sd_libpaw_fcflags="${LIBPAW_FCFLAGS}"
        test ! -z "${LIBPAW_LDFLAGS}" && sd_libpaw_ldflags="${LIBPAW_LDFLAGS}"
        test ! -z "${LIBPAW_LIBS}" && sd_libpaw_libs="${LIBPAW_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for LibPAW: '${sd_libpaw_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_libpaw_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_libpaw_init="esl"
    sd_libpaw_cppflags=""
    sd_libpaw_fcflags=""
    sd_libpaw_ldflags=""
    sd_libpaw_libs=""
  fi

  # Display configuration
  _SD_LIBPAW_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_libpaw_options)
  AC_SUBST(sd_libpaw_enable_def)
  AC_SUBST(sd_libpaw_policy)
  AC_SUBST(sd_libpaw_status)
  AC_SUBST(sd_libpaw_enable)
  AC_SUBST(sd_libpaw_init)
  AC_SUBST(sd_libpaw_ok)
  AC_SUBST(sd_libpaw_cppflags)
  AC_SUBST(sd_libpaw_fcflags)
  AC_SUBST(sd_libpaw_ldflags)
  AC_SUBST(sd_libpaw_libs)
  AC_SUBST(with_libpaw)

  # Clean-up
  unset tmp_libpaw_vars
]) # SD_LIBPAW_INIT


AC_DEFUN([SD_LIBPAW_DETECT], [
  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_libpaw_enable}" = "auto" -o "${sd_libpaw_enable}" = "yes"; then
    _SD_LIBPAW_CHECK_USE

    if test "${sd_libpaw_ok}" = "yes"; then
      if test "${sd_libpaw_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_libpaw_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_libpaw_fcflags}"
        LIBS="${sd_libpaw_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_libpaw_ldflags}"

      AC_DEFINE([HAVE_LIBPAW], 1,
        [Define to 1 if you have the LibPAW library.])

      _SD_LIBPAW_DUMP_CONFIG
    else
      if test "${sd_libpaw_status}" = "optional" -a \
              "${sd_libpaw_init}" = "def"; then
        sd_libpaw_enable="no"
        sd_libpaw_cppflags=""
        sd_libpaw_fcflags=""
        sd_libpaw_ldflags=""
        sd_libpaw_libs=""
      else
        AC_MSG_FAILURE([invalid LibPAW configuration])
      fi
    fi
  else
    sd_libpaw_enable="no"
    sd_libpaw_cppflags=""
    sd_libpaw_fcflags=""
    sd_libpaw_ldflags=""
    sd_libpaw_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_LIBPAW_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_libpaw_init}" = "esl"; then
    AC_MSG_NOTICE([will look for LibPAW in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_libpaw_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_libpaw_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_libpaw_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_libpaw_ldflags}"
    LIBS="${sd_libpaw_libs} ${sd_linalg_libs} ${LIBS}"
  fi

  # Check LibPAW API
  AC_MSG_CHECKING([whether the LibPAW Fortran interface works])
  for tmp_incs in "" "-I/usr/include"; do
    FCFLAGS="${FCFLAGS} ${tmp_incs}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use m_libpaw_defs
      ]])], [sd_libpaw_ok="yes"], [sd_libpaw_ok="no"])
    AC_LANG_POP([Fortran])
    if test "${sd_libpaw_ok}" = "yes"; then
      test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
      break
    fi
  done
  AC_MSG_RESULT([${sd_libpaw_ok}])
  unset tmp_incs

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_LIBPAW_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_LIBPAW_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_libpaw_invalid="no"
  if test "${sd_libpaw_enable_def}" != "auto" -a \
          "${sd_libpaw_enable_def}" != "no" -a \
          "${sd_libpaw_enable_def}" != "yes"; then
    case "${sd_libpaw_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_libpaw_enable_def = '${sd_libpaw_enable_def}'])
        ;;
      skip)
        tmp_libpaw_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_libpaw_enable_def = '${sd_libpaw_enable_def}'])
        tmp_libpaw_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_libpaw_invalid}" = "yes"; then
    if test "${sd_libpaw_status}" = "required"; then
      sd_libpaw_enable_def="yes"
    else
      sd_libpaw_enable_def="no"
    fi
    tmp_libpaw_invalid="no"
    AC_MSG_NOTICE([setting sd_libpaw_enable_def to '${sd_libpaw_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_libpaw_invalid="no"
  if test "${sd_libpaw_status}" = "implicit" -o \
          "${sd_libpaw_status}" = "required"; then
    if test "${sd_libpaw_enable}" = "no"; then
      case "${sd_libpaw_policy}" in
        fail)
          AC_MSG_ERROR([The LibPAW package is required and cannot be disabled
                  See https://launchpad.net/libpaw for details on how to
                  install it.])
          ;;
        skip)
          tmp_libpaw_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The LibPAW package is required and cannot be disabled])
          tmp_libpaw_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_libpaw_enable}" = "auto"; then
      AC_MSG_NOTICE([setting LibPAW trigger to yes])
      sd_libpaw_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_libpaw_invalid}" = "yes"; then
    case "${sd_libpaw_status}" in
      implicit|required)
        sd_libpaw_enable="yes"
        ;;
      optional)
        sd_libpaw_enable="no"
        ;;
    esac
    tmp_libpaw_invalid="no"
    AC_MSG_NOTICE([setting sd_libpaw_enable to '${sd_libpaw_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_libpaw_vars="${LIBPAW_CPPFLAGS}${LIBPAW_FCFLAGS}${LIBPAW_LDFLAGS}${LIBPAW_LIBS}"
  tmp_libpaw_invalid="no"
  if test ! -z "${tmp_libpaw_vars}" -a ! -z "${with_libpaw}"; then
    case "${sd_libpaw_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for LibPAW
                  Please use LIBPAW_FCFLAGS + LIBPAW_LIBS or --with-libpaw,
                  not both.])
        ;;
      skip)
        tmp_libpaw_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for LibPAW])
        tmp_libpaw_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_libpaw_vars}"; then
    sd_libpaw_enable="yes"
    sd_libpaw_init="env"
    if test "${tmp_libpaw_invalid}" = "yes"; then
      tmp_libpaw_invalid="no"
      AC_MSG_NOTICE([overriding --with-libpaw with LIBPAW_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_libpaw_status}" = "implicit"; then
    if test "${sd_libpaw_fcflags}" != ""; then
      sd_libpaw_fcflags=""
      AC_MSG_NOTICE([resetting LibPAW Fortran flags (implicit package)])
    fi
    if test "${sd_libpaw_ldflags}" != ""; then
      sd_libpaw_ldflags=""
      AC_MSG_NOTICE([resetting LibPAW linker flags (implicit package)])
    fi
    if test "${sd_libpaw_libs}" != ""; then
      sd_libpaw_libs=""
      AC_MSG_NOTICE([resetting LibPAW library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_libpaw_enable}" = "implicit"; then
    sd_libpaw_fcflags=""
    sd_libpaw_ldflags=""
    sd_libpaw_libs=""
  fi

  # Clean-up
  unset tmp_libpaw_invalid
  unset tmp_libpaw_vars
]) # _SD_LIBPAW_CHECK_CONFIG


AC_DEFUN([_SD_LIBPAW_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable LibPAW])
  AC_MSG_RESULT([${sd_libpaw_enable}])
  if test "${sd_libpaw_enable}" != "no"; then
    AC_MSG_CHECKING([how LibPAW parameters have been set])
    AC_MSG_RESULT([${sd_libpaw_init}])
    AC_MSG_CHECKING([for LibPAW C preprocessing flags])
    if test "${sd_libpaw_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpaw_cppflags}])
    fi
    AC_MSG_CHECKING([for LibPAW Fortran flags])
    if test "${sd_libpaw_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpaw_fcflags}])
    fi
    AC_MSG_CHECKING([for LibPAW linker flags])
    if test "${sd_libpaw_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpaw_ldflags}])
    fi
    AC_MSG_CHECKING([for LibPAW library flags])
    if test "${sd_libpaw_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpaw_libs}])
    fi
  fi
]) # _SD_LIBPAW_DUMP_CONFIG

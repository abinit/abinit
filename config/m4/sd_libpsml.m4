## Copyright (C) 2019 Yann Pouillon

#
# PSeudopotential Markup Language I/O library (LibPSML)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_LIBPSML_INIT], [
  # Init
  sd_libpsml_cppflags=""
  sd_libpsml_cflags=""
  sd_libpsml_cxxflags=""
  sd_libpsml_fcflags=""
  sd_libpsml_ldflags=""
  sd_libpsml_libs=""
  sd_libpsml_enable=""
  sd_libpsml_init="unknown"
  sd_libpsml_ok="unknown"

  # Set adjustable parameters
  sd_libpsml_options="$1"
  sd_libpsml_libs_def="$2"
  sd_libpsml_cppflags_def="$3"
  sd_libpsml_cflags_def="$4"
  sd_libpsml_cxxflags_def="$5"
  sd_libpsml_fcflags_def="$6"
  sd_libpsml_ldflags_def="$7"
  sd_libpsml_enable_def=""
  sd_libpsml_policy=""
  sd_libpsml_status=""

  # Process options
  for kwd in ${sd_libpsml_options}; do
    case "${kwd}" in
      auto)
        sd_libpsml_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_libpsml_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_libpsml_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg LibPSML option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_libpsml_libs_def}" && sd_libpsml_libs_def="-lpsml"
  test -z "${sd_libpsml_policy}" && sd_libpsml_policy="fail"
  test -z "${sd_libpsml_status}" && sd_libpsml_status="optional"
  test -z "${sd_libpsml_enable_def}" && sd_libpsml_enable_def="no"
  case "${sd_libpsml_status}" in
    implicit|required)
      sd_libpsml_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([libpsml],
    [AS_HELP_STRING([--with-libpsml],
      [Install prefix of the PSML I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_libpsml_enable="${withval}"
        sd_libpsml_init="yon"
      else
        sd_libpsml_enable="yes"
        sd_libpsml_init="dir"
      fi],
    [ sd_libpsml_enable="${sd_libpsml_enable_def}"; sd_libpsml_init="def"])

  # Declare environment variables
  AC_ARG_VAR([LIBPSML_CPPFLAGS], [C preprocessing flags for LibPSML.])
  AC_ARG_VAR([LIBPSML_FCFLAGS], [Fortran flags for LibPSML.])
  AC_ARG_VAR([LIBPSML_LDFLAGS], [Linker flags for LibPSML.])
  AC_ARG_VAR([LIBPSML_LIBS], [Library flags for LibPSML.])

  # Detect use of environment variables
  if test "${sd_libpsml_enable}" = "yes" -o "${sd_libpsml_enable}" = "auto"; then
    tmp_libpsml_vars="${LIBPSML_CPPFLAGS}${LIBPSML_FCFLAGS}${LIBPSML_LDFLAGS}${LIBPSML_LIBS}"
    if test "${sd_libpsml_init}" = "def" -a ! -z "${tmp_libpsml_vars}"; then
      sd_libpsml_enable="yes"
      sd_libpsml_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_LIBPSML_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_libpsml_enable}" = "yes" -o "${sd_libpsml_enable}" = "auto"; then

    # Set LibPSML-specific flags
    case "${sd_libpsml_init}" in

      def|yon)
        sd_libpsml_cppflags="${sd_libpsml_cppflags_def}"
        sd_libpsml_fcflags="${sd_libpsml_fcflags_def}"
        sd_libpsml_ldflags="${sd_libpsml_ldflags_def}"
        sd_libpsml_libs="${sd_libpsml_libs_def}"
        ;;

      dir)
        sd_libpsml_cppflags="-I${with_libpsml}/include"
        sd_libpsml_fcflags="${sd_libpsml_fcflags_def} -I${with_libpsml}/include"
        sd_libpsml_ldflags="${sd_libpsml_ldflags_def}"
        sd_libpsml_libs="-L${with_libpsml}/lib ${sd_libpsml_libs_def}"
        ;;

      env)
        sd_libpsml_cppflags="${sd_libpsml_cppflags_def}"
        sd_libpsml_fcflags="${sd_libpsml_fcflags_def}"
        sd_libpsml_ldflags="${sd_libpsml_ldflags_def}"
        sd_libpsml_libs="${sd_libpsml_libs_def}"
        test ! -z "${LIBPSML_CPPFLAGS}" && sd_libpsml_cppflags="${LIBPSML_CPPFLAGS}"
        test ! -z "${LIBPSML_FCFLAGS}" && sd_libpsml_fcflags="${LIBPSML_FCFLAGS}"
        test ! -z "${LIBPSML_LDFLAGS}" && sd_libpsml_ldflags="${LIBPSML_LDFLAGS}"
        test ! -z "${LIBPSML_LIBS}" && sd_libpsml_libs="${LIBPSML_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for LibPSML: '${sd_libpsml_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_libpsml_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_libpsml_init="esl"
    sd_libpsml_cppflags=""
    sd_libpsml_fcflags=""
    sd_libpsml_ldflags=""
    sd_libpsml_libs=""
  fi

  # Display configuration
  _SD_LIBPSML_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_libpsml_options)
  AC_SUBST(sd_libpsml_enable_def)
  AC_SUBST(sd_libpsml_policy)
  AC_SUBST(sd_libpsml_status)
  AC_SUBST(sd_libpsml_enable)
  AC_SUBST(sd_libpsml_init)
  AC_SUBST(sd_libpsml_ok)
  AC_SUBST(sd_libpsml_cppflags)
  AC_SUBST(sd_libpsml_fcflags)
  AC_SUBST(sd_libpsml_ldflags)
  AC_SUBST(sd_libpsml_libs)
  AC_SUBST(with_libpsml)

  # Clean-up
  unset tmp_libpsml_vars
]) # SD_LIBPSML_INIT


AC_DEFUN([SD_LIBPSML_DETECT], [
  # Display configuration
  _SD_LIBPSML_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_libpsml_enable}" = "auto" -o "${sd_libpsml_enable}" = "yes"; then
    _SD_LIBPSML_CHECK_USE

    if test "${sd_libpsml_ok}" = "yes"; then
      if test "${sd_libpsml_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_libpsml_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_libpsml_fcflags}"
        LIBS="${sd_libpsml_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_libpsml_ldflags}"

      AC_DEFINE([HAVE_LIBPSML], 1,
        [Define to 1 if you have the LibPSML library.])
    else
      if test "${sd_libpsml_status}" = "optional" -a \
              "${sd_libpsml_init}" = "def"; then
        sd_libpsml_enable="no"
        sd_libpsml_cppflags=""
        sd_libpsml_fcflags=""
        sd_libpsml_ldflags=""
        sd_libpsml_libs=""
      else
        AC_MSG_FAILURE([invalid LibPSML configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_LIBPSML_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_libpsml_init}" = "esl"; then
    AC_MSG_NOTICE([will look for LibPSML in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_libpsml_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_xmlf90_cppflags} ${sd_libpsml_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_xmlf90_fcflags} ${sd_libpsml_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_xmlf90_ldflags} ${sd_libpsml_ldflags}"
    LIBS="${sd_libpsml_libs} ${sd_xmlf90_libs} ${LIBS}"
  fi

  # Check LibPSML Fortran API
  AC_MSG_CHECKING([whether the LibPSML Fortran interface works])
  for tmp_incs in "" "-I/usr/include"; do
    FCFLAGS="${FCFLAGS} ${tmp_incs}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use m_psml
        use m_psml_api
        type(ps_t) :: psxml
        call ps_destroy(psxml)
      ]])], [sd_libpsml_ok="yes"], [sd_libpsml_ok="no"])
    AC_LANG_POP([Fortran])
    if test "${sd_libpsml_ok}" = "yes"; then
      test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
      break
    fi
  done
  AC_MSG_RESULT([${sd_libpsml_ok}])
  unset tmp_incs

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_LIBPSML_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_LIBPSML_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_libpsml_invalid="no"
  if test "${sd_libpsml_enable_def}" != "auto" -a \
          "${sd_libpsml_enable_def}" != "no" -a \
          "${sd_libpsml_enable_def}" != "yes"; then
    case "${sd_libpsml_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_libpsml_enable_def = '${sd_libpsml_enable_def}'])
        ;;
      skip)
        tmp_libpsml_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_libpsml_enable_def = '${sd_libpsml_enable_def}'])
        tmp_libpsml_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_libpsml_invalid}" = "yes"; then
    if test "${sd_libpsml_status}" = "required"; then
      sd_libpsml_enable_def="yes"
    else
      sd_libpsml_enable_def="no"
    fi
    tmp_libpsml_invalid="no"
    AC_MSG_NOTICE([setting sd_libpsml_enable_def to '${sd_libpsml_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_libpsml_invalid="no"
  if test "${sd_libpsml_status}" = "implicit" -o \
          "${sd_libpsml_status}" = "required"; then
    if test "${sd_libpsml_enable}" = "no"; then
      case "${sd_libpsml_policy}" in
        fail)
          AC_MSG_ERROR([The LibPSML package is required and cannot be disabled
                  See https://launchpad.net/libpsml for details on how to
                  install it.])
          ;;
        skip)
          tmp_libpsml_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The LibPSML package is required and cannot be disabled])
          tmp_libpsml_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_libpsml_enable}" = "auto"; then
      AC_MSG_NOTICE([setting LibPSML trigger to yes])
      sd_libpsml_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_libpsml_invalid}" = "yes"; then
    case "${sd_libpsml_status}" in
      implicit|required)
        sd_libpsml_enable="yes"
        ;;
      optional)
        sd_libpsml_enable="no"
        ;;
    esac
    tmp_libpsml_invalid="no"
    AC_MSG_NOTICE([setting sd_libpsml_enable to '${sd_libpsml_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_libpsml_vars="${LIBPSML_CPPFLAGS}${LIBPSML_FCFLAGS}${LIBPSML_LDFLAGS}${LIBPSML_LIBS}"
  tmp_libpsml_invalid="no"
  if test ! -z "${tmp_libpsml_vars}" -a ! -z "${with_libpsml}"; then
    case "${sd_libpsml_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for LibPSML
                  Please use LIBPSML_FCFLAGS + LIBPSML_LIBS or --with-libpsml,
                  not both.])
        ;;
      skip)
        tmp_libpsml_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for LibPSML])
        tmp_libpsml_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_libpsml_vars}"; then
    sd_libpsml_enable="yes"
    sd_libpsml_init="env"
    if test "${tmp_libpsml_invalid}" = "yes"; then
      tmp_libpsml_invalid="no"
      AC_MSG_NOTICE([overriding --with-libpsml with LIBPSML_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_libpsml_status}" = "implicit"; then
    if test "${sd_libpsml_fcflags}" != ""; then
      sd_libpsml_fcflags=""
      AC_MSG_NOTICE([resetting LibPSML Fortran flags (implicit package)])
    fi
    if test "${sd_libpsml_ldflags}" != ""; then
      sd_libpsml_ldflags=""
      AC_MSG_NOTICE([resetting LibPSML linker flags (implicit package)])
    fi
    if test "${sd_libpsml_libs}" != ""; then
      sd_libpsml_libs=""
      AC_MSG_NOTICE([resetting LibPSML library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_libpsml_enable}" = "implicit"; then
    sd_libpsml_fcflags=""
    sd_libpsml_ldflags=""
    sd_libpsml_libs=""
  fi

  # Clean-up
  unset tmp_libpsml_invalid
  unset tmp_libpsml_vars
]) # _SD_LIBPSML_CHECK_CONFIG


AC_DEFUN([_SD_LIBPSML_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable LibPSML])
  AC_MSG_RESULT([${sd_libpsml_enable}])
  if test "${sd_libpsml_enable}" != "no"; then
    AC_MSG_CHECKING([how LibPSML parameters have been set])
    AC_MSG_RESULT([${sd_libpsml_init}])
    AC_MSG_CHECKING([for LibPSML C preprocessing flags])
    if test "${sd_libpsml_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpsml_cppflags}])
    fi
    AC_MSG_CHECKING([for LibPSML Fortran flags])
    if test "${sd_libpsml_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpsml_fcflags}])
    fi
    AC_MSG_CHECKING([for LibPSML linker flags])
    if test "${sd_libpsml_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpsml_ldflags}])
    fi
    AC_MSG_CHECKING([for LibPSML library flags])
    if test "${sd_libpsml_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_libpsml_libs}])
    fi
  fi
]) # _SD_LIBPSML_DUMP_CONFIG

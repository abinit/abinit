## Copyright (C) 2019 Yann Pouillon

#
# XML Fortran I/O library (XMLF90)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_XMLF90_INIT], [
  # Init
  sd_xmlf90_cppflags=""
  sd_xmlf90_fcflags=""
  sd_xmlf90_ldflags=""
  sd_xmlf90_libs=""
  sd_xmlf90_enable=""
  sd_xmlf90_init="unknown"
  sd_xmlf90_ok="unknown"

  # Set adjustable parameters
  sd_xmlf90_options="$1"
  sd_xmlf90_libs_def="$2"
  sd_xmlf90_cppflags_def="$3"
  sd_xmlf90_cflags_def="$4"
  sd_xmlf90_cxxflags_def="$5"
  sd_xmlf90_fcflags_def="$6"
  sd_xmlf90_ldflags_def="$7"
  sd_xmlf90_enable_def=""
  sd_xmlf90_policy=""
  sd_xmlf90_status=""

  # Process options
  for kwd in ${sd_xmlf90_options}; do
    case "${kwd}" in
      auto)
        sd_xmlf90_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_xmlf90_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_xmlf90_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg XMLF90 option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_xmlf90_libs_def}" && sd_xmlf90_libs_def="-lxmlf90"
  test -z "${sd_xmlf90_policy}" && sd_xmlf90_policy="fail"
  test -z "${sd_xmlf90_status}" && sd_xmlf90_status="optional"
  test -z "${sd_xmlf90_enable_def}" && sd_xmlf90_enable_def="no"
  case "${sd_xmlf90_status}" in
    implicit|required)
      sd_xmlf90_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([xmlf90],
    [AS_HELP_STRING([--with-xmlf90],
      [Install prefix of the PSML I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_xmlf90_enable="${withval}"
        sd_xmlf90_init="yon"
      else
        sd_xmlf90_enable="yes"
        sd_xmlf90_init="dir"
      fi],
    [ sd_xmlf90_enable="${sd_xmlf90_enable_def}"; sd_xmlf90_init="def"])

  # Declare environment variables
  AC_ARG_VAR([XMLF90_CPPFLAGS], [C preprocessing flags for XMLF90.])
  AC_ARG_VAR([XMLF90_CFLAGS], [C flags for XMLF90.])
  AC_ARG_VAR([XMLF90_FCFLAGS], [Fortran flags for XMLF90.])
  AC_ARG_VAR([XMLF90_LDFLAGS], [Linker flags for XMLF90.])
  AC_ARG_VAR([XMLF90_LIBS], [Library flags for XMLF90.])

  # Detect use of environment variables
  if test "${sd_xmlf90_enable}" = "yes" -o "${sd_xmlf90_enable}" = "auto"; then
    tmp_xmlf90_vars="${XMLF90_CPPFLAGS}${XMLF90_CFLAGS}${XMLF90_FCFLAGS}${XMLF90_LDFLAGS}${XMLF90_LIBS}"
    if test "${sd_xmlf90_init}" = "def" -a ! -z "${tmp_xmlf90_vars}"; then
      sd_xmlf90_enable="yes"
      sd_xmlf90_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_XMLF90_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_xmlf90_enable}" = "yes" -o "${sd_xmlf90_enable}" = "auto"; then

    # Set XMLF90-specific flags
    case "${sd_xmlf90_init}" in

      def|yon)
        sd_xmlf90_cppflags="${sd_xmlf90_cppflags_def}"
        sd_xmlf90_fcflags="${sd_xmlf90_fcflags_def}"
        sd_xmlf90_ldflags="${sd_xmlf90_ldflags_def}"
        sd_xmlf90_libs="${sd_xmlf90_libs_def}"
        ;;

      dir)
        sd_xmlf90_cppflags="-I${with_xmlf90}/include"
        sd_xmlf90_fcflags="${sd_xmlf90_fcflags_def} -I${with_xmlf90}/include"
        sd_xmlf90_ldflags="${sd_xmlf90_ldflags_def}"
        sd_xmlf90_libs="-L${with_xmlf90}/lib ${sd_xmlf90_libs_def}"
        ;;

      env)
        sd_xmlf90_cppflags="${sd_xmlf90_cppflags_def}"
        sd_xmlf90_fcflags="${sd_xmlf90_fcflags_def}"
        sd_xmlf90_ldflags="${sd_xmlf90_ldflags_def}"
        sd_xmlf90_libs="${sd_xmlf90_libs_def}"
        test ! -z "${XMLF90_CPPFLAGS}" && sd_xmlf90_cppflags="${XMLF90_CPPFLAGS}"
        test ! -z "${XMLF90_FCFLAGS}" && sd_xmlf90_fcflags="${XMLF90_FCFLAGS}"
        test ! -z "${XMLF90_LDFLAGS}" && sd_xmlf90_ldflags="${XMLF90_LDFLAGS}"
        test ! -z "${XMLF90_LIBS}" && sd_xmlf90_libs="${XMLF90_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for XMLF90: '${sd_xmlf90_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_xmlf90_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_xmlf90_init="esl"
    sd_xmlf90_cppflags=""
    sd_xmlf90_fcflags=""
    sd_xmlf90_ldflags=""
    sd_xmlf90_libs=""
  fi

  # Display configuration
  _SD_XMLF90_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_xmlf90_options)
  AC_SUBST(sd_xmlf90_enable_def)
  AC_SUBST(sd_xmlf90_policy)
  AC_SUBST(sd_xmlf90_status)
  AC_SUBST(sd_xmlf90_enable)
  AC_SUBST(sd_xmlf90_init)
  AC_SUBST(sd_xmlf90_ok)
  AC_SUBST(sd_xmlf90_cppflags)
  AC_SUBST(sd_xmlf90_fcflags)
  AC_SUBST(sd_xmlf90_ldflags)
  AC_SUBST(sd_xmlf90_libs)
  AC_SUBST(with_xmlf90)

  # Clean-up
  unset tmp_xmlf90_vars
]) # SD_XMLF90_INIT


AC_DEFUN([SD_XMLF90_DETECT], [
  # Display configuration
  _SD_XMLF90_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_xmlf90_enable}" = "auto" -o "${sd_xmlf90_enable}" = "yes"; then
    _SD_XMLF90_CHECK_USE

    if test "${sd_xmlf90_ok}" = "yes"; then
      if test "${sd_xmlf90_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_xmlf90_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_xmlf90_fcflags}"
        LIBS="${sd_xmlf90_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_xmlf90_ldflags}"

      AC_DEFINE([HAVE_XMLF90], 1,
        [Define to 1 if you have the XMLF90 library.])
    else
      if test "${sd_xmlf90_status}" = "optional" -a \
              "${sd_xmlf90_init}" = "def"; then
        sd_xmlf90_enable="no"
        sd_xmlf90_cppflags=""
        sd_xmlf90_fcflags=""
        sd_xmlf90_ldflags=""
        sd_xmlf90_libs=""
      else
        AC_MSG_FAILURE([invalid XMLF90 configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_XMLF90_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_xmlf90_init}" = "esl"; then
    AC_MSG_NOTICE([will look for XMLF90 in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_xmlf90_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_xmlf90_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_xmlf90_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_xmlf90_ldflags}"
    LIBS="${sd_xmlf90_libs} ${LIBS}"
  fi

  # Check XMLF90 Fortran API
  AC_MSG_CHECKING([whether the XMLF90 Fortran interface works])
  for tmp_incs in "" "-I/usr/include"; do
    FCFLAGS="${FCFLAGS} ${tmp_incs}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use xmlf90_sax
        character(len=80) :: fname
        type(xml_t) :: fxml
        integer :: iostat
        call open_xmlfile(fname,fxml,iostat)
        call close_xmlfile(fxml)
      ]])], [sd_xmlf90_ok="yes"], [sd_xmlf90_ok="no"])
    AC_LANG_POP([Fortran])
    if test "${sd_xmlf90_ok}" = "yes"; then
      test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
      break
    fi
  done
  AC_MSG_RESULT([${sd_xmlf90_ok}])
  unset tmp_incs

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_XMLF90_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_XMLF90_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_xmlf90_invalid="no"
  if test "${sd_xmlf90_enable_def}" != "auto" -a \
          "${sd_xmlf90_enable_def}" != "no" -a \
          "${sd_xmlf90_enable_def}" != "yes"; then
    case "${sd_xmlf90_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_xmlf90_enable_def = '${sd_xmlf90_enable_def}'])
        ;;
      skip)
        tmp_xmlf90_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_xmlf90_enable_def = '${sd_xmlf90_enable_def}'])
        tmp_xmlf90_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_xmlf90_invalid}" = "yes"; then
    if test "${sd_xmlf90_status}" = "required"; then
      sd_xmlf90_enable_def="yes"
    else
      sd_xmlf90_enable_def="no"
    fi
    tmp_xmlf90_invalid="no"
    AC_MSG_NOTICE([setting sd_xmlf90_enable_def to '${sd_xmlf90_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_xmlf90_invalid="no"
  if test "${sd_xmlf90_status}" = "implicit" -o \
          "${sd_xmlf90_status}" = "required"; then
    if test "${sd_xmlf90_enable}" = "no"; then
      case "${sd_xmlf90_policy}" in
        fail)
          AC_MSG_ERROR([The XMLF90 package is required and cannot be disabled
                  See https://launchpad.net/xmlf90 for details on how to
                  install it.])
          ;;
        skip)
          tmp_xmlf90_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The XMLF90 package is required and cannot be disabled])
          tmp_xmlf90_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_xmlf90_enable}" = "auto"; then
      AC_MSG_NOTICE([setting XMLF90 trigger to yes])
      sd_xmlf90_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_xmlf90_invalid}" = "yes"; then
    case "${sd_xmlf90_status}" in
      implicit|required)
        sd_xmlf90_enable="yes"
        ;;
      optional)
        sd_xmlf90_enable="no"
        ;;
    esac
    tmp_xmlf90_invalid="no"
    AC_MSG_NOTICE([setting sd_xmlf90_enable to '${sd_xmlf90_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_xmlf90_vars="${XMLF90_FCFLAGS}${XMLF90_LDFLAGS}${XMLF90_LIBS}"
  tmp_xmlf90_invalid="no"
  if test ! -z "${tmp_xmlf90_vars}" -a ! -z "${with_xmlf90}"; then
    case "${sd_xmlf90_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for XMLF90
                  Please use XMLF90_FCFLAGS + XMLF90_LIBS or --with-xmlf90,
                  not both.])
        ;;
      skip)
        tmp_xmlf90_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for XMLF90])
        tmp_xmlf90_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_xmlf90_vars}"; then
    sd_xmlf90_enable="yes"
    sd_xmlf90_init="env"
    if test "${tmp_xmlf90_invalid}" = "yes"; then
      tmp_xmlf90_invalid="no"
      AC_MSG_NOTICE([overriding --with-xmlf90 with XMLF90_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_xmlf90_status}" = "implicit"; then
    if test "${sd_xmlf90_fcflags}" != ""; then
      sd_xmlf90_fcflags=""
      AC_MSG_NOTICE([resetting XMLF90 Fortran flags (implicit package)])
    fi
    if test "${sd_xmlf90_ldflags}" != ""; then
      sd_xmlf90_ldflags=""
      AC_MSG_NOTICE([resetting XMLF90 linker flags (implicit package)])
    fi
    if test "${sd_xmlf90_libs}" != ""; then
      sd_xmlf90_libs=""
      AC_MSG_NOTICE([resetting XMLF90 library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_xmlf90_enable}" = "implicit"; then
    sd_xmlf90_fcflags=""
    sd_xmlf90_ldflags=""
    sd_xmlf90_libs=""
  fi

  # Clean-up
  unset tmp_xmlf90_invalid
  unset tmp_xmlf90_vars
]) # _SD_XMLF90_CHECK_CONFIG


AC_DEFUN([_SD_XMLF90_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable XMLF90])
  AC_MSG_RESULT([${sd_xmlf90_enable}])
  if test "${sd_xmlf90_enable}" != "no"; then
    AC_MSG_CHECKING([how XMLF90 parameters have been set])
    AC_MSG_RESULT([${sd_xmlf90_init}])
    AC_MSG_CHECKING([for XMLF90 C preprocessing flags])
    if test "${sd_xmlf90_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_xmlf90_cppflags}])
    fi
    AC_MSG_CHECKING([for XMLF90 Fortran flags])
    if test "${sd_xmlf90_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_xmlf90_fcflags}])
    fi
    AC_MSG_CHECKING([for XMLF90 linker flags])
    if test "${sd_xmlf90_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_xmlf90_ldflags}])
    fi
    AC_MSG_CHECKING([for XMLF90 library flags])
    if test "${sd_xmlf90_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_xmlf90_libs}])
    fi
  fi
]) # _SD_XMLF90_DUMP_CONFIG

## Copyright (C) 2019-2022 ABINIT group (Yann Pouillon)

#
# YAKL library (Yet Another Kernel Launcher)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_YAKL_INIT], [
  # Init
  sd_yakl_cppflags=""
  sd_yakl_fcflags=""
  sd_yakl_ldflags=""
  sd_yakl_libs=""
  sd_yakl_enable=""
  sd_yakl_init="unknown"
  sd_yakl_fortran_ok="unknown"
  sd_yakl_ok="unknown"

  # Set adjustable parameters
  sd_yakl_options="$1"
  sd_yakl_libs_def="$2"
  sd_yakl_cppflags_def="$3"
  sd_yakl_fcflags_def="$4"
  sd_yakl_ldflags_def="$5"
  sd_yakl_enable_def=""
  sd_yakl_policy=""
  sd_yakl_status=""

  # Process options
  for kwd in ${sd_yakl_options}; do
    case "${kwd}" in
      auto)
        sd_yakl_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_yakl_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_yakl_policy="${kwd}"
        ;;
      mandatory)
        sd_yakl_enable="yes"
        sd_yakl_enable_def="yes"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg YAKL option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_yakl_libs_def}" && sd_yakl_libs_def="-lyakl_fortran_interface"
  test -z "${sd_yakl_policy}" && sd_yakl_policy="fail"
  test -z "${sd_yakl_status}" && sd_yakl_status="optional"
  test -z "${sd_yakl_enable_def}" && sd_yakl_enable_def="no"
  case "${sd_yakl_status}" in
    implicit|required)
      sd_yakl_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([yakl],
    [AS_HELP_STRING([--with-yakl],
      [Install prefix of the YAKL library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_yakl_enable="${withval}"
        sd_yakl_init="yon"
      else
        sd_yakl_enable="yes"
        sd_yakl_init="dir"
      fi],
    [ sd_yakl_enable="${sd_yakl_enable_def}"; sd_yakl_init="def"])

  # Declare environment variables
  AC_ARG_VAR([YAKL_CPPFLAGS], [C preprocessing flags for YAKL.])
  AC_ARG_VAR([YAKL_FCFLAGS], [Fortran flags for YAKL.])
  AC_ARG_VAR([YAKL_LDFLAGS], [Linker flags for YAKL.])
  AC_ARG_VAR([YAKL_LIBS], [Library flags for YAKL.])

  # Detect use of environment variables
  if test "${sd_yakl_enable}" = "yes" -o "${sd_yakl_enable}" = "auto"; then
    tmp_yakl_vars="${YAKL_CPPFLAGS}${YAKL_FCFLAGS}${YAKL_LDFLAGS}${YAKL_LIBS}"
    if test "${sd_yakl_init}" = "def" -a ! -z "${tmp_yakl_vars}"; then
      sd_yakl_enable="yes"
      sd_yakl_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_YAKL_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_yakl_enable}" = "yes" -o "${sd_yakl_enable}" = "auto"; then

    # Set YAKL-specific flags
    case "${sd_yakl_init}" in

      def|yon)
        sd_yakl_cppflags="${sd_yakl_cppflags_def}"
        sd_yakl_fcflags="${sd_yakl_fcflags_def}"
        sd_yakl_ldflags="${sd_yakl_ldflags_def}"
        sd_yakl_libs="${sd_yakl_libs_def}"
        ;;

      dir)
        sd_yakl_cppflags="-I${with_yakl}/include"
        sd_yakl_fcflags="${sd_yakl_fcflags_def} -I${with_yakl}/include"
        sd_yakl_ldflags="${sd_yakl_ldflags_def}"
        sd_yakl_libs="-L${with_yakl}/lib64 ${sd_yakl_libs_def} -lstdc++ -ldl -L${with_gpu}/lib64 -lcudart"
        ;;

      env)
        sd_yakl_cppflags="${sd_yakl_cppflags_def}"
        sd_yakl_fcflags="${sd_yakl_fcflags_def}"
        sd_yakl_ldflags="${sd_yakl_ldflags_def}"
        sd_yakl_libs="${sd_yakl_libs_def}"
        test ! -z "${YAKL_CPPFLAGS}" && sd_yakl_cppflags="${YAKL_CPPFLAGS}"
        test ! -z "${YAKL_LDFLAGS}" && sd_yakl_ldflags="${YAKL_LDFLAGS}"
        test ! -z "${YAKL_LIBS}" && sd_yakl_libs="${YAKL_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for YAKL: '${sd_yakl_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_yakl_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_yakl_init="esl"
    sd_yakl_cppflags=""
    sd_yakl_fcflags=""
    sd_yakl_ldflags=""
    sd_yakl_libs=""
  fi

  # Display configuration
  _SD_YAKL_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_yakl_options)
  AC_SUBST(sd_yakl_enable_def)
  AC_SUBST(sd_yakl_policy)
  AC_SUBST(sd_yakl_status)
  AC_SUBST(sd_yakl_enable)
  AC_SUBST(sd_yakl_init)
  AC_SUBST(sd_yakl_ok)
  AC_SUBST(sd_yakl_cppflags)
  AC_SUBST(sd_yakl_fcflags)
  AC_SUBST(sd_yakl_ldflags)
  AC_SUBST(sd_yakl_libs)
  AC_SUBST(with_yakl)

  # Clean-up
  unset tmp_yakl_vars
]) # SD_YAKL_INIT


AC_DEFUN([SD_YAKL_DETECT], [
  # Display configuration
  _SD_YAKL_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_yakl_enable}" = "auto" -o "${sd_yakl_enable}" = "yes"; then
    _SD_YAKL_CHECK_USE

    if test "${sd_yakl_ok}" = "yes"; then
      if test "${sd_yakl_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_yakl_libs_def} ${sd_esl_bundle_libs}"
      else
        CPPFLAGS="${CPPFLAGS} ${sd_yakl_cppflags}"
        FCFLAGS="${FCFLAGS} ${sd_yakl_fcflags}"
        FCFLAGS="${FCFLAGS} ${sd_yakl_fcflags}"
        LIBS="${sd_yakl_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_yakl_ldflags}"

      AC_DEFINE([HAVE_YAKL], 1,
        [Define to 1 if you have the YAKL library.])
    else
      if test "${sd_yakl_policy}" = "fail"; then
        AC_MSG_FAILURE([invalid YAKL configuration])
      fi
      if test "${sd_yakl_status}" = "optional" -a \
              "${sd_yakl_init}" = "def"; then
        sd_yakl_enable="no"
        sd_yakl_cppflags=""
        sd_yakl_fcflags=""
        sd_yakl_ldflags=""
        sd_yakl_libs=""
      else
        AC_MSG_WARN([invalid YAKL configuration])
      fi
    fi
  else
    sd_yakl_enable="no"
    sd_yakl_cppflags=""
    sd_yakl_fcflags=""
    sd_yakl_ldflags=""
    sd_yakl_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_YAKL_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_yakl_init}" = "esl"; then
    AC_MSG_NOTICE([will look for YAKL in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_yakl_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_yakl_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_yakl_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_yakl_ldflags}"
    LIBS="${sd_yakl_libs} ${LIBS}"
  fi

  # Check YAKL CXX API
  AC_MSG_CHECKING([whether the YAKL Fortran module works])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
       use gator_mod
       call gator_init()
    ]])], [sd_yakl_fortran_ok="yes"], [sd_yakl_fortran_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_yakl_fortran_ok}])

  # Combine the available results
  sd_yakl_ok="no"
  if test "${sd_yakl_fortran_ok}" = "yes"; then
    sd_yakl_ok="yes"
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_YAKL_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_YAKL_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_yakl_invalid="no"
  if test "${sd_yakl_enable_def}" != "auto" -a \
          "${sd_yakl_enable_def}" != "no" -a \
          "${sd_yakl_enable_def}" != "yes"; then
    case "${sd_yakl_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_yakl_enable_def = '${sd_yakl_enable_def}'])
        ;;
      skip)
        tmp_yakl_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_yakl_enable_def = '${sd_yakl_enable_def}'])
        tmp_yakl_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_yakl_invalid}" = "yes"; then
    if test "${sd_yakl_status}" = "required"; then
      sd_yakl_enable_def="yes"
    else
      sd_yakl_enable_def="no"
    fi
    tmp_yakl_invalid="no"
    AC_MSG_NOTICE([setting sd_yakl_enable_def to '${sd_yakl_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_yakl_invalid="no"
  if test "${sd_yakl_status}" = "implicit" -o \
          "${sd_yakl_status}" = "required"; then
    if test "${sd_yakl_enable}" = "no"; then
      case "${sd_yakl_policy}" in
        fail)
          AC_MSG_ERROR([The YAKL package is required and cannot be disabled
                  See https://tddft.org/programs/yakl/ for details on how to
                  install it.])
          ;;
        skip)
          tmp_yakl_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The YAKL package is required and cannot be disabled])
          tmp_yakl_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_yakl_enable}" = "auto"; then
      AC_MSG_NOTICE([setting YAKL trigger to yes])
      sd_yakl_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_yakl_invalid}" = "yes"; then
    case "${sd_yakl_status}" in
      implicit|required)
        sd_yakl_enable="yes"
        ;;
      optional)
        sd_yakl_enable="no"
        ;;
    esac
    tmp_yakl_invalid="no"
    AC_MSG_NOTICE([setting sd_yakl_enable to '${sd_yakl_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_yakl_vars="${YAKL_CPPFLAGS}${YAKL_CFLAGS}${YAKL_FCFLAGS}${YAKL_LDFLAGS}${YAKL_LIBS}"
  tmp_yakl_invalid="no"
  if test ! -z "${tmp_yakl_vars}" -a ! -z "${with_yakl}"; then
    case "${sd_yakl_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for YAKL
                  Please use YAKL_FCFLAGS + YAKL_LIBS or --with-yakl,
                  not both.])
        ;;
      skip)
        tmp_yakl_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for YAKL])
        tmp_yakl_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_yakl_vars}"; then
    sd_yakl_enable="yes"
    sd_yakl_init="env"
    if test "${tmp_yakl_invalid}" = "yes"; then
      tmp_yakl_invalid="no"
      AC_MSG_NOTICE([overriding --with-yakl with YAKL_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_yakl_status}" = "implicit"; then
    if test "${sd_yakl_fcflags}" != ""; then
      sd_yakl_fcflags=""
      AC_MSG_NOTICE([resetting YAKL Fortran flags (implicit package)])
    fi
    if test "${sd_yakl_ldflags}" != ""; then
      sd_yakl_ldflags=""
      AC_MSG_NOTICE([resetting YAKL linker flags (implicit package)])
    fi
    if test "${sd_yakl_libs}" != ""; then
      sd_yakl_libs=""
      AC_MSG_NOTICE([resetting YAKL library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_yakl_enable}" = "implicit"; then
    sd_yakl_fcflags=""
    sd_yakl_ldflags=""
    sd_yakl_libs=""
  fi

  # Clean-up
  unset tmp_yakl_invalid
  unset tmp_yakl_vars
]) # _SD_YAKL_CHECK_CONFIG


AC_DEFUN([_SD_YAKL_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable YAKL])
  AC_MSG_RESULT([${sd_yakl_enable}])
  if test "${sd_yakl_enable}" != "no"; then
    AC_MSG_CHECKING([how YAKL parameters have been set])
    AC_MSG_RESULT([${sd_yakl_init}])
    AC_MSG_CHECKING([for YAKL C preprocessing flags])
    if test "${sd_yakl_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_yakl_cppflags}])
    fi
    AC_MSG_CHECKING([for YAKL Fortran flags])
    if test "${sd_yakl_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_yakl_fcflags}])
    fi
    AC_MSG_CHECKING([for YAKL linker flags])
    if test "${sd_yakl_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_yakl_ldflags}])
    fi
    AC_MSG_CHECKING([for YAKL library flags])
    if test "${sd_yakl_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_yakl_libs}])
    fi
  fi
]) # _SD_YAKL_DUMP_CONFIG

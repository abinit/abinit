## Copyright (C) 2019 Yann Pouillon

#
# Fastest Fourier Transform in the West library (DFTI)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_DFTI_INIT], [
  # Init
  sd_dfti_cppflags=""
  sd_dfti_cflags=""
  sd_dfti_cxxflags=""
  sd_dfti_fcflags=""
  sd_dfti_ldflags=""
  sd_dfti_libs=""
  sd_dfti_enable=""
  sd_dfti_init="unknown"
  sd_dfti_mpi_ok="unknown"
  sd_dfti_ok="unknown"

  # Set adjustable parameters
  sd_dfti_options="$1"
  sd_dfti_libs_def="$2"
  sd_dfti_cppflags_def="$3"
  sd_dfti_cflags_def="$4"
  sd_dfti_cxxflags_def="$5"
  sd_dfti_fcflags_def="$6"
  sd_dfti_ldflags_def="$7"
  sd_dfti_enable_def=""
  sd_dfti_policy=""
  sd_dfti_status=""

  # Process options
  for kwd in ${sd_dfti_options}; do
    case "${kwd}" in
      auto)
        sd_dfti_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_dfti_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_dfti_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg DFTI option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  if test "${sd_mpi_enable}" = "yes"; then
    test -z "${sd_dfti_libs_def}" && sd_dfti_libs_def="-ldfti_mpi -ldfti"
  else
    test -z "${sd_dfti_libs_def}" && sd_dfti_libs_def="-ldfti"
  fi
  test -z "${sd_dfti_policy}" && sd_dfti_policy="fail"
  test -z "${sd_dfti_status}" && sd_dfti_status="optional"
  test -z "${sd_dfti_enable_def}" && sd_dfti_enable_def="no"
  case "${sd_dfti_status}" in
    implicit|required)
      sd_dfti_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([dfti],
    [AS_HELP_STRING([--with-dfti],
      [Install prefix of the DFTI library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_dfti_enable="${withval}"
        sd_dfti_init="yon"
      else
        sd_dfti_enable="yes"
        sd_dfti_init="dir"
      fi],
    [ sd_dfti_enable="${sd_dfti_enable_def}"; sd_dfti_init="def"])

  # Declare environment variables
  AC_ARG_VAR([DFTI_CPPFLAGS], [C preprocessing flags for DFTI.])
  AC_ARG_VAR([DFTI_CFLAGS], [C flags for DFTI.])
  AC_ARG_VAR([DFTI_FCFLAGS], [Fortran flags for DFTI.])
  AC_ARG_VAR([DFTI_LDFLAGS], [Linker flags for DFTI.])
  AC_ARG_VAR([DFTI_LIBS], [Library flags for DFTI.])

  # Detect use of environment variables
  if test "${sd_dfti_enable}" = "yes" -o "${sd_dfti_enable}" = "auto"; then
    tmp_dfti_vars="${DFTI_CPPFLAGS}${DFTI_CFLAGS}${DFTI_FCFLAGS}${DFTI_LDFLAGS}${DFTI_LIBS}"
    if test "${sd_dfti_init}" = "def" -a ! -z "${tmp_dfti_vars}"; then
      sd_dfti_enable="yes"
      sd_dfti_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_DFTI_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_dfti_enable}" = "yes" -o "${sd_dfti_enable}" = "auto"; then

    # Set DFTI-specific flags
    case "${sd_dfti_init}" in

      def|yon)
        sd_dfti_cppflags="${sd_dfti_cppflags_def}"
        sd_dfti_cflags="${sd_dfti_cflags_def}"
        sd_dfti_fcflags="${sd_dfti_fcflags_def}"
        sd_dfti_ldflags="${sd_dfti_ldflags_def}"
        sd_dfti_libs="${sd_dfti_libs_def}"
        ;;

      dir)
        sd_dfti_cppflags="-I${with_dfti}/include"
        sd_dfti_cflags="${sd_dfti_cflags_def}"
        sd_dfti_fcflags="${sd_dfti_fcflags_def}"
        sd_dfti_ldflags="${sd_dfti_ldflags_def}"
        sd_dfti_libs="-L${with_dfti}/lib ${sd_dfti_libs_def}"
        ;;

      env)
        sd_dfti_cppflags="${sd_dfti_cppflags_def}"
        sd_dfti_cflags="${sd_dfti_cflags_def}"
        sd_dfti_fcflags="${sd_dfti_fcflags_def}"
        sd_dfti_ldflags="${sd_dfti_ldflags_def}"
        sd_dfti_libs="${sd_dfti_libs_def}"
        test ! -z "${DFTI_CPPFLAGS}" && sd_dfti_cppflags="${DFTI_CPPFLAGS}"
        test ! -z "${DFTI_CFLAGS}" && sd_dfti_cflags="${DFTI_CFLAGS}"
        test ! -z "${DFTI_FCFLAGS}" && sd_dfti_fcflags="${DFTI_FCFLAGS}"
        test ! -z "${DFTI_LDFLAGS}" && sd_dfti_ldflags="${DFTI_LDFLAGS}"
        test ! -z "${DFTI_LIBS}" && sd_dfti_libs="${DFTI_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for DFTI: '${sd_dfti_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_dfti_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_dfti_init="esl"
    sd_dfti_cppflags=""
    sd_dfti_cflags=""
    sd_dfti_fcflags=""
    sd_dfti_ldflags=""
    sd_dfti_libs=""
  fi

  # Display configuration
  _SD_DFTI_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_dfti_options)
  AC_SUBST(sd_dfti_enable_def)
  AC_SUBST(sd_dfti_policy)
  AC_SUBST(sd_dfti_status)
  AC_SUBST(sd_dfti_enable)
  AC_SUBST(sd_dfti_init)
  AC_SUBST(sd_dfti_ok)
  AC_SUBST(sd_dfti_cppflags)
  AC_SUBST(sd_dfti_cflags)
  AC_SUBST(sd_dfti_fcflags)
  AC_SUBST(sd_dfti_ldflags)
  AC_SUBST(sd_dfti_libs)
  AC_SUBST(with_dfti)

  # Clean-up
  unset tmp_dfti_vars
]) # SD_DFTI_INIT


AC_DEFUN([SD_DFTI_DETECT], [
  # Display configuration
  _SD_DFTI_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_dfti_enable}" = "auto" -o "${sd_dfti_enable}" = "yes"; then
    _SD_DFTI_CHECK_USE

    if test "${sd_dfti_ok}" = "yes"; then
      if test "${sd_dfti_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_dfti_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_dfti_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_dfti_ldflags}"

      AC_DEFINE([HAVE_DFTI], 1,
        [Define to 1 if you have the DFTI library.])
    else
      if test "${sd_dfti_status}" = "optional" -a \
              "${sd_dfti_init}" = "def"; then
        sd_dfti_enable="no"
        sd_dfti_cppflags=""
        sd_dfti_cflags=""
        sd_dfti_fcflags=""
        sd_dfti_ldflags=""
        sd_dfti_libs=""
      else
        AC_MSG_FAILURE([invalid DFTI configuration])
      fi
    fi
  else
    sd_dfti_enable="no"
    sd_dfti_cppflags=""
    sd_dfti_cflags=""
    sd_dfti_fcflags=""
    sd_dfti_ldflags=""
    sd_dfti_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_DFTI_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_dfti_init}" = "esl"; then
    AC_MSG_NOTICE([will look for DFTI in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_dfti_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_dfti_cppflags}"
    CFLAGS="${CFLAGS} ${sd_dfti_cflags}"
    FCFLAGS="${FCFLAGS} ${sd_dfti_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_dfti_ldflags}"
    LIBS="${sd_dfti_libs} ${LIBS}"
  fi

  # Check DFTI C API
  # FIXME: Very complex to have it work properly, would need a replacement
  #        of AC_LINK_IFELSE accepting prologues, because of the included
  #        mkl_dfti.f90 file.
  AC_MSG_CHECKING([whether the DFTI library works])
  sd_dfti_ok="yes"
  AC_MSG_RESULT([${sd_dfti_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_DFTI_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_DFTI_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_dfti_invalid="no"
  if test "${sd_dfti_enable_def}" != "auto" -a \
          "${sd_dfti_enable_def}" != "no" -a \
          "${sd_dfti_enable_def}" != "yes"; then
    case "${sd_dfti_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_dfti_enable_def = '${sd_dfti_enable_def}'])
        ;;
      skip)
        tmp_dfti_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_dfti_enable_def = '${sd_dfti_enable_def}'])
        tmp_dfti_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_dfti_invalid}" = "yes"; then
    if test "${sd_dfti_status}" = "required"; then
      sd_dfti_enable_def="yes"
    else
      sd_dfti_enable_def="no"
    fi
    tmp_dfti_invalid="no"
    AC_MSG_NOTICE([setting sd_dfti_enable_def to '${sd_dfti_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_dfti_invalid="no"
  if test "${sd_dfti_status}" = "implicit" -o \
          "${sd_dfti_status}" = "required"; then
    if test "${sd_dfti_enable}" = "no"; then
      case "${sd_dfti_policy}" in
        fail)
          AC_MSG_ERROR([The DFTI package is required and cannot be disabled
                  See http://www.fftw.org/ for details on how to install it.])
          ;;
        skip)
          tmp_dfti_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The DFTI package is required and cannot be disabled])
          tmp_dfti_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_dfti_enable}" = "auto"; then
      AC_MSG_NOTICE([setting DFTI trigger to yes])
      sd_dfti_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_dfti_invalid}" = "yes"; then
    case "${sd_dfti_status}" in
      implicit|required)
        sd_dfti_enable="yes"
        ;;
      optional)
        sd_dfti_enable="no"
        ;;
    esac
    tmp_dfti_invalid="no"
    AC_MSG_NOTICE([setting sd_dfti_enable to '${sd_dfti_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_dfti_vars="${DFTI_CPPFLAGS}${DFTI_CFLAGS}${DFTI_LDFLAGS}${DFTI_LIBS}"
  tmp_dfti_invalid="no"
  if test ! -z "${tmp_dfti_vars}" -a ! -z "${with_dfti}"; then
    case "${sd_dfti_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for DFTI
                  Please use DFTI_CFLAGS + DFTI_LIBS or --with-dfti,
                  not both.])
        ;;
      skip)
        tmp_dfti_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for DFTI])
        tmp_dfti_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_dfti_vars}"; then
    sd_dfti_enable="yes"
    sd_dfti_init="env"
    if test "${tmp_dfti_invalid}" = "yes"; then
      tmp_dfti_invalid="no"
      AC_MSG_NOTICE([overriding --with-dfti with DFTI_{CPPFLAGS,CFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_dfti_status}" = "implicit"; then
    if test "${sd_dfti_ldflags}" != ""; then
      sd_dfti_ldflags=""
      AC_MSG_NOTICE([resetting DFTI linker flags (implicit package)])
    fi
    if test "${sd_dfti_libs}" != ""; then
      sd_dfti_libs=""
      AC_MSG_NOTICE([resetting DFTI library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_dfti_enable}" = "implicit"; then
    sd_dfti_cppflags=""
    sd_dfti_cflags=""
    sd_dfti_fcflags=""
    sd_dfti_ldflags=""
    sd_dfti_libs=""
  fi

  # Clean-up
  unset tmp_dfti_invalid
  unset tmp_dfti_vars
]) # _SD_DFTI_CHECK_CONFIG


AC_DEFUN([_SD_DFTI_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable DFTI])
  AC_MSG_RESULT([${sd_dfti_enable}])
  if test "${sd_dfti_enable}" != "no"; then
    AC_MSG_CHECKING([how DFTI parameters have been set])
    AC_MSG_RESULT([${sd_dfti_init}])
    AC_MSG_CHECKING([for DFTI C preprocessing flags])
    if test "${sd_dfti_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_cppflags}])
    fi
    AC_MSG_CHECKING([for DFTI C flags])
    if test "${sd_dfti_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_cflags}])
    fi
    AC_MSG_CHECKING([for DFTI Fortran flags])
    if test "${sd_dfti_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_fcflags}])
    fi
    AC_MSG_CHECKING([for DFTI linker flags])
    if test "${sd_dfti_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_ldflags}])
    fi
    AC_MSG_CHECKING([for DFTI library flags])
    if test "${sd_dfti_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_libs}])
    fi
  fi
]) # _SD_DFTI_DUMP_CONFIG

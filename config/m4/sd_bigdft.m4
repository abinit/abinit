## Copyright (C) 2019 Yann Pouillon

#
# Maximally-Localized Wannier Functions library (BigDFT)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_BIGDFT_INIT], [
  # Init
  sd_bigdft_cppflags=""
  sd_bigdft_cflags=""
  sd_bigdft_cxxflags=""
  sd_bigdft_fcflags=""
  sd_bigdft_ldflags=""
  sd_bigdft_libs=""
  sd_bigdft_enable=""
  sd_bigdft_init="unknown"
  sd_bigdft_ok="unknown"

  # Set adjustable parameters
  sd_bigdft_options="$1"
  sd_bigdft_libs_def="$2"
  sd_bigdft_cppflags_def="$3"
  sd_bigdft_cflags_def="$4"
  sd_bigdft_cxxflags_def="$5"
  sd_bigdft_fcflags_def="$6"
  sd_bigdft_ldflags_def="$7"
  sd_bigdft_enable_def=""
  sd_bigdft_policy=""
  sd_bigdft_status=""

  # Process options
  for kwd in ${sd_bigdft_options}; do
    case "${kwd}" in
      auto)
        sd_bigdft_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_bigdft_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_bigdft_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg BigDFT option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_bigdft_libs_def}" && sd_bigdft_libs_def="-lbigdft-1 -labinit -lpaw_bigdft -lyaml -lrt"
  test -z "${sd_bigdft_policy}" && sd_bigdft_policy="fail"
  test -z "${sd_bigdft_status}" && sd_bigdft_status="optional"
  test -z "${sd_bigdft_enable_def}" && sd_bigdft_enable_def="no"
  case "${sd_bigdft_status}" in
    implicit|required)
      sd_bigdft_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([bigdft],
    [AS_HELP_STRING([--with-bigdft],
      [Install prefix of the PSML I/O library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_bigdft_enable="${withval}"
        sd_bigdft_init="yon"
      else
        sd_bigdft_enable="yes"
        sd_bigdft_init="dir"
      fi],
    [ sd_bigdft_enable="${sd_bigdft_enable_def}"; sd_bigdft_init="def"])

  # Declare environment variables
  AC_ARG_VAR([BIGDFT_CPPFLAGS], [C preprocessing flags for BigDFT.])
  AC_ARG_VAR([BIGDFT_FCFLAGS], [Fortran flags for BigDFT.])
  AC_ARG_VAR([BIGDFT_LDFLAGS], [Linker flags for BigDFT.])
  AC_ARG_VAR([BIGDFT_LIBS], [Library flags for BigDFT.])

  # Detect use of environment variables
  if test "${sd_bigdft_enable}" = "yes" -o "${sd_bigdft_enable}" = "auto"; then
    tmp_bigdft_vars="${BIGDFT_CPPFLAGS}${BIGDFT_FCFLAGS}${BIGDFT_LDFLAGS}${BIGDFT_LIBS}"
    if test "${sd_bigdft_init}" = "def" -a ! -z "${tmp_bigdft_vars}"; then
      sd_bigdft_enable="yes"
      sd_bigdft_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_BIGDFT_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_bigdft_enable}" = "yes" -o "${sd_bigdft_enable}" = "auto"; then

    # Set BigDFT-specific flags
    case "${sd_bigdft_init}" in

      def|yon)
        sd_bigdft_cppflags="${sd_bigdft_cppflags_def}"
        sd_bigdft_fcflags="${sd_bigdft_fcflags_def}"
        sd_bigdft_ldflags="${sd_bigdft_ldflags_def}"
        sd_bigdft_libs="${sd_bigdft_libs_def}"
        ;;

      dir)
        sd_bigdft_cppflags="-I${with_bigdft}/include"
        sd_bigdft_fcflags="${sd_bigdft_fcflags_def} -I${with_bigdft}/include"
        sd_bigdft_ldflags="${sd_bigdft_ldflags_def}"
        sd_bigdft_libs="-L${with_bigdft}/lib ${sd_bigdft_libs_def}"
        ;;

      env)
        sd_bigdft_cppflags="${sd_bigdft_cppflags_def}"
        sd_bigdft_fcflags="${sd_bigdft_fcflags_def}"
        sd_bigdft_ldflags="${sd_bigdft_ldflags_def}"
        sd_bigdft_libs="${sd_bigdft_libs_def}"
        test ! -z "${BIGDFT_CPPFLAGS}" && sd_bigdft_cppflags="${BIGDFT_CPPFLAGS}"
        test ! -z "${BIGDFT_FCFLAGS}" && sd_bigdft_fcflags="${BIGDFT_FCFLAGS}"
        test ! -z "${BIGDFT_LDFLAGS}" && sd_bigdft_ldflags="${BIGDFT_LDFLAGS}"
        test ! -z "${BIGDFT_LIBS}" && sd_bigdft_libs="${BIGDFT_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for BigDFT: '${sd_bigdft_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_bigdft_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_bigdft_init="esl"
    sd_bigdft_cppflags=""
    sd_bigdft_fcflags=""
    sd_bigdft_ldflags=""
    sd_bigdft_libs=""
  fi

  # Display configuration
  _SD_BIGDFT_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_bigdft_options)
  AC_SUBST(sd_bigdft_enable_def)
  AC_SUBST(sd_bigdft_policy)
  AC_SUBST(sd_bigdft_status)
  AC_SUBST(sd_bigdft_enable)
  AC_SUBST(sd_bigdft_init)
  AC_SUBST(sd_bigdft_ok)
  AC_SUBST(sd_bigdft_cppflags)
  AC_SUBST(sd_bigdft_fcflags)
  AC_SUBST(sd_bigdft_ldflags)
  AC_SUBST(sd_bigdft_libs)
  AC_SUBST(with_bigdft)

  # Clean-up
  unset tmp_bigdft_vars
]) # SD_BIGDFT_INIT


AC_DEFUN([SD_BIGDFT_DETECT], [
  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_bigdft_enable}" = "auto" -o "${sd_bigdft_enable}" = "yes"; then
    _SD_BIGDFT_CHECK_USE

    if test "${sd_bigdft_ok}" = "yes"; then
      if test "${sd_bigdft_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_bigdft_libs_def} ${sd_esl_bundle_libs}"
      else
        FCFLAGS="${FCFLAGS} ${sd_bigdft_fcflags}"
        LIBS="${sd_bigdft_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_bigdft_ldflags}"

      AC_DEFINE([HAVE_BIGDFT], 1,
        [Define to 1 if you have the BigDFT library.])

      _SD_BIGDFT_DUMP_CONFIG
    else
      if test "${sd_bigdft_status}" = "optional" -a \
              "${sd_bigdft_init}" = "def"; then
        sd_bigdft_enable="no"
        sd_bigdft_cppflags=""
        sd_bigdft_fcflags=""
        sd_bigdft_ldflags=""
        sd_bigdft_libs=""
      else
        AC_MSG_FAILURE([invalid BigDFT configuration])
      fi
    fi
  else
    sd_bigdft_enable="no"
    sd_bigdft_cppflags=""
    sd_bigdft_fcflags=""
    sd_bigdft_ldflags=""
    sd_bigdft_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_BIGDFT_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_bigdft_init}" = "esl"; then
    AC_MSG_NOTICE([will look for BigDFT in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_bigdft_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_linalg_cppflags} ${sd_libxc_cppflags} ${sd_bigdft_cppflags}"
    FCFLAGS="${FCFLAGS} ${sd_linalg_fcflags} ${sd_libxc_fcflags} ${sd_bigdft_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_linalg_ldflags} ${sd_libxc_ldflags} ${sd_bigdft_ldflags}"
    LIBS="${sd_bigdft_libs} ${sd_libxc_libs} ${sd_linalg_libs} ${LIBS}"
  fi

  # Check BigDFT Fortran API
  AC_MSG_CHECKING([whether the BigDFT Fortran interface works])
  for tmp_incs in "" "-I/usr/include"; do
    FCFLAGS="${FCFLAGS} ${tmp_incs}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
          use bigdft_api
          implicit none
          integer iproc
          type(input_variables) :: inputs
          type(atoms_data) :: atoms
          type(restart_objects) :: rst
          call init_restart_objects(iproc, inputs, atoms, rst)
      ]])], [sd_bigdft_ok="yes"], [sd_bigdft_ok="no"])
    AC_LANG_POP([Fortran])
    if test "${sd_bigdft_ok}" = "yes"; then
      test "${sd_sys_fcflags}" = "" && sd_sys_fcflags="${tmp_incs}"
      break
    fi
  done
  AC_MSG_RESULT([${sd_bigdft_ok}])
  unset tmp_incs

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_BIGDFT_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_BIGDFT_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_bigdft_invalid="no"
  if test "${sd_bigdft_enable_def}" != "auto" -a \
          "${sd_bigdft_enable_def}" != "no" -a \
          "${sd_bigdft_enable_def}" != "yes"; then
    case "${sd_bigdft_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_bigdft_enable_def = '${sd_bigdft_enable_def}'])
        ;;
      skip)
        tmp_bigdft_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_bigdft_enable_def = '${sd_bigdft_enable_def}'])
        tmp_bigdft_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_bigdft_invalid}" = "yes"; then
    if test "${sd_bigdft_status}" = "required"; then
      sd_bigdft_enable_def="yes"
    else
      sd_bigdft_enable_def="no"
    fi
    tmp_bigdft_invalid="no"
    AC_MSG_NOTICE([setting sd_bigdft_enable_def to '${sd_bigdft_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_bigdft_invalid="no"
  if test "${sd_bigdft_status}" = "implicit" -o \
          "${sd_bigdft_status}" = "required"; then
    if test "${sd_bigdft_enable}" = "no"; then
      case "${sd_bigdft_policy}" in
        fail)
          AC_MSG_ERROR([The BigDFT package is required and cannot be disabled
                  See https://launchpad.net/bigdft for details on how to
                  install it.])
          ;;
        skip)
          tmp_bigdft_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The BigDFT package is required and cannot be disabled])
          tmp_bigdft_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_bigdft_enable}" = "auto"; then
      AC_MSG_NOTICE([setting BigDFT trigger to yes])
      sd_bigdft_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_bigdft_invalid}" = "yes"; then
    case "${sd_bigdft_status}" in
      implicit|required)
        sd_bigdft_enable="yes"
        ;;
      optional)
        sd_bigdft_enable="no"
        ;;
    esac
    tmp_bigdft_invalid="no"
    AC_MSG_NOTICE([setting sd_bigdft_enable to '${sd_bigdft_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_bigdft_vars="${BIGDFT_CPPFLAGS}${BIGDFT_FCFLAGS}${BIGDFT_LDFLAGS}${BIGDFT_LIBS}"
  tmp_bigdft_invalid="no"
  if test ! -z "${tmp_bigdft_vars}" -a ! -z "${with_bigdft}"; then
    case "${sd_bigdft_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for BigDFT
                  Please use BIGDFT_FCFLAGS + BIGDFT_LIBS or --with-bigdft,
                  not both.])
        ;;
      skip)
        tmp_bigdft_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for BigDFT])
        tmp_bigdft_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_bigdft_vars}"; then
    sd_bigdft_enable="yes"
    sd_bigdft_init="env"
    if test "${tmp_bigdft_invalid}" = "yes"; then
      tmp_bigdft_invalid="no"
      AC_MSG_NOTICE([overriding --with-bigdft with BIGDFT_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_bigdft_status}" = "implicit"; then
    if test "${sd_bigdft_fcflags}" != ""; then
      sd_bigdft_fcflags=""
      AC_MSG_NOTICE([resetting BigDFT Fortran flags (implicit package)])
    fi
    if test "${sd_bigdft_ldflags}" != ""; then
      sd_bigdft_ldflags=""
      AC_MSG_NOTICE([resetting BigDFT linker flags (implicit package)])
    fi
    if test "${sd_bigdft_libs}" != ""; then
      sd_bigdft_libs=""
      AC_MSG_NOTICE([resetting BigDFT library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_bigdft_enable}" = "implicit"; then
    sd_bigdft_fcflags=""
    sd_bigdft_ldflags=""
    sd_bigdft_libs=""
  fi

  # Clean-up
  unset tmp_bigdft_invalid
  unset tmp_bigdft_vars
]) # _SD_BIGDFT_CHECK_CONFIG


AC_DEFUN([_SD_BIGDFT_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable BigDFT])
  AC_MSG_RESULT([${sd_bigdft_enable}])
  if test "${sd_bigdft_enable}" != "no"; then
    AC_MSG_CHECKING([how BigDFT parameters have been set])
    AC_MSG_RESULT([${sd_bigdft_init}])
    AC_MSG_CHECKING([for BigDFT C preprocessing flags])
    if test "${sd_bigdft_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_bigdft_cppflags}])
    fi
    AC_MSG_CHECKING([for BigDFT Fortran flags])
    if test "${sd_bigdft_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_bigdft_fcflags}])
    fi
    AC_MSG_CHECKING([for BigDFT linker flags])
    if test "${sd_bigdft_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_bigdft_ldflags}])
    fi
    AC_MSG_CHECKING([for BigDFT library flags])
    if test "${sd_bigdft_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_bigdft_libs}])
    fi
  fi
]) # _SD_BIGDFT_DUMP_CONFIG

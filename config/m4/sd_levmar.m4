## Copyright (C) 2019 Yann Pouillon

#
# Exchange-Correlation functionals library (Levmar)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_LEVMAR_INIT], [
  # Init
  sd_levmar_cppflags=""
  sd_levmar_cflags=""
  sd_levmar_ldflags=""
  sd_levmar_libs=""
  sd_levmar_enable=""
  sd_levmar_init="unknown"
  sd_levmar_ok="unknown"

  # Set adjustable parameters
  sd_levmar_options="$1"
  sd_levmar_libs_def="$2"
  sd_levmar_cppflags_def="$3"
  sd_levmar_cflags_def="$4"
  sd_levmar_cxxflags_def="$5"
  sd_levmar_fcflags_def="$6"
  sd_levmar_ldflags_def="$7"
  sd_levmar_enable_def=""
  sd_levmar_policy=""
  sd_levmar_status=""

  # Process options
  for kwd in ${sd_levmar_options}; do
    case "${kwd}" in
      auto)
        sd_levmar_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_levmar_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_levmar_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg Levmar option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_levmar_libs_def}" && sd_levmar_libs_def="-llevmar"
  test -z "${sd_levmar_policy}" && sd_levmar_policy="fail"
  test -z "${sd_levmar_status}" && sd_levmar_status="optional"
  test -z "${sd_levmar_enable_def}" && sd_levmar_enable_def="no"
  case "${sd_levmar_status}" in
    implicit|required)
      sd_levmar_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([levmar],
    [AS_HELP_STRING([--with-levmar],
      [Install prefix of the Levmar library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_levmar_enable="${withval}"
        sd_levmar_init="yon"
      else
        sd_levmar_enable="yes"
        sd_levmar_init="dir"
      fi],
    [ sd_levmar_enable="${sd_levmar_enable_def}"; sd_levmar_init="def"])

  # Declare environment variables
  AC_ARG_VAR([LEVMAR_CPPFLAGS], [C preprocessing flags for Levmar.])
  AC_ARG_VAR([LEVMAR_CFLAGS], [C flags for Levmar.])
  AC_ARG_VAR([LEVMAR_LDFLAGS], [Linker flags for Levmar.])
  AC_ARG_VAR([LEVMAR_LIBS], [Library flags for Levmar.])

  # Detect use of environment variables
  if test "${sd_levmar_enable}" = "yes" -o "${sd_levmar_enable}" = "auto"; then
    tmp_levmar_vars="${LEVMAR_CPPFLAGS}${LEVMAR_CFLAGS}${LEVMAR_LDFLAGS}${LEVMAR_LIBS}"
    if test "${sd_levmar_init}" = "def" -a ! -z "${tmp_levmar_vars}"; then
      sd_levmar_enable="yes"
      sd_levmar_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_LEVMAR_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_levmar_enable}" = "yes" -o "${sd_levmar_enable}" = "auto"; then

    # Set Levmar-specific flags
    case "${sd_levmar_init}" in

      def|yon)
        sd_levmar_cppflags="${sd_levmar_cppflags_def}"
        sd_levmar_cflags="${sd_levmar_cflags_def}"
        sd_levmar_ldflags="${sd_levmar_ldflags_def}"
        sd_levmar_libs="${sd_levmar_libs_def}"
        ;;

      dir)
        sd_levmar_cppflags="-I${with_levmar}/include"
        sd_levmar_cflags="${sd_levmar_cflags_def}"
        sd_levmar_ldflags="${sd_levmar_ldflags_def}"
        sd_levmar_libs="-L${with_levmar}/lib ${sd_levmar_libs_def}"
        ;;

      env)
        sd_levmar_cppflags="${sd_levmar_cppflags_def}"
        sd_levmar_cflags="${sd_levmar_cflags_def}"
        sd_levmar_ldflags="${sd_levmar_ldflags_def}"
        sd_levmar_libs="${sd_levmar_libs_def}"
        test ! -z "${LEVMAR_CPPFLAGS}" && sd_levmar_cppflags="${LEVMAR_CPPFLAGS}"
        test ! -z "${LEVMAR_CFLAGS}" && sd_levmar_cflags="${LEVMAR_CFLAGS}"
        test ! -z "${LEVMAR_LDFLAGS}" && sd_levmar_ldflags="${LEVMAR_LDFLAGS}"
        test ! -z "${LEVMAR_LIBS}" && sd_levmar_libs="${LEVMAR_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for Levmar: '${sd_levmar_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_levmar_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_levmar_init="esl"
    sd_levmar_cppflags=""
    sd_levmar_cflags=""
    sd_levmar_ldflags=""
    sd_levmar_libs=""
  fi

  # Display configuration
  _SD_LEVMAR_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_levmar_options)
  AC_SUBST(sd_levmar_enable_def)
  AC_SUBST(sd_levmar_policy)
  AC_SUBST(sd_levmar_status)
  AC_SUBST(sd_levmar_enable)
  AC_SUBST(sd_levmar_init)
  AC_SUBST(sd_levmar_ok)
  AC_SUBST(sd_levmar_cppflags)
  AC_SUBST(sd_levmar_cflags)
  AC_SUBST(sd_levmar_ldflags)
  AC_SUBST(sd_levmar_libs)
  AC_SUBST(with_levmar)

  # Clean-up
  unset tmp_levmar_vars
]) # SD_LEVMAR_INIT


AC_DEFUN([SD_LEVMAR_DETECT], [
  # Display configuration
  _SD_LEVMAR_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_levmar_enable}" = "auto" -o "${sd_levmar_enable}" = "yes"; then
    _SD_LEVMAR_CHECK_USE

    if test "${sd_levmar_ok}" = "yes"; then
      if test "${sd_levmar_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_levmar_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_levmar_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_levmar_ldflags}"

      AC_DEFINE([HAVE_LEVMAR], 1,
        [Define to 1 if you have the Levmar library.])
    else
      if test "${sd_levmar_status}" = "optional" -a \
              "${sd_levmar_init}" = "def"; then
        sd_levmar_enable="no"
        sd_levmar_cppflags=""
        sd_levmar_cflags=""
        sd_levmar_ldflags=""
        sd_levmar_libs=""
      else
        AC_MSG_FAILURE([invalid Levmar configuration])
      fi
    fi
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_LEVMAR_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_levmar_init}" = "esl"; then
    AC_MSG_NOTICE([will look for Levmar in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_levmar_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_levmar_cppflags}"
    CFLAGS="${CFLAGS} ${sd_levmar_cflags}"
    LDFLAGS="${LDFLAGS} ${sd_levmar_ldflags}"
    LIBS="${sd_levmar_libs} ${LIBS}"
  fi

  # Check Levmar API
  AC_MSG_CHECKING([whether the Levmar library works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#     include <levmar.h>

      void dfit_function(double *p, double *y, int m, int n, void *adata)
      {
        p = 0;
      }
    ]],
    [[
      int ret;
      int c_npoles = 1;
      int c_nvals  = 1;
      int nparam   = 1;

      double adata[1];
      double p[1];
      double yvals[1];

      double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

      ret=dlevmar_dif(dfit_function, p, yvals, nparam, c_nvals, 5000, \
        opts, info, NULL, NULL, (void *)&adata);
    ]])], [sd_levmar_ok="yes"], [sd_levmar_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_levmar_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_LEVMAR_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_LEVMAR_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_levmar_invalid="no"
  if test "${sd_levmar_enable_def}" != "auto" -a \
          "${sd_levmar_enable_def}" != "no" -a \
          "${sd_levmar_enable_def}" != "yes"; then
    case "${sd_levmar_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_levmar_enable_def = '${sd_levmar_enable_def}'])
        ;;
      skip)
        tmp_levmar_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_levmar_enable_def = '${sd_levmar_enable_def}'])
        tmp_levmar_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_levmar_invalid}" = "yes"; then
    if test "${sd_levmar_status}" = "required"; then
      sd_levmar_enable_def="yes"
    else
      sd_levmar_enable_def="no"
    fi
    tmp_levmar_invalid="no"
    AC_MSG_NOTICE([setting sd_levmar_enable_def to '${sd_levmar_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_levmar_invalid="no"
  if test "${sd_levmar_status}" = "implicit" -o \
          "${sd_levmar_status}" = "required"; then
    if test "${sd_levmar_enable}" = "no"; then
      case "${sd_levmar_policy}" in
        fail)
          AC_MSG_ERROR([The Levmar package is required and cannot be disabled
                  See http://users.ics.forth.gr/~lourakis/levmar/ for details
                  on how to install it.])
          ;;
        skip)
          tmp_levmar_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The Levmar package is required and cannot be disabled])
          tmp_levmar_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_levmar_enable}" = "auto"; then
      AC_MSG_NOTICE([setting Levmar trigger to yes])
      sd_levmar_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_levmar_invalid}" = "yes"; then
    case "${sd_levmar_status}" in
      implicit|required)
        sd_levmar_enable="yes"
        ;;
      optional)
        sd_levmar_enable="no"
        ;;
    esac
    tmp_levmar_invalid="no"
    AC_MSG_NOTICE([setting sd_levmar_enable to '${sd_levmar_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_levmar_vars="${LEVMAR_CFLAGS}${LEVMAR_LDFLAGS}${LEVMAR_LIBS}"
  tmp_levmar_invalid="no"
  if test ! -z "${tmp_levmar_vars}" -a ! -z "${with_levmar}"; then
    case "${sd_levmar_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for Levmar
                  Please use LEVMAR_CFLAGS + LEVMAR_LIBS or --with-levmar,
                  not both.])
        ;;
      skip)
        tmp_levmar_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for Levmar])
        tmp_levmar_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_levmar_vars}"; then
    sd_levmar_enable="yes"
    sd_levmar_init="env"
    if test "${tmp_levmar_invalid}" = "yes"; then
      tmp_levmar_invalid="no"
      AC_MSG_NOTICE([overriding --with-levmar with LEVMAR_{CPPFLAGS,CFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_levmar_status}" = "implicit"; then
    if test "${sd_levmar_ldflags}" != ""; then
      sd_levmar_ldflags=""
      AC_MSG_NOTICE([resetting Levmar linker flags (implicit package)])
    fi
    if test "${sd_levmar_libs}" != ""; then
      sd_levmar_libs=""
      AC_MSG_NOTICE([resetting Levmar library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_levmar_enable}" = "implicit"; then
    sd_levmar_cppflags=""
    sd_levmar_cflags=""
    sd_levmar_ldflags=""
    sd_levmar_libs=""
  fi

  # Clean-up
  unset tmp_levmar_invalid
  unset tmp_levmar_vars
]) # _SD_LEVMAR_CHECK_CONFIG


AC_DEFUN([_SD_LEVMAR_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable Levmar])
  AC_MSG_RESULT([${sd_levmar_enable}])
  if test "${sd_levmar_enable}" != "no"; then
    AC_MSG_CHECKING([how Levmar parameters have been set])
    AC_MSG_RESULT([${sd_levmar_init}])
    AC_MSG_CHECKING([for Levmar C preprocessing flags])
    if test "${sd_levmar_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_levmar_cppflags}])
    fi
    AC_MSG_CHECKING([for Levmar C flags])
    if test "${sd_levmar_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_levmar_cflags}])
    fi
    AC_MSG_CHECKING([for Levmar linker flags])
    if test "${sd_levmar_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_levmar_ldflags}])
    fi
    AC_MSG_CHECKING([for Levmar library flags])
    if test "${sd_levmar_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_levmar_libs}])
    fi
  fi
]) # _SD_LEVMAR_DUMP_CONFIG

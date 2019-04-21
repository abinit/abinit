## Copyright (C) 2019 Yann Pouillon

#
# Toolbox for Research on Interacting Quantum Systems (TRIQS)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_TRIQS_INIT], [
  # Init
  sd_triqs_cppflags=""
  sd_triqs_cflags=""
  sd_triqs_cxxflags=""
  sd_triqs_ldflags=""
  sd_triqs_libs=""
  sd_triqs_enable=""
  sd_triqs_init="unknown"
  sd_triqs_ok="unknown"
  sd_triqs_api_version="unknown"

  # Set adjustable parameters
  sd_triqs_options="$1"
  sd_triqs_libs_def="$2"
  sd_triqs_cppflags_def="$3"
  sd_triqs_cflags_def="$4"
  sd_triqs_cxxflags_def="$5"
  sd_triqs_fcflags_def="$6"
  sd_triqs_ldflags_def="$7"
  sd_triqs_enable_def=""
  sd_triqs_policy=""
  sd_triqs_status=""

  # Process options
  for kwd in ${sd_triqs_options}; do
    case "${kwd}" in
      auto)
        sd_triqs_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_triqs_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_triqs_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg TRIQS option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_triqs_libs_def}" && sd_triqs_libs_def="-ltriqs -lcthyb_c"
  test -z "${sd_triqs_policy}" && sd_triqs_policy="fail"
  test -z "${sd_triqs_status}" && sd_triqs_status="optional"
  test -z "${sd_triqs_enable_def}" && sd_triqs_enable_def="no"
  case "${sd_triqs_status}" in
    implicit|required)
      sd_triqs_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([triqs],
    [AS_HELP_STRING([--with-triqs],
      [Install prefix of the TRIQS library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_triqs_enable="${withval}"
        sd_triqs_init="yon"
      else
        sd_triqs_enable="yes"
        sd_triqs_init="dir"
      fi],
    [ sd_triqs_enable="${sd_triqs_enable_def}"; sd_triqs_init="def"])

  # Declare environment variables
  AC_ARG_VAR([TRIQS_CPPFLAGS], [C preprocessing flags for TRIQS.])
  AC_ARG_VAR([TRIQS_CFLAGS], [C flags for TRIQS.])
  AC_ARG_VAR([TRIQS_CXXFLAGS], [C++ flags for TRIQS.])
  AC_ARG_VAR([TRIQS_LDFLAGS], [Linker flags for TRIQS.])
  AC_ARG_VAR([TRIQS_LIBS], [Library flags for TRIQS.])

  # Detect use of environment variables
  if test "${sd_triqs_enable}" = "yes" -o "${sd_triqs_enable}" = "auto"; then
    tmp_triqs_vars="${TRIQS_CPPFLAGS}${TRIQS_CFLAGS}${TRIQS_LDFLAGS}${TRIQS_LIBS}"
    if test "${sd_triqs_enable_cxx}" = "yes"; then
      tmp_triqs_vars="${tmp_triqs_vars}${TRIQS_CXXFLAGS}"
    fi
    if test "${sd_triqs_init}" = "def" -a ! -z "${tmp_triqs_vars}"; then
      sd_triqs_enable="yes"
      sd_triqs_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_TRIQS_CHECK_CONFIG
  fi
  # Adjust configuration depending on init type
  if test "${sd_triqs_enable}" = "yes" -o "${sd_triqs_enable}" = "auto"; then

    # Set TRIQS-specific flags
    case "${sd_triqs_init}" in

      def|yon)
        sd_triqs_cppflags="${sd_triqs_cppflags_def}"
        sd_triqs_cflags="${sd_triqs_cflags_def}"
        test "${sd_triqs_enable_cxx}" = "yes" && \
          sd_triqs_cxxflags="${sd_triqs_cxxflags_def}"
        sd_triqs_ldflags="${sd_triqs_ldflags_def}"
        sd_triqs_libs="${sd_triqs_libs_def}"
        ;;

      dir)
        sd_triqs_cppflags="-I${with_triqs}/include"
        sd_triqs_cflags="${sd_triqs_cflags_def}"
        test "${sd_triqs_enable_cxx}" = "yes" && \
          sd_triqs_cxxflags="${sd_triqs_cxxflags_def}"
        sd_triqs_ldflags="${sd_triqs_ldflags_def}"
        sd_triqs_libs="-L${with_triqs}/lib ${sd_triqs_libs_def}"
        ;;

      env)
        sd_triqs_cppflags="${sd_triqs_cppflags_def}"
        sd_triqs_cflags="${sd_triqs_cflags_def}"
        test "${sd_triqs_enable_cxx}" = "yes" && \
          sd_triqs_cxxflags="${sd_triqs_cxxflags_def}"
        sd_triqs_ldflags="${sd_triqs_ldflags_def}"
        sd_triqs_libs="${sd_triqs_libs_def}"
        test ! -z "${TRIQS_CPPFLAGS}" && sd_triqs_cppflags="${TRIQS_CPPFLAGS}"
        test ! -z "${TRIQS_CFLAGS}" && sd_triqs_cflags="${TRIQS_CFLAGS}"
        if test "${sd_triqs_enable_cxx}" = "yes"; then
          test ! -z "${TRIQS_CXXFLAGS}" && sd_triqs_cxxflags="${TRIQS_CXXFLAGS}"
        fi
        test ! -z "${TRIQS_LDFLAGS}" && sd_triqs_ldflags="${TRIQS_LDFLAGS}"
        test ! -z "${TRIQS_LIBS}" && sd_triqs_libs="${TRIQS_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for TRIQS: '${sd_triqs_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_triqs_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_triqs_init="esl"
    sd_triqs_cppflags=""
    sd_triqs_cflags=""
    sd_triqs_cxxflags=""
    sd_triqs_ldflags=""
    sd_triqs_libs=""
  fi

  # Display configuration
  _SD_TRIQS_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_triqs_options)
  AC_SUBST(sd_triqs_enable_def)
  AC_SUBST(sd_triqs_enable_cxx)
  AC_SUBST(sd_triqs_policy)
  AC_SUBST(sd_triqs_status)
  AC_SUBST(sd_triqs_enable)
  AC_SUBST(sd_triqs_init)
  AC_SUBST(sd_triqs_ok)
  AC_SUBST(sd_triqs_cppflags)
  AC_SUBST(sd_triqs_cflags)
  AC_SUBST(sd_triqs_ldflags)
  AC_SUBST(sd_triqs_libs)
  AC_SUBST(with_triqs)

  # Clean-up
  unset tmp_triqs_vars
]) # SD_TRIQS_INIT


AC_DEFUN([SD_TRIQS_DETECT], [
  # Display configuration
  _SD_TRIQS_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_triqs_enable}" = "auto" -o "${sd_triqs_enable}" = "yes"; then
    _SD_TRIQS_CHECK_USE

    if test "${sd_triqs_ok}" = "yes"; then
      if test "${sd_triqs_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_triqs_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_triqs_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_triqs_ldflags}"

      AC_DEFINE([HAVE_TRIQS], 1,
        [Define to 1 if you have the TRIQS library.])

      case "${sd_triqs_api_version}" in
        1.4)
          AC_DEFINE([HAVE_TRIQS_v1_4], 1,
            [Define to 1 if you have the TRIQS 1.4 libraries.])
          ;;
        2.0)
          AC_DEFINE([HAVE_TRIQS_v2_0], 1,
            [Define to 1 if you have the TRIQS 1.4 libraries.])
          ;;
        *)
          AC_MSG_ERROR([TRIQS API version ${sd_triqs_api_version} not implemented in the build system])
          ;;
      esac
    else
      if test "${sd_triqs_status}" = "optional" -a \
              "${sd_triqs_init}" = "def"; then
        sd_triqs_enable="no"
        sd_triqs_cppflags=""
        sd_triqs_cflags=""
        sd_triqs_cxxflags=""
        sd_triqs_ldflags=""
        sd_triqs_libs=""
      else
        AC_MSG_FAILURE([invalid TRIQS configuration])
      fi
    fi


  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_TRIQS_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_triqs_init}" = "esl"; then
    AC_MSG_NOTICE([will look for TRIQS in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_triqs_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_triqs_cppflags}"
    CFLAGS="${CFLAGS} ${sd_triqs_cflags}"
    LDFLAGS="${LDFLAGS} ${sd_triqs_ldflags}"
    LIBS="${sd_triqs_libs} ${LIBS}"
  fi

  # Check TRIQS C++ API
  AC_MSG_CHECKING([whether the TRIQS library works])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#       include <triqs/gfs.hpp>
        using namespace triqs::gfs;
    ]],
    [[
      triqs::clef::placeholder<0> iw_;
      double beta  = 10;
      auto iw_mesh = gf_mesh<imfreq>{beta, Fermion, 100};
      auto G = gf<imfreq, scalar_valued>{iw_mesh, make_shape()};
      G(iw_) << 1.0 / iw_ + 2.0 / iw_ / iw_ + 3.0 / iw_ / iw_ / iw_;
      auto known_moments = array<dcomplex, 1>{0.0, 1.0};
      auto [tail, err] = fit_tail(G, known_moments);
    ]])], [sd_triqs_ok="yes"; sd_triqs_api_version="2.0"], [sd_triqs_ok="no"])
  AC_LANG_POP([C++])
  AC_MSG_RESULT([${sd_triqs_ok}])

  # Check old TRIQS C++ API
  if test "${sd_triqs_ok}" != "yes"; then
    AC_MSG_CHECKING([whether the TRIQS library has an old API])
    AC_LANG_PUSH([C++])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#       include <triqs/gfs.hpp>
        using namespace triqs::gfs;
        using triqs::clef::placeholder;
      ]],
      [[
        double beta = 1;
        int nw      = 100;
        auto g      = gf<imfreq>{{beta, Fermion, nw}, {1, 1}};
        placeholder<0> w_;
        g(w_) << 1 / (w_ - 3);
      ]])], [sd_triqs_ok="yes"; sd_triqs_api_version="1.4"], [sd_triqs_ok="no"])
    AC_LANG_POP([C++])
    AC_MSG_RESULT([${sd_triqs_ok}])
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_TRIQS_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_TRIQS_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_triqs_invalid="no"
  if test "${sd_triqs_enable_def}" != "auto" -a \
          "${sd_triqs_enable_def}" != "no" -a \
          "${sd_triqs_enable_def}" != "yes"; then
    case "${sd_triqs_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_triqs_enable_def = '${sd_triqs_enable_def}'])
        ;;
      skip)
        tmp_triqs_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_triqs_enable_def = '${sd_triqs_enable_def}'])
        tmp_triqs_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_triqs_invalid}" = "yes"; then
    if test "${sd_triqs_status}" = "required"; then
      sd_triqs_enable_def="yes"
    else
      sd_triqs_enable_def="no"
    fi
    tmp_triqs_invalid="no"
    AC_MSG_NOTICE([setting sd_triqs_enable_def to '${sd_triqs_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_triqs_invalid="no"
  if test "${sd_triqs_status}" = "implicit" -o \
          "${sd_triqs_status}" = "required"; then
    if test "${sd_triqs_enable}" = "no"; then
      case "${sd_triqs_policy}" in
        fail)
          AC_MSG_ERROR([The TRIQS package is required and cannot be disabled
                  See https://triqs.github.io/ for details on how to
                  install it.])
          ;;
        skip)
          tmp_triqs_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The TRIQS package is required and cannot be disabled])
          tmp_triqs_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_triqs_enable}" = "auto"; then
      AC_MSG_NOTICE([setting TRIQS trigger to yes])
      sd_triqs_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_triqs_invalid}" = "yes"; then
    case "${sd_triqs_status}" in
      implicit|required)
        sd_triqs_enable="yes"
        ;;
      optional)
        sd_triqs_enable="no"
        ;;
    esac
    tmp_triqs_invalid="no"
    AC_MSG_NOTICE([setting sd_triqs_enable to '${sd_triqs_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_triqs_vars="${TRIQS_CXXFLAGS}${TRIQS_LDFLAGS}${TRIQS_LIBS}"
  tmp_triqs_invalid="no"
  if test ! -z "${tmp_triqs_vars}" -a ! -z "${with_triqs}"; then
    case "${sd_triqs_policy}" in
      fail)
        AC_MSG_ERROR([conflicting option settings for TRIQS
                  Please use TRIQS_CXXFLAGS + TRIQS_LIBS or --with-triqs,
                  not both.])
        ;;
      skip)
        tmp_triqs_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for TRIQS])
        tmp_triqs_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_triqs_vars}"; then
    sd_triqs_enable="yes"
    sd_triqs_init="env"
    if test "${tmp_triqs_invalid}" = "yes"; then
      tmp_triqs_invalid="no"
      AC_MSG_NOTICE([overriding --with-triqs with TRIQS_{CXXFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_triqs_status}" = "implicit"; then
    if test "${sd_triqs_ldflags}" != ""; then
      sd_triqs_ldflags=""
      AC_MSG_NOTICE([resetting TRIQS linker flags (implicit package)])
    fi
    if test "${sd_triqs_libs}" != ""; then
      sd_triqs_libs=""
      AC_MSG_NOTICE([resetting TRIQS library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_triqs_enable}" = "implicit"; then
    sd_triqs_ldflags=""
    sd_triqs_libs=""
  fi

  # Clean-up
  unset tmp_triqs_invalid
  unset tmp_triqs_vars
]) # _SD_TRIQS_CHECK_CONFIG


AC_DEFUN([_SD_TRIQS_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable TRIQS])
  AC_MSG_RESULT([${sd_triqs_enable}])
  if test "${sd_triqs_enable}" != "no"; then
    AC_MSG_CHECKING([how TRIQS parameters have been set])
    AC_MSG_RESULT([${sd_triqs_init}])
    AC_MSG_CHECKING([for TRIQS C preprocessing flags])
    if test "${sd_triqs_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_triqs_cppflags}])
    fi
    AC_MSG_CHECKING([for TRIQS C flags])
    if test "${sd_triqs_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_triqs_cflags}])
    fi
    if test "${sd_triqs_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([for TRIQS C++ flags])
      if test "${sd_triqs_cxxflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_triqs_cxxflags}])
      fi
    fi
    AC_MSG_CHECKING([for TRIQS linker flags])
    if test "${sd_triqs_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_triqs_ldflags}])
    fi
    AC_MSG_CHECKING([for TRIQS library flags])
    if test "${sd_triqs_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_triqs_libs}])
    fi
  fi
]) # _SD_TRIQS_DUMP_CONFIG

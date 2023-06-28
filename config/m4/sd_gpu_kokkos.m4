## Copyright (C) 2019-2022 ABINIT group (Yann Pouillon)

#
# Kokkos (Kokkos "core" libraries)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_KOKKOS_INIT], [
  # Init
  sd_kokkos_cppflags=""
  sd_kokkos_cxxflags=""
  sd_kokkos_ldflags=""
  sd_kokkos_libs=""
  sd_kokkos_enable=""
  sd_kokkos_init="unknown"
  sd_kokkos_cxx_ok="unknown"
  sd_kokkos_ok="unknown"

  # Set adjustable parameters
  sd_kokkos_options="$1"
  sd_kokkos_libs_def="$2"
  sd_kokkos_cppflags_def="$3"
  sd_kokkos_cxxflags_def="$4"
  sd_kokkos_ldflags_def="$5"
  sd_kokkos_enable_def=""
  sd_kokkos_policy=""
  sd_kokkos_status=""

  # Process options
  for kwd in ${sd_kokkos_options}; do
    case "${kwd}" in
      auto)
        sd_kokkos_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_kokkos_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_kokkos_policy="${kwd}"
        ;;
      mandatory)
        sd_kokkos_enable="yes"
        sd_kokkos_enable_def="yes"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg Kokkos core libraries option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_kokkos_libs_def}" && sd_kokkos_libs_def="-lkokkoscore"
  test -z "${sd_kokkos_policy}" && sd_kokkos_policy="fail"
  test -z "${sd_kokkos_status}" && sd_kokkos_status="optional"
  test -z "${sd_kokkos_enable_def}" && sd_kokkos_enable_def="no"
  case "${sd_kokkos_status}" in
    implicit|required)
      sd_kokkos_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([kokkos],
    [AS_HELP_STRING([--with-kokkos],
      [Install prefix of the Kokkos core libraries library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_kokkos_enable="${withval}"
        sd_kokkos_init="yon"
      else
        sd_kokkos_enable="yes"
        sd_kokkos_init="dir"
      fi],
    [ sd_kokkos_enable="${sd_kokkos_enable_def}"; sd_kokkos_init="def"])

  # Declare environment variables
  AC_ARG_VAR([KOKKOS_CPPFLAGS], [C preprocessing flags for Kokkos core libraries.])
  AC_ARG_VAR([KOKKOS_CXXFLAGS], [C++ flags for Kokkos core libraries.])
  AC_ARG_VAR([KOKKOS_LDFLAGS], [Linker flags for Kokkos core libraries.])
  AC_ARG_VAR([KOKKOS_LIBS], [Library flags for Kokkos core libraries.])

  # Detect use of environment variables
  if test "${sd_kokkos_enable}" = "yes" -o "${sd_kokkos_enable}" = "auto"; then
    tmp_kokkos_vars="${KOKKOS_CPPFLAGS}${KOKKOS_CXXFLAGS}${KOKKOS_LDFLAGS}${KOKKOS_LIBS}"
    if test "${sd_kokkos_init}" = "def" -a ! -z "${tmp_kokkos_vars}"; then
      sd_kokkos_enable="yes"
      sd_kokkos_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_KOKKOS_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_kokkos_enable}" = "yes" -o "${sd_kokkos_enable}" = "auto"; then

    # Set Kokkos core libraries-specific flags
    case "${sd_kokkos_init}" in

      def|yon)
        sd_kokkos_cppflags="${sd_kokkos_cppflags_def}"
        sd_kokkos_cxxflags="${sd_kokkos_cxxflags_def}"
        sd_kokkos_ldflags="${sd_kokkos_ldflags_def}"
        sd_kokkos_libs="${sd_kokkos_libs_def}"
        ;;

      dir)
        sd_kokkos_cppflags="-I${with_kokkos}/include"
        sd_kokkos_cxxflags="${sd_kokkos_cxxflags_def} -I${with_kokkos}/include"
        sd_kokkos_ldflags="${sd_kokkos_ldflags_def}"
        sd_kokkos_libs="-L${with_kokkos}/lib64 ${sd_kokkos_libs_def}"
        ;;

      env)
        sd_kokkos_cppflags="${sd_kokkos_cppflags_def}"
        sd_kokkos_cxxflags="${sd_kokkos_cxxflags_def}"
        sd_kokkos_ldflags="${sd_kokkos_ldflags_def}"
        sd_kokkos_libs="${sd_kokkos_libs_def}"
        test ! -z "${KOKKOS_CPPFLAGS}" && sd_kokkos_cppflags="${KOKKOS_CPPFLAGS}"
        test ! -z "${KOKKOS_LDFLAGS}" && sd_kokkos_ldflags="${KOKKOS_LDFLAGS}"
        test ! -z "${KOKKOS_LIBS}" && sd_kokkos_libs="${KOKKOS_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for Kokkos core libraries: '${sd_kokkos_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_kokkos_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_kokkos_init="esl"
    sd_kokkos_cppflags=""
    sd_kokkos_cxxflags=""
    sd_kokkos_ldflags=""
    sd_kokkos_libs=""
  fi

  AM_CONDITIONAL(DO_BUILD_16_KOKKOS_TOOLBOX,[test "${sd_kokkos_enable}" = "yes"])
  AM_CONDITIONAL(DO_BUILD_44_MANAGE_KOKKOS,[test "${sd_kokkos_enable}" = "yes"])

  # Display configuration
  _SD_KOKKOS_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_kokkos_options)
  AC_SUBST(sd_kokkos_enable_def)
  AC_SUBST(sd_kokkos_policy)
  AC_SUBST(sd_kokkos_status)
  AC_SUBST(sd_kokkos_enable)
  AC_SUBST(sd_kokkos_init)
  AC_SUBST(sd_kokkos_ok)
  AC_SUBST(sd_kokkos_cppflags)
  AC_SUBST(sd_kokkos_cxxflags)
  AC_SUBST(sd_kokkos_ldflags)
  AC_SUBST(sd_kokkos_libs)
  AC_SUBST(with_kokkos)

  # Clean-up
  unset tmp_kokkos_vars
]) # SD_KOKKOS_INIT


AC_DEFUN([SD_KOKKOS_DETECT], [
  # Display configuration
  _SD_KOKKOS_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_kokkos_enable}" = "auto" -o "${sd_kokkos_enable}" = "yes"; then
    _SD_KOKKOS_CHECK_USE

    if test "${sd_kokkos_ok}" = "yes"; then
      if test "${sd_kokkos_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_kokkos_libs_def} ${sd_esl_bundle_libs}"
      else
        CPPFLAGS="${CPPFLAGS} ${sd_kokkos_cppflags}"
        CXXFLAGS="${CXXFLAGS} ${sd_kokkos_cxxflags}"
        CXXFLAGS="${CXXFLAGS} ${sd_kokkos_cxxflags}"
        LIBS="${sd_kokkos_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_kokkos_ldflags}"

      AC_DEFINE([HAVE_KOKKOS], 1,
        [Define to 1 if you have the Kokkos core libraries library.])
    else
      if test "${sd_kokkos_status}" = "optional" -a \
              "${sd_kokkos_init}" = "def"; then
        sd_kokkos_enable="no"
        sd_kokkos_cppflags=""
        sd_kokkos_cxxflags=""
        sd_kokkos_ldflags=""
        sd_kokkos_libs=""
      else
        if test "${sd_kokkos_policy}" = "fail"; then
              AC_MSG_FAILURE([invalid Kokkos core libraries configuration])
        else
              AC_MSG_WARN([invalid Kokkos core libraries configuration])
        fi
      fi
    fi
  else
    sd_kokkos_enable="no"
    sd_kokkos_cppflags=""
    sd_kokkos_cxxflags=""
    sd_kokkos_ldflags=""
    sd_kokkos_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_KOKKOS_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_kokkos_init}" = "esl"; then
    AC_MSG_NOTICE([will look for Kokkos core libraries in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_kokkos_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_kokkos_cppflags}"
    CXXFLAGS="${CXXFLAGS} ${sd_kokkos_cxxflags}"
    LDFLAGS="${LDFLAGS} ${sd_kokkos_ldflags}"
    LIBS="${sd_kokkos_libs} ${LIBS}"
  fi

  # Check Kokkos core libraries C++ API
  AC_MSG_CHECKING([whether the Kokkos core library works])
#  AC_LANG_PUSH([C++])
#  AC_LINK_IFELSE([AC_LANG_PROGRAM(
#    [[
##include <Kokkos_Core.hpp>
#    ]],
#    [[
#      Kokkos::abort("Testing purposes");
#    ]])], [sd_kokkos_cxx_ok="yes"], [sd_kokkos_cxx_ok="no"])
#  AC_LANG_POP([C++])
#  FIXME Kokkos is actually used with a wrapper on NVCC and C++ compiler.
#  It's complicate, hence we only check for library presence and run no link check (WCGW?).
  sd_kokkos_cxx_ok="no"
  if test -e "${with_kokkos}/lib64/libkokkoscore.${abi_so_ext}"; then
    sd_kokkos_cxx_ok="yes"
  fi
  AC_MSG_RESULT([${sd_kokkos_cxx_ok}])

  # Combine the available results
  sd_kokkos_ok="no"
  if test "${sd_kokkos_cxx_ok}" = "yes"; then
    sd_kokkos_ok="yes"
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_KOKKOS_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_KOKKOS_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_kokkos_invalid="no"
  if test "${sd_kokkos_enable_def}" != "auto" -a \
          "${sd_kokkos_enable_def}" != "no" -a \
          "${sd_kokkos_enable_def}" != "yes"; then
    case "${sd_kokkos_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_kokkos_enable_def = '${sd_kokkos_enable_def}'])
        ;;
      skip)
        tmp_kokkos_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_kokkos_enable_def = '${sd_kokkos_enable_def}'])
        tmp_kokkos_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_kokkos_invalid}" = "yes"; then
    if test "${sd_kokkos_status}" = "required"; then
      sd_kokkos_enable_def="yes"
    else
      sd_kokkos_enable_def="no"
    fi
    tmp_kokkos_invalid="no"
    AC_MSG_NOTICE([setting sd_kokkos_enable_def to '${sd_kokkos_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_kokkos_invalid="no"
  if test "${sd_kokkos_status}" = "implicit" -o \
          "${sd_kokkos_status}" = "required"; then
    if test "${sd_kokkos_enable}" = "no"; then
      case "${sd_kokkos_policy}" in
        fail)
          AC_MSG_ERROR([The Kokkos core libraries package is required and cannot be disabled
                  See https://tddft.org/programs/kokkos/ for details on how to
                  install it.])
          ;;
        skip)
          tmp_kokkos_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The Kokkos core libraries package is required and cannot be disabled])
          tmp_kokkos_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_kokkos_enable}" = "auto"; then
      AC_MSG_NOTICE([setting Kokkos core libraries trigger to yes])
      sd_kokkos_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_kokkos_invalid}" = "yes"; then
    case "${sd_kokkos_status}" in
      implicit|required)
        sd_kokkos_enable="yes"
        ;;
      optional)
        sd_kokkos_enable="no"
        ;;
    esac
    tmp_kokkos_invalid="no"
    AC_MSG_NOTICE([setting sd_kokkos_enable to '${sd_kokkos_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_kokkos_vars="${KOKKOS_CPPFLAGS}${KOKKOS_CFLAGS}${KOKKOS_CXXFLAGS}${KOKKOS_LDFLAGS}${KOKKOS_LIBS}"
  tmp_kokkos_invalid="no"
  if test ! -z "${tmp_kokkos_vars}" -a ! -z "${with_kokkos}"; then
    case "${sd_kokkos_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for Kokkos core libraries
                  Please use KOKKOS_CXXFLAGS + KOKKOS_LIBS or --with-kokkos,
                  not both.])
        ;;
      skip)
        tmp_kokkos_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for Kokkos core libraries])
        tmp_kokkos_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_kokkos_vars}"; then
    sd_kokkos_enable="yes"
    sd_kokkos_init="env"
    if test "${tmp_kokkos_invalid}" = "yes"; then
      tmp_kokkos_invalid="no"
      AC_MSG_NOTICE([overriding --with-kokkos with KOKKOS_{CXXFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_kokkos_status}" = "implicit"; then
    if test "${sd_kokkos_cxxflags}" != ""; then
      sd_kokkos_cxxflags=""
      AC_MSG_NOTICE([resetting Kokkos core libraries C++ flags (implicit package)])
    fi
    if test "${sd_kokkos_ldflags}" != ""; then
      sd_kokkos_ldflags=""
      AC_MSG_NOTICE([resetting Kokkos core libraries linker flags (implicit package)])
    fi
    if test "${sd_kokkos_libs}" != ""; then
      sd_kokkos_libs=""
      AC_MSG_NOTICE([resetting Kokkos core libraries library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_kokkos_enable}" = "implicit"; then
    sd_kokkos_cxxflags=""
    sd_kokkos_ldflags=""
    sd_kokkos_libs=""
  fi

  # Clean-up
  unset tmp_kokkos_invalid
  unset tmp_kokkos_vars
]) # _SD_KOKKOS_CHECK_CONFIG


AC_DEFUN([_SD_KOKKOS_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable Kokkos core libraries])
  AC_MSG_RESULT([${sd_kokkos_enable}])
  if test "${sd_kokkos_enable}" != "no"; then
    AC_MSG_CHECKING([how Kokkos core libraries parameters have been set])
    AC_MSG_RESULT([${sd_kokkos_init}])
    AC_MSG_CHECKING([for Kokkos core libraries C preprocessing flags])
    if test "${sd_kokkos_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_kokkos_cppflags}])
    fi
    AC_MSG_CHECKING([for Kokkos core libraries C++ flags])
    if test "${sd_kokkos_cxxflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_kokkos_cxxflags}])
    fi
    AC_MSG_CHECKING([for Kokkos core libraries linker flags])
    if test "${sd_kokkos_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_kokkos_ldflags}])
    fi
    AC_MSG_CHECKING([for Kokkos core libraries library flags])
    if test "${sd_kokkos_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_kokkos_libs}])
    fi
  fi
]) # _SD_KOKKOS_DUMP_CONFIG

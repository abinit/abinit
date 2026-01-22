## Copyright (C) 2019-2026 ABINIT group (Yann Pouillon. MTorrent)

#
# openMP support for Steredeg
#


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_OPENMP_INIT], [
  # Init
  sd_openmp_cppflags=""
  sd_openmp_cflags=""
  sd_openmp_fcflags=""
  sd_openmp_ldflags=""
  sd_openmp_libs=""
  sd_openmp_cc=""
  sd_openmp_cc_ok="unknown"
  sd_openmp_cxx=""
  sd_openmp_cxx_ok="unknown"
  sd_openmp_fc=""
  sd_openmp_fc_ok="unknown"
  sd_openmp_enable=""
  sd_openmp_init="unknown"
  sd_openmp_ok="unknown"

  # Set adjustable parameters
  sd_openmp_options="$1"
  sd_openmp_libs_def="$2"
  sd_openmp_cppflags_def="$3"
  sd_openmp_cflags_def="$4"
  sd_openmp_cxxflags_def="$5"
  sd_openmp_fcflags_def="$6"
  sd_openmp_ldflags_def="$7"
  sd_openmp_enable_def=""
  sd_openmp_enable_cc=""
  sd_openmp_enable_cxx=""
  sd_openmp_enable_fc=""
  sd_openmp_policy=""
  sd_openmp_status=""

  # Process options
  for kwd in ${sd_openmp_options}; do
    case "${kwd}" in
      auto|no|yes)
        sd_openmp_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_openmp_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_openmp_policy="${kwd}"
        ;;
      no-cc)
        sd_openmp_enable_cc="no"
        ;;
      no-cxx)
        sd_openmp_enable_cxx="no"
        ;;
      no-fortran)
        sd_openmp_enable_fc="no"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg openMP option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_openmp_policy}" && sd_openmp_policy="warn"
  test -z "${sd_openmp_status}" && sd_openmp_status="optional"
  test "${sd_openmp_enable_cc}" = "" && sd_openmp_enable_cc="yes"
  test "${sd_openmp_enable_cxx}" = "" && sd_openmp_enable_cxx="yes"
  test "${sd_openmp_enable_fc}" = "" && sd_openmp_enable_fc="yes"
  test -z "${sd_openmp_enable_def}" && sd_openmp_enable_def="no"
  case "${sd_openmp_status}" in
    implicit|required)
      sd_openmp_enable_def="yes"
      ;;
  esac

  # Declare configure option
  AC_ARG_ENABLE(openmp,
    AS_HELP_STRING([--enable-openmp],
      [Activate support for OpenMP (default: no)]),
    [sd_openmp_enable="${enableval}"; sd_openmp_init="yon"],
    [sd_openmp_enable="no"; sd_openmp_init="def"])

  # Declare environment variables
  #AC_ARG_VAR([OPENMP_CPPFLAGS], [C preprocessing flags for openMP.])
  #AC_ARG_VAR([OPENMP_CFLAGS], [C flags for openMP.])
  #AC_ARG_VAR([OPENMP_CXXFLAGS], [C++ flags for openMP.])
  #AC_ARG_VAR([OPENMP_FCFLAGS], [Fortran flags for openMP.])
  #AC_ARG_VAR([OPENMP_FFLAGS], [Fortran flags for openMP (better use OPENMP_FCFLAGS).])
  #AC_ARG_VAR([OPENMP_LDFLAGS], [Linker flags for openMP.])
  #AC_ARG_VAR([OPENMP_LIBS], [Library flags for openMP.])

  # Detect use of environment variables
  if test "${sd_openmp_enable}" = "yes" -o "${sd_openmp_enable}" = "auto"; then
    tmp_openmp_vars="${OPENMP_CPPFLAGS}${OPENMP_LDFLAGS}${OPENMP_LIBS}"
    if test "${sd_openmp_enable_cc}" = "yes"; then
      test ! -z "${OPENMP_CFLAGS}" && tmp_openmp_vars="${tmp_openmp_vars}${OPENMP_CFLAGS}"
    fi
    if test "${sd_openmp_enable_cxx}" = "yes"; then
      test ! -z "${OPENMP_CXXFLAGS}" && tmp_openmp_vars="${tmp_openmp_vars}${OPENMP_CXXFLAGS}"
    fi
    if test "${sd_openmp_enable_fc}" = "yes"; then
      test ! -z "${OPENMP_FFLAGS}" && tmp_openmp_vars="${tmp_openmp_vars}${OPENMP_FFLAGS}"
      test ! -z "${OPENMP_FCFLAGS}" && tmp_openmp_vars="${tmp_openmp_vars}${OPENMP_FCFLAGS}"
    fi
    if test "${sd_openmp_init}" = "def" -o "${sd_openmp_init}" = "yon"; then
      if test "${tmp_openmp_vars}" != ""; then
        sd_openmp_enable="yes"
        sd_openmp_init="env"
      fi
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
    _SD_OPENMP_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_openmp_enable}" = "yes" -o "${sd_openmp_enable}" = "auto"; then

    # Set openMP-specific flags
    case "${sd_openmp_init}" in

      def|yon)
        sd_openmp_cppflags="${sd_openmp_cppflags_def}"
        test "${sd_openmp_enable_cc}" = "yes" && sd_openmp_cflags="${sd_openmp_cflags_def}"
        test "${sd_openmp_enable_cxx}" = "yes" && sd_openmp_cxxflags="${sd_openmp_cxxflags_def}"
        test "${sd_openmp_enable_fc}" = "yes" && sd_openmp_fcflags="${sd_openmp_fcflags_def}"
        sd_openmp_ldflags="${sd_openmp_ldflags_def}"
        sd_openmp_libs="${sd_openmp_libs_def}"
        ;;

      env)
        sd_openmp_cppflags="${sd_openmp_cppflags_def}"
        test "${sd_openmp_enable_cc}" = "yes" && sd_openmp_cflags="${sd_openmp_cflags_def}"
        test "${sd_openmp_enable_cxx}" = "yes" && sd_openmp_cxxflags="${sd_openmp_cxxflags_def}"
        test "${sd_openmp_enable_fc}" = "yes" && sd_openmp_fcflags="${sd_openmp_fcflags_def}"
        sd_openmp_ldflags="${sd_openmp_ldflags_def}"
        sd_openmp_libs="${sd_openmp_libs_def}"
        test "${OPENMP_CPPFLAGS}" != "" && sd_openmp_cppflags="${OPENMP_CPPFLAGS}"
        if test "${sd_openmp_enable_cc}" = "yes"; then
          test "${OPENMP_CFLAGS}" != "" && sd_openmp_cflags="${OPENMP_CFLAGS}"
        fi
        if test "${sd_openmp_enable_cxx}" = "yes"; then
          test "${OPENMP_CXXFLAGS}" != "" && sd_openmp_cxxflags="${OPENMP_CXXFLAGS}"
        fi
        if test "${sd_openmp_enable_fc}" = "yes"; then
          test "${OPENMP_FFLAGS}" != "" && sd_openmp_fcflags="${OPENMP_FFLAGS}"
          test "${OPENMP_FCFLAGS}" != "" && sd_openmp_fcflags="${OPENMP_FCFLAGS}"
        fi
        test "${OPENMP_LDFLAGS}" != "" && sd_openmp_ldflags="${OPENMP_LDFLAGS}"
        test "${OPENMP_LIBS}" != "" && sd_openmp_libs="${OPENMP_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for openMP: '${sd_openmp_init}'])
        ;;

    esac

  fi

  # Display configuration
  #if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
  #  _SD_OPENMP_DUMP_CONFIG
  #fi

  # Export configuration
  AC_SUBST(sd_openmp_options)
  AC_SUBST(sd_openmp_enable_def)
  AC_SUBST(sd_openmp_enable_cc)
  AC_SUBST(sd_openmp_enable_cxx)
  AC_SUBST(sd_openmp_enable_fc)
  AC_SUBST(sd_openmp_policy)
  AC_SUBST(sd_openmp_status)
  AC_SUBST(sd_openmp_enable)
  AC_SUBST(sd_openmp_init)
  AC_SUBST(sd_openmp_ok)
  AC_SUBST(sd_openmp_cc_ok)
  AC_SUBST(sd_openmp_cxx_ok)
  AC_SUBST(sd_openmp_fc_ok)
  AC_SUBST(sd_openmp_cppflags)
  AC_SUBST(sd_openmp_cflags)
  AC_SUBST(sd_openmp_cxxflags)
  AC_SUBST(sd_openmp_fcflags)
  AC_SUBST(sd_openmp_ldflags)
  AC_SUBST(sd_openmp_libs)
  AC_SUBST(enable_openmp)

  # Clean-up
  unset tmp_openmp_vars
]) # SD_OPENMP_INIT


AC_DEFUN([SD_OPENMP_DETECT], [
  # Display configuration
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
  _SD_OPENMP_DUMP_CONFIG
  fi

  # Check whether we can compile and run a simple programs
  # and update build flags if successful
  if test "${sd_openmp_enable}" = "yes" -o "${sd_openmp_enable}" = "auto"; then

    sd_openmp_ok="yes"
    if test "${sd_openmp_enable_cc}" = "yes"; then
      _SD_OPENMP_CHECK_CC
      if test "${sd_openmp_cc_ok}" = "no"; then
        sd_openmp_ok="no"
      fi
    fi
    if test "${sd_openmp_enable_cxx}" = "yes"; then
      _SD_OPENMP_CHECK_CXX
      if test "${sd_openmp_cxx_ok}" = "no"; then
        sd_openmp_ok="no"
      fi
    fi
    if test "${sd_openmp_enable_fc}" = "yes"; then
      _SD_OPENMP_CHECK_FC
      if test "${sd_openmp_fc_ok}" = "no"; then
        sd_openmp_ok="no"
      fi
    fi

    # Take decision according to policy
    if test "${sd_openmp_ok}" = "yes"; then
      sd_openmp_enable="yes"
    else
      if test "${sd_openmp_enable}" = "yes"; then
        case "${sd_openmp_policy}" in
          fail)
            AC_MSG_FAILURE([openMP support does not work])
            ;;
          skip)
            sd_openmp_enable="no"
            ;;
          warn)
            sd_openmp_enable="no"
            AC_MSG_WARN([openMP support does not work and has been disabled])
	    ;;
        esac
      else
        sd_openmp_enable="no"
      fi
    fi

    # Update build flags
    if test "${sd_openmp_ok}" = "yes"; then
      CPPFLAGS="${CPPFLAGS} ${sd_openmp_cppflags}"
      if test "${sd_openmp_enable_cc}" = "yes"; then
        CFLAGS="${CFLAGS} ${sd_openmp_cflags}"
      fi
      if test "${sd_openmp_enable_cxx}" = "yes"; then
        CXXFLAGS="${CXXFLAGS} ${sd_openmp_cxxflags}"
      fi
      if test "${sd_openmp_enable_fc}" = "yes"; then
        FCFLAGS="${FCFLAGS} ${sd_openmp_fcflags}"
      fi
      LIBS="${sd_openmp_libs} ${LIBS}"
      LDFLAGS="${LDFLAGS} ${sd_openmp_ldflags}"
    else
      if test "${sd_openmp_status}" = "optional" -a \
              "${sd_openmp_init}" = "def"; then
        sd_openmp_enable="no"
        sd_openmp_cppflags=""
        sd_openmp_cflags=""
        sd_openmp_fcflags=""
        sd_openmp_ldflags=""
        sd_openmp_libs=""
      else
        AC_MSG_WARN([invalid openMP configuration])
      fi
    fi

  else
    sd_openmp_enable="no"
    sd_openmp_cppflags=""
    sd_openmp_cflags=""
    sd_openmp_fcflags=""
    sd_openmp_ldflags=""
    sd_openmp_libs=""
  fi

  # Make openMP status available to the source code
  if test "${sd_openmp_enable}" = "yes" -a "${sd_openmp_ok}" = "yes"; then
    AC_DEFINE([HAVE_OPENMP], 1,
      [Define to 1 if you have a working openMP installation.])
  fi

]) # _SD_OPENMP_DETECT


                    # ------------------------------------ #
                    # ------------------------------------ #

#
# Private macros for C
#

AC_DEFUN([_SD_OPENMP_CHECK_CC], [
  # Init
  sd_openmp_cc_ok="unknown"
  sd_openmp_cc_api_ok="unknown"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_openmp_cppflags}"
  CFLAGS="${CFLAGS} ${sd_openmp_cflags}"
  LDFLAGS="${LDFLAGS} ${sd_openmp_ldflags}"
  LIBS="${sd_openmp_libs} ${LIBS}"

  # Check the openMP implementation
  _SD_OPENMP_CHECK_CC_API

  # Validate C support
  AC_MSG_CHECKING([whether the openMP C environment works])
  if test "${sd_openmp_cc_api_ok}" = "yes"; then
    sd_openmp_cc_ok="yes"
  else
    sd_openmp_cc_ok="no"
  fi
  AC_MSG_RESULT([${sd_openmp_cc_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_OPENMP_CHECK_CC

AC_DEFUN([_SD_OPENMP_CHECK_CC_API], [
  # Check if openMP directives work with C
  AC_MSG_CHECKING([whether the openMP directives work with C])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <omp.h>
#include <math.h>
    ]],
    [[
    int i, n, m;
    int array[10];
    n = 10;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        array[i] = sqrt((float)i);
    }
    m = omp_get_max_threads();
    return 0;
    ]])], [sd_openmp_cc_api_ok="yes"], [sd_openmp_cc_api_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_openmp_cc_api_ok}])
])

                    # ------------------------------------ #

#
# Private macros for C++
#

AC_DEFUN([_SD_OPENMP_CHECK_CXX], [
  # Init
  sd_openmp_cxx_ok="unknown"
  sd_openmp_cxx_api_ok="unknown"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_openmp_cppflags}"
  CXXFLAGS="${CXXFLAGS} ${sd_openmp_cxxflags}"
  LDFLAGS="${LDFLAGS} ${sd_openmp_ldflags}"
  LIBS="${sd_openmp_libs} ${LIBS}"

  # Check the openMP implementation
  _SD_OPENMP_CHECK_CXX_API

  # Validate C++ support
  AC_MSG_CHECKING([whether the openMP C++ environment works])
  if test "${sd_openmp_cxx_api_ok}" = "yes"; then
    sd_openmp_cxx_ok="yes"
  else
    sd_openmp_cxx_ok="no"
  fi
  AC_MSG_RESULT([${sd_openmp_cxx_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_OPENMP_CHECK_CXX

AC_DEFUN([_SD_OPENMP_CHECK_CXX_API], [
  # Check if openMP directives work with C++
  AC_MSG_CHECKING([whether the openMP directives work with C++])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <omp.h>
#include <cmath>
    ]],
    [[
    int i, n, m;
    int array[10];
    n = 10;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
        array[i] = std::sqrt(static_cast<float>(i));
    }
    m = omp_get_max_threads();
    return 0;
    ]])], [sd_openmp_cxx_api_ok="yes"], [sd_openmp_cxx_api_ok="no"])
  AC_LANG_POP([C++])
  AC_MSG_RESULT([${sd_openmp_cxx_api_ok}])
])

                    # ------------------------------------ #

#
# Private macros for Fortran 
#

AC_DEFUN([_SD_OPENMP_CHECK_FC], [
  # Init
  sd_openmp_fc_ok="unknown"
  sd_openmp_fc_api_ok="unknown"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_openmp_cppflags}"
  FCFLAGS="${FCFLAGS} ${sd_openmp_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_openmp_ldflags}"
  LIBS="${sd_openmp_libs} ${LIBS}"

  # Check the openMP implementation
  _SD_OPENMP_CHECK_FC_API

  # Validate Fortran support
  AC_MSG_CHECKING([whether the openMP Fortran environment works])
  if test "${sd_openmp_fc_api_ok}" = "yes"; then
    sd_openmp_fc_ok="yes"
  else
    sd_openmp_fc_ok="no"
  fi
  AC_MSG_RESULT([${sd_openmp_fc_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_OPENMP_CHECK_FC


AC_DEFUN([_SD_OPENMP_CHECK_FC_API], [
  # Check if openMP directives work with Fortran
  AC_MSG_CHECKING([whether the openMP directives work with Fortran])
  AC_LANG_PUSH([Fortran])
  AC_RUN_IFELSE([AC_LANG_PROGRAM([],
    [[
      use omp_lib 
      integer :: i,n,m
      integer :: array(10)
      n = 10 
!$OMP PARALLEL DO
      do i = 1, n
        array(i) = sqrt(real(i))
      enddo
!$OMP END PARALLEL DO
      m=omp_get_max_threads()
    ]])], [sd_openmp_fc_api_ok="yes"], [sd_openmp_fc_api_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_openmp_fc_api_ok}])
])


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_OPENMP_CHECK_CONFIG], [
  # Main trigger must be yes, no, or auto
  tmp_openmp_invalid="no"
  if test "${sd_openmp_enable}" != "auto" -a \
          "${sd_openmp_enable}" != "no" -a \
          "${sd_openmp_enable}" != "yes"; then
    case "${sd_openmp_policy}" in
      fail)
        AC_MSG_ERROR([invalid trigger value: sd_openmp_enable = '${sd_openmp_enable}'])
        ;;
      skip)
        tmp_openmp_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid trigger value: sd_openmp_enable = '${sd_openmp_enable}'])
        tmp_openmp_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_openmp_invalid}" = "yes"; then
    if test "${sd_openmp_status}" = "required"; then
      sd_openmp_enable="yes"
    else
      sd_openmp_enable="no"
    fi
    tmp_openmp_invalid="no"
    AC_MSG_NOTICE([setting sd_openmp_enable to '${sd_openmp_enable}'])
  fi

  # Check consistency between trigger value and package status
  tmp_openmp_invalid="no"
  if test "${sd_openmp_status}" = "implicit" -o \
          "${sd_openmp_status}" = "required"; then
    if test "${sd_openmp_enable}" = "no"; then
      case "${sd_openmp_policy}" in
        fail)
          AC_MSG_ERROR([The openMP option is required and cannot be disabled.])
          ;;
        skip)
          tmp_openmp_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The openMP option is recommended and cannot be disabled.])
          tmp_openmp_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_openmp_enable}" = "auto"; then
      AC_MSG_NOTICE([setting openMP trigger to yes])
      sd_openmp_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_openmp_invalid}" = "yes"; then
    case "${sd_openmp_status}" in
      implicit|required)
        sd_openmp_enable="yes"
        ;;
      optional)
        sd_openmp_enable="no"
        ;;
    esac
    tmp_openmp_invalid="no"
    AC_MSG_NOTICE([setting sd_openmp_enable to '${sd_openmp_enable}'])
  fi

  # When using environment variables, triggers must be set to yes
  if test "${sd_openmp_init}" = "env" -a "${sd_openmp_enable}" = "no"; then
    if test "${sd_openmp_policy}" != "skip"; then
      AC_MSG_WARN([openMP environment variables will be ignored])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_openmp_status}" = "implicit"; then
    sd_openmp_enable="yes"
    if test "${sd_openmp_cppflags}" != ""; then
      sd_openmp_cppflags=""
      AC_MSG_NOTICE([resetting openMP CPP flags (implicit package)])
    fi
    if test "${sd_openmp_cflags}" != ""; then
      sd_openmp_cflags=""
      AC_MSG_NOTICE([resetting openMP C flags (implicit package)])
    fi
    if test "${sd_openmp_cxxflags}" != ""; then
      sd_openmp_cxxflags=""
      AC_MSG_NOTICE([resetting openMP C++ flags (implicit package)])
    fi
    if test "${sd_openmp_fcflags}" != ""; then
      sd_openmp_fcflags=""
      AC_MSG_NOTICE([resetting openMP Fortran flags (implicit package)])
    fi
    if test "${sd_openmp_ldflags}" != ""; then
      sd_openmp_ldflags=""
      AC_MSG_NOTICE([resetting openMP linker flags (implicit package)])
    fi
    if test "${sd_openmp_libs}" != ""; then
      sd_openmp_libs=""
      AC_MSG_NOTICE([resetting openMP library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_openmp_enable}" = "no"; then
    sd_openmp_cppflags=""
    sd_openmp_cflags=""
    sd_openmp_cxxflags=""
    sd_openmp_fcflags=""
    sd_openmp_ldflags=""
    sd_openmp_libs=""
    sd_openmp_ok="no"
  fi

  # Clean-up
  unset tmp_openmp_invalid
]) # _SD_OPENMP_CHECK_CONFIG


AC_DEFUN([_SD_OPENMP_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable openMP])
  AC_MSG_RESULT([${sd_openmp_enable}])
  if test "${sd_openmp_enable}" != "no"; then
    AC_MSG_CHECKING([how openMP parameters have been set])
    AC_MSG_RESULT([${sd_openmp_init}])
    # AC_MSG_CHECKING([for openMP C preprocessing flags])
    # AC_MSG_RESULT([${sd_openmp_cppflags}])
    if test "${sd_openmp_enable_cc}" = "yes"; then
      AC_MSG_CHECKING([for openMP C flags])
      AC_MSG_RESULT([${sd_openmp_cflags}])
    fi
    if test "${sd_openmp_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([for openMP C++ flags])
      AC_MSG_RESULT([${sd_openmp_cxxflags}])
    fi
    if test "${sd_openmp_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([for openMP Fortran flags])
      AC_MSG_RESULT([${sd_openmp_fcflags}])
    fi
    AC_MSG_CHECKING([for openMP linker flags])
    AC_MSG_RESULT([${sd_openmp_ldflags}])
    AC_MSG_CHECKING([for openMP library flags])
    AC_MSG_RESULT([${sd_openmp_libs}])
  fi
]) # _SD_OPENMP_DUMP_CONFIG

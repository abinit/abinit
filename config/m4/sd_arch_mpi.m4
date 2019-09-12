## Copyright (C) 2019 Yann Pouillon

#
# MPI detection with Steredeg
#


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_MPI_INIT], [
  # Init
  sd_mpi_cppflags=""
  sd_mpi_cflags=""
  sd_mpi_fcflags=""
  sd_mpi_ldflags=""
  sd_mpi_libs=""
  sd_mpi_cc=""
  sd_mpi_cc_ok="unknown"
  sd_mpi_cc_set="no"
  sd_mpi_cxx=""
  sd_mpi_cxx_ok="unknown"
  sd_mpi_cxx_set="no"
  sd_mpi_fc=""
  sd_mpi_fc_ok="unknown"
  sd_mpi_fc_set="no"
  sd_mpi_enable=""
  sd_mpi_init="unknown"
  sd_mpi_ok="unknown"

  # Set adjustable parameters
  sd_mpi_options="$1"
  sd_mpi_libs_def="$2"
  sd_mpi_cppflags_def="$3"
  sd_mpi_cflags_def="$4"
  sd_mpi_cxxflags_def="$5"
  sd_mpi_fcflags_def="$6"
  sd_mpi_ldflags_def="$7"
  sd_mpi_enable_cxx=""
  sd_mpi_enable_def=""
  sd_mpi_enable_fc=""
  sd_mpi_policy=""
  sd_mpi_status=""

  # Process options
  for kwd in ${sd_mpi_options}; do
    case "${kwd}" in
      auto)
        sd_mpi_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_mpi_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_mpi_policy="${kwd}"
        ;;
      no-cxx)
        sd_mpi_enable_cxx="no"
        ;;
      no-fortran)
        sd_mpi_enable_fc="no"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg MPI option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_mpi_libs_def}" && sd_mpi_libs_def="-lmpi"
  test -z "${sd_mpi_policy}" && sd_mpi_policy="fail"
  test -z "${sd_mpi_status}" && sd_mpi_status="optional"
  test -z "${sd_mpi_enable_def}" && sd_mpi_enable_def="no"
  case "${sd_mpi_status}" in
    implicit|required)
      sd_mpi_enable_def="yes"
      ;;
  esac
  test -z "${sd_mpi_enable_cxx}" && sd_mpi_enable_cxx="yes"
  test -z "${sd_mpi_enable_fc}" && sd_mpi_enable_fc="yes"

  # Declare configure option
  AC_ARG_WITH([mpi],
    [AS_HELP_STRING(
      [--with-mpi],
      [Install prefix of MPI (e.g. /usr/local). The default behaviour is to detect whether the specified compilers properly support MPI and to fall back to serial mode if not. You may use --with-mpi without argument to force MPI detection, in which case detection failures will result in errors, and --without-mpi to disable MPI support completely.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_mpi_enable="${withval}"
        sd_mpi_init="yon"
      else
        sd_mpi_enable="yes"
        sd_mpi_init="dir"
      fi],
    [ sd_mpi_enable="${sd_mpi_enable_def}"; sd_mpi_init="def"])

  # Declare environment variables
  AC_ARG_VAR([MPI_CPPFLAGS], [C preprocessing flags for MPI.])
  AC_ARG_VAR([MPI_CFLAGS], [C flags for MPI.])
  AC_ARG_VAR([MPI_CXXFLAGS], [C++ flags for MPI.])
  AC_ARG_VAR([MPI_FCFLAGS], [Fortran flags for MPI.])
  AC_ARG_VAR([MPI_LDFLAGS], [Linker flags for MPI.])
  AC_ARG_VAR([MPI_LIBS], [Library flags for MPI.])

  # Detect use of environment variables
  if test "${sd_mpi_enable}" = "yes" -o "${sd_mpi_enable}" = "auto"; then
    tmp_compil_vars="${CC}"
    tmp_mpi_vars="${MPI_CPPFLAGS}${MPI_CFLAGS}${MPI_LDFLAGS}${MPI_LIBS}"
    if test "${sd_mpi_enable_cxx}" = "yes"; then
      tmp_compil_vars="${tmp_compil_vars}${CXX}"
      tmp_mpi_vars="${tmp_mpi_vars}${MPI_CXXFLAGS}"
    fi
    if test "${sd_mpi_enable_fc}" = "yes"; then
      tmp_compil_vars="${tmp_compil_vars}${FC}"
      tmp_mpi_vars="${tmp_mpi_vars}${MPI_FCFLAGS}"
    fi
    if test "${sd_mpi_init}" = "def" -a ! -z "${tmp_mpi_vars}"; then
      sd_mpi_enable="yes"
      sd_mpi_init="env"
    fi
    if test "${sd_mpi_init}" = "def" -a ! -z "${tmp_compil_vars}"; then
      sd_mpi_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
    _SD_MPI_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_mpi_enable}" = "yes" -o "${sd_mpi_enable}" = "auto"; then

    # Set MPI-specific flags
    case "${sd_mpi_init}" in

      def|yon)
        sd_mpi_cppflags="${sd_mpi_cppflags_def}"
        sd_mpi_cflags="${sd_mpi_cflags_def}"
        test "${sd_mpi_enable_cxx}" = "yes" && \
          sd_mpi_cxxflags="${sd_mpi_cxxflags_def}"
        test "${sd_mpi_enable_fc}" = "yes" && \
          sd_mpi_fcflags="${sd_mpi_fcflags_def}"
        sd_mpi_ldflags="${sd_mpi_ldflags_def}"
        sd_mpi_libs="${sd_mpi_libs_def}"
        ;;

      dir)
        sd_mpi_cppflags="-I${with_mpi}/include"
        sd_mpi_cflags="${sd_mpi_cflags_def}"
        test "${sd_mpi_enable_cxx}" = "yes" && \
          sd_mpi_cxxflags="${sd_mpi_cxxflags_def}"
        test "${sd_mpi_enable_fc}" = "yes" && \
          sd_mpi_fcflags="${sd_mpi_fcflags_def} -I${with_mpi}/include"
        sd_mpi_ldflags="${sd_mpi_ldflags_def}"
        sd_mpi_libs="-L${with_mpi}/lib ${sd_mpi_libs_def}"
        ;;

      env)
        sd_mpi_cppflags="${sd_mpi_cppflags_def}"
        sd_mpi_cflags="${sd_mpi_cflags_def}"
        test "${sd_mpi_enable_cxx}" = "yes" && \
          sd_mpi_cxxflags="${sd_mpi_cxxflags_def}"
        test "${sd_mpi_enable_fc}" = "yes" && \
          sd_mpi_fcflags="${sd_mpi_fcflags_def}"
        sd_mpi_ldflags="${sd_mpi_ldflags_def}"
        sd_mpi_libs="${sd_mpi_libs_def}"
        test ! -z "${MPI_CPPFLAGS}" && sd_mpi_cppflags="${MPI_CPPFLAGS}"
        test ! -z "${MPI_CFLAGS}" && sd_mpi_cflags="${MPI_CFLAGS}"
        if test "${sd_mpi_enable_cxx}" = "yes"; then
          test ! -z "${MPI_CXXFLAGS}" && sd_mpi_cxxflags="${MPI_CXXFLAGS}"
        fi
        if test "${sd_mpi_enable_fc}" = "yes"; then
          test ! -z "${MPI_FCFLAGS}" && sd_mpi_fcflags="${MPI_FCFLAGS}"
        fi
        test ! -z "${MPI_LDFLAGS}" && sd_mpi_ldflags="${MPI_LDFLAGS}"
        test ! -z "${MPI_LIBS}" && sd_mpi_libs="${MPI_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for MPI: '${sd_mpi_init}'])
        ;;

    esac

    # Look for MPI compilers
    if test "${sd_mpi_init}" != "env"; then
      tmp_mpi_keep_libs=""
      _SD_MPI_INIT_CC
      tmp_mpi_keep_libs="${tmp_mpi_keep_libs}${sd_mpi_cc_set}"
      if test "${sd_mpi_enable_cxx}" = "yes"; then
        _SD_MPI_INIT_CXX
        tmp_mpi_keep_libs="${tmp_mpi_keep_libs}${sd_mpi_cxx_set}"
      fi
      if test "${sd_mpi_enable_fc}" = "yes"; then
        _SD_MPI_INIT_FC
        tmp_mpi_keep_libs="${tmp_mpi_keep_libs}${sd_mpi_fc_set}"
      fi
      tmp_mpi_keep_libs=`echo "${tmp_mpi_keep_libs}" | grep "no"`
      if test "${tmp_mpi_keep_libs}" = ""; then
        sd_mpi_libs=""
      fi
      unset tmp_mpi_keep_libs
    fi

  fi

  # Display configuration
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
    _SD_MPI_DUMP_CONFIG
  fi

  # Export configuration
  AC_SUBST(sd_mpi_options)
  AC_SUBST(sd_mpi_enable_def)
  AC_SUBST(sd_mpi_enable_cxx)
  AC_SUBST(sd_mpi_enable_fc)
  AC_SUBST(sd_mpi_policy)
  AC_SUBST(sd_mpi_status)
  AC_SUBST(sd_mpi_enable)
  AC_SUBST(sd_mpi_init)
  AC_SUBST(sd_mpi_ok)
  AC_SUBST(sd_mpi_cppflags)
  AC_SUBST(sd_mpi_cflags)
  AC_SUBST(sd_mpi_fcflags)
  AC_SUBST(sd_mpi_ldflags)
  AC_SUBST(sd_mpi_libs)
  AC_SUBST(with_mpi)

  # Clean-up
  unset tmp_compil_vars
  unset tmp_mpi_vars
]) # SD_MPI_INIT


AC_DEFUN([SD_MPI_DETECT], [
  # Display configuration
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
  _SD_MPI_DUMP_CONFIG
  fi

  # Check compilers or APIs
  if test "${sd_mpi_enable}" = "yes" -o "${sd_mpi_enable}" = "auto"; then

    _SD_MPI_CHECK_CC
    if test "${sd_mpi_cc_ok}" = "yes"; then
      if test "${sd_mpi_enable_cxx}" = "yes"; then
        _SD_MPI_CHECK_CXX
      fi
      if test "${sd_mpi_enable_fc}" = "yes"; then
        _SD_MPI_CHECK_FC
      fi
    fi

    # Validate implementation status
    tmp_mpi_ok="yes"
    if test "${sd_mpi_cc_ok}" = "yes"; then
      if test "${sd_mpi_enable_cxx}" = "yes" -a "${sd_mpi_cxx_ok}" != "yes"; then
        tmp_mpi_ok="no"
      fi
      if test "${sd_mpi_enable_fc}" = "yes" -a "${sd_mpi_fc_ok}" != "yes"; then
        tmp_mpi_ok="no"
      fi
      sd_mpi_ok="${tmp_mpi_ok}"
    else
      sd_mpi_ok="no"
    fi

    # Take decision according to policy
    if test "${sd_mpi_ok}" = "yes"; then
      sd_mpi_enable="yes"
    else
      if test "${sd_mpi_enable}" = "yes"; then
        case "${sd_mpi_policy}" in
          fail)
            AC_MSG_FAILURE([MPI support does not work])
            ;;
          skip)
            sd_mpi_enable="no"
            ;;
          warn)
            sd_mpi_enable="no"
            AC_MSG_WARN([MPI support does not work and has been disabled])
	    ;;
        esac
      else
        sd_mpi_enable="no"
      fi
    fi

  fi

  # Make MPI status available to the source code
  if test "${sd_mpi_enable}" = "yes" -a "${sd_mpi_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI], 1,
      [Define to 1 if you have a working MPI installation.])
  fi
]) # SD_MPI_DETECT


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Private macros for C
#


AC_DEFUN([_SD_MPI_INIT_CC], [
  dnl Look for a MPI C compiler
  case "${sd_mpi_init}" in

    dir)
      if test -z "${CC}"; then
        sd_mpi_cc="${with_mpi}/bin/mpicc"
      else
        tmp_cc_has_path=`echo "${CC}" | grep '/'`
        if test -z "${tmp_cc_has_path}"; then
          sd_mpi_cc="${with_mpi}/bin/${CC}"
          test -x "${sd_mpi_cc}" || sd_mpi_cc="${with_mpi}/bin/mpicc"
        else
          sd_mpi_libs="-L${with_mpi}/lib ${sd_mpi_libs_def}"
          AC_MSG_NOTICE([user-defined MPI library flags: ${sd_mpi_libs}])
        fi
      fi

      if test ! -z "${sd_mpi_cc}"; then
        AC_MSG_NOTICE([user-defined MPI C compiler: ${sd_mpi_cc}])
        AC_MSG_CHECKING([for an executable MPI C compiler])
        if test -x "${sd_mpi_cc}"; then
          AC_MSG_RESULT([yes])
          AC_MSG_NOTICE([setting CC to ${sd_mpi_cc}])
          CC="${sd_mpi_cc}"
          sd_mpi_cc_set="yes"
        else
          AC_MSG_RESULT([not found])
          AC_MSG_ERROR([invalid MPI settings
                  Please adjust --with-mpi and/or CC and re-run configure])
        fi
      fi
      ;;

    def|env|yon)
      sd_mpi_cc="${CC}"
      if test ! -z "${sd_mpi_cc}"; then
        sd_mpi_cc_set="yes"
      else
        AC_CHECK_PROGS([sd_mpi_cc], [mpicc])
        if test ! -z "${sd_mpi_cc}"; then
          AC_MSG_NOTICE([setting CC to '${sd_mpi_cc}'])
          CC="${sd_mpi_cc}"
          sd_mpi_cc_set="yes"
        fi
      fi
      ;;

  esac
]) # _SD_MPI_INIT_CC


AC_DEFUN([_SD_MPI_CHECK_CC], [
  # Init
  sd_mpi_cc_ok="unknown"
  sd_mpi_cc_api_ok="unknown"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_mpi_cppflags}"
  CFLAGS="${CFLAGS} ${sd_mpi_cflags}"
  LDFLAGS="${LDFLAGS} ${sd_mpi_ldflags}"
  LIBS="${sd_mpi_libs} ${LIBS}"

  # Check the MPI implementation
  _SD_MPI_CHECK_CC_API

  # Validate C support
  AC_MSG_CHECKING([whether the MPI C environment works])
  if test "${sd_mpi_cc_api_ok}" = "yes"; then
    sd_mpi_cc_ok="yes"
  else
    sd_mpi_cc_ok="no"
  fi
  AC_MSG_RESULT([${sd_mpi_cc_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_MPI_CHECK_CC


AC_DEFUN([_SD_MPI_CHECK_CC_API], [
  # Check MPI C API
  AC_MSG_CHECKING([whether the MPI C API works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <mpi.h>
    ]],
    [[
      int *argc;
      char **argv;
      MPI_Init(argc, argv);
    ]])], [sd_mpi_cc_api_ok="yes"], [sd_mpi_cc_api_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_mpi_cc_api_ok}])
])


                    # ------------------------------------ #


#
# Private macros for C++
#


AC_DEFUN([_SD_MPI_INIT_CXX], [
  # Look for a MPI C compiler
  case "${sd_mpi_init}" in

    dir)
      if test -z "${CXX}"; then
        sd_mpi_cxx="${with_mpi}/bin/mpic++"
      else
        tmp_cxx_has_path=`echo "${CXX}" | grep '/'`
        if test -z "${tmp_cxx_has_path}"; then
          sd_mpi_cxx="${with_mpi}/bin/${CXX}"
          test -x "${sd_mpi_cxx}" || sd_mpi_cxx="${with_mpi}/bin/mpic++"
        else
          sd_mpi_libs="-L${with_mpi}/lib ${sd_mpi_libs_def}"
          AC_MSG_NOTICE([user-defined MPI library flags: ${sd_mpi_libs}])
        fi
      fi

      if test ! -z "${sd_mpi_cxx}"; then
        AC_MSG_NOTICE([user-defined MPI C++ compiler: ${sd_mpi_cxx}])
        AC_MSG_CHECKING([for an executable MPI C compiler])
        if test -x "${sd_mpi_cxx}"; then
          AC_MSG_RESULT([yes])
          AC_MSG_NOTICE([setting CXX to ${sd_mpi_cxx}])
          CXX="${sd_mpi_cxx}"
          sd_mpi_cxx_set="yes"
        else
          AC_MSG_RESULT([not found])
          AC_MSG_ERROR([invalid MPI settings
                  Please adjust --with-mpi and/or CXX and re-run configure])
        fi
      fi
      ;;

    def|env|yon)
      sd_mpi_cxx="${CXX}"
      if test ! -z "${sd_mpi_cxx}"; then
        sd_mpi_cxx_set="yes"
      else
        AC_CHECK_PROGS([sd_mpi_cxx], [mpic++ mpicxx])
        if test ! -z "${sd_mpi_cxx}"; then
          AC_MSG_NOTICE([setting CXX to '${sd_mpi_cxx}'])
          CXX="${sd_mpi_cxx}"
          sd_mpi_cxx_set="yes"
        fi
      fi
      ;;

  esac
]) # _SD_MPI_INIT_CXX


AC_DEFUN([_SD_MPI_CHECK_CXX], [
  # Init
  sd_mpi_cxx_ok="unknown"
  sd_mpi_cxx_api_ok="unknown"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_mpi_cppflags}"
  CXXFLAGS="${CXXFLAGS} ${sd_mpi_cxxflags}"
  LDFLAGS="${LDFLAGS} ${sd_mpi_ldflags}"
  LIBS="${sd_mpi_libs} ${LIBS}"

  # Check the MPI implementation
  _SD_MPI_CHECK_CXX_API

  # Validate C++ support
  AC_MSG_CHECKING([whether the MPI C++ environment works])
  if test "${sd_mpi_cxx_api_ok}" = "yes"; then
    sd_mpi_cxx_ok="yes"
  else
    sd_mpi_cxx_ok="no"
  fi
  AC_MSG_RESULT([${sd_mpi_cxx_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_MPI_CHECK_CXX


AC_DEFUN([_SD_MPI_CHECK_CXX_API], [
  # Check MPI C++ API
  AC_MSG_CHECKING([whether the MPI C++ API works])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <mpi.h>
    ]],
    [[
      MPI::Init()
    ]])], [sd_mpi_cxx_api_ok="yes"], [sd_mpi_cxx_api_ok="no"])
  AC_LANG_POP([C++])
  AC_MSG_RESULT([${sd_mpi_cxx_api_ok}])
])


                    # ------------------------------------ #


#
# Private macros for Fortran
#


AC_DEFUN([_SD_MPI_INIT_FC], [
  # Look for a MPI Fortran compiler
  case "${sd_mpi_init}" in

    dir)
      if test -z "${FC}"; then
        sd_mpi_fc="${with_mpi}/bin/mpif90"
      else
        tmp_fc_has_path=`echo "${FC}" | grep '/'`
        if test -z "${tmp_fc_has_path}"; then
          sd_mpi_fc="${with_mpi}/bin/${FC}"
          test -x "${sd_mpi_fc}" || sd_mpi_fc="${with_mpi}/bin/mpif90"
        else
          sd_mpi_libs="-L${with_mpi}/lib ${sd_mpi_libs_def}"
          AC_MSG_NOTICE([user-defined MPI library flags: ${sd_mpi_libs}])
        fi
      fi

      if test ! -z "${sd_mpi_fc}"; then
        AC_MSG_NOTICE([user-defined MPI Fortran compiler: ${sd_mpi_fc}])
        AC_MSG_CHECKING([for an executable MPI Fortran compiler])
        if test -x "${sd_mpi_fc}"; then
          AC_MSG_RESULT([yes])
          AC_MSG_NOTICE([setting FC to ${sd_mpi_fc}])
          FC="${sd_mpi_fc}"
          sd_mpi_fc_set="yes"
        else
          AC_MSG_RESULT([not found])
          AC_MSG_ERROR([invalid MPI settings for Fortran
                  Please adjust --with-mpi and/or FC and re-run configure])
        fi
      fi
      ;;

    def|env|yon)
      sd_mpi_fc="${FC}"
      if test ! -z "${sd_mpi_fc}"; then
        sd_mpi_fc_set="yes"
      else
        AC_CHECK_PROGS([sd_mpi_fc], [mpifort mpif90 mpif95])
        if test ! -z "${sd_mpi_fc}"; then
          AC_MSG_NOTICE([setting FC to '${sd_mpi_fc}'])
          FC="${sd_mpi_fc}"
          sd_mpi_fc_set="yes"
        fi
      fi
      ;;

  esac
]) # _SD_MPI_INIT_FC


AC_DEFUN([_SD_MPI_CHECK_FC], [
  # Init
  sd_mpi_fc_ok="unknown"
  sd_mpi_fc_api_ok="unknown"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_mpi_cppflags}"
  FCFLAGS="${FCFLAGS} ${sd_mpi_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_mpi_ldflags}"
  LIBS="${sd_mpi_libs} ${LIBS}"

  # Check the MPI implementation
  _SD_MPI_CHECK_FC_API

  # Validate Fortran support
  AC_MSG_CHECKING([whether the MPI Fortran environment works])
  if test "${sd_mpi_fc_api_ok}" = "yes"; then
    sd_mpi_fc_ok="yes"
  else
    sd_mpi_fc_ok="no"
  fi
  AC_MSG_RESULT([${sd_mpi_fc_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_MPI_CHECK_FC


AC_DEFUN([_SD_MPI_CHECK_FC_API], [
  # Check MPI C API
  AC_MSG_CHECKING([whether the MPI Fortran API works])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use mpi
      integer :: ierr
      call mpi_init(ierr)
    ]])], [sd_mpi_fc_api_ok="yes"], [sd_mpi_fc_api_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_mpi_fc_api_ok}])
])


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_MPI_CHECK_CONFIG], [
  # Main trigger must be yes, no, or auto
  tmp_mpi_invalid="no"
  if test "${sd_mpi_enable}" != "auto" -a \
          "${sd_mpi_enable}" != "no" -a \
          "${sd_mpi_enable}" != "yes"; then
    case "${sd_mpi_policy}" in
      fail)
        AC_MSG_ERROR([invalid trigger value: sd_mpi_enable = '${sd_mpi_enable}'])
        ;;
      skip)
        tmp_mpi_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid trigger value: sd_mpi_enable = '${sd_mpi_enable}'])
        tmp_mpi_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_mpi_invalid}" = "yes"; then
    if test "${sd_mpi_status}" = "required"; then
      sd_mpi_enable="yes"
    else
      sd_mpi_enable="no"
    fi
    tmp_mpi_invalid="no"
    AC_MSG_NOTICE([setting sd_mpi_enable to '${sd_mpi_enable}'])
  fi

  # Check consistency between trigger value and package status
  tmp_mpi_invalid="no"
  if test "${sd_mpi_status}" = "implicit" -o \
          "${sd_mpi_status}" = "required"; then
    if test "${sd_mpi_enable}" = "no"; then
      case "${sd_mpi_policy}" in
        fail)
          AC_MSG_ERROR([The MPI package is required and cannot be disabled.])
          ;;
        skip)
          tmp_mpi_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The MPI package is required and cannot be disabled.])
          tmp_mpi_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_mpi_enable}" = "auto"; then
      AC_MSG_NOTICE([setting MPI trigger to yes])
      sd_mpi_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_mpi_invalid}" = "yes"; then
    case "${sd_mpi_status}" in
      implicit|required)
        sd_mpi_enable="yes"
        ;;
      optional)
        sd_mpi_enable="no"
        ;;
    esac
    tmp_mpi_invalid="no"
    AC_MSG_NOTICE([setting sd_mpi_enable to '${sd_mpi_enable}'])
  fi

  # When using environment variables, triggers must be set to yes
  if test "${sd_mpi_init}" = "env" -a "${sd_mpi_enable}" = "no"; then
    if test "${sd_mpi_policy}" != "skip"; then
      AC_MSG_WARN([MPI environment variables will be ignored])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_mpi_status}" = "implicit"; then
    sd_mpi_enable="yes"
    if test "${sd_mpi_cppflags}" != ""; then
      sd_mpi_cppflags=""
      AC_MSG_NOTICE([resetting MPI CPP flags (implicit package)])
    fi
    if test "${sd_mpi_cflags}" != ""; then
      sd_mpi_cflags=""
      AC_MSG_NOTICE([resetting MPI C flags (implicit package)])
    fi
    if test "${sd_mpi_cxxflags}" != ""; then
      sd_mpi_cxxflags=""
      AC_MSG_NOTICE([resetting MPI C++ flags (implicit package)])
    fi
    if test "${sd_mpi_fcflags}" != ""; then
      sd_mpi_fcflags=""
      AC_MSG_NOTICE([resetting MPI Fortran flags (implicit package)])
    fi
    if test "${sd_mpi_ldflags}" != ""; then
      sd_mpi_ldflags=""
      AC_MSG_NOTICE([resetting MPI linker flags (implicit package)])
    fi
    if test "${sd_mpi_libs}" != ""; then
      sd_mpi_libs=""
      AC_MSG_NOTICE([resetting MPI library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_mpi_enable}" = "no"; then
    sd_mpi_cppflags=""
    sd_mpi_cflags=""
    sd_mpi_cxxflags=""
    sd_mpi_fcflags=""
    sd_mpi_ldflags=""
    sd_mpi_libs=""
    sd_mpi_ok="no"
  fi

  # Clean-up
  unset tmp_mpi_invalid
]) # _SD_MPI_CHECK_CONFIG


AC_DEFUN([_SD_MPI_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable MPI])
  AC_MSG_RESULT([${sd_mpi_enable}])
  if test "${sd_mpi_enable}" != "no"; then
    AC_MSG_CHECKING([how MPI parameters have been set])
    AC_MSG_RESULT([${sd_mpi_init}])
    AC_MSG_CHECKING([whether the MPI C compiler is set])
    AC_MSG_RESULT([${sd_mpi_cc_set}])
    if test "${sd_mpi_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([whether the MPI C++ compiler is set])
      AC_MSG_RESULT([${sd_mpi_cxx_set}])
    fi
    if test "${sd_mpi_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([whether the MPI Fortran compiler is set])
      AC_MSG_RESULT([${sd_mpi_fc_set}])
    fi
    AC_MSG_CHECKING([for MPI C preprocessing flags])
    AC_MSG_RESULT([${sd_mpi_cppflags}])
    AC_MSG_CHECKING([for MPI C flags])
    AC_MSG_RESULT([${sd_mpi_cflags}])
    if test "${sd_mpi_enable_cxx}" = "yes"; then
      AC_MSG_CHECKING([for MPI C++ flags])
      AC_MSG_RESULT([${sd_mpi_cxxflags}])
    fi
    if test "${sd_mpi_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([for MPI Fortran flags])
      AC_MSG_RESULT([${sd_mpi_fcflags}])
    fi
    AC_MSG_CHECKING([for MPI linker flags])
    AC_MSG_RESULT([${sd_mpi_ldflags}])
    AC_MSG_CHECKING([for MPI library flags])
    AC_MSG_RESULT([${sd_mpi_libs}])
  fi
]) # _SD_MPI_DUMP_CONFIG

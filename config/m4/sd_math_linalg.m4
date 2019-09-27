# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 Yann Pouillon, Marc Torrent
#
# This file is part of the Steredeg software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Support for external linear algebra libraries
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_LINALG_INIT], [
  # FIXME: Importing parameters from ABINIT
  sd_fc_vendor="${abi_fc_vendor}"
  sd_fc_version="${abi_fc_version}"

  # Init build flags
  sd_linalg_cppflags=""
  sd_linalg_cflags=""
  sd_linalg_cxxflags=""
  sd_linalg_fcflags=""
  sd_linalg_ldflags=""
  sd_linalg_libdir=""
  sd_linalg_libs=""
  sd_linalg_enable=""
  sd_linalg_init="unknown"
  sd_linalg_ok="unknown"

  # Init internal variables
  sd_linalg_chk_gpu=""
  sd_linalg_chk_mpi=""
  sd_linalg_chk_mpiacc=""
  sd_linalg_chk_serial=""
  sd_linalg_has_blas="unknown"
  sd_linalg_has_lapack="unknown"
  sd_linalg_has_lapacke="unknown"
  sd_linalg_has_blacs="unknown"
  sd_linalg_has_scalapack="unknown"
  sd_linalg_has_elpa="unknown"
  sd_linalg_has_elpa_2013="unknown"
  sd_linalg_has_elpa_2014="unknown"
  sd_linalg_has_elpa_2015="unknown"
  sd_linalg_has_elpa_2016="unknown"
  sd_linalg_has_plasma="unknown"
  sd_linalg_has_magma="unknown"
  sd_linalg_gpu_ok="unknown"
  sd_linalg_mpi_ok="unknown"
  sd_linalg_mpiacc_ok="unknown"
  sd_linalg_serial_ok="unknown"
  sd_linalg_provided=""

  # Set adjustable parameters
  sd_linalg_options="$1"
  sd_linalg_libs_def="$2"
  sd_linalg_cppflags_def="$3"
  sd_linalg_cflags_def="$4"
  sd_linalg_cxxflags_def="$5"
  sd_linalg_fcflags_def="$6"
  sd_linalg_ldflags_def="$7"
  sd_linalg_enable_def=""
  sd_linalg_policy=""
  sd_linalg_status=""

  # Process options
  for kwd in ${sd_linalg_options}; do
    case "${kwd}" in
      auto)
        sd_linalg_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_linalg_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_linalg_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg linear algebra option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_linalg_libs}" && sd_linalg_libs="-llapack -lblas"
  test -z "${sd_linalg_policy}" && sd_linalg_policy="fail"
  test -z "${sd_linalg_status}" && sd_linalg_status="optional"
  test -z "${sd_linalg_enable_def}" && sd_linalg_enable_def="no"
  case "${sd_linalg_status}" in
    implicit|required)
      sd_linalg_enable_def="yes"
      ;;
  esac

  # Declare configure option
  AC_ARG_WITH([linalg],
    [AS_HELP_STRING([--with-linalg],
      [Install prefix of the linear algebra libraries (e.g. /usr/local), or automatic configuration if set to 'yes' or 'no'.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_linalg_enable="${withval}"
        sd_linalg_init="yon"
      else
        sd_linalg_enable="yes"
        sd_linalg_init="dir"
        test -d "${withval}/lib" && sd_linalg_libdir="${withval}/lib"
        test -d "${withval}/lib64" && sd_linalg_libdir="${withval}/lib64"
      fi],
    [ sd_linalg_enable="${sd_linalg_enable_def}"; sd_linalg_init="def"])

  # Declare environment variables
  AC_ARG_VAR([LINALG_CPPFLAGS], [C preprocessing flags for linear algebra.])
  AC_ARG_VAR([LINALG_CFLAGS], [C flags for linear algebra.])
  AC_ARG_VAR([LINALG_CXXFLAGS], [C++ flags for linear algebra.])
  AC_ARG_VAR([LINALG_FCFLAGS], [Fortran flags for linear algebra.])
  AC_ARG_VAR([LINALG_LDFLAGS], [Linker flags for linear algebra.])
  AC_ARG_VAR([LINALG_LIBS], [Library flags for linear algebra.])

  # Detect use of environment variables
  if test "${sd_linalg_enable}" = "yes" -o "${sd_linalg_enable}" = "auto"; then
    tmp_linalg_vars="${LINALG_CPPFLAGS}${LINALG_CFLAGS}${LINALG_CXXFLAGS}${LINALG_FCFLAGS}${LINALG_LDFLAGS}${LINALG_LIBS}"
    if test "${sd_linalg_init}" = "def" -a ! -z "${tmp_linalg_vars}"; then
      sd_linalg_enable="yes"
      sd_linalg_init="env"
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_LINALG_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_linalg_enable}" = "yes" -o "${sd_linalg_enable}" = "auto"; then

    # Set LINALG-specific flags
    case "${sd_linalg_init}" in

      def|yon)
        if test "${sd_linalg_enable}" != "no"; then
          sd_linalg_cppflags="${sd_linalg_cppflags_def}"
          sd_linalg_cflags="${sd_linalg_cflags_def}"
          sd_linalg_cxxflags="${sd_linalg_cxxflags_def}"
          sd_linalg_fcflags="${sd_linalg_fcflags_def}"
          sd_linalg_ldflags="${sd_linalg_ldflags_def}"
          sd_linalg_libs="${sd_linalg_libs_def}"
        fi
        ;;

      dir)
        sd_linalg_cppflags="${sd_linalg_cppflags_def} -I${with_linalg}/include"
        sd_linalg_cflags="${sd_linalg_cflags_def}"
        sd_linalg_cxxflags="${sd_linalg_cxxflags_def}"
        sd_linalg_fcflags="${sd_linalg_fcflags_def} -I${with_linalg}/include"
        sd_linalg_ldflags="${sd_linalg_ldflags_def}"
        sd_linalg_libs="-L${sd_linalg_libdir} ${sd_linalg_libs_def}"
        ;;

      env)
        sd_linalg_cppflags="${sd_linalg_cppflags_def}"
        sd_linalg_cflags="${sd_linalg_cflags_def}"
        sd_linalg_cxxflags="${sd_linalg_cxxflags_def}"
        sd_linalg_fcflags="${sd_linalg_fcflags_def}"
        sd_linalg_ldflags="${sd_linalg_ldflags_def}"
        sd_linalg_libs="${sd_linalg_libs_def}"
        test ! -z "${LINALG_CPPFLAGS}" && sd_linalg_cppflags="${LINALG_CPPFLAGS}"
        test ! -z "${LINALG_CFLAGS}" && sd_linalg_cflags="${LINALG_CFLAGS}"
        test ! -z "${LINALG_CXXFLAGS}" && sd_linalg_cxxflags="${LINALG_CXXFLAGS}"
        test ! -z "${LINALG_FCFLAGS}" && sd_linalg_fcflags="${LINALG_FCFLAGS}"
        test ! -z "${LINALG_LDFLAGS}" && sd_linalg_ldflags="${LINALG_LDFLAGS}"
        test ! -z "${LINALG_LIBS}" && sd_linalg_libs="${LINALG_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for linear algebra: '${sd_linalg_init}'])
        ;;

    esac

  fi

  # Export configuration
  AC_SUBST(sd_linalg_cppflags)
  AC_SUBST(sd_linalg_cflags)
  AC_SUBST(sd_linalg_cxxflags)
  AC_SUBST(sd_linalg_fcflags)
  AC_SUBST(sd_linalg_ldflags)
  AC_SUBST(sd_linalg_libs)
  AC_SUBST(sd_linalg_choices)
  AC_SUBST(sd_linalg_enable)
  AC_SUBST(sd_linalg_init)
  AC_SUBST(sd_linalg_ok)
]) # SD_LINALG_INIT


AC_DEFUN([SD_LINALG_INIT_FLAVOR], [
  # Init internal parameters
  sd_linalg_valid_flavors="auto acml asl atlas easybuild elpa essl magma mkl netlib none openblas plasma"
  sd_linalg_flavor=""
  sd_linalg_flavor_init="unknown"
  sd_linalg_flavor_gpu=""
  sd_linalg_flavor_mpi=""
  sd_linalg_flavor_mpiacc=""
  sd_linalg_flavor_serial=""

  # Set adjustable parameters
  sd_linalg_flavor_def="$1"

  # Set reasonable defaults if not provided
  test -z "${sd_linalg_flavor_def}" && sd_linalg_flavor_def="auto"

  # Declare configure option
  AC_ARG_WITH([linalg-flavor],
    [AS_HELP_STRING([--with-linalg-flavor],
      [Linear algebra flavor to select])],
    [ tmp_linalg_iter=`echo "${withval}" | sed -e 's/+/ /g'`
      for tmp_linalg_flavor in ${tmp_linalg_iter}; do
        tmp_linalg_flavor_ok=`echo "${sd_linalg_valid_flavors}" | grep "${tmp_linalg_flavor}"`
        if test "${tmp_linalg_flavor_ok}" = ""; then
          AC_MSG_ERROR([invalid linear algebra flavor: '${tmp_linalg_flavor}'])
        fi
      done
      sd_linalg_flavor="${withval}"
      sd_linalg_enable="yes"
      sd_linalg_flavor_init="kwd"
      unset tmp_linalg_iter
      unset tmp_linalg_flavor
      unset tmp_linalg_flavor_ok],
    [ sd_linalg_flavor="${sd_linalg_flavor_def}"
      sd_linalg_flavor_init="def"])

  # Check that the specified flavor is consistent
  _SD_LINALG_CHECK_FLAVOR

  # Export configuration
  AC_SUBST(sd_linalg_flavor)
  AC_SUBST(with_linalg_flavor)
])


# SD_LINALG_DETECT()
# ------------------
#
# Sets all variables needed to handle the optimized linear algebra
# libraries.
#
AC_DEFUN([SD_LINALG_DETECT], [
  AC_MSG_CHECKING([how to detect linear algebra libraries])
  case "${sd_linalg_init}" in
    def|yon)
      AC_MSG_RESULT([explore])
      _SD_LINALG_EXPLORE
      ;;
    dir|env)
      AC_MSG_RESULT([verify])
      _SD_LINALG_CHECK_LIBS
      ;;
    *)
      AC_MSG_ERROR([unsupported linear algebra init type: '${sd_linalg_init}'])
      ;;
  esac
  _SD_LINALG_DUMP_CONFIG

  # Define vendor-specific C preprocessing options
  tmp_linalg_iter=`echo "${sd_linalg_flavor}" | sed -e 's/+/ /g'`
  for tmp_linalg_vendor in ${tmp_linalg_iter}; do
    case "${tmp_linalg_vendor}" in
      asl)
        AC_DEFINE([HAVE_LINALG_ASL], 1,
          [Define to 1 if you have the ASL linear algebra library.])
        ;;
      essl)
        AC_DEFINE([HAVE_LINALG_ESSL], 1,
          [Define to 1 if you have the ESSL linear algebra library.])
        ;;
    esac
  done
  unset tmp_linalg_iter
  unset tmp_linalg_vendor
]) # SD_LINALG_DETECT


                    # ------------------------------------ #


#
# Internal macros
#


AC_DEFUN([_SD_LINALG_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_linalg_invalid="no"
  if test "${sd_linalg_enable_def}" != "auto" -a \
          "${sd_linalg_enable_def}" != "no" -a \
          "${sd_linalg_enable_def}" != "yes"; then
    case "${sd_linalg_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_linalg_enable_def = '${sd_linalg_enable_def}'])
        ;;
      skip)
        tmp_linalg_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_linalg_enable_def = '${sd_linalg_enable_def}'])
        tmp_linalg_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_linalg_invalid}" = "yes"; then
    if test "${sd_linalg_status}" = "required"; then
      sd_linalg_enable_def="yes"
    else
      sd_linalg_enable_def="no"
    fi
    tmp_linalg_invalid="no"
    AC_MSG_NOTICE([setting sd_linalg_enable_def to '${sd_linalg_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_linalg_invalid="no"
  if test "${sd_linalg_status}" = "implicit" -o \
          "${sd_linalg_status}" = "required"; then
    if test "${sd_linalg_enable}" = "no"; then
      case "${sd_linalg_policy}" in
        fail)
          AC_MSG_ERROR([The linear algebra package is required and cannot be disabled
                  See https://www.netlib.org/ for details on how to
                  install it.])
          ;;
        skip)
          tmp_linalg_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The linear algebra package is required and cannot be disabled])
          tmp_linalg_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_linalg_enable}" = "auto"; then
      AC_MSG_NOTICE([setting linear algebra trigger to yes])
      sd_linalg_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_linalg_invalid}" = "yes"; then
    case "${sd_linalg_status}" in
      implicit|required)
        sd_linalg_enable="yes"
        ;;
      optional)
        sd_linalg_enable="no"
        ;;
    esac
    tmp_linalg_invalid="no"
    AC_MSG_NOTICE([setting sd_linalg_enable to '${sd_linalg_enable}'])
  fi

  # Implicit status overrides everything
  if test "${sd_linalg_status}" = "implicit"; then
    if test "${sd_linalg_ldflags}" != ""; then
      sd_linalg_ldflags=""
      AC_MSG_NOTICE([resetting linear algebra linker flags (implicit package)])
    fi
    if test "${sd_linalg_libs}" != ""; then
      sd_linalg_libs=""
      AC_MSG_NOTICE([resetting linear algebra library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_linalg_enable}" = "implicit"; then
    sd_linalg_cppflags=""
    sd_linalg_cflags=""
    sd_linalg_fcflags=""
    sd_linalg_ldflags=""
    sd_linalg_libs=""
  fi

  # Clean-up
  unset tmp_linalg_invalid
  unset tmp_linalg_vars
]) # _SD_LINALG_CHECK_CONFIG


AC_DEFUN([_SD_LINALG_CHECK_FLAVOR], [
  # Display configured flavor
  AC_MSG_CHECKING([for the requested linear algebra flavor])
  AC_MSG_RESULT([${sd_linalg_flavor}])

  if test "${sd_linalg_flavor}" = "auto"; then

    # Make sure linear algebra is enabled
    sd_linalg_enable="yes"

    # Set generic flavors first
    sd_linalg_chk_serial="netlib"
    if test "${sd_mpi_enable}" = "yes"; then
      sd_linalg_chk_mpi="netlib"
      sd_linalg_chk_mpiacc="elpa"
    fi
    if test "${sd_gpu_enable}" = "yes"; then
      sd_linalg_chk_gpu="magma"
    fi

    # Refine with vendor-specific flavors
    case "${sd_fc_vendor}" in
      gnu)
        if test "${MKLROOT}" != ""; then
          sd_linalg_chk_serial="mkl ${sd_linalg_chk_serial}"
        fi
        sd_linalg_chk_serial="openblas atlas ${sd_linalg_chk_serial}"
        ;;
      intel)
        sd_linalg_chk_serial="mkl atlas ${sd_linalg_chk_serial}"
        sd_linalg_chk_mpi="mkl netlib"
        ;;
    esac

  else

    # Reformat flavor
    tmp_linalg_netlib_explicit=`echo "${sd_linalg_flavor}" | grep "netlib"`
    sd_linalg_iter=`echo "${sd_linalg_flavor}" | tr '+' '\n' | grep -v "netlib" | sort -u | awk '{printf " %s", [$]1} END{printf "\n"}'`
    sd_linalg_iter=`echo "${sd_linalg_iter}" | sed -e 's/^[ ]*//'`

    # The 'none' flavor excludes any other one
    tmp_linalg_none=`echo "${sd_linalg_iter}" | grep "none"`
    if test "${tmp_linalg_none}" != ""; then
      if test "${sd_linalg_iter}" != "none"; then
        AC_MSG_ERROR([the linear algebra flavor cannot be 'none' and
                  something else at the same time])
      fi
    fi
    unset tmp_linalg_none

    # Check flavor unicity for each detection sequence
    for tmp_linalg_flavor in ${sd_linalg_iter}; do
      case "${tmp_linalg_flavor}" in
        easybuild|mkl)
          if test "${sd_linalg_chk_serial}" != ""; then
            AC_MSG_ERROR([only one serial linear algebra flavor is permitted])
          fi
          sd_linalg_chk_serial="${tmp_linalg_flavor}"
          if test "${sd_mpi_enable}" = "yes"; then
            if test "${sd_linalg_chk_mpi}" != ""; then
              AC_MSG_ERROR([only one MPI linear algebra flavor is permitted])
            fi
            sd_linalg_chk_mpi="${tmp_linalg_flavor}"
          fi
          ;;
        elpa)
          if test "${sd_mpi_enable}" = "yes"; then
            if test "${sd_linalg_chk_mpiacc}" != ""; then
              AC_MSG_ERROR([only one MPI acceleration linear algebra flavor is permitted])
            fi
            sd_linalg_chk_mpiacc="${tmp_linalg_flavor}"
          else
            AC_MSG_NOTICE([ignoring '${tmp_linalg_flavor}', since MPI is disabled])
          fi
          ;;
        magma)
          if test "${sd_gpu_enable}" = "yes"; then
            if test "${sd_linalg_chk_gpu}" != ""; then
              AC_MSG_ERROR([only one GPU linear algebra flavor is permitted])
            fi
            sd_linalg_chk_gpu="${tmp_linalg_flavor}"
          else
            AC_MSG_NOTICE([ignoring '${tmp_linalg_flavor}', since GPU is disabled])
          fi
          ;;
        plasma)
          if test "${sd_mpi_enable}" = "yes"; then
            if test "${sd_linalg_chk_mpi}" != ""; then
              AC_MSG_ERROR([only one MPI linear algebra flavor is permitted])
            fi
            sd_linalg_chk_mpi="${tmp_linalg_flavor}"
          else
            AC_MSG_NOTICE([ignoring '${tmp_linalg_flavor}', since MPI is disabled])
          fi
          ;;
        *)
          if test "${sd_linalg_chk_serial}" != ""; then
            AC_MSG_ERROR([only one serial linear algebra flavor is permitted])
          fi
          sd_linalg_chk_serial="${tmp_linalg_flavor}"
          ;;
      esac
    done

    # Some vendors only provide a partial serial linear algebra support
    case "${sd_linalg_chk_serial}" in
      atlas|openblas)
        if test "${tmp_linalg_netlib_explicit}" = ""; then
          sd_linalg_chk_serial="${sd_linalg_chk_serial} netlib"
        fi
        ;;
    esac

    # Handle the explicit 'netlib' case
    if test "${tmp_linalg_netlib_explicit}" != ""; then
      case "${sd_linalg_chk_serial}" in
        atlas|openblas)
          sd_linalg_chk_serial="${sd_linalg_chk_serial} netlib"
          if test "${sd_mpi_enable}" = "yes"; then
            if test "${sd_linalg_chk_mpi}" = ""; then
              sd_linalg_chk_mpi="netlib"
            fi
          fi
          ;;
        *)
          if test "${sd_linalg_chk_serial}" = ""; then
            sd_linalg_chk_serial="netlib"
          fi
          ;;
      esac
    fi
    unset tmp_linalg_netlib_explicit

    # At this point we can assume that the serial detection sequence is set,
    # unless we want to bypass linear algebra tests
    if test "${sd_linalg_flavor}" != "none"; then
      if test "${sd_linalg_chk_serial}" = ""; then
        sd_linalg_chk_serial="netlib"
      fi
    fi

  fi   # sd_linalg_flavor = auto

  # Display detection sequences
  AC_MSG_CHECKING([for the serial linear algebra detection sequence])
  if test "${sd_linalg_chk_serial}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${sd_linalg_chk_serial}])
  fi
  if test "${sd_mpi_enable}" = "yes"; then
    AC_MSG_CHECKING([for the MPI linear algebra detection sequence])
    if test "${sd_linalg_chk_mpi}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_chk_mpi}])
    fi
    AC_MSG_CHECKING([for the MPI acceleration linear algebra detection sequence])
    if test "${sd_linalg_chk_mpiacc}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_chk_mpiacc}])
    fi
  fi
  if test "${sd_gpu_enable}" = "yes"; then
    AC_MSG_CHECKING([for the GPU linear algebra detection sequence])
    if test "${sd_linalg_chk_gpu}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_chk_gpu}])
    fi
  fi
]) # _SD_LINALG_CHECK_FLAVOR


AC_DEFUN([_SD_LINALG_DUMP_CONFIG], [
  if test "${sd_linalg_enable}" != "no"; then
    AC_MSG_CHECKING([how linear algebra parameters have been set])
    AC_MSG_RESULT([${sd_linalg_init} (flavor: ${sd_linalg_flavor_init})])
    AC_MSG_CHECKING([for the actual linear algebra flavor])
    AC_MSG_RESULT([${sd_linalg_flavor}])
    AC_MSG_CHECKING([for linear algebra C preprocessing flags])
    if test "${sd_linalg_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_cppflags}])
    fi
    AC_MSG_CHECKING([for linear algebra C flags])
    if test "${sd_linalg_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_cflags}])
    fi
    AC_MSG_CHECKING([for linear algebra C++ flags])
    if test "${sd_linalg_cxxflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_cxxflags}])
    fi
    AC_MSG_CHECKING([for linear algebra Fortran flags])
    if test "${sd_linalg_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_fcflags}])
    fi
    AC_MSG_CHECKING([for linear algebra linker flags])
    if test "${sd_linalg_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_ldflags}])
    fi
    AC_MSG_CHECKING([for linear algebra library flags])
    if test "${sd_linalg_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_linalg_libs}])
    fi
  fi
]) # _SD_LINALG_DUMP_CONFIG

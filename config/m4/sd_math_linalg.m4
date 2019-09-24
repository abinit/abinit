# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon, Marc Torrent)
#
# This file is part of the Steredeg software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Support for external linear algebra libraries
#


# SD_LINALG_INIT()
# ----------------
#
# Select a linear algebra flavor according to the available system information.
#
AC_DEFUN([SD_LINALG_INIT],[
  # Init build flags
  sd_linalg_cppflags=""
  sd_linalg_cflags=""
  sd_linalg_cxxflags=""
  sd_linalg_fcflags=""
  sd_linalg_ldflags=""
  sd_linalg_libs=""
  sd_linalg_enable=""
  sd_linalg_flavor=""
  sd_linalg_init="unknown"
  sd_linalg_ok="unknown"

  # Init internal variables
  sd_linalg_chk_gpu=""
  sd_linalg_chk_mpi=""
  sd_linalg_chk_mpiext=""
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
  sd_linalg_gpu="unknown"
  sd_linalg_mpi="unknown"
  sd_linalg_provided=""
  sd_linalg_valid_flavors="auto acml atlas custom easybuild elpa essl magma mkl netlib openblas plasma scalapack"

  # Set adjustable parameters
  sd_linalg_options="$1"
  sd_linalg_libs_def="$2"
  sd_linalg_cppflags_def="$3"
  sd_linalg_cflags_def="$4"
  sd_linalg_cxxflags_def="$5"
  sd_linalg_fcflags_def="$6"
  sd_linalg_ldflags_def="$7"
  sd_linalg_enable_def=""
  sd_linalg_flavor_def=""
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
  test -z "${sd_linalg_policy}" && sd_linalg_policy="fail"
  test -z "${sd_linalg_status}" && sd_linalg_status="optional"
  test -z "${sd_linalg_enable_def}" && sd_linalg_enable_def="no"
  test -z "${sd_linalg_flavor_def}" && sd_linalg_flavor_def="auto"
  case "${sd_linalg_status}" in
    implicit|required)
      sd_linalg_enable_def="yes"
      ;;
  esac

  # Declare configure option
  AC_ARG_WITH([linalg],
    AC_HELP_STRING([--with-linalg],
      [Install prefix of the linear algebra libraries (e.g. /usr/local), or automatic configuration if set to 'yes' or 'no' (default: yes)]),
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_linalg_enable="${withval}"
        sd_linalg_init="yon"
      else
        sd_linalg_enable="yes"
        sd_linalg_init="dir"
      fi],
    [ sd_linalg_enable="${sd_linalg_enable_def}"; sd_linalg_init="def"])

  # Declare flavoring configure option
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
      sd_linalg_init="kwd"
      unset tmp_linalg_iter
      unset tmp_linalg_flavor
      unset tmp_linalg_flavor_ok],
    [ sd_linalg_enable="${sd_linalg_enable_def}"
      sd_linalg_flavor="${sd_linalg_flavor_def}"
      sd_linalg_init="def"])

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

      def|yon|kwd)
        sd_linalg_cppflags="${sd_linalg_cppflags_def}"
        sd_linalg_cflags="${sd_linalg_cflags_def}"
        sd_linalg_cxxflags="${sd_linalg_cxxflags_def}"
        sd_linalg_fcflags="${sd_linalg_fcflags_def}"
        sd_linalg_ldflags="${sd_linalg_ldflags_def}"
        sd_linalg_libs="${sd_linalg_libs_def}"
        ;;

      dir)
        sd_linalg_cppflags="-I${with_linalg}/include"
        sd_linalg_cflags="${sd_linalg_cflags_def}"
        sd_linalg_cxxflags="${sd_linalg_cxxflags_def}"
        sd_linalg_fcflags="${sd_linalg_fcflags_def}"
        sd_linalg_ldflags="${sd_linalg_ldflags_def}"
        sd_linalg_libs="-L${with_linalg}/lib ${sd_linalg_libs_def}"
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

  # Check that the specified flavor is consistent
  if test "${sd_linalg_flavor}" != "auto"; then

    # Reformat flavor
    sd_linalg_iter=`echo "${sd_linalg_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s",[$]1}'`

    # Check serial and parallel flavor unicity
    for tmp_linalg_flavor in ${sd_linalg_iter}; do
      case "${tmp_linalg_flavor}" in
        magma)
          if test "${sd_linalg_chk_gpu}" != ""; then
            AC_MSG_ERROR([only one GPU linear algebra flavor is permitted])
          fi
          sd_linalg_chk_gpu="${tmp_linalg_flavor}"
          ;;
        scalapack|plasma)
          if test "${sd_linalg_chk_mpi}" != ""; then
            AC_MSG_ERROR([only one MPI linear algebra flavor is permitted])
          fi
          sd_linalg_chk_mpi="${tmp_linalg_flavor}"
          ;;
        elpa)
          sd_linalg_chk_mpiext="${tmp_linalg_flavor}"
          ;;
        *)
          if test "${sd_linalg_chk_serial}" != ""; then
            AC_MSG_ERROR([only one serial linear algebra flavor is permitted])
          fi
          sd_linalg_chk_serial="${tmp_linalg_flavor}"
          ;;
      esac
      _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_flavor}])
    done
    if test "${sd_linalg_chk_serial}" = ""; then
      AC_MSG_ERROR([you must choose a serial linear algebra flavor])
    fi

  fi   # sd_linalg_flavor != auto

  # Init specific flavors
  AC_MSG_CHECKING([for the requested linear algebra support])
  AC_MSG_RESULT([${sd_linalg_flavor}])
  case "${sd_linalg_flavor}" in

    auto)
      # Set generic flavors first
      sd_linalg_chk_serial="netlib"
      if test "${sd_mpi_enable}" = "yes"; then
        sd_linalg_chk_mpi="elpa scalapack"
      fi
      if test "${sd_gpu_enable}" = "yes"; then
        sd_linalg_chk_gpu="magma"
      fi

      # Refine with vendor-specific flavors
      case "${sd_fc_vendor}" in
        gnu)
          sd_linalg_chk_serial="atlas netlib"
          ;;
        intel)
          sd_linalg_chk_serial="mkl atlas netlib"
          ;;
        *)
          sd_linalg_chk_serial="netlib"
          ;;
      esac
      ;;

    none)
      AC_MSG_WARN([bypassing linear algebra tests])
      sd_linalg_has_blas="yes"
      sd_linalg_has_lapack="yes"
      sd_linalg_serial="yes"
      sd_linalg_mpi="no"
      sd_linalg_gpu="no"
      ;;

  esac

  # Display configuration
  _SD_LINALG_DUMP_CONFIG

  # Display detection sequences
  AC_MSG_CHECKING([for the serial linear algebra detection sequence])
  if test "${sd_linalg_chk_serial}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${sd_linalg_chk_serial}])
  fi
  AC_MSG_CHECKING([for the MPI linear algebra detection sequence])
  if test "${sd_linalg_chk_mpi}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${sd_linalg_chk_mpi}])
  fi
  AC_MSG_CHECKING([for the GPU linear algebra detection sequence])
  if test "${sd_linalg_chk_gpu}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${sd_linalg_chk_gpu}])
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
  AC_SUBST(sd_linalg_flavor)
  AC_SUBST(sd_linalg_init)
  AC_SUBST(sd_linalg_ok)
  AC_SUBST(with_linalg_flavor)

]) # SD_LINALG_INIT


# SD_LINALG_DETECT()
# ------------------
#
# Sets all variables needed to handle the optimized linear algebra
# libraries.
#
AC_DEFUN([SD_LINALG_DETECT], [
  # Explore the system for if the user did not specify anything
  if test "${sd_linalg_init}" = "def" -o "${sd_linalg_init}" = "yon"; then
    _SD_LINALG_EXPLORE
  else
    _SD_LINALG_CHECK_LIBS
  fi
]) # SD_LINALG_DETECT


                    # ------------------------------------ #


#
# Private macros
#


# _SD_LINALG_CHECK_LIBS()
# -----------------------
#
# Check whether the specified libraries are BLAS and LAPACK
# implementations.
#
AC_DEFUN([_SD_LINALG_CHECK_LIBS], [
  # Init
  sd_linalg_has_blas="no"
  sd_linalg_has_lapack="no"
  sd_linalg_has_lapacke="no"
  sd_linalg_has_blacs="no"
  sd_linalg_has_scalapack="no"
  sd_linalg_has_elpa="no"
  sd_linalg_has_elpa_2013="no"
  sd_linalg_has_elpa_2014="no"
  sd_linalg_has_elpa_2015="no"
  sd_linalg_has_elpa_2016="no"
  sd_linalg_has_plasma="no"
  sd_linalg_has_magma="no"

  # Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${sd_linalg_cppflags}"
  FCFLAGS="${FCFLAGS} ${sd_linalg_fcflags}"
  LIBS="${sd_linalg_libs} ${sd_gpu_libs} ${sd_mpi_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  # BLAS?
  AC_MSG_CHECKING([for BLAS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
    ]])], [sd_linalg_has_blas="yes"], [sd_linalg_has_blas="no"])
  AC_MSG_RESULT([${sd_linalg_has_blas}])

  # BLAS extensions?
  _SD_LINALG_CHECK_BLAS_EXTS()

  # MKL BLAS extensions?
  _SD_LINALG_CHECK_BLAS_MKL_EXTS()

  # LAPACK?
  AC_MSG_CHECKING([for LAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zhpev
    ]])], [sd_linalg_has_lapack="yes"], [sd_linalg_has_lapack="no"])
  AC_MSG_RESULT([${sd_linalg_has_lapack}])

  # LAPACKE?
  AC_MSG_CHECKING([for LAPACKE C API support in specified libraries])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [#include <lapacke.h>],
    [[
      zhpev_;
    ]])],[sd_linalg_has_lapacke="yes"], [sd_linalg_has_lapacke="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_linalg_has_lapacke}])

  # BLACS?
  AC_MSG_CHECKING([for BLACS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call blacs_gridinit
    ]])], [sd_linalg_has_blacs="yes"], [sd_linalg_has_blacs="no"])
  AC_MSG_RESULT([${sd_linalg_has_blacs}])

  # ScaLAPACK?
  AC_MSG_CHECKING([for ScaLAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [sd_linalg_has_scalapack="yes"], [sd_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${sd_linalg_has_scalapack}])

  # ELPA
  AC_MSG_CHECKING([for ELPA support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer,parameter :: n=1, comm=1
      integer :: comm1, comm2, success
      success = get_elpa_communicators(comm, n, n, comm1, comm2)
    ]])], [sd_linalg_has_elpa="yes"], [sd_linalg_has_elpa="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa}])
  if test "${sd_linalg_has_elpa}" = "yes"; then
    _SD_LINALG_CHECK_ELPA_2017()
    _SD_LINALG_CHECK_ELPA_2016()
    _SD_LINALG_CHECK_ELPA_2015()
    _SD_LINALG_CHECK_ELPA_2014()
    _SD_LINALG_CHECK_ELPA_2013()
  fi

  # PLASMA?
  AC_MSG_CHECKING([for PLASMA support in specified libraries])
  sd_linalg_chk_plasma="${sd_linalg_has_lapacke}"
  if test "${sd_linalg_chk_plasma}" = "yes"; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use plasma
        call plasma_zhegv
      ]])], [sd_linalg_has_plasma="yes"], [sd_linalg_has_plasma="no"])
  fi
  AC_MSG_RESULT([${sd_linalg_has_plasma}])

  # MAGMA?
  AC_MSG_CHECKING([for MAGMA (version>=1.1.0) support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call magmaf_zhegvd
    ]])], [sd_linalg_has_magma="yes"], [sd_linalg_has_magma="no"])
  AC_MSG_RESULT([${sd_linalg_has_magma}])
  if test "${sd_linalg_has_magma}" = "yes"; then
    _SD_LINALG_CHECK_MAGMA_15()
  fi

  # Restore environment
  AC_LANG_POP([Fortran])
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # _SD_LINALG_CHECK_LIBS


                    # ------------------------------------ #


# _SD_LINALG_EXPLORE()
# --------------------
#
# Looks for linear algebra components by going through all the selected
# vendor sequences.
#
AC_DEFUN([_SD_LINALG_EXPLORE], [
  # Prepare environment
  SD_ENV_BACKUP
  LDFLAGS="${FC_LDFLAGS}"
  sd_saved_CPPFLAGS="${CPPFLAGS}"
  sd_saved_CXXFLAGS="${CXXFLAGS}"
  sd_saved_CFLAGS="${CFLAGS}"
  sd_saved_FCFLAGS="${FCFLAGS}"
  sd_saved_LDFLAGS="${LDFLAGS}"
  sd_saved_LIBS="${LIBS}"

  # Look for serial linear algebra support
  for tmp_linalg_vendor in ${sd_linalg_chk_serial}; do

    # Configure vendor libraries
    _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])
    CPPFLAGS="${sd_saved_CPPFLAGS} ${sd_linalg_vendor_cppflags}"
    CFLAGS="${sd_saved_CFLAGS} ${sd_linalg_vendor_cflags}"
    CXXFLAGS="${sd_saved_CXXFLAGS} ${sd_linalg_vendor_cxxflags}"
    FCFLAGS="${sd_saved_FCFLAGS} ${sd_linalg_vendor_fcflags}"
    LDFLAGS="${sd_saved_LDFLAGS} ${sd_linalg_vendor_ldflags}"

    # Look for BLAS
    tmp_linalg_blas_proceed=`echo "${sd_linalg_vendor_provided}" | grep "blas"`
    if test "${tmp_linalg_blas_proceed}" != "" -a \
            "${sd_linalg_has_blas}" != "yes"; then
      LIBS="${sd_linalg_vendor_blas_libs} ${sd_linalg_vendor_blas_prqs} ${sd_saved_LIBS}"
      AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for BLAS])
      if test "${sd_linalg_vendor_blas_libs}${sd_linalg_vendor_blas_prqs}" = ""; then
        AC_MSG_RESULT([none required])
      else
        AC_MSG_RESULT([${sd_linalg_vendor_blas_libs} ${sd_linalg_vendor_blas_prqs}])
      fi
      _SD_LINALG_CHECK_BLAS
      if test "${sd_linalg_has_blas}" = "yes"; then
         #sd_linalg_blas_cppflags="${sd_linalg_vendor_blas_cppflags}"
         #sd_linalg_blas_cflags="${sd_linalg_vendor_blas_cflags}"
         #sd_linalg_blas_fcflags="${sd_linalg_vendor_blas_fcflags}"
         #sd_linalg_blas_ldflags="${sd_linalg_vendor_blas_ldflags}"
         sd_linalg_blas_libs="${sd_linalg_vendor_blas_libs} ${sd_linalg_vendor_blas_prqs}"
         sd_linalg_blas_vendor="${tmp_linalg_vendor}"
         sd_linalg_provided="${sd_linalg_provided} blas"
         _SD_LINALG_CHECK_BLAS_EXTS
         if test "${tmp_linalg_vendor}" = "mkl"; then
           _SD_LINALG_CHECK_BLAS_MKL_EXTS
         fi
      fi
    fi

    # Look for LAPACK
    tmp_linalg_lapack_proceed=`echo "${sd_linalg_vendor_provided}" | grep "lapack"`
    if test "${tmp_linalg_lapack_proceed}" != "" -a \
            "${sd_linalg_has_blas}" = "yes" -a \
            "${sd_linalg_has_lapack}" != "yes"; then

      AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for LAPACK])
      if test "${sd_linalg_vendor_lapack_libs}${sd_linalg_vendor_lapack_prqs}" = ""; then
        AC_MSG_RESULT([none required])
      else
       AC_MSG_RESULT([${sd_linalg_vendor_lapack_libs} ${sd_linalg_vendor_lapack_prqs}])
      fi
      LIBS="${sd_linalg_vendor_lapack_libs} ${sd_linalg_vendor_lapack_prqs} ${sd_linalg_blas_libs} ${sd_saved_LIBS}"
      _SD_LINALG_CHECK_LAPACK
      if test "${sd_linalg_has_lapack}" = "yes"; then
         sd_linalg_lapack_libs="${sd_linalg_vendor_lapack_libs} ${sd_linalg_vendor_lapack_prqs}"
         sd_linalg_lapack_vendor="${tmp_linalg_vendor}"
         sd_linalg_provided="${sd_linalg_provided} lapack"
      fi
    fi

  done

  # FIXME: Look for MPI linear algebra support
  if test "${sd_mpi_enable}" = "yes"; then
    for tmp_linalg_vendor in ${sd_linalg_chk_mpi}; do
      AC_MSG_WARN([MPI linear algebra exploration not implemented!])
    done
  fi

  # FIXME: Look for GPU linear algebra support
  if test "${sd_gpu_enable}" = "yes"; then
    for tmp_linalg_vendor in ${sd_linalg_chk_gpu}; do
      AC_MSG_WARN([GPU linear algebra exploration not implemented!])
    done
  fi

  # Transmit linear algebra parameters found
  sd_linalg_cppflags="${sd_linalg_vendor_cppflags}"
  sd_linalg_cflags="${sd_linalg_vendor_cflags}"
  sd_linalg_cxxflags="${sd_linalg_vendor_cxxflags}"
  sd_linalg_fcflags="${sd_linalg_vendor_fcflags}"
  sd_linalg_ldflags="${sd_linalg_vendor_ldflags}"
  if test "${sd_mpi_enable}" = "yes"; then
    sd_linalg_libs="${sd_linalg_libs} ${sd_linalg_vendor_scalapack_libs}"
  fi
  sd_linalg_libs="${sd_linalg_libs} ${sd_linalg_vendor_lapack_libs} ${sd_linalg_vendor_blas_libs}"

  SD_ENV_RESTORE
  LIBS="${sd_saved_LIBS}"
]) # SD_LINALG_EXPLORE


                    # ------------------------------------ #
                    # ------------------------------------ #


# _SD_LINALG_CHECK_BLAS()
# -----------------------
#
# Check whether the build environment provides BLAS.
#
AC_DEFUN([_SD_LINALG_CHECK_BLAS], [
  sd_linalg_has_blas="unknown"

  AC_MSG_CHECKING([for BLAS support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
    ]])], [sd_linalg_has_blas="yes"], [sd_linalg_has_blas="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_blas}])
]) # _SD_LINALG_CHECK_BLAS


# _SD_LINALG_CHECK_BLAS_EXTS()
# ----------------------------
#
# Check whether the specified BLAS implementation provides useful extensions.
#
AC_DEFUN([_SD_LINALG_CHECK_BLAS_EXTS], [
  # AXPBY family?
  AC_MSG_CHECKING([for AXPBY support in the BLAS libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call saxpby
      call daxpby
      call caxpby
      call zaxpby
    ]])], [sd_linalg_has_axpby="yes"], [sd_linalg_has_axpby="no"])
  AC_MSG_RESULT([${sd_linalg_has_axpby}])

  if test "${sd_linalg_has_axpby}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_AXPBY], 1,
      [Define to 1 if you have an AXPBY BLAS1 extensions.])
  fi

  # gemm3m family
  AC_MSG_CHECKING([for GEMM3M in the BLAS libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call cgemm3m
      call zgemm3m
    ]])], [sd_linalg_has_gemm3m="yes"], [sd_linalg_has_gemm3m="no"])
  AC_MSG_RESULT([${sd_linalg_has_gemm3m}])

  if test "${sd_linalg_has_gemm3m}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_GEMM3M], 1,
      [Define to 1 if you have ?GEMM3M BLAS3 extensions.])
  fi
]) # _SD_LINALG_CHECK_BLAS_EXTS


# _SD_LINALG_CHECK_BLAS_MKL_EXTS()
# --------------------------------
#
# Check whether the specified MKL implementation provides BLAS extensions.
#
AC_DEFUN([_SD_LINALG_CHECK_BLAS_MKL_EXTS], [
  # mkl_imatcopy family
  AC_MSG_CHECKING([for mkl_imatcopy in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call mkl_simatcopy
      call mkl_dimatcopy
      call mkl_cimatcopy
      call mkl_zimatcopy
    ]])], [sd_linalg_mkl_has_imatcopy="yes"], [sd_linalg_mkl_has_imatcopy="no"])
  AC_MSG_RESULT([${sd_linalg_mkl_has_imatcopy}])

  if test "${sd_linalg_mkl_has_imatcopy}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_IMATCOPY], 1,
      [Define to 1 if you have mkl_?imatcopy extensions.])
  fi

  # mkl_omatcopy family
  AC_MSG_CHECKING([for mkl_omatcopy in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call mkl_somatcopy
      call mkl_domatcopy
      call mkl_comatcopy
      call mkl_zomatcopy
    ]])], [sd_linalg_mkl_has_omatcopy="yes"], [sd_linalg_mkl_has_omatcopy="no"])
  AC_MSG_RESULT([${sd_linalg_mkl_has_omatcopy}])

  if test "${sd_linalg_mkl_has_omatcopy}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_OMATCOPY], 1,
      [Define to 1 if you have mkl_?omatcopy extensions.])
  fi

  # mkl_omatadd family
  AC_MSG_CHECKING([for mkl_omatadd in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call mkl_somatadd
      call mkl_comatadd
      call mkl_domatadd
      call mkl_zomatadd
    ]])], [sd_linalg_mkl_has_omatadd="yes"], [sd_linalg_mkl_has_omatadd="no"])
  AC_MSG_RESULT([${sd_linalg_mkl_has_omatadd}])

  if test "${sd_linalg_mkl_has_omatadd}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_OMATADD], 1,
      [Define to 1 if you have mkl_?omatadd extensions.])
  fi

  # mkl_threads support functions
  AC_MSG_CHECKING([for mkl_set/get_threads in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      integer :: a
      a = mkl_get_max_threads()
      call mkl_set_num_threads
    ]])], [sd_linalg_mkl_has_threads="yes"], [sd_linalg_mkl_has_threads="no"])
  AC_MSG_RESULT([${sd_linalg_mkl_has_threads}])

  if test "${sd_linalg_mkl_has_threads}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_THREADS], 1,
      [Define to 1 if you have mkl_*threads extensions.])
  fi
]) # _SD_LINALG_CHECK_BLAS_MKL_EXTS


                    # ------------------------------------ #


# _SD_LINALG_CHECK_LAPACK()
# -------------------------
#
# Check whether the build environment provides LAPACK.
#
AC_DEFUN([_SD_LINALG_CHECK_LAPACK], [
  sd_linalg_has_lapack="unknown"

  AC_MSG_CHECKING([for LAPACK support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zhpev
    ]])], [sd_linalg_has_lapack="yes"], [sd_linalg_has_lapack="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_lapack}])
]) # _SD_LINALG_CHECK_LAPACK


                    # ------------------------------------ #


# _SD_LINALG_CHECK_ELPA_2013()
# ----------------------------
#
# Look for a ELPA 2013 API.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA_2013], [
  # Init
  sd_linalg_has_elpa_2013="no"

  AC_MSG_CHECKING([for ELPA 2013 API support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer, parameter :: na=1, lda=1, ldq=1, nev=1, nblk=1
      integer, parameter :: comm_r=1, comm_c=1
      real*8 :: a(lda,na), ev(na), q(ldq,na)
      complex*16 :: ac(lda,na)
      call solve_evp_real(na, nev, a, lda, ev, q, ldq, nblk, comm_r, comm_c)
      call invert_trm_complex(na, ac, lda, nblk, comm_r, comm_c)
    ]])], [sd_linalg_has_elpa_2013="yes"], [sd_linalg_has_elpa_2013="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa_2013}])

  if test "${sd_linalg_has_elpa_2013}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2013], 1,
      [Define to 1 if you have ELPA 2013 API support.])
  fi
]) # _SD_LINALG_CHECK_ELPA_2013


# _SD_LINALG_CHECK_ELPA_2014()
# ----------------------------
#
# Look for a ELPA 2014 API.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA_2014], [
  # Init
  sd_linalg_has_elpa_2014="no"

  AC_MSG_CHECKING([for ELPA 2014 API support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      logical :: success
      integer, parameter :: na=1, lda=1, ldq=1, nev=1, nblk=1
      integer, parameter :: comm_r=1, comm_c=1
      real*8 :: a(lda,na), ev(na), q(ldq,na)
      complex*16 :: ac(lda,na)
      success = solve_evp_real(na, nev, a, lda, ev, q, ldq, nblk, &
&       comm_r, comm_c)
      call invert_trm_complex(na, ac, lda, nblk, comm_r, comm_c, success)
    ]])], [sd_linalg_has_elpa_2014="yes"], [sd_linalg_has_elpa_2014="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa_2014}])

  if test "${sd_linalg_has_elpa_2014}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2014], 1,
      [Define to 1 if you have ELPA 2014 API support.])
  fi
]) # _SD_LINALG_CHECK_ELPA_2014


# _SD_LINALG_CHECK_ELPA_2015()
# ----------------------------
#
# Look for a ELPA 2015 API.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA_2015], [
  # Init
  sd_linalg_has_elpa_2015="no"

  AC_MSG_CHECKING([for ELPA 2015 API support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      logical :: success1, debug
      integer, parameter :: na=1, lda=1, ldq=1, nev=1, nblk=1
      integer :: comm_g=1, comm_r=1, comm_c=1, success2
      real*8 :: a(lda,na), ev(na), q(ldq,na)
      complex*16 :: ac(lda,na)
      success1 = solve_evp_real(na, nev, a, lda, ev, q, ldq, nblk, &
&       comm_r, comm_c)
      call cholesky_complex(na, ac, lda, nblk, comm_r, comm_c, debug, success1)
      success2 = get_elpa_row_col_comms(comm_g, na, na, comm_r, comm_c)
    ]])], [sd_linalg_has_elpa_2015="yes"], [sd_linalg_has_elpa_2015="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa_2015}])

  if test "${sd_linalg_has_elpa_2015}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2015], 1,
      [Define to 1 if you have ELPA 2015 API support.])
  fi
]) # _SD_LINALG_CHECK_ELPA_2015


# _SD_LINALG_CHECK_ELPA_2016()
# ----------------------------
#
# Look for a ELPA 2016 API.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA_2016], [
  # Init
  sd_linalg_has_elpa_2016="no"

  AC_MSG_CHECKING([for ELPA 2016 API support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      logical :: success1, debug
      integer, parameter :: na=1, lda=1, ldq=1, nev=1, nblk=1, nrow=1
      integer :: comm_g=1, comm_r=1, comm_c=1, success2
      real*8 :: a(lda,nrow), ev(na), q(ldq,nrow)
      complex*16 :: ac(lda,nrow)
      success1 = solve_evp_real_1stage(na, nev, a, lda, ev, q, ldq, &
&       nblk, nrow, comm_r, comm_c)
      success1 = cholesky_complex(na, ac, lda, nblk, nrow, &
&       comm_r, comm_c, debug)
      success2 = get_elpa_communicators(comm_g, na, na, comm_r, comm_c)
    ]])], [sd_linalg_has_elpa_2016="yes"], [sd_linalg_has_elpa_2016="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa_2016}])

  if test "${sd_linalg_has_elpa_2016}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2016], 1,
      [Define to 1 if you have ELPA 2016 API support.])
  fi
]) # _SD_LINALG_CHECK_ELPA_2016


# _SD_LINALG_CHECK_ELPA_2017()
# ----------------------------
#
# Look for a ELPA 2017+ API.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA_2017], [
  # Init
  sd_linalg_has_elpa_2017="no"

  AC_MSG_CHECKING([for ELPA 2017 API support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa
      class(elpa_t),pointer :: e
      integer,parameter :: na=1,ncol=1,nrow=1 ; integer :: err
      real*8 :: a(ncol,nrow),ev(na),q(ncol,nrow)
      call e%eigenvectors(a,ev,q,err)
      call e%cholesky(a,err)
    ]])], [sd_linalg_has_elpa_2017="yes"], [sd_linalg_has_elpa_2017="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa_2017}])

  if test "${sd_linalg_has_elpa_2017}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2017], 1,
      [Define to 1 if you have ELPA 2017 API support.])
    AC_DEFINE([HAVE_ELPA_FORTRAN2008], 1,
      [Define to 1 if you have ELPA Fortran 2008 API support.])
  fi
]) # _SD_LINALG_CHECK_ELPA_2017


                    # ------------------------------------ #


# _SD_LINALG_CHECK_MAGMA_15()
# ---------------------------
#
# Look for MAGMA >=1.5 (requires magma_init and magma_finalize).
#
AC_DEFUN([_SD_LINALG_CHECK_MAGMA_15], [
  # Init
  sd_linalg_has_magma_15="no"

  AC_MSG_CHECKING([for magma_init/magma_finalize support in specified MAGMA libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call magmaf_init
      call magma_finalize
    ]])], [sd_linalg_has_magma_15="yes"], [sd_linalg_has_magma_15="no"])
  AC_MSG_RESULT([${sd_linalg_has_magma_15}])

  if test "${sd_linalg_has_magma_15}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MAGMA_15], 1,
      [Define to 1 if you have MAGMA >=1.5 API support])
  fi
]) # _SD_LINALG_CHECK_MAGMA_15


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_BLACS(BLACS, EXTRA_LIBS)
# ------------------------------------------
#
# Look for a BLACS implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_BLACS], [
  # Init
  sd_linalg_has_blacs="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([blacs_gridinit], $1,
    [sd_linalg_has_blacs="yes"], [sd_linalg_has_blacs="no"],
    [$2 ${sd_linalg_libs}])
  if test "${sd_linalg_has_blacs}" = "yes"; then
    if test "${ac_cv_search_blacs_gridinit}" != "none required"; then
      sd_linalg_libs="${ac_cv_search_blacs_gridinit} $2 ${sd_linalg_libs}"
    fi
  fi
]) # _SD_LINALG_SEARCH_BLACS


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_BLAS(BLAS, EXTRA_LIBS)
# ----------------------------------------
#
# Look for a BLAS implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_BLAS], [
  # Init
  sd_linalg_has_blas="unknown"

  # Look for libraries and routines
  AC_MSG_CHECKING([for libraries that may contain BLAS]) 
  AC_MSG_RESULT([$1])
  AC_LANG_PUSH([Fortran])
  AC_SEARCH_LIBS([zgemm], $1,
    [sd_linalg_has_blas="yes"], [sd_linalg_has_blas="no"],
    [$2 ${sd_linalg_libs}])
  AC_LANG_POP([Fortran])
  if test "${sd_linalg_has_blas}" = "yes"; then
    if test "${ac_cv_search_zgemm}" != "none required"; then
      sd_linalg_libs="${ac_cv_search_zgemm} $2 ${sd_linalg_libs}"
    fi
  fi
]) # _SD_LINALG_SEARCH_BLAS


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_ELPA(ELPA, EXTRA_LIBS)
# ----------------------------------------
#
# Look for a ELPA implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_ELPA], [
  # Init
  sd_linalg_has_elpa="no"

  # Look for libraries and routines
  # Has to rewrite AC_SEARCH_LIBS because of mandatory F90 module
  AC_MSG_CHECKING([for the ELPA library])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer, parameter :: n=1, comm=1
      integer :: comm1, comm2, success
      success = get_elpa_communicators(comm, n, n, comm1, comm2)
    ]])], [sd_linalg_has_elpa="yes"], [sd_linalg_has_elpa="no"])
  if test "${sd_linalg_has_elpa}" = "no"; then
    tmp_saved_LIBS="${LIBS}"
    for test_lib in $1; do
      LIBS="-l${test_lib} $2 ${tmp_saved_LIBS}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use elpa1
          integer, parameter :: n=1, comm=1
          integer :: comm1, comm2, success
          success = get_elpa_communicators(comm, n, n, comm1, comm2)
        ]])], [sd_linalg_has_elpa="yes"], [sd_linalg_has_elpa="no"])
      if test "${sd_linalg_has_elpa}" = "yes"; then
        sd_linalg_libs="-l${test_lib} $2 ${sd_linalg_libs}"
        break
      fi
      LIBS="${tmp_saved_LIBS}"
    done
  fi
  AC_MSG_RESULT([${sd_linalg_has_elpa}])

  if test "${sd_linalg_has_elpa}" = "yes"; then
    _SD_LINALG_CHECK_ELPA_2017()
    _SD_LINALG_CHECK_ELPA_2016()
    _SD_LINALG_CHECK_ELPA_2015()
    _SD_LINALG_CHECK_ELPA_2014()
    _SD_LINALG_CHECK_ELPA_2013()
  fi
]) # _SD_LINALG_SEARCH_ELPA


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_LAPACK(LAPACK, EXTRA_LIBS)
# --------------------------------------------
#
# Look for a LAPACK implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_LAPACK], [
  # Init
  sd_linalg_has_lapack="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([zhpev], $1,
    [sd_linalg_has_lapack="yes"], [sd_linalg_has_lapack="no"],
    [$2 ${sd_linalg_libs}])
  if test "${sd_linalg_has_lapack}" = "yes"; then
    if test "${ac_cv_search_zhpev}" != "none required"; then
      sd_linalg_libs="${ac_cv_search_zhpev} $2 ${sd_linalg_libs}"
    fi
  fi
]) # _SD_LINALG_SEARCH_LAPACK


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_LAPACKE(LAPACKE, EXTRA_LIBS)
# ----------------------------------------------
#
# Look for a LAPACKE C API implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_LAPACKE], [
  # Init
  sd_linalg_has_lapacke="no"

  # Look for libraries and routines
  AC_MSG_CHECKING([for library containing zhpev_ C API])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [#include <lapacke.h>],
    [[
      zhpev_;
    ]])], [sd_linalg_has_lapacke="yes"], [sd_linalg_has_lapacke="no"])
  if test "${sd_linalg_has_lapacke}" = "no"; then
    tmp_saved_LIBS="${LIBS}"
    for test_lib in $1; do
      LIBS="-l${test_lib} $2 ${tmp_saved_LIBS}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [#include <lapacke.h>],
        [[
          zhpev_;
        ]])], [sd_linalg_has_lapacke="yes"], [sd_linalg_has_lapacke="no"])
      if test "${sd_linalg_has_lapacke}" = "yes"; then
        sd_linalg_libs="-l${test_lib} $2 ${sd_linalg_libs}"
        break
      fi  
    done
    if test "${sd_linalg_has_lapacke}" = "no"; then
      LIBS="${tmp_saved_LIBS}"
    fi
  fi
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_linalg_has_lapacke}])
]) # _SD_LINALG_SEARCH_LAPACKE


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_MAGMA(MAGMA, EXTRA_LIBS)
# ------------------------------------------
#
# Look for a MAGMA implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_MAGMA], [
  # Init
  sd_linalg_has_magma="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([magmaf_zheevd], $1,
    [sd_linalg_has_magma="yes"], [sd_linalg_has_magma="no"],
    [$2 ${sd_linalg_libs}])
  if test "${sd_linalg_has_magma}" = "yes"; then
    if test "${ac_cv_search_magmaf_zheevd}" != "none required"; then
      sd_linalg_libs="${ac_cv_search_magmaf_zheevd} $2 ${sd_linalg_libs}"
    fi
    _SD_LINALG_CHECK_MAGMA_15()
  fi
]) # _SD_LINALG_SEARCH_MAGMA


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_PLASMA(PLASMA, EXTRA_LIBS)
# --------------------------------------------
#
# Look for a PLASMA implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_PLASMA], [
  # Init
  sd_linalg_has_plasma="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([plasma_zhegv], $1,
    [sd_linalg_has_plasma="yes"], [sd_linalg_has_plasma="no"],
    [$2 ${sd_linalg_libs}])
  if test "${sd_linalg_has_plasma}" = "yes"; then
    if test "${ac_cv_search_plasma_zhegv}" != "none required"; then
      sd_linalg_libs="${ac_cv_search_plasma_zhegv} $2 ${sd_linalg_libs}"
    fi
  fi
]) # _SD_LINALG_SEARCH_PLASMA


                    # ------------------------------------ #


# _SD_LINALG_SEARCH_SCALAPACK(SCALAPACK, EXTRA_LIBS)
# --------------------------------------------------
#
# Look for a ScaLAPACK implementation.
#
AC_DEFUN([_SD_LINALG_SEARCH_SCALAPACK], [
  # Init
  sd_linalg_has_scalapack="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([pzheevx], $1,
    [sd_linalg_has_scalapack="yes"], [sd_linalg_has_scalapack="no"],
    [$2 ${sd_linalg_libs}])
  if test "${sd_linalg_has_scalapack}" = "yes"; then
    if test "${ac_cv_search_pzheevx}" != "none required"; then
      sd_linalg_libs="${ac_cv_search_pzheevx} $2 ${sd_linalg_libs}"
    fi
  fi
]) # _SD_LINALG_SEARCH_SCALAPACK


                    # ------------------------------------ #


# _SD_LINALG_SET_VENDOR_FLAGS(VENDOR)
# -----------------------------------
#
# Set libraries to look for depending on the specified flavor.
#
AC_DEFUN([_SD_LINALG_SET_VENDOR_FLAGS], [
  # Reset components
  sd_linalg_vendor_provided=""
  sd_linalg_vendor_cppflags=""
  sd_linalg_vendor_cflags=""
  sd_linalg_vendor_cxxflags=""
  sd_linalg_vendor_fcflags=""
  sd_linalg_vendor_ldflags=""
  sd_linalg_vendor_blas_libs=""
  sd_linalg_vendor_blas_prqs=""
  sd_linalg_vendor_lapack_libs=""
  sd_linalg_vendor_lapack_prqs=""
  sd_linalg_vendor_lapacke_libs=""
  sd_linalg_vendor_lapacke_prqs=""
  sd_linalg_vendor_scalapack_libs=""
  sd_linalg_vendor_scalapack_prqs=""
  sd_linalg_vendor_elpa_libs=""
  sd_linalg_vendor_elpa_prqs=""
  sd_linalg_vendor_magma_libs=""
  sd_linalg_vendor_magma_prqs=""
  sd_linalg_vendor_plasma_libs=""
  sd_linalg_vendor_plasma_prqs=""

  # Update components according to specified vendor
  case "$1" in

    acml)
      sd_linalg_vendor_provided="blas lapack lapacke scalapack"
      sd_linalg_vendor_blas_libs="-lacml"
      sd_linalg_vendor_blas_prqs="-lacml_mv"
      sd_linalg_vendor_scalapack_prqs="${sd_mpi_libs}"
      ;;

    atlas)
      sd_linalg_vendor_provided="blas"
      sd_linalg_vendor_blas_libs="-lf77blas"
      sd_linalg_vendor_blas_prqs="-lcblas -latlas"
      ;;

    custom)
      sd_linalg_vendor_provided="blas lapack"
      if test "${sd_mpi_enable}" = "yes"; then
        sd_linalg_vendor_provided="${sd_linalg_vendor_provided} scalapack"
      fi
      sd_linalg_vendor_blas_libs="${LINALG_LIBS}"
      ;;

    debian-mpich)
      sd_linalg_vendor_provided="blas lapack scalapack"
      sd_linalg_vendor_blas_libs="-lblas"
      sd_linalg_vendor_lapack_libs="-llapack"
      sd_linalg_vendor_scalapack_libs="-lscalapack-mpich -lblacs-mpich -lblacsCinit-mpich -lblacsF77init-mpich"
      ;;

    debian-openmpi)
      sd_linalg_vendor_provided="blas lapack scalapack"
      sd_linalg_vendor_blas_libs="-lblas"
      sd_linalg_vendor_lapack_libs="-llapack"
      sd_linalg_vendor_scalapack_libs="-lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -lblacsF77init-openmpi"
      ;;

    easybuild)
      sd_linalg_vendor_provided="blas lapack scalapack"
      sd_linalg_vendor_blas_libs="-lopenblas"
      sd_linalg_vendor_scalapack_libs="-lscalapack"
      ;;

    elpa)
      sd_linalg_vendor_provided="elpa"
      sd_linalg_vendor_elpa_libs="-lelpa"
      ;;

    essl)
      sd_linalg_vendor_provided="blas lapack lapacke scalapack"
      sd_linalg_vendor_fcflags="-qessl"
      sd_linalg_vendor_ldflags="-qessl"
      sd_linalg_vendor_blas_libs="-lessl"
      sd_linalg_vendor_scalapack_prqs="${sd_mpi_libs}"
      ;;

    magma)
      sd_linalg_vendor_provided="magma"
      sd_linalg_vendor_magma_libs="-lmagma"
      sd_linalg_vendor_magma_prqs="${sd_gpu_libs}"
      ;;

    mkl)
      if test "${MKLROOT}" = ""; then
        AC_MSG_ERROR([MKLROOT is not set, which means that MKL is not
                  properly configured])
      fi
      sd_linalg_vendor_provided="blas lapack scalapack"
      sd_linalg_vendor_cppflags="-I${MKLROOT}/include"
      sd_linalg_vendor_fcflags="-I${MKLROOT}/include"
      if test "${sd_mpi_enable}" = "yes"; then
        sd_linalg_vendor_ldflags="-mkl=cluster"
      else
        sd_linalg_vendor_ldflags="-mkl"
      fi
      ;;

    netlib)
      sd_linalg_vendor_provided="blas lapack lapacke scalapack"
      sd_linalg_vendor_blas_libs="-lblas"
      sd_linalg_vendor_lapack_libs="-llapack"
      sd_linalg_vendor_lapacke_libs="-llapacke"
      sd_linalg_vendor_scalapack_libs="-lscalapack -lblacs -lblacsCinit -lblacsF77init"
      sd_linalg_vendor_scalapack_prqs="${sd_mpi_libs}"
      ;;

    openblas)
      sd_linalg_vendor_provided="blas"
      sd_linalg_vendor_blas_libs="-lopenblas"
      ;;

    plasma)
      sd_linalg_vendor_provided="plasma"
      sd_linalg_vendor_plasma_libs="-lplasma"
      sd_linalg_vendor_plasma_prqs="-lcoreblas -lcorelapack"
      ;;

    *)
      AC_MSG_ERROR([no library settings for linear algebra flavor '$1'])
      ;;

  esac
]) # _SD_LINALG_SET_VENDOR_FLAGS


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
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


AC_DEFUN([_SD_LINALG_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable linear algebra])
  AC_MSG_RESULT([${sd_linalg_enable}])
  if test "${sd_linalg_enable}" != "no"; then
    AC_MSG_CHECKING([how linear algebra parameters have been set])
    AC_MSG_RESULT([${sd_linalg_init}])
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

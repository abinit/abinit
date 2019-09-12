# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon, Marc Torrent)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for external linear algebra libraries
#


# ABI_LINALG_INIT()
# -----------------
#
# Select a linear algebra flavor according to the available system information.
#
AC_DEFUN([ABI_LINALG_INIT],[
  # Init internal variables
  abi_linalg_chk_gpu=""
  abi_linalg_chk_mpi=""
  abi_linalg_chk_mpiext=""
  abi_linalg_chk_serial=""
  abi_linalg_has_blas="unknown"
  abi_linalg_has_lapack="unknown"
  abi_linalg_has_lapacke="unknown"
  abi_linalg_has_blacs="unknown"
  abi_linalg_has_scalapack="unknown"
  abi_linalg_has_elpa="unknown"
  abi_linalg_has_elpa_2013="unknown"
  abi_linalg_has_elpa_2014="unknown"
  abi_linalg_has_elpa_2015="unknown"
  abi_linalg_has_elpa_2016="unknown"
  abi_linalg_has_plasma="unknown"
  abi_linalg_has_magma="unknown"
  abi_linalg_gpu="unknown"
  abi_linalg_mpi="unknown"
  abi_linalg_provided=""

  # Init build flags
  abi_linalg_cppflags=""
  abi_linalg_cflags=""
  abi_linalg_cxxflags=""
  abi_linalg_fcflags=""
  abi_linalg_ldflags=""
  abi_linalg_libdirs=""
  abi_linalg_libs=""
  case "${abi_linalg_init}" in
    dir)
      abi_linalg_cppflags="-I${abi_linalg_prefix}/include"
      abi_linalg_fcflags="-I${abi_linalg_prefix}/include"
      abi_linalg_libdirs="-L${abi_linalg_prefix}/lib"
      ;;
    env)
      abi_linalg_cppflags="${LINALG_CPPFLAGS}"
      abi_linalg_cflags="${LINALG_CFLAGS}"
      abi_linalg_fcflags="${LINALG_FCFLAGS}"
      abi_linalg_ldflags="${LINALG_LDFLAGS}"
      abi_linalg_libs="${LINALG_LIBS}"
      ;;
  esac

  # Make sure a flavor is selected
  if test "${abi_linalg_flavor}" = ""; then
    abi_linalg_flavor="auto"
  fi

  # Check that the specified flavor is consistent
  if test "${abi_linalg_flavor}" != "auto"; then

    # Reformat flavor
    abi_linalg_iter=`echo "${abi_linalg_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s",[$]1}'`

    # Check serial and parallel flavor unicity
    for tmp_linalg_flavor in ${abi_linalg_iter}; do
      case "${tmp_linalg_flavor}" in
        magma)
          if test "${abi_linalg_chk_gpu}" != ""; then
            AC_MSG_ERROR([only one GPU linear algebra flavor is permitted])
          fi
          abi_linalg_chk_gpu="${tmp_linalg_flavor}"
          ;;
        scalapack|plasma)
          if test "${abi_linalg_chk_mpi}" != ""; then
            AC_MSG_ERROR([only one MPI linear algebra flavor is permitted])
          fi
          abi_linalg_chk_mpi="${tmp_linalg_flavor}"
          ;;
        elpa)
          abi_linalg_chk_mpiext="${tmp_linalg_flavor}"
          ;;
        *)
          if test "${abi_linalg_chk_serial}" != ""; then
            AC_MSG_ERROR([only one serial linear algebra flavor is permitted])
          fi
          abi_linalg_chk_serial="${tmp_linalg_flavor}"
          ;;
      esac
      _ABI_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_flavor}])
    done
    if test "${abi_linalg_chk_serial}" = ""; then
      AC_MSG_ERROR([you must choose a serial linear algebra flavor])
    fi

  fi   # abi_linalg_flavor != auto

  # Init specific flavors
  AC_MSG_CHECKING([for the requested linear algebra support])
  AC_MSG_RESULT([${abi_linalg_flavor}])
  case "${abi_linalg_flavor}" in

    auto)
      # Set generic flavors first
      abi_linalg_chk_serial="netlib"
      if test "${abi_mpi_enable}" = "yes"; then
        abi_linalg_chk_mpi="elpa scalapack"
      fi
      if test "${abi_gpu_enable}" = "yes"; then
        abi_linalg_chk_gpu="magma"
      fi

      # Refine with vendor-specific flavors
      case "${abi_fc_vendor}" in
        gnu)
          abi_linalg_chk_serial="atlas netlib"
          ;;
        intel)
          abi_linalg_chk_serial="mkl atlas netlib"
          ;;
        *)
          abi_linalg_chk_serial="netlib"
          ;;
      esac
      ;;

    none)
      AC_MSG_WARN([bypassing linear algebra tests])
      abi_linalg_has_blas="yes"
      abi_linalg_has_lapack="yes"
      abi_linalg_serial="yes"
      abi_linalg_mpi="no"
      abi_linalg_gpu="no"
      ;;

  esac

  # Display detection sequences
  AC_MSG_CHECKING([for the serial linear algebra detection sequence])
  if test "${abi_linalg_chk_serial}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${abi_linalg_chk_serial}])
  fi
  AC_MSG_CHECKING([for the MPI linear algebra detection sequence])
  if test "${abi_linalg_chk_mpi}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${abi_linalg_chk_mpi}])
  fi
  AC_MSG_CHECKING([for the GPU linear algebra detection sequence])
  if test "${abi_linalg_chk_gpu}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${abi_linalg_chk_gpu}])
  fi

  # Substitute variables needed for the use of the library
  AC_SUBST(abi_linalg_flavor)
  AC_SUBST(abi_linalg_cppflags)
  AC_SUBST(abi_linalg_cflags)
  AC_SUBST(abi_linalg_cxxflags)
  AC_SUBST(abi_linalg_fcflags)
  AC_SUBST(abi_linalg_ldflags)
  AC_SUBST(abi_linalg_libs)
]) # ABI_LINALG_INIT


# ABI_LINALG_DETECT()
# -------------------
#
# Sets all variables needed to handle the optimized linear algebra
# libraries.
#
AC_DEFUN([ABI_LINALG_DETECT], [
  # Explore the system for if the user did not specify anything
  if test "${abi_linalg_cppflags}${abi_linalg_cflags}${abi_linalg_fcflags}${abi_linalg_ldflags}${abi_linalg_libs}" = ""; then
    _ABI_LINALG_EXPLORE
  else
    _ABI_LINALG_CHECK_LIBS
  fi
]) # ABI_LINALG_DETECT


# ABI_LINALG_DETECT_OLD()
# -----------------------
#
# Sets all variables needed to handle the optimized linear algebra
# libraries.
#
AC_DEFUN([ABI_LINALG_DETECT_OLD], [
  # Prepare environment
  ABI_ENV_BACKUP
  LDFLAGS="${FC_LDFLAGS}"
  abi_saved_FCFLAGS="${FCFLAGS}"
  abi_saved_LDFLAGS="${LDFLAGS}"
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_linalg_cppflags}"
  FCFLAGS="${FCFLAGS} ${abi_linalg_fcflags}"
  LIBS="${abi_linalg_libs} ${LIBS}"

  # Look for linear algebra libraries
  if test "${abi_linalg_libs}" != "" -o \
          "${abi_linalg_flavor}" = "custom"; then

    _ABI_LINALG_CHECK_LIBS

  elif test "${abi_linalg_flavor}" != "none"; then

    # BLAS extension?
    _ABI_LINALG_CHECK_BLAS_EXTS()

    # MKL extensions?
    if test "${abi_linalg_chk_serial}" = "mkl"; then
      _ABI_LINALG_CHECK_BLAS_MKL_EXTS()
    fi

    # Look for the selected libraries
    if test "${abi_linalg_fallback}" = "no"; then
      FCFLAGS="${abi_saved_FCFLAGS} ${abi_linalg_fcflags}"
      LDFLAGS="${abi_saved_LDFLAGS} ${abi_linalg_ldflags}"
      _ABI_LINALG_SEARCH_BLAS([${abi_linalg_blas_libs}],
        [${abi_linalg_blas_prqs}])
      _ABI_LINALG_SEARCH_LAPACK([${abi_linalg_lapack_libs}],
        [${abi_linalg_lapack_prqs}])

      # MPI libraries
      case "${abi_linalg_chk_mpi}" in
        scalapack)
          if test "${abi_mpi_enable}" != "yes"; then
            AC_MSG_ERROR([ScaLAPACK support requires MPI])
          fi
          _ABI_LINALG_SEARCH_BLACS([${abi_linalg_blacs_libs}],
            [${abi_linalg_blacs_prqs}])
          _ABI_LINALG_SEARCH_SCALAPACK([${abi_linalg_scalapack_libs}],
            [${abi_linalg_scalapack_prqs}])
          ;;
        plasma)
          if test "${abi_mpi_enable}" != "yes"; then
            AC_MSG_ERROR([PLASMA support requires MPI])
          fi
          if test "${abi_openmp_enable}" != "yes"; then
            AC_MSG_ERROR([PLASMA support requires openMP])
          fi
          if test "${fc_has_iso_c_binding}" != "yes"; then
            AC_MSG_ERROR([PLASMA support requires Fortran 2003 ISO C bindings])
          fi
          _ABI_LINALG_SEARCH_LAPACKE([${abi_linalg_lapacke_libs}],[${abi_linalg_lapacke_prqs}])
          if test "${abi_linalg_has_lapacke}" != ""; then
            _ABI_LINALG_SEARCH_PLASMA([${abi_linalg_plasma_libs}],[${abi_linalg_plasma_prqs}])
          fi
          ;;
        *)
          if test "${abi_linalg_chk_mpi}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_mpi} not implemented])
          fi
          ;;
      esac

      # MPI extension libraries
      case "${abi_linalg_chk_mpiext}" in
        elpa)
          if test "${abi_mpi_enable}" != "yes"; then
            AC_MSG_ERROR([ELPA support requires MPI])
          fi
          if test "${abi_linalg_has_scalapack}" != "yes"; then
            AC_MSG_ERROR([ELPA support requires ScaLAPACK])
          fi
          _ABI_LINALG_SEARCH_ELPA([${abi_linalg_elpa_libs}],
            [${abi_linalg_elpa_prqs}])
          ;;
        *)
          if test "${abi_linalg_chk_mpiext}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_mpiext} not implemented])
          fi
          ;;
      esac

      # GPU extension libraries
      case "${abi_linalg_chk_gpu}" in
        magma)
          if test "${abi_gpu_enable}" != "yes"; then
            AC_MSG_ERROR([MAGMA requires GPU support])
          fi
          _ABI_LINALG_SEARCH_MAGMA([${abi_linalg_magma_libs}],
            [${abi_linalg_magma_prqs}])
          ;;
        *)
          if test "${abi_linalg_chk_gpu}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_gpu} not implemented])
          fi
          ;;
      esac
    fi
  fi

  # Set serial, MPI and GPU status
  if test "${abi_linalg_has_blas}" = "yes" -a \
          "${abi_linalg_has_lapack}" = "yes"; then
    abi_linalg_serial="yes"
    if test "${abi_linalg_has_blacs}" = "yes" -a \
            "${abi_linalg_has_scalapack}" = "yes"; then
      abi_linalg_mpi="yes"
    fi
    if test "${abi_linalg_has_plasma}" = "yes"; then
      abi_linalg_mpi="yes"
    fi
    if test "${abi_linalg_has_magma}" = "yes"; then
      abi_linalg_gpu="yes"
    fi
  fi

  # Transmit serial status to the source code
  AC_MSG_CHECKING([whether we have a serial linear algebra support])
  AC_MSG_RESULT([${abi_linalg_serial}])
  if test "${abi_linalg_serial}" = "yes"; then
    AC_DEFINE([HAVE_LINALG], 1,
      [Define to 1 if you have an optimized linear algebra library.])
    AC_DEFINE([HAVE_LINALG_SERIAL], 1,
      [Define to 1 if you have an optimized serial linear algebra library.])

    case "${abi_linalg_chk_serial}" in
      asl)
        AC_DEFINE([HAVE_LINALG_ASL], 1,
          [Define to 1 if you have the ASL linear algebra library.])
        ;;
      essl)
        AC_DEFINE([HAVE_LINALG_ESSL], 1,
          [Define to 1 if you have the ESSL linear algebra library.])
        ;;
    esac
  else
    abi_linalg_flavor="broken"
    AC_MSG_WARN([falling back to internal linear algebra libraries])
    abi_fallbacks="${abi_fallbacks} linalg"
    abi_linalg_flavor="netlib-fallback"
    abi_dft_linalg_fallback="yes"
  fi

  # Transmit MPI status to the source code
  AC_MSG_CHECKING([whether we have a MPI linear algebra support])
  AC_MSG_RESULT([${abi_linalg_mpi}])
  if test "${abi_linalg_mpi}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MPI], 1,
      [Define to 1 if you have an optimized MPI-parallel linear algebra library.])
    case "${abi_linalg_chk_mpi}" in
      plasma)
        AC_DEFINE([HAVE_LINALG_PLASMA], 1,
          [Define to 1 if you have an optimized PLASMA linear algebra library.])
        ;;
      scalapack)
        AC_DEFINE([HAVE_LINALG_SCALAPACK], 1,
          [Define to 1 if you have an optimized ScaLAPACK linear algebra library.])
        ;;
    esac
    case "${abi_linalg_chk_mpiext}" in
      elpa)
        AC_DEFINE([HAVE_LINALG_ELPA], 1,
          [Define to 1 if you have an optimized ELPA linear algebra library.])
        ;;
    esac
  elif test "${abi_linalg_chk_mpi}" != ""; then
    abi_linalg_flavor="broken"
  fi

  # Transmit GPU status to the source code
  AC_MSG_CHECKING([whether we have a GPU linear algebra support])
  AC_MSG_RESULT([${abi_linalg_gpu}])
  if test "${abi_linalg_gpu}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_GPU], 1,
      [Define to 1 if you have an optimized GPU-compatible linear algebra library.])
    case "${abi_linalg_chk_gpu}" in
      magma)
        AC_DEFINE([HAVE_LINALG_MAGMA], 1,
          [Define to 1 if you have the MAGMA linear algebra library.])
        ;;
    esac
  elif test "${abi_linalg_chk_gpu}" != ""; then
    abi_linalg_flavor="broken"
  fi

  # Restore build environment
  FCFLAGS="${abi_saved_FCFLAGS}"
  LDFLAGS="${abi_saved_LDFLAGS}"
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  # Output final flavor
  AC_MSG_CHECKING([for the actual linear algebra support])
  AC_MSG_RESULT([${abi_linalg_flavor}])
  if test "${abi_linalg_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Linear algebra detection failure])
    AC_MSG_ERROR([the requested ${abi_linalg_flavor} linear algebra flavor is not supported in this environment])
  fi
]) # ABI_LINALG_DETECT_OLD


                    # ------------------------------------ #


#
# Private macros
#


# _ABI_LINALG_CHECK_LIBS()
# ------------------------
#
# Check whether the specified libraries are BLAS and LAPACK
# implementations.
#
AC_DEFUN([_ABI_LINALG_CHECK_LIBS], [
  # Init
  abi_linalg_has_blas="no"
  abi_linalg_has_lapack="no"
  abi_linalg_has_lapacke="no"
  abi_linalg_has_blacs="no"
  abi_linalg_has_scalapack="no"
  abi_linalg_has_elpa="no"
  abi_linalg_has_elpa_2013="no"
  abi_linalg_has_elpa_2014="no"
  abi_linalg_has_elpa_2015="no"
  abi_linalg_has_elpa_2016="no"
  abi_linalg_has_plasma="no"
  abi_linalg_has_magma="no"

  # Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_linalg_cppflags}"
  FCFLAGS="${FCFLAGS} ${abi_linalg_fcflags}"
  LIBS="${abi_linalg_libs} ${abi_gpu_libs} ${abi_mpi_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  # BLAS?
  AC_MSG_CHECKING([for BLAS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
    ]])], [abi_linalg_has_blas="yes"], [abi_linalg_has_blas="no"])
  AC_MSG_RESULT([${abi_linalg_has_blas}])

  # BLAS extensions?
  _ABI_LINALG_CHECK_BLAS_EXTS()

  # MKL BLAS extensions?
  _ABI_LINALG_CHECK_BLAS_MKL_EXTS()

  # LAPACK?
  AC_MSG_CHECKING([for LAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zhpev
    ]])], [abi_linalg_has_lapack="yes"], [abi_linalg_has_lapack="no"])
  AC_MSG_RESULT([${abi_linalg_has_lapack}])

  # LAPACKE?
  AC_MSG_CHECKING([for LAPACKE C API support in specified libraries])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [#include <lapacke.h>],
    [[
      zhpev_;
    ]])],[abi_linalg_has_lapacke="yes"], [abi_linalg_has_lapacke="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${abi_linalg_has_lapacke}])

  # BLACS?
  AC_MSG_CHECKING([for BLACS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call blacs_gridinit
    ]])], [abi_linalg_has_blacs="yes"], [abi_linalg_has_blacs="no"])
  AC_MSG_RESULT([${abi_linalg_has_blacs}])

  # ScaLAPACK?
  AC_MSG_CHECKING([for ScaLAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [abi_linalg_has_scalapack="yes"], [abi_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${abi_linalg_has_scalapack}])

  # ELPA
  AC_MSG_CHECKING([for ELPA support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer,parameter :: n=1, comm=1
      integer :: comm1, comm2, success
      success = get_elpa_communicators(comm, n, n, comm1, comm2)
    ]])], [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
  AC_MSG_RESULT([${abi_linalg_has_elpa}])
  if test "${abi_linalg_has_elpa}" = "yes"; then
    _ABI_LINALG_CHECK_ELPA_2016()
    _ABI_LINALG_CHECK_ELPA_2015()
    _ABI_LINALG_CHECK_ELPA_2014()
    _ABI_LINALG_CHECK_ELPA_2013()
  fi

  # PLASMA?
  AC_MSG_CHECKING([for PLASMA support in specified libraries])
  abi_linalg_chk_plasma="${abi_linalg_has_lapacke}"
  if test "${abi_linalg_chk_plasma}" = "yes"; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use plasma
        call plasma_zhegv
      ]])], [abi_linalg_has_plasma="yes"], [abi_linalg_has_plasma="no"])
  fi
  AC_MSG_RESULT([${abi_linalg_has_plasma}])

  # MAGMA?
  AC_MSG_CHECKING([for MAGMA (version>=1.1.0) support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call magmaf_zhegvd
    ]])], [abi_linalg_has_magma="yes"], [abi_linalg_has_magma="no"])
  AC_MSG_RESULT([${abi_linalg_has_magma}])
  if test "${abi_linalg_has_magma}" = "yes"; then
    _ABI_LINALG_CHECK_MAGMA_15()
  fi

  # Restore environment
  AC_LANG_POP([Fortran])
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_LINALG_CHECK_LIBS


                    # ------------------------------------ #


# _ABI_LINALG_CHECK_BLAS()
# ------------------------
#
# Check whether the build environment provides BLAS.
#
AC_DEFUN([_ABI_LINALG_CHECK_BLAS], [
  abi_linalg_has_blas="unknown"

  AC_MSG_CHECKING([for BLAS support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
    ]])], [abi_linalg_has_blas="yes"], [abi_linalg_has_blas="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_linalg_has_blas}])
]) # _ABI_LINALG_CHECK_BLAS


# _ABI_LINALG_CHECK_BLAS_EXTS()
# -----------------------------
#
# Check whether the specified BLAS implementation provides useful extensions.
#
AC_DEFUN([_ABI_LINALG_CHECK_BLAS_EXTS], [
  # AXPBY family?
  AC_MSG_CHECKING([for AXPBY support in the BLAS libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call saxpby
      call daxpby
      call caxpby
      call zaxpby
    ]])], [abi_linalg_has_axpby="yes"], [abi_linalg_has_axpby="no"])
  AC_MSG_RESULT([${abi_linalg_has_axpby}])

  if test "${abi_linalg_has_axpby}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_AXPBY], 1,
      [Define to 1 if you have an AXPBY BLAS1 extensions.])
  fi

  # gemm3m family
  AC_MSG_CHECKING([for GEMM3M in the BLAS libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call cgemm3m
      call zgemm3m
    ]])], [abi_linalg_has_gemm3m="yes"], [abi_linalg_has_gemm3m="no"])
  AC_MSG_RESULT([${abi_linalg_has_gemm3m}])

  if test "${abi_linalg_has_gemm3m}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_GEMM3M], 1,
      [Define to 1 if you have ?GEMM3M BLAS3 extensions.])
  fi
]) # _ABI_LINALG_CHECK_BLAS_EXTS


# _ABI_LINALG_CHECK_BLAS_MKL_EXTS()
# ---------------------------------
#
# Check whether the specified MKL implementation provides BLAS extensions.
#
AC_DEFUN([_ABI_LINALG_CHECK_BLAS_MKL_EXTS], [
  # mkl_imatcopy family
  AC_MSG_CHECKING([for mkl_imatcopy in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call mkl_simatcopy
      call mkl_dimatcopy
      call mkl_cimatcopy
      call mkl_zimatcopy
    ]])], [abi_linalg_mkl_has_imatcopy="yes"], [abi_linalg_mkl_has_imatcopy="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_imatcopy}])

  if test "${abi_linalg_mkl_has_imatcopy}" = "yes"; then
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
    ]])], [abi_linalg_mkl_has_omatcopy="yes"], [abi_linalg_mkl_has_omatcopy="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_omatcopy}])

  if test "${abi_linalg_mkl_has_omatcopy}" = "yes"; then
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
    ]])], [abi_linalg_mkl_has_omatadd="yes"], [abi_linalg_mkl_has_omatadd="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_omatadd}])

  if test "${abi_linalg_mkl_has_omatadd}" = "yes"; then
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
    ]])], [abi_linalg_mkl_has_threads="yes"], [abi_linalg_mkl_has_threads="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_threads}])

  if test "${abi_linalg_mkl_has_threads}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_THREADS], 1,
      [Define to 1 if you have mkl_*threads extensions.])
  fi
]) # _ABI_LINALG_CHECK_BLAS_MKL_EXTS


# _ABI_LINALG_CHECK_LAPACK()
# --------------------------
#
# Check whether the build environment provides LAPACK.
#
AC_DEFUN([_ABI_LINALG_CHECK_LAPACK], [
  abi_linalg_has_lapack="unknown"

  AC_MSG_CHECKING([for LAPACK support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zhpev
    ]])], [abi_linalg_has_lapack="yes"], [abi_linalg_has_lapack="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_linalg_has_lapack}])
]) # _ABI_LINALG_CHECK_LAPACK


# _ABI_LINALG_CHECK_ELPA_2013()
# -----------------------------
#
# Look for a ELPA 2013 API.
#
AC_DEFUN([_ABI_LINALG_CHECK_ELPA_2013], [
  # Init
  abi_linalg_has_elpa_2013="no"

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
    ]])], [abi_linalg_has_elpa_2013="yes"], [abi_linalg_has_elpa_2013="no"])
  AC_MSG_RESULT([${abi_linalg_has_elpa_2013}])

  if test "${abi_linalg_has_elpa_2013}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2013], 1,
      [Define to 1 if you have ELPA 2013 API support.])
  fi
]) # _ABI_LINALG_CHECK_ELPA_2013


# _ABI_LINALG_CHECK_ELPA_2014()
# -----------------------------
#
# Look for a ELPA 2014 API.
#
AC_DEFUN([_ABI_LINALG_CHECK_ELPA_2014], [
  # Init
  abi_linalg_has_elpa_2014="no"

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
    ]])], [abi_linalg_has_elpa_2014="yes"], [abi_linalg_has_elpa_2014="no"])
  AC_MSG_RESULT([${abi_linalg_has_elpa_2014}])

  if test "${abi_linalg_has_elpa_2014}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2014], 1,
      [Define to 1 if you have ELPA 2014 API support.])
  fi
]) # _ABI_LINALG_CHECK_ELPA_2014


# _ABI_LINALG_CHECK_ELPA_2015()
# -----------------------------
#
# Look for a ELPA 2015 API.
#
AC_DEFUN([_ABI_LINALG_CHECK_ELPA_2015], [
  # Init
  abi_linalg_has_elpa_2015="no"

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
    ]])], [abi_linalg_has_elpa_2015="yes"], [abi_linalg_has_elpa_2015="no"])
  AC_MSG_RESULT([${abi_linalg_has_elpa_2015}])

  if test "${abi_linalg_has_elpa_2015}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2015], 1,
      [Define to 1 if you have ELPA 2015 API support.])
  fi
]) # _ABI_LINALG_CHECK_ELPA_2015


# _ABI_LINALG_CHECK_ELPA_2016()
# -----------------------------
#
# Look for a ELPA 2016 API.
#
AC_DEFUN([_ABI_LINALG_CHECK_ELPA_2016], [
  # Init
  abi_linalg_has_elpa_2016="no"

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
    ]])], [abi_linalg_has_elpa_2016="yes"], [abi_linalg_has_elpa_2016="no"])
  AC_MSG_RESULT([${abi_linalg_has_elpa_2016}])

  if test "${abi_linalg_has_elpa_2016}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA_2016], 1,
      [Define to 1 if you have ELPA 2016 API support.])
  fi
]) # _ABI_LINALG_CHECK_ELPA_2016


# _ABI_LINALG_CHECK_MAGMA_15()
# ----------------------------
#
# Look for MAGMA >=1.5 (requires magma_init and magma_finalize).
#
AC_DEFUN([_ABI_LINALG_CHECK_MAGMA_15], [
  # Init
  abi_linalg_has_magma_15="no"

  AC_MSG_CHECKING([for magma_init/magma_finalize support in specified MAGMA libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call magmaf_init
      call magma_finalize
    ]])], [abi_linalg_has_magma_15="yes"], [abi_linalg_has_magma_15="no"])
  AC_MSG_RESULT([${abi_linalg_has_magma_15}])

  if test "${abi_linalg_has_magma_15}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MAGMA_15], 1,
      [Define to 1 if you have MAGMA >=1.5 API support])
  fi
]) # _ABI_LINALG_CHECK_MAGMA_15


                    # ------------------------------------ #


# _ABI_LINALG_SEARCH_BLACS(BLACS, EXTRA_LIBS)
# -------------------------------------------
#
# Look for a BLACS implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_BLACS], [
  # Init
  abi_linalg_has_blacs="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([blacs_gridinit], $1,
    [abi_linalg_has_blacs="yes"], [abi_linalg_has_blacs="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_blacs}" = "yes"; then
    if test "${ac_cv_search_blacs_gridinit}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_blacs_gridinit} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_BLACS


# _ABI_LINALG_SEARCH_BLAS(BLAS, EXTRA_LIBS)
# -----------------------------------------
#
# Look for a BLAS implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_BLAS], [
  # Init
  abi_linalg_has_blas="unknown"

  # Look for libraries and routines
  AC_MSG_CHECKING([for libraries that may contain BLAS]) 
  AC_MSG_RESULT([$1])
  AC_LANG_PUSH([Fortran])
  AC_SEARCH_LIBS([zgemm], $1,
    [abi_linalg_has_blas="yes"], [abi_linalg_has_blas="no"],
    [$2 ${abi_linalg_libs}])
  AC_LANG_POP([Fortran])
  if test "${abi_linalg_has_blas}" = "yes"; then
    if test "${ac_cv_search_zgemm}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_zgemm} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_BLAS


# _ABI_LINALG_SEARCH_ELPA(ELPA, EXTRA_LIBS)
# -----------------------------------------
#
# Look for a ELPA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_ELPA], [
  # Init
  abi_linalg_has_elpa="no"

  # Look for libraries and routines
  # Has to rewrite AC_SEARCH_LIBS because of mandatory F90 module
  AC_MSG_CHECKING([for the ELPA library])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer, parameter :: n=1, comm=1
      integer :: comm1, comm2, success
      success = get_elpa_communicators(comm, n, n, comm1, comm2)
    ]])], [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
  if test "${abi_linalg_has_elpa}" = "no"; then
    tmp_saved_LIBS="${LIBS}"
    for test_lib in $1; do
      LIBS="-l${test_lib} $2 ${tmp_saved_LIBS}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use elpa1
          integer, parameter :: n=1, comm=1
          integer :: comm1, comm2, success
          success = get_elpa_communicators(comm, n, n, comm1, comm2)
        ]])], [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
      if test "${abi_linalg_has_elpa}" = "yes"; then
        abi_linalg_libs="-l${test_lib} $2 ${abi_linalg_libs}"
        break
      fi
      LIBS="${tmp_saved_LIBS}"
    done
  fi
  AC_MSG_RESULT([${abi_linalg_has_elpa}])

  if test "${abi_linalg_has_elpa}" = "yes"; then
    _ABI_LINALG_CHECK_ELPA_2016()
    _ABI_LINALG_CHECK_ELPA_2015()
    _ABI_LINALG_CHECK_ELPA_2014()
    _ABI_LINALG_CHECK_ELPA_2013()
  fi
]) # _ABI_LINALG_SEARCH_ELPA


# _ABI_LINALG_SEARCH_LAPACK(LAPACK, EXTRA_LIBS)
# ---------------------------------------------
#
# Look for a LAPACK implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_LAPACK], [
  # Init
  abi_linalg_has_lapack="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([zhpev], $1,
    [abi_linalg_has_lapack="yes"], [abi_linalg_has_lapack="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_lapack}" = "yes"; then
    if test "${ac_cv_search_zhpev}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_zhpev} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_LAPACK


# _ABI_LINALG_SEARCH_LAPACKE(LAPACKE, EXTRA_LIBS)
# -----------------------------------------------
#
# Look for a LAPACKE C API implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_LAPACKE], [
  # Init
  abi_linalg_has_lapacke="no"

  # Look for libraries and routines
  AC_MSG_CHECKING([for library containing zhpev_ C API])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [#include <lapacke.h>],
    [[
      zhpev_;
    ]])], [abi_linalg_has_lapacke="yes"], [abi_linalg_has_lapacke="no"])
  if test "${abi_linalg_has_lapacke}" = "no"; then
    tmp_saved_LIBS="${LIBS}"
    for test_lib in $1; do
      LIBS="-l${test_lib} $2 ${tmp_saved_LIBS}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [#include <lapacke.h>],
        [[
          zhpev_;
        ]])], [abi_linalg_has_lapacke="yes"], [abi_linalg_has_lapacke="no"])
      if test "${abi_linalg_has_lapacke}" = "yes"; then
        abi_linalg_libs="-l${test_lib} $2 ${abi_linalg_libs}"
        break
      fi  
    done
    if test "${abi_linalg_has_lapacke}" = "no"; then
      LIBS="${tmp_saved_LIBS}"
    fi
  fi
  AC_LANG_POP([C])
  AC_MSG_RESULT([${abi_linalg_has_lapacke}])
]) # _ABI_LINALG_SEARCH_LAPACKE


# _ABI_LINALG_SEARCH_MAGMA(MAGMA, EXTRA_LIBS)
# -------------------------------------------
#
# Look for a MAGMA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_MAGMA], [
  # Init
  abi_linalg_has_magma="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([magmaf_zheevd], $1,
    [abi_linalg_has_magma="yes"], [abi_linalg_has_magma="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_magma}" = "yes"; then
    if test "${ac_cv_search_magmaf_zheevd}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_magmaf_zheevd} $2 ${abi_linalg_libs}"
    fi
    _ABI_LINALG_CHECK_MAGMA_15()
  fi
]) # _ABI_LINALG_SEARCH_MAGMA


# _ABI_LINALG_SEARCH_PLASMA(PLASMA, EXTRA_LIBS)
# ---------------------------------------------
#
# Look for a PLASMA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_PLASMA], [
  # Init
  abi_linalg_has_plasma="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([plasma_zhegv], $1,
    [abi_linalg_has_plasma="yes"], [abi_linalg_has_plasma="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_plasma}" = "yes"; then
    if test "${ac_cv_search_plasma_zhegv}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_plasma_zhegv} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_PLASMA


# _ABI_LINALG_SEARCH_SCALAPACK(SCALAPACK, EXTRA_LIBS)
# ---------------------------------------------------
#
# Look for a ScaLAPACK implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_SCALAPACK], [
  # Init
  abi_linalg_has_scalapack="no"

  # Look for libraries and routines
  AC_SEARCH_LIBS([pzheevx], $1,
    [abi_linalg_has_scalapack="yes"], [abi_linalg_has_scalapack="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_scalapack}" = "yes"; then
    if test "${ac_cv_search_pzheevx}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_pzheevx} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_SCALAPACK


                    # ------------------------------------ #


# _ABI_LINALG_SET_VENDOR_FLAGS(VENDOR)
# ------------------------------------
#
# Set libraries to look for depending on the specified flavor.
#
AC_DEFUN([_ABI_LINALG_SET_VENDOR_FLAGS], [
  # Reset components
  abi_linalg_vendor_provided=""
  abi_linalg_vendor_fcflags=""
  abi_linalg_vendor_ldflags=""
  abi_linalg_vendor_blas_libs=""
  abi_linalg_vendor_blas_prqs=""
  abi_linalg_vendor_lapack_libs=""
  abi_linalg_vendor_lapack_prqs=""
  abi_linalg_vendor_lapacke_libs=""
  abi_linalg_vendor_lapacke_prqs=""
  abi_linalg_vendor_scalapack_libs=""
  abi_linalg_vendor_scalapack_prqs=""
  abi_linalg_vendor_elpa_libs=""
  abi_linalg_vendor_elpa_prqs=""
  abi_linalg_vendor_magma_libs=""
  abi_linalg_vendor_magma_prqs=""
  abi_linalg_vendor_plasma_libs=""
  abi_linalg_vendor_plasma_prqs=""

  # Update components according to specified vendor
  case "$1" in

    acml)
      abi_linalg_vendor_provided="blas lapack lapacke scalapack"
      abi_linalg_vendor_blas_libs="-lacml"
      abi_linalg_vendor_blas_prqs="-lacml_mv"
      abi_linalg_vendor_scalapack_prqs="${abi_mpi_libs}"
      ;;

    atlas)
      abi_linalg_vendor_provided="blas"
      abi_linalg_vendor_blas_libs="-lf77blas"
      abi_linalg_vendor_blas_prqs="-lcblas -latlas"
      ;;

    debian-mpich)
      abi_linalg_vendor_provided="blas lapack scalapack"
      abi_linalg_vendor_blas_libs="-lblas"
      abi_linalg_vendor_lapack_libs="-llapack"
      abi_linalg_vendor_scalapack_libs="-lscalapack-mpich -lblacs-mpich -lblacsCinit-mpich -lblacsF77init-mpich"
      ;;

    debian-openmpi)
      abi_linalg_vendor_provided="blas lapack scalapack"
      abi_linalg_vendor_blas_libs="-lblas"
      abi_linalg_vendor_lapack_libs="-llapack"
      abi_linalg_vendor_scalapack_libs="-lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -lblacsF77init-openmpi"
      ;;

    easybuild)
      abi_linalg_vendor_provided="blas lapack scalapack"
      abi_linalg_vendor_blas_libs="-lopenblas"
      abi_linalg_vendor_scalapack_libs="-lscalapack"
      ;;

    elpa)
      abi_linalg_vendor_provided="elpa"
      abi_linalg_vendor_elpa_libs="-lelpa"
      ;;

    essl)
      abi_linalg_vendor_provided="blas lapack lapacke scalapack"
      abi_linalg_vendor_fcflags="-qessl"
      abi_linalg_vendor_ldflags="-qessl"
      abi_linalg_vendor_blas_libs="-lessl"
      abi_linalg_vendor_scalapack_prqs="${abi_mpi_libs}"
      ;;

    magma)
      abi_linalg_vendor_provided="magma"
      abi_linalg_vendor_magma_libs="-lmagma"
      abi_linalg_vendor_magma_prqs="${abi_gpu_libs}"
      ;;

    mkl)
      abi_linalg_vendor_provided="blas lapack scalapack"
      if test "${abi_mpi_enable}" = "yes"; then
        abi_linalg_vendor_ldflags="-mkl=cluster"
      else
        abi_linalg_vendor_ldflags="-mkl"
      fi
      ;;

    netlib)
      abi_linalg_vendor_provided="blas lapack lapacke scalapack"
      abi_linalg_vendor_blas_libs="-lblas"
      abi_linalg_vendor_lapack_libs="-llapack"
      abi_linalg_vendor_lapacke_libs="-llapacke"
      abi_linalg_vendor_scalapack_libs="-lscalapack -lblacs -lblacsCinit -lblacsF77init"
      abi_linalg_vendor_scalapack_prqs="${abi_mpi_libs}"
      ;;

    openblas)
      abi_linalg_vendor_provided="blas"
      abi_linalg_vendor_blas_libs="-lopenblas"
      ;;

    plasma)
      abi_linalg_vendor_provided="plasma"
      abi_linalg_vendor_plasma_libs="-lplasma"
      abi_linalg_vendor_plasma_prqs="-lcoreblas -lcorelapack"
      ;;

    *)
      AC_MSG_ERROR([no library settings for linear algebra flavor '$1'])
      ;;

  esac
]) # _ABI_LINALG_SET_VENDOR_FLAGS


                    # ------------------------------------ #


# _ABI_LINALG_EXPLORE()
# ---------------------
#
# Looks for linear algebra components by going through all the selected
# vendor sequences.
#
AC_DEFUN([_ABI_LINALG_EXPLORE], [
  # Prepare environment
  ABI_ENV_BACKUP
  LDFLAGS="${FC_LDFLAGS}"
  abi_saved_CPPFLAGS="${CPPFLAGS}"
  abi_saved_CFLAGS="${CFLAGS}"
  abi_saved_FCFLAGS="${FCFLAGS}"
  abi_saved_LDFLAGS="${LDFLAGS}"
  abi_saved_LIBS="${LIBS}"

  # Look for serial linear algebra support
  for tmp_linalg_vendor in ${abi_linalg_chk_serial}; do

    # Configure vendor libraries
    _ABI_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])
    CPPFLAGS="${abi_saved_CPPFLAGS} ${abi_linalg_vendor_cppflags}"
    CFLAGS="${abi_saved_CFLAGS} ${abi_linalg_vendor_cflags}"
    FCFLAGS="${abi_saved_FCFLAGS} ${abi_linalg_vendor_fcflags}"
    LDFLAGS="${abi_saved_LDFLAGS} ${abi_linalg_vendor_ldflags}"

    # Look for BLAS
    tmp_linalg_blas_proceed=`echo "${abi_linalg_vendor_provided}" | grep "blas"`
    if test "${tmp_linalg_blas_proceed}" != "" -a \
            "${abi_linalg_has_blas}" != "yes"; then
      LIBS="${abi_linalg_vendor_blas_libs} ${abi_linalg_vendor_blas_prqs} ${abi_saved_LIBS}"
      AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for BLAS])
      if test "${abi_linalg_vendor_blas_libs}${abi_linalg_vendor_blas_prqs}" = ""; then
        AC_MSG_RESULT([none required])
      else
        AC_MSG_RESULT([${abi_linalg_vendor_blas_libs} ${abi_linalg_vendor_blas_prqs}])
      fi
      _ABI_LINALG_CHECK_BLAS
      if test "${abi_linalg_has_blas}" = "yes"; then
         #abi_linalg_blas_cppflags="${abi_linalg_vendor_blas_cppflags}"
         #abi_linalg_blas_cflags="${abi_linalg_vendor_blas_cflags}"
         #abi_linalg_blas_fcflags="${abi_linalg_vendor_blas_fcflags}"
         #abi_linalg_blas_ldflags="${abi_linalg_vendor_blas_ldflags}"
         abi_linalg_blas_libs="${abi_linalg_vendor_blas_libs} ${abi_linalg_vendor_blas_prqs}"
         abi_linalg_blas_vendor="${tmp_linalg_vendor}"
         abi_linalg_provided="${abi_linalg_provided} blas"
         _ABI_LINALG_CHECK_BLAS_EXTS
         if test "${tmp_linalg_vendor}" = "mkl"; then
           _ABI_LINALG_CHECK_BLAS_MKL_EXTS
         fi
      fi
    fi

    # Look for LAPACK
    tmp_linalg_lapack_proceed=`echo "${abi_linalg_vendor_provided}" | grep "lapack"`
    if test "${tmp_linalg_lapack_proceed}" != "" -a \
            "${abi_linalg_has_blas}" = "yes" -a \
            "${abi_linalg_has_lapack}" != "yes"; then

      AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for LAPACK])
      if test "${abi_linalg_vendor_lapack_libs}${abi_linalg_vendor_lapack_prqs}" = ""; then
        AC_MSG_RESULT([none required])
      else
       AC_MSG_RESULT([${abi_linalg_vendor_lapack_libs} ${abi_linalg_vendor_lapack_prqs}])
      fi
      LIBS="${abi_linalg_vendor_lapack_libs} ${abi_linalg_vendor_lapack_prqs} ${abi_linalg_blas_libs} ${abi_saved_LIBS}"
      _ABI_LINALG_CHECK_LAPACK
      if test "${abi_linalg_has_lapack}" = "yes"; then
         abi_linalg_lapack_libs="${abi_linalg_vendor_lapack_libs} ${abi_linalg_vendor_lapack_prqs}"
         abi_linalg_lapack_vendor="${tmp_linalg_vendor}"
         abi_linalg_provided="${abi_linalg_provided} lapack"
      fi
    fi

  done

  # FIXME: Look for MPI linear algebra support
  if test "${abi_mpi_enable}" = "yes"; then
    for tmp_linalg_vendor in ${abi_linalg_chk_mpi}; do
      AC_MSG_WARN([MPI linear algebra exploration not implemented!])
    done
  fi

  # FIXME: Look for GPU linear algebra support
  if test "${abi_gpu_enable}" = "yes"; then
    for tmp_linalg_vendor in ${abi_linalg_chk_gpu}; do
      AC_MSG_WARN([GPU linear algebra exploration not implemented!])
    done
  fi

  # Transmit linear algebra libraries found
  if test "${abi_mpi_enable}" = "yes"; then
    abi_linalg_libs="${abi_linalg_libs} ${abi_linalg_vendor_scalapack_libs}"
  fi
  abi_linalg_libs="${abi_linalg_libs} ${abi_linalg_vendor_lapack_libs} ${abi_linalg_vendor_blas_libs}"

  ABI_ENV_RESTORE
  LIBS="${abi_saved_LIBS}"
]) # ABI_LINALG_EXPLORE

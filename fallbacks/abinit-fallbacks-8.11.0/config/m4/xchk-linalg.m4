# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_LINALG()
# ----------------------------
#
# Check whether the specified libraries are BLAS and LAPACK
# implementations.
#
AC_DEFUN([AFB_CHECK_LINALG],[
  dnl Init
  afb_linalg_default_libs="-llapack -lblas"
  afb_linalg_has_incs="unknown"
  afb_linalg_has_libs="unknown"
  afb_linalg_has_blas="unknown"
  afb_linalg_has_blas_axbpy="unknown"
  afb_linalg_has_blas_gemm3m="unknown"
  afb_linalg_has_blas_mkl_imatcopy="unknown"
  afb_linalg_has_blas_mkl_omatcopy="unknown"
  afb_linalg_has_blas_mkl_omatadd="unknown"
  afb_linalg_has_lapack="unknown"
  afb_linalg_has_lapacke="unknown"
  afb_linalg_has_blacs="unknown"
  afb_linalg_has_scalapack="unknown"
  afb_linalg_has_elpa="unknown"
  afb_linalg_has_elpa_new="unknown"
  afb_linalg_has_plasma="unknown"
  afb_linalg_has_magma="unknown"
  afb_linalg_ext_ok="unknown"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${with_linalg_incs}"
  FCFLAGS="${FCFLAGS} ${with_linalg_incs}"
  if test "${afb_linalg_libs}" = ""; then
    AC_MSG_CHECKING([for linear algebra libraries to try])
    LIBS="${afb_linalg_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_linalg_default_libs}])
  else
    LIBS="${afb_linalg_libs} ${LIBS}"
  fi

  dnl BLAS?
  AC_MSG_CHECKING([for BLAS support in specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
    ]])], [afb_linalg_has_blas="yes"], [afb_linalg_has_blas="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${afb_linalg_has_blas}])

  dnl BLAS AXPBY extensions?
  if test "${afb_linalg_has_blas}" = "yes"; then
    AC_MSG_CHECKING([for BLAS AXPBY support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call saxpby
        call daxpby
        call caxpby
        call zaxpby
      ]])], [afb_linalg_has_blas_axpby="yes"], [afb_linalg_has_blas_axpby="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_blas_axpby}])
  fi

  dnl BLAS GEMM3M extensions?
  if test "${afb_linalg_has_blas}" = "yes"; then
    AC_MSG_CHECKING([for GEMM3M support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call cgemm3m
        call zgemm3m
      ]])], [afb_linalg_has_blas_gemm3m="yes"], [afb_linalg_has_blas_gemm3m="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_blas_gemm3m}])
  fi

  dnl BLAS MKL IMATCOPY extensions?
  if test "${afb_linalg_has_blas}" = "yes"; then
    AC_MSG_CHECKING([for MKL IMATCOPY support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call mkl_simatcopy
        call mkl_dimatcopy
        call mkl_cimatcopy
        call mkl_zimatcopy
      ]])], [afb_linalg_has_blas_mkl_imatcopy="yes"], [afb_linalg_has_blas_mkl_imatcopy="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_blas_mkl_imatcopy}])
  fi

  dnl BLAS MKL OMATCOPY extensions?
  if test "${afb_linalg_has_blas}" = "yes"; then
    AC_MSG_CHECKING([for MKL OMATCOPY support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call mkl_somatcopy
        call mkl_domatcopy
        call mkl_comatcopy
        call mkl_zomatcopy
      ]])], [afb_linalg_has_blas_mkl_omatcopy="yes"], [afb_linalg_has_blas_mkl_omatcopy="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_blas_mkl_omatcopy}])
  fi

  dnl BLAS MKL OMATADD extensions?
  if test "${afb_linalg_has_blas}" = "yes"; then
    AC_MSG_CHECKING([for MKL OMATADD support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call mkl_somatadd
        call mkl_domatadd
        call mkl_comatadd
        call mkl_zomatadd
      ]])], [afb_linalg_has_blas_mkl_omatadd="yes"], [afb_linalg_has_blas_mkl_omatadd="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_blas_mkl_omatadd}])
  fi

  dnl LAPACK?
  if test "${afb_linalg_has_blas}" = "yes"; then
    AC_MSG_CHECKING([for LAPACK support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call zhpev
      ]])], [afb_linalg_has_lapack="yes"], [afb_linalg_has_lapack="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_lapack}])
  fi

  dnl LAPACKE?
  if test "${afb_linalg_has_lapack}" = "yes"; then
    AC_MSG_CHECKING([for LAPACKE C API support in specified libraries])
    AC_LANG_PUSH([C])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [
#include <lapacke.h>
      ],[[
        zhpev_;
      ]])],[afb_linalg_has_lapacke="yes"], [afb_linalg_has_lapacke="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${afb_linalg_has_lapacke}])
  fi

  dnl BLACS?
  if test "${afb_linalg_has_lapack}" = "yes"; then
    AC_MSG_CHECKING([for BLACS support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call blacs_gridinit
      ]])], [afb_linalg_has_blacs="yes"], [afb_linalg_has_blacs="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_blacs}])
  fi

  dnl ScaLAPACK?
  if test "${afb_linalg_has_blacs}" = "yes"; then
    AC_MSG_CHECKING([for ScaLAPACK support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call pzheevx
      ]])], [afb_linalg_has_scalapack="yes"], [afb_linalg_has_scalapack="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_scalapack}])
  fi

  dnl ELPA?
  if test "${afb_linalg_has_scalapack}" = "yes"; then
    AC_MSG_CHECKING([for ELPA support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use elpa1
        logical :: success
        real*8 :: a(lda,na),ev(na),q(ldq,na)
        success=solve_evp_real(na,nev,a,lda,ev,q,ldq,nblk,comm_r,comm_c)
      ]])], [afb_linalg_has_elpa="yes"], [afb_linalg_has_elpa="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_elpa}])
  fi

  dnl PLASMA?
  if test "${afb_linalg_has_lapacke}" = "yes"; then
    AC_MSG_CHECKING([for PLASMA support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [[
        use plasma
        call plasma_zhegv
      ]])], [afb_linalg_has_plasma="yes"], [afb_linalg_has_plasma="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_plasma}])
  fi

  dnl MAGMA?
  if test "${afb_linalg_has_lapack}" = "yes"; then
    AC_MSG_CHECKING([for MAGMA (version>=1.1.0) support in specified libraries])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call magmaf_zhegvd
      ]])], [afb_linalg_has_magma="yes"], [afb_linalg_has_magma="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_linalg_has_magma}])
  fi

  dnl Final adjustments
  if test "${afb_linalg_has_lapacke}" = "yes"; then
    afb_linalg_has_incs="yes"
  fi
  if test "${afb_linalg_has_blas}" -a \
          "${afb_linalg_has_lapack}" = "yes"; then
    afb_linalg_has_libs="yes"
  fi
  if test "${afb_linalg_has_libs}" = "yes"; then
    afb_linalg_ext_ok="yes"
  else
    afb_linalg_ext_ok="no"
  fi

  dnl Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_LINALG

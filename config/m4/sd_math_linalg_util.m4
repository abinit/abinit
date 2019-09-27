# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 Yann Pouillon, Marc Torrent
#
# This file is part of the Steredeg software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Support for external linear algebra libraries - Utility macros
#


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
    ]])],
    [sd_linalg_has_blas="yes"; sd_linalg_provided="${sd_linalg_provided} blas"],
    [sd_linalg_has_blas="no"])
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
      [Define to 1 if you have GEMM3M BLAS3 extensions.])
  fi
]) # _SD_LINALG_CHECK_BLAS_EXTS


# _SD_LINALG_CHECK_BLAS_MKL_EXTS()
# --------------------------------
#
# Check whether the specified MKL implementation provides BLAS extensions.
#
AC_DEFUN([_SD_LINALG_CHECK_BLAS_MKL_EXTS], [
  # mkl_imatcopy family
  AC_MSG_CHECKING([for mkl_imatcopy in the specified libraries])
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
      [Define to 1 if you have MKL imatcopy extensions.])
  fi

  # mkl_omatcopy family
  AC_MSG_CHECKING([for mkl_omatcopy in the specified libraries])
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
      [Define to 1 if you have MKL omatcopy extensions.])
  fi

  # mkl_omatadd family
  AC_MSG_CHECKING([for mkl_omatadd in the specified libraries])
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
      [Define to 1 if you have MKL omatadd extensions.])
  fi

  # mkl_threads support functions
  AC_MSG_CHECKING([for mkl_set/get_threads in the specified libraries])
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
    ]])],
    [sd_linalg_has_lapack="yes"; sd_linalg_provided="${sd_linalg_provided} lapack"],
    [sd_linalg_has_lapack="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_lapack}])
]) # _SD_LINALG_CHECK_LAPACK



# _SD_LINALG_CHECK_LAPACKE()
# --------------------------
#
# Check whether the build environment provides LAPACKE.
#
AC_DEFUN([_SD_LINALG_CHECK_LAPACKE], [
  sd_linalg_has_lapacke="unknown"

  AC_MSG_CHECKING([for LAPACKE C API support in the specified libraries])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [#include <lapacke.h>],
    [[
      zhpev_;
    ]])],
    [sd_linalg_has_lapacke="yes"; sd_linalg_provided="${sd_linalg_provided} lapacke"],
    [sd_linalg_has_lapacke="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_linalg_has_lapacke}])
]) # _SD_LINALG_CHECK_LAPACKE


                    # ------------------------------------ #
                    # ------------------------------------ #


# _SD_LINALG_CHECK_BLACS()
# ------------------------
#
# Check whether the build environment provides BLACS.
#
AC_DEFUN([_SD_LINALG_CHECK_BLACS], [
  sd_linalg_has_blacs="unknown"

  AC_MSG_CHECKING([for BLACS support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call blacs_gridinit
    ]])],
    [sd_linalg_has_blacs="yes"; sd_linalg_provided="${sd_linalg_provided} blacs"],
    [sd_linalg_has_blacs="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_blacs}])
]) # _SD_LINALG_CHECK_BLACS


                    # ------------------------------------ #


# _SD_LINALG_CHECK_SCALAPACK()
# ----------------------------
#
# Check whether the build environment provides ScaLAPACK.
#
AC_DEFUN([_SD_LINALG_CHECK_SCALAPACK], [
  sd_linalg_has_scalapack="unknown"

  AC_MSG_CHECKING([for ScaLAPACK support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])],
    [sd_linalg_has_scalapack="yes"; sd_linalg_provided="${sd_linalg_provided} scalapack"],
    [sd_linalg_has_scalapack="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_scalapack}])

  if test "${sd_linalg_has_scalapack}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_SCALAPACK], 1,
      [Define to 1 if you have the ScaLAPACK linear algebra library.])
  fi
]) # _SD_LINALG_CHECK_SCALAPACK


                    # ------------------------------------ #
                    # ------------------------------------ #


# _SD_LINALG_CHECK_PLASMA()
# -------------------------
#
# Look for PLASMA (requires LAPACKE).
#
AC_DEFUN([_SD_LINALG_CHECK_PLASMA], [
  sd_linalg_has_plasma="unknown"

  AC_MSG_CHECKING([for PLASMA support in the specified libraries])
  sd_linalg_chk_plasma="${sd_linalg_has_lapacke}"
  if test "${sd_linalg_chk_plasma}" = "yes"; then
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use plasma
        call plasma_zhegv
      ]])],
      [sd_linalg_has_plasma="yes"; sd_linalg_provided="${sd_linalg_provided} plasma"],
      [sd_linalg_has_plasma="no"])
    AC_LANG_POP([Fortran])
  else
    sd_linalg_has_plasma="no"
  fi
  AC_MSG_RESULT([${sd_linalg_has_plasma}])

  if test "${sd_linalg_has_plasma}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_PLASMA], 1,
      [Define to 1 if you have the PLASMA linear algebra library.])
  fi
]) # _SD_LINALG_CHECK_PLASMA


                    # ------------------------------------ #
                    # ------------------------------------ #


# _SD_LINALG_CHECK_ELPA()
# -----------------------
#
# Look for the ELPA library.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA], [
  sd_linalg_has_elpa="unknown"

  AC_MSG_CHECKING([for ELPA support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer,parameter :: n=1, comm=1
      integer :: comm1, comm2, success
      success = get_elpa_communicators(comm, n, n, comm1, comm2)
    ]])],
    [sd_linalg_has_elpa="yes"; sd_linalg_provided="${sd_linalg_provided} elpa"],
    [sd_linalg_has_elpa="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_elpa}])

  if test "${sd_linalg_has_elpa}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_ELPA], 1,
      [Define to 1 if you have the ELPA linear algebra library.])
    _SD_LINALG_CHECK_ELPA_2017
    if test "${sd_linalg_has_elpa_2017}" != "yes"; then
      _SD_LINALG_CHECK_ELPA_2016
    fi
    if test "${sd_linalg_has_elpa_2016}" != "yes"; then
      _SD_LINALG_CHECK_ELPA_2015
    fi
    if test "${sd_linalg_has_elpa_2015}" != "yes"; then
      _SD_LINALG_CHECK_ELPA_2014
    fi
    if test "${sd_linalg_has_elpa_2014}" != "yes"; then
      _SD_LINALG_CHECK_ELPA_2013
    fi
  fi
]) # _SD_LINALG_CHECK_ELPA


# _SD_LINALG_CHECK_ELPA_2013()
# ----------------------------
#
# Look for a ELPA 2013 API.
#
AC_DEFUN([_SD_LINALG_CHECK_ELPA_2013], [
  sd_linalg_has_elpa_2013="unknown"

  AC_MSG_CHECKING([for ELPA 2013 API support in the specified libraries])
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
  sd_linalg_has_elpa_2014="unknown"

  AC_MSG_CHECKING([for ELPA 2014 API support in the specified libraries])
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
  sd_linalg_has_elpa_2015="unknown"

  AC_MSG_CHECKING([for ELPA 2015 API support in the specified libraries])
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
  sd_linalg_has_elpa_2016="unknown"

  AC_MSG_CHECKING([for ELPA 2016 API support in the specified libraries])
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
  sd_linalg_has_elpa_2017="unknown"

  AC_MSG_CHECKING([for ELPA 2017 API support in the specified libraries])
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
                    # ------------------------------------ #


# _SD_LINALG_CHECK_MAGMA()
# ------------------------
#
# Look for the MAGMA library (requires GPU support).
#
AC_DEFUN([_SD_LINALG_CHECK_MAGMA], [
  sd_linalg_has_magma_15="unknown"

  AC_MSG_CHECKING([for MAGMA (version>=1.1.0) support in the specified libraries])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call magmaf_zhegvd
    ]])],
    [sd_linalg_has_magma="yes"; sd_linalg_provided="${sd_linalg_provided} magma"],
    [sd_linalg_has_magma="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${sd_linalg_has_magma}])
  if test "${sd_linalg_has_magma}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MAGMA], 1,
      [Define to 1 if you have the MAGMA linear algebra library.])
    _SD_LINALG_CHECK_MAGMA_15
  fi
]) # _SD_LINALG_CHECK_MAGMA


# _SD_LINALG_CHECK_MAGMA_15()
# ---------------------------
#
# Look for MAGMA >=1.5 (requires magma_init and magma_finalize).
#
AC_DEFUN([_SD_LINALG_CHECK_MAGMA_15], [
  sd_linalg_has_magma_15="unknown"

  AC_MSG_CHECKING([for MAGMA 1.5 support in the specified libraries])
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

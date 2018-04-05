# -*- Autoconf -*-
#
# Copyright (C) 2005-2018 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for external linear algebra libraries
#


# _ABI_LINALG_CHECK_LIBS()
# ------------------------
#
# Check whether the specified libraries are BLAS and LAPACK
# implementations.
#
AC_DEFUN([_ABI_LINALG_CHECK_LIBS],[
  dnl Init
  abi_linalg_has_blas="no"
  abi_linalg_has_lapack="no"
  abi_linalg_has_lapacke="no"
  abi_linalg_has_blacs="no"
  abi_linalg_has_scalapack="no"
  abi_linalg_has_elpa="no"
  abi_linalg_has_plasma="no"
  abi_linalg_has_magma="no"

  dnl Prepare environment
  tmp_saved_LIBS="${LIBS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  LIBS="${LIBS} ${lib_gpu_libs} ${lib_mpi_libs}"
  FCFLAGS="${FCFLAGS} ${with_linalg_incs}"

  dnl BLAS?
  AC_MSG_CHECKING([for BLAS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call zgemm
    ])], [abi_linalg_has_blas="yes"], [abi_linalg_has_blas="no"])
  AC_MSG_RESULT([${abi_linalg_has_blas}])

  dnl BLAS extension?
  _ABI_LINALG_CHECK_BLAS_EXTS()

  dnl MKL BLAS extensions?
  _ABI_LINALG_CHECK_BLAS_MKL_EXTS()

  dnl LAPACK?
  AC_MSG_CHECKING([for LAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call zhpev
    ])], [abi_linalg_has_lapack="yes"], [abi_linalg_has_lapack="no"])
  AC_MSG_RESULT([${abi_linalg_has_lapack}])

  dnl LAPACKE?
  AC_MSG_CHECKING([for LAPACKE C API support in specified libraries])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([
     #include <lapacke.h>],[zhpev_;
     ])],[abi_linalg_has_lapacke="yes"], [abi_linalg_has_lapacke="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${abi_linalg_has_lapacke}])

  dnl BLACS?
  AC_MSG_CHECKING([for BLACS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call blacs_gridinit
    ])], [abi_linalg_has_blacs="yes"], [abi_linalg_has_blacs="no"])
  AC_MSG_RESULT([${abi_linalg_has_blacs}])

  dnl ScaLAPACK?
  AC_MSG_CHECKING([for ScaLAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call pzheevx
    ])], [abi_linalg_has_scalapack="yes"], [abi_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${abi_linalg_has_scalapack}])

  dnl ELPA?
  abi_linalg_has_elpa="${abi_linalg_has_scalapack}"
  if test "${abi_linalg_has_elpa}" = "yes"; then
    AC_MSG_CHECKING([for ELPA support in specified libraries])
    _ABI_LINALG_TEST_ELPA()
    AC_MSG_RESULT([${abi_linalg_has_elpa}])
    if test "${abi_linalg_has_elpa}" = "yes"; then
      _ABI_LINALG_FIND_ELPA_VERSION()
    fi
  fi

  dnl PLASMA?
  AC_MSG_CHECKING([for PLASMA support in specified libraries])
  abi_linalg_has_plasma="${abi_linalg_has_lapacke}"
  if test "${abi_linalg_has_plasma}" = "yes"; then
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [
        use plasma
      ])], [abi_linalg_has_plasma="yes"], [abi_linalg_has_plasma="no"])
    if test "${abi_linalg_has_plasma}" = "yes"; then
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [
          call plasma_zhegv
        ])], [abi_linalg_has_plasma="yes"], [abi_linalg_has_plasma="no"])
    fi
  fi
  AC_MSG_RESULT([${abi_linalg_has_plasma}])

  dnl MAGMA?
  AC_MSG_CHECKING([for MAGMA (version>=1.1.0) support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call magmaf_zhegvd
    ])], [abi_linalg_has_magma="yes"], [abi_linalg_has_magma="no"])
  AC_MSG_RESULT([${abi_linalg_has_magma}])
  if test "${abi_linalg_has_magma}" = "yes"; then
    _ABI_LINALG_CHECK_MAGMA_15()
  fi

  dnl Restore environment
  LIBS="${tmp_saved_LIBS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
]) # _ABI_LINALG_CHECK_LIBS


# _ABI_LINALG_SEARCH_BLAS(BLAS,EXTRA_LIBS)
# ----------------------------------------
#
# Look for a BLAS implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_BLAS],[
  dnl Init
  abi_linalg_has_blas="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([zgemm],$1,
    [abi_linalg_has_blas="yes"],[abi_linalg_has_blas="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_blas}" = "yes"; then
    if test "${ac_cv_search_zgemm}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_zgemm} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_BLAS


# _ABI_LINALG_CHECK_BLAS_EXTS()
# -----------------------------
#
# Check whether the specified BLAS implementation provides useful extensions.
#

AC_DEFUN([_ABI_LINALG_CHECK_BLAS_EXTS],[

  dnl AXPBY family?
  AC_MSG_CHECKING([for AXPBY support in specified BLAS libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call saxpby
      call daxpby
      call caxpby
      call zaxpby
    ])], [abi_linalg_has_axpby="yes"], [abi_linalg_has_axpby="no"])
  AC_MSG_RESULT([${abi_linalg_has_axpby}])

  if test "${abi_linalg_has_axpby}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_AXPBY],1,[Define to 1 if you have an AXPBY BLAS1 extensions.])
  fi

  dnl gemm3m family
  AC_MSG_CHECKING([for gemm3m in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
     call cgemm3m
     call zgemm3m
    ])], [abi_linalg_has_gemm3m="yes"], [abi_linalg_has_gemm3m="no"])
  AC_MSG_RESULT([${abi_linalg_has_gemm3m}])

  if test "${abi_linalg_has_gemm3m}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_GEMM3M],1,[Define to 1 if you have ?GEMM3M BLAS3 extensions.])
  fi


]) # _ABI_LINALG_CHECK_BLAS_EXTS


# _ABI_LINALG_CHECK_BLAS_MKL_EXTS()
# ---------------------------------
#
# Check whether the specified MKL implementation provides BLAS extensions.
#
AC_DEFUN([_ABI_LINALG_CHECK_BLAS_MKL_EXTS],[

  dnl mkl_imatcopy family
  AC_MSG_CHECKING([for mkl_imatcopy in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call mkl_simatcopy
      call mkl_dimatcopy
      call mkl_cimatcopy
      call mkl_zimatcopy
    ])], [abi_linalg_mkl_has_imatcopy="yes"], [abi_linalg_mkl_has_imatcopy="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_imatcopy}])

  if test "${abi_linalg_mkl_has_imatcopy}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_IMATCOPY],1,[Define to 1 if you have mkl_?imatcopy extensions.])
  fi

  dnl mkl_omatcopy family
  AC_MSG_CHECKING([for mkl_omatcopy in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call mkl_somatcopy
      call mkl_domatcopy
      call mkl_comatcopy
      call mkl_zomatcopy
    ])], [abi_linalg_mkl_has_omatcopy="yes"], [abi_linalg_mkl_has_omatcopy="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_omatcopy}])

  if test "${abi_linalg_mkl_has_omatcopy}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_OMATCOPY],1,[Define to 1 if you have mkl_?omatcopy extensions.])
  fi

  dnl mkl_omatadd family
  AC_MSG_CHECKING([for mkl_omatadd in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call mkl_somatadd
      call mkl_comatadd
      call mkl_domatadd
      call mkl_zomatadd
    ])], [abi_linalg_mkl_has_omatadd="yes"], [abi_linalg_mkl_has_omatadd="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_omatadd}])

  if test "${abi_linalg_mkl_has_omatadd}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_OMATADD],1,[Define to 1 if you have mkl_?omatadd extensions.])
  fi

  dnl mkl_threads support functions
  AC_MSG_CHECKING([for mkl_set/get_threads in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      integer :: a
      a = mkl_get_max_threads()
      call mkl_set_num_threads
    ])], [abi_linalg_mkl_has_threads="yes"], [abi_linalg_mkl_has_threads="no"])
  AC_MSG_RESULT([${abi_linalg_mkl_has_threads}])

  if test "${abi_linalg_mkl_has_threads}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MKL_THREADS],1,[Define to 1 if you have mkl_*threads extensions.])
  fi

]) # _ABI_LINALG_CHECK_BLAS_MKL_EXTS


# _ABI_LINALG_SEARCH_LAPACK(LAPACK,EXTRA_LIBS)
# --------------------------------------------
#
# Look for a LAPACK implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_LAPACK],[
  dnl Init
  abi_linalg_has_lapack="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([zhpev],$1,
    [abi_linalg_has_lapack="yes"],[abi_linalg_has_lapack="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_lapack}" = "yes"; then
    if test "${ac_cv_search_zhpev}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_zhpev} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_LAPACK


# _ABI_LINALG_SEARCH_LAPACKE(LAPACKE,EXTRA_LIBS)
# ----------------------------------------------
#
# Look for a LAPACKE C API implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_LAPACKE],[
  dnl Init
  abi_linalg_has_lapacke="no"

  dnl Look for libraries and routines
  dnl Has to rewrite AC_SEARCH_LIBS because of mandatory C header
  AC_MSG_CHECKING([for library containing zhpev_ C API])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <lapacke.h>],[zhpev_;])],
                 [abi_linalg_has_lapacke="yes"], [])
  if test "${abi_linalg_has_lapacke}" = "no"; then
    tmp_saved_LIBS="${LIBS}"
    for test_lib in $1; do
      LIBS="-l${test_lib} $2 ${tmp_saved_LIBS}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <lapacke.h>],[zhpev_;])],
                     [abi_linalg_has_lapacke="yes"], [])
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


# _ABI_LINALG_SEARCH_BLACS(BLACS,EXTRA_LIBS)
# ------------------------------------------
#
# Look for a BLACS implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_BLACS],[
  dnl Init
  abi_linalg_has_blacs="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([blacs_gridinit],$1,
    [abi_linalg_has_blacs="yes"],[abi_linalg_has_blacs="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_blacs}" = "yes"; then
    if test "${ac_cv_search_blacs_gridinit}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_blacs_gridinit} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_BLACS


# _ABI_LINALG_SEARCH_SCALAPACK(SCALAPACK,EXTRA_LIBS)
# --------------------------------------------------
#
# Look for a ScaLAPACK implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_SCALAPACK],[
  dnl Init
  abi_linalg_has_scalapack="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([pzheevx],$1,
    [abi_linalg_has_scalapack="yes"],[abi_linalg_has_scalapack="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_scalapack}" = "yes"; then
    if test "${ac_cv_search_pzheevx}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_pzheevx} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_SCALAPACK


# _ABI_LINALG_SEARCH_ELPA(ELPA,EXTRA_LIBS)
# ----------------------------------------
#
# Look for a ELPA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_ELPA],[
  dnl Init
  abi_linalg_has_elpa="no"

  dnl Look for libraries and routines
  dnl Has to rewrite AC_SEARCH_LIBS because of mandatory F90 module
  AC_MSG_CHECKING([for ELPA library])
  _ABI_LINALG_TEST_ELPA()
  if test "${abi_linalg_has_elpa}" = "no"; then
    tmp_saved_LIBS="${LIBS}"
    for test_lib in $1; do
      LIBS="-l${test_lib} $2 ${tmp_saved_LIBS}"
      _ABI_LINALG_TEST_ELPA()
      if test "${abi_linalg_has_elpa}" = "yes"; then
        abi_linalg_libs="-l${test_lib} $2 ${abi_linalg_libs}"
        break
      fi
    done
    if test "${abi_linalg_has_elpa}" = "no"; then
      LIBS="${tmp_saved_LIBS}"
    fi
  fi

  AC_MSG_RESULT([${abi_linalg_has_elpa}])

  if test "${abi_linalg_has_elpa}" = "yes"; then
    _ABI_LINALG_FIND_ELPA_VERSION()
  fi
]) # _ABI_LINALG_SEARCH_ELPA


# _ABI_LINALG_FIND_ELPA_VERSION(ELPA,EXTRA_LIBS)
# --------------------------------------------
#
# Try to determine ELPA version.
#
AC_DEFUN([_ABI_LINALG_FIND_ELPA_VERSION],[

  abi_linalg_elpa_version="none"
  AC_MSG_CHECKING([for ELPA library version])

# Check for ELPA 2017
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
    [
    use elpa
    class(elpa_t),pointer :: e
    integer,parameter :: na=1,ncol=1,nrow=1 ; integer :: err
    real*8 :: a(ncol,nrow),ev(na),q(ncol,nrow)
    call e%eigenvectors(a,ev,q,err)
    call e%cholesky(a,err)
    ])], [abi_linalg_has_elpa_2017="yes"], [abi_linalg_has_elpa_2017="no"])
  if test "${abi_linalg_has_elpa_2017}" = "yes"; then
    abi_linalg_elpa_version="2017"
    AC_DEFINE([HAVE_LINALG_ELPA_2017],1,[Define to 1 if you have ELPA 2017 API support])
    AC_DEFINE([HAVE_ELPA_FORTRAN2008],1,[Define to 1 if you have ELPA Fortran 2008 API support])
  else

# Check for ELPA 2016
   AC_LINK_IFELSE([AC_LANG_PROGRAM([],
     [
     use elpa1
     logical :: success1,debug
     integer,parameter :: na=1,lda=1,ldq=1,mcol=1,nev=1,nblk=1,nrow=1
     integer :: comm_g=1,comm_r=1,comm_c=1,success2
     real*8 :: a(lda,nrow),ev(na),q(ldq,nrow)
     complex*16 :: ac(lda,nrow)
     success1=elpa_solve_evp_real_1stage(na,nev,a,lda,ev,q,ldq,nblk,mcol,comm_r,comm_c)
     success1=elpa_cholesky_complex(na,ac,lda,nblk,nrow,comm_r,comm_c,debug)
     success2=elpa_get_communicators(comm_g,na,na,comm_r,comm_c)
     ])], [abi_linalg_has_elpa_2016="yes"], [abi_linalg_has_elpa_2016="no"])
   if test "${abi_linalg_has_elpa_2016}" = "yes"; then
     abi_linalg_elpa_version="2016"
     AC_DEFINE([HAVE_LINALG_ELPA_2016],1,[Define to 1 if you have ELPA 2016 API support])
   else

# Check for ELPA 2015
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [
      use elpa1
      logical :: success1,debug
      integer,parameter :: na=1,lda=1,ldq=1,nev=1,nblk=1
      integer :: comm_g=1,comm_r=1,comm_c=1,success2
      real*8 :: a(lda,na),ev(na),q(ldq,na)
      complex*16 :: ac(lda,na)
      success1=solve_evp_real(na,nev,a,lda,ev,q,ldq,nblk,comm_r,comm_c)
      call cholesky_complex(na,ac,lda,nblk,comm_r,comm_c,debug,success1)
      success2=get_elpa_row_col_comms(comm_g,na,na,comm_r,comm_c)
      ])], [abi_linalg_has_elpa_2015="yes"], [abi_linalg_has_elpa_2015="no"])
    if test "${abi_linalg_has_elpa_2015}" = "yes"; then
      abi_linalg_elpa_version="2015"
       AC_DEFINE([HAVE_LINALG_ELPA_2015],1,[Define to 1 if you have ELPA 2015 API support])
    else

# Check for ELPA 2014
     AC_LINK_IFELSE([AC_LANG_PROGRAM([],
       [
       use elpa1
       logical :: success
       integer,parameter :: na=1,lda=1,ldq=1,nev=1,nblk=1,comm_r=1,comm_c=1
       real*8 :: a(lda,na),ev(na),q(ldq,na)
       complex*16 :: ac(lda,na)
       success=solve_evp_real(na,nev,a,lda,ev,q,ldq,nblk,comm_r,comm_c)
       call invert_trm_complex(na,ac,lda,nblk,comm_r,comm_c,success)
       ])], [abi_linalg_has_elpa_2014="yes"], [abi_linalg_has_elpa_2014="no"])
     if test "${abi_linalg_has_elpa_2014}" = "yes"; then
       abi_linalg_elpa_version="2014"
       AC_DEFINE([HAVE_LINALG_ELPA_2014],1,[Define to 1 if you have ELPA 2014 API support])
     else

# Check for ELPA 2011-2013
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [
        use elpa1
        integer,parameter :: na=1,lda=1,ldq=1,nev=1,nblk=1,comm_r=1,comm_c=1
        real*8 :: a(lda,na),ev(na),q(ldq,na)
        complex*16 :: ac(lda,na)
        call solve_evp_real(na,nev,a,lda,ev,q,ldq,nblk,comm_r,comm_c)
        call invert_trm_complex(na,ac,lda,nblk,comm_r,comm_c)
        ])], [abi_linalg_has_elpa_2013="yes"], [abi_linalg_has_elpa_2013="no"])
      if test "${abi_linalg_has_elpa_2013}" = "yes"; then
        abi_linalg_elpa_version="2011-13"
        AC_DEFINE([HAVE_LINALG_ELPA_2013],1,[Define to 1 if you have ELPA 2013 API support])
      fi
     fi
    fi
   fi
  fi

  AC_MSG_RESULT([${abi_linalg_elpa_version}])
  if test "${abi_linalg_elpa_version}" = "none"; then
    AC_MSG_ERROR([ELPA version was not recognized!])
  fi
]) # _ABI_LINALG_FIND_ELPA_VERSION



# _ABI_LINALG_TEST_ELPA(ELPA,EXTRA_LIBS)
# ----------------------------------------
#
# Test if ELPA can be compiled or not
#
AC_DEFUN([_ABI_LINALG_TEST_ELPA],[
  dnl Init
  abi_linalg_has_elpa="no"

# Check ELPA v2017+
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[use elpa])],
    [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
  if test "${abi_linalg_has_elpa}" = "yes"; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [use elpa
       integer :: nrows=1,err
       class(elpa_t),pointer :: e
       call e%set("local_nrows",nrows,err)
      ])], [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
  fi

# Check ELPA v2017-
  if test "${abi_linalg_has_elpa}" = "no"; then
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[use elpa1])],
      [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
    if test "${abi_linalg_has_elpa}" = "yes"; then
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [use elpa1
         integer,parameter :: n=1,comm=1 ; integer :: comm1,comm2,success
         success=get_elpa_communicators(comm,n,n,comm1,comm2)
        ])], [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
    fi
  fi

# Check ELPA v2013-
  if test "${abi_linalg_has_elpa}" = "no"; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [call elpa_transpose_vectors
      ])], [abi_linalg_has_elpa="yes"], [abi_linalg_has_elpa="no"])
  fi
]) # _ABI_LINALG_TEST_ELPA


# _ABI_LINALG_CHECK_MAGMA_15(MAGMA,EXTRA_LIBS)
# --------------------------------------------
#
# Look for MAGMA >=1.5 (requires magma_init and magma_finalize).
#
AC_DEFUN([_ABI_LINALG_CHECK_MAGMA_15],[
  dnl Init
  abi_linalg_has_magma_15="no"

  AC_MSG_CHECKING([for magma_init/magma_finalize support in specified MAGMA libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call magmaf_init
    ])], [abi_linalg_has_magma_15="yes"], [abi_linalg_has_magma_15="no"])
  AC_MSG_RESULT([${abi_linalg_has_magma_15}])

  if test "${abi_linalg_has_magma_15}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MAGMA_15],1,[Define to 1 if you have MAGMA >=1.5 API support])
  fi
]) # _ABI_LINALG_CHECK_MAGMA_15


# _ABI_LINALG_SEARCH_PLASMA(PLASMA,EXTRA_LIBS)
# --------------------------------------------
#
# Look for a PLASMA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_PLASMA],[
  dnl Init
  abi_linalg_has_plasma="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([plasma_zhegv],$1,
    [abi_linalg_has_plasma="yes"],[abi_linalg_has_plasma="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_plasma}" = "yes"; then
    if test "${ac_cv_search_plasma_zhegv}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_plasma_zhegv} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_PLASMA


# _ABI_LINALG_SEARCH_MAGMA(MAGMA,EXTRA_LIBS)
# ------------------------------------------
#
# Look for a MAGMA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_MAGMA],[
  dnl Init
  abi_linalg_has_magma="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([magmaf_zheevd],$1,
    [abi_linalg_has_magma="yes"],[abi_linalg_has_magma="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_magma}" = "yes"; then
    if test "${ac_cv_search_magmaf_zheevd}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_magmaf_zheevd} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_MAGMA



# ABI_LINALG_DETECT()
# -------------------
#
# Sets all variables needed to handle the optimized linear algebra
# libraries.
#
AC_DEFUN([ABI_LINALG_DETECT],[
  dnl Initial setup
  abi_linalg_chk_gpu=""
  abi_linalg_chk_mpi=""
  abi_linalg_chk_mpiext=""
  abi_linalg_chk_serial=""
  abi_linalg_gpu="no"
  abi_linalg_mpi="no"
  abi_linalg_serial="no"
  abi_linalg_has_blas="no"
  abi_linalg_has_lapack="no"
  abi_linalg_has_lapacke="no"
  abi_linalg_has_blacs="no"
  abi_linalg_has_scalapack="no"
  abi_linalg_has_elpa="no"
  abi_linalg_has_plasma="no"
  abi_linalg_has_magma="no"
  abi_linalg_incs="${with_linalg_incs}"
  abi_linalg_libs="${with_linalg_libs}"
  lib_linalg_flavor="${with_linalg_flavor}"
  lib_linalg_fcflags=""
  lib_linalg_incs=""
  lib_linalg_ldflags=""
  lib_linalg_libs=""

  dnl Prepare environment
  ABI_ENV_BACKUP
  LDFLAGS="${FC_LDFLAGS}"
  abi_saved_FCFLAGS="${FCFLAGS}"
  abi_saved_LDFLAGS="${LDFLAGS}"
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_linalg_incs} ${CPPFLAGS}"
  LIBS="${with_linalg_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Make sure the 'none' flavor is not overriden
  if test "${with_linalg_flavor}" = "none"; then
    if test "${with_linalg_incs}" != "" -o \
            "${with_linalg_libs}" != ""; then
      AC_MSG_ERROR([user-defined linear algebra includes and libraries
                  are not allowed when the flavor is set to 'none'
           solution: use consistent linear algebra options])
    fi
  fi

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested linear algebra support])
  AC_MSG_RESULT([${lib_linalg_flavor}])

  dnl Reformat flavor
  abi_linalg_iter=`echo "${lib_linalg_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s",[$]1}'`

  dnl Check serial and parallel flavor unicity
  for abi_linalg_flavor in ${abi_linalg_iter}; do
    case "${abi_linalg_flavor}" in
      magma)
        if test "${abi_linalg_chk_gpu}" != ""; then
          AC_MSG_ERROR([only one GPU linear algebra flavor is permitted])
        fi
        abi_linalg_chk_gpu="${abi_linalg_flavor}"
        ;;
      scalapack|plasma)
        if test "${abi_linalg_chk_mpi}" != ""; then
          AC_MSG_ERROR([only one MPI linear algebra flavor is permitted])
        fi
        abi_linalg_chk_mpi="${abi_linalg_flavor}"
        ;;
      elpa)
        abi_linalg_chk_mpiext="${abi_linalg_flavor}"
        ;;
      *)
        if test "${abi_linalg_chk_serial}" != ""; then
          AC_MSG_ERROR([only one serial linear algebra flavor is permitted])
        fi
        abi_linalg_chk_serial="${abi_linalg_flavor}"
        ;;
    esac
  done
  if test "${abi_linalg_chk_serial}" = ""; then
    AC_MSG_ERROR([you must choose a serial linear algebra flavor])
  fi

  dnl Check if the user has requested a fallback
  AC_MSG_CHECKING([whether to select a fallback for linear algebra])
  abi_linalg_fallback=`echo "${abi_linalg_chk_serial}" | cut -s -d- -f2`
  if test "${abi_linalg_fallback}" = "fallback"; then
    abi_linalg_fallback="yes"
  else
    abi_linalg_fallback="no"
  fi
  AC_MSG_RESULT([${abi_linalg_fallback}])
  if test "${abi_linalg_fallback}" = "yes"; then
    if test "${enable_fallbacks}" = "no"; then
      AC_MSG_ERROR([fallback requested while fallbacks have been globally disabled])
    fi
    if test "${with_linalg_incs}" != "" -o "${with_linalg_libs}" != ""; then
      AC_MSG_ERROR([you may not specify include or link flags when requesting
                  a fallback (--with-linalg-incs and --with-linalg-libs)])
    fi
  fi

  dnl Look for linear algebra libraries
  if test "${with_linalg_libs}" != "" -o \
          "${lib_linalg_flavor}" = "custom"; then
    _ABI_LINALG_CHECK_LIBS

  elif test "${lib_linalg_flavor}" != "none"; then
    case "${abi_linalg_chk_serial}" in

      acml)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="acml"
        abi_linalg_blas_prqs="-lacml_mv"
        abi_linalg_lapack_libs="acml"
        abi_linalg_lapack_prqs=""
        abi_linalg_lapacke_libs="acml"
        abi_linalg_lapacke_prqs=""
        abi_linalg_blacs_libs="acml"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="acml"
        abi_linalg_scalapack_prqs=""
        ;;

      asl)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="asl"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="asl"
        abi_linalg_lapack_prqs=""
        abi_linalg_lapacke_libs="asl"
        abi_linalg_lapacke_prqs=""
        abi_linalg_blacs_libs="asl"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="asl"
        abi_linalg_scalapack_prqs=""
        ;;

      atlas)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="f77blas"
        abi_linalg_blas_prqs="-lcblas -latlas"
        abi_linalg_lapack_libs="lapack"
        abi_linalg_lapack_prqs=""
        abi_linalg_lapacke_libs=""
        abi_linalg_lapacke_prqs=""
        abi_linalg_blacs_libs="atlas"
        abi_linalg_blacs_prqs="-lscalapack"
        abi_linalg_scalapack_libs="scalapack"
        abi_linalg_scalapack_prqs=""
        ;;

      essl)
        abi_linalg_fcflags="-qessl"
        abi_linalg_ldflags="-qessl"
        abi_linalg_blas_libs="essl"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="essl"
        abi_linalg_lapack_prqs=""
        abi_linalg_lapacke_libs="essl"
        abi_linalg_lapacke_prqs=""
        abi_linalg_blacs_libs="essl"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="essl"
        abi_linalg_scalapack_prqs=""
        ;;

      mkl)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        if test "${abi_cpu_bits}" = "64"; then
          abi_linalg_blas_libs="mkl_intel_lp64"
          if test "${MKLROOT}" != ""; then
            abi_linalg_blas_prqs="-L${MKLROOT}/lib/intel64 -lmkl_sequential -lmkl_core -lpthread -lm"
          else
            abi_linalg_blas_prqs="-lmkl_sequential -lmkl_core -lpthread -lm"
          fi
          abi_linalg_lapack_libs="mkl_intel_lp64"
          abi_linalg_lapack_prqs=""
          abi_linalg_lapacke_libs=""
          abi_linalg_lapacke_prqs=""
          abi_linalg_blacs_libs="mkl_blacs_intelmpi_lp64"
          abi_linalg_blacs_prqs=""
          abi_linalg_scalapack_libs="mkl_scalapack_lp64"
          abi_linalg_scalapack_prqs=""
        else
          abi_linalg_blas_libs="mkl_intel"
          if test "${MKLROOT}" != ""; then
            abi_linalg_blas_prqs="-L${MKLROOT}/lib/ia32 -lmkl_sequential -lmkl_core -lpthread -lm"
          else
            abi_linalg_blas_prqs="-lmkl_sequential -lmkl_core -lpthread -lm"
          fi
          abi_linalg_lapack_libs="mkl_intel"
          abi_linalg_lapack_prqs=""
          abi_linalg_lapacke_libs=""
          abi_linalg_lapacke_prqs=""
          abi_linalg_blacs_libs="mkl_blacs_intelmpi"
          abi_linalg_blacs_prqs=""
          abi_linalg_scalapack_libs="mkl_scalapack_core"
          abi_linalg_scalapack_prqs=""
        fi
        ;;

      netlib|goto)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        if test "${abi_linalg_chk_serial}" = "goto"; then
          abi_linalg_blas_libs="goto"
          abi_linalg_blas_prqs=""
        else
          abi_linalg_blas_libs="blas"
          abi_linalg_blas_prqs=""
        fi
        abi_linalg_lapack_libs="lapack"
        abi_linalg_lapack_prqs=""
        abi_linalg_lapacke_libs="lapacke"
        abi_linalg_lapacke_prqs=""
        abi_linalg_blacs_libs="blacs"
        abi_linalg_blacs_prqs="-lblacsCinit -lblacsF77init"
        abi_linalg_scalapack_libs="scalapack"
        abi_linalg_scalapack_prqs=""
        ;;

      *)
        if test "${abi_linalg_fallback}" = "no"; then
          AC_MSG_ERROR([unknown linear algebra flavor '${lib_linalg_flavor}'])
        fi
        ;;

    esac

    dnl ELPA support is always separate
    abi_linalg_elpa_libs="elpa"
    abi_linalg_elpa_prqs=""

    dnl PLASMA support is always separate
    abi_linalg_plasma_libs="plasma"
    abi_linalg_plasma_prqs="-lcoreblas -lcorelapack"

    dnl MAGMA support is always separate
    abi_linalg_magma_libs="magma"
    abi_linalg_magma_prqs="${lib_gpu_libs}"

    dnl BLAS extension?
    _ABI_LINALG_CHECK_BLAS_EXTS()

    dnl MKL extensions?
    if test "${abi_linalg_chk_serial}" = "mkl"; then
      _ABI_LINALG_CHECK_BLAS_MKL_EXTS()
    fi

    dnl Look for the selected libraries
    if test "${abi_linalg_fallback}" = "no"; then
      FCFLAGS="${abi_saved_FCFLAGS} ${abi_linalg_fcflags}"
      LDFLAGS="${abi_saved_LDFLAGS} ${abi_linalg_ldflags}"
      _ABI_LINALG_SEARCH_BLAS([${abi_linalg_blas_libs}],[${abi_linalg_blas_prqs}])
      _ABI_LINALG_SEARCH_LAPACK([${abi_linalg_lapack_libs}],[${abi_linalg_lapack_prqs}])

      dnl MPI libraries
      case "${abi_linalg_chk_mpi}" in
        scalapack)
          if test "${enable_mpi}" != "yes"; then
            AC_MSG_ERROR([ScaLAPACK support requires MPI])
          fi
          _ABI_LINALG_SEARCH_BLACS([${abi_linalg_blacs_libs}],[${abi_linalg_blacs_prqs}])
          _ABI_LINALG_SEARCH_SCALAPACK([${abi_linalg_scalapack_libs}],[${abi_linalg_scalapack_prqs}])
          ;;
        plasma)
          if test "${enable_mpi}" != "yes"; then
            AC_MSG_ERROR([PLASMA support requires MPI])
          fi
          if test "${enable_openmp}" != "yes"; then
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

      dnl MPI extension libraries
      case "${abi_linalg_chk_mpiext}" in
        elpa)
          if test "${enable_mpi}" != "yes"; then
            AC_MSG_ERROR([ELPA support requires MPI])
          fi
          if test "${abi_linalg_has_scalapack}" != "yes"; then
            AC_MSG_ERROR([ELPA support requires ScaLAPACK])
          fi
          _ABI_LINALG_SEARCH_ELPA([${abi_linalg_elpa_libs}],[${abi_linalg_elpa_prqs}])
          if test "${abi_linalg_has_elpa}" == "no"; then
            AC_MSG_ERROR([ELPA library requested but not found (libelpa.x and/or elpa1.mod missing)])
          fi
          ;;
        *)
          if test "${abi_linalg_chk_mpiext}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_mpiext} not implemented])
          fi
          ;;
      esac

      dnl GPU libraries
      case "${abi_linalg_chk_gpu}" in
        magma)
          if test "${enable_gpu}" != "yes"; then
            AC_MSG_ERROR([MAGMA requires GPU support])
          fi
          _ABI_LINALG_SEARCH_MAGMA([${abi_linalg_magma_libs}],[${abi_linalg_magma_prqs}])
          if test "${abi_linalg_has_magma}" == "no"; then
            AC_MSG_ERROR([Magma library requested but not found])
          fi
          ;;
        *)
          if test "${abi_linalg_chk_gpu}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_gpu} not implemented])
          fi
          ;;
      esac
    fi
  fi

  dnl Set serial, MPI and GPU status
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

  dnl Transmit serial status to the source code
  AC_MSG_CHECKING([whether we have a serial linear algebra support])
  AC_MSG_RESULT([${abi_linalg_serial}])
  if test "${abi_linalg_serial}" = "yes"; then
    AC_DEFINE([HAVE_LINALG],1,[Define to 1 if you have an optimized linear algebra library.])
    AC_DEFINE([HAVE_LINALG_SERIAL],1,[Define to 1 if you have an optimized serial linear algebra library.])

    case "${abi_linalg_chk_serial}" in
      asl)
        AC_DEFINE([HAVE_LINALG_ASL],1,[Define to 1 if you have the ASL linear algebra library.])
        ;;
      essl)
        AC_DEFINE([HAVE_LINALG_ESSL],1,[Define to 1 if you have the ESSL linear algebra library.])
        ;;
    esac

    lib_linalg_fcflags="${abi_linalg_fcflags}"
    lib_linalg_ldflags="${abi_linalg_ldflags}"
    lib_linalg_incs="${abi_linalg_incs}"
    lib_linalg_libs="${abi_linalg_libs}"
  else
    lib_linalg_flavor="broken"
    AC_MSG_WARN([falling back to internal linear algebra libraries])
    abi_fallbacks="${abi_fallbacks} linalg"
    lib_linalg_flavor="netlib-fallback"
    abi_dft_linalg_fallback="yes"
  fi

  dnl Transmit MPI status to the source code
  AC_MSG_CHECKING([whether we have a MPI linear algebra support])
  AC_MSG_RESULT([${abi_linalg_mpi}])
  if test "${abi_linalg_mpi}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MPI],1,[Define to 1 if you have an optimized MPI-parallel linear algebra library.])
    case "${abi_linalg_chk_mpi}" in
      plasma)
        AC_DEFINE([HAVE_LINALG_PLASMA],1,[Define to 1 if you have an optimized PLASMA linear algebra library.])
        ;;
      scalapack)
        AC_DEFINE([HAVE_LINALG_SCALAPACK],1,[Define to 1 if you have an optimized ScaLAPACK linear algebra library.])
        ;;
    esac
    case "${abi_linalg_chk_mpiext}" in
      elpa)
        AC_DEFINE([HAVE_LINALG_ELPA],1,[Define to 1 if you have an optimized ELPA linear algebra library.])
        ;;
    esac
  elif test "${abi_linalg_chk_mpi}" != ""; then
    lib_linalg_flavor="broken"
  fi

  dnl Transmit GPU status to the source code
  AC_MSG_CHECKING([whether we have a GPU linear algebra support])
  AC_MSG_RESULT([${abi_linalg_gpu}])
  if test "${abi_linalg_gpu}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_GPU],1,[Define to 1 if you have an optimized GPU-compatible linear algebra library.])
    case "${abi_linalg_chk_gpu}" in
      magma)
        AC_DEFINE([HAVE_LINALG_MAGMA],1,[Define to 1 if you have the MAGMA linear algebra library.])
        ;;
    esac
  elif test "${abi_linalg_chk_gpu}" != ""; then
    lib_linalg_flavor="broken"
  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  FCFLAGS="${abi_saved_FCFLAGS}"
  LDFLAGS="${abi_saved_LDFLAGS}"
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual linear algebra support])
  AC_MSG_RESULT([${lib_linalg_flavor}])
  if test "${lib_linalg_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    AC_MSG_ERROR([the requested ${with_linalg_flavor} linear algebra flavor is not supported on this architecture])
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_linalg_flavor)
  AC_SUBST(lib_linalg_fcflags)
  AC_SUBST(lib_linalg_ldflags)
  AC_SUBST(lib_linalg_incs)
  AC_SUBST(lib_linalg_libs)
]) # ABI_LINALG_DETECT

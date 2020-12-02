!!****m* ABINIT/m_hide_blas
!! NAME
!!  m_hide_blas
!!
!! FUNCTION
!! This module defines interfaces for overloading BLAS routines.
!! whose goal is twofold. On one hand, using generic interfaces renders
!! the code more readable, especially when the routine can be compiled with
!! different precision type (single-precision or double precision as done for example in the GW code)
!! On the other hand, the generic interfaces defined here introduce a programming
!! layer that can be exploited for interfacing non-standard libraries such as for
!! example CUBLAS routines for GPU computations.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! The convention about names of interfaced routine is: x<name>,
!! where <name> is usually equal to the name of the standard routine
!! without the first character specifying the type and kind.
!! The full list of names is reported below.
!! BLAS procedures interfaced in this module are marked with an asterisk.
!! A complete list of possible overloaded interfaces is provided as guide for future additions.
!!
!! ================
!! ==== BLAS 1 ====
!! ================
!! FUNCTION idamax isamax icamax izamax  ---> XIAMAX(n,dx,incx)
!! * FUNCTION  snrm2  dnrm2 scnrm2 dznmr2  ---> XNRM2(n,x,incx)
!! FUNCTION  sasum  dasum scasum dzasum  ---> XASUM(n,x,incx)
!! * FUNCTION               cdotu  zdotu ---> XDOTU(n,x,incx,y,incy)
!! * FUNCTION               cdotc  zdotc ---> XDOTC(n,x,incx,y,incy)
!! FUNCTION  sdot   ddot                 ---> XDOT(n,x,incx,y,incy)
!! FUNCTION  sdsdot sdot                 ---> XDSDOT(n,x,incx,y,incy)
!! SUBROUTINE saxpy daxpy caxpy  zaxpy   ---> XAXPY(n,ca,cx,incx,cy,incy)
!! * SUBROUTINE scopy dcopy ccopy  zcopy   ---> XCOPY(n,cx,incx,cy,incy)
!! SUBROUTINE srotg drotg crotg  zrotg   ---> XROTG(a,b,c,s)
!! SUBROUTINE srot  drot  csrot  zdrot   ---> XROT(n,x,incx,y,incy,c,s)
!! * SUBROUTINE sscal dscal cscal  zscal
!!                        csscal zdscal  ---> XSCAL(n,a,x,incx)
!! SUBROUTINE sswap dswap cswap  zswap   ---> XSWAP(n,x,incx,y,incy)
!!
!! ================
!! ==== BLAS 2 ====
!! ================
!! SUBROUTINE sgbmv dgbmv cgbmv zgbmv    ---> XGBMV(trans,m,kl,ku,n,alpha,A,lda,X,incx,beta,Y,incy)
!! * SUBROUTINE sgemv dgemv cgemv zgemv    ---> XGEMV(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy)
!! * SUBROUTINE             cgerc zgerc    ---> XGERC(m,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE             cgeru zgeru    ---> XGERU(m,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE             chbmv zhbmv    ---> XHBMV(uplo,n,k,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE             chemv zhemv    ---> XHEMV(uplo,n,alpha,A,lda,X,incx,beta,Y,incy)
!! * SUBROUTINE             cher  zher     ---> XHER(uplo,n,alpha,x,incx,A,lda)
!! SUBROUTINE             cher2 zher2    ---> XHER2(uplo,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE             chpr  zhpr     ---> XHPR(uplo,n,alpha,x,incx,AP)
!! SUBROUTINE             chpr2 zhpr2    ---> XHPR2(uplo,n,alpha,x,incx,y,incy,AP)
!! SUBROUTINE             chpmv zhpmv    ---> XHPMV(uplo,n,alpha,AP,X,incx,beta,Y,incy)
!! SUBROUTINE stbmv dtbmv ctbmv ztbmv    ---> XTBMV(uplo,trans,diag,n,k,A,lda,X,incx)
!! SUBROUTINE stpmv dtpmv ctpmv ztpmv    ---> XTPMV(uplo,trans,diag,n,AP,X,incx)
!! SUBROUTINE strmv dtrmv ctrmv ztrmv    ---> XTRMV(uplo,trans,diag,n,A,lda,X,incx)
!! SUBROUTINE ssymv dsymv                ---> XSYMV(uplo,n,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE ssbmv dsbmv                ---> XSBMV(uplo,n,k,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE sspmv dspmv                ---> XSPMV(uplo,n,alpha,AP,X,incx,beta,Y,incy)
!! SUBROUTINE stbsv dtbsv ctbsv ztbsv    ---> XTBSV(uplo,trans,diag,n,k,A,lda,X,incx)
!! SUBROUTINE stpsv dtpsv ctpsv ztpsv    ---> XTPSV(uplo,trans,diag,n,AP,X,incx)
!! SUBROUTINE strsv dtrsv ctrsv ztrsv    ---> XTRSV(uplo,trans,diag,n,A,lda,X,incx)
!! SUBROUTINE  sger  dger                ---> XGER(m,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE  sspr  dspr                ---> XSPR(uplo,n,alpha,x,incx,AP)
!! SUBROUTINE sspr2 dspr2                ---> XSPR2(uplo,n,alpha,x,incx,y,incy,AP)
!! SUBROUTINE  ssyr  dsyr                ---> XSYR(uplo,n,alpha,x,incx,A,lda)
!! SUBROUTINE ssyr2 dsyr2                ---> XSYR2(uplo,n,alpha,x,incx,y,incy,A,lda)
!!
!! ================
!! ==== BLAS 3 ====
!! ================
!! * SUBROUTINE sgemm dgemm cgemm zgemm      ---> XGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE             chemm zhemm      ---> XHEMM(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE            cher2k zher2k     ---> XHER2K(uplo,trans,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
!! * SUBROUTINE             cherk zherk      ---> XHERK(uplo,trans,n,k,alpha,A,lda,beta,C,ldc)
!! SUBROUTINE ssymm dsymm csymm zsymm      ---> XSYMM(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE ssyr2k dsyr2k csyr2k zsyr2k  ---> XSYR2K(uplo,trans,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE ssyrk dsyrk csyrk zsyrk      ---> XSYRK(uplo,trans,n,k,alpha,A,lda,beta,C,ldc)
!! SUBROUTINE strmm dtrmm ctrmm ztrmm      ---> XTRMM(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb)
!! SUBROUTINE strsm dtrsm ctrsm ztrsm      ---> XTRSM(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb)
!!-------------------------------------------------------------------------------
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_hide_blas

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private

!BLAS1
 public :: xnrm2
 public :: xscal
 public :: xdotu
 public :: xdotc
 public :: xcopy

!BLAS2
 public :: xgemv
 public :: xgerc
 public :: xher

!BLAS3
 public :: xgemm
 public :: xherk

! Helper functions
 public :: blas_cholesky_ortho     ! Cholesky orthogonalization.

 public :: sqmat_itranspose        ! In-place transposition of a square matrix.
 public :: sqmat_otranspose        ! out-of-place transposition of a square matrix.

 public :: sqmat_iconjgtrans       ! in-place conjugate transpose of a square matrix.
 public :: sqmat_oconjgtrans       ! out-of-place conjugate transpose of a square matrix.

!----------------------------------------------------------------------

interface xnrm2
  !
  function snrm2 ( n, x, incx )
    use defs_basis
    real(sp) ::  snrm2
    integer,intent(in) :: incx, n
    real(sp),intent(in) ::  x( * )
  end function snrm2
  !
  function dnrm2 ( n, x, incx )
    use defs_basis
    real(dp) :: dnrm2
    integer,intent(in) :: incx, n
    real(dp),intent(in) ::  x( * )
  end function dnrm2
  !
  function scnrm2( n, x, incx )
    use defs_basis
    real(sp) :: scnrm2
    integer,intent(in) :: incx, n
    complex(spc),intent(in) :: x( * )
  end function scnrm2
  !
  function dznrm2( n, x, incx )
    use defs_basis
    real(dp) :: dznrm2
    integer,intent(in) :: incx, n
    complex(dpc),intent(in) :: x( * )
  end function dznrm2
  !
end interface xnrm2

!-------------------------------------------------------------------------------

interface xscal
  !
  subroutine sscal(n,sa,sx,incx)
    use defs_basis
    implicit none
    integer :: incx
    integer :: n
    real(sp) :: sa
    real(sp) :: sx(*)
  end subroutine sscal
  !
  subroutine  dscal(n,da,dx,incx)
    use defs_basis
    implicit none
    integer :: incx
    integer :: n
    real(dp):: da
    real(dp):: dx(*)
  end subroutine dscal
  !
  subroutine  cscal(n,ca,cx,incx)
    use defs_basis
    implicit none
    integer :: incx
    integer :: n
    complex(spc) :: ca
    complex(spc) :: cx(*)
  end subroutine cscal
  !
  subroutine  zscal(n,za,zx,incx)
    use defs_basis
    implicit none
    integer :: incx
    integer :: n
    complex(dpc) :: za
    complex(dpc) :: zx(*)
  end subroutine zscal
  !
  subroutine  csscal(n,sa,cx,incx)
    use defs_basis
    implicit none
    integer :: incx
    integer :: n
    real(sp) :: sa
    complex(spc) :: cx(*)
  end subroutine csscal
  !
  subroutine  zdscal(n,da,zx,incx)
    use defs_basis
    implicit none
    integer :: incx
    integer :: n
    real(dp) :: da
    complex(dpc) :: zx(*)
  end subroutine zdscal
  !
end interface xscal

!-------------------------------------------------------------------------------

interface xdotu
  !
#ifdef HAVE_LINALG_ZDOTU_BUG
  module procedure cdotu
  module procedure zdotu
#else
  function cdotu(n,cx,incx,cy,incy)
    use defs_basis
    complex(spc) :: cdotu
    complex(spc),intent(in) :: cx(*),cy(*)
    integer,intent(in) :: incx,incy,n
  end function cdotu
  !
  function zdotu(n,zx,incx,zy,incy)
    use defs_basis
    complex(dpc) :: zdotu
    complex(dpc),intent(in) :: zx(*),zy(*)
    integer,intent(in) :: incx,incy,n
  end function zdotu
#endif
  !
end interface xdotu

!-------------------------------------------------------------------------------


! CDOTC, CDOTU, ZDOTC, and ZDOTU are problematic if Mac OS X's Vec lib is used.
! See http://developer.apple.com/hardwaredrivers/ve/errata.html.
! If needed, we replace them with plain Fortran code.

interface xdotc
  !
#ifdef HAVE_LINALG_ZDOTC_BUG
   module procedure cdotc
   module procedure zdotc
#else
  function cdotc(n,cx,incx,cy,incy)
    use defs_basis
    complex(spc) :: cdotc
    complex(spc),intent(in) :: cx(*),cy(*)
    integer,intent(in) :: incx,incy,n
  end function cdotc
  !
  function zdotc(n,zx,incx,zy,incy)
    use defs_basis
    complex(dpc) :: zdotc
    complex(dpc),intent(in) :: zx(*),zy(*)
    integer,intent(in) :: incx,incy,n
  end function zdotc
#endif
  !
end interface xdotc

!-------------------------------------------------------------------------------

interface xcopy
   !module procedure ABI_xcopy
 !
 subroutine scopy(n,sx,incx,sy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   real(sp),intent(in) ::  sx(*)
   real(sp),intent(inout) :: sy(*)
 end subroutine scopy
 !
 subroutine  dcopy(n,dx,incx,dy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   real(dp),intent(in) :: dx(*)
   real(dp),intent(inout) :: dy(*)
 end subroutine dcopy
 !
 subroutine  ccopy(n,cx,incx,cy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   complex(spc),intent(in) :: cx(*)
   complex(spc),intent(inout) :: cy(*)
 end subroutine ccopy
 !
 subroutine  zcopy(n,cx,incx,cy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   complex(dpc),intent(in) :: cx(*)
   complex(dpc),intent(inout) :: cy(*)
 end subroutine zcopy
 !
end interface xcopy

!-------------------------------------------------------------------------------

interface xgemv
  !
  subroutine sgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    real(sp),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    real(sp),intent(in) :: a( lda, * ), x( * )
    real(sp),intent(inout) :: y( * )
  end subroutine sgemv
  !
  subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    real(dp),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    real(dp),intent(in) :: a( lda, * ), x( * )
    real(dp),intent(inout) :: y( * )
  end subroutine dgemv
  !
  subroutine cgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    complex(spc),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    complex(spc),intent(in) :: a( lda, * ), x( * )
    complex(spc),intent(inout) :: y( * )
  end subroutine cgemv
  !
  subroutine zgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    complex(dpc),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    complex(dpc),intent(in) :: a( lda, * ), x( * )
    complex(dpc),intent(inout) :: y( * )
  end subroutine zgemv
  !
end interface xgemv

!-------------------------------------------------------------------------------

interface xgerc
  !
  subroutine cgerc ( m, n, alpha, x, incx, y, incy, a, lda )
    use defs_basis
    complex(spc),intent(in) :: alpha
    integer,intent(in) :: incx, incy, lda, m, n
    complex(spc),intent(inout) :: a( lda, * )
    complex(spc),intent(in) :: x( * ), y( * )
  end subroutine cgerc
  !
  subroutine zgerc ( m, n, alpha, x, incx, y, incy, a, lda )
    use defs_basis
    complex(dpc),intent(in) :: alpha
    integer,intent(in) :: incx, incy, lda, m, n
    complex(dpc),intent(inout) :: a( lda, * )
    complex(dpc),intent(in) :: x( * ), y( * )
  end subroutine zgerc
  !
end interface xgerc

!-------------------------------------------------------------------------------

interface xher
  !
  subroutine cher ( uplo, n, alpha, x, incx, a, lda )
    use defs_basis
    character(len=1),intent(in) :: uplo
    real(spc),intent(in) :: alpha
    integer,intent(in) :: incx, lda, n
    complex(spc),intent(inout) :: a( lda, * )
    complex(spc),intent(in) :: x( * )
  end subroutine cher
  !
  subroutine zher ( uplo, n, alpha, x, incx, a, lda )
    use defs_basis
    character(len=1),intent(in) :: uplo
    real(dpc),intent(in) :: alpha
    integer,intent(in) :: incx, lda, n
    complex(dpc),intent(inout) :: a( lda, * )
    complex(dpc),intent(in) :: x( * )
  end subroutine zher
  !
end interface xher

!-------------------------------------------------------------------------------

interface xgemm
  !
  subroutine sgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    real(sp),intent(in) :: alpha, beta
    real(sp),intent(in) :: a( lda, * ), b( ldb, * )
    real(sp),intent(inout) :: c( ldc, * )
  end subroutine sgemm
  !
  subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    real(dp),intent(in) :: alpha, beta
    real(dp),intent(in) :: a( lda, * ), b( ldb, * )
    real(dp),intent(inout) :: c( ldc, * )
  end subroutine dgemm
  !
  subroutine cgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    complex(spc),intent(in) :: alpha, beta
    complex(spc),intent(in) :: a( lda, * ), b( ldb, * )
    complex(spc),intent(inout) :: c( ldc, * )
  end subroutine cgemm
  !
  subroutine zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    complex(dpc),intent(in) :: alpha, beta
    complex(dpc),intent(in) :: a( lda, * ), b( ldb, * )
    complex(dpc),intent(inout) :: c( ldc, * )
  end subroutine zgemm
  !
end interface xgemm

interface xherk
  !
  subroutine cherk( uplo, trans, n, k, alpha, a, lda, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: uplo
    character(len=1),intent(in) :: trans
    integer,intent(in) :: n,k,lda,ldc
    real(sp),intent(in) :: alpha
    complex(spc),intent(in) :: a( lda, * )
    real(sp),intent(in) :: beta
    complex(spc),intent(inout) :: c( ldc, * )
  end subroutine cherk
  !
  subroutine zherk( uplo, trans, n, k, alpha, a, lda, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: uplo
    character(len=1),intent(in) :: trans
    integer,intent(in) :: n,k,lda,ldc
    real(dp),intent(in) :: alpha
    complex(dpc),intent(in) :: a( lda, * )
    real(dp),intent(in) :: beta
    complex(dpc),intent(inout) :: c( ldc, * )
  end subroutine zherk
  !
end interface xherk

!-------------------------------------------------------------------------------

interface blas_cholesky_ortho
  module procedure blas_cholesky_ortho_spc
  module procedure blas_cholesky_ortho_dpc
end interface blas_cholesky_ortho

interface sqmat_itranspose
  module procedure sqmat_itranspose_sp
  module procedure sqmat_itranspose_dp
  module procedure sqmat_itranspose_spc
  module procedure sqmat_itranspose_dpc
end interface sqmat_itranspose

interface sqmat_otranspose
  module procedure sqmat_otranspose_sp
  module procedure sqmat_otranspose_dp
  module procedure sqmat_otranspose_spc
  module procedure sqmat_otranspose_dpc
end interface sqmat_otranspose

interface sqmat_iconjgtrans
  module procedure sqmat_iconjgtrans_spc
  module procedure sqmat_iconjgtrans_dpc
end interface sqmat_iconjgtrans

interface sqmat_oconjgtrans
  module procedure sqmat_oconjgtrans_spc
  module procedure sqmat_oconjgtrans_dpc
end interface sqmat_oconjgtrans

 real(dp),private,parameter ::  zero_dp = 0._dp
 real(dp),private,parameter ::  one_dp  = 1._dp

 complex(dpc),private,parameter :: czero_dpc = (0._dp,0._dp)
 complex(dpc),private,parameter :: cone_dpc  = (1._dp,0._dp)

CONTAINS  !========================================================================================

! CDOTC, CDOTU, ZDOTC, and ZDOTU are problematic if Mac OS X's Vec lib is used.
! See http://developer.apple.com/hardwaredrivers/ve/errata.html.
! Here we replace them with plain Fortran code.

#ifdef HAVE_LINALG_ZDOTC_BUG
!#warning "Using internal replacement for zdotc. External library cannot be used"
#include "replacements/cdotc.f"
#include "replacements/zdotc.f"
#endif

#ifdef HAVE_LINALG_ZDOTU_BUG
!#warning "Using internal replacement for zdotu. External library cannot be used"
#include "replacements/cdotu.f"
#include "replacements/zdotu.f"
#endif

!----------------------------------------------------------------------

!!***

!!****f* m_hide_blas/blas_cholesky_ortho_spc
!! NAME
!!  blas_cholesky_ortho_spc
!!
!! FUNCTION
!!  Performs the Cholesky orthonormalization of the vectors stored in iomat.
!!
!! INPUTS
!!  vec_size=Size of each vector.
!!  nvec=Number of vectors in iomat
!!
!! OUTPUT
!!  cf_ovlp=Cholesky factorization of the overlap matrix. ovlp = U^H U with U upper triangle matrix returned in cf_ovlp
!!
!! SIDE EFFECTS
!!  iomat(vec_size,nvec)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine blas_cholesky_ortho_spc(vec_size,nvec,iomat,cf_ovlp,use_gemm)

!Arguments ------------------------------------
 integer,intent(in) :: vec_size,nvec
 logical,optional,intent(in) :: use_gemm
 complex(spc),intent(inout) :: iomat(vec_size,nvec)
 complex(spc),intent(out) :: cf_ovlp(nvec,nvec)

!Local variables ------------------------------
!scalars
 integer :: ierr
 logical :: my_usegemm
 character(len=500) :: msg

! *************************************************************************

 ! 1) Calculate overlap_ij =  <phi_i|phi_j>
 ! TODO: use dsyrk
 my_usegemm = .FALSE.; if (PRESENT(use_gemm)) my_usegemm = use_gemm

 if (my_usegemm) then
   call xgemm("Conjugate","Normal",nvec,nvec,vec_size,cone_sp,iomat,vec_size,iomat,vec_size,czero_sp,cf_ovlp,nvec)
 else
   call xherk("U","C", nvec, vec_size, one_sp, iomat, vec_size, zero_sp, cf_ovlp, nvec)
 end if
 !
 ! 2) Cholesky factorization: ovlp = U^H U with U upper triangle matrix.
 call CPOTRF('U',nvec,cf_ovlp,nvec,ierr)
 if (ierr/=0)  then
   write(msg,'(a,i0)')' ZPOTRF returned info= ',ierr
   MSG_ERROR(msg)
 end if
 !
 ! 3) Solve X U = io_mat. On exit iomat is orthonormalized.
 call CTRSM('Right','Upper','Normal','Normal',vec_size,nvec,cone_sp,cf_ovlp,nvec,iomat,vec_size)

end subroutine blas_cholesky_ortho_spc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/blas_cholesky_ortho_dpc
!! NAME
!!  blas_cholesky_ortho_dpc
!!
!! FUNCTION
!!  Performs the Cholesky orthonormalization of the vectors stored in iomat.
!!
!! INPUTS
!!  vec_size=Size of each vector.
!!  nvec=Number of vectors in iomat
!!
!! OUTPUT
!!  cf_ovlp=Cholesky factorization of the overlap matrix. ovlp = U^H U with U upper triangle matrix returned in cf_ovlp
!!
!! SIDE EFFECTS
!!  iomat(vec_size,nvec)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine blas_cholesky_ortho_dpc(vec_size,nvec,iomat,cf_ovlp,use_gemm)

!Arguments ------------------------------------
 integer,intent(in) :: vec_size,nvec
 logical,optional,intent(in) :: use_gemm
 complex(dpc),intent(inout) :: iomat(vec_size,nvec)
 complex(dpc),intent(out) :: cf_ovlp(nvec,nvec)

!Local variables ------------------------------
!scalars
 integer :: ierr
 logical :: my_usegemm
 character(len=500) :: msg

! *************************************************************************

 ! 1) Calculate overlap_ij =  <phi_i|phi_j>
 my_usegemm = .FALSE.; if (PRESENT(use_gemm)) my_usegemm = use_gemm

 if (my_usegemm) then
   call xgemm("Conjugate","Normal",nvec,nvec,vec_size,cone_dpc,iomat,vec_size,iomat,vec_size,czero_dpc,cf_ovlp,nvec)
 else
   call xherk("U","C", nvec, vec_size, one_dp, iomat, vec_size, zero_dp, cf_ovlp, nvec)
 end if
 !
 ! 2) Cholesky factorization: ovlp = U^H U with U upper triangle matrix.
 call ZPOTRF('U',nvec,cf_ovlp,nvec,ierr)
 if (ierr/=0)  then
   write(msg,'(a,i0)')' ZPOTRF returned info= ',ierr
   MSG_ERROR(msg)
 end if
 !
 ! 3) Solve X U = io_mat. On exit io_mat is orthonormalized.
 call ZTRSM('Right','Upper','Normal','Normal',vec_size,nvec,cone_dpc,cf_ovlp,nvec,iomat,vec_size)

end subroutine blas_cholesky_ortho_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_itranspose_sp
!! NAME
!!  sqmat_itranspose_sp
!!
!! FUNCTION
!!  Compute alpha * mat^T in place. target: single precision real matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!
!! SIDE EFFECTS
!!   mat(n,n)=in output, it contains alpha * mat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_itranspose_sp(n,mat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(sp),optional,intent(in) :: alpha
!arrays
 real(sp),intent(inout) :: mat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_IMATCOPY
  if (PRESENT(alpha)) then
    call mkl_simatcopy("Column", "Trans", n, n, alpha, mat, n, n)
  else
    call mkl_simatcopy("Column", "Trans", n, n, one_sp, mat, n, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    mat = alpha * TRANSPOSE(mat)
  else
    mat = TRANSPOSE(mat)
  end if
#endif

end subroutine sqmat_itranspose_sp
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_itranspose_dp
!! NAME
!!  sqmat_itranspose_dp
!!
!! FUNCTION
!!  Compute alpha * mat^T in place. target: double precision real matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!
!! SIDE EFFECTS
!!   mat(n,n)=in output, it contains alpha * mat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_itranspose_dp(n,mat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),optional,intent(in) :: alpha
!arrays
 real(dp),intent(inout) :: mat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_IMATCOPY
  if (PRESENT(alpha)) then
    call mkl_dimatcopy("Column", "Trans", n, n, alpha, mat, n, n)
  else
    call mkl_dimatcopy("Column", "Trans", n, n, one_dp, mat, n, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    mat = alpha * TRANSPOSE(mat)
  else
    mat = TRANSPOSE(mat)
  end if
#endif

end subroutine sqmat_itranspose_dp
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_itranspose_spc
!! NAME
!!  sqmat_itranspose_spc
!!
!! FUNCTION
!!  Compute alpha * mat^T in place. target: single precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!
!! SIDE EFFECTS
!!   mat(n,n)=in output, it contains alpha * mat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_itranspose_spc(n,mat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(spc),optional,intent(in) :: alpha
!arrays
 complex(spc),intent(inout) :: mat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_IMATCOPY
  if (PRESENT(alpha)) then
    call mkl_cimatcopy("Column", "Trans", n, n, alpha, mat, n, n)
  else
    call mkl_cimatcopy("Column", "Trans", n, n, cone_sp, mat, n, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    mat = alpha * TRANSPOSE(mat)
  else
    mat = TRANSPOSE(mat)
  end if
#endif

end subroutine sqmat_itranspose_spc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_itranspose_dpc
!! NAME
!!  sqmat_itranspose_dpc
!!
!! FUNCTION
!!  Compute alpha * mat^T in place. target: double precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!
!! SIDE EFFECTS
!!   mat(n,n)=in output, it contains alpha * mat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_itranspose_dpc(n,mat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),optional,intent(in) :: alpha
!arrays
 complex(dpc),intent(inout) :: mat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_IMATCOPY
  if (PRESENT(alpha)) then
    call mkl_zimatcopy("Column", "Trans", n, n, alpha, mat, n, n)
  else
    call mkl_zimatcopy("Column", "Trans", n, n, cone_dpc, mat, n, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    mat = alpha * TRANSPOSE(mat)
  else
    mat = TRANSPOSE(mat)
  end if
#endif

end subroutine sqmat_itranspose_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_otranspose_sp
!! NAME
!!  sqmat_otranspose_sp
!!
!! FUNCTION
!!  Compute alpha * mat^T out-of-place. target: single precision real matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!  imat(n,n)=Input matrix.
!!
!! OUTPUT
!!  omat(n,n)=contains alpha * imat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_otranspose_sp(n,imat,omat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(sp),optional,intent(in) :: alpha
!arrays
 real(sp),intent(in) :: imat(n,n)
 real(sp),intent(out) :: omat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_OMATCOPY
  if (PRESENT(alpha)) then
    call mkl_somatcopy("Column", "Transpose", n, n, alpha, imat, n, omat, n)
  else
    call mkl_somatcopy("Column", "Transpose", n, n, one_sp, imat, n, omat, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    omat = alpha * TRANSPOSE(imat)
  else
    omat = TRANSPOSE(imat)
  end if
#endif

end subroutine sqmat_otranspose_sp
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_otranspose_dp
!! NAME
!!  sqmat_otranspose_dp
!!
!! FUNCTION
!!  Compute alpha * mat^T out-of-place. target: double precision real matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!  imat(n,n)=Input matrix.
!!
!! OUTPUT
!!  omat(n,n)=contains alpha * imat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_otranspose_dp(n,imat,omat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),optional,intent(in) :: alpha
!arrays
 real(dp),intent(in) :: imat(n,n)
 real(dp),intent(out) :: omat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_OMATCOPY
  if (PRESENT(alpha)) then
    call mkl_domatcopy("Column", "Transpose", n, n, alpha, imat, n, omat, n)
  else
    call mkl_domatcopy("Column", "Transpose", n, n, one_dp, imat, n, omat, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    omat = alpha * TRANSPOSE(imat)
  else
    omat = TRANSPOSE(imat)
  end if
#endif

end subroutine sqmat_otranspose_dp
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_otranspose_spc
!! NAME
!!  sqmat_otranspose_spc
!!
!! FUNCTION
!!  Compute alpha * mat^T out-of-place. target: single precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!  imat(n,n)=Input matrix.
!!
!! OUTPUT
!!  omat(n,n)=contains alpha * imat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_otranspose_spc(n,imat,omat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(spc),optional,intent(in) :: alpha
!arrays
 complex(spc),intent(in) :: imat(n,n)
 complex(spc),intent(out) :: omat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_OMATCOPY
  if (PRESENT(alpha)) then
    call mkl_comatcopy("Column", "Transpose", n, n, alpha, imat, n, omat, n)
  else
    call mkl_comatcopy("Column", "Transpose", n, n, cone_sp, imat, n, omat, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    omat = alpha * TRANSPOSE(imat)
  else
    omat = TRANSPOSE(imat)
  end if
#endif

end subroutine sqmat_otranspose_spc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_otranspose_dpc
!! NAME
!!  sqmat_otranspose_dpc
!!
!! FUNCTION
!!  Compute alpha * mat^T out-of-place. target: double precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!  imat(n,n)=Input matrix.
!!
!! OUTPUT
!!  omat(n,n)=contains alpha * imat^T.
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_otranspose_dpc(n,imat,omat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),optional,intent(in) :: alpha
!arrays
 complex(dpc),intent(in) :: imat(n,n)
 complex(dpc),intent(out) :: omat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_OMATCOPY
  if (PRESENT(alpha)) then
    call mkl_zomatcopy("Column", "Transpose", n, n, alpha, imat, n, omat, n)
  else
    call mkl_zomatcopy("Column", "Transpose", n, n, cone_dpc, imat, n, omat, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    omat = alpha * TRANSPOSE(imat)
  else
    omat = TRANSPOSE(imat)
  end if
#endif

end subroutine sqmat_otranspose_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_iconjgtrans_spc
!! NAME
!!  sqmat_iconjgtrans_spc
!!
!! FUNCTION
!!  Compute alpha * CONJG(mat^T) in place. target: single precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!
!! SIDE EFFECTS
!!   mat(n,n)=in output, it contains alpha * CONJG(mat^T).
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_iconjgtrans_spc(n,mat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(spc),optional,intent(in) :: alpha
!arrays
 complex(spc),intent(inout) :: mat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_IMATCOPY
  if (PRESENT(alpha)) then
    call mkl_cimatcopy("Column", "C", n, n, alpha, mat, n, n)
  else
    call mkl_cimatcopy("Column", "C", n, n, cone_sp, mat, n, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    mat = alpha * TRANSPOSE(CONJG(mat))
  else
    mat = TRANSPOSE(CONJG(mat))
  end if
#endif

end subroutine sqmat_iconjgtrans_spc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_iconjgtrans_dpc
!! NAME
!!  sqmat_iconjgtrans_dpc
!!
!! FUNCTION
!!  Compute alpha * CONJG(mat^T) in place. target: double precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!
!! SIDE EFFECTS
!!   mat(n,n)=in output, it contains alpha * CONJG(mat^T).
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_iconjgtrans_dpc(n, mat, alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),optional,intent(in) :: alpha
!arrays
 complex(dpc),intent(inout) :: mat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_IMATCOPY
  if (PRESENT(alpha)) then
    call mkl_zimatcopy("Column", "C", n, n, alpha, mat, n, n)
  else
    call mkl_zimatcopy("Column", "C", n, n, cone_dpc, mat, n, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    mat = alpha * TRANSPOSE(CONJG(mat))
  else
    mat = TRANSPOSE(CONJG(mat))
  end if
#endif

end subroutine sqmat_iconjgtrans_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_oconjgtrans_spc
!! NAME
!!  sqmat_oconjgtrans_spc
!!
!! FUNCTION
!!  Compute alpha * CONJG(mat^T) out-of-place. target: single precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!  imat(n,n)=Input matrix.
!!
!! OUTPUT
!!  omat(n,n)=contains alpha * CONJG(imat^T).
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_oconjgtrans_spc(n, imat, omat, alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(spc),optional,intent(in) :: alpha
!arrays
 complex(spc),intent(in) :: imat(n,n)
 complex(spc),intent(out) :: omat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_OMATCOPY
  if (PRESENT(alpha)) then
    call mkl_comatcopy("Column", "C", n, n, alpha, imat, n, omat, n)
  else
    call mkl_comatcopy("Column", "C", n, n, cone_sp, imat, n, omat, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    omat = alpha * TRANSPOSE(CONJG(imat))
  else
    omat = TRANSPOSE(CONJG(imat))
  end if
#endif

end subroutine sqmat_oconjgtrans_spc
!!***

!----------------------------------------------------------------------

!!****f* m_hide_blas/sqmat_oconjgtrans_dpc
!! NAME
!!  sqmat_oconjgtrans_dpc
!!
!! FUNCTION
!!  Compute alpha * CONJG(mat^T) out-of-place. target: double precision complex matrix.
!!
!! INPUTS
!!  n=size of the matrix
!!  [alpha]=scalar, set to 1.0 if not present
!!  imat(n,n)=Input matrix.
!!
!! OUTPUT
!!  omat(n,n)=contains alpha * CONJG(imat^T).
!!
!! PARENTS
!!
!! CHILDREN
!!      mkl_zomatcopy
!!
!! SOURCE

subroutine sqmat_oconjgtrans_dpc(n,imat,omat,alpha)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),optional,intent(in) :: alpha
!arrays
 complex(dpc),intent(in) :: imat(n,n)
 complex(dpc),intent(out) :: omat(n,n)

! *************************************************************************

#ifdef HAVE_LINALG_MKL_OMATCOPY
  if (PRESENT(alpha)) then
    call mkl_zomatcopy("Column", "C", n, n, alpha, imat, n, omat, n)
  else
    call mkl_zomatcopy("Column", "C", n, n, cone_dpc, imat, n, omat, n)
  end if
#else
  ! Fallback to Fortran.
  if (PRESENT(alpha)) then
    omat = alpha * TRANSPOSE(CONJG(imat))
  else
    omat = TRANSPOSE(CONJG(imat))
  end if
#endif

end subroutine sqmat_oconjgtrans_dpc
!!***

!----------------------------------------------------------------------

END MODULE m_hide_blas
!!***

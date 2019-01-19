!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_linalg_interfaces
!! NAME
!!  m_linalg_interfaces
!!
!! FUNCTION
!!  Interfaces for the BLAS and LAPACK linear algebra routines.
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2019 ABINIT group (Yann Pouillon)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! WARNING
!!  These routines are used both by real and complex arrays
!!  and are commented (no interface):
!!  - ztrsm, zgemm, zgemv, zhemm, zherk, zher, zgerc
!!  - zcopy, zaxpy, zdscal, zscal, zdotc
!!  - zhpev, zgsev, zheev, zgetrf, zpotrf, zhegv, zhpevx
!!  - zhpgv, zhegst
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_linalg_interfaces

 implicit none

 interface
  subroutine caxpy(n,ca,cx,incx,cy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   complex :: ca
   complex :: cx(*)
   complex :: cy(*)
  end subroutine caxpy
 end interface

 interface
  subroutine ccopy(n,cx,incx,cy,incy)
   implicit none
   integer, intent(in) :: incx, incy, n    !vz_i
   complex, intent(in) :: cx(*)    !vz_i
   complex, intent(inout) :: cy(*)    !vz_i
  end subroutine ccopy
 end interface

 interface
  complex function cdotc(n,cx,incx,cy,incy)
   implicit none
   integer, intent(in) :: incx, incy, n    !vz_i
   complex, intent(in) :: cx(*), cy(*)    !vz_i
  end function cdotc
 end interface

 interface
  complex function cdotu(n,cx,incx,cy,incy)
   implicit none
   integer, intent(in) :: incx, incy, n    !vz_i
   complex, intent(in) :: cx(*), cy(*)    !vz_i
  end function cdotu
 end interface

 interface
  subroutine cgbmv ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: KL
   integer :: KU
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: BETA
   character*1 :: TRANS
   complex :: X( * )
   complex :: Y( * )
  end subroutine cgbmv
 end interface

 interface
  subroutine cgemm ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer,intent(in) :: K,lda,ldb,ldc,m,n    !vz_i
   complex,intent(in) :: A( LDA, * )    !vz_i
   complex,intent(in) :: ALPHA    !vz_i
   complex,intent(in) :: B( LDB, * )    !vz_i
   complex,intent(in) :: BETA    !vz_i
   complex,intent(inout) :: C( LDC, * )    !vz_i
   character*1,intent(in) :: TRANSA    !vz_i
   character*1,intent(in) :: TRANSB    !vz_i
  end subroutine cgemm
 end interface

 interface
  subroutine cgemv ( TRANS, M, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer, intent(in) :: INCX, incy, lda, m, n    !vz_i
   complex, intent(in) :: A( LDA, * )    !vz_i
   complex, intent(in) :: ALPHA, beta    !vz_i
   character*1, intent(in) :: TRANS    !vz_i
   complex, intent(in) :: X( * )    !vz_i
   complex, intent(inout) :: Y( * )    !vz_i
  end subroutine cgemv
 end interface

 interface
  subroutine cgerc ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer, intent(in) :: INCX, incy, lda, m, n    !vz_i
   complex, intent(inout) :: A( LDA, * )    !vz_i
   complex, intent(in) :: ALPHA    !vz_i
   complex, intent(in) :: X( * ), Y( * )    !vz_i
  end subroutine cgerc
 end interface

 interface
  subroutine cgeru ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: X( * )
   complex :: Y( * )
  end subroutine cgeru
 end interface

 interface
  subroutine chbmv ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: K
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: BETA
   character*1 :: UPLO
   complex :: X( * )
   complex :: Y( * )
  end subroutine chbmv
 end interface

 interface
  subroutine chemm ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: B( LDB, * )
   complex :: BETA
   complex :: C( LDC, * )
   character*1 :: SIDE
   character*1 :: UPLO
  end subroutine chemm
 end interface

 interface
  subroutine chemv ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: BETA
   character*1 :: UPLO
   complex :: X( * )
   complex :: Y( * )
  end subroutine chemv
 end interface

 interface
  subroutine cher  ( UPLO, N, ALPHA, X, INCX, A, LDA )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   real :: ALPHA
   character*1 :: UPLO
   complex :: X( * )
  end subroutine cher
 end interface

 interface
  subroutine cher2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   character*1 :: UPLO
   complex :: X( * )
   complex :: Y( * )
  end subroutine cher2
 end interface

 interface
  subroutine cher2k( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: B( LDB, * )
   real :: BETA
   complex :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine cher2k
 end interface

 interface
  subroutine cherk ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
 BETA, C, LDC )
   implicit none
   integer, intent(in) :: K,lda,ldc,n    !vz_i
   complex,intent(in) :: A( LDA, * )    !vz_i
   real,intent(in) :: ALPHA    !vz_i
   real,intent(in) :: BETA    !vz_i
   complex,intent(inout) :: C( LDC, * )    !vz_i
   character*1,intent(in) :: TRANS    !vz_i
   character*1,intent(in) :: UPLO    !vz_i
  end subroutine cherk
 end interface

 interface
  subroutine chpmv ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   complex :: ALPHA
   complex :: AP( * )
   complex :: BETA
   character*1 :: UPLO
   complex :: X( * )
   complex :: Y( * )
  end subroutine chpmv
 end interface

 interface
  subroutine chpr  ( UPLO, N, ALPHA, X, INCX, AP )
   implicit none
   integer :: INCX
   integer :: N
   real :: ALPHA
   complex :: AP( * )
   character*1 :: UPLO
   complex :: X( * )
  end subroutine chpr
 end interface

 interface
  subroutine chpr2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   complex :: ALPHA
   complex :: AP( * )
   character*1 :: UPLO
   complex :: X( * )
   complex :: Y( * )
  end subroutine chpr2
 end interface

 interface
  subroutine crotg(ca,cb,c,s)
   implicit none
   real :: c
   complex :: ca
   complex :: cb
   complex :: s
  end subroutine crotg
 end interface

 interface
  subroutine  cscal(n,ca,cx,incx)
   implicit none
   integer :: incx
   integer :: n
   complex :: ca
   complex :: cx(*)
  end subroutine cscal
 end interface

 interface
  subroutine  csrot (n,cx,incx,cy,incy,c,s)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   real :: c
   real :: s
   complex :: cx(1)
   complex :: cy(1)
  end subroutine csrot
 end interface

 interface
  subroutine  csscal(n,sa,cx,incx)
   implicit none
   integer :: incx
   integer :: n
   real :: sa
   complex :: cx(*)
  end subroutine csscal
 end interface

 interface
  subroutine  cswap (n,cx,incx,cy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   complex :: cx(*)
   complex :: cy(*)
  end subroutine cswap
 end interface

 interface
  subroutine csymm ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: B( LDB, * )
   complex :: BETA
   complex :: C( LDC, * )
   character*1 :: SIDE
   character*1 :: UPLO
  end subroutine csymm
 end interface

 interface
  subroutine csyr2k( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: B( LDB, * )
   complex :: BETA
   complex :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine csyr2k
 end interface

 interface
  subroutine csyrk ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: BETA
   complex :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine csyrk
 end interface

 interface
  subroutine ctbmv ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: K
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex :: X( * )
  end subroutine ctbmv
 end interface

 interface
  subroutine ctbsv ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: K
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex :: X( * )
  end subroutine ctbsv
 end interface

 interface
  subroutine ctpmv ( UPLO, TRANS, DIAG, N, AP, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex :: AP( * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex :: X( * )
  end subroutine ctpmv
 end interface

 interface
  subroutine ctpsv ( UPLO, TRANS, DIAG, N, AP, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex :: AP( * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex :: X( * )
  end subroutine ctpsv
 end interface

 interface
  subroutine ctrmm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine ctrmm
 end interface

 interface
  subroutine ctrmv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex :: X( * )
  end subroutine ctrmv
 end interface

 interface
  subroutine ctrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine ctrsm
 end interface

 interface
  subroutine ctrsv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex :: X( * )
  end subroutine ctrsv
 end interface

 interface
  double precision function dasum(n,dx,incx)
   implicit none
   integer :: incx
   integer :: n
   double precision :: dx(*)
  end function dasum
 end interface

 interface
  subroutine daxpy(n,da,dx,incx,dy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   double precision :: da
   double precision :: dx(*)
   double precision :: dy(*)
  end subroutine daxpy
 end interface

 interface
  double precision function dcabs1(z)
   implicit none
   double complex :: z
  end function dcabs1
 end interface

 !interface
 ! subroutine  dcopy(n,dx,incx,dy,incy)
 !  implicit none
 !  integer :: incx
 !  integer :: incy
 !  integer :: n
 !  double precision :: dx(*)
 !  double precision :: dy(*)
 ! end subroutine dcopy
 !end interface

 interface
  double precision function ddot(n,dx,incx,dy,incy)
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   double precision,intent(in) :: dx(*)
   double precision,intent(in) :: dy(*)
  end function ddot
 end interface

 interface
  subroutine dgbmv ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: KL
   integer :: KU
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: BETA
   character*1 :: TRANS
   double precision :: X( * )
   double precision :: Y( * )
  end subroutine dgbmv
 end interface

 interface
  subroutine dgemm ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer,intent(in) :: K,lda,ldb,ldc,m,n    !vz_i
   double precision, intent(in) :: A( LDA, * )    !vz_i
   double precision,intent(in) :: ALPHA    !vz_i
   double precision,intent(in) :: B( LDB, * )    !vz_i
   double precision,intent(in) :: BETA    !vz_i
   double precision,intent(inout) :: C( LDC, * )    !vz_i
   character*1,intent(in) :: TRANSA    !vz_i
   character*1,intent(in) :: TRANSB    !vz_i
  end subroutine dgemm
 end interface

 interface
  subroutine dgemv ( TRANS, M, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer, intent(in) :: INCX,incy,lda,m,n    !vz_i
   double precision, intent(in) :: A( LDA, * )    !vz_i
   double precision, intent(in) :: ALPHA    !vz_i
   double precision, intent(in) :: BETA    !vz_i
   character*1, intent(in) :: TRANS    !vz_i
   double precision, intent(in) :: X( * )    !vz_i
   double precision, intent(inout) :: Y( * )    !vz_i
  end subroutine dgemv
 end interface

 interface
  subroutine dger  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: X( * )
   double precision :: Y( * )
  end subroutine dger
 end interface

 interface
  double precision function dnrm2 ( N, X, INCX )
   implicit none
   integer, intent(in) :: INCX, n    !vz_i
   double precision,intent(in) :: X( * )    !vz_i
  end function dnrm2
 end interface

 interface
  subroutine  drot (n,dx,incx,dy,incy,c,s)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   double precision :: c
   double precision :: s
   double precision :: dx(*)
   double precision :: dy(*)
  end subroutine drot
 end interface

 interface
  subroutine drotg(da,db,c,s)
   implicit none
   double precision :: c
   double precision :: da
   double precision :: db
   double precision :: s
  end subroutine drotg
 end interface

 interface
  subroutine drotm (N,DX,INCX,DY,INCY,DPARAM)
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   double precision :: DPARAM(5)
   double precision :: DX(1)
   double precision :: DY(1)
  end subroutine drotm
 end interface

 interface
  subroutine drotmg (DD1,DD2,DX1,DY1,DPARAM)
   implicit none
   double precision :: DD1
   double precision :: DD2
   double precision :: DPARAM(5)
   double precision :: DX1
   double precision :: DY1
  end subroutine drotmg
 end interface

 interface
  subroutine dsbmv ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: K
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: BETA
   character*1 :: UPLO
   double precision :: X( * )
   double precision :: Y( * )
  end subroutine dsbmv
 end interface

 interface
  subroutine  dscal(n,da,dx,incx)
   implicit none
   integer :: incx
   integer :: n
   double precision :: da
   double precision :: dx(*)
  end subroutine dscal
 end interface

 interface
  double precision function dsdot (N, SX, INCX, SY, INCY)
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   real :: SX(*)
   real :: SY(*)
  end function dsdot
 end interface

 interface
  subroutine  dswap (n,dx,incx,dy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   double precision :: dx(*)
   double precision :: dy(*)
  end subroutine dswap
 end interface

 interface
  subroutine dsymm ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: B( LDB, * )
   double precision :: BETA
   double precision :: C( LDC, * )
   character*1 :: SIDE
   character*1 :: UPLO
  end subroutine dsymm
 end interface

 interface
  subroutine dsymv ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: BETA
   character*1 :: UPLO
   double precision :: X( * )
   double precision :: Y( * )
  end subroutine dsymv
 end interface

 interface
  subroutine dtrmm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine dtrmm
 end interface

 interface
  subroutine dtrmv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   double precision :: X( * )
  end subroutine dtrmv
 end interface

 interface
  subroutine dtrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: ALPHA
   double precision :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine dtrsm
 end interface

 interface
  subroutine dtrsv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   double precision :: X( * )
  end subroutine dtrsv
 end interface

 interface
  double precision function dzasum(n,zx,incx)
   implicit none
   integer :: incx
   integer :: n
   double complex :: zx(*)
  end function dzasum
 end interface

 interface
  double precision function dznrm2( N, X, INCX )
   implicit none
   integer, intent(in) :: INCX, n    !vz_i
   complex*16,intent(in) :: X( * )    !vz_i
  end function dznrm2
 end interface

 interface
  integer function icamax(n,cx,incx)
   implicit none
   integer :: incx
   integer :: n
   complex :: cx(*)
  end function icamax
 end interface

 interface
  integer function idamax(n,dx,incx)
   implicit none
   integer :: incx
   integer :: n
   double precision :: dx(*)
  end function idamax
 end interface

 interface
  integer function isamax(n,sx,incx)
   implicit none
   integer :: incx
   integer :: n
   real :: sx(*)
  end function isamax
 end interface

 interface
  integer function izamax(n,zx,incx)
   implicit none
   integer :: incx
   integer :: n
   double complex :: zx(*)
  end function izamax
 end interface

 interface
  real function sasum(n,sx,incx)
   implicit none
   integer :: incx
   integer :: n
   real :: sx(*)
  end function sasum
 end interface

 interface
  subroutine saxpy(n,sa,sx,incx,sy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   real :: sa
   real :: sx(*)
   real :: sy(*)
  end subroutine saxpy
 end interface

 interface
  real function scasum(n,cx,incx)
   implicit none
   integer :: incx
   integer :: n
   complex :: cx(*)
  end function scasum
 end interface

 interface
  real function scnrm2( N, X, INCX )
   implicit none
   integer, intent(in) :: INCX, n    !vz_i
   complex, intent(in) :: X( * )    !vz_i
  end function scnrm2
 end interface

 interface
  subroutine scopy(n,sx,incx,sy,incy)
   implicit none
   integer, intent(in) :: incx,incy,n    !vz_i
   real, intent(in) :: sx(*)    !vz_i
   real, intent(inout) :: sy(*)    !vz_i
  end subroutine scopy
 end interface

 interface
  real function sdot(n,sx,incx,sy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   real :: sx(*)
   real :: sy(*)
  end function sdot
 end interface

 interface
  real function sdsdot (N, SB, SX, INCX, SY, INCY)
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   real :: SB
   real :: SX(*)
   real :: SY(*)
  end function sdsdot
 end interface

 interface
  subroutine sgbmv ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: KL
   integer :: KU
   integer :: LDA
   integer :: M
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: BETA
   character*1 :: TRANS
   real :: X( * )
   real :: Y( * )
  end subroutine sgbmv
 end interface

 interface
  subroutine sgemm ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer, intent(in) :: K,lda,ldb,ldc,m,n    !vz_i
   real,intent(in) :: A( LDA, * )    !vz_i
   real,intent(in) :: ALPHA    !vz_i
   real,intent(in) :: B( LDB, * )    !vz_i
   real,intent(in) :: BETA    !vz_i
   real,intent(inout) :: C( LDC, * )    !vz_i
   character*1,intent(in) :: TRANSA    !vz_i
   character*1,intent(in) :: TRANSB    !vz_i
  end subroutine sgemm
 end interface

 interface
  subroutine sgemv ( TRANS, M, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer,intent(in) :: INCX, incy, lda,m,n    !vz_i
   real,intent(in) :: A( LDA, * )    !vz_i
   real,intent(in) :: ALPHA    !vz_i
   real,intent(in) :: BETA    !vz_i
   character*1,intent(in) :: TRANS    !vz_i
   real,intent(in) :: X( * )    !vz_i
   real,intent(inout) :: Y( * )    !vz_i
  end subroutine sgemv
 end interface

 interface
  subroutine sger  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: M
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: X( * )
   real :: Y( * )
  end subroutine sger
 end interface

 interface
  real function snrm2 ( N, X, INCX )
   implicit none
   integer,intent(in) :: INCX,n    !vz_i
   real,intent(in) :: X( * )    !vz_i
  end function snrm2
 end interface

 interface
  subroutine srot (n,sx,incx,sy,incy,c,s)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   real :: c
   real :: s
   real :: sx(*)
   real :: sy(*)
  end subroutine srot
 end interface

 interface
  subroutine srotg(sa,sb,c,s)
   implicit none
   real :: c
   real :: s
   real :: sa
   real :: sb
  end subroutine srotg
 end interface

 interface
  subroutine srotm (N,SX,INCX,SY,INCY,SPARAM)
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   real :: SPARAM(5)
   real :: SX(1)
   real :: SY(1)
  end subroutine srotm
 end interface

 interface
  subroutine srotmg (SD1,SD2,SX1,SY1,SPARAM)
   implicit none
   real :: SD1
   real :: SD2
   real :: SPARAM(5)
   real :: SX1
   real :: SY1
  end subroutine srotmg
 end interface

 interface
  subroutine ssbmv ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: K
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: BETA
   character*1 :: UPLO
   real :: X( * )
   real :: Y( * )
  end subroutine ssbmv
 end interface

 interface
  subroutine sscal(n,sa,sx,incx)
   implicit none
   integer :: incx
   integer :: n
   real :: sa
   real :: sx(*)
  end subroutine sscal
 end interface

 interface
  subroutine sspmv ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   real :: ALPHA
   real :: AP( * )
   real :: BETA
   character*1 :: UPLO
   real :: X( * )
   real :: Y( * )
  end subroutine sspmv
 end interface

 interface
  subroutine sspr  ( UPLO, N, ALPHA, X, INCX, AP )
   implicit none
   integer :: INCX
   integer :: N
   real :: ALPHA
   real :: AP( * )
   character*1 :: UPLO
   real :: X( * )
  end subroutine sspr
 end interface

 interface
  subroutine sspr2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   real :: ALPHA
   real :: AP( * )
   character*1 :: UPLO
   real :: X( * )
   real :: Y( * )
  end subroutine sspr2
 end interface

 interface
  subroutine sswap (n,sx,incx,sy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   real :: sx(*)
   real :: sy(*)
  end subroutine sswap
 end interface

 interface
  subroutine ssymm ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: M
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: B( LDB, * )
   real :: BETA
   real :: C( LDC, * )
   character*1 :: SIDE
   character*1 :: UPLO
  end subroutine ssymm
 end interface

 interface
  subroutine ssymv ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: BETA
   character*1 :: UPLO
   real :: X( * )
   real :: Y( * )
  end subroutine ssymv
 end interface

 interface
  subroutine ssyr  ( UPLO, N, ALPHA, X, INCX, A, LDA )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   character*1 :: UPLO
   real :: X( * )
  end subroutine ssyr
 end interface

 interface
  subroutine ssyr2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   character*1 :: UPLO
   real :: X( * )
   real :: Y( * )
  end subroutine ssyr2
 end interface

 interface
  subroutine ssyr2k( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: B( LDB, * )
   real :: BETA
   real :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine ssyr2k
 end interface

 interface
  subroutine ssyrk ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: BETA
   real :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine ssyrk
 end interface

 interface
  subroutine stbmv ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: K
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   real :: X( * )
  end subroutine stbmv
 end interface

 interface
  subroutine stbsv ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: K
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   real :: X( * )
  end subroutine stbsv
 end interface

 interface
  subroutine stpmv ( UPLO, TRANS, DIAG, N, AP, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   real :: AP( * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   real :: X( * )
  end subroutine stpmv
 end interface

 interface
  subroutine stpsv ( UPLO, TRANS, DIAG, N, AP, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   real :: AP( * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   real :: X( * )
  end subroutine stpsv
 end interface

 interface
  subroutine strmm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine strmm
 end interface

 interface
  subroutine strmv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   real :: X( * )
  end subroutine strmv
 end interface

 interface
  subroutine strsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   real :: A( LDA, * )
   real :: ALPHA
   real :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine strsm
 end interface

 interface
  subroutine strsv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   real :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   real :: X( * )
  end subroutine strsv
 end interface

 !interface
 ! subroutine zaxpy(n,za,zx,incx,zy,incy)
 !  implicit none
 !  integer :: incx
 !  integer :: incy
 !  integer :: n
 !  double complex :: za
 !  double complex :: zx(*)
 !  double complex :: zy(*)
 ! end subroutine zaxpy
 !end interface

 !interface
 ! subroutine  zcopy(n,zx,incx,zy,incy)
 !  implicit none
 !  integer :: incx
 !  integer :: incy
 !  integer :: n
 !  double complex :: zx(*)
 !  double complex :: zy(*)
 ! end subroutine zcopy
 !end interface

 !interface
 ! double complex function zdotc(n,zx,incx,zy,incy)
 !  implicit none
 !  integer :: incx
 !  integer :: incy
 !  integer :: n
 !  double complex :: zx(*)
 !  double complex :: zy(*)
 ! end function zdotc
 !end interface

 interface
  double complex function zdotu(n,zx,incx,zy,incy)
   implicit none
   integer, intent(in) :: incx, incy, n    !vz_i
   double complex, intent(in) :: zx(*), zy(*)    !vz_i
  end function zdotu
 end interface

 interface
  subroutine zdrot( N, CX, INCX, CY, INCY, C, S )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   double precision :: C
   complex*16 :: CX( * )
   complex*16 :: CY( * )
   double precision :: S
  end subroutine zdrot
 end interface

 !interface
 ! subroutine  zdscal(n,da,zx,incx)
 !  implicit none
 !  integer :: incx
 !  integer :: n
 !  double precision :: da
 !  double complex :: zx(*)
 ! end subroutine zdscal
 !end interface

 interface
  subroutine zgbmv ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: KL
   integer :: KU
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: BETA
   character*1 :: TRANS
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zgbmv
 end interface

 !interface
 ! subroutine zgemm ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
 !  implicit none
 !  integer :: K
 !  integer :: LDA
 !  integer :: LDB
 !  integer :: LDC
 !  integer :: M
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: ALPHA
 !  complex*16 :: B( LDB, * )
 !  complex*16 :: BETA
 !  complex*16 :: C( LDC, * )
 !  character*1 :: TRANSA
 !  character*1 :: TRANSB
 ! end subroutine zgemm
 !end interface

 !interface
 ! subroutine zgemv ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
 !  implicit none
 !  integer :: INCX
 !  integer :: INCY
 !  integer :: LDA
 !  integer :: M
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: ALPHA
 !  complex*16 :: BETA
 !  character*1 :: TRANS
 !  complex*16 :: X( * )
 !  complex*16 :: Y( * )
 ! end subroutine zgemv
 !end interface

 !interface
 ! subroutine zgerc ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
 !  implicit none
 !  integer :: INCX
 !  integer :: INCY
 !  integer :: LDA
 !  integer :: M
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: ALPHA
 !  complex*16 :: X( * )
 !  complex*16 :: Y( * )
 ! end subroutine zgerc
 !end interface

 interface
  subroutine zgeru ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zgeru
 end interface

 interface
  subroutine zhbmv ( UPLO, N, K, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: K
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: BETA
   character*1 :: UPLO
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zhbmv
 end interface

 !interface
 ! subroutine zhemm ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
 !  implicit none
 !  integer :: LDA
 !  integer :: LDB
 !  integer :: LDC
 !  integer :: M
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: ALPHA
 !  complex*16 :: B( LDB, * )
 !  complex*16 :: BETA
 !  complex*16 :: C( LDC, * )
 !  character*1 :: SIDE
 !  character*1 :: UPLO
 ! end subroutine zhemm
 !end interface

 interface
  subroutine zhemv ( UPLO, N, ALPHA, A, LDA, X, INCX,&  
 BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: BETA
   character*1 :: UPLO
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zhemv
 end interface

 !interface
 ! subroutine zher  ( UPLO, N, ALPHA, X, INCX, A, LDA )
 !  implicit none
 !  integer :: INCX
 !  integer :: LDA
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  double precision :: ALPHA
 !  character*1 :: UPLO
 !  complex*16 :: X( * )
 ! end subroutine zher
 !end interface

 interface
  subroutine zher2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   character*1 :: UPLO
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zher2
 end interface

 interface
  subroutine zher2k( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,&  
 C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: B( LDB, * )
   double precision :: BETA
   complex*16 :: C( LDC, * )
   character :: TRANS
   character :: UPLO
  end subroutine zher2k
 end interface

 !interface
 ! subroutine zherk( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
 !  implicit none
 !  integer :: K
 !  integer :: LDA
 !  integer :: LDC
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  double precision :: ALPHA
 !  double precision :: BETA
 !  complex*16 :: C( LDC, * )
 !  character :: TRANS
 !  character :: UPLO
 ! end subroutine zherk
 !end interface

 interface
  subroutine zhpmv ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   complex*16 :: ALPHA
   complex*16 :: AP( * )
   complex*16 :: BETA
   character*1 :: UPLO
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zhpmv
 end interface

 interface
  subroutine zhpr  ( UPLO, N, ALPHA, X, INCX, AP )
   implicit none
   integer :: INCX
   integer :: N
   double precision :: ALPHA
   complex*16 :: AP( * )
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine zhpr
 end interface

 interface
  subroutine zhpr2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   complex*16 :: ALPHA
   complex*16 :: AP( * )
   character*1 :: UPLO
   complex*16 :: X( * )
   complex*16 :: Y( * )
  end subroutine zhpr2
 end interface

 interface
  subroutine zrotg(ca,cb,c,s)
   implicit none
   double precision :: c
   double complex :: ca
   double complex :: cb
   double complex :: s
  end subroutine zrotg
 end interface

 !interface
 ! subroutine  zscal(n,za,zx,incx)
 !  implicit none
 !  integer :: incx
 !  integer :: n
 !  double complex :: za
 !  double complex :: zx(*)
 ! end subroutine zscal
 !end interface

 interface
  subroutine  zswap (n,zx,incx,zy,incy)
   implicit none
   integer :: incx
   integer :: incy
   integer :: n
   double complex :: zx(*)
   double complex :: zy(*)
  end subroutine zswap
 end interface

 interface
  subroutine zsymm ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: B( LDB, * )
   complex*16 :: BETA
   complex*16 :: C( LDC, * )
   character*1 :: SIDE
   character*1 :: UPLO
  end subroutine zsymm
 end interface

 interface
  subroutine zsyr2k( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: B( LDB, * )
   complex*16 :: BETA
   complex*16 :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine zsyr2k
 end interface

 interface
  subroutine zsyrk ( UPLO, TRANS, N, K, ALPHA, A, LDA,&  
 BETA, C, LDC )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: BETA
   complex*16 :: C( LDC, * )
   character*1 :: TRANS
   character*1 :: UPLO
  end subroutine zsyrk
 end interface

 interface
  subroutine ztbmv ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: K
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine ztbmv
 end interface

 interface
  subroutine ztbsv ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: K
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine ztbsv
 end interface

 interface
  subroutine ztpmv ( UPLO, TRANS, DIAG, N, AP, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex*16 :: AP( * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine ztpmv
 end interface

 interface
  subroutine ztpsv ( UPLO, TRANS, DIAG, N, AP, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex*16 :: AP( * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine ztpsv
 end interface

 interface
  subroutine ztrmm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&  
 B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: B( LDB, * )
   character*1 :: DIAG
   character*1 :: SIDE
   character*1 :: TRANSA
   character*1 :: UPLO
  end subroutine ztrmm
 end interface

 interface
  subroutine ztrmv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine ztrmv
 end interface

 !interface
 ! subroutine ztrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
 !  implicit none
 !  integer :: LDA
 !  integer :: LDB
 !  integer :: M
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: ALPHA
 !  complex*16 :: B( LDB, * )
 !  character*1 :: DIAG
 !  character*1 :: SIDE
 !  character*1 :: TRANSA
 !  character*1 :: UPLO
 ! end subroutine ztrsm
 !end interface

 interface
  subroutine ztrsv ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
   implicit none
   integer :: INCX
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character*1 :: DIAG
   character*1 :: TRANS
   character*1 :: UPLO
   complex*16 :: X( * )
  end subroutine ztrsv
 end interface

 interface
  subroutine cgetf2( M, N, A, LDA, IPIV, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
  end subroutine cgetf2
 end interface

 interface
  subroutine cgetrf( M, N, A, LDA, IPIV, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
  end subroutine cgetrf
 end interface

 interface
  subroutine cgetri( N, A, LDA, IPIV, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: LWORK
   integer :: N
   complex :: A( LDA, * )
   complex :: WORK( * )
  end subroutine cgetri
 end interface

 interface
  subroutine chpev( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,&  
 INFO )
   implicit none
   integer :: INFO
   integer :: LDZ
   integer :: N
   complex :: AP( * )
   character :: JOBZ
   real :: RWORK( * )
   character :: UPLO
   real :: W( * )
   complex :: WORK( * )
   complex :: Z( LDZ, * )
  end subroutine chpev
 end interface

 interface
  subroutine chptrd( UPLO, N, AP, D, E, TAU, INFO )
   implicit none
   integer :: INFO
   integer :: N
   complex :: AP( * )
   real :: D( * )
   real :: E( * )
   complex :: TAU( * )
   character :: UPLO
  end subroutine chptrd
 end interface

 interface
  complex function cladiv( X, Y )
   implicit none
   complex :: X
   complex :: Y
  end function cladiv
 end interface

 interface
  real function clanhp( NORM, UPLO, N, AP, WORK )
   implicit none
   integer :: N
   complex :: AP( * )
   character :: NORM
   character :: UPLO
   real :: WORK( * )
  end function clanhp
 end interface

 interface
  subroutine clarf( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
   implicit none
   integer :: INCV
   integer :: LDC
   integer :: M
   integer :: N
   complex :: C( LDC, * )
   character :: SIDE
   complex :: TAU
   complex :: V( * )
   complex :: WORK( * )
  end subroutine clarf
 end interface

 interface
  subroutine clarfg( N, ALPHA, X, INCX, TAU )
   implicit none
   integer :: INCX
   integer :: N
   complex :: ALPHA
   complex :: TAU
   complex :: X( * )
  end subroutine clarfg
 end interface

 interface
  subroutine clasr( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
   implicit none
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   real :: C( * )
   character :: DIRECT
   character :: PIVOT
   real :: S( * )
   character :: SIDE
  end subroutine clasr
 end interface

 interface
  subroutine classq( N, X, INCX, SCALE, SUMSQ )
   implicit none
   integer :: INCX
   integer :: N
   real :: SCALE
   real :: SUMSQ
   complex :: X( * )
  end subroutine classq
 end interface

 interface
  subroutine claswp( N, A, LDA, K1, K2, IPIV, INCX )
   implicit none
   integer :: INCX
   integer :: IPIV( * )
   integer :: K1
   integer :: K2
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
  end subroutine claswp
 end interface

 interface
  subroutine clazro( M, N, ALPHA, BETA, A, LDA )
   implicit none
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: ALPHA
   complex :: BETA
  end subroutine clazro
 end interface

 interface
  subroutine csteqr( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDZ
   integer :: N
   character :: COMPZ
   real :: D( * )
   real :: E( * )
   real :: WORK( * )
   complex :: Z( LDZ, * )
  end subroutine csteqr
 end interface

 interface
  subroutine ctrtri( UPLO, DIAG, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   complex :: A( LDA, * )
   character :: DIAG
   character :: UPLO
  end subroutine ctrtri
 end interface

 interface
  subroutine cung2l( M, N, K, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: TAU( * )
   complex :: WORK( * )
  end subroutine cung2l
 end interface

 interface
  subroutine cung2r( M, N, K, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: M
   integer :: N
   complex :: A( LDA, * )
   complex :: TAU( * )
   complex :: WORK( * )
  end subroutine cung2r
 end interface

 interface
  subroutine cupgtr( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDQ
   integer :: N
   complex :: AP( * )
   complex :: Q( LDQ, * )
   complex :: TAU( * )
   character :: UPLO
   complex :: WORK( * )
  end subroutine cupgtr
 end interface

 interface
  subroutine dbdsqr( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,&  
 LDU, C, LDC, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDC
   integer :: LDU
   integer :: LDVT
   integer :: N
   integer :: NCC
   integer :: NCVT
   integer :: NRU
   double precision :: C( LDC, * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: U( LDU, * )
   character :: UPLO
   double precision :: VT( LDVT, * )
   double precision :: WORK( * )
  end subroutine dbdsqr
 end interface

 interface
  subroutine dgebd2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: TAUP( * )
   double precision :: TAUQ( * )
   double precision :: WORK( * )
  end subroutine dgebd2
 end interface

 interface
  subroutine dgebrd( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,&  
 INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: TAUP( * )
   double precision :: TAUQ( * )
   double precision :: WORK( * )
  end subroutine dgebrd
 end interface

 interface
  subroutine dgelq2( M, N, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dgelq2
 end interface

 interface
  subroutine dgelqf( M, N, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dgelqf
 end interface

 interface
  subroutine dgelss( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,&  
 WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LDB
   integer :: LWORK
   integer :: M
   integer :: N
   integer :: NRHS
   integer :: RANK
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   double precision :: RCOND
   double precision :: S( * )
   double precision :: WORK( * )
  end subroutine dgelss
 end interface

 interface
  subroutine dgeqr2( M, N, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dgeqr2
 end interface

 interface
  subroutine dgeqrf( M, N, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dgeqrf
 end interface

 interface
  subroutine dgesvd( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&  
 WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LDU
   integer :: LDVT
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   character :: JOBU
   character :: JOBVT
   double precision :: S( * )
   double precision :: U( LDU, * )
   double precision :: VT( LDVT, * )
   double precision :: WORK( * )
  end subroutine dgesvd
 end interface

 interface
  subroutine dgetf2( M, N, A, LDA, IPIV, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
  end subroutine dgetf2
 end interface

 interface
  subroutine dgetrf( M, N, A, LDA, IPIV, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
  end subroutine dgetrf
 end interface

 interface
  subroutine dgetri( N, A, LDA, IPIV, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: LWORK
   integer :: N
   double precision :: A( LDA, * )
   double precision :: WORK( * )
  end subroutine dgetri
 end interface

 interface
  subroutine dopgtr( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDQ
   integer :: N
   double precision :: AP( * )
   double precision :: Q( LDQ, * )
   double precision :: TAU( * )
   character :: UPLO
   double precision :: WORK( * )
  end subroutine dopgtr
 end interface

 interface
  subroutine dorg2l( M, N, K, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dorg2l
 end interface

 interface
  subroutine dorg2r( M, N, K, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dorg2r
 end interface

 interface
  subroutine dorgbr( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   character :: VECT
   double precision :: WORK( * )
  end subroutine dorgbr
 end interface

 interface
  subroutine dorgl2( M, N, K, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dorgl2
 end interface

 interface
  subroutine dorglq( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dorglq
 end interface

 interface
  subroutine dorgql( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dorgql
 end interface

 interface
  subroutine dorgqr( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   double precision :: WORK( * )
  end subroutine dorgqr
 end interface

 interface
  subroutine dorgtr( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: N
   double precision :: A( LDA, * )
   double precision :: TAU( * )
   character :: UPLO
   double precision :: WORK( LWORK )
  end subroutine dorgtr
 end interface

 interface
  subroutine dorm2r( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
 WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: C( LDC, * )
   character :: SIDE
   double precision :: TAU( * )
   character :: TRANS
   double precision :: WORK( * )
  end subroutine dorm2r
 end interface

 interface
  subroutine dormbr( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,&  
 LDC, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: C( LDC, * )
   character :: SIDE
   double precision :: TAU( * )
   character :: TRANS
   character :: VECT
   double precision :: WORK( * )
  end subroutine dormbr
 end interface

 interface
  subroutine dorml2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
 WORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: C( LDC, * )
   character :: SIDE
   double precision :: TAU( * )
   character :: TRANS
   double precision :: WORK( * )
  end subroutine dorml2
 end interface

 interface
  subroutine dormlq( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
 WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: C( LDC, * )
   character :: SIDE
   double precision :: TAU( * )
   character :: TRANS
   double precision :: WORK( * )
  end subroutine dormlq
 end interface

 interface
  subroutine dormqr( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,&  
 WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: K
   integer :: LDA
   integer :: LDC
   integer :: LWORK
   integer :: M
   integer :: N
   double precision :: A( LDA, * )
   double precision :: C( LDC, * )
   character :: SIDE
   double precision :: TAU( * )
   character :: TRANS
   double precision :: WORK( * )
  end subroutine dormqr
 end interface

 interface
  subroutine dposv( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LDB
   integer :: N
   integer :: NRHS
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: UPLO
  end subroutine dposv
 end interface

 interface
  subroutine dpotf2( UPLO, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character :: UPLO
  end subroutine dpotf2
 end interface

 interface
  subroutine dpotrf( UPLO, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character :: UPLO
  end subroutine dpotrf
 end interface

 interface
  subroutine dpotrs( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LDB
   integer :: N
   integer :: NRHS
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: UPLO
  end subroutine dpotrs
 end interface

 interface
  subroutine dpptrf( UPLO, N, AP, INFO )
   implicit none
   integer :: INFO
   integer :: N
   double precision :: AP( * )
   character :: UPLO
  end subroutine dpptrf
 end interface

 interface
  subroutine dspgst( ITYPE, UPLO, N, AP, BP, INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: N
   double precision :: AP( * )
   double precision :: BP( * )
   character :: UPLO
  end subroutine dspgst
 end interface

 interface
  subroutine dspgv( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,&
 INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: LDZ
   integer :: N
   double precision :: AP( * )
   double precision :: BP( * )
   double precision :: W( * )
   double precision :: WORK( * )
   double precision :: Z( LDZ, * )
   character :: JOBZ
   character :: UPLO
  end subroutine dspgv
 end interface

 interface
  subroutine drscl( N, SA, SX, INCX )
   implicit none
   integer :: INCX
   integer :: N
   double precision :: SA
   double precision :: SX( * )
  end subroutine drscl
 end interface

 interface
  subroutine dspev( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDZ
   integer :: N
   double precision :: AP( * )
   character :: JOBZ
   character :: UPLO
   double precision :: W( * )
   double precision :: WORK( * )
   double precision :: Z( LDZ, * )
  end subroutine dspev
 end interface

 interface
  subroutine dsptrd( UPLO, N, AP, D, E, TAU, INFO )
   implicit none
   integer :: INFO
   integer :: N
   double precision :: AP( * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: TAU( * )
   character :: UPLO
  end subroutine dsptrd
 end interface

 interface
  subroutine dstebz( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E,&  
 M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,&  
 INFO )
   implicit none
   integer :: IBLOCK( * )
   integer :: IL
   integer :: INFO
   integer :: ISPLIT( * )
   integer :: IU
   integer :: IWORK( * )
   integer :: M
   integer :: N
   integer :: NSPLIT
   double precision :: ABSTOL
   double precision :: D( * )
   double precision :: E( * )
   character :: ORDER
   character :: RANGE
   double precision :: VL
   double precision :: VU
   double precision :: W( * )
   double precision :: WORK( * )
  end subroutine dstebz
 end interface

 interface
  subroutine dsteqr( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDZ
   integer :: N
   character :: COMPZ
   double precision :: D( * )
   double precision :: E( * )
   double precision :: WORK( * )
   double precision :: Z( LDZ, * )
  end subroutine dsteqr
 end interface

 interface
  subroutine dsterf( N, D, E, INFO )
   implicit none
   integer :: INFO
   integer :: N
   double precision :: D( * )
   double precision :: E( * )
  end subroutine dsterf
 end interface

 interface
  subroutine dsyev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: N
   double precision :: A( LDA, * )
   character :: JOBZ
   character :: UPLO
   double precision :: W( * )
   double precision :: WORK( * )
  end subroutine dsyev
 end interface

 interface
  subroutine dsygs2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: LDA
   integer :: LDB
   integer :: N
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: UPLO
  end subroutine dsygs2
 end interface

 interface
  subroutine dsygst( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: LDA
   integer :: LDB
   integer :: N
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: UPLO
  end subroutine dsygst
 end interface

 interface
  subroutine dsygv( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,&  
 LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: LDA
   integer :: LDB
   integer :: LWORK
   integer :: N
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: JOBZ
   character :: UPLO
   double precision :: W( * )
   double precision :: WORK( * )
  end subroutine dsygv
 end interface

 interface
  subroutine dsysv( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,&  
 LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: LDB
   integer :: LWORK
   integer :: N
   integer :: NRHS
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: UPLO
   double precision :: WORK( * )
  end subroutine dsysv
 end interface

 interface
  subroutine dsytd2( UPLO, N, A, LDA, D, E, TAU, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: TAU( * )
   character :: UPLO
  end subroutine dsytd2
 end interface

 interface
  subroutine dsytf2( UPLO, N, A, LDA, IPIV, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character :: UPLO
  end subroutine dsytf2
 end interface

 interface
  subroutine dsytrd( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: N
   double precision :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: TAU( * )
   character :: UPLO
   double precision :: WORK( * )
  end subroutine dsytrd
 end interface

 interface
  subroutine dsytrf( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: LWORK
   integer :: N
   double precision :: A( LDA, * )
   character :: UPLO
   double precision :: WORK( * )
  end subroutine dsytrf
 end interface

 interface
  subroutine dsytrs( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: LDB
   integer :: N
   integer :: NRHS
   double precision :: A( LDA, * )
   double precision :: B( LDB, * )
   character :: UPLO
  end subroutine dsytrs
 end interface

 interface
  subroutine dtrti2( UPLO, DIAG, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character :: DIAG
   character :: UPLO
  end subroutine dtrti2
 end interface

 interface
  subroutine dtrtri( UPLO, DIAG, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   double precision :: A( LDA, * )
   character :: DIAG
   character :: UPLO
  end subroutine dtrtri
 end interface

 interface
  double precision function dzsum1( N, CX, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex*16 :: CX( * )
  end function dzsum1
 end interface

 interface
  integer function izmax1( N, CX, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex*16 :: CX( * )
  end function izmax1
 end interface

 interface
  subroutine ssterf( N, D, E, INFO )
   implicit none
   integer :: INFO
   integer :: N
   real :: D( * )
   real :: E( * )
  end subroutine ssterf
 end interface

 interface
  subroutine zbdsqr( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,&  
 LDU, C, LDC, RWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDC
   integer :: LDU
   integer :: LDVT
   integer :: N
   integer :: NCC
   integer :: NCVT
   integer :: NRU
   complex*16 :: C( LDC, * )
   double precision :: D( * )
   double precision :: E( * )
   double precision :: RWORK( * )
   complex*16 :: U( LDU, * )
   character :: UPLO
   complex*16 :: VT( LDVT, * )
  end subroutine zbdsqr
 end interface

 interface
  subroutine zgebak( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,&  
 INFO )
   implicit none
   integer :: IHI
   integer :: ILO
   integer :: INFO
   integer :: LDV
   integer :: M
   integer :: N
   character :: JOB
   double precision :: SCALE( * )
   character :: SIDE
   complex*16 :: V( LDV, * )
  end subroutine zgebak
 end interface

 interface
  subroutine zgebal( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
   implicit none
   integer :: IHI
   integer :: ILO
   integer :: INFO
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character :: JOB
   double precision :: SCALE( * )
  end subroutine zgebal
 end interface

 interface
  subroutine zgebd2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   complex*16 :: TAUP( * )
   complex*16 :: TAUQ( * )
   complex*16 :: WORK( * )
  end subroutine zgebd2
 end interface

 interface
  subroutine zgebrd( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,&  
 INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   complex*16 :: TAUP( * )
   complex*16 :: TAUQ( * )
   complex*16 :: WORK( * )
  end subroutine zgebrd
 end interface

 interface
  subroutine zgees( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,&  
 LDVS, WORK, LWORK, RWORK, BWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LDVS
   integer :: LWORK
   integer :: N
   integer :: SDIM
   complex*16 :: A( LDA, * )
   logical :: BWORK( * )
   character :: JOBVS
   double precision :: RWORK( * )
   logical :: SELECT
   character :: SORT
   complex*16 :: VS( LDVS, * )
   complex*16 :: W( * )
   complex*16 :: WORK( * )
  end subroutine zgees
 end interface

 interface
  subroutine zgeev( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,&
                    WORK, LWORK, RWORK, INFO )
   implicit none
   character :: JOBVL
   character :: JOBVR
   integer :: INFO 
   integer :: LDA
   integer :: LDVL 
   integer :: LDVR  
   integer :: LWORK
   integer :: N
   double precision :: RWORK( * )
   complex*16 :: A( LDA, * )
   complex*16 :: VL( LDVL, * ) 
   complex*16 :: VR( LDVR, * )
   complex*16 :: W( * )
   complex*16 :: WORK( * )
  end subroutine zgeev
 end interface

 interface
  subroutine zgehd2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: IHI
   integer :: ILO
   integer :: INFO
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: TAU( * )
   complex*16 :: WORK( * )
  end subroutine zgehd2
 end interface

 interface
  subroutine zgehrd( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: IHI
   integer :: ILO
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: TAU( * )
   complex*16 :: WORK( LWORK )
  end subroutine zgehrd
 end interface

 interface
  subroutine zgelq2( M, N, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: TAU( * )
   complex*16 :: WORK( * )
  end subroutine zgelq2
 end interface

 interface
  subroutine zgelqf( M, N, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: TAU( * )
   complex*16 :: WORK( * )
  end subroutine zgelqf
 end interface

 interface
  subroutine zgeqr2( M, N, A, LDA, TAU, WORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: TAU( * )
   complex*16 :: WORK( * )
  end subroutine zgeqr2
 end interface

 interface
  subroutine zgeqrf( M, N, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: TAU( * )
   complex*16 :: WORK( * )
  end subroutine zgeqrf
 end interface

 !interface
 ! subroutine zgesv( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: IPIV( * )
 !  integer :: LDA
 !  integer :: LDB
 !  integer :: N
 !  integer :: NRHS
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: B( LDB, * )
 ! end subroutine zgesv
 !end interface

 interface
  subroutine zgesvd( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&  
 WORK, LWORK, RWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LDU
   integer :: LDVT
   integer :: LWORK
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   character :: JOBU
   character :: JOBVT
   double precision :: RWORK( * )
   double precision :: S( * )
   complex*16 :: U( LDU, * )
   complex*16 :: VT( LDVT, * )
   complex*16 :: WORK( * )
  end subroutine zgesvd
 end interface

 interface
  subroutine zgetf2( M, N, A, LDA, IPIV, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
  end subroutine zgetf2
 end interface

 interface
  subroutine zgetri( N, A, LDA, IPIV, WORK, LWORK, INFO )
   implicit none
   integer :: INFO 
   integer :: LDA
   integer :: LWORK
   integer :: N
   integer :: IPIV( * )
   complex*16 :: A( LDA, * )
   complex*16 :: WORK( * )
  end subroutine zgetri
 end interface

 !interface 
 ! subroutine zgetrf( M, N, A, LDA, IPIV, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: IPIV( * )
 !  integer :: LDA
 !  integer :: M
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 ! end subroutine zgetrf
 !end interface

 interface
  subroutine zgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: IPIV( * )
   integer :: LDA
   integer :: LDB
   integer :: N
   integer :: NRHS
   complex*16 :: A( LDA, * )
   complex*16 :: B( LDB, * )
   character :: TRANS
  end subroutine zgetrs
 end interface

 !interface
 ! subroutine zheev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: LDA
 !  integer :: LWORK
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  character :: JOBZ
 !  double precision :: RWORK( * )
 !  character :: UPLO
 !  double precision :: W( * )
 !  complex*16 :: WORK( * )
 ! end subroutine zheev
 !end interface

 interface
  subroutine zhegs2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: LDA
   integer :: LDB
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: B( LDB, * )
   character :: UPLO
  end subroutine zhegs2
 end interface

 !interface
 ! subroutine zhegst( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: ITYPE
 !  integer :: LDA
 !  integer :: LDB
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: B( LDB, * )
 !  character :: UPLO
 ! end subroutine zhegst
 !end interface

 !interface
 ! subroutine zhegv( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,&  
 !  LWORK, RWORK, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: ITYPE
 !  integer :: LDA
 !  integer :: LDB
 !  integer :: LWORK
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  complex*16 :: B( LDB, * )
 !  character :: JOBZ
 !  double precision :: RWORK( * )
 !  character :: UPLO
 !  double precision :: W( * )
 !  complex*16 :: WORK( * )
 ! end subroutine zhegv
 !end interface

 interface
  subroutine zhetd2( UPLO, N, A, LDA, D, E, TAU, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   complex*16 :: TAU( * )
   character :: UPLO
  end subroutine zhetd2
 end interface

 interface
  subroutine zhetrd( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: LWORK
   integer :: N
   complex*16 :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   complex*16 :: TAU( * )
   character :: UPLO
   complex*16 :: WORK( * )
  end subroutine zhetrd
 end interface

 !interface
 ! subroutine zhpev( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: LDZ
 !  integer :: N
 !  complex*16 :: AP( * )
 !  character :: JOBZ
 !  double precision :: RWORK( * )
 !  character :: UPLO
 !  double precision :: W( * )
 !  complex*16 :: WORK( * )
 !  complex*16 :: Z( LDZ, * )
 ! end subroutine zhpev
 !end interface

 !interface
 ! subroutine zhpevx( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,&  
 !  ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK,&  
 !  IFAIL, INFO )
 !  implicit none
 !  integer :: IFAIL( * )
 !  integer :: IL
 !  integer :: INFO
 !  integer :: IU
 !  integer :: IWORK( * )
 !  integer :: LDZ
 !  integer :: M
 !  integer :: N
 !  double precision :: ABSTOL
 !  complex*16 :: AP( * )
 !  character :: JOBZ
 !  character :: RANGE
 !  double precision :: RWORK( * )
 !  character :: UPLO
 !  double precision :: VL
 !  double precision :: VU
 !  double precision :: W( * )
 !  complex*16 :: WORK( * )
 !  complex*16 :: Z( LDZ, * )
 ! end subroutine zhpevx
 !end interface

 interface
  subroutine zhpgst( ITYPE, UPLO, N, AP, BP, INFO )
   implicit none
   integer :: INFO
   integer :: ITYPE
   integer :: N
   complex*16 :: AP( * )
   complex*16 :: BP( * )
   character :: UPLO
  end subroutine zhpgst
 end interface

 !interface
 ! subroutine zhpgv( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,&  
 !  RWORK, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: ITYPE
 !  integer :: LDZ
 !  integer :: N
 !  complex*16 :: AP( * )
 !  complex*16 :: BP( * )
 !  character :: JOBZ
 !  double precision :: RWORK( * )
 !  character :: UPLO
 !  double precision :: W( * )
 !  complex*16 :: WORK( * )
 !  complex*16 :: Z( LDZ, * )
 ! end subroutine zhpgv
 !end interface

 interface
  subroutine zhptrd( UPLO, N, AP, D, E, TAU, INFO )
   implicit none
   integer :: INFO
   integer :: N
   complex*16 :: AP( * )
   double precision :: D( * )
   double precision :: E( * )
   complex*16 :: TAU( * )
   character :: UPLO
  end subroutine zhptrd
 end interface

 interface
  subroutine zhseqr( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,&  
 WORK, LWORK, INFO )
   implicit none
   integer :: IHI
   integer :: ILO
   integer :: INFO
   integer :: LDH
   integer :: LDZ
   integer :: LWORK
   integer :: N
   character :: COMPZ
   complex*16 :: H( LDH, * )
   character :: JOB
   complex*16 :: W( * )
   complex*16 :: WORK( * )
   complex*16 :: Z( LDZ, * )
  end subroutine zhseqr
 end interface

 interface
  subroutine zlabrd( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,&  
 LDY )
   implicit none
   integer :: LDA
   integer :: LDX
   integer :: LDY
   integer :: M
   integer :: N
   integer :: NB
   complex*16 :: A( LDA, * )
   double precision :: D( * )
   double precision :: E( * )
   complex*16 :: TAUP( * )
   complex*16 :: TAUQ( * )
   complex*16 :: X( LDX, * )
   complex*16 :: Y( LDY, * )
  end subroutine zlabrd
 end interface

 interface
  subroutine zlacgv( N, X, INCX )
   implicit none
   integer :: INCX
   integer :: N
   complex*16 :: X( * )
  end subroutine zlacgv
 end interface

 interface
  subroutine zlacon( N, V, X, EST, KASE )
   implicit none
   integer :: KASE
   integer :: N
   double precision :: EST
   complex*16 :: V( N )
   complex*16 :: X( N )
  end subroutine zlacon
 end interface

 interface
  subroutine zlacpy( UPLO, M, N, A, LDA, B, LDB )
   implicit none
   integer :: LDA
   integer :: LDB
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: B( LDB, * )
   character :: UPLO
  end subroutine zlacpy
 end interface

 interface
  double complex function zladiv( X, Y )
   implicit none
   complex*16 :: X
   complex*16 :: Y
  end function zladiv
 end interface

 interface
  subroutine zlahqr( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,&  
 IHIZ, Z, LDZ, INFO )
   implicit none
   integer :: IHI
   integer :: IHIZ
   integer :: ILO
   integer :: ILOZ
   integer :: INFO
   integer :: LDH
   integer :: LDZ
   integer :: N
   complex*16 :: H( LDH, * )
   complex*16 :: W( * )
   logical :: WANTT
   logical :: WANTZ
   complex*16 :: Z( LDZ, * )
  end subroutine zlahqr
 end interface

 interface
  subroutine zlahrd( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDT
   integer :: LDY
   integer :: N
   integer :: NB
   complex*16 :: A( LDA, * )
   complex*16 :: T( LDT, NB )
   complex*16 :: TAU( NB )
   complex*16 :: Y( LDY, NB )
  end subroutine zlahrd
 end interface

 interface
  subroutine zlahr2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
   implicit none
   integer :: K
   integer :: LDA
   integer :: LDT
   integer :: LDY
   integer :: N
   integer :: NB
   complex*16 :: A( LDA, * )
   complex*16 :: T( LDT, NB )
   complex*16 :: TAU( NB )
   complex*16 :: Y( LDY, NB )
  end subroutine zlahr2
 end interface

 interface
  double precision function zlange( NORM, M, N, A, LDA, WORK )
   implicit none
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   character :: NORM
   double precision :: WORK( * )
  end function zlange
 end interface

 interface
  double precision function zlanhe( NORM, UPLO, N, A, LDA, WORK )
   implicit none
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character :: NORM
   character :: UPLO
   double precision :: WORK( * )
  end function zlanhe
 end interface

 interface
  double precision function zlanhp( NORM, UPLO, N, AP, WORK )
   implicit none
   integer :: N
   complex*16 :: AP( * )
   character :: NORM
   character :: UPLO
   double precision :: WORK( * )
  end function zlanhp
 end interface

 interface
  double precision function zlanhs( NORM, N, A, LDA, WORK )
   implicit none
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character :: NORM
   double precision :: WORK( * )
  end function zlanhs
 end interface

 interface
  subroutine zlaqr0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,&
 IHIZ, Z, LDZ, WORK, LWORK, INFO )
   implicit none
   integer :: IHI
   integer :: IHIZ
   integer :: ILO
   integer :: ILOZ
   integer :: INFO
   integer :: LDH
   integer :: LDZ
   integer :: LWORK
   integer :: N
   
   complex*16 :: H( LDH, * )
   complex*16 :: W( * )
   complex*16 :: WORK( * )
   complex*16 :: Z( LDZ, * )

   logical :: WANTT
   logical :: WANTZ
  end subroutine zlaqr0
 end interface

 interface
  subroutine zlaqr1( N, H, LDH, S1, S2, V )
   implicit none
   integer :: LDH
   integer :: N
   complex*16 :: S1 
   complex*16 :: S2
   complex*16 :: H( LDH, * )
   complex*16 :: V( * )
  end subroutine zlaqr1
 end interface

 interface
  subroutine zlaqr2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,&
 IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,&
 NV, WV, LDWV, WORK, LWORK )
   implicit none
   integer :: IHIZ, ILOZ, KBOT, KTOP, LDH, LDT
   integer :: LDV, LDWV, LDZ, LWORK, N, ND, NH
   integer :: NS, NV, NW
   logical :: WANTT, WANTZ
   complex*16 :: H( LDH, * ), SH( * ), T( LDT, * )
   complex*16 :: V( LDV, * ),WORK( * ), WV( LDWV, * ), Z( LDZ, * )
  end subroutine zlaqr2
 end interface

 interface
  subroutine zlaqr3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,&
 IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,&
 NV, WV, LDWV, WORK, LWORK )
   implicit none
   integer :: IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV
   integer :: LDZ, LWORK, N, ND, NH, NS, NV, NW
   logical :: WANTT, WANTZ
   complex*16 :: H( LDH, * ), SH( * ), T( LDT, * )
   complex*16 :: V( LDV, * ),WORK( * ), WV( LDWV, * ), Z( LDZ, * )
  end subroutine zlaqr3
 end interface

 interface
  subroutine zlaqr4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,&
 IHIZ, Z, LDZ, WORK, LWORK, INFO )
   implicit none
   integer :: IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
   logical :: WANTT, WANTZ
   complex*16 :: H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
  end subroutine zlaqr4
 end interface

 interface
  subroutine zlaqr5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,&
 H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,&
 WV, LDWV, NH, WH, LDWH )
   implicit none
   integer :: IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV
   integer :: LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
   logical :: WANTT, WANTZ
   complex*16 :: H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * )
   complex*16 :: WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
  end subroutine zlaqr5
 end interface

 interface
  subroutine zlarf( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
   implicit none
   integer :: INCV
   integer :: LDC
   integer :: M
   integer :: N
   complex*16 :: C( LDC, * )
   character :: SIDE
   complex*16 :: TAU
   complex*16 :: V( * )
   complex*16 :: WORK( * )
  end subroutine zlarf
 end interface

 interface
  subroutine zlarfb( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,&  
 T, LDT, C, LDC, WORK, LDWORK )
   implicit none
   integer :: K
   integer :: LDC
   integer :: LDT
   integer :: LDV
   integer :: LDWORK
   integer :: M
   integer :: N
   complex*16 :: C( LDC, * )
   character :: DIRECT
   character :: SIDE
   character :: STOREV
   complex*16 :: T( LDT, * )
   character :: TRANS
   complex*16 :: V( LDV, * )
   complex*16 :: WORK( LDWORK, * )
  end subroutine zlarfb
 end interface

 interface
  subroutine zlarfg( N, ALPHA, X, INCX, TAU )
   implicit none
   integer :: INCX
   integer :: N
   complex*16 :: ALPHA
   complex*16 :: TAU
   complex*16 :: X( * )
  end subroutine zlarfg
 end interface

 interface
  subroutine zlarft( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
   implicit none
   integer :: K
   integer :: LDT
   integer :: LDV
   integer :: N
   character :: DIRECT
   character :: STOREV
   complex*16 :: T( LDT, * )
   complex*16 :: TAU( * )
   complex*16 :: V( LDV, * )
  end subroutine zlarft
 end interface

 interface
  subroutine zlarfx( SIDE, M, N, V, TAU, C, LDC, WORK )
   implicit none
   integer :: LDC
   integer :: M
   integer :: N
   complex*16 :: C( LDC, * )
   character :: SIDE
   complex*16 :: TAU
   complex*16 :: V( * )
   complex*16 :: WORK( * )
  end subroutine zlarfx
 end interface

 interface
  subroutine zlartg( F, G, CS, SN, R )
   implicit none
   double precision :: CS
   complex*16 :: F
   complex*16 :: G
   complex*16 :: R
   complex*16 :: SN
  end subroutine zlartg
 end interface

 interface
  subroutine zlascl( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: KL
   integer :: KU
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   double precision :: CFROM
   double precision :: CTO
   character :: TYPE
  end subroutine zlascl
 end interface

 interface
  subroutine zlaset( UPLO, M, N, ALPHA, BETA, A, LDA )
   implicit none
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: BETA
   character :: UPLO
  end subroutine zlaset
 end interface

 interface
  subroutine zlasr( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
   implicit none
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   double precision :: C( * )
   character :: DIRECT
   character :: PIVOT
   double precision :: S( * )
   character :: SIDE
  end subroutine zlasr
 end interface

 interface
  subroutine zlassq( N, X, INCX, SCALE, SUMSQ )
   implicit none
   integer :: INCX
   integer :: N
   double precision :: SCALE
   double precision :: SUMSQ
   complex*16 :: X( * )
  end subroutine zlassq
 end interface

 interface
  subroutine zlaswp( N, A, LDA, K1, K2, IPIV, INCX )
   implicit none
   integer :: INCX
   integer :: IPIV( * )
   integer :: K1
   integer :: K2
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
  end subroutine zlaswp
 end interface

 interface
  subroutine zlatrd( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
   implicit none
   integer :: LDA
   integer :: LDW
   integer :: N
   integer :: NB
   complex*16 :: A( LDA, * )
   double precision :: E( * )
   complex*16 :: TAU( * )
   character :: UPLO
   complex*16 :: W( LDW, * )
  end subroutine zlatrd
 end interface

 interface
  subroutine zlatrs( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,&
 CNORM, INFO )
   implicit none
   character :: DIAG, NORMIN, TRANS, UPLO
   integer :: INFO, LDA, N
   double precision :: SCALE      
   double precision :: CNORM( * )
   complex*16 :: A( LDA, * ), X( * )
  end subroutine zlatrs
 end interface

 interface
  subroutine zlazro( M, N, ALPHA, BETA, A, LDA )
   implicit none
   integer :: LDA
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: ALPHA
   complex*16 :: BETA
  end subroutine zlazro
 end interface

 interface
  subroutine zpotf2( UPLO, N, A, LDA, INFO )
   implicit none
   integer :: INFO
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
   character :: UPLO
  end subroutine zpotf2
 end interface

 !interface
 ! subroutine zpotrf( UPLO, N, A, LDA, INFO )
 !  implicit none
 !  integer :: INFO
 !  integer :: LDA
 !  integer :: N
 !  complex*16 :: A( LDA, * )
 !  character :: UPLO
 ! end subroutine zpotrf
 !end interface

 interface
  subroutine zpptrf( UPLO, N, AP, INFO )
   implicit none
   integer :: INFO
   integer :: N
   complex*16 :: AP( * )
   character :: UPLO
  end subroutine zpptrf
 end interface

 interface
  subroutine zrot( N, CX, INCX, CY, INCY, C, S )
   implicit none
   integer :: INCX
   integer :: INCY
   integer :: N
   double precision :: C
   complex*16 :: CX( * )
   complex*16 :: CY( * )
   complex*16 :: S
  end subroutine zrot
 end interface

 interface
  subroutine ztrevc( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,&
                     LDVR, MM, M, WORK, RWORK, INFO )
   implicit none
   character :: HOWMNY, SIDE
   integer :: INFO, LDT, LDVL, LDVR, M, MM, N
   logical :: SELECT( * )
   double precision :: RWORK( * )
   complex*16 :: T( LDT, * ), VL( LDVL, * ), VR( LDVR, * )
   complex*16 :: WORK( * )
  end subroutine ztrevc
 end interface

 interface
  subroutine ztrexc( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
   implicit none
   integer :: IFST
   integer :: ILST
   integer :: INFO
   integer :: LDQ
   integer :: LDT
   integer :: N
   character :: COMPQ
   complex*16 :: Q( LDQ, * )
   complex*16 :: T( LDT, * )
  end subroutine ztrexc
 end interface

 interface
  subroutine ztrsen( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,&  
 SEP, WORK, LWORK, INFO )
   implicit none
   integer :: INFO
   integer :: LDQ
   integer :: LDT
   integer :: LWORK
   integer :: M
   integer :: N
   character :: COMPQ
   character :: JOB
   complex*16 :: Q( LDQ, * )
   double precision :: S
   logical :: SELECT( * )
   double precision :: SEP
   complex*16 :: T( LDT, * )
   complex*16 :: W( * )
   complex*16 :: WORK( * )
  end subroutine ztrsen
 end interface

 interface
  subroutine ztrsyl( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,&  
 LDC, SCALE, INFO )
   implicit none
   integer :: INFO
   integer :: ISGN
   integer :: LDA
   integer :: LDB
   integer :: LDC
   integer :: M
   integer :: N
   complex*16 :: A( LDA, * )
   complex*16 :: B( LDB, * )
   complex*16 :: C( LDC, * )
   double precision :: SCALE
   character :: TRANA
   character :: TRANB
  end subroutine ztrsyl
 end interface

 interface
  subroutine ztrti2( UPLO, DIAG, N, A, LDA, INFO )
   implicit none
   character :: DIAG
   character :: UPLO
   integer :: INFO 
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
  end subroutine ztrti2
 end interface

 interface
  subroutine ztrtri( UPLO, DIAG, N, A, LDA, INFO )
   implicit none
   character :: DIAG
   character :: UPLO
   integer :: INFO
   integer :: LDA
   integer :: N
   complex*16 :: A( LDA, * )
  end subroutine ztrtri
 end interface

 interface
  subroutine zungqr( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
   implicit none
   integer :: INFO, K, LDA, LWORK, M, N
   complex*16 :: A( LDA, * ), TAU( * ), WORK( * )
  end subroutine zungqr
 end interface

end module m_linalg_interfaces
!!***

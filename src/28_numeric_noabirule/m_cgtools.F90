!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_cgtools
!! NAME
!!  m_cgtools
!!
!! FUNCTION
!! This module defines wrappers for BLAS routines. The arguments are stored
!! using the "cg" convention, namely real array of shape cg(2,...)
!!
!! COPYRIGHT
!! Copyright (C) 1992-2018 ABINIT group (MG, MT, XG, DCA, GZ, FB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! 1) The convention about names of interfaced routine is: cg_<name>,
!!    where <name> is equal to the name of the standard BLAS routine
!!
!! 2) Blas routines are called without an explicit interface on purpose since
!!
!!    a) The compiler should pass the base address of the array to the F77 BLAS
!!
!!    b) Any compiler would complain about type mismatch (REAL,COMPLEX)
!!       if an explicit interface is given.
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

#if defined HAVE_LINALG_GEMM3M
#define ABI_ZGEMM ZGEMM3M
#else
#define ABI_ZGEMM ZGEMM
#endif

MODULE m_cgtools

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_fstrings,   only : toupper

 implicit none

 private

 real(dp),public,parameter :: cg_czero(2) = (/0._dp,0._dp/)
 real(dp),public,parameter :: cg_cone(2)  = (/1._dp,0._dp/)


 ! Helper functions.
 public :: cg_setval
 public :: cg_tocplx
 public :: cg_fromcplx
 public :: cg_filter
 public :: cg_setaug_zero
 public :: cg_to_reim
 public :: cg_from_reim
 !public :: cg_times_eigr

 ! Blas1
 public :: cg_zcopy
 public :: cg_zscal
 public :: cg_dznrm2
 public :: cg_zdotc
 public :: cg_real_zdotc
 public :: cg_zdotu
 public :: cg_zaxpy
 public :: cg_zaxpby

 ! Blas2
 public :: cg_zgemv

 ! Blas3
 public :: cg_zgemm

 ! Helper functions for DFT calculations.
 public :: set_istwfk               ! Returns the value of istwfk associated to the input k-point.
 public :: sqnorm_g                 ! Square of the norm in reciprocal space.
 public :: dotprod_g                ! Scalar product <vec1|vect2> of complex vectors vect1 and vect2 (can be the same)
 public :: matrixelmt_g             ! matrix element <wf1|O|wf2> of two wavefunctions, in reciprocal space, for an operator diagonal in G-space.
 public :: dotprod_v                ! Dot product of two potentials (integral over FFT grid).
 public :: sqnorm_v                 ! Compute square of the norm of a potential (integral over FFT grid).
 public :: mean_fftr                ! Compute the mean of an arraysp(nfft,nspden), over the FFT grid.
 public :: cg_getspin               ! Sandwich a single wave function on the Pauli matrices
 public :: cg_gsph2box              ! Transfer data from the G-sphere to the FFT box.
 public :: cg_box2gsph              ! Transfer data from the FFT box to the G-sphere
 public :: cg_addtorho              ! Add |ur|**2 to the ground-states density rho.
 public :: cg_vlocpsi               ! Apply the local part of the potential to the wavefunction in real space.
 public :: cgnc_cholesky            ! Cholesky orthonormalization (version optimized for NC wavefunctions).
 public :: cgpaw_cholesky           ! Cholesky orthonormalization of PAW wavefunctions.
 public :: cgnc_normalize           ! Normalize NC wavefunctions.
 public :: cgnc_gramschmidt         ! Gram-Schmidt orthogonalization for NC wavefunctions.
 public :: cgpaw_normalize          ! Normalize PAW wavefunctions.
 public :: cgpaw_gramschmidt        ! Gram-Schmidt orthogonalization for PAW wavefuncion
 public :: projbd                   ! Project out vector "direc" onto the bands i.e. direc=direc-$sum_{j/=i} { <cg_{j}|direc>.|cg_{j}> }$
 public :: cg_envlop                ! Multiply random number values in cg by envelope function to lower initial kinetic energy.
 public :: cg_normev                ! Normalize a set of num eigenvectors of complex length ndim
 public :: cg_precon                ! precondition $<G|(H-e_{n,k})|C_{n,k}>$
 public :: cg_precon_block          ! precondition $<G|(H-e_{n,k})|C_{n,k}>$ for a block of band in the case of real WFs (istwfk/=1)
 public :: cg_zprecon_block         ! precondition $<G|(H-e_{n,k})|C_{n,k}>$ for a block of band
 public :: fxphas_seq               ! Fix phase of all bands. Keep normalization but maximize real part
 public :: overlap_g                ! Compute the scalar product between WF at two different k-points
!***

 !integer,parameter,private :: MIN_SIZE = 5000
 !complex(spc),private,parameter :: czero_spc =(0._sp,0._sp)
 !complex(spc),private,parameter :: cone_spc  =(1._sp,0._sp)
 !complex(dpc),private,parameter :: czero_dpc =(0._dp,0._dp)
 !complex(dpc),private,parameter :: cone_dpc  =(1._dp,0._dp)

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_setval
!! NAME
!!  cg_setval
!!
!! FUNCTION
!!  Set cg=alpha.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_setval(n,cg,alpha)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_setval'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),optional,intent(in) :: alpha(2)
 real(dp),intent(inout) :: cg(2,n)

! *************************************************************************

 if (PRESENT(alpha)) then
!$OMP PARALLEL
!$OMP WORKSHARE
   cg(1,:)=alpha(1)
   cg(2,:)=alpha(2)
!$OMP END WORKSHARE
!$OMP END PARALLEL
 else
!$OMP PARALLEL
!$OMP WORKSHARE
   cg(:,:)=zero
!$OMP END WORKSHARE
!$OMP END PARALLEL
 end if

end subroutine cg_setval
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_tocplx
!! NAME
!!  cg_tocplx
!!
!! FUNCTION
!!  Convert a real array with (real,imag) part to complex.
!!
!! INPUTS
!!  n = Specifies the number of elements in cg and ocplx
!!  cg(2*n)=Input array with real and imaginary part.
!!
!! OUTPUT
!!  ocplx(n)=Output complex array.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_tocplx(n, cg, ocplx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_tocplx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(in) :: cg(2*n)
 complex(dpc),intent(out) :: ocplx(n)

!Local variables ------------------------------
!scalars
 integer :: ii,idx

! *************************************************************************

!$OMP PARALLEL DO PRIVATE(ii,idx)
 do ii=1,n
   idx = 2*ii-1
   ocplx(ii) = DCMPLX(cg(idx),cg(idx+1))
 end do

end subroutine cg_tocplx
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_fromcplx
!! NAME
!!  cg_fromcplx
!!
!! FUNCTION
!!  Convert a complex array to a real array with (real,imag) part
!!
!! INPUTS
!!  n = Specifies the number of elements in icplx and ocg.
!!  icplx(n)=Input complex array.
!!
!! OUTPUT
!!  ocg(2*n)=Output array with real and imaginary part.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_fromcplx(n,icplx,ocg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_fromcplx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(out) :: ocg(2*n)
 complex(dpc),intent(in) :: icplx(n)

!Local variables ------------------------------
!scalars
 integer :: ii,idx

! *************************************************************************

!$OMP PARALLEL DO PRIVATE(ii,idx)
 do ii=1,n
   idx = 2*ii-1
   ocg(idx  ) = DBLE (icplx(ii))
   ocg(idx+1) = AIMAG(icplx(ii))
 end do

end subroutine cg_fromcplx
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_filter
!! NAME
!!  cg_filter
!!
!! FUNCTION
!!  Set all the elements of x to zero where mask is .TRUE.
!!
!! INPUTS
!!  n=Specifies the number of elements in vectors x and y.
!!  mask(n)=Logical array.
!!
!! SIDE EFFECTS
!!  x(n)=See description.
!!
!! PARENTS
!!
!! SOURCE

pure subroutine cg_filter(n, x, mask)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_filter'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(inout) :: x(2,n)
 logical,intent(in) :: mask(n)

! *************************************************************************

 where (mask)
   x(1,:) = zero
   x(2,:) = zero
 end where

end subroutine cg_filter
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_setaug_zero
!! NAME
!!  cg_setaug_zero
!!
!! FUNCTION
!!  Set to zero all elements of the array that are not in the FFT box.
!!
!! INPUTS
!! nx,ny,nz=physical dimensions of the FFT box
!! ldx,ldy,ldx=memory dimension of arr
!! ndat=number of FFTs
!!
!! SIDE EFFECT
!!  arr(2,ldx,ldy,ldz*ndat)= all entries in the augmented region are set to zero
!!
!! PARENTS
!!
!! SOURCE

pure subroutine cg_setaug_zero(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,arr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_setaug_zero'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 real(dp),intent(inout) :: arr(cplex,ldx,ldy,ldz*ndat)

!Local variables-------------------------------
 integer :: iy,iz,dat,padat

! *************************************************************************

 if (nx /= ldx) then
   do iz=1,ldz*ndat
     do iy=1,ldy
       arr(:,nx+1:ldx,iy,iz) = zero
     end do
   end do
 end if

 if (ny /= ldy) then
   do iz=1,ldz*ndat
     arr(:,:,ny+1:ldy,iz) = zero
   end do
 end if

 if (nz /= ldz) then
   do dat=1,ndat
     padat = ldz*(dat-1)
     do iz=nz+1,ldz
       arr(:,:,:,iz+padat) = zero
     end do
   end do
 end if

end subroutine cg_setaug_zero
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_to_reim
!! NAME
!!  cg_to_reim
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_to_reim(npw,ndat,cg,factor,reim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_to_reim'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,ndat
 real(dp),intent(in) :: factor
!arrays
 real(dp),intent(in) :: cg(2*npw*ndat)
 real(dp),intent(out) :: reim(npw*ndat*2)

! *************************************************************************

 ! Pack real and imaginary part of the wavefunctions.
 ! and multiply by scale factor if factor /= one. Could block but oh well
 call dcopy(npw*ndat, cg(1), 2, reim(1), 1)
 if (factor /= one) call dscal(npw*ndat,factor,reim(1),1)

 call dcopy(npw*ndat, cg(2), 2, reim(npw*ndat+1), 1)
 if (factor /= one) call dscal(npw*ndat,factor,reim(npw*ndat+1),1)

end subroutine cg_to_reim
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_from_reim
!! NAME
!!  cg_from_reim
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_from_reim(npw,ndat,reim,factor,cg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_from_reim'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,ndat
 real(dp),intent(in) :: factor
!arrays
 real(dp),intent(in) :: reim(npw*ndat*2)
 real(dp),intent(out) :: cg(2*npw*ndat)

! *************************************************************************

 ! UnPack real and imaginary part and multiply by scale factor if /= one.
 ! Could use blocking but oh well
 call dcopy(npw*ndat, reim(1), 1, cg(1), 2)
 call dcopy(npw*ndat, reim(npw*ndat+1), 1, cg(2), 2)

 if (factor /= one) call dscal(2*npw*ndat,factor, cg(1), 1)

end subroutine cg_from_reim
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zcopy
!! NAME
!!  cg_zcopy
!!
!! FUNCTION
!!  Perform y = x, where x and y are vectors.
!!
!! INPUTS
!!  n = Specifies the number of elements in vectors x and y.
!!  x = Input Array
!!
!! OUTPUT
!!  y = In output, y contains a copy of the values of x.
!!
!! PARENTS
!!      cgwf,corrmetalwf1,dfpt_cgwf,dfpt_mkrho,dfpt_vtowfk,lapackprof,m_iowf
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zcopy(n, x, y)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zcopy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(in) :: x(2*n)
 real(dp),intent(out) :: y(2*n)

! *************************************************************************

 call zcopy(n,x,1,y,1)

end subroutine cg_zcopy
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zscal
!! NAME
!!  cg_zscal
!!
!! FUNCTION
!!  Perform x = a*x
!!
!! INPUTS
!!  n = Specifies the number of elements in vector x.
!!  a(2)= The scalar a. If a(2) is zero, x = a*x is computed via zdscal
!!
!! OUTPUT
!!  x = Updated vector.
!!
!! OUTPUT
!!
!! PARENTS
!!      cgwf,m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zscal(n, a, x)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zscal'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: a(2)
!arrays
 real(dp),intent(inout) :: x(2*n)

! *************************************************************************

 if (a(2) == zero) then
   call zdscal(n, a, x, 1)
 else
   call zscal(n, a, x, 1)
 end if

end subroutine cg_zscal
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_dznrm2
!! NAME
!!  cg_dznrm2
!!
!! FUNCTION
!!   returns the euclidean norm of a vector via the function name, so that
!!   DZNRM2 := sqrt( x**H*x )
!!
!! INPUTS
!!  n = Specifies the number of elements in vector x.
!!  x(2*x) = Input array.
!!
!! OUTPUT
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

function cg_dznrm2(n, x) result(res)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_dznrm2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp) :: res
!arrays
 real(dp),intent(in) :: x(2*n)
 real(dp),external :: dznrm2

! *************************************************************************

 res = dznrm2(n, x, 1)

end function cg_dznrm2
!!***
!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zdotc
!! NAME
!!  cg_zdotc
!!
!! FUNCTION
!!   Perform a vector-vector operation defined as res = \Sigma (conjg(x)*y) where x and y are n-element vectors.
!!
!! INPUTS
!!  n = Specifies the number of elements in vector x and y
!!  x,y = Input arrays.
!!
!! OUTPUT
!!  res(2)=Real and Imaginary part of the scalar product.
!!
!! PARENTS
!!
!! SOURCE

function cg_zdotc(n,x,y) result(res)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zdotc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(in) :: x(2,n)
 real(dp),intent(in) :: y(2,n)
 real(dp) :: res(2)

!Local variables-------------------------------
#ifdef HAVE_LINALG_ZDOTC_BUG
 integer :: ii
#endif
 complex(dpc) :: cres
 complex(dpc),external :: zdotc

! *************************************************************************

#ifdef HAVE_LINALG_ZDOTC_BUG
 ! Workaround for veclib on MacOSx
 res = zero
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:res)
 do ii=1,n
   res(1) = res(1) + x(1,ii)*y(1,ii) + x(2,ii)*y(2,ii)
   res(2) = res(2) + x(1,ii)*y(2,ii) - x(2,ii)*y(1,ii)
 end do

#else
 cres = zdotc(n, x, 1, y, 1)
 res(1) = REAL(cres)
 res(2) = AIMAG(cres)
#endif

end function cg_zdotc
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_real_zdotc
!! NAME
!!  cg_real_zdotc
!!
!! FUNCTION
!!   Perform a vector-vector operation defined as res = REAL (\Sigma (conjg(x)*y)) where x and y are n-element vectors.
!!
!! INPUTS
!!  n = Specifies the number of elements in vector x and y
!!  x,y = Input arrays.
!!
!! OUTPUT
!!  res=Real part of the scalar product.
!!
!! PARENTS
!!
!! SOURCE

function cg_real_zdotc(n,x,y) result(res)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_real_zdotc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(in) :: x(2,n)
 real(dp),intent(in) :: y(2,n)
 real(dp) :: res

!Local variables-------------------------------
 real(dp),external :: ddot

! *************************************************************************

 res = ddot(2*n,x,1,y,1)

end function cg_real_zdotc
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zdotu
!! NAME
!!  cg_zdotu
!!
!! FUNCTION
!!   Perform a vector-vector operation defined as res = \Sigma (x*y) where x and y are n-element vectors.
!!   Note that x is unconjugated.
!!
!! INPUTS
!!  n = Specifies the number of elements in vector x and y
!!  x,y = Input arrays.
!!
!! OUTPUT
!!  res(2)=Real and Imaginary part of the scalar product.
!!
!! PARENTS
!!
!! SOURCE

function cg_zdotu(n, x, y) result(res)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zdotu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(in) :: x(2,n)
 real(dp),intent(in) :: y(2,n)
 real(dp) :: res(2)

!Local variables-------------------------------
#ifdef HAVE_LINALG_ZDOTU_BUG
 integer :: ii
#endif
 complex(dpc) :: cres
 complex(dpc),external :: zdotu

! *************************************************************************

#ifdef HAVE_LINALG_ZDOTU_BUG
 ! Workaround for veclib on MacOSx
 res = zero
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:res)
 do ii=1,n
   res(1) = res(1) + x(1,ii)*y(1,ii) - x(2,ii)*y(2,ii)
   res(2) = res(2) + x(1,ii)*y(2,ii) + x(2,ii)*y(1,ii)
 end do
#else
 cres = zdotu(n, x, 1, y, 1)
 res(1) = REAL(cres)
 res(2) = AIMAG(cres)
#endif

end function cg_zdotu
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zaxpy
!! NAME
!!  cg_zaxpy
!!
!! FUNCTION
!!  Computes y = alpha*x + y
!!
!! INPUTS
!!  n = Specifies the number of elements in vectors x and y.
!!  alpha = Specifies the scalar alpha.
!!  x = Array
!!
!! SIDE EFFECTS
!!  y = Array. In output, y contains the updated vector.
!!
!! PARENTS
!!      cgwf,dfpt_cgwf,lapackprof,rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zaxpy(n,alpha,x,y)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zaxpy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: alpha(2)
!arrays
 real(dp),intent(in) :: x(2*n)
 real(dp),intent(inout) :: y(2*n)

!local variables
! integer :: ii

! *************************************************************************

 if (alpha(2) == zero) then
   call daxpy(2*n,alpha(1),x,1,y,1)
 else
   call zaxpy(n,alpha,x,1,y,1)
 end if

end subroutine cg_zaxpy
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zaxpby
!! NAME
!!  cg_zaxpby
!!
!! FUNCTION
!!  Scales two vectors, adds them to one another and stores result in the vector.
!!  y := a*x + b*y
!!
!! INPUTS
!! n = the number of elements in vectors x and y.
!! a = Specifies the scalar a.
!! x = Array.
!! b = Specifies the scalar b.
!! y = Array
!!
!! OUTPUT
!! y Contains the updated vector y.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zaxpby(n,a,x,b,y)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zaxpby'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: a(2),b(2)
!arrays
 real(dp),intent(in) :: x(2*n)
 real(dp),intent(inout) :: y(2*n)

! *************************************************************************

#ifdef HAVE_LINALG_AXPBY
 call zaxpby(n, a, x, 1, b, y, 1)
#else
 call zscal(n, b, y, 1)
 call zaxpy(n, a, x, 1, y,1)
#endif

end subroutine cg_zaxpby
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zgemv
!! NAME
!!  cg_zgemv
!!
!! FUNCTION
!! The ?gemv routines perform a matrix-vector operation defined as
!!
!! y := alpha*A*x + beta*y,
!! or
!! y := alpha*A'*x + beta*y,
!! or
!! y := alpha*conjg(A')*x + beta*y,
!!
!! where: alpha and beta are scalars, x and y are vectors, A is an m-by-n matrix.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      lapackprof,m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zgemv(trans,nrows,ncols,cgmat,vec,matvec,alpha,beta)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zgemv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nrows,ncols
 real(dp),optional,intent(in) :: alpha(2),beta(2)
 character(len=1),intent(in) :: trans
!arrays
 real(dp),intent(in) :: cgmat(2,nrows*ncols)
 real(dp),intent(in) :: vec(2,*)
 real(dp),intent(inout) :: matvec(2,*)

!Local variables-------------------------------
!scalars
 integer :: mm,nn,kk,lda,ldb,ldc
 real(dp) :: my_alpha(2),my_beta(2)

! *************************************************************************

 lda = nrows
 mm  = nrows
 nn  = 1
 kk  = ncols

 if (toupper(trans) /= 'N') then
   mm = ncols
   kk = nrows
 end if

 ldb = kk
 ldc = mm

 my_alpha = cg_cone;  if (PRESENT(alpha)) my_alpha = alpha
 my_beta  = cg_czero; if (PRESENT(beta))  my_beta  = beta

 call ZGEMM(trans,"N",mm,nn,kk,my_alpha,cgmat,lda,vec,ldb,my_beta,matvec,ldc)
 ! ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

 !call ZGEMV(trans,mm,nn,my_alpha,cgmat,lda,vec,1,my_beta,matvec,1)

end subroutine cg_zgemv
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zgemm
!! NAME
!!  cg_zgemm
!!
!! FUNCTION
!!  The ?gemm routines perform a matrix-matrix operation with general matrices.
!!  The operation is defined as C := alpha*op(A)*op(B) + beta*C,
!!  where:
!!
!!  op(x) is one of op(x) = x, or op(x) = x', or op(x) = conjg(x'),
!!
!!  alpha and beta are scalars,
!!  A, B and C are matrices:
!!  op(A) is an m-by-k matrix,
!!  op(B) is a k-by-n matrix,
!!  C is an m-by-n matrix.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      lapackprof,m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zgemm(transa,transb,npws,ncola,ncolb,cg_a,cg_b,cg_c,alpha,beta)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zgemm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,ncola,ncolb
 real(dp),optional,intent(in) :: alpha(2),beta(2)
 character(len=1),intent(in) :: transa,transb
!arrays
 real(dp),intent(in) :: cg_a(2,npws*ncola)
 real(dp),intent(in) :: cg_b(2,npws*ncolb)
 real(dp),intent(inout) :: cg_c(2,*)

!Local variables-------------------------------
!scalars
 integer :: mm,nn,kk,lda,ldb,ldc
 real(dp) :: my_alpha(2),my_beta(2)

! *************************************************************************

 lda = npws
 ldb = npws

 mm  = npws
 nn  = ncolb
 kk  = ncola

 if (toupper(transa) /= 'N') then
   mm = ncola
   kk = npws
 end if
 if (toupper(transb) /= 'N') nn = npws

 ldc = mm

 my_alpha = cg_cone;  if (PRESENT(alpha)) my_alpha = alpha
 my_beta  = cg_czero; if (PRESENT(beta))  my_beta  = beta

 call ZGEMM(transa,transb,mm,nn,kk,my_alpha,cg_a,lda,cg_b,ldb,my_beta,cg_c,ldc)

end subroutine cg_zgemm
!!***


!!****f* m_cgtools/set_istwfk
!! NAME
!!  set_istwfk
!!
!! FUNCTION
!!  Returns the value of istwfk associated to the input k-point.
!!
!! INPUTS
!!  kpoint(3)=The k-point in reduced coordinates.
!!
!! OUTPUT
!!  istwfk= Integer flag internally used in the code to define the storage mode of the wavefunctions.
!!  It also define the algorithm used to apply an operator in reciprocal space as well as the FFT
!!  algorithm used to go from G- to r-space and vice versa.
!!
!!   1 => time-reversal cannot be used
!!   2 => use time-reversal at the Gamma point.
!!   3 => use time-reversal symmetry for k=(1/2, 0 , 0 )
!!   4 => use time-reversal symmetry for k=( 0 , 0 ,1/2)
!!   5 => use time-reversal symmetry for k=(1/2, 0 ,1/2)
!!   6 => use time-reversal symmetry for k=( 0 ,1/2, 0 )
!!   7 => use time-reversal symmetry for k=(1/2,1/2, 0 )
!!   8 => use time-reversal symmetry for k=( 0 ,1/2,1/2)
!!   9 => use time-reversal symmetry for k=(1/2,1/2,1/2)
!!
!!  Useful relations:
!!   u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^*
!!  and therefore:
!!   u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function set_istwfk(kpoint) result(istwfk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'set_istwfk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer :: bit0,ii
!arrays
 integer :: bit(3)

! *************************************************************************

 bit0=1

 do ii=1,3
   if (DABS(kpoint(ii))<tol10) then
     bit(ii)=0
   else if (DABS(kpoint(ii)-half)<tol10 ) then
     bit(ii)=1
   else
     bit0=0
   end if
 end do

 if (bit0==0) then
   istwfk=1
 else
   istwfk=2+bit(1)+4*bit(2)+2*bit(3) ! Note the inversion between bit(2) and bit(3)
 end if

end function set_istwfk
!!***

!!****f* m_cgtools/sqnorm_g
!! NAME
!! sqnorm_g
!!
!! FUNCTION
!! Compute the square of the norm of one complex vector vecti, in reciprocal space
!! Take into account the storage mode of the vector (istwf_k)
!!
!! INPUTS
!!  istwf_k=option parameter that describes the storage of wfs
!!  npwsp= (effective) number of planewaves at this k point.
!!  vect(2,npwsp)=the vector in reciprocal space (npw*nspinor, usually)
!!  me_g0=1 if this processors treats G=0, 0 otherwise.
!!  comm=MPI communicator for MPI sum.
!!
!! OUTPUT
!!  dotr= <vect|vect>
!!
!! PARENTS
!!      cgwf,dfpt_cgwf,dfpt_vtowfk,m_epjdos,mkresi,rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine sqnorm_g(dotr,istwf_k,npwsp,vect,me_g0,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sqnorm_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,npwsp,me_g0,comm
 real(dp),intent(out) :: dotr
!arrays
 real(dp),intent(in) :: vect(2,npwsp)

!Local variables-------------------------------
!scalars
 integer :: ierr

! *************************************************************************

 if (istwf_k==1) then ! General k-point
   !dotr = cg_real_zdotc(npwsp,vect,vect)
   dotr = cg_dznrm2(npwsp, vect)
   dotr = dotr * dotr

 else
   if (istwf_k==2 .and. me_g0==1) then
     ! Gamma k-point and I have G=0
     dotr=half*vect(1,1)**2
     dotr = dotr + cg_real_zdotc(npwsp-1,vect(1,2),vect(1,2))
   else
     ! Other TR k-points
     dotr = cg_real_zdotc(npwsp,vect,vect)
   end if
   dotr=two*dotr
 end if

 if (xmpi_comm_size(comm)>1) then
   call xmpi_sum(dotr,comm,ierr)
 end if

end subroutine sqnorm_g
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/dotprod_g
!! NAME
!! dotprod_g
!!
!! FUNCTION
!! Compute scalar product <vec1|vect2> of complex vectors vect1 and vect2 (can be the same)
!! Take into account the storage mode of the vectors (istwf_k)
!! If option=1, compute only real part, if option=2 compute also imaginary part.
!! If the number of calls to the dot product scales quadratically
!! with the volume of the system, it is preferable not to
!! call the present routine, but but to write a specially
!! optimized routine, that will avoid many branches related to
!! the existence of 'option' and 'istwf_k'.
!!
!! INPUTS
!!  istwf_k=option parameter that describes the storage of wfs
!!  vect1(2,npw)=first vector (one should take its complex conjugate)
!!  vect2(2,npw)=second vector
!!  npw= (effective) number of planewaves at this k point (including spinorial level)
!!  option= 1 if only real part to be computed,
!!          2 if both real and imaginary.
!!          3 if in case istwf_k==1 must compute real and imaginary parts,
!!               but if  istwf_k >1 must compute only real part
!!  me_g0=1 if this processor treats G=0, 0 otherwise
!!  comm=MPI communicator used to reduce the results.
!!
!! OUTPUT
!!  $doti=\Im ( <vect1|vect2> )$ , output only if option=2 and eventually option=3.
!!  $dotr=\Re ( <vect1|vect2> )$
!!
!! PARENTS
!!      cgwf,chebfi,corrmetalwf1,d2frnl,dfpt_cgwf,dfpt_nsteltwf,dfpt_nstpaw
!!      dfpt_nstwf,dfpt_vtowfk,dfpt_wfkfermi,dfptnl_resp,dotprod_set_cgcprj
!!      dotprodm_sumdiag_cgcprj,eig2stern,extrapwf,fock2ACE,fock_ACE_getghc
!!      fock_getghc,m_efmas,m_gkk,m_phgamma,m_phpi,m_rf2,m_sigmaph,mkresi
!!      nonlop_gpu,nonlop_test,rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine dotprod_g(dotr,doti,istwf_k,npw,option,vect1,vect2,me_g0,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprod_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,npw,option,me_g0,comm
 real(dp),intent(out) :: doti,dotr
!arrays
 real(dp),intent(in) :: vect1(2,npw),vect2(2,npw)

!Local variables-------------------------------
!scalars
 integer :: ierr
!arrays
 real(dp) :: dotarr(2)

! *************************************************************************

 if (istwf_k==1) then ! General k-point

   if(option==1)then
     dotr = cg_real_zdotc(npw,vect1,vect2)
   else
     dotarr = cg_zdotc(npw,vect1,vect2)
     dotr = dotarr(1)
     doti = dotarr(2)
   end if

 else if (istwf_k==2 .and. me_g0==1) then ! Gamma k-point and I have G=0

   dotr=half*vect1(1,1)*vect2(1,1)
   dotr = dotr + cg_real_zdotc(npw-1,vect1(1,2),vect2(1,2))
   dotr = two*dotr
   if (option==2) doti=zero

 else ! Other TR k-points

   dotr = cg_real_zdotc(npw,vect1,vect2)
   dotr=two*dotr
   if (option==2) doti=zero
 end if

!Reduction in case of parallelism
 if (xmpi_comm_size(comm)>1) then
   if (option==1.or.istwf_k/=1) then
     call xmpi_sum(dotr,comm,ierr)
   else
     dotarr(1)=dotr ; dotarr(2)=doti
     call xmpi_sum(dotarr,comm,ierr)
     dotr=dotarr(1) ; doti=dotarr(2)
   end if
 end if

end subroutine dotprod_g
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/matrixelmt_g
!! NAME
!! matrixelmt_g
!!
!! FUNCTION
!!  Compute a matrix element of two wavefunctions, in reciprocal space,
!!  for an operator that is diagonal in reciprocal space: <wf1|op|wf2>
!!  For the time being, only spin-independent operators are treated.
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent !)
!!  istwf_k=storage mode of the vectors
!!  needimag=0 if the imaginary part is not needed ; 1 if the imaginary part is needed
!!  npw=number of planewaves of the first vector
!!  nspinor=number of spinor components
!!  vect1(2,npw*nspinor)=first vector
!!  vect2(2,npw*nspinor)=second vector
!!  comm_fft=MPI communicator for the FFT
!!  me_g0=1 if this processors treats the G=0 component.
!!
!! OUTPUT
!!  ai=imaginary part of the matrix element
!!  ar=real part of the matrix element
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine matrixelmt_g(ai,ar,diag,istwf_k,needimag,npw,nspinor,vect1,vect2,me_g0,comm_fft)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matrixelmt_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,needimag,npw,nspinor,me_g0,comm_fft
 real(dp),intent(out) :: ai,ar
!arrays
 real(dp),intent(in) :: diag(npw),vect1(2,npw*nspinor),vect2(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: i1,ierr,ipw
 character(len=500) :: message
!arrays
 real(dp) :: buffer2(2)
 !real(dp),allocatable :: re_prod(:), im_prod(:)

! *************************************************************************

 if (nspinor==2 .and. istwf_k/=1) then
   write(message,'(a,a,a,i6,a,i6)')&
&   'When istwf_k/=1, nspinor must be 1,',ch10,&
&   'however, nspinor=',nspinor,', and istwf_k=',istwf_k
   MSG_BUG(message)
 end if

#if 0
 !TODO
 ABI_MALLOC(re_prod,(npw*nspinor))
 do ipw=1,npw*nspinor
  re_prod(ipw) = vect1(1,ipw)*vect2(1,ipw) + vect1(2,ipw)*vect2(2,ipw)
 end do

 if (needimag == 1) then
   ABI_MALLOC(im_prod,(npw*nspinor))
   do ipw=1,npw*nspinor
     im_prod(ipw) = vect1(1,ipw)*vect2(2,ipw) - vect1(2,ipw)*vect2(1,ipw)
   end do
 end if
#endif

 ar=zero
 if(needimag==1)ai=zero

!Normal storage mode
 if(istwf_k==1)then

!  Need only real part
   if(needimag==0)then

     do ipw=1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
     end do

     if(nspinor==2)then
       do ipw=1+npw,2*npw
         ar=ar+diag(ipw-npw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
       end do
     end if

   else ! Need also the imaginary part

     do ipw=1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
       ai=ai+diag(ipw)*(vect1(1,ipw)*vect2(2,ipw)-vect1(2,ipw)*vect2(1,ipw))
     end do

     if(nspinor==2)then
       do ipw=1+npw,2*npw
         ar=ar+diag(ipw-npw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
         ai=ai+diag(ipw-npw)*(vect1(1,ipw)*vect2(2,ipw)-vect1(2,ipw)*vect2(1,ipw))
       end do
     end if

   end if ! needimag

 else if(istwf_k>=2)then

!  XG030513 : MPIWF need to know which proc has G=0

   i1=1
   if(istwf_k==2 .and. me_g0==1)then
     ar=half*diag(1)*vect1(1,1)*vect2(1,1) ; i1=2
   end if

!  Need only real part
   if(needimag==0)then

     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
     end do
     ar=two*ar

   else ! Need also the imaginary part

     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
       ai=ai+diag(ipw)*(vect1(1,ipw)*vect2(2,ipw)-vect1(2,ipw)*vect2(1,ipw))
     end do
     ar=two*ar ; ai=two*ai

   end if

 end if ! istwf_k

#if 0
 ABI_FREE(re_prod)
 if (needimag == 1) then
   ABI_FREE(im_prod)
 end if
#endif

!MPIWF need to make reduction on ar and ai .
 if (xmpi_comm_size(comm_fft)>1) then
   buffer2(1)=ai
   buffer2(2)=ar
   call xmpi_sum(buffer2,comm_fft ,ierr)
   ai=buffer2(1)
   ar=buffer2(2)
 end if

end subroutine matrixelmt_g
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/dotprod_v
!! NAME
!! dotprod_v
!!
!! FUNCTION
!! Compute dot product of two potentials (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the potentials (nspden), and sum over them.
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potentials are stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potentials are stored as V, B_x, B_y, Bz  (B=magn. field)
!!  pot1(cplex*nfft,nspden)=first real space potential on FFT grid
!!  pot2(cplex*nfft,nspden)=second real space potential on FFT grid
!!  comm=MPI communicator in which results will be reduced.
!!
!! OUTPUT
!!  dotr= value of the dot product
!!
!! PARENTS
!!      m_epjdos
!!
!! CHILDREN
!!
!! SOURCE


subroutine dotprod_v(cplex,dotr,nfft,nspden,opt_storage,pot1,pot2,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprod_v'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,opt_storage,comm
 real(dp),intent(out) :: dotr
!arrays
 real(dp),intent(in) :: pot1(cplex*nfft,nspden),pot2(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ispden
 real(dp) :: ar
!arrays

! *************************************************************************

!Real or complex inputs are coded

 dotr=zero
!$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:dotr)
 do ispden=1,min(nspden,2)
   do ifft=1,cplex*nfft
     dotr =dotr + pot1(ifft,ispden)*pot2(ifft,ispden)
   end do
 end do

 if (nspden==4) then
   ar=zero
!$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:ar)
   do ispden=3,4
     do ifft=1,cplex*nfft
       ar = ar + pot1(ifft,ispden)*pot2(ifft,ispden)
     end do
   end do

   if (opt_storage==0) then
     if (cplex==1) then
       dotr = dotr+two*ar
     else
       dotr = dotr+ar
     end if
   else
     dotr = half*(dotr+ar)
   end if
 end if

!MPIWF reduction (addition) on dotr is needed here
 if (xmpi_comm_size(comm)>1) then
   call xmpi_sum(dotr,comm,ierr)
 end if

end subroutine dotprod_v
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/sqnorm_v
!! NAME
!! sqnorm_v
!!
!! FUNCTION
!! Compute square of the norm of a potential (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the potentials (nspden),
!! and sum over them.
!!
!! INPUTS
!!  cplex=if 1, real space function on FFT grid is REAL, if 2, COMPLEX
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potential is stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potential is stored as V, B_x, B_y, Bz  (B=magn. field)
!!  pot(cplex*nfft,nspden)=real space potential on FFT grid
!!
!! OUTPUT
!!  norm2= value of the square of the norm
!!
!! PARENTS
!!      dfpt_rhotov,dfpt_vtorho,rhotov,vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine sqnorm_v(cplex,nfft,norm2,nspden,opt_storage,pot,mpi_comm_sphgrid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sqnorm_v'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,opt_storage
 integer,intent(in),optional :: mpi_comm_sphgrid
 real(dp),intent(out) :: norm2
!arrays
 real(dp),intent(in) :: pot(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ispden,nproc_sphgrid
 real(dp) :: ar

! *************************************************************************

!Real or complex inputs are coded

 norm2=zero
 do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ispden,nfft,pot) REDUCTION(+:norm2)
   do ifft=1,cplex*nfft
     norm2=norm2 + pot(ifft,ispden)**2
   end do
 end do
 if (nspden==4) then
   ar=zero
   do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ispden,nfft,pot) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
       ar=ar + pot(ifft,ispden)**2
     end do
   end do
   if (opt_storage==0) then
     if (cplex==1) then
       norm2=norm2+two*ar
     else
       norm2=norm2+ar
     end if
   else
     norm2=half*(norm2+ar)
   end if
 end if

!MPIWF reduction (addition) on norm2 is needed here
 if(present(mpi_comm_sphgrid)) then
   nproc_sphgrid=xmpi_comm_size(mpi_comm_sphgrid)
   if(nproc_sphgrid>1)then
     call xmpi_sum(norm2,mpi_comm_sphgrid,ierr)
   end if
 end if

end subroutine sqnorm_v
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/mean_fftr
!! NAME
!! mean_fftr
!!
!! FUNCTION
!!  Compute the mean of an arraysp(nfft,nspden), over the FFT grid, for each component nspden,
!!  and return it in meansp(nspden).
!!  Take into account the spread of the array due to parallelism: the actual number of fft
!!  points is nfftot, but the number of points on this proc is nfft only.
!!  So : for ispden from 1 to nspden
!!       meansp(ispden) = sum(ifft=1,nfftot) arraysp(ifft,ispden) / nfftot
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  nspden=number of spin-density components
!!
!! OUTPUT
!!  meansp(nspden)=mean value for each nspden component
!!
!! PARENTS
!!      fresid,newvtr,pawmknhat,prcref,prcref_PMA,psolver_rhohxc,rhohxcpositron
!!      rhotov,rhotoxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine mean_fftr(arraysp,meansp,nfft,nfftot,nspden,mpi_comm_sphgrid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mean_fftr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nfftot,nspden
 integer,intent(in),optional:: mpi_comm_sphgrid
!arrays
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(out) :: meansp(nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ispden,nproc_sphgrid
 real(dp) :: invnfftot,tmean

! *************************************************************************

 invnfftot=one/(dble(nfftot))

 do ispden=1,nspden
   tmean=zero
!$OMP PARALLEL DO REDUCTION(+:tmean)
   do ifft=1,nfft
     tmean=tmean+arraysp(ifft,ispden)
   end do
   meansp(ispden)=tmean*invnfftot
 end do

!XG030514 : MPIWF The values of meansp(ispden) should
!now be summed accross processors in the same WF group, and spread on all procs.
 if(present(mpi_comm_sphgrid)) then
   nproc_sphgrid=xmpi_comm_size(mpi_comm_sphgrid)
   if(nproc_sphgrid>1) then
     call xmpi_sum(meansp,nspden,mpi_comm_sphgrid,ierr)
   end if
 end if

end subroutine mean_fftr
!!***

!!****f* m_cgtools/cg_getspin
!! NAME
!! cg_getspin
!!
!! FUNCTION
!!  Sandwich a single wave function on the Pauli matrices
!!
!! INPUTS
!!  npw_k = number of plane waves
!!  cgcband = coefficients of spinorial wave function
!!
!! OUTPUT
!!  spin = 3-vector of spin components for this state
!!  cgcmat = outer spin product of spinorial wf with itself
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!
!! SOURCE

subroutine  cg_getspin(cgcband, npw_k, spin, cgcmat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_getspin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: npw_k
 real(dp), intent(in) :: cgcband(2,2*npw_k)
 complex(dpc), intent(out),optional :: cgcmat(2,2)
 real(dp), intent(out) :: spin(3)

!Local variables-------------------------------
!scalars
 complex(dpc),parameter :: pauli_0(2,2) = reshape([cone,czero,czero,cone], [2,2])
 complex(dpc),parameter :: pauli_x(2,2) = reshape([czero,cone,cone,czero], [2,2])
 complex(dpc),parameter :: pauli_y(2,2) = reshape([czero,j_dpc,-j_dpc,czero], [2,2])
 complex(dpc),parameter :: pauli_z(2,2) = reshape([cone,czero,czero,-cone], [2,2])
 complex(dpc) :: cspin(0:3), cgcmat_(2,2)
! ***********************************************************************

! cgcmat_ = cgcband * cgcband^T*  i.e. 2x2 matrix of spin components (dpcomplex)
 cgcmat_ = czero
 call zgemm('n','c',2,2,npw_k,cone,cgcband,2,cgcband,2,czero,cgcmat_,2)

! spin(*)  = sum_{si sj pi} cgcband(si,pi)^* pauli_*(si,sj) cgcband(sj,pi)
 cspin(0) = cgcmat_(1,1)*pauli_0(1,1) + cgcmat_(2,1)*pauli_0(2,1) &
& + cgcmat_(1,2)*pauli_0(1,2) + cgcmat_(2,2)*pauli_0(2,2)
 cspin(1) = cgcmat_(1,1)*pauli_x(1,1) + cgcmat_(2,1)*pauli_x(2,1) &
& + cgcmat_(1,2)*pauli_x(1,2) + cgcmat_(2,2)*pauli_x(2,2)
 cspin(2) = cgcmat_(1,1)*pauli_y(1,1) + cgcmat_(2,1)*pauli_y(2,1) &
& + cgcmat_(1,2)*pauli_y(1,2) + cgcmat_(2,2)*pauli_y(2,2)
 cspin(3) = cgcmat_(1,1)*pauli_z(1,1) + cgcmat_(2,1)*pauli_z(2,1) &
& + cgcmat_(1,2)*pauli_z(1,2) + cgcmat_(2,2)*pauli_z(2,2)
!write(std_out,*) 'cgmat: ', cgcmat_
!write(std_out,*) 'real(spin): ', real(cspin)
!write(std_out,*) 'aimag(spin): ', aimag(cspin)

 spin = real(cspin(1:3))
 if (present(cgcmat)) cgcmat = cgcmat_

end subroutine cg_getspin
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_gsph2box
!! NAME
!! cg_gsph2box
!!
!! FUNCTION
!! Array iarrsph is defined in sphere with npw_k points. Insert iarrsph inside box
!! of nx*ny*nz points to define array oarrbox for fft box. rest of oarrbox is filled with 0 s.
!!
!! INPUTS
!! iarrsph(2,npw_k*ndat)= contains values for npw_k G vectors in basis sphere
!! ndat=number of FFT to perform.
!! npw_k=number of G vectors in basis at this k point
!! oarrbox(2,ldx*ldy*ldz*ndat) = fft box
!! nx,ny,nz=physical dimension of the box (oarrbox)
!! ldx,ldy,ldz=memory dimension of oarrbox
!! kg_k(3,npw_k)=integer coordinates of G vectors in basis sphere
!! istwf_k=option parameter that describes the storage of wfs
!!
!! OUTPUT
!!   oarrbox(ldx*ldy*ldz*ndat)
!!
!! NOTES
!! If istwf_k differs from 1, then special storage modes must be taken
!! into account, for symmetric wavefunctions coming from k=(0 0 0) or other
!! special k points.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_gsph2box(nx,ny,nz,ldx,ldy,ldz,ndat,npw_k,istwf_k,kg_k,iarrsph,oarrbox)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_gsph2box'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,nx,ny,nz,ldx,ldy,ldz,ndat,npw_k
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: iarrsph(2,npw_k*ndat)
 real(dp),intent(out) :: oarrbox(2,ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: me_g0=1
 integer :: ix,ixinv,iy,iyinv,iz,izinv,dat,ipw,npwmin,pad_box,pad_sph,ifft,ifft_inv,ldxyz
 character(len=500) :: msg
!arrays
 integer,allocatable :: ixinver(:),iyinver(:),izinver(:)

! *************************************************************************

!In the case of special k-points, invariant under time-reversal,
!but not Gamma, initialize the inverse coordinates
!Remember indeed that
!u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^*
!and therefore:
!u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.
 if (istwf_k>=2) then
   ABI_MALLOC(ixinver,(nx))
   ABI_MALLOC(iyinver,(ny))
   ABI_MALLOC(izinver,(nz))
   if ( ANY(istwf_k==(/2,4,6,8/)) ) then
     ixinver(1)=1
     do ix=2,nx
       ixinver(ix)=nx+2-ix
     end do
   else
     do ix=1,nx
       ixinver(ix)=nx+1-ix
     end do
   end if
   if (istwf_k>=2 .and. istwf_k<=5) then
     iyinver(1)=1
     do iy=2,ny
       iyinver(iy)=ny+2-iy
     end do
   else
     do iy=1,ny
       iyinver(iy)=ny+1-iy
     end do
   end if
   if ( ANY(istwf_k==(/2,3,6,7/)) ) then
     izinver(1)=1
     do iz=2,nz
       izinver(iz)=nz+2-iz
     end do
   else
     do iz=1,nz
       izinver(iz)=nz+1-iz
     end do
   end if
 end if

 ldxyz = ldx*ldy*ldz

 if (istwf_k==1) then

!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,iy,iz,ifft)
   do dat=1,ndat
     pad_sph = (dat-1)*npw_k
     pad_box = (dat-1)*ldxyz
     oarrbox(:,1+pad_box:ldxyz+pad_box) = zero ! zero the sub-array
     do ipw=1,npw_k
       ix=kg_k(1,ipw); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ipw); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ipw); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
#if defined __INTEL_COMPILER && defined HAVE_OPENMP
       if (ifft==0) then
         MSG_ERROR("prevent ifort+OMP from miscompiling this section on cronos")
       end if
#endif
       oarrbox(1,ifft+pad_box) = iarrsph(1,ipw+pad_sph)
       oarrbox(2,ifft+pad_box) = iarrsph(2,ipw+pad_sph)
     end do
   end do

 else if (istwf_k>=2) then
   !
   npwmin=1
   if(istwf_k==2 .and. me_g0==1) then ! If gamma point, then oarrbox must be completed
     do dat=1,ndat
       pad_sph = (dat-1)*npw_k
       pad_box = (dat-1)*ldxyz
       oarrbox(1,1+pad_box) = iarrsph(1,1+pad_sph)
       oarrbox(2,1+pad_box) = zero
     end do
     npwmin=2
   end if

!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,iy,iz,ixinv,iyinv,izinv,ifft,ifft_inv)
   do dat=1,ndat
     pad_sph = (dat-1)*npw_k
     pad_box = (dat-1)*ldxyz
     oarrbox(:,npwmin+pad_box:ldxyz+pad_box) = zero
     do ipw=npwmin,npw_k
       ix=kg_k(1,ipw); if(ix<0)ix=ix+nx; ix=ix+1
       iy=kg_k(2,ipw); if(iy<0)iy=iy+ny; iy=iy+1
       iz=kg_k(3,ipw); if(iz<0)iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       ! Construct the coordinates of -k-G
       ixinv=ixinver(ix); iyinv=iyinver(iy); izinv=izinver(iz)
       ifft_inv = ixinv + (iyinv-1)*ldx + (izinv-1)*ldx*ldy

       oarrbox(:,ifft    +pad_box) =  iarrsph(:,ipw+pad_sph)
       oarrbox(1,ifft_inv+pad_box) =  iarrsph(1,ipw+pad_sph)
       oarrbox(2,ifft_inv+pad_box) = -iarrsph(2,ipw+pad_sph)
     end do
   end do
   !
 else
   write(msg,'(a,i0)')"Wrong istwfk ",istwf_k
   MSG_ERROR(msg)
 end if

 if (istwf_k>=2) then
   ABI_FREE(ixinver)
   ABI_FREE(iyinver)
   ABI_FREE(izinver)
 end if

end subroutine cg_gsph2box
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_box2gsph
!! NAME
!!  cg_box2gsph
!!
!! FUNCTION
!!
!! INPUTS
!!  nx,ny,nz=physical dimension of the FFT box.
!!  ldx,ldy,ldz=Logical dimensions of the arrays.
!!  ndat=number of data in iarrbox
!!  npw_k=Number of planewaves in the G-sphere.
!!  kg_k(3,npw_k)=Reduced coordinates of the G-vectoes.
!!  iarrbox(2,ldx,ldy,ldz*ndat)=Input arrays on the FFT box.
!!  [rscal] = Scaling factor
!!
!! OUTPUT
!!  oarrsph(2,npw_k*ndat)=Data defined on the G-sphere.
!!
!! PARENTS
!!      fourwf,m_dfti,m_fftw3
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_box2gsph(nx,ny,nz,ldx,ldy,ldz,ndat,npw_k,kg_k,iarrbox,oarrsph,rscal)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_box2gsph'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat
 real(dp),optional,intent(in) :: rscal
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: iarrbox(2,ldx*ldy*ldz*ndat)
 real(dp),intent(out) :: oarrsph(2,npw_k*ndat)

!Local variables-------------------------------
!scalars
 integer :: ig,ix,iy,iz,idat,sph_pad,box_pad,ifft

! *************************************************************************

 if (.not. PRESENT(rscal)) then
   !
   if (ndat==1) then
!$OMP PARALLEL DO PRIVATE(ix,iy,iz,ifft)
     do ig=1,npw_k
       ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       oarrsph(1,ig) = iarrbox(1,ifft)
       oarrsph(2,ig) = iarrbox(2,ifft)
     end do
   else
!$OMP PARALLEL DO PRIVATE(sph_pad,box_pad,ix,iy,iz,ifft)
     do idat=1,ndat
       sph_pad = (idat-1)*npw_k
       box_pad = (idat-1)*ldx*ldy*ldz
       do ig=1,npw_k
         ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
         iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
         iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
         ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
         oarrsph(1,ig+sph_pad) = iarrbox(1,ifft+box_pad)
         oarrsph(2,ig+sph_pad) = iarrbox(2,ifft+box_pad)
       end do
     end do
   end if
   !
 else
   if (ndat==1) then
!$OMP PARALLEL DO PRIVATE(ix,iy,iz,ifft)
     do ig=1,npw_k
       ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       oarrsph(1,ig) = iarrbox(1,ifft) * rscal
       oarrsph(2,ig) = iarrbox(2,ifft) * rscal
     end do
   else
!$OMP PARALLEL DO PRIVATE(sph_pad,box_pad,ix,iy,iz,ifft)
     do idat=1,ndat
       sph_pad = (idat-1)*npw_k
       box_pad = (idat-1)*ldx*ldy*ldz
       do ig=1,npw_k
         ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
         iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
         iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
         ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
         oarrsph(1,ig+sph_pad) = iarrbox(1,ifft+box_pad) * rscal
         oarrsph(2,ig+sph_pad) = iarrbox(2,ifft+box_pad) * rscal
       end do
     end do
   end if
 end if

end subroutine cg_box2gsph
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_addtorho
!! NAME
!!  cg_addtorho
!!
!! FUNCTION
!!  Add |ur|**2 to the ground-states density rho.
!!    rho = rho + weight_r * Re[ur]**2 + weight_i * Im[ur]**2
!!
!! INPUTS
!!  nx,ny,nz=physical dimension of the FFT box.
!!  ldx,ldy,ldz=leading dimensions of the arrays.
!!  ndat=number of contributions to accumulate.
!!  weight_r=weight used for the accumulation of the density in real space
!!  weight_i=weight used for the accumulation of the density in real space
!!  ur(2,ldx,ldy,ldz*ndat)=wavefunctions in real space
!!
!! SIDE EFFECTS
!!  rho(ldx,ldy,ldz) = contains the input density at input,
!!                  modified in input with the contribution gived by ur.
!!
!! PARENTS
!!      fourwf,m_dfti,m_fftw3
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_addtorho(nx,ny,nz,ldx,ldy,ldz,ndat,weight_r,weight_i,ur,rho)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_addtorho'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
 real(dp),intent(in) :: weight_i,weight_r
!arrays
 real(dp),intent(in) :: ur(2,ldx,ldy,ldz*ndat)
 real(dp),intent(inout) :: rho(ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,idat,izdat

! *************************************************************************

 if (ndat==1) then
!$OMP PARALLEL DO
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         rho(ix,iy,iz) = rho(ix,iy,iz) + weight_r * ur(1,ix,iy,iz)**2 &
&                                      + weight_i * ur(2,ix,iy,iz)**2
       end do
     end do
   end do

 else
! It would be nice to use $OMP PARALLEL DO PRIVATE(izdat) REDUCTION(+:rho)
! but it's risky as the private rho is allocated on the stack of the thread.
!$OMP PARALLEL PRIVATE(izdat)
   do idat=1,ndat
!$OMP DO
     do iz=1,nz
       izdat = iz + (idat-1)*ldz
       do iy=1,ny
         do ix=1,nx
           rho(ix,iy,iz) = rho(ix,iy,iz) + weight_r * ur(1,ix,iy,izdat)**2 &
&                                        + weight_i * ur(2,ix,iy,izdat)**2
         end do
       end do
     end do
!$OMP END DO NOWAIT
   end do
!$OMP END PARALLEL
 end if

end subroutine cg_addtorho
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_vlocpsi
!! NAME
!!  cg_vlocpsi
!!
!! FUNCTION
!!  Apply the local part of the potentatil to the wavefunction in real space.
!!
!! INPUTS
!!  nx,ny,nz=physical dimension of the FFT box.
!!  ldx,ldy,ldz=leading dimensions of the arrays.
!!  ndat=number of wavefunctions.
!!  cplex=  1 if vloc is real, 2 for complex
!!  vloc(cplex*ldx,ldy,ldz)=Local potential on the FFT box.
!!
!! SIDE EFFECTS
!!  ur(2,ldx,ldy,ldz*ndat)=
!!    Input = wavefunctions in real space.
!!    Output= vloc |ur>
!!
!! PARENTS
!!      m_dfti,m_fftw3
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_vlocpsi(nx,ny,nz,ldx,ldy,ldz,ndat,cplex,vloc,ur)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_vlocpsi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,cplex
!arrays
 real(dp),intent(in) :: vloc(cplex*ldx,ldy,ldz)
 real(dp),intent(inout) :: ur(2,ldx,ldy,ldz*ndat)

!Local variables-------------------------------
!scalars
 integer :: idat,ix,iy,iz,padat
 real(dp) :: fim,fre

! *************************************************************************

 if (cplex==1) then
   !
   if (ndat==1) then
!$OMP PARALLEL DO
     do iz=1,nz
       do iy=1,ny
         do ix=1,nx
           ur(1,ix,iy,iz) = vloc(ix,iy,iz) * ur(1,ix,iy,iz)
           ur(2,ix,iy,iz) = vloc(ix,iy,iz) * ur(2,ix,iy,iz)
         end do
       end do
     end do
     !
   else
     !
!$OMP PARALLEL DO PRIVATE(padat)
     do idat=1,ndat
       padat = ldz*(idat-1)
       do iz=1,nz
         do iy=1,ny
           do ix=1,nx
             ur(1,ix,iy,iz+padat) = vloc(ix,iy,iz) * ur(1,ix,iy,iz+padat)
             ur(2,ix,iy,iz+padat) = vloc(ix,iy,iz) * ur(2,ix,iy,iz+padat)
           end do
         end do
       end do
     end do
     !
   end if
   !
 else if (cplex==2)then
   !
   if (ndat==1) then
!$OMP PARALLEL DO PRIVATE(fre,fim)
     do iz=1,nz
       do iy=1,ny
         do ix=1,nx
           fre = ur(1,ix,iy,iz)
           fim = ur(2,ix,iy,iz)
           ur(1,ix,iy,iz) = vloc(2*ix-1,iy,iz)*fre - vloc(2*ix,iy,iz)*fim
           ur(2,ix,iy,iz) = vloc(2*ix-1,iy,iz)*fim + vloc(2*ix,iy,iz)*fre
         end do
       end do
     end do
   else
!$OMP PARALLEL DO PRIVATE(padat,fre,fim)
     do idat=1,ndat
       padat = ldz*(idat-1)
       do iz=1,nz
         do iy=1,ny
           do ix=1,nx
             fre = ur(1,ix,iy,iz+padat)
             fim = ur(2,ix,iy,iz+padat)
             ur(1,ix,iy,iz+padat) = vloc(2*ix-1,iy,iz)*fre - vloc(2*ix,iy,iz)*fim
             ur(2,ix,iy,iz+padat) = vloc(2*ix-1,iy,iz)*fim + vloc(2*ix,iy,iz)*fre
           end do
         end do
       end do
     end do
   end if
   !
 else
   ur = huge(one)
   !MSG_BUG("Wrong cplex")
 end if

end subroutine cg_vlocpsi
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgnc_cholesky
!! NAME
!!  cgnc_cholesky
!!
!! FUNCTION
!!  Cholesky orthonormalization of the vectors stored in cgblock (version optimized for NC wavefunctions).
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband=Number of band in cgblock
!!  istwfk=Storage mode for the wavefunctions. 1 for standard full mode
!!  me_g0=1 if this node has G=0.
!!  comm_pw=MPI communicator for the planewave group. Set to xmpi_comm_self for sequential mode.
!!
!! SIDE EFFECTS
!!  cgblock(2*npws*nband)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!      pw_orthon
!!
!! CHILDREN
!!
!! SOURCE

subroutine cgnc_cholesky(npws,nband,cgblock,istwfk,me_g0,comm_pw,use_gemm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgnc_cholesky'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband,istwfk
 integer,intent(in) :: comm_pw,me_g0
 logical,optional,intent(in) :: use_gemm
!arrays
 real(dp),intent(inout) :: cgblock(2*npws*nband)

!Local variables ------------------------------
!scalars
 integer :: ierr,b1,b2
#ifdef DEBUG_MODE
 integer :: ptr
#endif
 logical :: my_usegemm
 character(len=500) :: msg
!arrays
 real(dp) :: rcg0(nband)
 real(dp),allocatable :: rovlp(:,:)
 complex(dpc),allocatable :: cf_ovlp(:,:)

! *************************************************************************

#ifdef DEBUG_MODE
 if (istwfk==2) then
   ierr = 0
   do b1=1,nband
     ptr = 2 + 2*(b1-1)*npws
     if (ABS(cgblock(ptr)) > zero) then
       ierr = ierr + 1
       write(msg,'(a,i0,es13.6)')" Input b1, Im u(g=0) should be zero ",b1,cgblock(ptr)
       call wrtout(std_out,msg,"COLL")
       !cgblock(ptr) = zero
     end if
   end do
   ABI_CHECK(ierr==0,"Non zero imag part")
 end if
#endif

 my_usegemm=.FALSE.; if (PRESENT(use_gemm)) my_usegemm = use_gemm

 ABI_MALLOC(cf_ovlp,(nband,nband))
 !
 ! 1) Calculate O_ij = <phi_i|phi_j>
 if (my_usegemm) then
   call ABI_ZGEMM("Conjugate","Normal",nband,nband,npws,cone,cgblock,npws,cgblock,npws,czero,cf_ovlp,nband)
 else
   call ZHERK("U","C",nband,npws,one,cgblock,npws,zero,cf_ovlp,nband)
 end if

 if (istwfk==1) then
   !
   ! Sum the overlap if PW are distributed.
   if (comm_pw /= xmpi_comm_self) then
     call xmpi_sum(cf_ovlp,comm_pw,ierr)
   end if
   !
   ! 2) Cholesky factorization: O = U^H U with U upper triangle matrix.
   call ZPOTRF('U',nband,cf_ovlp,nband,ierr)

   if (ierr/=0)  then
     write(msg,'(a,i0)')' ZPOTRF returned info = ',ierr
     MSG_ERROR(msg)
   end if

 else
   ! overlap is real. Note that nspinor is always 1 in this case.
   ABI_MALLOC(rovlp,(nband,nband))
   rovlp = two * REAL(cf_ovlp)

   if (istwfk==2 .and. me_g0==1) then
     ! Extract the real part at G=0 and subtract its contribution to the overlap.
     call dcopy(nband,cgblock,2*npws,rcg0,1)
     do b2=1,nband
       do b1=1,b2
         rovlp(b1,b2) = rovlp(b1,b2) - rcg0(b1)*rcg0(b2)
       end do
     end do
   end if

   ! Sum the overlap if PW are distributed.
   if (comm_pw /= xmpi_comm_self) then
     call xmpi_sum(rovlp,comm_pw,ierr)
   end if
   !
   ! 2) Cholesky factorization: O = U^H U with U upper triangle matrix.
   call DPOTRF('U',nband,rovlp,nband,ierr)

   if (ierr/=0)  then
     write(msg,'(a,i0)')' DPOTRF returned info = ',ierr
     MSG_ERROR(msg)
   end if

   cf_ovlp = DCMPLX(rovlp)
   ABI_FREE(rovlp)
 end if
 !
 ! 3) Solve X U = cgblock. On exit cgblock is orthonormalized.
 call ZTRSM('Right','Upper','Normal','Normal',npws,nband,cone,cf_ovlp,nband,cgblock,npws)

#ifdef DEBUG_MODE
 if (istwfk==2) then
   ierr = 0
   do b1=1,nband
     ptr = 2 + 2*(b1-1)*npws
     if (ABS(cgblock(ptr)) > zero) then
       ierr = ierr + 1
       write(msg,'(a,i0,es13.6)')" Output b1, Im u(g=0) should be zero ",b1,cgblock(ptr)
     end if
   end do
   ABI_CHECK(ierr==0,"Non zero imag part")
 end if
#endif

 ABI_FREE(cf_ovlp)

end subroutine cgnc_cholesky
!!***

!----------------------------------------------------------------------

!!****f* m_blas/cgpaw_cholesky
!! NAME
!!  cgpaw_cholesky
!!
!! FUNCTION
!!  Cholesky orthonormalization of the vectors stored in cgblock. (version for PAW wavefunctions).
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband=Number of band in cgblock and gsc
!!  istwfk=Storage mode for the wavefunctions. 1 for standard full mode
!!  me_g0=1 if this node has G=0.
!!  comm_pw=MPI communicator for the planewave group. Set to xmpi_comm_self for sequential mode.
!!
!! SIDE EFFECTS
!!  cgblock(2*npws*nband)
!!    input: Input set of vectors |C>, S|C>
!!    output: Orthonormalized set such as  <C|S|C> = 1
!!  gsc(2*npws*nband)
!!   destroyed in output.
!!
!! PARENTS
!!      pw_orthon
!!
!! CHILDREN
!!
!! SOURCE

subroutine cgpaw_cholesky(npws,nband,cgblock,gsc,istwfk,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgpaw_cholesky'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband,istwfk
 integer,intent(in) :: me_g0,comm_pw
!arrays
 real(dp),intent(inout) :: cgblock(2*npws*nband)
 real(dp),intent(inout) :: gsc(2*npws*nband)

!Local variables ------------------------------
!scalars
 integer :: ierr,b1,b2
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rovlp(:,:)
 real(dp) :: rcg0(nband),rg0sc(nband)
 complex(dpc),allocatable :: cf_ovlp(:,:)

! *************************************************************************

 ! 1) Calculate O_ij =  <phi_i|S|phi_j>
 ABI_MALLOC(cf_ovlp,(nband,nband))

 call ABI_ZGEMM("C","N",nband,nband,npws,cone,cgblock,npws,gsc,npws,czero,cf_ovlp,nband)

 if (istwfk==1) then
   !
   ! Sum the overlap if PW are distributed.
   if (comm_pw /= xmpi_comm_self) then
     call xmpi_sum(cf_ovlp,comm_pw,ierr)
   end if
   !
   ! 2) Cholesky factorization: O = U^H U with U upper triangle matrix.
   call ZPOTRF('U',nband,cf_ovlp,nband,ierr)

   if (ierr/=0)  then
     write(msg,'(a,i0)')' ZPOTRF returned info= ',ierr
     MSG_ERROR(msg)
   end if

 else
   ! overlap is real. Note that nspinor is always 1 in this case.
   ABI_MALLOC(rovlp,(nband,nband))
   rovlp = two * REAL(cf_ovlp)

   if (istwfk==2 .and. me_g0==1) then
     ! Extract the real part at G=0 and subtract its contribution to the overlap.
     call dcopy(nband,cgblock,2*npws,rcg0,1)
     call dcopy(nband,gsc,2*npws,rg0sc,1)
     do b2=1,nband
       do b1=1,b2
        rovlp(b1,b2) = rovlp(b1,b2) - rcg0(b1)*rg0sc(b2)
       end do
     end do
   end if
   !
   ! Sum the overlap if PW are distributed.
   if (comm_pw /= xmpi_comm_self) then
     call xmpi_sum(rovlp,comm_pw,ierr)
   end if
  !
   ! 2) Cholesky factorization: O = U^H U with U upper triangle matrix.
   call DPOTRF('U',nband,rovlp,nband,ierr)

   if (ierr/=0)  then
     write(msg,'(a,i0)')' DPOTRF returned info= ',ierr
     MSG_ERROR(msg)
   end if

   cf_ovlp = DCMPLX(rovlp)
   ABI_FREE(rovlp)
 end if
 !
 ! 3) Solve X U = cgblock.
 call ZTRSM('Right','Upper','Normal','Normal',npws,nband,cone,cf_ovlp,nband,cgblock,npws)

 ! 4) Solve Y U = gsc. On exit <cgblock|gsc> = 1
 call ZTRSM('Right','Upper','Normal','Normal',npws,nband,cone,cf_ovlp,nband,gsc,npws)

 ABI_FREE(cf_ovlp)

end subroutine cgpaw_cholesky
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgnc_normalize
!! NAME
!!  cgnc_normalize
!!
!! FUNCTION
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband=Number of vectors in icg1
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cgnc_normalize(npws,nband,cg,istwfk,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgnc_normalize'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband,istwfk,me_g0,comm_pw
!arrays
 real(dp),intent(inout) :: cg(2*npws*nband)

!Local variables ------------------------------
!scalars
 integer :: ptr,ierr,band
 character(len=500) :: msg
!arrays
 real(dp) :: norm(nband),alpha(2)

! *************************************************************************

!$OMP PARALLEL DO PRIVATE(ptr) IF (nband > 1)
 do band=1,nband
   ptr = 1 + 2*npws*(band-1)
   norm(band) = cg_dznrm2(npws, cg(ptr))
   norm(band) = norm(band) ** 2
   !norm(band) = cg_real_zdotc(npws, cg(ptr), cg(ptr))
 end do

 if (istwfk>1) then
   norm = two * norm
   if (istwfk==2 .and. me_g0==1) then
!$OMP PARALLEL DO PRIVATE(ptr) IF (nband >1)
     do band=1,nband
       ptr = 1 + 2*npws*(band-1)
       norm(band) = norm(band) - cg(ptr)**2
     end do
   end if
 end if

 if (comm_pw /= xmpi_comm_self) then
   call xmpi_sum(norm,comm_pw,ierr)
 end if

 ierr = 0
 do band=1,nband
   if (norm(band) > zero) then
     norm(band) = SQRT(norm(band))
   else
     ierr = ierr + 1
   end if
 end do

 if (ierr/=0) then
   write(msg,'(a,i0,a)')" Found ",ierr," vectors with norm <= zero!"
   MSG_ERROR(msg)
 end if

!$OMP PARALLEL DO PRIVATE(ptr,alpha) IF (nband > 1)
 do band=1,nband
   ptr = 1 + 2*npws*(band-1)
   alpha = (/one/norm(band), zero/)
   call cg_zscal(npws,alpha,cg(ptr))
 end do

end subroutine cgnc_normalize
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgnc_gsortho
!! NAME
!!  cgnc_gsortho
!!
!! FUNCTION
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband1=Number of vectors in icg1
!!  nband1=Number of vectors in cg2
!!  comm_pw=MPI communicator.
!!
!! SIDE EFFECTS
!!  cg2(2*npws*nband2)
!!  icg1(2*npws*nband1)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!      m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cgnc_gsortho(npws,nband1,icg1,nband2,iocg2,istwfk,normalize,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgnc_gsortho'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband1,nband2,istwfk,me_g0
 integer,optional,intent(in) :: comm_pw
 logical,intent(in) :: normalize
!arrays
 real(dp),intent(in) :: icg1(2*npws*nband1)
 real(dp),intent(inout) :: iocg2(2*npws*nband2)

!Local variables ------------------------------
!scalars
 integer :: ierr,b1,b2
!arrays
 real(dp) :: r_icg1(nband1),r_iocg2(nband2)
 real(dp),allocatable :: proj(:,:,:)

! *************************************************************************

 ABI_MALLOC(proj,(2,nband1,nband2))
 !proj = zero

 ! 1) Calculate <cg1|cg2>
 call cg_zgemm("C","N",npws,nband1,nband2,icg1,iocg2,proj)

 if (istwfk>1) then
   ! nspinor is always 1 in this case.
   ! Account for the missing G and set the imaginary part to zero since wavefunctions are real.
   proj(1,:,:) = two * proj(1,:,:)
   proj(2,:,:) = zero
   !
   if (istwfk==2 .and. me_g0==1) then
     ! Extract the real part at G=0 and subtract its contribution.
     call dcopy(nband1,icg1, 2*npws,r_icg1, 1)
     call dcopy(nband2,iocg2,2*npws,r_iocg2,1)
     do b2=1,nband2
       do b1=1,nband1
         proj(1,b1,b2) = proj(1,b1,b2) - r_icg1(b1) * r_iocg2(b2)
       end do
     end do
   end if
   !
 end if
 !
 ! This is for the MPI version
 if (comm_pw /= xmpi_comm_self) then
   call xmpi_sum(proj,comm_pw,ierr)
 end if

 ! 2) cg2 = cg2 - <cg1|cg2> cg1
 call cg_zgemm("N","N",npws,nband1,nband2,icg1,proj,iocg2,alpha=-cg_cone,beta=cg_cone)

 ABI_FREE(proj)

 ! 3) Normalize iocg2 if required.
 if (normalize) then
   call cgnc_normalize(npws,nband2,iocg2,istwfk,me_g0,comm_pw)
 end if

end subroutine cgnc_gsortho
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgnc_gramschmidt
!! NAME
!!  cgnc_grortho
!!
!! FUNCTION
!!  Gram-Schmidt orthonormalization of the vectors stored in cgblock
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband=Number of band in cgblock
!!  istwfk=Storage mode for the wavefunctions. 1 for standard full mode
!!  me_g0=1 if this node has G=0.
!!  comm_pw=MPI communicator for the planewave group. Set to xmpi_comm_self for sequential mode.
!!
!! SIDE EFFECTS
!!  cgblock(2*npws*nband)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!
!! SOURCE

subroutine cgnc_gramschmidt(npws,nband,cgblock,istwfk,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgnc_gramschmidt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband,istwfk
 integer,intent(in) :: comm_pw,me_g0
!arrays
 real(dp),intent(inout) :: cgblock(2*npws*nband)

!Local variables ------------------------------
!scalars
 integer :: b1,nb2,opt
 logical :: normalize

! *************************************************************************

 ! Normalize the first vector.
 call cgnc_normalize(npws,1,cgblock(1),istwfk,me_g0,comm_pw)
 if (nband == 1) RETURN

 ! Orthogonaluze b1 wrt to the bands in [1,b1-1].
 normalize = .TRUE.
 do b1=2,nband
    opt = 1 + 2*npws*(b1-1)
    nb2=b1-1
    call cgnc_gsortho(npws,nb2,cgblock(1),1,cgblock(opt),istwfk,normalize,me_g0,comm_pw)
 end do

end subroutine cgnc_gramschmidt
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgpaw_normalize
!! NAME
!!  cgpaw_normalize
!!
!! FUNCTION
!!  Normalize a set of PAW pseudo wavefunctions.
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband=Number of band in cgblock and gsc
!!  istwfk=Storage mode for the wavefunctions. 1 for standard full mode
!!  me_g0=1 if this node has G=0.
!!  comm_pw=MPI communicator for the planewave group. Set to xmpi_comm_self for sequential mode.
!!
!! SIDE EFFECTS
!!  cg(2*npws*nband)
!!    input: Input set of vectors |C>
!!    output: Normalized set such as  <C|S|C> = 1
!!  gsc(2*npws*nband)
!!    input: Input set of vectors S|C>
!!    output: New S|C> compute with the new |C>
!!
!! PARENTS
!!      m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cgpaw_normalize(npws,nband,cg,gsc,istwfk,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgpaw_normalize'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband,istwfk,me_g0,comm_pw
!arrays
 real(dp),intent(inout) :: cg(2*npws*nband),gsc(2*npws*nband)

!Local variables ------------------------------
!scalars
 integer :: ptr,ierr,band
 character(len=500) :: msg
!arrays
 real(dp) :: norm(nband),alpha(2)

! *************************************************************************

!$OMP PARALLEL DO PRIVATE(ptr) IF (nband > 1)
 do band=1,nband
   ptr = 1 + 2*npws*(band-1)
   norm(band) = cg_real_zdotc(npws,gsc(ptr),cg(ptr))
 end do

 if (istwfk>1) then
   norm = two * norm
   if (istwfk==2 .and. me_g0==1) then
!$OMP PARALLEL DO PRIVATE(ptr) IF (nband > 1)
     do band=1,nband
       ptr = 1 + 2*npws*(band-1)
       norm(band) = norm(band) - gsc(ptr) * cg(ptr)
     end do
   end if
 end if

 if (comm_pw /= xmpi_comm_self) then
   call xmpi_sum(norm,comm_pw,ierr)
 end if

 ierr = 0
 do band=1,nband
   if (norm(band) > zero) then
     norm(band) = SQRT(norm(band))
   else
     ierr = ierr + 1
   end if
 end do

 if (ierr/=0) then
   write(msg,'(a,i0,a)')" Found ",ierr," vectors with norm <= zero!"
   MSG_ERROR(msg)
 end if

 ! Scale |C> and S|C>.
!$OMP PARALLEL DO PRIVATE(ptr,alpha) IF (nband > 1)
 do band=1,nband
   ptr = 1 + 2*npws*(band-1)
   alpha = (/one/norm(band), zero/)
   call cg_zscal(npws,alpha,cg(ptr))
   call cg_zscal(npws,alpha,gsc(ptr))
 end do

end subroutine cgpaw_normalize
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgpaw_gsortho
!! NAME
!!  cgpaw_gsortho
!!
!! FUNCTION
!!  This routine uses the Gram-Schmidt method to orthogonalize a set of PAW wavefunctions.
!!  with respect to an input block of states.
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband1=Number of vectors in the input block icg1
!!  icg1(2*npws*nband1)=Input block of vectors.
!!  igsc1(2*npws*nband1)= S|C> for C in icg1.
!!  nband2=Number of vectors to orthogonalize
!!  normalize=True if output wavefunction must be normalized.
!!  istwfk=Storage mode for the wavefunctions. 1 for standard full mode
!!  me_g0=1 if this node has G=0.
!!  comm_pw=MPI communicator for the planewave group. Set to xmpi_comm_self for sequential mode.
!!
!! SIDE EFFECTS
!!  iocg2(2*npws*nband2), iogsc2(2*npws*nband1)
!!    input: set of |C> and S|C> wher |C> is the set of states to orthogonalize
!!    output: Orthonormalized set.
!!
!! PARENTS
!!      m_cgtools
!!
!! CHILDREN
!!
!! SOURCE

subroutine cgpaw_gsortho(npws,nband1,icg1,igsc1,nband2,iocg2,iogsc2,istwfk,normalize,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgpaw_gsortho'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband1,nband2,istwfk,me_g0
 integer,optional,intent(in) :: comm_pw
 logical,intent(in) :: normalize
!arrays
 real(dp),intent(in) :: icg1(2*npws*nband1),igsc1(2*npws*nband1)
 real(dp),intent(inout) :: iocg2(2*npws*nband2),iogsc2(2*npws*nband2)

!Local variables ------------------------------
!scalars
 integer :: ierr,b1,b2
!arrays
 real(dp) :: r_icg1(nband1),r_iocg2(nband2)
 real(dp),allocatable :: proj(:,:,:)

! *************************************************************************

 ABI_MALLOC(proj,(2,nband1,nband2))

 ! 1) Calculate <cg1|cg2>
 call cg_zgemm("C","N",npws,nband1,nband2,igsc1,iocg2,proj)

 if (istwfk>1) then
   ! nspinor is always 1 in this case.
   ! Account for the missing G and set the imaginary part to zero since wavefunctions are real.
   proj(1,:,:) = two * proj(1,:,:)
   proj(2,:,:) = zero
   !
   if (istwfk==2 .and. me_g0==1) then
     ! Extract the real part at G=0 and subtract its contribution.
     call dcopy(nband1,igsc1,2*npws,r_icg1, 1)
     call dcopy(nband2,iocg2,2*npws,r_iocg2,1)
     do b2=1,nband2
       do b1=1,nband1
         proj(1,b1,b2) = proj(1,b1,b2) - r_icg1(b1) * r_iocg2(b2)
       end do
     end do
   end if
   !
 end if
 !
 ! This is for the MPI version
 if (comm_pw /= xmpi_comm_self) then
   call xmpi_sum(proj,comm_pw,ierr)
 end if

 ! 2)
 !   cg2 = cg2 - <Scg1|cg2> cg1
 ! S cg2 = S cg2 - <Scg1|cg2> S cg1
 call cg_zgemm("N","N",npws,nband1,nband2,icg1,proj,iocg2,alpha=-cg_cone,beta=cg_cone)
 call cg_zgemm("N","N",npws,nband1,nband2,igsc1,proj,iogsc2,alpha=-cg_cone,beta=cg_cone)

 ABI_FREE(proj)

 ! 3) Normalize iocg2 and iogsc2 if required.
 if (normalize) then
   call cgpaw_normalize(npws,nband2,iocg2,iogsc2,istwfk,me_g0,comm_pw)
 end if

end subroutine cgpaw_gsortho
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cgpaw_gramschmidt
!! NAME
!!  cgpaw_grortho
!!
!! FUNCTION
!!  Gram-Schmidt orthonormalization of the vectors stored in cg
!!
!! INPUTS
!!  npws=Size of each vector (usually npw*nspinor)
!!  nband=Number of bands in cg
!!  istwfk=Storage mode for the wavefunctions. 1 for standard full mode
!!  me_g0=1 if this node has G=0.
!!  comm_pw=MPI communicator for the planewave group. Set to xmpi_comm_self for sequential mode.
!!
!! SIDE EFFECTS
!!  cg(2*npws*nband), gsc(2*npws*nband)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!
!! SOURCE

subroutine cgpaw_gramschmidt(npws,nband,cg,gsc,istwfk,me_g0,comm_pw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgpaw_gramschmidt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,nband,istwfk,comm_pw,me_g0
!arrays
 real(dp),intent(inout) :: cg(2*npws*nband),gsc(2*npws*nband)

!Local variables ------------------------------
!scalars
 integer :: b1,nb2,opt
 logical :: normalize

! *************************************************************************

 ! Normalize the first vector.
 call cgpaw_normalize(npws,1,cg(1),gsc(1),istwfk,me_g0,comm_pw)
 if (nband == 1) RETURN

 ! Orthogonalize b1 wrt to the bands in [1,b1-1].
 normalize = .TRUE.
 do b1=2,nband
    opt = 1 + 2*npws*(b1-1)
    nb2=b1-1
    call cgpaw_gsortho(npws,nb2,cg(1),gsc(1),1,cg(opt),gsc(opt),istwfk,normalize,me_g0,comm_pw)
 end do

end subroutine cgpaw_gramschmidt
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/projbd
!!
!! NAME
!! projbd
!!
!! FUNCTION
!! Project out vector "direc" onto the bands contained in "cg".
!! if useoverlap==0
!!  New direc=direc-$sum_{j/=i} { <cg_{j}|direc>.|cg_{j}> }$
!! if useoverlap==1 (use of overlap matrix S)
!!  New direc=direc-$sum_{j/=i} { <cg_{j}|S|direc>.|cg_{j}> }$
!! (index i can be set to -1 to sum over all bands)
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction coefficients for ALL bands
!!  iband0=which particular band we are interested in ("i" in the above formula)
!!         Can be set to -1 to sum over all bands...
!!  icg=shift to be given to the location of the data in cg
!!  iscg=shift to be given to the location of the data in cg
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcg=maximum size of second dimension of cg
!!  mscg=maximum size of second dimension of scg
!!  nband=number of bands
!!  npw=number of planewaves
!!  nspinor=number of spinorial components (on current proc)
!!  scg(2,mscg*useoverlap)=<G|S|band> for ALL bands,
!!                        where S is an overlap matrix
!!  scprod_io=0 if scprod array has to be computed; 1 if it is input (already in memory)
!!  tim_projbd=timing code of the calling subroutine(can be set to 0 if not attributed)
!!  useoverlap=describe the overlap of wavefunctions:
!!               0: no overlap (S=Identity_matrix)
!!               1: wavefunctions are overlapping
!! me_g0=1 if this processors treats G=0, 0 otherwise.
!! comm=MPI communicator (used if G vectors are distributed.
!!
!! SIDE EFFECTS
!!  direc(2,npw)= input: vector to be orthogonalised with respect to cg (and S)
!!                output: vector that has been orthogonalized wrt cg (and S)
!!
!!  scprod(2,nband)=scalar_product
!!        if useoverlap==0: scalar_product_i=$<cg_{j}|direc_{i}>$
!!        if useoverlap==1: scalar_product_i=$<cg_{j}|S|direc_{i}>$
!!    if scprod_io=0, scprod is output
!!    if scprod_io=1, scprod is input
!!
!! NOTES
!!  1) MPIWF Might have to be recoded for efficient paralellism
!!
!!  2) The new version employs BLAS2 routine so that the OMP parallelism is delegated to BLAS library.
!!
!!  3) Note for PAW: ref.= PRB 73, 235101 (2006), equations (71) and (72):
!!     in normal use, projbd applies P_c projector
!!     if cg and scg are inverted, projbd applies P_c+ projector
!!
!!  4) cg_zgemv wraps ZGEMM whose implementation is more efficient, especially in the threaded case.
!!
!! PARENTS
!!      cgwf,dfpt_cgwf,dfpt_nstpaw,getdc1,lapackprof
!!
!! CHILDREN
!!
!! SOURCE

subroutine projbd(cg,direc,iband0,icg,iscg,istwf_k,mcg,mscg,nband,&
&                 npw,nspinor,scg,scprod,scprod_io,tim_projbd,useoverlap,me_g0,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'projbd'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband0,icg,iscg,istwf_k,mcg,mscg,nband,npw,nspinor
 integer,intent(in) :: scprod_io,tim_projbd,useoverlap,me_g0,comm
!arrays
 real(dp),intent(in) :: cg(2,mcg),scg(2,mscg*useoverlap)
 real(dp),intent(inout) :: direc(2,npw*nspinor)
 real(dp),intent(inout) :: scprod(2,nband)

!Local variables-------------------------------
!scalars
 integer :: nbandm,npw_sp,ierr
!arrays
 real(dp) :: tsec(2),bkp_scprod(2),bkp_dirg0(2)

! *************************************************************************

 call timab(210+tim_projbd,1,tsec)

 npw_sp=npw*nspinor

 nbandm=nband

 if (istwf_k==1) then

   if (scprod_io==0) then
     if (useoverlap==1) then
       call cg_zgemv("C",npw_sp,nbandm,scg(1,iscg+1),direc,scprod)
     else
       call cg_zgemv("C",npw_sp,nbandm,cg(1,icg+1),  direc,scprod)
     end if
     call xmpi_sum(scprod,comm,ierr)
   end if

   if (iband0>0.and.iband0<=nbandm) then
     bkp_scprod = scprod(:,iband0)
     scprod(:,iband0) = zero
   end if

   call cg_zgemv("N",npw_sp,nbandm,cg(1,icg+1),scprod,direc,alpha=-cg_cone,beta=cg_cone)

   if (iband0>0.and.iband0<=nbandm) scprod(:,iband0) = bkp_scprod ! Restore previous value as scprod is output.

 else if (istwf_k>=2) then
   !
   !  u_{G0/2}(G) = u_{G0/2}(-G-G0)^*; k = G0/2
   !  hence:
   !  sum_G f*(G) g(G) = 2 REAL sum_G^{IZONE} w(G) f*(G)g(G)
   !  where
   !  w(G) = 1 except for k=0 and G=0 where w(G=0) = 1/2.
   !
   if (scprod_io==0) then

     if (useoverlap==1) then

       if (istwf_k==2 .and. me_g0==1) then
         bkp_dirg0 = direc(:,1)
         direc(1,1) = half * direc(1,1)
         direc(2,1) = zero
       end if

       call cg_zgemv("C",npw_sp,nbandm,scg(1,iscg+1),direc,scprod)
       scprod = two * scprod
       scprod(2,:) = zero

       if(istwf_k==2 .and. me_g0==1) direc(:,1) = bkp_dirg0

     else

       if (istwf_k==2 .and. me_g0==1) then
         bkp_dirg0 = direc(:,1)
         direc(1,1) = half * direc(1,1)
         direc(2,1) = zero
       end if

       call cg_zgemv("C",npw_sp,nbandm,cg(1,icg+1),direc,scprod)
       scprod = two * scprod
       scprod(2,:) = zero

       if (istwf_k==2 .and. me_g0==1) direc(:,1) = bkp_dirg0
     end if ! useoverlap

     call xmpi_sum(scprod,comm,ierr)
   end if

   if (iband0>0.and.iband0<=nbandm) then
     bkp_scprod = scprod(:,iband0)
     scprod(:,iband0) = zero
   end if

   call cg_zgemv("N",npw_sp,nbandm,cg(1,icg+1),scprod,direc,alpha=-cg_cone,beta=cg_cone)

   if (iband0>0.and.iband0<=nbandm) scprod(:,iband0) = bkp_scprod ! Restore previous value as scprod is output.

 end if ! Test on istwf_k

 call timab(210+tim_projbd,2,tsec)

end subroutine projbd
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_envlop
!!
!! NAME
!! cg_envlop
!!
!! FUNCTION
!! Multiply random number values in cg by envelope function to lower initial kinetic energy.
!! Envelope  $\left( 1-\left( G/G_{\max }\right) ^2\right) ^{power}$ for |G|<= Gmax.
!! Near G=0, little scaling, and goes to zero flatly near Gmax.Loop over perturbations
!!
!! INPUTS
!! cg(2,npw*nband)=initial random number wavefunctions
!! ecut=kinetic energy cutoff in Ha
!! gmet(3,3)=reciprocal space metric (bohr^-2)
!! icgmod=shift to be given to the location of data in cg
!! kg(3,npw)=reduced coordinates of G vectors in basis sphere
!! kpoint(3)=reduced coordinates of k point
!! mcg=maximum second dimension of cg (at least npw*nband*nspinor)
!! nband=number of bands being considered
!! npw=number of planewaves in basis sphere
!! nspinor=number of spinorial components of the wavefunctions
!!
!! OUTPUT
!!  cg(2,mcg)=revised values (not orthonormalized)
!!
!! PARENTS
!!      wfconv
!!
!! CHILDREN
!!
!! SOURCE


subroutine cg_envlop(cg,ecut,gmet,icgmod,kg,kpoint,mcg,nband,npw,nspinor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_envlop'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icgmod,mcg,nband,npw,nspinor
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),kpoint(3)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
 integer,parameter :: re=1,im=2,power=12
 integer :: i1,i2,i3,ig,igs,ispinor,nn,spad
 real(dp) :: cutoff,gs,kpgsqc
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: cut_pws(:)

! *************************************************************************

!$(k+G)^2$ cutoff from $(1/2)(2 Pi (k+G))^2 = ecut$
 kpgsqc=ecut/(2.0_dp*pi**2)
 cutoff = kpgsqc

 ABI_MALLOC(cut_pws,(npw))

!Run through G vectors in basis
!$OMP PARALLEL DO PRIVATE(gs)
 do ig=1,npw
   i1=kg(1,ig) ; i2=kg(2,ig) ; i3=kg(3,ig)
!(k+G)^2 evaluated using metric and kpoint
   gs = gmet(1,1)*(kpoint(1)+dble(i1))**2+&
&    gmet(2,2)*(kpoint(2)+dble(i2))**2+&
&    gmet(3,3)*(kpoint(3)+dble(i3))**2+&
&    2.0_dp*(gmet(2,1)*(kpoint(2)+dble(i2))*(kpoint(1)+dble(i1))+&
&    gmet(3,2)*(kpoint(3)+dble(i3))*(kpoint(2)+dble(i2))+&
&    gmet(1,3)*(kpoint(1)+dble(i1))*(kpoint(3)+dble(i3)))
   if (gs>cutoff) then
     cut_pws(ig) = zero
   else
     cut_pws(ig) = (1.0_dp-(gs/cutoff))**power
   end if
 end do

!Run through bands (real and imaginary components)
!$OMP PARALLEL DO PRIVATE(spad,igs)
 do nn=1,nband
   spad = (nn-1)*npw*nspinor+icgmod
   do ispinor=1,nspinor
     do ig=1,npw
       igs=ig+(ispinor-1)*npw
       cg(1,igs+spad) = cg(1,igs+spad) * cut_pws(ig)
       cg(2,igs+spad) = cg(2,igs+spad) * cut_pws(ig)
     end do
   end do
 end do

 ABI_FREE(cut_pws)

end subroutine cg_envlop
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_normev
!! NAME
!! cg_normev
!!
!! FUNCTION
!! Normalize a set of nband eigenvectors of complex length npw
!! (real length 2*npw) and set phases to make cg(i,i) real and positive.
!! Near convergence, cg(i,j) approaches delta(i,j).
!!
!! INPUTS
!!  cg(2*npw,nband)=unnormalized eigenvectors
!!  npw=dimension of cg as shown
!!  nband=number of eigenvectors and complex length thereof.
!!
!! OUTPUT
!!  cg(2*npw,nband)=nband normalized eigenvectors
!!
!! PARENTS
!!      subdiago
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_normev(cg,npw,nband)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_normev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nband
!arrays
 real(dp),intent(inout) :: cg(2*npw,nband)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: den,evim,evre,phim,phre,xnorm
 character(len=500) :: message

! *************************************************************************

!Loop over vectors
 do ii=1,nband
   ! find norm
   xnorm=0.0d0
   do jj=1,2*npw
     xnorm=xnorm+cg(jj,ii)**2
   end do

   if((xnorm-one)**2>tol6)then
     write(message,'(6a,i6,a,es16.6,3a)' )ch10,&
&     'normev: ',ch10,&
&     'Starting xnorm should be close to one (tol is tol6).',ch10,&
&     'However, for state number',ii,', xnorm=',xnorm,ch10,&
&     'It might be that your LAPACK library has not been correctly installed.'
     MSG_BUG(message)
   end if

   xnorm=1.0d0/sqrt(xnorm)
!  Set up phase
   phre=cg(2*ii-1,ii)
   phim=cg(2*ii,ii)
   if (phim/=0.0d0) then
     den=1.0d0/sqrt(phre**2+phim**2)
     phre=phre*xnorm*den
     phim=phim*xnorm*den
   else
!    give xnorm the same sign as phre (negate if negative)
     phre=sign(xnorm,phre)
     phim=0.0d0
   end if
!  normalize with phase change
   do jj=1,2*npw,2
     evre=cg(jj,ii)
     evim=cg(jj+1,ii)
     cg(jj,ii)=phre*evre+phim*evim
     cg(jj+1,ii)=phre*evim-phim*evre
   end do
 end do

end subroutine cg_normev
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_precon
!!
!! NAME
!! cg_precon
!!
!! FUNCTION
!! precondition $<G|(H-e_{n,k})|C_{n,k}>$
!!
!! INPUTS
!!  $cg(2,npw)=<G|C_{n,k}>$.
!!  $eval=current band eigenvalue=<C_{n,k}|H|C_{n,k}>$.
!!  istwf_k=option parameter that describes the storage of wfs
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  nspinor=number of spinorial components of the wavefunctions
!!  $vect(2,npw)=<G|H|C_{n,k}>$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!  mg_g0=1 if the node treats G0.
!!  comm=MPI communicator
!!
!! OUTPUT
!!  pcon(npw)=preconditionning matrix
!!  vect(2,npw*nspinor)=<G|(H-eval)|C_{n,k}>*(polynomial ratio)
!!
!! PARENTS
!!      cgwf,dfpt_cgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_precon(cg,eval,istwf_k,kinpw,npw,nspinor,me_g0,optekin,pcon,vect,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_precon'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,npw,nspinor,optekin,me_g0,comm
 real(dp),intent(in) :: eval
!arrays
 real(dp),intent(in) :: cg(2,npw*nspinor),kinpw(npw)
 real(dp),intent(inout) :: pcon(npw),vect(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ierr,ig,igs,ipw1,ispinor
 real(dp) :: ek0,ek0_inv,fac,poly,xx
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Compute mean kinetic energy of band
 if(istwf_k==1)then
   ek0=zero
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw
!    !$OMP PARALLEL DO PRIVATE(ig) REDUCTION(+:ek0) SHARED(cg,igs,kinpw,npw)
     do ig=1+igs,npw+igs
       if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
         ek0=ek0+kinpw(ig-igs)*(cg(1,ig)**2+cg(2,ig)**2)
       end if
     end do
   end do

 else if (istwf_k>=2)then
   if (istwf_k==2 .and. me_g0 == 1)then
     ek0=zero ; ipw1=2
     if(kinpw(1)<huge(0.0_dp)*1.d-11)ek0=0.5_dp*kinpw(1)*cg(1,1)**2
   else
     ek0=zero ; ipw1=1
   end if
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw
!    !$OMP PARALLEL DO PRIVATE(ig) REDUCTION(+:ek0) SHARED(cg,ipw1,kinpw,npw)
     do ig=ipw1+igs,npw+igs
       if(kinpw(ig)<huge(0.0_dp)*1.d-11)then
         ek0=ek0+kinpw(ig)*(cg(1,ig)**2+cg(2,ig)**2)
       end if
     end do
   end do
   ek0=2.0_dp*ek0
 end if

 call timab(48,1,tsec)
 call xmpi_sum(ek0,comm,ierr)
 call timab(48,2,tsec)

 if(ek0<1.0d-10)then
   write(message,'(3a)')&
&   'The mean kinetic energy of a wavefunction vanishes.',ch10,&
&   'It is reset to 0.1Ha.'
   MSG_WARNING(message)
   ek0=0.1_dp
 end if

 if (optekin==1) then
   ek0_inv=2.0_dp/(3._dp*ek0)
 else
   ek0_inv=1.0_dp/ek0
 end if

!Carry out preconditioning
 do ispinor=1,nspinor
   igs=(ispinor-1)*npw
!$OMP PARALLEL DO PRIVATE(fac,ig,poly,xx) SHARED(cg,ek0_inv,eval,kinpw,igs,npw,vect,pcon)
   do ig=1+igs,npw+igs
     if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
       xx=kinpw(ig-igs)*ek0_inv
!      Teter polynomial ratio
       poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
       fac=poly/(poly+16._dp*xx**4)
       if (optekin==1) fac=two*fac
       pcon(ig-igs)=fac
       vect(1,ig)=( vect(1,ig)-eval*cg(1,ig) )*fac
       vect(2,ig)=( vect(2,ig)-eval*cg(2,ig) )*fac
     else
       pcon(ig-igs)=zero
       vect(1,ig)=zero
       vect(2,ig)=zero
     end if
   end do
 end do

end subroutine cg_precon
!!***

!!****f* m_cgtools/cg_precon_block
!!
!! NAME
!! cg_precon_block
!!
!! FUNCTION
!! precondition $<G|(H-e_{n,k})|C_{n,k}>$ for a block of band (band-FFT parallelisation)
!! in the case of real WFs (istwfk/=1)
!!
!! INPUTS
!!  blocksize= size of blocks of bands
!!  $cg(vectsize,blocksize)=<G|C_{n,k}> for a block of bands$.
!!  $eval(blocksize,blocksize)=current block of bands eigenvalues=<C_{n,k}|H|C_{n,k}>$.
!!  $ghc(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  iterationnumber=number of iterative minimizations in LOBPCG
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  $vect(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!  optpcon= 0 the TPA preconditionning matrix does not depend on band
!!           1 the TPA preconditionning matrix (not modified)
!!           2 the TPA preconditionning matrix is independant of iterationnumber
!!  vectsize= size of vectors
!!  mg_g0=1 if this node has Gamma, 0 otherwise.
!!
!! OUTPUT
!!  vect(2,npw)=<g|(h-eval)|c_{n,k}>*(polynomial ratio)
!!
!! SIDE EFFECTS
!!  pcon(npw,blocksize)=preconditionning matrix
!!            input  if optpcon=0,2 and iterationnumber/=1
!!            output if optpcon=0,2 and iterationnumber==1
!!
!! PARENTS
!!      m_lobpcg
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_precon_block(cg,eval,blocksize,iterationnumber,kinpw,&
& npw,nspinor,me_g0,optekin,optpcon,pcon,ghc,vect,vectsize,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_precon_block'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,iterationnumber,npw,nspinor,optekin
 integer,intent(in) :: optpcon,vectsize,me_g0,comm
!arrays
 real(dp),intent(in) :: cg(vectsize,blocksize),eval(blocksize,blocksize)
 real(dp),intent(in) :: ghc(vectsize,blocksize),kinpw(npw)
 real(dp),intent(inout) :: pcon(npw,blocksize),vect(vectsize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: iblocksize,ierr,ig,igs,ipw1,ispinor
 real(dp) :: fac,poly,xx
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: ek0(:),ek0_inv(:)

! *************************************************************************

 call timab(536,1,tsec)

!In this case, the Teter, Allan and Payne preconditioner is approximated:
!the factor xx=Ekin(G) and no more Ekin(G)/Ekin(iband)
 if (optpcon==0) then
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw
     if (me_g0 == 1) then
       do ig=1+igs,1+igs !g=0
         if (iterationnumber==1) then
           if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
             xx=kinpw(ig-igs)
!            teter polynomial ratio
             poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
             fac=poly/(poly+16._dp*xx**4)
             if (optekin==1) fac=two*fac
             pcon(ig-igs,1)=fac
             do iblocksize=1,blocksize
               vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
             end do
           else
             pcon(ig-igs,1)=zero
             vect(ig,:)=0.0_dp
           end if
         else
           do iblocksize=1,blocksize
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
           end do
         end if
       end do
       do ig=2+igs,npw+igs
         if (iterationnumber==1) then
           if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
             xx=kinpw(ig-igs)
!            teter polynomial ratio
             poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
             fac=poly/(poly+16._dp*xx**4)
             if (optekin==1) fac=two*fac
             pcon(ig-igs,1)=fac
             do iblocksize=1,blocksize
               vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
               vect(ig+npw-1,iblocksize)=(ghc(ig+npw-1,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig+npw-1,iblocksize))*pcon(ig-igs,1)
             end do
           else
             pcon(ig-igs,1)=zero
             vect(ig,:)=zero
             vect(ig+npw-1,:)=zero
           end if
         else
           do iblocksize=1,blocksize
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
             vect(ig+npw-1,iblocksize)=(ghc(ig+npw-1,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig+npw-1,iblocksize))*pcon(ig-igs,1)
           end do
         end if
       end do
     else
       do ig=1+igs,npw+igs
         if (iterationnumber==1) then
           if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
             xx=kinpw(ig-igs)
!            teter polynomial ratio
             poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
             fac=poly/(poly+16._dp*xx**4)
             if (optekin==1) fac=two*fac
             pcon(ig-igs,1)=fac
             do iblocksize=1,blocksize
               vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
               vect(ig+npw,iblocksize)=(ghc(ig+npw,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig+npw,iblocksize))*pcon(ig-igs,1)
             end do
           else
             pcon(ig-igs,:)=zero
             vect(ig,:)=zero
             vect(ig+npw,:)=zero
           end if
         else
           do iblocksize=1,blocksize
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
             vect(ig+npw,iblocksize)=(ghc(ig+npw,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig+npw,iblocksize))*pcon(ig-igs,1)
           end do
         end if
       end do
     end if
   end do

 else if (optpcon>0) then
!  Compute mean kinetic energy of all bands
   ABI_ALLOCATE(ek0,(blocksize))
   ABI_ALLOCATE(ek0_inv,(blocksize))
   if (iterationnumber==1.or.optpcon==1) then
     do iblocksize=1,blocksize
       if (me_g0 == 1)then
         ek0(iblocksize)=0.0_dp ; ipw1=2
         if(kinpw(1)<huge(0.0_dp)*1.d-11)ek0(iblocksize)=0.5_dp*kinpw(1)*cg(1,iblocksize)**2
         do ig=ipw1,npw
           if(kinpw(ig)<huge(0.0_dp)*1.d-11)then
             ek0(iblocksize)=ek0(iblocksize)+&
&             kinpw(ig)*(cg(ig,iblocksize)**2+cg(ig+npw-1,iblocksize)**2)
           end if
         end do
       else
         ek0(iblocksize)=0.0_dp ; ipw1=1
         do ig=ipw1,npw
           if(kinpw(ig)<huge(0.0_dp)*1.d-11)then
             ek0(iblocksize)=ek0(iblocksize)+&
&             kinpw(ig)*(cg(ig,iblocksize)**2+cg(ig+npw,iblocksize)**2)
           end if
         end do
       end if
     end do

     call xmpi_sum(ek0,comm,ierr)

     do iblocksize=1,blocksize
       if(ek0(iblocksize)<1.0d-10)then
         write(message, '(4a)' )ch10,&
&         'cg_precon_block: the mean kinetic energy of a wavefunction vanishes.',ch10,&
&         'it is reset to 0.1ha.'
         MSG_WARNING(message)
         ek0(iblocksize)=0.1_dp
       end if
     end do
     if (optekin==1) then
       ek0_inv(:)=2.0_dp/(3._dp*ek0(:))
     else
       ek0_inv(:)=1.0_dp/ek0(:)
     end if
   end if !iterationnumber==1.or.optpcon==1

!  Carry out preconditioning
   do iblocksize=1,blocksize
     do ispinor=1,nspinor
       igs=(ispinor-1)*npw
       if (me_g0 == 1) then
         do ig=1+igs,1+igs !g=0
           if (iterationnumber==1.or.optpcon==1) then
             if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
               xx=kinpw(ig-igs)*ek0_inv(iblocksize)
!              teter polynomial ratio
               poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
               fac=poly/(poly+16._dp*xx**4)
               if (optekin==1) fac=two*fac
               pcon(ig-igs,iblocksize)=fac
               vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig,iblocksize))*fac
             else
               pcon(ig-igs,iblocksize)=zero
               vect(ig,iblocksize)=0.0_dp
             end if
           else
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,iblocksize)
           end if
         end do
         do ig=2+igs,npw+igs
           if (iterationnumber==1.or.optpcon==1) then
             if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
               xx=kinpw(ig-igs)*ek0_inv(iblocksize)
!              teter polynomial ratio
               poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
               fac=poly/(poly+16._dp*xx**4)
               if (optekin==1) fac=two*fac
               pcon(ig-igs,iblocksize)=fac
               vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig,iblocksize))*fac
               vect(ig+npw-1,iblocksize)=(ghc(ig+npw-1,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig+npw-1,iblocksize))*fac
             else
               pcon(ig-igs,iblocksize)=zero
               vect(ig,iblocksize)=zero
               vect(ig+npw-1,iblocksize)=zero
             end if
           else
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,iblocksize)
             vect(ig+npw-1,iblocksize)=(ghc(ig+npw-1,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig+npw-1,iblocksize))*pcon(ig-igs,iblocksize)
           end if
         end do
       else
         do ig=1+igs,npw+igs
           if (iterationnumber==1.or.optpcon==1) then
             if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
               xx=kinpw(ig-igs)*ek0_inv(iblocksize)
!              teter polynomial ratio
               poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
               fac=poly/(poly+16._dp*xx**4)
               if (optekin==1) fac=two*fac
               pcon(ig-igs,iblocksize)=fac
               vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig,iblocksize))*fac
               vect(ig+npw,iblocksize)=(ghc(ig+npw,iblocksize)-&
&               eval(iblocksize,iblocksize)*cg(ig+npw,iblocksize))*fac
             else
               pcon(ig-igs,iblocksize)=zero
               vect(ig,iblocksize)=zero
               vect(ig+npw,iblocksize)=zero
             end if
           else
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,iblocksize)
             vect(ig+npw,iblocksize)=(ghc(ig+npw,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig+npw,iblocksize))*pcon(ig-igs,iblocksize)
           end if
         end do
       end if
     end do
   end do
   ABI_DEALLOCATE(ek0)
   ABI_DEALLOCATE(ek0_inv)
 end if !optpcon

 call timab(536,2,tsec)

end subroutine cg_precon_block
!!***

!!****f* m_cgtools/cz_zprecon_block
!!
!! NAME
!! cg_zprecon_block
!!
!! FUNCTION
!! precondition $<G|(H-e_{n,k})|C_{n,k}>$ for a block of band (band-FFT parallelisation)
!!
!! INPUTS
!!  blocksize= size of blocks of bands
!!  $cg(vectsize,blocksize)=<G|C_{n,k}> for a block of bands$.
!!  $eval(blocksize,blocksize)=current block of bands eigenvalues=<C_{n,k}|H|C_{n,k}>$.
!!  $ghc(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  iterationnumber=number of iterative minimizations in LOBPCG
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  $vect(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!  optpcon= 0 the TPA preconditionning matrix does not depend on band
!!           1 the TPA preconditionning matrix (not modified)
!!           2 the TPA preconditionning matrix is independant of iterationnumber
!!  vectsize= size of vectors
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  vect(2,npw)=<g|(h-eval)|c_{n,k}>*(polynomial ratio)
!!
!! SIDE EFFECTS
!!  pcon(npw,blocksize)=preconditionning matrix
!!            input  if optpcon=0,2 and iterationnumber/=1
!!            output if optpcon=0,2 and iterationnumber==1
!!
!! PARENTS
!!      m_lobpcg
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zprecon_block(cg,eval,blocksize,iterationnumber,kinpw,&
&  npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_zprecon_block'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,iterationnumber,npw,nspinor,optekin
 integer,intent(in) :: optpcon,vectsize,comm
!arrays
 real(dp),intent(in) :: kinpw(npw)
 real(dp),intent(inout) :: pcon(npw,blocksize)
 complex(dpc),intent(in) :: cg(vectsize,blocksize),eval(blocksize,blocksize)
 complex(dpc),intent(in) :: ghc(vectsize,blocksize)
 complex(dpc),intent(inout) :: vect(vectsize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: iblocksize,ierr,ig,igs,ispinor
 real(dp) :: fac,poly,xx
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: ek0(:),ek0_inv(:)

! *************************************************************************

 call timab(536,1,tsec)

!In this case, the Teter, Allan and Payne preconditioner is approximated:
!the factor xx=Ekin(G) and no more Ekin(G)/Ekin(iband)
 if (optpcon==0) then
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw
     do ig=1+igs,npw+igs
       if (iterationnumber==1) then
         if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
           xx=kinpw(ig-igs)
!          teter polynomial ratio
           poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
           fac=poly/(poly+16._dp*xx**4)
           if (optekin==1) fac=two*fac
           pcon(ig-igs,1)=fac
           do iblocksize=1,blocksize
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
           end do
         else
           pcon(ig-igs,1)=zero
           vect(ig,:)=dcmplx(0.0_dp,0.0_dp)
         end if
       else
         do iblocksize=1,blocksize
           vect(ig,iblocksize)=(ghc(ig,iblocksize)-eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,1)
         end do
       end if
     end do
   end do

 else if (optpcon>0) then
!  Compute mean kinetic energy of all bands
   ABI_ALLOCATE(ek0,(blocksize))
   ABI_ALLOCATE(ek0_inv,(blocksize))
   if (iterationnumber==1.or.optpcon==1) then
     do iblocksize=1,blocksize
       ek0(iblocksize)=0.0_dp
       do ispinor=1,nspinor
         igs=(ispinor-1)*npw
         do ig=1+igs,npw+igs
           if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
             ek0(iblocksize)=ek0(iblocksize)+kinpw(ig-igs)*&
&             (real(cg(ig,iblocksize))**2+aimag(cg(ig,iblocksize))**2)
           end if
         end do
       end do
     end do

     call xmpi_sum(ek0,comm,ierr)

     do iblocksize=1,blocksize
       if(ek0(iblocksize)<1.0d-10)then
         message = 'the mean kinetic energy of a wavefunction vanishes. it is reset to 0.1ha.'
         MSG_WARNING(message)
         ek0(iblocksize)=0.1_dp
       end if
     end do
     if (optekin==1) then
       ek0_inv(:)=2.0_dp/(3._dp*ek0(:))
     else
       ek0_inv(:)=1.0_dp/ek0(:)
     end if
   end if !iterationnumber==1.or.optpcon==1

!  Carry out preconditioning
   do iblocksize=1,blocksize
     do ispinor=1,nspinor
       igs=(ispinor-1)*npw
       do ig=1+igs,npw+igs
         if (iterationnumber==1.or.optpcon==1) then
           if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
             xx=kinpw(ig-igs)*ek0_inv(iblocksize)
!            teter polynomial ratio
             poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
             fac=poly/(poly+16._dp*xx**4)
             if (optekin==1) fac=two*fac
             pcon(ig-igs,iblocksize)=fac
             vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&             eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,iblocksize)
           else
             pcon(ig-igs,iblocksize)=zero
             vect(ig,iblocksize)=dcmplx(0.0_dp,0.0_dp)
           end if
         else
           vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&           eval(iblocksize,iblocksize)*cg(ig,iblocksize))*pcon(ig-igs,iblocksize)
         end if
       end do
     end do
   end do
   ABI_DEALLOCATE(ek0)
   ABI_DEALLOCATE(ek0_inv)
 end if !optpcon

 call timab(536,2,tsec)

end subroutine cg_zprecon_block
!!***

!!****f* m_cgtools/fxphas_seq
!!
!! NAME
!! fxphas_seq
!!
!! FUNCTION
!! Fix phase of all bands. Keep normalization but maximize real part
!! (minimize imag part). Also fix the sign of real part
!! by setting the first non-zero element to be positive.
!!
!! This version has been stripped of all the mpi_enreg junk by MJV
!!
!! INPUTS
!!  cg(2,mcg)= contains the wavefunction |c> coefficients.
!!  gsc(2,mgsc)= if useoverlap==1, contains the S|c> coefficients,
!!               where S is an overlap matrix.
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  istwfk=input option parameter that describes the storage of wfs
!!    (set to 1 if usual complex vectors)
!!  mcg=size of second dimension of cg
!!  mgsc=size of second dimension of gsc
!!  nband_k=number of bands
!!  npw_k=number of planewaves
!!  useoverlap=describe the overlap of wavefunctions:
!!               0: no overlap (S=Identi0,ty_matrix)
!!               1: wavefunctions are overlapping
!!
!! OUTPUT
!!  cg(2,mcg)=same array with altered phase.
!!  gsc(2,mgsc)= same array with altered phase.
!!
!! NOTES
!! When the sign of the real part was fixed (modif v3.1.3g.6), the
!! test Tv3#5 , dataset 5, behaved differently than previously.
!! This should be cleared up.
!!
!! PARENTS
!!      m_dynmat,rayleigh_ritz
!!
!! CHILDREN
!!
!! SOURCE

subroutine fxphas_seq(cg,gsc,icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fxphas_seq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap
!arrays
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc*useoverlap)

!Local variables-------------------------------
!scalars
 integer :: iband,ii,indx
 real(dp) :: cim,cre,gscim,gscre,quotient,root1,root2,saa,sab,sbb,theta
 real(dp) :: thppi,xx,yy
 character(len=500) :: message
!arrays
 real(dp),allocatable :: cimb(:),creb(:),saab(:),sabb(:),sbbb(:) !,sarr(:,:)

! *************************************************************************

!The general case, where a complex phase indeterminacy is present
 if(istwfk==1)then

   ABI_ALLOCATE(cimb,(nband_k))
   ABI_ALLOCATE(creb,(nband_k))
   ABI_ALLOCATE(saab,(nband_k))
   ABI_ALLOCATE(sabb,(nband_k))
   ABI_ALLOCATE(sbbb,(nband_k))
   cimb(:)=zero ; creb(:)=zero

!  Loop over bands
!  TODO: MG store saa arrays in sarr(3,nband_k) to reduce false sharing.
   do iband=1,nband_k
     indx=icg+(iband-1)*npw_k

!    Compute several sums over Re, Im parts of c
     saa=0.0_dp ; sbb=0.0_dp ; sab=0.0_dp
     do ii=1+indx,npw_k+indx
       saa=saa+cg(1,ii)*cg(1,ii)
       sbb=sbb+cg(2,ii)*cg(2,ii)
       sab=sab+cg(1,ii)*cg(2,ii)
     end do
     saab(iband)=saa
     sbbb(iband)=sbb
     sabb(iband)=sab
   end do ! iband


   do iband=1,nband_k

     indx=icg+(iband-1)*npw_k

     saa=saab(iband)
     sbb=sbbb(iband)
     sab=sabb(iband)

!    Get phase angle theta
     if (sbb+saa>tol8)then
       if(abs(sbb-saa)>tol8*(sbb+saa) .or. 2*abs(sab)>tol8*(sbb+saa))then
         if (abs(sbb-saa)>tol8*abs(sab)) then
           quotient=sab/(sbb-saa)
           theta=0.5_dp*atan(2.0_dp*quotient)
         else
!          Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
           theta=0.25_dp*(pi-(sbb-saa)/sab)
         end if
!        Check roots to get theta for max Re part
         root1=cos(theta)**2*saa+sin(theta)**2*sbb-2.0_dp*cos(theta)*sin(theta)*sab
         thppi=theta+0.5_dp*pi
         root2=cos(thppi)**2*saa+sin(thppi)**2*sbb-2.0_dp*cos(thppi)*sin(thppi)*sab
         if (root2>root1) theta=thppi
       else
!        The real part vector and the imaginary part vector are orthogonal, and of same norm. Strong indeterminacy.
!        Will determine the first non-zero coefficient, and fix its phase
         do ii=1+indx,npw_k+indx
           cre=cg(1,ii)
           cim=cg(2,ii)
           if(cre**2+cim**2>tol8**2*(saa+sbb))then
             if(cre**2>tol8**2**cim**2)then
               theta=atan(cim/cre)
             else
!              Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
               theta=pi/2-cre/cim
             end if
             exit
           end if
         end do
       end if
     else
       write(message,'(a,i0,5a)')&
&       'The eigenvector with band ',iband,' has zero norm.',ch10,&
&       'This usually happens when the number of bands (nband) is comparable to the number of planewaves (mpw)',ch10,&
&       'Action: Check the parameters of the calculation. If nband ~ mpw, then decrease nband or, alternatively, increase ecut'
       MSG_ERROR(message)
     end if

     xx=cos(theta)
     yy=sin(theta)

!    Here, set the first non-zero element to be positive
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       cim=cg(2,ii)
       cre=xx*cre-yy*cim
       if(abs(cre)>tol8)exit
     end do
     if(cre<zero)then
       xx=-xx ; yy=-yy
     end if

     creb(iband)=xx
     cimb(iband)=yy

   end do

   do iband=1,nband_k

     indx=icg+(iband-1)*npw_k

     xx=creb(iband)
     yy=cimb(iband)
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       cim=cg(2,ii)
       cg(1,ii)=xx*cre-yy*cim
       cg(2,ii)=xx*cim+yy*cre
     end do

!    Alter phase of array S|cg>
     if (useoverlap==1) then
       indx=igsc+(iband-1)*npw_k
       do ii=1+indx,npw_k+indx
         gscre=gsc(1,ii)
         gscim=gsc(2,ii)
         gsc(1,ii)=xx*gscre-yy*gscim
         gsc(2,ii)=xx*gscim+yy*gscre
       end do
     end if

   end do ! iband

   ABI_DEALLOCATE(cimb)
   ABI_DEALLOCATE(creb)
   ABI_DEALLOCATE(saab)
   ABI_DEALLOCATE(sabb)
   ABI_DEALLOCATE(sbbb)

!  ====================================================================

!  Storages that take into account the time-reversal symmetry : the freedom is only a sign freedom
 else  ! if istwfk/=1

   ABI_ALLOCATE(creb,(nband_k))
   creb(:)=zero
!  Loop over bands
   do iband=1,nband_k

     indx=icg+(iband-1)*npw_k

!    Here, set the first non-zero real element to be positive
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       if(abs(cre)>tol8)exit
     end do
     creb(iband)=cre

   end do ! iband

   do iband=1,nband_k

     cre=creb(iband)
     if(cre<zero)then
       indx=icg+(iband-1)*npw_k
       do ii=1+indx,npw_k+indx
         cg(1,ii)=-cg(1,ii)
         cg(2,ii)=-cg(2,ii)
       end do
       if(useoverlap==1)then
         indx=igsc+(iband-1)*npw_k
         do ii=1+indx,npw_k+indx
           gsc(1,ii)=-gsc(1,ii)
           gsc(2,ii)=-gsc(2,ii)
         end do
       end if
     end if

   end do ! iband

   ABI_DEALLOCATE(creb)

 end if ! istwfk

end subroutine fxphas_seq
!!***

!!****f* m_cgtools/overlap_g
!! NAME
!! overlap_g
!!
!! FUNCTION
!! Compute the scalar product between WF at two different k-points
!! < u_{n,k1} | u_{n,k2}>
!!
!! INPUTS
!! mpw = maximum dimensioned size of npw
!! npw_k1 = number of plane waves at k1
!! npw_k2 = number of plane waves at k2
!! nspinor = 1 for scalar, 2 for spinor wavefunctions
!! pwind_k = array required to compute the scalar product (see initberry.f)
!! vect1 = wavefunction at k1: | u_{n,k1} >
!! vect2 = wavefunction at k1: | u_{n,k2} >
!!
!! OUTPUT
!! doti = imaginary part of the scalarproduct
!! dotr = real part of the scalarproduct
!!
!! NOTES
!! In case a G-vector of the basis sphere of plane waves at k1
!! does not belong to the basis sphere of plane waves at k2,
!! pwind = 0. Therefore, the dimensions of vect1 &
!! vect2 are (1:2,0:mpw) and the element (1:2,0) MUST be set to zero.
!!
!! The current implementation if not compatible with TR-symmetry (i.e. istwfk/=1) !
!!
!! PARENTS
!!      dfptff_bec,dfptff_die,dfptff_ebp,dfptff_edie,dfptff_gbefd
!!      dfptff_gradberry,qmatrix,smatrix
!!
!! CHILDREN
!!
!! SOURCE

subroutine overlap_g(doti,dotr,mpw,npw_k1,npw_k2,nspinor,pwind_k,vect1,vect2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'overlap_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpw,npw_k1,npw_k2,nspinor
 real(dp),intent(out) :: doti,dotr
!arrays
 integer,intent(in) :: pwind_k(mpw)
 real(dp),intent(in) :: vect1(1:2,0:mpw*nspinor),vect2(1:2,0:mpw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ipw,ispinor,jpw,spnshft1,spnshft2
 character(len=500) :: message

! *************************************************************************

!Check if vect1(:,0) = 0 and vect2(:,0) = 0
 if ((abs(vect1(1,0)) > tol12).or.(abs(vect1(2,0)) > tol12).or. &
& (abs(vect2(1,0)) > tol12).or.(abs(vect2(2,0)) > tol12)) then
   message = ' vect1(:,0) and/or vect2(:,0) are not equal to zero'
   MSG_BUG(message)
 end if

!Compute the scalar product
 dotr = zero; doti = zero
 do ispinor = 1, nspinor
   spnshft1 = (ispinor-1)*npw_k1
   spnshft2 = (ispinor-1)*npw_k2
!$OMP PARALLEL DO PRIVATE(jpw) REDUCTION(+:doti,dotr)
   do ipw = 1, npw_k1
     jpw = pwind_k(ipw)
     dotr = dotr + vect1(1,spnshft1+ipw)*vect2(1,spnshft2+jpw) + vect1(2,spnshft1+ipw)*vect2(2,spnshft2+jpw)
     doti = doti + vect1(1,spnshft1+ipw)*vect2(2,spnshft2+jpw) - vect1(2,spnshft1+ipw)*vect2(1,spnshft2+jpw)
   end do
 end do

end subroutine overlap_g
!!***

END MODULE m_cgtools
!!***

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_cplxtools
!! NAME
!!  m_cplxtools
!!
!! FUNCTION
!! This module defines helper functions to operate on complex arrays (mainly used in the GW code)
!!
!! COPYRIGHT
!! Copyright (C) 1992-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! 1) The convention about names of interfaced routine is: cplx_<name>,
!!    where <name> is equal to the name of the standard BLAS routine
!!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_cplxtools

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use m_fstrings,  only : toupper

 implicit none

 private

 ! Helper functions.
 public :: cplx_fromreal
 public :: cplx_filter

 ! Blas1
 public :: cplx_real_zdotc
 public :: cplx_zaxpby

 ! Blas2
 public :: cplx_zgemv

 !Blas3
 public :: cplx_zgemm

 ! Helper functions for DFT calculations.
 public :: cplx_box2gsph
 public :: cplx_gsph2box
 public :: cplx_setaug_zero
 public :: cplx_setaug_zero_dpc
 public :: cplx_setaug_zero_spc
 public :: cplx_addtorho
!***

 ! Interfaces
 interface cplx_box2gsph
   module procedure cplx_box2gsph_spc
   module procedure cplx_box2gsph_dpc
 end interface cplx_box2gsph

 interface cplx_gsph2box
   module procedure cplx_gsph2box_spc
   module procedure cplx_gsph2box_dpc
 end interface cplx_gsph2box

 interface cplx_setaug_zero
   module procedure cplx_setaug_zero_spc
   module procedure cplx_setaug_zero_dpc
 end interface cplx_setaug_zero

 interface cplx_addtorho
   module procedure cplx_addtorho_dpc
 end interface cplx_addtorho

 !integer,parameter,private :: MIN_SIZE = 5000

 complex(spc),private,parameter :: czero_spc =  (0._sp,0._sp)
 complex(spc),private,parameter :: cone_spc  =  (1._sp,0._sp)
 !complex(spc) ,parameter :: j_spc=(0._sp,1.0_sp)

 complex(dpc),private,parameter :: czero_dpc =  (0._dp,0._dp)
 complex(dpc),private,parameter :: cone_dpc  =  (1._dp,0._dp)
 !complex(dpc) ,parameter :: j_dpc=(0._dp,1.0_dp)

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_fromreal
!! NAME
!!  cplx_fromreal
!!
!! FUNCTION
!!  Convert a real array with (real,imag) part to a complex array
!!
!! INPUTS
!!  n = Specifies the number of elements in ocplx
!!  ireal(2*n)=Input real array.
!!
!! OUTPUT
!!  ocplx(n)=Output complex array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_fromreal(n,ireal,ocplx)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 real(dp),intent(in) :: ireal(2,n)
 complex(dpc),intent(out) :: ocplx(n)

!Local variables ------------------------------
!scalars
 integer :: ii

! *************************************************************************

!$OMP PARALLEL DO PRIVATE(ii)
 do ii=1,n
   ocplx(ii) = DCMPLX(ireal(1,ii),ireal(2,ii))
 end do

end subroutine cplx_fromreal
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_filter
!! NAME
!!  cplx_filter
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
!! CHILDREN
!!
!! SOURCE

subroutine cplx_filter(n, x, mask)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 complex(dpc),intent(inout) :: x(n)
 logical,intent(in) :: mask(n)

! *************************************************************************

 where (mask)
   x = czero
 end where

end subroutine cplx_filter
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_real_zdotc
!! NAME
!!  cplx_real_zdotc
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
!! CHILDREN
!!
!! SOURCE

function cplx_real_zdotc(n,x,y) result(res)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 complex(dpc),intent(in) :: x(n)
 complex(dpc),intent(in) :: y(n)
 real(dp) :: res

!Local variables-------------------------------
 real(dp),external :: ddot

! *************************************************************************

 res = ddot(2*n,x,1,y,1)

end function cplx_real_zdotc
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_zaxpby
!! NAME
!!  cplx_zaxpby
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

subroutine cplx_zaxpby(n,a,x,b,y)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: a,b
!arrays
 complex(dpc),intent(in) :: x(n)
 complex(dpc),intent(inout) :: y(n)

! *************************************************************************

#ifdef HAVE_LINALG_AXPBY
 call zaxpby(n, a, x, 1, b, y, 1)
#else
 call zscal(n, b, y, 1)
 call zaxpy(n, a, x, 1, y,1)
#endif

end subroutine cplx_zaxpby
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_zgemv
!! NAME
!!  cplx_zgemv
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
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_zgemv(trans,nrows,ncols,mat,vec,matvec,alpha,beta)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nrows,ncols
 complex(dpc),optional,intent(in) :: alpha,beta
 character(len=1),intent(in) :: trans
!arrays
 complex(dpc),intent(in) :: mat(nrows*ncols)
 complex(dpc),intent(in) :: vec(*)
 complex(dpc),intent(inout) :: matvec(*)

!Local variables-------------------------------
!scalars
 integer :: mm,nn,kk,lda,ldb,ldc
 complex(dpc) :: my_alpha,my_beta

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

 my_alpha = cone_dpc;  if (PRESENT(alpha)) my_alpha = alpha
 my_beta  = czero_dpc; if (PRESENT(beta))  my_beta  = beta

 call ZGEMM(trans,"N",mm,nn,kk,my_alpha,mat,lda,vec,ldb,my_beta,matvec,ldc)

 !call ZGEMV(trans,mm,nn,my_alpha,mat,lda,vec,1,my_beta,matvec,1)

end subroutine cplx_zgemv
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_zgemm
!! NAME
!!  cplx_zgemm
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
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_zgemm(transa,transb,npws,ncola,ncolb,amat,bmat,cmat,alpha,beta)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npws,ncola,ncolb
 complex(dpc),optional,intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
!arrays
 complex(dpc),intent(in) :: amat(npws*ncola)
 complex(dpc),intent(in) :: bmat(npws*ncolb)
 complex(dpc),intent(inout) :: cmat(*)

!Local variables-------------------------------
!scalars
 integer :: mm,nn,kk,lda,ldb,ldc
 complex(dpc) :: my_alpha,my_beta

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

 my_alpha = cone_dpc;  if (PRESENT(alpha)) my_alpha = alpha
 my_beta  = czero_dpc; if (PRESENT(beta))  my_beta  = beta

 call ZGEMM(transa,transb,mm,nn,kk,my_alpha,amat,lda,bmat,ldb,my_beta,cmat,ldc)

end subroutine cplx_zgemm
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_box2gsph_spc
!! NAME
!!  cplx_box2gsph_spc
!!
!! FUNCTION
!!   Transfer data from the FFT box to the G-sphere. Target SPC complex array.
!!
!! INPUTS
!!  nx,ny,nz=physical dimension of the FFT box.
!!  ldx,ldy,ldz=Logical dimensions of the arrays.
!!  ndat=number of data in iarrbox
!!  npw_k=Number of planewaves in the G-sphere.
!!  kg_k(3,npw_k)=Reduced coordinates of the G-vectoes.
!!  iarrbox(ldx*ldy*ldz*ndat)=Complex Input arrays on the FFT box.
!!  [rscal] = Scaling factor
!!
!! OUTPUT
!!  oarrsph(npw_k*ndat)=Complex Data defined on the G-sphere.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_box2gsph_spc(nx,ny,nz,ldx,ldy,ldz,ndat,npw_k,kg_k,iarrbox,oarrsph,rscal)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat
 real(sp),optional,intent(in) :: rscal
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 complex(spc),intent(in) :: iarrbox(ldx*ldy*ldz*ndat)
 complex(spc),intent(out) :: oarrsph(npw_k*ndat)

!Local variables-------------------------------
!scalars
 integer :: ig,ix,iy,iz,dat,pad_sph,pad_box,ifft,ldxyz

! *************************************************************************

 ldxyz = ldx*ldy*ldz
 if (.not. PRESENT(rscal)) then
   !
   if (ndat==1) then
!$OMP PARALLEL DO PRIVATE(ix,iy,iz,ifft)
     do ig=1,npw_k
       ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       oarrsph(ig) = iarrbox(ifft)
     end do
   else
!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,iy,iz,ifft)
     do dat=1,ndat
       pad_sph = (dat-1)*npw_k
       pad_box = (dat-1)*ldxyz
       do ig=1,npw_k
         ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
         iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
         iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
         ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
         oarrsph(ig+pad_sph) = iarrbox(ifft+pad_box)
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
       oarrsph(ig) = iarrbox(ifft) * rscal
     end do
   else
!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,iy,iz,ifft)
     do dat=1,ndat
       pad_sph = (dat-1)*npw_k
       pad_box = (dat-1)*ldxyz
       do ig=1,npw_k
         ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
         iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
         iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
         ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
         oarrsph(ig+pad_sph) = iarrbox(ifft+pad_box) * rscal
       end do
     end do
   end if
 end if

end subroutine cplx_box2gsph_spc
!!***

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_box2gsph_dpc
!! NAME
!!  cplx_box2gsph_dpc
!!
!! FUNCTION
!!   Transfer data from the FFT box to the G-sphere. Target DPC complex array.
!!
!! INPUTS
!!  nx,ny,nz=physical dimension of the FFT box.
!!  ldx,ldy,ldz=Logical dimensions of the arrays.
!!  ndat=number of data in iarrbox
!!  npw_k=Number of planewaves in the G-sphere.
!!  kg_k(3,npw_k)=Reduced coordinates of the G-vectoes.
!!  iarrbox(ldx*ldy*ldz*ndat)=Complex Input arrays on the FFT box.
!!  [rscal] = Scaling factor
!!
!! OUTPUT
!!  oarrsph(npw_k*ndat)=Complex Data defined on the G-sphere.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_box2gsph_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,npw_k,kg_k,iarrbox,oarrsph,rscal)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nx,ny,nz,ldx,ldy,ldz,ndat
 real(dp),optional,intent(in) :: rscal
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 complex(dpc),intent(in) :: iarrbox(ldx*ldy*ldz*ndat)
 complex(dpc),intent(out) :: oarrsph(npw_k*ndat)

!Local variables-------------------------------
!scalars
 integer :: ig,ix,iy,iz,dat,pad_sph,pad_box,ifft,ldxyz

! *************************************************************************

 ldxyz = ldx*ldy*ldz
 if (.not. PRESENT(rscal)) then
   !
   if (ndat==1) then
!$OMP PARALLEL DO PRIVATE(ix,iy,iz,ifft)
     do ig=1,npw_k
       ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       oarrsph(ig) = iarrbox(ifft)
     end do
   else
!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,iy,iz,ifft)
     do dat=1,ndat
       pad_sph = (dat-1)*npw_k
       pad_box = (dat-1)*ldxyz
       do ig=1,npw_k
         ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
         iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
         iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
         ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
         oarrsph(ig+pad_sph) = iarrbox(ifft+pad_box)
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
       oarrsph(ig) = iarrbox(ifft) * rscal
     end do
   else
!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,iy,iz,ifft)
     do dat=1,ndat
       pad_sph = (dat-1)*npw_k
       pad_box = (dat-1)*ldxyz
       do ig=1,npw_k
         ix=kg_k(1,ig); if (ix<0) ix=ix+nx; ix=ix+1
         iy=kg_k(2,ig); if (iy<0) iy=iy+ny; iy=iy+1
         iz=kg_k(3,ig); if (iz<0) iz=iz+nz; iz=iz+1
         ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
         oarrsph(ig+pad_sph) = iarrbox(ifft+pad_box) * rscal
       end do
     end do
   end if
 end if

end subroutine cplx_box2gsph_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_gsph2box_spc
!! NAME
!! cplx_gsph2box_spc
!!
!! FUNCTION
!! Array iarrsph is defined in sphere with npw points. Insert iarrsph inside box
!! of nx*ny*nz points to define array oarrbox for fft box. rest of oarrbox is filled with 0 s.
!! targer: SPC complex arrays
!!
!! INPUTS
!! iarrsph(2,npw*ndat)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to perform.
!! npw=number of G vectors in basis at this k point
!! oarrbox(2,ldx*ldy*ldz*ndat) = fft box
!! nx,ny,nz=physical dimension of the box (oarrbox)
!! ldx,ldy,ldz=memory dimension of oarrbox
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
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

subroutine cplx_gsph2box_spc(nx,ny,nz,ldx,ldy,ldz,ndat,npw,istwf_k,kg_k,iarrsph,oarrbox)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,nx,ny,nz,ldx,ldy,ldz,ndat,npw
!arrays
 integer,intent(in) :: kg_k(3,npw)
 complex(spc),intent(in) :: iarrsph(npw*ndat)
 complex(spc),intent(out) :: oarrbox(ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: me_g0=1
 integer :: ix,ixinv,iy,iyinv,iz,izinv,dat,ipw,npwmin,pad_box,pad_sph,ifft,ifft_inv,ldxyz
 !character(len=500) :: msg
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
     pad_sph = (dat-1)*npw
     pad_box = (dat-1)*ldxyz
     oarrbox(1+pad_box:ldxyz+pad_box) = czero_spc ! zero the sub-array
     do ipw=1,npw
       ix=kg_k(1,ipw); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ipw); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ipw); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
#if defined __INTEL_COMPILER && defined HAVE_OPENMP
       if (ifft==0) then
         MSG_ERROR("prevent ifort+OMP from miscompiling this section on cronos")
       end if
#endif
       oarrbox(ifft+pad_box) = iarrsph(ipw+pad_sph)
     end do
   end do

 else if (istwf_k>=2) then
   !
   npwmin=1
   if(istwf_k==2 .and. me_g0==1) then ! If gamma point, then oarrbox must be completed
     do dat=1,ndat
       pad_sph = (dat-1)*npw
       pad_box = (dat-1)*ldxyz
       oarrbox(1+pad_box) = REAL(iarrsph(1+pad_sph))
     end do
     npwmin=2
   end if

!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,ixinv,iy,iyinv,iz,izinv,ifft)
   do dat=1,ndat
     pad_sph = (dat-1)*npw
     pad_box = (dat-1)*ldxyz
     oarrbox(npwmin+pad_box:ldxyz+pad_box) = czero_spc
     do ipw=npwmin,npw
       ix=kg_k(1,ipw); if(ix<0)ix=ix+nx; ix=ix+1
       iy=kg_k(2,ipw); if(iy<0)iy=iy+ny; iy=iy+1
       iz=kg_k(3,ipw); if(iz<0)iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       ! Construct the coordinates of -k-G
       ixinv=ixinver(ix); iyinv=iyinver(iy); izinv=izinver(iz)
       ifft_inv = ixinv + (iyinv-1)*ldx + (izinv-1)*ldx*ldy
#if defined __INTEL_COMPILER && defined HAVE_OPENMP
       if (ifft==0 .or. ifft_inv==0) then
         MSG_ERROR("prevent ifort+OMP from miscompiling this section on cronos")
       end if
#endif
       oarrbox(ifft    +pad_box) =       iarrsph(ipw+pad_sph)
       oarrbox(ifft_inv+pad_box) = CONJG(iarrsph(ipw+pad_sph))
     end do
   end do
   !
 else
   MSG_ERROR("Wrong istwfk")
 end if

 if (istwf_k>=2) then
   ABI_FREE(ixinver)
   ABI_FREE(iyinver)
   ABI_FREE(izinver)
 end if

end subroutine cplx_gsph2box_spc
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_gsph2box_dpc
!! NAME
!! cplx_gsph2box_dpc
!!
!! FUNCTION
!! Array iarrsph is defined in sphere with npw points. Insert iarrsph inside box
!! of nx*ny*nz points to define array oarrbox for fft box. rest of oarrbox is filled with 0 s.
!! targer: DPC complex arrays
!!
!! INPUTS
!! iarrsph(2,npw*ndat)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to perform.
!! npw=number of G vectors in basis at this k point
!! oarrbox(2,ldx*ldy*ldz*ndat) = fft box
!! nx,ny,nz=physical dimension of the box (oarrbox)
!! ldx,ldy,ldz=memory dimension of oarrbox
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
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

subroutine cplx_gsph2box_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,npw,istwf_k,kg_k,iarrsph,oarrbox)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,nx,ny,nz,ldx,ldy,ldz,ndat,npw
!arrays
 integer,intent(in) :: kg_k(3,npw)
 complex(dpc),intent(in) :: iarrsph(npw*ndat)
 complex(dpc),intent(out) :: oarrbox(ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: me_g0=1
 integer :: ix,ixinv,iy,iyinv,iz,izinv,dat,ipw,npwmin,pad_box,pad_sph,ifft,ifft_inv,ldxyz
 !character(len=500) :: msg
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
     pad_sph = (dat-1)*npw
     pad_box = (dat-1)*ldxyz
     oarrbox(1+pad_box:ldxyz+pad_box) = czero_dpc ! zero the sub-array
     do ipw=1,npw
       ix=kg_k(1,ipw); if (ix<0) ix=ix+nx; ix=ix+1
       iy=kg_k(2,ipw); if (iy<0) iy=iy+ny; iy=iy+1
       iz=kg_k(3,ipw); if (iz<0) iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
#if defined __INTEL_COMPILER && defined HAVE_OPENMP
       if (ifft==0) then
         MSG_ERROR("prevent ifort+OMP from miscompiling this section on cronos")
       end if
#endif
       oarrbox(ifft+pad_box) = iarrsph(ipw+pad_sph)
     end do
   end do

 else if (istwf_k>=2) then
   !
   npwmin=1
   if(istwf_k==2 .and. me_g0==1) then ! If gamma point, then oarrbox must be completed
     do dat=1,ndat
       pad_sph = (dat-1)*npw
       pad_box = (dat-1)*ldxyz
       oarrbox(1+pad_box) = REAL(iarrsph(1+pad_sph))
     end do
     npwmin=2
   end if

!$OMP PARALLEL DO PRIVATE(pad_sph,pad_box,ix,ixinv,iy,iyinv,iz,izinv,ifft)
   do dat=1,ndat
     pad_sph = (dat-1)*npw
     pad_box = (dat-1)*ldxyz
     oarrbox(npwmin+pad_box:ldxyz+pad_box) = czero_dpc
     do ipw=npwmin,npw
       ix=kg_k(1,ipw); if(ix<0)ix=ix+nx; ix=ix+1
       iy=kg_k(2,ipw); if(iy<0)iy=iy+ny; iy=iy+1
       iz=kg_k(3,ipw); if(iz<0)iz=iz+nz; iz=iz+1
       ifft = ix + (iy-1)*ldx + (iz-1)*ldx*ldy
       ! Construct the coordinates of -k-G
       ixinv=ixinver(ix); iyinv=iyinver(iy); izinv=izinver(iz)
       ifft_inv = ixinv + (iyinv-1)*ldx + (izinv-1)*ldx*ldy
#if defined __INTEL_COMPILER && defined HAVE_OPENMP
       if (ifft==0 .or. ifft_inv==0) then
         MSG_ERROR("prevent ifort+OMP from miscompiling this section on cronos")
       end if
#endif
       oarrbox(ifft    +pad_box) =        iarrsph(ipw+pad_sph)
       oarrbox(ifft_inv+pad_box) = DCONJG(iarrsph(ipw+pad_sph))
     end do
   end do
   !
 else
   MSG_ERROR("Wrong istwfk")
 end if

 if (istwf_k>=2) then
   ABI_FREE(ixinver)
   ABI_FREE(iyinver)
   ABI_FREE(izinver)
 end if

end subroutine cplx_gsph2box_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_setaug_zero_spc
!! NAME
!!  cplx_setaug_zero_spc
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
!!  arr(ldx,ldy,ldz*ndat)= all entries in the augmented region are set to zero
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_setaug_zero_spc(nx,ny,nz,ldx,ldy,ldz,ndat,arr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 complex(spc),intent(inout) :: arr(ldx,ldy,ldz*ndat)

!Local variables-------------------------------
 integer :: iy,iz,dat,padat

! *************************************************************************

 if (nx /= ldx) then
   do iz=1,ldz*ndat
     do iy=1,ldy
       arr(nx+1:ldx,iy,iz) = czero_spc
     end do
   end do
 end if

 if (ny /= ldy) then
   do iz=1,ldz*ndat
     arr(:,ny+1:ldy,iz) = czero_spc
   end do
 end if

 if (nz /= ldz) then
   do dat=1,ndat
     padat = ldz*(dat-1)
     do iz=nz+1,ldz
       arr(:,:,iz+padat) = czero_spc
     end do
   end do
 end if

end subroutine cplx_setaug_zero_spc
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_setaug_zero_dpc
!! NAME
!!  cplx_setaug_zero_dpc
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
!!  arr(ldx,ldy,ldz*ndat)= all entries in the augmented region are set to zero
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_setaug_zero_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,arr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
!arrays
 complex(dpc),intent(inout) :: arr(ldx,ldy,ldz*ndat)

!Local variables-------------------------------
 integer :: iy,iz,dat,padat

! *************************************************************************

 if (nx /= ldx) then
   do iz=1,ldz*ndat
     do iy=1,ldy
       arr(nx+1:ldx,iy,iz) = czero_dpc
     end do
   end do
 end if

 if (ny /= ldy) then
   do iz=1,ldz*ndat
     arr(:,ny+1:ldy,iz) = czero_dpc
   end do
 end if

 if (nz /= ldz) then
   do dat=1,ndat
     padat = ldz*(dat-1)
     do iz=nz+1,ldz
       arr(:,:,iz+padat) = czero_dpc
     end do
   end do
 end if

end subroutine cplx_setaug_zero_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_cplxtools/cplx_addtorho_dpc
!! NAME
!!  cplx_addtorho_dpc
!!
!! FUNCTION
!!  Add |ur|**2 to the ground-states density rho.
!!    rho = rho + weight_r * |ur|**2
!!
!! INPUTS
!!  nx,ny,nz=physical dimension of the FFT box.
!!  ldx,ldy,ldz=leading dimensions of the arrays.
!!  ndat=number of contributions to accumulate.
!!  weight_r=weight used for the accumulation of the density in real space
!!  ur(ldx,ldy,ldz*ndat)=wavefunctions in real space
!!
!! SIDE EFFECTS
!!  rho(ldx,ldy,ldz) = contains the input density at input,
!!                     modified in input with the contribution gived by ur.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cplx_addtorho_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,weight_r,ur,rho)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat
 real(dp),intent(in) :: weight_r
!arrays
 complex(dpc),intent(in) :: ur(ldx*ldy*ldz*ndat)
 real(dp),intent(inout) :: rho(ldx*ldy*ldz)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,dat,ifft,ldxyz,pad_box,padz,pady

! *************************************************************************

 ldxyz = ldx*ldy*ldz

 if (ndat==1) then
!$OMP PARALLEL DO PRIVATE(padz, pady, ifft)
   do iz=1,nz
     padz = (iz-1)*ldx*ldy
     do iy=1,ny
       pady = (iy-1)*ldx
       do ix=1,nx
         ifft = ix + pady + padz
         rho(ifft) = rho(ifft) + weight_r * (REAL(ur(ifft))**2 + AIMAG(ur(ifft))**2)
       end do
     end do
   end do

 else
! It would be nice to use $OMP PARALLEL DO REDUCTION(+:rho)
! but it's risky as the private rho is allocated on the stack of the thread.
!$OMP PARALLEL PRIVATE(pad_box, padz, pady, ifft)
   do dat=1,ndat
     pad_box = (dat-1)*ldxyz
!$OMP DO
     do iz=1,nz
       padz = (iz-1)*ldx*ldy
       do iy=1,ny
         pady = (iy-1)*ldx
         do ix=1,nx
           ifft = ix + pady + padz
           rho(ifft) = rho(ifft) + weight_r * (REAL(ur(ifft+pad_box)**2 + AIMAG(ur(ifft+pad_box))**2))
         end do
       end do
     end do
!$OMP END DO NOWAIT
   end do
!$OMP END PARALLEL
 end if

end subroutine cplx_addtorho_dpc
!!***

!----------------------------------------------------------------------

END MODULE m_cplxtools
!!***

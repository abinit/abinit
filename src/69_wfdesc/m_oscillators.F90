!****m* ABINIT/m_oscillators
!! NAME
!!  m_oscillators
!!
!! FUNCTION
!!  This module contains procedures to calculate the oscillator matrix elements used in the GW code.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_oscillators

 use defs_basis
 use m_errors
 use m_abicore
 use m_fft

 use m_gwdefs,    only : czero_gw
 use m_fstrings,  only : toupper, itoa, sjoin
 use m_geometry,  only : spinrot_cmat
 use m_hide_blas, only : xcopy
 use m_gsphere,   only : gsphere_t

 implicit none

 private
!!***

!----------------------------------------------------------------------

 public :: rho_tw_g           ! Calculate rhotwg(G) = <wfn1| exp(-i(q+G).r) |wfn2>
 public :: calc_wfwfg         ! Calculate the Fourier transform of the product u_{bk}^*(r).u_{b"k}(r) at an arbitrary k in the BZ.
 public :: sym_rhotwgq0       ! Symmetrize the oscillator matrix elements in the BZ in the special case of q=0.
!!***

!----------------------------------------------------------------------

CONTAINS

!!****f* m_oscillators/rho_tw_g
!! NAME
!! rho_tw_g
!!
!! FUNCTION
!! Calculate rhotwg(G) = <wfn1| exp(-i(q+G).r) |wfn2>
!!
!! INPUTS
!! dim_rtwg=Define the size of the output array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=1 if the sum of the matrix elements is wanted.
!!    dim_rtwg=2 if <up|up>, <dwn|dwn> matrix elements are required
!! map2sphere= 1 to retrieve Fourier components indexed according to igfftg0.
!!             0 to retrieve Fourier components indexed according to the FFT box.
!!            NOTE: If map2sphere==0 npwvec must be equal to nr
!! use_padfft= Only compatible with map2sphere 1.
!!             1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfftg0(npwvec*map2sphere)=index of G-G_o in the FFT array for each G in the sphere.
!! i1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! i2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau}
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! npwvec=number of plane waves (in the sphere if map2sphere==1, in the FFT box if map2sphere==1)
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! nspinor=number of spinorial components.
!! spinrot1(4),spinrot2(4)=components of the spinor rotation matrix. See getspinrot
!! wfn1(nr*nspinor*ndat),wfn2(nr*nspinor*ndat)=the two wavefunctions (periodic part)
!! [nhat12(2,nr,nspinor**2*ndat)]=Compensation charge in real space to be added to \Psi_1^*\Psi_2 -- Only for PAW.
!!
!! OUTPUT
!! rhotwg(npwvec)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      m_chi0,m_cohsex,m_exc_build,m_fft_prof,m_prep_calc_ucrpa,m_sigc,m_sigx
!!
!! CHILDREN
!!
!! SOURCE

subroutine rho_tw_g(nspinor,npwvec,nr,ndat,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
& wfn1,i1,ktabr1,ktabp1,spinrot1,&
& wfn2,i2,ktabr2,ktabp2,spinrot2,&
& dim_rtwg,rhotwg) !& nhat12)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: i1,i2,npwvec,nr,nspinor,dim_rtwg,map2sphere,use_padfft,ndat
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: gbound(:,:) !gbound(2*mgfft+8,2)
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 real(dp),intent(in) :: spinrot1(4),spinrot2(4)
 complex(gwpc),intent(in) :: wfn1(nr*nspinor*ndat),wfn2(nr*nspinor*ndat)
 complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg*ndat)
! real(dp),optional,intent(in) :: nhat12(2,nr,nspinor**2*ndat)

!Local variables-------------------------------
!scalars
 integer :: ig,igfft,iab,spad1,spad2,spad0,nx,ny,nz,ldx,ldy,ldz,mgfft
 type(fftbox_plan3_t) :: plan
!arrays
 integer :: spinor_pad(2,4)
 complex(gwpc),allocatable :: u12prod(:),cwavef1(:),cwavef2(:),cwork(:)

! *************************************************************************

 SELECT CASE (nspinor)
 CASE (1)
   ! Collinear case.
   call ts_usug_kkp_bz(npwvec,nr,ndat,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
                       wfn1,i1,ktabr1,ktabp1,&
                       wfn2,i2,ktabr2,ktabp2,rhotwg)

 CASE (2)
   ! Spinorial case.
   ABI_CHECK(ndat==1,"ndat != 1 not coded")
   ABI_MALLOC(cwavef1, (nr * nspinor * ndat))
   ABI_MALLOC(cwavef2, (nr * nspinor * ndat))
   ABI_MALLOC(cwork, (nr * nspinor * ndat))

   call rotate_spinor(i1, ktabr1, ktabp1, spinrot1, nr, nspinor, ndat, wfn1, cwork, cwavef1)
   call rotate_spinor(i2, ktabr2, ktabp2, spinrot2, nr, nspinor, ndat, wfn2, cwork, cwavef2)

   ABI_FREE(cwork)
   ABI_MALLOC(u12prod, (nr))

   spinor_pad = reshape([0, 0, nr, nr, 0, nr, nr, 0], [2, 4])
   rhotwg = czero_gw
   do iab=1,2
     spad1 = spinor_pad(1,iab); spad2=spinor_pad(2,iab)

     u12prod = GWPC_CONJG(cwavef1(spad1+1:spad1+nr)) * cwavef2(spad2+1:spad2+nr)
     ! Add compensation charge.
     !if (PRESENT(nhat12)) u12prod = u12prod + CMPLX(nhat12(1,:,iab),nhat12(2,:,iab))

     spad0 = (iab-1)*npwvec
     SELECT CASE (map2sphere)
     CASE (0)
       ! Need results on the full FFT box thus cannot use zero-padded FFT.
       call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
       call fftbox_execute(plan,u12prod)
       if (dim_rtwg == 1) then
         rhotwg(1:npwvec) = rhotwg(1:npwvec) + u12prod
       else
         rhotwg(spad0+1:spad0+npwvec) = u12prod
       end if

     CASE (1)
       ! Need results on the G-sphere. Call zero-padded FFT routines if required.
       if (use_padfft == 1) then
         nx = ngfft(1); ny = ngfft(2); nz = ngfft(3); mgfft = maxval(ngfft(1:3))
         ldx = nx; ldy = ny; ldz = nz
         call fftpad(u12prod, ngfft, nx, ny, nz, ldx, ldy, ldz, ndat, mgfft, -1, gbound)
       else
         call fftbox_plan3_many(plan, ndat, ngfft(1:3), ngfft(1:3), ngfft(7), -1)
         call fftbox_execute(plan,u12prod)
       end if

       ! Have to map FFT to G-sphere.
       if (dim_rtwg == 1) then
         do ig=1,npwvec
           igfft = igfftg0(ig)
           ! G-G0 belong to the FFT mesh.
           if (igfft /= 0) rhotwg(ig) = rhotwg(ig) + u12prod(igfft)
         end do
       else
         do ig=1,npwvec
           igfft = igfftg0(ig)
           ! G-G0 belong to the FFT mesh.
           if (igfft /= 0) rhotwg(ig+spad0) = u12prod(igfft)
         end do
       end if

     CASE DEFAULT
       ABI_BUG("Wrong map2sphere")
     END SELECT
   end do !iab

   ABI_FREE(u12prod)
   ABI_FREE(cwavef1)
   ABI_FREE(cwavef2)

 CASE DEFAULT
   ABI_BUG('Wrong nspinor')
 END SELECT

end subroutine rho_tw_g
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/ts_usug_kkp_bz
!! NAME
!! ts_usug_kkp_bz
!!
!! FUNCTION
!! Calculate usug(G) = <u1|exp(-i(q+G).r)|u2> for ndat pair of wavefunctions
!! TODO: The routine is thread-safe hence it can be called within an OMP parallel region.
!!
!! INPUTS
!! map2sphere= 1 to retrieve Fourier components indexed according to igfftg0.
!!             0 to retrieve Fourier components indexed according to the FFT box.
!!               NOTE: If map2sphere==0 npw must be equal to nr
!! use_padfft= Only compatible with map2sphere 1.
!!             1 if matrix elements are calculated via zero-padded FFT.
!!             0 R-->G Transform in done on the full FFT box.
!! igfftg0(npw*map2sphere)=index of G-G_o in the FFT array for each G in the sphere.
!! time1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! time2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau}
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! npw=number of plane waves (in the sphere if map2sphere==1, in the FFT box if map2sphere==1)
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! u1(nr*ndat),u2(nr*ndat)=the two wavefunctions (periodic part)
!!
!! OUTPUT
!! usug(npw*ndat)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!
!! SOURCE

subroutine ts_usug_kkp_bz(npw,nr,ndat,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
& u1,time1,ktabr1,ktabp1,&
& u2,time2,ktabr2,ktabp2,usug) !& nhat12)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: time1,time2,npw,nr,map2sphere,use_padfft,ndat
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: gbound(:,:) !gbound(2*mgfft+8,2)
 integer,intent(in) :: igfftg0(npw*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 complex(gwpc),intent(in) :: u1(nr*ndat),u2(nr*ndat)
 complex(gwpc),intent(out) :: usug(npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: nx,ny,nz,ldx,ldy,ldz,mgfft
 type(fftbox_plan3_t) :: plan
!arrays
 complex(gwpc),allocatable :: u12prod(:)

! *************************************************************************

 ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for symmetries:
 ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
 !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 !
 ABI_MALLOC(u12prod,(nr*ndat))
 call usur_kkp_bz(nr,ndat,time1,ktabr1,ktabp1,u1,time2,ktabr2,ktabp2,u2,u12prod)

 ! Add compensation charge.
 !if (PRESENT(nhat12)) u12prod = u1prod + CMPLX(nhat12(1,:,1),nhat12(2,:,1))

 SELECT CASE (map2sphere)
 CASE (0)
   ! Need results on the full FFT box thus cannot use zero-padded FFT.
   call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
   call fftbox_execute(plan,u12prod)
   call xcopy(nr*ndat,u12prod,1,usug,1)

 CASE (1)
   ! Need results on the G-sphere. Call zero-padded FFT routines if required.
   if (use_padfft==1) then
     nx = ngfft(1); ny = ngfft(2); nz = ngfft(3); mgfft = MAXVAL(ngfft(1:3))
     ldx=nx; ldy=ny; ldz=nz
     call fftpad(u12prod,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,-1,gbound)
   else
     call fftbox_plan3_many(plan,ndat,ngfft(1:3),ngfft(1:3),ngfft(7),-1)
     call fftbox_execute(plan,u12prod)
   end if

   ! From the FFT to the G-sphere.
   call gw_box2gsph(nr,ndat,npw,igfftg0,u12prod,usug)

 CASE DEFAULT
   ABI_BUG("Wrong map2sphere")
 END SELECT

 ABI_FREE(u12prod)

end subroutine ts_usug_kkp_bz
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/usur_kkp_bz
!! NAME
!! usur_kkp_bz
!!
!! FUNCTION
!! Calculate u1_kbz^*(r) u2_kbz(r) in real space from the symmetric images in the IBZ.
!! Does not support spinor wavefunctions.
!!
!! INPUTS
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! u1(nr*ndat),u2(nr*ndat)=the two wavefunctions iin the IBZ (periodic part)
!! time1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! time2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau}
!!
!! OUTPUT
!!  u12prod(nr*dat) = u1_kbz^*(r) u2_kbz(r) for the ndat pairs.
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!
!! SOURCE

subroutine usur_kkp_bz(nr,ndat,time1,ktabr1,ktabp1,u1,time2,ktabr2,ktabp2,u2,u12prod)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nr,ndat,time1,time2
 complex(dpc),intent(in) :: ktabp1,ktabp2
!arrays
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 complex(gwpc),intent(in) :: u1(nr*ndat),u2(nr*ndat)
 complex(gwpc),intent(out) :: u12prod(nr*ndat)

!Local variables-------------------------------
!scalars
 integer :: ir,dat,padat
 complex(gwpc) :: my_ktabp1,my_ktabp2
!arrays
 complex(gwpc),allocatable :: u1_bz(:),u2_bz(:)

! *************************************************************************

 ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for symmetries:
 ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
 !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 !
 ABI_MALLOC(u1_bz,(nr*ndat))
 ABI_MALLOC(u2_bz,(nr*ndat))

 my_ktabp1 = ktabp1
 my_ktabp2 = ktabp2

 if (ndat==1) then
   do ir=1,nr
     u1_bz(ir) = u1(ktabr1(ir))*my_ktabp1
   end do
   do ir=1,nr
     u2_bz(ir) = u2(ktabr2(ir))*my_ktabp2
   end do
 else
!$OMP PARALLEL PRIVATE(padat)
!$OMP DO
   do dat=1,ndat
     padat = (dat-1)*nr
     do ir=1,nr
       u1_bz(ir+padat) = u1(ktabr1(ir)+padat)*my_ktabp1
     end do
   end do
!$OMP END DO NOWAIT
!$OMP DO
   do dat=1,ndat
     padat = (dat-1)*nr
     do ir=1,nr
       u2_bz(ir+padat) = u2(ktabr2(ir)+padat)*my_ktabp2
     end do
   end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
 end if

 ! Treat time-reversal.
 SELECT CASE (time1)
 CASE (1)
   if (ndat==1) then
     if (time2==1) then
       do ir=1,nr
         u12prod(ir) = GWPC_CONJG(u1_bz(ir)) * u2_bz(ir)
       end do
     else if (time2==2) then
       do ir=1,nr
         u12prod(ir) = GWPC_CONJG(u1_bz(ir)) * GWPC_CONJG(u2_bz(ir))
       end do
     else
       ABI_ERROR("Wrong time2")
     end if
   else
     if (time2==1) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = GWPC_CONJG(u1_bz(ir+padat)) * u2_bz(ir+padat)
         end do
       end do
     else if (time2==2) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = GWPC_CONJG(u1_bz(ir+padat)) * GWPC_CONJG(u2_bz(ir+padat))
         end do
       end do
     else
       ABI_ERROR("Wrong time2")
     end if
   end if

 CASE (2)
   if (ndat==1) then
     if (time2==1) then
       do ir=1,nr
         u12prod(ir) = u1_bz(ir) * u2_bz(ir)
       end do
     else if (time2==2) then
       do ir=1,nr
         u12prod(ir) = u1_bz(ir) * GWPC_CONJG(u2_bz(ir))
       end do
     else
       ABI_ERROR("Wrong time2")
     end if
   else
     if (time2==1) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = u1_bz(ir+padat) * u2_bz(ir+padat)
         end do
       end do
     else if (time2==2) then
!$OMP PARALLEL DO PRIVATE(padat)
       do dat=1,ndat
         padat = (dat-1)*nr
         do ir=1,nr
           u12prod(ir+padat) = u1_bz(ir+padat) * GWPC_CONJG(u2_bz(ir+padat))
         end do
       end do
     else
       ABI_ERROR("Wrong time2")
     end if
   end if
 CASE DEFAULT
   ABI_ERROR("Wrong time1")
 END SELECT

 ABI_FREE(u1_bz)
 ABI_FREE(u2_bz)

end subroutine usur_kkp_bz
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/gw_box2gsph
!! NAME
!! gw_box2gsph
!!
!! FUNCTION
!! Trasnfer data from the FFT box to the G-sphere.
!!
!! INPUTS
!! nr=number of FFT grid points
!! ndat=Number of wavefunctions to transform.
!! npw=number of plane waves in the sphere
!! igfftg0(npw)=index of G-G_o in the FFT array for each G in the sphere.
!! iarrbox(nr*ndat)=Input array on the FFT mesh
!!
!! OUTPUT
!! oarrsph(npw*ndat)=output array on the sphere.
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!
!! SOURCE

subroutine gw_box2gsph(nr,ndat,npw,igfftg0,iarrbox,oarrsph)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nr,ndat,npw
!arrays
 integer,intent(in) :: igfftg0(npw)
 complex(gwpc),intent(in) :: iarrbox(nr*ndat)
 complex(gwpc),intent(out) :: oarrsph(npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: ig,igfft,dat,pgsp,pfft

! *************************************************************************

 if (ndat==1) then
   do ig=1,npw
     igfft=igfftg0(ig)
     if (igfft/=0) then
       ! G-G0 belong to the FFT mesh.
       oarrsph(ig) = iarrbox(igfft)
     else
       ! Set this component to zero.
       oarrsph(ig) = czero_gw
     end if
   end do
 else
!$OMP PARALLEL DO PRIVATE(pgsp,pfft,igfft)
   do dat=1,ndat
     pgsp = (dat-1)*npw
     pfft = (dat-1)*nr
     do ig=1,npw
       igfft=igfftg0(ig)
       if (igfft/=0) then
         ! G-G0 belong to the FFT mesh.
         oarrsph(ig+pgsp) = iarrbox(igfft+pfft)
       else
         ! Set this component to zero.
         oarrsph(ig+pgsp) = czero_gw
       end if
     end do
   end do
 end if

end subroutine gw_box2gsph
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/calc_wfwfg
!! NAME
!! calc_wfwfg
!!
!! FUNCTION
!!  Calculate the Fourier transform of the product u_{bk}^*(r) u_{b"k}(r)
!!  Return values on the FFT box.
!!
!! INPUTS
!! nspinor=number of spinorial components.
!! spinrot(4)=components of the spinor rotation matrix
!!
!! OUTPUT
!!
!! PARENTS
!!      m_chi0,m_cohsex,m_sigc
!!
!! CHILDREN
!!
!! SOURCE

subroutine calc_wfwfg(ktabr_k,ktabi_k,spinrot,nr,nspinor,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nr,nspinor
!arrays
 integer,intent(in) :: ktabr_k(nr),ngfft_gw(18)
 real(dp),intent(in) :: spinrot(4)
 complex(gwpc),intent(in) :: wfr_jb(nr*nspinor),wfr_kb(nr*nspinor)
 complex(gwpc),intent(out) :: wfg2_jk(nr*nspinor)

!Local variables-------------------------------
 integer,parameter :: ndat1=1
 type(fftbox_plan3_t) :: plan
!arrays
 complex(gwpc),allocatable :: wfr2_dpcplx(:),ujb_bz(:),ukb_bz(:)

! *************************************************************************

 ! There is no need to take into account phases arising from non-symmorphic
 ! operations since the wavefunctions are evaluated at the same k-point.
 ABI_MALLOC(wfr2_dpcplx, (nr * nspinor * ndat1))

 if (nspinor == 1) then
   select case (ktabi_k)
   case (1)
     wfr2_dpcplx = GWPC_CONJG(wfr_jb(ktabr_k)) * wfr_kb(ktabr_k)
   case (2)
     ! Conjugate the product if time-reversal is used to reconstruct this k-point
     wfr2_dpcplx = wfr_jb(ktabr_k) * GWPC_CONJG(wfr_kb(ktabr_k))
   case default
     ABI_ERROR(sjoin("Wrong ktabi_k:", itoa(ktabi_k)))
   end select

 else if (nspinor == 2) then
   ABI_MALLOC(ujb_bz, (nr * nspinor * ndat1))
   ABI_MALLOC(ukb_bz, (nr * nspinor * ndat1))
   ! Use wfr2_dpcplx as workspace array
   call rotate_spinor(ktabi_k, ktabr_k, cone, spinrot, nr, nspinor, ndat1, wfr_jb, wfr2_dpcplx, ujb_bz)
   call rotate_spinor(ktabi_k, ktabr_k, cone, spinrot, nr, nspinor, ndat1, wfr_kb, wfr2_dpcplx, ukb_bz)
   wfr2_dpcplx = GWPC_CONJG(ujb_bz) * ukb_bz
   ABI_FREE(ujb_bz)
   ABI_FREE(ukb_bz)

 else
   ABI_ERROR(sjoin("Wrong nspinor:", itoa(nspinor)))
 end if

 ! Transform to Fourier space (result in wfg2_jk)
 call fftbox_plan3_many(plan, nspinor, ngfft_gw(1:3), ngfft_gw(1:3), ngfft_gw(7), -1)
 call fftbox_execute(plan, wfr2_dpcplx, wfg2_jk)
 ABI_FREE(wfr2_dpcplx)

end subroutine calc_wfwfg
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/sym_rhotwgq0
!! NAME
!!  sym_rhotwgq0
!!
!! FUNCTION
!!  Symmetrization of the oscillator matrix elements <k-q,b1|exp(-i(q+G).r)|k,b2> in the special case of q=0.
!!  The matrix elements in the full BZ is obtained from the matrix elements in the IBZ by
!!  rotating the wavefunctions and taking into account time reversal symmetry.
!!  strictly speaking the symmetrization can be performed only for non-degenerate states.
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).
!!  npw=Number of G-vectors
!!  dim_rtwg=Number of spin-spin combinations, 1 for collinear spin, 4 is nspinor==2 (TODO NOT CODED)
!!  itim_k=2 if time reversal is used to reconstruct the k in the BZ, 1 otherwise.
!!  isym_k=The index of the symmetry symrec rotains k_IBZ onto k_BZ.
!!  rhxtwg_in(dim_rtwg*npw)=The input matrix elements in the IBZ.
!!
!! OUTPUT
!!  rhxtwg_sym(dim_rtwg*npw)=The symmetrized matrix elements in the BZ.
!!
!! NOTES
!! Let M_{G}(k,q) =<k-q,b1|exp(-i(q+G).r)|k,b2>
!!  At q ==0, supposing non-degenerate bands, one obtains:
!!
!!  1) M_{ SG}( Sk) = e^{-iSG.t} M_{G}   (k)
!!  2) M_{-SG}(-Sk) = e^{+iSG.t} M_{G}^* (k)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sym_rhotwgq0(itim_k,isym_k,dim_rtwg,npw,rhxtwg_in,Gsph) result(rhxtwg_sym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,dim_rtwg,itim_k,isym_k
 type(gsphere_t),intent(in) :: Gsph
!arrays
 complex(gwpc),intent(in) :: rhxtwg_in(dim_rtwg*npw)
 complex(gwpc) :: rhxtwg_sym(dim_rtwg*npw)

!Local variables ------------------------------
!scalars
 integer :: ig

!************************************************************************

 ABI_CHECK(dim_rtwg == 1, "dim_rtwg/=1 not coded")

 SELECT CASE (isym_k)
 CASE (1)
   ! Fractional translation associated to E is assumed to be (zero,zero,zero).
   SELECT CASE (itim_k)
   CASE (1)
     ! Identity, no time-reversal. No symmetrization is needed.
     rhxtwg_sym(:) = rhxtwg_in(:)
   CASE (2)
     ! Identity + Time-reversal.
     do ig=1,npw
       rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = GWPC_CONJG(rhxtwg_in(ig))
     end do
   CASE DEFAULT
     ABI_ERROR(sjoin("Wrong value of itim_k:", itoa(itim_k)))
   END SELECT

 CASE DEFAULT
   ! Rotate wavefunctions.
   SELECT CASE (itim_k)
   CASE (1)
     ! no time-reversal, only rotation.
     do ig=1,npw
       rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = rhxtwg_in(ig) * Gsph%phmSGt(ig,isym_k)
     end do
   CASE (2)
     ! time-reversal + spatial rotation.
     do ig=1,npw
       rhxtwg_sym( Gsph%rottb(ig,itim_k,isym_k) ) = GWPC_CONJG( rhxtwg_in(ig) * Gsph%phmSGt(ig,isym_k) )
     end do
   CASE DEFAULT
     ABI_ERROR(sjoin("Wrong value of itim_k:", itoa(itim_k)))
   END SELECT
 END SELECT

end function sym_rhotwgq0
!!***

!----------------------------------------------------------------------

!!****f* m_oscillators/rotate_spinor
!! NAME
!!  rotate_spinor
!!
!! FUNCTION
!!  Return a spinor in the full BZ from its symmetrical image in the IBZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_oscillators
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_spinor(itim_kbz, ktabr_kbz, ktabp_kbz, spinrot, nr, nspinor, ndat, ug_ibz, cwork, oug_bz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itim_kbz, nr, nspinor, ndat
 complex(dpc),intent(in) :: ktabp_kbz
!arrays
 integer,intent(in) :: ktabr_kbz(nr)
 real(dp),intent(in) :: spinrot(4)
 complex(gwpc),intent(in) :: ug_ibz(nr*nspinor*ndat)
 complex(gwpc),intent(out) :: cwork(nr*nspinor*ndat), oug_bz(nr*nspinor*ndat)

!Local variables ------------------------------
!scalars
 integer :: ir,ir1,spad0,ispinor
 complex(gwpc) :: u1a,u1b
!arrays
 complex(dpc) :: spinrot_cmat1(2,2)

!************************************************************************

 ABI_CHECK(ndat == 1, "ndat > 1 not coded")
 ABI_CHECK(nspinor == 2, "nspinor should be 1")

 ! Apply Time-reversal if required.
 ! \psi_{-k}^1 =  (\psi_k^2)^*
 ! \psi_{-k}^2 = -(\psi_k^1)^*
 if (itim_kbz == 1) then
   cwork(:) = ug_ibz(:)
 else if (itim_kbz == 2) then
   cwork(1:nr) = GWPC_CONJG(ug_ibz(nr+1:2*nr))
   cwork(nr+1:2*nr) = -GWPC_CONJG(ug_ibz(1:nr))
 else
   ABI_ERROR(sjoin('Wrong itim_kbz:', itoa(itim_kbz)))
 end if

 do ispinor=1,nspinor
   spad0 = (ispinor-1) * nr
   do ir=1,nr
     ir1 = ktabr_kbz(ir)
     oug_bz(ir+spad0) = cwork(ir1+spad0) * ktabp_kbz
   end do
 end do

 ! Rotation in spinor space (same equations as in wfconv)
 spinrot_cmat1 = spinrot_cmat(spinrot)
 cwork = oug_bz
 do ir=1,nr
   u1a = cwork(ir); u1b = cwork(ir+nr)
   oug_bz(ir)    = spinrot_cmat1(1, 1) * u1a + spinrot_cmat1(1, 2) * u1b
   oug_bz(ir+nr) = spinrot_cmat1(2, 1) * u1a + spinrot_cmat1(2, 2) * u1b
 end do

end subroutine rotate_spinor
!!***

END MODULE m_oscillators
!!***

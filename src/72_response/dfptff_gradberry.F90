!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptff_gradberry
!! NAME
!! dfptff_gradberry
!!
!! FUNCTION
!! Calculation of the gradient of Berry-phase term in finite electric field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (XW).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to finite electric field calculation
!! ikpt = the index of the current k point 
!! isppol=1 for unpolarized, 2 for spin-polarized
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) = 
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>  
!!
!! OUTPUT
!! grad_berry(2,mpw1,dtefield%mband_occ) = the gradient of the Berry phase term
!!
!! PARENTS
!!      dfpt_vtorho
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfptff_gradberry(cg,cg1,dtefield,grad_berry,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&
&                     npwarr,npwar1,nspinor,nsppol,qmat,pwindall)


 use defs_basis
 use m_profiling_abi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptff_gradberry'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: ikpt,isppol,mband,mk1mem,mkmem,mpw,mpw1,nkpt,nspinor
 integer,intent(in) :: nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%mband_occ)

!Local variables -------------------------
!scalars
 integer :: iband,icg,icg1,idir,ikpt1
 integer :: ikpt1m,ikptn,ikptnm,ikptnp1,ipw,jband,jpw,kband
 integer :: mpw_tmp,npw_k1,npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,fac,wfi,wfr
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: z1(2),z2(2)
 real(dp),allocatable :: Amat(:,:,:),Bmat(:,:,:),s1mat(:,:,:),vect1(:,:)
 real(dp),allocatable :: vect2(:,:)

! *************************************************************************

 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(s1mat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 ABI_ALLOCATE(Amat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(Bmat,(2,dtefield%mband_occ,dtefield%mband_occ))
 vect1(:,0) = zero ; vect2(:,0) = zero
 s1mat(:,:,:)=zero
 grad_berry(:,:,:) = zero

 do idir=1,3
   fac = dtefield%efield_dot(idir)*dble(nkpt)/&
&   (dble(dtefield%nstr(idir))*four_pi)

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,3,idir)

   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)

     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg1(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg1(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)

         do  jband = 1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfr - fac*qmat(2,iband,jband,ikpt,1,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfi + fac*qmat(2,iband,jband,ikpt,1,idir)*wfr
           
         end do
       end do
     end if

   end do

!  compute <u^(0)_{k_j+n}|u^(1)_{k_j+1,q}> matrix----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg = dtefield%cgindex(ikptn,isppol)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwarr(ikptn)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,7,idir)


   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

     do iband = 1, dtefield%mband_occ

       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,iband,jband) = dotr
       s1mat(2,iband,jband) = doti

     end do    ! iband
   end do    !jband

!  compute <u^(0)_{-k_j+1}|u^(1)_{-k_j+n},q> matrix--------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikptnm= dtefield%ikpt_dk(ikptn,9,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   ikpt1m= dtefield%ikpt_dk(ikpt1,9,idir)

   icg = dtefield%cgindex(ikpt1m,isppol)
   icg1 = dtefield%cgindex(ikptnm,isppol+nsppol)
   npw_k1 = npwarr(ikpt1m)
   npw_k2 = npwar1(ikptnm)
   pwind_tmp(1:npw_k1) = pwindall((ikpt1m-1)*mpw_tmp+1:(ikpt1m-1)*mpw_tmp+npw_k1,7,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

     do iband = 1, dtefield%mband_occ

       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)

       s1mat(1,jband,iband) = s1mat(1,jband,iband) + dotr
       s1mat(2,jband,iband) = s1mat(2,jband,iband) + doti

     end do    ! iband
   end do    !jband

   Amat(:,:,:)=zero

!  calculate Amat
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Amat(1,iband,jband) = Amat(1,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(1,kband,jband,ikpt,1,idir)&
&         - s1mat(2,iband,kband)*qmat(2,kband,jband,ikpt,1,idir)
         Amat(2,iband,jband) = Amat(2,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(2,kband,jband,ikpt,1,idir)&
&         + s1mat(2,iband,kband)*qmat(1,kband,jband,ikpt,1,idir)
       end do
     end do
   end do

   Bmat(:,:,:)=zero

!  calculate Bmat
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Bmat(1,jband,kband) = Bmat(1,jband,kband) + Amat(1,iband,kband)*&
&         qmat(1,jband,iband,ikptn,1,idir)&
&         - Amat(2,iband,kband)*qmat(2,jband,iband,ikptn,1,idir)
         Bmat(2,jband,kband) = Bmat(2,jband,kband) + Amat(1,iband,kband)*&
&         qmat(2,jband,iband,ikptn,1,idir)&
&         + Amat(2,iband,kband)*qmat(1,jband,iband,ikptn,1,idir)
       end do
     end do
   end do

!  calc. the second term of gradient------------------------------

!  preparation

   ikptnp1 = dtefield%ikpt_dk(ikpt,3,idir)
   icg = dtefield%cgindex(ikptnp1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikptnp1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)

   z1(:) = zero
   z2(:) = zero

   do ipw = 1, npw_k1

     jpw = pwind_tmp(ipw)

     if (jpw > 0) then

       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg + (iband - 1)*npw_k2*nspinor + jpw)

         do jband=1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = grad_berry(1,ipw,jband) &
&           - fac*(Bmat(1,iband,jband)*wfr - Bmat(2,iband,jband)*wfi)
           grad_berry(2,ipw,jband) = grad_berry(2,ipw,jband) &
&           - fac*(Bmat(1,iband,jband)*wfi + Bmat(2,iband,jband)*wfr)

         end do
       end do
     end if
   end do

!  Second part of gradient of Berry phase++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   vect1(:,0) = zero ; vect2(:,0) = zero

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,4,idir)

!  write(std_out,*)'dfpt_cgwf:pwind_tmp',pwind_tmp
!  stop

   do ipw = 1, npw_k1

     jpw = pwind_tmp(ipw)

     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg1(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg1(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)

         do  jband = 1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfr + fac*qmat(2,iband,jband,ikpt,2,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfi - fac*qmat(2,iband,jband,ikpt,2,idir)*wfr

         end do
       end do
     end if
   end do




!  compute <u^(0)_{k_j+n}|u^(1)_{k_j-1,q}> matrix----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg = dtefield%cgindex(ikptn,isppol)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwarr(ikptn)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) =pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,8,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

     do iband = 1, dtefield%mband_occ

       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,iband,jband) = dotr
       s1mat(2,iband,jband) = doti

     end do    ! iband
   end do    !jband

!  compute <u^(0)_{-k_j-1}|u^(1)_{-k_j+n,q}> matrix-----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikptnm= dtefield%ikpt_dk(ikptn,9,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   ikpt1m= dtefield%ikpt_dk(ikpt1,9,idir)
   icg = dtefield%cgindex(ikpt1m,isppol)
   icg1 = dtefield%cgindex(ikptnm,isppol+nsppol)
   npw_k1 = npwarr(ikpt1m)
   npw_k2 = npwar1(ikptnm)
   pwind_tmp(1:npw_k1) =pwindall((ikpt1m-1)*mpw_tmp+1:(ikpt1m-1)*mpw_tmp+npw_k1,8,idir)



   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
     do iband = 1, dtefield%mband_occ
       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,jband,iband) = s1mat(1,jband,iband) + dotr
       s1mat(2,jband,iband) = s1mat(2,jband,iband) + doti

     end do    ! iband
   end do    !jband

   Amat(:,:,:)=zero

!  calculate Amat
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Amat(1,iband,jband) = Amat(1,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(1,kband,jband,ikpt,2,idir)&
&         - s1mat(2,iband,kband)*qmat(2,kband,jband,ikpt,2,idir)
         Amat(2,iband,jband) = Amat(2,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(2,kband,jband,ikpt,2,idir)&
&         + s1mat(2,iband,kband)*qmat(1,kband,jband,ikpt,2,idir)
       end do
     end do
   end do

   Bmat(:,:,:)=zero

!  calculate Bmat
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Bmat(1,jband,kband) = Bmat(1,jband,kband) + Amat(1,iband,kband)*&
&         qmat(1,jband,iband,ikptn,2,idir)&
&         - Amat(2,iband,kband)*qmat(2,jband,iband,ikptn,2,idir)
         Bmat(2,jband,kband) = Bmat(2,jband,kband) + Amat(1,iband,kband)*&
&         qmat(2,jband,iband,ikptn,2,idir)&
         + Amat(2,iband,kband)*qmat(1,jband,iband,ikptn,2,idir)
       end do
     end do
   end do

!  calc. the second term of gradient------------------------------

!  preparation

   ikptnp1 = dtefield%ikpt_dk(ikpt,4,idir)
   icg = dtefield%cgindex(ikptnp1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikptnp1)
   pwind_tmp(1:npw_k1) =pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,6,idir)

   z1(:) = zero
   z2(:) = zero
   do ipw = 1, npw_k1

     jpw = pwind_tmp(ipw)

     if (jpw > 0) then

       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg + (iband - 1)*npw_k2*nspinor + jpw)

         do jband=1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = grad_berry(1,ipw,jband) + fac*(Bmat(1,iband,jband)*wfr - Bmat(2,iband,jband)*wfi)
           grad_berry(2,ipw,jband) = grad_berry(2,ipw,jband) + fac*(Bmat(1,iband,jband)*wfi + Bmat(2,iband,jband)*wfr)

         end do
       end do
     end if
   end do

 end do !idir

 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(Amat)
 ABI_DEALLOCATE(Bmat)
 ABI_DEALLOCATE(pwind_tmp)

end subroutine dfptff_gradberry
!!***

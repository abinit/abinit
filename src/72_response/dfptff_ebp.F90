!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptff_ebp
!! NAME
!! dfptff_ebp
!!
!! FUNCTION
!! calculation of the energy from the term \Omega E \cdot P
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (XW).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! ikpt = the index of the current k point
!! isppol = the index of the spin component
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
!!      dfpt_scfcv
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfptff_ebp(cg,cg1,dtefield,eberry,mband,mkmem,&
&               mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat)

 use defs_basis
 use m_profiling_abi
 use m_efield

 use m_cgtools,   only : overlap_g

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptff_ebp'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,mpw1,nkpt,nspinor,nsppol
 real(dp),intent(out) :: eberry
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)

!Local variables ----------------------------------
!scalars
 integer :: iband,icg,icg1,idir
 integer :: ikpt,ikpt1,ikptn,ikptnm
 integer :: jband,kband,mpw_tmp,npw_k1,npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,e0,fac
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: z1(2)
 real(dp),allocatable :: Amat(:,:,:),umat(:,:,:,:),vect1(:,:),vect2(:,:)

! *************************************************************************

!calculate 4 matrices -----------------------------
 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(umat,(2,dtefield%mband_occ,dtefield%mband_occ,4))
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 ABI_ALLOCATE(Amat,(2,dtefield%mband_occ,dtefield%mband_occ))
 vect1(:,0) = zero ; vect2(:,0) = zero
 eberry=zero

 do ikpt=1,nkpt

   do idir=1,3

     fac = dtefield%efield_dot(idir)/&
&     (dble(dtefield%nstr(idir))*four_pi)

!    compute <u^(1)_{k_j,q}|u^(1)_{k_j+1,q}> matrix---------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwar1(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,3,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)

       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,1) = dotr
         umat(2,iband,jband,1) = doti

       end do    ! iband
     end do    !jband

!    compute <u^(0)_{k_j}|u^(1)_{k_j-n+1,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,5,idir)
     icg = dtefield%cgindex(ikpt,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,7,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,2) = dotr
         umat(2,iband,jband,2) = doti

       end do    ! iband
     end do    !jband

!    compute <u^(1)_{k_j-n,q}|u^(0)_{k_j+1}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn  = dtefield%ikpt_dk(ikpt,8,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikptn,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikptn)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,5,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,3) = dotr
         umat(2,iband,jband,3) = doti

       end do    ! iband
     end do    !jband

!    compute <u^(0)_{-k_j-n+1}|u^(1)_{-k_j,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn = dtefield%ikpt_dk(ikpt,5,idir)
     ikptnm = dtefield%ikpt_dk(ikptn,9,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,9,idir)
     icg = dtefield%cgindex(ikptnm,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikptnm)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptnm-1)*mpw_tmp+1:(ikptnm-1)*mpw_tmp+npw_k1,7,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,4) = dotr
         umat(2,iband,jband,4) = doti

       end do    ! iband
     end do    !jband

!    sum over the whole------------------------------------------------------------

     e0=zero
     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 + 4_dp*(umat(1,iband,jband,1)*qmat(2,jband,iband,ikpt,1,idir)&
&         +       umat(2,iband,jband,1)*qmat(1,jband,iband,ikpt,1,idir))

       end do
     end do

     eberry = eberry - e0*fac

     e0=zero

     ikptn=dtefield%ikpt_dk(ikpt,8,idir)

     Amat(:,:,:)=zero

!    calculate Amat
     do iband=1, dtefield%mband_occ
       do jband=1, dtefield%mband_occ
         do kband=1, dtefield%mband_occ
           Amat(1,iband,jband) = Amat(1,iband,jband) + (umat(1,iband,kband,3))*&
&           qmat(1,kband,jband,ikpt,1,idir)&
&           - (umat(2,iband,kband,3))*qmat(2,kband,jband,ikpt,1,idir)
           Amat(2,iband,jband) = Amat(2,iband,jband) + (umat(1,iband,kband,3))*&
&           qmat(2,kband,jband,ikpt,1,idir)&
&           + (umat(2,iband,kband,3))*qmat(1,kband,jband,ikpt,1,idir)
         end do
       end do
     end do

     do iband=1, dtefield%mband_occ
       do jband=1, dtefield%mband_occ
         do kband=1, dtefield%mband_occ

           z1(1) = (umat(1,jband,iband,4)+umat(1,iband,jband,2))*&
&           qmat(1,jband,kband,ikptn,1,idir)&
&           -    (umat(2,jband,iband,4)+umat(2,iband,jband,2))*&
&           qmat(2,jband,kband,ikptn,1,idir)
           z1(2) = (umat(1,jband,iband,4)+umat(1,iband,jband,2))*&
&           qmat(2,jband,kband,ikptn,1,idir)&
&           +    (umat(2,jband,iband,4)+umat(2,iband,jband,2))*&
&           qmat(1,jband,kband,ikptn,1,idir)

           e0 = e0 - 4_dp*(z1(1)*Amat(2,kband,iband)+z1(2)*Amat(1,kband,iband))

         end do
       end do
     end do

     eberry = eberry - e0*fac

   end do !end idir
 end do !end ikpt

 ABI_DEALLOCATE(umat)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(pwind_tmp)
 ABI_DEALLOCATE(Amat)
end subroutine dfptff_ebp
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptff_die
!! NAME
!! dfptff_die
!!
!! FUNCTION
!! calculate electric susceptibility tensor in Eq.(28) in PRB 75, 115116(2007). 
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (XW).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! idirpert = the current coloumn of the dielectric permittivity tensor
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
!! rprimd(3,3)=dimensional primitive translations in real space (bohr) 
!!
!! OUTPUT
!! diet = electric susceptibility tensor 
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

subroutine dfptff_die(cg,cg1,dtefield,d2lo,idirpert,ipert,mband,mkmem,&
&               mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)

 use defs_basis
 use m_profiling_abi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptff_die'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: idirpert,ipert,mband,mkmem,mpert,mpw,mpw1,nkpt,nspinor
 integer,intent(in) :: nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert) !vz_i

!Local variables ----------------------------------
!scalars
 integer :: ialpha,iband,icg,icg1,idir,ikpt,ikpt1,jband,mpw_tmp,npw_k1
 integer :: npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,e0,fac
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: edir(3)
 real(dp),allocatable :: s1mat(:,:,:),vect1(:,:),vect2(:,:)

! *************************************************************************

!calculate s1 matrices -----------------------------
 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(s1mat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 vect1(:,0) = zero ; vect2(:,0) = zero

 edir(:)=zero

 do ikpt=1,nkpt
   do idir=1,3
!    compute <u^(0)_{k_j}|u^(1)_{k_j+1,q}> matrix--- q=0 ----------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
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
         s1mat(1,iband,jband) = dotr
         s1mat(2,iband,jband) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(1)_{k_j,q}|u^(0)_{k_j+1}> matrix-- q=0 -------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikpt)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)
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
         s1mat(1,iband,jband) = s1mat(1,iband,jband) + dotr
         s1mat(2,iband,jband) = s1mat(2,iband,jband) + doti
       end do    ! iband
     end do    !jband

!    sum over the whole------------------------------------------------------------

     e0=zero

     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 + (s1mat(1,iband,jband)*qmat(2,jband,iband,ikpt,1,idir)&
&         +    s1mat(2,iband,jband)*qmat(1,jband,iband,ikpt,1,idir))
         
       end do
     end do

     do ialpha=1,3
       fac = rprimd(ialpha,idir)/&
&       (dble(dtefield%nstr(idir))*pi)
       edir(ialpha)=edir(ialpha)+ e0*fac
     end do

   end do !idir
 end do !ikpt

 d2lo(1,1:3,ipert,idirpert,ipert)=edir(:)

 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(pwind_tmp)
end subroutine dfptff_die
!!***

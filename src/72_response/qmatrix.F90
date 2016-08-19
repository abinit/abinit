!{\src2tex{textfont=tt}}
!!****f* ABINIT/qmatrix
!! NAME
!! qmatrix
!!
!! FUNCTION
!! calculation of the inverse of the overlap matrix 
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (XW).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! ikpt = the index of the current k point
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!!
!! OUTPUT
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      dzgedi,dzgefa,overlap_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,mband,npwarr,nkpt,nspinor,nsppol,pwindall)

 use defs_basis
 use m_profiling_abi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'qmatrix'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,mpw1,nkpt,nspinor,nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwarr(nkpt),pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(out) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)

!Local variables -------------------------
!scalars
 integer :: iband,icg,icg1,idir,ifor,ikpt,ikpt2,info,jband,job
 integer :: npw_k1,npw_k2,pwmax,pwmin
 integer :: isppol
 real(dp) :: doti,dotr
!arrays
 integer,allocatable :: ipvt(:),pwind_k(:)
 real(dp) :: det(2,2)
 real(dp),allocatable :: sinv(:,:,:),smat_k(:,:,:),vect1(:,:),vect2(:,:)
 real(dp),allocatable :: zgwork(:,:)

! *************************************************************************

 ABI_ALLOCATE(ipvt,(dtefield%mband_occ))
 ABI_ALLOCATE(sinv,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(zgwork,(2,dtefield%mband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw))
 ABI_ALLOCATE(vect2,(2,0:mpw))
 ABI_ALLOCATE(smat_k,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(pwind_k,(max(mpw,mpw1)))
 vect1(:,0) = zero ; vect2(:,0) = zero

 job = 11

!**************
!loop over k points
 do isppol = 1, nsppol
   do ikpt = 1, nkpt
     npw_k1 = npwarr(ikpt)
     icg  = dtefield%cgindex(ikpt,1)
  
     do idir = 1, 3
       
       do ifor = 1, 2
  
         ikpt2 = dtefield%ikpt_dk(ikpt,ifor,idir)
         npw_k2 = npwarr(ikpt2)
         icg1 = dtefield%cgindex(ikpt2,1)
         pwind_k(1:npw_k1) = pwindall((ikpt-1)*max(mpw,mpw1)+1:(ikpt-1)*max(mpw,mpw1)+npw_k1,ifor,idir)
  
         do jband = 1, dtefield%mband_occ
           vect2(:,1:npw_k2) = &
&           cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
           if (npw_k2 < mpw) vect2(:,npw_k2+1:mpw) = zero
  
           do iband = 1, dtefield%mband_occ
  
             pwmin = (iband-1)*npw_k1*nspinor
             pwmax = pwmin + npw_k1*nspinor
             vect1(:,1:npw_k1) = &
&             cg(:,icg + 1 + pwmin:icg + pwmax)
             if (npw_k1 < mpw) vect1(:,npw_k1+1:mpw) = zero
             call overlap_g(doti,dotr,mpw,npw_k1,npw_k2,nspinor,pwind_k,&
&             vect1,vect2)
             smat_k(1,iband,jband) = dotr
             smat_k(2,iband,jband) = doti
  
           end do    ! iband
  
         end do    !jband
  
         sinv(:,:,:) = smat_k(:,:,:)
  
         call dzgefa(sinv,dtefield%nband_occ(isppol),dtefield%nband_occ(isppol),ipvt,info)
         call dzgedi(sinv,dtefield%nband_occ(isppol),dtefield%nband_occ(isppol),ipvt,det,zgwork,job)
  
         qmat(:,:,:,ikpt,ifor,idir) = sinv(:,:,:)
  
  
       end do
  
     end do

   end do  !end loop over k 
 end do

 ABI_DEALLOCATE(ipvt)
 ABI_DEALLOCATE(sinv)
 ABI_DEALLOCATE(zgwork)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(smat_k)
 ABI_DEALLOCATE(pwind_k)

end subroutine qmatrix
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix_k_paw
!! NAME
!! smatrix_k_paw
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj_k (pawcprj_type) :: cprj for occupied bands at point k
!!  cprj_kb :: cprj for occupied bands at point k+b
!!  dtefield :: structure referring to all efield and berry's phase variables
!!  kdir :: integer giving direction along which overlap is computed for ket
!!  kfor :: integer indicating whether to compute forward (1) or backward (2)
!!    along kpt string
!!  natom :: number of atoms in cell
!!  typat :: typat(natom) type of each atom
!!
!! OUTPUT
!! smat_k_paw :: array of the on-site PAW parts of the overlaps between Bloch states at points
!!   k and k+b, for the various pairs of bands, that is, the on-site part of 
!!   <u_nk|u_mk+b>
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine assumes that the cprj are not explicitly ordered by 
!! atom type.
!!
!! PARENTS
!!      berryphase_new,cgwf,make_grad_berry
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine smatrix_k_paw(cprj_k,cprj_kb,dtefield,kdir,kfor,mband,natom,smat_k_paw,typat)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_efield
 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smatrix_k_paw'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: kdir,kfor,mband,natom
 type(efield_type),intent(in) :: dtefield
 type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%nspinor*mband)
 type(pawcprj_type),intent(in) :: cprj_kb(natom,dtefield%nspinor*mband)

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: smat_k_paw(2,dtefield%mband_occ,dtefield%mband_occ)

!Local variables---------------------------
!scalars
 integer :: iatom,iband,ibs,ilmn,ispinor,itypat
 integer :: jband,jbs,jlmn,klmn,nspinor
 complex(dpc) :: cpk,cpkb,cterm,paw_onsite

! *************************************************************************

!initialize smat_k_paw
 smat_k_paw(:,:,:) = zero

 nspinor = dtefield%nspinor

 do iatom = 1, natom
   itypat = typat(iatom)

   do ilmn=1,dtefield%lmn_size(itypat)
     do jlmn=1,dtefield%lmn_size(itypat)
       klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
       paw_onsite = cmplx(dtefield%qijb_kk(1,klmn,iatom,kdir),&
&       dtefield%qijb_kk(2,klmn,iatom,kdir))
       if (kfor > 1) paw_onsite = conjg(paw_onsite)
       do iband = 1, dtefield%mband_occ
         do jband = 1, dtefield%mband_occ
           do ispinor = 1, nspinor
             ibs = nspinor*(iband-1) + ispinor
             jbs = nspinor*(jband-1) + ispinor
             cpk=cmplx(cprj_k(iatom,ibs)%cp(1,ilmn),cprj_k(iatom,ibs)%cp(2,ilmn))
             cpkb=cmplx(cprj_kb(iatom,jbs)%cp(1,jlmn),cprj_kb(iatom,jbs)%cp(2,jlmn))
             cterm = conjg(cpk)*paw_onsite*cpkb
             smat_k_paw(1,iband,jband) = smat_k_paw(1,iband,jband)+dreal(cterm)
             smat_k_paw(2,iband,jband) = smat_k_paw(2,iband,jband)+dimag(cterm)
           end do ! end loop over ispinor
         end do ! end loop over jband
       end do ! end loop over iband
     end do ! end loop over ilmn
   end do ! end loop over jlmn

 end do ! end loop over atoms

 end subroutine    smatrix_k_paw
!!***


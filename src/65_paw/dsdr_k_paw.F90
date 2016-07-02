!{\src2tex{textfont=tt}}
!!****f* ABINIT/dsdr_k_paw
!! NAME
!! dsdr_k_paw
!!
!! FUNCTION
!! compute on-site terms for forces and stresses for finite electric fields with PAW
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
!!
!! SIDE EFFECTS
!! dsdr :: array of the on-site PAW parts of the derivatives with respect to atm
!!               positions and/or strains of the overlaps between Bloch states at points
!!               k and k+b, for the various pairs of bands
!!
!! NOTES
!! This routine assumes that the cprj are not explicitly ordered by 
!! atom type.
!!
!! PARENTS
!!      berryphase_new
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine dsdr_k_paw(cprj_k,cprj_kb,dsdr,dtefield,kdir,kfor,mband,natom,ncpgr,typat)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_efield
 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dsdr_k_paw'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: kdir,kfor,mband,natom,ncpgr
 character(len=500) :: message
 type(efield_type),intent(in) :: dtefield
 type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%nspinor*mband)
 type(pawcprj_type),intent(in) :: cprj_kb(natom,dtefield%nspinor*mband)

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(inout) :: dsdr(2,natom,ncpgr,dtefield%nband_occ,dtefield%nband_occ)

!Local variables---------------------------
!scalars
 integer :: iatom,iband,ibs,icpgr,ilmn,ispinor,itypat
 integer :: jband,jbs,jlmn,klmn,nspinor
 complex(dpc) :: cpk,cpkb,dcpk,dcpkb,cterm,paw_onsite

! *************************************************************************

!initialize dsdr
 dsdr(:,:,:,:,:) = zero

! if 3 gradients we are in the ctocprj choice 2 case
! and the 3 gradients are due to the atomic displacements
! if 6 gradients we are in the ctocprj choice 3 case
! and the 6 gradients are due to the strains
! if 9 gradients we are in the ctocprj choice 23 case
! and the first six are due to strain, last three due to displacements
 if (ncpgr /= 3 .and. ncpgr /= 6 .and. ncpgr /= 9) then
   message = ' dsdr_k_paw called with ncpgr /= 3, 6, or 9 (no gradients) '
   MSG_BUG(message)
 end if
 
 nspinor = dtefield%nspinor

 do iatom = 1, natom
   itypat = typat(iatom)

   do ilmn=1,dtefield%lmn_size(itypat)
     do jlmn=1,dtefield%lmn_size(itypat)
       klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
       paw_onsite = cmplx(dtefield%qijb_kk(1,klmn,iatom,kdir),&
&       dtefield%qijb_kk(2,klmn,iatom,kdir))
       if (kfor > 1) paw_onsite = conjg(paw_onsite)
       do iband = 1, dtefield%nband_occ
         do jband = 1, dtefield%nband_occ
           do ispinor = 1, nspinor
             do icpgr = 1, ncpgr
               ibs = nspinor*(iband-1) + ispinor
               jbs = nspinor*(jband-1) + ispinor
               cpk=cmplx(cprj_k(iatom,ibs)%cp(1,ilmn),cprj_k(iatom,ibs)%cp(2,ilmn))
               dcpk=cmplx(cprj_k(iatom,ibs)%dcp(1,icpgr,ilmn),cprj_k(iatom,ibs)%dcp(2,icpgr,ilmn))
               cpkb=cmplx(cprj_kb(iatom,jbs)%cp(1,jlmn),cprj_kb(iatom,jbs)%cp(2,jlmn))
               dcpkb=cmplx(cprj_kb(iatom,jbs)%dcp(1,icpgr,jlmn),cprj_kb(iatom,jbs)%dcp(2,icpgr,jlmn))
               cterm=paw_onsite*(conjg(dcpk)*cpkb+conjg(cpk)*dcpkb)
               dsdr(1,iatom,icpgr,iband,jband) = dsdr(1,iatom,icpgr,iband,jband)+real(cterm)
               dsdr(2,iatom,icpgr,iband,jband) = dsdr(2,iatom,icpgr,iband,jband)+aimag(cterm)
             end do ! end loop over icpgr
           end do ! end loop over ispinor
         end do ! end loop over jband
       end do ! end loop over iband
     end do ! end loop over ilmn
   end do ! end loop over jlmn

 end do ! end loop over atoms

 end subroutine    dsdr_k_paw
!!***


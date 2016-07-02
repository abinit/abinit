!{\src2tex{textfont=tt}}
!!****f* ABINIT/fxgkkphase
!!
!! NAME
!! fxgkkphase
!!
!! FUNCTION
!!   Set phase factors to eliminate gauge variable in gkk matrix elements
!!    (comes from phase factors in wavefunctions)
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = datastructure for elph data (dimensions and eventually data)
!!   gkk_flag = flags for presence of gkk matrix elements
!!   h1_mat_el = irreducible matrix elements to be completed
!!   iqptfull = qpoint number in full zone
!!
!! OUTPUT
!!   h1_mat_el = changed on output
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


subroutine fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)

 use defs_basis
 use defs_elphon
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fxgkkphase'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: gkk_flag(elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%k_phon%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ib1,ib2,ibranch
 real(dp) :: imphase,norm,rephase,tmpim
 real(dp) :: tmpre
!arrays

! *************************************************************************

 do ikpt_phon=1,elph_ds%k_phon%nkpt
   do ib1=1,elph_ds%nFSband
     do ib2=1,elph_ds%nFSband
!      Determine phase for first perturbation
       rephase =  h1_mat_el(1,ib2,ib1,1,ikpt_phon)
       imphase = -h1_mat_el(2,ib2,ib1,1,ikpt_phon)
       norm = sqrt(rephase**2+imphase**2)
       rephase = rephase / norm
       imphase = imphase / norm

!      DEBUG
!      if (ikpt_phon == 1) then
!      write(std_out,*) 'fxgkkphase : rephase,imphase = ', rephase,imphase
!      end if
!      ENDDEBUG

!      apply same phase factor to all perturbations
!      ----------------------------------------------------------
!      Very important ! Otherwise the scalar product with the
!      displacement vector will not be preserved.
!      ----------------------------------------------------------
       do ibranch=1,elph_ds%nbranch
!        if we already have data
         if (gkk_flag(ibranch,1,iqptfull) /= -1) then
           tmpre =    rephase * h1_mat_el(1,ib2,ib1,ibranch,ikpt_phon)&
&           -imphase * h1_mat_el(2,ib2,ib1,ibranch,ikpt_phon)
           tmpim =    rephase * h1_mat_el(2,ib2,ib1,ibranch,ikpt_phon)&
&           +imphase * h1_mat_el(1,ib2,ib1,ibranch,ikpt_phon)
           h1_mat_el(1,ib2,ib1,ibranch,ikpt_phon) = tmpre
           h1_mat_el(2,ib2,ib1,ibranch,ikpt_phon) = tmpim
         end if
!        end if
       end do
     end do
   end do
 end do
!end loop over FS kpt

end subroutine fxgkkphase
!!***

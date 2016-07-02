!{\src2tex{textfont=tt}}
!!****f* ABINIT/test_ftgkk
!! NAME
!! test_ftgkk
!!
!! FUNCTION
!!  Test the fourier transform routine ftgkk for the el-phon matrix elements
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVerstra)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with matrix elements
!!   gprim = reciprocal lattice vectors
!!   natom = number of atoms
!!   nrpt = number of real space points for FT interpolation
!!   rpt = coordinates of real space points for FT interpolation
!!   qpt_full = qpoint coordinates
!!   wghatm = weights for pairs of atoms in FT interpolation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!  MJV 18/5/2008 reverted to old syntax/use for ftgkk, with all ft being done
!!   in a batch. Might come back to 5.5 version with atomic FT in ftgkk, but later.
!!
!! PARENTS
!!
!! CHILDREN
!!      ftgkk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,qpt_full,wghatm)

 use defs_basis
 use defs_elphon
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'test_ftgkk'
 use interfaces_77_ddb, except_this_one => test_ftgkk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt_full(3,elph_ds%nqpt_full)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,iqpt,isppol,qtor,sz1,sz2
!arrays
 real(dp),allocatable :: gkq_disk(:,:,:,:,:),tmp_gkq(:,:,:,:,:)

! *************************************************************************

!for each qpt do FT to recuperate original values

 isppol = 1
 qtor = 0
 sz1=elph_ds%ngkkband*elph_ds%ngkkband
 sz2=elph_ds%nbranch*elph_ds%nbranch
 ABI_ALLOCATE(gkq_disk,(2,sz1,sz2,elph_ds%k_phon%nkpt,elph_ds%nsppol))
 ABI_ALLOCATE(tmp_gkq,(2,sz1,sz2,elph_ds%k_phon%nkpt,elph_ds%nsppol))

 do iqpt=1,elph_ds%nqpt_full
   tmp_gkq(:,:,:,:,:) = zero

   call ftgkk (wghatm,tmp_gkq,elph_ds%gkk_rpt,elph_ds%gkqwrite,&
&   elph_ds%gkk_rptwrite,gprim,1,natom,&
&   elph_ds%k_phon%nkpt,elph_ds%ngkkband,elph_ds%k_phon%nkpt,1,&
&   nrpt,elph_ds%nsppol,qtor,rpt,qpt_full,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

   if (elph_ds%gkqwrite == 0) then
     do ikpt_phon=1,10
       write (93,*) tmp_gkq(:,:,:,ikpt_phon,isppol)-elph_ds%gkk_qpt(:,:,:,ikpt_phon,isppol,iqpt)
     end do
   else
     do ikpt_phon=1, elph_ds%k_phon%nkpt
       read (elph_ds%unitgkq,REC=((iqpt-1)*elph_ds%k_phon%nkpt+ikpt_phon)) gkq_disk(:,:,:,ikpt_phon,:)
     end do
     do ikpt_phon=1,10
       write (93,*) tmp_gkq(:,:,:,ikpt_phon,isppol)-gkq_disk(:,:,:,ikpt_phon,isppol)
     end do
   end if
 end do

end subroutine test_ftgkk
!!***

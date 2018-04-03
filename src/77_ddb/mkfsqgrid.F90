!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkfsqgrid
!!
!! NAME
!! mkfsqgrid
!!
!! FUNCTION
!! This routine sets up the qpoints between the full FS kpt grid points
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kpt_phon = kpoints on the FS
!!  nkpt_phon = number of kpts on the FS
!!
!! OUTPUT
!!  nFSqpt = full number of qpoints between points on the FS
!!  FStoqpt = qpoint index for each pair of kpt on the FS
!!  tmpFSqpt = temp array with coordinates of the
!!    qpoints between points on the FS
!!
!! NOTES
!!   might need the inverse indexing qpt -> (kpt1,kpt2) (kpt3,kpt4) ...
!!
!! PARENTS
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkfsqgrid(kpt_phon,FStoqpt,nkpt_phon,nFSqpt,tmpFSqpt)

 use defs_basis
 use m_profiling_abi
 use m_numeric_tools, only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkfsqgrid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_phon
 integer,intent(out) :: nFSqpt
!arrays
 integer,intent(out) :: FStoqpt(nkpt_phon,nkpt_phon)
 real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
 real(dp),intent(out) :: tmpFSqpt(3,nkpt_phon*nkpt_phon)

!Local variables-------------------------------
!scalars
 integer :: iFSqpt,ikpt1,ikpt2,new
 real(dp) :: shift,ss
!arrays
 real(dp) :: k1(3),kpt(3)

! *************************************************************************

 nFSqpt=0
 tmpFSqpt(:,:)=zero
 do ikpt1=1,nkpt_phon
   do ikpt2=1,nkpt_phon
     k1(:) = kpt_phon(:,ikpt1) - kpt_phon(:,ikpt2)
     call wrap2_pmhalf(k1(1),kpt(1),shift)
     call wrap2_pmhalf(k1(2),kpt(2),shift)
     call wrap2_pmhalf(k1(3),kpt(3),shift)

     new=1
!    is kpt among the FS qpts found already?
     do iFSqpt=1,nFSqpt
       ss=(kpt(1)-tmpFSqpt(1,iFSqpt))**2 + &
&       (kpt(2)-tmpFSqpt(2,iFSqpt))**2 + &
&       (kpt(3)-tmpFSqpt(3,iFSqpt))**2
       if (ss < tol6) then
         FStoqpt(ikpt1,ikpt2) = iFSqpt
         new=0
         exit
       end if
     end do
     if (new == 1) then
       nFSqpt=nFSqpt+1
       tmpFSqpt(:,nFSqpt) = kpt(:)
       FStoqpt(ikpt1,ikpt2) = nFSqpt
     end if

   end do
 end do

!got nFSqpt,tmpFSqpt,FStoqpt

end subroutine mkfsqgrid
!!***

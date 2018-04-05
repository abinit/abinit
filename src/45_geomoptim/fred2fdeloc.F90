!{\src2tex{textfont=tt}}
!!****f* ABINIT/fred2fdeloc
!! NAME
!! fred2fdeloc
!!
!! FUNCTION
!!  calculate delocalized forces from reduced coordinate ones
!!
!! COPYRIGHT
!! Copyright (C) 2003-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! btinv(3*(natom-1),3*natom)= inverse transpose of B matrix (see delocint)
!! natom = number of atoms
!! gprimd(3,3)=dimensional translations in reciprocal space (bohr-1)
!!
!! OUTPUT
!! deloc_force(3*(natom-1))=delocalized forces from reduced coordinate ones
!! fred(3,natom)=delocalized forces in reduced coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      pred_delocint,xfh_recover_deloc
!!
!! CHILDREN
!!      dgemm,dgemv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fred2fdeloc(btinv,deloc_force,fred,natom,gprimd)

 use defs_basis
 use m_profiling_abi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fred2fdeloc'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 real(dp),intent(in) :: btinv(3*(natom-1),3*natom),gprimd(3,3),fred(3,natom)
 real(dp),intent(out) :: deloc_force(3*(natom-1))

!Local variables-------------------------------
 integer :: ii
!arrays
 real(dp) :: fcart(3,natom)
 character(len=500) :: message

! ******************************************************************

!make cartesian forces

 call dgemm('N','N',3,natom,3,one,&
& gprimd,3,fred,3,zero,fcart,3)

!turn cartesian to delocalized forces
 call dgemv('N',3*(natom-1),3*natom,one,&
& btinv,3*(natom-1),fcart,1,zero,deloc_force,1)

 write (message,'(a)') 'fred2fdeloc : deloc_force = '
 call wrtout(std_out,message,'COLL')

 do ii = 1, 3*(natom-1)
   write (message,'(I6,E16.6)') ii, deloc_force(ii)
   call wrtout(std_out,message,'COLL')
 end do

end subroutine fred2fdeloc
!!***

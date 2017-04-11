!{\src2tex{textfont=tt}}
!!****f* ABINIT/order_fs_kpts
!!
!! NAME
!! order_fs_kpts
!!
!! FUNCTION
!! This routine re-orders the kpoints on the standard grid which belong
!!  to the Fermi surface: put them in increasing z, then y,  then x
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   nkptirr = number of irreducible FS kpoints
!!   nkpt = input nkpt from header
!!   kptns = input kpt from header
!!
!! OUTPUT
!!   FSirredtoGS = mapping of irreducible kpoints to GS set
!!   kptirr = irreducible FS kpoint coordinates
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      destroy_kptrank,mkkptrank,wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine order_fs_kpts(kptns, nkpt, kptirr,nkptirr,FSirredtoGS)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_kptrank

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'order_fs_kpts'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkptirr
 integer,intent(in) :: nkpt

!arrays
 integer,intent(out) :: FSirredtoGS(nkptirr)
 real(dp),intent(in) :: kptns(3,nkpt)
 real(dp),intent(out) :: kptirr(3,nkptirr)

!Local variables-------------------------------
!scalars
 integer :: ikpt,jkpt,kkpt,new, ik
 real(dp) :: res
 type(kptrank_type) :: kptrank_t
!arrays
 integer :: kptirrank(nkptirr)

! *************************************************************************

!rank is used to order kpoints
 call mkkptrank (kptns,nkpt,kptrank_t)

 ik=1
 do ikpt=1,nkpt
!  add kpt to FS kpts, in order, increasing z, then y, then x !
   new = 1
!  look for position to insert kpt ikpt among irredkpts already found
   do jkpt=1,ik-1
     if (kptirrank(jkpt) > kptrank_t%rank(ikpt)) then
!      shift all the others up
       do kkpt=ik-1,jkpt,-1
         kptirr(:,kkpt+1) = kptirr(:,kkpt)
         kptirrank(kkpt+1) = kptirrank(kkpt)
         FSirredtoGS(kkpt+1) = FSirredtoGS(kkpt)
       end do
!      insert kpoint ikpt
       call wrap2_pmhalf(kptns(1,ikpt),kptirr(1,jkpt),res)
       call wrap2_pmhalf(kptns(2,ikpt),kptirr(2,jkpt),res)
       call wrap2_pmhalf(kptns(3,ikpt),kptirr(3,jkpt),res)

       kptirrank(jkpt) = kptrank_t%rank(ikpt)
       FSirredtoGS(jkpt) = ikpt
       new=0
       exit
     end if
   end do
!  ikpt not counted yet and higher rank than all previous
   if (new == 1) then
     call wrap2_pmhalf(kptns(1,ikpt),kptirr(1,ikpt),res)
     call wrap2_pmhalf(kptns(2,ikpt),kptirr(2,ikpt),res)
     call wrap2_pmhalf(kptns(3,ikpt),kptirr(3,ikpt),res)
     kptirrank(ik) = kptrank_t%rank(ikpt)
     FSirredtoGS(ik) = ikpt
   end if
   ik=ik+1
 end do

 call destroy_kptrank (kptrank_t)

end subroutine order_fs_kpts
!!***

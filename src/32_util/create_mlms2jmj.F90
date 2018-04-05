!{\src2tex{textfont=tt}}
!!****f* ABINIT/create_mlms2jmj
!! NAME
!! create_mlms2jmj
!!
!! FUNCTION
!! For a given angular momentum lcor, give the rotation matrix msml2jmj
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum
!!
!! SIDE EFFECTS
!!  mlms2jmj= rotation matrix
!!
!! NOTES
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

subroutine create_mlms2jmj(lcor,mlmstwojmj)

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'create_mlms2jmj'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
!arrays
 complex(dpc),intent(out) :: mlmstwojmj(2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 integer :: jc1,jj,jm,ll,ml1,ms1
 real(dp) :: invsqrt2lp1,xj,xmj
 character(len=500) :: message
!arrays
 integer, allocatable :: ind_msml(:,:)
 complex(dpc),allocatable :: mat_mlms2(:,:)

!*********************************************************************

!--------------- Built indices + allocations
 ll=lcor
 mlmstwojmj=czero
 ABI_ALLOCATE(ind_msml,(2,-ll:ll))
 ABI_ALLOCATE(mat_mlms2,(2*(2*lcor+1),2*(2*lcor+1)))
 mlmstwojmj=czero
 jc1=0
 do ms1=1,2
   do ml1=-ll,ll
     jc1=jc1+1
     ind_msml(ms1,ml1)=jc1
   end do
 end do
!--------------- built mlmstwojmj
!do jj=ll,ll+1    ! the physical value of j are ll-0.5,ll+0.5
!xj(jj)=jj-0.5
 if(ll==0)then
   message=' ll should not be equal to zero !'
   MSG_BUG(message)
 end if
 jc1=0
 invsqrt2lp1=one/sqrt(float(2*lcor+1))
 do jj=ll,ll+1
   xj=float(jj)-half
   do jm=-jj,jj-1
     xmj=float(jm)+half
     jc1=jc1+1
     if(nint(xj+0.5)==ll+1) then
       if(nint(xmj+0.5)==ll+1)  then
         mlmstwojmj(ind_msml(2,ll),jc1)=1.0   !  J=L+0.5 and m_J=L+0.5
       else if(nint(xmj-0.5)==-ll-1) then
         mlmstwojmj(ind_msml(1,-ll),jc1)=1.0   !  J=L+0.5 and m_J=-L-0.5
       else
         mlmstwojmj(ind_msml(2,nint(xmj-0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
         mlmstwojmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
       end if
     end if
     if(nint(xj+0.5)==ll) then
       mlmstwojmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
       mlmstwojmj(ind_msml(2,nint(xmj-0.5)),jc1)=-invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
     end if
   end do
 end do

 ABI_DEALLOCATE(ind_msml)
 ABI_DEALLOCATE(mat_mlms2)

 end subroutine create_mlms2jmj
!!***

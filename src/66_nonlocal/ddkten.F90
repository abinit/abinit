!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddkten
!! NAME
!! ddkten
!!
!! FUNCTION
!! Compact or decompact the tensors related to the ffnl(:,1,...)
!! part of the ddk operator, taking into account the direction
!! of the ddk perturbation.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  compact= if 1, compact from tmpfac
!!  idir=direction of the ddk perturbation
!!  rank=0,1,2, or 3 = rank of tmpfac tensor, also angular momentum (=l)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  temp(2,(rank*(rank+1))/2)=compacted tensor
!!    for l=1, just a scalar
!!    for l=2, a vector
!!  tmpfac(2,(rank+1)*(rank+2)/2)=decompacted tensor
!!    for l=1, a vector
!!    for l=2, a symmetric matrix, stored as
!!     (1 . .)
!!     (6 2 .)
!!     (5 4 3)
!!
!! NOTES
!! For l=0, there is no contribution.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ddkten(compact,idir,rank,temp,tmpfac)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddkten'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: compact,idir,rank
!arrays
 real(dp),intent(inout) :: temp(2,(rank*(rank+1))/2)
 real(dp),intent(inout) :: tmpfac(2,((rank+1)*(rank+2))/2)

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

 if(rank/=1 .and. rank/=2 .and. rank/=3)then
   write(message, '(a,i10,a,a,a)' )&
&   'Input rank=',rank,' not allowed.',ch10,&
&   'Possible values are 1,2,3 only.'
   MSG_BUG(message)
 end if

!Take care of p angular momentum
 if(rank==1)then

!  Compaction tmpfac -> temp
   if(compact==1)then
     temp(:,1)=tmpfac(:,idir)

!    Decompaction temp -> tmpfac
   else
     tmpfac(:,1:3)=0.0d0
     tmpfac(:,idir)=temp(:,1)
   end if

!  Take care of d angular momentum
!  rank=2 11->1 22->2 33->3 32->4 31->5 21->6

 else if(rank==2)then

!  Compaction tmpfac -> temp
   if(compact==1)then
     if(idir==1)then
!      Count the number of non-zero derivatives with respect to k(idir)
!      The factor of 2 on the diagonal comes from the derivative with
!      respect to the first K then to the second K
       temp(:,1)=2.0d0*tmpfac(:,1); temp(:,2)=tmpfac(:,6); temp(:,3)=tmpfac(:,5)
     else if(idir==2)then
       temp(:,2)=2.0d0*tmpfac(:,2); temp(:,1)=tmpfac(:,6); temp(:,3)=tmpfac(:,4)
     else if(idir==3)then
       temp(:,3)=2.0d0*tmpfac(:,3); temp(:,1)=tmpfac(:,5); temp(:,2)=tmpfac(:,4)
     end if
!    Decompaction temp -> tmpfac
   else
     tmpfac(:,1:6)=0.0d0
     tmpfac(:,idir)=2.0d0*temp(:,idir)
     if(idir==1)then
       tmpfac(:,5)=temp(:,3); tmpfac(:,6)=temp(:,2)
     else if(idir==2)then
       tmpfac(:,4)=temp(:,3); tmpfac(:,6)=temp(:,1)
     else if(idir==3)then
       tmpfac(:,4)=temp(:,2); tmpfac(:,5)=temp(:,1)
     end if
   end if

!  Take care of f angular momentum
 else if(rank==3)then
!  rank=3 111->1 221->2 331->3 321->4 311->5 211->6 222->7 332->8 322->9 333->10
!  rank=2 11->1 22->2 33->3 32->4 31->5 21->6

!  Compaction tmpfac -> temp
   if(compact==1)then
     if(idir==1)then
!      Count the number of non-zero derivatives with respect to k(idir)
       temp(:,1)=3.0d0*tmpfac(:,1)
       temp(:,2:4)=tmpfac(:,2:4)
       temp(:,5:6)=2.0d0*tmpfac(:,5:6)
     else if(idir==2)then
       temp(:,6)=2.0d0*tmpfac(:,2)
       temp(:,4)=2.0d0*tmpfac(:,9)
       temp(:,5)=tmpfac(:,4)
       temp(:,1)=tmpfac(:,6)
       temp(:,3)=tmpfac(:,8)
       temp(:,2)=3.0d0*tmpfac(:,7)
     else if(idir==3)then
       temp(:,3)=3.0d0*tmpfac(:,10)
       temp(:,5)=2.0d0*tmpfac(:,3)
       temp(:,4)=2.0d0*tmpfac(:,8)
       temp(:,6)=tmpfac(:,4)
       temp(:,1)=tmpfac(:,5)
       temp(:,2)=tmpfac(:,9)
     end if
!    Decompaction temp -> tmpfac
   else
     tmpfac(:,1:10)=0.0d0
     if(idir==1)then
       tmpfac(:,1)=3.0d0*temp(:,1)
       tmpfac(:,2:4)=temp(:,2:4)
       tmpfac(:,5:6)=2.0d0*temp(:,5:6)
     else if(idir==2)then
       tmpfac(:,2)=2.0d0*temp(:,6)
       tmpfac(:,9)=2.0d0*temp(:,4)
       tmpfac(:,4)=temp(:,5)
       tmpfac(:,6)=temp(:,1)
       tmpfac(:,8)=temp(:,3)
       tmpfac(:,7)=3.0d0*temp(:,2)
     else if(idir==3)then
       tmpfac(:,10)=3.0d0*temp(:,3)
       tmpfac(:,3)=2.0d0*temp(:,5)
       tmpfac(:,8)=2.0d0*temp(:,4)
       tmpfac(:,4)=temp(:,6)
       tmpfac(:,5)=temp(:,1)
       tmpfac(:,9)=temp(:,2)
     end if
   end if

 end if

end subroutine ddkten
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/dist2
!! NAME
!!  dist2 
!!
!! FUNCTION
!!  dist2(v1,v2,rprimd,option) calculates the distance of v1 and v2 in a crystal by
!!  repeating the unit cell 
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2017 ABINIT group (DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  v1,v2
!!  rprimd: dimensions of the unit cell. if not given 1,0,0/0,1,0/0,0,1 is assumed
!!  option: 0 v1, v2 given in cartesian coordinates (default) / 1 v1,v2 given in reduced coordinates
!!
!! OUTPUT
!!  dist2
!!
!! PARENTS
!!  ioniondist
!!
!! CHILDREN
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function dist2(v1,v2,rprimd,option)

 use defs_basis

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dist2'
 use interfaces_41_geometry, except_this_one => dist2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional                      :: option
 real(dp)                                         :: dist2
!arrays
 real(dp),intent(in),optional      :: rprimd(3,3)
 real(dp),intent(in)                  :: v1(3),v2(3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,opt,s1,s2,s3
 real(dp):: min2,norm2,ucvol
!arrays
 integer :: limits(3)
 real(dp) :: corner(3),dred(3),dtot(3),dv(3),dwrap(3),sh(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: vprimd(3,3)

! *************************************************************************

 if (.not.PRESENT(rprimd)) then
   vprimd=reshape((/1,0,0,  0,1,0,  0,0,1/),(/3,3/))
 else
   vprimd=rprimd
 end if

 call metric(gmet,gprimd,-1,rmet,vprimd,ucvol)

 dv(:)=v2(:)-v1(:)

!If in cartesian coordinates, need to be transformed to reduced coordinates.
 opt=0
 if(present(option))then
   opt=option
 end if
 if(opt==0)then
   dred(:)=gprimd(1,:)*dv(1)+gprimd(2,:)*dv(2)+gprimd(3,:)*dv(3)
 else
   dred(:)=dv(:)
 end if

!Wrap in the ]-1/2,1/2] interval
 call wrap2_pmhalf(dred(1),dwrap(1),sh(1))
 call wrap2_pmhalf(dred(2),dwrap(2),sh(2))
 call wrap2_pmhalf(dred(3),dwrap(3),sh(3))

!Compute the limits of the parallelipipedic box that contains the Wigner-Seitz cell
!The reduced coordinates of the corners of the Wigner-Seitz cell are computed (multiplied by two)
!Then, the maximal values of these reduced coordinates are stored.
 limits(:)=0
 do s1=-1,1,2
   do s2=-1,1,2
     do s3=-1,1,2
       corner(:)=gmet(:,1)*s1*rmet(1,1)+gmet(:,2)*s2*rmet(2,2)+gmet(:,3)*s3*rmet(3,3) 
       limits(1)=max(limits(1),ceiling(abs(corner(1))+tol14))
       limits(2)=max(limits(2),ceiling(abs(corner(2))+tol14))
       limits(3)=max(limits(3),ceiling(abs(corner(3))+tol14))
     end do
   end do
 end do

!Use all relevant primitive real space lattice vectors to find the minimal difference vector
 min2=huge(zero)
 do i1=-limits(1),limits(1)
   do i2=-limits(2),limits(2)
     do i3=-limits(3),limits(3)
       dtot(1)=dwrap(1)+i1
       dtot(2)=dwrap(2)+i2
       dtot(3)=dwrap(3)+i3
       norm2=dtot(1)*rmet(1,1)*dtot(1)+dtot(2)*rmet(2,2)*dtot(2)+dtot(3)*rmet(3,3)*dtot(3)+&
&       2*(dtot(1)*rmet(1,2)*dtot(2)+dtot(2)*rmet(2,3)*dtot(3)+dtot(3)*rmet(3,1)*dtot(1))
       min2=min(norm2,min2)
     end do
   end do
 end do
 dist2=sqrt(min2)

 end function dist2
!!***

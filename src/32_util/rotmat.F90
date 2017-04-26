!{\src2tex{textfont=tt}}
!!****f* ABINIT/rotmat
!! NAME
!! rotmat
!!
!! FUNCTION
!! Finds the rotation matrix.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (TRangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  xaxis(3)= vectors defining the x axis
!!  zaxis(3)= vectors defining the z axis
!! OUTPUT
!!  inversion_flag = flag that indicates that an inversion operation
!!   on the coordinate system should be done
!!  umat(3,3)= matrix that rotates the x=(1 0 0) and z=(0 0 1) to the new
!!   values defined in xaxis and zaxis
!! SIDE EFFECTS
!!
!! NOTES
!! Here I set that the axe x is originally at the 1 0 0 direction
!! and z is originally 0 0 1.
!! So calling rotmat(x',z') will find the rotation 
!! matrix for the case in which we rotate the x and z
!! axes from their default values to x' and z'.
!!
!! PARENTS
!!      mlwfovlp_ylmfac,mlwfovlp_ylmfar
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rotmat(xaxis,zaxis,inversion_flag,umat)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotmat'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: inversion_flag
!arrays
 real(dp),intent(in) :: xaxis(3),zaxis(3)
 real(dp),intent(out) :: umat(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: cosine,xmod,zmod
 character(len=500) :: message
!arrays
 real(dp) :: yaxis(3)

! *************************************************************************
 
 xmod = sqrt(xaxis(1)**2 + xaxis(2)**2 + xaxis(3)**2)
 zmod = sqrt(zaxis(1)**2 + zaxis(2)**2 + zaxis(3)**2)

 if(xmod < 1.d-8)then
   write(message,'(a,a,a,i6)')&
&   'The module of the xaxis should be greater than 1.d-8,',ch10,&
&   'however, |xaxis|=',xmod
   MSG_BUG(message)
 end if


 if(zmod < 1.d-8)then
   write(message,'(a,a,a,i6)')&
&   'The module of the zaxis should be greater than 1.d-8,',ch10,&
&   'however, |zaxis|=',zmod
   MSG_ERROR(message)
 end if

!verify that both axis are perpendicular
 cosine = (xaxis(1)*zaxis(1) + xaxis(2)*zaxis(2) &
& + xaxis(3)*zaxis(3))/(xmod*zmod)

 if(abs(cosine) > 1.d-8)then
   write(message,'(a,a,a,i6)')&
&   'xaxis and zaxis should be perpendicular,',ch10,&
&   'however, cosine=',cosine
   MSG_BUG(message)
 end if

!new y axis as cross product
 yaxis(1) = (zaxis(2)*xaxis(3) - xaxis(2)*zaxis(3))/(xmod*zmod)
 yaxis(2) = (zaxis(3)*xaxis(1) - xaxis(3)*zaxis(1))/(xmod*zmod)
 yaxis(3) = (zaxis(1)*xaxis(2) - xaxis(1)*zaxis(2))/(xmod*zmod)

!hack to allow inversion operation on coordinate transformation
!uses unlikely large but legal values of proj_x and/or proj_z
!to flag inversion
 inversion_flag=0
 if(xmod>10._dp .or. zmod>10._dp) then
   inversion_flag=1
   write(message, '(4a)' )&
&   'inversion operation will be appended to axis transformation',ch10,&
&   'Action : If you did not intend this, make |z|<10 and |x|<10 ',ch10
   call wrtout(std_out,message,'COLL')
 end if


 umat(1,:) = xaxis(:)/xmod
 umat(2,:) = yaxis(:)
 umat(3,:) = zaxis(:)/zmod

end subroutine rotmat
!!***

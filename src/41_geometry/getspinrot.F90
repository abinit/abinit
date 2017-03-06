!{\src2tex{textfont=tt}}
!!****f* ABINIT/getspinrot
!! NAME
!! getspinrot
!!
!! FUNCTION
!! From the symmetry matrix symrel_conv expressed
!! in the coordinate system rprimd,
!! compute the components of the spinor rotation matrix
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel_conv(3,3)=symmetry operation in real space in terms
!!  of primitive translations rprimd
!!
!! OUTPUT
!! spinrot(4)=components of the spinor rotation matrix :
!!  spinrot(1)=$\cos \phi / 2$
!!  spinrot(2)=$\sin \phi / 2 \times u_x$
!!  spinrot(3)=$\sin \phi / 2 \times u_y$
!!  spinrot(4)=$\sin \phi / 2 \times u_z$
!!  where $\phi$ is the angle of rotation, and
!!  $(u_x,u_y,u_z)$ is the normalized direction of the rotation axis
!!
!! NOTES
!! Only the proper part of the symmetry operation is taken into account :
!! pure rotations, while the inversion part is taken away, if present.
!!
!! The whole collection of symmetry matrices is call symrel(3,3,nsym)
!! symrel1 contains just one of those matrices symrel1(3,3)
!!
!! PARENTS
!!      cg_rotate,m_crystal,wfconv
!!
!! CHILDREN
!!      mati3det,matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getspinrot(rprimd,spinrot,symrel_conv)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getspinrot'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: symrel_conv(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: spinrot(4)

!Local variables-------------------------------
!scalars
 integer :: det,ii
 real(dp) :: cos_phi,norminv,phi,scprod,sin_phi
 !character(len=500) :: message
!arrays
 integer :: identity(3,3),symrel1(3,3)
 real(dp) :: axis(3),coord(3,3),coordinvt(3,3),matr1(3,3),matr2(3,3)
 real(dp) :: rprimd_invt(3,3),vecta(3),vectb(3),vectc(3)

!**************************************************************************

 symrel1(:,:)=symrel_conv(:,:)

!Compute determinant of the matrix
 call mati3det(symrel1,det)

!Produce a rotation from an improper symmetry
 if(det==-1)symrel1(:,:)=-symrel1(:,:)

!Test the possibility of the unit matrix
 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1

 if( sum((symrel1(:,:)-identity(:,:))**2)/=0)then

!  Transform symmetry matrix in the system defined by rprimd
   call matr3inv(rprimd,rprimd_invt)
   do ii=1,3
     coord(:,ii)=rprimd_invt(ii,:)
   end do
   call matr3inv(coord,coordinvt)
   do ii=1,3
     matr1(:,ii)=symrel1(:,1)*coord(1,ii)+&
&     symrel1(:,2)*coord(2,ii)+&
&     symrel1(:,3)*coord(3,ii)
   end do
   do ii=1,3
     matr2(:,ii)=coordinvt(1,:)*matr1(1,ii)+&
&     coordinvt(2,:)*matr1(2,ii)+&
&     coordinvt(3,:)*matr1(3,ii)
   end do

!  Find the eigenvector with unit eigenvalue of the
!  rotation matrix in cartesian coordinate, matr2

   matr1(:,:)=matr2(:,:)
   matr1(1,1)=matr1(1,1)-one
   matr1(2,2)=matr1(2,2)-one
   matr1(3,3)=matr1(3,3)-one

!  Compute the axis of rotation and the cos and sin of rotation angle
   if(matr1(1,1)**2 + matr1(2,1)**2 + matr1(3,1)**2 < tol8 )then
!    The first direction is the axis
     axis(1)=one ; axis(2)=zero ; axis(3)=zero
     cos_phi=matr2(2,2)
     sin_phi=matr2(3,2)
   else if(matr1(1,2)**2 + matr1(2,2)**2 + matr1(3,2)**2 < tol8 )then
!    The second direction is the axis
     axis(1)=zero ; axis(2)=one ; axis(3)=zero
     cos_phi=matr2(3,3)
     sin_phi=matr2(1,3)
   else
!    In this case, try use the first and second vector to build the
!    rotation axis : compute their cross product
     axis(1)=matr1(2,1)*matr1(3,2)-matr1(2,2)*matr1(3,1)
     axis(2)=matr1(3,1)*matr1(1,2)-matr1(3,2)*matr1(1,1)
     axis(3)=matr1(1,1)*matr1(2,2)-matr1(1,2)*matr1(2,1)
!    Then, try to normalize it
     scprod=sum(axis(:)**2)
     if(scprod<tol8)then
!      The first and second vectors were linearly dependent
!      Thus, use the first and third vectors
       axis(1)=matr1(2,1)*matr1(3,3)-matr1(2,3)*matr1(3,1)
       axis(2)=matr1(3,1)*matr1(1,3)-matr1(3,3)*matr1(1,1)
       axis(3)=matr1(1,1)*matr1(2,3)-matr1(1,3)*matr1(2,1)
!      Normalize the vector
       scprod=sum(axis(:)**2)
       if(scprod<tol8)then
         MSG_BUG('Cannot find the rotation axis.')
       end if
     end if
     norminv=one/sqrt(scprod)
     axis(:)=axis(:)*norminv

!    Project the axis vector out of the first unit vector,
!    and renormalize the projected vector
!    (the first vector cannot be the axis, as tested before)
     vecta(1)=one-axis(1)**2
     vecta(2)=-axis(1)*axis(2)
     vecta(3)=-axis(1)*axis(3)
     scprod=sum(vecta(:)**2)
     norminv=one/sqrt(scprod)
     vecta(:)=vecta(:)*norminv
!    Rotate the vector A, to get vector B
     vectb(:)=matr2(:,1)*vecta(1)+matr2(:,2)*vecta(2)+matr2(:,3)*vecta(3)
!    Get dot product of vectors A and B, giving cos of the rotation angle
     cos_phi=sum(vecta(:)*vectb(:))
!    Compute the cross product of the axis and vector A
     vectc(1)=axis(2)*vecta(3)-axis(3)*vecta(2)
     vectc(2)=axis(3)*vecta(1)-axis(1)*vecta(3)
     vectc(3)=axis(1)*vecta(2)-axis(2)*vecta(1)
!    Get dot product of vectors B and C, giving sin of the rotation angle
     sin_phi=sum(vectb(:)*vectc(:))
   end if

!  Get the rotation angle, then the parameters of the spinor rotation
!  Here, treat possible inaccurate values of the cosine of phi
   if(cos_phi>  one-tol8 )cos_phi=  one-tol8
   if(cos_phi<-(one-tol8))cos_phi=-(one-tol8)
   phi=acos(cos_phi)
   if(sin_phi<zero)phi=-phi
!  Rectify the angle, such that its absolute values corresponds to
!  180, 120, 90, 60, or 0 degrees
   phi=(nint(six*phi/pi))/six*pi
!  Compute components of the spinor matrix
   spinrot(1)=cos(half*phi)
   spinrot(2)=axis(1)*sin(half*phi)
   spinrot(3)=axis(2)*sin(half*phi)
   spinrot(4)=axis(3)*sin(half*phi)

 else

!  Here, the case of the unit matrix
   axis(:)=zero
   phi=zero
   spinrot(1)=one
   spinrot(2)=zero
   spinrot(3)=zero
   spinrot(4)=zero

 end if ! the case of the identity matrix

!DEBUG
!write(std_out,*)' getspinrot :'
!write(std_out,*)' symrel_conv =',symrel_conv(:,:)
!write(std_out,*)' symrel =',symrel1(:,:)
!write(std_out,*)' rprimd =',rprimd(:,:)
!write(std_out,*)' matr2 =',matr2(:,:)
!write(std_out,*)' matr1 =',matr1(:,:)
!write(std_out,*)' phi (degree)=',phi*180._dp/pi
!write(std_out,'(a,3d16.6)' )' axis=',axis(:)
!write(std_out,*)' vecta=',vecta(:)
!stop
!ENDDEBUG

end subroutine getspinrot
!!***

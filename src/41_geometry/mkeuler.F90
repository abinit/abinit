!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkeuler
!! NAME
!! mkeuler
!!
!! FUNCTION
!! For a given symmetry operation, determines the corresponding Euler angles
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rot(3,3)= symmetry matrix
!!
!! OUTPUT
!!  cosalp=  cos(alpha) with alpha=Euler angle 1
!!  cosbeta= cos(beta)  with beta =Euler angle 2
!!  cosgam=  cos(gamma) with gamma=Euler angle 3
!!  isn= error code (0 if the routine exit normally)
!!  sinalp= sin(alpha) with alpha=Euler angle 1
!!  singam= sin(gamma) with gamma=Euler angle 3
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!      setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkeuler(rot,cosbeta,cosalp,sinalp,cosgam,singam,isn)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkeuler'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(out) :: isn
 real(dp),intent(out) :: cosalp,cosbeta,cosgam,sinalp,singam
!arrays
 real(dp),intent(in) :: rot(3,3)

!Local variables ---------------------------------------
!scalars
 integer :: ier
 real(dp) :: check,sinbeta
 character(len=500) :: message

! *********************************************************************

 do isn= -1,1,2
   cosbeta=real(isn)*rot(3,3)
   if(abs(1._dp-cosbeta*cosbeta)<tol10) then
     sinbeta=zero
   else
     sinbeta=sqrt(1._dp-cosbeta*cosbeta)
   end if
   if (abs(sinbeta).gt.tol10)  then
     cosalp=isn*rot(3,1)/sinbeta
     sinalp=isn*rot(3,2)/sinbeta
     cosgam=-isn*rot(1,3)/sinbeta
     singam=isn*rot(2,3)/sinbeta
   else
     cosalp=isn*rot(1,1)/cosbeta
     sinalp=isn*rot(1,2)/cosbeta
     cosgam=one
     singam=zero
   end if

!  Check matrix:
   ier=0
   check=cosalp*cosbeta*cosgam-sinalp*singam
   if (abs(check-isn*rot(1,1))>tol8) ier=ier+1
   check=sinalp*cosbeta*cosgam+cosalp*singam
   if (abs(check-isn*rot(1,2))>tol8) ier=ier+1
   check=-sinbeta*cosgam
   if (abs(check-isn*rot(1,3))>tol8) ier=ier+1
   check=-cosalp*cosbeta*singam-sinalp*cosgam
   if (abs(check-isn*rot(2,1))>tol8) ier=ier+1
   check=-sinalp*cosbeta*singam+cosalp*cosgam
   if (abs(check-isn*rot(2,2))>tol8) ier=ier+1
   check=sinbeta*singam
   if (abs(check-isn*rot(2,3))>tol8) ier=ier+1
   check=cosalp*sinbeta
   if (abs(check-isn*rot(3,1))>tol8) ier=ier+1
   check=sinalp*sinbeta
   if (abs(check-isn*rot(3,2))>tol8) ier=ier+1
   if (ier.eq.0) return
 end do

 isn=0
 write(message, '(7a)' )&
& 'Error during determination of symetries !',ch10,&
& 'Action: check your input file:',ch10,&
& '        unit cell vectors and/or atoms positions',ch10,&
& '        have to be given with a better precision...'
 MSG_ERROR(message)

end subroutine mkeuler
!!***

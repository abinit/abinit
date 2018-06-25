!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_angles
!! NAME
!! m_angles
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (BA, NH, FJ, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_angles

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_special_funcs, only : factorial

 implicit none

 private
!!***

 public :: mkeuler             ! For a given symmetry operation, determines the corresponding Euler angles
 public :: dbeta               ! Calculate the rotation matrix d^l_{m{\prim}m}(beta)
 public :: create_slm2ylm      ! For a given angular momentum lcor, compute slm2ylm
 public :: create_mlms2jmj     ! For a given angular momentum lcor, give the rotation matrix msml2jmj


contains
!!***

!!****f* m_angles/mkeuler
!! NAME
!! mkeuler
!!
!! FUNCTION
!! For a given symmetry operation, determines the corresponding Euler angles
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
!!  This file comes from the file crystal_symmetry.f
!!  by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!      setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkeuler(rot,cosbeta,cosalp,sinalp,cosgam,singam,isn)


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

!!****f* m_angles/dbeta
!! NAME
!! dbeta
!!
!! FUNCTION
!!  Calculate the rotation matrix d^l_{m{\prim}m}(beta) using Eq. 4.14 of
!!  M.E. Rose, Elementary Theory of Angular Momentum,
!!             John Wiley & Sons, New-York, 1957
!!
!! INPUTS
!!  cosbeta= cosinus of beta (=Euler angle)
!!  ll= index l
!!  mm= index m
!!  mp= index m_prime
!!
!! OUTPUT
!!  dbeta= rotation matrix
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!  - Assume l relatively small so that factorials do not cause
!!    roundoff error
!!
!! PARENTS
!!     setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

function dbeta(cosbeta,ll,mp,mm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dbeta'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,mm,mp
 real(dp) :: dbeta
 real(dp),intent(in) :: cosbeta

!Local variables ------------------------------
!scalars
 integer,parameter :: mxterms=200
 integer :: ii,ina,inb,inc,ml,ms
 real(dp) :: arg,cosbetab2,pref,sinbetab2,sum,tt

!************************************************************************
 dbeta=zero

!Special cases
 if (abs(cosbeta-1._dp).lt.tol10) then
   if (mp.eq.mm) dbeta=1
 else if (abs(cosbeta+1._dp).lt.tol10) then
   if (mp.eq.-mm) dbeta=(-1)**(ll+mm)
 else

!  General case
   cosbetab2=sqrt((1+cosbeta)*0.5_dp)
   sinbetab2=sqrt((1-cosbeta)*0.5_dp)
   ml=max(mp,mm)
   ms=min(mp,mm)
   if (ml.ne.mp) sinbetab2=-sinbetab2
   tt=-(sinbetab2/cosbetab2)**2
   pref=sqrt((factorial(ll-ms)*factorial(ll+ml))&
&   /(factorial(ll+ms)*factorial(ll-ml)))&
&   /factorial(ml-ms)*(cosbetab2**(2*ll+ms-ml))&
&   *((-sinbetab2)**(ml-ms))
   sum=1._dp
   arg=1._dp
   ina=ml-ll
   inb=-ms-ll
   inc=ml-ms+1
   do ii=1,mxterms
     if (ina.eq.0.or.inb.eq.0) exit
     arg=(arg*ina*inb*tt)/(ii*inc)
     sum=sum+arg
     ina=ina+1
     inb=inb+1
     inc=inc+1
   end do
   dbeta=pref*sum
 end if

end function dbeta
!!***

!!****f* m_angles/create_slm2ylm
!! NAME
!! create_slm2ylm
!!
!! FUNCTION
!! For a given angular momentum lcor, compute slm2ylm
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!
!! OUTPUT
!!  slm2ylm(2lcor+1,2lcor+1) = rotation matrix.
!!
!! NOTES
!!  useful only in ndij==4
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_slm2ylm(lcor,slmtwoylm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'create_slm2ylm'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
!arrays
 complex(dpc),intent(out) :: slmtwoylm(2*lcor+1,2*lcor+1)

!Local variables ---------------------------------------
!scalars
 integer :: jm,ll,mm,im
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
!arrays

! *********************************************************************

 ll=lcor
 slmtwoylm=czero
 do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1
   onem=dble((-1)**mm)
   if (mm> 0) then
     slmtwoylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
     slmtwoylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
   end if
   if (mm==0) then
     slmtwoylm(im,im)=cone
   end if
   if (mm< 0) then
     slmtwoylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
     slmtwoylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
   end if
 end do

end subroutine create_slm2ylm
!!***

!!****f* m_angles/create_mlms2jmj
!! NAME
!! create_mlms2jmj
!!
!! FUNCTION
!! For a given angular momentum lcor, give the rotation matrix msml2jmj
!!
!! INPUTS
!!  lcor= angular momentum
!!
!! SIDE EFFECTS
!!  mlms2jmj= rotation matrix
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_mlms2jmj(lcor,mlmstwojmj)


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

end module m_angles
!!***

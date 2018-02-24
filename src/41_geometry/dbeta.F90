!{\src2tex{textfont=tt}}
!!****f* ABINIT/dbeta
!! NAME
!! dbeta
!!
!! FUNCTION
!!  Calculate the rotation matrix d^l_{m{\prim}m}(beta) using Eq. 4.14 of
!!  M.E. Rose, Elementary Theory of Angular Momentum,
!!             John Wiley & Sons, New-York, 1957
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FJ, MT, NH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function dbeta(cosbeta,ll,mp,mm)

 use defs_basis

 use m_special_funcs, only : factorial

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

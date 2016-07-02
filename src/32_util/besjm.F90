!{\src2tex{textfont=tt}}
!!****f* ABINIT/besjm
!! NAME
!! besjm
!!
!! FUNCTION
!! Spherical bessel function of order nn. Handles nn=0,1,2,3,4, or 5 only.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCG, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  arg= scaling to be applied to xx(nx)
!!  nn=order of spherical bessel function (only 0 through 5 allowed)
!!  cosx(1:nx)=cosines of arg*xx(1:nx)
!!  xx(1:nx)=set of dimensionless arguments of function
!!  nx=number of arguments
!!  sinx(1:nx)=sines of arg*xx(1:nx)
!!
!! OUTPUT
!!  besjx(1:nx)=returned values
!!
!! NOTES
!! besj(nn,y)=$ j_{nn}(y) =(\frac{\pi}{2y})^{\frac{1}{2}}J(nn+\frac{1}{2},y)$
!! where J=Bessel function of the first kind.
!! besjm compute multiple values, and relies on precomputed values of sin and cos of y.
!! The argument y is arg*xx(ix), for ix from 1 to nx
!! The values of xx must be positive, and ordered by increasing order
!! At small arg, the higher orders have so much cancellation that the
!! analytic expression is very poor computationally.  In that case we
!! use a rational polynomial approximation.
!!
!! PARENTS
!!      init_bess_spl,mlwfovlp_radial,psp1nl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine besjm(arg,besjx,cosx,nn,nx,sinx,xx)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'besjm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn,nx
 real(dp),intent(in) :: arg
!arrays
 real(dp),intent(in) :: cosx(nx),sinx(nx),xx(nx)
 real(dp),intent(out) :: besjx(nx)

!Local variables-------------------------------
!scalars
 integer :: ix,switchx
!no_abirules
!Series or rational polynomial coefficients
 real(dp),parameter :: b01=1.d0/6.d0,b02=1.d0/120.d0,b03=1.d0/5040.d0
 real(dp),parameter :: b04=1.d0/362880.d0,b11=0.8331251468724171d-1
 real(dp),parameter :: b12=0.2036961284395412d-2,b13=0.1932970379901801d-4
 real(dp),parameter :: b14=0.6526053169009489d-7,b21=0.5867824627555163d-1
 real(dp),parameter :: b22=0.1152501878595934d-2,b23=0.1011071389414764d-4
 real(dp),parameter :: b24=0.4172322111421287d-7,b25=0.6790616688656543d-10
 real(dp),parameter :: b31=0.439131885807176d-1,b32=0.6813139609887099d-3
 real(dp),parameter :: b33=0.4899103784264755d-5,b34=0.17025590795625d-7
 real(dp),parameter :: b35=0.2382642910613347d-10,b41=0.3587477991030971d-1
 real(dp),parameter :: b42=0.4833719855268907d-3,b43=0.3238388977796242d-5
 real(dp),parameter :: b44=0.1171802513125112d-7,b45=0.223261650431992d-10
 real(dp),parameter :: b46=.1800045587335951d-13,b51=0.295232406376567d-1
 real(dp),parameter :: b52=0.3359864457080573d-3,b53=0.19394750603618d-5
 real(dp),parameter :: b54=0.6143166228216219d-8,b55=0.10378501636108d-10
 real(dp),parameter :: b56=.749975122872713d-14
 real(dp),parameter :: c11=0.1668748531275829d-1,c12=0.1342812442426702d-3
 real(dp),parameter :: c13=0.6378249315355233d-6,c14=0.1573564527360138d-8
 real(dp),parameter :: c21=0.127503251530198d-1,c22=0.7911240539893565d-4
 real(dp),parameter :: c23=0.3044380758068054d-6,c24=0.7439837832363479d-9
 real(dp),parameter :: c25=0.9515065658793124d-12,c31=0.1164236697483795d-1
 real(dp),parameter :: c32=0.654858636312224d-4,c33=0.2265576367562734d-6
 real(dp),parameter :: c34=0.4929905563217352d-9,c35=0.555120465710914d-12
 real(dp),parameter :: c41=0.9579765544235745d-2,c42=0.4468999977536864d-4
 real(dp),parameter :: c43=0.1315634305905896d-6,c44=0.2615492488301639d-9
 real(dp),parameter :: c45=0.3387473312408129d-12,c46=.2280866204624012d-15
 real(dp),parameter :: c51=0.8938297823881763d-2,c52=0.3874149021633025d-4
 real(dp),parameter :: c53=0.1054692715135225d-6,c54=0.192879620987602d-9
 real(dp),parameter :: c55=0.2284469423833734d-12,c56=0.139729234332572d-15
 real(dp),parameter :: ffnth=1.d0/15.d0,o10395=1.d0/10395d0,oo105=1.d0/105.d0
 real(dp),parameter :: oo945=1.d0/945.d0
 real(dp) :: bot,rr,rsq,top
 character(len=500) :: message

! *************************************************************************

 if (nn==0) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=1.d-1) then
       rsq=rr*rr
       besjx(ix)=1.d0-rsq*(b01-rsq*(b02-rsq*(b03-rsq*b04)))
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     besjx(ix)=sinx(ix)/rr
   end do

 else if (nn==1) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=1.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b11-rsq*(b12-rsq*(b13-rsq*b14)))
       bot=1.d0+rsq*(c11+rsq*(c12+rsq*(c13+rsq*c14)))
       besjx(ix)=third*rr*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=(sinx(ix)-rr*cosx(ix))/rsq
   end do

 else if (nn==2) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=2.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b21-rsq*(b22-rsq*(b23-rsq*(b24-rsq*b25))))
       bot=1.d0+rsq*(c21+rsq*(c22+rsq*(c23+rsq*(c24+rsq*c25))))
       besjx(ix)=ffnth*rsq*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=((3.d0-rsq)*sinx(ix)-3.d0*rr*cosx(ix))/(rr*rsq)
   end do

 else if (nn==3) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=2.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b31-rsq*(b32-rsq*(b33-rsq*(b34-rsq*b35))))
       bot=1.d0+rsq*(c31+rsq*(c32+rsq*(c33+rsq*(c34+rsq*c35))))
       besjx(ix)=rr*rsq*oo105*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=( (15.d0-6.d0*rsq)*sinx(ix)&
&     + rr*(rsq-15.d0)  *cosx(ix) ) /(rsq*rsq)
   end do

 else if (nn==4) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=4.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b41-rsq*(b42-rsq*(b43-rsq*(b44-rsq*(b45-rsq*b46)))))
       bot=1.d0+rsq*(c41+rsq*(c42+rsq*(c43+rsq*(c44+rsq*(c45+rsq*c46)))))
       besjx(ix)=rsq*rsq*oo945*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=( (105.d0-rsq*(45.d0-rsq)) *sinx(ix)&
&     + rr * (10.d0*rsq-105.d0)  *cosx(ix) ) /(rsq*rsq*rr)
   end do

 else if (nn==5) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=4.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b51-rsq*(b52-rsq*(b53-rsq*(b54-rsq*(b55-rsq*b56)))))
       bot=1.d0+rsq*(c51+rsq*(c52+rsq*(c53+rsq*(c54+rsq*(c55+rsq*c56)))))
       besjx(ix)=rsq*rsq*rr*o10395*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=( (945.d0-rsq*(420.d0-rsq*15.d0)) *sinx(ix)&
&     + rr * (945.d0-rsq*(105.d0-rsq))  *cosx(ix) ) /(rsq*rsq*rr)
   end do

 else
   write(message, '(a,i0,a)' )' besjm only defined for nn in [0,5]; input was nn=',nn,'.'
   MSG_BUG(message)
 end if

end subroutine besjm
!!***

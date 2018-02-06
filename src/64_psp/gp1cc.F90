!{\src2tex{textfont=tt}}
!!****f* ABINIT/gp1cc
!! NAME
!! gp1cc
!!
!! FUNCTION
!! Derivative of gg(xx) wrt xx.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  xx=abscisse to which gp1cc_xx is calculated
!!
!! OUTPUT
!!  gp1cc_xx=derivative of gg(xx) wrt xx.
!!
!! NOTES
!! $ phi(x) = \frac{\sin(2\pi x)}{(2\pi x)(1-4x^2)(1-x^2)}$
!! $ gg(x)= phi(x)^2$
!! $ gp(x)= 2 * phi(x) * phi''(x)$
!! $ phi''(x)=\frac{\cos(2\pi x)-(1-15x^2+20x^4) phi(x)}{x(1-4x^2)(1-x^2)}$
!!
!!
!! PARENTS
!!      psp1cc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gp1cc(gp1cc_xx,xx)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gp1cc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: gp1cc_xx

!Local variables -------------------------------------------
!scalars
 real(dp),parameter :: c11=20.d0-8.d0*pi**2/3.d0
 real(dp),parameter :: c12=268.d0-160.d0/3.d0*pi**2+128.d0/45.d0*pi**4
 real(dp),parameter :: c21=-40.d0/27.d0,c22=40.d0/3.d0-32.d0*pi**2/27.d0
 real(dp),parameter :: c23=-4160.d0/81.d0+160.d0*pi**2/27.d0
 real(dp),parameter :: c24=157712.d0/729.d0-320.d0*pi**2/9.d0+512.d0*pi**4/405.d0
 real(dp),parameter :: c25=-452200.d0/729.d0+83200.d0*pi**2/729.d0-1280.d0*pi**4/243.d0
 real(dp),parameter :: c31=-25.d0/108.d0,c32=485.d0/216.d0-2.d0*pi**2/27.d0
 real(dp),parameter :: c33=-4055.d0/324.d0+25.d0*pi**2/27.d0
 real(dp),parameter :: c34=616697.d0/11664.d0-485.d0*pi**2/81.d0+32.d0*pi**4/405.d0
 real(dp),parameter :: c35=-2933875.d0/15552.d0+20275.d0*pi**2/729.d0-200.d0*pi**4/243.d0
 real(dp),parameter :: two_pim1=1.0d0/two_pi
 real(dp) :: denom,phi,phip

! *************************************************************************

!Cut off beyond r=3*xcccrc is already done at the calling level
 if (xx>1.001d0) then
!  The part that follows will be repeated later, but written in this way,
!  only one "if" condition is tested in most of the cases (1.001 < x < 3.0)
   denom=1.d0/(xx*(1.d0-4.d0*xx**2)*(1.d0-xx**2))
   phi=denom*sin(two_pi*xx)*two_pim1
   phip=denom*(cos(two_pi*xx)-(1.d0-xx**2*(15.d0-xx**2*20))*phi)
   gp1cc_xx=2.d0*phi*phip
!  Handle limits where denominator vanishes
 else if (abs(xx)<1.d-03) then
   gp1cc_xx=xx*(c11+xx**2*c12)
 else if (abs(xx-0.5d0)<=1.d-03) then
   gp1cc_xx=c21+(xx-0.5d0)*(c22+(xx-0.5d0)*(c23+(xx-0.5d0)*(c24+(xx-0.5d0)*c25)))
 else if (abs(xx-1.d0)<=1.d-03) then
   gp1cc_xx=c31+(xx-1.0d0)*(c32+(xx-1.0d0)*(c33+(xx-1.0d0)*(c34+(xx-1.0d0)*c35)))
 else
!  Here is the repeated part ...
   denom=1.d0/(xx*(1.d0-4.d0*xx**2)*(1.d0-xx**2))
   phi=denom*sin(two_pi*xx)*two_pim1
   phip=denom*(cos(two_pi*xx)-(1.d0-xx**2*(15.d0-xx**2*20))*phi)
   gp1cc_xx=2.d0*phi*phip
 end if

end subroutine gp1cc
!!***

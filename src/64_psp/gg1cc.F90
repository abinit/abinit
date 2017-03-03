!{\src2tex{textfont=tt}}
!!****f* ABINIT/gg1cc
!! NAME
!! gg1cc
!!
!! FUNCTION
!! gg1cc_xx=$(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4xx^2)(1-xx^2)})^2$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  xx= abscisse to which gg1cc_xx is calculated
!!
!! OUTPUT
!!  gg1cc_xx= gg1cc_x(xx)
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


subroutine gg1cc(gg1cc_xx,xx)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gg1cc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: gg1cc_xx

!Local variables -------------------------------------------
!The c s are coefficients for Taylor expansion of the analytic
!form near xx=0, 1/2, and 1.
!scalars
 real(dp) :: c21=4.d0/9.d0,c22=-40.d0/27.d0,c23=20.d0/3.d0-16.d0*pi**2/27.d0
 real(dp) :: c24=-4160.d0/243.d0+160.d0*pi**2/81.d0,c31=1.d0/36.d0
 real(dp) :: c32=-25.d0/108.d0,c33=485.d0/432.d0-pi**2/27.d0
 real(dp) :: c34=-4055.d0/972.d0+25.d0*pi**2/81.d0

! *************************************************************************

!Cut off beyond 3/gcut=xcccrc
 if (xx>3.0d0) then
   gg1cc_xx=0.0d0
!  Take care of difficult limits near x=0, 1/2, and 1
 else if (abs(xx)<=1.d-09) then
   gg1cc_xx=1.d0
 else if (abs(xx-0.5d0)<=1.d-04) then
!  (this limit and next are more troublesome for numerical cancellation)
   gg1cc_xx=c21+(xx-0.5d0)*(c22+(xx-0.5d0)*(c23+(xx-0.5d0)*c24))
 else if (abs(xx-1.d0)<=1.d-04) then
   gg1cc_xx=c31+(xx-1.0d0)*(c32+(xx-1.0d0)*(c33+(xx-1.0d0)*c34))
 else
!  The following is the square of the Fourier transform of a
!  function built out of two spherical bessel functions in G
!  space and cut off absolutely beyond gcut
   gg1cc_xx=(sin(2.0d0*pi*xx)/( (2.0d0*pi*xx) * &
&   (1.d0-4.0d0*xx**2)*(1.d0-xx**2) )  )**2
 end if

end subroutine gg1cc
!!***

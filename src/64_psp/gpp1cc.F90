!{\src2tex{textfont=tt}}
!!****f* ABINIT/gpp1cc
!! NAME
!! gpp1cc
!!
!! FUNCTION
!! Second derivative of gg wrt xx.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  xx= abscisse to which gpp1cc_xx is calculated
!!
!! OUTPUT
!!  gpp1cc_xx=second derivative of gg wrt xx.
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


subroutine gpp1cc(gpp1cc_xx,xx)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpp1cc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: gpp1cc_xx

!Local variables -------------------------------------------
!scalars
 real(dp),parameter :: c1=20.d0-8.d0*pi**2/3.d00
 real(dp),parameter :: c2=40.d0/3.d0-32.d0*pi**2/27.d0
 real(dp),parameter :: c3=-8320.d0/81.d0+320.d0*pi**2/27.d0
 real(dp),parameter :: c4=157712.d0/243.d0-320.d0*pi**2/3.d0+512.d0*pi**4/135.d0
 real(dp),parameter :: c5=-18088.d2/729.d0+3328.d2*pi**2/729.d0-5120.d0*pi**4/243.d0
 real(dp),parameter :: c6=485.d0/216.d0-2.d0*pi**2/27.d0
 real(dp),parameter :: c7=-4055.d0/162.d0+50.d0*pi**2/27.d0
 real(dp),parameter :: c8=616697.d0/3888.d0-485.d0*pi**2/27.d0+32.d0*pi**4/135.d0
 real(dp),parameter :: c9=-2933875.d0/3888.d0+81100.d0*pi**2/729.d0-800.d0*pi**4/243.d0
 real(dp) :: t1,t10,t100,t11,t12,t120,t121,t122,t127,t138,t14,t140,t15,t152
 real(dp) :: t157,t16,t160,t17,t174,t175,t18,t19,t2,t20,t21,t23,t24,t3,t31,t33
 real(dp) :: t34,t4,t41,t42,t44,t45,t46,t5,t54,t55,t56,t57,t6,t62,t64,t65,t7
 real(dp) :: t72,t78,t79,t8,t85,t9,t93

! *************************************************************************

 if (xx>3.0d0) then
!  Cut off beyond 3/gcut=3*xcccrc
   gpp1cc_xx=0.0d0
!  Take care of difficult limits near xx=0, 1/2, and 1
 else if (abs(xx)<=1.d-09) then
   gpp1cc_xx=c1
 else if (abs(xx-0.5d0)<=1.d-04) then
!  (this limit and next are more troublesome for numerical cancellation)
   gpp1cc_xx=c2+(xx-0.5d0)*(c3+(xx-0.5d0)*(c4+(xx-0.5d0)*c5))
 else if (abs(xx-1.d0)<=1.d-04) then
   gpp1cc_xx=c6+(xx-1.0d0)*(c7+(xx-1.0d0)*(c8+(xx-1.0d0)*c9))
 else

!  Should fix up this Maple fortran later
   t1 = xx**2
   t2 = 1/t1
   t3 = 1/Pi
   t4 = 2*xx
   t5 = t4-1
   t6 = t5**2
   t7 = 1/t6
   t8 = t4+1
   t9 = t8**2
   t10 = 1/t9
   t11 = xx-1
   t12 = t11**2
   t14 = 1/t12/t11
   t15 = xx+1
   t16 = t15**2
   t17 = 1/t16
   t18 = Pi*xx
   t19 = sin(t18)
   t20 = cos(t18)
   t21 = t20**2
   t23 = t19*t21*t20
   t24 = t17*t23
   t31 = t19**2
   t33 = t31*t19*t20
   t34 = t17*t33
   t41 = Pi**2
   t42 = 1/t41
   t44 = 1/t16/t15
   t45 = t31*t21
   t46 = t44*t45
   t54 = 1/t1/xx
   t55 = 1/t12
   t56 = t55*t46
   t57 = t10*t56
   t62 = t9**2
   t64 = t17*t45
   t65 = t55*t64
   t72 = 1/t9/t8
   t78 = t14*t64
   t79 = t10*t78
   t85 = t12**2
   t93 = t21**2
   t100 = t31**2
   t120 = 1/t6/t5
   t121 = t55*t34
   t122 = t10*t121
   t127 = t16**2
   t138 = t6**2
   t140 = t10*t65
   t152 = t72*t65
   t157 = t7*t140
   t160 = t1**2
   t174 = t55*t24
   t175 = t10*t174
   gpp1cc_xx = 8*t2*t3*t7*t10*t14*t34+8*t2*t42*t7*t10*t14*t46&
&   -8*t2*t3*t7*t10*t14*t24+8*t2*t3*t7*t10*t55*t44*t33+&
&   6*t2*t42*t7*t10*t55/t127*t45+24*t2*t42/t138*t140+&
&   16*t54*t42*t120*t140+16*t2*t3*t120*t122+16*t2&
&   *t42*t7*t72*t78-8*t2*t3*t7*t10*t55*t44*t23-8*t54*t3*t7*t175&
&   +2*t2*t7*t10*t55*t17*t100+2*t2*t7*t10*t55*t17*t93+&
&   8*t54*t42*t7*t79+16*t2*t42*t7*t72*t56+6*t2*t42*t7*t10/t85&
&   *t64+24*t2*t42*t7/t62*t65+8*t54*t42*t7*t57-&
&   16*t2*t3*t7*t72*t174+8*t54*t3*t7*t122-16*t2*t3*t120*t175&
&   +16*t2*t42*t120*t79+16*t2*t42*t120*t57+16*t54*t42*t7*t152+&
&   32*t2*t42*t120*t152+16*t2*t3*t7*t72*t121-12*t2*t157+&
&   6/t160*t42*t157
 end if

end subroutine gpp1cc
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp1cc
!! NAME
!! psp1cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 1 and 5 of pseudopotentials.
!! WARNING : the fifth derivate is actually set to zero
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fchrg=magnitude of the core charge correction
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! NOTES
!! This is a revised expression for core density (5 Nov 1992) :
!! density(r)=fchrg*gg(xx)
!! with
!! $ gg(xx)=(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4 xx^2)(1-xx^2)})^2 $
!! and
!! $ xx=\frac{r}{rchrg}=\frac{r}{xcccrc/3.0d0}=3*\frac{r}{xcccrc}=3*yy $
!!
!! Code for gg(xx), gp(xx), and gpp(xx) has been tested by numerical
!! derivatives--looks ok. gpp(x) should still be rewritten.
!! The argument of xccc1d is assumed to be normalized, and to vary
!! from yy=0 to 1 (from r=0 to r=xcccrc, or from xx=0 to 3)
!! Thus :
!!{{\ \begin{equation}
!! xccc1d(yy)=fchrg*[\frac{\sin(2*\pi*(3yy))}
!! {(6*\pi*(3yy))(1-4*(3yy)^2)(1-(3yy)^2)}]^2
!!\end{equation} }}
!!
!! WARNINGS
!! Warning : the fifth derivative is not yet delivered.
!!
!! PARENTS
!!      psp1in,psp5in
!!
!! CHILDREN
!!      gg1cc,gp1cc,gpp1cc,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp1cc(fchrg,n1xccc,xccc1d)

 use defs_basis
 use m_splines
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp1cc'
 use interfaces_64_psp, except_this_one => psp1cc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1xccc
 real(dp),intent(in) :: fchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1xccc,ider
 real(dp) :: der1,dern,factor,gg1cc_xx,gp1cc_xx,gpp1cc_xx,xx
 character(len=500) :: message
!arrays
 real(dp),allocatable :: ff(:),ff2(:),work(:),yy(:)

! *************************************************************************

 ABI_ALLOCATE(ff,(n1xccc))
 ABI_ALLOCATE(ff2,(n1xccc))
 ABI_ALLOCATE(work,(n1xccc))
 ABI_ALLOCATE(yy,(n1xccc))

 if(n1xccc > 1)then
   factor=one/dble(n1xccc-1)
   do i1xccc=1,n1xccc
     yy(i1xccc)=(i1xccc-1)*factor
   end do
 else
   write(message, '(a,i0)' )' n1xccc should larger than 1, while it is n1xccc=',n1xccc
   MSG_BUG(message)
 end if

!Initialization, to avoid some problem with some compilers
 xccc1d(1,:)=zero ; xccc1d(n1xccc,:)=zero

!Take care of each derivative separately
 do ider=0,2

   if(ider==0)then
!    Generate spline fitting for the function gg
     do i1xccc=1,n1xccc
       xx=three*yy(i1xccc)
       call gg1cc(gg1cc_xx,xx)
       ff(i1xccc)=fchrg*gg1cc_xx
     end do
!    Complete with derivatives at end points
     der1=zero
     call gp1cc(gp1cc_xx,three)
     dern=three*fchrg*gp1cc_xx
   else if(ider==1)then
!    Generate spline fitting for the function gp
     do i1xccc=1,n1xccc
       xx=three*yy(i1xccc)
       call gp1cc(gp1cc_xx,xx)
       ff(i1xccc)=three*fchrg*gp1cc_xx
     end do
!    Complete with derivatives at end points, already estimated
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   else if(ider==2)then
!    Generate spline fitting for the function gpp
!    (note : the function gpp has already been estimated, for the spline
!    fitting of the function gg, but it is replaced here by the more
!    accurate analytic derivative)
     do i1xccc=1,n1xccc
       xx=three*yy(i1xccc)
       call gpp1cc(gpp1cc_xx,xx)
       ff(i1xccc)=9.0_dp*fchrg*gpp1cc_xx
     end do
!    Complete with derivatives of end points
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   end if

!  Produce second derivative numerically, for use with splines
   call spline(yy,ff,n1xccc,der1,dern,ff2)
   xccc1d(:,ider+1)=ff(:)
   xccc1d(:,ider+3)=ff2(:)

 end do

 xccc1d(:,6)=zero

!DEBUG
!write(std_out,*)' psp1cc : output of core charge density and derivatives '
!write(std_out,*)'   yy          gg           gp  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) yy(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   yy          gpp          gg2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) yy(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   yy          gp2          gpp2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) yy(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(yy)

end subroutine psp1cc
!!***

!!****f* ABINIT/gg1cc
!! NAME
!! gg1cc
!!
!! FUNCTION
!! gg1cc_xx=$(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4xx^2)(1-xx^2)})^2$
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

!!****f* ABINIT/gp1cc
!! NAME
!! gp1cc
!!
!! FUNCTION
!! Derivative of gg(xx) wrt xx.
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
!! PARENTS
!!      psp1cc
!!
!! CHILDREN
!!
!! SOURCE

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

!!****f* ABINIT/gpp1cc
!! NAME
!! gpp1cc
!!
!! FUNCTION
!! Second derivative of gg wrt xx.
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

!!****m* ABINIT/m_psptk
!! NAME
!!  m_psptk
!!
!! FUNCTION
!!  This module collects low-level procedures used by the other psp modules
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (XG, DCA, MM, DRH, FrD, GZ, AF)
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

module m_psptk

 use defs_basis
 use m_errors
 use m_abicore
 use m_splines

 use m_numeric_tools,  only : ctrap
 use m_special_funcs,  only : sbf8

 implicit none

 private
!!***

 public :: psp1cc
 public :: psp5lo
 public :: psp5nl
 public :: psp8lo
 public :: psp8nl
 public :: cc_derivatives
!!***

contains
!!***

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
!! Warning: the fifth derivative is not yet delivered.
!!
!! PARENTS
!!      m_psp1,m_psp5
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine psp1cc(fchrg,n1xccc,xccc1d)

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
!!      m_psptk
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine gg1cc(gg1cc_xx,xx)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: gg1cc_xx

!Local variables -------------------------------------------
!The c s are coefficients for Taylor expansion of the analytic form near xx=0, 1/2, and 1.
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
!!      m_psptk
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine gp1cc(gp1cc_xx,xx)

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
!!      m_psptk
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine gpp1cc(gpp1cc_xx,xx)

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

!!****f* ABINIT/psp5lo
!! NAME
!! psp5lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on logarithmic grid using related uniform
!! grid in exponent and corrected trapezoidal integration.
!!
!! INPUTS
!!  al=spacing in exponent for radial atomic grid.
!!  mmax=number of radial r grid points (logarithmic atomic grid).
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rad(mmax)=r grid values (bohr).
!!  vloc(mmax)=V(r) on radial grid.
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm=$ 4\pi\int[r^2 (V(r)+\frac{Zv}{r}dr]$.
!!{{\\ \begin{equation}
!!  q2vq(mqgrid)
!!   =q^2 V(q)
!!   = -\frac{Zv}{\pi}
!!     + q^2 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 V(r)+r Zv)dr].
!!\end{equation} }}
!!  yp1,ypn=derivative of q^2 V(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      m_psp5,m_psp6
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine psp5lo(al,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&
&                  vloc,yp1,ypn,zion)

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mmax,mqgrid
 real(dp),intent(in) :: al,zion
 real(dp),intent(out) :: epsatm,yp1,ypn
!arrays
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables-------------------------------
!scalars
 integer :: iq,ir
 real(dp),parameter :: scale=10.0d0
 real(dp) :: arg,result,rmtoin,test,ztor1
!arrays
 real(dp),allocatable :: work(:)

! *************************************************************************

 ABI_ALLOCATE(work,(mmax))

!Do q=0 separately (compute epsatm)
!Do integral from 0 to r1
 ztor1=(zion/2.0d0+rad(1)*vloc(1)/3.d0)*rad(1)**2

!Set up integrand for q=0: $ \int[r^2 (V(r)+\frac{Zv}{r}) dr]$
!with extra factor of r to convert to uniform grid in exponent
 do ir=1,mmax
!  First handle tail region
   test=vloc(ir)+zion/rad(ir)
!  DEBUG
!  write(std_out,*)ir,rad(ir),test
!  ENDDEBUG
!  Ignore small contributions, or impose a cut-off in the case
!  the pseudopotential data are in single precision.
!  (it is indeed expected that vloc is very close to zero beyond 20,
!  so a value larger than 2.0d-8 is considered anomalous)
   if (abs(test)<1.0d-20 .or. (rad(ir)>20.0d0 .and. abs(test)>2.0d-8) ) then
     work(ir)=zero
   else
     work(ir)=(rad(ir)*rad(ir))*(rad(ir)*vloc(ir)+zion)
   end if
 end do
!DEBUG
!write(std_out,*)' psp5lo : stop '
!stop
!ENDDEBUG

!Do integral from r(1) to r(max)
 call ctrap(mmax,work,al,result)
!Do integral from r(mmax) to infinity
!compute decay length lambda at r(mmax)
!$\lambda=-\log((rad(im1)*vloc(im1)+zion)$/ &
!$(rad(imat)*vloc(imat)+zion))/(rad(im1)-rad(imat))$
!rmtoin=$(rad(mmax)*vloc(mmax)+zion)*(rad(mmax)+1.d0/\lambda)/\lambda$
!Due to inability to fit exponential decay to r*V(r)+Zv
!in tail, NO TAIL CORRECTION IS APPLIED
!(numerical trouble might be removed if atomic code is
!cleaned up in tail region)
 rmtoin=0.0d0

 epsatm=4.d0*pi*(result+ztor1+rmtoin)

 q2vq(1)=-zion/pi

!Loop over q values
 do iq=2,mqgrid
   arg=2.d0*pi*qgrid(iq)
!  ztor1=$ -Zv/\pi+2q \int_0^{r1}[\sin(2\pi q r)(rV(r)+Zv) dr]$
   ztor1=(vloc(1)*sin(arg*rad(1))/arg-(rad(1)*vloc(1)+zion)* &
&   cos(arg*rad(1)) )/pi

!  set up integrand
   do  ir=1,mmax
     test=vloc(ir)+zion/rad(ir)
!    Ignore contributions within decade of machine precision
     if ((scale+abs(test)).eq.scale) then
       work(ir)=zero
     else
       work(ir)=rad(ir)*sin(arg*rad(ir))*(rad(ir)*vloc(ir)+zion)
     end if
   end do
!  do integral from r(1) to r(mmax)
   call ctrap(mmax,work,al,result)

!  do integral from r(mmax) to infinity
!  rmtoin=(r(mmax)*vr(mmax)+zion)*(lambda*sin(arg*r(mmax))+
!  arg*cos(arg*r(mmax)))/(arg**2+lambda**2)
!  See comment above; no tail correction
   rmtoin=0.0d0

!  store q^2 v(q)
   q2vq(iq)=ztor1+2.d0*qgrid(iq)*(result+rmtoin)

 end do

!Compute derivatives of q^2 v(q) at ends of interval
 yp1=0.0d0
!ypn=$ 2\int_0^\infty[(\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r)(r V(r)+Z) dr]$
!integral from 0 to r1
 arg=2.0d0*pi*qgrid(mqgrid)
 ztor1=zion*rad(1)*sin(arg*rad(1))
 ztor1=ztor1+ 3.d0*rad(1)*vloc(1)*cos(arg*rad(1))/arg + &
& (rad(1)**2-1.0d0/arg**2)*vloc(1)*sin(arg*rad(1))
!integral from r(mmax) to infinity is overkill; ignore
!set up integrand
 do ir=1,mmax
   test=vloc(ir)+zion/rad(ir)
!  Ignore contributions within decade of machine precision
   if ((scale+abs(test)).eq.scale) then
     work(ir)=0.0d0
   else
     work(ir)=rad(ir)*(sin(arg*rad(ir))+arg*rad(ir)*cos(arg*rad(ir))) * &
&     (rad(ir)*vloc(ir)+zion)
   end if
 end do
 call ctrap(mmax,work,al,result)
 ypn=2.0d0 * (ztor1 + result)

 ABI_DEALLOCATE(work)

end subroutine psp5lo
!!***

!!****f* ABINIT/psp5nl
!! NAME
!! psp5nl
!!
!! FUNCTION
!! Make Kleinman-Bylander form factors f_l(q) for each l from
!! 0 to lmax; Vloc is assumed local potential.
!!
!! INPUTS
!!  al=grid spacing in exponent for radial grid
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!   Usually lmax=1, sometimes = 0 (e.g. for oxygen); lmax <= 2 allowed.
!!  mmax=number of radial grid points for atomic grid
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for q grid
!!  qgrid(mqgrid)=values at which form factors are returned
!!  rad(mmax)=radial grid values
!!  vloc(mmax)=local pseudopotential on radial grid
!!  vpspll(mmax,3)=nonlocal pseudopotentials for each l on radial grid
!!  wfll(mmax,3)=reference state wavefunctions on radial grid
!!                mmax and mqgrid
!!
!! OUTPUT
!!  ekb(mpsang)=Kleinman-Bylander energy,
!!             {{\\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!               {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!               \end{equation} }}
!!                for each l
!!  ffspl(mqgrid,2,mpsang)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!
!! NOTES
!! u_l(r) is reference state wavefunction (input as wf);
!! j_l(q) is a spherical Bessel function;
!! dV_l(r) = vpsp_l(r)-vloc(r) for angular momentum l;
!! f_l(q) = $ \int_0^{rmax}[j_l(2\pi q r) u_l(r) dV_l(r) r dr]/\sqrt{dvms}$
!! where dvms = $\int_0^{rmax} [(u_l(r) dV_l(r))^2 dr]$ is the mean
!! square value of the nonlocal correction for angular momentum l.
!! Xavier Gonze s E_KB = $ dvms/\int_0^{rmax}[(u_l(r))^2 dV_l(r) dr]$.
!! This is the eigenvalue of the Kleinman-Bylander operator and sets
!! the energy scale of the nonlocal psp corrections.
!!
!! PARENTS
!!      m_psp5,m_psp6
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine psp5nl(al,ekb,ffspl,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,mmax,mpsang,mqgrid
 real(dp),intent(in) :: al
!arrays
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vloc(mmax),vpspll(mmax,mpsang)
 real(dp),intent(in) :: wfll(mmax,mpsang)
 real(dp),intent(out) :: ekb(mpsang),ffspl(mqgrid,2,mpsang)

!Local variables-------------------------------
!scalars
 integer,parameter :: dpsang=5
 integer :: iq,ir,lp1
 real(dp) :: arg,bessel,dvwf,qr,result,yp1,ypn,ztor1
 character(len=500) :: message
!arrays
 real(dp) :: ckb(dpsang),dvms(dpsang),eta(dpsang),renorm(dpsang)
 real(dp),allocatable :: work1(:),work2(:),work3(:),work4(:)

!*************************************************************************

!l=0,1,2 and 3 spherical Bessel functions
!The accuracy of the bes1, bes2, bes3 functions for small arguments
!may be insufficient. In the present version
!of the routines, some care is taken with the value of the argument.
!If smaller than 1.d-3, a two terms
!Taylor series expansion is prefered.
! bes0(arg)=sin(arg)/arg
! bes1(arg)=(sin(arg)-arg*cos(arg))/arg**2
! bes2(arg)=( (3.0d0-arg**2)*sin(arg)-&
!& 3.0d0*arg*cos(arg) )      /arg**3

! bes3(arg)=(15.d0*sin(arg)-15.d0*arg*cos(arg) &
!& -6.d0*arg**2*sin(arg)+arg**3*cos(arg) )/arg**4

!Zero out Kleinman-Bylander energies ekb
 ekb(:)=0.0d0

 ABI_ALLOCATE(work1,(mmax))
 ABI_ALLOCATE(work2,(mmax))
 ABI_ALLOCATE(work3,(mmax))
 ABI_ALLOCATE(work4,(mmax))

!Allow for no nonlocal correction (lmax=-1)
 if (lmax/=-1) then

!  Check that lmax is within allowed range
   if (lmax<0.or.lmax>3) then
     write(message, '(a,i12,a,a,a,a,a,a,a)' )&
&     'lmax=',lmax,' is not an allowed value.',ch10,&
&     'Allowed values are -1 for no nonlocal correction or else',ch10,&
&     '0, 1,2 or 3 for maximum l nonlocal correction.',ch10,&
&     'Action: check the input atomic psp data file for lmax.'
     MSG_ERROR(message)
   end if

!  Compute normalizing integrals eta=<dV> and mean square
!  nonlocal psp correction dvms=<dV^2>
!  "dvwf" consistently refers to dV(r)*wf(r) where dV=nonlocal correction
   do lp1=1,lmax+1

!    integral from 0 to r1
     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)
     ztor1=(wfll(1,lp1)*dvwf)*rad(1)/dble(2*(lp1-1)+3)
!    integrand for r1 to r(mmax) (incl extra factor of r)
     do ir=1,mmax
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)
       work1(ir)=rad(ir)*(wfll(ir,lp1)*dvwf)
     end do
!    do integral by corrected trapezoidal integration
     call ctrap(mmax,work1,al,result)
     eta(lp1)=ztor1+result

     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)
     ztor1=dvwf**2*rad(1)/dble(2*(lp1-1)+3)
     do ir=1,mmax
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)
       work1(ir)=rad(ir)*(dvwf**2)
     end do
     call ctrap(mmax,work1,al,result)
     dvms(lp1)=ztor1+result

!    DEBUG
!    Compute the norm of wfll
!    wf=wfll(1,lp1)
!    ztor1=wf**2*rad(1)/dble(2*(lp1-1)+3)
!    do ir=1,mmax
!    wf=wfll(ir,lp1)
!    work1(ir)=rad(ir)*(wf**2)
!    end do
!    call ctrap(mmax,work1,al,result)
!    norm=ztor1+result
!    write(std_out,*)' lp1, norm',lp1,norm
!    ENDDEBUG

!    If dvms is not 0 for any given angular momentum l,
!    compute Xavier Gonze s definition of the Kleinman-Bylander
!    energy E_KB = dvms/eta.  In this case also renormalize
!    the projection operator to u_KB(r)=$u_l(r)*dV(r)/\sqrt{dvms}$.
!    This means dvwf gets multiplied by the normalization factor
!    "renorm"=$1/\sqrt{dvms}$ as seen below.
     if (dvms(lp1)/=0.0d0) then
       ekb(lp1)=dvms(lp1)/eta(lp1)
       renorm(lp1)=1.0d0/sqrt(dvms(lp1))
!      ckb is Kleinman-Bylander "cosine" (Xavier Gonze)
       ckb(lp1)=eta(lp1)/sqrt(dvms(lp1))
     else
       ekb(lp1)=0.0d0
     end if

   end do

!  l=0 form factor if ekb(1) not 0 (lmax always at least 0)
   if (ekb(1)/=0.0d0) then

!    do q=0 separately
     lp1=1
!    0 to r1 integral
     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
     ztor1=(rad(1)*dvwf)*rad(1)/3.0d0
!    integrand
     do ir=1,mmax
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
       work1(ir)=rad(ir)*(rad(ir)*dvwf)
     end do
     call ctrap(mmax,work1,al,result)
     ffspl(1,1,1)=ztor1+result

!    do rest of q points
     do iq=2,mqgrid
       arg=two_pi*qgrid(iq)
!      0 to r1 integral
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       ztor1=(bes0_psp5(arg*rad(1))*rad(1)*dvwf)*rad(1)/3.0d0
!      integrand
       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
         work1(ir)=rad(ir)*(rad(ir)*bes0_psp5(arg*rad(ir))*dvwf)
       end do
       call ctrap(mmax,work1,al,result)
       ffspl(iq,1,1)=ztor1+result
     end do

!    Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
!    yp1=0 for l=0
     yp1=0.0d0
!    ypn=$ \int [2\pi r (-bes1(2\pi r q)) wf(r) dV(r) r dr]$
     arg=two_pi*qgrid(mqgrid)
     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
     qr=arg*rad(1)
     if(qr<1.d-3)then
       bessel=(10.d0-qr*qr)*qr/30.0d0
     else
       bessel=bes1_psp5(qr)
     end if
!    ztor1=(-bes1(arg*rad(1))*two_pi*rad(1)*r(1)*dvwf)*rad(1)/5.0d0
     ztor1=(-bessel*two_pi*rad(1)*rad(1)*dvwf)*rad(1)/5.0d0
     do ir=1,mmax
       qr=arg*rad(ir)
       if(qr<1.d-3)then
         bessel=(10.d0-qr*qr)*qr/30.0d0
       else
         bessel=bes1_psp5(qr)
       end if
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!      work(ir)=rad(ir)*(-bes1(arg*rad(ir))*two_pi*rad(ir)*rad(ir)*dvwf)
       work1(ir)=rad(ir)*(-bessel*two_pi*rad(ir)*rad(ir)*dvwf)
     end do
     call ctrap(mmax,work1,al,result)
     ypn=ztor1+result

!    Fit spline to get second derivatives by spline fit
     call spline(qgrid,ffspl(1,1,1),mqgrid,yp1,ypn,ffspl(1,2,1))

   else
!    or else put nonlocal correction at l=0 to 0
     ffspl(:,:,1)=0.0d0
   end if

!  Finished if lmax=0 (highest nonlocal correction)
!  Do l=1 form factor if ekb(2) not 0 and lmax>=1
   if (lmax>0)then
     if(ekb(2)/=0.0d0) then

       lp1=2
!      do q=0 separately: f_1(q=0) vanishes !
       ffspl(1,1,2)=0.0d0

!      do rest of q points
       do iq=2,mqgrid
         arg=two_pi*qgrid(iq)
         dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
         qr=arg*rad(1)
         if(qr<1.d-3)then
           bessel=(10.d0-qr*qr)*qr/30.0d0
         else
           bessel=bes1_psp5(qr)
         end if
!        ztor1=(bes1(arg*rad(1))*rad(1)*dvwf)*rad(1)/5.0d0
         ztor1=(bessel*rad(1)*dvwf)*rad(1)/5.0d0

         do ir=1,mmax
           qr=arg*rad(ir)
           if(qr<1.d-3)then
             bessel=(10.d0-qr*qr)*qr/30.0d0
           else
             bessel=bes1_psp5(qr)
           end if
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
           work2(ir)=rad(ir)*(rad(ir)*bessel*dvwf)
         end do

         call ctrap(mmax,work2,al,result)
         ffspl(iq,1,2)=ztor1+result
       end do

!      Compute yp1,ypn for l=1
!      yp1=$\displaystyle \int [2\pi r^2 wf(r) dV(r)]/3$
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       ztor1=((two_pi*rad(1)**2)*dvwf)*rad(1)/(3.0d0*5.0d0)
       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
         work2(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf/3.0d0)
       end do
       call ctrap(mmax,work2,al,result)
       yp1=ztor1+result
!      ypn=$\int [2\pi r^2 wf(r) dV(r) (j_0(x)-(2/x)j_1(x)) dr]$
!      where x=2 Pi qgrid(mqgrid) r
       arg=two_pi*qgrid(mqgrid)
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       qr=arg*rad(1)
       if(qr<1.d-3)then
         bessel=(10.d0-3.0d0*qr*qr)/30.0d0
       else
         bessel=bes0_psp5(qr)-2.d0*bes1_psp5(qr)/qr
       end if
!      ztor1=( (two_pi*rad(1)**2)*dvwf* (bes0(arg*rad(1))-
!      2.0d0*bes1(arg*rad(1))/(arg*rad(1))) ) * rad(1)/5.0d0
       ztor1=( (two_pi*rad(1)**2)*dvwf*bessel)*  rad(1)/5.0d0

       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
         qr=arg*rad(ir)
         if(qr<1.d-3)then
           bessel=(10.d0-3.0d0*qr*qr)/30.0d0
         else
           bessel=bes0_psp5(qr)-2.d0*bes1_psp5(qr)/qr
         end if
!        work(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*
!        (bes0(arg*rad(ir))-2.d0*bes1(arg*rad(ir))/(arg*rad(ir))) )
         work2(ir)=rad(ir)*(two_pi*rad(ir)**2)*dvwf*bessel
       end do
       call ctrap(mmax,work2,al,result)
       ypn=ztor1+result

!      Fit spline for l=1 Kleinman-Bylander form factor
       call spline(qgrid,ffspl(1,1,2),mqgrid,yp1,ypn,ffspl(1,2,2))

     else
!      or else put form factor to 0 for l=1
       ffspl(:,:,2)=0.0d0
     end if
!    Endif condition of lmax>0
   end if

!  Finished if lmax=1 (highest nonlocal correction)
!  Do l=2 nonlocal form factor if eta(3) not 0 and lmax>=2
   if (lmax>1)then
     if(ekb(3)/=0.0d0) then

       lp1=3
!      do q=0 separately; f_2(q=0) vanishes
       ffspl(1,1,3)=0.0d0

!      do rest of q points
       do iq=2,mqgrid
         arg=two_pi*qgrid(iq)
         dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
         qr=arg*rad(1)
         if(qr<1.d-3)then
           bessel=qr*qr/15.0d0-qr**4/210.0d0
         else
           bessel=bes2_psp5(qr)
         end if
!        ztor1=(bes2(arg*rad(1))*rad(1)*dvwf)*rad(1)/7.0d0
         ztor1=(bessel*rad(1)*dvwf)*rad(1)/7.0d0
         do ir=1,mmax
           qr=arg*rad(ir)
           if(qr<1.d-3)then
             bessel=qr*qr/15.0d0-qr**4/210.0d0
           else
             bessel=bes2_psp5(qr)
           end if
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!          work(ir)=rad(ir)*(r(ir)*bes2(arg*rad(ir))*dvwf)
           work3(ir)=rad(ir)*(rad(ir)*bessel*dvwf)
         end do
         call ctrap(mmax,work3,al,result)
         ffspl(iq,1,3)=ztor1+result
       end do

!      Compute yp1,ypn for l=2
!      yp1=0 for l=2
       yp1=0.0d0
!      ypn=$\int [2 \pi r^2 wf(r) dV(r) (j_1(x)-(3/x)j_2(x)) dr]$
!      where x=2 Pi qgrid(mqgrid) r
       arg=two_pi*qgrid(mqgrid)
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       qr=arg*rad(1)
       if(qr<1.d-3)then
         bessel=qr*2.0d0/15.0d0-qr**3*4.0d0/210.0d0
       else
         bessel=bes1_psp5(qr)-3.0d0*bes2_psp5(qr)/qr
       end if
!      ztor1=( (two_pi*rad(1)**2)*dvwf* (bes1(arg*rad(1))-
!      3.0d0*bes2(arg*rad(1))/(arg*rad(1))) ) * rad(1)/7.0d0
       ztor1=( (two_pi*rad(1)**2)*dvwf* bessel ) * rad(1)/7.0d0
       do ir=1,mmax
         qr=arg*rad(ir)
         if(qr<1.d-3)then
           bessel=qr*2.0d0/15.0d0-qr**3*4.0d0/210.0d0
         else
           bessel=bes1_psp5(qr)-3.0d0*bes2_psp5(qr)/qr
         end if
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!        work3(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*
!        (bes1(arg*rad(ir))-3.d0*bes2(arg*rad(ir))/(arg*rad(ir))) )
         work3(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*bessel)
       end do
       call ctrap(mmax,work3,al,result)
       ypn=ztor1+result

!      Fit spline for l=2 Kleinman-Bylander form factor
       call spline(qgrid,ffspl(1,1,3),mqgrid,yp1,ypn,ffspl(1,2,3))

     else
!      or else put form factor to 0 for l=1
       ffspl(:,:,3)=0.0d0
     end if
!    Endif condition of lmax>1
   end if

!  Finished if lmax=2 (highest nonlocal correction)
!  Do l=3 nonlocal form factor if eta(4) not 0 and lmax>=3
   if (lmax>2)then
     if(ekb(4)/=0.0d0) then

       lp1=4
!      do q=0 separately; f_3(q=0) vanishes
       ffspl(1,1,4)=0.0d0

!      do rest of q points
       do iq=2,mqgrid
         arg=two_pi*qgrid(iq)
         dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
         qr=arg*rad(1)
         if(qr<1.d-3)then
           bessel=qr*qr*qr/105.0d0-qr**5/1890.0d0+qr**7/83160.0d0
         else
           bessel=bes3_psp5(qr)
         end if
!        ztor1=(bes3(arg*rad(1))*rad(1)*dvwf)*rad(1)/9.0d0
         ztor1=(bessel*rad(1)*dvwf)*rad(1)/9.0d0
         do ir=1,mmax
           qr=arg*rad(ir)
           if(qr<1.d-3)then
             bessel=qr*qr*qr/105.0d0-qr**5/1890.0d0+qr**7/83160.0d0
           else
             bessel=bes3_psp5(qr)
           end if
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!          work(ir)=rad(ir)*(rad(ir)*bes3(arg*rad(ir))*dvwf)
           work4(ir)=rad(ir)*(rad(ir)*bessel*dvwf)
         end do
         call ctrap(mmax,work4,al,result)
         ffspl(iq,1,4)=ztor1+result
       end do

!      Compute yp1,ypn for l=3
!      yp1=0 for l=3
       yp1=0.0d0
!      ypn=$\int [2\pi r^2 wf(r) dV(r) (j_2(x)-(4/x)j_3(x)) dr]$
!      where x=2 Pi qgrid(mqgrid) r
       arg=two_pi*qgrid(mqgrid)
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       qr=arg*rad(1)
       if(qr<1.d-3)then
         bessel=3.d0*qr**2/105.0d0-5.d0*qr**4/1890.0d0+7.d0*qr**6/83160.0d0
       else
         bessel=bes2_psp5(qr)-4.0d0*bes3_psp5(qr)/qr
       end if
!      ztor1=( (two_pi*rad(1)**2)*dvwf* (bes2(arg*rad(1))-
!      3.0d0*bes3(arg*rad(1))/(arg*rad(1))) ) * rad(1)/9.0d0
       ztor1=( (two_pi*rad(1)**2)*dvwf* bessel ) * rad(1)/9.0d0
       do ir=1,mmax
         qr=arg*rad(ir)
         if(qr<1.d-3)then
           bessel=3.d0*qr**2/105.0d0-5.d0*qr**4/1890.0d0+7.d0*qr**6/83160.0d0
         else
           bessel=bes2_psp5(qr)-4.0d0*bes3_psp5(qr)/qr
         end if
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!        work4(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*
!        (bes2(arg*rad(ir))-4.d0*bes3(arg*rad(ir))/(arg*rad(ir))) )
         work4(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*bessel)
       end do
       call ctrap(mmax,work4,al,result)
       ypn=ztor1+result

!      Fit spline for l=3 Kleinman-Bylander form factor
       call spline(qgrid,ffspl(1,1,4),mqgrid,yp1,ypn,ffspl(1,2,4))

     else
!      or else put form factor to 0 for l=3
       ffspl(:,:,4)=0.0d0
     end if
!    Endif condition of lmax>2
   end if

!  Endif condition lmax/=-1
 end if

!DEBUG
!write(std_out,*) 'EKB=',(ekb(iq),iq=1,3)
!write(std_out,*) 'COSKB=',(ckb(iq),iq=1,3)
!ENDDEBUG

 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)
 ABI_DEALLOCATE(work3)
 ABI_DEALLOCATE(work4)

 contains

   function  bes0_psp5(arg)

   real(dp) :: bes0_psp5,arg
   bes0_psp5=sin(arg)/arg
 end function bes0_psp5

   function bes1_psp5(arg)

   real(dp) :: bes1_psp5,arg
   bes1_psp5=(sin(arg)-arg*cos(arg))/arg**2
 end function bes1_psp5

   function bes2_psp5(arg)

   real(dp) :: bes2_psp5,arg
   bes2_psp5=( (3.0d0-arg**2)*sin(arg)-&
&   3.0d0*arg*cos(arg) )      /arg**3
 end function bes2_psp5

   function bes3_psp5(arg)

   real(dp) :: bes3_psp5, arg
   bes3_psp5=(15.d0*sin(arg)-15.d0*arg*cos(arg) &
&   -6.d0*arg**2*sin(arg)+arg**3*cos(arg) )/arg**4
 end function bes3_psp5

end subroutine psp5nl
!!***

!!****f* ABINIT/psp8lo
!! NAME
!! psp8lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on linear grid interpolated from the linear input
!! grid with a spacing adjusted to ensure convergence at the maximum
!! wavevector using corrected trapezoidal integration.
!!
!! INPUTS
!!  amesh=spacing for linear radial atomic grid.
!!  mmax=number of radial r grid points
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rad(mmax)=r grid values (bohr).
!!  vloc(mmax)=V(r) on radial grid.
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm=$ 4\pi\int[r^2 (V(r)+\frac{Zv}{r}dr]$.
!!{{\\ \begin{equation}
!!  q2vq(mqgrid)
!!   =q^2 V(q)
!!   = -\frac{Zv}{\pi}
!!     + q^2 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 V(r)+r Zv)dr].
!!\end{equation} }}
!!  yp1,ypn=derivative of q^2 V(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      m_psp8,m_psp9
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine psp8lo(amesh,epsatm,mmax,mqgrid,qgrid,q2vq,rad,vloc,yp1,ypn,zion)

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mmax,mqgrid
 real(dp),intent(in) :: amesh,zion
 real(dp),intent(out) :: epsatm,yp1,ypn
!arrays
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables-------------------------------
!Following parameter controls accuracy of Fourier transform based on qmax
!and represents the minimun number of integration points in one period.
!scalars
 integer,parameter :: NPT_IN_2PI=200
 integer :: ider,iq,ir,irmu,irn,mesh_mult,mmax_new
 real(dp) :: amesh_new,arg,fp1,fpn,qmesh,result,ztor1
!arrays
 real(dp),allocatable :: rad_new(:),rvlpz(:),rvlpz_new(:),sprvlpz(:,:),work(:)

! *************************************************************************

 ABI_ALLOCATE(work,(mmax))
 ABI_ALLOCATE(rvlpz,(mmax))

!Do q=0 separately (compute epsatm)
 ztor1=(zion/2.0d0+rad(1)*vloc(1)/3.d0)*rad(1)**2
!Set up integrand for q=0: $ \int[r^2 (V(r)+\frac{Zv}{r}) dr]$
 do ir=1,mmax
   rvlpz(ir)=rad(ir)*vloc(ir)+zion
   work(ir)=rad(ir)*rvlpz(ir)
 end do

!Do integral from zero to r(max)
 call ctrap(mmax,work,amesh,result)

 epsatm=4.d0*pi*result
 q2vq(1)=-zion/pi

!Find r mesh spacing necessary for accurate integration at qmax
 amesh_new=2.d0*pi/(NPT_IN_2PI*qgrid(mqgrid))

!Choose submultiple of input mesh
 mesh_mult=int(amesh/amesh_new) + 1
 mmax_new=mesh_mult*(mmax-1)+1
 amesh_new=amesh/dble(mesh_mult)

 ABI_ALLOCATE(rad_new,(mmax_new))
 ABI_ALLOCATE(rvlpz_new,(mmax_new))

 if(mesh_mult==1) then
   rad_new(:)=rad(:)
   rvlpz_new(:)=rvlpz(:)
 else
!  Set up spline and interpolate to finer mesh.
!  First, compute derivatives at end points
   fp1=(-50.d0*rvlpz(1)+96.d0*rvlpz(2)-72.d0*rvlpz(3)+32.d0*rvlpz(4)&
&   -6.d0*rvlpz(5))/(24.d0*amesh)
   fpn=(6.d0*rvlpz(mmax-4)-32.d0*rvlpz(mmax-3)+72.d0*rvlpz(mmax-2)&
&   -96.d0*rvlpz(mmax-1)+50.d0*rvlpz(mmax))/(24.d0*amesh)
   ABI_ALLOCATE(sprvlpz,(mmax,2))
   work(:)=zero

!  Spline fit
   call spline(rad, rvlpz,mmax,fp1,fpn,sprvlpz(:,2))
   sprvlpz(:,1)=rvlpz(:)

!  Set up new radial mesh
   irn=1
   do ir=1,mmax-1
     do irmu=0,mesh_mult-1
       rad_new(irn)=rad(ir)+dble(irmu)*amesh_new
       irn=irn+1
     end do
   end do
   rad_new(mmax_new)=rad(mmax)

   ider=0
   call splfit(rad,work,sprvlpz,ider,rad_new,rvlpz_new,mmax,mmax_new)

   ABI_DEALLOCATE(sprvlpz)
   ABI_DEALLOCATE(work)
   ABI_ALLOCATE(work,(mmax_new))
 end if

!Loop over q values
 do iq=2,mqgrid
   arg=2.d0*pi*qgrid(iq)

!  Set up integrand
   do  ir=1,mmax_new
     work(ir)=sin(arg*rad_new(ir))*rvlpz_new(ir)
   end do

!  Do integral from zero to rad(mmax)
   call ctrap(mmax_new,work,amesh_new,result)

!  Store q^2 v(q)
   q2vq(iq)=q2vq(1)+2.d0*qgrid(iq)*result

 end do

!Compute derivatives of q^2 v(q) at ends of interval
 qmesh=qgrid(2)-qgrid(1)
 yp1=(-50.d0*q2vq(1)+96.d0*q2vq(2)-72.d0*q2vq(3)+32.d0*q2vq(4)&
& -6.d0*q2vq(5))/(24.d0*qmesh)
 ypn=(6.d0*q2vq(mqgrid-4)-32.d0*q2vq(mqgrid-3)+72.d0*q2vq(mqgrid-2)&
& -96.d0*q2vq(mqgrid-1)+50.d0*q2vq(mqgrid))/(24.d0*qmesh)

 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(rad_new)
 ABI_DEALLOCATE(rvlpz_new)
 ABI_DEALLOCATE(rvlpz)

end subroutine psp8lo
!!***

!!****f* ABINIT/psp8nl
!! NAME
!! psp8nl
!!
!! FUNCTION
!! Make Kleinman-Bylander/Bloechl form factors f_ln(q) for each
!!  projector n for each angular momentum l excepting an l corresponding
!!  to the local potential.
!! Note that an arbitrary local potential can be used, so all l from
!!  0 to lmax may be represented.
!!
!! INPUTS
!!  amesh=grid spacing for uniform (linear) radial grid
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!    lmax <= 2 allowed.
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=number of radial grid points for atomic grid
!!  mqgrid=number of grid points for q grid
!!  pspso=spin-orbit characteristics, govern the content of ffspl and ekb
!!   if =0 : this input requires NO spin-orbit characteristics of the psp
!!   if =2 : this input requires HGH or psp8 characteristics of the psp
!!   if =3 : this input requires HFN characteristics of the psp
!!  qgrid(mqgrid)=values at which form factors are returned
!!  rad(mmax)=radial grid values
!!  vpspll(mmax,lnmax)=nonlocal projectors for each (l,n) on linear
!!   radial grid.  Here, these are the  product of the reference
!!   wave functions and (v(l,n)-vloc), calculated in the psp generation
!!   program and normalized so that integral(0,rc(l)) vpsll^2 dr = 1,
!!   which leads to the the usual convention for the energies ekb(l,n)
!!   also calculated in the psp generation program.
!!
!! OUTPUT
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_ln(q) and
!!   second derivative from spline fit for each (l,n).
!!
!! NOTES
!! u_l(r) is reference state wavefunction (input as wf);
!! j_l(q) is a spherical Bessel function;
!! dV_l(r) = vpsp_l(r)-vloc(r) for angular momentum l;
!! f_l(q) = $ \int_0^{rmax}[j_l(2\pi q r) u_l(r) dV_l(r) r dr]/\sqrt{dvms}$
!! where dvms = $\int_0^{rmax} [(u_l(r) dV_l(r))^2 dr]$ is the mean
!! square value of the nonlocal correction for angular momentum l.
!! Xavier Gonze s E_KB = $ dvms/\int_0^{rmax}[(u_l(r))^2 dV_l(r) dr]$.
!! This is the eigenvalue of the Kleinman-Bylander operator and sets
!! the energy scale of the nonlocal psp corrections.
!!
!! PARENTS
!!      m_psp8,m_psp9
!!
!! CHILDREN
!!      spline,splint
!!
!! SOURCE

subroutine psp8nl(amesh,ffspl,indlmn,lmax,lmnmax,lnmax,mmax,mqgrid,qgrid,rad,vpspll)

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: lmax,lmnmax,lnmax,mmax,mqgrid
 real(dp),intent(in) :: amesh
!arrays
 integer,intent(in) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vpspll(mmax,lnmax)
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax) !vz_i

!Local variables-------------------------------
!Following parameter controls accuracy of Fourier transform based on qmax
!and represents the minimun number of integration points in one period.
!scalars
 integer,parameter :: NPT_IN_2PI=200
 integer :: iln,iln0,ilmn,iq,ir,irmu,irn,ll,mesh_mult,mmax_new,mvpspll
 real(dp) :: amesh_new,arg,c1,c2,c3,c4,dri,qmesh,result,tv,xp,xpm1,xpm2,xpp1
 real(dp) :: yp1,ypn
!arrays
 real(dp) :: sb_out(4)
 real(dp),allocatable :: rad_new(:),vpspll_new(:,:),work(:,:),work2(:)

! *************************************************************************

!Find r mesh spacing necessary for accurate integration at qmax
 amesh_new=2.d0*pi/(NPT_IN_2PI*qgrid(mqgrid))

!Choose submultiple of input mesh
 mesh_mult=int(amesh/amesh_new) + 1
 mmax_new=mesh_mult*(mmax-1)+1
 amesh_new=amesh/dble(mesh_mult)

 ABI_ALLOCATE(rad_new,(mmax_new))
 ABI_ALLOCATE(vpspll_new,(mmax_new,lnmax))

 if(mesh_mult==1) then
   rad_new(:)=rad(:)
 else
!  Set up new radial mesh
   irn=1
   do ir=1,mmax-1
     do irmu=0,mesh_mult-1
       rad_new(irn)=rad(ir)+dble(irmu)*amesh_new
       irn=irn+1
     end do
   end do
   rad_new(mmax_new)=rad(mmax)
 end if

!Interpolate projectors onto new grid if called for
!Cubic polynomial interpolation is used which is consistent
!with the original interpolation of these functions from
!a log grid to the input linear grid.
 dri = one/amesh
 do irn=1,mmax_new
!  index to find bracketing input mesh points
   if(mesh_mult>1) then
     ir = irn/mesh_mult + 1
     ir = max(ir,2)
     ir = min(ir,mmax-2)
!    interpolation coefficients
     xp = dri * (rad_new(irn) - rad(ir))
     xpp1 = xp + one
     xpm1 = xp - one
     xpm2 = xp - two
     c1 = -xp * xpm1 * xpm2 * sixth
     c2 = xpp1 * xpm1 * xpm2 * half
     c3 = - xp * xpp1 * xpm2 * half
     c4 = xp * xpp1 * xpm1 * sixth
!    Now do the interpolation on all projectors for this grid point

     iln0=0
     do ilmn=1,lmnmax
       iln=indlmn(5,ilmn)
       if (iln>iln0) then
         iln0=iln
         tv =  c1 * vpspll(ir - 1, iln) &
&         + c2 * vpspll(ir    , iln) &
&         + c3 * vpspll(ir + 1, iln) &
&         + c4 * vpspll(ir + 2, iln)
         if(abs(tv)>tol10) then
           vpspll_new(irn,iln)=tv
           mvpspll=irn
         else
           vpspll_new(irn,iln)=zero
         end if
       end if
     end do

   else
!    With no mesh multiplication, just copy projectors
     ir=irn
     iln0=0
     do ilmn=1,lmnmax
       iln=indlmn(5,ilmn)
       if (iln>iln0) then
         iln0=iln
         tv = vpspll(ir,iln)
         if(abs(tv)>tol10) then
           vpspll_new(irn,iln)=tv
           mvpspll=irn
         else
           vpspll_new(irn,iln)=zero
         end if
       end if
     end do

   end if
 end do !irn

 ABI_ALLOCATE(work,(mvpspll,lnmax))

!Loop over q values
 do iq=1,mqgrid
   arg=2.d0*pi*qgrid(iq)

!  Set up integrands
   do  ir=1,mvpspll
     call sbf8(lmax+1,arg*rad_new(ir),sb_out)
     iln0=0
     do ilmn=1,lmnmax
       iln=indlmn(5,ilmn)
       if (iln>iln0) then
         iln0=iln
         ll=indlmn(1,ilmn)
         work(ir,iln)=sb_out(ll+1)*vpspll_new(ir,iln)*rad_new(ir)
       end if
     end do
   end do !ir

!  Do integral from zero to rad_new(mvpspll)
   iln0=0
   do ilmn=1,lmnmax
     iln=indlmn(5,ilmn)
     if (iln>iln0) then
       iln0=iln
       call ctrap(mvpspll,work(1,iln),amesh_new,result)
       ffspl(iq,1,iln)=result
     end if
   end do

!  End loop over q mesh
 end do !iq

!Fit splines for form factors
 ABI_ALLOCATE(work2,(mqgrid))
 qmesh=qgrid(2)-qgrid(1)

 iln0=0
 do ilmn=1,lmnmax
   iln=indlmn(5,ilmn)
   if (iln>iln0) then
     iln0=iln
!    Compute derivatives of form factors at ends of interval
     yp1=(-50.d0*ffspl(1,1,iln)+96.d0*ffspl(2,1,iln)-72.d0*ffspl(3,1,iln)&
&     +32.d0*ffspl(4,1,iln)- 6.d0*ffspl(5,1,iln))/(24.d0*qmesh)
     ypn=(6.d0*ffspl(mqgrid-4,1,iln)-32.d0*ffspl(mqgrid-3,1,iln)&
&     +72.d0*ffspl(mqgrid-2,1,iln)-96.d0*ffspl(mqgrid-1,1,iln)&
&     +50.d0*ffspl(mqgrid,1,iln))/(24.d0*qmesh)

     call spline(qgrid,ffspl(1,1,iln),mqgrid,yp1,ypn,ffspl(1,2,iln))
   end if
 end do

 ABI_DEALLOCATE(rad_new)
 ABI_DEALLOCATE(vpspll_new)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(work2)

end subroutine psp8nl
!!***

!!****f* ABINIT/cc_derivatives
!! NAME
!! cc_derivatives
!!
!! FUNCTION
!! subroutine to spline the core charge and get derivatives
!!   extracted from previous version of psp6cc_drh
!! input on log grid, and splined to regular grid between 0 and rchrg
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!  rad=radial grid points
!!  ff=core charge at points in rad
!!  ff1=first derivative of ff on log grid
!!  ff2=second derivative of ff on log grid
!!
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! PARENTS
!!      m_psp6,m_upf2abinit
!!
!! CHILDREN
!!      spline,splint
!!
!! NOTES
!! Test version by DRH - requires very smooth model core charge
!!
!! SOURCE

subroutine cc_derivatives(rad,ff,ff1,ff2,mmax,n1xccc,rchrg,xccc1d)

!Arguments ------------------------------------
! scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg
!arrays
 real(dp),intent(in) :: rad(mmax),ff(mmax),ff1(mmax),ff2(mmax)
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
! scalars
 integer :: i1xccc
 real(dp) :: der1,dern
!arrays
 real(dp),allocatable :: ff3(:),ff4(:),gg(:),gg1(:),gg2(:)
 real(dp),allocatable :: gg3(:),gg4(:),work(:),xx(:)

! *************************************************************************
 ABI_ALLOCATE(ff3,(mmax))
 ABI_ALLOCATE(ff4,(mmax))
 ABI_ALLOCATE(gg,(n1xccc))
 ABI_ALLOCATE(gg1,(n1xccc))
 ABI_ALLOCATE(gg2,(n1xccc))
 ABI_ALLOCATE(gg3,(n1xccc))
 ABI_ALLOCATE(gg4,(n1xccc))
 ABI_ALLOCATE(work,(mmax))
 ABI_ALLOCATE(xx,(n1xccc))

 !write(std_out,*) 'cc_derivatives : enter'

!calculate third derivative ff3 on logarithmic grid
 der1=ff2(1)
 dern=ff2(mmax)
 call spline(rad,ff1,mmax,der1,dern,ff3)

!calculate fourth derivative ff4 on logarithmic grid
 der1=0.d0
 dern=0.d0
 call spline(rad,ff2,mmax,der1,dern,ff4)

!generate uniform mesh xx in the box cut by rchrg:

 do i1xccc=1,n1xccc
   xx(i1xccc)=(i1xccc-1)* rchrg/dble(n1xccc-1)
 end do
!
!now interpolate core charge and derivatives on the uniform grid
!
!core charge, input=ff,  output=gg
 call splint(mmax,rad,ff,ff2,n1xccc,xx,gg)

!first derivative input=ff1, output=gg1
 call splint(mmax,rad,ff1,ff3,n1xccc,xx,gg1)

!normalize gg1
!gg1(:)=gg1(:)*rchrg

!second derivative input=ff2, output=gg2
 call splint(mmax,rad,ff2,ff4,n1xccc,xx,gg2)

!normalize gg2
!gg2(:)=gg2(:)*rchrg**2

!reallocate work otherwise the calls to spline crash (n1xccc /= mmax)
 ABI_DEALLOCATE(work)
 ABI_ALLOCATE(work,(n1xccc))

!recalculate 3rd derivative consistent with spline fit to first derivative
!on linear grid
 der1=gg2(1)
 dern=gg2(n1xccc)
 call spline(xx,gg1,n1xccc,der1,dern,gg3)

!calculate 4th derivative consistent with spline fit to second derivative
!on linear grid
 der1=0.0d0
 dern=0.0d0
 call spline(xx,gg2,n1xccc,der1,dern,gg4)

!now calculate second to fourth derivative by forward differences
!to avoid numerical noise uses a smoothing function
!
!call smooth(gg1,n1xccc,10)

!gg2(n1xccc)=0.0
!do i1xccc=1,n1xccc-1
!gg2(i1xccc)=(gg1(i1xccc+1)-gg1(i1xccc))*dble(n1xccc-1)
!end do

!call smooth(gg2,n1xccc,10)

!gg3(n1xccc)=0.0
!do i1xccc=1,n1xccc-1
!gg3(i1xccc)=(gg2(i1xccc+1)-gg2(i1xccc))*dble(n1xccc-1)
!end do

!call smooth(gg3,n1xccc,10)

!gg4(n1xccc)=0.0
!do i1xccc=1,n1xccc-1
!gg4(i1xccc)=(gg3(i1xccc+1)-gg3(i1xccc))*dble(n1xccc-1)
!end do

!call smooth(gg4,n1xccc,10)

!write on xcc1d
!normalize to unit range usage later in program
 xccc1d(:,1)=gg(:)
 xccc1d(:,2)=gg1(:)*rchrg
 xccc1d(:,3)=gg2(:)*rchrg**2
 xccc1d(:,4)=gg3(:)*rchrg**3
 xccc1d(:,5)=gg4(:)*rchrg**4
!***drh test
!write(std_out,'(a,2i6)') 'drh:psp6cc_drh - mmax,n1xccc',mmax,n1xccc
!***end drh test


!DEBUG
!note: the normalization condition is the following:
!4pi rchrg /dble(n1xccc-1) sum xx^2 xccc1d(:,1) = qchrg
!
!norm=0.d0
!do i1xccc=1,n1xccc
!norm = norm + 4.d0*pi*rchrg/dble(n1xccc-1)*&
!&             xx(i1xccc)**2*xccc1d(i1xccc,1)
!end do
!write(std_out,*) ' norm=',norm
!
!write(std_out,*)' psp6cc_drh : output of core charge density and derivatives '
!write(std_out,*)'   xx          gg           gg1  '
!do i1xccc=1,n1xccc
!write(10, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   xx          gg2          gg3  '
!do i1xccc=1,n1xccc
!write(11, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   xx          gg4          gg5  '
!do i1xccc=1,n1xccc
!write(12, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 ABI_DEALLOCATE(ff3)
 ABI_DEALLOCATE(ff4)
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(gg1)
 ABI_DEALLOCATE(gg2)
 ABI_DEALLOCATE(gg3)
 ABI_DEALLOCATE(gg4)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(xx)

end subroutine cc_derivatives
!!***

end module m_psptk
!!***

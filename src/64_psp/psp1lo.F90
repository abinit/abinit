!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp1lo
!! NAME
!! psp1lo
!!
!! FUNCTION
!! Compute sine transform to transform from v(r) to q^2 v(q)
!! using subroutines related to Teter atomic structure grid.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  drad(mmax)=inverse of r grid spacing at each point
!!  mmax=number of radial r grid points (Teter grid)
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rad(mmax)=r grid values (bohr).
!!  vloc(mmax)=v(r) on radial grid.
!!  wksincos(mmax,2,2)=contains sine and cosine of 2*pi*r(:)*dq and 2*pi*r(:)*q
!!    at input :  wksincos(:,1,1)=sine of 2*pi*r(:)*dq
!!                wksincos(:,2,1)=cosine of 2*pi*r(:)*dq
!!    wksincos(:,:,2) is not initialized, will be used inside the routine
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm= $4\pi \int[r^2 (v(r)+Zv/r) dr]$
!!  q2vq(mqgrid)=$q^2 v(q)$
!!  =$\displaystyle -Zv/\pi+q^2 4\pi\int(\frac{\sin(2\pi q r)}{2 \pi q r})(r^2 v(r)+r Zv)dr$.
!!  yp1,ypn=derivative of q^2 v(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      psp1in
!!
!! CHILDREN
!!      der_int,sincos
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psp1lo(drad,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&
&  vloc,wksincos,yp1,ypn,zion)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp1lo'
 use interfaces_64_psp, except_this_one => psp1lo
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,mqgrid
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm,yp1,ypn
!arrays
 real(dp),intent(in) :: drad(mmax),qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(inout) :: wksincos(mmax,2,2)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables-------------------------------
!scalars
 integer,parameter :: mma0=2001
 integer :: iq,ir,irmax
 real(dp),parameter :: scale=10.0d0
 real(dp) :: result,test,tpiq
!arrays
 real(dp) :: wk(mma0),wk1(mma0),wk2(mma0)

! *************************************************************************

!DEBUG
!write(std_out,*)' psp1lo : enter '
!if(.true.)stop
!ENDDEBUG

!Do q=0 separately (compute epsatm)
!Set up integrand for q=0: Int[r^2 (V(r)+Zv/r) dr]
!Treat r=0 by itself
 wk(1)=0.0d0

 do ir=2,mmax
!  (at large r do not want prefactor of r^2 and should see
!  V(r)+Zv/r go to 0 at large r)
   test=vloc(ir)+zion/rad(ir)
!  DEBUG
!  write(std_out,'(i4,3es20.10)' )ir,rad(ir),test,rad(ir)*test
!  ENDDEBUG
!  In this routine, NO cut-off radius is imposed : the input
!  vloc MUST be in real(dp) to obtain numerically
!  accurate values. The error can be on the order of 0.001 Ha !
   if (abs(test)<1.0d-20) then
     wk(ir)=0.0d0
   else
     wk(ir)=rad(ir)*(rad(ir)*vloc(ir)+zion)
   end if
 end do
!Do integral from 0 to r(max) (disregard contrib beyond r(max)
!(need numerical derivatives to do integral)
!Use mmax-1 to convert to Teter s dimensioning starting at 0
 call der_int(wk,wk2,rad,drad,mmax-1,result)

!DEBUG
!write(std_out,*)' psp1lo : result ',result
!stop
!ENDDEBUG

 epsatm=4.d0*pi*(result)
!q=0 value of integral is -zion/Pi + q^2 * epsatm = -zion/Pi
 q2vq(1)=-zion/pi

!Prepare loop over q values
 irmax=mmax+1
 do ir=mmax,2,-1
   test=vloc(ir)+zion/rad(ir)
   wk1(ir)=test*rad(ir)
!  Will ignore tail within decade of machine precision
   if ((scale+abs(test))==scale .and. irmax==ir+1) then
     irmax=ir
   end if
 end do
!Increase irmax a bit : this is copied from psp1nl
 irmax=irmax+4
 if(irmax>mmax)irmax=mmax

!Loop over q values
 do iq=2,mqgrid
   tpiq=two_pi*qgrid(iq)
   call sincos(iq,irmax,mmax,wksincos,rad,tpiq)
!  set up integrand Sin(2Pi q r)(rV(r)+Zv) for integral
!$\displaystyle -Zv/\pi + q^2 4\pi \int[\frac{\sin(2\pi q r)}{2\pi q r}(r^2 v(r)+r Zv)dr]$.
!  Handle r=0 separately
   wk(1)=0.0d0
   do ir=2,irmax
     wk(ir)=wksincos(ir,1,2)*wk1(ir)
   end do
!  do integral from 0 to r(max)
   if(irmax>mmax-1)irmax=mmax-1

   call der_int(wk,wk2,rad,drad,irmax,result)
!  store q^2 v(q)
   q2vq(iq)=-zion/pi+2.d0*qgrid(iq)*result
 end do

!Compute derivatives of q^2 v(q) at ends of interval
 yp1=0.0d0
!ypn=$\displaystyle 2\int_0^\infty (\sin (2\pi qmax r)+(2\pi qmax r)\cos (2\pi qmax r)(r V(r)+Z)dr]$
!integral from r(mmax) to infinity is overkill; ignore
!set up integrand
!Handle r=0 separately
 wk(1)=0.0d0
 tpiq=two_pi*qgrid(mqgrid)
 do ir=2,mmax
   test=vloc(ir)+zion/rad(ir)
!  Ignore contributions within decade of machine precision
   if ((scale+abs(test))==scale) then
     wk(ir)=0.0d0
   else
     wk(ir)=(sin(tpiq*rad(ir))+tpiq*rad(ir)*cos(tpiq*rad(ir))) * &
&     (rad(ir)*vloc(ir)+zion)
   end if
 end do
 call der_int(wk,wk2,rad,drad,mmax-1,result)

 ypn=2.0d0*result

end subroutine psp1lo
!!***

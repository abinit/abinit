!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp5lo
!! NAME
!! psp5lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on logarithmic grid using related uniform
!! grid in exponent and corrected trapezoidal integration.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, FrD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      psp5in,psp6in
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp5lo(al,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&
&                  vloc,yp1,ypn,zion)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp5lo'
 use interfaces_32_util
!End of the abilint section

 implicit none

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

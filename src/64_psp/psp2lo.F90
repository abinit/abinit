!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp2lo
!! NAME
!! psp2lo
!!
!! FUNCTION
!! Treat local part of Goedecker-Teter-Hutter pseudopotentials (pspcod=2),
!! as well as Hartwigsen-Goedecker-Hutter pseudopotentials (pspcod=3)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cc1,2,3,4=parameters from analytic pseudopotential form
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=values of q (or G) on grid from 0 to qmax (bohr^-1)
!!                if vlspl_recipSpace is .true. else values of r on grid from
!!                0 to 2pi / qmax * mqgrid_ff (bohr).
!!  rloc=local pseudopotential core radius (bohr)
!!  vlspl_recipSpace= .true. if computation of vlspl is done in reciprocal space
!!  zion=valence charge of atom
!!  parameters for local potential: rloc,c1,c2,c3,c4
!!
!! OUTPUT
!!  dvloc(mqgrid)=dVloc(r)/dr (only allocated if vlspl_recipSpace is false).
!!  epsatm=$4\pi\int[r^2 (v(r)+\frac{Zv}{r} dr]$
!!{{\ \begin{eqnarray}
!!  q2vq(mqgrid)&=&q^2 v(q) \nonumber \\
!!  &=&-Zv/\pi
!!   +q^2 4\pi\int[(\frac{\sin(2\pi qr)}{2\pi qr})(r^2 v(r)+r Zv)dr]\nonumber\\
!!  &=&\exp(-K^2*rloc^2/2) \nonumber \\
!!  &&   *(-\frac{zion}{\pi}+(\frac{K^2*rloc^3}{\sqrt{2*\pi}}*
!!       (c1+c2*(3-(rloc*K)^2) \nonumber \\
!!  &&    +c3*(15-10(rloc*K)^2+(rloc*K)^4) \nonumber \\
!!  &&    +c4*(105-105*(rloc*K)^2+21*(rloc*K)^4-(rloc*K)^6)) \nonumber
!!\end{eqnarray} }}
!! for GTH vloc with $K=(2\pi q)$.
!!  yp1,ypn=derivative of q^2 v(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      psp10in,psp2in,psp3in
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psp2lo(cc1,cc2,cc3,cc4,dvloc,epsatm,mqgrid,qgrid,q2vq,&
&  rloc,vlspl_recipSpace,yp1,ypn,zion)

 use defs_basis
 use m_profiling_abi

 use m_special_funcs,  only : abi_derfc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2lo'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mqgrid
 real(dp),intent(in) :: cc1,cc2,cc3,cc4,rloc,zion
 real(dp),intent(out) :: epsatm,yp1,ypn
 logical,intent(in) :: vlspl_recipSpace
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: dvloc(mqgrid),q2vq(mqgrid)

!Local variables-------------------------------
!scalars
 integer :: iqgrid
 real(dp) :: erfValue,gaussValue,polyValue,qmax,rq,rq2
 character(len=500) :: message

! *************************************************************************

!Compute epsatm = lim(q->0) [Vloc(q) + zion/(Pi*q^2)]
 epsatm=2.d0*pi*rloc**2*zion+(2.d0*pi)**(1.5d0)*rloc**3*&
& (cc1+3.d0*cc2+15.d0*cc3+105.d0*cc4)

!If vlspl_recipSpace is .true., we compute V(q)*q^2 in reciprocal space,
!else we compute V(r) in real space.
 if (vlspl_recipSpace) then
   write(message, '(a)' ) '-  Local part computed in reciprocal space.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  d(q^2*V(q))/d(q) at q=0 and q=qmax
   qmax=qgrid(mqgrid)
   rq2=(2.d0*pi*qmax*rloc)**2
   yp1=0.d0
   ypn= (2.d0*pi*qmax*rloc**2)*exp(-0.5d0*rq2)* &
&   (2.d0*zion + sqrt(2.d0*pi)*rloc*&
&   (cc1*(2.d0-rq2) + cc2*(6.d0-7.d0*rq2+rq2**2) +&
&   cc3*(30.d0-55.d0*rq2+16.d0*rq2**2-rq2**3) +&
&   cc4*(210.d0-525.d0*rq2+231.d0*rq2**2-29.d0*rq2**3+rq2**4)))
!  ypn has been tested against Maple-derived expression.

!  Compute q^2*vloc(q) on uniform grid
   do iqgrid=1,mqgrid
     rq2=(2.d0*pi*qgrid(iqgrid)*rloc)**2
     q2vq(iqgrid)=exp(-0.5d0*rq2)*(-zion/pi+rq2*(rloc/sqrt(2.d0*pi)) *&
&     ( cc1 + cc2*(3.d0-rq2) + cc3*(15.d0-10.d0*rq2+rq2**2) +&
&     cc4*(105.d0-rq2*(105.d0-rq2*(21.d0-rq2)))  ))
   end do
 else
   write(message, '(a)' ) '-  Local part computed in real space.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   
!  Compute derivatives for splines computations
   yp1 = 0.d0
   rq2 = (qgrid(mqgrid) / rloc) ** 2
   erfValue = abi_derfc(sqrt(0.5d0 * rq2))
   ypn = - 2.0d0 * zion / sqrt(2.d0 * pi) / qgrid(mqgrid) / rloc
   ypn = ypn - rq2 * (cc1 + cc2 * rq2 + cc3 * rq2 ** 2 + cc4 * rq2 ** 3) / qgrid(mqgrid)
   ypn = ypn + (2.d0 * cc2 * rq2 + 4.d0 * cc3 * rq2 ** 2 + 6.d0 * cc4 * rq2 ** 3) / qgrid(mqgrid)
   ypn = ypn * exp(-0.5d0 * rq2)
   ypn = ypn + zion / qgrid(mqgrid) ** 2 * erfValue
!  Note that ypn has been calculated on a full-proof a4 paper sheet.

!  Compute local potential and its first derivatives.
   do iqgrid = 1, mqgrid, 1
     rq2 = (qgrid(iqgrid) / rloc) ** 2
!    Compute erf() part
!    Case r = 0
     gaussValue = exp(-0.5d0 * rq2)
     if (qgrid(iqgrid) == 0.d0) then
       q2vq(iqgrid) = -zion / rloc * sqrt(2.d0 / pi)
       dvloc(iqgrid) = 0.d0
     else
       erfValue = abi_derfc(sqrt(0.5d0 * rq2))
       q2vq(iqgrid) = -zion / qgrid(iqgrid) * (1.0d0 - erfValue)
       dvloc(iqgrid) = - sqrt(2.d0 / pi) * zion * gaussValue / (qgrid(iqgrid) * rloc) - &
&       q2vq(iqgrid) / qgrid(iqgrid)
     end if
!    Add the gaussian part
     polyValue = cc1 + cc2 * rq2 + cc3 * rq2 ** 2 + cc4 * rq2 ** 3
     q2vq(iqgrid) = q2vq(iqgrid) + gaussValue * polyValue
     rq = qgrid(iqgrid) / rloc
     dvloc(iqgrid) = dvloc(iqgrid) - qgrid(iqgrid) / rloc ** 2 * gaussValue * polyValue + &
&     gaussValue * (2.0d0 * cc2 * rq / rloc + 3.0d0 * cc3 * rq ** 3 / rloc + &
&     6.0d0 * cc4 * rq ** 5 / rloc)
   end do

   write(message, '(a,f12.7,a,a,f12.7,a,a,a,f12.7)' ) &
&   '  | dr spline step is : ', qgrid(2), ch10, &
&   '  | r > ', qgrid(mqgrid) ,' is set to 0.', ch10, &
&   '  | last non-nul potential value is : ', q2vq(mqgrid)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

end subroutine psp2lo
!!***

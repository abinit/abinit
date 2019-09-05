!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xciit
!! NAME
!!  m_xciit
!!
!! FUNCTION
!! Exchange-correlation at finite temperature of an electron gas
!! Ichimaru S., Iyetomi H., Tanaka S., Phys. Rep. 149, 91-205 (1987) [[cite:Ichimaru1987]]
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2019 ABINIT group (JFD,LK)
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

module m_xciit

 use defs_basis
 use m_errors

 implicit none

 private
!!***

 public :: xciit
!!***

contains
!!***

!!****f* ABINIT/xciit
!! NAME
!!  xciit
!!
!! FUNCTION
!! Exchange-correlation at finite temperature of an electron gas
!! Ichimaru S., Iyetomi H., Tanaka S., Phys. Rep. 149, 91-205 (1987) [[cite:Ichimaru1987]]
!!
!! INPUTS
!!  temp= (electronic) temperature
!!  npt=number of real space points
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=Wigner-Seitz radii at each point
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density
!!  fxc(npt)=exchange-correlation free energy at finite temperature
!!  vxc(npt)=exchange-correlation potential
!!  --- optional output ---
!!  [dvxc(npt)]=partial second derivatives of the xc energy
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine xciit(exc,fxc,npt,order,rspts,temp,vxc, &
&                dvxc)!Optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
 real(dp),intent(in) :: temp
!arrays
 real(dp),intent(in) :: rspts(npt)
 real(dp),intent(out) :: exc(npt),fxc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!scalars
 integer :: ipt
 real(dp) :: ef,deltavxc,Gamma,rs,rsm1,tt
 character(len=500) :: msg

! *************************************************************************

!Checks the values of order
 if(order<0.or.order>2)then
   write(msg, '(4a,i3,a)' ) ch10,&
&   'With Ishimaru-Iyetomi-Tanka xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be ',order,'.'
   MSG_BUG(msg)
 end if

!Loop over grid points
 do ipt=1,npt

   rs=rspts(ipt)
   rsm1=one/rs
!  Step for the Vxc computation
   deltavxc=0.01_dp
!  Compute ef
   ef=0.5_dp*(9.0_dp*pi/4.0_dp)**(2.0_dp/3.0_dp)*rsm1**2
!  Compute temperature
   tt=max(temp/ef,tol12)
!  Compute Gamma
   Gamma=one/(tt*ef)/rs

!  Exchange-correlation of Ichimaru functional
   fxc(ipt)= fexsGamma(Gamma,tt)*rsm1
   exc(ipt)=fxc(ipt) - tdexcsdtiit(rs,tt);
   vxc(ipt)=(8.0_dp*(Fxc_iit(rs,tt,deltavxc)-Fxc_iit(rs,tt,-deltavxc)) &
&   -(Fxc_iit(rs,tt,two*deltavxc)-Fxc_iit(rs,tt,-two*deltavxc)))/ &
&   (12.0_dp*deltavxc*3.0_dp/(4.0_dp*pi)/rs**3)
   if (order==2) then
     dvxc(ipt)=(-30.0_dp*Fxc_iit(rs,tt,zero)+16.0_dp*(Fxc_iit(rs,tt,deltavxc)+Fxc_iit(rs,tt,-deltavxc)) &
&     -(Fxc_iit(rs,tt,two*deltavxc)+Fxc_iit(rs,tt,two*deltavxc)))/ &
&     (12.0_dp*(deltavxc*3.0_dp/(4.0_dp*pi)/rs**3)**2)
   end if
 end do

 CONTAINS
!!***

!!****f* ABINIT/fexsGamma
!!
!! NAME
!! fexsGamma
!!
!! FUNCTION
!! Free energy for the IIT finite temperautre XC functional
!!
!! INPUTS
!! Gamma= ?
!! t=temperature
!!
!! OUTPUT
!! fexsGamma=free energy
!!
!! SOURCE

 function fexsGamma(Gamma,t)

!Arguments ------------------------------------
 real(dp) :: fexsGamma
 real(dp),intent(in) :: Gamma,t
!Local variables-------------------------------
 real(dp) :: lambda
 real(dp) :: tanht,tanhst
 real(dp) :: a,b,c,d,e
 real(dp) :: bmcdse,amcse,sqrt4emd2

! *************************************************************************

   lambda=(4.0_dp/(9.0_dp*pi))**(one/3.0_dp)
   tanht=tanh(one/t)
   tanhst=tanh(one/sqrt(t))

   a=one/(pi*lambda)*(0.75_dp+3.04363_dp*t**2-0.09227_dp*t**3+1.7035_dp*t**4)/(one+8.31051_dp*t**2+5.1105_dp*t**4)*tanht
   b=(0.341308_dp+12.070873_dp*t**2+1.148889_dp*t**4)/(one+10.495346_dp*t**2+1.326623_dp*t**4)*sqrt(t)*tanhst
   e=(0.539409_dp+2.522206_dp*t**2+0.178484_dp*t**4)/(one+2.555501_dp*t**2+0.146319_dp*t**4)*t*tanht
   c=(0.872496_dp+0.025248_dp*exp(-1./t))*e
   d=(0.614925_dp+16.996055_dp*t**2+1.489056_dp*t**4)/(one+10.10935_dp*t**2+1.22184_dp*t**4)*sqrt(t)*tanhst

   bmcdse=b-c*d/e
   amcse=a-c/e
   sqrt4emd2=sqrt(4.0_dp*e-d**2)

   fexsGamma=-one/Gamma*(c/e*Gamma+2.0_dp/e*bmcdse*sqrt(Gamma)+one/e*(amcse-d/e*bmcdse)*log(e*Gamma+d*sqrt(Gamma)+one)- &
&   2.0_dp/(e*sqrt4emd2)*(d*amcse+(2.0_dp-d**2/e)*bmcdse)*(atan((2.0_dp*e*sqrt(Gamma)+d)/sqrt4emd2)-atan(d/sqrt4emd2)))

 end function fexsGamma
!!***

!!****f* ABINIT/Fxc_iit
!!
!! NAME
!! Fxc_iit
!!
!! FUNCTION
!! Auxiliary function for the IIT finite temperature XC functional
!!
!! INPUTS
!! deltavxc= ?
!! rs=Wigner-Seitz radius
!! t=temperature
!!
!! OUTPUT
!! Fxc_iit=auxiliary function
!!
!! SOURCE

 function Fxc_iit(rs,t,deltavxc)

!Arguments ------------------------------------
 real(dp) :: Fxc_iit
 real(dp),intent(in) :: rs,t,deltavxc
!Local variables-------------------------------
 real(dp) :: newrs,newt,newGamma

! *************************************************************************

   newrs=rs/(one+deltavxc)**(one/3.0_dp)
   newt=t/(one+deltavxc)**(2.0_dp/3.0_dp)
   newGamma=2.0_dp*(4.0_dp/(9.0_dp*pi))**(2.0_dp/3.0_dp)*newrs/newt
   Fxc_iit=3.0_dp/(4.0_dp*pi)*fexsGamma(newGamma,newt)/newrs**4

 end function Fxc_iit
!!***

!!****f* ABINIT/tdexcsdtiit
!!
!! NAME
!! tdexcsdtiit
!!
!! FUNCTION
!! Auxiliary function for the IIT finite temperature XC functional
!!
!! INPUTS
!! rs=Wigner-Seitz radius
!! t=temperature
!!
!! OUTPUT
!! tdexcsdtiit=auxiliary function
!!
!! SOURCE

 function tdexcsdtiit(rs,t)

!Arguments ------------------------------------
 real(dp) :: tdexcsdtiit
 real(dp),intent(in) :: rs,t
!Local variables-------------------------------
 real(dp) :: ef,Gamma
 real(dp) :: deltat=1.0d-2

! *************************************************************************

   ef=half*(9.0_dp*pi/4.0_dp)**(2.0_dp/3.0_dp)/rs**2
   Gamma=one/(t*ef)/rs
   tdexcsdtiit=8.0_dp*(fexsGamma(Gamma/(one+deltat),(one+deltat)*t) &
&   -fexsGamma(Gamma/(one-deltat),(one-deltat)*t)) &
&   -(fexsGamma(Gamma/(one+2.0_dp*deltat),(one+2.0_dp*deltat)*t) &
&   -fexsGamma(Gamma/(one-2.0_dp*deltat),(one-2.0_dp*deltat)*t))
   tdexcsdtiit=t*tdexcsdtiit/(12.0_dp*deltat*t)/rs

 end function tdexcsdtiit
!!***

end subroutine xciit
!!***

end module m_xciit
!!***

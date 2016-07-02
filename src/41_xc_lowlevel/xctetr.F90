!{\src2tex{textfont=tt}}
!!****f* ABINIT/xctetr
!! NAME
!! xctetr
!!
!! FUNCTION
!! Returns exc, vxc, and d(vxc)/d($\rho$) from input $\rho$.
!! Also returns $d^2(Vxc)/d(\rho)^2$ as needed for third-order DFPT
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rhor(npt)=electron number density (bohr^-3)
!!  rspts(npt)=corresponding Wigner-Seitz radii, precomputed
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d(rho*exc)/d(rho)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!  if(order>2) d2vxc(npt)=derivative d$^2$(Vxc)/d$(\rho)^2$ (hartree*bohr^6)
!!
!! NOTES
!! Teter exchange and correlation (xc)--Mike Teter s fit
!! to Ceperly-Alder electron gas energy data.  Data from
!! D.M. Ceperley and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
!! and private communication from authors.
!! This form is based on Mike Teter s rational polynomial
!! exc=-(a0+a1*rs+a2*rs**2+a3*rs**3)/(b1*rs+b2*rs**2+b3*rs**3+b4*rs**4)
!! where the parameters of the fit are fit to reproduce
!! Ceperley-Alder data and the high density limit (rs->0)
!! of the electron gas (pure exchange).
!! rs = $(3/(4\pi))^{1/3} * \rho(r)^{-1/3}$.
!! b1 must be 1 and a0 must be $(3/4)(3/(2\pi))^{2/3}$.
!! Fit is by Mike Teter, Corning Incorporated.
!! Note that d(vxc)/d($\rho$) gets a little wild at small rho.
!! d$^2$(Vxc)/d$(\rho)^2$ is probably wilder.
!!
!! Some notation:  (XG 990224, sign convention should be changed, see xcspol.f)
!!  $Exc = N1/D1$ with $N1=-(a0+a1*rs+...)$ given above and
!!              $D1= (b1*rs+b2*rs^2+...)$ also given above.
!!  $Vxc = N2/D1^2$ with $N2=d(N1)/d(rs)$.
!!  $d(Vxc)/d(rs)=(N3-D3*(2*N2/D1))/D1^2 with N3=d(N2)/d(rs)$ and
!!              $D3=d(D1)/d(rs)$.
!!  $d(Vxc)/d(\rho) = (-rs/(3*\rho))* d(Vxc)/d(rs)$.
!!  $d^2(Vxc)/d(rs)^2 = (N4-2*(2*N3*D3+N2*D4-3*N2*D3^2/D1)/D1)/D1^2$
!!   with $N4=d(N3)/d(rs), D4=d(D3)/d(rs)$.
!!  $d^2(Vxc)/d(\rho)^2= rs/(3*\rho)^2)*(4*d(Vxc)/d(rs)+rs*d^2(Vxc)/d(rs)^2)$.
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xctetr(exc,npt,order,rhor,rspts,vxc,& !Mandatory arguments
&                 d2vxc,dvxc)                    !Optional arguments

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xctetr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rhor(npt),rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: d2vxc(npt),dvxc(npt)

!Local variables-------------------------------
!rsfac=(3/(4 Pi))^(1/3)
!Mike Teter s parameters: (keep 8 digits after decimal)
!(a0=(3/4)(3/(2 Pi))^(2/3)
!scalars
 integer :: ipt
 real(dp),parameter :: a0=.4581652932831429_dp,a1=2.40875407_dp,a2=.88642404_dp
 real(dp),parameter :: a3=.02600342_dp,b1=1.0_dp,b2=4.91962865_dp
 real(dp),parameter :: b3=1.34799453_dp,b4=.03120453_dp,c1=4._dp*a0*b1/3.0_dp
 real(dp),parameter :: c2=5.0_dp*a0*b2/3.0_dp+a1*b1
 real(dp),parameter :: c3=2.0_dp*a0*b3+4.0_dp*a1*b2/3.0_dp+2.0_dp*a2*b1/3.0_dp
 real(dp),parameter :: c4=7.0_dp*a0*b4/3.0_dp+5.0_dp*a1*b3/3.0_dp+a2*b2+a3*b1/3.0_dp
 real(dp),parameter :: c5=2.0_dp*a1*b4+4.0_dp*a2*b3/3.0_dp+2.0_dp*a3*b2/3.0_dp
 real(dp),parameter :: c6=5.0_dp*a2*b4/3.0_dp+a3*b3,c7=4.0_dp*a3*b4/3.0_dp
 real(dp),parameter :: rsfac=0.6203504908994000_dp
 real(dp) :: d1,d1m1,d2vxcr,d3,d4,dvxcdr,n1,n2,n3,n4,rhom1,rs
 character(len=500) :: message

! *************************************************************************
!
!Checks the values of order
 if(order<0 .or. order>3)then
   write(message, '(a,a,a,i6)' )&
&   'With Teter 91 Ceperley-Alder xc functional, the only',ch10,&
&   'allowed values for order are 0, 1, 2 or 3, while it is found to be',order
   MSG_BUG(message)
 end if

!Checks the compatibility between the order and the presence of the optional arguments
 if(order /=3 .and. present(d2vxc))then
   write(message, '(a,a,a,i6)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector d2vxc, that is needed only with order=3, while we have',order
   MSG_BUG(message)
 end if

 if(order <= 1 .and. present(dvxc))then
   write(message, '(a,a,a,i6)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector dvxc, that is needed with order > 1, while we have',order
   MSG_BUG(message)
 end if

!separated cases with respect to order

 if (order<=1) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     n1=-(a0+rs*(a1+rs*(a2+rs*a3)))
     d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
     d1m1=1.0_dp/d1
     n2=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
!    
!    Exchange-correlation energy
     exc(ipt)=n1*d1m1
!    
!    Exchange-correlation potential
     vxc(ipt)=n2*d1m1**2
   end do
 else if (order>2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     n1=-(a0+rs*(a1+rs*(a2+rs*a3)))
     d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
     d1m1=1.0_dp/d1
     n2=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
!    
!    Exchange-correlation energy
     exc(ipt)=n1*d1m1
!    
!    Exchange-correlation potential
     vxc(ipt)=n2*d1m1**2
!    Assemble derivative of vxc wrt rs
     n3=-(c1+rs*(2._dp*c2+rs*(3._dp*c3+rs*(4._dp*c4+rs*(5._dp*c5+&
&     rs*(6._dp*c6+rs*(7._dp*c7)))))))
     d3=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
     dvxcdr=(n3-d3*(2._dp*n2*d1m1))*d1m1**2
     rhom1=1.0_dp/rhor(ipt)
!    
!    derivative of vxc wrt rho
     dvxc(ipt)=-dvxcdr*rs*third*rhom1
!    
     
!    Assemble derivative d^2(Vxc)/d(rs)^2
     n4=-(2.0_dp*c2+rs*(6.0_dp*c3+rs*(12.0_dp*c4+rs*(20.0_dp*c5+&
&     rs*(30.0_dp*c6+rs*(42.0_dp*c7))))))
     d4=2.0_dp*b2+rs*(6.0_dp*b3+rs*(12.0_dp*b4))
     d2vxcr=(n4-2.0_dp*(2.0_dp*n3*d3+n2*d4-3.0_dp*n2*d3**2*d1m1)*d1m1)*d1m1**2
     
!    Derivative d^2(Vxc)/d(rho)^2
     d2vxc(ipt)=(rs*third*rhom1)*(4.0_dp*dvxcdr+rs*d2vxcr)*third*rhom1
     
   end do
 else if (order>1) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     n1=-(a0+rs*(a1+rs*(a2+rs*a3)))
     d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
     d1m1=1.0_dp/d1
     n2=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
!    
!    Exchange-correlation energy
     exc(ipt)=n1*d1m1
!    
!    Exchange-correlation potential
     vxc(ipt)=n2*d1m1**2
!    Assemble derivative of vxc wrt rs
     n3=-(c1+rs*(2._dp*c2+rs*(3._dp*c3+rs*(4._dp*c4+rs*(5._dp*c5+&
&     rs*(6._dp*c6+rs*(7._dp*c7)))))))
     d3=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
     dvxcdr=(n3-d3*(2._dp*n2*d1m1))*d1m1**2
     rhom1=1.0_dp/rhor(ipt)
!    
!    derivative of vxc wrt rho
     dvxc(ipt)=-dvxcdr*rs*third*rhom1
!    
   end do
 end if
end subroutine xctetr
!!***

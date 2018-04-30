!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcwign
!! NAME
!! xcwign
!!
!! FUNCTION
!! Returns exc, vxc, and eventually d(vxc)/d($\rho$) from input $\rho$.
!! Wigner exchange and correlation (xc)--see e.g. David Pines,
!! Elementary Excitations in Solids, p. 94, NY 1964.
!! Expression is exc=-(0.44)/(rs+7.8)-efac/rs (hartree), efac below.
!! rs = $(3/(4\pi))^{1/3}* \rho (r)^{-1/3}$.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=corresponding Wigner-Seitz radii, precomputed
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d($\rho$*exc)/d($\rho$)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xcwign(exc,npt,order,rspts,vxc,& !Mandatory arguments
&                dvxc)                           !Optional arguments

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcwign'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!c1 and c2 are the Wigner parameters in hartree and bohr resp.
!scalars
 integer :: ipt
 real(dp),parameter :: c1=0.44_dp,c2=7.8_dp,c4_3=4.0_dp/3.0_dp
 real(dp),parameter :: c8_27=8.0_dp/27.0_dp
 real(dp) :: dfac,efac,rs,rsc2m1,rsm1,vfac,vxcnum
 character(len=500) :: message

! *************************************************************************

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,i0)' )&
&   'With Wigner xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be',order
   MSG_BUG(message)
 end if

!Checks the compatibility between the order and the presence of the optional arguments
 if(order <= 1 .and. present(dvxc))then
   write(message, '(a,a,a,i3)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector dvxc, that is needed only with order=2 , while we have',order
   MSG_BUG(message)
 end if

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac

!separate cases with respect to order
 if (order==2) then
   
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
     rsc2m1=1.0_dp/(rs+c2)
!    compute energy density (hartree)
     exc(ipt)=-c1*rsc2m1-efac*rsm1
     vxcnum=-(c4_3*rs+c2)*c1
!    compute potential (hartree)
     vxc(ipt)=vxcnum*rsc2m1**2-vfac*rsm1
!    compute d(vxc)/d(rho) (hartree*bohr^3)
     dvxc(ipt)=-(c8_27*pi)*(c1*rs**4)*(rs+rs+c2)*rsc2m1**3-dfac*rs**2
   end do
 else
   
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
     rsc2m1=1.0_dp/(rs+c2)
!    compute energy density (hartree)
     exc(ipt)=-c1*rsc2m1-efac*rsm1
     vxcnum=-(c4_3*rs+c2)*c1
!    compute potential (hartree)
     vxc(ipt)=vxcnum*rsc2m1**2-vfac*rsm1
   end do
   
 end if
!
end subroutine xcwign
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcxalp
!! NAME
!! xcxalp
!!
!! FUNCTION
!! Returns exc, vxc, and eventually d(vxc)/d($\rho$) from input $\rho$.
!! "X$\alpha$" method is used in this subroutine:
!! a single fixed value is chosen for "alpha", set below.
!! Expression is exc=-alpha*efac/rs (hartree), efac below.
!! rs = $(3/(4\pi))^{1/3}* \rho (r)^{-1/3}$.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=Wigner-Seitz radii, at each point
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


subroutine xcxalp(exc,npt,order,rspts,vxc, dvxc)  ! dvxc is optional

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcxalp'
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
!Set value of alpha in "X-alpha" method
!scalars
 integer :: ipt
 real(dp),parameter :: alpha=1.0_dp
 real(dp) :: dfac,efac,rs,rsm1,vfac
 character(len=500) :: message

! *************************************************************************

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,i3)' )&
&   'With X-alpha xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be',order
   MSG_BUG(message)
 end if

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac

!separate cases with respect to order
 if(order==2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
!    compute energy density (hartree)
     exc(ipt)=-alpha*efac*rsm1
!    compute potential (hartree)
     vxc(ipt)=-alpha*vfac*rsm1
!    compute d(vxc)/d(rho) (hartree*bohr^3)
     dvxc(ipt)=-alpha*dfac*rs**2
   end do
 else
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
!    compute energy density (hartree)
     exc(ipt)=-alpha*efac*rsm1
!    compute potential (hartree)
     vxc(ipt)=-alpha*vfac*rsm1
   end do
 end if
!
end subroutine xcxalp
!!***

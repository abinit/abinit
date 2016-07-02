!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcpzca
!! NAME
!! xcpzca
!!
!! FUNCTION
!! Returns exc, vxc, and d(vxc)/d($\rho$) from input rho.
!!
!! NOTE
!! Perdew-Zunger parameterization of Ceperly-Alder electron gas
!! energy data--
!! J. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981).
!! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980).
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


subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
&                dvxc)                            !Optional arguments

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcpzca'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rhor(npt),rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!Perdew-Zunger parameters a, b, b1, b2, c, d, gamma
!scalars
 integer :: ipt
 real(dp),parameter :: aa=0.0311_dp,b1=1.0529_dp,b2=0.3334_dp,bb=-0.048_dp
 real(dp),parameter :: c4_3=4.0_dp/3.0_dp,c7_6=7.0_dp/6.0_dp,cc=0.0020_dp
 real(dp),parameter :: dd=-0.0116_dp,ga=-0.1423_dp
 real(dp) :: den,den3,dfac,efac,logrs,rs,rsm1,t1,t2,vfac
 character(len=500) :: message

! *************************************************************************

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,i0)' )&
&   'With Perdew-Zunger Ceperley-Alder xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be',order
   MSG_BUG(message)
 end if

!Checks the compatibility between the order and the presence of the optional arguments
 if(order <= 1 .and. present(dvxc))then
   write(message, '(a,a,a,i0)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector dvxc, that is needed only with order=2 , while we have',order
   MSG_BUG(message)
 end if

!separate cases with respect to order
 if(order==2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
!    Consider two regimes: rs<1 or rs>=1
     if (rs<1._dp) then
       logrs=log(rs)
!      compute energy density exc (hartree)
       exc(ipt)=(aa+cc*rs)*logrs+dd*rs+bb-efac*rsm1
!      compute potential vxc=d(rho*exc)/d(rho) (hartree)
       vxc(ipt)=(aa+two_thirds*cc*rs)*logrs+(dd+dd-cc)*rs*third+&
&       (bb-aa*third)-vfac*rsm1
!      compute d(vxc)/d(rho) (hartree*bohr^3)
       dvxc(ipt)=-(3._dp*aa+(cc+dd+dd)*rs+2._dp*cc*rs*logrs)&
&       /(9._dp*rhor(ipt))-dfac*rs**2
     else if (rs<1000._dp) then
       t1=b1*sqrt(rs)
       t2=b2*rs
       den=1._dp/(1._dp+t1+t2)
       exc(ipt)=ga*den-efac*rsm1
       vxc(ipt)=ga*(1._dp+c7_6*t1+c4_3*t2)*den**2-vfac*rsm1
       den3=den**3
       dvxc(ipt)=(ga*den3/(36._dp*rhor(ipt)))*(5._dp*t1+8._dp*t2+&
&       7._dp*t1**2+16._dp*t2**2+21._dp*t1*t2)-dfac*rs**2
     else
       t1=b1*sqrt(rs)
       t2=b2*rs
       den=1._dp/(1._dp+t1+t2)
       exc(ipt)=ga*den-efac*rsm1
       vxc(ipt)=ga*(1._dp+c7_6*t1+c4_3*t2)*den**2-vfac*rsm1
       dvxc(ipt)=0._dp
     end if
   end do
 else
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
!    Consider two regimes: rs<1 or rs>=1
     if (rs<1._dp) then
       logrs=log(rs)
!      compute energy density exc (hartree)
       exc(ipt)=(aa+cc*rs)*logrs+dd*rs+bb-efac*rsm1
!      compute potential vxc=d(rho*exc)/d(rho) (hartree)
       vxc(ipt)=(aa+two_thirds*cc*rs)*logrs+(dd+dd-cc)*rs*third+&
&       (bb-aa*third)-vfac*rsm1
!      compute d(vxc)/d(rho) (hartree*bohr^3)
     else
       t1=b1*sqrt(rs)
       t2=b2*rs
       den=1._dp/(1._dp+t1+t2)
       exc(ipt)=ga*den-efac*rsm1
       vxc(ipt)=ga*(1._dp+c7_6*t1+c4_3*t2)*den**2-vfac*rsm1
     end if
   end do
 end if
!
end subroutine xcpzca
!!***

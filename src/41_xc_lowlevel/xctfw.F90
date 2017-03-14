!{\src2tex{textfont=tt}}
!!****f* ABINIT/xctfw
!! NAME
!! xctfw
!!
!! FUNCTION
!! Add gradient part of the Thomas-Fermi-Weizsacker functional
!! Perrot F.,Phys. Rev. A20,586-594 (1979)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (JFD,LK)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  ndvxcdgr= size of dvxcdgr(npts,ndvxcdgr)
!!  ngr2= size of grho2_updn(npts,ngr2)
!!  npts= number of points to be computed
!!  nspden=number if spin density component (necessarily 1 here)
!!  grho2_updn(npts,ngr2)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2),
!!     only used if gradient corrected functional (option=2,-2,-4 and 4 or beyond)
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!  temp= electronic temperature
!!  usefxc=1 if free energy fxc is used
!!
!! SIDE EFFECTS
!!  The following arrays are modified (gradient correction added):
!!  dvxcdgr(npts,3)=partial derivative of the XC energy divided by the norm of the gradient
!!  exci(npts)=exchange-correlation energy density
!!  fxci(npts)=free energy energy density
!!  vxci(npts,nspden)=exchange-correlation potential
!!
!! PARENTS
!!      rhohxc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE
!!$#if defined HAVE_CONFIG_H
!!$#include "config.h"
!!$#endif

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xctfw(temp,exci,fxci,usefxc,rho_updn,vxci,npts,nspden,dvxcdgr,ndvxcdgr,grho2_updn,ngr2)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xctfw'
 use interfaces_41_xc_lowlevel, except_this_one => xctfw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndvxcdgr,ngr2,npts,nspden,usefxc
 real(dp),intent(in) :: temp
!arrays
 real(dp),intent(in) :: grho2_updn(npts,ngr2),rho_updn(npts,nspden)
 real(dp),intent(inout) :: dvxcdgr(npts,ndvxcdgr),fxci(npts*usefxc),exci(npts),vxci(npts,nspden)

!Local variables-------------------------------
!scalars
 integer :: iperrot,ipts
 logical :: has_fxc,has_dvxcdgr
 real(dp) :: etfw,rho,rho_inv,rhomot,yperrot0,vtfw
 real(dp) :: yperrot,uperrot,dyperrotdn,duperrotdyperrot
 real(dp) :: hperrot,dhperrotdyperrot,dhperrotdn,dhperrotduperrot
!arrays
 real(dp) :: wpy(0:7), wpu(0:7)
 real(dp),allocatable :: rho_updnm1_3(:,:)

! *************************************************************************

!DBG_ENTER('COLL')

 has_fxc=(usefxc/=0)
 has_dvxcdgr=(ndvxcdgr/=0)

 yperrot0=1.666081101_dp

 wpy(0)=0.5_dp; wpy(1)=-0.1999176316_dp
 wpy(2)=0.09765615709_dp; wpy(3)=-0.06237609924_dp
 wpy(4)=0.05801466322_dp; wpy(5)=-0.04449287774_dp
 wpy(6)=0.01903211697_dp; wpy(7)=-0.003284096926_dp

 wpu(0)=one/6._dp; wpu(1)=0.311590799_dp
 wpu(2)=3.295662439_dp; wpu(3)=-29.22038326_dp
 wpu(4)=116.1084531_dp; wpu(5)=-250.4543147_dp
 wpu(6)=281.433688_dp; wpu(7)=-128.8784806_dp

 ABI_ALLOCATE(rho_updnm1_3,(npts,2))

 call invcb(rho_updn(:,1),rho_updnm1_3(:,1),npts)

 do ipts=1,npts

   rho   =rho_updn(ipts,1)
   rhomot=rho_updnm1_3(ipts,1)
   rho_inv=rhomot*rhomot*rhomot

   yperrot=pi*pi/sqrt2/temp**1.5*two*rho
   uperrot=yperrot**(2./3.)

   dyperrotdn=pi*pi/sqrt2/temp**1.5*2.0_dp

   hperrot=zero
   dhperrotdyperrot=zero
   dhperrotduperrot=zero
   if(yperrot<=yperrot0)then
     do iperrot=0,7
       hperrot=hperrot+wpy(iperrot)*yperrot**iperrot
       dhperrotdyperrot=dhperrotdyperrot+iperrot*wpy(iperrot)*yperrot**(iperrot-1)
     end do
     hperrot=one/12.0_dp*hperrot
     dhperrotdyperrot=one/12.0_dp*dhperrotdyperrot
     dhperrotdn=dhperrotdyperrot*dyperrotdn
   else
     do iperrot=0,7
       hperrot=hperrot+wpu(iperrot)/uperrot**(2*iperrot)
       dhperrotduperrot=dhperrotduperrot-2.*iperrot*wpu(iperrot)/uperrot**(2*iperrot+1)
     end do
     hperrot=one/12.0_dp*hperrot
     dhperrotduperrot=one/12.0_dp*dhperrotduperrot
     duperrotdyperrot=two/3._dp/yperrot**(1./3.)
     dhperrotdn=dhperrotduperrot*duperrotdyperrot*dyperrotdn
   end if

   etfw=hperrot*grho2_updn(ipts,1)*rho_inv*rho_inv
   vtfw=-etfw + rho/hperrot*dhperrotdn*etfw

   if(yperrot<=yperrot0)then
     exci(ipts)   = exci(ipts) + etfw + 1.5_dp*yperrot*dhperrotdyperrot*grho2_updn(ipts,1)*rho_inv*rho_inv
   else
     exci(ipts)   = exci(ipts) + etfw + uperrot*dhperrotduperrot*grho2_updn(ipts,1)*rho_inv*rho_inv
   end if
   vxci(ipts,1) = vxci(ipts,1)  + vtfw
   if (has_fxc) fxci(ipts)   = fxci(ipts) + etfw
   if (has_dvxcdgr) dvxcdgr(ipts,1)= dvxcdgr(ipts,1)+two*hperrot*rho_inv

 end do

 ABI_DEALLOCATE(rho_updnm1_3)

!DBG_EXIT('COLL')

end subroutine xctfw
!!***

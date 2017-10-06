!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcmult
!! NAME
!! xcmult
!!
!! FUNCTION
!! In the case of GGA, multiply the different gradient of spin-density
!! by the derivative of the XC functional with respect
!! to the norm of the gradient, then divide it by the
!! norm of the gradient
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  depsxc(nfft,nspgrad)=derivative of Exc with respect to the (spin-)density,
!!    or to the norm of the gradient of the (spin-)density,
!!    further divided by the norm of the gradient of the (spin-)density
!!   The different components of depsxc will be
!!   for nspden=1,         depsxc(:,1)=d(rho.exc)/d(rho)
!!         and if ngrad=2, depsxc(:,2)=1/2*1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!!                                      +   1/|grad rho|*d(rho.exc)/d(|grad rho|)
!!         (do not forget : |grad rho| /= |grad rho_up| + |grad rho_down|
!!   for nspden=2,         depsxc(:,1)=d(rho.exc)/d(rho_up)
!!                         depsxc(:,2)=d(rho.exc)/d(rho_down)
!!         and if ngrad=2, depsxc(:,3)=1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!!                         depsxc(:,4)=1/|grad rho_down|*d(rho.exc)/d(|grad rho_down|)
!!                         depsxc(:,5)=1/|grad rho|*d(rho.exc)/d(|grad rho|)
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngrad = must be 2
!!  nspden=number of spin-density components
!!  nspgrad=number of spin-density and spin-density-gradient components
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  rhonow(nfft,nspden,ngrad*ngrad)=
!!   at input :
!!    electron (spin)-density in real space and its gradient,
!!    either on the unshifted grid (if ishift==0,
!!      then equal to rhor), or on the shifted grid
!!     rhonow(:,:,1)=electron density in electrons/bohr**3
!!     rhonow(:,:,2:4)=gradient of electron density in el./bohr**4
!!   at output :
!!    rhonow(:,:,2:4) has been multiplied by the proper factor,
!!    described above.
!!
!! PARENTS
!!      m_pawxc,rhotoxc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xcmult (depsxc,nfft,ngrad,nspden,nspgrad,rhonow)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcmult'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ngrad,nspden,nspgrad
!arrays
 real(dp),intent(in) :: depsxc(nfft,nspgrad)
 real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)

!Local variables-------------------------------
!scalars
 integer :: idir,ifft
 real(dp) :: rho_tot,rho_up

! *************************************************************************

 do idir=1,3

   if(nspden==1)then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(depsxc,idir,nfft,rhonow)
     do ifft=1,nfft
       rhonow(ifft,1,1+idir)=rhonow(ifft,1,1+idir)*depsxc(ifft,2)
     end do

   else

!    In the spin-polarized case, there are more factors to take into account
!$OMP PARALLEL DO PRIVATE(ifft,rho_tot,rho_up) SHARED(depsxc,idir,nfft,rhonow)
     do ifft=1,nfft
       rho_tot=rhonow(ifft,1,1+idir)
       rho_up =rhonow(ifft,2,1+idir)
       rhonow(ifft,1,1+idir)=rho_up *depsxc(ifft,3)         + rho_tot*depsxc(ifft,5)
       rhonow(ifft,2,1+idir)=(rho_tot-rho_up)*depsxc(ifft,4)+ rho_tot*depsxc(ifft,5)
     end do

   end if ! nspden==1

 end do ! End loop on directions

end subroutine xcmult
!!***

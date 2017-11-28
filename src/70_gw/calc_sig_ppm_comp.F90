!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sig_ppm_comp
!!
!! NAME
!! calc_sig_ppm_comp
!!
!! FUNCTION
!! Calculating contributions to self-energy operator using a plasmon-pole model
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (FB, GMR, VO, LR, RWG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nomega=number of frequencies to consider
!!  npwc= number of G vectors in the plasmon pole
!!  npwc1= 1 if ppmodel==3, =npwc if ppmodel== 4, 1 for all the other cases
!!  npwc2= 1 if ppmodel==3, =1    if ppmodel== 4, 1 for all the other cases
!!  npwx=number of G vectors in rhotwgp
!!  ppmodel=plasmon pole model
!!  theta_mu_minus_e0i= $\theta(\mu-\epsilon_{k-q,b1,s}), defines if the state is occupied or not
!!  zcut=small imaginary part to avoid the divergence. (see related input variable)
!!  omegame0i(nomega)=frequencies where evaluate \Sigma_c ($\omega$ - $\epsilon_i$
!!  otq(npwc,npwc2)=plasmon pole parameters for this q-point
!!  botsq(npwc,npwc1)=plasmon pole parameters for this q-point
!!  eig(npwc,npwc)=the eigvectors of the symmetrized inverse dielectric matrix for this q point
!!   (first index for G, second index for bands)
!!  rhotwgp(npwx)=oscillator matrix elements divided by |q+G| i.e
!!    $\frac{\langle b1 k-q s | e^{-i(q+G)r | b2 k s \rangle}{|q+G|}$
!!
!! OUTPUT
!!  sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!
!!  ket(npwc,nomega):
!!
!!  In case of ppmodel==1,2 it contains
!!
!!   ket(G,omega) = Sum_G2       conjg(rhotw(G)) * Omega(G,G2) * rhotw(G2)
!!                          ---------------------------------------------------
!!                            2 omegatw(G,G2) (omega-E_i + omegatw(G,G2)(2f-1))
!!
!! NOTES
!! Taken from old routine
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_sig_ppm_comp(npwc,nomega,rhotwgp,botsq,otq,omegame0i_io,zcut,theta_mu_minus_e0i,ket,ppmodel,npwx,npwc1,npwc2)

 use m_profiling_abi
 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sig_ppm_comp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwc1,npwc2,npwx,ppmodel
 real(dp),intent(in) :: omegame0i_io,theta_mu_minus_e0i,zcut
!arrays
 complex(gwpc),intent(in) :: botsq(npwc,npwc1),rhotwgp(npwx),otq(npwc,npwc2)
 complex(gwpc),intent(inout) :: ket(npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,io
 real(dp) :: den,otw,twofm1_zcut
 complex(gwpc) :: num,rhotwgdp_igp
 logical :: fully_occupied,totally_empty
 character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: ket_comp(:)

!*************************************************************************

 if (ppmodel/=1.and.ppmodel/=2) then
   write(msg,'(a,i0,a)')' The completeness trick cannot be used when ppmodel is ',ppmodel,' It should be set to 1 or 2. '
   MSG_ERROR(msg)
 end if

 ABI_ALLOCATE(ket_comp,(npwc))
 ket_comp(:)=0.d0

 fully_occupied=(abs(theta_mu_minus_e0i-1.)<0.001)
 totally_empty=(abs(theta_mu_minus_e0i)<0.001)

 if(.not.(totally_empty)) then ! not totally empty
   twofm1_zcut=zcut
   do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
       otw=DBLE(otq(ig,igp)) ! in principle otw -> otw - ieta
       num = botsq(ig,igp)*rhotwgdp_igp

       den = omegame0i_io-otw
       if (den**2>zcut**2) then
         ket_comp(ig) = ket_comp(ig) - num/(den*otw)*theta_mu_minus_e0i
       end if
     end do !ig
   end do !igp
 end if ! not totally empty

 if(.not.(fully_occupied)) then ! not fully occupied
   twofm1_zcut=-zcut

   do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
       otw=DBLE(otq(ig,igp)) ! in principle otw -> otw - ieta
       num = botsq(ig,igp)*rhotwgdp_igp

       den = omegame0i_io-otw
       if (den**2>zcut**2) then
         ket_comp(ig) = ket_comp(ig) - num/(den*otw)*(1.-theta_mu_minus_e0i)
       end if
     end do !ig
   end do !igp
 end if ! not fully occupied

 do io=1,nomega
   ket(:,io)=ket(:,io)+0.5*ket_comp(:)
 end do

 ABI_DEALLOCATE(ket_comp)

end subroutine calc_sig_ppm_comp
!!***

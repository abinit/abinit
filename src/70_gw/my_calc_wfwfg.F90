!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_coh
!! NAME
!! calc_coh
!!
!! FUNCTION
!!  Calculates the partial contribution to the COH part of the COHSEX self-energy for a given q-point.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (FB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! iqibz=index of the irreducible q-point in the array qibz, point which is
!!  related by a symmetry operation to the point q summed over (see csigme).
!!  This index is also used to treat the integrable coulombian singularity at q=0
!! ngfft(18)=contain all needed information about 3D FFT for GW wavefuntions,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nsig_ab=Number of components in the self-energy operator (1 for collinear magnetism)
!! npwc=number of plane waves in $\tilde epsilon^{-1}$
!! nspinor=Number of spinorial components.
!! nfftot=number of points in real space
!! i_sz=contribution arising from the integrable coulombian singularity at q==0
!! (see csigme for the method used), note that in case of 3-D systems the factor
!! 4pi in the coulombian potential is included in the definition of i_sz
!! gvec(3,npwc)=G vectors in reduced coordinates
!! vc_sqrt(npwc)= square root of the coulombian matrix elements for this q-point
!! epsm1q_o(npwc,npwc)= contains $\tilde epsilon^{-1}(q,w=0) - \delta_{G Gp}$ for
!!  the particular q-point considered in the sum
!! wfg2_jk(nsig_ab*nfftot)= Fourier Transform of $\u_{jb k}^*(r) u_{kb k}$
!!  jb,kb=left and righ band indeces definining the left and right states where the
!!  partial contribution to the matrix element of $\Sigma_{COH}$ is evaluated
!!
!! OUTPUT
!! sigcohme=partial contribution to the matrix element of $<jb k \sigma|\Sigma_{COH} | kb k \sigma>$
!!  coming from this single q-point
!!
!! PARENTS
!!      cohsex_me
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_coh(nspinor,nsig_ab,nfftot,ngfft,npwc,gvec,wfg2_jk,epsm1q_o,vc_sqrt,i_sz,iqibz,same_band,sigcohme)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_gwdefs, only : czero_gw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_coh'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nfftot,npwc,nsig_ab,nspinor
 real(dp),intent(in) :: i_sz
 logical,intent(in) :: same_band
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfft(18)
 complex(gwpc),intent(in) :: epsm1q_o(npwc,npwc),vc_sqrt(npwc)
 complex(gwpc),intent(in) :: wfg2_jk(nfftot*nsig_ab)
 complex(gwpc),intent(out) :: sigcohme(nsig_ab)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ig,ig4,ig4x,ig4y,ig4z,igp,igmin,ispinor
 integer :: spad,outofbox
 character(len=500) :: msg
!arrays
 integer :: g2mg1(3)

! *************************************************************************

 DBG_ENTER("COLL")

 ! === Partial contribution to the matrix element of Sigma_c ===
 ! * For nspinor==2, the closure relation reads:
 !  $\sum_s \psi_a^*(1)\psi_b(2) = \delta_{ab} \delta(1-2)$
 !  where a,b are the spinor components. As a consequence, Sigma_{COH} is always
 !  diagonal in spin-space and only diagonal matrix elements have to be calculated.
 ! MG  TODO wfg2_jk should be calculated on an augmented FFT box to avoid spurious wrapping of G1-G2.
 ! MG: One has to make sure G1-G2 is still in the FFT mesh for each G1 and G2 in chi0 (not always true)
 ! MODULO wraps G1-G2 in the FFT box but the Fourier components are not periodic!

 ! * Treat the case q --> 0 adequately.
 ! TODO Better treatment of wings, check cutoff in the coulombian interaction.
 igmin=1; if (iqibz==1) igmin=2

 sigcohme(:)=czero_gw

 do ispinor=1,nspinor
   spad=(ispinor-1)*nfftot
   outofbox=0

   do igp=igmin,npwc
     do ig=igmin,npwc

      g2mg1 = gvec(:,igp)-gvec(:,ig)
      if (ANY(g2mg1(:)>ngfft(1:3)/2) .or. ANY(g2mg1(:)<-(ngfft(1:3)-1)/2)) then
        outofbox = outofbox+1; CYCLE
      end if

      ig4x=MODULO(g2mg1(1),ngfft(1))
      ig4y=MODULO(g2mg1(2),ngfft(2))
      ig4z=MODULO(g2mg1(3),ngfft(3))
      ig4= 1+ig4x+ig4y*ngfft(1)+ig4z*ngfft(1)*ngfft(2)

      sigcohme(ispinor) = sigcohme(ispinor) + &
&                         half*wfg2_jk(spad+ig4)*epsm1q_o(ig,igp)*vc_sqrt(ig)*vc_sqrt(igp)
     end do !ig
   end do !igp

   if (iqibz==1.and.same_band) then
     sigcohme(ispinor) = sigcohme(ispinor) + half*wfg2_jk(spad+ig4)*epsm1q_o(1,1)*i_sz
   end if
 end do !ispinor

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=50) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==50) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine calc_coh
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/calc_coh_comp
!! NAME
!! calc_coh_comp
!!
!! FUNCTION
!!  Calculates the COH-like contribution to the self-energy when
!!  the extrapolar technique and the closure relation is used to
!!  reduce the number of empty states to be summed over in the Green
!!  function entering the definition of the GW self-energy.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (FB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! iqibz=index of the irreducible q-point in the array qibz, point which is
!!  related by a symmetry operation to the point q summed over (see csigme).
!!  This index is also used to treat the integrable coulombian singularity at q=0
!! ngfft(18)=contain all needed information about 3D FFT for GW wavefuntions,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nsig_ab=Number of components in the self-energy operator (1 for collinear magnetism)
!! npwc=number of plane waves in $\tilde epsilon^{-1}$
!! nspinor=Number of spinorial components.
!! i_sz=contribution arising from the integrable coulombian singularity at q==0
!! (see csigme for the method used), note that in case of 3-D systems the factor
!! 4pi in the coulombian potential is included in the definition of i_sz
!! gvec(3,npwc)=G vectors in reduced coordinates
!! vc_sqrt(npwc)= square root of the coulombian matrix elements for this q-point
!! botsq = Plasmon-pole parameters
!! otq  = PPm parameters
!!
!! OUTPUT
!! sigcohme=partial contribution to the matrix element of $<jb k|\Sigma_{COH}| kb k>$
!!  coming from this single q-point for completeness trick
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_coh_comp(iqibz,i_sz,same_band,nspinor,nsig_ab,ediff,npwc,gvec,&
&  ngfft,nfftot,wfg2_jk,vc_sqrt,botsq,otq,sigcohme)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_gwdefs, only : czero_gw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_coh_comp'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,npwc,nsig_ab,nspinor,nfftot
 real(dp),intent(in) :: i_sz,ediff
 logical,intent(in) :: same_band
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfft(18)
 complex(gwpc),intent(in) :: botsq(npwc,npwc),otq(npwc,npwc)
 complex(gwpc),intent(in) :: vc_sqrt(npwc)
 complex(gwpc),intent(in) :: wfg2_jk(nfftot*nsig_ab)
 complex(gwpc),intent(out) :: sigcohme(nsig_ab)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ig,ig4,ig4x,ig4y,ig4z,igp,igmin,ispinor,ngfft1,ngfft2,ngfft3
 integer :: spad,outofbox
 character(len=500) :: msg
!arrays
 integer :: g2mg1(3)

! *************************************************************************

 DBG_ENTER("COLL")

 ! === Treat the case q --> 0 adequately ===
 ! TODO Better treatment of wings
 igmin=1 ; if (iqibz==1) igmin=2
 !
 ! === Partial contribution to the matrix element of Sigma_c ===
 ! * For nspinor==2, the closure relation reads:
 !  $\sum_s \psi_a^*(1)\psi_b(2) = \delta_{ab} \delta(1-2)$
 !  where a,b are the spinor components. As a consequence, Sigma_{COH} is always
 !  diagonal in spin-space and only diagonal matrix elements have to be calculated.
 ! MG  TODO wfg2_jk should be calculated on an augmented FFT box to avoid spurious wrapping of G1-G2.
 !
 ngfft1 =ngfft(1)
 ngfft2 =ngfft(2)
 ngfft3 =ngfft(3)
 sigcohme(:)=czero_gw

 do ispinor=1,nspinor
  spad=(ispinor-1)*nfftot
  outofbox=0

   do igp=igmin,npwc
     do ig=igmin,npwc

      g2mg1 = gvec(:,igp)-gvec(:,ig)
      if (ANY(g2mg1(:)>ngfft(1:3)/2) .or. ANY(g2mg1(:)<-(ngfft(1:3)-1)/2)) then
        outofbox = outofbox+1; CYCLE
      end if

      ig4x=MODULO(g2mg1(1),ngfft1)
      ig4y=MODULO(g2mg1(2),ngfft2)
      ig4z=MODULO(g2mg1(3),ngfft3)
      ig4= 1+ig4x+ig4y*ngfft1+ig4z*ngfft1*ngfft2

      !MG where is neta here, ediff, otq might be close to zero depending on gwecomp
      sigcohme(ispinor) = sigcohme(ispinor) + &
&       half*wfg2_jk(spad+ig4)*vc_sqrt(ig)*vc_sqrt(igp) * botsq(ig,igp) / ( otq(ig,igp) * ( ediff -otq(ig,igp) ) )
     end do
   end do

   if (iqibz==1.and.same_band) then
     sigcohme(ispinor) = sigcohme(ispinor) + half*wfg2_jk(spad+ig4)*i_sz*botsq(1,1) / ( otq(1,1) * (ediff -otq(1,1)) )
   end if
 end do !ispinor

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=50) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==50) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine calc_coh_comp
!!***

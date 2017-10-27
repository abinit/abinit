!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2kindstr2
!! NAME
!! d2kindstr2
!!
!! FUNCTION
!! compute expectation value of the second derivatives of the kinetic energy
!! wrt strain for one band and kpoint
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cwavef(2,npw*nspinor)=wavefunction for current band
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  gprimd(3,3)=primitive vectors in reciprocal space
!!  istwfk=information about wavefunction storage
!!  kg_k(3,npw)=integer coordinates of planewaves in basis sphere.
!!  kpt(3)=reduced coordinates of k point
!!  npw=number of plane waves at kpt.
!!  nspinor=number of spinorial components of the wavefunction
!!
!! OUTPUT
!!  ekinout(36)=expectation values of the second strain derivatives
!!   of the (modified) kinetic energy
!!
!! NOTES
!! Usually, the kinetic energy expression is $(1/2) (2 \pi)^2 (k+G)^2 $
!! However, the present implementation allows for a modification
!! of this kinetic energy, in order to obtain smooth total energy
!! curves with respect to the cut-off energy or the cell size and shape.
!! Thus the usual expression is kept if it is lower then ecut-ecutsm,
!! zero is returned beyond ecut, and in between, the kinetic
!! energy is DIVIDED by a smearing factor (to make it infinite at the
!! cut-off energy). The smearing factor is $x^2 (3-2x)$, where
!! x = (ecut- unmodified energy)/ecutsm.
!!
!! PARENTS
!!      dfpt_eltfrkin
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine d2kindstr2(cwavef,ecut,ecutsm,effmass_free,ekinout,gmet,gprimd,&
&            istwfk,kg_k,kpt,npw,nspinor)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2kindstr2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,npw,nspinor
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: cwavef(2,npw*nspinor),gmet(3,3),gprimd(3,3),kpt(3)
 real(dp),intent(inout) :: ekinout(36) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: ig,igs,ii,ispinor,istr1,istr2,ka,kb,kd,kg
 real(dp) :: d2fkin,d2fsm,d2kinacc,d2kpg2,dfkin,dfsm,dkpg21,dkpg22,ecutsm_inv
 real(dp) :: fsm,gpk1,gpk2,gpk3,htpisq,kpg2,term,xx
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: d2gm(3,3),dgm01(3,3),dgm10(3,3)

! *************************************************************************
!
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

!Loop over 2nd strain index
 do istr2=1,6
!  Loop over 1st strain index, upper triangle only
   do istr1=1,istr2

     ka=idx(2*istr1-1);kb=idx(2*istr1);kg=idx(2*istr2-1);kd=idx(2*istr2)

     do ii = 1,3
       dgm01(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
       dgm10(:,ii)=-(gprimd(kg,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kg,ii))
     end do

     d2gm(:,:)=0._dp
     do ii = 1,3
       if(ka==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kb,ii)
       if(ka==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(kb,ii)
       if(kb==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(ka,ii)
       if(kb==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(ka,ii)
     end do
     d2gm(:,:)=0.5_dp*d2gm(:,:)

     d2kinacc=0._dp

!    loop on spinor index
     do ispinor=1,nspinor
       igs=(ispinor-1)*npw
!      loop on plane waves
       do ig=1,npw
         gpk1=dble(kg_k(1,ig))+kpt(1)
         gpk2=dble(kg_k(2,ig))+kpt(2)
         gpk3=dble(kg_k(3,ig))+kpt(3)
         kpg2=htpisq*&
&         ( gmet(1,1)*gpk1**2+         &
&         gmet(2,2)*gpk2**2+         &
&         gmet(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*gmet(1,2)*gpk2+  &
&         gpk1*gmet(1,3)*gpk3+  &
&         gpk2*gmet(2,3)*gpk3 )  )
         dkpg21=htpisq*&
&         ( dgm01(1,1)*gpk1**2+         &
&         dgm01(2,2)*gpk2**2+         &
&         dgm01(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*dgm01(1,2)*gpk2+  &
&         gpk1*dgm01(1,3)*gpk3+  &
&         gpk2*dgm01(2,3)*gpk3 )  )
         dkpg22=htpisq*&
&         ( dgm10(1,1)*gpk1**2+         &
&         dgm10(2,2)*gpk2**2+         &
&         dgm10(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*dgm10(1,2)*gpk2+  &
&         gpk1*dgm10(1,3)*gpk3+  &
&         gpk2*dgm10(2,3)*gpk3 )  )
         d2kpg2=htpisq*&
&         ( d2gm(1,1)*gpk1**2+         &
&         d2gm(2,2)*gpk2**2+         &
&         d2gm(3,3)*gpk3**2          &
&         +2.0_dp*(gpk1*d2gm(1,2)*gpk2+  &
&         gpk1*d2gm(1,3)*gpk3+  &
&         gpk2*d2gm(2,3)*gpk3 )  )

         if(kpg2>ecut-tol12)then
           dfkin=0._dp
           d2fkin=0._dp
         elseif(kpg2>ecut-ecutsm)then
!          This kinetic cutoff smoothing function and its xx derivatives
!          were produced with Mathematica and the fortran code has been
!          numerically checked against Mathematica.
           xx=(ecut-kpg2)*ecutsm_inv
           fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
           dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
           d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
&           (-144+45*xx))))))*fsm**3
           dfkin=fsm-ecutsm_inv*kpg2*dfsm
           d2fkin=ecutsm_inv*(-2.0_dp*dfsm+ecutsm_inv*kpg2*d2fsm)
         else
           dfkin=1._dp
           d2fkin=0._dp
         end if

!        accumulate kinetic energy 2nd derivative with wavefunction components
         term=d2fkin*dkpg21*dkpg22 + dfkin*d2kpg2
         if(istwfk==2 .and. ig/=1)term=2.0_dp*term
         if(istwfk>2)term=2.0_dp*term
         d2kinacc=d2kinacc + term*(cwavef(re,ig+igs)**2 + cwavef(im,ig+igs)**2)

       end do  !ig
     end do !ispinor

     ekinout(istr1+6*(istr2-1))=d2kinacc/effmass_free

   end do !istr1
 end do !istr2

end subroutine d2kindstr2
!!***

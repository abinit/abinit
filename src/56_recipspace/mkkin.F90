!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkkin
!! NAME
!! mkkin
!!
!! FUNCTION
!! compute elements of kinetic energy operator in reciprocal
!! space at a given k point
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  idir1 = 1st direction of the derivative (if 1 <= idir1 <= 3, not used otherwise)
!!  idir2 = 2st direction of the derivative (if 1 <= idir1,idir2 <= 3, not used otherwise))
!!  kg(3,npw)=integer coordinates of planewaves in basis sphere.
!!  kpt(3)=reduced coordinates of k point
!!  npw=number of plane waves at kpt.
!!
!! OUTPUT
!!  kinpw(npw)=(modified) kinetic energy (or derivative) for each plane wave (Hartree)
!!
!! SIDE EFFECTS
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
!! This smearing factor is also used to derived a modified kinetic
!! contribution to stress, in another routine (forstrnps.f)
!!
!! Also, in order to break slightly the symmetry between axes, that causes
!! sometimes a degeneracy of eigenvalues and do not allow to obtain
!! the same results on different machines, there is a modification
!! by one part over 1.0e12 of the metric tensor elements (1,1) and (3,3)
!!
!! PARENTS
!!      calc_vhxc_me,d2frnl,dfpt_nsteltwf,dfpt_nstpaw,dfpt_rhofermi,energy
!!      getgh1c,ks_ddiago,m_io_kss,m_vkbr,mkffnl,vtorho
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkkin (ecut,ecutsm,effmass,gmet,kg,kinpw,kpt,npw,idir1,idir2)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkkin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 integer,intent(in) :: idir1,idir2
 real(dp),intent(in) :: ecut,ecutsm,effmass

!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),kpt(3)
 real(dp),intent(out) :: kinpw(npw)

!Local variables-------------------------------
!scalars
 integer :: ig,order
 real(dp),parameter :: break_symm=1.0d-11
 real(dp) :: ecutsm_inv,fsm,gpk1,gpk2,gpk3,htpisq,kinetic,kpg2,dkpg2,xx
 real(dp) :: d1kpg2,d2kpg2,ddfsm, dfsm
!arrays
 real(dp) :: gmet_break(3,3)

! *************************************************************************
!
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

 gmet_break(:,:)=gmet(:,:)
 gmet_break(1,1)=(1.0_dp+break_symm)*gmet(1,1)
 gmet_break(3,3)=(1.0_dp-break_symm)*gmet(3,3)

 order=0 ! Compute the kinetic operator
 if (idir1>0.and.idir1<4) then
   order=1 ! Compute the 1st derivative of the kinetic operator
   if (idir2>0.and.idir2<4) then
     order=2 ! Compute the 2nd derivative of the kinetic operator
   end if
 end if

!$OMP PARALLEL DO PRIVATE(dkpg2,d1kpg2,d2kpg2,gpk1,gpk2,gpk3,ig,kinetic,kpg2,xx,fsm,dfsm,ddfsm) &
!$OMP SHARED(kinpw,ecut,ecutsm,ecutsm_inv) &
!$OMP SHARED(gmet_break,htpisq,idir1,idir2,kg,kpt,npw)
 do ig=1,npw
   gpk1=dble(kg(1,ig))+kpt(1)
   gpk2=dble(kg(2,ig))+kpt(2)
   gpk3=dble(kg(3,ig))+kpt(3)
   kpg2=htpisq*&
&   ( gmet_break(1,1)*gpk1**2+         &
&   gmet_break(2,2)*gpk2**2+         &
&   gmet_break(3,3)*gpk3**2          &
&   +2.0_dp*(gpk1*gmet_break(1,2)*gpk2+  &
&   gpk1*gmet_break(1,3)*gpk3+  &
&   gpk2*gmet_break(2,3)*gpk3 )  )
   select case (order)
   case(0)
     kinetic=kpg2
   case(1)
     dkpg2=htpisq*2.0_dp*&
&     (gmet_break(idir1,1)*gpk1+gmet_break(idir1,2)*gpk2+gmet_break(idir1,3)*gpk3)
     kinetic=dkpg2
   case(2)
     dkpg2=htpisq*2.0_dp*gmet_break(idir1,idir2)
     kinetic=dkpg2
   end select
   if(kpg2>ecut-ecutsm)then
     if(kpg2>ecut-tol12)then
       if(order==0) then
!        Will filter the wavefunction, based on this value, in cgwf.f, getghc.f and precon.f
         kinetic=huge(0.0_dp)*1.d-10
       else
!        The wavefunction has been filtered : no derivative
         kinetic=0
       end if
     else
       if(order==0) then
         xx=max( (ecut-kpg2)*ecutsm_inv , 1.0d-20)
       else
         xx=(ecut-kpg2)*ecutsm_inv
       end if
       if(order==2) then
         d1kpg2=htpisq*2.0_dp*&
&         (gmet_break(idir1,1)*gpk1+gmet_break(idir1,2)*gpk2+gmet_break(idir1,3)*gpk3)
         d2kpg2=htpisq*2.0_dp*&
&         (gmet_break(idir2,1)*gpk1+gmet_break(idir2,2)*gpk2+gmet_break(idir2,3)*gpk3)
       end if
!      This kinetic cutoff smoothing function and its xx derivatives
!      were produced with Mathematica and the fortran code has been
!      numerically checked against Mathematica.
       fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
       if(order>0) dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
       if(order>1) ddfsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*(-144+45*xx))))))*fsm**3
       select case (order)
       case(0)
         kinetic=kpg2*fsm
       case(1)
         kinetic=dkpg2*(fsm-ecutsm_inv*kpg2*dfsm)
       case(2)
         kinetic=dkpg2*fsm&
&         -2.0_dp*d1kpg2*dfsm*ecutsm_inv*d2kpg2&
&         +kpg2*ddfsm*(ecutsm_inv**2)*d1kpg2*d2kpg2&
&         -kpg2*dfsm*ecutsm_inv*dkpg2
       end select
     end if
   end if
   kinpw(ig)=kinetic/effmass
 end do
!$OMP END PARALLEL DO

end subroutine mkkin
!!***

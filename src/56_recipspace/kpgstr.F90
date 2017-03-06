!{\src2tex{textfont=tt}}
!!****f* ABINIT/kpgstr
!! NAME
!! kpgstr
!!
!! FUNCTION
!! Compute elements of the derivative the kinetic energy operator in reciprocal
!! space at given k point wrt a single cartesian strain component
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (DRH, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  gmet(3,3) = reciprocal lattice metric tensor (Bohr**-2)
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  istr=1,...6 specifies cartesian strain component 11,22,33,32,31,21
!!  kg(3,npw) = integer coordinates of planewaves in basis sphere.
!!  kpt(3)    = reduced coordinates of k point
!!  npw       = number of plane waves at kpt.
!!
!! OUTPUT
!!  dkinpw(npw)=d/deps(istr) ( (1/2)*(2 pi)**2 * (k+G)**2 )
!!
!! NOTES
!!  Src_6response/kpg3.f
!!
!! PARENTS
!!      dfpt_nsteltwf,dfpt_nstpaw,dfpt_rhofermi,getgh1c
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine kpgstr(dkinpw,ecut,ecutsm,effmass,gmet,gprimd,istr,kg,kpt,npw)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpgstr'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: istr,npw
 real(dp),intent(in) :: ecut,ecutsm,effmass
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpt(3)
 real(dp),intent(out) :: dkinpw(npw)

!Local variables -------------------------
!scalars
 integer :: ig,ii,ka,kb
 real(dp) :: dfsm,dkinetic,dkpg2,ecutsm_inv,fsm,gpk1,gpk2,gpk3,htpisq
! real(dp) :: d2fsm ! used in commented section below
 real(dp) :: kpg2,xx
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: dgmetds(3,3)

! *********************************************************************

!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

!Compute derivative of metric tensor wrt strain component istr
 if(istr<1 .or. istr>6)then
   write(message, '(a,i10,a,a,a)' )&
&   'Input istr=',istr,' not allowed.',ch10,&
&   'Possible values are 1,2,3,4,5,6 only.'
   MSG_BUG(message)
 end if

 ka=idx(2*istr-1);kb=idx(2*istr)
 do ii = 1,3
   dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
 end do
!For historical reasons:
 dgmetds(:,:)=0.5_dp*dgmetds(:,:)

 do ig=1,npw
   gpk1=dble(kg(1,ig))+kpt(1)
   gpk2=dble(kg(2,ig))+kpt(2)
   gpk3=dble(kg(3,ig))+kpt(3)
   kpg2=htpisq*&
&   ( gmet(1,1)*gpk1**2+         &
&   gmet(2,2)*gpk2**2+         &
&   gmet(3,3)*gpk3**2          &
&   +2.0_dp*(gpk1*gmet(1,2)*gpk2+  &
&   gpk1*gmet(1,3)*gpk3+  &
&   gpk2*gmet(2,3)*gpk3 )  )
   dkpg2=htpisq*2.0_dp*&
&   (gpk1*(dgmetds(1,1)*gpk1+dgmetds(1,2)*gpk2+dgmetds(1,3)*gpk3)+  &
&   gpk2*(dgmetds(2,1)*gpk1+dgmetds(2,2)*gpk2+dgmetds(2,3)*gpk3)+  &
&   gpk3*(dgmetds(3,1)*gpk1+dgmetds(3,2)*gpk2+dgmetds(3,3)*gpk3) )
   dkinetic=dkpg2
   if(kpg2>ecut-ecutsm)then
     if(kpg2>ecut-tol12)then
!      The wavefunction has been filtered : no derivative
       dkinetic=0.0_dp
     else
       xx=(ecut-kpg2)*ecutsm_inv
!      This kinetic cutoff smoothing function and its xx derivatives
!      were produced with Mathematica and the fortran code has been
!      numerically checked against Mathematica.
       fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
       dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
!      d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
!      &                        (-144+45*xx))))))*fsm**3
       dkinetic=dkpg2*(fsm-ecutsm_inv*kpg2*dfsm)
     end if
   end if
   dkinpw(ig)=dkinetic/effmass
 end do

end subroutine kpgstr
!!***

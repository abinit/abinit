!!****f* ABINIT/size_dvxc
!! NAME
!! size_dvxc
!!
!! FUNCTION
!! Give the size of the array dvxc(npts,ndvxc) and the second dimension of the d2vxc(npts,nd2vxc)
!! needed for the allocations depending on the routine which is called from the drivexc routine
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc (DCA, XG, GMR, MF, GZ)
!!
!! INPUTS
!!  [add_tfw]= optional flag controling the addition of Weiszacker gradient correction to Thomas-Fermi XC energy
!!  ixc= choice of exchange-correlation scheme
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!
!! OUTPUT
!!  ndvxc size of the array dvxc(npts,ndvxc) for allocation
!!  ngr2 size of the array grho2_updn(npts,ngr2) for allocation
!!  nd2vxc size of the array d2vxc(npts,nd2vxc) for allocation
!!  nvxcdgr size of the array dvxcdgr(npts,nvxcdgr) for allocation
!!
!! PARENTS
!!      m_pawxc,rhohxc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order,add_tfw)

 use defs_basis
 use m_profiling_abi
 use libxc_functionals

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'size_dvxc'
!End of the abilint section

 implicit none

!Arguments----------------------
 integer, intent(in) :: ixc,nspden,order
 integer, intent(out) :: ndvxc,nd2vxc,ngr2,nvxcdgr
 logical, intent(in), optional :: add_tfw

!Local variables----------------
 logical :: add_tfw_

! *************************************************************************

 add_tfw_=.false.;if (present(add_tfw)) add_tfw_=add_tfw

 ngr2=0
 nvxcdgr=0
 ndvxc=0
 nd2vxc=0

!Dimension for the gradient of the density (only allocated for GGA or mGGA)
 if ((ixc>=11.and.ixc<=17).or.(ixc>=23.and.ixc<=24).or.ixc==26.or.ixc==27.or. &
& (ixc>=31.and.ixc<=34).or.(ixc==41.or.ixc==42).or.(add_tfw_)) ngr2=2*min(nspden,2)-1
 if (ixc<0.and.(libxc_functionals_isgga().or.libxc_functionals_ismgga())) &
& ngr2=2*min(nspden,2)-1

!A-Only Exc and Vxc
!=======================================================================================
 if (order**2 <= 1) then
   if (((ixc>=11 .and. ixc<=15) .and. ixc/=13) .or. (ixc>=23 .and. ixc<=24) .or. &
&   (ixc==41 .or. ixc==42)) nvxcdgr=3
   if (ixc==16.or.ixc==17.or.ixc==26.or.ixc==27) nvxcdgr=2
   if (ixc<0) nvxcdgr=3
   if (ixc>=31 .and. ixc<=34) nvxcdgr=3 !Native fake metaGGA functionals (for testing purpose only)
   if (add_tfw_) nvxcdgr=3

!  B- Exc+Vxc and other derivatives
!  =======================================================================================
 else

!  Definition of ndvxc and nvxcdgr, 2nd dimension of the arrays of 2nd-order derivatives
!  -------------------------------------------------------------------------------------
   if (ixc==1 .or. ixc==21 .or. ixc==22 .or. (ixc>=7 .and. ixc<=10) .or. ixc==13) then
!    Routine xcspol: new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option routine xcspol
!    Routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=min(nspden,2)+1
   else if (ixc>=2 .and. ixc<=6) then
!    Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)     !routine xcpzca
!    Teter fit (4/91) to Ceperley-Alder values (no spin-pol)   !routine xctetr
!    Wigner xc (no spin-pol)                                   !routine xcwign
!    Hedin-Lundqvist xc (no spin-pol)                          !routine xchelu
!    X-alpha (no spin-pol)                                     !routine xcxalp
     ndvxc=1
   else if (ixc==12 .or. ixc==24) then
!    Routine xcpbe, with optpbe=-2 and different orders (order)
     ndvxc=8
     nvxcdgr=3
   else if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. ixc==23 .or. ixc==41 .or. ixc==42) then
!    Routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=15
     nvxcdgr=3
   else if(ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27 ) then
     nvxcdgr=2
   else if (ixc==50) then
     ndvxc=1 !  IIT xc (no spin-pol)
   else if (ixc<0) then
     if(libxc_functionals_isgga() .or. libxc_functionals_ismgga()) then
       ndvxc=15
     else
       ndvxc=3
     end if
     nvxcdgr=3
   end if
   if (add_tfw_) nvxcdgr=3

!  Definition of nd2vxc, 2nd dimension of the array of 3rd-order derivatives
!  -------------------------------------------------------------------------------------
   if (order==3) then
     if (ixc==3) nd2vxc=1 ! Non spin polarized LDA case
     if ((ixc>=7 .and. ixc<=10) .or. (ixc==13)) nd2vxc=3*min(nspden,2)-2
!    Following line to be corrected when the calculation of d2vxcar is implemented for these functionals
     if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. (ixc==23.and.ixc<=24) .or. (ixc==41.or.ixc==42)) nd2vxc=1
     if ((ixc<0.and.(.not.(libxc_functionals_isgga().or.libxc_functionals_ismgga())))) &
&     nd2vxc=3*min(nspden,2)-2
   end if

 end if

end subroutine size_dvxc
!!***

!!****m* ABINIT/m_xcpositron
!! NAME
!!  m_xcpositron
!!
!! FUNCTION
!! Compute electron-positron correlation potentials and energy density.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (GJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xcpositron

 use defs_basis
 use m_errors
 use m_abicore

 use m_numeric_tools,      only : invcb

 implicit none

 private
!!***

 public :: xcpositron
!!***

contains
!!***

!!****f* ABINIT/xcpositron
!! NAME
!! xcpositron
!!
!! FUNCTION
!! Compute electron-positron correlation potentials and energy density.
!! Used electron-positron correlation functional is controlled by ixcpositron argument.
!! Returns Fxc, Vxc_pos, Vxc_el from input rhor_pos and rhor_el for positron and electrons.
!!
!! INPUTS
!!  grhoe2(ngr)=square of the gradient of electronic density rhoe (needed for GGA)
!!  ixcpositron=type of electron-positron correlation functional:
!!     1:  LDA zero positron density limit parametrized by Arponen & Pajanne
!!         and provided by Boronski & Nieminen [1,2]
!!     11: LDA zero positron density limit parametrized by Arponen & Pajanne
!!         and fitted by Sterne & Kaiser [1,3]
!!     2:  LDA electron-positron correlation
!!         provided by Puska, Seitsonen, and Nieminen [1,4]
!!     3:  GGA zero positron density limit parametrized by Arponen & Pajanne
!!         and provided by Boronski & Nieminen [1,2,5]
!!     31: GGA zero positron density limit parametrized by Arponen & Pajanne
!!         and fitted by Sterne & Kaiser [1,3,5]
!!     See references below
!!  ngr=size of grho2 array (0 if LDA, npt if GGA)
!!  npt=number of real space points on which density is provided
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rhoer(npt)=electron density (bohr^-3)
!!  rhopr(npt)=positron density (bohr^-3)
!!
!! OUTPUT
!!  fnxc(npt)=correlation energy per unit volume fxc
!!  vxce(npt)=correlation potential for electron dfxc/drhoe (hartree)
!!  vxcp(npt)=correlation potential for positron dfxc/drhop (hartree)
!!  vxcegr(ngr)= 1/|gradRhoe| dfxc/d|gradRhoe| (empty if LDA, i.e. ngr=0)
!!  Optional outputs:
!!    dvxce(npt)=partial second derivatives of the xc energy wr to the electronic density
!!               dvxce(:)=dVxce/dRhoe
!!    dvxcp(npt)=partial second derivatives of the xc energy wr to the positronic density
!!               dvxcp(:)=dVxcp/drhop
!!
!! NOTES
!!   References for electron-positron correlation functionals:
!!         [1] J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979) [[cite:Arponen1979a]].
!!         [2] E. Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986) [[cite:Boronski1986]].
!!         [3] P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991) [[cite:Sterne1991]].
!!         [4] M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994) [[cite:Puska1994]].
!!         [5] B. Barbiellini, M.J. Puska, T. Torsti and R.M.Nieminen, Phys. Rev. B 51, 7341 (1995) [[cite:Barbiellini1995]]
!!
!! PARENTS
!!      m_electronpositron,m_pawxc,m_rhotoxc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp,&
&                     dvxce,dvxcp) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixcpositron,ngr,npt
 logical,intent(in) :: posdensity0_limit
!arrays
 real(dp),intent(in) :: grhoe2(ngr),rhoer(npt),rhopr(npt)
 real(dp),intent(out) :: fnxc(npt),vxce(npt),vxcegr(ngr),vxcp(npt)
 real(dp),intent(out),optional :: dvxce(npt),dvxcp(npt)

!Local variables-------------------------------
!scalars
 integer,parameter :: idebug=0
 integer :: ipt
 logical :: gga,need_dvxce,need_dvxcp
 real(dp),parameter :: alpha_gga=0.22_dp
 real(dp),parameter :: ap_a1=-1.56_dp,ap_b1=0.051_dp,ap_c1=-0.081_dp,ap_d1=1.14_dp
 real(dp),parameter :: ap_a2=-0.92305_dp,ap_b2=-0.05459_dp
 real(dp),parameter :: ap_a3=-0.6298_dp,ap_b3=-13.15111_dp,ap_c3=2.8655_dp
 real(dp),parameter :: ap_a4=-179856.2768_dp,ap_b4=186.4207_dp,ap_c4=-0.524_dp
 real(dp),parameter :: ap_psn_limit=0.7_dp
 real(dp),parameter :: ap_psn_1=0.9_dp*ap_psn_limit,ap_psn_2=1.1_dp*ap_psn_limit
 real(dp),parameter :: fpi3=third*four_pi
 real(dp),parameter :: psn_aa=69.7029_dp,psn_ba=-107.4927_dp,psn_bb=141.8458_dp
 real(dp),parameter :: psn_ca=23.7182_dp,psn_cb=-33.6472_dp ,psn_cc=5.21152_dp
 real(dp),parameter :: sk_a=-1.56_dp,sk_b=0.1324_dp,sk_c=-4.092_dp,sk_d=51.96_dp,sk_e=0.7207_dp
 real(dp),parameter :: rsfac=0.6203504908994000_dp
 real(dp) :: arse,brse,crse,darse,dbrse,dcrse,d2arse,d2brse,d2crse
 real(dp) :: d2eps,deps,dexc,dexcdg,dexc_p,d2expgga,dexpgga,d2invrs,dinvrs,d2kf,dkf,d2nqtf2,dnqtf2
 real(dp) :: drse,drsp,d2exc,d2exc_p,d2rse,d2rsp,d2sqr,dsqr
 real(dp) :: eexp,eps,exc,exc_p,expgga,invf,dinvf,d2invf,invrhoe,invrhop,invrs,invsqr
 real(dp) :: kf,logrs,nqtf2,opr2,ratio_ap,ratio_psn,rhoe,rhop,rse,rsp,sqr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: rsepts(:),rsppts(:)

! *************************************************************************

 gga=(ngr==npt)
 need_dvxce=present(dvxce)
 need_dvxcp=present(dvxcp)

 if (gga.and.ixcpositron==2) then
   msg = 'xcpositron: GGA not yet implemented for ixcpositron=2 !'
   MSG_ERROR(msg)
 end if
 if (posdensity0_limit.and.ixcpositron==2) then
   msg = 'xcpositron: ixcpositron=2 cannot be treated in the zero positron density limit !'
   MSG_ERROR(msg)
 end if
 if (abs(ixcpositron)/=1.and.ixcpositron/=11.and.ixcpositron/=2.and.ixcpositron/=3.and.ixcpositron/=31) then
   msg = 'xcpositron: unknown electron-positron correlation !'
   MSG_ERROR(msg)
 end if

!Compute density radii for rhor_el, rhor_pos
 ABI_ALLOCATE(rsepts,(npt))
 call invcb(rhoer(:),rsepts,npt)
 rsepts(:)=rsfac*rsepts(:)
 if (ixcpositron==2) then
   ABI_ALLOCATE(rsppts,(npt))
   call invcb(rhopr(:),rsppts,npt)
   rsppts(:)=rsfac*rsppts(:)
 end if

!Loop over grid points
!----------------------------------------------------
 do ipt=1,npt

   rhoe=rhoer(ipt)
   rhop=rhopr(ipt)
   exc=zero;dexc=zero;d2exc=zero;dexcdg=zero

   rse=rsepts(ipt)
   invrhoe=one/rhoe
   drse=-third*rse*invrhoe
   if (need_dvxce) d2rse= four/nine*rse*invrhoe**2

!  Arponen & Pajane parametrization for electron
   if (ixcpositron/=11.and.ixcpositron/=31) then
     if (rse<0.302_dp) then
       invrs=one/rse;invsqr=sqrt(invrs);logrs=log(rse)
       exc =ap_a1*invsqr+(ap_b1*logrs+ap_c1)*logrs+ap_d1
       dexc=drse*invrs*(-half*ap_a1*invsqr+two*ap_b1*logrs+ap_c1)
       if (need_dvxce) d2exc=(d2rse/drse-drse*invrs)*dexc+drse**2*invrs**2*(quarter*ap_a1*invsqr+two*ap_b1)
     else if (rse>=0.302_dp.and.rse<=0.56_dp) then
       invrs=one/rse
       exc =ap_a2+ap_b2*invrs**2
       dexc=-drse*ap_b2*two*invrs**3
       if (need_dvxce) d2exc=d2rse/drse*dexc+six*drse**2*ap_b2*invrs**4
     else if (rse>0.56_dp.and.rse<=8.0_dp) then
       invrs=one/(rse+2.5_dp)
       dinvrs=-drse*invrs**2
! jmb : d2rse initialized only if need_dvxce = .True.
       if (need_dvxce) d2invrs=-d2rse*invrs**2-two*invrs*drse**2
       exc =ap_a3+ap_b3*invrs**2+ap_c3*invrs
       dexc=two*ap_b3*invrs*dinvrs+ap_c3*dinvrs
       if (need_dvxce) d2exc=two*ap_b3*dinvrs**2+(two*ap_b3*invrs+ap_c3)*d2invrs
     else
       exc  =ap_a4*rhoe**2+ap_b4*rhoe+ap_c4
       dexc =two*ap_a4*rhoe+ap_b4
       if (need_dvxce) d2exc=two*ap_a4
     end if

!    Sterne & Kaiser parametrization for electron
   else
     eexp=exp(-(rse+sk_c)**2/sk_d)
     opr2=(one+rse**2)
     arse=atan(rse)
     exc = sk_a/sqrt(arse)+sk_b*eexp+sk_e
     dexc= -(two*sk_b*eexp*(sk_c+rse)/sk_d + sk_a/(two*opr2*sqrt(arse)**3))*drse
     if (need_dvxce) d2exc=-(two*sk_b*eexp*(sk_c+rse)/sk_d + sk_a/(two*opr2*arse**1.5_dp))*d2rse &
&     +(two*sk_b*eexp*(two*sk_c**2-sk_d+four*sk_c*rse+two*rse**2)/sk_d**2 &
&     +sk_a*(three+four*rse*arse)/(four*opr2**2*sqrt(arse)**5))*drse**2
   end if

!  Puska, Seitsonen and Nieminen parametrization for positron
   if (ixcpositron==2.and.rse>=ap_psn_1) then
     rsp=rsppts(ipt)
     invrhop=one/rhop
     drsp=-third*rsp*invrhop
     if (need_dvxcp) d2rsp= four/nine*rsp*invrhop**2
     exc_p=zero;dexc_p=zero;d2exc_p=zero
     if (rsp<0.302_dp) then
       invrs=one/rsp;invsqr=sqrt(invrs);logrs=log(rsp)
       exc_p =ap_a1*invsqr+(ap_b1*logrs+ap_c1)*logrs+ap_d1
       dexc_p=drsp*invrs*(-half*ap_a1*invsqr+two*ap_b1*logrs+ap_c1)
       if (need_dvxcp) d2exc_p=(d2rsp/drsp-drsp*invrs)*dexc_p+drsp**2*invrs**2*(quarter*ap_a1*invsqr+two*ap_b1)
     else if (rsp>=0.302_dp.and.rsp<=0.56_dp) then
       invrs=one/rsp
       exc_p =ap_a2+ap_b2*invrs**2
       dexc_p=-drsp*ap_b2*two*invrs**3
       if (need_dvxcp) d2exc_p=d2rsp/drsp*dexc_p+six*drsp**2*ap_b2*invrs**4
     else if (rsp>0.56_dp.and.rsp<=8.0_dp) then
       invrs=one/(rsp+2.5_dp)
       dinvrs=-drsp*invrs**2
! jmb : d2rsp initialized only if need_dvxcp = .True.*
       if (need_dvxcp) d2invrs=-d2rsp*invrs**2-two*invrs*drsp**2
       exc_p =ap_a3+ap_b3*invrs**2+ap_c3*invrs
       dexc_p=two*ap_b3*invrs*dinvrs+ap_c3*dinvrs
       if (need_dvxcp) d2exc_p=two*ap_b3*dinvrs**2+(two*ap_b3*invrs+ap_c3)*d2invrs
     else
       exc_p  =ap_a4*rhop**2+ap_b4*rhop+ap_c4
       dexc_p =two*ap_a4*rhop+ap_b4
       if (need_dvxcp) d2exc_p=two*ap_a4
     end if
   end if

!  GGA correction
   if (gga) then
     kf=(three*pi*pi*rhoe)**third
     nqtf2=(rhoe*sqrt(four*kf/pi))**2
     eps=grhoe2(ipt)/nqtf2
     if (eps<zero) then
       MSG_ERROR('xcpositron: problem, negative GGA espilon !')
     end if
     expgga=exp(-alpha_gga*eps*third)

     dkf=pi*pi/(sqrt(three*pi*pi*rhoe)**third)
     d2kf=-two*pi*pi*pi*pi*(three*pi*pi*rhoe)**(-5.0_dp/3.0_dp)
     sqr=sqrt(four*kf/pi)
     dsqr=(four*dkf/pi)/(two*sqr)
     d2sqr=two/(pi*sqr*dkf)*(d2kf*sqr-dsqr*dkf)
     nqtf2=(rhoe*sqr)**two
     dnqtf2=two*(sqr+rhoe*dsqr)*rhoe*sqr
     d2nqtf2=two*(rhoe*sqr*(two*dsqr+rhoe*d2sqr) &
&     +sqr*(sqr+rhoe*dsqr) &
&     +rhoe*(sqr+rhoe*dsqr) )
     deps=-grhoe2(ipt)*dnqtf2/(nqtf2**two)
     d2eps=-grhoe2(ipt)/(nqtf2*nqtf2*dnqtf2)*(d2nqtf2*nqtf2*nqtf2-two*nqtf2*dnqtf2*dnqtf2)
     dexpgga=-alpha_gga*third*deps*expgga
     d2expgga=-alpha_gga*third*(d2eps*expgga+deps*dexpgga)

     exc   = exc  *expgga
     dexc=(dexc*expgga+exc*dexpgga)
     if (need_dvxce) d2exc=d2exc*expgga+two*dexc*dexpgga+exc*d2expgga
     if (abs(grhoe2(ipt))<1.e24_dp) dexcdg=-exc*alpha_gga*two_thirds/nqtf2
   end if

!  Computation of XC energy, potentials and kernels
!  Zero positron density limit
   if (ixcpositron/=2.or.rse<ap_psn_1) then
     fnxc(ipt)=rhop*exc
     vxce(ipt)=rhop*dexc
     vxcp(ipt)=exc
     if (need_dvxce) dvxce(ipt)=rhop*d2exc
     if (need_dvxcp) dvxcp(ipt)=zero
     if (gga)       vxcegr(ipt)=rhop*dexcdg
   else
!    Puska, Seitsonen and Nieminen functional
     arse=psn_aa+psn_ba*rse+psn_ca*rse**2
     brse=psn_ba+psn_bb*rse+psn_cb*rse**2
     crse=psn_ca+psn_cb*rse+psn_cc*rse**2
     darse=(psn_ba+two*psn_ca*rse)*drse
     dbrse=(psn_bb+two*psn_cb*rse)*drse
     dcrse=(psn_cb+two*psn_cc*rse)*drse
     invf=arse+brse*rsp+crse*rsp**2+invrhop/exc+invrhoe/exc_p
     fnxc(ipt)=one/invf
     dinvf=darse+dbrse*rsp+dcrse*rsp**2-invrhop*dexc/exc**2-invrhoe**2/exc_p
     vxce(ipt)=-dinvf/invf**2
     if (need_dvxce) then
       d2arse=darse*d2rse/drse+two*psn_ca*drse**2
       d2brse=dbrse*d2rse/drse+two*psn_cb*drse**2
       d2crse=dcrse*d2rse/drse+two*psn_cc*drse**2
       d2invf=d2arse+d2brse*rsp+d2crse*rsp**2 &
&       +invrhop*(two*dexc**2/exc-d2exc)/exc**2+two*invrhoe**3/exc_p
       dvxce(ipt)=(two*dinvf**2/invf-d2invf)/invf**2
     end if
     dinvf=(brse+two*crse*rsp)*drsp-invrhop**2/exc-invrhoe*dexc_p/exc_p**2
     vxcp(ipt)=-dinvf/invf**2
     if (need_dvxcp) then
       d2invf=two*crse*drsp+(brse+two*crse*rsp)*d2rsp &
&       +two*invrhop**3/exc+invrhoe*(two*dexc_p**2/exc_p-d2exc_p)/exc_p**2
       dvxcp(ipt)=(two*dinvf**2/invf-d2invf)/invf**2
     end if
!    For small rse, use pure Arponen/Pajanne functional
!    Around the limit (rse=0.7, see PSN paper), switch smoothly from PSN to AP
     if (rse>=ap_psn_1.and.rse<=ap_psn_2) then
       ratio_psn=(rse-ap_psn_1)/(ap_psn_2-ap_psn_1);ratio_ap=one-ratio_psn
       fnxc(ipt)=ratio_psn*fnxc(ipt)+ratio_ap*rhop*exc
       vxce(ipt)=ratio_psn*vxce(ipt)+ratio_ap*rhop*dexc
       vxcp(ipt)=ratio_psn*vxcp(ipt)+ratio_ap*exc
       if (need_dvxce) dvxce(ipt)=ratio_psn*dvxce(ipt)+ratio_ap*rhop*d2exc
       if (need_dvxcp) dvxcp(ipt)=ratio_psn*dvxcp(ipt)
     end if
   end if

!  Debug statements: use polynomial functionals
   if (idebug>0) then
     if (idebug==4) then ! order 4
       fnxc(ipt)=tol3*((rhop**4+rhoe**4)/12._dp+(rhop**3*rhoe+rhop*rhoe**3)/3._dp+rhop**2*rhoe**2)
       vxce(ipt)=tol3*((rhop**3*rhoe+rhop*rhoe**3)/3._dp+rhop**2*rhoe+rhop*rhoe**2)
       vxcp(ipt)=tol3*((rhop**3*rhoe+rhop*rhoe**3)/3._dp+rhop**2*rhoe+rhop*rhoe**2)
       if (need_dvxce) dvxce(ipt)=tol3*(rhop**3/3._dp+rhop*rhoe**2+rhop**2+two*rhop*rhoe)
       if (need_dvxcp) dvxcp(ipt)=tol3*(rhoe**3/3._dp+rhoe*rhop**2+rhoe**2+two*rhop*rhoe)
     end if
     if (idebug==3) then ! order 3
       fnxc(ipt)=tol3*((rhop**3+rhoe**3)*third+rhop**2*rhoe+rhop*rhoe**2)
       vxce(ipt)=tol3*(rhop+rhoe)**2
       vxcp(ipt)=tol3*(rhop+rhoe)**2
       if (need_dvxce) dvxce(ipt)=tol3*two*rhoe
       if (need_dvxcp) dvxcp(ipt)=tol3*two*rhop
     end if
     if (idebug==2) then ! order 2
       fnxc(ipt)=tol3*(rhop+rhoe)**2
       vxce(ipt)=tol3*two*(rhop+rhoe)
       vxcp(ipt)=tol3*two*(rhop+rhoe)
       if (need_dvxce) dvxce(ipt)=tol3*two
       if (need_dvxcp) dvxcp(ipt)=tol3*two
     end if
     if (idebug==1) then ! order 1
       fnxc(ipt)=tol3*(rhop+rhoe)
       vxce(ipt)=tol3
       vxcp(ipt)=tol3
       if (need_dvxce) dvxce(ipt)=zero
       if (need_dvxcp) dvxcp(ipt)=zero
     end if
   end if

 end do ! ipt

 ABI_DEALLOCATE(rsepts)
 if (ixcpositron==2) then
   ABI_DEALLOCATE(rsppts)
 end if

!Convert everything in Hartree units
 fnxc(:)=half*fnxc(:)
 vxce(:)=half*vxce(:)
 vxcp(:)=half*vxcp(:)
 if (need_dvxce) dvxce(:)=half*dvxce(:)
 if (need_dvxcp) dvxcp(:)=half*dvxcp(:)
 if (gga)       vxcegr(:)=half*vxcegr(:)

end subroutine xcpositron
!!***

end module m_xcpositron
!!***

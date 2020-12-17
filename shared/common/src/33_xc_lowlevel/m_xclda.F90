!!****m* ABINIT/m_xclda
!! NAME
!!  m_xclda
!!
!! FUNCTION
!!  LDA or LDA-like XC functionals.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LG, MF, JFD, LK)
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

module m_xclda

 use defs_basis
 use m_errors
 use m_abicore

 use m_numeric_tools,      only : invcb

 implicit none

 private
!!***

 public :: xcpzca     ! Perdew-Zunger parameterization of Ceperly-Alder electron gas energy data.
 public :: xcspol     ! Spin-polarized exchange and correlation, parameterized by Mike Teter
 public :: xctetr     ! Teter exchange and correlation --Mike Teter s fit
 public :: xcwign     ! Wigner exchange and correlation.
 public :: xchelu     ! Hedin-Lundqvist exchange and correlation
 public :: xcxalp     ! X$\alpha$ method.
 public :: xclb       ! GGA like part (vx_lb) of the Leeuwen-Baerends XC potential.
 public :: xctfw      ! Thomas-Fermi-Weizsacker functional
!!***

contains
!!***

!!****f* ABINIT/xcpzca
!! NAME
!! xcpzca
!!
!! FUNCTION
!! Returns exc, vxc, and d(vxc)/d($\rho$) from input rho.
!!
!! NOTE
!! Perdew-Zunger parameterization of Ceperly-Alder electron gas energy data.
!! J. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981) [[cite:Perdew1981]]
!! D.M. Ceperley and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980) [[cite:Ceperley1980]]
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
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
&                dvxc)                            !Optional arguments

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

!!****f* ABINIT/xcspol
!! NAME
!! xcspol
!!
!! FUNCTION
!! Spin-polarized exchange and correlation, parameterized by Mike Teter of Corning Incorporated.
!!
!! INPUTS
!!  nspden=number of spin-density components
!!  npts= number of points to be computed
!!  order=its absolute value gives the maximal derivative of Exc to be computed.
!!  rspts(npts)=Seitz electron radius (bohr)
!!  zeta(npts)=$(\rho\uparrow-\rho\downarrow)/(\rho\uparrow+\rho\downarrow)$=degree of polarization
!!  (ignored if nspden=1, in which case zeta should be 0)
!!
!! OUTPUT
!!  if(abs(order)>1) dvxc(npts,1+nspden)=              (Hartree*bohr^3)
!!   if(nspden=1 .and. order==2): dvxc(:,1)=dvxc/d$\rho$ , dvxc(:,2) empty
!!   if(nspden=1 .and. order==-2): also compute dvxc(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$
!!   if(nspden=2): dvxc(:,1)=dvxc($\uparrow$)/d$\rho(\uparrow)$,
!!       dvxc(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$, dvxc(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!
!!  exc(npts)=exchange-correlation energy density (hartree)
!!  vxc(npts,nspden)=xc potent. (d($\rho$*exc)/d($\rho\uparrow$)) and d/d($\rho\downarrow$) (ha)
!!  (only overall potential d($\rho$*exc)/d($\rho$) returned in vxc(1) for nspden=1)
!!  ndvxc= size of dvxc(npts,ndvxc)
!!
!! Normalization: Exc=$\int(exc(r)*\rho(r) d^3 r)$ for $\rho$(r)=electron density.
!!
!! TODO
!! To be added later
!!  d2vxc=derivative $d^2 (Vxc)/d(rho)^2$ (hartree*bohr^6)
!!
!! NOTES
!! This form is based on Mike Teter s rational polynomial
!! exc=-(a0+a1*rs+a2*rs**2+a3*rs**3)/(b1*rs+b2*rs**2+b3*rs**3+b4*rs**4)
!! where the parameters are fit to reproduce
!! (in this case) the Perdew-Wang parameterization of the correlation
!! energy given in Phys. Rev. B 45, 13244-13249 (1992) [[cite:Perdew1992]].
!!
!! Each parameter is interpolated between zeta=0 and 1 by
!! a_i(zeta)=a_i(0)+(a_i(1)-a_i(0))*f_x(zeta) and
!! f_x(zeta)=[(1+zeta)$^{4/3}$+(1-zeta)$^{4/3}$-2]/(2*(2$^{1/3}$-1)).
!!
!! Beware : in this expression, zeta is actually replaced by zeta*alpha_zeta,
!! where alpha_zeta is very close to 1, but slightly lower.
!! This is to remove the singularity in the derivatives when abs(zeta) is 1
!! Below,  a_i(1)-a_i(0) is called "da" for delta a, same for b s.
!!
!! rs = $(3/(4\pi))^{1/3} * \rho(r)^{-1/3}$
!! zeta = $(\rho\uparrow-\rho\downarrow)/(\rho\uparrow+\rho\downarrow)$
!! b1 must be 1 and a0 must be $(3/4)(3/(2\pi))^{2/3}$.
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,& !Mandatory arguments
&                 dvxc)                            !Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndvxc,npts,nspden,order
!arrays
 real(dp),intent(in) :: rspts(npts),zeta(npts)
 real(dp),intent(out) :: exc(npts),vxc(npts,nspden)
 real(dp),intent(out),optional :: dvxc(npts,ndvxc)

!Local variables-------------------------------
!The generation of density from rs needs rsfac and rsfac^(-3) :
!rsfac=(3/(4 Pi))^(1/3) ; rsfacm3=4pi/3
!Mike Teter s parameters of 8 April 1993.
!New parameters which accomodate spin polarization (fit to P-W)
!Paramagnetic limit:a0p,...b4p
!(a0=(3/4)(3/(2 Pi))^(2/3)) (note that b1=1 is fixed)
!Differences, ferromagnetic - paramagnetic (delta params):da1,da2,da3,db1,db2,db3,db4
!scalars
 integer :: ipts
 real(dp),parameter :: a0p=.4581652932831429_dp,a1p=2.217058676663745_dp
 real(dp),parameter :: a2p=0.7405551735357053_dp,a3p=0.01968227878617998_dp
 real(dp),parameter :: alpha_zeta=one-1.0d-6,b1p=one,b2p=4.504130959426697_dp
 real(dp),parameter :: b3p=1.110667363742916_dp,b4p=0.02359291751427506_dp
 real(dp),parameter :: da0=.119086804055547_dp,da1=0.6157402568883345_dp
 real(dp),parameter :: da2=0.1574201515892867_dp,da3=0.003532336663397157_dp
 real(dp),parameter :: db1=zero,db2=0.2673612973836267_dp
 real(dp),parameter :: db3=0.2052004607777787_dp,db4=0.004200005045691381_dp
 real(dp),parameter :: ft=4._dp/3._dp,rsfac=0.6203504908994000_dp
 real(dp),parameter :: rsfacm3=rsfac**(-3)
 real(dp) :: a0,a1,a2,a3,b1,b2,b3,b4,d1,d1m1,d2d1drs2,d2d1drsdf,d2excdf2
 real(dp) :: d2excdrs2,d2excdrsdf,d2excdz2,d2fxcdz2,d2n1drs2,d2n1drsdf,dd1df
 real(dp) :: dd1drs,dexcdf,dexcdrs,dexcdz,dfxcdz,dn1df,dn1drs,dvxcdrs
 real(dp) :: dvxcpdrho,dvxcpdz,excipt,fact,fxc,n1
 real(dp) :: rhom1,rs,vxcp,zet,zetm,zetm_third
 real(dp) :: zetp,zetp_third
 character(len=500) :: message
!no_abirules
!Set a minimum rho below which terms are 0
 real(dp),parameter :: rhotol=1.d-28
!real(dp) :: delta,rho,rho_dn,rho_dnm,rho_dnp,rho_up,rho_upm,rho_upp,zeta_mean

! *************************************************************************

!Checks the compatibility between the presence of dvxc and ndvxc
 if(ndvxc /=0 .neqv. present(dvxc))then
   message = 'If ndvxc/=0 there must be the optional argument dvxc'
   MSG_BUG(message)
 end if

!Checks the compatibility between the inputs and the presence of the optional arguments
 if(abs(order) <= 1 .and. ndvxc /= 0)then
   write(message, '(4a,i0)' )ch10,&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector dvxc, that is needed only with |order|>1 , while we have',order
   MSG_BUG(message)
 end if

 if(nspden == 1 .and. ndvxc /=0 .and. ndvxc /= 2)then
   write(message,'(a,i0)')' Once nspden=1 we must have ndvxc=2, while we have',ndvxc
   MSG_BUG(message)
 end if

 if(nspden == 2 .and. ndvxc /=0 .and. ndvxc /= 3)then
   write(message, '(a,i0)' )' Once nspden=2 we must have ndvxc=3, while we have',ndvxc
   MSG_BUG(message)
 end if


!Although fact is parameter value, some compilers are not able to evaluate
!it at compile time.
 fact=one/(two**(four*third)-two)

!DEBUG
!Finite-difference debugging, do not take away
!debug=1
!zeta_mean=0.1_dp
!delta=0.0001
!if(debug==1)then
!do ipts=1,npts,5
!rho=ipts*0.01_dp
!rho_up=rho*(one+zeta_mean)*half
!rho_dn=rho*(one-zeta_mean)*half
!rho_upp=rho_up+delta
!rho_upm=rho_up-delta
!rho_dnp=rho_dn+delta
!rho_dnm=rho_dn-delta
!First possibility : vary rho up , and then rho down
!zeta(ipts  )=(rho_up -rho_dn )/(rho_up +rho_dn )
!zeta(ipts+1)=(rho_upp-rho_dn )/(rho_upp+rho_dn )
!zeta(ipts+2)=(rho_upm-rho_dn )/(rho_upm+rho_dn )
!zeta(ipts+3)=(rho_up -rho_dnp)/(rho_up +rho_dnp)
!zeta(ipts+4)=(rho_up -rho_dnm)/(rho_up +rho_dnm)
!rspts(ipts  )=rsfac*(rho_up +rho_dn )**(-third)
!rspts(ipts+1)=rsfac*(rho_upp+rho_dn )**(-third)
!rspts(ipts+2)=rsfac*(rho_upm+rho_dn )**(-third)
!rspts(ipts+3)=rsfac*(rho_up +rho_dnp)**(-third)
!rspts(ipts+4)=rsfac*(rho_up +rho_dnm)**(-third)
!DEBUGBUG : another possibility : vary rho and zeta
!zeta(ipts+1)=zeta(ipts  )
!zeta(ipts+2)=zeta(ipts  )
!zeta(ipts+3)=zeta(ipts  )+delta
!zeta(ipts+4)=zeta(ipts  )-delta
!rspts(ipts+1)=rsfac*(rho+delta)**(-third)
!rspts(ipts+2)=rsfac*(rho-delta )**(-third)
!rspts(ipts+3)=rspts(ipts  )
!rspts(ipts+4)=rspts(ipts  )
!ENDDEBUGBUG
!end do
!end if
!nspden=2
!order=2
!ENDDEBUG

 if (nspden==1) then
!  separate cases with respect to order
   if(order==-2) then
!    No spin-polarization so skip steps related to zeta not 0
     do ipts=1,npts

       rs=rspts(ipts)
       n1=a0p+rs*(a1p+rs*(a2p+rs*a3p))
       d1=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
       d1m1=one/d1

!      Exchange-correlation energy
       excipt=-n1*d1m1
       exc(ipts)=excipt

!      Exchange-correlation potential
       dn1drs=a1p+rs*(2._dp*a2p+rs*(3._dp*a3p))
       dd1drs=b1p+rs*(2._dp*b2p+rs*(3._dp*b3p+rs*(4._dp*b4p)))

!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+excipt*dd1drs)*d1m1
       vxc(ipts,1)=excipt-third*rs*dexcdrs

!      If the exchange-correlation kernel is needed

       d2n1drs2=2._dp*a2p+rs*(6._dp*a3p)
       d2d1drs2=2._dp*b2p+rs*(6._dp*b3p+rs*(12._dp*b4p))
!      d2excdrs2 is d2(exc)/d(rs)2
       d2excdrs2=-(d2n1drs2+2._dp*dexcdrs*dd1drs+excipt*d2d1drs2)*d1m1
       dvxcdrs=third*(2.0_dp*dexcdrs-rs*d2excdrs2)
!      And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
       dvxc(ipts,1)= -rs**4*rsfacm3*third*dvxcdrs

!      dn1df=d(n1)/d(fxc) and dd1df=d(d1)/d(fxc)
       dn1df=da0+rs*(da1+rs*(da2+rs*da3))
       dd1df=rs*(db1+rs*(db2+rs*(db3+rs*db4)))
       dexcdf=-(dn1df+excipt*dd1df)*d1m1
!      d2(fxc)/d(zeta)2
       d2fxcdz2=ft*third*(alpha_zeta**2)*2._dp*fact
!      d2(exc)/d(zeta)2
       d2excdz2=d2fxcdz2*dexcdf
       rhom1=rsfacm3*rs**3
       dvxc(ipts,2)= dvxc(ipts,1) - d2excdz2*rhom1
     end do
   else if(order**2>1) then
!    No spin-polarization so skip steps related to zeta not 0
     do ipts=1,npts

       rs=rspts(ipts)
       n1=a0p+rs*(a1p+rs*(a2p+rs*a3p))
       d1=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
       d1m1=one/d1

!      Exchange-correlation energy
       excipt=-n1*d1m1
       exc(ipts)=excipt

!      Exchange-correlation potential
       dn1drs=a1p+rs*(2._dp*a2p+rs*(3._dp*a3p))
       dd1drs=b1p+rs*(2._dp*b2p+rs*(3._dp*b3p+rs*(4._dp*b4p)))

!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+excipt*dd1drs)*d1m1
       vxc(ipts,1)=excipt-third*rs*dexcdrs

!      If the exchange-correlation kernel is needed
       d2n1drs2=2._dp*a2p+rs*(6._dp*a3p)
       d2d1drs2=2._dp*b2p+rs*(6._dp*b3p+rs*(12._dp*b4p))
!      d2excdrs2 is d2(exc)/d(rs)2
       d2excdrs2=-(d2n1drs2+2._dp*dexcdrs*dd1drs+excipt*d2d1drs2)*d1m1
       dvxcdrs=third*(2.0_dp*dexcdrs-rs*d2excdrs2)
!      And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
       dvxc(ipts,1)= -rs**4*rsfacm3*third*dvxcdrs

     end do
   else
!    No spin-polarization so skip steps related to zeta not 0
     do ipts=1,npts

       rs=rspts(ipts)
       n1=a0p+rs*(a1p+rs*(a2p+rs*a3p))
       d1=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
       d1m1=one/d1

!      Exchange-correlation energy
       excipt=-n1*d1m1
       exc(ipts)=excipt

!      Exchange-correlation potential
       dn1drs=a1p+rs*(2._dp*a2p+rs*(3._dp*a3p))
       dd1drs=b1p+rs*(2._dp*b2p+rs*(3._dp*b3p+rs*(4._dp*b4p)))

!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+excipt*dd1drs)*d1m1
       vxc(ipts,1)=excipt-third*rs*dexcdrs
     end do

   end if


!  Allows for nspden==1, in the case of testing nspden=1 against nspden=2
 else if (nspden<=2) then


!  DEBUG
!  do not take away : allows to compare nspden=1 and nspden=2 coding
!  if (nspden==1)then
!  zeta(:)=zero
!  end if
!  ENDDEBUG
!  separate cases with respect to order
   if(abs(order)>1) then
!    Allow for spin polarization. This part could be optimized for speed.
     do ipts=1,npts

       rs=rspts(ipts)
       zet=zeta(ipts)
       zetp=one+zet*alpha_zeta
       zetm=one-zet*alpha_zeta
       zetp_third=zetp**third
       zetm_third=zetm**third
!      Exchange energy spin interpolation function f(zeta)
       fxc=( zetp*zetp_third + zetm*zetm_third - two ) *fact

       a0=a0p+fxc*da0
       a1=a1p+fxc*da1
       a2=a2p+fxc*da2
       a3=a3p+fxc*da3
       b1=b1p+fxc*db1
       b2=b2p+fxc*db2
       b3=b3p+fxc*db3
       b4=b4p+fxc*db4

       n1= a0+rs*(a1+rs*(a2+rs*a3))
       d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
       d1m1=one/d1

!      Exchange-correlation energy
       excipt=-n1*d1m1
       exc(ipts)=excipt

!      Exchange-correlation potential
       dn1drs=a1+rs*(2._dp*a2+rs*(3._dp*a3))
       dd1drs=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+excipt*dd1drs)*d1m1

!      Only vxcp contributes when paramagnetic
       vxcp=excipt-third*rs*dexcdrs

!      d(fxc)/d(zeta)  (which is 0 at zeta=0)
       dfxcdz=ft*alpha_zeta*(zetp_third-zetm_third)*fact

!      dn1df=d(n1)/d(fxc) and dd1df=d(d1)/d(fxc)
       dn1df=da0+rs*(da1+rs*(da2+rs*da3))
       dd1df=rs*(db1+rs*(db2+rs*(db3+rs*db4)))

!      dexcdz is d(exc)/d(zeta)
       dexcdf=-(dn1df+excipt*dd1df)*d1m1
       dexcdz=dfxcdz*dexcdf

!      Compute Vxc for both spin channels

       vxc(ipts,1)=vxcp - (zet-one)*dexcdz
       vxc(ipts,2)=vxcp - (zet+one)*dexcdz

!      DEBUG Allow to check the variation of rho and zeta
!      vxc(ipts,1)=vxcp
!      vxc(ipts,2)=dexcdz
!      ENDDEBUG
!      Compute second derivative with respect to rho
       d2n1drs2=2._dp*a2+rs*(6._dp*a3)
       d2d1drs2=2._dp*b2+rs*(6._dp*b3+rs*(12._dp*b4))
!      d2excdrs2 is d2(exc)/d(rs)2
       d2excdrs2=-(d2n1drs2+two*dexcdrs*dd1drs+excipt*d2d1drs2)*d1m1
       dvxcdrs=third*(two*dexcdrs-rs*d2excdrs2)
!      And d(vxc)/d(rho) paramagnetic =(-rs/(3*rho))*d(vxcp)/d(rs)
!      remember : 1/rho=(4pi/3)*rs**3=rsfacm3*rs**3
       rhom1=rsfacm3*rs**3
       dvxcpdrho= -rs*rhom1*third * dvxcdrs

!      Compute mixed second derivative with respect to rho and zeta
       d2n1drsdf=da1+rs*(2._dp*da2+rs*(3._dp*da3))
       d2d1drsdf=db1+rs*(2._dp*db2+rs*(3._dp*db3+rs*(4._dp*db4)))
!      d2excdrsdf is d2(exc)/d(rs)df
       d2excdrsdf=-(d2n1drsdf+dexcdrs*dd1df+dexcdf*dd1drs+excipt*d2d1drsdf)*d1m1
!      d(vxc)/d(zeta) paramagnetic
       dvxcpdz=dexcdz-third*rs*dfxcdz*d2excdrsdf

!      Compute second derivative with respect to zeta
!      the second derivative of n1 and d1 wrt f vanishes
       d2excdf2=-(two*dexcdf*dd1df)*d1m1
!      d2(fxc)/d(zeta)2
       d2fxcdz2=ft*third*(alpha_zeta**2)*(zetp_third**(-2)+zetm_third**(-2))*fact
!      d2(exc)/d(zeta)2
       d2excdz2=d2fxcdz2*dexcdf+dfxcdz**2*d2excdf2

!      Compute now the three second derivatives of the Exc energy with respect
!      to : wrt twice spin-up ; wrt spin-up and spin-dn ; wrt twice spin-down
       dvxc(ipts,1)= dvxcpdrho   &
&       +two*rhom1*( one-zet)*(dvxcpdz-dexcdz) &
&       +d2excdz2*rhom1*(one-zet)**2
       dvxc(ipts,2)= dvxcpdrho   &
&       +two*rhom1*(    -zet)*(dvxcpdz-dexcdz) &
&       +d2excdz2*rhom1*(one-zet)*(-one-zet)
!      if(nspden==2)then
       dvxc(ipts,3)= dvxcpdrho   &
&       +two*rhom1*(-one-zet)*(dvxcpdz-dexcdz) &
&       +d2excdz2*rhom1*(-one-zet)**2
!      else
!      !    For testing purposes, need the spin-averaged quantity
!      dvxc(ipts,1)= ( dvxc(ipts,1) + dvxc(ipts,2) ) * half
!      end if

!      DEBUG Allow to check the variation of rho and zeta
!      dvxc(ipts,1)=dvxcpdrho
!      dvxc(ipts,2)=d2excdz2
!      dvxc(ipts,3)=dvxcpdz
!      ENDDEBUG
     end do
   else
!    Allow for spin polarization. This part could be optimized for speed.
     do ipts=1,npts

       rs=rspts(ipts)
       zet=zeta(ipts)
       zetp=one+zet*alpha_zeta
       zetm=one-zet*alpha_zeta
       zetp_third=zetp**third
       zetm_third=zetm**third
!      Exchange energy spin interpolation function f(zeta)
       fxc=( zetp*zetp_third + zetm*zetm_third - two ) *fact

       a0=a0p+fxc*da0
       a1=a1p+fxc*da1
       a2=a2p+fxc*da2
       a3=a3p+fxc*da3
       b1=b1p+fxc*db1
       b2=b2p+fxc*db2
       b3=b3p+fxc*db3
       b4=b4p+fxc*db4

       n1= a0+rs*(a1+rs*(a2+rs*a3))
       d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
       d1m1=one/d1

!      Exchange-correlation energy
       excipt=-n1*d1m1
       exc(ipts)=excipt

!      Exchange-correlation potential
       dn1drs=a1+rs*(2._dp*a2+rs*(3._dp*a3))
       dd1drs=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
!      dexcdrs is d(exc)/d(rs)
       dexcdrs=-(dn1drs+excipt*dd1drs)*d1m1

!      Only vxcp contributes when paramagnetic
       vxcp=excipt-third*rs*dexcdrs

!      d(fxc)/d(zeta)  (which is 0 at zeta=0)
       dfxcdz=ft*alpha_zeta*(zetp_third-zetm_third)*fact

!      dn1df=d(n1)/d(fxc) and dd1df=d(d1)/d(fxc)
       dn1df=da0+rs*(da1+rs*(da2+rs*da3))
       dd1df=rs*(db1+rs*(db2+rs*(db3+rs*db4)))

!      dexcdz is d(exc)/d(zeta)
       dexcdf=-(dn1df+excipt*dd1df)*d1m1
       dexcdz=dfxcdz*dexcdf

!      Compute Vxc for both spin channels

       vxc(ipts,1)=vxcp - (zet-one)*dexcdz
       vxc(ipts,2)=vxcp - (zet+one)*dexcdz

!      DEBUG Allow to check the variation of rho and zeta
!      vxc(ipts,1)=vxcp
!      vxc(ipts,2)=dexcdz
!      ENDDEBUG
     end do
   end if
 else

!  Disallowed value for nspden
   write(message, '(3a,i0)' )&
&   ' Argument nspden must be 1 or 2; ',ch10,&
&   ' Value provided as argument was ',nspden
   MSG_BUG(message)
 end if

!DEBUG
!Finite-difference debugging, do not take away
!if(debug==1)then
!write(std_out,*)' delta =',delta
!do ipts=1,npts,5
!rho=(rspts(ipts)/rsfac)**(-3)
!write(std_out,'(a,i5,a,2es16.8)' ) ' Point number',ipts,' with rho,zeta=',rho,zeta(ipts)
!write(std_out,'(3es16.8)' )exc(ipts)*rho,vxc(ipts,1),vxc(ipts,2)
!write(std_out,'(3es16.8)' )dvxc(ipts,1),dvxc(ipts,3),dvxc(ipts,2)
!write(std_out,'(3es16.8)' )exc(ipts)*rho,&
!&      ( exc(ipts+1)*(rho+delta) - exc(ipts+2)*(rho-delta) )/2._dp/delta,&
!&      ( exc(ipts+3)*(rho+delta) - exc(ipts+4)*(rho-delta) )/2._dp/delta
!write(std_out,'(4es16.8)' )&
!&    ( vxc(ipts+1,1) - vxc(ipts+2,1) )/2._dp/delta,&
!&    ( vxc(ipts+3,2) - vxc(ipts+4,2) )/2._dp/delta,&
!&    ( vxc(ipts+3,1) - vxc(ipts+4,1) )/2._dp/delta,&
!&    ( vxc(ipts+1,2) - vxc(ipts+2,2) )/2._dp/delta
!end do
!stop
!end if
!ENDDEBUG

!DEBUG
!if(order==-2)then
!write(std_out,*)' xcspol : ipts,npts ',ipts,npts
!write(std_out,*)dvxcdrs,d2excdz2,d2fxcdz2,dexcdf
!write(std_out,*)rhom1
!write(std_out,*)dvxc(1000,1),dvxc(1000,2)
!stop
!end if
!ENDDEBUG

end subroutine xcspol
!!***


!!****f* ABINIT/xctetr
!! NAME
!! xctetr
!!
!! FUNCTION
!! Returns exc, vxc, and d(vxc)/d($\rho$) from input $\rho$.
!! Also returns $d^2(Vxc)/d(\rho)^2$ as needed for third-order DFPT
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rhor(npt)=electron number density (bohr^-3)
!!  rspts(npt)=corresponding Wigner-Seitz radii, precomputed
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d(rho*exc)/d(rho)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!  if(order>2) d2vxc(npt)=derivative d$^2$(Vxc)/d$(\rho)^2$ (hartree*bohr^6)
!!
!! NOTES
!! Teter exchange and correlation (xc)--Mike Teter s fit
!! to Ceperly-Alder electron gas energy data.  Data from
!! D.M. Ceperley and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980) [[cite:Ceperley1980]]
!! and private communication from authors.
!! This form is based on Mike Teter s rational polynomial
!! exc=-(a0+a1*rs+a2*rs**2+a3*rs**3)/(b1*rs+b2*rs**2+b3*rs**3+b4*rs**4)
!! where the parameters of the fit are fit to reproduce
!! Ceperley-Alder data and the high density limit (rs->0)
!! of the electron gas (pure exchange).
!! rs = $(3/(4\pi))^{1/3} * \rho(r)^{-1/3}$.
!! b1 must be 1 and a0 must be $(3/4)(3/(2\pi))^{2/3}$.
!! Fit is by Mike Teter, Corning Incorporated.
!! Note that d(vxc)/d($\rho$) gets a little wild at small rho.
!! d$^2$(Vxc)/d$(\rho)^2$ is probably wilder.
!!
!! Some notation:  (XG 990224, sign convention should be changed, see xcspol.f)
!!  $Exc = N1/D1$ with $N1=-(a0+a1*rs+...)$ given above and
!!              $D1= (b1*rs+b2*rs^2+...)$ also given above.
!!  $Vxc = N2/D1^2$ with $N2=d(N1)/d(rs)$.
!!  $d(Vxc)/d(rs)=(N3-D3*(2*N2/D1))/D1^2 with N3=d(N2)/d(rs)$ and
!!              $D3=d(D1)/d(rs)$.
!!  $d(Vxc)/d(\rho) = (-rs/(3*\rho))* d(Vxc)/d(rs)$.
!!  $d^2(Vxc)/d(rs)^2 = (N4-2*(2*N3*D3+N2*D4-3*N2*D3^2/D1)/D1)/D1^2$
!!   with $N4=d(N3)/d(rs), D4=d(D3)/d(rs)$.
!!  $d^2(Vxc)/d(\rho)^2= rs/(3*\rho)^2)*(4*d(Vxc)/d(rs)+rs*d^2(Vxc)/d(rs)^2)$.
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xctetr(exc,npt,order,rhor,rspts,vxc,& !Mandatory arguments
&                 d2vxc,dvxc)                    !Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rhor(npt),rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: d2vxc(npt),dvxc(npt)

!Local variables-------------------------------
!rsfac=(3/(4 Pi))^(1/3)
!Mike Teter s parameters: (keep 8 digits after decimal)
!(a0=(3/4)(3/(2 Pi))^(2/3)
!scalars
 integer :: ipt
 real(dp),parameter :: a0=.4581652932831429_dp,a1=2.40875407_dp,a2=.88642404_dp
 real(dp),parameter :: a3=.02600342_dp,b1=1.0_dp,b2=4.91962865_dp
 real(dp),parameter :: b3=1.34799453_dp,b4=.03120453_dp,c1=4._dp*a0*b1/3.0_dp
 real(dp),parameter :: c2=5.0_dp*a0*b2/3.0_dp+a1*b1
 real(dp),parameter :: c3=2.0_dp*a0*b3+4.0_dp*a1*b2/3.0_dp+2.0_dp*a2*b1/3.0_dp
 real(dp),parameter :: c4=7.0_dp*a0*b4/3.0_dp+5.0_dp*a1*b3/3.0_dp+a2*b2+a3*b1/3.0_dp
 real(dp),parameter :: c5=2.0_dp*a1*b4+4.0_dp*a2*b3/3.0_dp+2.0_dp*a3*b2/3.0_dp
 real(dp),parameter :: c6=5.0_dp*a2*b4/3.0_dp+a3*b3,c7=4.0_dp*a3*b4/3.0_dp
 real(dp),parameter :: rsfac=0.6203504908994000_dp
 real(dp) :: d1,d1m1,d2vxcr,d3,d4,dvxcdr,n1,n2,n3,n4,rhom1,rs
 character(len=500) :: message

! *************************************************************************
!
!Checks the values of order
 if(order<0 .or. order>3)then
   write(message, '(a,a,a,i6)' )&
&   'With Teter 91 Ceperley-Alder xc functional, the only',ch10,&
&   'allowed values for order are 0, 1, 2 or 3, while it is found to be',order
   MSG_BUG(message)
 end if

!Checks the compatibility between the order and the presence of the optional arguments
 if(order /=3 .and. present(d2vxc))then
   write(message, '(a,a,a,i6)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector d2vxc, that is needed only with order=3, while we have',order
   MSG_BUG(message)
 end if

 if(order <= 1 .and. present(dvxc))then
   write(message, '(a,a,a,i6)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector dvxc, that is needed with order > 1, while we have',order
   MSG_BUG(message)
 end if

!separated cases with respect to order

 if (order<=1) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     n1=-(a0+rs*(a1+rs*(a2+rs*a3)))
     d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
     d1m1=1.0_dp/d1
     n2=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
!
!    Exchange-correlation energy
     exc(ipt)=n1*d1m1
!
!    Exchange-correlation potential
     vxc(ipt)=n2*d1m1**2
   end do
 else if (order>2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     n1=-(a0+rs*(a1+rs*(a2+rs*a3)))
     d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
     d1m1=1.0_dp/d1
     n2=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
!
!    Exchange-correlation energy
     exc(ipt)=n1*d1m1
!
!    Exchange-correlation potential
     vxc(ipt)=n2*d1m1**2
!    Assemble derivative of vxc wrt rs
     n3=-(c1+rs*(2._dp*c2+rs*(3._dp*c3+rs*(4._dp*c4+rs*(5._dp*c5+&
&     rs*(6._dp*c6+rs*(7._dp*c7)))))))
     d3=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
     dvxcdr=(n3-d3*(2._dp*n2*d1m1))*d1m1**2
     rhom1=1.0_dp/rhor(ipt)
!
!    derivative of vxc wrt rho
     dvxc(ipt)=-dvxcdr*rs*third*rhom1
!

!    Assemble derivative d^2(Vxc)/d(rs)^2
     n4=-(2.0_dp*c2+rs*(6.0_dp*c3+rs*(12.0_dp*c4+rs*(20.0_dp*c5+&
&     rs*(30.0_dp*c6+rs*(42.0_dp*c7))))))
     d4=2.0_dp*b2+rs*(6.0_dp*b3+rs*(12.0_dp*b4))
     d2vxcr=(n4-2.0_dp*(2.0_dp*n3*d3+n2*d4-3.0_dp*n2*d3**2*d1m1)*d1m1)*d1m1**2

!    Derivative d^2(Vxc)/d(rho)^2
     d2vxc(ipt)=(rs*third*rhom1)*(4.0_dp*dvxcdr+rs*d2vxcr)*third*rhom1

   end do
 else if (order>1) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     n1=-(a0+rs*(a1+rs*(a2+rs*a3)))
     d1=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
     d1m1=1.0_dp/d1
     n2=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
!
!    Exchange-correlation energy
     exc(ipt)=n1*d1m1
!
!    Exchange-correlation potential
     vxc(ipt)=n2*d1m1**2
!    Assemble derivative of vxc wrt rs
     n3=-(c1+rs*(2._dp*c2+rs*(3._dp*c3+rs*(4._dp*c4+rs*(5._dp*c5+&
&     rs*(6._dp*c6+rs*(7._dp*c7)))))))
     d3=b1+rs*(2._dp*b2+rs*(3._dp*b3+rs*(4._dp*b4)))
     dvxcdr=(n3-d3*(2._dp*n2*d1m1))*d1m1**2
     rhom1=1.0_dp/rhor(ipt)
!
!    derivative of vxc wrt rho
     dvxc(ipt)=-dvxcdr*rs*third*rhom1
!
   end do
 end if
end subroutine xctetr
!!***

!!****f* ABINIT/xcwign
!! NAME
!! xcwign
!!
!! FUNCTION
!! Returns exc, vxc, and eventually d(vxc)/d($\rho$) from input $\rho$.
!! Wigner exchange and correlation (xc)--see e.g. David Pines,
!! Elementary Excitations in Solids, p. 94, NY 1964.
!! Expression is exc=-(0.44)/(rs+7.8)-efac/rs (hartree), efac below.
!! rs = $(3/(4\pi))^{1/3}* \rho (r)^{-1/3}$.
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=corresponding Wigner-Seitz radii, precomputed
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d($\rho$*exc)/d($\rho$)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xcwign(exc,npt,order,rspts,vxc,& !Mandatory arguments
&                dvxc)                           !Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!c1 and c2 are the Wigner parameters in hartree and bohr resp.
!scalars
 integer :: ipt
 real(dp),parameter :: c1=0.44_dp,c2=7.8_dp,c4_3=4.0_dp/3.0_dp
 real(dp),parameter :: c8_27=8.0_dp/27.0_dp
 real(dp) :: dfac,efac,rs,rsc2m1,rsm1,vfac,vxcnum
 character(len=500) :: message

! *************************************************************************

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,i0)' )&
&   'With Wigner xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be',order
   MSG_BUG(message)
 end if

!Checks the compatibility between the order and the presence of the optional arguments
 if(order <= 1 .and. present(dvxc))then
   write(message, '(a,a,a,i3)' )&
&   'The order chosen does not need the presence',ch10,&
&   'of the vector dvxc, that is needed only with order=2 , while we have',order
   MSG_BUG(message)
 end if

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac

!separate cases with respect to order
 if (order==2) then

!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
     rsc2m1=1.0_dp/(rs+c2)
!    compute energy density (hartree)
     exc(ipt)=-c1*rsc2m1-efac*rsm1
     vxcnum=-(c4_3*rs+c2)*c1
!    compute potential (hartree)
     vxc(ipt)=vxcnum*rsc2m1**2-vfac*rsm1
!    compute d(vxc)/d(rho) (hartree*bohr^3)
     dvxc(ipt)=-(c8_27*pi)*(c1*rs**4)*(rs+rs+c2)*rsc2m1**3-dfac*rs**2
   end do
 else

!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
     rsc2m1=1.0_dp/(rs+c2)
!    compute energy density (hartree)
     exc(ipt)=-c1*rsc2m1-efac*rsm1
     vxcnum=-(c4_3*rs+c2)*c1
!    compute potential (hartree)
     vxc(ipt)=vxcnum*rsc2m1**2-vfac*rsm1
   end do

 end if
!
end subroutine xcwign
!!***


!!****f* ABINIT/xchelu
!! NAME
!! xchelu
!!
!! FUNCTION
!! Returns exc, vxc, and eventually d(vxc)/d($\rho$) from input rho.
!!
!! NOTES
!! Hedin-Lundqvist exchange and correlation (xc)--
!! L. Hedin and B.I. Lundqvist, J. Phys. C. 4, 2064 (1971) [[cite:Hedin1971]]
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=Wigner-Seitz radii at each point
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d($\rho$*exc)/d($\rho$)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!aa and cc are H-L fitting parameters A and C (C in hartree)
!rs = (3/(4 Pi))**(1/3) * rho(r)**(-1/3).
!scalars
 integer :: ipt
 real(dp),parameter :: aa=21_dp,c1_21=one/21_dp,c4_9=4.0_dp/9.0_dp,cc=0.0225_dp
 real(dp) :: dfac,efac,rs,rsm1,vfac,xx
 character(len=500) :: message

! *************************************************************************

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,i0)' )&
&   'With Hedin-Lundqvist xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be',order
   MSG_BUG(message)
 end if

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac
!separate cases with respect to order
 if (order==2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=one/rs
!    compute energy density exc (hartree)
     xx=rs*c1_21
     exc(ipt)=-cc*((one+xx**3)*log(one+one/xx)+&
&     half*xx-xx*xx-third) - efac*rsm1
!    compute xc potential d(rho*exc)/d(rho) (hartree)
     vxc(ipt)=-cc*log(one+aa*rsm1)-vfac*rsm1
!    compute d(vxc)/d(rho) (hartree*bohr^3)
     dvxc(ipt)=-(rs**2)*((c4_9*pi)*cc*rs/(one+xx) + dfac)
   end do
 else
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=one/rs
!    compute energy density exc (hartree)
     xx=rs*c1_21
     exc(ipt)=-cc*((one+xx**3)*log(one+one/xx)+&
&     half*xx-xx*xx-third) - efac*rsm1
!    compute xc potential d(rho*exc)/d(rho) (hartree)
     vxc(ipt)=-cc*log(one+aa*rsm1)-vfac*rsm1
   end do
 end if
!
end subroutine xchelu
!!***

!!****f* ABINIT/xcxalp
!! NAME
!! xcxalp
!!
!! FUNCTION
!! Returns exc, vxc, and eventually d(vxc)/d($\rho$) from input $\rho$.
!! "X$\alpha$" method is used in this subroutine:
!! a single fixed value is chosen for "alpha", set below.
!! Expression is exc=-alpha*efac/rs (hartree), efac below.
!! rs = $(3/(4\pi))^{1/3}* \rho (r)^{-1/3}$.
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=Wigner-Seitz radii, at each point
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d($\rho$*exc)/d($\rho$)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xcxalp(exc,npt,order,rspts,vxc, dvxc)  ! dvxc is optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!Set value of alpha in "X-alpha" method
!scalars
 integer :: ipt
 real(dp),parameter :: alpha=1.0_dp
 real(dp) :: dfac,efac,rs,rsm1,vfac
 character(len=500) :: message

! *************************************************************************

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,i3)' )&
&   'With X-alpha xc functional, the only',ch10,&
&   'allowed values for order are 0, 1 or 2, while it is found to be',order
   MSG_BUG(message)
 end if

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac

!separate cases with respect to order
 if(order==2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
!    compute energy density (hartree)
     exc(ipt)=-alpha*efac*rsm1
!    compute potential (hartree)
     vxc(ipt)=-alpha*vfac*rsm1
!    compute d(vxc)/d(rho) (hartree*bohr^3)
     dvxc(ipt)=-alpha*dfac*rs**2
   end do
 else
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=1.0_dp/rs
!    compute energy density (hartree)
     exc(ipt)=-alpha*efac*rsm1
!    compute potential (hartree)
     vxc(ipt)=-alpha*vfac*rsm1
   end do
 end if
!
end subroutine xcxalp
!!***

!!****f* ABINIT/xclb
!! NAME
!! xclb
!!
!! FUNCTION
!! Computes the GGA like part (vx_lb) of the Leeuwen-Baerends
!! exchange-correlation potential (vxc_lb) and adds it to the
!! lda exchange-correlation potential (vxc_lda) which
!! must be provided as input,
!!            vxci  <--  vxc_lb =: vxc_lda + vx_lb
!!
!! R van Leeuwen and EJ Baerends, Phys Rev A 49, 2421 (1994) [[cite:VanLeeuwen1994]]
!!
!! With respect to spin, the van Leeuwen-Baerends
!! potential is "exchange-like" : separate contributions from
!! spin up and spin down.
!!
!! INPUTS
!!  npts= number of points to be computed
!!  nspden=1 for unpolarized, 2 for spin-polarized
!!  grho2_updn(npts,2*nspden-1)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2)
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  vxci(npts,nspden)=input xc potential to which Leeuwen-Baerends correction
!!   is added at output.
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xclb(grho2_updn,npts,nspden,rho_updn,vxci)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts,nspden
!arrays
 real(dp),intent(in) :: grho2_updn(npts,2*nspden-1),rho_updn(npts,nspden)
 real(dp),intent(inout) :: vxci(npts,nspden)

!Local variables-------------------------------
!scalars
 integer :: ipts,ispden
 real(dp),parameter :: beta=0.05_dp
 real(dp) :: density,density_gradient,density_t13,s_g_sq,scaled_gradient
 real(dp) :: scaling_factor,vx_lb

! *************************************************************************

!DEBUG
!write(std_out,*) ' %xclb: enter'
!ENDDEBUG

!scale the spin densities for evaluating spin up or down exchange
 scaling_factor=one
 if(nspden == 2) scaling_factor=two

 do ispden=1,nspden

   do ipts=1,npts

     density= scaling_factor * rho_updn(ipts,ispden)
     density_gradient= scaling_factor * sqrt(grho2_updn(ipts,ispden))

     density_t13= density**third
     scaled_gradient= density_gradient/max(density*density_t13,1.e-12_dp)

     s_g_sq= scaled_gradient*scaled_gradient

     vx_lb= -beta*density_t13 * s_g_sq/ &
&     (one+3.d0*beta* scaled_gradient*log(scaled_gradient+sqrt(one+s_g_sq*s_g_sq)))

     vxci(ipts,ispden)=vxci(ipts,ispden)+vx_lb
   end do

 end do

end subroutine xclb
!!***

!!****f* ABINIT/xctfw
!! NAME
!! xctfw
!!
!! FUNCTION
!! Add gradient part of the Thomas-Fermi-Weizsacker functional
!! Perrot F., Phys. Rev. A20, 586-594 (1979) [[cite:Perrot1979]]
!!
!! INPUTS
!!  ndvxcdgr= size of dvxcdgr(npts,ndvxcdgr)
!!  npts= number of points to be computed
!!  nspden=number if spin density component (necessarily 1 here)
!!  grho2_updn(npts,2*nspden-1)=square of the gradient of the spin-up,
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
!!      m_rhotoxc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xctfw(temp,exci,fxci,usefxc,rho_updn,vxci,npts,nspden,dvxcdgr,ndvxcdgr,grho2_updn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndvxcdgr,npts,nspden,usefxc
 real(dp),intent(in) :: temp
!arrays
 real(dp),intent(in) :: grho2_updn(npts,2*nspden-1),rho_updn(npts,nspden)
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

end subroutine xctfw
!!***

end module m_xclda
!!***

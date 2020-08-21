!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xchcth
!! NAME
!!  m_xchcth
!!
!! FUNCTION
!! Treat XC GGA functional of HCTH type.
!! See Hamprecht, Cohen, Tozer and Handy, J. Chem. Phys. 109, 6264 (1998) [[cite:Hamprecht1998]] for HCTH-93
!!  Boese, Doltsinis, Handy and Sprik, J. Chem. Phys. 112, 1670 (2000) [[cite:Boese2000]] for HCTH-120 and HCTH-147.
!!  Boese and Handy , J. Chem. Phys. 114, 5497 (2001) [[cite:Boese2001]] for HCTH-407.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2020 ABINIT group (XG,LG)
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

module m_xchcth

 use defs_basis
 use m_errors
 use m_abicore

 use m_numeric_tools,      only : invcb

 implicit none

 private
!!***

 public :: xchcth
!!***

contains
!!***

!!****f* ABINIT/xchcth
!! NAME
!! xchcth
!!
!! FUNCTION
!! Treat XC GGA functional of HCTH type.
!! See Hamprecht, Cohen, Tozer and Handy, J. Chem. Phys. 109, 6264 (1998) [[cite:Hamprecht1998]] for HCTH-93.
!!  Boese, Doltsinis, Handy and Sprik, J. Chem. Phys. 112, 1670 (2000) [[cite:Boese2000]] for HCTH-120 and HCTH-147.
!!  Boese and Handy , J. Chem. Phys. 114, 5497 (2001) [[cite:Boese2001]] for HCTH-407.
!!
!! For a series of values of the density and the square of the
!! gradient of the density, return the associated Exc energy,
!! potential, and, in case of response-function, functions needed
!! to build the XC kernel.
!!
!! INPUTS
!!  ixc=number of the XC functional : 16 for HCTH-93, 17 for HCTH-120, 26 for HCTH-147 and 27 for HCTH-407.
!!  npts= number of points to be computed
!!  nspden=1 for unpolarized, 2 for spin-polarized
!!  grho2_updn(npts,2*nspden-1)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2)
!!  order=maximal derivative of Exc to be computed
!!   (1 => energy and potential, or 2 => also XC kernel )
!!   Warning : order=2 not yet available
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!
!! OUTPUT
!!
!!  dvxcdgr(npts,3)=partial derivative of the exchange-correlation
!!    energy (exci*$\rho$) with respect to the spin-up (dvxcdgr(:,1)),
!!    spin-down (dvxcdgr(:,2)) square of gradients of the density
!!
!!  exci(npts)=exchange-correlation energy density (hartree)
!!  vxci(npts,nspden)=partial derivative of the exchange-correlation energy (exci*$\rho$)
!!    with respect to the spin-down (vxci(:,1)) and spin-up (vxci(:,2) densities
!! Normalization: Exc=$\int (exc(r)*\rho (r) d^3 r)$ for $\rho$(r)=electron density.
!!
!! TODO
!! Response function not coded yet, but part of it are already present
!!
!! PARENTS
!!      m_drivexc
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine xchcth(dvxcdgr,exci,grho2_updn,ixc,npts,nspden,order,rho_updn,vxci)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,npts,nspden,order
!arrays
 real(dp),intent(in) :: grho2_updn(npts,2*nspden-1),rho_updn(npts,nspden)
 real(dp),intent(out) :: dvxcdgr(npts,2),exci(npts),vxci(npts,nspden)

!Local variables-------------------------------
!scalars
 integer,save :: initialized=0
 integer :: ipts,ispden
 real(dp),parameter :: alpha_zeta2=one-1e-6_dp,alpha_zeta=one-1e-6_dp
 real(dp),parameter :: fsec_inv=one/1.709921_dp,gammacab=0.006_dp
 real(dp),parameter :: gammacsig=0.2_dp,gammax=0.004_dp
 real(dp),parameter :: rsfac=0.6203504908994000_dp
 real(dp),save :: factf_zeta,factfp_zeta,sixpi2_1_3,sixpi2m1_3,sq_rsfac
 real(dp),save :: sq_rsfac_inv,threefourth_divpi,twom1_3
 real(dp) :: ccab0,ccab1,ccab2,ccab3,ccab4,ccsig0,ccsig1,ccsig2,ccsig3,ccsig4
 real(dp) :: coeffss,cxsig0,cxsig1,cxsig2,cxsig3,cxsig4,d2ecrs0_drs2
 real(dp) :: d2ecrs1_drs2,d2ecrs_drs2,d2ecrs_drsdzeta,d2ecrs_dzeta2
 real(dp) :: d2fzeta4_dzeta2,d2gcrs_drs2,d2macrs_drs2,decrs0_drs,decrs1_drho
 real(dp) :: decrs1_drs,decrs_drs,decrs_dzeta,dfzeta4_dzeta,dgcabdss
 real(dp) :: dgcrs_drs,dgcsigdss,dgxsigdss,divcab,divcsig,divx,dmacrs_drs
 real(dp) :: drhoecab_drhodn,drhoecab_drhoup,drhoecrs1_drhodn,drhoecrs1_drhoup
 real(dp) :: drsdrho,dssdg,dssdndg,dssdndrho,dssdrho,dssupdg,dssupdrho,ducabdss
 real(dp) :: ducsigdss,duxsigdss,ec0_a1,ec0_aa,ec0_b1,ec0_b2,ec0_b3,ec0_b4
 real(dp) :: ec0_den,ec0_log,ec0_q0,ec0_q1,ec0_q1p,ec0_q1pp,ec1_a1,ec1_aa
 real(dp) :: ec1_b1,ec1_b2,ec1_b3,ec1_b4,ec1_den,ec1_log,ec1_q0,ec1_q1,ec1_q1p
 real(dp) :: ec1_q1pp,ecrs,ecrs0,ecrs1,ex_lsd,exc,f_zeta,factfpp_zeta
 real(dp) :: fp_zeta,fpp_zeta,gcab,gcrs,gcsig,gxsig,mac_a1,mac_aa,mac_b1
 real(dp) :: mac_b2,mac_b3,mac_b4,mac_den,mac_log,mac_q0,mac_q1,mac_q1p
 real(dp) :: mac_q1pp,macrs,rho,rho_inv
 real(dp) :: rhoecab,rhoecrs1_dn,rhoecrs1_up,rhomo6,rhomot,rhoo6
 real(dp) :: rhotmo6,rhotmot,rhoto6,rhotot,rhotot_inv,rs,rsm1_2,sqr_rs,ss,ss_dn
 real(dp) :: ss_up,ssavg,ucab,ucsig,uxsig,vxcadd,zeta,zeta4,zetm_1_3
 real(dp) :: zetp_1_3
 character(len=500) :: message
!arrays
 real(dp),allocatable :: rho_updnm1_3(:,:),rhoarr(:),rhom1_3(:),zetm(:)
 real(dp),allocatable :: zetmm1_3(:),zetp(:),zetpm1_3(:)
!no_abirules
!real(dp) :: delta,factor,grr,rho_dn,rho_dnm,rho_dnp,rho_up,rho_upm,rho_upp,zeta_mean

! *************************************************************************

!DEBUG
!write(std_out,*)' xchcth : enter'
!write(std_out,*)' nspden=',nspden
!ENDDEBUG

 if (order/=1) then
   write(message, '(a,i0)' )' Order must be 1 ; argument was ',order
   MSG_BUG(message)
 end if

 if(initialized==0)then
   twom1_3=two**(-third)
   sixpi2_1_3=(six*pi**2)**third
   sixpi2m1_3=one/sixpi2_1_3
   threefourth_divpi=three_quarters*piinv
   factf_zeta= one / ( two**(four/three)-two )
   factfp_zeta= four_thirds * factf_zeta * alpha_zeta2
   sq_rsfac=sqrt(rsfac)
   sq_rsfac_inv=one/sq_rsfac
   initialized=1
 end if

 if(ixc==16)then
!  HCTH-93 parameters from table II of JCP 109, 6264 (1998). Note that there is
!  an error in tables III of JCP 112, 1670 (2000) and JCP 114, 5497 (2001) for the coefficient ccab1
   cxsig0=1.09320_dp ; cxsig1=-0.744056_dp ; cxsig2=5.59920_dp   ; cxsig3=-6.78549_dp ; cxsig4=4.49357_dp
   ccsig0=0.222601_dp; ccsig1=-0.0338622_dp; ccsig2=-0.0125170_dp; ccsig3=-0.802496_dp; ccsig4=1.55396_dp
   ccab0=0.729974_dp ; ccab1=3.35287_dp   ; ccab2=-11.543_dp    ; ccab3=8.08564_dp   ; ccab4=-4.47857_dp
 else if(ixc==17)then
!  HCTH-120 parameters from table III of JCP 112, 1670 (2000)
!  Note the correction of the sign of cxsig1 and ccsig1, as there is a misprint in the HCTH paper !
!  see the exchange of mail with MMarques 23 Dec 2008, and 6 Jan 2009.
!  Now, the functional agrees with the libXC results
   cxsig0=1.09163_dp; cxsig1=-0.7472_dp; cxsig2=5.0783_dp; cxsig3=-4.1075_dp; cxsig4=1.1717_dp
   ccsig0=0.48951_dp; ccsig1=-0.2607_dp; ccsig2=0.4329_dp; ccsig3=-1.9925_dp; ccsig4=2.4853_dp
   ccab0=0.51473_dp ; ccab1=6.9298_dp  ; ccab2=-24.707_dp; ccab3=23.110_dp  ; ccab4=-11.323_dp
!  Exactly same values than in lib_xc
   cxsig0=1.09163_dp; cxsig1=-0.747215_dp; cxsig2=5.07833_dp; cxsig3=-4.10746_dp; cxsig4=1.17173_dp
   ccsig0=0.489508_dp; ccsig1=-0.260699_dp; ccsig2=0.432917_dp; ccsig3=-1.99247_dp; ccsig4=2.48531_dp
   ccab0=0.51473_dp ; ccab1=6.92982_dp  ; ccab2=-24.7073_dp; ccab3=23.1098_dp  ; ccab4=-11.3234_dp
 else if(ixc==26)then
!  HCTH-147 parameters from table III of JCP 112, 1670 (2000)
   cxsig0=1.09025_dp; cxsig1=-0.7992_dp; cxsig2=5.5721_dp ; cxsig3=-5.8676_dp; cxsig4=3.0454_dp
   ccsig0=0.56258_dp; ccsig1=-0.0171_dp; ccsig2=-1.3064_dp; ccsig3= 1.0575_dp; ccsig4=0.8854_dp
   ccab0=0.54235_dp ; ccab1=7.0146_dp  ; ccab2=-28.382_dp ; ccab3=35.033_dp  ; ccab4=-20.428_dp
!  Exactly same values than in lib_xc
   cxsig0=1.09025_dp; cxsig1=-0.799194_dp; cxsig2=5.57212_dp ; cxsig3=-5.86760_dp; cxsig4=3.04544_dp
   ccsig0=0.562576_dp; ccsig1= 0.0171436_dp; ccsig2=-1.30636_dp; ccsig3= 1.05747_dp; ccsig4=0.885429_dp
   ccab0=0.542352_dp ; ccab1=7.01464_dp  ; ccab2=-28.3822_dp ; ccab3=35.0329_dp  ; ccab4=-20.4284_dp
 else if(ixc==27)then
!  HCTH-407 parameters from table IV of JCP 114, 5497 (2001)
   cxsig0=1.08184_dp; cxsig1=-0.5183_dp; cxsig2=3.4256_dp; cxsig3=-2.6290_dp; cxsig4=2.2886_dp
   ccsig0=1.18777_dp; ccsig1=-2.4029_dp; ccsig2=5.6174_dp; ccsig3=-9.1792_dp; ccsig4=6.2480_dp
   ccab0=0.58908_dp ; ccab1=4.4237_dp  ; ccab2=-19.222_dp; ccab3=42.572_dp  ; ccab4=-42.005_dp
!  Exactly same values than in lib_xc
   cxsig0=1.08184_dp; cxsig1=-0.518339_dp; cxsig2=3.42562_dp; cxsig3=-2.62901_dp; cxsig4=2.28855_dp
   ccsig0=1.18777_dp; ccsig1=-2.40292_dp; ccsig2=5.61741_dp; ccsig3=-9.17923_dp; ccsig4=6.24798_dp
   ccab0=0.589076_dp ; ccab1=4.42374_dp  ; ccab2=-19.2218_dp; ccab3=42.5721_dp  ; ccab4=-42.0052_dp
 else
   write(message, '(a,i0)' )' xchcth : ixc must be 16, 17, 26, or 27 ; argument was ',ixc
   MSG_BUG(message)
 end if

!Parameters for the Perdew-Wang 92 LSD as well as LSD-RPA,
!see Table I of Phys.Rev.B 45,13244 (1992) [[cite:Perdew1992a]]
 ec0_aa=0.031091_dp ; ec1_aa=0.015545_dp ; mac_aa=0.016887_dp
 ec0_a1=0.21370_dp  ; ec1_a1=0.20548_dp  ; mac_a1=0.11125_dp
 ec0_b1=7.5957_dp   ; ec1_b1=14.1189_dp  ; mac_b1=10.357_dp
 ec0_b2=3.5876_dp   ; ec1_b2=6.1977_dp   ; mac_b2=3.6231_dp
 ec0_b3=1.6382_dp   ; ec1_b3=3.3662_dp   ; mac_b3=0.88026_dp
 ec0_b4=0.49294_dp  ; ec1_b4=0.62517_dp  ; mac_b4=0.49671_dp

!DEBUG
!Finite-difference debugging, do not take away
!Note : here work with collinear gradients. Might be generalized ...
!debug=2  ! Choose 1 (rho grads) or 2 (grho grads)
!factor=one
!zeta_mean=0.98_dp
!zeta_mean=zero
!delta=0.000025*factor
!delta=0.0000125*factor
!if(debug/=0)then
!do ipts=1,npts,5
!rho=ipts*0.zero*factor
!rho_up=rho*(one+zeta_mean)*half
!rho_dn=rho*(one-zeta_mean)*half
!rho_upp=rho_up+delta
!rho_upm=rho_up-delta
!rho_dnp=rho_dn+delta
!rho_dnm=rho_dn-delta
!! Here, vary rho
!if(debug==1)then
!rho_updn(ipts  ,1)=rho_up ; rho_updn(ipts  ,2)=rho_dn
!rho_updn(ipts+1,1)=rho_upp; rho_updn(ipts+1,2)=rho_dn
!rho_updn(ipts+2,1)=rho_upm; rho_updn(ipts+2,2)=rho_dn
!rho_updn(ipts+3,1)=rho_up ; rho_updn(ipts+3,2)=rho_dnp
!rho_updn(ipts+4,1)=rho_up ; rho_updn(ipts+4,2)=rho_dnm
!grho2_updn(ipts:ipts+4,1)=(0.2_dp*factor)**2     ! grad2 of spin up density
!grho2_updn(ipts:ipts+4,2)=(0.2_dp*factor)**2     ! grad2 of spin down density
!grho2_updn(ipts:ipts+4,3)=(0.3_dp*factor)**2     ! grad2 of total density
!else
!!  Here, vary grho (interchange rho and grho)
!grho2_updn(ipts  ,1)=rho_up**2 ; grho2_updn(ipts  ,2)=rho_dn**2
!grho2_updn(ipts+1,1)=rho_upp**2; grho2_updn(ipts+1,2)=rho_dn**2
!grho2_updn(ipts+2,1)=rho_upm**2; grho2_updn(ipts+2,2)=rho_dn**2
!grho2_updn(ipts+3,1)=rho_up**2 ; grho2_updn(ipts+3,2)=rho_dnp**2
!grho2_updn(ipts+4,1)=rho_up**2 ; grho2_updn(ipts+4,2)=rho_dnm**2
!grho2_updn(ipts  ,3)=(ipts*0.zero*factor)**2
!grho2_updn(ipts+1,3)=(ipts*0.zero*factor+delta)**2
!grho2_updn(ipts+2,3)=(ipts*0.zero*factor-delta)**2
!grho2_updn(ipts+3,3)=(ipts*0.zero*factor+delta)**2   ! identical to ipts+1
!grho2_updn(ipts+4,3)=(ipts*0.zero*factor-delta)**2   ! identical to ipts+2
!rho_updn(ipts:ipts+4,1)=0.2_dp*factor*(one+zeta_mean)*half    ! spin up density
!rho_updn(ipts:ipts+4,2)=0.2_dp*factor*(one-zeta_mean)*half    ! spin down density
!end if
!end do
!end if
!Usual option :
!nspden=2 ; order=2
!GGA
!nspden=2 ; order=1
!Might take also, although finite difference later is meaningless
!nspden=1 ; order=-2
!ENDDEBUG

 if(order**2 >1)then
   factfpp_zeta= third * factfp_zeta * alpha_zeta2
 end if

 ABI_ALLOCATE(rhoarr,(npts))
 ABI_ALLOCATE(rhom1_3,(npts))
 ABI_ALLOCATE(rho_updnm1_3,(npts,2))
 do ispden=1,nspden
   call invcb(rho_updn(:,ispden),rho_updnm1_3(:,ispden),npts)
 end do
 if(nspden==1)then
   rhoarr(:)=two*rho_updn(:,1)
   rhom1_3(:)=twom1_3*rho_updnm1_3(:,1)
   rho_updnm1_3(:,2)=rho_updnm1_3(:,1)
 else
   rhoarr(:)=rho_updn(:,1)+rho_updn(:,2)
   call invcb(rhoarr,rhom1_3,npts)
   ABI_ALLOCATE(zetm,(npts))
   ABI_ALLOCATE(zetmm1_3,(npts))
   ABI_ALLOCATE(zetp,(npts))
   ABI_ALLOCATE(zetpm1_3,(npts))
   do ipts=1,npts
     rhotmot=rhom1_3(ipts)
     rhotot_inv=rhotmot*rhotmot*rhotmot
     zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
     zetp(ipts)=one+zeta*alpha_zeta
     zetm(ipts)=one-zeta*alpha_zeta
   end do
   call invcb(zetp,zetpm1_3,npts)
   call invcb(zetm,zetmm1_3,npts)
 end if

 if (nspden==1) then

   if(order==-2) then

     do ipts=1,npts

!      -----------------------------------------------------------------------
!      First take care of the spin-split part of the functional
       exc=zero
       ispden=1
       rho   =rho_updn(ipts,ispden)
       rhomot=rho_updnm1_3(ipts,ispden)
       rho_inv=rhomot*rhomot*rhomot

!      Exchange part
       ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!      Note that this definition differs from the PBE one
       coeffss=rho_inv*rho_inv*rhomot*rhomot
       ss=grho2_updn(ipts,ispden)*coeffss
       dssdrho=-eight*third*ss*rho_inv
       dssdg=two*coeffss

       divx=one/(one+gammax*ss)
       uxsig=gammax*ss*divx
       duxsigdss=gammax*divx*(one-ss*gammax*divx)

       gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
       dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
&       *duxsigdss

       exc=exc+ex_lsd*rho*gxsig
       vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
       dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

!      Spin parallel correlation part
!      Note that this definition is for fully spin-polarized quantities
       rs=rsfac*rhomot
       rhomo6=sqrt(rhomot)
       sqr_rs=sq_rsfac*rhomo6
       rhoo6=rho*rhomot*rhomot*rhomo6
       rsm1_2=sq_rsfac_inv*rhoo6
       drsdrho=-third*rs*rho_inv

       ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
       ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
       ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
       ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!      ec1_log=log( one + one / ec1_q1 )
       ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
       ecrs1=ec1_q0*ec1_log
       decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
       decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

!      Store the LSDA corr energy and density derivative
       rhoecrs1_up=rho*ecrs1
       drhoecrs1_drhoup=decrs1_drho
       ss_up=ss
       dssupdrho=dssdrho
       dssupdg=dssdg

       divcsig=one/(one+gammacsig*ss)
       ucsig=gammacsig*ss*divcsig
       ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

       gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
       dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
&       *ducsigdss

!      DEBUG
!      gcsig=zero
!      dgcsigdss=zero
!      ENDDEBUG
       exc=exc+ecrs1*rho*gcsig
       vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
&       ecrs1*rho*dgcsigdss*dssdrho
       dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

       rhoecrs1_dn=rhoecrs1_up
       drhoecrs1_drhodn=drhoecrs1_drhoup
       ss_dn=ss_up
       dssdndrho=dssupdrho
       dssdndg=dssupdg
       exc=exc*2

!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

!      From now, the coding of the PW92 functional start. It is identical in xcpbe.f
       rs=rsfac*rhotmot
       sqr_rs=sq_rsfac*rhotmo6
       rsm1_2=sq_rsfac_inv*rhoto6

!      Formulas A6-A8 of PW92LSD
       ec0_q0=-2.0d0*ec0_aa*(one+ec0_a1*rs)
       ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
       ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
       ec0_den=one/(ec0_q1*ec0_q1+ec0_q1)
!      ec0_log=log( one + one / ec0_q1 )
       ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
       ecrs0=ec0_q0*ec0_log
       decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

       ec0_q1pp=half*ec0_aa*(-ec0_b1*rsm1_2**3+3.d0*ec0_b3*rsm1_2+8.d0*ec0_b4)
       d2ecrs0_drs2= 4.0d0*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&       -ec0_q0*ec0_q1pp*ec0_den                        &
&       +ec0_q0*ec0_q1p**2*ec0_den**2*(2.d0*ec0_q1+one)



       mac_q0=-2.0d0*mac_aa*(one+mac_a1*rs)
       mac_q1=2.0d0*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
       mac_q1p=mac_aa*(mac_b1*rsm1_2+2.d0*mac_b2+3.d0*mac_b3*sqr_rs+4.d0*mac_b4*rs)
       mac_den=one/(mac_q1*mac_q1+mac_q1)
       mac_log=-log( mac_q1*mac_q1*mac_den )
       macrs=mac_q0*mac_log
       dmacrs_drs= -2.0d0*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den


       ecrs=ecrs0
       decrs_drs=decrs0_drs
       decrs_dzeta=zero

       d2ecrs_drs2=d2ecrs0_drs2

       d2ecrs_dzeta2=alpha_zeta**2*(-macrs)

       zeta=zero

!      At this point, the coding of the PW92 functional finishes.

!      The correlation between different spin from HCTH is now computed
!      First, the part without gradient correction factor
       rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
       vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
       drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
       drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

!      Now, the gradient correction factor
       ssavg=half*(ss_up+ss_dn)
       divcab=one/(one+gammacab*ssavg)
       ucab=gammacab*ssavg*divcab
       ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

       gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
       dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
&       *ducabdss

       exc=exc+rhoecab*gcab

       vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
       dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
!      If non spin-polarized, treat spin down contribution now, similar to spin up
!      vxci(ipts,2)=vxci(ipts,1)
       dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

!      Final division by the total density, to give the energy density
       exci(ipts)=exc*rhotot_inv

     end do ! ipts=1,npts

   else if(order**2>1) then

     do ipts=1,npts

!      -----------------------------------------------------------------------
!      First take care of the spin-split part of the functional
       exc=zero
       ispden=1
       rho   =rho_updn(ipts,ispden)
       rhomot=rho_updnm1_3(ipts,ispden)
       rho_inv=rhomot*rhomot*rhomot

!      Exchange part
       ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!      Note that this definition differs from the PBE one
       coeffss=rho_inv*rho_inv*rhomot*rhomot
       ss=grho2_updn(ipts,ispden)*coeffss
       dssdrho=-eight*third*ss*rho_inv
       dssdg=two*coeffss

       divx=one/(one+gammax*ss)
       uxsig=gammax*ss*divx
       duxsigdss=gammax*divx*(one-ss*gammax*divx)

       gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
       dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
&       *duxsigdss

       exc=exc+ex_lsd*rho*gxsig
       vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
       dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

!      Spin parallel correlation part
!      Note that this definition is for fully spin-polarized quantities
       rs=rsfac*rhomot
       rhomo6=sqrt(rhomot)
       sqr_rs=sq_rsfac*rhomo6
       rhoo6=rho*rhomot*rhomot*rhomo6
       rsm1_2=sq_rsfac_inv*rhoo6
       drsdrho=-third*rs*rho_inv

       ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
       ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
       ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
       ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!      ec1_log=log( one + one / ec1_q1 )
       ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
       ecrs1=ec1_q0*ec1_log
       decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
       decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

!      Store the LSDA corr energy and density derivative
       rhoecrs1_up=rho*ecrs1
       drhoecrs1_drhoup=decrs1_drho
       ss_up=ss
       dssupdrho=dssdrho
       dssupdg=dssdg

       divcsig=one/(one+gammacsig*ss)
       ucsig=gammacsig*ss*divcsig
       ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

       gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
       dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
&       *ducsigdss

!      DEBUG
!      gcsig=zero
!      dgcsigdss=zero
!      ENDDEBUG
       exc=exc+ecrs1*rho*gcsig
       vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
&       ecrs1*rho*dgcsigdss*dssdrho
       dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

       rhoecrs1_dn=rhoecrs1_up
       drhoecrs1_drhodn=drhoecrs1_drhoup
       ss_dn=ss_up
       dssdndrho=dssupdrho
       dssdndg=dssupdg
       exc=exc*2

!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

!      From now, the coding of the PW92 functional start. It is identical in xcpbe.f
       rs=rsfac*rhotmot
       sqr_rs=sq_rsfac*rhotmo6
       rsm1_2=sq_rsfac_inv*rhoto6

!      Formulas A6-A8 of PW92LSD
       ec0_q0=-2.0d0*ec0_aa*(one+ec0_a1*rs)
       ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
       ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
       ec0_den=one/(ec0_q1*ec0_q1+ec0_q1)
!      ec0_log=log( one + one / ec0_q1 )
       ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
       ecrs0=ec0_q0*ec0_log
       decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
       ec0_q1pp=half*ec0_aa*(-ec0_b1*rsm1_2**3+3.d0*ec0_b3*rsm1_2+8.d0*ec0_b4)
       d2ecrs0_drs2= 4.0d0*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&       -ec0_q0*ec0_q1pp*ec0_den                        &
&       +ec0_q0*ec0_q1p**2*ec0_den**2*(2.d0*ec0_q1+one)


       ecrs=ecrs0
       decrs_drs=decrs0_drs
       decrs_dzeta=zero
       d2ecrs_drs2=d2ecrs0_drs2
       zeta=zero

!      At this point, the coding of the PW92 functional finishes.

!      The correlation between different spin from HCTH is now computed
!      First, the part without gradient correction factor
       rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
       vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
       drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
       drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

!      Now, the gradient correction factor
       ssavg=half*(ss_up+ss_dn)
       divcab=one/(one+gammacab*ssavg)
       ucab=gammacab*ssavg*divcab
       ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

       gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
       dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
&       *ducabdss

       exc=exc+rhoecab*gcab

       vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
       dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
!      If non spin-polarized, treat spin down contribution now, similar to spin up
!      vxci(ipts,2)=vxci(ipts,1)
       dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

!      Final division by the total density, to give the energy density
       exci(ipts)=exc*rhotot_inv

     end do ! ipts=1,npts
   else

     do ipts=1,npts

!      -----------------------------------------------------------------------
!      First take care of the spin-split part of the functional
       exc=zero
       ispden=1
       rho   =rho_updn(ipts,ispden)
       rhomot=rho_updnm1_3(ipts,ispden)
       rho_inv=rhomot*rhomot*rhomot

!      Exchange part
       ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!      Note that this definition differs from the PBE one
       coeffss=rho_inv*rho_inv*rhomot*rhomot
       ss=grho2_updn(ipts,ispden)*coeffss
       dssdrho=-eight*third*ss*rho_inv
       dssdg=two*coeffss

       divx=one/(one+gammax*ss)
       uxsig=gammax*ss*divx
       duxsigdss=gammax*divx*(one-ss*gammax*divx)

       gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
       dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
&       *duxsigdss

       exc=exc+ex_lsd*rho*gxsig
       vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
       dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

!      Spin parallel correlation part
!      Note that this definition is for fully spin-polarized quantities
       rs=rsfac*rhomot
       rhomo6=sqrt(rhomot)
       sqr_rs=sq_rsfac*rhomo6
       rhoo6=rho*rhomot*rhomot*rhomo6
       rsm1_2=sq_rsfac_inv*rhoo6
       drsdrho=-third*rs*rho_inv

       ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
       ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
       ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
       ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!      ec1_log=log( one + one / ec1_q1 )
       ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
       ecrs1=ec1_q0*ec1_log
       decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
       decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

!      Store the LSDA corr energy and density derivative
       rhoecrs1_up=rho*ecrs1
       drhoecrs1_drhoup=decrs1_drho
       ss_up=ss
       dssupdrho=dssdrho
       dssupdg=dssdg

       divcsig=one/(one+gammacsig*ss)
       ucsig=gammacsig*ss*divcsig
       ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

       gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
       dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
&       *ducsigdss

!      DEBUG
!      gcsig=zero
!      dgcsigdss=zero
!      ENDDEBUG
       exc=exc+ecrs1*rho*gcsig
       vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
&       ecrs1*rho*dgcsigdss*dssdrho
       dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

       rhoecrs1_dn=rhoecrs1_up
       drhoecrs1_drhodn=drhoecrs1_drhoup
       ss_dn=ss_up
       dssdndrho=dssupdrho
       dssdndg=dssupdg
       exc=exc*2

!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

!      From now, the coding of the PW92 functional start. It is identical in xcpbe.f
       rs=rsfac*rhotmot
       sqr_rs=sq_rsfac*rhotmo6
       rsm1_2=sq_rsfac_inv*rhoto6

!      Formulas A6-A8 of PW92LSD
       ec0_q0=-2.0d0*ec0_aa*(one+ec0_a1*rs)
       ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
       ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
       ec0_den=one/(ec0_q1*ec0_q1+ec0_q1)
!      ec0_log=log( one + one / ec0_q1 )
       ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
       ecrs0=ec0_q0*ec0_log
       decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

       ecrs=ecrs0
       decrs_drs=decrs0_drs
       decrs_dzeta=zero
       zeta=zero

!      At this point, the coding of the PW92 functional finishes.

!      The correlation between different spin from HCTH is now computed
!      First, the part without gradient correction factor
       rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
       vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
       drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
       drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

!      Now, the gradient correction factor
       ssavg=half*(ss_up+ss_dn)
       divcab=one/(one+gammacab*ssavg)
       ucab=gammacab*ssavg*divcab
       ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

       gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
       dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
&       *ducabdss

       exc=exc+rhoecab*gcab

       vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
       dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
!      If non spin-polarized, treat spin down contribution now, similar to spin up
!      vxci(ipts,2)=vxci(ipts,1)
       dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

!      Final division by the total density, to give the energy density
       exci(ipts)=exc*rhotot_inv

     end do ! ipts=1,npts
   end if

 else if(nspden==2) then

   if(order**2>1) then

     do ipts=1,npts

!      -----------------------------------------------------------------------
!      First take care of the spin-split part of the functional
       exc=zero
       do ispden=1,nspden
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         rho_inv=rhomot*rhomot*rhomot

!        Exchange part
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Note that this definition differs from the PBE one
         coeffss=rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         dssdrho=-eight*third*ss*rho_inv
         dssdg=two*coeffss

         divx=one/(one+gammax*ss)
         uxsig=gammax*ss*divx
         duxsigdss=gammax*divx*(one-ss*gammax*divx)

         gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
         dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
&         *duxsigdss

         exc=exc+ex_lsd*rho*gxsig
         vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

!        Spin parallel correlation part
!        Note that this definition is for fully spin-polarized quantities
         rs=rsfac*rhomot
         rhomo6=sqrt(rhomot)
         sqr_rs=sq_rsfac*rhomo6
         rhoo6=rho*rhomot*rhomot*rhomo6
         rsm1_2=sq_rsfac_inv*rhoo6
         drsdrho=-third*rs*rho_inv

         ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
         ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
         ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
         ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!        ec1_log=log( one + one / ec1_q1 )
         ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
         ecrs1=ec1_q0*ec1_log
         decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
         decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

!        Store the LSDA corr energy and density derivative
         if(ispden==1)then
           rhoecrs1_up=rho*ecrs1
           drhoecrs1_drhoup=decrs1_drho
           ss_up=ss
           dssupdrho=dssdrho
           dssupdg=dssdg
         else
           rhoecrs1_dn=rho*ecrs1
           drhoecrs1_drhodn=decrs1_drho
           ss_dn=ss
           dssdndrho=dssdrho
           dssdndg=dssdg
         end if

         divcsig=one/(one+gammacsig*ss)
         ucsig=gammacsig*ss*divcsig
         ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

         gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
         dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
&         *ducsigdss

!        DEBUG
!        gcsig=zero
!        dgcsigdss=zero
!        ENDDEBUG
         exc=exc+ecrs1*rho*gcsig
         vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
&         ecrs1*rho*dgcsigdss*dssdrho
         dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

       end do

!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

!      From now, the coding of the PW92 functional start. It is identical in xcpbe.f
       rs=rsfac*rhotmot
       sqr_rs=sq_rsfac*rhotmo6
       rsm1_2=sq_rsfac_inv*rhoto6

!      Formulas A6-A8 of PW92LSD
       ec0_q0=-2.0d0*ec0_aa*(one+ec0_a1*rs)
       ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
       ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
       ec0_den=one/(ec0_q1*ec0_q1+ec0_q1)
!      ec0_log=log( one + one / ec0_q1 )
       ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
       ecrs0=ec0_q0*ec0_log
       decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
       ec0_q1pp=half*ec0_aa*(-ec0_b1*rsm1_2**3+3.d0*ec0_b3*rsm1_2+8.d0*ec0_b4)
       d2ecrs0_drs2= 4.0d0*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&       -ec0_q0*ec0_q1pp*ec0_den                        &
&       +ec0_q0*ec0_q1p**2*ec0_den**2*(2.d0*ec0_q1+one)

       mac_q0=-2.0d0*mac_aa*(one+mac_a1*rs)
       mac_q1=2.0d0*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
       mac_q1p=mac_aa*(mac_b1*rsm1_2+2.d0*mac_b2+3.d0*mac_b3*sqr_rs+4.d0*mac_b4*rs)
       mac_den=one/(mac_q1*mac_q1+mac_q1)
       mac_log=-log( mac_q1*mac_q1*mac_den )
       macrs=mac_q0*mac_log
       dmacrs_drs= -2.0d0*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

       zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
       ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
       ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
       ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
       ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!      ec1_log=log( one + one / ec1_q1 )
       ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
       ecrs1=ec1_q0*ec1_log
       decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

!      alpha_zeta is introduced in order to remove singularities for fully
!      polarized systems.
       zetp_1_3=(one+zeta*alpha_zeta)*zetpm1_3(ipts)**2
       zetm_1_3=(one-zeta*alpha_zeta)*zetmm1_3(ipts)**2

       f_zeta=( (one+zeta*alpha_zeta2)*zetp_1_3 +                      &
&       (one-zeta*alpha_zeta2)*zetm_1_3 - 2.0d0 ) * factf_zeta
       fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
       zeta4=zeta**4
       gcrs=ecrs1-ecrs0+macrs*fsec_inv
!      ecrs=ecrs0+f_zeta*(-macrs*(one-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
       ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

       dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
!      decrs_drs=decrs0_drs+f_zeta*&
!      &       (-dmacrs_drs*(one-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
       decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
       dfzeta4_dzeta=4.0d0*zeta**3*f_zeta+fp_zeta*zeta4
       decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

       ec1_q1pp=half*ec1_aa*(-ec1_b1*rsm1_2**3+3.d0*ec1_b3*rsm1_2+8.d0*ec1_b4)
       d2ecrs1_drs2= 4.0d0*ec1_aa*ec1_a1*ec1_q1p*ec1_den            &
&       -ec1_q0*ec1_q1pp*ec1_den                        &
&       +ec1_q0*ec1_q1p**2*ec1_den**2*(2.d0*ec1_q1+one)

       mac_q1pp=half*mac_aa*(-mac_b1*rsm1_2**3+3.d0*mac_b3*rsm1_2+8.d0*mac_b4)
       d2macrs_drs2= 4.0d0*mac_aa*mac_a1*mac_q1p*mac_den            &
&       -mac_q0*mac_q1pp*mac_den                        &
&       +mac_q0*mac_q1p**2*mac_den**2*(2.d0*mac_q1+one)

       d2gcrs_drs2=d2ecrs1_drs2-d2ecrs0_drs2+d2macrs_drs2*fsec_inv
       fpp_zeta=(zetpm1_3(ipts)**2+zetmm1_3(ipts)**2) * factfpp_zeta
       d2fzeta4_dzeta2=12.0d0*zeta**2*f_zeta  &
&       + 8.0d0*zeta**3*fp_zeta &
&       +       zeta4  *fpp_zeta

       d2ecrs_drs2=d2ecrs0_drs2+&
&       f_zeta*(zeta4*d2gcrs_drs2-d2macrs_drs2*fsec_inv)
       d2ecrs_drsdzeta=dfzeta4_dzeta*dgcrs_drs-fp_zeta*dmacrs_drs*fsec_inv
       d2ecrs_dzeta2=d2fzeta4_dzeta2*gcrs-fpp_zeta*macrs*fsec_inv


!      At this point, the coding of the PW92 functional finishes.

!      The correlation between different spin from HCTH is now computed
!      First, the part without gradient correction factor
       rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
       vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
       drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
       drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

!      Now, the gradient correction factor
       ssavg=half*(ss_up+ss_dn)
       divcab=one/(one+gammacab*ssavg)
       ucab=gammacab*ssavg*divcab
       ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

       gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
       dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
&       *ducabdss

       exc=exc+rhoecab*gcab

       vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
       dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
       vxci(ipts,2)=vxci(ipts,2)+drhoecab_drhodn*gcab+rhoecab*dgcabdss*half*dssdndrho
       dvxcdgr(ipts,2)=dvxcdgr(ipts,2)+rhoecab*dgcabdss*half*dssdndg

!      Final division by the total density, to give the energy density
       exci(ipts)=exc*rhotot_inv

     end do ! ipts=1,npts

   else

     do ipts=1,npts

!      -----------------------------------------------------------------------
!      First take care of the spin-split part of the functional
       exc=zero
       do ispden=1,nspden
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         rho_inv=rhomot*rhomot*rhomot

!        Exchange part
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Note that this definition differs from the PBE one
         coeffss=rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         dssdrho=-eight*third*ss*rho_inv
         dssdg=two*coeffss

         divx=one/(one+gammax*ss)
         uxsig=gammax*ss*divx
         duxsigdss=gammax*divx*(one-ss*gammax*divx)

         gxsig=cxsig0+uxsig*(cxsig1+uxsig*(cxsig2+uxsig*(cxsig3+uxsig*cxsig4)))
         dgxsigdss=(cxsig1+uxsig*(two*cxsig2+uxsig*(three*cxsig3+uxsig*four*cxsig4)))&
&         *duxsigdss

         exc=exc+ex_lsd*rho*gxsig
         vxci(ipts,ispden)=ex_lsd*(four_thirds*gxsig+rho*dgxsigdss*dssdrho)
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dgxsigdss*dssdg

!        Spin parallel correlation part
!        Note that this definition is for fully spin-polarized quantities
         rs=rsfac*rhomot
         rhomo6=sqrt(rhomot)
         sqr_rs=sq_rsfac*rhomo6
         rhoo6=rho*rhomot*rhomot*rhomo6
         rsm1_2=sq_rsfac_inv*rhoo6
         drsdrho=-third*rs*rho_inv

         ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
         ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
         ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
         ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!        ec1_log=log( one + one / ec1_q1 )
         ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
         ecrs1=ec1_q0*ec1_log
         decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den
         decrs1_drho=ecrs1+decrs1_drs*drsdrho*rho

!        Store the LSDA corr energy and density derivative
         if(ispden==1)then
           rhoecrs1_up=rho*ecrs1
           drhoecrs1_drhoup=decrs1_drho
           ss_up=ss
           dssupdrho=dssdrho
           dssupdg=dssdg
         else
           rhoecrs1_dn=rho*ecrs1
           drhoecrs1_drhodn=decrs1_drho
           ss_dn=ss
           dssdndrho=dssdrho
           dssdndg=dssdg
         end if

         divcsig=one/(one+gammacsig*ss)
         ucsig=gammacsig*ss*divcsig
         ducsigdss=gammacsig*divcsig*(one-ss*gammacsig*divcsig)

         gcsig=ccsig0+ucsig*(ccsig1+ucsig*(ccsig2+ucsig*(ccsig3+ucsig*ccsig4)))
         dgcsigdss=(ccsig1+ucsig*(two*ccsig2+ucsig*(three*ccsig3+ucsig*four*ccsig4)))&
&         *ducsigdss

!        DEBUG
!        gcsig=zero
!        dgcsigdss=zero
!        ENDDEBUG
         exc=exc+ecrs1*rho*gcsig
         vxci(ipts,ispden)=vxci(ipts,ispden)+decrs1_drho*gcsig+&
&         ecrs1*rho*dgcsigdss*dssdrho
         dvxcdgr(ipts,ispden)=dvxcdgr(ipts,ispden)+ecrs1*rho*dgcsigdss*dssdg

       end do

!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6

!      From now, the coding of the PW92 functional start. It is identical in xcpbe.f
       rs=rsfac*rhotmot
       sqr_rs=sq_rsfac*rhotmo6
       rsm1_2=sq_rsfac_inv*rhoto6

!      Formulas A6-A8 of PW92LSD
       ec0_q0=-2.0d0*ec0_aa*(one+ec0_a1*rs)
       ec0_q1=2.0d0*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
       ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2.d0*ec0_b2+3.d0*ec0_b3*sqr_rs+4.d0*ec0_b4*rs)
       ec0_den=one/(ec0_q1*ec0_q1+ec0_q1)
!      ec0_log=log( one + one / ec0_q1 )
       ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
       ecrs0=ec0_q0*ec0_log
       decrs0_drs= -2.0d0*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

       mac_q0=-2.0d0*mac_aa*(one+mac_a1*rs)
       mac_q1=2.0d0*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
       mac_q1p=mac_aa*(mac_b1*rsm1_2+2.d0*mac_b2+3.d0*mac_b3*sqr_rs+4.d0*mac_b4*rs)
       mac_den=one/(mac_q1*mac_q1+mac_q1)
       mac_log=-log( mac_q1*mac_q1*mac_den )
       macrs=mac_q0*mac_log
       dmacrs_drs= -2.0d0*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

       zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
       ec1_q0=-two*ec1_aa*(one+ec1_a1*rs)
       ec1_q1=two*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
       ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+two*ec1_b2+three*ec1_b3*sqr_rs+four*ec1_b4*rs)
       ec1_den=one/(ec1_q1*ec1_q1+ec1_q1)
!      ec1_log=log( one + one / ec1_q1 )
       ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
       ecrs1=ec1_q0*ec1_log
       decrs1_drs= -two*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

!      alpha_zeta is introduced in order to remove singularities for fully
!      polarized systems.
       zetp_1_3=(one+zeta*alpha_zeta)*zetpm1_3(ipts)**2
       zetm_1_3=(one-zeta*alpha_zeta)*zetmm1_3(ipts)**2

       f_zeta=( (one+zeta*alpha_zeta2)*zetp_1_3 +                      &
&       (one-zeta*alpha_zeta2)*zetm_1_3 - 2.0d0 ) * factf_zeta
       fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
       zeta4=zeta**4
       gcrs=ecrs1-ecrs0+macrs*fsec_inv
!      ecrs=ecrs0+f_zeta*(-macrs*(one-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
       ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

       dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
!      decrs_drs=decrs0_drs+f_zeta*&
!      &       (-dmacrs_drs*(one-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
       decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
       dfzeta4_dzeta=4.0d0*zeta**3*f_zeta+fp_zeta*zeta4
       decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

!      At this point, the coding of the PW92 functional finishes.

!      The correlation between different spin from HCTH is now computed
!      First, the part without gradient correction factor
       rhoecab=ecrs*rhotot-rhoecrs1_dn-rhoecrs1_up
       vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
       drhoecab_drhoup=vxcadd+decrs_dzeta-drhoecrs1_drhoup
       drhoecab_drhodn=vxcadd-decrs_dzeta-drhoecrs1_drhodn

!      Now, the gradient correction factor
       ssavg=half*(ss_up+ss_dn)
       divcab=one/(one+gammacab*ssavg)
       ucab=gammacab*ssavg*divcab
       ducabdss=gammacab*divcab*(one-ssavg*gammacab*divcab)

       gcab=ccab0+ucab*(ccab1+ucab*(ccab2+ucab*(ccab3+ucab*ccab4)))
       dgcabdss=(ccab1+ucab*(two*ccab2+ucab*(three*ccab3+ucab*four*ccab4)))&
&       *ducabdss

       exc=exc+rhoecab*gcab

       vxci(ipts,1)=vxci(ipts,1)+drhoecab_drhoup*gcab+rhoecab*dgcabdss*half*dssupdrho
       dvxcdgr(ipts,1)=dvxcdgr(ipts,1)+rhoecab*dgcabdss*half*dssupdg
       vxci(ipts,2)=vxci(ipts,2)+drhoecab_drhodn*gcab+rhoecab*dgcabdss*half*dssdndrho
       dvxcdgr(ipts,2)=dvxcdgr(ipts,2)+rhoecab*dgcabdss*half*dssdndg

!      Final division by the total density, to give the energy density
       exci(ipts)=exc*rhotot_inv

     end do ! ipts=1,npts

   end if

 else
!  Disallowed value for nspden
   write(message, '(3a,i0)' )&
&   '  Argument nspden must be 1 or 2; ',ch10,&
&   '  Value provided as argument was ',nspden
   MSG_BUG(message)
 end if

!DEBUG
!Finite-difference debugging, do not take away
!Beware that dvxcdgr(:,3) no longer exists
!if(debug/=0)then
!do ipts=1,npts,5

!rho=rho_updn(ipts,1)+rho_updn(ipts,2)
!write(std_out,'(a,i5,a,es16.8)' ) ' Point number',ipts,' with rho=',rho

!! For rho
!if(debug==1)then
!write(std_out,'(3es16.8)' )exci(ipts)*rho,vxci(ipts,1),vxci(ipts,2)
!else
!!  For grho2
!write(std_out,'(4es16.8)' )exci(ipts)*rho,dvxcdgr(ipts,1),&
!&  dvxcdgr(ipts,2),dvxcdgr(ipts,3)
!end if

!if(debug==1)then
!!  For rho
!write(std_out,'(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*(rho+delta) - exci(ipts+2)*(rho-delta) )/2.d0/delta,&
!&      ( exci(ipts+3)*(rho+delta) - exci(ipts+4)*(rho-delta) )/2.d0/delta
!write(std_out,'(3es16.8)' )&
!&    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2.d0/delta,&
!&    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2.d0/delta,&
!&    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2.d0/delta
!write(std_out,'(4es16.8)' )&
!&    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2.d0/delta,&
!&    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2.d0/delta,&
!&    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2.d0/delta,&
!&    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2.d0/delta
!else
!!  For grho2  (should distinguish exchange and correlation ...)
!grr=sqrt(grho2_updn(ipts,1)) ! Analysis of exchange
!grr=sqrt(grho2_updn(ipts,3)) ! Analysis of correlation
!write(std_out,'(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*rho - exci(ipts+2)*rho )/2.d0/delta/grr,&
!&      ( exci(ipts+3)*rho - exci(ipts+4)*rho )/2.d0/delta/grr
!write(std_out,'(3es16.8)' )&
!&    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2.d0/delta/grr,&
!&    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2.d0/delta/grr,&
!&    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2.d0/delta/grr
!write(std_out,'(4es16.8)' )&
!&    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2.d0/delta/grr,&
!&    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2.d0/delta/grr,&
!&    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2.d0/delta/grr,&
!&    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2.d0/delta/grr
!end if
!end do
!stop
!end if
!ENDDEBUG

 ABI_DEALLOCATE(rhoarr)
 ABI_DEALLOCATE(rhom1_3)
 ABI_DEALLOCATE(rho_updnm1_3)
 if(nspden==2) then
   ABI_DEALLOCATE(zetm)
   ABI_DEALLOCATE(zetmm1_3)
   ABI_DEALLOCATE(zetp)
   ABI_DEALLOCATE(zetpm1_3)
 end if

!DEBUG
!write(std_out,*)' xchcth : exit'
!write(std_out,*)' nspden=',nspden
!if(order==2)stop
!ENDDEBUG

end subroutine xchcth
!!***

end module m_xchcth
!!***

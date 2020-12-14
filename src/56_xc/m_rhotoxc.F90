!!****m* ABINIT/m_rhotoxc
!! NAME
!!  m_rhotox
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MF, GZ, DRH, MT, SPr)
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

module m_rhotoxc

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_errors
 use m_cgtools
 use m_xcdata
 use m_xc_vdw
 use libxc_functionals

 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_geometry,         only : metric
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_xcpositron,       only : xcpositron
 use m_drivexc,          only : size_dvxc,drivexc,xcmult,mkdenpos
 use m_xclda,            only : xctfw
 use m_xctk,             only : xcden, xcpot

 implicit none

 private
!!***

 public :: rhotoxc
!!***

contains
!!***

!!****f* ABINIT/rhotoxc
!! NAME
!! rhotoxc
!!
!! FUNCTION
!! Start from the density or spin-density, and
!! compute xc correlation potential and energies.
!! Eventually compute xc kernel (if option=-2, 2, 3, 10 or 12).
!! Cannot be used with wavelets.
!!
!! INPUTS
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,xcdata%nspden*nhatdim)= -PAW only- compensation density
!!  nhatdim= -PAW only- 0 if nhat array is not used ; 1 otherwise
!!  nhatgr(nfft,xcdata%nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  n3xccc=dimension of the xccc3d array (0 or nfft or cplx*nfft).
!!  option=0 or 1 for xc only (exc, vxc, strsxc),
!!         2 for xc and kxc (no paramagnetic part if xcdata%nspden=1)
!!        10 for xc  and kxc with only partial derivatives wrt density part (d2Exc/drho^2)
!!        12 for xc and kxc with only partial derivatives wrt density part (d2Exc/drho^2)
!!              and, in the case of hybrid functionals, substitution of the hybrid functional
!!              by the related auxiliary GGA functional for the computation of the xc kernel (not for other quantities)
!!         3 for xc, kxc and k3xc
!!        -2 for xc and kxc (with paramagnetic part if xcdata%nspden=1)
!!  rhor(nfft,xcdata%nspden)=electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if xcdata%nspden=2)
!!   (total in first comp. and magnetization in comp. 2 to 4 if xcdata%nspden=4)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  [vhartr(nfft)=Hartree potential (only needed for Fermi-Amaldi functional)]
!!  xcdata <type(xcdata_type)>=storage for different input variables and derived parameters needed to compute the XC functional
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!
!!  === optional inputs ===
!!  [add_tfw]=flag controling the addition of Weiszacker gradient correction to Thomas-Fermi kin energy
!!  [taur(nfftf,xcdata%nspden*xcdata%usekden)]=array for kinetic energy density
!!  [xc_funcs(2)]= <type(libxc_functional_type)>, optional : libxc XC functionals. Must be coherent with xcdata.
!!  [xccctau3d(n3xccc)]=3D core electron kinetic energy density for XC core correction (bohr^-3)
!!
!! OUTPUT
!!  enxc=returned exchange and correlation energy (hartree).
!!  strsxc(6)= contribution of xc to stress tensor (hartree/bohr^3),
!!   given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!   (note: fxc is rho*exc in the following)
!!   Explicitely : strsxc(mu,nu) = (1/N) Sum(i=1,N)
!!    ( delta(mu,nu) * [  exc(i)rhotot(i)
!!               - depsxc_drho(up,i)*rhor(up,i)-depsxc_drho(dn,i)*rhor(dn,i)]
!!     - gradrho(up,mu)*gradrho(up,nu) * depsxc_dgradrho(up,i) / gradrho(up,i)
!!     - gradrho(dn,mu)*gradrho(dn,nu) * depsxc_dgradrho(dn,i) / gradrho(dn,i) )
!!  vxc(nfft,xcdata%nspden)=xc potential
!!    (spin up in first half and spin down in second half if xcdata%nspden=2)
!!    (v^11, v^22, Re[V^12], Im[V^12] if xcdata%nspden=4)
!!  vxcavg=<Vxc>=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!
!!  === Only if abs(option)=2, -2, 3, 10, 12 (in case 12, for hybrids, substitution of the related GGA) ===
!!  kxc(nfft,nkxc)=exchange and correlation kernel (returned only if nkxc/=0)
!!    Content of Kxc array:
!!   ===== if LDA
!!    if xcdata%nspden==1: kxc(:,1)= d2Exc/drho2
!!              that is 1/2 ( d2Exc/drho_up drho_up + d2Exc/drho_up drho_dn )
!!                         kxc(:,2)= d2Exc/drho_up drho_dn
!!    if xcdata%nspden>=2: kxc(:,1)= d2Exc/drho_up drho_up
!!                         kxc(:,2)= d2Exc/drho_up drho_dn
!!                         kxc(:,3)= d2Exc/drho_dn drho_dn
!!   ===== if GGA
!!    if xcdata%nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if xcdata%nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!
!!  === Only if abs(option)=3 ===
!!  [k3xc(nfft,nk3xc)]= -- optional -- third derivative of the XC energy functional of the density,
!!    at each point of the real space grid (only in the LDA or LSDA)
!!    Content of K3xc array:
!!    ===== if LDA
!!    if xcdata%nspden==1: k3xc(:,1)= d3Exc/drho3
!!    if xcdata%nspden>=2, k3xc(:,1)= d3Exc/drho_up drho_up drho_up
!!                  k3xc(:,2)= d3Exc/drho_up drho_up drho_dn
!!                  k3xc(:,3)= d3Exc/drho_up drho_dn drho_dn
!!                  k3xc(:,4)= d3Exc/drho_dn drho_dn drho_dn

!! === Additional optional output ===
!!  [exc_vdw_out]= vdW-DF contribution to enxc (hartree)
!!  [vxctau(nfft,xcdata%nspden,4*xcdata%usekden)]=(only for meta-GGA)=
!!    vxctau(:,:,1): derivative of XC energy density with respect to kinetic energy density (depsxcdtau).
!!    vxctau(:,:,2:4): gradient of vxctau (gvxctau)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>= -- optional argument -- quantities for the electron-positron annihilation
!!
!! NOTES
!! Start from the density, and compute Hartree (if option>=1) and xc correlation potential and energies.
!! Eventually compute xc kernel (if option=-2, 2, 3, 10 or 12 - in the latter case, substitution by the related GGA kernel).
!! Allows a variety of exchange-correlation functionals
!! according to ixc. Here is a list of allowed values.
!!                                                    subroutine name
!!   <0 means use of libxc
!!    0 means no xc applied (usually for testing)
!! *LDA,LSD
!!    1 means new Teter (4/93) with spin-pol option        xcspol
!!    2 means Perdew-Zunger-Ceperley-Alder                 xcpzca
!!    3 means old Teter (4/91) fit to Ceperley-Alder data  xctetr
!!    4 means Wigner                                       xcwign
!!    5 means Hedin-Lundqvist                              xchelu
!!    6 means "X-alpha" xc                                 xcxalp
!!    7 mean Perdew-Wang 92 LSD fit to Ceperley-Alder data xcpbe
!!    8 mean Perdew-Wang 92 LSD , exchange-only            xcpbe
!!    9 mean Perdew-Wang 92 Ex+Ec_RPA  energy              xcpbe
!!   10 means RPA LSD energy (only the energy !!)          xcpbe
!! *GGA
!!   11 means Perdew-Burke-Ernzerhof GGA functional        xcpbe
!!   12 means x-only Perdew-Burke-Ernzerhof GGA functional xcpbe
!!   13 means LDA (ixc==7), except that the xc potential
!!      is given within the van Leeuwen-Baerends GGA       xclb
!!   14 means revPBE GGA functional                        xcpbe
!!   15 means RPBE GGA functional                          xcpbe
!!   16 means HCTH GGA functional                          xchcth
!!   23 means WC GGA functional                            xcpbe
!!   24 means C09x GGA exchange functional                 xcpbe
!! *Fermi-Amaldi
!!   20 means Fermi-Amaldi correction
!!   21 means Fermi-Amaldi correction with LDA(ixc=1) kernel
!!   22 means Fermi-Amaldi correction with hybrid BPG kernel
!! *Hybrid GGA
!!   41 means PBE0-1/4                                     xcpbe
!!   42 means PBE0-1/3                                     xcpbe
!! *Other
!!   50 means IIT xc                                       xciit
!!
!! NOTE: please update echo_xc_name.F90 if you add new functional (apart from libxc)
!!
!! Allow for improved xc quadrature (intxc=1) by using the usual FFT grid
!! as well as another, shifted, grid, and combining both results.
!! Spin-polarization is allowed only with ixc=0, 1, and GGAs until now.
!!
!! To make the variable names easier to understand, a rule notation is tentatively proposed here:
!!   rho ---> means density
!!   tau ---> means kinetic energy density
!!   exc ---> means exchange-correlation energy density per particle
!!   epsxc ---> means rho*exc == exchange-correlation energy density
!!   vxc ---> means exchange-correlation potential
!!   bigexc ---> means exchange-correlation energy E_xc (for the moment it is named "enxc")
!!   m_norm ---> means norm of magnetization
!!
!!   g... --> means gradient of something (e.g. : grho --> means gradient of electron density)
!!   g...2 -> means square norm of gradient of something (e.g. : grho2 -> means square norm of gradient of electron density)
!!   l... --> means laplacian of something (e.g. : lrho --> means laplacian of electron density)
!!   d...d... --> means derivative of something with regards to something else.
!!   (d2...d...d...  ---> means second derivative of ... with regards to ... and to ...) etc...
!!   d... --> without the occurence of the second "d" means that this is an array of
!!            several derivative of the same quantity (e.g. : depsxc)
!!
!!   ..._b ----> means a block of the quantity "..." (use in mpi loops which treat the data block by block)
!!   ..._updn -> means that spin up and spin down is available in that array
!!               as (..,1) and (..,2). (if xcdata%nspden >=2 of course).
!!   ..._apn --> in case of positrons are concerned.
!!
!!   for more details about notations please see pdf in /doc/theory/MGGA/
!!
!! PARENTS
!!      m_dft_energy,m_forstr,m_kxc,m_longwave,m_nonlinear,m_odamix,m_prcref
!!      m_respfn_driver,m_rhotov,m_scfcv_core,m_setvtr,m_vhxc_me,m_xchybrid
!!
!! CHILDREN
!!      dotprod_vn,drivexc,libxc_functionals_end
!!      libxc_functionals_get_hybridparams,libxc_functionals_init,mean_fftr
!!      metric,mkdenpos,size_dvxc,timab,xc_vdw_aggregate,xcden,xcmult
!!      xcpositron,xcpot,xctfw,xmpi_sum
!!
!! SOURCE

subroutine rhotoxc(enxc,kxc,mpi_enreg,nfft,ngfft, &
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,option, &
& rhor,rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata, &
& add_tfw,exc_vdw_out,electronpositron,k3xc,taur,vhartr,vxctau,xc_funcs,xcctau3d) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nk3xc,n3xccc,nfft,nhatdim,nhatgrdim,nkxc,option
 integer,intent(in) :: usexcnhat
 logical,intent(in) :: non_magnetic_xc
 logical,intent(in),optional :: add_tfw
 real(dp),intent(out) :: enxc,vxcavg
 real(dp),intent(out),optional :: exc_vdw_out
 type(MPI_type),intent(in) :: mpi_enreg
 type(electronpositron_type),pointer,optional :: electronpositron
 type(xcdata_type), intent(in) :: xcdata
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat(nfft,xcdata%nspden*nhatdim)
 real(dp),intent(in) :: nhatgr(nfft,xcdata%nspden,3*nhatgrdim)
 real(dp),intent(in),target :: rhor(nfft,xcdata%nspden)
 real(dp),intent(in) :: rprimd(3,3),xccc3d(n3xccc)
 real(dp),intent(in),optional :: xcctau3d(:)
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vxc(nfft,xcdata%nspden)
 real(dp),intent(in),optional :: vhartr(nfft)
 real(dp),intent(in),target,optional :: taur(:,:)
 real(dp),intent(out),optional :: k3xc(1:nfft,1:nk3xc),vxctau(:,:,:)
 type(libxc_functional_type),intent(inout),optional :: xc_funcs(2)

!Local variables-------------------------------
!scalars
 integer :: auxc_ixc,cplex,ierr,ifft,ii,ixc,ixc_from_lib,indx,ipositron,ipts,ishift,ispden,iwarn,iwarnp
 integer :: jj,mpts,ndvxc,nd2vxc,nfftot,ngr,ngrad,ngrad_apn,nkxc_eff,npts
 integer :: nspden,nspden_apn,nspden_eff,nspden_updn,nspgrad,nvxcgrho,nvxclrho,nvxctau
 integer :: n3xctau,order,usefxc,nproc_fft,comm_fft,usegradient,usekden,uselaplacian
 logical :: my_add_tfw
 real(dp),parameter :: mot=-one/3.0_dp
 real(dp) :: coeff,divshft,doti,dstrsxc,dvdn,dvdz,epsxc,exc_str,factor,m_norm_min,s1,s2,s3
 real(dp) :: strdiag,strsxc1_tot,strsxc2_tot,strsxc3_tot,strsxc4_tot
 real(dp) :: strsxc5_tot,strsxc6_tot,ucvol
 logical :: test_nhat,need_nhat,need_nhatgr,with_vxctau
 character(len=500) :: message
 real(dp) :: hyb_mixing, hyb_mixing_sr, hyb_range
!arrays
 real(dp) :: gm_norm(3),grho(3),gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3)
 real(dp) :: tsec(2),vxcmean(4)
 real(dp),allocatable :: d2vxc_b(:,:),depsxc(:,:),depsxc_apn(:,:),dvxc_apn(:),dvxc_b(:,:)
 real(dp),allocatable :: exc_b(:),fxc_b(:),fxc_apn(:),grho2_apn(:),grho2_b_updn(:,:),lrhonow(:,:),lrho_b_updn(:,:)
 real(dp),allocatable :: m_norm(:),nhat_up(:),rho_b_updn(:,:),rho_b(:),rhonow_apn(:,:,:)
 real(dp),allocatable :: tau_b_updn(:,:),vxc_apn(:,:),vxcgr_apn(:),vxcgrho_b_updn(:,:),vxcrho_b_updn(:,:)
 real(dp),allocatable :: vxc_b_apn(:),vxc_ep(:),vxctau_b_updn(:,:),vxclrho_b_updn(:,:)
 real(dp),allocatable,target :: rhonow(:,:,:),taunow(:,:,:)
 real(dp),pointer :: rhocorval(:,:),rhor_(:,:),taucorval(:,:),taur_(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: rhonow_ptr(:,:,:)
 real(dp) :: deltae_vdw,exc_vdw
 real(dp) :: decdrho_vdw(xcdata%nspden),decdgrho_vdw(3,xcdata%nspden)
 real(dp) :: strsxc_vdw(3,3)
 type(libxc_functional_type) :: xc_funcs_auxc(2)

! *************************************************************************

! Note: the following cases seem to never be tested (should be fixed)
!      - ipositron==2 and ngrad_apn==2
!      - usewvl/=0
!      - test_nhat and usexcnhat==1 and nspden==4

 call timab(81,1,tsec)

!Optional arguments
 my_add_tfw=.false.;if (present(add_tfw)) my_add_tfw=add_tfw

!Useful scalars
 nspden=xcdata%nspden
 ixc=xcdata%ixc
 auxc_ixc=xcdata%auxc_ixc
 n3xctau=0

!nspden_updn: 1 for non-polarized, 2 for polarized
 nspden_updn=min(nspden,2)

!The variable order indicates to which derivative of the energy
!the computation must be done. Computing exc and vxc needs order=1 .
!Meaningful values are 1, 2, 3. Lower than 1 is the same as 1, and larger
!than 3 is the same as 3.
!order=1 or 2 supported for all LSD and GGA ixc
!order=3 supported only for ixc=3 and ixc=7
 order=1
 if(option==2.or.option==10.or.option==12)order=2
 if(option==-2)order=-2
 if(option==3)order=3

!Sizes of local arrays
 if (present(xc_funcs)) then
   call size_dvxc(ixc,order,nspden_updn,&
&            usegradient=usegradient,uselaplacian=uselaplacian,usekden=usekden,&
&            nvxcgrho=nvxcgrho,nvxclrho=nvxclrho,nvxctau=nvxctau,&
&            ndvxc=ndvxc,nd2vxc=nd2vxc,add_tfw=my_add_tfw,xc_funcs=xc_funcs)
 else
   call size_dvxc(ixc,order,nspden_updn,&
&            usegradient=usegradient,uselaplacian=uselaplacian,usekden=usekden,&
&            nvxcgrho=nvxcgrho,nvxclrho=nvxclrho,nvxctau=nvxctau,&
&            ndvxc=ndvxc,nd2vxc=nd2vxc,add_tfw=my_add_tfw)
 end if

!ngrad=1 is for LDAs or LSDs, ngrad=2 is for GGAs/mGGAs
 ngrad=1;if(xcdata%xclevel==2.or.usegradient==1) ngrad=2

!nspden_eff: effective value of nspden used to compute gradients of density:
!  1 for non-polarized system,
!  2 for collinear polarized system or LDA (can be reduced to a collinear system)
!  4 for non-collinear polarized system and GGA
 nspden_eff=nspden_updn;if (nspden==4.and.ngrad==2) nspden_eff=4

!Number of kcxc components depends on option (force LDA type if option==10 or 12)
 nkxc_eff=nkxc;if (option==10.or.option==12) nkxc_eff=min(nkxc,3)

!Check options
 if(option==3.and.nd2vxc==0.and.ixc/=0)then
   write(message, '(3a,i0)' )&
&   'Third-order xc kernel can only be computed for ixc = 0, 3, 7 or 8,',ch10,&
&   'while it is found to be ',ixc
   MSG_ERROR(message)
 end if
 if(nspden==4.and.xcdata%xclevel==2.and.(abs(option)==2))then
   MSG_BUG('When nspden==4 and GGA, the absolute value of option cannot be 2 !')
 end if
 if(ixc<0) then
   if (present(xc_funcs)) then
     ixc_from_lib=libxc_functionals_ixc(xc_functionals=xc_funcs)
   else
     ixc_from_lib=libxc_functionals_ixc()
   end if
!  Check consistency between ixc passed in input and the one used to initialize the library.
   if (ixc /= ixc_from_lib) then
     write(message, '(a,i0,2a,i0,2a)')&
&     'The value of ixc specified in input, ixc = ',ixc,ch10,&
&     'differs from the one used to initialize the functional ',ixc_from_lib,ch10,&
&     'Action: reinitialize the global structure funcs, see NOTES in m_libxc_functionals'
     MSG_BUG(message)
   end if
 end if

!Handling of mGGA functionals
 with_vxctau=(present(vxctau))
 if (with_vxctau) with_vxctau=(size(vxctau)>0)
 if (usekden==1) then
   if (.not.present(taur)) then
     message=' For mGGA functionals, kinetic energy density is needed. Set input variable usekden to 1.' 
     message=trim(message)//' Also use NC pseudopotentials without non-linear XC core correction.'
     MSG_BUG(message)
   else if (size(taur)/=nfft*nspden) then
     message=' Invalid size for taur!'
     MSG_BUG(message)
   end if
   if (present(xcctau3d)) then
     n3xctau=size(xcctau3d)
     if (n3xctau/=0.and.n3xctau/=nfft) then
       message=' Invalid size for xccctau3d!'
       MSG_BUG(message)
     end if
   end if
   if (with_vxctau) then
     if (size(vxctau)/=nfft*nspden*4) then
       message=' Invalid size for vxctau!'
       MSG_BUG(message)
     end if
   end if
 end if
 if((usekden==1.or.uselaplacian==1).and.nspden==4)then
   !mGGA en NC-magnetism: how do we rotate tau kinetic energy density?
   message=' At present, meta-GGA (usekden=1 or uselaplacian=1)  is not comptatible with non-collinear magnetism (nspden=4).'
   MSG_ERROR(message)
 end if

!MPI FFT communicator
 comm_fft = mpi_enreg%comm_fft; nproc_fft = mpi_enreg%nproc_fft

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!In this routine, hartre, xcden and xcpot are called for real
!densities and potentials, corresponding to zero wavevector
 cplex=1
 qphon(:)=zero
 iwarn=0
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 usefxc=0;if (ixc==50) usefxc=1

!Initializations
 enxc=zero
 epsxc=zero
 vxc(:,:)=zero
 vxcavg=zero
 strsxc(:)=zero
 strsxc1_tot=zero;strsxc2_tot=zero;strsxc3_tot=zero
 strsxc4_tot=zero;strsxc5_tot=zero;strsxc6_tot=zero
 if (with_vxctau) vxctau(:,:,:)=zero
 if (nkxc/=0) kxc(:,:)=zero
 if(abs(option)==3.and.nk3xc/=0) k3xc(:,:)=zero
 ipositron=0
 if (present(electronpositron)) then
   ipositron=electronpositron_calctype(electronpositron)
   if (ipositron==2) then
     electronpositron%e_xc  =zero
     electronpositron%e_xcdc=zero
   end if
 end if
 deltae_vdw = zero
 exc_vdw = zero
 decdrho_vdw(:) = zero
 decdgrho_vdw(:,:) = zero
 strsxc_vdw(:,:) = zero


 if ((xcdata%xclevel==0.or.ixc==0).and.(.not.my_add_tfw)) then
!  No xc at all is applied (usually for testing)
   MSG_WARNING('Note that no xc is applied (ixc=0).')

 else if (ixc/=20) then

!  Test: has a compensation density to be added/substracted (PAW) ?
   need_nhat=(nhatdim==1.and.usexcnhat==0)
   need_nhatgr=(nhatdim==1.and.nhatgrdim==1.and.ngrad==2.and.xcdata%intxc==0)
   test_nhat=(need_nhat.or.need_nhatgr)

!  The different components of depsxc will be
!  for nspden=1,   depsxc(:,1)=d(rho.exc)/d(rho) == (depsxcdrho) == (vxcrho)
!  and if ngrad=2, depsxc(:,2)=1/2*1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!  +1/|grad rho|*d(rho.exc)/d(|grad rho|)
!  == (1/2 * 1/|grho_up| * depsxcd|grho_up|) +  1/|grho| * depsxcd|grho|
!  (vxcgrho=1/|grho| * depsxcd|grho|)
!  (do not forget : |grad rho| /= |grad rho_up| + |grad rho_down|
!  and if use_laplacian, depsxc(:,3)=d(rho.exc)/d(lapl rho) == (depsxcdlrho) == (vxclrho)
!
!  for nspden>=2,  depsxc(:,1)=d(rho.exc)/d(rho_up) == (depsxcdrho_up) == (vxcrho_up)
!  depsxc(:,2)=d(rho.exc)/d(rho_down) == (depsxcdrho_dn) == (vxcrho_dn)
!  and if ngrad=2, depsxc(:,3)=1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|) == (1/|grho_up| * depsxcd|grho_up|) == (vxcgrho_up)
!  depsxc(:,4)=1/|grad rho_down|*d(rho.exc)/d(|grad rho_down|) == (1/|grho_dn| * depsxcd|grho_dn|) == (vxcgrho_dn)
!  depsxc(:,5)=1/|grad rho|*d(rho.exc)/d(|grad rho|) == (1/|grho| * depsxcd|grho|) == (vxcgrho)
!  and if use_laplacian, depsxc(:,6)=d(rho.exc)/d(lapl rho_up) == (depsxcdlrho_up) == (vxclrho_up)
!  depsxc(:,7)=d(rho.exc)/d(lapl rho_dn) == (depsxcdlrho_dn) == (vxclrho_dn)
!  Note: if nspden=4, rho_up=(rho+|m|)/2, rho_down=(rho-|m|)/2
   nspgrad=nspden_updn*ngrad;if(nspden_updn==2.and.ngrad==2)nspgrad=5
   if(uselaplacian==1) nspgrad=nspgrad+nspden_updn
   ABI_ALLOCATE(depsxc,(nfft,nspgrad))
   depsxc(:,:)=zero

!  PAW: select the valence density (and magnetization) to use:
!  link the correct density, according to usexcnhat option
   if ((.not.need_nhat).and.(.not.non_magnetic_xc)) then
     rhor_ => rhor
   else
     ABI_ALLOCATE(rhor_,(nfft,nspden))
     if (need_nhat) then
       do ispden=1,nspden
         do ifft=1,nfft
           rhor_(ifft,ispden)=rhor(ifft,ispden)-nhat(ifft,ispden)
         end do
       end do
     else
       do ispden=1,nspden
         do ifft=1,nfft
           rhor_(ifft,ispden)=rhor(ifft,ispden)
         end do
       end do
     end if
     if(non_magnetic_xc) then
       if(nspden==2) rhor_(:,2)=rhor_(:,1)*half
       if(nspden==4) rhor_(:,2:4)=zero
     endif
   end if
   if (usekden==1) then
     if(non_magnetic_xc) then
       ABI_ALLOCATE(taur_,(nfft,nspden))
       if(nspden==2) taur_(:,2)=taur_(:,1)*half
       if(nspden==4) taur_(:,2:4)=zero
     else
       taur_ => taur
     end if
   end if

!  Some initializations for the electron-positron correlation
   if (ipositron==2) then
     nspden_apn=1;ngrad_apn=1;iwarnp=1
     if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad_apn=2
     if (ngrad_apn==2.and.xcdata%xclevel<2) then
       message = 'GGA for the positron can only be performed with GGA pseudopotentials for the electron !'
       MSG_ERROR(message)
     end if
     if (ngrad_apn>1.and.option/=0.and.option/=1.and.option/=10.and.option/=12) then
       message = 'You cannot compute full GGA XC kernel for electrons-positron systems !'
       MSG_ERROR(message)
     end if
     ABI_ALLOCATE(depsxc_apn,(nfft,ngrad_apn))
   end if

!  Non-collinear magnetism: store norm of magnetization
!   m_norm_min= EPSILON(0.0_dp)**2 ! EB: TOO SMALL!!!
   m_norm_min=tol8  ! EB: tol14 is still too small, tests are underway
   if (nspden==4) then
     ABI_ALLOCATE(m_norm,(nfft))
     m_norm(:)=sqrt(rhor_(:,2)**2+rhor_(:,3)**2+rhor_(:,4)**2)
   end if

!  rhocorval will contain effective density used to compute gradients:
!  - with core density (if NLCC)
!  - without compensation density (if PAW under certain conditions)
!  - in (up+dn,up) or (n,mx,my,mz) format according to collinearity
!  of polarization and use of gradients (GGA)
   if (n3xccc>0.or.test_nhat.or.nspden_eff/=nspden) then
     ABI_ALLOCATE(rhocorval,(nfft,nspden_eff))
     if (nspden==nspden_eff) then
       rhocorval(:,1:nspden)=rhor_(:,1:nspden)
     else if (nspden==4) then
       rhocorval(:,1)=rhor_(:,1)
       rhocorval(:,2)=half*(rhor_(:,1)+m_norm(:))
     else
       rhocorval=zero
     end if
   else
     rhocorval => rhor_
   end if
   if (usekden==1.and.(n3xctau>0.or.nspden_eff/=nspden)) then
     ABI_ALLOCATE(taucorval,(nfft,nspden_eff))
     if (nspden==nspden_eff) then
       taucorval(:,1:nspden)=taur_(:,1:nspden)
     else
       taucorval=zero
     end if
   else
     taucorval => taur_
   end if

!  Add core electron density to effective density
   if (n3xccc>0) then
     rhocorval(:,1)=rhocorval(:,1)+xccc3d(:)
     if(nspden_eff==2) then
       rhocorval(:,2)=rhocorval(:,2)+half*xccc3d(:)
     end if
   end if
   if (n3xctau>0) then
     taucorval(:,1)=taucorval(:,1)+xcctau3d(:)
     if(nspden_eff==2) then
       taucorval(:,2)=taucorval(:,2)+half*xcctau3d(:)
     end if
   end if

!  If PAW, substract compensation density from effective density:
!  - if GGA, because nhat gradients are computed separately
   if (test_nhat.and.usexcnhat==1) then
     if (nspden==nspden_eff) then
       rhocorval(:,1:nspden)=rhocorval(:,1:nspden)-nhat(:,1:nspden)
     else if (nspden==4) then
       ABI_ALLOCATE(nhat_up,(nfft))
       do ifft=1,nfft
         if (m_norm(ifft)>m_norm_min) then
           nhat_up(ifft)=half*(nhat(ifft,1) &
&           +(rhor_(ifft,2)*nhat(ifft,2) &
&            +rhor_(ifft,3)*nhat(ifft,3) &
&            +rhor_(ifft,4)*nhat(ifft,4))/m_norm(ifft))
         else
           nhat_up(ifft)=half*(nhat(ifft,1) &
&           +sqrt(nhat(ifft,2)**2+nhat(ifft,3)**2+nhat(ifft,4)**2))
         end if
       end do
       rhocorval(:,1)=rhocorval(:,1)-nhat(:,1)
       rhocorval(:,2)=rhocorval(:,2)-nhat_up(:)
     end if
   end if

!  rhonow will contain effective density (and gradients if GGA)
!  taunow will contain effective kinetic energy density (if MetaGGA)
!  lrhonow will contain the laplacian if we have a MGGA
   ABI_ALLOCATE(rhonow,(nfft,nspden_eff,ngrad*ngrad))
   ABI_ALLOCATE(lrhonow,(nfft,nspden_eff*uselaplacian))
   ABI_ALLOCATE(taunow,(nfft,nspden_eff,usekden))

!  ====================================================================
!  ====================================================================
!  Loop on unshifted or shifted grids
   do ishift=0,xcdata%intxc

!    Set up density on unshifted or shifted grid (will be in rhonow(:,:,1)),
!    as well as the gradient of the density, also on the unshifted
!    or shifted grid (will be in rhonow(:,:,2:4)), if needed.
     if (uselaplacian==1) then
       call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,&
&                 qphon,rhocorval,rhonow,lrhonow=lrhonow)
     else
       call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,&
&                 qphon,rhocorval,rhonow)
     end if
     if (usekden==1) then
       call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,1,nspden_eff,&
&                 qphon,taucorval,taunow)
     end if

!    PAW+GGA: add "exact" gradients of compensation density
     !if (test_nhat.and.usexcnhat==1.and.ishift==0) then
     if (test_nhat.and.usexcnhat==1) then
       if (nspden==nspden_eff) then
         rhonow(:,1:nspden,1)=rhocorval(:,1:nspden)+nhat(:,1:nspden)
       else if (nspden==4) then
         rhonow(:,1,1)=rhocorval(:,1)+nhat(:,1)
         rhonow(:,2,1)=rhocorval(:,2)+nhat_up(:)
       end if
       if (ngrad==2.and.nhatgrdim==1.and.nspden==nspden_eff) then
         do ii=1,3
           jj=ii+1
           do ispden=1,nspden
             do ifft=1,nfft
               rhonow(ifft,ispden,jj)=rhonow(ifft,ispden,jj)+nhatgr(ifft,ispden,ii)
             end do
           end do
         end do
       end if
     end if

!    Deallocate temporary arrays
     if (ishift==xcdata%intxc) then
       if (n3xccc>0.or.test_nhat.or.nspden_eff/=nspden)  then
         ABI_DEALLOCATE(rhocorval)
       end if
       if (usekden==1.and.(n3xccc>0.or.nspden_eff/=nspden))  then
         ABI_DEALLOCATE(taucorval)
       end if
       if (test_nhat.and.nspden/=nspden_eff.and.usexcnhat==1)  then
         ABI_DEALLOCATE(nhat_up)
       end if
     end if

!    In case of non-collinear magnetism, extract up and down density and gradients (if GGA)
     if (nspden==4.and.nspden_eff==nspden) then
       if (ngrad==2) then
         do ifft=1,nfft
           gm_norm(1:3)=zero
           if(m_norm(ifft)>m_norm_min) then
!            if(m_norm(ifft)>rhonow(ifft,1,1)*tol10+tol14) then
             do jj=1,3  ! Compute here nabla(|m|)=(m.nabla(m))/|m| == (g|m| = m/|m| * gm)
               do ii=2,4
                 gm_norm(jj)=gm_norm(jj)+rhonow(ifft,ii,1+jj)*rhonow(ifft,ii,1)
               end do
             end do
             gm_norm(1:3)=gm_norm(1:3)/m_norm(ifft)
           end if
           rhonow(ifft,2,2)=half*(rhonow(ifft,1,2)+gm_norm(1))
           rhonow(ifft,2,3)=half*(rhonow(ifft,1,3)+gm_norm(2))
           rhonow(ifft,2,4)=half*(rhonow(ifft,1,4)+gm_norm(3))
         end do
       end if
       rhonow(:,2,1)=half*(rhonow(:,1,1)+m_norm(:))
       if (usekden==1) taunow(:,2,1)=half*(taunow(:,1,1)+m_norm(:))
     end if
!    Make the density positive everywhere (but do not care about gradients)
     call mkdenpos(iwarn,nfft,nspden_updn,1,rhonow(:,1:nspden_updn,1),xcdata%xc_denpos)

!    write(std_out,*) 'rhonow',rhonow

!    Uses a block formulation, in order to save simultaneously
!    CPU time and memory : xc routines
!    are called only once over mpts times, while the amount of allocated
!    space is kept at a low value, even if a lot of different
!    arrays are allocated, for use in different xc functionals.

     mpts=4000
     if (usekden==1) mpts=nfft   ! Why?

     do ifft=1,nfft,mpts
!      npts=mpts
!      npts is the number of points to be treated in this bunch
       npts=min(nfft-ifft+1,mpts)

!      Allocation of mandatory arguments of drivexc
       ABI_ALLOCATE(exc_b,(npts))
       ABI_ALLOCATE(rho_b,(npts))
       ABI_ALLOCATE(rho_b_updn,(npts,nspden_updn))
       ABI_ALLOCATE(vxcrho_b_updn,(npts,nspden_updn))
       vxcrho_b_updn(:,:)=zero

!      Allocation of optional arguments of drivexc
       ABI_ALLOCATE(grho2_b_updn,(npts,(2*nspden_updn-1)*usegradient))
       ABI_ALLOCATE(lrho_b_updn,(npts,nspden_updn*uselaplacian))
       ABI_ALLOCATE(tau_b_updn,(npts,nspden_updn*usekden))
       ABI_ALLOCATE(vxcgrho_b_updn,(npts,nvxcgrho))
       ABI_ALLOCATE(vxclrho_b_updn,(npts,nvxclrho))
       ABI_ALLOCATE(vxctau_b_updn,(npts,nvxctau))
       ABI_ALLOCATE(dvxc_b,(npts,ndvxc))
       ABI_ALLOCATE(d2vxc_b,(npts,nd2vxc))
       ABI_ALLOCATE(fxc_b,(npts*usefxc))
       if (nvxcgrho>0) vxcgrho_b_updn(:,:)=zero
       if (nvxclrho>0) vxclrho_b_updn(:,:)=zero
       if (nvxctau>0) vxctau_b_updn(:,:)=zero

       do ipts=ifft,ifft+npts-1
!        indx=ipts-ifft+1 varies from 1 to npts
         indx=ipts-ifft+1
         rho_b(indx)=rhonow(ipts,1,1)
         if(nspden_updn==1)then
           rho_b_updn(indx,1)=rhonow(ipts,1,1)*half
           if (usegradient==1) grho2_b_updn(indx,1)=quarter*(rhonow(ipts,1,2)**2 &
&                                       +rhonow(ipts,1,3)**2+rhonow(ipts,1,4)**2)
           if (usekden==1) tau_b_updn(indx,1)=taunow(ipts,1,1)*half
           if (uselaplacian==1) lrho_b_updn(indx,1)=lrhonow(ipts,1)*half
         else
           rho_b_updn(indx,1)=rhonow(ipts,2,1)
           rho_b_updn(indx,2)=rhonow(ipts,1,1)-rhonow(ipts,2,1)
           if(usegradient==1)then
             grho2_b_updn(indx,1)=rhonow(ipts,2,2)**2+   &
&             rhonow(ipts,2,3)**2+   &
&             rhonow(ipts,2,4)**2
             grho2_b_updn(indx,2)=(rhonow(ipts,1,2)-rhonow(ipts,2,2))**2 +   &
&             (rhonow(ipts,1,3)-rhonow(ipts,2,3))**2 +   &
&             (rhonow(ipts,1,4)-rhonow(ipts,2,4))**2
             grho2_b_updn(indx,3)=rhonow(ipts,1,2)**2+   &
&             rhonow(ipts,1,3)**2+   &
&             rhonow(ipts,1,4)**2
           end if
           if (usekden==1) then
             tau_b_updn(indx,1)=taunow(ipts,2,1)
             tau_b_updn(indx,2)=taunow(ipts,1,1)-taunow(ipts,2,1)
           end if
           if (uselaplacian==1) then
             lrho_b_updn(indx,1)=lrhonow(ipts,2)
             lrho_b_updn(indx,2)=lrhonow(ipts,1)-lrhonow(ipts,2)
           end if
         end if
       end do
!      In case of a hybrid functional, if one needs to compute the auxiliary GGA Kxc,
!      a separate call to drivexc is first needed to compute Kxc using such auxiliary GGA,
!      before calling again drivexc using the correct functional for Exc and Vxc.

       if(xcdata%usefock==1 .and. auxc_ixc/=0)then
         if (auxc_ixc<0) then
           call libxc_functionals_init(auxc_ixc,nspden,xc_functionals=xc_funcs_auxc)
         end if
         call drivexc(auxc_ixc,order,npts,nspden_updn,usegradient,0,0,&
&          rho_b_updn,exc_b,vxcrho_b_updn,nvxcgrho,0,0,ndvxc,nd2vxc, &
&          grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b_updn,dvxc=dvxc_b, &
&          fxcT=fxc_b,hyb_mixing=xcdata%hyb_mixing,el_temp=xcdata%tphysel,&
&          xc_funcs=xc_funcs_auxc)
!        Transfer the xc kernel
         if (nkxc_eff==1.and.ndvxc==15) then
           kxc(ifft:ifft+npts-1,1)=half*(dvxc_b(1:npts,1)+dvxc_b(1:npts,9)+dvxc_b(1:npts,10))
         else if (nkxc_eff==3.and.ndvxc==15) then
           kxc(ifft:ifft+npts-1,1)=dvxc_b(1:npts,1)+dvxc_b(1:npts,9)
           kxc(ifft:ifft+npts-1,2)=dvxc_b(1:npts,10)
           kxc(ifft:ifft+npts-1,3)=dvxc_b(1:npts,2)+dvxc_b(1:npts,11)
         end if
         if (auxc_ixc<0) then
           call libxc_functionals_end(xc_functionals=xc_funcs_auxc)
         end if
       end if
       if (present(xc_funcs)) then
         call libxc_functionals_get_hybridparams(hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&                                                hyb_range=hyb_range,xc_functionals=xc_funcs)
       else
         call libxc_functionals_get_hybridparams(hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&                                                hyb_range=hyb_range)
       end if

!      Call to main XC driver
       if (present(xc_funcs)) then
         call drivexc(ixc,order,npts,nspden_updn,&
&          usegradient,uselaplacian,usekden,&
&          rho_b_updn,exc_b,vxcrho_b_updn,&
&          nvxcgrho,nvxclrho,nvxctau,ndvxc,nd2vxc, &
&          grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b_updn,&
&          lrho_updn=lrho_b_updn,vxclrho=vxclrho_b_updn,&
&          tau_updn=tau_b_updn,vxctau=vxctau_b_updn,&
&          dvxc=dvxc_b,d2vxc=d2vxc_b,el_temp=xcdata%tphysel,fxcT=fxc_b,&
&          hyb_mixing=xcdata%hyb_mixing,&
&          xc_funcs=xc_funcs)
       else
         call drivexc(ixc,order,npts,nspden_updn,&
&          usegradient,uselaplacian,usekden,&
&          rho_b_updn,exc_b,vxcrho_b_updn,&
&          nvxcgrho,nvxclrho,nvxctau,ndvxc,nd2vxc, &
&          grho2_updn=grho2_b_updn,vxcgrho=vxcgrho_b_updn,&
&          lrho_updn=lrho_b_updn,vxclrho=vxclrho_b_updn,&
&          tau_updn=tau_b_updn,vxctau=vxctau_b_updn,&
&          dvxc=dvxc_b,d2vxc=d2vxc_b,el_temp=xcdata%tphysel,fxcT=fxc_b,&
&          hyb_mixing=xcdata%hyb_mixing)
       end if

!      If fake meta-GGA, has to remove the core contribution
!        when electronic effective mass has been modified
       if (n3xccc>0.and.(ixc==31.or.ixc==34.or.ixc==35)) then
         if (ixc==31.or.ixc==35) then
           coeff=one-(one/1.01_dp)
           if (nspden_updn==1) then
             coeff=half*coeff
             do ipts=1,npts
               exc_b(ipts)=exc_b(ipts)-coeff*xcctau3d(ifft+ipts-1) &
&                         /rho_b_updn(ipts,1)
             end do
           else
             do ipts=1,npts
               exc_b(ipts)=exc_b(ipts)-coeff*xcctau3d(ifft+ipts-1) &
&                         /(rho_b_updn(ipts,1)+rho_b_updn(ipts,2))
             end do
           end if
         else
           message = 'MetaGGA ixc=34 is not yet allowed with a core kinetic energy density!'
           MSG_ERROR(message)
         end if
       end if

!      Gradient Weiszacker correction to a Thomas-Fermi functional
       if (my_add_tfw) then
         vxcgrho_b_updn(:,:)=zero
         call xctfw(xcdata%tphysel,exc_b,fxc_b,usefxc,rho_b_updn,vxcrho_b_updn,npts,nspden_updn, &
&                   vxcgrho_b_updn,nvxcgrho,grho2_b_updn)
       end if

!      Accumulate enxc, strsxc and store vxc (and eventually kxc)
       dstrsxc=zero
       do ipts=ifft,ifft+npts-1
         indx=ipts-ifft+1
         epsxc=epsxc+rho_b(indx)*exc_b(indx)  !will be normalized with respect to the volume later to get enxc ("bigexc").
         depsxc(ipts,1)=vxcrho_b_updn(indx,1)
         exc_str=exc_b(indx);if(usefxc==1) exc_str=fxc_b(indx)
         if(nspden_updn==1)then
           strdiag=rho_b(indx)*(exc_str-vxcrho_b_updn(indx,1))
         else if(nspden_updn==2)then
           depsxc(ipts,2)=vxcrho_b_updn(indx,2)
!          Note : this is not the complete Vxc in the GGA case
           strdiag=rho_b(indx)*exc_str &
&           -rho_b_updn(indx,1)*vxcrho_b_updn(indx,1)&
&           -(rho_b(indx)-rho_b_updn(indx,1))*vxcrho_b_updn(indx,2)
         end if
         dstrsxc=dstrsxc+strdiag

!        For GGAs, additional terms appear
!        (the LB functional does not lead to additional terms)
         if(ngrad==2 .and. ixc/=13)then

!          Treat explicitely spin up, spin down and total spin for spin-polarized
!          Will exit when ispden=1 is finished if non-spin-polarized
           do ispden=1,3

             if(nspden_updn==1 .and. ispden>=2)exit

!            If the norm of the gradient vanishes, then the different terms vanishes,
!            but the inverse of the gradient diverges, so skip the update.
             if(grho2_b_updn(indx,ispden) < 1.0d-24) then
               depsxc(ipts,ispden+nspden_updn)=zero
               cycle
             end if

!            Compute the derivative of n.e_xc wrt the
!            spin up, spin down, or total density. In the non-spin-polarized
!            case take the coefficient that will be multiplied by the
!            gradient of the total density
             if(nspden_updn==1)then
!              !              Definition of vxcgrho_b_updn changed in v3.3
               if (nvxcgrho == 3) then
                 coeff=half*vxcgrho_b_updn(indx,1) + vxcgrho_b_updn(indx,3)
               else
                 coeff=half*vxcgrho_b_updn(indx,1)
               end if
             else if(nspden_updn==2)then
               if (nvxcgrho == 3) then
                 coeff=vxcgrho_b_updn(indx,ispden)
               else if (ispden /= 3) then
                 coeff=vxcgrho_b_updn(indx,ispden)
               else if (ispden == 3) then
                 coeff=zero
               end if
             end if
             depsxc(ipts,ispden+nspden_updn)=coeff

!            Store the gradient of up, down or total density, depending on ispden and nspden, at point ipts
             if(nspden_updn==1)then
               grho(1:3)=rhonow(ipts,1,2:4)
             else if(ispden==1 .and. nspden_updn==2)then
               grho(1:3)=rhonow(ipts,2,2:4)
             else if(ispden==2 .and. nspden_updn==2)then
               grho(1:3)=rhonow(ipts,1,2:4)-rhonow(ipts,2,2:4)
             else if(ispden==3 .and. nspden_updn==2)then
               grho(1:3)=rhonow(ipts,1,2:4)
             end if

!            In case of ixc 31 (mGGA functional fake 1),
!            skip the stress tensor to follow a LDA scheme (see doc/theory/MGGA/report_MGGA.pdf)
             if(ixc==31) cycle

!            Compute the contribution to the stress tensor
             s1=-grho(1)*grho(1)*coeff
             s2=-grho(2)*grho(2)*coeff
             s3=-grho(3)*grho(3)*coeff
!            The contribution of the next line comes from the part of Vxc
!            obtained from the derivative wrt the gradient
             dstrsxc=dstrsxc+s1+s2+s3
             strsxc1_tot=strsxc1_tot+s1
             strsxc2_tot=strsxc2_tot+s2
             strsxc3_tot=strsxc3_tot+s3
             strsxc4_tot=strsxc4_tot-grho(3)*grho(2)*coeff
             strsxc5_tot=strsxc5_tot-grho(3)*grho(1)*coeff
             strsxc6_tot=strsxc6_tot-grho(2)*grho(1)*coeff

           end do
         end if

!        For meta-GGAs, add the laplacian term (vxclrho) and/or kinetic energy density term (vxctau)
         if (usekden==1.and.with_vxctau) then
           if (nspden_updn==1)then
             vxctau(ipts,1,1) = vxctau_b_updn(indx,1)
           else if (nspden_updn==2)then
             vxctau(ipts,1,1) = vxctau_b_updn(indx,1)
             vxctau(ipts,2,1) = vxctau_b_updn(indx,2)
           end if
         end if
         if (uselaplacian==1) then
           if (nspden_updn==1)then
             depsxc(ipts,3) = vxclrho_b_updn(indx,1)
           else if (nspden_updn==2)then
             depsxc(ipts,6)   = vxclrho_b_updn(indx,1)
             depsxc(ipts,7)   = vxclrho_b_updn(indx,2)
           end if
         end if

       end do

!      Additional electron-positron correlation terms
       if (ipositron==2) then
!        Compute electron-positron XC energy per unit volume, potentials and derivatives
         ngr=0;if (ngrad_apn==2) ngr=npts
         ABI_ALLOCATE(fxc_apn,(npts))
         ABI_ALLOCATE(vxc_b_apn,(npts))
         ABI_ALLOCATE(vxcgr_apn,(ngr))
         ABI_ALLOCATE(vxc_ep,(npts))
         ABI_ALLOCATE(rhonow_apn,(npts,nspden_apn,1))
         ABI_ALLOCATE(grho2_apn,(ngr))
         rhonow_apn(1:npts,1,1)=electronpositron%rhor_ep(ifft:ifft+npts-1,1)
         if (usexcnhat==0) rhonow_apn(1:npts,1,1)=rhonow_apn(1:npts,1,1)-electronpositron%nhat_ep(ifft:ifft+npts-1,1)
         if (.not.electronpositron%posdensity0_limit) then
           call mkdenpos(iwarnp,npts,nspden_apn,1,rhonow_apn(:,1,1),xcdata%xc_denpos)
         end if
         if (ngrad_apn==2.and.usegradient==1) then
           if (nspden_apn==1) grho2_apn(:)=four*grho2_b_updn(:,1)
           if (nspden_apn==2) grho2_apn(:)=grho2_b_updn(:,3)
         end if
         if (ndvxc==0) then
           call xcpositron(fxc_apn,grho2_apn,electronpositron%ixcpositron,ngr,npts,&
&           electronpositron%posdensity0_limit,rho_b,&
&           rhonow_apn(:,1,1),vxc_b_apn,vxcgr_apn,vxc_ep)
         else
           ABI_ALLOCATE(dvxc_apn,(npts))
           call xcpositron(fxc_apn,grho2_apn,electronpositron%ixcpositron,ngr,npts,&
&           electronpositron%posdensity0_limit,rho_b,&
&           rhonow_apn(:,1,1),vxc_b_apn,vxcgr_apn,vxc_ep,dvxce=dvxc_apn)
         end if
         ABI_DEALLOCATE(vxc_ep)
         ABI_DEALLOCATE(rhonow_apn)
         ABI_DEALLOCATE(grho2_apn)
!        Accumulate electron-positron XC energies
         s1=zero
         do ipts=1,npts
           s1=s1+fxc_apn(ipts)
         end do
         electronpositron%e_xc=electronpositron%e_xc+s1*ucvol/dble(nfftot)
!        Add electron-positron dVxc_el/dRho_el to electron-electron one
         if (ndvxc==1) dvxc_b(:,1)=dvxc_b(:,1)+dvxc_apn(:)
         if (ndvxc==3) then
           dvxc_b(:,1)=dvxc_b(:,1)+four*dvxc_apn(:)
           dvxc_b(:,2)=dvxc_b(:,2)+four*dvxc_apn(:)
           dvxc_b(:,3)=dvxc_b(:,3)+four*dvxc_apn(:)
         end if
         if (ndvxc==15) then
           dvxc_b(:, 9)=dvxc_b(:, 9)+four*dvxc_apn(:)
           dvxc_b(:,10)=dvxc_b(:,10)+four*dvxc_apn(:)
           dvxc_b(:,11)=dvxc_b(:,11)+four*dvxc_apn(:)
         end if
!        Modify stresses - Compute factors for GGA
         do ipts=ifft,ifft+npts-1
           indx=ipts-ifft+1
           depsxc_apn(ipts,1)=vxc_b_apn(indx)
           dstrsxc=dstrsxc+fxc_apn(indx)-rho_b(indx)*vxc_b_apn(indx)
           if (ngrad_apn==2) then
             depsxc_apn(ipts,2)=vxcgr_apn(indx)
             s1=-grho(1)*grho(1)*vxcgr_apn(indx)
             s2=-grho(2)*grho(2)*vxcgr_apn(indx)
             s3=-grho(3)*grho(3)*vxcgr_apn(indx)
             dstrsxc=dstrsxc+s1+s2+s3
             strsxc1_tot=strsxc1_tot+s1
             strsxc2_tot=strsxc2_tot+s2
             strsxc3_tot=strsxc3_tot+s3
             strsxc4_tot=strsxc4_tot-grho(3)*grho(2)*vxcgr_apn(indx)
             strsxc5_tot=strsxc5_tot-grho(3)*grho(1)*vxcgr_apn(indx)
             strsxc6_tot=strsxc6_tot-grho(2)*grho(1)*vxcgr_apn(indx)
           end if ! GGA
         end do ! ipts
!        Deallocations
         ABI_DEALLOCATE(fxc_apn)
         ABI_DEALLOCATE(vxc_b_apn)
         ABI_DEALLOCATE(vxcgr_apn)
         if (ndvxc>0) then
           ABI_DEALLOCATE(dvxc_apn)
         end if
       end if

!      Transfer the xc kernel (if this must be done, and has not yet been done)
       if (nkxc_eff>0.and.ndvxc>0 .and. (xcdata%usefock==0 .or. auxc_ixc==0)) then
         if (nkxc_eff==1.and.ndvxc==15) then
           kxc(ifft:ifft+npts-1,1)=half*(dvxc_b(1:npts,1)+dvxc_b(1:npts,9)+dvxc_b(1:npts,10))
         else if (nkxc_eff==3.and.ndvxc==15) then
           kxc(ifft:ifft+npts-1,1)=dvxc_b(1:npts,1)+dvxc_b(1:npts,9)
           kxc(ifft:ifft+npts-1,2)=dvxc_b(1:npts,10)
           kxc(ifft:ifft+npts-1,3)=dvxc_b(1:npts,2)+dvxc_b(1:npts,11)
         else if (nkxc_eff==7.and.ndvxc==8) then
           kxc(ifft:ifft+npts-1,1)=half*dvxc_b(1:npts,1)
           kxc(ifft:ifft+npts-1,2)=half*dvxc_b(1:npts,3)
           kxc(ifft:ifft+npts-1,3)=quarter*dvxc_b(1:npts,5)
           kxc(ifft:ifft+npts-1,4)=eighth*dvxc_b(1:npts,7)
         else if (nkxc_eff==7.and.ndvxc==15) then
           kxc(ifft:ifft+npts-1,1)=half*(dvxc_b(1:npts,1)+dvxc_b(1:npts,9)+dvxc_b(1:npts,10))
           kxc(ifft:ifft+npts-1,2)=half*dvxc_b(1:npts,3)+dvxc_b(1:npts,12)
           kxc(ifft:ifft+npts-1,3)=quarter*dvxc_b(1:npts,5)+dvxc_b(1:npts,13)
           kxc(ifft:ifft+npts-1,4)=eighth*dvxc_b(1:npts,7)+dvxc_b(1:npts,15)
         else if (nkxc_eff==19.and.ndvxc==15) then
           kxc(ifft:ifft+npts-1,1)=dvxc_b(1:npts,1)+dvxc_b(1:npts,9)
           kxc(ifft:ifft+npts-1,2)=dvxc_b(1:npts,10)
           kxc(ifft:ifft+npts-1,3)=dvxc_b(1:npts,2)+dvxc_b(1:npts,11)
           kxc(ifft:ifft+npts-1,4)=dvxc_b(1:npts,3)
           kxc(ifft:ifft+npts-1,5)=dvxc_b(1:npts,4)
           kxc(ifft:ifft+npts-1,6)=dvxc_b(1:npts,5)
           kxc(ifft:ifft+npts-1,7)=dvxc_b(1:npts,6)
           kxc(ifft:ifft+npts-1,8)=dvxc_b(1:npts,7)
           kxc(ifft:ifft+npts-1,9)=dvxc_b(1:npts,8)
           kxc(ifft:ifft+npts-1,10)=dvxc_b(1:npts,12)
           kxc(ifft:ifft+npts-1,11)=dvxc_b(1:npts,13)
           kxc(ifft:ifft+npts-1,12)=dvxc_b(1:npts,14)
           kxc(ifft:ifft+npts-1,13)=dvxc_b(1:npts,15)
         else ! All other cases
           kxc(ifft:ifft+npts-1,1:nkxc_eff)=zero
           kxc(ifft:ifft+npts-1,1:min(nkxc_eff,ndvxc))=dvxc_b(1:npts,1:min(nkxc_eff,ndvxc))
         end if
         if (nkxc_eff==7) then
           kxc(ifft:ifft+npts-1,5)=rhonow(ifft:ifft+npts-1,1,2)
           kxc(ifft:ifft+npts-1,6)=rhonow(ifft:ifft+npts-1,1,3)
           kxc(ifft:ifft+npts-1,7)=rhonow(ifft:ifft+npts-1,1,4)
         else if (nkxc_eff==19) then
           kxc(ifft:ifft+npts-1,14)=rhonow(ifft:ifft+npts-1,1,2)
           kxc(ifft:ifft+npts-1,15)=rhonow(ifft:ifft+npts-1,2,2)
           kxc(ifft:ifft+npts-1,16)=rhonow(ifft:ifft+npts-1,1,3)
           kxc(ifft:ifft+npts-1,17)=rhonow(ifft:ifft+npts-1,2,3)
           kxc(ifft:ifft+npts-1,18)=rhonow(ifft:ifft+npts-1,1,4)
           kxc(ifft:ifft+npts-1,19)=rhonow(ifft:ifft+npts-1,2,4)
         end if
       end if

!      Transfer the XC 3rd-derivative
       if (abs(option)==3.and.order==3.and.nd2vxc>0) then
         k3xc(ifft:ifft+npts-1,1:nd2vxc)=d2vxc_b(1:npts,1:nd2vxc)
       end if

!      Add the diagonal part to the xc stress
       strsxc1_tot=strsxc1_tot+dstrsxc
       strsxc2_tot=strsxc2_tot+dstrsxc
       strsxc3_tot=strsxc3_tot+dstrsxc

       ABI_DEALLOCATE(exc_b)
       ABI_DEALLOCATE(rho_b)
       ABI_DEALLOCATE(rho_b_updn)
       ABI_DEALLOCATE(grho2_b_updn)
       ABI_DEALLOCATE(vxcrho_b_updn)
       ABI_DEALLOCATE(dvxc_b)
       ABI_DEALLOCATE(d2vxc_b)
       ABI_DEALLOCATE(vxcgrho_b_updn)
       ABI_DEALLOCATE(fxc_b)
       ABI_DEALLOCATE(vxclrho_b_updn)
       ABI_DEALLOCATE(lrho_b_updn)
       ABI_DEALLOCATE(tau_b_updn)
       ABI_DEALLOCATE(vxctau_b_updn)

!      End of the loop on blocks of data
     end do

     strsxc(1)=strsxc1_tot
     strsxc(2)=strsxc2_tot
     strsxc(3)=strsxc3_tot
     strsxc(4)=strsxc4_tot
     strsxc(5)=strsxc5_tot
     strsxc(6)=strsxc6_tot

!    If GGA, multiply the gradient of the density by the proper
!    local partial derivatives of the XC functional
     rhonow_ptr => rhonow
     if (ipositron==2) then
       ABI_ALLOCATE(rhonow_ptr,(nfft,nspden_eff,ngrad*ngrad))
       rhonow_ptr=rhonow
     end if
     if(ngrad==2 .and. ixc/=13)then
       call xcmult(depsxc,nfft,ngrad,nspden_eff,nspgrad,rhonow_ptr)
     end if

!    Compute contribution from this grid to vxc, and ADD to existing vxc
     if (nspden/=4) then
       if(with_vxctau)then
         call xcpot(cplex,gprimd,ishift,uselaplacian,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&         qphon,depsxc=depsxc,rhonow=rhonow_ptr,vxc=vxc,vxctau=vxctau)
       else
         call xcpot(cplex,gprimd,ishift,uselaplacian,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&         qphon,depsxc=depsxc,rhonow=rhonow_ptr,vxc=vxc)
       end if

     else

!      If non-collinear magnetism, restore potential in proper axis before adding it
       ABI_ALLOCATE(vxcrho_b_updn,(nfft,4))
       vxcrho_b_updn=zero
       call xcpot(cplex,gprimd,ishift,uselaplacian,mpi_enreg,nfft,ngfft,ngrad,nspden_eff,nspgrad,&
&       qphon,depsxc=depsxc,rhonow=rhonow_ptr,vxc=vxcrho_b_updn)
       do ifft=1,nfft
         dvdn=half*(vxcrho_b_updn(ifft,1)+vxcrho_b_updn(ifft,2))
         if(m_norm(ifft)>m_norm_min) then
!          if(m_norm(ifft)>rhor_(ifft,1)*tol10+tol14) then
           dvdz=half*(vxcrho_b_updn(ifft,1)-vxcrho_b_updn(ifft,2))/m_norm(ifft)
           vxc(ifft,1)=vxc(ifft,1)+dvdn+rhor_(ifft,4)*dvdz
           vxc(ifft,2)=vxc(ifft,2)+dvdn-rhor_(ifft,4)*dvdz
           vxc(ifft,3)=vxc(ifft,3)+rhor_(ifft,2)*dvdz
           vxc(ifft,4)=vxc(ifft,4)-rhor_(ifft,3)*dvdz
         else
           vxc(ifft,1:2)=vxc(ifft,1:2)+dvdn
         end if
       end do
       ABI_DEALLOCATE(vxcrho_b_updn)
     end if
     if (ipositron==2)  then
       ABI_DEALLOCATE(rhonow_ptr)
     end if
     nullify(rhonow_ptr)

!    Add electron-positron XC potential to electron-electron one
!    Eventually compute GGA contribution
     if (ipositron==2) then
       ABI_ALLOCATE(rhonow_apn,(nfft,nspden_apn,ngrad_apn**2))
       rhonow_apn(1:nfft,1,1:ngrad_apn**2)=rhonow(1:nfft,1,1:ngrad_apn**2)
       if (ngrad_apn==2) then
         call xcmult(depsxc_apn,nfft,ngrad_apn,nspden_apn,ngrad_apn,rhonow_apn)
       end if
       ABI_ALLOCATE(vxc_apn,(nfft,nspden_apn))
       vxc_apn=zero
       call xcpot(cplex,gprimd,ishift,0,mpi_enreg,nfft,ngfft,ngrad_apn,&
&       nspden_apn,ngrad_apn,qphon,depsxc=depsxc_apn,rhonow=rhonow_apn,vxc=vxc_apn)
       vxc(:,1)=vxc(:,1)+vxc_apn(:,1)
       if (nspden_updn==2) vxc(:,2)=vxc(:,2)+vxc_apn(:,1)
       s1=zero
       do ipts=1,nfft
         s1=s1+vxc_apn(ipts,1)*rhonow(ipts,1,1)
       end do
       electronpositron%e_xcdc=electronpositron%e_xcdc+s1*ucvol/dble(nfftot)
       ABI_DEALLOCATE(rhonow_apn)
       ABI_DEALLOCATE(vxc_apn)
       ABI_DEALLOCATE(depsxc_apn)
     end if

!    End loop on unshifted or shifted grids
   end do

!  Calculate van der Waals correction when requested
#if defined DEV_YP_VDWXC
   if ( (xcdata%vdw_xc > 0) .and. (xcdata%vdw_xc < 3) .and. (xc_vdw_status()) ) then
     call xc_vdw_aggregate(ucvol,gprimd,nfft,nspden_updn,ngrad*ngrad, &
&     ngfft(1),ngfft(2),ngfft(3),rhonow, &
&     deltae_vdw,exc_vdw,decdrho_vdw,decdgrho_vdw,strsxc_vdw)
   end if
#else
   if ( (xcdata%vdw_xc > 0) .and. (xcdata%vdw_xc < 3) ) then
     write(message,'(3a)')&
&     'vdW-DF functionals are not fully operational yet.',ch10,&
&     'Action : modify vdw_xc'
     MSG_ERROR(message)
   end if
#endif
!  Normalize enxc, strsxc and vxc
   divshft=one/dble(xcdata%intxc+1)
   strsxc(:)=strsxc(:)/dble(nfftot)*divshft
   enxc=epsxc*ucvol/dble(nfftot)*divshft
   vxc=vxc*divshft
   if (with_vxctau) vxctau=vxctau*divshft

!  Reduction in case of FFT distribution
   if (nproc_fft>1)then
     call timab(48,1,tsec)
     call xmpi_sum(strsxc,comm_fft ,ierr)
     call xmpi_sum(enxc  ,comm_fft ,ierr)
     if (ipositron==2) then
       s1=electronpositron%e_xc;s2=electronpositron%e_xcdc
       call xmpi_sum(s1,comm_fft ,ierr)
       call xmpi_sum(s2,comm_fft ,ierr)
       electronpositron%e_xc=s1;electronpositron%e_xcdc=s2
     end if
     call timab(48,2,tsec)
   end if

!  Compute vxcavg
   call mean_fftr(vxc,vxcmean,nfft,nfftot,min(nspden,2),mpi_comm_sphgrid=comm_fft)
   if(nspden==1)then
     vxcavg=vxcmean(1)
   else
     vxcavg=half*(vxcmean(1)+vxcmean(2))
   end if

   ABI_DEALLOCATE(depsxc)
   ABI_DEALLOCATE(rhonow)
   ABI_DEALLOCATE(lrhonow)
   ABI_DEALLOCATE(taunow)
   if (need_nhat.or.non_magnetic_xc) then
     ABI_DEALLOCATE(rhor_)
   end if
   if ((usekden==1).and.(non_magnetic_xc)) then
     ABI_DEALLOCATE(taur_)
   end if
   if (allocated(m_norm))  then
     ABI_DEALLOCATE(m_norm)
   end if

 end if

!Treat separately the Fermi-Amaldi correction.
 if (ixc==20 .or. ixc==21 .or. ixc==22) then
   if(present(vhartr))then

!    Fermi-Amaldi correction : minus Hartree divided by the
!    number of electrons per unit cell. This is not size consistent, but
!    interesting for isolated systems with a few electrons.
!    nelect=ucvol*rhog(1,1)
     factor=-one/xcdata%nelect
     vxc(:,1)=factor*vhartr(:)
     if(nspden>=2) vxc(:,2)=factor*vhartr(:)

!    Compute corresponding xc energy and stress as well as vxcavg
     call dotprod_vn(1,rhor,enxc,doti,nfft,nfftot,1,1,vxc,ucvol,mpi_comm_sphgrid=comm_fft)
     enxc=half*enxc
     strsxc(1:3)=-enxc/ucvol

!    Compute average of vxc (one component only).
     call mean_fftr(vxc,vxcmean,nfft,nfftot,1,mpi_comm_sphgrid=comm_fft)
     vxcavg = vxcmean(1)
!    For ixc=20, the local exchange-correlation kernel is zero, but the Hartree
!    kernel will be modified in tddft. No other use of kxc should be made with ixc==20
     if(nkxc/=0 .and. ixc==20) kxc(:,:)=zero
!    For ixc=21 or 22, the LDA (ixc=1) kernel has been computed previously.

   else

     MSG_BUG('When ixc=20,21 or 22, vhartr needs to be present in the call to rhotoxc !')

   end if

 end if

!Add van der Waals terms
#if defined DEV_YP_VDWXC
 if ( (xcdata%vdw_xc > 0) .and. (xcdata%vdw_xc < 10) .and. (xc_vdw_status()) ) then
   enxc = enxc + exc_vdw
   do ispden=1,nspden
     vxc(:,ispden) = vxc(:,ispden) + decdrho_vdw(ispden)
   end do
   strsxc(1) = strsxc(1) + strsxc_vdw(1,1)
   strsxc(2) = strsxc(2) + strsxc_vdw(2,2)
   strsxc(3) = strsxc(3) + strsxc_vdw(3,3)
   strsxc(4) = strsxc(4) + strsxc_vdw(3,2)
   strsxc(5) = strsxc(5) + strsxc_vdw(3,1)
   strsxc(6) = strsxc(6) + strsxc_vdw(2,1)
 end if
#endif
 if ( present(exc_vdw_out) ) exc_vdw_out = exc_vdw

 call timab(81,2,tsec)

 DBG_EXIT("COLL")

end subroutine rhotoxc
!!***

end module m_rhotoxc
!!***

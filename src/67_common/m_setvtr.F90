!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_setvtr
!! NAME
!!  m_setvtr
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (XG, GMR, FJ, MT, EB, SPr)
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

module m_setvtr

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_errors
 use m_abi2big
 use m_xmpi
 use m_xcdata
 use m_dtset

 use defs_datatypes,      only : pseudopotential_type
 use defs_abitypes,       only : MPI_type
 use m_time,              only : timab
 use m_geometry,          only : xred2xcart
 use m_cgtools,           only : dotprod_vn
 use m_ewald,             only : ewald
 use m_energies,          only : energies_type
 use m_electronpositron,  only : electronpositron_type, electronpositron_calctype, rhohxcpositron
 use libxc_functionals,   only : libxc_functionals_is_hybrid
 use m_pawrad,            only : pawrad_type
 use m_pawtab,            only : pawtab_type
 use m_jellium,           only : jellium
 use m_spacepar,          only : hartre
 use m_dens,              only : constrained_dft_t,constrained_dft_ini,constrained_dft_free,mag_penalty
 use m_vdw_dftd2,         only : vdw_dftd2
 use m_vdw_dftd3,         only : vdw_dftd3
 use m_atm2fft,           only : atm2fft
 use m_rhotoxc,           only : rhotoxc
 use m_mklocl,            only : mklocl
 use m_xchybrid,          only : xchybrid_ncpp_cc
 use m_mkcore,            only : mkcore, mkcore_alt
 use m_psolver,           only : psolver_rhohxc
 use m_wvl_psi,          only : wvl_psitohpsi
 use m_mkcore_wvl,       only : mkcore_wvl

#if defined HAVE_BIGDFT
 use BigDFT_API, only: denspot_set_history
#endif

 implicit none

 private
!!***

 public :: setvtr
!!***

contains
!!***

!!****f* m_setvtr/setvtr
!! NAME
!! setvtr
!!
!! FUNCTION
!! Set up the trial potential and some energy terms
!!
!! INPUTS
!!  [add_tfw]=flag controling the addition of Weiszacker gradient correction to Thomas-Fermi kin energy
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   |       if =0,1 no xc kernel, =2 spin-averaged (LDA) kernel
!!   | densfor_pred=govern the choice of preconditioner for the SCF cycle
!!   | iscf=determines the way the SCF cycle is handled
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | qprtrb(3)= integer wavevector of possible perturbing potential
!!   |            in basis of reciprocal lattice translations
!!   | typat(natom)=type integer for each atom in cell
!!   | vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!   |  perturbing potential is added of the form
!!   |  V(G)=(vprtrb(1)+I*vprtrb(2))/2 at the values G=qprtrb and
!!   |  (vprtrb(1)-I*vprtrb(2))/2 at G=-qprtrb (integers)
!!   |  for each type of atom, from psp (used in norm-conserving only)
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2) (sphere for density and potential)
!!  istep=step number in the main loop of scfcv
!!  mgfft=maximum size of 1D FFTs
!!  moved_rhor=1 if the density was moved just before
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc
!!  ntypat=number of types of atoms in unit cell.
!!  n1xccc=dimension of xccc1d; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=>0 if some additional energies have to be computed
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase (structure factor) information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of electron density
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!   | definition for spin components:
!!   | case of nspden = 2
!!   |      rhor(:,1) => rho_up + rho_dwn
!!   |      rhor(:,2) => rho_up
!!   | case of nspden = 4
!!   |      rhor(:,1)   => rho_upup + rho_dwndwn
!!   |      rhor(:,2:4) => {m_x,m_y,m_z}
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  [taur(nfftf,nspden*dtset%usekden)]=array for kinetic energy density
!!  ucvol = unit cell volume (bohr^3)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_xc=exchange-correlation energy (hartree)
!!   | In case of hybrid compensation algorithm:
!!   | e_hybcomp_E0=energy compensation term for hybrid exchange-correlation energy (hartree) at fixed density
!!   | e_hybcomp_v0=potential compensation term for hybrid exchange-correlation energy (hartree) at fixed density
!!   | e_hybcomp_v=potential compensation term for hybrid exchange-correlation energy (hartree) at self-consistent density
!!  ==== if optene==2 or 4
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if dtset%icoulomb == 0
!!   | e_ewald=Ewald energy (hartree)
!!  ==== if optene>=1
!!   | e_hartree=Hartree part of total energy (hartree)
!!  ==== if optene==3 or 4
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!  ==== if dtset%vdw_xc == 5 or 6 or 7
!!   | e_vdw_dftd=Dispersion energy from DFT-D Van der Waals correction (hartree)
!!  grchempottn(3,natom)=grads of spatially-varying chemical energy (hartree)
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D2 dispersion (hartree)
!!  kxc(nfft,nkxc)=exchange-correlation kernel, will be computed if nkxc/=0 .
!!                 see routine rhotoxc for a more complete description
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxcavg=mean of the vxc potential
!!
!! SIDE EFFECTS
!!  [electronpositron <type(electronpositron_type)>]=quantities for the electron-positron annihilation (optional argument)
!!  moved_atm_inside=1 if the atomic positions were moved inside the SCF loop.
!!  vhartr(nfft)=Hartree potential (Hartree)
!!  vpsp(nfft)=local psp (Hartree)
!!  vtrial(nfft,nspden)= trial potential (Hartree)
!!  vxc(nfft,nspden)= xc potential (Hartree)
!!  [vxc_hybcomp(nfft,nspden)= compensation xc potential (Hartree) in case of hybrids] Optional output
!!       i.e. difference between the hybrid Vxc at fixed density and the auxiliary Vxc at fixed density
!!  [vxctau(nfftf,dtset%nspden,4*dtset%usekden)]=derivative of XC energy density with respect to
!!    kinetic energy density (metaGGA cases) (optional output)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  [xccctau3d(n3xccc)]=3D core electron kinetic energy density for XC core correction, bohr^-3
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  Developers have to be careful when introducing others arrays: they have to be stored on the fine FFT grid.

!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      bethe_salpeter,scfcv,screening,sigma
!!
!! CHILDREN
!!      atm2fft,denspot_set_history,dotprod_vn,ewald,hartre,ionion_realspace
!!      ionion_surface,jellium,mag_penalty,mkcore,mkcore_alt,mkcore_wvl,mklocl
!!      psolver_rhohxc,rhohxcpositron,rhotoxc,spatialchempot,timab,vdw_dftd2
!!      vdw_dftd3,wvl_psitohpsi,wvl_vtrial_abi2big,xcdata_init,xchybrid_ncpp_cc
!!      xred2xcart
!!
!! SOURCE

subroutine setvtr(atindx1,dtset,energies,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcut,&
&  istep,kxc,mgfft,moved_atm_inside,moved_rhor,mpi_enreg,&
&  nattyp,nfft,ngfft,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,ntypat,n1xccc,n3xccc,&
&  optene,pawrad,pawtab,ph1d,psps,rhog,rhor,rmet,rprimd,strsxc,&
&  ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,xccc3d,xred,&
&  electronpositron,taur,vxc_hybcomp,vxctau,add_tfw,xcctau3d) ! optionals arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mgfft,n1xccc,n3xccc,nfft,ngrvdw,nhatgrdim,nkxc,ntypat
 integer,intent(in) :: optene,usexcnhat
 integer,intent(inout) :: moved_atm_inside,moved_rhor
 logical,intent(in),optional :: add_tfw
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vxcavg
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data), intent(inout) :: wvl
!arrays
 integer, intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(:,:,:) !(nfft,dtset%nspden,3*nhatgrdim)
 real(dp),intent(in) :: rhog(2,nfft)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout),optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
 real(dp),intent(out),optional :: vxc_hybcomp(:,:) ! (nfft,nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(inout),optional ::xcctau3d(n3xccc*dtset%usekden)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(out) :: grchempottn(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom),grvdw(3,ngrvdw),kxc(nfft,nkxc),strsxc(6)
 type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: coredens_method,coretau_method,mpi_comm_sphgrid,nk3xc
 integer :: iatom,ifft,ipositron,ispden,nfftot
 integer :: optatm,optdyfr,opteltfr,optgr,option,option_eff,optn,optn2,optstr,optv,vloc_method
 real(dp) :: doti,e_xcdc_vxctau,ebb,ebn,evxc,ucvol_local,rpnrm
 logical :: add_tfw_,is_hybrid_ncpp,non_magnetic_xc,with_vxctau,wvlbigdft
 real(dp), allocatable :: xcart(:,:)
 character(len=500) :: message
 type(constrained_dft_t) :: constrained_dft
 type(xcdata_type) :: xcdata,xcdatahyb
!arrays
 real(dp),parameter :: identity(1:4)=(/1._dp,1._dp,0._dp,0._dp/)
 real(dp) :: dummy6(6),tsec(2)
 real(dp) :: grewtn_fake(3,1)
 real(dp) :: dummy_in(0)
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0)
 real(dp) :: strn_dummy6(6), strv_dummy6(6)
 real(dp) :: vzeeman(4)
 real(dp),allocatable :: grtn(:,:),dyfr_dum(:,:,:),gr_dum(:,:)
 real(dp),allocatable :: rhojellg(:,:),rhojellr(:),rhowk(:,:),vjell(:)
 real(dp),allocatable :: v_constr_dft_r(:,:),rhog_dum(:,:)

! *********************************************************************

 call timab(91,1,tsec)

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = (present(vxctau).and.present(taur).and.(dtset%usekden/=0))

!Check if we're in hybrid norm conserving pseudopotential with a core correction
 is_hybrid_ncpp=(dtset%usepaw==0 .and. n3xccc/=0 .and. &
& (dtset%ixc==41.or.dtset%ixc==42.or.libxc_functionals_is_hybrid()))

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)

!Get size of FFT grid
 nfftot=PRODUCT(ngfft(1:3))

!mpi
 mpi_comm_sphgrid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) mpi_comm_sphgrid=mpi_enreg%comm_wvl

!Test electron-positron case
 ipositron=0;if (present(electronpositron)) ipositron=electronpositron_calctype(electronpositron)

!Test addition of Weiszacker gradient correction to Thomas-Fermi kin energy
 add_tfw_=.false.;if (present(add_tfw)) add_tfw_=add_tfw

!Get Ewald energy and Ewald forces, as well as vdW-DFTD energy and forces, and chemical potential energy and forces.
!-------------------------------------------------------------------------------------------------------------------
 call timab(5,1,tsec)
 if (ipositron/=1) then
   if (dtset%icoulomb == 0 .or. (dtset%usewvl == 0 .and. dtset%icoulomb == 2)) then
!    Periodic system, need to compute energy and forces due to replica and
!    to correct the shift in potential calculation.
     call ewald(energies%e_ewald,gmet,grewtn,dtset%natom,ntypat,rmet,dtset%typat,ucvol,xred,psps%ziontypat)
!    For a periodic system bearing a finite charge, the monopole correction to the
!    energy is relevant.
!    See Leslie and Gillan, JOURNAL OF PHYSICS C-SOLID STATE PHYSICS 18, 973 (1985)
     if(abs(dtset%charge)>tol8) then
       call ewald(energies%e_monopole,gmet,grewtn_fake,1,1,rmet,(/1/),ucvol,(/0.0_dp,0.0_dp,0.0_dp/),(/dtset%charge/))
       energies%e_monopole=-energies%e_monopole
     end if
   else if (dtset%icoulomb == 1) then
!    In a non periodic system (real space computation), the G=0 divergence
!    doesn't occur and ewald is not needed. Only the ion/ion interaction
!    energy is relevant and used as Ewald energy and gradient.
     call ionion_realSpace(dtset, energies%e_ewald, grewtn, rprimd, xred, psps%ziontypat)
   else if (dtset%icoulomb == 2) then
     call ionion_surface(dtset, energies%e_ewald, grewtn, mpi_enreg%me_wvl, mpi_enreg%nproc_wvl, rprimd, &
&     wvl%descr, wvl%den, xred)
   end if
   if (dtset%nzchempot>0) then
     call spatialchempot(energies%e_chempot,dtset%chempot,grchempottn,dtset%natom,ntypat,dtset%nzchempot,dtset%typat,xred)
   end if
   if (dtset%vdw_xc==5.and.ngrvdw==dtset%natom) then
     call vdw_dftd2(energies%e_vdw_dftd,dtset%ixc,dtset%natom,ntypat,1,dtset%typat,rprimd,&
&     dtset%vdw_tol,xred,psps%znucltypat,fred_vdw_dftd2=grvdw)
   end if
   if ((dtset%vdw_xc==6.or.dtset%vdw_xc==7).and.ngrvdw==dtset%natom) then
     call vdw_dftd3(energies%e_vdw_dftd,dtset%ixc,dtset%natom,&
&     ntypat,1,dtset%typat,rprimd,dtset%vdw_xc,dtset%vdw_tol,dtset%vdw_tol_3bt,&
&     xred,psps%znucltypat,fred_vdw_dftd3=grvdw)
   end if
 else
   energies%e_ewald=zero
   energies%e_chempot=zero
   grchempottn=zero
   grewtn=zero
   energies%e_vdw_dftd=zero
   if (ngrvdw>0) grvdw=zero
 end if
 call timab(5,2,tsec)

!Compute parts of total energy depending on potentials
!--------------------------------------------------------------
 if (dtset%usewvl == 0) then
   ucvol_local = ucvol
#if defined HAVE_BIGDFT
 else
!  We need to tune the volume when wavelets are used because, not all FFT points are used.
!  ucvol_local = (half * dtset%wvl_hgrid) ** 3 * ngfft(1)*ngfft(2)*ngfft(3)
   ucvol_local = product(wvl%den%denspot%dpbox%hgrids) * real(product(wvl%den%denspot%dpbox%ndims), dp)
#endif
 end if

!Determine by which method the local ionic potential and/or the pseudo core charge density
! have to be computed
!Local ionic potential:
! Method 1: PAW
! Method 2: Norm-conserving PP, icoulomb>0, wavelets
 vloc_method=1;if (psps%usepaw==0) vloc_method=2
 if (dtset%icoulomb>0) vloc_method=2
 if (psps%usewvl==1) vloc_method=2
!Pseudo core charge density:
! Method 1: PAW, nc_xccc_gspace
! Method 2: Norm-conserving PP, wavelets
 coredens_method=1;if (psps%usepaw==0) coredens_method=2
 if (psps%nc_xccc_gspace==1) coredens_method=1
 if (psps%nc_xccc_gspace==0) coredens_method=2
 if (psps%usewvl==1) coredens_method=2
 coretau_method=0
 if (dtset%usekden==1.and.psps%usepaw==1) then
   coretau_method=1;if (psps%nc_xccc_gspace==0) coretau_method=2
 end if
!In some specific cases, XC has to be handled as non-magnetic
 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

!Local ionic potential and/or pseudo core charge by method 1
 if (vloc_method==1.or.coredens_method==1) then
   call timab(552,1,tsec)
   optv=0;if (vloc_method==1) optv=1
   optn=0;if (coredens_method==1) optn=n3xccc/nfft
   optatm=1;optdyfr=0;opteltfr=0;optgr=0;optstr=0;optn2=1
   call atm2fft(atindx1,xccc3d,vpsp,dummy_out1,dummy_out2,dummy_out3,dummy_in,&
&   gmet,gprimd,dummy_out4,dummy_out5,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,&
&   dummy_in,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,dummy_in,dummy_in,dummy_in,dtset%vprtrb,psps%vlspl,&
&   comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&   paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
   call timab(552,2,tsec)
 end if
 if (coretau_method==1) then
   call timab(552,1,tsec)
   optv=0;optn=1
   optatm=1;optdyfr=0;opteltfr=0;optgr=0;optstr=0;optn2=4
   call atm2fft(atindx1,xcctau3d,dummy_out6,dummy_out1,dummy_out2,dummy_out3,dummy_in,&
&   gmet,gprimd,dummy_out4,dummy_out5,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,&
&   dummy_in,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,dummy_in,dummy_in,dummy_in,dtset%vprtrb,psps%vlspl,&
&   comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&   paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
   call timab(552,2,tsec)
 end if

!Local ionic potential by method 2
 if (vloc_method==2) then
   option=1
   ABI_ALLOCATE(gr_dum,(3,dtset%natom))
   ABI_ALLOCATE(dyfr_dum,(3,3,dtset%natom))
   ABI_ALLOCATE(rhog_dum,(2,nfft))
   call mklocl(dtset,dyfr_dum,energies%e_localpsp,gmet,gprimd,&
&   gr_dum,gsqcut,dummy6,mgfft,mpi_enreg,dtset%natom,nattyp,&
&   nfft,ngfft,dtset%nspden,ntypat,option,pawtab,ph1d,psps,&
&   dtset%qprtrb,rhog_dum,rhor,rprimd,ucvol,dtset%vprtrb,vpsp,wvl%descr,wvl%den,xred)
   ABI_DEALLOCATE(gr_dum)
   ABI_DEALLOCATE(dyfr_dum)
   ABI_DEALLOCATE(rhog_dum)
 end if

!3D pseudo core electron density xccc3d by method 2
 if (coredens_method==2.and.n1xccc/=0) then
   call timab(91,2,tsec)
   call timab(92,1,tsec)
   option=1
   ABI_ALLOCATE(gr_dum,(3,dtset%natom))
   ABI_ALLOCATE(dyfr_dum,(3,3,dtset%natom))
   if (psps%usewvl==0.and.psps%usepaw==0.and.dtset%icoulomb==0) then
     call mkcore(dummy6,dyfr_dum,gr_dum,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&     ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,dtset%typat,ucvol,&
&     vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)
   else if (psps%usewvl==0.and.(psps%usepaw==1.or.dtset%icoulomb==1)) then
     call mkcore_alt(atindx1,dummy6,dyfr_dum,gr_dum,dtset%icoulomb,mpi_enreg,dtset%natom,&
&     nfft,dtset%nspden,nattyp,ntypat,ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,&
&     ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred,pawrad,pawtab,psps%usepaw)
   else if (psps%usewvl==1.and.psps%usepaw==1) then
#if defined HAVE_BIGDFT
!      call mkcore_wvl_old(atindx1,dummy6,dyfr_dum,wvl%descr%atoms%astruct%geocode,gr_dum,wvl%descr%h,&
! &         dtset%natom,nattyp,nfft,wvl%den%denspot%dpbox%nscatterarr(mpi_enreg%me_wvl,:),&
! &         dtset%nspden,ntypat,wvl%descr%Glr%d%n1,wvl%descr%Glr%d%n1i,wvl%descr%Glr%d%n2,&
! &         wvl%descr%Glr%d%n2i,wvl%descr%Glr%d%n3,wvl%den%denspot%dpbox%n3pi,n3xccc,option,&
! &         pawrad,pawtab,psps%gth_params%psppar,rprimd,ucvol_local,vxc,xccc3d,xred,&
! &         mpi_comm_wvl=mpi_enreg%comm_wvl)
     call mkcore_wvl(atindx1,dummy6,gr_dum,dtset%natom,nattyp,nfft,dtset%nspden,ntypat,&
&     n1xccc,n3xccc,option,pawrad,pawtab,rprimd,vxc,psps%xccc1d,xccc3d,&
&     psps%xcccrc,xred,wvl%den,wvl%descr,mpi_comm_wvl=mpi_enreg%comm_wvl)
#endif
   end if
   ABI_DEALLOCATE(gr_dum)
   ABI_DEALLOCATE(dyfr_dum)
   call timab(92,2,tsec)
   call timab(91,1,tsec)
 end if
 if (coretau_method==2) then
   call timab(91,2,tsec)
   call timab(92,1,tsec)
   option=1
   ABI_ALLOCATE(gr_dum,(3,dtset%natom))
   ABI_ALLOCATE(dyfr_dum,(3,3,dtset%natom))
   call mkcore_alt(atindx1,dummy6,dyfr_dum,gr_dum,dtset%icoulomb,mpi_enreg,dtset%natom,&
&   nfft,dtset%nspden,nattyp,ntypat,ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,&
&   ucvol,vxc,psps%xcccrc,psps%xccc1d,xcctau3d,xred,pawrad,pawtab,psps%usepaw,&
&   usekden=.true.)
   ABI_DEALLOCATE(gr_dum)
   ABI_DEALLOCATE(dyfr_dum)
   call timab(92,2,tsec)
   call timab(91,1,tsec)
 end if

!Adds the jellium potential to the local part of ionic potential
 if (dtset%jellslab/=0) then
   ABI_ALLOCATE(vjell,(nfft))
   ABI_ALLOCATE(rhojellg,(2,nfft))
   ABI_ALLOCATE(rhojellr,(nfft))
   option=1
   call jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,dtset%nspden,option,&
&   dtset%slabwsrad,rhojellg,rhojellr,rprimd,vjell,dtset%slabzbeg,dtset%slabzend)
!  Compute background-background energy
   call dotprod_vn(1,rhojellr,ebb,doti,nfft,nfftot,1,1,vjell,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
   ebb=half*ebb
!  Compute electrostatic energy between background and nuclei before adding vjell to vpsp
   call dotprod_vn(1,rhojellr,ebn,doti,nfft,nfftot,1,1,vpsp,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
!  Update e_ewald with ebb and ebn
   energies%e_ewald=energies%e_ewald+ebb+ebn
!  Compute gradient of ebn wrt tn
!  This is not yet coded for usewvl or icoulomb=1
   if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
     write(message,'(3a)')&
&     'The computation of forces due to jellium background',ch10,&
&     'has to be verified in the PAW formalism.'
     MSG_WARNING(message)

     ABI_ALLOCATE(grtn,(3,dtset%natom))
     optatm=0;optdyfr=0;opteltfr=0;optgr=1;optstr=0;optv=1;optn=0;optn2=1
     call atm2fft(atindx1,dummy_out1,vpsp,dummy_out2,dummy_out3,dummy_out4,dummy_in,&
&     gmet,gprimd,dummy_out5,grtn,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&     optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,&
&     rhojellg,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,dummy_in,dummy_in,dummy_in,dtset%vprtrb,psps%vlspl,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&     paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)

!    Update grewtn with gradient of ebn wrt tn
     do iatom=1,dtset%natom
       grewtn(1:3,iatom)=grewtn(1:3,iatom)+grtn(1:3,iatom)
     end do
     ABI_DEALLOCATE(grtn)
   else ! of usepaw==1
     option=2
     ABI_ALLOCATE(dyfr_dum,(3,3,dtset%natom))
     ABI_ALLOCATE(grtn,(3,dtset%natom))
     call mklocl(dtset,dyfr_dum,energies%e_localpsp,gmet,gprimd,&
&     grtn,gsqcut,dummy6,mgfft,mpi_enreg,dtset%natom,nattyp,&
&     nfft,ngfft,1,ntypat,option,pawtab,ph1d,psps,dtset%qprtrb,rhojellg,&
&     rhojellr,rprimd,ucvol,dtset%vprtrb,vpsp,wvl%descr,wvl%den,xred)
!    Update grewtn with gradient of ebn wrt tn (reestablish order of atoms)
     do iatom=1,dtset%natom
       grewtn(1:3,atindx1(iatom))=grewtn(1:3,atindx1(iatom))+grtn(1:3,iatom)
     end do
     ABI_DEALLOCATE(dyfr_dum)
     ABI_DEALLOCATE(grtn)
   end if ! of usepaw==1
   vpsp(:)=vpsp(:)+vjell(:)
   ABI_DEALLOCATE(vjell)
   ABI_DEALLOCATE(rhojellg)
   ABI_DEALLOCATE(rhojellr)
 end if

!Additional stuff for electron-positron calculation
!Compute the electronic/positronic local (Hartree) potential
 if (ipositron==1) vpsp=-vpsp

!If we are at the initialisation, or
!if the atom positions has changed and the non-linear core correction
!is included, or the rhor has changed, one needs to compute the xc stuff.
!One needs also to compute the Hartree stuff if the density changed,
!or at initialisation.
!--------------------------------------------------------------

 if(istep==1 .or. n1xccc/=0 .or. moved_rhor==1 .or. dtset%positron<0 .or. mod(dtset%fockoptmix,100)==11) then

   option=0
   if(istep==1 .or. moved_rhor==1 .or. dtset%positron<0 .or. mod(dtset%fockoptmix,100)==11) option=1
   if (nkxc>0) option=2
   if (dtset%iscf==-1) option=-2
!  Not yet able to deal fully with the full XC kernel in case of GGA + spin
   option_eff=option
   if (option==2.and.xcdata%xclevel==2.and.(nkxc==3-2*mod(xcdata%nspden,2))) option_eff=12

   if (ipositron/=1) then
     if (dtset%icoulomb == 0 .and. dtset%usewvl == 0) then
       if(option/=0 .and. option/=10)then
         call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfft,ngfft,rhog,rprimd,vhartr)
       end if
       call xcdata_init(xcdata,dtset=dtset)
       if(mod(dtset%fockoptmix,100)==11)then
         xcdatahyb=xcdata
!        Setup the auxiliary xc functional information
         call xcdata_init(xcdata,dtset=dtset,auxc_ixc=0,ixc=dtset%auxc_ixc)
       end if
!      Use the periodic solver to compute Hxc
       nk3xc=1
       if (ipositron==0) then

!        Compute energies%e_xc and associated quantities
         if(.not.is_hybrid_ncpp .or. mod(dtset%fockoptmix,100)==11)then
           call rhotoxc(energies%e_xc,kxc,mpi_enreg,nfft,ngfft,&
&           nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,&
&           option_eff,rhor,rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata,&
&           taur=taur,vhartr=vhartr,vxctau=vxctau,add_tfw=add_tfw_,xcctau3d=xcctau3d)
         else
!          Only when is_hybrid_ncpp, and moreover, the xc functional is not the auxiliary xc functional, then call xchybrid_ncpp_cc
           call xchybrid_ncpp_cc(dtset,energies%e_xc,mpi_enreg,nfft,ngfft,n3xccc,rhor,rprimd,&
&           strsxc,vxcavg,xccc3d,vxc=vxc)
         end if

!        Possibly compute energies%e_hybcomp_E0
         if(mod(dtset%fockoptmix,100)==11)then
!          This call to rhotoxc uses the hybrid xc functional
           if(.not.is_hybrid_ncpp)then
             call rhotoxc(energies%e_hybcomp_E0,kxc,mpi_enreg,nfft,ngfft,&
&             nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,&
&             option,rhor,rprimd,strsxc,usexcnhat,vxc_hybcomp,vxcavg,xccc3d,xcdatahyb,&
&             taur=taur,vhartr=vhartr,vxctau=vxctau,add_tfw=add_tfw_)
           else
             call xchybrid_ncpp_cc(dtset,energies%e_hybcomp_E0,mpi_enreg,nfft,ngfft,n3xccc,rhor,rprimd,&
&             strsxc,vxcavg,xccc3d,vxc=vxc_hybcomp)
           end if

!          Combine hybrid and auxiliary quantities
           energies%e_xc=energies%e_xc*dtset%auxc_scal
           energies%e_hybcomp_E0=energies%e_hybcomp_E0-energies%e_xc
           vxc(:,:)=vxc(:,:)*dtset%auxc_scal
           vxc_hybcomp(:,:)=vxc_hybcomp(:,:)-vxc(:,:)

         end if

       else if (ipositron==2) then
         call rhotoxc(energies%e_xc,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,&
&         option,rhor,rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata,&
&         taur=taur,vhartr=vhartr,vxctau=vxctau,add_tfw=add_tfw_,&
&         electronpositron=electronpositron)
       end if

     elseif(.not. wvlbigdft) then
!      Use the free boundary solver
       call psolver_rhohxc(energies%e_hartree, energies%e_xc, evxc, &
&       dtset%icoulomb, dtset%ixc, &
&       mpi_enreg, nfft, ngfft,&
&       nhat,psps%usepaw,&
&       dtset%nscforder,dtset%nspden,n3xccc,rhor,rprimd, &
&       usexcnhat,psps%usepaw,dtset%usewvl,vhartr, vxc, &
&       vxcavg,wvl%descr,wvl%den,wvl%e,&
&       xccc3d,dtset%xclevel,dtset%xc_denpos)
     end if
   else
     energies%e_xc=zero
     call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat,nkxc,dtset%nspden,n3xccc,&
&     dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,psps%usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
   end if
   if (ipositron/=0) then
     if (optene>=1) then
       call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,&
&       nfft,nfftot,1,1,electronpositron%vha_ep,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
     end if
     vhartr=vhartr+electronpositron%vha_ep
   end if
 end if

!Compute the trial potential
!-------------------------------------------------------------
 if (.not. wvlbigdft) then
!  Now, compute trial Hxc potential. Local psp potential will be added later.
   if(moved_atm_inside==0 .or.dtset%iscf>=10) then

!    Compute starting Hxc potential.
!    Multiply by identity, should not change anything if nspden /= 4
     do ispden=1,dtset%nspden
       vtrial(:,ispden)=vhartr(:)*identity(ispden)+vxc(:,ispden)
     end do

   else

!    One should be here only when moved_atm_inside==1
!    The (H)xc now added corrects the previous one.
     if(dtset%densfor_pred==1)then
!      xc was substracted off. This should be rationalized later
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+vxc(:,ispden)
       end do
     else if(abs(dtset%densfor_pred)==2.or.abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)then
!      Hxc was substracted off. This should be rationalized later
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+vhartr(:)*identity(ispden)+vxc(:,ispden)
       end do
     end if
   end if

!  Adds the local part of the potential
   if ((moved_atm_inside==0).or.(dtset%densfor_pred/=3)) then
     do ispden=1,min(2,dtset%nspden)
       do ifft=1,nfft
         vtrial(ifft,ispden)=vtrial(ifft,ispden)+vpsp(ifft)
       end do
     end do
   end if

!  Adds the compensating vxc for hybrids
   if(mod(dtset%fockoptmix,100)==11)then
     vtrial(:,:)=vtrial(:,:)+vxc_hybcomp(:,:)
   end if

   if(dtset%usewvl==1) then
     call wvl_vtrial_abi2big(1,vtrial,wvl%den)
   end if

 else

!  Compute with covering comms the different part of the potential.
#if defined HAVE_BIGDFT
!  Copy e_ewald.
   wvl%e%energs%eion = energies%e_ewald
!  Setup the mixing, if necessary
   call denspot_set_history(wvl%den%denspot,dtset%iscf,dtset%nsppol, &
&   wvl%den%denspot%dpbox%ndims(1),wvl%den%denspot%dpbox%ndims(2))
#endif
   ABI_ALLOCATE(xcart,(3, dtset%natom))
   call xred2xcart(dtset%natom, rprimd, xcart, xred)
   call wvl_psitohpsi(dtset%diemix,energies%e_exactX, energies%e_xc, energies%e_hartree, &
&   energies%e_kinetic, energies%e_localpsp, energies%e_nlpsp_vfock, energies%e_sicdc, &
&   istep, 1, dtset%iscf, mpi_enreg%me_wvl, dtset%natom, dtset%nfft, mpi_enreg%nproc_wvl, dtset%nspden, &
&   rpnrm, .true.,evxc, wvl,.true., xcart, strsxc,&
&   vtrial, vxc)
   if (optene==3.or.optene==4) energies%e_xcdc=evxc
   ABI_DEALLOCATE(xcart)

 end if

!Add the zeeman field to vtrial
 if (any(abs(dtset%zeemanfield(:))>tol8)) then
   vzeeman(:) = zero                            ! vzeeman_ij = -1/2*sigma_ij^alpha*B_alpha
   if(dtset%nspden==2)then
     vzeeman(1) = -half*dtset%zeemanfield(3)   ! v_dwndwn = -1/2*B_z
     vzeeman(2) =  half*dtset%zeemanfield(3)   ! v_upup   =  1/2*B_z
     do ifft=1,nfft
       vtrial(ifft,1) = vtrial(ifft,1) + vzeeman(1) !SPr: added 1st component
       vtrial(ifft,2) = vtrial(ifft,2) + vzeeman(2)
     end do !ifft
   end if
   if(dtset%nspden==4)then
     vzeeman(1)=-half*dtset%zeemanfield(3)     ! v_dwndwn                  => v_11
     vzeeman(2)= half*dtset%zeemanfield(3)     ! v_upup                    => v_22
     vzeeman(3)=-half*dtset%zeemanfield(1)     ! Re(v_dwnup) = Re(v_updwn) => Re(v_12)
     vzeeman(4)= half*dtset%zeemanfield(2)     ! Im(v_dwnup) =-Im(v_dwnup) => Im(v_12)
     do ispden=1,dtset%nspden
       do ifft=1,nfft
         vtrial(ifft,ispden) = vtrial(ifft,ispden) + vzeeman(ispden)
       end do
     end do
   end if
 end if

!Compute the constrained potential for the magnetic moments
 if (dtset%magconon==1.or.dtset%magconon==2) then
!  Initialize the datastructure constrained_dft, for penalty function constrained magnetization
   call constrained_dft_ini(dtset%chrgat,constrained_dft,dtset%constraint_kind,dtset%magconon,dtset%magcon_lambda,&
&    mpi_enreg,dtset%natom,dtset%nfft,dtset%ngfft,dtset%nspden,dtset%ntypat,dtset%ratsm,&
&    dtset%ratsph,rprimd,dtset%spinat,dtset%typat,xred,dtset%ziontypat)
   ABI_ALLOCATE(v_constr_dft_r, (nfft,dtset%nspden))
   v_constr_dft_r = zero
   call mag_penalty(constrained_dft,mpi_enreg,rhor,v_constr_dft_r,xred)
   if(dtset%nspden==4)then
     do ispden=1,dtset%nspden ! (SPr: both components should be used? EB: Yes it should be the case, corrected now)
       do ifft=1,nfft
         vtrial(ifft,ispden) = vtrial(ifft,ispden) + v_constr_dft_r(ifft,ispden)
       end do !ifft
     end do !ispden
   else if(dtset%nspden==2)then
     do ifft=1,nfft
!      TODO : MJV: check that magnetic constraint works also for nspden 2 or add input variable condition
!              EB: ispden=2 is rho_up only: to be tested
!             SPr: for ispden=2, both components should be used (e.g. see definition for vzeeman)?
       vtrial(ifft,1) = vtrial(ifft,1) + v_constr_dft_r(ifft,1) !SPr: added the first component here
       vtrial(ifft,2) = vtrial(ifft,2) + v_constr_dft_r(ifft,2)
     end do !ifft
   end if
   ABI_DEALLOCATE(v_constr_dft_r)
   call constrained_dft_free(constrained_dft)
 end if

!Compute parts of total energy depending on potentials
!--------------------------------------------------------------

!For icoulomb==0 or usewvl Ehartree is calculated in psolver_rhohxc().
!For PAW we recalculate this since nhat was not taken into account
!in psolver_rhohxc: E_H= int v_H (n+nhat) dr

 if (optene>=1 .and. .not. wvlbigdft .and. (dtset%icoulomb==0 .or. dtset%usepaw==1 ) ) then
!  Compute Hartree energy ehart
!  Already available in the Psolver case through psolver_rhohxc().
   if (ipositron/=1) then
     call dotprod_vn(1,rhor,energies%e_hartree,doti,nfft,nfftot,1,1,vhartr,ucvol_local,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
     if (ipositron==0) energies%e_hartree = half * energies%e_hartree
     if (ipositron==2) energies%e_hartree = half * (energies%e_hartree-electronpositron%e_hartree)
   else
     energies%e_hartree=zero
   end if
 end if

 if(mod(dtset%fockoptmix,100)==11)then
   if (.not. wvlbigdft) then
     call dotprod_vn(1,rhor,energies%e_hybcomp_v0,doti,nfft,nfftot,1,1,vxc_hybcomp,ucvol_local,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
     energies%e_hybcomp_v=energies%e_hybcomp_v0
   end if
 end if

 if (optene==2.or.optene==4 .and. .not. wvlbigdft) then
!  Compute local psp energy eei
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,nfft,nfftot,1,1,vpsp,ucvol_local,&
&   mpi_comm_sphgrid=mpi_comm_sphgrid)
 end if

 if (optene==3.or.optene==4 .and. .not. wvlbigdft) then
!  Compute double-counting XC energy enxcdc
   if (ipositron/=1) then
     if (dtset%usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local,&
&       mpi_comm_sphgrid=mpi_comm_sphgrid)
     else
       ABI_ALLOCATE(rhowk,(nfft,dtset%nspden))
       rhowk=rhor-nhat
       call dotprod_vn(1,rhowk,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local,&
&                      mpi_comm_sphgrid=mpi_comm_sphgrid)
       ABI_DEALLOCATE(rhowk)
     end if
     if (with_vxctau) then
       call dotprod_vn(1,taur,e_xcdc_vxctau,doti,nfft,nfftot,dtset%nspden,1,vxctau(:,:,1),&
&                      ucvol_local,mpi_comm_sphgrid=mpi_comm_sphgrid)
       energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if
 end if

!--------------------------------------------------------------

!The initialisation for the new atomic positions has been done
 moved_atm_inside=0

 call timab(91,2,tsec)

end subroutine setvtr
!!***

!!****m* m_setvtr/spatialchempot
!! NAME
!!  spatialchempot
!! FUNCTION
!!  Treat spatially varying chemical potential.
!!  Compute energy and derivative with respect to dimensionless reduced atom coordinates of the
!!  spatially varying chemical potential. No contribution to stresses.
!!
!! INPUTS
!! chempot(3,nzchempot,ntypat)=input array with information about the chemical potential (see input variable description)
!! natom=number of atoms in unit cell
!! ntypat=number of type of atoms
!! nzchempot=number of limiting planes for chemical potential
!! typat(natom)=integer label of each type of atom (1,2,...)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!
!! OUTPUT
!! e_chempot=chemical potential energy in hartrees
!! grchempottn(3,natom)=grads of e_chempot wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!
!! SOURCE

subroutine spatialchempot(e_chempot,chempot,grchempottn,natom,ntypat,nzchempot,typat,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,nzchempot
 real(dp),intent(out) :: e_chempot
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: chempot(3,nzchempot,ntypat),xred(3,natom)
 real(dp),intent(out) :: grchempottn(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,itypat,iz
 real(dp) :: a_2,a_3,cp0,cp1,dcp0,dcp1,ddz,deltaz,deltaziz
 real(dp) :: dqz,dz1,qz,zred,z0
!character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,'(a)')' spatialchempot : enter '
!write(std_out,'(a,2i6)')' nzchempot,ntypat',nzchempot,ntypat
!write(std_out,'(a,6es13.3)') ' chempot(1:3,1:2,1)=',chempot(1:3,1:2,1)
!write(std_out,'(a,6es13.3)') ' chempot(1:3,1:2,2)=',chempot(1:3,1:2,2)
!ENDDEBUG

 e_chempot=zero
 grchempottn(:,:)=zero

!Loop on the different atoms
 do iatom=1,natom

   itypat=typat(iatom)
   zred=xred(3,iatom)

!  Determine the delimiting plane just lower to zred
!  First compute the relative zred with respect to the first delimiting plane
!  Take into account a tolerance :
   deltaz=zred-chempot(1,1,itypat)
   deltaz=modulo(deltaz+tol12,1.0d0)-tol12
!  deltaz is positive (or higher than -tol12), and lower than one-tol12.
   do iz=2,nzchempot+1
     if(iz/=nzchempot+1)then
       deltaziz=chempot(1,iz,itypat)-chempot(1,1,itypat)
     else
       deltaziz=one
     end if
     if(deltaziz>deltaz)exit
   end do

!  Defines coordinates and values inside the delimiting interval,
!  with respect to the lower delimiting plane
   z0=chempot(1,iz-1,itypat)-chempot(1,1,itypat) ; cp0=chempot(2,iz-1,itypat) ; dcp0=chempot(3,iz-1,itypat)
   if(iz/=nzchempot+1)then
     dz1=chempot(1,iz,itypat)-chempot(1,iz-1,itypat) ; cp1=chempot(2,iz,itypat) ; dcp1=chempot(3,iz,itypat)
   else
     dz1=(chempot(1,1,itypat)+one)-chempot(1,nzchempot,itypat) ; cp1=chempot(2,1,itypat) ; dcp1=chempot(3,1,itypat)
   end if
   ddz=deltaz-z0

!DEBUG
!  write(std_out,'(a,2i5)')' Delimiting planes, iz-1 and iz=',iz-1,iz
!  write(std_out,'(a,2es13.3)')' z0,  dz1= :',z0,dz1
!  write(std_out,'(a,2es13.3)')' cp0, cp1= :',cp0,cp1
!  write(std_out,'(a,2es13.3)')' dcp0, dcp1= :',dcp0,dcp1
!  write(std_out,'(a,2es13.3)')' deltaz,ddz=',deltaz,ddz
!ENDDEBUG

!  Determine the coefficient of the third-order polynomial taking z0 as origin
!  P(dz=z-z0)= a_3*dz**3 + a_2*dz**2 + a_1*dz + a_0 ; obviously a_0=cp0 and a_1=dcp0
!  Define qz=a_3*dz + a_2 and dqz=3*a_3*dz + 2*a_2
   qz=((cp1-cp0)-dcp0*dz1)/dz1**2
   dqz=(dcp1-dcp0)/dz1
   a_3=(dqz-two*qz)/dz1
   a_2=three*qz-dqz

!  Compute value and gradient of the chemical potential, at ddz wrt to lower delimiting plane
   e_chempot=e_chempot+(a_3*ddz**3 + a_2*ddz**2 + dcp0*ddz + cp0)
   grchempottn(3,iatom)=three*a_3*ddz**2 + two*a_2*ddz + dcp0

!DEBUG
!  write(std_out,'(a,4es16.6)')' qz,dqz=',qz,dqz
!  write(std_out,'(a,4es16.6)')' cp0,dcp0,a_2,a_3=',cp0,dcp0,a_2,a_3
!  write(std_out,'(a,2es13.3)')' dcp0*ddz + cp0=',dcp0*ddz + cp0
!  write(std_out,'(a,2es13.3)')' a_2*ddz**2=',a_2*ddz**2
!  write(std_out,'(a,2es13.3)')' a_3*ddz**3=',a_3*ddz**3
!  write(std_out,'(a,2es13.3)')' contrib=',a_3*ddz**3 + a_2*ddz**2 + dcp0*ddz + cp0
!  write(std_out,'(a,2es13.3)')' e_chempot=',e_chempot
!  write(std_out,'(a,3es20.10)')' grchempottn=',grchempottn(:,iatom)
!ENDDEBUG

 end do

!DEBUG
!write(std_out,'(a)')' spatialchempot : exit '
!write(std_out,'(a,es16.6)') ' e_chempot=',e_chempot
!ENDDEBUG

end subroutine spatialchempot
!!***

!!****f* ABINIT/ionion_realspace
!!
!! NAME
!! ionion_realspace
!!
!! FUNCTION
!! Compute the ion/ion interaction energies and forces in real space
!! case. Use ewald() instead if computations are done in reciprocal
!! space since it also includes the correction for the shift done in
!! potentials calculations and includes replica interactions.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  rmet(3,3)=metric tensor in real space (bohr^2)
!!  xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!  zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!!  eew=final ion/ion energy in hartrees
!!  grewtn(3,natom)=grads of ion/ion wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine ionion_realSpace(dtset, eew, grewtn, rprimd, xred, zion)

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: eew
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: rprimd(3,3),zion(dtset%ntypat)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ia1,ia2,iatom,igeo
 real(dp) :: r
!arrays
 real(dp),allocatable :: grew_cart(:,:),xcart(:,:)

! *************************************************************************

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

!Summing the interaction between ions.
 eew = 0._dp
 do ia1 = 1, dtset%natom, 1
   do ia2 = ia1 + 1, dtset%natom, 1
     r = sqrt((xcart(1, ia1) - xcart(1, ia2)) ** 2 + &
&     (xcart(2, ia1) - xcart(2, ia2)) ** 2 + &
&     (xcart(3, ia1) - xcart(3, ia2)) ** 2)
     eew = eew + zion(dtset%typat(ia1)) * zion(dtset%typat(ia2)) / r
   end do
 end do

!Allocate temporary array to store cartesian gradients.
 ABI_ALLOCATE(grew_cart,(3, dtset%natom))

!Summing the forces for each atom
 do ia1 = 1, dtset%natom, 1
   grew_cart(:, ia1) = 0._dp
   do ia2 = 1, dtset%natom, 1
     if (ia1 /= ia2) then
       r = (xcart(1, ia1) - xcart(1, ia2)) ** 2 + &
&       (xcart(2, ia1) - xcart(2, ia2)) ** 2 + &
&       (xcart(3, ia1) - xcart(3, ia2)) ** 2
       do igeo = 1, 3, 1
         grew_cart(igeo, ia1) = grew_cart(igeo, ia1) - (xcart(igeo, ia1) - xcart(igeo, ia2)) * &
&         zion(dtset%typat(ia1)) * zion(dtset%typat(ia2)) / (r ** 1.5_dp)
       end do
     end if
   end do
 end do

 ABI_DEALLOCATE(xcart)

!Transform cartesian gradients to reduced gradients.
 do iatom = 1, dtset%natom, 1
   do igeo = 1, 3, 1
     grewtn(igeo, iatom) = rprimd(1, igeo) * grew_cart(1, iatom) + &
&     rprimd(2, igeo) * grew_cart(2, iatom) + &
&     rprimd(3, igeo) * grew_cart(3, iatom)
   end do
 end do
 ABI_DEALLOCATE(grew_cart)

end subroutine ionion_realSpace
!!***

!!****f* ABINIT/ionion_surface
!!
!! NAME
!! ionion_surface
!!
!! FUNCTION
!! Compute the ion/ion interaction energies and forces in real space
!! case. Use ewald() instead if computations are done in reciprocal
!! space since it also includes the correction for the shift done in
!! potentials calculations and includes replica interactions.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  rmet(3,3)=metric tensor in real space (bohr^2)
!!  xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!  zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!!  eew=final ion/ion energy in hartrees
!!  grewtn(3,natom)=grads of ion/ion wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!      ionicenergyandforces,xred2xcart
!!
!! SOURCE

subroutine ionion_surface(dtset, eew, grewtn, me, nproc, rprimd, wvl, wvl_den, xred)

#if defined HAVE_BIGDFT
 use BigDFT_API, only: IonicEnergyandForces
#endif

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: me, nproc
 real(dp),intent(out) :: eew
 type(dataset_type),intent(in) :: dtset
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: dispersion, iatom, igeo
 real(dp) :: psoffset
!arrays
 real(dp),allocatable :: xcart(:,:)
 real(dp),pointer :: grew_cart(:,:),fdisp(:,:)
#if defined HAVE_BIGDFT
 real(dp) :: edisp
 real(dp) :: ewaldstr(6)
#endif

! *************************************************************************

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

 nullify(fdisp)
 nullify(grew_cart)
 dispersion = 0
 psoffset = 0._dp
#if defined HAVE_BIGDFT
 call IonicEnergyandForces(me, nproc, wvl_den%denspot%dpbox,&
& wvl%atoms, dtset%efield, xcart, &
& eew, grew_cart, dispersion, edisp, fdisp,&
& ewaldstr,wvl%Glr%d%n1,wvl%Glr%d%n2,wvl%Glr%d%n3,&
& wvl_den%denspot%V_ext, wvl_den%denspot%pkernel,psoffset)

 if (associated(fdisp)) then
   ABI_DEALLOCATE(fdisp)
 end if
#endif

 ABI_DEALLOCATE(xcart)

!Transform cartesian gradients to reduced gradients.
 do iatom = 1, dtset%natom, 1
   do igeo = 1, 3, 1
     grewtn(igeo, iatom) = -rprimd(1, igeo) * grew_cart(1, iatom) - &
&     rprimd(2, igeo) * grew_cart(2, iatom) - &
&     rprimd(3, igeo) * grew_cart(3, iatom)
   end do
 end do
 if (associated(grew_cart)) then
   ABI_DEALLOCATE(grew_cart)
 end if

#if !defined HAVE_BIGDFT
 if (.false.) write(std_out,*) me,nproc,wvl%h(1),wvl_den%symObj
#endif

end subroutine ionion_surface
!!***

end module m_setvtr
!!***

!!****m* ABINIT/m_dft_energy
!! NAME
!!  m_dft_energy
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, AR, MB, MT, EB)
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

module m_dft_energy

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_hamiltonian
 use m_errors
 use m_xmpi
 use m_gemm_nonlop
 use m_xcdata
 use m_cgtools
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_geometry,         only : metric
 use m_kg,               only : mkkin
 use m_energies,         only : energies_type
 use m_electronpositron, only : electronpositron_type, electronpositron_calctype, rhohxcpositron
 use m_bandfft_kpt,      only : bandfft_kpt, bandfft_kpt_type, prep_bandfft_tabs, &
                                bandfft_kpt_savetabs, bandfft_kpt_restoretabs
 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_init_unpacked, &
                                pawrhoij_mpisum_unpacked, pawrhoij_free_unpacked, pawrhoij_inquire_dim, &
&                               pawrhoij_symrhoij
 use m_pawcprj,          only : pawcprj_type,pawcprj_alloc,pawcprj_free,pawcprj_gather_spin
 use m_pawfgr,           only : pawfgr_type
 use m_paw_dmft,         only : paw_dmft_type
 use m_paw_nhat,         only : pawmknhat
 use m_paw_occupancies,  only : pawaccrhoij
 use m_fft,              only : fftpac, fourdp
 use m_spacepar,         only : meanvalue_g, hartre
 use m_dens,             only : constrained_dft_t,mag_penalty
 use m_mkrho,            only : mkrho
 use m_mkffnl,           only : mkffnl
 use m_getghc,           only : getghc
 use m_rhotoxc,          only : rhotoxc
 use m_mpinfo,           only : proc_distrb_cycle
 use m_nonlop,           only : nonlop
 use m_fourier_interpol, only : transgrid
 use m_prep_kgb,         only : prep_getghc, prep_nonlop
 use m_psolver,          only : psolver_rhohxc

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

 implicit none

 private
!!***

 public :: energy
!!***

contains
!!***

!!****f* ABINIT/energy
!! NAME
!! energy
!!
!! FUNCTION
!! Compute electronic energy terms
!! energies%e_eigenvalues, ek and enl from arbitrary (orthonormal) provided wf,
!! ehart, enxc, and eei from provided density and potential,
!! energies%e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!! energies%e_zeeman=Zeeman spin energy from applied magnetic field -m.B
!! ek=kinetic energy, ehart=Hartree electron-electron energy,
!! enxc,enxcdc=exchange-correlation energies, eei=local pseudopotential energy,
!! enl=nonlocal pseudopotential energy
!! Also, compute new density from provided wfs, after the evaluation
!! of ehart, enxc, and eei.
!! WARNING XG180913 : At present, Fock energy not computed !
!!
!! INPUTS
!!  [add_tfw]=flag controling the addition of Weiszacker gradient correction to Thomas-Fermi kin energy
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=number of k points treated by this node.
!!   | mpw=maximum dimension for number of planewaves
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for polarized
!!   | nspinor=number of spinorial components
!!   | nsym=number of symmetry elements in space group (at least 1)
!!   | occopt=option for occupancies
!!   | tsmear=smearing energy or temperature (if metal)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut=G^2 cutoff from gsqcut=ecut/(2 Pi^2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  optene=option for the computation of total energy (direct scheme or double-counting scheme)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | ntypat=number of types of atoms in cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfftf)=local pseudopotential in real space (hartree)
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!  wvl <type(wvl_internal_type)>=wavelets internal data
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  compch_fft=-PAW only- compensation charge inside spheres computed over fine fft grid
!!  etotal=total energy (hartree):
!!    - computed by direct scheme if optene=0 or 2
!!    - computed by double-counting scheme if optene=1 or 3
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points (hartree^2)
!!  strsxc(6)=exchange-correlation contribution to stress tensor
!!  vhartr(nfftf)=work space to hold Hartree potential in real space (hartree)
!!  vtrial(nfftf,nspden)=total local potential (hartree)
!!  vxc(nfftf,nspden)=work space to hold Vxc(r) in real space (hartree)
!!  [vxctau(nfft,nspden,4*usekden)]=(only for meta-GGA): derivative of XC energy density
!!    with respect to kinetic energy density (depsxcdtau). The arrays vxctau contains also
!!    the gradient of vxctau (gvxctau) in vxctau(:,:,2:4)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_chempot(IN)=energy from spatially varying chemical potential (hartree)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_vdw_dftd(IN)=VdW DFT-D energy
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_eigenvalues(OUT)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_hartree(OUT)=Hartree part of total energy (hartree units)
!!   | e_kinetic(OUT)=kinetic energy part of total energy.
!!   | e_nlpsp_vfock(OUT)=nonlocal psp + potential Fock ACE part of total energy.
!!   | e_xc(OUT)=exchange-correlation energy (hartree)
!!  ==== if optene==0, 2 or 3
!!   | e_localpsp(OUT)=local psp energy (hartree)
!!  ==== if optene==1, 2 or 3
!!   | e_xcdc(OUT)=exchange-correlation double-counting energy (hartree)
!!  rhog(2,nfftf)=work space for rho(G); save intact on return (? MT 08-12-2008: is that true now ?)
!!  rhor(nfftf,nspden)=work space for rho(r); save intact on return (? MT 08-12-2008: is that true now ?)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  nspinor should not be modified in the call of rdnpw
!!  === if psps%usepaw==1 ===
!!    nhat(nfftf,nspden*usepaw)= compensation charge density
!!    pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!!  There is a large amount of overhead in the way this routine do the computation of the energy !
!!  For example, the density has already been precomputed, so why to compute it again here ??
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      bandfft_kpt_restoretabs,bandfft_kpt_savetabs,destroy_hamiltonian
!!      dotprod_vn,fftpac,fourdp,gpu_finalize_ffnl_ph3d,gpu_update_ffnl_ph3d
!!      hartre,init_hamiltonian,load_k_hamiltonian,load_spin_hamiltonian
!!      mag_penalty,make_gemm_nonlop,meanvalue_g,metric,mkffnl,mkkin,mkresi
!!      mkrho,nonlop,pawaccrhoij,pawcprj_alloc,pawcprj_free,pawcprj_gather_spin
!!      pawmknhat,pawrhoij_alloc,pawrhoij_free,pawrhoij_free_unpacked
!!      pawrhoij_init_unpacked,pawrhoij_mpisum_unpacked,prep_bandfft_tabs
!!      prep_nonlop,psolver_rhohxc,rhohxcpositron,rhotoxc,pawrhoij_symrhoij,timab
!!      transgrid,xcdata_init,xmpi_sum
!!
!! SOURCE

subroutine energy(cg,compch_fft,constrained_dft,dtset,electronpositron,&
& energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,my_natom,nfftf,ngfftf,nhat,&
& nhatgr,nhatgrdim,npwarr,n3xccc,occ,optene,paw_dmft,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrhoij,pawtab,phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,&
& taug,taur,usexcnhat,vhartr,vtrial,vpsp,vxc,wfs,wvl,wvl_den,wvl_e,xccc3d,xred,ylm,&
& add_tfw,vxctau) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,my_natom,n3xccc,nfftf,nhatgrdim,optene,usexcnhat
 logical,intent(in),optional :: add_tfw
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: compch_fft,etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(constrained_dft_t),intent(in) :: constrained_dft
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_wf_type),intent(inout) :: wfs
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(wvl_energy_terms),intent(inout) ::wvl_e
!arrays
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)
 integer :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)),kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: ngfftf(18),npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: cg(2,mcg),eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp), intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden),taur(nfftf,dtset%nspden*dtset%usekden)
 real(dp), intent(out) :: strsxc(6)
 real(dp), intent(in) :: rprimd(3,3),vpsp(nfftf),xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp), intent(out) :: vhartr(nfftf),vtrial(nfftf,dtset%nspden),vxc(nfftf,dtset%nspden)
 real(dp),intent(out),optional,target :: vxctau(nfftf,dtset%nspden,4*dtset%usekden)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 type(paw_ij_type), intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,blocksize,choice,cplex,cplex_rhoij,cpopt,dimffnl
 integer :: iband,iband_last,iblock,iblocksize,icg,ider,idir,ierr,ifft,ikg,ikpt,ilm
 integer :: ipert,ipositron,iresid,ispden,isppol,istwf_k,izero
 integer :: me_distrb,mpi_comm_sphgrid,my_ikpt,my_nspinor,n1,n2,n3,n4,n5,n6
 integer :: nband_k,nblockbd,nfftotf,nkpg,nkxc,nk3xc,nnlout,npw_k,nspden_rhoij,option
 integer :: option_rhoij,paw_opt,signs,spaceComm,tim_mkrho,tim_nonlop
 logical :: add_tfw_,paral_atom,usetimerev,with_vxctau
 logical :: non_magnetic_xc,wvlbigdft=.false.
 real(dp) :: dotr,doti,eeigk,ekk,enlk,evxc,e_xcdc_vxctau,ucvol,ucvol_local,vxcavg
 !character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
 type(xcdata_type) :: xcdata
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpg_dum(0,0),kpoint(3),nonlop_out(1,1)
 real(dp) :: qpt(3),rhodum(1),rmet(3,3),tsec(2),ylmgr_dum(1,1,1),vzeeman(4)
 real(dp) :: magvec(dtset%nspden)
 real(dp),target :: vxctau_dum(0,0,0)
 real(dp),allocatable :: buffer(:),cgrvtrial(:,:)
 real(dp),allocatable :: cwavef(:,:),eig_k(:),enlout(:),ffnl(:,:,:,:),ffnl_sav(:,:,:,:)
 real(dp),allocatable :: kinpw(:),kinpw_sav(:),kxc(:,:),occ_k(:),occblock(:)
 real(dp),allocatable :: ph3d(:,:,:),ph3d_sav(:,:,:)
 real(dp),allocatable :: resid_k(:),rhowfg(:,:),rhowfr(:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: vlocal_tmp(:,:,:),vxctaulocal(:,:,:,:,:),ylm_k(:,:),v_constr_dft_r(:,:)
 real(dp),pointer :: vxctau_(:,:,:)
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 type(pawcprj_type),target,allocatable :: cwaveprj(:,:)
 type(pawcprj_type),pointer :: cwaveprj_gat(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = (present(vxctau).and.dtset%usekden/=0)
 vxctau_ => vxctau_dum ; if (with_vxctau) vxctau_ => vxctau

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 nfftotf=PRODUCT(ngfftf(1:3))
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   MSG_BUG('wrong values for nfft, nfftf!')
 end if

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1 .and. dtset%wvl_bigdft_comp==1)

 call timab(59,1,tsec)

!Data for parallelism
 spaceComm=mpi_enreg%comm_cell
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 if(mpi_enreg%paral_hf==1) spaceComm=mpi_enreg%comm_kpt
 mpi_comm_sphgrid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) then
   spaceComm=mpi_enreg%comm_wvl
   mpi_comm_sphgrid=mpi_enreg%comm_wvl
 end if
 me_distrb=mpi_enreg%me_kpt
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 paral_atom=(my_natom/=dtset%natom)

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 if (dtset%usewvl == 0) then
   ucvol_local = ucvol
#if defined HAVE_BIGDFT
 else
!  We need to tune the volume when wavelets are used because, not
!  all FFT points are used.
!  ucvol_local = (half * dtset%wvl_hgrid) ** 3 * ngfft(1)*ngfft(2)*ngfft(3)
   ucvol_local = product(wvl_den%denspot%dpbox%hgrids) * real(product(wvl_den%denspot%dpbox%ndims), dp)
#endif
 end if

!Compute Hxc potential from density
 option=1;nkxc=0
 ipositron=electronpositron_calctype(electronpositron)
 add_tfw_=.false.;if (present(add_tfw)) add_tfw_=add_tfw
 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

 if (ipositron/=1) then

   if (dtset%icoulomb == 0) then
!    Use the periodic solver to compute Hxc.
     call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfftf,ngfftf,rhog,rprimd,vhartr)
     call xcdata_init(xcdata,dtset=dtset)
     ABI_ALLOCATE(kxc,(1,nkxc))
!    to be adjusted for the call to rhotoxc
     nk3xc=1
     if (ipositron==0) then
       call rhotoxc(energies%e_xc,kxc, &
&       mpi_enreg,nfftf,ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim, &
&       nkxc,nk3xc,non_magnetic_xc,n3xccc,option,rhor,rprimd,strsxc, &
&       usexcnhat,vxc,vxcavg,xccc3d,xcdata,taur=taur,vhartr=vhartr, &
&       vxctau=vxctau_,exc_vdw_out=energies%e_xc_vdw,add_tfw=add_tfw_)
     else
       call rhotoxc(energies%e_xc,kxc, &
&       mpi_enreg,nfftf,ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim, &
&       nkxc,nk3xc,non_magnetic_xc,n3xccc,option,rhor,rprimd,strsxc, &
&       usexcnhat,vxc,vxcavg,xccc3d,xcdata, &
&       electronpositron=electronpositron,taur=taur,vhartr=vhartr, &
&       vxctau=vxctau_,exc_vdw_out=energies%e_xc_vdw,add_tfw=add_tfw_)
     end if
     ABI_DEALLOCATE(kxc)
   else if (dtset%usewvl == 0) then
!    Use the free boundary solver.
     call psolver_rhohxc(energies%e_hartree, energies%e_xc, evxc, &
&     dtset%icoulomb, dtset%ixc, mpi_enreg, nfftf, &
&     ngfftf,nhat,psps%usepaw,&
&     dtset%nscforder,dtset%nspden,n3xccc,rhor,rprimd, &
&     usexcnhat,psps%usepaw,dtset%usewvl,&
&     vhartr, vxc, vxcavg, wvl,wvl_den,wvl_e,&
&     xccc3d,dtset%xclevel,dtset%xc_denpos)
   end if
 else
   energies%e_xc=zero
   call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfftf,ngfftf,nhat,nkxc,dtset%nspden,n3xccc,&
&   dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,psps%usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
 end if
 if (ipositron/=0) then
   call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,nfftf,nfftotf,1,1,electronpositron%vha_ep,&
&   ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
   vhartr=vhartr+electronpositron%vha_ep
 end if

!Total local potential (for either spin channel) is
!Hartree + local psp + Vxc(spin), minus its mean
!(Note : this potential should agree with the input vtrial)
 do ispden=1,min(dtset%nspden,2)
   do ifft=1,nfftf
     vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
   end do
 end do
 if (dtset%nspden==4) then
   do ifft=1,nfftf
     vtrial(ifft,3:4)=vxc(ifft,3:4)
   end do
 end if

!Add the vzeeman pot in the trial pot
!Vzeeman might have to be allocated correctly --> to be checked
 if (any(abs(dtset%zeemanfield(:))>tol8)) then
   vzeeman(:) = zero
   if(dtset%nspden==2)then
!TODO: check this against rhotov and setvtr, where the potential is -1/2 and +1/2 for the 2 spin components.
! see comment by SPr in rhotov
! TODO: check this 1/2 factor is for the electron spin magnetic moment.
     vzeeman(1) = -half*dtset%zeemanfield(3) ! For collinear ispden=1 potential is v_upup
     vzeeman(2) = +half*dtset%zeemanfield(3) ! For collinear ispden=2 potential is v_dndn
   end if
   if(dtset%nspden==4)then
     vzeeman(1)=-half*dtset%zeemanfield(3)
     vzeeman(2)= half*dtset%zeemanfield(3)
     vzeeman(3)=-half*dtset%zeemanfield(1)
     vzeeman(4)= half*dtset%zeemanfield(2)
   end if
   magvec = zero
   do ispden=1,dtset%nspden
     do ifft=1,nfftf
!TODO: the full cell magnetization will need extra PAW terms, and is certainly calculated elsewhere.
!The calculation of the zeeman energy can be moved there
       magvec(ispden) = magvec(ispden) + rhor(ifft,ispden)
       vtrial(ifft,ispden)=vtrial(ifft,ispden)+vzeeman(ispden)
     end do
   end do
   if(dtset%nspden==2)then
     energies%e_zeeman = -half*dtset%zeemanfield(3)*(two*magvec(2)-magvec(1)) !  diff rho = rhoup-rhodown = 2 rhoup - rho
   else if(dtset%nspden==4)then
     energies%e_zeeman = -half * (dtset%zeemanfield(1)*magvec(2)& ! x
&                                +dtset%zeemanfield(2)*magvec(3)& ! y
&                                +dtset%zeemanfield(3)*magvec(4)) ! z
   end if
 end if

!Compute the constrained potential for the magnetic moments
!NOTE: here in energy.F90 rhor and vtrial are given on nfftf grid
!the values coming from mag_penalty may be different from those calculated
!calling mag_penalty with nfft in setvtr and rhotov
 if (dtset%magconon==1.or.dtset%magconon==2) then
   ABI_ALLOCATE(v_constr_dft_r, (nfftf,dtset%nspden))
   v_constr_dft_r = zero
   call mag_penalty(constrained_dft,mpi_enreg,rhor,v_constr_dft_r,xred)
!   call mag_penalty(dtset%natom, dtset%spinat, dtset%nspden, dtset%magconon, dtset%magcon_lambda, rprimd, &
!&   mpi_enreg, nfftf, dtset%ngfft, dtset%ntypat, dtset%ratsph, rhor, &
!&   dtset%typat, v_constr_dft_r, xred)
   do ispden=1,dtset%nspden
     do ifft=1,nfftf
       vtrial(ifft,ispden)=vtrial(ifft,ispden)+v_constr_dft_r(ifft,ispden)
     end do
   end do
   ABI_DEALLOCATE(v_constr_dft_r)
 end if

!Compute Hartree energy - use up+down rhor
 if (ipositron/=1) then
   call dotprod_vn(1,rhor,energies%e_hartree ,doti,nfftf,nfftotf,1,1,vhartr,&
&   ucvol_local,mpi_comm_sphgrid=mpi_comm_sphgrid)
   if (ipositron==0) energies%e_hartree=half*energies%e_hartree
   if (ipositron==2) energies%e_hartree = half *(energies%e_hartree-electronpositron%e_hartree)
 else
   energies%e_hartree=zero
 end if

!Compute local psp energy - use up+down rhor
 if (optene/=1) then
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,nfftf,nfftotf,1,1,vpsp,&
&   ucvol_local,mpi_comm_sphgrid=mpi_comm_sphgrid)
 end if

!Compute DC-xc energy - use up+down rhor
 if (optene>0) then
   if (ipositron/=1) then
     if (psps%usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,nfftf,nfftotf,dtset%nspden,1,vxc,&
&       ucvol_local,mpi_comm_sphgrid=mpi_comm_sphgrid)
       if (with_vxctau)then
         call dotprod_vn(1,taur,e_xcdc_vxctau,doti,nfftf,nfftotf,dtset%nspden,1,vxctau(:,:,1),&
&         ucvol_local,mpi_comm_sphgrid=mpi_comm_sphgrid)
         energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
       end if
     else
       ABI_ALLOCATE(rhowfr,(nfftf,dtset%nspden))
       rhowfr=rhor-nhat
       call dotprod_vn(1,rhowfr,energies%e_xcdc,doti,nfftf,nfftotf,dtset%nspden,1,vxc,&
&       ucvol_local,mpi_comm_sphgrid=mpi_comm_sphgrid)
       ABI_DEALLOCATE(rhowfr)
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if
 end if

 energies%e_eigenvalues=zero
 energies%e_kinetic=zero
 energies%e_nlpsp_vfock=zero
 energies%e_fock0=zero
 bdtot_index=0
 icg=0

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

!============================================
!==== Initialize most of the Hamiltonian ====
!============================================
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
!* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,&
& dtset%natom,dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& paw_ij=paw_ij,ph1d=ph1d,electronpositron=electronpositron,&
& nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamk%nvloc))
 if (with_vxctau) then
   ABI_ALLOCATE(vxctaulocal,(n4,n5,n6,gs_hamk%nvloc,4))
 end if

!PAW: additional initializations
 if (psps%usepaw==1) then
   ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,my_nspinor))
   call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
   if (mpi_enreg%paral_spinor==1) then
     ABI_DATATYPE_ALLOCATE(cwaveprj_gat,(dtset%natom,dtset%nspinor))
     call pawcprj_alloc(cwaveprj_gat,0,gs_hamk%dimcprj)
   else
     cwaveprj_gat => cwaveprj
   end if
   if (paral_atom) then
     ABI_DATATYPE_ALLOCATE(pawrhoij_unsym,(dtset%natom))
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&                nspden=dtset%nspden,spnorb=dtset%pawspnorb,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&     dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
   else
     pawrhoij_unsym => pawrhoij
     call pawrhoij_init_unpacked(pawrhoij_unsym)
   end if
   option_rhoij=1
   usetimerev=(dtset%kptopt>0.and.dtset%kptopt<3)
 else
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
 end if

!LOOP OVER SPINS
 do isppol=1,dtset%nsppol
   ikg=0

!  Set up local potential vlocal with proper dimensioning, from vtrial
!  Also take into account the spin.
   if(dtset%nspden/=4)then
     if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
       call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
       if(with_vxctau) then
         do ispden=1,4
           call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,&
&           vxctau(:,:,ispden),vxctaulocal(:,:,:,:,ispden),2)
         end do
       end if
     else
       ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
       if(with_vxctau) then
         do ispden=1,4
           call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,&
&                         rhodum,rhodum,cgrvtrial,vxctau(:,:,ispden))
           call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,&
&                      cgrvtrial,vxctaulocal(:,:,:,:,ispden),2)
         end do
       end if
       ABI_DEALLOCATE(cgrvtrial)
     end if
   else
     ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
     if (psps%usepaw==0) then
       do ispden=1,dtset%nspden
         call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal_tmp,2)
         vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       end do
     else
       ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       do ispden=1,dtset%nspden
         call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal_tmp,2)
         vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       end do
       ABI_DEALLOCATE(cgrvtrial)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
   end if

!  Continue Hamiltonian initialization
   call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)
   if (with_vxctau) then
     call gs_hamk%load_spin(isppol,vxctaulocal=vxctaulocal)
   end if

!  Loop over k points
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     istwf_k=dtset%istwfk(ikpt)
     npw_k=npwarr(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then
       resid(1+bdtot_index : nband_k+bdtot_index) = zero
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

!    Parallelism over FFT and/or bands: define sizes and tabs
     if (mpi_enreg%paral_kgb==1) then
       my_ikpt=mpi_enreg%my_kpttab(ikpt)
       nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
       my_bandfft_kpt => bandfft_kpt(my_ikpt)
     else
       my_ikpt=ikpt
       nblockbd=nband_k/mpi_enreg%bandpp
       !if (nband_k/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
     end if
     blocksize=nband_k/nblockbd

     ABI_ALLOCATE(eig_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(resid_k,(nband_k))
     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
     resid_k(:)=zero
     kpoint(:)=dtset%kptns(:,ikpt)
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
     eig_k(:)=eigen(1+bdtot_index:nband_k+bdtot_index)
     if (minval(eig_k)>1.d100) eig_k=zero
     eeigk=zero ; ekk=zero ; enlk=zero

     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Compute kinetic energy
     ABI_ALLOCATE(kinpw,(npw_k))
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

!    Compute kinetic energy of each band
     do iblock=1,nblockbd
       do iblocksize=1,blocksize
         iband=(iblock-1)*blocksize+iblocksize
         if (abs(occ_k(iband))>tol8) then
           cwavef(1:2,1:npw_k*my_nspinor)= &
&           cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg)
           call meanvalue_g(dotr,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,cwavef,cwavef,0)
           energies%e_kinetic=energies%e_kinetic+dtset%wtk(ikpt)*occ_k(iband)*dotr
         end if
       end do
     end do

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;dimffnl=1;nkpg=0
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,ider,psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&     npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&     psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)

!    Load k-dependent part in the Hamiltonian datastructure
!     - Compute 3D phase factors
!     - Prepare various tabs in case of band-FFT parallelism
!     - Load k-dependent quantities in the Hamiltonian
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
     call gs_hamk%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,&
&     kinpw_k=kinpw,kg_k=kg_k,ffnl_k=ffnl,ph3d_k=ph3d,&
&     compute_ph3d=.true.,compute_gbound=(mpi_enreg%paral_kgb/=1))

!    Load band-FFT tabs (transposed k-dependent arrays)
     if (mpi_enreg%paral_kgb==1) then
       call bandfft_kpt_savetabs(my_bandfft_kpt,ffnl=ffnl_sav,ph3d=ph3d_sav,kinpw=kinpw_sav)
       call prep_bandfft_tabs(gs_hamk,ikpt,dtset%mkmem,mpi_enreg)
       call gs_hamk%load_k(npw_fft_k=my_bandfft_kpt%ndatarecv, &
&       gbound_k =my_bandfft_kpt%gbound, &
&       kinpw_k  =my_bandfft_kpt%kinpw_gather, &
&       kg_k     =my_bandfft_kpt%kg_k_gather, &
&       ffnl_k   =my_bandfft_kpt%ffnl_gather, &
&       ph3d_k   =my_bandfft_kpt%ph3d_gather)
     end if

!    Setup gemm_nonlop
     if (gemm_nonlop_use_gemm) then
       gemm_nonlop_ikpt_this_proc_being_treated = my_ikpt
       call make_gemm_nonlop(my_ikpt,gs_hamk%npw_fft_k,gs_hamk%lmnmax, &
&       gs_hamk%ntypat, gs_hamk%indlmn, gs_hamk%nattyp, gs_hamk%istwf_k, gs_hamk%ucvol, gs_hamk%ffnl_k,&
&       gs_hamk%ph3d_k)
     end if

#if defined HAVE_GPU_CUDA
     if (dtset%use_gpu_cuda==1) then
       call gpu_update_ffnl_ph3d(ph3d,size(ph3d),ffnl,size(ffnl))
     end if
#endif

!    Compute nonlocal psp energy (NCPP) or Rhoij (PAW)
     ABI_ALLOCATE(enlout,(blocksize))
     ABI_ALLOCATE(occblock,(blocksize))
     do iblock=1,nblockbd
       iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband_k)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,me_distrb)) cycle

!      Select occupied bands
       occblock(:)=occ(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
       if(abs(maxval(occblock))>=tol8 ) then
         cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&         cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)

         choice=1-gs_hamk%usepaw ; signs=1 ; idir=0 ; nnlout=blocksize
         paw_opt=gs_hamk%usepaw;cpopt=gs_hamk%usepaw-1

         if (mpi_enreg%paral_kgb/=1) then
           tim_nonlop=3
           call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir,(/zero/),mpi_enreg,blocksize,nnlout,&
&           paw_opt,signs,nonlop_out,tim_nonlop,cwavef,cwavef)
         else
           tim_nonlop=14
           call prep_nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir,(/zero/),blocksize,&
&           mpi_enreg,nnlout,paw_opt,signs,nonlop_out,tim_nonlop,cwavef,cwavef)
         end if

         do iblocksize=1,blocksize
           iband=(iblock-1)*blocksize+iblocksize
           energies%e_eigenvalues=energies%e_eigenvalues+dtset%wtk(ikpt)*occ_k(iband)*eig_k(iband)
!          WARNING : the Fock contribution is NOT computed !!!
           energies%e_nlpsp_vfock=energies%e_nlpsp_vfock+dtset%wtk(ikpt)*occ_k(iband)*enlout(iblocksize)
         end do

!        PAW: accumulate rhoij
         if (psps%usepaw==1) then
           cplex=merge(1,2,istwf_k>1)
           if (mpi_enreg%paral_spinor==1) then
             call pawcprj_gather_spin(cwaveprj,cwaveprj_gat,dtset%natom,1,my_nspinor,dtset%nspinor,&
&             mpi_enreg%comm_spinor,ierr)
             call pawaccrhoij(gs_hamk%atindx,cplex,cwaveprj_gat,cwaveprj_gat,0,isppol,dtset%natom,dtset%natom,&
&             dtset%nspinor,occ_k(iband),option_rhoij,pawrhoij_unsym,usetimerev,dtset%wtk(ikpt))
           else
             call pawaccrhoij(gs_hamk%atindx,cplex,cwaveprj,cwaveprj,0,isppol,dtset%natom,dtset%natom,&
&             dtset%nspinor,occ_k(iband),option_rhoij,pawrhoij_unsym,usetimerev,dtset%wtk(ikpt))
           end if
         end if

!        End loop on bands
       end if
     end do

!    Compute residual of each band (for informative purposes)
     call mkresi(cg,eig_k,gs_hamk,icg,ikpt,isppol,mcg,mpi_enreg,nband_k,dtset%prtvol,resid_k)
     resid(1+bdtot_index : nband_k+bdtot_index) = resid_k(:)

!    Restore the bandfft tabs
     if (mpi_enreg%paral_kgb==1) then
       call bandfft_kpt_restoretabs(my_bandfft_kpt,ffnl=ffnl_sav,ph3d=ph3d_sav,kinpw=kinpw_sav)
     end if

!    Incremente indexes
     bdtot_index=bdtot_index+nband_k
     if (dtset%mkmem/=0) then
       icg=icg+npw_k*my_nspinor*nband_k
       ikg=ikg+npw_k
     end if

#if defined HAVE_GPU_CUDA
     if(dtset%use_gpu_cuda==1) then
       call gpu_finalize_ffnl_ph3d()
     end if
#endif

     ABI_DEALLOCATE(eig_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(resid_k)
     ABI_DEALLOCATE(enlout)
     ABI_DEALLOCATE(occblock)
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kinpw)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(ylm_k)

!    End loops on isppol and ikpt
   end do
 end do

 call gs_hamk%free()

 if(xmpi_paral==1)then
!  Accumulate enl eeig and ek on all proc.
   ABI_ALLOCATE(buffer,(3+dtset%mband*dtset%nkpt*dtset%nsppol))
   buffer(1)=energies%e_nlpsp_vfock ; buffer(2)=energies%e_kinetic ; buffer(3)=energies%e_eigenvalues
   do iresid=1,dtset%mband*dtset%nkpt*dtset%nsppol
     buffer(iresid+3)=resid(iresid)
   end do
   call timab(48,1,tsec)
   call xmpi_sum(buffer,spaceComm,ierr)
   call timab(48,2,tsec)
   energies%e_nlpsp_vfock=buffer(1) ; energies%e_kinetic=buffer(2) ; energies%e_eigenvalues=buffer(3)
   do iresid=1,dtset%mband*dtset%nkpt*dtset%nsppol
     resid(iresid)=buffer(iresid+3)
   end do
   ABI_DEALLOCATE(buffer)
!  Accumulate rhoij_
   if (psps%usepaw==1) then
     call pawrhoij_mpisum_unpacked(pawrhoij_unsym,spaceComm,comm2=mpi_enreg%comm_band)
   end if
 end if

!Compute total (free) energy
 if (optene==0.or.optene==2) then
   etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
!&   energies%e_nlpsp_vfock - energies%e_fock0 +
!   Should compute the e_fock0 energy !! Also, the Fock contribution to e_nlpsp_vfock
&   energies%e_nlpsp_vfock + energies%e_localpsp + energies%e_corepsp
   if (psps%usepaw==1) etotal=etotal + energies%e_paw
 else if (optene==1.or.optene==3) then
   etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&   energies%e_xcdc + energies%e_corepsp - energies%e_corepspdc
   if (psps%usepaw==1) etotal=etotal + energies%e_pawdc
 end if
 etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
 if(dtset%occopt>=3 .and. dtset%occopt<=8) etotal=etotal-dtset%tsmear*energies%entropy

!Additional stuff for electron-positron
 if (dtset%positron/=0) then
   if (ipositron==0) then
     energies%e_electronpositron  =zero
     energies%edc_electronpositron=zero
   else
     energies%e_electronpositron  =electronpositron%e_hartree+electronpositron%e_xc
     energies%edc_electronpositron=electronpositron%e_hartree+electronpositron%e_xcdc
     if (psps%usepaw==1) then
       energies%e_electronpositron  =energies%e_electronpositron  +electronpositron%e_paw
       energies%edc_electronpositron=energies%edc_electronpositron+electronpositron%e_pawdc
     end if
   end if
   if (optene==0.or.optene==2) electronpositron%e0=etotal
   if (optene==1.or.optene==3) electronpositron%e0=etotal-energies%edc_electronpositron
   etotal=electronpositron%e0+energies%e0_electronpositron+energies%e_electronpositron
 end if

!Compute new charge density based on incoming wf
!Keep rhor and rhog intact for later use e.g. in stress. (? MT 08-12-2008: is that true now ?)
!=== Norm-conserving psps: simply compute rho from WFs
 !paw_dmft%use_dmft=0 ! dmft not used here
 !paw_dmft%use_sc_dmft=0 ! dmft not used here
 if (psps%usepaw==0) then
   tim_mkrho=3
   call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&   npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl_den,wfs)
   if(dtset%usekden==1)then
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&     npwarr,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,wvl_den,wfs,option=1)
   end if
 else
!  === PAW case: symmetrize rhoij and add compensation charge density
   tim_mkrho=3;option=1;choice=1
   call pawrhoij_symrhoij(pawrhoij,pawrhoij_unsym,choice,gprimd,indsym,0,dtset%natom,dtset%nsym,&
&   dtset%ntypat,option,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call pawrhoij_free_unpacked(pawrhoij_unsym)
   if (paral_atom) then
     call pawrhoij_free(pawrhoij_unsym)
     ABI_DATATYPE_DEALLOCATE(pawrhoij_unsym)
   end if
   ider=0;izero=0;cplex=1;ipert=0;idir=0;qpt(:)=zero
   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,&
&   my_natom,dtset%natom,nfftf,ngfftf,&
&   0,dtset%nspden,dtset%ntypat,pawang,pawfgrtab,rhodum,nhat,pawrhoij,pawrhoij,&
&   pawtab,qpt,rprimd,ucvol_local,dtset%usewvl,xred,&
&   comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&   distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)

   ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
   ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
   rhowfr(:,:)=zero
   call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&   npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol_local,wvl_den,wfs)
   call transgrid(1,mpi_enreg,dtset%nspden,+1,1,0,dtset%paral_kgb,pawfgr,rhowfg,rhodum,rhowfr,rhor)
   rhor(:,:)=rhor(:,:)+nhat(:,:)
   call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
   if(dtset%usekden==1)then
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&               rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,wvl_den,wfs,option=1)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,taug,rhowfr,taur)
   end if
   ABI_DEALLOCATE(rhowfr)
   ABI_DEALLOCATE(rhowfg)
 end if

 MSG_COMMENT('New density rho(r) made from input wfs')

 call timab(59,2,tsec)

 ABI_DEALLOCATE(vlocal)
 if (with_vxctau) then
   ABI_DEALLOCATE(vxctaulocal)
 end if

 if (psps%usepaw==1) then
   call pawcprj_free(cwaveprj)
   ABI_DATATYPE_DEALLOCATE(cwaveprj)
   if (mpi_enreg%paral_spinor==1) then
     call pawcprj_free(cwaveprj_gat)
     ABI_DATATYPE_DEALLOCATE(cwaveprj_gat)
   else
     nullify(cwaveprj_gat)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine energy
!!***

!!****f* ABINIT/mkresi
!! NAME
!! mkresi
!!
!! FUNCTION
!! Make residuals from knowledge of wf in G space and application of Hamiltonian.
!!
!! INPUTS
!!  cg(2,mcg)=<G|Cnk>=Fourier coefficients of wavefunction
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=index of k-point
!!  isppol=index of spin
!!  mcg=second dimension of the cg array
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands involved in subspace matrix.
!!  npw=number of planewaves in basis sphere at this k point.
!!  prtvol=control print volume and debugging output
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  eig_k(nband)$= \langle C_n \mid H \mid C_n \rangle $ for each band.
!!  resid_k(nband)=residual for each band
!!   $= \langle C_n \mid H H \mid C_n \rangle- \langle C_n \mid H \mid C_n \rangle^2 $.
!!
!! PARENTS
!!      energy
!!
!! CHILDREN
!!      dotprod_g,getghc,prep_getghc,sqnorm_g,timab
!!
!! SOURCE

subroutine mkresi(cg,eig_k,gs_hamk,icg,ikpt,isppol,mcg,mpi_enreg,nband,prtvol,resid_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,ikpt,isppol,mcg,nband,prtvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: eig_k(nband),resid_k(nband)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_getghc=3
 integer :: blocksize,cpopt,iband,iband_last,iblock,iblocksize,ipw,ipw_shift
 integer :: my_nspinor,nblockbd,npw_k
 real(dp) :: doti,dotr
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable,target :: cwavef(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: cwavef_ptr(:,:),ghc_ptr(:,:),gsc_ptr(:,:)
 type(pawcprj_type) :: cwaveprj(0,0)

! *************************************************************************

!Keep track of total time spent in mkresi
 call timab(13,1,tsec)

!Parallelism over FFT and/or bands: define sizes and tabs
 my_nspinor=max(1,gs_hamk%nspinor/mpi_enreg%nproc_spinor)
 if (mpi_enreg%paral_kgb==1) then
   nblockbd=nband/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
 else
   nblockbd=nband/mpi_enreg%nproc_fft
   if (nband/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
 end if
 blocksize=nband/nblockbd

 npw_k=gs_hamk%npw_k
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(ghc,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(gvnlxc,(2,npw_k*my_nspinor))
 if (gs_hamk%usepaw==1)  then
   ABI_ALLOCATE(gsc,(2,npw_k*my_nspinor))
 else
   ABI_ALLOCATE(gsc,(0,0))
 end if

!Loop over (blocks of) bands
 do iblock=1,nblockbd
   iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband)
   if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,mpi_enreg%me_kpt)) cycle

!  Load |Cn>
   ipw_shift=(iblock-1)*npw_k*my_nspinor*blocksize+icg
!$OMP PARALLEL DO
   do ipw=1,npw_k*my_nspinor*blocksize
     cwavef(1,ipw)=cg(1,ipw+ipw_shift)
     cwavef(2,ipw)=cg(2,ipw+ipw_shift)
   end do

!  Compute H|Cn>
   cpopt=-1
   if (mpi_enreg%paral_kgb==0) then
     call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlxc,zero,mpi_enreg,1,&
&     prtvol,gs_hamk%usepaw,tim_getghc,0)
   else
     call prep_getghc(cwavef,gs_hamk,gvnlxc,ghc,gsc,zero,nband,mpi_enreg,&
&     prtvol,gs_hamk%usepaw,cpopt,cwaveprj,&
&     already_transposed=.false.)
   end if

!  Compute the residual, <Cn|(H-<Cn|H|Cn>)**2|Cn>:
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize
     ipw_shift=(iblocksize-1)*npw_k*my_nspinor
     cwavef_ptr => cwavef(:,1+ipw_shift:npw_k*my_nspinor+ipw_shift)
     ghc_ptr    => ghc   (:,1+ipw_shift:npw_k*my_nspinor+ipw_shift)

!    First get eigenvalue <Cn|H|Cn>:
     call dotprod_g(dotr,doti,gs_hamk%istwf_k,npw_k*my_nspinor,1,cwavef_ptr,ghc_ptr,&
&     mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     eig_k(iband)=dotr

!    Next need <G|(H-S<Cn|H|Cn>)|Cn> (in ghc):
     if (gs_hamk%usepaw==0) then
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(cwavef_ptr,ghc_ptr,eig_k,iband,npw_k,my_nspinor)
       do ipw=1,npw_k*my_nspinor
         ghc_ptr(1,ipw)=ghc_ptr(1,ipw)-eig_k(iband)*cwavef_ptr(1,ipw)
         ghc_ptr(2,ipw)=ghc_ptr(2,ipw)-eig_k(iband)*cwavef_ptr(2,ipw)
       end do
     else
       gsc_ptr => gsc(:,1+ipw_shift:npw_k*my_nspinor+ipw_shift)
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(gsc_ptr,ghc_ptr,eig_k,iband,npw_k,my_nspinor)
       do ipw=1,npw_k*my_nspinor
         ghc_ptr(1,ipw)=ghc_ptr(1,ipw)-eig_k(iband)*gsc_ptr(1,ipw)
         ghc_ptr(2,ipw)=ghc_ptr(2,ipw)-eig_k(iband)*gsc_ptr(2,ipw)
       end do
     end if

!    Then simply square the result:
     call sqnorm_g(dotr,gs_hamk%istwf_k,npw_k*my_nspinor,ghc_ptr,&
&     mpi_enreg%me_g0,mpi_enreg%comm_fft)
     resid_k(iband)=dotr

   end do ! iblocksize

 end do ! iblock

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlxc)
 ABI_DEALLOCATE(gsc)

 call timab(13,2,tsec)

end subroutine mkresi
!!***

end module m_dft_energy
!!***

!!****m* ABINIT/m_vtorho
!! NAME
!!  m_vtorho
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MF, AR, MM, MT, FJ, MB, MT, TR)
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

module m_vtorho

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_xmpi
 use m_ab7_mixing
 use m_errors
 use m_wffile
 use m_efield
 use m_cgtools
 use m_gemm_nonlop
 use m_hdr
 use m_dtset
 use m_dtfil

 use defs_datatypes,       only : pseudopotential_type
 use defs_abitypes,        only : MPI_type
 use m_time,               only : timab
 use m_geometry,           only : xred2xcart
 use m_occ,                only : newocc
 use m_pawang,             only : pawang_type
 use m_pawtab,             only : pawtab_type
 use m_paw_ij,             only : paw_ij_type
 use m_pawfgrtab,          only : pawfgrtab_type
 use m_pawrhoij,           only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_io, pawrhoij_inquire_dim
 use m_pawcprj,            only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_getdim
 use m_pawfgr,             only : pawfgr_type
 use m_energies,           only : energies_type
 use m_hamiltonian,        only : init_hamiltonian, gs_hamiltonian_type
 use m_bandfft_kpt,        only : bandfft_kpt, bandfft_kpt_type, bandfft_kpt_set_ikpt, &
                                  bandfft_kpt_savetabs, bandfft_kpt_restoretabs, prep_bandfft_tabs
 use m_electronpositron,   only : electronpositron_type,electronpositron_calctype
 use m_paw_dmft,           only : paw_dmft_type,init_dmft,destroy_dmft,print_dmft,saveocc_dmft
 use m_paw_correlations,   only : setnoccmmp
 use m_paw_occupancies,   only : pawmkrhoij
 use m_paw_mkrho,          only : pawmkrho
 use m_crystal,            only : crystal_init, crystal_t
 use m_oper,               only : oper_type,init_oper,destroy_oper
 use m_io_tools,           only : flush_unit
 use m_abi2big,            only : wvl_occ_abi2big, wvl_rho_abi2big, wvl_occopt_abi2big, wvl_eigen_abi2big
 use m_fock,               only : fock_type, fock_ACE_type, fock_updateikpt, fock_calc_ene
 use m_invovl,             only : make_invovl
 use m_tddft,              only : tddft
 use m_kg,                 only : mkkin, mkkpg
 use m_suscep_stat,        only : suscep_stat
 use m_fft,                only : fftpac
 use m_spacepar,           only : symrhg
 use m_vtowfk,             only : vtowfk
 use m_mkrho,              only : mkrho, prtrhomxmn
 use m_mkffnl,             only : mkffnl
 use m_mpinfo,             only : proc_distrb_cycle
 use m_common,             only : prteigrs
 use m_dmft,               only : dmft_solve
 use m_datafordmft,        only : datafordmft
 use m_fourier_interpol,   only : transgrid
 use m_cgprj,              only : ctocprj
 use m_wvl_rho,            only : wvl_mkrho
 use m_wvl_psi,            only : wvl_hpsitopsi, wvl_psitohpsi, wvl_nl_gradient

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

#if defined HAVE_BIGDFT
 use BigDFT_API,           only : last_orthon, evaltoocc, write_energies, eigensystem_info
#endif

 implicit none

 private
!!***

 public :: vtorho
!!***

contains
!!***

!!****f* ABINIT/vtorho
!! NAME
!! vtorho
!!
!! FUNCTION
!! This routine compute the new density from a fixed potential (vtrial)
!! but might also simply compute eigenvectors and eigenvalues.
!! The main part of it is a wf update over all k points.
!!
!! INPUTS
!!  afford=used to dimension susmat
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cpus= cpu time limit in seconds
!!  dbl_nnsclo=if 1, will double the value of dtset%nnsclo
!!  dielop= if positive, the dielectric matrix must be computed.
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dmatpawu= fixed occupation matrix of correlated orbitals (DFT+U or DMFT only)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points treated by this node.
!!   | mpw=maximum dimensioned size of npw
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!   | typat= array of types of the natoms
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  etotal=total energy (Ha) - only needed for tddft
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for diel matrix
!!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfftf,nkxc)=exchange-correlation kernel, needed only if nkxc/=0 .
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mpi_enreg=informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!                see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  npwdiel=size of the susmat array.
!!  ntypat=number of types of atoms in unit cell.
!!  optforces=option for the computation of forces (0: no force;1: forces)
!!  optres=0: the new value of the density is computed in place of the input value
!!         1: only the density residual is computed ; the input density is kept
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!                                    nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  phnonsdiel(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases,
!!   for diel matr
!!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1ddiel(2,3*(2*mgfftdiel+1)*natom*usepaw)=one-dimensional structure factor information
!!                                             for the dielectric matrix
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive vectors
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume in bohr**3.
!!  usecprj=1 if cprj datastructure is stored in memory
!!  wffnew,unit numbers for wf disk files.
!!  with_vectornd = 1 if vectornd allocated
!!  vectornd(with_vectornd*nfftf,3)=nuclear dipole moment vector potential
!!  vtrial(nfftf,nspden)=INPUT potential Vtrial(r).
!!  [vxctau(nfft,nspden,4*usekden)]=(only for meta-GGA): derivative of XC energy density
!!    with respect to kinetic energy density (depsxcdtau). The arrays vxctau contains also
!!    the gradient of vxctau (gvxctau) in vxctau(:,:,2:4)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!  ylmdiel(npwdiel,lmax_diel**2)= real spherical harmonics for each G and k point
!!                                 for the dielectric matrix
!!
!! OUTPUT
!!  compch_fft=-PAW only- compensation charge inside spheres computed over fine fft grid
!!  dphase(3) : dphase(idir) = accumulated change in the string-averaged
!!     Zak phase along the idir-th direction caused by the update of all
!!     the occupied Bloch states at all the k-points (only if finite electric field)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  susmat(2,npwdiel*afford,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  === if optforces>0 ===
!!    grnl(3*natom)=stores grads of nonlocal energy wrt length scales
!!  ==== if optres==1
!!    nres2=square of the norm of the residual
!!    nvresid(nfftf,nspden)=density residual
!!  ==== if psps%usepaw==1
!!    cprj(natom,mcprj*usecprj)= wave functions projected with non-local projectors:
!!                               cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!    nhat(nfftf,nspden*psps%usepaw)=compensation charge density on rectangular grid in real space
!!
!! SIDE EFFECTS
!!  cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!   At output contains updated wavefunctions coefficients;
!!    if nkpt>1, these are kept in a disk file.
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_kinetic=kinetic energy part of total energy
!!   | e_nlpsp_vfock=nonlocal psp + potential Fock ACE part of total energy
!!   | e_fermie=fermi energy (Hartree)
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k.
!!      (input if insulator - occopt<3 - ; output if metallic)
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfftf)=Fourier transform of total electron density
!!  rhor(nfftf,nspden)=total electron density (el/bohr**3)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  tauresid(nfftf,nspden*dtset%usekden)=array for kinetic energy density residual
!!  wvl <type(wvl_data)>=wavelets structures in case of wavelets basis.
!!  ==== if (usepaw==1) ====
!!    cprj(natom,mcprj*usecprj)= wave functions projected with non-local projectors:
!!                               cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!!  The total electronic density (rhor,rhog) is divided into two terms:
!!   - The density related to WFs =Sum[Psi**2]
!!   - The compensation density (nhat) - only in PAW
!!
!!  The parallelisation needed for the electric field should be
!!  made an independent subroutine, so that this routine could be put
!!  back in the 95_drive directory.
!!
!! SOURCE

subroutine vtorho(afford,atindx,atindx1,cg,compch_fft,cprj,cpus,dbl_nnsclo,&
&           dielop,dielstrt,dmatpawu,dphase,dtefield,dtfil,dtset,&
&           eigen,electronpositron,energies,etotal,gbound_diel,&
&           gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&           istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mcprj,mgfftdiel,mpi_enreg,&
&           my_natom,natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&           npwarr,npwdiel,nres2,ntypat,nvresid,occ,optforces,&
&           optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
&           phnons,phnonsdiel,ph1d,ph1ddiel,psps,fock,&
&           pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,&
&           rmet,rprimd,susmat,symrec,taug,taur,tauresid,&
&           ucvol,usecprj,wffnew,with_vectornd,vectornd,vtrial,vxctau,wvl,xred,&
&           ylm,ylmgr,ylmdiel)

!Arguments -------------------------------
 integer, intent(in) :: afford,dbl_nnsclo,dielop,dielstrt,istep,istep_mix,lmax_diel,mcg,mcprj,mgfftdiel
 integer, intent(in) :: my_natom,natom,nfftf,nfftdiel,nkxc,npwdiel
 integer, intent(in) :: ntypat,optforces,optres,pwind_alloc,usecprj,with_vectornd
 real(dp), intent(in) :: cpus,etotal,gsqcut,ucvol
 real(dp), intent(out) :: compch_fft,nres2,residm
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(inout) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type), intent(inout) :: energies
 type(hdr_type), intent(inout) :: hdr
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 type(pawfgr_type), intent(in) :: pawfgr
 type(pseudopotential_type), intent(in) :: psps
 type(fock_type),pointer, intent(inout) :: fock
 type(wffile_type), intent(inout) :: wffnew
 type(wvl_data), intent(inout) :: wvl
 integer, intent(in) :: atindx(natom),atindx1(natom),gbound_diel(2*mgfftdiel+8,2)
 integer, intent(in) :: indsym(4,dtset%nsym,natom)
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),kg_diel(3,npwdiel),nattyp(ntypat),ngfftdiel(18),npwarr(dtset%nkpt)
 integer, intent(in) :: pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: dmatpawu(:,:,:,:),gmet(3,3),gprimd(3,3),ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp), intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*psps%usepaw)
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc),rmet(3,3),rprimd(3,3)
 real(dp), intent(inout) :: vectornd(with_vectornd*nfftf,3),vtrial(nfftf,dtset%nspden)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 real(dp), intent(out) :: dphase(3),grnl(3*natom)
 real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(out) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(out) :: nvresid(nfftf,dtset%nspden)
 real(dp), intent(out) :: susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(inout) :: kxc(nfftf,nkxc),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden),taur(nfftf,dtset%nspden*dtset%usekden)
 real(dp), intent(inout) :: tauresid(nfftf,dtset%nspden*dtset%usekden)
 real(dp), intent(inout),optional :: vxctau(nfftf,dtset%nspden,4*dtset%usekden)
 type(pawcprj_type),pointer,intent(inout) :: cprj(:,:)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=111,tim_mkrho=2
 integer,save :: nwarning=0
 integer :: bdtot_index,counter,cplex,cplex_rhoij,dimffnl,enunit,iband,iband1,ibdkpt
 integer :: ibg,icg,ider,idir,ierr,ifft,ifor,ifor1,ii,ikg,ikpt
 integer :: ikpt_loc,ikpt1,my_ikpt,ikxc,ilm,imagn,index1,iorder_cprj,ipert,iplex
 integer :: iscf,ispden,isppol,istwf_k,mband_cprj,mbdkpsp,mb2dkpsp
 integer :: mcgq,mcprj_local,mcprj_tmp,me_distrb,mkgq,mpi_comm_sphgrid
 integer :: mwarning,my_nspinor,n1,n2,n3,n4,n5,n6,nband_eff
 integer :: nband_k,nband_cprj_k,nbuf,neglect_pawhat,nfftot,nkpg,nkpt1,nnn,nnsclo_now
 integer :: nproc_distrb,npw_k,nspden_rhoij,option,prtvol
 integer :: spaceComm_distrb,usecprj_local,usefock_ACE,usetimerev
 logical :: berryflag,computesusmat,fixed_occ,has_vectornd
 logical :: locc_test,paral_atom,remove_inv,usefock,with_vxctau
 logical :: do_last_ortho,wvlbigdft=.false.
 real(dp) :: dmft_dftocc
 real(dp) :: edmft,ebandlda,ebanddmft,ebandldatot,ekindmft,ekindmft2,ekinlda
 real(dp) :: min_occ,vxcavg_dum,strsxc(6)
 character(len=500) :: message
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: dielar(7),dphase_k(3),kpoint(3),qpt(3),rhodum(1),tsec(2),ylmgr_dum(0,0,0)
 real(dp),allocatable :: EigMin(:,:),buffer1(:),buffer2(:),cgq(:,:)
 real(dp),allocatable :: cgrkxc(:,:),cgrvtrial(:,:),doccde(:)
 real(dp),allocatable :: dphasek(:,:),eig_k(:),ek_k(:),ek_k_nd(:,:,:),eknk(:),eknk_nd(:,:,:,:,:)
 real(dp),allocatable :: enlx_k(:),enlxnk(:),focknk(:),fockfornk(:,:,:),ffnl(:,:,:,:),grnl_k(:,:), xcart(:,:)
 real(dp),allocatable :: grnlnk(:,:),kinpw(:),kpg_k(:,:),occ_k(:),ph3d(:,:,:)
 real(dp),allocatable :: pwnsfacq(:,:),resid_k(:),rhoaug(:,:,:,:)
 real(dp),allocatable :: rhowfg(:,:),rhowfr(:,:),tauwfg(:,:),tauwfr(:,:)
 real(dp),allocatable :: vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:),vlocal_tmp(:,:,:)
 real(dp),allocatable :: vxctaulocal(:,:,:,:,:),ylm_k(:,:),zshift(:)
 complex(dpc),target,allocatable :: nucdipmom_k(:)
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)
 type(pawcprj_type),allocatable,target:: cprj_local(:,:)
 type(oper_type) :: dft_occup
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)

 type(crystal_t) :: cryst_struc
 integer :: idum1(0),idum3(0,0,0)
 real(dp) :: rdum2(0,0),rdum4(0,0,0,0)

!Variables for BigDFT
#if defined HAVE_BIGDFT
 integer :: occopt_bigdft
#endif

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in vtorho
 call timab(980,1,tsec)
 call timab(981,1,tsec)

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)

!Several inits
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 usecprj_local=0;if (psps%usepaw==1) usecprj_local=1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 paral_atom=(my_natom/=natom)
 compch_fft=-1.d5

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = (present(vxctau).and.dtset%usekden/=0)

!Check that fock is present if want to use fock option
 usefock = (dtset%usefock==1 .and. associated(fock))
 usefock_ACE=0
 if (usefock) usefock_ACE=fock%fock_common%use_ACE

!Init MPI
 spaceComm_distrb=mpi_enreg%comm_cell
 if (mpi_enreg%paral_kgb==1) spaceComm_distrb=mpi_enreg%comm_kpt
 if (mpi_enreg%paral_hf ==1) spaceComm_distrb=mpi_enreg%comm_kpt
 nproc_distrb=xmpi_comm_size(spaceComm_distrb)
 me_distrb=xmpi_comm_rank(spaceComm_distrb)
 mpi_comm_sphgrid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) mpi_comm_sphgrid=mpi_enreg%comm_wvl
 if (mpi_enreg%me_img/=0) nwarning=nwarning+1

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   MSG_BUG('wrong values for nfft, nfftf!')
 end if

!Test optforces (to prevent memory overflow)
 if (optforces/=0.and.optforces/=1) then
   write(message,'(a,i0)')' wrong value for optforces = ',optforces
   MSG_BUG(message)
 end if

 iscf=dtset%iscf
 fixed_occ=(dtset%occopt<3.or.electronpositron_calctype(electronpositron)==1)
 if(.not. wvlbigdft) then
   energies%e_eigenvalues = zero
   energies%e_kinetic     = zero
   energies%e_nlpsp_vfock = zero
   if (usefock) then
     energies%e_fock=zero
     energies%e_fockdc=zero
   end if
   grnl(:)=zero
   resid(:) = zero ! JWZ 13 May 2010. resid and eigen need to be fully zeroed each time before use
   eigen(:) = zero
   bdtot_index=0
   ibg=0;icg=0
   mbdkpsp=dtset%mband*dtset%nkpt*dtset%nsppol
   if(paw_dmft%use_dmft==1) mb2dkpsp=2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol
 end if

 if(dtset%usewvl==0) then
   ABI_ALLOCATE(eknk,(mbdkpsp))
   ABI_ALLOCATE(enlxnk,(mbdkpsp))
   ABI_ALLOCATE(eknk_nd,(dtset%nsppol,dtset%nkpt,2,dtset%mband,dtset%mband*paw_dmft%use_dmft))
   ABI_ALLOCATE(EigMin,(2,dtset%mband))
   ABI_ALLOCATE(grnlnk,(3*natom,mbdkpsp*optforces))
   if (usefock) then
     ABI_ALLOCATE(focknk,(mbdkpsp))
     focknk=zero
     if (optforces>0)then
       ABI_ALLOCATE(fockfornk,(3,natom,mbdkpsp))
       fockfornk=zero
     end if
   end if
   eknk(:)=zero;enlxnk(:)=zero
   if (optforces>0) grnlnk(:,:)=zero
   if(paw_dmft%use_dmft==1) eknk_nd=zero
 end if !usewvl==0

!Initialize rhor if needed; store old rhor
 if(iscf>=0 .or. iscf==-3) then
   if (optres==1) then
     nvresid=rhor ; tauresid=taur
   end if
!  NC and plane waves
   if (psps%usepaw==0 .and. dtset%usewvl==0) then
     rhor=zero ; taur=zero
!  PAW
   else if(psps%usepaw==1) then
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     ABI_ALLOCATE(tauwfr,(dtset%nfft,dtset%nspden*dtset%usekden))
     ABI_ALLOCATE(tauwfg,(2,dtset%nfft*dtset%usekden))
     rhowfr(:,:)=zero ; tauwfr(:,:)=zero
   end if
 end if

!Set max number of non-self-consistent loops nnsclo_now for use in vtowfk
 if(iscf<0)then
   ! ===== Non self-consistent =====
   nnsclo_now=dtset%nstep
 else
   ! ===== Self-consistent =====
   if(dtset%nnsclo>0) then
   ! ===== Self-consistent + imposed =====
     nnsclo_now=dtset%nnsclo
   else if (dtset%nnsclo<0) then
   ! ===== Self-consistent + imposed during abs(nnsclo) steps =====
     nnsclo_now=1
     if (istep<=abs(dtset%nnsclo)) nnsclo_now=merge(5,dtset%useria,dtset%useria==0)
   else
   ! ===== Self-consistent + default =====
     nnsclo_now=1
     if (dtset%usewvl==0) then
     ! ----- Plane waves -----
       if (istep<=2.and.iscf/=0) nnsclo_now=2
     else
     ! ----- Wavelets -----
       if (iscf==0) then
         nnsclo_now=0
       else if (istep<=2) then
         nnsclo_now=3
       else if (istep<=4) then
         nnsclo_now=2
       end if
     end if
   end if
   ! ===== Double is required =====
   if(dbl_nnsclo==1)then
     nnsclo_now=nnsclo_now*2
   end if
 end if
 if(dtset%wfoptalg==2)nnsclo_now=40  ! UNDER DEVELOPMENT

 write(message, '(a,i0,a,3(i0,1x))' ) ' vtorho: nnsclo_now = ',nnsclo_now,&
   ', note that nnsclo, dbl_nnsclo, istep= ',dtset%nnsclo,dbl_nnsclo,istep
 call wrtout(std_out,message)

!==== Initialize most of the Hamiltonian ====
!Allocate all arrays and initialize quantities that do not depend on k and spin.
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& paw_ij=paw_ij,ph1d=ph1d,usecprj=usecprj_local,electronpositron=electronpositron,fock=fock,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

!Initializations for PAW (projected wave functions)
 mcprj_local=0 ; mband_cprj=0
 if (psps%usepaw==1) then
   mband_cprj=dtset%mband
   if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
   iorder_cprj=0 ; mcprj_local=mcprj
   if (usecprj==0) then
     mcprj_local=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
     !This is a check but should always be true since scfcv allocated cprj anyway
     if (allocated(cprj_local)) then
       !Was allocated in scfcv so we just destroy and reconstruct it as desired
       call pawcprj_free(cprj_local)
       ABI_DATATYPE_DEALLOCATE(cprj_local)
     end if
     ABI_DATATYPE_ALLOCATE(cprj_local,(dtset%natom,mcprj_local))
     call pawcprj_alloc(cprj_local,0,gs_hamk%dimcprj)
     cprj=> cprj_local
   end if
 end if

 call timab(981,2,tsec)

!===================================================================
! WAVELETS - Branching with a separate VTORHO procedure
!===================================================================

 if (dtset%usewvl == 1) then
#ifndef HAVE_BIGDFT
   BIGDFT_NOTENABLED_ERROR()
#else

!  do_last_ortho in case of diagonalization scheme
   if (     wvlbigdft) do_last_ortho=(dtset%iscf/=0)
   if (.not.wvlbigdft) do_last_ortho=(.true.)

   ABI_ALLOCATE(xcart,(3, dtset%natom))
   call xred2xcart(dtset%natom, rprimd, xcart, xred)

   if(wvlbigdft) then
!    NSCF loop for wvlbigdt:
     call wvl_nscf_loop_bigdft()
   else
!    NSCF loop for WVL: (not wvlbigdft)
     call wvl_nscf_loop()
   end if

!  Eventually orthogonalize WFs now
   if (do_last_ortho) then
     call write_energies(ii,0,wvl%e%energs,0.d0,0.d0,"final")
     call last_orthon(me_distrb, nproc_distrb, ii, wvl%wfs%ks, wvl%e%energs%evsum, .true.)
     if(wvlbigdft) energies%e_xcdc = wvl%e%energs%evxc
!    If occupation numbers are not changed...
     if (fixed_occ .or. (iscf<0 .and. iscf/=-3)) then
       call wvl_comm_eigen()
     end if
!    ... or update occupations:
     if( ( .not.fixed_occ) .and. (iscf>0.or.iscf==-3)) then
       if(wvlbigdft) then
         call wvl_occ_bigdft()
       else
!        Communicate eigenvalues:
         call wvl_comm_eigen()
!        Update occ and Fermi level
         call wvl_occ()
       end if
     end if
!    This might accelerate convergence:
     wvl%wfs%ks%diis%energy_min=one
     wvl%wfs%ks%diis%alpha=two
   end if !do_last_ortho

!  Compute eigenvalues energy
   if(.not. wvlbigdft .and. nnsclo_now>0) then
     call e_eigen(eigen,energies%e_eigenvalues,dtset%mband,dtset%nband,dtset%nkpt,&
&     dtset%nsppol,occ,dtset%wtk)
   else
     energies%e_eigenvalues = energies%e_kinetic + energies%e_localpsp &
&     + energies%e_xcdc  + two*energies%e_hartree +energies%e_nlpsp_vfock
   end if

   if (optforces == 1) then ! not compatible with iscf=0 and wvlbigdftcomp=1 + PAW
     call wvl_nl_gradient(grnl, mpi_enreg, dtset%natom, rprimd, wvl, xcart)
   end if

!  For iscf<0 we do not update the density
   if (dtset%iscf>=0) then !(dtset%iscf>=0 .and. .not. wvlbigdft ) then
     call wvl_mkrho(dtset,irrzon,mpi_enreg,phnons,rhor,wvl%wfs,wvl%den)
   end if
   ABI_DEALLOCATE(xcart)

!  Note in WVL+NC: the rest will be skipped.
!  For PAW: we will compute Rho_ij at the end.
   !if(wvlbigdft) return
#endif
 else

!===================================================================
! PLANE WAVES - Standard VTORHO procedure
!===================================================================

!  Electric fields: set flag to turn on various behaviors
   berryflag = (dtset%berryopt == 4 .or. dtset%berryopt == 14 .or. &
&   dtset%berryopt == 6 .or. dtset%berryopt == 16 .or. &
&   dtset%berryopt == 7 .or. dtset%berryopt == 17)

!  Electric field: allocate dphasek
   nkpt1 = dtset%nkpt
   if ( berryflag ) then
     ABI_ALLOCATE(dphasek,(3,dtset%nkpt*dtset%nsppol))
     dphasek(:,:) = zero
     nkpt1 = dtefield%mkmem_max
   end if

   ABI_ALLOCATE(rhoaug,(n4,n5,n6,gs_hamk%nvloc))
   ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamk%nvloc))
   if(with_vxctau) then
     ABI_ALLOCATE(vxctaulocal,(n4,n5,n6,gs_hamk%nvloc,4))
   end if

   has_vectornd = (with_vectornd .EQ. 1)
   if(has_vectornd) then
      ABI_ALLOCATE(vectornd_pac,(n4,n5,n6,gs_hamk%nvloc,3))
   end if

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol
     call timab(982,1,tsec)

     ikpt_loc = 0
     ikg=0

!    Set up local potential vlocal with proper dimensioning, from vtrial
!    Also take into account the spin.
     if(dtset%nspden/=4)then
       if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
         call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
         if(with_vxctau) then
           do ii=1,4
             call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,&
                         vxctau(:,:,ii),vxctaulocal(:,:,:,:,ii),2)
           end do
         end if
       else
         ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
         call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,&
                        rhodum,rhodum,cgrvtrial,vtrial)
         call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,&
                     cgrvtrial,vlocal,2)
         if(with_vxctau) then
           do ii=1,4
             call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,&
                           rhodum,rhodum,cgrvtrial,vxctau(:,:,ii))
             call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,&
                         cgrvtrial,vxctaulocal(:,:,:,:,ii),2)
           end do
         end if
         ABI_DEALLOCATE(cgrvtrial)
       end if
     else
       ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
       if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
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
     end if ! nspden

     rhoaug(:,:,:,:)=zero

     ! if vectornd is present, set it up for addition to gs_hamk similarly to how it's done for
     ! vtrial. Note that it must be done for the three Cartesian directions. Also, the following
     ! code assumes explicitly and implicitly that nvloc = 1. This should eventually be generalized.
     if(has_vectornd) then
        do idir = 1, 3
           ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
           call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,idir))
           call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
           ABI_DEALLOCATE(cgrvtrial)
        end do
     end if

!    Continue to initialize the Hamiltonian
     call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)
     if (with_vxctau) then
       call gs_hamk%load_spin(isppol,vxctaulocal=vxctaulocal)
     end if
     if (has_vectornd) then
       call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
     end if

     call timab(982,2,tsec)

!    BIG FAT k POINT LOOP
!    MVeithen: I had to modify the structure of this loop in order to implement MPI // of the electric field
!    Note that the loop here differs from the similar one in berryphase_new.F90.
!    here, ikpt_loc numbers the kpts treated by the current processor.
!    in berryphase_new.F90, ikpt_loc ALSO includes info about value of isppol.

     ikpt = 0
     do while (ikpt_loc < nkpt1)

       call timab(997,1,tsec)

       if ( .not.berryflag ) then
         ikpt_loc = ikpt_loc + 1
         ikpt = ikpt_loc
         my_ikpt = mpi_enreg%my_kpttab(ikpt)
       else
         if (ikpt_loc < dtset%mkmem) ikpt = ikpt + 1
         if ((ikpt > dtset%nkpt).and.(ikpt_loc < dtset%mkmem)) exit
         my_ikpt=ikpt
       end if

       dphase_k(:) = zero
       counter=100*ikpt+isppol
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       nband_cprj_k=nband_k/mpi_enreg%nproc_band
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)

       mcgq=1 ; mkgq=1
       if (.not.berryflag) then
         if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then
           eigen(1+bdtot_index : nband_k+bdtot_index) = zero
           resid(1+bdtot_index : nband_k+bdtot_index) = zero
           bdtot_index=bdtot_index+nband_k
           cycle
         end if
       else
         if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) .and.(ikpt_loc <= dtset%mkmem)) then
           eigen(1+bdtot_index : nband_k+bdtot_index) = zero
           resid(1+bdtot_index : nband_k+bdtot_index) = zero
           bdtot_index = bdtot_index + nband_k
           cycle
         end if
         ikpt_loc = ikpt_loc + 1
         mcgq = dtset%mpw*my_nspinor*nband_k*dtefield%nneigh(ikpt)
         ikg = dtefield%kgindex(ikpt)
         mkgq = 6*dtset%mpw
       end if ! berryflag

       call timab(997,2,tsec)

!      In case of MPI // of a finite field calculation
!      build the cgq array that stores the wavefunctions for the
!      neighbours of ikpt, and the pwnsfacq array that stores the
!      corresponding phase factors (in case of tnons)
       ABI_ALLOCATE(cgq,(2,mcgq))
       ABI_ALLOCATE(pwnsfacq,(2,mkgq))
       if ( berryflag ) then
         call cgq_builder(berryflag,cg,cgq,dtefield,dtset,ikpt,ikpt_loc,isppol,mcg,mcgq,&
&         me_distrb,mkgq,mpi_enreg,my_nspinor,nband_k,nproc_distrb,&
&         npwarr,pwnsfac,pwnsfacq,pwind_alloc,spaceComm_distrb)
         if (ikpt_loc > dtset%mkmem) then
           ABI_DEALLOCATE(cgq)
           ABI_DEALLOCATE(pwnsfacq)
           cycle
         end if
       end if !berryopt

       call timab(984,1,tsec)

       if (mpi_enreg%paral_kgb==1) my_bandfft_kpt => bandfft_kpt(my_ikpt)
       call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)
!      my_ikpt = ikpt
!       if (mpi_enreg%paral_kgb==1) then
!        my_ikpt = mpi_enreg%my_kpttab(ikpt)
!         my_bandfft_kpt => bandfft_kpt(my_ikpt)
!         call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)
!       end if

       ABI_ALLOCATE(eig_k,(nband_k))
       ABI_ALLOCATE(ek_k,(nband_k))
       ABI_ALLOCATE(enlx_k,(nband_k))
       ABI_ALLOCATE(ek_k_nd,(2,nband_k,nband_k*paw_dmft%use_dmft))
       ABI_ALLOCATE(occ_k,(nband_k))
       ABI_ALLOCATE(resid_k,(nband_k))
       ABI_ALLOCATE(zshift,(nband_k))
       ABI_ALLOCATE(grnl_k,(3*natom,nband_k*optforces))

       eig_k(:)=zero
       ek_k(:)=zero
       enlx_k(:)=zero
       if(paw_dmft%use_dmft==1) ek_k_nd(:,:,:)=zero
       if (optforces>0) grnl_k(:,:)=zero
       kpoint(:)=dtset%kptns(:,ikpt)
       occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
       resid_k(:)=zero
       zshift(:)=dtset%eshift

       ABI_ALLOCATE(kg_k,(3,npw_k))
       ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       if (psps%useylm==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
         end do
       end if

!      Set up remaining of the Hamiltonian

!      Compute (1/2) (2 Pi)**2 (k+G)**2:
       ABI_ALLOCATE(kinpw,(npw_k))
!       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k)
       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

!      Compute (k+G) vectors (only if useylm=1)
       nkpg=3*optforces*dtset%nloalg(3)
       ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
       if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
         call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
       end if

!      Compute nonlocal form factors ffnl at all (k+G):
       ider=0;idir=0;dimffnl=1
       ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
       if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
         call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&         npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&         psps%usepaw,psps%useylm,ylm_k,ylmgr)
       end if

!      Load k-dependent part in the Hamiltonian datastructure
!       - Compute 3D phase factors
!       - Prepare various tabs in case of band-FFT parallelism
!       - Load k-dependent quantities in the Hamiltonian
       ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))

       if(usefock_ACE/=0) then
         call gs_hamk%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,&
&         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,fockACE_k=fock%fockACE(ikpt,isppol),ph3d_k=ph3d,&
&         compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1),&
&         compute_gbound=(mpi_enreg%paral_kgb/=1))
       else
         call gs_hamk%load_k(kpt_k=dtset%kptns(:,ikpt),istwf_k=istwf_k,npw_k=npw_k,&
&         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,&
&         compute_ph3d=(mpi_enreg%paral_kgb/=1.or.istep<=1),&
&         compute_gbound=(mpi_enreg%paral_kgb/=1))
       end if

!      Load band-FFT tabs (transposed k-dependent arrays)
       if (mpi_enreg%paral_kgb==1) then
         if (istep<=1) then
           call prep_bandfft_tabs(gs_hamk,ikpt,dtset%mkmem,mpi_enreg)
         end if
         call gs_hamk%load_k(npw_fft_k=my_bandfft_kpt%ndatarecv, &
&         gbound_k =my_bandfft_kpt%gbound, &
&         kinpw_k  =my_bandfft_kpt%kinpw_gather, &
&         kg_k     =my_bandfft_kpt%kg_k_gather, &
&         kpg_k    =my_bandfft_kpt%kpg_k_gather, &
         ffnl_k   =my_bandfft_kpt%ffnl_gather, &
         ph3d_k   =my_bandfft_kpt%ph3d_gather)
       end if

!      Build inverse of overlap matrix for chebfi
       if(psps%usepaw == 1 .and. dtset%wfoptalg == 1 .and. istep <= 1) then
         call make_invovl(gs_hamk, dimffnl, ffnl, ph3d, mpi_enreg)
       end if

       ! Setup gemm_nonlop
       if (gemm_nonlop_use_gemm) then
         !set the global variable indicating to gemm_nonlop where to get its data from
         gemm_nonlop_ikpt_this_proc_being_treated = my_ikpt
         if (istep <= 1) then
           !Init the arrays
           call make_gemm_nonlop(my_ikpt,gs_hamk%npw_fft_k,gs_hamk%lmnmax, &
&           gs_hamk%ntypat, gs_hamk%indlmn, gs_hamk%nattyp, gs_hamk%istwf_k, gs_hamk%ucvol, gs_hamk%ffnl_k,&
&           gs_hamk%ph3d_k)
         end if
       end if

#if defined HAVE_GPU_CUDA
       if (dtset%use_gpu_cuda==1) then
         call gpu_update_ffnl_ph3d(ph3d,size(ph3d),ffnl,size(ffnl))
       end if
#endif

       call timab(984,2,tsec)

!      Update the value of ikpt,isppol in fock_exchange and allocate the memory space to perform HF calculation.
       if (usefock) then
         call fock_updateikpt(fock%fock_common,ikpt,isppol)
       end if
       if ((psps%usepaw==1).and.(usefock)) then
         if ((fock%fock_common%optfor).and.(usefock_ACE==0)) then
           fock%fock_common%forces_ikpt=zero
         end if
       end if

!      Compute the eigenvalues, wavefunction, residuals,
!      contributions to kinetic energy, nonlocal energy, forces,
!      and update of rhor to this k-point and this spin polarization.
       call vtowfk(cg,cgq,cprj,cpus,dphase_k,dtefield,dtfil,&
&       dtset,eig_k,ek_k,ek_k_nd,enlx_k,fixed_occ,grnl_k,gs_hamk,&
&       ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,mband_cprj,mcg,mcgq,mcprj_local,mkgq,&
&       mpi_enreg,dtset%mpw,natom,nband_k,dtset%nkpt,nnsclo_now,npw_k,npwarr,&
&       occ_k,optforces,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&       rhoaug,paw_dmft,dtset%wtk(ikpt),zshift)

       call timab(985,1,tsec)

#if defined HAVE_GPU_CUDA
       if(dtset%use_gpu_cuda==1) call gpu_finalize_ffnl_ph3d()
#endif
       ABI_DEALLOCATE(ffnl)
       ABI_DEALLOCATE(kinpw)
       ABI_DEALLOCATE(kg_k)
       ABI_DEALLOCATE(kpg_k)
       if(allocated(nucdipmom_k)) then
         ABI_DEALLOCATE(nucdipmom_k)
       end if
       ABI_DEALLOCATE(ylm_k)
       ABI_DEALLOCATE(ph3d)
       ABI_DEALLOCATE(cgq)
       ABI_DEALLOCATE(pwnsfacq)

!      electric field
       if (berryflag) then
         dphasek(:,ikpt + (isppol - 1)*dtset%nkpt) = dphase_k(:)

!        The overlap matrices for all first neighbours of ikpt are no more up to date
         do idir = 1, 3
           do ifor = 1, 2
             ikpt1 = dtefield%ikpt_dk(dtefield%i2fbz(ikpt),ifor,idir)
             ikpt1 = dtefield%indkk_f2ibz(ikpt1,1)
             ifor1 = -1*ifor + 3   ! ifor = 1 -> ifor1 = 2 & ifor = 2 -> ifor1 = 1
             dtefield%sflag(:,ikpt1+(isppol-1)*dtset%nkpt,ifor1,idir) = 0
           end do
         end do
       end if  ! berryflag

!      Save eigenvalues (hartree), residuals (hartree**2)
       eigen(1+bdtot_index : nband_k+bdtot_index) = eig_k(:)
       eknk (1+bdtot_index : nband_k+bdtot_index) = ek_k (:)
       if(usefock) then
         focknk (1+bdtot_index : nband_k+bdtot_index) = fock%fock_common%eigen_ikpt (:)
         if (optforces>0) fockfornk(:,:,1+bdtot_index : nband_k+bdtot_index) = fock%fock_common%forces_ikpt(:,:,:)
       end if
       if(paw_dmft%use_dmft==1) eknk_nd(isppol,ikpt,:,:,:) = ek_k_nd(:,:,:)
       resid(1+bdtot_index : nband_k+bdtot_index) = resid_k(:)
       if (optforces>0) grnlnk(:,1+bdtot_index : nband_k+bdtot_index) = grnl_k(:,:)
       enlxnk(1+bdtot_index : nband_k+bdtot_index) = enlx_k(:)

       if(iscf>0 .or. iscf==-3)then
!        Accumulate sum over k points for band, nonlocal and kinetic energies,
!        also accumulate gradients of Enonlocal:
         do iband=1,nband_k
           if (abs(occ_k(iband))>tol8) then
             energies%e_kinetic     = energies%e_kinetic     + dtset%wtk(ikpt)*occ_k(iband)*ek_k(iband)
             energies%e_eigenvalues = energies%e_eigenvalues + dtset%wtk(ikpt)*occ_k(iband)*eig_k(iband)
             energies%e_nlpsp_vfock = energies%e_nlpsp_vfock + dtset%wtk(ikpt)*occ_k(iband)*enlx_k(iband)
             if (optforces>0) grnl(:)=grnl(:)+dtset%wtk(ikpt)*occ_k(iband)*grnl_k(:,iband)
             if (usefock) then
               energies%e_fock=energies%e_fock + half*fock%fock_common%eigen_ikpt(iband)*occ_k(iband)*dtset%wtk(ikpt)
               if (usefock_ACE==0)then
                 energies%e_fock0=energies%e_fock
               endif
             endif
           end if
         end do

!        Calculate Fock contribution to the total energy if required
         if ((psps%usepaw==1).and.(usefock)) then
           if ((fock%fock_common%optfor).and.(usefock_ACE==0)) then
             !WARNING : this routine actually does NOT compute the Fock contrib to total energy, but modifies the force ONLY.
             call fock_calc_ene(dtset,fock%fock_common,energies%e_exactX,ikpt,nband_k,occ_k)
           end if
         end if
       end if

       call timab(985,2,tsec)

       ABI_DEALLOCATE(eig_k)
       ABI_DEALLOCATE(ek_k)
       ABI_DEALLOCATE(ek_k_nd)
       ABI_DEALLOCATE(grnl_k)
       ABI_DEALLOCATE(occ_k)
       ABI_DEALLOCATE(resid_k)
       ABI_DEALLOCATE(zshift)
       ABI_DEALLOCATE(enlx_k)

!      Keep track of total number of bands (all k points so far, even for k points not treated by me)
       bdtot_index=bdtot_index+nband_k

!      Also shift array memory if dtset%mkmem/=0
       if (dtset%mkmem/=0) then
         ibg=ibg+my_nspinor*nband_cprj_k
         icg=icg+npw_k*my_nspinor*nband_k
         ikg=ikg+npw_k
       end if

     end do ! End big k point loop

     call timab(986,1,tsec)

     if (fixed_occ .and. mpi_enreg%paral_kgb==1) then
       call xmpi_sum(rhoaug,mpi_enreg%comm_bandspinorfft,ierr) !Sum the contributions over bands/FFT/spinors
     end if

!    Transfer density on augmented fft grid to normal fft grid in real space
!    Also take into account the spin.
     if(iscf>0.or.iscf==-3)then
       if (psps%usepaw==0) then
         call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug(:,:,:,1),1)
         if(dtset%nspden==4)then
           do imagn=2,4
             call fftpac(imagn,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug(:,:,:,imagn),1)
           end do
         end if
       else
         call fftpac(isppol,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhowfr,rhoaug(:,:,:,1),1)
         if(dtset%nspden==4)then
           do imagn=2,4
             call fftpac(imagn,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhowfr,rhoaug(:,:,:,imagn),1)
           end do
         end if
       end if
     end if

     call timab(986,2,tsec)

   end do ! End loop over spins

   call timab(988,1,tsec)
   if (usefock) then
     if(fock%fock_common%optfor) then
       call xmpi_sum(fock%fock_common%forces,mpi_enreg%comm_kpt,ierr)
     end if
   end if
!  Electric field: compute string-averaged change in Zak phase
!  along each direction, store it in dphase(idir)
!  ji: it is not convenient to do this anymore. Remove. Set dphase(idir)=0.0_dp.
!  eventually, dphase(idir) will have to go...
   if (berryflag) then
     dphase(:) = zero
!    In case of MPI // of a finite field calculation, send dphasek to all cpus
     call xmpi_sum(dphasek,spaceComm_distrb,ierr)
     ABI_DEALLOCATE(dphasek)
   end if ! berryflag
   ABI_DEALLOCATE(rhoaug)
   ABI_DEALLOCATE(vlocal)
   if(with_vxctau) then
     ABI_DEALLOCATE(vxctaulocal)
   end if
   if(has_vectornd) then
      ABI_DEALLOCATE(vectornd_pac)
   end if

   call timab(988,2,tsec)

   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
   doccde(:)=zero

!  Treat now varying occupation numbers, in the self-consistent case
   if((.not.fixed_occ) .and. (iscf>0.or.iscf==-3)) then

!    Parallel case
     if (mpi_enreg%nproc_kpt>1) then

       call timab(989,1,tsec)

!      If needed, exchange the values of eigen,resid,eknk,enlxnk,grnlnk
       ABI_ALLOCATE(buffer1,((4+3*natom*optforces+dtset%usefock+3*natom*dtset%usefock*optforces)*mbdkpsp))
       if(paw_dmft%use_dmft==1) then
         ABI_ALLOCATE(buffer2,(mb2dkpsp*paw_dmft%use_dmft))
       end if
!      Pack eigen,resid,eknk,enlxnk,grnlnk in buffer1
       buffer1(1          :  mbdkpsp)=eigen(:)
       buffer1(1+  mbdkpsp:2*mbdkpsp)=resid(:)
       buffer1(1+2*mbdkpsp:3*mbdkpsp)=eknk(:)
       buffer1(1+3*mbdkpsp:4*mbdkpsp)=enlxnk(:)
       index1=4*mbdkpsp
       if (optforces>0) then
         buffer1(index1+1:index1+3*natom*mbdkpsp)=reshape(grnlnk,(/(3*natom)*mbdkpsp/) )
         index1=index1+3*natom*mbdkpsp
       end if
       if (usefock) then
         buffer1(1+index1:index1+mbdkpsp)=focknk(:)
         if (optforces>0) then
           index1=index1+mbdkpsp
           buffer1(index1+1:index1+3*natom*mbdkpsp)=reshape(fockfornk,(/(3*natom)*mbdkpsp/) )
         end if
       end if
       if(paw_dmft%use_dmft==1) then
         nnn=0
         do ikpt=1,dtset%nkpt
           do isppol=1,dtset%nsppol
             do iband=1,dtset%mband
               do iband1=1,dtset%mband
                 do iplex=1,2
                   nnn=nnn+1
                   buffer2(nnn)=eknk_nd(isppol,ikpt,iplex,iband,iband1)
                 end do
               end do
             end do
           end do
         end do
         if(nnn.ne.mb2dkpsp)  then
           write(message,*)' BUG in vtorho2, buffer2',nnn,mb2dkpsp
           MSG_BUG(message)
         end if
       end if
!      Build sum of everything
       call timab(48,1,tsec)
       call xmpi_sum(buffer1,mpi_enreg%comm_kpt,ierr)
!      if(mpi_enreg%paral_kgb/=1.and.paw_dmft%use_dmft==1) then
       if(paw_dmft%use_dmft==1) then
         call xmpi_sum(buffer2,mpi_enreg%comm_kpt,ierr)
       end if
       call timab(48,2,tsec)

!      Unpack eigen,resid,eknk,enlxnk,grnlnk in buffer1
       eigen(:) =buffer1(1          :  mbdkpsp)
       resid(:) =buffer1(1+  mbdkpsp:2*mbdkpsp)
       eknk(:)  =buffer1(1+2*mbdkpsp:3*mbdkpsp)
       enlxnk(:) =buffer1(1+3*mbdkpsp:4*mbdkpsp)
       index1=4*mbdkpsp
       if (optforces>0) then
         grnlnk(:,:)=reshape(buffer1(index1+1:index1+3*natom*mbdkpsp),(/3*natom,mbdkpsp/) )
       end if
       if (usefock) then
         focknk(:)=buffer1(1+index1:index1+mbdkpsp)
         if (optforces>0) then
           index1=index1+mbdkpsp
           fockfornk(:,:,:)=reshape(buffer1(index1+1:index1+3*natom*mbdkpsp),(/3,natom,mbdkpsp/) )
         end if
       end if
       if(paw_dmft%use_dmft==1) then
         nnn=0
         do ikpt=1,dtset%nkpt
           do isppol=1,dtset%nsppol
             do iband=1,dtset%mband
               do iband1=1,dtset%mband
                 do iplex=1,2
                   nnn=nnn+1
                   eknk_nd(isppol,ikpt,iplex,iband,iband1)=buffer2(nnn)
                 end do
               end do
             end do
           end do
         end do
       end if
       if(allocated(buffer2))  then
         ABI_DEALLOCATE(buffer2)
       end if
       ABI_DEALLOCATE(buffer1)
       call timab(989,2,tsec)

     end if ! nproc_kpt>1

!    Compute the new occupation numbers from eigen
     call timab(990,1,tsec)
     call newocc(doccde,eigen,energies%entropy,energies%e_fermie,dtset%spinmagntarget,&
&     dtset%mband,dtset%nband,dtset%nelect,dtset%nkpt,dtset%nspinor,&
&     dtset%nsppol,occ,dtset%occopt,prtvol,dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)
     call timab(990,2,tsec)

!    !=========  DMFT call begin ============================================
     dmft_dftocc=0
     if(paw_dmft%use_dmft==1.and.psps%usepaw==1.and.dtset%nbandkss==0) then
       call timab(991,1,tsec)

       if (dtset%dmftcheck>=0.and.dtset%usedmft>=1.and.(sum(pawtab(:)%upawu)>=tol8.or.  &
&       sum(pawtab(:)%jpawu)>tol8).and.dtset%dmft_entropy==0) energies%entropy=zero

!      ==  0 to a dmft calculation and do not use lda occupations
!      ==  1 to a lda calculation with the dmft loop
       if(dtset%dmftcheck==-1) dmft_dftocc=1

!      ==  initialise occnd
       paw_dmft%occnd=zero

       bdtot_index = 1
       do isppol=1,dtset%nsppol
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
             paw_dmft%occnd(1,iband,iband,ikpt,isppol)=occ(bdtot_index)
             bdtot_index = bdtot_index + 1
           end do
         end do
       end do


       if(dmft_dftocc==0) then
         if(dtset%occopt/=3) then
           MSG_ERROR('occopt should be equal to 3 in dmft')
         end if
!        ==  initialise edmft
         if(paw_dmft%use_dmft>=1) edmft = zero

         !  Compute residm to check the value
         ibdkpt=1
         residm=zero
         do isppol=1,dtset%nsppol
           do ikpt=1,dtset%nkpt
             nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
             nband_eff=max(1,nband_k-dtset%nbdbuf)
             residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_eff-1)))
             ibdkpt=ibdkpt+nband_k
           end do
         end do

         !  Test residm
         if (paw_dmft%use_dmft>0 .and. residm>tol4 .and. dtset%dmftcheck>=0) then
           if(dtset%dmft_entropy>0)  then
             write(message,'(a,e12.3)')&
               ' WARNING: Wavefunctions not converged: DFT+DMFT calculation cannot be carried out safely ',residm
             call wrtout(std_out,message)
           else
             write(message,'(a,e12.3)')&
              ' ERROR: Wavefunctions not converged: DFT+DMFT calculation cannot be carried out safely ',residm
             call wrtout(std_out,message)
             write(message,'(a,i0)')'  Action: increase nline and nnsclo',dtset%nstep
             MSG_ERROR(message)
           end if

         else if (paw_dmft%use_dmft>0 .and. residm>tol10.and. dtset%dmftcheck>=0) then
           write(message,'(3a)')ch10,&
            '  Wavefunctions not converged: DFT+DMFT calculation might not be carried out safely ',ch10
           MSG_WARNING(message)
         end if

!        ==  allocate paw_dmft%psichi and paw_dmft%eigen_dft
         call init_dmft(dmatpawu,dtset,energies%e_fermie,dtfil%fnameabo_app,&
           dtfil%filnam_ds(3),dtset%nspinor,paw_dmft,pawtab,psps,dtset%typat)
         call print_dmft(paw_dmft,dtset%pawprtvol)

!        ==  gather crystal structure date into data "cryst_struc"
         remove_inv=.false.
         if(dtset%nspden==4) remove_inv=.true.
         call crystal_init(dtset%amu_orig(:,1),cryst_struc,dtset%spgroup,natom,dtset%npsp,ntypat, &
          dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
          dtset%nspden==2.and.dtset%nsppol==1,remove_inv,hdr%title,&
          dtset%symrel,dtset%tnons,dtset%symafm)

!        ==  compute psichi
         call xmpi_barrier(spaceComm_distrb)
         call init_oper(paw_dmft,dft_occup)
         call flush_unit(std_out)
         call timab(620,1,tsec)
         call datafordmft(cryst_struc,cprj,gs_hamk%dimcprj,dtset,eigen,energies%e_fermie,&
&         dft_occup,dtset%mband,mband_cprj,dtset%mkmem,mpi_enreg,&
&         dtset%nkpt,my_nspinor,dtset%nsppol,occ,&
&         paw_dmft,paw_ij,pawang,pawtab,psps,usecprj_local,dtfil%unpaw)
         call timab(620,2,tsec)
         call flush_unit(std_out)


!        ==  solve dmft loop
         call xmpi_barrier(spaceComm_distrb)

         call dmft_solve(cryst_struc,istep,dft_occup,paw_dmft,pawang,pawtab,dtset%pawprtvol)
         edmft=paw_dmft%edmft
         energies%e_paw=energies%e_paw+edmft
         energies%e_pawdc=energies%e_pawdc+edmft
         call flush_unit(std_out)
!        paw_dmft%occnd(:,:,:,:,:)=0.5_dp

!        For compatibility with old test, do not use for calculation
         if(dtset%dmft_occnd_imag==0) paw_dmft%occnd(2,:,:,:,:)=zero

!        call print_dmft(paw_dmft,dtset%pawprtvol)
!         if(dtset%paral_kgb==1) then
!           write(message,'(5a)')ch10,&
!&           ' Parallelization over bands is not yet compatible with self-consistency in DMFT ',ch10,&
!&           ' Calculation of density does not taken into account non diagonal occupations',ch10
!           call wrtout(std_out,message)
!           call wrtout(ab_out,message)
!!          MSG_ERROR(message)
!           if(dtset%nstep>1) then
!             write(message,'(a,i0)')'  Action: use nstep=1 instead of nstep=',dtset%nstep
!             MSG_ERROR(message)
!           end if
!           residm=zero
!         end if
!        if(dtset%nspinor==2) then
!          call flush_unit(ab_out)
!          write(message,'(3a)')&
!          &         ' Self consistent DFT+DMFT with nspinor==2 is not possible yet ',ch10,&
!          &         ' Calculation are restricted to nstep =1'
!          !         MSG_ERROR(message)
!          if(dtset%nstep>1) then
!          write(message,'(a,i0)')' Action: use nstep=1 instead of nstep=',dtset%nstep
!          !           MSG_ERROR(message)
!          endif
!        end if

         if(me_distrb==0) then
           call saveocc_dmft(paw_dmft)
         end if
         call destroy_dmft(paw_dmft)

!        ==  destroy crystal_t cryst_struc
         call cryst_struc%free()
         call destroy_oper(dft_occup)
       end if ! dmft_dftocc
       call timab(991,2,tsec)

     end if ! usedmft

     if(dtset%nbandkss/=0) then
       write(message,'(a,i3,2a,i3,4a)') &
        " dtset%nbandkss = ",dtset%nbandkss,ch10,&
        " and dtset%usedmft = ",dtset%usedmft,ch10,&
        " a DFT loop is carried out without DMFT.",ch10,&
        " Only psichi's will be written at convergence of the DFT loop."
       call wrtout(std_out,message)
     end if
!    !=========  DMFT call end   ============================================

     call timab(992,1,tsec)

!    Compute eeig, ek,enl and grnl from the new occ, and the shared eknk,enlxnk,grnlnk
     energies%e_eigenvalues = zero
     energies%e_kinetic     = zero
     energies%e_nlpsp_vfock = zero
     if (usefock) then
       energies%e_fock     = zero
       if (optforces>0) fock%fock_common%forces=zero
     end if
     if (optforces>0) grnl(:)=zero
     if(paw_dmft%use_dmft>=1) then
       ebandlda               = zero
       ebanddmft              = zero
       ebandldatot            = zero
       ekindmft               = zero
       ekindmft2              = zero
       ekinlda                = zero
     end if

!    Compute new energy terms due to non diagonal occupations and DMFT.
     bdtot_index=1
     do isppol=1,dtset%nsppol
       do ikpt=1,dtset%nkpt
         nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
         do iband=1,nband_k

           locc_test = abs(occ(bdtot_index))>tol8
!          dmft
           if(paw_dmft%use_dmft>=1.and.dtset%nbandkss==0) then
             if(paw_dmft%band_in(iband)) then
               if( paw_dmft%use_dmft == 1 .and. dmft_dftocc == 1 ) then ! test of the code
                 paw_dmft%occnd(1,iband,iband,ikpt,isppol)=occ(bdtot_index)
               end if
               locc_test = abs(paw_dmft%occnd(1,iband,iband,ikpt,isppol))+&
&               abs(paw_dmft%occnd(2,iband,iband,ikpt,isppol))>tol8
             end if
           end if

           if (locc_test) then
!            dmft
             if(paw_dmft%use_dmft==1.and.dtset%nbandkss==0) then
               ebandldatot=ebandldatot+dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
               if(paw_dmft%band_in(iband)) then
                 ebandlda=ebandlda+dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
                 ekinlda=ekinlda+dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
                 occ(bdtot_index)=paw_dmft%occnd(1,iband,iband,ikpt,isppol)
                 ebanddmft=ebanddmft+dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
                 ekindmft=ekindmft+dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
               end if
             end if

             energies%e_eigenvalues = energies%e_eigenvalues + &
&             dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
             energies%e_kinetic = energies%e_kinetic + &
&             dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
             energies%e_nlpsp_vfock = energies%e_nlpsp_vfock + &
&             dtset%wtk(ikpt)*occ(bdtot_index)*enlxnk(bdtot_index)
             if (usefock) then
               energies%e_fock=energies%e_fock + half*focknk(bdtot_index)*occ(bdtot_index)*dtset%wtk(ikpt)
               if (optforces>0) fock%fock_common%forces(:,:)=fock%fock_common%forces(:,:)+&
&               dtset%wtk(ikpt)*occ(bdtot_index)*fockfornk(:,:,bdtot_index)
             end if
             if (optforces>0) grnl(:)=grnl(:)+dtset%wtk(ikpt)*occ(bdtot_index)*grnlnk(:,bdtot_index)
           end if
           bdtot_index=bdtot_index+1
           if(paw_dmft%use_dmft==1.and.dtset%nbandkss==0) then
             do iband1=1,nband_k
               if(paw_dmft%band_in(iband).and.paw_dmft%band_in(iband1)) then
!                write(std_out,*) "II+", isppol,ikpt,iband,iband1
                 ekindmft2=ekindmft2  +  dtset%wtk(ikpt)*paw_dmft%occnd(1,iband,iband1,ikpt,isppol)*&
&                 eknk_nd(isppol,ikpt,1,iband,iband1)
                 ekindmft2=ekindmft2  -  dtset%wtk(ikpt)*paw_dmft%occnd(2,iband,iband1,ikpt,isppol)*&
&                 eknk_nd(isppol,ikpt,2,iband,iband1)
!                write(std_out,*) "II", occnd(1,iband,iband1,ikpt,isppol),eknk_nd(isppol,ikpt,iband,iband1)
               end if
             end do
           end if
         end do
       end do
     end do

     if(paw_dmft%use_dmft==1) then
       energies%e_kinetic = energies%e_kinetic -ekindmft+ekindmft2
       if(abs(dtset%pawprtvol)>=2) then
         write(message,'(4a,7(2x,a,2x,e14.7,a),a)') &
           "-----------------------------------------------",ch10,&
           "--- Energy for DMFT and tests (in Ha)  ",ch10,&
           "--- Ebandldatot    (Ha.) = ",ebandldatot,ch10,&
           "--- Ebandlda       (Ha.) = ",ebandlda,ch10,&
           "--- Ebanddmft      (Ha.) = ",ebanddmft,ch10,&
           "--- Ekinlda        (Ha.) = ",ekinlda,ch10, &
           "--- Ekindmftdiag   (Ha.) = ",ekindmft,ch10,&
           "--- Ekindmftnondiag(Ha.) = ",ekindmft2,ch10,&
           "--- Edmft=         (Ha.) = ",edmft,ch10,&
           "-----------------------------------------------"
         call wrtout(std_out,message)
       end if
!       if(paw_dmft%use_dmft==1.and.mpi_enreg%paral_kgb==1) paw_dmft%use_dmft=0
     end if

     if (psps%usepaw==0) then
       call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
         rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
     else
       call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
         rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
     end if
     call timab(992,2,tsec)

!    Treat fixed occupation numbers or non-self-consistent case
   else

     if (mpi_enreg%nproc_kpt>1) then

       call timab(989,1,tsec)

       nbuf=2*mbdkpsp+dtset%nfft*dtset%nspden+3+3*natom*optforces
       ! If Hartree-Fock calculation, the exact exchange energy is k-dependent.
       if(dtset%usefock==1) then
         nbuf=nbuf+1
         if (optforces>0) nbuf=nbuf+3*natom
       end if
       if(iscf==-1 .or. iscf==-2)nbuf=2*mbdkpsp
       ABI_ALLOCATE(buffer1,(nbuf))
       ! Pack eigen,resid,rho[wf]r,grnl,enl,ek
       buffer1(1:mbdkpsp)=eigen(:)
       buffer1(1+mbdkpsp:2*mbdkpsp)=resid(:)
       index1=2*mbdkpsp
       if(iscf/=-1 .and. iscf/=-2)then
         if (psps%usepaw==0) then
           buffer1(index1+1:index1+dtset%nfft*dtset%nspden)=reshape(rhor,(/dtset%nfft*dtset%nspden/))
         else
           buffer1(index1+1:index1+dtset%nfft*dtset%nspden)=reshape(rhowfr,(/dtset%nfft*dtset%nspden/))
         end if
         index1=index1+dtset%nfft*dtset%nspden
         buffer1(index1+1) = energies%e_kinetic
         buffer1(index1+2) = energies%e_eigenvalues
         buffer1(index1+3) = energies%e_nlpsp_vfock
         index1=index1+3
         ! If Hartree-Fock calculation, save e_fock in buffer1
         if (dtset%usefock==1) then
           buffer1(index1+1) = energies%e_fock
           index1=index1+1
           if (optforces>0)then
             buffer1(index1+1:index1+3*natom)=reshape(fock%fock_common%forces,(/3*natom/))
             index1=index1+3*natom
           end if
         end if
         if (optforces>0) buffer1(index1+1:index1+3*natom)=grnl(1:3*natom)
       end if

       ! Build sum of everything
       call timab(48,1,tsec)
       call xmpi_sum(buffer1,nbuf,mpi_enreg%comm_kpt ,ierr)
       call timab(48,2,tsec)

       ! Unpack the final result
       eigen(:)=buffer1(1:mbdkpsp)
       resid(:)=buffer1(1+mbdkpsp:2*mbdkpsp)
       index1=2*mbdkpsp
       if(iscf/=-1 .and. iscf/=-2)then
         if (psps%usepaw==0) then
           ii=1
           do ispden=1,dtset%nspden
             do ifft=1,dtset%nfft
               rhor(ifft,ispden)=buffer1(index1+ii)
               ii=ii+1
             end do
           end do
         else
           ii=1
           do ispden=1,dtset%nspden
             do ifft=1,dtset%nfft
               rhowfr(ifft,ispden)=buffer1(index1+ii)
               ii=ii+1
             end do
           end do
         end if
         index1=index1+dtset%nfft*dtset%nspden
         energies%e_kinetic = buffer1(index1+1)
         energies%e_eigenvalues = buffer1(index1+2)
         energies%e_nlpsp_vfock = buffer1(index1+3)
         index1=index1+3
         ! If Hartree-Fock calculation, save e_fock in buffer1
         if (dtset%usefock==1) then
           energies%e_fock = buffer1(index1+1)
           index1=index1+1
           if (optforces>0) then
             fock%fock_common%forces(:,:)=reshape(buffer1(index1+1:index1+3*natom),(/3,natom/))
             index1=index1+3*natom
           end if
         end if
         if (optforces>0) grnl(1:3*natom)=buffer1(index1+1:index1+3*natom)
       end if
       ABI_DEALLOCATE(buffer1)
       call timab(989,2,tsec)

     end if ! nproc_kpt>1


!    Compute the highest occupied eigenenergy
     if(iscf/=-1 .and. iscf/=-2)then
       call timab(993,1,tsec)
       energies%e_fermie = -huge(one)
       bdtot_index=1
       do isppol=1,dtset%nsppol
         do ikpt=1,dtset%nkpt
           nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
           do iband=1,nband_k
             if(abs(occ(bdtot_index))>tol8 .and. eigen(bdtot_index)>energies%e_fermie+tol10) then
               energies%e_fermie=eigen(bdtot_index)
             end if
             bdtot_index=bdtot_index+1
           end do
         end do
       end do
       call xmpi_max(energies%e_fermie,spaceComm_distrb,ierr)
       call timab(993,2,tsec)
     end if

!    If needed, compute rhog, and symmetrizes the density
     if (iscf > 0 .or. iscf==-3 ) then
!      energies%e_fermie=zero  ! Actually, should determine the maximum of the valence band XG20020802
       nfftot=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)

       call timab(994,1,tsec)
       if (psps%usepaw==0) then
         call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,&
&         dtset%nsppol,dtset%nsym,phnons,rhog,rhor,rprimd,dtset%symafm,dtset%symrel,dtset%tnons)
       else
         call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,&
&         dtset%nsppol,dtset%nsym,phnons,rhowfg,rhowfr,rprimd,dtset%symafm,dtset%symrel,dtset%tnons)
       end if
       call timab(994,2,tsec)
!      We now have both rho(r) and rho(G), symmetrized, and if dtset%nsppol=2
!      we also have the spin-up density, symmetrized, in rhor(:,2).
     end if

   end if !  End of test on varying or fixed occupation numbers

   call timab(994,1,tsec)

!  Compute the kinetic energy density
   if(dtset%usekden==1 .and. (iscf > 0 .or. iscf==-3 ) )then
     if (psps%usepaw==0) then
       call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&       taug,taur,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
     else
       call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&      tauwfg,tauwfr,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
     end if
   end if

   ABI_DEALLOCATE(eknk)
   if (usefock) then
     ABI_DEALLOCATE(focknk)
     if (optforces>0)then
       ABI_DEALLOCATE(fockfornk)
     end if
   end if
   ABI_DEALLOCATE(eknk_nd)
   ABI_DEALLOCATE(grnlnk)
   ABI_DEALLOCATE(enlxnk)

!  In the non-self-consistent case, print eigenvalues and residuals
   if(iscf<=0 .and. me_distrb==0)then
     option=2 ; enunit=1 ; vxcavg_dum=zero
     call prteigrs(eigen,enunit,energies%e_fermie,dtfil%fnameabo_app_eig,&
&     ab_out,iscf,dtset%kptns,dtset%kptopt,dtset%mband,dtset%nband,&
&     dtset%nkpt,nnsclo_now,dtset%nsppol,occ,dtset%occopt,option,&
&     dtset%prteig,prtvol,resid,dtset%tolwfr,vxcavg_dum,dtset%wtk)
   end if

!  Find largest residual over bands, k points, and spins, except for nbdbuf highest bands
   ibdkpt=1
   residm=zero
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       nband_eff=max(1,nband_k-dtset%nbdbuf)
       residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_eff-1)))
       ibdkpt=ibdkpt+nband_k
     end do
   end do

 end if !usewvl==0

!===================================================================
! End of PLANE WAVES section
!===================================================================

!In the self-consistent case, diagnose lack of unoccupied state (for each spin and k-point).
!Print a warning if the number of such messages already written does not exceed mwarning.
 mwarning=5
 if(nwarning<mwarning .and. iscf>=0)then
   nwarning=nwarning+1
   bdtot_index=1
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       min_occ=two
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       do iband=1,nband_k
         if(occ(bdtot_index)<min_occ)min_occ=occ(bdtot_index)
         bdtot_index=bdtot_index+1
       end do
       if(min_occ>0.01_dp)then
         if(dtset%nsppol==1)then
           write(message, '(a,i0,3a,f7.3,5a)' )&
             'For k-point number ',ikpt,',',ch10,&
             'The minimal occupation factor is',min_occ,'.',ch10,&
             'An adequate monitoring of convergence requires it to be  at most 0.01_dp.',ch10,&
             'Action: increase slightly the number of bands.'
         else
           write(message, '(a,i0,3a,i0,a,f7.3,5a)' )&
             'For k-point number ',ikpt,', and',ch10,&
             'for spin polarization ',isppol, ' the minimal occupation factor is',min_occ,'.',ch10,&
             'An adequate monitoring of convergence requires it to be at most 0.01_dp.',ch10,&
             'Action: increase slightly the number of bands.'
         end if
         MSG_WARNING(message)
         exit ! It is enough if one lack of adequate occupation is identified, so exit.
       end if
     end do
   end do
 end if

 if (iscf>0.or.iscf==-3 .or. (dtset%usewvl==1 .and. iscf==0)) then

!  PAW: Build new rhoij quantities from new occ then symetrize them
!  Compute and add the compensation density to rhowfr to get the total density
   if (psps%usepaw==1) then
     call timab(555,1,tsec)
     if (paral_atom) then
       ABI_DATATYPE_ALLOCATE(pawrhoij_unsym,(natom))
       call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&                nspden=dtset%nspden,spnorb=dtset%pawspnorb,cpxocc=dtset%pawcpxocc)
       call pawrhoij_alloc(pawrhoij_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&       dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0)
     else
       pawrhoij_unsym => pawrhoij
     end if
     if (usecprj_local==1) then
       call pawmkrhoij(atindx,atindx1,cprj,gs_hamk%dimcprj,dtset%istwfk,dtset%kptopt,dtset%mband,mband_cprj,&
&       mcprj_local,dtset%mkmem,mpi_enreg,natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
&       occ,dtset%paral_kgb,paw_dmft,pawrhoij_unsym,dtfil%unpaw,dtset%usewvl,dtset%wtk)
     else
       mcprj_tmp=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
       ABI_DATATYPE_ALLOCATE(cprj_tmp,(natom,mcprj_tmp))
       call pawcprj_alloc(cprj_tmp,0,gs_hamk%dimcprj)
       call ctocprj(atindx,cg,1,cprj_tmp,gmet,gprimd,0,0,0,dtset%istwfk,kg,dtset%kptns,&
&       mcg,mcprj_tmp,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&       dtset%natom,nattyp,dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg,&
&       npwarr,dtset%nspinor,dtset%nsppol,ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&       ucvol,dtfil%unpaw,xred,ylm,ylmgr_dum)
       call pawmkrhoij(atindx,atindx1,cprj_tmp,gs_hamk%dimcprj,dtset%istwfk,dtset%kptopt,&
&       dtset%mband,mband_cprj,mcprj_tmp,dtset%mkmem,mpi_enreg,natom,dtset%nband,dtset%nkpt,&
&       dtset%nspinor,dtset%nsppol,occ,dtset%paral_kgb,paw_dmft,pawrhoij_unsym,dtfil%unpaw,&
&       dtset%usewvl,dtset%wtk)
       call pawcprj_free(cprj_tmp)
       ABI_DATATYPE_DEALLOCATE(cprj_tmp)
     end if
     call timab(555,2,tsec)
!    Build symetrized packed rhoij and compensated pseudo density
     cplex=1;ipert=0;idir=0;qpt(:)=zero
     if(dtset%usewvl==0) then
       call pawmkrho(1,compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&       my_natom,natom,dtset%nspden,dtset%nsym,ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&       dtset%pawprtvol,pawrhoij,pawrhoij_unsym,pawtab,qpt,rhowfg,rhowfr,rhor,rprimd,dtset%symafm,&
&       symrec,dtset%typat,ucvol,dtset%usewvl,xred,pawnhat=nhat,rhog=rhog)
       if (dtset%usekden==1) then
!        DO WE NEED TAUG?
         call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,tauwfg,taug,tauwfr,taur)
       end if
     else
!      here do not pass rhog, we do not use it
       call pawmkrho(1,compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&       my_natom,natom,dtset%nspden,dtset%nsym,ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&       dtset%pawprtvol,pawrhoij,pawrhoij_unsym,pawtab,qpt,rhowfg,rhowfr,rhor,rprimd,dtset%symafm,&
&       symrec,dtset%typat,ucvol,dtset%usewvl,xred,pawnhat=nhat)
!      In WVL: copy density to BigDFT object:
       call wvl_rho_abi2big(1,rhor,wvl%den)
     end if
     if (paral_atom) then
       call pawrhoij_free(pawrhoij_unsym)
       ABI_DATATYPE_DEALLOCATE(pawrhoij_unsym)
     end if
   end if

   if(paw_dmft%use_dmft==1) then
!    == check noccmmp
     call setnoccmmp(1,0,rdum4,0,0,idum3,my_natom,natom,0,1,dtset%nsppol,0,ntypat,&
&     paw_ij,pawang,dtset%pawprtvol,pawrhoij,pawtab,rdum2,idum1,dtset%typat,0,dtset%usepawu,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if

!  Find and print minimum and maximum total electron density and locations
!  Compute density residual (if required) and its squared norm
   if (iscf>=0) then
     if (psps%usepaw==0) then
       call prtrhomxmn(std_out,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor,ucvol=ucvol)
     else
       call prtrhomxmn(std_out,mpi_enreg,nfftf,pawfgr%ngfft,dtset%nspden,1,rhor,ucvol=ucvol)
     end if
     if (optres==1) then
       nvresid=rhor-nvresid
       call sqnorm_v(1,nfftf,nres2,dtset%nspden,optres,nvresid,mpi_comm_sphgrid=mpi_comm_sphgrid)
       if (dtset%usekden==1) then
         if (optres==1) tauresid=taur-tauresid
       endif
     end if
   end if

 end if ! iscf>0 or iscf=-3

 if(psps%usepaw==1.and.(iscf>=0.or.iscf==-3))  then
   ABI_DEALLOCATE(rhowfr)
   ABI_DEALLOCATE(rhowfg)
   if (dtset%usekden==1) then
     ABI_DEALLOCATE(tauwfr)
     ABI_DEALLOCATE(tauwfg)
   end if
 end if

 call timab(994,2,tsec)

 if (iscf==-1) then
!  Eventually compute the excited states within tddft
   call timab(995,1,tsec)
   if (psps%usepaw==1) then
!    In case of PAW calculation, have to transfer kxc from the fine to the coarse grid:
     ABI_ALLOCATE(cgrkxc,(dtset%nfft,nkxc))
     do ikxc=1,nkxc
       call transgrid(1,mpi_enreg,1,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrkxc(:,ikxc),kxc(:,ikxc))
     end do
     call tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&     kg,cgrkxc,dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nfft,&
&     ngfftdiel,dtset%nkpt,nkxc,npwarr,dtset%nspinor,dtset%nsppol,occ,ucvol,wffnew)
     ABI_DEALLOCATE(cgrkxc)
   else
     call tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&     kg,kxc,dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nfft,&
&     ngfftdiel,dtset%nkpt,nkxc,npwarr,dtset%nspinor,dtset%nsppol,occ,ucvol,wffnew)
   end if
   call timab(995,2,tsec)

 else
!  Eventually compute the susceptibility matrix and the
!  dielectric matrix when istep_mix is equal to 1 or dielstrt
   call timab(996,1,tsec)
   computesusmat = dtset%testsusmat(dielop, dielstrt, istep_mix) !test if the matrix is to be computed
   if(computesusmat) then
     dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
     dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
     dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
     dielar(7)=dtset%diemix;if (iscf>=10) dielar(7)=dtset%diemixmag
     usetimerev=1
     if (psps%usepaw==1.and.dtset%pawspnorb>0.and.dtset%kptopt/=1.and.dtset%kptopt/=2) usetimerev=0
     neglect_pawhat=1-dtset%pawsushat
     call suscep_stat(atindx,atindx1,cg,cprj,dielar,&
&     gs_hamk%dimcprj,doccde,eigen,gbound_diel,gprimd,&
&     irrzondiel,dtset%istwfk,kg,kg_diel,lmax_diel,&
&     dtset%mband,mcg,mcprj_local,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,natom,dtset%nband,&
&     neglect_pawhat,nfftdiel,ngfftdiel,&
&     dtset%nkpt,npwarr,npwdiel,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,ntypat,&
&     occ,dtset%occopt,pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&
&     susmat,dtset%symafm,dtset%symrel,dtset%tnons,dtset%typat,ucvol,&
&     dtfil%unpaw,usecprj_local,psps%usepaw,usetimerev,dtset%wtk,ylmdiel)
   end if
   call timab(996,2,tsec)

 end if ! end condition on iscf

 call gs_hamk%free()

 if (psps%usepaw==1) then
   if (usecprj==0) then
     call pawcprj_free(cprj_local)
     ABI_DATATYPE_DEALLOCATE(cprj_local)
   end if
 end if

 if(dtset%usewvl==0) then
   ABI_DEALLOCATE(EigMin)
   ABI_DEALLOCATE(doccde)
#if defined HAVE_GPU_CUDA
   if(dtset%use_gpu_cuda==1) call gpu_finalize_ham_data()
#endif
 end if

 call timab(980,2,tsec)

 DBG_EXIT("COLL")

 contains
!!***

!!****f* ABINIT/wvl_nscf_loop
!! NAME
!!  wvl_nscf_loop
!!
!! FUNCTION
!!  Non-self-consistent field cycle in Wavelets
!!  See also "wvl_nscf_loop_bigdft"
!!
!! INPUTS
!!  nnsclo= number of non-self consistent field iterations
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine wvl_nscf_loop()

!Arguments ------------------------------------
! integer, intent(in)                    :: istep,mcprj,nfft,nnsclo
! real(dp), intent(inout)                :: residm
! type(dataset_type), intent(in)         :: dtset
! type(MPI_type), intent(in)             :: mpi_enreg
! type(energies_type), intent(inout)     :: energies
! type(wvl_data), intent(inout)          :: wvl
! !arrays
! real(dp), intent(inout)                :: xcart(3, dtset%natom)
! real(dp), dimension(6), intent(out)    :: strsxc
! type(pawcprj_type),dimension(dtset%natom,mcprj),intent(out)::cprj

!Local variables-------------------------------
 integer :: inonsc,ii
 integer,parameter :: iscf_=-1       !do not do a SCF cycle
 logical,parameter :: do_scf=.false. !do not do a SCF cycle
 logical,parameter :: wvlbigdft=.false.
 real(dp) :: dum,eexctx,eh,ekin,eloc,enl,esicdc,evxc,exc

! *************************************************************************

   DBG_ENTER("COLL")

   if(nnsclo_now>0) then
     do inonsc=1,nnsclo_now
       call wvl_psitohpsi(dtset%diemix,eexctx,exc,eh,ekin,eloc,enl,esicdc,&
&       istep,inonsc,iscf_,mpi_enreg%me_wvl,dtset%natom,&
&       nfftf,mpi_enreg%nproc_wvl,dtset%nspden,&
&       dum,do_scf,evxc,wvl,wvlbigdft,xcart,strsxc)
       call wvl_hpsitopsi(cprj,dtset,energies,inonsc,mcprj_local,mpi_enreg, &
&       residm,wvl,xcart)
       if(residm<dtset%tolwfr) exit !Exit loop if converged
     end do

   else
     do ii=1, dtset%nline
!      Direct minimization technique: no diagonalization
       call wvl_psitohpsi(dtset%diemix,eexctx,exc,eh,ekin,eloc,enl,esicdc,&
&       istep,ii,iscf_,mpi_enreg%me_wvl,dtset%natom,&
&       nfftf,mpi_enreg%nproc_wvl,dtset%nspden,&
&       dum,do_scf,evxc,wvl,wvlbigdft,xcart,strsxc)
       call wvl_hpsitopsi(cprj,dtset,energies,ii,mcprj_local,mpi_enreg, &
&       residm,wvl,xcart)
       if(residm<dtset%tolwfr) exit !Exit loop if converged
     end do
   end if

!  Update energies depending on new WF
   energies%e_kinetic=ekin
   energies%e_nlpsp_vfock=enl
   energies%e_exactX=eexctx
   energies%e_sicdc=esicdc

!  Eventually update energies depending on density
   if (dtset%iscf<10) then
     energies%e_localpsp=eloc
     energies%e_hartree=eh
     energies%e_xc=exc ; energies%e_xcdc=evxc
   else if (nnsclo_now==0) then
     energies%e_localpsp=eloc
   end if

!  End of nscf iterations
   if (do_last_ortho) then
!    !Don't update energies (nscf cycle has been done); just recompute potential
     inonsc=nnsclo_now;if (nnsclo_now==0) inonsc=dtset%nline
     call wvl_psitohpsi(dtset%diemix,eexctx,exc,eh,ekin,eloc,enl,esicdc,&
&     istep,inonsc,iscf_,mpi_enreg%me_wvl,dtset%natom,&
&     nfftf,mpi_enreg%nproc_wvl,dtset%nspden,&
&     dum,do_scf,evxc,wvl,wvlbigdft,xcart,strsxc)
   end if

   DBG_EXIT("COLL")

 end subroutine wvl_nscf_loop
!!***

!!****f* ABINIT/wvl_nscf_loop_bigdft
!! NAME
!!  wvl_nscf_loop_bigdft
!!
!! FUNCTION
!!  Non-self-consistent field cycle in Wavelets
!!  It follows the BigDFT scheme.
!!  See also "wvl_nscf_loop"
!!
!! INPUTS
!!  nnsclo= number of non-self consistent field iterations
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine wvl_nscf_loop_bigdft()

!Arguments ------------------------------------
! integer, intent(in)                    :: istep,mcprj,nfft,nnsclo
! real(dp), intent(inout)                :: residm
! real(dp), intent(out)                  :: nres2
! type(dataset_type), intent(in)         :: dtset
! type(MPI_type), intent(in)             :: mpi_enreg
! type(energies_type), intent(inout)     :: energies
! type(wvl_data), intent(inout)          :: wvl
 !arrays
! real(dp), intent(inout)                :: xcart(3, dtset%natom)
! real(dp), dimension(6), intent(out)    :: strsxc
! type(pawcprj_type),dimension(dtset%natom,mcprj),intent(out)::cprj

!Local variables-------------------------------
 integer :: inonsc
 integer,parameter :: iscf_=-1       !do not do a SCF cycle
 logical,parameter :: do_scf=.false. !do not do a SCF cycle
 logical,parameter :: wvlbigdft=.true.
 real(dp) :: eexctx,eh,ekin,eloc,enl,esicdc,evxc,exc

! *************************************************************************

   DBG_ENTER("COLL")

   call wvl_hpsitopsi(cprj,dtset, energies, istep, mcprj_local,mpi_enreg, &
&   residm, wvl,xcart)

   if (nnsclo_now>2) then
     do inonsc = 2, nnsclo_now-1
       call wvl_psitohpsi(dtset%diemix, energies%e_exactX, energies%e_xc, &
&       energies%e_hartree, energies%e_kinetic, energies%e_localpsp, &
&       energies%e_nlpsp_vfock, energies%e_sicdc, istep, inonsc, iscf_, &
&       mpi_enreg%me_wvl, dtset%natom, nfftf, mpi_enreg%nproc_wvl,&
&       dtset%nspden, nres2, do_scf,energies%e_xcdc, &
&       wvl, wvlbigdft, xcart, strsxc)
       call wvl_hpsitopsi(cprj,dtset, energies, inonsc, mcprj_local,mpi_enreg, &
&       residm, wvl,xcart)
     end do
   end if

!  End of nscf iterations
   if (do_last_ortho.and.nnsclo_now<=1) then
!    !Don't update energies (nscf cycle has been done); just recompute potential
     call wvl_psitohpsi(dtset%diemix,eexctx,exc,eh,ekin,eloc,enl,esicdc, &
&     istep, 1, iscf_, mpi_enreg%me_wvl, dtset%natom, nfftf, &
&     mpi_enreg%nproc_wvl,dtset%nspden, nres2, do_scf,evxc, &
&     wvl, wvlbigdft, xcart, strsxc)
   else if (do_last_ortho.and.nnsclo_now>1) then
!    !Update energies and potential (nscf cycles are not finished)
     call wvl_psitohpsi(dtset%diemix, energies%e_exactX, energies%e_xc, &
&     energies%e_hartree,energies%e_kinetic, energies%e_localpsp, &
&     energies%e_nlpsp_vfock, energies%e_sicdc, istep, nnsclo_now, iscf_, &
&     mpi_enreg%me_wvl, dtset%natom, nfftf, mpi_enreg%nproc_wvl,&
&     dtset%nspden, nres2, do_scf,energies%e_xcdc, &
&     wvl, wvlbigdft, xcart, strsxc)
   end if

   DBG_EXIT("COLL")

 end subroutine wvl_nscf_loop_bigdft
!!***

!!****f* ABINIT/e_eigen
!! NAME
!!  e_eigen
!!
!! FUNCTION
!!  Computes eigenvalues energy from eigen, occ, kpt, wtk
!!
!! INPUTS
!!  eigen(nkpt*nsppol)=eigenvalues
!!  mband= maximum number of bands
!!  nband(nkpt*nsppol)= number of bands for each k-point and spin
!!  nkpt= number of k-points
!!  nsppol= number of spin polarization
!!  occ(mband*nkpt*nsppol)=occupations
!!  wtk(nkpt)= k-point weights
!!
!! OUTPUT
!!  e_eigenvalues= eigenvalues energy
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine e_eigen(eigen,e_eigenvalues,mband,nband,nkpt,nsppol,occ,wtk)

!Arguments ------------------------------------
 integer , intent(in)  :: mband,nkpt,nsppol
 integer , intent(in)  :: nband(nkpt*nsppol)
 real(dp) , intent(in)  :: eigen(mband*nkpt*nsppol)
 real(dp) , intent(in)  :: occ(mband*nkpt*nsppol)
 real(dp) , intent(in)  :: wtk(nkpt)
 real(dp) , intent(out) :: e_eigenvalues

!Local variables-------------------------------
 integer :: ib,iband,ii,ikpt,isppol,nband_k
 real(dp) :: wtk_k

! *************************************************************************

   DBG_ENTER("COLL")
   ii=0;ib=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       ii=ii+1
       nband_k=nband(ii) ;  wtk_k=wtk(ii)
       do iband=1,nband_k
         ib=ib+1
         if(abs(occ(ib)) > tol8) then
           e_eigenvalues = e_eigenvalues + wtk_k*occ(ib)*eigen(ib)
         end if
       end do
     end do
   end do

   DBG_EXIT("COLL")

 end subroutine e_eigen
!!***

!!****f* ABINIT/wvl_occ
!! NAME
!!  wvl_occ
!!
!! FUNCTION
!!  Computes occupations for the wavelet case
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!! for the wvlbigdft case, see the routine 'wvl_occ_bigdft'
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine wvl_occ()

!Local variables-------------------------------
 real(dp):: doccde_(dtset%mband*dtset%nkpt*dtset%nsppol)
! *************************************************************************

   DBG_ENTER("COLL")

!  Compute the new occupation numbers from eigen
   call newocc(doccde_,eigen,energies%entropy,energies%e_fermie,dtset%spinmagntarget,&
&   dtset%mband,dtset%nband,dtset%nelect,dtset%nkpt,dtset%nspinor,&
&   dtset%nsppol,occ,dtset%occopt,prtvol,&
&   dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)

! Copy occupations and efermi to BigDFT variables
   call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,1,wvl%wfs)

#if defined HAVE_BIGDFT
!  Copy Fermi level to BigDFT variable:
   wvl%wfs%ks%orbs%efermi=energies%e_fermie
#endif

   DBG_EXIT("COLL")

 end subroutine wvl_occ
!!***

!!****f* ABINIT/wvl_occ_bigdft
!! NAME
!!  wvl_occ_bigdft
!!
!! FUNCTION
!!  Computes occupations for the wavelet case
!!  Using BigDFT routines
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! for the wvlbigdft case, see the routine 'wvl_occ_bigdft'
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine wvl_occ_bigdft()

! *************************************************************************

   DBG_ENTER("COLL")

! Transfer occopt from ABINIT to BigDFT
#if defined HAVE_BIGDFT
   occopt_bigdft=dtset%occopt
   call wvl_occopt_abi2big(occopt_bigdft,occopt_bigdft,1)

!Calculate occupations using BigDFT routine
   call evaltoocc(mpi_enreg%me_wvl, mpi_enreg%nproc_wvl, .false., &
&   dtset%tsmear, wvl%wfs%ks%orbs,  occopt_bigdft)

!Pass e_fermi from BigDFT object to ABINIT variable:
   energies%e_fermie = wvl%wfs%ks%orbs%efermi

!Copy occupations from BigDFT to ABINIT variables
   call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,2,wvl%wfs)
#endif

   DBG_EXIT("COLL")

 end subroutine wvl_occ_bigdft
!!***

!!****f* ABINIT/wvl_comm_eigen
!! NAME
!!  wvl_comm_eigen
!!
!! FUNCTION
!!  Computes occupations for the wavelet case
!!  Using BigDFT routines
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! for the wvlbigdft case, see the routine 'wvl_occ_bigdft'
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine wvl_comm_eigen()

!Arguments ------------------------------------

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer:: ikpt,norb,shift
#endif

! *************************************************************************

   DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
   if(wvlbigdft) then
!  Communicates eigenvalues to all procs.
!  This will print out the eigenvalues and Fermi level.
     call eigensystem_info(mpi_enreg%me_wvl, mpi_enreg%nproc_wvl,0.d0,&
&     wvl%wfs%ks%Lzd%Glr%wfd%nvctr_c+7*wvl%wfs%ks%Lzd%Glr%wfd%nvctr_f,&
&     wvl%wfs%ks%orbs,wvl%wfs%ks%psi)
   else
!  Send all eigenvalues to all procs.
!  I simply communicate eigenvalues: I do not print them into screen, nor calculate Fermi-level.
     norb=wvl%wfs%ks%orbs%norb
     if (mpi_enreg%nproc_wvl > 1) then
       shift=1
       do ikpt = 1, wvl%wfs%ks%orbs%nkpts
         call xmpi_bcast(wvl%wfs%ks%orbs%eval(shift:shift+norb-1),wvl%wfs%ks%orbs%ikptproc(ikpt),mpi_enreg%comm_wvl,ierr)
         shift=shift+norb
       end do
     end if
   end if

!Copy eigenvalues from BigDFT object to "eigen"
   call wvl_eigen_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,eigen,2,wvl%wfs)

#else
   BIGDFT_NOTENABLED_ERROR()
#endif

   DBG_EXIT("COLL")

 end subroutine wvl_comm_eigen

end subroutine vtorho
!!***

!!****f* ABINIT/cgq_builder
!! NAME
!! cgq_builder
!!
!! FUNCTION
!! This routine locates cgq for efield calculations, especially for parallel case
!!
!! INPUTS
!!  berryflag = logical flag determining use of electric field variables
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ikpt=index of current k kpt
!!  ikpt_loc=index of k point on current processor (see vtorho.F90)
!!  isspol=value of spin polarization currently treated
!!  me_distrb=current value from spaceComm_distrb (see vtorho.F90)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcgq=size of cgq array (see vtorho.F90)
!!  mkgq=size of pwnsfacq array (see vtorho.F90)
!!  my_nspinor=nspinor value determined by current // set up
!!  nband_k=number of bands at each k point
!!  nproc_distrb=nproc from spaceComm_distrb (see vtorho.F90)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  spaceComm_distrb=comm_cell from mpi_enreg
!!
!! OUTPUT
!!  cgq(2,mcgq)=planewave coefficients of wavenfunctions adjacent to cg at ikpt
!!  pwnsfacq(2,mkgq)=phase factors for non-symmorphic translations for cg's adjacent to cg(ikpt)
!!
!! SIDE EFFECTS
!! Input/Output
!!   dtefield <type(efield_type)> = efield variables
!!   mpi_enreg=information about MPI parallelization
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_recv,xmpi_send
!!
!! SOURCE

subroutine cgq_builder(berryflag,cg,cgq,dtefield,dtset,ikpt,ikpt_loc,isppol,mcg,mcgq,&
&                      me_distrb,mkgq,mpi_enreg,my_nspinor,nband_k,nproc_distrb,&
&                      npwarr,pwnsfac,pwnsfacq,pwind_alloc,spaceComm_distrb)

!Arguments ------------------------------------
 integer,intent(in) :: ikpt,ikpt_loc,isppol,me_distrb,mcg,mcgq,mkgq,my_nspinor,nband_k
 integer,intent(in) :: nproc_distrb,pwind_alloc,spaceComm_distrb
 logical,intent(in) :: berryflag
 type(dataset_type), intent(in) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(MPI_type), intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),pwnsfac(2,pwind_alloc)
 real(dp),intent(out) :: cgq(2,mcgq),pwnsfacq(2,mkgq)

!Local variables -------------------------
!scalars
 integer :: count,count1,icg1,icg2,dest,his_source
 integer :: idir,ierr,ifor,ikg1,ikg2,ikptf,ikpt1f,ikpt1i
 integer :: jkpt,jkpt1i,jkptf,jkpt1f,jsppol,my_source,npw_k1,tag
!arrays
 integer,allocatable :: flag_send(:,:), flag_receive(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer(:,:)

! *************************************************************************

 if (mcgq==0.or.mkgq==0) return

 call timab(983,1,tsec)

!Test compatbility of berryflag
 if (berryflag) then
   ABI_ALLOCATE(flag_send,(0:nproc_distrb-1,dtefield%fnkpt))
 end if
 ABI_ALLOCATE(flag_receive,(dtset%nkpt))
 flag_send(:,:) = 0
 flag_receive(:) = 0

 if (berryflag) ikptf = dtefield%i2fbz(ikpt)

 do idir = 1, 3

!  skip idir values for which efield_dot(idir) = 0
   if (berryflag .and. abs(dtefield%efield_dot(idir)) < tol12 ) cycle

   do ifor = 1, 2

     if(berryflag) then
       dtefield%sflag(:,ikpt + dtset%nkpt*(isppol - 1),ifor,idir) = 0
       ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
       ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)
     end if

     npw_k1 = npwarr(ikpt1i)
     count = npw_k1*my_nspinor*nband_k
     my_source = mpi_enreg%proc_distrb(ikpt1i,1,isppol)

     do dest = 0, nproc_distrb-1

       if ((dest==me_distrb).and.(ikpt_loc <= dtset%mkmem)) then
!        I am dest and have something to do

         if ( my_source == me_distrb ) then
!          I am destination and source

           if(berryflag) then
             ikg1 = dtefield%fkgindex(ikpt1f)
             ikg2 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
             icg1 = dtefield%cgindex(ikpt1i,isppol)
             icg2 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
           end if

           pwnsfacq(:,ikg2 + 1:ikg2 + npw_k1) = pwnsfac(:,ikg1 + 1:ikg1 + npw_k1)
           cgq(:,icg2 + 1:icg2 + count) = cg(:,icg1 + 1:icg1 + count)

         else !  I am the destination but not the source -> receive
!          receive pwnsfacq
           if(berryflag) then
             tag = ikpt1f + (isppol - 1)*dtefield%fnkpt
             ikg1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
           end if
           ABI_ALLOCATE(buffer,(2,npw_k1))
           call xmpi_recv(buffer,my_source,tag,spaceComm_distrb,ierr)
           pwnsfacq(:,ikg1+1:ikg1+npw_k1) = buffer(:,1:npw_k1)
           ABI_DEALLOCATE(buffer)

!          receive cgq if necessary
           if(flag_receive(ikpt1i) == 0) then
             ABI_ALLOCATE(buffer,(2,count))
             tag = ikpt1i + (isppol - 1)*dtset%nkpt
             call xmpi_recv(buffer,my_source,tag,spaceComm_distrb,ierr)
             if(berryflag) icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
             cgq(:,icg1+1:icg1+count) = buffer(:,1:count)
             ABI_DEALLOCATE(buffer)
             flag_receive(ikpt1i) = 1
           end if ! end if flag_receive == 0
         end if ! end tasks if I am the destination

       else if (ikpt_loc <= mpi_enreg%mkmem(dest)) then  ! dest != me and the dest has a k-point to treat

!        jkpt is the kpt which is being treated by dest (in ibz)
!        jsppol is his isppol
         jkpt = mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,1)
         jsppol = mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,2)

         if(jkpt > 0 .and. jsppol > 0) then

           if(berryflag) then
             jkptf = dtefield%i2fbz(jkpt)
             jkpt1f = dtefield%ikpt_dk(jkptf,ifor,idir)
             jkpt1i = dtefield%indkk_f2ibz(jkpt1f,1)
           end if
           his_source = mpi_enreg%proc_distrb(jkpt1i,1,jsppol)

           if (his_source == me_distrb) then

!            send
!            pwnsfacq
             if(berryflag) then
               ikg1 = dtefield%fkgindex(jkpt1f)
               tag = jkpt1f + (jsppol - 1)*dtefield%fnkpt
             end if
             count1 = npwarr(jkpt1i)
             ABI_ALLOCATE(buffer,(2,count1))
             buffer(:,1:count1)  = pwnsfac(:,ikg1+1:ikg1+count1)
             call xmpi_send(buffer,dest,tag,spaceComm_distrb,ierr)
             ABI_DEALLOCATE(buffer)

!            send cgq if necessary
             if(flag_send(dest, jkpt1i)==0) then
               if(berryflag) icg1 = dtefield%cgindex(jkpt1i,jsppol)
               tag = jkpt1i + (jsppol - 1)*dtset%nkpt
               count1 = npwarr(jkpt1i)*nband_k*my_nspinor
               ABI_ALLOCATE(buffer,(2,count1))
               buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)
               call xmpi_send(buffer,dest,tag,spaceComm_distrb,ierr)
               ABI_DEALLOCATE(buffer)
               flag_send(dest, jkpt1i)=1
             end if ! if send cgq

           end if ! end check that his_source == me

         end if ! end check on jkpt > 0 and jsppol > 0

       end if ! end check on me = dest else if me != dest

     end do ! end loop over dest = 0, nproc-1

   end do !end loop over ifor

 end do !end loop over idir

 call timab(983,2,tsec)

 ABI_DEALLOCATE(flag_send)
 ABI_DEALLOCATE(flag_receive)

end subroutine cgq_builder
!!***

end module m_vtorho
!!***

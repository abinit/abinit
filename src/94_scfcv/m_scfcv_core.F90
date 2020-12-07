!!****m* ABINIT/m_scfcv_core
!! NAME
!!  m_scfcv_core
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
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

module m_scfcv_core

 use defs_basis
 use defs_wvltypes
 use defs_rectypes
 use m_xmpi
 use m_abicore
 use m_wffile
 use m_rec
 use m_ab7_mixing
 use m_errors
 use m_efield
 use mod_prc_memory
 use m_nctk
 use m_hdr
 use m_xcdata
 use m_cgtools
 use m_dtfil
 use m_distribfft
 use m_nonlop,           only : nonlop_counter


 use defs_datatypes,     only : pseudopotential_type
 use defs_abitypes,      only : MPI_type
 use m_berryphase_new,   only : update_e_field_vars
 use m_dens,             only : constrained_dft_t, constrained_dft_ini, constrained_dft_free
 use m_time,             only : timab
 use m_fstrings,         only : int2char4, sjoin, itoa
 use m_symtk,            only : symmetrize_xred
 use m_geometry,         only : metric
 use m_fftcore,          only : getng, sphereboundary
 use m_time,             only : abi_wtime, sec2str
 use m_exit,             only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_mpinfo,           only : destroy_mpi_enreg, iwrite_fftdatar, initmpi_seq, proc_distrb_cycle
 use m_ioarr,            only : fftdatar_write_from_hdr
 use m_results_gs ,      only : results_gs_type
 use m_scf_history,      only : scf_history_type, scf_history_init, scf_history_free
 use m_energies,         only : energies_type, energies_init, energies_copy
 use m_electronpositron, only : electronpositron_type, electronpositron_calctype
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type,pawtab_get_lsize
 use m_paw_an,           only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify, paw_an_reset_flags
 use m_pawfgrtab,        only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_get, pawcprj_lincom, &
&                               pawcprj_free, pawcprj_axpby, pawcprj_put, pawcprj_getdim, pawcprj_reorder
 use m_pawdij,           only : pawdij, symdij
 use m_pawfgr,           only : pawfgr_type
 use m_paw_ij,           only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_reset_flags
 use m_paw_dmft,         only : paw_dmft_type
 use m_paw_nhat,         only : nhatgrid,wvl_nhatgrid,pawmknhat
 use m_paw_tools,        only : chkpawovlp
 use m_paw_denpot,       only : pawdenpot
 use m_paw_occupancies,  only : pawmkrhoij
 use m_paw_correlations, only : setnoccmmp,setrhoijpbe0
 use m_orbmag,           only : orbmag_type
 use m_paw_mkrho,        only : pawmkrho
 use m_paw_uj,           only : pawuj_red, macro_uj_type
 use m_paw_dfpt,         only : pawgrnl
 use m_fock,             only : fock_type, fock_init, fock_destroy, fock_ACE_destroy, fock_common_destroy, &
                                fock_BZ_destroy, fock_update_exc, fock_updatecwaveocc
 use m_gwls_hamiltonian, only : build_vxc
#if defined HAVE_BIGDFT
 use BigDFT_API,         only : cprj_clean,cprj_paw_alloc
#endif
 use m_outxml,           only : out_resultsgs_XML, out_geometry_XML
 use m_kg,               only : getcut, getmpw, kpgio, getph
 use m_fft,              only : fourdp
 use m_vtorhorec,        only : first_rec, vtorhorec
 use m_vtorhotf,         only : vtorhotf
 use m_outscfcv,         only : outscfcv
 use m_afterscfloop,     only : afterscfloop
 use m_extraprho,        only : extraprho
 use m_spacepar,         only : make_vectornd,setsym
 use m_newrho,           only : newrho
 use m_newvtr,           only : newvtr
 use m_vtorho,           only : vtorho
 use m_setvtr,           only : setvtr
 use m_mkrho,            only : mkrho
 use m_rhotov,           only : rhotov
 use m_forces,           only : fresid, forces
 use m_dft_energy,       only : energy
 use m_initylmg,         only : initylmg
 use m_rhotoxc,          only : rhotoxc
 use m_drivexc,          only : check_kxc
 use m_odamix,           only : odamix
 use m_common,           only : scprqt, prtene
 use m_fourier_interpol, only : transgrid
 use m_fock_getghc,      only : fock2ACE
 use m_forstr,           only : nres2vres
 use m_positron,         only : setup_positron
 use m_cgprj,            only : ctocprj
 use m_psolver,          only : psolver_rhohxc
 use m_paw2wvl,          only : paw2wvl_ij, wvl_cprjreorder

 implicit none

 private
!!***

 public :: scfcv_core
!!***

contains
!!***

!!****f* ABINIT/scfcv_core
!! NAME
!! scfcv_core
!!
!! FUNCTION
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cpus= cpu time limit in seconds
!!  dmatpawu= fixed occupation matrix of correlated orbitals (DFT+U or DMFT only)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points treated by this node.
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |    for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  fatvshift=factor to multiply dtset%atvshift
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mcg=size of wave-functions array (cg) =mpw*my_nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)= # atoms of each type.
!!  ndtpawuj=size of dtpawuj
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and
!!     related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for
!!     each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real
!!     spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions; if mkmem>=nkpt, these are kept in a disk file.
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each NL proj |p_lmn>
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!  dtpawuj(ndtpawuj)= data used for the automatic determination of U
!!     (relevant only for PAW+U) calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for
!!     the electron-positron annihilation
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!     for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!     forces and its components, the stress tensor) of a ground-state
!!     computation (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous
!!     SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic
!!     energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wffnew=struct info for wf disk files
!!  wvl <type(wvl_data)>=all wavelets data
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic
!!     coordinates at output, current xred is transferred to xred_old
!!  conv_retcode=return code, 0 if convergence was achieved
!!
!! NOTES
!! It is worth to explain THE USE OF FFT GRIDS:
!! ============================================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj
!!      lincom_cgcprj,pawcprj_alloc,pawcprj_axpby,pawcprj_free,pawcprj_get
!!      pawcprj_getdim,pawcprj_lincom,pawcprj_put,timab,xmpi_sum,zgesv
!!
!! SOURCE

subroutine scfcv_core(itime, atindx,atindx1,cg,cprj,cpus,dmatpawu,dtefield,dtfil,dtorbmag,dtpawuj,&
&  dtset,ecore,eigen,electronpositron,fatvshift,hdr,indsym,&
&  initialized,irrzon,kg,mcg,mcprj,mpi_enreg,my_natom,nattyp,ndtpawuj,nfftf,npwarr,occ,&
&  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,&
&  pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,taug,taur,wffnew,wvl,xred,xred_old,ylm,ylmgr,conv_retcode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime, mcg,my_natom,ndtpawuj,pwind_alloc
 integer,intent(inout) :: initialized,nfftf,mcprj
 integer,intent(out) :: conv_retcode
 real(dp),intent(in) :: cpus,ecore,fatvshift
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(orbmag_type),intent(inout) :: dtorbmag
 type(electronpositron_type),pointer:: electronpositron
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(recursion_type),intent(inout) :: rec_set
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
 type(wffile_type),intent(inout) :: wffnew
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in) :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout) :: cg(2,mcg),dmatpawu(:,:,:,:)
 real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), pointer :: taug(:,:),taur(:,:)
 real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: xred(3,dtset%natom)
 real(dp), intent(inout) :: xred_old(3,dtset%natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
 type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawcprj_type),pointer, intent(inout) :: cprj(:,:)
!Local variables -------------------------
!scalars
 integer,parameter :: level=110,response=0,cplex1=1
 integer :: afford,bantot,choice
 integer :: computed_forces,cplex,cplex_hf,ctocprj_choice,dbl_nnsclo,dielop,dielstrt,dimdmat
 integer :: forces_needed,errid,has_vxctau,has_dijhat,has_dijnd,has_dijU,has_vhartree,has_dijfock
 integer :: history_size,usefock
 integer :: iatom,ider,idir,ierr,ii,ikpt,impose_dmat,denpot
 integer :: initialized0,iorder_cprj,ipert,ipositron,isave_den,isave_kden,iscf10,ispden
 integer :: ispmix,istep,istep_fock_outer,istep_mix,istep_updatedfock,itypat,izero,lmax_diel,lpawumax,mband_cprj
#if defined HAVE_BIGDFT
 integer :: mcprj_wvl
#endif
 integer :: me,me_wvl,mgfftdiel,mgfftf,moved_atm_inside,moved_rhor,my_nspinor,n1xccc
 integer :: n3xccc,ncpgr,nfftdiel,nfftmix,nfftmix_per_nfft,nfftotf,ngrcondft,ngrvdw,nhatgrdim,nk3xc,nkxc
 integer :: npawmix,npwdiel,nstep,nzlmopt,optcut,optcut_hf,optene,optgr0,optgr0_hf
 integer :: optgr1,optgr2,optgr1_hf,optgr2_hf,option,optrad,optrad_hf,optres,optxc,prtfor,prtxml,quit
 integer :: quit_sum,rdwrpaw,shft,spaceComm,spaceComm_fft,spaceComm_wvl,spaceComm_grid
 integer :: spare_mem
 integer :: stress_needed,sz1,sz2,tim_mkrho,unit_out
 integer :: usecprj,usexcnhat,use_hybcomp
 integer :: my_quit,quitsum_request,timelimit_exit,usecg,wfmixalg,with_vectornd
 integer ABI_ASYNC :: quitsum_async
 real(dp) :: boxcut,compch_fft,compch_sph,deltae,diecut,diffor,ecut
 real(dp) :: ecutf,ecutsus,edum,elast,etotal,evxc,fermie,gsqcut,hyb_mixing,hyb_mixing_sr
 real(dp) :: maxfor,res2,residm,ucvol,ucvol_local,val_max
 real(dp) :: val_min,vxcavg,vxcavg_dum
 real(dp) :: zion,wtime_step,now,prev
 character(len=10) :: tag
 character(len=500) :: MY_NAME = "scfcv_core"
 character(len=1500) :: msg
 !character(len=500) :: dilatmx_errmsg
 character(len=fnlen) :: fildata
 type(MPI_type) :: mpi_enreg_diel
 type(xcdata_type) :: xcdata
 type(energies_type) :: energies
 type(ab7_mixing_object) :: mix,mix_mgga
 logical,parameter :: VERBOSE=.FALSE.
 logical :: dummy_nhatgr
 logical :: finite_efield_flag=.false.
 logical :: non_magnetic_xc=.false.
 logical :: recompute_cprj=.false.,reset_mixing=.false.
 logical,save :: tfw_activated=.false.
 logical :: wvlbigdft=.false.
!type(energies_type),pointer :: energies_wvl  ! TO BE ACTIVATED LATER
!arrays
 integer :: ngfft(18),ngfftdiel(18),ngfftf(18),ngfftmix(18),npwarr_diel(1)
 integer :: npwtot_diel(1)
 integer, save :: scfcv_jdtset = 0 ! To simulate iapp behavior
 integer, save :: scfcv_itime = 1 ! To simulate iapp behavior
 integer,allocatable :: dimcprj(:),dimcprj_srt(:)
 integer,allocatable :: gbound_diel(:,:),irrzondiel(:,:,:),kg_diel(:,:)
 integer,allocatable :: l_size_atm(:)
 integer,allocatable :: indsym_dum(:,:,:),symrec_dum(:,:,:), rmm_diis_status(:,:,:)
 logical,pointer :: lmselect_ep(:,:)
 real(dp) :: dielar(7),dphase(3),dummy2(6),favg(3),gmet(3,3),gprimd(3,3)
 real(dp) :: kpt_diel(3),pel(3),pel_cg(3),pelev(3),pion(3),ptot(3),qpt(3),red_ptot(3) !!REC
 real(dp) :: rhodum(1),rmet(3,3),strsxc(6),strten(6),tollist(12)
 real(dp) :: tsec(2),vnew_mean(dtset%nspden),vres_mean(dtset%nspden)
 real(dp) :: efield_old_cart(3), ptot_cart(3)
 real(dp) :: red_efield2(3),red_efield2_old(3)
 real(dp) :: vpotzero(2)
! red_efield1(3),red_efield2(3) is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009) [[cite:Stengel2009]]
! red_efield1(3) for fixed ebar calculation, red_efield2(3) for fixed reduced d calculation, in mixed BC
! red_efieldbar_lc(3) is local reduced electric field, defined by Eq.(28) of Nat. Phys. suppl. (2009) [[cite:Stengel2009]]
! pbar(3) and dbar(3) are reduced polarization and displacement field,
!    defined by Eq.(27) and (29) Nat. Phys. suppl. (2009) [[cite:Stengel2009]]
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:)
 real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
 real(dp),allocatable :: grchempottn(:,:),grcondft(:,:),grewtn(:,:)
 real(dp),allocatable :: grhf(:,:),grnl(:),grvdw(:,:),grxc(:,:)
 real(dp),allocatable :: intgres(:,:),kxc(:,:),nhat(:,:),nhatgr(:,:,:),nvresid(:,:),nvtauresid(:,:)
 real(dp),allocatable :: ph1d(:,:),ph1ddiel(:,:),ph1df(:,:)
 real(dp),allocatable :: phnonsdiel(:,:,:),rhowfg(:,:),rhowfr(:,:),shiftvector(:)
 real(dp),allocatable :: susmat(:,:,:,:,:),synlgr(:,:)
 real(dp),allocatable :: vectornd(:,:),vhartr(:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),vxc_hybcomp(:,:),vxctau(:,:,:),workr(:,:),xccc3d(:),xcctau3d(:),ylmdiel(:,:)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:)
 type(scf_history_type) :: scf_history_wf
 type(constrained_dft_t) :: constrained_dft
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(pawrhoij_type),pointer :: pawrhoij_ep(:)
 type(fock_type),pointer :: fock
 type(pawcprj_type),allocatable, target :: cprj_local(:,:)

! *********************************************************************

 _IBM6("Hello, I'm running on IBM6")

 DBG_ENTER("COLL")

 call timab(238,1,tsec)
 call timab(54,1,tsec)

 ! enable time limit handler if not done in callers.
 if (enable_timelimit_in(MY_NAME) == MY_NAME) then
   write(std_out,*)"Enabling timelimit check in function: ",trim(MY_NAME)," with timelimit: ",trim(sec2str(get_timelimit()))
 end if

!######################################################################
!Initializations - Memory allocations
!----------------------------------------------------------------------
 lmax_diel = 0

!MPI communicators
 if (xmpi_paral==1.and.mpi_enreg%paral_hf==1) then
   spaceComm=mpi_enreg%comm_kpt
 else
   spaceComm=mpi_enreg%comm_cell
 end if
 me=xmpi_comm_rank(spaceComm)
 spaceComm_fft=mpi_enreg%comm_fft
 spaceComm_wvl=mpi_enreg%comm_wvl
 me_wvl=mpi_enreg%me_wvl
 spaceComm_grid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) spaceComm_grid=mpi_enreg%comm_wvl
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!Save some variables from dataset definition
 nstep=dtset%nstep
 ecut=dtset%ecut
 ecutf=ecut; if (psps%usepaw==1) ecutf=dtset%pawecutdg
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 iscf10=mod(dtset%iscf,10)
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff
 tollist(8)=dtset%vdw_df_threshold
 dielstrt=0
 finite_efield_flag=(dtset%berryopt == 4  .or. &
& dtset%berryopt == 6  .or. &
& dtset%berryopt == 7  .or. &
& dtset%berryopt == 14 .or. &
& dtset%berryopt == 16 .or. &
& dtset%berryopt == 17)

!Get FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
   mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
   mgfftf=dtset%mgfft;ngfftf(:)=ngfft(:)
 end if

!Calculate zion: the total positive charge acting on the valence electrons
 zion=zero
 do iatom=1,dtset%natom
   zion=zion+psps%ziontypat(dtset%typat(iatom))
 end do

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Fock: be sure that the pointer is initialized to Null.
 usefock=dtset%usefock
 nullify(fock)

!Special care in case of WVL
!wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)
!if (wvlbigdft) then   ! TO BE ACTIVATED LATER
!  ABI_ALLOCATE(energies_wvl,)
!end if
 ucvol_local = ucvol
#if defined HAVE_BIGDFT
 if (dtset%usewvl == 1) then
!  We need to tune the volume when wavelets are used because, not
!  all FFT points are used.
!  ucvol_local = (half * dtset%wvl_hgrid) ** 3 * ngfft(1)*ngfft(2)*ngfft(3)
   ucvol_local = product(wvl%den%denspot%dpbox%hgrids) * real(product(wvl%den%denspot%dpbox%ndims), dp)
 end if
#endif

!Some variables need to be initialized/nullified at start
 nullify(grhor,lrhor,elfr)
 quit=0 ; dbl_nnsclo=0 ; conv_retcode=0
 dielop=0 ; strsxc=zero
 deltae=zero ; elast=zero ;
 vpotzero(:)=zero
 ! JWZ April 12 2018: Intel 18 compiler seems to require maxfor initialized,
 ! else it dies in scprqt in some scenarios
 maxfor=zero
 !
 results_gs%residm=zero;results_gs%res2=zero
 results_gs%deltae=zero;results_gs%diffor=zero
 call energies_init(energies)
 if (dtset%positron/=0.and.initialized/=0) then
   energies%e0_electronpositron =results_gs%energies%e0_electronpositron
   energies%e_electronpositron  =results_gs%energies%e_electronpositron
   energies%edc_electronpositron=results_gs%energies%edc_electronpositron
   maxfor=zero
 end if

 ! Initialize fermi level.
 if (dtset%nstep==0 .or. dtset%iscf < 0) then
 !if ((dtset%nstep==0 .or. dtset%iscf < 0) .and. dtset%plowan_compute==0) then
   energies%e_fermie = results_gs%energies%e_fermie
   results_gs%fermie = results_gs%energies%e_fermie
   !write(std_out,*)"in scfcv_core: results_gs%fermie: ",results_gs%fermie
 end if

 select case(dtset%usepotzero)
 case(0,1)
   energies%e_corepsp   = ecore / ucvol
   energies%e_corepspdc = zero
 case(2)
   ! No need to include the PspCore energy since it is already included in the
   ! local pseudopotential  (vpsp)
   energies%e_corepsp   = zero
   energies%e_corepspdc = zero
 end select
 if(wvlbigdft) energies%e_corepsp = zero

 fermie=energies%e_fermie
 isave_den=0; isave_kden=0 !initial index of density protection file
 optres=merge(0,1,dtset%iscf<10)
 usexcnhat=0!;mcprj=0
 initialized0=initialized
 if (dtset%tfkinfunc==12) tfw_activated=.true.
 ipert=0;idir=0;cplex=1
 istep_mix=1
 istep_fock_outer=1
 ipositron=electronpositron_calctype(electronpositron)
 unit_out=0;if (dtset%prtvol >= 10) unit_out=ab_out
 nfftotf=product(ngfftf(1:3))

 usecprj=0
 if (mcprj>0) then
  usecprj=1
 end if

!Stresses and forces flags
 forces_needed=0;prtfor=0
 if ((dtset%optforces==1.or.dtset%ionmov==4.or.dtset%ionmov==5.or.abs(tollist(3))>tiny(0._dp))) then
   if (dtset%iscf>0.and.nstep>0) forces_needed=1
   if (nstep==0) forces_needed=2
   prtfor=1
 else if (dtset%iscf>0.and.dtset%optforces==2) then
   forces_needed=2
 end if

 stress_needed=0
 if (dtset%optstress>0.and.dtset%iscf>0.and.dtset%prtstm==0.and. (nstep>0.or.dtfil%ireadwf==1)) stress_needed=1
 if (dtset%optstress>0.and.dtset%iscf>0.and.psps%usepaw==1 &
& .and.finite_efield_flag.and.(nstep>0.or.dtfil%ireadwf==1)) stress_needed=1

!This is only needed for the tddft routine, and does not
!correspond to the intented use of results_gs (should be only
!for output of scfcv_core
 etotal = results_gs%etotal

!Entering a scfcv_core loop, printing data to XML file if required.
 prtxml=0;if (me==0.and.dtset%prtxml==1) prtxml=1
 if (prtxml == 1) then
!  scfcv_core() will handle a scf loop, so we output the scfcv markup.
   write(ab_xml_out, "(A)") '    <scfcvLoop>'
   write(ab_xml_out, "(A)") '      <initialConditions>'
!  We output the geometry of the dataset given in argument.
!  xred and rprimd are given independently since dtset only
!  stores original and final values.
   call out_geometry_XML(dtset, 4, dtset%natom, rprimd, xred)
   write(ab_xml_out, "(A)") '      </initialConditions>'
 end if

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor and res2 are here initialized to 0)
 choice=1 ; diffor=zero ; res2=zero
 ABI_ALLOCATE(fcart,(3,dtset%natom))
 ABI_ALLOCATE(fred,(3,dtset%natom))
 fred(:,:)=zero
 fcart(:,:)=results_gs%fcart(:,:) ! This is a side effect ...
!results_gs should not be used as input of scfcv_core
!HERE IS PRINTED THE FIRST LINE OF SCFCV

 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
& dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,istep_fock_outer,istep_mix,dtset%kptns,&
& maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
& occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
& psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode)

!Various allocations (potentials, gradients, ...)
 ABI_ALLOCATE(forold,(3,dtset%natom))
 ABI_ALLOCATE(grchempottn,(3,dtset%natom))
 ABI_ALLOCATE(grcondft,(3,dtset%natom))
 ABI_ALLOCATE(gresid,(3,dtset%natom))
 ABI_ALLOCATE(grewtn,(3,dtset%natom))
 ABI_ALLOCATE(grnl,(3*dtset%natom))
 ABI_ALLOCATE(grxc,(3,dtset%natom))
 ABI_ALLOCATE(synlgr,(3,dtset%natom))
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
 ABI_ALLOCATE(vhartr,(nfftf))
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vpsp,(nfftf))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))
 ABI_ALLOCATE(vxctau,(nfftf,dtset%nspden,4*dtset%usekden))

 wfmixalg=dtset%fockoptmix/100
 use_hybcomp=0
 if(mod(dtset%fockoptmix,100)==11)use_hybcomp=1
 ABI_ALLOCATE(vxc_hybcomp,(nfftf,dtset%nspden*use_hybcomp))

 ngrvdw=0;if (dtset%vdw_xc>=5.and.dtset%vdw_xc<=7) ngrvdw=dtset%natom
 ABI_ALLOCATE(grvdw,(3,ngrvdw))

 ngrcondft=0
 if(any(dtset%constraint_kind(:)/=0)) ngrcondft=dtset%natom
 ABI_ALLOCATE(intgres,(dtset%nspden,ngrcondft))
 if(ngrcondft/=0)then
   intgres(:,:)=zero
 endif

 grchempottn(:,:)=zero
 grcondft(:,:)=zero
 forold(:,:)=zero ; gresid(:,:)=zero ; pel(:)=zero
 vtrial(:,:)=zero; vxc(:,:)=zero
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))

!Allocations/initializations for PAW only
 lpawumax=-1
 if(psps%usepaw==1) then
!  Variables/arrays related to the fine FFT grid
   ABI_ALLOCATE(xcctau3d,(nfftf*dtset%usekden))
   ABI_ALLOCATE(nhat,(nfftf,dtset%nspden*psps%usepaw))
   if (nstep==0) nhat=zero
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(my_natom))
   if (my_natom>0) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,dtset%typat,&
&     mpi_atmtab=mpi_enreg%my_atmtab)
     call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat,&
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     ABI_DEALLOCATE(l_size_atm)
   end if
   compch_fft=-1.d5
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   if (usexcnhat==0.and.dtset%ionmov==4.and.dtset%iscf<10) then
     MSG_ERROR('You cannot simultaneously use ionmov=4 and such a PAW psp file !')
   end if

!  Variables/arrays related to the PAW spheres
   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   call paw_an_nullify(paw_an)
   call paw_ij_nullify(paw_ij)
   has_dijhat=0;if (dtset%iscf==22) has_dijhat=1
   has_vhartree=0; if (dtset%prtvha > 0 .or. dtset%prtvclmb > 0) has_vhartree=1
   has_dijfock=0; if (usefock==1) has_dijfock=1
   has_dijnd=0;if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
   has_dijU=merge(0,1,dtset%usepawu>0) !Be careful on this!
   has_vxctau=dtset%usekden
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,0,0,dtset%nspden,&
&   cplex,dtset%pawxcdev,dtset%typat,pawang,pawtab,has_vxc=1,&
&   has_vxctau=has_vxctau,has_vxc_ex=1,has_vhartree=has_vhartree,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,&
&   dtset%pawspnorb,dtset%natom,dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijfock=has_dijfock,has_dijhartree=1,has_dijnd=has_dijnd,has_dijso=1,has_dijhat=has_dijhat,&
&   has_dijU=has_dijU,has_pawu_occ=1,has_exexch_pot=1,nucdipmom=dtset%nucdipmom,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   if(dtset%usewvl==1) then
     call paw2wvl_ij(1,paw_ij,wvl%descr)
   end if
   compch_sph=-1.d5
   ABI_ALLOCATE(dimcprj,(dtset%natom))
   ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
   call pawcprj_getdim(dimcprj    ,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')
   call pawcprj_getdim(dimcprj_srt,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
   do itypat=1,dtset%ntypat
     if (pawtab(itypat)%usepawu/=0) lpawumax=max(pawtab(itypat)%lpawu,lpawumax)
   end do
   if (dtset%usedmatpu/=0.and.lpawumax>0) then
     if (2*lpawumax+1/=size(dmatpawu,1).or.2*lpawumax+1/=size(dmatpawu,2)) then
       MSG_BUG('Incorrect size for dmatpawu!')
     end if
   end if

!  Allocation of projected WF (optional)
   if (usecprj==1) then
     iorder_cprj=0
     if (usefock==1) then
       ctocprj_choice = 1
       if (dtset%optforces == 1) then
        ctocprj_choice = 2; ! ncpgr = 3
       end if
!       if (dtset%optstress /= 0) then
!         ncpgr = 6 ; ctocprj_choice = 3
!       end if
     end if

#if defined HAVE_BIGDFT
     if (dtset%usewvl==1) then
       mband_cprj=dtset%mband;if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
       mcprj_wvl=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
       ABI_DATATYPE_ALLOCATE(wvl%descr%paw%cprj,(dtset%natom,mcprj_wvl))
       call cprj_paw_alloc(wvl%descr%paw%cprj,0,dimcprj_srt)
     end if
#endif
   end if

!  Other variables for PAW
   nullify(pawrhoij_ep);if(associated(electronpositron))pawrhoij_ep=>electronpositron%pawrhoij_ep
   nullify(lmselect_ep);if(associated(electronpositron))lmselect_ep=>electronpositron%lmselect_ep
 else
   ABI_ALLOCATE(dimcprj,(0))
   ABI_ALLOCATE(dimcprj_srt,(0))
   ABI_ALLOCATE(nhat,(0,0))
   ABI_ALLOCATE(xcctau3d,(0))
   ABI_DATATYPE_ALLOCATE(paw_ij,(0))
   ABI_DATATYPE_ALLOCATE(paw_an,(0))
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(0))
 end if ! PAW

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (dtset%iscf>=0) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
   dielar(7)=dtset%diemix;if (dtset%iscf>=10) dielar(7)=dtset%diemixmag
   ABI_ALLOCATE(nvresid,(nfftf,dtset%nspden))
   ABI_ALLOCATE(nvtauresid,(nfftf,dtset%nspden*dtset%usekden))
   if (nstep==0) then
    nvresid=zero
    nvtauresid=zero
   end if
   ABI_ALLOCATE(dtn_pc,(3,dtset%natom))
!  The next arrays are needed if iscf==5 and ionmov==4,
!  but for the time being, they are always allocated
   ABI_ALLOCATE(grhf,(3,dtset%natom))
!  Additional allocation for mixing within PAW
   npawmix=0
   if(psps%usepaw==1) then
     do iatom=1,my_natom
       itypat=pawrhoij(iatom)%itypat
       pawrhoij(iatom)%use_rhoijres=1
       sz1=pawrhoij(iatom)%cplex_rhoij*pawtab(itypat)%lmn2_size
       sz2=pawrhoij(iatom)%nspden
       ABI_ALLOCATE(pawrhoij(iatom)%rhoijres,(sz1,sz2))
       do ispden=1,pawrhoij(iatom)%nspden
         pawrhoij(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(pawrhoij(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       pawrhoij(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij(iatom)%nspden*pawtab(itypat)%lmnmix_sz &
&                     *pawrhoij(iatom)%cplex_rhoij*pawrhoij(iatom)%qphase
     end do
   end if
   if (dtset%iscf > 0) then
     denpot = AB7_MIXING_POTENTIAL
     if (dtset%iscf > 10) denpot = AB7_MIXING_DENSITY
     if (psps%usepaw==1.and.dtset%pawmixdg==0 .and. dtset%usewvl==0) then
       ispmix=AB7_MIXING_FOURRIER_SPACE;nfftmix=dtset%nfft;ngfftmix(:)=ngfft(:)
     else
       ispmix=AB7_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
     end if
     !TRangel: added to avoid segfaults with Wavelets
     nfftmix_per_nfft=0;if(nfftf>0) nfftmix_per_nfft=(1-nfftmix/nfftf)
     call ab7_mixing_new(mix, iscf10, denpot, ispmix, nfftmix, dtset%nspden, npawmix, errid, msg, dtset%npulayit)
     if (errid /= AB7_NO_ERROR) then
       MSG_ERROR(msg)
     end if
     if (dtset%usekden/=0) then
       if (dtset%useria==12345) then  ! This is temporary
         call ab7_mixing_new(mix_mgga, iscf10, denpot, ispmix, nfftmix, dtset%nspden, 0, errid, msg, dtset%npulayit)
       else
         call ab7_mixing_new(mix_mgga, 0, denpot, ispmix, nfftmix, dtset%nspden, 0, errid, msg, dtset%npulayit)
       end if
       if (errid /= AB7_NO_ERROR) then
         MSG_ERROR(msg)
       end if
     end if
     if (dtset%mffmem == 0) then
       call ab7_mixing_use_disk_cache(mix, dtfil%fnametmp_fft)
       if (dtset%usekden/=0.and.denpot==AB7_MIXING_DENSITY) &
&        call ab7_mixing_use_disk_cache(mix, dtfil%fnametmp_fft_mgga)
     end if
!   else if (dtset%iscf==0.and.dtset%usewvl==1) then
!     ispmix=AB7_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
   end if
 else
   ABI_ALLOCATE(nvresid,(0,0))
   ABI_ALLOCATE(nvtauresid,(0,0))
   ABI_ALLOCATE(dtn_pc,(0,0))
   ABI_ALLOCATE(grhf,(0,0))
 end if ! iscf>0

! Here initialize the datastructure constrained_dft, for constrained DFT calculations
! as well as penalty function constrained magnetization
 if(any(dtset%constraint_kind(:)/=0).or.dtset%magconon/=0)then
   call constrained_dft_ini(dtset%chrgat,constrained_dft,dtset%constraint_kind,dtset%magconon,dtset%magcon_lambda,&
&    mpi_enreg,dtset%natom,nfftf,ngfftf,dtset%nspden,dtset%ntypat,&
&    dtset%ratsm,dtset%ratsph,rprimd,dtset%spinat,dtset%typat,xred,dtset%ziontypat)
 endif

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT

 if( (nstep>0 .and. dtset%iscf>=0) .or. dtset%iscf==-1 ) then !MF

!  Here, for TDDFT, artificially set iprcel . Also set a variable to reduce the memory needs.
   afford=1
   if(dtset%iscf==-1) then
!    dtset%iprcel=21
     afford=0
   end if

!  First compute dimensions
   if(dtset%iprcel>=21 .or. dtset%iscf==-1)then
!    With dielop=1, the matrices will be computed when istep=dielstrt
!    With dielop=2, the matrices will be computed when istep=dielstrt and 1
     dielop=1
     if(dtset%iprcel>=41)dielop=2
     if((dtset%iprcel >= 71).and.(dtset%iprcel<=79)) dielop=0 !RSkerker preconditioner do not need the susceptibility matrix
!    Immediate computation of dielectric matrix
     dielstrt=1
!    Or delayed computation
     if(modulo(dtset%iprcel,100)>21 .and. modulo(dtset%iprcel,100)<=29)dielstrt=modulo(dtset%iprcel,100)-20
     if(modulo(dtset%iprcel,100)>31 .and. modulo(dtset%iprcel,100)<=39)dielstrt=modulo(dtset%iprcel,100)-30
     if(modulo(dtset%iprcel,100)>41 .and. modulo(dtset%iprcel,100)<=49)dielstrt=modulo(dtset%iprcel,100)-40
     if(modulo(dtset%iprcel,100)>51 .and. modulo(dtset%iprcel,100)<=59)dielstrt=modulo(dtset%iprcel,100)-50
     if(modulo(dtset%iprcel,100)>61 .and. modulo(dtset%iprcel,100)<=69)dielstrt=modulo(dtset%iprcel,100)-60
!    Get diecut, and the fft grid to be used for the susceptibility computation
     diecut=abs(dtset%diecut)
     if( dtset%diecut<0.0_dp )then
       ecutsus=ecut
     else
       ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
     end if
!    Impose sequential calculation
     ngfftdiel(1:3)=0 ; ngfftdiel(7)=100 ; ngfftdiel(9)=0; ngfftdiel(8)=dtset%ngfft(8);ngfftdiel(10:18)=0
     if(dtset%iscf==-1)ngfftdiel(7)=102

!    The dielectric stuff is performed in sequential mode; set mpi_enreg_diel accordingly
     call initmpi_seq(mpi_enreg_diel)
     call getng(dtset%boxcutmin,ecutsus,gmet,k0,mpi_enreg_diel%me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
&     mpi_enreg_diel%nproc_fft,dtset%nsym,mpi_enreg_diel%paral_kgb,dtset%symrel,&
&     use_gpu_cuda=dtset%use_gpu_cuda)
!    Update the fft distribution
     call init_distribfft_seq(mpi_enreg_diel%distribfft,'c',ngfftdiel(2),ngfftdiel(3),'all')

!    Compute the size of the dielectric matrix
     kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
     call getmpw(diecut,dtset%exchn2n3d,gmet,(/1/),kpt_diel,mpi_enreg_diel,npwdiel,1)
     lmax_diel=0
     if (psps%usepaw==1) then
       do ii=1,dtset%ntypat
         lmax_diel=max(lmax_diel,pawtab(ii)%lcut_size)
       end do
     end if
   else
     npwdiel=1
     mgfftdiel=1
     nfftdiel=1
     lmax_diel=0
     afford=0
   end if

!  Now, performs allocation
   ABI_ALLOCATE(dielinv,(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
   ABI_ALLOCATE(susmat,(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
   ABI_ALLOCATE(kg_diel,(3,npwdiel))
   ABI_ALLOCATE(gbound_diel,(2*mgfftdiel+8,2))
   ABI_ALLOCATE(irrzondiel,(nfftdiel**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
   ABI_ALLOCATE(phnonsdiel,(2,nfftdiel**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
   ABI_ALLOCATE(ph1ddiel,(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw))
   ABI_ALLOCATE(ylmdiel,(npwdiel,lmax_diel**2))
!  Then, compute the values of different arrays
   if(dielop>=1)then
!    Note : npwarr_diel is dummy, npwtot_diel is dummy
!    This kpgio call for going from the suscep FFT grid to the diel sphere
     npwarr_diel(1)=npwdiel

     call kpgio(diecut,dtset%exchn2n3d,gmet,(/1/),kg_diel,&
&     kpt_diel,1,(/1/),1,'COLL',mpi_enreg_diel,npwdiel,&
&     npwarr_diel,npwtot_diel,dtset%nsppol)
     call sphereboundary(gbound_diel,1,kg_diel,mgfftdiel,npwdiel)

     if (dtset%nsym>1 .and. dtset%iscf>=0 ) then
!      Should replace this initialization of irrzondiel and phnonsdiel through setsym by a direct call to irrzg
       ABI_ALLOCATE(indsym_dum,(4,dtset%nsym,dtset%natom))
       ABI_ALLOCATE(symrec_dum,(3,3,dtset%nsym))
       call setsym(indsym_dum,irrzondiel,dtset%iscf,dtset%natom,&
&       nfftdiel,ngfftdiel,dtset%nspden,dtset%nsppol,dtset%nsym,phnonsdiel,&
&       dtset%symafm,symrec_dum,dtset%symrel,dtset%tnons,dtset%typat,xred)
       ABI_DEALLOCATE(indsym_dum)
       ABI_DEALLOCATE(symrec_dum)
     end if
     if (psps%usepaw==1) then
       call getph(atindx,dtset%natom,ngfftdiel(1),ngfftdiel(2),&
&       ngfftdiel(3),ph1ddiel,xred)
       call initylmg(gprimd,kg_diel,kpt_diel,1,mpi_enreg_diel,&
&       lmax_diel,npwdiel,dtset%nband,1,npwarr_diel,dtset%nsppol,0,&
&       rprimd,ylmdiel,rhodum)
     end if
   end if

   if(dtset%iprcel>=21 .or. dtset%iscf==-1)then
     call destroy_mpi_enreg(mpi_enreg_diel)
   end if

 else
   npwdiel=1
   mgfftdiel=1
   nfftdiel=1
   afford = 0
   ABI_ALLOCATE(susmat,(0,0,0,0,0))
   ABI_ALLOCATE(kg_diel,(0,0))
   ABI_ALLOCATE(gbound_diel,(0,0))
   ABI_ALLOCATE(irrzondiel,(0,0,0))
   ABI_ALLOCATE(phnonsdiel,(0,0,0))
   ABI_ALLOCATE(ph1ddiel,(0,0))
   ABI_ALLOCATE(ylmdiel,(0,0))
 end if

 nkxc=0
!TDDFT - For a first coding
 if (dtset%iscf==-1 .and. dtset%nspden==1) nkxc=2
 if (dtset%iscf==-1 .and. dtset%nspden==2) nkxc=3
!Eventually need kxc-LDA when susceptibility matrix has to be computed
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=2*min(dtset%nspden,2)-1
!Eventually need kxc-LDA for residual forces (when density mixing is selected)
 if (dtset%iscf>=10.and.dtset%usewvl==0.and.forces_needed>0 .and. &
& abs(dtset%densfor_pred)>=1.and.abs(dtset%densfor_pred)<=6.and.abs(dtset%densfor_pred)/=5) then
   if (dtset%xclevel==1.or.dtset%densfor_pred>=0) nkxc=2*min(dtset%nspden,2)-1
   if (dtset%xclevel==2.and.dtset%nspden==1.and.dtset%densfor_pred<0) nkxc=7
   if (dtset%xclevel==2.and.dtset%nspden==2.and.dtset%densfor_pred<0) nkxc=19
 end if
 if (nkxc>0) then
   call check_kxc(dtset%ixc,dtset%optdriver)
 end if
 ABI_ALLOCATE(kxc,(nfftf,nkxc))

!This flag will be set to 1 just before an eventual change of atomic
!positions inside the iteration, and set to zero when the consequences
!of this change are taken into account.
 moved_atm_inside=0
!This flag will be set to 1 if the forces are computed inside the iteration.
 computed_forces=0

 if(dtset%wfoptalg==2)then
   ABI_ALLOCATE(shiftvector,((dtset%mband+2)*dtset%nkpt))
   val_min=-1.0_dp
   val_max=zero
 else
   ABI_ALLOCATE(shiftvector,(1))
 end if

!!PAW+DMFT: allocate structured datatype paw_dmft if dtset%usedmft=1
!call init_sc_dmft(dtset%dmftbandi,dtset%dmftbandf,dtset%mband,dtset%nkpt,&
!&  dtset%nsppol,dtset%usedmft,paw_dmft,dtset%usedmft)
!call print_sc_dmft(paw_dmft)

!!Electric field initializations: initialize pel_cg(:) and p_ion(:)
 call update_e_field_vars(atindx,atindx1,cg,dimcprj,dtefield,dtfil,dtset,&
& efield_old_cart,gmet,gprimd,hdr,idir,kg,mcg,&
& dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,dtset%natom,nattyp,ngfft,dtset%nkpt,npwarr,&
& dtset%ntypat,pawrhoij,pawtab,pel_cg,pelev,pion,psps,ptot,ptot_cart,&
& pwind,pwind_alloc,pwnsfac,red_efield2,red_efield2_old,red_ptot,rmet,rprimd,&
& 0,quit,istep,ucvol,unit_out,psps%usepaw,xred,ylm,ylmgr)

 if (dtset%iscf==22) energies%h0=zero

 call timab(54,2,tsec)

!##################################################################
!PERFORM ELECTRONIC ITERATIONS
!##################################################################

!Offer option of computing total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to vtorho
!Pass through the first routines even when nstep==0

 quitsum_request = xmpi_request_null; timelimit_exit = 0
 istep_updatedfock=0

 ABI_ICALLOC(rmm_diis_status, (2, dtset%nkpt, dtset%nsppol))

! start SCF loop
 do istep=1,max(1,nstep)

   ! Handle time limit condition.
   if (istep == 1) prev = abi_wtime()
   if (istep  > 1) then
     now = abi_wtime()
     wtime_step = now - prev
     prev = now
     call wrtout(std_out, sjoin("{SCF_istep:", itoa(istep-1), ", Vnl|psi>:", itoa(nonlop_counter), &
                  ", wall_time: '", sec2str(wtime_step), "'} <<< TIME"))
     nonlop_counter = 0

     if (have_timelimit_in(MY_NAME)) then
       if (istep > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then
           write(msg,"(3a)")"Approaching time limit ",trim(sec2str(get_timelimit())),". Will exit istep loop in scfcv_core."
           MSG_COMMENT(msg)
           call wrtout(ab_out, msg, "COLL")
           timelimit_exit = 1
           exit
         end if
       end if

       my_quit = 0; if (now - get_start_time() + 2.15 * wtime_step > get_timelimit()) my_quit = 1
       call xmpi_isum(my_quit,quitsum_async,spacecomm,quitsum_request,ierr)
     end if
   end if

   call timab(240,1,tsec)
   if (moved_atm_inside==1 .or. istep==1) then
!    ##############################################################
!    The following steps are done once for a given set of atomic
!    coordinates or for the nstep=1 case
!    --------------------------------------------------------------

!    Eventually symmetrize atomic coordinates over space group elements:
     call symmetrize_xred(dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred,indsym=indsym)

     if (dtset%usewvl == 0) then
!      Get cut-off for g-vectors
       if (psps%usepaw==1) then
         call wrtout(std_out,' FFT (fine) grid used in SCF cycle:','COLL')
       end if
       call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!      Compute structure factor phases and large sphere cut-off (gsqcut):
       call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)

       if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
         call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
       else
         ph1df(:,:)=ph1d(:,:)
       end if
     end if

!    Initialization of atomic data for PAW
     if (psps%usepaw==1) then
!      Check for non-overlapping spheres
       call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)

!      Identify parts of the rectangular grid where the density has to be calculated
       optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
       if (forces_needed==1.or.(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)) then
         optgr1=dtset%pawstgylm;if (stress_needed==1) optrad=1; if (dtset%pawprtwf==1) optrad=1
       end if

       if(dtset%usewvl==0) then
         call nhatgrid(atindx1,gmet,my_natom,dtset%natom,&
&         nattyp,ngfftf,psps%ntypat,optcut,optgr0,optgr1,optgr2,optrad,&
&         pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred,&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&         comm_fft=spaceComm_fft,distribfft=mpi_enreg%distribfft)
       else
         shft=0
#if defined HAVE_BIGDFT
         shft=wvl%descr%Glr%d%n1i*wvl%descr%Glr%d%n2i*wvl%den%denspot%dpbox%nscatterarr(me_wvl,4)
         call wvl_nhatgrid(atindx1,wvl%descr%atoms%astruct%geocode,&
&         wvl%descr%h,wvl%den%denspot%dpbox%i3s,dtset%natom,dtset%natom,&
&         nattyp,psps%ntypat,wvl%descr%Glr%d%n1,wvl%descr%Glr%d%n1i,&
&         wvl%descr%Glr%d%n2,wvl%descr%Glr%d%n2i,wvl%descr%Glr%d%n3,&
&         wvl%den%denspot%dpbox%n3pi,optcut,optgr0,optgr1,optgr2,optrad,&
&         pawfgrtab,pawtab,psps%gth_params%psppar,rprimd,shft,xred)
#endif
       end if
     end if

!    If we are inside SCF cycle or inside dynamics over ions,
!    we have to translate the density of previous iteration
     moved_rhor=0

     if (initialized/=0.and.dtset%usewvl == 0.and.ipositron/=1.and. &
&     (abs(dtset%densfor_pred)==2.or.abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)) then
       moved_rhor=1
       if (abs(dtset%densfor_pred)==2) then
         option=2
         ABI_ALLOCATE(workr,(nfftf,dtset%nspden))
         call fresid(dtset,gresid,mpi_enreg,nfftf,ngfftf,&
&         psps%ntypat,option,pawtab,rhor,rprimd,&
&         ucvol,workr,xred,xred_old,psps%znuclpsp)
         rhor=workr
         ABI_DEALLOCATE(workr)
       else if (abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6) then
         scf_history%icall=scf_history%icall+1
         call extraprho(atindx,atindx1,cg,cprj,dtset,gmet,gprimd,gsqcut,&
&         scf_history%icall,kg,mcg,mcprj,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&         my_natom,nattyp,nfftf,ngfftf,npwarr,psps%ntypat,pawrhoij,pawtab,&
&         ph1df,psps,psps%qgrid_vl,rhor,rprimd,scf_history,ucvol,&
&         psps%usepaw,xred,xred_old,ylm,psps%ziontypat,psps%znuclpsp)
       end if
       call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
     end if

     ! if any nuclear dipoles are nonzero, compute the vector potential in real space (depends on
     ! atomic position so should be done for nstep = 1 and for updated ion positions
     if ( any(abs(dtset%nucdipmom(:,:))>tol8) ) then
        with_vectornd = 1
     else
        with_vectornd = 0
     end if
     if(allocated(vectornd)) then
        ABI_DEALLOCATE(vectornd)
     end if
     ABI_ALLOCATE(vectornd,(with_vectornd*nfftf,3))
     if(with_vectornd .EQ. 1) then
        call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,dtset%nucdipmom,&
             & rprimd,vectornd,xred)
     endif

   end if ! moved_atm_inside==1 .or. istep==1

   !Initialize/Update data in the case of an Exact-exchange (Hartree-Fock) or hybrid XC calculation
   hyb_mixing=zero;hyb_mixing_sr=zero
   if (usefock==1) then
     if (istep==1) then
       ! Initialize data_type fock for the calculation
       cplex_hf=cplex
       if (psps%usepaw==1) cplex_hf=dtset%pawcpxocc
       call fock_init(atindx,cplex_hf,dtset,fock,gsqcut,kg,mpi_enreg,nattyp,npwarr,pawang,pawfgr,pawtab,rprimd)
       if (fock%fock_common%usepaw==1) then
         optcut_hf = 0 ! use rpaw to construct local_pawfgrtab
         optgr0_hf = 0; optgr1_hf = 0; optgr2_hf = 0 ! dont need gY terms locally
         optrad_hf = 1 ! do store r-R
         call nhatgrid(atindx1,gmet,dtset%natom,dtset%natom,nattyp,ngfftf,psps%ntypat,&
&         optcut_hf,optgr0_hf,optgr1_hf,optgr2_hf,optrad_hf,fock%fock_common%pawfgrtab,pawtab,&
&         rprimd,dtset%typat,ucvol,xred,typord=1)
         iatom=-1;idir=0
         call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,iatom,idir,&
&         iorder_cprj,dtset%istwfk,kg,dtset%kptns,mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft, dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,&
&         dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,dtfil%unpaw,&
&         xred,ylm,ylmgr)
       end if
       if(wfmixalg/=0)then
         spare_mem=0
         if(spare_mem==1)history_size=wfmixalg ! Not yet coded
         if(spare_mem==0)history_size=2*(wfmixalg-1)+1
!        Specific case of simple mixing : always history_size=1
         if(wfmixalg==2)history_size=1
         scf_history_wf%history_size=history_size
         usecg=2
         call scf_history_init(dtset,mpi_enreg,usecg,scf_history_wf)
       end if
     end if

     !Fock energy
     energies%e_exactX=zero
     if (fock%fock_common%optfor) fock%fock_common%forces=zero

     if (istep==1 .or. istep_updatedfock==fock%fock_common%nnsclo_hf .or. &
&        (fock%fock_common%nnsclo_hf>1 .and. fock%fock_common%scf_converged) ) then

       istep_updatedfock=1
       fock%fock_common%scf_converged=.false.

       !Possibly mix the wavefunctions from different steps before computing the Fock operator
       if(wfmixalg/=0 .and. .not. (wfmixalg==2 .and. abs(scf_history_wf%alpha-one)<tol8) )then
         call wf_mixing(atindx1,cg,cprj,dtset,istep_fock_outer,mcg,mcprj,mpi_enreg,&
&         nattyp,npwarr,pawtab,scf_history_wf)
         istep_fock_outer=istep_fock_outer+1

!DEBUG
         if(.false.)then
           !Update the density, from the newly mixed cg and cprj.
           !Be careful: in PAW, rho does not include the compensation density (added later) !
           tim_mkrho=6
           if (psps%usepaw==1) then
             ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
             ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
!          1-Compute density from WFs
             call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,&
&                       rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
             call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
!          2-Compute rhoij
             call pawmkrhoij(atindx,atindx1,cprj,dimcprj,dtset%istwfk,dtset%kptopt,dtset%mband,mband_cprj,&
&             mcprj,dtset%mkmem,mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
&             occ,dtset%paral_kgb,paw_dmft,pawrhoij,dtfil%unpaw,dtset%usewvl,dtset%wtk)
!          3-Symetrize rhoij, compute nhat and add it to rhor
!            Note pawrhoij_unsym and pawrhoij are the same, which means that pawrhoij
!            cannot be distributed over different atomic sites.
             cplex=1;ipert=0;idir=0;qpt(:)=zero
             call pawmkrho(1,compch_fft,cplex,gprimd,idir,indsym,ipert,mpi_enreg,&
&             my_natom,dtset%natom,dtset%nspden,dtset%nsym,dtset%ntypat,&
&             dtset%paral_kgb,pawang,pawfgr,pawfgrtab,dtset%pawprtvol,pawrhoij,pawrhoij,&
&             pawtab,qpt,rhowfg,rhowfr,rhor,rprimd,dtset%symafm,symrec,dtset%typat,ucvol,&
&             dtset%usewvl,xred,pawnhat=nhat,rhog=rhog)
!          2-Take care of kinetic energy density
             if(dtset%usekden==1)then
               call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,&
&                         rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
               call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,taug,rhowfr,taur)
             end if
             ABI_DEALLOCATE(rhowfg)
             ABI_DEALLOCATE(rhowfr)
           else
             write(std_out,*)' scfcv_core : recompute the density after the wf mixing '
             call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&             mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
             if(dtset%usekden==1)then
               call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,taug,taur,&
&                         rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
             end if
           end if
         end if
       end if
!ENDDEBUG

       ! Update data relative to the occupied states in fock
       call fock_updatecwaveocc(cg,cprj,dtset,fock,indsym,mcg,mcprj,mpi_enreg,nattyp,npwarr,occ,ucvol)
       ! Possibly (re)compute the ACE operator
       if(fock%fock_common%use_ACE/=0) then
         call fock2ACE(cg,cprj,fock,dtset%istwfk,kg,dtset%kptns,dtset%mband,mcg,mcprj,dtset%mgfft,&
&         dtset%mkmem,mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,dtset%natom,dtset%nband,dtset%nfft,ngfft,dtset%nkpt,dtset%nloalg,npwarr,dtset%nspden,&
&         dtset%nspinor,dtset%nsppol,dtset%ntypat,occ,dtset%optforces,paw_ij,pawtab,ph1d,psps,rprimd,&
&         dtset%typat,usecprj,dtset%use_gpu_cuda,dtset%wtk,xred,ylm)
         energies%e_fock0=fock%fock_common%e_fock0
       end if

       !Should place a test on whether there should be the final exit of the istep loop.
       !This test should use focktoldfe.
       !This should update the info in fock%fock_common%fock_converged.
       !For the time being, fock%fock_common%fock_converged=.false., so the loop end with the maximal value of nstep always,
       !except when nnsclo_hf==1 (so the Fock operator is always updated), in which case, the usual exit tests (toldfe, tolvrs, etc)
       !work fine.
       !if(fock%fock_common%nnsclo_hf==1 .and. fock%fock_common%use_ACE==0)then
       if(fock%fock_common%nnsclo_hf==1) fock%fock_common%fock_converged=.TRUE.

       !Depending on fockoptmix, possibly restart the mixing procedure for the potential
       if(mod(dtset%fockoptmix,10)==1) istep_mix=1
     else
       istep_updatedfock=istep_updatedfock+1
     end if

     !Used locally
     hyb_mixing=fock%fock_common%hyb_mixing ; hyb_mixing_sr=fock%fock_common%hyb_mixing_sr
   end if ! usefock

!  Initialize/update data in the electron-positron case
   if (dtset%positron<0.or.(dtset%positron>0.and.istep==1)) then
     call setup_positron(atindx,atindx1,cg,cprj,dtefield,dtfil,dtset,ecore,eigen,&
&     etotal,electronpositron,energies,fock,forces_needed,fred,gmet,gprimd,&
&     grchempottn,grcondft,grewtn,grvdw,gsqcut,hdr,initialized0,indsym,istep,istep_mix,kg,&
&     kxc,maxfor,mcg,mcprj,mgfftf,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfftf,ngrvdw,nhat,&
&     nkxc,npwarr,nvresid,occ,optres,paw_ij,pawang,pawfgr,pawfgrtab,&
&     pawrad,pawrhoij,pawtab,ph1df,ph1d,psps,rhog,rhor,rprimd,&
&     stress_needed,strsxc,symrec,ucvol,usecprj,vhartr,vpsp,vxc,vxctau,&
&     xccc3d,xcctau3d,xred,ylm,ylmgr)
     ipositron=electronpositron_calctype(electronpositron)
   end if

   if ((moved_atm_inside==1 .or. istep==1).or.&
&   (dtset%positron<0.and.istep_mix==1).or.&
&   (mod(dtset%fockoptmix,100)==11 .and. istep_updatedfock==1)) then
!    PAW only: we sometimes have to compute compensation density
!    and eventually add it to density from WFs
     nhatgrdim=0
     dummy_nhatgr = .False.
!    This is executed only in the positron case.
     if (psps%usepaw==1.and.(dtset%positron>=0.or.ipositron/=1) &
&     .and.((usexcnhat==0) &
&     .or.(dtset%xclevel==2.and.(dtfil%ireadwf/=0.or.dtfil%ireadden/=0.or.initialized/=0)) &
&     .or.(dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0))) then
       call timab(558,1,tsec)
       nhatgrdim=0;if (dtset%xclevel==2) nhatgrdim=usexcnhat*dtset%pawnhatxc
       ider=2*nhatgrdim;izero=0
       if (nhatgrdim>0)   then
         ABI_ALLOCATE(nhatgr,(cplex*nfftf,dtset%nspden,3*nhatgrdim))
       else
         ABI_ALLOCATE(nhatgr,(0,0,0))
         dummy_nhatgr = .True.
       end if
       call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
&       nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,&
&       nhatgr,nhat,pawrhoij,pawrhoij,pawtab,k0,rprimd,ucvol_local,dtset%usewvl,xred,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       comm_fft=spaceComm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&       distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
       if (dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0) then
         rhor(:,:)=rhor(:,:)+nhat(:,:)
         if(dtset%usewvl==0) then
           call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
         end if
       end if
       call timab(558,2,tsec)
     end if

!    The following steps have been gathered in the setvtr routine:
!    - get Ewald energy and Ewald forces
!    - compute local ionic pseudopotential vpsp
!    - possibly compute 3D core electron density xccc3d
!    - possibly compute 3D core kinetic energy density
!    - possibly compute vxc and vhartr
!    - set up vtrial

     optene = 4 * optres
     if(dtset%iscf==-3) optene=4
     if (wvlbigdft) optene = 1 ! VH needed for the WF mixing

     if (.not.allocated(nhatgr))  then
       ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3*nhatgrdim))
       dummy_nhatgr = .True.
     end if

     call setvtr(atindx1,dtset,energies,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcut,&
&     istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&     nattyp,nfftf,ngfftf,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,&
&     n1xccc,n3xccc,optene,pawrad,pawtab,ph1df,psps,rhog,rhor,rmet,rprimd,&
&     strsxc,ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,&
&     xccc3d,xred,electronpositron=electronpositron,&
&     taur=taur,vxc_hybcomp=vxc_hybcomp,vxctau=vxctau,add_tfw=tfw_activated,xcctau3d=xcctau3d)

     ! set the zero of the potentials here
     if(dtset%usepotzero==2) vpsp(:) = vpsp(:) + ecore / ( zion * ucvol )

     if(dtset%optdriver==RUNL_GWLS) call build_vxc(vxc,nfftf,dtset%nspden)

     if ((nhatgrdim>0.and.nstep>0).or.dummy_nhatgr) then
       ABI_DEALLOCATE(nhatgr)
     end if

!    Recursion Initialisation
     if(dtset%userec==1 .and. istep==1)  then
       rec_set%quitrec = 0
!      --At any step calculate the metric
       call Init_MetricRec(rec_set%inf,rec_set%nl%nlpsp,rmet,ucvol,rprimd,xred,dtset%ngfft(1:3),dtset%natom,rec_set%debug)
       call destroy_distribfft(rec_set%mpi%distribfft)
       call init_distribfft(rec_set%mpi%distribfft,'c',rec_set%mpi%nproc_fft,rec_set%ngfftrec(2),rec_set%ngfftrec(3))
       call init_distribfft(rec_set%mpi%distribfft,'f',rec_set%mpi%nproc_fft,dtset%ngfft(2),dtset%ngfft(3))
       if(initialized==0) call first_rec(dtset,psps,rec_set)
     end if

!    End the condition of atomic position change or istep==1
   end if
   call timab(240,2,tsec)
   call timab(241,1,tsec)

!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------
!  PAW: Compute energies and potentials in the augmentation regions (spheres)
!  Compute pseudopotential strengths (Dij quantities)
   if (psps%usepaw==1)then

!    Local exact exch.: impose occ. matrix if required
     if (dtset%useexexch/=0) then
       call setrhoijpbe0(dtset,initialized0,istep,istep_mix,&
&       spaceComm,my_natom,dtset%natom,dtset%ntypat,pawrhoij,pawtab,dtset%typat,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if

!    Computation of on-site densities/potentials/energies
     nzlmopt=0;if (istep_mix==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep_mix>2) nzlmopt=dtset%pawnzlm
     call paw_an_reset_flags(paw_an) ! Force the recomputation of on-site potentials
     call paw_ij_reset_flags(paw_ij,self_consistent=.true.) ! Force the recomputation of Dij
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,dtset%ixc,my_natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,&
&     pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&     electronpositron=electronpositron,vpotzero=vpotzero)

!    Correct the average potential with the calculated constant vpotzero
!    Correct the total energies accordingly
!    vpotzero(1) = -beta/ucvol
!    vpotzero(2) = -1/ucvol sum_ij rho_ij gamma_ij
     write(msg,'(a,f14.6,2x,f14.6)') &
&     ' average electrostatic smooth potential [Ha] , [eV]',SUM(vpotzero(:)),SUM(vpotzero(:))*Ha_eV
     call wrtout(std_out,msg,'COLL')
     vtrial(:,:)=vtrial(:,:)+SUM(vpotzero(:))
     if(option/=1)then
!      Fix the direct total energy (non-zero only for charged systems)
       energies%e_paw=energies%e_paw-SUM(vpotzero(:))*dtset%charge
!      Fix the double counting total energy accordingly (for both charged AND
!      neutral systems)
       energies%e_pawdc=energies%e_pawdc-SUM(vpotzero(:))*zion+vpotzero(2)*dtset%charge
     end if

!    PAW+U: impose density matrix if required
!           not available if usepawu<0 (PAW+U without occupation matrix)
     if (dtset%usepawu>0.and.(ipositron/=1)) then
       impose_dmat=0
       if ((istep<=abs(dtset%usedmatpu)).and.(dtset%usedmatpu<0.or.initialized0==0)) impose_dmat=1
       if (impose_dmat==1.or.dtset%dmatudiag/=0) then
         dimdmat=0;if (impose_dmat==1) dimdmat=2*lpawumax+1
         call setnoccmmp(0,dimdmat,&
&         dmatpawu(1:dimdmat,1:dimdmat,1:dtset%nsppol*dtset%nspinor,1:dtset%natpawu*impose_dmat),&
&         dtset%dmatudiag,impose_dmat,indsym,my_natom,dtset%natom,dtset%natpawu,&
&         dtset%nspinor,dtset%nsppol,dtset%nsym,dtset%ntypat,paw_ij,pawang,dtset%pawprtvol,&
&         pawrhoij,pawtab,dtset%spinat,dtset%symafm,dtset%typat,0,dtset%usepawu,&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!        Reinitalize mixing if PAW+U and occupation matrix now allowed to change
!        For experimental purpose...
         if ((dtset%userib==1234).and.(istep==abs(dtset%usedmatpu)).and. &
&         (dtset%usedmatpu<0.or.initialized0==0)) reset_mixing=.true.
       end if
     end if

!    Dij computation
     call timab(561,1,tsec)

     call pawdij(cplex,dtset%enunit,gprimd,ipert,my_natom,dtset%natom,nfftf,nfftotf,&
&     dtset%nspden,psps%ntypat,paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,&
&     pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,k0,dtset%spnorbscl,&
&     ucvol_local,dtset%charge,vtrial,vxc,xred,&
&     natvshift=dtset%natvshift,atvshift=dtset%atvshift,fatvshift=fatvshift,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     mpi_comm_grid=spaceComm_grid,&
&     hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,&
&     electronpositron_calctype=ipositron,&
&     electronpositron_pawrhoij=pawrhoij_ep,&
&     electronpositron_lmselect=lmselect_ep,&
&     nucdipmom=dtset%nucdipmom)

!    Symetrize Dij
     call symdij(gprimd,indsym,ipert,my_natom,dtset%natom,dtset%nsym,&
&     psps%ntypat,0,paw_ij,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     if (has_dijhat==1) then
       call symdij(gprimd,indsym,ipert,my_natom,dtset%natom,dtset%nsym,&
&       psps%ntypat,1,paw_ij,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if
     if(dtset%usewvl==1) then
       call paw2wvl_ij(3,paw_ij,wvl%descr)
     end if

     call timab(561,2,tsec)
   end if

!  Write out occupancies to dtpawuj-dataset
   if (dtset%usepawu/=0.and.dtset%macro_uj>0.and.istep>1.and.ipositron/=1) then
     call pawuj_red(dtset,dtpawuj,fatvshift,my_natom,dtset%natom,dtset%ntypat,&
     paw_ij,pawrad,pawtab,ndtpawuj,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if

   call timab(241,2,tsec)

!  No need to continue and call vtorho, when nstep==0
   if(nstep==0)exit

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   call timab(56,1,tsec)

   if(dtset%iscf>=0)then
     write(msg, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,msg,'COLL')
   end if

!  The next flag says whether the xred have to be changed in the current iteration
   moved_atm_inside=0
   ! /< Hack to remove iapp from scfcv_core
   ! for ionmov 4|5 ncycle=1
   ! Hence iapp = itime
   if ( dtset%jdtset /= scfcv_jdtset ) then
     ! new dtset -> reinitialize
     scfcv_jdtset = dtset%jdtset
     scfcv_itime = 0
   end if
   if ( istep .eq. 1 ) scfcv_itime = scfcv_itime + 1
   if(dtset%ionmov==4 .and. mod(scfcv_itime,2)/=1 .and. dtset%iscf>=0 ) moved_atm_inside=1
   if(dtset%ionmov==5 .and. scfcv_itime/=1 .and. istep==1 .and. dtset%iscf>=0) moved_atm_inside=1
   ! /< Hack to remove iapp from scfcv_core

!  Thomas-Fermi scheme might use a different toldfe criterion
   if (dtset%tfkinfunc>0.and.dtset%tfkinfunc/=2) then
     tollist(4)=dtset%toldfe;if (.not.tfw_activated) tollist(4)=dtset%tfw_toldfe
   end if

!  The next flag says whether the forces have to be computed in the current iteration
   computed_forces=0
   if ((dtset%optforces==1 .and. dtset%usewvl == 0).or.(moved_atm_inside==1)) computed_forces=1
   if (abs(tollist(3))>tiny(0._dp)) computed_forces=1
   if (dtset%iscf<0) computed_forces=0
   if ((istep==1).and.(dtset%optforces/=1)) then
     if (moved_atm_inside==1) then
       write(msg,'(5a)')&
&       'Although the computation of forces during electronic iterations',ch10,&
&       'was not required by user, it is done (required by the',ch10,&
&       'choice of ionmov input parameter).'
       MSG_WARNING(msg)
     end if
     if (abs(tollist(3))+abs(tollist(7))>tiny(0._dp)) then
       write(msg,'(5a)')&
&       'Although the computation of forces during electronic iterations',ch10,&
&       'was not required by user, it is done (required by the',ch10,&
&       '"toldff" or "tolrff" tolerance criteria).'
       MSG_WARNING(msg)
     end if
   end if
   if ((istep==1).and.(dtset%optforces==1).and. dtset%usewvl == 1) then
     write(msg,'(5a)')&
&     'Although the computation of forces during electronic iterations',ch10,&
&     'was required by user, it has been disable since the tolerence',ch10,&
&     'is not on forces (force computation is expensive in wavelets).'
     MSG_WARNING(msg)
   end if

   call timab(56,2,tsec)

!  ######################################################################
!  Compute the density rho from the trial potential
!  ----------------------------------------------------------------------
   call timab(242,1,tsec)
!  Compute the density from the trial potential
   if (dtset%tfkinfunc==0) then
     if(VERBOSE) call wrtout(std_out,'*. Compute the density from the trial potential (vtorho)',"COLL")

     call vtorho(itime,afford,atindx,atindx1,cg,compch_fft,cprj,cpus,dbl_nnsclo,&
&     dielop,dielstrt,dmatpawu,dphase,dtefield,dtfil,dtset,&
&     eigen,electronpositron,energies,etotal,gbound_diel,&
&     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&     istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mcprj,mgfftdiel,mpi_enreg,&
&     my_natom,dtset%natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&     npwarr,npwdiel,res2,psps%ntypat,nvresid,occ,&
&     computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
&     pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,fock,&
&     pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,rmet,rprimd,&
&     susmat,symrec,taug,taur,nvtauresid,ucvol_local,usecprj,wffnew,with_vectornd,&
&     vectornd,vtrial,vxctau,wvl,xred,ylm,ylmgr,ylmdiel, rmm_diis_status)

   else if (dtset%tfkinfunc==1.or.dtset%tfkinfunc==11.or.dtset%tfkinfunc==12) then
     MSG_WARNING('THOMAS FERMI')
     call vtorhotf(dtset,energies%e_kinetic,energies%e_nlpsp_vfock,&
&     energies%entropy,energies%e_fermie,gprimd,grnl,irrzon,mpi_enreg,&
&     dtset%natom,nfftf,dtset%nspden,dtset%nsppol,dtset%nsym,phnons,&
&     rhog,rhor,rprimd,ucvol,vtrial)

     residm=zero
     energies%e_eigenvalues=zero
   end if

!  Recursion method
   if(dtset%userec==1)then
     call vtorhorec(dtset,&
&     energies%e_kinetic,energies%e_nlpsp_vfock,energies%entropy,energies%e_eigenvalues,&
&     energies%e_fermie,grnl,initialized,irrzon,nfftf,phnons,&
&     rhog,rhor,vtrial,rec_set,istep-nstep,rprimd,gprimd)
     residm=zero
   end if

   ! Update Fermi level in energies
   results_gs%fermie = energies%e_fermie

   if(dtset%wfoptalg==2)then
     do ikpt=1,dtset%nkpt
       shiftvector(1+(ikpt-1)*(dtset%mband+2))=val_min
       shiftvector(2+(ikpt-1)*(dtset%mband+2):ikpt*(dtset%mband+2)-1)=&
&       eigen((ikpt-1)*dtset%mband+1:ikpt*dtset%mband)
       shiftvector(ikpt*(dtset%mband+2))=val_max
     end do
   end if

   call timab(242,2,tsec)

!  ######################################################################
!  Skip out of step loop if non-SCF (completed)
!  ----------------------------------------------------------------------

!  Indeed, nstep loops have been done inside vtorho
   if (dtset%iscf<0) exit

!  ######################################################################
!  In case of density mixing or wavelet handling, compute the total energy
!  ----------------------------------------------------------------------
   call timab(60,1,tsec)
   if (dtset%iscf>=10 .or. wvlbigdft) then
     optene = 1  ! use double counting scheme (default)
     if (wvlbigdft.and.dtset%iscf==0) optene = 0 ! use direct scheme
     if (dtset%iscf==22) optene = -1

!    Add the Fock contribution to E_xc and E_xcdc if required
     if (usefock==1) then
       energies%e_fockdc=two*energies%e_fock
     end if

!    if the mixing is the ODA mixing, compute energy and new density here
     if (dtset%iscf==22) then
       call odamix(deltae,dtset,&
&       elast,energies,etotal,gprimd,gsqcut,kxc,mpi_enreg,&
&       my_natom,nfftf,ngfftf,nhat,nkxc,psps%ntypat,nvresid,n3xccc,optres,&
&       paw_ij,paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&       red_ptot,psps,rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,&
&       usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred,&
&       taur=taur,vxctau=vxctau,add_tfw=tfw_activated)
     end if
!    If the density mixing is required, compute the total energy here
!    TODO: add nvtauresid if needed (for forces?)
     call etotfor(atindx1,deltae,diffor,dtefield,dtset,&
&     elast,electronpositron,energies,&
&     etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,grcondft,gresid,grewtn,grhf,grnl,grvdw,&
&     grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,my_natom,&
&     nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,psps%ntypat,nvresid,n1xccc,n3xccc,&
&     optene,computed_forces,optres,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&     ph1df,red_ptot,psps,rhog,rhor,rmet,rprimd,symrec,synlgr,ucvol,&
&     psps%usepaw,vhartr,vpsp,vxc,vxctau,wvl%descr,wvl%den,xccc3d,xred)
    !if (wvlbigdft) energies_copy(energies,energies_wvl) ! TO BE ACTIVATED LATER
   end if
   call timab(60,2,tsec)

!  ######################################################################
!  In case of density mixing, check the exit criterion
!  ----------------------------------------------------------------------
   if (dtset%iscf>=10.or.(wvlbigdft.and.dtset%iscf>0)) then
!    Check exit criteria
     call timab(52,1,tsec)
     choice=2
     if(paw_dmft%use_dmft==1) then
       call prtene(dtset,energies,std_out,psps%usepaw)
     end if
     call scprqt(choice,cpus,deltae,diffor,dtset,&
&     eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
&     dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,istep_fock_outer,istep_mix,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode,&
&     electronpositron=electronpositron,fock=fock)
     call timab(52,2,tsec)

!    Check if we need to exit the loop
     call timab(244,1,tsec)
     if (dtset%tfkinfunc>10.and.(.not.tfw_activated).and.quit==1) then
       quit=0;tfw_activated=.true.;reset_mixing=.true.
     end if
     if(dtset%userec==1.and.rec_set%quitrec==2)quit=1
     if (istep==nstep) quit=1
     quit_sum=quit
     call xmpi_sum(quit_sum,spaceComm,ierr)
     if (quit_sum>0) quit=1
     call timab(244,2,tsec)

!    If criteria in scprqt say to quit, then exit the loop over istep.
     if (quit==1) exit
   end if

!  ######################################################################
!  Mix the total density (if required)
!  ----------------------------------------------------------------------
   call timab(68,1,tsec)

   if (dtset%iscf>=10 .and.dtset%iscf/=22.and. .not. wvlbigdft ) then

!    If LDA dielectric matrix is used for preconditionning, has to update here Kxc
     if (nkxc>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79) &
&     .and.((istep==1.or.istep==dielstrt).or.(dtset%iprcel>=100))) then
       optxc=10
       call xcdata_init(xcdata,dtset=dtset)
!      to be adjusted for the call to rhotoxc
       nk3xc=1
       if(dtset%icoulomb==0 .and. dtset%usewvl==0) then
         non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
         call rhotoxc(edum,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,0,nkxc,nk3xc,non_magnetic_xc,n3xccc,&
&         optxc,rhor,rprimd,dummy2,0,vxc,vxcavg_dum,xccc3d,xcdata,&
&         add_tfw=tfw_activated,taur=taur,vhartr=vhartr,vxctau=vxctau,xcctau3d=xcctau3d)
       else if(.not. wvlbigdft) then
!        WVL case:
         call psolver_rhohxc(energies%e_hartree, energies%e_xc, evxc, &
&         dtset%icoulomb, dtset%ixc, &
&         mpi_enreg, nfftf, ngfftf,&
&         nhat,psps%usepaw,&
&         dtset%nscforder,dtset%nspden,n3xccc,rhor,rprimd, &
&         usexcnhat,psps%usepaw,dtset%usewvl,vhartr, vxc, vxcavg,&
&         wvl%descr,wvl%den,&
&         wvl%e,xccc3d,dtset%xclevel,dtset%xc_denpos)
       end if
     end if

     call newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,&
&     dtset,etotal,fcart,pawfgr%fintocoa,&
&     gmet,grhf,gsqcut,initialized,ispmix,istep_mix,kg_diel,kxc,&
&     mgfftf,mix,pawfgr%coatofin,moved_atm_inside,mpi_enreg,my_natom,nattyp,nfftf,&
&     nfftmix,nfftmix_per_nfft,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,nvresid,psps%ntypat,&
&     n1xccc,pawrhoij,pawtab,ph1df,psps,rhog,rhor,&
&     rprimd,susmat,psps%usepaw,vtrial,wvl%descr,wvl%den,xred,&
&     mix_mgga=mix_mgga,taug=taug,taur=taur,tauresid=nvtauresid)
   end if   ! iscf>=10

   call timab(68,2,tsec)

!  ######################################################################
!  Additional computation in case of an electric field or electric displacement field
!  ----------------------------------------------------------------------

   call timab(239,1,tsec)

   call update_e_field_vars(atindx,atindx1,cg,dimcprj,dtefield,dtfil,dtset,&
&   efield_old_cart,gmet,gprimd,hdr,idir,kg,mcg,&
&   dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,dtset%natom,nattyp,ngfft,dtset%nkpt,npwarr,&
&   dtset%ntypat,pawrhoij,pawtab,pel_cg,pelev,pion,psps,ptot,ptot_cart,&
&   pwind,pwind_alloc,pwnsfac,red_efield2,red_efield2_old,red_ptot,rmet,rprimd,&
&   1,quit,istep,ucvol,unit_out,psps%usepaw,xred,ylm,ylmgr)

   call timab(239,2,tsec)

!  ######################################################################
!  Compute the new potential from the trial density
!  ----------------------------------------------------------------------

   call timab(243,1,tsec)
   if(VERBOSE) call wrtout(std_out,'*. Compute the new potential from the trial density',"COLL")

!  Set XC computation flag
   optxc=1
   if (nkxc>0) then
! MJV 2017 May 25: you should not be able to get here with iscf < 0
     if (dtset%iscf<0) optxc=2
     if (modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79).and. &
&     dtset%iscf<10.and. &
&     (dtset%iprcel>=100.or.istep==1.or.istep==dielstrt)) optxc=2
     if (dtset%iscf>=10.and.dtset%densfor_pred/=0.and.abs(dtset%densfor_pred)/=5) optxc=2
     if (optxc==2.and.dtset%xclevel==2.and.nkxc==2*min(dtset%nspden,2)-1) optxc=12
   end if

   if (dtset%iscf/=22) then
!    PAW: eventually recompute compensation density (and gradients)
     nhatgrdim=0
     if ( allocated(nhatgr) ) then
       ABI_DEALLOCATE(nhatgr)
     end if
     if (psps%usepaw==1) then
       ider=-1;if (dtset%iscf>=10.and.((dtset%xclevel==2.and.dtset%pawnhatxc>0).or.usexcnhat==0)) ider=0
       if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) ider=ider+2
       if (ipositron==1) ider=-1
       if (ider>0) then
         nhatgrdim=1
         ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3))
       else
         ABI_ALLOCATE(nhatgr,(0,0,0))
       end if
       if (ider>=0) then
         call timab(558,1,tsec)
         izero=0

         call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,nfftf,ngfftf,&
&         nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,nhatgr,nhat,&
&         pawrhoij,pawrhoij,pawtab,k0,rprimd,ucvol_local,dtset%usewvl,xred,&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&         comm_fft=spaceComm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&         distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)

         call timab(558,2,tsec)
       end if
     else
       ABI_ALLOCATE(nhatgr,(0,0,0))
     end if

!    Compute new potential from the trial density

     optene=2*optres;if(psps%usepaw==1) optene=2
     call rhotov(constrained_dft,dtset,energies,gprimd,grcondft,gsqcut,intgres,istep,kxc,mpi_enreg,nfftf,ngfftf, &
&     nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,optene,optres,optxc,&
&     rhog,rhor,rprimd,strsxc,ucvol_local,psps%usepaw,usexcnhat,&
&     vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,wvl,xccc3d,xred,&
&     electronpositron=electronpositron,taur=taur,vxctau=vxctau,vtauresid=nvtauresid,&
&     vxc_hybcomp=vxc_hybcomp,add_tfw=tfw_activated,xcctau3d=xcctau3d)

   end if

   call timab(243,2,tsec)
   call timab(60,1,tsec)

!  This is inside the loop, its not equivalent to the line 1821
   if(moved_atm_inside==1) xred_old(:,:)=xred(:,:)

   if (dtset%iscf<10) then

     if(VERBOSE) call wrtout(std_out,'Check exit criteria in case of potential mixing',"COLL")

!    If the potential mixing is required, compute the total energy here
!    PAW: has to compute here spherical terms
     if (psps%usepaw==1) then
       nzlmopt=0;if (istep_mix==1.and.dtset%pawnzlm>0) nzlmopt=-1
       if (istep_mix>1) nzlmopt=dtset%pawnzlm
       call paw_an_reset_flags(paw_an) ! Force the recomputation of on-site potentials
       option=2
       call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,&
&       dtset%ixc,my_natom,dtset%natom,dtset%nspden,&
&       psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,&
&       paw_ij,pawang,dtset%pawprtvol,pawrad,pawrhoij,dtset%pawspnorb,&
&       pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&       hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       electronpositron=electronpositron)

     end if

!    Add the Fock contribution to E_xc and E_xcdc if required
     if (usefock==1) energies%e_fockdc=two*energies%e_fock

     if (.not.wvlbigdft) then
! TODO: add nvtauresid if needed (for forces?)
       call etotfor(atindx1,deltae,diffor,dtefield,dtset,&
&       elast,electronpositron,energies,&
&       etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,grcondft,gresid,grewtn,grhf,grnl,grvdw,&
&       grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,my_natom,&
&       nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,dtset%ntypat,nvresid,n1xccc, &
&       n3xccc,0,computed_forces,optres,pawang,pawfgrtab,pawrad,pawrhoij,&
&       pawtab,ph1df,red_ptot,psps,rhog,rhor,rmet,rprimd,symrec,synlgr,ucvol,&
&       psps%usepaw,vhartr,vpsp,vxc,vxctau,wvl%descr,wvl%den,xccc3d,xred)
     end if

   end if
   call timab(60,2,tsec)

!  ######################################################################
!  Check exit criteria in case of potential mixing or direct minimization
!  ----------------------------------------------------------------------
   if ((dtset%iscf<10.and.(.not.wvlbigdft)) .or. dtset%iscf == 0) then
!    Check exit criteria
     call timab(52,1,tsec)
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,&
&     eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
&     dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,istep_fock_outer,istep_mix,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode,&
&     electronpositron=electronpositron,fock=fock)
     call timab(52,2,tsec)

!    Check if we need to exit the loop
     call timab(244,1,tsec)
     if (dtset%tfkinfunc>10.and.(.not.tfw_activated).and.quit==1) then
       quit=0;tfw_activated=.true.;reset_mixing=.true.
     end if
     if (istep==nstep.and.psps%usepaw==1) quit=1
     if(dtset%userec==1 .and. rec_set%quitrec==2) quit=1
     quit_sum=quit
     call xmpi_sum(quit_sum,spaceComm,ierr)
     if (quit_sum > 0) quit=1

!    If criteria in scprqt say to quit, then exit the loop over istep.
     if (quit==1) then
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+nvresid(:,ispden)+vres_mean(ispden)
       end do
       !if (dtset%usekden==1) then
       !  do ispden=1,dtset%nspden
       !    vxctau(:,ispden,1)=vxctau(:,ispden,1)+nvtauresid(:,ispden)
       !  end do
       !end if
       call timab(244,2,tsec) ! Due to the exit instruction, two timab calls are needed
       exit ! exit the loop over istep
     end if
     call timab(244,2,tsec) ! Due to the exit instruction, two timab calls are needed
   end if

!  ######################################################################
!  Mix the potential (if required) - Check exit criteria
!  ----------------------------------------------------------------------

   call timab(245,1,tsec)
   if (dtset%iscf<10 .and. dtset%iscf>0 .and. .not. wvlbigdft) then

     if(VERBOSE) call wrtout(std_out,'*. Mix the potential (if required) - Check exit criteria',"COLL")

!    Precondition the residual and forces, then determine the new vtrial
!    (Warning: the (H)xc potential may have been subtracted from vtrial)

     call newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
&     dtn_pc,dtset,etotal,fcart,pawfgr%fintocoa,&
&     gmet,grhf,gsqcut,initialized,ispmix,&
&     istep_mix,kg_diel,kxc,mgfftf,mix,pawfgr%coatofin,&
&     moved_atm_inside,mpi_enreg,my_natom,nattyp,nfftf,nfftmix,&
&     ngfftf,ngfftmix,nkxc,npawmix,npwdiel,&
&     nstep,psps%ntypat,n1xccc,&
&     pawrhoij,ph1df,psps,rhor,rprimd,susmat,psps%usepaw,&
&     vhartr,vnew_mean,vpsp,nvresid,vres_mean,vtrial,vxc,xred,&
&     nfftf,pawtab,rhog,wvl,&
&     mix_mgga=mix_mgga,vtau=vxctau,vtauresid=nvtauresid)

   end if   ! iscf<10

!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  ######################################################################

   if(VERBOSE) call wrtout(std_out,'*. END MINIMIZATION ITERATIONS',"COLL")

!  The initialisation of the gstate run should be done when this point is reached
   initialized=1

!  This is to save the density for restart.
   if (iwrite_fftdatar(mpi_enreg)) then

     if(dtset%prtden<0.or.dtset%prtkden<0) then
!      Update the content of the header (evolving variables)
!      Don't use parallelism over atoms because only me=0 accesses here
       bantot=hdr%bantot
       if (dtset%positron==0) then
         call hdr%update(bantot,etotal,energies%e_fermie,residm,&
&         rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1))
       else
         call hdr%update(bantot,electronpositron%e0,energies%e_fermie,residm,&
&         rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1))
       end if
     end if

     if (dtset%prtden<0) then
       if (mod(istep-1,abs(dtset%prtden))==0) then
         isave_den=isave_den+1
         rdwrpaw=0
         call int2char4(mod(isave_den,2),tag)
         ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
         fildata=trim(dtfil%fnametmp_app_den)//'_'//trim(tag)
         if (dtset%iomode == IO_MODE_ETSF) fildata = nctk_ncify(fildata)
         call fftdatar_write_from_hdr("density",fildata,dtset%iomode,hdr,ngfftf,cplex1,nfftf,&
&         dtset%nspden,rhor,mpi_enreg,eigen=eigen)
       end if
     end if

     if (dtset%prtkden<0) then
       if (mod(istep-1,abs(dtset%prtkden))==0) then
         isave_kden=isave_kden+1
         rdwrpaw=0
         call int2char4(mod(isave_kden,2),tag)
         ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
         fildata=trim(dtfil%fnametmp_app_kden)//'_'//trim(tag)
         if (dtset%iomode == IO_MODE_ETSF) fildata = nctk_ncify(fildata)
         ! output the Laplacian of density
         call fftdatar_write_from_hdr("kinedr",fildata,dtset%iomode,hdr,ngfftf,cplex1,nfftf,&
&         dtset%nspden,taur,mpi_enreg,eigen=eigen)
       end if
     end if

   end if

   ABI_DEALLOCATE(nhatgr)

   istep_mix=istep_mix+1
   if (reset_mixing) then
     istep_mix=1;reset_mixing=.false.
   end if
   if (ipositron/=0) electronpositron%istep_scf=electronpositron%istep_scf+1

   call timab(245,2,tsec)
 end do ! istep

 ABI_FREE(rmm_diis_status)
 ABI_SFREE(nhatgr)

 ! Avoid pending requests if itime == ntime.
 call xmpi_wait(quitsum_request,ierr)
 if (timelimit_exit == 1) istep = istep - 1

 call timab(246,1,tsec)

 if (dtset%iscf > 0) then
   call ab7_mixing_deallocate(mix)
   if (dtset%usekden/=0) call ab7_mixing_deallocate(mix_mgga)
 end if

 if (usefock==1)then
   if(wfmixalg/=0) call scf_history_free(scf_history_wf)
 end if

 if (quit==1.and.nstep==1) initialized=1

!######################################################################
!Case nstep==0: compute energy based on incoming wf
!----------------------------------------------------------------------

 if(nstep==0) then
   optene=2*psps%usepaw+optres
   energies%entropy=results_gs%energies%entropy  !MT20070219: entropy is not recomputed in routine energy
   if (.not.allocated(nhatgr) ) then
     ABI_ALLOCATE(nhatgr,(0,0,0))
   end if

   call energy(cg,compch_fft,constrained_dft,dtset,electronpositron,&
&   energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,my_natom,&
&   nfftf,ngfftf,nhat,nhatgr,nhatgrdim,npwarr,n3xccc,&
&   occ,optene,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
&   phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,taug,taur,usexcnhat,&
&   vhartr,vtrial,vpsp,vxc,wvl%wfs,wvl%descr,wvl%den,wvl%e,xccc3d,xred,ylm,&
&   add_tfw=tfw_activated,vxctau=vxctau)

   if (nhatgrdim>0)  then
     ABI_DEALLOCATE(nhatgr)
   end if

 end if ! nstep==0

!######################################################################
!Additional steps after SC iterations, including force, stress, polarization calculation
!----------------------------------------------------------------------

 if (dtset%userec==1) then
   call prtene(dtset,energies,ab_out,psps%usepaw)
   call prtene(dtset,energies,std_out,psps%usepaw)
 end if

!if (wvlbigdft) call energies_copy(energies_wvl,energies) ! TO BE ACTIVATED LATER

!PAW: if cprj=<p_lmn|Cnk> are in memory,
!need to reorder them (from atom-sorted to unsorted)
 if (psps%usepaw==1.and.usecprj==1) then
   iorder_cprj=1
   call pawcprj_reorder(cprj,atindx1)
   if (dtset%positron/=0) then
     if (electronpositron%dimcprj>0) then
       call pawcprj_reorder(electronpositron%cprj_ep,atindx1)
     end if
   end if
   if (dtset%usewvl==1) then
     call wvl_cprjreorder(wvl%descr,atindx1)
   end if
 end if

!PAW: if cprj=<p_lmn|Cnk> are not in memory,need to compute them in some cases
 recompute_cprj = psps%usepaw ==1 .and. usecprj==0 .and. &
& (dtset%prtwant  ==2 .or. &
& dtset%prtwant  ==3  .or. &
& dtset%prtnabla > 0  .or. &
& dtset%prtdos   ==3  .or. &
& dtset%berryopt /=0  .or. &
& dtset%orbmag /=0    .or. &
& dtset%kssform  ==3  .or. &
& dtset%pawfatbnd> 0  .or. &
& dtset%pawprtwf > 0  )

 if(dtset%orbmag .NE. 0) recompute_cprj=.TRUE.
 if (recompute_cprj) then
   usecprj=1
   mband_cprj=dtset%mband/mpi_enreg%nproc_band
   mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
   ABI_DATATYPE_ALLOCATE(cprj_local,(dtset%natom,mcprj))
   ncpgr = 0 ; ctocprj_choice = 1
   if (finite_efield_flag) then
     if (forces_needed /= 0 .and. stress_needed == 0) then
       ncpgr = 3 ; ctocprj_choice = 2
     else if (forces_needed /= 0 .and. stress_needed /= 0) then
       ncpgr = 9 ; ctocprj_choice = 23
     else if (forces_needed == 0 .and. stress_needed /= 0) then
       ncpgr = 6 ; ctocprj_choice = 3
     end if
   end if
   if ((dtset%orbmag .GT. 1) .OR. (dtset%orbmag .LT. -1) ) then
      ncpgr=3; ctocprj_choice=5; idir=0
   end if
   call pawcprj_alloc(cprj_local,ncpgr,dimcprj)
   cprj=> cprj_local
   iatom=0 ; iorder_cprj=1 ! cprj are not ordered
   call ctocprj(atindx,cg,ctocprj_choice,cprj_local,gmet,gprimd,&
&   iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&   mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&   dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,&
&   dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&   ucvol,dtfil%unpaw,xred,ylm,ylmgr)
 end if

 call timab(246,2,tsec)
 call timab(247,1,tsec)

!SHOULD CLEAN THE ARGS OF THIS ROUTINE
 call afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&
& deltae,diffor,dtefield,dtfil,dtorbmag,dtset,eigen,electronpositron,elfr,&
& energies,etotal,favg,fcart,fock,forold,fred,grchempottn,grcondft,&
& gresid,grewtn,grhf,grhor,grvdw,&
& grxc,gsqcut,hdr,indsym,intgres,irrzon,istep,istep_fock_outer,istep_mix,&
& kg,kxc,lrhor,maxfor,mcg,mcprj,mgfftf,&
& moved_atm_inside,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfft,ngfftf,ngrvdw,nhat,&
& nkxc,npwarr,nvresid,occ,optres,paw_an,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,ph1d,ph1df,phnons,pion,prtfor,&
& prtxml,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,&
& taur,tollist,usecprj,vectornd,vhartr,vpsp,vtrial,vxc,vxctau,vxcavg,with_vectornd,wvl,&
& xccc3d,xcctau3d,xred,ylm,ylmgr,dtset%charge*SUM(vpotzero(:)),conv_retcode)

!Before leaving the present routine, save the current value of xred.
 xred_old(:,:)=xred(:,:)

 call timab(247,2,tsec)

!######################################################################
!All calculations in scfcv_core are finished. Printing section
!----------------------------------------------------------------------

 call timab(248,1,tsec)

 call outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dmatpawu,dtfil,&
& dtset,ecut,eigen,electronpositron,elfr,etotal,&
& gmet,gprimd,grhor,hdr,intgres,kg,lrhor,dtset%mband,mcg,mcprj,dtset%mgfft,&
& dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,my_natom,dtset%natom,nattyp,&
& nfftf,ngfftf,nhat,dtset%nkpt,npwarr,dtset%nspden,&
& dtset%nsppol,dtset%nsym,psps%ntypat,n3xccc,occ,paw_dmft,pawang,pawfgr,pawfgrtab,&
& pawrad,pawrhoij,pawtab,paw_an,paw_ij,dtset%prtvol,psps,results_gs,&
& rhor,rprimd,taur,ucvol,usecprj,vhartr,vpsp,vtrial,vxc,wvl%den,xccc3d,xred)

 call timab(248,2,tsec)
 call timab(249,1,tsec)

!Transfer eigenvalues and occupation computed by BigDFT in afterscfloop to eigen.
#if defined HAVE_BIGDFT
 if (dtset%usewvl == 1) then
   if (dtset%nsppol == 1) then
     eigen = wvl%wfs%ks%orbs%eval
     occ = wvl%wfs%ks%orbs%occup
   else
     eigen(1:wvl%wfs%ks%orbs%norbu) = wvl%wfs%ks%orbs%eval(1:wvl%wfs%ks%orbs%norbu)
     eigen(dtset%mband + 1:dtset%mband + wvl%wfs%ks%orbs%norbd) = &
&     wvl%wfs%ks%orbs%eval(wvl%wfs%ks%orbs%norbu + 1:wvl%wfs%ks%orbs%norb)
     occ(1:wvl%wfs%ks%orbs%norbu) = wvl%wfs%ks%orbs%occup(1:wvl%wfs%ks%orbs%norbu)
     occ(dtset%mband + 1:dtset%mband + wvl%wfs%ks%orbs%norbd) = &
&     wvl%wfs%ks%orbs%occup(wvl%wfs%ks%orbs%norbu + 1:wvl%wfs%ks%orbs%norb)
   end if
 end if
#endif
!need to reorder cprj (from unsorted to atom-sorted)
 if (psps%usepaw==1.and.usecprj==1) then
   iorder_cprj=0
   call pawcprj_reorder(cprj,atindx)
   if (dtset%positron/=0) then
     if (electronpositron%dimcprj>0) then
       call pawcprj_reorder(electronpositron%cprj_ep,atindx)
     end if
   end if
 end if
!######################################################################
!Deallocate memory and save results
!----------------------------------------------------------------------

 call prc_mem_free()

 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(fred)
 ABI_DEALLOCATE(forold)
 ABI_DEALLOCATE(grchempottn)
 ABI_DEALLOCATE(grcondft)
 ABI_DEALLOCATE(gresid)
 ABI_DEALLOCATE(grewtn)
 ABI_DEALLOCATE(grnl)
 ABI_DEALLOCATE(grvdw)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(intgres)
 ABI_DEALLOCATE(synlgr)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(vhartr)
 ABI_DEALLOCATE(vtrial)
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(vxc)
 ABI_DEALLOCATE(vxc_hybcomp)
 ABI_DEALLOCATE(vxctau)
 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(shiftvector)
 ABI_DEALLOCATE(dtn_pc)
 ABI_DEALLOCATE(grhf)
 ABI_DEALLOCATE(nvresid)
 ABI_DEALLOCATE(nvtauresid)

 if(allocated(vectornd)) then
    ABI_DEALLOCATE(vectornd)
 end if

 if((nstep>0.and.dtset%iscf>0).or.dtset%iscf==-1) then
   ABI_DEALLOCATE(dielinv)
 end if
 ABI_DEALLOCATE(gbound_diel)
 ABI_DEALLOCATE(irrzondiel)
 ABI_DEALLOCATE(kg_diel)
 ABI_DEALLOCATE(phnonsdiel)
 ABI_DEALLOCATE(susmat)
 ABI_DEALLOCATE(ph1ddiel)
 ABI_DEALLOCATE(ylmdiel)

 if (psps%usepaw==1) then
   if (dtset%iscf>0) then
     do iatom=1,my_natom
       pawrhoij(iatom)%lmnmix_sz=0
       pawrhoij(iatom)%use_rhoijres=0
       ABI_DEALLOCATE(pawrhoij(iatom)%kpawmix)
       ABI_DEALLOCATE(pawrhoij(iatom)%rhoijres)
     end do
   end if
!   if (recompute_cprj.or.usecprj==1) then
   if (recompute_cprj) then
     usecprj=0;mcprj=0
     call pawcprj_free(cprj)
     ABI_DATATYPE_DEALLOCATE(cprj_local)
   end if
   call paw_an_free(paw_an)
   call paw_ij_free(paw_ij)
   call pawfgrtab_free(pawfgrtab)
   if(dtset%usewvl==1) then
#if defined HAVE_BIGDFT
     call cprj_clean(wvl%descr%paw%cprj)
     ABI_DATATYPE_DEALLOCATE(wvl%descr%paw%cprj)
#endif
     call paw2wvl_ij(2,paw_ij,wvl%descr)
   end if
 end if
 ABI_DATATYPE_DEALLOCATE(pawfgrtab)
 ABI_DATATYPE_DEALLOCATE(paw_an)
 ABI_DATATYPE_DEALLOCATE(paw_ij)
 ABI_DEALLOCATE(nhat)
 ABI_DEALLOCATE(xcctau3d)
 ABI_DEALLOCATE(dimcprj_srt)
 ABI_DEALLOCATE(dimcprj)


! Deallocate exact exchange data at the end of the calculation
 if (usefock==1) then
   if (fock%fock_common%use_ACE/=0) call fock_ACE_destroy(fock%fockACE)
   call fock_common_destroy(fock%fock_common)
   call fock_BZ_destroy(fock%fock_BZ)
   call fock_destroy(fock)
   nullify(fock)
 end if

 if (prtxml == 1) then
!  We output the final result given in results_gs
   write(ab_xml_out, "(A)") '      <finalConditions>'
   call out_resultsgs_XML(dtset, 4, results_gs, psps%usepaw)
   write(ab_xml_out, "(A)") '      </finalConditions>'
   write(ab_xml_out, "(A)") '    </scfcvLoop>'
 end if

!Free the datastructure constrained_dft
 call constrained_dft_free(constrained_dft)

 call timab(249,2,tsec)
 call timab(238,2,tsec)

 DBG_EXIT("COLL")

end subroutine scfcv_core
!!***

!!****f* ABINIT/etotfor
!! NAME
!! etotfor
!!
!! FUNCTION
!! This routine is called to compute the total energy and various parts of it.
!! The routine computes -if requested- the forces.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   | bfield = cartesian coordinates of magnetic field in atomic units
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along some direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | densfor_pred=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell.
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetry elements in space group
!!   | occopt=option for occupancies
!!   | prtvol=integer controlling volume of printed output
!!   | tsmear=smearing energy or temperature (if metal)
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  grchempottn(3,natom)=grads of spatially-varying chemical potential energy (hartree)
!!  grcondft(3,natom)=grads of constrained DFT energy (hartree)
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D dispersion (hartree)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  ntypat=number of types of atoms in unit cell.
!!  nvresid(nfft,nspden)=potential or density residual
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of total energy
!!         (-1=no computation; 0=direct scheme; 1=double-counting scheme)
!!  optforces=option for the computation of forces
!!  optres=0 if residual array (nvresid) contains the potential residual
!!        =1 if residual array (nvresid) contains the density residual
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=array for holding XC potential
!!  vxctau(nfftf,dtset%nspden,4*usekden)]=derivative of XC energy density wrt
!!      kinetic energy density (metaGGA cases)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  deltae=change in total energy between the previous and present SCF cycle
!!  etotal=total energy (hartree)
!!  ===== if optforces==1
!!   diffor=maximum absolute change in component of forces between present and previous SCF cycle.
!!   favg(3)=mean of fcart before correction for translational symmetry
!!   fcart(3,natom)=cartesian forces from fred (hartree/bohr)
!!   fred(3,natom)=symmetrized form of grtn (grads of Etot) (hartree)
!!   gresid(3,natom)=forces due to the residual of the density/potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!   maxfor=maximum absolute value of force
!!   synlgr(3,natom)=symmetrized form of grads of Enl (hartree)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  elast=previous value of the energy,
!!        needed to compute deltae, then updated.
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_eigenvalues(IN)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_chempot(IN)=energy from spatially varying chemical potential (hartree)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_vdw_dftd(IN)=VdW DFT-D energy
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_hybcomp_E0(IN)=energy compensation energy for the hybrid functionals at frozen density
!!   | e_hybcomp_v0(IN)=potential compensation energy for the hybrid functionals at frozen density
!!   | e_hybcomp_v (IN)=potential compensation energy for the hybrid functionals at self-consistent density
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!   | e_nlpsp_vfock(IN)=nonlocal psp + potential Fock ACE part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the electric field:  enefield = -ucvol*E*P
!!   | e_magfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the magnetic field:  e_magfield = -ucvol*E*P
!!   | e_entropy(OUT)=entropy energy due to the occupation number smearing (if metal)
!!   |                this value is %entropy * dtset%tsmear (hartree).
!!  ===== if optforces==1
!!   forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!   grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!                 Input for norm-conserving psps, output for PAW
!!  ===== if psps%usepaw==1
!!   pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj
!!      lincom_cgcprj,pawcprj_alloc,pawcprj_axpby,pawcprj_free,pawcprj_get
!!      pawcprj_getdim,pawcprj_lincom,pawcprj_put,timab,xmpi_sum,zgesv
!!
!! SOURCE

subroutine etotfor(atindx1,deltae,diffor,dtefield,dtset,&
&  elast,electronpositron,energies,&
&  etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,grcondft,gresid,grewtn,grhf,grnl,grvdw,&
&  grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,my_natom,nattyp,&
&  nfft,ngfft,ngrvdw,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&
&  pawang,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,red_ptot,psps,rhog,rhor,rmet,rprimd,&
&  symrec,synlgr,ucvol,usepaw,vhartr,vpsp,vxc,vxctau,wvl,wvl_den,xccc3d,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,mgfft,n1xccc,n3xccc,nfft,ngrvdw,nkxc,ntypat,optene,optforces
 integer,intent(in) :: optres,usepaw
 real(dp),intent(in) :: gsqcut
 real(dp),intent(inout) :: elast,ucvol
 real(dp),intent(out) :: deltae,diffor,etotal,maxfor
 type(MPI_type),intent(in) :: mpi_enreg
 type(efield_type),intent(in) :: dtefield
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: gmet(3,3),grchempottn(3,dtset%natom),grcondft(3,dtset%natom)
 real(dp),intent(in) :: grewtn(3,dtset%natom),grvdw(3,ngrvdw),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),red_ptot(3)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rmet(3,3)
 real(dp),intent(in) :: vhartr(nfft),vpsp(nfft),vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
 real(dp),intent(in) :: xccc3d(n3xccc)
 real(dp),intent(inout) :: forold(3,dtset%natom),grnl(3*dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(inout),target :: nvresid(nfft,dtset%nspden)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fred(3,dtset%natom)
 real(dp),intent(inout) :: fcart(3,dtset%natom)
 real(dp),intent(inout) :: rprimd(3,3)
 real(dp),intent(out) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(inout) :: grxc(3,dtset%natom)
 real(dp),intent(out) :: synlgr(3,dtset%natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: comm_grid,dimnhat,ifft,ipositron,ispden,optgr,optgr2,option,optnc,optstr,optstr2,iir,jjr,kkr
 logical :: apply_residual
 real(dp) :: eenth,ucvol_
!arrays
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: tsec(2),A(3,3),A1(3,3),A_new(3,3),efield_new(3)
 real(dp) :: dummy(0),nhat_dum(0,0)
 real(dp),allocatable :: vlocal(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: resid(:,:)

! *********************************************************************

 call timab(80,1,tsec)

 ipositron=electronpositron_calctype(electronpositron)

 if (optene>-1) then

!  When the finite-temperature VG broadening scheme is used,
!  the total entropy contribution "tsmear*entropy" has a meaning,
!  and gather the two last terms of Eq.8 of VG paper
!  Warning : might have to be changed for fixed moment calculations
   if(dtset%occopt>=3 .and. dtset%occopt<=8) then
     if (abs(dtset%tphysel) < tol10) then
       energies%e_entropy = - dtset%tsmear * energies%entropy
     else
       energies%e_entropy = - dtset%tphysel * energies%entropy
     end if
   else
     energies%e_entropy = zero
   end if

!  Turn it into an electric enthalpy, refer to Eq.(33) of Suppl. of Nat. Phys. paper (5,304,2009) [[cite:Stengel2009]],
!    the missing volume is added here
   energies%e_elecfield=zero
   if ((dtset%berryopt==4.or.dtset%berryopt==14).and.ipositron/=1) then
     energies%e_elecfield=-dot_product(dtset%red_efieldbar,red_ptot)  !!ebar_i p_i
     eenth=zero
     do iir=1,3
       do jjr=1,3
         eenth=eenth+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)  !!g^{-1})_ij ebar_i ebar_j
       end do
     end do
     energies%e_elecfield=energies%e_elecfield-eenth*ucvol/(8._dp*pi)
   end if

!  Turn it into an internal energy, refer to Eq.(36) of Suppl. of Nat. Phys. paper (5,304,2009) [[cite:Stengel2009]],
!    but a little different: U=E_ks + (vol/8*pi) *  g^{-1})_ij ebar_i ebar_j
   if ((dtset%berryopt==6.or.dtset%berryopt==16).and.ipositron/=1) then
     energies%e_elecfield=zero
     eenth=zero
     do iir=1,3
       do jjr=1,3
         eenth=eenth+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)  !! g^{-1})_ij ebar_i ebar_j
       end do
     end do
     energies%e_elecfield=energies%e_elecfield+eenth*ucvol/(8._dp*pi)
   end if

!  Calculate internal energy and electric enthalpy for mixed BC case.
   if (dtset%berryopt==17.and.ipositron/=1) then
     energies%e_elecfield=zero
     A(:,:)=(four_pi/ucvol)*rmet(:,:)
     A1(:,:)=A(:,:) ; A_new(:,:)=A(:,:)
     efield_new(:)=dtset%red_efield(:)
     eenth=zero
     do kkr=1,3
       if (dtset%jfielddir(kkr)==1) then    ! fixed ebar direction
!        step 1 add -ebar*p
         eenth=eenth-dtset%red_efieldbar(kkr)*red_ptot(kkr)
!        step 2  chang to e_new (change e to ebar)
         efield_new(kkr)=dtset%red_efieldbar(kkr)
!        step 3  chang matrix A to A1
         do iir=1,3
           do jjr=1,3
             if (iir==kkr .and. jjr==kkr) A1(iir,jjr)=-1.0/A(kkr,kkr)
             if ((iir==kkr .and. jjr/=kkr) .or.  (iir/=kkr .and.  jjr==kkr)) &
&             A1(iir,jjr)=-1.0*A(iir,jjr)/A(kkr,kkr)
             if (iir/=kkr .and. jjr/=kkr) A1(iir,jjr)=A(iir,jjr)-A(iir,kkr)*A(kkr,jjr)/A(kkr,kkr)
           end do
         end do
         A(:,:)=A1(:,:) ; A_new(:,:)=A1(:,:)
       end if
     end do  ! end fo kkr
     do iir=1,3
       do jjr=1,3
         eenth= eenth+half*A_new(iir,jjr)*efield_new(iir)*efield_new(jjr)
       end do
     end do
     energies%e_elecfield=energies%e_elecfield+eenth
   end if   ! berryopt==17

!  Turn it into a magnetic enthalpy, by adding orbital electronic contribution
   energies%e_magfield = zero
!  if (dtset%berryopt == 5 .and. ipositron/=1) then
!  emag = dot_product(mag_cart,dtset%bfield)
!  energies%e_magfield = emag
!  end if

!  Compute total (free)- energy by direct scheme
   if (optene==0) then
     etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
      energies%e_localpsp + energies%e_corepsp +&
!&    two*energies%e_fock - energies%e_fock0 +&   ! The Fock energy is already included in the non-local one
!&     energies%e_nlpsp_vfock - energies%e_fock0 +&
&     energies%e_entropy + energies%e_elecfield + energies%e_magfield+&
&     energies%e_hybcomp_E0 - energies%e_hybcomp_v0 + energies%e_hybcomp_v + energies%e_constrained_dft
     etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd

!    See similar section in m_energies.F90
!    XG 20181025 This gives a variational energy in case of NCPP with all bands occupied - not yet for metals.
     if (usepaw==0) etotal = etotal + energies%e_nlpsp_vfock - energies%e_fock0
!    XG 20181025 I was expecting the following to give also a variational energy in case of PAW, but this is not true.
!    if (usepaw==1) etotal = etotal + energies%e_paw + energies%e_nlpsp_vfock - energies%e_fock0
!    XG 20181025 So, the following is giving a non-variational expression ...
     if (usepaw==1) etotal = etotal + energies%e_paw + energies%e_fock

   end if

!  Compute total (free) energy by double-counting scheme
   if (optene==1) then
     etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc &
&     - energies%e_xcdc + energies%e_corepsp - energies%e_corepspdc- energies%e_fock0 &
&     + energies%e_entropy + energies%e_elecfield + energies%e_magfield &
&     + energies%e_hybcomp_E0 - energies%e_hybcomp_v0 + energies%e_constrained_dft
     etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
     if (usepaw/=0) etotal = etotal + energies%e_pawdc
   end if

!  Additional stuff for electron-positron
   if (dtset%positron/=0) then
     if (ipositron==0) then
       energies%e_electronpositron  =zero
       energies%edc_electronpositron=zero
     else
       energies%e_electronpositron  =electronpositron%e_hartree+electronpositron%e_xc
       energies%edc_electronpositron=electronpositron%e_hartree+electronpositron%e_xcdc
       if (usepaw==1) then
         energies%e_electronpositron  =energies%e_electronpositron  +electronpositron%e_paw
         energies%edc_electronpositron=energies%edc_electronpositron+electronpositron%e_pawdc
       end if
     end if
     if (optene==0) electronpositron%e0=etotal
     if (optene==1) electronpositron%e0=etotal-energies%edc_electronpositron
     etotal=electronpositron%e0+energies%e0_electronpositron+energies%e_electronpositron
   end if

!  Compute energy residual
   deltae=etotal-elast
   elast=etotal
 end if !optene/=-1

 call timab(80,2,tsec)

!------Compute forces-----------------------------------------------------

 if (optforces==1) then

!  PAW: add gradients due to Dij derivatives to non-local term
   if (usepaw==1) then
     ABI_ALLOCATE(vlocal,(nfft,dtset%nspden))
     do ispden=1,min(dtset%nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,vhartr,vlocal,vpsp,vxc)
       do ifft=1,nfft
         vlocal(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
       end do
     end do

     if(dtset%nspden==4)then
       do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,vlocal,vxc)
         do ifft=1,nfft
           vlocal(ifft,ispden)=vxc(ifft,ispden)
         end do
       end do
     end if
     ucvol_=ucvol
#if defined HAVE_BIGDFT
     if (dtset%usewvl==1) ucvol_=product(wvl_den%denspot%dpbox%hgrids)*real(product(wvl_den%denspot%dpbox%ndims),dp)
#endif
     dimnhat=0;optgr=1;optgr2=0;optstr=0;optstr2=0
     comm_grid=mpi_enreg%comm_fft;if(dtset%usewvl==1) comm_grid=mpi_enreg%comm_wvl
     call pawgrnl(atindx1,dimnhat,dummy,1,dummy,grnl,gsqcut,mgfft,my_natom, &
&     dtset%natom, nattyp,nfft,ngfft,nhat_dum,dummy,dtset%nspden,dtset%nsym,ntypat,optgr,optgr2,optstr,optstr2,&
&     pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,k0,rprimd,symrec,dtset%typat,ucvol_,vlocal,vxc,xred, &
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,mpi_comm_grid=mpi_enreg%comm_fft,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,paral_kgb=mpi_enreg%paral_kgb)
     ABI_DEALLOCATE(vlocal)
   end if

   apply_residual=(optres==1 .and. dtset%usewvl==0.and.abs(dtset%densfor_pred)>=1 .and. &
&   abs(dtset%densfor_pred)<=6.and.abs(dtset%densfor_pred)/=5)

!  If residual is a density residual (and forces from residual asked),
!  has to convert it into a potential residual before calling forces routine
   if (apply_residual) then
     ABI_ALLOCATE(resid,(nfft,dtset%nspden))
     option=0; if (dtset%densfor_pred<0) option=1
     optnc=1;if (dtset%nspden==4.and.(abs(dtset%densfor_pred)==4.or.abs(dtset%densfor_pred)==6)) optnc=2
     call nres2vres(dtset,gsqcut,usepaw,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&
&     nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&     rhor,rprimd,usepaw,resid,xccc3d,xred,vxc)
   else
     resid => nvresid
   end if
   call forces(atindx1,diffor,dtefield,dtset,favg,fcart,fock,forold,fred,grchempottn,grcondft,gresid,grewtn,&
&   grhf,grnl,grvdw,grxc,gsqcut,indsym,maxfor,mgfft,mpi_enreg,&
&   n1xccc,n3xccc,nattyp,nfft,ngfft,ngrvdw,ntypat,pawrad,pawtab,&
&   ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,dtset%usefock,resid,vxc,vxctau,wvl,wvl_den,xred,&
&   electronpositron=electronpositron)
   if (apply_residual) then
     ABI_DEALLOCATE(resid)
   end if

!  Returned fred are full symmetrized gradients of Etotal
!  wrt reduced coordinates xred, d(Etotal)/d(xred)
!  Forces are contained in array fcart

 else   ! if optforces==0
   fcart=zero
   fred=zero
   favg=zero
   diffor=zero
   gresid=zero
   grhf=zero
   maxfor=zero
   synlgr=zero
 end if

 call timab(80,2,tsec)

end subroutine etotfor
!!***

!!****f* ABINIT/wf_mixing
!! NAME
!! wf_mixing
!!
!! FUNCTION
!! Mixing of wavefunctions in the outer loop of a double loop SCF approach.
!! Different algorithms are implemented, depending on the value of wfmixalg.
!!
!! INPUTS
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  istep=number of call the routine (usually the outer loop in the SCF double loop)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of cprj array
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(dtset%ntypat)=number of atoms of each type in cell.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  scf_history_wf <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj
!!      lincom_cgcprj,pawcprj_alloc,pawcprj_axpby,pawcprj_free,pawcprj_get
!!      pawcprj_getdim,pawcprj_lincom,pawcprj_put,timab,xmpi_sum,zgesv
!!
!! SOURCE

subroutine wf_mixing(atindx1,cg,cprj,dtset,istep,mcg,mcprj,mpi_enreg,&
& nattyp,npwarr,pawtab,scf_history_wf)

 use m_cgcprj,  only : dotprod_set_cgcprj, dotprodm_sumdiag_cgcprj, lincom_cgcprj, cgcprj_cholesky

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history_wf
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: hermitian
 integer :: ibdmix,ibdsp,ibg,ibg_hist,icg,icg_hist
 integer :: ierr,ikpt,indh,ind_biorthog,ind_biorthog_eff,ind_newwf,ind_residual,inplace
 integer :: iorder,iset2,isppol,istep_cycle,istep_new,istwf_k,kk,me_distrb,my_nspinor
 integer :: nband_k,nbdmix,npw_k,nset1,nset2,ntypat
 integer :: shift_set1,shift_set2,spaceComm_band,spare_mem,usepaw,wfmixalg
 real(dp) :: alpha,beta
 complex(dpc) :: sum_coeffs
!arrays
 integer,allocatable :: ipiv(:),dimcprj(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: al(:,:),mmn(:,:,:)
 real(dp),allocatable :: dotprod_res(:,:,:),dotprod_res_k(:,:,:),res_mn(:,:,:),smn(:,:,:)
 complex(dpc),allocatable :: coeffs(:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kh(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)
!write(std_out,*)' wf_mixing : enter, istep= ',istep
!call flush(std_out)
!write(std_out,*)' istep,scf_history_wf%alpha=',istep,scf_history_wf%alpha
!write(std_out,*)' cg(1,1)=',cg(1,1)
!write(std_out,*)' scf_history_wf%cg(1,1,1:5)=',scf_history_wf%cg(1,1,1:5)
!ABI_ALLOCATE(cg_ref,(2,mcg))
!cg_ref(:,:)=cg(:,:)
!ABI_DATATYPE_ALLOCATE(cprj_ref,(dtset%natom,mcprj))
!cprj_ref(:,:)=cprj(:,:)
!      write(std_out,*)' scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)=',&
!&       scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)
!       call flush(std_out)
!ENDDEBUG

 if (istep==0) return

 ntypat=dtset%ntypat
 usepaw=dtset%usepaw
 wfmixalg=scf_history_wf%wfmixalg
 nbdmix=dtset%nbandhf
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 me_distrb=mpi_enreg%me_kpt
 spaceComm_band=xmpi_comm_self

 spare_mem=0
 if(scf_history_wf%history_size==wfmixalg-1)spare_mem=1

!scf_history_wf%alpha contains dtset%wfmix
 alpha=scf_history_wf%alpha
 beta=one-scf_history_wf%alpha
 icg=0
 icg_hist=0
 ibg=0
 ibg_hist=0

!Useful array
 ABI_ALLOCATE(dimcprj,(dtset%natom))
 if (usepaw==1) then
   iorder=0 ! There is no change of ordering in the mixing of wavefunctions
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

 if(istep==1)then
   do indh=1,scf_history_wf%history_size
     call pawcprj_alloc(scf_history_wf%cprj(:,:,indh),0,dimcprj)
   end do
 end if

 ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,my_nspinor*nbdmix))
 ABI_DATATYPE_ALLOCATE(cprj_kh,(dtset%natom,my_nspinor*nbdmix))
 if(usepaw==1) then
   call pawcprj_alloc(cprj_k,0,dimcprj)
   call pawcprj_alloc(cprj_kh,0,dimcprj)
 end if
 ABI_ALLOCATE(smn,(2,nbdmix,nbdmix))
 ABI_ALLOCATE(mmn,(2,nbdmix,nbdmix))

 if(wfmixalg>2)then
   nset1=1
   nset2=min(istep-1,wfmixalg-1)
   ABI_ALLOCATE(dotprod_res_k,(2,1,nset2))
   ABI_ALLOCATE(dotprod_res,(2,1,nset2))
   ABI_ALLOCATE(res_mn,(2,wfmixalg-1,wfmixalg-1))
   dotprod_res=zero
   if(istep==1)then
     scf_history_wf%dotprod_sumdiag_cgcprj_ij=zero
   end if
 end if

!Explanation for the index for the wavefunction stored in scf_history_wf
!The reference is the cg+cprj output after the wf optimization at istep 1.
!It comes as input to the present routine as cgcprj input at step 2, and is usually found at indh=1.

!In the simple mixing case (wfmixalg==2), the reference is never stored, because it is used "on-the-fly" to biothogonalize the
!previous input (that was stored in indh=1), then generate the next input, which is stored again in indh=1

!When the storage is not spared:
!- the values of indh from 2 to wfmixalg store the (computed here) biorthogonalized input cgcprj, then the residual
!- the values of indh from wfmixalg+1 to 2*wfmixalg-1 store the biorthogonalized output cgcprj (coming as argument)

!First step
 if (istep==1 .or. (wfmixalg==2 .and. abs(scf_history_wf%alpha-one)<tol8) ) then

   indh=2   ! This input wavefunction is NOT the reference
   if(wfmixalg==2)indh=1 ! But this does not matter in the simple mixing case that has history_size=1

!  Simply store the wavefunctions and cprj. However, nband_k might be different from nbandhf...
!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       npw_k=npwarr(ikpt)

       scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
       if(usepaw==1) then
!        scf_history_wf%cprj(:,ibg_hist+1:ibg_hist+my_nspinor*nbdmix,1)=cprj(:,ibg+1:ibg+my_nspinor*nbdmix)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,iorder,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nbdmix,nband_k,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,iorder,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

!      Update the counters
       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

     end do
   end do

 else
!  From istep==2

!  First part of the computation : biorthogonalization, and computation of the residual (possibly, prediction of the next input in the case of simple mixing)
!  Index for the wavefunctions stored in scf_history_wf whose scalar products with the argument cgcprj will have to be computed.
   indh=1   ! This input wavefunction is the reference
   if(wfmixalg/=2 .and. istep==2)indh=2 ! except for istep=2 in the rmm-diis

   if(wfmixalg>2)then
!    istep inside the cycle defined by wfmixalg, and next index. Then, indices of the wavefunction sets.
     istep_cycle=mod((istep-2),wfmixalg-1)
     istep_new=mod((istep-1),wfmixalg-1)
     ind_biorthog=1+wfmixalg+istep_cycle
     ind_residual=2+istep_cycle
     ind_newwf=2+istep_new
     shift_set1=ind_residual-1
     shift_set2=1
   end if

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)

!      Biorthogonalization

       if(usepaw==1) then
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,iorder,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nbdmix,nband_k,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_get(atindx1,cprj_kh,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,iorder,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       hermitian=0
       if(wfmixalg==2 .or. istep==2)then
         call dotprod_set_cgcprj(atindx1,cg,scf_history_wf%cg(:,:,indh),cprj_k,cprj_kh,dimcprj,hermitian,&
&         0,0,icg,icg_hist,ikpt,isppol,istwf_k,nbdmix,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       else
         call dotprod_set_cgcprj(atindx1,scf_history_wf%cg(:,:,indh),cg,cprj_kh,cprj_k,dimcprj,hermitian,&
&         0,0,icg_hist,icg,ikpt,isppol,istwf_k,nbdmix,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)
       end if

!      Invert S matrix, that is NOT hermitian.
!      Calculate M=S^-1
       mmn=zero
       do kk=1,nbdmix
         mmn(1,kk,kk)=one
       end do

       ABI_ALLOCATE(ipiv,(nbdmix))
!      The smn is destroyed by the following inverse call
       call zgesv(nbdmix,nbdmix,smn,nbdmix,ipiv,mmn,nbdmix,ierr)
       ABI_DEALLOCATE(ipiv)
!DEBUG
       if(ierr/=0)then
         MSG_ERROR(' The call to cgesv general inversion routine failed')
       end if
!ENDDEBUG

!      The M matrix is used to compute the biorthogonalized set of wavefunctions, and to store it at the proper place
       if(wfmixalg==2 .or. istep==2)then
         inplace=1
         call lincom_cgcprj(mmn,scf_history_wf%cg(:,:,indh),cprj_kh,dimcprj,&
&         icg_hist,inplace,mcg,my_nspinor*nbdmix,dtset%natom,nbdmix,nbdmix,npw_k,my_nspinor,usepaw)
       else
         inplace=0
         call lincom_cgcprj(mmn,cg,cprj_k,dimcprj,&
&         icg,inplace,mcg,my_nspinor*nbdmix,dtset%natom,nbdmix,nbdmix,npw_k,my_nspinor,usepaw,&
&         cgout=scf_history_wf%cg(:,:,ind_biorthog),cprjout=scf_history_wf%cprj(:,:,ind_biorthog),icgout=icg_hist)
       end if

!      The biorthogonalised set of wavefunctions is now stored at the proper place

!      Finalize this first part of the computation, depending on the algorithm and the step.

       if(wfmixalg==2)then

!        Wavefunction extrapolation, simple mixing case
!        alpha contains dtset%wfmix, beta contains one-alpha
         cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)=&
&         alpha*cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)&
&         +beta*scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_axpby(beta,alpha,cprj_kh(:,ibdmix:ibdmix),cprj_k(:,ibdmix:ibdmix))
           end do ! end loop on ibdmix
         end if

!        Back to usual orthonormalization
         call cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf_k,mcg,my_nspinor*nband_k,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,usepaw)

!        Store the newly extrapolated wavefunctions, orthonormalized, in scf_history_wf
         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,iorder,isppol,&
&             nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&             mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
           end do ! end loop on ibdmix
         end if

       else  !  wfmixalg/=2
!        RMM-DIIS

         if (istep==2)then
!          Store the argument wf as the reference for all future steps, in scf_history_wf with index 1.
           scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,1)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
           if(usepaw==1) then
             do ibdmix=1,nbdmix
               call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,1),dtset%natom,1,ibg_hist,ikpt,iorder,isppol,&
&               nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&               mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
             end do ! end loop on ibdmix
           end if
         end if

         ind_biorthog_eff=ind_biorthog
         if(istep==2)ind_biorthog_eff=1 ! The argument wf has not been stored in ind_biorthog
!        Compute the residual of the wavefunctions for this istep,
!        that replaces the previously stored set of (biorthogonalized) input wavefunctions
         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_residual)=&
&         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_biorthog_eff)&
&         -scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_residual)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_axpby(one,-one,scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog_eff),&
&             scf_history_wf%cprj(:,ibdmix:ibdmix,ind_residual))
           end do ! end loop on ibdmix
         end if

!        Compute the new scalar products to fill the res_mn matrix
         call dotprodm_sumdiag_cgcprj(atindx1,scf_history_wf%cg,scf_history_wf%cprj,dimcprj,&
&         ibg_hist,icg_hist,ikpt,isppol,istwf_k,nbdmix,mcg,mcprj,dtset%mkmem,&
&         mpi_enreg,scf_history_wf%history_size,dtset%natom,nattyp,nbdmix,npw_k,nset1,nset2,my_nspinor,dtset%nsppol,ntypat,&
&         shift_set1,shift_set2,pawtab,dotprod_res_k,usepaw)

         dotprod_res=dotprod_res+dotprod_res_k

!        scf_history_wf for index ind_biorthog will contain the extrapolated wavefunctions (and no more the output of the SCF loop).
         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_biorthog)=&
&         scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_biorthog_eff)+&
&         (alpha-one)*scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_residual)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             if(ind_biorthog/=ind_biorthog_eff)then
               scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog)=scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog_eff)
             end if
             call pawcprj_axpby((alpha-one),one,scf_history_wf%cprj(:,ibdmix:ibdmix,ind_residual),&
&             scf_history_wf%cprj(:,ibdmix:ibdmix,ind_biorthog))
           end do ! end loop on ibdmix
         end if

       end if

       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

     end do ! End big k point loop
   end do ! End loop over spins

 end if ! istep>=2

 if(wfmixalg>2 .and. istep>1)then

!DEBUG
!  write(std_out,*)' '
!  write(std_out,*)' Entering the residual minimisation part '
!  write(std_out,*)' '
!  call flush(std_out)
!ENDDEBUG

   call timab(48,1,tsec)
   call xmpi_sum(dotprod_res,mpi_enreg%comm_kpt,ierr)
   call timab(48,2,tsec)

   scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,1+shift_set1,1+shift_set2:nset2+shift_set2)=dotprod_res(:,1,1:nset2)
   scf_history_wf%dotprod_sumdiag_cgcprj_ij(1,1+shift_set2:nset2+shift_set2,1+shift_set1)=dotprod_res(1,1,1:nset2)
   scf_history_wf%dotprod_sumdiag_cgcprj_ij(2,1+shift_set2:nset2+shift_set2,1+shift_set1)=-dotprod_res(2,1,1:nset2)

 end if ! wfmixalg>2 and istep>1

 if(wfmixalg>2 .and. istep>2)then

!  Extract the relevant matrix R_mn
   res_mn(:,1:nset2,1:nset2)=&
&   scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,1+shift_set2:nset2+shift_set2,1+shift_set2:nset2+shift_set2)

!DEBUG
!      write(std_out,*)' The matrix res_mn(:,1:nset2,1:nset2) is :'
!      write(std_out,*)res_mn(:,1:nset2,1:nset2)
!      call flush(std_out)
!ENDDEBUG

!  Solve R_mn \alpha_n = 1_m
   ABI_ALLOCATE(ipiv,(nset2))
   ABI_ALLOCATE(coeffs,(nset2))
   coeffs(:)=cone
!  The res_mn is destroyed by the following inverse call
   call zgesv(nset2,1,res_mn,wfmixalg-1,ipiv,coeffs,nset2,ierr)
   ABI_DEALLOCATE(ipiv)
!  The coefficients must sum to one
   sum_coeffs=sum(coeffs)
   coeffs=coeffs/sum_coeffs

!DEBUG
!      write(std_out,*)' The coefficients that minimize the residual have been found'
!      write(std_out,*)' coeffs =',coeffs
!      call flush(std_out)
!ENDDEBUG
 end if ! wfmixalg>2 and istep>2

 if(wfmixalg>2 .and. istep>1)then

!  Find the new "input" wavefunction, bi-orthogonalized, and store it replacing the adequate "old" input wavefunction.

   icg=0
   icg_hist=0
   ibg=0
   ibg_hist=0
   ABI_ALLOCATE(al,(2,nset2))
   if(istep>2)then
     do iset2=1,nset2
       al(1,iset2)=real(coeffs(iset2)) ; al(2,iset2)=aimag(coeffs(iset2))
     end do
   else
     al(1,1)=one ; al(2,1)=zero
   end if

!DEBUG
!      write(std_out,*)' Overload the coefficients, in order to simulate a simple mixing with wfmix '
!      write(std_out,*)' Set al(1,ind_biorthog-3)=one, for ind_biorthog=',ind_biorthog
!      write(std_out,*)' This will feed scf_history for set ind_biorthog-3+wfmixalg=',ind_biorthog-3+wfmixalg
!      al(:,:)=zero
!      al(1,ind_biorthog-3)=one
!      call flush(std_out)
!ENDDEBUG

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)

       if(istep>2)then
!        Make the appropriate linear combination (from the extrapolated wfs)
         cg(:,icg+1:icg+my_nspinor*npw_k*nband_k)=zero
         do iset2=1,nset2
           cg(1,icg+1:icg+my_nspinor*npw_k*nband_k)=cg(1,icg+1:icg+my_nspinor*npw_k*nband_k)&
&           +al(1,iset2)*scf_history_wf%cg(1,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)&
&           -al(2,iset2)*scf_history_wf%cg(2,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)
           cg(2,icg+1:icg+my_nspinor*npw_k*nband_k)=cg(2,icg+1:icg+my_nspinor*npw_k*nband_k)&
&           +al(1,iset2)*scf_history_wf%cg(2,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)&
&           +al(2,iset2)*scf_history_wf%cg(1,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,iset2+wfmixalg)
         end do
       else ! One needs a simple copy from the extrapolated wavefunctions
         cg(:,icg+1:icg+my_nspinor*npw_k*nband_k)=scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,1+wfmixalg)
       end if
!      Note the storage in cprj_k. By the way, a simple copy might also be used in case istep=2.
       if(usepaw==1) then
         do ibdsp=1,my_nspinor*nbdmix
           call pawcprj_lincom(al,scf_history_wf%cprj(:,ibdsp,1+wfmixalg:nset2+wfmixalg),cprj_k(:,ibdsp:ibdsp),nset2)
         end do
       end if

!      Store the newly extrapolated wavefunctions for this k point, still bi-orthonormalized, in scf_history_wf
       scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind_newwf)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)
       if(usepaw==1) then
         call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,ind_newwf),dtset%natom,1,ibg_hist,ikpt,iorder,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

!      Back to usual orthonormalization for the cg and cprj_k
       call cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf_k,mcg,my_nspinor*nband_k,dtset%mkmem,&
&       mpi_enreg,dtset%natom,nattyp,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,usepaw)

!      Need to transfer cprj_k to cprj
       if(usepaw==1) then
         call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,iorder,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

     end do ! End big k point loop
   end do ! End loop over spins

   if(istep>2)then
     ABI_DEALLOCATE(coeffs)
   end if
   ABI_DEALLOCATE(al)

 end if ! wfmixalg>2 and istep>1

!DEBUG
! write(std_out,*)' wf_mixing : exit '
!      write(std_out,*)' scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)=',&
!&       scf_history_wf%dotprod_sumdiag_cgcprj_ij(:,2,2)
! write(std_out,*)' cg(1:2,1:2)=',cg(1:2,1:2)
! write(std_out,*)' scf_history_wf%cg(1:2,1:2,1)=',scf_history_wf%cg(1:2,1:2,1)
! ABI_DEALLOCATE(cg_ref)
! ABI_DATATYPE_DEALLOCATE(cprj_ref)
!ENDDEBUG

 if(usepaw==1) then
   call pawcprj_free(cprj_k)
   call pawcprj_free(cprj_kh)
 end if
 ABI_DATATYPE_DEALLOCATE(cprj_k)
 ABI_DATATYPE_DEALLOCATE(cprj_kh)
 ABI_DEALLOCATE(dimcprj)
 ABI_DEALLOCATE(mmn)
 ABI_DEALLOCATE(smn)
 if(wfmixalg>2)then
   ABI_DEALLOCATE(dotprod_res_k)
   ABI_DEALLOCATE(dotprod_res)
   ABI_DEALLOCATE(res_mn)
 end if

end subroutine wf_mixing
!!***

end module m_scfcv_core
!!***

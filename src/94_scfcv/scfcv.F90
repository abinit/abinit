!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv
!! NAME
!! scfcv
!!
!! FUNCTION
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
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
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
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
!!  dtefield <type(efield_type)> = variables related to Berry phase
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
!!      ab7_mixing_deallocate,ab7_mixing_new,ab7_mixing_use_disk_cache
!!      afterscfloop,build_vxc,check_kxc,chkpawovlp,cprj_clean,cprj_paw_alloc
!!      ctocprj,destroy_distribfft,destroy_mpi_enreg,energies_init,energy
!!      etotfor,extraprho,fftdatar_write_from_hdr,first_rec,fock_destroy
!!      fock_init,fock_update_exc,fock_updatecwaveocc,fourdp,fresid,getcut
!!      getmpw,getng,getph,gshgg_mkncwrite,hdr_update,init_distribfft
!!      init_distribfft_seq,init_metricrec,initmpi_seq,initylmg,int2char4,kpgio
!!      metric,newrho,newvtr,nhatgrid,odamix,out_geometry_xml,out_resultsgs_xml
!!      outscfcv,paw2wvl_ij,paw_an_free,paw_an_init,paw_an_nullify
!!      paw_an_reset_flags,paw_ij_free,paw_ij_init,paw_ij_nullify
!!      paw_ij_reset_flags,pawcprj_alloc,pawcprj_free,pawcprj_getdim
!!      pawcprj_reorder,pawdenpot,pawdij,pawfgrtab_free,pawfgrtab_init
!!      pawmknhat,pawtab_get_lsize,pawuj_red,prc_mem_free,prtene,psolver_rhohxc
!!      rhohxc,rhotov,scprqt,setnoccmmp,setrhoijpbe0,setsym,setup_positron
!!      setvtr,sphereboundary,status,symdij,symmetrize_xred,timab
!!      update_e_field_vars,vtorho,vtorhorec,vtorhotf,wrtout,wvl_cprjreorder
!!      wvl_nhatgrid,xmpi_isum,xmpi_sum,xmpi_wait
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine scfcv(atindx,atindx1,cg,cpus,dmatpawu,dtefield,dtfil,dtpawuj,&
&  dtset,ecore,eigen,electronpositron,fatvshift,hdr,indsym,&
&  initialized,irrzon,kg,mcg,mpi_enreg,my_natom,nattyp,ndtpawuj,nfftf,npwarr,occ,&
&  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,&
&  pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,taug,taur,wffnew,wvl,xred,xred_old,ylm,ylmgr,conv_retcode)


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_xmpi
 use m_profiling_abi
 use m_wffile
 use m_rec
 use m_ab7_mixing
 use m_errors
 use m_efield
 use mod_prc_memory
 use m_nctk
 use m_hdr

 use m_fstrings,         only : int2char4, sjoin
 use m_time,             only : abi_wtime, sec2str
 use m_exit,             only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_abi_etsf,         only : abi_etsf_init
 use m_mpinfo,           only : destroy_mpi_enreg, iwrite_fftdatar
 use m_ioarr,            only : ioarr, fftdatar_write_from_hdr
 use m_results_gs ,      only : results_gs_type
 use m_scf_history,      only : scf_history_type
 use m_energies,         only : energies_type, energies_init, energies_copy
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type,pawtab_get_lsize
 use m_paw_an,           only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify,&
&                               paw_an_reset_flags
 use m_pawfgrtab,        only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_reorder, pawcprj_getdim
 use m_pawdij,           only : pawdij, symdij,pawdijhat
 use m_pawfgr,           only : pawfgr_type
 use m_paw_ij,           only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify,&
&                               paw_ij_reset_flags
 use m_paw_dmft,         only : paw_dmft_type
 use m_fock,             only : fock_type,fock_init,fock_destroy,fock_update_exc,fock_updatecwaveocc
 use gwls_hamiltonian,   only : build_vxc
#if defined HAVE_BIGDFT
 use BigDFT_API,         only : cprj_clean,cprj_paw_alloc
#endif
 use m_io_kss,             only : gshgg_mkncwrite

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcv'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
 use interfaces_43_wvl_wrappers
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_57_iovars
 use interfaces_62_poisson
 use interfaces_65_paw
 use interfaces_66_nonlocal
 use interfaces_67_common
 use interfaces_68_recursion
 use interfaces_68_rsprc
 use interfaces_79_seqpar_mpi
 use interfaces_94_scfcv, except_this_one => scfcv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,my_natom,ndtpawuj,pwind_alloc
 integer,intent(inout) :: initialized,nfftf
 integer,intent(out) :: conv_retcode
 real(dp),intent(in) :: cpus,ecore,fatvshift
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
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
 real(dp), intent(in) :: rprimd(3,3)
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
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(inout) :: paw_dmft

!Local variables -------------------------
!scalars
 integer,parameter :: level=110,response=0,cplex1=1
 integer :: afford,bantot,choice
 integer :: computed_forces,cplex,cplex_hf,ctocprj_choice,dbl_nnsclo,dielop,dielstrt,dimdmat
 integer :: forces_needed,errid,has_dijhat,has_dijnd,has_vhartree,has_dijfock,usefock
!integer :: dtset_iprcel
 integer :: iatom,ider,idir,ierr,iexit,ii,ikpt,impose_dmat,denpot
 integer :: initialized0,iorder_cprj,ipert,ipositron,isave_den,isave_kden,iscf10,ispden
 integer :: ispmix,istep,istep_mix,itypat,izero,lmax_diel,lpawumax,mband_cprj,mcprj,me,me_wvl
 integer :: mgfftdiel,mgfftf,moved_atm_inside,moved_rhor,my_nspinor,n1xccc
 integer :: n3xccc,ncpgr,nfftdiel,nfftmix,nfftmix_per_nfft,nfftotf,ngrvdw,nhatgrdim,nk3xc,nkxc
 integer :: npawmix,npwdiel,nstep,nzlmopt,optcut,optcut_hf,optene,optgr0,optgr0_hf
 integer :: optgr1,optgr2,optgr1_hf,optgr2_hf,option,optrad,optrad_hf,optres,optxc,prtfor,prtxml,quit
 integer :: quit_sum,req_cplex_dij,rdwrpaw,shft,spaceComm,spaceComm_fft,spaceComm_wvl,spaceComm_grid
 integer :: stress_needed,sz1,sz2,unit_out
 integer :: usecprj,usexcnhat
 integer :: my_quit,quitsum_request,timelimit_exit
 integer ABI_ASYNC :: quitsum_async
 real(dp) :: boxcut,compch_fft,compch_sph,deltae,diecut,diffor,ecut
 real(dp) :: ecutf,ecutsus,edum,elast,etotal,evxc,fermie,gsqcut,hybrid_mixing,hybrid_mixing_sr
 real(dp) :: maxfor,res2,residm,ucvol,ucvol_local,val_max
 real(dp) :: val_min,vxcavg,vxcavg_dum
 real(dp) :: zion,wtime_step,now,prev
 character(len=10) :: tag
 character(len=1500) :: message
 !character(len=500) :: dilatmx_errmsg
 character(len=fnlen) :: fildata
 type(MPI_type) :: mpi_enreg_diel
 type(energies_type) :: energies
 type(ab7_mixing_object) :: mix
 logical,parameter :: VERBOSE=.FALSE.
 logical :: dummy_nhatgr
 logical :: finite_efield_flag=.false.
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
 integer,allocatable :: indsym_dum(:,:,:),symrec_dum(:,:,:)
 logical,pointer :: lmselect_ep(:,:)
 real(dp) :: dielar(7),dphase(3),dummy2(6),favg(3),gmet(3,3),gprimd(3,3)
 real(dp) :: kpt_diel(3),pel(3),pel_cg(3),pelev(3),pion(3),ptot(3),red_ptot(3) !!REC
 real(dp) :: rhodum(1),rmet(3,3),strsxc(6),strten(6),tollist(12)
 real(dp) :: tsec(2),vnew_mean(dtset%nspden),vres_mean(dtset%nspden)
 real(dp) :: efield_old_cart(3), ptot_cart(3)
 real(dp) :: red_efield2(3),red_efield2_old(3)
 real(dp) :: vpotzero(2)
! red_efield1(3),red_efield2(3) is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009)
! red_efield1(3) for fixed ebar calculation, red_efield2(3) for fixed reduced d calculation, in mixed BC
! red_efieldbar_lc(3) is local reduced electric field, defined by Eq.(28) of Nat. Phys. suppl. (2009)
! pbar(3) and dbar(3) are reduced polarization and displacement field, defined by Eq.(27) and (29) Nat. Phys. suppl. (2009)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:)
 real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
 real(dp),allocatable :: grchempottn(:,:),grewtn(:,:),grhf(:,:),grnl(:),grvdw(:,:),grxc(:,:)
 real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),nvresid(:,:)
 real(dp),allocatable :: ph1d(:,:),ph1ddiel(:,:),ph1df(:,:)
 real(dp),allocatable :: phnonsdiel(:,:,:),shiftvector(:)
 real(dp),allocatable :: susmat(:,:,:,:,:),synlgr(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),vxctau(:,:,:),workr(:,:),xccc3d(:),ylmdiel(:,:)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:)
 type(pawcprj_type),allocatable :: cprj(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(pawrhoij_type),pointer :: pawrhoij_ep(:)
 type(fock_type),pointer :: fock

! *********************************************************************

 _IBM6("Hello, I'm running on IBM6")

 DBG_ENTER("COLL")

 call timab(238,1,tsec)
 call timab(54,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter')

 ! enable time limit handler if not done in callers.
 if (enable_timelimit_in(ABI_FUNC) == ABI_FUNC) then
   write(std_out,*)"Enabling timelimit check in function: ",trim(ABI_FUNC)," with timelimit: ",trim(sec2str(get_timelimit()))
 end if

!######################################################################
!Initializations - Memory allocations
!----------------------------------------------------------------------
 lmax_diel = 0

!MPI communicators
 if ((xmpi_paral==1).and.(mpi_enreg%paral_hf==1)) then
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
 !if (dtset%nstep==0 .or. dtset%iscf < 0) then
 if ((dtset%nstep==0 .or. dtset%iscf < 0) .and. dtset%plowan_compute==0) then
   energies%e_fermie = results_gs%energies%e_fermie
   results_gs%fermie = results_gs%energies%e_fermie
   write(std_out,*)"in scfcv: results_gs%fermie: ",results_gs%fermie
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
 usexcnhat=0;mcprj=0
 initialized0=initialized
 if (dtset%tfkinfunc==12) tfw_activated=.true.
 ipert=0;idir=0;cplex=1
 istep_mix=1
 ipositron=electronpositron_calctype(electronpositron)
 unit_out=0;if (dtset%prtvol >= 10) unit_out=ab_out
 nfftotf=product(ngfftf(1:3))

 usecprj=0 ; iorder_cprj=0
 if (psps%usepaw==1) then
   if (associated(electronpositron)) then
     if (dtset%positron/=0.and.electronpositron%dimcprj>0) usecprj=1
   end if
   if (dtset%prtnabla>0) usecprj=1
   if (dtset%extrapwf>0) usecprj=1
   if (dtset%usewvl==1)  usecprj=1
   if (dtset%pawfatbnd>0)usecprj=1
   if (dtset%prtdos==3)  usecprj=1
   if (nstep==0) usecprj=0
   if (usefock==1)  usecprj=1
 end if

!Stresses and forces flags
 forces_needed=0;prtfor=0
 if ((dtset%optforces==1.or.dtset%ionmov==4.or.abs(tollist(3))>tiny(0._dp))) then
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
!for output of scfcv
 etotal = results_gs%etotal

!Entering a scfcv loop, printing data to XML file if required.
 prtxml=0;if (me==0.and.dtset%prtxml==1) prtxml=1
 if (prtxml == 1) then
!  scfcv() will handle a scf loop, so we output the scfcv markup.
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
!results_gs should not be used as input of scfcv
!HERE IS PRINTED THE FIRST LINE OF SCFCV

 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
& dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
& maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
& occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
& psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode)

!Various allocations (potentials, gradients, ...)
 ABI_ALLOCATE(forold,(3,dtset%natom))
 ABI_ALLOCATE(grchempottn,(3,dtset%natom))
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
 ngrvdw=0;if (dtset%vdw_xc>=5.and.dtset%vdw_xc<=7) ngrvdw=dtset%natom
 ABI_ALLOCATE(grvdw,(3,ngrvdw))
 grchempottn(:,:)=zero
 forold(:,:)=zero ; gresid(:,:)=zero ; pel(:)=zero
 vtrial(:,:)=zero; vxc(:,:)=zero
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))

!Allocations/initializations for PAW only
 lpawumax=-1

 !Allocate fake cprj -> valgrind is happier and so am I
 ABI_DATATYPE_ALLOCATE(cprj,(1,1))
 ABI_ALLOCATE(dimcprj,(1))
 dimcprj(1) = 1
 call pawcprj_alloc(cprj,0,dimcprj)
 ABI_DEALLOCATE(dimcprj)
 if(psps%usepaw==1) then
!  Variables/arrays related to the fine FFT grid
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
     message = 'You cannot simultaneously use ionmov=4 and such a PAW psp file !'
     MSG_ERROR(message)
   end if

!  Variables/arrays related to the PAW spheres
   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   call paw_an_nullify(paw_an)
   call paw_ij_nullify(paw_ij)
   has_dijhat=0;if (dtset%iscf==22) has_dijhat=1
   has_vhartree=0; if (dtset%prtvha > 0) has_vhartree=1
   has_dijfock=0; if (usefock==1) has_dijfock=1
   has_dijnd=0; req_cplex_dij=1
   if(any(abs(dtset%nucdipmom)>tol8)) then
     has_dijnd=1; req_cplex_dij=2
   end if
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,0,dtset%nspden,cplex,dtset%pawxcdev,&
&   dtset%typat,pawang,pawtab,has_vxc=1,has_vxc_ex=1,has_vhartree=has_vhartree,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,&
&   dtset%pawspnorb,dtset%natom,dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijfock=has_dijfock,has_dijhartree=1,has_dijnd=has_dijnd,has_dijso=1,has_dijhat=has_dijhat,&
&   has_pawu_occ=1,has_exexch_pot=1,req_cplex_dij=req_cplex_dij,&
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
     if (pawtab(itypat)%usepawu>0) lpawumax=max(pawtab(itypat)%lpawu,lpawumax)
   end do
   if (dtset%usedmatpu/=0.and.lpawumax>0) then
     if (2*lpawumax+1/=size(dmatpawu,1).or.2*lpawumax+1/=size(dmatpawu,2)) then
       message = 'Incorrect size for dmatpawu!'
       MSG_BUG(message)
     end if
   end if

!  Allocation of projected WF (optional)
   if (usecprj==1) then
     iorder_cprj=0
     mband_cprj=dtset%mband;if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
     mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
     !Was allocated above for valgrind sake so should always be true (safety)
     if (allocated(cprj)) then
       call pawcprj_free(cprj)
       ABI_DATATYPE_DEALLOCATE(cprj)
     end if
     ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
     ncpgr=0
     if (usefock==1) then
       ctocprj_choice = 1
       if (dtset%optforces /= 0 .and. dtset%optstress == 0) then
         ncpgr = 3 ; ctocprj_choice = 2
 !        else if (dtset%optstress /= 0) then
 !       ncpgr = 9 ; ctocprj_choice = 23
       end if
     end if
     call pawcprj_alloc(cprj,ncpgr,dimcprj_srt)
#if defined HAVE_BIGDFT
     if (dtset%usewvl==1) then
       ABI_DATATYPE_ALLOCATE(wvl%descr%paw%cprj,(dtset%natom,mcprj))
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
   if (nstep==0) nvresid=zero
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
       sz1=pawrhoij(iatom)%cplex*pawtab(itypat)%lmn2_size
       sz2=pawrhoij(iatom)%nspden
       ABI_ALLOCATE(pawrhoij(iatom)%rhoijres,(sz1,sz2))
       do ispden=1,pawrhoij(iatom)%nspden
         pawrhoij(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(pawrhoij(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       pawrhoij(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij(iatom)%nspden*pawtab(itypat)%lmnmix_sz*pawrhoij(iatom)%cplex
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
     call ab7_mixing_new(mix, iscf10, denpot, ispmix, nfftmix, dtset%nspden, npawmix, errid, message, dtset%npulayit)
     if (errid /= AB7_NO_ERROR) then
       MSG_ERROR(message)
     end if
     if (dtset%mffmem == 0) then
       call ab7_mixing_use_disk_cache(mix, dtfil%fnametmp_fft)
     end if
!   else if (dtset%iscf==0.and.dtset%usewvl==1) then
!     ispmix=AB7_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
   end if
 else
   ABI_ALLOCATE(nvresid,(0,0))
   ABI_ALLOCATE(dtn_pc,(0,0))
   ABI_ALLOCATE(grhf,(0,0))
 end if ! iscf>0

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
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=min(2*dtset%nspden-1,3)
!Eventually need kxc-LDA for residual forces (when density mixing is selected)
 if (dtset%iscf>=10.and.dtset%usewvl==0.and.forces_needed>0 .and. &
& abs(dtset%densfor_pred)>=1.and.abs(dtset%densfor_pred)<=6.and.abs(dtset%densfor_pred)/=5) then
   if (dtset%xclevel==1.or.dtset%densfor_pred>=0) nkxc=min(2*dtset%nspden-1,3)
   if (dtset%xclevel==2.and.dtset%nspden==2.and.dtset%densfor_pred<0) nkxc=23
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

 call status(0,dtfil%filstat,iexit,level,'berryphase    ')

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

! start SCF loop
 do istep=1,max(1,nstep)
   call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

   ! Handle time limit condition.
   if (istep == 1) prev = abi_wtime()
   if (istep  > 1) then
     now = abi_wtime()
     wtime_step = now - prev
     prev = now
     call wrtout(std_out, sjoin(" scfcv: previous iteration took ",sec2str(wtime_step)))

     if (have_timelimit_in(ABI_FUNC)) then
       if (istep > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then
           write(message,"(3a)")"Approaching time limit ",trim(sec2str(get_timelimit())),". Will exit istep loop in scfcv."
           MSG_COMMENT(message)
           call wrtout(ab_out, message, "COLL")
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
     call symmetrize_xred(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

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
         call extraprho(atindx,atindx1,cg,dtset,gmet,gprimd,gsqcut,&
&         scf_history%icall,kg,mcg,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&         my_natom,nattyp,nfftf,ngfftf,npwarr,psps%ntypat,pawrhoij,pawtab,&
&         ph1df,psps,psps%qgrid_vl,rhor,rprimd,scf_history,ucvol,&
&         psps%usepaw,xred,xred_old,ylm,psps%ziontypat,psps%znuclpsp)
       end if
       call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     end if

   end if ! moved_atm_inside==1 .or. istep==1

   !Initialize/Update data in the case of an Exact-exchange (Hartree-Fock) or hybrid XC calculation
   hybrid_mixing=zero;hybrid_mixing_sr=zero
   if (usefock==1) then
     if (istep==1) then
       ! Initialize data_type fock for the calculation
       cplex_hf=cplex
       if (psps%usepaw==1) cplex_hf=dtset%pawcpxocc
       call fock_init(atindx,cplex_hf,dtset,fock,gsqcut,indsym,kg,mpi_enreg,nattyp,npwarr,pawang,pawfgr,pawtab,rprimd)
       if (fock%usepaw==1) then
         optcut_hf = 0 ! use rpaw to construct local_pawfgrtab
         optgr0_hf = 0; optgr1_hf = 0; optgr2_hf = 0 ! dont need gY terms locally
         optrad_hf = 1 ! do store r-R
         call nhatgrid(atindx1,gmet,dtset%natom,dtset%natom,nattyp,ngfftf,psps%ntypat,&
&         optcut_hf,optgr0_hf,optgr1_hf,optgr2_hf,optrad_hf,fock%pawfgrtab,pawtab,&
&         rprimd,dtset%typat,ucvol,xred,typord=1)
         iatom=-1;idir=0
         call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,iatom,idir,&
&         iorder_cprj,dtset%istwfk,kg,dtset%kptns,dtset%mband,mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft, dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,&
&         dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,dtfil%unpaw,&
&         xred,ylm,ylmgr)
       end if
     end if

     if (fock%optfor) then
       fock%forces=zero
     end if
     hybrid_mixing=fock%hybrid_mixing ; hybrid_mixing_sr=fock%hybrid_mixing_sr
     ! Update data relative to the occupied states in fock
     call fock_updatecwaveocc(cg,cprj,dtset,fock,energies%e_exactX,indsym,istep,mcg,mcprj,mpi_enreg,nattyp,npwarr,occ,ucvol)
   end if ! usefock

!  Initialize/update data in the electron-positron case
   if (dtset%positron<0.or.(dtset%positron>0.and.istep==1)) then
     call setup_positron(atindx,atindx1,cg,cprj,dtefield,dtfil,dtset,ecore,eigen,&
&     etotal,electronpositron,energies,fock,forces_needed,fred,gmet,gprimd,&
&     grchempottn,grewtn,grvdw,gsqcut,hdr,initialized0,indsym,istep,istep_mix,kg,&
&     kxc,maxfor,mcg,mcprj,mgfftf,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfftf,ngrvdw,nhat,&
&     nkxc,npwarr,nvresid,occ,optres,paw_ij,pawang,pawfgr,pawfgrtab,&
&     pawrad,pawrhoij,pawtab,ph1df,ph1d,psps,rhog,rhor,rprimd,&
&     stress_needed,strsxc,symrec,ucvol,usecprj,vhartr,vpsp,vxc,&
&     xccc3d,xred,ylm,ylmgr)
     ipositron=electronpositron_calctype(electronpositron)
   end if

   if ((moved_atm_inside==1 .or. istep==1).or. (dtset%positron<0.and.istep_mix==1)) then
!    PAW only: we sometimes have to compute compensation density
!    and eventually add it to density from WFs
     nhatgrdim=0
     dummy_nhatgr = .False.
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
           call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
         end if
       end if
       call timab(558,2,tsec)
     end if

!    The following steps have been gathered in the setvtr routine:
!    - get Ewald energy and Ewald forces
!    - compute local ionic pseudopotential vpsp
!    - eventually compute 3D core electron density xccc3d
!    - eventually compute vxc and vhartr
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
&     taug=taug,taur=taur,vxctau=vxctau,add_tfw=tfw_activated)

     ! set the zero of the potentials here
     if(dtset%usepotzero==2) then
       vpsp(:) = vpsp(:) + ecore / ( zion * ucvol )
     end if

     if(dtset%optdriver==RUNL_GWLS) then
       call build_vxc(vxc,nfftf,dtset%nspden)
     end if

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
       if(initialized==0) then
         call first_rec(dtset,psps,rec_set)
       end if
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
     if (dtset%useexexch>0) then
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
&     hyb_mixing=hybrid_mixing,hyb_mixing_sr=hybrid_mixing_sr,&
&     electronpositron=electronpositron,vpotzero=vpotzero)

!    Correct the average potential with the calculated constant vpotzero
!    Correct the total energies accordingly
!    vpotzero(1) = -beta/ucvol
!    vpotzero(2) = -1/ucvol sum_ij rho_ij gamma_ij
     write(message,'(a,f14.6,2x,f14.6)') &
&     ' average electrostatic smooth potential [Ha] , [eV]',SUM(vpotzero(:)),SUM(vpotzero(:))*Ha_eV
     call wrtout(std_out,message,'COLL')
     vtrial(:,:)=vtrial(:,:)+SUM(vpotzero(:))
     if(option/=1)then
!      Fix the direct total energy (non-zero only for charged systems)
       energies%e_paw=energies%e_paw-SUM(vpotzero(:))*dtset%charge
!      Fix the double counting total energy accordingly (for both charged AND
!      neutral systems)
       energies%e_pawdc=energies%e_pawdc-SUM(vpotzero(:))*zion+vpotzero(2)*dtset%charge
     end if

!    PAW+U: impose density matrix if required
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
&     hyb_mixing=hybrid_mixing,hyb_mixing_sr=hybrid_mixing_sr,&
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
   if (dtset%usepawu>0.and.dtset%macro_uj>0.and.istep>1.and.ipositron/=1) then
     call pawuj_red(dtset,dtpawuj,fatvshift,my_natom,dtset%natom,dtset%ntypat,&
     paw_ij,pawrad,pawtab,ndtpawuj,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if

   call timab(241,2,tsec)


!  No need to continue and call vtorho, when nstep==0
   if(nstep==0)exit

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   call timab(56,1,tsec)

   if(dtset%iscf>=0)then
     write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,message,'COLL')
   end if

   if (dtset%useria == -4242) then
     call gshgg_mkncwrite(istep, dtset, dtfil, psps, hdr, pawtab, pawfgr, paw_ij, mpi_enreg, &
     rprimd, xred, eigen, npwarr, kg, ylm, ngfft, dtset%nfft, ngfftf, nfftf, vtrial) !electronpositron) ! Optional arguments
   end if

!  The next flag says whether the xred have to be changed in the current iteration
   moved_atm_inside=0
   ! /< Hack to remove iapp from scfcv
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
   ! /< Hack to remove iapp from scfcv

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
       write(message,'(5a)')&
&       'Although the computation of forces during electronic iterations',ch10,&
&       'was not required by user, it is done (required by the',ch10,&
&       'choice of ionmov input parameter).'
       MSG_WARNING(message)
     end if
     if (abs(tollist(3))+abs(tollist(7))>tiny(0._dp)) then
       write(message,'(5a)')&
&       'Although the computation of forces during electronic iterations',ch10,&
&       'was not required by user, it is done (required by the',ch10,&
&       '"toldff" or "tolrff" tolerance criteria).'
       MSG_WARNING(message)
     end if
   end if
   if ((istep==1).and.(dtset%optforces==1).and. dtset%usewvl == 1) then
     write(message,'(5a)')&
&     'Although the computation of forces during electronic iterations',ch10,&
&     'was required by user, it has been disable since the tolerence',ch10,&
&     'is not on forces (force computation is expensive in wavelets).'
     MSG_WARNING(message)
   end if

   call timab(56,2,tsec)

!  ######################################################################
!  Compute the density rho from the trial potential
!  ----------------------------------------------------------------------

   call timab(242,1,tsec)
!  Compute the density from the trial potential
   if (dtset%tfkinfunc==0) then

     if(VERBOSE)then
       call wrtout(std_out,'*. Compute the density from the trial potential (vtorho)',"COLL")
     end if

     call vtorho(afford,atindx,atindx1,cg,compch_fft,cprj,cpus,dbl_nnsclo,&
&     dielop,dielstrt,dmatpawu,dphase,dtefield,dtfil,dtset,&
&     eigen,electronpositron,energies,etotal,gbound_diel,&
&     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&     istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mcprj,mgfftdiel,mpi_enreg,&
&     my_natom,dtset%natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&     npwarr,npwdiel,res2,psps%ntypat,nvresid,occ,&
&     computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
&     pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,fock,&
&     pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,rmet,rprimd,&
&     susmat,symrec,taug,taur,ucvol_local,usecprj,wffnew,vtrial,vxctau,wvl,xred,&
&     ylm,ylmgr,ylmdiel)

   else if (dtset%tfkinfunc==1.or.dtset%tfkinfunc==11.or.dtset%tfkinfunc==12) then
     MSG_WARNING('THOMAS FERMI')
     call vtorhotf(dtfil,dtset,energies%e_kinetic,energies%e_nonlocalpsp,&
&     energies%entropy,energies%e_fermie,gprimd,grnl,irrzon,mpi_enreg,&
&     dtset%natom,nfftf,dtset%nspden,dtset%nsppol,dtset%nsym,phnons,&
&     rhog,rhor,rprimd,ucvol,vtrial)
     residm=zero
     energies%e_eigenvalues=zero
   end if

!  Recursion method
   if(dtset%userec==1)then
     call vtorhorec(dtset,&
&     energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,energies%e_eigenvalues,&
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

    !if (dtset%useria == -4242) then
    !  call gshgg_mkncwrite(istep, dtset, dtfil, psps, hdr, pawtab, pawfgr, paw_ij, mpi_enreg, &
    !     rprimd, xred, eigen, npwarr, kg, ylm, ngfft, dtset%nfft, ngfftf, nfftf, vtrial) !electronpositron) ! Optional arguments
    !end if

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
&       taug=taug,taur=taur,vxctau=vxctau,add_tfw=tfw_activated)
     end if
!    If the density mixing is required, compute the total energy here
     call etotfor(atindx1,deltae,diffor,dtefield,dtset,&
&     elast,electronpositron,energies,&
&     etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,gresid,grewtn,grhf,grnl,grvdw,&
&     grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,my_natom,&
&     nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,psps%ntypat,nvresid,n1xccc,n3xccc,&
&     optene,computed_forces,optres,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&     ph1df,red_ptot,psps,rhog,rhor,rmet,rprimd,symrec,synlgr,ucvol,&
&     psps%usepaw,vhartr,vpsp,vxc,wvl%descr,wvl%den,xccc3d,xred)
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
&     dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode,electronpositron=electronpositron)
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
!      to be adjusted for the call to rhohxc
       nk3xc=1
       if(dtset%icoulomb==0 .and. dtset%usewvl==0) then
         call rhohxc(dtset,edum,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,0,nkxc,nk3xc,dtset%nspden,n3xccc,&
&         optxc,rhog,rhor,rprimd,dummy2,0,vhartr,vxc,vxcavg_dum,xccc3d,&
&         add_tfw=tfw_activated,taug=taug,taur=taur,vxctau=vxctau)
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
&     rprimd,susmat,psps%usepaw,vtrial,wvl%descr,wvl%den,xred)
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
   if(VERBOSE)then
     call wrtout(std_out,'*. Compute the new potential from the trial density',"COLL")
   end if

!  Set XC computation flag
   optxc=1
   if (nkxc>0) then
     if (dtset%iscf<0) optxc=2
     if (modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79).and. &
&     dtset%iscf<10.and. &
&     (dtset%iprcel>=100.or.istep==1.or.istep==dielstrt)) optxc=2
     if (dtset%iscf>=10.and.dtset%densfor_pred/=0.and.abs(dtset%densfor_pred)/=5) optxc=2
     if (optxc==2.and.dtset%xclevel==2.and.nkxc==3-2*mod(dtset%nspden,2)) optxc=12
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

     call rhotov(dtset,energies,gprimd,gsqcut,istep,kxc,mpi_enreg,nfftf,ngfftf, &
&     nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,optene,optres,optxc,&
&     rhog,rhor,rprimd,strsxc,ucvol_local,psps%usepaw,usexcnhat,&
&     vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,wvl,xccc3d,xred,&
&     electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau,add_tfw=tfw_activated)

   end if

   call timab(243,2,tsec)
   call timab(60,1,tsec)

!  This is inside the loop, its not equivalent to the line 1821
   if(moved_atm_inside==1) xred_old(:,:)=xred(:,:)

   if (dtset%iscf<10) then

     if(VERBOSE)then
       call wrtout(std_out,'Check exit criteria in case of potential mixing',"COLL")
     end if

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
&       hyb_mixing=hybrid_mixing,hyb_mixing_sr=hybrid_mixing_sr,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       electronpositron=electronpositron)

     end if

!    Add the Fock contribution to E_xc and E_xcdc if required
     if (usefock==1) then
       energies%e_fockdc=two*energies%e_fock
     end if

     if (.not.wvlbigdft) then
       call etotfor(atindx1,deltae,diffor,dtefield,dtset,&
&       elast,electronpositron,energies,&
&       etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,gresid,grewtn,grhf,grnl,grvdw,&
&       grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,my_natom,&
&       nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,dtset%ntypat,nvresid,n1xccc, &
&       n3xccc,0,computed_forces,optres,pawang,pawfgrtab,pawrad,pawrhoij,&
&       pawtab,ph1df,red_ptot,psps,rhog,rhor,rmet,rprimd,symrec,synlgr,ucvol,&
&       psps%usepaw,vhartr,vpsp,vxc,wvl%descr,wvl%den,xccc3d,xred)
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
&     dtfil%filnam_ds(1),initialized0,dtset%iscf,istep,dtset%kptns,&
&     maxfor,moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,nstep,&
&     occ,optres,prtfor,prtxml,quit,res2,resid,residm,response,tollist,&
&     psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode,electronpositron=electronpositron)
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

     if(VERBOSE)then
       call wrtout(std_out,'*. Mix the potential (if required) - Check exit criteria',"COLL")
     end if

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
&     vhartr,vnew_mean,vpsp,nvresid,vtrial,vxc,xred,&
&     nfftf,&
&     pawtab,rhog,wvl)

   end if   ! iscf<10

!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  ######################################################################

   if(VERBOSE)then
     call wrtout(std_out,'*. END MINIMIZATION ITERATIONS',"COLL")
   end if

!  The initialisation of the gstate run should be done when this point is reached
   initialized=1

   !if (dtset%useria == -4242) then
   !  call gshgg_mkncwrite(istep, dtset, dtfil, psps, hdr, pawtab, pawfgr, paw_ij, mpi_enreg, &
   !     rprimd, xred, eigen, npwarr, kg, ylm, ngfft, dtset%nfft, ngfftf, nfftf, vtrial) !electronpositron) ! Optional arguments
   !end if

!  This is to save the density for restart.
   if (iwrite_fftdatar(mpi_enreg)) then

     if(dtset%prtden<0.or.dtset%prtkden<0) then
!      Update the content of the header (evolving variables)
       bantot=hdr%bantot
       if (dtset%positron==0) then
         call hdr_update(hdr,bantot,etotal,energies%e_fermie,residm,&
&         rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1),&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       else
         call hdr_update(hdr,bantot,electronpositron%e0,energies%e_fermie,residm,&
&         rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1),&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
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
         dtset%nspden,rhor,mpi_enreg,eigen=eigen)
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
         dtset%nspden,taur,mpi_enreg,eigen=eigen)
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

 ! Avoid pending requests if itime == ntime.
 call xmpi_wait(quitsum_request,ierr)
 if (timelimit_exit == 1) istep = istep - 1

 call timab(246,1,tsec)

 if (dtset%iscf > 0) then
   call ab7_mixing_deallocate(mix)
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
   call energy(cg,compch_fft,dtset,electronpositron,&
&   energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,my_natom,&
&   nfftf,ngfftf,nhat,nhatgr,nhatgrdim,npwarr,n3xccc,&
&   occ,optene,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
&   phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,taug,taur,usexcnhat,&
&   vhartr,vtrial,vpsp,vxc,vxctau,wvl%wfs,wvl%descr,wvl%den,wvl%e,xccc3d,xred,ylm,&
&   add_tfw=tfw_activated)
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
!   call pawcprj_reorder(cprj,atindx1)
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
& dtset%kssform  ==3  .or. &
& dtset%pawfatbnd> 0  .or. &
& dtset%pawprtwf > 0  )

 if (recompute_cprj) then
   usecprj=1
   mband_cprj=dtset%mband/mpi_enreg%nproc_band
   mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
   call pawcprj_free(cprj)
   ABI_DATATYPE_DEALLOCATE(cprj) ! Was previously allocated (event if size = 0,0)
   ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
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
   call pawcprj_alloc(cprj,ncpgr,dimcprj)
   iatom=0 ; iorder_cprj=1 ! cprj are not ordered
   call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,&
&   iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&   dtset%mband,mcg,mcprj,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,&
&   dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,&
&   dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&   ucvol,dtfil%unpaw,xred,ylm,ylmgr)
 end if

 call timab(246,2,tsec)
 call timab(247,1,tsec)

!SHOULD CLEAN THE ARGS OF THIS ROUTINE
 call afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&
& deltae,diffor,dtefield,dtfil,dtset,eigen,electronpositron,elfr,&
& energies,etotal,favg,fcart,fock,forold,fred,grchempottn,&
& gresid,grewtn,grhf,grhor,grvdw,&
& grxc,gsqcut,hdr,indsym,irrzon,istep,kg,kxc,lrhor,maxfor,mcg,mcprj,mgfftf,&
& moved_atm_inside,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfft,ngfftf,ngrvdw,nhat,&
& nkxc,npwarr,nvresid,occ,optres,paw_an,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,ph1d,ph1df,phnons,pion,prtfor,&
& prtxml,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,&
& taur,tollist,usecprj,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,&
& xccc3d,xred,ylm,ylmgr,dtset%charge*SUM(vpotzero(:)),conv_retcode)

!Before leaving the present routine, save the current value of xred.
 xred_old(:,:)=xred(:,:)

 call timab(247,2,tsec)

!######################################################################
!All calculations in scfcv are finished. Printing section
!----------------------------------------------------------------------

 call timab(248,1,tsec)

 call outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dmatpawu,dtfil,&
& dtset,ecut,eigen,electronpositron,elfr,etotal,energies%e_fermie,&
& gmet,gprimd,grhor,hdr,kg,lrhor,dtset%mband,mcg,mcprj,dtset%mgfft,&
& dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,my_natom,dtset%natom,nattyp,&
& nfftf,ngfftf,nhat,dtset%nkpt,npwarr,dtset%nspden,&
& dtset%nsppol,dtset%nsym,psps%ntypat,n3xccc,occ,paw_dmft,pawang,pawfgr,pawfgrtab,&
& pawrad,pawrhoij,pawtab,paw_an,paw_ij,dtset%prtvol,psps,results_gs,&
& rhor,rprimd,taur,ucvol,usecprj,vhartr,vpsp,vtrial,vxc,wvl%den,xccc3d,xred)

 if(associated(elfr))then
   ABI_DEALLOCATE(elfr)
   nullify(elfr)
 end if

 if(associated(grhor))then
   ABI_DEALLOCATE(grhor)
   nullify(grhor)
 end if

 if(associated(lrhor))then
   ABI_DEALLOCATE(lrhor)
   nullify(lrhor)
 end if

 if(dtset%prtkden/=0 .or. dtset%prtelf/=0)then
   ABI_DEALLOCATE(taur)
 end if

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

!######################################################################
!Deallocate memory and save results
!----------------------------------------------------------------------

 call prc_mem_free()

 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(fred)
 ABI_DEALLOCATE(forold)
 ABI_DEALLOCATE(grchempottn)
 ABI_DEALLOCATE(gresid)
 ABI_DEALLOCATE(grewtn)
 ABI_DEALLOCATE(grnl)
 ABI_DEALLOCATE(grvdw)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(synlgr)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(vhartr)
 ABI_DEALLOCATE(vtrial)
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(vxc)
 ABI_DEALLOCATE(vxctau)
 ABI_DEALLOCATE(xccc3d)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(shiftvector)
 ABI_DEALLOCATE(dtn_pc)
 ABI_DEALLOCATE(grhf)
 ABI_DEALLOCATE(nvresid)

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
 !end if

 if (psps%usepaw==1) then
   if (dtset%iscf>0) then
     do iatom=1,my_natom
       pawrhoij(iatom)%lmnmix_sz=0
       pawrhoij(iatom)%use_rhoijres=0
       ABI_DEALLOCATE(pawrhoij(iatom)%kpawmix)
       ABI_DEALLOCATE(pawrhoij(iatom)%rhoijres)
     end do
   end if
   if (recompute_cprj.or.usecprj==1) then
     usecprj=0;mcprj=0
     if (recompute_cprj.or.dtset%mkmem/=0) then
       call pawcprj_free(cprj)
     end if
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
 ABI_DEALLOCATE(dimcprj_srt)
 ABI_DEALLOCATE(dimcprj)
 call pawcprj_free(cprj)
 ABI_DATATYPE_DEALLOCATE(cprj)

! Deallocate exact exchange data at the end of the calculation
 if (usefock==1) then
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

 call status(0,dtfil%filstat,iexit,level,'exit')

 call timab(249,2,tsec)
 call timab(238,2,tsec)

 DBG_EXIT("COLL")

end subroutine scfcv
!!***

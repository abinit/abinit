!!****m* ABINIT/m_ksdiago
!! NAME
!!  m_ksdiago
!!
!! FUNCTION
!!  Direct diagonalization of the KS Hamiltonian H_k(G,G')
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ksdiago

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_hamiltonian
 use m_distribfft
 use libxc_functionals
 use m_ebands
 use m_nctk
 use m_dtfil
 use m_hdr
 use m_wfk

 use defs_datatypes,      only : pseudopotential_type
 use defs_abitypes,       only : MPI_type
 use m_gwdefs,            only : GW_TOLQ0, GW_Q0_DEFAULT !, cone_gw, czero_gw, j_gw
 use m_dtset,             only : dataset_type
 use m_fstrings,          only : toupper, ktoa, itoa, sjoin, ftoa, ltoa
 use m_io_tools,          only : iomode_from_fname, get_unit
 use m_yaml,              only : yamldoc_t, yamldoc_open
 use m_numeric_tools,     only : blocked_loop
 use m_time,              only : cwtime, cwtime_report, timab
 use m_geometry,          only : metric, normv
 use m_hide_lapack,       only : xhegv_cplex, xheev_cplex, xheevx_cplex, xhegvx_cplex
 use m_slk,               only : slkmat_dp_t, slk_processor_t, block_dist_1d, &
                                 compute_eigen_problem, compute_generalized_eigen_problem
 use m_bz_mesh,           only : findnq, findq, findqg0, identk
 use m_kg,                only : mkkin, mkkpg
 use m_crystal,           only : crystal_t
 use m_fftcore,           only : kpgsph, get_kg
 use m_fft_mesh,          only : calc_ceigr, get_gfft
 use m_fft,               only : fftpac, uplan_t, fftbox_plan3_t, zerosym
 use m_cgtools,           only : set_istwfk
 use m_electronpositron,  only : electronpositron_type
 use m_mpinfo,            only : destroy_mpi_enreg, initmpi_seq
 use m_pawtab,            only : pawtab_type
 use m_paw_ij,            only : paw_ij_type
 use m_pawcprj,           only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_reorder, &
                                 pawcprj_set_zero, pawcprj_mpi_sum, pawcprj_copy
 use m_cgprj,             only : getcprj
 use m_pawfgr,            only : pawfgr_type
 use m_initylmg,          only : initylmg
 use m_mkffnl,            only : mkffnl
 use m_getghc,            only : getghc, multithreaded_getghc
 use m_wfd,               only : wfd_t
 use m_vcoul,             only : vcgen_t
 use m_occ,               only : get_fact_spin_tol_empty
 use m_pstat,             only : pstat_proc

 implicit none

 private
!!***

!!****t* m_ksdiago/ugb_t
!! NAME
!!  ugb_t
!!
!! FUNCTION
!!  This object stores the wavefunctions for given (k-point, spin) in a
!!  PBLAS matrix distributed over bands (columns)
!!
!! SOURCE

 type, public :: ugb_t

   integer :: istwf_k = -1
   ! Storage mode of cg_k.

   integer :: nspinor = -1
   ! Number of spinors.

   integer :: npw_k = -1
   ! Number of planewaves.

   integer :: npwsp = -1
   ! nnpw_k * nspinor.

   integer :: nband_k = - 1
   ! Total number of bands (global)

   integer :: my_bstart = -1, my_bstop = - 1, my_nband = - 1
   ! 1) Initial band
   ! 2) Last band
   ! 3) Number of bands treated by this proc. 0 if idle proc.

   logical :: has_idle_procs
   ! True if there are procs in comm who don't own any column.

   integer,pointer :: comm
   ! pointer to MPI communicator in mat

   type(slk_processor_t) :: processor

   type(slkmat_dp_t) :: mat
   ! PBLAS matrix with MPI-distributed Fourier components (double precision)
   ! Local buffer: (2, npwsp * my_nband)
   ! Global matrix: (npwsp, nband_k)

   integer, allocatable :: kg_k(:,:)
   ! (3, npw_k)
   ! G-vectors in reduced coordinates.

   real(dp), contiguous, pointer :: cg_k(:,:,:)
   ! (2, npwsp * my_nband)
   ! NB: This is a pointer to mat%buffer_cplx

   type(pawcprj_type),allocatable :: cprj_k(:,:)
   ! (natom, nspinor * my_nband))
   ! PAW projections ordered according to natom and NOT according to typat.
   ! NOTE my_nband

 contains

   procedure :: from_diago  => ugb_from_diago
    ! Build object by direct diagonalization of the KS Hamiltonian.

   procedure :: from_wfk_file  => ugb_from_wfk_file
    ! Build object from WFK file.

   procedure :: free => ugb_free
    ! Free memory.

   procedure :: print => ugb_print
    ! Print info on object.

   procedure :: collect_cprj => ugb_collect_cprj
    ! Collect a subset of PAW cprj on all processors.

 end type ugb_t
!!***

!!****t* m_ksdiago/hyb_t
!! NAME
!!  hyb_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: hyb_t

   integer :: nkibz = -1, nkbz = -1
   integer :: nqibz = -1, nqbz = -1

   integer :: mg0(3) = [2, 2, 2]
   ! Max shifts to account for umklapps.

   type(wfd_t) :: wfd
   type(vcgen_t) :: vcgen
   type(ebands_t) :: ebands

   real(dp),allocatable :: kibz(:,:), kbz(:,:)
   real(dp),allocatable :: qibz(:,:), qbz(:,:), wtq(:)

   integer,allocatable :: kbz2ibz(:,:)
    ! kbz2ibz(6, nkbz))

   integer,allocatable :: kbz2ibz_symrel(:,:)
    ! kbz2ibz_symrel(6, nkbz))

   integer,allocatable :: qbz2ibz(:,:)
    ! qbz2ibz(6, nqbz))

 contains

   procedure :: from_wfk_file  => hyb_from_wfk_file
    ! Build object from WFK file

   !procedure :: print => hyb_print
    ! Print info on object.

   procedure :: free => hyb_free
    ! Free memory.

 end type hyb_t
!!***

!!****t* m_ksdiago/ddiago_ctl_type
!! NAME
!!  ddiago_ctl_type
!!
!! FUNCTION
!!  Structure storing the variables controlling the direct diagonalization of the Kohn-Sham Hamiltonian.
!!  Mainly used for debugging (and in the KSS code!)
!!
!! SOURCE

 type, public :: ddiago_ctl_type

  integer :: spin
   ! The spin component of the Hamiltonian (1 if nspinor==1 or nsppol==1).

  integer :: istwf_k
   ! Option defining whether time-reversal symmetry is used at particular k-points
   ! If 0, the code will automatically use TR symmetry if possible (depending on the k-point)

  integer :: nband_k
   ! Number of bands to be calculated.

  integer :: npw_k
  ! The number of planes waves for the wavefunctions taking into account time-reversal symmetry.

  integer :: npwtot
  ! The number of planes waves in the Hamiltonian without taking into account istwf_k

  integer :: nspinor
  ! Number of spinorial components.

  integer :: prtvol
   ! Flag controlling the verbosity.

  integer :: use_scalapack
  ! 0 if diagonalization is done in sequential on each node.
  ! 1 to use scalapack
  ! TODO Not implemented

  real(dp) :: abstol
   ! used fro RANGE= "V", "I", and "A" when do_full_diago=.FALSE.
   ! The absolute error tolerance for the eigenvalues. An approximate eigenvalue is accepted
   ! as converged when it is determined to lie in an interval [a,b] of width less than or equal to
   !
   !         ABSTOL + EPS *   max( |a|,|b| ) ,
   !
   ! where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then  EPS*|T|  will be used in its place,
   ! where |T| is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.
   !
   ! Eigenvalues will be computed most accurately when ABSTOL is
   ! set to twice the underflow threshold 2*DLAMCH('S'), not zero.
   ! If this routine returns with INFO>0, indicating that some
   ! eigenvectors did not converge, try setting ABSTOL to 2*DLAMCH('S').

  real(dp) :: ecut
   ! The cutoff energy for the plane wave basis set.

  real(dp) :: ecutsm
   ! Smearing energy for plane wave kinetic energy (Ha)

  real(dp) :: effmass_free
   ! Effective mass for electrons (usually one).

  logical :: do_full_diago
  ! Specifies whether direct or partial diagonalization will be performed.
  ! Meaningful only if RANGE='A'.

  integer :: ilu(2)
   ! If RANGE='I', the indices (in ascending order) of the smallest and largest eigenvalues to be returned.
   ! il=ilu(1), iu=ilu(2) where
   ! 1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. NOT used if RANGE = 'A' or 'V'.

  integer :: nloalg(3)

  real(dp) :: kpoint(3)
   ! The k-point in reduced coordinates at which the Hamiltonian is diagonalized.

  real(dp) :: vlu(2)
   ! If RANGE='V', the lower and upper bounds of the interval to
   ! be searched for eigenvalues. vl=vlu(1) and vu=vlu(2) with VL < VU.
   ! Not referenced if RANGE = 'A' or 'I'.

  character(len=1) :: jobz
   ! character defining whether wavefunctions are required (lapack option).
   ! "N":  Compute eigenvalues only;
   ! "V":  Compute eigenvalues and eigenvectors.

  character(len=1) :: range
   ! character defining the subset of eigenstates that will be calculated (lapack option).
   ! "A": all eigenvalues will be found.
   ! "V": all eigenvalues in the half-open interval (VL,VU] will be found.
   ! "I": the IL-th through IU-th eigenvalues will be found.

  !$character(len=fnlen) :: fname
  ! The name of the file storing the eigenvectors and eigenvalues (only if jobz="V")

 end type ddiago_ctl_type
!!***


!!****t* m_ksdiago/psbands_t
!! NAME
!!  psbands_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: psbands_t

   integer :: nb_tot = -1
   ! Total number of states (protected + pseudo bands)

   integer :: nb_protected = -1
   ! Number of protected bands.

   integer :: nslices = -1
   ! Number of slices.

   integer :: maxsto_per_slice = -1
   ! Max number of pseudo bands per slice.

   real(dp) :: efrac

   integer,allocatable :: subspace(:,:)
   ! (3, nslices)
   ! For each slice, the first and last band index and the number of pseudo bands in the slice.

   real(dp),allocatable :: ps_eig(:)
   ! (nb_tot)
   ! eigenvalues (KS + pseudo energies)

 contains

   procedure :: init => psbands_init
    ! Initialize the object

   procedure :: band2slice => psbands_band2slice
   ! Return the slice index from the band index.

   procedure :: free => psbands_free
    ! Free memory.

 end type psbands_t
!!***

 public :: ksdiago
 public :: init_ddiago_ctl

!!***

contains
!!***

!!****f* m_ksdiago/ksdiago
!! NAME
!! ksdiago
!!
!! FUNCTION
!!  This routine performs the direct diagonalization of the Kohn-Sham Hamiltonian
!!  for a given k-point and spin. The routine drives the following operations:
!!
!!    1) Re-computing <G|H|G_prim> matrix elements for all (G, G_prim).
!!       starting from the knowledge of the local potential on the real-space FFT mesh.
!!
!!    2) Diagonalizing H in the plane-wave basis.
!!
!!  It is called in outkss.F90 during the generation of the KSS file
!!  needed for a GW post-treatment. Since many-body calculations usually
!!  require a large number of eigenstates eigen-functions, a direct
!!  diagonalization of the Hamiltonian might reveal more stable than iterative
!!  techniques that might be problematic when several high energy states are required.
!!  The main drawback of the direct diagonalization is the bad scaling with the size
!!  of the basis set (npw**3) and the large memory requirements.
!!
!! INPUTS
!!  kpoint(3)
!!  prtvol=Integer Flags  defining verbosity level
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  mgfftc=maximum size of 1D FFTs (coarse mesh).
!!  natom=number of atoms in cell.
!!  nfftf=(effective) number of FFT grid points in the dense FFT mesh (for this processor)
!!         (nfftf=nfft for norm-conserving potential runs)
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=number of density components
!!  pawtab(psps%ntypat*psps%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pawfgr<pawfgr_type>=fine grid parameters and related data
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  vtrial(nfftf,nspden)=the trial potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  comm=MPI communicator.
!!  [electronpositron] <electronpositron_type>=quantities for the electron-positron annihilation.
!!  nfftc=Number of points in the coarse FFT mesh.
!!  ngfftc(18)=Info about 3D FFT for the coarse mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  Diago_ctl<ddiago_ctl_type>=Datatype storing variables and options controlling the direct diagonalization.
!!
!! OUTPUT
!!  ierr=Status error.
!!  onband_diago
!!
!! SIDE EFFECTS
!!  eig_ene(:)=Pointer used for allocating and storing the eigenvalues (hartree)
!!    input: pointer to NULL
!!     output: eig_ene(onband_diago)=The calculatated eigenvalues in ascending order.
!!
!!  eig_vec(:,:,:)=Pointer used for allocating and holding the wave functions at this k-point and spin.
!!    input: pointer to NULL
!!    output: eig_vec(2,npw_k*nspinor,onband_diago)=The calculated eigenvectors.
!!
!!  cprj_k(natom,nspinor*onband_diago) PAW only===
!!   input: pointer to NULL
!!   output: Projected eigenstates <Proj_i|Cnk> from output eigenstates.
!!
!! NOTES
!! * The routine can be time consuming (in particular when computing <G1|H|G2> elements for all (G1,G2)).
!!   So, it is recommended to call it once per run.
!!
!! * The routine RE-compute all Hamiltonian terms. So it is equivalent to an additional electronic SCF cycle.
!!   (This has no effect is convergence was reached.
!!   If not, eigenvalues/vectors may differs from the conjugate gradient ones)
!!
!! * Please, do NOT pass Dtset% to this routine. Either use a local variable properly initialized
!!   or add the additional variable to ddiago_ctl_type and change the creation method accordingly.
!!   ksdiago is designed such that it is possible to diagonalize the Hamiltonian at an arbitrary k-point
!!   or spin (not efficient but easy to code). Therefore ksdiago is useful non only for
!!   the KSS generation but also for testing more advanced iterative algorithms as well as interpolation techniques.
!!
!! SOURCE

subroutine ksdiago(Diago_ctl, nband_k, nfftc, mgfftc, ngfftc, natom, &
                   typat, nfftf, nspinor, nspden, nsppol, pawtab, pawfgr, paw_ij,&
                   psps, rprimd, vtrial, xred, onband_diago, eig_ene, eig_vec, cprj_k, comm, ierr,&
                   electronpositron) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfftc,natom,comm,nband_k,nfftf,nsppol,nspden,nspinor,nfftc
 integer,intent(out) :: ierr, onband_diago
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ddiago_ctl_type),intent(in) :: Diago_ctl
!arrays
 integer,intent(in) :: typat(natom), ngfftc(18)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: vtrial(nfftf,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),pointer :: eig_ene(:),eig_vec(:,:,:)
 type(pawcprj_type),pointer :: cprj_k(:,:)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
 type(electronpositron_type),optional,pointer :: Electronpositron

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem1 = 1, tim_getghc = 4, paral_kgb0 = 0, master = 0, ndat1 = 1, ncomp1 = 1
 integer :: cprj_choice,cpopt,dimffnl,ib,ider,idir,spin,npw_k
 integer :: ikg,istwf_k,exchn2n3d,prtvol
 integer :: jj,n1,n2,n3,n4,n5,n6,negv,nkpg,nproc,npw_k_test,my_rank,optder
 integer :: type_calc,sij_opt,igsp2,cplex_ghg,iband,ibs1,ibs2
 real(dp),parameter :: lambda0 = zero
 real(dp) :: ucvol,ecutsm,effmass_free,size_mat,ecut
 logical :: do_full_diago
 character(len=50) :: jobz,range
 character(len=80) :: frmt1,frmt2
 character(len=10) :: stag(2)
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg_seq
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer :: nloalg(3)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),kptns_(3,1),kpoint(3),ylmgr_dum(1,1,1)
 real(dp),allocatable :: ph3d(:,:,:),bras(:,:),ffnl(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),ylm_k(:,:),dum_ylm_gr_k(:,:,:), vxctaulocal(:,:,:,:,:)
 real(dp),allocatable :: ghc(:,:),gvnlxc(:,:),gsc(:,:),ghg_mat(:,:,:),gsg_mat(:,:,:)
 real(dp),pointer :: cwavef(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
! *********************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 if (nproc > 1) then
   ABI_WARNING("ksdiago not supported in parallel. Running in sequential.")
 end if

 call initmpi_seq(mpi_enreg_seq) ! Fake MPI_type for sequential part.
 call init_distribfft_seq(mpi_enreg_seq%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 if (pawfgr%usefinegrid /= 0) then
   call init_distribfft_seq(mpi_enreg_seq%distribfft, 'f', pawfgr%ngfft(2), pawfgr%ngfft(3), 'all')
 end if

 spin  = Diago_ctl%spin
 kpoint  = Diago_ctl%kpoint
 istwf_k = Diago_ctl%istwf_k
 !% nband_k = Diago_ctl%nband_k
 npw_k   = Diago_ctl%npw_k
 nloalg  = Diago_ctl%nloalg
 ecut    = Diago_ctl%ecut
 ecutsm  = Diago_ctl%ecutsm
 effmass_free = Diago_ctl%effmass_free
 prtvol  = Diago_ctl%prtvol

 call metric(gmet, gprimd, -1, rmet, rprimd, ucvol)

 if (nsppol == 1) stag = ['          ','          ']
 if (nsppol == 2) stag = ['SPIN UP:  ','SPIN DOWN:']

 ! The coarse FFT mesh.
 n1 = ngfftc(1); n2 = ngfftc(2); n3 = ngfftc(3)
 n4 = ngfftc(4); n5 = ngfftc(5); n6 = ngfftc(6)

 !====================
 !=== Check input ====
 !====================
 ierr = 0

 ! istwfk must be 1 for each k-point
 if (istwf_k/=1) then
   write(msg,'(7a)')&
   ' istwfk /= 1 not allowed:',ch10,&
   ' States output not programmed for time-reversal symmetry.',ch10,&
   ' Action: change istwfk in input file (put it to 1 for all kpt).',ch10,&
   ' Program does not stop but _KSS file will not be created...'
   ABI_WARNING(msg)
   ierr = ierr + 1
 end if

 if (ierr /= 0) RETURN ! Houston we have a problem!

 ! Initialize the Hamiltonian datatype on the coarse FFT mesh.
 if (present(electronpositron)) then
   call gs_hamk%init(psps, pawtab, nspinor, nsppol, nspden, natom, typat, xred, nfftc, &
    mgfftc, ngfftc, rprimd, nloalg, paw_ij=paw_ij, usecprj=0, electronpositron=electronpositron)
 else
   call gs_hamk%init(psps, pawtab, nspinor, nsppol, nspden, natom, typat, xred, nfftc, &
    mgfftc, ngfftc, rprimd, nloalg, paw_ij=paw_ij, usecprj=0)
 end if

 ! Check on the number of stored bands.
 onband_diago = nband_k
 if (nband_k==-1 .or. nband_k >= npw_k*nspinor) then
   onband_diago = npw_k*nspinor
   write(msg,'(4a,i0)')ch10,&
    ' Since the number of bands to be computed was -1 or',ch10,&
    ' too large, it has been set to the maximum value npw_k*nspinor: ',npw_k*nspinor
   call wrtout(std_out, msg)
 end if

 !do_full_diago = (onband_diago==npw_k*nspinor)
 do_full_diago = Diago_ctl%do_full_diago

 if (do_full_diago) then
   write(msg,'(6a)')ch10,&
   ' Since the number of bands to be computed',ch10,&
   ' is equal to the number of G-vectors found for this kpt,',ch10,&
   ' the program will perform complete diagonalization.'
 else
   write(msg,'(6a)')ch10,&
   ' Since the number of bands to be computed',ch10,&
   ' is less than the number of G-vectors found,',ch10,&
   ' the program will perform partial diagonalization.'
 end if
 if (prtvol > 0) call wrtout(std_out, msg)

 ! Set up local potential vlocal with proper dimensioning, from vtrial.
 ! Select spin component of interest if nspden<=2 as nvloc==1, for nspden==4, nvloc==4
 ! option=2: vtrial(n1*n2*n3,ispden) --> vlocal(nd1,nd2,nd3) real case

 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamk%nvloc))
 !if (with_vxctau) then
 !  ABI_MALLOC(vxctaulocal,(n4,n5,n6,gs_hamk%nvloc,4))
 !end if

 ! Set up local potential vlocal on the coarse FFT mesh from vtrial taking into account the spin.

 call gspot_transgrid_and_pack(spin, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
                               nspden, gs_hamk%nvloc, ncomp1, pawfgr, mpi_enreg_seq, vtrial, vlocal)
 call gs_hamk%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

 ! This for meta-gga.
 !if (with_vxctau) then
 !  call gspot_transgrid_and_pack(spin, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
 !                                nspden, gs_hamk%nvloc, 4, pawfgr, mpi_enreg, vxctau, vxctaulocal)
 !  call gs_hamk%load_spin(spin, vxctaulocal=vxctaulocal)
 !end if

 ! Calculate G-vectors, for this k-point. Count also the number of planewaves as a check.
 exchn2n3d = 0; ikg = 0
 ABI_MALLOC(kg_k, (3, npw_k))

 call kpgsph(ecut, exchn2n3d, gmet, ikg, 0, istwf_k, kg_k, kpoint, 0, mpi_enreg_seq, 0, npw_k_test)
 ABI_CHECK(npw_k_test == npw_k, "npw_k_test/=npw_k")
 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,mkmem1,mpi_enreg_seq,npw_k,npw_k_test)

 !========================
 !==== Kinetic energy ====
 !========================
 ABI_MALLOC(kinpw, (npw_k))
 call mkkin(ecut, ecutsm, effmass_free, gmet, kg_k, kinpw, kpoint, npw_k, 0, 0)

 !================================
 !==== Non-local form factors ====
 !================================
 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2*psps%useylm))

 if (psps%useylm == 1) then
   optder = 0
   ABI_MALLOC(dum_ylm_gr_k, (npw_k, 3+6*(optder/2),psps%mpsang**2))
   kptns_(:,1) = kpoint

   ! Here mband is not used if paral_compil_kpt=0
   call initylmg(gprimd, kg_k, kptns_, mkmem1, mpi_enreg_seq, psps%mpsang, npw_k, [nband_k], 1, &
     [npw_k], 1, optder, rprimd, ylm_k, dum_ylm_gr_k)

   ABI_FREE(dum_ylm_gr_k)
 end if

 ! Compute (k+G) vectors (only if useylm=1)
 nkpg = 3 * nloalg(3)
 ABI_MALLOC(kpg_k, (npw_k, nkpg))
 if (nkpg > 0) call mkkpg(kg_k, kpg_k, kpoint, nkpg, npw_k)

 ! Compute nonlocal form factors ffnl at all (k+G):
 idir=0; ider=0; dimffnl=1+ider ! Now the derivative is not needed anymore.
 ABI_MALLOC(ffnl, (npw_k, dimffnl, psps%lmnmax, psps%ntypat))

 call mkffnl(psps%dimekb, dimffnl, psps%ekb, ffnl, psps%ffspl, gmet, gprimd, ider, idir, psps%indlmn, &
   kg_k, kpg_k, kpoint, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, &
   psps%ntypat, psps%pspso, psps%qgrid_ff, rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum)

 ABI_FREE(ylm_k)

 ! Load k-dependent part in the Hamiltonian datastructure
 ABI_MALLOC(ph3d, (2, npw_k, gs_hamk%matblk))
 call gs_hamk%load_k(kpt_k=kpoint, istwf_k=istwf_k, npw_k=npw_k, kinpw_k=kinpw, &
                     kg_k=kg_k, kpg_k=kpg_k, ffnl_k=ffnl, ph3d_k=ph3d, compute_ph3d=.true., compute_gbound=.true.)

 ! Prepare call to getghc.
 type_calc = 0                                ! For applying the whole Hamiltonian
 sij_opt = 0; if (psps%usepaw==1) sij_opt = 1 ! For PAW, <k+G|S|k+G"> is also needed.

 cpopt = -1    ! If cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
 if (psps%usepaw==1.and..FALSE.) then ! TODO Calculate <p_lmn|k+G>.
   cpopt = 0  ! <p_lmn|in> are computed here and saved
 end if

 ABI_MALLOC(ghc, (2, npw_k*nspinor*ndat1))
 ABI_MALLOC(gvnlxc, (2, npw_k*nspinor*ndat1))
 ABI_MALLOC(gsc, (2, npw_k*nspinor*ndat1*(sij_opt+1)/2))

 cplex_ghg = 2
 size_mat = cplex_ghg*(npw_k*nspinor)**2*dp*b2Mb
 write(msg,'(a,f0.3,a)')" Out-of-memory in ghg_mat. Memory required by the Hamiltonian matrix: ",size_mat," [Mb]."
 ABI_STAT_MALLOC(ghg_mat, (cplex_ghg, npw_k*nspinor, npw_k*nspinor), ierr)
 ABI_CHECK(ierr == 0, msg)
 write(msg,'(a,f0.3,a)')" Out-of-memory in gsg_mat. Memory required by the PAW overlap operator: ",size_mat," [Mb]."
 ABI_STAT_MALLOC(gsg_mat, (cplex_ghg, npw_k*nspinor, npw_k*nspinor*psps%usepaw), ierr)
 ABI_CHECK(ierr == 0, msg)

 ! cwaveprj is ordered by atom type, see nonlop_ylm.
 ABI_MALLOC(cwaveprj, (natom, nspinor*(1+cpopt)*gs_hamk%usepaw))
 if (cpopt == 0) call pawcprj_alloc(cwaveprj, 0, gs_hamk%dimcprj)

 ! Initialize plane-wave array with zeros
 ABI_CALLOC(bras, (2, npw_k*nspinor))
 if (prtvol > 0) call wrtout(std_out, ' Calculating <G|H|G''> elements')

 ! Loop over the |beta,G''> component.
 do igsp2=1,npw_k*nspinor
   bras(1, igsp2) = one

   ! Get <:|H|beta,G''> and <:|S_{PAW}|beta,G''>
   call getghc(cpopt, bras, cwaveprj, ghc, gsc, gs_hamk, gvnlxc, lambda0, mpi_enreg_seq, ndat1, &
               prtvol, sij_opt, tim_getghc, type_calc)

   ! Fill the upper triangle.
   ghg_mat(:,1:igsp2,igsp2) = ghc(:,1:igsp2)
   if (psps%usepaw == 1) gsg_mat(:,1:igsp2,igsp2) = gsc(:,1:igsp2)

   ! Reset the |G,beta> component that has been treated.
   bras(1, igsp2) = zero
 end do

 ! Free workspace memory allocated so far.
 ABI_FREE(bras)
 ABI_FREE(kinpw)
 ABI_FREE(vlocal)
 ABI_FREE(ghc)
 ABI_FREE(gvnlxc)
 ABI_FREE(gsc)
 ABI_SFREE(vxctaulocal)

 if (psps%usepaw == 1 .and. cpopt == 0) call pawcprj_free(Cwaveprj)
 ABI_FREE(cwaveprj)

 !===========================================
 !=== Diagonalization of <G|H|G''> matrix ===
 !===========================================
 ABI_MALLOC(eig_ene, (onband_diago))
 ABI_MALLOC(eig_vec, (cplex_ghg, npw_k*nspinor, onband_diago))

 jobz = Diago_ctl%jobz  !jobz="Vectors"

 if (do_full_diago) then
   ! Full diagonalization
   write(msg,'(6a,i0)')ch10,&
     ' Begin full diagonalization for kpt: ',trim(ktoa(kpoint)), stag(spin), ch10,&
     ' Matrix size: ', npw_k*nspinor
   call wrtout(std_out, msg)

   if (psps%usepaw == 0) then
     call xheev_cplex(jobz, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, eig_ene, msg, ierr)
   else
     call xhegv_cplex(1, jobz, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, gsg_mat, eig_ene, msg, ierr)
   end if
   ABI_CHECK(ierr == 0, msg)
   eig_vec(:,:,:)=  ghg_mat

 else
   ! Partial diagonalization
   range = Diago_ctl%range !range="Irange"

   write(msg,'(2a,3es16.8,3a,i0,a,i0)')ch10,&
     ' Begin partial diagonalization for kpt= ',kpoint, stag(spin),ch10,&
     ' - Size of mat.=',npw_k*nspinor,' - # out_nband: ',onband_diago
   call wrtout(std_out, msg)

   if (psps%usepaw == 0) then
     call xheevx_cplex(jobz, range, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, zero, zero,&
       1, onband_diago, -tol8, negv, eig_ene, eig_vec, npw_k*nspinor, msg, ierr)
   else
     call xhegvx_cplex(1, jobz, range, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, gsg_mat, zero, zero,&
       1, onband_diago, -tol8, negv, eig_ene, eig_vec, npw_k*nspinor, msg, ierr)
   end if
   ABI_CHECK(ierr == 0, msg)
 end if

 ABI_FREE(ghg_mat)
 ABI_FREE(gsg_mat)

 if (prtvol > 0 .and. my_rank == master) then
   ! Write eigenvalues.
   frmt1 = '(8x,9(1x,f7.2))'; frmt2 = '(8x,9(1x,f7.2))'
   write(msg,'(2a,3x,a)')' Eigenvalues in eV for kpt: ', trim(ktoa(kpoint)), stag(spin)
   call wrtout(std_out, msg)

   write(msg,frmt1)(eig_ene(ib)*Ha_eV,ib=1,MIN(9,onband_diago))
   call wrtout(std_out, msg)
   if (onband_diago >9 ) then
     do jj=10,onband_diago,9
       write(msg, frmt2) (eig_ene(ib)*Ha_eV,ib=jj,MIN(jj+8,onband_diago)); call wrtout(std_out, msg)
     end do
   end if
 end if

 !========================================================
 !==== Calculate <Proj_i|Cnk> from output eigenstates ====
 !========================================================
 if (psps%usepaw == 1) then

   ABI_MALLOC(cprj_k,(natom, nspinor*onband_diago))
   call pawcprj_alloc(cprj_k, 0, gs_hamk%dimcprj)

   idir = 0; cprj_choice = 1  ! Only projected wave functions.

   do iband=1,onband_diago
     ibs1 = nspinor * (iband - 1) + 1
     ibs2 = ibs1; if (nspinor == 2) ibs2=ibs2+1
     cwavef => eig_vec(1:2,1:npw_k,iband)

     call getcprj(cprj_choice, 0, cwavef, cprj_k(:,ibs1:ibs2), &
       gs_hamk%ffnl_k, idir, gs_hamk%indlmn, gs_hamk%istwf_k, gs_hamk%kg_k, &
       gs_hamk%kpg_k, gs_hamk%kpt_k, gs_hamk%lmnmax, gs_hamk%mgfft, mpi_enreg_seq, 1, &
       gs_hamk%natom, gs_hamk%nattyp, gs_hamk%ngfft, gs_hamk%nloalg, gs_hamk%npw_k, gs_hamk%nspinor, &
       gs_hamk%ntypat, gs_hamk%phkxred, gs_hamk%ph1d, gs_hamk%ph3d_k, gs_hamk%ucvol, gs_hamk%useylm)
   end do

   !  Reorder the cprj (order is now the same as in input file)
   call pawcprj_reorder(cprj_k, gs_hamk%atindx1)
 end if ! usepaw

 ! Free memory.
 ABI_FREE(kpg_k)
 ABI_FREE(kg_k)
 ABI_FREE(ph3d)
 ABI_FREE(ffnl)

 call destroy_mpi_enreg(mpi_enreg_seq)
 call gs_hamk%free()
 call xmpi_barrier(comm)

end subroutine ksdiago
!!***
!----------------------------------------------------------------------

!!****f* m_ksdiago/init_ddiago_ctl
!! NAME
!!  init_ddiago_ctl
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine init_ddiago_ctl(Dctl, jobz, spin, nspinor, ecut, kpoint, nloalg, gmet, &
  nband_k, istwf_k, ecutsm, effmass_free, abstol, range, ilu, vlu, use_scalapack, prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,nspinor
 integer,optional,intent(in) :: istwf_k,prtvol,use_scalapack,nband_k
 real(dp),intent(in) :: ecut
 real(dp),optional,intent(in) :: ecutsm,effmass_free
 real(dp),optional,intent(in) :: abstol
 character(len=*),intent(in) :: jobz
 character(len=*),optional,intent(in) :: range
 type(ddiago_ctl_type),intent(out) :: Dctl
!arrays
 integer,intent(in) :: nloalg(3)
 integer,optional,intent(in) :: ilu(2)
 real(dp),intent(in) :: kpoint(3)
 real(dp),optional,intent(in) :: vlu(2)
 real(dp),intent(in) :: gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: npw_k
 logical :: ltest
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg_seq
!arrays
 integer,allocatable :: kg_k(:,:)
! *************************************************************************

 call initmpi_seq(mpi_enreg_seq) ! Fake MPI_type.

 Dctl%spin  = spin
 Dctl%nspinor = nspinor
 Dctl%kpoint  = kpoint

 if (PRESENT(istwf_k)) then
  Dctl%istwf_k = istwf_k
 else
  Dctl%istwf_k = set_istwfk(kpoint)
 end if

 ABI_CHECK(Dctl%istwf_k == 1, "istwf_k/=1 not coded")

 Dctl%jobz   = toupper(jobz(1:1))
 Dctl%range  = "A"
 if (PRESENT(range)) Dctl%range = toupper(range)

 Dctl%ecut = ecut
 Dctl%ecutsm = zero; if (PRESENT(ecutsm)) Dctl%ecutsm = ecutsm
 Dctl%effmass_free = one; if (PRESENT(effmass_free)) Dctl%effmass_free = effmass_free
 Dctl%nloalg  = nloalg
 Dctl%prtvol = 0; if (PRESENT(prtvol)) Dctl%prtvol = prtvol
 Dctl%abstol = -tol8; if (PRESENT(abstol)) Dctl%abstol = abstol

 ABI_MALLOC(kg_k,(3,0))

 ! Total number of G-vectors for this k-point with istwf_k=1.
 call kpgsph(ecut,0,gmet,0,0,1,kg_k,kpoint,0,mpi_enreg_seq,0,Dctl%npwtot)

 ! G-vectors taking into account time-reversal symmetry.
 call kpgsph(ecut,0,gmet,0,0,istwf_k,kg_k,kpoint,0,mpi_enreg_seq,0,npw_k)

 Dctl%npw_k = npw_k
 ABI_FREE(kg_k)

 Dctl%do_full_diago = .FALSE.

 SELECT CASE (Dctl%range)
 CASE ("A")

  ! Check on the number of stored bands.
  Dctl%nband_k=-1
  if (PRESENT(nband_k)) Dctl%nband_k=nband_k

  if (Dctl%nband_k==-1.or.Dctl%nband_k>=npw_k*nspinor) then
    Dctl%nband_k=npw_k*nspinor
    write(msg,'(4a)')ch10,&
    'Since the number of bands to be computed was (-1) or',ch10,&
    'too large, it has been set to the max. value npw_k*nspinor. '
    if (Dctl%prtvol>0) call wrtout(std_out, msg)
  else
    Dctl%nband_k=nband_k
  end if

  Dctl%do_full_diago = (Dctl%nband_k==npw_k*nspinor)

  if (Dctl%do_full_diago) then
    write(msg,'(6a)')ch10,&
     'Since the number of bands to be computed',ch10,&
     'is equal to the number of G-vectors found for this k-point,',ch10,&
     'the program will perform complete diagonalization.'
  else
    write(msg,'(6a)')ch10,&
      'Since the number of bands to be computed',ch10,&
      'is less than the number of G-vectors found,',ch10,&
      'the program will perform partial diagonalization.'
  end if
  if (Dctl%prtvol>0) call wrtout(std_out, msg)

 CASE ("I")
  if (.not.PRESENT(ilu)) then
    ABI_ERROR(" ilu must be specified when range=I ")
  end if
  Dctl%ilu = ilu

  ltest = ( ( ilu(2)>=ilu(1) ) .and. ilu(1)>=1 .and. ilu(2)<=Dctl%npwtot )
  write(msg,'(a,2i0)')" Illegal value for ilu: ",ilu
  ABI_CHECK(ltest,msg)
  Dctl%nband_k= ilu(2)-ilu(1)+1

 CASE ("V")
  if (.not.PRESENT(vlu)) then
    ABI_ERROR(" vlu must be specified when range=V ")
  end if
  Dctl%vlu = vlu

  Dctl%nband_k=-1 !??

  ltest = (vlu(2)>vlu(1))
  write(msg,'(a,2f0.3)')" Illegal value for vlu: ",vlu
  ABI_CHECK(ltest,msg)

 CASE DEFAULT
   ABI_ERROR(" Unknown value for range: "//TRIM(Dctl%range))
 END SELECT

 ! Consider the case in which we asked for the entire set of eigenvectors
 ! but the number of bands is less that npw_k. Therefore have to prepare the call to ZHEEVX.
 ! TODO this has to be done in a cleaner way.
 if (Dctl%range == "A" .and. .not. dctl%do_full_diago) then
   Dctl%range="I"
   Dctl%ilu(1) = 1
   Dctl%ilu(2) = npw_k*nspinor
   Dctl%nband_k= npw_k*nspinor
 end if

 Dctl%use_scalapack=0
 if (PRESENT(use_scalapack)) Dctl%use_scalapack=use_scalapack
 ABI_CHECK(Dctl%use_scalapack==0," scalapack mode not coded yet")

 call destroy_mpi_enreg(mpi_enreg_seq)

end subroutine init_ddiago_ctl
!!***

!!****f* m_ksdiago/ugb_from_diago
!! NAME
!! ugb_from_diago
!!
!! FUNCTION
!!  This routine performs the direct diagonalization of the Kohn-Sham Hamiltonian
!!  for a given k-point and spin using Scalapack/ELPA.
!!
!! INPUTS
!!  spin= spin index.
!!  istwf_k= Storage mode for wavefunctions.
!!  kpoint(3)= k-point in reduced coordinates
!!  ecut= Cutoff energy
!!  gs_fermie=Fermi level as computed from the previous GS run.
!!  nband_k=Number of bands
!!  prtvol=Verbosity level
!!  nfftf=(effective) number of FFT grid points in the dense FFT mesh (for this processor)
!!         (nfftf=nfft for norm-conserving potential runs)
!!  pawtab(psps%ntypat*psps%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pawfgr<pawfgr_type>=fine grid parameters and related data
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  vtrial(nfftf,nspden)=the trial potential
!!  comm=MPI communicator.
!!  nfftc=Number of points in the coarse FFT mesh.
!!  ngfftc(18)=Info about 3D FFT for the coarse mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  [electronpositron] <electronpositron_type>=quantities for the electron-positron annihilation.
!!
!! OUTPUT
!!  eig_k(1:nband_k)=The calculatated eigenvalues in ascending order.
!!
!! SOURCE

subroutine ugb_from_diago(ugb, spin, istwf_k, kpoint, ecut, gs_fermie, nband_k, ngfftc, nfftf, &
                          dtset, pawtab, pawfgr, paw_ij, cryst, psps, vtrial, eig_k, hyb, comm, &
                          electronpositron) ! Optional arguments

!Arguments ------------------------------------
!scalars
 class(ugb_t),target,intent(out) :: ugb
 integer,intent(in) :: spin, istwf_k
 real(dp),intent(in) :: kpoint(3), ecut, gs_fermie
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: comm,nfftf
 integer,intent(inout) :: nband_k
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(hyb_t),intent(inout) :: hyb
!arrays
 integer,intent(in) :: ngfftc(18)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 !real(dp),intent(inout) :: vxctau(nfftf, dtset%nspden, 4*usevxctau)
 real(dp),allocatable,intent(out) :: eig_k(:)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij(cryst%natom*psps%usepaw)
 type(electronpositron_type),optional,pointer :: electronpositron

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem1 = 1, tim_getghc = 4, paral_kgb0 = 0, master = 0, ncomp1 = 1
 integer :: cprj_choice,cpopt,dimffnl,ib,ider,idir,npw_k,nfftc,mgfftc, igs, ige, omp_nt
 integer :: jj,n1,n2,n3,n4,n5,n6,nkpg,nproc,my_rank,optder, ib_glob, nb_glob, ib_loc
 integer :: type_calc,sij_opt,igsp2_start,ig, my_ib, ibs1, ipwsp, islice, igsp_loc, my_npwsp
 integer :: npwsp, col_bsize, nsppol, nspinor, nspden, loc2_size, il_g1, il_g2, ig1, ig2, ierr, min_my_nband, band_sum
 integer :: idat, ndat, batch_size, h_size !, mene_found
 integer :: ik_ibz, ik_bz, isym_k, trev_k, g0_k(3), g0(3)
 integer :: iq_bz !, isym_q, trev_q !, g0_q(3) iq_ibz
 real(dp),parameter :: lambda0 = zero
 real(dp) :: cpu, wall, gflops, mem_mb, f_bsum, fact_spin, tol_empty_in, tol_empty, gsq_max, inv_sqrt_ucvol, rcut
 logical :: do_full_diago, haveit, isirr_k, q_is_gamma ! isirr_q,
 character(len=80) :: frmt1
 character(len=10) :: stag(2)
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg_seq
 type(gs_hamiltonian_type) :: gs_hamk
 type(slkmat_dp_t) :: ghg_mat, gsg_mat, ghg_4diag, gsg_4diag, eigvec
 type(slk_processor_t) :: proc_1d, proc_4diag
 type(uplan_t) :: uplan_k
 type(fftbox_plan3_t) :: box_plan
 type(psbands_t) :: psb
!arrays
 integer,allocatable :: gfft(:,:)
 real(dp) :: kptns_(3,1), ylmgr_dum(1,1,1), tsec(2), ksum(3), kk_ibz(3), kgw_m_ksum(3), qq_bz(3), my_gw_qlwl(3) ! q0(3),
 real(dp),allocatable :: ph3d(:,:,:), ffnl(:,:,:,:), kinpw(:), kpg_k(:,:), thetas(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:), ylm_k(:,:), dum_ylm_gr_k(:,:,:), eig_ene(:), ghc(:,:), gvnlxc(:,:), gsc(:,:), vcg_qbz(:,:)
 real(dp),target,allocatable :: bras(:,:)
 complex(dp),allocatable :: ps_ug(:,:,:)
 complex(gwpc),allocatable :: cbras_box(:,:), cbras_g(:,:), vc_sqrt(:), ur(:), rfg_box(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
! *********************************************************************

 call timab(1919, 1, tsec)
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! See sequence of calls in vtorho.
 ! Check that usekden is not 0 if want to use vxctau
 !with_vxctau = dtset%usekden/=0

 if (dtset%usekden/=0) then
   ABI_ERROR("nscf_init with mgga not yet coded")
 end if
 ! Check if want to use vxctau
 !with_vxctau = (present(vxctau).and.usevxctau/=0)

 !====================
 !=== Check input ====
 !====================
 if (all(istwf_k /= [1, 2])) then
   ABI_ERROR(sjoin("istwfk:", itoa(istwf_k), "not allowed:"))
 end if

 if (istwf_k == 2) then
   ABI_ERROR("istwfk == 2 with direct diago is still under development")
   !ABI_WARNING("istwfk == 2 with direct diago is still under development")
 end if

 if (dtset%ixc < 0) then
   if (libxc_functionals_ismgga() .and. .not. libxc_functionals_is_potential_only()) then
     ABI_ERROR("meta-gga functionals are not compatible with direct diagonalization!")
   end if
 end if

 ! MPI_type for sequential part.
 call initmpi_seq(mpi_enreg_seq)
 call init_distribfft_seq(mpi_enreg_seq%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 if (pawfgr%usefinegrid /= 0) then
   call init_distribfft_seq(mpi_enreg_seq%distribfft, 'f', pawfgr%ngfft(2), pawfgr%ngfft(3), 'all')
 end if

 nspinor = dtset%nspinor; nsppol = dtset%nsppol; nspden = dtset%nspden
 if (nsppol == 1) stag = ['          ','          ']
 if (nsppol == 2) stag = ['SPIN UP:  ','SPIN DOWN:']

 ! Get g-vectors from kpt and ecut.
 call get_kg(kpoint, istwf_k, ecut, cryst%gmet, npw_k, ugb%kg_k)
 npwsp = npw_k * nspinor

 ! The coarse FFT mesh for the application of the Hamiltonian.
 n1 = ngfftc(1); n2 = ngfftc(2); n3 = ngfftc(3)
 n4 = ngfftc(4); n5 = ngfftc(5); n6 = ngfftc(6)
 nfftc = product(ngfftc(1:3)); mgfftc = maxval(ngfftc(1:3))

 ! Initialize the Hamiltonian on the coarse FFT mesh.
 if (present(electronpositron)) then
   call gs_hamk%init(psps, pawtab, nspinor, nsppol, nspden, cryst%natom, cryst%typat, cryst%xred, nfftc, &
    mgfftc, ngfftc, cryst%rprimd, dtset%nloalg, paw_ij=paw_ij, usecprj=0, electronpositron=electronpositron)
 else
   call gs_hamk%init(psps, pawtab, nspinor, nsppol, nspden, cryst%natom, cryst%typat, cryst%xred, nfftc, &
    mgfftc, ngfftc, cryst%rprimd, dtset%nloalg, paw_ij=paw_ij, usecprj=0)
 end if

 ! Check on the number of stored bands.
 if (nband_k == -1 .or. nband_k >= npwsp) then
   nband_k = npwsp
   write(msg,'(4a, i0)')ch10,&
    ' Since the number of bands to be computed was -1 or',ch10,&
    ' too large, it has been set to the maximum value. npw_k*nspinor: ',npwsp
   call wrtout(std_out, msg)
 end if

 do_full_diago = nband_k == npwsp

 ! Set up local potential vlocal with proper dimensioning, from vtrial.
 ! Select spin component of interest if nspden<=2 as nvloc==1, for nspden==4, nvloc==4
 ! option=2: vtrial(n1*n2*n3,ispden) --> vlocal(nd1,nd2,nd3) real case

 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamk%nvloc))
 call gspot_transgrid_and_pack(spin, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
                               nspden, gs_hamk%nvloc, ncomp1, pawfgr, mpi_enreg_seq, vtrial, vlocal)
 call gs_hamk%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

 ! TODO: This for meta-gga.
 !if (with_vxctau) then
 !  call gspot_transgrid_and_pack(spin, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
 !                                nspden, gs_hamk%nvloc, 4, pawfgr, mpi_enreg, vxctau, vxctaulocal)
 !  call gs_hamk%load_spin(spin, vxctaulocal=vxctaulocal)
 !end if

 !========================
 !==== Kinetic energy ====
 !========================
 ABI_MALLOC(kinpw, (npw_k))
 call mkkin(ecut, dtset%ecutsm, dtset%effmass_free, cryst%gmet, ugb%kg_k, kinpw, kpoint, npw_k, 0, 0)

 !================================
 !==== Non-local form factors ====
 !================================
 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2*psps%useylm))

 if (psps%useylm == 1) then
   optder = 0
   ABI_MALLOC(dum_ylm_gr_k, (npw_k, 3+6*(optder/2),psps%mpsang**2))
   kptns_(:,1) = kpoint

   ! NB: Here mband is not used if paral_compil_kpt = 0
   call initylmg(cryst%gprimd, ugb%kg_k, kptns_, mkmem1, mpi_enreg_seq, psps%mpsang, npw_k, [nband_k], 1, &
     [npw_k], 1, optder, cryst%rprimd, ylm_k, dum_ylm_gr_k)

   ABI_FREE(dum_ylm_gr_k)
 end if

 ! Compute (k+G) vectors (only if useylm=1)
 nkpg = 3 * dtset%nloalg(3)
 ABI_MALLOC(kpg_k, (npw_k, nkpg))
 if (nkpg > 0) call mkkpg(ugb%kg_k, kpg_k, kpoint, nkpg, npw_k)

 ! Compute nonlocal form factors ffnl at all (k+G):
 idir=0; ider=0; dimffnl=1+ider ! Now the derivative is not needed anymore.
 ABI_MALLOC(ffnl, (npw_k, dimffnl, psps%lmnmax, psps%ntypat))

 call mkffnl(psps%dimekb, dimffnl, psps%ekb, ffnl, psps%ffspl, cryst%gmet, cryst%gprimd, ider, idir, psps%indlmn, &
             ugb%kg_k, kpg_k, kpoint, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, &
             psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum)

 ABI_FREE(ylm_k)

 ! Load k-dependent part of the Hamiltonian.
 ABI_MALLOC(ph3d, (2, npw_k, gs_hamk%matblk))
 call gs_hamk%load_k(kpt_k=kpoint, istwf_k=istwf_k, npw_k=npw_k, kinpw_k=kinpw, &
                     kg_k=ugb%kg_k, kpg_k=kpg_k, ffnl_k=ffnl, ph3d_k=ph3d, compute_ph3d=.true., compute_gbound=.true.)

 ! Prepare call to getghc.
 type_calc = 0                                ! For applying the whole Hamiltonian
 sij_opt = 0; if (psps%usepaw==1) sij_opt = 1 ! For PAW, <k+G|S|k+G"> is also needed.

 cpopt = -1
 if (psps%usepaw == 1) cpopt = 0  ! <p_lmn|in> are computed here and saved

 ! Init 1D PBLAS grid to block-distribute H along columns.
 call proc_1d%init(comm, grid_dims=[1, nproc])
 h_size = npwsp; if (istwf_k == 2) h_size = 2*npwsp - 1

 ABI_CHECK(block_dist_1d(h_size, nproc, col_bsize, msg), msg)
 call ghg_mat%init(h_size, h_size, proc_1d, istwf_k, size_blocs=[h_size, col_bsize])
 if (psps%usepaw == 1) call gsg_mat%init(h_size, h_size, proc_1d, istwf_k, size_blocs=[h_size, col_bsize])

 ! Estimate memory
 mem_mb = ghg_mat%locmem_mb()
 mem_mb = two * (psps%usepaw + 1) * mem_mb + mem_mb  ! last term for eigvec matrix
 call wrtout(std_out, sjoin(" Local memory for scalapack matrices:", ftoa(mem_mb, fmt="(f8.1)"), ' [Mb] <<< MEM'))

 ! Define batch size for the application of the Hamiltonian
 ! This is useful if OpenMP is activated thus we use multiples of omp_nt.
 omp_nt = xomp_get_num_threads(open_parallel=.True.)
 batch_size = 8 * omp_nt
 if (istwf_k == 2) batch_size = 1  ! FIXME
 !batch_size = 1
 call wrtout(std_out, sjoin(" Building H^KS with batch_size:", itoa(batch_size)))

 ABI_MALLOC(bras, (2, npwsp * batch_size))
 ! cwaveprj is ordered by atom type, see nonlop_ylm.
 ABI_MALLOC(cwaveprj, (cryst%natom, nspinor*(1+cpopt)*gs_hamk%usepaw*batch_size))
 if (cpopt == 0) call pawcprj_alloc(cwaveprj, 0, gs_hamk%dimcprj)
 ABI_MALLOC(ghc, (2, npwsp * batch_size))
 ABI_MALLOC(gvnlxc, (2, npwsp * batch_size))
 ABI_MALLOC(gsc, (2, npwsp * batch_size*(sij_opt+1)/2))

 ! Loop over the |beta,G''> component.
 call cwtime(cpu, wall, gflops, "start")
 loc2_size = ghg_mat%size_local(2)

 if (my_rank == master) call pstat_proc%print(_PSTAT_ARGS_)

 do il_g2=1, loc2_size, batch_size
   ! Operate on ndat g-vectors starting at the igsp2_start global index.
   igsp2_start = ghg_mat%loc2gcol(il_g2)
   ndat = blocked_loop(il_g2, loc2_size, batch_size)

   bras = zero
   if (istwf_k == 1) then
     do idat=0,ndat-1
       bras(1, igsp2_start + idat * npwsp + idat) = one
     end do
   else
     ! only istwf_k == 2 is coded here. NB: there's a check at the beginning of this routine.
     do idat=0,ndat-1
       if (igsp2_start + idat <= npwsp) then
         ! Cosine term
         bras(1, igsp2_start + idat*npwsp + idat) = half
         if (igsp2_start == 1) bras(1, igsp2_start + idat*npwsp + idat) = one
       else
         ! Sine term
         !ig = igsp2_start - npwsp + 1
         ig = igsp2_start - npwsp + 1 + 1  ! This should be OK
         bras(2, ig + idat*npwsp + idat) = half
       end if
     end do
   end if

   ! Get <:|H|beta,G''> and <:|S_{PAW}|beta,G''>
   call multithreaded_getghc(cpopt, bras, cwaveprj, ghc, gsc, gs_hamk, gvnlxc, lambda0, mpi_enreg_seq, ndat, &
                             dtset%prtvol, sij_opt, tim_getghc, type_calc)

   ! Now fill my local buffer of ghg/gsg.
   if (istwf_k == 1) then
     ! Complex wavefunctions.
     do idat=0,ndat-1
       igs = 1 + idat * npwsp; ige = igs + npwsp - 1
       ghg_mat%buffer_cplx(:, il_g2+idat) = cmplx(ghc(1, igs:ige), ghc(2, igs:ige), kind=dp)
     end do
     if (psps%usepaw == 1) then
       do idat=0,ndat-1
         igs = 1 + idat * npwsp; ige = igs + npwsp - 1
         gsg_mat%buffer_cplx(:, il_g2+idat) = cmplx(gsc(1,igs:ige), gsc(2,igs:ige), kind=dp)
       end do
     end if

   else
     ! Real wavefunctions.
     do idat=0,ndat-1
       igs = 1 + idat*npwsp; ige = igs + npwsp - 1
       !if (igsp2_start == 1 .or. igsp2_start == npwsp + 1 .and. idat == 0) then
       !  ghc(:, igs:ige) = tol3 !; print *, ghc(:, igs:ige)
       !end if
       ghg_mat%buffer_real(1:npwsp,  il_g2+idat) =  ghc(1, igs:ige)     ! CC or CS
       ghg_mat%buffer_real(npwsp+1:, il_g2+idat) = -ghc(2, igs+1:ige)   ! SC or SS. Note igs+1
     end do
     if (psps%usepaw == 1) then
       NOT_IMPLEMENTED_ERROR()
       !gsg_mat%buffer_real(...)
       !gsg_mat%buffer_real(...)
     end if
   end if ! istwf_k
 end do ! il_g2

 ! MG: DEBUG
 !call wrtout(std_out, " WARNING: Setting H_KS to zero for debugging purposes!"); ghg_mat%buffer_cplx = czero
 call cwtime_report(" build H^KS_g1g2", cpu, wall, gflops)

 ! Free workspace memory allocated so far.
 ABI_FREE(bras)
 ABI_FREE(kinpw)
 ABI_FREE(vlocal)
 ABI_FREE(ghc)
 ABI_FREE(gvnlxc)
 ABI_FREE(gsc)
 if (psps%usepaw == 1 .and. cpopt == 0) call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 ! ==================================
 ! Compute Fock operator F^k_{g1,g2}
 ! ==================================
 if (dtset%usefock == 1) then
 !if (.False.) then
   call cwtime(cpu, wall, gflops, "start")
   call wrtout(std_out, sjoin(" Building Fock operator F^k_{g1,g2}  with batch_size:", itoa(batch_size)))
   ABI_CHECK(dtset%usepaw == 0, "DIRECT DIAGO OF FOCK OPERATOR WITH PAW IS NOT CODED!")
   inv_sqrt_ucvol = one / sqrt(cryst%ucvol)

   call hyb%wfd%change_ngfft(cryst, psps, ngfftc)

   ABI_MALLOC(ur, (nfftc*nspinor))
   ABI_MALLOC(cbras_g, (npw_k*nspinor, batch_size))
   ABI_MALLOC(cbras_box, (nfftc*nspinor, batch_size))
   ABI_MALLOC(rfg_box, (nfftc*nspinor, batch_size))
   ABI_MALLOC(vc_sqrt, (nfftc))
   ABI_MALLOC(vcg_qbz, (nfftc, hyb%nqbz))

   ! Set tolerance used to decide if a band is empty.
   tol_empty_in = 0.01_dp
   call get_fact_spin_tol_empty(nsppol, nspinor, tol_empty_in, fact_spin, tol_empty)

   ! Precompute the Coulomb term here to avoid tons of calls inside the loop over ig2.
   ! Get g-vectors in the FFT box for vcoul.
   ABI_MALLOC(gfft, (3, nfftc))
   call get_gfft(ngfftc, kpoint, cryst%gmet, gsq_max, gfft)

   my_gw_qlwl(:) = GW_Q0_DEFAULT; if (dtset%gw_nqlwl > 0) my_gw_qlwl = dtset%gw_qlwl(:,1)
   !my_gw_qlwl = zero
   do ik_bz=1,hyb%nkbz
     ksum = hyb%kbz(:, ik_bz)
     kgw_m_ksum = kpoint - ksum
     !print *, "kpoint", kpoint, "ksum:", ksum
     call findqg0(iq_bz, g0, kgw_m_ksum, hyb%nqbz, hyb%qbz, hyb%mG0)
     ABI_CHECK(all(g0 == 0), sjoin("g0 = ", ltoa(g0)))
     qq_bz = hyb%qbz(:,iq_bz)
     q_is_gamma = normv(qq_bz, cryst%gmet, "G") < GW_TOLQ0
     call hyb%vcgen%get_vc_sqrt(qq_bz, nfftc, gfft, my_gw_qlwl, cryst, vc_sqrt, comm, vc=vcg_qbz(:,iq_bz))
     ! A non-positive value of rcut activates the recipe of Spencer & Alavi, PRB 77, 193110 (2008) [[cite:Spencer2008]].
     rcut = (cryst%ucvol * hyb%nkbz * 3.d0 / four_pi) ** third
     !vcgen%i_sz = two_pi * rcut**2
     if (q_is_gamma) then
       !vcg_qbz(1,iq_bz) = two_pi * rcut**2      ! FIXME: This is used in GW
       !vcg_qbz(1,iq_bz) = hyb%vcgen%i_sz
       vcg_qbz(1,iq_bz) = two_pi/three * rcut**2 ! FIXME: This is used in m_fock
       !vcg_qbz(1,iq_bz) = zero
     end if
     !vcg_qbz(2:,iq_bz) = vcg_qbz(2:,iq_bz) * (inv_sqrt_ucvol**2)
     vcg_qbz(:,iq_bz) = vcg_qbz(:,iq_bz) * (inv_sqrt_ucvol**2)
     vcg_qbz(:,iq_bz) = one
     call zerosym(vcg_qbz(:,iq_bz), 1, n1, n2, n3)
   end do ! ik_bz
   ABI_FREE(gfft)

   ! Build plans for (dense, g-sphere) FFTs.
   call box_plan%from_ngfft(ngfftc, nspinor*batch_size, dtset%gpu_option)
   call uplan_k%init(npw_k, nspinor, batch_size, ngfftc, istwf_k, ugb%kg_k, gwpc, dtset%gpu_option)

   ! Blocked loop over the columns of F^k_{g1,g2}.
   do ig2=1, npwsp, batch_size
     ndat = blocked_loop(ig2, npwsp, batch_size)
     ! Fill cbras_box(r) with e^{ig2.r}.
     do idat=1,ndat
       call calc_ceigr(ugb%kg_k(:,ig2+idat-1), nfftc, nspinor, ngfftc, cbras_box(:,idat))
       cbras_box(:,idat) = cbras_box(:,idat) * inv_sqrt_ucvol
     end do

     ! ==============================
     ! ==== Sum over k in the BZ ====
     ! ==============================
     rfg_box = zero
     do ik_bz=1,hyb%nkbz
       ksum = hyb%kbz(:, ik_bz)
       ! Parallelism over k-points.
       !if (.not. hyb%wfd%ihave_ug(0, ik_ibz, spin)) cycle

       ! Find the symmetrical image of ksum in the IBZ
       ! FIXME: Be careful with the symmetry conventions here and the interplay between umklapp in q and FFT
       ik_ibz = hyb%kbz2ibz_symrel(1, ik_bz); isym_k = hyb%kbz2ibz_symrel(2, ik_bz)
       trev_k = hyb%kbz2ibz_symrel(6, ik_bz); g0_k = hyb%kbz2ibz_symrel(3:5, ik_bz)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       kk_ibz = hyb%kibz(:, ik_ibz)

       ! Identify q and G0 where q + G0 = k_GW - ksum
       !kgw_m_ksum = kpoint - ksum
       !call findqg0(iq_bz, g0, kgw_m_ksum, hyb%nqbz, hyb%qbz, hyb%mG0)
       !ABI_CHECK(all(g0 == 0), sjoin("g0 = ", ltoa(g0)))

       !qq_bz = hyb%qbz(:, iq_bz)
       !iq_ibz = hyb%qbz2ibz(1, iq_bz); isym_q = hyb%qbz2ibz(2, iq_bz)
       !trev_q = hyb%qbz2ibz(6, iq_bz); g0_q = hyb%qbz2ibz(3:5, iq_bz)
       !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))

       !! Find the corresponding irreducible q-point.
       !! NB: non-zero umklapp G_o is not allowed. There's a check in setup_sigma
       !!call qmesh%get_BZ_item(iq_bz, qbz, iq_ibz, isym_q, itim_q)
       !q_is_gamma = normv(qq_bz, cryst%gmet, "G") < GW_TOLQ0

       ! ==========================
       ! Sum over (occupied) bands
       ! ==========================
       do band_sum=1, hyb%wfd%nband(ik_ibz, spin)
         ! MPI parallelism over bands.
         if (.not. hyb%wfd%ihave_ug(band_sum, ik_ibz, spin)) cycle
         f_bsum = hyb%ebands%occ(band_sum, ik_ibz, spin) * fact_spin; if (abs(f_bsum) <= tol_empty) cycle

         !print *, "band_sum, ik_ibz, spin", band_sum, ik_ibz, spin
         call hyb%wfd%get_ur(band_sum, ik_ibz, spin, ur)
         ur = ur * inv_sqrt_ucvol
         do idat=1,ndat
           cbras_box(:,idat) = f_bsum * conjg(ur) * cbras_box(:,idat)
         end do

         ! FFT r --> g and multiply by v(g,q) on the FFT box.
         call box_plan%execute(cbras_box(:,1), -1, ndat=ndat)
         do idat=1,ndat
           cbras_box(:,idat) = cbras_box(:,idat) * vcg_qbz(:, iq_bz)
         end do
         ! FFT g --> r, multiply by u(r) and accumulate in rfg_box
         call box_plan%execute(cbras_box(:,1), +1, ndat=ndat)
         do idat=1,ndat
           rfg_box(:,idat) = rfg_box(:,idat) + cbras_box(:,idat) * ur
         end do
       end do ! band_sum
     end do ! ik_bz

     ! FFT r --> g_sphere and MPI sum partial contributions.
     call uplan_k%execute_rg(ndat, rfg_box(:,1), cbras_g(:,1))
     call xmpi_sum(cbras_g, hyb%wfd%comm_spin(spin), ierr)
     cbras_g = - cbras_g * sqrt(cryst%ucvol)
     !cbras_g = - half * cbras_g * sqrt(cryst%ucvol)
     !cbras_g = - 10000 * cbras_g / (hyb%nkbz*cryst%ucvol) ! * alpha_hyb
     !cbras_g = - sqrt(cryst%ucvol) * cbras_g ! / (hyb%nkbz*cryst%ucvol) ! * alpha_hyb
     !cbras_g = -half * inv_sqrt_ucvol * cbras_g / (hyb%nkbz*cryst%ucvol) ! * alpha_hyb
     !cbras_g = - cbras_g / (hyb%nkbz * cryst%ucvol) ! * alpha_hyb
     !cbras_g = -half * sqrt(cryst%ucvol) * cbras_g  !/ (hyb%nkbz*cryst%ucvol) ! * alpha_hyb

     ! Update my local buffer of ghg_mat.
     do idat=1,ndat
       do ig1=1,npwsp
         ! From global to local indices.
         call ghg_mat%glob2loc(ig1, ig2+idat-1, il_g1, il_g2, haveit); if (.not. haveit) cycle
         !print *, "ig1, idat, cbras_g", ig1, idat, cbras_g(ig1, idat)
         if (istwf_k == 1) then
           ! Complex wavefunctions.
           ghg_mat%buffer_cplx(il_g1, il_g2) = ghg_mat%buffer_cplx(il_g1, il_g2) + cbras_g(ig1,idat)
         else
          ! Real wavefunctions.
          NOT_IMPLEMENTED_ERROR()
         end if ! istwf_k
       end do ! ig1
     end do ! idat

   end do ! ig2

   ABI_FREE(vc_sqrt)
   ABI_FREE(vcg_qbz)
   ABI_FREE(cbras_box)
   ABI_FREE(rfg_box)
   ABI_FREE(cbras_g)
   ABI_FREE(ur)
   call uplan_k%free(); call box_plan%free()
   call cwtime_report(" build Fock_g1g2", cpu, wall, gflops)
 end if ! usefock

 !===========================================
 !=== Diagonalization of <G|H|G''> matrix ===
 !===========================================
 ABI_MALLOC(eig_ene, (h_size))

 ! Change size block. Use 2D rectangular grid of processors for diagonalization, if possible.
 call proc_4diag%init(comm)
 call ghg_mat%change_size_blocs(ghg_4diag, processor=proc_4diag, free=.True.)
 if (psps%usepaw == 1) call gsg_mat%change_size_blocs(gsg_4diag, processor=proc_4diag, free=.True.)
 !call ghg_mat%copy(ghg_4diag); call ghg_mat%free()

 ! NB: global H shape is (h_size, h_size) even for partial diago.
 ! then one extracts the (hsize, nband_k) sub-matrix before returning.
 call ghg_4diag%copy(eigvec)
 if (my_rank == master) call pstat_proc%print(_PSTAT_ARGS_)

#ifndef HAVE_LINALG_ELPA
 call wrtout([std_out, ab_out], &
 "- WARNING: Using ScaLAPACK for diagonalization, but ELPA library is highly recommended for both efficiency and memory reasons.")
#endif

 if (do_full_diago) then
   write(msg,'(5a, (a,i0), 2a)')ch10,&
     ' Begin full diagonalization for kpt: ',trim(ktoa(kpoint)), stag(spin), ch10,&
     " H_gg' Matrix size: ",npwsp, ", Scalapack grid: ", trim(ltoa(ghg_4diag%processor%grid%dims))
   call wrtout(std_out, msg)
   call cwtime(cpu, wall, gflops, "start")
   if (psps%usepaw == 0) then
     !call ghg_4diag%pzheev("V", "U", eigvec, eig_ene)
     call compute_eigen_problem(ghg_4diag%processor, ghg_4diag, eigvec, eig_ene, comm, istwf_k)
   else
     call compute_generalized_eigen_problem(ghg_4diag%processor, ghg_4diag, gsg_4diag, eigvec, eig_ene, comm, istwf_k)
   end if
   call cwtime_report(" full_diago", cpu, wall, gflops)

 else
   write(msg,'(6a,i0,(a,i0), 2a)') ch10,&
     ' Begin partial diagonalization for kpt: ',trim(ktoa(kpoint)), stag(spin), ch10,&
     " H_gg' Matrix size: ",npwsp,', nband_k: ', nband_k,", Scalapack grid: ", trim(ltoa(ghg_4diag%processor%grid%dims))
   call wrtout(std_out, msg)

   call cwtime(cpu, wall, gflops, "start")
   if (psps%usepaw == 0) then
     !call ghg_4diag%pzheevx("V", "I", "U", zero, zero, 1, nband_k, -tol8, eigvec, mene_found, eig_ene)
     call compute_eigen_problem(ghg_4diag%processor, ghg_4diag, eigvec, eig_ene, comm, istwf_k, nev=nband_k)
   else
     !call ghg_4diag%pzhegvx(1, "V", "I", "U", gsg_4diag, zero, zero, 1, nband_k, -tol8, eigvec, mene_found, eig_ene)
     call compute_generalized_eigen_problem(ghg_4diag%processor, ghg_4diag, gsg_4diag, eigvec, eig_ene, comm, istwf_k, &
                                            nev=nband_k)
   end if
   call cwtime_report(" partial_diago", cpu, wall, gflops)
 end if

 if (my_rank == master) then
   call pstat_proc%print(_PSTAT_ARGS_)
   ! Write eigenvalues.
   frmt1 = '(8x,*(1x,f7.3))'
   write(msg, '(2a,3x,a)')' Eigenvalues in eV for kpt: ', trim(ktoa(kpoint)), stag(spin); call wrtout(std_out, msg)
   write(msg, frmt1)(eig_ene(ib)*Ha_eV,ib=1,min(9,nband_k)); call wrtout(std_out, msg)
   ! HYB DEBUG
   !call wrtout(std_out, "hyb%ebands")
   !write(msg, frmt1)(hyb%ebands%eig(ib,1,spin)*Ha_eV, ib=1,min(9,hyb%ebands%mband)); call wrtout(std_out, msg)
   if (nband_k > 9 .and. dtset%prtvol > 0) then
     do jj=10,nband_k,9
       write(msg, frmt1) (eig_ene(ib)*Ha_eV,ib=jj,min(jj+8,nband_k)); call wrtout(std_out, msg)
     end do
   end if
 end if

 ! Free memory
 call ghg_4diag%free(); call gsg_4diag%free(); call proc_1d%free()

 ! ================
 ! Stochastic bands
 ! ================
 if (dtset%nb_protected /= 0) then
   call wrtout(std_out, " Generating stochastic bands...")
   ! Initial setup.
   call psb%init(dtset, h_size, eig_ene, gs_fermie) !, nband_k)
   my_npwsp = eigvec%size_local(1)
   nb_glob = eigvec%size_global(2)
   ABI_CALLOC(ps_ug, (my_npwsp, psb%maxsto_per_slice, psb%nslices))
   ABI_MALLOC(thetas, (nb_glob, psb%maxsto_per_slice))

   ! Loop over global bands.
   do ib_glob=1, nb_glob
    ! Need the same random phases on all MPI procs.
    if (eigvec%processor%my_rank == master) call random_number(thetas)
    call xmpi_bcast(thetas, master, eigvec%processor%comm, ierr)

     ! Get slice index from ib_glob.
     islice = psb%band2slice(ib_glob); if (islice == -1) cycle
     !band_block = psb%subspace(1:2, islice)
     !nb_in_slice  = psb%subspace(3,islice)

     ! Loop over global PW index.
     do ipwsp=1,npwsp
       call eigvec%glob2loc(ipwsp, ib_glob, igsp_loc, ib_loc, haveit); if (.not. haveit) cycle
       do ib=1,psb%subspace(3,islice)
         ps_ug(igsp_loc, ib, islice) = ps_ug(igsp_loc, ib, islice) + &
           eigvec%buffer_cplx(igsp_loc, ib_loc) * exp(j_dpc*two_pi*thetas(ib_glob,ib))
       end do
     end do
   end do ! ib_glob
   ABI_FREE(thetas)

   ! Normalize
   ! TODO: Need MPI communicator over columns here.
   !call xmpi_sum(ps_ug, eigvec%column_comm, ierr)
   do islice=1,psb%nslices
     do ib=1,psb%subspace(3,islice)
       if (psb%subspace(3,islice) == 1) cycle
       ps_ug(:,ib,islice) = ps_ug(:,ib,islice) / sqrt(one * psb%subspace(3,islice))
     end do
   end do

   ! Now insert ps_ug in the right position in eigevec
   do ib_glob=1, nb_glob
     islice = psb%band2slice(ib_glob); if (islice == -1) cycle
     !band_start = 1 + (islice - 1) * psb%nb_per_slice
     ! Loop over global PW index.
     do ipwsp=1,npwsp
       call eigvec%glob2loc(ipwsp, ib_glob, igsp_loc, ib_loc, haveit); if (.not. haveit) cycle
       !do ib=1,psb%nb_per_slice
       !  eigvec%buffer_cplx(igsp_loc, ib_loc) = ps_ug(igsp_loc, ib, islice)
       !end do
     end do
   end do

   ABI_FREE(ps_ug)

   ! here we change the value of nband_k and eig_k.
   nband_k = psb%nb_tot
   ABI_MALLOC(eig_k, (nband_k))
   eig_k = psb%ps_eig

 else
   ! No pseudo bands.
   ABI_MALLOC(eig_k, (nband_k))
   eig_k(:) = eig_ene(1:nband_k)
 end if

 ! Now transfer eigvec to the ugb datastructure using 1d grid (block column distribution).
 call wrtout(std_out, " Moving to PBLAS block column distribution...")
 call cwtime(cpu, wall, gflops, "start")

 call ugb%processor%init(comm, grid_dims=[1, nproc])
 ABI_CHECK(block_dist_1d(nband_k, nproc, col_bsize, msg), msg)
 call eigvec%cut(h_size, nband_k, ugb%mat, size_blocs=[h_size, col_bsize], processor=ugb%processor, free=.True.)
 call proc_4diag%free()

 ! =================
 ! Build ugb object
 ! =================
 ugb%istwf_k = istwf_k
 ugb%nspinor = nspinor
 ugb%npw_k = npw_k
 ugb%npwsp = npwsp
 ugb%nband_k = nband_k
 ugb%comm => ugb%mat%processor%comm

 ugb%my_bstart = ugb%mat%loc2gcol(1)
 ugb%my_bstop = ugb%mat%loc2gcol(ugb%mat%size_local(2))
 ugb%my_nband = ugb%my_bstop - ugb%my_bstart + 1

 if (ugb%my_nband > 0) then
   call c_f_pointer(c_loc(ugb%mat%buffer_cplx), ugb%cg_k, shape=[2, npwsp, ugb%my_nband])
 else
   ugb%my_nband = 0; ugb%cg_k => null()
 end if

 call xmpi_min(ugb%my_nband, min_my_nband, comm, ierr)
 ugb%has_idle_procs = min_my_nband == 0

 if (psps%usepaw == 1 .and. ugb%my_nband > 0) then
   ! Calculate <Proj_i|Cnk> from output eigenstates. Note array allocated with ugb%my_nband
   ABI_MALLOC(ugb%cprj_k, (cryst%natom, nspinor * ugb%my_nband))
   call pawcprj_alloc(ugb%cprj_k, 0, gs_hamk%dimcprj)
   idir = 0; cprj_choice = 1  ! Only projected wave functions.

   do my_ib=1,ugb%my_nband
     ibs1 = nspinor * (my_ib - 1) + 1
     call getcprj(cprj_choice, 0, ugb%cg_k(:,:,my_ib), ugb%cprj_k(:,ibs1), &
                  gs_hamk%ffnl_k, idir, gs_hamk%indlmn, gs_hamk%istwf_k, gs_hamk%kg_k, &
                  gs_hamk%kpg_k, gs_hamk%kpt_k, gs_hamk%lmnmax, gs_hamk%mgfft, mpi_enreg_seq, 1, &
                  gs_hamk%natom, gs_hamk%nattyp, gs_hamk%ngfft, gs_hamk%nloalg, gs_hamk%npw_k, gs_hamk%nspinor, &
                  gs_hamk%ntypat, gs_hamk%phkxred, gs_hamk%ph1d, gs_hamk%ph3d_k, gs_hamk%ucvol, gs_hamk%useylm)
   end do

   !  Reorder the cprj (order is now the same as in the input file)
   call pawcprj_reorder(ugb%cprj_k, gs_hamk%atindx1)
 end if ! usepaw

 call cwtime_report(" block column distribution completed", cpu, wall, gflops)

 ! Free memory.

 ABI_FREE(eig_ene)
 ABI_FREE(kpg_k)
 ABI_FREE(ph3d)
 ABI_FREE(ffnl)
 call destroy_mpi_enreg(mpi_enreg_seq); call gs_hamk%free(); call psb%free()

 if (my_rank == master) call pstat_proc%print(_PSTAT_ARGS_)

 call timab(1919, 2, tsec)

end subroutine ugb_from_diago
!!***

!!****f* m_ksdiago/ugb_from_wfk_file
!! NAME
!! ugb_from_wfk_file
!!
!! FUNCTION
!!  Initialize an ugb_t instance from a WFK file.
!!
!! INPUTS
!!  spin: spin index.
!!  kpoint(3)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  eig_k(1:nband_k)=The calculatated eigenvalues in ascending order.
!!
!! SOURCE

subroutine ugb_from_wfk_file(ugb, ik_ibz, spin, istwf_k, kpoint, nband_k, &
                             dtset, dtfil, cryst, eig_k, comm)

 use m_wfk

!Arguments ------------------------------------
!scalars
 class(ugb_t),target,intent(out) :: ugb
 integer,intent(in) :: ik_ibz, spin, istwf_k
 real(dp),intent(in) :: kpoint(3)
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 integer,intent(in) :: nband_k, comm
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),allocatable,intent(out) :: eig_k(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, formeig0 = 0
 integer :: ierr, bcast_comm, color, min_my_nband
 integer :: nprocs, my_rank, nbsum, npwsp, bstart, bstop, band_step, nb, npw_k, col_bsize, band, ib, il_b, iloc
 logical :: have_band
 type(ebands_t) :: wfk_ebands
 type(wfk_t) :: wfk
 type(hdr_type) :: wfk_hdr
 character(len=fnlen) :: wfk_path
 character(len=500) :: msg
 !type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer :: units(2)
 real(dp),target,allocatable :: cg_work(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: cg_k(:,:)
! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 units(:) = [std_out, ab_out]

 wfk_path = dtfil%fnamewffk
 if (my_rank == master) then
   if (nctk_try_fort_or_ncfile(wfk_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find HYBRYD WFK file:", wfk_path, ". Error:", msg))
   end if
   call wrtout(units, sjoin("- Reading HYBRID orbitals from WFK file:", wfk_path), pre_newlines=2)
 end if

 ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
 call xmpi_bcast(wfk_path, master, comm, ierr)

 ! Read energies and performs some basic consistency checks.
 wfk_ebands = wfk_read_ebands(wfk_path, comm, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ABI_CHECK_IEQ(dtset%ixc, wfk_hdr%ixc, "dtset%ixc /= wfk_hdr%ixc")
 ABI_CHECK(all(abs(wfk_hdr%kptns(:,ik_ibz) - kpoint) < tol6), "Different kpoint")
 ABI_CHECK_IRANGE(nband_k, 1, wfk_ebands%mband, "nband_k > mband.")

 ABI_MALLOC(eig_k, (nband_k))
 eig_k = wfk_ebands%eig(1:nband_k, ik_ibz, spin)

 npw_k = wfk_hdr%npwarr(ik_ibz)
 npwsp = npw_k * wfk_hdr%nspinor
 ABI_CHECK_IEQ(istwf_k, wfk_hdr%istwfk(ik_ibz), "different istwfk_k")

 ! Init scalapack matrix
 call ugb%processor%init(comm, grid_dims=[1, nprocs])
 ABI_CHECK(block_dist_1d(nband_k, nprocs, col_bsize, msg), msg)
 call ugb%mat%init(npwsp, nband_k, ugb%processor, 1, size_blocs=[-1, col_bsize])

 ABI_MALLOC(ugb%kg_k, (3, npw_k))

 ! Master reads and broadcasts. Much faster on lumi
 if (my_rank == master) then
   call wfk%open_read(wfk_path, formeig0, iomode_from_fname(wfk_path), get_unit(), xmpi_comm_self)
 end if

 ! TODO: Optimize this part
 ! Find band_step that gives good compromise between memory and efficiency.
 !band_step = memb_limited_step(1, nbsum, 2*npwsp, xmpi_bsize_dp, 1024.0_dp)
 band_step = 200
 !band_step = 100
 nbsum = nband_k
 do bstart=1, nbsum, band_step
   bstop = min(bstart + band_step - 1, nbsum); nb = bstop - bstart + 1

   ABI_MALLOC(cg_work, (2, npwsp, nb)) ! This array is always dp
   if (my_rank == master) then
     call c_f_pointer(c_loc(cg_work), cg_k, shape=[2, npwsp * nb])
     call wfk%read_band_block([bstart, bstop], ik_ibz, spin, xmpio_single, kg_k=ugb%kg_k, cg_k=cg_k)
   end if

   call xmpi_bcast(ugb%kg_k, master, comm, ierr)

   ! Create communicator with master and all procs requiring this set of bands block (color == 1)
   color = 0
   do band=bstart, bstop
     call ugb%mat%glob2loc(1, band, iloc, il_b, have_band)
     if (have_band) then
       color = 1; exit
     end if
   end do
   if (my_rank == master) color = 1
   call xmpi_comm_split(comm, color, my_rank, bcast_comm, ierr)

   if (color == 1) then
     call xmpi_bcast(cg_work, master, bcast_comm, ierr)
   endif
   call xmpi_comm_free(bcast_comm)

   ! Copy my portion of cg_work to buffer_cplx (here we have dp --> sp conversion).
   if (color == 1) then
     do band=bstart, bstop
       ib = band - bstart + 1
       call ugb%mat%glob2loc(1, band, iloc, il_b, have_band); if (.not. have_band) cycle
       ugb%mat%buffer_cplx(:, il_b) = cmplx(cg_work(1,:,ib), cg_work(2,:,ib), kind=gwpc)
     end do
   end if
   ABI_FREE(cg_work)
 end do ! bstart

 if (my_rank == master) call wfk%close()

 ! =================
 ! Build ugb object
 ! =================
 ugb%istwf_k = istwf_k
 ugb%nspinor = wfk_hdr%nspinor
 ugb%npw_k = npw_k
 ugb%npwsp = npwsp
 ugb%nband_k = nband_k
 ugb%comm => ugb%mat%processor%comm

 ugb%my_bstart = ugb%mat%loc2gcol(1)
 ugb%my_bstop = ugb%mat%loc2gcol(ugb%mat%size_local(2))
 ugb%my_nband = ugb%my_bstop - ugb%my_bstart + 1

 if (ugb%my_nband > 0) then
   call c_f_pointer(c_loc(ugb%mat%buffer_cplx), ugb%cg_k, shape=[2, ugb%npwsp, ugb%my_nband])
 else
   ugb%my_nband = 0
   ugb%cg_k => null()
 end if

 call xmpi_min(ugb%my_nband, min_my_nband, comm, ierr)
 ugb%has_idle_procs = min_my_nband == 0

 ! TODO
 if (dtset%usepaw == 1 .and. ugb%my_nband > 0) then
    ABI_ERROR("ugb_from_wfk does not support PAW")
    ! Calculate <Proj_i|Cnk> from output eigenstates. Note array allocated with ugb%my_nband
    ABI_MALLOC(ugb%cprj_k, (cryst%natom, ugb%nspinor * ugb%my_nband))
    !call pawcprj_alloc(ugb%cprj_k, 0, gs_hamk%dimcprj)
    !idir = 0; cprj_choice = 1  ! Only projected wave functions.

    !do my_ib=1,ugb%my_nband
    !  ibs1 = nspinor * (my_ib - 1) + 1
    !  call getcprj(cprj_choice, 0, ugb%cg_k(:,:,my_ib), ugb%cprj_k(:,ibs1), &
    !               gs_hamk%ffnl_k, idir, gs_hamk%indlmn, gs_hamk%istwf_k, gs_hamk%kg_k, &
    !               gs_hamk%kpg_k, gs_hamk%kpt_k, gs_hamk%lmnmax, gs_hamk%mgfft, mpi_enreg_seq, &
    !               gs_hamk%natom, gs_hamk%nattyp, gs_hamk%ngfft, gs_hamk%nloalg, gs_hamk%npw_k, gs_hamk%nspinor, &
    !               gs_hamk%ntypat, gs_hamk%phkxred, gs_hamk%ph1d, gs_hamk%ph3d_k, gs_hamk%ucvol, gs_hamk%useylm)
    !end do

    !!  Reorder the cprj (order is now the same as in the input file)
    !call pawcprj_reorder(ugb%cprj_k, gs_hamk%atindx1)
 end if ! usepaw

 call ugb%print(units, dtset%prtvol)

 call wfk_hdr%free(); call wfk_ebands%free()

end subroutine ugb_from_wfk_file
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/ugb_free
!! NAME
!!  ugb_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! SOURCE

subroutine ugb_free(ugb)

!Arguments ------------------------------------
 class(ugb_t),intent(inout) :: ugb
! *************************************************************************

 call ugb%mat%free()
 call ugb%processor%free()
 ABI_SFREE(ugb%kg_k)
 ugb%cg_k => null()
 ugb%comm => null()

 if (allocated(ugb%cprj_k)) then
   call pawcprj_free(ugb%cprj_k)
   ABI_FREE(ugb%cprj_k)
 end if

end subroutine ugb_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/ugb_print
!! NAME
!!  ugb_print
!!
!! FUNCTION
!!  Print info on the object.
!!
!! SOURCE

subroutine ugb_print(ugb, units, prtvol, header)

!Arguments ------------------------------------
 class(ugb_t),intent(in) :: ugb
 integer,intent(in) :: units(:), prtvol
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 character(len=500) :: msg
 type(yamldoc_t) :: ydoc
! *************************************************************************

 ABI_UNUSED(prtvol)

 msg = ' ==== Info on the ugb_t object ==== '; if (present(header)) msg = ' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(units, msg)

 ydoc = yamldoc_open('ugb_t') !, width=11, real_fmt='(3f8.3)')
 call ydoc%add_int("istwf_k", ugb%istwf_k)
 call ydoc%add_int("nspinor", ugb%nspinor)
 call ydoc%add_int("npw_k", ugb%npw_k)
 call ydoc%add_int("nband_k", ugb%nband_k)
 call ydoc%add_int("my_bstart", ugb%my_bstart)
 call ydoc%add_int("my_bstop", ugb%my_bstop)
 call ydoc%add_int("my_nband", ugb%my_nband)
 call ydoc%write_units_and_free(units)

end subroutine ugb_print
!!***
!----------------------------------------------------------------------

!!****f* m_ksdiago/ugb_collect_cprj
!! NAME
!!  ugb_collect_cprj
!!
!! FUNCTION
!!  This is a collective routine that returns in `out_cprj` the PAW projections
!!  for `nb` bands starting at `band_start` NB: `out_cprj` is supposed to be allocated in the parent
!!
!! SOURCE

subroutine ugb_collect_cprj(ugb, nspinor, nb, band_start, out_cprj)

!Arguments ------------------------------------
 class(ugb_t),intent(in) :: ugb
 integer,intent(in) :: nspinor, nb, band_start
 type(pawcprj_type),intent(inout) :: out_cprj(:,:)

!Local variables-------------------------------
 integer :: ierr, my_ibs, out_ibs, band, cnt
! *************************************************************************

 ABI_CHECK_IEQ(size(ugb%cprj_k, dim=1), size(out_cprj, dim=1), "size1 should be the same")
 ABI_CHECK_IGEQ(size(out_cprj, dim=2), nb*nspinor, "size2 too small!")

 ! TODO: Numb algorithm based on xmpi_sum. Might be optimized.
 call pawcprj_set_zero(out_cprj)

 cnt = nspinor - 1
 do band=band_start, band_start+nb-1
   if (band >= ugb%my_bstart .and. band <= ugb%my_bstop) then
     my_ibs = 1 + (band - ugb%my_bstart) * nspinor
     out_ibs = 1 + (band - band_start) * nspinor
     call pawcprj_copy(ugb%cprj_k(:,my_ibs:my_ibs+cnt), out_cprj(:,out_ibs:out_ibs+cnt))
   end if
 end do

 call pawcprj_mpi_sum(out_cprj, ugb%comm, ierr)

end subroutine ugb_collect_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/hyb_from_wfk_file
!! NAME
!!  hyb_from_wfk_file
!!
!! FUNCTION
!!  Read the WFK file compute with HYBRID functional
!!
!! SOURCE

subroutine hyb_from_wfk_file(hyb, cryst, dtfil, dtset, psps, pawtab, ngfftc, diago_pool, comm)

 use m_krank
 use m_kpts

!Arguments ------------------------------------
 class(hyb_t),intent(out) :: hyb
 type(crystal_t),intent(in) :: cryst
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(xmpi_pool2d_t),intent(in) :: diago_pool
 integer,intent(in) :: ngfftc(18), comm

!Local variables ------------------------------
 integer,parameter :: master = 0
 integer :: nprocs, my_rank, ierr, mband, nkibz, nsppol, spin, ik_ibz, ebands_kptopt ! b1, b2,
 real(dp) :: vc_ecut
 character(len=5000) :: msg
 type(hdr_type) :: wfk_hdr
 type(crystal_t) :: wfk_cryst
 type(krank_t) :: krank_ibz ! qrank,
 character(len=fnlen) :: wfk_path
 integer :: units(2)
 integer :: nqbzX
 integer,allocatable :: nband(:,:), wfd_istwfk(:), qtab(:), qtabi(:), qtabo(:)
 real(dp),allocatable :: qbz(:,:), wtk(:), wtq(:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
!************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 units(:) = [std_out, ab_out]

 wfk_path = dtfil%fnamewffk
 if (my_rank == master) then
   if (nctk_try_fort_or_ncfile(wfk_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find HYBRYD WFK file:", wfk_path, ". Error:", msg))
   end if
   call wrtout(units, sjoin("- Reading HYBRID orbitals from WFK file:", wfk_path), pre_newlines=2)
 end if

 ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
 call xmpi_bcast(wfk_path, master, comm, ierr)

 ! Construct crystal and hyb%ebands from the GS WFK file.
 hyb%ebands = wfk_read_ebands(wfk_path, comm, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ABI_CHECK_IEQ(dtset%ixc, wfk_hdr%ixc, "dtset%ixc /= wfk_hdr%ixc")

 wfk_cryst = wfk_hdr%get_crystal()
 if (cryst%compare(wfk_cryst, header=" Comparing input crystal with WFK crystal") /= 0) then
   ABI_ERROR("Crystal structure from input and from WFK file do not agree! Check messages above!")
 end if
 !call wfk_cryst%print(header="crystal structure from WFK file")
 call wfk_cryst%free()
 ! TODO: Add more consistency checks e.g. nkibz,...
 !cryst = wfk_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 nkibz = hyb%ebands%nkpt; nsppol = hyb%ebands%nsppol
 mband = hyb%ebands%mband

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k wavevectors
 ! treated by this MPI rank are stored.
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 ! Set tolerance used to decide if a band is empty
 !tol_empty_in = 0.01_dp
 !call get_fact_spin_tol_empty(nsppol, nspinor, tol_empty_in, fact_spin, tol_empty)

 !do hyb_ik_ibz=1,nkibz
 do spin=1,nsppol
   if (all(.not. diago_pool%treats(:, spin))) cycle ! MPI distribution of collinear spins.
   do ik_ibz=1,nkibz
     bks_mask(:, ik_ibz, spin) = .True.
     !bks_mask(b1:b2, ik_ibz, spin) = .True.
   end do
 end do
 !end do

 ! Impose istwfk = 1 for all k-points.
 ! wfd_read_wfk will handle a possible conversion if the WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1 !; wfd_istwfk = wfk_hdr%istwf_k

 call hyb%wfd%init(cryst, pawtab, psps, keep_ur, mband, nband, nkibz, dtset%nsppol, bks_mask, &
               dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, hyb%ebands%kptns, ngfftc, &
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call hyb%wfd%print([std_out], header="Wavefunctions for Hybrid WKF file")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 call wfk_hdr%free()

 ! Read wavefunctions.
 call hyb%wfd%read_wfk(wfk_path, iomode_from_fname(wfk_path))

 ! This piece of code is taken from m_gwr.

 ! =======================
 ! Setup k-mesh and q-mesh
 ! =======================
 ! Get full kBZ associated to hyb%ebands
 call kpts_ibz_from_kptrlatt(cryst, hyb%ebands%kptrlatt, hyb%ebands%kptopt, hyb%ebands%nshiftk, hyb%ebands%shiftk, &
                             hyb%nkibz, hyb%kibz, wtk, hyb%nkbz, hyb%kbz) !, bz2ibz=bz2ibz)
                             !new_kptrlatt=gwr%kptrlatt, new_shiftk=gwr%kshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME
 ABI_FREE(wtk)

 ! In principle kibz should be equal to hyb%ebands%kptns.
 ABI_CHECK_IEQ(hyb%nkibz, hyb%ebands%nkpt, "nkibz != hyb%ebands%nkpt")
 ABI_CHECK(all(abs(hyb%ebands%kptns - hyb%kibz) < tol12), "hyb%ebands%kibz != hyb%kibz")

 ! Note symrec convention.
 ebands_kptopt = hyb%ebands%kptopt
 krank_ibz = krank_from_kptrlatt(hyb%nkibz, hyb%kibz, hyb%ebands%kptrlatt, compute_invrank=.False.)

 ABI_MALLOC(hyb%kbz2ibz, (6, hyb%nkbz))
 if (kpts_map("symrec", ebands_kptopt, cryst, krank_ibz, hyb%nkbz, hyb%kbz, hyb%kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 ! Order kbz by stars and rearrange entries in kbz2ibz table.
 call kpts_pack_in_stars(hyb%nkbz, hyb%kbz, hyb%kbz2ibz)

 if (my_rank == master) then
   call kpts_map_print(units, " Mapping kBZ --> kIBZ", "symrec", hyb%kbz, hyb%kibz, hyb%kbz2ibz, dtset%prtvol)
 end if

 ! Table with symrel conventions for the symmetrization of the wfs.
 ABI_MALLOC(hyb%kbz2ibz_symrel, (6, hyb%nkbz))
 if (kpts_map("symrel", ebands_kptopt, cryst, krank_ibz, hyb%nkbz, hyb%kbz, hyb%kbz2ibz_symrel) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if
 call krank_ibz%free()

 ! Setup qIBZ, weights and BZ.
 ! Always use q --> -q symmetry even in systems without inversion
 ! TODO: Might add input variable to rescale the q-mesh.

 ! Find the number of q-points such that q = k1-k2.
 call findnq(hyb%nkbz, hyb%kbz, cryst%nsym, cryst%symrec, cryst%symafm, hyb%nqibz, cryst%timrev)

 ! Find the coordinates of the q-points in the IBZ.
 ABI_MALLOC(hyb%qibz, (3, hyb%nqibz))
 call findq(hyb%nkbz, hyb%kbz, cryst%nsym, cryst%symrec, cryst%symafm, cryst%gprimd, hyb%nqibz, hyb%qibz, cryst%timrev)
 ABI_CHECK(all(abs(hyb%qibz(:,1)) < tol16), "First qpoint in qibz should be Gamma!")

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 !ABI_MALLOC(hyb%qbz2ibz, (6, hyb%nqbz))
 !qrank = krank_from_kptrlatt(hyb%nqibz, hyb%qibz, qptrlatt, compute_invrank=.False.)

 !if (kpts_map("symrec", qtimrev1, cryst, qrank, hyb%nqbz, hyb%qbz, hyb%qbz2ibz) /= 0) then
 !  ABI_ERROR("Cannot map qBZ to IBZ!")
 !end if
 !call qrank%free()

 ! Order qbz by stars and rearrange entries in qbz2ibz table.
 !call kpts_pack_in_stars(hyb%nqbz, hyb%qbz, hyb%qbz2ibz)
 !if (my_rank == master) then
 !  call kpts_map_print(units, " Mapping qBZ --> qIBZ", "symrec", hyb%qbz, hyb%qibz, hyb%qbz2ibz, dtset%prtvol)
 !end if

 nqbzX = hyb%nqibz*cryst%nsym*cryst%timrev ! Maximum possible number
 ABI_MALLOC(qbz, (3, nqbzX))
 ABI_MALLOC(wtq, (hyb%nqibz))
 ABI_MALLOC(qtab, (nqbzX))
 ABI_MALLOC(qtabi, (nqbzX))
 ABI_MALLOC(qtabo, (nqbzX))

 call identk(hyb%qibz, hyb%nqibz, nqbzX, cryst%nsym, cryst%timrev, cryst%symrec, cryst%symafm, qbz, qtab, qtabi, qtabo, hyb%nqbz, wtq)

 ABI_MALLOC(hyb%qbz, (3, hyb%nqibz))
 hyb%qbz = qbz(:,1:hyb%nqibz)

 ABI_FREE(qbz)
 ABI_FREE(wtq)
 ABI_FREE(qtab)
 ABI_FREE(qtabi)
 ABI_FREE(qtabo)

 ! TODO: MC technique does not seem to work as expected, even in the legacy code.
 vc_ecut = dtset%ecut ! * four
 call hyb%vcgen%init(cryst, hyb%ebands%kptrlatt, hyb%nkbz, hyb%nqibz, hyb%nqbz, hyb%qbz, &
                     dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, vc_ecut, comm)

end subroutine hyb_from_wfk_file
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/hyb_free
!! NAME
!!  hyb_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! SOURCE

subroutine hyb_free(hyb)

!Arguments ------------------------------------
 class(hyb_t),intent(inout) :: hyb
! *************************************************************************

 ABI_SFREE(hyb%kibz)
 ABI_SFREE(hyb%kbz)
 ABI_SFREE(hyb%qibz)
 ABI_SFREE(hyb%qbz)
 ABI_SFREE(hyb%wtq)
 ABI_SFREE(hyb%kbz2ibz)
 ABI_SFREE(hyb%kbz2ibz_symrel)
 ABI_SFREE(hyb%qbz2ibz)

 ! Free datatypes
 call hyb%wfd%free(); call hyb%vcgen%free(); call hyb%ebands%free()

end subroutine hyb_free
!!***

!!****f* m_ksdiago/psbands_init
!! NAME
!! psbands_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine psbands_init(psb, dtset, eig_size, eig_k, gs_fermie)

!Arguments ------------------------------------
 class(psbands_t),intent(out) :: psb
 class(dataset_type),target,intent(in) :: dtset
 integer,intent(in) :: eig_size
 real(dp),intent(in) :: gs_fermie
!arrays
 real(dp),intent(in) :: eig_k(eig_size)

!Local variables-------------------------------
!scalars
 integer :: islice, ib, cnt, units(2), first_band, last_band, nb
 real(dp) :: first_eig, last_eig
 real(dp),allocatable :: tmp_eig_k(:)
! *********************************************************************

 ! Shift energies wrt the input Fermi level.
 ABI_MALLOC(tmp_eig_k, (eig_size))
 tmp_eig_k = eig_k - gs_fermie

 psb%nb_protected = dtset%nb_protected
 psb%maxsto_per_slice = dtset%nb_per_slice
 ! TODO
 psb%efrac = 0.02_dp   ! dtset%efrac

 ! Compute nslices and subspace
 ! TODO: Add possibility of treating occupied states as well?
 ABI_MALLOC(psb%subspace, (3, eig_size))
 first_band = psb%nb_protected + 1
 psb%nslices = 0

 do while (first_band > 0)
   first_eig = tmp_eig_k(first_band)
   last_eig = first_eig + (first_eig * psb%efrac)
   last_band = get_band_with_energy_small_than(first_band+1, eig_size, last_eig)
   psb%nslices = psb%nslices + 1
   psb%subspace(1, psb%nslices) = first_band
   if (last_band == -1) then
     psb%subspace(2, psb%nslices) = eig_size
   else
     psb%subspace(2, psb%nslices) = last_band
   end if
   nb = psb%subspace(2, psb%nslices) - psb%subspace(1, psb%nslices)  + 1
   if (last_band == first_band) then
     ! Won't use pseudo bands in this case.
     psb%subspace(3, psb%nslices) = 1
   else
     psb%subspace(3, psb%nslices) = min(dtset%nb_per_slice, nb)
   end if
   first_band = last_band + 1
   !write(std_out,'(a,i0,a,*(1x,i0))')" islice: ", psb%nslices, " subspace:", psb%subspace(:, psb%nslices)
 end do

 ! Copy eigenvalues of the protected states.
 psb%nb_tot = psb%nb_protected + sum(psb%subspace(3,1:psb%nslices))
 ABI_MALLOC(psb%ps_eig, (psb%nb_tot))
 psb%ps_eig(1:psb%nb_protected) = tmp_eig_k(1:psb%nb_protected)

 cnt = 0
 do islice=1,psb%nslices
   first_band = psb%subspace(1, islice)
   last_band = psb%subspace(2, islice)
   ! Take average of eigenvalues inside the slice.
   do ib=1,psb%subspace(3, islice)
     cnt = cnt + 1
     psb%ps_eig(psb%nb_protected + cnt) = sum(tmp_eig_k(first_band:last_band)) / dble(last_band - first_band + 1)
   end do
 end do
 ABI_FREE(tmp_eig_k)

 psb%ps_eig = psb%ps_eig + gs_fermie

 units = [std_out, ab_out]
 call wrtout(units, ' Stochastic pseudobands setup:', pre_newlines=1)
 call wrtout(units, sjoin('     Number of stochastic subspaces: ', itoa(psb%nslices)))
 call wrtout(units, sjoin('     Number of stochastic pseudobands per subspace: ', itoa(dtset%nb_per_slice)))
 call wrtout(units, sjoin('     Original number of bands: ', itoa(eig_size)))
 call wrtout(units, sjoin('     Number of bands in the protection window: ', itoa(psb%nb_protected)))
 call wrtout(units, sjoin('     Final number of bands: ', itoa(psb%nb_tot)), newlines=1)

 !if (dtset%prtvol > 5) then
 !  do islice=1,psb%nslices
 !    write(msg,'(a,i0,a,*(1x,i0))')" islice: ", psb%nslices, " subspace:", psb%subspace(:, psb%nslices)
 !    call wrtout(units, msg)
 !  end do
 !end if

contains

integer function get_band_with_energy_small_than(idx_start, idx_end, energy) result(band)
  integer, intent(in) :: idx_start, idx_end
  integer :: ib
  real(dp), intent(in) :: energy

  band = -1
  do ib=idx_start,idx_end
    if (tmp_eig_k(ib) > energy) then
      band = ib - 1; return
    end if
  end do
end function get_band_with_energy_small_than

end subroutine psbands_init
!!***

!!****f* m_ksdiago/psbands_band2slice
!! NAME
!! psbands_band2slice
!!
!! FUNCTION
!!  Return the slice index from the band index. -1 if band is protected.
!!
!! SOURCE

integer function psbands_band2slice(psb, band) result(islice)

!Arguments ------------------------------------
 class(psbands_t),intent(in) :: psb
 integer,intent(in) :: band
! *********************************************************************

 do islice=1,psb%nslices
   if (band >= psb%subspace(1,islice) .and. &
       band <= psb%subspace(2,islice)) return
 end do
 islice = -1

end function psbands_band2slice
!!***

!!****f* m_ksdiago/psbands_free
!! NAME
!! psbands_free
!!
!! FUNCTION
!! Free memory
!!
!! SOURCE

subroutine psbands_free(psb)

!Arguments ------------------------------------
 class(psbands_t),intent(inout) :: psb
! *********************************************************************

 ABI_SFREE(psb%ps_eig)
 ABI_SFREE(psb%subspace)

end subroutine psbands_free
!!***

end module m_ksdiago
!!***

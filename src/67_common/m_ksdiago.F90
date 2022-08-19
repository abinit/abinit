!!****m* ABINIT/m_ksdiago
!! NAME
!!  m_ksdiago
!!
!! FUNCTION
!!  Direct diagonalization of the KS Hamiltonian H_k(G,G')
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MG)
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

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_hamiltonian
 use m_distribfft

 use defs_datatypes,      only : pseudopotential_type
 use defs_abitypes,       only : MPI_type
 use m_dtset,             only : dataset_type
 use m_fstrings,          only : toupper, ktoa, itoa, sjoin
 use m_time,              only : cwtime, cwtime_report
 use m_geometry,          only : metric
 use m_hide_lapack,       only : xhegv_cplex, xheev_cplex, xheevx_cplex, xhegvx_cplex
 use m_slk,               only : matrix_scalapack, processor_scalapack, &
                                 compute_eigen_problem, compute_generalized_eigen_problem
 use m_kg,                only : mkkin, mkkpg
 use m_crystal,           only : crystal_t
 use m_fftcore,           only : kpgsph, get_kg
 use m_fft,               only : fftpac
 use m_cgtools,           only : set_istwfk
 use m_electronpositron,  only : electronpositron_type
 use m_mpinfo,            only : destroy_mpi_enreg, initmpi_seq
 use m_pawtab,            only : pawtab_type
 use m_paw_ij,            only : paw_ij_type
 use m_pawcprj,           only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_reorder
 use m_pawfgr,            only : pawfgr_type
 use m_initylmg,          only : initylmg
 use m_mkffnl,            only : mkffnl
 use m_getghc,            only : getghc, multithreaded_getghc
 use m_cgprj,             only : getcprj

 implicit none

 private
!!***


!!****t* m_ksdiago/ugb_t
!! NAME
!!  ugb_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: ugb_t
   integer :: istwf_k = -1
   integer :: nband_k = - 1
   integer :: npw_k = -1
   integer :: nspinor = -1
   type(processor_scalapack) :: processor
   type(matrix_scalapack) :: mat
   integer, allocatable :: kg_k(:,:)
   real(dp),allocatable :: eig_k(:)
   !type(pawcprj_type),allocatable :: cprj_k(:,:)
   ! (natom, nspinor*onband_diago))

 !contains
 !  procedure :: free =>  ugb_free
 end type ugb_t
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

  integer :: isppol
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
   ! used fro RANGE="V","I", and "A" when do_full_diago=.FALSE.
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

 public :: ksdiago
 public :: init_ddiago_ctl
 public :: ksdiago_slk
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
 integer,parameter :: mkmem_ = 1, tim_getghc = 4, paral_kgb0 = 0, master = 0, ndat1 = 1, ncomp1 = 1
 integer :: cprj_choice,cpopt,dimffnl,ib,ider,idir,isppol,npw_k
 integer :: ikg,istwf_k,exchn2n3d,prtvol
 integer :: jj,n1,n2,n3,n4,n5,n6,negv,nkpg,nproc,npw_k_test,my_rank,optder
 integer :: type_calc,sij_opt,igsp2,cplex_ghg,iband,iorder_cprj,ibs1,ibs2
 real(dp),parameter :: lambda0 = zero
 real(dp) :: ucvol,ecutsm,effmass_free,size_mat,ecut
 logical :: do_full_diago
 character(len=50) :: jobz,range
 character(len=80) :: frmt1,frmt2
 character(len=10) :: stag(2)
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer :: nloalg(3)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),kptns_(3,1),kpoint(3),ylmgr_dum(1,1,1)
 real(dp),allocatable :: ph3d(:,:,:),bras(:,:),ffnl(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),ylm_k(:,:),dum_ylm_gr_k(:,:,:)
 real(dp),allocatable :: ghc(:,:),gvnlxc(:,:),gsc(:,:),ghg_mat(:,:,:),gsg_mat(:,:,:)
 real(dp),pointer :: cwavef(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 if (nproc > 1) then
   ABI_WARNING("ksdiago not supported in parallel. Running in sequential.")
 end if

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type for sequential part.
 call init_distribfft_seq(MPI_enreg_seq%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 if (pawfgr%usefinegrid /= 0) then
   call init_distribfft_seq(MPI_enreg_seq%distribfft, 'f', pawfgr%ngfft(2), pawfgr%ngfft(3), 'all')
 end if

 isppol  = Diago_ctl%isppol
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

 if (nsppol == 1) stag= ['          ','          ']
 if (nsppol == 2) stag= ['SPIN UP:  ','SPIN DOWN:']

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
   call init_hamiltonian(gs_hamk, psps, pawtab, nspinor, nsppol, nspden, natom, typat, xred, nfftc, &
    mgfftc, ngfftc, rprimd, nloalg, paw_ij=paw_ij, usecprj=0, electronpositron=electronpositron)
 else
   call init_hamiltonian(gs_hamk, psps, pawtab, nspinor, nsppol, nspden, natom, typat, xred, nfftc, &
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

 ! Set up local potential vlocal on the coarse FFT mesh from vtrial taking into account the spin.

 call gspot_transgrid_and_pack(isppol, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
                               nspden, gs_hamk%nvloc, ncomp1, pawfgr, mpi_enreg_seq, vtrial, vlocal)
 call gs_hamk%load_spin(isppol, vlocal=vlocal, with_nonlocal=.true.)

 ! This for meta-gga.
 !if (with_vxctau) then
 !  call gspot_transgrid_and_pack(isppol, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
 !                                nspden, gs_hamk%nvloc, 4, pawfgr, mpi_enreg, vxctau, vxctaulocal)
 !  call gs_hamk%load_spin(isppol, vxctaulocal=vxctaulocal)
 !end if

 ! Calculate G-vectors, for this k-point. Count also the number of planewaves as a check.
 exchn2n3d = 0; ikg = 0
 ABI_MALLOC(kg_k, (3, npw_k))

 call kpgsph(ecut, exchn2n3d, gmet, ikg, 0, istwf_k, kg_k, kpoint, 0, MPI_enreg_seq, 0, npw_k_test)
 ABI_CHECK(npw_k_test == npw_k, "npw_k_test/=npw_k")
 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,mkmem_,MPI_enreg_seq,npw_k,npw_k_test)

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
   call initylmg(gprimd, kg_k, kptns_, mkmem_, MPI_enreg_seq, psps%mpsang, npw_k, [nband_k], 1, &
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

 ! cwaveprj is ordered, see nonlop_ylm.
 ABI_MALLOC(cwaveprj, (natom, nspinor*(1+cpopt)*gs_hamk%usepaw))
 if (cpopt == 0) call pawcprj_alloc(cwaveprj, 0, gs_hamk%dimcprj)

 ! Initialize plane-wave array with zeros
 ABI_CALLOC(bras, (2, npw_k*nspinor))
 if (prtvol > 0) call wrtout(std_out, ' Calculating <G|H|G''> elements')

 ! Loop over the |beta,G''> component.
 do igsp2=1,npw_k*nspinor

   bras(1, igsp2) = one

   ! Get <:|H|beta,G''> and <:|S_{PAW}|beta,G''>
   call getghc(cpopt, bras, cwaveprj, ghc, gsc, gs_hamk, gvnlxc, lambda0, MPI_enreg_seq, ndat1, &
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
     ' Begin full diagonalization for kpt: ',trim(ktoa(kpoint)), stag(isppol), ch10,&
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
     ' Begin partial diagonalization for kpt= ',kpoint,stag(isppol),ch10,&
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
   frmt1='(8x,9(1x,f7.2))' ; frmt2='(8x,9(1x,f7.2))'
   write(msg,'(2a,3x,a)')' Eigenvalues in eV for kpt: ', trim(ktoa(kpoint)), stag(isppol)
   call wrtout(std_out, msg)

   write(msg,frmt1)(eig_ene(ib)*Ha_eV,ib=1,MIN(9,onband_diago))
   call wrtout(std_out, msg)
   if (onband_diago >9 ) then
     do jj=10,onband_diago,9
       write(msg,frmt2) (eig_ene(ib)*Ha_eV,ib=jj,MIN(jj+8,onband_diago))
       call wrtout(std_out, msg)
     end do
   end if
 end if

 !========================================================
 !==== Calculate <Proj_i|Cnk> from output eigenstates ====
 !========================================================
 if (psps%usepaw == 1) then

   iorder_cprj = 1 !  Ordered (order does change wrt input file); will be changed later
   ABI_MALLOC(cprj_k,(natom, nspinor*onband_diago))
   call pawcprj_alloc(cprj_k, 0, gs_hamk%dimcprj)

   idir = 0; cprj_choice = 1  ! Only projected wave functions.

   do iband=1,onband_diago
     ibs1=nspinor*(iband-1)+1
     ibs2=ibs1; if (nspinor==2) ibs2=ibs2+1
     cwavef => eig_vec(1:2,1:npw_k,iband)

     call getcprj(cprj_choice,0,cwavef,cprj_k(:,ibs1:ibs2),&
      gs_hamk%ffnl_k,idir,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,&
      gs_hamk%kpg_k,gs_hamk%kpt_k,gs_hamk%lmnmax,gs_hamk%mgfft,MPI_enreg_seq,&
      gs_hamk%natom,gs_hamk%nattyp,gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,&
      gs_hamk%ntypat,gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
   end do

   !  Reorder the cprj (order is now the same as in input file)
   call pawcprj_reorder(cprj_k, gs_hamk%atindx1)
 end if ! usepaw

 ! Free memory.
 ABI_FREE(kpg_k)
 ABI_FREE(kg_k)
 ABI_FREE(ph3d)
 ABI_FREE(ffnl)

 call destroy_mpi_enreg(MPI_enreg_seq)
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

subroutine init_ddiago_ctl(Dctl, jobz, isppol, nspinor, ecut, kpoint, nloalg, gmet, &
  nband_k, istwf_k, ecutsm, effmass_free, abstol, range, ilu, vlu, use_scalapack, prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,nspinor
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
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: kg_k(:,:)

! *************************************************************************

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type.

 Dctl%isppol  = isppol
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
 call kpgsph(ecut,0,gmet,0,0,1,kg_k,kpoint,0,MPI_enreg_seq,0,Dctl%npwtot)

 ! G-vectors taking into account time-reversal symmetry.
 call kpgsph(ecut,0,gmet,0,0,istwf_k,kg_k,kpoint,0,MPI_enreg_seq,0,npw_k)

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

 call destroy_mpi_enreg(MPI_enreg_seq)

end subroutine init_ddiago_ctl
!!***

!!****f* m_ksdiago/ksdiago_slk
!! NAME
!! ksdiago_slk
!!
!! FUNCTION
!!  This routine performs the direct diagonalization of the Kohn-Sham Hamiltonian
!!  for a given k-point and spin using Scalapack/ELPA
!!
!! INPUTS
!!  kpoint(3)
!!  prtvol=Integer Flags  defining verbosity level
!!  mgfftc=maximum size of 1D FFTs (coarse mesh).
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
!!  onband_diago
!!  eig_ene(onband_diago)=The calculatated eigenvalues in ascending order.
!!  cprj_k(natom,nspinor*onband_diago): Projected eigenstates <Proj_i|Cnk> from output eigenstates. PAW only
!!
!! SOURCE

subroutine ksdiago_slk(isppol, istwf_k, kpoint, ecut, nband_k, ngfftc, nfftf, &
                       dtset, pawtab, pawfgr, paw_ij, cryst, psps, vtrial, &
                       onband_diago, eig_ene, ugb, cprj_k, comm, &
                       electronpositron) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol, istwf_k
 real(dp),intent(in) :: kpoint(3), ecut
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: comm,nband_k,nfftf
 integer,intent(inout) :: onband_diago
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: ngfftc(18) !, ngfftf(18)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),allocatable,intent(inout) :: eig_ene(:)
 type(matrix_scalapack),intent(inout) :: ugb
 type(pawcprj_type),allocatable,intent(inout) :: cprj_k(:,:)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij(cryst%natom*psps%usepaw)
 type(electronpositron_type),optional,pointer :: electronpositron

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_ = 1, tim_getghc = 4, paral_kgb0 = 0, master = 0, ncomp1 = 1
 integer :: cprj_choice,cpopt,dimffnl,ib,ider,idir,npw_k,nfftc,mgfftc, igs, ige, omp_nt
 integer :: jj,n1,n2,n3,n4,n5,n6,negv,nkpg,nproc,npw_k_test,my_rank,optder
 integer :: type_calc,sij_opt,igsp2,ig, cplex_ghg,iband,iorder_cprj,ibs1,ibs2
 integer :: npwsp, col_bsize, nsppol, nspinor, nspden, loc2_size, il_g2
 integer :: idat, ndat, batch_size, mat_size !, ierr
 real(dp),parameter :: lambda0 = zero
 real(dp) :: size_mat, cpu, wall, gflops !, mem_mb
 logical :: do_full_diago
 character(len=80) :: frmt1,frmt2
 character(len=10) :: stag(2)
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
 type(gs_hamiltonian_type) :: gs_hamk
 type(matrix_scalapack) :: ghg_mat, gsg_mat, ghg_4diag, gsg_4diag
 type(processor_scalapack) :: slkproc, slkproc_4diag
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kptns_(3,1), ylmgr_dum(1,1,1)
 real(dp),allocatable :: ph3d(:,:,:),ffnl(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),ylm_k(:,:),dum_ylm_gr_k(:,:,:)
 real(dp),target,allocatable :: bras(:,:), ghc(:,:), gvnlxc(:,:), gsc(:,:) !,ghg_mat(:,:,:),gsg_mat(:,:,:)
 real(dp),pointer :: cwavef(:,:)
 !real(dp), ABI_CONTIGUOUS pointer :: gsc_bk(:,:), cg_bk(:,:), ghc_bk(:,:), gvnlxc_bk(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 !====================
 !=== Check input ====
 !====================
 if (all(istwf_k /= [1, 2])) then
   ABI_ERROR(sjoin("istwfk:", itoa(istwf_k), "not allowed:"))
 end if

 if (istwf_k == 2) then
   !ABI_WARNING("istwfk == 2 is still under development")
   ABI_ERROR("istwfk == 2 is still under development")
 end if

 ! MPI_type for sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 if (pawfgr%usefinegrid /= 0) then
   call init_distribfft_seq(MPI_enreg_seq%distribfft, 'f', pawfgr%ngfft(2), pawfgr%ngfft(3), 'all')
 end if

 nspinor = dtset%nspinor
 nsppol = dtset%nsppol
 nspden = dtset%nspden
 if (nsppol == 1) stag= ['          ','          ']
 if (nsppol == 2) stag= ['SPIN UP:  ','SPIN DOWN:']

 call get_kg(kpoint, istwf_k, ecut, cryst%gmet, npw_k, kg_k)
 npwsp = npw_k * nspinor

 ! The coarse FFT mesh.
 n1 = ngfftc(1); n2 = ngfftc(2); n3 = ngfftc(3)
 n4 = ngfftc(4); n5 = ngfftc(5); n6 = ngfftc(6)
 nfftc = product(ngfftc(1:3))
 mgfftc = maxval(ngfftc(1:3))

 ! Initialize the Hamiltonian datatype on the coarse FFT mesh.
 if (present(electronpositron)) then
   call init_hamiltonian(gs_hamk, psps, pawtab, nspinor, nsppol, nspden, cryst%natom, cryst%typat, cryst%xred, nfftc, &
    mgfftc, ngfftc, cryst%rprimd, dtset%nloalg, paw_ij=paw_ij, usecprj=0, electronpositron=electronpositron)
 else
   call init_hamiltonian(gs_hamk, psps, pawtab, nspinor, nsppol, nspden, cryst%natom, cryst%typat, cryst%xred, nfftc, &
    mgfftc, ngfftc, cryst%rprimd, dtset%nloalg, paw_ij=paw_ij, usecprj=0)
 end if

 ! Check on the number of stored bands.
 !onband_diago = nband_k
 !if (nband_k == -1 .or. nband_k >= npwsp) then
 !  onband_diago = npwsp
 !  write(msg,'(4a, i0)')ch10,&
 !   ' Since the number of bands to be computed was -1 or',ch10,&
 !   ' too large, it has been set to the maximum value npw_k*nspinor: ',npwsp
 !  call wrtout(std_out, msg)
 !end if

 !do_full_diago = (onband_diago == npw_k*nspinor)
 do_full_diago = .True.

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
 if (dtset%prtvol > 0) call wrtout(std_out, msg)

 ! Set up local potential vlocal with proper dimensioning, from vtrial.
 ! Select spin component of interest if nspden<=2 as nvloc==1, for nspden==4, nvloc==4
 ! option=2: vtrial(n1*n2*n3,ispden) --> vlocal(nd1,nd2,nd3) real case

 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamk%nvloc))

 ! Set up local potential vlocal on the coarse FFT mesh from vtrial taking into account the spin.

 call gspot_transgrid_and_pack(isppol, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
                               nspden, gs_hamk%nvloc, ncomp1, pawfgr, mpi_enreg_seq, vtrial, vlocal)
 call gs_hamk%load_spin(isppol, vlocal=vlocal, with_nonlocal=.true.)

 ! This for meta-gga.
 !if (with_vxctau) then
 !  call gspot_transgrid_and_pack(isppol, psps%usepaw, paral_kgb0, nfftc, ngfftc, nfftf, &
 !                                nspden, gs_hamk%nvloc, 4, pawfgr, mpi_enreg, vxctau, vxctaulocal)
 !  call gs_hamk%load_spin(isppol, vxctaulocal=vxctaulocal)
 !end if

 !========================
 !==== Kinetic energy ====
 !========================
 ABI_MALLOC(kinpw, (npw_k))
 call mkkin(ecut, dtset%ecutsm, dtset%effmass_free, cryst%gmet, kg_k, kinpw, kpoint, npw_k, 0, 0)

 !================================
 !==== Non-local form factors ====
 !================================
 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2*psps%useylm))

 if (psps%useylm == 1) then
   optder = 0
   ABI_MALLOC(dum_ylm_gr_k, (npw_k, 3+6*(optder/2),psps%mpsang**2))
   kptns_(:,1) = kpoint

   ! Here mband is not used if paral_compil_kpt=0
   call initylmg(cryst%gprimd, kg_k, kptns_, mkmem_, MPI_enreg_seq, psps%mpsang, npw_k, [nband_k], 1, &
     [npw_k], 1, optder, cryst%rprimd, ylm_k, dum_ylm_gr_k)

   ABI_FREE(dum_ylm_gr_k)
 end if

 ! Compute (k+G) vectors (only if useylm=1)
 nkpg = 3 * dtset%nloalg(3)
 ABI_MALLOC(kpg_k, (npw_k, nkpg))
 if (nkpg > 0) call mkkpg(kg_k, kpg_k, kpoint, nkpg, npw_k)

 ! Compute nonlocal form factors ffnl at all (k+G):
 idir=0; ider=0; dimffnl=1+ider ! Now the derivative is not needed anymore.
 ABI_MALLOC(ffnl, (npw_k, dimffnl, psps%lmnmax, psps%ntypat))

 call mkffnl(psps%dimekb, dimffnl, psps%ekb, ffnl, psps%ffspl, cryst%gmet, cryst%gprimd, ider, idir, psps%indlmn, &
   kg_k, kpg_k, kpoint, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, &
   psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum)

 ABI_FREE(ylm_k)

 ! Load k-dependent part in the Hamiltonian datastructure
 ABI_MALLOC(ph3d, (2, npw_k, gs_hamk%matblk))
 call gs_hamk%load_k(kpt_k=kpoint, istwf_k=istwf_k, npw_k=npw_k, kinpw_k=kinpw, &
                     kg_k=kg_k, kpg_k=kpg_k, ffnl_k=ffnl, ph3d_k=ph3d, compute_ph3d=.true., compute_gbound=.true.)

 ! Prepare call to getghc.
 type_calc = 0                                ! For applying the whole Hamiltonian
 sij_opt = 0; if (psps%usepaw==1) sij_opt = 1 ! For PAW, <k+G|S|k+G"> is also needed.

 cpopt = -1    ! If cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
 if (psps%usepaw == 1 .and. .FALSE.) then ! TODO Calculate <p_lmn|k+G>.
   cpopt = 0  ! <p_lmn|in> are computed here and saved
 end if

 ! Define batch size for the application of the Hamiltonian
 ! This is useful if OpenMP is activated thus we use multiple of omp_nt.
 omp_nt = xomp_get_num_threads(open_parallel=.True.)
 batch_size = 1
 !batch_size = 4 * omp_nt

 ABI_MALLOC(ghc, (2, npwsp * batch_size))
 ABI_MALLOC(gvnlxc, (2, npwsp * batch_size))
 ABI_MALLOC(gsc, (2, npwsp * batch_size*(sij_opt+1)/2))

 ! Init 1D PBLAS grid to block-distribute matrices along columns.
 call slkproc%init(comm, grid_dims=[1, nproc])

 mat_size = npwsp; if (istwf_k == 2) mat_size = 2 * npwsp - 1

 col_bsize = mat_size / nproc; if (mod(mat_size, nproc) /= 0) col_bsize = col_bsize + 1
 call ghg_mat%init(mat_size, mat_size, slkproc, istwf_k, size_blocs=[mat_size, col_bsize])
 if (psps%usepaw == 1) call gsg_mat%init(mat_size, mat_size, slkproc, istwf_k, size_blocs=[mat_size, col_bsize])

 ! cwaveprj is ordered, see nonlop_ylm.
 ABI_MALLOC(cwaveprj, (cryst%natom, nspinor*(1+cpopt)*gs_hamk%usepaw*batch_size))
 if (cpopt == 0) call pawcprj_alloc(cwaveprj, 0, gs_hamk%dimcprj)

 ! Initialize plane-wave array with zeros
 ABI_CALLOC(bras, (2, npwsp * batch_size))
 if (dtset%prtvol > 0) call wrtout(std_out, ' Calculating <G|H|G''> elements')

 ! Loop over the |beta,G''> component.
 loc2_size = ghg_mat%sizeb_local(2)

 do il_g2=1, ghg_mat%sizeb_local(2), batch_size
   ndat = merge(batch_size, loc2_size-il_g2+1, il_g2+batch_size-1 <= loc2_size)
   igsp2 = ghg_mat%loc2gcol(il_g2)

   bras = zero
   if (istwf_k == 1) then
     do idat=0,ndat-1
       bras(1, igsp2 + idat * npwsp + idat) = one
     end do
   else
     do idat=0,ndat-1
       if (igsp2 + idat <= npwsp) then
         bras(1, igsp2 + idat * npwsp + idat) = half   ! Cosine
         if (igsp2 == 1) bras(1, igsp2 + idat * npwsp + idat) = one
       else
         ig = igsp2 - npwsp + 1 + 1
         bras(2, ig + idat * npwsp + idat) = +half  ! Sine
       end if
     end do
   end if

   ! Get <:|H|beta,G''> and <:|S_{PAW}|beta,G''>
   call multithreaded_getghc(cpopt, bras, cwaveprj, ghc, gsc, gs_hamk, gvnlxc, lambda0, mpi_enreg_seq, ndat, &
                             dtset%prtvol, sij_opt, tim_getghc, type_calc)

   ! Fill local buffer.
   if (istwf_k == 1) then
     do idat=0,ndat-1
       igs = 1 + idat * npwsp; ige = igs + npwsp - 1
       ghg_mat%buffer_cplx(:, il_g2 + idat) = cmplx(ghc(1, igs:ige), ghc(2, igs:ige), kind=dp)
     end do
     if (psps%usepaw == 1) then
       do idat=0,ndat-1
         igs = 1 + idat * npwsp; ige = igs + npwsp - 1
         gsg_mat%buffer_cplx(:,il_g2+idat) = cmplx(gsc(1,igs:ige), gsc(2,igs:ige), kind=dp)
       end do
     end if

   else

     do idat=0,ndat-1
       igs = 1 + idat * npwsp; ige = igs + npwsp - 1
       !if (igsp2 == 1 .or. igsp2 == npwsp + 1 .and. idat == 0) then
       !  ghc(:, igs:ige) = tol3
       !  !print *, ghc(:, igs:ige)
       !end if
       ghg_mat%buffer_real(1:npwsp, il_g2 + idat) = ghc(1, igs:ige)    ! CC or CS
       !ghg_mat%buffer_real(npwsp+1:, il_g2 + idat) = -ghc(2, igs:ige)  ! SC or SS
       ghg_mat%buffer_real(npwsp+1:, il_g2 + idat) = -ghc(2, igs+1:ige)  ! SC or SS
     end do
   end if

   ! Reset the |G,beta> components that has been treated.
   !do idat=0,ndat-1
   !  bras(1, igsp2 + idat * npwsp + idat) = zero
   !end do
 end do ! il_g2

 ! Free workspace memory allocated so far.
 ABI_FREE(bras)
 ABI_FREE(kinpw)
 ABI_FREE(vlocal)
 ABI_FREE(ghc)
 ABI_FREE(gvnlxc)
 ABI_FREE(gsc)

 if (psps%usepaw == 1 .and. cpopt == 0) call pawcprj_free(Cwaveprj)
 ABI_FREE(cwaveprj)

 !===========================================
 !=== Diagonalization of <G|H|G''> matrix ===
 !===========================================
 onband_diago = mat_size
 ABI_MALLOC(eig_ene, (mat_size))

 ! If possible, use 2D rectangular grid of processors for diagonalization.
 call slkproc_4diag%init(comm)

 call ghg_mat%change_size_blocs(ghg_4diag, processor=slkproc_4diag)
 call ghg_mat%free()
 if (psps%usepaw == 1) then
   call gsg_mat%change_size_blocs(gsg_4diag, processor=slkproc_4diag)
   call gsg_mat%free()
 end if
 call ghg_4diag%copy(ugb)

 if (do_full_diago) then
   ! Full diagonalization
   write(msg,'(6a, i0)')ch10,&
     ' Begin full diagonalization for kpt: ',trim(ktoa(kpoint)), stag(isppol), ch10,&
     ' Matrix size: ',npwsp
   call wrtout(std_out, msg)
   call cwtime(cpu, wall, gflops, "start")

   if (psps%usepaw == 0) then
     !call ghg_4diag%pzheev("V", "U", ugb, eig_ene)
     call compute_eigen_problem(ghg_4diag%processor, ghg_4diag, ugb, eig_ene, comm, istwf_k)
   else
     call compute_generalized_eigen_problem(ghg_4diag%processor, ghg_4diag, gsg_4diag, ugb, eig_ene, comm, istwf_k)
   end if
   call cwtime_report(" full_diago", cpu, wall, gflops)

 else
   ! Partial diagonalization
   ! range = Diago_ctl%range !range="Irange"
   ABI_ERROR("Not Implemented Error")
   call cwtime(cpu, wall, gflops, "start")

   write(msg,'(2a,3es16.8,3a,i0,a,i0)')ch10,&
     ' Begin partial diagonalization for kpt= ',kpoint,stag(isppol),ch10,&
     ' - Size of mat.=',npw_k*nspinor,' - # out_nband: ', onband_diago
   call wrtout(std_out, msg)

   !range = "I"
   !ilu =
   if (psps%usepaw == 0) then
     !call slk_pzheevx(Slk_mat, "V", range, "U", vl, vu, il, iu, abstol, Slk_vec, mene_found, eigen)
     !call xheevx_cplex(jobz, range, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, zero, zero,&
     !  1, onband_diago, -tol8, negv, eig_ene, eig_vec, npw_k*nspinor, msg, ierr)
   else
     !call slk_pzhegvx(Slk_matA, ibtype, "V", range, "U", Slk_matB, vl, vu, il, iu, abstol, Slk_vec, mene_found, eigen)
     !call xhegvx_cplex(1, jobz, range, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, gsg_mat, zero, zero,&
     !  1, onband_diago, -tol8, negv, eig_ene, eig_vec, npw_k*nspinor, msg, ierr)
   end if
   call cwtime_report(" partial_diago", cpu, wall, gflops)
 end if

 call ghg_4diag%free()
 call gsg_4diag%free()
 call slkproc%free()
 !call ugb%free()
 !call slkproc_4diag%free()

 if (dtset%prtvol > 0 .and. my_rank == master) then
   ! Write eigenvalues.
   frmt1='(8x,9(1x,f7.2))' ; frmt2='(8x,9(1x,f7.2))'
   write(msg,'(2a,3x,a)')' Eigenvalues in eV for kpt: ', trim(ktoa(kpoint)), stag(isppol)
   call wrtout(std_out, msg)

   write(msg,frmt1)(eig_ene(ib)*Ha_eV,ib=1,MIN(9,onband_diago))
   call wrtout(std_out, msg)
   if (onband_diago > 9 ) then
     do jj=10,onband_diago,9
       write(msg,frmt2) (eig_ene(ib)*Ha_eV,ib=jj,MIN(jj+8,onband_diago))
       call wrtout(std_out, msg)
     end do
   end if
 end if

 !========================================================
 !==== Calculate <Proj_i|Cnk> from output eigenstates ====
 !========================================================
 !if (psps%usepaw == 1) then

 !  iorder_cprj = 1 !  Ordered (order does change wrt input file); will be changed later
 !  ABI_MALLOC(cprj_k, (cryst%natom, nspinor*onband_diago))
 !  call pawcprj_alloc(cprj_k, 0, gs_hamk%dimcprj)
 !  idir = 0; cprj_choice = 1  ! Only projected wave functions.

 !  do iband=1,onband_diago
 !    ibs1=nspinor*(iband-1)+1
 !    ibs2=ibs1; if (nspinor==2) ibs2=ibs2+1
 !    cwavef => eig_vec(1:2,1:npw_k,iband)

 !    call getcprj(cprj_choice,0,cwavef,cprj_k(:,ibs1:ibs2),&
 !     gs_hamk%ffnl_k,idir,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,&
 !     gs_hamk%kpg_k,gs_hamk%kpt_k,gs_hamk%lmnmax,gs_hamk%mgfft,MPI_enreg_seq,&
 !     gs_hamk%natom,gs_hamk%nattyp,gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,&
 !     gs_hamk%ntypat,gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
 !  end do

 !  !  Reorder the cprj (order is now the same as in input file)
 !  call pawcprj_reorder(cprj_k, gs_hamk%atindx1)
 !end if ! usepaw

 ! Free memory.
 ABI_FREE(kpg_k)
 ABI_FREE(ph3d)
 ABI_FREE(ffnl)
 ABI_FREE(kg_k)

 call destroy_mpi_enreg(MPI_enreg_seq)
 call gs_hamk%free()

 call xmpi_barrier(comm)

end subroutine ksdiago_slk
!!***

end module m_ksdiago
!!***

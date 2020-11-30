!!****m* ABINIT/m_ksdiago
!! NAME
!!  m_ksdiago
!!
!! FUNCTION
!!  Direct diagonalization of the KS Hamiltonian H_k(G,G')
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
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

module m_ksdiago

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_hamiltonian
 use m_distribfft

 use defs_datatypes,      only : pseudopotential_type
 use defs_abitypes,       only : MPI_type
 use m_fstrings,          only : toupper
 use m_geometry,          only : metric
 use m_hide_lapack,       only : xhegv_cplex, xheev_cplex, xheevx_cplex, xhegvx_cplex
 use m_kg,                only : mkkin, mkkpg
 use m_fftcore,           only : kpgsph
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
 use m_getghc,            only : getghc
 use m_fourier_interpol,  only : transgrid
 use m_cgprj,             only : getcprj

 implicit none

 private
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

 public ::  ksdiago
 public ::  init_ddiago_ctl
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
!!  At present, only norm-conserving pseudopotentials are implemented.
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
!!  Pawtab(Psps%ntypat*Psps%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  Pawfgr<pawfgr_type>=fine grid parameters and related data
!!  Paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  vtrial(nfftf,nspden)=the trial potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  comm=MPI communicator.
!!  [Electronpositron] <electronpositron_type>=quantities for the electron-positron annihilation.
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
!!  Cprj_k(natom,nspinor*onband_diago) PAW only===
!!   input: pointer to NULL
!!   output: Projected eigenstates <Proj_i|Cnk> from output eigenstates.
!!
!! NOTES
!! * The routine can be time consuming (in particular when computing <G1|H|G2> elements for all (G1,G2)).
!!   So, it is recommended to call it once per run.
!!
!! * The routine RE-compute all Hamiltonian terms. So it is equivalent to an additional electronic SC cycle.
!!   (This has no effect is convergence was reach. If not, eigenvalues/vectors may differs from the conjugate gradient ones)
!!
!! * Please, do NOT pass Dtset% to this routine. Either use a local variable properly initialized
!!   or add the additional variable to ddiago_ctl_type and change the creation method accordingly.
!!   ksdiago is designed such that it is possible to diagonalize the Hamiltonian at an arbitrary k-point
!!   or spin (not efficient but easy to code). Therefore ksdiago is useful non only for
!!   the KSS generation but also for testing more advanced iterative algorithms as well as interpolation techniques.
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine ksdiago(Diago_ctl,nband_k,nfftc,mgfftc,ngfftc,natom,&
& typat,nfftf,nspinor,nspden,nsppol,Pawtab,Pawfgr,Paw_ij,&
& Psps,rprimd,vtrial,xred,onband_diago,eig_ene,eig_vec,Cprj_k,comm,ierr,&
& Electronpositron) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfftc,natom,comm,nband_k,nfftf,nsppol,nspden,nspinor,nfftc
 integer,intent(out) :: ierr, onband_diago
 type(pseudopotential_type),intent(in) :: Psps
 type(pawfgr_type),intent(in) :: Pawfgr
 type(ddiago_ctl_type),intent(in) :: Diago_ctl
!arrays
 integer,intent(in) :: typat(natom), ngfftc(18)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: vtrial(nfftf,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),pointer :: eig_ene(:),eig_vec(:,:,:)
 type(pawcprj_type),pointer :: Cprj_k(:,:)
 type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(paw_ij_type),intent(in) :: Paw_ij(natom*Psps%usepaw)
 type(electronpositron_type),optional,pointer :: Electronpositron

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_=1,tim_getghc=4,paral_kgb=0,master=0,ndat1=1
 integer :: cprj_choice,cpopt,dimffnl,ib,ider,idir,isppol,npw_k
 integer :: ikg,istwf_k,exchn2n3d,prtvol
 integer :: jj,n1,n2,n3,n4,n5,n6,negv,nkpg,nprocs,npw_k_test,my_rank,optder
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
 real(dp),allocatable :: ph3d(:,:,:),pwave(:,:),ffnl(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),ylm_k(:,:),dum_ylm_gr_k(:,:,:)
 real(dp),allocatable :: ghc(:,:),gvnlxc(:,:),gsc(:,:),ghg_mat(:,:,:),gsg_mat(:,:,:)
 real(dp),pointer :: cwavef(:,:)
 type(pawcprj_type),allocatable :: Cwaveprj(:,:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 if (nprocs > 1) then
   MSG_WARNING("ksdiago not supported in parallel. Running in sequential.")
 end if

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type for sequential part.
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 if (Pawfgr%usefinegrid /= 0) then
   call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',Pawfgr%ngfft(2),pawfgr%ngfft(3),'all')
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
 ierr=0

 ! istwfk must be 1 for each k-point
 if (istwf_k/=1) then
   write(msg,'(7a)')&
   ' istwfk /= 1 not allowed:',ch10,&
   ' States output not programmed for time-reversal symmetry.',ch10,&
   ' Action: change istwfk in input file (put it to 1 for all kpt).',ch10,&
   ' Program does not stop but _KSS file will not be created...'
   MSG_WARNING(msg)
   ierr=ierr+1
 end if

 if (paral_kgb /= 0) then
   write(msg,'(3a)')&
   ' paral_kgb /= 0 not allowed:',ch10,&
   ' Program does not stop but _KSS file will not be created...'
   MSG_WARNING(msg)
   ierr=ierr+1
 end if

 if (ierr /= 0) RETURN ! Houston we have a problem!

 ! Initialize the Hamiltonian datatype on the coarse FFT mesh.
 if (present(Electronpositron)) then
   call init_hamiltonian(gs_hamk,Psps,pawtab,nspinor,nsppol,nspden,natom,typat,xred,nfftc,&
    mgfftc,ngfftc,rprimd,nloalg,paw_ij=Paw_ij,usecprj=0,Electronpositron=Electronpositron)
 else
   call init_hamiltonian(gs_hamk,Psps,pawtab,nspinor,nsppol,nspden,natom,typat,xred,nfftc,&
    mgfftc,ngfftc,rprimd,nloalg,paw_ij=Paw_ij,usecprj=0)
 end if

 ! Check on the number of stored bands.
 onband_diago = nband_k
 if (nband_k==-1 .or. nband_k >= npw_k*nspinor) then
   onband_diago = npw_k*nspinor
   write(msg,'(4a)')ch10,&
    ' Since the number of bands to be computed was (-1) or',ch10,&
    ' too large, it has been set to the max value npw_k*nspinor. '
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

 !nvloc=1; if (nspden==4) nvloc=4
 ABI_MALLOC(vlocal, (n4,n5,n6,gs_hamk%nvloc))

 ! Set up local potential vlocal on the coarse FFT mesh from vtrial taking into account the spin.
 ! Also take into account the spin.

 call gspot_transgrid_and_pack(isppol, psps%usepaw, paral_kgb, nfftc, ngfftc, nfftf, &
                               nspden, gs_hamk%nvloc, 1, pawfgr, mpi_enreg_seq, vtrial, vlocal)
 call gs_hamk%load_spin(isppol, vlocal=vlocal, with_nonlocal=.true.)

 !if (with_vxctau) then
 !  call gspot_transgrid_and_pack(isppol, psps%usepaw, paral_kgb, nfftc, ngfftc, nfftf, &
 !                                nspden, gs_hamk%nvloc, 4, pawfgr, mpi_enreg, vxctau, vxctaulocal)
 !  call gs_hamk%load_spin(isppol, vxctaulocal=vxctaulocal)
 !end if

 ! Calculate G-vectors, for this k-point. Count also the number of planewaves as a check.
 exchn2n3d=0; ikg=0
 ABI_MALLOC(kg_k, (3,npw_k))

 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,0,MPI_enreg_seq,0,npw_k_test)
 ABI_CHECK(npw_k_test == npw_k,"npw_k_test/=npw_k")
 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,mkmem_,MPI_enreg_seq,npw_k,npw_k_test)

 !========================
 !==== Kinetic energy ====
 !========================
 ABI_MALLOC(kinpw, (npw_k))
 call mkkin(ecut, ecutsm, effmass_free, gmet, kg_k, kinpw, kpoint, npw_k, 0, 0)

 !================================
 !==== Non-local form factors ====
 !================================
 ABI_MALLOC(ylm_k, (npw_k,Psps%mpsang**2*Psps%useylm))

 if (Psps%useylm==1) then
   optder=0
   ABI_MALLOC(dum_ylm_gr_k,(npw_k,3+6*(optder/2),Psps%mpsang**2))
   kptns_(:,1) = kpoint

   ! Here mband is not used if paral_compil_kpt=0
   call initylmg(gprimd,kg_k,kptns_,mkmem_,MPI_enreg_seq,Psps%mpsang,npw_k,[nband_k],1,&
     [npw_k],1,optder,rprimd,ylm_k,dum_ylm_gr_k)

   ABI_FREE(dum_ylm_gr_k)
 end if

 ! Compute (k+G) vectors (only if useylm=1)
 nkpg = 3*nloalg(3)
 ABI_MALLOC(kpg_k, (npw_k, nkpg))
 if (nkpg > 0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

 ! Compute nonlocal form factors ffnl at all (k+G):
 idir=0; ider=0; dimffnl=1+ider ! Now the derivative is not needed anymore.
 ABI_MALLOC(ffnl, (npw_k, dimffnl, Psps%lmnmax, Psps%ntypat))

 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,gmet,gprimd,ider,idir,Psps%indlmn,&
   kg_k,kpg_k,kpoint,Psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,&
   Psps%ntypat,Psps%pspso,Psps%qgrid_ff,rmet,Psps%usepaw,Psps%useylm,ylm_k,ylmgr_dum)

 ABI_FREE(ylm_k)

 ! Load k-dependent part in the Hamiltonian datastructure
 ABI_MALLOC(ph3d, (2, npw_k, gs_hamk%matblk))
 call gs_hamk%load_k(kpt_k=kpoint,istwf_k=istwf_k,npw_k=npw_k,kinpw_k=kinpw,&
                     kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)

 ! Prepare call to getghc.
 type_calc=0         ! For applying the whole Hamiltonian
 sij_opt=0; if (Psps%usepaw==1) sij_opt=1 ! For PAW, <k+G|S|k+G"> is also needed.

 cpopt = -1    ! If cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
 if (Psps%usepaw==1.and..FALSE.) then ! TODO Calculate <p_lmn|k+G>.
   cpopt = 0  ! <p_lmn|in> are computed here and saved
 end if

 ABI_MALLOC(ghc, (2, npw_k*nspinor*ndat1))
 ABI_MALLOC(gvnlxc, (2, npw_k*nspinor*ndat1))
 ABI_MALLOC(gsc, (2, npw_k*nspinor*ndat1*(sij_opt+1)/2))

 cplex_ghg = 2
 size_mat = cplex_ghg*(npw_k*nspinor)**2*dp*b2Mb
 write(msg,'(a,f0.3,a)')" Out-of-memory in ghg_mat. Memory required by the Hamiltonian matrix: ",size_mat," [Mb]."
 ABI_STAT_MALLOC(ghg_mat, (cplex_ghg,npw_k*nspinor,npw_k*nspinor), ierr)
 ABI_CHECK(ierr == 0, msg)
 write(msg,'(a,f0.3,a)')" Out-of-memory in gsg_mat. Memory required by the PAW overlap operator: ",size_mat," [Mb]."
 ABI_STAT_MALLOC(gsg_mat, (cplex_ghg,npw_k*nspinor,npw_k*nspinor*Psps%usepaw), ierr)
 ABI_CHECK(ierr == 0, msg)

 ! Cwaveprj is ordered, see nonlop_ylm.
 ABI_MALLOC(Cwaveprj, (natom,nspinor*(1+cpopt)*gs_hamk%usepaw))
 if (cpopt == 0) call pawcprj_alloc(Cwaveprj,0,gs_hamk%dimcprj)

 ! Initialize plane-wave array with zeros
 ABI_CALLOC(pwave, (2, npw_k*nspinor))
 if (prtvol > 0) call wrtout(std_out,' Calculating <G|H|G''> elements')

 ! Loop over the |beta,G''> component.
 do igsp2=1,npw_k*nspinor
   ! Get <:|H|beta,G''> and <:|S_{PAW}|beta,G''>
   pwave(1,igsp2) = one

   call getghc(cpopt,pwave,Cwaveprj,ghc,gsc,gs_hamk,gvnlxc,lambda0,MPI_enreg_seq,ndat1,&
               prtvol,sij_opt,tim_getghc,type_calc)

   ! Fill the upper triangle.
   ghg_mat(:,1:igsp2,igsp2) = ghc(:,1:igsp2)
   if (Psps%usepaw == 1) gsg_mat(:,1:igsp2,igsp2) = gsc(:,1:igsp2)

   ! Reset the |G,beta> component that has been treated.
   pwave(1,igsp2) = zero
 end do

 ! Free workspace memory allocated so far.
 ABI_FREE(pwave)
 ABI_FREE(kinpw)
 ABI_FREE(vlocal)
 ABI_FREE(ghc)
 ABI_FREE(gvnlxc)
 ABI_FREE(gsc)

 if (Psps%usepaw==1.and.cpopt==0) call pawcprj_free(Cwaveprj)
 ABI_FREE(Cwaveprj)

 !===========================================
 !=== Diagonalization of <G|H|G''> matrix ===
 !===========================================
 ABI_MALLOC(eig_ene, (onband_diago))
 ABI_MALLOC(eig_vec, (cplex_ghg, npw_k*nspinor, onband_diago))

 jobz = Diago_ctl%jobz  !jobz="Vectors"

 if (do_full_diago) then
   ! Full diagonalization
   write(msg,'(2a,3es16.8,3x,3a,i5)')ch10,&
     ' Begin complete diagonalization for kpt= ',kpoint(:),stag(isppol),ch10,&
     ' - Size of mat.=',npw_k*nspinor
   if (prtvol>0) call wrtout(std_out, msg)

   if (Psps%usepaw==0) then
     call xheev_cplex(jobz, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, eig_ene, msg, ierr)
   else
     call xhegv_cplex(1, jobz, "Upper", cplex_ghg, npw_k*nspinor, ghg_mat, gsg_mat, eig_ene, msg, ierr)
   end if
   ABI_CHECK(ierr == 0, msg)
   eig_vec(:,:,:)=  ghg_mat

 else
   ! Partial diagonalization
   range = Diago_ctl%range !range="Irange"

   write(msg,'(2a,3es16.8,3a,i5,a,i5)')ch10,&
     ' Begin partial diagonalization for kpt= ',kpoint,stag(isppol),ch10,&
     ' - Size of mat.=',npw_k*nspinor,' - # bnds=',onband_diago
   if (prtvol>0) call wrtout(std_out, msg)

   if (Psps%usepaw==0) then
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
   write(msg,'(a,3es16.8,3x,a)')' Eigenvalues in eV for kpt= ',kpoint,stag(isppol)
   call wrtout(std_out,msg)

   write(msg,frmt1)(eig_ene(ib)*Ha_eV,ib=1,MIN(9,onband_diago))
   call wrtout(std_out,msg)
   if (onband_diago >9 ) then
     do jj=10,onband_diago,9
       write(msg,frmt2) (eig_ene(ib)*Ha_eV,ib=jj,MIN(jj+8,onband_diago))
       call wrtout(std_out,msg)
     end do
   end if
 end if

 !========================================================
 !==== Calculate <Proj_i|Cnk> from output eigenstates ====
 !========================================================
 if (Psps%usepaw == 1) then

   iorder_cprj=1 !  Ordered (order does change wrt input file); will be changed later
   ABI_MALLOC(Cprj_k,(natom,nspinor*onband_diago))
   call pawcprj_alloc(Cprj_k,0,gs_hamk%dimcprj)

   idir=0; cprj_choice=1  ! Only projected wave functions.

   do iband=1,onband_diago
     ibs1=nspinor*(iband-1)+1
     ibs2=ibs1; if (nspinor==2) ibs2=ibs2+1
     cwavef => eig_vec(1:2,1:npw_k,iband)

     call getcprj(cprj_choice,0,cwavef,Cprj_k(:,ibs1:ibs2),&
      gs_hamk%ffnl_k,idir,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,&
      gs_hamk%kpg_k,gs_hamk%kpt_k,gs_hamk%lmnmax,gs_hamk%mgfft,MPI_enreg_seq,&
      gs_hamk%natom,gs_hamk%nattyp,gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,&
      gs_hamk%ntypat,gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
   end do

   !  Reorder the cprj (order is now the same as in input file)
   call pawcprj_reorder(Cprj_k,gs_hamk%atindx1)

 end if !usepaw==1

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
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      destroy_mpi_enreg,initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_ddiago_ctl(Dctl,jobz,isppol,nspinor,ecut,kpoint,nloalg,gmet,&
& nband_k,istwf_k,ecutsm,effmass_free,abstol,range,ilu,vlu,use_scalapack,prtvol)

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

 ABI_CHECK(Dctl%istwf_k==1,"istwf_k/=1 not coded")

 Dctl%jobz   = toupper(jobz(1:1))
 Dctl%range  = "A"
 if (PRESENT(range)) Dctl%range = toupper(range)

 Dctl%ecut = ecut
 Dctl%ecutsm = zero; if (PRESENT(ecutsm)) Dctl%ecutsm = ecutsm
 Dctl%effmass_free = one; if (PRESENT(effmass_free)) Dctl%effmass_free = effmass_free
 Dctl%nloalg  = nloalg
 Dctl%prtvol = 0; if (PRESENT(prtvol)) Dctl%prtvol = prtvol
 Dctl%abstol = -tol8; if (PRESENT(abstol)) Dctl%abstol = abstol

 ABI_ALLOCATE(kg_k,(3,0))

 ! Total number of G-vectors for this k-point with istwf_k=1.
 call kpgsph(ecut,0,gmet,0,0,1,kg_k,kpoint,0,MPI_enreg_seq,0,Dctl%npwtot)

 ! G-vectors taking into account time-reversal symmetry.
 call kpgsph(ecut,0,gmet,0,0,istwf_k,kg_k,kpoint,0,MPI_enreg_seq,0,npw_k)

 Dctl%npw_k = npw_k
 ABI_DEALLOCATE(kg_k)

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
    if (Dctl%prtvol>0) call wrtout(std_out,msg)
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
  if (Dctl%prtvol>0) call wrtout(std_out,msg)

 CASE ("I")
  if (.not.PRESENT(ilu)) then
    MSG_ERROR(" ilu must be specified when range=I ")
  end if
  Dctl%ilu = ilu

  ltest = ( ( ilu(2)>=ilu(1) ) .and. ilu(1)>=1 .and. ilu(2)<=Dctl%npwtot )
  write(msg,'(a,2i0)')" Illegal value for ilu: ",ilu
  ABI_CHECK(ltest,msg)
  Dctl%nband_k= ilu(2)-ilu(1)+1

 CASE ("V")
  if (.not.PRESENT(vlu)) then
    MSG_ERROR(" vlu must be specified when range=V ")
  end if
  Dctl%vlu = vlu

  Dctl%nband_k=-1 !??

  ltest = (vlu(2)>vlu(1))
  write(msg,'(a,2f0.3)')" Illegal value for vlu: ",vlu
  ABI_CHECK(ltest,msg)

 CASE DEFAULT
   MSG_ERROR(" Unknown value for range: "//TRIM(Dctl%range))
 END SELECT

 ! Consider the case in which we asked for the entire set of eigenvectors
 ! but the number of bands is less that npw_k. Therefore have to prepare the call to ZHEEVX.
 ! TODO this has to be done in a cleaner way.
 if (Dctl%range=="A".and. (.not.Dctl%do_full_diago)) then
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

end module m_ksdiago
!!***

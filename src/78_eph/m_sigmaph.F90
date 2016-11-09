!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigmaph
!! NAME
!!  m_sigmaph
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sigmaph

 use defs_basis
 use defs_abitypes
 use iso_c_binding
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_ifc
 use m_ebands
 use m_wfk
 use m_ddb
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfd
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime
 use m_fstrings,       only : itoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth
 use m_io_tools,       only : iomode_from_fname
 use m_special_funcs,  only : dirac_delta
 use m_fftcore,        only : ngfft_seq
 use m_cgtools,        only : dotprod_g !set_istwfk
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_fftcore,        only : get_kg
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
! use m_paw_an,		     only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
! use m_paw_ij,		     only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
! use m_pawfgrtab,	   only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
! use m_pawrhoij,	     only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
! use m_pawdij,		     only : pawdij, symdij

 implicit none

 private
!!***

 public :: sigmaph_driver
!!***

!----------------------------------------------------------------------

!!****t* m_sigmaph/sigmaph_t
!! NAME
!! sigmaph_t
!!
!! FUNCTION
!! Container for the (diagonal) matrix elements of the electronic self-energy (phonon contribution)
!! computed in the Bloch-state representation i.e. Sigma_ph(omega, T, band, k, spin).
!! Provides methods to compute the QP corrections, the spectral functions and save the results to file.
!!
!! TODO
!!  1) Store (q, mu) or q-resolved contributions? More memory but could be useful for post-processing
!!  2) Use list of i.delta values instead of single shift? Useful for convergence studies.
!!  3) Write python code to analyze/merge results.
!!
!! SOURCE

 type,public :: sigmaph_t

  integer :: nkcalc
   ! Number of k-points computed.

  integer :: max_nbcalc
   ! Maximum number of bands computed (max over nkcalc and spin).

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: nspinor
   ! Number of spinorial components.

  !integer :: nsab
   ! Number of spin components in self-energy matrix elements.
   !   1 if (nspinor=1, nsppol=1),
   !   2 if (nspinor=1, nsppol=2),
   !   4 if (nspinor=2)

  integer :: nwr
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Odd number so that the mesh is centered on the KS energy.
   ! The spectral function is computed only if nwr > 1.

  integer :: ntemp
   ! Number of temperatures.

  integer :: symsigma
   ! 1 if matrix elements should be symmetrized.
   ! Rrequired when the sum over q in the BZ is replaced by IBZ(k).

  integer :: timrev
   !  timrev=1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nqbz
  ! Number of q-points in the BZ.

  integer :: ncid
  ! Netcdf file used to save results.

  integer :: mpw
  integer :: gmax(3)

  real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  integer :: ngqpt(3)
  ! Number of divisions in the Q mesh in the BZ.

  integer,allocatable :: bstart_ks(:,:)
   ! bstart_ks(nkcalc,nsppol)
   ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used

  integer,allocatable :: nbcalc_ks(:,:)
   ! nbcalc_ks(nkcalc,nsppol)
   ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used

  integer,allocatable :: kcalc2ibz(:,:)
   ! Mapping kcalc --> ibz as reported in listkk

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.

  real(dp),allocatable :: kTmesh(:)
   ! kTmesh(ntemp)
   ! List of temperatures (kT units).

  real(dp),allocatable :: mu_e(:)
   ! mu_e(ntemp)
   ! chemical potential of electrons for the different temperatures.

  real(dp),allocatable :: wrmesh_b(:,:)
  ! wrmesh_b(nwr, max_nbcalc)
  ! Frequency mesh along the real axis (Ha units) used for the different bands

  complex(dpc),allocatable :: vals_wr(:,:,:)
   ! vals_wr(nwr, ntemp, max_nbcalc)
   ! Sigma_ph(omega, kT, band) for given (k, spin).
   ! enk_KS corresponds to nwr/2 + 1.

  complex(dpc),allocatable :: dvals_dwr(:,:)
   ! vals_wr(ntemp, max_nbcalc)
   ! The derivative of Sigma_ph(omega, kT, band) for given (k, spin) wrt omega computed at the KS energy.
   ! enk_KS corresponds to nwr/2 + 1.

  !complex(dpc),allocatable :: vals_qwr(:,:,:)
   !  vals(nwr, ntemp, max_nbcalc, nqpt, natom*3, nkcalc)
   !  (q, mu)-resolved self-energy for given (k, spin)

  !complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  !integer :: nqibz_k
  ! Number of q-points in the IBZ(k)

  !real(dp),allocatable :: qibz_k(:,:)
  ! qibz(3,nqibz)
  ! Reduced coordinates of the q-points in the IBZ.

  !real(dp),allocatable :: wtq_k(:)
  ! wtq(nqibz)
  ! Weights of the q-points in the IBZ (normalized to one)

  real(dp),allocatable :: qbz(:,:)
  ! qbz(3,nqbz)
  ! Reduced coordinates of the q-points in the BZ.

 end type sigmaph_t

 private :: sigmaph_new             ! Creation method (allocates memory, initialize data from input vars).
 private :: sigmaph_free            ! Free memory.
 !private :: sigmaph_setup_kcalc    ! Return tables used to perform the sum over q-points for given (k-point, spin)
 !private :: sigmaph_setup_kcalc_spin   ! Return tables used to perform the sum over q-points for given (k-point, spin)
 !private :: sigmaph_solve          ! Compute the QP corrections.

 !private :: sigmaph_print          ! Print results to main output file.
 !private :: sigmaph_write          ! Write results to formatted file.
 !private :: sigmaph_symvals        ! Symmetrize self-energy matrix elements (symsigma == 1).

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_driver
!! NAME
!!  sigmaph_driver
!!
!! FUNCTION
!!  Compute phonon-contribution to the electron self-energy.
!!
!! INPUTS
!! wk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_driver(wfk0_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
                          pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_driver'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: dummy_npw=1,tim_getgh1c=1,berryopt0=0
 integer,parameter :: useylmgr=0,useylmgr1=0,master=0,ndat1=1,copt0=0
 integer :: my_rank,nproc,iomode,mband,my_minb,my_maxb,nsppol,nkpt,idir,ipert,iq_ibz
 integer :: cplex,db_iqpt,natom,natom3,ipc,ipc1,ipc2,nspinor
 integer :: ib_kq,ib_k,band,num_smallw
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq !,!timerev_q,
 integer :: spin,istwf_k,istwf_kq,istwf_kirr,npw_k,npw_kq,npw_kirr
 integer :: mpw,ierr,iw,it !ipw
 integer :: n1,n2,n3,n4,n5,n6,nspden,do_ftv1q,imode
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnl1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,nqibz_k
 integer :: nbcalc_ks,bstart_ks,ikcalc
 real(dp) :: cpu,wall,gflops
 real(dp) :: ecut,eshift,dotr,doti,dksqmax
 complex(dpc) :: ieta
 logical,parameter :: have_ktimerev=.True.
 logical :: isirr_k,isirr_kq,gen_eigenpb
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(sigmaph_t) :: sigma
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0bz_kq(3),g0_kq(3),symq(4,2,cryst%nsym),dummy_gvec(3,dummy_npw)
 integer :: work_ngfft(18),gmax(3) !,my_gmax(3) !g0ibz_kq(3),
 integer :: indkk_kq(1,6)
 integer,allocatable :: gtmp(:,:),kg_k(:,:),kg_kq(:,:),nband(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),phfrq(3*cryst%natom),lf(2),rg(2),res(2)
 real(dp) :: f_nk,f_mkq,wqnu,nqnu,gkk2,eig0nk,eig0mkq
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom),displ_red(2,3,cryst%natom,3*cryst%natom)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:)
 real(dp),allocatable :: gkk_atm(:,:,:),gkk_nu(:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:),kets_k(:,:,:),h1kets_kq(:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnl1(:,:),work(:,:,:,:)
 real(dp),allocatable ::  gs1c(:,:)
 !real(dp),allocatable :: wt_kq(:,:)
 real(dp),allocatable :: qibz_k(:,:)
 complex(dpc),allocatable :: fact_wr(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)
 !real(dp),allocatable :: cwave0(:,:),gvnl1(:,:)

!************************************************************************

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkpt = ebands%nkpt; mband=ebands%mband

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx,natom,n1,n2,n3,ph1d,cryst%xred)

 ! Broadening parameter
 if (dtset%elph2_imagden > tol12) then
   ieta = j_dpc * dtset%elph2_imagden
 else
   ieta = j_dpc * 0.0001_dp
 end if

 ! Construct object to store final results.
 ecut = dtset%ecut ! dtset%dilatmx
 sigma = sigmaph_new(dtset, ecut, cryst, ebands, dtfil, comm)

 ! This is the maximum number of PWs for all possible k+q
 mpw = sigma%mpw; gmax = sigma%gmax

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))

 ! Initialize the wave function descriptor.
 ! For the time being, no memory distribution, each node has the full set of states.
 my_minb = 1; my_maxb = mband
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask,(mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur,(mband, nkpt ,nsppol))
 nband=mband; bks_mask=.True.; keep_ur=.False.

 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,dtset%paral_kgb,dummy_npw,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,dtset%ecutsm,dtset%dilatmx,ebands%istwfk,ebands%kptns,ngfft,&
   dummy_gvec,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm,opt_ecut=ecut)

 call wfd_print(wfd,header="Wavefunctions on the Fermi Surface",mode_paral='PERS')
 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)

 iomode = iomode_from_fname(wfk0_path)
 call wfd_read_wfk(wfd,wfk0_path,iomode)
 if (.False.) call wfd_test_ortho(wfd,cryst,pawtab,unit=std_out,mode_paral="PERS")

 ! TODO FOR PAW
 usecprj = 0
 ABI_DT_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1  ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2     ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnl1 = 0 ! gvnl1 is output
 ABI_MALLOC(gvnl1, (2,usevnl))
 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nspden,natom,&
   dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
   usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 !PAW:allocate memory for non-symetrized 1st-order occupancies matrix (pawrhoij1)
 ! pawrhoij1_unsym => pawrhoij1
 ! if (psps%usepaw==1.and.iscf_mod>0) then
 !   if (paral_atom) then
 !     ABI_DATATYPE_ALLOCATE(pawrhoij1_unsym,(natom))
 !     cplex_rhoij=max(cplex,dtset%pawcpxocc);nspden_rhoij=nspden
 !     call pawrhoij_alloc(pawrhoij1_unsym,cplex_rhoij,nspden_rhoij,nspinor,&
 !       dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
 !   else
 !     pawrhoij1_unsym => pawrhoij1
 !     call pawrhoij_init_unpacked(pawrhoij1_unsym)
 !   end if
 ! end if

 ! Allocate work space arrays.
 ! Note nvloc in vlocal.
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 vlocal = huge(one)
 ABI_MALLOC(fact_wr, (sigma%nwr))

 ! Open the DVDB file
 call dvdb_open_read(dvdb, ngfftf, xmpi_comm_self)

 ! Loop over k-points for QP corrections.
 do ikcalc=1,sigma%nkcalc

   ! Symmetry indices for k.
   kk = sigma%kcalc(:, ikcalc)
   ik_ibz = sigma%kcalc2ibz(ikcalc,1); isym_k = sigma%kcalc2ibz(ikcalc,2)
   trev_k = sigma%kcalc2ibz(ikcalc,6); g0_k = sigma%kcalc2ibz(ikcalc, 3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   ABI_CHECK(isirr_k, "k-point must be in IBZ")
   kk_ibz = ebands%kptns(:,ik_ibz)
   npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)

   !call sigmaph_setup_kpoint(sigma, ebands, ikcalc)

   ! For the time being we do not use symmetries --> nqibz_k == nqbz
   nqibz_k = sigma%nqbz
   ABI_MALLOC(qibz_k, (3, nqibz_k))
   qibz_k = sigma%qbz(:,1:sigma%nqbz)
   call wrtout(std_out, sjoin("Will compute", itoa(nqibz_k), "q-points in the IBZ(k)"))

   ! Activate Fourier interpolation if irred q-points are not in the DVDB file.
   ! TODO: This should be done only once
   do_ftv1q = 0
   do iq_ibz=1,nqibz_k
     if (dvdb_findq(dvdb, qibz_k(:,iq_ibz)) == -1) do_ftv1q = do_ftv1q + 1
   end do
   if (do_ftv1q /= 0) then
     write(msg, "(2(a,i0),a)")"Will use Fourier interpolation of DFPT potentials [",do_ftv1q,"/",nqibz_k,"]"
     call wrtout(std_out, msg)
     call wrtout(std_out, sjoin("From ngqpt", ltoa(ifc%ngqpt), "to", ltoa(sigma%ngqpt)))
     call dvdb_ftinterp_setup(dvdb, ifc%ngqpt, 1, [zero,zero,zero], nfft, ngfft, comm)
   end if

   ! Allow PW-arrays dimensioned with mpw
   ABI_MALLOC(kg_k, (3, npw_k))
   kg_k = wfd%kdata(ik_ibz)%kg_k
   ABI_MALLOC(kg_kq, (3, mpw))

   ! Spherical Harmonics for useylm==1.
   ABI_MALLOC(ylm_k,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
   ABI_MALLOC(ylm_kq,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
   ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

   do spin=1,nsppol
     ! Bands to compute.
     bstart_ks = sigma%bstart_ks(ikcalc,spin)
     nbcalc_ks = sigma%nbcalc_ks(ikcalc, spin)

     ! Zero self-energy matrix elements. Build frequency mesh for nk states.
     sigma%vals_wr = zero; sigma%dvals_dwr = zero
     do ib_k=1,nbcalc_ks
       band = ib_k + bstart_ks - 1
       eig0nk = ebands%eig(band,ik_ibz,spin) - sigma%wr_step * (sigma%nwr/2)
       sigma%wrmesh_b(:,ib_k) = arth(eig0nk, sigma%wr_step, sigma%nwr)
     end do
     !call sigmaph_setup_kcalc_spin(sigma, ikcalc, spin, ebands) ??

     ABI_MALLOC(gkk_atm, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkk_nu, (2, nbcalc_ks, natom3))

     ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
     ! TODO: symmetrize them if kk is not irred
     ABI_MALLOC(kets_k, (2, npw_k*nspinor, nbcalc_ks))
     do ib_k=1,nbcalc_ks
       band = ib_k + bstart_ks - 1
       call wfd_copy_cg(wfd, band, ik_ibz, spin, kets_k(1,1,ib_k))
     end do

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal, &
       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

     ! Integrations over q-points.
     do iq_ibz=1,nqibz_k
       !write(std_out,*)"in q-point loopq",iq_ibz,nqibz_k
       call cwtime(cpu,wall,gflops,"start")
       qpt = qibz_k(:,iq_ibz)

       ! Find the index of the q-point in the DVDB.
       db_iqpt = dvdb_findq(dvdb, qpt)

       ! TODO: handle q_bz = S q_ibz case by symmetrizing the potentials in the DVDB.
       if (db_iqpt /= -1) then
         if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
         ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
         ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
         call dvdb_readsym_allv1(dvdb, db_iqpt, cplex, nfftf, ngfftf, v1scf, xmpi_comm_self)
       else
         if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Could not find: ",ktoa(qpt), "in DVDB - interpolating"))
         ! Fourier interpolate of the potential
         ABI_CHECK(any(abs(qpt) > tol12), "qpt cannot be zero if Fourier interpolation is used")
         cplex = 2
         ABI_MALLOC(v1scf, (cplex,nfftf,nspden,natom3))
         call dvdb_ftinterp_qpt(dvdb, qpt, nfftf, ngfftf, v1scf, xmpi_comm_self)
       end if

       ! Examine the symmetries of the q wavevector
       !call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

       ! TODO: Make sure that symmetries in Q-space are preserved.
       ! Avoid fourq if qpt is in ddb

       ! Get phonon frequencies and displacements in reduced coordinates for this q-point
       call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

       ! Find k+q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qpt

       call listkk(dksqmax,cryst%gmet,indkk_kq,ebands%kptns,kq,ebands%nkpt,1,cryst%nsym,&
         1,cryst%symafm,cryst%symrel,sigma%timrev,use_symrec=.False.)

       if (dksqmax > tol12) then
         write(msg, '(7a,es16.6,4a)' )&
          'The WFK file cannot be used to start thee present calculation ',ch10,&
          'It was asked that the wavefunctions be accurate, but',ch10,&
          'at least one of the k points could not be generated from a symmetrical one.',ch10,&
          'dksqmax=',dksqmax,ch10,&
          'Action: check your WFK file and k point input variables',ch10,&
          '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
         MSG_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1,1); isym_kq = indkk_kq(1,2)
       trev_kq = indkk_kq(1, 6); g0_kq = indkk_kq(1, 3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:,ikq_ibz)

       ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
       ! Be careful with time-reversal symmetry.
       if (isirr_kq) then
         ! Copy u_kq(G)
         !write(std_out,*)"before isirr"
         istwf_kq = wfd%istwfk(ikq_ibz); npw_kq = wfd%npwarr(ikq_ibz)
         ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
         kg_kq(:,1:npw_kq) = wfd%kdata(ikq_ibz)%kg_k
       else
         ! Will Reconstruct u_kq(G) from the IBZ image.
         !istwf_kq = set_istwfk(kq); if (.not. have_ktimerev) istwf_kq = 1
         !call change_istwfk(from_npw,from_kg,from_istwfk,to_npw,to_kg,to_istwfk,n1,n2,n3,ndat1,from_cg,to_cg)
         istwf_kq = 1
         call get_kg(kq,istwf_kq,ecut,cryst%gmet,npw_kq,gtmp)
         ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
         kg_kq(:,1:npw_kq) = gtmp(:,:npw_kq)
         ABI_FREE(gtmp)
       end if

       ! This array will store H1 |psi_nk>
       ABI_MALLOC(h1kets_kq, (2, npw_kq*nspinor, nbcalc_ks))

       ! Allocate vlocal1 with correct cplex. Note nvloc
       ABI_STAT_MALLOC(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)
       ABI_CHECK(ierr==0, "oom vlocal1")

       ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
       do ipc=1,natom3
         call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
           pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))
       end do
       !write(std_out,*)"before load spin"

       ! if PAW, one has to solve a generalized eigenproblem
       ! BE careful here because I will need sij_opt==-1
       gen_eigenpb = (psps%usepaw==1)
       sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm==1) then
       !   optder=0; if (useylmgr==1) optder=1
       !   call initylmg(cryst%gprimd,kg_k,kk,mkmem1,mpi_enreg,psps%mpsang,mpw,nband,mkmem1,&
       !     [npw_k],dtset%nsppol,optder,cryst%rprimd,ylm_k,ylmgr)
       !   call initylmg(cryst%gprimd,kg_kq,kq,mkmem1,mpi_enreg,psps%mpsang,mpw,nband,mkmem1,&
       !     [npw_kq],dtset%nsppol,optder,cryst%rprimd,ylm_kq,ylmgr_kq)
       !end if

       ! Loop over all 3*natom perturbations.
       do ipc=1,natom3
         idir = mod(ipc-1, 3) + 1
         ipert = (ipc - idir) / 3 + 1

         ! Prepare application of the NL part.
         !write(std_out,*)"before init_rf spin"
         call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=1)
         call load_spin_rf_hamiltonian(rf_hamkq,gs_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc), &
           comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
         !write(std_out,*)"after load_rf spin"

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&  ! In
           cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&          ! In
           npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&         ! In
           dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)     ! Out
         !write(std_out,*)"after getgh1c_setup"

         ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0)
         do ib_k=1,nbcalc_ks
           band = ib_k + bstart_ks - 1
           eig0nk = ebands%eig(band,ik_ibz,spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0,copt0,kets_k(:,:,ib_k),cwaveprj0,h1kets_kq(:,:,ib_k),&
             grad_berry,gs1c,gs_hamkq,gvnl1,idir,ipert,eshift,mpi_enreg,optlocal,&
             optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)

           ! TODO: Solve Sternheimer equations non-self-consistently (??)
         end do

         call destroy_rf_hamiltonian(rf_hamkq)

         ABI_FREE(kinpw1)
         ABI_FREE(kpg1_k)
         ABI_FREE(kpg_k)
         ABI_FREE(dkinpw)
         ABI_FREE(ffnlk)
         ABI_FREE(ffnl1)
         ABI_FREE(ph3d)
         if (allocated(ph3d1)) then
           ABI_FREE(ph3d1)
         end if

       end do ! ipc (loop over 3*natom atomic perturbations)

       ABI_FREE(gs1c)

       ! Sum over bands (all bands are taken for the time being)
       istwf_kirr = wfd%istwfk(ikq_ibz); npw_kirr = wfd%npwarr(ikq_ibz)
       ABI_MALLOC(bras_kq, (2, npw_kq*nspinor, 1))
       ABI_MALLOC(cgwork, (2, npw_kirr*nspinor))

       do ib_kq=1,0
       !do ib_kq=1,wfd%nband(ikq_ibz, spin)

         ! symmetrize wavefunctions from IBZ (if needed).
         ! Be careful with time-reversal symmetry.
         if (isirr_kq) then
           ! Copy u_kq(G)
           call wfd_copy_cg(wfd, ib_kq, ikq_ibz, spin, bras_kq(1,1,1))
         else
           ! Reconstruct u_kq(G) from the IBZ image.
           !istwf_kq = set_istwfk(kq); if (.not. have_ktimerev) istwf_kq = 1
           !call change_istwfk(from_npw,from_kg,from_istwfk,to_npw,to_kg,to_istwfk,n1,n2,n3,ndat1,from_cg,to_cg)

           ! Use cgwork as workspace array, results stored in bras_kq
           !g0_kq =  g0ibz_kq + g0bz_kq
           call wfd_copy_cg(wfd, ib_kq, ikq_ibz, spin, cgwork)
           call cg_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1,&
                          npw_kirr, wfd%kdata(ikq_ibz)%kg_k,&
                          npw_kq, kg_kq, istwf_kirr, istwf_kq, cgwork, bras_kq(:,:,1), work_ngfft, work)
         end if

         do ipc=1,natom3
           ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
           !The array eig1_k contains:
           !
           ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
           ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
           do ib_k=1,nbcalc_ks
             call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras_kq(1,1,1),h1kets_kq(1,1,ib_k),&
                  mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
             gkk_atm(:,ib_k,ipc) = [dotr, doti]
           end do
         end do

         ! Get gkk(k, q, nu)
         call gkknu_from_atm(1, nbcalc_ks, 1, natom, gkk_atm, phfrq, displ_red, gkk_nu, num_smallw)

         ! Accumulate contribution to phonon self-energy
         eig0mkq = ebands%eig(ib_kq,ikq_ibz,spin)
         !f_mkq = ebands%occ(ib_kq,ikq_ibz,spin)

         do imode=1,natom3
           wqnu = phfrq(imode)
           do ib_k=1,nbcalc_ks
             eig0nk = ebands%eig(ib_k,ik_ibz,spin)
             !f_nk = ebands%occ(ib_k,ik_ibz,spin)
             gkk2 = gkk_nu(1,ib_k,imode) ** 2 + gkk_nu(2,ib_k,imode) ** 2

             do it=1,sigma%ntemp
               nqnu = nbe(wqnu, sigma%kTmesh(it), zero)
               f_mkq = nfd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))
               f_nk = nfd(eig0nk, sigma%kTmesh(it), sigma%mu_e(it))

               ! Accumulate Sigma(w) for state ib_k
               fact_wr(:) = (nqnu + f_mkq      ) / (sigma%wrmesh_b(:,ib_k) - eig0mkq + wqnu - ieta) + &
                            (nqnu - f_mkq + one) / (sigma%wrmesh_b(:,ib_k) - eig0mkq - wqnu - ieta)
               !sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + gkk2 * fact_wr(:) !* wq

               ! Accumulate d(Sigma)/dw for state ib_k
               !fact_wr(1) = (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu - ieta) + &
               !             (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu - ieta)
               !sigma%dvals_dwr(it, ib_k) = sigma%dvals_dwr(it, ib_k) + gkk2 * fact_wr(1) !* wq
             end do

           end do
         end do ! imode

       end do ! ib_kq (sum over bands at k+q)

       ABI_FREE(bras_kq)
       ABI_FREE(cgwork)
       ABI_FREE(h1kets_kq)
       ABI_FREE(v1scf)
       ABI_FREE(vlocal1)

       call cwtime(cpu,wall,gflops,"stop")
       write(msg,'(a,i0,2(a,f8.2))')"q-point [",iq_ibz,"] completed. cpu:",cpu,", wall:",wall
       call wrtout(std_out,msg,"COLL",do_flush=.True.)
     end do ! iq_ibz

     ABI_FREE(kets_k)
     ABI_FREE(gkk_atm)
     ABI_FREE(gkk_nu)

     ! Collect results and divide by the total number of q-points in the full mesh.
     call xmpi_sum(sigma%vals_wr, comm, ierr); sigma%vals_wr = sigma%vals_wr / sigma%nqbz
     call xmpi_sum(sigma%dvals_dwr, comm, ierr); sigma%dvals_dwr = sigma%dvals_dwr / sigma%nqbz
     ! Writes the results for a single (k-point, spin) to NETCDF file
     call sigmaph_solve(sigma, ikcalc, spin, ebands)
   end do ! spin

   ABI_FREE(qibz_k)
   ABI_FREE(kg_k)
   ABI_FREE(kg_kq)
   ABI_FREE(ylm_k)
   ABI_FREE(ylm_kq)
   ABI_FREE(ylmgr_kq)
 end do !ikcalc

 call wrtout(std_out, "Computation of Sigma_ph completed", "COLL", do_flush=.True.)

 ! Free memory
 ABI_FREE(gvnl1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(fact_wr)

 call destroy_hamiltonian(gs_hamkq)
 call sigmaph_free(sigma)
 call wfd_free(wfd)
 call pawcprj_free(cwaveprj0)
 ABI_DT_FREE(cwaveprj0)

end subroutine sigmaph_driver
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/gkknu_from_atm
!! NAME
!!  gkknu_from_atm
!!
!! FUNCTION
!!  Transform the gkk matrix elements from (atom, red_direction) basis to phonon-mode basis
!!
!! INPUTS
!!  nb1,nb2=Number of bands in gkk_atm matrix.
!!  nk=Number of k-points (usually 1)
!!  natom=Number of atoms.
!!  gkk_atm(2,nb1,nb2,3*natom)=EPH matrix elements in the atomic basis.
!!  phfrq(3*natom)=Phonon frequencies in Ha
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduded coordinates.
!!
!! OUTPUT
!!  gkk_nu(2,nb1,nb2,3*natom)=EPH matrix elements in the phonon-mode basis.
!!  num_smallw=Number of negative/too small frequencies that have been ignored
!!    by setting the corresponding gkk_nu to zero.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkknu_from_atm(nb1, nb2, nk, natom, gkk_atm, phfrq, displ_red, gkk_nu, num_smallw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkknu_from_atm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb1, nb2, nk, natom
 integer,intent(out) :: num_smallw
!arrays
 real(dp),intent(in) :: phfrq(3*natom),displ_red(2,3*natom,3*natom)
 real(dp),intent(in) :: gkk_atm(2,nb1,nb2,nk,3*natom)
 real(dp),intent(out) :: gkk_nu(2,nb1,nb2,nk,3*natom)

!Local variables-------------------------
!scalars
 integer :: nu,ipc

! *************************************************************************

 gkk_nu = zero; num_smallw = 0

 ! Loop over phonon branches.
 do nu=1,3*natom
   ! Ignore negative or too small frequencies
   if (phfrq(nu) < tol6) then
     num_smallw = num_smallw + 1; cycle
   end if

   ! Transform the gkk from atom,red_dir basis to phonon mode basis
   do ipc=1,3*natom
     gkk_nu(1,:,:,:,nu) = gkk_nu(1,:,:,:,nu) &
       + gkk_atm(1,:,:,:,ipc) * displ_red(1,ipc,nu) &
       - gkk_atm(2,:,:,:,ipc) * displ_red(2,ipc,nu)
     gkk_nu(2,:,:,:,nu) = gkk_nu(2,:,:,:,nu) &
       + gkk_atm(1,:,:,:,ipc) * displ_red(2,ipc,nu) &
       + gkk_atm(2,:,:,:,ipc) * displ_red(1,ipc,nu)
   end do

   gkk_nu(:,:,:,:,nu) = gkk_nu(:,:,:,:,nu) / sqrt(two * phfrq(nu))
 end do

end subroutine gkknu_from_atm
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/nfd
!! NAME
!!  nfd
!!
!! FUNCTION
!!  Fermi-Dirac statistic
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   kT=Value of K_Boltzmann x T in Ha.
!!   mu=Chemical potential in Ha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function nfd(ee, kT, mu)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nfd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

! *************************************************************************

 !TODO: Find decent value.
 if (kT > tol16) then
   nfd = one / (exp((ee-mu) / kT) + one)
 else
   ! Heaviside
   if (ee > mu) then
     nfd = zero
   else if (ee < mu) then
     nfd = one
   else
     nfd = half
   end if
 end if

end function nfd
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/nbe
!! NAME
!!  nbe
!!
!! FUNCTION
!!   Bose-Einstein statistic
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   kT=Value of K_Boltzmann x T in Ha.
!!   mu=Chemical potential in Ha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function nbe(ee, kT, mu)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nbe'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

! *************************************************************************

 if (kT > tol16) then
   nbe = one / (exp((ee-mu) / kT) - one)
 else
   ! No condensate for T-->0
   nbe = zero
 end if

end function nbe
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_new
!! NAME
!!  sigmaph_new
!!
!! FUNCTION
!!  Creation method (allocates memory, initialize data from input vars).
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst(crystal_t)=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type (sigmaph_t) function sigmaph_new(dtset, ecut, cryst, ebands, dtfil, comm) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_new'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 real(dp),intent(in) :: ecut
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(datafiles_type),intent(in) :: dtfil

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0,brav1=1,option1=1,occopt3=3
 integer :: my_rank,nqpt_max,ik,my_nshiftq,nct,my_mpw,cnt,nproc,iq_ibz
 integer :: onpw,ii,ipw,ierr,it
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp),parameter :: spinmagntarget=-99.99_dp
 real(dp) :: tstep,dksqmax
 character(len=500) :: msg
 type(ebands_t) :: tmp_ebands
!arrays
 integer :: qptrlatt(3,3),indkk_k(1,6),gmax(3),my_gmax(3)
 integer,allocatable :: gtmp(:,:)
 real(dp) :: my_shiftq(3,1),temp_range(2),kk(3),kq(3)
 real(dp),allocatable :: qibz_k(:,:),tmp_qbz(:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Copy important dimensions.
 new%nsppol = ebands%nsppol; new%nspinor = ebands%nspinor
 new%timrev = 1

 ! Build (linear) mesh of temperatures.
 ! TODO: Will become input eph_tsmesh start stop num
 new%ntemp = 10
 ABI_CHECK(new%ntemp > 1, "ntemp cannot be 1")
 temp_range = [0.6_dp, 1.2_dp]
 tstep = (temp_range(2) - temp_range(1)) / (new%ntemp - 1)
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 new%kTmesh = arth(temp_range(1), tstep, new%ntemp)

 ! Compute the chemical potential at the different physical temperatures with Fermi-Dirac.
 ABI_MALLOC(new%mu_e, (new%ntemp))
 new%mu_e = ebands%fermie

 call ebands_copy(ebands, tmp_ebands)
 do it=1,new%ntemp
   call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(it), spinmagntarget, dtset%prtvol)
   new%mu_e(it) = tmp_ebands%fermie
 end do
 call ebands_free(tmp_ebands)

 ! Frequency mesh for spectral function.
 new%nwr = 11
 ABI_CHECK(mod(new%nwr, 2) == 1, "nwr should be odd!")
 new%wr_step = 0.01

 ! Define q-mesh for integration.
 ! Either q-mesh from DDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation if q not in DDB)
 new%ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   new%ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if

 qptrlatt = 0
 qptrlatt(1,1) = new%ngqpt(1)
 qptrlatt(2,2) = new%ngqpt(2)
 qptrlatt(3,3) = new%ngqpt(3)
 nqpt_max = (product(new%ngqpt) * my_nshiftq)

 ABI_MALLOC(tmp_qbz, (3,nqpt_max))
 call smpbz(brav1,std_out,qptrlatt,nqpt_max,new%nqbz,my_nshiftq,option1,my_shiftq,tmp_qbz)
 ABI_MALLOC(new%qbz, (3, new%nqbz))
 new%qbz = tmp_qbz(:,:new%nqbz)
 ABI_FREE(tmp_qbz)

 ! These parameters will be passed via the input file (similar to kptgw, bdgw)
 ! k-point and bands where corrections are wanted
 ! if symsigma == 1, we have to include all degenerate states in the set
 ! because the final QP corrections will be obtained by averaging the results in the degenerate subspace.
 ! k-point and bands where corrections are wanted
 ! These parameters will be passed via the input file (similar to kptgw, bdgw)
 new%symsigma = 0
 new%nkcalc = ebands%nkpt
 !new%nkcalc = 1

 ABI_MALLOC(new%kcalc, (3, new%nkcalc))
 ABI_MALLOC(new%bstart_ks, (new%nkcalc, new%nsppol))
 ABI_MALLOC(new%nbcalc_ks, (new%nkcalc, new%nsppol))

 new%kcalc = ebands%kptns
 !new%kcalc = zero
 new%bstart_ks = 1; new%nbcalc_ks = 2
 new%max_nbcalc = maxval(new%nbcalc_ks)

 ! The k-point and the symmetries relating the BZ point to the IBZ.
 ABI_MALLOC(new%kcalc2ibz, (new%nkcalc, 6))

 do ik=1,new%nkcalc
   kk = new%kcalc(:,ik)
   call listkk(dksqmax,cryst%gmet,indkk_k,ebands%kptns,kk,ebands%nkpt,1,cryst%nsym,&
      1,cryst%symafm,cryst%symrel,new%timrev,use_symrec=.False.)

   new%kcalc2ibz(ik, :) = indkk_k(1, :)
   if (dksqmax > tol12) then
     write(msg, '(7a,es16.6,4a)' )&
      'The WFK file cannot be used to start thee present calculation ',ch10,&
      'It was asked that the wavefunctions be accurate, but',ch10,&
      'at least one of the k points could not be generated from a symmetrical one.',ch10,&
      'dksqmax=',dksqmax,ch10,&
      'Action: check your WFK file and k point input variables',ch10,&
      '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
      MSG_ERROR(msg)
   end if
 end do

 !call wrtout(std_out, sjoin("Will compute", itoa(nqibz_k), "q-points in the IBZ(k)"))

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! used to symmetrize the wavefunctions in G-space.
 ! TODO: Should loop over IBZ(k)
 new%mpw = 0; new%gmax = 0; cnt = 0
 do ik=1,new%nkcalc
   kk = new%kcalc(:, ik)
   do iq_ibz=1,new%nqbz
   !do iq_ibz=1,nqibz_k
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     !kq = kk + qibz_k(:,iq_ibz)
     kq = kk + new%qbz(:,iq_ibz)

     ! TODO: g0 umklapp here can enter into play!
     ! fstab should contains the max of the umlapp G-vectors.
     ! gmax could not be large enough!
     call get_kg(kq,1,ecut,cryst%gmet,onpw,gtmp)
     new%mpw = max(new%mpw, onpw)
     do ipw=1,onpw
       do ii=1,3
        new%gmax(ii) = max(new%gmax(ii), abs(gtmp(ii,ipw)))
       end do
     end do
     ABI_FREE(gtmp)
   end do
 end do

 my_mpw = new%mpw; call xmpi_max(my_mpw, new%mpw, comm, ierr)
 my_gmax = new%gmax; call xmpi_max(my_gmax, new%gmax, comm, ierr)
 call wrtout(std_out, sjoin('optimal value of mpw=', itoa(new%mpw)))

 ! Allocate arrays used to store final results and set them to zero.
 ABI_CALLOC(new%vals_wr, (new%nwr, new%ntemp, new%max_nbcalc))
 ABI_CALLOC(new%dvals_dwr, (new%ntemp, new%max_nbcalc))
 ABI_CALLOC(new%wrmesh_b, (new%nwr, new%max_nbcalc))

 ! Open netcdf file within MPI communicator comm.
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   ! Master creates the netcdf file used to store the results of the calculation.
   NCF_CHECK(nctk_open_create(new%ncid, strcat(dtfil%filnam_ds(4), "_SIGMAPH.nc"), comm))
   ncid = new%ncid

   NCF_CHECK(crystal_ncwrite(cryst, ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))

   ! Add sigma_ph dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nkcalc", new%nkcalc), &
     nctkdim_t("max_nbcalc", new%max_nbcalc), &
     nctkdim_t("nsppol", new%nsppol), &
     nctkdim_t("nwr", new%nwr), &
     nctkdim_t("ntemp", new%ntemp)], &
     !nctkdim_t("nqbz", new%nqbz)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "symsigma"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "wr_step"])
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("bstart_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("nbcalc_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("kcalc", "int", "three, nkcalc"), &
     nctkarr_t("kTmesh", "dp", "ntemp"), &
     nctkarr_t("mu_e", "dp", "ntemp"), &
     ! These arrays acquire two extra dimensions on file (nkcalc, nsppol).
     nctkarr_t("vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("dvals_dwr", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol") &
   ])
   NCF_CHECK(ncerr)

   if (new%nwr > 1) then
      ! Make room for the spectral function on file.
      ncerr = nctk_def_arrays(ncid, [ &
         nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol")])
      NCF_CHECK(ncerr)
   end if

   ! Write data that do not depend on the (kpt, spin) loop.
   NCF_CHECK(nctk_set_datamode(ncid))
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: "symsigma"], [new%symsigma])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "wr_step"], [new%wr_step])
   NCF_CHECK(ncerr)

   NCF_CHECK(nf90_put_var(ncid, vid("bstart_ks"), new%bstart_ks))
   NCF_CHECK(nf90_put_var(ncid, vid("nbcalc_ks"), new%nbcalc_ks))
   NCF_CHECK(nf90_put_var(ncid, vid("kcalc"), new%kcalc))
   NCF_CHECK(nf90_close(new%ncid))
 end if ! master

 call xmpi_barrier(comm)
 ! Now reopen the file in parallel.
 NCF_CHECK(nctk_open_modify(new%ncid, strcat(dtfil%filnam_ds(4), "_SIGMAPH.nc"), comm))
 NCF_CHECK(nctk_set_datamode(new%ncid))

contains

 integer function vid(vname)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid
!!***
#endif

end function sigmaph_new
!!***

!!****f* m_sigmaph/sigmaph_free
!! NAME
!!  sigmaph_free
!!
!! FUNCTION
!!  Deallocate dynamic memory
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_free(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(sigmaph_t),intent(inout) :: self

! *************************************************************************

 ! integer
 if (allocated(self%bstart_ks)) then
   ABI_FREE(self%bstart_ks)
 end if
 if (allocated(self%nbcalc_ks)) then
   ABI_FREE(self%nbcalc_ks)
 end if
 if (allocated(self%kcalc2ibz)) then
   ABI_FREE(self%kcalc2ibz)
 end if

 ! real
 if (allocated(self%kcalc)) then
   ABI_FREE(self%kcalc)
 end if
 if (allocated(self%kTmesh)) then
   ABI_FREE(self%kTmesh)
 end if
 if (allocated(self%mu_e)) then
   ABI_FREE(self%mu_e)
 end if
 if (allocated(self%wrmesh_b)) then
   ABI_FREE(self%wrmesh_b)
 end if
 if (allocated(self%qbz)) then
   ABI_FREE(self%qbz)
 end if

 ! complex
 if (allocated(self%vals_wr)) then
   ABI_FREE(self%vals_wr)
 end if
 if (allocated(self%dvals_dwr)) then
   ABI_FREE(self%dvals_dwr)
 end if
 !if (allocated(self%vals_qwr)) then
 !  ABI_FREE(self%vals_qwr)
 !end if

 ! Close files.
#ifdef HAVE_NETCDF
 NCF_CHECK(nf90_close(self%ncid))
#endif

end subroutine sigmaph_free
!!***

!!****f* m_sigmaph/sigmaph_solve
!! NAME
!!  sigmaph_solve
!!
!! FUNCTION
!!  Compute QP energies, Z factor and spectral function (if required).
!!  Save results to file.
!!
!! INPUTS
!!  ikcalc=Index of the computed k-point
!!  spin=Spin index.
!!  ebands<ebands_t>=KS band energies.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_solve(self, ikcalc, spin, ebands)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nbe'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, spin
 type(sigmaph_t),target,intent(in) :: self
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!arrays
 integer :: shape3(3),shape4(4)
 real(dp), ABI_CONTIGUOUS pointer :: rdata3(:,:,:),rdata4(:,:,:,:)
 !read(dp),allocatable :: spfunc(:,:,:,:)

! *************************************************************************

 ! Compute QP corrections.

 ! Compute spectral functions.
 if (self%nwr > 1) then
 end if

#ifdef HAVE_NETCDF
 shape3(1) = 2; shape4(1) = 2

 ! Write data (use iso_c_binding to associate a real pointer to complex data because
 ! netcdf does not support complex types.
 shape4(2:) = shape(self%vals_wr)
 call c_f_pointer(c_loc(self%vals_wr), rdata4, shape4)
 NCF_CHECK(nf90_put_var(self%ncid, vid("vals_wr"), rdata4, start=[1,1,1,1,ikcalc,spin]))

 shape3(2:) = shape(self%dvals_dwr)
 call c_f_pointer(c_loc(self%dvals_dwr), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, vid("dvals_dwr"), rdata3, start=[1,1,1,ikcalc,spin]))

 !nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol")])
 !shape4(2:) = shape(self%spfunc_wr)
 !call c_f_pointer(c_loc(self%spfunc_wr), rdata4, shape4)
 !NCF_CHECK(nf90_put_var(self%ncid, vid("spfunc_dwr"), rdata4, start=[1,1,ikcalc,spin]))

 !NCF_CHECK(nf90_put_var(self%ncid, vid("vals_qwr"), self%vals_qwr))

contains
 integer function vid(vname)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(self%ncid, vname)
 end function vid
#endif

end subroutine sigmaph_solve
!!***

end module m_sigmaph
!!***

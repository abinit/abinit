!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigmaph
!! NAME
!!  m_sigmaph
!!
!! FUNCTION
!!  Compute the matrix elements of the Fan-Migdal Debye-Waller self-energy in the KS basis set.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
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
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
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
 use m_time,           only : cwtime, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit
 use m_io_tools,       only : iomode_from_fname
 use m_special_funcs,  only : dirac_delta, gspline_t, gspline_new, gspline_eval, gspline_free
 use m_fftcore,        only : ngfft_seq
 use m_cgtools,        only : dotprod_g !set_istwfk
 use m_cgtk,           only : cgtk_rotate
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, listkk
 use m_lgroup,         only : lgroup_t, lgroup_new, lgroup_free
 use m_fftcore,        only : get_kg
 use m_kg,             only : getph
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_fourier_interpol, only : transgrid
 use m_symkpt,     only : symkpt
! use m_paw_an,	       only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
! use m_paw_ij,	       only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
! use m_pawfgrtab,     only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
! use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
! use m_pawdij,	       only : pawdij, symdij

 implicit none

 private
!!***

 public :: sigmaph   ! Main entry point to compute self-energy matrix elements
!!***

 ! Tables for degenerated KS states.
 type bids_t
   integer, allocatable :: vals(:)
 end type bids_t

 type degtab_t
   type(bids_t), allocatable :: bids(:)
 end type degtab_t

!----------------------------------------------------------------------

!!****t* m_sigmaph/sigmaph_t
!! NAME
!! sigmaph_t
!!
!! FUNCTION
!! Container for the (diagonal) matrix elements of the electronic self-energy (phonon contribution)
!! computed in the KS representation i.e. Sigma_eph(omega, T, band, k, spin).
!! Provides methods to compute the QP corrections, the spectral functions and save the results to file.
!!
!! TODO
!!  1) Store (q, mu) or q-resolved contributions? More memory but could be useful for post-processing
!!  2) Use list of i.delta values instead of single shift? Useful for convergence studies.
!!
!! SOURCE

 type,private :: sigmaph_t

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
   ! Required when the sum over q in the BZ is replaced by IBZ(k).

  integer :: timrev
   !  timrev=1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nbsum
   ! Number of bands used in sum over states.

  integer :: nqbz
  ! Number of q-points in the BZ.

  integer :: nqibz
  ! Number of q-points in the IBZ.

  integer :: nqibz_k
  ! Number of q-points in the IBZ(k). Depends on ikcalc.

  integer :: ncid = nctk_noid
  ! Netcdf file used to save results.

  integer :: mpw
  ! Maximum number of PWs for all possible k+q

  complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  integer :: gmax(3)

  integer :: ngqpt(3)
   ! Number of divisions in the Q mesh in the BZ.

  integer,allocatable :: bstart_ks(:,:)
   ! bstart_ks(nkcalc,nsppol)
   ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used.

  integer,allocatable :: nbcalc_ks(:,:)
   ! nbcalc_ks(nkcalc,nsppol)
   ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used.

  integer,allocatable :: kcalc2ibz(:,:)
    !kcalc2ibz(nkcalc, 6))
    ! Mapping kcalc --> ibz as reported by listkk.

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.

  real(dp),allocatable :: qbz(:,:)
  ! qbz(3,nqbz)
  ! Reduced coordinates of the q-points in the full BZ.

  real(dp),allocatable :: qibz(:,:)
  ! qibz(3,nqibz)
  ! Reduced coordinates of the q-points in the IBZ.

  real(dp),allocatable :: wtq(:)
  ! wtq(nqibz)
  ! Weights of the q-points in the IBZ (normalized to one).

  real(dp),allocatable :: qibz_k(:,:)
  ! qibz(3,nqibz_k)
  ! Reduced coordinates of the q-points in the IBZ(k). Depends on ikcalc.

  real(dp),allocatable :: wtq_k(:)
  ! wtq(nqibz_k)
  ! Weights of the q-points in the IBZ(k) (normalized to one). Depends on ikcalc.

  real(dp),allocatable :: kTmesh(:)
   ! kTmesh(ntemp)
   ! List of temperatures (kT units).

  real(dp),allocatable :: mu_e(:)
   ! mu_e(ntemp)
   ! chemical potential of electrons for the different temperatures.

  real(dp),allocatable :: wrmesh_b(:,:)
  ! wrmesh_b(nwr, max_nbcalc)
  ! Frequency mesh along the real axis (Ha units) used for the different bands
  ! This array depends on (ikcalc, spin)

  complex(dpc),allocatable :: vals_e0ks(:,:)
   ! vals_e0ks(ntemp, max_nbcalc)
   ! Sigma_eph(omega=eKS, kT, band) for fixed (kcalc, spin).

  complex(dpc),allocatable :: dvals_de0ks(:,:)
   ! dvals_de0ks(ntemp, max_nbcalc) for fixed (kcalc, spin)
   ! d Sigma_eph(omega, kT, band, kcalc, spin) / d omega (omega=eKS)

  real(dp),allocatable :: dw_vals(:,:)
  !  dw_vals(ntemp, max_nbcalc) for fixed (kcalc, spin)
  !  Debye-Waller term (static).

  !real(dp),allocatable :: qpadb_enes(:,:)
  ! qp_adbenes(ntemp, max_nbcalc)
  ! (Real) QP energies computed with the adiabatic formalism.

  !complex(dpc),allocatable :: qp_enes(:,:)
  ! qp_enes(ntemp, max_nbcalc)
  ! (Complex) QP energies computed with the non-adiabatic formalism.

  complex(dpc),allocatable :: vals_wr(:,:,:)
   ! vals_wr(nwr, ntemp, max_nbcalc)
   ! Sigma_eph(omega, kT, band) for given (k, spin).
   ! enk_KS corresponds to nwr/2 + 1.
   ! This array depends on (ikcalc, spin)

  integer :: gfw_nomega

  real(dp),allocatable :: gfw_mesh(:)
   ! gfw_mesh(gfw_nomega)

  complex(dpc),allocatable :: gfw_vals(:,:,:,:)
   !gfw_vals(gfw_nomega, ntemp, 2, max_nbcalc)

  logical :: has_nuq_terms

  complex(dpc),allocatable :: vals_nuq(:,:,:,:,:)
   ! vals(ntemp, max_nbcalc, natom*3, nqbz, 2)
   ! (nu, q)-resolved self-energy for given (k, spin)
   ! First slice: Fan term at omega=e0_KS
   ! Second slice: DW term

  type(degtab_t),allocatable :: degtab(:,:)
   ! (nkcalc, nsppol)
   ! Table used to average QP results in the degenerate subspace if symsigma == 1

 end type sigmaph_t
!!***

 private :: sigmaph_new               ! Creation method (allocates memory, initialize data from input vars).
 private :: sigmaph_free              ! Free memory.
 private :: sigmaph_setup_kcalc       ! Return tables used to perform the sum over q-points for given k-point.
 private :: sigmaph_gather_and_write  ! Compute the QP corrections.
 private :: sigmaph_print             ! Print results to main output file.

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph
!! NAME
!!  sigmaph
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

subroutine sigmaph(wfk0_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
                   pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph'
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
 integer,parameter :: useylmgr=0,useylmgr1=0,master=0,ndat1=1
 integer :: my_rank,mband,my_minb,my_maxb,nsppol,nkpt,iq_ibz
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,nprocs
 integer :: ibsum_kq,ib_k,band,num_smallw,ibsum,ii,im,in,ndeg !,ib,nstates
 integer :: idir,ipert,ip1,ip2,idir1,ipert1,idir2,ipert2
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq !,!timerev_q,
 integer :: spin,istwf_k,istwf_kq,istwf_kqirr,npw_k,npw_kq,npw_kqirr
 integer :: mpw,ierr,it !ipw
 integer :: n1,n2,n3,n4,n5,n6,nspden,do_ftv1q,nu
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnl1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,nq
 integer :: nbcalc_ks,nbsum,bstart_ks,ikcalc
 real(dp),parameter :: tol_enediff=0.001_dp*eV_Ha
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,gflops_all
 real(dp) :: ecut,eshift,dotr,doti,dksqmax,weigth_q,rfact,alpha,beta,gmod2,hmod2
 complex(dpc) :: cfact,dka,dkap,dkpa,dkpap,my_ieta,cplx_ediff
 logical,parameter :: have_ktimerev=.True.
 logical :: isirr_k,isirr_kq,gen_eigenpb,isqzero
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(sigmaph_t) :: sigma
 type(gspline_t) :: gspl
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0_kq(3),dummy_gvec(3,dummy_npw)
 integer :: work_ngfft(18),gmax(3) !!g0ibz_kq(3),
 integer :: indkk_kq(1,6)
 integer,allocatable :: gtmp(:,:),kg_k(:,:),kg_kq(:,:),nband(:,:),distrib_bq(:,:),deg_ibk(:) !,degblock(:,:),
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),phfrq(3*cryst%natom),sqrt_phfrq0(3*cryst%natom)
 real(dp) :: lf(2),rg(2),res(2)
 real(dp) :: wqnu,nqnu,gkk2,eig0nk,eig0mk,eig0mkq,f_mkq
 !real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),phfrq(3*cryst%natom)
 !real(dp) :: wqnu,nqnu,gkk2,eig0nk,eig0mk,eig0mkq,ediff,f_mkq !,f_nk
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom),displ_red(2,3,cryst%natom,3*cryst%natom)
 !real(dp) :: ucart(2,3,cryst%natom,3*cryst%natom)
 real(dp) :: d0mat(2,3*cryst%natom,3*cryst%natom)
 complex(dpc) :: cmat(3*cryst%natom,3*cryst%natom),cvec1(3*cryst%natom),cvec2(3*cryst%natom)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkk_atm(:,:,:),gkk_nu(:,:,:),dbwl_nu(:,:,:,:),gdw2_mn(:,:),gkk0_atm(:,:,:,:)
 complex(dpc),allocatable :: tpp(:,:),hka_mn(:,:,:),wmat1(:,:,:)
 real(dp),allocatable :: bra_kq(:,:),kets_k(:,:,:),h1kets_kq(:,:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnl1(:,:),work(:,:,:,:)
 real(dp),allocatable ::  gs1c(:,:),nqnu_tlist(:),dt_weights(:,:)
 complex(dpc),allocatable :: cfact_wr(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)

!************************************************************************

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor
 nspden = dtset%nspden; nkpt = ebands%nkpt

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx,natom,n1,n2,n3,ph1d,cryst%xred)

 ! Construct object to store final results.
 ecut = dtset%ecut ! dtset%dilatmx
 sigma = sigmaph_new(dtset, ecut, cryst, ebands, ifc, dtfil, comm)
 if (my_rank == master) call sigmaph_print(sigma, dtset, ab_out)

 ! This is the maximum number of PWs for all possible k+q treated.
 mpw = sigma%mpw; gmax = sigma%gmax

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))

 ! Initialize the wave function descriptor (read up to mband states where mband is defined by dtset%nband)
 ! For the time being, no memory distribution, each node has the full set of states.
 ! TODO: Distribute memory
 mband = dtset%mband
 my_minb = 1; my_maxb = mband
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask,(mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur,(mband, nkpt ,nsppol))
 nband=mband; bks_mask=.True.; keep_ur=.False.

 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,dtset%paral_kgb,dummy_npw,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,dtset%ecutsm,dtset%dilatmx,ebands%istwfk,ebands%kptns,ngfft,&
   dummy_gvec,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm,opt_ecut=ecut)

 call wfd_print(wfd,header="Wavefunctions for self-energy calculation.",mode_paral='PERS')
 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)

 call wfd_read_wfk(wfd, wfk0_path, iomode_from_fname(wfk0_path))
 if (.False.) call wfd_test_ortho(wfd, cryst, pawtab, unit=std_out, mode_paral="PERS")

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

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nsppol,nspden,natom,&
&  dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
&  comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
&  usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

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
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 vlocal = huge(one)
 if (sigma%nwr > 0) then
   ABI_MALLOC(cfact_wr, (sigma%nwr))
 end if
 ABI_MALLOC(nqnu_tlist, (sigma%ntemp))

 ! Open the DVDB file
 call dvdb_open_read(dvdb, ngfftf, xmpi_comm_self)
 if (dtset%prtvol > 10) dvdb%debug = .True.
 ! This to symmetrize the DFPT potentials.
 if (dtset%symdynmat == 1) dvdb%symv1 = .True.
 call dvdb_print(dvdb, prtvol=dtset%prtvol)

 ! Compute gaussian spline.
 if (sigma%gfw_nomega > 0) then
   gspl = gspline_new(three * (sigma%gfw_mesh(2) - sigma%gfw_mesh(1)))
   ABI_MALLOC(dt_weights, (sigma%gfw_nomega, 2))
 end if

 ! Loop over k-points in Sigma_nk.
 do ikcalc=1,sigma%nkcalc
   kk = sigma%kcalc(:, ikcalc)

   ! Find IBZ(k) for q-point integration.
   call sigmaph_setup_kcalc(sigma, cryst, ikcalc)
   call wrtout(std_out, sjoin("Computing self-energy matrix elements for k-point:", ktoa(kk)))
   call wrtout(std_out, sjoin("Number of q-points in the IBZ(k):", itoa(sigma%nqibz_k)))

   ! Symmetry indices for kk.
   ik_ibz = sigma%kcalc2ibz(ikcalc,1); isym_k = sigma%kcalc2ibz(ikcalc,2)
   trev_k = sigma%kcalc2ibz(ikcalc,6); g0_k = sigma%kcalc2ibz(ikcalc,3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   ABI_CHECK(isirr_k, "For the time being the k-point must be in the IBZ")
   kk_ibz = ebands%kptns(:,ik_ibz)
   npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)

   ! Activate Fourier interpolation if q-points are not in the DVDB file.
   ! TODO: handle q_bz = S q_ibz case by symmetrizing the potentials already available in the DVDB.
   ! without performing FT interpolation.
   do_ftv1q = 0
   do iq_ibz=1,sigma%nqibz_k
     if (dvdb_findq(dvdb, sigma%qibz_k(:,iq_ibz)) == -1) do_ftv1q = do_ftv1q + 1
   end do
   if (do_ftv1q /= 0) then
     write(msg, "(2(a,i0),a)")"Will use Fourier interpolation of DFPT potentials [",do_ftv1q,"/",sigma%nqibz_k,"]"
     call wrtout(std_out, msg)
     call wrtout(std_out, sjoin("From ngqpt", ltoa(ifc%ngqpt), "to", ltoa(sigma%ngqpt)))
     call dvdb_ftinterp_setup(dvdb, ifc%ngqpt, 1, [zero,zero,zero], nfftf, ngfftf, comm)
   end if

   ! Allocate PW-arrays. Note mpw in kg_kq
   ABI_MALLOC(kg_k, (3, npw_k))
   kg_k = wfd%kdata(ik_ibz)%kg_k
   ABI_MALLOC(kg_kq, (3, mpw))

   ! Spherical Harmonics for useylm==1.
   ABI_MALLOC(ylm_k,(mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylm_kq,(mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

   do spin=1,nsppol
     ! Bands in Sigma_nk to compute and number of bands in sum over states.
     bstart_ks = sigma%bstart_ks(ikcalc, spin)
     nbcalc_ks = sigma%nbcalc_ks(ikcalc, spin)
     nbsum = sigma%nbsum

     ! Zero self-energy matrix elements. Build frequency mesh for nk states.
     sigma%vals_e0ks = zero; sigma%dvals_de0ks = zero; sigma%dw_vals = zero

     if (sigma%nwr > 0) then
       sigma%vals_wr = zero
       do ib_k=1,nbcalc_ks
         band = ib_k + bstart_ks - 1
         eig0nk = ebands%eig(band,ik_ibz,spin) - sigma%wr_step * (sigma%nwr / 2)
         sigma%wrmesh_b(:,ib_k) = arth(eig0nk, sigma%wr_step, sigma%nwr)
       end do
     end if
     if (sigma%gfw_nomega > 0) sigma%gfw_vals = zero
     if (sigma%has_nuq_terms) sigma%vals_nuq = zero

     ! Allocate eph matrix elements.
     ABI_MALLOC(gkk_atm, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkk_nu, (2, nbcalc_ks, natom3))
     ABI_STAT_MALLOC(dbwl_nu, (2, nbcalc_ks, nbsum, natom3), ierr)
     ABI_CHECK(ierr == 0, "oom in dbwl_nu")
     dbwl_nu = zero
     ABI_MALLOC(gkk0_atm, (2, nbcalc_ks, nbsum, natom3))
     ABI_CHECK(ierr == 0, "oom in gkk0_atm")
     gkk0_atm = zero

     ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
     ! TODO: symmetrize them if kk is not irred
     ABI_MALLOC(kets_k, (2, npw_k*nspinor, nbcalc_ks))
     do ib_k=1,nbcalc_ks
       band = ib_k + bstart_ks - 1
       call wfd_copy_cg(wfd, band, ik_ibz, spin, kets_k(1,1,ib_k))
     end do

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Distribute q-points and bands.
     ABI_MALLOC(distrib_bq, (nbsum, sigma%nqibz_k))
     if (sigma%nqibz_k >= nprocs .and. mod(sigma%nqibz_k, nprocs) == 0) then
       do iq_ibz=1,sigma%nqibz_k
         distrib_bq(:, iq_ibz) = mod(iq_ibz, nprocs)
       end do
     else
       do it=1,nbsum*sigma%nqibz_k
         ibsum_kq = mod(it-1, nbsum) + 1; iq_ibz = (it - ibsum_kq) / nbsum + 1
         distrib_bq(ibsum_kq, iq_ibz) = mod(it, nprocs)
       end do
     end if

     ! Integrations over q-points in the IBZ(k)
     do iq_ibz=1,sigma%nqibz_k
       ! Quick-parallelization over q-points
       if (all(distrib_bq(1:nbsum, iq_ibz) /= my_rank)) cycle

       qpt = sigma%qibz_k(:,iq_ibz)
       isqzero = (sum(qpt**2) < tol14) !; if (isqzero) cycle
       call cwtime(cpu,wall,gflops,"start")

       ! Find the index of the q-point in the DVDB.
       db_iqpt = dvdb_findq(dvdb, qpt)

       ! TODO: handle q_bz = S q_ibz case by symmetrizing the potentials already available in the DVDB.
       if (db_iqpt /= -1) then
         if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found:", ktoa(qpt), "in DVDB with index", itoa(db_iqpt)))
         ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
         ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
         call dvdb_readsym_allv1(dvdb, db_iqpt, cplex, nfftf, ngfftf, v1scf, xmpi_comm_self)
       else
         if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Could not find:", ktoa(qpt), "in DVDB - interpolating"))
         ! Fourier interpolation of the potential
         ABI_CHECK(any(abs(qpt) > tol12), "qpt cannot be zero if Fourier interpolation is used")
         cplex = 2
         ABI_MALLOC(v1scf, (cplex,nfftf,nspden,natom3))
         call dvdb_ftinterp_qpt(dvdb, qpt, nfftf, ngfftf, v1scf, xmpi_comm_self)
         !v1scf = zero
       end if

       ! TODO: Make sure that symmetries in Q-space are preserved.
       ! Avoid fourq if qpt is in ddb
       ! Examine the symmetries of the q wavevector
       !call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

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
         write(msg, '(3a,es16.6,7a)' )&
          "The WFK file cannot be used to compute self-energy corrections.",ch10,&
          "At least one of the k-points could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
          "Q-mesh: ",ltoa(sigma%ngqpt),", K-mesh (from kptrlatt) ",ltoa(get_diag(dtset%kptrlatt)), ch10, &
          'Action: check your WFK file and (k,q) point input variables'
         MSG_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1,1); isym_kq = indkk_kq(1,2)
       trev_kq = indkk_kq(1, 6); g0_kq = indkk_kq(1, 3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0)) !; isirr_kq = .True.
       kq_ibz = ebands%kptns(:,ikq_ibz)

       ! Get npw_kq, kg_kq for k+q
       ! Be careful with time-reversal symmetry and istwf_kq
       if (isirr_kq) then
         ! Copy u_kq(G)
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

       ! Allocate array to store H1 |psi_nk> for all 3*natom perturbations
       ABI_STAT_MALLOC(h1kets_kq, (2, npw_kq*nspinor, nbcalc_ks, natom3), ierr)
       ABI_CHECK(ierr==0, "oom in h1kets_kq")

       ! Allocate vlocal1 with correct cplex. Note nvloc
       ABI_STAT_MALLOC(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)
       ABI_CHECK(ierr==0, "oom in vlocal1")

       ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
       do ipc=1,natom3
         call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
           pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))
       end do
       !vlocal1 = one

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
       ! In the inner loop, I calculate H1 * psi_k, stored in h1kets_kq on the k+q sphere.
       do ipc=1,natom3
         idir = mod(ipc-1, 3) + 1; ipert = (ipc - idir) / 3 + 1

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=.true.)
             !&paw_ij1=paw_ij1,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
             !&mpi_spintab=mpi_enreg%my_isppoltab)
         call load_spin_rf_hamiltonian(rf_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&  ! In
           cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&          ! In
           npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&         ! In
           dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)     ! Out

         ! Compute H(1) applied to GS wavefunction Psi_nk(0)
         do ib_k=1,nbcalc_ks
           band = ib_k + bstart_ks - 1
           eig0nk = ebands%eig(band,ik_ibz,spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0,kets_k(:,:,ib_k),cwaveprj0,h1kets_kq(:,:,ib_k,ipc),&
             grad_berry,gs1c,gs_hamkq,gvnl1,idir,ipert,eshift,mpi_enreg,optlocal,&
             optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
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

       ! ==============
       ! Sum over bands
       ! ==============
       istwf_kqirr = wfd%istwfk(ikq_ibz); npw_kqirr = wfd%npwarr(ikq_ibz)
       ABI_MALLOC(bra_kq, (2, npw_kq*nspinor))
       ABI_MALLOC(cgwork, (2, npw_kqirr*nspinor))

       do ibsum_kq=1,nbsum
         if (distrib_bq(ibsum_kq, iq_ibz) /= my_rank) cycle
         ! This is to check whether the gkk elements in the degenerate subspace break symmetry
         !if (ibsum_kq >= bstart_ks .and. ibsum_kq <= bstart_ks + nbcalc_ks - 1) cycle

         ! symmetrize wavefunctions from IBZ (if needed).
         ! Be careful with time-reversal symmetry.
         if (isirr_kq) then
           ! Copy u_kq(G)
           call wfd_copy_cg(wfd, ibsum_kq, ikq_ibz, spin, bra_kq)
         else
           ! Reconstruct u_kq(G) from the IBZ image.
           !istwf_kq = set_istwfk(kq); if (.not. have_ktimerev) istwf_kq = 1
           !call change_istwfk(from_npw,from_kg,from_istwfk,to_npw,to_kg,to_istwfk,n1,n2,n3,ndat1,from_cg,to_cg)

           ! Use cgwork as workspace array, results stored in bra_kq
           !g0_kq =  g0ibz_kq + g0bz_kq
           call wfd_copy_cg(wfd, ibsum_kq, ikq_ibz, spin, cgwork)
           call cgtk_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1,&
                            npw_kqirr, wfd%kdata(ikq_ibz)%kg_k,&
                            npw_kq, kg_kq, istwf_kqirr, istwf_kq, cgwork, bra_kq, work_ngfft, work)
           !bra_kq = zero
         end if

         ! h1kets_kq is the main responsible for the symmetry breaking
         !h1kets_kq = one
         do ipc=1,natom3
           ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
           !The array eig1_k contains:
           !
           ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
           ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
           do ib_k=1,nbcalc_ks
             call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bra_kq,h1kets_kq(:,:,ib_k,ipc),&
                  mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
             gkk_atm(:,ib_k,ipc) = [dotr, doti]
           end do
         end do

         ! Get gkk(kcalc, q, nu)
         call gkknu_from_atm(1, nbcalc_ks, 1, natom, gkk_atm, phfrq, displ_red, gkk_nu, num_smallw)

         ! Save data for Debye-Waller computation (performed outside the q-loop)
         ! dbwl_nu(2, nbcalc_ks, nbsum, natom3), gkk_nu(2, nbcalc_ks, natom3)
         if (isqzero) then
           gkk0_atm(:, :, ibsum_kq, :) = gkk_atm
           dbwl_nu(:, :, ibsum_kq, :) = gkk_nu
         end if

         ! Accumulate contribution to self-energy
         eig0mkq = ebands%eig(ibsum_kq,ikq_ibz,spin)
         f_mkq = ebands%occ(ibsum_kq,ikq_ibz,spin)
         if (nsppol == 1 .and. nspinor == 1 .and. nspden == 1) f_mkq = f_mkq * half
         weigth_q = sigma%wtq_k(iq_ibz)

         do nu=1,natom3
           ! Ignore acoustic or unstable modes.
           wqnu = phfrq(nu); if (wqnu < tol6) cycle

           ! Compute gaussian weights for Eliashberg function.
           if (sigma%gfw_nomega > 0) call gspline_eval(gspl, wqnu, sigma%gfw_nomega, sigma%gfw_mesh, dt_weights)

           do ib_k=1,nbcalc_ks
             band = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band,ik_ibz,spin)
             gkk2 = weigth_q * (gkk_nu(1,ib_k,nu) ** 2 + gkk_nu(2,ib_k,nu) ** 2)

             do it=1,sigma%ntemp
               !nqnu = one; f_mkq = one
               !write(std_out,*)wqnu,sigma%kTmesh(it)
               !write(std_out,*)eig0mkq,sigma%kTmesh(it),sigma%mu_e(it),eig0mkq-sigma%mu_e(it) / sigma%kTmesh(it)
               nqnu = nbe(wqnu, sigma%kTmesh(it), zero)
               ! TODO: Find good tolerances to treat limits
               !f_mkq = nfd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))
               !if (nsppol == 1 .and. nspinor == 1 .and. nspden == 1) f_mkq = f_mkq * half

               ! Accumulate Sigma(w=eKS) for state ib_k
               !my_ieta = sigma%ieta * sign(one, ebands%fermie - eig0nk)
               my_ieta = sigma%ieta
               cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + my_ieta) + &
                        (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + my_ieta)

               sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * cfact

               ! Accumulate contribution to Eliashberg functions (Fan term)
               if (sigma%gfw_nomega > 0) then
                 sigma%gfw_vals(:, it, 1, ib_k) = sigma%gfw_vals(:, it, 1, ib_k) + &
                   gkk2 * cfact * dt_weights(:, 1)
               end if
               if (sigma%has_nuq_terms) then
                 !TODO: in principle iq_ibz --> iq_bz
                 sigma%vals_nuq(it, ib_k, nu, iq_ibz, 1) = sigma%vals_nuq(it, ib_k, nu, iq_ibz, 1) + &
                   gkk2 * cfact / weigth_q
               end if

               ! Accumulate dSigma(w)/dw(w=eKS) derivative for state ib_k
#if 0
!old derivative
               ! This to avoid numerical instability
               !my_ieta = five * sigma%ieta
               cfact = -(nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + my_ieta) ** 2 &
                       -(nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + my_ieta) ** 2

               ! Try to avoid numerical instability
               !dka = (eig0nk - eig0mkq + wqnu + my_ieta) ** 2 *  (eig0nk - eig0mkq - wqnu + my_ieta) ** 2
               !cfact = -(nqnu + f_mkq) * (eig0nk - eig0mkq - wqnu + my_ieta) ** 2 &
               !        -(nqnu - f_mkq + one) * (eig0nk - eig0mkq + wqnu + my_ieta) ** 2
               !cfact = cfact / dka

               sigma%dvals_de0ks(it, ib_k) = sigma%dvals_de0ks(it, ib_k) + gkk2 * cfact
#else
               !cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + my_ieta) + &
               !         (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + my_ieta)
               !sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * cfact
               cfact = (eig0nk - eig0mkq + wqnu + my_ieta)
               gmod2 = cfact * dconjg(cfact)
               cfact = (eig0nk - eig0mkq - wqnu + my_ieta)
               hmod2 = cfact * dconjg(cfact)

               sigma%dvals_de0ks(it, ib_k) = sigma%dvals_de0ks(it, ib_k) + gkk2 * ( &
                 (nqnu + f_mkq)        * (gmod2 - two * (eig0nk - eig0mkq + wqnu) ** 2) / gmod2 ** 2 + &
                 (nqnu - f_mkq + one)  * (hmod2 - two * (eig0nk - eig0mkq - wqnu) ** 2) / hmod2 ** 2   &
               )
#endif

               ! Accumulate Sigma(w) for state ib_k
               if (sigma%nwr > 0) then
                 my_ieta = sigma%ieta
                 cfact_wr(:) = (nqnu + f_mkq      ) / (sigma%wrmesh_b(:,ib_k) - eig0mkq + wqnu + my_ieta) + &
                               (nqnu - f_mkq + one) / (sigma%wrmesh_b(:,ib_k) - eig0mkq - wqnu + my_ieta)
                 sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + gkk2 * cfact_wr(:)
               end if
             end do ! it

           end do ! ib_k
         end do ! nu

       end do ! ibsum_kq (sum over bands at k+q)

       ABI_FREE(bra_kq)
       ABI_FREE(cgwork)
       ABI_FREE(h1kets_kq)
       ABI_FREE(v1scf)
       ABI_FREE(vlocal1)

       call cwtime(cpu,wall,gflops,"stop")
       write(msg,'(2(a,i0),2(a,f8.2))')"q-point [",iq_ibz,"/",sigma%nqibz_k,"] completed. cpu:",cpu,", wall:",wall
       call wrtout(std_out, msg, do_flush=.True.)
     end do ! iq_ibz (sum over q-points)

     ABI_FREE(kets_k)
     ABI_FREE(gkk_atm)
     ABI_FREE(gkk_nu)
     ABI_FREE(distrib_bq)

     ! =========================
     ! Compute Debye-Waller term
     ! =========================
     call xmpi_sum(dbwl_nu, comm, ierr)
     call xmpi_sum(gkk0_atm, comm, ierr)

     ABI_MALLOC(deg_ibk, (nbcalc_ks))
     deg_ibk = 1
     ndeg = 1
#if 0
     ! Count number of degeneracies.
     do ib_k=2,nbcalc_ks
       ib = ib_k + bstart_ks - 1
       if ( abs(ebands%eig(ib,ik_ibz,spin) - ebands%eig(ib-1,ik_ibz,spin) ) > tol_enediff) ndeg = ndeg + 1
     end do

     ! Build degblock table.
     ABI_MALLOC(degblock, (2, ndeg))
     ndeg = 1; degblock(1, 1) = 1
     do ib_k=2,nbcalc_ks
       ib = ib_k + bstart_ks - 1
       if ( abs(ebands%eig(ib,ik_ibz,spin) - ebands%eig(ib-1,ik_ibz,spin) ) > tol_enediff) then
         degblock(2, ndeg) = ib_k - 1
         ndeg = ndeg + 1
         degblock(1, ndeg) = ib_k
       end if
     end do
     degblock(2, ndeg) = nbcalc_ks

     !ABI_MALLOC(gkk0_atm, (2, nbcalc_ks, nbsum, natom3))
     do ii=1,ndeg
       !write(std_out,*)"ideg, ", ii, ", degblock", degblock(:,ii)
       nstates = degblock(2, ii) - degblock(1, ii) + 1
       do ib_k=degblock(1, ii), degblock(2, ii)
         deg_ibk(ib_k) = nstates
       end do
       if (nstates == 1) continue
       ! Use dbwl_nu as workspace
       dbwl_nu(:, 1, :, :) = zero
       do ib_k=degblock(1, ii), degblock(2, ii)
         dbwl_nu(:, 1, :, :) = dbwl_nu(:, 1, :, :) + gkk0_atm(:, ib_k, :, :)
       end do
       !dbwl_nu(:, 1, :, :) = dbwl_nu(:, 1, :, :) / nstates
       do ib_k=degblock(1, ii), degblock(2, ii)
         gkk0_atm(:, ib_k, :, :) = dbwl_nu(:, 1, :, :)
       end do
     end do
     ABI_FREE(degblock)
     deg_ibk(:) = 1
#endif

     ABI_MALLOC(gdw2_mn, (nbsum, nbcalc_ks))
     ABI_MALLOC(tpp, (natom3, natom3))
     ABI_MALLOC(hka_mn, (natom3, nbsum, nbcalc_ks))

     qpt = zero
     call ifc_fourq(ifc, cryst, qpt, sqrt_phfrq0, displ_cart)
     where (sqrt_phfrq0 >= zero)
       sqrt_phfrq0 = sqrt(sqrt_phfrq0)
     else where
       sqrt_phfrq0 = zero
     end where
     d0mat = reshape(displ_cart, [2, natom3, natom3])
     ! cmat contains the displament vectors as complex array
     cmat = dcmplx(d0mat(1,:,:), d0mat(2,:,:))
     ! Multiply d by M to get e * M^{1/2}
     do nu=1,natom3
       do ip1=1,natom3
         idir1 = mod(ip1-1, 3) + 1; ipert1 = (ip1 - idir1) / 3 + 1
         rfact = ifc%amu(cryst%typat(ipert1)) * amu_emass !; rfact = one
         cmat(ip1, nu) = cmat(ip1, nu) * rfact
       end do
     end do

     ! Integral over IBZ. Note that here we can use IBZ(k=0).
     ! symsigma == 0 activates sum over full BZ for debugging purpose.
     ! For the time being, integrate over full BZ
     ! TODO Bug if symsigma 0 ?
     nq = sigma%nqibz
     if (sigma%symsigma == 0) nq = sigma%nqbz
     !if (sigma%symsigma == 0) nq = sigma%nqibz_k
     !do iq_ibz=1,sigma%nqibz
     !nq = sigma%nqbz
     do iq_ibz=1,nq
       if (mod(iq_ibz, nprocs) /= my_rank) cycle  ! MPI parallelism
#if 0
       qpt = sigma%qbz(:,iq_ibz); weigth_q = one / sigma%nqbz
#else
       if (abs(sigma%symsigma) == 1) then
         qpt = sigma%qibz(:,iq_ibz); weigth_q = sigma%wtq(iq_ibz)
         !qpt = sigma%qibz_k(:,iq_ibz); weigth_q = sigma%wtq_k(iq_ibz)
       else
         qpt = sigma%qbz(:,iq_ibz); weigth_q = one / sigma%nqbz
         !qpt = sigma%qibz_k(:,iq_ibz); weigth_q = sigma%wtq_k(iq_ibz)
       end if
#endif

       ! Get phonons for this q-point.
       call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart)

       ! Compute hka_mn matrix with shape: (natom3, nbsum, nbcalc_ks))
       ! Needed for Giustino's equation.
       ! dbwl_nu(2, nbcalc_ks, nbsum, natom3)
       do in=1,nbcalc_ks
         do im=1,nbsum
           cvec1 = dcmplx(sqrt_phfrq0(:) * dbwl_nu(1, in, im, :), &
                          sqrt_phfrq0(:) * dbwl_nu(2, in, im, :))
           do ii=1,natom3
             cvec2 = cmat(ii,:)
             !hka_mn(ii, im, in) = xdotu(natom3, cvec2, 1, cvec1, 1)
             hka_mn(ii, im, in) = dot_product(dconjg(cvec2), cvec1)
             !hka_mn(ii, im, in) = dconjg(hka_mn(ii, im, in))
             !write(std_out,*)"hka_mn: ",hka_mn(ii, im, in)
           end do
         end do
       end do

       ! Sum over modes for this q-point.
       do nu=1,natom3
         ! Ignore acoustic or unstable modes.
         wqnu = phfrq(nu); if (wqnu < tol6) cycle

         ! Compute gaussian weights for Eliashberg function.
         if (sigma%gfw_nomega > 0) call gspline_eval(gspl, wqnu, sigma%gfw_nomega, sigma%gfw_mesh, dt_weights)

         ! Compute T_pp'(q,nu) matrix in cartesian coordinates.
         do ip2=1,natom3
           idir2 = mod(ip2-1, 3) + 1; ipert2 = (ip2 - idir2) / 3 + 1
           do ip1=1,natom3
             idir1 = mod(ip1-1, 3) + 1; ipert1 = (ip1 - idir1) / 3 + 1
             ! (k,a) (k,a')* + (k',a) (k',a')*
             dka   = dcmplx(displ_cart(1, idir1, ipert1, nu), displ_cart(2, idir1, ipert1, nu))
             dkap  = dcmplx(displ_cart(1, idir2, ipert1, nu), displ_cart(2, idir2, ipert1, nu))
             dkpa  = dcmplx(displ_cart(1, idir1, ipert2, nu), displ_cart(2, idir1, ipert2, nu))
             dkpap = dcmplx(displ_cart(1, idir2, ipert2, nu), displ_cart(2, idir2, ipert2, nu))
             tpp(ip1,ip2) = dka * dconjg(dkap) + dkpa * dconjg(dkpap)
             !write(std_out,*)"tpp: ",tpp(ip1, ip2)
           end do
         end do

         ! Giustino's equation in RMP
         ! gdw2_mn = diag(conjg(H) T H) / (2 w(qu,nu))
         ! tpp(natom3, natom3), hka_mn(natom3, nbsum, nbcalc_ks)
         ABI_MALLOC(wmat1, (natom3, nbsum, nbcalc_ks))
         call zgemm("N", "N", natom3, nbsum*nbcalc_ks, natom3, cone, tpp, natom3, hka_mn, natom3, czero, wmat1, natom3)

         do ibsum=1,nbsum
           do ib_k=1,nbcalc_ks
             !write(std_out,*)"maxval(wmat)",maxval(abs(wmat1(:, ibsum, ib_k)))
             cfact = dot_product(hka_mn(:, ibsum, ib_k), wmat1(:, ibsum, ib_k))
             !cfact = dot_product(wmat1(:, ibsum, ib_k), wmat1(:, ibsum, ib_k))
             !cfact = xdotc(natom3, hka_mn(:, ibsum, ib_k), 1, wmat1(:, ibsum, ib_k), 1)
             !cfact = xdotu(natom3, hka_mn(:, ibsum, ib_k), 1, wmat1(:, ibsum, ib_k), 1)
             gdw2_mn(ibsum, ib_k) = real(cfact) / (two * wqnu)
             !write(std_out, *)"gdw2_mn: ", gdw2_mn(ibsum, ib_k)
           end do
         end do
         ABI_FREE(wmat1)
         !gdw2_mn = zero

         ! Get phonon occupation for all temperatures.
         nqnu_tlist = nbe(wqnu, sigma%kTmesh(:), zero)

         ! Sum over bands and add (static) DW contribution for the different temperatures.
         do ibsum=1,nbsum
           eig0mk = ebands%eig(ibsum, ik_ibz, spin)
           do ib_k=1,nbcalc_ks
             band = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band, ik_ibz, spin)
             ! Handle n == m and degenerate states (either ignore or add broadening)
             cplx_ediff = (eig0nk - eig0mk)
             !if (abs(cplx_ediff) < tol6) cplx_ediff = cplx_ediff + sigma%ieta
             if (abs(cplx_ediff) < tol6) cycle

#if 1
             ! Compute DW term following XG paper. Check prefactor.
             ! gkk0_atm(2, nbcalc_ks, nbsum, natom3)
             ! previous version
             gdw2_mn(ibsum, ib_k) = zero
             do ip2=1,natom3
               do ip1=1,natom3
                 gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) + tpp(ip1,ip2) * ( &
                   gkk0_atm(1, ib_k, ibsum, ip1) * gkk0_atm(1, ib_k, ibsum, ip2) + &
                   gkk0_atm(2, ib_k, ibsum, ip1) * gkk0_atm(2, ib_k, ibsum, ip2) + &
                   gkk0_atm(1, ib_k, ibsum, ip2) * gkk0_atm(1, ib_k, ibsum, ip1) + &
                   gkk0_atm(2, ib_k, ibsum, ip2) * gkk0_atm(2, ib_k, ibsum, ip1) &
                  )
               end do

               !do ip1=1,natom3
               !  gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) + tpp(ip1,ip2) * ( &
               !    gkk0_atm(1, ib_k, ibsum, ip1) * gkk0_atm(1, ib_k, ibsum, ip2) + &
               !    gkk0_atm(2, ib_k, ibsum, ip1) * gkk0_atm(2, ib_k, ibsum, ip2)  ) + &
               !    tpp(ip2, ip1) * ( &
               !    gkk0_atm(1, ib_k, ibsum, ip2) * gkk0_atm(1, ib_k, ibsum, ip1) + &
               !    gkk0_atm(2, ib_k, ibsum, ip2) * gkk0_atm(2, ib_k, ibsum, ip1) &
               !   )
               !end do

             end do
             !gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) * two
             gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) * two / deg_ibk(ib_k)
             !write(std_out,*)"gdw2_mn: ",gdw2_mn(ibsum, ib_k)
#endif
             ! dbwl_nu(2, nbcalc_ks, nbsum, natom3), gkk_nu(2, nbcalc_ks, natom3)

             ! accumulate DW for each T, add it to Sigma(e0) and Sigma(w) as well
             do it=1,sigma%ntemp
               cfact = - weigth_q * gdw2_mn(ibsum, ib_k) * (two * nqnu_tlist(it) + one)  / cplx_ediff
               rfact = real(cfact)
               sigma%dw_vals(it, ib_k) = sigma%dw_vals(it, ib_k) + rfact
               sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + rfact
               if (sigma%nwr > 0) sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + rfact

               ! Optionally, accumulate contribution to Eliashberg functions (DW term)
               if (sigma%gfw_nomega > 0) then
                 sigma%gfw_vals(:, it, 2, ib_k) = sigma%gfw_vals(:, it, 2, ib_k) + &
                   rfact * dt_weights(:, 1)
               end if

               if (sigma%has_nuq_terms) then
                 ! TODO: in principle iq_ibz --> iq_bz if symsigma == 0
                 sigma%vals_nuq(it, ib_k, nu, iq_ibz, 2) = sigma%vals_nuq(it, ib_k, nu, iq_ibz, 2) + &
                   gkk2 * rfact / weigth_q
               end if
             end do

           end do
         end do

       end do ! nu
     end do ! iq_ibz

     ABI_FREE(deg_ibk)
     ABI_FREE(gdw2_mn)
     ABI_FREE(tpp)
     ABI_FREE(hka_mn)
     ABI_FREE(dbwl_nu)
     ABI_FREE(gkk0_atm)

     ! Collect results inside comm and write results for this (k-point, spin) to NETCDF file.
     call sigmaph_gather_and_write(sigma, ebands, ikcalc, spin, comm)
   end do ! spin

   ABI_FREE(kg_k)
   ABI_FREE(kg_kq)
   ABI_FREE(ylm_k)
   ABI_FREE(ylm_kq)
   ABI_FREE(ylmgr_kq)
 end do !ikcalc

 call cwtime(cpu_all, wall_all, gflops_all, "stop")
 call wrtout(std_out, "Computation of Sigma_eph completed", do_flush=.True.)
 call wrtout(std_out, sjoin("Total wall-time:", sec2str(cpu_all), ", Total cpu time:", sec2str(wall_all), ch10, ch10))

 ! Free memory
 ABI_FREE(gvnl1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(nqnu_tlist)
 if (sigma%nwr > 0) then
   ABI_FREE(cfact_wr)
 end if
 if (sigma%gfw_nomega > 0) then
   ABI_FREE(dt_weights)
 end if

 call gspline_free(gspl)
 call destroy_hamiltonian(gs_hamkq)
 call sigmaph_free(sigma)
 call wfd_free(wfd)
 call pawcprj_free(cwaveprj0)
 ABI_DT_FREE(cwaveprj0)

end subroutine sigmaph
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/gkknu_from_atm
!! NAME
!!  gkknu_from_atm
!!
!! FUNCTION
!!  Transform the gkk matrix elements from (atom, red_direction) basis to phonon-mode basis.
!!
!! INPUTS
!!  nb1,nb2=Number of bands in gkk_atm matrix.
!!  nk=Number of k-points (usually 1)
!!  natom=Number of atoms.
!!  gkk_atm(2,nb1,nb2,3*natom)=EPH matrix elements in the atomic basis.
!!  phfrq(3*natom)=Phonon frequencies in Ha
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduced coordinates.
!!
!! OUTPUT
!!  gkk_nu(2,nb1,nb2,3*natom)=EPH matrix elements in the phonon-mode basis.
!!  num_smallw=Number of negative/too small frequencies that have been ignored
!!    by setting the corresponding gkk_nu to zero.
!!
!! PARENTS
!!      m_sigmaph
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

   ! Transform the gkk from atom, red_dir basis to phonon mode basis
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
!!  Fermi-Dirac statistic 1 / [(exp((e - mu)/ KT) + 1]
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

!Local variables ------------------------------
 real(dp) :: ee_mu,arg
! *************************************************************************

 ee_mu = ee - mu

 !TODO: Find good tols.
 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol6) then
   arg = ee_mu / kT
   nfd = one / (exp(arg) + one)
 else
   ! Heaviside
   if (ee_mu > zero) then
     nfd = zero
   else if (ee < zero) then
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
!!   Bose-Einstein statistic  1 / [(exp((e - mu)/ KT) - 1]
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

!Local variables ------------------------------
 real(dp) :: ee_mu, arg
! *************************************************************************

 ee_mu = ee - mu

 !TODO: Find good tols.
 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol12) then
   arg = ee_mu / kT
   if (arg > tol12 .and. arg < 600._dp) then
     nbe = one / (exp(arg) - one)
   else
     nbe = zero
   end if
 else
   ! No condensate for T --> 0
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
!!  ecut=Cutoff energy for wavefunctions.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  dtfil<datafiles_type>=variables related to files.
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type (sigmaph_t) function sigmaph_new(dtset, ecut, cryst, ebands, ifc, dtfil, comm) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 real(dp),intent(in) :: ecut
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),intent(in) :: ifc
 type(datafiles_type),intent(in) :: dtfil

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0,brav1=1,occopt3=3,qptopt1=1
 integer :: my_rank,ik,my_nshiftq,my_mpw,cnt,nprocs,iq_ibz,ik_ibz,ndeg
 integer :: onpw,ii,ipw,ierr,it,spin,gap_err,ikcalc,gw_qprange,ibstop
 integer :: nk_found,ifo,jj
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 real(dp),parameter :: spinmagntarget=-99.99_dp,tol_enediff=0.001_dp*eV_Ha
 real(dp) :: dksqmax,ph_wstep
 character(len=500) :: msg
 logical :: changed,found
 type(ebands_t) :: tmp_ebands
 type(gaps_t) :: gaps
!arrays
 integer :: qptrlatt(3,3),indkk_k(1,6),my_gmax(3),kpos(6)
 integer :: val_indeces(ebands%nkpt, ebands%nsppol)
 integer,allocatable :: gtmp(:,:),degblock(:,:)
 real(dp) :: my_shiftq(3,1),kk(3),kq(3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Copy important dimensions.
 new%nsppol = ebands%nsppol; new%nspinor = ebands%nspinor

 ! Broadening parameter from zcut
 new%ieta = + j_dpc * dtset%zcut

 ! Build (linear) mesh of K * temperatures. tsmesh(1:3) = [start, step, num]
 new%ntemp = nint(dtset%tmesh(3))
 ABI_CHECK(new%ntemp > 0, "ntemp <= 0")
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 !write(std_out, *)"dtset%tmesh", dtset%tmesh
 new%kTmesh = arth(dtset%tmesh(1), dtset%tmesh(2), new%ntemp) * kb_HaK

 ! Number of bands included in self-energy summation.
 new%nbsum = dtset%mband

 gap_err = get_gaps(ebands, gaps)
 call gaps_print(gaps, unit=std_out)
 call ebands_report_gap(ebands, unit=std_out)
 val_indeces = get_valence_idx(ebands)

 ! Compute the chemical potential at the different physical temperatures with Fermi-Dirac.
 ABI_MALLOC(new%mu_e, (new%ntemp))
 new%mu_e = ebands%fermie

 ! TODO T == 0 >> SIGSEGV
#if 0
 call ebands_copy(ebands, tmp_ebands)
 do it=1,new%ntemp
   call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(it), spinmagntarget, dtset%prtvol)
   new%mu_e(it) = tmp_ebands%fermie
 end do
 call ebands_free(tmp_ebands)
#endif

 ! Frequency mesh for sigma(w) and spectral function.
 ! TODO: Use GW variables but change default
 !dtset%freqspmin
 !dtset%freqspmax
 !dtset%nfreqsp
 new%nwr = 201
 ABI_CHECK(mod(new%nwr, 2) == 1, "nwr should be odd!")
 new%wr_step = 0.05 * eV_Ha

 ! Define q-mesh for integration.
 ! Either q-mesh from DDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation if q not in DDB)
 new%ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   new%ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems wo inversion
 qptrlatt = 0; qptrlatt(1,1) = new%ngqpt(1); qptrlatt(2,2) = new%ngqpt(2); qptrlatt(3,3) = new%ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
   new%nqibz, new%qibz, new%wtq, new%nqbz, new%qbz)

 ! Select k-point and bands where corrections are wanted
 ! if symsigma == 1, we have to include all degenerate states in the set
 ! because the final QP corrections will be obtained by averaging the results in the degenerate subspace.
 ! k-point and bands where corrections are wanted
 ! We initialize IBZ(k) here so that we have all the basic dimensions of the run and it's possible
 ! to distribuite the calculations among processors.
 new%symsigma = dtset%symsigma; new%timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ! TODO: debug gw_qprange > 0. Rename variable
 if (dtset%nkptgw /= 0) then
   ! Treat the k-points and bands specified in the input file.
   call wrtout(std_out, "Getting list of k-points for self-energy from kptgw and bdgw.")
   new%nkcalc = dtset%nkptgw
   ABI_MALLOC(new%kcalc, (3, new%nkcalc))
   ABI_MALLOC(new%bstart_ks, (new%nkcalc, new%nsppol))
   ABI_MALLOC(new%nbcalc_ks, (new%nkcalc, new%nsppol))

   new%kcalc = dtset%kptgw(:,1:new%nkcalc)
   do spin=1,new%nsppol
     new%bstart_ks(:,spin) = dtset%bdgw(1,1:new%nkcalc,spin)
     new%nbcalc_ks(:,spin) = dtset%bdgw(2,1:new%nkcalc,spin) - dtset%bdgw(1,1:new%nkcalc,spin) + 1
   end do

   ! Consistency check on bdgw and dtset%mband
   ierr = 0
   do spin=1,new%nsppol
     do ikcalc=1,new%nkcalc
       if (dtset%bdgw(2,ikcalc,spin) > dtset%mband) then
         ierr = ierr + 1
         write(msg,'(a,2i0,2(a,i0))')&
          "For (k, s) ",ikcalc,spin," bdgw= ",dtset%bdgw(2,ikcalc,spin), " > mband=",dtset%mband
         MSG_WARNING(msg)
       end if
     end do
   end do
   ABI_CHECK(ierr == 0, "Not enough bands in WFK file. See messages above. Aborting now.")

 else
   ! Use qp_range to select the interesting k-points and the corresponing bands.
   !
   !    0 --> Compute the QP corrections only for the fundamental and the optical gap.
   ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
   !           bands above and below the Fermi level.
   ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
   !          Include all occupied states and `num` empty states.

   call wrtout(std_out, "nkptgw == 0 ==> Automatic selection of k-points and bands for the corrections.")
   gw_qprange = dtset%gw_qprange

   if (gap_err /=0 .and. gw_qprange == 0) then
     msg = "Problem while computing the fundamental and optical gap (likely metal). Will replace gw_qprange=0 with gw_qprange=1"
     MSG_WARNING(msg)
     gw_qprange = 1
   end if

   if (gw_qprange /= 0) then
     ! Include all the k-points in the IBZ.
     ! Note that kcalc == ebands%kptns so we can use a single ik index in the loop over k-points.
     ! No need to map kcalc onto ebands%kptns.
     new%nkcalc = ebands%nkpt
     ABI_MALLOC(new%kcalc, (3, new%nkcalc))
     ABI_MALLOC(new%bstart_ks, (new%nkcalc, new%nsppol))
     ABI_MALLOC(new%nbcalc_ks, (new%nkcalc, new%nsppol))

     new%kcalc = ebands%kptns
     new%bstart_ks = 1; new%nbcalc_ks = 6

     if (gw_qprange > 0) then
       ! All k-points: Add buffer of bands above and below the Fermi level.
       do spin=1,new%nsppol
         do ik=1,new%nkcalc
           new%bstart_ks(ik,spin) = max(val_indeces(ik,spin) - gw_qprange, 1)
           ii = max(val_indeces(ik,spin) + gw_qprange + 1, dtset%mband)
           new%nbcalc_ks(ik,spin) = ii - new%bstart_ks(ik,spin) + 1
         end do
       end do

     else
       ! All k-points: include all occupied states and -gw_qprange empty states.
       new%bstart_ks = 1
       do spin=1,new%nsppol
         do ik=1,new%nkcalc
           new%nbcalc_ks(ik,spin) = min(val_indeces(ik,spin) - gw_qprange, dtset%mband)
         end do
       end do
     end if

   else
     ! gw_qprange is not specified in the input.
     ! Include the optical and the fundamental KS gap.
     ! The main problem here is that kptgw and nkptgw do not depend on the spin and therefore
     ! we have compute the union of the k-points where the fundamental and the optical gaps are located.
     !
     ! Find the list of `interesting` kpoints.
     ABI_CHECK(gap_err == 0, "gw_qprange 0 cannot be used because I cannot find the gap (gap_err !=0)")
     nk_found = 1; kpos(1) = gaps%fo_kpos(1,1)

     do spin=1,new%nsppol
       do ifo=1,3
         ik = gaps%fo_kpos(ifo, spin)
         found = .False.; jj = 0
         do while (.not. found .and. jj < nk_found)
           jj = jj + 1; found = (kpos(jj) == ik)
         end do
         if (.not. found) then
           nk_found = nk_found + 1; kpos(nk_found) = ik
         end if
       end do
     end do

     ! Now we can define the list of k-points and the bands range.
     new%nkcalc = nk_found
     ABI_MALLOC(new%kcalc, (3, new%nkcalc))
     ABI_MALLOC(new%bstart_ks, (new%nkcalc, new%nsppol))
     ABI_MALLOC(new%nbcalc_ks, (new%nkcalc, new%nsppol))

     do ii=1,new%nkcalc
       ik = kpos(ii)
       new%kcalc(:,ii) = ebands%kptns(:,ik)
       do spin=1,new%nsppol
         new%bstart_ks(ii,spin) = val_indeces(ik,spin)
         new%nbcalc_ks(ii,spin) = 2
       end do
     end do
   end if ! gw_qprange

 end if ! nkptgw /= 0

 call gaps_free(gaps)

 ! The k-point and the symmetries relating the BZ k-point to the IBZ.
 ABI_MALLOC(new%kcalc2ibz, (new%nkcalc, 6))
 if (new%symsigma == 1) then
   ABI_DT_MALLOC(new%degtab, (new%nkcalc, new%nsppol))
 end if

 do ikcalc=1,new%nkcalc
   kk = new%kcalc(:,ikcalc)
   call listkk(dksqmax,cryst%gmet,indkk_k,ebands%kptns,kk,ebands%nkpt,1,cryst%nsym,&
      1,cryst%symafm,cryst%symrel,new%timrev,use_symrec=.False.)

   new%kcalc2ibz(ikcalc, :) = indkk_k(1, :)

   ik_ibz = indkk_k(1,1) !; isym_k = indkk_k(1,2)
   !trev_k = indkk_k(1, 6); g0_k = indkk_k(1, 3:5)
   !isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !kk_ibz = ebands%kptns(:,ik_ibz)

   if (dksqmax > tol12) then
      write(msg, '(4a,es16.6,7a)' )&
       "The WFK file cannot be used to compute self-energy corrections at kpoint: ",ktoa(kk),ch10,&
       "the k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
       "Q-mesh: ",ltoa(new%ngqpt),", K-mesh (from kptrlatt) ",ltoa(get_diag(dtset%kptrlatt)), ch10, &
       'Action: check your WFK file and (k,q) point input variables'
      MSG_ERROR(msg)
   end if

   ! Find the little-group of the k-point and initialize the weights for BZ integration.
   ! Examine the symmetries of the k-point
   ! TODO: See timerev_k
   !call littlegroup_q(cryst%nsym,kk,symk,cryst%symrec,cryst%symafm,timerev_k,prtvol=dtset%prtvol)
   !call wrtout(std_out, sjoin("Will compute", itoa(sigma%nqibz_k), "q-points in the IBZ(k)"))

   ! We will have to average the QP corrections over degenerate states if symsigma=1 is used.
   ! Here we make sure that all the degenerate states are included.
   ! Store also band indices of the degenerate sets, used to average final results.
   if (abs(new%symsigma) == 1) then
     do spin=1,new%nsppol
       ibstop = new%bstart_ks(ikcalc,spin) + new%nbcalc_ks(ikcalc,spin) - 1
       call enclose_degbands(ebands, ik_ibz, spin, new%bstart_ks(ikcalc,spin), ibstop, &
         changed, tol_enediff, degblock=degblock)
       if (changed) then
         new%nbcalc_ks(ikcalc,spin) = ibstop - new%bstart_ks(ikcalc,spin) + 1
         write(msg,'(2(a,i0),2a,2(1x,i0))')&
           "Not all the degenerate states at ikcalc= ",ikcalc,", spin= ",spin,ch10,&
           "were included in the bdgw set. bdgw has been changed to: ",new%bstart_ks(ikcalc,spin),ibstop
         MSG_COMMENT(msg)
       end if
       ! Store band indices used for averaging.
       ndeg = size(degblock, dim=2)
       ABI_DT_MALLOC(new%degtab(ikcalc, spin)%bids, (ndeg))
       do ii=1,ndeg
         cnt = degblock(2, ii) - degblock(1, ii) + 1
         ABI_DT_MALLOC(new%degtab(ikcalc, spin)%bids(ii)%vals, (cnt))
         new%degtab(ikcalc, spin)%bids(ii)%vals = [(jj, jj=degblock(1, ii), degblock(2, ii))]
       end do
       ABI_FREE(degblock)
     end do
   end if ! symsigma

 end do

 ! Now we can finally compute max_nbcalc
 new%max_nbcalc = maxval(new%nbcalc_ks)

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! used to symmetrize the wavefunctions in G-space.
 ! TODO: Should loop over IBZ(k)
 new%mpw = 0; new%gmax = 0; cnt = 0
 do ik=1,new%nkcalc
   kk = new%kcalc(:, ik)
   !call sigmaph_setup_kcalc(self, crystal, ik)
   do iq_ibz=1,new%nqbz
   !do iq_ibz=1,new%nqibz_k
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle
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

 ! ================================================================
 ! Allocate arrays used to store final results and set them to zero
 ! ================================================================
 ABI_CALLOC(new%vals_e0ks, (new%ntemp, new%max_nbcalc))
 ABI_CALLOC(new%dvals_de0ks, (new%ntemp, new%max_nbcalc))
 ABI_CALLOC(new%dw_vals, (new%ntemp, new%max_nbcalc))

 ! Frequency dependent stuff
 if (new%nwr > 0) then
   ABI_CALLOC(new%vals_wr, (new%nwr, new%ntemp, new%max_nbcalc))
   ABI_CALLOC(new%wrmesh_b, (new%nwr, new%max_nbcalc))
 end if

 ! Prepare calculation of generalized Eliashberg functions.
 new%gfw_nomega = 0; ph_wstep = one/Ha_cmm1

 if (new%gfw_nomega /= 0 .and. ph_wstep > zero) then
   new%gfw_nomega = nint((ifc%omega_minmax(2) - ifc%omega_minmax(1) ) / ph_wstep) + 1
   ABI_MALLOC(new%gfw_mesh, (new%gfw_nomega))
   new%gfw_mesh = arth(ifc%omega_minmax(1), ph_wstep, new%gfw_nomega)
   ABI_MALLOC(new%gfw_vals, (new%gfw_nomega, new%ntemp, 2, new%max_nbcalc))
 end if

 new%has_nuq_terms = .False.
 if (new%has_nuq_terms) then
   ABI_CALLOC(new%vals_nuq, (new%ntemp, new%max_nbcalc, 3*cryst%natom, new%nqbz, 2))
 end if

 ! Open netcdf file (only master work for the time being because cannot assume HDF5 + MPI-IO)
 ! This could create problems if MPI parallelism over (spin, nkptgw) ...
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   ! Master creates the netcdf file used to store the results of the calculation.
   NCF_CHECK(nctk_open_create(new%ncid, strcat(dtfil%filnam_ds(4), "_SIGEPH.nc"), xmpi_comm_self))
   ncid = new%ncid

   NCF_CHECK(crystal_ncwrite(cryst, ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))

   ! Add sigma_eph dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nkcalc", new%nkcalc), nctkdim_t("max_nbcalc", new%max_nbcalc), &
     nctkdim_t("nsppol", new%nsppol), nctkdim_t("ntemp", new%ntemp), nctkdim_t("natom3", 3 * cryst%natom), &
     nctkdim_t("nqibz", new%nqibz), nctkdim_t("nqbz", new%nqbz)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   if (new%nwr > 0) then
     NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nwr", new%nwr)]))
   end if
   if (new%gfw_nomega > 0) then
     NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("gfw_nomega", new%gfw_nomega)]))
   end if

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "symsigma", "nbsum"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "eph_intmeth", "eph_transport", "symdynmat"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "ph_intmeth"], defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "eta", "wr_step"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "ph_wstep", "ph_smear"])
   NCF_CHECK(ncerr)

   ! Define arrays for results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngqpt", "int", "three"), &
     nctkarr_t("eph_ngqpt_fine", "int", "three"), &
     nctkarr_t("ddb_ngqpt", "int", "three"), &
     nctkarr_t("ph_ngqpt", "int", "three"), &
     nctkarr_t("bstart_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("nbcalc_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("kcalc", "dp", "three, nkcalc"), &
     nctkarr_t("kcalc2ibz", "int", "nkcalc, six"), &
     nctkarr_t("kTmesh", "dp", "ntemp"), &
     nctkarr_t("mu_e", "dp", "ntemp"), &
     nctkarr_t("vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("dvals_de0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("dw_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("qpadb_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ze0_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ks_gaps", "dp", "nkcalc, nsppol"), &
     nctkarr_t("qpadb_gaps", "dp", "ntemp, nkcalc, nsppol"), &
     nctkarr_t("qp_gaps", "dp", "ntemp, nkcalc, nsppol") &
   ])
   NCF_CHECK(ncerr)

   ! TODO: Check the mesh
   if (new%nwr > 0) then
     ! These arrays get two extra dimensions on file (nkcalc, nsppol).
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("wrmesh_b", "dp", "nwr, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
   end if
   if (new%gfw_nomega > 0) then
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("gfw_mesh", "dp", "gfw_nomega"), &
       nctkarr_t("gfw_vals", "dp", "two, gfw_nomega, ntemp, two, max_nbcalc, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
   end if
   if (new%has_nuq_terms) then
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("vals_nuq", "dp", "two, ntemp, max_nbcalc, natom3, nqbz, two, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
   end if

   if (new%nwr > 1) then
     ! Make room for the spectral function.
     ncerr = nctk_def_arrays(ncid, [nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol")])
     NCF_CHECK(ncerr)
   end if

   ! ======================================================
   ! Write data that do not depend on the (kpt, spin) loop.
   ! ======================================================
   NCF_CHECK(nctk_set_datamode(ncid))
   ncerr = nctk_write_iscalars(ncid, &
       [character(len=nctk_slen) :: "symsigma", "nbsum", "symdynmat", "ph_intmeth", "eph_intmeth", "eph_transport"], &
       [new%symsigma, new%nbsum, dtset%symdynmat, dtset%ph_intmeth, dtset%eph_intmeth, dtset%eph_transport])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, &
     [character(len=nctk_slen) :: "eta", "wr_step"], &
     [aimag(new%ieta), new%wr_step])
   NCF_CHECK(ncerr)

   ncerr = nctk_write_dpscalars(ncid, &
     [character(len=nctk_slen) :: "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"], &
     [dtset%eph_fsewin, dtset%eph_fsmear, dtset%eph_extrael, dtset%eph_fermie, dtset%ph_wstep, dtset%ph_smear])
   NCF_CHECK(ncerr)

   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt"), new%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_ngqpt_fine"), dtset%eph_ngqpt_fine))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ddb_ngqpt"), dtset%ddb_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ph_ngqpt"), dtset%ph_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "bstart_ks"), new%bstart_ks))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nbcalc_ks"), new%nbcalc_ks))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc"), new%kcalc))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ibz"), new%kcalc2ibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), new%kTmesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mu_e"), new%mu_e))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eta"), aimag(new%ieta)))
   if (new%gfw_nomega > 0) then
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gfw_mesh"), new%gfw_mesh))
   end if
   NCF_CHECK(nf90_close(ncid))
 end if ! master

 ! Now reopen the file in parallel.
 call xmpi_barrier(comm)
 NCF_CHECK(nctk_open_modify(new%ncid, strcat(dtfil%filnam_ds(4), "_SIGEPH.nc"), xmpi_comm_self))
 !NCF_CHECK(nctk_open_modify(new%ncid, strcat(dtfil%filnam_ds(4), "_SIGEPH.nc"), comm))
 NCF_CHECK(nctk_set_datamode(new%ncid))
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
!!      m_sigmaph
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

!Local variables-------------------------------
 integer :: ii,jj,ideg

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
 if (allocated(self%qibz)) then
   ABI_FREE(self%qibz)
 end if
 if (allocated(self%wtq)) then
   ABI_FREE(self%wtq)
 end if
 if (allocated(self%qibz_k)) then
   ABI_FREE(self%qibz_k)
 end if
 if (allocated(self%wtq_k)) then
   ABI_FREE(self%wtq_k)
 end if
 if (allocated(self%gfw_mesh)) then
   ABI_FREE(self%gfw_mesh)
 end if

 ! complex
 if (allocated(self%vals_e0ks)) then
   ABI_FREE(self%vals_e0ks)
 end if
 if (allocated(self%dvals_de0ks)) then
   ABI_FREE(self%dvals_de0ks)
 end if
 if (allocated(self%dw_vals)) then
   ABI_FREE(self%dw_vals)
 end if
 if (allocated(self%vals_wr)) then
   ABI_FREE(self%vals_wr)
 end if
 if (allocated(self%gfw_vals)) then
   ABI_FREE(self%gfw_vals)
 end if
 if (allocated(self%vals_nuq)) then
   ABI_FREE(self%vals_nuq)
 end if

 ! types.
 if (allocated(self%degtab)) then
    do jj=1,size(self%degtab, dim=2)
      do ii=1,size(self%degtab, dim=1)
         do ideg=1,size(self%degtab(ii, jj)%bids)
           ABI_FREE(self%degtab(ii, jj)%bids(ideg)%vals)
         end do
         ABI_DT_FREE(self%degtab(ii, jj)%bids)
      end do
    end do
    ABI_DT_FREE(self%degtab)
 end if

 ! Close netcdf file.
#ifdef HAVE_NETCDF
 if (self%ncid /= nctk_noid) then
   NCF_CHECK(nf90_close(self%ncid))
 end if
#endif

end subroutine sigmaph_free
!!***

!!****f* m_sigmaph/sigmaph_setup_kcalc
!! NAME
!!  sigmaph_setup_kcalc
!!
!! FUNCTION
!!  Prepare calculations of self-energy matrix elements for ikcalc index.
!!
!! INPUTS
!!  crystal<crystal_t> = Crystal structure.
!!  ikcalc=Index of the k-point to compute.
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_setup_kcalc(self, cryst, ikcalc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_setup_kcalc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc
 type(crystal_t),intent(in) :: cryst
 type(sigmaph_t),intent(inout) :: self

!Local variables-------------------------------
 type(lgroup_t) :: lgk

! *************************************************************************

 if (allocated(self%qibz_k)) then
   ABI_FREE(self%qibz_k)
 end if
 if (allocated(self%wtq_k)) then
   ABI_FREE(self%wtq_k)
 end if

 if (self%symsigma == 0) then
   ! Do not use symmetries in BZ sum_q --> nqibz_k == nqbz
   self%nqibz_k = self%nqbz
   ABI_MALLOC(self%qibz_k, (3, self%nqibz_k))
   ABI_MALLOC(self%wtq_k, (self%nqibz_k))
   self%qibz_k = self%qbz; self%wtq_k = one / self%nqbz

 else if (abs(self%symsigma) == 1) then
   ! Use the symmetries of the little group
   lgk = lgroup_new(cryst, self%kcalc(:, ikcalc), self%timrev, self%nqbz, self%qbz, self%nqibz, self%qibz)

   self%nqibz_k = lgk%nibz
   ABI_MALLOC(self%qibz_k, (3, self%nqibz_k))
   ABI_MALLOC(self%wtq_k, (self%nqibz_k))
   self%qibz_k = lgk%ibz; self%wtq_k = lgk%weights
   call lgroup_free(lgk)
 else
   MSG_ERROR(sjoin("Wrong symsigma:", itoa(self%symsigma)))
 end if

end subroutine sigmaph_setup_kcalc
!!***

!!****f* m_sigmaph/sigmaph_gather_and_write
!! NAME
!!  sigmaph_gather_and_write
!!
!! FUNCTION
!!  Gather results from the MPI processes. Then master:
!!
!!      1. Computes QP energies, Z factor and spectral function (if required).
!!      2. Saves results to file.
!!
!! INPUTS
!!  ebands<ebands_t>=KS band energies.
!!  ikcalc=Index of the computed k-point
!!  spin=Spin index.
!!  comm=MPI communicator.
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_gather_and_write(self, ebands, ikcalc, spin, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_gather_and_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, spin, comm
 type(sigmaph_t),target,intent(inout) :: self
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 integer,parameter :: master=0
 integer :: ideg,ib,it,ii,iw,nstates,ierr,my_rank,band,ik_ibz,ibc,ib_val,ib_cond
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_adb,qpe_adb_val,qpe_adb_cond
 real(dp) :: smrt,alpha,beta,e0pde(9)
 complex(dpc) :: sig0c,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2
 character(len=500) :: msg
!arrays
 integer :: shape3(3),shape4(4),shape5(5),shape6(6)
 real(dp), ABI_CONTIGUOUS pointer :: rdata3(:,:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:),rdata6(:,:,:,:,:,:)
 integer, ABI_CONTIGUOUS pointer :: bids(:)
 real(dp) :: qp_gaps(self%ntemp),qpadb_gaps(self%ntemp)
 real(dp),allocatable :: aw(:,:,:)
 complex(dpc),target :: qpadb_enes(self%ntemp, self%max_nbcalc),qp_enes(self%ntemp, self%max_nbcalc)
 real(dp) :: ks_enes(self%max_nbcalc),ze0_vals(self%ntemp, self%max_nbcalc)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call xmpi_sum_master(self%vals_e0ks, master, comm, ierr)
 call xmpi_sum_master(self%dvals_de0ks, master, comm, ierr)
 call xmpi_sum_master(self%dw_vals, master, comm, ierr)
 if (self%nwr > 0) call xmpi_sum_master(self%vals_wr, master, comm, ierr)
 if (self%gfw_nomega > 0) call xmpi_sum_master(self%gfw_vals, master, comm, ierr)
 if (self%has_nuq_terms) call xmpi_sum_master(self%vals_nuq, master, comm, ierr)

 ! Only master writes
 if (my_rank /= master) return

 ik_ibz = self%kcalc2ibz(ikcalc, 1)

#if 0
 ! This to compute derivative of Re Sigma with 9 points.
 write(ab_out, *)"Using finite derivative"
 iw = 1 + (self%nwr / 2) - 4
 e0pde = arth(zero, self%wr_step, 9)
 do ibc=1,self%nbcalc_ks(ikcalc,spin)
   do it=1,self%ntemp
     smrt = linfit(9, e0pde, real(self%vals_wr(iw:iw+8, it, ibc)), alpha, beta)
     self%dvals_de0ks(it, ibc) = alpha
     if (smrt > 0.1 / Ha_eV) then
       write(msg,'(3a,i0,a,i0,2a,2(f22.15,2a))')&
         'Values of Re Sig_c are not linear ',ch10,&
         'band index ibc = ',ibc,' spin component = ',spin,ch10,&
         'root mean square= ',smrt,ch10,&
         'estimated slope = ',alpha,ch10,&
         'Omega [eV] SigC [eV]'
       MSG_WARNING(msg)
     end if
   end do
 end do
#endif

 if (self%symsigma == +1) then
   ! Average self-energy matrix elements in the degenerate subspace.
   do ideg=1,size(self%degtab(ikcalc, spin)%bids)
     bids => self%degtab(ikcalc, spin)%bids(ideg)%vals
     nstates = size(bids)

     do it=1,self%ntemp
       ! Average QP(T) and Z(T).
       cavg1 = sum(self%vals_e0ks(it, bids(:))) / nstates
       cavg2 = sum(self%dvals_de0ks(it, bids(:))) / nstates
       ravg = sum(self%dw_vals(it, bids(:))) / nstates
       do ii=1,nstates
         self%vals_e0ks(it, bids(ii)) = cavg1
         self%dvals_de0ks(it, bids(ii)) = cavg2
         self%dw_vals(it, bids(ii)) = ravg
       end do

       if (self%nwr > 0) then
         ! Average Sigma(omega, T)
         do iw=1,self%nwr
           cavg1 = sum(self%vals_wr(iw, it, bids(:))) / nstates
           do ii=1,nstates
             self%vals_wr(iw, it, bids(ii)) = cavg1
           end do
         end do
       end if
       !if (self%has_nuq_terms) then
       !  self%vals_nuq(it, bids(:), :, :, :) = sum(self%vals_nuq(it, bids(:), :, :, :)) / nstates
       !end if
     end do ! it
   end do ! ideg
 end if ! symsigma == +1

 ! Compute QP energies and Gaps.
 ib_val = nint(ebands%nelect / two); ib_cond = ib_val + 1
 kse_val = huge(one) * tol6; kse_cond = huge(one) * tol6
 qp_enes = huge(one) * tol6; qpadb_enes = huge(one) * tol6
 ks_enes = huge(one) * tol6; ze0_vals = huge(one) * tol6
 ks_gap = -one; qpadb_gaps = -one; qp_gaps = -one

 ! Write legend.
 if (ikcalc == 1 .and. spin == 1) then
   write(ab_out,"(a)")repeat("=", 80)
   write(ab_out,"(a)")"Final results in eV."
   write(ab_out,"(a)")"Notations:"
   write(ab_out,"(a)")"   eKS: Kohn-Sham energy. eQP: quasi-particle energy."
   write(ab_out,"(a)")"   eQP-eKS: Difference between the QP and the KS energy."
   write(ab_out,"(a)")"   SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part."
   write(ab_out,"(a)")"   Z(eKS): Renormalization factor."
   write(ab_out,"(a)")"   FAN: Real part of the Fan term at eKS. DW: Debye-Waller term."
   write(ab_out,"(a)")"   DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP."
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
 end if

 ! FIXME
 !do it=1,1
 do it=1,self%ntemp
   if (it == 1) then
     if (self%nsppol == 1) then
       write(ab_out,"(a)")sjoin("K-point:", ktoa(self%kcalc(:,ikcalc)))
     else
       write(ab_out,"(a)")sjoin("K-point:", ktoa(self%kcalc(:,ikcalc)), ", spin:", itoa(spin))
     end if
     write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
   end if

   do ibc=1,self%nbcalc_ks(ikcalc,spin)
     band = self%bstart_ks(ikcalc,spin) + ibc - 1
     kse = ebands%eig(band,ik_ibz,spin)
     ks_enes(ibc) = kse
     sig0c = self%vals_e0ks(it, ibc)
     dw = self%dw_vals(it, ibc)
     fan0 = real(sig0c) - dw
     ! TODO: Note that here I use the full Sigma including the imaginary part
     !zc = one / (one - self%dvals_de0ks(it, ibc))
     zc = one / (one - real(self%dvals_de0ks(it, ibc)))
     ze0_vals(it, ibc) = real(zc)
     qpe = kse + real(zc) * real(sig0c)
     qpe_adb = kse + real(sig0c)
     if (ibc == 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     if (band == ib_val) then
       kse_val = kse; qpe_val = qpe; qpe_adb_val = qpe_adb
     end if
     if (band == ib_cond) then
       kse_cond = kse; qpe_cond = qpe; qpe_adb_cond = qpe_adb
     end if
     if (it == 1) then
       !   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
       write(ab_out, "(i4,10(f8.3,1x))") &
         band, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
         real(sig0c) * Ha_eV, aimag(sig0c) * Ha_eV, real(zc), &
         fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
     end if
     if (ibc > 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     qpadb_enes(it, ibc) = qpe_adb
     qp_enes(it, ibc) = qpe
     if (kse_val /= huge(one) .and. kse_cond /= huge(one)) then
       ! We have enough states to compute the gap.
       if (it == 1) ks_gap = kse_cond - kse_val
       qpadb_gaps(it) = qpe_adb_cond - qpe_adb_val
       qp_gaps(it) = real(qpe_cond - qpe_val)
     end if
   end do ! ibc

   ! Print KS and QP gap
   if (it == 1 .and. kse_val /= huge(one) .and. kse_cond /= huge(one)) then
     write(ab_out, "(a)")" "
     write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, &
       "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
     write(ab_out, "(2(a,f8.3),a)")" QP gap: ",qp_gaps(it) * Ha_eV," (adiabatic: ",qpadb_gaps(it) * Ha_eV, ")"
     write(ab_out, "(2(a,f8.3),a)")" QP_gap - KS_gap: ",(qp_gaps(it) - ks_gap) * Ha_eV,&
         " (adiabatic: ",(qpadb_gaps(it) - ks_gap) * Ha_eV, ")"
     write(ab_out, "(a)")" "
   end if
 end do ! it

 ! Write data to netcdf file **Only master writes**
 ! (use iso_c_binding to associate a real pointer to complex data because netcdf does not support complex types).

#ifdef HAVE_NETCDF
 ! Write self-energy matrix elements for this (kpt, spin)
 ! Cannot use c_loc with gcc <= 4.8 due to internal compiler error so use c2r and stack memory.
 !shape3(1) = 2; shape4(1) = 2; shape5(1) = 2; shape6(1) = 2

 !shape3(2:) = shape(self%vals_e0ks); call c_f_pointer(c_loc(self%vals_e0ks), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "vals_e0ks"), c2r(self%vals_e0ks), start=[1,1,1,ikcalc,spin]))

 !shape3(2:) = shape(self%dvals_de0ks); call c_f_pointer(c_loc(self%dvals_de0ks), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "dvals_de0ks"), c2r(self%dvals_de0ks), start=[1,1,1,ikcalc,spin]))

 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "dw_vals"), self%dw_vals, start=[1,1,ikcalc,spin]))

 ! Dump QP energies and gaps for this (kpt, spin)
 !shape3(2:) = shape(qpadb_enes); call c_f_pointer(c_loc(qpadb_enes), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qpadb_enes"), c2r(qpadb_enes), start=[1,1,1,ikcalc,spin]))

 !shape3(2:) = shape(qp_enes); call c_f_pointer(c_loc(qp_enes), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qp_enes"), c2r(qp_enes), start=[1,1,1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "ze0_vals"), ze0_vals, start=[1,1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "ks_enes"), ks_enes, start=[1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "ks_gaps"), ks_gap, start=[ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qpadb_gaps"), qpadb_gaps, start=[1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qp_gaps"), qp_gaps, start=[1,ikcalc,spin]))

 ! Write frequency dependent data.
 if (self%nwr > 0) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "wrmesh_b"), self%wrmesh_b, start=[1,1,ikcalc,spin]))

   !shape4(2:) = shape(self%vals_wr); call c_f_pointer(c_loc(self%vals_wr), rdata4, shape4)
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "vals_wr"), c2r(self%vals_wr), start=[1,1,1,1,ikcalc,spin]))

   ! Compute spectral function.
   ! A = 1/pi [Im Sigma(ww)] / ([ww - ee - Re Sigma(ww)] ** 2 + Im Sigma(ww) ** 2])
   ABI_MALLOC(aw, (self%nwr, self%ntemp, self%max_nbcalc))

   do ib=1,self%nbcalc_ks(ikcalc,spin)
     band = self%bstart_ks(ikcalc, spin) + ib - 1
     kse = ebands%eig(band, ik_ibz, spin)
     do it=1,self%ntemp
       aw(:, it, ib) = piinv * abs(aimag(self%vals_wr(:, it, ib))) / &
         ((self%wrmesh_b(:, ib) - kse - real(self%vals_wr(:, it, ib))) ** 2 + aimag(self%vals_wr(:, it, ib)) ** 2)
     end do
   end do
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "spfunc_wr"), aw, start=[1,1,1,ikcalc,spin]))
   ABI_FREE(aw)
 end if

 ! Write Eliashberg functions
 if (self%gfw_nomega > 0) then
   !shape5(2:) = shape(self%gfw_vals); call c_f_pointer(c_loc(self%gfw_vals), rdata5, shape5)
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "gfw_vals"), c2r(self%gfw_vals), start=[1,1,1,1,1,ikcalc,spin]))
 end if
 if (self%has_nuq_terms) then
   !shape6(2:) = shape(self%vals_nuq); call c_f_pointer(c_loc(self%vals_nuq), rdata6, shape6)
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "vals_nuq"), c2r(self%vals_nuq), start=[1,1,1,1,1,1,ikcalc,spin]))
 end if
#endif

end subroutine sigmaph_gather_and_write
!!***

!!****f* m_sigmaph/sigmaph_print
!! NAME
!!  sigmaph_print
!!
!! FUNCTION
!!  Print self-energy and QP corrections for given (k-point, spin).
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  unt=Fortran unit number
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_print(self, dtset, unt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unt
 type(dataset_type),intent(in) :: dtset
 type(sigmaph_t),intent(in) :: self

!Local variables-------------------------------
!integer
 integer :: ikc,is

! *************************************************************************

 ! Write dimensions
 write(unt,"(a)")sjoin("Number of bands in e-ph self-energy:", itoa(self%nbsum))
 write(unt,"(a)")sjoin("Symsigma: ",itoa(self%symsigma), "Timrev:",itoa(self%timrev))
 write(unt,"(a)")sjoin("Imaginary shift in the denominator (zcut): ", ftoa(aimag(self%ieta) * Ha_eV, fmt="f5.3"), "[eV]")
 write(unt,"(a)")sjoin("Number of frequencies along the real axis:", itoa(self%nwr), ", Step:", ftoa(self%wr_step*Ha_eV), "[eV]")
 write(unt,"(a)")sjoin("Number of temperatures:", itoa(self%ntemp), &
   "From:", ftoa(self%kTmesh(1) / kb_HaK), "to", ftoa(self%kTmesh(self%ntemp) / kb_HaK), "[K]")
 write(unt,"(a)")sjoin("Ab-initio q-mesh from DDB file:", ltoa(dtset%ddb_ngqpt))
 write(unt,"(a)")sjoin("Q-mesh used for self-energy integration [ngqpt]:", ltoa(self%ngqpt))
 write(unt,"(a)")sjoin("Number of q-points in the IBZ:", itoa(self%nqibz))
 write(unt,"(a)")"List of K-points for self-energy corrections:"
 do ikc=1,self%nkcalc
   do is=1,self%nsppol
   if (self%nsppol == 2) write(unt,"(a,i1)")"... For spin: ",is
     write(unt, "(2(i4,2x),a,2(i4,1x))")&
       ikc, is, trim(ktoa(self%kcalc(:,ikc))), self%bstart_ks(ikc,is), self%bstart_ks(ikc,is) + self%nbcalc_ks(ikc,is) - 1
   end do
 end do

end subroutine sigmaph_print
!!***

end module m_sigmaph
!!***

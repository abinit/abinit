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
 use m_lgroup
 use m_ephwg
 use m_nctk
 use m_hdr
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven
 use m_io_tools,       only : iomode_from_fname
 use m_special_funcs,  only : dirac_delta
 use m_fftcore,        only : ngfft_seq
 use m_cgtools,        only : dotprod_g
 use m_cgtk,           only : cgtk_rotate
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, listkk
 use m_occ,            only : occ_fd, occ_be
 use m_double_grid,    only : double_grid_t
 use m_fftcore,        only : get_kg
 use m_kg,             only : getph
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 implicit none

 private
!!***

 public :: sigmaph   ! Main entry point to compute self-energy matrix elements
!!***

 ! Double grid datatype for electron-phonon
 type,public :: eph_double_grid_t
   type(ebands_t) :: ebands_dense

   real(dp),allocatable :: kpts_coarse(:,:)
   real(dp),allocatable :: kpts_dense(:,:)
   real(dp),allocatable :: kpts_dense_ibz(:,:)
   integer :: coarse_nbz, dense_nbz, dense_nibz

   real(dp),allocatable :: weights_dense(:)
   !weights in the dense grid

   integer,allocatable  :: bz2ibz_coarse(:)
   ! map full Brillouin zone to ibz (in the coarse grid)

   integer,allocatable  :: bz2ibz_dense(:)
   ! map full Brillouin zone to ibz (in the dense grid)

   integer :: nkpt_coarse(3), nkpt_dense(3)
   ! size of the coarse and dense meshes

   integer :: interp_kmult(3)
   ! multiplicity of the meshes

   integer :: ndiv
   ! interp_kmult(1)*interp_kmult(2)*interp_kmult(3)

   ! the integer indexes for the coarse grid are calculated for kk(3) with:
   ! [mod(nint((kpt(1)+1)*self%nkpt_coarse(1)),self%nkpt_coarse(1))+1,
   !  mod(nint((kpt(2)+1)*self%nkpt_coarse(2)),self%nkpt_coarse(2))+1,
   !  mod(nint((kpt(3)+1)*self%nkpt_coarse(3)),self%nkpt_coarse(3))+1]

   integer,allocatable :: indexes_to_coarse(:,:,:)
   ! given integer indexes get the array index of the kpoint (coarse)
   integer,allocatable :: coarse_to_indexes(:,:)
   ! given the array index get the integer indexes of the kpoints (coarse)

   integer,allocatable :: indexes_to_dense(:,:,:)
   ! given integer indexes get the array index of the kpoint (dense)
   integer,allocatable :: dense_to_indexes(:,:)
   ! given the array index get the integer indexes of the kpoint (dense)

   integer,allocatable :: coarse_to_dense(:,:)
   ! map coarse to dense mesh (nbz_coarse,mult(interp_kmult))

 end type eph_double_grid_t

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
   ! The spectral function is computed only if nwr > 0 (taken from dtser%nfreqsp)

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

  integer :: qint_method
   ! Defines the method used to integrate in q-space
   ! 0 --> Standard quadrature (one point per small box)
   ! 1 --> Use tetrahedron method

  logical :: use_doublegrid
   ! whether to use double grid or not

  type(eph_double_grid_t) :: eph_doublegrid
   ! store the double grid related object

  logical :: imag_only
   ! True if only the imaginary part of the self-energy must be computed

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
  ! Reduced coordinates of the q-points in the IBZ (full simmetry of the system).

  real(dp),allocatable :: wtq(:)
  ! wtq(nqibz)
  ! Weights of the q-points in the IBZ (normalized to one).

  real(dp),allocatable :: qibz_k(:,:)
  ! qibz(3,nqibz_k)
  ! Reduced coordinates of the q-points in the IBZ(k). Depends on ikcalc.

  integer,allocatable:: indkk(:, :)
   ! Mapping IBZ_k --> initial IBZ (self%lgrp%ibz --> self%ibz)

  real(dp),allocatable :: wtq_k(:)
  ! wtq(nqibz_k)
  ! Weights of the q-points in the IBZ(k) (normalized to one). Depends on ikcalc.

  real(dp),allocatable :: kTmesh(:)
   ! kTmesh(ntemp)
   ! List of temperatures (kT units).

  real(dp),allocatable :: mu_e(:)
   ! mu_e(ntemp)
   ! chemical potential of electrons for the different temperatures.

  real(dp),allocatable :: e0vals(:)
  ! (nbcalc_ks)
  ! KS energies where QP corrections are wantend
  ! This array is initialized inside the (ikcalc, spin) loop

  complex(dpc),allocatable :: cweights(:,:,:,:,:,:)
  ! (nz, 2, nbcalc_ks, natom3, nq, ibsum))
  ! Weights for the q-integration of 1 / (e1 - e2 \pm w_{q, nu} + i.eta)
  ! This array is initialized inside the (ikcalc, spin) loop

  real(dp),allocatable :: deltaw_pm(:,:,:,:,:)
  ! (nbcalc_ks, 2, natom3, nq, ibsum))
  ! Weights for the q-integration of the two delta (abs/emission) if imag_only
  ! This array is initialized inside the (ikcalc, spin) loop

  real(dp),allocatable :: wrmesh_b(:,:)
  ! wrmesh_b(nwr, max_nbcalc)
  ! Frequency mesh along the real axis (Ha units) used for the different bands
  ! This array depends on (ikcalc, spin)

  complex(dpc),allocatable :: vals_e0ks(:,:)
   ! vals_e0ks(ntemp, max_nbcalc)
   ! Sigma_eph(omega=eKS, kT, band) for fixed (kcalc, spin).

  complex(dpc),allocatable :: dvals_de0ks(:,:)
   ! dvals_de0ks(ntemp, max_nbcalc) for fixed (kcalc, spin)
   ! d Re Sigma_eph(omega, kT, band, kcalc, spin) / d omega (omega=eKS)

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
   ! Number of frequencies in Eliashberg function.
   ! Set to 0 to deactivate this part.
   ! TODO: Integration in q-space is done according to eph_intmeth.

  real(dp),allocatable :: gfw_mesh(:)
   ! Frequency mesh for Eliashberg function (phonons)
   ! gfw_mesh(gfw_nomega)

   real(dp),allocatable :: gf_nnuq(:,:,:,:)
   ! (nbcalc_ks, natom3, %nqibz_k, 3)
   ! Quantities needed to compute generalized Eliashberg functions  (gkk2/Fan-Migdal/DW terms)
   ! This array depends on (ikcalc, spin)

  ! TODO: Can be removed now.
  real(dp),allocatable :: gfw_vals(:,:,:)
   !gfw_vals(gfw_nomega, 3, max_nbcalc)
   ! Generalized Eliashberg function a2F_{n,k,spin}(w)
   ! 1:   gkk^2 with cutoff on energy difference (en - em)
   ! 2:3 (Fan-Migdal/DW contribution)
   ! This array depends on (ikcalc, spin)

  type(ephwg_t) :: ephwg
   ! This object compute the weights for the BZ integration in q-space if qint_method > 0

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
 integer,parameter :: dummy_npw=1,tim_getgh1c=1,berryopt0=0,bcorr0=0,timrev0=0
 integer,parameter :: useylmgr=0,useylmgr1=0,master=0,ndat1=1,nz=1
 integer,parameter :: sppoldbl1=1,timrev1=1
 integer :: my_rank,mband,my_minb,my_maxb,nsppol,nkpt,iq_ibz
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,nprocs
 integer :: ibsum_kq,ib_k,band,num_smallw,ibsum,ii,jj,im,in !,ib,nstates
 integer :: isym, this_calc
 integer :: idir,ipert,ip1,ip2,idir1,ipert1,idir2,ipert2
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq !,!timerev_q,
 integer :: ik_ibz_fine,iq_ibz_fine,ikq_ibz_fine,ik_bz_fine,ikq_bz_fine,iq_bz_fine
 integer :: ik_bz, ikq_bz, iq_bz
 integer :: spin,istwf_k,istwf_kq,istwf_kqirr,npw_k,npw_kq,npw_kqirr
 integer :: mpw,ierr,it !ipw
 integer :: n1,n2,n3,n4,n5,n6,nspden,do_ftv1q,nu,ndiv
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnl1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,nq
 integer :: nbcalc_ks,nbsum,bstart_ks,ikcalc
 real(dp),parameter :: tol_enediff=0.001_dp*eV_Ha
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,cpu_dvscf,wall_dvscf,gflops_all
 real(dp) :: ecut,eshift,dotr,doti,dksqmax,weigth_q,rfact,alpha,beta,gmod2,hmod2,ediff,weight
 complex(dpc) :: cfact,dka,dkap,dkpa,dkpap,cplx_ediff
 logical,parameter :: have_ktimerev=.True.
 logical :: isirr_k,isirr_kq,gen_eigenpb,isqzero
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(sigmaph_t) :: sigma
 type(eph_double_grid_t) :: eph_dg
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0_kq(3),dummy_gvec(3,dummy_npw)
 integer :: work_ngfft(18),gmax(3) !!g0ibz_kq(3),
 integer :: indkk_kq(1,6)
 integer,allocatable :: gtmp(:,:),kg_k(:,:),kg_kq(:,:),nband(:,:),distrib_bq(:,:) !,degblock(:,:),
 integer,allocatable :: eph_dg_mapping(:,:), iqlk(:), indkk(:)
 integer,allocatable :: indq2dvdb(:,:),wfd_istwfk(:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),phfrq(3*cryst%natom),sqrt_phfrq0(3*cryst%natom)
 real(dp) :: kq_sym(3)
 real(dp) :: wqnu,nqnu,gkk2,eig0nk,eig0mk,eig0mkq,f_mkq,eminmax(2)
 !real(dp) :: wqnu,nqnu,gkk2,eig0nk,eig0mk,eig0mkq,ediff,f_mkq !,f_nk
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom),displ_red(2,3,cryst%natom,3*cryst%natom)
 !real(dp) :: ucart(2,3,cryst%natom,3*cryst%natom)
 real(dp) :: d0mat(2,3*cryst%natom,3*cryst%natom)
 complex(dpc) :: cmat(3*cryst%natom,3*cryst%natom),cvec1(3*cryst%natom),cvec2(3*cryst%natom)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkk_atm(:,:,:),gkk_nu(:,:,:),dbwl_nu(:,:,:,:),gdw2_mn(:,:),gkk0_atm(:,:,:,:)
 complex(dpc),allocatable :: tpp(:,:),tpp_red(:,:),hka_mn(:,:,:),wmat1(:,:,:),zvals(:,:)
 real(dp),allocatable :: bra_kq(:,:),kets_k(:,:,:),h1kets_kq(:,:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnl1(:,:),work(:,:,:,:)
 real(dp),allocatable ::  gs1c(:,:),nqnu_tlist(:),dt_weights(:,:),dargs(:)
 real(dp),allocatable :: tmp_deltaw_pm(:,:,:)
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
 cpu_dvscf = 0; wall_dvscf = 0

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
 if (my_rank == master) then
   call sigmaph_print(sigma, dtset, ab_out)
   call sigmaph_print(sigma, dtset, std_out)
 end if

 ! This is the maximum number of PWs for all possible k+q treated.
 mpw = sigma%mpw; gmax = sigma%gmax

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
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

 ! Distribute bands
 !if (sigma%imag_only) then
 !  bks_mask = .False.
 !  bks_mask(bstart:bstop, : ,:) = .True.
 !else
 !  bks_mask = .False.
 !  call xmpi_split_work(sigma%nbsum, comm, ip1, ip2, msg, ierr)
 !  if (ip1 == sigma%nbsum + 1) then
 !    MSG_ERROR("sigmaph with idle processes should be tested!")
 !  end if
 !  bks_mask(ip1:ip2, : ,:) = .True.
 !end if

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,dtset%paral_kgb,dummy_npw,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands%kptns,ngfft,&
   dummy_gvec,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm,opt_ecut=ecut)

 call wfd_print(wfd,header="Wavefunctions for self-energy calculation.",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 call wfd_read_wfk(wfd, wfk0_path, iomode_from_fname(wfk0_path))
 if (.False.) call wfd_test_ortho(wfd, cryst, pawtab, unit=std_out, mode_paral="PERS")

!ENCAPSULATE
 ! Double grid stuff
 ndiv = 1
 if (sigma%use_doublegrid) then
   eph_dg       = sigma%eph_doublegrid
   ABI_MALLOC(eph_dg_mapping,(6,eph_dg%ndiv))
   ABI_MALLOC(indkk,(eph_dg%dense_nbz))
   ndiv = eph_dg%ndiv
 end if

 ! Mapping of q point to IBZ
 ABI_MALLOC(iqlk,(ndiv))
!ENCAPSULATE

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
 ! Set cache in Mb for q-points.
 ! When we ask for a qpt in the BZ the code checks whether the IBZ image is in the cache.
 ! and if we have a cache hit we "rotate" the qibz in the cache to return V(qbz) else the code reads qibz and rotates.
 ! When we need to remove a cache item because we've reached the dvdb_qcache_mb limit, we select the q-point
 ! which the largest number of operations in the little group (e.g. Gamma) while trying to keep the previous qibz in cache
 ! The cache is built dynamically so it depends on the way we loop over q-points in the caller.
 ! This is also the reason why we reorder the q-points in ibz_k to pack the points in *shells*  to minimise cache misses.
 call dvdb_set_qcache_mb(dvdb, dtset%dvdb_qcache_mb)
 call dvdb_print(dvdb, prtvol=dtset%prtvol)

 ! Loop over k-points in Sigma_nk.
 do ikcalc=1,sigma%nkcalc
   kk = sigma%kcalc(:, ikcalc)

   ! Find IBZ(k) for q-point integration.
   call sigmaph_setup_kcalc(sigma, cryst, ikcalc, dtset%prtvol)
   call wrtout(std_out, sjoin(ch10, repeat("=", 92)))
   call wrtout(std_out, sjoin("Computing self-energy matrix elements for k-point:", ktoa(kk)))
   call wrtout(std_out, sjoin("Number of q-points in the IBZ(k):", itoa(sigma%nqibz_k)))

   ! Symmetry indices for kk.
   ik_ibz = sigma%kcalc2ibz(ikcalc,1); isym_k = sigma%kcalc2ibz(ikcalc,2)
   trev_k = sigma%kcalc2ibz(ikcalc,6); g0_k = sigma%kcalc2ibz(ikcalc,3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   ABI_CHECK(isirr_k, "For the time being the k-point must be in the IBZ")
   kk_ibz = ebands%kptns(:,ik_ibz)
   npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)

   ! Find correspondence IBZ_k --> set of q-points in DVDB.
   ! Need to handle q_bz = S q_ibz by symmetrizing the potentials already available in the DVDB.
   ! Activate Fourier interpolation only if q-points cannot be reconstructed from the DVDB file.
   !
   ! Note:
   !   * q --> -q symmetry is always used for phonons.
   !   * we use symrec instead of symrel (see also m_dvdb)
   !
   ABI_MALLOC(indq2dvdb, (6, sigma%nqibz_k))
   do_ftv1q = 0
   do iq_ibz=1,sigma%nqibz_k
     call listkk(dksqmax, cryst%gmet, indkk_kq, dvdb%qpts, sigma%qibz_k(:,iq_ibz), dvdb%nqpt, 1, cryst%nsym, &
        1, cryst%symafm, cryst%symrec, timrev1, use_symrec=.True.)
     indq2dvdb(:, iq_ibz) = indkk_kq(1, :)
     if (dksqmax > tol12) then
       ! Cannot recostruct this qpt by symmetry. Set entry to -1 and activate Fourier interpolation.
       indq2dvdb(:, iq_ibz) = -1
       do_ftv1q = do_ftv1q + 1
     end if
   end do

   if (do_ftv1q /= 0) then
     call cwtime(cpu,wall,gflops,"start")
     write(msg, "(2(a,i0),a)")"Will use Fourier interpolation of DFPT potentials [",do_ftv1q,"/",sigma%nqibz_k,"]"
     call wrtout(std_out, msg)
     call wrtout(std_out, sjoin("From ngqpt", ltoa(ifc%ngqpt), "to", ltoa(sigma%ngqpt)))
     call dvdb_ftinterp_setup(dvdb, ifc%ngqpt, 1, [zero,zero,zero], nfftf, ngfftf, comm)
     call cwtime(cpu,wall,gflops,"stop")
     wall_dvscf = wall_dvscf + wall
     cpu_dvscf = cpu_dvscf + cpu
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

     ! Prepare spectral function and/or Eliasberg function.
     if (sigma%nwr > 0) then
       sigma%vals_wr = zero
       do ib_k=1,nbcalc_ks
         band = ib_k + bstart_ks - 1
         eig0nk = ebands%eig(band,ik_ibz,spin) - sigma%wr_step * (sigma%nwr / 2)
         sigma%wrmesh_b(:,ib_k) = arth(eig0nk, sigma%wr_step, sigma%nwr)
       end do
     end if
     if (sigma%gfw_nomega > 0) then
       if (allocated(sigma%gf_nnuq)) then
         ABI_FREE(sigma%gf_nnuq)
       end if
       ABI_CALLOC(sigma%gf_nnuq, (nbcalc_ks, natom3, sigma%nqibz_k, 3))
     end if

     ! Allocate eph matrix elements.
     ABI_MALLOC(gkk_atm, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkk_nu, (2, nbcalc_ks, natom3))

     ! Arrays for Debye-Waller
     if (.not. sigma%imag_only) then
       ABI_STAT_MALLOC(dbwl_nu, (2, nbcalc_ks, nbsum, natom3), ierr)
       ABI_CHECK(ierr == 0, "oom in dbwl_nu")
       dbwl_nu = zero
       ABI_MALLOC(gkk0_atm, (2, nbcalc_ks, nbsum, natom3))
       ABI_CHECK(ierr == 0, "oom in gkk0_atm")
       gkk0_atm = zero
     end if

     ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
     ! and save KS energies in sigma%e0vals
     ! TODO: symmetrize them if kk is not irred
     ABI_MALLOC(kets_k, (2, npw_k*nspinor, nbcalc_ks))
     ABI_MALLOC(sigma%e0vals, (nbcalc_ks))
     do ib_k=1,nbcalc_ks
       band = ib_k + bstart_ks - 1
       call wfd_copy_cg(wfd, band, ik_ibz, spin, kets_k(1,1,ib_k))
       sigma%e0vals(ib_k) = ebands%eig(band, ik_ibz, spin)
     end do

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Distribute q-points and bands.
     ABI_MALLOC(distrib_bq, (nbsum, sigma%nqibz_k))
     !if (sigma%nqibz_k >= nprocs .and. mod(sigma%nqibz_k, nprocs) == 0) then
     if (sigma%nqibz_k >= nprocs) then
       do iq_ibz=1,sigma%nqibz_k
         distrib_bq(:, iq_ibz) = mod(iq_ibz, nprocs)
       end do
     else
       do it=1,nbsum*sigma%nqibz_k
         ibsum_kq = mod(it-1, nbsum) + 1; iq_ibz = (it - ibsum_kq) / nbsum + 1
         distrib_bq(ibsum_kq, iq_ibz) = mod(it, nprocs)
       end do
     end if

!ENCAPSULATE
     if (sigma%qint_method > 0) then
       ! Weights for Re-Im with i.eta shift.
       ABI_MALLOC(sigma%cweights, (nz, 2, nbcalc_ks, natom3, nbsum, sigma%ephwg%nq_k))
       ! Weights for Im (tethraedron, eta --> 0)
       ABI_CALLOC(sigma%deltaw_pm, (2 ,nbcalc_ks, natom3, nbsum, sigma%ephwg%nq_k))
     endif
     ABI_MALLOC(zvals, (nz, nbcalc_ks))

     if (sigma%qint_method > 0) then
       ! Map eph_dg%dense -> ephwg%lgk%ibz
       if (sigma%use_doublegrid) then
         call cwtime(cpu,wall,gflops,"start")
         indkk = 0
         do ikq_ibz=1,sigma%ephwg%lgk%nibz
           ! get coordinates of ephwg%lgk%ibz point
           kq(:) = sigma%ephwg%lgk%ibz(:,ikq_ibz)
           ! Loop over the star of q
           do isym=1,sigma%ephwg%lgk%nsym_lg
             ! Get the symmetric of q
             do ii=1,3
               kq_sym(ii)=kq(1)*sigma%ephwg%lgk%symrec_lg(ii,1,isym)&
                         +kq(2)*sigma%ephwg%lgk%symrec_lg(ii,2,isym)&
                         +kq(3)*sigma%ephwg%lgk%symrec_lg(ii,3,isym)
             end do
             ! get the index of the ibz point in bz
             ikq_bz = eph_double_grid_get_index(eph_dg,kq_sym,2)
             indkk(ikq_bz) = ikq_ibz
           end do
         end do
         call cwtime(cpu,wall,gflops,"stop")
         write(msg,'(2(a,f8.2))') "little group of k mapping cpu:",cpu,", wall:",wall
         call wrtout(std_out, msg, do_flush=.True.)
       endif

       ! Precompute the weights for tetrahedron
       call cwtime(cpu,wall,gflops,"start")
       ABI_MALLOC(tmp_deltaw_pm,(3,sigma%ephwg%nq_k, 2))
       ! loop over bands to sum
       do ibsum_kq=1,nbsum
         ! loop over phonon modes
         do nu=1,natom3
           ! loop over bands in the self-energy
           do ib_k=1,nbcalc_ks
             this_calc = (ibsum_kq-1)*natom3*nbcalc_ks + (nu-1)*nbcalc_ks + ib_k
             if (mod(this_calc,nprocs) /= my_rank) cycle
             band = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band,ik_ibz,spin)
             eminmax(1) = eig0nk - 0.01 
             eminmax(2) = eig0nk + 0.01
             call ephwg_get_deltas(sigma%ephwg, ibsum_kq, spin, nu, 3, eminmax, bcorr0, tmp_deltaw_pm, xmpi_comm_self)
             ! we pay the efficiency here
             sigma%deltaw_pm(1,ib_k,nu,ibsum_kq,:) = tmp_deltaw_pm(2, :, 1) / ( sigma%ephwg%lgk%weights(:) )
             sigma%deltaw_pm(2,ib_k,nu,ibsum_kq,:) = tmp_deltaw_pm(2, :, 2) / ( sigma%ephwg%lgk%weights(:) )
           enddo
         enddo
       enddo
       ABI_FREE(tmp_deltaw_pm)
       call xmpi_sum(sigma%deltaw_pm, comm, ierr)
       call cwtime(cpu,wall,gflops,"stop")
       write(msg,'(2(a,f8.2))') "weights with tetrahedron  cpu:",cpu,", wall:",wall
       call wrtout(std_out, msg, do_flush=.True.)
     endif
!ENCAPSULATE

     ! Integrations over q-points in the IBZ(k)
     do iq_ibz=1,sigma%nqibz_k
       ! Quick-parallelization over q-points
       if (all(distrib_bq(1:nbsum, iq_ibz) /= my_rank)) cycle

       qpt = sigma%qibz_k(:,iq_ibz)
       isqzero = (sum(qpt**2) < tol14) !; if (isqzero) cycle
       call cwtime(cpu,wall,gflops,"start")

       db_iqpt = indq2dvdb(1, iq_ibz)

       if (db_iqpt /= -1) then
         if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found:", ktoa(qpt), "in DVDB with index", itoa(db_iqpt)))
         ! Read and reconstruct the dvscf potentials for qpt and all 3*natom perturbations.
         ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
         call dvdb_readsym_qbz(dvdb, cryst, qpt, indq2dvdb(:,iq_ibz), cplex, nfftf, ngfftf, v1scf, xmpi_comm_self)
       else
         ! Fourier interpolation of the potential
         if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Could not find:", ktoa(qpt), "in DVDB - interpolating"))
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
       ! DEBUG
       !call ifc_fourq(ifc, cryst, sigma%qibz(:, sigma%indkk(iq_ibz, 1)), &
       !   phfrq_ibz, displ_cart, out_displ_red=displ_red)
       call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)
       !phfrq = phfrq_ibz

       ! Find k+q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qpt
       call listkk(dksqmax,cryst%gmet,indkk_kq,ebands%kptns,kq,ebands%nkpt,1,cryst%nsym,&
         sppoldbl1,cryst%symafm,cryst%symrel,sigma%timrev,use_symrec=.False.)

       if (dksqmax > tol12) then
         write(msg, '(4a,es16.6,7a)' )&
          "The WFK file cannot be used to compute self-energy corrections at k:.", trim(ktoa(kk)), ch10,&
          "At least one of the k+q points could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
          "Q-mesh: ",trim(ltoa(sigma%ngqpt)),", K-mesh (from kptrlatt) ",trim(ltoa(get_diag(dtset%kptrlatt))),ch10, &
          'Action: check your WFK file and (k,q) point input variables'
         MSG_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1,1); isym_kq = indkk_kq(1,2)
       trev_kq = indkk_kq(1, 6); g0_kq = indkk_kq(1, 3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0)) !; isirr_kq = .True.
       kq_ibz = ebands%kptns(:,ikq_ibz)

!ENCAPSULATE
       ! Double grid stuff
       if (sigma%use_doublegrid) then
         ik_bz  = eph_double_grid_get_index(eph_dg,kk,1)
         ikq_bz = eph_double_grid_get_index(eph_dg,kq,1)
         iq_bz  = eph_double_grid_get_index(eph_dg,qpt,1)

         ik_bz_fine  = eph_dg%coarse_to_dense(ik_bz,1)
         ik_ibz_fine = eph_dg%bz2ibz_dense(ik_bz_fine)
         !fine grid around kq
         do jj=1,eph_dg%ndiv

           !kq
           ikq_bz_fine = eph_dg%coarse_to_dense(ikq_bz,jj)
           ikq_ibz_fine = eph_dg%bz2ibz_dense(ikq_bz_fine)

           !qq
           iq_bz_fine = eph_dg%coarse_to_dense(iq_bz,jj)
           iq_ibz_fine = eph_dg%bz2ibz_dense(iq_bz_fine)

           eph_dg_mapping(:, jj) = &
             [ik_bz_fine,  ikq_bz_fine,  iq_bz_fine,&
              ik_ibz_fine, ikq_ibz_fine, iq_ibz_fine]
         enddo
       endif
 
       ! Map q to qibz for tetrahedron
       if (sigma%qint_method > 0) then
         if (sigma%use_doublegrid) then
           ! set the qpoints to be mapped
           do jj=1,eph_dg%ndiv
               iq_bz_fine = eph_dg_mapping(3,jj)
               iqlk(jj) = indkk(iq_bz_fine)
           end do
         else
           iqlk(1) = iq_ibz
           if (sigma%symsigma == 0) iqlk(1) = lgroup_find_ibzimage(sigma%ephwg%lgk, qpt)
           ABI_CHECK(iqlk(1) /= -1, sjoin("Cannot find q-point in IBZ(k)", ktoa(qpt)))
           if (sigma%symsigma == 1) then
             ABI_CHECK(all(abs(sigma%qibz_k(:, iq_ibz) - sigma%ephwg%lgk%ibz(:, iqlk(1))) < tol12), "Mismatch in qpoints.")
           end if
         endif
       end if
!ENCAPSULATE
       
       ! Get npw_kq, kg_kq for k+q
       ! Be careful with time-reversal symmetry and istwf_kq
       if (isirr_kq) then
         ! Copy u_kq(G)
         istwf_kq = wfd%istwfk(ikq_ibz); npw_kq = wfd%npwarr(ikq_ibz)
         ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
         kg_kq(:,1:npw_kq) = wfd%kdata(ikq_ibz)%kg_k
       else
         ! Will Reconstruct u_kq(G) from the IBZ image.
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
         if (isqzero .and. .not. sigma%imag_only) then
           gkk0_atm(:, :, ibsum_kq, :) = gkk_atm
           dbwl_nu(:, :, ibsum_kq, :) = gkk_nu
           !dbwl_nu(:, :, ibsum_kq, :) = one
         end if

         ! Accumulate contribution to self-energy
         eig0mkq = ebands%eig(ibsum_kq,ikq_ibz,spin)
         ! q-weigths for naive integration
         weigth_q = sigma%wtq_k(iq_ibz)

         do nu=1,natom3
           ! Ignore acoustic or unstable modes.
           wqnu = phfrq(nu); if (wqnu < tol6) cycle

           ! For each band in Sigma_{bk}
           do ib_k=1,nbcalc_ks
             band = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band,ik_ibz,spin)
             gkk2 = weigth_q * (gkk_nu(1,ib_k,nu) ** 2 + gkk_nu(2,ib_k,nu) ** 2)

             ! Optionally, accumulate contribution to Eliashberg functions
             ediff = eig0nk - eig0mkq
             !if (abs(cplx_ediff) < tol6) cplx_ediff = cplx_ediff + sigma%ieta
             if (sigma%gfw_nomega > 0 .and. abs(eig0nk - eig0mkq) > tol6) then
               ! eph strength with "energy" cutoff.
               if ((abs(ediff) > wqnu -  three * dtset%ph_smear) .and. (abs(ediff) < wqnu + three * dtset%ph_smear)) then
                 sigma%gf_nnuq(ib_k, nu, iq_ibz, 1) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 1) + &
                    (gkk_nu(1,ib_k,nu) ** 2 + gkk_nu(2,ib_k,nu) ** 2)
               end if
               ! Fan term.
               if (ediff > wqnu) then
                  rfact = one / ediff
               else
                 ! Non adiabatic regime --> Add complex shift.
                 ! Note however that the expression for the Eliashberg function relies on adiabaticity.
                 rfact = real(one/(ediff + sigma%ieta))
               end if
               sigma%gf_nnuq(ib_k, nu, iq_ibz, 2) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 2) + &
                  (gkk_nu(1,ib_k,nu) ** 2 + gkk_nu(2,ib_k,nu) ** 2) * rfact
             end if

             do it=1,sigma%ntemp
               nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)
               ! TODO: In principle mu_e(T) (important if semimetal) but need to treat T --> 0 in new_occ
               f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

!ENCAPSULATE
               ! Here we have to handle 3 different logical values
               ! leading to 9 different cases:
               !
               ! qint_method         0      1
               !   use_doublegrid   .true. .false.
               !     imag_only      .true. .false.
               !
               ! we will write this with nested conditionals using the order above

               if (sigma%qint_method == 0) then
                 if (sigma%use_doublegrid) then
                   cfact = 0
                   do jj=1,eph_dg%ndiv
                     ! Double Grid shared points weights
                     ikq_bz_fine  = eph_dg_mapping(2, jj)
                     weight = eph_dg%weights_dense(ikq_bz_fine)

                     ! Electronic eigenvalue
                     ikq_ibz_fine = eph_dg_mapping(5, jj)
                     eig0mkq = eph_dg%ebands_dense%eig(ibsum_kq,ikq_ibz_fine,spin)
                     f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

                     ! Phonon frequency
                     iq_ibz_fine = eph_dg_mapping(6, jj)
                     wqnu = sigma%ephwg%phfrq_ibz(iq_ibz_fine,nu)
                     !if (wqnu < tol6) cycle
                     nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)

                     cfact = cfact + &
                            ((nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                             (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta) )*weight
                   enddo
                 else
                   cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                            (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)
                 endif
                 if (sigma%imag_only) then
                   sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * j_dpc * aimag(cfact)
                 else
                   sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * cfact
                 end if
               else
                 if (sigma%use_doublegrid) then
                   do jj=1,eph_dg%ndiv
                     ! Double Grid shared points weights
                     ikq_bz_fine  = eph_dg_mapping(2, jj)
                     weight = eph_dg%weights_dense(ikq_bz_fine)

                     ! Electronic eigenvalue
                     ikq_ibz_fine = eph_dg_mapping(5, jj)
                     eig0mkq = eph_dg%ebands_dense%eig(ibsum_kq,ikq_ibz_fine,spin)
                     f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

                     ! Phonon frequency
                     iq_ibz_fine = eph_dg_mapping(6, jj)
                     wqnu = sigma%ephwg%phfrq_ibz(iq_ibz_fine,nu)
                     !if (wqnu < tol6) cycle
                     nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)

                     if (sigma%imag_only) then
                       sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * j_dpc * pi * ( &
                         (nqnu + f_mkq      ) * sigma%deltaw_pm(1, ib_k, nu, ibsum_kq, iqlk(jj)) +  &
                         (nqnu - f_mkq + one) * sigma%deltaw_pm(2, ib_k, nu, ibsum_kq, iqlk(jj)) )*weight
                     else
                       sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * ( &
                         (nqnu + f_mkq      ) * sigma%cweights(1, 1, ib_k, nu, ibsum_kq, ikq_ibz_fine) +  &
                         (nqnu - f_mkq + one) * sigma%cweights(1, 2, ib_k, nu, ibsum_kq, ikq_ibz_fine) )*weight
                     end if
                   end do
                 else
                   if (sigma%imag_only) then
                     sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * j_dpc * pi * ( &
                       (nqnu + f_mkq      ) * sigma%deltaw_pm(1, ib_k, nu, ibsum_kq, iqlk(1)) +  &
                       (nqnu - f_mkq + one) * sigma%deltaw_pm(2, ib_k, nu, ibsum_kq, iqlk(1)) )
                   else
                     sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * ( &
                       (nqnu + f_mkq      ) * sigma%cweights(1, 1, ib_k, nu, ibsum_kq, ikq_ibz_fine) +  &
                       (nqnu - f_mkq + one) * sigma%cweights(1, 2, ib_k, nu, ibsum_kq, ikq_ibz_fine) )
                   endif
                 end if
               end if
!ENCAPSULATE

               ! Derivative of sigma
               ! TODO: should calculate this with the double grid as well
               if (.not. sigma%imag_only) then
                 !if (sigma%qint_method == 1) then
                   ! Have to rescale gkk2 before computing derivative (HM: why?)
                 !  gkk2 = gkk2 * sigma%wtq_k(iq_ibz)
                 !end if

                 ! Accumulate d(Re Sigma) / dw(w=eKS) for state ib_k
                 !cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
                 !            (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
                 gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
                 hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
                 rfact = (nqnu + f_mkq      ) * (-gmod2 + aimag(sigma%ieta)**2) / (gmod2 + aimag(sigma%ieta)**2) ** 2 + &
                         (nqnu - f_mkq + one) * (-hmod2 + aimag(sigma%ieta)**2) / (hmod2 + aimag(sigma%ieta)**2) ** 2
                 sigma%dvals_de0ks(it, ib_k) = sigma%dvals_de0ks(it, ib_k) + gkk2 * rfact
                 !cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                 !         (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)
                 !sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkk2 * cfact

                 !cfact = (eig0nk - eig0mkq + wqnu + sigma%ieta)
                 !gmod2 = cfact * dconjg(cfact)
                 !cfact = (eig0nk - eig0mkq - wqnu + sigma%ieta)
                 !hmod2 = cfact * dconjg(cfact)
                 !sigma%dvals_de0ks(it, ib_k) = sigma%dvals_de0ks(it, ib_k) + gkk2 * ( &
                 !  (nqnu + f_mkq)        * (gmod2 - two * (eig0nk - eig0mkq + wqnu) ** 2) / gmod2 ** 2 + &
                 !  (nqnu - f_mkq + one)  * (hmod2 - two * (eig0nk - eig0mkq - wqnu) ** 2) / hmod2 ** 2   &
                 !)

                 ! Accumulate Sigma(w) for state ib_k if spectral function is wanted.
                 ! TODO: weigths
                 if (sigma%nwr > 0) then
                   cfact_wr(:) = (nqnu + f_mkq      ) / (sigma%wrmesh_b(:,ib_k) - eig0mkq + wqnu + sigma%ieta) + &
                                 (nqnu - f_mkq + one) / (sigma%wrmesh_b(:,ib_k) - eig0mkq - wqnu + sigma%ieta)
                   sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + gkk2 * cfact_wr(:)
                 end if

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
       write(msg,'(4(a,i0),2(a,f8.2))') "k-point [",ikcalc,"/",sigma%nkcalc, &
                                        "] q-point [",iq_ibz,"/",sigma%nqibz_k,"] completed. cpu:",cpu,", wall:",wall
       call wrtout(std_out, msg, do_flush=.True.)
     end do ! iq_ibz (sum over q-points)

     if (sigma%qint_method > 0) then
       ABI_FREE(sigma%deltaw_pm)
       ABI_FREE(sigma%cweights)
     end if
     ! Print cache stats. The first k-point is expected to have lots of misses
     ! especially if it's the Gamma point with symsigma = 1.
     if (dvdb%qcache_size > 0) then
       write(std_out, "(a)")"Qcache stats"
       write(std_out, "(a,i0)")"Total Number of calls: ", dvdb%qcache_stats(1)
       write(std_out, "(a,i6,f6.1,a)")&
         "Cache hit:  ", dvdb%qcache_stats(2), (100.0_dp * dvdb%qcache_stats(2)) / dvdb%qcache_stats(1), "%"
       write(std_out, "(a,i6,f6.1,a)")&
         "Cache miss: ", dvdb%qcache_stats(3), (100.0_dp * dvdb%qcache_stats(3)) / dvdb%qcache_stats(1), "%"
       dvdb%qcache_stats = 0
     end if

     ABI_FREE(sigma%e0vals)
     ABI_FREE(zvals)
     ABI_FREE(kets_k)
     ABI_FREE(gkk_atm)
     ABI_FREE(gkk_nu)
     ABI_FREE(distrib_bq)

     ! =========================
     ! Compute Debye-Waller term
     ! =========================
     if (.not. sigma%imag_only) then
       call xmpi_sum(dbwl_nu, comm, ierr)
       call xmpi_sum(gkk0_atm, comm, ierr)

       ABI_MALLOC(gdw2_mn, (nbsum, nbcalc_ks))
       ABI_MALLOC(tpp, (natom3, natom3))
       ABI_MALLOC(tpp_red, (natom3, natom3))
       ABI_MALLOC(hka_mn, (natom3, nbsum, nbcalc_ks))

       qpt = zero
       call ifc_fourq(ifc, cryst, qpt, sqrt_phfrq0, displ_cart)
       where (sqrt_phfrq0 >= zero)
         sqrt_phfrq0 = sqrt(sqrt_phfrq0)
       else where
         sqrt_phfrq0 = zero
       end where
       ! cmat contains the displament vectors as complex array
       d0mat = reshape(displ_cart, [2, natom3, natom3])
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
       ! TODO Bug if symsigma 0 ?
       nq = sigma%nqibz; if (sigma%symsigma == 0) nq = sigma%nqbz
       if (sigma%symsigma == +1) nq = sigma%nqibz_k

       do iq_ibz=1,nq
         if (mod(iq_ibz, nprocs) /= my_rank) cycle  ! MPI parallelism

         if (abs(sigma%symsigma) == 1) then
           !qpt = sigma%qibz(:,iq_ibz); weigth_q = sigma%wtq(iq_ibz)
           qpt = sigma%qibz_k(:,iq_ibz); weigth_q = sigma%wtq_k(iq_ibz)
         else
           qpt = sigma%qbz(:,iq_ibz); weigth_q = one / sigma%nqbz
         end if

         ! Get phonons for this q-point.
         ! DEBUG
         !call ifc_fourq(ifc, cryst, sigma%qibz(:, sigma%indkk(iq_ibz, 1)), phfrq_ibz, displ_cart)
         call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)
         !phfrq = phfrq_ibz
         !if (all(abs(qpt) < tol12)) cycle

         ! Compute hka_mn matrix with shape: (natom3, nbsum, nbcalc_ks))
         ! Needed for Giustino's equation.
         ! dbwl_nu(2, nbcalc_ks, nbsum, natom3)
         !dbwl_nu = one !; dbwl_nu(:,:,:,1:3) = zero
         !sqrt_phfrq0 = one
         do in=1,nbcalc_ks
           do im=1,nbsum
             cvec1 = dcmplx(sqrt_phfrq0(:) * dbwl_nu(1, in, im, :), &
                            sqrt_phfrq0(:) * dbwl_nu(2, in, im, :))
             !cvec1 = cone * tol6
             do ii=1,natom3
               !cvec2 = cone * tol6
               cvec2 = cmat(ii,:)
               !cvec2 = real(cmat(ii,:)
               hka_mn(ii, im, in) = dot_product(dconjg(cvec2), cvec1)
               !hka_mn(ii, im, in) = dconjg(hka_mn(ii, im, in))
               !hka_mn(ii, im, in) = cone
               !write(std_out, *)"hka_mn(ii, im, in)", hka_mn(ii, im, in)
             end do
           end do
         end do

         ! Sum over modes for this q-point.
         do nu=1,natom3
           ! Ignore acoustic or unstable modes.
           wqnu = phfrq(nu); if (wqnu < tol6) cycle

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
               !tpp(ip1,ip2) = dconjg(dka) * dkap + dconjg(dkpa) * dkpap
               !write(std_out,*)"tpp: ",tpp(ip1, ip2)

               dka   = dcmplx(displ_red(1, idir1, ipert1, nu), displ_red(2, idir1, ipert1, nu))
               dkap  = dcmplx(displ_red(1, idir2, ipert1, nu), displ_red(2, idir2, ipert1, nu))
               dkpa  = dcmplx(displ_red(1, idir1, ipert2, nu), displ_red(2, idir1, ipert2, nu))
               dkpap = dcmplx(displ_red(1, idir2, ipert2, nu), displ_red(2, idir2, ipert2, nu))
               tpp_red(ip1,ip2) = dka * dconjg(dkap) + dkpa * dconjg(dkpap)
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

           ! Get phonon occupation for all temperatures.
           nqnu_tlist = occ_be(wqnu, sigma%kTmesh(:), zero)

           ! Sum over bands and add (static) DW contribution for the different temperatures.
if (dtset%useria == 0) then
           do ibsum=1,nbsum
             do ib_k=1,nbcalc_ks
               ! Compute DW term following XG paper. Check prefactor.
               ! gkk0_atm(2, nbcalc_ks, nbsum, natom3)
               ! previous version
               gdw2_mn(ibsum, ib_k) = zero
               do ip2=1,natom3
                 do ip1=1,natom3
                   cfact = ( &
                     + gkk0_atm(1, ib_k, ibsum, ip1) * gkk0_atm(1, ib_k, ibsum, ip2) &
                     + gkk0_atm(2, ib_k, ibsum, ip1) * gkk0_atm(2, ib_k, ibsum, ip2) &
                     + gkk0_atm(1, ib_k, ibsum, ip2) * gkk0_atm(1, ib_k, ibsum, ip1) &
                     + gkk0_atm(2, ib_k, ibsum, ip2) * gkk0_atm(2, ib_k, ibsum, ip1) &
                   )

                   gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) + real(tpp_red(ip1,ip2) * cfact)
                 end do
               end do

               !gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) * two
               gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) /  (four * two * wqnu)
               !write(std_out,*)"gdw2_mn: ",gdw2_mn(ibsum, ib_k)
               ! dbwl_nu(2, nbcalc_ks, nbsum, natom3), gkk_nu(2, nbcalc_ks, natom3)
             end do ! ibsum
           end do ! ib_k
end if

           do ibsum=1,nbsum
             eig0mk = ebands%eig(ibsum, ik_ibz, spin)

             do ib_k=1,nbcalc_ks
               band = ib_k + bstart_ks - 1
               eig0nk = ebands%eig(band, ik_ibz, spin)
               ! Handle n == m and degenerate states (either ignore or add broadening)
               cplx_ediff = (eig0nk - eig0mk)
               if (abs(cplx_ediff) < tol6) cycle
               !if (abs(cplx_ediff) < tol6) cplx_ediff = cplx_ediff + sigma%ieta

               ! Optionally, accumulate contribution to Eliashberg functions (DW term)
               if (sigma%gfw_nomega > 0) then
                  sigma%gf_nnuq(ib_k, nu, iq_ibz, 3) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 3) &
                      - gdw2_mn(ibsum, ib_k) / real(cplx_ediff)
               end if

               ! Accumulate DW for each T, add it to Sigma(e0) and Sigma(w) as well
               ! - (2 n_{q\nu} + 1) * gdw2 / (e_{nk} - e_{mk})
               do it=1,sigma%ntemp
                 cfact = - weigth_q * gdw2_mn(ibsum, ib_k) * (two * nqnu_tlist(it) + one)  / cplx_ediff
                 rfact = real(cfact)
                 sigma%dw_vals(it, ib_k) = sigma%dw_vals(it, ib_k) + rfact
                 sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + rfact
                 if (sigma%nwr > 0) sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + rfact
               end do
             end do ! ib_k
           end do ! ib_sum

         end do ! nu
       end do ! iq_ibz

       ABI_FREE(gdw2_mn)
       ABI_FREE(tpp)
       ABI_FREE(tpp_red)
       ABI_FREE(hka_mn)
       ABI_FREE(dbwl_nu)
       ABI_FREE(gkk0_atm)
     end if ! not %imag_only

     if (sigma%gfw_nomega /= 0) then
       ! Accumulate DW Eliashberg function with gaussian.
       call xmpi_sum(sigma%gf_nnuq, comm, ierr)
       ABI_MALLOC(dargs, (sigma%gfw_nomega))
       ABI_MALLOC(dt_weights, (sigma%gfw_nomega, 2))
       sigma%gfw_vals = zero
       do iq_ibz=1,sigma%nqibz_k
         if (mod(iq_ibz, nprocs) /= my_rank) cycle  ! MPI parallelism
         call ifc_fourq(ifc, cryst, sigma%qibz_k(:,iq_ibz), phfrq, displ_cart) ! TODO: Use phfreqs in ephgw
         do nu=1,natom3
           dargs = sigma%gfw_mesh - phfrq(nu)
           dt_weights(:,1) = dirac_delta(dargs, dtset%ph_smear)
           do ib_k=1,nbcalc_ks
             do ii=1,3
               sigma%gfw_vals(:, ii, ib_k) = sigma%gfw_vals(:, ii, ib_k) +  &
                 sigma%wtq_k(iq_ibz) * sigma%gf_nnuq(ib_k, nu, iq_ibz, ii) * dt_weights(:, 1)
             end do
           end do
         end do
       end do
       ABI_FREE(dargs)
       ABI_FREE(dt_weights)
       call xmpi_sum(sigma%gfw_vals, comm, ierr)

       ! For tetrahedron method.
       !do nu=1,natom3
       !  do iq_ibz=1,sigma%nqibz_k
       !    do ib_k=1,nbcalc_ks
       !        vals_ibz_k = sigma%gf_nnuq(ib_k, nu, :, ii)
       !        call tetra_get_onewk(sigma%ephwg%tetra_k, iq_ibz, bcorr0, sigma%gfw_nomega, sigma%nqibz_k, vals_ibz_k, &
       !          sigma%gfw_mesh(1), sigma%gfw_mesh(sigma%gfw_nomega), max_occ1, dt_weights)

       !        sigma%gfw_vals(:, ii, ib_k) = sigma%gfw_vals(:, ii, ib_k) +  &
       !          sigma%gf_nnuq(ib_k, nu, iq_ibz, ii) * dt_weights(:, 1)  ! sigma%wtq_k(iq_ibz) *
       !    end do
       !  end do
       !end do
     end if

     ! Collect results inside comm and write results for this (k-point, spin) to NETCDF file.
     call sigmaph_gather_and_write(sigma, ebands, ikcalc, spin, comm)
   end do ! spin

   ABI_FREE(indq2dvdb)
   ABI_FREE(kg_k)
   ABI_FREE(kg_kq)
   ABI_FREE(ylm_k)
   ABI_FREE(ylm_kq)
   ABI_FREE(ylmgr_kq)
 end do !ikcalc

 call cwtime(cpu_all, wall_all, gflops_all, "stop")
 call wrtout(std_out, "Computation of Sigma_eph completed", do_flush=.True.)
 call wrtout(std_out, sjoin("Interpolate dVscf:", sec2str(cpu_dvscf), ", Total cpu time:", sec2str(wall_dvscf)))
 call wrtout(std_out, sjoin("Total wall-time:", sec2str(cpu_all), ", Total cpu time:", sec2str(wall_all), ch10, ch10))

 ! Free memory
 ABI_FREE(gvnl1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(nqnu_tlist)
 ABI_FREE(iqlk)
 if (sigma%nwr > 0) then
   ABI_FREE(cfact_wr)
 end if

 if (sigma%use_doublegrid) then
   ABI_FREE(eph_dg_mapping)
   ABI_FREE(indkk)
 endif

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

   ! Transform the gkk from (atom, reduced direction) basis to phonon mode representation
   do ipc=1,3*natom
     gkk_nu(1,:,:,:,nu) = gkk_nu(1,:,:,:,nu) &
       + gkk_atm(1,:,:,:,ipc) * displ_red(1,ipc,nu) &
       - gkk_atm(2,:,:,:,ipc) * displ_red(2,ipc,nu)
       !+ gkk_atm(2,:,:,:,ipc) * displ_red(2,ipc,nu)
     gkk_nu(2,:,:,:,nu) = gkk_nu(2,:,:,:,nu) &
       + gkk_atm(1,:,:,:,ipc) * displ_red(2,ipc,nu) &
       + gkk_atm(2,:,:,:,ipc) * displ_red(1,ipc,nu)
       !- gkk_atm(2,:,:,:,ipc) * displ_red(1,ipc,nu)
   end do

   gkk_nu(:,:,:,:,nu) = gkk_nu(:,:,:,:,nu) / sqrt(two * phfrq(nu))
 end do

end subroutine gkknu_from_atm
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
 integer,parameter :: master=0,occopt3=3,qptopt1=1,sppoldbl1=1
 integer :: my_rank,ik,my_nshiftq,my_mpw,cnt,nprocs,iq_ibz,ik_ibz,ndeg
 integer :: onpw,ii,ipw,ierr,it,spin,gap_err,ikcalc,gw_qprange,bstop,band
 integer :: nk_found,ifo,jj,bstart
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 logical :: downsample
 real(dp),parameter :: spinmagntarget=-99.99_dp,tol_enediff=0.001_dp*eV_Ha
 character(len=500) :: wfk_fname_dense
 real(dp) :: dksqmax,ph_wstep,elow,ehigh,wmax
 character(len=500) :: msg
 logical :: changed,found
 type(ebands_t) :: tmp_ebands, ebands_dense
 type(gaps_t) :: gaps
 type(hdr_type) :: hdr_wfk_dense
!arrays
 integer :: intp_kptrlatt(3,3)
 integer :: qptrlatt(3,3),indkk_k(1,6),my_gmax(3),kpos(6),band_block(2)
 integer :: val_indeces(ebands%nkpt, ebands%nsppol), intp_nshiftk
 real(dp):: params(3), nelect
 integer,allocatable :: gtmp(:,:),degblock(:,:)
 real(dp) :: my_shiftq(3,1),kk(3),kq(3),intp_shiftk(3)
 real(dp),pointer :: energies_dense(:,:,:)

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
 new%kTmesh = arth(dtset%tmesh(1), dtset%tmesh(2), new%ntemp) * kb_HaK

 ! Number of bands included in self-energy summation.
 ! This value depends on the kind of calculation as imag_only can use
 ! an energy window aroud the Fermi level.
 new%nbsum = dtset%mband

 gap_err = get_gaps(ebands, gaps)
 call gaps_print(gaps, unit=std_out)
 call ebands_report_gap(ebands, unit=std_out)
 val_indeces = get_valence_idx(ebands)

 ! Frequency mesh for sigma(w) and spectral function.
 ! TODO: Use GW variables but change default
 new%nwr = dtset%nfreqsp
 if (new%nwr > 0) then
   if (mod(new%nwr, 2) == 0) new%nwr = new%nwr + 1
 end if
 !dtset%freqspmin; dtset%freqspmax
 new%wr_step = 0.05 * eV_Ha

 ! Define q-mesh for integration of the self-energy.
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
 ! if symsigma == +1, we have to include all degenerate states in the set
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
     !new%bstart_ks = 1; new%nbcalc_ks = 6

     if (gw_qprange > 0) then
       ! All k-points: Add buffer of bands above and below the Fermi level.
       do spin=1,new%nsppol
         do ik=1,new%nkcalc
           new%bstart_ks(ik,spin) = max(val_indeces(ik,spin) - gw_qprange, 1)
           new%nbcalc_ks(ik,spin) = min(val_indeces(ik,spin) + gw_qprange, dtset%mband)
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
 if (abs(new%symsigma) == 1) then
   ABI_DT_MALLOC(new%degtab, (new%nkcalc, new%nsppol))
 end if

 do ikcalc=1,new%nkcalc
   kk = new%kcalc(:,ikcalc)
   call listkk(dksqmax,cryst%gmet,indkk_k,ebands%kptns,kk,ebands%nkpt,1,cryst%nsym,&
      sppoldbl1,cryst%symafm,cryst%symrel,new%timrev,use_symrec=.False.)

   new%kcalc2ibz(ikcalc, :) = indkk_k(1, :)

   ik_ibz = indkk_k(1,1) !; isym_k = indkk_k(1,2)
   !trev_k = indkk_k(1, 6); g0_k = indkk_k(1, 3:5)
   !isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !kk_ibz = ebands%kptns(:,ik_ibz)

   if (dksqmax > tol12) then
      write(msg, '(4a,es16.6,7a)' )&
       "The WFK file cannot be used to compute self-energy corrections at kpoint: ",ktoa(kk),ch10,&
       "the k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
       "Q-mesh: ",trim(ltoa(new%ngqpt)),", K-mesh (from kptrlatt) ",trim(ltoa(get_diag(dtset%kptrlatt))),ch10, &
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
       bstop = new%bstart_ks(ikcalc,spin) + new%nbcalc_ks(ikcalc,spin) - 1
       call enclose_degbands(ebands, ik_ibz, spin, new%bstart_ks(ikcalc,spin), bstop, &
         changed, tol_enediff, degblock=degblock)
       if (changed) then
         new%nbcalc_ks(ikcalc,spin) = bstop - new%bstart_ks(ikcalc,spin) + 1
         write(msg,'(2(a,i0),2a,2(1x,i0))')&
           "Not all the degenerate states at ikcalc= ",ikcalc,", spin= ",spin,ch10,&
           "were included in the bdgw set. bdgw has been changed to: ",new%bstart_ks(ikcalc,spin),bstop
         MSG_COMMENT(msg)
         write(msg,'(2(a,i0),2a)') "The number of included states ",bstop,&
                      " is larger than the number of bands in the input ",dtset%nband(new%nkcalc*(spin-1)+ikcalc),ch10,&
                      "Please increase nband."
         ABI_CHECK( bstop <= dtset%nband(new%nkcalc*(spin-1)+ikcalc), msg )
       end if
       ! Store band indices used for averaging (shifted by bstart_ks)
       ndeg = size(degblock, dim=2)
       ABI_DT_MALLOC(new%degtab(ikcalc, spin)%bids, (ndeg))
       do ii=1,ndeg
         cnt = degblock(2, ii) - degblock(1, ii) + 1
         ABI_DT_MALLOC(new%degtab(ikcalc, spin)%bids(ii)%vals, (cnt))
         new%degtab(ikcalc, spin)%bids(ii)%vals = [(jj, jj=degblock(1, ii) - new%bstart_ks(ikcalc, spin) + 1, &
                                                           degblock(2, ii) - new%bstart_ks(ikcalc, spin) + 1)]
       end do
       ABI_FREE(degblock)
     end do
   end if ! symsigma

 end do

 ! Now we can finally compute max_nbcalc
 new%max_nbcalc = maxval(new%nbcalc_ks)

 ! Now compute the min/max KS energy to be included in the imaginary part.
 ! ifc%omega_minmax(2) comes froms the coarse Q-mesh of the DDB so increase it by 10%.
 ! Also take into account Lorentzian shape if zcut is used.
 elow = huge(one); ehigh = - huge(one)
 wmax = 1.1_dp * ifc%omega_minmax(2) + five * dtset%zcut
 do ikcalc=1,new%nkcalc
   ik_ibz = new%kcalc2ibz(ikcalc, 1)
   do spin=1,new%nsppol
     bstart = new%bstart_ks(ikcalc,spin)
     bstop = new%bstart_ks(ikcalc,spin) + new%nbcalc_ks(ikcalc,spin) - 1
     ehigh = max(ehigh, maxval(ebands%eig(bstart:bstop, ik_ibz, spin)) + wmax)
     elow = min(elow, minval(ebands%eig(bstart:bstop, ik_ibz, spin)) - wmax)
   end do
 end do

 !call ebands_get_binds_from_range(ebands, elow, ehigh, bstart, bstop)
 bstart = huge(1); bstop = -huge(1)
 do spin=1,ebands%nsppol
   do ik=1,ebands%nkpt
     do band=1,ebands%nband(ik+(spin-1)*ebands%nkpt)
       if (ebands%eig(band, ik , spin) >= elow .and. ebands%eig(band, ik , spin) <= ehigh) then
          bstart = min(bstart, band)
          bstop = max(bstop, band)
       end if
     end do
   end do
 end do

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! used to symmetrize the wavefunctions in G-space.
 ! TODO: Should loop over IBZ(k)
 new%mpw = 0; new%gmax = 0; cnt = 0
 do ik=1,new%nkcalc
   kk = new%kcalc(:, ik)
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
 call wrtout(std_out, sjoin('Optimal value of mpw=', itoa(new%mpw)))

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
 ! TODO: Input variable or activate it by default?
 new%gfw_nomega = 0
 new%gfw_nomega = 200; ph_wstep = dtset%ph_wstep
 if (new%gfw_nomega /= 0) then
   new%gfw_nomega = nint((ifc%omega_minmax(2) - ifc%omega_minmax(1) ) / ph_wstep) + 1
   ABI_MALLOC(new%gfw_mesh, (new%gfw_nomega))
   new%gfw_mesh = arth(ifc%omega_minmax(1), ph_wstep, new%gfw_nomega)
   ABI_MALLOC(new%gfw_vals, (new%gfw_nomega, 3, new%max_nbcalc))
 end if

 new%imag_only = .False.; if (dtset%eph_task == -4) new%imag_only = .True.
 ! TODO: Remove qint_method, use eph_intmeth or perhaps dtset%qint_method dtset%kint_method
 ! DEcide default behaviour for Re-Im/Im
 new%qint_method = dtset%eph_intmeth - 1
 !new%qint_method = 0; if (dtset%userib == 23) new%qint_method = 1
 !write(std_out, *)"imag_only:", new%imag_only, ", ;qint_method:", new%qint_method

 ! Initialize object for the computation of integration weights (integration in q-space).
 ! Weights can be obtained in different ways:
 !
 !  1. Computed from eigens on the same coarse q-mesh as the one used for the self-energy.
 !  2. Obtained from eigens on a denser q-mesh and then transfered to the coarse q-mesh.
 !     In this case the eigens on the dense mesh are either read from an external file (ab-initio)
 !     or interpolated on the fly with star-functions.
 !
 !  NB: The routines assume that the k-mesh for electrons and the q-mesh for phonons are the same.
 !  Thus we need to downsample the k-mesh if it's denser that the q-mesh.


 new%use_doublegrid = .False.
 if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
     (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then

   wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)//'FINE'
   call wrtout(std_out,"EPH Interpolation: will read energies from: "//trim(wfk_fname_dense),"COLL")

   if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
     MSG_ERROR(msg)
   end if

   call wfk_read_eigenvalues(wfk_fname_dense,energies_dense,hdr_wfk_dense,comm)
   ebands_dense = ebands_from_hdr(hdr_wfk_dense,maxval(hdr_wfk_dense%nband),energies_dense)
   ABI_FREE(energies_dense)

   !TODO add a check for consistency
   ! number of bands and kpoints (comensurability)
   ABI_CHECK(hdr_wfk_dense%mband == ebands%mband, 'Inconsistent number of bands for the fine and dense grid')
   call hdr_free(hdr_wfk_dense)

   new%use_doublegrid = .True.

 !read bs_interpmult
 else if (dtset%bs_interp_kmult(1) /= 0 .or.&
          dtset%bs_interp_kmult(2) /= 0 .or.&
          dtset%bs_interp_kmult(3) /= 0 ) then

   call wrtout(std_out,"EPH Interpolation: will use star functions interpolation","COLL")

   ! Interpolate band energies with star-functions
   params = 0; params(1) = 1; params(2) = 5
   !TODO: mband should be min of nband
   band_block = [1, ebands%mband]
   intp_kptrlatt(:,1) = [ebands%kptrlatt(1,1)*dtset%bs_interp_kmult(1), 0, 0]
   intp_kptrlatt(:,2) = [0, ebands%kptrlatt(2,2)*dtset%bs_interp_kmult(2), 0]
   intp_kptrlatt(:,3) = [0, 0, ebands%kptrlatt(3,3)*dtset%bs_interp_kmult(3)]

   intp_shiftk = zero
   intp_nshiftk = 1
   ebands_dense = ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt,&
                                      intp_nshiftk, intp_shiftk, band_block, comm)
   new%use_doublegrid = .True.
 end if

 if (new%qint_method > 0) then
   ! bstart and new%bsum select the band range.
   bstart = 1
   if (new%use_doublegrid) then
     ! Double-grid technique from ab-initio energies or star-function interpolation.
     !if (dtset% dtftil% ....) then
     !  tmp_ebands = wfk_read_ebands(filepath, comm)
     !else
     !  do ii = 1, 3
     !    intp_kptrlatt(:, ii) = dtset%bs_interp_kmult(ii) * ebands%kptrlatt(:, ii)
     !  end do
     !  tmp_ebands = ebands_interp_kmesh(ebands, cryst, dtset%einterp, &
     !    intp_kptrlatt, my_nshiftq, my_shiftq, [bstart, bstart + new%nbsum - 1], comm)
     !  call ebands_interpolate_kpath(ebands, dtset, cryst, band_block, prefix, comm)
     !end if
     new%ephwg          = ephwg_from_ebands(cryst, ifc, ebands_dense, bstart, new%nbsum, comm)
     new%eph_doublegrid = eph_double_grid_new(cryst, ebands_dense, ebands%kptrlatt, ebands_dense%kptrlatt)
     !call ebands_free(tmp_ebands)
   else
     downsample = any(ebands%kptrlatt /= qptrlatt) .or. ebands%nshiftk /= my_nshiftq
     if (ebands%nshiftk == my_nshiftq) downsample = downsample .or. any(ebands%shiftk /= my_shiftq)
     if (downsample) then
       MSG_COMMENT("K-mesh != Q-mesh for self-energy. Will downsample electron energies.")
       tmp_ebands = ebands_downsample(ebands, cryst, qptrlatt, my_nshiftq, my_shiftq)
       new%ephwg = ephwg_from_ebands(cryst, ifc, tmp_ebands, bstart, new%nbsum, comm)
       call ebands_free(tmp_ebands)
     else
       new%ephwg = ephwg_from_ebands(cryst, ifc, ebands, bstart, new%nbsum, comm)
     end if
   end if
 else
   if (new%use_doublegrid) then
     new%eph_doublegrid = eph_double_grid_new(cryst, ebands_dense, ebands%kptrlatt, ebands_dense%kptrlatt)
     new%ephwg          = ephwg_from_ebands(cryst, ifc, ebands_dense, bstart, new%nbsum, comm)
   endif
 end if

 ! Compute the chemical potential at the different physical temperatures with Fermi-Dirac.
 ABI_MALLOC(new%mu_e, (new%ntemp))
 new%mu_e = ebands%fermie

 if (dtset%eph_fermie == zero) then
   if (new%use_doublegrid) then
     call ebands_copy(ebands_dense, tmp_ebands)
   else
     call ebands_copy(ebands, tmp_ebands)
   endif
   !
   do it=1,new%ntemp
     call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(it), spinmagntarget, dtset%prtvol)
     call ebands_set_nelect(tmp_ebands, ebands%nelect, dtset%spinmagntarget, msg)
     new%mu_e(it) = tmp_ebands%fermie
   end do
   !
   ! Check that the total number of electrons is correct
   ! This is to trigger problems as the routines that calculate the occupations in ebands_set_nelect
   ! are different from the occ_fd that will be used in the rest of the subroutine
   !
   do it=1,new%ntemp
     nelect = 0
     !loop over spin
     do spin=1,tmp_ebands%nsppol
       !loop over kpoints
       do ik=1,tmp_ebands%nkpt
         !loop over bands
         do ii=1,tmp_ebands%nband(ik)
           nelect = nelect + 2*tmp_ebands%wtk(ik)*occ_fd(tmp_ebands%eig(ii,ik,spin),new%kTmesh(it),new%mu_e(it))
         end do
       end do
     end do
     !
     write(msg,'(3(a,f10.6))')&
       'Calculated number of electrons nelect = ',nelect,&
       ' does not correspond with ebands%nelect = ',tmp_ebands%nelect,&
       ' for T = ',new%kTmesh(it)
     ! For T = 0 the number of occupied states goes in discrete steps (according to the k-point sampling)
     ! for finite doping its hard to find nelect that exactly matches ebands%nelect.
     ! in this case we print a warning
     if (new%kTmesh(it) == 0) then
       MSG_WARNING(msg)
     else
       ABI_CHECK(abs(nelect-ebands%nelect) < tol6,msg)
     end if
     !
   end do
   call ebands_free(tmp_ebands)
 else
   new%mu_e(:) = ebands%fermie
 endif

 call ebands_free(ebands_dense)

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

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
       "symsigma", "nbsum", "symdynmat", "ph_intmeth", "eph_intmeth", "qint_method", "eph_transport", "imag_only"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
       "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"])
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
       nctkarr_t("gfw_vals", "dp", "gfw_nomega, three, max_nbcalc, nkcalc, nsppol") &
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
   ii = 0; if (new%imag_only) ii = 1
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "symsigma", "nbsum", "symdynmat", "ph_intmeth", "eph_intmeth", "qint_method", "eph_transport", "imag_only"], &
     [new%symsigma, new%nbsum, dtset%symdynmat, dtset%ph_intmeth, dtset%eph_intmeth, new%qint_method, &
     dtset%eph_transport, ii])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"], &
     [aimag(new%ieta), new%wr_step, dtset%eph_fsewin, dtset%eph_fsmear, dtset%eph_extrael, dtset%eph_fermie, &
     dtset%ph_wstep, dtset%ph_smear])
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
 if (allocated(self%e0vals)) then
   ABI_FREE(self%e0vals)
 end if
 if (allocated(self%cweights)) then
   ABI_FREE(self%cweights)
 end if
 if (allocated(self%deltaw_pm)) then
   ABI_FREE(self%deltaw_pm)
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
 if (allocated(self%indkk)) then
   ABI_FREE(self%indkk)
 end if
 if (allocated(self%wtq_k)) then
   ABI_FREE(self%wtq_k)
 end if
 if (allocated(self%gfw_mesh)) then
   ABI_FREE(self%gfw_mesh)
 end if
 if (allocated(self%gf_nnuq)) then
   ABI_FREE(self%gf_nnuq)
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
 call ephwg_free(self%ephwg)

 if (self%use_doublegrid) then
   call eph_double_grid_free(self%eph_doublegrid)
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
!!  prtbol= Verbosity level
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_setup_kcalc(self, cryst, ikcalc, prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_setup_kcalc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, prtvol
 type(crystal_t),intent(in) :: cryst
 type(sigmaph_t),intent(inout) :: self

!Local variables-------------------------------
 integer,parameter :: sppoldbl1 = 1
 real(dp) :: dksqmax
 type(lgroup_t) :: lgk

! *************************************************************************

 if (allocated(self%qibz_k)) then
   ABI_FREE(self%qibz_k)
 end if
 if (allocated(self%wtq_k)) then
   ABI_FREE(self%wtq_k)
 end if

 ! Prepare weights for BZ(k) integration
 if (self%qint_method > 0) call ephwg_setup_kpoint(self%ephwg, self%kcalc(:, ikcalc), prtvol)

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

 ! DEBUGGING
 ! Get mapping IBZ_k --> initial IBZ (self%lgrp%ibz --> self%ibz)
 if (allocated(self%indkk)) then
   ABI_FREE(self%indkk)
 end if
 ABI_MALLOC(self%indkk, (self%nqibz_k, 6))
 call listkk(dksqmax, cryst%gmet, self%indkk, self%qibz, self%qibz_k, self%nqibz, self%nqibz_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)
 if (dksqmax > tol12) MSG_ERROR("Wrong mapping")

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
 real(dp) :: smrt,e0pde(9)
 complex(dpc) :: sig0c,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2
 character(len=500) :: msg
!arrays
 integer :: shape3(3),shape4(4),shape5(5),shape6(6)
 integer, ABI_CONTIGUOUS pointer :: bids(:)
 real(dp), ABI_CONTIGUOUS pointer :: rdata3(:,:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:),rdata6(:,:,:,:,:,:)
 real(dp) :: qp_gaps(self%ntemp),qpadb_gaps(self%ntemp)
 real(dp),allocatable :: aw(:,:,:)
 real(dp) :: ks_enes(self%max_nbcalc),ze0_vals(self%ntemp, self%max_nbcalc)
 real(dp) :: gfw_avg(self%gfw_nomega, 3)
 complex(dpc),target :: qpadb_enes(self%ntemp, self%max_nbcalc),qp_enes(self%ntemp, self%max_nbcalc)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call xmpi_sum_master(self%vals_e0ks, master, comm, ierr)
 call xmpi_sum_master(self%dvals_de0ks, master, comm, ierr)
 call xmpi_sum_master(self%dw_vals, master, comm, ierr)
 if (self%nwr > 0) call xmpi_sum_master(self%vals_wr, master, comm, ierr)

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

     ! Symmetrize Eliashberg function
     if (self%gfw_nomega > 0) then
       gfw_avg = sum(self%gfw_vals(:, :, bids(:)), dim=3) / nstates
       do ii=1,nstates
         self%gfw_vals(:, :, bids(ii)) = gfw_avg
       end do
     end if

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
     ! FIXME
     if (it == 1) then
       !   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP
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

#ifdef HAVE_NETCDF
 ! **Only master writes**
 ! Write self-energy matrix elements for this (kpt, spin)
 ! (use iso_c_binding to associate a real pointer to complex data because netcdf does not support complex types).
 ! Well, cannot use c_loc with gcc <= 4.8 due to internal compiler error so use c2r and stack memory.
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
   ! A = -1/pi [Im Sigma(ww)] / ([ww - ee - Re Sigma(ww)] ** 2 + Im Sigma(ww) ** 2])
   ABI_MALLOC(aw, (self%nwr, self%ntemp, self%max_nbcalc))
   do ib=1,self%nbcalc_ks(ikcalc,spin)
     band = self%bstart_ks(ikcalc, spin) + ib - 1
     kse = ebands%eig(band, ik_ibz, spin)
     do it=1,self%ntemp
       aw(:, it, ib) = piinv * abs(aimag(self%vals_wr(:, it, ib))) / &
         ((self%wrmesh_b(:, ib) - kse - real(self%vals_wr(:, it, ib))) ** 2 + aimag(self%vals_wr(:, it, ib)) ** 2)
     end do
   end do
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "spfunc_wr"), aw, start=[1, 1, 1, ikcalc, spin]))
   ABI_FREE(aw)
 end if

 ! Write Eliashberg functions
 if (self%gfw_nomega > 0) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "gfw_vals"), self%gfw_vals, start=[1, 1, 1, ikcalc, spin]))
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
 integer :: ikc,is
 character(len=500) :: msg

! *************************************************************************

 if (unt == dev_null) return

 ! Write dimensions
 write(unt,"(a)")sjoin("Number of bands in e-ph self-energy:", itoa(self%nbsum))
 write(unt,"(a)")sjoin("Symsigma: ",itoa(self%symsigma), "Timrev:", itoa(self%timrev))
 write(unt,"(a)")sjoin("Imaginary shift in the denominator (zcut): ", ftoa(aimag(self%ieta) * Ha_eV, fmt="f5.3"), "[eV]")
 msg = "Standard quadrature"; if (self%qint_method == 1) msg = "tetrahedron method"
 write(unt, "(2a)")sjoin("Method for q-space integration:", msg)
 if (self%imag_only) write(unt, "(a)")"Only the Imaginary part of Sigma will be computed."
 if (.not. self%imag_only) write(unt, "(a)")"Both Real and Imaginary part of Sigma will be computed."
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




!!****f* m_sigmaph/eph_double_grid_new
!! NAME
!!
!! FUNCTION
!!   Prepare Double grid integration
!!
!!   double grid:
!!   ----------------- interp_kmult 2
!!   |. .|.|.|.|. . .| side 1
!!   |. x|.|x|.|. x .| size 3 (2*side+1)
!!   |. .|.|.|.|. . .|
!!   -----------------
!!
!!   triple grid:
!!   ------------------- interp_kmult 3
!!   |. . .|. . .|. . .| side 1
!!   |. x .|. x .|. x .| size 3 (2*side+1)
!!   |. . .|. . .|. . .|
!!   -------------------
!!
!!   quadruple grid:
!!   --------------------------- interp_kmult 4
!!   |. . . .|.|. . .|.|. . . .| side 2
!!   |. . . .|.|. . .|.|. . . .| size 5  (2*side+1)
!!   |. . x .|.|. x .|.|. x . .|
!!   |. . . .|.|. . .|.|. . . .|
!!   |. . . .|.|. . .|.|. . . .|
!!   ---------------------------
!!
!!   . = double grid
!!   x = coarse grid
!!
!!   The fine grid is used to evaluate the weights on the coarse grid
!!   The points of the fine grid are associated to the points of the
!!   coarse grid according to proximity
!!   The integration weights are returned on the coarse grid.
!!
!!   Steps of the implementation
!!   1. Get the fine k grid from file or interpolation (ebands_dense)
!!   2. Find the matching between the k_coarse and k_dense using the double_grid object
!!   3. Calculate the phonon frequencies on the dense mesh and store them on a array
!!   4. Create an array to bring the points in the full brillouin zone to the irreducible brillouin zone
!!   5. Create a scatter array between the points in the fine grid
!!
!! INPUTS
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

type (eph_double_grid_t) function eph_double_grid_new(cryst, ebands_dense, kptrlatt_coarse, kptrlatt_dense) result(eph_dg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph_double_grid_new'
!End of the abilint section

 type(crystal_t), intent(in) :: cryst
 type(ebands_t), intent(in) :: ebands_dense
 integer, intent(in) :: kptrlatt_coarse(3,3), kptrlatt_dense(3,3)

 integer,parameter :: sppoldbl1=1,timrev1=1
 integer :: i_dense,i_coarse,this_dense,i_subdense,i1,i2,i3,ii,jj,kk
 integer :: nkpt_coarse(3), nkpt_dense(3), interp_kmult(3), interp_side(3)
 integer,allocatable :: indqq(:,:)
 real(dp) :: dksqmax

 nkpt_coarse(1) = kptrlatt_coarse(1,1)
 nkpt_coarse(2) = kptrlatt_coarse(2,2)
 nkpt_coarse(3) = kptrlatt_coarse(3,3)
 nkpt_dense(1)  = kptrlatt_dense(1,1)
 nkpt_dense(2)  = kptrlatt_dense(2,2)
 nkpt_dense(3)  = kptrlatt_dense(3,3)

 eph_dg%dense_nbz = nkpt_dense(1)*nkpt_dense(2)*nkpt_dense(3)
 eph_dg%coarse_nbz = nkpt_coarse(1)*nkpt_coarse(2)*nkpt_coarse(3)
 interp_kmult = nkpt_dense/nkpt_coarse
 eph_dg%interp_kmult = interp_kmult
 eph_dg%nkpt_coarse = nkpt_coarse
 eph_dg%nkpt_dense = nkpt_dense
 eph_dg%ebands_dense = ebands_dense

 ! A microzone is the set of points in the fine grid belonging to a certain coarse point
 ! we have to consider a side of a certain size around the coarse point
 ! to make sure the microzone is centered around its point.
 ! The fine points shared by multiple microzones should have weights
 ! according to in how many microzones they appear
 !
 ! this is integer division
 interp_side = interp_kmult/2

 eph_dg%ndiv = (2*interp_side(1)+1)*&
               (2*interp_side(2)+1)*&
               (2*interp_side(3)+1)


 write(std_out,*) 'coarse:      ', nkpt_coarse
 write(std_out,*) 'dense:       ', nkpt_dense
 write(std_out,*) 'interp_kmult:', interp_kmult
 write(std_out,*) 'ndiv:        ', eph_dg%ndiv
 ABI_CHECK(all(nkpt_dense(:) >= nkpt_coarse(:)), 'dense mesh is smaller than coarse mesh.')

 ABI_MALLOC(eph_dg%kpts_coarse,(3,eph_dg%coarse_nbz))
 ABI_MALLOC(eph_dg%kpts_dense,(3,eph_dg%dense_nbz))
 ABI_MALLOC(eph_dg%coarse_to_dense,(eph_dg%coarse_nbz,eph_dg%ndiv))

 ABI_MALLOC(eph_dg%dense_to_indexes,(3,eph_dg%dense_nbz))
 ABI_MALLOC(eph_dg%indexes_to_dense,(nkpt_dense(1),nkpt_dense(2),nkpt_dense(3)))

 ABI_MALLOC(eph_dg%coarse_to_indexes,(3,eph_dg%dense_nbz))
 ABI_MALLOC(eph_dg%indexes_to_coarse,(nkpt_coarse(1),nkpt_coarse(2),nkpt_coarse(3)))

 ABI_MALLOC(eph_dg%weights_dense,(eph_dg%dense_nbz))

 write(std_out,*) 'create dense to coarse mapping'
 ! generate mapping of points in dense bz to the dense bz
 ! coarse loop
 i_dense = 0
 i_coarse = 0
 do kk=1,nkpt_coarse(3)
   do jj=1,nkpt_coarse(2)
     do ii=1,nkpt_coarse(1)
       i_coarse = i_coarse + 1
       !calculate reduced coordinates of point in coarse mesh
       eph_dg%kpts_coarse(:,i_coarse) = [dble(ii-1)/nkpt_coarse(1),&
                                         dble(jj-1)/nkpt_coarse(2),&
                                         dble(kk-1)/nkpt_coarse(3)]
       !create the fine mesh
       do i3=1,interp_kmult(3)
         do i2=1,interp_kmult(2)
           do i1=1,interp_kmult(1)
             i_dense = i_dense + 1
             !calculate reduced coordinates of point in dense mesh
             eph_dg%kpts_dense(:,i_dense) =  &
                  [dble((ii-1)*interp_kmult(1)+i1-1)/(nkpt_coarse(1)*interp_kmult(1)),&
                   dble((jj-1)*interp_kmult(2)+i2-1)/(nkpt_coarse(2)*interp_kmult(2)),&
                   dble((kk-1)*interp_kmult(3)+i3-1)/(nkpt_coarse(3)*interp_kmult(3))]
             !integer indexes mapping
             eph_dg%indexes_to_dense((ii-1)*interp_kmult(1)+i1,&
                                     (jj-1)*interp_kmult(2)+i2,&
                                     (kk-1)*interp_kmult(3)+i3) = i_dense
             eph_dg%dense_to_indexes(:,i_dense) = [(ii-1)*interp_kmult(1)+i1,&
                                                   (jj-1)*interp_kmult(2)+i2,&
                                                   (kk-1)*interp_kmult(3)+i3]
           enddo
         enddo
       enddo
       eph_dg%indexes_to_coarse(ii,jj,kk) = i_coarse
       eph_dg%coarse_to_indexes(:,i_coarse) = [ii,jj,kk]
     enddo
   enddo
 enddo

 ! here we need to iterate again because we can have points of the dense grid
 ! belonging to multiple coarse points
 i_coarse = 0
 do kk=1,nkpt_coarse(3)
   do jj=1,nkpt_coarse(2)
     do ii=1,nkpt_coarse(1)
       i_coarse = i_coarse + 1

       !create a mapping from coarse to dense
       i_subdense = 0
       do i3=-interp_side(3),interp_side(3)
         do i2=-interp_side(2),interp_side(2)
           do i1=-interp_side(1),interp_side(1)
             i_subdense = i_subdense + 1
             !integer indexes mapping
             this_dense = eph_dg%indexes_to_dense(&
                    mod((ii-1)*interp_kmult(1)+i1+nkpt_dense(1),nkpt_dense(1))+1,&
                    mod((jj-1)*interp_kmult(2)+i2+nkpt_dense(2),nkpt_dense(2))+1,&
                    mod((kk-1)*interp_kmult(3)+i3+nkpt_dense(3),nkpt_dense(3))+1)

             !array indexes mapping
             eph_dg%coarse_to_dense(i_coarse,i_subdense) = this_dense
           enddo
         enddo
       enddo
     enddo
   enddo
 enddo

 ABI_CHECK(i_dense == eph_dg%dense_nbz, 'dense mesh mapping is incomplete')

 !calculate the weights of each fine point
 !different methods to distribute the weights might lead to better convergence
 !loop over coarse points
 eph_dg%weights_dense = 0
 do ii=1,eph_dg%coarse_nbz
   !loop over points in the microzone
   do jj=1,eph_dg%ndiv
     i_dense = eph_dg%coarse_to_dense(ii,jj)
     eph_dg%weights_dense(i_dense) = eph_dg%weights_dense(i_dense) + 1
   end do
 end do
 !weights_dense is array, ndiv is scalar
 eph_dg%weights_dense = 1/eph_dg%weights_dense/(interp_kmult(1)*interp_kmult(2)*interp_kmult(3))

 !3.
 ABI_MALLOC(eph_dg%kpts_dense_ibz,(3,ebands_dense%nkpt))
 eph_dg%kpts_dense_ibz = ebands_dense%kptns
 eph_dg%dense_nibz = ebands_dense%nkpt

 !4.
 write(std_out,*) 'map bz -> ibz'
 ! TODO: this is slow for large grids, improve!
 ABI_MALLOC(eph_dg%bz2ibz_dense,(eph_dg%dense_nbz))
 ABI_MALLOC(indqq,(eph_dg%dense_nbz,6))
 call listkk(dksqmax,cryst%gmet,    indqq,&
             eph_dg%kpts_dense_ibz, eph_dg%kpts_dense,&
             eph_dg%dense_nibz,     eph_dg%dense_nbz, &
             cryst%nsym,sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)
 ABI_CHECK(dksqmax < tol6, 'Problem creating a bz to ikbz kpoint mapping')
 eph_dg%bz2ibz_dense(:) = indqq(:,1)
 ABI_FREE(indqq)


end function eph_double_grid_new
!!***

!!****f* m_sigmaph/eph_double_grid_free
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine eph_double_grid_free(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph_double_grid_free'
!End of the abilint section

 type(eph_double_grid_t) :: self

 ABI_FREE(self%weights_dense)
 ABI_FREE(self%bz2ibz_dense)
 ABI_FREE(self%kpts_coarse)
 ABI_FREE(self%kpts_dense)
 ABI_FREE(self%kpts_dense_ibz)
 ABI_FREE(self%coarse_to_dense)
 ABI_FREE(self%dense_to_indexes)
 ABI_FREE(self%indexes_to_dense)
 ABI_FREE(self%coarse_to_indexes)
 ABI_FREE(self%indexes_to_coarse)

 !ABI_FREE(self%bz2ibz_coarse)

 !call ebands_free(ebands_dense)

end subroutine eph_double_grid_free

!------------------------------------------------------------------------

!!****f* m_sigmaph/eph_double_grid_get_index
!! NAME
!!
!! FUNCTION
!!   Get the indeex of a certain k-point in the double grid
!!
!! INPUTS
!!   kpt=kpoint to be mapped (reduced coordinates)
!!   opt=Map to the coarse (1) or dense grid (2)
!! 
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

type (integer) function eph_double_grid_get_index(self,kpt,opt) result(ikpt)

 type(eph_double_grid_t) :: self
 integer :: opt
 real(dp) :: kpt(3)

 if (opt==1) then
   ikpt = self%indexes_to_coarse(&
             mod(nint((kpt(1)+1)*self%nkpt_coarse(1)),self%nkpt_coarse(1))+1,&
             mod(nint((kpt(2)+1)*self%nkpt_coarse(2)),self%nkpt_coarse(2))+1,&
             mod(nint((kpt(3)+1)*self%nkpt_coarse(3)),self%nkpt_coarse(3))+1)
 else if (opt==2) then
   ikpt = self%indexes_to_dense(&
             mod(nint((kpt(1)+1)*self%nkpt_dense(1)),self%nkpt_dense(1))+1,&
             mod(nint((kpt(2)+1)*self%nkpt_dense(2)),self%nkpt_dense(2))+1,&
             mod(nint((kpt(3)+1)*self%nkpt_dense(3)),self%nkpt_dense(3))+1)
 else
   MSG_ERROR(sjoin("Error in eph_double_grid_get_index opt. Possible values are 1 or 2. Got", itoa(opt)))
 endif    

end function eph_double_grid_get_index

!!****f* m_sigmaph/eph_double_grid_debug
!! NAME
!!
!! FUNCTION
!!   Write different files to disk for debugging of the double grid calculation
!!
!! INPUTS
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine eph_double_grid_debug(self,cryst,ifc)

 type(eph_double_grid_t) :: self
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: cryst
 integer  :: ii, jj, i_dense
 real(dp) :: maxfreq, error, dksqmin, dksqmean, dksqmax
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom)
 real(dp),allocatable :: phfreq_bz(:), phfreq_ibz(:), phfrq_dense(:,:)
 real(dp) :: qpt(3)

 write(std_out,*) 'calculate phonon frequencies'
 !calculate the phonon frequencies at the q-points on the ibz of the dense q-grid
 !ibz version
 ! HM: I noticed that the fourier interpolation sometimes breaks the symmetries
 !     for low q-point sampling
 ABI_MALLOC(phfrq_dense,(3*cryst%natom,self%dense_nibz))
 do ii=1,self%dense_nibz
   qpt = self%kpts_dense_ibz(:,ii)
   ! Get phonon frequencies and displacements in reduced coordinates for this q-point
   call ifc_fourq(ifc, cryst, qpt, phfrq_dense(:,ii), displ_cart )
 enddo
 !bz version
 ABI_MALLOC(phfrq_dense,(3*cryst%natom,self%dense_nbz))
 do ii=1,self%dense_nbz
   qpt = self%kpts_dense(:,ii)
   ! Get phonon frequencies and displacements in reduced coordinates for this q-point
   call ifc_fourq(ifc, cryst, qpt, phfrq_dense(:,ii), displ_cart )
 enddo

 open (unit = 2, file = "ibz.dat")
 do ii=1,self%dense_nibz
   write(2,*) self%kpts_dense_ibz(:,ii)
 end do

 open (unit = 2, file = "bz.dat")
 do ii=1,self%dense_nbz
   write(2,*) self%kpts_dense(:,ii), self%bz2ibz_dense(ii)
 end do

 ABI_MALLOC(phfreq_bz,(cryst%natom*3))
 ABI_MALLOC(phfreq_ibz,(cryst%natom*3))
 open (unit = 2, file = "phbz.dat")
 open (unit = 3, file = "phibz.dat")
 dksqmax = 0
 dksqmin = 1e8
 maxfreq = 0
 do ii=1,self%dense_nbz
   call ifc_fourq(ifc, cryst, self%kpts_dense(:,ii), phfreq_bz, displ_cart )
   call ifc_fourq(ifc, cryst, self%kpts_dense_ibz(:,self%bz2ibz_dense(ii)), phfreq_ibz, displ_cart )
   write(2,*) phfreq_bz
   write(3,*) phfreq_ibz
   do jj=1,cryst%natom*3
     error = abs(phfreq_bz(jj)-phfreq_ibz(jj))
     if (dksqmax < error) dksqmax = error
     if (dksqmin > error) dksqmin = error
     if (maxfreq < phfreq_bz(jj)) maxfreq = phfreq_bz(jj)
     dksqmean = dksqmean + error
   enddo
 end do
 write(std_out,*) 'bz2ibz phonon error min: ', dksqmin
 write(std_out,*) 'bz2ibz phonon error max: ', dksqmax, dksqmax/maxfreq
 write(std_out,*) 'bz2ibz phonon error mean:', dksqmean/self%dense_nbz, &
                                         dksqmean/self%dense_nbz/maxfreq

 open (unit = 2, file = "coarse2dense.dat")
 do ii=1,self%coarse_nbz
   write(2,*)
   write(2,*)
   write(2,*) self%kpts_coarse(:,ii)
   do jj=1,self%ndiv
     i_dense = self%coarse_to_dense(ii,jj)
     write(2,*) self%kpts_dense(:,i_dense), self%weights_dense(i_dense)
   end do
 end do

 open (unit = 2, file = "coarse.dat")
 do ii=1,self%coarse_nbz
   write(2,*) self%kpts_coarse(:,ii)
 end do
 close(2)

 open (unit = 2, file = "dense.dat")
 do ii=1,self%dense_nbz
   write(2,*) self%kpts_dense(:,ii)
 end do
 close(2)

end subroutine eph_double_grid_debug


end module m_sigmaph
!!***

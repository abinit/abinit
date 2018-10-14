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
 use m_skw
 use m_lgroup
 use m_ephwg
 use m_nctk
 use m_hdr
 use m_eph_double_grid
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, simpson
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
 use m_pawang,         only : pawang_type, gauleg
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

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

  integer :: nwr
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Odd number so that the mesh is centered on the KS energy.
   ! The spectral function is computed only if nwr > 0 (taken from dtset%nfreqsp)

  integer :: ntemp
   ! Number of temperatures.

  integer :: symsigma
   ! 1 if matrix elements should be symmetrized.
   ! Required when the sum over q in the BZ is replaced by IBZ(k).

  integer :: timrev
   !  timrev=1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nbsum
   ! Total number of bands used in sum over states without taking into account MPI distribution.

  !integer :: my_nbsum
   ! Number of bands sum over states treated by this MPI rank

  !integer :: my_bstart, my_bstop
   ! Initial KS band index included in self-energy sum
   ! 1 if Re-Im
   ! computed at runtime on the basis of the nk states in Sigma_{nk} if imag_only

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

  integer :: ntheta, nphi
  !integer :: ndvis(3)=3

  integer :: nqr = 0
   ! Number of points on the radial mesh for spherical integration of the Frohlich self-energy

  integer :: angl_size = 0
   ! Dimension of angular mesh for spherical integration of the Frohlich self-energy
   ! angl_size = ntheta * nphi

  complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  real(dp) :: qrad
   ! Radius of the sphere for the numerical integration of the Frohlich self-energy

  !real(dp) :: qdamp
   ! Exponential damping added to Frohlich model.

  integer :: qint_method
   ! Defines the method used to integrate in q-space
   ! 0 --> Standard quadrature (one point per small box)
   ! 1 --> Use tetrahedron method

  integer :: frohl_model
   ! 1 to activate the computation of the Frohlich self-energy
   ! to treat the q-->0 divergence and accelerate convergence in polar semiconductors.

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

  !integer,allocatable:: indkk(:, :)
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
  ! (nz, 2, nbcalc_ks, natom3, nbsum, nq_k))
  ! Weights for the q-integration of 1 / (e1 - e2 \pm w_{q, nu} + i.eta)
  ! This array is initialized inside the (ikcalc, spin) loop

  real(dp),allocatable :: deltaw_pm(:,:,:,:,:)
  ! (2, nbcalc_ks, natom3, nbsum, nq_k))
  ! Weights for the q-integration of the two delta (abs/emission) if imag_only
  ! This array is initialized inside the (ikcalc, spin) loop

  real(dp),allocatable :: wrmesh_b(:,:)
  ! wrmesh_b(nwr, max_nbcalc)
  ! Frequency mesh along the real axis (Ha units) used for the different bands
  ! Each mesh is **centered** on the corresponding KS energy.
  ! This array depends on (ikcalc, spin)

  !real(dp),allocatable :: frohl_gkq2(:,:)
  ! frohl_gkq2(max_nbcalc, natom3)
  ! Stores the long-range Frohlich matrix element (squared)

  real(dp), allocatable :: qvers_cart(:,:)
   ! qvers_cart(3, angl_size)
   ! For each point of the angular mesh, gives the Cartesian coordinates
   ! of the corresponding point on an unitary sphere (Frohlich self-energy)

  real(dp), allocatable :: angwgth(:)
   ! angwgth(angl_size)
   ! For each point of the angular mesh, gives the weight
   ! of the corresponding point on an unitary sphere (Frohlich self-energy)

  complex(dpc),allocatable :: vals_e0ks(:,:)
   ! vals_e0ks(ntemp, max_nbcalc)
   ! Sigma_eph(omega=eKS, kT, band) for fixed (kcalc, spin).

  complex(dpc),allocatable :: frohl_vals_e0ks(:,:)
   ! frohl_vals_e0ks(ntemp, max_nbcalc)
   ! Sigma_frohl(omega=eKS, kT, band) for fixed (kcalc, spin).

  complex(dpc),allocatable :: dvals_de0ks(:,:)
   ! dvals_de0ks(ntemp, max_nbcalc) for fixed (kcalc, spin)
   ! d Re Sigma_eph(omega, kT, band, kcalc, spin) / d omega (omega=eKS)

  complex(dpc),allocatable :: frohl_dvals_de0ks(:,:)
   ! frohl_dvals_de0ks(ntemp, max_nbcalc) for fixed (kcalc, spin)
   ! d Re Sigma_frohl(omega, kT, band, kcalc, spin) / d omega (omega=eKS)

  real(dp),allocatable :: dw_vals(:,:)
  !  dw_vals(ntemp, max_nbcalc) for fixed (kcalc, spin)
  !  Debye-Waller term (static).

  !real(dp),allocatable :: qpoms_enes(:,:)
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

  complex(dpc),allocatable :: frohl_vals_wr(:,:,:)
   ! frohl_vals_wr(nwr, ntemp, max_nbcalc)
   ! Sigma_frohl(omega, kT, band) for given (k, spin).
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
   ! Quantities needed to compute generalized Eliashberg functions  (gkq2/Fan-Migdal/DW terms)
   ! This array depends on (ikcalc, spin)
   ! NB: q-weights for integration are not included.

  ! TODO: Can be removed now.
  real(dp),allocatable :: gfw_vals(:,:,:)
   !gfw_vals(gfw_nomega, 3, max_nbcalc)
   ! Generalized Eliashberg function a2F_{n,k,spin}(w)
   ! 1:   gkk^2 with delta(en - em)
   ! 2:3 (Fan-Migdal/DW contribution)
   ! This array depends on (ikcalc, spin)

  type(ephwg_t) :: ephwg
   ! This object compute the weights for the BZ integration in q-space if qint_method > 0

  type(skw_t) :: frohl_skw
   ! Star-function interpolator.
   ! Used to interpolate \epsilon_{nk+q) in the numerical integration of the Frohlich self-energy.

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


 real(dp),private,parameter :: EPH_WTOL = tol6
 ! Tolerance for phonon frequencies.

!----------------------------------------------------------------------

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
 integer,parameter :: dummy_npw=1,tim_getgh1c=1,berryopt0=0,timrev0=0
 integer,parameter :: useylmgr=0,useylmgr1=0,master=0,ndat1=1,nz=1
 integer,parameter :: sppoldbl1=1,timrev1=1
 integer :: my_rank,mband,my_minb,my_maxb,nsppol,nkpt,iq_ibz
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,nprocs
 integer :: ibsum_kq,ib_k,band_ks,num_smallw,ibsum,ii,jj,im,in
 integer :: idir,ipert,ip1,ip2,idir1,ipert1,idir2,ipert2
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq
 integer :: iq_ibz_fine,ikq_ibz_fine,ikq_bz_fine,iq_bz_fine
 integer :: spin,istwf_k,istwf_kq,istwf_kqirr,npw_k,npw_kq,npw_kqirr
 integer :: mpw,ierr,it
 integer :: n1,n2,n3,n4,n5,n6,nspden,nu
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnl1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,nq
 integer :: nbcalc_ks,nbsum,bstart_ks,ikcalc,bstart,bstop,my_bstart,my_bstop,iatom
 real(dp),parameter :: tol_enediff=0.001_dp*eV_Ha
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,gflops_all,cpu_ks,wall_ks,gflops_ks,cpu_dw,wall_dw,gflops_dw
 real(dp) :: cpu_setk, wall_setk, gflops_setk
 real(dp) :: ecut,eshift,dotr,doti,dksqmax,weigth_q,rfact,gmod2,hmod2,ediff,weight, inv_qepsq
 real(dp) :: elow,ehigh,wmax
 complex(dpc) :: cfact,dka,dkap,dkpa,dkpap,cplx_ediff, cnum
 logical :: isirr_k,isirr_kq,gen_eigenpb,isqzero
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(sigmaph_t) :: sigma
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0_kq(3),dummy_gvec(3,dummy_npw)
 integer :: work_ngfft(18),gmax(3), indkk_kq(1,6) !g0ibz_kq(3),
 integer,allocatable :: gtmp(:,:),kg_k(:,:),kg_kq(:,:),nband(:,:),distrib_bq(:,:)
 integer,allocatable :: indq2dvdb(:,:),wfd_istwfk(:),iqk2dvdb(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),qpt_cart(3),phfrq(3*cryst%natom)
 real(dp) :: wqnu,nqnu,gkq2,eig0nk,eig0mk,eig0mkq,f_mkq
 real(dp),allocatable :: displ_cart(:,:,:,:),displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:),gdw2_mn(:,:),gkq0_atm(:,:,:,:),gkq2_lr(:,:)
 complex(dpc),allocatable :: tpp_red(:,:) !,zvals(:,:)
 complex(dpc) :: cdd(3)
 real(dp),allocatable :: bra_kq(:,:),kets_k(:,:,:),h1kets_kq(:,:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnl1(:,:),work(:,:,:,:)
 real(dp),allocatable ::  gs1c(:,:),nqnu_tlist(:),dt_weights(:,:),dargs(:)
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

 ! Notes about MPI version.
 ! If eph_task == -4:
 !    Loops are MPI parallelized over bands so that we can distribute wavefunctions over nband.
 !
 ! If eph_task == -4:
 !    Loops are MPI parallelized over q-points
 !    wavefunctions are not distributed but only states between my_bstart
 !    and my_bstop are allocated and read from file

 if (sigma%imag_only) then
   ! Compute the min/max KS energy to be included in the imaginary part.
   ! ifc%omega_minmax(2) comes froms the coarse Q-mesh of the DDB so increase it by 10%.
   ! Also take into account Lorentzian shape if zcut is used.
   elow = huge(one); ehigh = - huge(one)
   wmax = 1.1_dp * ifc%omega_minmax(2) + five * dtset%zcut
   do ikcalc=1,sigma%nkcalc
     ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
     do spin=1,sigma%nsppol
       bstart = sigma%bstart_ks(ikcalc,spin)
       bstop = sigma%bstart_ks(ikcalc,spin) + sigma%nbcalc_ks(ikcalc,spin) - 1
       ehigh = max(ehigh, maxval(ebands%eig(bstart:bstop, ik_ibz, spin)) + wmax)
       elow = min(elow, minval(ebands%eig(bstart:bstop, ik_ibz, spin)) - wmax)
     end do
   end do

   ! Select indices for energy window.
   call get_bands_from_erange(ebands, elow, ehigh, my_bstart, my_bstop)
   call wrtout(std_out, sjoin("Allocating bands for imag part between bstart: ", itoa(my_bstart), &
     " and my_bstop:", itoa(my_bstop)))
   call wrtout(std_out, sjoin("elow:", ftoa(elow), "ehigh:", ftoa(ehigh), "[Ha]"))
   ABI_CHECK(my_bstart <= my_bstop, "my_bstart > my_bstop")
   bks_mask = .False.; bks_mask(my_bstart:my_bstop, : ,:) = .True.
   bks_mask = .True. ! TODO: Have to redefine nbsum and reshift band index in sum
   !sigma%nbsum
 else
   ! eph_task == 4 (Re + Im)
   call xmpi_split_work(sigma%nbsum, comm, my_bstart, my_bstop, msg, ierr)
   call wrtout(std_out, sjoin("Allocating and treating bands from bstart: ", itoa(my_bstart), &
     " up to my_bstop:", itoa(my_bstop)))
   if (my_bstart == sigma%nbsum + 1) then
     MSG_ERROR("sigmaph with idle processes should be tested! Decrease ncpus or increase nband")
   end if
   bks_mask = .False.; bks_mask(my_bstart:my_bstop, : ,:) = .True.
 end if

 ! Each node needs the wavefunctions for Sigma_{nk}
 do spin=1,sigma%nsppol
   do ikcalc=1,sigma%nkcalc
     ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
     bstart = sigma%bstart_ks(ikcalc,spin)
     bstop = bstart + sigma%nbcalc_ks(ikcalc,spin) - 1
     bks_mask(bstart:bstop, ik_ibz, spin) = .True.
   end do
 end do
 !bks_mask = .True. ! Uncomment this line to have all states on each MPI rank.

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
 ! This is also the reason why we reorder the q-points in ibz_k to pack the points in *shells* to minimise cache misses.
 call dvdb_set_qcache_mb(dvdb, dtset%dvdb_qcache_mb)
 call dvdb_print(dvdb, prtvol=dtset%prtvol)
 call dvdb_qcache_read(dvdb, nfftf, ngfftf, comm)

 ABI_MALLOC(displ_cart, (2,3,cryst%natom,3*cryst%natom))
 ABI_MALLOC(displ_red, (2,3,cryst%natom,3*cryst%natom))

 ! Loop over k-points in Sigma_nk.
 do ikcalc=1,sigma%nkcalc
   call cwtime(cpu_ks, wall_ks, gflops_ks, "start")
   kk = sigma%kcalc(:, ikcalc)

   ! Find IBZ(k) for q-point integration.
   call cwtime(cpu_setk, wall_setk, gflops_setk, "start")
   call sigmaph_setup_kcalc(sigma, cryst, ikcalc, dtset%prtvol)
   call wrtout(std_out, sjoin(ch10, repeat("=", 92)))
   msg = sjoin("[", itoa(ikcalc), "/", itoa(sigma%nqibz_k), "]")
   call wrtout(std_out, sjoin("Computing self-energy matrix elements for k-point:", ktoa(kk), msg))
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

   ABI_MALLOC(iqk2dvdb, (sigma%nqibz_k, 6))
   call listkk(dksqmax, cryst%gmet, iqk2dvdb, dvdb%qpts, sigma%qibz_k, dvdb%nqpt, sigma%nqibz_k, cryst%nsym, &
        1, cryst%symafm, cryst%symrec, timrev1, use_symrec=.True.)
   if (dksqmax > tol12) then
     write(msg, '(a,es16.6,2a)' )&
       "At least one of the q points could not be generated from a symmetrical one in the DVDB. dksqmax: ",dksqmax, ch10,&
       'Action: check your DVDB file and use eph_task to interpolate the potentials on a denser q-mesh.'
     MSG_ERROR(msg)
   end if
   do iq_ibz=1,sigma%nqibz_k
     indq2dvdb(:, iq_ibz) = iqk2dvdb(iq_ibz, :)
   end do
   ABI_FREE(iqk2dvdb)

   ! Allocate PW-arrays. Note mpw in kg_kq
   ABI_MALLOC(kg_k, (3, npw_k))
   kg_k = wfd%kdata(ik_ibz)%kg_k
   ABI_MALLOC(kg_kq, (3, mpw))

   ! Spherical Harmonics for useylm==1.
   ABI_MALLOC(ylm_k,(mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylm_kq,(mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

   call cwtime(cpu_setk, wall_setk, gflops_setk, "stop")
   call wrtout(std_out, sjoin("Setup kcalc completed. cpu-time:", sec2str(cpu_setk), ",wall-time:", sec2str(wall_setk)))

   do spin=1,nsppol
     ! Bands in Sigma_nk to compute and number of bands in sum over states.
     bstart_ks = sigma%bstart_ks(ikcalc, spin)
     nbcalc_ks = sigma%nbcalc_ks(ikcalc, spin)
     nbsum = sigma%nbsum

     ! Zero self-energy matrix elements. Build frequency mesh for nk states.
     sigma%vals_e0ks = zero; sigma%dvals_de0ks = zero; sigma%dw_vals = zero

     ! Prepare computation of Sigma_{nk}(w) and spectral function.
     if (sigma%nwr > 0) then
       sigma%vals_wr = zero
       do ib_k=1,nbcalc_ks
         band_ks = ib_k + bstart_ks - 1
         ! Build linear mesh **centered** on KS energy.
         eig0nk = ebands%eig(band_ks, ik_ibz, spin) - sigma%wr_step * (sigma%nwr / 2)
         sigma%wrmesh_b(:,ib_k) = arth(eig0nk, sigma%wr_step, sigma%nwr)
       end do
     end if

     ! Prepare Eliasberg function.
     if (sigma%gfw_nomega > 0) then
       if (allocated(sigma%gf_nnuq)) then
         ABI_FREE(sigma%gf_nnuq)
       end if
       ABI_CALLOC(sigma%gf_nnuq, (nbcalc_ks, natom3, sigma%nqibz_k, 3))
     end if

     ! Allocate eph matrix elements.
     ABI_MALLOC(gkq_atm, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkq_nu, (2, nbcalc_ks, natom3))

     ! Arrays for Debye-Waller
     if (.not. sigma%imag_only) then
       ABI_STAT_ALLOCATE(gkq0_atm, (2, nbcalc_ks, nbsum, natom3), ierr)
       ABI_CHECK(ierr == 0, "oom in gkq0_atm")
       gkq0_atm = zero
     end if

     ! Compute self-energy with Frohlich model for the gkq to handle divergence for q-->0
     ! and improve convergence with respect to nqpt.
     ABI_CALLOC(gkq2_lr, (nbcalc_ks, natom3))
     if (sigma%frohl_model /= 0) then
       call eval_sigfrohl(sigma, cryst, ifc, ebands, ikcalc, spin, comm)
     end if

     ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
     ! and save KS energies in sigma%e0vals
     ! TODO: symmetrize them if kk is not irred
     ABI_MALLOC(kets_k, (2, npw_k*nspinor, nbcalc_ks))
     ABI_MALLOC(sigma%e0vals, (nbcalc_ks))
     do ib_k=1,nbcalc_ks
       band_ks = ib_k + bstart_ks - 1
       call wfd_copy_cg(wfd, band_ks, ik_ibz, spin, kets_k(1,1,ib_k))
       sigma%e0vals(ib_k) = ebands%eig(band_ks, ik_ibz, spin)
     end do

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Distribute q-points and bands.
     ! The distribution must be consistent with the WF distribution done with bks_mask
     ABI_MALLOC(distrib_bq, (nbsum, sigma%nqibz_k))
     distrib_bq = -1

     if (dtset%eph_task == 4) then
       distrib_bq(my_bstart:my_bstop, :) = my_rank
     else if (dtset%eph_task == -4) then
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
     else
       MSG_ERROR("Invalid eph_task")
     end if

     !ABI_MALLOC(zvals, (nz, nbcalc_ks))

     if (sigma%qint_method > 0) then
       ! Weights for Re-Im with i.eta shift.
       ! FIXME: This part is broken now. Lot of memory allocated here!
       ABI_MALLOC(sigma%cweights, (nz, 2, nbcalc_ks, natom3, nbsum, sigma%ephwg%nq_k))
       ! Weights for Im (tethraedron, eta --> 0)
       ABI_MALLOC(sigma%deltaw_pm, (2 ,nbcalc_ks, natom3, nbsum, sigma%ephwg%nq_k))

       ! Map sigma%eph_doublegrid%dense -> ephwg%lgk%ibz
       if (sigma%use_doublegrid) then
         call eph_double_grid_bz2ibz(sigma%eph_doublegrid, sigma%ephwg%lgk%ibz, sigma%ephwg%lgk%nibz,&
                                     sigma%ephwg%lgk%symrec_lg, sigma%ephwg%lgk%nsym_lg, &
                                     sigma%eph_doublegrid%bz2lgkibz)
       endif

       ! Precompute the weights for tetrahedron
       call sigmaph_get_all_qweights(sigma,cryst,ebands,spin,ikcalc,comm)
     endif

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
         MSG_ERROR(sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB"))
       end if

       ! TODO: Make sure that symmetries in Q-space are preserved.
       ! Avoid fourq if qpt is in ddb
       ! Examine the symmetries of the q wavevector
       !call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

       ! Get phonon frequencies and displacements in reduced coordinates for this q-point
       call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

       if (sigma%frohl_model /= 0) then
         ! TODO: Recheck this part
         qpt_cart = matmul(cryst%rprimd, qpt)
         inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))

         ! Compute gkq_{LR} without (i 4pi/ucvol)
         do nu=1,natom3
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
           cnum = zero
           do iatom=1,cryst%natom
             cdd = cmplx(displ_cart(1,:, iatom, nu), displ_cart(2,:, iatom, nu)) * exp(-j_dpc * two_pi * cryst%xred(:, iatom))
             !cdd = cdd * exp(-(q/sigma%qdamp)**2)
             cnum = cnum + dot_product(qpt_cart, matmul(ifc%zeff(:, :, iatom), cdd))
           end do
           do ib_k=1,nbcalc_ks
             gkq2_lr(ib_k, nu) = (real(cnum) ** 2 + aimag(cnum) ** 2) / (two * wqnu) * inv_qepsq
           end do
         end do
       end if

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

       ! Double grid stuff
       if (sigma%use_doublegrid) then
         call eph_double_grid_get_mapping(sigma%eph_doublegrid,kk,kq,qpt)
       endif

       ! Map q to qibz for tetrahedron
       if (sigma%qint_method > 0) then
         if (.not.sigma%use_doublegrid) then
           iq_ibz_fine = iq_ibz
           if (sigma%symsigma == 0) iq_ibz_fine = lgroup_find_ibzimage(sigma%ephwg%lgk, qpt)
           ABI_CHECK(iq_ibz_fine /= -1, sjoin("Cannot find q-point in IBZ(k)", ktoa(qpt)))
           if (sigma%symsigma == 1) then
             ABI_CHECK(all(abs(sigma%qibz_k(:, iq_ibz_fine) - sigma%ephwg%lgk%ibz(:, iq_ibz_fine)) < tol12), "Mismatch in qpoints.")
           end if
         endif
       end if

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
           band_ks = ib_k + bstart_ks - 1
           eig0nk = ebands%eig(band_ks, ik_ibz, spin)
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
         !band_sum =
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
             gkq_atm(:,ib_k,ipc) = [dotr, doti]
           end do
         end do

         ! Get gkk(kcalc, q, nu)
         call gkknu_from_atm(1, nbcalc_ks, 1, natom, gkq_atm, phfrq, displ_red, gkq_nu, num_smallw)

         ! DEBUG: Write gkk matrix elements to file. DO NOT USE IN PRODUCTION
         ! nctkarr_t("gkk", "dp", "two, max_nbcalc, nbsum, natom3, nqbz, nkcalc, nsppol"), &
         !ii = nf90_put_var(sigma%ncid, nctk_idname(sigma%ncid, "gkk"), gkq_nu, &
         !  start=[1, 1, ibsum, 1, iq_ibz, ikcalc, spin], count=[2, nbcalc_ks, 1, natom3, 1, 1, 1])
         !NCF_CHECK(ii)
         !ii = nf90_put_var(sigma%ncid, nctk_idname(ncid, "frohl_gkq2"), gkq_nu, &
         !  start=[1, 1, ibsum, 1, iq_ibz, ikcalc, spin], count=[])
         !NCF_CHECK(ii)

         ! Save data for Debye-Waller computation (performed outside the q-loop)
         ! gkq_nu(2, nbcalc_ks, nbsum, natom3)
         if (isqzero .and. .not. sigma%imag_only) then
           gkq0_atm(:, :, ibsum_kq, :) = gkq_atm
         end if

         ! Accumulate contribution to self-energy
         eig0mkq = ebands%eig(ibsum_kq,ikq_ibz,spin)
         ! q-weigths for naive integration
         weigth_q = sigma%wtq_k(iq_ibz)

         do nu=1,natom3
           ! Ignore acoustic or unstable modes.
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle

           ! For each band in Sigma_{bk}
           do ib_k=1,nbcalc_ks
             band_ks = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band_ks, ik_ibz, spin)
             gkq2 = weigth_q * (gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2)

             ! Handle intra-band term in polar materials
             ! TODO: be careful when reshifting nbsum in imag_only case.
             ! DFPT gkk for q == 0 should not diverge (should plot gkq2(q) though)
             if (sigma%frohl_model /= 0 .and. ibsum_kq == band_ks .and. .not. isqzero) then
               gkq2 = gkq2 - weigth_q * gkq2_lr(ib_k, nu)
             end if

             ! Optionally, accumulate contribution to Eliashberg functions
             ediff = eig0nk - eig0mkq
             !if (abs(cplx_ediff) < EPH_WTOL) cplx_ediff = cplx_ediff + sigma%ieta
             if (sigma%gfw_nomega > 0 .and. abs(eig0nk - eig0mkq) > EPH_WTOL) then
               ! eph strength with delta(en - emkq)
               rfact = dirac_delta(eig0nk - eig0mkq, dtset%tsmear)
               sigma%gf_nnuq(ib_k, nu, iq_ibz, 1) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 1) + &
                    rfact * (gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2)

               ! Fan term.
               if (ediff > wqnu) then
                  rfact = one / ediff
               else
                 ! Non adiabatic regime --> Add complex shift.
                 ! Note however that the expression for the Eliashberg function relies on adiabaticity.
                 rfact = real(one/(ediff + sigma%ieta))
               end if
               sigma%gf_nnuq(ib_k, nu, iq_ibz, 2) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 2) + &
                  (gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2) * rfact
             end if

             do it=1,sigma%ntemp
               ! Compute occ for this T (note mu_e(it) Fermi level)
               nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)
               f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

               ! Here we have to handle 3 different logical values
               ! leading to 9 different cases:
               !
               ! qint_method         0      1
               !   use_doublegrid   .true. .false.
               !     imag_only      .true. .false.
               !
               ! we will write this with nested conditionals using the order above

               if (sigma%qint_method == 0) then
                 ! zcut mode
                 if (sigma%use_doublegrid) then
                   cfact = zero
                   do jj=1,sigma%eph_doublegrid%ndiv
                     ! Double Grid shared points weights
                     ikq_bz_fine  = sigma%eph_doublegrid%mapping(2, jj)
                     weight = sigma%eph_doublegrid%weights_dense(ikq_bz_fine)

                     ! Electronic eigenvalue
                     ikq_ibz_fine = sigma%eph_doublegrid%mapping(5, jj)
                     eig0mkq = sigma%eph_doublegrid%ebands_dense%eig(ibsum_kq,ikq_ibz_fine,spin)
                     f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

                     ! Phonon frequency
                     iq_ibz_fine = sigma%eph_doublegrid%mapping(6, jj)
                     wqnu = sigma%ephwg%phfrq_ibz(iq_ibz_fine,nu)
                     !if (wqnu < EPH_WTOL) cycle
                     nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)

                     cfact = cfact + &
                            ((nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                             (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta) ) * weight
                   enddo
                 else
                   cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                            (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)
                 endif
                 if (sigma%imag_only) then
                   sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * j_dpc * aimag(cfact)
                 else
                   sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * cfact
                 end if

               else
                 ! Tetrahedron method
                 if (sigma%use_doublegrid) then
                   do jj=1,sigma%eph_doublegrid%ndiv
                     ! Double Grid shared points weights
                     ikq_bz_fine  = sigma%eph_doublegrid%mapping(2, jj)
                     weight = sigma%eph_doublegrid%weights_dense(ikq_bz_fine)

                     ! Electronic eigenvalue
                     ikq_ibz_fine = sigma%eph_doublegrid%mapping(5, jj)
                     eig0mkq = sigma%eph_doublegrid%ebands_dense%eig(ibsum_kq,ikq_ibz_fine,spin)
                     f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

                     ! Phonon frequency
                     iq_ibz_fine = sigma%eph_doublegrid%mapping(6, jj)
                     wqnu = sigma%ephwg%phfrq_ibz(iq_ibz_fine,nu)
                     !if (wqnu < EPH_WTOL) cycle
                     nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)

                     ! Map bz to lgk
                     iq_bz_fine = sigma%eph_doublegrid%mapping(3,jj)
                     iq_ibz_fine = sigma%eph_doublegrid%bz2lgkibz(iq_bz_fine)

                     if (sigma%imag_only) then
                       ! note pi factor (Sokhotski–Plemelj theorem)
                       sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * j_dpc * pi * ( &
                         (nqnu + f_mkq      ) * sigma%deltaw_pm(1, ib_k, nu, ibsum_kq, iq_ibz_fine) +  &
                         (nqnu - f_mkq + one) * sigma%deltaw_pm(2, ib_k, nu, ibsum_kq, iq_ibz_fine) ) * weight
                     else
                       sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * ( &
                         (nqnu + f_mkq      ) * sigma%cweights(1, 1, ib_k, nu, ibsum_kq, iq_ibz_fine) +  &
                         (nqnu - f_mkq + one) * sigma%cweights(1, 2, ib_k, nu, ibsum_kq, iq_ibz_fine) ) * weight
                     end if
                   end do
                 else
                   if (sigma%imag_only) then
                     sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * j_dpc * pi * ( &
                       (nqnu + f_mkq      ) * sigma%deltaw_pm(1, ib_k, nu, ibsum_kq, iq_ibz_fine) +  &
                       (nqnu - f_mkq + one) * sigma%deltaw_pm(2, ib_k, nu, ibsum_kq, iq_ibz_fine) )
                   else
                     sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * ( &
                       (nqnu + f_mkq      ) * sigma%cweights(1, 1, ib_k, nu, ibsum_kq, iq_ibz_fine) +  &
                       (nqnu - f_mkq + one) * sigma%cweights(1, 2, ib_k, nu, ibsum_kq, iq_ibz_fine) )
                   endif
                 end if
               end if

               ! Derivative of sigma
               ! TODO: should calculate this with the double grid as well
               if (.not. sigma%imag_only) then
                 !if (sigma%qint_method == 1) then
                   ! Have to rescale gkq2 before computing derivative (HM: why?)
                 !  gkq2 = gkq2 * sigma%wtq_k(iq_ibz)
                 !end if

                 ! Accumulate d(Re Sigma) / dw(w=eKS) for state ib_k
                 !cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
                 !            (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
                 gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
                 hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
                 rfact = (nqnu + f_mkq      ) * (-gmod2 + aimag(sigma%ieta)**2) / (gmod2 + aimag(sigma%ieta)**2) ** 2 + &
                         (nqnu - f_mkq + one) * (-hmod2 + aimag(sigma%ieta)**2) / (hmod2 + aimag(sigma%ieta)**2) ** 2
                 sigma%dvals_de0ks(it, ib_k) = sigma%dvals_de0ks(it, ib_k) + gkq2 * rfact
                 !cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                 !         (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)
                 !sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2 * cfact

                 !cfact = (eig0nk - eig0mkq + wqnu + sigma%ieta)
                 !gmod2 = cfact * dconjg(cfact)
                 !cfact = (eig0nk - eig0mkq - wqnu + sigma%ieta)
                 !hmod2 = cfact * dconjg(cfact)
                 !sigma%dvals_de0ks(it, ib_k) = sigma%dvals_de0ks(it, ib_k) + gkq2 * ( &
                 !  (nqnu + f_mkq)        * (gmod2 - two * (eig0nk - eig0mkq + wqnu) ** 2) / gmod2 ** 2 + &
                 !  (nqnu - f_mkq + one)  * (hmod2 - two * (eig0nk - eig0mkq - wqnu) ** 2) / hmod2 ** 2   &
                 !)

                 ! Accumulate Sigma(w) for state ib_k if spectral function is wanted.
                 ! TODO: weigths
                 if (sigma%nwr > 0) then
                   cfact_wr(:) = (nqnu + f_mkq      ) / (sigma%wrmesh_b(:,ib_k) - eig0mkq + wqnu + sigma%ieta) + &
                                 (nqnu - f_mkq + one) / (sigma%wrmesh_b(:,ib_k) - eig0mkq - wqnu + sigma%ieta)
                   sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + gkq2 * cfact_wr(:)
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

       if (sigma%nqibz_k < 1000 .or. (sigma%nqibz_k > 1000 .and. mod(iq_ibz, 200) == 0) .or. iq_ibz <= nprocs) then
         call cwtime(cpu, wall, gflops, "stop")
         write(msg,'(4(a,i0),2(a,f8.2))') "k-point [",ikcalc,"/",sigma%nkcalc, &
                                          "] q-point [",iq_ibz,"/",sigma%nqibz_k,"] completed. cpu:",cpu,", wall:",wall
         call wrtout(std_out, msg, do_flush=.True.)
       end if
     end do ! iq_ibz (sum over q-points)

     if (sigma%qint_method > 0) then
       ABI_FREE(sigma%deltaw_pm)
       ABI_FREE(sigma%cweights)
     end if

     ! Print cache stats. The first k-point is expected to have lots of misses
     ! especially if it's the Gamma point and symsigma = 1.
     if (dvdb%qcache_size > 0 .and. dvdb%qcache_stats(1) /= 0) then
       write(std_out, "(a)")"Qcache stats"
       write(std_out, "(a,i0)")"Total Number of calls: ", dvdb%qcache_stats(1)
       write(std_out, "(a,i0,2x,f5.1,a)")&
         "Cache hit: ", dvdb%qcache_stats(2), (100.0_dp * dvdb%qcache_stats(2)) / dvdb%qcache_stats(1), "%"
       write(std_out, "(a,i0,2x,f5.1,a)")&
         "Cache miss: ", dvdb%qcache_stats(3), (100.0_dp * dvdb%qcache_stats(3)) / dvdb%qcache_stats(1), "%"
       dvdb%qcache_stats = 0
     end if

     ABI_FREE(sigma%e0vals)
     !ABI_FREE(zvals)
     ABI_FREE(kets_k)
     ABI_FREE(gkq_atm)
     ABI_FREE(gkq_nu)
     ABI_FREE(gkq2_lr)
     ABI_FREE(distrib_bq)

     ! =========================
     ! Compute Debye-Waller term
     ! =========================
     if (.not. sigma%imag_only) then
       call wrtout(std_out, "Computing Debye-Waller term...")
       call cwtime(cpu_dw, wall_dw, gflops_dw, "start")
       call xmpi_sum(gkq0_atm, comm, ierr)

       ABI_MALLOC(gdw2_mn, (nbsum, nbcalc_ks))
       ABI_MALLOC(tpp_red, (natom3, natom3))

       ! Integral over IBZ(k).
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
         call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

         ! Sum over modes for this q-point.
         do nu=1,natom3
           ! Ignore acoustic or unstable modes.
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle

           ! Compute T_pp'(q,nu) matrix in cartesian coordinates.
           do ip2=1,natom3
             idir2 = mod(ip2-1, 3) + 1; ipert2 = (ip2 - idir2) / 3 + 1
             do ip1=1,natom3
               idir1 = mod(ip1-1, 3) + 1; ipert1 = (ip1 - idir1) / 3 + 1
               ! (k,a) (k,a')* + (k',a) (k',a')*
               !dka   = dcmplx(displ_cart(1, idir1, ipert1, nu), displ_cart(2, idir1, ipert1, nu))
               !dkap  = dcmplx(displ_cart(1, idir2, ipert1, nu), displ_cart(2, idir2, ipert1, nu))
               !dkpa  = dcmplx(displ_cart(1, idir1, ipert2, nu), displ_cart(2, idir1, ipert2, nu))
               !dkpap = dcmplx(displ_cart(1, idir2, ipert2, nu), displ_cart(2, idir2, ipert2, nu))
               !tpp(ip1,ip2) = dka * dconjg(dkap) + dkpa * dconjg(dkpap)

               dka   = dcmplx(displ_red(1, idir1, ipert1, nu), displ_red(2, idir1, ipert1, nu))
               dkap  = dcmplx(displ_red(1, idir2, ipert1, nu), displ_red(2, idir2, ipert1, nu))
               dkpa  = dcmplx(displ_red(1, idir1, ipert2, nu), displ_red(2, idir1, ipert2, nu))
               dkpap = dcmplx(displ_red(1, idir2, ipert2, nu), displ_red(2, idir2, ipert2, nu))
               tpp_red(ip1,ip2) = dka * dconjg(dkap) + dkpa * dconjg(dkpap)
             end do
           end do

           ! Get phonon occupation for all temperatures.
           nqnu_tlist = occ_be(wqnu, sigma%kTmesh(:), zero)

           ! Sum over bands and add (static) DW contribution for the different temperatures.
           do ibsum=1,nbsum
             do ib_k=1,nbcalc_ks
               ! Compute DW term following XG paper. Check prefactor.
               ! gkq0_atm(2, nbcalc_ks, nbsum, natom3)
               ! previous version
               gdw2_mn(ibsum, ib_k) = zero
               do ip2=1,natom3
                 do ip1=1,natom3
                   cfact = ( &
                     + gkq0_atm(1, ib_k, ibsum, ip1) * gkq0_atm(1, ib_k, ibsum, ip2) &
                     + gkq0_atm(2, ib_k, ibsum, ip1) * gkq0_atm(2, ib_k, ibsum, ip2) &
                     + gkq0_atm(1, ib_k, ibsum, ip2) * gkq0_atm(1, ib_k, ibsum, ip1) &
                     + gkq0_atm(2, ib_k, ibsum, ip2) * gkq0_atm(2, ib_k, ibsum, ip1) &
                   )

                   gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) + real(tpp_red(ip1,ip2) * cfact)
                 end do
               end do

               gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) /  (four * two * wqnu)
             end do ! ibsum
           end do ! ib_k

           do ibsum=1,nbsum
             eig0mk = ebands%eig(ibsum, ik_ibz, spin)

             do ib_k=1,nbcalc_ks
               band_ks = ib_k + bstart_ks - 1
               eig0nk = ebands%eig(band_ks, ik_ibz, spin)
               ! Handle n == m and degenerate states (either ignore or add broadening)
               cplx_ediff = (eig0nk - eig0mk)
               if (abs(cplx_ediff) < EPH_WTOL) cycle
               !if (abs(cplx_ediff) < EPH_WTOL) cplx_ediff = cplx_ediff + sigma%ieta

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
       ABI_FREE(tpp_red)
       ABI_FREE(gkq0_atm)

       call cwtime(cpu_dw, wall_dw, gflops_dw, "stop")
       call wrtout(std_out, sjoin("DW completed. cpu-time:", sec2str(cpu_dw), ",wall-time:", &
          sec2str(wall_dw)), do_flush=.True.)
     end if ! not %imag_only

     if (sigma%gfw_nomega /= 0) then
       ! Compute Eliashberg function (useful but cost is not negligible.
       ! May need to deactivate this part for HTC.
       call wrtout(std_out, sjoin("Computing Eliashberg function with nomega:", itoa(sigma%gfw_nomega), &
           "Use prtphdos 0 to disable this part"))
       call cwtime(cpu, wall, gflops, "start")
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

       call cwtime(cpu, wall, gflops, "stop")
       call wrtout(std_out, sjoin("Eliashberg function completed. cpu-time:", sec2str(cpu), &
           ",wall time:", sec2str(wall)), do_flush=.True.)

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
     call sigmaph_gather_and_write(sigma, ebands, ikcalc, spin, dtset%prtvol, comm)
   end do ! spin

   ABI_FREE(indq2dvdb)
   ABI_FREE(kg_k)
   ABI_FREE(kg_kq)
   ABI_FREE(ylm_k)
   ABI_FREE(ylm_kq)
   ABI_FREE(ylmgr_kq)

   call cwtime(cpu_ks, wall_ks, gflops_ks, "stop")
   call wrtout(std_out, sjoin("Sigma_nk completed. cpu-time:", sec2str(cpu_ks), &
           ",wall-time:", sec2str(wall_ks)), do_flush=.True.)
 end do ! ikcalc

 call cwtime(cpu_all, wall_all, gflops_all, "stop")
 call wrtout(std_out, "Sigma_eph completed")
 call wrtout(std_out, sjoin("Total cpu-time:", sec2str(cpu_all), ", Total wall-time:", sec2str(wall_all), ch10))

 ! Free memory
 ABI_FREE(gvnl1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(nqnu_tlist)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 if (sigma%nwr > 0) then
   ABI_FREE(cfact_wr)
 end if

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
!!  nb1,nb2=Number of bands in gkq_atm matrix.
!!  nk=Number of k-points (usually 1)
!!  natom=Number of atoms.
!!  gkq_atm(2,nb1,nb2,3*natom)=EPH matrix elements in the atomic basis.
!!  phfrq(3*natom)=Phonon frequencies in Ha
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduced coordinates.
!!
!! OUTPUT
!!  gkq_nu(2,nb1,nb2,3*natom)=EPH matrix elements in the phonon-mode basis.
!!  num_smallw=Number of negative/too small frequencies that have been ignored
!!    by setting the corresponding gkq_nu to zero.
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkknu_from_atm(nb1, nb2, nk, natom, gkq_atm, phfrq, displ_red, gkq_nu, num_smallw)


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
 real(dp),intent(in) :: gkq_atm(2,nb1,nb2,nk,3*natom)
 real(dp),intent(out) :: gkq_nu(2,nb1,nb2,nk,3*natom)

!Local variables-------------------------
!scalars
 integer :: nu,ipc

! *************************************************************************

 gkq_nu = zero; num_smallw = 0

 ! Loop over phonon branches.
 do nu=1,3*natom
   ! Ignore negative or too small frequencies
   if (phfrq(nu) < EPH_WTOL) then
     num_smallw = num_smallw + 1; cycle
   end if

   ! Transform the gkk from (atom, reduced direction) basis to phonon mode representation
   do ipc=1,3*natom
     gkq_nu(1,:,:,:,nu) = gkq_nu(1,:,:,:,nu) &
       + gkq_atm(1,:,:,:,ipc) * displ_red(1,ipc,nu) &
       - gkq_atm(2,:,:,:,ipc) * displ_red(2,ipc,nu)
     gkq_nu(2,:,:,:,nu) = gkq_nu(2,:,:,:,nu) &
       + gkq_atm(1,:,:,:,ipc) * displ_red(2,ipc,nu) &
       + gkq_atm(2,:,:,:,ipc) * displ_red(1,ipc,nu)
   end do

   gkq_nu(:,:,:,:,nu) = gkq_nu(:,:,:,:,nu) / sqrt(two * phfrq(nu))
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
 integer :: onpw,ii,ipw,ierr,it,spin,gap_err,ikcalc,gw_qprange,bstop !,band_ks
 integer :: nk_found,ifo,jj,bstart,nbcount,sigma_nkbz
 integer :: isym_k, trev_k
 integer :: ip,npoints,skw_cplex
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 logical :: downsample
 real(dp),parameter :: spinmagntarget=-99.99_dp,tol_enediff=0.001_dp*eV_Ha
 character(len=500) :: wfk_fname_dense
 real(dp) :: dksqmax,ang,con,cos_phi,cos_theta,sin_phi,sin_theta,nelect
 character(len=500) :: msg
 logical :: changed,found,isirr_k
 type(ebands_t) :: tmp_ebands, ebands_dense
 type(gaps_t) :: gaps
!arrays
 integer :: intp_kptrlatt(3,3), g0_k(3), skw_band_block(2)
 integer :: qptrlatt(3,3),indkk_k(1,6),my_gmax(3),kpos(6),band_block(2),kptrlatt(3,3)
 integer :: val_indeces(ebands%nkpt, ebands%nsppol), intp_nshiftk
 integer,allocatable :: gtmp(:,:),degblock(:,:)
 real(dp):: params(4), my_shiftq(3,1),kk(3),kq(3),intp_shiftk(3)
 real(dp),allocatable :: sigma_wtk(:),sigma_kbz(:,:),th(:),wth(:)

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
 !dtset%freqspmin;
 new%nwr = dtset%nfreqsp
 new%wr_step = zero
 if (new%nwr > 0) then
   if (mod(new%nwr, 2) == 0) new%nwr = new%nwr + 1
   new%wr_step = two * eV_Ha / (new%nwr - 1)
   if (dtset%freqspmax /= zero) new%wr_step = dtset%freqspmax / (new%nwr - 1)
 end if

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

 ! TODO: Rename variable
 if (dtset%nkptgw /= 0) then
   ! Treat the k-points and bands specified in the input file.
   call wrtout(std_out, "Generating list of k-points for self-energy from kptgw and bdgw.")
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
          "For (k, s) ",ikcalc,spin," bdgw= ",dtset%bdgw(2,ikcalc,spin), " > mband = ",dtset%mband
         MSG_WARNING(msg)
       end if
     end do
   end do
   ABI_CHECK(ierr == 0, "Not enough bands in WFK file. See messages above. Aborting now.")

 else
   ! Use qp_range to select the interesting k-points and the corresponing bands.
   !
   !    0 --> Compute the QP corrections only for the fundamental and the direct gap.
   ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
   !           bands above and below the Fermi level.
   ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
   !          Include all occupied states and `num` empty states.

   gw_qprange = dtset%gw_qprange
   if (gap_err /=0 .and. gw_qprange == 0) then
     msg = "Problem while computing the fundamental and direct gap (likely metal). Will replace gw_qprange=0 with gw_qprange=1"
     MSG_WARNING(msg)
     gw_qprange = 1
   end if

   if (gw_qprange /= 0) then

     if (any(dtset%sigma_ngkpt /= 0)) then
        call wrtout(std_out, "Generating list of k-points for self-energy from sigma_nkpt and qprange.")
        ABI_CHECK(gw_qprange /= 0, "gw_qprange must be != 0")
        ! Get %kcalc from sigma_ngkpt
        kptrlatt = 0
        kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)
        !write(std_out,*)"kptrlatt", kptrlatt
        !write(std_out,*)"sigma_nshiftk", dtset%sigma_nshiftk
        !write(std_out,*)"sigma_shiftk", dtset%sigma_shiftk
        call kpts_ibz_from_kptrlatt(cryst, kptrlatt, dtset%kptopt, dtset%sigma_nshiftk, dtset%sigma_shiftk, &
          new%nkcalc, new%kcalc, sigma_wtk, sigma_nkbz, sigma_kbz)
        ABI_FREE(sigma_kbz)
        ABI_FREE(sigma_wtk)
     else
        ! Include all the k-points in the IBZ.
        ! Note that kcalc == ebands%kptns so we can use a single ik index in the loop over k-points.
        ! No need to map kcalc onto ebands%kptns.
        call wrtout(std_out, "nkptgw set to 0 ==> Automatic selection of k-points and bands for sigma_{nk}.")
        new%nkcalc = ebands%nkpt
        ABI_MALLOC(new%kcalc, (3, new%nkcalc))
        new%kcalc = ebands%kptns
     end if

     ABI_MALLOC(new%bstart_ks, (new%nkcalc, new%nsppol))
     ABI_MALLOC(new%nbcalc_ks, (new%nkcalc, new%nsppol))

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
     ! Include the direct and the fundamental KS gap.
     ! The main problem here is that kptgw and nkptgw do not depend on the spin and therefore
     ! we have compute the union of the k-points where the fundamental and the direct gaps are located.
     !
     ! Find the list of `interesting` kpoints.
     call wrtout(std_out, "qprange not specified in input --> Include direct and fundamental KS gap in Sigma_{nk}")
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

 ierr = 0
 do ikcalc=1,new%nkcalc
   kk = new%kcalc(:,ikcalc)
   call listkk(dksqmax,cryst%gmet,indkk_k,ebands%kptns,kk,ebands%nkpt,1,cryst%nsym,&
      sppoldbl1,cryst%symafm,cryst%symrel,new%timrev,use_symrec=.False.)

   new%kcalc2ibz(ikcalc, :) = indkk_k(1, :)

   ik_ibz = indkk_k(1,1); isym_k = indkk_k(1,2)
   trev_k = indkk_k(1, 6); g0_k = indkk_k(1, 3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !kk_ibz = ebands%kptns(:,ik_ibz)
   if (.not. isirr_k) then
     MSG_WARNING(sjoin("For the time being the k-point must be in the IBZ but got", ktoa(kk)))
     ierr = ierr + 1
   end if

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
                      "Action: Increase nband."
         ABI_CHECK(bstop <= dtset%nband(new%nkcalc*(spin-1)+ikcalc), msg)
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
 ABI_CHECK(ierr == 0, "Fatal error, kptgw list must be in the IBZ")

 ! Now we can finally compute max_nbcalc
 new%max_nbcalc = maxval(new%nbcalc_ks)

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
 ! Allow users to deactivate this part with prtphdos == 0
 new%gfw_nomega = 0
 if (dtset%prtphdos == 1) then
   ! TODO: Use phdos min/max?
   new%gfw_nomega = nint((ifc%omega_minmax(2) - ifc%omega_minmax(1) ) / dtset%ph_wstep) + 1
   ABI_MALLOC(new%gfw_mesh, (new%gfw_nomega))
   new%gfw_mesh = arth(ifc%omega_minmax(1), dtset%ph_wstep, new%gfw_nomega)
   ABI_MALLOC(new%gfw_vals, (new%gfw_nomega, 3, new%max_nbcalc))
 end if

 new%imag_only = .False.; if (dtset%eph_task == -4) new%imag_only = .True.
 ! TODO: Remove qint_method, use eph_intmeth or perhaps dtset%qint_method dtset%kint_method
 ! Decide default behaviour for Re-Im/Im
 new%qint_method = dtset%eph_intmeth - 1

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
 if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine == 0) .or.&
     (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /= 0) )  then

   wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)//'FINE'
   if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
     MSG_ERROR(msg)
   end if
   call wrtout(std_out,"EPH Interpolation: will read energies from: "//trim(wfk_fname_dense))

   ebands_dense = wfk_read_ebands(wfk_fname_dense, comm)

   !TODO add a check for consistency
   ! number of bands and kpoints (comensurability)
   ABI_CHECK(ebands_dense%mband == ebands%mband, 'Inconsistent number of bands for the fine and dense grid')

   new%use_doublegrid = .True.

 ! read bs_interpmult
 else if (any(dtset%bs_interp_kmult /= 0)) then

   call wrtout(std_out,"EPH Interpolation: will use star functions interpolation")
   ! Interpolate band energies with star-functions
   params = 0; params(1) = 1; params(2) = 5
   if (nint(dtset%einterp(1)) == 1) params = dtset%einterp
   write(std_out, "(a, 4(f5.2, 2x))")"SKW parameters for double-grid:", params

   !TODO: mband should be min of nband
   band_block = [1, ebands%mband]
   intp_kptrlatt(:,1) = [ebands%kptrlatt(1,1)*dtset%bs_interp_kmult(1), 0, 0]
   intp_kptrlatt(:,2) = [0, ebands%kptrlatt(2,2)*dtset%bs_interp_kmult(2), 0]
   intp_kptrlatt(:,3) = [0, 0, ebands%kptrlatt(3,3)*dtset%bs_interp_kmult(3)]

   intp_nshiftk = 1
   intp_shiftk = zero
   ebands_dense = ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt, &
                                      intp_nshiftk, intp_shiftk, band_block, comm)
   new%use_doublegrid = .True.
 end if

 ! bstart and new%bsum select the band range.
 bstart = 1
 if (new%qint_method > 0) then
   if (new%use_doublegrid) then
     ! Double-grid technique from ab-initio energies or star-function interpolation.
     new%ephwg          = ephwg_from_ebands(cryst, ifc, ebands_dense, bstart, new%nbsum, comm)
     new%eph_doublegrid = eph_double_grid_new(cryst, ebands_dense, ebands%kptrlatt, ebands_dense%kptrlatt)
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
 new%mu_e(:) = ebands%fermie

 if (dtset%eph_fermie == zero) then
   if (new%use_doublegrid) then
     call ebands_copy(ebands_dense, tmp_ebands)
   else
     call ebands_copy(ebands, tmp_ebands)
   endif

   ! We only need mu_e so MPI parallelize the T-loop.
   new%mu_e = zero
   do it=1,new%ntemp
     if (mod(it, nprocs) /= my_rank) cycle
     ! Use Fermi-Dirac occopt
     call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(it), spinmagntarget, dtset%prtvol)
     call ebands_set_nelect(tmp_ebands, ebands%nelect, dtset%spinmagntarget, msg)
     new%mu_e(it) = tmp_ebands%fermie
     !
     ! Check that the total number of electrons is correct
     ! This is to trigger problems as the routines that calculate the occupations in ebands_set_nelect
     ! are different from the occ_fd that will be used in the rest of the subroutine
     !
     nelect = ebands_calc_nelect(tmp_ebands, new%kTmesh(it), new%mu_e(it))

     if (abs(nelect - ebands%nelect) > tol6) then
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
         MSG_ERROR(msg)
       end if
     end if
   end do ! it

   call ebands_free(tmp_ebands)
   call xmpi_sum(new%mu_e, comm, ierr)
 endif

 call ebands_free(ebands_dense)

 ! Prepare computation of Frohlich self-energy
 new%frohl_model = 0
 new%frohl_model = dtset%useria
 if (new%frohl_model /= 0) then
   ! Init parameters for numerical integration inside sphere.
   ! Set sphere radius to a fraction of the smallest reciprocal lattice vector.
   !ABI_MALLOC(new%frohl_gkq2, (new%max_nbcalc, natom3))
   new%ntheta = 10
   !new%ntheta = dtser%efmas_ntheta
   new%nphi = 2 * new%ntheta
   new%qrad = tol6
   !new%qrad = half * min(norm2(cryst%gprimd(:, 1)), norm2(cryst%gprimd(:, 2)), norm2(cryst%gprimd(:, 3)))
   new%qrad = new%qrad / 2.0_dp
   !new%qdamp = one
   new%nqr = 10
   write(std_out,"(a)")"Activating computation of Frohlich self-energy:"
   write(std_out,"(2(a,i0,1x))")"ntheta:", new%ntheta, "nphi:", new%nphi
   write(std_out,"((a,i0,1x,a,f5.3,1x,a))")"nqr points:", new%nqr, "qrad:", new%qrad, "[Bohr^-1]"

   ! Initialize angular mesh qvers_cart and angwgth (inspired to initang in m_pawang)
   ABI_MALLOC(th, (new%ntheta))
   ABI_MALLOC(wth, (new%ntheta))
   con = two_pi / new%nphi
   call gauleg(-one, one, th, wth, new%ntheta)

   ! Initialize qvers_cart and angwgth
   new%angl_size = new%ntheta * new%nphi
   ABI_MALLOC(new%qvers_cart, (3, new%angl_size))
   ABI_MALLOC(new%angwgth, (new%angl_size))
   npoints = 0
   do it = 1, new%ntheta
     cos_theta = th(it)
     sin_theta = sqrt(one - cos_theta*cos_theta)
     do ip = 1, new%nphi
       ang = con * (ip-1)
       cos_phi = cos(ang)
       sin_phi = sin(ang)
       npoints = npoints + 1
       new%qvers_cart(1, npoints) = sin_theta * cos_phi
       new%qvers_cart(2, npoints) = sin_theta * sin_phi
       new%qvers_cart(3, npoints) = cos_theta
       ! Normalization required
       new%angwgth(npoints) = wth(it) / (two * new%nphi)
     end do
   end do

   ABI_FREE(th)
   ABI_FREE(wth)

   !if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine == 0) .or.&
   !    (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /= 0) )  then
   !  wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)//'FINE'
   !  if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
   !    MSG_ERROR(msg)
   !  end if
   !  call wrtout(std_out,"Will read energies from: "//trim(wfk_fname_dense))
   !  ebands_dense = wfk_read_ebands(path, comm) result(ebands)
   !  ebands_ptr => ebands_dense
   !  call ebands_free(ebands_dense)
   !else
   !  ebands_ptr => ebands
   !end if

   ! Build SKW object for all bands where QP corrections are wanted.
   ! TODO: check band_block because I got weird results (don't remember if with AbiPy or Abinit)
   skw_cplex = 1; if (kpts_timrev_from_kptopt(ebands%kptopt) == 0) skw_cplex = 2
   skw_band_block = [minval(new%bstart_ks), maxval(new%bstart_ks + new%nbcalc_ks - 1)]
   params = 0; params(1) = 1; params(2) = 5
   if (nint(dtset%einterp(1)) == 1) params = dtset%einterp
   write(std_out, "(a, 4(f5.2, 2x))")"SKW parameters used to interpolate e_{nk+q} in Frohlich self-energy:", params
   new%frohl_skw = skw_new(cryst, params(2:), skw_cplex, ebands%mband, ebands%nkpt, ebands%nsppol, &
     ebands%kptns, ebands%eig, skw_band_block, comm)

   ! Allocate arrays to store Frohlich results for given (ikcalc, spin)
   ABI_CALLOC(new%frohl_vals_e0ks, (new%ntemp, new%max_nbcalc))
   ABI_CALLOC(new%frohl_dvals_de0ks, (new%ntemp, new%max_nbcalc))
   if (new%nwr > 0) then
     ABI_CALLOC(new%frohl_vals_wr, (new%nwr, new%ntemp, new%max_nbcalc))
   end if
 end if

 ! Open netcdf file (only master works for the time being because I cannot assume HDF5 + MPI-IO)
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
     "eph_task", "symsigma", "nbsum", "symdynmat", "ph_intmeth", "eph_intmeth", "qint_method", "eph_transport", &
     "imag_only", "frohl_model"])
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
     nctkarr_t("qpoms_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ze0_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ks_gaps", "dp", "nkcalc, nsppol"), &
     nctkarr_t("qpoms_gaps", "dp", "ntemp, nkcalc, nsppol"), &
     nctkarr_t("qp_gaps", "dp", "ntemp, nkcalc, nsppol") &
   ])
   NCF_CHECK(ncerr)

   ! DEBUG: Write gkq matrix elements to file. DO NOT USE IN PRODUCTION
   !ncerr = nctk_def_arrays(ncid, [ &
   !  nctkarr_t("gkq", "dp", "two, max_nbcalc, nbsum, natom3, nqbz, nkcalc, nsppol"), &
   !  nctkarr_t("gkq_frohl", "dp", "two, max_nbcalc, nbsum, natom3, nqbz, nkcalc, nsppol") &
   !])
   !NCF_CHECK(ncerr)

   if (new%frohl_model == 1) then
     ! Arrays storing the Frohlich self-energy.
     ! These arrays get two extra dimensions on file (nkcalc, nsppol).
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("frohl_vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("frohl_dvals_de0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol") ])
     NCF_CHECK(ncerr)
     if (new%nwr > 0) then
       ncerr = nctk_def_arrays(ncid, [ &
         nctkarr_t("frohl_vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
         nctkarr_t("frohl_spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol") &
       ])
       NCF_CHECK(ncerr)
     end if
   end if

   if (new%nwr > 0) then
     ! Make room for the spectral function.
     ! These arrays get two extra dimensions on file (nkcalc, nsppol).
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("wrmesh_b", "dp", "nwr, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol") &
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

   ! ======================================================
   ! Write data that do not depend on the (kpt, spin) loop.
   ! ======================================================
   NCF_CHECK(nctk_set_datamode(ncid))
   ii = 0; if (new%imag_only) ii = 1
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "symsigma", "nbsum", "symdynmat", "ph_intmeth", "eph_intmeth", "qint_method", &
     "eph_transport", "imag_only", "frohl_model"], &
     [dtset%eph_task, new%symsigma, new%nbsum, dtset%symdynmat, dtset%ph_intmeth, dtset%eph_intmeth, new%qint_method, &
     dtset%eph_transport, ii, new%frohl_model])
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

 ! Now reopen the file (not xmpi_comm_self --> only master writes)
 call xmpi_barrier(comm)
 NCF_CHECK(nctk_open_modify(new%ncid, strcat(dtfil%filnam_ds(4), "_SIGEPH.nc"), xmpi_comm_self))
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
 !if (allocated(self%bstart_ks)) then
 !  ABI_FREE(self%bstart_ks)
 !end if
 ABI_SFREE(self%bstart_ks)
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
 !if (allocated(self%frohl_gkq2)) then
 !  ABI_FREE(self%frohl_gkq2)
 !end if
 if (allocated(self%qvers_cart)) then
   ABI_FREE(self%qvers_cart)
 end if
 if (allocated(self%angwgth)) then
   ABI_FREE(self%angwgth)
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
 !if (allocated(self%indkk)) then
 !  ABI_FREE(self%indkk)
 !end if
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
 if (allocated(self%frohl_vals_e0ks)) then
   ABI_FREE(self%frohl_vals_e0ks)
 end if
 if (allocated(self%dvals_de0ks)) then
   ABI_FREE(self%dvals_de0ks)
 end if
 if (allocated(self%frohl_dvals_de0ks)) then
   ABI_FREE(self%frohl_dvals_de0ks)
 end if
 if (allocated(self%dw_vals)) then
   ABI_FREE(self%dw_vals)
 end if
 if (allocated(self%vals_wr)) then
   ABI_FREE(self%vals_wr)
 end if
 if (allocated(self%frohl_vals_wr)) then
   ABI_FREE(self%frohl_vals_wr)
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

 call skw_free(self%frohl_skw)

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
 !integer,parameter :: sppoldbl1 = 1
 !real(dp) :: dksqmax
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
 !if (allocated(self%indkk)) then
 !  ABI_FREE(self%indkk)
 !end if
 !ABI_MALLOC(self%indkk, (self%nqibz_k, 6))
 !call listkk(dksqmax, cryst%gmet, self%indkk, self%qibz, self%qibz_k, self%nqibz, self%nqibz_k, cryst%nsym,&
 !   sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)
 !if (dksqmax > tol12) MSG_ERROR("Wrong mapping")

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
!!  prtvol= Verbosity level
!!  comm=MPI communicator.
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_gather_and_write(self, ebands, ikcalc, spin, prtvol, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_gather_and_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, spin, prtvol, comm
 type(sigmaph_t),target,intent(inout) :: self
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 integer,parameter :: master=0, max_ntemp=1
 integer :: ideg,ib,it,ii,iw,nstates,ierr,my_rank,band_ks,ik_ibz,ibc,ib_val,ib_cond,jj
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_oms,qpe_oms_val,qpe_oms_cond
 real(dp) :: cpu, wall, gflops, invsig2fmts, tau
 complex(dpc) :: sig0c,sig0fr,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2
 character(len=500) :: msg
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
!arrays
 integer :: shape3(3),shape4(4),shape5(5),shape6(6)
 integer, ABI_CONTIGUOUS pointer :: bids(:)
 !real(dp), ABI_CONTIGUOUS pointer :: rdata3(:,:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:),rdata6(:,:,:,:,:,:)
 real(dp) :: qp_gaps(self%ntemp),qpoms_gaps(self%ntemp)
 real(dp),allocatable :: aw(:,:,:), frohl_aw(:,:,:)
 real(dp) :: ks_enes(self%max_nbcalc),ze0_vals(self%ntemp, self%max_nbcalc)
 real(dp) :: gfw_avg(self%gfw_nomega, 3)
 complex(dpc) :: qpoms_enes(self%ntemp, self%max_nbcalc),qp_enes(self%ntemp, self%max_nbcalc)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call cwtime(cpu, wall, gflops, "start")
 call xmpi_sum_master(self%vals_e0ks, master, comm, ierr)
 call xmpi_sum_master(self%dvals_de0ks, master, comm, ierr)
 call xmpi_sum_master(self%dw_vals, master, comm, ierr)
 if (self%nwr > 0) call xmpi_sum_master(self%vals_wr, master, comm, ierr)
 call cwtime(cpu, wall, gflops, "stop", comm=comm)
 call wrtout(std_out, sjoin("Sigma_{nk} gather completed. Average cpu-time:", sec2str(cpu), &
     ", Total Average wall-time:", sec2str(wall), ch10), do_flush=.True.)

 ! Only master writes
 if (my_rank /= master) return

 ik_ibz = self%kcalc2ibz(ikcalc, 1)

 ! Add Frohlich results to final results (frohl arrays are not symmetrized if symsigma == 1)
 if (self%frohl_model == 1) then
   self%vals_e0ks = self%vals_e0ks + self%frohl_vals_e0ks
   self%dvals_de0ks = self%dvals_de0ks + self%frohl_dvals_de0ks
   if (self%nwr > 0) self%vals_wr = self%vals_wr + self%frohl_vals_wr
 end if

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
 qp_enes = huge(one) * tol6; qpoms_enes = huge(one) * tol6
 ks_enes = huge(one) * tol6; ze0_vals = huge(one) * tol6
 ks_gap = -one; qpoms_gaps = -one; qp_gaps = -one

 ! Write legend.
 if (ikcalc == 1 .and. spin == 1) then
   write(ab_out,"(a)")repeat("=", 80)
   write(ab_out,"(a)")"Final results in eV."
   write(ab_out,"(a)")"Notations:"
   write(ab_out,"(a)")"   eKS: Kohn-Sham energy. eQP: quasi-particle energy."
   write(ab_out,"(a)")"   eQP - eKS: Difference between the QP and the KS energy."
   write(ab_out,"(a)")"   SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part."
   if (self%frohl_model == 1) then
     write(ab_out,"(a)")"   SF1(eKS): Real part of the (model) Frohlich self-energy at the KS energy, SF2 for imaginary part."
   end if
   write(ab_out,"(a)")"   Z(eKS): Renormalization factor."
   write(ab_out,"(a)")"   FAN: Real part of the Fan term at eKS. DW: Debye-Waller term."
   write(ab_out,"(a)")"   DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP."
   write(ab_out,"(a)")"   OTMS: On-the-mass-shell approximation with eQP ~= eKS + Sigma(omega=eKS)"
   write(ab_out,"(a)")"   TAU(eKS): Lifetime in femtoseconds computed at the KS energy."
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
 end if

 do it=1,min(self%ntemp, max_ntemp)
   ! Write header.
   if (self%nsppol == 1) then
     write(ab_out,"(3a,f4.1,a)") &
       "K-point: ", trim(ktoa(self%kcalc(:,ikcalc))), ", T= ", self%kTmesh(it) / kb_HaK, " [K]"
   else
     write(ab_out,"(3a,i1,a,f4.1,a)") &
       "K-point: ", trim(ktoa(self%kcalc(:,ikcalc))), ", spin: ", spin, ", T= ",self%kTmesh(it) / kb_HaK, " [K]"
   end if
   if (self%imag_only) then
     if (self%frohl_model == 0) then
       write(ab_out,"(a)")"   B    eKS    SE2(eKS)  TAU(eKS)  DeKS"
     else
       write(ab_out,"(a)")"   B    eKS    SE2(eKS)  SF2(eKS)  TAU(eKS)  DeKS"
     end if
   else
     if (self%frohl_model == 0) then
       write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
     else
       write(ab_out,"(a)")&
         "   B    eKS     eQP    eQP-eKS   SE1(eKS)  SF1(eKS)  SE2(eKS)  SF2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
     end if
   end if

   do ibc=1,self%nbcalc_ks(ikcalc,spin)
     band_ks = self%bstart_ks(ikcalc,spin) + ibc - 1
     kse = ebands%eig(band_ks, ik_ibz, spin)
     ks_enes(ibc) = kse
     sig0c = self%vals_e0ks(it, ibc)
     dw = self%dw_vals(it, ibc)
     fan0 = real(sig0c) - dw
     if (self%frohl_model == 1) sig0fr = self%frohl_vals_e0ks(it, ibc)
     ! Compute QP energies with On-the-Mass-Shell approximation and first renormalization i.e. Z(eKS)
     ! TODO: Note that here I use the full Sigma including the imaginary part
     !zc = one / (one - self%dvals_de0ks(it, ibc))
     zc = one / (one - real(self%dvals_de0ks(it, ibc)))
     ze0_vals(it, ibc) = real(zc)
     qpe = kse + real(zc) * real(sig0c)
     qpe_oms = kse + real(sig0c)
     if (ibc == 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     if (band_ks == ib_val) then
       kse_val = kse; qpe_val = qpe; qpe_oms_val = qpe_oms
     end if
     if (band_ks == ib_cond) then
       kse_cond = kse; qpe_cond = qpe; qpe_oms_cond = qpe_oms
     end if
     ! FIXME
     if (self%imag_only) then
       invsig2fmts = Time_Sec * 1e+15 / two
       tau = 999999.0_dp
       if (abs(aimag(sig0c)) > tol16) tau = invsig2fmts / abs(aimag(sig0c))
       tau = min(tau, 999999.0_dp)
       if (self%frohl_model == 0) then
         ! B    eKS     SE2(eKS)  TAU(eKS)  DeKS
         write(ab_out, "(i4,2(f8.3,1x),f8.1,1x,f8.3)") &
           band_ks, kse * Ha_eV, aimag(sig0c) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
       else
         ! B    eKS     SE2(eKS)  SF2(eKS)  TAU(eKS)  DeKS
         write(ab_out, "(i4,2(f8.3,1x),2(f8.1,1x),f8.3)") &
           band_ks, kse * Ha_eV, aimag(sig0c) * Ha_eV, aimag(sig0fr) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
       end if
     else
       if (self%frohl_model == 0) then
         ! B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP
         write(ab_out, "(i4, 10(f8.3,1x))") &
           band_ks, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
           real(sig0c) * Ha_eV, aimag(sig0c) * Ha_eV, real(zc), &
           fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
       else
         ! B    eKS     eQP    eQP-eKS   SE1(eKS)  SF1(eKS) SE2(eKS) SF2(eKS) Z(eKS)  FAN(eKS)   DW      DeKS     DeQP
         write(ab_out, "(i4, 12(f8.3,1x))") &
           band_ks, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
           real(sig0c) * Ha_eV, real(sig0fr) * Ha_eV, aimag(sig0c) * Ha_eV, aimag(sig0fr) * Ha_eV, real(zc), &
           fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
       end if
     end if
     if (ibc > 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     qpoms_enes(it, ibc) = qpe_oms
     qp_enes(it, ibc) = qpe
     if (kse_val /= huge(one) .and. kse_cond /= huge(one)) then
       ! We have enough states to compute the gap.
       if (it == 1) ks_gap = kse_cond - kse_val
       qpoms_gaps(it) = qpe_oms_cond - qpe_oms_val
       qp_gaps(it) = real(qpe_cond - qpe_val)
     end if
   end do ! ibc

   ! Print KS and QP gaps.
   if (.not. self%imag_only) then
     if (kse_val /= huge(one) .and. kse_cond /= huge(one)) then
       write(ab_out, "(a)")" "
       write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, &
         "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
       write(ab_out, "(2(a,f8.3),a)")" QP gap: ",qp_gaps(it) * Ha_eV," (OTMS: ",qpoms_gaps(it) * Ha_eV, ")"
       write(ab_out, "(2(a,f8.3),a)")" QP_gap - KS_gap: ",(qp_gaps(it) - ks_gap) * Ha_eV,&
           " (OTMS: ",(qpoms_gaps(it) - ks_gap) * Ha_eV, ")"
       write(ab_out, "(a)")" "
     end if
   else
     if (kse_val /= huge(one) .and. kse_cond /= huge(one)) then
       write(ab_out, "(a)")" "
       write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
       write(ab_out, "(a)")" "
     end if
   end if

   write(ab_out, "(a)")repeat("=", 92)
 end do ! it

 if (self%ntemp > max_ntemp .and. (ikcalc == 1 .and. spin == 1)) then
   write(ab_out, "(a,i0,a)")"No more than ", max_ntemp, " temperatures are written to the main output file."
   write(ab_out, "(2a)")"Please use SIGEPH.nc file and AbiPy to analyze the results.",ch10
 end if

 if (prtvol > 0 .and. (ikcalc == 1 .and. spin == 1)) then
   if (self%gfw_nomega > 0) then
     write(ab_out, "(2a)")"omega and Eliashberg function gf_{nk}(omega) for testing purposes:"
     iw = (self%gfw_nomega / 2)
     do ib=1,min(self%nbcalc_ks(ikcalc, spin), 5)
       band_ks = self%bstart_ks(ikcalc, spin) + ib - 1
       write(ab_out, "(a, i0)")"For band:", band_ks
       do jj=0,1
         write(ab_out, "(4(f8.3,2x))")self%gfw_mesh(iw+jj), (self%gfw_vals(iw+jj, ii, ib), ii=1,3)
       end do
     end do
     write(ab_out, "(a)")ch10
   end if

   if (self%nwr >= 3) then
     write(ab_out, "(2a)")ch10,"omega and Sigma_nk(omega, T=1) in eV for testing purposes:"
     it = 1; iw = (self%nwr / 2)
     do ib=1,min(self%nbcalc_ks(ikcalc, spin), 5)
       band_ks = self%bstart_ks(ikcalc, spin) + ib - 1
       write(ab_out, "(a, i0)")"For band:", band_ks
       do ii=0,1
         write(ab_out, "(3(f8.3,2x))")self%wrmesh_b(iw+ii, ib) * Ha_eV, self%vals_wr(iw+ii, it, ib) * Ha_eV
       end do
     end do
     write(ab_out, "(a)")ch10
   end if
 end if

 call cwtime(cpu, wall, gflops, "start")
#ifdef HAVE_NETCDF

 ! Write self-energy matrix elements for this (kpt, spin)
 ! NB: Only master writes
 ! (use iso_c_binding to associate a real pointer to complex data because netcdf does not support complex types).
 ! Well, cannot use c_loc with gcc <= 4.8 due to internal compiler error so use c2r and stack memory.
 !shape3(1) = 2; shape4(1) = 2; shape5(1) = 2; shape6(1) = 2

 !shape3(2:) = shape(self%vals_e0ks); call c_f_pointer(c_loc(self%vals_e0ks), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "vals_e0ks"), c2r(self%vals_e0ks), start=[1,1,1,ikcalc,spin]))

 !shape3(2:) = shape(self%dvals_de0ks); call c_f_pointer(c_loc(self%dvals_de0ks), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "dvals_de0ks"), c2r(self%dvals_de0ks), start=[1,1,1,ikcalc,spin]))

 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "dw_vals"), self%dw_vals, start=[1,1,ikcalc,spin]))

 ! Dump QP energies and gaps for this (kpt, spin)
 !shape3(2:) = shape(qpoms_enes); call c_f_pointer(c_loc(qpoms_enes), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qpoms_enes"), c2r(qpoms_enes), start=[1,1,1,ikcalc,spin]))

 !shape3(2:) = shape(qp_enes); call c_f_pointer(c_loc(qp_enes), rdata3, shape3)
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qp_enes"), c2r(qp_enes), start=[1,1,1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "ze0_vals"), ze0_vals, start=[1,1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "ks_enes"), ks_enes, start=[1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "ks_gaps"), ks_gap, start=[ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qpoms_gaps"), qpoms_gaps, start=[1,ikcalc,spin]))
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qp_gaps"), qp_gaps, start=[1,ikcalc,spin]))

 if (self%frohl_model == 1) then
   ! Write Forhlich self-energy to file.
   ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_vals_e0ks"), &
      c2r(self%frohl_vals_e0ks), start=[1,1,1,ikcalc,spin])
   NCF_CHECK(ncerr)
   ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_dvals_de0ks"), &
      c2r(self%frohl_dvals_de0ks), start=[1,1,1,ikcalc,spin])
   NCF_CHECK(ncerr)
   if (self%nwr > 0) then
     ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_vals_wr"), &
       c2r(self%frohl_vals_wr), start=[1,1,1,1,ikcalc,spin])
     NCF_CHECK(ncerr)
   end if
 end if

 ! Write frequency dependent data.
 if (self%nwr > 0) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "wrmesh_b"), self%wrmesh_b, start=[1,1,ikcalc,spin]))

   !shape4(2:) = shape(self%vals_wr); call c_f_pointer(c_loc(self%vals_wr), rdata4, shape4)
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "vals_wr"), c2r(self%vals_wr), start=[1,1,1,1,ikcalc,spin]))

   ! Compute spectral function.
   ! A = -1/pi [Im Sigma(ww)] / ([ww - ee - Re Sigma(ww)] ** 2 + Im Sigma(ww) ** 2])
   ABI_MALLOC(aw, (self%nwr, self%ntemp, self%max_nbcalc))
   ABI_MALLOC(frohl_aw, (self%nwr, self%ntemp, self%max_nbcalc))
   do ib=1,self%nbcalc_ks(ikcalc,spin)
     band_ks = self%bstart_ks(ikcalc, spin) + ib - 1
     kse = ebands%eig(band_ks, ik_ibz, spin)
     do it=1,self%ntemp
       aw(:, it, ib) = -piinv * aimag(self%vals_wr(:, it, ib)) / &
         ((self%wrmesh_b(:, ib) - kse - real(self%vals_wr(:, it, ib))) ** 2 + aimag(self%vals_wr(:, it, ib)) ** 2)
       if (self%frohl_model /= 0) then
         ! Spectral function associated to Frohlich model.
         frohl_aw(:, it, ib) = -piinv * aimag(self%frohl_vals_wr(:, it, ib)) / &
           ((self%wrmesh_b(:, ib) - kse - real(self%frohl_vals_wr(:, it, ib))) ** 2 + aimag(self%frohl_vals_wr(:, it, ib)) ** 2)
       end if
     end do
   end do
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "spfunc_wr"), aw, start=[1, 1, 1, ikcalc, spin]))
   if (self%frohl_model /= 0) then
     NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_spfunc_wr"), frohl_aw, start=[1, 1, 1, ikcalc, spin]))
   end if
   ABI_FREE(aw)
   ABI_FREE(frohl_aw)
 end if

 ! Write Eliashberg functions
 if (self%gfw_nomega > 0) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "gfw_vals"), self%gfw_vals, start=[1, 1, 1, ikcalc, spin]))
 end if
#endif

 call cwtime(cpu, wall, gflops, "stop")
 call wrtout(std_out, sjoin("Sigma_{nk} netcdf output completed. cpu-time:", sec2str(cpu), &
     ", Total wall-time:", sec2str(wall), ch10), do_flush=.True.)

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
 !if (self%imag_only) then write(unt,"(a)")sjoin("From bmin:", itoa(), "to bmax:", itoa())
 write(unt,"(a)")sjoin("Symsigma: ",itoa(self%symsigma), "Timrev:", itoa(self%timrev))
 write(unt,"(a)")sjoin("Imaginary shift in the denominator (zcut): ", ftoa(aimag(self%ieta) * Ha_eV, fmt="f5.3"), "[eV]")
 msg = "Standard quadrature"; if (self%qint_method == 1) msg = "tetrahedron method"
 write(unt, "(2a)")sjoin("Method for q-space integration:", msg)
 if (self%imag_only) write(unt, "(a)")"Only the Imaginary part of Sigma will be computed."
 if (.not. self%imag_only) write(unt, "(a)")"Both Real and Imaginary part of Sigma will be computed."
 write(unt,"(a)")sjoin("Number of frequencies along the real axis:", itoa(self%nwr), &
    ", Step:", ftoa(self%wr_step*Ha_eV, fmt="f5.3"), "[eV]")
 write(unt, "(a)")sjoin("Number of frequency in generalized Eliashberg functions:", itoa(self%gfw_nomega))
 write(unt,"(a)")sjoin("Number of temperatures:", itoa(self%ntemp), &
   "From:", ftoa(self%kTmesh(1) / kb_HaK), "to", ftoa(self%kTmesh(self%ntemp) / kb_HaK), "[K]")
 write(unt,"(a)")sjoin("Ab-initio q-mesh from DDB file:", ltoa(dtset%ddb_ngqpt))
 write(unt,"(a)")sjoin("Q-mesh used for self-energy integration [ngqpt]:", ltoa(self%ngqpt))
 write(unt,"(a)")sjoin("Number of q-points in the IBZ:", itoa(self%nqibz))
 write(unt,"(a)")sjoin("asr:", itoa(dtset%asr), "dipdip:", itoa(dtset%dipdip), "symdynmat:", itoa(dtset%symdynmat))
 if (self%frohl_model == 0) then
   write(unt,"(a)")"No special treatment of Frohlich divergence in gkq for q --> 0"
 else
   write(unt,"(a)")"Activating computation of Frohlich self-energy to treat divergence in gkq for q --> 0"
   write(unt,"(2(a,i0,1x))")"ntheta:", self%ntheta, "nphi:", self%nphi
   write(unt,"((a,i0,1x,a,f5.3,1x,a))")"nr points:", self%nqr, "qrad:", self%qrad, "[Bohr^-1]"
 end if
 write(unt,"(a)")"List of K-points for self-energy corrections:"
 do ikc=1,self%nkcalc
   if (ikc > 20) then
     write(unt, "(a)")"nkcalc > 20. Stop printing more k-point information."
     exit
   end if
   do is=1,self%nsppol
     if (self%nsppol == 2) write(unt,"(a,i1)")"... For spin: ",is
     write(unt, "(2(i4,2x),a,2(i4,1x))") &
       ikc, is, trim(ktoa(self%kcalc(:,ikc))), self%bstart_ks(ikc,is), self%bstart_ks(ikc,is) + self%nbcalc_ks(ikc,is) - 1
     end do
 end do

end subroutine sigmaph_print
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_get_all_qweights
!! NAME
!!  sigmaph_get_all_qweights
!!
!! FUNCTION
!!  Compute all the weights for q-space integration using the tetrahedron method
!!
!! INPUTS
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  spin: Spin index
!!  ikcalc: Index of the self-energy k-point in the kcalc array.
!!  comm: MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_get_all_qweights(sigma,cryst,ebands,spin,ikcalc,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_get_all_qweights'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigmaph_t),intent(inout) :: sigma
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: ikcalc, spin, comm

!Local variables ------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: nu, band_ks, ibsum_kq, ik_ibz, ib_k, bstart_ks, nbcalc_ks, my_rank, natom3, ierr
 integer :: nprocs,this_calc
 real(dp) :: eig0nk, eminmax(2)
 real(dp) :: cpu,wall,gflops
 real(dp),allocatable :: tmp_deltaw_pm(:,:,:)
 character(len=500) :: msg

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 ik_ibz = sigma%kcalc2ibz(ikcalc,1)
 nbcalc_ks = sigma%nbcalc_ks(ikcalc,spin)
 bstart_ks = sigma%bstart_ks(ikcalc,spin)
 natom3 = 3 * cryst%natom

 sigma%deltaw_pm = zero

 call cwtime(cpu,wall,gflops,"start")
 ABI_MALLOC(tmp_deltaw_pm,(3,sigma%ephwg%nq_k, 2))
 ! loop over bands to sum
 do ibsum_kq=1,sigma%nbsum
  ! loop over phonon modes
  do nu=1,natom3
    ! loop over bands in the self-energy
    do ib_k=1,nbcalc_ks
      this_calc = (ibsum_kq-1)*natom3*nbcalc_ks + (nu-1)*nbcalc_ks + ib_k
      if (mod(this_calc,nprocs) /= my_rank) cycle
      band_ks = ib_k + bstart_ks - 1
      eig0nk = ebands%eig(band_ks, ik_ibz, spin)
      eminmax(1) = eig0nk - 0.01
      eminmax(2) = eig0nk + 0.01
      call ephwg_get_deltas(sigma%ephwg, ibsum_kq, spin, nu, 3, eminmax, bcorr0, tmp_deltaw_pm, xmpi_comm_self)
      ! we pay the efficiency here
      sigma%deltaw_pm(1,ib_k,nu,ibsum_kq,:) = tmp_deltaw_pm(2, :, 1) / ( sigma%ephwg%lgk%weights(:) )
      sigma%deltaw_pm(2,ib_k,nu,ibsum_kq,:) = tmp_deltaw_pm(2, :, 2) / ( sigma%ephwg%lgk%weights(:) )
    enddo
  enddo
 enddo

 ! TODO: Reintegrate cweights
 ! Compute \int 1/z with tetrahedron if both real and imag part of sigma are wanted.
 !ABI_MALLOC(zvals, (nz, nbc))
 !zvals(1, 1:nbc) = sigma%e0vals + sigma%ieta
 !call ephwg_zinv_weights(sigma%ephwg, iqlk, nz, nbc, zvals, ibsum_kq, spin, sigma%cweights, xmpi_comm_self, &
 !  use_bzsum=sigma%symsigma == 0)
 !ABI_FREE(zvals)

 ABI_FREE(tmp_deltaw_pm)
 call xmpi_sum(sigma%deltaw_pm, comm, ierr)

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2))') "weights with tetrahedron  cpu:",cpu,", wall:",wall
 call wrtout(std_out, msg, do_flush=.True.)

end subroutine sigmaph_get_all_qweights
!!***

!!****f* m_sigmaph/eval_sigfrohl
!! NAME
!!  eval_sigfrohl
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine eval_sigfrohl(sigma, cryst, ifc, ebands, ikcalc, spin, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eval_sigfrohl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(sigmaph_t),intent(inout) :: sigma
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: ikcalc, spin, comm

!Local variables ------------------------------
!scalars
 integer :: my_rank, nprocs, nbcalc_ks, iq_ibz, i1, i2, i3, nu, it, natom3, iatom, bstart_ks
 integer :: ib_k, band_ks, ik_ibz, iang, iqr, ierr, iw
 real(dp) :: wqnu, nqnu, eig0nk, eig0mkq, f_mkq ,gkq2, qrad2, weigth_q, vol_fact
 real(dp) :: inv_qepsq, qstep, gmod2, hmod2, rfact
 complex(dpc) :: cfact, cnum
!arrays
 integer :: ndivs(3)
 real(dp) :: q0(3), qpt(3), qpt_cart(3), kq(3), kk(3) !, qvers_cart(3)
 real(dp) :: qrmesh(sigma%nqr), gmod2r(sigma%nqr), hmod2r(sigma%nqr), rfactr(sigma%nqr), nqr(sigma%nqr)
 real(dp) :: phfrq(cryst%natom*3)
 !real(dp) :: oder1(3),oder2(3,3)
 real(dp),allocatable :: displ_cart(:,:,:,:), eigs_kq(:), eigs_kqr(:,:), wqr(:,:)
 real(dp),allocatable ::  f_mkqr(:), eig0mkqr(:), gkqr2(:,:)
 complex(dpc) :: cdd(3), cfqr(sigma%nqr)
 complex(dpc),allocatable :: cfact_wr(:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 natom3 = 3 * cryst%natom
 nbcalc_ks = sigma%nbcalc_ks(ikcalc, spin)
 bstart_ks = sigma%bstart_ks(ikcalc, spin)
 ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
 kk = sigma%kcalc(:, ikcalc)

 ndivs = 1
 qrad2 = sigma%qrad ** 2

 ! Allocate workspace arrays.
 ABI_MALLOC(displ_cart, (2,3,cryst%natom,3*cryst%natom))
 ABI_MALLOC(eigs_kq, (nbcalc_ks))
 if (sigma%nwr > 0) then
   ABI_MALLOC(cfact_wr, (sigma%nwr))
 end if

 ! Integrations over q-points in the IBZ(k).
 ! Use densified mesh in each q-microzone.
 ! Use spherical cutoff of radius qrad around q = 0 when computing gkq^{LR}

 sigma%frohl_vals_e0ks = zero
 sigma%frohl_dvals_de0ks = zero
 if (sigma%nwr > 0) sigma%frohl_vals_wr = zero

 do iq_ibz=1,sigma%nqibz_k
    if (mod(iq_ibz, nprocs) /= my_rank) cycle ! MPI parallelism
    !write(std_out,*)"iq_ibz:", iq_ibz
    q0 = sigma%qibz_k(:,iq_ibz)
    ! q-weigths for IBZ(k) integration
    weigth_q = sigma%wtq_k(iq_ibz) ! * (four_pi / cryst%ucvol)

    ! Build densified grid inside microzone centered around iq_ibz
    do i1=-ndivs(1),ndivs(1),1
      qpt(1) = q0(1) + real(i1) / (ndivs(1) * sigma%ngqpt(1))
      do i2=-ndivs(2),ndivs(2),1
        qpt(2) = q0(2) + real(i2) / (ndivs(2) * sigma%ngqpt(2))
        do i3=-ndivs(3),ndivs(3),1
          qpt(3) = q0(3) + real(i2) / (ndivs(3) * sigma%ngqpt(3))

          qpt_cart = two * pi * matmul(cryst%gprimd, qpt)
          ! TODO: In principle one should introduce weigths for points that are close to the border
          ! Should write routine to compute weights given kptrlatt, kibz and radius.
          if (dot_product(qpt_cart, qpt_cart) <= qrad2) cycle
          vol_fact = one
          call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart)
          inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))

          ! Interpolate e_{n k+q}
          kq = kk + qpt
          do ib_k=1,nbcalc_ks
            band_ks = ib_k + bstart_ks - 1
            call skw_eval_bks(sigma%frohl_skw, band_ks, kq, spin,  eigs_kq(ib_k))
          end do

          ! Sum over modes for this q-point.
          do nu=1,natom3
            wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
            ! Compute gkq_{LR} without (i 4pi/ucvol)
            ! note that at this level it does not dependend on ib_k.
            cnum = zero
            do iatom=1,cryst%natom
              cdd = cmplx(displ_cart(1,:, iatom, nu), displ_cart(2,:, iatom, nu)) * exp(-j_dpc * two_pi * cryst%xred(:, iatom))
              !cdd = cdd * exp(-q/self%qdamp**2)
              cnum = cnum + dot_product(qpt_cart, matmul(ifc%zeff(:, :, iatom), cdd))
            end do
            gkq2 = (real(cnum) ** 2 + aimag(cnum) ** 2) / (two * wqnu) * inv_qepsq * weigth_q * vol_fact

            ! For each band in Sigma_{bk}
            do ib_k=1,nbcalc_ks
              band_ks = ib_k + bstart_ks - 1
              eig0nk = ebands%eig(band_ks, ik_ibz, spin)
              eig0mkq = eigs_kq(ib_k)

              do it=1,sigma%ntemp
                nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)
                f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))
                cfact =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                         (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)

                if (sigma%imag_only) then
                  sigma%frohl_vals_e0ks(it, ib_k) = sigma%frohl_vals_e0ks(it, ib_k) + gkq2 * j_dpc * aimag(cfact)
                else
                  sigma%frohl_vals_e0ks(it, ib_k) = sigma%frohl_vals_e0ks(it, ib_k) + gkq2 * cfact

                  ! Accumulate d(Re Sigma) / dw(w=eKS) for state ib_k
                  !cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
                  !            (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
                  gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
                  hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
                  rfact = (nqnu + f_mkq      ) * (-gmod2 + aimag(sigma%ieta)**2) / (gmod2 + aimag(sigma%ieta)**2) ** 2 + &
                          (nqnu - f_mkq + one) * (-hmod2 + aimag(sigma%ieta)**2) / (hmod2 + aimag(sigma%ieta)**2) ** 2
                  sigma%frohl_dvals_de0ks(it, ib_k) = sigma%frohl_dvals_de0ks(it, ib_k) + gkq2 * rfact

                  ! Accumulate Sigma(w) for state ib_k if spectral function is wanted.
                  ! TODO: weigths
                  if (sigma%nwr > 0) then
                    cfact_wr(:) = (nqnu + f_mkq      ) / (sigma%wrmesh_b(:,ib_k) - eig0mkq + wqnu + sigma%ieta) + &
                                  (nqnu - f_mkq + one) / (sigma%wrmesh_b(:,ib_k) - eig0mkq - wqnu + sigma%ieta)
                    sigma%frohl_vals_wr(:, it, ib_k) = sigma%frohl_vals_wr(:, it, ib_k) + gkq2 * cfact_wr(:)
                  end if

                end if
              end do ! it
            end do ! ib_k
          end do ! nu

        end do ! i3
      end do ! i2
    end do ! i1
 end do ! iq_ibz

 !write(std_out,*)"IBZ(k) done"

 ! Now integrate inside sphere of radius qcut using spherical coordinates.
 ABI_MALLOC(eigs_kqr, (sigma%nqr, nbcalc_ks))
 ABI_MALLOC(wqr, (sigma%nqr, natom3))
 ABI_MALLOC(f_mkqr, (sigma%nqr))
 ABI_MALLOC(eig0mkqr, (sigma%nqr))
 ABI_CALLOC(gkqr2, (sigma%nqr, natom3))

 qstep = sigma%qrad / (sigma%nqr - 1)
 qrmesh = arth(zero, qstep, sigma%nqr)

 ! Weighted summation over angles
 do iang=1, sigma%angl_size
    if (mod(iang, nprocs) /= my_rank) cycle
    !write(std_out, *)"iang=", iang
    !qvers_cart = matmul(cryst%rprimd, sigma%qvers_cart(:, iang))
    ! Precompute denominator
    inv_qepsq = one / dot_product(sigma%qvers_cart(:, iang), matmul(ifc%dielt, sigma%qvers_cart(:, iang)))

    ! Get all omega_{q nu} and e_{n k+q} along the q-line for this (theta, phi)
    do iqr=1,sigma%nqr
      ! TODO: Decide how to treat the first point
      qpt_cart = sigma%qvers_cart(:, iang) * qrmesh(iqr)
      !qpt = matmul(..., qpt_cart)
      qpt = zero
      if (iqr == 1) then
        ! Include non-analytical part even for q==0.
        call ifc_fourq(ifc, cryst, qpt_cart, phfrq, displ_cart, nanaqdir="cart")
      else
        call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart)
      end if
      wqr(iqr, :) = phfrq

      ! Compute gkq_{LR}. Note that in our approx it does not dependend on ib_k.
      do nu=1,natom3
        wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
        cnum = zero
        do iatom=1,cryst%natom
          ! This is complex
          cdd = cmplx(displ_cart(1,:, iatom, nu), displ_cart(2,:, iatom, nu))
          !cdd = cdd * exp(-q**2)
          cnum = cnum + dot_product(qpt_cart, matmul(ifc%zeff(:, :, iatom), cdd))
        end do
        gkqr2(iqr, nu) = (real(cnum) ** 2 + aimag(cnum) ** 2) / (two * wqnu) * inv_qepsq
      end do

      ! Interpolate e_{n k+q} and store results.
      kq = kk + qpt
      do ib_k=1,nbcalc_ks
        band_ks = ib_k + bstart_ks - 1
        call skw_eval_bks(sigma%frohl_skw, band_ks, kq, spin, eigs_kqr(iqr, ib_k))
      end do
    end do ! iqr

    ! Sum over phonon modes.
    do nu=1,natom3
      wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle

      ! For each band in Sigma_{bk}
      do ib_k=1,nbcalc_ks
        band_ks = ib_k + bstart_ks - 1
        eig0nk = ebands%eig(band_ks, ik_ibz, spin)
        eig0mkqr(:) = eigs_kqr(:, ib_k)

        ! For each temperature and corresponding Fermi level and phonon occ.
        do it=1,sigma%ntemp
          nqr = occ_be(wqr(:, nu), sigma%kTmesh(it), zero)
          f_mkqr = occ_fd(eig0mkqr, sigma%kTmesh(it), sigma%mu_e(it))

          ! Radial integration for Sigma(w=e_{nk})
          do iqr=1,sigma%nqr
            cfqr(iqr) = gkqr2(iqr, nu) * ( &
               (nqr(iqr) + f_mkqr(iqr)      ) / (eig0nk - eig0mkqr(iqr) + wqr(iqr, nu) + sigma%ieta) + &
               (nqr(iqr) - f_mkqr(iqr) + one) / (eig0nk - eig0mkqr(iqr) - wqr(iqr, nu) + sigma%ieta) )
          end do
          cfact = simpson_cplx(sigma%nqr, qstep, cfqr) * sigma%angwgth(iang)

          if (sigma%imag_only) then
            ! NB: Here we are not completely consistent if the integration is done with tetrahedra.
            ! In principle one should avoid ieta and use some kind of tessellation for the sphere.
            sigma%frohl_vals_e0ks(it, ib_k) = sigma%frohl_vals_e0ks(it, ib_k) + j_dpc * aimag(cfact)
          else
            sigma%frohl_vals_e0ks(it, ib_k) = sigma%frohl_vals_e0ks(it, ib_k) + cfact

            ! Accumulate d(Re Sigma) / dw(w=eKS) for state ib_k
            ! cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
            !             (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
            gmod2r(:) = (eig0nk - eig0mkqr(:) + wqr(:, nu)) ** 2
            hmod2r(:) = (eig0nk - eig0mkqr(:) - wqr(:, nu)) ** 2
            rfactr(:) = gkqr2(:, nu) * &
              (nqr(:) + f_mkqr(:)      ) * (-gmod2r(:) + aimag(sigma%ieta)**2) / (gmod2r(:) + aimag(sigma%ieta)**2) ** 2 + &
              (nqr(:) - f_mkqr(:) + one) * (-hmod2r(:) + aimag(sigma%ieta)**2) / (hmod2r(:) + aimag(sigma%ieta)**2) ** 2

            sigma%frohl_dvals_de0ks(it, ib_k) = sigma%frohl_dvals_de0ks(it, ib_k) + simpson(qstep, rfactr) * sigma%angwgth(iang)

            ! Accumulate Sigma(w) for state ib_k if spectral function is wanted.
            ! This is gonna be costly ...
            ! TODO: weigths
            if (sigma%nwr > 0) then
              do iw=1,sigma%nwr
                cfqr(:) = gkqr2(:, nu) * ( &
                  (nqr(:) + f_mkqr(:)      ) / (sigma%wrmesh_b(iw, ib_k) - eig0mkqr(:) + wqr(:, nu) + sigma%ieta) + &
                  (nqr(:) - f_mkqr(:) + one) / (sigma%wrmesh_b(iw, ib_k) - eig0mkqr(:) - wqr(:, nu) + sigma%ieta) )
                cfact = simpson_cplx(sigma%nqr, qstep, cfqr)
                sigma%frohl_vals_wr(iw, it, ib_k) = sigma%frohl_vals_wr(iw, it, ib_k) * cfact * sigma%angwgth(iang)
              end do
            end if
          end if

        end do
      end do
    end do

 end do ! idir

 ABI_FREE(displ_cart)
 ABI_FREE(eigs_kq)
 ABI_FREE(eigs_kqr)
 ABI_FREE(wqr)
 ABI_FREE(eig0mkqr)
 ABI_FREE(f_mkqr)
 ABI_FREE(gkqr2)

 if (sigma%nwr > 0) then
   ABI_FREE(cfact_wr)
 end if

 ! Accumulate final results
 call xmpi_sum(sigma%frohl_vals_e0ks, comm, ierr)
 call xmpi_sum(sigma%frohl_dvals_de0ks, comm, ierr)
 if (sigma%nwr > 0) call xmpi_sum(sigma%frohl_vals_wr, comm, ierr)

end subroutine eval_sigfrohl
!!***

end module m_sigmaph
!!***

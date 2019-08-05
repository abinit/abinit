!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigmaph
!! NAME
!!  m_sigmaph
!!
!! FUNCTION
!!  Compute the matrix elements of the Fan-Migdal Debye-Waller self-energy in the KS basis set.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG, HM)
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
 use m_copy
 use m_ifc
 use m_ebands
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfd
 use m_skw
 use m_krank
 use m_lgroup
 use m_ephwg
 use m_sort
 use m_hdr
 use m_sigtk
 use m_eph_double_grid
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_rf2

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime, cwtime_report, timab
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, simpson
 use m_io_tools,       only : iomode_from_fname, file_exists
 use m_special_funcs,  only : gaussian
 use m_fftcore,        only : ngfft_seq
 use m_cgtk,           only : cgtk_rotate
 use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm
 use m_crystal,        only : crystal_t
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, listkk
 use m_occ,            only : occ_fd, occ_dfd, occ_be
 use m_fftcore,        only : get_kg
 use m_kg,             only : getph
 use m_bz_mesh,        only : isamek
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_ioarr,          only : read_rhor
 use m_pawang,         only : pawang_type, gauleg
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
 use m_dfpt_cgwf,      only : dfpt_cgwf

 implicit none

 private
!!***

 ! Tables for degenerated KS states.
 type bids_t
   integer, allocatable :: vals(:)
 end type bids_t

 type degtab_t
   type(bids_t), allocatable :: bids(:)
 end type degtab_t

 ! Store the weights in single or double precision
 integer,private,parameter :: DELTAW_KIND = dp
 !integer,private,parameter :: DELTAW_KIND = sp


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

 type,public :: sigmaph_t

  integer :: nkcalc
   ! Number of k-points computed.

  integer :: max_nbcalc
   ! Maximum number of bands computed (max over nkcalc and spin).

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

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
   ! timrev=1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nbsum
   ! Total number of bands used in sum over states without taking into account MPI distribution.

  integer :: bsum_start, bsum_stop
   ! Fist and last band included in self-energy sum without taking into account MPI distribution.
   ! nbsum = bsum_stop - bsum_start + 1

  integer :: my_bsum_start, my_bsum_stop
   ! Initial and final band index included in self-energy sum
   ! Processor-dependent if Re-Im calculation.
   ! Processor-independent and computed at runtime on the basis of the nk states in Sigma_{nk} if imag_only

  integer :: my_npert
   ! Number of atomic perturbations or phonon modes treated by this MPI rank.

  integer :: comm_pert = xmpi_comm_null
   ! MPI communicator for parallelism over atomic perturbations.

  integer :: nprocs_pert
   ! Number of cpus for parallelism over atomic perturbations.

  integer :: me_pert
   ! My rank in comm_pert

  integer :: comm_bq = xmpi_comm_null
   ! MPI communicator used to distribute bsum, q-points (high-level parallelization)

  integer :: me_bq
  integer :: nprocs_bq

  integer :: comm_qpt = xmpi_comm_null
  ! MPI communicator for q-points
  integer :: me_qpt
  integer :: nprocs_qpt

  integer :: comm_bsum = xmpi_comm_null
  ! MPI communicator for bands in sum
  integer :: me_bsum
  integer :: nprocs_bsum

  integer :: coords(3)

  integer :: nqbz
  ! Number of q-points in the (dense) BZ for sigma integration

  integer :: nqibz
  ! Number of q-points in the (dense) IBZ for sigma integration

  integer :: nqibz_k
  ! Number of q-points in the IBZ(k). Depends on ikcalc.

  integer :: my_nqibz_k
  ! Number of q-points in the IBZ(k) treated by this MPI proc. Depends on ikcalc.
  ! Differs from nqibz_k only if imag with tetra because in this case we can introduce a cutoff on the weights

  integer :: ncid = nctk_noid
  ! Netcdf file used to save results.

  integer :: mpw
  ! Maximum number of PWs for all possible k+q

  integer :: bcorr = 0
  ! 1 to include Blochl correction in tetrahedron method else 0.

  integer :: ntheta, nphi
  ! Number of division for spherical integration of Frohlich term.

  integer :: nqr = 0
   ! Number of points on the radial mesh for spherical integration of the Frohlich self-energy

  integer :: angl_size = 0
   ! Dimension of angular mesh for spherical integration of the Frohlich self-energy
   ! angl_size = ntheta * nphi

  complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  real(dp) :: elow, ehigh
   ! min and Max KS energy treated in self-energy +- max phonon energy
   ! Used to select bands in self-energy sum if imag only and select q-points in qpoints_oracle

  real(dp) :: winfact = four
   ! winfact * wmax is used to define the energy window for filtering electronic states
   ! in the computation of electron lifetimes.

  real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  real(dp) :: qrad = zero
   ! Radius of the sphere for the numerical integration of the Frohlich self-energy

  real(dp) :: qdamp = one
   ! Exponential damping e-(qmod**2/qdamp**2) added to the Frohlich matrix elements in q-space.

  real(dp) :: wmax
  ! Max phonon energy + buffer. Used to select the bands to sum for the imaginary part
  ! and filter q-points on the basis of electron energy difference.

  integer :: qint_method
   ! Defines the method used to integrate in q-space
   ! 0 -> Standard quadrature (one point per micro zone).
   ! 1 -> Use tetrahedron method.

  integer :: frohl_model
   ! > 0 to treat the q-->0 divergence and accelerate convergence in polar semiconductors.
   !   1: Use spherical integration inside the micro zone around the Gamma point
   !   2: More advanced (and expensive) treatmentent inside the BZ (NOT IMPLEMENTED YET)

  logical :: use_doublegrid = .False.
   ! whether to use double grid or not

  logical :: use_ftinterp = .False.
   ! whether DFPT potentials should be read from the DVDB or Fourier-interpolated on the fly.

  type(eph_double_grid_t) :: eph_doublegrid
   ! store the double grid related object

  logical :: imag_only
   ! True if only the imaginary part of the self-energy must be computed

  logical :: calc_mrta
   ! True if the self-energy in the energy-momentum relaxation time approximation should be computed

  integer :: gmax(3)

  integer :: ngqpt(3)
   ! Number of divisions in the Q mesh in the BZ.

  integer,allocatable :: bstart_ks(:,:)
   ! bstart_ks(nkcalc,nsppol)
   ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used.

  integer,allocatable :: bstop_ks(:,:)
   ! bstop_ks(nkcalc,nsppol)

  integer,allocatable :: nbcalc_ks(:,:)
   ! nbcalc_ks(nkcalc,nsppol)
   ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used.

  integer,allocatable :: kcalc2ibz(:,:)
    !kcalc2ibz(nkcalc, 6))
    ! Mapping kcalc --> ibz as reported by listkk.

  integer,allocatable :: myq2ibz_k(:)
   ! myq2ibz_k(my_nqibz_k)
   ! Mapping my q-point index --> index in nqibz_k arrays (IBZ_k)
   ! Differs from nqibz_k only if imag with tetra because in this case we can introduce a cutoff.

  integer(i1b),allocatable :: itreat_qibz(:)
   ! itreat_qibz(nqibz)
   ! Table used to distribute potentials over q-points in the IBZ.
   ! The loop over qpts in the IBZ(k) is MPI distributed inside comm_qpt accordinging to this table.
   ! 0 if this IBZ point is not treated by the procs inside comm_pert_bsum
   ! 1 or 2-9 if this IBZ is treated.

  integer,allocatable :: my_pinfo(:,:)
    ! my_pinfo(3, my_npert)
    ! my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
    ! my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
    ! my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3

  integer,allocatable :: pert_table(:,:)
   ! pert_table(2, natom3)
   ! pert_table(1, npert): rank of the processor treating this atomic perturbation.
   ! pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.

  integer,allocatable :: phmodes_skip(:)
   ! (natom3) A mask to skip accumulating the contribution of certain phonon modes

  integer,allocatable:: ind_bz2ibz(:,:)
   ! (6, %nqibz)
   ! Mapping q IBZ to BZ

  integer,allocatable:: indkk_kq(:, :)
   ! (6, %nqibz_k))
   ! Mapping k+q --> initial IBZ. Depends on ikcalc.

  !integer,allocatable :: ind_ibz2dvdb(:,:)

  integer,allocatable :: indq2dvdb_k(:,:)
   ! (6, %nqibz_k))
   ! Mapping qibz_k --> IBZ found in DVDB file. Used when DFPT potentials are read from DVDB file
   ! so that we know how to access/symmetrize v1scf
   ! Depends on ikcalc.

  integer,allocatable :: ind_ibzk2ibz(:,:)
   ! (6, %nqibz_k))
   ! Mapping qibz_k --> IBZ defined by eph_ngqpt_fine.
   ! Depends on ikcalc.

  integer,allocatable :: ibz2dvdb(:)
   ! (6, %nqpts))
   ! Mapping dvdb%ibz --> %ibz

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

  real(dp),allocatable :: vcar_calc(:,:,:,:)
  ! (3, max_nbcalc, nkcalc, nsppol))
  ! Diagonal elements of velocity operator in cartesian coordinates for all states in Sigma_nk.

  real(dp),allocatable :: linewidth_mrta(:,:)
   ! linewidth_mrta(ntemp, max_nbcalc)
   ! Linewidths computed withing the momentum relaxation time approximation
   ! for fixed (kcalc, spin). Only if imag_only

  complex(dpc),allocatable :: cweights(:,:,:,:,:,:,:)
  ! (nz, 2, nbcalc_ks, my_npert, my_bsum_start:my_bsum_stop, my_nqibz_k, ndiv))
  ! Weights for the q-integration of 1 / (e1 - e2 \pm w_{q, nu} + i.eta)
  ! This array is initialized inside the (ikcalc, spin) loop

  real(kind=DELTAW_KIND),allocatable :: deltaw_pm(:,:,:,:,:,:)
  ! (2, nbcalc_ks, my_npert, bsum_start:bsum_start, my_nqibz_k, ndiv))
  ! Weights for the q-integration of the two delta (abs/emission) if imag_only
  ! This array is initialized inside the (ikcalc, spin) loop

  real(dp),allocatable :: wrmesh_b(:,:)
  ! wrmesh_b(nwr, max_nbcalc)
  ! Frequency mesh along the real axis (Ha units) used for the different bands
  ! Each mesh is **centered** on the corresponding KS energy.
  ! This array depends on (ikcalc, spin)

  real(dp), allocatable :: qvers_cart(:,:)
   ! qvers_cart(3, angl_size)
   ! For each point of the angular mesh, gives the Cartesian coordinates
   ! of the corresponding point on an unitary sphere (Frohlich self-energy)

  real(dp), allocatable :: angwgth(:)
   ! angwgth(angl_size)
   ! For each point of the angular mesh, gives the weight
   ! of the corresponding point on an unitary sphere (Frohlich self-energy)

  real(dp),allocatable :: frohl_deltas_sphcorr(:, :, :, :)
   ! (2, ntemp, max_nbcalc, natom3))
   ! Integration of the imaginary part inside the small sphere around Gamma
   ! computed numerically with the Frohlich model by Verdi and angular integration.
   ! The first dimension stores the contributions due to +/- omega_qn
   ! Used if frohl_model == 1 and imag_only. This array depend on (ikcalc, spin)

  integer, allocatable :: qp_done(:,:)
   ! qp_done(kcalc,spin)
   ! Keep track of the QP states already computed for restart of the calculation

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

  integer :: gfw_nomega = 0
   ! Number of phonon frequencies in Eliashberg function.
   ! Set to 0 to deactivate this part.
   ! Integration in q-space is done according to eph_intmeth.

  real(dp),allocatable :: gfw_mesh(:)
   ! gfw_mesh(gfw_nomega)
   ! Frequency mesh for Eliashberg function (Allen-Cardona adiabatic, phonons only)

   real(dp),allocatable :: gf_nnuq(:,:,:,:)
   ! (nbcalc_ks, natom3, %nqibz_k, 3)
   ! Quantities needed to compute generalized Eliashberg functions (gkq2/Fan-Migdal/DW terms)
   ! This array depends on (ikcalc, spin)
   ! NB: q-weights for integration are not included.

  real(dp),allocatable :: gfw_vals(:,:,:)
   !gfw_vals(gfw_nomega, 3, max_nbcalc)
   ! Generalized Eliashberg function a2F_{n,k,spin}(w)
   ! 1:   gkk^2 with delta(en - em)
   ! 2:3 (Fan-Migdal/DW contribution in adiabatic approximation)
   ! This array depends on (ikcalc, spin)

  integer :: a2f_ne = 0
   ! Number of points in a2f_emesh

  real(dp),allocatable :: a2f_emesh(:)
  ! a2f_emesh(a2f_ne)

  real(dp),allocatable :: a2few(:,:,:)
   ! a2few(a2f_ne, gwf_nomega, max_nbcalc)

  type(ephwg_t) :: ephwg
   ! This object compute the weights for the BZ integration in q-space if qint_method > 0

  type(skw_t) :: frohl_skw
   ! Star-function interpolator.
   ! Used to interpolate \epsilon_{nk+q) in the numerical integration of the Frohlich self-energy.

  type(degtab_t),allocatable :: degtab(:,:)
   ! (nkcalc, nsppol)
   ! Table used to average QP results in the degenerate subspace if symsigma == 1

  contains

    procedure :: write => sigmaph_write
     ! Write main dimensions and header of sigmaph on a netcdf file.

    procedure :: comp  => sigmaph_comp
     ! Compare two instances of sigmaph raise error if different

    procedure :: setup_kcalc => sigmaph_setup_kcalc
     ! Return tables used to perform the sum over q-points for given k-point.

    procedure :: gather_and_write => sigmaph_gather_and_write
     ! Compute the QP corrections.

    procedure :: print => sigmaph_print
     ! Print results to main output file.

    procedure :: free => sigmaph_free
      ! Free sigmaph object

    procedure :: get_ebands => sigmaph_get_ebands
      ! Fill in values in ebands from the sigmaph structure and netcdf file

 end type sigmaph_t
!!***

 public :: sigmaph         ! Main entry point to compute self-energy matrix elements
 public :: sigmaph_read    ! Read main dimensions and header of sigmaph from a netcdf file.
 private :: sigmaph_new   ! Creation method (allocates memory, initialize data from input vars).

 real(dp),private,parameter :: EPH_WTOL = tol6
 ! Tolerance for phonon frequencies to be ignored.

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha

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
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! wfk_hdr=Header of the WFK file.
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

subroutine sigmaph(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, wfk_hdr, &
                   pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

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
 type(hdr_type),intent(in) :: wfk_hdr
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_getgh1c=1,berryopt0=0,timrev0=0
 integer,parameter :: useylmgr=0,useylmgr1=0,master=0,ndat1=1,sppoldbl1=1,timrev1=1
 integer,parameter :: igscq0=0, icgq0 = 0, usedcwavef0 = 0, nbdbuf0 = 0, quit0 = 0, cplex1 = 1, pawread0 = 0
 integer :: my_rank,nsppol,nkpt,iq_ibz,iq_ibz_frohl,iq_bz_frohl, my_npert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,nprocs
 integer :: ibsum_kq,ib_k,band_ks,ibsum,ii,jj, iw
 integer :: mcgq, mgscq, nband_kq
 integer :: idir,ipert,ip1,ip2,idir1,ipert1,idir2,ipert2
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq
 integer :: iq_ibz_fine,ikq_ibz_fine,ikq_bz_fine
 integer :: spin,istwf_k,istwf_kq,istwf_kqirr,npw_k,npw_kq,npw_kqirr
 integer :: mpw,ierr,it,imyq,band, ignore_kq, ignore_ibsum_kq
 integer :: n1,n2,n3,n4,n5,n6,nspden,nu, iang
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,nq,cnt,imyp, q_start, q_stop, restart
 integer :: nlines_done, nline_in, grad_berry_size_mpw1
 integer :: nbcalc_ks,nbsum,bsum_start, bsum_stop, bstart_ks,ikcalc,bstart,bstop,iatom
 integer :: comm_rpt, method
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,gflops_all,cpu_ks,wall_ks,gflops_ks,cpu_dw,wall_dw,gflops_dw
 real(dp) :: cpu_setk, wall_setk, gflops_setk, cpu_qloop, wall_qloop, gflops_qloop, gf_val
 real(dp) :: ecut,eshift,weight_q,rfact,gmod2,hmod2,ediff,weight, inv_qepsq, qmod, fqdamp, simag, q0rad, out_resid
 real(dp) :: vkk_norm2
 complex(dpc) :: cfact,dka,dkap,dkpa,dkpap,cplx_ediff, cnum, sig_cplx
 logical :: isirr_k,isirr_kq,gen_eigenpb,isqzero
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(sigmaph_t) :: sigma, sigma_restart
 type(ddkop_t) :: ddkop
 type(rf2_t) :: rf2
 type(hdr_type) :: pot_hdr
 character(len=500) :: msg
 character(len=fnlen) :: sigeph_path
!arrays
 integer :: g0_k(3),g0_kq(3)
 integer :: work_ngfft(18),gmax(3) !g0ibz_kq(3),
 integer(i1b),allocatable :: itreatq_dvdb(:)
 integer,allocatable :: gtmp(:,:),kg_k(:,:),kg_kq(:,:),nband(:,:), qselect(:), wfd_istwfk(:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),qpt_cart(3),phfrq(3*cryst%natom), dotri(2),qq_ibz(3)
 real(dp) :: vk(2, 3), vk_red(2,3), vkq(2,3), vkq_red(2,3), tsec(2), eminmax(2)
 real(dp) :: frohl_sphcorr(3*cryst%natom), vec_natom3(2, 3*cryst%natom)
 real(dp) :: wqnu,nqnu,gkq2,gkq2_pf,gkq2_dfrohl,eig0nk,eig0mk,eig0mkq,f_mkq
 real(dp),allocatable :: displ_cart(:,:,:,:),displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:),gdw2_mn(:,:),gkq0_atm(:,:,:,:),gkq2_lr(:,:)
 real(dp),allocatable :: cgq(:,:,:), gscq(:,:,:), out_eig1_k(:), cg1s_kq(:,:,:,:), h1kets_kq_allperts(:,:,:,:)
 real(dp),allocatable :: dcwavef(:, :), gh1c_n(:, :), ghc(:,:), gsc(:,:), stern_ppb(:,:,:,:), stern_dw(:,:,:,:)
 logical,allocatable :: ihave_ikibz_spin(:,:)
 complex(dpc),allocatable :: tpp_red(:,:)
 complex(dpc) :: cdd(3), cp3(3)
 real(dp),allocatable :: bra_kq(:,:),kets_k(:,:,:),h1kets_kq(:,:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: vtrial(:,:),gvnlx1(:,:),gvnlxc(:,:),work(:,:,:,:)
 real(dp),allocatable :: gs1c(:,:),nqnu_tlist(:),dtw_weights(:,:),dt_tetra_weights(:,:,:),dwargs(:),alpha_mrta(:)
 real(dp),allocatable :: delta_e_minus_emkq(:)
 complex(dpc),allocatable :: cfact_wr(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:), cwaveprj(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)

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
 ! Check if a previous netcdf file is present and restart the calculation
 ! Here we try to read an existing SIGEPH file if eph_restart
 ! we compare the variables with the state of the code (i.e. new sigmaph generated in sigmaph_new)
 restart = 0; ierr = 1
 sigeph_path = strcat(dtfil%filnam_ds(4), "_SIGEPH.nc")
 if (my_rank == master .and. dtset%eph_restart == 1) then
   sigma_restart = sigmaph_read(sigeph_path, dtset, xmpi_comm_self, msg, ierr)
 end if

 sigma = sigmaph_new(dtset, ecut, cryst, ebands, ifc, dtfil, comm)

 if (my_rank == master .and. dtset%eph_restart == 1) then
   if (ierr == 0) then
     if (any(sigma_restart%qp_done /= 1)) then
       call sigma%comp(sigma_restart)
       sigma%qp_done = sigma_restart%qp_done
       restart = 1
       call wrtout([std_out, ab_out], "- Restarting from previous SIGEPH.nc file")
     else
       restart = 0; sigma%qp_done = 0
       MSG_COMMENT("Found SIGEPH.nc file with all QP entries already computed. Will overwrite file.")
     end if
   end if
   call sigma_restart%free()
 end if

 call xmpi_bcast(restart, master, comm, ierr)
 call xmpi_bcast(sigma%qp_done, master, comm, ierr)

 if (restart == 0) then
   call sigma%write(dtset, cryst, ebands, wfk_hdr, dtfil, comm)
 else
   if (my_rank == master) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_modify(sigma%ncid, sigeph_path, xmpi_comm_self))
#endif
   end if
 end if

 if (my_rank == master) then
   call sigma%print(dtset, ab_out)
   call sigma%print(dtset, std_out)
 end if
 my_npert = sigma%my_npert

 ! This is the maximum number of PWs for all possible k+q treated.
 mpw = sigma%mpw; gmax = sigma%gmax

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))

 ! Initialize the wave function descriptor.
 ! Each node has all k-points and spins and bands between my_bsum_start and my_bsum_stop
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, nkpt ,nsppol))

 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.

 ! Each node needs the wavefunctions for Sigma_{nk}
 do spin=1,sigma%nsppol
   do ikcalc=1,sigma%nkcalc
     ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
     bstart = sigma%bstart_ks(ikcalc,spin)
     bstop = bstart + sigma%nbcalc_ks(ikcalc,spin) - 1
     bks_mask(bstart:bstop, ik_ibz, spin) = .True.
   end do
 end do

 if (sigma%imag_only .and. sigma%qint_method == 1) then
   call wrtout(std_out, " Including restricted set of states within energy window around bdgw states.")
   do spin=1,sigma%nsppol
     do ik_ibz=1,ebands%nkpt
       do band=sigma%my_bsum_start, sigma%my_bsum_stop
         eig0mk = ebands%eig(band, ik_ibz, spin)
         if (eig0mk >= sigma%elow  - sigma%winfact * sigma%wmax .and. &
             eig0mk <= sigma%ehigh + sigma%winfact * sigma%wmax) then
            bks_mask(band, ik_ibz ,spin) = .True.
         end if
       end do
     end do
   end do
   ! Uncomment these lines to disable energy window trick and allocate all bands.
   if (dtset%userie == 123) then
     call wrtout(std_out, " Storing all bands between my_bsum_start and my_bsum_stop.")
     bks_mask(sigma%my_bsum_start:sigma%my_bsum_stop, : ,:) = .True.
   end if
 else
   bks_mask(sigma%my_bsum_start:sigma%my_bsum_stop, : ,:) = .True.
 endif

 if (dtset%userie == 124) then
   ! Uncomment this line to have all states on each MPI rank.
   bks_mask = .True.
   call wrtout(std_out, " Storing all bands for debugging purposes.")
 end if

 ! This table is needed when computing tau:
 ! k+q states outside energy window are not read hence their contribution won't be included.
 ! Error is small provided calculation is close to convergence.
 ! TODO: Store min/max energy difference ...
 ABI_MALLOC(ihave_ikibz_spin, (nkpt, nsppol))
 ihave_ikibz_spin = .False.
 do spin=1,sigma%nsppol
   do ik_ibz=1,ebands%nkpt
     if (any(bks_mask(:, ik_ibz, spin))) ihave_ikibz_spin(ik_ibz, spin) = .True.
   end do
 end do

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,dtset%mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)

 call wfd%print(header="Wavefunctions for self-energy calculation.",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))
 if (.False.) call wfd%test_ortho(cryst, pawtab, unit=std_out, mode_paral="PERS")

 ! if PAW, one has to solve a generalized eigenproblem
 ! BE careful here because I will need sij_opt == -1
 usecprj = 0
 gen_eigenpb = psps%usepaw == 1
 sij_opt = 0; if (gen_eigenpb) sij_opt = 1

 ABI_DT_MALLOC(cwaveprj0, (natom, nspinor*usecprj))
 ABI_DT_MALLOC(cwaveprj, (natom, nspinor*usecprj))
 ABI_MALLOC(displ_cart, (2,3,cryst%natom,3*cryst%natom))
 ABI_MALLOC(displ_red, (2,3,cryst%natom,3*cryst%natom))
 ABI_MALLOC(tpp_red, (natom3, natom3))

 !if (sigma%calc_velocity == 1) then
 call cwtime(cpu_ks, wall_ks, gflops_ks, "start", &
     msg=" Computing diagonal elements of velocity operator for all states in Sigma_nk...")
 ABI_MALLOC(cgwork,   (2, mpw*wfd%nspinor))
 ABI_CALLOC(sigma%vcar_calc, (3, sigma%max_nbcalc, sigma%nkcalc, nsppol))

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
 !if (my_rank == master) call ddkop%print(ab_out)
 ! All sigma_nk states are available on each node so parallelization is easy.
 cnt = 0
 do spin=1,nsppol
   do ikcalc=1,sigma%nkcalc
     kk = sigma%kcalc(:, ikcalc)
     bstart_ks = sigma%bstart_ks(ikcalc, spin)
     ik_ibz = sigma%kcalc2ibz(ikcalc,1)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

     do ib_k=1,sigma%nbcalc_ks(ikcalc, spin)
       cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
       band_ks = ib_k + bstart_ks - 1
       call wfd%copy_cg(band_ks, ik_ibz, spin, cgwork)
       eig0nk = ebands%eig(band_ks, ik_ibz, spin)
       call ddkop%apply(eig0nk, npw_k, wfd%nspinor, cgwork, cwaveprj0, wfd%mpi_enreg)
       vk_red = ddkop%get_velocity(eig0nk, istwf_k, npw_k, wfd%nspinor, wfd%mpi_enreg%me_g0, cgwork)
       call ddk_red2car(cryst%rprimd,vk_red,vk)
       sigma%vcar_calc(:, ib_k, ikcalc, spin) = vk(1, :)
     end do

   end do
 end do

 call xmpi_sum(sigma%vcar_calc, comm, ierr)
 call cwtime_report(" Velocities", cpu_ks, wall_ks, gflops_ks)

 ! Write results to disk.
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK(nctk_set_defmode(sigma%ncid))
   NCF_CHECK(nctk_def_arrays(sigma%ncid, [nctkarr_t("vcar_calc", "dp", "three, max_nbcalc, nkcalc, nsppol")]))
   NCF_CHECK(nctk_set_datamode(sigma%ncid))
   NCF_CHECK(nf90_put_var(sigma%ncid, nctk_idname(sigma%ncid, "vcar_calc"), sigma%vcar_calc))
 end if
#endif

 ABI_FREE(cgwork)

 if (.not.sigma%calc_mrta) call ddkop%free()

 ! Radius of sphere with volume equivalent to the micro zone.
 q0rad = two_pi * (three / (four_pi * cryst%ucvol * sigma%nqbz)) ** third

 frohl_sphcorr = zero
 if (sigma%frohl_model == 1 .and. .not. sigma%imag_only) then
   ! Prepare treatment of Frohlich divergence with spherical integration in the microzone around Gamma.
   ! Correction does not depend on (nk) so we precompute values at this level.
   call wrtout(std_out, " Computing spherical average to treat Frohlich divergence ...")
   do iang=1,sigma%angl_size
     if (mod(iang, nprocs) /= my_rank) cycle ! MPI parallelism
     qpt_cart = sigma%qvers_cart(:, iang)
     inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))

     call ifc%fourq(cryst, qpt_cart, phfrq, displ_cart, nanaqdir="cart")

     ! Note that acoustic modes are ignored.
     do nu=4,natom3
       wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
       cp3 = czero
       do iatom=1, natom
         cp3 = cp3 + matmul(ifc%zeff(:, :, iatom), cmplx(displ_cart(1,:,iatom, nu), displ_cart(2,:,iatom, nu), kind=dpc))
       end do
       cnum = dot_product(qpt_cart, cp3)
       ! Compute spherical average.
       frohl_sphcorr(nu) = frohl_sphcorr(nu) + sigma%angwgth(iang) * abs(cnum) ** 2 * inv_qepsq ** 2 / wqnu ** 2
     end do
   end do
   call xmpi_sum(frohl_sphcorr, comm, ierr)

   frohl_sphcorr = frohl_sphcorr * eight * pi / cryst%ucvol * (three / (four_pi * cryst%ucvol * sigma%nqbz)) ** third
   if (my_rank == master) then
     write(std_out, "(/,a)")" Frohlich model integrated inside the small q-sphere around Gamma: "
     write(std_out, "(a)")" This correction is used to accelerate the convergence of the ZPR with the q-point sampling "
     write(std_out, "(a)")" Note that this term tends to zero for N_q --> oo "
     write(std_out, "(a,/)")" so it's different from the integral of the Frohlich potential in the full BZ."
     do nu=1,natom3
       write(std_out, "(a,i0,a,f8.3,a)")" For phonon mode: ", nu," sphcorr:", frohl_sphcorr(nu) * Ha_eV, " [eV]"
     end do
     write(std_out, "(a)")ch10
   end if
 end if

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1  ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2     ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0 ! gvnlx1 is output

 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nsppol,nspden,natom,&
  dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
  comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
  usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 ! Allocate work space arrays.
 ! vtrial and vlocal are required for Sternheimer (H0). DFPT routines do not need it.
 ! Note nvloc in vlocal (we will select one/four spin components afterwards)
 ABI_CALLOC(vtrial, (nfftf, nspden))
 ABI_CALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))

 if (dtset%eph_stern == 1) then
   ! Read GS POT (vtrial) from input POT file
   ! In principle one may store vtrial in the DVDB but getpot_path is simpler to implement.
   call wrtout([std_out, ab_out], sjoin(" Reading GS KS potential for Sternheimer from: ", dtfil%filpotin))
   call read_rhor(dtfil%filpotin, cplex1, nspden, nfftf, ngfftf, pawread0, mpi_enreg, vtrial, pot_hdr, pawrhoij, comm, &
                  allow_interp=.True.)
   call hdr_free(pot_hdr)
 end if

 if (sigma%nwr > 0) then
   ABI_MALLOC(cfact_wr, (sigma%nwr))
 end if
 ABI_MALLOC(nqnu_tlist, (sigma%ntemp))

 ! Allocate workspace arrays for Eliashberg calculation.
 if (dtset%prteliash /= 0) then
   ABI_MALLOC(dtw_weights, (sigma%gfw_nomega, 2))
   ABI_MALLOC(dwargs, (sigma%gfw_nomega))
   if (sigma%a2f_ne > 0) then
     ABI_MALLOC(delta_e_minus_emkq, (sigma%a2f_ne))
   end if
 end if

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 ! This to symmetrize the DFPT potentials.
 dvdb%symv1 = dtset%symv1scf
 if (sigma%nprocs_pert > 1) then
   call dvdb%set_pert_distrib(sigma%my_npert, natom3, sigma%my_pinfo, sigma%pert_table, sigma%comm_pert)
 end if

 sigma%use_ftinterp = .False.

 ! Find correspondence IBZ --> set of q-points in DVDB.
 ABI_MALLOC(sigma%ibz2dvdb, (sigma%nqibz))
 if (dvdb%find_qpts(sigma%nqibz, sigma%qibz, sigma%ibz2dvdb, comm) /= 0) then
   call wrtout([std_out, ab_out], " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
   sigma%use_ftinterp = .True.
 else
   call wrtout([std_out, ab_out], " DVDB file contains all q-points in the IBZ --> Reading DFPT potentials from file.")
   sigma%use_ftinterp = .False.
 end if

 if (sigma%use_ftinterp) then
   ! Use ddb_ngqpt q-mesh to compute real-space represention of DFPT v1scf potentials to prepare Fourier interpolation.
   ! R-points are distributed inside comm_rpt
   ! Note that when R-points are distributed inside comm_qpt we cannot interpolate potentials on-the-fly
   ! inside the loop over q-points.
   ! In this case, indeed, the interpolation must be done in sigma_setup_qloop once we know the q-points contributing
   ! to the integral and the potentials must be cached.
   comm_rpt = sigma%comm_qpt
   if (.not. sigma%imag_only) comm_rpt = sigma%comm_bsum
   !FIXME: comm_qpt is buggy.
   !if (sigma%imag_only) comm_rpt = xmpi_comm_self
   !comm_rpt = sigma%comm_bsum
   comm_rpt = xmpi_comm_self
   method = dtset%userid
   call dvdb%ftinterp_setup(dtset%dvdb_ngqpt, dtset%ddb_qrefine, 1, dtset%ddb_shiftq, &
                            nfftf, ngfftf, method, comm_rpt)

   ! Build q-cache in the *dense* IBZ using the global mask qselect and itreat_qibz.
   ABI_MALLOC(qselect, (sigma%nqibz))
   qselect = 1
   if (sigma%imag_only .and. sigma%qint_method == 1) then
     call qpoints_oracle(sigma, cryst, ebands, sigma%qibz, sigma%nqibz, sigma%nqbz, sigma%qbz, qselect, comm)
     !qselect = 1
   end if
   call dvdb%ftqcache_build(nfftf, ngfftf, sigma%nqibz, sigma%qibz, dtset%dvdb_qcache_mb, qselect, sigma%itreat_qibz, comm)

 else
   ABI_MALLOC(qselect, (dvdb%nqpt))
   qselect = 1
   ! Try to predict the q-points required to compute tau.
   if (sigma%imag_only .and. sigma%qint_method == 1) then
     call qpoints_oracle(sigma, cryst, ebands, dvdb%qpts, dvdb%nqpt, sigma%nqbz, sigma%qbz, qselect, comm)
   end if
 end if

 call dvdb%print(prtvol=dtset%prtvol)

 if (.not. sigma%use_ftinterp) then
   ! Need to translate itreat_qibz into itreatq_dvdb.
   ABI_ICALLOC(itreatq_dvdb, (dvdb%nqpt))
   do iq_ibz=1,sigma%nqibz
     if (sigma%itreat_qibz(iq_ibz) == 0) cycle
     db_iqpt = sigma%ibz2dvdb(iq_ibz)
     ABI_CHECK(db_iqpt /= -1, sjoin("Could not find IBZ q-point:", ktoa(sigma%qibz(:, iq_ibz)), "in DVDB file."))
     itreatq_dvdb(db_iqpt) = 1
   end do
   !itreatq_dvdb = 1
   call dvdb%qcache_read(nfftf, ngfftf, dtset%dvdb_qcache_mb, qselect, itreatq_dvdb, comm)
   ABI_FREE(itreatq_dvdb)
 end if

 ABI_FREE(qselect)

 ! Loop over k-points in Sigma_nk. Loop over spin is internal we operate on nspden components at once.
 do ikcalc=1,sigma%nkcalc
   ! Check if this (kpoint, spin) was already calculated
   if (all(sigma%qp_done(ikcalc, :) == 1)) cycle
   call cwtime(cpu_ks, wall_ks, gflops_ks, "start")
   !call abimem_report(std_out, with_mallinfo=.False.)
   !write(std_out, *)"xmpi_count_requests", xmpi_count_requests

   ! Find IBZ(k) for q-point integration.
   call cwtime(cpu_setk, wall_setk, gflops_setk, "start")
   call sigma%setup_kcalc(dtset, cryst, dvdb, ebands, ikcalc, dtset%prtvol, comm)

   ! Symmetry indices for kk.
   kk = sigma%kcalc(:, ikcalc)
   ik_ibz = sigma%kcalc2ibz(ikcalc, 1); isym_k = sigma%kcalc2ibz(ikcalc, 2)
   trev_k = sigma%kcalc2ibz(ikcalc, 6); g0_k = sigma%kcalc2ibz(ikcalc, 3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   ABI_CHECK(isirr_k, "For the time being the k-point must be in the IBZ")
   kk_ibz = ebands%kptns(:,ik_ibz)
   npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)

   ! Allocate PW-arrays. Note mpw in kg_kq
   ABI_MALLOC(kg_k, (3, npw_k))
   kg_k = wfd%kdata(ik_ibz)%kg_k
   ABI_MALLOC(kg_kq, (3, mpw))

   ! Spherical Harmonics for useylm==1.
   ABI_MALLOC(ylm_k,(mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylm_kq,(mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

   call cwtime_report(" Setup kcalc", cpu_setk, wall_setk, gflops_setk)

   ! TODO: Spin should be treated in a more flexible and scalable way --> kcalc and bdgw should depend on spin.
   ! Introduce other comm and cartesian dimension for spin
   do spin=1,nsppol
     ! Check if this kpoint and spin was already calculated
     if (sigma%qp_done(ikcalc, spin) == 1) cycle

     call timab(1900, 1, tsec)
     ! Bands in Sigma_nk to compute and number of bands in sum over states.
     bstart_ks = sigma%bstart_ks(ikcalc, spin)
     nbcalc_ks = sigma%nbcalc_ks(ikcalc, spin)
     bsum_start = sigma%bsum_start; bsum_stop = sigma%bsum_stop
     nbsum = sigma%nbsum
     !nband_k = ebands%nband(ik_ibz + (spin-1) * ebands%nkpt)

     ! Zero self-energy matrix elements. Build frequency mesh for nk states.
     sigma%vals_e0ks = zero; sigma%dvals_de0ks = zero; sigma%dw_vals = zero
     if (sigma%calc_mrta) then
       sigma%linewidth_mrta = zero
       ABI_MALLOC(alpha_mrta, (nbcalc_ks))
     end if

     ! Prepare computation of Sigma_{nk}(w) and spectral function.
     if (sigma%nwr > 0) then
       sigma%vals_wr = zero
       do ib_k=1,nbcalc_ks
         band_ks = ib_k + bstart_ks - 1
         ! Build linear mesh **centered** around the KS energy.
         eig0nk = ebands%eig(band_ks, ik_ibz, spin) - sigma%wr_step * (sigma%nwr / 2)
         sigma%wrmesh_b(:,ib_k) = arth(eig0nk, sigma%wr_step, sigma%nwr)
       end do
     end if

     ! Prepare Eliasberg function.
     if (sigma%gfw_nomega > 0) then
       ABI_SFREE(sigma%gf_nnuq)
       ABI_CALLOC(sigma%gf_nnuq, (nbcalc_ks, natom3, sigma%nqibz_k, 3))
     end if
     if (dtset%prteliash == 3) sigma%a2few = zero

     ! Allocate eph matrix elements.
     ABI_MALLOC(gkq_atm, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkq_nu, (2, nbcalc_ks, natom3))

     ! Arrays for Debye-Waller
     if (.not. sigma%imag_only) then
       ABI_CALLOC_OR_DIE(gkq0_atm, (2, nbcalc_ks, bsum_start:bsum_stop, natom3), ierr)
     end if

     ! Integrate delta functions inside miniBZ around Gamma.
     if (sigma%frohl_model == 1 .and. sigma%imag_only) then
       call eval_sigfrohl_deltas(sigma, cryst, ifc, ebands, ikcalc, spin, dtset%prtvol, comm)
     end if

     ABI_CALLOC(gkq2_lr, (nbcalc_ks, natom3))
     if (sigma%frohl_model == 2) then
       ! Compute self-energy with Frohlich model for the gkq to handle divergence for q-->0
       ! and improve convergence with respect to nqpt.
       call eval_sigfrohl2(sigma, cryst, ifc, ebands, ikcalc, spin, comm)
     end if

     ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
     ! and save KS energies in sigma%e0vals
     ! Note: One should rotate the wavefunctions if kk is not irred (not implemented)
     ABI_MALLOC(kets_k, (2, npw_k*nspinor, nbcalc_ks))
     ABI_MALLOC(sigma%e0vals, (nbcalc_ks))
     do ib_k=1,nbcalc_ks
       band_ks = ib_k + bstart_ks - 1
       call wfd%copy_cg(band_ks, ik_ibz, spin, kets_k(1, 1, ib_k))
       sigma%e0vals(ib_k) = ebands%eig(band_ks, ik_ibz, spin)
     end do

     ! Distribute q-points, compute tetra weigths.
     call sigmaph_setup_qloop(sigma, dtset, cryst, ebands, dvdb, spin, ikcalc, nfftf, ngfftf, comm)

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq, spin, vlocal=vlocal, with_nonlocal=.true.)
     call timab(1900, 2, tsec)

     ! =======================================
     ! Integration over q-points in the IBZ(k)
     ! =======================================
     call cwtime(cpu_qloop, wall_qloop, gflops_qloop, "start")
     ignore_kq = 0; ignore_ibsum_kq = 0

     do imyq=1,sigma%my_nqibz_k
       iq_ibz = sigma%myq2ibz_k(imyq)
       qpt = sigma%qibz_k(:, iq_ibz)
       !qibz = sigma%qibz(:, sigma%ind_ibzk2ibz(1, iq_ibz))
       isqzero = sum(qpt**2) < tol14 !; if (isqzero) cycle

       call cwtime(cpu, wall, gflops, "start")

       ! Find k+q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       kq = kk + qpt
       ikq_ibz = sigma%indkk_kq(1, iq_ibz); isym_kq = sigma%indkk_kq(2, iq_ibz)
       trev_kq = sigma%indkk_kq(6, iq_ibz); g0_kq = sigma%indkk_kq(3:5, iq_ibz)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0)) !; isirr_kq = .True.
       kq_ibz = ebands%kptns(:, ikq_ibz)
       nband_kq = ebands%nband(ikq_ibz + (spin-1) * ebands%nkpt)
       if (.not. ihave_ikibz_spin(ikq_ibz, spin)) then
         ignore_kq = ignore_kq + 1; cycle
       end if

       if (sigma%use_ftinterp) then
         ! Use Fourier interpolation to get DFPT potentials for this qpt (hopefully in cache).
         db_iqpt = sigma%ind_ibzk2ibz(1, iq_ibz)
         qq_ibz = sigma%qibz(:, db_iqpt)
         call dvdb%get_ftqbz(cryst, qpt, qq_ibz, sigma%ind_ibzk2ibz(:, iq_ibz), cplex, nfftf, ngfftf, v1scf, sigma%comm_pert)
       else
         ! Read and reconstruct the dvscf potentials for qpt and my_npert perturbations.
         ! This call allocates v1scf(cplex, nfftf, nspden, my_npert))
         db_iqpt = sigma%indq2dvdb_k(1, iq_ibz)
         ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB file."))
         call dvdb%readsym_qbz(cryst, qpt, sigma%indq2dvdb_k(:,iq_ibz), cplex, nfftf, ngfftf, v1scf, sigma%comm_pert)
       end if

       ! Get phonon frequencies and displacements in reduced coordinates for this q-point
       call timab(1901, 1, tsec)
       call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red, comm=sigma%comm_pert)

       if (sigma%frohl_model == 2) then
         ! TODO: Recheck this part
         qpt_cart = matmul(cryst%rprimd, qpt)
         qmod = sqrt(sum(qpt_cart ** 2))
         inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))
         fqdamp = (two / cryst%ucvol) ** 2 * inv_qepsq ** 2 !* exp(-(qmod/sigma%qdamp) ** 2)

         ! Compute gkq_{LR} without
         gkq2_lr = zero
         do nu=1,natom3
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
           cnum = zero
           do iatom=1,cryst%natom
             cdd = cmplx(displ_cart(1,:,iatom, nu), displ_cart(2,:,iatom, nu), kind=dpc) * &
                 exp(-j_dpc * two_pi * dot_product(qpt, cryst%xred(:, iatom)))
             cnum = cnum + dot_product(qpt_cart, matmul(ifc%zeff(:, :, iatom), cdd))
           end do
           do ib_k=1,nbcalc_ks
             gkq2_lr(ib_k, nu) = (real(cnum) ** 2 + aimag(cnum) ** 2) / (two * wqnu)
           end do
         end do ! nu
         gkq2_lr = gkq2_lr * fqdamp
       end if ! frohl_model

       ! Double grid stuff
       if (sigma%use_doublegrid) then
         call sigma%eph_doublegrid%get_mapping(kk, kq, qpt)
         iq_bz_frohl = sigma%eph_doublegrid%get_index(qpt,2)
         iq_ibz_frohl = sigma%eph_doublegrid%bz2ibz_dense(iq_bz_frohl)
       end if

       ! Map q to qibz for tetrahedron
       if (sigma%qint_method > 0) then
         if (.not. sigma%use_doublegrid) then
           iq_ibz_fine = iq_ibz
           if (sigma%symsigma == 0) iq_ibz_fine = sigma%ephwg%lgk%find_ibzimage(qpt)
           ABI_CHECK(iq_ibz_fine /= -1, sjoin("Cannot find q-point in IBZ(k)", ktoa(qpt)))
           if (sigma%symsigma == 1) then
              if (.not. all(abs(sigma%qibz_k(:, iq_ibz_fine) - sigma%ephwg%lgk%ibz(:, iq_ibz_fine)) < tol12)) then
                MSG_ERROR("Mismatch in qpoints.")
              end if
           end if
         endif
       end if

       ! Get npw_kq, kg_kq for k+q. Be careful with time-reversal symmetry and istwf_kq
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
       call timab(1901, 2, tsec)
       call timab(1902, 1, tsec)

       istwf_kqirr = wfd%istwfk(ikq_ibz); npw_kqirr = wfd%npwarr(ikq_ibz)
       ABI_MALLOC(bra_kq, (2, npw_kq*nspinor))
       ABI_MALLOC(cgwork, (2, npw_kqirr*nspinor))

       ! Allocate array to store H1 |psi_nk> for all 3*natom perturbations
       ABI_MALLOC_OR_DIE(h1kets_kq, (2, npw_kq*nspinor, my_npert, nbcalc_ks), ierr)

       ! Allocate vlocal1 with correct cplex. Note nvloc
       ABI_MALLOC_OR_DIE(vlocal1,(cplex*n4,n5,n6, gs_hamkq%nvloc, my_npert), ierr)

       ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))
       ABI_MALLOC(gvnlx1, (2, npw_kq*nspinor))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm==1) then
       !   optder = 0; if (useylmgr == 1) optder = 1
       !   call initylmg(cryst%gprimd, kg_k, kk, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1, &
       !     [npw_k], dtset%nsppol, optder, cryst%rprimd, ylm_k, ylmgr)
       !   call initylmg(cryst%gprimd, kg_kq, kq, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1, &
       !     [npw_kq], dtset%nsppol, optder, cryst%rprimd, ylm_kq, ylmgr_kq)
       !end if

       if (sigma%calc_mrta) then
         call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kq, istwf_kq, npw_kq, kg_kq)
       end if

       if (dtset%eph_stern == 1 .and. .not. sigma%imag_only) then
         ABI_CALLOC(cg1s_kq, (2, npw_kq*nspinor, natom3, nbcalc_ks))
       end if

       ! Loop over all 3*natom perturbations (Each CPU prepares its own potentials)
       ! In the inner loop, I calculate H1 * psi_k, stored in h1kets_kq on the k+q sphere.
       do imyp=1,my_npert
         idir = sigma%my_pinfo(1, imyp); ipert = sigma%my_pinfo(2, imyp); ipc = sigma%my_pinfo(3, imyp)

         ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
         ! Each CPU prepares its own potentials.
         call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc, &
           pawfgr,mpi_enreg,vtrial,v1scf(:,:,:,imyp),vlocal,vlocal1(:,:,:,:,imyp))

         ! This is needed to call dfpt_cgwf (Sternheimer).
         !if (dtset%eph_stern == 1 .and. .not. sigma%imag_only) then
         call load_spin_hamiltonian(gs_hamkq, spin, vlocal=vlocal, with_nonlocal=.true.)
         !end if

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)
         call load_spin_rf_hamiltonian(rf_hamkq, spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert, &  ! In
           cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k, &          ! In
           npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq, &         ! In
           dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)      ! Out

         ! Compute H(1) applied to GS wavefunction Psi_nk(0)
         do ib_k=1,nbcalc_ks
           band_ks = ib_k + bstart_ks - 1
           eig0nk = ebands%eig(band_ks, ik_ibz, spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0,kets_k(:,:,ib_k),cwaveprj0,h1kets_kq(:, :, imyp, ib_k),&
             grad_berry,gs1c,gs_hamkq,gvnlx1,idir,ipert,eshift,mpi_enreg,optlocal,&
             optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
         end do

         ! Activate Sternheimer (we are still inside the MPI loop over my_npert).
         ! NB: Assume adiabatic AHC expression to compute the contribution of states above nband_kq.
         if (dtset%eph_stern == 1 .and. .not. sigma%imag_only) then
           mcgq = npw_kq * nspinor * nband_kq
           mgscq = npw_kq * nspinor * nband_kq * psps%usepaw
           ABI_CALLOC(cgq, (2, npw_kq * nspinor, nband_kq))
           ABI_MALLOC(gscq, (2, npw_kq * nspinor, nband_kq*psps%usepaw))

           ! Build global array with cg_kq wavefunctions to prepare call to dfpt_cgwf.
           ! TODO: Ideally, dfpt_cgwf should be modified so that we can pass cgq that is MPI distributed over nband_kq
           ! bsum_range is not compatible with Sternheimer. There's a check at the level of the parser in chkinp.
           do ibsum_kq=sigma%my_bsum_start, sigma%my_bsum_stop
             if (isirr_kq) then
                call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, bra_kq)
              else
                ! Reconstruct u_kq(G) from the IBZ image.
                call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, cgwork)
                call cgtk_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1, &
                                 npw_kqirr, wfd%kdata(ikq_ibz)%kg_k, &
                                 npw_kq, kg_kq, istwf_kqirr, istwf_kq, cgwork, bra_kq, work_ngfft, work)
              end if
              cgq(:, :, ibsum_kq) = bra_kq
           end do
           call xmpi_sum(cgq, sigma%comm_bsum, ierr)

           ! In principle, one can provide an improved initial guess if the WFK contains > nband states.
           ABI_CALLOC(out_eig1_k, (2*nband_kq**2))
           ABI_MALLOC(dcwavef, (2, npw_kq*nspinor*usedcwavef0))
           ABI_MALLOC(gh1c_n, (2, npw_kq*nspinor))
           ABI_MALLOC(ghc, (2, npw_kq*nspinor))
           ABI_MALLOC(gsc, (2, npw_kq*nspinor))
           ABI_MALLOC(gvnlxc, (2, npw_kq*nspinor))

           ! TODO: grad_berry is problematic because in dfpt_cgwf, the array is declared with
           !
           !  real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband)
           !
           ! and
           !
           !  npw1_k=number of plane waves at this k+q point
           !
           ! So in principle we should allocate lot of memory to avoid bound checking error!
           ! For the time being use mpw1 = 0 because mpw1 is not used in this call to dfpt_cgwf
           ! still it's clear that the treatment of this array must be completely refactored in the DFPT code.
           !
           grad_berry_size_mpw1 = 0
           do ib_k=1,nbcalc_ks
             ! MPI parallelism inside comm_bsum
             ! TODO: To be replaced by MPI parallellism in projbd inside dfpt_cgwf
             if (mod(ib_k, sigma%nprocs_bsum) /= sigma%me_bsum) cycle
             band_ks = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band_ks, ik_ibz, spin)

             nline_in = min(100, npw_kq); if (dtset%nline > nline_in) nline_in = min(dtset%nline, npw_kq)
             nlines_done = 0

             ! Save u1 in cg1s_kq
             cg1s_kq(:, :, ipc, ib_k) = zero

             call dfpt_cgwf(band_ks, berryopt0, cgq, cg1s_kq(:,:,ipc, ib_k), kets_k(:,:,ib_k), &
               cwaveprj, cwaveprj0, rf2, dcwavef, &
               eig0nk, ebands%eig(:, ikq_ibz, spin), out_eig1_k, ghc, gh1c_n, grad_berry, gsc, gscq, &
               gs_hamkq, gvnlxc, gvnlx1, icgq0, idir, ipert, igscq0, &
               mcgq, mgscq, mpi_enreg, grad_berry_size_mpw1, cryst%natom, nband_kq, nbdbuf0, nline_in, npw_k, npw_kq, nspinor, &
               opt_gvnlx1, dtset%prtvol, quit0, out_resid, rf_hamkq, dtset%dfpt_sciss, -one, dtset%tolwfr, &
               !opt_gvnlx1, -15, quit0, out_resid, rf_hamkq, dtset%dfpt_sciss, dtset%tolrde, dtset%tolwfr, &
               usedcwavef0, dtset%wfoptalg, nlines_done)

             ! Handle possible convergence issues.
             if (out_resid > dtset%tolwfr) then
               write(msg, "(2(a, es13.5), a, i0, a)") &
                 "Sternheimer solver didn't convergence: out_resid:", out_resid, " >= tolwfr: ", dtset%tolwfr, &
                 " after nline: ", nlines_done, " iterations. Increase nline and/or tolwfr."
               MSG_ERROR(msg)
             else if (out_resid < zero) then
               MSG_ERROR(sjoin(" out_resid: ", ftoa(out_resid), ", nlines_done:", itoa(nlines_done)))
             else
               if (my_rank == master .and. dtset%prtvol > 10) then
                 write(std_out,*)" Converged with out_resid: ", out_resid, " after ", nlines_done, " iterations."
                 write(std_out,*)" |psi1|^2", cg_real_zdotc(npw_kq*nspinor, cg1s_kq(:, :, ipc, ib_k), cg1s_kq(:, :, ipc, ib_k))
               end if
             end if

           end do ! ib_k

           ABI_FREE(cgq)
           ABI_FREE(gscq)
           ABI_FREE(out_eig1_k)
           ABI_FREE(dcwavef)
           ABI_FREE(gh1c_n)
           ABI_FREE(ghc)
           ABI_FREE(gsc)
           ABI_FREE(gvnlxc)
         end if ! sternheimer

         call destroy_rf_hamiltonian(rf_hamkq)

         ABI_FREE(kinpw1)
         ABI_FREE(kpg1_k)
         ABI_FREE(kpg_k)
         ABI_FREE(dkinpw)
         ABI_FREE(ffnlk)
         ABI_FREE(ffnl1)
         ABI_FREE(ph3d)
         ABI_SFREE(ph3d1)
       end do ! imyp
       call timab(1902, 2, tsec)

       ABI_FREE(gs1c)
       ABI_FREE(gvnlx1)
       ABI_FREE(vlocal1)
       ABI_FREE(v1scf)

       ! Add contribution to Fan-Migdal self-energy coming from Sternheimer.
       if (dtset%eph_stern == 1 .and. .not. sigma%imag_only) then
         call xmpi_sum(cg1s_kq, comm, ierr)
         ! h1kets_kq are MPI distributed inside comm_pert but we need off-diagonal pp' terms --> collect results.
         ABI_CALLOC(h1kets_kq_allperts, (2, npw_kq*nspinor, natom3, nbcalc_ks))
         do imyp=1,my_npert
           ipc = sigma%my_pinfo(3, imyp)
           h1kets_kq_allperts(:, :, ipc, :) = h1kets_kq(:, :, imyp, :)
         end do
         call xmpi_sum(h1kets_kq_allperts, sigma%comm_pert, ierr)

         if (isqzero) then
           ABI_CALLOC(stern_dw, (2, natom3, natom3, nbcalc_ks))
         end if

         ! Compute S_pp' = <D_{qp} vscf u_nk|u'_{nk+q p'}>
         ABI_CALLOC(stern_ppb, (2, natom3, natom3, nbcalc_ks))
         do ib_k=1,nbcalc_ks
           call cg_zgemm("C", "N", npw_kq*nspinor, natom3, natom3, &
             h1kets_kq_allperts(:,:,:,ib_k), cg1s_kq(:,:,:,ib_k), stern_ppb(:,:,:,ib_k))
             !cg1s_kq(:,:,:,ib_k), h1kets_kq_allperts(:,:,:,ib_k), stern_ppb(:,:,:,ib_k))

             ! Save data for Debye-Waller computation (performed outside the q-loop)
             if (isqzero) stern_dw(:,:,:,ib_k) = stern_ppb(:,:,:,ib_k)
         end do

         ! Compute contribution for M > sigma%nbsum using static + adiabatic approximation for self-energy.
         do nu=1,natom3
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
           ! Get phonon occupation for all temperatures.
           nqnu_tlist = occ_be(wqnu, sigma%kTmesh(:), zero)

           do ib_k=1,nbcalc_ks
             ! sum_{pp'} d_p* Stern_{pp'} d_p' with d = displ_red(:, :, :, nu), S = stern_ppb(:, :, :, ib_k)
             vec_natom3 = zero
             call cg_zgemm("N", "N", natom3, natom3, 1, stern_ppb(:,:,:,ib_k), displ_red(:,:,:,nu), vec_natom3)
             dotri = cg_zdotc(natom3, displ_red(:,:,:,nu), vec_natom3)
             !write(std_out, *)"dotri:", dotri
             rfact = dotri(1)
             !rfact = cg_real_zdotc(natom3, displ_red(:,:,:,nu), vec_natom3)
             rfact = rfact * sigma%wtq_k(iq_ibz) / (two * wqnu)
             ! Need to rescale by nprocs because for the time being all procs enter this part (DOH!)
             rfact = rfact / nprocs

             do it=1,sigma%ntemp
               sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + (two * nqnu_tlist(it) + one) * rfact
               !if (sigma%nwr > 0) sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + gkq2 * cfact_wr(:)

               ! Add DW contribution for this T
               if (isqzero) then
               !  gdw2_mn(ibsum, ib_k) = - gdw2_mn(ibsum, ib_k) /  (four * two * wqnu)
               !  sigma%dw_vals(it, ib_k) = sigma%dw_vals(it, ib_k) + rfact
               !  sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + rfact
               !  !if (sigma%nwr > 0) sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + rfact
               end if

               ! TODO Eliashberg functions.
               !if (dtset%prteliash /= 0) then
               !end if
             end do
           end do
         end do

         ABI_FREE(cg1s_kq)
         ABI_FREE(h1kets_kq_allperts)
         ABI_FREE(stern_ppb)
       end if ! eph_stern == 1

       ! ==============
       ! Sum over bands
       ! ==============
       call timab(1903, 1, tsec)
       do ibsum_kq=sigma%my_bsum_start, sigma%my_bsum_stop

         ! This can happen if we have loaded the wavefunctions inside the energy range.
         if (sigma%imag_only .and. sigma%qint_method == 1) then
           if (.not. wfd%ihave_ug(ibsum_kq, ikq_ibz, spin)) then
             ignore_ibsum_kq = ignore_ibsum_kq + 1; cycle
           end if
         end if

         ! Symmetrize wavefunctions from IBZ (if needed). Be careful with time-reversal symmetry.
         if (isirr_kq) then
           ! Copy u_kq(G)
           call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, bra_kq)
         else
           ! Reconstruct u_kq(G) from the IBZ image.
           ! Use cgwork as workspace array, results stored in bra_kq
           ! g0_kq = g0ibz_kq + g0bz_kq
           call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, cgwork)
           call cgtk_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1, &
                            npw_kqirr, wfd%kdata(ikq_ibz)%kg_k, &
                            npw_kq, kg_kq, istwf_kqirr, istwf_kq, cgwork, bra_kq, work_ngfft, work)
         end if

         ! Get gkk(kcalc, q, idir-ipert) (atom representation)
         ! No need to handle istwf_kq because it's always 1.
         gkq_atm = zero
         do imyp=1,my_npert
           ipc = sigma%my_pinfo(3, imyp)
           ! Calculate <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)> for this pert (NC psps) istwf_k always 1 ...
           do ib_k=1,nbcalc_ks
             gkq_atm(:, ib_k, ipc) = cg_zdotc(npw_kq*nspinor, bra_kq, h1kets_kq(:,:,imyp,ib_k))
           end do
           !call cg_zgemv("C", npw_kq*nspinor, nbcalc_ks, h1kets_kq(:,:,:,imyp), bra_kq, gkq_atm(:,:,ipc))
         end do
         !ii = nbcalc_ks * my_npert
         !call cg_zgemm("H", "N", npw_kq*nspinor, ii, ii, h1kets_kq, bra_kq, gkq_atm)
         !call cg_zgemm("H", "N", npw_kq*nspinor, ii, ii, bra_kq, h1kets_kq, gkq_atm)

         ! Get gkk(kcalc, q, nu)
         if (sigma%nprocs_pert > 1) call xmpi_sum(gkq_atm, sigma%comm_pert, ierr)
         call gkknu_from_atm(1, nbcalc_ks, 1, natom, gkq_atm, phfrq, displ_red, gkq_nu)

         ! Save data for Debye-Waller computation (performed outside the q-loop)
         ! gkq_nu(2, nbcalc_ks, bsum_start:bsum_stop, natom3)
         if (isqzero .and. .not. sigma%imag_only) then
           gkq0_atm(:, :, ibsum_kq, :) = gkq_atm
         end if

         ! Accumulate contribution to self-energy
         eig0mkq = ebands%eig(ibsum_kq, ikq_ibz, spin)
         ! q-weight for naive integration
         weight_q = sigma%wtq_k(iq_ibz)

         if (sigma%calc_mrta) then
           call ddkop%apply(eig0mkq, npw_kq, wfd%nspinor, bra_kq, cwaveprj0, wfd%mpi_enreg)
           vkq_red = ddkop%get_velocity(eig0mkq, istwf_kq, npw_kq, wfd%nspinor, wfd%mpi_enreg%me_g0, bra_kq)
           call ddk_red2car(cryst%rprimd,vkq_red,vkq)
           do ib_k=1,nbcalc_ks
             vk(1,:) = sigma%vcar_calc(:, ib_k, ikcalc, spin)
             vkk_norm2 = dot_product(vk(1,:), vk(1,:))
             alpha_mrta(ib_k) = zero
             if (vkk_norm2 > tol6) alpha_mrta(ib_k) = one - dot_product(vkq(1,:), vk(1,:)) / vkk_norm2**2
           end do
         end if

         do imyp=1,my_npert
           ! Ignore acoustic or unstable modes.
           nu = sigma%my_pinfo(3, imyp)
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
           if (sigma%phmodes_skip(nu)==1) cycle

           ! For each band in Sigma_{bk}
           do ib_k=1,nbcalc_ks
             band_ks = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band_ks, ik_ibz, spin)
             gkq2 = weight_q * (gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2)
             ediff = eig0nk - eig0mkq

             ! Handle intra-band term in polar materials
             ! DFPT gkk for q == 0 should not diverge (should plot gkq2(q) though)
             !if (sigma%frohl_model == 2 .and. abs(ediff) <= TOL_EDIFF) then
             !  gkq2 = gkq2 - weight_q * gkq2_lr(ib_k, nu)
             !end if

             ! Optionally, accumulate contribution to Eliashberg functions
             if (dtset%prteliash /= 0) then
               ! EPH strength with delta(en - emkq)
               rfact = gaussian(eig0nk - eig0mkq, dtset%tsmear)
               sigma%gf_nnuq(ib_k, nu, iq_ibz, 1) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 1) + &
                    rfact * (gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2)

               ! Treat Fan term.
               if (ediff > wqnu) then
                  rfact = one / ediff
               else
                 ! Non adiabatic regime --> Add complex shift.
                 ! Note however that the expression for this flavor of Eliashberg function relies on adiabaticity.
                 rfact = real(one/(ediff + sigma%ieta))
               end if

               gf_val = gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2
               if (sigma%frohl_model == 1 .and. isqzero .and. ibsum_kq == band_ks) then
                 gf_val = frohl_sphcorr(nu) * (four_pi / three * q0rad ** 3)
               end if
               sigma%gf_nnuq(ib_k, nu, iq_ibz, 2) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 2) + gf_val * rfact

               if (dtset%prteliash == 3) then
                 delta_e_minus_emkq = gaussian(sigma%a2f_emesh - eig0mkq, dtset%tsmear)
                 dwargs = sigma%gfw_mesh - phfrq(nu)
                 dtw_weights(:, 1) = gaussian(dwargs, dtset%ph_smear)
                 do iw=1,sigma%gfw_nomega
                   sigma%a2few(:, iw, ib_k) = sigma%a2few(:, iw, ib_k) + &
                      delta_e_minus_emkq(:) * dtw_weights(iw, 1) * gf_val * sigma%wtq_k(iq_ibz)
                 end do
               end if
             end if  ! dtset%prteliash /= 0

             do it=1,sigma%ntemp
               ! Compute occ for this T (note mu_e(it) Fermi level)
               nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)
               !f_nk = occ_fd(eig0nk, sigma%kTmesh(it), sigma%mu_e(it))
               f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

               ! Here we have to handle 3 different logical values leading to 9 different cases:
               !
               ! qint_method         0      1
               !   use_doublegrid   .true. .false.
               !     imag_only      .true. .false.
               !
               ! We will write this with nested conditionals using the order above

               if (sigma%qint_method == 0) then
                 ! =========
                 ! zcut mode
                 ! =========
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
                   simag = gkq2 * aimag(cfact)
                   sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + j_dpc * simag
                   if (sigma%calc_mrta) then
                     sigma%linewidth_mrta(it, ib_k) = sigma%linewidth_mrta(it, ib_k) + simag * alpha_mrta(ib_k)
                   end if

                 else
                   sig_cplx = gkq2 * cfact
                   if (sigma%frohl_model == 1 .and. isqzero .and. ediff <= TOL_EDIFF) then
                     ! Treat Frohlich divergence with spherical integration around the Gamma point.
                     ! In principle one should rescale by the number of degenerate states but it's
                     ! easier to move all the weight to a single band.
                     sig_cplx = czero
                     if (ibsum_kq == band_ks) sig_cplx = frohl_sphcorr(nu) * (two * f_mkq - one)
                   end if

                   sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + sig_cplx
                 end if

               else
                 ! ===================
                 ! Tetrahedron method
                 ! ===================
                 if (sigma%use_doublegrid) then
                   do jj=1,sigma%eph_doublegrid%ndiv
                     ! Double Grid shared points weights
                     ikq_bz_fine  = sigma%eph_doublegrid%mapping(2, jj)
                     weight = sigma%eph_doublegrid%weights_dense(ikq_bz_fine)

                     ! Electronic eigenvalue
                     ikq_ibz_fine = sigma%eph_doublegrid%mapping(5, jj)
                     eig0mkq = sigma%eph_doublegrid%ebands_dense%eig(ibsum_kq, ikq_ibz_fine, spin)
                     f_mkq = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))

                     ! Phonon frequency
                     iq_ibz_fine = sigma%eph_doublegrid%mapping(6, jj)
                     wqnu = sigma%ephwg%phfrq_ibz(iq_ibz_fine,nu)
                     nqnu = occ_be(wqnu, sigma%kTmesh(it), zero)

                     ! Add Frohlich contribution
                     gkq2_pf = gkq2
                     if (ediff <= TOL_EDIFF .and. sigma%frohl_model == 3) then
                       gkq2_dfrohl = sigma%ephwg%frohl_ibz(iq_ibz_fine,nu) - sigma%ephwg%frohl_ibz(iq_ibz_frohl,nu)
                       gkq2_pf = gkq2_pf + weight_q*gkq2_dfrohl
                     end if

                     if (sigma%imag_only) then
                       ! Note pi factor from Sokhotski-Plemelj theorem.
                       simag = gkq2_pf * pi * ( &
                         (nqnu + f_mkq      ) * sigma%deltaw_pm(1, ib_k, imyp, ibsum_kq, imyq, jj) +  &
                         (nqnu - f_mkq + one) * sigma%deltaw_pm(2, ib_k, imyp, ibsum_kq, imyq, jj) ) * weight
                       sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + j_dpc * simag
                       if (sigma%calc_mrta) then
                         sigma%linewidth_mrta(it, ib_k) = sigma%linewidth_mrta(it, ib_k) + simag * alpha_mrta(ib_k)
                       end if
                     else
                       sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + gkq2_pf * ( &
                         (nqnu + f_mkq      ) * sigma%cweights(1, 1, ib_k, imyp, ibsum_kq, imyq, jj) +  &
                         (nqnu - f_mkq + one) * sigma%cweights(1, 2, ib_k, imyp, ibsum_kq, imyq, jj) ) * weight
                     end if
                   end do

                 else
                   ! Tetrahedron method WITHOUT double grid.
                   if (sigma%imag_only) then
                     simag = gkq2 * pi * ( &
                       (nqnu + f_mkq      ) * sigma%deltaw_pm(1, ib_k, imyp, ibsum_kq, imyq, 1) +  &
                       (nqnu - f_mkq + one) * sigma%deltaw_pm(2, ib_k, imyp, ibsum_kq, imyq, 1) )

                     if (sigma%frohl_model == 1 .and. isqzero .and. ediff <= TOL_EDIFF) then
                       ! Treat Frohlich divergence with spherical integration of deltas around the Gamma point.
                       ! In principle one should rescale by the number of degenerate states but it's
                       ! easier to move all the weight to a single band
                       simag = zero
                       ! TODO: Check the sign, use retarded convention.
                       if (ibsum_kq == band_ks) simag = -pi * sum(sigma%frohl_deltas_sphcorr(1:2, it, ib_k, nu), dim=1)
                     end if

                     sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + j_dpc * simag
                     if (sigma%calc_mrta) then
                       sigma%linewidth_mrta(it, ib_k) = sigma%linewidth_mrta(it, ib_k) + simag * alpha_mrta(ib_k)
                     end if

                   else
                     sig_cplx = gkq2 * ( &
                       (nqnu + f_mkq      ) * sigma%cweights(1, 1, ib_k, imyp, ibsum_kq, imyq, 1) +  &
                       (nqnu - f_mkq + one) * sigma%cweights(1, 2, ib_k, imyp, ibsum_kq, imyq, 1) )
                     if (sigma%frohl_model == 1 .and. isqzero .and. ediff <= TOL_EDIFF) then
                       ! Treat Frohlich divergence with spherical integration around the Gamma point.
                       ! In principle one should rescale by the number of degenerate states but it's
                       ! easier to move all the weight to a single band
                       sig_cplx = czero
                       if (ibsum_kq == band_ks) sig_cplx = frohl_sphcorr(nu) * (two * f_mkq - one)
                     end if

                     sigma%vals_e0ks(it, ib_k) = sigma%vals_e0ks(it, ib_k) + sig_cplx
                   endif
                 end if
               end if

               ! Derivative of sigma
               ! TODO: should calculate this with the double grid as well
               if (.not. sigma%imag_only) then
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
                 if (sigma%nwr > 0) then
                   if (sigma%qint_method == 1) then
                     cfact_wr(:) = (nqnu + f_mkq      ) * sigma%cweights(2:, 1, ib_k, imyp, ibsum_kq, imyq, 1) + &
                                   (nqnu - f_mkq + one) * sigma%cweights(2:, 2, ib_k, imyp, ibsum_kq, imyq, 1)
                   else
                     cfact_wr(:) = (nqnu + f_mkq      ) / (sigma%wrmesh_b(:,ib_k) - eig0mkq + wqnu + sigma%ieta) + &
                                   (nqnu - f_mkq + one) / (sigma%wrmesh_b(:,ib_k) - eig0mkq - wqnu + sigma%ieta)
                   end if
                   sigma%vals_wr(:, it, ib_k) = sigma%vals_wr(:, it, ib_k) + gkq2 * cfact_wr(:)
                 end if

               end if

             end do ! it
           end do ! ib_k
         end do ! nu

       end do ! ibsum_kq (sum over bands at k+q)
       call timab(1903, 2, tsec)

       ABI_FREE(bra_kq)
       ABI_FREE(cgwork)
       ABI_FREE(h1kets_kq)

       if (sigma%nqibz_k < 1000 .or. (sigma%nqibz_k > 1000 .and. mod(iq_ibz, 200) == 0) .or. iq_ibz <= nprocs) then
         write(msg,'(4(a,i0),a,f8.2)') " k-point [",ikcalc,"/",sigma%nkcalc, "] q-point [",iq_ibz,"/",sigma%nqibz_k,"]"
         call cwtime_report(msg, cpu, wall, gflops)
       end if
     end do ! iq_ibz (sum over q-points in IBZ_k)

     call cwtime_report(" Q-loop", cpu_qloop, wall_qloop, gflops_qloop)

     ! Print cache stats.
     if (sigma%use_ftinterp) then
       call dvdb%ft_qcache%report_stats()
     else
       call dvdb%qcache%report_stats()
     end if

     ABI_FREE(sigma%e0vals)
     ABI_FREE(kets_k)
     ABI_FREE(gkq_atm)
     ABI_FREE(gkq_nu)
     ABI_FREE(gkq2_lr)

     ! =========================
     ! Compute Debye-Waller term
     ! =========================
     if (.not. sigma%imag_only) then
       call cwtime(cpu_dw, wall_dw, gflops_dw, "start", msg=" Computing Debye-Waller term...")
       ! Collect gkq0_atm inside comm_bq so that we can parallelize CPU easily inside comm_bq
       call xmpi_sum(gkq0_atm, sigma%comm_bq, ierr)
       ! if I have treated q = Gamma
       !call xmpi_sum(gkq0_atm, sigma%comm_bsum, ierr)
       !call xmpi_sum(gkq0_atm, sigma%comm_qpt, ierr)
       if (dtset%eph_stern == 1) call xmpi_sum(stern_dw, sigma%comm_bq, ierr)

       ABI_CALLOC(gdw2_mn, (bsum_start:bsum_stop+1, nbcalc_ks))

       ! Integral over IBZ(k). Distributed inside comm_bq
       nq = sigma%nqibz; if (sigma%symsigma == 0) nq = sigma%nqbz
       if (sigma%symsigma == +1) nq = sigma%nqibz_k
       call xmpi_split_work(nq, sigma%comm_bq, q_start, q_stop)

       do iq_ibz=q_start,q_stop
         if (abs(sigma%symsigma) == 1) then
           !qpt = sigma%qibz(:,iq_ibz); weight_q = sigma%wtq(iq_ibz)
           qpt = sigma%qibz_k(:,iq_ibz); weight_q = sigma%wtq_k(iq_ibz)
         else
           ! Use full BZ for q-integration.
           qpt = sigma%qbz(:,iq_ibz); weight_q = one / sigma%nqbz
         end if

         ! Get phonons for this q-point.
         call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red, comm=sigma%comm_pert)

         ! Sum over my modes for this q-point.
         do imyp=1,my_npert
           nu = sigma%my_pinfo(3, imyp)
           ! Ignore acoustic or unstable modes.
           wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle

           ! Get phonon occupation for all temperatures.
           nqnu_tlist = occ_be(wqnu, sigma%kTmesh(:), zero)

           ! Compute T_pp'(q,nu) matrix in reduced coordinates.
           do ip2=1,natom3
             idir2 = mod(ip2-1, 3) + 1; ipert2 = (ip2 - idir2) / 3 + 1
             do ip1=1,natom3
               idir1 = mod(ip1-1, 3) + 1; ipert1 = (ip1 - idir1) / 3 + 1
               ! (k,a) (k,a')* + (k',a) (k',a')*
               dka   = dcmplx(displ_red(1, idir1, ipert1, nu), displ_red(2, idir1, ipert1, nu))
               dkap  = dcmplx(displ_red(1, idir2, ipert1, nu), displ_red(2, idir2, ipert1, nu))
               dkpa  = dcmplx(displ_red(1, idir1, ipert2, nu), displ_red(2, idir1, ipert2, nu))
               dkpap = dcmplx(displ_red(1, idir2, ipert2, nu), displ_red(2, idir2, ipert2, nu))
               tpp_red(ip1, ip2) = dka * dconjg(dkap) + dkpa * dconjg(dkpap)
             end do
           end do

           ! Sum over bands and add (static) DW contribution for the different temperatures.
           do ibsum=bsum_start, bsum_stop
             do ib_k=1,nbcalc_ks
               ! Compute DW term following XG paper. Check prefactor.
               ! gkq0_atm(2, nbcalc_ks, bsum_start:bsum_stop, natom3)
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

               gdw2_mn(ibsum, ib_k) = gdw2_mn(ibsum, ib_k) / (four * two * wqnu)

               !if (dtset%eph_stern == 1 .and. ibsum == bsum_stop) then
               !  ! Compute DW term for m > nband
               !  cfact = zero
               !  do ip2=1,natom3
               !    do ip1=1,natom3
               !      cfact = cfact + tpp_red(ip1, ip2) * cmplx(stern_dw(1,ip1,ip2,ib_k), stern_dw(2,ip1,ip2,ib_k), kind=dpc)
               !    end do
               !  end do
               !  gdw2_mn(ibsum+1, ib_k) = real(cfact) / (four * two * wqnu)
               !end if
             end do ! ib_k
           end do ! ibsum

           do ibsum=bsum_start,bsum_stop
             eig0mk = ebands%eig(ibsum, ik_ibz, spin)

             do ib_k=1,nbcalc_ks
               band_ks = ib_k + bstart_ks - 1
               eig0nk = ebands%eig(band_ks, ik_ibz, spin)
               ! Handle n == m and degenerate states (either ignore or add broadening)
               cplx_ediff = (eig0nk - eig0mk)
               if (abs(cplx_ediff) < EPH_WTOL) cycle
               !if (abs(cplx_ediff) < EPH_WTOL) cplx_ediff = cplx_ediff + sigma%ieta

               ! Optionally, accumulate DW contribution to Eliashberg functions.
               if (dtset%prteliash /= 0) then
                  sigma%gf_nnuq(ib_k, nu, iq_ibz, 3) = sigma%gf_nnuq(ib_k, nu, iq_ibz, 3) &
                      - gdw2_mn(ibsum, ib_k) / real(cplx_ediff)
               end if

               !if (dtset%prteliash == 3) then
               !  delta_e_minus_emkq = gaussian(sigma%a2f_emesh - eig0mk, dtset%tsmear)
               !  dwargs = sigma%gfw_mesh - phfrq(nu)
               !  dtw_weights(:, 1) = gaussian(dwargs, dtset%ph_smear)
               !  do ie=1,sigma%a2f_ne
               !    sigma%a2few(:, ie, ib_k, 2) = sigma%a2few(:, ie, ib_k, 2) + &
               !         delta_e_minus_emkq(ie) * dtw_weights(:, 1) * gdw2_mn(ibsum, ib_k) / (enk - e) * sigma%wtq_k(iq_ibz)
               !end if

               ! Accumulate DW for each T, add it to Sigma(e0) and Sigma(w) as well
               ! - (2 n_{q\nu} + 1) * gdw2 / (e_nk - e_mk)
               do it=1,sigma%ntemp
                 cfact = - weight_q * gdw2_mn(ibsum, ib_k) * (two * nqnu_tlist(it) + one)  / cplx_ediff
                 if (dtset%eph_stern == 1 .and. ibsum == bsum_stop) then
                   cfact = cfact - weight_q * gdw2_mn(ibsum+1, ib_k) * (two * nqnu_tlist(it) + one)
                 end if
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
       ABI_FREE(gkq0_atm)
       ABI_SFREE(stern_dw)

       call cwtime_report(" DW completed", cpu_dw, wall_dw, gflops_dw)
     end if ! not %imag_only

     if (dtset%prteliash /= 0) then
       ! Compute Eliashberg function (useful but cost is not negligible so it must be activated by the user).
       call cwtime(cpu, wall, gflops, "start", msg=sjoin(" Computing Eliashberg function with nomega: ", &
           itoa(sigma%gfw_nomega)))

       if (dtset%prteliash == 3) call xmpi_sum(sigma%a2few, comm, ierr)

       ! Collect all terms on each node so that we can MPI-parallelize easily.
       ! NB: gf_nnuq does not include the q-weights from integration.
       call xmpi_sum(sigma%gf_nnuq, comm, ierr)
       sigma%gfw_vals = zero

       if (sigma%qint_method == 0 .or. sigma%symsigma == 0) then
         ! Eliashberg function with gaussian method (ph_smear)

         do iq_ibz=1,sigma%nqibz_k
           if (mod(iq_ibz, nprocs) /= my_rank) cycle ! MPI parallelism
           ! Recompute phonons (cannot use sigma%ephwg in this case)
           call ifc%fourq(cryst, sigma%qibz_k(:,iq_ibz), phfrq, displ_cart)
           do nu=1,natom3
             dwargs = sigma%gfw_mesh - phfrq(nu)
             dtw_weights(:, 1) = gaussian(dwargs, dtset%ph_smear)
             do ib_k=1,nbcalc_ks
               do ii=1,3
                 sigma%gfw_vals(:, ii, ib_k) = sigma%gfw_vals(:, ii, ib_k) +  &
                   sigma%gf_nnuq(ib_k, nu, iq_ibz, ii) * dtw_weights(:, 1) * sigma%wtq_k(iq_ibz)
               end do

             end do
           end do
         end do

       else
         ! Eliashberg function with tetrahedron method.
         eminmax = [sigma%gfw_mesh(1), sigma%gfw_mesh(sigma%gfw_nomega)]
         ABI_MALLOC(dt_tetra_weights, (sigma%gfw_nomega, sigma%nqibz_k, 2))
         do nu=1,natom3
           ! All procs compute weights.
           call sigma%ephwg%get_deltas_qibzk(nu, sigma%gfw_nomega, eminmax, sigma%bcorr, dt_tetra_weights, comm, &
                                             with_qweights=.True.)
           do iq_ibz=1,sigma%nqibz_k
             if (mod(iq_ibz, nprocs) /= my_rank) cycle ! MPI parallelism
             do ib_k=1,nbcalc_ks
               do ii=1,3
                 sigma%gfw_vals(:, ii, ib_k) = sigma%gfw_vals(:, ii, ib_k) +  &
                   sigma%gf_nnuq(ib_k, nu, iq_ibz, ii) * dt_tetra_weights(:, iq_ibz, 1)
               end do
             end do
           end do
         end do
         ABI_FREE(dt_tetra_weights)
       end if

       ! Collect final results.
       call xmpi_sum(sigma%gfw_vals, comm, ierr)
       call cwtime_report(" Eliashberg function", cpu, wall, gflops)
     end if

     if (my_rank == master) then
       if (ignore_kq /= 0) write(std_out, "(a, 1x, i0)")" Number of ignored k+q points:", ignore_kq
       if (ignore_ibsum_kq /= 0) write(std_out, "(a, 1x, i0)")" Number of ignored (k+q, m) states:", ignore_ibsum_kq
     end if

     ! Collect results inside comm and write results for this (k-point, spin) to NETCDF file.
     call sigma%gather_and_write(ebands, ikcalc, spin, dtset%prtvol, comm)

     ABI_SFREE(alpha_mrta)
   end do ! spin

   ABI_FREE(kg_k)
   ABI_FREE(kg_kq)
   ABI_FREE(ylm_k)
   ABI_FREE(ylm_kq)
   ABI_FREE(ylmgr_kq)

   call cwtime_report(" One ikcalc k-point", cpu_ks, wall_ks, gflops_ks)
 end do ! ikcalc

 call cwtime_report(" Sigma_eph full calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Free memory
 ABI_FREE(ihave_ikibz_spin)
 ABI_FREE(grad_berry)
 ABI_FREE(vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(nqnu_tlist)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_FREE(tpp_red)
 ABI_SFREE(cfact_wr)
 ABI_SFREE(dwargs)
 ABI_SFREE(dtw_weights)
 ABI_SFREE(delta_e_minus_emkq)

 call destroy_hamiltonian(gs_hamkq)
 call wfd%free()
 call pawcprj_free(cwaveprj0)
 ABI_DT_FREE(cwaveprj0)
 call pawcprj_free(cwaveprj)
 ABI_DT_FREE(cwaveprj)
 call ddkop%free()
 call sigma%free()

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
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkknu_from_atm(nb1, nb2, nk, natom, gkq_atm, phfrq, displ_red, gkq_nu)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb1, nb2, nk, natom
!arrays
 real(dp),intent(in) :: phfrq(3*natom),displ_red(2,3*natom,3*natom)
 real(dp),intent(in) :: gkq_atm(2,nb1,nb2,nk,3*natom)
 real(dp),intent(out) :: gkq_nu(2,nb1,nb2,nk,3*natom)

!Local variables-------------------------
!scalars
 integer :: nu,ipc

! *************************************************************************

 gkq_nu = zero

 ! Loop over phonon branches.
 do nu=1,3*natom
   ! Ignore negative or too small frequencies
   if (phfrq(nu) < EPH_WTOL) cycle

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

type(sigmaph_t) function sigmaph_new(dtset, ecut, cryst, ebands, ifc, dtfil, comm) result(new)

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
 integer,parameter :: master=0, occopt3=3, qptopt1=1, sppoldbl1=1, istwfk1=1
 integer :: my_rank,ik,my_nshiftq,my_mpw,cnt,nprocs,ik_ibz,ndeg, iq_ibz
 integer :: onpw,ii,ipw,ierr,it,spin,gap_err,ikcalc,qprange_,bstop
 integer :: jj,bstart,natom,natom3
 integer :: isym_k, trev_k, mband, i1,i2,i3
 integer :: idir, iatom, pertcase, nrest
 integer :: ip, npoints !,skw_cplex !, edos_intmeth
 logical :: downsample
 character(len=fnlen) :: wfk_fname_dense
 character(len=500) :: msg
 real(dp) :: dksqmax,ang,con,cos_phi,cos_theta,sin_phi,sin_theta,nelect, estep
 real(dp) :: cpu_all, wall_all, gflops_all, cpu, wall, gflops
 logical :: changed, isirr_k
 type(ebands_t) :: tmp_ebands, ebands_dense
 type(gaps_t) :: gaps
!arrays
 integer :: intp_kptrlatt(3,3), g0_k(3) !, skw_band_block(2)
 integer :: qptrlatt(3,3),indkk_k(1,6),my_gmax(3),band_block(2) !,kptrlatt(3,3)
 integer :: val_indeces(ebands%nkpt, ebands%nsppol), intp_nshiftk
 integer :: all_pinfo(3, cryst%natom * 3)
 integer,allocatable :: gtmp(:,:),degblock(:,:), degblock_all(:,:,:,:), ndeg_all(:,:), iperm(:)
 real(dp):: params(4), my_shiftq(3,1),kk(3),kq(3),intp_shiftk(3)
 real(dp),allocatable :: th(:),wth(:)
 integer,allocatable :: temp(:,:)
#ifdef HAVE_MPI
 integer :: ndims, comm_cart, me_cart
 logical :: reorder
 integer,allocatable :: dims(:)
 logical,allocatable :: periods(:), keepdim(:)
#endif

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call cwtime(cpu, wall, gflops, "start")

 ! Copy important dimensions.
 new%nsppol = ebands%nsppol; new%nspinor = ebands%nspinor; mband = dtset%mband
 natom = cryst%natom; natom3 = cryst%natom * 3

 ! Re-Im or Im only?
 new%imag_only = .False.; if (dtset%eph_task == -4) new%imag_only = .True.

 ! TODO: Activate by default?
 ! Compute lifetimes in the MRTA approximation
 new%calc_mrta = .False.; if (dtset%userib == 11) new%calc_mrta = .True.

 ! TODO: Remove qint_method, use eph_intmeth or perhaps dtset%qint_method dtset%kint_method
 ! FIXME: Tetra gives positive SIGE2 while zcut gives negative (retarded)
 ! Decide default behaviour for Re-Im/Im
 new%qint_method = dtset%eph_intmeth - 1

 ! Broadening parameter from zcut
 new%ieta = + j_dpc * dtset%zcut

 ! Define q-mesh for integration of the self-energy.
 ! Either q-mesh from DDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation if q not in DDB)
 new%ngqpt = dtset%dvdb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   new%ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if

 ! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems without inversion
 qptrlatt = 0; qptrlatt(1,1) = new%ngqpt(1); qptrlatt(2,2) = new%ngqpt(2); qptrlatt(3,3) = new%ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
   new%nqibz, new%qibz, new%wtq, new%nqbz, new%qbz, bz2ibz=new%ind_bz2ibz)

!HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(temp,(new%nqbz, 6))
 call listkk(dksqmax, cryst%gmet, temp, new%qibz, new%qbz, new%nqibz, new%nqbz, cryst%nsym, &
      1, cryst%symafm, cryst%symrec, 1, comm, exit_loop=.True., use_symrec=.True.)

 new%ind_bz2ibz(1,:) = temp(:,1)
 new%ind_bz2ibz(2,:) = temp(:,2)
 new%ind_bz2ibz(3,:) = temp(:,6)
 new%ind_bz2ibz(4,:) = temp(:,3)
 new%ind_bz2ibz(5,:) = temp(:,4)
 new%ind_bz2ibz(6,:) = temp(:,5)
 ABI_FREE(temp)
!END DEBUG

 ! Build (linear) mesh of K * temperatures. tsmesh(1:3) = [start, step, num]
 new%ntemp = nint(dtset%tmesh(3))
 ABI_CHECK(new%ntemp > 0, "ntemp <= 0")
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 new%kTmesh = arth(dtset%tmesh(1), dtset%tmesh(2), new%ntemp) * kb_HaK

 gap_err = get_gaps(ebands, gaps)
 if (gap_err/=0) then
   ! In case of error try to enforce semiconductor occupations before calling get_gaps
   ! This might still fail though...
   call gaps%free()
   call ebands_copy(ebands,tmp_ebands)
   call ebands_set_scheme(tmp_ebands, occopt3, dtset%tsmear, dtset%spinmagntarget, prtvol=dtset%prtvol)
   call ebands_set_nelect(tmp_ebands, tmp_ebands%nelect-dtset%eph_extrael, dtset%spinmagntarget, msg)
   gap_err = get_gaps(tmp_ebands,gaps)
   call ebands_free(tmp_ebands)
 end if
 call gaps%print(unit=std_out)
 val_indeces = get_valence_idx(ebands)

 ! Frequency mesh for sigma(w) and spectral function.
 ! TODO: Use GW variables but change default
 !dtset%freqspmin
 new%nwr = dtset%nfreqsp; new%wr_step = zero
 if (new%nwr > 0) then
   if (mod(new%nwr, 2) == 0) new%nwr = new%nwr + 1
   new%wr_step = two * eV_Ha / (new%nwr - 1)
   if (dtset%freqspmax /= zero) new%wr_step = dtset%freqspmax / (new%nwr - 1)
 end if

 ! Select k-point and bands where corrections are wanted
 ! if symsigma == +1, we have to include all degenerate states in the set
 ! because the final QP corrections will be obtained by averaging the results in the degenerate subspace.
 ! We initialize IBZ(k) here so that we have all the basic dimensions of the run and it's possible
 ! to distribuite the calculations among processors.
 new%symsigma = dtset%symsigma; new%timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 call cwtime_report(" sigmaph_new: k-points", cpu, wall, gflops)

 ! TODO: Rename variable
 ! Reorder k-points by lengths to improve cache locality?
 if (dtset%nkptgw /= 0) then
   ! Treat the k-points and bands specified in the input file via kptgw and bdgw.
   call sigtk_kcalc_from_nkptgw(dtset, mband, new%nkcalc, new%kcalc, new%bstart_ks, new%nbcalc_ks)

 else

   if (any(dtset%sigma_erange > zero)) then
     call sigtk_kcalc_from_erange(dtset, cryst, ebands, gaps, new%nkcalc, new%kcalc, new%bstart_ks, new%nbcalc_ks, comm)

   else
     ! Use qp_range to select the interesting k-points and the corresponding bands.
     !
     !    0 --> Compute the QP corrections only for the fundamental and the direct gap.
     ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
     !           bands above and below the Fermi level.
     ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
     !          Include all occupied states and `num` empty states.

     qprange_ = dtset%gw_qprange
     if (gap_err /= 0 .and. qprange_ == 0) then
       MSG_WARNING("Cannot compute fundamental and direct gap (likely metal). Will replace qprange 0 with qprange 1")
       qprange_ = 1
     end if

     if (qprange_ /= 0) then
       call sigtk_kcalc_from_qprange(dtset, cryst, ebands, qprange_, new%nkcalc, new%kcalc, new%bstart_ks, new%nbcalc_ks)
     else
       ! qprange is not specified in the input.
       ! Include direct and fundamental KS gap or include states depending on the position wrt band edges.
       call sigtk_kcalc_from_gaps(dtset, ebands, gaps, new%nkcalc, new%kcalc, new%bstart_ks, new%nbcalc_ks)
     end if
    end if

 end if ! nkptgw /= 0

 call gaps%free()

 ! The k-point and the symmetries relating the BZ k-point to the IBZ.
 ABI_MALLOC(new%kcalc2ibz, (new%nkcalc, 6))
 if (abs(new%symsigma) == 1) then
   ABI_DT_MALLOC(new%degtab, (new%nkcalc, new%nsppol))
 end if

 ! Workspace arrays used to compute degeneracy tables.
 ABI_ICALLOC(degblock_all, (2, mband, new%nkcalc, new%nsppol))
 ABI_ICALLOC(ndeg_all, (new%nkcalc, new%nsppol))
 ierr = 0

 do ikcalc=1,new%nkcalc
   if (mod(ikcalc, nprocs) /= my_rank) then
     new%kcalc2ibz(ikcalc, :) = 0
     new%bstart_ks(ikcalc, :) = 0
     new%nbcalc_ks(ikcalc, :) = 0
     cycle ! MPI parallelism.
   end if
   ! Note symrel and use_symrel.
   kk = new%kcalc(:,ikcalc)
   call listkk(dksqmax,cryst%gmet,indkk_k,ebands%kptns,kk,ebands%nkpt,1,cryst%nsym,&
      sppoldbl1,cryst%symafm,cryst%symrel,new%timrev,xmpi_comm_self,exit_loop=.True., use_symrec=.False.)

   if (dksqmax > tol12) then
      write(msg, '(4a,es16.6,7a)' )&
       "The WFK file cannot be used to compute self-energy corrections at kpoint: ",ktoa(kk),ch10,&
       "the k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
       "Q-mesh: ",trim(ltoa(new%ngqpt)),", K-mesh (from kptrlatt) ",trim(ltoa(get_diag(dtset%kptrlatt))),ch10, &
       'Action: check your WFK file and (k, q) point input variables'
      MSG_ERROR(msg)
   end if

   new%kcalc2ibz(ikcalc, :) = indkk_k(1, :)

   ik_ibz = indkk_k(1,1); isym_k = indkk_k(1,2)
   trev_k = indkk_k(1, 6); g0_k = indkk_k(1, 3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !kk_ibz = ebands%kptns(:,ik_ibz)
   if (.not. isirr_k) then
     MSG_WARNING(sjoin("For the time being the k-point must be in the IBZ but got:", ktoa(kk)))
     ierr = ierr + 1
   end if

   ! We will have to average the QP corrections over degenerate states if symsigma=1 is used.
   ! Here we make sure that all the degenerate states are included.
   ! Store also band indices of the degenerate sets, used to average final results.
   if (abs(new%symsigma) == 1) then
     cnt = 0
     do spin=1,new%nsppol
       bstop = new%bstart_ks(ikcalc,spin) + new%nbcalc_ks(ikcalc,spin) - 1
       call enclose_degbands(ebands, ik_ibz, spin, new%bstart_ks(ikcalc,spin), bstop, changed, TOL_EDIFF, degblock=degblock)
       if (changed) then
         new%nbcalc_ks(ikcalc,spin) = bstop - new%bstart_ks(ikcalc,spin) + 1
         cnt = cnt + 1
         if (cnt < 5) then
           write(msg,'(2(a,i0),2a,2(1x,i0))') &
             "Not all the degenerate states for ikcalc: ",ikcalc,", spin: ",spin,ch10, &
             "were included in the bdgw set. bdgw has been automatically changed to: ",new%bstart_ks(ikcalc,spin),bstop
           MSG_COMMENT(msg)
         end if
         write(msg,'(2(a,i0),2a)') &
           "The number of included states: ", bstop, &
           " is larger than the number of bands in the input ",dtset%nband(new%nkcalc*(spin-1)+ikcalc),ch10,&
           "Action: Increase nband."
         ABI_CHECK(bstop <= dtset%nband(new%nkcalc*(spin-1)+ikcalc), msg)
       end if

       ! Store band indices used for averaging (shifted by bstart_ks)
       ndeg = size(degblock, dim=2)
       ndeg_all(ikcalc, spin) = ndeg
       degblock_all(:, 1:ndeg, ikcalc, spin) = degblock(:, 1:ndeg)

       ABI_FREE(degblock)
     end do
   end if ! symsigma
 end do ! ikcalc
 ABI_CHECK(ierr == 0, "Fatal error, kptgw list must be in the IBZ.")

 ! Collect data
 call xmpi_sum(new%kcalc2ibz, comm, ierr)
 call xmpi_sum(new%bstart_ks, comm, ierr)
 call xmpi_sum(new%nbcalc_ks, comm, ierr)

 ! Build degtab tables.
 if (abs(new%symsigma) == 1) then
   call xmpi_sum(ndeg_all, comm, ierr)
   call xmpi_sum(degblock_all, comm, ierr)
   do ikcalc=1,new%nkcalc
     do spin=1,new%nsppol
       ndeg = ndeg_all(ikcalc, spin)
       ABI_DT_MALLOC(new%degtab(ikcalc, spin)%bids, (ndeg))
       do ii=1,ndeg
           cnt = degblock_all(2, ii, ikcalc, spin) - degblock_all(1, ii, ikcalc, spin) + 1
           ABI_DT_MALLOC(new%degtab(ikcalc, spin)%bids(ii)%vals, (cnt))
           new%degtab(ikcalc, spin)%bids(ii)%vals = [(jj, jj= &
             degblock_all(1, ii, ikcalc, spin) - new%bstart_ks(ikcalc, spin) + 1, &
             degblock_all(2, ii, ikcalc, spin) - new%bstart_ks(ikcalc, spin) + 1)]
       end do
     end do
   end do
 end if
 ABI_FREE(degblock_all)
 ABI_FREE(ndeg_all)

 call cwtime_report(" sigmaph_new: kptgw", cpu, wall, gflops)

 ! Now we can finally compute max_nbcalc
 new%max_nbcalc = maxval(new%nbcalc_ks)

 ABI_MALLOC(new%bstop_ks, (new%nkcalc, new%nsppol))
 new%bstop_ks = new%bstart_ks + new%nbcalc_ks - 1

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! used to symmetrize the wavefunctions in G-space.
 ! Note that we loop over the full BZ instead of the IBZ(k)
 ! This part is slow for very dense meshes, should try to use a geometrical approach...

 new%mpw = 0; new%gmax = 0; cnt = 0
 do ik=1,new%nkcalc
   kk = new%kcalc(:, ik)
   do i3=-1,1
     do i2=-1,1
       do i1=-1,1
         cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
         kq = kk + half * [i1, i2, i3]
         ! TODO: g0 umklapp here can enter into play gmax could not be large enough!
         call get_kg(kq, istwfk1, 1.1_dp * ecut, cryst%gmet, onpw, gtmp)
         new%mpw = max(new%mpw, onpw)
         do ipw=1,onpw
           do ii=1,3
            new%gmax(ii) = max(new%gmax(ii), abs(gtmp(ii,ipw)))
           end do
         end do
         ABI_FREE(gtmp)
       end do
     end do
   end do
 end do

 my_mpw = new%mpw; call xmpi_max(my_mpw, new%mpw, comm, ierr)
 my_gmax = new%gmax; call xmpi_max(my_gmax, new%gmax, comm, ierr)

 call wrtout(std_out, sjoin(' Optimal value of mpw:', itoa(new%mpw), "gmax:", ltoa(new%gmax)))
 call cwtime_report(" sigmaph_new: mpw", cpu, wall, gflops)

 ! Define number of bands included in self-energy summation as well as band range.
 ! This value depends on the kind of calculation as imag_only can take advantage of
 ! an energy window aroud the Fermi level.
 !
 ! Notes about MPI version.
 ! If eph_task == -4:
 !    Loops are MPI parallelized over bands so that we can distribute memory for wavefunctions over nband.
 !
 ! If eph_task == -4:
 !    Loops are MPI parallelized over q-points
 !    wavefunctions are not distributed but only states between my_bsum_start and my_bsum_stop
 !    are allocated and read from file.

 new%wmax = 1.1_dp * abs(ifc%omega_minmax(2))
 if (new%qint_method == 0) new%wmax = new%wmax + five * dtset%zcut
 ! TODO: One should be consistent with tolerances when using tetra + q-point filtering.
 !if (new%qint_method == 1) new%wmax = new%wmax + five * dtset%zcut
 !new%wmax = new%wmax + five * dtset%zcut
 !write(std_out,*)"wmax:", new%wmax * Ha_eV, " (eV)"

 new%elow = huge(one); new%ehigh = - huge(one)
 do spin=1,new%nsppol
   do ikcalc=1,new%nkcalc
     ik_ibz = new%kcalc2ibz(ikcalc, 1)
     bstart = new%bstart_ks(ikcalc, spin)
     bstop = new%bstart_ks(ikcalc, spin) + new%nbcalc_ks(ikcalc, spin) - 1
     new%ehigh = max(new%ehigh, maxval(ebands%eig(bstart:bstop, ik_ibz, spin)) + new%wmax)
     new%elow = min(new%elow, minval(ebands%eig(bstart:bstop, ik_ibz, spin)) - new%wmax)
   end do
 end do
 !call wrtout(std_out, sjoin("elow:", ftoa(elow), "ehigh:", ftoa(ehigh), "[Ha]"))

 if (new%imag_only) then

   if (all(dtset%sigma_bsum_range /= 0)) then
     new%bsum_start = max(dtset%sigma_bsum_range(1), 1)
     new%bsum_stop = min(dtset%sigma_bsum_range(2), mband)
     new%nbsum = new%bsum_stop - new%bsum_start + 1
     new%my_bsum_start = new%bsum_start; new%my_bsum_stop = new%bsum_stop
   else
     ! Compute the min/max KS energy to be included in the imaginary part.
     ! ifc%omega_minmax(2) comes froms the coarse Q-mesh of the DDB so increase it by 10%.
     ! Also take into account the Lorentzian function if zcut is used.
     ! In principle this should be large enough but it seems that the linewidths in v8[160] are slightly affected.
     ! Select indices for energy window.

     call get_bands_from_erange(ebands, new%elow, new%ehigh, new%bsum_start, new%bsum_stop)
     new%bsum_stop = min(new%bsum_stop, mband)
     ABI_CHECK(new%bsum_start <= new%bsum_stop, "bsum_start > bsum_bstop")
     new%nbsum = new%bsum_stop - new%bsum_start + 1
     new%my_bsum_start = new%bsum_start; new%my_bsum_stop = new%bsum_stop
   end if

   if (dtset%useria == 567) then
     ! Uncomment this part to use all states to debug.
     call wrtout([std_out, ab_out], "- Setting bstart to 1 and bstop to nband for debugging purposes")
     new%nbsum = mband; new%bsum_start = 1; new%bsum_stop = new%bsum_start + new%nbsum - 1
     new%my_bsum_start = new%bsum_start; new%my_bsum_stop = new%bsum_stop
   end if

 else
   ! Re + Im
   new%bsum_start = 1; new%bsum_stop = mband
   if (all(dtset%sigma_bsum_range /= 0)) then
     new%bsum_start = max(dtset%sigma_bsum_range(1), 1)
     new%bsum_stop = min(dtset%sigma_bsum_range(2), mband)
   end if
   new%nbsum = new%bsum_stop - new%bsum_start + 1
 end if

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! Init for sequential execution.
 new%comm_pert = xmpi_comm_self; new%my_npert = natom3; new%me_pert = 0; new%nprocs_pert = 1
 new%comm_bq = xmpi_comm_self; new%me_bq = 0; new%nprocs_bq = 1
 new%comm_qpt = xmpi_comm_self; new%me_qpt = 0; new%nprocs_qpt = 1
 new%comm_bsum = xmpi_comm_self; new%me_bsum = 0; new%nprocs_bsum = 1

 if (any(dtset%eph_np_pqbks /= 0)) then
   ! Use parameters from input file.
   new%nprocs_pert = dtset%eph_np_pqbks(1)
   new%nprocs_qpt  = dtset%eph_np_pqbks(2)
   new%nprocs_bsum = dtset%eph_np_pqbks(3)
   new%my_npert = natom3 / new%nprocs_pert
   ABI_CHECK(new%my_npert > 0, "nprocs_pert cannot be greater than 3 * natom.")
   ABI_CHECK(mod(natom3, new%nprocs_pert) == 0, "nprocs_pert must divide 3 * natom.")
   if (new%imag_only .and. new%nprocs_bsum /= 1) then
     MSG_ERROR("nprocs_bsum should be 1 when computing Imag(Sigma)")
   end if
 else
   ! Automatic grid generation.
   !
   ! Handle parallelism over perturbations first.
   ! Use MPI communicator to distribute 3 * natom perturbations to reduce memory requirements for DFPT potentials.
   ! Ideally, perturbations are equally distributed --> total number of CPUs should be divisible by 3 * natom.
   ! or at least, divisible by one integer i for i in [2, 3 * natom - 1].
   do cnt=natom3,2,-1
     if (mod(nprocs, cnt) == 0 .and. mod(natom3, cnt) == 0) then
       new%nprocs_pert = cnt; new%my_npert = natom3 / cnt; exit
     end if
   end do
   if (new%my_npert == natom3 .and. nprocs > 1) then
     MSG_WARNING("The number of MPI procs should be divisible by 3*natom to reduce memory requirements!")
   end if

   !do cnt=natom,2,-1
   !  if (mod(nprocs, cnt) == 0 .and. mod(natom, cnt) == 0) then
   !    new%nprocs_pert = cnt; new%my_npert = natom / cnt; exit
   !  end if
   !end do
   !if (new%my_npert == natom .and. nprocs > 1) then
   !  MSG_WARNING("The number of MPI procs should be divisible by natom to reduce memory requirements!")
   !end if

   ! Define number of procs for q-points and bands. nprocs is divisible by nprocs_pert.
   if (new%imag_only) then
     ! Just one extra MPI level for q-points.
     new%nprocs_qpt = nprocs / new%nprocs_pert
   else
     ! Try to distribute equally nbsum first.
     nrest = nprocs / new%nprocs_pert
     do bstop=nrest,1,-1
       if (mod(new%nbsum, bstop) == 0 .and. mod(nprocs, new%nprocs_pert * bstop) == 0) then
         new%nprocs_bsum = bstop; new%nprocs_qpt = nrest / new%nprocs_bsum
         exit
       end if
     end do
     ! This to recover the previous behaviour in Re-Im that is OK
     !new%nprocs_bsum = nrest; new%nprocs_qpt = 1
   end if
 end if

 ! Consistency check.
 if (new%nprocs_pert * new%nprocs_qpt * new%nprocs_bsum /= nprocs) then
   write(msg, "(a,i0,3a, 4(a,1x,i0))") &
     "Cannot create 3d Cartesian grid with nprocs: ", nprocs, ch10, &
     "Idle processes are not supported. The product of `nprocs_*` should be equal to nprocs.", ch10, &
     "nprocs_pert (", new%nprocs_pert, ") x nprocs_qpt (", new%nprocs_qpt, ") x nprocs_bsum (", new%nprocs_bsum, &
     ") = ", new%nprocs_pert * new%nprocs_qpt * new%nprocs_bsum
   MSG_ERROR(msg)
 end if

#ifdef HAVE_MPI
 ! Create 3d cartesian communicator: 3*natom perturbations, q-points in IBZ, bands in Sigma sum.
 ndims = 3
 ABI_MALLOC(dims, (ndims))
 ABI_MALLOC(periods, (ndims))
 ABI_MALLOC(keepdim, (ndims))
 periods(:) = .False.; reorder = .False.
 dims = [new%nprocs_pert, new%nprocs_qpt, new%nprocs_bsum]

 call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, comm_cart, ierr)
 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, new%coords, ierr)

 ! Create communicator to distribute perturbations.
 keepdim = .False.; keepdim(1) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, new%comm_pert, ierr)
 call MPI_COMM_RANK(new%comm_pert, new%me_pert, ierr)

 ! Create communicator for qpoints.
 keepdim = .False.; keepdim(2) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, new%comm_qpt, ierr)
 call MPI_COMM_RANK(new%comm_qpt, new%me_qpt, ierr)

 ! Create communicator for bands.
 keepdim = .False.; keepdim(3) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, new%comm_bsum, ierr)
 call MPI_COMM_RANK(new%comm_bsum, new%me_bsum, ierr)

 ! Create communicator for the (band_sum, qpoint_sum) loops
 keepdim = .False.; keepdim(2:3) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, new%comm_bq, ierr)
 call MPI_COMM_RANK(new%comm_bq, new%me_bq, ierr)

 ABI_FREE(dims)
 ABI_FREE(periods)
 ABI_FREE(keepdim)
 call xmpi_comm_free(comm_cart)
#endif

 ! Build table with list of perturbations treated by this CPU.
 ABI_MALLOC(new%my_pinfo, (3, new%my_npert))
 ABI_MALLOC(new%pert_table, (2, natom3))

 do iatom=1,cryst%natom
   do idir=1,3
     pertcase = idir + (iatom-1) * 3
     all_pinfo(:, pertcase) = [idir, iatom, pertcase]
     new%pert_table(1, pertcase) = (pertcase - 1) / (natom3 / new%nprocs_pert)
   end do
 end do
 bstart = (natom3 / new%nprocs_pert) * new%me_pert + 1
 bstop = bstart + new%my_npert - 1
 new%my_pinfo = all_pinfo(:, bstart:bstop)

 new%pert_table(2, :) = -1
 do ii=1,new%my_npert
   ip = new%my_pinfo(3, ii)
   new%pert_table(2, ip) = ii
 end do
 !write(std_out,*)"my_npert", new%my_npert, "nprocs_pert", new%nprocs_pert; write(std_out,*)"my_pinfo", new%my_pinfo

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 ! By default do not skip, if set skip all but specified
 ABI_MALLOC(new%phmodes_skip, (natom3))
 new%phmodes_skip = 0
 if (all(dtset%eph_phrange /= 0)) then
   if (minval(dtset%eph_phrange) < 1 .or. &
       maxval(dtset%eph_phrange) > natom3 .or. &
       dtset%eph_phrange(2) < dtset%eph_phrange(1)) then
     MSG_ERROR('Invalid range for eph_phrange. Should be between [1, 3*natom] and eph_modes(2) > eph_modes(1)')
   end if
   call wrtout(std_out, sjoin(" Including phonon modes between [", &
               itoa(dtset%eph_phrange(1)), ',', itoa(dtset%eph_phrange(2)), "]"))
   new%phmodes_skip = 0
   new%phmodes_skip(dtset%eph_phrange(1):dtset%eph_phrange(2)) = 1
 end if

 if (.not. new%imag_only) then
   ! Split bands among the procs inside comm_bsum.
   call xmpi_split_work(new%nbsum, new%comm_bsum, new%my_bsum_start, new%my_bsum_stop)
   if (new%my_bsum_start == new%nbsum + 1) then
     MSG_ERROR("sigmaph code does not support idle processes! Decrease ncpus or increase nband or use eph_np_pqbks input var.")
   end if
   new%my_bsum_start = new%bsum_start + new%my_bsum_start - 1
   new%my_bsum_stop = new%bsum_start + new%my_bsum_stop - 1
 end if

 call wrtout(std_out, sjoin(" Global bands for self-energy sum, bsum_start: ", itoa(new%bsum_start), &
   " bsum_bstop:", itoa(new%bsum_stop)))
 call wrtout(std_out, sjoin(" Allocating and treating bands from my_bsum_start: ", itoa(new%my_bsum_start), &
   " up to my_bsum_stop:", itoa(new%my_bsum_stop)))

 ! Distribute DFPT potentials (IBZ q-points) inside comm_qpt.
 ! Note that we distribute IBZ instead of the full BZ or the IBZ_k inside the loop over ikcalc.
 ! This means that the load won't be equally distributed but memory will scale with nprocs_qpt.
 ! To reduce load imbalance, we sort qibz points by norm and use cyclic distribution.
 ABI_ICALLOC(new%itreat_qibz, (new%nqibz))
 call sort_rpts(new%nqibz, new%qibz, cryst%gmet, iperm)
 do ii=1,new%nqibz
   iq_ibz = iperm(ii)
   if (mod(ii, new%nprocs_qpt) == new%me_qpt) new%itreat_qibz(iq_ibz) = 1
 end do
 !new%itreat_qibz = 1
 ABI_FREE(iperm)

 call wrtout(std_out, sjoin("P Number of q-points in the IBZ treated by this proc: " ,itoa(count(new%itreat_qibz == 1))))

 ! ================================================================
 ! Allocate arrays used to store final results and set them to zero
 ! ================================================================
 ABI_ICALLOC(new%qp_done, (new%nkcalc, new%nsppol))
 ABI_CALLOC(new%vals_e0ks, (new%ntemp, new%max_nbcalc))
 ABI_CALLOC(new%dvals_de0ks, (new%ntemp, new%max_nbcalc))
 ABI_CALLOC(new%dw_vals, (new%ntemp, new%max_nbcalc))

 ! Frequency dependent stuff
 if (new%nwr > 0) then
   ABI_CALLOC(new%vals_wr, (new%nwr, new%ntemp, new%max_nbcalc))
   ABI_CALLOC(new%wrmesh_b, (new%nwr, new%max_nbcalc))
 end if

 ! Prepare calculation of generalized Eliashberg functions by setting gfw_nomega
 ! prteliash == 0 deactivates computation (default).
 new%gfw_nomega = 0
 if (dtset%prteliash /= 0) then
   new%gfw_nomega = nint((ifc%omega_minmax(2) - ifc%omega_minmax(1) ) / dtset%ph_wstep) + 1
   ABI_MALLOC(new%gfw_mesh, (new%gfw_nomega))
   new%gfw_mesh = arth(ifc%omega_minmax(1), dtset%ph_wstep, new%gfw_nomega)
   ABI_MALLOC(new%gfw_vals, (new%gfw_nomega, 3, new%max_nbcalc))
 end if

 new%a2f_ne = 0
 if (dtset%prteliash == 3) then
   ! TODO: dosdeltae should have a default value.
   ! TODO: Use logmesh/double mesh for electrons?
   estep = dtset%dosdeltae; if (estep <= zero) estep = 0.05 * eV_Ha
   new%a2f_ne = nint((maxval(ebands%eig) - minval(ebands%eig)) / estep) + 1
   if (my_rank == master) then
     write(std_out, *)"Computing a2f with ", new%a2f_ne, " points for electrons and ", new%gfw_nomega, " points for phonons."
     write(std_out, *)"doseltae:", estep, ", tsmear:", dtset%tsmear
   end if
   ABI_MALLOC(new%a2f_emesh, (new%a2f_ne))
   new%a2f_emesh = arth(minval(ebands%eig), estep, new%a2f_ne)
   ABI_CALLOC(new%a2few, (new%a2f_ne, new%gfw_nomega, new%max_nbcalc))
 end if

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

 ! TODO: Should add support for getwfkfine_path
 new%use_doublegrid = .False.
 if (dtset%getwfkfine /= 0 .or. dtset%irdwfkfine /= 0 .or. dtset%getwfkfine_path /= ABI_NOFILE) then

   wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)
   if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
     MSG_ERROR(msg)
   end if
   call wrtout([std_out, ab_out], " EPH double grid interpolation: will read energies from: "//trim(wfk_fname_dense), newlines=1)
   ebands_dense = wfk_read_ebands(wfk_fname_dense, comm)

   !TODO add consistency check: number of bands and kpoints (commensurability)
   ABI_CHECK(ebands_dense%mband == ebands%mband, 'Inconsistent number of bands for the fine and dense grid')
   new%use_doublegrid = .True.

 else if (any(dtset%bs_interp_kmult /= 0)) then
   ! Read bs_interpmult
   call wrtout([std_out, ab_out]," EPH interpolation: will use star functions interpolation.", newlines=1)
   ! Interpolate band energies with star-functions
   params = 0; params(1) = 1; params(2) = 5
   if (nint(dtset%einterp(1)) == 1) params = dtset%einterp
   write(std_out, "(a, 4(f5.2, 2x))")"SKW parameters for double-grid:", params

   !TODO: mband should be min of nband
   band_block = [1, ebands%mband]
   ! TODO: Now we should use this band range.
   ! Note that we start from 1 because we are gonna use ebands_dense to compute the Fermi level.
   !band_block = [1, new%bsum_stop]
   intp_kptrlatt(:,1) = [ebands%kptrlatt(1,1)*dtset%bs_interp_kmult(1), 0, 0]
   intp_kptrlatt(:,2) = [0, ebands%kptrlatt(2,2)*dtset%bs_interp_kmult(2), 0]
   intp_kptrlatt(:,3) = [0, 0, ebands%kptrlatt(3,3)*dtset%bs_interp_kmult(3)]

   intp_nshiftk = 1
   intp_shiftk = zero
   ebands_dense = ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt, &
                                      intp_nshiftk, intp_shiftk, band_block, comm)
   new%use_doublegrid = .True.
 end if

 ! TODO: Should recheck the case bstart > 1 as I got weird results.
 ! bstart and new%bsum select the band range.
 ! Here we compute the weights only for the states included in the sum
 bstart = new%bsum_start
 new%frohl_model = nint(dtset%frohl_params(1))
 if (new%qint_method > 0) then
   if (new%use_doublegrid) then
     ! Double-grid technique from ab-initio energies or star-function interpolation.
     new%ephwg = ephwg_from_ebands(cryst, ifc, ebands_dense, bstart, new%nbsum, new%frohl_model, comm)
     new%eph_doublegrid = eph_double_grid_new(cryst, ebands_dense, ebands%kptrlatt, ebands_dense%kptrlatt)
   else
     downsample = any(ebands%kptrlatt /= qptrlatt) .or. ebands%nshiftk /= my_nshiftq
     if (ebands%nshiftk == my_nshiftq) downsample = downsample .or. any(ebands%shiftk /= my_shiftq)
     if (downsample) then
       MSG_COMMENT("K-mesh != Q-mesh for self-energy. Will downsample electron energies.")
       tmp_ebands = ebands_downsample(ebands, cryst, qptrlatt, my_nshiftq, my_shiftq)
       new%ephwg = ephwg_from_ebands(cryst, ifc, tmp_ebands, bstart, new%nbsum, new%frohl_model, comm)
       call ebands_free(tmp_ebands)
     else
       new%ephwg = ephwg_from_ebands(cryst, ifc, ebands, bstart, new%nbsum, new%frohl_model, comm)
     end if
   end if
 else
   if (new%use_doublegrid) then
     new%eph_doublegrid = eph_double_grid_new(cryst, ebands_dense, ebands%kptrlatt, ebands_dense%kptrlatt)
     new%ephwg = ephwg_from_ebands(cryst, ifc, ebands_dense, bstart, new%nbsum, new%frohl_model, comm)
   endif
 end if

 call cwtime_report(" sigmaph_new: doublegrid", cpu, wall, gflops)

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
     if (mod(it, nprocs) /= my_rank) cycle ! MPI parallelism.
     ! Use Fermi-Dirac occopt
     call ebands_set_scheme(tmp_ebands, occopt3, new%kTmesh(it), dtset%spinmagntarget, dtset%prtvol)
     call ebands_set_nelect(tmp_ebands, ebands%nelect, dtset%spinmagntarget, msg)
     new%mu_e(it) = tmp_ebands%fermie
     !
     ! Check that the total number of electrons is correct
     ! This is to trigger problems as the routines that calculate the occupations in ebands_set_nelect
     ! are different from the occ_fd that will be used in the rest of the code.
     nelect = ebands_calc_nelect(tmp_ebands, new%kTmesh(it), new%mu_e(it))

     if (abs(nelect - ebands%nelect) > tol6) then
       ! For T = 0 the number of occupied states goes in discrete steps (according to the k-point sampling)
       ! for finite doping its hard to find nelect that exactly matches ebands%nelect.
       ! in this case we print a warning
       write(msg,'(3(a,f10.6))')&
         'Calculated number of electrons nelect = ',nelect, &
         ' does not correspond with ebands%nelect = ',tmp_ebands%nelect,' for T = ',new%kTmesh(it)
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
 if (new%frohl_model /= 0) then
   ! Init parameters for numerical integration inside sphere.
   ! Set sphere radius to a fraction of the smallest reciprocal lattice vector.
   ! Use sphere whose radius is a fraction of the smallest rec lattice vector.
   new%qrad = huge(one)
   do ii=1,3
     new%qrad = min(new%qrad, sqrt(dot_product(cryst%gprimd(:, ii), cryst%gprimd(:, ii))))
   end do
   ! Set angular mesh.
   new%ntheta = int(dtset%frohl_params(2))
   if (new%ntheta <= 0) new%ntheta = 10
   new%nphi = 2 * new%ntheta
   new%nqr = 10  !int(dtset%frohl_params(3))
   new%qrad = new%qrad !/ 2.0_dp  dtset%frohl_params(2)
   new%qdamp = new%qrad !/ four
   write(std_out,"(a)")" Activating computation of Frohlich self-energy:"
   write(std_out,"(2(a,i0,1x))")" ntheta: ", new%ntheta, "nphi: ", new%nphi
   write(std_out,"((a,i0,1x,a,f6.3,1x,a))")" nqr points: ", new%nqr, "qrad: ", new%qrad, " [Bohr^-1]"

   ! Initialize angular mesh qvers_cart and angwgth (inspired to initang in m_pawang)
   ABI_MALLOC(th, (new%ntheta))
   ABI_MALLOC(wth, (new%ntheta))
   con = two_pi / new%nphi
   call gauleg(-one, one, th, wth, new%ntheta)

   ! Initialize qvers_cart and angular weights angwgth
   ! NB: summing over f * angwgth gives the spherical average 1/(4pi) \int domega f(omega)
   new%angl_size = new%ntheta * new%nphi
   ABI_MALLOC(new%qvers_cart, (3, new%angl_size))
   ABI_MALLOC(new%angwgth, (new%angl_size))
   npoints = 0
   do it = 1, new%ntheta
     cos_theta = th(it)
     sin_theta = sqrt(one - cos_theta*cos_theta)
     do ip = 1, new%nphi
       ang = con * (ip-1)
       cos_phi = cos(ang); sin_phi = sin(ang)
       npoints = npoints + 1
       new%qvers_cart(1, npoints) = sin_theta * cos_phi
       new%qvers_cart(2, npoints) = sin_theta * sin_phi
       new%qvers_cart(3, npoints) = cos_theta
       ! Normalization required
       new%angwgth(npoints) = wth(it) / (two * new%nphi)
     end do
   end do
   !write(std_out, *)"Sum angwgth: ", sum(new%angwgth)
   !write(std_out, *)"int sig(theta): ", sum(new%angwgth * sqrt(one - new%qvers_cart(3, :) **2))
   !write(std_out, *)pi/4
   !stop

   ABI_FREE(th)
   ABI_FREE(wth)

   ! Build SKW object for all bands where QP corrections are wanted.
   ! TODO: check band_block because I got weird results (don't remember if with AbiPy or Abinit)
   !skw_cplex = 1; if (kpts_timrev_from_kptopt(ebands%kptopt) == 0) skw_cplex = 2
   !skw_band_block = [minval(new%bstart_ks), maxval(new%bstart_ks + new%nbcalc_ks - 1)]
   !params = 0; params(1) = 5
   !if (nint(dtset%einterp(1)) == 1) params(1:3) = dtset%einterp(2:4)
   !write(std_out, "(a, 4(f5.2, 2x))")"SKW parameters used to interpolate e_{nk+q} in Frohlich self-energy:", params
   !new%frohl_skw = skw_new(cryst, params(2:), skw_cplex, ebands%mband, ebands%nkpt, ebands%nsppol, &
   !  ebands%kptns, ebands%eig, skw_band_block, comm)
   !call new%frohl_skw%print(std_out)

   ! Allocate arrays to store Frohlich results for given (ikcalc, spin)
   ABI_CALLOC(new%frohl_vals_e0ks, (new%ntemp, new%max_nbcalc))
   ABI_CALLOC(new%frohl_dvals_de0ks, (new%ntemp, new%max_nbcalc))
   if (new%nwr > 0) then
     ABI_CALLOC(new%frohl_vals_wr, (new%nwr, new%ntemp, new%max_nbcalc))
   end if
 end if

 if (new%calc_mrta) then
   ABI_CALLOC(new%linewidth_mrta, (new%ntemp, new%max_nbcalc))
 end if

 call cwtime_report(" sigmaph_new: all", cpu_all, wall_all, gflops_all)

end function sigmaph_new
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_write
!! NAME
!!  sigmaph_write
!!
!! FUNCTION
!!  Define dimensions and netcdf arrays in SIGEPH file.
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  wfk_hdr=Header of the WFK file.
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  dtfil<datafiles_type>=variables related to files.
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_write(self, dtset, cryst, ebands, wfk_hdr, dtfil, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 class(sigmaph_t),intent(inout) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(hdr_type),intent(in) :: wfk_hdr
 type(datafiles_type),intent(in) :: dtfil

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: my_rank, ii, edos_intmeth
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
 !character(len=500) :: msg
 real(dp) :: edos_broad, edos_step,  cpu_all, wall_all, gflops_all, cpu, wall, gflops
 character(len=fnlen) :: path
 type(edos_t) :: edos

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 if (dtset%prtdos /= 0) then
   call cwtime(cpu, wall, gflops, "start")
   ! Compute electron DOS.
   edos_intmeth = 2; if (self%bcorr == 1) edos_intmeth = 3
   if (dtset%prtdos == 1) edos_intmeth = 1
   !if (dtset%prtdos == -2) edos_intmeth = 3
   edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
   call wrtout(std_out, " Computing electron dos. Use prtdos 0 to disable this part...", do_flush=.True.)
   edos = ebands_get_edos(ebands, cryst, edos_intmeth, edos_step, edos_broad, comm)
   if (my_rank == master) then
     path = strcat(dtfil%filnam_ds(4), "_EDOS")
     call wrtout(ab_out, sjoin("- Writing electron DOS to file:", path))
     call edos%write(path)
     call edos%print(unit=std_out)
     !call edos%print(unit=ab_out)
   end if
   call cwtime_report(" sigmaph_new: ebands", cpu, wall, gflops)
 end if

 ! Open netcdf file (only master works for the time being because I cannot assume HDF5 + MPI-IO)
 ! This could create problems if MPI parallelism over (spin, nkptgw) ...
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   path = strcat(dtfil%filnam_ds(4), "_SIGEPH.nc")

   ! Master creates the netcdf file used to store the results of the calculation.
   NCF_CHECK(nctk_open_create(self%ncid, path, xmpi_comm_self))
   ncid = self%ncid

   NCF_CHECK(hdr_ncwrite(wfk_hdr, ncid, fform_from_ext("SIGEPH.nc"), nc_define=.True.))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))
   if (dtset%prtdos /= 0) then
     NCF_CHECK(edos%ncwrite(ncid))
   end if

   ! Add sigma_eph dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nkcalc", self%nkcalc), nctkdim_t("max_nbcalc", self%max_nbcalc), &
     nctkdim_t("nsppol", self%nsppol), nctkdim_t("ntemp", self%ntemp), nctkdim_t("natom3", 3 * cryst%natom), &
     nctkdim_t("nqibz", self%nqibz), nctkdim_t("nqbz", self%nqbz)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   if (self%nwr > 0) then
     NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nwr", self%nwr)]))
   end if
   if (dtset%prteliash /= 0) then
     NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("gfw_nomega", self%gfw_nomega)]))
     if (dtset%prteliash == 3) then
       NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("a2f_ne", self%a2f_ne)]))
     end if
   end if

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", "symdynmat", &
     "ph_intmeth", "eph_intmeth", "qint_method", "eph_transport", &
     "imag_only", "frohl_model", "symv1scf", "dvdb_add_lr"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"])
   NCF_CHECK(ncerr)

   ! Define arrays for results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngqpt", "int", "three"), &
     nctkarr_t("eph_ngqpt_fine", "int", "three"), &
     nctkarr_t("eph_phrange", "int", "two"), &
     nctkarr_t("ddb_ngqpt", "int", "three"), &
     nctkarr_t("dvdb_ngqpt", "int", "three"), &
     nctkarr_t("ph_ngqpt", "int", "three"), &
     nctkarr_t("sigma_ngkpt", "int", "three"), &
     nctkarr_t("sigma_erange", "dp", "two"), &
     nctkarr_t("frohl_params", "dp", "four"), &
     nctkarr_t("bstart_ks", "int", "nkcalc, nsppol"), &
     !nctkarr_t("bstop_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("nbcalc_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("kcalc", "dp", "three, nkcalc"), &
     nctkarr_t("kcalc2ibz", "int", "nkcalc, six"), &
     nctkarr_t("kTmesh", "dp", "ntemp"), &
     nctkarr_t("mu_e", "dp", "ntemp"), &
     nctkarr_t("qp_done", "int", "nkcalc, nsppol"), &
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

   if (self%calc_mrta) then
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("linewidth_mrta", "dp", "ntemp, max_nbcalc, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
   end if

   if (self%frohl_model == 1) then
     if (self%imag_only) then
       ncerr = nctk_def_arrays(ncid, [ &
         nctkarr_t("frohl_deltas_sphcorr", "dp", "two, ntemp, max_nbcalc, natom3, nkcalc, nsppol") &
       ])
       NCF_CHECK(ncerr)
     end if
   end if

   !elif (self%frohl_model == 2) then
   !  ! Arrays storing the Frohlich self-energy.
   !  ! These arrays get two extra dimensions on file (nkcalc, nsppol).
   !  ncerr = nctk_def_arrays(ncid, [ &
   !    nctkarr_t("frohl_vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
   !    nctkarr_t("frohl_dvals_de0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol") ])
   !  NCF_CHECK(ncerr)
   !  if (self%nwr > 0) then
   !    ncerr = nctk_def_arrays(ncid, [ &
   !      nctkarr_t("frohl_vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
   !      nctkarr_t("frohl_spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol") &
   !    ])
   !    NCF_CHECK(ncerr)
   !  end if
   !end if

   if (self%nwr > 0) then
     ! Make room for the spectral function. These arrays get two extra dimensions on file (nkcalc, nsppol).
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("wrmesh_b", "dp", "nwr, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("spfunc_wr", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
   end if

   if (dtset%prteliash /= 0) then
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("gfw_mesh", "dp", "gfw_nomega"), &
       nctkarr_t("gfw_vals", "dp", "gfw_nomega, three, max_nbcalc, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
     if (dtset%prteliash == 3) then
       ncerr = nctk_def_arrays(ncid, [ &
         nctkarr_t("a2f_emesh", "dp", "a2f_ne"), &
         nctkarr_t("a2few", "dp", "a2f_ne, gfw_nomega, max_nbcalc, nkcalc, nsppol") &
       ])
       NCF_CHECK(ncerr)
     end if
   end if

   ! ======================================================
   ! Write data that do not depend on the (kpt, spin) loop.
   ! ======================================================
   NCF_CHECK(nctk_set_datamode(ncid))
   ii = 0; if (self%imag_only) ii = 1
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", &
     "symdynmat", "ph_intmeth", "eph_intmeth", "qint_method", &
     "eph_transport", "imag_only", "frohl_model", "symv1scf", "dvdb_add_lr"], &
     [dtset%eph_task, self%symsigma, self%nbsum, self%bsum_start, self%bsum_stop, &
     dtset%symdynmat, dtset%ph_intmeth, dtset%eph_intmeth, self%qint_method, &
     dtset%eph_transport, ii, self%frohl_model, dtset%symv1scf, dtset%dvdb_add_lr])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"], &
     [aimag(self%ieta), self%wr_step, dtset%eph_fsewin, dtset%eph_fsmear, dtset%eph_extrael, dtset%eph_fermie, &
     dtset%ph_wstep, dtset%ph_smear])
   NCF_CHECK(ncerr)

   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt"), self%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_ngqpt_fine"), dtset%eph_ngqpt_fine))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ddb_ngqpt"), dtset%ddb_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "dvdb_ngqpt"), dtset%dvdb_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ph_ngqpt"), dtset%ph_ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma_ngkpt"), dtset%sigma_ngkpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma_erange"), dtset%sigma_erange))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "frohl_params"), dtset%frohl_params))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_phrange"), dtset%eph_phrange))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "bstart_ks"), self%bstart_ks))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nbcalc_ks"), self%nbcalc_ks))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc"), self%kcalc))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ibz"), self%kcalc2ibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), self%kTmesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mu_e"), self%mu_e))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eta"), aimag(self%ieta)))
   if (dtset%prteliash /= 0) then
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gfw_mesh"), self%gfw_mesh))
     if (dtset%prteliash == 3) then
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "a2f_emesh"), self%a2f_emesh))
     end if
   end if
   !NCF_CHECK(nf90_sync(ncid))
 end if ! master
#endif

 call edos%free()
 call cwtime_report(" sigmaph_new: netcdf", cpu_all, wall_all, gflops_all)

end subroutine sigmaph_write
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_read
!! NAME
!!  sigmaph_read
!!
!! FUNCTION
!!  Start a sigmaph instance from a netcdf file.
!!  This routine serves only to read some basic dimensions of sigmaph_t to verify if a restart is possible
!!
!! INPUTS
!!  path= SIGEPH Filename.
!!  dtset<dataset_type>=All input variables for this dataset.
!!  ecut=Cutoff energy for wavefunctions.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(sigmaph_t) function sigmaph_read(path, dtset, comm, msg, ierr, keep_open, extrael_fermie) result(new)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 integer,intent(out) :: ierr
 type(dataset_type),intent(in) :: dtset
 character(len=500),intent(out) :: msg
 real(dp), optional, intent(inout) :: extrael_fermie(2)
 logical,optional,intent(in) :: keep_open

!Local variables ------------------------------
!scalars
 integer :: imag_only
#ifdef HAVE_NETCDF
 integer :: ncid,varid,ncerr
#endif
 real(dp) :: eph_fermie, eph_fsewin, ph_wstep, ph_smear, eta, eph_extrael, eph_fsmear
 character(len=fnlen) :: path
!arrays
 integer :: eph_task, symdynmat, ph_intmeth, eph_intmeth, eph_transport
 integer :: eph_ngqpt_fine(3), ddb_ngqpt(3), ph_ngqpt(3), sigma_ngkpt(3), frohl_params(4)
 real(dp) :: sigma_erange(2)

! *************************************************************************

 ! Open netcdf file
 msg = "No message"
#ifdef HAVE_NETCDF
 ierr = 0

 if (.not. file_exists(path)) then
   ierr = 1; return
 end if
 NCF_CHECK(nctk_open_read(ncid, path, comm))

 !TODO?
 !NCF_CHECK(cryst%ncread(ncid))
 !NCF_CHECK(ebands_ncread(ebands, ncid))

 ! Read sigma_eph dimensions.
 NCF_CHECK(nctk_get_dim(ncid, "nkcalc", new%nkcalc))
 NCF_CHECK(nctk_get_dim(ncid, "max_nbcalc", new%max_nbcalc))
 NCF_CHECK(nctk_get_dim(ncid, "nsppol", new%nsppol))
 NCF_CHECK(nctk_get_dim(ncid, "ntemp", new%ntemp))
 !NCF_CHECK(nctk_get_dim(ncid, "natom3", natom3))
 NCF_CHECK(nctk_get_dim(ncid, "nqibz", new%nqibz))
 NCF_CHECK(nctk_get_dim(ncid, "nqbz", new%nqbz))

 !NCF_CHECK(nctk_get_dim(ncid,"nwr", new%nwr))
 !NCF_CHECK(nctk_get_dim(ncid,"gfw_nomega", new%gfw_nomega))

 ! ======================================================
 ! Read data that does not depend on the (kpt, spin) loop.
 ! ======================================================
 NCF_CHECK(nf90_get_var(ncid, vid("symsigma"), new%symsigma))
 NCF_CHECK(nf90_get_var(ncid, vid("nbsum"), new%nbsum))
 NCF_CHECK(nf90_get_var(ncid, vid("bsum_start"), new%bsum_start))
 NCF_CHECK(nf90_get_var(ncid, vid("bsum_stop"), new%bsum_stop))

 NCF_CHECK(nf90_get_var(ncid, vid("qint_method"), new%qint_method))
 NCF_CHECK(nf90_get_var(ncid, vid("frohl_model"), new%frohl_model))
 NCF_CHECK(nf90_get_var(ncid, vid("imag_only"), imag_only))
 new%imag_only = (imag_only == 1)

 ABI_MALLOC(new%kcalc, (3, new%nkcalc))
 ABI_MALLOC(new%bstart_ks, (new%nkcalc, new%nsppol))
 ABI_MALLOC(new%bstop_ks, (new%nkcalc, new%nsppol))
 ABI_MALLOC(new%nbcalc_ks, (new%nkcalc, new%nsppol))
 ABI_MALLOC(new%mu_e, (new%ntemp))
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 ABI_MALLOC(new%kcalc2ibz, (new%nkcalc, 6))

 NCF_CHECK(nf90_get_var(ncid, vid("ngqpt"), new%ngqpt))
 NCF_CHECK(nf90_get_var(ncid, vid("bstart_ks"), new%bstart_ks))
 NCF_CHECK(nf90_get_var(ncid, vid("nbcalc_ks"), new%nbcalc_ks))
 new%bstop_ks = new%bstart_ks + new%nbcalc_ks - 1
 NCF_CHECK(nf90_get_var(ncid, vid("kcalc"), new%kcalc))
 NCF_CHECK(nf90_get_var(ncid, vid("kcalc2ibz"), new%kcalc2ibz))
 NCF_CHECK(nf90_get_var(ncid, vid("kTmesh"), new%kTmesh))
 NCF_CHECK(nf90_get_var(ncid, vid("wr_step"), new%wr_step))
 NCF_CHECK(nf90_get_var(ncid, vid("mu_e"), new%mu_e))
 NCF_CHECK(nf90_get_var(ncid, vid("eta"), eta))
 new%ieta = j_dpc * eta

 ! Try to read the done array used to implement restart capabilities.
 ABI_ICALLOC(new%qp_done, (new%nkcalc, new%nsppol))
 ncerr = nf90_inq_varid(ncid, "qp_done", varid)
 if (ncerr == nf90_noerr) then
   NCF_CHECK(nf90_get_var(ncid, varid, new%qp_done))
 end if

 ! ============================================================
 ! Read and check consistency against dtset
 ! ============================================================
 NCF_CHECK(nf90_get_var(ncid, vid("eph_fsewin"), eph_fsewin))
 NCF_CHECK(nf90_get_var(ncid, vid("eph_fsmear"), eph_fsmear))
 NCF_CHECK(nf90_get_var(ncid, vid("eph_extrael"), eph_extrael))
 NCF_CHECK(nf90_get_var(ncid, vid("eph_fermie"), eph_fermie))
 NCF_CHECK(nf90_get_var(ncid, vid("ph_wstep"), ph_wstep))
 NCF_CHECK(nf90_get_var(ncid, vid("ph_smear"), ph_smear))
 ABI_CHECK(eph_fsewin  == dtset%eph_fsewin,  "netcdf eph_fsewin != input file")
 ABI_CHECK(eph_fsmear  == dtset%eph_fsmear,  "netcdf eph_fsmear != input file")
 ABI_CHECK(ph_wstep    == dtset%ph_wstep,    "netcdf ph_wstep != input file")
 ABI_CHECK(ph_smear    == dtset%ph_smear,    "netcdf ph_smear != input file")
 if (present(extrael_fermie)) then
   extrael_fermie = [eph_extrael, eph_fermie]
 else
   ABI_CHECK(eph_extrael == dtset%eph_extrael, "netcdf eph_extrael != input file")
   ABI_CHECK(eph_fermie  == dtset%eph_fermie,  "netcdf eph_feremie != input file")
 end if

 NCF_CHECK(nf90_get_var(ncid, vid("eph_task"), eph_task))
 NCF_CHECK(nf90_get_var(ncid, vid("symdynmat"), symdynmat))
 NCF_CHECK(nf90_get_var(ncid, vid("ph_intmeth"), ph_intmeth))
 NCF_CHECK(nf90_get_var(ncid, vid("eph_intmeth"), eph_intmeth))
 NCF_CHECK(nf90_get_var(ncid, vid("eph_transport"), eph_transport))
 ABI_CHECK(symdynmat     == dtset%symdynmat,    "netcdf symdynmat != input file")
 ABI_CHECK(ph_intmeth    == dtset%ph_intmeth,   "netcdf ph_intmeth != input file")
 ABI_CHECK(eph_intmeth   == dtset%eph_intmeth,  "netcdf eph_intmeth != input file")
 ABI_CHECK(eph_transport == dtset%eph_transport,"netcdf eph_transport != input file")

 NCF_CHECK(nf90_get_var(ncid, vid("eph_ngqpt_fine"), eph_ngqpt_fine))
 NCF_CHECK(nf90_get_var(ncid, vid("ddb_ngqpt"), ddb_ngqpt))
 NCF_CHECK(nf90_get_var(ncid, vid("ph_ngqpt"), ph_ngqpt))
 NCF_CHECK(nf90_get_var(ncid, vid("sigma_ngkpt"), sigma_ngkpt))
 NCF_CHECK(nf90_get_var(ncid, vid("frohl_params"), frohl_params))
 NCF_CHECK(nf90_get_var(ncid, vid("sigma_erange"), sigma_erange))

 if (present(keep_open)) then
   new%ncid = ncid
 else
   NCF_CHECK(nf90_close(ncid))
   ! so that the structure is properly freed
   new%ncid = nctk_noid
 end if

 ABI_CHECK(all(dtset%eph_ngqpt_fine == eph_ngqpt_fine),"netcdf eph_ngqpt_fine != input file")
 ABI_CHECK(all(dtset%ddb_ngqpt      == ddb_ngqpt),     "netcdf ddb_ngqpt != input file")
 ABI_CHECK(all(dtset%ph_ngqpt       == ph_ngqpt),      "netcdf ph_ngqpt != input file")
 ABI_CHECK(all(dtset%sigma_ngkpt    == sigma_ngkpt),   "netcdf sigma_ngkpt != input file")
 !ABI_CHECK(all(abs(dtset%frohl_params - frohl_params) < tol6), "netcdf frohl_params != input file")
 !ABI_CHECK(all(abs(dtset%sigma_erange - sigma_erange) < tol6),  "netcdf sigma_erange != input file")

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
end function vid

#endif

end function sigmaph_read
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_get_ebands
!! NAME
!!  sigmaph_get_ebands
!!
!! FUNCTION
!!  Read quantities from the sigmaph to a ebands structure and return mapping
!!
!! INPUTS
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  opt=integer option selecting what to read on the ebands object. 1-only mapping, 10+n-read n temperature linewidths
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function sigmaph_get_ebands(self, cryst, ebands, linewidth_serta, linewidth_mrta, &
                                       velocity, comm, ierr, indq2ebands) result(new)

!Arguments -----------------------------------------------
 integer,intent(in) :: comm
 class(sigmaph_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 integer, intent(out) :: ierr
 real(dp), allocatable, intent(out) :: linewidth_serta(:,:,:,:), linewidth_mrta(:,:,:,:), velocity(:,:,:,:)
 integer, allocatable, optional, intent(out) :: indq2ebands(:)
 character(len=500) :: msg

!Local variables -----------------------------------------
 integer,parameter :: sppoldbl1 = 1, timrev1 = 1
 integer :: spin, ikpt, ikcalc, iband, itemp, nkcalc, nsppol, nkpt
 integer :: band_ks, bstart_ks, nbcalc_ks, mband
 logical :: has_mrta, has_vel, has_car_vel, has_red_vel
#ifdef HAVE_NETCDF
 integer :: ncerr, varid
#endif
 !character(len=fnlen) :: path
 integer,allocatable :: indkk(:,:)
 real(dp) :: dksqmax, vk_red(2,3), vk_car(2,3)

! *************************************************************************

 ! copy useful dimensions
 nsppol = self%nsppol
 nkcalc = self%nkcalc
 nkpt = ebands%nkpt

 ! Map ebands kpoints to sigmaph
 ABI_MALLOC(indkk,(nkcalc, 6))
 call listkk(dksqmax, cryst%gmet, indkk, ebands%kptns, self%kcalc, ebands%nkpt, nkcalc, cryst%nsym, &
             sppoldbl1, cryst%symafm, cryst%symrec, timrev1, comm, exit_loop=.True., use_symrec=.True.)

 if (dksqmax > tol12) then
    write(msg, '(3a,es16.6,a)' ) &
     "Error mapping ebands to sigmaph",ch10,&
     "the k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10
    MSG_ERROR(msg)
 end if

 ! store mapping to return
 if (present(indq2ebands)) then
   ABI_MALLOC(indq2ebands,(nkcalc))
   indq2ebands(:) = indkk(:,1)
 end if

 ! Allocate using only the relevant bands for transport
 ! includin valence states to allow to compute different doping
 mband = maxval(self%bstop_ks)
 new = ebands_chop(ebands, 1, mband)
 !mband = ebands%mband
 !call ebands_copy(ebands,new)

#ifdef HAVE_NETCDF
 has_mrta = .true.
 ierr = nf90_inq_varid(self%ncid, "vals_e0k", varid)
 if (ierr /= nf90_noerr) has_mrta = .false.

 has_car_vel = .false.
 has_red_vel = .false.
 ierr = nf90_inq_varid(self%ncid, "vcar_calc", varid)
 if (ierr == nf90_noerr) has_car_vel = .true.
 ierr = nf90_inq_varid(self%ncid, "vred_calc", varid)
 if (ierr == nf90_noerr) has_red_vel = .true.
 has_vel = (has_car_vel .or. has_red_vel)

 ! Read linewidths from sigmaph
 ABI_CALLOC(linewidth_serta,(self%ntemp,mband,nkpt,nsppol))
 if (has_mrta) then
   ABI_CALLOC(linewidth_mrta,(self%ntemp,mband,nkpt,nsppol))
 end if
 if (has_vel) then
   ABI_CALLOC(velocity,(3,mband,nkpt,nsppol))
 end if

 do spin=1,nsppol
   do ikcalc=1,nkcalc
     bstart_ks = self%bstart_ks(ikcalc,spin)
     nbcalc_ks = self%nbcalc_ks(ikcalc,spin)
     do iband=1,nbcalc_ks
       band_ks = iband + bstart_ks - 1
       ikpt = indkk(ikcalc,1)

       do itemp=1,self%ntemp
         ! Read SERTA lifetimes
         ncerr = nf90_get_var(self%ncid, nctk_idname(self%ncid, "vals_e0ks"),&
                 linewidth_serta(itemp,band_ks,ikpt,spin), start=[2,itemp,iband,ikcalc,spin])
         NCF_CHECK(ncerr)
         ! Read MRTA lifetimes
         if (has_mrta) then
            ncerr = nf90_get_var(self%ncid, nctk_idname(self%ncid, "linewidth_mrta"),&
                                 linewidth_mrta(itemp,band_ks,ikpt,spin), start=[itemp,iband,ikcalc,spin])
            NCF_CHECK(ncerr)
         end if
       end do

       ! Read band velocities
       if (has_vel) then
         if (has_car_vel) then
           ncerr = nf90_get_var(self%ncid, nctk_idname(self%ncid, "vcar_calc"),&
                                velocity(:,band_ks,ikpt,spin), start=[1,iband,ikcalc,spin])
           NCF_CHECK(ncerr)
         end if
         if (has_red_vel) then
           ncerr = nf90_get_var(self%ncid, nctk_idname(self%ncid, "vred_calc"),&
                                vk_red(1,:), start=[1,iband,ikcalc,spin])
           NCF_CHECK(ncerr)
           call ddk_red2car(cryst%rprimd,vk_red,vk_car)
           velocity(:,band_ks,ikpt,spin) = vk_car(1,:)
         end if
       end if

     end do
   end do
 end do
#endif

 ABI_FREE(indkk)

end function sigmaph_get_ebands
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_comp
!! NAME
!!  sigmaph_comp
!!
!! FUNCTION
!!  Compare the headers of two sigmaph_t instances
!!
!! INPUTS
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_comp(self,comp)

!Arguments ------------------------------------
 class(sigmaph_t),intent(in) :: self, comp

! *************************************************************************

 ABI_CHECK(self%nkcalc == comp%nkcalc, "Difference found in nkcalc.")
 ABI_CHECK(self%max_nbcalc == comp%max_nbcalc, "Difference found in max_nbcalc.")
 ABI_CHECK(self%nsppol == self%nsppol, "Difference found in nsppol.")
 ABI_CHECK(self%ntemp == self%ntemp, "Difference found in ntemp.")
 !ABI_CHECK(natom3 == natom3, "")
 ABI_CHECK(self%nqibz == self%nqibz, "Difference found in nqibz.")
 ABI_CHECK(self%nqbz == self%nqbz, "Difference found in nqbz.")

 ! ======================================================
 ! Read data that does not depend on the (kpt, spin) loop.
 ! ======================================================
 ABI_CHECK(self%symsigma == comp%symsigma, "Different value found for symsigma.")
 ABI_CHECK(self%nbsum == comp%nbsum, "Different value found for nbsum.")
 ABI_CHECK(self%bsum_start == comp%bsum_start, "Different value found for bsum_start.")
 ABI_CHECK(self%bsum_stop == comp%bsum_stop, "Different value found for bsum_stop.")
 ABI_CHECK(self%qint_method == comp%qint_method, "Different value found for qint_method.")
 ABI_CHECK(self%frohl_model == comp%frohl_model, "Different value found for frohl_model.")
 ABI_CHECK(self%imag_only .eqv. comp%imag_only, "Difference found in imag_only")
 ABI_CHECK(self%wr_step==comp%wr_step, "Different value found for wr_step")
 ABI_CHECK(self%ieta == comp%ieta, "Different value found for zcut.")

 ABI_CHECK(all(self%ngqpt == comp%ngqpt), "Different value found for ngqpt")
 ABI_CHECK(all(self%bstart_ks == comp%bstart_ks), "Different value found for bstart_ks")
 ABI_CHECK(all(self%nbcalc_ks == comp%nbcalc_ks), "Different value found for bstop_ks")
 ABI_CHECK(all(self%kcalc == comp%kcalc), "Different value found for kcalc")
 ABI_CHECK(all(self%kcalc2ibz == comp%kcalc2ibz), "Different value found for kcalc2ibz")
 ABI_CHECK(all(self%kTmesh == comp%kTmesh), "Different value found for kTmesh")
 ABI_CHECK(all(self%mu_e == comp%mu_e), "Different value found for mu_e")

end subroutine sigmaph_comp
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

!Arguments ------------------------------------
 class(sigmaph_t),intent(inout) :: self

!Local variables-------------------------------
 integer :: ii,jj,ideg

! *************************************************************************

 ! integer
 ABI_SFREE(self%bstart_ks)
 ABI_SFREE(self%bstop_ks)
 ABI_SFREE(self%nbcalc_ks)
 ABI_SFREE(self%kcalc2ibz)
 ABI_SFREE(self%myq2ibz_k)
 ABI_SFREE(self%itreat_qibz)
 ABI_SFREE(self%my_pinfo)
 ABI_SFREE(self%pert_table)
 ABI_SFREE(self%phmodes_skip)
 ABI_SFREE(self%ind_bz2ibz)
 ABI_SFREE(self%indkk_kq)
 ABI_SFREE(self%indq2dvdb_k)
 ABI_SFREE(self%ind_ibzk2ibz)
 ABI_SFREE(self%ibz2dvdb)

 ! real
 ABI_SFREE(self%kcalc)
 ABI_SFREE(self%kTmesh)
 ABI_SFREE(self%mu_e)
 ABI_SFREE(self%e0vals)
 ABI_SFREE(self%vcar_calc)
 ABI_SFREE(self%linewidth_mrta)
 ABI_SFREE(self%cweights)
 ABI_SFREE(self%deltaw_pm)
 ABI_SFREE(self%wrmesh_b)
 ABI_SFREE(self%qvers_cart)
 ABI_SFREE(self%angwgth)
 ABI_SFREE(self%frohl_deltas_sphcorr)
 ABI_SFREE(self%qp_done)
 ABI_SFREE(self%qbz)
 ABI_SFREE(self%qibz)
 ABI_SFREE(self%wtq)
 ABI_SFREE(self%qibz_k)
 ABI_SFREE(self%wtq_k)
 ABI_SFREE(self%gfw_mesh)
 ABI_SFREE(self%gf_nnuq)

 ! complex
 ABI_SFREE(self%vals_e0ks)
 ABI_SFREE(self%frohl_vals_e0ks)
 ABI_SFREE(self%dvals_de0ks)
 ABI_SFREE(self%frohl_dvals_de0ks)
 ABI_SFREE(self%dw_vals)
 ABI_SFREE(self%vals_wr)
 ABI_SFREE(self%frohl_vals_wr)
 ABI_SFREE(self%gfw_vals)
 ABI_SFREE(self%a2f_emesh)
 ABI_SFREE(self%a2few)

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

 call self%ephwg%free()
 call self%eph_doublegrid%free()
 call self%frohl_skw%free()

 ! Deallocate MPI communicators
 call xmpi_comm_free(self%comm_pert)
 call xmpi_comm_free(self%comm_qpt)
 call xmpi_comm_free(self%comm_bsum)
 call xmpi_comm_free(self%comm_bq)

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
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t> = Crystal structure.
!!  dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  ikcalc=Index of the k-point to compute.
!!  prtvol= Verbosity level
!!  comm= MPI communicator
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_setup_kcalc(self, dtset, cryst, dvdb, ebands, ikcalc, prtvol, comm)

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, prtvol, comm
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 class(sigmaph_t),target,intent(inout) :: self
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(in) :: dvdb

!Local variables-------------------------------
 integer,parameter :: sppoldbl1 = 1, timrev1 = 1, master = 0
 integer :: spin, my_rank, iq_ibz, nprocs !, nbcalc_ks !, bstart_ks
 integer :: ikpt, ibz_k, isym_k, itim_k !isym_lgk,
 real(dp) :: dksqmax, cpu, wall, gflops
 character(len=500) :: msg
 logical :: compute_lgk
 type(lgroup_t),target :: lgk
 type(lgroup_t),pointer :: lgk_ptr
!arrays
 integer,allocatable :: iqk2dvdb(:,:)
 real(dp) :: kk(3)
 real(dp),allocatable :: kq_list(:,:)

! *************************************************************************

 ABI_UNUSED(dvdb%nqpt)

 ABI_SFREE(self%qibz_k)
 ABI_SFREE(self%wtq_k)

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 kk = self%kcalc(:, ikcalc)

 call wrtout(std_out, sjoin(ch10, repeat("=", 92)))
 msg = sjoin("[", itoa(ikcalc), "/", itoa(self%nkcalc), "]")
 call wrtout(std_out, sjoin(" Computing self-energy matrix elements for k-point:", ktoa(kk), msg))
 spin = 1
 write(msg, "(3(a, i0))")" Treating ", self%nbcalc_ks(ikcalc, spin), " band(s) between: ", &
   self%bstart_ks(ikcalc, spin)," and: ", self%bstart_ks(ikcalc, spin) + self%nbcalc_ks(ikcalc, spin) - 1
 call wrtout(std_out, msg)
 write(msg, "(2(a,i0))")"P Allocating and treating bands from my_bsum_start: ", self%my_bsum_start, &
     " up to my_bsum_stop: ", self%my_bsum_stop
 call wrtout(std_out, msg)

 ! Prepare weights for BZ(k) integration
 if (self%qint_method > 0) then
   if (self%use_doublegrid) then
     call self%ephwg%double_grid_setup_kpoint(self%eph_doublegrid, kk, prtvol, comm)
   else
     call self%ephwg%setup_kpoint(kk, prtvol, comm, skip_mapping=.true.)
   end if
   call self%ephwg%report_stats()
 endif

 call cwtime(cpu, wall, gflops, "start")

 if (self%symsigma == 0) then
   ! Do not use symmetries in BZ sum_q --> nqibz_k == nqbz
   self%nqibz_k = self%nqbz
   ABI_MALLOC(self%qibz_k, (3, self%nqibz_k))
   ABI_MALLOC(self%wtq_k, (self%nqibz_k))
   self%qibz_k = self%qbz; self%wtq_k = one / self%nqbz
   call wrtout(std_out, sjoin(" symsigma = 0 --> Integration done over full BZ with nqbz:", itoa(self%nqibz_k)))

 else if (abs(self%symsigma) == 1) then
   ! Use the symmetries of the little group of the k-point
   ! Pack points in *shells* to minimise cache misses.
   compute_lgk = .not. (self%qint_method > 0 .and. .not. self%use_doublegrid)
   if (compute_lgk) then
     lgk = lgroup_new(cryst, kk, self%timrev, self%nqbz, self%qbz, self%nqibz, self%qibz, comm)
     lgk_ptr => lgk
   else
     ! Avoid this call to lgroup new. Use lgk already computed in self%ephwg
     lgk_ptr => self%ephwg%lgk
   end if

   call wrtout(std_out, sjoin(" Number of operations in little group(k):", itoa(lgk_ptr%nsym_lg), &
     "(including time-reversal symmetry)"))
   call wrtout(std_out, sjoin(" Number of q-points in the IBZ(k):", itoa(lgk_ptr%nibz)))

   self%nqibz_k = lgk_ptr%nibz
   ABI_MALLOC(self%qibz_k, (3, self%nqibz_k))
   ABI_MALLOC(self%wtq_k, (self%nqibz_k))
   self%qibz_k = lgk_ptr%ibz; self%wtq_k = lgk_ptr%weights
   !if (compute_lgk) call lgk%free()
 else
   MSG_ERROR(sjoin("Wrong symsigma:", itoa(self%symsigma)))
 end if

 call cwtime_report(" lgroup_symsigma", cpu, wall, gflops)

 !if (self%symsigma == 0 .or. .true.) then
 if (self%symsigma == 0) then
   ! Find correspondence IBZ_k --> IBZ
   ABI_MALLOC(iqk2dvdb, (self%nqibz_k, 6))
   call listkk(dksqmax, cryst%gmet, iqk2dvdb, self%qibz, self%qibz_k, self%nqibz, self%nqibz_k, cryst%nsym, &
        1, cryst%symafm, cryst%symrec, timrev1, comm, exit_loop=.True., use_symrec=.True.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6,2a)' )&
       "At least one of the q points in the IBZ_k could not be generated from one in the IBZ. dksqmax: ",dksqmax, ch10,&
       "Action: check your DVDB file and use eph_task to interpolate the potentials on a denser q-mesh."
     MSG_ERROR(msg)
   end if
   ABI_REMALLOC(self%ind_ibzk2ibz, (6, self%nqibz_k))
   do iq_ibz=1,self%nqibz_k
     self%ind_ibzk2ibz(:, iq_ibz) = iqk2dvdb(iq_ibz, :)
   end do
   ABI_FREE(iqk2dvdb)
 elseif (abs(self%symsigma) == 1) then
   ! IBZ_k --> BZ --> IBZ
   ABI_REMALLOC(self%ind_ibzk2ibz, (6, self%nqibz_k))
   self%ind_ibzk2ibz = 0
   do ikpt=1,self%nqbz
     ibz_k    = lgk_ptr%bz2ibz_smap(1,ikpt)
     !isym_lgk = lgk_ptr%bz2ibz_smap(2,ikpt)
     !isym_k   = lgk_ptr%lgsym2glob(1,isym_lgk)
     !itim_k   = lgk_ptr%lgsym2glob(2,isym_lgk)
     isym_k   = lgk_ptr%bz2ibz_smap(2,ikpt)
     itim_k   = lgk_ptr%bz2ibz_smap(3,ikpt)
     ! I assume that isym=1 and itim_k=0 is identity but still verify the kpoint
     if (isym_k /= 1 .or. itim_k /= 1 .or. any(lgk_ptr%bz2ibz_smap(4:,ikpt)/=0)) cycle
     ! check IBZ_k --> BZ
     ABI_CHECK(sum(abs(self%qbz(:,ikpt)-self%qibz_k(:,ibz_k)))<tol8,'Wrong mapping')
     ! IBZ_k --> IBZ
     !self%ind_ibzk2ibz(:, ibz_k) = self%ind_bz2ibz(:,ikpt)
     self%ind_ibzk2ibz(1, ibz_k) = self%ind_bz2ibz(1,ikpt)
     self%ind_ibzk2ibz(2, ibz_k) = self%ind_bz2ibz(2,ikpt)
     self%ind_ibzk2ibz(6, ibz_k) = self%ind_bz2ibz(3,ikpt)
     self%ind_ibzk2ibz(3:5, ibz_k) = self%ind_bz2ibz(4:6,ikpt)
   end do
   do ikpt=1,self%nqibz_k
     ABI_CHECK(self%ind_ibzk2ibz(1,ikpt)/=0,'Did not find mapping')
   end do
   if (compute_lgk) call lgk%free()
 else
   MSG_ERROR(sjoin("Wrong symsigma:", itoa(self%symsigma)))
 endif

 call cwtime_report(" IBZ_k --> IBZ", cpu, wall, gflops)

 if (.not. self%use_ftinterp) then
   ! Find correspondence IBZ_k --> set of q-points in DVDB.
   ! Need to handle q_bz = S q_ibz by symmetrizing the potentials already available in the DVDB.
   !
   ! Note:
   !   q --> -q symmetry is always used for phonons.
   !   we use symrec instead of symrel (see also m_dvdb)
   ! IBZ_K -> BZ -> IBZ -> DVDB
   ABI_REMALLOC(self%indq2dvdb_k, (6, self%nqibz_k))
   self%indq2dvdb_k = self%ind_ibzk2ibz
   do ikpt=1,self%nqibz_k
     self%indq2dvdb_k(1,ikpt) = self%ibz2dvdb(self%ind_ibzk2ibz(1,ikpt))
   end do
   call cwtime_report(" IBZ_k --> DVDB", cpu, wall, gflops)
 end if

 ! Find k+q in the extended zone and extract symmetry info.
 ! Be careful here because there are two umklapp vectors to be considered:
 !
 !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
 !
 ABI_MALLOC(kq_list, (3, self%nqibz_k))
 do iq_ibz=1,self%nqibz_k
   kq_list(:, iq_ibz) = kk + self%qibz_k(:,iq_ibz)
 end do

 ! Use iqk2dvdb as workspace array.
 ABI_MALLOC(iqk2dvdb, (self%nqibz_k, 6))
 call listkk(dksqmax, cryst%gmet, iqk2dvdb, ebands%kptns, kq_list, ebands%nkpt, self%nqibz_k, cryst%nsym, &
   sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, comm, exit_loop=.True., use_symrec=.False.)
 if (dksqmax > tol12) then
   write(msg, '(4a,es16.6,7a)' )&
    "The WFK file cannot be used to compute self-energy corrections at k: ", trim(ktoa(kk)), ch10,&
    "At least one of the k+q points could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
    "Q-mesh: ",trim(ltoa(self%ngqpt)),", K-mesh (from kptrlatt) ",trim(ltoa(get_diag(dtset%kptrlatt))),ch10, &
    "Action: check your WFK file and (k, q) point input variables."
   MSG_ERROR(msg)
 end if
 ABI_FREE(kq_list)

 ABI_REMALLOC(self%indkk_kq, (6, self%nqibz_k))
 do iq_ibz=1,self%nqibz_k
   self%indkk_kq(:, iq_ibz) = iqk2dvdb(iq_ibz, :)
 end do
 ABI_FREE(iqk2dvdb)

 call cwtime_report(" k+q --> ebands", cpu, wall, gflops)

 if (self%qint_method > 0 .and. .not. self%use_doublegrid) then
   !self%ephwg%lgk2ibz
   ABI_REMALLOC(self%ephwg%lgk2ibz, (self%nqibz_k))
   self%ephwg%lgk2ibz = self%ind_ibzk2ibz(1, :)
   !self%ephwg%kq2ibz
   ABI_REMALLOC(self%ephwg%kq2ibz, (self%nqibz_k))
   self%ephwg%kq2ibz = self%indkk_kq(1, :)
 end if

end subroutine sigmaph_setup_kcalc
!!***

!!****f* m_sigmaph/sigmaph_setup_qloop
!! NAME
!!  sigmaph_setup_qloop
!!
!! FUNCTION
!!  Prepare integration of self-energy matrix in q-space for given (spin, ikcalc)
!!  Distribute q-points and precompute weights if tetrahedron method and imag_only
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t> = Crystal structure.
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!!  spin: spin index.
!!  ikcalc=Index of the k-point to compute.
!!  nfftf=Number of fft-points on the fine grid for interpolated potential
!!  ngfftf(18)=information on 3D FFT for interpolated potential
!!  comm= MPI communicator
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_setup_qloop(self, dtset, cryst, ebands, dvdb, spin, ikcalc, nfftf, ngfftf, comm)

!Arguments ------------------------------------
 integer,intent(in) :: spin, ikcalc, nfftf, comm
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 class(sigmaph_t),intent(inout) :: self
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb
!arrays
 integer,intent(in) :: ngfftf(18)

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: my_rank, iq_ibz_k, iq_ibz, ierr, nprocs, imyq, iq_dvdb, ii, cnt, itreat, iq
 integer :: nqeff, ndiv
 real(dp) :: cpu, wall, gflops !, weight_q
 logical :: qfilter
 character(len=500) :: msg
!arrays
 integer,allocatable :: mask_qibz_k(:), imask(:), qtab(:), ineed_qibz(:), ineed_qdvdb(:)
! real(dp) :: kk(3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 msg = "Standard quadrature"; if (self%qint_method == 1) msg = "tetrahedron method"
 call wrtout(std_out, sjoin(" Preparing q-loop with integration method:", msg))
 call cwtime(cpu, wall, gflops, "start")
 !kk = self%kcalc(:, ikcalc)

 select case (dtset%eph_task)
 case (4)
   call distribute_nqibz_k_nofilter()
   if (self%qint_method == 1) call sigmaph_get_all_qweights(self, cryst, ebands, spin, ikcalc, comm)

 case (-4)
   ! Computation of imaginary part
   if (self%qint_method == 0) then
     call distribute_nqibz_k_nofilter()

   else if (self%qint_method == 1) then
     ! Imag with tetra --> Precompute weights in IBZ_k.
     call distribute_nqibz_k_nofilter()
     call sigmaph_get_all_qweights(self, cryst, ebands, spin, ikcalc, comm)

     qfilter = any(dtset%eph_tols_idelta >= zero)
     if (qfilter) then
       ! Two-pass algorithm:
       ! Select q-points with significant contribution, recompute my_nqibz_k and myq2ibz_k.
       ! Finally, recompute integration weights with new distribution.
       ! %deltaw_pm(2, nbcalc_ks, my_npert, bsum_start:bsum_stop, my_nqibz_k, ndiv)
       ndiv = 1; if (self%use_doublegrid) ndiv = self%eph_doublegrid%ndiv
       ABI_ICALLOC(mask_qibz_k, (self%nqibz_k))
       do imyq=1,self%my_nqibz_k
         iq_ibz_k = self%myq2ibz_k(imyq)
         ! TODO: Weight should be reintroduced when Henrique's htetra is merged....
         !weight_q = self%wtq_k(iq_ibz_k)
         if (any(abs(self%deltaw_pm(1, :, :, :, imyq, :)) >= dtset%eph_tols_idelta(1) / ndiv)) mask_qibz_k(iq_ibz_k) = 1
         if (any(abs(self%deltaw_pm(2, :, :, :, imyq, :)) >= dtset%eph_tols_idelta(2) / ndiv)) mask_qibz_k(iq_ibz_k) = 1
       end do

       ! Take max inside comm.
       call alloc_copy(mask_qibz_k, imask)
       call xmpi_max(imask, mask_qibz_k, comm, ierr)
       ABI_FREE(imask)

       ! Find all qpts in the IBZ_k contributing to Im(Sigma).
       ABI_MALLOC(qtab, (self%nqibz_k))
       nqeff = 0
       do iq_ibz_k=1,self%nqibz_k
         if (mask_qibz_k(iq_ibz_k) == 1) then
           nqeff = nqeff + 1; qtab(nqeff) = iq_ibz_k
         end if
       end do
       ABI_FREE(mask_qibz_k)

       if (my_rank == master) then
         write(std_out, "(a, 2(es16.6,1x))")" Removing q-points with integration weights < ", dtset%eph_tols_idelta / ndiv
         write(std_out, "(a,i0,a,f5.1,a)")" Total number of q-points contributing to Im(Sigma): ", nqeff, &
           " (nqeff / nqibz_k): ", (100.0_dp * nqeff) / self%nqibz_k, " [%]"
       end if

       ! Redistribute relevant q-points inside comm_qpt taking into account itreat_qibz
       ! I may need to update qcache after this operation -->  build ineed_* table.
       ! Must handle two cases: potentials from DVDB or Fourier-interpolated.
       if (self%use_ftinterp) then
         ABI_ICALLOC(ineed_qibz, (self%nqibz))
       else
         ABI_ICALLOC(ineed_qdvdb, (dvdb%nqpt))
       end if

       self%my_nqibz_k = 0
       do ii=1,2
         if (ii == 2) then
           ABI_REMALLOC(self%myq2ibz_k, (self%my_nqibz_k))
         end if
         cnt = 0
         do iq=1,nqeff
           iq_ibz_k = qtab(iq)
           iq_ibz = self%ind_ibzk2ibz(1, iq_ibz_k)
           itreat = int(self%itreat_qibz(iq_ibz), kind=i4b)
           if (.not. self%use_ftinterp) iq_dvdb = self%indq2dvdb_k(1, iq_ibz_k)
           if (itreat /= 0) then
             !if (itreat == 1 .or. (itreat > 1 .and. itreat == 1 + iqpos(self%qibz_k(:, iq_ibz_k)))) then
             if (ii == 1) self%my_nqibz_k = self%my_nqibz_k + 1
             if (ii == 2) then
               cnt = cnt + 1
               self%myq2ibz_k(cnt) = qtab(iq)
               if (self%use_ftinterp) then
                 if (.not. allocated(dvdb%ft_qcache%key(iq_ibz)%v1scf)) ineed_qibz(iq_ibz) = 1
               else
                 if (.not. allocated(dvdb%qcache%key(iq_dvdb)%v1scf)) ineed_qdvdb(iq_dvdb) = 1
               end if
             end if
           end if
         end do
       end do
       ABI_FREE(qtab)

       ! Recompute weights with new q-point distribution.
       call sigmaph_get_all_qweights(self, cryst, ebands, spin, ikcalc, comm)

       write(msg, "(a,i0,2a,f7.3,a)") &
        " Number of q-points in IBZ(k) treated by this proc: ", self%my_nqibz_k, ch10, &
        " Load balance inside comm_qpt: ", (one * self%nqibz_k) / (self%my_nqibz_k * self%nprocs_qpt), " (should be ~1)"
       call wrtout(std_out, msg)
       MSG_WARNING_IF(self%my_nqibz_k == 0, "my_nqibz_k == 0")

       ! Make sure each node has the q-points we need. Try not to break qcache_size_mb contract!
       if (self%use_ftinterp) then
         ! Update cache by Fourier interpolating W(r,R)
         call dvdb%ftqcache_update_from_ft(nfftf, ngfftf, self%nqibz, self%qibz, ineed_qibz, comm)
       else
         ! Update cache. Perform collective IO inside comm if needed.
         call dvdb%qcache_update_from_file(nfftf, ngfftf, ineed_qdvdb, comm)
       end if

       ABI_SFREE(ineed_qibz)
       ABI_SFREE(ineed_qdvdb)
     end if ! qfilter

   else
     MSG_ERROR(sjoin("Invalid eph_intmeth:", itoa(dtset%eph_intmeth)))
   end if ! intmeth

 case default
   MSG_ERROR(sjoin("Invalid eph_task:", itoa(dtset%eph_task)))
 end select

 call cwtime_report(" Setup qloop", cpu, wall, gflops)

contains

subroutine distribute_nqibz_k_nofilter()
  ! Find number of q-points in IBZ(k) treated by this MPI rank
  ! taking into account itreat_qibz and build redirection table myq2ibz_k.
  ! The distribution must be consistent with the WF distribution done with bks_mask

  self%my_nqibz_k = 0
  do ii=1,2
    if (ii == 2) then
      ABI_REMALLOC(self%myq2ibz_k, (self%my_nqibz_k))
    end if
    cnt = 0
    do iq_ibz_k=1,self%nqibz_k
      iq_ibz = self%ind_ibzk2ibz(1, iq_ibz_k)
      if (self%itreat_qibz(iq_ibz) == 0) cycle
      if (ii == 1) self%my_nqibz_k = self%my_nqibz_k + 1
      if (ii == 2) then
        cnt = cnt + 1
        self%myq2ibz_k(cnt) = iq_ibz_k
      end if
    end do
  end do

end subroutine distribute_nqibz_k_nofilter

end subroutine sigmaph_setup_qloop
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

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, spin, prtvol, comm
 class(sigmaph_t),target,intent(inout) :: self
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 integer,parameter :: master=0, max_ntemp=1
 integer :: ideg,ib,it,ii,iw,nstates,ierr,my_rank,band_ks,ik_ibz,ibc,ib_val,ib_cond,jj
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_oms,qpe_oms_val,qpe_oms_cond
 real(dp) :: cpu, wall, gflops, invsig2fmts, tau
 complex(dpc) :: sig0c,sig0fr,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2
 !character(len=500) :: msg
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
!arrays
 !integer :: shape3(3),shape4(4),shape5(5),shape6(6)
 integer, ABI_CONTIGUOUS pointer :: bids(:)
 !real(dp), ABI_CONTIGUOUS pointer :: rdata3(:,:,:), rdata4(:,:,:,:), rdata5(:,:,:,:,:), rdata6(:,:,:,:,:,:)
 real(dp) :: qp_gaps(self%ntemp),qpoms_gaps(self%ntemp)
 real(dp),allocatable :: aw(:,:,:), frohl_aw(:,:,:), a2few_avg(:,:)
 real(dp) :: ks_enes(self%max_nbcalc),ze0_vals(self%ntemp, self%max_nbcalc)
 real(dp) :: gfw_avg(self%gfw_nomega, 3)
 complex(dpc) :: qpoms_enes(self%ntemp, self%max_nbcalc),qp_enes(self%ntemp, self%max_nbcalc)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 ! Could use non-blocking communications and double buffer technique to reduce synchronisation cost...
 call cwtime(cpu, wall, gflops, "start", msg=" Gathering results. Waiting for other processes...")

 call xmpi_sum_master(self%vals_e0ks, master, comm, ierr)
 call xmpi_sum_master(self%dvals_de0ks, master, comm, ierr)
 call xmpi_sum_master(self%dw_vals, master, comm, ierr)
 if (self%nwr > 0) call xmpi_sum_master(self%vals_wr, master, comm, ierr)
 if (self%calc_mrta) call xmpi_sum_master(self%linewidth_mrta, master, comm, ierr)

 call cwtime_report(" Sigma_nk gather completed", cpu, wall, gflops, comm=comm)
 ! Only master writes
 if (my_rank /= master) return

 ik_ibz = self%kcalc2ibz(ikcalc, 1)

 ! Add Frohlich results to final results (frohl arrays are not symmetrized if symsigma == 1)
 if (self%frohl_model == 2) then
   self%vals_e0ks = self%vals_e0ks + self%frohl_vals_e0ks
   if (.not. self%imag_only) then
     self%dvals_de0ks = self%dvals_de0ks + self%frohl_dvals_de0ks
   end if
   if (self%nwr > 0) self%vals_wr = self%vals_wr + self%frohl_vals_wr
 end if

 if (self%a2f_ne > 0) then
   ABI_MALLOC(a2few_avg, (self%a2f_ne, self%gfw_nomega))
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
       if (self%a2f_ne > 0) then
          a2few_avg = sum(self%a2few(:, :, bids(:)), dim=3) / nstates
          do ii=1,nstates
            self%a2few(:, :, bids(ii)) = a2few_avg
          end do
       end if
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

       ! Average TAU_MRTA
       if (self%calc_mrta) then
         ravg = sum(self%linewidth_mrta(it, bids(:))) / nstates
         do ii=1,nstates
           self%linewidth_mrta(it, bids(ii)) = ravg
         end do
       end if

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

 ABI_SFREE(a2few_avg)

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
   !if (self%frohl_model == 1) then
   !  write(ab_out,"(a)")"   SF1(eKS): Real part of the (model) Frohlich self-energy at the KS energy, SF2 for imaginary part."
   !end if
   write(ab_out,"(a)")"   Z(eKS): Renormalization factor."
   write(ab_out,"(a)")"   FAN: Real part of the Fan term at eKS. DW: Debye-Waller term."
   write(ab_out,"(a)")"   DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP."
   write(ab_out,"(a)")"   OTMS: On-the-mass-shell approximation with eQP ~= eKS + Sigma(omega=eKS)"
   write(ab_out,"(a)")"   TAU(eKS): Lifetime in femtoseconds computed at the KS energy."
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
 end if

 do it=1,self%ntemp
   ! Write header.
   if (it <= max_ntemp) then
     if (self%nsppol == 1) then
       write(ab_out,"(3a,f6.1,a)") &
         "K-point: ", trim(ktoa(self%kcalc(:,ikcalc))), ", T= ", self%kTmesh(it) / kb_HaK, " [K]"
     else
       write(ab_out,"(3a,i1,a,f6.1,a)") &
         "K-point: ", trim(ktoa(self%kcalc(:,ikcalc))), ", spin: ", spin, ", T= ",self%kTmesh(it) / kb_HaK, " [K]"
     end if
     if (self%imag_only) then
       !if (self%frohl_model == 0) then
         write(ab_out,"(a)")"   B    eKS    SE2(eKS)  TAU(eKS)  DeKS"
       !else
       !  write(ab_out,"(a)")"   B    eKS    SE2(eKS)  SF2(eKS)  TAU(eKS)  DeKS"
       !end if
     else
       !if (self%frohl_model == 0) then
         write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
       !else
       !  write(ab_out,"(a)")&
       !    "   B    eKS     eQP    eQP-eKS   SE1(eKS)  SF1(eKS)  SE2(eKS)  SF2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
       !end if
     end if
   end if

   do ibc=1,self%nbcalc_ks(ikcalc,spin)
     band_ks = self%bstart_ks(ikcalc,spin) + ibc - 1
     kse = ebands%eig(band_ks, ik_ibz, spin)
     ks_enes(ibc) = kse
     sig0c = self%vals_e0ks(it, ibc)
     dw = self%dw_vals(it, ibc)
     fan0 = real(sig0c) - dw
     if (self%frohl_model == 2) sig0fr = self%frohl_vals_e0ks(it, ibc)
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
     if (it <= max_ntemp) then
       if (self%imag_only) then
         invsig2fmts = Time_Sec * 1e+15 / two
         tau = 999999.0_dp
         if (abs(aimag(sig0c)) > tol16) tau = invsig2fmts / abs(aimag(sig0c))
         tau = min(tau, 999999.0_dp)
         !if (self%frohl_model == 0) then
           ! B    eKS     SE2(eKS)  TAU(eKS)  DeKS
           write(ab_out, "(i4,2(f8.3,1x),f8.1,1x,f8.3)") &
             band_ks, kse * Ha_eV, aimag(sig0c) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
         !else
         !  ! B    eKS     SE2(eKS)  SF2(eKS)  TAU(eKS)  DeKS
         !  write(ab_out, "(i4,2(f8.3,1x),2(f8.1,1x),f8.3)") &
         !    band_ks, kse * Ha_eV, aimag(sig0c) * Ha_eV, aimag(sig0fr) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
         !end if
       else
         !if (self%frohl_model == 0) then
           ! B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP
           write(ab_out, "(i4, 10(f8.3,1x))") &
             band_ks, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
             real(sig0c) * Ha_eV, aimag(sig0c) * Ha_eV, real(zc), &
             fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
         !else
         !  ! B    eKS     eQP    eQP-eKS   SE1(eKS)  SF1(eKS) SE2(eKS) SF2(eKS) Z(eKS)  FAN(eKS)   DW      DeKS     DeQP
         !  write(ab_out, "(i4, 12(f8.3,1x))") &
         !    band_ks, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
         !    real(sig0c) * Ha_eV, real(sig0fr) * Ha_eV, aimag(sig0c) * Ha_eV, aimag(sig0fr) * Ha_eV, real(zc), &
         !    fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
         !end if
       end if
     end if
     if (ibc > 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     qpoms_enes(it, ibc) = qpe_oms
     qp_enes(it, ibc) = qpe
     if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
       ! We have enough states to compute the gap.
       if (it == 1) ks_gap = kse_cond - kse_val
       qpoms_gaps(it) = qpe_oms_cond - qpe_oms_val
       qp_gaps(it) = real(qpe_cond - qpe_val)
     end if
   end do ! ibc

   ! Print KS and QP gaps.
   if (it <= max_ntemp) then
     if (.not. self%imag_only) then
       if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
         write(ab_out, "(a)")" "
         write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, &
           "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
         write(ab_out, "(2(a,f8.3),a)")" QP gap: ",qp_gaps(it) * Ha_eV," (OTMS: ",qpoms_gaps(it) * Ha_eV, ")"
         write(ab_out, "(2(a,f8.3),a)")" QP_gap - KS_gap: ",(qp_gaps(it) - ks_gap) * Ha_eV,&
             " (OTMS: ",(qpoms_gaps(it) - ks_gap) * Ha_eV, ")"
         write(ab_out, "(a)")" "
       end if
     else
       if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
         write(ab_out, "(a)")" "
         write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
         write(ab_out, "(a)")" "
       end if
     end if

   write(ab_out, "(a)")repeat("=", 92)
   end if

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

 if (self%calc_mrta) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "linewidth_mrta"), self%linewidth_mrta, start=[1,1,ikcalc,spin]))
 end if

 if (self%frohl_model == 1 .and. self%imag_only) then
   ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_deltas_sphcorr"), &
      self%frohl_deltas_sphcorr, start=[1,1,1,1, ikcalc, spin])
   NCF_CHECK(ncerr)
 end if

 !if (self%frohl_model == 2) then
 !  ! Write Frohlich self-energy to file.
 !  ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_vals_e0ks"), &
 !     c2r(self%frohl_vals_e0ks), start=[1,1,1,ikcalc,spin])
 !  NCF_CHECK(ncerr)
 !  ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_dvals_de0ks"), &
 !     c2r(self%frohl_dvals_de0ks), start=[1,1,1,ikcalc,spin])
 !  NCF_CHECK(ncerr)
 !  if (self%nwr > 0) then
 !    ncerr = nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_vals_wr"), &
 !      c2r(self%frohl_vals_wr), start=[1,1,1,1,ikcalc,spin])
 !    NCF_CHECK(ncerr)
 !  end if
 !end if

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
       if (self%frohl_model == 2) then
         ! Spectral function associated to Frohlich model.
         frohl_aw(:, it, ib) = -piinv * aimag(self%frohl_vals_wr(:, it, ib)) / &
           ((self%wrmesh_b(:, ib) - kse - real(self%frohl_vals_wr(:, it, ib))) ** 2 + aimag(self%frohl_vals_wr(:, it, ib)) ** 2)
       end if
     end do
   end do
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "spfunc_wr"), aw, start=[1, 1, 1, ikcalc, spin]))
   !if (self%frohl_model == 2) then
   !  NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "frohl_spfunc_wr"), frohl_aw, start=[1, 1, 1, ikcalc, spin]))
   !end if
   ABI_FREE(aw)
   ABI_FREE(frohl_aw)
 end if

 ! Write Eliashberg functions
 if (self%gfw_nomega > 0) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "gfw_vals"), self%gfw_vals, start=[1, 1, 1, ikcalc, spin]))
 end if
 if (self%gfw_nomega > 0 .and. self%a2f_ne > 0) then
   NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "a2few"), self%a2few, start=[1, 1, 1, ikcalc, spin]))
 end if

 ! Write array for restart.
 self%qp_done(ikcalc,spin) = 1
 NCF_CHECK(nf90_put_var(self%ncid, nctk_idname(self%ncid, "qp_done"), 1, start=[ikcalc, spin]))

 ! Dump the cache to file (this is necessary to ensure we can restart)
 NCF_CHECK(nf90_sync(self%ncid))
#endif

 call cwtime_report(" Sigma_nk netcdf output", cpu, wall, gflops)

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

!Arguments ------------------------------------
 integer,intent(in) :: unt
 type(dataset_type),intent(in) :: dtset
 class(sigmaph_t),intent(in) :: self

!Local variables-------------------------------
 integer :: ikc, is, ndiv
 character(len=500) :: msg

! *************************************************************************

 if (unt == dev_null) return

 ! Write dimensions
 !write(unt,"(a, 2(f5.3, 1x), a)")"Computing self-energy corrections for states inside energy window:", &
 !    self%elow * Ha_eV, self%ehigh * Ha_eV, "[eV]"
 write(unt,"(a)")sjoin(" Number of bands in e-ph self-energy sum:", itoa(self%nbsum))
 write(unt,"(a)")sjoin(" From bsum_start:", itoa(self%bsum_start), "to bsum_stop:", itoa(self%bsum_stop))
 if (dtset%eph_stern == 1 .and. .not. self%imag_only) then
   write(unt, "(a)")" Treating high-energy bands with Sternheimer and static self-energy."
   write(unt, "(a, es16.6, a, i0)")" Tolwfr:", dtset%tolwfr, ", nline:", dtset%nline
 end if
 write(unt,"(a)")sjoin(" Symsigma: ",itoa(self%symsigma), "Timrev:", itoa(self%timrev))
 write(unt,"(a)")sjoin(" Imaginary shift in the denominator (zcut): ", ftoa(aimag(self%ieta) * Ha_eV, fmt="f5.3"), "[eV]")
 msg = " Standard quadrature"; if (self%qint_method == 1) msg = " Tetrahedron method"
 write(unt, "(2a)")sjoin(" Method for q-space integration:", msg)
 if (self%qint_method == 1) then
   ndiv = 1; if (self%use_doublegrid) ndiv = self%eph_doublegrid%ndiv
   write(unt, "(a, 2(es16.6,1x))")" Tolerance for integration weights < ", dtset%eph_tols_idelta(:) / ndiv
 end if
 if (self%use_doublegrid) write(unt, "(a, i0)")" Using double grid technique with ndiv: ", self%eph_doublegrid%ndiv
 if (self%imag_only) write(unt, "(a)")" Only the Imaginary part of Sigma will be computed."
 if (.not. self%imag_only) write(unt, "(a)")" Both Real and Imaginary part of Sigma will be computed."
 write(unt,"(a)")sjoin(" Number of frequencies along the real axis:", itoa(self%nwr), &
    ", Step:", ftoa(self%wr_step * Ha_eV, fmt="f5.3"), "[eV]")
 write(unt, "(a)")sjoin(" Number of frequency in generalized Eliashberg functions:", itoa(self%gfw_nomega))
 write(unt,"(a)")sjoin(" Number of temperatures:", itoa(self%ntemp), &
   "From:", ftoa(self%kTmesh(1) / kb_HaK), "to", ftoa(self%kTmesh(self%ntemp) / kb_HaK), "[K]")
 write(unt,"(a)")sjoin(" Ab-initio q-mesh from DDB file:", ltoa(dtset%ddb_ngqpt))
 write(unt,"(a)")sjoin(" Ab-initio q-mesh from DVDB file:", ltoa(dtset%dvdb_ngqpt))
 write(unt,"(a)")sjoin(" Q-mesh used for self-energy integration [ngqpt]:", ltoa(self%ngqpt))
 write(unt,"(a)")sjoin(" Number of q-points in the IBZ:", itoa(self%nqibz))
 write(unt,"(a)")sjoin(" asr:", itoa(dtset%asr), "dipdip:", itoa(dtset%dipdip), "symdynmat:", itoa(dtset%symdynmat))
 if (self%frohl_model == 0) then
   write(unt,"(a)")" No special treatment of Frohlich divergence in gkq for q --> 0"
 else if (self%frohl_model == 1) then
   write(unt,"(a)")" Integrating Frohlich model in small sphere around Gamma to accelerate qpt convergence"
   write(unt,"(2(a,i0,1x))")" Sperical integration performed with: ntheta: ", self%ntheta, ", nphi: ", self%nphi
 else if (self%frohl_model == 2) then
   !write(unt,"((a,i0,1x,a,f5.3,1x,a))")"nr points:", self%nqr, "qrad:", self%qrad, "[Bohr^-1]"
 end if
 write(unt,"(a, i0)")" Number of k-points for self-energy corrections: ", self%nkcalc
 if (all(dtset%sigma_erange /= -one)) then
   write(unt, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
 end if
 write(unt,"(a)")" List of K-points for self-energy corrections:"
 do ikc=1,self%nkcalc
   if (ikc > 10) then
     write(unt, "(2a)")" nkcalc > 10. Stop printing more k-point information.",ch10
     exit
   end if
   do is=1,self%nsppol
     if (self%nsppol == 2) write(unt,"(a,i1)")" For spin: ",is
     write(unt, "(2(i4,2x),a,2(i4,1x))") &
       ikc, is, trim(ktoa(self%kcalc(:,ikc))), self%bstart_ks(ikc,is), self%bstart_ks(ikc,is) + self%nbcalc_ks(ikc,is) - 1
     end do
 end do

 write(unt, "(a)")" === MPI parallelism ==="
 write(unt, "(2(a,i0))")"P Allocating and summing bands from my_bsum_start: ", self%my_bsum_start, &
     " up to my_bsum_stop: ", self%my_bsum_stop
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over perturbations: ", self%nprocs_pert
 write(unt, "(a,i0)")"P Number of perturbations treated by this CPU: ", self%my_npert
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over q-points: ", self%nprocs_qpt
 write(unt, "(2(a,i0))")"P Number of q-points in the IBZ treated by this proc: " , &
     count(self%itreat_qibz == 1), " of ", self%nqibz
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over bands: ", self%nprocs_bsum

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

subroutine sigmaph_get_all_qweights(sigma, cryst, ebands, spin, ikcalc, comm)

!Arguments ------------------------------------
!scalars
 class(sigmaph_t),intent(inout) :: sigma
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: ikcalc, spin, comm

!Local variables ------------------------------
!scalars
 integer :: nu, ibsum_kq, ik_ibz, bstart_ks, nbcalc_ks, my_rank, natom3 !, ierr
 integer :: nprocs, imyp, imyq, ndiv, bsum_start, bsum_stop, ib_k, band_ks
 integer :: iq_ibz_fine,iq_bz_fine,iq_ibz,jj, nz
 real(dp) :: weight, cpu,wall, gflops
 real(dp) :: eig0nk
!arrays
 real(dp) :: kk(3), kq(3), qpt(3), dpm(2)
 real(dp),allocatable :: tmp_deltaw_pm(:,:,:)
 complex(dpc),allocatable :: zvals(:,:), tmp_cweights(:,:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
 nbcalc_ks = sigma%nbcalc_ks(ikcalc,spin)
 bstart_ks = sigma%bstart_ks(ikcalc,spin)
 bsum_start = sigma%bsum_start; bsum_stop = sigma%bsum_stop
 natom3 = 3 * cryst%natom
 ndiv = 1; if (sigma%use_doublegrid) ndiv = sigma%eph_doublegrid%ndiv

 ABI_CHECK(sigma%symsigma == 1, "symsigma 0 with tetra not implemented")

 if (sigma%imag_only) then
   ! Weights for Im (tethraedron, eta --> 0)
   ABI_REMALLOC(sigma%deltaw_pm, (2, nbcalc_ks, sigma%my_npert, bsum_start:bsum_stop, sigma%my_nqibz_k, ndiv))
   sigma%deltaw_pm = zero

   ! Temporary weights (on the fine IBZ_k mesh if double grid is used)
   ABI_MALLOC(tmp_deltaw_pm, (1, sigma%ephwg%nq_k, 2))

   ! loop over bands to sum
   do ibsum_kq=sigma%bsum_start, sigma%bsum_stop
    ! loop over my phonon modes
    do imyp=1,sigma%my_npert
      nu = sigma%my_pinfo(3, imyp)

      ! HM: This one should be faster but uses more memory, I compute for each ib instead
      ! Compute weights inside comm_bq
      !call sigma%ephwg%get_deltas_wvals(ibsum_kq, spin, nu, nbcalc_ks, &
      !                                  ebands%eig(bstart_ks:bstart_ks+nbcalc_ks, ik_ibz, spin), &
      !                                  sigma%bcorr, tmp_deltaw_pm, sigma%comm_bq)

      ! loop over bands in self-energy matrix elements.
      do ib_k=1,nbcalc_ks
        band_ks = ib_k + bstart_ks - 1
        eig0nk = ebands%eig(band_ks, ik_ibz, spin)

        ! Compute weights inside comm_bq
        call sigma%ephwg%get_deltas_wvals(ibsum_kq, spin, nu, 1, [eig0nk], &
                                          sigma%bcorr, tmp_deltaw_pm, sigma%comm_bq)

        ! For all the q-points that I am going to calculate
        do imyq=1,sigma%my_nqibz_k
          iq_ibz = sigma%myq2ibz_k(imyq)

          if (sigma%use_doublegrid) then
            ! For all the q-points in the microzone
            ! This is done again in the main sigmaph routine
            qpt = sigma%qibz_k(:,iq_ibz)
            call sigma%eph_doublegrid%get_mapping(kk, kq, qpt)
            do jj=1,sigma%eph_doublegrid%ndiv
              iq_bz_fine = sigma%eph_doublegrid%mapping(3,jj)
              iq_ibz_fine = sigma%eph_doublegrid%bz2lgkibz(iq_bz_fine)
              weight = sigma%ephwg%lgk%weights(iq_ibz_fine)
              !dpm = tmp_deltaw_pm(ib_k, iq_ibz_fine, :)
              dpm = tmp_deltaw_pm(1, iq_ibz_fine, :)
              sigma%deltaw_pm(:, ib_k, imyp, ibsum_kq, imyq, jj) = dpm / weight
            end do
          else
            weight = sigma%ephwg%lgk%weights(iq_ibz)
            !dpm = tmp_deltaw_pm(ib_k, iq_ibz, :)
            dpm = tmp_deltaw_pm(1, iq_ibz, :)
            sigma%deltaw_pm(:, ib_k, imyp, ibsum_kq, imyq, 1) = dpm / weight
          end if

        end do
      end do
    end do
   end do

   ABI_FREE(tmp_deltaw_pm)

 else
   ! Both real and imag part --> compute \int 1/z with tetrahedron.
   ! Note that we still need a finite i.eta in the expression (hopefully smaller than default value).
   ! Besides we have to take into account the case in which the spectral function in wanted.
   ! Derivative wrt omega is still computed with finite i.eta.
   ABI_CHECK(.not. sigma%use_doublegrid, "double grid for Re-Im not implemented")

   ! TODO: This part should be tested.
   nz = 1; if (sigma%nwr > 0) nz = 1 + sigma%nwr
   ABI_REMALLOC(sigma%cweights, (nz,2,nbcalc_ks,sigma%my_npert,sigma%my_bsum_start:sigma%my_bsum_stop,sigma%my_nqibz_k,ndiv))
   ABI_MALLOC(tmp_cweights, (nz, 2, nbcalc_ks, sigma%nqibz_k))
   ABI_MALLOC(zvals, (nz, nbcalc_ks))

   ! Initialize z-points for Sigma_{nk} for different n bands.
   zvals(1, :) = sigma%e0vals + sigma%ieta
   if (sigma%nwr > 0) zvals(2:sigma%nwr+1, :) = sigma%wrmesh_b(:, 1:nbcalc_ks) + sigma%ieta
   ! Loop over my bands in self-energy sum.
   ! TODO: Really slow if nz >> 1, reduce the number of ibsum_kq bands for which tetra must be used.
   do ibsum_kq=sigma%my_bsum_start, sigma%my_bsum_stop
     ! Loop over my phonon modes
     do imyp=1,sigma%my_npert
       nu = sigma%my_pinfo(3, imyp)

       ! cweights(nz, 2, nbsigma, self%nq_k)
       call sigma%ephwg%get_zinv_weights(nz, nbcalc_ks, zvals, ibsum_kq, spin, nu, tmp_cweights, xmpi_comm_self)
       !  use_bzsum=sigma%symsigma == 0)

       ! For all the q-points that I am going to calculate
       do imyq=1,sigma%my_nqibz_k
         iq_ibz = sigma%myq2ibz_k(imyq)
         weight = sigma%ephwg%lgk%weights(iq_ibz)
         sigma%cweights(:, :, :, imyp, ibsum_kq, imyq, 1) = tmp_cweights(:, :, :, iq_ibz) / weight
       end do
     end do
   end do
   ABI_FREE(zvals)
   ABI_FREE(tmp_cweights)
 end if

 call cwtime_report(" weights with tetrahedron", cpu, wall, gflops)

end subroutine sigmaph_get_all_qweights
!!***

!!****f* m_sigmaph/eval_sigfrohl_deltas
!! NAME
!!  eval_sigfrohl_deltas
!!
!! FUNCTION
!!  Compute frohl_deltas_sphcorr array with contributions to the imaginary part
!!  of the Fan-Migdal self-energy due the Frohlich divergence.
!!  Use Verdi's model, group velocities and radial integration.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine eval_sigfrohl_deltas(sigma, cryst, ifc, ebands, ikcalc, spin, prtvol, comm)

!Arguments ------------------------------------
 class(sigmaph_t),intent(inout) :: sigma
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: ikcalc, spin, prtvol, comm

!Local variables ------------------------------
!scalars
 real(dp),parameter :: mytol = tol12
 integer :: nu, it, iatom, my_rank, nprocs !nbcalc_ks,
 integer :: ib_k, band_ks, ik_ibz, iang, ierr
 real(dp) :: wqnu, nqnu, eig0nk, f_nk, dfde_nk
 real(dp) :: inv_qepsq, q0rad, vnk_mod, cos_theta_qvnk, qroot, fact_qvers !, den
! real(dp) :: cpu, wall, gflops
 complex(dpc) :: cnum
!arrays
 real(dp) :: qvers_cart(3), vnk(3), nqnu_tlist(sigma%ntemp)
 real(dp) :: phfrq(cryst%natom*3), displ_cart(2,3,cryst%natom,3*cryst%natom)
 complex(dpc) :: cp3(3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Integrate delta functions inside the small sphere around the Gamma point.
 ! Radius of sphere with volume equivalent to the micro zone.
 q0rad = two_pi * (three / (four_pi * cryst%ucvol * sigma%nqbz)) ** third

 ABI_MALLOC_IFNOT(sigma%frohl_deltas_sphcorr, (2, sigma%ntemp, sigma%max_nbcalc, 3 * cryst%natom))
 sigma%frohl_deltas_sphcorr = zero

 !kk = sigma%kcalc(:, ikcalc)
 ik_ibz = sigma%kcalc2ibz(ikcalc, 1)

 ! Compute angular average.
 do iang=1,sigma%angl_size
   if (mod(iang, nprocs) /= my_rank) cycle ! MPI parallelism
   qvers_cart = sigma%qvers_cart(:, iang)
   inv_qepsq = one / dot_product(qvers_cart, matmul(ifc%dielt, qvers_cart))

   ! Compute phonons with NA behaviour along qvers_cart.
   call ifc%fourq(cryst, qvers_cart, phfrq, displ_cart, nanaqdir="cart")

   ! Note that acoustic modes are ignored.
   do nu=4,3*cryst%natom
     wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
     ! Get phonon occupation for all temperatures.
     nqnu_tlist = occ_be(wqnu, sigma%kTmesh(:), zero)

     cp3 = czero
     do iatom=1, cryst%natom
      cp3 = cp3 + matmul(ifc%zeff(:, :, iatom), cmplx(displ_cart(1,:,iatom, nu), displ_cart(2,:,iatom, nu), kind=dpc))
     end do
     cnum = dot_product(qvers_cart, cp3)

     ! For each band in sigma_(nk)
     do ib_k=1,sigma%nbcalc_ks(ikcalc, spin)
       band_ks = ib_k + sigma%bstart_ks(ikcalc, spin) - 1
       vnk = sigma%vcar_calc(:, ib_k, ikcalc, spin)
       vnk_mod = sqrt(dot_product(vnk, vnk))
       if (abs(vnk_mod) < mytol) cycle
       cos_theta_qvnk = dot_product(qvers_cart, vnk) / vnk_mod
       if (abs(cos_theta_qvnk) < mytol) cycle
       ! The two possible roots (+- wqn / (v cos_theta_qvnk) must be inside the small sphere.
       qroot = wqnu / (vnk_mod * cos_theta_qvnk)
       if (abs(qroot) > q0rad) cycle

       ! Use group velocities to expand e_kq and f_kq around (n, k)
       ! NB: Don't know how to treat degeneracies correctly
       ! Phonon quantities are not Taylor expanded, we just take into accout the angular dependence.
       ! As a matter of fact, we have a non-analytical point and ph freqs are usually flat in this small region.
       eig0nk = ebands%eig(band_ks, ik_ibz, spin)

       ! Extract part of the integrand that does not depend on T.
       !den = wqnu ** 2 * vnk_mod * cos_theta_qvnk; if (abs(den) < mytol) cycle
       fact_qvers = sigma%angwgth(iang) * abs(cnum) ** 2 * inv_qepsq ** 2 * vnk_mod * cos_theta_qvnk / wqnu ** 3

       do it=1,sigma%ntemp
         nqnu = nqnu_tlist(it)
         ! f_{k+q} = df/dek vk.q
         f_nk = occ_fd(eig0nk, sigma%kTmesh(it), sigma%mu_e(it))
         !dfde_nk = occ_dfd(eig0nk, sigma%kTmesh(it), sigma%mu_e(it))
         ! TODO: I got NAN for T --> 0
         !dfde_nk = zero

         ! Different expressions for absorption and emission.
         ! TODO: Check that terms are always >= 0
         if (qroot >= zero) then
           ! cos_theta >= 0
           sigma%frohl_deltas_sphcorr(1, it, ib_k, nu) = sigma%frohl_deltas_sphcorr(1, it, ib_k, nu) + &
             fact_qvers * (nqnu + f_nk + dfde_nk * wqnu)
         else
           ! cos_theta < 0
           sigma%frohl_deltas_sphcorr(2, it, ib_k, nu) = sigma%frohl_deltas_sphcorr(2, it, ib_k, nu) + &
             fact_qvers * (nqnu + one - f_nk + dfde_nk * wqnu)
         end if
       end do

     end do ! ib_k
   end do ! nu
 end do ! iang

 sigma%frohl_deltas_sphcorr = - four_pi * sigma%frohl_deltas_sphcorr / (two * pi ** 2 * cryst%ucvol)
 call xmpi_sum(sigma%frohl_deltas_sphcorr, comm, ierr)

 if (my_rank == 0 .and. prtvol > 1) then
   write(std_out, *)" eval_sigfrohl_deltas: f+d+, f-d-, sum(f+d+, f-d-) for nu in 4, ... 3natom"
   do ib_k=1,sigma%nbcalc_ks(ikcalc, spin)
     do nu=4,3*cryst%natom
       write(std_out, *)sum(sigma%frohl_deltas_sphcorr(:, 1, ib_k, nu), dim=1), sigma%frohl_deltas_sphcorr(:, 1, ib_k, nu)
     end do
   end do
 end if

end subroutine eval_sigfrohl_deltas
!!***

!!****f* m_sigmaph/eval_sigfrohl2
!! NAME
!!  eval_sigfrohl2
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine eval_sigfrohl2(sigma, cryst, ifc, ebands, ikcalc, spin, comm)

!Arguments ------------------------------------
 class(sigmaph_t),intent(inout) :: sigma
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: ebands
 integer,intent(in) :: ikcalc, spin, comm

!Local variables ------------------------------
!scalars
 integer :: my_rank, nprocs, nbcalc_ks, nu, it, natom3, iatom, bstart_ks !i1, i2, i3,  iq_ibz,
 integer :: ib_k, band_ks, ik_ibz, iang, iqr, ierr, iw
 real(dp) :: wqnu, eig0nk, qrad2 !, weight_q !, vol_fact !eig0mkq, f_mkq , gkq2,  nqnu,
 real(dp) :: inv_qepsq, qstep , qmod, fqdamp !gmod2, hmod2,  rfact,
 real(dp) :: cpu_fr,wall_fr,gflops_fr
 complex(dpc) :: cfact, cnum
!arrays
 !integer :: ndivs(3)
 real(dp) :: qpt(3), qpt_cart(3), kq(3), kk(3) !q0(3),
 real(dp) :: qrmesh(sigma%nqr), gmod2r(sigma%nqr), hmod2r(sigma%nqr), rfactr(sigma%nqr), nqr(sigma%nqr)
 real(dp) :: phfrq(cryst%natom*3)
 real(dp),allocatable :: displ_cart(:,:,:,:), eigs_kq(:), eigs_kqr(:,:), wqr(:,:)
 real(dp),allocatable :: f_mkqr(:), eig0mkqr(:), gkqr2(:,:)
 complex(dpc) :: cdd(3), cfqr(sigma%nqr)
 complex(dpc),allocatable :: cfact_wr(:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call cwtime(cpu_fr, wall_fr, gflops_fr, "start")

 natom3 = 3 * cryst%natom
 nbcalc_ks = sigma%nbcalc_ks(ikcalc, spin)
 bstart_ks = sigma%bstart_ks(ikcalc, spin)
 ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
 kk = sigma%kcalc(:, ikcalc)

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
 !return

#if 0
 ndivs = 1
 do iq_ibz=1,sigma%nqibz_k
    if (mod(iq_ibz, nprocs) /= my_rank) cycle ! MPI parallelism
    !write(std_out,*)"iq_ibz:", iq_ibz
    q0 = sigma%qibz_k(:,iq_ibz)
    ! q-weigths for IBZ(k) integration
    weight_q = sigma%wtq_k(iq_ibz) ! * (four_pi / cryst%ucvol)

    ! Build densified grid inside microzone centered around iq_ibz
    do i1=-ndivs(1),ndivs(1),1
      qpt(1) = q0(1) + real(i1) / (ndivs(1) * sigma%ngqpt(1))
      do i2=-ndivs(2),ndivs(2),1
        qpt(2) = q0(2) + real(i2) / (ndivs(2) * sigma%ngqpt(2))
        do i3=-ndivs(3),ndivs(3),1
          qpt(3) = q0(3) + real(i2) / (ndivs(3) * sigma%ngqpt(3))

          qpt_cart = two * pi * matmul(cryst%gprimd, qpt)
          qmod = sqrt(sum(qpt_cart ** 2))
          fqdamp = exp((-qmod/sigma%qdamp) ** 2)
          ! TODO: In principle one should introduce weigths for points that are close to the border
          ! Should write routine to compute weights given kptrlatt, kibz and radius.
          if (dot_product(qpt_cart, qpt_cart) <= qrad2) cycle
          vol_fact = one
          call ifc%fourq(cryst, qpt, phfrq, displ_cart)
          inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))

          ! Interpolate e_{n k+q}
          kq = kk + qpt
          do ib_k=1,nbcalc_ks
            band_ks = ib_k + bstart_ks - 1
            call sigma%frohl_skw%eval_bks(band_ks, kq, spin,  eigs_kq(ib_k))
          end do

          ! Sum over modes for this q-point.
          do nu=1,natom3
            wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
            ! Compute gkq_{LR} without (i 4pi/ucvol)
            ! note that at this level it does not dependend on ib_k.
            cnum = zero
            do iatom=1,cryst%natom
              cdd = cmplx(displ_cart(1,:, iatom, nu), displ_cart(2,:, iatom, nu), kind=dpc) * &
                    exp(-j_dpc * two_pi * dot_product(qpt, cryst%xred(:, iatom)))
              cnum = cnum + dot_product(qpt_cart, matmul(ifc%zeff(:, :, iatom), cdd))
            end do
            gkq2 = (real(cnum) ** 2 + aimag(cnum) ** 2) / (two * wqnu) * inv_qepsq * weight_q * vol_fact !* fqdamp ** 2

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
#endif

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
    ! Precompute denominator
    inv_qepsq = one / dot_product(sigma%qvers_cart(:, iang), matmul(ifc%dielt, sigma%qvers_cart(:, iang)))

    ! Get all omega_{q nu} and e_{n k+q} along the q-line for this (theta, phi)
    do iqr=1,sigma%nqr
      qpt_cart = sigma%qvers_cart(:, iang) * qrmesh(iqr)
      qmod = sqrt(sum(qpt_cart ** 2))
      qpt = matmul(cryst%gprimd, qpt_cart)
      !if (iqr == 1) cycle

      if (iqr == 1) then
        ! Include non-analytical part for the first point
        call ifc%fourq(cryst, qpt_cart, phfrq, displ_cart, nanaqdir="cart")
      else
        call ifc%fourq(cryst, qpt, phfrq, displ_cart)
      end if
      wqr(iqr, :) = phfrq

      !inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))
      fqdamp = (two / cryst%ucvol) ** 2 * inv_qepsq ** 2 !* exp(-(qmod/sigma%qdamp) ** 2)

      ! Compute gkq_{LR}. Note that in our approx it does not dependend on ib_k.
      gkqr2(iqr, :) = zero
      do nu=1,natom3
        wqnu = phfrq(nu); if (wqnu < EPH_WTOL) cycle
        cnum = zero
        do iatom=1,cryst%natom
          ! This is complex
          cdd = cmplx(displ_cart(1,:, iatom, nu), displ_cart(2,:, iatom, nu), kind=dpc) * &
                exp(-j_dpc * two_pi * dot_product(qpt, cryst%xred(:, iatom)))
          cnum = cnum + dot_product(qpt_cart, matmul(ifc%zeff(:, :, iatom), cdd))
        end do
        gkqr2(iqr, nu) = (real(cnum) ** 2 + aimag(cnum) ** 2) / (two * wqnu)
      end do
      gkqr2 = gkqr2 * fqdamp

      ! Interpolate e_{n k+q} and store results.
      kq = kk + qpt
      do ib_k=1,nbcalc_ks
        band_ks = ib_k + bstart_ks - 1
        call sigma%frohl_skw%eval_bks(band_ks, kq, spin, eigs_kqr(iqr, ib_k))
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
          !write(666, *)"#iang, nu, ib_k, it, cfqr(1:nqr)"
          !write(666, *)iang, nu, ib_k, it, cfqr
          cfact = two_pi * simpson_cplx(sigma%nqr, qstep, cfqr) * sigma%angwgth(iang)

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

            sigma%frohl_dvals_de0ks(it, ib_k) = sigma%frohl_dvals_de0ks(it, ib_k) + &
                two_pi * simpson(qstep, rfactr) * sigma%angwgth(iang)

            ! Accumulate Sigma(w) for state ib_k if spectral function is wanted.
            ! This is gonna be costly ...
            ! TODO: weigths
            if (sigma%nwr > 0) then
              do iw=1,sigma%nwr
                cfqr(:) = gkqr2(:, nu) * ( &
                  (nqr(:) + f_mkqr(:)      ) / (sigma%wrmesh_b(iw, ib_k) - eig0mkqr(:) + wqr(:, nu) + sigma%ieta) + &
                  (nqr(:) - f_mkqr(:) + one) / (sigma%wrmesh_b(iw, ib_k) - eig0mkqr(:) - wqr(:, nu) + sigma%ieta) )
                cfact = two_pi * simpson_cplx(sigma%nqr, qstep, cfqr)
                sigma%frohl_vals_wr(iw, it, ib_k) = sigma%frohl_vals_wr(iw, it, ib_k) * cfact * sigma%angwgth(iang)
              end do
            end if
          end if

        end do ! it
      end do ! ib_k

    end do ! iqr
 end do ! iang

 ABI_FREE(displ_cart)
 ABI_FREE(eigs_kq)
 ABI_FREE(eigs_kqr)
 ABI_FREE(wqr)
 ABI_FREE(eig0mkqr)
 ABI_FREE(f_mkqr)
 ABI_FREE(gkqr2)
 ABI_SFREE(cfact_wr)

 ! Accumulate final results
 call xmpi_sum(sigma%frohl_vals_e0ks, comm, ierr)
 call xmpi_sum(sigma%frohl_dvals_de0ks, comm, ierr)
 if (sigma%nwr > 0) call xmpi_sum(sigma%frohl_vals_wr, comm, ierr)

 call cwtime_report(" Frohlich self-energy", cpu_fr, wall_fr, gflops_fr)

end subroutine eval_sigfrohl2
!!***

!!****f* m_sigmaph/qpoints_oracle
!! NAME
!!  qpoints_oracle
!!
!! FUNCTION
!!  This function tries to predict the **full** list of q-points needed to compute the lifetimes
!!  once we know sigma%nkcalc. It uses an energy window computed from the max phonon frequency
!!  multiplied by sigma%winfact.
!!
!! INPUT
!! cryst=Crystalline structure
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! qpts(3, nqpt)=
!! nqpt= Number of points in qpts
!! nqbz=Number of q-points in BZ.
!! qbz(3, nbz) = full BZ
!! comm=MPI communicator.
!!
!! OUTPUT
!!  qselect(nqpt)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine qpoints_oracle(sigma, cryst, ebands, qpts, nqpt, nqbz, qbz, qselect, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpt, nqbz, comm
 class(sigmaph_t),intent(in) :: sigma
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
!arrays
 real(dp),intent(in) :: qpts(3,nqpt), qbz(3,nqbz)
 integer,intent(out) :: qselect(nqpt)

!Local variables ------------------------------
!scalars
 integer,parameter :: timrev1=1
 integer :: spin, ikcalc, ik_ibz, iq_bz, ierr, db_iqpt, ibsum_kq, ikq_ibz, ikq_bz
 integer :: cnt, my_rank, nprocs, ib_k, band_ks, nkibz, nkbz, kq_rank
 real(dp) :: eig0nk, eig0mkq, dksqmax, ediff, cpu, wall, gflops
 character(len=500) :: msg
 type(krank_t) :: krank
!arrays
 integer :: g0(3)
 integer,allocatable :: qbz_count(:), qbz2qpt(:,:), bz2ibz(:,:)
 real(dp) :: kq(3), kk(3)
 real(dp),allocatable :: wtk(:), kibz(:,:), kbz(:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call cwtime(cpu, wall, gflops, "start")
 call wrtout(std_out, sjoin(" qpoints_oracle: predicting number q-points for tau with winfact:", ftoa(sigma%winfact)))

 ! Get full BZ associated to ebands
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
   nkibz, kibz, wtk, nkbz, kbz, bz2ibz=bz2ibz)
 call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 ABI_FREE(wtk)
 ABI_FREE(kibz)
 ABI_CHECK(nkibz == ebands%nkpt, "nkibz != ebands%nkpt")

 ! Build BZ --> IBZ mapping using symrec.
 !ABI_MALLOC(bz2ibz, (nkbz, 6))
 !call listkk(dksqmax, cryst%gmet, bz2ibz, ebands%kptns, kbz, ebands%nkpt, nkbz, cryst%nsym, &
 !     1, cryst%symafm, cryst%symrec, sigma%timrev, comm, exit_loop=.True., use_symrec=.True.)
 !if (dksqmax > tol12) then
 !  write(msg, '(a, es16.6)' ) &
 !   "At least one of the points in BZ could not be generated from a symmetrical one. dksqmax: ", dksqmax
 !  MSG_ERROR(msg)
 !end if
 !call cwtime_report(" qpoints_oracle_listkk1", cpu, wall, gflops)

 ! Make full k-point rank arrays
 krank = krank_new(nkbz, kbz)

 ABI_ICALLOC(qbz_count, (nqbz))
 cnt = 0
 do spin=1,sigma%nsppol
   do ikcalc=1,sigma%nkcalc
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism
     kk = sigma%kcalc(:, ikcalc)
     ik_ibz = sigma%kcalc2ibz(ikcalc, 1)
     do iq_bz=1,nqbz
       if (qbz_count(iq_bz) /= 0) cycle ! No need to check this q-point again.
       kq = kk + qbz(:, iq_bz)
       kq_rank = krank%get_rank(kq)
       ikq_bz = krank%invrank(kq_rank)
       ABI_CHECK(ikq_bz > 0, sjoin("Cannot find kq: ", ktoa(kq)))
       ABI_CHECK(isamek(kq, kbz(:, ikq_bz), g0), "Wrong invrank")
       !ikq_ibz = bz2ibz(ikq_bz,1)
       ikq_ibz = bz2ibz(1, ikq_bz)
       do ib_k=1,sigma%nbcalc_ks(ikcalc, spin)
         band_ks = ib_k + sigma%bstart_ks(ikcalc, spin) - 1
         eig0nk = ebands%eig(band_ks, ik_ibz, spin)
         do ibsum_kq=sigma%bsum_start, sigma%bsum_stop
           eig0mkq = ebands%eig(ibsum_kq, ikq_ibz, spin)
           ediff = eig0nk - eig0mkq
           ! Perform check on the energy difference to exclude this q-point.
           if (abs(ediff) <= sigma%winfact * sigma%wmax) qbz_count(iq_bz) = qbz_count(iq_bz) + 1
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(kbz)
 call krank%free()

 call xmpi_sum(qbz_count, comm, ierr)
 call cwtime_report(" qbz_count", cpu, wall, gflops)

 ! Get mapping QBZ --> QPTS q-points.
 ABI_MALLOC(qbz2qpt, (nqbz, 6))
 call listkk(dksqmax, cryst%gmet, qbz2qpt, qpts, qbz, nqpt, nqbz, cryst%nsym, &
      1, cryst%symafm, cryst%symrec, timrev1, comm, exit_loop=.True., use_symrec=.True.)
 if (dksqmax > tol12) then
   write(msg, '(a,es16.6,2a)' )&
     "At least one of the q points could not be generated from a symmetrical one in the DVDB. dksqmax: ",dksqmax, ch10, &
     "Action: check your DVDB file and use eph_task to interpolate the potentials on a denser q-mesh."
   MSG_ERROR(msg)
 end if
 call cwtime_report(" oracle_listkk_qbz_qpts", cpu, wall, gflops)

 ! Compute qselect using qbz2qpt.
 qselect = 0
 do iq_bz=1,nqbz
   if (qbz_count(iq_bz) == 0) cycle
   db_iqpt = qbz2qpt(iq_bz, 1)
   qselect(db_iqpt) = qselect(db_iqpt) + 1
 end do

 ABI_FREE(qbz_count)
 ABI_FREE(qbz2qpt)
 ABI_FREE(bz2ibz)

 cnt = count(qselect /= 0)
 write(std_out, "(a, i0, a, f5.1, a)")" qpoints_oracle: calculation of tau_nk will need: ", cnt, &
   " q-points in the IBZ. (nqibz_eff / nqibz): ", (100.0_dp * cnt) / sigma%nqibz, " [%]"
 !qselect = 0

end subroutine qpoints_oracle
!!***

end module m_sigmaph
!!***

!!****m* ABINIT/m_gwr
!! NAME
!!  m_gwr
!!
!! FUNCTION
!!  Objects and procedures to implement GW in real-space and imaginary time.
!!
!! NOTES
!!
!! Memory and workload are distributed using a 4D cartesian grid: (g/r, tau, k-points, spin).
!!
!! Inside the g/r communicator, we use PBLAS matrices to store G, tchi and W.
!! using a 1D processor grid with block distribution either along columns or rows.
!! A 2D grid would require MPI-FFT.
!!
!! Let's assume for simplicity that we have only two MPI procs in the g/r communicator.
!! Matrices in (g,g') space are distributed along columns so that the g-index is local
!! and we can use sequential zero-padded FFTs to transform from g to r in the unit cell:
!!
!!                     g'-axis
!!              |--------------------
!!              |         |         |
!!     g-axis   |   P0    |   P1    |
!!              |         |         |
!!              |--------------------
!!
!! The results of the FFT transform along g are stored in another PBLAS matrix with the same layout:
!!
!!                     g'-axis
!!              |--------------------
!!              |         |         |
!!     r-axis   |   P0    |   P1    |
!!              |         |         |
!!              |--------------------
!!
!! At this point, we use ptrans to MPI transpose the (r, g') matrix + complex conjugation (if needed)
!! and we end up with:
!!
!!                     r-axis
!!              |--------------------
!!              |         |         |
!!     g'-axis  |   P0    |   P1    |
!!              |         |         |
!!              |--------------------
!!
!! Differences with respect to the GW code in frequency-domain:
!!
!!  - in GWR, the k-mesh for G must be Gamma-centered.
!!  - All the two-point functions are defined on k/q-centered g-spheres while GW uses a single Gamma-centered sphere.
!!  - At present only one-shot GW is supported.
!!  - It is not clear if it makes sense to support all the options of the GW code, e.g. COHSEX, SEX, etc.
!!  - No plasmopole model in GWR. In principle it is possible but where is the point?
!!  - The frequency/tau meshes are automatically defined by ntau and the KS spectrum (minimax meshes)
!!
!! Technical properties:
!!
!!   - Integration of vcoul_t is not easy (q-centered gvec vs single sphere, memory is not MPI distributed)
!!     Solution: extract reusable components from vcoul_t that can be called inside the loop over q in IBZ.
!!     Also it's not clear to me that one can use vc(Sq, SG) when a cutoff is used as the cutoff breaks
!!     the spherical symmetri of vc(r). Besides, when symmetries are used to reconstruct the term for q in the BZ,
!!     one might have to take into account umklapps. Use cache?
!!
!!   - Computation of Sigma_x = iGv must be done in Fourier space using the Lehmann representation so
!!     that we can handle the long-range behavior in q-space. Unfortunately, we cannot simply call calc_sigx_me
!!     from GWR so we have to implement a new routine.
!!
!!   - Treatment of the anisotropic behaviour of Wc. This part is badly coded in GW, in the sense that
!!     we use a finite small q when computing Wc for q --> 0. This breaks the symmetry of the system
!!     and QP degeneracies. The equations needed to express the angular dependency of W(q) for q --> 0
!!     are well known but one has to pass through the Adler-Wiser expression.
!!     Solution: Compute heads and wings using a WFK_fine wavefunction file with dense k-mesh and less bands.
!!     The dipole matrix elements are computed with the DFPT routines, still we need to
!!     recode a lot of stuff that is already done in cchi0q0, especially symmetries.
!!     Note, however, that tchi is Hermitian along the imaginary axis, expect for omega = 0 in metals
!!     but I don't think the minmax grids contain omega = 0.
!!
!!  - In principle, it's possible to compute QP correction along a k-path is a new WFK file is provided.
!!    The correlated part is evaluated in real-space in the super-cell.
!!    For Sigma_x, we need a specialized routine that can handle arbitrary q, especially at the level of v(q, G)
!!    but I don't know if this approach will give smooth bands
!!    as we don't have q --> 0 when k does not belong to the k-mesh.
!!
!!  - New routine for direct diagonalization of the KS Hamiltonian based on PBLAS?
!!
!!  - New routine to compute oscillator matrix elements with NC/PAW and PBLAS matrices.
!!    It can be used to compute tchi head/wings as well as Sigma_x + interface with coupled-cluster codes.
!!
!!  - Decide whether we should use VASP conventions for G and the analytic continuation
!!    or the "standard" ones by Godby.
!!    The standard ones are consistent with Hedin's notations and correspond to the ones used in the legacy GW code.
!!    On the other hand, VASP notations make life easier if one has to implement PAW.
!!
!!  - Address nspinor = 2 and PBLAS distribution as MPI proc can have both spinors in memory
!!    In other words, we should store the first/last index in gvec for each spinor
!!
!!  - Optimization for Gamma-only. Memory and c -> r FFTs
!!
!! COPYRIGHT
!! Copyright (C) 1999-2021 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwr

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_hdr
 use m_ebands
 use netcdf
 use m_nctk
 use m_dtfil
 use m_yaml
 use m_sigtk
 use iso_c_binding

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use defs_abitypes,   only : mpi_type
 use m_gwdefs,        only : GW_TOLQ0, GW_Q0_DEFAULT
 use m_time,          only : cwtime, cwtime_report, sec2str
 use m_io_tools,      only : iomode_from_fname, get_unit !, file_exists, open_file
 use m_numeric_tools, only : get_diag, isdiagmat, arth, print_arr
 use m_copy,          only : alloc_copy
 use m_geometry,      only : normv
 use m_fstrings,      only : sjoin, itoa, strcat, ktoa, ltoa, ftoa
 use m_krank,         only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_crystal,       only : crystal_t
 use m_dtset,         only : dataset_type
 use m_fftcore,       only : get_kg, sphereboundary, ngfft_seq, getng, print_ngfft !, kgindex
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_distribfft,    only : init_distribfft_seq
 use m_fft,           only : fft_ug, fft_ur, fftbox_plan3_t, fourdp
 use m_fft_mesh,      only : calc_ceikr !, times_eikr !, times_eigr, ig2gfft, get_gftt, calc_eigr
 use m_kpts,          only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_map_print, kpts_pack_in_stars
 use m_melemts,       only : melements_t
 use m_ioarr,         only : read_rhor
 use m_slk,           only : matrix_scalapack, processor_scalapack, slk_array_free, slk_array_locmem_mb
 use m_wfk,           only : wfk_read_ebands, wfk_t, wfk_open_read
 use m_wfd,           only : wfd_init, wfd_t, wfdgw_t, wave_t, WFD_STORED
 use m_pawtab,        only : pawtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, &
                             pawrhoij_inquire_dim, pawrhoij_symrhoij, pawrhoij_unpack


!#define __HAVE_GREENX
#undef __HAVE_GREENX

#ifdef __HAVE_GREENX
 use mp2_grids,      only : gx_minimax_grid
#endif

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_gwr/desc_t
!! NAME
!! desc_t
!!
!! FUNCTION
!!  Parameters related to a two-point function such as
!!  gvectors, tables used for zero padded FFTs and matrix elements of the coulomb interaction.
!!
!! SOURCE

 type,public :: desc_t

   integer :: istwfk = 1
   ! Storage mode for this k point.

   integer :: npw = -1
   ! Total number of plane-waves for this k/q-point.

   integer,allocatable :: gvec(:,:)
   ! gvec(3, npw)
   ! G vectors in reduced coordinates.
   ! Note that this array is global i.e. it is not MPI-distributed in the PBLAS communicator.

   integer,allocatable :: gbound(:,:)
   ! gbound(2*mgfft+8, 2)
   ! sphere boundary info for zero-padded FFT

   complex(gwpc),allocatable :: vc_sqrt(:)
   ! (npw)
   ! Square root of the Coulomb interaction in reciprocal space.
   ! Allocated and computed for tchi/W descriptors
   ! A cutoff might be applied

   integer, allocatable :: gvec2gsort(:)
   ! (npw)
   ! Mapping the gvec array and the sorted one.
   ! Computed by calling desc_calc_gnorm_table.

   real(dp), allocatable :: sorted_gnorm(:)
   ! (npw)
   ! Sorted list with the norm of the gvector
   ! Ccomputed by calling desc_calc_gnorm_table.

 contains

   procedure :: copy => desc_copy
   ! Copy object.

   procedure :: free => desc_free
   ! Free memory.

   procedure :: calc_gnorm_table => desc_calc_gnorm_table
   ! Compute mapping used to loop over G-vectors ordered by norm.

 end type desc_t

 interface desc_array_free
   module procedure desc_array1_free
 end interface desc_array_free
!!***

!----------------------------------------------------------------------

!!****t* m_gwr/gwr_t
!! NAME
!! gwr_t
!!
!! FUNCTION
!!  This object provides the API high-level to perform the different steps of the GWR algorithm.
!!
!! SOURCE

 type, public :: gwr_t

   integer :: nsppol = 1, nspinor = -1, nspden = -1
   ! Number of independent spin polarizations, number of spinor components and spin density.

   integer :: natom = -1
    ! Number of atoms

   integer :: usepaw = -1

   integer :: my_nspins = -1
   ! Number of independent spin polarizations treated by this MPI proc

   integer :: nkbz = -1, nkibz = -1
   ! Number of k-points in the BZ/IBZ

   integer :: my_nkibz = -1, my_nkbz = -1
   ! Number of k-points in the IBZ/BZ stored by this MPI proc.

   integer,allocatable :: my_kbz_inds(:)
   ! (my_nkbz)
   ! List of k-BZ indices treated by this proc.

   integer,allocatable :: my_kibz_inds(:)
   ! (my_nkibz)
   ! List of k-IBZ indices treated by this proc.

   integer :: nqbz = -1, nqibz = -1
   ! Number of q-points in the BZ/IBZ

   integer :: my_nqibz = -1, my_nqbz = -1
   ! Number of q-points in the IBZ/BZ stored by this MPI proc.

   integer,allocatable :: my_qibz_inds(:)
   ! (my_nqibz)
   ! List of q-IBZ indices treated by this proc.

   integer,allocatable :: my_qbz_inds(:)
   ! (my_nqbz)
   ! List of q-IBZ indices treated by this proc.

   integer :: ntau = -1
   ! Total number of imaginary time points.

   integer :: my_ntau = -1
   ! Number of imaginary time/frequency points treated by this MPI rank.

   !integer :: nomega = -1, my_nomega = -1

  integer :: nkcalc
   ! Number of Sigma_nk k-points computed
   ! TODO: Should be spin dependent + max_nkcalc

  integer :: max_nbcalc
   ! Maximum number of bands computed (max over nkcalc and spin).

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.

  integer,allocatable :: bstart_ks(:,:)
   ! bstart_ks(nkcalc, nsppol)
   ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all degenerate states should be included when symmetries are used.

  integer,allocatable :: bstop_ks(:,:)
   ! bstop_ks(nkcalc, nsppol)

  integer,allocatable :: nbcalc_ks(:,:)
   ! nbcalc_ks(nkcalc, nsppol)
   ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all degenerate states should be included when symmetries are used.

   integer,allocatable :: kcalc2ibz(:,:)
   !kcalc2ibz(nkcalc, 6))
   ! Mapping ikcalc --> IBZ as reported by listkk.

   logical :: use_supercell = .True.
   ! True if we are using the supercell formalism.
   ! False if we are using the mixed-space approach with convolutions in k-space.

   integer :: ngkpt(3) = -1
   ! Number of divisions in k-point mesh.

   integer :: ngqpt(3) = -1
   ! Number of divisions in q-point mesh.

   integer,allocatable :: my_spins(:)
   ! (my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2.

   integer,allocatable :: my_itaus(:)
   ! (my_ntau)
   ! Indirect table giving the tau indices treated by this MPI rank.

   integer,allocatable :: tau_master(:)
   ! (ntau)
   ! The rank of the MPI proc in tau_comm treating itau.

   integer, allocatable :: np_qibz(:)
   ! (nqibz)
   ! Number of processors in kpt_comm treating iq_ibz

   integer, allocatable :: np_kibz(:)
   ! (nkibz)
   ! Number of processors in kpt_comm treating iq_ibz

   real(dp),allocatable :: tau_mesh(:), tau_wgs(:)
   ! (ntau)
   ! Imaginary tau mesh and integration weights.

   real(dp),allocatable :: iw_mesh(:), iw_wgs(:)
   ! (ntau)
   ! Imaginary frequency mesh and integration weights

   real(dp),allocatable :: t2w_cos_wgs(:,:)
   ! (ntau, ntau)
   ! weights for cosine transform. (i tau --> i omega)

   real(dp),allocatable :: w2t_cos_wgs(:,:)
   ! (ntau, ntau)
   ! weights for sine transform (i iomega --> i tau)

   real(dp),allocatable :: t2w_sin_wgs(:,:)
   ! (ntau, ntau)
   ! weights for sine transform (i tau --> i omega)

   real(dp) :: ft_max_error(3) = -one
   ! Max error due to inhomogenous FT.

   real(dp) :: cosft_duality_err = -one
   ! Max_{ij} |CT CT^{-1} - I|

   integer :: green_mpw = -1
   ! Max number of g-vectors for Green's function over k-points.

   integer :: tchi_mpw = -1
   ! Max number of g-vectors for tchi over q-points.

   !integer :: sigma_mpw = -1
   ! Max number of g-vectors for Sigma over q-points.

   integer :: g_ngfft(18) = -1, g_mgfft = -1, g_nfft = -1
   ! FFT mesh for the Green's function.

   !integer :: chi_ngfft(18) = -1, chi_mgfft = -1, chi_nfft = -1
   !integer :: sig_ngfft(18) = -1, sig_mgfft = -1, sig_nfft = -1

   type(desc_t),allocatable :: green_desc_kibz(:)
   ! (nkibz)
   ! Descriptor for Green's functions

   type(desc_t),allocatable :: tchi_desc_qibz(:)
   ! (nqibz)
   ! Descriptor for tchi

   !type(desc_t),allocatable :: sigma_desc_kibz(:)
   ! (nkibz)
   ! Descriptor for self-energy

   integer :: coords_gtks(4) = 0
   ! Cartesian coordinates of this processor in the Cartesian grid.

   integer :: comm
   ! Initial MPI communicator with all procs involved in the calculation

   type(xcomm_t) :: spin_comm
   ! MPI communicator over spins.

   type(xcomm_t) :: kpt_comm
   ! MPI communicator for k/q-point distribution.

   type(xcomm_t) :: g_comm
   ! MPI communicator for g/r distribution

   type(xcomm_t) :: tau_comm
   ! MPI communicator for imag time distribution

   type(xcomm_t) :: gtau_comm
   ! MPI communicator for g/tau 2D subgrid.

   type(xcomm_t) :: tks_comm
   ! MPI communicator for tau/kpoint/spin 3D grid

   type(xcomm_t) :: gtk_comm
   ! MPI communicator for /gtau/kpoint 3D grid

   type(dataset_type), pointer :: dtset => null()
   ! Input variables.

   type(datafiles_type), pointer :: dtfil => null()
   ! Names of input/output files and prefixes.

   type(crystal_t), pointer  :: cryst => null()
   ! Crystal structure.

   type(ebands_t), pointer  :: ks_ebands => null()
   ! KS energies

   type(ebands_t) :: qp_ebands
   ! QP energies

   type(pseudopotential_type), pointer :: psps => null()
   ! NC Pseudos data

   type(pawtab_type), pointer :: pawtab(:) => null()
   ! PAW data

   type(mpi_type),pointer :: mpi_enreg => null()

   type(processor_scalapack) :: g_slkproc
   !type(processor_scalapack) :: gtau_slkproc

   type(matrix_scalapack),allocatable :: gt_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)
   ! Occupied/Empty Green's function G_k(g,g')

   type(matrix_scalapack),allocatable :: tchi_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! Irreducible polarizability tchi_q(g,g')

   character(len=10) :: tchi_space = "none"
   ! "none", "itau", "iomega"

   type(matrix_scalapack),allocatable :: wc_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! Correlated screened Coulomb interaction (does not depend on the spin, though)

   character(len=10) :: wc_space = "none"
   ! "none", "itau", "iomega"

   type(matrix_scalapack),allocatable :: sigc_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)

   !character(len=10) :: sigc_space = "none"
   ! "none", "itau", "iomega"

   !character(len=fnlen) :: wfk_path = ABI_NOFILE
   ! Path to the WFK file with the KS wavefunctions.

   character(len=fnlen) :: gwrnc_path = ABI_NOFILE
   ! Path to the GWR.nc file with output results.

   real(dp),allocatable :: kbz(:,:)
   ! (3, nkbz)
   ! Reduced coordinates of the k-points in the full BZ.

   real(dp), contiguous, pointer :: kibz(:,:) => null()
    ! (3, nkibz)
    ! Reduced coordinates of the k-points in the IBZ

   integer,allocatable :: my_kbz2ibz(:,:)
    ! (6, my_nkbz)
    ! These table used the conventions for the symmetrization of the wavefunctions expected by cgtk_rotate.
    ! In this case listkk has been called with symrel and use_symrec=False

   integer,allocatable :: my_qbz2ibz(:,:)
    ! (6, my_nqbz)

   real(dp), contiguous, pointer :: wtk(:) => null()
    ! (nkibz)
    ! Weights of the k-points in the IBZ (normalized to one).

   real(dp),allocatable :: qbz(:,:)
    ! (3, nqbz)
    ! Reduced coordinates of the q-points in the full BZ.

   integer,allocatable :: qbz2ibz(:,:)
   ! (6, nqbz)
   ! Mapping qBZ to IBZ.

   real(dp),allocatable :: qibz(:,:)
   ! (3, nqibz)
   ! Reduced coordinates of the q-points in the IBZ (full simmetry of the system).

   real(dp),allocatable :: wtq(:)
   ! (nqibz)
   ! Weights of the q-points in the IBZ (normalized to one).

   type(wfdgw_t) :: kcalc_wfd
   ! wavefunction descriptors with the KS states where QP corrections are wanted.

   type(melements_t) :: ks_me !, qp_me
   ! Matrix elements of the different potentials in the KS basis set.

   type(degtab_t),allocatable :: degtab(:,:)
   ! (nkcalc, nsppol)
   ! Table used to average QP results in the degenerate subspace if symsigma == 1

 contains

   procedure :: init => gwr_init
   ! Initialize the object.

   procedure :: rotate_gt => gwr_rotate_gt
   ! Reconstruct the Green's functions in the kBZ from the IBZ.

   procedure :: get_green_gpr => gwr_get_green_gpr
    ! G_k(g,g') --> G_k(g',r) for each k in the BZ treated by this MPI proc for given spin and tau.

   procedure :: rotate_wct => gwr_rotate_wct
   ! Reconstruct Wc(q) in the BZ from the IBZ.

   procedure :: get_wc_gpr => gwr_get_wc_gpr
   ! W_q(g,g') --> W_q(g',r) for each q in the BZ treated by this MPI procs for given spin and tau.

   procedure :: cos_transform  => gwr_cos_transform
   ! Inhomogeneous cosine transform.

   procedure :: free => gwr_free
   ! Free memory.

   procedure :: print => gwr_print
   ! Print info on the object.

   procedure :: print_trace => gwr_print_trace
    ! Print trace of matrices for testing purposes.

   procedure :: load_kcalc_wfd => gwr_load_kcalc_wfd
   !  Load the KS states for Sigma_nk from the WFK file

   procedure :: build_gtau_from_wfk => gwr_build_gtau_from_wfk
   ! Build the Green's function in imaginary time for k-points in the IBZ from the WFK file

   procedure :: build_tchi => gwr_build_tchi
   ! Build the irreducible polarizability

   procedure :: build_wc => gwr_build_wc
   ! Build the correlated part of the screened interaction.

   procedure :: build_sigmac => gwr_build_sigmac
   ! Build the correlated part of the self-energy iGWc
   ! and compute matrix elements in the KS representation.

   procedure :: rpa_energy => gwr_rpa_energy
   ! Compute RPA energy.

   procedure :: run_g0w0 => gwr_run_g0w0
   ! Compute QP corrections with G0W0.

   procedure :: ncwrite_tchi_wc => gwr_ncwrite_tchi_wc
   ! Write tchi or wc to netcdf file

   !procedure :: ncread_tchi_wc => gwr_ncread_tchi_wc
   ! Read tchi or wc to netcdf file

 end type gwr_t
!!***

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_init
!! NAME
!! gwr_init
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_init(gwr, dtset, dtfil, cryst, psps, pawtab, ks_ebands, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 class(gwr_t),target,intent(out) :: gwr
 type(dataset_type),target,intent(in) :: dtset
 type(datafiles_type),target,intent(in) :: dtfil
 type(crystal_t),target,intent(in) :: cryst
 type(pseudopotential_type),target,intent(in) :: psps
 type(pawtab_type),target,intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(ebands_t),target,intent(in) :: ks_ebands
 type(mpi_type),target,intent(in) :: mpi_enreg
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: me_fft0 = 0, paral_fft0 = 0, nproc_fft1 = 1, istwfk1 = 1
 integer,parameter :: qptopt1 = 1, timrev1 = 1, master = 0, ndims = 4
 integer :: my_is, my_it, my_ikf, my_iqf, ii, npwsp, col_bsize, ebands_timrev, my_iki, my_iqi, itau, spin
 integer :: my_nshiftq, iq_bz, iq_ibz, npw_, ncid, ig, ig_start
 integer :: comm_cart, me_cart, ierr, all_nproc, nps, my_rank, qprange_, gap_err
 integer :: jj, cnt, ikcalc, ndeg, mband, bstop, nbsum, it, iw
 integer :: ik_ibz, ik_bz, isym_k, trev_k, g0_k(3)
 logical :: isirr_k, changed, q_is_gamma
 real(dp) :: ecut_eff, cpu, wall, gflops, mem_mb, te_min, te_max
 character(len=5000) :: msg
 logical :: reorder
 type(krank_t) :: qrank, krank_ibz
 type(gaps_t) :: gaps
!arrays
 integer :: qptrlatt(3,3), dims(ndims)
 integer :: indkk_k(6,1)
 integer,allocatable :: kbz2ibz(:,:), qbz2ibz(:,:), gvec_(:,:)
 integer,allocatable :: degblock(:,:), degblock_all(:,:,:,:), ndeg_all(:,:)
 real(dp) :: my_shiftq(3,1), kk_ibz(3), kk_bz(3), qq_bz(3), qq_ibz(3), kk(3)
 real(dp),allocatable :: wtk(:), kibz(:,:), mat(:,:)
 logical :: periods(ndims), keepdim(ndims)

! *************************************************************************

 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 call cwtime(cpu, wall, gflops, "start")

 gwr%dtset => dtset
 gwr%dtfil => dtfil
 gwr%cryst => cryst
 gwr%psps => psps
 gwr%pawtab => pawtab
 gwr%ks_ebands => ks_ebands
 ! TODO Make sure fermie is inside the gap if semiconductor.
 ! perhaps this shoud be delegated to gwr_driver
 gwr%kibz => ks_ebands%kptns
 gwr%wtk => ks_ebands%wtk
 gwr%mpi_enreg => mpi_enreg

 ! Initialize qp_ebands with KS values.
 call ebands_copy(ks_ebands, gwr%qp_ebands)

 gwr%comm = comm
 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%nspden = dtset%nspden
 gwr%natom = dtset%natom
 gwr%usepaw = dtset%usepaw
 gwr%use_supercell = .True.

 mband = ks_ebands%mband
 nbsum = dtset%nband(1)
 ABI_CHECK_IRANGE(nbsum, 1, mband, "Invalid nbsum")

 ! =======================
 ! Setup k-mesh and q-mesh
 ! =======================

 ! Get full kBZ associated to ks_ebands
 call kpts_ibz_from_kptrlatt(cryst, ks_ebands%kptrlatt, ks_ebands%kptopt, ks_ebands%nshiftk, ks_ebands%shiftk, &
                             gwr%nkibz, kibz, wtk, gwr%nkbz, gwr%kbz) !, bz2ibz=bz2ibz)
                             !new_kptrlatt=gwr%kptrlatt, new_shiftk=gwr%kshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ABI_FREE(wtk)

 ! In principle kibz should be equal to ks_ebands%kptns
 ABI_CHECK(gwr%nkibz == ks_ebands%nkpt, "nkibz != ks_ebands%nkpt")
 ABI_CHECK(all(abs(ks_ebands%kptns - kibz) < tol12), "ks_ebands%kibz != kibz")

 if (.not. (isdiagmat(ks_ebands%kptrlatt) .and. ks_ebands%nshiftk == 1)) then
   ABI_ERROR("GWR requires ngkpt with one shift!")
 end if
 gwr%ngkpt = get_diag(ks_ebands%kptrlatt)

 ! Note symrec convention
 ABI_MALLOC(kbz2ibz, (6, gwr%nkbz))
 ebands_timrev = kpts_timrev_from_kptopt(ks_ebands%kptopt)

 krank_ibz = krank_from_kptrlatt(gwr%nkibz, kibz, ks_ebands%kptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", ebands_timrev, cryst, krank_ibz, gwr%nkbz, gwr%kbz, kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 ! Order kbz by stars and rearrange entries in kbz2ibz table.
 call kpts_pack_in_stars(gwr%nkbz, gwr%kbz, kbz2ibz)

 if (my_rank == master) then
   call kpts_map_print([std_out, ab_out], "Mapping kBZ --> kIBZ", "symrec", gwr%kbz, kibz, kbz2ibz, gwr%dtset%prtvol)
 end if

 !call get_ibz2bz(gwr%nkibz, gwr%nkbz, kbz2ibz, kibz2bz, ierr)
 !ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for k-points")

 ! Setup qIBZ, weights and BZ.
 ! Always use q --> -q symmetry even in systems without inversion
 ! TODO: Might add input variable to rescale the q-mesh.
 my_nshiftq = 1; my_shiftq = zero
 qptrlatt = ks_ebands%kptrlatt
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
                             gwr%nqibz, gwr%qibz, gwr%wtq, gwr%nqbz, gwr%qbz)
                             !new_kptrlatt=gwr%qptrlatt, new_shiftk=gwr%qshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 gwr%ngqpt = get_diag(qptrlatt)

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(qbz2ibz, (6, gwr%nqbz))

 qrank = krank_from_kptrlatt(gwr%nqibz, gwr%qibz, qptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", timrev1, cryst, qrank, gwr%nqbz, gwr%qbz, qbz2ibz) /= 0) then
   ABI_ERROR("Cannot map qBZ to IBZ!")
 end if
 call qrank%free()

 ! Order qbz by stars and rearrange entries in qbz2ibz table.
 call kpts_pack_in_stars(gwr%nqbz, gwr%qbz, qbz2ibz)
 if (my_rank == master) then
   call kpts_map_print([std_out, ab_out], "Mapping qBZ --> qIBZ", "symrec", gwr%qbz, gwr%qibz, qbz2ibz, gwr%dtset%prtvol)
 end if

 ! ==========================
 ! Setup k-points in Sigma_nk
 ! ==========================
 gaps = ebands_get_gaps(ks_ebands, gap_err)
 if (my_rank == master) then
   !call ebands_print(ks_ebands, header="KS band structure", unit=std_out, prtvol=gwr%dtset%prtvol)
   msg = "Kohn-Sham gaps and band edges from IBZ mesh"
   call gaps%print(unit=std_out, header=msg)
   call gaps%print(unit=ab_out, header=msg)
 end if

 ! TODO: nkcalc should be spin dependent.
 ! This piece of code is taken from m_sigmaph.
 ! In principle one should use the same algorithm in setup_sigma (legacy GW code).
 if (dtset%nkptgw /= 0) then

   ! Treat the k-points and bands specified in the input file via kptgw and bdgw.
   call sigtk_kcalc_from_nkptgw(dtset, mband, gwr%nkcalc, gwr%kcalc, gwr%bstart_ks, gwr%nbcalc_ks)

 else

   if (any(abs(dtset%sigma_erange) > zero)) then
     ! Use sigma_erange and (optionally) sigma_ngkpt
     call sigtk_kcalc_from_erange(dtset, cryst, ks_ebands, gaps, gwr%nkcalc, gwr%kcalc, gwr%bstart_ks, gwr%nbcalc_ks, comm)

   else
     ! Use qp_range to select the interesting k-points and the corresponding bands.
     !
     !    0 --> Compute the QP corrections only for the fundamental and the direct gap.
     ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
     !          bands above and below the Fermi level.
     ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
     !          Include all occupied states and `num` empty states.

     qprange_ = dtset%gw_qprange
     if (gap_err /= 0 .and. qprange_ == 0) then
       ABI_WARNING("Cannot compute fundamental and direct gap (likely metal). Will replace qprange 0 with qprange 1")
       qprange_ = 1
     end if

     if (qprange_ /= 0) then
       call sigtk_kcalc_from_qprange(dtset, cryst, ks_ebands, qprange_, gwr%nkcalc, gwr%kcalc, gwr%bstart_ks, gwr%nbcalc_ks)
     else
       ! qprange is not specified in the input.
       ! Include direct and fundamental KS gap or include states depending on the position wrt band edges.
       call sigtk_kcalc_from_gaps(dtset, ks_ebands, gaps, gwr%nkcalc, gwr%kcalc, gwr%bstart_ks, gwr%nbcalc_ks)
     end if
   end if

 end if ! nkptgw /= 0

 ! Include all degenerate states and map kcalc to ibz.
 ! NB: This part is copied from sigmaph.

 ! The k-point and the symmetries connecting the BZ k-point to the IBZ.
 ABI_MALLOC(gwr%kcalc2ibz, (gwr%nkcalc, 6))
 if (abs(gwr%dtset%symsigma) == 1) then
   ABI_MALLOC(gwr%degtab, (gwr%nkcalc, gwr%nsppol))
 end if

 ! Workspace arrays used to compute degeneracy tables.
 ABI_ICALLOC(degblock_all, (2, mband, gwr%nkcalc, gwr%nsppol))
 ABI_ICALLOC(ndeg_all, (gwr%nkcalc, gwr%nsppol))

 ierr = 0

 do ikcalc=1,gwr%nkcalc
   !if (mod(ikcalc, nprocs) /= my_rank) then
   !  gwr%kcalc2ibz(ikcalc, :) = 0
   !  gwr%bstart_ks(ikcalc, :) = 0
   !  gwr%nbcalc_ks(ikcalc, :) = 0
   !  cycle ! MPI parallelism inside comm
   !end if

   ! Note symrel and use_symrel.
   ! These are the conventions for the symmetrization of the wavefunctions used in cgtk_rotate.
   kk = gwr%kcalc(:, ikcalc)

   if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, 1, kk, indkk_k) /= 0) then
      write(msg, '(5a)' ) &
       "The WFK file cannot be used to compute self-energy corrections at k-point: ",trim(ktoa(kk)),ch10,&
       "The k-point cannot be generated from a symmetrical one.", ch10
      ABI_ERROR(msg)
   end if

   ! TODO: Invert dims and update abipy
   gwr%kcalc2ibz(ikcalc, :) = indkk_k(:, 1)

   ik_ibz = indkk_k(1,1); isym_k = indkk_k(2,1)
   trev_k = indkk_k(6,1); g0_k = indkk_k(3:5,1)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !kk_ibz = ks_ebands%kptns(:,ik_ibz)
   if (.not. isirr_k) then
     ABI_WARNING(sjoin("The k-point in Sigma_{nk} must be in the IBZ but got:", ktoa(kk)))
     ierr = ierr + 1
   end if

   ! We will have to average the QP corrections over degenerate states if symsigma=1 is used.
   ! Here we make sure that all the degenerate states are included.
   ! Store also band indices of the degenerate sets, used to average final results.
   if (abs(gwr%dtset%symsigma) == 1) then
     cnt = 0
     do spin=1,gwr%nsppol
       bstop = gwr%bstart_ks(ikcalc, spin) + gwr%nbcalc_ks(ikcalc, spin) - 1
       call ebands_enclose_degbands(ks_ebands, ik_ibz, spin, gwr%bstart_ks(ikcalc, spin), bstop, changed, TOL_EDIFF, &
                                    degblock=degblock)
       if (changed) then
         gwr%nbcalc_ks(ikcalc, spin) = bstop - gwr%bstart_ks(ikcalc, spin) + 1
         cnt = cnt + 1
         if (cnt < 5) then
           write(msg,'(2(a,i0),2a,2(1x,i0))') &
             "Not all the degenerate states for ikcalc: ",ikcalc,", spin: ",spin,ch10, &
             "were included in the bdgw set. bdgw has been automatically changed to: ",gwr%bstart_ks(ikcalc, spin), bstop
           ABI_COMMENT(msg)
         end if
         write(msg,'(2(a,i0),2a)') &
           "The number of included states: ", bstop, &
           " is larger than the number of bands in the input ",dtset%nband(ik_ibz + (spin-1)*ks_ebands%nkpt),ch10,&
           "Action: Increase nband."
         ABI_CHECK(bstop <= dtset%nband(ik_ibz + (spin-1)*ks_ebands%nkpt), msg)
       end if

       ! Store band indices used for averaging (shifted by bstart_ks)
       ndeg = size(degblock, dim=2)
       ndeg_all(ikcalc, spin) = ndeg
       degblock_all(:, 1:ndeg, ikcalc, spin) = degblock(:, 1:ndeg)

       ABI_FREE(degblock)

     end do
   end if ! symsigma
 end do ! ikcalc

 ABI_CHECK(ierr == 0, "kptgw wavevectors must be in the IBZ read from the WFK file.")

 ! Collect data
 !call xmpi_sum(gwr%kcalc2ibz, comm, ierr)
 !call xmpi_sum(gwr%bstart_ks, comm, ierr)
 !call xmpi_sum(gwr%nbcalc_ks, comm, ierr)

 ! Build degtab tables.
 if (abs(gwr%dtset%symsigma) == 1) then
   call xmpi_sum(ndeg_all, comm, ierr)
   call xmpi_sum(degblock_all, comm, ierr)
   do ikcalc=1,gwr%nkcalc
     do spin=1,gwr%nsppol
       ndeg = ndeg_all(ikcalc, spin)
       ABI_MALLOC(gwr%degtab(ikcalc, spin)%bids, (ndeg))
       do ii=1,ndeg
         cnt = degblock_all(2, ii, ikcalc, spin) - degblock_all(1, ii, ikcalc, spin) + 1
         ABI_MALLOC(gwr%degtab(ikcalc, spin)%bids(ii)%vals, (cnt))
         gwr%degtab(ikcalc, spin)%bids(ii)%vals = [(jj, jj= &
           degblock_all(1, ii, ikcalc, spin) - gwr%bstart_ks(ikcalc, spin) + 1, &
           degblock_all(2, ii, ikcalc, spin) - gwr%bstart_ks(ikcalc, spin) + 1)]
       end do
     end do
   end do
 end if

 ABI_FREE(degblock_all)
 ABI_FREE(ndeg_all)

 ! Now we can finally compute max_nbcalc
 gwr%max_nbcalc = maxval(gwr%nbcalc_ks)

 ABI_MALLOC(gwr%bstop_ks, (gwr%nkcalc, gwr%nsppol))
 gwr%bstop_ks = gwr%bstart_ks + gwr%nbcalc_ks - 1

 call krank_ibz%free()
 ABI_FREE(kibz) ! Deallocate kibz here because krank_ibz keeps a reference to this array.

 ! ================================
 ! Setup tau/omega mesh and weights
 ! ================================
 ! Compute min/max transition energy taking into account nsppol if any.
 te_min = minval(gaps%cb_min - gaps%vb_max)
 te_max = maxval(ks_ebands%eig(nbsum,:,:) - ks_ebands%eig(1,:,:))
 if (te_min <= tol6) then
   te_min = tol6
   ABI_WARNING("System is metallic or with a very small fundamental gap!")
 end if

 gwr%ntau = dtset%gwr_ntau

#ifndef __HAVE_GREENX
 ABI_MALLOC(gwr%tau_mesh, (gwr%ntau))
 ABI_MALLOC(gwr%tau_wgs, (gwr%ntau))
 ABI_MALLOC(gwr%iw_mesh, (gwr%ntau))
 ABI_MALLOC(gwr%iw_wgs, (gwr%ntau))
 te_max = one
 gwr%iw_mesh = arth(zero, te_max, gwr%ntau)
 gwr%iw_wgs = te_max / gwr%ntau
 gwr%tau_mesh = arth(zero, te_max, gwr%ntau)
 gwr%tau_wgs = one / gwr%ntau

 ABI_MALLOC(gwr%w2t_cos_wgs, (gwr%ntau, gwr%ntau))
 ABI_MALLOC(gwr%t2w_cos_wgs, (gwr%ntau, gwr%ntau))
 ABI_MALLOC(gwr%t2w_sin_wgs, (gwr%ntau, gwr%ntau))
 gwr%w2t_cos_wgs = one; gwr%t2w_cos_wgs = one; gwr%t2w_sin_wgs = one

#else
 call gx_minimax_grid(gwr%ntau, te_min, te_max,  &  ! in
                      gwr%tau_mesh, gwr%tau_wgs, &  ! all these args are out and allocated by the routine.
                      gwr%iw_mesh, gwr%iw_wgs,   &
                      gwr%t2w_cos_wgs, gwr%w2t_cos_wgs, gwr%t2w_sin_wgs, &
                      gwr%ft_max_error)

 ! Compute the "real" weights used for the inhomogeneous cosine/sine FT and check whether
 ! the two matrices for the forward/backward cosine transforms are the inverse of each other.
 do it=1,gwr%ntau
   do iw=1,gwr%ntau
     gwr%t2w_cos_wgs(iw, it) = gwr%t2w_cos_wgs(iw, it) * cos(gwr%tau_mesh(it) * gwr%iw_mesh(iw))
     gwr%w2t_cos_wgs(it, iw) = gwr%w2t_cos_wgs(it, iw) * cos(gwr%tau_mesh(it) * gwr%iw_mesh(iw))

     !gwr%t2w_cos_wgs(it, iw) = gwr%t2w_cos_wgs(it, iw) * cos(gwr%tau_mesh(it) * gwr%iw_mesh(iw))
     !gwr%w2t_cos_wgs(iw, it) = gwr%w2t_cos_wgs(iw, it) * cos(gwr%tau_mesh(it) * gwr%iw_mesh(iw))
   end do
 end do

 ABI_MALLOC(mat, (gwr%ntau, gwr%ntau))
 mat = matmul(gwr%t2w_cos_wgs, gwr%w2t_cos_wgs)
 !mat = matmul(gwr%w2t_cos_wgs, gwr%t2w_cos_wgs)
 do it=1,gwr%ntau
   mat(it, it) = mat(it, it) - one
 end do
 !print *, "mat", mat
 !call print_arr(mat)
 !mat = matmul(transpose(gwr%t2w_cos_wgs), transpose(gwr%w2t_cos_wgs))
 !mat = matmul(gwr%t2w_cos_wgs, transpose(gwr%w2t_cos_wgs))

 gwr%cosft_duality_err = maxval(abs(mat))

 call wrtout(std_out, sjoin("Max_{ij} |CT CT^{-1} - I|", ftoa(gwr%cosft_duality_err)))
 call wrtout(std_out, sjoin("ft_max_error", ltoa(gwr%ft_max_error)))

 ABI_FREE(mat)
#endif

 call gaps%free()

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================

 if (any(dtset%gwr_np_gtks /= 0)) then
   ! Use MPI parameters from input file.
   gwr%g_comm%nproc    = dtset%gwr_np_gtks(1)
   gwr%tau_comm%nproc  = dtset%gwr_np_gtks(2)
   gwr%kpt_comm%nproc  = dtset%gwr_np_gtks(3)
   gwr%spin_comm%nproc = dtset%gwr_np_gtks(4)
 else
   ! Automatic grid generation.
   !   Priorities        |  MPI Scalability                | Memory
   ! ==================================================================================================
   !   spin (if any)     |  excellent                      | scales
   !   tau               |  excellent                      | scales
   !   kbz               |  good?, requires communication  | scales (depends on BZ -> IBZ mapping)
   !   g/r (PBLAS)       |  network intensive              ! scales
   !
   gwr%spin_comm%nproc = 1
   if (gwr%nsppol == 2 .and. all_nproc > 1) then
     ABI_CHECK(mod(all_nproc, 2) == 0, "when nsppol == 2, nprocs should be even!")
     gwr%spin_comm%nproc = 2
   end if

   nps = all_nproc / gwr%spin_comm%nproc
   do ii=nps,1,-1
     if (mod(gwr%ntau, ii) == 0 .and. mod(nps, ii) == 0) exit
   end do

   if (ii == 1 .and. nps > 1) then
     if (gwr%nkbz > 1) then
       call xmpi_distrib_2d(nps, "12", gwr%ntau, gwr%nkbz, gwr%tau_comm%nproc, gwr%kpt_comm%nproc, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(nps), " with priority: tau/kbz"))
     else
       call xmpi_distrib_2d(nps, "12", gwr%ntau, gwr%green_mpw, gwr%tau_comm%nproc, gwr%g_comm%nproc, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(nps), " with priority: tau/g"))
     end if
   else
     ! ii divides ntau and nps.
     gwr%tau_comm%nproc  = ii
     nps = nps / gwr%tau_comm%nproc

     gwr%kpt_comm%nproc = 1  ! Init values assuming Gamma-only sampling.
     gwr%g_comm%nproc = nps

     if (gwr%nkbz > 1) then
       do ii=nps,1,-1
         if (mod(gwr%nkbz, ii) == 0 .and. mod(nps, ii) == 0) exit
       end do
       if (ii == 1 .and. nps > 1) then
         call xmpi_distrib_2d(nps, "12", gwr%nkbz, gwr%green_mpw, gwr%kpt_comm%nproc, gwr%g_comm%nproc, ierr)
         ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(nps), " with priority: k/g"))
       else
         ! ii divides nkbz and nps.
         gwr%kpt_comm%nproc = ii
         gwr%g_comm%nproc = nps / ii
       end if
     end if

   end if
 end if

 ! Consistency check.
 if (gwr%g_comm%nproc * gwr%tau_comm%nproc * gwr%kpt_comm%nproc * gwr%spin_comm%nproc /= all_nproc) then
   write(msg, "(a,i0,3a, 5(a,1x,i0))") &
     "Cannot create 4d Cartesian grid with total nproc: ", all_nproc, ch10, &
     "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
     "g_nproc (", gwr%g_comm%nproc, ") x tau_nproc (", gwr%tau_comm%nproc, ")  x kpt_nproc (", gwr%kpt_comm%nproc, &
     ")  x spin_nproc (", gwr%spin_comm%nproc, ") != ", all_nproc
   ABI_ERROR(msg)
 end if

 ! ==========================
 ! MPI grid and communicators
 ! ==========================

 dims = [gwr%g_comm%nproc, gwr%tau_comm%nproc, gwr%kpt_comm%nproc, gwr%spin_comm%nproc]
 periods(:) = .False.; reorder = .False.

#ifdef HAVE_MPI
 call MPI_CART_CREATE(gwr%comm, ndims, dims, periods, reorder, comm_cart, ierr)

 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, gwr%coords_gtks, ierr)

 ! Create communicator for g-vectors
 keepdim = .False.; keepdim(1) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%g_comm%value, ierr); gwr%g_comm%me = xmpi_comm_rank(gwr%g_comm%value)

 ! Create communicator for tau
 keepdim = .False.; keepdim(2) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%tau_comm%value, ierr); gwr%tau_comm%me = xmpi_comm_rank(gwr%tau_comm%value)

 ! Create communicator for k-points
 keepdim = .False.; keepdim(3) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%kpt_comm%value, ierr); gwr%kpt_comm%me = xmpi_comm_rank(gwr%kpt_comm%value)

 ! Create communicator for spin
 keepdim = .False.; keepdim(4) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%spin_comm%value, ierr); gwr%spin_comm%me = xmpi_comm_rank(gwr%spin_comm%value)

 ! Create communicator for the (g, tau) 2D grid.
 keepdim = .False.; keepdim(1) = .True.; keepdim(2) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%gtau_comm%value, ierr); gwr%gtau_comm%me = xmpi_comm_rank(gwr%gtau_comm%value)

 ! Create communicator for the (g, tau, k) 3D subgrid.
 keepdim = .True.; keepdim(4) = .False.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%gtk_comm%value, ierr); gwr%gtk_comm%me = xmpi_comm_rank(gwr%gtk_comm%value)

 ! Create communicator for the (tau, k, spin) 3D subgrid.
 keepdim = .True.; keepdim(1) = .False.
 call MPI_CART_SUB(comm_cart, keepdim, gwr%tks_comm%value, ierr); gwr%tks_comm%me = xmpi_comm_rank(gwr%tks_comm%value)

 call xmpi_comm_free(comm_cart)
#endif

 ! Block-distribute dimensions and allocate redirection table: local index --> global index.
 call xmpi_split_block(gwr%ntau, gwr%tau_comm%value, gwr%my_ntau, gwr%my_itaus)
 ABI_CHECK(gwr%my_ntau > 0, "my_ntau == 0, decrease number of procs for tau level")

 ! Store the rank of the MPI proc in tau_comm treating itau.
 ABI_MALLOC(gwr%tau_master, (gwr%ntau))
 gwr%tau_master = -1
 do my_it=1,gwr%my_ntau
   itau = gwr%my_itaus(my_it)
   gwr%tau_master(itau) = gwr%tau_comm%me
 end do
 call xmpi_max_ip(gwr%tau_master, gwr%tau_comm%value, ierr)

 call xmpi_split_block(gwr%nsppol, gwr%spin_comm%value, gwr%my_nspins, gwr%my_spins)
 ABI_CHECK(gwr%my_nspins > 0, "my_nspins == 0, decrease number of procs for spin level")

 ! Distribute k-points in the full BZ, transfer symmetry tables.
 ! Finally, find the number of IBZ points to be stored in memory by this MPI rank.

 call xmpi_split_block(gwr%nkbz, gwr%kpt_comm%value, gwr%my_nkbz, gwr%my_kbz_inds)
 ABI_CHECK(gwr%my_nkbz > 0, "my_nkbz == 0, decrease number of procs for k-point level")

 ABI_MALLOC(gwr%my_kbz2ibz, (6, gwr%my_nkbz))
 gwr%my_kbz2ibz = kbz2ibz(:, gwr%my_kbz_inds(:))

 ABI_ICALLOC(gwr%np_kibz, (gwr%nkibz))
 do my_ikf=1,gwr%my_nkbz
   ik_ibz = gwr%my_kbz2ibz(1, my_ikf)
   gwr%np_kibz(ik_ibz) = 1
 end do

 gwr%my_nkibz = count(gwr%np_kibz > 0)
 ABI_MALLOC(gwr%my_kibz_inds, (gwr%my_nkibz))
 ii = 0
 do ik_ibz=1,gwr%nkibz
   if (gwr%np_kibz(ik_ibz) > 0) then
     ii = ii + 1; gwr%my_kibz_inds(ii) = ik_ibz
   end if
 end do

 call xmpi_sum(gwr%np_kibz, gwr%kpt_comm%value, ierr)

 ABI_FREE(kbz2ibz)

 ! Distribute q-points in the full BZ, transfer symmetry tables.
 ! Finally find the number of IBZ points that should be stored in memory.
 call xmpi_split_block(gwr%nqbz, gwr%kpt_comm%value, gwr%my_nqbz, gwr%my_qbz_inds)
 ABI_MALLOC(gwr%my_qbz2ibz, (6, gwr%my_nqbz))
 gwr%my_qbz2ibz = qbz2ibz(:, gwr%my_qbz_inds(:))

 ABI_ICALLOC(gwr%np_qibz, (gwr%nqibz))

 do my_iqf=1,gwr%my_nqbz
   iq_ibz = gwr%my_qbz2ibz(1, my_iqf)
   gwr%np_qibz(iq_ibz) = 1
 end do

 gwr%my_nqibz = count(gwr%np_qibz > 0)
 ABI_MALLOC(gwr%my_qibz_inds, (gwr%my_nqibz))
 ii = 0
 do iq_ibz=1,gwr%nqibz
   if (gwr%np_qibz(iq_ibz) > 0) then
     ii = ii + 1; gwr%my_qibz_inds(ii) = iq_ibz
   end if
 end do

 call xmpi_sum(gwr%np_qibz, gwr%kpt_comm%value, ierr)

 ABI_FREE(qbz2ibz)

 ! =========================================
 ! Find FFT mesh and max number of g-vectors
 ! =========================================
 ! Note the usage of gwr_boxcutmin and loops over the full BZ.

 ecut_eff = dtset%ecut * dtset%dilatmx ** 2
 call ngfft_seq(gwr%g_ngfft, [0, 0, 0])

 do ik_bz=1,gwr%nkbz
   kk_bz = gwr%kbz(:, ik_bz)
   call get_kg(kk_bz, istwfk1, ecut_eff, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%gwr_boxcutmin, dtset%chksymtnons, ecut_eff, cryst%gmet, &
              kk_bz, me_fft0, gwr%g_mgfft, gwr%g_nfft, gwr%g_ngfft, nproc_fft1, cryst%nsym, paral_fft0, &
              cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%green_mpw = max(gwr%green_mpw, npw_)
 end do

 do iq_bz=1,gwr%nqbz
   qq_bz = gwr%qbz(:, iq_bz)
   call get_kg(qq_bz, istwfk1, dtset%ecuteps, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%gwr_boxcutmin, dtset%chksymtnons, dtset%ecuteps, cryst%gmet, &
              qq_bz, me_fft0, gwr%g_mgfft, gwr%g_nfft, gwr%g_ngfft, nproc_fft1, cryst%nsym, &
              paral_fft0, cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%tchi_mpw = max(gwr%tchi_mpw, npw_)
 end do

 ! TODO: For the time being no augmentation
 gwr%g_ngfft(4:6) = gwr%g_ngfft(1:3)

 if (my_rank == master) then
   call print_ngfft(gwr%g_ngfft, header="FFT mesh for GWR", unit=std_out)
   call print_ngfft(gwr%g_ngfft, header="FFT mesh for GWR", unit=ab_out)
 end if

 ! Now we know the value of g_ngfft. Setup tables for zero-padded FFTs.
 ! Build descriptors for Green's functions and tchi and setup tables for zero-padded FFTs.
 ABI_MALLOC(gwr%green_desc_kibz, (gwr%nkibz))

 do my_iki=1,gwr%my_nkibz
   ik_ibz = gwr%my_kibz_inds(my_iki)
   kk_ibz = gwr%kibz(:, ik_ibz)
   associate (desc_k => gwr%green_desc_kibz(ik_ibz))
   desc_k%istwfk = istwfk1
   call get_kg(kk_ibz, desc_k%istwfk, ecut_eff, gwr%cryst%gmet, desc_k%npw, desc_k%gvec)
   ABI_MALLOC(desc_k%gbound, (2 * gwr%g_mgfft + 8, 2))
   call sphereboundary(desc_k%gbound, desc_k%istwfk, desc_k%gvec, gwr%g_mgfft, desc_k%npw)
   end associate
 end do

 ABI_MALLOC(gwr%tchi_desc_qibz, (gwr%nqibz))

 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   q_is_gamma = (normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0)
   associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))
   desc_q%istwfk = istwfk1
   call get_kg(qq_ibz, desc_q%istwfk, dtset%ecuteps, gwr%cryst%gmet, desc_q%npw, desc_q%gvec)
   ABI_MALLOC(desc_q%gbound, (2 * gwr%g_mgfft + 8, 2))
   call sphereboundary(desc_q%gbound, desc_q%istwfk, desc_q%gvec, gwr%g_mgfft, desc_q%npw)
   ! Compute v(q,G). FIXME: Implement cutoff in vc
   ABI_MALLOC(desc_q%vc_sqrt, (desc_q%npw))
   ig_start = 1
   if (q_is_gamma) then
     ig_start = 2
     desc_q%vc_sqrt(1) = sqrt(four_pi) / normv(GW_Q0_DEFAULT, gwr%cryst%gmet, "G")
   end if
   do ig=ig_start,desc_q%npw
     desc_q%vc_sqrt(ig) = sqrt(four_pi) / normv(qq_ibz + desc_q%gvec(:,ig), gwr%cryst%gmet, "G")
   end do
   end associate
 end do

 ! ==========================================================
 ! Allocate PBLAS arrays for G_kibz(g,g') and tchi_qibz(g,g')
 ! ==========================================================
 ABI_MALLOC(gwr%gt_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))
 ABI_MALLOC(gwr%tchi_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))

 ! 1D grid block-distributed along columns.
 call gwr%g_slkproc%init(gwr%g_comm%value, grid_dims=[1, gwr%g_comm%nproc])

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_it=1,gwr%my_ntau
     itau = gwr%my_itaus(my_it)

     ! Allocate G_k(g,g') for k in my IBZ, MPI distributed over g' in blocks
     do my_iki=1,gwr%my_nkibz
       ik_ibz = gwr%my_kibz_inds(my_iki)
       npwsp = gwr%green_desc_kibz(ik_ibz)%npw * gwr%nspinor
       col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
       associate (gt => gwr%gt_kibz(:, ik_ibz, itau, spin))
       call gt(1)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
       call gt(2)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
       !call gt(1)%print(header=sjoin("gt_kibz(1) for ik_ibz:", itoa(ik_ibz)))
       end associate
     end do

    ! Allocate tchi_q(g,g') for q in my IBZ, MPI distributed over g' in blocks
    do my_iqi=1,gwr%my_nqibz
      iq_ibz = gwr%my_qibz_inds(my_iqi)
      npwsp = gwr%tchi_desc_qibz(iq_ibz)%npw * gwr%nspinor
      col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
      associate (tchi => gwr%tchi_qibz(iq_ibz, itau, spin))
      call tchi%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
      !call tchi%print(header=sjoin("tchi_qibz for iq_ibz:", itoa(iq_ibz)))
      end associate
    end do

   end do ! my_it
 end do ! my_is

 mem_mb = slk_array_locmem_mb(gwr%gt_kibz) + slk_array_locmem_mb(gwr%tchi_qibz)
 write(msg,'(a,f8.1,a)')' Local memory needed for G + tchi PBLAS matrices: ', mem_mb, ' [Mb] <<< MEM'
 call wrtout(std_out, msg)

 ! Create netcdf file to store results.
 gwr%gwrnc_path = strcat(dtfil%filnam_ds(4), "_GWR.nc")

 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, gwr%gwrnc_path, xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ks_ebands, ncid))

   ! Add dimensions.
   !ncerr = nctk_def_dims(ncid, [ &
   !  nctkdim_t("nsppol", gwr%nsppol), &
   !  !nctkdim_t("nkcalc", gwr%nkcalc), nctkdim_t("max_nbcalc", gwr%max_nbcalc), &
   !  nctkdim_t("nqibz", gwr%nqibz), nctkdim_t("nqbz", gwr%nqbz)], &
   !  defmode=.True.)
   !NCF_CHECK(ncerr)

   !ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", "symdynmat", &
   !  "ph_intmeth", "eph_intmeth", "qint_method", "eph_transport", &
   !  "imag_only", "symv1scf", "dvdb_add_lr", "mrta", "ibte_prep", "eph_prtscratew"])
   !NCF_CHECK(ncerr)
   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", &
   !  "ph_wstep", "ph_smear", "eph_phwinfact"])
   !NCF_CHECK(ncerr)

   ! Define arrays with results.
   !ncerr = nctk_def_arrays(ncid, [ &
   !  nctkarr_t("ngqpt", "int", "three"), &
   !  nctkarr_t("ddb_ngqpt", "int", "three"), &
   !  nctkarr_t("ph_ngqpt", "int", "three"), &
   !  nctkarr_t("sigma_ngkpt", "int", "three"), &
   !  nctkarr_t("sigma_erange", "dp", "two"), &
   !  nctkarr_t("bstart_ks", "int", "nkcalc, nsppol"), &
   !  !nctkarr_t("bstop_ks", "int", "nkcalc, nsppol"), &
   !  nctkarr_t("nbcalc_ks", "int", "nkcalc, nsppol"), &
   !  nctkarr_t("kcalc", "dp", "three, nkcalc"), &
   !  nctkarr_t("kcalc2ibz", "int", "nkcalc, six"), &
   !  nctkarr_t("qp_done", "int", "nkcalc, nsppol"), &
   !  nctkarr_t("vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
   !  nctkarr_t("dvals_de0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
   !  nctkarr_t("qpoms_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
   !  nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
   !  nctkarr_t("ze0_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"), &
   !  nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol"), &
   !  nctkarr_t("ks_gaps", "dp", "nkcalc, nsppol"), &
   !  nctkarr_t("qpoms_gaps", "dp", "ntemp, nkcalc, nsppol"), &
   !  nctkarr_t("qp_gaps", "dp", "ntemp, nkcalc, nsppol"), &
   !  nctkarr_t("vcar_calc", "dp", "three, max_nbcalc, nkcalc, nsppol") &
   !])
   !NCF_CHECK(ncerr)

   ! ======================================================
   ! Write data that do not depend on the (kpt, spin) loop.
   ! ======================================================
   !NCF_CHECK(nctk_set_datamode(ncid))
   !ii = 0; if (gwr%imag_only) ii = 1
   !ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", &
   !  "symdynmat", "ph_intmeth", "eph_intmeth", "qint_method", &
   !  "eph_transport", "imag_only", "symv1scf", "dvdb_add_lr", "mrta", "ibte_prep", "eph_prtscratew"], &
   !  [dtset%eph_task, gwr%symsigma, gwr%nbsum, gwr%bsum_start, gwr%bsum_stop, &
   !  dtset%symdynmat, dtset%ph_intmeth, dtset%eph_intmeth, gwr%qint_method, &
   !  dtset%eph_transport, ii, dtset%symv1scf, dtset%dvdb_add_lr, gwr%mrta, dtset%ibte_prep, dtset%eph_prtscratew])
   !NCF_CHECK(ncerr)
   !ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear", "eph_phwinfact"], &
   !  [aimag(gwr%ieta), gwr%wr_step, dtset%eph_fsewin, dtset%eph_fsmear, dtset%eph_extrael, dtset%eph_fermie, &
   !  dtset%ph_wstep, dtset%ph_smear, dtset%eph_phwinfact])
   !NCF_CHECK(ncerr)

   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt"), gwr%ngqpt))
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ph_ngqpt"), dtset%ph_ngqpt))

   NCF_CHECK(nf90_close(ncid))
 end if

 ! Now reopen the file inside ncwrite_comm to perform parallel-IO
 !call xmpi_barrier(comm)
 !if (gwr%ncwrite_comm%value /= xmpi_comm_null) then
 !  NCF_CHECK(nctk_open_modify(gwr%ncid, path, gwr%ncwrite_comm%value))
 !  NCF_CHECK(nctk_set_datamode(gwr%ncid))
 !end if

 if (my_rank == master) then
   call gwr%print(unit=ab_out)
   call gwr%print(unit=std_out)
 end if

 call cwtime_report(" gwr_init:", cpu, wall, gflops)

end subroutine gwr_init
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_free
!! NAME
!! gwr_free
!!
!! FUNCTION
!!  Free memory
!!
!! SOURCE

subroutine gwr_free(gwr)

!Arguments ------------------------------------
 class(gwr_t), intent(inout) :: gwr

!Local variables-------------------------------
 !integer :: ii, jj

! *************************************************************************

 ABI_SFREE(gwr%kbz)
 ABI_SFREE(gwr%my_kbz2ibz)
 ABI_SFREE(gwr%my_qbz2ibz)
 ABI_SFREE(gwr%my_kbz_inds)
 ABI_SFREE(gwr%my_kibz_inds)
 ABI_SFREE(gwr%my_qbz_inds)
 ABI_SFREE(gwr%my_qibz_inds)
 ABI_SFREE(gwr%qbz)
 ABI_SFREE(gwr%qibz)
 ABI_SFREE(gwr%wtq)
 ABI_SFREE(gwr%qbz2ibz)
 ABI_SFREE(gwr%my_spins)
 ABI_SFREE(gwr%my_itaus)
 ABI_SFREE(gwr%tau_master)
 ABI_SFREE(gwr%np_kibz)
 ABI_SFREE(gwr%np_qibz)
 ABI_SFREE(gwr%tau_mesh)
 ABI_SFREE(gwr%tau_wgs)
 ABI_SFREE(gwr%iw_mesh)
 ABI_SFREE(gwr%iw_wgs)
 ABI_SFREE(gwr%w2t_cos_wgs)
 ABI_SFREE(gwr%t2w_cos_wgs)
 ABI_SFREE(gwr%t2w_sin_wgs)
 ABI_SFREE(gwr%kcalc)
 ABI_SFREE(gwr%bstart_ks)
 ABI_SFREE(gwr%bstop_ks)
 ABI_SFREE(gwr%nbcalc_ks)
 ABI_SFREE(gwr%kcalc2ibz)

 call ebands_free(gwr%qp_ebands)
 call gwr%kcalc_wfd%free()

 ! Free descriptors
 if (allocated(gwr%green_desc_kibz)) then
   call desc_array_free(gwr%green_desc_kibz)
   ABI_FREE(gwr%green_desc_kibz)
 end if
 if (allocated(gwr%tchi_desc_qibz)) then
   call desc_array_free(gwr%tchi_desc_qibz)
   ABI_FREE(gwr%tchi_desc_qibz)
 end if

 ! Free PBLAS matrices
 if (allocated(gwr%gt_kibz)) then
   call slk_array_free(gwr%gt_kibz)
   ABI_FREE(gwr%gt_kibz)
 end if
 if (allocated(gwr%tchi_qibz)) then
   call slk_array_free(gwr%tchi_qibz)
   ABI_FREE(gwr%tchi_qibz)
 end if
 if (allocated(gwr%wc_qibz)) then
   call slk_array_free(gwr%wc_qibz)
   ABI_FREE(gwr%wc_qibz)
 end if
 if (allocated(gwr%sigc_kibz)) then
   call slk_array_free(gwr%sigc_kibz)
   ABI_FREE(gwr%sigc_kibz)
 end if
 call gwr%g_slkproc%free()

 call gwr%ks_me%free()

 ! datatypes.
 if (allocated(gwr%degtab)) then
   call degtab_array_free(gwr%degtab)
   ABI_FREE(gwr%degtab)
 end if

 ! Free MPI communicators
 call gwr%spin_comm%free()
 call gwr%g_comm%free()
 call gwr%tau_comm%free()
 call gwr%kpt_comm%free()
 call gwr%gtau_comm%free()
 call gwr%gtk_comm%free()
 call gwr%tks_comm%free()

end subroutine gwr_free
!!***

subroutine desc_array1_free(desc_array1)
  type(desc_t),intent(inout) :: desc_array1(:)
  integer :: ii
  do ii=1,size(desc_array1, dim=1)
    call desc_array1(ii)%free()
  end do
end subroutine desc_array1_free

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_load_kcalc_wfd
!! NAME
!! gwr_load_kcalc_wfd
!!
!! FUNCTION
!!  Load the KS states for Sigma_nk from the WFK file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_load_kcalc_wfd(gwr, wfk_path, tmp_kstab)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 character(len=*),intent(in) :: wfk_path
 integer,allocatable,intent(out) :: tmp_kstab(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: mband, nkibz, nsppol, spin, ik_ibz, ikcalc
 real(dp) :: cpu, wall, gflops
 !character(len=5000) :: msg
 type(ebands_t) :: ks_ebands
 type(hdr_type) :: wfk_hdr
!arrays
 integer,allocatable :: nband(:,:), wfd_istwfk(:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 associate (wfd => gwr%kcalc_wfd, dtset => gwr%dtset)

 ks_ebands = wfk_read_ebands(wfk_path, gwr%comm, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ! TODO: More consistency checks e.g. nkibz,...

 !cryst = wfk_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 nkibz = ks_ebands%nkpt; nsppol = ks_ebands%nsppol !; mband = ks_ebands%mband

 ! Don't take mband from ks_ebands but compute it from gwr%bstop_ks
 mband = maxval(gwr%bstop_ks)

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k wavevectors
 ! treated by this MPI rank are stored.
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 ABI_ICALLOC(tmp_kstab, (2, nkibz, nsppol))

 do spin=1,gwr%nsppol
   do ikcalc=1,gwr%nkcalc ! TODO: Should be spin dependent!
     ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
     associate (b1 => gwr%bstart_ks(ikcalc, spin), b2 => gwr%bstop_ks(ikcalc, spin))
     tmp_kstab(:, ik_ibz, spin) = [b1, b2]
     bks_mask(b1:b2, ik_ibz, spin) = .True.
     end associate
   end do
 end do

 ! Impose istwfk = 1 for all k-points.
 ! wfd_read_wfk will handle a possible conversion if the WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 call wfd_init(wfd, gwr%cryst, gwr%pawtab, gwr%psps, keep_ur, mband, nband, nkibz, dtset%nsppol, bks_mask, &
   dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ks_ebands%kptns, gwr%g_ngfft, &
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, gwr%comm)

 call wfd%print(header="Wavefunctions for GWR calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 call ebands_free(ks_ebands)
 call wfk_hdr%free()

 ! Read KS wavefunctions.
 call wfd%read_wfk(wfk_path, iomode_from_fname(wfk_path))
 end associate

 call cwtime_report(" gwr_load_kcalc_from_wfk:", cpu, wall, gflops)

end subroutine gwr_load_kcalc_wfd
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_gtau_from_wfk
!! NAME
!!  gwr_build_gtau_from_wfk
!!
!! FUNCTION
!!  Build the Green's function in imaginary time from the WFK file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_gtau_from_wfk(gwr, wfk_path)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 character(len=*),intent(in) :: wfk_path

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0 = 0, istwfk1 = 1
 integer :: mband, nkibz, nsppol, my_is, my_iki, spin, ik_ibz ! my_it,
 integer :: my_ib, my_nb, band, npw_k, mpw, istwf_k, itau, ipm, ig1, ig2, il_g1, il_g2
 integer :: nbsum, npwsp, col_bsize, nband_k, my_bstart, my_bstop, nb, ierr
 real(dp) :: f_nk, eig_nk, ef, spin_fact ! tau,
 real(dp) :: cpu, wall, gflops
 complex(dp) :: gt_cfact
 character(len=5000) :: msg
 type(wfd_t),target :: wfd
 type(wave_t),pointer :: wave
 type(ebands_t) :: wfk_ebands
 type(hdr_type) :: wfk_hdr
 type(wfk_t) :: wfk
 type(processor_scalapack) :: gtau_slkproc
 type(matrix_scalapack), target :: cg_mat
!arrays
 integer :: lsize(2)
 integer,allocatable :: nband(:,:), wfd_istwfk(:), my_bands(:), kg_k(:,:)
 real(dp) :: kk_ibz(3)
 !real(dp),allocatable :: cgwork(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: cg_k(:,:)
 complex(dp),allocatable :: loc_cwork(:,:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 associate (dtset => gwr%dtset)

 wfk_ebands = wfk_read_ebands(wfk_path, gwr%comm, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ! TODO: More consistency checks e.g. nkibz,...

 !cryst = wfk_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k-points
 ! treated by this MPI rank are stored in memory.

 nkibz = wfk_ebands%nkpt; nsppol = wfk_ebands%nsppol; mband = wfk_ebands%mband

 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 nbsum = dtset%nband(1)
 if (nbsum > mband) then
   ABI_WARNING(sjoin("WFK file contains", itoa(mband), "states while you want to use:", itoa(nbsum)))
   nbsum = mband
 end if
 call wrtout(std_out, sjoin(" Computing Green's function with nbsum:", itoa(nbsum)))

 ! Init bks_mask:
 !
 !  - each MPI proc reads only the (k-points, spin) it needs.
 !  - For given (k, s), bands are distributed inside tau_comm
 !    This means that bands are replicated across g_comm.
 !
 ! TODO: Use (g_comm x tau_comm) 2D grid to distribute bands
 !       Requires PBLAS redistribuution though.

 call xmpi_split_block(nbsum, gwr%tau_comm%value, my_nb, my_bands)
 ABI_CHECK(my_nb > 0, "my_nb == 0")

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     bks_mask(my_bands(:), ik_ibz, spin) = .True.
   end do
 end do

 ! Impose istwfk = 1 for all k-points.
 ! wfd_read_wfk will handle a possible conversion if the WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 call wfd_init(wfd, gwr%cryst, gwr%pawtab, gwr%psps, keep_ur, mband, nband, nkibz, dtset%nsppol, bks_mask, &
   dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, wfk_ebands%kptns, gwr%g_ngfft, &
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, gwr%comm)

 call wfd%print(header="Wavefunctions for GWR calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 call ebands_free(wfk_ebands)
 call wfk_hdr%free()

 ! Read KS wavefunctions.
 call wfd%read_wfk(wfk_path, iomode_from_fname(wfk_path))

 ! ==============================================
 ! Build Green's functions in g-space for given k
 ! ==============================================

 ! for tau > 0:
 !
 !      G_k(r,r',itau) = i \sum_b^{occ} psi_b(r) \psi_b^*(r') exp(e_b tau)
 !
 ! for tau < 0:
 !
 !      G_k(r,r',itau) = -i \sum_b^{empty} psi_b(r) \psi_b^*(r') exp(e_b tau)


 ! NB: G_k is constructed for k in the IBZ, then we rotate the k-point to obtain G_k in the BZ.
 !
 ! TODO:
 !     1) Make sure that gvec in gwr and wfd agree with each other.
 !     2) May implement the following trick to accelerate convergence wrt nband:
 !         - Init nb states with random numbers and orthogonalize wrt the nband states found in the WFK file
 !         - Compute <i|H|j> in the subspace spanned by the perp states
 !         - Diagonalize H in the subspace to get variational approximation to the KS states
 !         - Add these approximated eigenstats to G.
 !
 ! TODO: Mv m_ddk in 72_response below 70_gw or move gwr to higher level.
 !ddkop = ddkop_new(dtset, gwr%cryst, gwr%pawtab, gwr%psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 mpw = maxval(wfd%npwarr)
 ef = gwr%ks_ebands%fermie
 spin_fact = two / (gwr%nsppol * gwr%nspinor)

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     kk_ibz = gwr%kibz(:, ik_ibz)
     npw_k = wfd%npwarr(ik_ibz)
     istwf_k = wfd%istwfk(ik_ibz)
     associate (desc_k => gwr%green_desc_kibz(ik_ibz))

     ABI_CHECK(npw_k == desc_k%npw, "npw_k != desc_k%npw")
     ABI_CHECK(all(wfd%kdata(ik_ibz)%kg_k == desc_k%gvec), "kg_k != desc_k%gvec")
     !ABI_MALLOC(cgwork, (2, npw_k * wfd%nspinor))
     !call ddkop%setup_spin_kpoint(dtset, gwr%cryst, gwr%psps, spin, kk_ibz, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)
     end associate

     ! Loop over global set of taus.
     ! For each tau, sum local set of bands distributed inside tau_comm.
     ! Each proc store partial results in loc_cwork dimensioned as PBLAS matrix.
     ! Finally, reduce PBLAS data on the processor treating tau inside tau_comm.

     !call wfd%copy_cg(band_ks, ik_ibz, spin, cgwork)
     !eig0nk = ebands%eig(band_ks, ik_ibz, spin)
     !sigma%vcar_calc(:, ib_k, ikcalc, spin) = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cgwork, cwaveprj0)
     ! TODO: spinor and use BLAS2 but take into account PBLAS distribution
     ! or rewrite everything with PBLAS routines.

     ! This just to get the dimension of the local buffer from the first PBLAS matrix I treat.
     itau = gwr%my_itaus(1)
     associate (gt0 => gwr%gt_kibz(1, ik_ibz, itau, spin))
     lsize = gt0%sizeb_local
     ABI_MALLOC(loc_cwork, (lsize(1), lsize(2), 2))

     do itau=1,gwr%ntau

       ! Sum over bands treated by me.
       loc_cwork = zero
       do my_ib=1,my_nb
         band = my_bands(my_ib)
         f_nk = gwr%ks_ebands%occ(band, ik_ibz, spin)
         eig_nk = gwr%ks_ebands%eig(band, ik_ibz, spin) - ef

         ! Select occupied or empty G.
         if (eig_nk < tol6) then
           ipm = 1
           gt_cfact = j_dpc * exp(gwr%tau_mesh(itau) * eig_nk)
           ! Vasp convention
           !ipm = 2
           !gt_cfact = exp(gwr%tau_mesh(itau) * eig_nk)
         else if (eig_nk > tol6) then
           ipm = 2
           gt_cfact = -j_dpc * exp(-gwr%tau_mesh(itau) * eig_nk)
           ! Vasp convention
           !ipm = 1
           !gt_cfact = -exp(-gwr%tau_mesh(itau) * eig_nk)
         else
           ABI_WARNING("Metallic system of semiconductor with Fermi level inside bands!!!!")
         end if

         ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

         if (.not. wave%has_ug == WFD_STORED) then
           write(msg,'(a,3(i0,1x),a)')" ug for (band, ik_ibz, spin): ",band, ik_ibz, spin," is not stored in memory!"
           ABI_ERROR(msg)
         end if

         do il_g2=1, gt0%sizeb_local(2)
           ig2 = gt0%loc2gcol(il_g2)
           do il_g1=1, gt0%sizeb_local(1)
             ig1 = gt0%loc2grow(il_g1)
             loc_cwork(il_g1, il_g2, ipm) = loc_cwork(il_g1, il_g2, ipm) &
               + gt_cfact * wave%ug(ig1) * GWPC_CONJG(wave%ug(ig2))
           end do
         end do
       end do  ! my_ib

       call xmpi_sum_master(loc_cwork, gwr%tau_master(itau), gwr%tau_comm%value, ierr)
       if (gwr%tau_comm%me == gwr%tau_master(itau)) then
         do ipm=1,2
           gwr%gt_kibz(ipm, ik_ibz, itau, spin)%buffer_cplx = loc_cwork(:,:,ipm)
         end do
       end if

     end do ! itau

     ABI_FREE(loc_cwork)
     end associate

   end do ! my_ikf
 end do ! my_is

 !call ddkop%free()
 call wfd%free()
 ABI_FREE(my_bands)

 if (gwr%dtset%prtvol > 0) call gwr_print_trace(gwr, "gt_kibz")

 call cwtime_report(" gwr_build_gtau_from_wfk:", cpu, wall, gflops)
 end associate

 RETURN

 ! PBLAS-based version.
 call wfk_open_read(wfk, wfk_path, formeig0, iomode_from_fname(wfk_path), get_unit(), gwr%comm, hdr_out=wfk_hdr)

 mpw = maxval(wfk_hdr%npwarr)
 ABI_MALLOC(kg_k, (3, mpw))

 call gtau_slkproc%init(gwr%gtau_comm%value, grid_dims=[1, gwr%gtau_comm%nproc])

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     npw_k = wfk_hdr%npwarr(ik_ibz)
     istwf_k = wfk_hdr%istwfk(ik_ibz)
     nband_k = wfk_hdr%nband(ik_ibz)
     npwsp = npw_k * gwr%nspinor

     ! Init CG(npw_k * nspinor, nband) PBLAS matrix within gtau communicator.
     ! and distribute it over bands so that each proc reads a subset of bands in read_band_block
     col_bsize = nband_k / gwr%gtau_comm%nproc; if (mod(nband_k, gwr%gtau_comm%nproc) /= 0) col_bsize = col_bsize + 1
     call cg_mat%init(npwsp, nband_k, gtau_slkproc, istwfk1, size_blocs=[npwsp, col_bsize])

     my_bstart = cg_mat%loc2gcol(1)
     my_bstop = cg_mat%loc2gcol(cg_mat%sizeb_local(2))
     nb = my_bstop - my_bstart + 1
     ABI_CHECK(nb > 0, "nb == 0, decrease number of procs for G and tau parallelism.")

     ! Read my bands
     call c_f_pointer(c_loc(cg_mat%buffer_cplx), cg_k, shape=[2, npwsp * nb])
     call wfk%read_band_block([my_bstart, my_bstop], ik_ibz, spin, xmpio_single, kg_k=kg_k, cg_k=cg_k)

     !call cg_mat%copy(cg_work)
     !call cg_work%free()

     do itau=1,gwr%ntau
       !gt_cfact = exp(-gwr%tau_mesh(itau)  * eig_nk)
       !call slk_pzgemm(transa, transb, matrix1, cone, matrix2, czero, results)
       !my_it = gwr%myit_from_itau(itau)
       !if (my_it == -1) cycle
       !associate (gt => gwr%gt_kibz(ipm, ik_ibz, itau, spin))
     end do

     call cg_mat%free()
   end do ! my_iki
 end do ! my_is

 call gtau_slkproc%free()
 ABI_FREE(kg_k)

 call wfk%close()
 call wfk_hdr%free()

end subroutine gwr_build_gtau_from_wfk
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_rotate_gt
!! NAME
!!  gwr_rotate_gt
!!
!! FUNCTION
!!  Reconstruct the Green's functions in the kBZ from the IBZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!
!!    If q_bz = S q_ibz + G0:
!!
!!      $\epsilon^{-1}_{SG1-G0, SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau} \epsilon^{-1}_{G1, G2)}(q)
!!
!!    If time-reversal symmetry can be used then:
!!
!!      $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau} \epsilon^{-1}_{-S^{-1}(G1+Go), -S^{-1}(G2+G0)}^*(q)
!!
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G - G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!! SOURCE

subroutine gwr_rotate_gt(gwr, my_ikf, my_it, my_is, desc_kbz, gt_kbz)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 integer,intent(in) :: my_ikf, my_it, my_is
 type(desc_t),intent(out) :: desc_kbz
 type(matrix_scalapack),intent(inout) :: gt_kbz(2)

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, ipm, spin, itau, ik_ibz, isym_k, trev_k, g0_k(3), ik_bz
 logical :: isirr_k
 real(dp) :: tsign
!arrays
 integer :: g1(3), g2(3)
 real(dp) :: tnon(3)
 complex(dp) :: ph2, ph1
 !real(dp) :: cpu, wall, gflops

! *************************************************************************

 !call cwtime(cpu, wall, gflops, "start")

 spin = gwr%my_spins(my_is)
 itau = gwr%my_itaus(my_it)
 ik_bz = gwr%my_kbz_inds(my_ikf)

 ik_ibz = gwr%my_kbz2ibz(1, my_ikf); isym_k = gwr%my_kbz2ibz(2, my_ikf)
 trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 tsign = merge(1, -1, trev_k == 0)
 !ABI_CHECK(all(g0_k == 0), sjoin("For kbz:", ktoa(gwr%kbz(:, ik_bz)), "g0_k:", ltoa(g0_k), " != 0"))

 ! Copy descriptor from IBZ
 associate (desc_kibz => gwr%green_desc_kibz(ik_ibz))
 call desc_kibz%copy(desc_kbz)

 if (isirr_k) then
   ! Copy the Green functions and we are done
   do ipm=1,2
     call gwr%gt_kibz(ipm, ik_ibz, itau, spin)%copy(gt_kbz(ipm))
   end do
   goto 10
 end if

 !ABI_ERROR("green: k should be irred for the time being")

 ! Rotate gvec, recompute gbound and rotate vc_sqrt
 ! TODO: 1) Handle TR and routine to rotate tchi/W including vc_sqrt
 !       2) Make sure that FFT box is large enough to accomodate umklapps
 do ig1=1,desc_kbz%npw
   desc_kbz%gvec(:,ig1) = matmul(gwr%cryst%symrec(:,:,isym_k), desc_kibz%gvec(:,ig1)) - g0_k
 end do

 call sphereboundary(desc_kbz%gbound, desc_kbz%istwfk, desc_kbz%gvec, gwr%g_mgfft, desc_kbz%npw)

 ! Get G_k in the BZ.
 tnon = gwr%cryst%tnons(:, isym_k)
 do ipm=1,2
   associate (gk_i => gwr%gt_kibz(ipm, ik_ibz, itau, spin), gk_f => gt_kbz(ipm))
   call gk_i%copy(gk_f)
   do il_g2=1, gk_f%sizeb_local(2)
     ig2 = mod(gk_f%loc2gcol(il_g2) -1, desc_kbz%npw) + 1
     !print *, "ig2, npw", ig2, desc_kbz%npw
     g2 = desc_kbz%gvec(:,ig2)
     ph2 = exp(-j_dpc * two_pi * dot_product(g2, tnon))
     do il_g1=1, gk_f%sizeb_local(1)
       ig1 = mod(gk_f%loc2grow(il_g1) -1, desc_kbz%npw) + 1
       g1 = desc_kbz%gvec(:,ig1)
       ph1 = exp(+j_dpc * two_pi * dot_product(g1, tnon))
       gk_f%buffer_cplx(il_g1, il_g2) = gk_i%buffer_cplx(il_g1, il_g2) * ph1 * ph2
     end do
   end do
   !if (trev_k == 1) then ! TR
   !call gkf%ptrans_ip("C", -cone)
   !call gk%print(header=sjoin("not isirr_k with gt_kibz:", itoa(ik_ibz)))
   end associate
 end do
 !stop

 end associate

10 continue
 !call cwtime_report(" gwr_rotate_gt:", cpu, wall, gflops)

end subroutine gwr_rotate_gt
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_green_gpr
!! NAME
!!  gwr_get_green_gpr
!!
!! FUNCTION
!!  Use FFTs to compute:
!!
!!      G_k(g,g') --> G_k(g',r)
!!
!!  for each k in the BZ treated by this MPI proc for given spin and tau.
!!
!!  1) FFT Transform the first index: G_k(g,g') --> G_k(r,g') and multiply by e^{ik.r}
!!     This is a local operation.
!!
!!  2) MPI transpose the matrix to go from (r,g') to (g',r) distribution.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_green_gpr(gwr, my_it, my_is, desc_kbz, gt_gpr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 integer,intent(in) :: my_it, my_is
 type(desc_t),target,intent(out) :: desc_kbz(gwr%my_nkbz)
 type(matrix_scalapack),intent(inout) :: gt_gpr(2, gwr%my_nkbz)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1
 integer :: my_ikf, ik_bz, ig2, ipm, npwsp, col_bsize
 real(dp) :: cpu, wall, gflops
 real(dp) :: kk_bz(3)
 type(matrix_scalapack) :: rgp
 type(matrix_scalapack),target :: gt_kbz(2)
 complex(dp),allocatable :: ceikr(:)

! *************************************************************************

 !call cwtime(cpu, wall, gflops, "start")

 ABI_MALLOC(ceikr, (gwr%g_nfft * gwr%nspinor))

 do my_ikf=1,gwr%my_nkbz
   ik_bz = gwr%my_kbz_inds(my_ikf)
   kk_bz = gwr%kbz(:, ik_bz)

   call calc_ceikr(kk_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ceikr)

   ! Get G_k in the BZ.
   call gwr%rotate_gt(my_ikf, my_it, my_is, desc_kbz(my_ikf), gt_kbz)

   do ipm=1,2
     associate (ggp => gt_kbz(ipm), desc_k => desc_kbz(my_ikf))

     ! Allocate rgp PBLAS matrix to store G(r, g')
     npwsp = desc_k%npw * gwr%nspinor
     col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
     call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_k%istwfk, &
                   size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

     ABI_CHECK(size(ggp%buffer_cplx, dim=2) == size(rgp%buffer_cplx, dim=2), "len2")

     do ig2=1,ggp%sizeb_local(2)

       ABI_CHECK(size(ggp%buffer_cplx(:, ig2)) == desc_k%npw * gwr%nspinor, "npw")
       ABI_CHECK(size(rgp%buffer_cplx(:, ig2)) == gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

       ! Perform FFT G_k(g,g') -> G_k(r,g') and store results in rgp.
       call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat1, &
                   gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                   ggp%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))

       ! Multiply by e^{ik.r}
       rgp%buffer_cplx(:, ig2) = ceikr * rgp%buffer_cplx(:, ig2)
     end do ! ig2

     ! MPI transpose: G_k(r,g') -> G_k(g',r)
     call rgp%ptrans("N", gt_gpr(ipm, my_ikf))
     call rgp%free()
     end associate
   end do ! ipm

   call slk_array_free(gt_kbz)
 end do ! my_ikf

 ABI_FREE(ceikr)

 !call cwtime_report(" gwr_get_green_gpr:", cpu, wall, gflops)

end subroutine gwr_get_green_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_rotate_wct
!! NAME
!!  gwr_rotate_wct
!!
!! FUNCTION
!!  Reconstruct Wc(q) in the BZ from the IBZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_rotate_wct(gwr, my_iqf, my_it, my_is, desc_qbz, wct_qbz)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 integer,intent(in) :: my_iqf, my_it, my_is
 type(desc_t),intent(out) :: desc_qbz
 type(matrix_scalapack),intent(inout) :: wct_qbz

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, spin, itau, iq_ibz, isym_q, trev_q, g0_q(3), iq_bz
 logical :: isirr_q
 real(dp) :: tsign
!arrays
 integer :: g1(3), g2(3)
 real(dp) :: tnon(3) !, qq_bz(3)
 complex(dp) :: ph2, ph1

! *************************************************************************

 ABI_CHECK(gwr%wc_space == "itau", sjoin("wc_space:", gwr%wc_space, " != itau"))

 spin = gwr%my_spins(my_is)
 itau = gwr%my_itaus(my_it)
 iq_bz = gwr%my_qbz_inds(my_iqf)

 iq_ibz = gwr%my_qbz2ibz(1, my_iqf); isym_q = gwr%my_qbz2ibz(2, my_iqf)
 trev_q = gwr%my_qbz2ibz(6, my_iqf); g0_q = gwr%my_qbz2ibz(3:5, my_iqf)
 isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
 tsign = merge(1, -1, trev_q == 0)
 ! TODO: Understand why legacy GW does not need umklapp
 !ABI_CHECK(all(g0_q == 0), sjoin("For qbz:", ktoa(gwr%qbz(:, iq_bz)), "g0_q:", ltoa(g0_q), " != 0"))

 ! Copy descriptor from IBZ
 associate (desc_qibz => gwr%tchi_desc_qibz(iq_ibz))
 call desc_qibz%copy(desc_qbz)

 if (isirr_q) then
   ! Copy the matrices and we are done.
   call gwr%wc_qibz(iq_ibz, itau, spin)%copy(wct_qbz)
   return
 end if

 ! rotate gvec, recompute gbound and rotate vc_sqrt.
 ! TODO: 1) Handle TR and routine to rotate tchi/W including vc_sqrt
 !       2) Make sure that FFT box is large enough to accomodate umklapps
 do ig1=1,desc_qbz%npw
   desc_qbz%gvec(:,ig1) = matmul(gwr%cryst%symrec(:,:,isym_q), desc_qibz%gvec(:,ig1)) - g0_q
 end do

 call sphereboundary(desc_qbz%gbound, desc_qbz%istwfk, desc_qbz%gvec, gwr%g_mgfft, desc_qbz%npw)

 ! TODO: rotate vc_sqrt
 ! vc(Sq, Sg) = vc(q, g)
 ! vc(-q, -g) = vc(q, g)
 do ig1=1,desc_qbz%npw
   !desc_qbz%vc_sqrt(ig1) = desc_qbz%vc_sqrt(ig1)
 end do

 ! Get Wc_q with q in the BZ.
 tnon = gwr%cryst%tnons(:, isym_q)
 associate (wqi => gwr%wc_qibz(iq_ibz, itau, spin), wqf => wct_qbz)
 call wqi%copy(wct_qbz)
 do il_g2=1, wqf%sizeb_local(2)
   ig2 = mod(wqf%loc2gcol(il_g2) -1, desc_qbz%npw) + 1
   g2 = desc_qbz%gvec(:,ig2)
   ph2 = exp(-j_dpc * two_pi * dot_product(g2, tnon))
   do il_g1=1, wqf%sizeb_local(1)
     ig1 = mod(wqf%loc2grow(il_g1) -1, desc_qbz%npw) + 1
     g1 = desc_qbz%gvec(:,ig1)
     ph1 = exp(+j_dpc * two_pi * dot_product(g1, tnon))
     wqf%buffer_cplx(il_g1, il_g2) = wqi%buffer_cplx(il_g1, il_g2) * ph1 * ph2
   end do
 end do
 !call wqf%print(header=sjoin("not isirr_q with gt_kibz:", itoa(iq_ibz)))
 end associate

 !if (trev_q == 1)
 end associate

end subroutine gwr_rotate_wct
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_wc_gpr
!! NAME
!!  gwr_get_wc_gpr
!!
!! FUNCTION
!!  Use FFTs to compute:
!!
!!      Wc_q(g,g') --> Wc_q(g',r)
!!
!!  for each q in the BZ treated by this MPI procs for given spin and tau.
!!
!!  1) FFT Transform the first index: Wc(g,g',it) --> Wc(r,g',it)  (local operation)
!!
!!  2) MPI transposition: Wc(r,g',it) --> Wc(g',r,it)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_wc_gpr(gwr, my_it, my_is, desc_qbz, wct_gpr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 integer,intent(in) :: my_it, my_is
 type(desc_t),target,intent(out) :: desc_qbz(gwr%my_nqbz)
 type(matrix_scalapack),intent(inout) :: wct_gpr(gwr%my_nqbz)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1
 integer :: my_iqf, ig2, npwsp, col_bsize
 !real(dp) :: cpu, wall, gflops
 type(matrix_scalapack) :: rgp, wct_qbz

! *************************************************************************

 !call cwtime(cpu, wall, gflops, "start")

 do my_iqf=1,gwr%my_nqbz

   ! Get W_q in the BZ.
   call gwr%rotate_wct(my_iqf, my_it, my_is, desc_qbz(my_iqf), wct_qbz)
   associate (desc_q => desc_qbz(my_iqf))

   ! Allocate rgp PBLAS matrix to store Wc(r, g')
   npwsp = desc_q%npw * gwr%nspinor
   col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
   call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_q%istwfk, &
                 size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

   ABI_CHECK(size(wct_qbz%buffer_cplx, dim=2) == size(rgp%buffer_cplx, dim=2), "len2")

   ! FFT. Results stored in rgp
   do ig2=1,wct_qbz%sizeb_local(2)
     ABI_CHECK(size(wct_qbz%buffer_cplx(:, ig2)) == desc_q%npw, "npw")
     ABI_CHECK(size(rgp%buffer_cplx(:, ig2)) == gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

     call fft_ug(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat1, &
                 gwr%g_mgfft, gwr%g_ngfft, desc_q%istwfk, desc_q%gvec, desc_q%gbound, &
                 wct_qbz%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))
   end do ! ig2

   ! MPI Transpose: Wc(r,g') -> Wc(g',r)
   call rgp%ptrans("N", wct_gpr(my_iqf))
   call rgp%free()
   end associate

   call wct_qbz%free()
 end do ! my_iqf

 !call cwtime_report(" gwr_get_wc_gpr:", cpu, wall, gflops)

end subroutine gwr_get_wc_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_cos_transform
!! NAME
!!  gwr_cos_transform
!!
!! FUNCTION
!!  Perform cosine transform.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_cos_transform(gwr, what, mode, sum_spins)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: what, mode
 logical,optional,intent(in) :: sum_spins

!Local variables-------------------------------
!scalars
 integer :: my_iqi, my_is, ig1, ig2, my_it, ierr, iq_ibz, itau, spin, it0, iw
 real(dp) :: cpu, wall, gflops !, tau
 logical :: sum_spins_
!arrays
 real(dp),pointer :: weights(:,:)
 real(dp),allocatable :: my_weights(:,:)
 complex(dp):: my_cwork(gwr%my_ntau), glob_cwork(gwr%ntau)
 type(matrix_scalapack), pointer :: mats(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 sum_spins_ = .False.; if (present(sum_spins)) sum_spins_ = sum_spins

 ! Target weights depending on mode.
 select case(mode)
 case ("iw2t")
   ! From omega to tau
   if (what == "tchi") then
     ABI_CHECK(gwr%tchi_space == "iomega", sjoin("mode:", mode, "with what:", what, "and tchi_space:", gwr%tchi_space))
     gwr%tchi_space = "itau"
   end if
   if (what == "wc") then
     ABI_CHECK(gwr%wc_space == "iomega", sjoin("mode:", mode, "with what:", what, "and wc_space:", gwr%wc_space))
     gwr%wc_space = "itau"
   end if
   weights => gwr%w2t_cos_wgs

 case ("it2w")
   ! From tau to omega
   if (what == "tchi") then
     ABI_CHECK(gwr%tchi_space == "itau", sjoin("mode:", mode, " with what:", what, "and tchi_space:", gwr%tchi_space))
     gwr%tchi_space = "iomega"
   end if
   if (what == "wc") then
     ABI_CHECK(gwr%wc_space == "itau", sjoin("mode:", mode, " with what:", what, "and wc_space:", gwr%wc_space))
     gwr%wc_space = "iomega"
   end if
   weights => gwr%t2w_cos_wgs

 case default
   ABI_ERROR(sjoin("Wrong mode:", mode))
 end select

 ! Extract my weights.
 ABI_MALLOC(my_weights, (gwr%ntau, gwr%my_ntau))

 do my_it=1,gwr%my_ntau
   itau = gwr%my_itaus(my_it)
   do iw=1,gwr%ntau
     my_weights(iw, my_it) = weights(iw, itau)
   end do
 end do

 ! And now perform inhomegenous FFT in parallel.
 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))

     mats => null()
     if (what == "tchi") mats => gwr%tchi_qibz(iq_ibz, :, spin)
     if (what =="wc")   mats => gwr%wc_qibz(iq_ibz, :, spin)
     ABI_CHECK(associated(mats), sjoin("Invalid value for what:", what))

     it0 = gwr%my_itaus(1)
     do ig2=1,mats(it0)%sizeb_local(2)
       do ig1=1,mats(it0)%sizeb_local(1)

         ! Extract matrix elements as a function of tau
         ! TODO: Here we can block over ig1 and call zgemv to reduce the number of MPI communications.
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           my_cwork(my_it) = mats(itau)%buffer_cplx(ig1, ig2)
         end do

         do itau=1,gwr%ntau
           glob_cwork(itau) = dot_product(my_weights(itau, :), my_cwork)
         end do

         call xmpi_sum(glob_cwork, gwr%tau_comm%value, ierr)

         ! Update my local (g1, g2) entry to have it in imaginary-frequency space.
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           mats(itau)%buffer_cplx(ig1, ig2) = glob_cwork(itau)
         end do

       end do ! ig1
     end do ! ig2

     end associate
   end do ! my_iqi
 end do ! my_is

 if (sum_spins_ .and. gwr%nspinor /= 2) then  ! gwr%nsppol == 2 .and.
   ! Sum over spins
   do my_iqi=1,gwr%my_nqibz
      iq_ibz = gwr%my_qibz_inds(my_iqi)
      do my_is=1,gwr%my_nspins
        spin = gwr%my_spins(my_is)
        mats => null()
        if (what == "tchi") mats => gwr%tchi_qibz(iq_ibz,:,spin)
        !if (what =="wc")   mats => gwr%wc_qibz(iq_ibz, :, spin)
        ABI_CHECK(associated(mats), sjoin("Invalid value for what:", what))

        do my_it=1,gwr%my_ntau
          itau = gwr%my_itaus(my_it)

          if (gwr%nsppol == 1) then
            mats(itau)%buffer_cplx = two * mats(itau)%buffer_cplx

          else if (gwr%nsppol == 2) then
            if (gwr%spin_comm%nproc > 1) then
              ! Spins are distributed and we have to sum them.
              call xmpi_sum(mats(itau)%buffer_cplx, gwr%spin_comm%value, ierr)
            else
              ! Spins are not distributed (this should happen only in sequential execution).
              if (spin == 1) then
                mats(itau)%buffer_cplx = mats(itau)%buffer_cplx + gwr%tchi_qibz(iq_ibz,itau,spin+1)%buffer_cplx
              end if
            end if
          end if

        end do ! my_it
      end do ! my_is
   end do ! my_iqi
 end if

 ABI_FREE(my_weights)

 call cwtime_report(" gwr_cos_transform:", cpu, wall, gflops)

end subroutine gwr_cos_transform
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_copy
!! NAME
!!  desc_copy
!!
!! FUNCTION
!!  Copy object
!!  NB: cannot use obj1 = obj2 syntax because ABINIT memory-leak detector
!!  won't see the allocation performed by the compiler.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine desc_copy(in_desc, new_desc)

!Arguments ------------------------------------
 class(desc_t),intent(in) :: in_desc
 class(desc_t),intent(out) :: new_desc

! *************************************************************************

 new_desc%istwfk = in_desc%istwfk
 new_desc%npw = in_desc%npw

 call alloc_copy(in_desc%gvec, new_desc%gvec)
 call alloc_copy(in_desc%gbound, new_desc%gbound)
 if (allocated(in_desc%vc_sqrt)) call alloc_copy(in_desc%vc_sqrt, new_desc%vc_sqrt)

end subroutine desc_copy
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_calc_gnorm_table
!! NAME
!!  desc_calc_gnorm_table
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine desc_calc_gnorm_table(desc, cryst)

 use m_sort, only : sort_dp

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc
 class(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
 integer :: ig, isort
 integer, allocatable :: gsort2gvec(:)

! *************************************************************************

 if (allocated(desc%sorted_gnorm)) return

 ABI_MALLOC(desc%sorted_gnorm, (desc%npw))
 ABI_MALLOC(desc%gvec2gsort, (desc%npw))
 ABI_MALLOC(gsort2gvec, (desc%npw))

 do ig=1,desc%npw
   desc%sorted_gnorm(ig) = normv(desc%gvec(:,ig), cryst%gmet, "G")
   gsort2gvec(ig) = ig
 end do

 call sort_dp(desc%npw, desc%sorted_gnorm, gsort2gvec, tol16)

 ! Invert the mapping.
 do isort=1,desc%npw
   ig = gsort2gvec(isort)
   desc%gvec2gsort(ig) = isort
 end do

 ABI_FREE(gsort2gvec)

end subroutine desc_calc_gnorm_table
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_free
!! NAME
!!  desc_free
!!
!! FUNCTION
!!  Free memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine desc_free(desc)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc

! *************************************************************************

 ABI_SFREE(desc%gvec)
 ABI_SFREE(desc%gbound)
 ABI_SFREE(desc%vc_sqrt)
 ABI_SFREE(desc%gvec2gsort)
 ABI_SFREE(desc%sorted_gnorm)

end subroutine desc_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/print_gsort_
!! NAME
!!  print_gsort_
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_gsort_(mat, desc, cryst)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: mat
 class(desc_t),intent(inout) :: desc
 class(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
 integer :: il_g1, il_g2, ig1, ig2, iglob1, iglob2, is1, is2
 complex(dp),allocatable :: cwork(:,:)

! *************************************************************************

 call desc%calc_gnorm_table(cryst)

 ! FIXME: It won't work if nspinor == 2
 ABI_MALLOC(cwork, (mat%sizeb_local(1), mat%sizeb_local(2)))
 cwork = huge(one)

 do il_g2=1,mat%sizeb_local(2)
   iglob2 = mat%loc2gcol(il_g2)
   ig2 = mod(iglob2 -1 , desc%npw) + 1
   is2 = desc%gvec2gsort(ig2)
   do il_g1=1,mat%sizeb_local(1)
     iglob1 = mat%loc2grow(il_g1)
     ig1 = mod(iglob1 - 1, desc%npw) + 1
     is1 = desc%gvec2gsort(ig1)
     !print *, "is1, is2", is1, is2
     cwork(is1, is2) = mat%buffer_cplx(il_g1, il_g2)
   end do
 end do
 !stop

 call print_arr(cwork)
 ABI_FREE(cwork)

end subroutine print_gsort_
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_print
!! NAME
!!  gwr_print
!!
!! FUNCTION
!!  Print info on the gwr object.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_print(gwr, header, unit)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 character(len=*),optional,intent(in) :: header
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
!scalars
 integer :: my_is, spin, my_it, itau, my_iki, ik_ibz, unt, ii
 character(len=500) :: msg
 type(yamldoc_t) :: ydoc

! *********************************************************************

 unt = std_out; if (present(unit)) unt =unit

 msg = ' ==== Info on the GWR object ==== '
 if (present(header)) msg=' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(unt, msg)

 ydoc = yamldoc_open('GWR_params') !, width=11, real_fmt='(3f8.3)')
 call ydoc%add_string("gwr_task", gwr%dtset%gwr_task)
 call ydoc%add_int("ntau", gwr%ntau)
 call ydoc%add_int1d("ngkpt", gwr%ngkpt)
 call ydoc%add_int("nkbz", gwr%nkbz)
 call ydoc%add_int("nkibz", gwr%nkibz)
 call ydoc%add_int1d("ngqpt", gwr%ngqpt)
 call ydoc%add_int("nqbz", gwr%nqbz)
 call ydoc%add_int("nqibz", gwr%nqibz)
 call ydoc%add_int("green_mpw", gwr%green_mpw)
 call ydoc%add_int("tchi_mpw", gwr%tchi_mpw)
 call ydoc%add_int1d("g_ngfft", gwr%g_ngfft(1:8))
 call ydoc%add_int1d("P gwr_np_gtks", [gwr%g_comm%nproc, gwr%tau_comm%nproc, gwr%kpt_comm%nproc, gwr%spin_comm%nproc])
 call ydoc%add_int1d("P np_kibz", gwr%np_kibz)
 call ydoc%add_int1d("P np_qibz", gwr%np_qibz)
 ! Print Max error due to inhomogeneous FT.
 call ydoc%add_real("ft_max_err_t2w_cos", gwr%ft_max_error(1))
 call ydoc%add_real("ft_max_err_w2t_cos", gwr%ft_max_error(2))
 call ydoc%add_real("ft_max_err_t2w_sin", gwr%ft_max_error(3))
 call ydoc%add_real("cosft_duality_err", gwr%cosft_duality_err)

 ! Print imaginary time/frequency mesh with weights.
 call ydoc%open_tabular("Minimax imaginary tau/omega mesh", comment="tau, weight(tau), omega, weight(omega)")
 do ii=1,gwr%ntau
   write(msg, "(i0, 4(es12.5,2x))")ii, gwr%tau_mesh(ii), gwr%tau_wgs(ii), gwr%iw_mesh(ii), gwr%iw_wgs(ii)
   call ydoc%add_tabular_line(msg)
 end do

 call ydoc%write_and_free(unt)

 if (gwr%dtset%prtvol > 10) then
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       do my_iki=1,gwr%my_nkibz
         ik_ibz = gwr%my_kibz_inds(my_iki)
         associate (gt => gwr%gt_kibz(:, ik_ibz, itau, spin))
         call gt(1)%print(header=sjoin("gt_kibz(1) for ik_ibz:", itoa(ik_ibz)), unit=unt)
         call gt(2)%print(header=sjoin("gt_kibz(2) for ik_ibz:", itoa(ik_ibz)), unit=unt)
         end associate
       end do
     end do
   end do
 end if

end subroutine gwr_print
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_tchi
!! NAME
!!  gwr_build_cthi
!!
!! FUNCTION
!!  Compute irreducible polarizability.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_tchi(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1, istwfk1 = 1
 integer :: my_is, my_it, my_ikf, ig, my_ir, my_nr, npwsp, ncol_glob, col_bsize, my_iqi !, my_iki ! my_iqf,
 integer :: sc_nfft, spin, ik_ibz, ik_bz, iq_ibz, ierr, ipm, itau, ig2 ! ig1
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all, spin_fact
 character(len=5000) :: msg
 type(desc_t),pointer :: desc_k, desc_q
 type(matrix_scalapack) :: chi_rgp
!arrays
 integer :: sc_ngfft(18)
 integer,allocatable :: green_scg(:,:), chi_scg(:,:)
 real(dp) :: kk_bz(3), kpq_bz(3), qq_ibz(3) ! kk_ibz(3), qq_bz(3)
 complex(dp),allocatable :: gt_scbox(:,:), chit_scbox(:)
 type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz), chiq_gpr(gwr%my_nqibz)
 type(desc_t), target :: desc_kbz(gwr%my_nkbz)
 type(fftbox_plan3_t) :: plan_gp2rp, plan_rp2gp
 complex(dp),allocatable :: cemiqr(:)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ABI_CHECK(gwr%tchi_space == "none", sjoin("tchi_space: ", gwr%tchi_space, " != none"))
 gwr%tchi_space = "itau"

 spin_fact = two / (gwr%nsppol * gwr%nspinor)

 if (gwr%use_supercell) then
   ! ============================
   ! GWr algorithm with supercell
   ! ============================

   ! Set FFT mesh in the supercell
   ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
   ! Perhaps the safest approach would be to generate the plan on the fly.
   sc_ngfft = gwr%g_ngfft
   sc_ngfft(1:3) = gwr%ngkpt * gwr%g_ngfft(1:3)
   sc_ngfft(4:6) = sc_ngfft(1:3)
   sc_nfft = product(sc_ngfft(1:3))
   !sc_augsize = product(sc_ngfft(4:6))
   call wrtout(std_out, sjoin(" Building chi0 = iG0G0 with FFTs in the supercell:", ltoa(sc_ngfft(1:3))))

   ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat1, 2))
   ABI_CALLOC(chit_scbox, (sc_nfft * gwr%nspinor * ndat1))

   ! This should be the correct sequence
   ! (ngfft, ndat, isign)
   call plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1 * 2, -1)
   call plan_rp2gp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1,     +1)

   !call plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1 * 2, +1)
   !call plan_rp2gp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1,     -1)

   ! The g-vectors in the supercell for G and tchi.
   ABI_MALLOC(green_scg, (3, gwr%green_mpw))
   ABI_MALLOC(chi_scg, (3, gwr%tchi_mpw))
   ! The phase e^{-iq.r}
   ABI_MALLOC(cemiqr, (gwr%g_nfft * gwr%nspinor))

   ! Allocate PBLAS arrays for tchi_q(g',r) for all q in the IBZ treated by this MPI rank.
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     npwsp = gwr%tchi_desc_qibz(iq_ibz)%npw * gwr%nspinor
     ncol_glob = gwr%g_nfft * gwr%nspinor
     col_bsize = ncol_glob / gwr%g_comm%nproc; if (mod(ncol_glob, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
     call chiq_gpr(my_iqi)%init(npwsp, gwr%g_nfft * gwr%nspinor, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
   end do

   ! Loop over my spins and my taus.
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
       itau = gwr%my_itaus(my_it)

       ! G_k(g,g') --> G_k(g',r) e^{ik.r} for each k in the BZ treated by me.
       call gwr%get_green_gpr(my_it, my_is, desc_kbz, gt_gpr)

       ! Loop over the second index r that is MPI-distributed in g_comm.
       my_nr = gt_gpr(1, 1)%sizeb_local(2)
       do my_ir=1,my_nr

         ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
         ! Note that we need to take the union of (k, g') for k in the BZ.
         gt_scbox = zero
         do my_ikf=1,gwr%my_nkbz
           ik_bz = gwr%my_kbz_inds(my_ikf)
           ik_ibz = gwr%my_kbz2ibz(1, my_ikf)

           ! Compute k+g
           desc_k => desc_kbz(my_ikf)
           do ig=1,desc_k%npw
             green_scg(:,ig) = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ig)
           end do

           do ipm=1,2
             call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat1, green_scg, &
                           gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                           gt_scbox(:, ipm))                             ! inout
           end do
         end do ! my_ikf

         ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
         call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)

         ! G(G',r) --> G(R',r)
         call plan_gp2rp%execute_ip_dpc(gt_scbox)
         !gt_scbox = gt_scbox * sc_nfft / gwr%cryst%ucvol

         ! Compute tchi(R',r) for this r
         !chit_scbox = -j_dpc * gt_scbox(:, 1) * gt_scbox(:, 2) * spin_fact
         chit_scbox = gt_scbox(:, 1) * conjg(gt_scbox(:, 2)) ! * spin_fact

         ! Back to tchi(G'=q+g',r) space immediately.
         call plan_rp2gp%execute(chit_scbox)

         ! Extract tchi_q(g', r) on the g-sphere from the FFT box in the supercell.
         ! not necessarly equal to the one used for the Green's function.
         ! NB: For the time being I'm assuming k-mesh == q-mesh.
         ! Only q-points in the IBZ are considered.
         do my_iqi=1,gwr%my_nqibz
           iq_ibz = gwr%my_qibz_inds(my_iqi)
           qq_ibz = gwr%qibz(:, iq_ibz)
           desc_q => gwr%tchi_desc_qibz(iq_ibz)

           ! Compute q + g
           do ig=1,desc_q%npw
             chi_scg(:,ig) = nint(qq_ibz * gwr%ngqpt) + gwr%ngqpt * desc_q%gvec(:,ig)
           end do

           call box2gsph(sc_ngfft, desc_q%npw, gwr%nspinor * ndat1, chi_scg, &
                         chit_scbox, &                            ! in
                         chiq_gpr(my_iqi)%buffer_cplx(:, my_ir))  ! out
         end do ! my_iqi
       end do ! my_ir

       ! Free descriptors and PBLAS matrices in kBZ.
       call desc_array_free(desc_kbz)
       call slk_array_free(gt_gpr)

       ! Now we have tchi_q(g',r).
       ! For each IBZ q-point treated by this MPI proc, do:
       !
       !     1) MPI transpose to have tchi_q(r,g')
       !     2) FFT along first dimension to have tchi_q(g,g') stored in gwr%tchi_qibz
       !
       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         desc_q => gwr%tchi_desc_qibz(iq_ibz)

         ! Note minus sign in q
         call calc_ceikr(-gwr%qibz(:,iq_ibz), gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, cemiqr)

         ! MPI-transposition
         call chiq_gpr(my_iqi)%ptrans("N", chi_rgp)

         ABI_CHECK(size(gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx, dim=2) == size(chi_rgp%buffer_cplx, dim=2), "len2")

         ! FFT r --> g along the first dimension: tchi_q(r,g') --> tchi_q(g,g')
         do ig2=1,chi_rgp%sizeb_local(2)

           ABI_CHECK(size(chi_rgp%buffer_cplx(:, ig2)) == gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")
           ABI_CHECK(size(gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2)) == desc_q%npw, "npw")

           call fft_ur(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat1, gwr%g_mgfft, gwr%g_ngfft, &
                       istwfk1, desc_q%gvec, desc_q%gbound, &
                       chi_rgp%buffer_cplx(:, ig2),         &                  ! ur(in)
                       gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2))  ! ug(out)

           gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) = &
           gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) * cemiqr  * gwr%cryst%ucvol
           !gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) * cemiqr  ! * gwr%g_nfft / gwr%cryst%ucvol
         end do
         call chi_rgp%free()
       end do

       write(msg,'(2(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "]"
       call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
     end do ! my_it
   end do ! spin

   ABI_FREE(gt_scbox)
   ABI_FREE(chit_scbox)
   ABI_FREE(green_scg)
   ABI_FREE(chi_scg)
   ABI_FREE(cemiqr)
   call slk_array_free(chiq_gpr)

  else
    ! ===================================================================
    ! Mixed-space algorithm in the unit cell with convolutions in k-space
    ! ===================================================================

    do my_is=1,gwr%my_nspins
      spin = gwr%my_spins(my_is)
      do my_it=1,gwr%my_ntau
        call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
        itau = gwr%my_itaus(my_it)

        ! Redistribute G_k in the IBZ so that each MPI proc can recostruct G in the BZ on the fly.
        !need_kibz(gwr%nkibz)
        !got_kibz(gwr%nkibz)
        !call gwr%green_redistribute(my_it, my_is, need_kibz, got_kibz, "allocate')
        !call gwr%green_redistribute(my_it, my_is, need_kibz, got_kibz, "free")

        do my_ikf=1,gwr%my_nkbz
          ik_bz = gwr%my_kbz_inds(my_ikf)
          kk_bz = gwr%kbz(:, ik_bz)

          ! Use symmetries to get G_k from the corresponding image in the IBZ.
          ! G_k(g,g') --> G_k(r,r')
          !call gwr%gk_rrp(my_ikf, my_it, my_is, "C", gk_rrp)
          !call gk_rrp%free()

          do my_iqi=1,gwr%my_nqibz
            iq_ibz = gwr%my_qibz_inds(my_iqi)
            qq_ibz = gwr%qibz(:, iq_ibz)
            !desc_q => gwr%tchi_desc_qibz(iq_ibz)
            !iq_bz = gwr%my_qbz_inds(my_iqf)
            !qq_bz = gwr%qbz(:, iq_bz)

            kpq_bz = kk_bz + qq_ibz
            !ikpq_ibz = ??

            ! Use symmetries to get G_kq from the corresponding image in the IBZ.
            ! G_kq(g,g') --> G_kq(r,r')
            !call gwr%gk_rrp(my_ikqf, my_it, my_is, "N", gkq_rrp)
            !call gkq_rrp%free()

            !tchi_rrp(my_iqi)%buffer_cplx = tchi(my_iqi)%buffer_cplx + gk_rrp%buffer_cplx * gkq_rrp%buffer_cplx

            !desc_k = gwr%green_desc_kibz(ik_ibz)%rotate(symtab_k)
            !desc_kpq = gwr%green_desc_kibz(ikpq_ibz)%rotate(symtab_kpq)

            !gwr%gt_kibz(1, ik_ibz, itau, spin)%buffer_cplx(:,ig2)

            !do ig2=1,ggp%sizeb_local(2)
            !  call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat1, &
            !              gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
            !              ggp%buffer_cplx(:, ig2),
            !              rgp%buffer_cplx(:, ig2))
            !end do ! ig2
            !call rgp%ptrans("C", gpr)
            !call rgp%free()

            !do ir1=1,ggp%sizeb_local(2)
            !  call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat1, &
            !              gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
            !              gpr%buffer_cplx(:, ir1),
            !              rpr%buffer_cplx(:, ir1))
            !end do ! ir1
            !call gpr%free()
            !call rpr%ptrans("C", rrp)
            !call rpr%free()

            !call desc_k%free()
            !call desc_kpq%free()

            !my_nr = go_k_gpr%sizeb_local(2)
            !do my_ir=1,my_nr
            !   !chi_fft = go_rpr * ge_rpr
            !   ! from r' to g'
            !   !chi_rrp%buffer_cplx(:, my_ir) = gtk_rrp(1)%buffer_cplx(:, my_ir) * gtkq_rrp(2)%buffer_cplx(:, my_ir)
            !end do

            !call chi_gp_r%ptrans("C", chi_r_gp)
            !call chi_gp_r%free()
            !do ig2=1, chi_rp_g%sizeb_local(2)
              !chi_r_gp%buffer_cplx(:, ig2)
              !!gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) = FFT_rg
            !end do
            !call chi_r_gp%free()
          end do ! my_ikf
        end do ! iq_ibz

       write(msg,'(2(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "]"
       call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
     end do ! my_it
   end do ! spin
 end if

 ! Print trace of chi_q(i tau) matrices for testing purposes.
 if (gwr%dtset%prtvol > 0) call gwr%print_trace("tchi_qibz")

 ! Transform irreducible tchi from imaginary tau to imaginary omega.
 ! Also sum over spins to get total tchi if collinear spin.
 call gwr%cos_transform("tchi", "it2w", sum_spins=.True.)

 ! Print trace of chi_q(i omega) matrices for testing purposes.
 if (gwr%dtset%prtvol > 0) call gwr%print_trace("tchi_qibz")

 ! Write file with chi0(i omega)
 if (gwr%dtset%prtsuscep > 0) then
   call gwr%ncwrite_tchi_wc("tchi", trim(gwr%dtfil%filnam_ds(4))//'_TCHI.nc')
 end if

 ! When reading, we have to perform consistency checks handle possible difference in:
 !   - input ecuteps and the value found on disk
 ! Changing q-mesh or time/imaginary mesh is not possible.

 !call gwr%ncread_tchi_wc("tchi", trim(gwr%dtfil%filnam_ds(4))//'_TCHI.nc'))

 call cwtime_report(" gwr_build_tchi:", cpu_all, wall_all, gflops_all)

end subroutine gwr_build_tchi
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_print_trace
!! NAME
!!  gwr_print_trace
!!
!! FUNCTION
!!  Print traces of PBLAS matrices to std_out and ab_out.
!!  NB: This is a global routine that should be called by all procs inside gwr%comm.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_print_trace(gwr, what)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: what

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: my_is, spin, my_it, itau, iq_ibz, ierr, my_iqi, my_iki, ik_ibz, ipm, my_rank
 character(len=5000) :: comment !, msg
 integer :: units(2)
 complex(dp),allocatable :: ctrace3(:,:,:), ctrace4(:,:,:,:)
 type(matrix_scalapack),pointer :: mats(:,:,:)

! *************************************************************************

 ! The same q/k point in the IBZ might be available on different procs in kpt_comm
 ! thus we have to rescale the trace before summing the results in gwr%comm.

 my_rank = xmpi_comm_rank(gwr%comm)

 comment = "Invalid space!"
 units = [std_out, ab_out]

 select case (what)
 case ("tchi_qibz", "wc_qibz")
   ! Trace of tchi or Wc
   ABI_CALLOC(ctrace3, (gwr%nqibz, gwr%ntau, gwr%nsppol))

   if (what == "tchi_qibz") then
     mats => gwr%tchi_qibz
     if (gwr%tchi_space == "iomega") comment = " (iq_ibz, iomega) table"
     if (gwr%tchi_space == "itau") comment = " (iq_ibz, itau) table"
   else if (what == "wc_qibz") then
     mats => gwr%wc_qibz
     if (gwr%wc_space == "iomega") comment = " (iq_ibz, iomega) table"
     if (gwr%wc_space == "itau") comment = " (iq_ibz, itau) table"
   end if

   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         ctrace3(iq_ibz, itau, spin) = mats(iq_ibz, itau, spin)%get_trace() / gwr%np_qibz(iq_ibz)
!DEBUG
         if (xmpi_comm_size(gwr%comm) == 1 .and. my_rank == master &
             .and. what == "tchi_qibz" .and. gwr%tchi_space == "iomega") then
           if (itau == 1) then
             call wrtout(std_out, sjoin(what, " in iomega space for w:", ftoa(gwr%iw_mesh(itau)*Ha_eV), "(eV)"))
             call print_arr(mats(iq_ibz, itau, spin)%buffer_cplx)
             call wrtout(std_out, " Printing sorted GG' matrix")
             call print_gsort_(mats(iq_ibz, itau, spin), gwr%tchi_desc_qibz(iq_ibz), gwr%cryst)
             !call xmpi_abort(msg="Aborting now")
           end if
         end if
!END DEBUG
       end do
     end do
   end do

   call xmpi_sum_master(ctrace3, 0, gwr%tks_comm%value, ierr)

   if (xmpi_comm_rank(gwr%comm) == master) then
     do spin=1,gwr%nsppol
       call wrtout(units, sjoin(" Trace of:", what, "for spin:", itoa(spin), "for testing purposes:"))
       call wrtout(units, comment, newlines=1)
       call print_arr(ctrace3(:,:,spin), unit=ab_out)
       call print_arr(ctrace3(:,:,spin), unit=std_out)
     end do
   end if
   ABI_FREE(ctrace3)

 case ("gt_kibz")
   ! Trace of Green's functions.
   ABI_CALLOC(ctrace4, (gwr%nkibz, gwr%ntau, 2, gwr%nsppol))

   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       do my_iki=1,gwr%my_nkibz
         ik_ibz = gwr%my_kibz_inds(my_iki)
         do ipm=1,2
           ctrace4(ik_ibz, itau, ipm, spin) = gwr%gt_kibz(ipm, ik_ibz, itau, spin)%get_trace() / gwr%np_kibz(ik_ibz)
         end do
       end do
     end do
   end do
   comment = "(ik_ibz, itau) table"

   call xmpi_sum_master(ctrace4, master, gwr%tks_comm%value, ierr)

   if (xmpi_comm_rank(gwr%comm) == master) then
     do spin=1,gwr%nsppol
       do ipm=1,2
         call wrtout(units, sjoin(" Trace of:", what, "for ipm:", itoa(ipm), ", spin:", itoa(spin), "for testing purposes:"))
         call wrtout(units, comment, newlines=1)
         call print_arr(ctrace4(:,:,ipm, spin), unit=ab_out)
         call print_arr(ctrace4(:,:,ipm, spin), unit=std_out)
       end do
     end do
   end if
   ABI_FREE(ctrace4)

 case default
   ABI_ERROR(sjoin("Invalid value of what:", what))
 end select

end subroutine gwr_print_trace
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_wc
!! NAME
!!  gwr_build_wc
!!
!! FUNCTION
!!  Compute Wc(i tau) from tchi(i omega)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_wc(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_iqi, my_it, my_is, iq_ibz, spin, itau, iw, ierr, il_g1, il_g2, ig1, ig2, iglob1, iglob2
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_q, wall_q, gflops_q
 logical :: q_is_gamma, free_chi
 character(len=5000) :: msg
 complex(dpc) :: vcs_g1, vcs_g2
 type(matrix_scalapack) :: tmp_mat
 type(yamldoc_t) :: ydoc
 !arrays
 real(dp) :: qq_ibz(3)
 complex(dpc) :: em1_wq(gwr%ntau, gwr%nqibz), eps_wq(gwr%ntau, gwr%nqibz)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call wrtout(std_out, " Building correlated screening Wc from irreducible chi ...")

 ABI_CHECK(gwr%tchi_space == "iomega", sjoin("tchi_space: ", gwr%tchi_space, " != iomega"))
 ABI_CHECK(gwr%wc_space == "none", sjoin("wc_space: ", gwr%wc_space, " != none"))
 gwr%wc_space = "iomega"

 ! =======================================
 ! Allocate PBLAS arrays for wc_qibz(g,g')
 ! =======================================
 ! Note that we have already summed tchi over spin.
 ABI_MALLOC(gwr%wc_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))
 free_chi = .False.

 em1_wq = zero; eps_wq = zero

 do my_iqi=1,gwr%my_nqibz
   call cwtime(cpu_q, wall_q, gflops_q, "start")
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   q_is_gamma = normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0
   associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))

   ! The spin loop is needed so that procs in different pools operate
   ! on their own matrix that has been already summed over (collinear) spins.
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       associate (wc => gwr%wc_qibz(iq_ibz, itau, spin))
       call gwr%tchi_qibz(iq_ibz, itau, spin)%copy(wc, free=free_chi)
       !call wc%print(header="wc")

       ! Build epsilon(q,iw) = delta_{g,g'} - v_q(g,g') tchi_q(g,g,iw).

       ! RPA: \tepsilon = 1 - Vc^{1/2} chi0 Vc^{1/2}
       ! vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff.

       !wc%buffer_cplx = zero
       do il_g2=1,wc%sizeb_local(2)
         iglob2 = wc%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,wc%sizeb_local(1)
           iglob1 = wc%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1
           vcs_g1 = desc_q%vc_sqrt(ig1)
           wc%buffer_cplx(il_g1, il_g2) = - wc%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2
           if (iglob1 == 1 .and. iglob2 == 1) then
             wc%buffer_cplx(il_g1, il_g2) = one + wc%buffer_cplx(il_g1, il_g2)
             ! Store epsilon_{iw, iq_ibz}(0, 0)
             eps_wq(itau, iq_ibz) = wc%buffer_cplx(il_g1, il_g2) / gwr%np_qibz(iq_ibz)
           end if
         end do ! il_g1
       end do ! il_g2

       ! Invert tilde epsilon.
       ! NB: PZGETRF requires square block cyclic decomposition i.e., MB_A = NB_A.
       ! hence we neeed to redistribute the data before calling zinvert.
       ! TODO: Can call zhdp
       !call wrtout(std_out, "Printing wc%buffer_cplex before inversion")
       !call print_arr(wc%buffer_cplx, unit=std_out)

       call wc%change_size_blocs(tmp_mat)
       call tmp_mat%zinvert()
       !call tmp_mat%zhdp_invert()
       call wc%take_from(tmp_mat)
       call tmp_mat%free()

       !call wrtout(std_out, "Printing wc%buffer_cplex after inversion")
       !call print_arr(wc%buffer_cplx, unit=std_out)

       ! Build W(q, iw) = e^{-1}_q(g,g',iw) v_q(g,g') and remove bare vc
       do il_g2=1,wc%sizeb_local(2)
         iglob2 = wc%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,wc%sizeb_local(1)
           iglob1 = wc%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1
           if (iglob1 == 1 .and. iglob2 == 1) then
             ! Store epsilon^{-1}_{iw, iq_ibz}(0, 0)
             em1_wq(itau, iq_ibz) = wc%buffer_cplx(il_g1, il_g2) / gwr%np_qibz(iq_ibz)
           end if
           vcs_g1 = desc_q%vc_sqrt(ig1)
           wc%buffer_cplx(il_g1, il_g2) = wc%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2
           if (iglob1 == 1 .and. iglob2 == 1) then
             ! Subtract exchange part.
             wc%buffer_cplx(il_g1, il_g2) = wc%buffer_cplx(il_g1, il_g2) - cone
           end if
         end do ! il_g1
       end do ! il_g2

       end associate
     end do  ! my_it
   end do ! my_is
   end associate

   write(msg,'(2(a,i0),a)')" My iqi [", my_iqi, "/", gwr%my_nqibz, "]"
   call cwtime_report(msg, cpu_q, wall_q, gflops_q)
 end do ! my_iqi

 call xmpi_sum_master(em1_wq, master, gwr%gtk_comm%value, ierr)
 call xmpi_sum_master(eps_wq, master, gwr%gtk_comm%value, ierr)

 if (xmpi_comm_rank(gwr%comm) == master) then
   ! TODO: Write data to file.
   ydoc = yamldoc_open('EMACRO_WITHOUT_LOCAL_FIELDS') !, width=11, real_fmt='(3f8.3)')
   call ydoc%open_tabular("epsilon_{iw, q -> Gamma}(0,0)") ! comment="(iomega, iq_ibz)")
   do iw=1,gwr%ntau
     write(msg, "(3(es16.8,2x))") gwr%iw_mesh(iw), real(eps_wq(iw, 1)), aimag(eps_wq(iw, 1))
     call ydoc%add_tabular_line(msg)
   end do
   call ydoc%write_units_and_free([ab_out, std_out])

   ydoc = yamldoc_open('EMACRO_WITH_LOCAL_FIELDS') !, width=11, real_fmt='(3f8.3)')
   call ydoc%open_tabular("epsilon_{iw, q -> Gamma}(0,0)") !, comment="(iomega, iq_ibz)")
   do iw=1,gwr%ntau
     write(msg, "(3(es16.8,2x))") gwr%iw_mesh(iw), real(em1_wq(iw, 1)), aimag(em1_wq(iw, 1))
     call ydoc%add_tabular_line(msg)
   end do
   call ydoc%write_units_and_free([ab_out, std_out])
 end if

 ! Print trace of wc_q(itau) matrices for testing purposes.
 !if (gwr%dtset%prtvol > 0) call gwr%print_trace("wc_qibz")

 ! Cosine transform from iomega to itau to get Wc(i tau)
 call gwr%cos_transform("wc", "iw2t")

 ! Print trace of wc_q(iomega) matrices for testing purposes.
 !if (gwr%dtset%prtvol > 0) call gwr%print_trace("wc_qibz")

 call cwtime_report(" gwr_build_wc:", cpu_all, wall_all, gflops_all)

end subroutine gwr_build_wc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_sigmac
!! NAME
!!  gwr_build_sigmac
!!
!! FUNCTION
!!  Build Sigma_c(i tau) and compute matrix elements in the KS basis set.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_sigmac(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer, parameter :: ndat1 = 1
 integer :: my_is, my_it, spin, ik_ibz, sc_nfft, my_ir, my_nr !my_iki,
 integer :: my_iqf, iq_ibz, iq_bz, itau, ierr, jb, ib1, ib2, bmin, bmax ! col_bsize, npwsp,
 integer :: my_ikf, ipm, ik_bz, ig, ikcalc ! my_iqi,
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all !, spin_fact
 !logical :: q_is_gamma
 type(desc_t), pointer :: desc_q, desc_k
 character(len=500) :: msg
!arrays
 integer :: sc_ngfft(18)
 integer,allocatable :: green_scg(:,:), wc_scg(:,:)
 !real(dp) :: vals(2) !, qq_ibz(3) ! kk_ibz(3),
 complex(dp),allocatable :: cit_me(:,:,:,:,:), ciw_me(:,:,:,:,:)
 complex(dp),allocatable :: gt_scbox(:,:), wct_scbox(:)
 complex(gwpc),allocatable :: ur_bk(:,:,:)
 type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz), wc_gpr(gwr%my_nqbz)
 type(desc_t), target :: desc_kbz(gwr%my_nkbz), desc_qbz(gwr%my_nqbz)
 type(fftbox_plan3_t) :: gt_plan_gp2rp, wt_plan_gp2rp

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ABI_CHECK(gwr%wc_space == "itau", sjoin("wc_space: ", gwr%wc_space, " != itau"))

 ! ===========================================
 ! Allocate PBLAS arrays for sigmac_kibz(g,g')
 ! ===========================================
 ! TODO: This is not needed. Perhaps for self-consistent without KS/HF representation.
 !ABI_MALLOC(gwr%sigc_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))
 !gwr%sigc_space = "itau"

 !do my_is=1,gwr%my_nspins
 !  spin = gwr%my_spins(my_is)
 !  do my_it=1,gwr%my_ntau
 !    itau = gwr%my_itaus(my_it)

 !    do my_iki=1,gwr%my_nkibz
 !      ik_ibz = gwr%my_kibz_inds(my_iki)
 !      npwsp = gwr%green_desc_kibz(ik_ibz)%npw * gwr%nspinor
 !      col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
 !      associate (sigc => gwr%sigc_kibz(:, ik_ibz, itau, spin))
 !      call sigc(1)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
 !      call sigc(2)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
 !      !call sigc(1)%print(header=sjoin("sigc_kibz(1) for ik_ibz:", itoa(ik_ibz)))
 !      end associate
 !    end do

 !  end do ! my_it
 !end do ! my_is

 !if (gwr%use_supercell) then

 ! Set FFT mesh in the supercell
 ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
 ! Perhaps the safest approach would be to generate the plan on the fly.
 sc_ngfft = gwr%g_ngfft
 sc_ngfft(1:3) = gwr%ngkpt * gwr%g_ngfft(1:3)
 sc_ngfft(4:6) = sc_ngfft(1:3)
 sc_nfft = product(sc_ngfft(1:3))
 !sc_augsize = product(sc_ngfft(4:6))
 call wrtout(std_out, sjoin(" Building Sigma_c = i GWc with FFTs in the supercell:", ltoa(sc_ngfft(1:3))))

 ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat1, 2))
 ABI_CALLOC(wct_scbox, (sc_nfft * gwr%nspinor * ndat1))

 ! (ngfft, ndat, isign)
 call gt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1 * 2, +1)
 call wt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1, +1)

 ! The g-vectors in the supercell for G and tchi.
 ABI_MALLOC(green_scg, (3, gwr%green_mpw))
 ABI_MALLOC(wc_scg, (3, gwr%tchi_mpw))

 ! Set FFT mesh used to compute u(r) in the unit cell.
 call gwr%kcalc_wfd%change_ngfft(gwr%cryst, gwr%psps, gwr%g_ngfft)

 ABI_CALLOC(cit_me, (2, gwr%ntau, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   ! Load u(r) for GW corrections in the unit cell.
   bmin = minval(gwr%bstart_ks(:, spin))
   bmax = maxval(gwr%bstop_ks(:, spin))
   ABI_MALLOC_OR_DIE(ur_bk, (gwr%g_nfft * gwr%nspinor, bmin:bmax, gwr%nkcalc), ierr)

   do ikcalc=1,gwr%nkcalc
     ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
     ib1 = gwr%bstart_ks(ikcalc, spin)
     ib2 = gwr%bstop_ks(ikcalc, spin)
     call gwr%kcalc_wfd%get_many_ur([(jb, jb=ib1, ib2)], ik_ibz, spin, ur_bk(:, ib1:ib2, ikcalc))
   end do

   do my_it=1,gwr%my_ntau
     call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
     itau = gwr%my_itaus(my_it)

     ! G_k(g,g') --> G_k(g',r) for each k in the BZ treated by me.
     call gwr%get_green_gpr(my_it, my_is, desc_kbz, gt_gpr)

     ! Wc_q(g,g') --> Wc_q(g',r) for each q in the BZ treated by me.
     call gwr%get_wc_gpr(my_it, my_is, desc_qbz, wc_gpr)

     my_nr = gt_gpr(1, 1)%sizeb_local(2)
     ABI_CHECK(my_nr == wc_gpr(1)%sizeb_local(2), "my_nr != wc_gpr(1)%sizeb_local(2)")

     ! Loop over the second index r (MPI-distributed in g_comm).
     do my_ir=1,my_nr

       ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
       ! Note that we need to take the union of (k,g') for k in the BZ.
       gt_scbox = zero
       do my_ikf=1,gwr%my_nkbz
         ik_bz = gwr%my_kbz_inds(my_ikf)
         ik_ibz = gwr%my_kbz2ibz(1, my_ikf)

         ! Compute k + g
         desc_k => desc_kbz(my_ikf)
         do ig=1,desc_k%npw
           green_scg(:,ig) = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ig)
         end do

         do ipm=1,2
           call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat1, green_scg, &
                         gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                         gt_scbox(:, ipm))                             ! inout
         end do
       end do ! my_ikf

       ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)

       ! G(G',r) --> G(R',r)
       call gt_plan_gp2rp%execute_ip_dpc(gt_scbox)

       ! Insert Wc_q(g',r) in G'-space in the supercell FFT box for fixed r.
       ! Note that we need to take the union of (q,g') for q in the BZ.
       ! Also, note ngpt instead of ngkpt.
       wct_scbox = zero
       do my_iqf=1,gwr%my_nqbz
         iq_bz = gwr%my_qbz_inds(my_iqf)
         iq_ibz = gwr%my_qbz2ibz(1, my_iqf)

         ! Compute q + g
         desc_q => desc_qbz(my_iqf)
         do ig=1,desc_q%npw
           wc_scg(:,ig) = nint(gwr%qbz(:, iq_bz) * gwr%ngqpt) + gwr%ngqpt * desc_q%gvec(:,ig)
         end do

         call gsph2box(sc_ngfft, desc_q%npw, gwr%nspinor * ndat1, wc_scg, &
                       wc_gpr(my_iqf)%buffer_cplx(:, my_ir), & ! in
                       wct_scbox)                              ! inout
       end do ! my_iqf

       ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       call xmpi_sum(wct_scbox, gwr%kpt_comm%value, ierr)

       ! Wc(G',r) --> Wc(R',r)
       call wt_plan_gp2rp%execute_ip_dpc(wct_scbox)

       ! Use gt_scbox to store GW (R',r, +/- i tau) for this r (imaginary unit is included outside the loops)
       gt_scbox(:, 1) = gt_scbox(:, 1) * wct_scbox
       gt_scbox(:, 2) = gt_scbox(:, 2) * wct_scbox

       ! Integrate self-energy matrix elements in the R-supercell at fixed r.
       do ikcalc=1,gwr%nkcalc
          ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
          do jb=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
            !ur_bk(:, jb, ikcalc))
            ! Multiply the periodic part by e^{ik.r}
            !call times_eikr(desc%point, gwr%g_ngfft, gwr%g_nfft, ndat1, rgp%buffer_cplx(:, ig2)
            !gt_scbox(:, 1)
            !gt_scbox(:, 2)
            !vals = zero
            !ib = jb - gwr%bstart_ks(ikcalc, spin) + 1
            !cit_me(:, itau, ib, ikcalc, spin) = cit_me(:, itau, ib, ikcalc, spin) + vals(:)
          end do
       end do ! ikcalc

     end do ! my_ir

     ! Free descriptors and PBLAS matrices in kBZ and qBZ.
     call desc_array_free(desc_kbz)
     call slk_array_free(gt_gpr)
     call desc_array_free(desc_qbz)
     call slk_array_free(wc_gpr)

     write(msg,'(2(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "]"
     call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
   end do ! my_it

   ABI_FREE(ur_bk)
 end do ! my_is

 ABI_FREE(gt_scbox)
 ABI_FREE(wct_scbox)
 ABI_FREE(green_scg)
 ABI_FREE(wc_scg)

 ! TODO:
 ! Store matrix elements of Sigma_c(it), separate even and odd part then
 ! use sine/cosine transform to get Sigma_c(i omega).
 ! Finally, perform analytic continuation with Pade' to go to the real-axis
 ! and compute QP corrections and spectral function.

 call xmpi_sum(cit_me, gwr%comm, ierr)
 ABI_CALLOC(ciw_me, (2, gwr%ntau, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol))

 do spin=1,gwr%nsppol
   do ikcalc=1,gwr%nkcalc
      ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
      do jb=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
        !call ft_t2w(cit_me(:,:, ib, ikcalc, spin), ciw_me(:,:, ib, ikcalc, spin))
      end do
   end do
 end do

 ABI_FREE(cit_me)
 ABI_FREE(ciw_me)

 call cwtime_report(" gwr_build_sigmac:", cpu_all, wall_all, gflops_all)

contains

subroutine ft_t2w(t_vals, w_vals)

!Arguments ------------------------------------
 complex(dp),intent(in) :: t_vals(2, gwr%ntau)
 complex(dp),intent(out) :: w_vals(2, gwr%ntau)

!Local variables-------------------------------
 integer :: iw ! itau,
 complex(dp) :: t_odd(gwr%ntau), t_even(gwr%ntau)

! *************************************************************************

 ! f(t) = O(t) + E(t) = (f(t) + f(-t)) / 2  + (f(t) - f(-t)) / 2

 t_odd = (t_vals(1,:) + t_vals(2,:)) / two
 t_even = (t_vals(1,:) - t_vals(2,:)) / two

 do iw=1,gwr%ntau
   !w_vals(1, iw)
   !w_vals(2, iw)
   !matmul(gwr%t2w_cos_wgs(iw,:), t_odd)
   !matmul(gwr%t2w_sin_wgs(iw,:), t_even)
 end do

end subroutine ft_t2w

end subroutine gwr_build_sigmac
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_rpa_energy
!! NAME
!!  gwr_rpa_energy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_rpa_energy(gwr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_is, my_iqi, my_it, itau, spin, iq_ibz
 integer :: il_g1, il_g2, ig1, ig2, iglob2, iglob1, npw_q !, ierr
 complex(dpc) :: vcs_g1, vcs_g2
!arrays
 type(matrix_scalapack) :: chi_tmp, dummy_vec
 real(dp),allocatable :: eig(:)

! *************************************************************************

 call gwr%build_tchi()

 ABI_CHECK(gwr%tchi_space == "iomega", sjoin("tchi_space: ", gwr%tchi_space, " != iomega"))

 ! TODO: Compute RPA total energy with different number of PWs and extrapolate.
 ! See also calc_rpa_functional in m_screening_driver
 !spin_fact = two / (gwr%nsppol * gwr%nspinor)
 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       associate (desc_q => gwr%tchi_desc_qibz(iq_ibz), tchi => gwr%tchi_qibz(iq_ibz, itau, spin))
       if (my_it == 1) then
         ! Allocate workspace. NB: npw_q is the total number of PWs for this q.
         call tchi%copy(chi_tmp)
         npw_q = tchi%sizeb_global(1)
         ABI_MALLOC(eig, (npw_q))
       end if

       do il_g2=1,tchi%sizeb_local(2)
         iglob2 = tchi%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,tchi%sizeb_local(1)
           iglob1 = tchi%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1
           vcs_g1 = desc_q%vc_sqrt(ig1)
           chi_tmp%buffer_cplx(il_g1, il_g2) = tchi%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2
         end do
       end do

       call chi_tmp%pzheev("N", "U", dummy_vec, eig)

       ! Perform integration in imaginary frequency.
       !gwr%iw_wgs(itau)

       if (my_it == gwr%my_ntau) then
         ! Free workspace
         call chi_tmp%free()
         ABI_FREE(eig)
       end if
       end associate

     end do ! my_it
   end do ! my_iqi
 end do ! my_is

 ! Print results to ab_out and nc file
 !call xmpi_sum_master(ec_rpa, master, gwr%comm, ierr)
 !call xmpi_sum_master(ec_gm, master, gwr%comm, ierr)

end subroutine gwr_rpa_energy
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_run_g0w0
!! NAME
!!  gwr_run_g0w0
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_run_g0w0(gwr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr

!Local variables-------------------------------
!scalars
!arrays

! *************************************************************************

 !call gwr%calc_sigx_mels()
 call gwr%build_tchi()
 call gwr%build_wc()
 call gwr%build_sigmac()
 call xmpi_barrier(gwr%comm)

end subroutine gwr_run_g0w0
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_ncwrite_tchi_wc
!! NAME
!!  gwr_ncwrite_tchi_wc
!!
!! FUNCTION
!!  Write tchi or wc to netcdf file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_ncwrite_tchi_wc(gwr, what, filepath)

!Arguments ------------------------------------
 class(gwr_t),target,intent(in) :: gwr
 character(len=*),intent(in) :: what, filepath

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_is, my_iqi, my_it, spin, iq_ibz, itau, npwtot_q, my_ncols, my_gcol_start
 integer :: all_nproc, my_rank, ncid, ncerr, ierr
 real(dp) :: cpu, wall, gflops
!arrays
 integer :: dist_qibz(gwr%nqibz)
 real(dp),ABI_CONTIGUOUS pointer :: fptr(:,:,:)
 type(matrix_scalapack), pointer :: mats(:)

! *************************************************************************

 ! TODO
 ! 1) hscr_new requires ep%
 ! 2) file format assumes Gamma-centered gvectors.

 all_nproc = xmpi_comm_size(gwr%comm); my_rank = xmpi_comm_rank(gwr%comm)
 call cwtime(cpu, wall, gflops, "start")

 if (my_rank == master) then
   call wrtout(std_out, sjoin(" Writing", what, "to:", filepath))
   NCF_CHECK(nctk_open_create(ncid, filepath, xmpi_comm_self))
   NCF_CHECK(gwr%cryst%ncwrite(ncid))

   ! Add dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nsppol", gwr%nsppol), nctkdim_t("mpw", gwr%tchi_mpw), &
     nctkdim_t("ntau", gwr%ntau), &
     nctkdim_t("nqibz", gwr%nqibz), nctkdim_t("nqbz", gwr%nqbz)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define arrays with results.
   ! TODO: Add metadata for mats: spin sum, vc cutoff, t/w mesh, handle nspinor 2
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngqpt", "int", "three"), &
     nctkarr_t("qibz", "dp", "three, nqibz"), &
     nctkarr_t("wtq", "dp", "nqibz"), &
     nctkarr_t("gvec", "int", "three, mpw, nqibz"), &
     nctkarr_t("mats", "dp", "two, mpw, mpw, ntau, nqibz, nsppol") &
   ])
   NCF_CHECK(ncerr)

   ! Write global arrays.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("qibz"), gwr%qibz))
   NCF_CHECK(nf90_put_var(ncid, vid("wtq"), gwr%wtq))
   NCF_CHECK(nf90_close(ncid))
 end if

 call xmpi_barrier(gwr%comm)

 ! The same q-point in the IBZ might be stored on different pools.
 ! To avoid writing the same array multiple times, we use dist_qibz
 ! to select the procs inside gwr%kpt_comm who are gonna write the iq_ibz q-point.
 dist_qibz = huge(1)
 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   dist_qibz(iq_ibz) = gwr%kpt_comm%me
 end do
 call xmpi_min_ip(dist_qibz, gwr%kpt_comm%value, ierr)
 ABI_CHECK(all(dist_qibz /= huge(1)), "Wrong distribution of qibz points in gwr%kpt_comm")

 ! Reopen the file in gwr%comm.
 NCF_CHECK(nctk_open_modify(ncid, filepath, gwr%comm))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     if (dist_qibz(iq_ibz) /= gwr%kpt_comm%me) cycle

     associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))
     npwtot_q = desc_q%npw

     if (spin == 1 .and. gwr%gtau_comm%me == 0) then
       ! Write all G-vectors for this q
       NCF_CHECK(nf90_put_var(ncid, vid("gvec"), desc_q%gvec, start=[1,1,iq_ibz], count=[3,npwtot_q,1]))
     end if

     mats => null()
     if (what == "tchi") mats => gwr%tchi_qibz(iq_ibz, :, spin)
     if (what == "wc")   mats => gwr%wc_qibz(iq_ibz, :, spin)
     ABI_CHECK(associated(mats), sjoin("Invalid value for what:", what))

     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       ! NB: Assuming matrix is distributed in contigous blocks along the column index.
       my_ncols = mats(itau)%sizeb_local(2)
       my_gcol_start = mats(itau)%loc2gcol(1)
       call c_f_pointer(c_loc(mats(itau)%buffer_cplx), fptr, shape=[2, npwtot_q, my_ncols])

       ncerr = nf90_put_var(ncid, vid("mats"), fptr, &
                            start=[1, 1, my_gcol_start, itau, iq_ibz, spin], &
                            count=[2, npwtot_q, my_ncols, 1, 1, 1])
       NCF_CHECK(ncerr)
     end do
     end associate
   end do ! my_iqi
 end do ! my_is

 NCF_CHECK(nf90_close(ncid))

 call cwtime_report(" gwr_ncwrite_tchi_wc:", cpu, wall, gflops)

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
end function vid

end subroutine gwr_ncwrite_tchi_wc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gsph2box
!! NAME
!! gsph2box
!!
!! FUNCTION
!! Insert cg_k array defined on the k-centered g-sphere with npw points inside the FFT box.
!!
!! INPUTS
!! ngfft:
!!   n1,n2,n3=physical dimension of the FFT box
!!   n4,n5,n6=memory dimension of cfft
!! npw=number of G vectors in basis at this k point
!! ndat=number of items to process
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! cg(2,npw*ndat)= contains values for npw G vectors in basis sphere
!!
!! OUTPUT
!! cfft(2,n4,n5,n6*ndat) = array on FFT box filled with cg data
!!      Note that cfft is intent(inout) so that we can add contributions from different k-points.
!!
!! SOURCE

subroutine gsph2box(ngfft, npw, ndat, kg_k, cg, cfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft(6), npw, ndat
!arrays
 integer,intent(in) :: kg_k(3, npw)
 complex(dp),intent(in) :: cg(npw * ndat)
 complex(dp),intent(inout) :: cfft(ngfft(4), ngfft(5), ngfft(6) * ndat)

!Local variables-------------------------------
 integer :: n1, n2, n3, n4, n5, n6, i1, i2, i3, idat, ipw

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Insert cg into cfft
!$OMP PARALLEL DO PRIVATE(i1, i2, i3) IF (ndat > 1)
 do idat=1,ndat
   do ipw=1,npw
     i1 = kg_k(1,ipw); if (i1 < 0) i1 = i1+n1; i1 = i1+1
     i2 = kg_k(2,ipw); if (i2 < 0) i2 = i2+n2; i2 = i2+1
     i3 = kg_k(3,ipw); if (i3 < 0) i3 = i3+n3; i3 = i3+1

     cfft(i1,i2,i3+n6*(idat-1)) = cg(ipw+npw*(idat-1))
   end do
 end do

end subroutine gsph2box
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/box2gsph
!! NAME
!! box2gsph
!!
!! FUNCTION
!! Extract cg_k array defined on the k-centered g-sphere with npw points from the FFT box.
!!
!! INPUTS
!! ngfft:
!!   n1,n2,n3=physical dimension of the FFT box
!!   n4,n5,n6=memory dimension of cfft
!! npw=number of G vectors in basis at this k point
!! ndat=number of items to process
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! cfft(2,n4,n5,n6*ndat) = array on FFT box
!!
!! OUTPUT
!! cg(2,npw*ndat)= contains values for npw G vectors in basis sphere
!!
!! SOURCE

subroutine box2gsph(ngfft, npw, ndat, kg_k, cfft, cg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft(6), npw, ndat
!arrays
 integer,intent(in) :: kg_k(3, npw)
 complex(dp),intent(in) :: cfft(ngfft(4), ngfft(5), ngfft(6) * ndat)
 complex(dp),intent(out) :: cg(npw * ndat)

!Local variables-------------------------------
 integer :: n1, n2, n3, n4, n5, n6, i1, i2, i3, idat, ipw, ipwdat, i3dat

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Extract cg from cfft, ignoring components outside range of cg sphere
!$OMP PARALLEL DO PRIVATE(i1, i2, i3, ipwdat, i3dat) IF (ndat > 1)
 do idat=1,ndat
   do ipw=1,npw
     i1 = kg_k(1,ipw); if (i1 < 0) i1 = i1+n1; i1 = i1+1
     i2 = kg_k(2,ipw); if (i2 < 0) i2 = i2+n2; i2 = i2+1
     i3 = kg_k(3,ipw); if (i3 < 0) i3 = i3+n3; i3 = i3+1
     ipwdat = ipw + (idat-1) * npw
     i3dat = i3 + (idat-1) * n6

     cg(ipwdat) = cfft(i1,i2,i3dat)
   end do
 end do

end subroutine box2gsph
!!***

end module m_gwr
!!***

!!****m* ABINIT/m_gwr
!! NAME
!!  m_gwr
!!
!! FUNCTION
!!  Objects and procedures implementing GW in real-space and imaginary time.
!!
!! NOTES
!!
!! Memory and workload are distributed using a 4D cartesian grid: (g/r, tau, k-points, spin).
!!
!! Inside the g/r communicator, we use PBLAS matrices to store G, tchi and W.
!! using a 1D processor grid with block distribution either along columns or rows.
!! A 2D grid, indeed, would require MPI-FFT or some communication before performing the FFTs.
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
!! At this point, we use ptrans to MPI transpose the (r, g') matrix, and we end up with:
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
!!     the spherical symmetry of vc(r). Besides, when symmetries are used to reconstruct the term for q in the BZ,
!!     one might have to take into account umklapps. Use cache?
!!
!!   - Computation of Sigma_x = Gv must be done in Fourier space using the Lehmann representation so
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
!!  - DONE: New routine for direct diagonalization of the KS Hamiltonian based on PBLAS?
!!
!!  - New routine to compute oscillator matrix elements with NC/PAW and PBLAS matrices.
!!    It can be used to compute tchi head/wings as well as Sigma_x + interface with coupled-cluster codes.
!!
!!  - Decide whether we should use VASP conventions for G and the analytic continuation or the "standard" ones by Godby.
!!    The standard ones are consistent with Hedin's notations and correspond to the ones used in the legacy GW code.
!!    On the other hand, VASP notations make life easier if one has to implement PAW.
!!
!!  - Address nspinor = 2 and PBLAS distribution as MPI proc can have both spinors in memory
!!    In other words, we should store the first/last index in gvec for each spinor
!!
!!  - Optimization for Gamma-only. Memory and c -> r FFTs
!!
!!  - Need to extend FFT API to avoid scaling if isign = -1. Also fft_ug and fft_ur should accept isign
!!    optional argument. Refactoring of all the FFT routines used in the GW code is needed
!!    in order to exploit R2C, C2R (e.g. chi0(q=0) and GPU version.
!!
!!  - Possible incompatibilities between gwpc, slk matrices that are always in dp and GW machinery
!!
!!  - Single precision for scalapack matrices?
!!
!!  - Use round-robin distribution instead of blocked-distribution to improve load balance.
!!
!!  - Memory peaks:
!!
!!      (env3.9) [magianto@uan01 /scratch/project_465000061/magianto/DDIAGO_ZnO]
!!      $~/git_repos/abinit/tests/Scripts/abimem.py peaks abimem_rank0.mocc
!!      [0] <var=gt_scbox, A@m_gwr.F90:3395, addr=0x14aa53673010, size_mb=379.688>
!!      [1] <var=xsum, A@xmpi_sum.finc:2551, addr=0x14aa2fce9010, size_mb=379.688>
!!      [2] <var=gt_scbox, A@m_gwr.F90:4338, addr=0x14aa4f64f010, size_mb=379.688>
!!      [3] <var=allcg_k, A@m_wfd.F90:4631, addr=0x14aa56b57010, size_mb=217.865>
!!      [4] <var=chit_scbox, A@m_gwr.F90:3396, addr=0x14aa4789a010, size_mb=189.844>
!!      [5] <var=wct_scbox, A@m_gwr.F90:4339, addr=0x14aa43876010, size_mb=189.844>
!!      [6] <var=xsum, A@xmpi_sum.finc:2476, addr=0x14aa31bb0010, size_mb=189.844>
!!      [7] <var=cg_k, A@m_wfd.F90:4623, addr=0x14aa64535010, size_mb=108.932>
!!      [8] <var=sc_psi_bk, A@m_gwr.F90:4356, addr=0x14aa3d989010, size_mb=94.922>
!!      [9] <var=sc_ceikr, A@m_gwr.F90:4358, addr=0x10dd28a0, size_mb=47.461>
!!      [10] <var=loc_cwork, A@m_gwr.F90:1894, addr=0x1f3b6fb0, size_mb=11.691>
!!      [11] <var=xsum, A@xmpi_sum_master.finc:1313, addr=0x1ff67d00, size_mb=11.691>
!!      [12] <var=matrix%buffer_cplx, A@m_slk.F90:582, addr=0x8eded70, size_mb=5.845>
!!      [13] <var=Wfd%irottb, A@m_wfd.F90:3721, addr=0x8b4f960, size_mb=3.560>
!!      [14] <var=rhonow, A@m_rhotoxc.F90:613, addr=0x14aa6afc4010, size_mb=2.373>
!!      [15] <var=rhor_file, A@m_ioarr.F90:1045, addr=0x14aa6b0df010, size_mb=1.266>
!!
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
 use m_hide_blas

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use defs_abitypes,   only : mpi_type
 use m_gwdefs,        only : GW_TOLQ0, GW_Q0_DEFAULT
 use m_time,          only : cwtime, cwtime_report, sec2str
 use m_io_tools,      only : iomode_from_fname, get_unit, file_exists, open_file
 use m_numeric_tools, only : blocked_loop, get_diag, isdiagmat, arth, print_arr, pade, dpade, newrap_step, c2r !, linfit
 use m_copy,          only : alloc_copy
 use m_geometry,      only : normv, vdotw
 use m_fstrings,      only : sjoin, itoa, strcat, ktoa, ltoa, ftoa
 use m_sort,          only : sort_dp, sort_rvals
 use m_krank,         only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_crystal,       only : crystal_t
 use m_dtset,         only : dataset_type
 use m_fftcore,       only : get_kg, sphereboundary, getng, print_ngfft, fftcore_set_mixprec
 use m_mpinfo,        only : initmpi_seq, destroy_mpi_enreg
 use m_distribfft,    only : init_distribfft_seq
 use m_fft,           only : fft_ug, fft_ur, fftbox_plan3_t, fourdp
 use m_fft_mesh,      only : calc_ceikr, ctimes_eikr
 use m_kpts,          only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_map_print, kpts_pack_in_stars
 use m_bz_mesh,       only : littlegroup_t
 use m_gsphere,       only : kg_map
 use m_melemts,       only : melements_t
 use m_ioarr,         only : read_rhor
 use m_slk,           only : matrix_scalapack, processor_scalapack, slk_array_free, slk_array_set, slk_array_locmem_mb, &
                             block_dist_1d, slk_pzgemm
 use m_wfk,           only : wfk_read_ebands, wfk_t, wfk_open_read
 use m_wfd,           only : wfd_init, wfd_t, wfdgw_t, wave_t, WFD_STORED
 use m_pawtab,        only : pawtab_type
#ifdef __HAVE_GREENX
 use gx_api,          only : gx_minimax_grid, gx_get_error_message
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
!!  gvectors, tables used for zero padded FFTs and matrix elements
!!  of the coulomb interaction.
!!
!! SOURCE

 type,public :: desc_t

   integer :: istwfk = 1
   ! Storage mode for this k point.

   integer :: npw = -1
   ! Total number of plane-waves for this k/q-point.

   integer,allocatable :: gvec(:,:)
   ! gvec(3, npw)
   ! G-vectors in reduced coordinates.
   ! Note that this array is global i.e. it is not MPI-distributed inside the PBLAS communicator.

   integer,allocatable :: gbound(:,:)
   ! gbound(2*mgfft+8, 2)
   ! sphere boundary info for zero-padded FFT

   complex(gwpc),allocatable :: vc_sqrt(:)
   ! (npw)
   ! Square root of the Coulomb interaction in reciprocal space.
   ! Allocated and computed for tchi/W descriptors.
   ! A cutoff might be applied.

 contains

   procedure :: init => desc_init
   ! Init object

   procedure :: copy => desc_copy
   ! Copy object.

   procedure :: free => desc_free
   ! Free memory.

   procedure :: calc_gnorm_table => desc_calc_gnorm_table
   ! Compute mapping used to loop over G-vectors ordered by norm.

   procedure :: calc_ginv_map => desc_calc_ginv_map
   ! Compute mapping g --> -g

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
!!  This object provides the high-level API to perform the different steps of the GWR algorithm.
!!
!! SOURCE

 type, public :: gwr_t

   integer :: nsppol = 1, nspinor = -1, nspden = -1
   ! Number of independent spin polarizations, number of spinor components and spin density.

   integer :: natom = -1
    ! Number of atoms

   integer :: usepaw = -1
    ! 0 if NC pseudos. 1 if PAW is used (not yet supported).

   integer :: my_nspins = -1
   ! Number of independent spin polarizations treated by this MPI proc

   integer :: nkbz = -1, nkibz = -1
   ! Number of k-points in the BZ/IBZ

   integer :: my_nkibz = -1, my_nkbz = -1
   ! Number of k-points in the IBZ/BZ stored by this MPI proc.

   integer :: uc_batch_size = -1
   ! Max number of unit cell FFT-transforms done in batch mode.

   integer :: sc_batch_size = -1
   ! Max number of supercell-cell FFT-transforms done in batch mode.

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

  integer :: nwr = -1
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Odd number so that the mesh is centered on the KS energy.

  real(dp) :: wr_step = -one
   ! Step of the linear mesh along the real axis (Ha units).

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.

  logical :: timeit = .False.
  ! Internal variable used to activate profiling in low-level routines

  logical :: idle_proc

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

   logical :: use_supercell_for_tchi = .True.
   ! True if we are using the supercell formalism for tchi
   ! False if we are using the mixed-space approach with convolutions in k-space.

   logical :: use_supercell_for_sigma = .True.
   ! True if we are using the supercell formalism for sigma
   ! False if we are using the mixed-space approach with convolutions in k-space.

   integer :: ngkpt(3) = -1, ngqpt(3) = -1
   ! Number of divisions in k/q meshes.

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

   !logical, allocatable :: distrib_kibz(:)

   real(dp),allocatable :: tau_mesh(:), tau_wgs(:)
   ! (ntau)
   ! Imaginary tau mesh and integration weights.

   real(dp),allocatable :: iw_mesh(:), iw_wgs(:)
   ! (ntau)
   ! Imaginary frequency mesh and integration weights

   real(dp),allocatable :: cosft_wt(:,:)
   ! (ntau, ntau)
   ! weights for cosine transform. (i tau --> i omega)

   real(dp),allocatable :: cosft_tw(:,:)
   ! (ntau, ntau)
   ! weights for sine transform (i iomega --> i tau)

   real(dp),allocatable :: sinft_wt(:,:)
   ! (ntau, ntau)
   ! weights for sine transform (i tau --> i omega)

   real(dp) :: ft_max_error(3) = -one
   ! Max error due to inhomogenous FT.

   real(dp) :: cosft_duality_error = -one
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

   integer :: mg0(3) = [2, 2, 2]

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
   ! NB: gwr%comm is not necessarly the same as the input_comm
   ! as we might have removed some cores in input_comm to have nproc divisible by ntau * nsppol.

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

   type(xcomm_t) :: gk_comm
   ! MPI communicator for g/k 2D subgrid.

   type(xcomm_t) :: tks_comm
   ! MPI communicator for tau/kpoint/spin 3D grid

   type(xcomm_t) :: gtk_comm
   ! MPI communicator for g/tau/kpoint 3D grid

   type(dataset_type), pointer :: dtset => null()
   ! Input variables.

   type(datafiles_type), pointer :: dtfil => null()
   ! Names of input/output files and prefixes.

   type(crystal_t), pointer :: cryst => null()
   ! Crystal structure.

   integer :: sc_iteration = 0
   ! Internal counter used to implement self-consistency (not yet implemented)

   type(ebands_t), pointer :: ks_ebands => null()
   ! initial KS energies

   type(ebands_t) :: qp_ebands
   ! QP energies

   type(pseudopotential_type), pointer :: psps => null()
   ! NC Pseudos data

   type(pawtab_type), pointer :: pawtab(:) => null()
   ! PAW data

   type(mpi_type),pointer :: mpi_enreg => null()

   type(processor_scalapack) :: g_slkproc

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
   ! Correlated screened Coulomb interaction summed over collinear spin
   ! Replicated across the spin comm if nsppol == 2.

   character(len=10) :: wc_space = "none"
   ! "none", "itau", "iomega"

   !type(matrix_scalapack),allocatable :: em1_qibz(:,:,:)
   ! Inverse dielectric matrix at omega = 0
   ! (nqibz, nsppol)
   ! Replicated across the tau comm and the spin comm if nsppol == 2.

   type(matrix_scalapack),allocatable :: sigc_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)

   character(len=10) :: sigc_space = "none"
   ! "none", "itau", "iomega"

   type(matrix_scalapack),allocatable :: ugb(:,:)
   ! (nkibz, nsppol)
   ! Fourier components of the KS wavefunctions.
   ! Bands are distributed inside the gtau_comm in round-robin fashion.

   type(processor_scalapack) :: gtau_slkproc

   integer :: ugb_nband = -1

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

   integer,allocatable :: kbz2ibz(:,:)
    ! (6, nkbz)
    ! Mapping kBZ to IBZ (symrec conventions)

   real(dp), contiguous, pointer :: wtk(:) => null()
    ! (nkibz)
    ! Weights of the k-points in the IBZ (normalized to one).

   real(dp),allocatable :: qbz(:,:)
    ! (3, nqbz)
    ! Reduced coordinates of the q-points in the full BZ.

   integer,allocatable :: qbz2ibz(:,:)
   ! (6, nqbz)
   ! Mapping qBZ to IBZ (symrec conventions)

   real(dp),allocatable :: qibz(:,:)
   ! (3, nqibz)
   ! Reduced coordinates of the q-points in the IBZ (full simmetry of the system).

   real(dp),allocatable :: wtq(:)
   ! (nqibz)
   ! Weights of the q-points in the IBZ (normalized to one).

   complex(dp),allocatable :: chi0_head_myw(:,:,:)
   ! (3, 3, my_ntau)
   ! Head of the irred polarizability in i omega space.
   ! Note that spins have been summed over.

   complex(dp),allocatable :: chi0_uwing_myw(:,:,:), chi0_lwing_myw(:,:,:)
   ! (3, npw_chi_gamma, my_ntau)
   ! Upper wings of the irred polarizability in i omega space.
   ! Note that spins have been summed over.

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

   procedure :: rotate_gpm => gwr_rotate_gpm
   ! Reconstruct the Green's functions in the kBZ from the IBZ.

   procedure :: get_myk_green_gpr => gwr_get_myk_green_gpr
    ! G_k(g,g') --> G_k(g',r) for each k in the BZ treated by this MPI proc for given spin and tau.

   procedure :: get_gk_rpr_pm => gwr_get_gk_rpr_pm
   ! Compute G_k(r',r) with (r, r') in the unit cell and k in the full BZ.

   procedure :: rotate_wc => gwr_rotate_wc
   ! Reconstruct Wc(q) in the BZ from the IBZ.

   procedure :: get_myq_wc_gpr => gwr_get_myq_wc_gpr
   ! W_q(g,g') --> W_q(g',r) for each q in the BZ treated by this MPI procs for given spin and tau.

   procedure :: get_wc_rpr_qbz => gwr_get_wc_rpr_qbz

   procedure :: cos_transform  => gwr_cos_transform
   ! Inhomogeneous cosine transform.

   procedure :: alloc_free_mats => gwr_alloc_free_mats
   ! Allocate/Deallocate matrices for G/tchi/Sigma

   procedure :: free => gwr_free
   ! Free memory.

   procedure :: print => gwr_print
   ! Print info on the object.

   procedure :: print_mem => gwr_print_mem

   procedure :: print_trace => gwr_print_trace
    ! Print trace of matrices for testing purposes.

   procedure :: load_kcalc_wfd => gwr_load_kcalc_wfd
   ! Load the KS states for Sigma_nk from the WFK file

   procedure :: read_ugb_from_wfk => gwr_read_ugb_from_wfk
   ! Read wavefunctions from WFK file.

   procedure :: build_green => gwr_build_green
   ! Build the Green's function in imaginary time for k-points in the IBZ from the WFK file

   procedure :: build_tchi => gwr_build_tchi
   ! Build the irreducible polarizability

   procedure :: distrib_gt_kibz => gwr_distrib_gt_kibz
   ! Redistribute/deallocate G_k

   procedure :: distrib_mats_qibz => gwr_distrib_mats_qibz
   ! Redistribute/deallocate tchi_q or Wc_q

   procedure :: build_wc => gwr_build_wc
   ! Build the correlated part of the screened interaction.

   procedure :: build_sigmac => gwr_build_sigmac
   ! Build the correlated part of the self-energy GWc
   ! and compute matrix elements in the KS representation.

   procedure :: rpa_energy => gwr_rpa_energy
   ! Compute RPA energy.

   procedure :: build_chi0_head_and_wings => gwr_build_chi0_head_and_wings

   procedure :: build_sigxme => gwr_build_sigxme

   procedure :: run_g0w0 => gwr_run_g0w0
   ! Compute QP corrections with G0W0.

   procedure :: ncwrite_tchi_wc => gwr_ncwrite_tchi_wc
   ! Write tchi or wc to netcdf file

 end type gwr_t
!!***

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha

 integer,private,parameter :: PRINT_MODR = 20
 integer,private,parameter :: istwfk1 = 1, ndat1 = 1
 integer,private,parameter :: me_fft0 = 0, paral_fft0 = 0, nproc_fft1 = 1

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

subroutine gwr_init(gwr, dtset, dtfil, cryst, psps, pawtab, ks_ebands, mpi_enreg, input_comm)

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
 integer,intent(in) :: input_comm

!Local variables-------------------------------
!scalars
 integer,parameter :: qptopt1 = 1, qtimrev1 = 1, master = 0, ndims = 4
 integer :: my_it, my_ikf, my_iqf, ii, ebands_timrev, my_iki, my_iqi, itau, spin
 integer :: my_nshiftq, iq_bz, iq_ibz, npw_, ncid, ig, ig_start
 integer :: comm_cart, me_cart, ierr, all_nproc, nps, my_rank, qprange_, gap_err, ncerr, omp_nt
 integer :: jj, cnt, ikcalc, ndeg, mband, bstop, nbsum, it, iw
 integer :: ik_ibz, ik_bz, isym_k, trev_k, g0_k(3)
 real(dp) :: cpu, wall, gflops, te_min, te_max, wmax
 logical :: isirr_k, changed, q_is_gamma, reorder
 character(len=5000) :: msg
 type(krank_t) :: qrank, krank_ibz
 type(gaps_t) :: ks_gaps
!arrays
 integer :: qptrlatt(3,3), dims(ndims), indkk_k(6,1)
 integer,allocatable :: gvec_(:,:),degblock(:,:), degblock_all(:,:,:,:), ndeg_all(:,:)
 real(dp) :: my_shiftq(3,1), kk_ibz(3), kk_bz(3), qq_bz(3), qq_ibz(3), kk(3)
 real(dp),allocatable :: wtk(:), kibz(:,:)
 logical :: periods(ndims), keepdim(ndims)

! *************************************************************************

 all_nproc = xmpi_comm_size(input_comm); my_rank = xmpi_comm_rank(input_comm)

 call cwtime(cpu, wall, gflops, "start")

 gwr%dtset => dtset
 gwr%dtfil => dtfil
 gwr%cryst => cryst
 gwr%psps => psps
 gwr%pawtab => pawtab
 gwr%ks_ebands => ks_ebands
 gwr%kibz => ks_ebands%kptns
 gwr%wtk => ks_ebands%wtk
 gwr%mpi_enreg => mpi_enreg
 gwr%timeit = dtset%prtvol > 20

 ! Initialize qp_ebands with KS values.
 call ebands_copy(ks_ebands, gwr%qp_ebands)

 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%nspden = dtset%nspden
 gwr%natom = dtset%natom
 gwr%usepaw = dtset%usepaw

 gwr%use_supercell_for_tchi = .True.
 if (gwr%dtset%useria == 1) gwr%use_supercell_for_tchi = .False.
 gwr%use_supercell_for_sigma = .True.
 if (gwr%dtset%userib == 1) gwr%use_supercell_for_sigma = .False.

 mband = ks_ebands%mband
 nbsum = dtset%nband(1)
 ABI_CHECK_IRANGE(nbsum, 1, mband, "Invalid nbsum")

 ! Define Frequency mesh for sigma(w_real) and spectral functions.
 ! Note that in GWR computing quantities on the real-axis is really cheap
 ! so we can use very dense meshes without affecting performance.
 wmax = dtset%freqspmax; if (abs(wmax) < tol6) wmax = 100 * eV_Ha
 gwr%nwr = dtset%nfreqsp
 if (gwr%nwr ==  0) gwr%nwr = nint(wmax / (0.05 * eV_Ha))
 if (mod(gwr%nwr, 2) == 0) gwr%nwr = gwr%nwr + 1
 gwr%wr_step = wmax / (gwr%nwr - 1)

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
 ABI_CHECK_IEQ(gwr%nkibz, ks_ebands%nkpt, "nkibz != ks_ebands%nkpt")
 ABI_CHECK(all(abs(ks_ebands%kptns - kibz) < tol12), "ks_ebands%kibz != kibz")

 if (.not. (isdiagmat(ks_ebands%kptrlatt) .and. ks_ebands%nshiftk == 1)) then
   ABI_ERROR("GWR requires ngkpt with one shift!")
 end if
 gwr%ngkpt = get_diag(ks_ebands%kptrlatt)

 ! Note symrec convention
 ABI_MALLOC(gwr%kbz2ibz, (6, gwr%nkbz))
 ebands_timrev = kpts_timrev_from_kptopt(ks_ebands%kptopt)

 krank_ibz = krank_from_kptrlatt(gwr%nkibz, kibz, ks_ebands%kptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", ebands_timrev, cryst, krank_ibz, gwr%nkbz, gwr%kbz, gwr%kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 ! Order kbz by stars and rearrange entries in kbz2ibz table.
 call kpts_pack_in_stars(gwr%nkbz, gwr%kbz, gwr%kbz2ibz)

 if (my_rank == master) then
   call kpts_map_print([std_out, ab_out], " Mapping kBZ --> kIBZ", "symrec", &
                       gwr%kbz, kibz, gwr%kbz2ibz, gwr%dtset%prtvol)
 end if

 !call get_ibz2bz(gwr%nkibz, gwr%nkbz, gwr%kbz2ibz, kibz2bz, ierr)
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
 ABI_MALLOC(gwr%qbz2ibz, (6, gwr%nqbz))

 qrank = krank_from_kptrlatt(gwr%nqibz, gwr%qibz, qptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", qtimrev1, cryst, qrank, gwr%nqbz, gwr%qbz, gwr%qbz2ibz) /= 0) then
   ABI_ERROR("Cannot map qBZ to IBZ!")
 end if
 call qrank%free()

 ! Order qbz by stars and rearrange entries in qbz2ibz table.
 call kpts_pack_in_stars(gwr%nqbz, gwr%qbz, gwr%qbz2ibz)
 if (my_rank == master) then
   call kpts_map_print([std_out, ab_out], " Mapping qBZ --> qIBZ", "symrec", &
                       gwr%qbz, gwr%qibz, gwr%qbz2ibz, gwr%dtset%prtvol)
 end if

 ! ==========================
 ! Setup k-points in Sigma_nk
 ! ==========================
 ks_gaps = ebands_get_gaps(ks_ebands, gap_err)
 if (my_rank == master) then
   !call ebands_print(ks_ebands, header="KS band structure", unit=std_out, prtvol=gwr%dtset%prtvol)
   msg = "Kohn-Sham gaps and band edges from IBZ mesh"
   call ks_gaps%print(unit=std_out, header=msg)
   call ks_gaps%print(unit=ab_out, header=msg)
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
     call sigtk_kcalc_from_erange(dtset, cryst, ks_ebands, ks_gaps, &
                                  gwr%nkcalc, gwr%kcalc, gwr%bstart_ks, gwr%nbcalc_ks, input_comm)

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
       call sigtk_kcalc_from_gaps(dtset, ks_ebands, ks_gaps, gwr%nkcalc, gwr%kcalc, gwr%bstart_ks, gwr%nbcalc_ks)
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

 ! FIXME: This trigger SIGSEGV on lumi, don't know why!
#if 0
 ! Build degtab tables.
 if (abs(gwr%dtset%symsigma) == 1) then
   call xmpi_sum(ndeg_all, input_comm, ierr)
   call xmpi_sum(degblock_all, input_comm, ierr)
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
#endif

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
 te_min = minval(ks_gaps%cb_min - ks_gaps%vb_max)
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
 gwr%iw_mesh = arth(zero, te_max, gwr%ntau)
 gwr%iw_wgs = te_max / (gwr%ntau - 1)
 gwr%tau_mesh = arth(zero, two_pi * gwr%iw_mesh(2), gwr%ntau)
 gwr%tau_wgs = one / gwr%ntau

 ABI_MALLOC(gwr%cosft_tw, (gwr%ntau, gwr%ntau))
 ABI_MALLOC(gwr%cosft_wt, (gwr%ntau, gwr%ntau))
 ABI_MALLOC(gwr%sinft_wt, (gwr%ntau, gwr%ntau))

 do it=1, gwr%ntau
   do iw=1, gwr%ntau
      gwr%cosft_tw(it,iw) = cos(gwr%iw_mesh(iw) * gwr%tau_mesh(it)) * two * gwr%iw_mesh(gwr%ntau) / (gwr%ntau - 1)
      gwr%cosft_wt(iw,it) = cos(gwr%iw_mesh(iw) * gwr%tau_mesh(it)) * two * gwr%tau_mesh(gwr%ntau) / (gwr%ntau - 1)
   end do
 end do
 gwr%sinft_wt = one

#else
 call gx_minimax_grid(gwr%ntau, te_min, te_max, &  ! in
                      gwr%tau_mesh, gwr%tau_wgs, gwr%iw_mesh, gwr%iw_wgs, & ! out args allocated by the routine.
                      gwr%cosft_wt, gwr%cosft_tw, gwr%sinft_wt, &
                      gwr%ft_max_error, gwr%cosft_duality_error, ierr)
 if (ierr /= 0) then
   call gx_get_error_message(msg)
   ABI_ERROR(msg)
 end if

 call wrtout(std_out, sjoin("Max_{ij} |CT CT^{-1} - I|", ftoa(gwr%cosft_duality_error)))
 call wrtout(std_out, sjoin("ft_max_error", ltoa(gwr%ft_max_error)))
#endif

 call ks_gaps%free()

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! NB: Here we define gwr%comm and gwr%idle_proc
 ! Do not use input_comm after this section as idle processors return immediately.

 if (any(dtset%gwr_np_gtks /= 0)) then
   ! Use MPI parameters from input file.
   gwr%g_comm%nproc    = dtset%gwr_np_gtks(1)
   gwr%tau_comm%nproc  = dtset%gwr_np_gtks(2)
   gwr%kpt_comm%nproc  = dtset%gwr_np_gtks(3)
   gwr%spin_comm%nproc = dtset%gwr_np_gtks(4)

   !call xmpi_comm_multiple_of(gwr%tau_comm%nproc * gwr%spin_comm%nproc, input_comm, gwr%idle_proc, gwr%comm)
   !if (gwr%idle_proc) return
   gwr%comm = input_comm
#ifdef HAVE_MPI
   call MPI_Comm_dup(input_comm, gwr%comm, ierr)
#endif
   all_nproc = xmpi_comm_size(gwr%comm)
 else
   ! Automatic grid generation.
   !
   !   Priorities        |  MPI Scalability                | Memory
   ! ==================================================================================================
   !   spin (if any)     |  excellent                      | scales
   !   tau               |  excellent                      | scales
   !   kbz               |  to be optimized!               | scales (depends on the BZ -> IBZ mapping)
   !   g/r (PBLAS)       |  network intensive              ! scales
   !

   call xmpi_comm_multiple_of(gwr%ntau * gwr%dtset%nsppol, input_comm, gwr%idle_proc, gwr%comm)
   if (gwr%idle_proc) return
   all_nproc = xmpi_comm_size(gwr%comm)

   gwr%spin_comm%nproc = 1
   !if (gwr%nsppol == 2 .and. all_nproc > 1) then
   !  ABI_CHECK(mod(all_nproc, 2) == 0, "when nsppol == 2, nprocs should be even!")
   !  gwr%spin_comm%nproc = 2
   !end if

   nps = all_nproc / gwr%spin_comm%nproc
   do ii=nps,1,-1
     if (mod(gwr%ntau, ii) == 0 .and. mod(nps, ii) == 0) exit
   end do

   if (ii == 1 .and. nps > 1) then
     if (gwr%nkbz > 1) then
       ! Give priority to tau/kbz
       call xmpi_distrib_2d(nps, "12", gwr%ntau, gwr%nkbz, gwr%tau_comm%nproc, gwr%kpt_comm%nproc, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(nps), " with priority: tau/kbz"))
     else
       ! Give priority to tau/g
       call xmpi_distrib_2d(nps, "12", gwr%ntau, gwr%green_mpw, gwr%tau_comm%nproc, gwr%g_comm%nproc, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(nps), " with priority: tau/g"))
     end if
   else
     ! ii divides ntau and nps.
     gwr%tau_comm%nproc = ii
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
         ! In this case, ii divides nkbz and nps.
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
 call gwr%g_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for tau
 keepdim = .False.; keepdim(2) = .True.
 call gwr%tau_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for k-points
 keepdim = .False.; keepdim(3) = .True.
 call gwr%kpt_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for spin
 keepdim = .False.; keepdim(4) = .True.
 call gwr%spin_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for the (g, tau) 2D grid.
 keepdim = .False.; keepdim(1) = .True.; keepdim(2) = .True.
 call gwr%gtau_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for the (g, tau) 2D grid.
 keepdim = .False.; keepdim(1) = .True.; keepdim(3) = .True.
 call gwr%gk_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for the (g, tau, k) 3D subgrid.
 keepdim = .True.; keepdim(4) = .False.
 call gwr%gtk_comm%from_cart_sub(comm_cart, keepdim)

 ! Create communicator for the (tau, k, spin) 3D subgrid.
 keepdim = .True.; keepdim(1) = .False.
 call gwr%tks_comm%from_cart_sub(comm_cart, keepdim)

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
 ABI_CHECK(all(gwr%tau_master > -1), "tau_master!")

 call xmpi_split_block(gwr%nsppol, gwr%spin_comm%value, gwr%my_nspins, gwr%my_spins)
 ABI_CHECK(gwr%my_nspins > 0, "my_nspins == 0, decrease number of procs for spin level")

 ! Distribute k-points in the full BZ, build redirection tables.
 ! Finally, find the number of IBZ k-points stored by this MPI rank.

 call xmpi_split_block(gwr%nkbz, gwr%kpt_comm%value, gwr%my_nkbz, gwr%my_kbz_inds)
 ABI_CHECK(gwr%my_nkbz > 0, "my_nkbz == 0, decrease number of procs for k-point level")

 ABI_ICALLOC(gwr%np_kibz, (gwr%nkibz))
 do my_ikf=1,gwr%my_nkbz
   ik_bz = gwr%my_kbz_inds(my_ikf)
   ik_ibz = gwr%kbz2ibz(1, ik_bz)
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

 ! Distribute q-points in the full BZ, transfer symmetry tables.
 ! Finally find the number of IBZ q-points that should be stored in memory.
 call xmpi_split_block(gwr%nqbz, gwr%kpt_comm%value, gwr%my_nqbz, gwr%my_qbz_inds)

 ABI_ICALLOC(gwr%np_qibz, (gwr%nqibz))

 do my_iqf=1,gwr%my_nqbz
   iq_bz = gwr%my_qbz_inds(my_iqf)
   iq_ibz = gwr%qbz2ibz(1, iq_bz)
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

 ! =========================================
 ! Find FFT mesh and max number of g-vectors
 ! =========================================
 ! Note the usage of gwr_boxcutmin and the loops over the full BZ.
 ! All the procs execute this part.

 gwr%g_ngfft = gwr%dtset%ngfft ! Allow user to specify fftalg
 gwr%g_ngfft(1:6) = 0

 gwr%green_mpw = -1
 do ik_bz=1,gwr%nkbz
   kk_bz = gwr%kbz(:, ik_bz)
   call get_kg(kk_bz, istwfk1, dtset%ecut, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%gwr_boxcutmin, dtset%chksymtnons, dtset%ecut, cryst%gmet, &
              kk_bz, me_fft0, gwr%g_mgfft, gwr%g_nfft, gwr%g_ngfft, nproc_fft1, cryst%nsym, paral_fft0, &
              cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%green_mpw = max(gwr%green_mpw, npw_)
 end do

 gwr%tchi_mpw = -1
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

 ! Define batch sizes for FFT transforms, use multiples of OpenMP threads.
 omp_nt = xomp_get_num_threads(open_parallel=.True.)
 gwr%uc_batch_size = max(1, gwr%dtset%userid * omp_nt)
 gwr%sc_batch_size = max(1, gwr%dtset%userie * omp_nt)

 !gwr%uc_batch_size = 4; gwr%sc_batch_size = 4

 if (my_rank == master) then
   call print_ngfft(gwr%g_ngfft, header="FFT mesh for Green's function", unit=std_out)
   call wrtout(std_out, sjoin(" FFT uc_batch_size:", itoa(gwr%uc_batch_size)))
   call wrtout(std_out, sjoin(" FFT sc_batch_size:", itoa(gwr%sc_batch_size)))
 end if

 ! Now we know the value of g_ngfft. Setup tables for zero-padded FFTs.
 ! Build descriptors for Green's functions and tchi and setup tables for zero-padded FFTs.
 ABI_MALLOC(gwr%green_desc_kibz, (gwr%nkibz))

 do my_iki=1,gwr%my_nkibz
   ik_ibz = gwr%my_kibz_inds(my_iki)
   kk_ibz = gwr%kibz(:, ik_ibz)
   call gwr%green_desc_kibz(ik_ibz)%init(kk_ibz, istwfk1, dtset%ecut, gwr)
 end do

 ABI_MALLOC(gwr%tchi_desc_qibz, (gwr%nqibz))

 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   q_is_gamma = (normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0)
   ! Note ecuteps instead of ecut.
   call gwr%tchi_desc_qibz(iq_ibz)%init(qq_ibz, istwfk1, dtset%ecuteps, gwr)
   associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))
   ! Compute sqrt(v(q,G))
   !TODO: Implement cutoff in vc
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

 ! Init 1D PBLAS grid to block-distribute matrices along columns.
 call gwr%g_slkproc%init(gwr%g_comm%value, grid_dims=[1, gwr%g_comm%nproc])
 call gwr%gtau_slkproc%init(gwr%gtau_comm%value, grid_dims=[1, gwr%gtau_comm%nproc])

 ! ==================================
 ! Allocate arrays of PBLAS matrices
 ! ==================================

 ABI_MALLOC(gwr%gt_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))
 ABI_MALLOC(gwr%tchi_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))
 ABI_MALLOC(gwr%sigc_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))

 ! Create netcdf file to store results.
 gwr%gwrnc_path = strcat(dtfil%filnam_ds(4), "_GWR.nc")

 if (my_rank == master) then
   call gwr%print(unit=ab_out)
   call gwr%print(unit=std_out)

   NCF_CHECK(nctk_open_create(ncid, gwr%gwrnc_path, xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ks_ebands, ncid))

   ! Add dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nsppol", gwr%nsppol), nctkdim_t("ntau", gwr%ntau), nctkdim_t("nwr", gwr%nwr), &
     nctkdim_t("nkcalc", gwr%nkcalc), nctkdim_t("max_nbcalc", gwr%max_nbcalc) &
     !nctkdim_t("nqibz", gwr%nqibz), nctkdim_t("nqbz", gwr%nqbz)
     ], defmode=.True.)
   NCF_CHECK(ncerr)

   !ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", "symdynmat", &
   !  "imag_only", "symv1scf", "dvdb_add_lr", "mrta", "ibte_prep", "eph_prtscratew"])
   !NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "wr_step" & ! "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", &
   ])
   NCF_CHECK(ncerr)

   ! Define arrays with results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("tau_mesh", "dp", "ntau"), &
     nctkarr_t("iw_mesh", "dp", "ntau"), &
     !nctkarr_t("ngqpt", "int", "three"), &
     !nctkarr_t("ddb_ngqpt", "int", "three"), &
     !nctkarr_t("sigma_ngkpt", "int", "three"), &
     !nctkarr_t("sigma_erange", "dp", "two"), &
     !nctkarr_t("bstart_ks", "int", "nkcalc, nsppol"), &
     !!nctkarr_t("bstop_ks", "int", "nkcalc, nsppol"), &
     !nctkarr_t("nbcalc_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("kcalc", "dp", "three, nkcalc") &
     !nctkarr_t("kcalc2ibz", "int", "nkcalc, six"), &
     !nctkarr_t("qp_done", "int", "nkcalc, nsppol"), &
     !nctkarr_t("vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("dvals_de0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("qpoms_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("qp_enes", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("ze0_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("ks_gaps", "dp", "nkcalc, nsppol"), &
     !nctkarr_t("qpoms_gaps", "dp", "ntemp, nkcalc, nsppol"), &
     !nctkarr_t("qp_gaps", "dp", "ntemp, nkcalc, nsppol"), &
     !nctkarr_t("vcar_calc", "dp", "three, max_nbcalc, nkcalc, nsppol") &
   ])
   NCF_CHECK(ncerr)

   ! ======================================================
   ! Write data that do not depend on the (kpt, spin) loop.
   ! ======================================================
   NCF_CHECK(nctk_set_datamode(ncid))
   !ii = 0; if (gwr%imag_only) ii = 1
   !ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", &
   !  "eph_transport", "imag_only", "symv1scf", "dvdb_add_lr", "mrta", "ibte_prep", "eph_prtscratew"], &
   !  [dtset%eph_task, gwr%symsigma, gwr%nbsum, gwr%bsum_start, gwr%bsum_stop, &
   !  dtset%symdynmat, dtset%ph_intmeth, dtset%eph_intmeth, gwr%qint_method, &
   !  dtset%eph_transport, ii, dtset%symv1scf, dtset%dvdb_add_lr, gwr%mrta, dtset%ibte_prep, dtset%eph_prtscratew])
   !NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "wr_step"], &
     [gwr%wr_step])
   NCF_CHECK(ncerr)

   NCF_CHECK(nf90_put_var(ncid, vid("tau_mesh"), gwr%tau_mesh))
   NCF_CHECK(nf90_put_var(ncid, vid("iw_mesh"), gwr%iw_mesh))
   NCF_CHECK(nf90_put_var(ncid, vid("kcalc"), gwr%kcalc))
   NCF_CHECK(nf90_close(ncid))
 end if ! master

 call cwtime_report(" gwr_init:", cpu, wall, gflops)

contains
integer function vid(vname)
  character(len=*),intent(in) :: vname
  vid = nctk_idname(ncid, vname)
end function vid

end subroutine gwr_init
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_alloc_free_mats
!! NAME
!! gwr_alloc_free_mats
!!
!! FUNCTION
!!
!! SOURCE

subroutine gwr_alloc_free_mats(gwr, mask_ibz, what, action)

!Arguments ------------------------------------
 class(gwr_t), target, intent(inout) :: gwr
 integer,intent(in) :: mask_ibz(:)
 character(len=*),intent(in) :: what, action

!Local variables-------------------------------
 integer :: my_is, my_it, ipm, npwsp, col_bsize, itau, spin, ik_ibz, iq_ibz
 type(matrix_scalapack), pointer :: mat
 character(len=500) :: msg

! *************************************************************************

 ABI_CHECK(action == "malloc" .or. action == "free", sjoin("Invalid action:", action))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_it=1,gwr%my_ntau
     itau = gwr%my_itaus(my_it)
     ! All the PBLAS matrices are MPI distributed over g' in blocks

     select case (what)
     case ("green")
       ! ==================================
       ! Allocate G_k(g,g') for k in my IBZ
       ! ==================================
       ABI_CHECK_IEQ(size(mask_ibz), gwr%nkibz, "wrong mask size")

       do ik_ibz=1,gwr%nkibz
         if (mask_ibz(ik_ibz) == 0) cycle
         npwsp = gwr%green_desc_kibz(ik_ibz)%npw * gwr%nspinor
         ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
         associate (gt => gwr%gt_kibz(:, ik_ibz, itau, spin))
         do ipm=1,2
           if (action == "malloc") call gt(ipm)%init(npwsp, npwsp, gwr%g_slkproc, istwfk1, size_blocs=[npwsp, col_bsize])
           if (action == "free") call gt(ipm)%free()
         end do
         end associate
       end do

     case ("tchi", "wc")
       ! =====================================
       ! Allocate tchi_q(g,g') for q in my IBZ
       ! =====================================
       ABI_CHECK_IEQ(size(mask_ibz), gwr%nqibz, "wrong mask size")

       do iq_ibz=1,gwr%nqibz
         if (mask_ibz(iq_ibz) == 0) cycle
         npwsp = gwr%tchi_desc_qibz(iq_ibz)%npw * gwr%nspinor
         ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
         if (what == "tchi") mat => gwr%tchi_qibz(iq_ibz, itau, spin)
         if (what == "wc") mat => gwr%wc_qibz(iq_ibz, itau, spin)
         if (action == "malloc") call mat%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
         if (action == "free") call mat%free()
       end do

    case ("sigma")
       ! ===========================================
       ! Allocate PBLAS arrays for sigmac_kibz(g,g')
       ! ===========================================
       ABI_CHECK_IEQ(size(mask_ibz), gwr%nkibz, "wrong mask size")
       do ik_ibz=1,gwr%nkibz
         if (mask_ibz(ik_ibz) == 0) cycle
         npwsp = gwr%tchi_desc_qibz(iq_ibz)%npw * gwr%nspinor
         ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
         associate (sigc => gwr%sigc_kibz(:, ik_ibz, itau, spin))
         do ipm=1,2
           if (action == "malloc") call sigc(ipm)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
           if (action == "free") call sigc(ipm)%free()
         end do
         end associate
       end do

     case default
       ABI_ERROR(sjoin("Invalid what:", what))
     end select

   end do ! my_it
 end do ! my_is

 call wrtout(std_out, "")
 call gwr%print_mem(unit=std_out)

end subroutine gwr_alloc_free_mats
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

! *************************************************************************

 ABI_SFREE(gwr%kbz)
 ABI_SFREE(gwr%kbz2ibz)
 ABI_SFREE(gwr%qbz2ibz)
 ABI_SFREE(gwr%my_kbz_inds)
 ABI_SFREE(gwr%my_kibz_inds)
 ABI_SFREE(gwr%my_qbz_inds)
 ABI_SFREE(gwr%my_qibz_inds)
 ABI_SFREE(gwr%qbz)
 ABI_SFREE(gwr%qibz)
 ABI_SFREE(gwr%wtq)
 ABI_SFREE(gwr%chi0_head_myw)
 ABI_SFREE(gwr%chi0_uwing_myw)
 ABI_SFREE(gwr%chi0_lwing_myw)
 ABI_SFREE(gwr%qbz2ibz)
 ABI_SFREE(gwr%my_spins)
 ABI_SFREE(gwr%my_itaus)
 ABI_SFREE(gwr%tau_master)
 ABI_SFREE(gwr%np_kibz)
 ABI_SFREE(gwr%np_qibz)
#ifndef __HAVE_GREENX
 ABI_SFREE(gwr%tau_mesh)
 ABI_SFREE(gwr%tau_wgs)
 ABI_SFREE(gwr%iw_mesh)
 ABI_SFREE(gwr%iw_wgs)
 ABI_SFREE(gwr%cosft_tw)
 ABI_SFREE(gwr%cosft_wt)
 ABI_SFREE(gwr%sinft_wt)
#else
 ABI_SFREE_NOCOUNT(gwr%tau_mesh)
 ABI_SFREE_NOCOUNT(gwr%tau_wgs)
 ABI_SFREE_NOCOUNT(gwr%iw_mesh)
 ABI_SFREE_NOCOUNT(gwr%iw_wgs)
 ABI_SFREE_NOCOUNT(gwr%cosft_tw)
 ABI_SFREE_NOCOUNT(gwr%cosft_wt)
 ABI_SFREE_NOCOUNT(gwr%sinft_wt)
#endif
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

 if (allocated(gwr%ugb)) then
   call slk_array_free(gwr%ugb)
   ABI_FREE(gwr%ugb)
 end if
 call gwr%gtau_slkproc%free()

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
 call gwr%gk_comm%free()
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

!!****f* m_gwr/gwr_read_ugb_from_wfk
!! NAME
!!  gwr_read_ugb_from_wfk
!!
!! FUNCTION

!!  Read wavefunctions from WFK file.
!!
!! SOURCE

subroutine gwr_read_ugb_from_wfk(gwr, wfk_path)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: wfk_path

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0 = 0
 integer :: mband, min_nband, nkibz, nsppol, my_is, my_iki, spin, ik_ibz
 integer :: my_ib, my_nb, band, npw_k, mpw, istwf_k !, itau !, il_b
 integer :: nbsum, npwsp, col_bsize, nband_k, my_bstart, my_bstop, my_nband, ierr
 real(dp) :: f_nk, eig_nk, ef, cpu, wall, gflops, cpu_green, wall_green, gflops_green
 character(len=5000) :: msg
 type(ebands_t) :: wfk_ebands
 type(hdr_type) :: wfk_hdr
 type(wfk_t) :: wfk
 type(dataset_type),pointer :: dtset
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kk_ibz(3)
 !real(dp),allocatable :: cgwork(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: cg_k(:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 dtset => gwr%dtset

 wfk_ebands = wfk_read_ebands(wfk_path, gwr%comm, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ! TODO: More consistency checks e.g. nkibz,...

 !cryst = wfk_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 nkibz = wfk_ebands%nkpt; nsppol = wfk_ebands%nsppol; mband = wfk_ebands%mband
 min_nband = minval(wfk_ebands%nband)

 nbsum = dtset%nband(1)
 if (nbsum > min_nband) then
   ABI_WARNING(sjoin("WFK file contains", itoa(min_nband), "states while you're asking for:", itoa(nbsum)))
   nbsum = min_nband
 end if
 call wrtout(std_out, sjoin(" Computing Green's function with nbsum:", itoa(nbsum)))
 call ebands_free(wfk_ebands)

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
 !
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

 ! Select occupied or empty G.
 ! if (eig_nk < -tol6) then
 !   !ipm = 1
 !   !gt_cfact = j_dpc * exp(gwr%tau_mesh(itau) * eig_nk)
 !   ! Vasp convention
 !   ipm = 2
 !   gt_cfact = exp(gwr%tau_mesh(itau) * eig_nk)
 ! else if (eig_nk > tol6) then
 !   !ipm = 2
 !   !gt_cfact = -j_dpc * exp(-gwr%tau_mesh(itau) * eig_nk)
 !   ! Vasp convention
 !   ipm = 1
 !   gt_cfact = -exp(-gwr%tau_mesh(itau) * eig_nk)
 ! else
 !   ABI_WARNING("Metallic system of semiconductor with Fermi level inside bands!!!!")
 ! end if

 call wrtout(std_out, sjoin(" Building Green's functions from KS states with nbsum", itoa(nbsum), "..."))

 ! TODO
 ! Compute contribution to head and wings for this (kpoint, spin)
 ! taking into account that the same ik_ibz might be replicated.

 !mpw = maxval(wfd%npwarr)
 ! TODO: Mv m_ddk in 72_response below 70_gw or move gwr to higher level.
 !ddkop = ddkop_new(dtset, gwr%cryst, gwr%pawtab, gwr%psps, wfd%mpi_enreg, mpw, wfd%ngfft)
 !call ddkop%free()

 !ABI_MALLOC(cgwork, (2, npw_k * gwr%nspinor))
 !call ddkop%setup_spin_kpoint(gwr%dtset, gwr%cryst, gwr%psps, spin, kk_ibz, istwf_k, npw_k, kg_k)
 !do il_b=1, ugb%sizeb_local(2)
 !  band = ugb%loc2gcol(il_b)
 !  eig_nk = gwr%ks_ebands%eig(band, ik_ibz, spin)
 !  call c_f_pointer(c_loc(ugb%buffer_cplx(:,il_b)), cwave, shape=[2, npwsp])
 !  call ddkop%apply(eig_nk, npw_k, nspinor, cwave, cwaveprj)
 !end do
 !call ddkop%free()
 !ABI_FREE(cgwork)

 ! Init ugb(npwsp, nbsum) PBLAS matrix within the gtau communicator.
 ! and distribute it over bands so that each proc reads a subset of bands in read_band_block

 ABI_MALLOC(gwr%ugb, (gwr%nkibz, gwr%nsppol))
 gwr%ugb_nband = nbsum

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     npw_k = gwr%green_desc_kibz(ik_ibz)%npw
     npwsp = npw_k * gwr%nspinor
     ABI_CHECK(block_dist_1d(gwr%ugb_nband, gwr%gtau_comm%nproc, col_bsize, msg), msg)
     call gwr%ugb(ik_ibz, spin)%init(npwsp, gwr%ugb_nband, gwr%gtau_slkproc, istwfk1, size_blocs=[npwsp, col_bsize])
   end do
 end do

 call wfk_open_read(wfk, wfk_path, formeig0, iomode_from_fname(wfk_path), get_unit(), gwr%comm) !  hdr_out=wfk_hdr)

 mpw = maxval(wfk_hdr%npwarr)
 ABI_MALLOC(kg_k, (3, mpw))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     call cwtime(cpu_green, wall_green, gflops_green, "start")
     ik_ibz = gwr%my_kibz_inds(my_iki)
     kk_ibz = gwr%kibz(:, ik_ibz)
     npw_k = wfk_hdr%npwarr(ik_ibz)
     istwf_k = wfk_hdr%istwfk(ik_ibz)
     ! TODO
     ABI_CHECK_IEQ(istwf_k, 1, "istwfk_k should be 1")
     !nband_k = wfk_hdr%nband(ik_ibz)
     npwsp = npw_k * gwr%nspinor

     associate (ugb => gwr%ugb(ik_ibz, spin), desc_k => gwr%green_desc_kibz(ik_ibz))

     my_bstart = ugb%loc2gcol(1)
     my_bstop = ugb%loc2gcol(ugb%sizeb_local(2))
     my_nband = my_bstop - my_bstart + 1
     ABI_CHECK(my_nband > 0, "nb == 0, decrease number of procs for G and tau parallelism.")
     !print *, "my_bstart, my_bstop, nsum",  my_bstart, my_bstop, nbsum

     ! Read my bands
     call c_f_pointer(c_loc(ugb%buffer_cplx), cg_k, shape=[2, npwsp * my_nband])

     call wfk%read_band_block([my_bstart, my_bstop], ik_ibz, spin, xmpio_single, & ! xmpio_collective,
                              kg_k=kg_k, cg_k=cg_k)

     ! TODO: use round-robing + use gtau_comm% for IO.
     !call wfk%read_bmask(bmask, ik_ibz, spin, xmpio_single, & ! xmpio_collective,
     !                    kg_k=kg_k, cg_k=cg_k) !, eig_k, occ_k)

     ABI_CHECK_IEQ(npw_k, desc_k%npw, "npw_k != desc_k%npw")
     ABI_CHECK(all(kg_k(:,1:npw_k) == desc_k%gvec), "kg_k != desc_k%gvec")

     write(msg,'(3(a,i0),a)')" My ik_ibz [", my_iki, "/", gwr%my_nkibz, "] (tot: ", gwr%nkibz, ")"
     call cwtime_report(msg, cpu_green, wall_green, gflops_green)
     end associate
   end do ! my_iki
 end do ! my_is

 ABI_FREE(kg_k)
 call wfk%close()
 call wfk_hdr%free()

 !call gwr%build_green(free_ugb=.True.)

 call cwtime_report(" gwr_read_ugb_from_wfk:", cpu, wall, gflops)

end subroutine gwr_read_ugb_from_wfk
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_green
!! NAME
!!  gwr_build_green
!!
!! FUNCTION
!!  Build the Green's function in imaginary time from the WFK file `wfk_path`
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_green(gwr, free_ugb)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 logical,intent(in) :: free_ugb

!Local variables-------------------------------
!scalars
 integer :: my_is, my_iki, spin, ik_ibz, band, itau, ipm, il_b, npwsp, isgn
 real(dp) :: f_nk, eig_nk, cpu, wall, gflops, cpu_green, wall_green, gflops_green
 character(len=500) :: msg
 complex(dp) :: gt_cfact
 type(matrix_scalapack), target :: work, green ! cg_mat,
 integer :: mask_kibz(gwr%nkibz)

! *************************************************************************

 ABI_CHECK(abs(gwr%ks_ebands%fermie) < tol12, "ef should have been set to zero!")

 call cwtime(cpu, wall, gflops, "start")

 ! Allocate Green's functions
 mask_kibz = 0; mask_kibz(gwr%my_kibz_inds(:)) = 1
 call gwr%alloc_free_mats(mask_kibz, "green", "malloc")

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     call cwtime(cpu_green, wall_green, gflops_green, "start")
     ik_ibz = gwr%my_kibz_inds(my_iki)
     associate (ugb => gwr%ugb(ik_ibz, spin), desc_k => gwr%green_desc_kibz(ik_ibz))
     npwsp = desc_k%npw * gwr%nspinor

     call ugb%copy(work)
     !call ugb%change_size_blocs(work, size_blocs=, processor=)
     !call work%copy(green, empty=.True.)

     call green%init(npwsp, npwsp, gwr%gtau_slkproc, istwfk1) ! size_blocs=[npwsp, col_bsize])

     do itau=1,gwr%ntau
       do ipm=1,2
         ! Multiply columns by exponentials in imaginary time.
         work%buffer_cplx = ugb%buffer_cplx

         !!$OMP PARALLEL DO PRIVATE(band, f_nk, eig_nk, gt_cfact)
         do il_b=1, work%sizeb_local(2)
          band = work%loc2gcol(il_b)
          f_nk = gwr%ks_ebands%occ(band, ik_ibz, spin)
          eig_nk = gwr%ks_ebands%eig(band, ik_ibz, spin)
          gt_cfact = zero
          if (ipm == 2) then
            if (eig_nk < -tol6) gt_cfact = exp(gwr%tau_mesh(itau) * eig_nk)
          else
            if (eig_nk > tol6) gt_cfact = exp(-gwr%tau_mesh(itau) * eig_nk)
          end if

          work%buffer_cplx(:,il_b) = work%buffer_cplx(:,il_b) * sqrt(real(gt_cfact))
          !call xscal(npwsp, work%buffer_cplx(:,il_b), cone * sqrt(real(gt_cfact)), 1)
         end do ! il_b

         ! Build G(g,g',ipm)
         isgn = merge(1, -1, ipm == 2)
         call slk_pzgemm("N", "C", work, isgn * cone, work, czero, green)

         ! Here we redistribute the data between two different communicators: gtau_comm -> g_comm.
         call gwr%gt_kibz(ipm, ik_ibz, itau, spin)%take_from(green)
       end do ! ipm
     end do ! itau

     call work%free()
     call green%free()
     if (free_ugb) call ugb%free()

     write(msg,'(3(a,i0),a)')" My ik_ibz [", my_iki, "/", gwr%my_nkibz, "] (tot: ", gwr%nkibz, ")"
     call cwtime_report(msg, cpu_green, wall_green, gflops_green)
     end associate
   end do ! my_iki
 end do ! my_is

 if (gwr%dtset%prtvol > 0) call gwr_print_trace(gwr, "gt_kibz")
 call cwtime_report(" gwr_build_green:", cpu, wall, gflops)

end subroutine gwr_build_green
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_rotate_gpm
!! NAME
!!  gwr_rotate_gpm
!!
!! FUNCTION
!!  Reconstruct the Green's functions in the kBZ from the IBZ.
!!
!! INPUTS
!!   ik_bz = Index of the k-point in the BZ
!!   itau = tau index (global index)
!!   spin = spin index
!!
!! OUTPUT
!!  desc_kbz =
!!  gt_pm(2) =
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

subroutine gwr_rotate_gpm(gwr, ik_bz, itau, spin, desc_kbz, gt_pm)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 integer,intent(in) :: ik_bz, spin, itau
 type(desc_t),intent(out) :: desc_kbz
 type(matrix_scalapack),intent(out) :: gt_pm(2)

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, ipm, ik_ibz, isym_k, trev_k, g0_k(3), tsign_k
 logical :: isirr_k
!arrays
 integer :: g1(3), g2(3)
 real(dp) :: tnon(3) !, cpu, wall, gflops
 complex(dp) :: ph2, ph1

! *************************************************************************

 !call cwtime(cpu, wall, gflops, "start")

 ik_ibz = gwr%kbz2ibz(1, ik_bz); isym_k = gwr%kbz2ibz(2, ik_bz)
 trev_k = gwr%kbz2ibz(6, ik_bz); g0_k = gwr%kbz2ibz(3:5, ik_bz)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 tsign_k = merge(1, -1, trev_k == 0)
 !ABI_CHECK(all(g0_k == 0), sjoin("For kbz:", ktoa(gwr%kbz(:, ik_bz)), "g0_k:", ltoa(g0_k), " != 0"))

 ! Copy descriptor from IBZ
 associate (desc_kibz => gwr%green_desc_kibz(ik_ibz))
 call desc_kibz%copy(desc_kbz)

 if (isirr_k) then
   ! Copy the PBLAS matrices with the two Green's functions and we are done.
   do ipm=1,2
     call gwr%gt_kibz(ipm, ik_ibz, itau, spin)%copy(gt_pm(ipm))
   end do
   goto 10
 end if

 ! From:
 !
 !      u_{Sk}(Sg) = e^{-i(Sk+g).tnon} u_k(g)
 !
 ! and
 !
 !      u_{k+g0}(g-g0) = u_k(g)
 !
 ! one obtains:
 !
 !      G_{Sk+g0}(Sg-g0,Sg'-g0) = e^{-i tnon.S(g-g')} G_k{g,g'}
 !
 ! For time-reversal, we have u_{-k}(g) = u_{k}{-g}^*
 !
 !      G_{-k}(-g,-g') = [G_{k}(g,g')]*

 !ABI_WARNING_IF(trev_k == 0, "green: trev_k /= 0 not yet coded")

 ! Rotate gvec, recompute gbound and rotate vc_sqrt
 ! TODO: 1) Handle TR and routine to rotate tchi/W including vc_sqrt
 !       2) Make sure that the FFT box is large enough to accomodate umklapps
 do ig1=1,desc_kbz%npw
   desc_kbz%gvec(:,ig1) = tsign_k * matmul(gwr%cryst%symrec(:,:,isym_k), desc_kibz%gvec(:,ig1)) - g0_k
 end do

 call sphereboundary(desc_kbz%gbound, desc_kbz%istwfk, desc_kbz%gvec, gwr%g_mgfft, desc_kbz%npw)

 ! Get G_k with k in the BZ.
 tnon = gwr%cryst%tnons(:, isym_k)
 !print *, "tnon:", tnon
 do ipm=1,2
   associate (gk_i => gwr%gt_kibz(ipm, ik_ibz, itau, spin), gk_f => gt_pm(ipm))
   call gk_i%copy(gk_f)
   do il_g2=1, gk_f%sizeb_local(2)
     ig2 = mod(gk_f%loc2gcol(il_g2) - 1, desc_kbz%npw) + 1
     g2 = desc_kbz%gvec(:,ig2)
     !g2 = desc_kibz%gvec(:,ig2)
     ph2 = exp(+j_dpc * two_pi * dot_product(g2, tnon))
     do il_g1=1, gk_f%sizeb_local(1)
       ig1 = mod(gk_f%loc2grow(il_g1) - 1, desc_kbz%npw) + 1
       g1 = desc_kbz%gvec(:,ig1)
       !g1 = desc_kibz%gvec(:,ig1)
       ph1 = exp(-j_dpc * two_pi * dot_product(g1, tnon))
       gk_f%buffer_cplx(il_g1, il_g2) = gk_i%buffer_cplx(il_g1, il_g2) * ph1 * ph2
       if (trev_k == 1) gk_f%buffer_cplx(il_g1, il_g2) = conjg(gk_f%buffer_cplx(il_g1, il_g2))
     end do
   end do
   !TODO
   !if (trev_k == 1) then ! TR
   !call gkf%ptrans_ip("C", -cone)
   end associate
 end do
 end associate

10 continue
 !call cwtime_report(" gwr_rotate_gpm:", cpu, wall, gflops)

end subroutine gwr_rotate_gpm
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_myk_green_gpr
!! NAME
!!  gwr_get_myk_green_gpr
!!
!! FUNCTION
!!  Use FFTs to compute G_k(g,g') --> G_k(g',r)
!!  for each k in the BZ treated by this MPI proc for given spin and tau.
!!
!!  1) FFT Transform the first index and multiply by e^{ik.r}:
!!
!!          G_k(g,g') --> G_k(r,g') = e^{ik.r} \sum_g e^{ig.r} G_k(g,g')
!!
!!     NB: This is a local operation.
!!
!!  2) MPI transpose the matrix to go from (r,g') to (g',r) distribution.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_myk_green_gpr(gwr, itau, spin, desc_mykbz, gt_gpr)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 integer,intent(in) :: itau, spin
 type(desc_t),intent(out) :: desc_mykbz(gwr%my_nkbz)
 type(matrix_scalapack),intent(inout) :: gt_gpr(2, gwr%my_nkbz)

!Local variables-------------------------------
!scalars
 integer :: my_ikf, ik_bz, ig2, ipm, npwsp, col_bsize, idat, ndat !, ii
 logical :: k_is_gamma
 real(dp) :: kk_bz(3), cpu, wall, gflops
 character(len=500) :: msg
 type(matrix_scalapack) :: rgp, gt_pm(2)
 complex(dp),allocatable :: ceikr(:)

! *************************************************************************

 if (gwr%timeit) call cwtime(cpu, wall, gflops, "start")

 ABI_MALLOC(ceikr, (gwr%g_nfft * gwr%nspinor))

 do my_ikf=1,gwr%my_nkbz
   ik_bz = gwr%my_kbz_inds(my_ikf)
   kk_bz = gwr%kbz(:, ik_bz)

   k_is_gamma = normv(kk_bz, gwr%cryst%gmet, "G") < GW_TOLQ0
   if (.not. k_is_gamma) call calc_ceikr(kk_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ceikr)

   ! Get G_k(+/- itau) in the BZ.
   call gwr%rotate_gpm(ik_bz, itau, spin, desc_mykbz(my_ikf), gt_pm)

   do ipm=1,2
     associate (ggp => gt_pm(ipm), desc_k => desc_mykbz(my_ikf))

     ! Allocate rgp PBLAS matrix to store G(r,g')
     npwsp = desc_k%npw * gwr%nspinor
     ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
     call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_k%istwfk, &
                   size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

     ABI_CHECK_IEQ(size(ggp%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

     do ig2=1, ggp%sizeb_local(2), gwr%uc_batch_size
       ndat = blocked_loop(ig2, ggp%sizeb_local(2), gwr%uc_batch_size)

       !ABI_CHECK_IEQ(size(ggp%buffer_cplx(:, ig2)), desc_k%npw * gwr%nspinor, "npw")
       !ABI_CHECK_IEQ(size(rgp%buffer_cplx(:, ig2)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

       ! Perform FFT G_k(g,g') -> G_k(r,g') and store results in rgp.
       call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                   gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                   ggp%buffer_cplx(:, ig2), &  ! in
                   rgp%buffer_cplx(:, ig2))    ! out

       ! Multiply by e^{ik.r}
       if (.not. k_is_gamma) then
         !$OMP PARALLEL DO
         do idat=0,ndat-1
           rgp%buffer_cplx(:, ig2 + idat) = ceikr(:) * rgp%buffer_cplx(:, ig2 + idat)
         end do
       end if
     end do ! ig2

     ! MPI transpose: G_k(r,g') -> G_k(g',r)
     call rgp%ptrans("N", gt_gpr(ipm, my_ikf))
     call rgp%free()
     end associate
   end do ! ipm

   call slk_array_free(gt_pm)
 end do ! my_ikf

 ABI_FREE(ceikr)

 if (gwr%timeit) call cwtime_report(" gwr_get_myk_green_gpr:", cpu, wall, gflops)

end subroutine gwr_get_myk_green_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_gk_rpr_pm
!! NAME
!!  gwr_get_gk_rpr_pm
!!
!! FUNCTION
!!  Use FFTs to compute G_k(r',r) from G_k(g,g') for k in the BZ and given spin and tau.
!!  Note that output matrix is transposed i.e. (r',r) instead of (r,r')
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_gk_rpr_pm(gwr, my_ikf, itau, spin, gk_rpr_pm)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 integer,intent(in) :: my_ikf, itau, spin
 type(matrix_scalapack),intent(inout) :: gk_rpr_pm(2)

!Local variables-------------------------------
!scalars
 integer :: ik_bz, ig2, ipm, npwsp, col_bsize, ir1, ndat
 real(dp) :: cpu, wall, gflops
 type(matrix_scalapack) :: rgp, gt_pm(2), gpr
 type(desc_t) :: desc_k
 character(len=500) :: msg

! *************************************************************************

 if (gwr%timeit) call cwtime(cpu, wall, gflops, "start")

 ik_bz = gwr%my_kbz_inds(my_ikf)
 !kk_bz = gwr%kbz(:, ik_bz)

 ! Get G_k(+/- itau) in the BZ.
 call gwr%rotate_gpm(ik_bz, itau, spin, desc_k, gt_pm)

 do ipm=1,2
   ! Allocate rgp PBLAS matrix to store G(r,g')
   npwsp = desc_k%npw * gwr%nspinor
   ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
   call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_k%istwfk, &
                 size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

   associate (ggp => gt_pm(ipm))
   !ABI_CHECK_IEQ(size(ggp%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

   do ig2=1,ggp%sizeb_local(2),gwr%uc_batch_size
     ndat = blocked_loop(ig2, ggp%sizeb_local(2), gwr%uc_batch_size)

     !ABI_CHECK_IEQ(size(ggp%buffer_cplx(:, ig2)), desc_k%npw * gwr%nspinor, "npw")
     !ABI_CHECK_IEQ(size(rgp%buffer_cplx(:, ig2)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

     ! Perform FFT G_k(g,g') -> G_k(r,g') and store results in rgp.
     call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                 gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                 ggp%buffer_cplx(:, ig2), &  ! in
                 rgp%buffer_cplx(:, ig2))    ! out

   end do ! ig2
   end associate

   ! MPI transpose: G_k(r,g') -> G_k(g',r)
   call rgp%ptrans("N", gpr)
   call rgp%free()

   do ir1=1,gpr%sizeb_local(2),gwr%uc_batch_size
     ndat = blocked_loop(ir1, gpr%sizeb_local(2), gwr%uc_batch_size)

     !ABI_CHECK_IEQ(size(gpr%buffer_cplx(:, ir1)), desc_k%npw * gwr%nspinor, "npw")
     !ABI_CHECK_IEQ(size(gk_rpr_pm(ipm)%buffer_cplx(:, ir1)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

     ! FIXME: FFT sign is wrong (should be - instead of + but I need to change the API)
     ! Perform FFT G_k(g',r) -> G_k(r',r) and store results in rgp.
     call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                 gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                 gpr%buffer_cplx(:, ir1), &             ! in
                 gk_rpr_pm(ipm)%buffer_cplx(:, ir1))    ! out
   end do
   call gpr%free()

   ! Rescale
   !gk_rpr_pm(ipm)%buffer_cplx = gk_rpr_pm(ipm)%buffer_cplx * gwr%g_nfft
 end do ! ipm

 call slk_array_free(gt_pm)
 call desc_k%free()

 if (gwr%timeit) call cwtime_report(" gwr_get_gk_rpr_pm:", cpu, wall, gflops)

end subroutine gwr_get_gk_rpr_pm
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_ggp_to_rpr
!! NAME
!!  gwr_ggp_to_rpr
!!
!! FUNCTION
!!  F_{g,g'} --> F_{r',r}
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_ggp_to_rpr(gwr, desc, ggp, rpr)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 type(desc_t),intent(in) :: desc
 type(matrix_scalapack),intent(in) :: ggp
 type(matrix_scalapack),intent(inout) :: rpr

!Local variables-------------------------------
!scalars
 integer :: ig2, npwsp, col_bsize, ir1, ndat
 type(matrix_scalapack) :: rgp, gpr
 character(len=500) :: msg

! *************************************************************************

 ! Allocate rgp PBLAS matrix to store F(r,g')
 npwsp = desc%npw * gwr%nspinor
 ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
 call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc%istwfk, &
               size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

 do ig2=1,ggp%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ig2, ggp%sizeb_local(2), gwr%uc_batch_size)
   call fft_ug(desc%npw, gwr%g_nfft, gwr%nspinor, ndat, &
               gwr%g_mgfft, gwr%g_ngfft, desc%istwfk, desc%gvec, desc%gbound, &
               ggp%buffer_cplx(:, ig2), & ! in
               rgp%buffer_cplx(:, ig2))   ! out
 end do ! ig2

 ! F(r,g') --> F(g',r)
 call rgp%ptrans("N", gpr)
 call rgp%free()

 ! F(g',r) -> F(r',r) and store results in rpr
 do ir1=1,gpr%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ir1, gpr%sizeb_local(2), gwr%uc_batch_size)
   call fft_ur(desc%npw, gwr%g_nfft, gwr%nspinor, ndat, &
               gwr%g_mgfft, gwr%g_ngfft, desc%istwfk, desc%gvec, desc%gbound, &
               gpr%buffer_cplx(:, ir1), &  ! in
               rpr%buffer_cplx(:, ir1))    ! out
 end do

 call gpr%free()

end subroutine gwr_ggp_to_rpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_rotate_wc
!! NAME
!!  gwr_rotate_wc
!!
!! FUNCTION
!!  Reconstruct Wc(q) in the BZ from the IBZ.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_rotate_wc(gwr, iq_bz, itau, spin, desc_qbz, wc_qbz)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 integer,intent(in) :: iq_bz, itau, spin
 type(desc_t),intent(out) :: desc_qbz
 type(matrix_scalapack),intent(inout) :: wc_qbz

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, iq_ibz, isym_q, trev_q, g0_q(3), tsign_q
 logical :: isirr_q
!arrays
 integer :: g1(3), g2(3)
 real(dp) :: tnon(3), qq_bz(3)
 complex(dp) :: ph2, ph1

! *************************************************************************

 ABI_CHECK(gwr%wc_space == "itau", sjoin("wc_space:", gwr%wc_space, " != itau"))

 !spin = gwr%my_spins(my_is)
 !itau = gwr%my_itaus(my_it)
 !iq_bz = gwr%my_qbz_inds(my_iqf)
 qq_bz = gwr%qbz(:, iq_bz)

 iq_ibz = gwr%qbz2ibz(1, iq_bz); isym_q = gwr%qbz2ibz(2, iq_bz)
 trev_q = gwr%qbz2ibz(6, iq_bz); g0_q = gwr%qbz2ibz(3:5, iq_bz)
 isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
 tsign_q = merge(1, -1, trev_q == 0)
 ! TODO: Understand why legacy GW does not need umklapp
 !ABI_CHECK(all(g0_q == 0), sjoin("For qbz:", ktoa(gwr%qbz(:, iq_bz)), "g0_q:", ltoa(g0_q), " != 0"))

 ! Copy descriptor from IBZ
 associate (desc_qibz => gwr%tchi_desc_qibz(iq_ibz))
 call desc_qibz%copy(desc_qbz)

 if (isirr_q) then
   ! Copy the matrices and we are done.
   call gwr%wc_qibz(iq_ibz, itau, spin)%copy(wc_qbz)
   return
 end if

 !ABI_WARNING_IF(trev_q == 0, "wc_rotate: trev_q /= 0 not yet coded")

 ! rotate gvec, recompute gbound and rotate vc_sqrt.
 ! TODO: 1) Handle TR and routine to rotate tchi/W including vc_sqrt
 !       2) Make sure that FFT box is large enough to accomodate umklapps
 do ig1=1,desc_qbz%npw
   desc_qbz%gvec(:,ig1) = tsign_q * matmul(gwr%cryst%symrec(:,:,isym_q), desc_qibz%gvec(:,ig1)) - g0_q
 end do

 call sphereboundary(desc_qbz%gbound, desc_qbz%istwfk, desc_qbz%gvec, gwr%g_mgfft, desc_qbz%npw)

 ! TODO: rotate vc_sqrt
 ! vc(Sq, Sg) = vc(q, g)
 ! vc(-q, -g) = vc(q, g)
 do ig1=1,desc_qbz%npw
   desc_qbz%vc_sqrt(ig1) = sqrt(four_pi) / normv(qq_bz + desc_qbz%gvec(:,ig1), gwr%cryst%gmet, "G")
 end do

 ! Get Wc_q with q in the BZ.
 tnon = gwr%cryst%tnons(:, isym_q)
 associate (wq_i => gwr%wc_qibz(iq_ibz, itau, spin), wq_f => wc_qbz)
 call wq_i%copy(wc_qbz)
 do il_g2=1, wq_f%sizeb_local(2)
   ig2 = mod(wq_f%loc2gcol(il_g2) - 1, desc_qbz%npw) + 1
   g2 = desc_qbz%gvec(:,ig2)
   !g2 = desc_qibz%gvec(:,ig2)
   !ph2 = exp(-j_dpc * two_pi * dot_product(g2, tnon))
   ph2 = exp(+j_dpc * two_pi * dot_product(g2, tnon))
   do il_g1=1, wq_f%sizeb_local(1)
     ig1 = mod(wq_f%loc2grow(il_g1) - 1, desc_qbz%npw) + 1
     g1 = desc_qbz%gvec(:,ig1)
     !g1 = desc_qibz%gvec(:,ig1)
     !ph1 = exp(+j_dpc * two_pi * dot_product(g1, tnon))
     ph1 = exp(-j_dpc * two_pi * dot_product(g1, tnon))
     wq_f%buffer_cplx(il_g1, il_g2) = wq_i%buffer_cplx(il_g1, il_g2) * ph1 * ph2
     if (trev_q == 1) wq_f%buffer_cplx(il_g1, il_g2) = conjg(wq_f%buffer_cplx(il_g1, il_g2))
   end do
 end do
 end associate
 end associate

end subroutine gwr_rotate_wc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_myq_wc_gpr
!! NAME
!!  gwr_get_myq_wc_gpr
!!
!! FUNCTION
!!  Use FFTs to compute: Wc_q(g,g') --> Wc_q(g',r)
!!  for each q in the BZ treated by this MPI proc for given spin and tau.
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

subroutine gwr_get_myq_wc_gpr(gwr, itau, spin, desc_myqbz, wc_gpr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 integer,intent(in) :: itau, spin
 type(desc_t),target,intent(out) :: desc_myqbz(gwr%my_nqbz)
 type(matrix_scalapack),intent(inout) :: wc_gpr(gwr%my_nqbz)

!Local variables-------------------------------
!scalars
 integer :: my_iqf, iq_bz, ig2, npwsp, col_bsize, idat, ndat
 real(dp) :: cpu, wall, gflops, qq_bz(3)
 logical :: q_is_gamma
 character(len=500) :: msg
 type(matrix_scalapack) :: rgp, wc_qbz
 complex(dp),allocatable :: ceiqr(:)

! *************************************************************************

 if (gwr%timeit) call cwtime(cpu, wall, gflops, "start")

 ABI_MALLOC(ceiqr, (gwr%g_nfft * gwr%nspinor))

 do my_iqf=1,gwr%my_nqbz
   iq_bz = gwr%my_qbz_inds(my_iqf)
   qq_bz = gwr%qbz(:, iq_bz)

   q_is_gamma = normv(qq_bz, gwr%cryst%gmet, "G") < GW_TOLQ0
   if (.not. q_is_gamma) call calc_ceikr(qq_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ceiqr)

   ! Get W_q in the BZ.
   call gwr%rotate_wc(iq_bz, itau, spin, desc_myqbz(my_iqf), wc_qbz)
   associate (desc_q => desc_myqbz(my_iqf))

   ! Allocate rgp PBLAS matrix to store Wc(r, g')
   npwsp = desc_q%npw * gwr%nspinor
   ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
   call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_q%istwfk, &
                 size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

   !ABI_CHECK_IEQ(size(wc_qbz%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

   ! FFT. Results stored in rgp
   do ig2=1,wc_qbz%sizeb_local(2), gwr%uc_batch_size
     ndat = blocked_loop(ig2, wc_qbz%sizeb_local(2), gwr%uc_batch_size)

     !ABI_CHECK_IEQ(size(wc_qbz%buffer_cplx(:, ig2)), desc_q%npw, "npw")
     !ABI_CHECK_IEQ(size(rgp%buffer_cplx(:, ig2)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

     call fft_ug(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                 gwr%g_mgfft, gwr%g_ngfft, desc_q%istwfk, desc_q%gvec, desc_q%gbound, &
                 wc_qbz%buffer_cplx(:, ig2), &  ! in
                 rgp%buffer_cplx(:, ig2))        ! out

     ! Multiply by e^{iq.r}
     if (.not. q_is_gamma) then
       !$OMP PARALLEL DO
       do idat=0,ndat-1
         rgp%buffer_cplx(:, ig2+idat) = ceiqr(:) * rgp%buffer_cplx(:, ig2+idat)
       end do
     end if
   end do ! ig2

   ! MPI transposition: Wc(r,g') -> Wc(g',r)
   call rgp%ptrans("N", wc_gpr(my_iqf))
   call rgp%free()
   end associate

   call wc_qbz%free()
 end do ! my_iqf

 ABI_FREE(ceiqr)

 if (gwr%timeit) call cwtime_report(" gwr_get_myq_wc_gpr:", cpu, wall, gflops)

end subroutine gwr_get_myq_wc_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_wc_rpr_qbz
!! NAME
!!  gwr_get_wc_rpr_qbz
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_wc_rpr_qbz(gwr, qq_bz, itau, spin, wc_rpr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 real(dp),intent(in) :: qq_bz(3)
 integer,intent(in) :: itau, spin
 type(matrix_scalapack),intent(inout) :: wc_rpr

!Local variables-------------------------------
!scalars
 integer :: iq_bz, ig2, npwsp, col_bsize, ir1, ndat !, idat, g0_q(3), my_iqf,
 !logical :: q_is_gamma
 character(len=500) :: msg
 type(desc_t) :: desc_q
 type(matrix_scalapack) :: wc_ggp, rgp, gpr
 complex(dp),allocatable :: ceiqr(:)

! *************************************************************************

 ABI_ERROR("Not Implemented error")

 call wrtout(std_out, trim(ktoa(qq_bz)))

 ! Identify q + g0 = k - kp with q in gwr%qbz.
 ! NB: Non-zero g0, requires the application of the phase.
 !call findqg0(iq_bz, g0_q, kmkp, gwr%nqbz, gwr%qbz, gwr%mg0)

 !iq_bz = gwr%my_qbz_inds(my_iqf)
 !qq_bz = gwr%qbz(:, iq_bz)
 !q_is_gamma = normv(qq_bz, gwr%cryst%gmet, "G") < GW_TOLQ0
 !if (.not. q_is_gamma) then
 !ABI_MALLOC(ceiqr, (gwr%g_nfft * gwr%nspinor))
 !call calc_ceikr(qq_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ceiqr)
 !end if

 ! Get W_q in the BZ.
 call gwr%rotate_wc(iq_bz, itau, spin, desc_q, wc_ggp)

 ! Allocate rgp PBLAS matrix to store Wc(r, g')
 npwsp = desc_q%npw * gwr%nspinor
 ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
 call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_q%istwfk, &
               size_blocs=[gwr%g_nfft * gwr%nspinor, col_bsize])

 !ABI_CHECK_IEQ(size(wc_ggp%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

 ! FFT. Results stored in rgp
 do ig2=1,wc_ggp%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ig2, wc_ggp%sizeb_local(2), gwr%uc_batch_size)

   !ABI_CHECK_IEQ(size(wc_ggp%buffer_cplx(:, ig2)), desc_q%npw, "npw")
   !ABI_CHECK_IEQ(size(rgp%buffer_cplx(:, ig2)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

   call fft_ug(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat, &
               gwr%g_mgfft, gwr%g_ngfft, desc_q%istwfk, desc_q%gvec, desc_q%gbound, &
               wc_ggp%buffer_cplx(:, ig2), &  ! in
               rgp%buffer_cplx(:, ig2))       ! out

   ! Multiply by e^{iq.r}
   !if (.not. q_is_gamma) then
   !  do idat=0,ndat-1
   !    !rgp%buffer_cplx(:, ig2+idat) = ceiqr(:) * rgp%buffer_cplx(:, ig2+idat)
   !  end do
   !end if
 end do ! ig2

 ABI_FREE(ceiqr)

 ! MPI transpose: Wc(r,g') -> Wc(g',r)
 call rgp%ptrans("N", gpr)
 call rgp%free()

 do ir1=1,gpr%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ir1, gpr%sizeb_local(2), gwr%uc_batch_size)

   !ABI_CHECK_IEQ(size(gpr%buffer_cplx(:, ir1)), desc_q%npw * gwr%nspinor, "npw")
   !ABI_CHECK_IEQ(size(gk_rpr_pm(ipm)%buffer_cplx(:, ir1)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

   ! FIXME: FFT sign is wrong (should be - instead of + but I need to change the FFT API)
   ! Perform FFT G_k(g',r) -> G_k(r',r) and store results in rgp.
   call fft_ug(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat, &
               gwr%g_mgfft, gwr%g_ngfft, desc_q%istwfk, desc_q%gvec, desc_q%gbound, &
               gpr%buffer_cplx(:, ir1), &     ! in
               wc_rpr%buffer_cplx(:, ir1))    ! out
 end do

 call gpr%free()
 call desc_q%free()
 call wc_ggp%free()

end subroutine gwr_get_wc_rpr_qbz
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
 integer :: ndat, idat, loc1_size, loc2_size, batch_size
 real(dp) :: cpu, wall, gflops !, tau
 logical :: sum_spins_
!arrays
 real(dp), pointer :: weights_ptr(:,:)
 complex(dp) :: wgt_globmy(gwr%ntau, gwr%my_ntau)  ! Use complex instead of real to be able to use ZGEMM.
 complex(dp),allocatable :: cwork_myit(:,:,:), glob_cwork(:,:,:)
 type(matrix_scalapack), pointer :: mats(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 sum_spins_ = .False.; if (present(sum_spins)) sum_spins_ = sum_spins

 call wrtout(std_out, sjoin(" Performing cosine transform. what:", what, "mode:", mode))

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
   weights_ptr => gwr%cosft_tw

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
   weights_ptr => gwr%cosft_wt

 case default
   ABI_ERROR(sjoin("Wrong mode:", mode))
 end select

 ! Extract my weights.
 do my_it=1,gwr%my_ntau
   itau = gwr%my_itaus(my_it)
   do iw=1,gwr%ntau
     wgt_globmy(iw, my_it) = weights_ptr(iw, itau)
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

     ! Use the first itau index to get the size of the local buffer.
     ! Block over ig2 to reduce the number of MPI calls and take advantage of ZGEMM.
     it0 = gwr%my_itaus(1)
     loc1_size = mats(it0)%sizeb_local(1)
     loc2_size = mats(it0)%sizeb_local(2)
     !
     ! batch_size in terms of columns
     ! TODO: How to determine batch_size automatically to avoid going OOM
     !batch_size = loc2_size
     batch_size = 4
     !batch_size = 1

     ABI_MALLOC(cwork_myit, (gwr%my_ntau, loc1_size, batch_size))
     ABI_MALLOC(glob_cwork, (gwr%ntau, loc1_size, batch_size))

     do ig2=1,mats(it0)%sizeb_local(2), batch_size
       ndat = blocked_loop(ig2, mats(it0)%sizeb_local(2), batch_size)

       ! Extract matrix elements as a function of tau.
       do idat=1,ndat
         do ig1=1,mats(it0)%sizeb_local(1)
           do my_it=1,gwr%my_ntau
             itau = gwr%my_itaus(my_it)
             cwork_myit(my_it, ig1, idat) = mats(itau)%buffer_cplx(ig1, ig2+idat-1)
           end do
         end do
       end do

       ! Compute contribution to itau matrix
       !!$OMP PARALLEL DO
       do idat=1,ndat
         !do ig1=1,mats(it0)%sizeb_local(1)
         ! do itau=1,gwr%ntau
         !   glob_cwork(itau, ig1, idat) = dot_product(wgt_globmy(itau, :), cwork_myit(:, ig1, idat))
         ! end do
         !end do
         call ZGEMM("N", "N", gwr%ntau, loc1_size, gwr%my_ntau, cone, &
                    wgt_globmy, gwr%ntau, cwork_myit(:,:,idat), gwr%my_ntau, czero, glob_cwork(:,:,idat), gwr%ntau)
       end do

       call xmpi_sum(glob_cwork, gwr%tau_comm%value, ierr)

       ! Update my local (g1,g2) entry to have it in imaginary-frequency space.
       do idat=1,ndat
         do ig1=1,mats(it0)%sizeb_local(1)
           do my_it=1,gwr%my_ntau
             itau = gwr%my_itaus(my_it)
             mats(itau)%buffer_cplx(ig1, ig2+idat-1) = glob_cwork(itau, ig1, idat)
           end do
         end do
       end do

     end do ! ig2

     ABI_FREE(cwork_myit)
     ABI_FREE(glob_cwork)

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
                gwr%tchi_qibz(iq_ibz,itau,spin+1)%buffer_cplx = mats(itau)%buffer_cplx
              end if
            end if
          end if

        end do ! my_it
      end do ! my_is
   end do ! my_iqi
 end if

 call cwtime_report(" gwr_cos_transform:", cpu, wall, gflops)

end subroutine gwr_cos_transform
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_init
!! NAME
!!  desc_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine desc_init(desc, kk, istwfk, ecut, gwr)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc
 real(dp),intent(in) :: kk(3)
 integer,intent(in) :: istwfk
 real(dp),intent(in) :: ecut
 class(gwr_t),intent(in) :: gwr

! *************************************************************************

 desc%istwfk = istwfk
 call get_kg(kk, desc%istwfk, ecut, gwr%cryst%gmet, desc%npw, desc%gvec)
 ABI_MALLOC(desc%gbound, (2 * gwr%g_mgfft + 8, 2))
 call sphereboundary(desc%gbound, desc%istwfk, desc%gvec, gwr%g_mgfft, desc%npw)

end subroutine desc_init
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
!!integer, allocatable :: gvec2gsort(:)
!! (npw)
!! Mapping the gvec array and the sorted one.
!! Computed by calling desc_calc_gnorm_table. Mainly used to compare data with legacy code
!!real(dp), allocatable :: sorted_gnorm(:)
!! (npw)
!! Sorted list with the norm of the gvector
!! Computed by calling desc_calc_gnorm_table.
!!
!! SOURCE

subroutine desc_calc_gnorm_table(desc, cryst, gvec2gsort, sorted_gnorm)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc
 class(crystal_t),intent(in) :: cryst
 integer, allocatable,intent(out) :: gvec2gsort(:)
 real(dp),allocatable,intent(out) :: sorted_gnorm(:)

!Local variables-------------------------------
 integer :: ig, isort
 integer, allocatable :: gsort2gvec(:)

! *************************************************************************

 ABI_MALLOC(sorted_gnorm, (desc%npw))
 ABI_MALLOC(gvec2gsort, (desc%npw))
 ABI_MALLOC(gsort2gvec, (desc%npw))

 do ig=1,desc%npw
   sorted_gnorm(ig) = normv(desc%gvec(:,ig), cryst%gmet, "G")
   gsort2gvec(ig) = ig
 end do

 call sort_dp(desc%npw, sorted_gnorm, gsort2gvec, tol16)

 ! Invert the mapping.
 do isort=1,desc%npw
   ig = gsort2gvec(isort)
   gvec2gsort(ig) = isort
 end do

 ABI_FREE(gsort2gvec)

end subroutine desc_calc_gnorm_table
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_calc_ginv_map
!! NAME
!!  desc_calc_ginv_map
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!  ginv: index of -g
!!
!! SOURCE

subroutine desc_calc_ginv_map(desc, ginv)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc
 integer, allocatable, intent(out) :: ginv(:)

!Local variables-------------------------------
 integer :: nmiss
 integer, allocatable :: inv_gvec(:,:)

! *************************************************************************

 ABI_MALLOC(inv_gvec, (3, desc%npw))
 inv_gvec = -desc%gvec

 ABI_MALLOC(ginv, (desc%npw))
 call kg_map(desc%npw, desc%gvec, desc%npw, inv_gvec, ginv, nmiss)
 ABI_FREE(inv_gvec)

 ABI_CHECK(nmiss == 0, "nmiss > 0, Cannot find -g in gvec, perhaps q != 0")

end subroutine desc_calc_ginv_map
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_free
!! NAME
!!  desc_free
!!
!! FUNCTION
!!  Free memory
!!
!! SOURCE

subroutine desc_free(desc)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc

! *************************************************************************

 ABI_SFREE(desc%gvec)
 ABI_SFREE(desc%gbound)
 ABI_SFREE(desc%vc_sqrt)

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

subroutine print_gsort_(mat, desc, cryst, times_ucvol)

!Arguments ------------------------------------
 class(matrix_scalapack),intent(inout) :: mat
 class(desc_t),intent(inout) :: desc
 class(crystal_t),intent(in) :: cryst
 logical,intent(in) :: times_ucvol

!Local variables-------------------------------
 integer :: il_g1, il_g2, ig1, ig2, iglob1, iglob2, is1, is2
 real(dp) :: fact
 integer, allocatable :: gvec2gsort(:)
 real(dp), allocatable :: sorted_gnorm(:)
 complex(dp),allocatable :: cwork(:,:)

! *************************************************************************

 call desc%calc_gnorm_table(cryst, gvec2gsort, sorted_gnorm)

 ! FIXME: It won't work if nspinor == 2
 ABI_MALLOC(cwork, (mat%sizeb_local(1), mat%sizeb_local(2)))
 cwork = huge(one)
 fact = one; if (times_ucvol) fact = cryst%ucvol

 do il_g2=1,mat%sizeb_local(2)
   iglob2 = mat%loc2gcol(il_g2)
   ig2 = mod(iglob2 -1 , desc%npw) + 1
   is2 = gvec2gsort(ig2)
   do il_g1=1,mat%sizeb_local(1)
     iglob1 = mat%loc2grow(il_g1)
     ig1 = mod(iglob1 - 1, desc%npw) + 1
     is1 = gvec2gsort(ig1)
     cwork(is1, is2) = mat%buffer_cplx(il_g1, il_g2) * fact
   end do
 end do

 ABI_FREE(gvec2gsort)
 ABI_FREE(sorted_gnorm)

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
 call ydoc%add_int1d("g_ngfft", gwr%g_ngfft(1:6))
 call ydoc%add_int1d("P gwr_np_gtks", [gwr%g_comm%nproc, gwr%tau_comm%nproc, gwr%kpt_comm%nproc, gwr%spin_comm%nproc])
 call ydoc%add_int1d("P np_kibz", gwr%np_kibz)
 call ydoc%add_int1d("P np_qibz", gwr%np_qibz)
 ! Print Max error due to inhomogeneous FT.
 call ydoc%add_real("ft_max_err_t2w_cos", gwr%ft_max_error(1))
 call ydoc%add_real("ft_max_err_w2t_cos", gwr%ft_max_error(2))
 call ydoc%add_real("ft_max_err_t2w_sin", gwr%ft_max_error(3))
 call ydoc%add_real("cosft_duality_error", gwr%cosft_duality_error)

 ! Print imaginary time/frequency mesh with weights.
 call ydoc%open_tabular("Minimax imaginary tau/omega mesh", comment="tau, weight(tau), omega, weight(omega)")
 do ii=1,gwr%ntau
   write(msg, "(i0, 4(es12.5,2x))")ii, gwr%tau_mesh(ii), gwr%tau_wgs(ii), gwr%iw_mesh(ii), gwr%iw_wgs(ii)
   call ydoc%add_tabular_line(msg)
 end do

 call ydoc%write_and_free(unt)

 if (gwr%dtset%prtvol >= 30) then
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

!!****f* m_gwr/gwr_print_mem
!! NAME
!!  gwr_print_mem
!!
!! FUNCTION
!!  Print memory allocated for matrices.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_print_mem(gwr, unit)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
!scalars
 integer :: unt
 real(dp) :: mem_mb
 !character(len=500) :: msg

! *********************************************************************

 unt = std_out; if (present(unit)) unt =unit

 if (allocated(gwr%gt_kibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%gt_kibz))
   call wrtout(std_out, sjoin("- Local memory for Green's functions: ", ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
 end if
 if (allocated(gwr%tchi_qibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%tchi_qibz))
   call wrtout(std_out, sjoin('- Local memory for tchi_qibz: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
 end if
 if (allocated(gwr%wc_qibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%wc_qibz))
   call wrtout(std_out, sjoin('- Local memory for wc_iqbz: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
 end if
 if (allocated(gwr%sigc_kibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%sigc_kibz))
   call wrtout(std_out, sjoin('- Local memory for sigc_kibz: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
 end if

end subroutine gwr_print_mem
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_tchi
!! NAME
!!  gwr_build_cthi
!!
!! FUNCTION
!!  Compute the irreducible polarizability.
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
 integer :: my_is, my_it, my_ikf, ig, my_ir, my_nr, npwsp, ncol_glob, col_bsize, my_iqi !, my_iki ! my_iqf,
 integer :: idat, ndat, sc_nfft, spin, ik_bz, iq_ibz, ierr, ipm, itau, ig2, my_rank ! ig1
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all, cpu_ir, wall_ir, gflops_ir
 real(dp) :: tchi_rfact, mem_mb, local_max, max_abs_imag_chit !, spin_fact, sc_ucvol, ik_ibz,
 logical :: q_is_gamma
 character(len=5000) :: msg
 type(desc_t),pointer :: desc_k, desc_q
 type(matrix_scalapack) :: chi_rgp
 !type(mpi_type) :: tchi_mpi_enreg
!arrays
 integer :: sc_ngfft(18)
 integer,allocatable :: green_scg(:,:), chi_scg(:,:)
 integer :: mask_qibz(gwr%nqibz)
 real(dp) :: kk_bz(3), kpq_bz(3), qq_ibz(3) !, qq_bz(3) ! kk_ibz(3),
 !logical :: need_gt_kibz(gwr%nkibz), got_gt_kibz(gwr%nkibz)
 complex(dp),allocatable :: gt_scbox(:,:), chit_scbox(:), gt_ucbox(:,:)
 !type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz), chiq_gpr(gwr%my_nqibz), gk_rpr_pm(2)
 !type(desc_t), target :: desc_mykbz(gwr%my_nkbz)
 type(matrix_scalapack),allocatable :: gt_gpr(:,:), chiq_gpr(:), gk_rpr_pm(:)
 type(desc_t), target, allocatable :: desc_mykbz(:)
 type(fftbox_plan3_t) :: plan_gp2rp, plan_rp2gp
 complex(dp),allocatable :: cemiqr(:), sc_ceimkr(:,:) !, uc_ceikr(:)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 my_rank = xmpi_comm_rank(gwr%comm)

 ABI_CHECK(gwr%tchi_space == "none", sjoin("tchi_space: ", gwr%tchi_space, " != none"))
 gwr%tchi_space = "itau"

 ! Allocate tchi matrices
 mask_qibz = 0; mask_qibz(gwr%my_qibz_inds(:)) = 1
 call gwr%alloc_free_mats(mask_qibz, "tchi", "malloc")

 max_abs_imag_chit = zero

 if (gwr%use_supercell_for_tchi) then
   ! ============================
   ! GWr algorithm with supercell
   ! ============================

   ! Set FFT mesh in the supercell
   sc_ngfft = gwr%g_ngfft
   sc_ngfft(1:3) = gwr%ngkpt * gwr%g_ngfft(1:3)
   sc_ngfft(4:6) = sc_ngfft(1:3)
   sc_nfft = product(sc_ngfft(1:3))
   !sck_ucvol = gwr%cryst%ucvol * product(gwr%ngkpt)
   !scq_ucvol = gwr%cryst%ucvol * product(gwr%ngqpt)
   !sc_augsize = product(sc_ngfft(4:6))

   !call initmpi_seq(tchi_mpi_enreg)
   !call init_distribfft_seq(tchi_mpi_enreg%distribfft, 'c', sc_ngfft(2), sc_ngfft(3), 'all')
   !call init_distribfft_seq(tchi_mpi_enreg%distribfft, 'f', sc_ngfft(2), sc_ngfft(3), 'all')

   call wrtout(std_out, " Building chi0 with FFTs in the supercell:", pre_newlines=2)
   call wrtout(std_out, sjoin(" ngkpt:", ltoa(gwr%ngkpt), " ngqpt:", ltoa(gwr%ngqpt)))
   call wrtout(std_out, sjoin(" gwr_boxcutmin:", ftoa(gwr%dtset%gwr_boxcutmin)))
   call wrtout(std_out, sjoin(" sc_ngfft:", ltoa(sc_ngfft(1:8))))
   call wrtout(std_out, sjoin(" my_ntau:", itoa(gwr%my_ntau), "ntau:", itoa(gwr%ntau)))
   call wrtout(std_out, sjoin(" my_nkbz:", itoa(gwr%my_nkbz), "nkibz:", itoa(gwr%nkibz)))
   call wrtout(std_out, sjoin(" FFT uc_batch_size:", itoa(gwr%uc_batch_size), &
                              " FFT sc_batch_size:", itoa(gwr%sc_batch_size)), do_flush=.True.)

   ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
   ! Perhaps the safest approach would be to generate the plan on the fly.
   ndat = gwr%sc_batch_size
   ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat, 2))
   ABI_CALLOC(chit_scbox, (sc_nfft * gwr%nspinor * ndat))

   ! This for the version with my_nkbz FFTs in the unit cell
   ABI_CALLOC(gt_ucbox, (gwr%g_nfft * gwr%nspinor * ndat, 2))

   ! (ngfft, ndat, isign)
   call plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat * 2, -1)
   call plan_rp2gp%from_ngfft(sc_ngfft, gwr%nspinor * ndat,     +1)

   !call plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat * 2, +1)
   !call plan_rp2gp%from_ngfft(sc_ngfft, gwr%nspinor * ndat,     -1)

   ! The g-vectors in the supercell for G and tchi.
   ABI_MALLOC(green_scg, (3, gwr%green_mpw))
   ABI_MALLOC(chi_scg, (3, gwr%tchi_mpw))

   ! The phase e^{-iq.r} in the unit cell.
   ABI_MALLOC(cemiqr, (gwr%g_nfft * gwr%nspinor))

   ! Precompute phase factors e^{-ik.R} in the super cell.
   !ABI_MALLOC(sc_ceimkr, (sc_nfft * gwr%nspinor, gwr%my_nkbz))
   !ABI_MALLOC(uc_ceimkr, (gwr%g_nfft * gwr%nspinor, gwr%my_nkbz))

   !do my_ikf=1,gwr%my_nkbz
   !  ik_bz = gwr%my_kbz_inds(my_ikf)
   !  !call calc_sc_ceikr(-gwr%kbz(:,ik_bz), gwr%ngkpt, gwr%g_ngfft, gwr%nspinor, sc_ceimkr(:,my_ikf))
   !  !call calc_ceikr(-gwr%kbz(:,ik_bz), gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, uc_ceimkr(:,my_ikf))
   !end do
   !ABI_FREE(uc_ceimkr)

   ABI_MALLOC(gt_gpr, (2, gwr%my_nkbz))
   ABI_MALLOC(chiq_gpr, (gwr%my_nqibz))
   ABI_MALLOC(gk_rpr_pm, (2))
   ABI_MALLOC(desc_mykbz, (gwr%my_nkbz))

   ! Allocate PBLAS arrays for tchi_q(g',r) for all q in the IBZ treated by this MPI rank.
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     npwsp = gwr%tchi_desc_qibz(iq_ibz)%npw * gwr%nspinor
     ncol_glob = gwr%g_nfft * gwr%nspinor
     ABI_CHECK(block_dist_1d(ncol_glob, gwr%g_comm%nproc, col_bsize, msg), msg)
     call chiq_gpr(my_iqi)%init(npwsp, gwr%g_nfft * gwr%nspinor, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
   end do
   mem_mb = sum(slk_array_locmem_mb(chiq_gpr))
   call wrtout(std_out, sjoin(' Local memory for chiq_gpr: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

   ! Loop over my spins and my taus.
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
       itau = gwr%my_itaus(my_it)

       ! G_k(g,g') --> G_k(g',r) e^{ik.r} for each k in the BZ treated by me.
       call gwr%get_myk_green_gpr(itau, spin, desc_mykbz, gt_gpr)

       ! Loop over r in the unit cell that is now MPI-distributed in g_comm.
       my_nr = gt_gpr(1,1)%sizeb_local(2)

       do my_ir=1, my_nr, gwr%sc_batch_size
         if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
           call cwtime(cpu_ir, wall_ir, gflops_ir, "start")
         end if
         ndat = blocked_loop(my_ir, my_nr, gwr%sc_batch_size)

         !call cwtime(cpu, wall, gflops, "start")
         ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
         ! Note that we need to take the union of (k, g') for k in the BZ.

         !call xscal(size(gt_scbox), czero, gt_scbox, 1)
         gt_scbox = zero

!if (.False.) then
if (.True.) then
         do my_ikf=1,gwr%my_nkbz
           ik_bz = gwr%my_kbz_inds(my_ikf)

           ! Compute k+g
           desc_k => desc_mykbz(my_ikf)
           do ig=1,desc_k%npw
             green_scg(:,ig) = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ig)
           end do

           do ipm=1,2
             !ABI_CHECK_IEQ(size(gt_gpr(ipm, my_ikf)%buffer_cplx, dim=2), my_nr, "my_nr!")
             !ABI_CHECK_IEQ(size(gt_gpr(ipm, my_ikf)%buffer_cplx, dim=1), desc_k%npw, "desc_k!")

             call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat, green_scg, &
                           gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                           gt_scbox(:, ipm))                             ! inout
           end do
         end do ! my_ikf

         ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
         if (gwr%kpt_comm%nproc > 1) call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)

         !call cwtime_report("G part", cpu, wall, gflops)

         ! G(G',r) --> G(R',r) = sum_{k,g'} e^{-i(k+g').R'} G_k(g',r)
         call plan_gp2rp%execute_ip_dpc(gt_scbox)
         gt_scbox = gt_scbox * sc_nfft

else
         ! This is just to check if replacing a single FFT in the supercell with
         ! nkpt FFTs in the unit cell + communication is faster.
         do ipm=1,2
           do my_ikf=1,gwr%my_nkbz
             desc_k => desc_mykbz(my_ikf)

             ! FIXME The sign is wrong as it should be -1.
             call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                         gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                         gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                         gt_ucbox(:,ipm))                              ! out

             ! Sketching the algo.
             !gt_ucbox(:,ipm) = gt_ucbox(:,ipm) * uc_ceimkr(:,my_ikf)
             !gt_scbox(:,ipm) = matmul(gt_ucbox(:,:,ipm), emiLk)

             !call ctimes_eikr(-kk, gwr%g_ngfft, gwr%g_nfft, ndat, gt_ucbox(:,ipm)))

             ! This becomes the bottleneck but perhaps one can take advantage of localization.
             ! Moreover the FFTs are distributed inside kpt_comm
             ! Also, one can save all the FFTs in a matrix G(mnfft * ndat, my_nkbz) multiply by the e^{-ikr} phase
             ! and then use zgemm to compute Out(r,L) = [e^{-ikr}G_k(r)] e^{-ikL} with precomputed e^{-iLk} phases.
             !call ur_to_scpsi(gwr%ngkpt, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ndat, &
             !                 gt_ucbox(:,ipm), cone, sc_ceimkr(:, my_ikf), gt_scbox(:,ipm))
           end do ! my_ikf
         end do ! ipm

         ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
         if (gwr%kpt_comm%nproc > 1) call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)
         gt_scbox = gt_scbox / sc_nfft
end if

         ! Compute tchi(R',r) for this r. Note that results are real so one might use r2c
         chit_scbox(:) = gt_scbox(:, 1) * conjg(gt_scbox(:, 2))
         !max_abs_imag_chit = max(max_abs_imag_chit, maxval(abs(aimag(chit_scbox))))

         ! Back to tchi(G'=q+g',r) space immediately with isign + 1.
         call plan_rp2gp%execute(chit_scbox)

         ! The GG part is the hotspot
         !call cwtime_report("GG part", cpu, wall, gflops)

         ! Extract tchi_q(g',r) on the ecuteps g-sphere from the FFT box in the supercell
         ! and save data in chiq_gpr PBLAS matrix. Only my q-points in the IBZ are considered.
         ! Alternatively, one can avoid the FFT above and
         ! use zero-padded to go from the supercell to the ecuteps g-sphere inside the my_iqi loop.
         ! This approach should play well with k-point parallelism.
         do my_iqi=1,gwr%my_nqibz
           iq_ibz = gwr%my_qibz_inds(my_iqi)
           qq_ibz = gwr%qibz(:, iq_ibz)
           desc_q => gwr%tchi_desc_qibz(iq_ibz)

           ! Compute q+g vectors.
           do ig=1,desc_q%npw
             chi_scg(:,ig) = nint(qq_ibz * gwr%ngqpt) + gwr%ngqpt(:) * desc_q%gvec(:,ig)
           end do

           !ABI_CHECK_IEQ(size(chiq_gpr(my_iqi)%buffer_cplx, dim=1), desc_q%npw, "desc_q!")
           call box2gsph(sc_ngfft, desc_q%npw, gwr%nspinor * ndat, chi_scg, &
                         chit_scbox, &                            ! in
                         chiq_gpr(my_iqi)%buffer_cplx(:, my_ir))  ! out
         end do ! my_iqi
         !call cwtime_report("chiq part", cpu, wall, gflops)

         if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
           write(msg,'(3(a,i0),a)')" My ir [", my_ir, "/", my_nr, "] (tot: ", gwr%g_nfft, ")"
           if (my_rank == 0) call cwtime_report(msg, cpu_ir, wall_ir, gflops_ir)
         end if
       end do ! my_ir
       call wrtout(std_out, "Computation of chiq_gpr completed")

       ! Free descriptors and PBLAS matrices in kBZ.
       call desc_array_free(desc_mykbz)
       call slk_array_free(gt_gpr)

       ! Now we have tchi_q(g',r).
       ! For each IBZ q-point treated by this MPI proc, do:
       !
       !     1) MPI transpose to have tchi_q(r,g')
       !     2) FFT along first dimension to have tchi_q(g,g') stored in gwr%tchi_qibz
       !
       tchi_rfact = one / gwr%g_nfft / gwr%cryst%ucvol / (gwr%nkbz * gwr%nqbz)
       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         desc_q => gwr%tchi_desc_qibz(iq_ibz)
         q_is_gamma = normv(gwr%qibz(:,iq_ibz), gwr%cryst%gmet, "G") < GW_TOLQ0

         ! Note minus sign in q.
         if (.not. q_is_gamma) call calc_ceikr(-gwr%qibz(:,iq_ibz), gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, cemiqr)

         ! MPI-transposition: tchi_q(g',r) => tchi_q(r,g')
         call chiq_gpr(my_iqi)%ptrans("N", chi_rgp)

         !ABI_CHECK_IEQ(size(gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx, dim=2), size(chi_rgp%buffer_cplx, dim=2), "len2")

         ! FFT r --> g along the first dimension: tchi_q(r,g') --> tchi_q(g,g')
         ! Results stored in gwr%tchi_qibz.
         do ig2=1,chi_rgp%sizeb_local(2), gwr%uc_batch_size
           ndat = blocked_loop(ig2, chi_rgp%sizeb_local(2), gwr%uc_batch_size)

           !ABI_CHECK_IEQ(size(chi_rgp%buffer_cplx(:, ig2)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")
           !ABI_CHECK_IEQ(size(gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2)), desc_q%npw, "npw")
           if (.not. q_is_gamma) then
             !$OMP PARALLEL DO
             do idat=0,ndat-1
               chi_rgp%buffer_cplx(:, ig2 + idat) = cemiqr(:) * chi_rgp%buffer_cplx(:, ig2 + idat)
             end do
           end if

           call fft_ur(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat, gwr%g_mgfft, gwr%g_ngfft, &
                       istwfk1, desc_q%gvec, desc_q%gbound, &
                       chi_rgp%buffer_cplx(:, ig2),         &                  ! ur(in)
                       gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2))  ! ug(out)

           !$OMP PARALLEL DO
           do idat=0,ndat-1
             gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2 + idat) = &
             gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2 + idat) * tchi_rfact
             !call xscal(size(gt_scbox), czero, gt_scbox, 1)
           end do
         end do

         call chi_rgp%free()
         call gwr%tchi_qibz(iq_ibz, itau, spin)%set_imag_diago_to_zero(local_max)
       end do ! my_iqi

       write(msg,'(3(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
       call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
     end do ! my_it
   end do ! my_is

   ABI_FREE(gt_ucbox)
   ABI_FREE(gt_scbox)
   ABI_FREE(chit_scbox)
   ABI_FREE(green_scg)
   ABI_FREE(chi_scg)
   ABI_FREE(cemiqr)
   !ABI_SFREE(sc_ceimkr)

   call slk_array_free(chiq_gpr)
   !call destroy_mpi_enreg(tchi_mpi_enreg)

   ABI_FREE(gt_gpr)
   ABI_FREE(chiq_gpr)
   ABI_FREE(gk_rpr_pm)
   ABI_FREE(desc_mykbz)

  else
    ! ===================================================================
    ! Mixed-space algorithm in the unit cell with convolutions in k-space
    ! ===================================================================

    do my_is=1,gwr%my_nspins
      spin = gwr%my_spins(my_is)
      do my_it=1,gwr%my_ntau
        call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
        itau = gwr%my_itaus(my_it)

        ! Redistribute G_k(g,g') in the IBZ so that each MPI proc can recostruct G in the BZ on the fly.
        !need_gt_kibz = ?
        !call gwr%distrib_gt_kibz(itau, spin, need_kibz, got_kibz, "communicate')

        do my_ikf=1,gwr%my_nkbz
          ik_bz = gwr%my_kbz_inds(my_ikf)
          kk_bz = gwr%kbz(:, ik_bz)

          ! Use symmetries to get G_k from the IBZ. G_k(g,g') --> G_k(r,r')
          call gwr%get_gk_rpr_pm(my_ikf, my_it, my_is, gk_rpr_pm)

          do iq_ibz=1,gwr%nqibz
            !do my_iqi=1,gwr%my_nqibz
            !iq_ibz = gwr%my_qibz_inds(my_iqi)
            qq_ibz = gwr%qibz(:, iq_ibz)
            !iq_bz = gwr%my_qbz_inds(my_iqf)
            !qq_bz = gwr%qbz(:, iq_bz)
            kpq_bz = kk_bz + qq_ibz
            !ikpq_ibz = ??
            !kmq_bz = kk_bz - qq_ibz
            !ikmq_ibz = ??
            !call gwr%get_gk_rpr_pm(my_ikf, my_it, my_is, gk_rpr_pm)
            !gk_rpr_pm(1)%buffer_cplx * conjg(gkq_rpr_pm(2)%buffer_cplx)
          end do ! iq_ibz
        end do ! my_ikf

        ! Deallocate extra G's
        !call gwr%distrib_gt_kibz(itau, spin, need_kibz, got_kibz, "free")

        call slk_array_free(gk_rpr_pm)
        !gwr%tchi_qibz(iq_ibz, itau, spin)

       write(msg,'(3(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
       call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
     end do ! my_it
   end do ! spin
 end if

 !call wrtout(std_out, sjoin(" max_abs_imag_chit", ftoa(max_abs_imag_chit)))

 ! Print trace of chi_q(i tau) matrices for testing purposes.
 if (gwr%dtset%prtvol > 0) call gwr%print_trace("tchi_qibz")

 ! Transform irreducible tchi from imaginary tau to imaginary omega.
 ! Also sum over spins to get total tchi if collinear spin.
 call gwr%cos_transform("tchi", "it2w", sum_spins=.True.)

 call load_head_wings_from_sus_file__(gwr, "AW_CD/runo_DS3_SUS.nc")

 ! Print trace of chi_q(i omega) matrices for testing purposes.
 if (gwr%dtset%prtvol > 0) call gwr%print_trace("tchi_qibz")

 ! Write file with chi0(i omega)
 if (gwr%dtset%prtsuscep > 0) call gwr%ncwrite_tchi_wc("tchi", trim(gwr%dtfil%filnam_ds(4))//'_TCHIM.nc')

 call cwtime_report(" gwr_build_tchi:", cpu_all, wall_all, gflops_all)

end subroutine gwr_build_tchi
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_distrib_gt_kibz
!! NAME
!!  gwr_distrib_gt_kibz
!!
!! FUNCTION
!!  Redistribute G_k for fixed (itau, spin) according to need_kibz
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_distrib_gt_kibz(gwr, itau, spin, need_kibz, got_kibz, action)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 integer,intent(in) :: itau, spin, need_kibz(gwr%nkibz)
 integer,intent(out) :: got_kibz(gwr%nkibz)
 character(len=*),intent(in) :: action

!Local variables-------------------------------
 integer :: ik_ibz, ipm, ierr, lsize(2) ! col_bsize, npwsp,
 integer :: do_mpi_kibz(gwr%nkibz), sender_kibz(gwr%nkibz)
 real(dp) :: kk_ibz(3), cpu, wall, gflops
 complex(dp),allocatable :: cbuf_k(:,:)

! *************************************************************************

 if (gwr%timeit) call cwtime(cpu, wall, gflops, "start")

 select case (action)
 case ("communicate")
   do_mpi_kibz = need_kibz
   call xmpi_sum(do_mpi_kibz, gwr%kpt_comm%value, ierr)

   got_kibz = 0
   do ik_ibz=1,gwr%nkibz
     kk_ibz = gwr%kibz(:, ik_ibz)
     sender_kibz(ik_ibz) = huge(1)
     if (allocated(gwr%green_desc_kibz(ik_ibz)%gvec)) sender_kibz(ik_ibz) = gwr%kpt_comm%me
     if (need_kibz(ik_ibz) /= 0 .and. .not. allocated(gwr%green_desc_kibz(ik_ibz)%gvec)) then
       got_kibz(ik_ibz) = 1
       call gwr%green_desc_kibz(ik_ibz)%init(kk_ibz, istwfk1, gwr%dtset%ecut, gwr)
     end if
   end do

   call xmpi_min_ip(sender_kibz, gwr%kpt_comm%value, ierr)

   ! Allocate memory
   call gwr%alloc_free_mats(got_kibz, "green", "malloc")

   ! MPI communication
   do ik_ibz=1,gwr%nkibz
     if (do_mpi_kibz(ik_ibz) == 0) cycle
     lsize = gwr%gt_kibz(1, ik_ibz, itau, spin)%sizeb_local(:)
     ABI_MALLOC(cbuf_k, (lsize(1), lsize(2)))
     do ipm=1,2
       if (gwr%kpt_comm%me == sender_kibz(ik_ibz)) cbuf_k = gwr%gt_kibz(ipm, ik_ibz, itau, spin)%buffer_cplx
       call xmpi_bcast(cbuf_k, sender_kibz(ik_ibz), gwr%kpt_comm%value, ierr)
       if (need_kibz(ik_ibz) == 1) gwr%gt_kibz(ipm, ik_ibz, itau, spin)%buffer_cplx = cbuf_k
     end do
     ABI_FREE(cbuf_k)
   end do

 case ("free")
   ! Use got_kibz to free previously allocated memory
   do ik_ibz=1,gwr%nkibz
     if (got_kibz(ik_ibz) == 1) call gwr%green_desc_kibz(ik_ibz)%free()
   end do
   call gwr%alloc_free_mats(got_kibz, "green", "free")

 case default
   ABI_ERROR(sjoin("Invalid action:", action))
 end select

 if (gwr%timeit) call cwtime_report(" gwr_distrib_gt_kibz:", cpu, wall, gflops)

end subroutine gwr_distrib_gt_kibz
!!***

!!****f* m_gwr/gwr_distrib_mats_qibz
!! NAME
!!  gwr_distrib_mats_qibz
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_distrib_mats_qibz(gwr, what, itau, spin, need_qibz, got_qibz, action)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: what
 integer,intent(in) :: itau, spin, need_qibz(gwr%nqibz)
 integer,intent(out) :: got_qibz(gwr%nqibz)
 character(len=*),intent(in) :: action

!Local variables-------------------------------
 integer :: iq_ibz, ierr, lsize(2)  ! col_bsize, npwsp,
 integer :: do_mpi_qibz(gwr%nqibz), sender_qibz(gwr%nqibz)
 real(dp) :: qq_ibz(3), cpu, wall, gflops
 complex(dp),allocatable :: cbuf_q(:,:)

! *************************************************************************

 ABI_CHECK(what == "tchi" .or. what == "wc", sjoin("Invalid what:", what))

 if (gwr%timeit) call cwtime(cpu, wall, gflops, "start")

 select case (action)
 case ("communicate")
   do_mpi_qibz = need_qibz
   call xmpi_sum(do_mpi_qibz, gwr%kpt_comm%value, ierr)

   got_qibz = 0
   do iq_ibz=1,gwr%nqibz
     qq_ibz = gwr%qibz(:, iq_ibz)
     sender_qibz(iq_ibz) = huge(1)
     if (allocated(gwr%tchi_desc_qibz(iq_ibz)%gvec)) sender_qibz(iq_ibz) = gwr%kpt_comm%me
     if (need_qibz(iq_ibz) /= 0 .and. .not. allocated(gwr%tchi_desc_qibz(iq_ibz)%gvec)) then
       got_qibz(iq_ibz) = 1
       call gwr%tchi_desc_qibz(iq_ibz)%init(qq_ibz, istwfk1, gwr%dtset%ecuteps, gwr)
     end if
   end do

   call xmpi_min_ip(sender_qibz, gwr%kpt_comm%value, ierr)

   ! Allocate memory
   call gwr%alloc_free_mats(got_qibz, what, "malloc")

   ! MPI communication
   do iq_ibz=1,gwr%nqibz
     if (do_mpi_qibz(iq_ibz) == 0) cycle
     lsize = gwr%gt_kibz(1, iq_ibz, itau, spin)%sizeb_local(:)
     ABI_MALLOC(cbuf_q, (lsize(1), lsize(2)))
     if (gwr%kpt_comm%me == sender_qibz(iq_ibz)) then
       if (what == "tchi") cbuf_q = gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx
       if (what == "wc")   cbuf_q = gwr%wc_qibz(iq_ibz, itau, spin)%buffer_cplx
     end if
     call xmpi_bcast(cbuf_q, sender_qibz(iq_ibz), gwr%kpt_comm%value, ierr)
     if (need_qibz(iq_ibz) == 1) then
       if (what == "tchi") gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx = cbuf_q
       if (what == "wc")   gwr%wc_qibz(iq_ibz, itau, spin)%buffer_cplx = cbuf_q
     end if
     ABI_FREE(cbuf_q)
   end do

 case ("free")
   ! Use got_qibz to free previously allocated memory
   do iq_ibz=1,gwr%nqibz
     if (got_qibz(iq_ibz) == 1) call gwr%tchi_desc_qibz(iq_ibz)%free()
   end do
   call gwr%alloc_free_mats(got_qibz, what, "free")

 case default
   ABI_ERROR(sjoin("Invalid action:", action))
 end select

 if (gwr%timeit) call cwtime_report(" gwr_distrib_mats_qibz:", cpu, wall, gflops)

end subroutine gwr_distrib_mats_qibz
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

 comment = "Invalid space!"; units = [std_out, ab_out]

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
       end do
     end do
   end do

   call xmpi_sum_master(ctrace3, 0, gwr%tks_comm%value, ierr)

   if (xmpi_comm_rank(gwr%comm) == master) then
     do spin=1,gwr%nsppol
       call wrtout(units, sjoin(" Trace of:", what, "for spin:", itoa(spin), "for testing purposes:"))
       call wrtout(units, comment, pre_newlines=2)
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
   comment = " (ik_ibz, itau) table"

   call xmpi_sum_master(ctrace4, master, gwr%tks_comm%value, ierr)

   if (xmpi_comm_rank(gwr%comm) == master) then
     do spin=1,gwr%nsppol
       do ipm=1,2
         call wrtout(units, sjoin(" Trace of:", what, "for ipm:", itoa(ipm), ", spin:", itoa(spin), "for testing purposes:"))
         call wrtout(units, comment, newlines=1)
         call print_arr(ctrace4(:,:, ipm, spin), unit=ab_out)
         call print_arr(ctrace4(:,:, ipm, spin), unit=std_out)
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
 integer :: my_iqi, my_it, my_is, iq_ibz, spin, itau, iw, ierr, il_g1, il_g2, ig1, ig2, iglob1, iglob2, my_rank
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_q, wall_q, gflops_q, q0_vol, i_sz, q0sph
 logical :: q_is_gamma, free_tchi
 character(len=5000) :: msg
 complex(dpc) :: vcs_g1, vcs_g2
 type(matrix_scalapack) :: em1
 type(yamldoc_t) :: ydoc
!arrays
 real(dp) :: qq_ibz(3)
 complex(dpc) :: em1_wq(gwr%ntau, gwr%nqibz), eps_wq(gwr%ntau, gwr%nqibz)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call wrtout(std_out, " Building correlated screening Wc from irreducible chi ...", pre_newlines=2)

 my_rank = xmpi_comm_rank(gwr%comm)

 ABI_CHECK(gwr%tchi_space == "iomega", sjoin("tchi_space: ", gwr%tchi_space, " != iomega"))
 ABI_CHECK(gwr%wc_space == "none", sjoin("wc_space: ", gwr%wc_space, " != none"))
 gwr%wc_space = "iomega"

 ! =======================================
 ! Allocate PBLAS arrays for wc_qibz(g,g')
 ! =======================================
 ! Note that we have already summed tchi over spin.
 ! Also, G=0 corresponds to iglob = 1 since only q-points in the IBZ are treated.
 ! This is not true for ther other q-points in the full BZ as we may have a non-zero umklapp g0_q

 ABI_MALLOC(gwr%wc_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))

 free_tchi = .True.
 if (free_tchi) gwr%tchi_space = "none"

 em1_wq = zero; eps_wq = zero

 ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2
 ! Very crude approximation, see also m_vcoul
 q0_vol = (two_pi)**3 / (gwr%nqbz * gwr%cryst%ucvol)
 i_sz = four_pi*7.44*q0_vol**(-two_thirds)

 q0sph = two_pi * (three / (four_pi * gwr%cryst%ucvol * gwr%nqbz)) ** third

 ! If possible, use 2d rectangular grid of processors for diagonalization.
 !call slkproc_4diag%init(gwr%g_comm%value)

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

       ! Build symmetrized RPA epsilon: 1 - Vc^{1/2} chi0 Vc^{1/2}
       associate (wc => gwr%wc_qibz(iq_ibz, itau, spin))
       call gwr%tchi_qibz(iq_ibz, itau, spin)%copy(wc)
       if (free_tchi) call gwr%tchi_qibz(iq_ibz, itau, spin)%free()

       do il_g2=1,wc%sizeb_local(2)
         iglob2 = wc%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,wc%sizeb_local(1)
           iglob1 = wc%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1
           vcs_g1 = desc_q%vc_sqrt(ig1)
           wc%buffer_cplx(il_g1, il_g2) = -wc%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2
           if (iglob1 == iglob2) then
             wc%buffer_cplx(il_g1, il_g2) = one + wc%buffer_cplx(il_g1, il_g2)
             if (iglob1 == 1 .and. iglob2 == 1) then
               ! Store epsilon_{iw, iq_ibz}(0, 0)
               ! Rescale by np_qibz because we will MPI reduce this array.
               eps_wq(itau, iq_ibz) = wc%buffer_cplx(il_g1, il_g2) / gwr%np_qibz(iq_ibz)
             end if
           end if
         end do ! il_g1
       end do ! il_g2

       ! Invert symmetrized epsilon.
       ! NB: PZGETRF requires square block cyclic decomposition along the two axis
       ! hence we neeed to redistribute the data before calling zinvert.
       ! TODO: Can call zhdp
       !call wrtout(std_out, "Printing wc%buffer_cplex before inversion")
       !call print_arr(wc%buffer_cplx, unit=std_out)

       call wc%change_size_blocs(em1) ! processor=slkproc_4diag
       call em1%zinvert()
       !call em1%zdhp_invert("U")
       call wc%take_from(em1)  ! processor=wc%processor)
       call em1%free()

       !call wrtout(std_out, sjoin(" e-1 at q:", ktoa(qq_ibz), "i omega:", ftoa(gwr%iw_mesh(itau) * Ha_eV), "eV"))
       !call print_arr(wc%buffer_cplx, unit=std_out)

       ! Build Wc(q, iw) = [e^{-1}_q(g,g',iw) - delta_{gg'} v_q(g,g') by removing bare vc
       do il_g2=1,wc%sizeb_local(2)
         iglob2 = wc%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,wc%sizeb_local(1)
           iglob1 = wc%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1
           vcs_g1 = desc_q%vc_sqrt(ig1)

           if (iglob1 == 1 .and. iglob2 == 1) then
             ! Store epsilon^{-1}_{iw, iq_ibz}(0, 0). Rescale by np_qibz because we will MPI reduce this array.
             em1_wq(itau, iq_ibz) = wc%buffer_cplx(il_g1, il_g2) / gwr%np_qibz(iq_ibz)
           end if

           ! Subtract exchange part.
           if (iglob1 == iglob2) wc%buffer_cplx(il_g1, il_g2) = wc%buffer_cplx(il_g1, il_g2) - one

           ! Handle divergence in Wc for q --> 0
           if (q_is_gamma .and. (iglob1 == 1 .or. iglob2 == 1)) then
             if (iglob1 == 1 .and. iglob2 == 1) then
               vcs_g1 = sqrt(i_sz); vcs_g2 = sqrt(i_sz)
             else if (iglob1 == 1) then
               !vcs_g1 = (four_pi) ** (three/two) * q0sph ** 2 / two
               vcs_g1 = sqrt(i_sz)
             else if (iglob2 == 1) then
               !vcs_g2 = (four_pi) ** (three/two) * q0sph ** 2 / two
               vcs_g2 = sqrt(i_sz)
             end if
           end if

           wc%buffer_cplx(il_g1, il_g2) = wc%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2 / gwr%cryst%ucvol
         end do ! il_g1
       end do ! il_g2
       end associate

     end do  ! my_it
   end do ! my_is
   end associate

   write(msg,'(2(a,i0),a)')" My iqi [", my_iqi, "/", gwr%my_nqibz, "]"
   call cwtime_report(msg, cpu_q, wall_q, gflops_q)
 end do ! my_iqi

 !call slkproc_4diag%free()

 call xmpi_sum_master(em1_wq, master, gwr%gtk_comm%value, ierr)
 call xmpi_sum_master(eps_wq, master, gwr%gtk_comm%value, ierr)

 if (my_rank == master) then
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
 if (gwr%dtset%prtvol > 0) call gwr%print_trace("wc_qibz")

 ! Write file with Wc(i omega)
 !if (gwr%dtset%prtsuscep > 0) call gwr%ncwrite_tchi_wc("wc", trim(gwr%dtfil%filnam_ds(4))//'_WCIMW.nc')

 ! Cosine transform from iomega to itau to get Wc(i tau)
 call gwr%cos_transform("wc", "iw2t")

 ! Write file with Wc(i tau)
 !if (gwr%dtset%prtsuscep > 0) call gwr%ncwrite_tchi_wc("wc", trim(gwr%dtfil%filnam_ds(4))//'_WCIMT.nc')

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
 integer :: my_is, my_it, spin, ik_ibz, ikcalc_ibz, sc_nfft, sc_size, my_ir, my_nr, iw, idat, ndat, my_rank !my_iki,
 integer :: my_iqf, iq_ibz, iq_bz, itau, ierr, ibc, bmin, bmax, band, nbc ! col_bsize, npwsp, ib1, ib2,
 integer :: my_ikf, ipm, ik_bz, ig, ikcalc, uc_ir, ir, ncid, col_bsize, npwsp, ibeg, iend ! my_iqi, sc_ir,
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all !, cpu, wall, gflops
 real(dp) :: mem_mb, cpu_ir, wall_ir, gflops_ir
 real(dp) :: max_abs_imag_wct, max_abs_re_wct, sck_ucvol, scq_ucvol  !, spin_fact
 !logical :: k_is_gamma !, q_is_gamma
 character(len=500) :: msg
 type(desc_t), pointer :: desc_q, desc_k
 type(yamldoc_t) :: ydoc
!arrays
 integer :: sc_ngfft(18), need_qibz(gwr%nqibz),  got_qibz(gwr%nqibz), units(4)
 integer,allocatable :: green_scg(:,:), wc_scg(:,:)
 real(dp) :: kk_bz(3), kcalc_bz(3), qq_bz(3)  !, qq_ibz(3)
 !real(dp),allocatable :: inv_rmR(:), rylm_rmR(:,:)
 complex(dp) :: sigc_pm(2), cpsi_r, odd_t(gwr%ntau), even_t(gwr%ntau)
 complex(dp),target,allocatable :: sigc_it_diag_kcalc(:,:,:,:,:)
 complex(dp),allocatable :: gt_scbox(:,:), wct_scbox(:), sc_ceikr(:), uc_ceikr(:)
 complex(dp),allocatable :: uc_psi_bk(:,:,:), sc_psi_bk(:,:,:)
 complex(gwpc),allocatable :: ur(:)
 type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz), gk_rpr_pm(2), sigc_rpr(2,gwr%nkcalc), wc_rpr, wc_gpr(gwr%my_nqbz)
 type(desc_t), target :: desc_mykbz(gwr%my_nkbz), desc_myqbz(gwr%my_nqbz)
 type(fftbox_plan3_t) :: gt_plan_gp2rp, wt_plan_gp2rp

 integer :: band_val, ibv, ncerr, unt_it, unt_iw, unt_rw, unt_aw
 real(dp) :: e0, ks_gap, qp_gap
 complex(dp) :: zz, sigc_e0, dsigc_de0, z_e0, sig_xc
 integer,allocatable :: ks_vbik(:,:), iperm(:)
 real(dp),allocatable :: sorted_qpe(:)
 real(dp) :: e0_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 real(dp) :: spfunc_diag_kcalc(gwr%nwr, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 real(dp) :: sigx_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 !real(dp) :: ks_gap_kcalc(gwr%nkcalc, gwr%nsppol), qp_gap_kcalc(gwr%nkcalc, gwr%nsppol)
 complex(dp) :: ze0_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol), rw_mesh(gwr%nwr)
 complex(dp) :: sigc_iw_diag_kcalc(gwr%ntau, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: sigc_e0_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: qpe_zlin_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol), imag_zmesh(gwr%ntau,2)
 complex(dp) :: sigc_rw_diag_kcalc(gwr%nwr, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 my_rank = xmpi_comm_rank(gwr%comm)

 ABI_CHECK(gwr%wc_space == "itau", sjoin("wc_space: ", gwr%wc_space, " != itau"))

 !mask_kibz = 0; mask_kibz(gwr%my_kibz_inds(:)) = 1
 !call gwr%alloc_free_mats(mask_kibz, "sigma" "malloc")

 ! Set FFT mesh in the supercell
 ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
 ! Perhaps the safest approach would be to generate the plan on the fly.

 sc_ngfft = gwr%g_ngfft
 sc_size = product(gwr%ngkpt)
 sc_ngfft(1:3) = gwr%ngkpt * gwr%g_ngfft(1:3)
 sc_ngfft(4:6) = sc_ngfft(1:3)
 sc_nfft = product(sc_ngfft(1:3))
 sck_ucvol = gwr%cryst%ucvol * product(gwr%ngkpt)
 scq_ucvol = gwr%cryst%ucvol * product(gwr%ngqpt)
 !sc_augsize = product(sc_ngfft(4:6))

 call wrtout(std_out, sjoin(" Building Sigma_c with FFTs in the supercell:", ltoa(sc_ngfft(1:3))), pre_newlines=2)
 call wrtout(std_out, sjoin(" ngkpt:", ltoa(gwr%ngkpt), " ngqpt:", ltoa(gwr%ngqpt)))
 call wrtout(std_out, sjoin(" gwr_boxcutmin:", ftoa(gwr%dtset%gwr_boxcutmin)))
 call wrtout(std_out, sjoin(" sc_ngfft:", ltoa(sc_ngfft(1:8))))
 call wrtout(std_out, sjoin(" my_ntau:", itoa(gwr%my_ntau), "ntau:", itoa(gwr%ntau)))
 call wrtout(std_out, sjoin(" my_nkbz:", itoa(gwr%my_nkbz), "nkibz:", itoa(gwr%nkibz)))
 call wrtout(std_out, sjoin(" FFT uc_batch_size:", itoa(gwr%uc_batch_size), &
                            " FFT sc_batch_size:", itoa(gwr%sc_batch_size)), do_flush=.True.)

 ! The g-vectors in the supercell for G and tchi.
 ABI_MALLOC(green_scg, (3, gwr%green_mpw))
 ABI_MALLOC(wc_scg, (3, gwr%tchi_mpw))

 ! Set FFT mesh used to compute u(r) in the unit cell.
 call gwr%kcalc_wfd%change_ngfft(gwr%cryst, gwr%psps, gwr%g_ngfft)

 ! Diagonal matrix elements Sigmac_(itau) in the KS basis set.
 ABI_CALLOC(sigc_it_diag_kcalc, (2, gwr%ntau, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol))
 max_abs_imag_wct = zero; max_abs_re_wct = zero

 call gwr%print_mem(unit=std_out)

if (gwr%use_supercell_for_sigma) then

 ! NOTE:
 ! There are two possibilities here:
 !
 ! 1) Compute the matrix elements of Sigma_c in the KS basis set by integrating over the real-space supercell.
 !
 ! 2) Compute and store the Sigma_c^k(g,g', i omega) and then compute the matrix elements in g-space.
 !
 ! The first option requires less memory provided we are interested in a small set of KS states.
 ! The second option is interesting if we need to compute several matrix elements, including off-diagonal terms.

 ndat = gwr%sc_batch_size
 ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat, 2))
 ABI_CALLOC(wct_scbox, (sc_nfft * gwr%nspinor * ndat))

 ! (ngfft, ndat, isign)
 !call gt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat * 2, +1)
 !call wt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat, +1)

 call gt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat * 2, -1)
 call wt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat, -1)

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   ! Load wavefunctions for GW corrections in the unit cell.
   bmin = minval(gwr%bstart_ks(:, spin))
   bmax = maxval(gwr%bstop_ks(:, spin))

   ABI_MALLOC_OR_DIE(uc_psi_bk, (gwr%g_nfft * gwr%nspinor, bmin:bmax, gwr%nkcalc), ierr)
   ABI_MALLOC_OR_DIE(sc_psi_bk, (sc_nfft * gwr%nspinor, bmin:bmax, gwr%nkcalc), ierr)
   ABI_MALLOC(ur, (gwr%g_nfft * gwr%nspinor))
   ABI_MALLOC(sc_ceikr, (sc_nfft * gwr%nspinor))
   ABI_MALLOC(uc_ceikr, (gwr%g_nfft * gwr%nspinor))

   do ikcalc=1,gwr%nkcalc
     kcalc_bz = gwr%kcalc(:, ikcalc)
     ikcalc_ibz = gwr%kcalc2ibz(ikcalc, 1)  ! TODO: Assuming wfs in IBZ

     ! Compute e^{ik.r} phases on the two FFT meshes: super cell and unit cell.
     call calc_sc_ceikr(kcalc_bz, gwr%ngkpt, gwr%g_ngfft, gwr%nspinor, sc_ceikr)
     call calc_ceikr(kcalc_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, uc_ceikr)

     do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
       call gwr%kcalc_wfd%get_ur(band, ikcalc_ibz, spin, ur)
       uc_psi_bk(:, band, ikcalc) = ur

       ! Multiply by the phases to get psi_nk(r) in the super cell.
       call ur_to_scpsi(gwr%ngkpt, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, 1, &
                        uc_psi_bk(:, band, ikcalc), czero, sc_ceikr, sc_psi_bk(:, band, ikcalc))
       uc_psi_bk(:, band, ikcalc) = uc_psi_bk(:, band, ikcalc) * uc_ceikr
     end do
   end do ! ikcalc

   ABI_FREE(ur)
   ABI_FREE(sc_ceikr)
   ABI_FREE(uc_ceikr)

   ! Construct Sigma(itau) in the supercell.
   do my_it=1,gwr%my_ntau
     call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
     itau = gwr%my_itaus(my_it)

     ! G_k(g,g') --> G_k(g',r) e^{ik.r} for each k in the BZ treated by me.
     call gwr%get_myk_green_gpr(itau, spin, desc_mykbz, gt_gpr)

     ! Wc_q(g,g') --> Wc_q(g',r) e^{iq.r} for each q in the BZ treated by me.
     call gwr%get_myq_wc_gpr(itau, spin, desc_myqbz, wc_gpr)

     my_nr = gt_gpr(1,1)%sizeb_local(2)
     ABI_CHECK(my_nr == wc_gpr(1)%sizeb_local(2), "my_nr != wc_gpr(1)%sizeb_local(2)")

     ! Loop over r in the unit cell that is now MPI-distributed in g_comm.
     do my_ir=1, my_nr, gwr%sc_batch_size
       if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
         call cwtime(cpu_ir, wall_ir, gflops_ir, "start")
       end if
       ndat = blocked_loop(my_ir, my_nr, gwr%sc_batch_size)

       uc_ir = gt_gpr(1,1)%loc2gcol(my_ir)  ! FIXME: This won't work if nspinor 2

       !sc_ir
       !mpsang = gwr%dtser%gwr_lmax + 1
       !ABI_MALLOC(inv_rmR, (sc_nfft))
       !ABI_MALLOC(rylm_rmR, (sc_nfft, mpsang))
       !call inv_distances_and_rylm(uc_pt, mpsang, inv_rmR, rylm_rmR)
       !ABI_FREE(inv_rmR)
       !ABI_FREE(rylm_rmR)

       ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
       ! Take the union of (k,g') for k in the BZ.
       !call cwtime(cpu, wall, gflops, "start")

       gt_scbox = zero
       do my_ikf=1,gwr%my_nkbz
         ik_bz = gwr%my_kbz_inds(my_ikf)
         ik_ibz = gwr%kbz2ibz(1, ik_bz)

         ! Compute k+g'
         desc_k => desc_mykbz(my_ikf)
         do ig=1,desc_k%npw
           green_scg(:,ig) = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ig)
         end do

         do ipm=1,2
           call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat, green_scg, &
                         gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                         gt_scbox(:, ipm))                             ! inout
         end do
       end do ! my_ikf

       ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)

       ! G(G',r) --> G(R',r)
       call gt_plan_gp2rp%execute_ip_dpc(gt_scbox)
       gt_scbox = gt_scbox * (sc_nfft / sck_ucvol)
       !call cwtime_report("G part", cpu, wall, gflops)

       ! Insert Wc_q(g',r) in G'-space in the supercell FFT box for fixed r.
       ! Take the union of (q,g') for q in the BZ. Also, note ngqpt instead of ngkpt.
       wct_scbox = zero
       do my_iqf=1,gwr%my_nqbz
         iq_bz = gwr%my_qbz_inds(my_iqf)
         iq_ibz = gwr%qbz2ibz(1, iq_bz)

         ! Compute q+g' vectors.
         desc_q => desc_myqbz(my_iqf)
         do ig=1,desc_q%npw
           wc_scg(:,ig) = nint(gwr%qbz(:, iq_bz) * gwr%ngqpt) + gwr%ngqpt * desc_q%gvec(:,ig)
         end do

         call gsph2box(sc_ngfft, desc_q%npw, gwr%nspinor * ndat, wc_scg, &
                       wc_gpr(my_iqf)%buffer_cplx(:, my_ir), & ! in
                       wct_scbox)                              ! inout
       end do ! my_iqf

       ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       call xmpi_sum(wct_scbox, gwr%kpt_comm%value, ierr)
       !call cwtime_report("W part", cpu, wall, gflops)

       ! Wc(G',r) --> Wc(R',r)
       call wt_plan_gp2rp%execute_ip_dpc(wct_scbox)
       wct_scbox = wct_scbox * (sc_nfft / scq_ucvol)
       !wct_scbox = wct_scbox * sc_nfft
       !if (itau == 1) print *, "W(r,r):", wct_scbox(uc_ir)

!DEBUG
       !max_abs_re_wct = max(max_abs_re_wct, maxval(abs(real(wct_scbox))))
       !max_abs_imag_wct = max(max_abs_imag_wct, maxval(abs(aimag(wct_scbox))))
       !wct_scbox(uc_ir) = zero
       !wct_scbox = real(wct_scbox)
!END DEBUG

       ! Use gt_scbox to store GW (R',r, +/- i tau) for this r.
       gt_scbox(:, 1) = gt_scbox(:, 1) * wct_scbox(:)
       gt_scbox(:, 2) = gt_scbox(:, 2) * wct_scbox(:)
       !print *, "Maxval abs imag G:", maxval(abs(aimag(gt_scbox)))
       !call cwtime_report("GW part", cpu, wall, gflops)

       ! Integrate self-energy matrix elements in the R-supercell for fixed r(s)
       do ikcalc=1,gwr%nkcalc
         ! FIXME: Temporary hack till I find a better MPI algo for k-points.
         if (gwr%kpt_comm%skip(ikcalc)) cycle
         do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
           ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
           do idat=1,ndat
             ir = uc_ir + idat - 1
             ibeg = 1 + (idat - 1) * sc_nfft * gwr%nspinor
             iend = ibeg + sc_nfft * gwr%nspinor - 1
             cpsi_r = conjg(uc_psi_bk(ir, band, ikcalc))
             sigc_pm(1) = sum(cpsi_r * gt_scbox(ibeg:iend, 1) * sc_psi_bk(:, band, ikcalc))
             sigc_pm(2) = sum(cpsi_r * gt_scbox(ibeg:iend, 2) * sc_psi_bk(:, band, ikcalc))
             sigc_it_diag_kcalc(:, itau, ibc, ikcalc, spin) = &
             sigc_it_diag_kcalc(:, itau, ibc, ikcalc, spin) + sigc_pm(:)
           end  do
          end do
       end do ! ikcalc
       !call cwtime_report("Sig Matrix part", cpu, wall, gflops)

       if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
         write(msg,'(3(a,i0),a)')" My ir [", my_ir, "/", my_nr, "] (tot: ", gwr%g_nfft, ")"
         if (my_rank == 0) call cwtime_report(msg, cpu_ir, wall_ir, gflops_ir)
       end if
     end do ! my_ir

     ! Free descriptors and PBLAS matrices in kBZ and qBZ.
     call desc_array_free(desc_mykbz)
     call slk_array_free(gt_gpr)
     call desc_array_free(desc_myqbz)
     call slk_array_free(wc_gpr)

     write(msg,'(3(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
     call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
   end do ! my_it

   ABI_FREE(sc_psi_bk)
   ABI_FREE(uc_psi_bk)
 end do ! my_is

 !call wrtout(std_out, sjoin(" Maxval abs re W:", ftoa(max_abs_re_wct)))
 !call wrtout(std_out, sjoin(" Maxval abs imag W:", ftoa(max_abs_imag_wct)))

 ABI_FREE(gt_scbox)
 ABI_FREE(wct_scbox)
else

 call wrtout(std_out, sjoin(" Building Sigma_c with convolutions in q-space:", ltoa(sc_ngfft(1:3))), pre_newlines=2)

 ! Allocate Sigma(r',r)
 npwsp = gwr%g_nfft * gwr%nspinor
 col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
 call wc_rpr%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
 do ipm=1,2
   call gk_rpr_pm(ipm)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
   do ikcalc=1,gwr%nkcalc
     call sigc_rpr(ipm, ikcalc)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
   end do
 end do
 mem_mb = slk_array_locmem_mb(wc_rpr) + sum(slk_array_locmem_mb(gk_rpr_pm)) + sum(slk_array_locmem_mb(sigc_rpr))
 call wrtout(std_out, sjoin(' Memory for local PBLAS matrices: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

 ! Define tables to account for symmetries:
 !  - when looping over the BZ, we only need to include the union of IBZ_x for x in kcalc.
 !  - when accumulating the self-energy, we have to include weights that depend on x.

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   ! Load wavefunctions for GW corrections in the unit cell.
   bmin = minval(gwr%bstart_ks(:, spin))
   bmax = maxval(gwr%bstop_ks(:, spin))

   ABI_MALLOC_OR_DIE(uc_psi_bk, (gwr%g_nfft * gwr%nspinor, bmin:bmax, gwr%nkcalc), ierr)
   ABI_MALLOC(ur, (gwr%g_nfft * gwr%nspinor))

   do ikcalc=1,gwr%nkcalc
     kcalc_bz = gwr%kcalc(:, ikcalc)
     ikcalc_ibz = gwr%kcalc2ibz(ikcalc, 1)  ! TODO: Assuming wfs in IBZ

     do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
       call gwr%kcalc_wfd%get_ur(band, ikcalc_ibz, spin, ur)
       uc_psi_bk(:, band, ikcalc) = ur
     end do
   end do

   ABI_FREE(ur)

   ! Construct Sigma(itau) using convolutions in k-space and real-space representation in the unit cell.
   do my_it=1,gwr%my_ntau
     call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
     itau = gwr%my_itaus(my_it)

     ! Redistribute W_q(g,g') in the IBZ so that each MPI proc can recostruct W in the BZ on the fly.
     need_qibz = 0
     call gwr%distrib_mats_qibz("wc", itau, spin, need_qibz, got_qibz, "communicate")
     call slk_array_set(sigc_rpr, czero)

     do my_ikf=1,gwr%my_nkbz
       ik_bz = gwr%my_kbz_inds(my_ikf)
       kk_bz = gwr%kbz(:, ik_bz)

       ! Use symmetries to get G_k from the IBZ. G_k(g,g') --> G_k(r,r')
       call gwr%get_gk_rpr_pm(my_ikf, my_it, my_is, gk_rpr_pm)

       do ikcalc=1,gwr%nkcalc
         qq_bz = gwr%kcalc(:, ikcalc) - kk_bz

         ! FIXME: Finalize implementation
         call gwr%get_wc_rpr_qbz(qq_bz, itau, spin, wc_rpr)

         do ipm=1,2
           sigc_rpr(ipm, ikcalc)%buffer_cplx = sigc_rpr(ipm, ikcalc)%buffer_cplx + &
             gk_rpr_pm(ipm)%buffer_cplx * wc_rpr%buffer_cplx
         end do
       end do ! ikcalc
     end do ! my_ikf

     ! Deallocate extra Wc matrices
     call gwr%distrib_mats_qibz("wc", itau, spin, need_qibz, got_qibz, "free")

     ! Integrate self-energy matrix elements in the unit cell
     ! Remember that sigma is stored as (r', r) with the second dimension MPI-distributed
     do ikcalc=1,gwr%nkcalc
       do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
         !call diag_braket("T", sigc_rpr(1, ikcalc), gwr%g_nfft * gwr%nspinor, uc_psi_bk(:, band, ikcalc), sigc_pm(1))
         !call diag_braket("T", sigc_rpr(2, ikcalc), gwr%g_nfft * gwr%nspinor, uc_psi_bk(:, band, ikcalc), sigc_pm(2))
         ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
         sigc_it_diag_kcalc(:, itau, ibc, ikcalc, spin) = sigc_pm(:)
        end do
     end do ! ikcalc

     write(msg,'(3(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
     call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
   end do ! my_it

   ABI_FREE(uc_psi_bk)
 end do ! my_is

 call wc_rpr%free()
 call slk_array_free(sigc_rpr)
 call slk_array_free(gk_rpr_pm)

end if

 ABI_FREE(green_scg)
 ABI_FREE(wc_scg)

 ! Store matrix elements of Sigma_c(it), separate even and odd part then
 ! use sine/cosine transform to get Sigma_c(i omega).
 ! Finally, perform analytic continuation with Pade' to go to the real-axis,
 ! compute QP corrections and spectral functions.

 sigc_it_diag_kcalc = -sigc_it_diag_kcalc * (gwr%cryst%ucvol / gwr%g_nfft) ** 2 ! / (gwr%nkbz ** 2)

 call xmpi_sum(sigc_it_diag_kcalc, gwr%comm, ierr)

 ! TODO: Fake values for the exchange part.
 sigx_kcalc = zero
 sigc_iw_diag_kcalc = zero

 imag_zmesh(:,1) = j_dpc * gwr%iw_mesh
 imag_zmesh(:,2) = -j_dpc * gwr%iw_mesh

 do spin=1,gwr%nsppol
   do ikcalc=1,gwr%nkcalc
      ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
      do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
        ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
        e0 = gwr%ks_ebands%eig(band, ik_ibz, spin)

        associate (vals_pmt => sigc_it_diag_kcalc(:,:, ibc, ikcalc, spin))
        ! f(t) = E(t) + O(t) = (f(t) + f(-t)) / 2  + (f(t) - f(-t)) / 2
        even_t = (vals_pmt(1,:) + vals_pmt(2,:)) / two
        odd_t = (vals_pmt(1,:) - vals_pmt(2,:)) / two
        sigc_iw_diag_kcalc(:, ibc, ikcalc, spin) = matmul(gwr%cosft_wt, even_t) + j_dpc * matmul(gwr%sinft_wt, odd_t)
        end associate

        ! if zz in 2 or 3 quadrant, avoid branch cut in the complex plane using Sigma(-iw) = Sigma(iw)*.
        zz = cmplx(e0, zero)
        if (real(zz) > zero) then
        !if (real(zz) < zero) then
          sigc_e0 = pade(gwr%ntau, imag_zmesh(:,1), sigc_iw_diag_kcalc(:, ibc, ikcalc, spin), zz)
          dsigc_de0 = dpade(gwr%ntau, imag_zmesh(:,1), sigc_iw_diag_kcalc(:, ibc, ikcalc, spin), zz)
        else
          sigc_e0 = pade(gwr%ntau, imag_zmesh(:,2), conjg(sigc_iw_diag_kcalc(:, ibc, ikcalc, spin)), zz)
          dsigc_de0 = dpade(gwr%ntau, imag_zmesh(:,2), conjg(sigc_iw_diag_kcalc(:, ibc, ikcalc, spin)), zz)
        end if

        ! Z = (1 - dSigma / domega(E0))^{-1}
        z_e0 = one / (one - dsigc_de0)

        ! Store results.
        e0_kcalc(ibc, ikcalc, spin) = e0
        sigc_e0_kcalc(ibc, ikcalc, spin) = sigc_e0
        ze0_kcalc(ibc, ikcalc, spin) = z_e0

        ! Linearized QP solution.
        qpe_zlin_kcalc(ibc, ikcalc, spin) = e0 + &
          z_e0 * (sigc_e0 + sigx_kcalc(ibc, ikcalc, spin) - gwr%ks_me%vxcval(band, band, ik_ibz, spin))

        ! Build linear mesh on the real axis **centered** around e0.
        rw_mesh = arth(e0 - gwr%wr_step * (gwr%nwr / 2), gwr%wr_step, gwr%nwr)
        do iw=1,gwr%nwr
          zz = rw_mesh(iw)
          if (real(zz) > zero) then
          !if (real(zz) < zero) then
            sigc_e0 = pade(gwr%ntau, imag_zmesh(:,1), sigc_iw_diag_kcalc(:, ibc, ikcalc, spin), zz)
          else
            sigc_e0 = pade(gwr%ntau, imag_zmesh(:,2), conjg(sigc_iw_diag_kcalc(:, ibc, ikcalc, spin)), zz)
          end if
          sigc_rw_diag_kcalc(iw, ibc, ikcalc, spin) = sigc_e0
          sig_xc = sigc_e0 + sigx_kcalc(ibc, ikcalc, spin)

          ! TODO: Spectral function
          spfunc_diag_kcalc(iw, ibc, ikcalc, spin) = zero
          !spfunc_diag_kcalc(iw, ibc, ikcalc, spin) = &
          !  one / pi * abs(aimag(sigc_e0) &
          !  /( (real(rw_mesh(iw) - Sr%hhartree(ib, ib, ik_ibz, spin) - sigx_xc)) ** 2 &
          !    +(aimag(sigc_e0)) ** 2) / Ha_eV
        end do

        ! TODO
        !call pade%init(gwr%ntau, imag_zmesh, sigc_iw_diag_kcalc(:, ibc, ikcalc, spin), fermie=zero)
        !call pade%eval(zz, sigc_e0)
        !call pade%eval_der1(zz, dsigc_de0)
        !call pade%eval(rw_mesh, sigc_rw_diag_kcalc(:, ibc, ikcalc, spin))
        !call pade%free()

      end do ! band
   end do ! ikcalc
 end do ! spin

 ! Master writes results to ab_out, std_out and GWR.nc
 if (xmpi_comm_rank(gwr%comm) == 0) then

   !is_metallic = ebands_has_metal_scheme(gwr%ks_ebands)
   !call write_notations([std_out, ab_out]
   ABI_MALLOC(ks_vbik, (gwr%ks_ebands%nkpt, gwr%ks_ebands%nsppol))
   ks_vbik(:,:) = ebands_get_valence_idx(gwr%ks_ebands)

   do spin=1,gwr%nsppol
     do ikcalc=1,gwr%nkcalc
        ik_ibz = gwr%kcalc2ibz(ikcalc, 1)

        ydoc = yamldoc_open('GWR_SelfEnergy_ee', width=11, real_fmt='(3f8.3)')
        call ydoc%add_real1d('kpoint', gwr%kcalc(:, ikcalc))
        call ydoc%add_int('spin', spin, int_fmt="(i1)")

        ! Compute Gaps
        band_val = ks_vbik(ik_ibz, spin)
        nbc = gwr%bstop_ks(ikcalc, spin) - gwr%bstart_ks(ikcalc, spin) + 1

        if (band_val >= gwr%bstart_ks(ikcalc, spin) .and. band_val + 1 <= gwr%bstop_ks(ikcalc, spin)) then
          ibv = band_val - gwr%bstart_ks(ikcalc, spin) + 1
          ks_gap = e0_kcalc(ibv+1, ikcalc, spin) - e0_kcalc(ibv, ikcalc, spin)

          ! This to detect possible band inversion in the QP energies.
          ! and compute qp_gap accordingly.
          call sort_rvals(nbc, real(qpe_zlin_kcalc(:, ikcalc, spin)), iperm, sorted_qpe, tol=tol12)

          if (iperm(ibv) /= ibv .or. iperm(ibv + 1) /= ibv + 1) then
            call ydoc%add_int('QP_VBM_band', iperm(ibv) + gwr%bstart_ks(ikcalc, spin) - 1)
            call ydoc%add_int('QP_CBM_band', iperm(ibv+1) + gwr%bstart_ks(ikcalc, spin) - 1)
            qp_gap = sorted_qpe(ibv+1) - sorted_qpe(ibv)
          else
            call ydoc%add_int('QP_VBM_band', ibv + gwr%bstart_ks(ikcalc, spin) - 1)
            call ydoc%add_int('QP_CBM_band', ibv+1 + gwr%bstart_ks(ikcalc, spin) - 1)
            qp_gap = real(qpe_zlin_kcalc(ibv+1, ikcalc, spin) - qpe_zlin_kcalc(ibv, ikcalc, spin))
          end if
          ABI_FREE(iperm)
          ABI_FREE(sorted_qpe)

          call ydoc%add_real('KS_gap', ks_gap * Ha_eV)
          call ydoc%add_real('QP_gap', qp_gap * Ha_eV)
          call ydoc%add_real('Delta_QP_KS', (qp_gap - ks_gap) * Ha_eV)
        end if

        call ydoc%open_tabular('data') !, tag='SigmaeeData')
        msg = " Band      E0 <VxcDFT>   SigX SigC(E0)     Z1      Z2  E-E0       E"
        call ydoc%add_tabular_line(msg)

        do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
          ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
          e0 = e0_kcalc(ibc, ikcalc, spin)
          write(msg,'(i5, *(f8.3))') &
            band, &                                                          ! Band
            e0 * Ha_eV, &                                                    ! E0
            real(gwr%ks_me%vxcval(band, band, ik_ibz, spin)) * Ha_eV, &      ! <VxcDFT>
            sigx_kcalc(ibc, ikcalc, spin) * Ha_eV, &                         ! SigX
            real(sigc_e0_kcalc(ibc, ikcalc, spin)) * Ha_eV, &                ! SigC(E0)
            real(ze0_kcalc(ibc, ikcalc, spin)), &                            ! Z1
            aimag(ze0_kcalc(ibc, ikcalc, spin)), &                           ! Z2
            (real(qpe_zlin_kcalc(ibc, ikcalc, spin)) - e0) * Ha_eV, &        ! E-E0
            real(qpe_zlin_kcalc(ibc, ikcalc, spin)) * Ha_eV                  ! E
          call ydoc%add_tabular_line(msg)
        end do

        call ydoc%write_units_and_free([std_out, ab_out])
     end do ! ikcalc
   end do ! spin

   ABI_FREE(ks_vbik)

   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_SIGCIT'), msg, newunit=unt_it, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_it, "(a)")"# Diagonal elements of Sigma_c(i tau, +/-) in atomic units"
   write(unt_it, "(a)")"# tau Re/Im Sigma_c(+itau) Re/Im Sigma_c(-itau)"
   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_SIGCIW'), msg, newunit=unt_iw, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_iw, "(a)")"# Diagonal elements of Sigma_c(i omega) in eV units"
   write(unt_iw, "(a)")"# omega Re/Im Sigma_c(i omega)"
   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_SIGCRW'), msg, newunit=unt_rw, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_rw, "(a)")"# Diagonal elements of Sigma_c(omega) in eV units"
   write(unt_rw, "(a)")"# omega Re/Im Sigma_c(omega)"
   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_AW'), msg, newunit=unt_aw, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_aw, "(a)")"# Spectral function A(omega)"

   units = [unt_it, unt_iw, unt_rw, unt_aw]

   ! TODO: Finalize the implementation and put output files under test.
   do spin=1,gwr%nsppol
     do ikcalc=1,gwr%nkcalc
       ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
       call write_units(units, sjoin("# kpt:", ktoa(gwr%kcalc(:, ikcalc)), "spin:", itoa(spin)))
       do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
         ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
         e0 = e0_kcalc(ibc, ikcalc, spin)
         call write_units(units, sjoin("# band:", itoa(band), ", spin:", itoa(spin)))
         ! Write Sigma_c(i tau) and Sigma_c(i omega)
         do itau=1,gwr%ntau
           write(unt_it, "(*(es16.8))") &
             gwr%tau_mesh(itau), &
             c2r(sigc_it_diag_kcalc(1, itau, ibc, ikcalc, spin)), &
             c2r(sigc_it_diag_kcalc(2, itau, ibc, ikcalc, spin))
           write(unt_iw, "(*(es16.8))") gwr%iw_mesh(itau) * Ha_eV, c2r(sigc_iw_diag_kcalc(itau, ibc, ikcalc, spin)) * Ha_eV
         end do
         ! Write Sigma_c(omega) and A(omega)
         rw_mesh = arth(e0 - gwr%wr_step * (gwr%nwr / 2), gwr%wr_step, gwr%nwr) * Ha_eV
         do iw=1,gwr%nwr
           write(unt_rw, "(*(es16.8))") rw_mesh(itau), c2r(sigc_rw_diag_kcalc(iw, ibc, ikcalc, spin)) * Ha_eV
           write(unt_aw, "(*(es16.8))") rw_mesh(itau), spfunc_diag_kcalc(iw, ibc, ikcalc, spin) / Ha_eV
         end do
       end do
     end do
   end do

   close(unt_it); close(unt_iw); close(unt_rw); close(unt_aw)

   ! Add results to the GWR.nc
   NCF_CHECK(nctk_open_modify(ncid, gwr%gwrnc_path, xmpi_comm_self))
   ! Define arrays with results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("e0_kcalc", "dp", "max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ze0_kcalc", "dp", "two, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("qpe_zlin_kcalc", "dp", "two, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("sigx_kcalc", "dp", "max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("sigc_it_diag_kcalc", "dp", "two, two, ntau, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("sigc_iw_diag_kcalc", "dp", "two, ntau, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("sigc_rw_diag_kcalc", "dp", "two, nwr, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("spfunc_diag_kcalc", "dp", "nwr, max_nbcalc, nkcalc, nsppol") &
   ])
   NCF_CHECK(ncerr)
   ! Write data.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("e0_kcalc"), e0_kcalc))
   NCF_CHECK(nf90_put_var(ncid, vid("ze0_kcalc"), c2r(ze0_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("sigx_kcalc"), sigx_kcalc))
   NCF_CHECK(nf90_put_var(ncid, vid("qpe_zlin_kcalc"), c2r(qpe_zlin_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("sigc_it_diag_kcalc"), c2r(sigc_it_diag_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("sigc_iw_diag_kcalc"), c2r(sigc_iw_diag_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("sigc_rw_diag_kcalc"), c2r(sigc_rw_diag_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("spfunc_diag_kcalc"), spfunc_diag_kcalc))
   NCF_CHECK(nf90_close(ncid))
 end if

 call xmpi_barrier(gwr%comm)
 !stop

 ABI_FREE(sigc_it_diag_kcalc)
 call cwtime_report(" gwr_build_sigmac:", cpu_all, wall_all, gflops_all)

contains
integer function vid(vname)
 character(len=*),intent(in) :: vname
 vid = nctk_idname(ncid, vname)
end function vid

subroutine write_units(units, str)
 character(len=*),intent(in) :: str
 integer,intent(in) :: units(:)
 integer :: ii

 do ii=1,size(units)
   write(units(ii), "(a)") trim(str)
 end do

end subroutine write_units

subroutine write_notations(units)
 integer,intent(in) :: units(:)
 integer :: ii, unt

 do ii=1,size(units)
   unt = units(ii)
   write(unt,"(a)")repeat("=", 80)
   write(unt,"(a)")" QP results in eV."
   write(unt,"(a)")" Notations:"
   !write(unt,"(a)")"     eKS: Kohn-Sham energy. eQP: quasi-particle energy."
   !write(unt,"(a)")"     eQP - eKS: Difference between the QP and the KS energy."
   !write(unt,"(a)")"     SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part."
   !write(unt,"(a)")"     Z(eKS): Renormalization factor."
   !write(unt,"(a)")"     FAN: Real part of the Fan term at eKS. DW: Debye-Waller term."
   !write(unt,"(a)")"     DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP."
   !write(unt,"(a)")"     OTMS: On-the-mass-shell approximation with eQP ~= eKS + Sigma(omega=eKS)"
   !write(unt,"(a)")"     TAU(eKS): Lifetime in femtoseconds computed at the KS energy."
   write(unt,"(a)")" "
   write(unt,"(a)")" "
 end do
end subroutine write_notations

end subroutine gwr_build_sigmac
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/diag_braket
!! NAME
!!  diag_braket
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if 0

subroutine diag_braket(trans, mat, nfftsp, u_glob, cout, do_mpi_sum)

!Arguments ------------------------------------
 character(len=*),intent(in) :: trans
 integer,intent(in) :: nfftsp
 class(matrix_scalapack),intent(in) :: mat
 complex(dp),intent(in) :: u_glob(nfftsp)
 complex(dp),intent(out) :: cout
 logical,optional,intent(in) :: do_mpi_sum

!Local variables-------------------------------
 integer :: ir1, il_r1, ierr
 logical :: do_mpi_sum__
 complex(dp),allocatable :: loc_cwork(:)

! *************************************************************************

 do_mpi_sum__ = .True.; if (present(do_mpi_sum)) do_mpi_sum__ = do_mpi_sum

 select case (trans)
 case ("T")
   ABI_CHECK_IEQ(nfftsp, mat%sizeb_local(1), "First dimension should be local to each MPI proc!")

   ! (r',r) with r' local and r-index PBLAS-distributed.
   ! Integrate over r'
   ABI_MALLOC(loc_cwork, (mat%sizeb_local(2)))
   loc_cwork(:) = matmul(transpose(mat%buffer_cplx), u_glob)

   ! Integrate over r. Note complex conjugate.
   cout = zero
   do il_r1=1, mat%sizeb_local(2)
     ir1 = mat%loc2gcol(il_r1)
     cout = cout + conjg(u_glob(ir1)) * loc_cwork(il_r1)
   end do

   ABI_FREE(loc_cwork)
   if (do_mpi_sum__) call xmpi_sum(cout, mat%processor%grid%ictxt, ierr)

 case ("N")
   ABI_ERROR(sjoin("Not implemented trans:", trans))

 case default
   ABI_ERROR(sjoin("Invalid trans:", trans))
 end select

end subroutine diag_braket
!!***

#endif

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
 integer :: my_is, my_iqi, my_it, itau, spin, iq_ibz, ii
 integer :: il_g1, il_g2, ig1, ig2, npw_q !, ierr, iglob2, iglob1,
 logical :: q_is_gamma
 real(dp) :: weight_q, qq_ibz(3)
 complex(dpc) :: vcs_g1, vcs_g2
!arrays
 type(matrix_scalapack) :: chi_tmp, dummy_vec
 real(dp),allocatable :: eig(:)

! *************************************************************************

 call gwr%build_chi0_head_and_wings()
 call gwr%build_green(free_ugb=.True.)
 call gwr%build_tchi()

 ABI_CHECK(gwr%tchi_space == "iomega", sjoin("tchi_space:", gwr%tchi_space, "!= iomega"))

 ! TODO: Compute RPA energy with different number of PWs and extrapolate.
 ! See also calc_rpa_functional in m_screening_driver
 !spin_fact = two / (gwr%nsppol * gwr%nspinor)
 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     qq_ibz = gwr%qibz(:, iq_ibz)
     weight_q = gwr%wtq(iq_ibz)
     q_is_gamma = normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0

     ! NB: Take into account that iq_ibz might be replicated across procs.
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       associate (desc_q => gwr%tchi_desc_qibz(iq_ibz), tchi => gwr%tchi_qibz(iq_ibz, itau, spin))
       if (my_it == 1) then
         ! Allocate workspace. NB: npw_q is the total number of PWs for this q.
         call tchi%copy(chi_tmp)
         npw_q = desc_q%npw
         ABI_CHECK_IEQ(npw_q, tchi%sizeb_global(1), "npw_q")
         ABI_MALLOC(eig, (npw_q))
       end if

       do il_g2=1,tchi%sizeb_local(2)
         ig2 = mod(tchi%loc2gcol(il_g2) - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,tchi%sizeb_local(1)
           ig1 = mod(tchi%loc2grow(il_g1) - 1, desc_q%npw) + 1
           vcs_g1 = desc_q%vc_sqrt(ig1)
           chi_tmp%buffer_cplx(il_g1, il_g2) = tchi%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2
         end do
       end do

       call chi_tmp%pzheev("N", "U", dummy_vec, eig)

       ! Perform integration in imaginary frequency.
       do ii=1,desc_q%npw
         !eig(ii)
         !gwr%iw_wgs(itau)
       end do

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

 call gwr%build_chi0_head_and_wings()
 call gwr%build_sigxme()
 call gwr%build_green(free_ugb=.True.)
 call gwr%build_tchi()
 call gwr%build_wc()
 call gwr%build_sigmac()

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
 real(dp), ABI_CONTIGUOUS pointer :: fptr(:,:,:)
 type(matrix_scalapack), pointer :: mats(:)

! *************************************************************************

 ! Cannot reuse SCR.nc/SUSC.nc fileformat as:
 !  - hscr_new requires ep%
 !  - these file formats assume Gamma-centered gvectors.

 all_nproc = xmpi_comm_size(gwr%comm); my_rank = xmpi_comm_rank(gwr%comm)
 call cwtime(cpu, wall, gflops, "start")

 if (my_rank == master) then
   call wrtout(std_out, sjoin(" Writing", what, "to:", filepath))
   NCF_CHECK(nctk_open_create(ncid, filepath, xmpi_comm_self))
   NCF_CHECK(gwr%cryst%ncwrite(ncid))

   ! Add dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nsppol", gwr%nsppol), nctkdim_t("ntau", gwr%ntau), nctkdim_t("mpw", gwr%tchi_mpw), &
     nctkdim_t("nqibz", gwr%nqibz), nctkdim_t("nqbz", gwr%nqbz)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define arrays with results.
   ! TODO: Add metadata for mats: spin sum, vc cutoff, t/w mesh, handle nspinor 2
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngkpt", "int", "three"), &
     nctkarr_t("ngqpt", "int", "three"), &
     nctkarr_t("qibz", "dp", "three, nqibz"), &
     nctkarr_t("wtq", "dp", "nqibz"), &
     nctkarr_t("tau_mesh", "dp", "ntau"), &
     nctkarr_t("tau_wgs", "dp", "ntau"), &
     nctkarr_t("iw_mesh", "dp", "ntau"), &
     nctkarr_t("iw_wgs", "dp", "ntau"), &
     nctkarr_t("gvecs", "int", "three, mpw, nqibz"), &
     nctkarr_t("mats", "dp", "two, mpw, mpw, ntau, nqibz, nsppol") &
   ])
   NCF_CHECK(ncerr)

   ! Write global arrays.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("ngkpt"), gwr%ngkpt))
   NCF_CHECK(nf90_put_var(ncid, vid("ngqpt"), gwr%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, vid("qibz"), gwr%qibz))
   NCF_CHECK(nf90_put_var(ncid, vid("wtq"), gwr%wtq))
   NCF_CHECK(nf90_put_var(ncid, vid("tau_mesh"), gwr%tau_mesh))
   NCF_CHECK(nf90_put_var(ncid, vid("tau_wgs"), gwr%tau_wgs))
   NCF_CHECK(nf90_put_var(ncid, vid("iw_mesh"), gwr%iw_mesh))
   NCF_CHECK(nf90_put_var(ncid, vid("iw_wgs"), gwr%iw_wgs))
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
       NCF_CHECK(nf90_put_var(ncid, vid("gvecs"), desc_q%gvec, start=[1,1,iq_ibz], count=[3,npwtot_q,1]))
     end if

     mats => null()
     if (what == "tchi") mats => gwr%tchi_qibz(iq_ibz, :, spin)
     if (what == "wc")   mats => gwr%wc_qibz(iq_ibz, :, spin)
     ABI_CHECK(associated(mats), sjoin("Invalid value for what:", what))

     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       ! NB: Assuming PBLAS matrix distributed in contigous blocks along the column index.
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
 !character(len=500) :: msg

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

     !ABI_CHECK(i1 == modulo(kg_k(1, ipw), n1) + 1, "i1")
     !ABI_CHECK(i2 == modulo(kg_k(2, ipw), n2) + 1, "i2")
     !ABI_CHECK(i3 == modulo(kg_k(3, ipw), n3) + 1, "i3")

     !if (any(kg_k(:,ipw) > ngfft(1:3)/2) .or. any(kg_k(:,ipw) < -(ngfft(1:3)-1)/2) ) then
     !  write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ipw)," falls outside the FFT box. Increase boxcutmin (?)"
     !  ABI_ERROR(msg)
     !end if

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
 !character(len=500) :: msg

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

     !ABI_CHECK(i1 == modulo(kg_k(1, ipw), n1) + 1, "i1")
     !ABI_CHECK(i2 == modulo(kg_k(2, ipw), n2) + 1, "i2")
     !ABI_CHECK(i3 == modulo(kg_k(3, ipw), n3) + 1, "i3")

     !if (any(kg_k(:,ipw) > ngfft(1:3)/2) .or. any(kg_k(:,ipw) < -(ngfft(1:3)-1)/2) ) then
     !  write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ipw)," falls outside the FFT box. Increase boxcutmin (?)"
     !  ABI_ERROR(msg)
     !end if

     ipwdat = ipw + (idat-1) * npw
     i3dat = i3 + (idat-1) * n6

     cg(ipwdat) = cfft(i1,i2,i3dat)
   end do
 end do

end subroutine box2gsph
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/ur_to_scpsi
!! NAME
!!  ur_to_scpsi
!!
!! FUNCTION
!!  Build psi(R) = psi(r+L) in the supercell from the periodic part u(r) given in the unit cell
!!  and the array `sc_eikr` with the phase in the supercell.
!!  Results in sc_psi are stored using the same layout as the one used for the FFT in the supercell
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure subroutine ur_to_scpsi(sc_shape, uc_ngfft, uc_nfft, nspinor, ndat, ur, alpha, sc_eikr, sc_psi)

!Arguments ------------------------------------
 integer,intent(in) :: sc_shape(3)
 integer,intent(in) :: uc_ngfft(18)
 integer,intent(in) :: uc_nfft, nspinor, ndat
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: ur(uc_nfft, nspinor, ndat)
 complex(dp),intent(out) :: &
   sc_psi(uc_ngfft(1)*sc_shape(1), uc_ngfft(2)*sc_shape(2), uc_ngfft(3)*sc_shape(3), nspinor, ndat)
 complex(dp),intent(in) :: &
   sc_eikr(uc_ngfft(1)*sc_shape(1), uc_ngfft(2)*sc_shape(2), uc_ngfft(3)*sc_shape(3), nspinor)

!Local variables-------------------------------
 integer :: ix, iy, iz, spinor, idat, uc_ir, i1, i2, i3

! *************************************************************************

 if (abs(alpha) < tol12) then
!!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i3, i2, i1, uc_ir)
   do idat=1,ndat
     do spinor=1,nspinor
       do iz=1,uc_ngfft(3) * sc_shape(3)
         i3 = mod(iz -1, uc_ngfft(3)) + 1
         do iy=1,uc_ngfft(2) * sc_shape(2)
           i2 = mod(iy -1, uc_ngfft(2)) + 1
           do ix=1,uc_ngfft(1) * sc_shape(1)
             i1 = mod(ix -1, uc_ngfft(1)) + 1
             uc_ir = i1 + (i2 - 1) * uc_ngfft(1) + (i3 - 1) * uc_ngfft(1) * uc_ngfft(2)
             sc_psi(ix, iy, iz, spinor, idat) = &
                 sc_eikr(ix, iy, iz, spinor) * ur(uc_ir, spinor, idat)
           end do
         end do
       end do
     end do
   end do

 else
   ! accumulate
!!$OMP PARALLEL DO COLLAPSE (3) PRIVATE(i3, i2, i1, uc_ir)
   do idat=1,ndat
     do spinor=1,nspinor
       do iz=1,uc_ngfft(3) * sc_shape(3)
         i3 = mod(iz -1, uc_ngfft(3)) + 1
         do iy=1,uc_ngfft(2) * sc_shape(2)
           i2 = mod(iy -1, uc_ngfft(2)) + 1
           do ix=1,uc_ngfft(1) * sc_shape(1)
             i1 = mod(ix -1, uc_ngfft(1)) + 1
             uc_ir = i1 + (i2 - 1) * uc_ngfft(1) + (i3 - 1) * uc_ngfft(1) * uc_ngfft(2)
             sc_psi(ix, iy, iz, spinor, idat) = sc_psi(ix, iy, iz, spinor, idat) + &
                 (alpha * sc_eikr(ix, iy, iz, spinor)) * ur(uc_ir, spinor, idat)
           end do
         end do
       end do
     end do
   end do
 end if

end subroutine ur_to_scpsi
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/calc_sc_ceikr
!! NAME
!! calc_sc_ceikr
!!
!! FUNCTION
!!  Calculate e^{ik.r} on the FFT mesh in the supercell.
!!  starting from the FFT mesh for the unit cell
!!  Results in sc_ceikr are stored using the same layout as the one used for the FFT in the supercell.
!!
!! INPUTS
!!  kk(3): k-point in reduced coordinates of the unit cell
!!  sc_shape: shape of the supercell.
!!  uc_ngfft(18)=: nformation about the 3d FFT for the unit cell.
!!  nspinor=number of spinor components.
!!
!! OUTPUT
!!  sc_ceikr(..., nspinor) = e^{ik.r} on the supercell fft mesh.
!!
!! SOURCE

pure subroutine calc_sc_ceikr(kk, sc_shape, uc_ngfft, nspinor, sc_ceikr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sc_shape(3), nspinor
!arrays
 real(dp),intent(in) :: kk(3)
 integer,intent(in) :: uc_ngfft(18)
 complex(dp),intent(out) :: &
    sc_ceikr(uc_ngfft(1)*sc_shape(1), uc_ngfft(2)*sc_shape(2), uc_ngfft(3)*sc_shape(3), nspinor)

!Local variables-------------------------------
 integer :: ix, iy, iz
 real(dp) :: kdotr

! *************************************************************************

 if (all(abs(kk) < tol12)) then
   sc_ceikr = cone; return
 end if

 do iz=0,uc_ngfft(3) * sc_shape(3) - 1
   do iy=0,uc_ngfft(2) * sc_shape(2) - 1
     do ix=0,uc_ngfft(1) * sc_shape(1) - 1
       kdotr = two_pi*( kk(1) * (ix / dble(uc_ngfft(1))) &
                       +kk(2) * (iy / dble(uc_ngfft(2))) &
                       +kk(3) * (iz / dble(uc_ngfft(3))) )
       sc_ceikr(ix+1, iy+1, iz+1, 1) = cmplx(cos(kdotr), sin(kdotr), kind=dp)
     end do
   end do
 end do

 if (nspinor > 1) then
   do ix=2,nspinor
     sc_ceikr(:,:,:, ix) = sc_ceikr(:,:,:, 1)
   end do
 end if

end subroutine calc_sc_ceikr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/load_head_wings_from_sus_file__
!! NAME
!!  load_head_wings_from_sus_file__
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine load_head_wings_from_sus_file__(gwr, filepath)

!Arguments ------------------------------------
!scalars
 class(gwr_t),intent(inout) :: gwr
 character(len=*),intent(in) :: filepath

!Local variables-------------------------------
!scalars
 integer, parameter :: master = 0
 integer :: my_rank, sus_npw, sus_nomega, ncerr, ncid, nmiss, ig_sus, ig, ierr
 integer :: my_iqi, my_it, my_is, iq_ibz, spin, itau, iw, il_g1, il_g2, ig1, ig2, iglob1, iglob2, ig1_inv, ig2_inv
 real(dp) :: max_err
 logical :: q_is_gamma
!arrays
 !integer :: g1(3), g2(3)
 integer,allocatable :: sus_gvec(:,:), g2sus(:), ginv(:)
 real(dp) :: qq_ibz(3)
 real(dp),ABI_CONTIGUOUS pointer :: fptr2(:,:), fptr4(:,:,:,:)
 complex(dp),target,allocatable :: sus_head(:,:,:), sus_lwing(:,:,:), sus_uwing(:,:,:), sus_freqs(:)

! *************************************************************************

 my_rank = xmpi_comm_rank(gwr%comm)
 call wrtout(std_out, sjoin(" Loading head and wings from:", filepath))

 if (.not. file_exists(filepath)) then
   ABI_WARNING("Cannot find file with head and wings. Results are WRONG. Returning.")
   return
 end if

 ! Read head and wings from SUS file as we are still not able to compute
 ! these quantities in GWR.
 if (my_rank == master) then
   NCF_CHECK(nctk_open_read(ncid, filepath, xmpi_comm_self))

   !NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nctk_get_dim(ncid, "number_of_coefficients_dielectric_function", sus_npw))
   NCF_CHECK(nctk_get_dim(ncid, "number_of_frequencies_dielectric_function", sus_nomega))
   ABI_CHECK_IEQ(sus_nomega, gwr%ntau + 1, sjoin("sus_nomega: ", itoa(sus_nomega), "should be equal to ntau + 1"))
   !print *, "sus_npw, sus_nomega", sus_npw, sus_nomega

   ABI_MALLOC(sus_freqs, (sus_nomega))
   call c_f_pointer(c_loc(sus_freqs), fptr2, shape=[2, sus_nomega])
   NCF_CHECK(nf90_get_var(ncid, vid('frequencies_dielectric_function'), fptr2))
   max_err = maxval(abs(gwr%iw_mesh - aimag(sus_freqs(2:))))
   if (max_err > tol3) then
     do iw=1,gwr%ntau
       write(std_out, *) gwr%iw_mesh(iw), aimag(sus_freqs(iw+1))
     end do
     ABI_ERROR("Mismatch in imaginary frequency mesh")
   end if
   write(std_out, *)"max_err in frequency meshes:", max_err
   ABI_FREE(sus_freqs)

   ! (number_of_reduced_dimensions, number_of_coefficients_dielectric_function, number_of_qpoints_dielectric_function')
   ABI_MALLOC(sus_gvec, (3, sus_npw))
   ncerr = nf90_get_var(ncid, vid("reduced_coordinates_plane_waves_dielectric_function"), sus_gvec, &
                         start=[1, 1, 1], count=[3, sus_npw, 1])
   NCF_CHECK(ncerr)

   ABI_MALLOC(sus_head, (3, 3, sus_nomega))
   call c_f_pointer(c_loc(sus_head), fptr4, shape=[2, 3, 3, sus_nomega])
   NCF_CHECK(nf90_get_var(ncid, vid("sus_head"), fptr4))

   ABI_MALLOC(sus_lwing, (3, sus_npw, sus_nomega))
   call c_f_pointer(c_loc(sus_lwing), fptr4, shape=[2, 3, sus_npw, sus_nomega])
   NCF_CHECK(nf90_get_var(ncid, vid("sus_lower_wing"), fptr4))

   ABI_MALLOC(sus_uwing, (3, sus_npw, sus_nomega))
   call c_f_pointer(c_loc(sus_uwing), fptr4, shape=[2, 3, sus_npw, sus_nomega])
   NCF_CHECK(nf90_get_var(ncid, vid("sus_upper_wing"), fptr4))

   NCF_CHECK(nf90_close(ncid))
 end if

 call xmpi_bcast(sus_npw, master, gwr%comm, ierr)
 call xmpi_bcast(sus_nomega, master, gwr%comm, ierr)

 if (my_rank /= master) then
   ABI_MALLOC(sus_gvec, (3, sus_npw))
   ABI_MALLOC(sus_head, (3, 3, sus_nomega))
   ABI_MALLOC(sus_lwing, (3, sus_npw, sus_nomega))
   ABI_MALLOC(sus_uwing, (3, sus_npw, sus_nomega))
 end if

 call xmpi_bcast(sus_gvec, master, gwr%comm, ierr)
 call xmpi_bcast(sus_head, master, gwr%comm, ierr)
 call xmpi_bcast(sus_lwing, master, gwr%comm, ierr)
 call xmpi_bcast(sus_uwing, master, gwr%comm, ierr)

 !print *, "head"
 !do ig=1,sus_nomega
 !   call print_arr(sus_head(:, :, ig))
 !   call print_arr(matmul(matmul(gwr%cryst%rprimd, sus_head(:,:,ig)), transpose(gwr%cryst%rprimd)) / two_pi**2)
 !   call print_arr(matmul(matmul(transpose(gwr%cryst%rprimd), sus_head(:,:,ig)), gwr%cryst%rprimd) / two_pi**2)
 !   print *, "next"
 !end do
 !stop

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     qq_ibz = gwr%qibz(:, iq_ibz)
     q_is_gamma = normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0
     if (.not. q_is_gamma) cycle

     associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))
     ABI_CALLOC(gwr%chi0_head_myw, (3, 3, gwr%my_ntau))
     ABI_CALLOC(gwr%chi0_uwing_myw, (3, desc_q%npw, gwr%my_ntau))
     ABI_CALLOC(gwr%chi0_lwing_myw, (3, desc_q%npw, gwr%my_ntau))

     ! NB: The imaginary frequency points in the SUS file starts at 2 when AC is used so use itau + 1
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       gwr%chi0_head_myw(:,:,my_it) = sus_head(:,:,itau + 1)
     end do

     ABI_MALLOC(g2sus, (desc_q%npw))
     call kg_map(sus_npw, sus_gvec, desc_q%npw, desc_q%gvec, g2sus, nmiss)
     call wrtout(std_out, sjoin(" nmiss:", itoa(nmiss), "out of:", itoa(desc_q%npw)))
     ABI_CHECK_IEQ(nmiss, 0, "nmiss > 0, increase ecuteps in legacy GW code to have more G-vectors!")

     do ig=1,desc_q%npw
       ig_sus = g2sus(ig)
       do my_it=1,gwr%my_ntau
         itau = gwr%my_itaus(my_it)
         gwr%chi0_uwing_myw(:, ig, my_it) = sus_uwing(:, ig_sus, itau + 1)
         gwr%chi0_lwing_myw(:, ig, my_it) = sus_lwing(:, ig_sus, itau + 1)
       end do
     end do

     ! Now use head and wings to obtain tchi(q) for finite q. See also cchi0q.
     ! Note that along the imaginary axis the matrix is hermitian.
     ! Moreover at q == 0, one has M(-g,-g') = M(g,g')^*
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       associate (chi => gwr%tchi_qibz(iq_ibz, itau, spin))
       do il_g2=1,chi%sizeb_local(2)
         iglob2 = chi%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         do il_g1=1,chi%sizeb_local(1)
           iglob1 = chi%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1

           if (iglob1 == 1 .and. iglob2 == 1) then
             chi%buffer_cplx(il_g1, il_g2) = &
               vdotw(GW_Q0_DEFAULT, matmul(gwr%chi0_head_myw(:,:,my_it), GW_Q0_DEFAULT), gwr%cryst%gmet, "G")
           else if (iglob2 == 1) then
             chi%buffer_cplx(il_g1, il_g2) = &
               vdotw(GW_Q0_DEFAULT, gwr%chi0_lwing_myw(:, ig1, my_it), gwr%cryst%gmet, "G")
               !vdotw(GW_Q0_DEFAULT, conjg(gwr%chi0_uwing_myw(:, ig1, my_it)), gwr%cryst%gmet, "G")
           else if (iglob1 == 1) then
             chi%buffer_cplx(il_g1, il_g2) = &
               vdotw(GW_Q0_DEFAULT, gwr%chi0_uwing_myw(:, ig2, my_it), gwr%cryst%gmet, "G")
           end if

         end do ! il_ig1
       end do ! il_ig2

       call desc_q%calc_ginv_map(ginv)
       ierr = 0
       do ig2=1,chi%sizeb_global(2)
         !g2 = desc_q%gvec(:, ig2)
         ig2_inv = ginv(ig2)
         do ig1=1,chi%sizeb_global(2)
           !g1 = desc_q%gvec(:, ig1)
           ig1_inv = ginv(ig1)
           !print *, ig1_inv, ig2_inv
           !write(std_out, *) real(chi%buffer_cplx(ig1, ig2)), real(chi%buffer_cplx(ig1_inv, ig2_inv))

           if (abs(aimag(chi%buffer_cplx(ig1, ig2)) + aimag(chi%buffer_cplx(ig1_inv, ig2_inv))) > tol6) then
             ierr = ierr + 1
#if 0
             write(std_out, *) aimag(chi%buffer_cplx(ig1, ig2)), aimag(chi%buffer_cplx(ig1_inv, ig2_inv))
             write(std_out, *) desc_q%gvec(:, ig1), desc_q%gvec(:, ig2)
#endif
           end if
         end do
       end do
       ABI_FREE(ginv)
       write(std_out, *)" Found: ", ierr, "pairs over: ", desc_q%npw ** 2, "that do not fulfill g --> -g symmetry"
       !stop
       end associate

     end do ! my_it
     end associate
     ABI_FREE(g2sus)
   end do ! my_iqi
 end do ! my_is

 ABI_FREE(sus_gvec)
 ABI_FREE(sus_head)
 ABI_FREE(sus_lwing)
 ABI_FREE(sus_uwing)

contains
integer function vid(vname)
  character(len=*),intent(in) :: vname
  vid = nctk_idname(ncid, vname)
end function vid

end subroutine load_head_wings_from_sus_file__
!!***

!subroutine r2c_fft(cplex, fofg, fofr, isign, mpi_enreg, nfft, ndat, ngfft, tim_fourdp)
!
!!Arguments ------------------------------------
!!scalars
! integer,intent(in) :: cplex,isign,nfft,ndat,tim_fourdp
! type(MPI_type),intent(in) :: mpi_enreg
!!arrays
! integer,intent(in) :: ngfft(18)
! real(dp),intent(inout) :: fofg(2,nfft,ndat),fofr(cplex*nfft,ndat)
!
!!Local variables-------------------------------
!
! !call c_f_pointer(c_loc(cg_mat%buffer_cplx), fofg, shape=[2, nfft, ndat])
! !call c_f_pointer(c_loc(cg_mat%buffer_cplx), fofr, shape=[cplex*nfft * ndat])
!
! !call fourdp(cplex, fofg, fofr, isign, mpi_enreg, nfft, ndat, ngfft)
!
!end subroutine r2c_fft

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_chi0_head_and_wings
!! NAME
!!  gwr_build_chi0_head_and_wings
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_chi0_head_and_wings(gwr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr

!Local variables-------------------------------
!scalars
!arrays

! *************************************************************************

 call wrtout(std_out, "TODO: Building head and wings")
 call wrtout(std_out, sjoin("nkbz:", itoa(gwr%nkbz)))

end subroutine gwr_build_chi0_head_and_wings
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_sigxme
!! NAME
!!  gwr_build_sigxme
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_sigxme(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer :: nsppol, nspinor, ierr, my_ikf, ib_sum, ii, iab, kb, il_b, ig
 integer :: my_is, ikcalc, ikcalc_ibz, bmin, bmax
 integer :: spin, isym, jb, is_idx
 integer :: spad, spadx1, spadx2, irow, npw_k, ndegs, wtqm, wtqp
 integer :: npwx, x_nfft, x_mgfft, x_fftalga, nsig_ab
 integer :: ik_bz, ik_ibz, isym_k, trev_k, g0_k(3), tsign_k
 integer :: iq_bz, iq_ibz, isym_q, trev_q, g0_q(3), tsign_q
 logical :: isirr_k, isirr_q
 real(dp) :: cpu, wall, gflops, fact_spin, theta_mu_minus_esum, theta_mu_minus_esum2, tol_empty, tol_empty_in
 character(len=5000) :: msg
 type(matrix_scalapack),pointer :: ugb_kibz
 type(crystal_t),pointer :: cryst
 type(littlegroup_t) :: ltg_k
!arrays
 integer :: g0(3), spinor_padx(2,4)
 integer :: x_ngfft(18)
 real(dp) :: ksum(3), kgw(3), kgw_m_ksum(3), qbz(3), q0(3), spinrot_kbz(4), spinrot_kgw(4) !, tsec(2)
 real(dp),contiguous, pointer :: qp_ene(:,:,:), qp_occ(:,:,:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:), rhotwg(:), rhotwgp(:), rhotwg_ki(:,:), ur_bdgw(:,:), ur_ibz(:), ur_ksum(:)
 complex(dpc),allocatable  :: sigxcme_tmp(:,:), sigxme_tmp(:,:,:), sym_sigx(:,:,:), sigx(:,:,:,:)

! *************************************************************************

 call wrtout(std_out, "TODO: Computing matrix elements of Sigma_x ...")

 nsppol = gwr%nsppol; nspinor = gwr%nspinor
 cryst => gwr%cryst

 if (gwr%sc_iteration == 0) then
   qp_ene => gwr%ks_ebands%eig; qp_occ => gwr%ks_ebands%occ
 else
   qp_ene => gwr%qp_ebands%eig; qp_occ => gwr%qp_ebands%occ
 end if

 ! MRM allow lower occ numbers
 ! Normalization of theta_mu_minus_esum. If nsppol==2, qp_occ $\in [0,1]$
 tol_empty_in = 0.01                            ! Initialize the tolerance used to decide if a band is empty (passed to m_sigx.F90)
 select case (nsppol)
 case (1)
   fact_spin = half; tol_empty = tol_empty_in          ! below this value the state is assumed empty
   if (nspinor == 2) then
     fact_spin = one; tol_empty = half * tol_empty_in  ! below this value the state is assumed empty
   end if
 case (2)
   fact_spin = one; tol_empty = half * tol_empty_in  ! to be consistent and obtain similar results if a metallic
 case default                                        ! spin unpolarized system is treated using nsppol==2
   ABI_BUG(sjoin('Wrong nsppol:', itoa(nsppol)))
 end select

 ! FFT mesh from ecutsigx
 ! =========================================
 ! Find FFT mesh and max number of g-vectors
 ! =========================================
 ! Note the usage of ecutsigx and gwr_boxcutmin and the loops over the full BZ.
 ! All the procs execute this part.

 !x_ngfft = gwr%dtset%ngfft ! Allow user to specify fftalg
 !x_ngfft(1:6) = 0

 !x_mpw = -1
 !do ik_bz=1,gwr%nkbz
 !  kk_bz = gwr%kbz(:, ik_bz)
 !  call get_kg(kk_bz, istwfk1, dtset%ecut, cryst%gmet, npw_, gvec_)
 !  ABI_FREE(gvec_)
 !  call getng(dtset%gwr_boxcutmin, dtset%chksymtnons, dtset%ecut, cryst%gmet, &
 !             kk_bz, me_fft0, x_mgfft, gwr%g_nfft, x_ngfft, nproc_fft1, cryst%nsym, paral_fft0, &
 !             cryst%symrel, cryst%tnons, unit=dev_null)
 !  x_mpw = max(x_mpw, npw_)
 !end do

 x_nfft = product(x_ngfft(1:3)); x_mgfft = maxval(x_ngfft(1:3)); x_fftalga = x_ngfft(7) / 100
 npwx = 1
 nsig_ab = nspinor ** 2; spinor_padx = reshape([0, 0, npwx, npwx, 0, npwx, npwx, 0], [2, 4])

 do my_is=1,gwr%my_nspins
 spin = gwr%my_spins(my_is)
 do ikcalc=1,gwr%nkcalc ! TODO: Should be spin dependent!

   ikcalc_ibz = gwr%kcalc2ibz(ikcalc, 1)
   kgw = gwr%kcalc(:, ikcalc); bmin =  gwr%bstart_ks(ikcalc, spin); bmax = gwr%bstop_ks(ikcalc, spin)

   ! ==============================================================
   ! ==== Find little group of the k-points for GW corrections ====
   ! ==============================================================
   ! * The little group is used only if symsigma == 1
   ! * If use_umklp == 1 then symmetries requiring an umklapp to preserve k_gw are included as well.
   !call ltg_k%init(kgw, Qmesh, cryst, use_umklp=1, npwe=0)

   write(msg,'(6a)') ch10, &
    ' Calculating <nk|Sigma_x|nk> at k: ',trim(ktoa(kgw)), ", for bands: ", trim(ltoa([bmin, bmax])),ch10
   call wrtout(std_out, msg)

   ! ===============================================
   ! Load wavefunctions for Sigma_x matrix elements
   ! ===============================================
   ! All procs need ur_bdgw but the IBZ is distributed and, possibly, replicated in gwr%kpt_comm.
   ! Here we select the right procs, fill the buffer with the FFT results and then use a dumb xmpi_sum
   ! to gather the results.
   ABI_MALLOC_OR_DIE(ur_bdgw, (x_nfft * nspinor, bmin:bmax), ierr)
   ur_bdgw = zero

   !ABI_MALLOC(gbound_kcalc, (2 * x_mgfft + 8, 2))
   !call sphereboundary(gbound_kcalc, desc_kbz%istwfk, desc_kbz%gvec, x_mgfft, desc_kbz%npw)

   !call fft_ug(npw_kcalc, x_nfft, nspinor, ndat1, &
   !            x_mgfft, x_ngfft, istwf_kcalc, desc_k%gvec, desc_k%gbound, &
   !            gwr%ugb(ikcalc_ibz, spin)%buffer_cplx(:, il_b), &
   !            ur_bdgw)    ! out
   !ABI_FREE(gbound_kcalc)

   !ABI_MALLOC(ur_ibz, (x_nfft * nspinor))
   ABI_MALLOC(ur_ksum, (x_nfft * nspinor))
   ABI_MALLOC(rhotwg_ki, (npwx * nspinor, bmin:bmax))
   ABI_MALLOC(rhotwg, (npwx * nspinor))
   ABI_MALLOC(rhotwgp, (npwx * nspinor))
   ABI_MALLOC(vc_sqrt_qbz, (npwx))

   ABI_CALLOC(sigxme_tmp, (bmin:bmax, bmin:bmax, nsppol * nsig_ab))
   ABI_CALLOC(sigxcme_tmp, (bmin:bmax, nsppol * nsig_ab))
   ABI_CALLOC(sigx, (2, bmin:bmax, bmin:bmax, nsppol * nsig_ab))

   ! ==============================
   ! ==== Sum over k in the BZ ====
   ! ==============================

   do my_ikf=1,gwr%my_nkbz
     ik_bz = gwr%my_kbz_inds(my_ikf)
     ksum = gwr%kbz(:, ik_bz)

     ! Find the symmetrical image of ksum in the IBZ
     !call kmesh%get_BZ_item(ik_bz, ksum, ik_ibz, isym_ki, iik, ph_mkt)

     ! FIXME: Be careful with the symmetry conventions here!
     ik_ibz = gwr%kbz2ibz(1, ik_bz); isym_k = gwr%kbz2ibz(2, ik_bz)
     trev_k = gwr%kbz2ibz(6, ik_bz); g0_k = gwr%kbz2ibz(3:5, ik_bz)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     !tsign_k = merge(1, -1, trev_k == 0)

     ! Identify q and G0 where q + G0 = k_GW - ksum
     kgw_m_ksum = kgw - ksum
     !call findqg0(iq_bz, g0, kgw_m_ksum, Qmesh%nbz, Qmesh%bz, Sigp%mG0)

     ! If symmetries are exploited only q-points in the IBZ_k are computed.
     ! In this case elements are weighted according to wtqp and wtqm.
     ! wtqm is for time-reversal.
     wtqp = 1; wtqm = 0
     !if (can_symmetrize(spin)) then
     !  if (ltg_k%ibzq(iq_bz) /= 1) cycle
     !  wtqp = 0; wtqm = 0
     !  do isym=1,ltg_k%nsym_sg
     !    wtqp = wtqp + ltg_k%wtksym(1, isym, iq_bz)
     !    wtqm = wtqm + ltg_k%wtksym(2, isym, iq_bz)
     !  end do
     !end if

     ! Find the corresponding irreducible q-point.
     ! NB: non-zero umklapp G_o is not allowed. There's a check in setup_sigma
     !call qmesh%get_BZ_item(iq_bz, qbz, iq_ibz, isym_q, itim_q)
     !q_is_gamma = normv(qbz, cryst%gmet, "G") < GW_TOLQ0

     !qq_bz = gwr%qbz(:, iq_bz)
     !iq_ibz = gwr%qbz2ibz(1, iq_bz); isym_q = gwr%qbz2ibz(2, iq_bz)
     !trev_q = gwr%qbz2ibz(6, iq_bz); g0_q = gwr%qbz2ibz(3:5, iq_bz)
     !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     !tsign_q = merge(1, -1, trev_q == 0)

     ! Get Fourier components of the Coulomb interaction in the BZ
     ! In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     do ig=1,npwx
       !ig_rot = Gsph_x%rottb(ig, itim_q, isym_q)
       !vc_sqrt_qbz(ig_rot) = Vcp%vc_sqrt_resid(ig, iq_ibz)
     end do

     ! ==========================
     ! Sum over (occupied) bands
     ! ==========================
     ugb_kibz => gwr%ugb(ik_ibz, spin)

     do il_b=1,ugb_kibz%sizeb_local(2)
       ib_sum = ugb_kibz%loc2gcol(il_b)

       ! Skip empty states. MRM: allow negative occ numbers.
       if (abs(qp_occ(ib_sum, ik_ibz, spin)) < tol_empty) CYCLE

       !call wfd%get_ur(ib_sum, ik_ibz, spin, ur_ibz)

       ! Compute ur_ksum(r) from the symmetrical image.
       !ugb_kibz%buffer_cplx(:, il_b)
       !call fft_ug(npw_k, x_nfft, nspinor, ndat1, &
       !            x_mgfft, x_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
       !            ggp%buffer_cplx(:, ig2), &  ! in
       !            ur_ksum)    ! out

       ! Get all <k-q,ib_sum,s|e^{-i(q+G).r}|s,jb,k>
       do jb=bmin,bmax

         ! Multiply by the square root of the Coulomb term
         ! In 3-D systems, the factor sqrt(4pi) is included
         do ii=1,nspinor
           spad = (ii-1) * npwx
           !rhotwg_ki(spad+1:spad+npwx,jb) = rhotwg_ki(spad+1:spad + npwx,jb) * vc_sqrt_qbz(1:npwx)
         end do

#if 0
         if (ik_bz == jk_bz) then
           ! Treat analytically the case q --> 0:
           !
           !   * The oscillator is evaluated at q = 0 as it is considered constant in the small cube around Gamma
           !     while the Colulomb term is integrated out.
           !   * If nspinor == 1, we have nonzero contribution only if ib_sum == jb
           !   * If nspinor == 2, we evaluate <ib_sum,up|jb,up> and <ib_sum,dwn|jb,dwn>,
           !     and impose orthonormalization since npwwfn might be < npwvec.
           !   * Note the use of i_sz_resid and not i_sz, to account for the possibility
           !     to have generalized KS basis set from hybrid

           if (nspinor == 1) then
             rhotwg_ki(1, jb) = czero_gw
             if (ib_sum == jb) rhotwg_ki(1,jb) = CMPLX(SQRT(Vcp%i_sz_resid), 0.0_gwp)
             !rhotwg_ki(1,jb) = czero_gw ! DEBUG

           else
             npw_k = wfd%npwarr(ik_ibz)
             rhotwg_ki(1, jb) = zero; rhotwg_ki(npwx+1, jb) = zero
             if (ib_sum == jb) then
               ABI_CHECK(wfd%get_wave_ptr(ib_sum, ik_ibz, spin, wave_sum, msg) == 0, msg)
               cg_sum => wave_sum%ug
               ABI_CHECK(wfd%get_wave_ptr(jb, jk_ibz, spin, wave_jb, msg) == 0, msg)
               cg_jb  => wave_jb%ug
               ctmp = xdotc(npw_k, cg_sum(1:), 1, cg_jb(1:), 1)
               rhotwg_ki(1, jb) = CMPLX(SQRT(Vcp%i_sz_resid), 0.0_gwp) * real(ctmp)
               ctmp = xdotc(npw_k, cg_sum(npw_k+1:), 1, cg_jb(npw_k+1:), 1)
               rhotwg_ki(npwx+1, jb) = CMPLX(SQRT(Vcp%i_sz_resid), 0.0_gwp) * real(ctmp)
             end if
             !rhotwg_ki(1, jb) = zero; rhotwg_ki(npwx+1, jb) = zero
             ! PAW is missing
           end if
         end if
#endif
       end do ! jb Got all matrix elements from bmin up to bmax.

       theta_mu_minus_esum  = fact_spin * qp_occ(ib_sum, ik_ibz, spin)
       theta_mu_minus_esum2 = sqrt(abs(fact_spin * qp_occ(ib_sum, ik_ibz, spin))) ! MBB Nat. orb. funct. approx. sqrt(occ)


       if (abs(theta_mu_minus_esum / fact_spin) >= tol_empty) then     ! MRM: allow negative occ numbers
         do kb=bmin,bmax
#if 0
           ! Copy the ket Sigma_x |phi_{k,kb}>.
           rhotwgp(:) = rhotwg_ki(:, kb)

           ! Loop over the non-zero row elements of this column.
           ! If gwcalctyp <  20: only diagonal elements since QP == KS.
           ! If gwcalctyp >= 20:
           !      * Only off-diagonal elements connecting states with same character.
           !      * Only the upper triangle if HF, SEX, or COHSEX.

           do irow=1,Sigxij_tab(spin)%col(kb)%size1
             jb = Sigxij_tab(spin)%col(kb)%bidx(irow)
             rhotwg(:) = rhotwg_ki(:,jb)

             ! Calculate bare exchange <phi_jb|Sigma_x|phi_kb>.
             ! Do the scalar product only if ib_sum is occupied.
             do iab=1,nsig_ab
               spadx1 = spinor_padx(1, iab); spadx2 = spinor_padx(2, iab)
               xdot_tmp = -XDOTC(npwx, rhotwg(spadx1+1:), 1, rhotwgp(spadx2+1:), 1)
               gwpc_sigxme  = xdot_tmp * theta_mu_minus_esum
               gwpc_sigxme2 = xdot_tmp * theta_mu_minus_esum2

               ! Accumulate and symmetrize Sigma_x matrix elements.
               ! -wtqm comes from time-reversal (exchange of band indeces)
               is_idx = spin; if (nspinor == 2) is_idx = iab
               sigxme_tmp(jb, kb, is_idx) = sigxme_tmp(jb, kb, is_idx) + &
                  (wtqp + wtqm) * DBLE(gwpc_sigxme) + (wtqp - wtqm) * j_gw * AIMAG(gwpc_sigxme)
               if (jb == kb) then
                 sigxcme_tmp(jb, is_idx) = sigxcme_tmp(jb, is_idx) + &
                   (wtqp + wtqm) * DBLE(gwpc_sigxme2) + (wtqp - wtqm) *j_gw * AIMAG(gwpc_sigxme2)
               end if

               sigx(1, jb, kb, is_idx) = sigx(1, jb, kb, is_idx) + wtqp *      gwpc_sigxme
               sigx(2, jb, kb, is_idx) = sigx(2, jb, kb, is_idx) + wtqm *CONJG(gwpc_sigxme)
             end do
           end do ! jb
#endif
         end do ! kb
       end if

     end do ! ib_sum
   end do ! my_ifk Got all diagonal (off-diagonal) matrix elements.

   ABI_FREE(ur_bdgw)
   !ABI_FREE(ur_ibz)
   ABI_FREE(ur_ksum)
   ABI_FREE(rhotwg_ki)
   ABI_FREE(rhotwg)
   ABI_FREE(rhotwgp)
   ABI_FREE(vc_sqrt_qbz)
   ABI_FREE(sigxme_tmp)
   ABI_FREE(sigxcme_tmp)
   ABI_FREE(sigx)

   call ltg_k%free()
 end do ! ikcalc
 end do ! my_is

end subroutine gwr_build_sigxme
!!***

!----------------------------------------------------------------------

end module m_gwr
!!***

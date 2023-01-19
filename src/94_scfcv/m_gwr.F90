!!****m* ABINIT/m_gwr
!! NAME
!!  m_gwr
!!
!! FUNCTION
!!  Objects and procedures implementing GW in real-space and imaginary time.
!!
!! NOTES
!!   Memory and workload are distributed using a 4D cartesian grid: (g/r, tau, k-points, spin).
!!
!!   Inside the g/r communicator, we use PBLAS matrices to store G, tchi and W.
!!   using a 1D processor grid with block distribution along columns.
!!   A 2D grid, indeed, would require MPI-FFT or some communication before performing the FFTs.
!!
!!   Let's assume for simplicity that we have only two MPI procs in the g/r communicator.
!!   Matrices in (g,g') space are distributed along columns so that the g-index is local
!!   and we can use sequential zero-padded FFTs to transform from g to r in the unit cell:
!!
!!                       g'-axis
!!                |--------------------
!!                |         |         |
!!       g-axis   |   P0    |   P1    |
!!                |         |         |
!!                |--------------------
!!
!!   The results of the FFT transform along g are stored in another PBLAS matrix with the same layout:
!!
!!                       g'-axis
!!                |--------------------
!!                |         |         |
!!       r-axis   |   P0    |   P1    |
!!                |         |         |
!!                |--------------------
!!
!!   At this point, we use ptrans to MPI transpose the (r, g') matrix, and we end up with:
!!
!!                       r-axis
!!                |--------------------
!!                |         |         |
!!       g'-axis  |   P0    |   P1    |
!!                |         |         |
!!                |--------------------
!!
!!   Differences with respect to the GW code in frequency-domain:
!!
!!    - in GWR, the k-mesh must be Gamma-centered.
!!    - All the two-point functions are defined on k/q-centered g-spheres while GW uses a single Gamma-centered sphere.
!!    - The frequency/tau meshes are automatically defined by ntau and the KS spectrum (minimax meshes)
!!
!!   Technical properties:
!!
!!     - it's not clear to me that one can use vc(Sq, SG) when a cutoff is used as the cutoff breaks
!!       the spherical symmetry of vc(r). Besides, when symmetries are used to reconstruct the term for q in the BZ,
!!       one might have to take into account umklapps. Use cache?
!!
!!     - Treatment of the anisotropic behaviour of Wc. This part is badly coded in GW, in the sense that
!!       we use a finite small q when computing Wc for q --> 0. This breaks the symmetry of the system
!!       and QP degeneracies. The equations needed to express the angular dependency of W(q) for q --> 0
!!       are well known but one has to pass through the Adler-Wiser expression.
!!       Solution: Compute heads and wings using a WFK_fine wavefunction file with dense k-mesh and less bands.
!!       The dipole matrix elements are computed with the DFPT routines, still we need to
!!       recode a lot of stuff that is already done in cchi0q0, especially symmetries.
!!       Note, however, that tchi is Hermitian along the imaginary axis, expect for omega = 0 in metals
!!       but I don't think the minmax grids contain omega = 0.
!!
!!    - In principle, it's possible to compute QP correction along a k-path if a new WFK file is provided.
!!      The correlated part is evaluated in real-space in the super-cell.
!!      For Sigma_x, we need a specialized routine that can handle arbitrary q, especially at the level of v(q, G)
!!      but I don't know if this approach will give smooth bands
!!      as we don't have q --> 0 when k does not belong to the k-mesh.
!!
!!    - New routine to compute oscillator matrix elements with NC/PAW and PBLAS matrices.
!!      It can be used to compute tchi head/wings as well as Sigma_x + interface with coupled-cluster codes.
!!
!!    - Decide whether we should use VASP conventions for G and the analytic continuation or the "standard" ones by Godby.
!!      The standard ones are consistent with Hedin's notations and correspond to the ones used in the legacy GW code.
!!      On the other hand, VASP notations make life easier if one has to implement PAW.
!!
!!    - Address nspinor = 2 and PBLAS distribution as MPI proc can have both spinors in memory
!!      In other words, we should store the first/last index in gvec for each spinor
!!
!!    - Optimization for Gamma-only. Memory and c -> r FFTs
!!
!!    - Need to extend FFT API to avoid scaling if isign = -1. Also fft_ug and fft_ur should accept isign
!!      optional argument. Refactoring of all the FFT routines used in the GW code is needed
!!      in order to exploit R2C, C2R (e.g. chi0(q=0) and GPU version.
!!
!!    - TODO: Possible incompatibilities between gwpc, slk matrices that are always in dp and GW machinery
!!      Single precision for scalapack matrices?
!!
!!    - Use round-robin distribution instead of blocked-distribution to improve load balance.
!!
!!    - Memory peaks:
!!
!!        (env3.9) [magianto@uan01 /scratch/project_465000061/magianto/DDIAGO_ZnO]
!!        $~/git_repos/abinit/tests/Scripts/abimem.py peaks abimem_rank0.mocc
!!        [0] <var=gt_scbox, A@m_gwr.F90:3395, addr=0x14aa53673010, size_mb=379.688>
!!        [1] <var=xsum, A@xmpi_sum.finc:2551, addr=0x14aa2fce9010, size_mb=379.688>
!!        [2] <var=gt_scbox, A@m_gwr.F90:4338, addr=0x14aa4f64f010, size_mb=379.688>
!!        [3] <var=allcg_k, A@m_wfd.F90:4631, addr=0x14aa56b57010, size_mb=217.865>
!!        [4] <var=chit_scbox, A@m_gwr.F90:3396, addr=0x14aa4789a010, size_mb=189.844>
!!        [5] <var=wct_scbox, A@m_gwr.F90:4339, addr=0x14aa43876010, size_mb=189.844>
!!        [6] <var=xsum, A@xmpi_sum.finc:2476, addr=0x14aa31bb0010, size_mb=189.844>
!!        [7] <var=cg_k, A@m_wfd.F90:4623, addr=0x14aa64535010, size_mb=108.932>
!!
!!  TODO
!!  - Remove cryst%timrev, use kptopt and qptopt
!!  - Sig_c breaks QP degeneracies due to q0.
!!
!! NOTES:
!!
!!  1) _slk_mat_t is a macro defined in abi_common.h that allows us to use PBLAS in single/double precision
!!     Be careful when using c_f_pointer because there's no type checking.
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
 use m_gwdefs,        only : GW_TOL_DOCC, GW_TOLQ0, GW_TOL_W0, GW_Q0_DEFAULT, cone_gw, czero_gw, sigijtab_t, &
                             sigijtab_free, g0g0w
 use m_time,          only : cwtime, cwtime_report, sec2str, timab
 use m_io_tools,      only : iomode_from_fname, get_unit, file_exists, open_file, write_units
 use m_pstat,         only : pstat_t
 use m_numeric_tools, only : blocked_loop, get_diag, isdiagmat, arth, print_arr, imin_loc, imax_loc, &
                             c2r, linfit, bisect, hermitianize
 use m_copy,          only : alloc_copy
 use m_geometry,      only : normv, vdotw
 use m_fstrings,      only : sjoin, itoa, strcat, ktoa, ltoa, ftoa
 use m_sort,          only : sort_dp, sort_rvals
 use m_krank,         only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_crystal,       only : crystal_t
 use m_dtset,         only : dataset_type
 use m_fftcore,       only : get_kg, sphereboundary, getng, print_ngfft, fftcore_set_mixprec, ngfft_seq
 use m_cgtk,          only : cgtk_rotate
 use m_mpinfo,        only : initmpi_seq, destroy_mpi_enreg
 use m_distribfft,    only : init_distribfft_seq
 use m_fft,           only : fftbox_plan3_t, uplan_t, fft_ug, fft_ur !, fourdp
 use m_fft_mesh,      only : calc_ceikr, calc_ceigr, ctimes_eikr
 use m_kpts,          only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_map_print, kpts_pack_in_stars
 use m_bz_mesh,       only : littlegroup_t, findqg0
 use m_gsphere,       only : kg_map, gsphere_t
 use m_melemts,       only : melements_t
 use m_ioarr,         only : read_rhor
 use m_slk,           only : matrix_scalapack, slkmat_sp_t, processor_scalapack, slk_array_free, slk_array_set, &
                             slk_array_locmem_mb, block_dist_1d, slk_pgemm
 use m_wfk,           only : wfk_read_ebands, wfk_t, wfk_open_read
 use m_wfd,           only : wfd_init, wfd_t, wfdgw_t
 use m_ddk,           only : ddkop_t, ddkop_new
 use m_pawtab,        only : pawtab_type
 use m_pawcprj,       only : pawcprj_type
 use m_vcoul,         only : vcgen_t
 use m_vkbr,          only : vkbr_t, vkbr_free, vkbr_init, nc_ihr_comm
 use m_chi0tk,        only : chi0_bbp_mask, accumulate_head_wings_imagw, symmetrize_afm_chi0
 use m_sigx,          only : sigx_symmetrize
 use m_dyson_solver,  only : sigma_pade_t
#ifdef __HAVE_GREENX
 use gx_minimax,      only : gx_minimax_grid, gx_get_error_message
#endif

 implicit none

 private
!!***


!!****t* m_gwr/desc_t
!! NAME
!! desc_t
!!
!! FUNCTION
!!  Parameters related to a two-point function such as
!!  gvectors, tables used for zero padded FFTs and matrix elements of the Coulomb interaction.
!!
!! SOURCE

 type,public :: desc_t

   integer :: istwfk = 1
   ! Storage mode for this k/q point.

   integer :: npw = -1
   ! Total number of plane-waves for this k/q-point.

   integer :: ig0 = -1
   ! Index of g=0 in gvec.

   logical :: kin_sorted
   ! True if gvec is sorted by |q+g|^2/2

   integer,allocatable :: gvec(:,:)
   ! (3, npw)
   ! G-vectors in reduced coordinates.
   ! Note that this array is global i.e. it is not MPI-distributed inside the PBLAS communicator.

   integer,allocatable :: gbound(:,:)
   ! (2*mgfft+8, 2)
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

   procedure :: get_vc_sqrt => desc_get_vc_sqrt
   ! Compute square root of vc(q,g).

   procedure :: free => desc_free
   ! Free memory.

 end type desc_t

 interface desc_array_free
   module procedure desc_array1_free
 end interface desc_array_free
!!***

!----------------------------------------------------------------------

!!****t* m_gwr/est_t
!! NAME
!! est_t
!!
!! FUNCTION
!! Memory is given in Mb
!!
!! SOURCE

 type, public :: est_t

   real(dp) :: mem_green_gg = zero
   real(dp) :: mem_green_rg = zero
   real(dp) :: mem_chi_gg = zero
   real(dp) :: mem_chi_rg = zero
   real(dp) :: mem_ugb = zero
   real(dp) :: mem_total = zero
   real(dp) :: efficiency = zero
   real(dp) :: speedup = zero

 contains
   procedure :: print => est_print
 end type est_t
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

  integer :: nkcalc
   ! Number of Sigma_nk k-points computed
   ! TODO: Should be spin dependent + max_nkcalc

  integer :: max_nbcalc
   ! Maximum number of bands computed (max over nkcalc and spin).

  integer :: nwr = -1
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Odd number so that the mesh is centered on the KS energy.

  !real(dp) :: i_sz = huge(one)
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  real(dp) :: wr_step = -one
   ! Step of the linear mesh along the real axis (Ha units).

  real(dp) :: q0(3) = GW_Q0_DEFAULT
   ! The small q for the treatment of q --> 0

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.

  logical :: idle_proc = .False.
  ! True if there are idle procs i.e. if processes in the input_comm have been excluded.

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
   ! Number of processors in kpt_comm treating ik_ibz

   logical, allocatable :: itreat_ikibz(:)
   ! (nkibz)

   logical, allocatable :: itreat_iqibz(:)
   ! (nqibz)

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
   ! Max shifts to account for umklapps.

   type(desc_t),allocatable :: green_desc_kibz(:)
   ! (nkibz)
   ! Descriptor for Green's functions

   type(desc_t),allocatable :: tchi_desc_qibz(:)
   ! (nqibz)
   ! Descriptor for tchi. NB: The g-vectors are sorted by |q+g|^2/2

   !type(desc_t),allocatable :: sigma_desc_kibz(:)
   ! (nkibz)
   ! Descriptor for self-energy

   integer :: coords_gtks(4) = 0
   ! Cartesian coordinates of this processor in the Cartesian grid.

   type(xcomm_t) :: comm
   ! Communicator with all MPI procs involved in the computation
   ! NB: gwr%comm%value is not necessarly the same as the input_comm
   ! we may decide to remove some procs from input_comm before createring the Cartesian grid.

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

   type(xcomm_t) :: kg_comm
   ! MPI communicator for g/g 2D subgrid.

   type(xcomm_t) :: kts_comm
   ! MPI communicator for tau/kpoint/spin 3D grid

   type(xcomm_t) :: kgt_comm
   ! MPI communicator for g/tau/kpoint 3D grid

   type(dataset_type), pointer :: dtset => null()
   ! Input variables.

   type(datafiles_type), pointer :: dtfil => null()
   ! Names of input/output files and prefixes.

   type(crystal_t), pointer :: cryst => null()
   ! Crystal structure.

   integer :: scf_iteration = 1
   ! Internal counter used to implement self-consistency
   ! For the time being, only self-consistency in energies is supported.

   integer,allocatable :: ks_vbik(:,:)
   ! (gwr%ks_ebands%nkpt, gwr%ks_ebands%nsppol)
   ! KS valence band indices.

   type(ebands_t), pointer :: ks_ebands => null()
   ! initial KS energies

   type(ebands_t) :: qp_ebands
   ! QP energies

   type(ebands_t) :: qp_ebands_prev
   ! QP energies of the previous iteration. Used if self-consistency.

   type(pseudopotential_type), pointer :: psps => null()
   ! NC Pseudos data

   type(pawtab_type), pointer :: pawtab(:) => null()
   ! PAW data

   type(mpi_type),pointer :: mpi_enreg => null()
   ! Sequential mpi_type needed to invoke routines requiring it.

   type(processor_scalapack) :: g_slkproc

   type(__slkmat_t),allocatable :: gt_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)
   ! Occupied/Empty Green's function G_k(g,g')

   type(__slkmat_t),allocatable :: tchi_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! Irreducible polarizability tchi_q(g,g')

   character(len=10) :: tchi_space = "none"
   ! "none", "itau", "iomega"

   type(__slkmat_t),allocatable :: wc_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! Correlated screened Coulomb interaction summed over collinear spins
   ! Replicated across spin_comm if nsppol == 2.

   character(len=10) :: wc_space = "none"
   ! "none", "itau", "iomega"

   !type(__slkmat_t),allocatable :: em1_qibz(:,:,:)
   ! Inverse dielectric matrix at omega = 0
   ! (nqibz, nsppol)
   ! Replicated across the tau comm and the spin comm if nsppol == 2.

   type(__slkmat_t),allocatable :: sigc_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)

   character(len=10) :: sigc_space = "none"
   ! "none", "itau", "iomega"

   type(__slkmat_t),allocatable :: ugb(:,:)
   ! (nkibz, nsppol)
   ! Fourier components of the KS wavefunctions stored in a PBLAS matrix
   ! Bands are distributed inside the gtau_comm in a round-robin fashion.

   type(processor_scalapack) :: gtau_slkproc

   integer :: ugb_nband = -1

   type(vcgen_t) :: vcgen
   ! Object used to compute vc(q,g)

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

   integer,allocatable :: kbz2ibz_symrel(:,:)
    ! (6, nkbz)
    ! Mapping kBZ to IBZ (symrel conventions) TODO: To be removed

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
   ! (3,3,my_ntau)
   ! Head of the irred. polarizability in i.omega space.
   ! Note that spins have been summed over.

   complex(dp),allocatable :: chi0_uwing_myw(:,:,:), chi0_lwing_myw(:,:,:)
   ! (3, npw_chi_gamma, my_ntau)
   ! Upper wings of the irred. polarizability in i omega space.
   ! Note that spins have been summed over.

   type(wfdgw_t) :: kcalc_wfd
   ! wavefunction descriptor with the KS states where QP corrections are wanted.

   type(melements_t) :: ks_me !, qp_me
   ! Matrix elements of the different potentials in the KS basis set.

   type(degtab_t),allocatable :: degtab(:,:)
   ! (nkcalc, nsppol)
   ! Table used to average QP results in the degenerate subspace if symsigma == 1

   complex(dp),allocatable :: x_mat(:,:,:,:)
   ! (b1gw:b2gw, b1gw:b2gw, nkcalc, nsppol*nsig_ab)
   ! Matrix elements of $\<nks|\Sigma_x|nk's\>$ with
   ! b1gw = minval(gwr%bstart_ks); b2gw = maxval(gwr%bstop_ks)

   !type(pstat_t) :: ps

 contains

   procedure :: init => gwr_init
   ! Initialize the object.

   procedure :: rotate_gpm => gwr_rotate_gpm
   ! Reconstruct the Green's functions in the BZ from the IBZ.

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

   procedure :: malloc_free_mats => gwr_malloc_free_mats
   ! Allocate/Deallocate matrices for G/tchi/Sigma

   procedure :: free => gwr_free
   ! Free memory.

   procedure :: print => gwr_print
   ! Print info on the object.

   procedure :: print_mem => gwr_print_mem
   ! Print memory required by PBLAS matrices.

   procedure :: print_trace => gwr_print_trace
   ! Print trace of matrices for testing purposes.

   procedure :: load_kcalc_wfd => gwr_load_kcalc_wfd
   ! Load the KS states for Sigma_nk from the WFK file

   procedure :: read_ugb_from_wfk => gwr_read_ugb_from_wfk
   ! Read wavefunctions from WFK file.

   procedure :: build_green => gwr_build_green
   ! Build Green's functions in imaginary time from the gwr%ugb matrices stored in memory.

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
   ! Compute head and wings of chi0

   procedure :: build_sigxme => gwr_build_sigxme
   ! Compute matrix elements of the exchange part.

   procedure :: get_u_ngfft => gwr_get_u_ngfft
   ! Compute FFT mesh from boxcutmin

   procedure :: run_g0w0 => gwr_run_g0w0
   ! Compute QP corrections with one-shot G0W0.

   procedure :: run_energy_scf => gwr_run_energy_scf
   ! Compute QP corrections with energy-only self-consistent GW

   procedure :: check_scf_cycle => gwr_check_scf_cycle
   ! Check SCF cycle for convergence.

   procedure :: ncwrite_tchi_wc => gwr_ncwrite_tchi_wc
   ! Write tchi or wc to netcdf file

 end type gwr_t
!!***

 interface gsph2box
   module procedure gsph2box_dpc
   module procedure gsph2box_spc
 end interface gsph2box

 interface box2gsph
   module procedure box2gsph_dpc
   module procedure box2gsph_spc
 end interface box2gsph

 interface get_1d_sc_phases
   module procedure get_1d_sc_phases_dpc
   module procedure get_1d_sc_phases_spc
 end interface get_1d_sc_phases

 interface sc_sum
   module procedure sc_sum_dpc
   module procedure sc_sum_spc
 end interface sc_sum

 ! Handy named costants (private stuff)
 integer,private,parameter :: PRINT_MODR = 500

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha

 integer,private,parameter :: istwfk1 = 1, ndat1 = 1, me_fft0 = 0, paral_fft0 = 0, nproc_fft1 = 1
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
 integer :: my_nshiftq, iq_bz, iq_ibz, npw_, ncid !, ig
 integer :: comm_cart, me_cart, ierr, all_nproc, my_rank, qprange_, gap_err, ncerr, omp_nt !np_work,
 integer :: cnt, ikcalc, ndeg, mband, bstop, nbsum !, it, iw ! jj,
 integer :: ik_ibz, ik_bz, isym_k, trev_k, g0_k(3)
 integer :: ip_g, ip_k, ip_t, ip_s, np_g, np_k, np_t, np_s
 real(dp) :: cpu, wall, gflops, te_min, te_max, wmax, vc_ecut, delta, abs_rerr, exact_int, eval_int
 real(dp) :: prev_efficiency, prev_speedup
 logical :: isirr_k, changed, q_is_gamma, reorder
 character(len=5000) :: msg
 type(krank_t) :: qrank, krank_ibz
 type(gaps_t) :: ks_gaps
 type(est_t) :: est
 !type(pstat_t) :: ps
!arrays
 integer :: qptrlatt(3,3), dims_kgts(ndims), try_dims_kgts(ndims), indkk_k(6,1), units(2)
 integer,allocatable :: gvec_(:,:),degblock(:,:), degblock_all(:,:,:,:), ndeg_all(:,:), iwork(:,:), got(:)
 real(dp) :: my_shiftq(3,1), kk_ibz(3), kk_bz(3), qq_bz(3), qq_ibz(3), kk(3), tsec(2)
 real(dp),allocatable :: wtk(:), kibz(:,:)
 logical :: periods(ndims), keepdim(ndims)

! *************************************************************************

 all_nproc = xmpi_comm_size(input_comm); my_rank = xmpi_comm_rank(input_comm)
 units = [std_out, ab_out]

 call cwtime(cpu, wall, gflops, "start")
 call timab(1920, 1, tsec)

 gwr%dtset => dtset
 gwr%dtfil => dtfil
 gwr%cryst => cryst
 gwr%psps => psps
 gwr%pawtab => pawtab
 gwr%ks_ebands => ks_ebands
 gwr%kibz => ks_ebands%kptns
 gwr%wtk => ks_ebands%wtk
 gwr%mpi_enreg => mpi_enreg

 ! Initialize qp_ebands with KS values.
 call ebands_copy(ks_ebands, gwr%qp_ebands)
 call ebands_copy(ks_ebands, gwr%qp_ebands_prev)

 ABI_MALLOC(gwr%ks_vbik, (gwr%ks_ebands%nkpt, gwr%ks_ebands%nsppol))
 gwr%ks_vbik(:,:) = ebands_get_valence_idx(gwr%ks_ebands)

 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%nspden = dtset%nspden
 gwr%natom = dtset%natom
 gwr%usepaw = dtset%usepaw

 gwr%use_supercell_for_tchi = .True.
 !if (gwr%dtset%useria == 1) gwr%use_supercell_for_tchi = .False.
 gwr%use_supercell_for_sigma = .True.
 !if (gwr%dtset%userib == 1) gwr%use_supercell_for_sigma = .False.

 if (dtset%gw_nqlwl /= 0) gwr%q0 = dtset%gw_qlwl(:, 1)
 call wrtout(std_out, sjoin(" Using q0:", ktoa(gwr%q0), "for long-wavelenght limit"))

 mband = ks_ebands%mband
 nbsum = dtset%nband(1)
 ABI_CHECK_IRANGE(nbsum, 1, mband, "Invalid nbsum")

 ! Define frequency mesh for sigma(w_real) and spectral functions.
 ! Note that in GWR computing quantities on the real-axis is really cheap
 ! so we can use very dense meshes without affecting performance.
 ! The default for nfresp and freqspmax is zero.
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
   ABI_ERROR("GWR code requires ngkpt with one shift!")
 end if
 gwr%ngkpt = get_diag(ks_ebands%kptrlatt)

 ! Note symrec convention
 ebands_timrev = kpts_timrev_from_kptopt(ks_ebands%kptopt)
 krank_ibz = krank_from_kptrlatt(gwr%nkibz, kibz, ks_ebands%kptrlatt, compute_invrank=.False.)

 ABI_MALLOC(gwr%kbz2ibz, (6, gwr%nkbz))
 if (kpts_map("symrec", ebands_timrev, cryst, krank_ibz, gwr%nkbz, gwr%kbz, gwr%kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 ! Order kbz by stars and rearrange entries in kbz2ibz table.
 call kpts_pack_in_stars(gwr%nkbz, gwr%kbz, gwr%kbz2ibz)

 if (my_rank == master) then
   call kpts_map_print(units, " Mapping kBZ --> kIBZ", "symrec", gwr%kbz, kibz, gwr%kbz2ibz, gwr%dtset%prtvol)
 end if

 !call get_ibz2bz(gwr%nkibz, gwr%nkbz, gwr%kbz2ibz, kibz2bz, ierr)
 !ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for k-points")

 ! Table with symrel conventions for the symmetrization of the wfs.
 ABI_MALLOC(gwr%kbz2ibz_symrel, (6, gwr%nkbz))
 if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, gwr%nkbz, gwr%kbz, gwr%kbz2ibz_symrel) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 ! Setup qIBZ, weights and BZ.
 ! Always use q --> -q symmetry even in systems without inversion
 ! TODO: Might add input variable to rescale the q-mesh.
 my_nshiftq = 1; my_shiftq = zero; qptrlatt = ks_ebands%kptrlatt
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
                             gwr%nqibz, gwr%qibz, gwr%wtq, gwr%nqbz, gwr%qbz)
                             !new_kptrlatt=gwr%qptrlatt, new_shiftk=gwr%qshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ABI_CHECK(all(abs(gwr%qibz(:,1)) < tol16), "First qpoint in qibz should be Gamma!")
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
   call kpts_map_print(units, " Mapping qBZ --> qIBZ", "symrec", gwr%qbz, gwr%qibz, gwr%qbz2ibz, gwr%dtset%prtvol)
 end if

 ! ==========================
 ! Setup k-points in Sigma_nk
 ! ==========================
 ks_gaps = ebands_get_gaps(ks_ebands, gap_err)
 if (my_rank == master) then
   !call ebands_print(ks_ebands, header="KS band structure", unit=std_out, prtvol=gwr%dtset%prtvol)
   !call ebands_print_gaps(ks_ebands, ab_out, header="KS gaps (Fermi energy set to zero)")
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

 ! Include all degenerate states and map kcalc to the ibz.
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

#ifdef __HAVE_GREENX
 call gx_minimax_grid(gwr%ntau, te_min, te_max, &  ! in
                      gwr%tau_mesh, gwr%tau_wgs, gwr%iw_mesh, gwr%iw_wgs, & ! out args allocated by the routine.
                      gwr%cosft_wt, gwr%cosft_tw, gwr%sinft_wt, &
                      gwr%ft_max_error, gwr%cosft_duality_error, ierr)

 if (ierr /= 0) then
   call gx_get_error_message(msg)
   ABI_ERROR(msg)
 end if

 ! FIXME: Here we need to rescale the weights because greenx convention is not what we expect!
 gwr%iw_wgs(:) = gwr%iw_wgs(:) / four

 call wrtout(std_out, sjoin(" Max_{ij} |CT CT^{-1} - I|", ftoa(gwr%cosft_duality_error)))
 call wrtout(std_out, sjoin(" ft_max_error", ltoa(gwr%ft_max_error)))
#else
 ABI_ERROR("GWR code requires Green-X library!")
#endif

 if (gwr%comm%me == 0) then
   write(std_out, "(3a)")ch10, " Computing F(delta) = \int_0^{\infty} dw / (w^2 + delta^2) = pi/2/delta ", ch10
   write(std_out, "(*(a12,2x))")"delta", "numeric", "exact", "abs_rerr (%)"
   do ii=1,10
     delta = (ii * te_min)
     eval_int = sum(gwr%iw_wgs(:) / (gwr%iw_mesh(:)**2 + delta**2))
     exact_int = pi / (two * delta)
     abs_rerr = 100 * abs(eval_int - exact_int) / exact_int
     write(std_out, "(*(es12.5,2x))") delta, eval_int, exact_int, abs_rerr
   end do

   write(std_out, "(3a)")ch10," Computing F(w) = \int_0^{\infty} e^{-w tau} dtau", ch10
   write(std_out, "(*(a12,2x))")"w", "numeric", "exact", "abs_rerr (%)"
   do itau=1,gwr%ntau
     eval_int = sum(gwr%tau_wgs(:) * exp(-gwr%tau_mesh(:) * gwr%iw_mesh(itau)))
     exact_int = one / gwr%iw_mesh(itau)
     abs_rerr = 100 * abs(eval_int - exact_int) / exact_int
     write(std_out, "(*(es12.5,2x))") gwr%iw_mesh(itau), eval_int, exact_int, abs_rerr
   end do
   write(std_out, "(a)")

 endif
 !stop

 call ks_gaps%free()

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
              cryst%symrel, cryst%tnons, use_gpu_cuda=gwr%dtset%use_gpu_cuda, unit=dev_null)
   gwr%green_mpw = max(gwr%green_mpw, npw_)
 end do

 gwr%tchi_mpw = -1
 do iq_bz=1,gwr%nqbz
   qq_bz = gwr%qbz(:, iq_bz)
   call get_kg(qq_bz, istwfk1, dtset%ecuteps, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%gwr_boxcutmin, dtset%chksymtnons, dtset%ecuteps, cryst%gmet, &
              qq_bz, me_fft0, gwr%g_mgfft, gwr%g_nfft, gwr%g_ngfft, nproc_fft1, cryst%nsym, &
              paral_fft0, cryst%symrel, cryst%tnons, use_gpu_cuda=gwr%dtset%use_gpu_cuda, unit=dev_null)
   gwr%tchi_mpw = max(gwr%tchi_mpw, npw_)
   if (iq_bz == 1) then
     ABI_CHECK(all(abs(qq_bz) < tol16), "First qpoint in qbz should be Gamma!")
   end if
 end do

 ! For the time being no augmentation
 gwr%g_ngfft(4:6) = gwr%g_ngfft(1:3)

 ! Define batch sizes for FFT transforms, use multiples of OpenMP threads.
 omp_nt = xomp_get_num_threads(open_parallel=.True.)
 gwr%uc_batch_size = omp_nt; gwr%sc_batch_size = omp_nt
 !gwr%uc_batch_size = 4; gwr%sc_batch_size = 2
 if (gwr%dtset%use_gpu_cuda /= 0) then
   !gwr%uc_batch_size = 4; gwr%sc_batch_size = 2
 end if

 if (my_rank == master) then
   call print_ngfft(gwr%g_ngfft, header="FFT mesh for Green's function", unit=std_out)
   call print_ngfft(gwr%g_ngfft, header="FFT mesh for Green's function", unit=ab_out)
   call wrtout(std_out, sjoin(" FFT uc_batch_size:", itoa(gwr%uc_batch_size)))
   call wrtout(std_out, sjoin(" FFT sc_batch_size:", itoa(gwr%sc_batch_size)))
 end if

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! Here we define the following quantities:
 !
 !  - np_k, np_g, np_t, np_s
 !  - gwr%comm and gwr%idle_proc
 !
 ! NB: Do not use input_comm after this section as idle processors return immediately.

 if (any(dtset%gwr_np_kgts /= 0)) then
   ! Use MPI parameters from input file.
   np_k = dtset%gwr_np_kgts(1); np_g = dtset%gwr_np_kgts(2); np_t = dtset%gwr_np_kgts(3); np_s = dtset%gwr_np_kgts(4)

   !call xmpi_comm_multiple_of(product(dtset%gwr_np_kgts), input_comm, gwr%idle_proc, gwr%comm)
   !if (gwr%idle_proc) return
   gwr%comm = xcomm_from_mpi_int(input_comm)
   all_nproc = gwr%comm%nproc

 else
   ! Automatic grid generation.
   !
   !   Priorities        |  MPI Scalability                | Memory
   ! ==================================================================================================
   !   spin (if any)     |  excellent                      | scales
   !   g/r (PBLAS)       |  network-intensive              ! scales
   !   tau               |  excellent                      | scales
   !   kbz               |  newtwork-intensive             | scales (depends on the BZ -> IBZ mapping)

   gwr%comm = xcomm_from_mpi_int(input_comm)
   all_nproc = gwr%comm%nproc
   !call xmpi_comm_multiple_of(gwr%ntau * gwr%dtset%nsppol, input_comm, gwr%idle_proc, gwr%comm)
   !if (gwr%idle_proc) return
   !all_nproc = xmpi_comm_size(gwr%comm)

#if 1
   ! Start from a configuration that minimizes memory i.e use all procs for g-parallelism,
   ! then check whether it's possible to move some procs to the other levels
   ! without spoiling parallel efficiency and/or increasing memory per MPI proc.
   ! Only master rank works here for consistency reasons.
   if (my_rank == master) then
     dims_kgts = [1, all_nproc, 1, 1]
     est = estimate(gwr, dims_kgts)
     prev_efficiency = est%efficiency; prev_speedup = est%speedup
     call wrtout(units, sjoin("- Optimizing MPI grid with mem_per_cpu_mb:", ftoa(mem_per_cpu_mb), "[Mb]"), pre_newlines=1)
     call wrtout(units, "- Use `abinit run.abi --mem-per-cpu=4G` to set mem_per_cpu_mb in the submission script")
     write(msg, "(a,4(a4,2x),3(a12,2x))") "- ", "np_k", "np_g", "np_t", "np_s", "memb_per_cpu", "efficiency", "speedup"
     call wrtout(units, msg)
     ip_k = dims_kgts(1); ip_g = dims_kgts(2); ip_t = dims_kgts(3); ip_s = dims_kgts(4)
     write(msg, "(a,4(i4,2x),3(es12.5,2x))") "- ", ip_k, ip_g, ip_t, ip_s, est%mem_total, est%efficiency, est%speedup
     call wrtout(units, msg)

     do ip_s=1,gwr%nsppol
       do ip_t=1,gwr%ntau
         if (mod(gwr%ntau, ip_t) /= 0) cycle ! ip_t should divide gwr%ntau.
         do ip_k=1,gwr%nkbz
           if (mod(gwr%nkbz, ip_k) /= 0) cycle ! ip_k is should divide gwr%nkbz.
           do ip_g=1,gwr%green_mpw
             try_dims_kgts = [ip_k, ip_g, ip_t, ip_s]
             if (product(try_dims_kgts) /= all_nproc .or. all(try_dims_kgts == dims_kgts)) cycle
             !npwps
             !ABI_CHECK(block_dist_1d(npwsp, ip_g, col_bsize, msg), msg)
             est = estimate(gwr, try_dims_kgts)
             !if (est%mem_total < mem_per_cpu_mb * 0.8_dp .and. est%efficiency > prev_efficiency) then
             if (est%mem_total < mem_per_cpu_mb * 0.8_dp .and. est%speedup > prev_speedup) then
               prev_efficiency = est%efficiency; prev_speedup = est%speedup; dims_kgts = try_dims_kgts
             end if
             write(msg,"(a,4(i4,2x),3(es12.5,2x))")"- ", ip_k, ip_g, ip_t, ip_s, est%mem_total, est%efficiency, est%speedup
             call wrtout(units, msg)
           end do
         end do
       end do
     end do

     est = estimate(gwr, dims_kgts)
     call wrtout(units, "-")
     call wrtout(units, "- Selected MPI grid:")
     ip_k = dims_kgts(1); ip_g = dims_kgts(2); ip_t = dims_kgts(3); ip_s = dims_kgts(4)
     write(msg, "(a,4(a4,2x),3(a12,2x))") "- ", "np_k", "np_g", "np_t", "np_s", "memb_per_cpu", "efficiency", "speedup"
     call wrtout(units, msg)
     write(msg, "(a,4(i4,2x),3(es12.5,2x))")"- ", ip_k, ip_g, ip_t, ip_s, est%mem_total, est%efficiency, est%speedup
     call wrtout(units, msg, newlines=1)
     call est%print(units)

     !call ps%from_pid()
     !call ps%print([std_out])
   end if ! master

   call xmpi_bcast(dims_kgts, master, gwr%comm%value, ierr)
   np_k = dims_kgts(1); np_g = dims_kgts(2); np_t = dims_kgts(3); np_s = dims_kgts(4)

#else
   ! Determine number of processors for the spin axis. if all_nproc is odd, spin is not distributed when nsppol == 2
   np_s = 1
   if (gwr%nsppol == 2 .and. all_nproc > 1) then
     if (mod(all_nproc, 2) == 0) then
       np_s = 2
     else
       ABI_WARNING("When nsppol == 2, it's reconmended to use an even number of MPI nprocs!")
     end if
   end if
   np_work = all_nproc / np_s

   ! Try to find divisor of ntau and np_work
   do ii=np_work,1,-1
     if (mod(gwr%ntau, ii) == 0 .and. mod(np_work, ii) == 0) exit
   end do

   if (ii == 1 .and. np_work > 1) then
     ! No divisor found.
     if (gwr%nkbz > 1) then
       ! Give priority to tau/kbz
       call xmpi_distrib_2d(np_work, "12", gwr%ntau, gwr%nkbz, np_t, np_k, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np_work), "with priority: tau/kbz"))
     else
       ! Give priority to tau/g
       call xmpi_distrib_2d(np_work, "12", gwr%ntau, gwr%green_mpw, np_t, np_k, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np_work), "with priority: tau/g"))
     end if

   else
     ! ii divides ntau and np_work.
     np_t = ii; np_work = np_work / np_t
     ! Init values assuming Gamma-only sampling.
     np_k = 1; np_g = np_work

     if (gwr%nkbz > 1) then
       do ii=np_work,1,-1
         if (mod(gwr%nkbz, ii) == 0 .and. mod(np_work, ii) == 0) exit
       end do
       if (ii == 1 .and. np_work > 1) then
         call xmpi_distrib_2d(np_work, "12", gwr%nkbz, gwr%green_mpw, np_k, np_k, ierr)
         ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np_work), "with priority: k/g"))
       else
         ! In this case, ii divides both nkbz and np_work.
         np_k = ii; np_g = np_work / ii
       end if
     end if

   end if
#endif

 end if

 ! ================================
 ! Build MPI grid and communicators
 ! ================================
 dims_kgts = [np_k, np_g, np_t, np_s]
 gwr%dtset%gwr_np_kgts = dims_kgts
 periods(:) = .False.; reorder = .True.

 ! Consistency check.
 if (product(dims_kgts) /= all_nproc) then
   write(msg, "(a,i0,3a, 5(a,1x,i0))") &
     "Cannot create 4D Cartesian grid with total nproc: ", all_nproc, ch10, &
     "Idle MPI processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
     "g_nproc (", np_g, ") x tau_nproc (", np_t, ") x kpt_nproc (", np_k,") x spin_nproc (", np_s, ") != ", all_nproc
   ABI_ERROR(msg)
 end if

#ifdef HAVE_MPI
 call MPI_CART_CREATE(gwr%comm%value, ndims, dims_kgts, periods, reorder, comm_cart, ierr)

 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, gwr%coords_gtks, ierr)

 ! k-point communicator
 keepdim = .False.; keepdim(1) = .True.; call gwr%kpt_comm%from_cart_sub(comm_cart, keepdim)
 ! g-communicator
 keepdim = .False.; keepdim(2) = .True.; call gwr%g_comm%from_cart_sub(comm_cart, keepdim)
 ! tau-communicator
 keepdim = .False.; keepdim(3) = .True.; call gwr%tau_comm%from_cart_sub(comm_cart, keepdim)
 ! spin-communicator
 keepdim = .False.; keepdim(4) = .True.; call gwr%spin_comm%from_cart_sub(comm_cart, keepdim)
 ! Communicator for the (X, g, tau, X) 2D grid.
 keepdim = .False.; keepdim(2) = .True.; keepdim(3) = .True.; call gwr%gtau_comm%from_cart_sub(comm_cart, keepdim)
 ! Communicator for the (k, g, X, X) 2D grid.
 keepdim = .False.; keepdim(1) = .True.; keepdim(2) = .True.; call gwr%kg_comm%from_cart_sub(comm_cart, keepdim)
 ! Communicator for the (k, g, tau, X) 3D subgrid.
 keepdim = .True.; keepdim(4) = .False.; call gwr%kgt_comm%from_cart_sub(comm_cart, keepdim)
 ! Communicator for the (k, X, tau, spin) 3D subgrid.
 keepdim = .True.; keepdim(2) = .False.; call gwr%kts_comm%from_cart_sub(comm_cart, keepdim)
 call xmpi_comm_free(comm_cart)
#endif

 !call gwr%kpt_comm%print_names()
 !call gwr%g_comm%print_names()

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
 ABI_CHECK(gwr%my_nspins > 0, "my_nspins == 0, decrease number of MPI procs for spin level")

 ! Distribute k-points in full BZ and build redirection tables.
 ! Finally, find the number of IBZ k-points stored by this MPI rank.

 call xmpi_split_block(gwr%nkbz, gwr%kpt_comm%value, gwr%my_nkbz, gwr%my_kbz_inds)
 ABI_CHECK(gwr%my_nkbz > 0, "my_nkbz == 0, decrease number of MPI procs for k-point level")

 ! Compute np_kibz
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

 ! Build table to distribute iterations over ik_ibz as kIBZ might be replicated across MPI procs.
 ABI_ICALLOC(iwork, (gwr%kpt_comm%nproc, gwr%nkibz))
 ABI_ICALLOC(got, (gwr%kpt_comm%nproc))
 do my_iki=1,gwr%my_nkibz
   ik_ibz = gwr%my_kibz_inds(my_iki)
   iwork(gwr%kpt_comm%me + 1, ik_ibz) = 1
 end do
 call xmpi_sum(iwork, gwr%kpt_comm%value, ierr)
 ABI_MALLOC(gwr%itreat_ikibz, (gwr%nkibz))
 gwr%itreat_ikibz = .False.
 do ik_ibz=1,gwr%nkibz
   ii = imin_loc(got, mask=iwork(:, ik_ibz) /= 0); got(ii) = got(ii) + 1
   if (ii == gwr%kpt_comm%me + 1) gwr%itreat_ikibz(ik_ibz) = .True.
 end do
 ABI_FREE(got)
 ABI_FREE(iwork)

 ! Distribute q-points in full BZ, transfer symmetry tables.
 ! Finally find the number of IBZ q-points that should be stored in memory.
 call xmpi_split_block(gwr%nqbz, gwr%kpt_comm%value, gwr%my_nqbz, gwr%my_qbz_inds)

 ! Compute np_qibz
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

 ! Build table to distribute iterations over iq_ibz as qIBZ might be replicated.
 ABI_ICALLOC(iwork, (gwr%kpt_comm%nproc, gwr%nqibz))
 ABI_ICALLOC(got, (gwr%kpt_comm%nproc))
 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   iwork(gwr%kpt_comm%me + 1, iq_ibz) = 1
 end do
 call xmpi_sum(iwork, gwr%kpt_comm%value, ierr)

 ABI_MALLOC(gwr%itreat_iqibz, (gwr%nqibz))
 gwr%itreat_iqibz = .False.
 do iq_ibz=1,gwr%nqibz
   ii = imin_loc(got, mask=iwork(:, iq_ibz) /= 0); got(ii) = got(ii) + 1
   if (ii == gwr%kpt_comm%me + 1) gwr%itreat_iqibz(iq_ibz) = .True.
 end do
 ABI_FREE(got)
 ABI_FREE(iwork)

 ! TODO: MC technique does not seem to work as expected, even in the legacy code.
 vc_ecut = max(dtset%ecutsigx, dtset%ecuteps)
 call gwr%vcgen%init(cryst, ks_ebands%kptrlatt, gwr%nkbz, gwr%nqibz, gwr%nqbz, gwr%qbz, &
                     dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, vc_ecut, gwr%comm%value)

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
   ! Also, sort the g-vectors by |q+g|^2/2 when q is in the IBZ to facilitate
   ! the extrapolation of the RPA energy as a function of ecut_chi
   call gwr%tchi_desc_qibz(iq_ibz)%init(qq_ibz, istwfk1, dtset%ecuteps, gwr, kin_sorted=.True.)

   ! Compute sqrt(vc(q,G))
   associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))
   call desc_q%get_vc_sqrt(qq_ibz, q_is_gamma, gwr, gwr%gtau_comm%value)
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
   call gwr%print(unit=ab_out); call gwr%print(unit=std_out)

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

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "symsigma", "scf_iteration"])
   NCF_CHECK(ncerr)

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
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "symsigma", "scf_iteration"], &
     [gwr%dtset%symsigma, gwr%scf_iteration])
   NCF_CHECK(ncerr)

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
 call timab(1920, 2, tsec)

contains
integer function vid(vname)
  character(len=*),intent(in) :: vname
  vid = nctk_idname(ncid, vname)
end function vid

end subroutine gwr_init
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/estimate
!! NAME
!! estimate
!!
!! FUNCTION
!!  Estimate memory requirements and the parallel speedup of a given `np_kgts` configuration.
!!
!! SOURCE

type(est_t) pure function estimate(gwr, np_kgts) result(est)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 integer,intent(in) :: np_kgts(4)

!Local variables-------------------------------
 real(dp) :: np_k, np_g, np_t, np_s, w_k, w_g, w_t, w_s, np_tot

! *************************************************************************

! Use real quantities to avoid integer division
 np_k = np_kgts(1); np_g = np_kgts(2); np_t = np_kgts(3); np_s = np_kgts(4)
 np_tot = product(real(np_kgts))

 ! NB: array dimensioned with nkibz and nqibz do not scale as 1/np_k as we distribute the BZ, IBZ points might be replicated.

 ! Resident memory in Mb for G(g,g',+/-tau) and chi(g,g',tau)
 est%mem_green_gg = two * two * (one*gwr%nspinor*gwr%green_mpw)**2 * two*gwr%ntau * gwr%nkibz * gwr%nsppol * gwp*b2Mb / np_tot
 est%mem_chi_gg = two * (one*gwr%tchi_mpw)**2 * gwr%ntau * gwr%nqibz * gwp*b2Mb / (np_g * np_t * np_k)
 est%mem_ugb = two * gwr%green_mpw * gwr%nspinor * gwr%dtset%nband(1) * gwr%nkibz * gwr%nsppol * gwp*b2Mb / np_tot

 ! Temporary memory allocated inside the tau loops.
 ! This is the chunck we have to minimize by increasing np_g and/or np_k to avoid going OOM.
 est%mem_green_rg = two * two * gwr%nspinor**2 * gwr%green_mpw * gwr%g_nfft * gwr%nkbz * gwr%nsppol * gwp*b2Mb / (np_g * np_k)
 est%mem_chi_rg = two * gwr%tchi_mpw * gwr%g_nfft * gwr%nqbz * gwp*b2Mb / (np_g * np_k)

 est%mem_total = est%mem_green_gg + est%mem_chi_gg + est%mem_ugb + est%mem_green_rg + est%mem_chi_rg

 ! Estimate speedup and parallel efficiency using heuristic weights. Note g_nfft instead of green_mpw.
 w_k = 0.799_dp; w_g = 0.899_dp; w_t = 1.1_dp; w_s = 1.2_dp
 ! Promote kpt parallelism under particular circustamnces.
 if (gwr%nkbz > 4**3) w_k = w_g + tol2 * merge(+1, -1, np_k < 5)

 est%speedup = speedup(gwr%nkbz, nint(np_k), w_k) * speedup(gwr%g_nfft, nint(np_g), w_g) * &
               speedup(gwr%ntau, nint(np_t), w_t) * speedup(gwr%nsppol, nint(np_s), w_s)
 est%efficiency = est%speedup / np_tot

contains

real(dp) pure function speedup(size, np, weight)
 ! Expected speedup for a `size` problem and `np` processes
 integer,intent(in) :: size, np
 real(dp),intent(in) :: weight
 if (np == 1) then
   speedup = one
 else
   speedup = (weight*size) / (one* ((size / np) + merge(0, 1, mod(size, np) == 0)))
 end if
end function speedup

end function estimate
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/est_print
!! NAME
!! est_print
!!
!! FUNCTION
!!
!! SOURCE

subroutine est_print(est, units)

!Arguments ------------------------------------
 class(est_t), intent(in) :: est
 integer,intent(in) :: units(:)

!Local variables-------------------------------
 !real(dp) :: np_k, np_g, np_t, np_s, w_k, w_g, w_t, w_s, np_tot

! *************************************************************************

 call wrtout(units, "- Resident memory in Mb for G(g,g',+/-tau) and chi(g,g',tau):")
 !est%mem_green_gg
 !est%mem_chi_gg
 !est%mem_ugb

 call wrtout(units, "- Temporary memory allocated inside the tau loops:")
 !est%mem_green_rg
 !est%mem_chi_rg

end subroutine est_print
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_malloc_free_mats
!! NAME
!! gwr_malloc_free_mats
!!
!! FUNCTION
!! Allocate/Free PBLAS matrices according to `what` for the set of k/q-points selected by `mask_ibz`.
!!
!! SOURCE

subroutine gwr_malloc_free_mats(gwr, mask_ibz, what, action)

!Arguments ------------------------------------
 class(gwr_t), target, intent(inout) :: gwr
 integer,intent(in) :: mask_ibz(:)
 character(len=*),intent(in) :: what, action

!Local variables-------------------------------
 integer :: my_is, my_it, ipm, npwsp, col_bsize, itau, spin, ik_ibz, iq_ibz
 type(__slkmat_t), pointer :: mat
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
           if (action == "malloc") call gt(ipm)%init(npwsp, npwsp, gwr%g_slkproc, istwfk1, size_blocs=[-1, col_bsize])
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
         if (action == "malloc") call mat%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[-1, col_bsize])
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
           if (action == "malloc") call sigc(ipm)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[-1, col_bsize])
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

end subroutine gwr_malloc_free_mats
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_free
!! NAME
!! gwr_free
!!
!! FUNCTION
!!  Free dynamic memory
!!
!! SOURCE

subroutine gwr_free(gwr)

!Arguments ------------------------------------
 class(gwr_t), intent(inout) :: gwr

! *************************************************************************

 ABI_SFREE(gwr%ks_vbik)
 ABI_SFREE(gwr%kbz)
 ABI_SFREE(gwr%kbz2ibz)
 ABI_SFREE(gwr%kbz2ibz_symrel)
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
 ABI_SFREE(gwr%itreat_ikibz)
 ABI_SFREE(gwr%np_qibz)
 ABI_SFREE(gwr%itreat_iqibz)
!#ifdef __HAVE_GREENX
 ABI_SFREE_NOCOUNT(gwr%tau_mesh)
 ABI_SFREE_NOCOUNT(gwr%tau_wgs)
 ABI_SFREE_NOCOUNT(gwr%iw_mesh)
 ABI_SFREE_NOCOUNT(gwr%iw_wgs)
 ABI_SFREE_NOCOUNT(gwr%cosft_tw)
 ABI_SFREE_NOCOUNT(gwr%cosft_wt)
 ABI_SFREE_NOCOUNT(gwr%sinft_wt)
!#endif
 ABI_SFREE(gwr%kcalc)
 ABI_SFREE(gwr%bstart_ks)
 ABI_SFREE(gwr%bstop_ks)
 ABI_SFREE(gwr%nbcalc_ks)
 ABI_SFREE(gwr%kcalc2ibz)
 ABI_SFREE(gwr%x_mat)

 call ebands_free(gwr%qp_ebands)
 call ebands_free(gwr%qp_ebands_prev)
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

 ! datatypes.
 call gwr%ks_me%free()
 call gwr%vcgen%free()

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
 call gwr%kg_comm%free()
 call gwr%kgt_comm%free()
 call gwr%kts_comm%free()
 call gwr%comm%free()

end subroutine gwr_free
!!***

! Free array of desc_t objects.
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
!!  Load the KS states to compute Sigma_nk from the WFK file
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

 ks_ebands = wfk_read_ebands(wfk_path, gwr%comm%value, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ! TODO: Add more consistency checks e.g. nkibz,...

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
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, gwr%comm%value)

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
!!  Read wavefunctions from the WFK file wfk_path and store them in gwr%ugb (MPI distributed).
!!
!! SOURCE

subroutine gwr_read_ugb_from_wfk(gwr, wfk_path)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=*),intent(in) :: wfk_path

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0 = 0, master = 0
 integer :: mband, min_nband, nkibz, nsppol, my_is, my_iki, spin, ik_ibz, ierr, io_algo
 integer :: npw_k, mpw, istwf_k, il_b, ib, band, iloc !, itau
 integer :: nbsum, npwsp, bstart, bstop, bstep, nb !my_nband ! nband_k,
 real(dp) :: cpu, wall, gflops, cpu_green, wall_green, gflops_green
 character(len=5000) :: msg
 logical :: haveit
 type(ebands_t) :: wfk_ebands
 type(hdr_type) :: wfk_hdr
 type(wfk_t) :: wfk
 type(dataset_type),pointer :: dtset
!arrays
 integer,allocatable :: kg_k(:,:)
 logical,allocatable :: bmask(:)
 real(dp) :: kk_ibz(3), tsec(2)
 real(dp),target,allocatable :: cg_work(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: cg_k(:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 call timab(1921, 1, tsec)

 dtset => gwr%dtset

 wfk_ebands = wfk_read_ebands(wfk_path, gwr%comm%value, out_hdr=wfk_hdr)
 call wfk_hdr%vs_dtset(dtset)
 ! TODO: Add more consistency checks e.g. nkibz,...

 !cryst = wfk_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 nkibz = wfk_ebands%nkpt; nsppol = wfk_ebands%nsppol; mband = wfk_ebands%mband
 min_nband = minval(wfk_ebands%nband)

 nbsum = dtset%nband(1)
 if (nbsum > min_nband) then
   ABI_WARNING(sjoin("WFK file contains", itoa(min_nband), "states while you're asking for:", itoa(nbsum)))
   nbsum = min_nband
 end if
 !call wrtout(std_out, sjoin(" Computing Green's function with nbsum:", itoa(nbsum)), do_flush=.True.)
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
 !     2) May implement trick used in gwst to add empty states approximated with LC of PWs.

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

 call wrtout(std_out, sjoin(" Reading KS states with nbsum", itoa(nbsum), "..."), do_flush=.True.)

 ! Init set of (npwsp, nbsum) PBLAS matrix distributed within the gtau communicator.
 ! and distribute it over bands so that each proc reads a subset of bands in read_band_block
 ! Note size_blocs below that corresponds to roud-robin distribution along band axis.

 ABI_MALLOC(gwr%ugb, (gwr%nkibz, gwr%nsppol))
 gwr%ugb_nband = nbsum

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     npw_k = gwr%green_desc_kibz(ik_ibz)%npw
     npwsp = npw_k * gwr%nspinor
     call gwr%ugb(ik_ibz, spin)%init(npwsp, gwr%ugb_nband, gwr%gtau_slkproc, istwfk1, size_blocs=[-1, 1])
   end do
 end do
 call gwr%print_mem(unit=std_out)

 mpw = maxval(wfk_hdr%npwarr)
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(bmask, (mband))

 io_algo = 2

 if (io_algo == 1) then
   ! This version is very bad on LUMI
   call wrtout(std_out, " Using collective MPI-IO with wfk%read_bmask ...")
   call wfk_open_read(wfk, wfk_path, formeig0, iomode_from_fname(wfk_path), get_unit(), gwr%gtau_comm%value)
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
       npwsp = npw_k * gwr%nspinor

       associate (ugb => gwr%ugb(ik_ibz, spin), desc_k => gwr%green_desc_kibz(ik_ibz))
       ABI_CHECK_IEQ(npw_k, desc_k%npw, "npw_k != desc_k%npw")

       ! use round-robin distribution inside gtau_comm% for IO.
       ! TODO: Optimize wfk_read_bmask and/or read WFK with all procs and master broadcasting.
       bmask = .False.
       do il_b=1, ugb%sizeb_local(2)
         band = ugb%loc2gcol(il_b)
         bmask(band) = .True.
       end do
       ! FIXME: This is wrong is spc
       call c_f_pointer(c_loc(ugb%buffer_cplx), cg_k, shape=[2, npwsp * ugb%sizeb_local(2)])
       call wfk%read_bmask(bmask, ik_ibz, spin, &
                           !xmpio_single, &
                           xmpio_collective, &
                           kg_k=kg_k, cg_k=cg_k)

       ABI_CHECK(all(kg_k(:,1:npw_k) == desc_k%gvec), "kg_k != desc_k%gvec")

       write(msg,'(4x, 3(a,i0),a)')" Read ugb_k: my_iki [", my_iki, "/", gwr%my_nkibz, "] (tot: ", gwr%nkibz, ")"
       call cwtime_report(msg, cpu_green, wall_green, gflops_green)
       end associate
     end do ! my_iki
   end do ! my_is
   call wfk%close()

 else
   ! Master reads and broadcasts. Much faster on lumi
   call wrtout(std_out, " Using IO version based on master reads and brodcasts ...")
   if (gwr%comm%me == master) then
     call wfk_open_read(wfk, wfk_path, formeig0, iomode_from_fname(wfk_path), get_unit(), xmpi_comm_self)
   end if

   ! TODO This to be able to maximize the size of cg_work
   !call xmpi_get_vmrss(vmem_mb, gwr%comm%value)

   do spin=1,gwr%nsppol
     do ik_ibz=1,gwr%nkibz
       call cwtime(cpu_green, wall_green, gflops_green, "start")
       kk_ibz = gwr%kibz(:, ik_ibz)
       npw_k = wfk_hdr%npwarr(ik_ibz)
       istwf_k = wfk_hdr%istwfk(ik_ibz)
       ! TODO
       ABI_CHECK_IEQ(istwf_k, 1, "istwfk_k should be 1")
       npwsp = npw_k * gwr%nspinor

       bstep = memlimited_step(1, nbsum, 2 * npwsp, xmpi_bsize_dp, 1024.0_dp)
       do bstart=1, nbsum, bstep
         bstop = min(bstart + bstep - 1, nbsum)
         nb = bstop - bstart + 1

         ABI_MALLOC(cg_work, (2, npwsp, nb))
         if (gwr%comm%me == master) then
           call c_f_pointer(c_loc(cg_work), cg_k, shape=[2, npwsp * nb])
           call wfk%read_band_block([bstart, bstop], ik_ibz, spin, xmpio_single, kg_k=kg_k, cg_k=cg_k)
         end if

         call xmpi_bcast(kg_k, master, gwr%comm%value, ierr)
         call xmpi_bcast(cg_work, master, gwr%comm%value, ierr)

         ! Copy my portion of cg_work to buffer_cplx if I treat this (spin, ik_ibz).
         if (any(gwr%my_spins == spin) .and. any(gwr%my_kibz_inds == ik_ibz)) then
           associate (ugb => gwr%ugb(ik_ibz, spin), desc_k => gwr%green_desc_kibz(ik_ibz))
           ABI_CHECK(all(kg_k(:,1:npw_k) == desc_k%gvec), "kg_k != desc_k%gvec")
           do band=bstart, bstop
             ib = band - bstart + 1
             call ugb%glob2loc(1, band, iloc, il_b, haveit)
             if (.not. haveit) cycle
             ABI_CHECK_IEQ(iloc, 1, "iloc should be 1")
             ugb%buffer_cplx(:, il_b) = cmplx(cg_work(1,:,ib), cg_work(2,:,ib), kind=gwpc)
           end do
           end associate
         end if

         ABI_FREE(cg_work)
       end do ! bstart

       write(msg,'(4x,2(a,i0),a)')" Read ugb_k: ik_ibz [", ik_ibz, "/", gwr%nkibz, "]"
       call cwtime_report(msg, cpu_green, wall_green, gflops_green)
     end do ! ik_ibz
   end do ! spin
   if (gwr%comm%me == master) call wfk%close()

 end if ! io_algo

 call cwtime_report(" gwr_read_ugb_from_wfk:", cpu, wall, gflops)
 call timab(1921, 2, tsec)

 ABI_FREE(kg_k)
 ABI_FREE(bmask)

 call wfk_hdr%free()
 call gwr%print_mem(unit=std_out)

end subroutine gwr_read_ugb_from_wfk
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_green
!! NAME
!!  gwr_build_green
!!
!! FUNCTION
!!  Build Green's functions in imaginary time from the gwr%ugb matrices stored in memory.
!!  Store only G_k for the IBZ k-points treated by this MPI proc.
!!
!! INPUTS
!!  free_ugb: True if gwr%ugb wavefunctions should be deallocated before returning.
!!
!! SOURCE

subroutine gwr_build_green(gwr, free_ugb)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 logical,intent(in) :: free_ugb

!Local variables-------------------------------
!scalars
 integer :: my_is, my_iki, spin, ik_ibz, band, itau, ipm, il_b, npwsp, isgn
 real(dp) :: f_nk, eig_nk, cpu, wall, gflops, cpu_green, wall_green, gflops_green
 character(len=500) :: msg
 complex(dp) :: gt_cfact
 type(__slkmat_t), target :: work, green
!arrays
 integer :: mask_kibz(gwr%nkibz)
 real(dp) :: tsec(2)
 real(dp),contiguous, pointer :: qp_eig(:,:,:), qp_occ(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 call timab(1922, 1, tsec)

 ! Use KS or QP energies depending on the iteration state.
 if (gwr%scf_iteration == 1) then
   call wrtout([std_out, ab_out], " Building Green's functions from KS orbitals and KS energies...", &
               pre_newlines=2, newlines=1, do_flush=.True.)

   qp_eig => gwr%ks_ebands%eig; qp_occ => gwr%ks_ebands%occ
   msg = sjoin("Fermi energy is not set to zero! fermie:", ftoa(gwr%ks_ebands%fermie))
   ABI_CHECK(abs(gwr%ks_ebands%fermie) < tol12, msg)

   ! Allocate Green's functions if this is the first iteration
   mask_kibz = 0; mask_kibz(gwr%my_kibz_inds(:)) = 1
   call gwr%malloc_free_mats(mask_kibz, "green", "malloc")

 else
   call wrtout([std_out, ab_out], " Building Green's functions from KS orbitals and QP energies...", &
               pre_newlines=2, newlines=1, do_flush=.True.)
   qp_eig => gwr%qp_ebands%eig; qp_occ => gwr%qp_ebands%occ
   msg = sjoin("Fermi energy is not set to zero! fermie:", ftoa(gwr%qp_ebands%fermie))
   ABI_CHECK(abs(gwr%qp_ebands%fermie) < tol12, msg)
 end if

 ABI_CHECK(allocated(gwr%ugb), "gwr%ugb array should be allocated!")

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz ! my k-points in IBZ
     call cwtime(cpu_green, wall_green, gflops_green, "start")
     ik_ibz = gwr%my_kibz_inds(my_iki)
     associate (ugb => gwr%ugb(ik_ibz, spin), desc_k => gwr%green_desc_kibz(ik_ibz))
     npwsp = desc_k%npw * gwr%nspinor

     call ugb%copy(work)
     !call ugb%change_size_blocs(work, size_blocs=, processor=)
     !call work%copy(green, empty=.True.)

     ! Init output of pzgemm.
     call green%init(npwsp, npwsp, gwr%gtau_slkproc, istwfk1) ! size_blocs=[-1, col_bsize])

     ! We loop over ntau instead of my_ntau as pzgemm is done inside gtau_comm.
     do itau=1,gwr%ntau
       do ipm=1,2
         ! Multiply columns by exponentials in imaginary time.
         work%buffer_cplx = ugb%buffer_cplx

         !$OMP PARALLEL DO PRIVATE(band, f_nk, eig_nk, gt_cfact)
         do il_b=1, work%sizeb_local(2)
           band = work%loc2gcol(il_b)
           f_nk = qp_occ(band, ik_ibz, spin)
           eig_nk = qp_eig(band, ik_ibz, spin)
           gt_cfact = zero
           if (ipm == 2) then
             if (eig_nk < -tol6) gt_cfact = exp(gwr%tau_mesh(itau) * eig_nk)
           else
             if (eig_nk > tol6) gt_cfact = exp(-gwr%tau_mesh(itau) * eig_nk)
           end if

           work%buffer_cplx(:,il_b) = work%buffer_cplx(:,il_b) * sqrt(real(gt_cfact))
         end do ! il_b

         ! Build G(g,g',ipm)
         isgn = merge(1, -1, ipm == 2)
         call slk_pgemm("N", "C", work, isgn * cone_gw, work, czero_gw, green) ! TODO ija=, ijb=

         ! Redistribute data from gtau_comm to g_comm.
         call gwr%gt_kibz(ipm, ik_ibz, itau, spin)%take_from(green)
       end do ! ipm
     end do ! itau

     call work%free()
     call green%free()
     if (free_ugb) call ugb%free()

     write(msg,'(4x, 3(a,i0),a)')" G_ikbz [", my_iki, "/", gwr%my_nkibz, "] (tot: ", gwr%nkibz, ")"
     call cwtime_report(msg, cpu_green, wall_green, gflops_green)
     end associate
   end do ! my_iki
 end do ! my_is

 if (gwr%dtset%prtvol > 0) call gwr_print_trace(gwr, "gt_kibz")
 call cwtime_report(" gwr_build_green:", cpu, wall, gflops)
 call timab(1922, 2, tsec)

 call gwr%print_mem(unit=std_out)

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
 type(__slkmat_t),intent(out) :: gt_pm(2)

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
 !      G_{-k}(-g,-g') = [G_k(g,g')]*

 !ABI_WARNING_IF(trev_k == 0, "green: trev_k /= 0 should be tested")

 ! Rotate gvec, recompute gbound and rotate vc_sqrt
 ! TODO: 1) Handle TR and routine to rotate tchi/W including vc_sqrt
 !       2) Make sure that the FFT box is large enough to accomodate umklapps

 desc_kbz%ig0 = -1
 do ig1=1,desc_kbz%npw
   desc_kbz%gvec(:,ig1) = tsign_k * matmul(gwr%cryst%symrec(:,:,isym_k), desc_kibz%gvec(:,ig1)) - g0_k
   if (all(desc_kbz%gvec(:,ig1) == 0)) desc_kbz%ig0 = ig1
 end do
 desc_kbz%kin_sorted = .False.
 ABI_CHECK(desc_kbz%ig0 /= -1, "Cannot find g=0 after rotation!")

 call sphereboundary(desc_kbz%gbound, desc_kbz%istwfk, desc_kbz%gvec, gwr%g_mgfft, desc_kbz%npw)

 ! Get G_k with k in the BZ.
 tnon = gwr%cryst%tnons(:, isym_k)
 do ipm=1,2
   associate (gk_i => gwr%gt_kibz(ipm, ik_ibz, itau, spin), gk_f => gt_pm(ipm))
   call gk_i%copy(gk_f)
   !!$OMP PARALLEL DO PRIVATE(ig1, g2, ph2, ig1, g2, ph1)
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
 type(__slkmat_t),intent(inout) :: gt_gpr(2, gwr%my_nkbz)

!Local variables-------------------------------
!scalars
 integer :: my_ikf, ik_bz, ig2, ipm, npwsp, col_bsize, idat, ndat
 logical :: k_is_gamma
 real(dp) :: kk_bz(3), cpu, wall, gflops !, mem_mb
 complex(gwpc),allocatable :: ceikr(:)
 character(len=500) :: msg
 type(__slkmat_t) :: rgp, gt_pm(2)
 type(uplan_t) :: uplan_k

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 !mem_mb = two * gwr%my_nkbz * two * gwpc * gwr%g_nfft * gwr%gree_mpw * b2Mb /  gwr%g_slkproc%nbprocs
 !call wrtout(std_out, sjoin("Local memory for Green's functions: ", ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

 ABI_MALLOC(ceikr, (gwr%g_nfft * gwr%nspinor))

 do my_ikf=1,gwr%my_nkbz
   ik_bz = gwr%my_kbz_inds(my_ikf)
   kk_bz = gwr%kbz(:, ik_bz)

   k_is_gamma = normv(kk_bz, gwr%cryst%gmet, "G") < GW_TOLQ0
   if (.not. k_is_gamma) call calc_ceikr(kk_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ceikr)

   ! Get G_kbz(+/- itau) in the BZ.
   call gwr%rotate_gpm(ik_bz, itau, spin, desc_mykbz(my_ikf), gt_pm)

   associate (desc_k => desc_mykbz(my_ikf))
   call uplan_k%init(desc_k%npw, gwr%nspinor, gwr%uc_batch_size, gwr%g_ngfft, desc_k%istwfk, &
                     desc_k%gvec, gwpc, gwr%dtset%use_gpu_cuda)

   do ipm=1,2
     associate (ggp => gt_pm(ipm))

     ! Allocate rgp PBLAS matrix to store G_kbz(r,g')
     npwsp = desc_k%npw * gwr%nspinor
     ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
     call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_k%istwfk, size_blocs=[-1, col_bsize])

     ABI_CHECK_IEQ(size(ggp%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

     ! Perform FFT G_k(g,g') -> G_k(r,g') and store results in rgp.
     do ig2=1, ggp%sizeb_local(2), gwr%uc_batch_size
       ndat = blocked_loop(ig2, ggp%sizeb_local(2), gwr%uc_batch_size)

#if 0
       call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                   gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                   ggp%buffer_cplx(:, ig2), &  ! in
                   rgp%buffer_cplx(:, ig2))    ! out

#else
       call uplan_k%execute_gr(ndat, ggp%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))
#endif

       if (.not. k_is_gamma) then
         ! Multiply by e^{ik.r}
         !$OMP PARALLEL DO
         do idat=0,ndat-1
           rgp%buffer_cplx(:, ig2 + idat) = ceikr(:) * rgp%buffer_cplx(:, ig2 + idat)
         end do
       end if
     end do ! ig2

     ! MPI transpose: G_k(r,g') -> G_k(g',r)
     call rgp%ptrans("N", gt_gpr(ipm, my_ikf), free=.True.)
     end associate
   end do ! ipm

   call uplan_k%free()
   call slk_array_free(gt_pm)
   end associate
 end do ! my_ikf

 ABI_FREE(ceikr)

 call cwtime_report(" gwr_get_myk_green_gpr:", cpu, wall, gflops)

end subroutine gwr_get_myk_green_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_gk_rpr_pm
!! NAME
!!  gwr_get_gk_rpr_pm
!!
!! FUNCTION
!!  Compute G_k(r',r) from G_k(g,g') for k in the BZ and given spin and tau.
!!  Note that output matrix is transposed i.e. (r',r) instead of (r,r').
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_gk_rpr_pm(gwr, ik_bz, itau, spin, gk_rpr_pm)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 integer,intent(in) :: ik_bz, itau, spin
 type(__slkmat_t),intent(inout) :: gk_rpr_pm(2)

!Local variables-------------------------------
!scalars
 integer :: ig2, ipm, npwsp, col_bsize, ir1, ndat
 real(dp) :: cpu, wall, gflops
 type(__slkmat_t) :: rgp, gt_pm(2), gpr
 type(desc_t) :: desc_kbz
 type(uplan_t) :: uplan_k
 character(len=500) :: msg

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ABI_ERROR("This should be tested")

 !kk_bz = gwr%kbz(:, ik_bz)

 ! Get G_k(+/- itau) in the BZ.
 call gwr%rotate_gpm(ik_bz, itau, spin, desc_kbz, gt_pm)

 call uplan_k%init(desc_kbz%npw, gwr%nspinor, gwr%uc_batch_size, gwr%g_ngfft, desc_kbz%istwfk, &
                   desc_kbz%gvec, gwpc, gwr%dtset%use_gpu_cuda)

 do ipm=1,2
   ! Allocate rgp PBLAS matrix to store G(r,g')
   npwsp = desc_kbz%npw * gwr%nspinor
   ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
   call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_kbz%istwfk, size_blocs=[-1, col_bsize])

   associate (ggp => gt_pm(ipm))
   !ABI_CHECK_IEQ(size(ggp%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")
   do ig2=1, ggp%sizeb_local(2), gwr%uc_batch_size
     ndat = blocked_loop(ig2, ggp%sizeb_local(2), gwr%uc_batch_size)

     ! Perform FFT G_k(g,g') -> G_k(r,g') and store results in rgp.
     !call fft_ug(desc_kbz%npw, gwr%g_nfft, gwr%nspinor, ndat, &
     !            gwr%g_mgfft, gwr%g_ngfft, desc_kbz%istwfk, desc_kbz%gvec, desc_kbz%gbound, &
     !            ggp%buffer_cplx(:, ig2), &  ! in
     !            rgp%buffer_cplx(:, ig2))    ! out

     call uplan_k%execute_gr(ndat, ggp%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))
   end do ! ig2
   end associate

   ! MPI transpose: G_k(r,g') -> G_k(g',r) and transform g' index.
   call rgp%ptrans("N", gpr, free=.True.)

   do ir1=1, gpr%sizeb_local(2), gwr%uc_batch_size
     ndat = blocked_loop(ir1, gpr%sizeb_local(2), gwr%uc_batch_size)

     ! Perform FFT G_k(g',r) -> G_k(r',r) and store results in rgp.

     ! FIXME: FFT sign is wrong (should be - instead of + but I need to change the API)
     !call fft_ug(desc_kbz%npw, gwr%g_nfft, gwr%nspinor, ndat, &
     !            gwr%g_mgfft, gwr%g_ngfft, desc_kbz%istwfk, desc_kbz%gvec, desc_kbz%gbound, &
     !            gpr%buffer_cplx(:, ir1), &             ! in
     !            gk_rpr_pm(ipm)%buffer_cplx(:, ir1))    ! out

     call uplan_k%execute_gr(ndat, gpr%buffer_cplx(:, ir1), gk_rpr_pm(ipm)%buffer_cplx(:, ir1))
   end do ! ir1
   call gpr%free()

   ! Rescale
   !gk_rpr_pm(ipm)%buffer_cplx = gk_rpr_pm(ipm)%buffer_cplx * gwr%g_nfft
 end do ! ipm

 call slk_array_free(gt_pm)
 call desc_kbz%free()
 call uplan_k%free()

 call cwtime_report(" gwr_get_gk_rpr_pm:", cpu, wall, gflops)

end subroutine gwr_get_gk_rpr_pm
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_ggp_to_rpr
!! NAME
!!  gwr_ggp_to_rpr
!!
!! FUNCTION
!!  Helper function to FFT transform a two-point function: F_{g,g'} --> F_{r',r}
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
 type(__slkmat_t),intent(in) :: ggp
 type(__slkmat_t),intent(inout) :: rpr

!Local variables-------------------------------
!scalars
 integer :: ig2, npwsp, col_bsize, ir1, ndat
 type(__slkmat_t) :: rgp, gpr
 character(len=500) :: msg
 type(uplan_t) :: uplan_k

! *************************************************************************

 ABI_ERROR("Not Implemented Error")

 ! Allocate intermediate rgp PBLAS matrix to store F(r,g')
 npwsp = desc%npw * gwr%nspinor
 ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
 call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc%istwfk, size_blocs=[-1, col_bsize])

 call uplan_k%init(desc%npw, gwr%nspinor, gwr%uc_batch_size, gwr%g_ngfft, desc%istwfk, &
                   desc%gvec, gwpc, gwr%dtset%use_gpu_cuda)

 do ig2=1, ggp%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ig2, ggp%sizeb_local(2), gwr%uc_batch_size)
   call uplan_k%execute_gr(ndat, ggp%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))
 end do ! ig2

 ! F(r,g') --> F(g',r)
 call rgp%ptrans("N", gpr, free=.True.)

 ! F(g',r) -> F(r',r) and store results in rpr
 ! FIXME: FFT sign is wrong (should be - instead of + but I need to change the API)
 do ir1=1, gpr%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ir1, gpr%sizeb_local(2), gwr%uc_batch_size)
   call uplan_k%execute_gr(ndat, gpr%buffer_cplx(:, ir1), rpr%buffer_cplx(:, ir1))
 end do ! ir1

 call uplan_k%free()
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
 type(__slkmat_t),intent(inout) :: wc_qbz

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, iq_ibz, isym_q, trev_q, g0_q(3), tsign_q
 logical :: isirr_q, q_is_gamma
!arrays
 integer :: g1(3), g2(3)
 real(dp) :: tnon(3), qq_bz(3)
 complex(dp) :: ph2, ph1

! *************************************************************************

 ABI_CHECK(gwr%wc_space == "itau", sjoin("wc_space:", gwr%wc_space, " != itau"))

 qq_bz = gwr%qbz(:, iq_bz)
 q_is_gamma = normv(qq_bz, gwr%cryst%gmet, "G") < GW_TOLQ0

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
   ! Copy the PBLAS matrix in wc_qbz and we are done.
   call gwr%wc_qibz(iq_ibz, itau, spin)%copy(wc_qbz); return
 end if

 !ABI_WARNING_IF(trev_q == 0, "trev_q should be tested")

 ! rotate gvec, recompute gbound and rotate vc_sqrt.
 ! TODO: 1) Handle TR and routine to rotate tchi/W including vc_sqrt
 !       2) Make sure that FFT box is large enough to accomodate umklapps
 desc_qbz%ig0 = -1
 do ig1=1,desc_qbz%npw
   desc_qbz%gvec(:,ig1) = tsign_q * matmul(gwr%cryst%symrec(:,:,isym_q), desc_qibz%gvec(:,ig1)) - g0_q
   if (all(desc_qbz%gvec(:,ig1) == 0)) desc_qbz%ig0 = ig1
 end do
 desc_qbz%kin_sorted = .False.
 ABI_CHECK(desc_qbz%ig0 /= -1, "Cannot find g=0 after rotation!")

 call sphereboundary(desc_qbz%gbound, desc_qbz%istwfk, desc_qbz%gvec, gwr%g_mgfft, desc_qbz%npw)

 ! Compute sqrt(vc(q,G))
 ! TODO: rotate vc_sqrt
 ! vc(Sq, Sg) = vc(q, g)
 ! vc(-q, -g) = vc(q, g)
 call desc_qbz%get_vc_sqrt(qq_bz, q_is_gamma, gwr, gwr%gtau_comm%value)

 ! Get Wc_q with q in the BZ.
 tnon = gwr%cryst%tnons(:, isym_q)
 associate (wq_i => gwr%wc_qibz(iq_ibz, itau, spin), wq_f => wc_qbz)
 call wq_i%copy(wc_qbz)

 !!!$OMP PARALLEL DO PRIVATE(ig2, g2, phs2, ig1, g2, ph1)
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
!!  for each q in the BZ treated by this MPI proc for given `spin` and `itau` index:
!!
!!      1) FFT Transform the first index: Wc(g,g',it) --> Wc(r,g',it)  (local operation)
!!      2) MPI transposition: Wc(r,g',it) --> Wc(g',r,it)
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
 type(__slkmat_t),intent(inout) :: wc_gpr(gwr%my_nqbz)

!Local variables-------------------------------
!scalars
 integer :: my_iqf, iq_bz, ig2, npwsp, col_bsize, idat, ndat
 real(dp) :: cpu, wall, gflops, qq_bz(3)
 logical :: q_is_gamma
 character(len=500) :: msg
 type(__slkmat_t) :: rgp, wc_qbz
 type(uplan_t) :: uplan_q
 complex(gwpc),allocatable :: ceiqr(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 ABI_MALLOC(ceiqr, (gwr%g_nfft * gwr%nspinor))

 do my_iqf=1,gwr%my_nqbz
   iq_bz = gwr%my_qbz_inds(my_iqf)
   qq_bz = gwr%qbz(:, iq_bz)

   q_is_gamma = normv(qq_bz, gwr%cryst%gmet, "G") < GW_TOLQ0
   if (.not. q_is_gamma) call calc_ceikr(qq_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, ceiqr)

   ! Get Wc_q for q in the BZ.
   call gwr%rotate_wc(iq_bz, itau, spin, desc_myqbz(my_iqf), wc_qbz)
   associate (desc_q => desc_myqbz(my_iqf))

   ! Allocate rgp PBLAS matrix to store Wc_q(r, g')
   npwsp = desc_q%npw * gwr%nspinor
   ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
   call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_q%istwfk, size_blocs=[-1, col_bsize])

   !ABI_CHECK_IEQ(size(wc_qbz%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

   call uplan_q%init(desc_q%npw, gwr%nspinor, gwr%uc_batch_size, gwr%g_ngfft, desc_q%istwfk, &
                     desc_q%gvec, gwpc, gwr%dtset%use_gpu_cuda)

   ! FFT. Results stored in rgp
   do ig2=1,wc_qbz%sizeb_local(2), gwr%uc_batch_size
     ndat = blocked_loop(ig2, wc_qbz%sizeb_local(2), gwr%uc_batch_size)

     !ABI_CHECK_IEQ(size(wc_qbz%buffer_cplx(:, ig2)), desc_q%npw, "npw")
     !ABI_CHECK_IEQ(size(rgp%buffer_cplx(:, ig2)), gwr%g_nfft * gwr%nspinor, "gwr%g_nfft * gwr%nspinor")

#if 0
     call fft_ug(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat, &
                 gwr%g_mgfft, gwr%g_ngfft, desc_q%istwfk, desc_q%gvec, desc_q%gbound, &
                 wc_qbz%buffer_cplx(:, ig2), &  ! in
                 rgp%buffer_cplx(:, ig2))       ! out
#else
     call uplan_q%execute_gr(ndat, wc_qbz%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))
#endif

     ! Multiply by e^{iq.r}
     if (.not. q_is_gamma) then
       !$OMP PARALLEL DO
       do idat=0,ndat-1
         rgp%buffer_cplx(:, ig2+idat) = ceiqr(:) * rgp%buffer_cplx(:, ig2+idat)
       end do
     end if
   end do ! ig2

   call uplan_q%free()

   ! MPI transposition: Wc(r,g') -> Wc(g',r)
   call rgp%ptrans("N", wc_gpr(my_iqf), free=.True.)
   end associate

   call wc_qbz%free()
 end do ! my_iqf

 ABI_FREE(ceiqr)

 call cwtime_report(" gwr_get_myq_wc_gpr:", cpu, wall, gflops)

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
 type(__slkmat_t),intent(inout) :: wc_rpr

!Local variables-------------------------------
!scalars
 integer :: iq_bz, ig2, npwsp, col_bsize, ir1, ndat !, idat, g0_q(3), my_iqf,
 !logical :: q_is_gamma
 character(len=500) :: msg
 type(desc_t) :: desc_qbz
 type(__slkmat_t) :: wc_ggp, rgp, gpr
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
 call gwr%rotate_wc(iq_bz, itau, spin, desc_qbz, wc_ggp)

 ! Allocate rgp PBLAS matrix to store Wc(r, g')
 npwsp = desc_qbz%npw * gwr%nspinor
 ABI_CHECK(block_dist_1d(npwsp, gwr%g_comm%nproc, col_bsize, msg), msg)
 call rgp%init(gwr%g_nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_qbz%istwfk, size_blocs=[-1, col_bsize])

 !ABI_CHECK_IEQ(size(wc_ggp%buffer_cplx, dim=2), size(rgp%buffer_cplx, dim=2), "len2")

 !call uplan_k%init(npw, nspinor, gwr%uc_batch_size, ngfft, istwfk, gvec, gwpc, gwr%dtset%use_gpu_cuda, gbound=desc_qbz%gbound)
 !call uplan_k%execute_gr(ndat, ug, ur, isign=+1, scale=.False.)
 !call uplan_k%execute_rg(ndat, ur, ug, isign=-1, scale=.True.)
 !call uplan_k%free()

 ! FFT. Results stored in rgp
 do ig2=1,wc_ggp%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ig2, wc_ggp%sizeb_local(2), gwr%uc_batch_size)
   call fft_ug(desc_qbz%npw, gwr%g_nfft, gwr%nspinor, ndat, &
               gwr%g_mgfft, gwr%g_ngfft, desc_qbz%istwfk, desc_qbz%gvec, desc_qbz%gbound, &
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
 call rgp%ptrans("N", gpr, free=.True.)

 ! Wc(g',r) -> Wc(r',r)
 do ir1=1,gpr%sizeb_local(2), gwr%uc_batch_size
   ndat = blocked_loop(ir1, gpr%sizeb_local(2), gwr%uc_batch_size)

   ! FIXME: FFT sign is wrong (should be - instead of + but I need to change the FFT API)
   ! Perform FFT Wc_q(g',r) -> Wc_q(r',r) and store results in wc_rgp.
   call fft_ug(desc_qbz%npw, gwr%g_nfft, gwr%nspinor, ndat, &
               gwr%g_mgfft, gwr%g_ngfft, desc_qbz%istwfk, desc_qbz%gvec, desc_qbz%gbound, &
               gpr%buffer_cplx(:, ir1), &     ! in
               wc_rpr%buffer_cplx(:, ir1))    ! out
 end do

 call gpr%free()
 call desc_qbz%free()
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
 integer,allocatable :: requests(:)
 real(dp), pointer :: weights_ptr(:,:)
 complex(dp) :: wgt_globmy(gwr%ntau, gwr%my_ntau)  ! Use complex instead of real to be able to call ZGEMM.
 complex(dp),allocatable :: cwork_myit(:,:,:), glob_cwork(:,:,:)
 type(__slkmat_t), pointer :: mats(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 sum_spins_ = .False.; if (present(sum_spins)) sum_spins_ = sum_spins

 call wrtout(std_out, sjoin(" Performing cosine transform. what:", what, ", mode:", mode))

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

 ! And now perform inhomogenous FT in parallel.
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
     ! TODO: Determine batch_size automatically to avoid going OOM
     !batch_size = 1
     !batch_size = 24
     batch_size = 48
     !batch_size = loc2_size

     ABI_MALLOC(cwork_myit, (gwr%my_ntau, loc1_size, batch_size))
     ABI_MALLOC(glob_cwork, (gwr%ntau, loc1_size, batch_size))
     ABI_MALLOC(requests, (batch_size))

     do ig2=1,mats(it0)%sizeb_local(2), batch_size
       ndat = blocked_loop(ig2, mats(it0)%sizeb_local(2), batch_size)

       ! Extract matrix elements as a function of tau.
       do idat=1,ndat
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           do ig1=1,mats(it0)%sizeb_local(1)
             cwork_myit(my_it, ig1, idat) = mats(itau)%buffer_cplx(ig1, ig2+idat-1)
           end do
         end do
       end do

       ! Compute contribution to itau matrix
       do idat=1,ndat
         !do ig1=1,mats(it0)%sizeb_local(1)
         ! do itau=1,gwr%ntau
         !   glob_cwork(itau, ig1, idat) = dot_product(wgt_globmy(itau, :), cwork_myit(:, ig1, idat))
         ! end do
         !end do
         call ZGEMM("N", "N", gwr%ntau, loc1_size, gwr%my_ntau, cone, &
                    wgt_globmy, gwr%ntau, cwork_myit(1,1,idat), gwr%my_ntau, czero, glob_cwork(1,1,idat), gwr%ntau)
         !call xmpi_isum_ip(glob_cwork(:,:,idat), gwr%tau_comm%value, requests(idat), ierr)
       end do

       !call xmpi_waitall_1d(requests(1:ndat), ierr)
       call xmpi_sum(glob_cwork, gwr%tau_comm%value, ierr)

       ! Update my local (g1, g2) entry to have it in imaginary-frequency.
       !!!$OMP PARALLEL DO PRIVATE(itau)
       do idat=1,ndat
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           do ig1=1,mats(it0)%sizeb_local(1)
             mats(itau)%buffer_cplx(ig1, ig2+idat-1) = glob_cwork(itau, ig1, idat)
           end do
         end do
       end do

     end do ! ig2

     ABI_FREE(cwork_myit)
     ABI_FREE(glob_cwork)
     ABI_FREE(requests)
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
              ! Spins are distributed thus we have to sum them.
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

subroutine desc_init(desc, kk, istwfk, ecut, gwr, kin_sorted)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc
 real(dp),intent(in) :: kk(3)
 integer,intent(in) :: istwfk
 real(dp),intent(in) :: ecut
 class(gwr_t),intent(in) :: gwr
 logical,optional,intent(in) :: kin_sorted

!Local variables-------------------------------
 integer :: ig

! *************************************************************************

 desc%kin_sorted = .False.; if (present(kin_sorted)) desc%kin_sorted = kin_sorted
 desc%istwfk = istwfk

 call get_kg(kk, desc%istwfk, ecut, gwr%cryst%gmet, desc%npw, desc%gvec, kin_sorted=desc%kin_sorted)

 ABI_MALLOC(desc%gbound, (2 * gwr%g_mgfft + 8, 2))
 call sphereboundary(desc%gbound, desc%istwfk, desc%gvec, gwr%g_mgfft, desc%npw)

 ! Find the index of g = 0.
 desc%ig0 = -1
 do ig=1,desc%npw
   if (all(desc%gvec(:,ig) == 0)) then
     desc%ig0 = ig; exit
   end if
 end do

end subroutine desc_init
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_get_vc_sqrt
!! NAME
!!  desc_get_vc_sqrt
!!
!! FUNCTION
!!  Compute square root of Coulomb interaction vc(q,g).
!!
!! SOURCE

subroutine desc_get_vc_sqrt(desc, qpt, q_is_gamma, gwr, comm)

!Arguments ------------------------------------
 class(desc_t),intent(inout) :: desc
 real(dp),intent(in) :: qpt(3)
 logical, intent(in) :: q_is_gamma
 class(gwr_t),intent(in) :: gwr
 integer,intent(in) :: comm

! *************************************************************************

 ABI_UNUSED([q_is_gamma])

 if (allocated(desc%vc_sqrt)) return
 ABI_MALLOC(desc%vc_sqrt, (desc%npw))

 call gwr%vcgen%get_vc_sqrt(qpt, desc%npw, desc%gvec, gwr%q0, gwr%cryst, desc%vc_sqrt, comm)

end subroutine desc_get_vc_sqrt
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/desc_copy
!! NAME
!!  desc_copy
!!
!! FUNCTION
!!  Copy object
!!  NB: cannot use obj1 = obj2 syntax because ABINIT memory-leak detector
!!  won't see the allocation automatically performed by the compiler.
!!
!! SOURCE

subroutine desc_copy(in_desc, new_desc)

!Arguments ------------------------------------
 class(desc_t),intent(in) :: in_desc
 class(desc_t),intent(out) :: new_desc

! *************************************************************************

 call new_desc%free()

 new_desc%istwfk = in_desc%istwfk
 new_desc%npw = in_desc%npw
 new_desc%ig0 = in_desc%ig0
 new_desc%kin_sorted = in_desc%kin_sorted

 call alloc_copy(in_desc%gvec, new_desc%gvec)
 call alloc_copy(in_desc%gbound, new_desc%gbound)
 if (allocated(in_desc%vc_sqrt)) call alloc_copy(in_desc%vc_sqrt, new_desc%vc_sqrt)

end subroutine desc_copy
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

!!****f* m_gwr/gwr_print
!! NAME
!!  gwr_print
!!
!! FUNCTION
!!  Print info on the gwr object.
!!
!! INPUTS
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
 !call ydoc%add_int("gw_icutcoul", gwr%dtset%gw_icutcoul)
 call ydoc%add_int("green_mpw", gwr%green_mpw)
 call ydoc%add_int("tchi_mpw", gwr%tchi_mpw)
 call ydoc%add_int1d("g_ngfft", gwr%g_ngfft(1:6))
 call ydoc%add_int1d("P gwr_np_kgts", gwr%dtset%gwr_np_kgts)
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
   if (mem_mb > zero) then
     call wrtout(std_out, sjoin("- Local memory for Green's functions: ", ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
   end if
 end if
 if (allocated(gwr%tchi_qibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%tchi_qibz))
   if (mem_mb > zero) then
     call wrtout(std_out, sjoin('- Local memory for tchi_qibz: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
   end if
 end if
 if (allocated(gwr%wc_qibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%wc_qibz))
   if (mem_mb > zero) then
     call wrtout(std_out, sjoin('- Local memory for wc_iqbz: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
   end if
 end if
 if (allocated(gwr%sigc_kibz)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%sigc_kibz))
   if (mem_mb > zero) then
     call wrtout(std_out, sjoin('- Local memory for sigc_kibz: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
   end if
 end if
 if (allocated(gwr%ugb)) then
   mem_mb = sum(slk_array_locmem_mb(gwr%ugb))
   if (mem_mb > zero) then
     call wrtout(std_out, sjoin('- Local memory for ugb wavefunctions: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))
   end if
 end if
 call wrtout(std_out, " ")

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
!! SOURCE

subroutine gwr_build_tchi(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer :: my_is, my_it, my_ikf, ig, my_ir, my_nr, npwsp, ncol_glob, col_bsize, my_iqi, ii !, my_iki ! my_iqf,
 integer :: idat, ndat, sc_nfft, spin, ik_bz, iq_ibz, ierr, ipm, itau, ig2 !, ig1
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all, cpu_ir, wall_ir, gflops_ir
 real(dp) :: tchi_rfact, mem_mb, local_max, max_abs_imag_chit !, spin_fact, sc_ucvol, ik_ibz,
 complex(gwpc) :: head_q
 complex(dp) :: chq(3), wng(3)
 logical :: q_is_gamma
 character(len=5000) :: msg
 type(desc_t),pointer :: desc_k, desc_q
 type(__slkmat_t) :: chi_rgp
!arrays
 integer :: sc_ngfft(18), gg(3)
 integer,allocatable :: green_scgvec(:,:), chi_scgvec(:,:)
 integer :: mask_qibz(gwr%nqibz)
 real(dp) :: kk_bz(3), kpq_bz(3), qq_ibz(3), tsec(2) !, qq_bz(3) ! kk_ibz(3),
 !logical :: need_gt_kibz(gwr%nkibz), got_gt_kibz(gwr%nkibz)
 complex(gwpc),allocatable :: gt_scbox(:), chit_scbox(:)
 complex(gwpc),allocatable :: low_wing_q(:), up_wing_q(:), cemiqr(:) !, sc_ceimkr(:,:) !, uc_ceikr(:)
 !type(__slkmat_t) :: gt_gpr(2, gwr%my_nkbz), chiq_gpr(gwr%my_nqibz), gk_rpr_pm(2)
 !type(desc_t), target :: desc_mykbz(gwr%my_nkbz)
 type(__slkmat_t),allocatable :: gt_gpr(:,:), chiq_gpr(:), gk_rpr_pm(:)
 type(desc_t), target, allocatable :: desc_mykbz(:)
 type(fftbox_plan3_t) :: plan_gp2rp, plan_rp2gp
 type(uplan_t) :: uplan_q

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call timab(1923, 1, tsec)

 ABI_CHECK(gwr%tchi_space == "none", sjoin("tchi_space: ", gwr%tchi_space, " != none"))
 gwr%tchi_space = "itau"

 ! Allocate tchi matrices
 mask_qibz = 0; mask_qibz(gwr%my_qibz_inds(:)) = 1
 call gwr%print_mem(unit=std_out)
 call gwr%malloc_free_mats(mask_qibz, "tchi", "malloc")

 max_abs_imag_chit = zero

 if (gwr%use_supercell_for_tchi) then
   ! ============================
   ! GWR algorithm with supercell
   ! ============================

   ! Set FFT mesh in the supercell
   sc_ngfft = gwr%g_ngfft
   sc_ngfft(1:3) = gwr%ngkpt * gwr%g_ngfft(1:3)
   sc_ngfft(4:6) = sc_ngfft(1:3)
   sc_nfft = product(sc_ngfft(1:3))
   !sck_ucvol = gwr%cryst%ucvol * product(gwr%ngkpt)
   !scq_ucvol = gwr%cryst%ucvol * product(gwr%ngqpt)

   call wrtout(std_out, " Building chi0 with FFTs in the supercell:", pre_newlines=2)
   call wrtout(std_out, sjoin(" gwr_np_kgts:", ltoa(gwr%dtset%gwr_np_kgts)))
   call wrtout(std_out, sjoin(" ngkpt:", ltoa(gwr%ngkpt), " ngqpt:", ltoa(gwr%ngqpt)))
   call wrtout(std_out, sjoin(" gwr_boxcutmin:", ftoa(gwr%dtset%gwr_boxcutmin)))
   call wrtout(std_out, sjoin(" sc_ngfft:", ltoa(sc_ngfft(1:8))))
   call wrtout(std_out, sjoin(" my_ntau:", itoa(gwr%my_ntau), "ntau:", itoa(gwr%ntau)))
   call wrtout(std_out, sjoin(" my_nkbz:", itoa(gwr%my_nkbz), "nkbz:", itoa(gwr%nkbz)))
   call wrtout(std_out, sjoin(" my_nkibz:", itoa(gwr%my_nkibz), "nkibz:", itoa(gwr%nkibz)))
   call wrtout(std_out, sjoin(" FFT uc_batch_size:", itoa(gwr%uc_batch_size), &
                              " FFT sc_batch_size:", itoa(gwr%sc_batch_size)), do_flush=.True.)

   ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
   ! Perhaps the safest approach would be to generate the plan on the fly.
   ndat = gwr%sc_batch_size
   ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat * 2))
   ABI_CALLOC(chit_scbox, (sc_nfft * gwr%nspinor * ndat))
   mem_mb = (sc_nfft * gwr%nspinor * ndat * 3 * 2 * dp) * b2Mb
   call wrtout(std_out, sjoin(" Memory for scbox arrays:", ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

   ! This for the version with my_nkbz FFTs in the unit cell
   !ABI_CALLOC(gt_ucbox, (gwr%g_nfft * gwr%nspinor * ndat, 2))

   ! (ngfft, ndat)
   call plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat * 2, gwr%dtset%use_gpu_cuda)
   call plan_rp2gp%from_ngfft(sc_ngfft, gwr%nspinor * ndat,     gwr%dtset%use_gpu_cuda)

   ! The g-vectors in the supercell for G and tchi.
   ABI_MALLOC(green_scgvec, (3, gwr%green_mpw))
   ABI_MALLOC(chi_scgvec, (3, gwr%tchi_mpw))

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
     call chiq_gpr(my_iqi)%init(npwsp, gwr%g_nfft * gwr%nspinor, gwr%g_slkproc, 1, size_blocs=[-1, col_bsize])
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

       mem_mb = sum(slk_array_locmem_mb(gt_gpr))
       call wrtout(std_out, sjoin(' Local memory for gt_gpr: ', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

       ! Loop over r in the unit cell that is now MPI-distributed inside g_comm.
       my_nr = gt_gpr(1, 1)%sizeb_local(2)

       do my_ir=1, my_nr, gwr%sc_batch_size
         ndat = blocked_loop(my_ir, my_nr, gwr%sc_batch_size)
         if (my_ir <= 6 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
           call cwtime(cpu_ir, wall_ir, gflops_ir, "start")
         end if

         ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
         ! Note that we need to take the union of (k, g') for k in the BZ.
         gt_scbox = zero

!if (.False.) then
if (.True.) then
         do my_ikf=1,gwr%my_nkbz
           ik_bz = gwr%my_kbz_inds(my_ikf)

           ! Compute k+g
           desc_k => desc_mykbz(my_ikf)
           gg = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt)
           do ig=1,desc_k%npw
             green_scgvec(:,ig) = gg + gwr%ngkpt * desc_k%gvec(:,ig)
           end do

           do ipm=1,2
             !ABI_CHECK_IEQ(size(gt_gpr(ipm, my_ikf)%buffer_cplx, dim=2), my_nr, "my_nr!")
             !ABI_CHECK_IEQ(size(gt_gpr(ipm, my_ikf)%buffer_cplx, dim=1), desc_k%npw, "desc_k!")
             ii = 1 + (ipm - 1) * sc_nfft * gwr%nspinor * ndat

             call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat, green_scgvec, &
                           gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                           gt_scbox(ii:))                                ! inout
                           !gt_scbox(:,ipm))                             ! inout
           end do
         end do ! my_ikf

         ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
         if (gwr%kpt_comm%nproc > 1) call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)

         !call xmpi_isum_ip(gt_scbox(:,ipm), gwr%kpt_comm%value, request_ipm(ipm), ierr)
         !call xmpi_wait_all(request_ipm(ipm), ierr)
         !call cwtime_report("G part", cpu, wall, gflops)

         ! G(G',r) --> G(R',r) = sum_{k,g'} e^{-i(k+g').R'} G_k(g',r)
         call plan_gp2rp%execute(gt_scbox, -1)
         gt_scbox = gt_scbox * sc_nfft

else
         ! This is just to check whether replacing a single FFT in the supercell with
         ! nkpt FFTs in the unit cell + communication is faster.
         do ipm=1,2
           do my_ikf=1,gwr%my_nkbz
             desc_k => desc_mykbz(my_ikf)

             ! FIXME The sign is wrong as it should be -1.
             !call fft_ug(desc_k%npw, gwr%g_nfft, gwr%nspinor, ndat, &
             !            gwr%g_mgfft, gwr%g_ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
             !            gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
             !            gt_ucbox(:,ipm))                              ! out

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

         ! The GG part is the hotspot
         ! Compute tchi(R',r) for this r. Note that results are real so one might use r2c
         !chit_scbox(:) = gt_scbox(:, 1) * conjg(gt_scbox(:, 2))
         chit_scbox(:) = gt_scbox(1:sc_nfft * gwr%nspinor * ndat) * conjg(gt_scbox(sc_nfft * gwr%nspinor * ndat + 1:))
         !max_abs_imag_chit = max(max_abs_imag_chit, maxval(abs(aimag(chit_scbox))))

         ! Back to tchi(G'=q+g',r) space immediately with isign + 1.
         call plan_rp2gp%execute(chit_scbox, +1)

         ! Extract tchi_q(g',r) on the ecuteps g-sphere from the FFT box in the supercell
         ! and save data in chiq_gpr PBLAS matrix. Only my q-points in the IBZ are considered.
         do my_iqi=1,gwr%my_nqibz
           iq_ibz = gwr%my_qibz_inds(my_iqi)
           qq_ibz = gwr%qibz(:, iq_ibz)
           desc_q => gwr%tchi_desc_qibz(iq_ibz)

           ! Compute q+g vectors.
           gg = nint(qq_ibz * gwr%ngqpt)
           do ig=1,desc_q%npw
             chi_scgvec(:,ig) = gg + gwr%ngqpt(:) * desc_q%gvec(:,ig)
           end do

           !ABI_CHECK_IEQ(size(chiq_gpr(my_iqi)%buffer_cplx, dim=1), desc_q%npw, "desc_q!")
           call box2gsph(sc_ngfft, desc_q%npw, gwr%nspinor * ndat, chi_scgvec, &
                         chit_scbox, &                            ! in
                         chiq_gpr(my_iqi)%buffer_cplx(:, my_ir))  ! out

           ! Alternatively, one can avoid the above FFT,
           ! use zero-padded to go from the supercell to the ecuteps g-sphere inside the my_iqi loop.
           ! This approach should play well with k-point parallelism.
           !call fft_ur(desc_q%npw, sc_nfft, gwr%nspinor, ndat, sc_mgfft, sc_ngfft, &
           !            istwfk1, desc_q%gvec, sc_gbound_q, &
           !            chit_scbox, &                            ! in
           !            chiq_gpr(my_iqi)%buffer_cplx(:, my_ir))  ! out

         end do ! my_iqi
         !call cwtime_report("chiq part", cpu, wall, gflops)

         if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
           write(msg,'(4x, 3(a,i0),a)')" tChi my_ir [", my_ir, "/", my_nr, "] (tot: ", gwr%g_nfft, ")"
           if (gwr%comm%me == 0) call cwtime_report(msg, cpu_ir, wall_ir, gflops_ir)
         end if
       end do ! my_ir

       ! Free descriptors and PBLAS matrices in kBZ.
       call desc_array_free(desc_mykbz)
       call slk_array_free(gt_gpr)

       ! Now we have tchi_q(g',r).
       ! For each IBZ q-point treated by this MPI proc, do:
       !
       !     1) MPI transpose to have tchi_q(r,g')
       !     2) FFT along the first dimension to get tchi_q(g,g') stored in gwr%tchi_qibz
       !
       tchi_rfact = one / gwr%g_nfft / gwr%cryst%ucvol / (gwr%nkbz * gwr%nqbz)
       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         q_is_gamma = normv(gwr%qibz(:,iq_ibz), gwr%cryst%gmet, "G") < GW_TOLQ0
         desc_q => gwr%tchi_desc_qibz(iq_ibz)

         ! Note minus sign in q.
         if (.not. q_is_gamma) call calc_ceikr(-gwr%qibz(:,iq_ibz), gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, cemiqr)

         ! MPI-transposition: tchi_q(g',r) => tchi_q(r,g')
         call chiq_gpr(my_iqi)%ptrans("N", chi_rgp)
         !ABI_CHECK_IEQ(size(gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx, dim=2), size(chi_rgp%buffer_cplx, dim=2), "len2")

         ! FFT r --> g along the first dimension: tchi_q(r,g') --> tchi_q(g,g').
         ! Results stored in gwr%tchi_qibz.
         call uplan_q%init(desc_q%npw, gwr%nspinor, gwr%uc_batch_size, gwr%g_ngfft, istwfk1, &
                         desc_q%gvec, gwpc, gwr%dtset%use_gpu_cuda)

         do ig2=1, chi_rgp%sizeb_local(2), gwr%uc_batch_size
           ndat = blocked_loop(ig2, chi_rgp%sizeb_local(2), gwr%uc_batch_size)

           if (.not. q_is_gamma) then
             !$OMP PARALLEL DO
             do idat=0,ndat-1
               chi_rgp%buffer_cplx(:, ig2 + idat) = cemiqr(:) * chi_rgp%buffer_cplx(:, ig2 + idat)
             end do
           end if

#if 0
           call fft_ur(desc_q%npw, gwr%g_nfft, gwr%nspinor, ndat, gwr%g_mgfft, gwr%g_ngfft, &
                       istwfk1, desc_q%gvec, desc_q%gbound, &
                       chi_rgp%buffer_cplx(:, ig2),         &                  ! ur(in)
                       gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2))  ! ug(out)
#else

           call uplan_q%execute_rg(ndat, chi_rgp%buffer_cplx(:, ig2), &
                                 gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2))
#endif

           !$OMP PARALLEL DO
           do idat=0,ndat-1
             gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2 + idat) = &
             gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2 + idat) * tchi_rfact
           end do
         end do ! ig2

         call uplan_q%free()
         call chi_rgp%free()

         call gwr%tchi_qibz(iq_ibz, itau, spin)%set_imag_diago_to_zero(local_max)
       end do ! my_iqi

       write(msg,'(3(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
       call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau, end_str=ch10)
     end do ! my_it
   end do ! my_is

   ABI_FREE(gt_scbox)
   ABI_FREE(chit_scbox)
   ABI_FREE(green_scgvec)
   ABI_FREE(chi_scgvec)
   ABI_FREE(cemiqr)
   ABI_FREE(gt_gpr)
   ABI_FREE(gk_rpr_pm)
   ABI_FREE(desc_mykbz)
   !ABI_SFREE(gt_ucbox)
   !ABI_SFREE(sc_ceimkr)

   call slk_array_free(chiq_gpr)
   ABI_FREE(chiq_gpr)

   call plan_gp2rp%free()
   call plan_rp2gp%free()

  else
    ! ===================================================================
    ! Mixed-space algorithm in the unit cell with convolutions in k-space
    ! ===================================================================

    do my_is=1,gwr%my_nspins
      spin = gwr%my_spins(my_is)
      do my_it=1,gwr%my_ntau
        call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
        itau = gwr%my_itaus(my_it)

        ! Redistribute G_k(g,g') in the IBZ so that each MPI proc can recostruct G_k in the BZ on the fly.
        !need_gt_kibz = ?
        !call gwr%distrib_gt_kibz(itau, spin, need_kibz, got_kibz, "communicate')

        do my_ikf=1,gwr%my_nkbz
          ik_bz = gwr%my_kbz_inds(my_ikf)
          kk_bz = gwr%kbz(:, ik_bz)

          ! Use symmetries to get G_k(g,g') from the IBZ, then G_k(g,g') -> G_k(r,r').
          call gwr%get_gk_rpr_pm(ik_bz, itau, spin, gk_rpr_pm)

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
            ! here I may need to take into account the umklapp
            !call gwr%get_gk_rpr_pm(ikq_bz, itau, spin, gkq_rpr_pm)
            !chiq_rpr%buffer_cplx = chiq_rpr%buffer_cplx + gk_rpr_pm(1)%buffer_cplx * conjg(gkq_rpr_pm(2)%buffer_cplx)
          end do ! iq_ibz
        end do ! my_ikf

        ! Deallocate extra G's
        !call gwr%distrib_gt_kibz(itau, spin, need_kibz, got_kibz, "free")

        call slk_array_free(gk_rpr_pm)
        !gwr%tchi_qibz(iq_ibz, itau, spin)%buffer_cplx =

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

 if (gwr%kpt_comm%me == 0) then
   ! ===================================================
   ! ==== Construct head and wings from the tensor =====
   ! ===================================================
   associate (desc_q0 => gwr%tchi_desc_qibz(1), mat_ts => gwr%tchi_qibz(1,:,:))
   ABI_CHECK_IEQ(desc_q0%ig0, 1, "ig0 should be 1")
   ABI_MALLOC(up_wing_q, (desc_q0%npw))
   ABI_MALLOC(low_wing_q, (desc_q0%npw))

   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       do ig=2,desc_q0%npw
         wng = gwr%chi0_uwing_myw(:,ig, my_it)
         up_wing_q(ig) = vdotw(gwr%q0, wng, gwr%cryst%gmet, "G")
         wng = gwr%chi0_lwing_myw(:,ig,my_it)
         low_wing_q(ig) = vdotw(gwr%q0, wng, gwr%cryst%gmet, "G")
       end do
       chq = matmul(gwr%chi0_head_myw(:,:,my_it), gwr%q0)
       head_q = vdotw(gwr%q0, chq, gwr%cryst%gmet, "G")

       call mat_ts(itau, spin)%set_head_and_wings(head_q, low_wing_q, up_wing_q)
     end do ! my_it
   end do ! my_is
   end associate
   ABI_FREE(up_wing_q)
   ABI_FREE(low_wing_q)
 end if

 ! Print trace of chi_q(i omega) matrices for testing purposes.
 if (gwr%dtset%prtvol > 0) call gwr%print_trace("tchi_qibz")

 ! Write file with chi0(i omega)
 if (gwr%dtset%prtsuscep > 0) call gwr%ncwrite_tchi_wc("tchi", trim(gwr%dtfil%filnam_ds(4))//'_TCHIM.nc')

 call cwtime_report(" gwr_build_tchi:", cpu_all, wall_all, gflops_all)
 call timab(1923, 2, tsec)

end subroutine gwr_build_tchi
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_distrib_gt_kibz
!! NAME
!!  gwr_distrib_gt_kibz
!!
!! FUNCTION
!!  If action == "communicate":
!!      Redistribute G_k for fixed (itau, spin) according to `need_kibz` table.
!!      Also, set got_kibz to 1 for each IBZ k-point that has been received.
!!  If action == "free":
!!      Use got_kibz to deallocate matrices.
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
 integer,intent(inout) :: got_kibz(gwr%nkibz)
 character(len=*),intent(in) :: action

!Local variables-------------------------------
 integer :: ik_ibz, ipm, ierr, lsize(2) ! col_bsize, npwsp,
 integer :: do_mpi_kibz(gwr%nkibz), sender_kibz(gwr%nkibz)
 real(dp) :: kk_ibz(3), cpu, wall, gflops
 complex(dp),allocatable :: cbuf_k(:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

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
   call gwr%malloc_free_mats(got_kibz, "green", "malloc")

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
   call gwr%malloc_free_mats(got_kibz, "green", "free")

 case default
   ABI_ERROR(sjoin("Invalid action:", action))
 end select

 call cwtime_report(" gwr_distrib_gt_kibz:", cpu, wall, gflops)

end subroutine gwr_distrib_gt_kibz
!!***

!!****f* m_gwr/gwr_distrib_mats_qibz
!! NAME
!!  gwr_distrib_mats_qibz
!!
!! FUNCTION
!!  If action == "communicate":
!!      Redistribute chi_q (wc_q) for fixed (itau, spin) according to `need_qibz` table.
!!      Also, set got_qibz to 1 for each IBZ q-point that has been received.
!!  If action == "free":
!!      Use got_qibz to deallocate matrices.
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
 integer,intent(inout) :: got_qibz(gwr%nqibz)
 character(len=*),intent(in) :: action

!Local variables-------------------------------
 integer :: iq_ibz, ierr, lsize(2)  ! col_bsize, npwsp,
 integer :: do_mpi_qibz(gwr%nqibz), sender_qibz(gwr%nqibz)
 real(dp) :: qq_ibz(3), cpu, wall, gflops
 complex(dp),allocatable :: cbuf_q(:,:)

! *************************************************************************

 ABI_CHECK(what == "tchi" .or. what == "wc", sjoin("Invalid what:", what))
 call cwtime(cpu, wall, gflops, "start")

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
   call gwr%malloc_free_mats(got_qibz, what, "malloc")

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
   call gwr%malloc_free_mats(got_qibz, what, "free")

 case default
   ABI_ERROR(sjoin("Invalid action:", action))
 end select

 call cwtime_report(" gwr_distrib_mats_qibz:", cpu, wall, gflops)

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
 integer :: my_is, spin, my_it, itau, iq_ibz, ierr, my_iqi, my_iki, ik_ibz, ipm
 character(len=5000) :: comment !, msg
 integer :: units(2)
 complex(dp),allocatable :: ctrace3(:,:,:), ctrace4(:,:,:,:)
 type(__slkmat_t),pointer :: mats(:,:,:)

! *************************************************************************

 ! The same q/k point in the IBZ might be available on different procs in kpt_comm
 ! thus we have to rescale the trace before summing the results in gwr%comm.

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

   call xmpi_sum_master(ctrace3, 0, gwr%kts_comm%value, ierr)

   if (gwr%comm%me == master) then
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

   call xmpi_sum_master(ctrace4, master, gwr%kts_comm%value, ierr)

   if (gwr%comm%me == master) then
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
 integer :: my_iqi, my_it, my_is, iq_ibz, spin, itau, iw, ierr
 integer :: il_g1, il_g2, ig1, ig2, iglob1, iglob2, ig0
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_q, wall_q, gflops_q
 logical :: q_is_gamma, free_tchi
 character(len=5000) :: msg
 complex(dpc) :: vcs_g1, vcs_g2
 type(__slkmat_t) :: em1
 type(yamldoc_t) :: ydoc
!arrays
 real(dp) :: qq_ibz(3), tsec(2)
 complex(dpc) :: em1_wq(gwr%ntau, gwr%nqibz), eps_wq(gwr%ntau, gwr%nqibz)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call timab(1924, 1, tsec)
 call wrtout([std_out, ab_out], " Building correlated screening Wc", pre_newlines=2)

 ABI_CHECK(gwr%tchi_space == "iomega", sjoin("tchi_space: ", gwr%tchi_space, " != iomega"))

 if (allocated(gwr%wc_qibz)) then
   call slk_array_free(gwr%wc_qibz)
   ABI_FREE(gwr%wc_qibz)
   gwr%wc_space = "none"
 end if

 ABI_CHECK(gwr%wc_space == "none", sjoin("wc_space: ", gwr%wc_space, " != none"))
 gwr%wc_space = "iomega"

 ! =======================================
 ! Allocate PBLAS arrays for wc_qibz(g,g')
 ! =======================================
 ! Note that we have already summed tchi over spin.
 ! Also, G=0 corresponds to iglob = 1 as only q-points in the IBZ are treated.
 ! This is not true for the other q-points in the full BZ as we may have a non-zero umklapp g0_q
 ABI_MALLOC(gwr%wc_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))

 free_tchi = .True.; if (free_tchi) gwr%tchi_space = "none"
 em1_wq = zero; eps_wq = zero

 ! If possible, use 2d rectangular grid of processors for diagonalization.
 !call slkproc_4diag%init(gwr%g_comm%value)

 do my_iqi=1,gwr%my_nqibz
   call cwtime(cpu_q, wall_q, gflops_q, "start")
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   q_is_gamma = normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0

   associate (desc_q => gwr%tchi_desc_qibz(iq_ibz))
   ig0 = desc_q%ig0

   ! The spin loop is needed so that procs in different pools can operate
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
             if (iglob1 == ig0 .and. iglob2 == ig0) then
               ! Store epsilon_{iw, iq_ibz}(0, 0)
               ! Rescale by np_qibz because we will MPI reduce this array.
               eps_wq(itau, iq_ibz) = wc%buffer_cplx(il_g1, il_g2) / gwr%np_qibz(iq_ibz)
             end if
           end if
         end do ! il_g1
       end do ! il_g2

       ! Invert symmetrized epsilon.
       ! NB: PZGETRF requires square block cyclic decomposition along the two axis
       ! hence we neeed to redistribute the data before calling invert.

       call wc%change_size_blocs(em1) ! processor=slkproc_4diag
       !call em1%invert()
       call em1%hpd_invert("U") ! TODO: Can call hpd_invert
       call wc%take_from(em1, free=.True.)  ! processor=wc%processor)

       !call wrtout(std_out, sjoin(" e-1 at q:", ktoa(qq_ibz), "i omega:", ftoa(gwr%iw_mesh(itau) * Ha_eV), "eV"))
       !call print_arr(wc%buffer_cplx, unit=std_out)

       ! Build Wc(q, iw) = e^{-1}_q(g,g',iw) - delta_{gg'} v_q(g,g') by removing bare vc
       do il_g2=1,wc%sizeb_local(2)
         iglob2 = wc%loc2gcol(il_g2)
         ig2 = mod(iglob2 - 1, desc_q%npw) + 1
         vcs_g2 = desc_q%vc_sqrt(ig2)
         do il_g1=1,wc%sizeb_local(1)
           iglob1 = wc%loc2grow(il_g1)
           ig1 = mod(iglob1 - 1, desc_q%npw) + 1
           vcs_g1 = desc_q%vc_sqrt(ig1)

           if (iglob1 == ig0 .and. iglob2 == ig0) then
             ! Store epsilon^{-1}_{iw, iq_ibz}(0, 0). Rescale by np_qibz because we will MPI reduce this array.
             em1_wq(itau, iq_ibz) = wc%buffer_cplx(il_g1, il_g2) / gwr%np_qibz(iq_ibz)
           end if

           ! Subtract exchange part.
           if (iglob1 == iglob2) wc%buffer_cplx(il_g1, il_g2) = wc%buffer_cplx(il_g1, il_g2) - one

           ! Handle divergence in Wc for q --> 0
           if (q_is_gamma .and. (iglob1 == ig0 .or. iglob2 == ig0)) then
             if (iglob1 == ig0 .and. iglob2 == ig0) then
               vcs_g1 = sqrt(gwr%vcgen%i_sz); vcs_g2 = sqrt(gwr%vcgen%i_sz)
             else if (iglob1 == ig0) then
               !vcs_g1 = (four_pi) ** (three/two) * q0sph ** 2 / two
               vcs_g1 = sqrt(gwr%vcgen%i_sz)
             else if (iglob2 == ig0) then
               !vcs_g2 = (four_pi) ** (three/two) * q0sph ** 2 / two
               vcs_g2 = sqrt(gwr%vcgen%i_sz)
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

 call xmpi_sum_master(em1_wq, master, gwr%kgt_comm%value, ierr)
 call xmpi_sum_master(eps_wq, master, gwr%kgt_comm%value, ierr)

 if (gwr%comm%me == master) then
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
 call timab(1924, 2, tsec)

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
 integer,parameter :: master = 0
 integer :: my_is, my_it, spin, ik_ibz, ikcalc_ibz, sc_nfft, sc_size, my_ir, my_nr, iw, idat, ndat, ii !sc_mgfft,
 integer :: my_iqf, iq_ibz, iq_bz, itau, ierr, ibc, bmin, bmax, band, nbc ! col_bsize, ib1, ib2,
 integer :: my_ikf, ipm, ik_bz, ig, ikcalc, uc_ir, ir, ncid, col_bsize, npwsp, ibeg, iend, nrsp ! my_iqi, sc_ir,
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all !, cpu, wall, gflops
 real(dp) :: mem_mb, cpu_ir, wall_ir, gflops_ir
 real(dp) :: max_abs_imag_wct, max_abs_re_wct, sck_ucvol, scq_ucvol  !, spin_fact
 logical :: k_is_gamma !, q_is_gamma
 character(len=500) :: msg
 type(desc_t), pointer :: desc_q, desc_k
 type(yamldoc_t) :: ydoc
!arrays
 integer :: sc_ngfft(18), need_qibz(gwr%nqibz),  got_qibz(gwr%nqibz), units(3), gg(3)
 integer,allocatable :: green_scgvec(:,:), wc_scgvec(:,:)
 real(dp) :: kk_bz(3), kcalc_bz(3), qq_bz(3), tsec(2)  !, qq_ibz(3)
 complex(gwpc) :: cpsi_r, sigc_pm(2)
 complex(dp) :: odd_t(gwr%ntau), even_t(gwr%ntau)
 complex(dp),target,allocatable :: sigc_it_diag_kcalc(:,:,:,:,:)
 complex(gwpc) ABI_ASYNC, allocatable :: gt_scbox(:), wct_scbox(:) !, gt_scbox(:,:),
 complex(gwpc),allocatable :: uc_psi_bk(:,:,:), scph1d_kcalc(:,:,:), uc_ceikr(:), ur(:)
 type(__slkmat_t) :: gt_gpr(2, gwr%my_nkbz), gk_rpr_pm(2), sigc_rpr(2,gwr%nkcalc), wc_rpr, wc_gpr(gwr%my_nqbz)
 type(desc_t), target :: desc_mykbz(gwr%my_nkbz), desc_myqbz(gwr%my_nqbz)
 type(fftbox_plan3_t) :: gt_plan_gp2rp, wt_plan_gp2rp
 integer :: band_val, ibv, ncerr, unt_it, unt_iw, unt_rw
 real(dp) :: e0, ks_gap, qp_gap, sigx, vxc_val, vu, v_meanf, eshift
 complex(dp) :: zz, zsc, sigc_e0, dsigc_de0, z_e0, sig_xc, hhartree_bk, qp_ene, qp_ene_prev
 integer,allocatable :: iperm(:)
 integer :: gt_request, wct_request
 real(dp),allocatable :: sorted_qpe(:)
 real(dp) :: e0_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 real(dp) :: spfunc_diag_kcalc(gwr%nwr, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 real(dp) :: rw_mesh(gwr%nwr)
 !real(dp) :: sigx_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 logical :: define
 integer :: qp_solver_ierr(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: ze0_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: sigc_iw_diag_kcalc(gwr%ntau, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: sigc_e0_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: qpe_zlin_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol), imag_zmesh(gwr%ntau)
 complex(dp) :: qpe_pade_kcalc(gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 complex(dp) :: sigxc_rw_diag_kcalc(gwr%nwr, gwr%max_nbcalc, gwr%nkcalc, gwr%nsppol)
 type(sigma_pade_t) :: spade

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call timab(1925, 1, tsec)

 ABI_CHECK(gwr%wc_space == "itau", sjoin("wc_space: ", gwr%wc_space, " != itau"))

 !mask_kibz = 0; mask_kibz(gwr%my_kibz_inds(:)) = 1
 !call gwr%malloc_free_mats(mask_kibz, "sigma" "malloc")

 !if (gwr%scf_iteration == 1) then
 !else
 !end if

 ! Set FFT mesh in the supercell.
 ! Be careful when using the FFT plan as ndat can change inside the loop if we start to block.
 ! Perhaps the safest approach would be to generate the plan on the fly.

 sc_ngfft = gwr%g_ngfft
 sc_size = product(gwr%ngkpt)
 sc_ngfft(1:3) = gwr%ngkpt * gwr%g_ngfft(1:3)
 sc_ngfft(4:6) = sc_ngfft(1:3)
 sc_nfft = product(sc_ngfft(1:3))
 !sc_mgfft = maxval(sc_ngfft(1:3))
 sck_ucvol = gwr%cryst%ucvol * product(gwr%ngkpt)
 scq_ucvol = gwr%cryst%ucvol * product(gwr%ngqpt)

 call wrtout(std_out, sjoin(" Building Sigma_c with FFTs in the supercell:", ltoa(sc_ngfft(1:3))), pre_newlines=2)
 call wrtout(std_out, sjoin(" gwr_np_kgts:", ltoa(gwr%dtset%gwr_np_kgts)))
 call wrtout(std_out, sjoin(" ngkpt:", ltoa(gwr%ngkpt), " ngqpt:", ltoa(gwr%ngqpt)))
 call wrtout(std_out, sjoin(" gwr_boxcutmin:", ftoa(gwr%dtset%gwr_boxcutmin)))
 call wrtout(std_out, sjoin(" sc_ngfft:", ltoa(sc_ngfft(1:8))))
 call wrtout(std_out, sjoin(" my_ntau:", itoa(gwr%my_ntau), "ntau:", itoa(gwr%ntau)))
 call wrtout(std_out, sjoin(" my_nkbz:", itoa(gwr%my_nkbz), "nkibz:", itoa(gwr%nkibz)))
 call wrtout(std_out, sjoin(" FFT uc_batch_size:", itoa(gwr%uc_batch_size), &
                            " FFT sc_batch_size:", itoa(gwr%sc_batch_size)), do_flush=.True.)

 ! The g-vectors in the supercell for G and tchi.
 ABI_MALLOC(green_scgvec, (3, gwr%green_mpw))
 ABI_MALLOC(wc_scgvec, (3, gwr%tchi_mpw))

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
 ! 2) Compute and store Sigma_c^k(g,g', i omega) and then compute the matrix elements in g-space.
 !
 ! The first option requires less memory provided we are interested in a small set of KS states.
 ! The second option is interesting if we need to compute several matrix elements, including off-diagonal terms.

 ndat = gwr%sc_batch_size
 ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat * 2))
 ABI_CALLOC(wct_scbox, (sc_nfft * gwr%nspinor * ndat))

 ! (ngfft, ndat)
 call gt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat * 2, gwr%dtset%use_gpu_cuda)
 call wt_plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat,     gwr%dtset%use_gpu_cuda)

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   ! Load wavefunctions for GW corrections in the unit cell.
   bmin = minval(gwr%bstart_ks(:, spin))
   bmax = maxval(gwr%bstop_ks(:, spin))

   ABI_MALLOC_OR_DIE(uc_psi_bk, (gwr%g_nfft * gwr%nspinor, bmin:bmax, gwr%nkcalc), ierr)
   ABI_MALLOC(ur, (gwr%g_nfft * gwr%nspinor))
   ABI_MALLOC(uc_ceikr, (gwr%g_nfft * gwr%nspinor))

   ! Precompute one-dimensional factors to get 3d e^{ik.L}
   call get_1d_sc_phases(gwr%ngkpt, gwr%nkcalc, gwr%kcalc, scph1d_kcalc)

   do ikcalc=1,gwr%nkcalc
     kcalc_bz = gwr%kcalc(:, ikcalc)
     ikcalc_ibz = gwr%kcalc2ibz(ikcalc, 1)  ! NB: Assuming wfs in the IBZ.

     ! Compute e^{ik.r} phases in the unit cell.
     call calc_ceikr(kcalc_bz, gwr%g_ngfft, gwr%g_nfft, gwr%nspinor, uc_ceikr)

     do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
       call gwr%kcalc_wfd%get_ur(band, ikcalc_ibz, spin, ur)
       uc_psi_bk(:, band, ikcalc) = ur * uc_ceikr
     end do
   end do ! ikcalc

   ABI_FREE(ur)
   ABI_FREE(uc_ceikr)

   ! Construct Sigma(itau) in the supercell.
   do my_it=1,gwr%my_ntau
     call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
     itau = gwr%my_itaus(my_it)

     ! G_k(g,g') --> G_k(g',r) e^{ik.r} for each k in the BZ treated by me.
     call gwr%get_myk_green_gpr(itau, spin, desc_mykbz, gt_gpr)
     mem_mb = sum(slk_array_locmem_mb(gt_gpr))
     call wrtout(std_out, sjoin(' Local memory for gt_gpr:', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

     ! Wc_q(g,g') --> Wc_q(g',r) e^{iq.r} for each q in the BZ treated by me.
     call gwr%get_myq_wc_gpr(itau, spin, desc_myqbz, wc_gpr)
     mem_mb = sum(slk_array_locmem_mb(wc_gpr))
     call wrtout(std_out, sjoin(' Local memory for wc_gpr:', ftoa(mem_mb, fmt="f8.1"), ' [Mb] <<< MEM'))

     my_nr = gt_gpr(1,1)%sizeb_local(2)
     ABI_CHECK(my_nr == wc_gpr(1)%sizeb_local(2), "my_nr != wc_gpr(1)%sizeb_local(2)")

     ! Loop over r in the unit cell that is now MPI-distributed in g_comm.
     do my_ir=1, my_nr, gwr%sc_batch_size
       if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
         call cwtime(cpu_ir, wall_ir, gflops_ir, "start")
       end if
       ndat = blocked_loop(my_ir, my_nr, gwr%sc_batch_size)

       uc_ir = gt_gpr(1,1)%loc2gcol(my_ir)  ! FIXME: This won't work if nspinor 2

       ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
       ! Take the union of (k,g') for k in the BZ.
       gt_scbox = zero
       do my_ikf=1,gwr%my_nkbz
         ik_bz = gwr%my_kbz_inds(my_ikf)
         ik_ibz = gwr%kbz2ibz(1, ik_bz)

         ! Compute k+g'
         desc_k => desc_mykbz(my_ikf)
         gg = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt)
         do ig=1,desc_k%npw
           green_scgvec(:,ig) = gg + gwr%ngkpt * desc_k%gvec(:,ig)
         end do

         do ipm=1,2
           ii = 1 + (ipm - 1) * sc_nfft * gwr%nspinor * ndat
           call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat, green_scgvec, &
                         gt_gpr(ipm, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                         gt_scbox(ii:))                                ! inout
                         !gt_scbox(:, ipm))                            ! inout
         end do
       end do ! my_ikf

       ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       if (gwr%kpt_comm%nproc > 1) call xmpi_isum_ip(gt_scbox, gwr%kpt_comm%value, gt_request, ierr)

       ! Insert Wc_q(g',r) in G'-space in the supercell FFT box for fixed r.
       ! Take the union of (q,g') for q in the BZ. Also, note ngqpt instead of ngkpt.
       wct_scbox = zero
       do my_iqf=1,gwr%my_nqbz
         iq_bz = gwr%my_qbz_inds(my_iqf)
         iq_ibz = gwr%qbz2ibz(1, iq_bz)

         ! Compute q+g' vectors.
         desc_q => desc_myqbz(my_iqf)
         gg = nint(gwr%qbz(:, iq_bz) * gwr%ngqpt)
         do ig=1,desc_q%npw
           wc_scgvec(:,ig) = gg + gwr%ngqpt * desc_q%gvec(:,ig)
         end do

         call gsph2box(sc_ngfft, desc_q%npw, gwr%nspinor * ndat, wc_scgvec, &
                       wc_gpr(my_iqf)%buffer_cplx(:, my_ir), & ! in
                       wct_scbox)                              ! inout
       end do ! my_iqf

       ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       if (gwr%kpt_comm%nproc > 1) call xmpi_isum_ip(wct_scbox, gwr%kpt_comm%value, wct_request, ierr)

       ! G(G',r) --> G(R',r)
       if (gwr%kpt_comm%nproc > 1) call xmpi_wait(gt_request, ierr)
       call gt_plan_gp2rp%execute(gt_scbox, -1)
       gt_scbox = gt_scbox * (sc_nfft / sck_ucvol)

       ! Wc(G',r) --> Wc(R',r)
       if (gwr%kpt_comm%nproc > 1) call xmpi_wait(wct_request, ierr)
       call wt_plan_gp2rp%execute(wct_scbox, -1)
       wct_scbox = wct_scbox * (sc_nfft / scq_ucvol)

       !DEBUG
       !if (itau == 1) print *, "W(r,r):", wct_scbox(uc_ir)
       !max_abs_re_wct = max(max_abs_re_wct, maxval(abs(real(wct_scbox))))
       !max_abs_imag_wct = max(max_abs_imag_wct, maxval(abs(aimag(wct_scbox))))
       !wct_scbox(uc_ir) = zero
       !wct_scbox = real(wct_scbox)
       !END DEBUG

       ! Use gt_scbox to store GW (R',r, +/- i tau) for this r.
       ii = sc_nfft * gwr%nspinor * ndat
       gt_scbox(1:ii) = gt_scbox(1:ii) * wct_scbox(:)
       gt_scbox(ii+1:) = gt_scbox(ii+1:) * wct_scbox(:)
       !print *, "Maxval abs imag G:", maxval(abs(aimag(gt_scbox)))

       ! Integrate self-energy matrix elements in the R-supercell for fixed set of ndat r-points.
       do ikcalc=1,gwr%nkcalc
         k_is_gamma = normv(gwr%kcalc(:,ikcalc), gwr%cryst%gmet, "G") < GW_TOLQ0
         ! FIXME: Temporary hack till I find a better MPI algo for k-points.
         if (gwr%kpt_comm%skip(ikcalc)) cycle

         do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
           ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
           do idat=1,ndat
             ir = uc_ir + idat - 1
             cpsi_r = conjg(uc_psi_bk(ir, band, ikcalc))
             do ipm=1,2
               ibeg = 1 + (idat - 1) * sc_nfft * gwr%nspinor + (ipm - 1) * sc_nfft * gwr%nspinor * ndat
               iend = ibeg + sc_nfft * gwr%nspinor - 1
               call sc_sum(gwr%ngkpt, gwr%g_ngfft, gwr%nspinor, scph1d_kcalc(:,:,ikcalc), k_is_gamma, &
                 cpsi_r, gt_scbox(ibeg:), uc_psi_bk(:, band, ikcalc), sigc_pm(ipm))
             end do

             sigc_it_diag_kcalc(:, itau, ibc, ikcalc, spin) = &
             sigc_it_diag_kcalc(:, itau, ibc, ikcalc, spin) + sigc_pm(:)
           end  do
          end do
       end do ! ikcalc
       !call cwtime_report("Sig Matrix part", cpu, wall, gflops)

       if (my_ir <= 3 * gwr%sc_batch_size .or. mod(my_ir, PRINT_MODR) == 0) then
         write(msg,'(4x,3(a,i0),a)')" Sigma_c my_ir [", my_ir, "/", my_nr, "] (tot: ", gwr%g_nfft, ")"
         if (gwr%comm%me == 0) call cwtime_report(msg, cpu_ir, wall_ir, gflops_ir)
       end if
     end do ! my_ir

     ! Free descriptors and PBLAS matrices in kBZ and qBZ.
     call desc_array_free(desc_mykbz)
     call slk_array_free(gt_gpr)
     call desc_array_free(desc_myqbz)
     call slk_array_free(wc_gpr)

     write(msg,'(3(a,i0),a)')" Sigma_c my_itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
     call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau, end_str=ch10)
   end do ! my_it

   ABI_FREE(scph1d_kcalc)
   ABI_FREE(uc_psi_bk)
 end do ! my_is

 !call wrtout(std_out, sjoin(" Maxval abs re W:", ftoa(max_abs_re_wct)))
 !call wrtout(std_out, sjoin(" Maxval abs imag W:", ftoa(max_abs_imag_wct)))
 ABI_FREE(gt_scbox)
 ABI_FREE(wct_scbox)

 call gt_plan_gp2rp%free()
 call wt_plan_gp2rp%free()

else
 call wrtout(std_out, sjoin(" Building Sigma_c with convolutions in q-space:", ltoa(sc_ngfft(1:3))), pre_newlines=2)

 ! Allocate workspace to store Wc(r' r) and Sigma_kcalc(r',r)
 nrsp = gwr%g_nfft * gwr%nspinor
 col_bsize = nrsp / gwr%g_comm%nproc; if (mod(nrsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
 call wc_rpr%init(nrsp, nrsp, gwr%g_slkproc, 1, size_blocs=[-1, col_bsize])

 do ipm=1,2
   call gk_rpr_pm(ipm)%init(nrsp, nrsp, gwr%g_slkproc, 1, size_blocs=[-1, col_bsize])
   do ikcalc=1,gwr%nkcalc
     call sigc_rpr(ipm, ikcalc)%init(nrsp, nrsp, gwr%g_slkproc, 1, size_blocs=[-1, col_bsize])
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
   bmin = minval(gwr%bstart_ks(:, spin)); bmax = maxval(gwr%bstop_ks(:, spin))
   ABI_MALLOC_OR_DIE(uc_psi_bk, (nrsp, bmin:bmax, gwr%nkcalc), ierr)
   ABI_MALLOC(ur, (nrsp))

   do ikcalc=1,gwr%nkcalc
     kcalc_bz = gwr%kcalc(:, ikcalc)
     ikcalc_ibz = gwr%kcalc2ibz(ikcalc, 1)  ! NB: Assuming wfs in IBZ

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

     ! Redistribute W_q(g,g') in the IBZ so that each MPI proc can reconstruct Wc_q in the BZ on the fly.
     need_qibz = 0
     call gwr%distrib_mats_qibz("wc", itau, spin, need_qibz, got_qibz, "communicate")
     call slk_array_set(sigc_rpr, czero)

     do my_ikf=1,gwr%my_nkbz
       ik_bz = gwr%my_kbz_inds(my_ikf)
       kk_bz = gwr%kbz(:, ik_bz)

       ! Use symmetries to get G_k from the IBZ. G_k(g,g') --> G_k(r,r')
       call gwr%get_gk_rpr_pm(ik_bz, itau, spin, gk_rpr_pm)

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

     write(msg,'(3(a,i0),a)')" Sigma_c my_itau [", my_it, "/", gwr%my_ntau, "] (tot: ", gwr%ntau, ")"
     call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
   end do ! my_it

   ABI_FREE(uc_psi_bk)
 end do ! my_is

 call wc_rpr%free()
 call slk_array_free(sigc_rpr)
 call slk_array_free(gk_rpr_pm)
end if

 ABI_FREE(green_scgvec)
 ABI_FREE(wc_scgvec)

 ! Store matrix elements of Sigma_c(it), separate even and odd part
 ! then use sine/cosine transform to get Sigma_c(i omega).
 ! Finally, perform analytic continuation with Pade' to go to the real-axis
 ! and compute QP corrections and spectral functions.
 ! All procs execute this part as it's very cheap.

 sigc_it_diag_kcalc = -sigc_it_diag_kcalc * (gwr%cryst%ucvol / gwr%g_nfft) ** 2
 call xmpi_sum(sigc_it_diag_kcalc, gwr%comm%value, ierr)

 imag_zmesh(:) = j_dpc * gwr%iw_mesh
 qp_solver_ierr = 0
 sigc_iw_diag_kcalc = zero

 ! Save previous QP bands in qp_ebands_prev (needed for self-consistency)
 ! In the loop below, we also update gwr%qp_ebands%eig with the QP results and recompute occ/fermie.
 gwr%qp_ebands_prev%eig = gwr%qp_ebands%eig
 gwr%qp_ebands_prev%occ = gwr%qp_ebands%occ

 do spin=1,gwr%nsppol
   do ikcalc=1,gwr%nkcalc
      ik_ibz = gwr%kcalc2ibz(ikcalc, 1)

      do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
        ibc = band - gwr%bstart_ks(ikcalc, spin) + 1

        ! NB: e0 is always set to the KS energy even in case of self-consistency.
        e0 = gwr%ks_ebands%eig(band, ik_ibz, spin)
        sigx = gwr%x_mat(band, band, ikcalc, spin)

        ! Note vxc[n_val] instead of vxc[n_val + n_nlcc] with the model core charge.
        vxc_val = gwr%ks_me%vxcval(band, band, ik_ibz, spin)
        vu = zero
        if (gwr%dtset%usepawu /= 0) vu = gwr%ks_me%vu(band, band, ik_ibz, spin)
        v_meanf = vxc_val + vu

        ! f(t) = E(t) + O(t) = (f(t) + f(-t)) / 2  + (f(t) - f(-t)) / 2
        associate (vals_pmt => sigc_it_diag_kcalc(:,:, ibc, ikcalc, spin))
        even_t = (vals_pmt(1,:) + vals_pmt(2,:)) / two
        odd_t = (vals_pmt(1,:) - vals_pmt(2,:)) / two
        sigc_iw_diag_kcalc(:, ibc, ikcalc, spin) = matmul(gwr%cosft_wt, even_t) + j_dpc * matmul(gwr%sinft_wt, odd_t)
        end associate

        zz = cmplx(e0, zero)
        call spade%init(gwr%ntau, imag_zmesh, sigc_iw_diag_kcalc(:, ibc, ikcalc, spin), branch_cut=">")

        ! Solve the QP equation with Newton-Rapson starting from e0
        call spade%qp_solve(e0, v_meanf, sigx, zz, zsc, msg, ierr)
        qpe_pade_kcalc(ibc, ikcalc, spin) = zsc
        qp_solver_ierr(ibc, ikcalc, spin) = ierr
        if (ierr /= 0) then
          ABI_WARNING(msg)
        end if

        call spade%eval(zz, sigc_e0, dzdval=dsigc_de0)
        ! Z = (1 - dSigma / domega(E0))^{-1}
        z_e0 = one / (one - dsigc_de0)

        ! Compute linearized QP solution and store results
        qp_ene = e0 + z_e0 * (sigc_e0 + sigx - v_meanf)
        qpe_zlin_kcalc(ibc, ikcalc, spin) = qp_ene
        e0_kcalc(ibc, ikcalc, spin) = e0
        sigc_e0_kcalc(ibc, ikcalc, spin) = sigc_e0
        ze0_kcalc(ibc, ikcalc, spin) = z_e0

        ! IMPORTANT: Here we update qp_ebands%eig with the new results.
        gwr%qp_ebands%eig(band, ik_ibz, spin) = real(qp_ene)

        ! Compute Spectral function using linear mesh **centered** around KS e0.
        rw_mesh = arth(e0 - gwr%wr_step * (gwr%nwr / 2), gwr%wr_step, gwr%nwr)
        hhartree_bk = gwr%ks_ebands%eig(band, ik_ibz, spin) - v_meanf
        do iw=1,gwr%nwr
          zz = rw_mesh(iw)
          call spade%eval(zz, sigc_e0)
          sig_xc = sigx + sigc_e0
          sigxc_rw_diag_kcalc(iw, ibc, ikcalc, spin) = sig_xc

          spfunc_diag_kcalc(iw, ibc, ikcalc, spin) = one / pi * abs(aimag(sigc_e0)) &
            /( (real(rw_mesh(iw) - hhartree_bk - sig_xc)) ** 2 + (aimag(sigc_e0)) ** 2) / Ha_eV

          !Sr%hhartree = hdft - KS_me%vxcval
          !spfunc_diag_kcalc(iw, ibc, ikcalc, spin) = &
          !  one / pi * abs(aimag(sigc_e0)) &
          !  /( (real(rw_mesh(iw) - Sr%hhartree(ib, ib, ik_ibz, spin) - sigx_xc)) ** 2 &
          !    +(aimag(sigc_e0)) ** 2) / Ha_eV
        end do ! iw

      end do ! band
   end do ! ikcalc
 end do ! spin

 if (gwr%nkcalc == gwr%nkibz) then

   ! Here I shift the bands that are not explictly included in the SCF calculation.
   ! using the correction evaluated at bstop_ks/bstart_ks to accelerate self-consistent calculations.
   do spin=1,gwr%nsppol
     do ikcalc=1,gwr%nkcalc
        ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
        band = gwr%bstop_ks(ikcalc, spin)
        if (band + 1 <= size(gwr%qp_ebands%eig, dim=1)) then
          eshift = gwr%qp_ebands%eig(band, ik_ibz, spin) - gwr%qp_ebands_prev%eig(band, ik_ibz, spin)
          call wrtout(std_out, sjoin( "Correcting bands >= ", itoa(band+1), "with eshift:", ftoa(eshift * Ha_meV), "(meV)"))
          gwr%qp_ebands%eig(band + 1:, ik_ibz, spin) = gwr%qp_ebands%eig(band + 1:, ik_ibz, spin) + eshift
        end if
        band = gwr%bstart_ks(ikcalc, spin)
        if (band > 1) then ! unlikely
          eshift = gwr%qp_ebands%eig(band, ik_ibz, spin) - gwr%qp_ebands_prev%eig(band, ik_ibz, spin)
          call wrtout(std_out, sjoin( "Correcting bands < ", itoa(band), "with eshift:", ftoa(eshift * Ha_meV), "(meV)"))
          gwr%qp_ebands%eig(:band - 1, ik_ibz, spin) = gwr%qp_ebands%eig(:band - 1, ik_ibz, spin) + eshift
        end if
     end do
   end do

   ! Recompute occupancies and set fermie to zero.
   ! FIXME: Possible problem here if the QP energies are not ordered!
   call ebands_update_occ(gwr%qp_ebands, gwr%dtset%spinmagntarget, prtvol=gwr%dtset%prtvol, fermie_to_zero=.True.)
 end if

 ! Master writes results to ab_out, std_out and GWR.nc
 if (gwr%comm%me == 0) then

   if (any(qp_solver_ierr /= 0)) then
     ! Write warning if QP solver failed.
     ierr = count(qp_solver_ierr /= 0)
     call wrtout([ab_out, std_out], sjoin("QP solver failed for:", itoa(ierr), "states"))
   end if

   call write_notations([std_out, ab_out])
   do spin=1,gwr%nsppol
     do ikcalc=1,gwr%nkcalc
       ik_ibz = gwr%kcalc2ibz(ikcalc, 1)

       ydoc = yamldoc_open('GWR_SelfEnergy_ee', width=11, real_fmt='(3f8.3)')
       call ydoc%add_real1d('kpoint', gwr%kcalc(:, ikcalc))
       call ydoc%add_int('spin', spin, int_fmt="(i1)")
       call ydoc%add_int('gwr_scf_iteration', gwr%scf_iteration)
       call ydoc%add_string('gwr_task', gwr%dtset%gwr_task)

       ! Compute Gaps using KS band indices.
       band_val = gwr%ks_vbik(ik_ibz, spin)
       nbc = gwr%bstop_ks(ikcalc, spin) - gwr%bstart_ks(ikcalc, spin) + 1

       if (band_val >= gwr%bstart_ks(ikcalc, spin) .and. band_val + 1 <= gwr%bstop_ks(ikcalc, spin)) then
         ibv = band_val - gwr%bstart_ks(ikcalc, spin) + 1
         ks_gap = gwr%ks_ebands%eig(band_val+1, ik_ibz, spin) - gwr%ks_ebands%eig(band_val, ik_ibz, spin)

         ! This to detect possible band inversion in the QP energies and compute qp_gaps accordingly.
         call sort_rvals(nbc, real(qpe_zlin_kcalc(:, ikcalc, spin)), iperm, sorted_qpe, tol=tol12)

         if (iperm(ibv) /= ibv .or. iperm(ibv + 1) /= ibv + 1) then
           call ydoc%add_int('QP_VBM_band', iperm(ibv) + gwr%bstart_ks(ikcalc, spin) - 1)
           call ydoc%add_int('QP_CBM_band', iperm(ibv+1) + gwr%bstart_ks(ikcalc, spin) - 1)
           qp_gap = sorted_qpe(ibv+1) - sorted_qpe(ibv)
         else
           call ydoc%add_int('QP_VBM_band', ibv + gwr%bstart_ks(ikcalc, spin) - 1)
           call ydoc%add_int('QP_CBM_band', ibv+1 + gwr%bstart_ks(ikcalc, spin) - 1)
           qp_gap = gwr%qp_ebands%eig(band_val+1, ik_ibz, spin) - gwr%qp_ebands%eig(band_val, ik_ibz, spin)
         end if
         ABI_FREE(iperm)
         ABI_FREE(sorted_qpe)

         call ydoc%add_real('KS_gap', ks_gap * Ha_eV)
         call ydoc%add_real('QP_gap', qp_gap * Ha_eV)
         call ydoc%add_real('Delta_QP_KS', (qp_gap - ks_gap) * Ha_eV)
       end if

       call ydoc%open_tabular('data') !, tag='SigmaeeData')
       write(msg, "(a5, *(a9))") "Band", "E0", "<VxcDFT>", "SigX", "SigC(E0)", "Z", "E-E0", "E-Eprev", "E", "Occ(E)"
       call ydoc%add_tabular_line(msg)

       do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
         ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
         e0 = gwr%ks_ebands%eig(band, ik_ibz, spin)
         qp_ene = gwr%qp_ebands%eig(band, ik_ibz, spin)
         qp_ene_prev = gwr%qp_ebands_prev%eig(band, ik_ibz, spin)

         write(msg,'(i5, *(f9.3))') &
           band, &                                                        ! Band
           e0 * Ha_eV, &                                                  ! E0
           real(gwr%ks_me%vxcval(band, band, ik_ibz, spin)) * Ha_eV, &    ! <VxcDFT>
           real(gwr%x_mat(band, band, ikcalc, spin)) * Ha_eV, &           ! SigX
           real(sigc_e0_kcalc(ibc, ikcalc, spin)) * Ha_eV, &              ! SigC(E0)
           real(ze0_kcalc(ibc, ikcalc, spin)), &                          ! Z
           !aimag(ze0_kcalc(ibc, ikcalc, spin)), &                        ! Z2
           (real(qp_ene - e0)) * Ha_eV, &                                 ! E-E0
           real(qp_ene - qp_ene_prev) * Ha_eV, &                          ! E-Eprev
           real(qp_ene) * Ha_eV, &                                        ! E
           gwr%qp_ebands%occ(band, ik_ibz, spin)                          ! Occ(E)
         call ydoc%add_tabular_line(msg)
       end do

       call ydoc%write_units_and_free([std_out, ab_out])
     end do ! ikcalc
   end do ! spin

   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_SIGC_IT'), msg, newunit=unt_it, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_it, "(a)")"# Diagonal elements of Sigma_c(i tau, +/-) in atomic units"
   write(unt_it, "(a)")"# tau Re/Im Sigma_c(+itau) Re/Im Sigma_c(-itau)"

   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_SIGXC_IW'), msg, newunit=unt_iw, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_iw, "(a)")"# Diagonal elements of Sigma_xc(i omega) in eV units"
   write(unt_iw, "(a)")"# omega Re/Im Sigma_c(i omega)"

   if (open_file(strcat(gwr%dtfil%filnam_ds(4), '_SIGXC_RW'), msg, newunit=unt_rw, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt_rw, "(a)")"# Diagonal elements of Sigma_xc(omega) in eV units and spectral function A(omega)"
   write(unt_rw, "(a)")"# omega Re/Im Sigma_xc(omega), A(omega)"

   units = [unt_it, unt_iw, unt_rw]
   call write_units(units, "# Fermi energy set to zero. Energies in eV")
   call write_units(units, sjoin("# nkcalc:", itoa(gwr%nkcalc), ", nsppol:", itoa(gwr%nsppol)))

   ! TODO: Improve file format. Add compatibility with gnuplot format for datasets?
   do spin=1,gwr%nsppol
     do ikcalc=1,gwr%nkcalc
       ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
       call write_units(units, sjoin("# kpt:", ktoa(gwr%kcalc(:, ikcalc)), "spin:", itoa(spin)))
       do band=gwr%bstart_ks(ikcalc, spin), gwr%bstop_ks(ikcalc, spin)
         ibc = band - gwr%bstart_ks(ikcalc, spin) + 1
         e0 = gwr%ks_ebands%eig(band, ik_ibz, spin)
         sigx = gwr%x_mat(band, band, ikcalc, spin)

         call write_units(units, sjoin("# band:", itoa(band), ", spin:", itoa(spin)))
         call write_units(units, sjoin("# sigx_ev:", ftoa(sigx * Ha_eV)))

         ! Write Sigma_c(i tau) and Sigma_c(i omega)
         do itau=1,gwr%ntau
           ! FIXME itau is not ordered
           write(unt_it, "(*(es16.8))") &
             gwr%tau_mesh(itau), &
             c2r(sigc_it_diag_kcalc(1, itau, ibc, ikcalc, spin)), &
             c2r(sigc_it_diag_kcalc(2, itau, ibc, ikcalc, spin))
           write(unt_iw, "(*(es16.8))") &
             gwr%iw_mesh(itau) * Ha_eV, &
             (c2r(sigc_iw_diag_kcalc(itau, ibc, ikcalc, spin) + sigx)) * Ha_eV
         end do

         ! Write Sigma_xc(omega) and A(omega)
         rw_mesh = arth(e0 - gwr%wr_step * (gwr%nwr / 2), gwr%wr_step, gwr%nwr) * Ha_eV
         do iw=1,gwr%nwr
           write(unt_rw, "(*(es16.8))") &
             rw_mesh(iw), &
             c2r(sigxc_rw_diag_kcalc(iw, ibc, ikcalc, spin)) * Ha_eV, &
             spfunc_diag_kcalc(iw, ibc, ikcalc, spin) / Ha_eV
         end do
       end do
     end do
   end do

   close(unt_it); close(unt_iw); close(unt_rw)

   ! ======================
   ! Add results to GWR.nc
   ! ======================
   NCF_CHECK(nctk_open_modify(ncid, gwr%gwrnc_path, xmpi_comm_self))

   ! Define arrays with results.
   define = .True.
   if (define) then
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t("e0_kcalc", "dp", "max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("ze0_kcalc", "dp", "two, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("qpe_zlin_kcalc", "dp", "two, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("qpe_pade_kcalc", "dp", "two, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("qp_solver_ierr", "int", "max_nbcalc, nkcalc, nsppol"), &
       ! TODO: Write exchange matrix?
       !nctkarr_t("sigx_kcalc", "dp", "max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("sigc_it_diag_kcalc", "dp", "two, two, ntau, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("sigc_iw_diag_kcalc", "dp", "two, ntau, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("sigxc_rw_diag_kcalc", "dp", "two, nwr, max_nbcalc, nkcalc, nsppol"), &
       nctkarr_t("spfunc_diag_kcalc", "dp", "nwr, max_nbcalc, nkcalc, nsppol") &
     ])
     NCF_CHECK(ncerr)
   end if

   ! Write data.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("e0_kcalc"), e0_kcalc))
   NCF_CHECK(nf90_put_var(ncid, vid("ze0_kcalc"), c2r(ze0_kcalc)))
   !NCF_CHECK(nf90_put_var(ncid, vid("sigx_kcalc"), sigx_kcalc))
   NCF_CHECK(nf90_put_var(ncid, vid("qpe_zlin_kcalc"), c2r(qpe_zlin_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("qpe_pade_kcalc"), c2r(qpe_pade_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("qp_solver_ierr"), qp_solver_ierr))
   NCF_CHECK(nf90_put_var(ncid, vid("sigc_it_diag_kcalc"), c2r(sigc_it_diag_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("sigc_iw_diag_kcalc"), c2r(sigc_iw_diag_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("sigxc_rw_diag_kcalc"), c2r(sigxc_rw_diag_kcalc)))
   NCF_CHECK(nf90_put_var(ncid, vid("spfunc_diag_kcalc"), spfunc_diag_kcalc))
   NCF_CHECK(nf90_close(ncid))
 end if ! master

 !call xmpi_barrier(gwr%comm%value)

 ABI_FREE(sigc_it_diag_kcalc)
 call cwtime_report(" gwr_build_sigmac:", cpu_all, wall_all, gflops_all)
 call timab(1925, 2, tsec)

contains
integer function vid(vname)
 character(len=*),intent(in) :: vname
 vid = nctk_idname(ncid, vname)
end function vid

end subroutine gwr_build_sigmac
!!***

!!****f* m_gwr/write_notations
!! NAME
!!  write_notations
!!
!! FUNCTION
!!  Write the meaning of the different columns.
!!
!! SOURCE

subroutine write_notations(units)
 integer,intent(in) :: units(:)
 integer :: ii, unt

 do ii=1,size(units)
   unt = units(ii)
   write(unt,"(a)")repeat("=", 80)
   write(unt,"(a)")" QP results (energies in eV)"
   write(unt,"(a)")" Notations:"
   write(unt,"(a)")"     E0: Kohn-Sham energy"
   write(unt,"(a)")"     <VxcDFT>: Matrix elements of Vxc[n_val] without non-linear core correction (if any)"
   write(unt,"(a)")"     SigX: Matrix elements of Sigma_x"
   write(unt,"(a)")"     SigC(E0): Matrix elements of Sigma_c at E0"
   write(unt,"(a)")"     Z: Renormalization factor"
   write(unt,"(a)")"     E-E0: Difference between the QP and the KS energy."
   write(unt,"(a)")"     E-Eprev: Difference between QP energy at iteration i and i-1"
   write(unt,"(a)")"     E: Quasi-particle energy"
   write(unt,"(a)")"     Occ(E): Occupancy of QP state"
   !write(unt,"(a)")"     SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part."
   !write(unt,"(a)")"     TAU(eKS): Lifetime in femtoseconds computed at the KS energy."
   write(unt,"(a)")" "
   write(unt,"(a)")" "
 end do
end subroutine write_notations
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
 class(__slkmat_t),intent(in) :: mat
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
!!  Compute the correlated part of the total energy within ACFDT.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_rpa_energy(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_is, my_iqi, my_it, itau, spin, iq_ibz, ii, ierr, ig, ncut, icut, mat_size
 integer :: il_g1, il_g2, ig1, ig2, npw_q, ig0
 logical :: q_is_gamma
 real(dp) :: weight, qq_ibz(3), estep, aa, bb, rmsq, ecut_soft, damp, tsec(2)
 complex(dpc) :: vcs_g1, vcs_g2
 type(desc_t),pointer :: desc_q
 !character(len=500) :: msg
!arrays
 type(__slkmat_t) :: chi_tmp, dummy_vec, chi_4diag
 type(processor_scalapack) :: proc_4diag
 real(gwp),allocatable :: eig(:)
 real(dp),allocatable :: kin_qg(:), ec_rpa(:), ec_mp2(:), ecut_chi(:)

! *************************************************************************

 call gwr%build_chi0_head_and_wings()
 call gwr%build_green(free_ugb=.True.)
 call gwr%build_tchi()

 ABI_CHECK(gwr%tchi_space == "iomega", sjoin("tchi_space:", gwr%tchi_space, "!= iomega"))
 call timab(1928, 1, tsec)

 ! See also calc_rpa_functional in m_screening_driver
 ! Compute RPA energy for ncut cutoff energies in order to extrapolate for ecuteps --> oo
 ncut = 5; estep = -gwr%dtset%ecuteps * 0.05_dp

 ABI_CALLOC(ec_rpa, (ncut))
 ABI_CALLOC(ec_mp2, (ncut))
 ABI_MALLOC(ecut_chi, (ncut))
 ecut_chi = arth(gwr%dtset%ecuteps + tol12, estep, ncut)

 ! Polarizability has been summed over spins inside build_tchi.
 ! The loop over spins is needed to parallelize the loop over my_iqi if nsppol == 2.
 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   if (gwr%spin_comm%nproc == 1 .and. spin == 2) cycle

   do my_iqi=1,gwr%my_nqibz
     if (gwr%spin_comm%skip(my_iqi)) cycle
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     qq_ibz = gwr%qibz(:, iq_ibz)
     q_is_gamma = normv(qq_ibz, gwr%cryst%gmet, "G") < GW_TOLQ0
     !if (q_is_gamma) then
     !  call wrtout([std_out, ab_out], "RPA: Ignoring q==0"); cycle
     !end if

     ! iq_ibz might be replicated inside gwr%kpt_comm.
     if (.not. gwr%itreat_iqibz(iq_ibz)) cycle

     desc_q => gwr%tchi_desc_qibz(iq_ibz)
     ABI_CHECK(desc_q%kin_sorted, "g-vectors are not sorted by |q+g|^2/2 !")
     npw_q = desc_q%npw; ig0 = desc_q%ig0

     ABI_MALLOC(kin_qg, (npw_q))
     do ig=1,npw_q
       kin_qg(ig) = half * normv(qq_ibz + desc_q%gvec(:,ig), gwr%cryst%gmet, "G") ** 2
     end do

     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       associate (tchi => gwr%tchi_qibz(iq_ibz, itau, spin))

       if (my_it == 1) then
         ! Allocate workspace. NB: npw_q is the total number of PWs for this q.
         call tchi%copy(chi_tmp)
         ABI_CHECK_IEQ(npw_q, tchi%sizeb_global(1), "npw_q")
         ABI_MALLOC(eig, (npw_q))
       end if

       do icut=1,ncut

         ! Damp Coulomb kernel in order to have smooth E(V).
         ! See also https://www.vasp.at/wiki/index.php/ENCUTGWSOFT
         ! and Harl's PhD thesis available at: https://utheses.univie.ac.at/detail/2259
         ecut_soft = 0.8_dp * ecut_chi(icut)

         ! TODO: Contribution due to head for q --> 0 is ignored.
         ! This is not optimal but consistent with calc_rpa_functional
         do il_g2=1,tchi%sizeb_local(2)
           !ig2 = mod(tchi%loc2gcol(il_g2) - 1, desc_q%npw) + 1
           ig2 = tchi%loc2gcol(il_g2)
           damp = one
           !if (kin_qg(ig2) > ecut_soft) then
           !  damp = sqrt(half * (one + cos(pi * (kin_qg(ig2) - ecut_soft) / (ecut_chi(icut) - ecut_soft))))
           !end if
           vcs_g2 = desc_q%vc_sqrt(ig2) * damp
           if (q_is_gamma .and. ig2 == ig0) vcs_g2 = zero

           do il_g1=1,tchi%sizeb_local(1)
             !ig1 = mod(tchi%loc2grow(il_g1) - 1, desc_q%npw) + 1
             ig1 = tchi%loc2grow(il_g1)
             damp = one
             !if (kin_qg(ig1) > ecut_soft) then
             !  damp = sqrt(half * (one + cos(pi * (kin_qg(ig1) - ecut_soft) / (ecut_chi(icut) - ecut_soft))))
             !end if
             vcs_g1 = desc_q%vc_sqrt(ig1) * damp
             if (q_is_gamma .and. ig1 == ig0) vcs_g1 = zero

             chi_tmp%buffer_cplx(il_g1, il_g2) = tchi%buffer_cplx(il_g1, il_g2) * vcs_g1 * vcs_g2
           end do
         end do

         ! Diagonalize sub-matrix and perform integration in imaginary frequency.
         ! Eq (6) in 10.1103/PhysRevB.81.115126
         ! NB: have to build chi_tmp inside loop over icut as matrix is destroyed by pzheev.
         mat_size = bisect(kin_qg, ecut_chi(icut))

         ! Change size block and, if possible, use 2D rectangular grid of processors for diagonalization
         call proc_4diag%init(chi_tmp%processor%comm)
         call chi_tmp%change_size_blocs(chi_4diag, processor=proc_4diag)
         call chi_4diag%heev("N", "U", dummy_vec, eig, mat_size=mat_size)
         call chi_4diag%free()
         call proc_4diag%free()

         ! TODO: ELPA
         !call compute_eigen_problem(processor, matrix, results, eigen, comm, istwf_k, nev)

         if (xmpi_comm_rank(chi_tmp%processor%comm) == 0) then
           weight = gwr%wtq(iq_ibz) * gwr%iw_wgs(itau) / two_pi
           do ii=1,mat_size
             ec_rpa(icut) = ec_rpa(icut) + weight * (log(one - eig(ii)) + eig(ii))
             ! second order Moeller Plesset.
             ec_mp2(icut) = ec_mp2(icut) - weight * eig(ii) ** 2 / two
             !if (eig(ii) > zero) then
             !  write(msg, "(a, es16.8)")"Positive eigenvalue:", eig(ii)
             !  ABI_ERROR(msg)
             !end if
           end do
         end if
       end do ! icut

       if (my_it == gwr%my_ntau) then
         ! Free workspace
         call chi_tmp%free()
         ABI_FREE(eig)
       end if
       end associate

     end do ! my_it

     ABI_FREE(kin_qg)
   end do ! my_iqi
 end do ! my_is

 ! Collect results on the master node.
 call xmpi_sum_master(ec_rpa, master, gwr%comm%value, ierr)
 call xmpi_sum_master(ec_mp2, master, gwr%comm%value, ierr)

 if (gwr%comm%me == master) then
   ! Print results to ab_out.
   ! TODO: Add metadata: nband, nqbz...
   write(ab_out, "(4a16)")"ecut_chi", "ecut_chi^(-3/2)", "RPA Ec (eV)", "RPA Ec (Ha)"
   do icut=ncut,1,-1
     write(ab_out, "(*(es16.8))") ecut_chi(icut), ecut_chi(icut) ** (-three/two), ec_rpa(icut) * Ha_eV, ec_rpa(icut)
   end do
   if (ncut > 1) then
     ! Add last line with extrapolated value.
     rmsq = linfit(ncut, ecut_chi(:) ** (-three/two), ec_rpa, aa, bb)
     write(ab_out, "(2a16,*(es16.8))") "oo", "0", bb * Ha_eV, bb
   end if
 end if

 ABI_FREE(ec_rpa)
 ABI_FREE(ec_mp2)
 ABI_FREE(ecut_chi)

 call timab(1928, 2, tsec)

end subroutine gwr_rpa_energy
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_run_g0w0
!! NAME
!!  gwr_run_g0w0
!!
!! FUNCTION
!!  Compute QP energies within the G0W0 approximation and minimax meshes along the imaginary axis.
!!
!! INPUTS
!!  [free_ugb]: True if array with empty KS states should freed as soon as possibile. Default: True
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_run_g0w0(gwr, free_ugb)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr
 logical,optional,intent(in) :: free_ugb

!Local variables-------------------------------
 logical :: free_ugb__

! *************************************************************************

 ! Use ugb wavefunctions and the Lehmann representation to compute head/wings and Sigma_x matrix elements.
 call gwr%build_chi0_head_and_wings()
 call gwr%build_sigxme()

 ! Now compute G(itau) from ugb and start the GWR algorithm.
 free_ugb__ = .True.; if (present(free_ugb)) free_ugb__ = free_ugb
 call gwr%build_green(free_ugb=free_ugb__)
 call gwr%build_tchi()
 call gwr%build_wc()
 call gwr%build_sigmac()

end subroutine gwr_run_g0w0
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_run_energy_scf
!! NAME
!!  gwr_run_energy_scf
!!
!! FUNCTION
!!  Compute QP energies within energy-only self-consistent GW approximation
!!  and minimax meshes along the imaginary axis.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_run_energy_scf(gwr)

!Arguments ------------------------------------
 class(gwr_t),intent(inout) :: gwr

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: gwr_nstep, units(2)
 logical :: converged
 character(len=500) :: msg

! *************************************************************************

 !  [free_ugb]: True if array with empty KS states should freed as soon as possibile. Default: True

 ! TODO:
 ! To implement restart capabilities we need to read scf_iteration, qp_ebands and gwr_task from GWR.nc
 ! build_sigmac should be responsible for writing checkpoint data with qp_ebands at each iteration.
 units = [std_out, ab_out]
 gwr_nstep = gwr%dtset%gwr_nstep

 if (gwr%nkcalc /= gwr%nkibz) then
   ABI_ERROR("For energy-only GW, one should include all k-points in the IBZ")
 end if

 select case (gwr%dtset%gwr_task)
 case ("EGEW")
   converged = .False.
   call wrtout(units, " Begin energy-only self-consistency in both G and W (EGEW)")
   do while (.not. converged .and. gwr%scf_iteration <= gwr_nstep)
     call gwr%run_g0w0(free_ugb=.False.)
     gwr%scf_iteration = gwr%scf_iteration + 1
     call gwr%check_scf_cycle(converged)
   end do

 case ("EGW0")
   call wrtout(units, " Begin energy-only self-consistency in G (EGW0)")
   call gwr%run_g0w0(free_ugb=.False.)
   converged = .False.
   do while (.not. converged .and. gwr%scf_iteration <= gwr_nstep)
     gwr%scf_iteration = gwr%scf_iteration + 1
     call gwr%build_green(free_ugb=.False.)
     call gwr%build_sigxme()  ! NB: This should not change in semiconductors
     call gwr%build_sigmac()
     call gwr%check_scf_cycle(converged)
   end do

 case ("G0EW")
   ! This is more difficult to implement as we need to store G0 and eG
   ! and then use G only for chi and not in Sigma
   call wrtout(units, " Begin energy-only self-consistency in W (G0EW)")
   ABI_ERROR("G0WE is not yet implemented")
   call gwr%run_g0w0(free_ugb=.False.)
   converged = .False.
   do while (.not. converged .and. gwr%scf_iteration <= gwr_nstep)
     gwr%scf_iteration = gwr%scf_iteration + 1
     !call gwr%build_green(free_ugb=.False.)
     call gwr%build_chi0_head_and_wings()
     call gwr%build_tchi()
     call gwr%build_wc()
     call gwr%build_sigmac()
     call gwr%check_scf_cycle(converged)
   end do

 case default
   ABI_ERROR(sjoin("Invalid gwr_task:", gwr%dtset%gwr_task))
 end select

 if (gwr%comm%me == master) then
   if (converged) then
     write(msg, "(1x,4a,i0,a,f8.3,a)") &
       trim(gwr%dtset%gwr_task), " self-consistent loop:", ch10, &
       " Convergence achieved at iteration: ", gwr%scf_iteration, &
       " with gwr_tolqpe: ",gwr%dtset%gwr_tolqpe * Ha_meV, " (meV)"
     call wrtout(units, msg)
   else
     write(msg, "(1x,4a,f8.3,3a,i0,a)") &
       trim(gwr%dtset%gwr_task), " self-consistent loop:", ch10, &
       " WARNING: Could not converge with gwr_tolqpe: ",gwr%dtset%gwr_tolqpe * Ha_meV, " (meV)", ch10, &
       " after: ", gwr_nstep, " steps"
     call wrtout(units, msg)
   end if
 end if

end subroutine gwr_run_energy_scf
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/check_scf_cyle
!! NAME
!!  check_scf_cycle
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_check_scf_cycle(gwr, converged)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 logical,intent(out) :: converged

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: spin, ikcalc, ik_ibz, band, ib, jb
 character(len=500) :: msg
 real(dp) :: max_adiff, adiff(gwr%qp_ebands%mband)
 integer :: units(2)

! *************************************************************************

 max_adiff = -one; converged = .True.; units = [std_out, ab_out]

 if (gwr%comm%me == master) then
   call wrtout(units, sjoin(" Checking for convergence at iteration:", itoa(gwr%scf_iteration)))
 end if

 associate (now => gwr%qp_ebands, prev => gwr%qp_ebands_prev)
 do spin=1,gwr%nsppol
   do ikcalc=1,gwr%nkcalc ! TODO: Should be spin dependent!
     ! Compute max abs difference between QP at iteration i and i-1.
     ik_ibz = gwr%kcalc2ibz(ikcalc, 1)
     ib = gwr%bstart_ks(ikcalc, spin); jb = gwr%bstop_ks(ikcalc, spin)
     adiff = zero; adiff(ib:jb) = abs(now%eig(ib:jb, ik_ibz, spin) - prev%eig(ib:jb, ik_ibz, spin))
     band = maxloc(adiff, dim=1)
     max_adiff = max(max_adiff, adiff(band))
     if (adiff(band) > gwr%dtset%gwr_tolqpe) converged = .False.
     if (gwr%comm%me == master) then
       ! Write info
       write(msg, "(a,i0,1x,2a,i0)") " For k-point: ", ik_ibz, trim(ktoa(now%kptns(:,ik_ibz))),", spin: ", spin
       call wrtout(units, msg)
       write(msg, "(4x,a,es12.5,a,i0)")" max(abs(E_i - E_{i-1})): ", adiff(band) * Ha_meV, " (meV) for band: ", band
       call wrtout(units, msg)
     end if
   end do
 end do
 end associate

 ! Just to make sure that all MPI procs agree on this!
 call xmpi_land(converged, gwr%comm%value)

 if (gwr%comm%me == master) then
   write(msg, "(a,i0,a)") "QP gaps at iteration: ",gwr%scf_iteration," (Fermi energy set to zero)"
   call ebands_print_gaps(gwr%qp_ebands, std_out, header=msg)
   call ebands_print_gaps(gwr%qp_ebands, ab_out, header=msg)
   if (.not. converged) then
     call wrtout(units," Not converged --> start new iteration ...")
   !else
   !  call wrtout(units, sjoin(" Convergence achieved at iteration", itoa(gwr%scf_iteration)))
   end if
 end if

end subroutine gwr_check_scf_cycle
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
 integer :: my_is, my_iqi, my_it, spin, iq_ibz, itau, npwtot_q, my_ncols, my_gcol_start, ncid, ncerr !, ierr
 real(dp) :: cpu, wall, gflops
!arrays
 real(dp), ABI_CONTIGUOUS pointer :: fptr(:,:,:)
 type(__slkmat_t), pointer :: mats(:)

! *************************************************************************

 ! Cannot reuse SCR.nc/SUSC.nc fileformat as:
 !  - hscr_new requires ep%
 !  - old file formats assume Gamma-centered g vectors.

 call cwtime(cpu, wall, gflops, "start")

 if (gwr%comm%me == master) then
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

 call xmpi_barrier(gwr%comm%value)

 ! Reopen the file in gwr%comm.
 NCF_CHECK(nctk_open_modify(ncid, filepath, gwr%comm%value))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)

     ! The same q-point in the IBZ might be stored on different pools.
     ! To avoid writing the same array multiple times, we use itreat_qibz
     ! to select the procs inside gwr%kpt_comm who are gonna write this iq_ibz q-point.
     if (.not. gwr%itreat_iqibz(iq_ibz)) cycle

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

       ! FIXME: Assuming PBLAS matrix distributed in contiguous blocks along the column index.
       ! This part must be changed if we use round robin distribution.
       my_ncols = mats(itau)%sizeb_local(2)
       my_gcol_start = mats(itau)%loc2gcol(1)

       ! FIXME: This is wrong if spc
       !call c_f_pointer(c_loc(mats(itau)%buffer_cplx), fptr, shape=[2, npwtot_q, my_ncols])
       ABI_MALLOC(fptr, (2, npwtot_q, my_ncols))
       fptr(1,:,:) = dble(mats(itau)%buffer_cplx)
       fptr(2,:,:) = aimag(mats(itau)%buffer_cplx)

       ncerr = nf90_put_var(ncid, vid("mats"), fptr, &
                            start=[1, 1, my_gcol_start, itau, iq_ibz, spin], &
                            count=[2, npwtot_q, my_ncols, 1, 1, 1])
                            !stride=[1, gwr%g_comm%nproc, 1, 1, 1])
       ABI_FREE(fptr)
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

!!****f* m_gwr/gsph2box_dpc
!! NAME
!! gsph2box_dpc
!!
!! FUNCTION
!! Insert cg_k array defined on the k-centered g-sphere with npw vectors inside the FFT box.
!! The main difference wrt to sphere is that cfft is not initialized to zero. See notes below.
!!
!! INPUTS
!! ngfft:
!!   n1,n2,n3=physical dimension of the FFT box
!!   n4,n5,n6=memory dimension of cfft
!! npw=number of G vectors in basis at this k point
!! ndat=number of items to process
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! cg(npw*ndat)= contains values for npw G vectors in basis sphere
!!
!! OUTPUT
!! cfft(n4,n5,n6*ndat) = array on FFT box filled with cg data
!!      Note that cfft is intent(inout) so that we can add contributions from different k-points.
!!
!! SOURCE

subroutine gsph2box_dpc(ngfft, npw, ndat, kg_k, cg, cfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft(6), npw, ndat
!arrays
 integer,intent(in) :: kg_k(3, npw)
 complex(dp),intent(in) :: cg(npw * ndat)
 complex(dp),target,intent(inout) :: cfft(ngfft(4)*ngfft(5)*ngfft(6)*ndat)

!Local variables-------------------------------
 integer :: n1, n2, n3, n4, n5, n6, i1, i2, i3, idat, ipw
 complex(dp),contiguous,pointer :: cfft_ptr(:,:,:,:)
 !character(len=500) :: msg

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)
 call c_f_pointer(c_loc(cfft), cfft_ptr, shape=[n4, n5, n6, ndat])

 ! Insert cg into cfft
!$OMP PARALLEL DO PRIVATE(i1, i2, i3) IF (ndat > 1)
 do idat=1,ndat
   do ipw=1,npw
     i1 = modulo(kg_k(1, ipw), n1) + 1
     i2 = modulo(kg_k(2, ipw), n2) + 1
     i3 = modulo(kg_k(3, ipw), n3) + 1

     !if (any(kg_k(:,ipw) > ngfft(1:3)/2) .or. any(kg_k(:,ipw) < -(ngfft(1:3)-1)/2) ) then
     !  write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ipw)," falls outside the FFT box. Increase boxcutmin (?)"
     !  ABI_ERROR(msg)
     !end if

     cfft_ptr(i1,i2,i3,idat) = cg(ipw+npw*(idat-1))
   end do
 end do

end subroutine gsph2box_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gsph2box_spc
!! NAME
!! gsph2box_spc
!!
!! FUNCTION
!! Insert cg_k array defined on the k-centered g-sphere with npw vectors inside the FFT box.
!! The main difference wrt to sphere is that cfft is not initialized to zero. See notes below.
!!
!! INPUTS
!! ngfft:
!!   n1,n2,n3=physical dimension of the FFT box
!!   n4,n5,n6=memory dimension of cfft
!! npw=number of G vectors in basis at this k point
!! ndat=number of items to process
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! cg(npw*ndat)= contains values for npw G vectors in basis sphere
!!
!! OUTPUT
!! cfft(n4,n5,n6*ndat) = array on FFT box filled with cg data
!!      Note that cfft is intent(inout) so that we can add contributions from different k-points.
!!
!! SOURCE

subroutine gsph2box_spc(ngfft, npw, ndat, kg_k, cg, cfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft(6), npw, ndat
!arrays
 integer,intent(in) :: kg_k(3, npw)
 complex(sp),intent(in) :: cg(npw * ndat)
 complex(sp),target,intent(inout) :: cfft(ngfft(4)*ngfft(5)*ngfft(6)*ndat)

!Local variables-------------------------------
 integer :: n1, n2, n3, n4, n5, n6, i1, i2, i3, idat, ipw
 complex(sp),contiguous,pointer :: cfft_ptr(:,:,:,:)
 !character(len=500) :: msg

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)
 call c_f_pointer(c_loc(cfft), cfft_ptr, shape=[n4, n5, n6, ndat])

 ! Insert cg into cfft
!$OMP PARALLEL DO PRIVATE(i1, i2, i3) IF (ndat > 1)
 do idat=1,ndat
   do ipw=1,npw
     i1 = modulo(kg_k(1, ipw), n1) + 1
     i2 = modulo(kg_k(2, ipw), n2) + 1
     i3 = modulo(kg_k(3, ipw), n3) + 1

     !if (any(kg_k(:,ipw) > ngfft(1:3)/2) .or. any(kg_k(:,ipw) < -(ngfft(1:3)-1)/2) ) then
     !  write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ipw)," falls outside the FFT box. Increase boxcutmin (?)"
     !  ABI_ERROR(msg)
     !end if

     cfft_ptr(i1,i2,i3,idat) = cg(ipw+npw*(idat-1))
   end do
 end do

end subroutine gsph2box_spc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/box2gsph_dpc
!! NAME
!! box2gsph_dpc
!!
!! FUNCTION
!! Extract cg_k array defined on the k-centered g-sphere with npw vectors from the FFT box.
!!
!! INPUTS
!! ngfft:
!!   n1,n2,n3=physical dimension of the FFT box
!!   n4,n5,n6=memory dimension of cfft
!! npw=number of G vectors in basis at this k point
!! ndat=number of items to process
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! cfft(n4,n5,n6, ndat) = array on FFT box
!!
!! OUTPUT
!! cg(npw*ndat)= contains values for npw G vectors in basis sphere
!!
!! SOURCE

subroutine box2gsph_dpc(ngfft, npw, ndat, kg_k, cfft, cg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft(6), npw, ndat
!arrays
 integer,intent(in) :: kg_k(3, npw)
 complex(dp),target,intent(in) :: cfft(ngfft(4)*ngfft(5)*ngfft(6)*ndat)
 complex(dp),intent(out) :: cg(npw*ndat)

!Local variables-------------------------------
 integer :: n1, n2, n3, n4, n5, n6, i1, i2, i3, idat, ipw, icg
 complex(dp),contiguous,pointer :: cfft_ptr(:,:,:,:)
 !character(len=500) :: msg

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)
 call c_f_pointer(c_loc(cfft), cfft_ptr, shape=[n4, n5, n6, ndat])

 ! Extract cg from cfft, ignoring components outside range of cg sphere
 !$OMP PARALLEL DO PRIVATE(i1, i2, i3, icg) IF (ndat > 1)
 do idat=1,ndat
   do ipw=1,npw
     i1 = modulo(kg_k(1, ipw), n1) + 1
     i2 = modulo(kg_k(2, ipw), n2) + 1
     i3 = modulo(kg_k(3, ipw), n3) + 1

     !if (any(kg_k(:,ipw) > ngfft(1:3)/2) .or. any(kg_k(:,ipw) < -(ngfft(1:3)-1)/2) ) then
     !  write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ipw)," falls outside the FFT box. Increase boxcutmin (?)"
     !  ABI_ERROR(msg)
     !end if

     icg = ipw + (idat - 1) * npw
     cg(icg) = cfft_ptr(i1, i2, i3, idat)
   end do
 end do

end subroutine box2gsph_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/box2gsph_spc
!! NAME
!! box2gsph_spc
!!
!! FUNCTION
!! Extract cg_k array defined on the k-centered g-sphere with npw vectors from the FFT box.
!!
!! INPUTS
!! ngfft:
!!   n1,n2,n3=physical dimension of the FFT box
!!   n4,n5,n6=memory dimension of cfft
!! npw=number of G vectors in basis at this k point
!! ndat=number of items to process
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! cfft(n4,n5,n6, ndat) = array on FFT box
!!
!! OUTPUT
!! cg(npw*ndat)= contains values for npw G vectors in basis sphere
!!
!! SOURCE

subroutine box2gsph_spc(ngfft, npw, ndat, kg_k, cfft, cg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft(6), npw, ndat
!arrays
 integer,intent(in) :: kg_k(3, npw)
 complex(sp),target,intent(in) :: cfft(ngfft(4)*ngfft(5)*ngfft(6)*ndat)
 complex(sp),intent(out) :: cg(npw*ndat)

!Local variables-------------------------------
 integer :: n1, n2, n3, n4, n5, n6, i1, i2, i3, idat, ipw, icg
 complex(sp),contiguous,pointer :: cfft_ptr(:,:,:,:)
 !character(len=500) :: msg

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)
 call c_f_pointer(c_loc(cfft), cfft_ptr, shape=[n4, n5, n6, ndat])

 ! Extract cg from cfft, ignoring components outside range of cg sphere
 !$OMP PARALLEL DO PRIVATE(i1, i2, i3, icg) IF (ndat > 1)
 do idat=1,ndat
   do ipw=1,npw
     i1 = modulo(kg_k(1, ipw), n1) + 1
     i2 = modulo(kg_k(2, ipw), n2) + 1
     i3 = modulo(kg_k(3, ipw), n3) + 1

     !if (any(kg_k(:,ipw) > ngfft(1:3)/2) .or. any(kg_k(:,ipw) < -(ngfft(1:3)-1)/2) ) then
     !  write(msg,'(a,3(i0,1x),a)')" The G-vector: ",kg_k(:, ipw)," falls outside the FFT box. Increase boxcutmin (?)"
     !  ABI_ERROR(msg)
     !end if

     icg = ipw + (idat - 1) * npw
     cg(icg) = cfft_ptr(i1, i2, i3, idat)
   end do
 end do

end subroutine box2gsph_spc
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
 integer,intent(in) :: sc_shape(3), uc_ngfft(18), uc_nfft, nspinor, ndat
 complex(dp),intent(in) :: alpha, ur(uc_nfft, nspinor, ndat)
 complex(dp),intent(in) :: &
   sc_eikr(uc_ngfft(1)*sc_shape(1), uc_ngfft(2)*sc_shape(2), uc_ngfft(3)*sc_shape(3), nspinor)
 complex(dp),intent(out) :: &
   sc_psi(uc_ngfft(1)*sc_shape(1), uc_ngfft(2)*sc_shape(2), uc_ngfft(3)*sc_shape(3), nspinor, ndat)

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
!!  kk(3)= k-point in reduced coordinates of the unit cell
!!  sc_shape= shape of the supercell.
!!  uc_ngfft(18)= information about the 3d FFT for the unit cell.
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

!!****f* m_gwr/gwr_build_chi0_head_and_wings
!! NAME
!!  gwr_build_chi0_head_and_wings
!!
!! FUNCTION
!!  Compute head and wings of chi0 on the minimax frequency grid.
!!
!! SOURCE

subroutine gwr_build_chi0_head_and_wings(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer,parameter :: two_poles = 2, one_pole = 1, gwcomp0 = 0, spmeth0 = 0
 integer :: nsppol, nspinor, ierr, my_is, spin, my_ikf, itau, my_it
 integer :: ik_bz, ik_ibz, isym_k, trev_k, g0_k(3)
 !integer :: iq_bz, iq_ibz, isym_q, trev_q, g0_q(3)
 integer :: nkpt_summed, use_umklp, band1, band2, band1_start, band1_stop
 integer :: ib, il_b1, il_b2, nb, block_size, ii, mband
 integer :: istwf_ki, npw_ki, nI, nJ, nomega, io, iq, nq, dim_rtwg !ig,
 integer :: npwe, u_nfft, u_mgfft, u_mpw
 logical :: isirr_k, use_tr, is_metallic
 real(dp) :: spin_fact, weight, deltaf_b1b2, deltaeGW_b1b2, gwr_boxcutmin_c, zcut, qlen, eig_nk
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_k, wall_k, gflops_k
 complex(dpc) :: deltaeKS_b1b2
 type(__slkmat_t),pointer :: ugb_kibz
 character(len=5000) :: msg
 type(crystal_t),pointer :: cryst
 type(dataset_type),pointer :: dtset
 type(ebands_t),pointer :: now_ebands
 type(littlegroup_t) :: ltg_q
 type(desc_t),pointer :: desc_ki
!arrays
 integer :: gmax(3), u_ngfft(18), work_ngfft(18) ! spinor_padx(2,4), g0(3),
 integer,contiguous, pointer :: kg_ki(:,:)
 integer,allocatable :: gvec_q0(:,:), gbound_q0(:,:), u_gbound(:,:)
 real(dp) :: kk_ibz(3), kk_bz(3), tsec(2)
 real(dp),contiguous, pointer :: qp_eig(:,:,:), qp_occ(:,:,:), ks_eig(:,:,:), cwave(:,:)
 real(dp),allocatable :: work(:,:,:,:), qdirs(:,:)
 logical :: gradk_not_done(gwr%nkibz)
 logical,allocatable :: bbp_mask(:,:)
 complex(dpc) :: chq(3) !, wng(3)
 !complex(dp),allocatable :: ug1_block(:,:)
 complex(gwpc) :: rhotwx(3, gwr%nspinor**2) !, new_rhotwx(3, gwr%nspinor**2)
 complex(gwpc),allocatable :: ug2(:), ur1_kibz(:), ur2_kibz(:), ur_prod(:), rhotwg(:), ug1_block(:,:), ug1(:)
 complex(dpc) :: green_w(gwr%ntau), omega(gwr%ntau)
 complex(dpc),allocatable :: chi0_lwing(:,:,:), chi0_uwing(:,:,:), chi0_head(:,:,:), head_qvals(:)
 real(dp), allocatable :: gh1c_block(:,:,:,:)
 type(vkbr_t),allocatable :: vkbr(:)
 type(gsphere_t) :: gsph
 !type(ddkop_t) :: ddkop
 !type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *************************************************************************

 call timab(1927, 1, tsec)
 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call wrtout([std_out, ab_out], sjoin(" Computing chi0 head and wings with inclvkb:", itoa(gwr%dtset%inclvkb)), &
             pre_newlines=1)

 nspinor = gwr%nspinor; nsppol = gwr%nsppol; dtset => gwr%dtset; cryst => gwr%cryst
 use_tr = gwr%dtset%awtr == 1; zcut = gwr%dtset%zcut ! well, it's not used in g0w0 when omega is complex.

 ! Use KS or QP energies depending on the iteration state.
 if (gwr%scf_iteration == 1) then
   call wrtout([std_out, ab_out], " Using KS orbitals and KS energies...", newlines=1, do_flush=.True.)
   qp_eig => gwr%ks_ebands%eig; qp_occ => gwr%ks_ebands%occ
   now_ebands => gwr%ks_ebands
 else
   call wrtout([std_out, ab_out], " Using KS orbitals and QP energies...", newlines=1, do_flush=.True.)
   qp_eig => gwr%qp_ebands%eig; qp_occ => gwr%qp_ebands%occ
   now_ebands => gwr%qp_ebands
 end if

 ks_eig => gwr%ks_ebands%eig
 mband = gwr%ks_ebands%mband

 is_metallic = ebands_has_metal_scheme(now_ebands)

 ! Setup weight (2 for spin unpolarized systems, 1 for polarized).
 ! spin_fact is used to normalize the occupation factors to one.
 ! Consider also the AFM case.
 select case (nsppol)
 case (1)
   weight = two / gwr%nkbz; spin_fact = half
   if (gwr%nspden == 2) then
     weight = one / gwr%nkbz; spin_fact = half
   end if
   if (nspinor == 2) then
     weight = one / gwr%nkbz; spin_fact = one
   end if
 case (2)
   weight = one / gwr%nkbz; spin_fact = one
 case default
   ABI_BUG(sjoin("Wrong nsppol:", itoa(nsppol)))
 end select

 ! TODO: Replace vkbr with ddk and factorize calls to DDK |bra>
 ABI_MALLOC(vkbr, (gwr%nkibz))
 gradk_not_done = .TRUE.

 ! TODO: Might become 1b
 ABI_MALLOC(bbp_mask, (mband, mband))

 ! =========================================
 ! Find FFT mesh and max number of g-vectors
 ! =========================================
 ! TODO: Can be decreased. Consider also fftgw
 gwr_boxcutmin_c = two
 !gwr_boxcutmin_c = one
 call gwr%get_u_ngfft(gwr_boxcutmin_c, u_ngfft, u_nfft, u_mgfft, u_mpw, gmax)

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
 gmax = 2 * gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 if (gwr%comm%me == 0) then
   call print_ngfft(u_ngfft, header="FFT mesh for chi0 head/wings computation", unit=std_out)
   call print_ngfft(u_ngfft, header="FFT mesh for chi0 head/wings computation", unit=ab_out)
 endif

 ! Need to broacast G-vectors at q = 0 if k/q-point parallelism is activated.
 if (gwr%kpt_comm%me == 0) then
   npwe = gwr%tchi_desc_qibz(1)%npw
   ABI_CHECK(gwr%tchi_desc_qibz(1)%kin_sorted, "g-vectors are not sorted by |q+g|^2/2 !")
 end if
 call xmpi_bcast(npwe, 0, gwr%kpt_comm%value, ierr)
 ABI_MALLOC(gvec_q0, (3, npwe))
 if (gwr%kpt_comm%me == 0) gvec_q0 = gwr%tchi_desc_qibz(1)%gvec
 call xmpi_bcast(gvec_q0, 0, gwr%kpt_comm%value, ierr)

 ! This is needed to call accumulate_head_wings_imagw
 call gsph%init(cryst, npwe, gvec_q0)

 ABI_MALLOC(gbound_q0, (2 * u_mgfft + 8, 2))
 call sphereboundary(gbound_q0, istwfk1, gvec_q0, u_mgfft, npwe)

 ! Init little group to find IBZ_q
 use_umklp = 0
 call ltg_q%init([zero, zero, zero], gwr%nkbz, gwr%kbz, cryst, use_umklp, npwe) !, gvec=gvec_kss)

 nkpt_summed = gwr%nkbz
 if (dtset%symchi /= 0) then
   nkpt_summed = ltg_q%nibz_ltg
   call ltg_q%print(std_out, Dtset%prtvol)
 end if
 !call wrtout(std_out, sjoin(' Calculation status: ', itoa(nkpt_summed), ' k-points to be completed'))

 ! ============================================
 ! === Begin big fat loop over transitions ====
 ! ============================================

 ! NB: One might reduce the number of bands as head and wings converge fast wrt nband and slow wrt k-mesh.
 ! Should introduce a tolerance on the frequency part computed at the first minimax frequency and
 ! compute max_nband from this.

 ! Loop on spin to calculate $\chi_{\up,\up} + \chi_{\down,\down}$
 ! TODO: Spinor
 nI = 1; nJ = 1; nomega = gwr%ntau
 omega(:) = j_dpc * gwr%iw_mesh(:)
 ABI_CALLOC(chi0_lwing, (npwe*nI, nomega, 3))
 ABI_CALLOC(chi0_uwing, (npwe*nJ, nomega, 3))
 ABI_CALLOC(chi0_head, (3, 3, nomega))

 ABI_MALLOC(u_gbound, (2 * u_mgfft + 8, 2))
 ABI_MALLOC(ur1_kibz, (u_nfft * nspinor))
 ABI_MALLOC(ur2_kibz, (u_nfft * nspinor))
 ABI_MALLOC(ur_prod, (u_nfft * nspinor))
 dim_rtwg = 1 !; if (nspinor==2) dim_rtwg=2 ! Can reduce size depending on Ep%nI and Ep%nj
 ABI_MALLOC(rhotwg, (npwe * dim_rtwg))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   ! TODO:
   !ddkop = ddkop_new(dtset, gwr%cryst, gwr%pawtab, gwr%psps, gwr%mpi_enreg, u_mpw, u_ngfft)

   ! Loop over my k-points in the BZ.
   do my_ikf=1,gwr%my_nkbz
     ik_bz = gwr%my_kbz_inds(my_ikf)
     kk_bz = gwr%kbz(:, ik_bz)

     if (dtset%symchi == 1 .and. ltg_q%ibzq(ik_bz) /= 1) CYCLE ! Only IBZ_q
     call cwtime(cpu_k, wall_k, gflops_k, "start")

     ! FIXME: Be careful with the symmetry conventions here!
     ! and the interplay between umklapp in q and FFT
     ! Also, the assembly_chi0 routines assume symrec and trev_k in [1, 2]
     ik_ibz = gwr%kbz2ibz_symrel(1, ik_bz); isym_k = gwr%kbz2ibz_symrel(2, ik_bz)
     trev_k = gwr%kbz2ibz_symrel(6, ik_bz); g0_k = gwr%kbz2ibz_symrel(3:5, ik_bz)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     kk_ibz = gwr%kibz(:, ik_ibz)

     ugb_kibz => gwr%ugb(ik_ibz, spin)
     desc_ki => gwr%green_desc_kibz(ik_ibz)
     npw_ki   =  desc_ki%npw
     istwf_ki =  desc_ki%istwfk
     kg_ki    => desc_ki%gvec
     !print *, "ik_ibz", ik_ibz, "npw_ki", npw_ki

     ABI_MALLOC(ug1, (npw_ki * nspinor))
     ABI_MALLOC(ug2, (npw_ki * nspinor))

     call sphereboundary(u_gbound, istwf_ki, kg_ki, u_mgfft, npw_ki)

     if (gwr%usepaw == 0 .and. dtset%inclvkb /= 0 .and. gradk_not_done(ik_ibz)) then
       ! Include term <n,k|[Vnl,iqr]|n"k>' for q -> 0.
       call vkbr_init(vkbr(ik_ibz), cryst, gwr%psps, dtset%inclvkb, istwf_ki, npw_ki, kk_ibz, kg_ki)
       gradk_not_done(ik_ibz) = .FALSE.
     end if

     !call ddkop%setup_spin_kpoint(gwr%dtset, gwr%cryst, gwr%psps, spin, kk_ibz, istwf_ki, npw_ki, kg_ki)

     call chi0_bbp_mask(ik_ibz, ik_ibz, spin, spin_fact, use_tr, &
                        gwcomp0, spmeth0, gwr%ugb_nband, mband, now_ebands, bbp_mask)
     !bbp_mask = .True.

     ! FIXME: This part should be tested with tau/g-para
     ! TODO:
     !  1) Logic to determine block_size from memory.
     !  2) Add support for symchi = 0
     !  3) Invert the loops

     block_size = min(48, gwr%ugb_nband)
     block_size = min(200, gwr%ugb_nband)
     !block_size = 1

     do band1_start=1, gwr%ugb_nband, block_size
       if (all(.not. bbp_mask(band1_start:, :))) then
         !print *, "exiting band1_start loop"
         exit
       end if
       !print *, "band1_start, gwr%ugb_nband, block_size", band1_start, gwr%ugb_nband, block_size
       nb = blocked_loop(band1_start, gwr%ugb_nband, block_size)
       band1_stop = band1_start + nb - 1

       ! Collect nb bands starting from band1_start on each proc.
       call ugb_kibz%collect_cplx(npw_ki * nspinor, nb, [1, band1_start], ug1_block)

       ! Dump algorithm based on xmpi_sum. To be replaced by an all_gather
       ! The advantage is that it works indipendently of the PBLAS distribution.
       !ABI_CALLOC(ug1_block, (npw_ki * nspinor, nb))
       !do il_b1=1, ugb_kibz%sizeb_local(2)
       !  band1 = ugb_kibz%loc2gcol(il_b1)
       !  ii = band1 - band1_start + 1
       !  if (band1 >= band1_start .and. band1 <= band1_stop) ug1_block(:, ii) = ugb_kibz%buffer_cplx(:, il_b1)
       !end do ! il_b1
       !call xmpi_sum(ug1_block, ugb_kibz%processor%comm, ierr)

       ABI_MALLOC(gh1c_block, (2, npw_ki*nspinor, 3, nb))
       do il_b1=1, ugb_kibz%sizeb_local(2)
         band1 = ugb_kibz%loc2gcol(il_b1)
         eig_nk = gwr%ks_ebands%eig(band1, ik_ibz, spin)

         ! FIXME: This is wrong if spc
         call c_f_pointer(c_loc(ugb_kibz%buffer_cplx(:,il_b1)), cwave, shape=[2, npw_ki*nspinor])
         !call ddkop%apply(eig_nk, npw_ki, nspinor, cwave, cwaveprj)
         !gh1c_block(:,:,:,xx_ib) = ddkop%gh1c(:, 1:npw_ki*nspinor,:)
       end do

       ! Loop over "conduction" states.
       !do band1=band1_start, band1_stop
       do ib=1,nb
         band1 = band1_start + ib - 1
         ug1 = ug1_block(:, ib)
         call fft_ug(npw_ki, u_nfft, nspinor, ndat1, u_mgfft, u_ngfft, istwf_ki, kg_ki, u_gbound, ug1, ur1_kibz)
         !call fft_ug(npw_ki, u_nfft, nspinor, ndat1, u_mgfft, u_ngfft, istwf_ki, kg_ki, u_gbound, ug1_block(:,ib), ur1_kibz)

         ! Loop over "valence" states.
         !do band2=1,gwr%ugb_nband
         do il_b2=1, ugb_kibz%sizeb_local(2)
           band2 = ugb_kibz%loc2gcol(il_b2)

           deltaeKS_b1b2 = ks_eig(band1, ik_ibz, spin) - ks_eig(band2, ik_ibz, spin)
           deltaf_b1b2  = spin_fact * (qp_occ(band1, ik_ibz, spin) - qp_occ(band2, ik_ibz, spin))
           deltaeGW_b1b2 = qp_eig(band1, ik_ibz, spin) - qp_eig(band2, ik_ibz, spin)

           ! Skip negligible transitions.
           if (abs(deltaf_b1b2) < GW_TOL_DOCC) CYCLE
           ! Adler-Wiser expression.
           ! Add small imaginary of the Time-Ordered response function but only for non-zero real omega
           ! FIXME What about metals?
           if (.not. use_tr) then
             ! Adler-Wiser without time-reversal.
             do io=1,nomega
               green_w(io) = g0g0w(omega(io), deltaf_b1b2, deltaeGW_b1b2, zcut, GW_TOL_W0, one_pole)
             end do

           else
             if (band1 < band2) CYCLE ! Here we GAIN a factor ~2

             do io=1,nomega
               ! Rangel: In metals, the intra-band transitions term does not contain the antiresonant part
               ! if(abs(deltaeGW_b1b2)>GW_TOL_W0) green_w(io) = g0g0w(omega(io),deltaf_b1b2,deltaeGW_b1b2,zcut,GW_TOL_W0)
               if (band1 == band2) green_w(io) = g0g0w(omega(io), deltaf_b1b2, deltaeGW_b1b2, zcut, GW_TOL_W0, one_pole)
               if (band1 /= band2) green_w(io) = g0g0w(omega(io), deltaf_b1b2, deltaeGW_b1b2, zcut, GW_TOL_W0, two_poles)
             end do
           end if

           ug2 = ugb_kibz%buffer_cplx(:, il_b2)
           call fft_ug(npw_ki, u_nfft, nspinor, ndat1, u_mgfft, u_ngfft, istwf_ki, kg_ki, u_gbound, ug2, ur2_kibz)

           ! FIXME: nspinor 2 is wrong as we have a 2x2 matrix
           ur_prod(:) = conjg(ur1_kibz(:)) * ur2_kibz
           call fft_ur(npwe, u_nfft, nspinor, ndat1, u_mgfft, u_ngfft, istwfk1, gvec_q0, gbound_q0, ur_prod, rhotwg)

           if (gwr%usepaw == 0) then
             ! Matrix elements of i[H,r] for NC pseudopotentials.
             ! NB ug1 and ug2 are kind=gwpc
             rhotwx = nc_ihr_comm(vkbr(ik_ibz), cryst, gwr%psps, npw_ki, nspinor, istwf_ki, gwr%dtset%inclvkb, &
                                  kk_ibz, ug1, ug2, kg_ki)
           end if

           ! Treat a possible degeneracy between v and c.
           ! Adler-Wiser expression, to be consistent here we use the KS eigenvalues (?)
           if (abs(deltaeKS_b1b2) > GW_TOL_W0) then
             rhotwx = -rhotwx / deltaeKS_b1b2
           else
             rhotwx = czero_gw
           end if

           !new_rhotwx = zero
           !gh1c_block(:,:,:,ib) = ddkop%gh1c(:, 1:npw_ki*nspinor,:)

           ! NB: Using symrec conventions here
           ik_ibz = gwr%kbz2ibz(1, ik_bz); isym_k = gwr%kbz2ibz(2, ik_bz)
           trev_k = gwr%kbz2ibz(6, ik_bz); g0_k = gwr%kbz2ibz(3:5, ik_bz)
           trev_k = trev_k + 1  ! NB: GW routines assume trev in [1, 2]

           ! TODO: Metals
           call accumulate_head_wings_imagw( &
                                        npwe, nomega, nI, nJ, dtset%symchi, &
                                        is_metallic, ik_bz, isym_k, trev_k, nspinor, cryst, ltg_q, gsph, &
                                        rhotwx, rhotwg, green_w, chi0_head, chi0_lwing, chi0_uwing)
         end do ! band2
       end do ! band1

       ABI_FREE(ug1_block)
       ABI_SFREE(gh1c_block)
     end do ! band1_start

     !if (gwr%usepaw == 0 .and. dtset%inclvkb /= 0 .and. dtset%symchi == 1) then
     !  call vkbr_free(vkbr(ik_ibz)) ! Not need anymore as we loop only over IBZ.
     !end if

     ABI_FREE(ug1)
     ABI_FREE(ug2)
     write(msg,'(3(a,i0),a)')" my_ikf [", my_ikf, "/", gwr%my_nkbz, "] (tot: ", gwr%nkbz, ")"
     call cwtime_report(msg, cpu_k, wall_k, gflops_k)
   end do ! my_ikf

   !call ddkop%free()
 end do ! my_is

 ABI_FREE(bbp_mask)
 ABI_FREE(gvec_q0)
 ABI_FREE(gbound_q0)
 ABI_FREE(work)
 ABI_FREE(ur1_kibz)
 ABI_FREE(ur2_kibz)
 ABI_FREE(ur_prod)
 ABI_FREE(rhotwg)
 ABI_FREE(u_gbound)

 call vkbr_free(vkbr)
 ABI_FREE(vkbr)

 ! Collect head and wings.
 call xmpi_sum(chi0_head, gwr%comm%value, ierr)
 call xmpi_sum(chi0_lwing, gwr%comm%value, ierr)
 call xmpi_sum(chi0_uwing, gwr%comm%value, ierr)

 chi0_head = chi0_head * weight / cryst%ucvol
 ! Tensor in terms of reciprocal lattice vectors.
 do io=1,nomega
   chi0_head(:,:,io) = matmul(chi0_head(:,:,io), cryst%gmet) * (two_pi**2)
 end do
 chi0_lwing = chi0_lwing * weight / cryst%ucvol
 chi0_uwing = chi0_uwing * weight / cryst%ucvol

 ! ===============================================
 ! ==== Symmetrize chi0 in case of AFM system ====
 ! ===============================================
 ! Reconstruct $chi0{\down,\down}$ from $chi0{\up,\up}$.
 ! Works only in the case of magnetic group Shubnikov type IV.
 if (cryst%use_antiferro) then
   call symmetrize_afm_chi0(Cryst, gsph, ltg_q, npwe, nomega, &
     chi0_head=chi0_head, chi0_lwing=chi0_lwing, chi0_uwing=chi0_uwing)
 end if

 if (gwr%comm%me == 0 .and. gwr%dtset%prtvol >= 1) then
   ! Output results to ab_out and std_out
   ! ===================================================
   ! ==== Construct head and wings from the tensor =====
   ! ===================================================
   qlen = tol3
   call cryst%get_redcart_qdirs(nq, qdirs, qlen=qlen)
   ABI_MALLOC(head_qvals, (nq))
   call wrtout([std_out, ab_out], " Head of the irreducible polarizability for q --> 0", pre_newlines=1)
   call wrtout([std_out, ab_out], sjoin(" q0_len:", ftoa(qlen), "(Bohr^-1)"))
   write(msg, "(*(a14))") "iomega (eV)", "[100]", "[010]", "[001]", "x", "y", "z"
   call wrtout([std_out, ab_out], msg)
   do io=1,nomega
     do iq=1,nq
       chq = matmul(chi0_head(:,:,io), qdirs(:,iq))
       head_qvals(iq) = vdotw(qdirs(:, iq), chq, cryst%gmet, "G")
     end do
     write(msg, "(*(es12.5,2x))") gwr%iw_mesh(io) * Ha_eV, real(head_qvals(:))
     call wrtout([std_out, ab_out], msg)
     ! Write imag part to std_out
     write(msg, "(*(es12.5,2x))") gwr%iw_mesh(io) * Ha_eV, aimag(head_qvals(:))
     call wrtout(std_out, msg)
   end do
   call wrtout([std_out, ab_out], " ")
   ABI_FREE(qdirs)
   ABI_FREE(head_qvals)
 end if

 ! Save quantities for later use as this routine must be called before build_tchi.
 if (gwr%kpt_comm%me == 0) then
   ABI_REMALLOC(gwr%chi0_head_myw, (3, 3, gwr%my_ntau))
   ABI_REMALLOC(gwr%chi0_uwing_myw, (3, npwe, gwr%my_ntau))
   ABI_REMALLOC(gwr%chi0_lwing_myw, (3, npwe, gwr%my_ntau))

   do my_it=1,gwr%my_ntau
     itau = gwr%my_itaus(my_it)
     gwr%chi0_head_myw(:,:,my_it) = chi0_head(:,:,itau)
     do ii=1,3
       gwr%chi0_uwing_myw(ii,:,my_it) = chi0_uwing(:,itau,ii)
       gwr%chi0_lwing_myw(ii,:,my_it) = chi0_lwing(:,itau,ii)
     end do
   end do
 end if

 ABI_FREE(chi0_lwing)
 ABI_FREE(chi0_uwing)
 ABI_FREE(chi0_head)

 call ltg_q%free()
 call gsph%free()

 call cwtime_report(" gwr_build_chi0_head_and_wings:", cpu_all, wall_all, gflops_all)
 call timab(1927, 2, tsec)

end subroutine gwr_build_chi0_head_and_wings
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_sigxme
!! NAME
!!  gwr_build_sigxme
!!
!! FUNCTION
!! Compute matrix elements of the exchange part.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_build_sigxme(gwr, compute_qp)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 logical,optional,intent(in) :: compute_qp

!Local variables-------------------------------
!scalars
 integer :: nsppol, nspinor, ierr, my_ikf, band_sum, ii, jj, kb, il_b, iab !ig_start, ig,
 integer :: my_is, ikcalc, ikcalc_ibz, bmin, bmax, band, istwf_k, npw_k
 integer :: spin, isym, jb, is_idx, use_umklp
 integer :: spad, wtqm, wtqp, irow, spadx1, spadx2
 integer :: npwx, u_nfft, u_mgfft, u_mpw, nsig_ab
 integer :: ik_bz, ik_ibz, isym_k, trev_k, g0_k(3)
 integer :: iq_bz, iq_ibz, isym_q, trev_q, g0_q(3)
 logical :: isirr_k, isirr_q, only_diago, sigc_is_herm, compute_qp__
 real(dp) :: fact_spin, theta_mu_minus_esum, theta_mu_minus_esum2, tol_empty, tol_empty_in, gwr_boxcutmin_x
 real(dp) :: cpu_k, wall_k, gflops_k, cpu_all, wall_all, gflops_all
 character(len=5000) :: msg
 logical :: q_is_gamma
 type(__slkmat_t),pointer :: ugb_kibz
 type(crystal_t),pointer :: cryst
 type(dataset_type),pointer :: dtset
 type(littlegroup_t) :: ltg_k
 type(desc_t),pointer :: desc_ki
!arrays
 integer :: g0(3), gmax(3), spinor_padx(2,4), u_ngfft(18), work_ngfft(18)
 integer,allocatable :: gbound_kcalc(:,:), gvec_x(:,:), gbound_x(:,:), kg_k(:,:), gbound_ksum(:,:)
 real(dp) :: ksum(3), kk_ibz(3), kgw(3), kgw_m_ksum(3), qq_bz(3), tsec(2) !, kk_bz(3), q0(3) !, spinrot_kbz(4), spinrot_kgw(4)
 real(dp),contiguous, pointer :: ks_eig(:,:,:), qp_eig(:,:,:), qp_occ(:,:,:), cg2_ptr(:,:) ! cg1_ptr(:,:),
 real(dp),allocatable :: work(:,:,:,:), cg1_ibz(:,:) !, cg2_bz(:,:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:)
 complex(dp),allocatable :: rhotwg(:), rhotwgp(:), rhotwg_ki(:,:)
 complex(gwpc),allocatable :: ur_bdgw(:,:)
 complex(dp),allocatable :: ur_ksum(:), ur_prod(:), eig0r(:)
 complex(dp),target,allocatable :: ug_ksum(:)
 complex(dp),allocatable  :: sigxcme_tmp(:,:), sigxme_tmp(:,:,:), sigx(:,:,:,:)
 complex(dp) :: gwpc_sigxme, gwpc_sigxme2, xdot_tmp
 type(sigijtab_t),allocatable :: Sigxij_tab(:,:), Sigcij_tab(:,:)

! *************************************************************************

 call timab(1920, 1, tsec)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 nsppol = gwr%nsppol; nspinor = gwr%nspinor; cryst => gwr%cryst; dtset => gwr%dtset
 nsig_ab = nspinor ** 2

 ! Allocate array with Sigma_x matrix elements.
 ii = minval(gwr%bstart_ks); jj = maxval(gwr%bstop_ks)
 ABI_RECALLOC(gwr%x_mat, (ii:jj, ii:jj, gwr%nkcalc, gwr%nsppol * nsig_ab))

 ! Table for \Sigmax_ij matrix elements.
 only_diago = .True.
 sigc_is_herm = .False.
 call sigtk_sigma_tables(gwr%nkcalc, gwr%nkibz, gwr%nsppol, gwr%bstart_ks, gwr%bstop_ks, gwr%kcalc2ibz(:,1), &
                         only_diago, sigc_is_herm, sigxij_tab, sigcij_tab)

 call sigijtab_free(Sigcij_tab)
 ABI_FREE(Sigcij_tab)

 if (only_diago) then
   call wrtout([std_out, ab_out], " Computing diagonal matrix elements of Sigma_x", pre_newlines=1)
 else
   call wrtout([std_out, ab_out], " Computing diagonal + off-diagonal matrix elements of Sigma_x", pre_newlines=1)
 end if

 ks_eig => gwr%ks_ebands%eig
 if (gwr%scf_iteration == 1) then
   call wrtout([std_out, ab_out], " Using KS orbitals and KS energies...", newlines=1, do_flush=.True.)
   qp_eig => gwr%ks_ebands%eig; qp_occ => gwr%ks_ebands%occ
 else
   call wrtout([std_out, ab_out], " Using KS orbitals and QP energies...", newlines=1, do_flush=.True.)
   qp_eig => gwr%qp_ebands%eig; qp_occ => gwr%qp_ebands%occ
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

 ! =========================================
 ! Find FFT mesh and max number of g-vectors
 ! =========================================
 gwr_boxcutmin_x = two
 call gwr%get_u_ngfft(gwr_boxcutmin_x, u_ngfft, u_nfft, u_mgfft, u_mpw, gmax)

 if (gwr%comm%me == 0) then
   call print_ngfft(u_ngfft, header="FFT mesh for Sigma_x", unit=std_out)
   call print_ngfft(u_ngfft, header="FFT mesh for Sigma_x", unit=ab_out)
 end if

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
 gmax = 2 * gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 do my_is=1,gwr%my_nspins
 spin = gwr%my_spins(my_is)
 do ikcalc=1,gwr%nkcalc ! TODO: Should be spin dependent!
   call cwtime(cpu_k, wall_k, gflops_k, "start")
   ikcalc_ibz = gwr%kcalc2ibz(ikcalc, 1)
   kgw = gwr%kcalc(:, ikcalc)
   bmin = gwr%bstart_ks(ikcalc, spin); bmax = gwr%bstop_ks(ikcalc, spin)

   ! ==============================================================
   ! ==== Find little group of the k-points for GW corrections ====
   ! ==============================================================
   ! * The little group is used only if symsigma == 1
   ! * If use_umklp == 1 then symmetries requiring an umklapp to preserve k_gw are included as well.
   use_umklp = 1
   call ltg_k%init(kgw, gwr%nqbz, gwr%qbz, cryst, use_umklp, npwe=0)

   write(msg,'(5a)') ch10, &
    ' Calculating <nk|Sigma_x|nk> at k: ',trim(ktoa(kgw)), ", for band range: ", trim(ltoa([bmin, bmax]))
   call wrtout(std_out, msg)

   ! ===============================================
   ! Load wavefunctions for Sigma_x matrix elements
   ! ===============================================
   ! All procs need ur_bdgw but the IBZ is distributed and, possibly, replicated in gwr%kpt_comm.
   ! Here we select the right procs, fill the buffer with the FFT results and then use
   ! a dumb xmpi_sum + rescaling to gather the results.
   ! FIXME: g-vectors from Green's descriptor or use another array to be able to deal with istwfk == 2?

   ABI_MALLOC_OR_DIE(ur_bdgw, (u_nfft * nspinor, bmin:bmax), ierr)
   ur_bdgw = zero

   if (any(ikcalc_ibz == gwr%my_kibz_inds)) then
     associate (desc_kcalc => gwr%green_desc_kibz(ikcalc_ibz), ugb_kcalc => gwr%ugb(ikcalc_ibz, spin))
     ABI_MALLOC(gbound_kcalc, (2 * u_mgfft + 8, 2))
     call sphereboundary(gbound_kcalc, desc_kcalc%istwfk, desc_kcalc%gvec, u_mgfft, desc_kcalc%npw)

     do il_b=1,ugb_kcalc%sizeb_local(2)
       band = ugb_kcalc%loc2gcol(il_b); if (band < bmin .or. band > bmax) CYCLE
       call fft_ug(desc_kcalc%npw, u_nfft, nspinor, ndat1, &
                   u_mgfft, u_ngfft, desc_kcalc%istwfk, desc_kcalc%gvec, gbound_kcalc, &
                   gwr%ugb(ikcalc_ibz, spin)%buffer_cplx(:, il_b), &  ! in
                   ur_bdgw(:, band))                                  ! out
     end do
     ABI_FREE(gbound_kcalc)
     end associate
   end if

   ! Collect and rescale
   call xmpi_sum(ur_bdgw, gwr%kgt_comm%value, ierr)
   ur_bdgw = ur_bdgw / gwr%np_kibz(ikcalc_ibz)

   ABI_MALLOC(ur_prod, (u_nfft * nspinor))
   ABI_MALLOC(ur_ksum, (u_nfft * nspinor))
   ABI_MALLOC(eig0r, (u_nfft * nspinor))

   ABI_CALLOC(sigxme_tmp, (bmin:bmax, bmin:bmax, nsppol * nsig_ab))
   ABI_CALLOC(sigxcme_tmp, (bmin:bmax, nsppol * nsig_ab))
   ABI_CALLOC(sigx, (2, bmin:bmax, bmin:bmax, nsppol * nsig_ab))

   ! ========================================
   ! ==== Sum over my k-points in the BZ ====
   ! ========================================

   do my_ikf=1,gwr%my_nkbz
     ik_bz = gwr%my_kbz_inds(my_ikf)
     ksum = gwr%kbz(:, ik_bz)

     ! Find the symmetrical image of ksum in the IBZ
     !call kmesh%get_BZ_item(ik_bz, ksum, ik_ibz, isym_ki, iik, ph_mkt)

     ! FIXME: Be careful with the symmetry conventions here and the interplay between umklapp in q and FFT
     ik_ibz = gwr%kbz2ibz_symrel(1, ik_bz); isym_k = gwr%kbz2ibz_symrel(2, ik_bz)
     trev_k = gwr%kbz2ibz_symrel(6, ik_bz); g0_k = gwr%kbz2ibz_symrel(3:5, ik_bz)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     kk_ibz = gwr%kibz(:, ik_ibz)

     ! Identify q and G0 where q + G0 = k_GW - ksum
     kgw_m_ksum = kgw - ksum
     call findqg0(iq_bz, g0, kgw_m_ksum, gwr%nqbz, gwr%qbz, gwr%mG0)
     !ABI_CHECK(all(g0 == 0), sjoin("g0 = ", ltoa(g0)))

     call calc_ceigr(g0, u_nfft, nspinor, u_ngfft, eig0r)

     ! If symmetries are exploited, only q-points in the IBZ_k are computed.
     ! In this case elements are weighted according to wtqp and wtqm. wtqm is for time-reversal.
     wtqp = 1; wtqm = 0
     !if (can_symmetrize(spin)) then
     if (gwr%dtset%symsigma == 1) then
       if (ltg_k%ibzq(iq_bz) /= 1) CYCLE
       wtqp = 0; wtqm = 0
       do isym=1,ltg_k%nsym_sg
         wtqp = wtqp + ltg_k%wtksym(1, isym, iq_bz)
         wtqm = wtqm + ltg_k%wtksym(2, isym, iq_bz)
       end do
     end if

     qq_bz = gwr%qbz(:, iq_bz)
     iq_ibz = gwr%qbz2ibz(1, iq_bz); isym_q = gwr%qbz2ibz(2, iq_bz)
     trev_q = gwr%qbz2ibz(6, iq_bz); g0_q = gwr%qbz2ibz(3:5, iq_bz)
     isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))

     ! Find the corresponding irreducible q-point.
     ! NB: non-zero umklapp G_o is not allowed. There's a check in setup_sigma
     !call qmesh%get_BZ_item(iq_bz, qbz, iq_ibz, isym_q, itim_q)
     q_is_gamma = normv(qq_bz, cryst%gmet, "G") < GW_TOLQ0
     call get_kg(qq_bz, istwfk1, dtset%ecutsigx, cryst%gmet, npwx, gvec_x)

     ABI_MALLOC(gbound_x, (2*u_mgfft + 8, 2))
     call sphereboundary(gbound_x, istwfk1, gvec_x, u_mgfft, npwx)

     ! Tables for the FFT of the oscillators.
     !  a) FFT index of G-G0.
     !  b) x_gbound table for the zero-padded FFT performed in rhotwg.
     !ABI_MALLOC(x_gbound, (2*u_mgfft+8, 2))
     !call Gsph_x%fft_tabs(g0, u_mgfft, u_ngfft, use_padfft, x_gbound, igfftxg0)

     ABI_MALLOC(rhotwg_ki, (npwx * nspinor, bmin:bmax))
     ABI_MALLOC(rhotwg, (npwx * nspinor))
     ABI_MALLOC(rhotwgp, (npwx * nspinor))
     ABI_MALLOC(vc_sqrt_qbz, (npwx))
     spinor_padx = reshape([0, 0, npwx, npwx, 0, npwx, npwx, 0], [2, 4])

     ! Get Fourier components of the Coulomb interaction in the BZ
     ! In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     call gwr%vcgen%get_vc_sqrt(qq_bz, npwx, gvec_x, gwr%q0, gwr%cryst, vc_sqrt_qbz, gwr%gtau_comm%value)

     desc_ki => gwr%green_desc_kibz(ik_ibz)

     ! Get npw_k, kg_k for this k.
     if (isirr_k) then
       istwf_k = desc_ki%istwfk; npw_k = desc_ki%npw
       ABI_MALLOC(kg_k, (3, npw_k))
       kg_k(:,:) = desc_ki%gvec
     else
       istwf_k = 1
       call get_kg(ksum, istwf_k, dtset%ecut, cryst%gmet, npw_k, kg_k)
     end if

     ABI_MALLOC(ug_ksum, (npw_k * nspinor))
     ABI_MALLOC(cg1_ibz, (2, desc_ki%npw * nspinor))
     !ABI_MALLOC(cg2_bz, (2, npw_k * nspinor))

     ABI_MALLOC(gbound_ksum, (2*u_mgfft+8, 2))
     call sphereboundary(gbound_ksum, istwf_k, kg_k, u_mgfft, npw_k)

     ! ==========================
     ! Sum over (occupied) bands
     ! ==========================
     ugb_kibz => gwr%ugb(ik_ibz, spin)

     do il_b=1,ugb_kibz%sizeb_local(2)
       band_sum = ugb_kibz%loc2gcol(il_b)

       ! Skip empty states. MRM: allow negative occ numbers.
       if (abs(qp_occ(band_sum, ik_ibz, spin)) < tol_empty) CYCLE

       !call wfd%get_ur(band_sum, ik_ibz, spin, ur_ibz)

       ! Compute ur_ksum(r) from the symmetrical image.
       ! I should rotate the g-vectors outside the loop and rotate ug here
       ! but at present I cannot use cgtk_rotate due to the symrel^T convention.

       if (isirr_k) then
         !call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, bra_kq)
         ug_ksum(:) = ugb_kibz%buffer_cplx(:, il_b)
       else
         ! Reconstruct u_kq(G) from the IBZ image.
         !call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, cgwork)

         ! FIXME: This is wrong if spc
         call c_f_pointer(c_loc(ug_ksum), cg2_ptr, shape=[2, npw_k * nspinor])

         !call c_f_pointer(c_loc(ugb_kibz%buffer_cplx(:, il_b)), cg1_ptr, shape=[2, desc_ki%npw * nspinor])
         !call cgtk_rotate(cryst, kk_ibz, isym_k, trev_k, g0_k, nspinor, ndat1, &
         !                 desc_ki%npw, desc_ki%gvec, &
         !                 npw_k, kg_k, desc_ki%istwfk, istwf_k, cg1_ptr, cg2_ptr, work_ngfft, work)

         cg1_ibz(1,:) = real(ugb_kibz%buffer_cplx(:, il_b))
         cg1_ibz(2,:) = aimag(ugb_kibz%buffer_cplx(:, il_b))
         call cgtk_rotate(cryst, kk_ibz, isym_k, trev_k, g0_k, nspinor, ndat1, &
                          desc_ki%npw, desc_ki%gvec, &
                          npw_k, kg_k, desc_ki%istwfk, istwf_k, cg1_ibz, cg2_ptr, work_ngfft, work)
       end if

       call fft_ug(npw_k, u_nfft, nspinor, ndat1, u_mgfft, u_ngfft, istwf_k, kg_k, gbound_ksum, &
                   ug_ksum, ur_ksum)

       if (any(g0 /= 0)) ur_ksum = ur_ksum * conjg(eig0r)

       ! Get all <k-q,band_sum,s|e^{-i(q+G).r}|s,jb,k>
       do jb=bmin,bmax

         ! FIXME: nspinor 2 is wrong as we have a 2x2 matrix
         ur_prod(:) = conjg(ur_ksum(:)) * ur_bdgw(:,jb)
         call fft_ur(npwx, u_nfft, nspinor, ndat1, u_mgfft, u_ngfft, istwfk1, gvec_x, gbound_x, &
                     ur_prod, rhotwg_ki(:,jb))

         ! Multiply by the square root of the Coulomb term
         ! In 3-D systems, the factor sqrt(4pi) is included
         do ii=1,nspinor
           spad = (ii-1) * npwx
           rhotwg_ki(spad+1:spad+npwx,jb) = rhotwg_ki(spad+1:spad + npwx,jb) * vc_sqrt_qbz(1:npwx)
         end do

         if (q_is_gamma) then
         !if (ik_bz == jk_bz) then
           ! Treat analytically the case q --> 0:
           !
           !   * The oscillator is evaluated at q = 0 as it is considered constant in the small cube around Gamma
           !     while the Colulomb term is integrated out.
           !   * If nspinor == 1, we have nonzero contribution only if band_sum == jb
           !   * If nspinor == 2, we evaluate <band_sum,up|jb,up> and <band_sum,dwn|jb,dwn>,
           !     and impose orthonormalization since npwwfn might be < npwvec.
           !   * Note the use of i_sz_resid and not i_sz, to account for the possibility
           !     to have generalized KS basis set from hybrid

           if (nspinor == 1) then
             rhotwg_ki(1, jb) = czero_gw
             if (band_sum == jb) rhotwg_ki(1,jb) = cmplx(sqrt(gwr%vcgen%i_sz), 0.0_gwp)
             !rhotwg_ki(1,jb) = czero_gw ! DEBUG

           else
             !ABI_ERROR("Not implemented Error")
             rhotwg_ki(1, jb) = zero; rhotwg_ki(npwx+1, jb) = zero
             if (band_sum == jb) then
               !ABI_CHECK(wfd%get_wave_ptr(band_sum, ik_ibz, spin, wave_sum, msg) == 0, msg)
               !cg_sum => wave_sum%ug
               !ABI_CHECK(wfd%get_wave_ptr(jb, jk_ibz, spin, wave_jb, msg) == 0, msg)
               !cg_jb  => wave_jb%ug
               !ctmp = xdotc(npw_k, cg_sum(1:), 1, cg_jb(1:), 1)
               rhotwg_ki(1, jb) = cmplx(sqrt(gwr%vcgen%i_sz), 0.0_gwp)  !* real(ctmp)
               !ctmp = xdotc(npw_k, cg_sum(npw_k+1:), 1, cg_jb(npw_k+1:), 1)
               rhotwg_ki(npwx+1, jb) = cmplx(sqrt(gwr%vcgen%i_sz), 0.0_gwp) ! * real(ctmp)
             end if
             !!!rhotwg_ki(1, jb) = zero; rhotwg_ki(npwx+1, jb) = zero
             !!! PAW is missing
           end if
         end if

       end do ! jb Got all matrix elements from bmin up to bmax.

       theta_mu_minus_esum  = fact_spin * qp_occ(band_sum, ik_ibz, spin)
       theta_mu_minus_esum2 = sqrt(abs(fact_spin * qp_occ(band_sum, ik_ibz, spin))) ! MBB Nat. orb. funct. approx. sqrt(occ)

       if (abs(theta_mu_minus_esum / fact_spin) >= tol_empty) then     ! MRM: allow negative occ numbers
         do kb=bmin,bmax

           ! Copy the ket Sigma_x |phi_{k,kb}>.
           rhotwgp(:) = rhotwg_ki(:, kb)

           ! Loop over the non-zero row elements of this column.
           ! If gwcalctyp <  20: only diagonal elements since QP == KS.
           ! If gwcalctyp >= 20:
           !      * Only off-diagonal elements connecting states with same character.
           !      * Only the upper triangle if HF, SEX, or COHSEX.

           do irow=1,Sigxij_tab(ikcalc, spin)%col(kb)%size1
             jb = Sigxij_tab(ikcalc, spin)%col(kb)%bidx(irow)
             rhotwg(:) = rhotwg_ki(:,jb)

             ! Calculate bare exchange <phi_jb|Sigma_x|phi_kb>.
             ! Do the scalar product only if band_sum is occupied.
             do iab=1,nsig_ab
               spadx1 = spinor_padx(1, iab); spadx2 = spinor_padx(2, iab)
               xdot_tmp = -XDOTC(npwx, rhotwg(spadx1+1:), 1, rhotwgp(spadx2+1:), 1)
               gwpc_sigxme  = xdot_tmp * theta_mu_minus_esum
               gwpc_sigxme2 = xdot_tmp * theta_mu_minus_esum2

               ! Accumulate and symmetrize Sigma_x matrix elements.
               ! -wtqm comes from time-reversal (exchange of band indeces)
               is_idx = spin; if (nspinor == 2) is_idx = iab
               sigxme_tmp(jb, kb, is_idx) = sigxme_tmp(jb, kb, is_idx) + &
                  (wtqp + wtqm) * DBLE(gwpc_sigxme) + (wtqp - wtqm) * j_dpc * AIMAG(gwpc_sigxme)
               if (jb == kb) then
                 sigxcme_tmp(jb, is_idx) = sigxcme_tmp(jb, is_idx) + &
                   (wtqp + wtqm) * DBLE(gwpc_sigxme2) + (wtqp - wtqm) *j_dpc * AIMAG(gwpc_sigxme2)
               end if

               sigx(1, jb, kb, is_idx) = sigx(1, jb, kb, is_idx) + wtqp *      gwpc_sigxme
               sigx(2, jb, kb, is_idx) = sigx(2, jb, kb, is_idx) + wtqm *CONJG(gwpc_sigxme)
             end do
           end do ! irow

         end do ! kb
       end if

     end do ! band_sum

     ABI_FREE(gbound_x)
     ABI_FREE(kg_k)
     ABI_FREE(ug_ksum)
     ABI_FREE(cg1_ibz)
     !ABI_FREE(cg2_bz)
     ABI_FREE(gbound_ksum)
     ABI_FREE(gvec_x)
     ABI_FREE(rhotwg_ki)
     ABI_FREE(rhotwg)
     ABI_FREE(rhotwgp)
     ABI_FREE(vc_sqrt_qbz)
   end do ! my_ikf Got all diagonal (off-diagonal) matrix elements.

   ! Gather contributions from all the CPUs.
   call xmpi_sum(sigxme_tmp, gwr%kgt_comm%value, ierr)
   call xmpi_sum(sigxcme_tmp, gwr%kgt_comm%value, ierr)
   call xmpi_sum(sigx, gwr%kgt_comm%value, ierr)

   ! Multiply by constants. For 3D systems sqrt(4pi) is included in vc_sqrt_qbz.
   sigxme_tmp  = (one / (cryst%ucvol * gwr%nkbz)) * sigxme_tmp  ! * Sigp%sigma_mixing
   sigxcme_tmp = (one / (cryst%ucvol * gwr%nkbz)) * sigxcme_tmp ! * Sigp%sigma_mixing
   sigx        = (one / (cryst%ucvol * gwr%nkbz)) * sigx        ! * Sigp%sigma_mixing

   ! If we have summed over the IBZ_q, we have to average over degenerate states.
   ! Presently only diagonal terms are considered
   ! Note that here we pass ks_eig to sigx_symmetrize instead of qp_eig.
   ! The reason is that we use the eigenvalues to detect degeneracies before averaging
   ! and qp_eig may break degeneracies while ks_eig are much more accurate.
   ! Most of the breaking comes from the correlated part, likey due to the treatment of q --> 0.

   ! TODO QP-SCGW required a more involved approach, there is a check in sigma
   ! TODO it does not work if nspinor == 2.

   if (gwr%dtset%symsigma == 1) then
     call sigx_symmetrize(ikcalc_ibz, spin, bmin, bmax, nsppol, nspinor, nsig_ab, ks_eig, sigx, sigxme_tmp)
     !do ii=bmin, bmax; print *, "qp_eig:", ii, qp_eig(ii, ikcalc_ibz, spin) * Ha_eV; end do
     !call sigx_symmetrize(ikcalc_ibz, spin, bmin, bmax, nsppol, nspinor, nsig_ab, qp_eig, sigx, sigxme_tmp)
   end if

   ! Reconstruct the full sigma_x matrix from the upper triangle.
   !if (nspinor == 1) then
   !  call hermitianize(sigxme_tmp, "Upper")
   !else
   !  ABI_WARNING("Should hermitianize non-collinear sigma!")
   !end if

   ! Save exchange matrix in gwr%x_mat
   if (gwr%nspinor == 1) then
     gwr%x_mat(bmin:bmax, bmin:bmax, ikcalc, spin) = sigxme_tmp(bmin:bmax, bmin:bmax, spin)
   else
     gwr%x_mat(bmin:bmax, bmin:bmax, ikcalc, :) = sigxme_tmp(bmin:bmax, bmin:bmax, :)
   end if

   ABI_FREE(ur_bdgw)
   ABI_FREE(ur_prod)
   ABI_FREE(ur_ksum)
   ABI_FREE(eig0r)
   ABI_FREE(sigxme_tmp)
   ABI_FREE(sigxcme_tmp)
   ABI_FREE(sigx)
   call ltg_k%free()
   call cwtime_report(" Sigx_nk:", cpu_k, wall_k, gflops_k)
 end do ! ikcalc
 end do ! my_is

 if (gwr%spin_comm%nproc > 1) call xmpi_sum(gwr%x_mat, gwr%spin_comm%value, ierr)

 ABI_FREE(work)
 call sigijtab_free(Sigxij_tab)
 ABI_FREE(Sigxij_tab)

 ! Compute QP results. Done usually when gwr_task == G0v i.e. Hartree-Fock with KS states.
 compute_qp__ = .False.; if (present(compute_qp)) compute_qp__ = compute_qp
 if (compute_qp__ .and. gwr%comm%me == 0) then
   call write_notations([std_out, ab_out])
 end if

 call cwtime_report(" gwr_build_sigxme:", cpu_all, wall_all, gflops_all)
 call timab(1920, 2, tsec)

end subroutine gwr_build_sigxme
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_get_u_ngfft
!! NAME
!!  gwr_get_u_ngfft
!!
!! FUNCTION
!!  Compute FFT mesh from boxcutmin.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwr_get_u_ngfft(gwr, boxcutmin, u_ngfft, u_nfft, u_mgfft, u_mpw, gmax)

!Arguments ------------------------------------
 class(gwr_t),intent(in) :: gwr
 real(dp),intent(in) :: boxcutmin
 integer,intent(out) :: u_ngfft(18), u_nfft, u_mgfft, u_mpw, gmax(3)

!Local variables-------------------------------
 integer :: ik_bz, npw_, ig, ii
 real(dp) :: kk_bz(3)
 integer,allocatable :: gvec_(:,:)

! *************************************************************************

 ! All MPI procs execute this part.
 ! Note the loops over the full BZ to compute u_mpw
 ! FIXME: umklapp, ecutsigx and q-centered G-sphere
 ! TODO: Write new routine to compute best FFT mesh for ecut1 + ecut1. See set_mesh from GW code.

 u_ngfft = gwr%dtset%ngfft ! This to allow users to specify fftalg

 u_mpw = -1; gmax = 0
 do ik_bz=1,gwr%nkbz
   kk_bz = gwr%kbz(:, ik_bz)
   call get_kg(kk_bz, istwfk1, gwr%dtset%ecut, gwr%cryst%gmet, npw_, gvec_)
   u_mpw = max(u_mpw, npw_)
   ! TODO: g0 umklapp here can enter into play gmax may not be large enough!
   do ig=1,npw_
     do ii=1,3
      gmax(ii) = max(gmax(ii), abs(gvec_(ii, ig)))
     end do
   end do
   ABI_FREE(gvec_)
   call getng(boxcutmin, gwr%dtset%chksymtnons, gwr%dtset%ecut, gwr%cryst%gmet, &
              kk_bz, me_fft0, u_mgfft, u_nfft, u_ngfft, nproc_fft1, gwr%cryst%nsym, paral_fft0, &
              gwr%cryst%symrel, gwr%cryst%tnons, use_gpu_cuda=gwr%dtset%use_gpu_cuda, unit=dev_null)
 end do

end subroutine gwr_get_u_ngfft
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/get_1d_sc_phases_dpc
!! NAME
!!  get_1d_sc_phases_dpc
!!
!! FUNCTION
!!  Compute one-dimensional factors in the supercell.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine get_1d_sc_phases_dpc(sc_shape, nkpt, kpts, ph1d)

!Arguments ------------------------------------
 integer,intent(in) :: sc_shape(3), nkpt
 real(dp),intent(in) :: kpts(3, nkpt)
 complex(dp),allocatable,intent(out) :: ph1d(:,:,:)

!Local variables-------------------------------
 integer :: ikpt, ix, iy, iz
 real(dp) :: arg, fact, kk(3)

! *************************************************************************

 ABI_MALLOC(ph1d, (maxval(sc_shape), 3, nkpt))

 do ikpt=1,nkpt
   kk = kpts(:, ikpt)
   fact = two_pi * kk(1)
   do ix=0,sc_shape(1) - 1
     arg = fact * ix
     ph1d(ix + 1, 1, ikpt) = cmplx(cos(arg), sin(arg), kind=dp)
   end do
   fact = two_pi * kk(2)
   do iy=0,sc_shape(2) - 1
     arg = fact * iy
     ph1d(iy + 1, 2, ikpt) = cmplx(cos(arg), sin(arg), kind=dp)
   end do
   fact = two_pi * kk(3)
   do iz=0,sc_shape(3) - 1
     arg = fact * iz
     ph1d(iz + 1, 3, ikpt) = cmplx(cos(arg), sin(arg), kind=dp)
   end do
 end do ! ikpt

end subroutine get_1d_sc_phases_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/get_1d_sc_phases_spc
!! NAME
!!  get_1d_sc_phases_spc
!!
!! FUNCTION
!!  Compute one-dimensional factors in the supercell.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine get_1d_sc_phases_spc(sc_shape, nkpt, kpts, ph1d)

!Arguments ------------------------------------
 integer,intent(in) :: sc_shape(3), nkpt
 real(dp),intent(in) :: kpts(3, nkpt)
 complex(spc),allocatable,intent(out) :: ph1d(:,:,:)

!Local variables-------------------------------
 integer :: ikpt, ix, iy, iz
 real(dp) :: arg, fact, kk(3)

! *************************************************************************

 ABI_MALLOC(ph1d, (maxval(sc_shape), 3, nkpt))

 do ikpt=1,nkpt
   kk = kpts(:, ikpt)
   fact = two_pi * kk(1)
   do ix=0,sc_shape(1) - 1
     arg = fact * ix
     ph1d(ix + 1, 1, ikpt) = cmplx(cos(arg), sin(arg), kind=sp)
   end do
   fact = two_pi * kk(2)
   do iy=0,sc_shape(2) - 1
     arg = fact * iy
     ph1d(iy + 1, 2, ikpt) = cmplx(cos(arg), sin(arg), kind=sp)
   end do
   fact = two_pi * kk(3)
   do iz=0,sc_shape(3) - 1
     arg = fact * iz
     ph1d(iz + 1, 3, ikpt) = cmplx(cos(arg), sin(arg), kind=sp)
   end do
 end do ! ikpt

end subroutine get_1d_sc_phases_spc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/sc_sum_dpc
!! NAME
!!  sc_sum_dpc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine sc_sum_dpc(sc_shape, uc_ngfft, nspinor, ph1d, k_is_gamma, alpha, sc_data, uc_psi, cout)

!Arguments ------------------------------------
 integer,intent(in) :: sc_shape(3), uc_ngfft(18), nspinor
 complex(dp),intent(in) :: ph1d(maxval(sc_shape), 3)
 logical,intent(in) :: k_is_gamma
 complex(dp),target,intent(in) :: alpha, uc_psi(uc_ngfft(1)*uc_ngfft(2)*uc_ngfft(3)*nspinor)
 complex(dp),target,intent(in) :: &
    sc_data(uc_ngfft(1)*sc_shape(1)*uc_ngfft(2)*sc_shape(2)*uc_ngfft(3)*sc_shape(3)*nspinor)
 complex(dp),intent(out) :: cout

!Local variables-------------------------------
 integer :: il1, il2, il3, spinor, uc_n1, uc_n2, uc_n3, ix, iy, iz !, idat
 complex(dp) :: cphase, phl32, phl3
 complex(dp),contiguous,pointer :: uc_psi_ptr(:,:,:,:), sc_data_ptr(:,:,:,:,:,:,:)

! *************************************************************************

 uc_n1 = uc_ngfft(1); uc_n2 = uc_ngfft(2); uc_n3 = uc_ngfft(3)

 call c_f_pointer(c_loc(uc_psi), uc_psi_ptr, shape=[uc_n1, uc_n2, uc_n3, nspinor])
 call c_f_pointer(c_loc(sc_data), sc_data_ptr, &
                  shape=[uc_n1, sc_shape(1), uc_n2, sc_shape(2), uc_n3, sc_shape(3), nspinor])

 !do spinor=1,nspinor
 !end do
 ABI_CHECK(nspinor == 1, "nspinor 2 not coded")
 spinor = 1
 cout = zero

 if (k_is_gamma) then
   ! Don't need to multiply by e^{ik.L}
   do il3=1,sc_shape(3)
     do iz=1,uc_n3
       do il2=1,sc_shape(2)
         do iy=1,uc_n2
           do il1=1,sc_shape(1)
             do ix=1,uc_n1
               cout = cout + uc_psi_ptr(ix, iy, iz, spinor) * sc_data_ptr(ix, il1, iy, il2, iz, il3, spinor)
             end do
           end do
         end do
       end do
     end do
   end do

 else

   do il3=1,sc_shape(3)
     phl3 = ph1d(il3, 3)
     do iz=1,uc_n3
       do il2=1,sc_shape(2)
         phl32 = phl3 * ph1d(il2, 2)
         do iy=1,uc_n2
           do il1=1,sc_shape(1)
             cphase = phl32 * ph1d(il1, 1)  ! e^{ik.L}
             do ix=1,uc_n1
               cout = cout + cphase * uc_psi_ptr(ix, iy, iz, spinor) * sc_data_ptr(ix, il1, iy, il2, iz, il3, spinor)
             end do
           end do
         end do
       end do
     end do
   end do

 end if

 cout = alpha * cout

end subroutine sc_sum_dpc
!!***

!!****f* m_gwr/sc_sum_spc
!! NAME
!!  sc_sum_spc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine sc_sum_spc(sc_shape, uc_ngfft, nspinor, ph1d, k_is_gamma, alpha, sc_data, uc_psi, cout)

!Arguments ------------------------------------
 integer,intent(in) :: sc_shape(3), uc_ngfft(18), nspinor
 complex(sp),intent(in) :: ph1d(maxval(sc_shape), 3)
 logical,intent(in) :: k_is_gamma
 complex(sp),target,intent(in) :: alpha, uc_psi(uc_ngfft(1)*uc_ngfft(2)*uc_ngfft(3)*nspinor)
 complex(sp),target,intent(in) :: &
    sc_data(uc_ngfft(1)*sc_shape(1)*uc_ngfft(2)*sc_shape(2)*uc_ngfft(3)*sc_shape(3)*nspinor)
 complex(sp),intent(out) :: cout

!Local variables-------------------------------
 integer :: il1, il2, il3, spinor, uc_n1, uc_n2, uc_n3, ix, iy, iz !, idat
 complex(sp) :: cphase, phl32, phl3
 complex(sp),contiguous,pointer :: uc_psi_ptr(:,:,:,:), sc_data_ptr(:,:,:,:,:,:,:)

! *************************************************************************

 uc_n1 = uc_ngfft(1); uc_n2 = uc_ngfft(2); uc_n3 = uc_ngfft(3)

 call c_f_pointer(c_loc(uc_psi), uc_psi_ptr, shape=[uc_n1, uc_n2, uc_n3, nspinor])
 call c_f_pointer(c_loc(sc_data), sc_data_ptr, &
                  shape=[uc_n1, sc_shape(1), uc_n2, sc_shape(2), uc_n3, sc_shape(3), nspinor])


 !do spinor=1,nspinor
 !end do
 ABI_CHECK(nspinor == 1, "nspinor 2 not coded")
 spinor = 1
 cout = zero

 if (k_is_gamma) then
   ! Don't need to multiply by ! e^{ik.L}

   do il3=1,sc_shape(3)
     do iz=1,uc_n3
       do il2=1,sc_shape(2)
         do iy=1,uc_n2
           do il1=1,sc_shape(1)
             do ix=1,uc_n1
               cout = cout + uc_psi_ptr(ix, iy, iz, spinor) * sc_data_ptr(ix, il1, iy, il2, iz, il3, spinor)
             end do
           end do
         end do
       end do
     end do
   end do

 else

   do il3=1,sc_shape(3)
     phl3 = ph1d(il3, 3)
     do iz=1,uc_n3
       do il2=1,sc_shape(2)
         phl32 = phl3 * ph1d(il2, 2)
         do iy=1,uc_n2
           do il1=1,sc_shape(1)
             cphase = phl32 * ph1d(il1, 1)  ! e^{ik.L}
             do ix=1,uc_n1
               cout = cout + cphase * uc_psi_ptr(ix, iy, iz, spinor) * sc_data_ptr(ix, il1, iy, il2, iz, il3, spinor)
             end do
           end do
         end do
       end do
     end do
   end do

 end if

 cout = alpha * cout

end subroutine sc_sum_spc
!!***

integer pure function memlimited_step(start, stop, num_items, bsize, maxmem_mb) result(step)
  integer,intent(in) :: start, stop, num_items, bsize
  real(dp),intent(in) :: maxmem_mb

  real(dp) :: totmem_mb
  totmem_mb = one * (stop - start + 1) * num_items * bsize
  step = stop - start + 1
  if (totmem_mb > maxmem_mb) step = floor(totmem_mb / maxmem_mb)

end function memlimited_step
!!***

end module m_gwr
!!***

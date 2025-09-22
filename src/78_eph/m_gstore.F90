!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!  This module implements the gstore_t object that allows one to **precompute"" the e-ph matrix elements g
!!  and store them in memory with a MPI-distributed data structure.
!!  This approach is the most CPU-efficient one when one has to deal with
!!  algorithms in which the same g(q, k) is required several (many) times.
!!  Typical examples are iterative solvers for non-linear equations that are called inside a loop over T.
!!  At each iteration, indeed, we need g(q, k) and computing these quantities from scratch
!!  would be very expensive.
!!
!!  Note that g depends on the two wave vectors (q, k), two electron band indices (m, n),
!!  phonon mode nu with momentum q and spin index if nsppol == 2 (collinear case).
!!
!!  g(q, k) is therefore a sloppy notation for:
!!
!!          g(q, k) = <k+q,m,spin| \Delta_{q,\nu} V^{spin}_{scf} |k,n,spin>
!!
!!  There are lots of technical details that should be discussed but, roughly speaking,
!!  the gstore API allows one to:
!!
!!   - select whether q or k should be in the IBZ or in the BZ.
!!     NB: It is not possible to use the IBZ both for q and k as g(Sk, q) = g(k, S^{-1}q)
!!     thus one has to select the appropriate zones beforehand.
!!
!!   - filter bands and/or k/q wavevectors according to some criterion.
!!     In superconductors, for instance, only k/k+q states on the Fermi surface are usually needed.
!!     In semiconductors, one can include only k, k+q inside an energy window around the band edge.
!!     for transport properties or just the |n,k,spin> states at the band edges while <k+q,m,spin|
!!     have q in the BZ and m=1,nband.
!!
!!  - whether the code should compute and store the complex valued g or |g|^2.
!!    Expression depending of |g|^2 are gauge-invariant provided that all degenerate states are summed over.
!!    On the contrary, the complex valued g is gauge-dependent and hic sunt leones.
!!    In Abinit, the g elements are computed within the same gauge by reconstructing Bloch states
!!    in the BZ from the IBZ by using a deterministic symmetrization algorithm
!!    Client code reading the e-ph matrix elements produced by ABINIT is expected to follow the
!!    same conventions, especially if one needs to mix g with wavefunctions in the BZ.
!!
!!  At the level of the API, we have three different routines.
!!
!!      1) gstore_new builds the object, defines the BZ sampling type (e.g. k in the IBZ, q in the BZ)
!!         and implements filtering techniques. The MPI grid is automatically generated at this level.
!!
!!      2) gstore_compute evaluates the e-ph matrix elements in parallel and dumps the results to GSTORE.nc.
!!
!!      3) gstore%from_ncpath reconstructs the object from a GSTORE.nc file.
!!
!!  In a typical scenario, one uses eph_task 11 to generate GSTORE.nc i.e. steps 1) and 2).
!!  Then one introduces a new value of eph_task in which we read the object from file and call
!!  a specialized routine that implements the "post-processing" steps needed
!!  to compute the physical properties of interest.
!!
!!  Now, let us discuss the MPI-distribution.
!!
!!  The (q, k) matrix is distributed inside a 2D cartesian grid using block distribution.
!!  This is schematic representation for MPI 4 procs with 2 procs for k and 2 procs for q:
!!
!!                 k-axis (kpt_comm)
!!              |--------------------
!!              |         |         |
!!              |   P00   |   P01   |
!!              |         |         |
!!    q-axis    |--------------------
!!  (qpt_comm)  |         |         |
!!              |   P10   |   P11   |
!!              |         |         |
!!              |--------------------
!!
!!  Each processor stores all the (band_kq, band_k) transitions for a given (q, k) pair.
!!
!!  Perturbations can be optionally distributed along a third axis (pert_comm).
!!  Note, however, that the parallelism over perturbations is not expected to be the most efficient
!!  although it allows one to reduce the memory required to store the scattering potential in the supercell
!!  as we can distribute W(r, R, 3 * natom) over the last dimension.
!!
!!  For electronic properties, one usually uses k-points in the IBZ and q-points in the BZ.
!!  hence the parallelism over q-points is the most efficient one in terms of wall-time.
!!  Keep in mind, however, that the k-point parallelism allows one to reduce the memory allocated for the
!!  wavefunctions. Using some procs for k-points is also beneficial in terms of performance
!!  as we can reduce load imbalance is the number of procs in qpt_comm does not divide nqbz.
!!
!!  NB: If nsppol == 2, we create two gqk objects, one for each spin.
!!  The reason is that dimensions such as the number of effective bands/q-points/k-points
!!  depends on the collinear spin when filters are employed.
!!
!! TODO
!!  1) Implement possibility of reading a subset of data from a larger gstore ?
!!  2) Optimize v^1_loc|psi_nk> by precomputing <r|psi_nk> before the loop over my_npert
!!     Big speedup is expected, especially if one loops first over k and then q, provided
!!     the interpolation of v^1_q does not start to dominate
!!  3) Use similar trick in dfpt_cgw for H^0 |psi_nk>
!!  4) Operate on multiple n states in getgh1c (new version of getgh1c allows it)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gstore

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_clib
 use m_xmpi
 use m_errors
 use m_htetra
 use libtetrabz
 use netcdf
 use m_nctk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_ephtk
 use m_mkffnl
 use m_sigtk

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : pseudopotential_type
 use m_dtset,          only : dataset_type
 use m_dtfil,          only : datafiles_type
 use m_time,           only : cwtime, cwtime_report, sec2str
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat, replace_ch0, yesno, string_in
 use m_numeric_tools,  only : arth, get_diag, isdiagmat
 use m_krank,          only : krank_t, get_ibz2bz, star_from_ibz_idx
 use m_io_tools,       only : iomode_from_fname, file_exists
 use m_special_funcs,  only : gaussian
 use m_copy,           only : alloc_copy
 use m_fftcore,        only : ngfft_seq, get_kg
 use m_cgtools,        only : cg_zdotc
 use m_kg,             only : getph
 use m_crystal,        only : crystal_t
 use m_hdr,            only : hdr_type, fform_from_ext
 use m_matrix,         only : matr3inv
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_sort, kpts_pack_in_stars, &
                              kptrlatt_from_ngkpt
 use m_ebands,         only : ebands_t, gaps_t
 use m_lgroup,         only : lgroup_t
 use m_bz_mesh,        only : kmesh_t
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack
 use m_ifc,            only : ifc_type
 use m_phonons,        only : pheigvec_rotate
 use m_wfd,            only : wfd_t
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_mlwfovlp,       only : wan_t
 use m_pstat,          only : pstat_proc

 implicit none

 private

 character(len=abi_slen),public,parameter :: GSTORE_GMODE_ATOM   = "atom"
 character(len=abi_slen),public,parameter :: GSTORE_GMODE_PHONON = "phonon"

 ! Rank of the MPI Cartesian grid.
 integer,private,parameter :: ndims = 6
!!***

!----------------------------------------------------------------------

!!****t* m_gstore/gqk_t
!! NAME
!! gqk_t
!!
!! FUNCTION
!!  This object stores MPI-distributed e-ph matrix elements for
!!  a given spin index (if collinear magnetism i.e. nsppol 2).
!!  Local dimensions and arrays start with `my_`, global dimensions start with `glob_`
!!
!! SOURCE

type, public :: gqk_t

  integer :: cplex = -1
  ! 1 if |g|^2 is stored
  ! 2 if complex-valued g are stored (mind the gauge)

  integer :: spin = -1
  ! Spin index.

  integer :: natom3 = -1
  ! 3 * natom
  ! Mainly used to dimension arrays

  integer :: nb = -1
  ! Number of bands included in the calculation for this spin
  ! Global as this dimension is not MPI-distributed due to (m, n) pairs.
  ! NB: nb is not necessarily equal to nband.

  integer :: bstart = -1, bstop = -1
  ! The first band starts at bstart.
  ! The last band is bstop (NB: These are global indices)

  ! TODO
  ! These new entries will be used to implement band distribution.
  integer :: nb_kq = -1, nb_k = -1
  integer :: bstart_k = -1, bstop_k = -1
  integer :: bstart_kq = -1, bstop_kq = -1

  integer :: my_npert = -1
  ! Number of perturbations treated by this MPI rank.

  integer :: my_pert_start = -1
  ! Initial perturbation treated by this MPI proc

  integer :: glob_nk = -1, glob_nq = -1
  ! Total number of k/q points in global matrix.
  ! Note that k-points/q-points can be filtered. Use kzone, qzone and kfilter to interpret these dimensions.

  integer :: my_nk = -1, my_nq = -1
  ! Number of k/q points treated by this MPI proc. Used to loop and allocate local arrays.

  integer :: my_kstart = -1, my_qstart = -1
  ! Index of the first k/q point in the global matrix treated by this MPI proc

  integer,allocatable :: my_k2ibz(:,:)
  ! (6, my_nk)
  ! Mapping my_kpoints --> kibz (symrel conventions)

  !integer,allocatable :: my_k2bz(:,:)
  ! (my_nk)
  ! Mapping my_kpoints --> ik_bz

  real(dp),allocatable :: my_kpts(:,:)
  ! (3, my_nkpt)
  ! k-points treated by this MPI proc.

  real(dp),allocatable :: my_wtk(:)
  ! (my_nkpt)
  ! Weights for the k-points treated by this MPI proc.

  integer,allocatable :: my_q2ibz(:,:)
  ! (6, my_nq)
  ! Mapping my_qpoints --> qibz
  ! symrel conventions

  integer,allocatable :: my_q2bz(:)
  ! (my_nq)
  ! Mapping my_iq index --> iq_bz index in the full BZ

  integer,allocatable :: my_k2glob(:)
  ! (my_nk)
  ! Mapping my_ik index --> global index in the g(q, k) matrix.

  integer,allocatable :: my_q2glob(:)
  ! (my_nq)
  ! Mapping my_iq index --> global index in the g(q, k) matrix.

  real(dp),allocatable :: vk_cart_ibz(:,:,:)
  ! (3, nb, nkibz)
  ! Diagonal v_{n,k} for k in the IBZ.
  ! Values in the BZ can be reconstructed by symmetry.
  ! Allocated if gstore%with_vk == 1
  ! TODO: Here I should decide how to treat nb_k, nk_kq

  real(dp),allocatable :: vkmat_cart_ibz(:,:,:,:,:)
  ! (3, nb, nb, nkibz)
  ! v_{m, n,k} for the k in the IBZ
  ! Allocated if gstore%with_vk in (1, 2)
  ! TODO: Here I should decide how to treat nb_k, nk_kq

  integer,allocatable :: my_pertcases(:)
  ! (my_npert)
  ! List of perturbation indices treated by this MPI proc.
  ! Contiguous indices.

  complex(dp), allocatable :: my_g(:,:,:,:,:)
  ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
  ! (       p, b1_kq,     q, b2_k, k)  -->  <k+q, b1| D_{q,p}H |k, b2>
  ! e-ph matrix elements g (local buffer). Allocated if cplex == 2

  ! FIXME: I don't remember why I decided to have my_npert as first dimension
  ! now it seems much more more natural to me to have:

  ! (nb_kq, nb_k, my_npert, my_nk, my_nq) or
  ! (nb_kq, nb_k, my_npert, my_nq, my_nk)

  complex(dp), allocatable :: my_gq0nm_atm(:,:,:,:)
  ! (nb_k, nb_kq, natom3, my_nk)
  ! e-ph matrix elements g(k,q=0) required for the RIA DW term.
  ! Note: Perturbations are not distributed. Only k-points
  ! Also, m, n bands are exchanged (first n then m)

  real(dp), allocatable :: my_g2(:,:,:,:,:)
  ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
  ! |g|^2 (local buffer). Allocated if cplex == 1

  integer :: coords_qkpb_sumbp(ndims) = 0
  ! Coordinates of this processor in the (q, k, pert, band, band_sum, pp_sum) Cartesian grid.

  type(xcomm_t) :: kpt_comm
   ! MPI communicator over k-points

  type(xcomm_t) :: qpt_comm
   ! MPI communicator over q-points

  type(xcomm_t) :: qpt_kpt_comm
   ! MPI communicator over k/q subgrid

  type(xcomm_t) :: pert_comm
   ! MPI communicator over atomic perturbations.

  type(xcomm_t) :: band_comm
   ! MPI communicator for band distribution.

  type(xcomm_t) :: bsum_comm
   ! MPI communicator over bands in summation. NB: It is not used to distribute
   ! the memory for the g but to distribute a possible sum over bands as done in the GWPT code.

  type(xcomm_t) :: pp_sum_comm
   ! MPI communicator over wavevector summation. NB: It not used to distribute
   ! the memory for the g but to distribute a possible sum over wavevectors as done in the GWPT code.

  type(xcomm_t) :: qpt_pert_comm
   ! MPI communicator over the 2d grid (qpt, atomic perturbations)

  type(xcomm_t) :: pert_ppsum_comm
   ! MPI communicator over the 2d grid (atomic perturbations, pp_sum) used in GWPT

  type(xcomm_t) :: pert_ppsum_bsum_comm
   ! MPI communicator over the 3d grid (atomic perturbations, pp_sum, band_sum) used in GWPT

  type(xcomm_t) :: comm
   ! MPI communicator for full grid of procs treating this spin.

  type(wan_t) :: wan
   ! Object used to interpolate the e-ph matrix elements with Wannier.

  real(dp),allocatable :: my_wnuq(:,:)
  ! (my_npert, my_nq)
  ! Phonon frequencies in Ha (MPI distributed)

  real(dp),allocatable :: my_displ_cart(:,:,:,:,:)
  ! (2, 3, natom, my_npert, my_nq))
  ! Phonon displacements (MPI distributed)

 contains

  procedure :: gather => gqk_gather
  ! Gather the MPI-distributed matrix elements for a given k/q-point index

  procedure :: get_erange_mask => gqk_get_erange_mask
  ! Compute MPI-distributed & global mask for electronic states allowed by energy filtering

  procedure :: filter_erange => gqk_filter_erange
  ! Nullify all matrix elements connecting electronic states outside of specified erange

  procedure :: myqpt => gqk_myqpt
  ! Return the q-point and the weight from my local index my_iq

  procedure :: dbldelta_qpt => gqk_dbldelta_qpt
  ! Compute weights for the double delta.

  procedure :: free => gqk_free
  ! Free memory

 end type gqk_t
!!***

!----------------------------------------------------------------------

!!****t* m_gstore/gstore_t
!! NAME
!! gstore_t
!!
!! FUNCTION
!! This object stores:
!!
!!    - pointers to the crystal structure, the KS bands, the IFCs.
!!    - arrays that do not depend on the spin such as the IBZ and weights for k/q-points.
!!    - metadata such as kzone, qzone and kfilter that are needed to interpret
!!      the storage mode used for the g(k,q).
!!
!! NB: the e-ph matrix element are stored in gstore%gqk(my_is) where my_is counts
!!     the number of collinear spins treated by this MPI processor.
!!
!! SOURCE

type, public :: gstore_t

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: my_nspins = 0
   ! Number of collinear spins treated by this MPI rank

  integer :: nkibz = -1, nqibz = -1
  ! Number of k/q points in the IBZ.

  integer :: nkbz = -1, nqbz = -1
  ! Number of k/q points in the BZ.

  integer :: comm
   ! Global communicator
   ! Inherited by the caller thus we don't free it in gstore_free.

  integer :: with_vk = 0
  ! 0 if group velocities should not be computed
  ! 1 to compute diagonal terms only
  ! 2 to compute diagonal and off-diagonal terms

  integer :: qptopt = -1
  ! option for the generation of q points (defines whether spatial symmetries and/or time-reversal can be used)

  integer :: has_used_lgk = 0
  ! value of use_lgk used to generate GSTORE.nc (read from file).

  integer :: has_used_lgq = 0
  ! value of use_lgq used to generate GSTORE.nc (read from file).

  character(len=fnlen) :: path = " "
  ! Path to the nc file associated to the gstore

  character(len=fnlen) :: wfk0_path = " "

  character(len=abi_slen) :: kzone = " ", qzone = " "
   ! Specifies whether k- or q-points are in the BZ or in the IBZ.
   ! Possible values are "ibz" or "bz".
   ! Note that the combination ("ibz", "ibz") is not allowed.

  character(len=abi_slen) :: kfilter = "none"
  ! Specifies the technique used to filter k-points.
  ! Possible values: "none", "fs_tetra", "erange", "qprange"

  character(len=abi_slen) :: gmode = "none"
  ! "phonon" or "atom"

  character(len=abi_slen) :: gtype = "KS"
  ! Formalism used to compute g(k,q). Either KS or GWPT

  real(dp),allocatable :: erange_spin(:, :)
  ! (2, nsppol)
  ! Energy window. zero if not used. Requires kfilter == "erange"

  type(crystal_t), pointer :: cryst => null()
  ! Crystalline structure

  type(ebands_t), pointer :: ebands => null()
  ! Electron bands

  logical :: ebands_owns_memory = .False.
  ! True if ebands pointer owns memory and should therefore be deallocated in gstore_free
  ! TODO: To my future Matteo, here we have to be very careful if we try to update the energies with GW
  ! or dope/changhe the Fermi level.

  type(ifc_type), pointer :: ifc => null()
  ! interatomic force constants.

  type(dataset_type), pointer :: dtset => null()
  ! Reference to the dataset.

  type(krank_t) :: krank_ibz, qrank_ibz
  ! Object used to find k-points or q-points in the IBZ and map BZ to IBZ.

  integer :: ngqpt(3) = 0
  ! Number of grid points for q-points (either from ddb_ngqpt or eph_ngqpt_fine)

  integer,allocatable :: my_spins(:)
   ! (%my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2 and nspinor == 1

  integer,allocatable :: brange_spin(:, :)
  ! (2, nsppol)
  ! Range of bands for each spin
  ! TODO: Should be replaced by brange_k_spin and brange_kq_spin

  !integer :: max_nb = -1
  ! Max number of bands over spin

  integer,allocatable :: glob_nk_spin(:), glob_nq_spin(:)
  ! Total number of k/q points for each spin after filtering (if any)
  ! (nsppol)

  integer,allocatable :: kglob2bz(:,:)
  ! (max_nk, nsppol))
  ! Mapping ik_glob to BZ index for k-points.

  !integer,allocatable qglob2bz(:,:)
  ! (max_nq, nsppol))
  ! Mapping iq_glob to BZ index for q-points.

  integer,allocatable :: kbz2ibz(:,:)
  ! (6, gstore%nkbz))
  ! Mapping BZ --> IBZ (symrel conventions that can be used to symmetrize wavefunctions)

  real(dp), contiguous, pointer :: kibz(:,:)
  ! k-points in the IBZ. Points to ebands%kptns
  ! (3, nkibz)

  real(dp),allocatable :: delta_ef_kibz_spin(:,:,:)
  ! (nb, gstore%nkibz, nsppol))
  ! Tetrahedron weights at eF in the IBZ.

  real(dp), allocatable :: qibz(:,:)
  ! (3, nqibz)
  ! q-points in the IBZ in reduced coordinates.

  real(dp), allocatable :: wtq(:)
  ! (nqibz)
  ! q-points weights in the IBZ

  ! TODO: Use MPI shared memory
  real(dp),allocatable :: qbz(:,:)
  ! q-points in the BZ.

 ! TODO: Use MPI shared memory
  real(dp),allocatable :: kbz(:,:)
  ! k-points in the BZ.

  !integer :: qptrlatt(3, 3) = -1  ! kptrlatt(3, 3) = -1,
   ! k-mesh and q-mesh

  !real(dp),allocatable :: kshift(:, :), qshift(:, :)
  ! k/q-mesh shift (well, q-mesh is usually gamma-centered)

  type(gqk_t), allocatable :: gqk(:)
  ! (my_nspins)
  ! Datastructure storing e-ph matrix elements for the collinear spins treated by this MPI proc.

contains

  procedure :: fill_bks_mask => gstore_fill_bks_mask
  ! Fill the table used to read (b, k, s) wavefunctions from the WFK file
  ! keeping into account the distribution of the e-ph matrix elements.

  procedure :: fill_bks_mask_pp_mesh => gstore_fill_bks_mask_pp_mesh
  ! Fill the table used to read (b, k, s) wavefunctions from the WFK file
  ! keeping into account the distribution of the e-ph matrix elements in the GWPT code
  ! and the parallel distribution of the pp momenta.

  procedure :: get_mpw_gmax => gstore_get_mpw_gmax
  ! Compute the maximum number of PWs for all possible k+q treated.

  procedure :: spin2my_is => gstore_spin2my_is
  !  Return the local spin index from the global spin index.
  !  0 if this spin is not treated by this MPI proc.

  procedure :: free => gstore_free
  ! Free memory

  procedure :: print => gstore_print
  ! Print info on the object

  procedure :: check_little_group => gstore_check_little_group
   !  Check consistency between little group options from file and from input.

  procedure, private :: distribute_spins__ => gstore_distribute_spins
  ! Distribute spins, create indirect mapping to spin index and init %brange_spin

  procedure, private :: set_mpi_grid__ => gstore_set_mpi_grid__
  ! Set the MPI cartesian grid

  procedure, private :: malloc__ => gstore_malloc__
  ! Allocate local buffers once the MPI grid has been initialized.

  procedure, private :: filter_fs_tetra__ => gstore_filter_fs_tetra__
  ! Select k-points on the FS using the tetrahedron method

  procedure, private :: filter_erange__ => gstore_filter_erange__
  ! Select k-points inside an energy window.

  procedure, private :: filter_qprange__ => gstore_filter_qprange__
  ! Select k-points according to gw_qprange

  procedure :: compute => gstore_compute
  ! Compute e-ph matrix elements.

  procedure :: get_lambda_iso_iw => gstore_get_lambda_iso_iw
  ! Compute isotropic lambda(iw) along the imaginary axis.

  procedure :: get_a2fw => gstore_get_a2fw
  ! Compute Eliashberg function a^2F(w).

  procedure :: from_ncpath => gstore_from_ncpath
  ! Reconstruct object from netcdf file.

  procedure :: init => gstore_init
  ! Build object from scratch

  procedure :: get_missing_qbz_spin => gstore_get_missing_qbz_spin
  ! Return the number of (q-points, spin) entries that have been computed

  procedure :: set_perts_distrib => gstore_set_perts_distrib
  ! Activate parallelism over perturbations at the level of the DVDB file.

  procedure :: print_for_abitests => gstore_print_for_abitests
  ! Print subset of results to ab_out for testing purposes.

  procedure :: check_cplex_qkzone_gmode => gstore_check_cplex_qkzone_gmode
  ! Perform consistency checks.

  procedure :: wannierize_and_write_gwan => gstore_wannierize_and_write_gwan
  ! Compute g(R_e,R_ph) from g(k,q) and save results to GWAN.nc file

end type gstore_t
!!***

public :: gstore_check_restart
 ! Check whether restart is possible.

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_init
!! NAME
!! gstore_init
!!
!! FUNCTION
!! Initialize the object
!!
!! INPUTS
!! path=Filename of the output GSTORE.nc file
!!
!! SOURCE

subroutine gstore_init(gstore, path, dtset, dtfil, wfk0_hdr, cryst, ebands, ifc, comm, &
                       gtype) ! optional

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(out) :: gstore
 character(len=*),intent(in) :: path
 type(dataset_type),target,intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(in) :: wfk0_hdr
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),target,intent(in) :: ebands
 class(ifc_type),target,intent(in) :: ifc
 integer,intent(in) :: comm
 character(len=*),optional,intent(in) :: gtype

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: all_nproc, my_rank, ierr, my_nshiftq, nsppol, spin, natom3, cnt, timrev_q, with_cplex
 integer :: ik_ibz, ik_bz, iq_bz, iq_ibz, max_nq, max_nk, ncid, spin_ncid, ncerr, gstore_fform
 integer :: my_is, my_ik, my_iq, nq
 logical :: keep_umats, has_abiwan, has_gwan, write_gstore
 real(dp) :: cpu, wall, gflops, weight_qq, gstore_fill_dp
 character(len=5000) :: msg
!arrays
 integer :: ngqpt(3), qptrlatt(3,3), comm_spin(ebands%nsppol), nproc_spin(ebands%nsppol), units(2)
 integer,allocatable :: qbz2ibz(:,:), kibz2bz(:), qibz2bz(:), qglob2bz(:,:)
 integer,allocatable :: select_qbz_spin(:,:), select_kbz_spin(:,:)
 real(dp):: my_shiftq(3,1), kpt(3), kq(3), qpt(3)
 real(dp),allocatable :: wtk(:), kibz(:,:)
 type(wan_t),target :: wan_spin(ebands%nsppol)
 complex(dp),allocatable :: intp_gatm(:,:,:,:)
!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom3 = 3 * cryst%natom; nsppol = ebands%nsppol
 units = [std_out, ab_out]

 call wrtout(std_out, " gstore_init: building gstore_t instance...")
 call pstat_proc%print(_PSTAT_ARGS_)

 ! Set basic parameters.
 gstore%comm = comm; gstore%nsppol = nsppol; gstore%path = path
 if (present(gtype)) gstore%gtype = gtype

 ! Get references to other data structures.
 gstore%dtset => dtset
 gstore%cryst => cryst; gstore%ebands => ebands; gstore%ifc => ifc
 gstore%ebands_owns_memory = .False.

 has_abiwan = .False.; has_gwan = .False.; keep_umats = .False.
 if (dtfil%filabiwanin /= ABI_NOFILE) then
   has_abiwan = .True.
   call wrtout(units, sjoin(" Reading set of bands to be included in gstore computation from ABIWAN file:", dtfil%filabiwanin))
   do spin=1,ebands%nsppol
     call wan_spin(spin)%from_abiwan(dtfil%filabiwanin, spin, ebands%nsppol, keep_umats, dtfil%filnam_ds(4), comm)
     call wan_spin(spin)%print(units)
   end do
   if (dtfil%filgwanin /= ABI_NOFILE) then
     has_gwan = .True.
     ! TODO: Use Wannier to interpolate band energies on the dense k-mesh
     !call wan_interp_ebands(wan_spin, cryst, ebands, intp_kptrlatt, intp_nshiftk, intp_shiftk, dense_ebands, comm)
     !gstore%ebands => dense_ebands; gstore%ebands_owns_memory = .True.
   end if
 end if

 ! Set metadata.
 gstore%kibz => gstore%ebands%kptns
 gstore%kzone = dtset%gstore_kzone; gstore%qzone = dtset%gstore_qzone; gstore%kfilter = dtset%gstore_kfilter
 gstore%with_vk = dtset%gstore_with_vk; gstore%gmode = dtset%gstore_gmode

 ABI_CALLOC(gstore%erange_spin, (2, nsppol))
 gstore%erange_spin = dtset%gstore_erange(:, 1:nsppol)

 if (any(gstore%erange_spin /= zero)) then
   ABI_CHECK(gstore%kfilter == "none", sjoin("kfilter should be none when erange is used while it is:", gstore%kfilter))
   gstore%kfilter = "erange"
 end if

 if (gstore%kzone == "ibz" .and. gstore%qzone == "ibz") then
   ABI_ERROR("The combination kzone = 'ibz' and qzone = 'ibz' is not allowed!")
 end if

 ! TODO
 !gstore%kptrlatt(3, 3); gstore%kshift(3, 1); gstore%qptrlatt(3, 3); gstore%qshift(3, 1)

 ! Distribute spins, create indirect mapping to spin index and init %brange_spin
 ABI_CHECK_ILEQ(dtset%mband, ebands%mband, "dtset%mband > ebands%mband")
 call gstore%distribute_spins__(dtset%mband, dtset%gstore_brange, nproc_spin, comm_spin, comm)

 if (has_abiwan) then
   ! Here we set brange_spin to be consistent with the wannierization step.
   do spin=1,gstore%nsppol
     gstore%brange_spin(1, spin) = wan_spin(spin)%bmin
     gstore%brange_spin(2, spin) = wan_spin(spin)%bmax
   end do
 end if

 ! Free wan_spin
 do spin=1,ebands%nsppol
   call wan_spin(spin)%free()
 end do

 ! Define q-mesh: either from DVDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation)
 ! Save it in gstore for future reference.
 ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if
 gstore%ngqpt(:) = ngqpt(:)

 ! TODO: Should fix bz2ibz to use the same conventions as krank and listkk
 ! NB: only sigmaph seems to be using this optional argument

 ! Setup qIBZ, weights and BZ.
 ! Assume qptopt == kptopt unless value is specified in input
 qptrlatt = 0; qptrlatt(1, 1) = ngqpt(1); qptrlatt(2, 2) = ngqpt(2); qptrlatt(3, 3) = ngqpt(3)
 gstore%qptopt = ebands%kptopt; if (dtset%qptopt /= 0) gstore%qptopt = dtset%qptopt
 timrev_q = kpts_timrev_from_kptopt(gstore%qptopt)

 call wrtout(std_out, sjoin(" Generating q-mesh with ngqpt:", ltoa(ngqpt), " and qptopt:", itoa(gstore%qptopt)))
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, gstore%qptopt, my_nshiftq, my_shiftq, &
                             gstore%nqibz, gstore%qibz, gstore%wtq, gstore%nqbz, gstore%qbz)
                             !new_kptrlatt=gstore%qptrlatt, new_shiftk=gstore%qshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(qbz2ibz, (6, gstore%nqbz))
 call gstore%qrank_ibz%from_kptrlatt(gstore%nqibz, gstore%qibz, qptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", gstore%qptopt, cryst, gstore%qrank_ibz, gstore%nqbz, gstore%qbz, qbz2ibz) /= 0) then
   ABI_ERROR("Cannot map qBZ to IBZ!")
 end if

 ! Order qbz by stars and rearrange entries in qbz2ibz table.
 call kpts_pack_in_stars(gstore%nqbz, gstore%qbz, qbz2ibz)

 !call kpts_print_kmap(std_out, qibz, gstore%qbz, qbz2ibz)
 !do iq_bz=1,gstore%nqbz
 !  print *, "iq_bz -> iq_ibz", qbz2ibz(1, iq_bz), gstore%qbz(:, iq_bz)
 !end do

 call get_ibz2bz(gstore%nqibz, gstore%nqbz, qbz2ibz, qibz2bz, msg, ierr)
 ABI_CHECK(ierr == 0, sjoin("Something wrong in symmetry tables for q-points!", ch10, msg))

 ! Get full BZ associated to ebands
 call wrtout(std_out, sjoin(" Generating k-mesh with ngkpt:", ltoa(get_diag(ebands%kptrlatt)), " and kptopt:", itoa(ebands%kptopt)))
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
                             gstore%nkibz, kibz, wtk, gstore%nkbz, gstore%kbz) !, bz2ibz=bz2ibz)
                             !new_kptrlatt=gstore%kptrlatt, new_shiftk=gstore%kshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ! In principle kibz should be equal to ebands%kptns
 ABI_CHECK(gstore%nkibz == ebands%nkpt, "nkibz != ebands%nkpt")
 ABI_CHECK(all(abs(gstore%kibz - kibz) < tol12), "ebands%kibz != kibz")
 ABI_FREE(kibz)

 ! Note symrel and use_symrec=.False. in get_mapping.
 ! This means that this table can be used to symmetrize wavefunctions in cgtk_rotate.
 ! TODO This ambiguity should be removed. Change cgtk_rotate so that we can use the symrec convention.

 ABI_MALLOC(gstore%kbz2ibz, (6, gstore%nkbz))
 call gstore%krank_ibz%from_kptrlatt(gstore%nkibz, gstore%kibz, ebands%kptrlatt, compute_invrank=.False.)
 if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, gstore%nkbz, gstore%kbz, gstore%kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 ! TODO:
 ! Order kbz by stars and rearrange entries in kbz2ibz table.
 !call kpts_pack_in_stars(gstore%nkbz, kbz, gstore%kbz2ibz)

 call get_ibz2bz(gstore%nkibz, gstore%nkbz, gstore%kbz2ibz, kibz2bz, msg, ierr)
 ABI_CHECK(ierr == 0, sjoin("Something wrong in symmetry tables for k-points", ch10, msg))

 ! These tables are used to exclude q/k points
 ! We use the full BZ because this mask can be also used when points are restricted to the IBZ
 ! provided we convert from ik_ibz to ik_bz. Note that both arrays are initialized with zeros.
 ABI_ICALLOC(select_qbz_spin, (gstore%nqbz, nsppol))
 ABI_ICALLOC(select_kbz_spin, (gstore%nkbz, nsppol))

 select case (gstore%kzone)
 case ("ibz")
   do spin=1,nsppol
     do ik_ibz=1,gstore%nkibz
       ik_bz = kibz2bz(ik_ibz); select_kbz_spin(ik_bz, spin) = 1
     end do
   end do

 case ("bz")
   select_kbz_spin = 1

 case default
   ABI_ERROR(sjoin("Invalid kzone:", gstore%kzone))
 end select

 select case (gstore%qzone)
 case ("ibz")
   do spin=1,nsppol
     do iq_ibz=1, gstore%nqibz
       iq_bz = qibz2bz(iq_ibz); select_qbz_spin(iq_bz, spin) = 1
     end do
   end do

 case ("bz")
   select_qbz_spin = 1

 case default
   ABI_ERROR(sjoin("Invalid qzone:", gstore%qzone))
 end select

 ! Here we filter the electronic wavevectors k and recompute select_qbz_spin and select_kbz_spin according to kfilter.
 select case (gstore%kfilter)
 case ("none")
   continue

 case ("erange")
   call gstore%filter_erange__(gstore%qbz, qbz2ibz, qibz2bz, gstore%kbz, gstore%kibz, gstore%kbz2ibz, kibz2bz, &
                               select_qbz_spin, select_kbz_spin)

 case ("qprange")
   call gstore%filter_qprange__(dtset, gstore%qbz, qbz2ibz, qibz2bz, gstore%kbz, gstore%kibz, gstore%kbz2ibz, &
                                kibz2bz, select_qbz_spin, select_kbz_spin)

 case ("fs_tetra")
   ! Use the tetrahedron method to filter k- and k+q points on the FS in metals
   ! and define gstore%brange_spin automatically.
   call gstore%filter_fs_tetra__(gstore%qbz, qbz2ibz, qibz2bz, gstore%kbz, gstore%kibz, gstore%kbz2ibz, &
                                 kibz2bz, select_qbz_spin, select_kbz_spin)

 case default
   ABI_ERROR(sjoin("Invalid gstore%kfilter:", gstore%kfilter))
 end select

 ! Total number of k/q points for each spin after filtering (if any)
 ABI_MALLOC(gstore%glob_nk_spin, (nsppol))
 ABI_MALLOC(gstore%glob_nq_spin, (nsppol))
 gstore%glob_nk_spin(:) = count(select_kbz_spin > 0, dim=1)
 gstore%glob_nq_spin(:) = count(select_qbz_spin > 0, dim=1)

 ! We need another table mapping the global index in the gqk matrix to the q/k index in the BZ
 ! so that one can extract the symmetry tables computed above.
 ! Again this is needed as the global sizes of the gqk matrix
 ! is not necessarily equal to the size of the BZ/IBZ if we have filtered the wave vectors.

 max_nq = maxval(gstore%glob_nq_spin) ! Max dim over spin
 max_nk = maxval(gstore%glob_nk_spin)
 ABI_ICALLOC(qglob2bz, (max_nq, nsppol))
 ABI_ICALLOC(gstore%kglob2bz, (max_nk, nsppol))

 do spin=1,nsppol
   cnt = 0
   do iq_bz=1,gstore%nqbz
     if (select_qbz_spin(iq_bz, spin) /= 0) then
       cnt = cnt + 1; qglob2bz(cnt, spin) = iq_bz
     end if
   end do

   cnt = 0
   do ik_bz=1,gstore%nkbz
     if (select_kbz_spin(ik_bz, spin) /= 0) then
       cnt = cnt + 1; gstore%kglob2bz(cnt, spin) = ik_bz
     end if
   end do
 end do

 ! =============================================
 ! Initialize gqk basic dimensions and MPI grid
 ! =============================================
 call gstore%set_mpi_grid__(dtset%gstore_cplex, nproc_spin, comm_spin)
 call xmpi_comm_free(comm_spin)

 ! At this point, we have the Cartesian grid (one per spin if any)
 ! and we can finally allocate and distribute other arrays.
 ! Note with_cplex = 0 --> matrix elements are not allocated here
 with_cplex = 0
 if (has_gwan) with_cplex = 2
 call gstore%malloc__(with_cplex, max_nq, qglob2bz, max_nk, gstore%kglob2bz, qbz2ibz, gstore%kbz2ibz)

 ! Initialize GSTORE.nc file i.e. define dimensions and arrays
 ! Entries such as the e-ph matrix elements will be filled afterwards in gstore_compute.
 ! Master node defines dimensions and variables.

 write_gstore = .True.
 if (has_gwan) write_gstore = .False.

 if (my_rank == master .and. write_gstore) then
   NCF_CHECK(nctk_open_create(ncid, gstore%path, xmpi_comm_self))

   ! Write the abinit header with metadata, structure and occupancies.
   gstore_fform = fform_from_ext("GSTORE.nc")
   NCF_CHECK(wfk0_hdr%ncwrite(ncid, gstore_fform, spinat=dtset%spinat, nc_define=.True.))

   ! Add crystalline structure.
   NCF_CHECK(gstore%cryst%ncwrite(ncid))
   ! Add eigenvalues and occupations.
   NCF_CHECK(gstore%ebands%ncwrite(ncid))

   ! Write gstore dimensions
   ncerr = nctk_def_dims(ncid, [ &
      nctkdim_t("gstore_nkibz", gstore%nkibz), &
      nctkdim_t("gstore_nkbz", gstore%nkbz), &
      nctkdim_t("gstore_nqibz", gstore%nqibz), &
      nctkdim_t("gstore_nqbz", gstore%nqbz), &
      nctkdim_t("gstore_max_nq", max_nq), &
      nctkdim_t("gstore_max_nk", max_nk), &
      nctkdim_t("gstore_max_nb", maxval(gstore%brange_spin(2, :) - gstore%brange_spin(1, :) + 1) ), &
      nctkdim_t("natom", gstore%cryst%natom), &
      nctkdim_t("natom3", 3 * gstore%cryst%natom), &
      nctkdim_t("gstore_cplex", dtset%gstore_cplex) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "gstore_with_vk", "gstore_qptopt", "gstore_completed", &
     "gstore_use_lgk", "gstore_use_lgq" &
   ])
   NCF_CHECK(ncerr)
   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
   !NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("gstore_qibz", "dp", "three, gstore_nqibz"), &
     nctkarr_t("gstore_qbz", "dp", "three, gstore_nqbz"), &
     nctkarr_t("gstore_wtq", "dp", "gstore_nqibz"), &
     nctkarr_t("gstore_kbz", "dp", "three, gstore_nkbz"), &
     nctkarr_t("gstore_kzone", "c", "character_string_length"), &
     nctkarr_t("gstore_qzone", "c", "character_string_length"), &
     nctkarr_t("gstore_kfilter", "c", "character_string_length"), &
     nctkarr_t("gstore_gmode", "c", "character_string_length"), &
     nctkarr_t("gstore_gtype", "c", "character_string_length"), &
     nctkarr_t("gstore_wfk0_path", "c", "fnlen"), &
     nctkarr_t("gstore_brange_spin", "i", "two, number_of_spins"), &
     nctkarr_t("gstore_erange_spin", "dp", "two, number_of_spins"), &
     nctkarr_t("gstore_ngqpt", "i", "three"), &
     nctkarr_t("phfreqs_ibz", "dp", "natom3, gstore_nqibz"), &
     nctkarr_t("pheigvec_cart_ibz", "dp", "two, three, natom, natom3, gstore_nqibz"), &
     nctkarr_t("gstore_glob_nq_spin", "i", "number_of_spins"), &
     nctkarr_t("gstore_glob_nk_spin", "i", "number_of_spins"), &
     nctkarr_t("gstore_done_qbz_spin", "i", "gstore_nqbz, number_of_spins"), &
     nctkarr_t("gstore_kbz2ibz", "i", "six, gstore_nkbz"), &
     nctkarr_t("gstore_qbz2ibz", "i", "six, gstore_nqbz"), &
     nctkarr_t("gstore_qglob2bz", "i", "gstore_max_nq, number_of_spins"), &
     nctkarr_t("gstore_kglob2bz", "i", "gstore_max_nk, number_of_spins") &
   ])
   NCF_CHECK(ncerr)

   ! Internal table used to restart computation. Initialized with zeros.
   !  0 --> (ib_bz, spin) has not been computed.
   !  1 --> (iq_bz, spin) has been computed.
   ! In order to check if the whole generation is completed, one should test if "gstore_completed" == 1
   !NCF_CHECK(nf90_def_var_fill(ncid, vid("gstore_done_qbz_spin"), NF90_FILL, 0))

   ! Optional arrays
   if (allocated(gstore%delta_ef_kibz_spin)) then
     ncerr = nctk_def_arrays(ncid, &
       nctkarr_t("gstore_delta_ef_kibz_spin", "dp", "gstore_max_nb, gstore_nkibz, number_of_spins"))
   end if
   NCF_CHECK(ncerr)

   ! Write data
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_with_vk"), gstore%with_vk))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qptopt"), gstore%qptopt))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_use_lgk"), dtset%gstore_use_lgk))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_use_lgq"), dtset%gstore_use_lgq))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_completed"), 0))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kzone"), trim(gstore%kzone)))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qzone"), trim(gstore%qzone)))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kfilter"), trim(gstore%kfilter)))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_gmode"), trim(gstore%gmode)))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_gtype"), trim(gstore%gtype)))

   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qibz"), gstore%qibz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qbz"), gstore%qbz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_wtq"), gstore%wtq))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kbz"), gstore%kbz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_brange_spin"), gstore%brange_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_erange_spin"), gstore%erange_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_ngqpt"), gstore%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_glob_nq_spin"), gstore%glob_nq_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_glob_nk_spin"), gstore%glob_nk_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kbz2ibz"), gstore%kbz2ibz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qbz2ibz"), qbz2ibz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qglob2bz"), qglob2bz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kglob2bz"), gstore%kglob2bz))
   ! NB: kibz has been already written by ebands%ncwrite

   if (allocated(gstore%delta_ef_kibz_spin)) then
     NCF_CHECK(nf90_put_var(ncid, vid("gstore_delta_ef_kibz_spin"), gstore%delta_ef_kibz_spin))
   end if

   do spin=1,gstore%nsppol
     ! Create group for this spin.
     NCF_CHECK(nf90_def_grp(ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

     ! Dimensions in gqk_spin group
     ncerr = nctk_def_dims(spin_ncid, [ &
        nctkdim_t("nb", gstore%brange_spin(2, spin) - gstore%brange_spin(1, spin) + 1), &
        !nctkdim_t("nb_k", gstore%brange_spin(2, spin) - gstore%brange_spin(1, spin) + 1), &
        !nctkdim_t("nb_kq", gstore%brange_spin(2, spin) - gstore%brange_spin(1, spin) + 1), &
        nctkdim_t("glob_nk", gstore%glob_nk_spin(spin)), &
        nctkdim_t("glob_nq", gstore%glob_nq_spin(spin))  &
     ], defmode=.True.)
     NCF_CHECK(ncerr)

     ! Define scalars
     ncerr = nctk_def_iscalars(spin_ncid, [character(len=nctk_slen) :: "bstart"])
     NCF_CHECK(ncerr)

     ! arrays in gqk_spin group with the precious stuff. Note global dimensions.
     ncerr = nctk_def_arrays(spin_ncid, [ &
       nctkarr_t("gvals", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq") &
     ])
     NCF_CHECK(ncerr)

     ! Compress gvals to reduce size on disk.
     !NCF_CHECK(nf90_def_var_deflate(spin_ncid, vid_spin("gvals"), shuffle=1, deflate=1, deflate_level=5))
     ! IMPORTANT: Init gvals with zeros.

    ! Default value for entries in gvals, vk_cart_ibz and vkmat_cart_ibz arrays that have not been written.
    ! This can happen only if we have filtered wavevectors.
     gstore_fill_dp = zero
     if (gstore%kfilter == "none") gstore_fill_dp = -huge(one)
     NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("gvals"), NF90_FILL, gstore_fill_dp))

     ! In GWPT gvals is used for g^Sigma so we declare another array to store g^KS
     if (dtset%eph_task == 17) then
       ncerr = nctk_def_arrays(spin_ncid, [ &
         nctkarr_t("gvals_ks", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq") &
       ])
       NCF_CHECK(ncerr)
       NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("gvals_ks"), NF90_FILL, gstore_fill_dp))
     end if

     select case(gstore%with_vk)
     case (1)
       ! Diagonal terms only
       NCF_CHECK(nctk_def_arrays(spin_ncid, nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz")))
       NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("vk_cart_ibz"), NF90_FILL, gstore_fill_dp))
     case (2)
       ! Full (nb x nb) matrix.
       NCF_CHECK(nctk_def_arrays(spin_ncid, nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz")))
       NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("vkmat_cart_ibz"), NF90_FILL, gstore_fill_dp))
     end select

     ! Write (small) data
     NCF_CHECK(nctk_set_datamode(spin_ncid))
     NCF_CHECK(nf90_put_var(spin_ncid, vid_spin("bstart"), gstore%brange_spin(1, spin)))
   end do ! spin

   NCF_CHECK(nf90_close(ncid))
 end if ! master

 ! Make sure GSTORE.nc has been written by master
 call xmpi_barrier(gstore%comm)

 ABI_FREE(wtk)
 ABI_FREE(qibz2bz)
 ABI_FREE(kibz2bz)
 ABI_FREE(select_qbz_spin)
 ABI_FREE(select_kbz_spin)
 ABI_FREE(qglob2bz)
 ABI_FREE(qbz2ibz)

 call cwtime_report(" gstore_init:", cpu, wall, gflops)
 call pstat_proc%print(_PSTAT_ARGS_)

 if (has_gwan .and. with_cplex /= 0) then
   call wrtout(units, " Using Wannier interpolation to compute and store e-ph matrix elements ...", pre_newlines=1)

   do my_is=1,gstore%my_nspins
     spin = gstore%my_spins(my_is)
     associate (gqk => gstore%gqk(my_is), wan => gstore%gqk(my_is)%wan)

     ! Here we build gqk%wan for this spin from the ABIWAN.nc file
     ! and set the communicator for perturbations from gstore.
     call wan%from_abiwan(dtfil%filabiwanin, spin, ebands%nsppol, keep_umats, "", gqk%comm%value)
     wan%my_pert_start = gqk%my_pert_start; wan%my_npert = gqk%my_npert; wan%pert_comm => gqk%pert_comm

     ! Now load g(R_e, R_p) for this spin from GWAN.nc
     call wan%load_gwan(dtfil%filgwanin, gstore%cryst, spin, ebands%nsppol, gqk%comm) ! gqk%pert_comm,

     ! Interpolate my e-ph matrix elements.
     ! NOTE: the interpolated values are in the atomic representation and distributed over perts.
     ! Then one should take into account the change from atom and phonon representation.
     nq = 1
     ABI_MALLOC(intp_gatm, (wan%nwan, wan%nwan, wan%my_npert, nq))
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_qq, qpt)
       !call ifc%fourq(cryst, qq_ibz, phfr_qq, displ_cart_qibz, out_displ_red=displ_red_qibz, out_eigvec=pheigvec_qibz)
       do my_ik=1,gqk%my_nk
         kpt = gqk%my_kpts(:,my_ik); kq = kpt + qpt
         call wan%interp_eph_manyq(1, qpt, kpt, intp_gatm)
       end do ! my_ik
     end do ! my_iq
     ABI_FREE(intp_gatm)
     end associate
   end do ! my_is
 end if

contains
 integer function vid(var_name)
   character(len=*),intent(in) :: var_name
   vid = nctk_idname(ncid, var_name)
 end function vid
 integer function vid_spin(var_name)
   character(len=*),intent(in) :: var_name
   vid_spin = nctk_idname(spin_ncid, var_name)
 end function vid_spin

end subroutine gstore_init
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_distribute_spins
!! NAME
!! gstore_distribute_spins
!!
!! FUNCTION
!!  Distribute spins. Also create and return indirect mapping to spin index and init %brange_spin
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_distribute_spins(gstore, mband, gstore_brange, nproc_spin, comm_spin, comm)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: mband, comm, gstore_brange(2, gstore%nsppol)
 integer,intent(out) :: nproc_spin(gstore%nsppol), comm_spin(gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer :: spin, my_rank, ierr, color, nsppol, nprocs
!arrays
 integer :: buff_spin(2)
!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 nsppol = gstore%nsppol

 gstore%my_nspins = 0
 ABI_MALLOC(gstore%brange_spin, (2, nsppol))

 do spin=1,nsppol
   ! NB: If MPI_UNDEFINED is passed as the colour value, the subgroup in which the calling MPI process will be placed is MPI_COMM_NULL
   color = 1
   if (nsppol == 2 .and. nprocs > 1) then
     color = xmpi_undefined
     if (spin == 1 .and. my_rank <= (nprocs - 1) / 2) color = 1
     if (spin == 2 .and. my_rank >  (nprocs - 1) / 2) color = 1
   end if

   call xmpi_comm_split(comm, color, my_rank, comm_spin(spin), ierr)
   if (comm_spin(spin) /= xmpi_comm_null) then
     gstore%my_nspins = gstore%my_nspins + 1
     buff_spin(gstore%my_nspins) = spin
   end if

   nproc_spin(spin) = xmpi_comm_size(comm_spin(spin))

   gstore%brange_spin(:, spin) = [1, mband]
   if (all(gstore_brange /= 0)) gstore%brange_spin(:, spin) = gstore_brange(:, spin)

   ABI_CHECK_IRANGE(gstore%brange_spin(1, spin), 1, mband, "gstore_brange(1, spin)")
   ABI_CHECK_IRANGE(gstore%brange_spin(2, spin), 1, mband, "gstore_brange(2, spin)")
   ABI_CHECK(gstore%brange_spin(1, spin) <= gstore%brange_spin(2, spin), "brange_spin(1, spin) <= brange_spin(2, spin)")
 end do

 ABI_MALLOC(gstore%my_spins, (gstore%my_nspins))
 gstore%my_spins = buff_spin(1:gstore%my_nspins)
 ABI_MALLOC(gstore%gqk, (gstore%my_nspins))

end subroutine gstore_distribute_spins
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_set_mpi_grid__
!! NAME
!! gstore_set_mpi_grid__
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_set_mpi_grid__(gstore, gstore_cplex, nproc_spin, comm_spin)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: gstore_cplex
 integer,intent(in) :: nproc_spin(gstore%nsppol)
 integer,intent(inout) :: comm_spin(gstore%nsppol)
!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: spin, my_is, np, my_rank, ierr, npp_bz, units(2)
 type(gqk_t),pointer :: gqk
 character(len=5000) :: msg
 character(len=10) :: order
 integer :: comm_cart, me_cart, dims(ndims)
 logical :: reorder, periods(ndims), keepdim(ndims)
 character(len=10) :: priority
!----------------------------------------------------------------------

 units = [std_out, ab_out]

 associate (dtset => gstore%dtset)
 my_rank = xmpi_comm_rank(gstore%comm)

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is)
   gqk%spin = spin; gqk%natom3 = 3 * gstore%cryst%natom; gqk%cplex = gstore_cplex
   ABI_CHECK_IRANGE(gqk%cplex, 1, 2, "gstore_cplex")

   ! Compute bstart and band size for this spin.
   gqk%bstart = gstore%brange_spin(1, spin)
   gqk%bstop = gstore%brange_spin(2, spin)
   gqk%nb = gstore%brange_spin(2, spin) - gstore%brange_spin(1, spin) + 1

   ! FIXME: for the time being, a direct copy of %nb
   gqk%nb_k = gqk%nb; gqk%nb_kq = gqk%nb
   gqk%bstart_k = gqk%bstart; gqk%bstop_k = gqk%bstop
   gqk%bstart_kq = gqk%bstart; gqk%bstop_kq = gqk%bstop

   ! Store global shape of the q/k matrix for this spin.
   gqk%glob_nq = gstore%glob_nq_spin(spin)
   gqk%glob_nk = gstore%glob_nk_spin(spin)
 end do

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is)

   ! Init for sequential execution.
   gqk%my_npert = gqk%natom3
   gqk%qpt_comm%nproc = 1; gqk%kpt_comm%nproc = 1; gqk%pert_comm%nproc = 1; gqk%band_comm%nproc = 1
   ! NB: The communicators below are only used in GWPT.
   gqk%bsum_comm%nproc = 1; gqk%pp_sum_comm%nproc = 1

   np = nproc_spin(spin)

   if (dtset%eph_task /= 17) then
     ! =============================
     ! For all eph_tasks except GWPT
     ! =============================

     if (any(dtset%eph_np_pqbks /= 0)) then
       ! Use parameters from input file. Need to perform sanity check though.
       gqk%pert_comm%nproc = dtset%eph_np_pqbks(1)
       gqk%qpt_comm%nproc  = dtset%eph_np_pqbks(2)
       gqk%band_comm%nproc  = dtset%eph_np_pqbks(3)
       ABI_CHECK(dtset%eph_np_pqbks(3) == 1, "Band parallelism not yet implemented in gstore")
       gqk%kpt_comm%nproc = dtset%eph_np_pqbks(4)
       !gqk%spin_comm%nproc = dtset%eph_np_pqbks(5)
       gqk%my_npert = gqk%natom3 / gqk%pert_comm%nproc
       ABI_CHECK(gqk%my_npert > 0, "pert_comm_nproc cannot be greater than 3*natom.")
       ABI_CHECK(mod(gqk%natom3, gqk%pert_comm%nproc) == 0, "pert_comm_nproc must divide 3*natom.")

     else
       ! Automatic grid generation (hopefully smart)
       ! Keep in mind that in gstore_compute, the first loop is over q-points
       ! in order to reduce the number of interpolations of the DFPT potentials in q-space
       ! hence the q-point parallelism is expected to be more efficient.
       ! On the other hand, the k-point parallelism and the perturbation parallelism
       ! allow one to reduce the memory requirements associated to the wavefunctions (kpt) and
       ! the scattering potentials in the supercell (perturbations).
       ! Here we try to optimize performance but it's clear that for large systems the user
       ! should specify dtset%eph_np_pqbks in the input file.

       select case (dtset%eph_task)
       case (11)
         priority = "qk"
       case (12, -12)
         priority = "q"
       case (13)
         priority = "q"
       case (14, 17)
         priority = "kq"
       case (24)
         priority = "qk"
       case default
         ABI_ERROR(sjoin("Please register default priority for eph_task:", itoa(dtset%eph_task)))
       end select

       select case (priority)
       case ("q")
         order = "1"
       case ("k")
         order = "2"
       case ("kq")
         order = "21"
       case ("qk")
         order = "12"
       case default
         ABI_ERROR(sjoin("Wrong priority:", priority))
       end select

       if (gqk%glob_nk == 1) order = "21"
       if (gqk%glob_nq == 1) order = "12"
       call xmpi_distrib_2d(np, order, gqk%glob_nq, gqk%glob_nk, gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, ierr)
       ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np), " with priority: ", priority))
     end if

   else if (dtset%eph_task == 17) then
     ! =========================
     ! MPI distribution for GWPT
     ! =========================

     if (any(dtset%gwpt_np_wpqbks /= 0)) then
       ! Use parameters from input file. Need to perform sanity check though.
       gqk%pp_sum_comm%nproc  = dtset%gwpt_np_wpqbks(1)
       gqk%pert_comm%nproc = dtset%gwpt_np_wpqbks(2)
       gqk%qpt_comm%nproc  = dtset%gwpt_np_wpqbks(3)
       gqk%bsum_comm%nproc  = dtset%gwpt_np_wpqbks(4)
       gqk%kpt_comm%nproc = dtset%gwpt_np_wpqbks(5)
       !gqk%spin_comm%nproc = dtset%gwpt_np_wpqbks(6)
       gqk%my_npert = gqk%natom3 / gqk%pert_comm%nproc
       ABI_CHECK(gqk%my_npert > 0, "pert_comm_nproc cannot be greater than 3*natom.")
       ABI_CHECK(mod(gqk%natom3, gqk%pert_comm%nproc) == 0, "pert_comm_nproc must divide 3*natom.")

     else
       ! Automatic grid generation for GWPT (hopefully smart)

       !if (gqk%glob_nk == 1 .and. gqk%glob_nq == 1) then
         ! This may happen in GWPT when only of e-ph matrix element is wanted.
         ! Here we activate the parallelism over pp_sum and perturbations.
         ! In principle we should distributed the
         npp_bz = product(get_diag(gstore%dtset%kptrlatt))

         !call kmesh%init(gstore%cryst, gstore%dtset%nkibz, gstore%dtset%kptns, dtset%kptopt)
         !call find_qmesh(qmesh, gstore%cryst, kmesh)
         !npp_bz = qmesh%nbz
         !call kmesh%free(); call qmesh%free()

         order = "12"
         !order = "21"
         call xmpi_distrib_2d(np, order, npp_bz, gqk%natom3, gqk%pp_sum_comm%nproc, gqk%pert_comm%nproc, ierr)
         !call xmpi_distrib_2d(np, order, gqk%natom3, gstore%dtset%mband, gqk%pert_comm%nproc, gqk%bsum_comm%nproc, ierr)
         ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np), " with priority: ", priority))
       !end if
     end if

   else
     ABI_ERROR(sjoin("Invalid eph_task", itoa(dtset%eph_task)))
   end if

   ! Consistency check.
   if (gqk%pert_comm%nproc * gqk%qpt_comm%nproc * gqk%kpt_comm%nproc * gqk%band_comm%nproc * &
       gqk%bsum_comm%nproc * gqk%pp_sum_comm%nproc /= nproc_spin(spin)) then
     write(msg, "(a,i0,3a, 7(a,1x,i0))") &
       "Cannot create Cartesian grid with total nproc: ", nproc_spin(spin), ch10, &
       "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
       "qpt_nproc (", gqk%qpt_comm%nproc, ") x kpt_nproc (", gqk%kpt_comm%nproc, ") x pert_nproc", gqk%pert_comm%nproc, &
       "x band_nproc (", gqk%band_comm%nproc, "x bsum_nproc (", gqk%bsum_comm%nproc, ") x psum_nproc (", gqk%pp_sum_comm%nproc, &
       ") != ", nproc_spin(spin)
     ABI_ERROR(msg)
   end if

 end do ! my_is

 ! For each spin treated by this rank, create Cartesian communicator of rank ndims.
 periods(:) = .False.; reorder = .False.

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is)

   ! TODO: Should change order for GWPT
   dims = [gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, gqk%pert_comm%nproc, gqk%band_comm%nproc, &
           gqk%bsum_comm%nproc, gqk%pp_sum_comm%nproc]

   ! Note comm_spin(spin)
   gqk%comm = xcomm_from_mpi_int(comm_spin(spin))
   gqk%coords_qkpb_sumbp = 0

#ifdef HAVE_MPI
   call MPI_CART_CREATE(gqk%comm, ndims, dims, periods, reorder, comm_cart, ierr)
   ! Find the index and coordinates of the current processor
   call MPI_COMM_RANK(comm_cart, me_cart, ierr)
   call MPI_CART_COORDS(comm_cart, me_cart, ndims, gqk%coords_qkpb_sumbp, ierr)

   ! Communicator for q-points in g(k,q)
   keepdim = .False.; keepdim(1) = .True.; call gqk%qpt_comm%from_cart_sub(comm_cart, keepdim)
   ! Communicator for k-points in g(k,q)
   keepdim = .False.; keepdim(2) = .True.; call gqk%kpt_comm%from_cart_sub(comm_cart, keepdim)
   ! Communicator for the (qpt, kpt) 2D grid
   keepdim = .False.; keepdim(1) = .True.; keepdim(2) = .True.; call gqk%qpt_kpt_comm%from_cart_sub(comm_cart, keepdim)
   ! Communicator for perturbations in g(k,q)
   keepdim = .False.; keepdim(3) = .True.; call gqk%pert_comm%from_cart_sub(comm_cart, keepdim)
   ! 2d Communicator for the (qpt, pert) grid
   keepdim = .False.; keepdim(1) = .True.; keepdim(3) = .True.; call gqk%qpt_pert_comm%from_cart_sub(comm_cart, keepdim)
   ! Communicator for band in g(k,q)
   keepdim = .False.; keepdim(4) = .True.; call gqk%band_comm%from_cart_sub(comm_cart, keepdim)
   ! Communicator for bsum (GWPT mode)
   keepdim = .False.; keepdim(5) = .True.; call gqk%bsum_comm%from_cart_sub(comm_cart, keepdim)
   ! Communicator for pp_sum (GWPT mode)
   keepdim = .False.; keepdim(6) = .True.; call gqk%pp_sum_comm%from_cart_sub(comm_cart, keepdim)
   ! 2d Communicator for the (pert, pp_sum) 2D grid
   keepdim = .False.; keepdim(3) = .True.; keepdim(6) = .True.; call gqk%pert_ppsum_comm%from_cart_sub(comm_cart, keepdim)
   ! 3d Communicator for the (pert, pp_sum, band_sum) 3D grid (GWPT mode)
   keepdim = .False.; keepdim(3) = .True.; keepdim(5) = .True.; keepdim(6) = .True.
   call gqk%pert_ppsum_bsum_comm%from_cart_sub(comm_cart, keepdim)
   call xmpi_comm_free(comm_cart)
#endif

   ! Distribute perturbations inside pert_comm using block distribution.
   call xmpi_split_block(gqk%natom3, gqk%pert_comm%value, gqk%my_npert, gqk%my_pertcases)
   gqk%my_pert_start = gqk%my_pertcases(1)

   call wrtout(units, sjoin("P qpt_comm can use shmem:", yesno(gqk%qpt_comm%can_use_shmem())))
   call wrtout(units, sjoin("P kpt_comm can use shmem:", yesno(gqk%kpt_comm%can_use_shmem())))
   call wrtout(units, sjoin("P bsum_comm can use shmem:", yesno(gqk%bsum_comm%can_use_shmem())))
   call wrtout(units, sjoin("P pp_sum_comm can use shmem:", yesno(gqk%pp_sum_comm%can_use_shmem())), newlines=1)
 end do ! my_is

 if (my_rank == master) call gstore%print([std_out, ab_out])

 end associate

end subroutine gstore_set_mpi_grid__
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_print
!! NAME
!! gstore_print
!!
!! FUNCTION
!!  Print info on the gstore object.
!!
!! INPUTS
!!  units: Unit numbers.
!!  [header]: header string
!!  [prtvol]: Verbosity level.
!!
!! SOURCE

subroutine gstore_print(gstore, units, header, prtvol)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 integer,intent(in) :: units(:)
 character(len=*),optional,intent(in) :: header
 integer,optional,intent(in) :: prtvol

!Local variables ------------------------------
 integer,parameter :: max_nk=10
 integer :: my_is, my_prtvol, ik_calc, ik_bz, ik_ibz !, iq_calc, iq_bz, iq_ibz
 character(len=500) :: msg
!----------------------------------------------------------------------

 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 if (present(header)) then
   call wrtout(units, header)
 else
   call wrtout(units, " === Gstore parameters ===")
 end if

 call wrtout(units, sjoin(" kzone:", gstore%kzone))
 call wrtout(units, sjoin(" kfilter:", gstore%kfilter))
 call wrtout(units, sjoin(" nkibz:", itoa(gstore%nkibz)))
 call wrtout(units, sjoin(" nkbz:", itoa(gstore%nkbz)))
 call wrtout(units, sjoin(" glob_nk_spin:", ltoa(gstore%glob_nk_spin)))
 call wrtout(units, sjoin(" qzone:", gstore%qzone))
 call wrtout(units, sjoin(" nqibz:", itoa(gstore%nqibz)))
 call wrtout(units, sjoin(" nqbz:", itoa(gstore%nqbz)))
 call wrtout(units, sjoin(" glob_nq_spin:", ltoa(gstore%glob_nq_spin)))
 call wrtout(units, sjoin(" kptopt:", itoa(gstore%ebands%kptopt)))
 call wrtout(units, sjoin(" qptopt:", itoa(gstore%qptopt)))
 call wrtout(units, sjoin(" has_used_lgk:", itoa(gstore%has_used_lgk)))
 call wrtout(units, sjoin(" has_used_lgq:", itoa(gstore%has_used_lgq)))
 call wrtout(units, sjoin(" with_vk:", itoa(gstore%with_vk)))

 do my_is=1,gstore%my_nspins
   associate (spin => gstore%my_spins(my_is), gqk => gstore%gqk(my_is))
   call wrtout(units, sjoin(" gqk_cplex:", itoa(gqk%cplex)), pre_newlines=1)
   call wrtout(units, sjoin(" gqk_bstart_k:", itoa(gqk%bstart_k)))
   call wrtout(units, sjoin(" gqk_bstart_kq:", itoa(gqk%bstart_kq)))
   call wrtout(units, sjoin(" gqk_bstop_k:", itoa(gqk%bstop_k)))
   call wrtout(units, sjoin(" gqk_bstop_kq:", itoa(gqk%bstop_kq)))
   call wrtout(units, sjoin(" gqk_nb_k:", itoa(gqk%nb_k)))
   call wrtout(units, sjoin(" gqk_nb_kq:", itoa(gqk%nb_kq)))
   call wrtout(units, sjoin(" gqk_my_npert:", itoa(gqk%my_npert)))
   call wrtout(units, sjoin("P gqk_my_nk:", itoa(gqk%my_nk)))
   call wrtout(units, sjoin("P gqk_my_nq:", itoa(gqk%my_nq)))
   !tot_num_g = one * qqk%nb_kq * gqk%nb_k * gqk%glob_nk * gqk%glob_nq * gqk%natom3
   call wrtout(units, sjoin(ch10, " === MPI distribution ==="))
   call wrtout(units, sjoin("P Number of CPUs for parallelism over perturbations: ", itoa(gqk%pert_comm%nproc)))
   call wrtout(units, sjoin("P Number of perturbations treated by this CPU: ",  itoa(gqk%my_npert)))
   call wrtout(units, sjoin("P Number of CPUs for parallelism over q-points: ", itoa(gqk%qpt_comm%nproc)))
   call wrtout(units, sjoin("P Number of CPUs for parallelism over k-points: ", itoa(gqk%kpt_comm%nproc)))
   ! This only for GWPT.
   if (gqk%bsum_comm%nproc /= 1) then
     call wrtout(units, sjoin("P Number of CPUs for parallelism over band summation: ", itoa(gqk%bsum_comm%nproc)))
   end if
   if (gqk%pp_sum_comm%nproc /= 1) then
     call wrtout(units, sjoin("P Number of CPUs for parallelism over wavevector summation: ", itoa(gqk%pp_sum_comm%nproc)))
   end if

   ! Print q-points
   call wrtout(units, " k-points included in gstore:")
   do ik_calc=1,gqk%glob_nk
     ik_bz = gstore%kglob2bz(ik_calc, spin)
     ik_ibz = gstore%kbz2ibz(1, ik_bz)
     call wrtout(units, sjoin(itoa(ik_calc), ":", ktoa(gstore%kbz(:, ik_bz))))
     if (ik_calc > max_nk .and. my_prtvol == 0) then
       call wrtout(units, sjoin(" Max", itoa(max_nk), " k-points will be written. Use prtvol > 0 to print all of them."))
       exit
     end if
   end do

   ! Print q-points
   !call wrtout(units, " q-points included in gstore:")
   !do iq_calc=1,gqk%glob_nq
   !  iq_bz = gstore%qglob2bz(iq_calc, spin)
   !  iq_ibz = gstore%qbz2ibz(1, iq_bz)
   !  call wrtout(units, sjoin(itoa(iq_calc), ":", ktoa(gstore%qbz(:, iq_bz))))
   !  if (iq_calc > max_nk .and. my_prtvol == 0) then
   !    call wrtout(units, sjoin(" Max", itoa(max_nk), " q-points will be written. Use prtvol > 0 to print all of them."))
   !    exit
   !  end if
   !end do

   ! Print memory
   if (allocated(gqk%my_g2)) then
     write(msg,'(a,f8.1,a)')'- Local memory allocated for |g|^2 array: ',ABI_MEM_MB(gqk%my_g2),' [Mb] <<< MEM'
     call wrtout(units, msg)
   end if
   if  (allocated(gqk%my_g)) then
     write(msg,'(a,f8.1,a)')'- Local memory allocated for g array: ',ABI_MEM_MB(gqk%my_g),' [Mb] <<< MEM'
     call wrtout(units, msg)
   end if
   if  (allocated(gqk%my_gq0nm_atm)) then
     write(msg,'(a,f8.1,a)')'- Local memory allocated for g0nm_atm array: ',ABI_MEM_MB(gqk%my_gq0nm_atm),' [Mb] <<< MEM'
     call wrtout(units, msg)
   end if
   if (allocated(gqk%vk_cart_ibz)) then
     write(msg,'(a,f8.1,a)')'- Local memory allocated for vnk_cart_ibz: ',ABI_MEM_MB(gqk%vk_cart_ibz),' [Mb] <<< MEM'
     call wrtout(units, msg)
   end if
   if (allocated(gqk%vkmat_cart_ibz)) then
     write(msg,'(a,f8.1,a)')'- Local memory allocated for vkmat_cart_ibz: ',ABI_MEM_MB(gqk%vkmat_cart_ibz),' [Mb] <<< MEM'
     call wrtout(units, msg)
   end if
   end associate
 end do

 if (my_prtvol > 0) then
   call gstore%ebands%print(units, header="Electron bands in GSTORE", prtvol=my_prtvol)
 end if

end subroutine gstore_print
!!***

!!****f* m_gstore/gstore_check_little_group
!! NAME
!! gstore_check_little_group
!!
!! FUNCTION
!!  Check consistency between little group options from file and from input.
!!
!! INPUTS
!!
!! SOURCE

integer function gstore_check_little_group(gstore, dtset, msg) result(ierr)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 type(dataset_type),intent(in) :: dtset
 character(len=*),intent(out) :: msg
!----------------------------------------------------------------------

 ierr = 0; msg = ""

 if (dtset%gstore_use_lgk /= 0) then
   ! Cannot use IBZ_k if we have used IBZ_q.
   if (gstore%has_used_lgq /= 0) then
     msg = sjoin("Input var gstore_use_lgq: ", itoa(dtset%gstore_use_lgq), ", but GSTORE file has:", itoa(gstore%has_used_lgq))
     ABI_ERROR_NOSTOP(msg, ierr)
   end if
   ! Using IBZ_k when GSTORE has full BZ if OK but inefficient.
   if (gstore%has_used_lgk == 0) then
     msg = sjoin("Input var gstore_use_lgk: ", itoa(dtset%gstore_use_lgk), ", but GSTORE file has:", itoa(gstore%has_used_lgk))
     ABI_COMMENT(msg)
   end if
 end if

 if (dtset%gstore_use_lgq /= 0) then
   ! Cannot use IBZ_q if we have used IBZ_k.
   if (gstore%has_used_lgk /= 0) then
     msg = sjoin("Input var gstore_use_lgk: ", itoa(dtset%gstore_use_lgk), ", but GSTORE file has:", itoa(gstore%has_used_lgk))
     ABI_ERROR_NOSTOP(msg, ierr)
   end if
   ! Using IBZ_q when GSTORE has full BZ if OK but inefficient.
   if (gstore%has_used_lgq == 0) then
     msg = sjoin("Input var gstore_use_lgq: ", itoa(dtset%gstore_use_lgq), ", but GSTORE file has:", itoa(gstore%has_used_lgq))
     ABI_COMMENT(msg)
   end if
 end if

end function gstore_check_little_group
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_malloc__
!! NAME
!! gstore_malloc__
!!
!! FUNCTION
!! Allocate local buffers once the MPI grid has been initialized.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_malloc__(gstore, with_cplex, max_nq, qglob2bz, max_nk, kglob2bz, qbz2ibz, kbz2ibz)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: with_cplex, max_nq, max_nk
 integer,intent(in) :: qglob2bz(max_nq, gstore%nsppol), kglob2bz(max_nk, gstore%nsppol)
 integer,intent(in) :: qbz2ibz(6, gstore%nqbz), kbz2ibz(6, gstore%nkbz)

!Local variables-------------------------------
!scalars
 integer :: my_is, ierr, my_iq, my_ik, iq_glob, iq_bz, ik_glob, ik_bz
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3), nb_k, nb_kq
 logical :: isirr_k
 real(dp) :: mem_mb
!----------------------------------------------------------------------

 !gstore%max_nb = maxval(gstore%brange_spin(2, :) - gstore%brange_spin(1, :) + 1)

 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), spin => gstore%my_spins(my_is))
   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq

   ! Split q-points and transfer symmetry tables.
   ! Note that glob_nq and glob_nk does not necessarily correspond to the size of the BZ
   ! First of all we have to consider kzone
   ! Even if kzone == "bz" we may have filtered the wavevectors e.g. Fermi surface.
   call xmpi_split_block(gqk%glob_nq, gqk%qpt_comm%value, gqk%my_nq, gqk%my_q2glob)
   ABI_CHECK(gqk%my_nq > 0, sjoin("glob_nq:", itoa(gqk%glob_nq), ", qpt_comm%nproc:", itoa(gqk%qpt_comm%nproc), " => my_nq == 0"))
   gqk%my_qstart = gqk%my_q2glob(1)

   ABI_MALLOC(gqk%my_q2ibz, (6, gqk%my_nq))
   ABI_MALLOC(gqk%my_q2bz, (gqk%my_nq))

   do my_iq=1,gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1
     iq_bz = qglob2bz(iq_glob, spin)
     gqk%my_q2ibz(:, my_iq) = qbz2ibz(:, iq_bz)
     gqk%my_q2bz(my_iq) = iq_bz
   end do

   ! Split k-points and transfer symmetry tables
   call xmpi_split_block(gqk%glob_nk, gqk%kpt_comm%value, gqk%my_nk, gqk%my_k2glob)
   ABI_CHECK(gqk%my_nk > 0, sjoin("glob_nk:", itoa(gqk%glob_nk), ", kpt_comm%nproc:", itoa(gqk%kpt_comm%nproc), " => my_nk == 0"))
   gqk%my_kstart = gqk%my_k2glob(1)

   ABI_MALLOC(gqk%my_k2ibz, (6, gqk%my_nk))
   !ABI_MALLOC(gqk%my_k2bz, (gqk%my_nk))
   ABI_MALLOC(gqk%my_kpts, (3, gqk%my_nk))
   ABI_MALLOC(gqk%my_wtk, (gqk%my_nk))

   do my_ik=1,gqk%my_nk
     ik_glob = my_ik + gqk%my_kstart - 1
     ik_bz = kglob2bz(ik_glob, spin)
     gqk%my_k2ibz(:, my_ik) = kbz2ibz(:, ik_bz)
     !gqk%my_k2bz(my_ik) = ik_bz

     ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
     trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5, my_ik)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     tsign_k = 1; if (trev_k == 1) tsign_k = -1

     ! Note symrel^T convention for k
     gqk%my_kpts(:, my_ik) = tsign_k * matmul(transpose(gstore%cryst%symrel(:,:,isym_k)), gstore%kibz(:, ik_ibz)) + g0_k

     select case(gstore%kzone)
     case ("ibz")
       gqk%my_wtk(my_ik) = gstore%ebands%wtk(ik_ibz)
     case ("bz")
       gqk%my_wtk(my_ik) = one / gstore%nkbz
     end select
   end do

   ! Allocate storage for MPI-distributed e-ph matrix elements.
   if (with_cplex > 0) then
     mem_mb = with_cplex * gqk%my_npert * nb_kq * nb_k * gqk%my_nq * gqk%my_nk * eight * b2Mb
     call wrtout(std_out, sjoin(" Local memory for e-ph matrix elements:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))

     ! The initialization with zero is important as not all the g are computed when we filter in k-space.
     ! Abinit postprocessing tools will operate of the full my_g array
     ! and we don't want to trigger floating point exceptions.
     select case (with_cplex)
     case (1)
       ABI_MALLOC_OR_DIE(gqk%my_g2, (gqk%my_npert, nb_kq, gqk%my_nq, nb_k, gqk%my_nk), ierr)
       gqk%my_g2 = zero
     case (2)
       ABI_MALLOC_OR_DIE(gqk%my_g, (gqk%my_npert, nb_kq, gqk%my_nq, nb_k, gqk%my_nk), ierr)
       gqk%my_g = zero
     case default
       ABI_ERROR(sjoin("Wrong with_cplex:", itoa(with_cplex)))
     end select

     ! Allocate storage for MPI-distributed dH/dk matrix elements.
     if (any(gstore%with_vk == [1, 2])) then
       mem_mb = 3 * gqk%nb * gqk%my_nk * eight * b2Mb
       call wrtout(std_out, sjoin(" Memory for diagonal vk_cart_ibz:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))
       ABI_MALLOC_OR_DIE(gqk%vk_cart_ibz, (3, gqk%nb, gstore%nkibz), ierr)
       gqk%vk_cart_ibz = zero
     end if

     if (gstore%with_vk == 2) then
       mem_mb = two * 3 * gqk%nb ** 2 * gqk%my_nk * eight * b2Mb
       call wrtout(std_out, sjoin(" Memory for vkmat_cart_ibz:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))
       ABI_MALLOC_OR_DIE(gqk%vkmat_cart_ibz, (2, 3, gqk%nb, gqk%nb, gstore%nkibz), ierr)
       gqk%vkmat_cart_ibz = zero
     end if
   end if

   end associate
 end do ! my_is

end subroutine gstore_malloc__
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_filter_fs_tetra__
!! NAME
!! gstore_filter_fs_tetra__
!!
!! FUNCTION
!!  Compute delta(e_k - e_F) with the tetrahedron method. Use weights to filter k-points.
!!  Include only those q-points such that there exists at least one k on the FS with k + q on the FS.
!!  Also, store and precompute gstore%delta_ef_kibz_spin(max_nb, gstore%nkibz, nsppol)
!!  to be used to filter inside gstore%compute.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_filter_fs_tetra__(gstore, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
                                    select_qbz_spin, select_kbz_spin)

!Arguments ------------------------------------
!scalars
 class(gstore_t),intent(inout) :: gstore
 real(dp),intent(in) :: qbz(3, gstore%nqbz)
 real(dp),target,intent(in) :: kbz(3, gstore%nqbz), kibz(3, gstore%nkibz)
 integer,intent(in) :: qbz2ibz(6,gstore%nqbz), qibz2bz(gstore%nqibz)
 integer,intent(in) :: kbz2ibz(6,gstore%nkbz), kibz2bz(gstore%nkibz)
 integer,intent(out) :: select_qbz_spin(gstore%nqbz, gstore%nsppol)
 integer,intent(out) :: select_kbz_spin(gstore%nkbz, gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: tetra_opt0 = 0
 integer :: nsppol, ierr, cnt, spin, band, ib, ii, max_nb, all_nproc, my_rank, comm
 integer :: ik_bz, ik_ibz, iflag, nk_in_star
 real(dp) :: max_occ
 character(len=80) :: error_string
 type(htetra_t) :: ktetra
!arrays
 integer,allocatable :: indkk(:), kstar_bz_inds(:)
 real(dp):: rlatt(3,3), klatt(3,3), delta_theta_ef(2) ! qpt(3),
 real(dp),allocatable :: eig_ibz(:)
!----------------------------------------------------------------------

 associate (cryst => gstore%cryst, ebands => gstore%ebands)

 comm = gstore%comm; all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 nsppol = gstore%nsppol

 ! Use the tetrahedron method to filter k- and k+q points on the FS in metals
 ! and define gstore%brange_spin automatically.
 call wrtout(std_out, sjoin(" Filtering k-points using:", gstore%kfilter))

 ! NB: here we recompute brange_spin
 call ebands%get_bands_e0(ebands%fermie, gstore%brange_spin, ierr)
 ABI_CHECK(ierr == 0, "Error in ebands_get_bands_e0")

 ABI_MALLOC(indkk, (gstore%nkbz))
 indkk(:) = kbz2ibz(1, :)

 rlatt = ebands%kptrlatt; call matr3inv(rlatt, klatt)
 call ktetra%init(indkk, gstore%cryst%gprimd, klatt, kbz, gstore%nkbz, gstore%kibz, gstore%nkibz, &
                  ierr, error_string, gstore%comm)
 ABI_CHECK(ierr == 0, error_string)

 ABI_MALLOC(eig_ibz, (gstore%nkibz))
 max_occ = two / (ebands%nspinor * nsppol)
 select_kbz_spin = 0

 max_nb = maxval(gstore%brange_spin(2, :) - gstore%brange_spin(1, :) + 1)
 ABI_CALLOC(gstore%delta_ef_kibz_spin, (max_nb, gstore%nkibz, gstore%nsppol))

 cnt = 0
 do spin=1,nsppol
   do band=gstore%brange_spin(1, spin), gstore%brange_spin(2, spin)
     ib = band - gstore%brange_spin(1, spin) + 1
     eig_ibz = ebands%eig(band, :, spin)

     do ik_ibz=1,gstore%nkibz
       cnt = cnt + 1; if (mod(cnt, all_nproc) /= my_rank) cycle ! MPI parallelism inside comm

       call ktetra%get_onewk_wvals(ik_ibz, tetra_opt0, 1, [ebands%fermie], max_occ, &
                                   gstore%nkibz, eig_ibz, delta_theta_ef)

       gstore%delta_ef_kibz_spin(ib, ik_ibz, spin) = delta_theta_ef(1)

       iflag = merge(1, 0, abs(delta_theta_ef(1)) > zero)

       ! Use iflag to filter k-points.
       select case (gstore%kzone)
       case ("ibz")
         ik_bz = kibz2bz(ik_ibz)
         select_kbz_spin(ik_bz, spin) = select_kbz_spin(ik_bz, spin) + iflag

       case ("bz")
         call star_from_ibz_idx(ik_ibz, gstore%nkbz, kbz2ibz, nk_in_star, kstar_bz_inds)
         ABI_CHECK(nk_in_star > 0, "Something wrong in star_from_ibz_idx")
         do ii=1,nk_in_star
           ik_bz = kstar_bz_inds(ii)
           select_kbz_spin(ik_bz, spin) = select_kbz_spin(ik_bz, spin) + iflag
         end do
         ABI_FREE(kstar_bz_inds)
       end select
     end do

   end do
 end do

 !call ktetra%print(std_out)
 call xmpi_sum(select_kbz_spin, comm, ierr)
 call xmpi_sum(gstore%delta_ef_kibz_spin, comm, ierr)

 ! Now the tricky part as we want to remove q-points that
 ! do not lead to any scattering process between two states on the FS
 ! Remember that k+q is always a sub-mesh of the input ebands k-mesh.

 call recompute_select_qbz_spin(gstore, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
                                select_kbz_spin, select_qbz_spin)

 ABI_FREE(eig_ibz)
 ABI_FREE(indkk)
 call ktetra%free()
 end associate

end subroutine gstore_filter_fs_tetra__
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_filter_erange__
!! NAME
!! gstore_filter_erange__
!!
!! FUNCTION
!! Filter k-points and q-points according to an energy range.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_filter_erange__(gstore, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
                                  select_qbz_spin, select_kbz_spin)

!Arguments ------------------------------------
!scalars
 class(gstore_t),intent(inout) :: gstore
 real(dp),intent(in) :: qbz(3, gstore%nqbz)
 real(dp),target,intent(in) :: kbz(3, gstore%nqbz), kibz(3, gstore%nkibz)
 integer,intent(in) :: qbz2ibz(6,gstore%nqbz), qibz2bz(gstore%nqibz)
 integer,intent(in) :: kbz2ibz(6,gstore%nkbz), kibz2bz(gstore%nkibz)
 integer,intent(out) :: select_qbz_spin(gstore%nqbz, gstore%nsppol)
 integer,intent(out) :: select_kbz_spin(gstore%nkbz, gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: tetra_opt0 = 0
 integer :: nsppol, cnt, spin, ii, all_nproc, my_rank, comm, band ! ierr,
 integer :: ik_bz, ik_ibz, iflag, nk_in_star, gap_err
 logical :: assume_gap
 real(dp) :: ee, abs_erange1, abs_erange2, vmax, cmin
 type(gaps_t) :: gaps
!arrays
 integer,allocatable :: kstar_bz_inds(:)
!----------------------------------------------------------------------

 associate (cryst => gstore%cryst, ebands => gstore%ebands)

 comm = gstore%comm; all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 nsppol = gstore%nsppol

 ! filter k and k + q points according to erange and define gstore%brange_spin automatically.
 ! NB: here we recompute brange_spin
 call wrtout(std_out, sjoin(" Filtering k-points using gstore_erange:", &
                            ltoa(reshape(gstore%erange_spin, [2 * gstore%nsppol]) * Ha_eV), "(eV)"))

 assume_gap = .not. all(gstore%erange_spin < zero)
 gaps = ebands%get_gaps(gap_err)
 if (assume_gap) call gaps%print([std_out]) !, header=msg)

 select_kbz_spin = 0; cnt = 0

 do spin=1,nsppol
   ! Init brange
   gstore%brange_spin(:, spin) = [huge(1), -huge(1)]
   abs_erange1 = abs(gstore%erange_spin(1, spin))
   abs_erange2 = abs(gstore%erange_spin(2, spin))

   if (assume_gap) then
     ! Get CBM and VBM with some tolerance
     vmax = gaps%vb_max(spin) + tol2 * eV_Ha
     cmin = gaps%cb_min(spin) - tol2 * eV_Ha
   else
     vmax = ebands%fermie
     cmin = ebands%fermie
   end if

   do band=1, ebands%mband
     do ik_ibz=1,gstore%nkibz
       !cnt = cnt + 1; if (mod(cnt, all_nproc) /= my_rank) cycle ! MPI parallelism inside comm

       ! Use iflag to filter k-points.
       ee = ebands%eig(band, ik_ibz, spin)
       iflag = 0

       if (abs_erange1 > zero) then
         ! Filter valence states.
         if (ee <= vmax .and. vmax - ee <= abs_erange1) then
           iflag = 1 !; write(std_out, *), "Adding valence band", band, " with ee [eV]: ", ee * Ha_eV
         end if
       end if
       if (abs_erange2 > zero) then
         ! Filter conduction states.
         if (ee >= cmin .and. ee - cmin <= abs_erange2) then
           iflag = 1 !; write(std_out, *)"Adding conduction band", band, " with ee [eV]: ", ee * Ha_eV
         end if
       end if

       if (iflag == 1) then
         gstore%brange_spin(1, spin) = min(gstore%brange_spin(1, spin), band)
         gstore%brange_spin(2, spin) = max(gstore%brange_spin(2, spin), band)

         select case (gstore%kzone)
         case ("ibz")
           ik_bz = kibz2bz(ik_ibz)
           select_kbz_spin(ik_bz, spin) = iflag

         case ("bz")
           call star_from_ibz_idx(ik_ibz, gstore%nkbz, kbz2ibz, nk_in_star, kstar_bz_inds)
           ABI_CHECK(nk_in_star > 0, "Something wrong in star_from_ibz_idx")
           do ii=1,nk_in_star
             ik_bz = kstar_bz_inds(ii)
             select_kbz_spin(ik_bz, spin) = iflag
           end do
           ABI_FREE(kstar_bz_inds)
         end select
       end if

     end do ! band
   end do ! ik_ibz

   !call wrtout(std_out, sjoin("brange_spin:", ltoa(gstore%brange_spin(:, spin))))
   !call wrtout(std_out, sjoin("count_select_kbz:", itoa(count(select_kbz_spin == 1))))

   if (any(gstore%brange_spin(:, spin) == [huge(1), -huge(1)])) then
     ABI_ERROR("Empty list of states inside gstore_erange")
   end if
 end do ! spin

 !call xmpi_sum(select_kbz_spin, comm, ierr)
 call recompute_select_qbz_spin(gstore, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
                                select_kbz_spin, select_qbz_spin)

 call gaps%free()
 end associate

end subroutine gstore_filter_erange__
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_filter_qprange__
!! NAME
!! gstore_filter_qprange__
!!
!! FUNCTION
!! Filter k-points according to qprange
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_filter_qprange__(gstore, dtset, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
                                   select_qbz_spin, select_kbz_spin)

!Arguments ------------------------------------
!scalars
 class(gstore_t),intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 real(dp),intent(in) :: qbz(3, gstore%nqbz)
 real(dp),target,intent(in) :: kbz(3, gstore%nqbz), kibz(3, gstore%nkibz)
 integer,intent(in) :: qbz2ibz(6,gstore%nqbz), qibz2bz(gstore%nqibz)
 integer,intent(in) :: kbz2ibz(6,gstore%nkbz), kibz2bz(gstore%nkibz)
 integer,intent(out) :: select_qbz_spin(gstore%nqbz, gstore%nsppol)
 integer,intent(out) :: select_kbz_spin(gstore%nkbz, gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer :: spin, ik_bz, ik_ibz, gap_err, qprange, ik_calc, nkcalc, mapl_kk(6), my_rank
 type(gaps_t) :: gaps
!arrays
 integer,allocatable :: bstart_ks(:,:), nbcalc_ks(:,:)
 real(dp),allocatable :: kcalc(:,:)
!----------------------------------------------------------------------

 ABI_UNUSED(kbz2ibz)
 ABI_UNUSED(kbz)
 ABI_UNUSED(kibz)
 ABI_UNUSED(qbz2ibz)
 ABI_UNUSED(qbz)
 ABI_UNUSED(qibz2bz)

 associate (cryst => gstore%cryst, ebands => gstore%ebands)

 my_rank = xmpi_comm_rank(gstore%comm)

 qprange = dtset%gw_qprange
 call wrtout(std_out, sjoin(" Filtering k-points using qprange:", itoa(qprange)))
 if (gstore%qzone /= "bz") then
   ABI_ERROR(sjoin('qprange filtering requires gstore_qzone = "bz" while it is: ', gstore%qzone))
 end if

 gaps = ebands%get_gaps(gap_err)
 if (my_rank == 0) call gaps%print([std_out])
 if (gap_err /= 0) then
   ABI_ERROR("Cannot compute fundamental and direct gap (likely metal).")
 end if

 ! Use qp_range to select the interesting k-points and the corresponding bands.
 !
 !    0 --> Compute the QP corrections only for the fundamental and the direct gap.
 ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
 !          bands above and below the Fermi level.
 ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
 !          Include all occupied states and `num` empty states.

 ! Compute nkcalc, kcalc, bstart_ks, nbcalc_ks
 if (qprange /= 0) then
   call sigtk_kcalc_from_qprange(dtset, gstore%cryst, ebands, qprange, nkcalc, kcalc, bstart_ks, nbcalc_ks)
 else
   ! qprange is not specified in the input.
   ! Include direct and fundamental KS gap or include states depending on the position wrt band edges.
   call sigtk_kcalc_from_gaps(dtset, ebands, gaps, nkcalc, kcalc, bstart_ks, nbcalc_ks)
 end if

 ! TODO: kcalc should be spin-dependent to handle magnetic semiconductors.
 select_kbz_spin = 0
 do spin=1,gstore%nsppol
   do ik_calc=1,nkcalc
     if (kpts_map("symrel", ebands%kptopt, gstore%cryst, gstore%krank_ibz, 1, kcalc(:,ik_calc), mapl_kk) /= 0) then
       ABI_ERROR(sjoin("Cannot map kcalc to IBZ with kcalc:", ktoa(kcalc(:,ik_calc))))
     end if
     ! Change select_kbz_spin
     ik_ibz = mapl_kk(1)
     ik_bz = kibz2bz(ik_ibz); select_kbz_spin(ik_bz, spin) = 1
   end do
   ! FIXME: This requires a more careful treatment of (gqk%nb, gqk%nb) matrix that should become (nb1, nb2)
   ! Set brange_spin from bstart_ks and nbcalc_ks. Arrays have shape (nkcalc, nsppol)
   !gstore%brange_spin(:, spin) = [minval(bstart_ks(:,spin)), maxval(bstart_ks(:,spin) + nbcalc_ks(:,spin) - 1)]
 end do

 !call recompute_select_qbz_spin(gstore, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
 !                               select_kbz_spin, select_qbz_spin)

 ABI_FREE(kcalc)
 ABI_FREE(bstart_ks)
 ABI_FREE(nbcalc_ks)

 call gaps%free()
 end associate

end subroutine gstore_filter_qprange__
!!***

!!****f* m_gstore/recompute_select_qbz_spin
!! NAME
!! recompute_select_qbz_spin
!!
!! FUNCTION
!! Recompute select_qbz_spin table after the filtering on the k-points.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine recompute_select_qbz_spin(gstore, qbz, qbz2ibz, qibz2bz, kbz, kibz, kbz2ibz, kibz2bz, &
                                     select_kbz_spin, select_qbz_spin)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
!arrays
 real(dp),intent(in) :: qbz(3, gstore%nqbz)
 integer,intent(in) :: qbz2ibz(6,gstore%nqbz), qibz2bz(gstore%nqibz)
 integer,intent(in) :: kbz2ibz(6,gstore%nkbz), kibz2bz(gstore%nkibz)
 real(dp),target,intent(in) :: kibz(3, gstore%nkibz), kbz(3, gstore%nkbz)
 integer,intent(in) :: select_kbz_spin(gstore%nkbz, gstore%nsppol)
 integer,intent(out) :: select_qbz_spin(gstore%nqbz, gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer :: all_nproc, my_rank, ierr, ii, ik_bz, iq_bz, iq_ibz, ikq_ibz, ikq_bz, len_kpts_ptr, ebands_kptopt, spin
!arrays
 integer,allocatable :: map_kq(:,:)
 real(dp) :: qpt(3)
 real(dp),contiguous, pointer :: kpts_ptr(:,:)
! *************************************************************************

 ABI_UNUSED(kbz2ibz)
 ABI_UNUSED(qbz2ibz)

 all_nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)
 ebands_kptopt = gstore%ebands%kptopt

 select_qbz_spin = 0

 if (gstore%kzone == "ibz") kpts_ptr => kibz
 if (gstore%kzone == "bz")  kpts_ptr => kbz
 len_kpts_ptr = size(kpts_ptr, dim=2)
 ABI_MALLOC(map_kq, (6, len_kpts_ptr))

 select case (gstore%qzone)
 case ("ibz")
   do iq_ibz=1,gstore%nqibz
     if (mod(iq_ibz, all_nproc) /= my_rank) cycle ! MPI parallelism.
     qpt = gstore%qibz(:, iq_ibz)
     iq_bz = qibz2bz(iq_ibz)
     ! k + q_ibz --> k IBZ --> k BZ

     if (kpts_map("symrel", ebands_kptopt, gstore%cryst, gstore%krank_ibz, len_kpts_ptr, kpts_ptr, map_kq, qpt=qpt) /= 0) then
       ABI_ERROR("Cannot map k+q to IBZ!")
     end if

     do ii=1,len_kpts_ptr
       ! get the k-index in BZ
       select case (gstore%kzone)
       case ("bz")
         ik_bz = ii
       case ("ibz")
         ik_bz = kibz2bz(ii)
       end select

       do spin=1,gstore%nsppol
         if (select_kbz_spin(ik_bz, spin) /= 0) then
           ! now, see if q-point connects k-points inside the filtered zone
           ikq_ibz = map_kq(1, ii)
           ikq_bz = kibz2bz(ikq_ibz)
           select_qbz_spin(iq_bz, :) = select_qbz_spin(iq_bz, :) + select_kbz_spin(ikq_bz, :)
         end if
       end do
     end do
   end do ! iq_ibz

 case ("bz")
   do iq_bz=1,gstore%nqbz
     if (mod(iq_bz, all_nproc) /= my_rank) cycle ! MPI parallelism.
     qpt = qbz(:, iq_bz)
     !iq_ibz = qbz2ibz(1, iq_bz)
     ! k + q_bz --> k IBZ --> k BZ

     if (kpts_map("symrel", ebands_kptopt, gstore%cryst, gstore%krank_ibz, len_kpts_ptr, kpts_ptr, map_kq, qpt=qpt) /= 0) then
       ABI_ERROR("Cannot map k+q to IBZ!")
     end if

     ! here we loop over all k-points in iBZ (kzone="ibz") or BZ (kzone="bz")
     do ii=1,len_kpts_ptr
       ! but for each spin, we have to loop only over e-range filtered kpts

       ! get the k-index in BZ
       select case (gstore%kzone)
       case ("bz")
         ik_bz = ii
       case ("ibz")
         ik_bz = kibz2bz(ii)
       end select

       do spin=1,gstore%nsppol
         ! now, see if q-point connects k-points inside the filtered zone
         if (select_kbz_spin(ik_bz, spin) /= 0) then
           ikq_ibz = map_kq(1, ii)
           ikq_bz = kibz2bz(ikq_ibz)
           select_qbz_spin(iq_bz, spin) = select_qbz_spin(iq_bz, spin) + select_kbz_spin(ikq_bz, spin)
         end if
       end do

     end do
   end do
 end select

 call xmpi_sum(select_qbz_spin, gstore%comm, ierr)

 ABI_FREE(map_kq)

end subroutine recompute_select_qbz_spin
!!***

!!****f* m_gstore/gstore_spin2my_is
!! NAME
!! gstore_spin2my_is
!!
!! FUNCTION
!!  Return the local spin index from the global spin index.
!!  0 if this spin is not treated by this MPI proc.
!!
!! SOURCE

integer pure function gstore_spin2my_is(gstore, spin) result(my_is)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 integer,intent(in) :: spin
!----------------------------------------------------------------------

 do my_is=1,gstore%my_nspins
   if (gstore%my_spins(my_is) == spin) return
 end do
 my_is = 0

end function gstore_spin2my_is
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_fill_bks_mask
!! NAME
!! gstore_fill_bks_mask
!!
!! FUNCTION
!!  Fills the bks_mask array defining the set of wavefunctoins that should be read from file
!!  by this MPI rank when computing the KS e-ph matrix elements.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_fill_bks_mask(gstore, mband, nkibz, nsppol, bks_mask)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 integer,intent(in) :: mband, nkibz, nsppol
 logical,intent(out) :: bks_mask(mband, nkibz, nsppol)

!Local variables-------------------------------
!scalars
 integer :: my_is, my_ik, my_iq, spin, ik_ibz, ikq_ibz, ebands_kptopt
 integer :: bstart_k, bstop_k, bstart_kq, bstop_kq
 real(dp) :: weight_q, cpu, wall, gflops
!arrays
 integer,allocatable :: map_kq(:,:)
 real(dp) :: qpt(3)
!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
 associate (cryst => gstore%cryst, ebands => gstore%ebands)

 bks_mask = .False.; ebands_kptopt = gstore%ebands%kptopt

 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   spin = gstore%my_spins(my_is)
   bstart_k = gqk%bstart_k; bstop_k = gqk%bstop_k
   bstart_kq = gqk%bstart_kq; bstop_kq = gqk%bstop_kq

   ! We need the image of this k-point in the IBZ.
   do my_ik=1,gqk%my_nk
     ik_ibz = gqk%my_k2ibz(1, my_ik)
     bks_mask(bstart_k:bstop_k, ik_ibz, spin) = .True.
   end do

   ! As well as the image of k+q in the IBZ.
   ABI_MALLOC(map_kq, (6, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     call gqk%myqpt(my_iq, gstore, weight_q, qpt)

     if (kpts_map("symrel", ebands_kptopt, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, map_kq, qpt=qpt) /= 0) then
       ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qpt)))
     end if

     do my_ik=1,gqk%my_nk
       ikq_ibz = map_kq(1, my_ik)
       bks_mask(bstart_kq:bstop_kq, ikq_ibz, spin) = .True.
     end do
   end do

   ABI_FREE(map_kq)
   end associate
 end do ! my_is

 call cwtime_report(" gstore_fill_bks_mask", cpu, wall, gflops)
 end associate

end subroutine gstore_fill_bks_mask
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_fill_bks_mask_with_pp
!! NAME
!! gstore_fill_bks_mask_with_pp
!!
!! FUNCTION
!!  Fill the bks_mask array defining the set of states that should be read from the WFK file
!!  by this MPI rank when computing the GWPT e-ph matrix elements in which we have
!!  to consider k+q, k-q and k+q-p as well as the sum over states (bsum)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_fill_bks_mask_pp_mesh(gstore, ecut, mband, nkibz, nsppol, my_pp_start_spin, my_pp_stop_spin, pp_mesh, &
                                        my_bsum_start, my_bsum_stop, bks_mask, mpw, gmax)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 real(dp),intent(in) :: ecut
 integer,intent(in) :: mband, nkibz, nsppol
 integer :: my_pp_start_spin(nsppol), my_pp_stop_spin(nsppol)
 type(kmesh_t),intent(in) :: pp_mesh
 integer,intent(in) :: my_bsum_start(nsppol), my_bsum_stop(nsppol)
 logical,intent(out) :: bks_mask(mband, nkibz, nsppol)
 integer,intent(out) :: mpw, gmax(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: istwfk1 = 1
 integer :: b1, b2, ebands_kptopt, ierr, ikq_ibz, ik_ibz, ipp_bz, my_ik, my_iq, spin, my_is, onpw, my_mpw, my_gmax(3)
 real(dp) :: weight_q, cpu, wall, gflops
 type(gqk_t),pointer :: gqk
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: map_kq(:,:), gtmp(:,:)
 real(dp) :: qpt(3), pp(3), kk(3)
!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 ! TODO: These loops can be parallelized using bsum_comm and pert_comm if needed.
 bks_mask = .False.; cryst => gstore%cryst
 ebands_kptopt = gstore%ebands%kptopt
 mpw = 0; gmax = 0

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is); spin = gstore%my_spins(my_is)
   ! These are the first and last band indices used in the sum over states (possibly MPI-distributed)
   b1 = my_bsum_start(spin); b2 = my_bsum_stop(spin)

   ABI_MALLOC(map_kq, (6, gqk%my_nk))

   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik)

     ! We need the image of this k-point in the IBZ for the incoming state |psi_nk>.
     bks_mask(gqk%bstart_k:gqk%bstop_k, ik_ibz, spin) = .True.  ! gqk%n_start, gqk%n_stop

     ! We also need the image of k+q in the IBZ for the outgoing state <psi_mkq|.
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_q, qpt)
       if (kpts_map("symrel", ebands_kptopt, cryst, gstore%krank_ibz, 1, kk, map_kq, qpt=qpt) /= 0) then
         ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qpt)))
       end if
       ikq_ibz = map_kq(1, 1)
       bks_mask(gqk%bstart_kq:gqk%bstop_kq, ikq_ibz, spin) = .True. ! gqk%m_start, gqk%m_stop
     end do ! my_iq

   end do ! my_ok

   ! We also need the image of k-p in the IBZ for the pp wavevectors treated by this MPI rank.
   ! These states are summed over so use b1 and b2.
   do ipp_bz=my_pp_start_spin(spin), my_pp_stop_spin(spin)
     pp = pp_mesh%bz(:,ipp_bz)
     if (kpts_map("symrel", ebands_kptopt, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, map_kq, qpt=-pp) /= 0) then
       ABI_ERROR(sjoin("Cannot map k-p to IBZ with qpt:", ktoa(qpt), "and pp:", ktoa(pp)))
     end if

     do my_ik=1,gqk%my_nk
       ikq_ibz = map_kq(1, my_ik)
       bks_mask(b1:b2, ikq_ibz, spin) = .True.
       ! Compute g-sphere, returns onpw. Note istwfk == 1.
       kk = gqk%my_kpts(:, my_ik)
       call get_kg(kk-pp, istwfk1, ecut, gstore%cryst%gmet, onpw, gtmp, mpw=mpw, gmax=gmax)
       ABI_FREE(gtmp)
     end do
   end do ! ipp_bz

   ! We also need the image of k+q-p in the IBZ for the pp wavevectors treated by this MPI rank.
   ! These states are summed over so use b1 and b2.
   do my_iq=1,gqk%my_nq
     call gqk%myqpt(my_iq, gstore, weight_q, qpt)

     do ipp_bz=my_pp_start_spin(spin), my_pp_stop_spin(spin)
       pp = pp_mesh%bz(:,ipp_bz)

       if (kpts_map("symrel", ebands_kptopt, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, map_kq, qpt=qpt-pp) /= 0) then
         ABI_ERROR(sjoin("Cannot map k+q-p to IBZ with qpt:", ktoa(qpt), "and pp:", ktoa(pp)))
       end if

       do my_ik=1,gqk%my_nk
         ikq_ibz = map_kq(1, my_ik)
         bks_mask(b1:b2, ikq_ibz, spin) = .True.
         ! Compute g-sphere, returns onpw. Note istwfk = 1.
         kk = gqk%my_kpts(:, my_ik)
         call get_kg(kk+qpt-pp, istwfk1, ecut, gstore%cryst%gmet, onpw, gtmp, mpw=mpw, gmax=gmax)
         ABI_FREE(gtmp)
       end do
     end do ! ipp_bz

   end do ! my_iq

   ABI_FREE(map_kq)
 end do ! my_is

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, gstore%comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, gstore%comm, ierr)

 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))
 call cwtime_report(" gstore_fill_bks_mask_pp_mesh", cpu, wall, gflops)

end subroutine gstore_fill_bks_mask_pp_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_get_mpw_gmax
!! NAME
!! gstore_get_mpw_gmax
!!
!! FUNCTION
!! Compute the maximum number of PWs for all possible k+q treated.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_get_mpw_gmax(gstore, ecut, mpw, gmax)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 real(dp),intent(in) :: ecut
 integer,intent(out) :: mpw, gmax(3)

!Local variables-------------------------------
 integer,parameter :: istwfk1 = 1
 integer :: my_is, my_ik, my_iq, spin, onpw, ierr, my_mpw, ipx, ipy, ipz, pp_max
 real(dp) :: weight_q, cpu, wall, gflops
!arrays
 integer :: my_gmax(3)
 integer,allocatable :: gtmp(:,:)
 real(dp) :: kk(3), qpt(3), pp(3)
!----------------------------------------------------------------------

 mpw = 0; gmax = 0

 ! TODO: This is an hotspot due to the double loop over k and q. Should use a geometrical approach to compute mpw and gmax.
 call wrtout(std_out, " Computing mpw. This may take some time for dense k/q meshes...", pre_newlines=1)
 call cwtime(cpu, wall, gflops, "start")

 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   spin = gstore%my_spins(my_is)
   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik)

     ! Compute g-sphere, returns onpw. Note istwfk == 1.
     call get_kg(kk, istwfk1, ecut, gstore%cryst%gmet, onpw, gtmp, mpw=mpw, gmax=gmax)
     ABI_FREE(gtmp)

     pp_max = 0
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_q, qpt)
       ! TODO: g0 umklapp here can enter into play! gmax could not be large enough!
       do ipz=-pp_max,pp_max
          do ipy=-pp_max,pp_max
            do ipx=-pp_max,pp_max
             pp = [ipx, ipy, ipz] * half
             call get_kg(kk + qpt - pp, 1, ecut, gstore%cryst%gmet, onpw, gtmp, mpw=mpw, gmax=gmax)
             ABI_FREE(gtmp)
           end do
         end do
       end do
     end do ! my_iq

   end do ! my_ik
   end associate
 end do ! my_is

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, gstore%comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, gstore%comm, ierr)

 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))
 call cwtime_report(" gstore_get_mpw_gmax", cpu, wall, gflops)

end subroutine gstore_get_mpw_gmax
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_get_lambda_iso_iw
!! NAME
!! gstore_get_lambda_iso_iw
!!
!! FUNCTION
!!  Compute isotropic lambda along the imaginary axis
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_get_lambda_iso_iw(gstore, nw, imag_w, lambda)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 integer,intent(in) :: nw
 real(dp),intent(in) :: imag_w(nw)
 real(dp),intent(out) :: lambda(nw)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, nb_k, nb_kq
 real(dp) :: g2, wqnu, weight_k, weight_q
!arrays
 real(dp) :: qpt(3)
 real(dp),allocatable :: dbl_delta_q(:,:,:), g2_pmnk(:,:,:,:)
!----------------------------------------------------------------------

 ABI_CHECK(gstore%qzone == "bz", "gstore_get_lambda_iso_iw assumes qzone == `bz`")
 !if (gstore%check_cplex_qkzone_gmode(cplex1, "bz", kzone, gmode, kfilter) result(ierr)

 lambda = zero
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq
   ABI_CHECK_IEQ(nb_k, nb_kq, "gqk_dbldelta_qpt does not support nb_k != nb_kq")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbl_delta_q, (nb_kq, nb_k, gqk%my_nk))
   ABI_MALLOC(g2_pmnk, (gqk%my_npert, nb_kq, nb_k, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     ! Compute integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, gstore%dtset%eph_intmeth, gstore%dtset%eph_fsmear, qpt, weight_q, dbl_delta_q)

     ! Copy data to improve memory access in the loops below.
     g2_pmnk = gqk%my_g2(:,:,my_iq,:,:)

     do my_ik=1,gqk%my_nk
       weight_k = gqk%my_wtk(my_ik)
       do in_k=1,nb_k
         do im_kq=1,nb_kq
           do my_ip=1,gqk%my_npert
             g2 = g2_pmnk(my_ip, im_kq, in_k, my_ik)
             ! TODO: handle wqnu ~ 0
             wqnu = gqk%my_wnuq(my_ip, my_iq)
             lambda(:) = lambda(:) + &
               two * wqnu / (imag_w(:) ** 2 + wqnu ** 2) * g2 * weight_k * weight_q * dbl_delta_q(im_kq, in_k, my_ik)
           end do
         end do
       end do
     end do
   end do ! my_iq

   ABI_FREE(dbl_delta_q)
   ABI_FREE(g2_pmnk)
   end associate
 end do ! my_is

 ! Take into account collinear spin
 lambda = lambda * (two / (gstore%nsppol * gstore%dtset%nspinor))
 call xmpi_sum(lambda, gstore%comm, ierr)

end subroutine gstore_get_lambda_iso_iw
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_get_a2fw
!! NAME
!! gstore_get_a2fw
!!
!! FUNCTION
!!  Compute Eliashberg function a^2F(omega).
!!
!! INPUTS
!!  nw: Number of frequencies.
!!  wmesh: Frequency mesh.
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_get_a2fw(gstore, dtset, nw, wmesh, a2fw)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: nw
 real(dp),intent(in) :: wmesh(nw)
 real(dp),intent(out) :: a2fw(nw)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, timrev_q, ii, ik_ibz, nb_k, nb_kq
 real(dp) :: g2, wqnu, weight_k, weight_q, cpu, wall, gflops
 type(lgroup_t) :: lg_myq
 character(len=500) :: msg, kk_string !, qq_bz_string,
!arrays
 real(dp) :: qpt(3), kk(3)
 real(dp),allocatable :: dbl_delta_q(:,:,:), g2_mnkp(:,:,:,:), deltaw_nuq(:)
!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
 call wrtout(std_out, sjoin(" Computing a^2F(w) with ph_smear:", ftoa(gstore%dtset%ph_smear * Ha_meV), "(meV)"), pre_newlines=1)
 ABI_CHECK(gstore%qzone == "bz", "gstore_get_lambda_iso_iw assumes qzone == `bz`")

 ! Check consistency of little group options.
 ABI_CHECK(gstore%check_little_group(dtset, msg) == 0, msg)
 ABI_MALLOC(deltaw_nuq, (nw))

 a2fw = zero
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq
   ABI_CHECK_IEQ(nb_k, nb_kq, "gqk_dbldelta_qpt does not support nb_k != nb_kq")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbl_delta_q, (nb_kq, nb_k, gqk%my_nk))
   ABI_MALLOC(g2_mnkp, (nb_kq, nb_k, gqk%my_nk, gqk%my_npert))

   do my_iq=1,gqk%my_nq
     ! Compute integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, gstore%dtset%eph_intmeth, gstore%dtset%eph_fsmear, qpt, weight_q, dbl_delta_q)

     ! Copy data to improve memory access in the loops below.
     do my_ip=1,gqk%my_npert
       g2_mnkp(:,:,:,my_ip) = gqk%my_g2(my_ip,:,my_iq,:,:)
     end do

     ! Compute the little group of the q-point so that we only need to sum g(k,q) for k in the IBZ_q
     if (dtset%gstore_use_lgq /= 0) then
       timrev_q = kpts_timrev_from_kptopt(gstore%qptopt)
       call lg_myq%init(cryst, qpt, timrev_q, gstore%nkbz, gstore%kbz, gstore%nkibz, gstore%kibz, xmpi_comm_self)
     end if

     do my_ip=1,gqk%my_npert
       wqnu = gqk%my_wnuq(my_ip, my_iq)
       deltaw_nuq = gaussian(wmesh - wqnu, gstore%dtset%ph_smear)
       do my_ik=1,gqk%my_nk
         kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik); weight_k = gqk%my_wtk(my_ik)

         ! Handle little group and weight
         if (dtset%gstore_use_lgq /= 0) then
           ii = lg_myq%findq_ibzk(kk)
           if (ii == -1) then
             kk_string = ktoa(kk)
             call wrtout(std_out, sjoin(" my_ik:", itoa(my_ik), kk_string, " not in IBZ_q --> skipping iteration"))
             cycle
             weight_k = lg_myq%weights(ii)
             ! TODO: Check fillvalue (should be zero)
           end if
         end if

         do in_k=1,nb_k
           do im_kq=1,nb_kq
             g2 = g2_mnkp(im_kq, in_k, my_ik, my_ip)
             a2fw(:) = a2fw(:) + deltaw_nuq(:) * g2 * weight_k * weight_q * dbl_delta_q(im_kq, in_k, my_ik)
           end do
         end do
       end do
     end do

     call lg_myq%free()
   end do ! my_iq

   ABI_FREE(dbl_delta_q)
   ABI_FREE(g2_mnkp)
   end associate
 end do ! my_is

 ABI_FREE(deltaw_nuq)

 ! Take into account collinear spin and N(eF) TODO
 a2fw = a2fw * (two / (gstore%nsppol * gstore%dtset%nspinor))
 call xmpi_sum(a2fw, gstore%comm, ierr)

 call cwtime_report(" gstore_get_a2fw", cpu, wall, gflops)

end subroutine gstore_get_a2fw
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_free
!! NAME
!! gstore_free
!!
!! FUNCTION
!!  Free dynamic memory in gstore_t.
!!
!! SOURCE

subroutine gstore_free(gstore)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore

!Local variables-------------------------------
 integer :: my_is
!----------------------------------------------------------------------

 do my_is=1,gstore%my_nspins
   call gstore%gqk(my_is)%free()
 end do
 ABI_SFREE(gstore%gqk)

 ABI_SFREE(gstore%delta_ef_kibz_spin)
 ABI_SFREE(gstore%qibz)
 ABI_SFREE(gstore%wtq)
 ABI_SFREE(gstore%my_spins)
 ABI_SFREE(gstore%brange_spin)
 ABI_SFREE(gstore%glob_nk_spin)
 ABI_SFREE(gstore%glob_nq_spin)
 ABI_SFREE(gstore%kglob2bz)
 ABI_SFREE(gstore%kbz2ibz)
 ABI_SFREE(gstore%erange_spin)
 ABI_SFREE(gstore%qbz)
 ABI_SFREE(gstore%kbz)

 if (gstore%ebands_owns_memory) call gstore%ebands%free()

 call gstore%krank_ibz%free()
 call gstore%qrank_ibz%free()

end subroutine gstore_free
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_myqpt
!! NAME
!! gqk_myqpt
!!
!! FUNCTION
!!  Return the weight and the reduced coordinates of the q-point from the local index my_iq.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure subroutine gqk_myqpt(gqk, my_iq, gstore, weight, qpt)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 class(gstore_t),intent(in) :: gstore
 integer,intent(in) :: my_iq
 real(dp),intent(out) :: weight, qpt(3)

!Local variables ------------------------------
 integer :: iq_ibz, isym_q, trev_q, tsign_q, g0_q(3)
 logical :: isirr_q
!----------------------------------------------------------------------

 iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
 trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
 isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
 tsign_q = 1; if (trev_q == 1) tsign_q = -1

 ! NB: Use symrec convention for q
 qpt = tsign_q * matmul(gstore%cryst%symrec(:,:,isym_q), gstore%qibz(:, iq_ibz)) + g0_q

 select case(gstore%qzone)
 case ("ibz")
   weight = gstore%wtq(iq_ibz)
 case ("bz")
   weight = one / gstore%nqbz
 end select

end subroutine gqk_myqpt
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_dbldelta_qpt
!! NAME
!! gqk_dbldelta_qpt
!!
!! FUNCTION
!!  Note that k/q weights are not included in dbl_delta_q
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gqk_dbldelta_qpt(gqk, my_iq, gstore, eph_intmeth, eph_fsmear, qpt, weight_q, dbl_delta_q)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: my_iq, eph_intmeth
 real(dp),intent(in) :: eph_fsmear
 real(dp),intent(out) :: qpt(3), weight_q, dbl_delta_q(gqk%nb, gqk%nb, gqk%my_nk)

!Local variables ------------------------------
!scalars
 real(dp), parameter :: min_smear = tol9
 integer :: nb_k, nb_kq, nkbz, spin, my_ik, ib, ib1, ib2, band1, band2, nesting
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 integer :: ikq_ibz, isym_kq, trev_kq, tsign_kq, g0_kq(3), ii, i1, i2, i3, cnt, ik_bz, ltetra
 real(dp) :: g1, g2, sigma !, weight_k !, cpu, wall, gflops
 logical :: isirr_k, isirr_kq, use_adaptive
 type(ebands_t), pointer :: ebands
 type(crystal_t), pointer :: cryst
 type(krank_t) :: my_krank
!arrays
 integer :: nge(3), ngw(3)
 integer,allocatable :: my_kqmap(:,:), kmesh_map(:,:)
 real(dp) :: kk(3), kmesh_cartvec(3,3), rlatt(3,3), klatt(3,3), vb_k(3, gqk%nb), vb_kq(3, gqk%nb)
 real(dp),allocatable :: eig_k(:,:), eig_kq(:,:), kmesh(:,:), wght_bz(:,:,:)
!----------------------------------------------------------------------

 nb_k = gqk%nb_k; nb_kq = gqk%nb_kq; nkbz = gstore%nkbz; spin = gqk%spin
 ebands => gstore%ebands; cryst => gstore%cryst

 ABI_CHECK_IEQ(nb_k, nb_kq, "gqk_dbldelta_qpt does not support nb_k != nb_kq")

 call gqk%myqpt(my_iq, gstore, weight_q, qpt)

 ! The double delta with tetra is ill-defined for q == 0. In this case we fall back to gaussian.
 nesting = merge(1, 0, abs(eph_intmeth) == 2 .and. all(abs(qpt) < tol12))

 rlatt = gstore%ebands%kptrlatt; call matr3inv(rlatt, klatt)
 kmesh_cartvec(:, 1) = cryst%gprimd(:,1)*klatt(1,1) + cryst%gprimd(:,2)*klatt(2,1) + cryst%gprimd(:,3)*klatt(3,1)
 kmesh_cartvec(:, 2) = cryst%gprimd(:,1)*klatt(1,2) + cryst%gprimd(:,2)*klatt(2,2) + cryst%gprimd(:,3)*klatt(3,2)
 kmesh_cartvec(:, 3) = cryst%gprimd(:,1)*klatt(1,3) + cryst%gprimd(:,2)*klatt(2,3) + cryst%gprimd(:,3)*klatt(3,3)
 ! TODO: It seems that two_pi is not needed here!

 if (abs(eph_intmeth) == 1 .or. nesting /= 0) then
   use_adaptive = eph_fsmear < zero .or. abs(eph_intmeth) == 2
   if (use_adaptive) then
     ABI_CHECK(allocated(gqk%vk_cart_ibz), "vk_cart_ibz should be allocated when use_adaptive is .True.")
   end if

   ! Find k + q in the IBZ for all my k-points.
   ABI_MALLOC(my_kqmap, (6, gqk%my_nk))
   if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, my_kqmap, qpt=qpt) /= 0) then
     ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qpt)))
   end if

   ! Init default sigma
   sigma = eph_fsmear

   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik)

     ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
     trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5, my_ik)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     tsign_k = 1; if (trev_k == 1) tsign_k = -1

     ikq_ibz = my_kqmap(1, my_ik); isym_kq = my_kqmap(2, my_ik)
     trev_kq = my_kqmap(6, my_ik); g0_kq = my_kqmap(3:5, my_ik)
     isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
     tsign_kq = 1; if (trev_kq == 1) tsign_kq = -1

     if (use_adaptive) then
       ! If k or k+q is not in the IBZ, we need to recostruct the value by symmetry using v(Sq) = S v(q).
       ! Use transpose(R) because we are using the tables for the wavefunctions
       ! In this case listkk has been called with symrec and use_symrec=False
       ! so q_bz = S^T q_ibz where S is the isym_kq symmetry
       vb_k = gqk%vk_cart_ibz(:,:,ik_ibz)
       vb_kq = gqk%vk_cart_ibz(:,:,ikq_ibz)

       if (.not. isirr_k) then
         do ib=1,nb_k
           vb_k(:,ib) = tsign_k * matmul(transpose(cryst%symrel_cart(:,:,isym_k)), vb_k(:,ib))
         end do
       end if
       if (.not. isirr_kq) then
         do ib=1,nb_kq
           vb_kq(:,ib) = tsign_kq * matmul(transpose(cryst%symrel_cart(:,:,isym_kq)), vb_kq(:,ib))
         end do
       end if
     end if

     do ib2=1,nb_k
       band2 = ib2 + gqk%bstart - 1
       if (use_adaptive) then
         sigma = max(maxval([(abs(dot_product(vb_k(:, ib2), kmesh_cartvec(:,ii))), ii=1,3)]), min_smear)
         !write(std_out, *)"sigma:", sigma * Ha_eV
       end if
       g2 = gaussian(ebands%eig(band2, ik_ibz, spin) - ebands%fermie, sigma)

       do ib1=1,nb_kq
         band1 = ib1 + gqk%bstart - 1
         if (use_adaptive) then
           sigma = max(maxval([(abs(dot_product(vb_kq(:, ib1), kmesh_cartvec(:,ii))), ii=1,3)]), min_smear)
         end if
         g1 = gaussian(ebands%eig(band1, ikq_ibz, spin) - ebands%fermie, sigma)
         dbl_delta_q(ib1, ib2, my_ik) = g1 * g2 ! / fs%nktot
       end do

     end do
   end do

   ABI_FREE(my_kqmap)

 else if (abs(eph_intmeth) == 2) then

   ABI_CHECK(isdiagmat(ebands%kptrlatt), "kptrlatt must be diagonal when tetra is used.")
   ABI_CHECK(ebands%nshiftk == 1, "nshiftk must be 1 when tetra is used")
   nge = get_diag(ebands%kptrlatt); ngw = nge
   ABI_CHECK_IEQ(nkbz, product(nge(1:3)), "Wrong nge")

   ! Compute eig_k and eig_kq in full BZ for the relevant bands around Ef.
   ABI_MALLOC(kmesh, (3, nkbz))
   ABI_MALLOC(eig_k, (nb_k, nkbz))
   ABI_MALLOC(eig_kq, (nb_kq, nkbz))

   ! Technical problems:
   !
   ! 1) libtetrabz works with the BZ and assumes a certaing ordering of the k-points (see below)
   !    so we have to fill the array with eig_k and eig_kq from the IBZ by remapping the libtetra kk
   !    to the Abinit IBZ

   ! 2) The dbldelta weights are given in the BZ, while the caller requires weights for k in the IBZ
   !    and moreover only for the IBZ k-point treated by this MPI proc.

   ik_bz = 0
   do i3=0,nge(3) - 1
     do i2=0,nge(2) - 1
       do i1=0,nge(1) - 1
         ik_bz = ik_bz + 1
         kk = ([i1, i2, i3] + ebands%shiftk(:, 1)) / nge(:)
         kmesh(:, ik_bz) = kk
       end do
     end do
   end do

   ! Map libtetra BZ mesh to IBZ and fill eig_k
   !call cwtime(cpu, wall, gflops, "start")
   ABI_MALLOC(kmesh_map, (6, nkbz))

   ! Find correspondence between libtetra mesh and the IBZ.
   if (kpts_map("symrec", ebands%kptopt, cryst, gstore%krank_ibz, nkbz, kmesh, kmesh_map) /= 0) then
     ABI_ERROR("Cannot map libtetra mesh to IBZ")
   end if

   do ik_bz=1,nkbz
     ik_ibz = kmesh_map(1, ik_bz)
     eig_k(:, ik_bz) = ebands%eig(gqk%bstart_k:gqk%bstop_k, ik_ibz, spin) - ebands%fermie
   end do

   ! Map libtetra BZ mesh + q to IBZ and fill eig_kq
   if (kpts_map("symrec", ebands%kptopt, cryst, gstore%krank_ibz, nkbz, kmesh, kmesh_map, qpt=qpt) /= 0) then
     ABI_ERROR(sjoin("Cannot map libtetra k+q to IBZ with qpt:", ktoa(qpt)))
   end if

   do ik_bz=1,nkbz
     ikq_ibz = kmesh_map(1, ik_bz)
     eig_kq(:, ik_bz) = ebands%eig(gqk%bstart_kq:gqk%bstop_kq, ikq_ibz, spin) - ebands%fermie
   end do

   ABI_FREE(kmesh_map)
   !call cwtime_report(" kmesh_map", cpu, wall, gflops)

   ! Call libtetra routine to compute weights for double delta integration.
   ! Note that libtetra assumes Ef set to zero.
   ! TODO: Average weights over degenerate states?
   ! NB: This is a botleneck, can pass comm_kp

   ! Select option for double delta with tetra.
   !  2 for the optimized tetrahedron method.
   ! -2 for the linear tetrahedron method.
   ltetra = 0
   if (eph_intmeth ==  2) ltetra = 2
   if (eph_intmeth == -2) ltetra = 1

   ABI_MALLOC(wght_bz, (nb_k, nb_kq, nkbz))
   call libtetrabz_dbldelta(ltetra, gstore%cryst%gprimd, nb_k, nge, eig_k, eig_kq, ngw, wght_bz) !, comm=comm)
   !call cwtime_report(" libtetrabz_dbldelta", cpu, wall, gflops)

   call my_krank%init(gqk%my_nk, gqk%my_kpts)

   ! Reindex from full BZ to my set of kpoints and rescale weights.
   cnt = 0
   do ik_bz=1,nkbz
     my_ik = my_krank%get_index(kmesh(:, ik_bz))
     if (my_ik /= -1) then
       dbl_delta_q(:,:,my_ik) = wght_bz(:,:,ik_bz) * gstore%nkbz
       cnt = cnt + 1
     end if
   end do

   ! FIXME: bug if k-point (and q-point) parallelism.
   ABI_CHECK_IEQ(cnt, gqk%my_nk, sjoin("cnt != my_nk, ", itoa(cnt), itoa(gqk%my_nk)))
   call my_krank%free()
   !call cwtime_report(" transfer", cpu, wall, gflops)

   ABI_FREE(wght_bz)
   ABI_FREE(kmesh)
   ABI_FREE(eig_k)
   ABI_FREE(eig_kq)

 else
   ABI_ERROR(sjoin("Invalid eph_intmeth:", itoa(eph_intmeth)))
 end if

end subroutine gqk_dbldelta_qpt
!!***
!----------------------------------------------------------------------

!!****f* m_gstore/gqk_free
!! NAME
!! gqk_free
!!
!! FUNCTION
!!  Free dynamic memory in gqk_t instance.
!!
!! SOURCE

subroutine gqk_free(gqk)

!Arguments ------------------------------------
 class(gqk_t),intent(inout) :: gqk
!----------------------------------------------------------------------

 ABI_SFREE(gqk%my_k2ibz)
 ABI_SFREE(gqk%my_kpts)
 ABI_SFREE(gqk%my_wtk)
 ABI_SFREE(gqk%my_q2ibz)
 ABI_SFREE(gqk%my_q2bz)
 ABI_SFREE(gqk%my_k2glob)
 ABI_SFREE(gqk%my_q2glob)
 ABI_SFREE(gqk%my_wnuq)
 ABI_SFREE(gqk%my_displ_cart)
 ABI_SFREE(gqk%my_g)
 ABI_SFREE(gqk%my_gq0nm_atm)
 ABI_SFREE(gqk%my_g2)
 ABI_SFREE(gqk%my_pertcases)
 ABI_SFREE(gqk%vk_cart_ibz)
 ABI_SFREE(gqk%vkmat_cart_ibz)

 call gqk%wan%free()

 ! Free MPI communicators
 call gqk%kpt_comm%free(); call gqk%qpt_comm%free(); call gqk%qpt_kpt_comm%free()
 call gqk%pert_comm%free(); call gqk%band_comm%free(); call gqk%qpt_pert_comm%free()
 call gqk%pert_ppsum_comm%free(); call gqk%pert_ppsum_bsum_comm%free()
 call gqk%comm%free(); call gqk%bsum_comm%free(); call gqk%pp_sum_comm%free()

end subroutine gqk_free
!!***

!!****f* m_gstore/gstore_get_missing_qbz_spin
!! NAME
!! gstore_get_missing_qbz_spin
!!
!! FUNCTION
!!  Return the number of (q-points, spin) entries that have been computed.
!!
!! SOURCE

subroutine gstore_get_missing_qbz_spin(gstore, done_qbz_spin, ndone, nmiss)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 integer,intent(in) :: done_qbz_spin(gstore%nqbz, gstore%nsppol)
 integer,intent(out) :: ndone, nmiss

!Local variables ------------------------------
 integer :: my_is, my_iq, iq_bz, spin, ierr, nscale
!----------------------------------------------------------------------

 nmiss = 0
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   spin = gstore%my_spins(my_is)
   do my_iq=1,gqk%my_nq
     iq_bz = gqk%my_q2bz(my_iq)
     if (done_qbz_spin(iq_bz, spin) == 0) nmiss = nmiss + 1
     if (done_qbz_spin(iq_bz, spin) == 1) ndone = ndone + 1
   end do ! my_iq
   ! Rescale to avoid overcounting in xmpi_summ
   nscale = (gqk%pert_comm%nproc * gqk%bsum_comm%nproc * gqk%pp_sum_comm%nproc)
   nmiss = nmiss / nscale
   ndone = ndone / nscale
   end associate
 end do ! my_is

 call xmpi_sum(ndone, gstore%comm, ierr)
 call xmpi_sum(nmiss, gstore%comm, ierr)

end subroutine gstore_get_missing_qbz_spin
!!***

!!****f* m_gstore/gstore_set_perts_distrib
!! NAME
!! gstore_set_perts_distrib
!!
!! FUNCTION
!! Activate parallelism over perturbations at the level of the DVDB file.
!!
!! SOURCE

subroutine gstore_set_perts_distrib(gstore, cryst, dvdb, my_npert)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 type(crystal_t),intent(in) :: cryst
 type(dvdb_t),intent(inout) :: dvdb
 integer,intent(out) :: my_npert

!Local variables ------------------------------
!scalars
 integer :: my_is, spin
!arrays
 integer,allocatable :: my_pinfo(:,:), pert_table(:,:)
!----------------------------------------------------------------------

 my_npert = cryst%natom * 3
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   spin = gstore%my_spins(my_is)
   if (gqk%pert_comm%nproc > 1) then
     ! Activate parallelism over perturbations
     ! Build table with list of perturbations treated by this MPI rank inside pert_comm.
     !ABI_WARNING("GSTORE with pert_comm%nproc > 1 not tested")
     my_npert = gqk%my_npert
     call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, gqk%pert_comm%value)
     call dvdb%set_pert_distrib(my_npert, cryst%natom * 3, my_pinfo, pert_table, gqk%pert_comm%value)
     ABI_CHECK(all(my_pinfo(3, :) == gqk%my_pertcases), "my_pinfo(3, :) != gqk%my_pertcases")

     ABI_FREE(my_pinfo)
     ABI_FREE(pert_table)
   end if
   end associate
 end do

end subroutine gstore_set_perts_distrib
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_compute
!! NAME
!!  gstore_compute
!!
!! FUNCTION
!!  Compute MPI-distributed e-ph matrix elements
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
!!  GSTORE.nc file
!!
!! SOURCE

subroutine gstore_compute(gstore, wfk0_path, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, &
                          pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
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
 integer,parameter :: tim_getgh1c = 1, berryopt0 = 0, ider0 = 0, idir0 = 0, LOG_MODQ = 5, master = 0, ndat1 = 1
 integer :: my_rank,nproc,nproc_lim,mband,nsppol,nkibz,idir,ipert, iq_bz
 integer :: cplex,natom,natom3,ipc,nspinor, nskip_tetra_kq, timrev_k, timrev_q
 integer :: band_k, in_k, im_kq, ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq, nb_k, nb_kq
 integer :: my_ik, my_is, comm_rpt, my_npert, my_ip, my_iq, spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: mpw, nb,ierr,cnt, n1,n2,n3,n4,n5,n6,nspden,ndone, db_iqpt
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf, nkpg_k, nkpg_kq, qbuf_size, iqbuf_cnt, root_ncid, spin_ncid, ncerr
 integer :: ii, my_nqibz, iq_start, iq_ibz, isym_q, trev_q, prev_iqbz
 real(dp) :: cpu, wall, gflops, cpu_q, wall_q, gflops_q, cpu_all, wall_all, gflops_all
 real(dp) :: ecut, eshift, eig0nk, weight_q, weight_k
 logical :: gen_eigenpb, isirr_k, isirr_kq, isirr_q, print_time, need_ftinterp, qq_is_gamma
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_ham_kq
 type(rf_hamiltonian_type) :: rf_ham_kq
 type(ddkop_t) :: ddkop
 type(gqk_t),pointer :: gqk
 type(lgroup_t) :: lg_myq
 character(len=5000) :: msg, qq_bz_string, kk_string
!arrays
 integer :: g0_k(3), g0_kq(3), g0_q(3), work_ngfft(18),gmax(3),indkk_kq(6,1), units(2), qbz2dvdb(6)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), nband(:,:), wfd_istwfk(:), qmap_symrec(:,:)
 integer,allocatable :: iq_buf(:,:), done_qbz_spin(:,:), my_iqibz_inds(:)
 !integer,allocatable :: qibz2dvdb(:) !, displs(:), recvcounts(:)
 real(dp) :: kk_bz(3),kq_bz(3),kk_ibz(3),kq_ibz(3), qq_bz(3), qq_ibz(3), vk(3), phfr_qq(3*cryst%natom)
 real(dp),allocatable :: displ_cart_qibz(:,:,:,:), displ_red_qibz(:,:,:,:), pheigvec_qibz(:,:,:,:)
 real(dp),allocatable :: displ_cart_qbz(:,:,:,:), displ_red_qbz(:,:,:,:), pheigvec_qbz(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:), kinpw_k(:), kinpw_kq(:), kpg_kq(:,:), kpg_k(:,:)
 real(dp),allocatable :: ffnl_k(:,:,:,:), ffnl_kq(:,:,:,:), ph3d_k(:,:,:), ph3d_kq(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:), gkq_atm(:,:,:,:), gkq_nu(:,:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:), kets_k(:,:,:), h1_kets_kq(:,:,:), cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:), vlocal(:,:,:,:), vlocal1(:,:,:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:), gvnlx1(:,:), work(:,:,:,:)
 real(dp),allocatable :: gs1c_kq(:,:), vk_cart_ibz(:,:,:) !, vkmat_cart_ibz(:,:,:,:)
 real(dp),allocatable :: my_gbuf(:,:,:,:,:,:), buf_wqnu(:,:), buf_eigvec_cart(:,:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)
 type(lgroup_t),allocatable :: lg_myk(:)
 !type(lgroup_t) :: lg_myk(:)
!************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 units = [std_out, ab_out]

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 ! Check MPI distribution: abort if there are any idle processes
 ! NOTE that we have to perform this check only at the level of gstore%compute:
 ! if gstore object is loaded from a *GSTORE.nc file, there are no calls to individual
 ! phonons/electrons distribution, hence k/q-iBZ is not a limiting factor anymore
 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   nproc = gqk%comm%nproc
   nproc_lim = min(gstore%nkibz, gstore%nqibz)
   if (nproc > nproc_lim) then
      write(msg, "(a,i0,a,i0,4a)") "gstore%compute: nproc=", nproc, " > min(nkibz, nqibz)=", &
        nproc_lim, ch10, "This will lead to idle processes, which are not supported.", ch10, &
        "Please decrease the total number of CPUs as the size of your problem is relatively small"
      ABI_ERROR(msg)
   endif
 enddo

 ! This parameter defines the size of the q-buffer used to store the g(k, q) e-ph matrix elements
 ! for all the k-point treated by this MPI rank.
 ! Increasing the buffer size increases the memory requirements
 ! but it leads to better performance as the number of IO operations is decreased.
 ! TODO: Should compute it on the basis of my_nkpt and my_nqpt
 qbuf_size = 16
 call wrtout(std_out, sjoin(" Begin computation of e-ph matrix elements with qbuf_size:", itoa(qbuf_size)), pre_newlines=1)
 call pstat_proc%print(_PSTAT_ARGS_)

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkibz = ebands%nkpt; mband = ebands%mband

 ! FFT meshes
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)
 ABI_CHECK(dvdb%has_fields("pot1", msg), msg)

 ! Activate parallelism over perturbations at the level of the DVDB
 call gstore%set_perts_distrib(cryst, dvdb, my_npert)

 ! Prepare Fourier interpolation of DFPT potentials.
 comm_rpt = xmpi_comm_self
 !comm_rpt = bqs_comm%value

 ! qmap_symrec gives the mapping gstore%ibz --> dvdb%ibz
 call dvdb%need_ftinterp(gstore%nqibz, gstore%qibz, gstore%qptopt, qmap_symrec, need_ftinterp)
 ABI_FREE(qmap_symrec)
 !need_ftinterp = .True.

 if (.not. need_ftinterp .and. dtset%eph_use_ftinterp /= 0) then
   ABI_WARNING("Enforcing FT interpolation for q-points even if it's not strictly needed.")
   need_ftinterp = .True.
 end if

 if (need_ftinterp) then
   call wrtout(units, " Cannot find all IBZ q-points in the DVDB --> Activating Fourier interpolation.")
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, gstore%qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)
 else
   call wrtout(units, " DVDB file contains all q-points in the IBZ --> Reading DFPT potentials from file.")
 end if

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical image of the k/k+q wavevectors treated by this MPI rank are stored.
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 call gstore%fill_bks_mask(mband, nkibz, nsppol, bks_mask)

 ! Impose istwfk = 1 for all k-points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1
 ecut = dtset%ecut

 call wfd%init(cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
               nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print([std_out], header="Wavefunctions for GSTORE calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))
 call pstat_proc%print(_PSTAT_ARGS_)

 ! one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2, 3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! that will be used to symmetrize the wavefunctions in G-space.
 call gstore%get_mpw_gmax(ecut, mpw, gmax)

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 usecprj = 0
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1    ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2       ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0  ! gvnlx1 is output
 ABI_MALLOC(gvnlx1, (2, usevnl))
 ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call gs_ham_kq%init(psps, pawtab, nspinor, nsppol, nspden, natom, &
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg, &
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab, &
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 ! Allocate vlocal. Note nvloc
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_MALLOC(vlocal, (n4, n5, n6, gs_ham_kq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(displ_cart_qbz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_cart_qibz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red_qbz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red_qibz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(pheigvec_qbz, (2, 3, cryst%natom, 3*cryst%natom))
 ABI_MALLOC(pheigvec_qibz, (2, 3, cryst%natom, 3*cryst%natom))
 ABI_CALLOC(dummy_vtrial, (nfftf, nspden))

 ! Open GSTORE file, and read table used for restarting.
 NCF_CHECK(nctk_open_modify(root_ncid, gstore%path, gstore%comm))

 ! integer scalars
 ncerr = nctk_def_iscalars(root_ncid, [character(len=nctk_slen) :: &
   "used_ftinterp" &
 ])
 NCF_CHECK(ncerr)

 ABI_MALLOC(done_qbz_spin, (gstore%nqbz, nsppol))
 NCF_CHECK(nf90_get_var(root_ncid, nctk_idname(root_ncid, "gstore_done_qbz_spin"), done_qbz_spin))

 gstore%wfk0_path = wfk0_path

 !if (my_rank == master) then
 ii = merge(1, 0, need_ftinterp)
 ncerr = nctk_write_iscalars(root_ncid, [character(len=nctk_slen) :: &
   "used_ftinterp"], &
   [ii &
 ])
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_wfk0_path"), trim(gstore%wfk0_path)))
 !end if

 if (my_rank == master) call gstore%print([std_out])

 ndone = count(done_qbz_spin == 1)

 ! NB: Write phonon data here as we are not guaranteed to have all the IBZ q-points
 ! inside the loop over my_iq if filtering has been used.
 if (ndone == 0) then
   call wrtout(std_out, " Computing phonon frequencies and displacements in the IBZ", pre_newlines=1)
   call cwtime(cpu, wall, gflops, "start")

   call xmpi_split_block(gstore%nqibz, gstore%comm, my_nqibz, my_iqibz_inds)
   ABI_MALLOC(buf_wqnu, (natom3, my_nqibz))
   ABI_MALLOC(buf_eigvec_cart, (2, 3, natom, natom3, my_nqibz))

   do ii=1,my_nqibz
     iq_ibz = my_iqibz_inds(ii)
     call ifc%fourq(cryst, gstore%qibz(:, iq_ibz), buf_wqnu(:,ii), displ_cart_qibz, &
                    out_eigvec=buf_eigvec_cart(:,:,:,:,ii))
   end do
   if (nproc > 1 .and. gstore%nqibz >= nproc) then
     NCF_CHECK(nctk_set_collective(root_ncid, root_vid("phfreqs_ibz")))
     NCF_CHECK(nctk_set_collective(root_ncid, root_vid("pheigvec_cart_ibz")))
   end if
   call xmpi_barrier(gstore%comm)

   if (my_nqibz > 0) then
     iq_start = my_iqibz_inds(1)
     ncerr = nf90_put_var(root_ncid, root_vid("phfreqs_ibz"), buf_wqnu, start=[1, iq_start], count=[natom3, my_nqibz])
     NCF_CHECK(ncerr)
     ncerr = nf90_put_var(root_ncid, root_vid("pheigvec_cart_ibz"), buf_eigvec_cart, &
                          start=[1,1,1,1,iq_start], count=[2, 3, natom, natom3, my_nqibz])
     NCF_CHECK(ncerr)
   end if

   ABI_FREE(my_iqibz_inds)
   ABI_FREE(buf_wqnu)
   ABI_FREE(buf_eigvec_cart)
   call cwtime_report(" Phonon computation + output", cpu, wall, gflops)
 else
   call wrtout(std_out, sjoin(" Restarting GSTORE calculation. Found: ", itoa(ndone), " (qpt, spin) entries already computed"))
 end if

 ! Create ddkop object to compute group velocities (if needed)
 call ddkop%init(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 if (gstore%with_vk /= 0 .and. ndone == 0) then
   call wrtout(std_out, " Computing and writing velocity operator matrix elements in the IBZ", pre_newlines=1)
   call wrtout(std_out, " Note that not all the k-points in the IBZ are computed when kfilter is activated!")
   call cwtime(cpu, wall, gflops, "start")

   ! On disk, we have:
   !    nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz"))
   !    nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz")))

   ABI_MALLOC(cgwork, (2, mpw*wfd%nspinor))

   do my_is=1,gstore%my_nspins
     gqk => gstore%gqk(my_is)
     spin = gstore%my_spins(my_is)

     if (gstore%with_vk == 1) then
       ABI_CALLOC(vk_cart_ibz, (3, gqk%nb_k, gstore%nkibz))
     else
       ABI_ERROR("gstore%with_vk 2 not implemented")
     end if

     cnt = 0
     do my_ik=1,gqk%my_nk
       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       kk_bz = gqk%my_kpts(:, my_ik)
       weight_k = gqk%my_wtk(my_ik)

       ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
       trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5,my_ik)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       if (.not. isirr_k) cycle

       ! parallelize inside (q, pert) so that only one proc in the 3D grid
       ! computes vk for this kpt in the BZ and we can use xmpi_sum_master.
       cnt = cnt + 1
       if (gqk%qpt_pert_comm%skip(cnt)) cycle

       kk_ibz = ebands%kptns(:,ik_ibz)
       npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
       call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk_ibz, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

       select case (gstore%with_vk)
       case (1)
         do in_k=1,gqk%nb_k
           band_k = in_k + gqk%bstart_k - 1
           call wfd%copy_cg(band_k, ik_ibz, spin, cgwork)
           vk = ddkop%get_vdiag(ebands%eig(band_k, ik_ibz, spin), istwf_k, npw_k, wfd%nspinor, cgwork, cwaveprj0)
           vk_cart_ibz(:, in_k, ik_ibz) = vk
         end do

       case (2)
         ABI_ERROR("with_vk 2")
         !do in_k=1,nb_k
         !  band_k = in_k + bstart_k - 1
         !end do
       end select
     end do ! my_ik

     call xmpi_sum_master(vk_cart_ibz, master, gqk%comm%value, ierr)
     !if (gqk%comm%me == master) then
       NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))
       NCF_CHECK(nf90_put_var(spin_ncid, spin_vid("vk_cart_ibz"), vk_cart_ibz))
     !end if
     ABI_SFREE(vk_cart_ibz)
   end do ! my_is

   ABI_FREE(cgwork)
   call cwtime_report(sjoin(" Computation of v_k group velocities with with_vk:", itoa(gstore%with_vk)), cpu, wall, gflops)
 end if

 call wrtout(std_out, " Begin computation of e-ph matrix elements...", pre_newlines=1)

 ! TODO: Exchange the q/k loops so that one can reduce the number of calls to _k dependente routines
 ! and precompute u_nk(r) so that we can save one FFT every time we compute <g|v^1_{kappa,a}(r)|u_nk(r)>
 ! The price to pay is an increase in the number of calls to get_ftqbz but for small systems this part does not dominate
 ! Alternatively, one cah have two versions that will be invoked depending on my_nk, my_nq

 ! Here we decide if the q-points can be reduced to the IBZ(k)
 if (dtset%gstore_use_lgk /= 0) then
   call wrtout(units, " Only q-points in the IBZ_k will be computed.")
 else if (dtset%gstore_use_lgq /= 0) then
   call wrtout(units, " Only k-points in the IBZ_q will be computed.")
 else
   call wrtout(units, " Little group operations won't be used")
 end if

 ! if PAW, one has to solve a generalized eigenproblem
 gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1

 ! Loop over my spins.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is)
   my_npert = gqk%my_npert
   NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

   ! Allocate workspace for wavefunctions using mpw and nb
   ! FIXME: Should be allocated with npw_k and npw_kw but one has to change wfd_sym_ug_kg to get rid of mpw
   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq

   ABI_MALLOC(kets_k, (2, mpw*nspinor, nb_k))
   ABI_MALLOC(bras_kq, (2, mpw*nspinor, nb_kq))
   ABI_MALLOC(h1_kets_kq, (2, mpw*nspinor, nb_kq))

   ABI_MALLOC(iq_buf, (2, qbuf_size))
   ABI_MALLOC(gkq_atm, (2, nb_kq, nb_k, natom3))
   ABI_MALLOC(gkq_nu, (2, nb_kq, nb_k, natom3))

   ! Inside the loops we compute gkq_nu(2, nb_kq, nb_kq, natom3)
   ABI_MALLOC_OR_DIE(my_gbuf, (gqk%cplex, nb_kq, nb_k, natom3, gqk%my_nk, qbuf_size), ierr)
   call pstat_proc%print(_PSTAT_ARGS_)

   ! Compute the little group of the k-point so that we can compute g(k,q) only for q in the IBZ_k
   if (dtset%gstore_use_lgk /= 0) then
     timrev_k = kpts_timrev_from_kptopt(ebands%kptopt)
     ABI_MALLOC(lg_myk, (gqk%my_nk))
     do my_ik=1,gqk%my_nk
       kk_bz = gqk%my_kpts(:, my_ik)
       call lg_myk(my_ik)%init(cryst, kk_bz, timrev_k, gstore%nqbz, gstore%qbz, gstore%nqibz, gstore%qibz, xmpi_comm_self)
     end do
   end if

   ! Loop over my set of q-points
   prev_iqbz = -1
   do my_iq=1,gqk%my_nq
     print_time = my_rank == 0 .and. (my_iq <= LOG_MODQ .or. mod(my_iq, LOG_MODQ) == 0)
     if (print_time) call cwtime(cpu_q, wall_q, gflops_q, "start")
     iq_bz = gqk%my_q2bz(my_iq)

     call gqk%myqpt(my_iq, gstore, weight_q, qq_bz)
     qq_is_gamma = sum(qq_bz**2) < tol14
     qq_bz_string = ktoa(qq_bz)

     ! Handle possible restart.
     if (done_qbz_spin(iq_bz, spin) == 1) then
       call wrtout(std_out, sjoin(" iq_bz:", itoa(iq_bz), ", spin: ", itoa(spin), " already computed --> skipping iteration"))
       cycle
     end if

     ! Compute the little group of the q-point so that we can compute g(k,q) only for k in the IBZ_q
     if (dtset%gstore_use_lgq /= 0) then
       timrev_q = kpts_timrev_from_kptopt(gstore%qptopt)
       call lg_myq%init(cryst, qq_bz, timrev_q, gstore%nkbz, gstore%kbz, gstore%nkibz, gstore%kibz, xmpi_comm_self)
     end if

     iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
     trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5,my_iq)
     ! Don't test if umklapp == 0 because we use the periodic gauge:
     !
     !      phfreq(q+G) = phfreq(q) and eigvec(q) = eigvec(q+G)
     !
     !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     isirr_q = (isym_q == 1 .and. trev_q == 0)
     qq_ibz = gstore%qibz(:, iq_ibz)

     nskip_tetra_kq = 0
     iqbuf_cnt = 1 + mod(my_iq - 1, qbuf_size)
     iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]

     if (iq_ibz /= prev_iqbz) then
       ! Get phonon frequencies and eigenvectors for the corresponding q-point in the IBZ.
       call ifc%fourq(cryst, qq_ibz, phfr_qq, displ_cart_qibz, out_displ_red=displ_red_qibz, out_eigvec=pheigvec_qibz)
       prev_iqbz = iq_ibz
     end if

     if (isirr_q) then
       displ_cart_qbz = displ_cart_qibz; displ_red_qbz = displ_red_qibz; pheigvec_qbz = pheigvec_qibz
     else
       ! Rotate phonon eigenvectors from q_ibz to q_bz.
       ! This part is needed to enforce the gauge in the ph eigenvectors, including e(-q) = e(q)^*
       call pheigvec_rotate(cryst, qq_ibz, isym_q, trev_q, pheigvec_qibz, pheigvec_qbz, displ_cart_qbz, &
                            displ_red_qbz=displ_red_qbz)
     end if

     if (need_ftinterp) then
       ! Fourier interpolation.
       call dvdb%get_ftqbz(qq_bz, cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)
     else
       ! Read and reconstruct the dvscf potentials for qpt and my_npert perturbations.
       db_iqpt = dvdb%findq(qq_ibz)
       ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qq_bz), "in DVDB file."))
       ! The first entry in qbz2dvdb gives the index in dvdb%qpts.
       ! The other entries in mapc_qq are OK as they refer to symmetries.
       qbz2dvdb = gqk%my_q2ibz(:, my_iq); qbz2dvdb(1) = db_iqpt
       call dvdb%readsym_qbz(cryst, qq_bz, qbz2dvdb, cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)
     end if

     ! Allocate vlocal1 with correct cplex. Note nvloc and my_npert.
     ABI_MALLOC(vlocal1, (cplex*n4, n5, n6, gs_ham_kq%nvloc, my_npert))

     ! Set up local potential vlocal1 with proper dimensioning from vtrial1 taking into account the spin.
     do my_ip=1,my_npert
       call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_ham_kq%nvloc,&
                                  pawfgr, mpi_enreg, dummy_vtrial, v1scf(:,:,:,my_ip), vlocal, vlocal1(:,:,:,:,my_ip))
     end do

     ! Continue to initialize the GS Hamiltonian
     call gs_ham_kq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

     ! Loop over my k-points
     do my_ik=1,gqk%my_nk
       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       kk_bz = gqk%my_kpts(:, my_ik)
       weight_k = gqk%my_wtk(my_ik)
       kk_string = ktoa(kk_bz)

       ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
       trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5,my_ik)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       kk_ibz = ebands%kptns(:,ik_ibz)

       ! Set entry to zero. Important as there are cycle instructions inside these loops
       ! and we don't want random numbers written to disk.
       my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = zero

       if (dtset%gstore_use_lgk /= 0) then
         ii = lg_myk(my_ik)%findq_ibzk(qq_bz)
         if (ii == -1) then
           call wrtout(std_out, sjoin(" iq_bz:", itoa(iq_bz), qq_bz_string, " not in IBZ_k --> skipping iteration"))
           cycle
           ! TODO: Check fillvalue (should be zero)
         end if
       end if

       if (dtset%gstore_use_lgq /= 0) then
         ii = lg_myq%findq_ibzk(kk_bz)
         if (ii == -1) then
           call wrtout(std_out, sjoin(" my_ik:", itoa(my_ik), kk_string, " not in IBZ_q --> skipping iteration"))
           cycle
           ! TODO: Check fillvalue (should be zero)
         end if
       end if

       ! =========================================
       ! Find symmetrical image of k+q in the kIBZ
       ! =========================================
       kq_bz = kk_bz + qq_bz
       if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kq_bz, indkk_kq) /= 0) then
         write(msg, '(3a)' ) &
          "Cannot find k+q in k-mesh", ch10, 'Check your WFK file and the (k,q) point input variables.'
          ABI_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1, 1); isym_kq = indkk_kq(2, 1)
       trev_kq = indkk_kq(6, 1); g0_kq = indkk_kq(3:5, 1)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)

       ! If we have used the KERANGE trick, we may have k or k+q points with just one G component set to zero
       ! so we skip this transition immediately. This should happen only if fsewin > sigma_erange.
       if (wfd%npwarr(ik_ibz) == 1 .or. wfd%npwarr(ikq_ibz) == 1) cycle

       if (gstore%kfilter == "fs_tetra") then
         ! Check tetra delta(e_{k+q}) and cycle if all the weights at k+q are zero.
         if (all(abs(gstore%delta_ef_kibz_spin(:, ikq_ibz, spin)) == zero)) then
           nskip_tetra_kq = nskip_tetra_kq + 1; cycle
         end if
       end if

       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       ! TODO: these routines now should allocate wavefunctions as
       !real(dp),intent(out) :: cgs_kbz(2, npw_k*self%nspinor, nband)
       call wfd%sym_ug_kg(ecut, kk_bz, kk_ibz, gqk%bstart_k, nb_k, spin, mpw, gqk%my_k2ibz(:, my_ik), cryst, &
                          work_ngfft, work, istwf_k, npw_k, kg_k, kets_k)

       ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kq_bz, kq_ibz, gqk%bstart_kq, nb_kq, spin, mpw, indkk_kq(:,1), cryst, &
                          work_ngfft, work, istwf_kq, npw_kq, kg_kq, bras_kq)

       call gs_ham_kq%eph_setup_k("k" , kk_bz, istwf_k, npw_k, kg_k, dtset, cryst, psps, &
                                  nkpg_k, kpg_k, ffnl_k, kinpw_k, ph3d_k, gqk%pert_comm%value)

       call gs_ham_kq%eph_setup_k("kq", kq_bz, istwf_k, npw_kq, kg_kq, dtset, cryst, psps, &
                                  nkpg_kq, kpg_kq, ffnl_kq, kinpw_kq, ph3d_kq, gqk%pert_comm%value)

       ABI_MALLOC(gs1c_kq, (2, npw_kq*nspinor*((sij_opt+1)/2)))

       ! Loop over my atomic perturbations and compute gkq_atm.
       gkq_atm = zero
       do my_ip=1,my_npert
         idir = dvdb%my_pinfo(1, my_ip); ipert = dvdb%my_pinfo(2, my_ip); ipc = dvdb%my_pinfo(3, my_ip)

         ! Prepare application of the NL part.
         call rf_ham_kq%init(cplex, gs_ham_kq, ipert, has_e1kbsc=.true.)
         call rf_ham_kq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,my_ip), with_nonlocal=.true.)

         ! Calculate dvscf * psi_k, results stored in h1_kets_kq on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0)
         do in_k=1,nb_k
           band_k = in_k + gqk%bstart_k - 1
           eig0nk = ebands%eig(band_k, ik_ibz, spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0, kets_k(:,:,in_k), cwaveprj0, h1_kets_kq(:,:,in_k), &
                        grad_berry, gs1c_kq, gs_ham_kq, gvnlx1, idir, ipert, [eshift], mpi_enreg, ndat1, optlocal, &
                        optnl, opt_gvnlx1, rf_ham_kq, sij_opt, tim_getgh1c, usevnl)
         end do

         call rf_ham_kq%free()

         ! Calculate <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation. No need to handle istwf_kq because it's always 1.
         !$OMP PARALLEL DO COLLAPSE(2)
         do in_k=1,nb_k
           do im_kq=1,nb_kq
             gkq_atm(:, im_kq, in_k, ipc) = cg_zdotc(npw_kq*nspinor, bras_kq(1,1,im_kq), h1_kets_kq(1,1,in_k))
           end do
         end do

       end do ! my_ip

       ABI_FREE(gs1c_kq)
       ABI_FREE(ffnl_k)
       ABI_FREE(ffnl_kq)
       ABI_FREE(kpg_k)
       ABI_FREE(kpg_kq)
       ABI_FREE(ph3d_k)
       ABI_FREE(ph3d_kq)
       ABI_FREE(kinpw_k)
       ABI_FREE(kinpw_kq)

       ! Collect gkq_atm inside pert_comm so that all procs can operate on the data.
       if (gqk%pert_comm%nproc > 1) call xmpi_sum(gkq_atm, gqk%pert_comm%value, ierr)

       ! Get g in the phonon representation.
       select case (gstore%gmode)
       case (GSTORE_GMODE_PHONON)
         call ephtk_gkknu_from_atm(nb_kq, nb_k, 1, natom, gkq_atm, phfr_qq, displ_red_qbz, gkq_nu)

         ! Save e-ph matrix elements in the buffer.
         select case (gqk%cplex)
         case (1)
           my_gbuf(1,:,:,:, my_ik, iqbuf_cnt) = gkq_nu(1,:,:,:) ** 2 + gkq_nu(2,:,:,:) ** 2
         case (2)
           my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = gkq_nu
         end select

       case (GSTORE_GMODE_ATOM)
         ! Save e-ph matrix elements in the buffer.
         select case (gqk%cplex)
         case (1)
           my_gbuf(1,:,:,:, my_ik, iqbuf_cnt) = gkq_atm(1,:,:,:) ** 2 + gkq_atm(2,:,:,:) ** 2
         case (2)
           my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = gkq_atm
         end select

       case default
         ABI_ERROR(sjoin("Invalid gstore%gmode:", gstore%gmode))
       end select

     end do ! my_ik

     ABI_FREE(v1scf)
     ABI_FREE(vlocal1)

     ! Dump buffer
     if (iqbuf_cnt == qbuf_size) then
       call dump_my_gbuf()
     end if

     if (print_time) then
       write(msg,'(2(a,i0),a)')" My q-point [", my_iq, "/", gqk%my_nq, "]"
       call cwtime_report(msg, cpu_q, wall_q, gflops_q); if (my_iq == LOG_MODQ) call wrtout(std_out, "...", do_flush=.True.)
     end if
     call lg_myq%free()
   end do ! my_iq

   ! Dump the remainder.
   if (iqbuf_cnt /= 0) then
     call dump_my_gbuf()
   end if

   ABI_FREE(iq_buf)
   ABI_FREE(my_gbuf)
   ABI_FREE(bras_kq)
   ABI_FREE(kets_k)
   ABI_FREE(h1_kets_kq)
   ABI_FREE(gkq_atm)
   ABI_FREE(gkq_nu)

   if (dtset%gstore_use_lgk /= 0) then
     do my_ik=1,gqk%my_nk
       call lg_myk(my_ik)%free()
     end do
     ABI_FREE(lg_myk)
   end if
 end do ! my_is

 call cwtime_report(" GSTORE computation done", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10) !, comm=gstore%comm)
 !call gstore%print([std_out], header="GSTORE at the end of gstore%compute")

 ! Set gstore_completed to 1 so that we can easily check if restarted is needed.
 !if (my_rank == master) then
   NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_completed"), 1))
 !end if
 NCF_CHECK(nf90_sync(root_ncid))
 NCF_CHECK(nf90_close(root_ncid))
 call xmpi_barrier(gstore%comm)

 ! Output some of the results to ab_out for testing purposes
 call gstore%print_for_abitests(dtset)

 ! Free memory
 ABI_FREE(gvnlx1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(displ_cart_qbz)
 ABI_FREE(displ_cart_qibz)
 ABI_FREE(displ_red_qbz)
 ABI_FREE(displ_red_qibz)
 ABI_FREE(pheigvec_qbz)
 ABI_FREE(pheigvec_qibz)
 ABI_FREE(done_qbz_spin)

 call ddkop%free(); call gs_ham_kq%free(); call wfd%free()
 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)

contains

subroutine dump_my_gbuf()

 ! This function is called inside the double loop over (my_is, my_iq) or when we exit
 ! from the my_iq loop to dump the remainder that is still in the q-buffer,
 ! All the MPI procs in the (kpt_comm x pert_comm) grid shall call this contained routine
 ! as we have side-effects i.e. iqbuf_cnt set to 0.

 ! On disk we have the global arrays:
 !
 !      nctkarr_t("gvals", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq")
 !
 ! while the local MPI buffers are dimensioned as follows:
 !
 !      my_gbuf(2, nb, nb, natom3, gqk%my_nk, qbuf_size)

 ! If parallelism over perturbation is activated, only the procs treating the first perturbation
 ! i.e. the procs treating different k-points for this q are involved in IO
 ! as all the local buffers store results for all natom3 perturbations.

 integer :: ii, iq_bz, iq_glob, my_iq

 if (gqk%coords_qkpb_sumbp(3) /= 0) goto 10 ! Yes, I'm very proud of this GOTO.

 !iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]
 my_iq = iq_buf(1, 1)
 iq_glob = my_iq + gqk%my_qstart - 1

 !print *, "in dump_my_gbuf with start: ", [1, 1, 1, 1, gqk%my_kstart, iq_glob]
 !print *, "                  count; ", [gqk%cplex, gqk%nb_kq, gqk%nb_k, gqk%natom3, gqk%my_nk, iqbuf_cnt]

 ! NB: this is an individual IO operation
 ncerr = nf90_put_var(spin_ncid, spin_vid("gvals"), my_gbuf, &
                      start=[1, 1, 1, 1, gqk%my_kstart, iq_glob], &
                      count=[gqk%cplex, gqk%nb_kq, gqk%nb_k, gqk%natom3, gqk%my_nk, iqbuf_cnt])
 NCF_CHECK(ncerr)

 ! Only one proc sets the entry in done_qbz_spin to 1 for all the q-points in the buffer.
 !if (all(gqk%coords_qkpb_sumbp(2:3) == [0, 0]))  then
   do ii=1,iqbuf_cnt
     iq_bz = iq_buf(2, ii)
     NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_done_qbz_spin"), 1, start=[iq_bz, spin]))
   end do
 !end if

 ! Zero the counter before returning
10 iqbuf_cnt = 0

 NCF_CHECK(nf90_sync(spin_ncid))
 NCF_CHECK(nf90_sync(root_ncid))

end subroutine dump_my_gbuf

integer function root_vid(var_name)
  character(len=*),intent(in) :: var_name
  root_vid = nctk_idname(root_ncid, var_name)
end function root_vid

integer function spin_vid(var_name)
  character(len=*),intent(in) :: var_name
  spin_vid = nctk_idname(spin_ncid, var_name)
end function spin_vid

end subroutine gstore_compute
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_check_qkzone_gmode
!! NAME
!! gstore_check_cplex_qkzone_gmode
!!
!! FUNCTION
!! Perform consistency checks
!!
!! INPUTS
!!
!! SOURCE

integer function gstore_check_cplex_qkzone_gmode(gstore, cplex, qzone, kzone, gmode, kfilter) result(ierr)

!Arguments ------------------------------------
 class(gstore_t),target,intent(in) :: gstore
 integer,intent(in) :: cplex
 character(len=*),intent(in) :: qzone, kzone, gmode
 character(len=*),optional,intent(in) :: kfilter

!Local variables-------------------------------
 integer :: my_is
 type(gqk_t),pointer :: gqk
! *************************************************************************

 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == qzone, sjoin("qzone: ", qzone, "required but got: ", gstore%qzone), ierr)
 ABI_CHECK_NOSTOP(gstore%kzone == kzone, sjoin("kzone: ", kzone, "required but got: ", gstore%kzone), ierr)
 ABI_CHECK_NOSTOP(gstore%gmode == gmode, sjoin("gmode: ", gmode, "required but got: ", gstore%gmode), ierr)
 if (present(kfilter)) then
   ABI_CHECK_NOSTOP(gstore%kfilter == kfilter, sjoin("kfilter: ", kfilter, "required but got: ", gstore%kfilter), ierr)
 end if

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   ABI_CHECK_NOSTOP(gqk%cplex == cplex, sjoin("cplex:", itoa(cplex), "required but got: ", itoa(gqk%cplex)), ierr)
   if (cplex == 1) then
     ABI_CHECK_NOSTOP(allocated(gqk%my_g2), "my_g2 array is not allocated", ierr)
   else if (cplex == 2) then
     ABI_CHECK_NOSTOP(allocated(gqk%my_g), "my_g array is not allocated", ierr)
   end if
 end do

end function gstore_check_cplex_qkzone_gmode
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_from_ncpath
!! NAME
!! gstore_from_ncpath
!!
!! FUNCTION
!!  Reconstruct a gstore_t instance from a netcdf file.
!!
!! INPUTS
!!  path: Path to the GSTORE file.
!!  with_cplex: 1 for |g|^2, 2 for complex g
!!  dtset: Input variables
!!  cryst: crystalline structure
!!  ebands: KS energies
!!  ifc: interatomic force constants
!!  comm: MPI communicator
!!  [with_gmode]
!!  [gvals_name]: "gvals" (default) or "gvals_ks" to read the KS g produced by the GWPT code.
!!  [read_dw]:
!!
!! SOURCE

subroutine gstore_from_ncpath(gstore, path, with_cplex, dtset, cryst, ebands, ifc, comm, &
                              with_gmode, gvals_name, read_dw)  ! optional

!Arguments ------------------------------------
 class(gstore_t),target,intent(out) :: gstore
 character(len=*),intent(in) :: path
 integer,intent(in) :: with_cplex
 type(dataset_type),target,intent(in) :: dtset
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),target,intent(in) :: ebands
 class(ifc_type),target,intent(in) :: ifc
 integer,intent(in) :: comm
 character(len=*),optional,intent(in) :: with_gmode, gvals_name
 logical,optional,intent(in) :: read_dw

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_rank, ncid, spin, spin_ncid, nproc, ierr, fform, max_nb, ib, natom, natom3, varid, ib_m, ib_n
 integer :: max_nq, max_nk, gstore_cplex, ncerr, my_is, my_iq, iq_glob, my_ik, ik_glob, nb_k, nb_kq
 integer :: my_ip, ipert, iq_ibz, iq_bz, isym_q, trev_q, tsign_q, ii
 real(dp) :: cpu, wall, gflops
 character(len=500) :: gvals_name__
 logical :: isirr_q, read_dw__, from_atm_to_nu
 type(hdr_type) :: wfk0_hdr
 type(crystal_t) :: gstore_cryst
 type(gqk_t),pointer :: gqk
!arrays
 integer :: units(2), ibuffer(9), nproc_spin(ebands%nsppol), comm_spin(ebands%nsppol), brange_spin(2, ebands%nsppol), g0_q(3)
 integer,allocatable :: qglob2bz(:,:), qbz2ibz(:,:), kbz2ibz(:,:)
 real(dp) :: qq_ibz(3)
 real(dp),allocatable :: gwork_q(:,:,:,:,:), slice_bb(:,:,:)
 real(dp),allocatable :: phfreqs_ibz(:,:), pheigvec_cart_ibz(:,:,:,:,:), pheigvec_cart_qbz(:,:,:,:)
 real(dp),allocatable :: displ_cart_qbz(:,:,:,:), displ_red_qbz(:,:,:,:), gmn_nu(:,:,:,:)
! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 gvals_name__ = "gvals"; if (present(gvals_name)) gvals_name__ = gvals_name
 read_dw__ = .False.; if (present(read_dw)) read_dw__ = read_dw

 call wrtout(units, sjoin("- Reading e-ph matrix elements from: ", path))

 ! Set basic parameters.
 gstore%comm = comm; gstore%nsppol = dtset%nsppol; gstore%path = path

 ! Get references to other data structures.
 gstore%dtset => dtset; gstore%cryst => cryst; gstore%ebands => ebands; gstore%ifc => ifc; gstore%kibz => ebands%kptns

 natom = cryst%natom; natom3 = cryst%natom * 3

 ABI_CALLOC(gstore%erange_spin, (2, gstore%nsppol))
 ABI_MALLOC(gstore%glob_nk_spin, (gstore%nsppol))
 ABI_MALLOC(gstore%glob_nq_spin, (gstore%nsppol))

 ! =======================================================
 ! Master node reads basic objects and gstore dimensions
 ! =======================================================

 if (my_rank == master) then
   call wrtout([std_out, ab_out], sjoin(" Initializing gstore object from:", path, ch10))
   ABI_CHECK(path /= ABI_NOFILE, "Use getgstore_filepath to specify the path to GSTORE.nc")
   NCF_CHECK(nctk_open_read(ncid, path, xmpi_comm_self))

   call wfk0_hdr%ncread(ncid, fform)
   ABI_CHECK(fform /= 0, sjoin("Error while reading:", path))

   ! Read gstore dimensions
   NCF_CHECK(nctk_get_dim(ncid, "gstore_cplex", gstore_cplex))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_nkibz", gstore%nkibz))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_nkbz", gstore%nkbz))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_nqibz", gstore%nqibz))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_nqbz", gstore%nqbz))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_max_nq", max_nq))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_max_nk", max_nk))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_max_nb", max_nb))

   ! Read gstore variables
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_with_vk"), gstore%with_vk))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qptopt"), gstore%qptopt))

   ! little group variables were added in Abinit v10.5.6.
   gstore%has_used_lgk = 0; gstore%has_used_lgq = 0
   ncerr = nf90_inq_varid(ncid, "gstore_has_used_lgk", varid)
   if (ncerr == nf90_noerr) then
     NCF_CHECK(nf90_get_var(ncid, vid("gstore_has_used_lgk"), gstore%has_used_lgk))
   end if
   ncerr = nf90_inq_varid(ncid, "gstore_has_used_lgq", varid)
   if (ncerr == nf90_noerr) then
     NCF_CHECK(nf90_get_var(ncid, vid("gstore_has_used_lgq"), gstore%has_used_lgq))
   end if
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kzone"), gstore%kzone))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qzone"), gstore%qzone))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kfilter"), gstore%kfilter))
   call replace_ch0(gstore%kzone); call replace_ch0(gstore%qzone); call replace_ch0(gstore%kfilter)

   ! gstore_gtype was added in Abinit v10.5.6
   gstore%gtype = "KS"
   ncerr = nf90_inq_varid(ncid, "gstore_gtype", varid)
   if (ncerr == nf90_noerr) then
     NCF_CHECK(nf90_get_var(ncid, vid("gstore_gtype"), gstore%gtype))
   end if
   call replace_ch0(gstore%gtype)

   if (gvals_name__ == "gvals_ks") then
     call wrtout(units, " Reading KS e-ph matrix elements")
   else if (gvals_name__ == "gvals") then
     if (gstore%gtype == "GWPT") then
       call wrtout(units, " Reading GWPT e-ph matrix elements")
     else
       call wrtout(units, " Reading KS e-ph matrix elements")
     end if
   end if

   ! gstore_gmode was added in Abinit v10.1.2
   gstore%gmode = GSTORE_GMODE_PHONON
   ncerr = nf90_inq_varid(ncid, "gstore_gmode", varid)
   if (ncerr == nf90_noerr) then
     NCF_CHECK(nf90_get_var(ncid, vid("gstore_gmode"), gstore%gmode))
     call replace_ch0(gstore%gmode)
   end if

   ! gstore_ngqpt was added during the 78_eph/m_varpeq.f90 module development
   gstore%ngqpt(:) = 0
   ncerr = nf90_inq_varid(ncid, "gstore_ngqpt", varid)
   if (ncerr == nf90_noerr) then
     NCF_CHECK(nf90_get_var(ncid, vid("gstore_ngqpt"), gstore%ngqpt))
   endif

   NCF_CHECK(nf90_get_var(ncid, vid("gstore_wfk0_path"), gstore%wfk0_path))
   call replace_ch0(gstore%wfk0_path)

   ABI_MALLOC(gstore%qibz, (3, gstore%nqibz))
   ABI_MALLOC(gstore%wtq, (gstore%nqibz))
   ABI_MALLOC(gstore%qbz, (3, gstore%nqbz))
   ABI_MALLOC(gstore%kbz, (3, gstore%nkbz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_brange_spin"), brange_spin))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_erange_spin"), gstore%erange_spin))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qibz"), gstore%qibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_wtq"), gstore%wtq))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qbz"), gstore%qbz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kbz"), gstore%kbz))

   ABI_MALLOC(qbz2ibz, (6, gstore%nqbz))
   ABI_MALLOC(gstore%kbz2ibz, (6, gstore%nkbz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qbz2ibz"), qbz2ibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kbz2ibz"), gstore%kbz2ibz))

   ABI_MALLOC(qglob2bz, (max_nq, gstore%nsppol))
   ABI_MALLOC(gstore%kglob2bz, (max_nk, gstore%nsppol))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qglob2bz"), qglob2bz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kglob2bz"), gstore%kglob2bz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_glob_nq_spin"), gstore%glob_nq_spin))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_glob_nk_spin"), gstore%glob_nk_spin))

   ! Read optional variables:
   if (gstore%kfilter == "fs_tetra") then
     ABI_MALLOC(gstore%delta_ef_kibz_spin, (max_nb, gstore%nkibz, gstore%nsppol))
     NCF_CHECK(nf90_get_var(ncid, vid("gstore_delta_ef_kibz_spin"), gstore%delta_ef_kibz_spin))
   end if

   NCF_CHECK(nf90_close(ncid))

   ! =================
   ! Consistency check
   ! =================
   gstore_cryst = wfk0_hdr%get_crystal()
   if (cryst%compare(gstore_cryst, header=" Comparing input crystal with the one from GSTORE file") /= 0) then
     ABI_ERROR("Crystal structure from input and GSTORE do not agree! Check messages above!")
   end if
   call gstore_cryst%free()

   if (gstore_cplex == 1 .and. with_cplex == 2) then
     ABI_ERROR("GSTORE file contains |g|^2 while ABINIT needs complex g. Recompute GSTORE file with gstore_cplex 2")
   end if
    if (read_dw__ .and. gstore_cplex == 1) then
       ABI_ERROR("read_dw requires complex g. Recompute GSTORE file with gstore_cplex 2")
    end if
    if (read_dw__ .and. gstore%gmode == GSTORE_GMODE_PHONON) then
       ABI_ERROR('read_dw requires g in the atom representation. Recompute GSTORE file with gstore_mode "atom"')
    end if
 end if ! master

 if (nproc > 1) then
   ! Broadcast the header.
   call wfk0_hdr%bcast(master, my_rank, comm)

   ! Broadcast dimensions.
   if (my_rank == master) then
     ibuffer = [gstore_cplex, gstore%nkibz, gstore%nkbz, gstore%nqibz, gstore%nqbz, gstore%with_vk, max_nq, max_nk, max_nb]
   end if
   call xmpi_bcast(ibuffer, master, comm, ierr)

   if (my_rank /= master) then
     ! Other MPI procs need to store dims and allocate memory before the bcast.
     gstore_cplex = ibuffer(1)
     gstore%nkibz = ibuffer(2)
     gstore%nkbz = ibuffer(3)
     gstore%nqibz = ibuffer(4)
     gstore%nqbz = ibuffer(5)
     gstore%with_vk = ibuffer(6)
     max_nq = ibuffer(7)
     max_nk = ibuffer(8)
     max_nb = ibuffer(9)

     ABI_MALLOC(gstore%qibz, (3, gstore%nqibz))
     ABI_MALLOC(gstore%wtq, (gstore%nqibz))
     ABI_MALLOC(gstore%qbz, (3, gstore%nqbz))
     ABI_MALLOC(gstore%kbz, (3, gstore%nkbz))
     ABI_MALLOC(qbz2ibz, (6, gstore%nqbz))
     ABI_MALLOC(gstore%kbz2ibz, (6, gstore%nkbz))
     ABI_MALLOC(qglob2bz, (max_nq, gstore%nsppol))
     ABI_MALLOC(gstore%kglob2bz, (max_nk, gstore%nsppol))
   end if

   call xmpi_bcast(gstore%has_used_lgk, master, comm, ierr)
   call xmpi_bcast(gstore%has_used_lgq, master, comm, ierr)
   call xmpi_bcast(gstore%qptopt, master, comm, ierr)
   call xmpi_bcast(gstore%has_used_lgk, master, comm, ierr)
   call xmpi_bcast(gstore%has_used_lgq, master, comm, ierr)
   call xmpi_bcast(gstore%kzone, master, comm, ierr)
   call xmpi_bcast(gstore%qzone, master, comm, ierr)
   call xmpi_bcast(gstore%kfilter, master, comm, ierr)
   call xmpi_bcast(gstore%gmode, master, comm, ierr)
   call xmpi_bcast(gstore%gtype, master, comm, ierr)
   call xmpi_bcast(gstore%ngqpt, master, comm, ierr)
   call xmpi_bcast(gstore%wfk0_path, master, comm, ierr)
   call xmpi_bcast(brange_spin, master, comm, ierr)
   call xmpi_bcast(gstore%erange_spin, master, comm, ierr)
   call xmpi_bcast(gstore%qibz, master, comm, ierr)
   call xmpi_bcast(gstore%wtq, master, comm, ierr)
   call xmpi_bcast(gstore%qbz, master, comm, ierr)
   call xmpi_bcast(gstore%kbz, master, comm, ierr)
   call xmpi_bcast(qbz2ibz, master, comm, ierr)
   call xmpi_bcast(gstore%kbz2ibz, master, comm, ierr)
   call xmpi_bcast(qglob2bz, master, comm, ierr)
   call xmpi_bcast(gstore%kglob2bz, master, comm, ierr)
   call xmpi_bcast(gstore%glob_nq_spin, master, comm, ierr)
   call xmpi_bcast(gstore%glob_nk_spin, master, comm, ierr)

   if (gstore%kfilter == "fs_tetra") then
     if (my_rank /= master) then
       ABI_MALLOC(gstore%delta_ef_kibz_spin, (max_nb, gstore%nkibz, gstore%nsppol))
     end if
     call xmpi_bcast(gstore%delta_ef_kibz_spin, master, comm, ierr)
   end if
 end if

 ! Consistency checl
 call wfk0_hdr%vs_dtset(dtset); call wfk0_hdr%free()

 ! Distribute spins, create indirect mapping to spin index and init gstore%brange_spin
 call gstore%distribute_spins__(ebands%mband, brange_spin, nproc_spin, comm_spin, comm)

 ! Compute krank
 call gstore%krank_ibz%from_kptrlatt(gstore%nkibz, gstore%kibz, ebands%kptrlatt, compute_invrank=.False.)

 call gstore%set_mpi_grid__(with_cplex, nproc_spin, comm_spin)

 ! At this point, we have the Cartesian grid (one per spin if any) and we can finally allocate and distribute other arrays.
 call gstore%malloc__(with_cplex, max_nq, qglob2bz, max_nk, gstore%kglob2bz, qbz2ibz, gstore%kbz2ibz)

 call xmpi_comm_free(comm_spin)

 ! Now we read the big arrays with MPI-IO and hdf5 groups.
 ! Note the loop over spin as each gqk has its own dimensions.
 ! Recall the shape of the arrays:
 !
 ! In memory, we have allocated:
 !
 !    my_g(my_npert, nb_kq, my_nq, nb_k, my_nk) if with_cplex == 2 (complex array)
 !
 ! or
 !
 !    my_g2(my_npert, nb, my_nq, nb, my_nk) if with_cplex == 1
 !
 ! On disk, we have:
 !
 !   nctkarr_t("gvals", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq")
 !
 ! so gstore_cplex == 1 and with_cplex == 2 are not compatible.
 !
 call cwtime(cpu, wall, gflops, "start")

 ! ===================================================
 ! Load phonon frequencies and eigenvectors in the IBZ
 ! ===================================================

 ! nctkarr_t("phfreqs_ibz", "dp", "natom3, gstore_nqibz")
 ! nctkarr_t("pheigvec_cart_ibz", "dp", "two, three, natom, natom3, gstore_nqibz")
 ABI_MALLOC(phfreqs_ibz, (natom3, gstore%nqibz))
 ABI_MALLOC(pheigvec_cart_ibz, (2, 3, cryst%natom, cryst%natom * 3, gstore%nqibz))
 ABI_MALLOC(pheigvec_cart_qbz, (2, 3, cryst%natom, cryst%natom * 3))
 ABI_MALLOC(displ_cart_qbz, (2, 3, cryst%natom, cryst%natom * 3))
 ABI_MALLOC(displ_red_qbz, (2, 3, cryst%natom, natom3))

 NCF_CHECK(nctk_open_read(ncid, gstore%path, gstore%comm))

 if (nproc > 1) then
   NCF_CHECK(nctk_set_collective(ncid, vid("phfreqs_ibz")))
   NCF_CHECK(nctk_set_collective(ncid, vid("pheigvec_cart_ibz")))
 end if
 NCF_CHECK(nf90_get_var(ncid, vid("phfreqs_ibz"), phfreqs_ibz))
 NCF_CHECK(nf90_get_var(ncid, vid("pheigvec_cart_ibz"), pheigvec_cart_ibz))
 NCF_CHECK(nf90_close(ncid))

 ! =========================
 ! Load e-ph matrix elements
 ! =========================

 from_atm_to_nu = .False.
 if (present(with_gmode)) then
   ! The only conversion I can think of is: atom --> phonon.
   if (gstore%gmode /= with_gmode) then
     if (gstore%gmode == GSTORE_GMODE_ATOM .and. with_gmode == GSTORE_GMODE_PHONON) then
       from_atm_to_nu = .True.; gstore%gmode = GSTORE_GMODE_PHONON ! Change gstore%gmode here
     else
       ABI_ERROR(sjoin("Conversion from gstore%gmode: ", gstore%gmode, "to:", with_gmode, " is not yet supported"))
     end if
   end if
 end if

 do spin=1,gstore%nsppol
   my_is = gstore%spin2my_is(spin)

   if (my_is /= 0) then
     gqk => gstore%gqk(my_is)
     nb_k = gqk%nb_k; nb_kq = gqk%nb_kq

     ABI_MALLOC(gqk%my_wnuq, (gqk%my_npert, gqk%my_nq))
     ABI_MALLOC(gqk%my_displ_cart, (2, 3, cryst%natom, gqk%my_npert, gqk%my_nq))

     NCF_CHECK(nctk_open_read(ncid, gstore%path, gqk%comm%value))
     NCF_CHECK(nf90_inq_ncid(ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

     ! gstore_cplex defines the data on disk while cplex defines what we want to store in memory
     ABI_MALLOC_OR_DIE(gwork_q, (gstore_cplex, nb_kq, nb_k, gqk%natom3, gqk%glob_nk), ierr)
     ABI_MALLOC(slice_bb, (gstore_cplex, nb_kq, nb_k))

     ! Read my_gq0nm_atm matrix elements for DW.
     if (read_dw__) then
        ! Find the index of q = 0.
        iq_glob = -1
        do ii=1, gstore%glob_nq_spin(spin)
          iq_bz = qglob2bz(ii, spin)
          if (sum(gstore%qbz(:, iq_bz)**2) < tol4) then
            iq_glob = ii; exit
          end if
        end do
        ABI_CHECK_INEQ(iq_glob, -1, "Cannot finq q=0 in g(k,q)!")
        call wrtout(std_out, sjoin(" Reading g_atm(k,q=0) for Debye-Waller with iq_glob:", itoa(iq_glob)))

       ! Read q-slice (individual IO)
       ncerr = nf90_get_var(spin_ncid, spin_vid(gvals_name__), gwork_q, start=[1, 1, 1, 1, 1, iq_glob])
       NCF_CHECK(ncerr)

       ! Allocate my_gq0nm_atm and transfer data. Note TRANSPOSITION in (m, n) indices.
       ABI_MALLOC(gqk%my_gq0nm_atm, (nb_k, nb_kq, natom3, gqk%my_nk))  ! nb_k, nb_kq
        do my_ik=1,gqk%my_nk
          ik_glob = my_ik + gqk%my_kstart - 1
          do ib_m=1,nb_kq
            do ib_n=1,nb_k
              gqk%my_gq0nm_atm(ib_n,ib_m,:,my_ik) = gwork_q(1,ib_m,ib_n,:,ik_glob) + j_dpc * gwork_q(2,ib_m,ib_n,:,ik_glob)
            end do
          end do
        end do
     end if ! read_dw__

     if (from_atm_to_nu) then
       ABI_MALLOC(gmn_nu, (2, nb_kq, nb_k, 3*natom))
     end if

     ! Read my e-ph matrix elements.
     do my_iq=1,gqk%my_nq
        iq_glob = my_iq + gqk%my_qstart - 1

        !call wrtout(std_out, " Computing and storing phonons in the full BZ by rotating the data in the IBZ...")
        iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
        trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
        !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
        isirr_q = (isym_q == 1 .and. trev_q == 0)
        tsign_q = 1; if (trev_q == 1) tsign_q = -1
        qq_ibz = gstore%qibz(:, iq_ibz)
        call pheigvec_rotate(cryst, qq_ibz, isym_q, trev_q, pheigvec_cart_ibz(:,:,:,:,iq_ibz), pheigvec_cart_qbz, displ_cart_qbz, &
                             displ_red_qbz=displ_red_qbz)

        ! Save my frequencies and my phonon displacements.
        gqk%my_wnuq(:, my_iq) = phfreqs_ibz(gqk%my_pertcases(:), iq_ibz)
        gqk%my_displ_cart(:,:,:,:,my_iq) = displ_cart_qbz(:,:,:,gqk%my_pertcases(:))

        ! Read q-slice (individual IO).
        ncerr = nf90_get_var(spin_ncid, spin_vid(gvals_name__), gwork_q, start=[1, 1, 1, 1, 1, iq_glob])
        NCF_CHECK(ncerr)

        do my_ik=1,gqk%my_nk
          ik_glob = my_ik + gqk%my_kstart - 1

          if (from_atm_to_nu) then
            ! Here we convert from g(k,q)_atm to g(k,q)_phonon and replace data in gwork_q at ik_glob
            call ephtk_gkknu_from_atm(nb_kq, nb_k, 1, natom, gwork_q(:,:,:,:, ik_glob), &
                                      phfreqs_ibz(:, iq_ibz), displ_red_qbz, gmn_nu)
            gwork_q(:,:,:,:, ik_glob) = gmn_nu
          end if

          do my_ip=1,gqk%my_npert
            ipert = gqk%my_pertcases(my_ip)
            slice_bb = gwork_q(:,:,:, ipert, ik_glob)

            ! Put data in the right place and handle conversion g --> |g|^2.
            if (with_cplex == gstore_cplex) then
              if (with_cplex == 1) gqk%my_g2(my_ip,:,my_iq,:,my_ik) = slice_bb(1,:,:)
              if (with_cplex == 2) gqk%my_g(my_ip,:,my_iq,:,my_ik) = slice_bb(1,:,:) + j_dpc * slice_bb(2,:,:)
            else
              if (with_cplex == 1 .and. gstore_cplex == 2) then
                gqk%my_g2(my_ip, :, my_iq, :, my_ik) = slice_bb(1,:,:) ** 2 + slice_bb(2,:,:) ** 2
              else
                ABI_ERROR("Conversion from g2 on file to g_cplx in memory is not possible!")
              end if
            end if

          end do
        end do

     end do ! my_iq

     ABI_FREE(gwork_q)
     ABI_FREE(slice_bb)
     ABI_SFREE(gmn_nu)

     ! ==============================================
     ! Read matrix elements of the velocity operator
     ! ==============================================
     if (gstore%with_vk == 1) then
       if (gqk%comm%nproc > 1) then
         NCF_CHECK(nctk_set_collective(spin_ncid, spin_vid("vk_cart_ibz")))
       end if
       NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("vk_cart_ibz"), gqk%vk_cart_ibz))

     else if (gstore%with_vk == 2) then
       if (gqk%comm%nproc > 1) then
         NCF_CHECK(nctk_set_collective(spin_ncid, spin_vid("vkmat_cart_ibz")))
       end if
       NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("vkmat_cart_ibz"), gqk%vkmat_cart_ibz))

       ! Transfer diagonal terms to vk_cart_ibz.
       do ib=1,gqk%nb_k
         gqk%vk_cart_ibz(:, ib, :) = gqk%vkmat_cart_ibz(1, :, ib, ib, :)
       end do
     end if

     NCF_CHECK(nf90_close(ncid))
   end if
 end do ! spin

 ABI_FREE(qglob2bz)
 ABI_FREE(qbz2ibz)
 ABI_FREE(phfreqs_ibz)
 ABI_FREE(pheigvec_cart_ibz)
 ABI_FREE(displ_cart_qbz)
 ABI_FREE(displ_red_qbz)
 ABI_FREE(pheigvec_cart_qbz)

 call xmpi_barrier(gstore%comm)
 call cwtime_report(" gstore_from_ncpath", cpu, wall, gflops)

contains
integer function vid(var_name)
  character(len=*),intent(in) :: var_name
  vid = nctk_idname(ncid, var_name)
end function vid

integer function spin_vid(var_name)
  character(len=*),intent(in) :: var_name
  spin_vid = nctk_idname(spin_ncid, var_name)
end function spin_vid

end subroutine gstore_from_ncpath
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_check_restart
!! NAME
!! gstore_check_restart
!!
!! FUNCTION
!!  Check whether restart from a previous GSTORE.nc is possible.
!!
!! INPUTS
!!
!! SOURCE

subroutine gstore_check_restart(filepath, dtset, nqbz, done_qbz_spin, restart, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: filepath
 type(dataset_type),intent(in) :: dtset
 integer,intent(out) :: nqbz, restart
 integer,allocatable,intent(out) :: done_qbz_spin(:,:)
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_rank, root_ncid, ierr, gstore_completed, gstore_fform, units(2) ! fform
 character(len=500) :: msg
 type(hdr_type) :: gstore_hdr
! *************************************************************************

 my_rank = xmpi_comm_rank(comm); units = [std_out, ab_out]

 restart = 0; nqbz = 0
 if (my_rank == master .and. dtset%eph_restart == 1) then
    if (file_exists(filepath)) then
      ! Use gstore_completed to understand if the previous GSTORE run completed else we need to restart.
      NCF_CHECK(nctk_open_read(root_ncid, filepath, xmpi_comm_self))
      NCF_CHECK(nf90_get_var(root_ncid, root_vid("gstore_completed"), gstore_completed))
      call gstore_hdr%ncread(root_ncid, gstore_fform)
      ABI_CHECK_INEQ(gstore_fform, 0, "Wrong gstore_fform")
      call gstore_hdr%vs_dtset(dtset); call gstore_hdr%free()

      NCF_CHECK(nctk_get_dim(root_ncid, "gstore_nqbz", nqbz))
      ABI_MALLOC(done_qbz_spin, (nqbz, dtset%nsppol))
      NCF_CHECK(nf90_get_var(root_ncid, root_vid("gstore_done_qbz_spin"), done_qbz_spin))
      NCF_CHECK(nf90_close(root_ncid))
      !print *, "done_qbz_spin:", done_qbz_spin

      ! NOTE: done_qbz_spin is dimensioned with the q-points in the BZ but we
      ! FIXME: should set to 1 the q-points in the BZ else we never restart
      ! Perhaps can can set to -1 if q = TS q_ibz if q_ibz is done.
      if (gstore_completed /= 0) then
        ! Previous computation completed, keep a backup of the file and start from scratch.
        restart = 0; done_qbz_spin = 0
        msg = sjoin("Found GSTORE.nc file with all entries already computed.", ch10, &
                    "Will overwrite:", trim(filepath), ch10, "Keeping backup copy in:", strcat(filepath, ".bkp"))
        call wrtout(ab_out, sjoin("WARNING: ", msg))
        ABI_WARNING(msg)
        ! Keep backup copy
        ABI_CHECK(clib_rename(trim(filepath), strcat(filepath, ".bkp")) == 0, "Failed to rename GSTORE file.")
      else
        restart = 1
        call wrtout(units, "- Restarting from a previous GSTORE.nc file")
      end if
    end if
 end if

 call xmpi_bcast(restart, master, comm, ierr)
 call xmpi_bcast(nqbz, master, comm, ierr)
 if (my_rank /= master) then
   ABI_MALLOC(done_qbz_spin, (nqbz, dtset%nsppol))
 end if
 if (nqbz /= 0) call xmpi_bcast(done_qbz_spin, master, comm, ierr)

contains
integer function root_vid(var_name)
  character(len=*),intent(in) :: var_name
  root_vid = nctk_idname(root_ncid, var_name)
end function root_vid

end subroutine gstore_check_restart
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_print_for_abitests
!! NAME
!! gstore_print_for_abitests
!!
!! FUNCTION
!!  Print subset of results to ab_out for testing purposes.
!!  This routine should be called when the output of the GSTORE.nc file is completed
!!
!! INPUTS
!!  dtset<type(dataset_type)>=all input variables for this dataset
!!  [with_ks]: if True, read "gvals" as well as "gvals_ks". This options is used in GWPT
!!   in which "gvals" are the GWPT matrix elements and ""gvals_ks" are the KS ones.
!!
!! SOURCE

subroutine gstore_print_for_abitests(gstore, dtset, with_ks)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 type(dataset_type),intent(in) :: dtset
 logical,optional,intent(in) :: with_ks

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: root_ncid, spin_ncid, gstore_completed, spin, nb, ik_glob, iq_glob, ipc, cplex, ncerr, natom3, varid
 integer :: glob_nq, glob_nk, im_kq, in_k, nb_k, nb_kq ! ib, ik_ibz,
 logical :: with_ks__
 real(dp) :: g2, g2_ks
 character(len=abi_slen) :: gstore_gmode
!arrays
 integer,allocatable :: done_qbz_spin(:,:)
 real(dp),allocatable :: gslice_mn(:,:,:), gslice_ks_mn(:,:,:)  !,vk_cart_ibz(:,:,:), vkmat_cart_ibz(:,:,:,:)
! *************************************************************************

 ! Only master prints to ab_out
 if (xmpi_comm_rank(gstore%comm) /= master) return

 with_ks__ = .False.; if (present(with_ks)) with_ks__ = with_ks

 natom3 = dtset%natom * 3
 NCF_CHECK(nctk_open_read(root_ncid, gstore%path, xmpi_comm_self))

 NCF_CHECK(nf90_get_var(root_ncid, root_vid("gstore_completed"), gstore_completed))
 write(ab_out, *)""
 write(ab_out, "(a,i0)")" gstore_completed: ", gstore_completed

 ABI_MALLOC(done_qbz_spin, (gstore%nqbz, dtset%nsppol))
 NCF_CHECK(nf90_get_var(root_ncid, root_vid("gstore_done_qbz_spin"), done_qbz_spin))
 write(ab_out, "(a,*(i0,1x))")" gstore_done_qbz_spin: ", count(done_qbz_spin == 1)
 ABI_FREE(done_qbz_spin)

 ! gstore_gmode was added in Abinit v10.1.2
 gstore_gmode = GSTORE_GMODE_PHONON
 ncerr = nf90_inq_varid(root_ncid, "gstore_gmode", varid)
 if (ncerr == nf90_noerr) then
   NCF_CHECK(nf90_get_var(root_ncid, root_vid("gstore_gmode"), gstore_gmode))
   call replace_ch0(gstore_gmode)
 end if

 do spin=1,gstore%nsppol
   NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))
   NCF_CHECK(nctk_get_dim(spin_ncid, "glob_nq", glob_nq))
   NCF_CHECK(nctk_get_dim(spin_ncid, "glob_nk", glob_nk))
   NCF_CHECK(nctk_get_dim(spin_ncid, "gstore_cplex", cplex))

   NCF_CHECK(nctk_get_dim(spin_ncid, "nb", nb))
   write(ab_out, "(a,i0)")" gqk%nb: ", nb
   nb_k = nb; nb_kq = nb

   !NCF_CHECK(nctk_get_dim(spin_ncid, "nb_k", nb_k))
   !NCF_CHECK(nctk_get_dim(spin_ncid, "nb_kq", nb_kq))
   !write(ab_out, "(a,i0)")" gqk%nb_kq: ", nb_kq
   !write(ab_out, "(a,i0)")" gqk%nb_k: ", nb_k
   write(ab_out, "(a,i0)")" gqk%glob_nq: ", glob_nq
   write(ab_out, "(a,i0)")" gqk%glob_nk: ", glob_nk

   ! FIXME
   ! Handle the output of group velocities. On disk, we have:
   !
   !    nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz"))
   ! or
   !    nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz")))

   !select case (gstore%with_vk)
   !case (0)
   !  continue

   !case (1)
   !  ABI_MALLOC(vk_cart_ibz, (2, 3, nb))
   !  write(ab_out,"(a)") " Group velocities v_nk in Cartesian coordinates and atomic units:"
   !  do ik_ibz=1,gstore%nkibz
   !    ! Only the first and the last k-points are written.
   !    if (ik_ibz /= 1 .and. ik_ibz /= gstore%nkibz) cycle
   !    NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("vk_cart_ibz"), vk_cart_ibz, start=[1,1,ik_ibz], count=[3,nb,1]))
   !    write(ab_out, "(a)")sjoin(" For k-point:", ktoa(gstore%kibz(:,ik_ibz)), ", spin", itoa(spin))
   !    do ib=1,nb
   !      write(ab_out, "(6es16.6)") vk_cart_ibz(:,:,ib)
   !    end do
   !  end do
   !  write(ab_out, "(a)")" "
   !  ABI_FREE(vk_cart_ibz)

   !case (2)
   !  write(ab_out, "(a)")" TEXT Output of vkmat is not coded yet!"
   !end select

   ! Handle the output of the e-ph matrix elements. On disk we have the global array:
   !
   !    nctkarr_t("gvals", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq")

   !cplex = dtset%gstore_cplex
   ABI_MALLOC(gslice_mn, (cplex, nb_kq, nb_k))
   ABI_MALLOC(gslice_ks_mn, (cplex, nb_kq, nb_k))

   write(ab_out,"(a)") " E-PH matrix elements:"
   if (with_ks__) then
     write(ab_out, "(1x,5(a5,1x),2(a16))") "iq","ik", "mode", "im_kq", "in_k", "|g^SE|^2 in Ha^2", "|g^KS|^2 in Ha^2"
   else
     write(ab_out, "(1x,5(a5,1x),a16)") "iq","ik", "mode", "im_kq", "in_k", "|g|^2 in Ha^2"
   end if

   do iq_glob=1,glob_nq
     if (iq_glob /= 1 .and. iq_glob /= glob_nq) cycle  ! Write the first and the last q-point.
     do ik_glob=1,glob_nk
       if (ik_glob /= 1 .and. ik_glob /= glob_nk) cycle ! Write the first and the last k-point.
       do ipc=1,natom3
         if (ipc /= 4 .and. ipc /= natom3) cycle ! Write the 4th and the last perturbation.
         ncerr = nf90_get_var(spin_ncid, spin_vid("gvals"), gslice_mn, &
                              start=[1,1,1,ipc,ik_glob,iq_glob], count=[cplex,nb_kq,nb_k,1,1,1])
         NCF_CHECK(ncerr)

         ! TODO: Get rid of cplex, write everything using complex and atom representation
         write(ab_out, "(3(a,1x,i0,1x))")" |g(k,q)|^2 in Ha^2 for iq:", iq_glob, "ik:", ik_glob, "mode:", ipc

         if (.not. with_ks__) then
          ! gvals only.
           write(ab_out, "(1x,5(a5,1x),a16)")"iq","ik", "mode", "im_kq", "in_k", "|g|^2 in Ha^2"
           do im_kq=1,nb_kq
             do in_k=1,nb_k
               if (cplex == 1) g2 = gslice_mn(1, im_kq, in_k)
               if (cplex == 2) g2 = gslice_mn(1, im_kq, in_k)**2 + gslice_mn(2, im_kq, in_k)**2
               write(ab_out, "(1x,5(i5,1x),es16.6)") iq_glob, ik_glob, ipc, im_kq, in_k, g2
             end do
           end do
        else
          ! g^SE and g^KS
          ncerr = nf90_get_var(spin_ncid, spin_vid("gvals_ks"), gslice_ks_mn, &
                               start=[1,1,1,ipc,ik_glob,iq_glob], count=[cplex,nb,nb,1,1,1])
          NCF_CHECK(ncerr)
          write(ab_out, "(1x,5(a5,1x),2a16)")"iq","ik", "mode", "im_kq", "in_k", "|g^SE|^2 in Ha^2", "|g^KS|^2 in Ha^2"
          do im_kq=1,nb_kq
            do in_k=1,nb_k
              if (cplex == 1) then
                g2 = gslice_mn(1, im_kq, in_k)
                g2_ks = gslice_ks_mn(1, im_kq, in_k)
              end if
              if (cplex == 2) then
                g2 = gslice_mn(1, im_kq, in_k)**2 + gslice_mn(2, im_kq, in_k)**2
                g2_ks = gslice_ks_mn(1, im_kq, in_k)**2 + gslice_ks_mn(2, im_kq, in_k)**2
              end if
              write(ab_out, "(1x,5(i5,1x),2(es16.6))") iq_glob, ik_glob, ipc, im_kq, in_k, g2, g2_ks
            end do
          end do
        end if

       end do
     end do
   end do

   ABI_FREE(gslice_mn)
   ABI_FREE(gslice_ks_mn)
 end do

 NCF_CHECK(nf90_close(root_ncid))

contains
integer function root_vid(var_name)
  character(len=*),intent(in) :: var_name
  root_vid = nctk_idname(root_ncid, var_name)
end function root_vid

integer function spin_vid(var_name)
  character(len=*),intent(in) :: var_name
  spin_vid = nctk_idname(spin_ncid, var_name)
end function spin_vid

end subroutine gstore_print_for_abitests
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_gather
!! NAME
!!  gqk_gather
!!
!! FUNCTION
!!  Gather the MPI-distributed matrix elements for a given k/q-point index.
!!  Once one reciprocal dimension is fixed, the gathering is performed across the other one.
!!
!! INPUTS
!!  mode = String controlling the choice of the fixed dimension ("k" or "q")
!!  fixed_pt = Index of the fixed k/q-point
!!
!! OUTPUT
!!  g_gathered(:,:,:,:) = Matrix elements for the given k/q-point, ordered as
!!      (my_npert, nb, nb, glob_nq) - k-point is fixed
!!      (my_npert, nb, nb, glob_nk) - q-point is fixed
!!
!! SOURCE

subroutine gqk_gather(gqk, mode, fixed_pt, g_gathered)

!Arguments ------------------------------------
!scalars
 class(gqk_t), target, intent(in) :: gqk
 character(len=*),intent(in) :: mode
 integer,intent(in) :: fixed_pt
!arrays
 complex(dp), allocatable, intent(out) :: g_gathered(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: comm, ierr, ipt_glob, ngather, my_ipt, my_ngather, my_ptstart, my_pert, ib, jb
!arrays
 complex(dp), pointer :: my_g(:,:,:,:)
!----------------------------------------------------------------------

 select case (mode)
 case ("k")
   comm = gqk%qpt_comm%value
   ngather = gqk%glob_nq
   my_ngather = gqk%my_nq
   my_ptstart = gqk%my_qstart
   my_g => gqk%my_g(:,:,:,:,fixed_pt)
   ABI_MALLOC(g_gathered, (gqk%my_npert, gqk%nb, gqk%nb, ngather))

 case ("q")
   comm = gqk%kpt_comm%value
   ngather = gqk%glob_nk
   my_ngather = gqk%my_nk
   my_ptstart = gqk%my_kstart
   my_g => gqk%my_g(:,:,fixed_pt,:,:)
   ABI_MALLOC(g_gathered, (gqk%my_npert, gqk%nb, gqk%nb, ngather))

 case default
   ABI_ERROR(sjoin("Gathering MPI-distributed matrix elements, unsupported mode: ", mode))
 end select

 g_gathered(:,:,:,:) = zero
 do my_ipt=1,my_ngather
   ipt_glob = my_ipt + my_ptstart - 1

   ! FIXME: can the rearrangement be done in the select case statement?
   do ib=1,gqk%nb
     do jb=1,gqk%nb
       do my_pert=1,gqk%my_npert
         if (mode == "k") then
           g_gathered(my_pert, jb, ib, ipt_glob) = my_g(my_pert, jb, my_ipt, ib)
         else
           g_gathered(my_pert, jb, ib, ipt_glob) = my_g(my_pert, jb, ib, my_ipt)
         endif
       enddo
     enddo
   enddo
 enddo

 call xmpi_sum(g_gathered, comm, ierr)

end subroutine gqk_gather
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_wannierize_and_write_gwan
!! NAME
!! gstore_wannierize_and_write_gwan
!!
!! FUNCTION
!!  Compute g(R_e,R_ph) from g(k,q) and save results to GWAN.nc file
!!
!! INPUTS
!!
!! SOURCE

subroutine gstore_wannierize_and_write_gwan(gstore, dvdb, dtfil)

!Arguments ------------------------------------
 class(gstore_t),target, intent(in) :: gstore
 type(dvdb_t),intent(in) :: dvdb
 type(datafiles_type),intent(in) :: dtfil

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: nr_e, nr_p, nwan, iwan, jwan, spin, my_is, my_ip, ir, my_ik, my_iq, nb_k, nb_kq
 integer :: my_nk, my_nq, ierr, ik, ikq, my_npert, nwin_k, nwin_kq, ii, jj, band_kq, band_k, ib_k, ib_kq
 !character(len=500) :: msg
 logical :: keep_umats
 type(wan_t),pointer :: wan
 type(gqk_t),pointer :: gqk
!arrays
 integer :: qptrlatt_(3,3), units(2)
 real(dp) :: weight_qq, qpt(3), kpt(3), kq(3), cpu, wall, gflops
 complex(dp),allocatable :: emikr(:), emiqr(:), u_k(:,:), u_kq(:,:), gww_epq(:,:,:,:,:), gww_pk(:,:,:,:), g_bb(:,:), tmp_mat(:,:)
! *************************************************************************

 units = [std_out, ab_out]
 call wrtout(units, " Computing e-ph matrix elements in the Wannier representation...", pre_newlines=1)
 call cwtime(cpu, wall, gflops, "start")

 if (gstore%check_cplex_qkzone_gmode(2, "bz", "bz", "atom", kfilter="none") /= 0) then
   ABI_ERROR("The gstore object is inconsistent with gstore_wannierize_and_write_gwan. See messages above.")
 end if

 ! TODO: Handle long-range part.
 if (dvdb%has_zeff .or. dvdb%has_quadrupoles) then
   ABI_WARNING("Treatment of long-range part not yet coded in gstore_wannierize_and_write_gwan!")
 end if

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is); my_nq = gqk%my_nq; my_nk = gqk%my_nk; my_npert = gqk%my_npert
   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq

   ! Initialize gkq%wan from ABIWAN.nc for this spin.
   keep_umats = .False.
   call gqk%wan%from_abiwan(dtfil%filabiwanin, spin, gstore%nsppol, keep_umats, dtfil%filnam_ds(4), gqk%comm%value)
   wan => gqk%wan

   ! Compute WS lattice vectors. Also allocate grpe_wwp with shape: (nr_p, nr_e, nwan, nwan, my_npert)
   call kptrlatt_from_ngkpt(gstore%ngqpt, qptrlatt_)
   call wan%setup_eph_ws_kq(gstore%cryst, gstore%ebands%shiftk(:,1), gstore%ebands%kptrlatt, qptrlatt_, &
                            gqk%my_pert_start, my_npert, gqk%pert_comm)

   nr_p = wan%nr_p; nr_e = wan%nr_e; nwan = wan%nwan
   !if (gqk%comm%me == master) call wan%print(units)

   ABI_MALLOC(emikr, (nr_e))
   ABI_MALLOC(emiqr, (nr_p))
   ABI_MALLOC(gww_pk, (nwan, nwan, my_npert, my_nk))
   ! Intermediate buffer to store the sum over k-points. Note my_nq.
   ABI_CALLOC(gww_epq, (nwan, nwan, nr_e, my_npert, my_nq))

   ! Loop over my q-points (partial sum over q)
   do my_iq=1,my_nq
     call gqk%myqpt(my_iq, gstore, weight_qq, qpt)

     !call get_kg(qpt, 1, ecut_lr, gstore%cryst%gmet, ng_q, gvec_q)
     !ABI_FREE(gvec_q)

     ! Loop over my k-points (partial sum over k)
     do my_ik=1,my_nk
       kpt = gqk%my_kpts(:,my_ik); kq = kpt + qpt
       ik = wan%krank%get_index(kpt); ikq = wan%krank%get_index(kq)
       ABI_CHECK(ik  /= -1, sjoin("Cannot find kpt: ", ktoa(kpt)))
       ABI_CHECK(ikq /= -1, sjoin("Cannot find k+q: ", ktoa(kq)))

       ! Get rotation matrices at k and k+q.
       nwin_k = wan%dimwin(ik); nwin_kq = wan%dimwin(ikq)

       ABI_MALLOC(g_bb, (nwin_kq, nwin_k))
       ABI_MALLOC(u_k, (1:nwin_k, 1:nwan))
       ABI_MALLOC(u_kq, (1:nwin_kq, 1:nwan))
       ABI_MALLOC(tmp_mat, (nwan, nwin_k))

       u_k = wan%u_k(1:nwin_k, 1:nwan, ik)
       u_kq = wan%u_k(1:nwin_kq, 1:nwan, ikq)

       do my_ip=1,gqk%my_npert
         !----------------------------------------------------------
         !  STEP 1: rotation to optimally smooth Bloch states
         !----------------------------------------------------------
         ! [Eqn. 24 of PRB 76, 165108 (2007)]
         ! g~(k,q) = U(k+q)^\dagger * g(k,q) * U(k)

         ! gqk%nb = wan%bmax - wan%bmin + 1
         ! (my_npert, nb, my_nq, nb, my_nk)
         ! (       p, b1_kq,     q, b2_k, k)  -->  <k+q, b1| D_{q,p}H |k, b2>

         ! Extract e-ph matrix elements from my_g buffer to align bands with the U matrices at k and k+q.
         jj = 0
         do ib_k=1,nb_k
           band_k = ib_k - wan%bmin + 1; if (.not. wan%lwindow(band_k, ik)) cycle
           jj = jj + 1; ii = 0
           do ib_kq=1,nb_kq
             band_kq = ib_kq - wan%bmin + 1; if (.not. wan%lwindow(band_kq, ikq)) cycle
             ii = ii + 1
             g_bb(ii, jj) = gqk%my_g(my_ip, ib_kq, my_iq, ib_k, my_ik)
           end do
         end do

         ! TODO: Remove LR part.
         !call handle_lr_term(cryst, qpt, ng_q, gvec_q, nwin_kq, nwin_k, nwan, u_kq, u_k, dvdb%zeff, dvdb%qstar, -1, g_bb)

         ! the two zgemm calls perform: epmats  = [ cu(ikq)^\dagger * epmatk ] * cu(ikk)
         ! [here we have a size-reduction from nbnd*nbnd to nwan*nwan]
         ! output stored in gww_pk(:,:, my_ip, my_ik)

         !gww_pk(:,:, my_ip, my_ik) = MATMUL(CONJG(TRANSPOSE(u_kq)), MATMUL(g_bb, u_k))
         call ZGEMM('N', 'N', nwan, nwin_k, nwin_kq, cone, g_bb, nwin_kq, u_k, nwin_k, czero, tmp_mat, nwan)
         call ZGEMM('C', 'N', nwin_kq, nwan, nwin_k, cone, u_kq, nwan, tmp_mat, nwan, czero, gww_pk(:,:, my_ip, my_ik), nwin_kq)

         !call ZGEMM('C', 'N', nwan, nbnd, nbnd, cone, u_kq, nbnd, epmatk(:, :, ik, imode), nbnd, czero, eptmp, nwan)
         !call ZGEMM('N', 'N', nwan, nwan, nbnd, cone, eptmp, nwan, u_k, nbnd, czero, epmats(:, :, ik, imode), nwan)
       end do ! my_ip

       ABI_FREE(g_bb)
       ABI_FREE(tmp_mat)

       !----------------------------------------------------------------------
       !  STEP 3: Fourier transform to obtain matrix elements in electron wannier basis
       !----------------------------------------------------------------------
       ! [Eqn. 24 of PRB 76, 165108 (2007)]
       ! g(R_e,q) = (1/nkc) sum_k e^{-ikR_e} g~(k,q)
       ! g(R_e,q) is epmatw (nwan,nwan,ir)

       do ir=1,nr_e
         emikr(ir) = exp(-j_dpc * two_pi * dot_product(kpt, wan%r_e(:, ir))) / dble(gstore%nkbz)
       end do
       ! gww_pk(nwan, nwan, my_npert,my_nk)
       do my_ip=1,gqk%my_npert
         do ir=1,nr_e
           gww_epq(:,:,ir, my_ip, my_iq) = gww_epq(:,:,ir, my_ip, my_iq) + emikr(ir) * gww_pk(:,:,my_ip, my_ik)
         end do
       end do

       ABI_FREE(u_k)
       ABI_FREE(u_kq)
     end do ! my_ik
   end do ! my_iq

   call xmpi_sum(gww_epq, gqk%kpt_comm%value, ierr)

   !----------------------------------------------------------
   !  Fourier transform to go into Wannier basis
   !----------------------------------------------------------
   ! [Eqn. 24 of PRB 76, 165108 (2007)]
   ! g(R_e,R_p) = (1/nq) sum_q e^{-iqR_p} g(R_e,q)

   ! Loop over my q-points (partial sum over q)
   do my_iq=1,my_nq
     call gqk%myqpt(my_iq, gstore, weight_qq, qpt)
     do ir=1,nr_p
       emiqr(ir) = exp(-j_dpc * two_pi * dot_product(qpt, wan%r_p(:, ir))) / dble(gstore%nqbz)
     end do

      do my_ip=1,my_npert
      do jwan=1,nwan
      do iwan=1,nwan
      do ir=1,nr_e
         wan%grpe_wwp(:, ir, iwan, jwan, my_ip) = wan%grpe_wwp(:, ir, iwan, jwan, my_ip) + &
           gww_epq(iwan, jwan, :, my_ip, my_iq) * emiqr(:)
      end do
      end do
      end do
      end do ! my_ip
   end do ! my_iq

   call xmpi_sum(wan%grpe_wwp, gqk%qpt_kpt_comm%value, ierr)

   if (gqk%comm%me == master) then
     write(std_out, '(a)') '#   R_e [Bohr]    max_{m,n,nu} |g(m,n,nu R_e,:)|  min_{m,n,nu} |g(m,n,nu R_e,:)|[Ha/Bohr] '
     do ir=1,nr_e
       write(std_out, *) wan%rmod_e(ir), maxval(abs(wan%grpe_wwp(:,ir,:,:,:))), sum(abs(wan%grpe_wwp(:,ir,:,:,:))) / size(wan%grpe_wwp(:,ir,:,:,:))
     end do
   end if

   ! Free memory for this spin.
   ABI_FREE(emikr)
   ABI_FREE(emiqr)
   ABI_FREE(gww_pk)
   ABI_FREE(gww_epq)
 end do ! my_is

 ! =====================
 ! Write data to GWAN.nc
 ! =====================
 do spin=1,gstore%nsppol
   my_is = gstore%spin2my_is(spin)
   if (my_is /= 0) then
     gqk => gstore%gqk(my_is)
     ! TODO: and I'm the in the first slice of gqk%comm ...
     !gqk%coords_qkpb_sumbp(ndims)
     call gqk%wan%ncwrite_gwan(dtfil, gstore%cryst, gstore%ebands, gqk%pert_comm)
   end if
   call xmpi_barrier(gstore%comm)
 end do ! spin

 call cwtime_report(" gstore_wannierize_and_write_gwan:", cpu, wall, gflops)

end subroutine gstore_wannierize_and_write_gwan
!!***

!!****f* m_gstore/handle_lr_term
!! NAME
!! handle_lr_term
!!
!! FUNCTION
!!  Add/Remove the long range term to/from the e-ph matrix elements.
!!
!! INPUTS
!!
!! SOURCE

!subroutine handle_lr_term(cryst, qpt, ng, gvec, nwin_kq, nwin_k, nwan, u_kq, u_k, zeff, qstar, isgn, g_bb)
!
!!Arguments ------------------------------------
! type(crystal_t),intent(in) :: cryst
! real(dp),intent(in) :: qpt(3)
! integer,intent(in) :: ng, nwin_kq, nwin_k, nwan, isgn, gvec(3,ng)
! complex(dp),intent(in) :: u_kq(nwin_kq, nwan), u_k(nwin_k,nwan)
! real(dp),intent(in) :: zeff(3,3,cryst%natom), qstar(3,3,3,cryst%natom)
! complex(dp),intent(inout) :: g_bb(nwin_kq, nwin_k)
!
!!Local variables-------------------------------
!!scalars
! !integer :: ig
! !character(len=500) :: msg
!!arrays
!! *************************************************************************
!
!end subroutine handle_lr_term
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_get_erange_mask
!! NAME
!!  gqk_get_erange_mask
!!
!! FUNCTION
!!  Compute MPI-distributed and global masks for electronic states allowed by erange.
!!
!! INPUTS
!!  gstore<gstore_t>=Electron-phonon object containing dimensions and related quantities.
!!  erange=Energy range:
!!    -- if both entries are negative, assume metal and include states within the
!! [efermi-abs(erange(1)), efermi+abs(erange(2))] window;
!!    -- otherwise, erange(1) & erange(2) select window wrt VBM & CBM, respectively.
!!
!! OUTPUT
!!  my_states(gqk%nb, gqk%my_nk)=Mask for selected states at this MPI proc.
!!  glob_states(gqk%nb, gqk%my_nk)=Global mask for selected states.
!!
!! SOURCE

subroutine gqk_get_erange_mask(gqk, gstore, erange, my_states, glob_states)

!Arguments ------------------------------------
!scalars
 class(gqk_t), target, intent(inout) :: gqk
 class(gstore_t), target, intent(in) :: gstore
!arrays
 real(dp), intent(in) :: erange(2)
 integer, intent(out) :: my_states(gqk%nb, gqk%my_nk), glob_states(gqk%nb, gqk%glob_nk)

!Local variables-------------------------------
!scalars
 class(ebands_t), pointer :: ebands
 type(gaps_t) :: gaps
 integer :: my_ik, ik_ibz, ik_glob, ib, bstart, gap_err, ierr
 real(dp) :: vmax, cmin, eig
 logical :: assume_gap
!----------------------------------------------------------------------

 ebands => gstore%ebands

 assume_gap = .not. all(erange < tol12)
 gaps = ebands%get_gaps(gap_err)

 if (assume_gap) then
   call gaps%print([std_out])
   vmax = gaps%vb_max(gqk%spin) + tol2 * eV_Ha
   cmin = gaps%cb_min(gqk%spin) - tol2 * eV_Ha
 else
   vmax = ebands%fermie
   cmin = ebands%fermie
 end if

 ! Fill the mask for allowed states
 my_states(:,:) = 0
 glob_states(:,:) = 0
 bstart = gstore%brange_spin(1, gqk%spin)

 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   ik_glob = my_ik + gqk%my_kstart - 1

   do ib=1,gqk%nb
     eig = ebands%eig(bstart + ib - 1, ik_ibz, gqk%spin)

     if (abs(erange(1)) > tol12) then
       ! Filter valence states.
       if (eig <= vmax .and. vmax - eig <= abs(erange(1))) then
         my_states(ib, my_ik) = 1
         glob_states(ib, ik_glob) = 1
       end if
     end if

     if (abs(erange(2)) > tol12) then
       ! Filter conduction states.
       if (eig >= cmin .and. eig - cmin <= abs(erange(2))) then
         my_states(ib, my_ik) = 1
         glob_states(ib, ik_glob) = 1
       end if
     end if

   enddo
 enddo

 call xmpi_sum(glob_states, gqk%kpt_comm%value, ierr)
 call gaps%free()

end subroutine gqk_get_erange_mask
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_filter_erange
!! NAME
!!  gqk_filter_erange
!!
!! FUNCTION
!!  Nullify all matrix elements connecting electronic states excluded by energy range.
!!
!! INPUTS
!!  gstore<gstore_t>=Electron-phonon object containing dimensions and related quantities.
!!  erange=Energy range:
!!    -- if both entries are negative, assume metal and include states within the
!! [efermi-abs(erange(1)), efermi+abs(erange(2))] window;
!!    -- otherwise, erange(1) & erange(2) select window wrt VBM & CBM, respectively.
!!
!! OUTPUT
!!
!! SOURCE

subroutine gqk_filter_erange(gqk, gstore, erange)

!Arguments ------------------------------------
!scalars
 class(gqk_t), target, intent(inout) :: gqk
 class(gstore_t), target, intent(in) :: gstore
!arrays
 real(dp), intent(in) :: erange(2)

!Local variables-------------------------------
!scalars
 integer :: my_ik, ik_glob, my_iq, ikq, ipert, ierr, ib, jb
 real(dp) :: wtq
 logical :: skip_nk, skip_mkq, skip_q
 type(krank_t) :: krank_kpts
 type(ebands_t), pointer :: ebands
!arrays
 integer :: my_states(gqk%nb, gqk%my_nk), glob_states(gqk%nb, gqk%glob_nk)
 real(dp) :: kpt(3), qpt(3), kpq(3), kpts(3, gqk%glob_nk), my_qpts(3, gqk%my_nq)
!----------------------------------------------------------------------

 ebands => gstore%ebands

 ! Compute masks
 call gqk%get_erange_mask(gstore, erange, my_states, glob_states)

 ! Get global krank for k+q transitions
 kpts(:, :) = zero
 do my_ik=1,gqk%my_nk
   ik_glob = my_ik + gqk%my_kstart - 1
   kpts(:, ik_glob) = gqk%my_kpts(:, my_ik)
 enddo
 call xmpi_sum(kpts, gqk%kpt_comm%value, ierr)

 call krank_kpts%from_kptrlatt(gqk%glob_nk, kpts, ebands%kptrlatt, compute_invrank=.True.)

 ! Get all q-points for this proc
 do my_iq=1,gqk%my_nq
   call gqk%myqpt(my_iq, gstore, wtq, my_qpts(:, my_iq))
 enddo

 ! Nullify matrix elements connecting the excluded states
 do my_ik=1,gqk%my_nk
   kpt(:) = gqk%my_kpts(:, my_ik)

   do ib=1,gqk%nb
     ! |nk> is forbidden
     skip_nk = .false.
     if (my_states(ib, my_ik) == 0) skip_nk = .true.

     do my_iq=1,gqk%my_nq
       qpt(:) = my_qpts(:, my_iq)

       ! Find k+q --> k' index in krank_kpts
       kpq(:) = kpt(:) + qpt(:)
       ikq = krank_kpts%get_index(kpq)

       ! k+q falls outside the filtered kpts pool
       skip_q = .false.
       if (ikq == -1) skip_q = .true.

       do jb=1,gqk%nb
         ! |mk+q> is forbidden
         skip_mkq = .false.
         if (glob_states(jb, ikq) == 0) skip_mkq = .true.

         do ipert=1,gqk%my_npert

           if (skip_nk .or. skip_q .or. skip_mkq) then
             select case (gqk%cplex)
             case (1)
               gqk%my_g2(ipert, jb, my_iq, ib, my_ik) = zero
             case (2)
               gqk%my_g(ipert, jb, my_iq, ib, my_ik) = czero
             case default
               ABI_ERROR(sjoin("Invalid gqk%cplex:", itoa(gqk%cplex)))
             end select
           end if

         enddo
       enddo
     enddo
   enddo
 enddo

 call krank_kpts%free()

end subroutine gqk_filter_erange
!!***

end module m_gstore
!!***

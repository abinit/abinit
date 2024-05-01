!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!
!!  This module implements the gstore_t object that allows one to **precompute"" the e-ph matrix elements g
!!  and store them in memory with a MPI-distributed data structure.
!!
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
!!          g(q, k) = <k + q, m, spin| \Delta_{q, \nu} V^{spin}_{scf} | k, n, spin>
!!
!!  There are lots of technical details that should be discussed but, roughly speaking,
!!  the gqk API allows one to:
!!
!!   - select whether q or k should be in the IBZ or in the BZ.
!!     NB: It is not possible to use the IBZ both for q and k as g(Sk, q) = g(k, S^{-1}q)
!!     thus one has to select the appropriate zones beforehand.
!!
!!   - filter bands and/or k/q wavevectors according to some criterion.
!!     In superconductors, for instance, only k/k+q states on the Fermi surface are needed.
!!     In semiconductors, one can include only k, k+q inside an energy window around the band edge.
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
!!         and implements filtering techniques. The MPI grid is automatically generated at this level
!!
!!      2) gstore_compute evaluates the e-ph matrix elements in parallel and dumps the results to GSTORE.nc
!!
!!      3) gstore%from_ncpath reconstructs the object from GSTORE.nc
!!
!!  In a typical scenario, one uses eph_task 11 to generate GSTORE.nc i.e. steps 1) and 2).
!!  Then one introduces a new value of eph_task in which we read the object from file and call
!!  a specialized routine that implements the "post-processing" steps needed to compute the physical properties of interest.
!!
!!  Now, let's discuss the MPI-distribution.
!!
!!  The (q, k) matrix is distributed inside a 2D cartesian grid using block distribution.
!!  This is schematic representation for 4 procs with 2 procs for k and 2 procs for q:
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
!!  Each processor stores all the (band_1, band_2) transitions for a given (q, k) pair.
!!
!!  Perturbations can be optionally distributed along a third axis (pert_comm).
!!  Note, however, that the parallelism over perturbations is not expected to be the most efficient
!!  although it allows one to reduce the memory required to store the scattering potential in the supercell
!!  as we can distribute W(r, R, 3 * natom) over the last dimension.
!!
!!  For electronic properties, one usually uses k-points in the IBZ and q-points in the BZ.
!!  hence the parallelism over q-points is the most efficient one in terms of wall-time.
!!  Keep in mind, however, that the k-point parallelism allows one to reduce the memory allocated for the
!!  wavefunctions. Using some procs for k-point is also beneficial in terms of performance
!!  as we reduce load imbalance with the number of procs in qpt_comm does not divide nqbz.
!!
!!  NB: If nsppol == 2, we create two gqk objects, one for each spin.
!!  The reason is that dimensions such as the number of effective bands/q-points/k-points
!!  depends on spin if filters are employed.
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2024 ABINIT group (MG)
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

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_htetra
 use libtetrabz
 use m_ebands
 use, intrinsic :: iso_c_binding
 use netcdf
 use m_nctk
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_crystal
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_dtset
 use m_dtfil
 use m_wfd
 use m_ephtk
 use m_mkffnl

 use defs_abitypes,    only : mpi_type
 use m_time,           only : cwtime, cwtime_report, sec2str
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat, replace_ch0
 !use m_yaml,           only : yamldoc_t
 use m_numeric_tools,  only : arth, get_diag, isdiagmat
 use m_krank,          only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_io_tools,       only : iomode_from_fname
 use m_special_funcs,  only : gaussian
 use m_copy,           only : alloc_copy
 use m_fftcore,        only : ngfft_seq, get_kg
 use m_cgtools,        only : cg_zdotc
 use m_kg,             only : getph, mkkpg
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_hdr,            only : hdr_type, fform_from_ext, hdr_ncread
 use m_symtk,          only : matr3inv
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_sort, kpts_pack_in_stars
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_ifc,            only : ifc_type
 use m_phonons,        only : pheigvec_rotate
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_gstore/gqk_t
!! NAME
!! gqk_t
!!
!! FUNCTION
!!  This object stores MPI-distributed e-ph matrix elements for
!!  a given spin index (if collinear magnetism i.e. nsppol 2).
!!  local dimensions and arrays start with `my_`,
!!  global dimensions start with `glob_
!!
!! SOURCE

type, public :: gqk_t

  integer :: cplex = -1
  ! 1 if |g|^2 should be stored
  ! 2 if complex valued g (mind the gauge)

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
  ! The last band is bstop (global indices)

  integer :: my_npert = -1
  ! Number of perturbations treated by this MPI rank.

  integer :: glob_nk = -1, glob_nq = -1
  ! Total number of k/q points in global matrix.
  ! Note that k-points/q-points can be filtered.
  ! Use kzone, qzone and kfilter to interpret these dimensions.

  integer :: my_nk = -1, my_nq = -1
  ! Number of k/q points treated by this MPI proc.
  ! Used to loop and allocate local arrays.

  integer :: my_kstart = -1, my_qstart = -1
  ! Index of the first k/q point in the global matrix treated by this MPI proc

  integer,allocatable :: my_k2ibz(:,:)
  ! (6, my_nk)
  ! Mapping my_kpoints --> kibz

  !integer,allocatable :: my_k2bz(:,:)
  ! (my_nk)
  ! Mapping my_kpoints --> ik_bz

  real(dp),allocatable :: my_kpts(:,:)
  ! (3, my_nkpt)
  ! k-points treated by this MPI proc.

  real(dp),allocatable :: my_wtk(:)
  ! (my_nkpt)
  ! k-point weights

  integer,allocatable :: my_q2ibz(:,:)
  ! (6, my_nq)
  ! Mapping my_qpoints --> qibz

  integer,allocatable :: my_q2bz(:)
  ! (my_nq)
  ! Mapping my_iq index --> iq_bz index in the full BZ

  !integer,allocatable :: my_q2glob(:)
  ! (my_nq)
  ! Mapping my_iq index --> global index in the g(q, k) matrix.

  real(dp),allocatable :: vk_cart_ibz(:,:,:)
  ! (3, nb, nkibz)
  ! Diagonal v_{m, m,k} for k in the IBZ.
  ! Values in the BZ can be reconstructed by symmetry.
  ! Allocated if gstore%with_vk == 1

  real(dp),allocatable :: vkmat_cart_ibz(:,:,:,:,:)
  ! (3, nb, nb, nkibz)
  ! v_{m, n,k} for the k in the IBZ
  ! Allocated if gstore%with_vk in (1, 2)

  integer,allocatable :: my_iperts(:)
  ! (my_npert)
  ! List of perturbation indices treated by this MPI proc.
  ! Contiguous indices.

  complex(dp), allocatable :: my_g(:,:,:,:,:)
  ! (my_npert, nb, my_nq, nb, my_nk)
  ! (       p, b1,     q, b2, k)  -->  <k+q, b1| D_{q,p}H |k, b2>
  ! e-ph matrix elements g (local buffer). Allocated if cplex == 2

  real(dp), allocatable :: my_g2(:,:,:,:,:)
  ! (my_npert, nb, my_nq, nb, my_nk)
  ! |g|^2 (local buffer). Allocated if cplex == 1

  integer :: coords_qkp(3)
  ! Coordinates of this processor in the (q, k, pert) Cartesian grid.

  type(xcomm_t) :: kpt_comm
   ! MPI communicator over k-points

  type(xcomm_t) :: qpt_comm
   ! MPI communicator over q-points

  type(xcomm_t) :: pert_comm
   ! MPI communicator over atomic perturbations.

  type(xcomm_t) :: qpt_pert_comm
   ! MPI communicator over the 2d grid (qpt, atomic perturbations)

  type(xcomm_t) :: grid_comm
   ! MPI communicator for full (q,k,pert) grid.

  real(dp),allocatable :: my_wnuq(:,:)
  ! (my_npert, my_nq)
  ! Phonon frequencies in Ha (MPI distributed)

  real(dp),allocatable :: my_displ_cart(:,:,:,:,:)
  ! (2, 3, cryst%natom, my_npert, my_nq))
  ! Phonon displacements (MPI distributed)

 contains

  procedure :: myqpt => gqk_myqpt
  ! Return the q-point and the weight from my local index my_iq

  !procedure :: my_qweight => gqk_my_qweight
  ! Return the q-point weight from my local index my_iq

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
!!    - pointers to the cyrstalline structure, the KS bands, the IFCs.
!!    - arrays that do not depend on the spin such as the IBZ and weights for k/q-points.
!!    - metadata such as kzone, qzone and kfilter that are needed to interpret
!!      the storage mode used for the g(k, q)
!!
!! NB: the e-ph matrix element are stored in gstore%qqk(my_is)
!!     where my_is counts the number of spins treated by this MPI processor.
!!
!! NOTES
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
  !  qptopt=option for the generation of q points (defines whether spatial symmetries and/or time-reversal can be used)

  character(len=fnlen) :: path = " "
  ! Path to the nc file associated to the gstore

  character(len=fnlen) :: wfk0_path = " "

  character(len=fnlen) :: kzone = " ", qzone = " "
   ! Specifies whether k- or q-points are in the BZ or in the IBZ.
   ! Possible values are "ibz" or "bz".
   ! Note that the combination ("ibz", "ibz") is not allowed

  character(len=fnlen) :: kfilter = "none"
  ! Specifies the tecnique used to filter k-points.
  ! Possible values: "none", "fs_tetra", "erange".

  real(dp),allocatable :: erange_spin(:, :)
  ! (2, nsppol)
  ! Energy window. zero if not used. Requires kfilter == "erange"

  type(crystal_t), pointer :: cryst => null()

  type(ebands_t), pointer :: ebands => null()

  type(ifc_type), pointer :: ifc => null()

  type(krank_t) :: krank_ibz !, qrank
  ! Object used to find k-points in the IBZ and map BZ to IBZ.

  integer,allocatable :: my_spins(:)
   ! (%my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2.

  integer,allocatable :: brange_spin(:, :)
  ! (2, nsppol)

  !integer :: max_nb = -1
  ! Max number of bands over spin

  integer,allocatable :: glob_nk_spin(:), glob_nq_spin(:)
  ! (nsppol)

  real(dp), contiguous, pointer :: kibz(:,:)
  ! k-points in the IBZ. Points to ebands%kptns
  ! (3, nkibz)

  real(dp),allocatable :: delta_ef_kibz_spin(:,:,:)
  ! (nb, gstore%nkibz, nsppol))
  ! Tetrahedron weights at eF in the IBZ.

  real(dp), allocatable :: qibz(:,:)
  ! (3, nqibz)
  ! q-points in the IBZ

  real(dp), allocatable :: wtq(:)
  ! (nqibz)
  ! q-points weights in the IBZ

  !integer :: qptrlatt(3, 3) = -1  ! kptrlatt(3, 3) = -1,
   ! k-mesh and q-mesh

  !real(dp),allocatable :: kshift(:, :), qshift(:, :)
  ! k/q-mesh shift (well, q-mesh is usually gamma-centered)

  type(gqk_t), allocatable :: gqk(:)
  ! (my_nspins)
  ! Datastructure storing e-ph matrix elements.

  !type(htetra_t) :: ktetra
  ! Used to evaluate integrals in k-space with the tetrahedron method.

  !type(htetra_t) :: qtetra
  ! Used to evaluate integrals in q-space with the tetrahedron method.

contains

  procedure :: fill_bks_mask => gstore_fill_bks_mask
  ! Fill the table used to read (b, k, s) wavefunctions from the WFK file
  ! keeping into account the distribution of the e-ph matrix elements.

  procedure :: get_mpw_gmax => gstore_get_mpw_gmax

  procedure :: spin2my_is => gstore_spin2my_is
  !  Return the local spin index from the global spin index.
  !  0 if this spin is not treated by this MPI proc.

  procedure :: free => gstore_free
  ! Free memory

  procedure :: print => gstore_print
  ! Print info on the object

  procedure, private :: distribute_spins__ => gstore_distribute_spins
  ! Distribute spins (nsppol = 2) and create indirect mapping to spin index.

  procedure, private :: set_mpi_grid__ => gstore_set_mpi_grid__
  ! Set the MPI cartesian grid

  procedure, private :: malloc__ => gstore_malloc__
  ! Allocate local buffers once the MPI grid has been initialized.

  procedure, private :: filter_fs_tetra__ => gstore_filter_fs_tetra__
  ! Select k-points on the FS using the tetrahedron method

  procedure, private :: filter_erange__ => gstore_filter_erange__
  ! Select k-points inside an energy window

  procedure :: compute => gstore_compute
  ! Compute e-ph matrix elements.

  procedure :: calc_my_phonons => gstore_calc_my_phonons
  ! Helper function to compute ph quantities for all q-points treated by the MPI proc.

  procedure :: get_lambda_iso_iw => gstore_get_lambda_iso_iw
  ! Compute isotropic lamdda(iw) along the imaginary axis.

  procedure :: get_a2fw => gstore_get_a2fw
  ! Compute Eliashberg function a^2F(w).

  procedure :: from_ncpath => gstore_from_ncpath
  ! Reconstruct object from netcdf file

  procedure :: init => gstore_init
  ! Build object

end type gstore_t
!!***

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
!! OUTPUT
!!
!! SOURCE

subroutine gstore_init(gstore, path, dtset, wfk0_hdr, cryst, ebands, ifc, comm)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(out) :: gstore
 integer,intent(in) :: comm
 character(len=*),intent(in) :: path
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(in) :: wfk0_hdr
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),target,intent(in) :: ebands
 class(ifc_type),target,intent(in) :: ifc

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: all_nproc, my_rank, ierr, my_nshiftq, nsppol, spin, natom3, cnt, qtimrev
 integer :: ik_ibz, ik_bz, iq_bz, iq_ibz, ebands_timrev, max_nq, max_nk, ncid, spin_ncid, ncerr, gstore_fform
 real(dp) :: cpu, wall, gflops
 character(len=10) :: priority
 !character(len=5000) :: msg
 type(krank_t) :: qrank
!arrays
 integer :: ngqpt(3), qptrlatt(3,3), comm_spin(ebands%nsppol), nproc_spin(ebands%nsppol)
 !integer :: glob_nk_spin(ebands%nsppol), glob_nq_spin(ebands%nsppol)
 integer,allocatable :: qbz2ibz(:,:), kbz2ibz(:,:), kibz2bz(:), qibz2bz(:), qglob2bz(:,:), kglob2bz(:,:)
 integer,allocatable :: select_qbz_spin(:,:), select_kbz_spin(:,:) !, done_qbz_spin(:,:)
 real(dp):: my_shiftq(3,1)
 real(dp),allocatable :: qbz(:,:), wtk(:), kibz(:,:), kbz(:,:)
 !integer :: out_kptrlatt(3,3)
 !real(dp),allocatable :: out_kibz(:,:), out_wtk(:)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Set basic parameters.
 gstore%comm = comm
 gstore%nsppol = nsppol
 gstore%path = path

 ! Get references to other data structures.
 gstore%cryst => cryst
 gstore%ebands => ebands
 gstore%ifc => ifc
 gstore%kibz => ebands%kptns

 ! Set metadata.
 gstore%kzone =  dtset%gstore_kzone; gstore%qzone = dtset%gstore_qzone; gstore%kfilter = dtset%gstore_kfilter
 gstore%with_vk = dtset%gstore_with_vk

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
 !gstore%kptrlatt(3, 3)
 !gstore%kshift(3, 1)
 !gstore%qptrlatt(3, 3)
 !gstore%qshift(3, 1)

 ! Distribute spins, create indirect mapping to spin index and init gstore%brange_spin
 ABI_CHECK_ILEQ(dtset%mband, ebands%mband, "dtset%mband > ebands%mband")
 call gstore%distribute_spins__(dtset%mband, dtset%gstore_brange, nproc_spin, comm_spin, comm)

 ! Define q-mesh: either from DVDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation)
 ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if

 ! TODO: Should fix bz2ibz to use the same conventions as krank and listkk
 ! NB: only sigmaph seems to be using this optional argument

 ! Setup qIBZ, weights and BZ.
 ! Assume qptopt == kptopt unless value is specified in input
 qptrlatt = 0; qptrlatt(1, 1) = ngqpt(1); qptrlatt(2, 2) = ngqpt(2); qptrlatt(3, 3) = ngqpt(3)
 gstore%qptopt = ebands%kptopt; if (dtset%qptopt /= 0) gstore%qptopt = dtset%qptopt
 qtimrev = kpts_timrev_from_kptopt(gstore%qptopt)

 !call wrtout(std_out, sjoin(" Generating q-IBZ for with qptopt:", itoa(gstore%qptopt)))
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, gstore%qptopt, my_nshiftq, my_shiftq, &
                             gstore%nqibz, gstore%qibz, gstore%wtq, gstore%nqbz, qbz)
                             !new_kptrlatt=gstore%qptrlatt, new_shiftk=gstore%qshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(qbz2ibz, (6, gstore%nqbz))

 qrank = krank_from_kptrlatt(gstore%nqibz, gstore%qibz, qptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", qtimrev, cryst, qrank, gstore%nqbz, qbz, qbz2ibz) /= 0) then
   ABI_ERROR("Cannot map qBZ to IBZ!")
 end if
 call qrank%free()

 ! Order qbz by stars and rearrange entries in qbz2ibz table.
 call kpts_pack_in_stars(gstore%nqbz, qbz, qbz2ibz)

 !call kpts_print_kmap(std_out, qibz, qbz, qbz2ibz)
 !do iq_bz=1,gstore%nqbz
 !  print *, "iq_bz -> iq_ibz", qbz2ibz(1, iq_bz), qbz(:, iq_bz)
 !end do

 call get_ibz2bz(gstore%nqibz, gstore%nqbz, qbz2ibz, qibz2bz, ierr)
 ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for q-points!")

 ! Get full BZ associated to ebands
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
                             gstore%nkibz, kibz, wtk, gstore%nkbz, kbz) !, bz2ibz=bz2ibz)
                             !new_kptrlatt=gstore%kptrlatt, new_shiftk=gstore%kshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ! In principle kibz should be equal to ebands%kptns
 ABI_CHECK(gstore%nkibz == ebands%nkpt, "nkibz != ebands%nkpt")
 ABI_CHECK(all(abs(gstore%kibz - kibz) < tol12), "ebands%kibz != kibz")
 ABI_FREE(kibz)

 ! Note symrel and use_symrec=.False. in get_mapping.
 ! This means that this table can be used to symmetrize wavefunctions in cgtk_rotate.
 ! TODO This ambiguity should be removed. Change cgtk_rotate so that we can use the symrec convention.

 ABI_MALLOC(kbz2ibz, (6, gstore%nkbz))
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 gstore%krank_ibz = krank_from_kptrlatt(gstore%nkibz, gstore%kibz, ebands%kptrlatt, compute_invrank=.False.)
 if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, gstore%nkbz, kbz, kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 call get_ibz2bz(gstore%nkibz, gstore%nkbz, kbz2ibz, kibz2bz, ierr)
 ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for k-points")

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

 ! Here we filter the electronic wavevectors k and recompute select_qbz_spin and select_kbz_spin.

 select case (gstore%kfilter)
 case ("none")
   continue

 case ("erange")
   call gstore%filter_erange__(qbz, qbz2ibz, qibz2bz, kbz, gstore%kibz, kbz2ibz, kibz2bz, &
                               select_qbz_spin, select_kbz_spin)

 case ("fs_tetra")
   ! Use the tetrahedron method to filter k- and k+q points on the FS in metals
   ! and define gstore%brange_spin automatically.
   call gstore%filter_fs_tetra__(qbz, qbz2ibz, qibz2bz, kbz, gstore%kibz, kbz2ibz, kibz2bz, &
                                 select_qbz_spin, select_kbz_spin)
 case default
   ABI_ERROR(sjoin("Invalid kfilter:", gstore%kfilter))
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
 ABI_ICALLOC(kglob2bz, (max_nk, nsppol))

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
       cnt = cnt + 1; kglob2bz(cnt, spin) = ik_bz
     end if
   end do
 end do

 ! =============================================
 ! Initialize gqk basic dimensions and MPI grid
 ! =============================================
 !priority = "qk"
 call priority_from_eph_task(dtset%eph_task, priority)

 call gstore%set_mpi_grid__(dtset%gstore_cplex, dtset%eph_np_pqbks, priority, nproc_spin, comm_spin)
 call xmpi_comm_free(comm_spin)

 ! At this point, we have the Cartesian grid (one per spin if any)
 ! and we can finally allocate and distribute other arrays.
 call gstore%malloc__(0, max_nq, qglob2bz, max_nk, kglob2bz, qbz2ibz, kbz2ibz)

 ! Initialize GSTORE.nc file i.e. define dimensions and arrays
 ! Entries such as the e-ph matrix elements will be filled afterwards in gstore_compute.
 ! Master node writes basic objects and gstore dimensions

 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, gstore%path, xmpi_comm_self))

   ! Write the abinit header with metadata, structure and occupancies.
   gstore_fform = fform_from_ext("GSTORE.nc")
   NCF_CHECK(wfk0_hdr%ncwrite(ncid, gstore_fform, spinat=dtset%spinat, nc_define=.True.))
   !NCF_CHECK(gstore%cryst%ncwrite(ncid))

   ! Add eigenvalues and occupations.
   NCF_CHECK(ebands_ncwrite(gstore%ebands, ncid))

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

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "gstore_with_vk", "gstore_qptopt"])
   NCF_CHECK(ncerr)
   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
   !NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("gstore_qibz", "dp", "three, gstore_nqibz"), &
     nctkarr_t("gstore_qbz", "dp", "three, gstore_nqbz"), &
     nctkarr_t("gstore_wtq", "dp", "gstore_nqibz"), &
     nctkarr_t("gstore_kbz", "dp", "three, gstore_nkbz"), &
     nctkarr_t("gstore_kzone", "c", "fnlen"), &
     nctkarr_t("gstore_qzone", "c", "fnlen"), &
     nctkarr_t("gstore_kfilter", "c", "fnlen"), &
     nctkarr_t("gstore_wfk0_path", "c", "fnlen"), &
     nctkarr_t("gstore_brange_spin", "i", "two, number_of_spins"), &
     nctkarr_t("gstore_erange_spin", "dp", "two, number_of_spins"), &
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

   ! internal table used to restart computation. Init with zeros.
   NCF_CHECK(nf90_def_var_fill(ncid, vid("gstore_done_qbz_spin"), NF90_FILL, zero))

   ! Optional arrays
   if (allocated(gstore%delta_ef_kibz_spin)) then
     ncerr = nctk_def_arrays(ncid, &
       nctkarr_t("gstore_delta_ef_kibz_spin", "dp", "gstore_max_nb, gstore_nkibz, number_of_spins"))
   end if
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_with_vk"), gstore%with_vk))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qptopt"), gstore%qptopt))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kzone"), trim(gstore%kzone)))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qzone"), trim(gstore%qzone)))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kfilter"), trim(gstore%kfilter)))

   ! Write arrays
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qibz"), gstore%qibz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qbz"), qbz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_wtq"), gstore%wtq))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kbz"), kbz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_brange_spin"), gstore%brange_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_erange_spin"), gstore%erange_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_glob_nq_spin"), gstore%glob_nq_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_glob_nk_spin"), gstore%glob_nk_spin))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kbz2ibz"), kbz2ibz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qbz2ibz"), qbz2ibz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_qglob2bz"), qglob2bz))
   NCF_CHECK(nf90_put_var(ncid, vid("gstore_kglob2bz"), kglob2bz))
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
        nctkdim_t("glob_nk", gstore%glob_nk_spin(spin)), &
        nctkdim_t("glob_nq", gstore%glob_nq_spin(spin))  &
     ], defmode=.True.)
     NCF_CHECK(ncerr)

     ! Define scalars
     ncerr = nctk_def_iscalars(spin_ncid, [character(len=nctk_slen) :: "bstart"])
     NCF_CHECK(ncerr)

     ! arrays in gqk_spin group with the precious stuff. Note glob dimensions
     ncerr = nctk_def_arrays(spin_ncid, [ &
       nctkarr_t("gvals", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq") &
     ])
     NCF_CHECK(ncerr)

     ! IMPORTANT: Init with zeros.
     NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("gvals"), NF90_FILL, zero))

     select case(gstore%with_vk)
     case (1)
       ! Diagonal terms only
       NCF_CHECK(nctk_def_arrays(spin_ncid, nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz")))
       NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("vk_cart_ibz"), NF90_FILL, zero))
     case (2)
       ! Full (nb x nb) matrix.
       NCF_CHECK(nctk_def_arrays(spin_ncid, nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz")))
       NCF_CHECK(nf90_def_var_fill(spin_ncid, vid_spin("vkmat_cart_ibz"), NF90_FILL, zero))
     end select

     ! Write (small) data
     NCF_CHECK(nctk_set_datamode(spin_ncid))
     NCF_CHECK(nf90_put_var(spin_ncid, vid_spin("bstart"), gstore%brange_spin(1, spin)))
   end do ! spin

   NCF_CHECK(nf90_close(ncid))
 end if ! master

 call xmpi_barrier(gstore%comm)

 ABI_FREE(wtk)
 ABI_FREE(kbz)
 ABI_FREE(qbz)
 ABI_FREE(qibz2bz)
 ABI_FREE(kibz2bz)
 ABI_FREE(select_qbz_spin)
 ABI_FREE(select_kbz_spin)
 ABI_FREE(qglob2bz)
 ABI_FREE(kglob2bz)
 ABI_FREE(qbz2ibz)
 ABI_FREE(kbz2ibz)

 call cwtime_report(" gstore_init:", cpu, wall, gflops)

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(NCID, vname)
 end function vid
 integer function vid_spin(vname)
   character(len=*),intent(in) :: vname
   vid_spin = nctk_idname(SPIN_NCID, vname)
 end function vid_spin

end subroutine gstore_init
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/priority_from_eph_task
!! NAME
!! priority_from_eph_task
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine priority_from_eph_task(eph_task, priority)

!Arguments ------------------------------------
 integer,intent(in) :: eph_task
 character(len=*),intent(out) :: priority
!----------------------------------------------------------------------

 select case (eph_task)
 case (11)
   priority = "qk"
 case (12, -12)
   priority = "q"
 case (14)
   priority = "kq"
 case default
   ABI_ERROR(sjoin("Please register default priority for eph_task:", itoa(eph_task)))
 end select

end subroutine priority_from_eph_task
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_distribute_spins
!! NAME
!! gstore_distribute_spins
!!
!! FUNCTION
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
 integer :: spin, my_rank, ierr, color, nsppol, all_nproc
!arrays
 integer :: buff_spin(2)
!----------------------------------------------------------------------

 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 nsppol = gstore%nsppol

 gstore%my_nspins = 0
 ABI_MALLOC(gstore%brange_spin, (2, nsppol))

 do spin=1,nsppol
   ! NB: If MPI_UNDEFINED is passed as the colour value, the subgroup in which the calling
   ! MPI process will be placed is MPI_COMM_NULL
   color = 1
   if (nsppol == 2 .and. all_nproc > 1) then
     color = xmpi_undefined
     if (spin == 1 .and. my_rank <= (all_nproc - 1) / 2) color = 1
     if (spin == 2 .and. my_rank > (all_nproc - 1) / 2) color = 1
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

subroutine gstore_set_mpi_grid__(gstore, gstore_cplex, eph_np_pqbks, priority, nproc_spin, comm_spin)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: gstore_cplex
 integer,intent(in) :: eph_np_pqbks(5), nproc_spin(gstore%nsppol), comm_spin(gstore%nsppol)
 character(len=*),intent(in) :: priority

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, ndims = 3
 integer :: spin, my_is, np, my_rank, ierr !, ii
 type(gqk_t),pointer :: gqk
 character(len=5000) :: msg
 character(len=10) :: order
 integer :: comm_cart, me_cart, dims(ndims)
 logical :: reorder, periods(ndims), keepdim(ndims)
!----------------------------------------------------------------------

 my_rank = xmpi_comm_rank(gstore%comm)

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   gqk%spin = spin; gqk%natom3 = 3 * gstore%cryst%natom; gqk%cplex = gstore_cplex
   ABI_CHECK_IRANGE(gqk%cplex, 1, 2, "gstore_cplex")

   ! Compute bstart and band size for this spin.
   gqk%bstart = gstore%brange_spin(1, spin)
   gqk%bstop = gstore%brange_spin(2, spin)
   gqk%nb = gstore%brange_spin(2, spin) - gstore%brange_spin(1, spin) + 1

   ! Store global shape of the q/k matrix for this spin.
   gqk%glob_nq = gstore%glob_nq_spin(spin)
   gqk%glob_nk = gstore%glob_nk_spin(spin)
 end do

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   ! Init for sequential execution.
   gqk%qpt_comm%nproc = 1; gqk%kpt_comm%nproc = 1; gqk%pert_comm%nproc = 1
   gqk%my_npert = gqk%natom3

   if (any(eph_np_pqbks /= 0)) then
     ! Use parameters from input file. Need to perform sanity check though.
     gqk%pert_comm%nproc = eph_np_pqbks(1)
     gqk%qpt_comm%nproc  = eph_np_pqbks(2)
     ABI_CHECK(eph_np_pqbks(3) == 1, "Band parallelism not implemented in gstore")
     gqk%kpt_comm%nproc = eph_np_pqbks(4)
     !gqk%spin_comm%nproc = eph_np_pqbks(5)
     gqk%my_npert = gqk%natom3 / gqk%pert_comm%nproc
     ABI_CHECK(gqk%my_npert > 0, "pert_comm_nproc cannot be greater than 3 * natom.")
     ABI_CHECK(mod(gqk%natom3, gqk%pert_comm%nproc) == 0, "pert_comm_nproc must divide 3 * natom.")

   else
     ! Automatic grid generation (hopefully smart)
     ! Keep in mind that in gstore_build, the first loop is over q-points
     ! in order to reduce the number of interpolations of the DFPT potentials in q-space
     ! hence the q-point parallelism is expected to be more efficient.
     ! On the other hand, the k-point parallelism and the perturbation parallelism
     ! allow one to reduce the memory requirements associated to the wavefunctions (kpt) and
     ! the scattering potentials in the supercell (perturbations).
     ! Here we try to optimize performance but it's clear that for large systems the user
     ! should specify eph_np_pqbks in the input file.

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

     np = nproc_spin(spin)
     call xmpi_distrib_2d(np, order, gqk%glob_nq, gqk%glob_nk, gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, ierr)
     ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np), " with priority: ", priority))
     !if (ierr /= 0) call xmpi_distrib_2d_extra(np, gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, gqk%pert_comm%nproc, ierr)
   end if

   ! Consistency check.
   if (gqk%pert_comm%nproc * gqk%qpt_comm%nproc * gqk%kpt_comm%nproc /= nproc_spin(spin)) then
     write(msg, "(a,i0,3a, 4(a,1x,i0))") &
       "Cannot create 3d Cartesian grid with total nproc: ", nproc_spin(spin), ch10, &
       "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
       "qpt_nproc (", gqk%qpt_comm%nproc, ") x kpt_nproc (", gqk%kpt_comm%nproc, ")  x kpt_nproc", gqk%pert_comm%nproc, &
       ") != ", nproc_spin(spin)
     ABI_ERROR(msg)
   end if

 end do ! my_is

 ! For each spin treated by this rank, create 3D cartesian communicator: (q-points, k-points, perturbations)
 periods(:) = .False.; reorder = .False.

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   dims = [gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, gqk%pert_comm%nproc]

   ! Note comm_spin(spin)
   gqk%grid_comm = xcomm_from_mpi_int(comm_spin(spin))

#ifdef HAVE_MPI
   call MPI_CART_CREATE(gqk%grid_comm, ndims, dims, periods, reorder, comm_cart, ierr)

   ! Find the index and coordinates of the current processor
   call MPI_COMM_RANK(comm_cart, me_cart, ierr)
   call MPI_CART_COORDS(comm_cart, me_cart, ndims, gqk%coords_qkp, ierr)

   ! Create communicator for q-points
   keepdim = .False.; keepdim(1) = .True.; call gqk%qpt_comm%from_cart_sub(comm_cart, keepdim)
   ! Create communicator for k-points
   keepdim = .False.; keepdim(2) = .True.; call gqk%kpt_comm%from_cart_sub(comm_cart, keepdim)
   ! Create communicator for perturbations.
   keepdim = .False.; keepdim(3) = .True.; call gqk%pert_comm%from_cart_sub(comm_cart, keepdim)
   ! Create communicator for the (qpt, pert) 2D grid
   keepdim = .False.; keepdim(1) = .True.; keepdim(3) = .True.; call gqk%qpt_pert_comm%from_cart_sub(comm_cart, keepdim)
   call xmpi_comm_free(comm_cart)
#endif

   ! Distribute perturbations inside per_comm using block distribution.
   call xmpi_split_block(gqk%natom3, gqk%pert_comm%value, gqk%my_npert, gqk%my_iperts)
 end do ! my_is

 if (my_rank == master) then
   call gstore%print(std_out)
   call gstore%print(ab_out)
 end if

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
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_print(gstore, unit, header, prtvol)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: unit
 character(len=*),optional,intent(in) :: header
 integer,optional,intent(in) :: prtvol

!Local variables ------------------------------
 integer :: my_is, spin, my_prtvol
 type(gqk_t),pointer :: gqk
 !character(len=500) :: msg
 !type(yamldoc_t) :: ydoc
!----------------------------------------------------------------------

 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 if (present(header)) then
   call wrtout(unit, header)
 else
   call wrtout(unit, " === Gstore parameters ===")
 end if
 !call ebands_print(gstore%ebands, header="Electron bands", unit=unit, prtvol=my_prtvol)
 call wrtout(unit, sjoin(" kzone:", gstore%kzone))
 call wrtout(unit, sjoin(" kfilter:", gstore%kfilter))
 call wrtout(unit, sjoin(" nkibz:", itoa(gstore%nkibz)))
 call wrtout(unit, sjoin(" nkbz:", itoa(gstore%nkbz)))
 call wrtout(unit, sjoin(" glob_nk_spin:", ltoa(gstore%glob_nk_spin)))
 call wrtout(unit, sjoin(" qzone:", gstore%qzone))
 call wrtout(unit, sjoin(" nqibz:", itoa(gstore%nqibz)))
 call wrtout(unit, sjoin(" nqbz:", itoa(gstore%nqbz)))
 call wrtout(unit, sjoin(" glob_nq_spin:", ltoa(gstore%glob_nq_spin)))
 call wrtout(unit, sjoin(" kptopt:", itoa(gstore%ebands%kptopt)))
 call wrtout(unit, sjoin(" qptopt:", itoa(gstore%qptopt)))
 call wrtout(unit, sjoin(" with_vk:", itoa(gstore%with_vk)))

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   call wrtout(unit, sjoin(ch10, " === MPI distribution ==="))
   call wrtout(unit, sjoin("P Number of CPUs for parallelism over perturbations: ", itoa(gqk%pert_comm%nproc)))
   call wrtout(unit, sjoin("P Number of perturbations treated by this CPU: ",  itoa(gqk%my_npert)))
   call wrtout(unit, sjoin("P Number of CPUs for parallelism over q-points: ", itoa(gqk%qpt_comm%nproc)))
   call wrtout(unit, sjoin("P Number of CPUs for parallelism over k-points: ", itoa(gqk%kpt_comm%nproc)))
   call wrtout(unit, sjoin(" gqk_cplex:", itoa(gqk%cplex)), pre_newlines=1)
   call wrtout(unit, sjoin(" gqk_bstart:", itoa(gqk%bstart)))
   call wrtout(unit, sjoin(" gqk_bstop:", itoa(gqk%bstop)))
   call wrtout(unit, sjoin(" gqk_nb:", itoa(gqk%nb)))
   call wrtout(unit, sjoin(" gqk_my_npert:", itoa(gqk%my_npert)))
   call wrtout(unit, sjoin(" gqk_my_nk:", itoa(gqk%my_nk)))
   call wrtout(unit, sjoin(" gqk_my_nq:", itoa(gqk%my_nq)))
   !if (gqk%cplex == 1) then
   !  write(msg,'(a,f8.1,a)')'- Local memory allocated for |g|^2 array: ',ABI_MEM_MB(gqk%my_g2),' [Mb] <<< MEM'
   !  call wrtout(unit, msg) !; print *, "my_g2 shape:", shape(gqk%my_g2)
   !else
   !  write(msg,'(a,f8.1,a)')'- Local memory allocated for g array: ',ABI_MEM_MB(gqk%my_g),' [Mb] <<< MEM'
   !  call wrtout(unit, msg) !; print *, "my_g shape:", shape(gqk%my_g)
   !end if
 end do

end subroutine gstore_print
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_malloc__
!! NAME
!! gstore_malloc__
!!
!! FUNCTION
!! Allocate memory.
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
 integer,intent(in) :: qglob2bz(max_nq, gstore%nsppol)
 integer,intent(in) :: kglob2bz(max_nk, gstore%nsppol)
 integer,intent(in) :: qbz2ibz(6, gstore%nqbz), kbz2ibz(6, gstore%nkbz)

!Local variables-------------------------------
!scalars
 integer :: spin, my_is, ierr, my_iq, my_ik, iq_glob, iq_bz, ik_glob, ik_bz
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 logical :: isirr_k
 real(dp) :: mem_mb
 type(gqk_t),pointer :: gqk
!arrays
 integer,allocatable :: myq2glob(:), myk2glob(:)
!----------------------------------------------------------------------

 !gstore%max_nb = maxval(gstore%brange_spin(2, :) - gstore%brange_spin(1, :) + 1)

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   ! Split q-points and transfer symmetry tables.
   ! Note that glob_nq and glob_nk does not necessarily correspond to the size of the BZ
   ! First of all we have to consider kzone
   ! Even if kzone == "bz" we may have filtered the wavevectors e.g. Fermi surface.
   call xmpi_split_block(gqk%glob_nq, gqk%qpt_comm%value, gqk%my_nq, myq2glob)
   ABI_CHECK(gqk%my_nq > 0, "my_nq == 0")
   gqk%my_qstart = myq2glob(1)
   ABI_FREE(myq2glob)

   ABI_MALLOC(gqk%my_q2ibz, (6, gqk%my_nq))
   ABI_MALLOC(gqk%my_q2bz, (gqk%my_nq))

   do my_iq=1,gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1
     iq_bz = qglob2bz(iq_glob, spin)
     gqk%my_q2ibz(:, my_iq) = qbz2ibz(:, iq_bz)
     gqk%my_q2bz(my_iq) = iq_bz
   end do

   ! Split k-points and transfer symmetry tables
   call xmpi_split_block(gqk%glob_nk, gqk%kpt_comm%value, gqk%my_nk, myk2glob)
   ABI_CHECK(gqk%my_nk > 0, "my_nk == 0")
   gqk%my_kstart = myk2glob(1)
   ABI_FREE(myk2glob)

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

     !! symrel^T convention for k
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
     mem_mb = with_cplex * gqk%my_npert * gqk%nb ** 2 * gqk%my_nq * gqk%my_nk * eight * b2Mb
     call wrtout(std_out, sjoin(" Local memory for e-ph matrix elements:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))

     ! The initialization with zero is important as not all the g are computed when we filter in k-space.
     ! Abinit postprocessing tools will operate of the full my_g array
     ! and we don't want to trigger floating point exceptions.
     select case (with_cplex)
     case (1)
       ABI_MALLOC_OR_DIE(gqk%my_g2, (gqk%my_npert, gqk%nb, gqk%my_nq, gqk%nb, gqk%my_nk), ierr)
       gqk%my_g2 = zero
     case (2)
        ABI_MALLOC_OR_DIE(gqk%my_g, (gqk%my_npert, gqk%nb, gqk%my_nq, gqk%nb, gqk%my_nk), ierr)
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
 class(gstore_t),target,intent(inout) :: gstore
 real(dp),intent(in) :: qbz(3, gstore%nqbz)
 real(dp),target,intent(in) :: kbz(3, gstore%nqbz), kibz(3, gstore%nkibz)
 integer,intent(in) :: qbz2ibz(6,gstore%nqbz), qibz2bz(gstore%nqibz)
 integer,intent(in) :: kbz2ibz(6,gstore%nkbz), kibz2bz(gstore%nkibz)
 integer,intent(out) :: select_qbz_spin(gstore%nqbz, gstore%nsppol)
 integer,intent(out) :: select_kbz_spin(gstore%nkbz, gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: tetra_opt0 = 0
 integer :: nsppol, ierr, cnt, spin, band, ib, ii, max_nb, all_nproc, my_rank, comm, ebands_timrev
 integer :: ik_bz, ik_ibz, iflag, nk_in_star
 real(dp) :: max_occ
 character(len=80) :: errorstring
 type(ebands_t),pointer :: ebands
 type(crystal_t),pointer :: cryst
 type(htetra_t) :: ktetra
!arrays
 integer,allocatable :: indkk(:), kstar_bz_inds(:)
 real(dp):: rlatt(3,3), klatt(3,3), delta_theta_ef(2) ! qpt(3),
 real(dp),allocatable :: eig_ibz(:)
!----------------------------------------------------------------------

 comm = gstore%comm
 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Use the tetrahedron method to filter k- and k+q points on the FS in metals
 ! and define gstore%brange_spin automatically.
 call wrtout(std_out, sjoin(" Filtering k-points using:", gstore%kfilter))

 ! NB: here we recompute brange_spin
 cryst => gstore%cryst
 ebands => gstore%ebands
 nsppol = gstore%nsppol

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 call ebands_get_bands_e0(ebands, ebands%fermie, gstore%brange_spin, ierr)
 ABI_CHECK(ierr == 0, "Error in ebands_get_bands_e0")

 ABI_MALLOC(indkk, (gstore%nkbz))
 indkk(:) = kbz2ibz(1, :)

 rlatt = ebands%kptrlatt; call matr3inv(rlatt, klatt)
 call htetra_init(ktetra, indkk, gstore%cryst%gprimd, klatt, kbz, gstore%nkbz, gstore%kibz, gstore%nkibz, &
                  ierr, errorstring, gstore%comm)
 ABI_CHECK(ierr == 0, errorstring)

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
 class(gstore_t),target,intent(inout) :: gstore
 real(dp),intent(in) :: qbz(3, gstore%nqbz)
 real(dp),target,intent(in) :: kbz(3, gstore%nqbz), kibz(3, gstore%nkibz)
 integer,intent(in) :: qbz2ibz(6,gstore%nqbz), qibz2bz(gstore%nqibz)
 integer,intent(in) :: kbz2ibz(6,gstore%nkbz), kibz2bz(gstore%nkibz)
 integer,intent(out) :: select_qbz_spin(gstore%nqbz, gstore%nsppol)
 integer,intent(out) :: select_kbz_spin(gstore%nkbz, gstore%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: tetra_opt0 = 0
 integer :: nsppol, cnt, spin, ii, all_nproc, my_rank, comm, ebands_timrev, band ! ierr,
 integer :: ik_bz, ik_ibz, iflag, nk_in_star, gap_err
 logical :: assume_gap
 real(dp) :: ee, abs_erange1, abs_erange2, vmax, cmin
 type(ebands_t),pointer :: ebands
 type(crystal_t),pointer :: cryst
 type(gaps_t) :: gaps
!arrays
 integer,allocatable :: kstar_bz_inds(:)
!----------------------------------------------------------------------

 comm = gstore%comm
 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! filter k and k + q points according to erange and define gstore%brange_spin automatically.
 ! NB: here we recompute brange_spin
 call wrtout(std_out, sjoin(" Filtering k-points using gstore_erange:", &
                            ltoa(reshape(gstore%erange_spin, [2 * gstore%nsppol]) * Ha_eV), "(eV)"))

 cryst => gstore%cryst; ebands => gstore%ebands; nsppol = gstore%nsppol

 assume_gap = .not. all(gstore%erange_spin < zero)
 gaps = ebands_get_gaps(ebands, gap_err)
 if (assume_gap) call gaps%print(unit=std_out) !, header=msg)

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 select_kbz_spin = 0; cnt = 0

 do spin=1,nsppol
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
         if (ee <= vmax .and. vmax - ee <= abs_erange1) then
           iflag = 1 !; write(std_out, *), "Adding valence band", band, " with ee [eV]: ", ee * Ha_eV
         end if
       end if
       if (abs_erange2 > zero) then
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

end subroutine gstore_filter_erange__
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
 integer :: all_nproc, my_rank, ierr
 integer :: ii, iq_bz, iq_ibz, ikq_ibz, ikq_bz, len_kpts_ptr, ebands_timrev
!arrays
 integer,allocatable :: map_kq(:,:)
 real(dp):: qpt(3)
 real(dp),contiguous, pointer :: kpts_ptr(:,:)

! *************************************************************************

 ABI_UNUSED(kbz2ibz)
 ABI_UNUSED(qbz2ibz)

 all_nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)
 ebands_timrev = kpts_timrev_from_kptopt(gstore%ebands%kptopt)

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

     if (kpts_map("symrel", ebands_timrev, gstore%cryst, gstore%krank_ibz, len_kpts_ptr, kpts_ptr, map_kq, qpt=qpt) /= 0) then
       ABI_ERROR("Cannot map k+q to IBZ!")
     end if

     do ii=1,len_kpts_ptr
       ikq_ibz = map_kq(1, ii)
       ikq_bz = kibz2bz(ikq_ibz)
       select_qbz_spin(iq_bz, :) = select_qbz_spin(iq_bz, :) + select_kbz_spin(ikq_bz, :)
     end do
   end do ! iq_ibz

 case ("bz")
   do iq_bz=1,gstore%nqbz
     if (mod(iq_bz, all_nproc) /= my_rank) cycle ! MPI parallelism.
     qpt = qbz(:, iq_bz)
     !iq_ibz = qbz2ibz(1, iq_bz)
     ! k + q_bz --> k IBZ --> k BZ

     if (kpts_map("symrel", ebands_timrev, gstore%cryst, gstore%krank_ibz, len_kpts_ptr, kpts_ptr, map_kq, qpt=qpt) /= 0) then
       ABI_ERROR("Cannot map k+q to IBZ!")
     end if

     do ii=1,len_kpts_ptr
       ikq_ibz = map_kq(1, ii)
       ikq_bz = kibz2bz(ikq_ibz)
       select_qbz_spin(iq_bz, :) = select_qbz_spin(iq_bz, :) + select_kbz_spin(ikq_bz, :)
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
!! INPUTS
!!
!! OUTPUT
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
!!  Fills the bks_mask array defining the set of states that should be read from file
!!  by this MPI rank when computing the e-ph matrix elements.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_fill_bks_mask(gstore, mband, nkibz, nsppol, bks_mask)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: mband, nkibz, nsppol
 logical,intent(out) :: bks_mask(mband, nkibz, nsppol)

!Local variables-------------------------------
!scalars
 integer :: my_is, my_ik, my_iq, spin, ik_ibz, iqk_ibz, ebands_timrev
 real(dp) :: weight_q, cpu, wall, gflops
 type(gqk_t),pointer :: gqk
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: map_kq(:,:)
 real(dp) :: qpt(3)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 bks_mask = .False.
 cryst => gstore%cryst

 ebands_timrev = kpts_timrev_from_kptopt(gstore%ebands%kptopt)

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   ! We need the image of this k in the IBZ.
   do my_ik=1,gqk%my_nk
     ik_ibz = gqk%my_k2ibz(1, my_ik)
     bks_mask(gqk%bstart:gqk%bstart + gqk%nb - 1, ik_ibz, spin) = .True.
   end do

   ! as well as the image of k+q in the IBZ.
   ABI_MALLOC(map_kq, (6, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     call gqk%myqpt(my_iq, gstore, weight_q, qpt)

     if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, map_kq, qpt=qpt) /= 0) then
       ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qpt)))
     end if

     do my_ik=1,gqk%my_nk
       iqk_ibz = map_kq(1, my_ik)
       bks_mask(gqk%bstart:gqk%bstop, iqk_ibz, spin) = .True.
     end do
   end do

   ABI_FREE(map_kq)
 end do

 call cwtime_report(" gstore_fill_bks_mask", cpu, wall, gflops)

end subroutine gstore_fill_bks_mask
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_get_mpw_gmax
!! NAME
!! gstore_get_mpw_gmax
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_get_mpw_gmax(gstore, ecut, mpw, gmax)

!Arguments ------------------------------------
 class(gstore_t),target,intent(in) :: gstore
 real(dp),intent(in) :: ecut
 integer,intent(out) :: mpw, gmax(3)

!Local variables-------------------------------
 integer,parameter :: istwfk1 = 1
 integer :: my_is, my_ik, my_iq, spin, onpw, ii, ipw, ierr, my_mpw
 real(dp) :: weight_q, cpu, wall, gflops !weight_k,
 type(gqk_t),pointer :: gqk
!arrays
 integer :: my_gmax(3)
 integer,allocatable :: gtmp(:,:)
 real(dp) :: kk(3), qpt(3)

!----------------------------------------------------------------------

 mpw = 0; gmax = 0

 ! TODO: This is an hotspot due to the double loop over k and q.
 ! Should use a geometrical approach to compute mpw and gmax.

 call wrtout(std_out, " Computing mpw. This may take some time for dense k/q meshes...")
 call cwtime(cpu, wall, gflops, "start")

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik)

     ! Compute G sphere, returning npw. Note istwfk == 1.
     call get_kg(kk, istwfk1, ecut, gstore%cryst%gmet, onpw, gtmp)
     mpw = max(mpw, onpw)
     do ipw=1,onpw
       do ii=1,3
         gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
       end do
     end do
     ABI_FREE(gtmp)

     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_q, qpt)

       ! TODO: g0 umklapp here can enter into play!
       ! gmax could not be large enough!
       call get_kg(kk + qpt, 1, ecut, gstore%cryst%gmet, onpw, gtmp)
       mpw = max(mpw, onpw)
       do ipw=1,onpw
         do ii=1,3
          gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
         end do
       end do
       ABI_FREE(gtmp)
     end do

   end do ! my_ik
 end do ! my_is

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, gstore%comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, gstore%comm, ierr)

 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))
 call cwtime_report(" gmax and mpw", cpu, wall, gflops)

end subroutine gstore_get_mpw_gmax
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_calc_my_phonons
!! NAME
!! gstore_calc_my_phonons
!!
!! FUNCTION
!!  Compute and store ph frequencies and optionally the displacements for all
!!  the q-points treated by this MPI rank.
!!  TODO: Rewrite this part, Read stuff from the IBZ and symmetrize.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_calc_my_phonons(gstore, store_phdispl)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 logical,intent(in) :: store_phdispl

!Local variables-------------------------------
 integer :: my_is, my_iq, ierr
 real(dp) :: weight_q, cpu, wall, gflops, qpt(3)
 type(crystal_t),pointer :: cryst
 type(gqk_t),pointer :: gqk
 real(dp),allocatable :: phfrq(:), displ_cart(:,:,:,:)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
 cryst => gstore%cryst

 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, 3 * cryst%natom))
 ABI_MALLOC(phfrq, (3 * cryst%natom))

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)

   ABI_RECALLOC(gqk%my_wnuq, (gqk%my_npert, gqk%my_nq))
   if (store_phdispl) then
     ABI_RECALLOC(gqk%my_displ_cart, (2, 3, cryst%natom, gqk%my_npert, gqk%my_nq))
   end if

   do my_iq=1,gqk%my_nq
     if (gqk%kpt_comm%skip(my_iq)) cycle ! MPI parallelism inside kpt_comm

     ! Get phonon frequencies and eigenvectors for this q-point.
     ! FIXME: Perhaps it's better if we compute everything at q_ibz and then rotate
     ! so that we are guaranteed to have e(-q) == e(q)^*
     call gqk%myqpt(my_iq, gstore, weight_q, qpt)
     call gstore%ifc%fourq(cryst, qpt, phfrq, displ_cart)

     gqk%my_wnuq(:, my_iq) = phfrq(gqk%my_iperts(:))
     if (store_phdispl) gqk%my_displ_cart(:,:,:,:,my_iq) = displ_cart(:,:,:,gqk%my_iperts(:))
   end do

   ! Note that the computation is replicated across the perturbation communicator.
   ! This should not represent a problem since the phase in the displacement is fixed in dfpt_phfrq
   ! hence results should be coherent across different procs.
   call xmpi_sum(gqk%my_wnuq, gqk%kpt_comm%value, ierr)
   if (store_phdispl) call xmpi_sum(gqk%my_displ_cart, gqk%kpt_comm%value, ierr)
 end do ! my_is

 ABI_FREE(displ_cart)
 ABI_FREE(phfrq)

 call cwtime_report(" gstore_calc_my_phonons", cpu, wall, gflops)

end subroutine gstore_calc_my_phonons
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

subroutine gstore_get_lambda_iso_iw(gstore, dtset, nw, imag_w, lambda)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: nw
 real(dp),intent(in) :: imag_w(nw)
 real(dp),intent(out) :: lambda(nw)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, ib_k, ib_kq, ierr
 real(dp) :: g2, w_nuq, weight_k, weight_q
 type(gqk_t), pointer :: gqk
!arrays
 real(dp) :: qpt(3)
 real(dp),allocatable :: dbldelta_q(:,:,:), g2_pmnk(:,:,:,:)

!----------------------------------------------------------------------

 ABI_CHECK(gstore%qzone == "bz", "gstore_get_lambda_iso_iw assumes qzone == `bz`")
 lambda = zero

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)

   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbldelta_q, (gqk%nb, gqk%nb, gqk%my_nk))
   ABI_MALLOC(g2_pmnk, (gqk%my_npert, gqk%nb, gqk%nb, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     ! Compute integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, dtset%eph_intmeth, dtset%eph_fsmear, qpt, weight_q, dbldelta_q)

     ! Copy data to improve memory access in the loops below.
     g2_pmnk = gqk%my_g2(:,:,my_iq,:,:)

     do my_ik=1,gqk%my_nk
       !weight_k = gqk%my_kweight(my_ik, gstore)
       weight_k = gqk%my_wtk(my_ik)
       do ib_k=1,gqk%nb
         do ib_kq=1,gqk%nb
           do my_ip=1,gqk%my_npert
             g2 = g2_pmnk(my_ip, ib_kq, ib_k, my_ik)
             ! TODO: handle w_nuq ~ 0
             w_nuq = gqk%my_wnuq(my_ip, my_iq)
             lambda(:) = lambda(:) + &
               two * w_nuq / (imag_w(:) ** 2 + w_nuq ** 2) * g2 * weight_k * weight_q * dbldelta_q(ib_kq, ib_k, my_ik)
           end do
         end do
       end do
     end do
   end do ! my_iq

   ABI_FREE(dbldelta_q)
   ABI_FREE(g2_pmnk)
 end do ! my_is

 ! Take into account collinear spin
 lambda = lambda * (two / (gstore%nsppol * dtset%nspinor))

 call xmpi_sum(lambda, gstore%comm, ierr)

end subroutine gstore_get_lambda_iso_iw
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_get_a2fw
!! NAME
!! gstore_get_a2fw
!!
!! FUNCTION
!!  Compute a^2F(omega)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_get_a2fw(gstore, dtset, nw, wmesh, a2fw)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: nw
 real(dp),intent(in) :: wmesh(nw)
 real(dp),intent(out) :: a2fw(nw)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, ib_k, ib_kq, ierr
 real(dp) :: g2, w_nuq, weight_k, weight_q, cpu, wall, gflops
 type(gqk_t), pointer :: gqk
!arrays
 real(dp) :: qpt(3)
 real(dp),allocatable :: dbldelta_q(:,:,:), g2_mnkp(:,:,:,:), deltaw_nuq(:)

!----------------------------------------------------------------------

 call wrtout(std_out, sjoin(" Computing a^2F(w) with ph_smear:", ftoa(dtset%ph_smear * Ha_meV), "(meV)"))
 ABI_CHECK(gstore%qzone == "bz", "gstore_get_lambda_iso_iw assumes qzone == `bz`")
 call cwtime(cpu, wall, gflops, "start")

 ABI_MALLOC(deltaw_nuq, (nw))
 a2fw = zero

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)

   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbldelta_q, (gqk%nb, gqk%nb, gqk%my_nk))
   ABI_MALLOC(g2_mnkp, (gqk%nb, gqk%nb, gqk%my_nk, gqk%my_npert))

   do my_iq=1,gqk%my_nq
     ! Compute integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, dtset%eph_intmeth, dtset%eph_fsmear, qpt, weight_q, dbldelta_q)

     ! Copy data to improve memory access in the loops below.
     do my_ip=1,gqk%my_npert
       g2_mnkp(:,:,:,my_ip) = gqk%my_g2(my_ip,:,my_iq,:,:)
     end do

     do my_ip=1,gqk%my_npert
       w_nuq = gqk%my_wnuq(my_ip, my_iq)
       deltaw_nuq = gaussian(wmesh - w_nuq, dtset%ph_smear)
       do my_ik=1,gqk%my_nk
         !weight_k = gqk%my_kweight(my_ik, gstore)
         weight_k = gqk%my_wtk(my_ik)
         do ib_k=1,gqk%nb
           do ib_kq=1,gqk%nb
             g2 = g2_mnkp(ib_kq, ib_k, my_ik, my_ip)
             a2fw(:) = a2fw(:) + deltaw_nuq(:) * g2 * weight_k * weight_q * dbldelta_q(ib_kq, ib_k, my_ik)
           end do
         end do
       end do
     end do
   end do ! my_iq

   ABI_FREE(dbldelta_q)
   ABI_FREE(g2_mnkp)
 end do ! my_is

 ABI_FREE(deltaw_nuq)

 ! Take into account collinear spin and N(eF) TODO
 a2fw = a2fw * (two / (gstore%nsppol * dtset%nspinor))
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
!!  Free dynamic memory.
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
 ABI_SFREE(gstore%erange_spin)

 call gstore%krank_ibz%free()

end subroutine gstore_free
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_myqpt
!! NAME
!! gqk_myqpt
!!
!! FUNCTION
!!  Return the weight and the reduced coordinates of the q-point from the local index my_iq
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
!scalars
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
!!  Note that k/q weights are not included in dbldelta_q
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gqk_dbldelta_qpt(gqk, my_iq, gstore, eph_intmeth, eph_fsmear, qpt, weight_q, dbldelta_q)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: my_iq, eph_intmeth
 real(dp),intent(in) :: eph_fsmear
 real(dp),intent(out) :: qpt(3), weight_q, dbldelta_q(gqk%nb, gqk%nb, gqk%my_nk)

!Local variables ------------------------------
!scalars
 real(dp), parameter :: min_smear = tol9
 integer :: nb, nkbz, spin, my_ik, ib, ib1, ib2, band1, band2, nesting, ebands_timrev
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 integer :: ikq_ibz, isym_kq, trev_kq, tsign_kq, g0_kq(3)
 integer :: ii, i1, i2, i3, cnt, ik_bz, ltetra
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

 nb = gqk%nb; nkbz = gstore%nkbz; spin = gqk%spin

 ebands => gstore%ebands; cryst => gstore%cryst
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

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
   if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, my_kqmap, qpt=qpt) /= 0) then
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
         do ib=1,nb
           vb_k(:,ib) = tsign_k * matmul(transpose(cryst%symrel_cart(:,:,isym_k)), vb_k(:,ib))
         end do
       end if
       if (.not. isirr_kq) then
         do ib=1,nb
           vb_kq(:,ib) = tsign_kq * matmul(transpose(cryst%symrel_cart(:,:,isym_kq)), vb_kq(:,ib))
         end do
       end if
     end if

     do ib2=1,nb
       band2 = ib2 + gqk%bstart - 1
       if (use_adaptive) then
         sigma = max(maxval([(abs(dot_product(vb_k(:, ib2), kmesh_cartvec(:,ii))), ii=1,3)]), min_smear)
         !write(std_out, *)"sigma:", sigma * Ha_eV
       end if
       g2 = gaussian(ebands%eig(band2, ik_ibz, spin) - ebands%fermie, sigma)

       do ib1=1,nb
         band1 = ib1 + gqk%bstart - 1
         if (use_adaptive) then
           sigma = max(maxval([(abs(dot_product(vb_kq(:, ib1), kmesh_cartvec(:,ii))), ii=1,3)]), min_smear)
         end if
         g1 = gaussian(ebands%eig(band1, ikq_ibz, spin) - ebands%fermie, sigma)
         dbldelta_q(ib1, ib2, my_ik) = g1 * g2 ! / fs%nktot
       end do

     end do
   end do

   ABI_FREE(my_kqmap)

 else if (abs(eph_intmeth) == 2) then

   ABI_CHECK(isdiagmat(ebands%kptrlatt), "kptrlatt must be diagonal when tetra is used.")
   ABI_CHECK(ebands%nshiftk == 1, "nshiftk must be 1 when tetra is used")
   nge = get_diag(ebands%kptrlatt); ngw = nge
   ABI_CHECK(nkbz == product(nge(1:3)), "Wrong nge")

   ! Compute eig_k and eig_kq in full BZ for the relevant bands around Ef.
   ABI_MALLOC(kmesh, (3, nkbz))
   ABI_MALLOC(eig_k, (nb, nkbz))
   ABI_MALLOC(eig_kq, (nb, nkbz))

   ! Technical problems:
   !
   ! 1) libtetrabz works with the BZ and assume a certaing ordering of the k-points (see below)
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
   if (kpts_map("symrec", ebands_timrev, cryst, gstore%krank_ibz, nkbz, kmesh, kmesh_map) /= 0) then
     ABI_ERROR("Cannot map libtetra mesh to IBZ")
   end if

   do ik_bz=1,nkbz
     ik_ibz = kmesh_map(1, ik_bz)
     eig_k(:, ik_bz) = ebands%eig(gqk%bstart:gqk%bstop, ik_ibz, spin) - ebands%fermie
   end do

   ! Map libtetra BZ mesh + q to IBZ and fill eig_kq
   if (kpts_map("symrec", ebands_timrev, cryst, gstore%krank_ibz, nkbz, kmesh, kmesh_map, qpt=qpt) /= 0) then
     ABI_ERROR(sjoin("Cannot map libtetra k+q to IBZ with qpt:", ktoa(qpt)))
   end if

   do ik_bz=1,nkbz
     ikq_ibz = kmesh_map(1, ik_bz)
     eig_kq(:, ik_bz) = ebands%eig(gqk%bstart:gqk%bstop, ikq_ibz, spin) - ebands%fermie
   end do

   ABI_FREE(kmesh_map)
   !call cwtime_report(" kmesh_map", cpu, wall, gflops)

   ! Call libtetra routine to compute weights for double delta integration.
   ! Note that libtetra assumes Ef set to zero.
   ! TODO: Average weights over degenerate states?
   ! NB: This is a bootleneck, can pass comm_kp

   ! Select option for double delta with tetra.
   !  2 for the optimized tetrahedron method.
   ! -2 for the linear tetrahedron method.
   ltetra = 0
   if (eph_intmeth ==  2) ltetra = 2
   if (eph_intmeth == -2) ltetra = 1

   ABI_MALLOC(wght_bz, (nb, nb, nkbz))
   call libtetrabz_dbldelta(ltetra, gstore%cryst%gprimd, nb, nge, eig_k, eig_kq, ngw, wght_bz) !, comm=comm)
   !call cwtime_report(" libtetrabz_dbldelta", cpu, wall, gflops)

   my_krank = krank_new(gqk%my_nk, gqk%my_kpts)

   ! Reindex from full BZ to my set of kpoints and rescale weights.
   cnt = 0
   do ik_bz=1,nkbz
     my_ik = my_krank%get_index(kmesh(:, ik_bz))
     if (my_ik /= -1) then
       dbldelta_q(:,:,my_ik) = wght_bz(:,:,ik_bz) * gstore%nkbz
       cnt = cnt + 1
     end if
   end do

   ! FIXME: bug if k-point (and q-point) parallelism.
   ABI_CHECK(cnt == gqk%my_nk, sjoin("cnt != my_nk, ", itoa(cnt), itoa(gqk%my_nk)))
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
!!  Free dynamic memory.
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
 ABI_SFREE(gqk%my_wnuq)
 ABI_SFREE(gqk%my_displ_cart)
 ABI_SFREE(gqk%my_g)
 ABI_SFREE(gqk%my_g2)
 ABI_SFREE(gqk%my_iperts)
 ABI_SFREE(gqk%vk_cart_ibz)
 ABI_SFREE(gqk%vkmat_cart_ibz)

 ! Free communicators
 call gqk%kpt_comm%free()
 call gqk%qpt_comm%free()
 call gqk%pert_comm%free()
 call gqk%qpt_pert_comm%free()
 call gqk%grid_comm%free()

end subroutine gqk_free
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
 integer,parameter :: tim_getgh1c = 1, berryopt0 = 0, ider0 = 0, idir0 = 0, LOG_MODQ = 5
 integer,parameter :: useylmgr = 0, useylmgr1 = 0, master = 0, ndat1 = 1
 integer :: my_rank,nproc,mband,nsppol,nkibz,idir,ipert,ebands_timrev, iq_bz
 integer :: cplex,natom,natom3,ipc,nspinor, nskip_tetra_kq
 integer :: bstart_k,bstart_kq,nband_k,nband_kq,band_k, ib_k, ib_kq !ib1,ib2, band_kq,
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq
 integer :: my_ik, my_is, comm_rpt, my_npert, my_ip, my_iq, spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: mpw, nb,ierr,cnt, n1,n2,n3,n4,n5,n6,nspden,ndone
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf, nkpg, nkpg1, qbuf_size, iqbuf_cnt, root_ncid, spin_ncid, ncerr
 integer :: ii, my_nqibz, iq_start, iq_ibz, isym_q, trev_q, prev_iqbz
 real(dp) :: cpu, wall, gflops, cpu_q, wall_q, gflops_q, cpu_all, wall_all, gflops_all
 real(dp) :: ecut, eshift, eig0nk, weight_q, weight_k
 logical :: gen_eigenpb, isirr_k, isirr_kq, isirr_q, print_time
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(ddkop_t) :: ddkop
 type(gqk_t),pointer :: gqk
 character(len=500) :: msg
!arrays
 integer :: g0_k(3), g0_kq(3), g0_q(3), work_ngfft(18),gmax(3),indkk_kq(6,1)
 integer(i1b),allocatable :: itreat_qibz(:)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), nband(:,:), wfd_istwfk(:), qselect(:)
 integer,allocatable :: my_pinfo(:,:), pert_table(:,:)
 integer,allocatable :: iq_buf(:,:), done_qbz_spin(:,:), my_iqibz_inds(:)
 real(dp) :: kk_bz(3),kq_bz(3),kk_ibz(3),kq_ibz(3), qq_bz(3), qq_ibz(3), vk(3)
 real(dp) :: phfrq(3*cryst%natom), ylmgr_dum(1,1,1)
 real(dp),allocatable :: displ_cart_qibz(:,:,:,:), displ_red_qibz(:,:,:,:), pheigvec_qibz(:,:,:,:)
 real(dp),allocatable :: displ_cart_qbz(:,:,:,:), displ_red_qbz(:,:,:,:), pheigvec_qbz(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:), kinpw1(:), kpg1_k(:,:), kpg_k(:,:), dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:), ffnl1(:,:,:,:), ph3d(:,:,:), ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:), gkk_atm(:,:,:,:),gkq_nu(:,:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:), kets_k(:,:,:), h1kets_kq(:,:,:), cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:), vlocal(:,:,:,:), vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:), ylm_k(:,:), ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:), gvnlx1(:,:), work(:,:,:,:)
 real(dp),allocatable :: gs1c(:,:), vk_cart_ibz(:,:,:) !, vkmat_cart_ibz(:,:,:,:)
 real(dp),allocatable :: my_gbuf(:,:,:,:,:,:), buf_wqnu(:,:), buf_eigvec_cart(:,:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)

!************************************************************************

 ! This parameter defines the size of the q-buffer used to store the g(k, q) e-ph matrix elements
 ! for all the k-point treated by this MPI rank.
 ! Increasing the buffer size increases the memory requirements
 ! but it leads to better performance as the number of IO operations is decreased.
 ! TODO: Should compute it on the basis of my_nkpt and my_nqpt
 qbuf_size = 4
 call wrtout(std_out, sjoin(" Begin computation of e-ph matrix elements with qbuf_size:", itoa(qbuf_size)))

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkibz = ebands%nkpt; mband = ebands%mband
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ! FFT meshes
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   if (gqk%pert_comm%nproc > 1) then
     ! Activate parallelism over perturbations
     ! Build table with list of perturbations treated by this CPU inside pert_comm
     ABI_WARNING("pert_comm%nproc > 1 not tested")
     call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, gqk%pert_comm%value)
     call dvdb%set_pert_distrib(my_npert, natom3, my_pinfo, pert_table, gqk%pert_comm%value)
     ABI_CHECK(all(my_pinfo(3, :) == gqk%my_iperts), "my_pinfo(3, :) != gqk%my_iperts")
     ABI_FREE(my_pinfo)
     ABI_FREE(pert_table)
   end if
 end do

 !call wrtout([std_out, ab_out], " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
 ! Prepare Fourier interpolation of DFPT potentials.
 comm_rpt = xmpi_comm_self
 !comm_rpt = bqs_comm%value
 call dvdb%ftinterp_setup(dtset%ddb_ngqpt, gstore%qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)

 ! Build q-cache in the *dense* IBZ using the global mask qselect and itreat_qibz.
 ABI_MALLOC(itreat_qibz, (gstore%nqibz))
 ABI_MALLOC(qselect, (gstore%nqibz))
 qselect = 0; itreat_qibz = 0
 call dvdb%ftqcache_build(nfftf, ngfftf, gstore%nqibz, gstore%qibz, dtset%dvdb_qcache_mb, qselect, itreat_qibz, gstore%comm)
 ABI_FREE(itreat_qibz)
 ABI_FREE(qselect)

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k/k+q wavevectors treated by this MPI rank are stored.

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

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
   nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print(header="Wavefunctions for GSTORE calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

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

 ! Spherical Harmonics for useylm == 1.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

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

 call init_hamiltonian(gs_hamkq, psps, pawtab, nspinor, nsppol, nspden, natom, &
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg, &
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab, &
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 ! Allocate vlocal. Note nvloc
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(displ_cart_qbz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_cart_qibz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red_qbz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red_qibz, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(pheigvec_qbz, (2, 3, cryst%natom, 3*cryst%natom))
 ABI_MALLOC(pheigvec_qibz, (2, 3, cryst%natom, 3*cryst%natom))

 ABI_CALLOC(dummy_vtrial, (nfftf, nspden))

 ! Create ddkop object to compute group velocities (if needed)
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ! Open GSTORE file, and read table used for restarting.
 NCF_CHECK(nctk_open_modify(root_ncid, gstore%path, gstore%comm))

 ABI_MALLOC(done_qbz_spin, (gstore%nqbz, nsppol))
 NCF_CHECK(nf90_get_var(root_ncid, nctk_idname(root_ncid, "gstore_done_qbz_spin"), done_qbz_spin))
 gstore%wfk0_path = wfk0_path
 if (my_rank == master) then
   NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_wfk0_path"), trim(gstore%wfk0_path)))
 end if

 if (my_rank == master) call gstore%print(std_out)

 ndone = count(done_qbz_spin == 1)

 ! NB: Write phonon data here as we are not guaranteed to have all the IBZ q-points
 ! inside the loop over my_iq if filtering has been used.
 ! TODO: Write phdata only if eph_task == -11
 if (ndone == 0) then
   call wrtout(std_out, " Computing phonon frequencies and displacements in the IBZ")
   call cwtime(cpu, wall, gflops, "start")

   call xmpi_split_block(gstore%nqibz, gstore%comm, my_nqibz, my_iqibz_inds)
   ABI_MALLOC(buf_wqnu, (natom3, my_nqibz))
   ABI_MALLOC(buf_eigvec_cart, (2, 3, natom, natom3, my_nqibz))

   do ii=1,my_nqibz
     iq_ibz = my_iqibz_inds(ii)
     !call ifc%fourq(cryst, gstore%qibz(:, iq_ibz), buf_wqnu(:,ii), buf_displ_cart(:,:,:,:,ii))
     call ifc%fourq(cryst, gstore%qibz(:, iq_ibz), buf_wqnu(:,ii), displ_cart_qibz, &
                    out_eigvec=buf_eigvec_cart(:,:,:,:,ii))
   end do
   if (nproc > 1 .and. gstore%nqibz >= nproc) then
     NCF_CHECK(nctk_set_collective(root_ncid, root_vid("phfreqs_ibz")))
     NCF_CHECK(nctk_set_collective(root_ncid, root_vid("pheigvec_cart_ibz")))
   end if

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
   call cwtime_report(" phonon computation + output", cpu, wall, gflops)
 else
   call wrtout(std_out, &
               sjoin(" Restarting GSTORE calculation. Found: ", itoa(ndone), " (qpt, spin) entries already computed"))
 end if

 if (gstore%with_vk /= 0 .and. ndone == 0) then
   call wrtout(std_out, " Computing and writing velocity operator matrix elements in the IBZ")
   call wrtout(std_out, " Note that not all the k-points in the IBZ are computed when kfilter is activated!")
   call cwtime(cpu, wall, gflops, "start")

   ! On disk, we have:
   !
   !    nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz"))
   !    nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz")))

   ABI_MALLOC(cgwork, (2, mpw*wfd%nspinor))

   do my_is=1,gstore%my_nspins
     spin = gstore%my_spins(my_is)
     gqk => gstore%gqk(my_is)

     if (gstore%with_vk == 1) then
       ABI_CALLOC(vk_cart_ibz, (3, gqk%nb, gstore%nkibz))
     else
       ABI_ERROR("with_vk 2")
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
         do ib_k=1,gqk%nb
           band_k = ib_k + gqk%bstart - 1
           call wfd%copy_cg(band_k, ik_ibz, spin, cgwork)
           vk = ddkop%get_vdiag(ebands%eig(band_k, ik_ibz, spin), istwf_k, npw_k, wfd%nspinor, cgwork, cwaveprj0)
           vk_cart_ibz(:, ib_k, ik_ibz) = vk
         end do

       case (2)
         ABI_ERROR("with_vk 2")
         !do ib_k=1,nband_k
         !  band_k = ib_k + bstart_k - 1
         !end do
       end select
     end do ! my_ik

     call xmpi_sum_master(vk_cart_ibz, master, gqk%grid_comm%value, ierr)
     if (gqk%grid_comm%me == master) then
       NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))
       NCF_CHECK(nf90_put_var(spin_ncid, spin_vid("vk_cart_ibz"), vk_cart_ibz))
     end if
     ABI_SFREE(vk_cart_ibz)
   end do ! my_is

   ABI_FREE(cgwork)
   call cwtime_report(sjoin(" Computation of v_k group velocities with with_vk:", itoa(gstore%with_vk)), cpu, wall, gflops)
 end if

 ! Loop over my spins.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   my_npert = gqk%my_npert
   NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

   ! Allocate workspace for wavefunctions using mpw and nb
   nb = gqk%nb
   ABI_MALLOC(bras_kq, (2, mpw*nspinor, nb))
   ABI_MALLOC(kets_k, (2, mpw*nspinor, nb))
   ABI_MALLOC(h1kets_kq, (2, mpw*nspinor, nb))
   ABI_MALLOC(gkk_atm, (2, nb, nb, natom3))
   ABI_MALLOC(gkq_nu, (2, nb, nb, natom3))
   ABI_MALLOC(iq_buf, (2, qbuf_size))

   ! Inside the loops we compute gkq_nu(2, nb, nb, natom3)
   ABI_MALLOC_OR_DIE(my_gbuf, (gqk%cplex, nb, nb, natom3, gqk%my_nk, qbuf_size), ierr)

   ! Loop over my set of q-points
   prev_iqbz = -1
   do my_iq=1,gqk%my_nq
     print_time = my_rank == 0 .and. (my_iq <= LOG_MODQ .or. mod(my_iq, LOG_MODQ) == 0)
     if (print_time) call cwtime(cpu_q, wall_q, gflops_q, "start")
     iq_bz = gqk%my_q2bz(my_iq)
     if (done_qbz_spin(iq_bz, spin) == 1) cycle

     call gqk%myqpt(my_iq, gstore, weight_q, qq_bz)

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
       call ifc%fourq(cryst, qq_ibz, phfrq, displ_cart_qibz, out_displ_red=displ_red_qibz, out_eigvec=pheigvec_qibz)
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

     !call ifc%fourq(cryst, qq_bz, phfrq, displ_cart_qbz, out_displ_red=displ_red_qbz, out_eigvec=pheigvec_qbz))
     ! Use Fourier interpolation of DFPT potentials to get my_npert potentials.
     !cplex = 2
     !ABI_MALLOC(v1scf, (cplex, nfft, nspden, dvdb%my_npert))
     !call dvdb%ftinterp_qpt(qq_bz, nfftf, ngfftf, v1scf, dvdb%comm_rpt)

     ! Version with qcache.
     call dvdb%get_ftqbz(cryst, qq_bz, qq_ibz, gqk%my_q2ibz(:, my_iq), cplex, nfftf, ngfftf, v1scf, &
                         gqk%pert_comm%value)

     ! Allocate vlocal1 with correct cplex. Note nvloc and my_npert.
     ABI_MALLOC(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert))

     ! Set up local potential vlocal1 with proper dimensioning from vtrial1 taking into account the spin.
     do my_ip=1,my_npert
       call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc,&
                 pawfgr, mpi_enreg, dummy_vtrial, v1scf(:,:,:,my_ip), vlocal, vlocal1(:,:,:,:,my_ip))
     end do

     ! Continue to initialize the GS Hamiltonian
     call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

     ! Loop over my k-points
     do my_ik=1,gqk%my_nk
       ! Set entry to zero. Important as there are cycle instructions inside these loops
       ! and we don't want to write random numbers to disk.
       my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = zero

       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       kk_bz = gqk%my_kpts(:, my_ik)
       weight_k = gqk%my_wtk(my_ik)

       ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
       trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5,my_ik)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       kk_ibz = ebands%kptns(:,ik_ibz)

       ! Number of bands crossing the Fermi level at k
       !bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)
       bstart_k = gqk%bstart; nband_k = gqk%nb

       ! ============================================
       ! Find symmetrical image of k+q in in the kIBZ
       ! ============================================
       kq_bz = kk_bz + qq_bz

       if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, 1, kq_bz, indkk_kq) /= 0) then
         write(msg, '(3a)' ) &
          "Cannot find k+q in kmesh", ch10, 'Action: check your WFK file and the (k, q) point input variables.'
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

       ! Number of bands crossing the Fermi level at k+q
       !bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz)
       bstart_kq = gqk%bstart; nband_kq = gqk%nb
       ABI_CHECK(nband_k <= nb .and. nband_kq <= nb, "wrong nband")

       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kk_bz, kk_ibz, bstart_k, nband_k, spin, mpw, gqk%my_k2ibz(:, my_ik), cryst, &
                          work_ngfft, work, istwf_k, npw_k, kg_k, kets_k)

       ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kq_bz, kq_ibz, bstart_kq, nband_kq, spin, mpw, indkk_kq(:,1), cryst, &
                          work_ngfft, work, istwf_kq, npw_kq, kg_kq, bras_kq)

       ! if PAW, one has to solve a generalized eigenproblem
       ! Be careful here because I will need sij_opt==-1
       gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2, npw_kq*nspinor*((sij_opt+1)/2)))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm == 1) then
       !   optder = 0; if (useylmgr == 1) optder = 1
       !   call initylmg(cryst%gprimd, kg_k, kk_bz, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_k], dtset%nsppol, optder, cryst%rprimd, ylm_k, ylmgr)
       !   call initylmg(cryst%gprimd, kg_kq, kq_bz, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_kq], dtset%nsppol, optder, cryst%rprimd, ylm_kq, ylmgr_kq)
       !end if

       ! Compute k+G vectors
       nkpg = 3 * dtset%nloalg(3)
       ABI_MALLOC(kpg_k, (npw_k, nkpg))
       if (nkpg > 0) call mkkpg(kg_k, kpg_k, kk_bz, nkpg, npw_k)

       ! Compute nonlocal form factors ffnlk at (k+G)
       ABI_MALLOC(ffnlk, (npw_k, 1, psps%lmnmax, psps%ntypat))
       call mkffnl_objs(cryst, psps, 1, ffnlk, ider0, idir0, kg_k, kpg_k, kk_bz, nkpg, npw_k, ylm_k, ylmgr_dum, &
                        comm=gqk%pert_comm%value) !, request=ffnlk_request)

       ! Compute k+q+G vectors
       nkpg1 = 3 * dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq_bz, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G)
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       call mkffnl_objs(cryst, psps, 1, ffnl1, ider0, idir0, kg_kq, kpg1_k, kq_bz, nkpg1, npw_kq, ylm_kq, ylmgr_kq, &
                        comm=gqk%pert_comm%value) ! request=ffnl1_request)

       ! Loop over my atomic perturbations and compute gkk_atm.
       gkk_atm = zero
       do my_ip=1,my_npert
         idir = dvdb%my_pinfo(1, my_ip); ipert = dvdb%my_pinfo(2, my_ip); ipc = dvdb%my_pinfo(3, my_ip)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)

         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,my_ip), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk_bz, kq_bz, idir, ipert, &             ! In
                            cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
                            npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
                            dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1, &  ! Out
                            reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)              ! Reuse some arrays

         ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0)
         do ib_k=1,nband_k
           band_k = ib_k + bstart_k - 1
           eig0nk = ebands%eig(band_k, ik_ibz, spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0, kets_k(:,:,ib_k), cwaveprj0, h1kets_kq(:,:,ib_k), &
                        grad_berry, gs1c, gs_hamkq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
                        optnl, opt_gvnlx1, rf_hamkq, sij_opt, tim_getgh1c, usevnl)
         end do

         call rf_hamkq%free()
         ABI_FREE(kinpw1)
         ABI_FREE(dkinpw)
         ABI_FREE(ph3d)
         ABI_SFREE(ph3d1)

         ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
         ! No need to handle istwf_kq because it's always 1.
         do ib_k=1,nband_k
           do ib_kq=1,nband_kq
             gkk_atm(:, ib_kq, ib_k, ipc) = cg_zdotc(npw_kq*nspinor, bras_kq(1,1,ib_kq), h1kets_kq(1,1,ib_k))
           end do
         end do

       end do ! my_ip

       ABI_FREE(gs1c)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)

       ! Collect gkk_atm inside pert_comm so that all procs can operate on the data.
       if (gqk%pert_comm%nproc > 1) call xmpi_sum(gkk_atm, gqk%pert_comm%value, ierr)

       ! Get g in the phonon representation.
       call ephtk_gkknu_from_atm(nb, nb, 1, natom, gkk_atm, phfrq, displ_red_qbz, gkq_nu)

       ! Save e-ph matrix elements in the buffer.
       select case (gqk%cplex)
       case (1)
         my_gbuf(1,:,:,:, my_ik, iqbuf_cnt) = gkq_nu(1,:,:,:) ** 2 + gkq_nu(2,:,:,:) ** 2
       case (2)
         my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = gkq_nu
       end select
     end do ! my_ik

     ABI_FREE(v1scf)
     ABI_FREE(vlocal1)

     ! Dump buffers
     if (iqbuf_cnt == qbuf_size) call dump_data()

     if (print_time) then
       write(msg,'(2(a,i0),a)')" My q-point [", my_iq, "/", gqk%my_nq, "]"
       call cwtime_report(msg, cpu_q, wall_q, gflops_q); if (my_iq == LOG_MODQ) call wrtout(std_out, "...", do_flush=.True.)
     end if
     !if (my_rank == master) then
     !  ! Print cache stats.
     !  !call dvdb%ft_qcache%report_stats()
     !end if
   end do ! my_iq

   ! Dump the remainder.
   if (iqbuf_cnt /= 0) call dump_data()

   ABI_FREE(iq_buf)
   ABI_FREE(my_gbuf)
   ABI_FREE(bras_kq)
   ABI_FREE(kets_k)
   ABI_FREE(h1kets_kq)
   ABI_FREE(gkk_atm)
   ABI_FREE(gkq_nu)
 end do ! my_is

 call cwtime_report(" GSTORE computation done", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10) !, comm=gstore%comm)
 call gstore%print(std_out, header="GSTORE at the end of gstore%compute")

 !call xmpi_barrier(gstore%comm)
 NCF_CHECK(nf90_close(root_ncid))

 ! Free memory
 ABI_FREE(gvnlx1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(displ_cart_qbz)
 ABI_FREE(displ_cart_qibz)
 ABI_FREE(displ_red_qbz)
 ABI_FREE(displ_red_qibz)
 ABI_FREE(pheigvec_qbz)
 ABI_FREE(pheigvec_qibz)
 ABI_FREE(done_qbz_spin)

 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 call ddkop%free()
 call gs_hamkq%free()
 call wfd%free()

contains

subroutine dump_data()

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

 ! If parallelism over pertubation is activated, only the procs treating the first perturbation
 ! i.e. the procs treating different k-points for this q are involved in IO
 ! as all the local buffers store results for all natom3 pertubations.

 integer :: ii, iq_bz, iq_glob, my_iq

 if (gqk%coords_qkp(3) /= 0) goto 10 ! Yes, I'm very proud of this GOTO.

 !iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]
 my_iq = iq_buf(1, 1)
 iq_glob = my_iq + gqk%my_qstart - 1

 ! NB: this is an individual IO operation
 ncerr = nf90_put_var(spin_ncid, spin_vid("gvals"), my_gbuf, &
                      start=[1, 1, 1, 1, gqk%my_kstart, iq_glob], &
                      count=[gqk%cplex, gqk%nb, gqk%nb, gqk%natom3, gqk%my_nk, iqbuf_cnt] &
                     )
 NCF_CHECK(ncerr)

 ! Only one proc sets the entry in done_qbz_spin to 1 for all the q-points in the buffer.
 if (all(gqk%coords_qkp(2:3) == [0, 0]))  then
   do ii=1,iqbuf_cnt
     iq_bz = iq_buf(2, ii)
     NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_done_qbz_spin"), 1, start=[iq_bz, spin]))
   end do
 end if

 ! Zero the counter before returning
10 iqbuf_cnt = 0

end subroutine dump_data

integer function root_vid(vname)
  character(len=*),intent(in) :: vname
  root_vid = nctk_idname(root_ncid, vname)
end function root_vid

integer function spin_vid(vname)
  character(len=*),intent(in) :: vname
  spin_vid = nctk_idname(spin_ncid, vname)
end function spin_vid

end subroutine gstore_compute
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_from_ncpath
!! NAME
!! gstore_from_ncpath
!!
!! FUNCTION
!!  Reconstruct a gstore object from a GSTORE.nc netcdf file.
!!
!! INPUTS
!!
!! SOURCE

subroutine gstore_from_ncpath(gstore, path, with_cplex, dtset, cryst, ebands, ifc, comm)

!Arguments ------------------------------------
 class(gstore_t),target,intent(out) :: gstore
 character(len=*),intent(in) :: path
 integer,intent(in) :: with_cplex
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: comm
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),target,intent(in) :: ebands
 class(ifc_type),target,intent(in) :: ifc

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_rank, ncid, spin, spin_ncid, nproc, ierr, fform, max_nb, ib, natom, natom3
 integer :: max_nq, max_nk, gstore_cplex, ncerr, my_is, my_iq, iq_glob, my_ik, ik_glob, my_ip, ipert
 integer :: iq_ibz, isym_q, trev_q, tsign_q, g0_q(3)
 real(dp) :: cpu, wall, gflops
 character(len=10) :: priority
 logical :: store_phdispl, isirr_q
 type(hdr_type) :: wfk0_hdr
 type(crystal_t) :: gstore_cryst
 type(gqk_t),pointer :: gqk
!arrays
 integer :: ibuffer(9), nproc_spin(ebands%nsppol), comm_spin(ebands%nsppol), brange_spin(2, ebands%nsppol)
 integer,allocatable :: qglob2bz(:,:), kglob2bz(:,:), qbz2ibz(:,:), kbz2ibz(:,:)
 real(dp) :: qq_ibz(3)
 real(dp),allocatable :: gwork_q(:,:,:,:,:), slice_bb(:,:,:)
 real(dp),allocatable :: phfreqs_ibz(:,:), pheigvec_cart_ibz(:,:,:,:,:)
 real(dp),allocatable :: pheigvec_cart_qbz(:,:,:,:), displ_cart_qbz(:,:,:,:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Set basic parameters.
 gstore%comm = comm
 gstore%nsppol = dtset%nsppol
 gstore%path = path

 ! Get references to other data structures.
 gstore%cryst => cryst
 gstore%ebands => ebands
 gstore%ifc => ifc
 gstore%kibz => ebands%kptns
 natom = cryst%natom
 natom3 = cryst%natom * 3

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
   !NCF_CHECK(cryst%ncwrite(ncid))

   ! TODO Should compare ebands_file with input ebands
   !NCF_CHECK(ebands_ncwrite(ebands, ncid))

   call hdr_ncread(wfk0_hdr, ncid, fform)
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
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kzone"), gstore%kzone))
   call replace_ch0(gstore%kzone)
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qzone"), gstore%qzone))
   call replace_ch0(gstore%qzone)
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kfilter"), gstore%kfilter))
   call replace_ch0(gstore%kfilter)
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_wfk0_path"), gstore%wfk0_path))
   call replace_ch0(gstore%wfk0_path)

   ABI_MALLOC(gstore%qibz, (3, gstore%nqibz))
   ABI_MALLOC(gstore%wtq, (gstore%nqibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_brange_spin"), brange_spin))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_erange_spin"), gstore%erange_spin))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qibz"), gstore%qibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_wtq"), gstore%wtq))

   ABI_MALLOC(qbz2ibz, (6, gstore%nqbz))
   ABI_MALLOC(kbz2ibz, (6, gstore%nkbz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qbz2ibz"), qbz2ibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kbz2ibz"), kbz2ibz))

   ABI_MALLOC(qglob2bz, (max_nq, gstore%nsppol))
   ABI_MALLOC(kglob2bz, (max_nk, gstore%nsppol))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qglob2bz"), qglob2bz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_kglob2bz"), kglob2bz))
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
   !call gstore_cryst%print(header="crystal structure from WFK file")

   if (cryst%compare(gstore_cryst, header=" Comparing input crystal with the one from GSTORE file") /= 0) then
     ABI_ERROR("Crystal structure from input and GSTORE do not agree! Check messages above!")
   end if
   call gstore_cryst%free()

   if (gstore_cplex == 1 .and. with_cplex == 2) then
     ABI_ERROR("GSTORE file contains |g| while Abinit needs complex g. Regenerate GSTORE file with gstore_cplex 2")
   end if
 end if ! master

 if (nproc > 1) then
   ! Broadcast the header.
   call wfk0_hdr%bcast(master, my_rank, comm)

   ! Broadcast dimensions.
   if (my_rank == master) then
     ibuffer = [gstore_cplex, gstore%nkibz, gstore%nkbz, gstore%nqibz, gstore%nqbz, &
                gstore%with_vk, max_nq, max_nk, max_nb]
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
     ABI_MALLOC(qbz2ibz, (6, gstore%nqbz))
     ABI_MALLOC(kbz2ibz, (6, gstore%nkbz))
     ABI_MALLOC(qglob2bz, (max_nq, gstore%nsppol))
     ABI_MALLOC(kglob2bz, (max_nk, gstore%nsppol))
   end if

   call xmpi_bcast(gstore%qptopt, master, comm, ierr)
   call xmpi_bcast(gstore%kzone, master, comm, ierr)
   call xmpi_bcast(gstore%qzone, master, comm, ierr)
   call xmpi_bcast(gstore%kfilter, master, comm, ierr)
   call xmpi_bcast(gstore%wfk0_path, master, comm, ierr)
   call xmpi_bcast(brange_spin, master, comm, ierr)
   call xmpi_bcast(gstore%erange_spin, master, comm, ierr)
   call xmpi_bcast(gstore%qibz, master, comm, ierr)
   call xmpi_bcast(gstore%wtq, master, comm, ierr)
   call xmpi_bcast(qbz2ibz, master, comm, ierr)
   call xmpi_bcast(kbz2ibz, master, comm, ierr)
   call xmpi_bcast(qglob2bz, master, comm, ierr)
   call xmpi_bcast(kglob2bz, master, comm, ierr)
   call xmpi_bcast(gstore%glob_nq_spin, master, comm, ierr)
   call xmpi_bcast(gstore%glob_nk_spin, master, comm, ierr)

   if (gstore%kfilter == "fs_tetra") then
     if (my_rank /= master) then
       ABI_MALLOC(gstore%delta_ef_kibz_spin, (max_nb, gstore%nkibz, gstore%nsppol))
     end if
     call xmpi_bcast(gstore%delta_ef_kibz_spin, master, comm, ierr)
   end if

 end if

 ! Construct crystal and ebands from the GS WFK file.
 !ebands_file = wfk_read_ebands(wfk0_path, comm, out_hdr=wfk0_hdr)
 call wfk0_hdr%vs_dtset(dtset)
 call wfk0_hdr%free()

 ! Distribute spins, create indirect mapping to spin index and init gstore%brange_spin
 call gstore%distribute_spins__(ebands%mband, brange_spin, nproc_spin, comm_spin, comm)

 ! Compute krank
 gstore%krank_ibz = krank_from_kptrlatt(gstore%nkibz, gstore%kibz, ebands%kptrlatt, compute_invrank=.False.)

 ! Set MPI grid. Also define priority according to eph_task.
 call priority_from_eph_task(dtset%eph_task, priority)

 call gstore%set_mpi_grid__(with_cplex, dtset%eph_np_pqbks, priority, nproc_spin, comm_spin)

 ! At this point, we have the Cartesian grid (one per spin if any)
 ! and we can finally allocate and distribute other arrays.
 call gstore%malloc__(with_cplex, max_nq, qglob2bz, max_nk, kglob2bz, qbz2ibz, kbz2ibz)

 ABI_FREE(qglob2bz)
 ABI_FREE(kglob2bz)
 ABI_FREE(qbz2ibz)
 ABI_FREE(kbz2ibz)

 call xmpi_comm_free(comm_spin)

 ! Now we read the big arrays with MPI-IO and hdf5 groups.
 ! Note the loop over spin as each gqk has its own dimensions.
 ! Recall the shape on the arrays:
 !
 ! In memory, we have allocated:
 !
 !    my_g(my_npert, nb, my_nq, nb, my_nk) if with_cplex == 2 (complex array)
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

 ! ========================================
 ! Load phonon frequencies and eigenvectors
 ! ========================================

 ! nctkarr_t("phfreqs_ibz", "dp", "natom3, gstore_nqibz")
 ! nctkarr_t("pheigvec_cart_ibz", "dp", "two, three, natom, natom3, gstore_nqibz")
 ABI_MALLOC(phfreqs_ibz, (natom3, gstore%nqibz))
 ABI_MALLOC(pheigvec_cart_ibz, (2, 3, cryst%natom, cryst%natom * 3, gstore%nqibz))
 ABI_MALLOC(pheigvec_cart_qbz, (2, 3, cryst%natom, cryst%natom * 3))
 ABI_MALLOC(displ_cart_qbz, (2, 3, cryst%natom, cryst%natom * 3))

 NCF_CHECK(nctk_open_read(ncid, gstore%path, gstore%comm))

 if (nproc > 1) then
   NCF_CHECK(nctk_set_collective(ncid, vid("phfreqs_ibz")))
 end if
 NCF_CHECK(nf90_get_var(ncid, vid("phfreqs_ibz"), phfreqs_ibz))
 if (nproc > 1) then
   NCF_CHECK(nctk_set_collective(ncid, vid("pheigvec_cart_ibz")))
 end if
 NCF_CHECK(nf90_get_var(ncid, vid("pheigvec_cart_ibz"), pheigvec_cart_ibz))
 NCF_CHECK(nf90_close(ncid))

 store_phdispl = .True.
 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)

   ABI_MALLOC(gqk%my_wnuq, (gqk%my_npert, gqk%my_nq))
   if (store_phdispl) then
     ABI_MALLOC(gqk%my_displ_cart, (2, 3, cryst%natom, gqk%my_npert, gqk%my_nq))
   end if

   do my_iq=1,gqk%my_nq
     iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
     trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
     !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     isirr_q = (isym_q == 1 .and. trev_q == 0)
     tsign_q = 1; if (trev_q == 1) tsign_q = -1
     qq_ibz = gstore%qibz(:, iq_ibz)
     call pheigvec_rotate(cryst, qq_ibz, isym_q, trev_q, pheigvec_cart_ibz(:,:,:,:,iq_ibz), &
                          pheigvec_cart_qbz, displ_cart_qbz)

     gqk%my_wnuq(:, my_iq) = phfreqs_ibz(gqk%my_iperts(:), iq_ibz)
     if (store_phdispl) gqk%my_displ_cart(:,:,:,:,my_iq) = displ_cart_qbz(:,:,:,gqk%my_iperts(:))
   end do ! my_iq
 end do ! my_is

 ABI_FREE(phfreqs_ibz)
 ABI_SFREE(pheigvec_cart_ibz)
 ABI_SFREE(displ_cart_qbz)
 ABI_FREE(pheigvec_cart_qbz)

 do spin=1,gstore%nsppol
   my_is = gstore%spin2my_is(spin)

   if (my_is /= 0) then
     gqk => gstore%gqk(my_is)

     NCF_CHECK(nctk_open_read(ncid, gstore%path, gqk%grid_comm%value))
     NCF_CHECK(nf90_inq_ncid(ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

     !NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("bstart"), gqk%bstart))
     !NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("bstop"), gqk%bstop))

     ! gstore_cplex defines the data on disk while
     ! cplex defines what we want to store in memory
     ABI_MALLOC_OR_DIE(gwork_q, (gstore_cplex, gqk%nb, gqk%nb, gqk%natom3, gqk%glob_nk), ierr)
     ABI_MALLOC(slice_bb, (gstore_cplex, gqk%nb, gqk%nb))

     do my_iq=1,gqk%my_nq
        iq_glob = my_iq + gqk%my_qstart - 1

        ! Read q-slice (individual IO)
        ncerr = nf90_get_var(spin_ncid, spin_vid("gvals"), gwork_q, start=[1, 1, 1, 1, 1, iq_glob]) ! count=[])
        NCF_CHECK(ncerr)

        do my_ik=1,gqk%my_nk
          ik_glob = my_ik + gqk%my_kstart - 1
          do my_ip=1,gqk%my_npert
            ipert = gqk%my_iperts(my_ip)
            slice_bb = gwork_q(:,:,:, ipert, ik_glob)

            ! Put data in the right place and handle conversion g --> |g|^2
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

     ! ==============================================
     ! Read matrix elements of the velocity operator
     ! ==============================================
     if (gstore%with_vk == 1) then
       if (gqk%grid_comm%nproc > 1) then
         NCF_CHECK(nctk_set_collective(spin_ncid, spin_vid("vk_cart_ibz")))
       end if
       NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("vk_cart_ibz"), gqk%vk_cart_ibz))

     else if (gstore%with_vk == 2) then
       if (gqk%grid_comm%nproc > 1) then
         NCF_CHECK(nctk_set_collective(spin_ncid, spin_vid("vkmat_cart_ibz")))
       end if
       NCF_CHECK(nf90_get_var(spin_ncid, spin_vid("vkmat_cart_ibz"), gqk%vkmat_cart_ibz))

       ! Tranfer diagonal terms to vk_cart_ibz.
       do ib=1,gqk%nb
         gqk%vk_cart_ibz(:, ib, :) = gqk%vkmat_cart_ibz(1, :, ib, ib, :)
       end do
     end if

     NCF_CHECK(nf90_close(ncid))
   end if

   call xmpi_barrier(gstore%comm)
 end do ! spin

 call cwtime_report(" gstore_from_ncpath", cpu, wall, gflops)

contains
integer function vid(vname)
  character(len=*),intent(in) :: vname
  vid = nctk_idname(ncid, vname)
end function vid

integer function spin_vid(vname)
  character(len=*),intent(in) :: vname
  spin_vid = nctk_idname(spin_ncid, vname)
end function spin_vid

end subroutine gstore_from_ncpath
!!***

end module m_gstore
!!***

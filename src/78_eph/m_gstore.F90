!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!
!!  This module implements the gstore object that allows one to
!!  **precompute"" the e-ph matrix elements g and store them in memory
!!  with a MPI-distributed data structure.
!!
!!  This approach is the most CPU-efficient when one has to deal with
!!  algorithms in which the same g(q, k) is required several (many) times.
!!  Typical examples are iterative solvers for non-linear equations that are called inside a loop over T.
!!  At each iteration, indeed, we need g(q, k) and computing these quantities from scratch
!!  would be very expensive.
!!
!!  Note that g depends on the two wave vectors (q, k), two electron band indices (m, n),
!!  phonon mode nu with momentum q and spin if nsppol == 2 (collinear case).
!!
!!  g(q, k) is therefore a sloppy notation for:
!!
!!          g(q, k) = <k + q, m, spin| \Delta_{q, \nu} V^{spin}_{scf} | k, n, spin>
!!
!!  There are lots of technical details that should be discussed but, roughly speaking,
!!  the gqk API allows one to:
!!
!!   - select whether q or k should be in the IBZ or in the BZ.
!!     NB: It is not possible to use the IBZ both for q and k as g(Sk, q) = g(k, Sq)
!!     thus one has to select the appropriate sampling beforehand
!!
!!   - filter bands and/or k/q wavevectors according to some criterion.
!!     In superconductors, for instance, only states on the Fermi surface are needed
!!
!!  - whether the code should compute and store the complex valued g or |g|^2.
!!    Expression depending of |g|^2 are gauge-invariant provided that all degenerate states are summed over.
!!    On the contrary, the complex valued g is gauge-dependent and hic sunt leones.
!!    In Abinit, the g elements are computed within the same gauge by reconstructing Bloch states
!!    in the BZ from the IBZ by using a deterministic symmetrization algorithm
!!    Client code reading the e-ph matrix elements produced by ABINIT is expected to follow the
!!    same conventions, especially if they need to mix g matrix elements with wavefunctions in the BZ.
!!
!!  At the level of the API, we have three different steps.
!!
!!      1) call gstore_new to define the BZ sampling type (e.g. k in the IBZ, q in the BZ)
!!         and additional filtering techniques. The MPI grid is automatically generated at this level
!!
!!      2) call gstore_compute to compute the e-ph matrix elements
!!
!!      3) use the gstore object to implement your equations or write the results to nc file.
!!
!!  Now, let's discuss the MPI-distribution.
!!
!!  The (k, q) matrix is distributed inside a 2D cartesian grid using block distribution.
!!  This is schematic representation for 4 procs with 2 procs for k and 2 procs for q:
!!
!!                      k-axis
!!              |--------------------
!!              |         |         |
!!     q-axis   |   P00   |   P01   |
!!              |         |         |
!!              |--------------------
!!              |         |         |
!!              |   P10   |   P11   |
!!              |         |         |
!!              |--------------------
!!
!!  Each processor stores all the (band_1, band_2) transitions for a given (q, k) pair.
!!  Perturbations can be optionally distributed along a third axis (pert_comm).
!!  The parallelism over perturbations is not expected to be the most efficient but it allows
!!  one to reduce the memory required to store the scattering potential in the supercell
!!  as we can distribute W(r, R, 3*natom) over the last dimension.
!!
!!  NB: If nsppol == 2, we create two gqk objects, one for each spin.
!!  The reason is that dimensions such as the number of effective bands/q-points/k-points
!!  depends on the spin if we start to filter.
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
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
 use m_krank
 use m_htetra

 use m_ebands
 use iso_c_binding
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
 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, get_diag !, isdiagmat
 use m_io_tools,       only : iomode_from_fname
 use m_copy,           only : alloc_copy
 use m_fftcore,        only : ngfft_seq, get_kg
 use m_cgtools,        only : cg_zdotc
 use m_kg,             only : getph, mkkpg
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_symtk,          only : matr3inv
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_ifc,            only : ifc_type
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
!!  This object Stores MPI-distributed e-p matrix elements for
!!  a given spin index (collinear case).
!!
!! NOTES
!!  local dimensions start with `my_`,
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

  integer :: bstart = 1
  ! The first band starts at bstart. (global index)

  integer :: my_npert = -1
  ! Number of perturbations treated by this MPI rank.

  integer :: glob_nk = -1, glob_nq = -1
  ! Total number of k/q points in global matrix.
  ! Note that k-points/q-points can be filtered.
  ! Use kzone and qzone, kfilter to interpret these dimensions

  integer :: my_nk = -1, my_nq = -1
  ! Number of k/q points treated by this MPI proc.
  ! Used to loop and allocate local arrays.

  integer :: my_kstart = -1, my_qstart = -1
  ! Index of the first k/q point in the global matrix treated by this MPI proc

  !real(dp),allocatable :: my_wtk(:)
  ! (my_nk)
  ! Weights for k-points treated by this MPI rank

  !real(dp),allocatable :: my_wtq(:)
  ! (my_nq)
  ! Weights for q-points treated by this MPI rank

  integer,allocatable :: my_k2ibz(:,:)
  ! (6, my_nk)
  ! Mapping my_kpoints --> kibz

  integer,allocatable :: my_q2ibz(:,:)
  ! (6, my_nq)
  ! Mapping my_qpoints --> qibz

  complex(dp),allocatable :: my_vk_cart(:,:,:)
  ! (3, nb, my_nk)
  ! Diagonal v_{m, m,k} for the k-points treated by this MPI proc.
  ! Allocated if gstore%with_vk == 1

  complex(dp),allocatable :: my_vkmat_cart(:,:,:,:)
  ! (3, nb, nb, my_nk)
  ! v_{m, n,k} for the k-points treated by this MPI proc.
  ! Allocated if gstore%with_vk == 2

  !complex(dp),allocatable :: my_vcart_qk(:,:,:,:)
  ! group velocities v_{nm,q+k} for the k/q-points treated by this MPI proc.
  ! (3, nb, nb, my_nq, my_nk)

  ! or alternatively:
  !complex(dp),allocatable :: my_vcart_kibz(:,:,:,:)

  integer,allocatable :: my_iperts(:)
  ! (my_npert)
  ! List of perturbation indices treated by this MPI proc.

  real(dp), allocatable :: my_vals(:,:,:,:,:,:)
  ! E-ph matrix elements
  ! (cplex, my_npert, nb, my_nq, nb, my_nk)

  integer :: coords_qkp(3)
  ! coordinates of this processor in the Cartesian grid.

  type(xcomm_t) :: kpt_comm
   ! MPI communicator over k-points

  type(xcomm_t) :: qpt_comm
   ! MPI communicator over q-points

  type(xcomm_t) :: pert_comm
   ! MPI communicator over atomic perturbations.

  type(xcomm_t) :: grid_comm
   ! MPI communicator for full (q,k,pert) grid.

  real(dp),allocatable :: my_wqnu(:,:)
  ! (my_npert, my_nq)
  ! Phonon frequencies (MPI distributed)

  real(dp),allocatable :: my_displ_cart(:,:,:,:,:)
  ! (2, 3, cryst%natom, my_npert, my_nq))
  ! Phonon displacement (MPI distributed)
  ! Seldom needed because e-ph matrix elements are already in the phonon representation.

 contains

  procedure :: get_mykpt => gqk_get_mykpt
  ! Return the k-point from my local index my_ik

  procedure :: get_all_mykpts => gqk_get_all_mykpts
  ! Allocate and return array with the all reduced coordinates of the k-point
  ! treated by this MPI rank

  procedure :: get_myqpt => gqk_get_myqpt
  ! Return the q-point from my local index my_iq

  procedure :: ncwrite_path => gqk_ncwrite_path

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
!!
!! NOTES
!!
!! SOURCE

type, public :: gstore_t

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: my_nspins
   ! Number of collinear spins treated by this MPI rank

  integer :: nkibz = -1, nqibz = -1
  ! Number of k/q points in the IBZ.

  integer :: nkbz = -1, nqbz = -1
  ! Number of k/q points in the BZ.

  integer :: comm
   ! Global communicator
   ! Inherited by the caller so we don't free it in gstore_free.

  integer :: with_vk = 0
  ! 0 if group velocities should not be computed
  ! 1 to compute diagonal terms only
  ! 2 to compute diagonal + off-diagonal terms
  !
  ! NB: We compute v_k for all the k-points treated by this MPI proc.
  ! as v_kq can be reconstructed by symmetry although this step requires an all_gather along the k-axis.

  character(len=nctk_slen) :: kzone = "", qzone = ""
    ! Specifies whether k- or q-points are in the BZ or in the IBZ.
    ! Possible values are "ibz" or "bz".
    ! Note that combination ("ibz", "ibz") is not allowed

  character(len=nctk_slen) :: kfilter = "none"
  ! Specifies the tecnique used to filter k-points.
  ! Possible values:
  !     "None"
  !     "fs_tetra"
  !     "fs_erange"
  !     "erange" ! referred to CBM, VBM in semiconductors.
  !     ...

  !real(dp) :: erange(2) = zero
  ! Energy window

  type(crystal_t), pointer :: cryst  => null()

  type(ebands_t), pointer :: ebands  => null()

  type(ifc_type), pointer :: ifc => null()

  type(krank_t) :: krank !, qrank
  ! Object used to find the correspondence with the kibz array with the IBZ.

  integer,allocatable :: my_spins(:)
   ! (%my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2.

  real(dp), contiguous, pointer :: kibz(:,:)
  ! k-points in the IBZ. Points to ebands%kptns
  ! (3, nkibz)

  real(dp), allocatable :: qibz(:,:)
  ! (3, nqibz)
  ! q-points in the IBZ

  real(dp), allocatable :: wtq(:)
  ! (nqibz)
  ! q-points weights in the IBZ

  !integer :: kptrlatt(3, 3), qptrlatt(3, 3)
   ! k-mesh and q-mesh

  !real(dp) :: kshift(3, 1), qshift(3, 1)
  ! k/q-mesh shift (well, q-mesh is usually gamma-centered)

  type(xcomm_t) :: spin_comm
    ! MPI communicator over spins (nsppol = 2)

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

  procedure :: free => gstore_free
  ! Free memory

  procedure :: ncwrite_path => gstore_ncwrite_path
  ! Write object to file

  !procedure :: ncread => gstore_ncread_path
  ! Reconstruct object from file

  procedure :: compute => gstore_compute
  ! Compute e-ph matrix elements.

  procedure :: calc_my_phonons => gstore_calc_my_phonons
  ! Helper function to compute ph quantities for all q-points treated by the MPI proc.

end type gstore_t

public :: gstore_new
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_new
!! NAME
!! gstore_new
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function gstore_new(dtset, cryst, ebands, ifc, comm) result (gstore)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),target,intent(in) :: ebands
 class(ifc_type),target,intent(in) :: ifc
 type(gstore_t), target :: gstore

!Local variables-------------------------------
!scalars
 integer,parameter :: qptopt1 = 1, timrev1 = 1, tetra_opt0 = 0, master = 0
 integer :: all_nproc, my_rank, ierr, my_nshiftq, nsppol, iq_glob, ik_glob, ii ! out_nkibz,
 integer :: my_is, my_ik, my_iq, spin, natom3, cnt, np, band, iflag ! bstart, bstop, nw,
 integer :: ik_ibz, ik_bz, ebands_timrev, nk_in_star  ! isym_k, trev_k,
 integer :: iq_bz, iq_ibz, ikq_ibz, ikq_bz, len_kpts_ptr, color  !isym_q, trev_q,
 real(dp) :: cpu, wall, gflops, dksqmax, max_occ, mem_mb ! , elow, ehigh, estep
 !logical :: isirr_k, isirr_q
 character(len=5000) :: msg
 character(len=80) :: errorstring
 type(gqk_t),pointer :: gqk
 type(krank_t),target :: qrank
 type(htetra_t) :: ktetra ! , qtetra
!arrays
 integer :: ngqpt(3), qptrlatt(3,3)
 integer :: comm_spin(ebands%nsppol), nproc_spin(ebands%nsppol), unts(2)
 integer :: glob_nk_spin(ebands%nsppol), glob_nq_spin(ebands%nsppol)
 integer :: bands_spin(2, ebands%nsppol), fs_bands_spin(2, ebands%nsppol)
 integer,allocatable :: myq2glob(:), myk2glob(:)
 integer,allocatable :: qbz2ibz_(:,:), kbz2ibz_(:,:), kibz2bz(:), qibz2bz(:), indkk(:), kstar_bz_inds(:)
 integer,allocatable :: qglob2bz_idx(:,:), kglob2bz_idx(:,:)
 integer,allocatable :: select_qbz_spin(:,:), select_kbz_spin(:,:), map_kq(:,:)
 real(dp):: my_shiftq(3,1), qpt(3), rlatt(3,3), klatt(3,3), delta_theta_ef(2)  ! kpt(3), qlatt(3,3),
 real(dp),allocatable :: qbz_(:,:), wtk_(:) !, my_kpts(:,:) !wtq_(:),
 real(dp),target,allocatable :: kibz_(:,:), kbz_(:,:)
 real(dp),allocatable :: eig_ibz(:)  !, wvals(:)
 real(dp),contiguous, pointer :: kpts_ptr(:,:)
 !integer :: out_kptrlatt(3,3)
 !real(dp),allocatable :: out_kibz(:,:), out_wtk(:)
#ifdef HAVE_MPI
 integer,parameter :: ndims = 3
 integer :: comm_cart, me_cart
 logical :: reorder
 integer :: dims(ndims)
 logical :: periods(ndims), keepdim(ndims)
#endif

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Set basic parameters.
 gstore%nsppol = nsppol
 gstore%comm = comm

 ! Get references to other data structures.
 gstore%cryst => cryst
 gstore%ebands => ebands
 gstore%ifc => ifc
 gstore%kibz => ebands%kptns

 ! Set metadata.
 gstore%kzone = "ibz"    ! = dtset%gstore_kzone; if (present(kzone)) gstore%kzone = kzone
 gstore%qzone = "bz"     ! = dtset%gstore_qzone; if (present(qzone)) gstore%qzone = qzone
 gstore%kfilter = "none" ! = dtset%gstore_kfilter; if (present(kfilter)) gstore%kfilter = kfilter
 gstore%with_vk = 0      ! = dtset%gstore_with_vk; if (present(with_vk)) gstore%with_vk = with_vk

 !gstore%kptrlatt(3, 3)
 !gstore%kshift(3, 1)
 !gstore%qptrlatt(3, 3)
 !gstore%qshift(3, 1)
 !gstore%spin_comm

 ! Distribute spins and create mapping to spin index.
 !gstore%spin_comm%nproc = 1
 !if (gstore%nsppol == 2 .and. mod(all_nproc, 2) == 0) then
 !  gstore%spin_comm%nproc = 2
 !end if
 !if (any(dtset%eph_np_pqbks /= 0)) gstore%spin_comm%nproc = dtset%eph_np_pqbks(5)

 if (gstore%nsppol == 2) then
   call xmpi_split_block(gstore%nsppol, gstore%comm, gstore%my_nspins, gstore%my_spins)
   !msg = sjoin("nsppol (", itoa(gstore%nsppol), ") < spin_comm_nproc (", itoa(gstore%spin_comm%nproc), ")")
   !ABI_CHECK(gstore%my_nspins > 0, msg)
 else
   ! No nsppol parallelism DOH!
   gstore%my_nspins = 1
   ABI_MALLOC(gstore%my_spins, (gstore%my_nspins))
   gstore%my_spins = 1
 end if

 !comm_spin(:) = comm

 do spin=1,nsppol
   ! NB: If MPI_UNDEFINED is passed as the colour value, the subgroup in which the calling
   ! MPI process will be placed is MPI_COMM_NULL
   color = xmpi_undefined; if (any(gstore%my_spins == spin)) color = 1
   call xmpi_comm_split(comm, color, my_rank, comm_spin(spin), ierr)
 end do

 ABI_MALLOC(gstore%gqk, (gstore%my_nspins))
 do spin=1,nsppol
   nproc_spin(spin) = xmpi_comm_size(comm_spin(spin))
   bands_spin(:, spin) = [1, ebands%mband]
   !if (all(dtsegt%gstore_brange /= 0)) bands_spin(:, spin) = dtset%gstore_brange(:, spin)

   ABI_CHECK_IRANGE(bands_spin(1, spin), 1, ebands%mband, "gstore_brange(1,spin)")
   ABI_CHECK_IRANGE(bands_spin(2, spin), 1, ebands%mband, "gstore_brange(2,spin)")
   ABI_CHECK(bands_spin(1, spin) <= bands_spin(2, spin), "bands_spin(1, spin) <= bands_spin(2, spin)")
 end do

 ! Define q-mesh
 ! Either q-mesh from DVDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation if q not in DDB)
 ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if

 ! TODO: Should fix bz2ibz to use the same conventions as krank and listkk
 ! NB: only sigmaph seems to be using this optional argument

 ! Setup qIBZ, weights and BZ.
 ! Always use q --> -q symmetry for phonons even in systems without inversion
 qptrlatt = 0; qptrlatt(1, 1) = ngqpt(1); qptrlatt(2, 2) = ngqpt(2); qptrlatt(3, 3) = ngqpt(3)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
                             gstore%nqibz, gstore%qibz, gstore%wtq, gstore%nqbz, qbz_)
                             !new_kptrlatt, new_shiftk,
                             !, bz2ibz=new%ind_qbz2ibz)

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(qbz2ibz_, (6, gstore%nqbz))

 qrank = krank_from_kptrlatt(gstore%nqibz, gstore%qibz, qptrlatt, compute_invrank=.False.)
 call qrank%get_mapping(gstore%nqbz, qbz_, dksqmax, cryst%gmet, qbz2ibz_, &
                        cryst%nsym, cryst%symafm, cryst%symrec, timrev1, use_symrec=.True.)

 if (dksqmax > tol12) then
   ABI_ERROR("Cannot map qBZ to IBZ!")
 end if
 call qrank%free()

 call get_ibz2bz(gstore%nqibz, gstore%nqbz, qbz2ibz_, qibz2bz, ierr)
 ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for q-points!")

 ! Get full BZ associated to ebands
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
   gstore%nkibz, kibz_, wtk_, gstore%nkbz, kbz_) !, bz2ibz=bz2ibz)

 ! In principle kibz_ should be equal to ebands%kptns
 ABI_CHECK(gstore%nkibz == ebands%nkpt, "nkibz != ebands%nkpt")
 ABI_FREE(kibz_)

 ! Note symrel and use_symrec=.False. in get_mapping.
 ! This means that this table can be used to symmetrize wavefunctions in cgtk_rotate.
 ! TODO This ambiguity should be removed. Change cgtk_rotate so that we can use the symrec convention.

 ABI_MALLOC(kbz2ibz_, (6, gstore%nkbz))
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 gstore%kibz => ebands%kptns

 gstore%krank = krank_from_kptrlatt(gstore%nkibz, gstore%kibz, ebands%kptrlatt, compute_invrank=.False.)
 call gstore%krank%get_mapping(gstore%nkbz, kbz_, dksqmax, cryst%gmet, kbz2ibz_, &
                        cryst%nsym, cryst%symafm, cryst%symrel, ebands_timrev, &
                        use_symrec=.False.)

 if (dksqmax > tol12) then
    ABI_ERROR("Cannot map kBZ to IBZ!")
 end if

 call get_ibz2bz(gstore%nkibz, gstore%nkbz, kbz2ibz_, kibz2bz, ierr)
 ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for k-points")

 ! These tables are used to exclude q/k points
 ! We use the full BZ because this mask can be also used when points are restricted to the IBZ
 ! Note that both arrays are initialized with zeros.
 ABI_ICALLOC(select_qbz_spin, (gstore%nqbz, nsppol))
 ABI_ICALLOC(select_kbz_spin, (gstore%nkbz, nsppol))

 if (gstore%kzone == "ibz" .and. gstore%qzone == "ibz") then
   ABI_ERROR("The combination kzone = 'ibz' and qzone = 'ibz' is not allowed")
 end if

 select case (gstore%kzone)
 case ("ibz")
   do spin=1,nsppol
     do ik_ibz=1,gstore%nkibz
       ik_bz = kibz2bz(ik_ibz)
       select_kbz_spin(ik_bz, spin) = 1
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
     cnt = 0
     do iq_ibz=1, gstore%nqibz
       iq_bz = qibz2bz(iq_ibz)
       select_qbz_spin(iq_bz, spin) = 1
     end do
     ABI_CHECK(cnt == gstore%nqibz, "cnt != nqibz. Something wrong with symmetry tables.")
   end do

 case ("bz")
    select_qbz_spin = 1

 case default
   ABI_ERROR(sjoin("Invalid qzone:", gstore%qzone))
 end select

 select case (gstore%kfilter)
 case ("none")
   continue

 !case ("fs_erange")
 !case ("erange") ! referred to CBM, VBM in semiconductors.

 case ("fs_tetra")
   ! This is just to show how to use the tetrahedron method
   ! to filter k- and k+q points on the FS in metals and define bands_spin automatically
   call wrtout(std_out, sjoin("Filtering k-points using:", gstore%kfilter))
   !
   call ebands_get_bands_e0(ebands, ebands%fermie, fs_bands_spin, ierr)
   ABI_CHECK(ierr == 0, "Error in ebands_get_bands_e0")
   !print *, "after ebands_get_bands_e0 with fs_bands_spin:", fs_bands_spin

   bands_spin = fs_bands_spin
   rlatt = ebands%kptrlatt; call matr3inv(rlatt, klatt)

   ABI_MALLOC(indkk, (gstore%nkbz))
   indkk(:) = kbz2ibz_(1, :) ! TODO: Decide whether it makes sense to store ktetra or indkk in gstore.

   call htetra_init(ktetra, indkk, cryst%gprimd, klatt, kbz_, gstore%nkbz, gstore%kibz, gstore%nkibz, ierr, errorstring, comm)
   ABI_CHECK(ierr == 0, errorstring)

   !np = 5; estep = one / Ha_meV
   !elow  = ebands%fermie - np * estep
   !ehigh = ebands%fermie + np * estep
   !nw = 2 * np + 1
   !nw = 1
   !ABI_MALLOC(wvals, (nw))
   !wvals = arth(elow, estep, nw)
   !wvals = [ebands%fermie]

   ABI_MALLOC(eig_ibz, (gstore%nkibz))
   max_occ = two / (ebands%nspinor * ebands%nsppol)
   select_kbz_spin = 0

   cnt = 0
   do spin=1,nsppol
     do band=fs_bands_spin(1, spin), fs_bands_spin(2, spin)
       eig_ibz = ebands%eig(band, :, spin)
       do ik_ibz=1,gstore%nkibz
         cnt = cnt + 1; if (mod(cnt, all_nproc) /= my_rank) cycle ! MPI parallelism inside comm

         call ktetra%get_onewk_wvals(ik_ibz, tetra_opt0, 1, [ebands%fermie], max_occ, &
                                     gstore%nkibz, eig_ibz, delta_theta_ef)
         iflag = merge(1, 0, abs(delta_theta_ef(1)) > zero)

         ! Use iflag to filter k-points.
         select case (gstore%kzone)
         case ("ibz")
           ik_bz = kibz2bz(ik_ibz)
           select_kbz_spin(ik_bz, spin) = iflag

         case ("bz")
           call get_star_from_ibz(ik_ibz, gstore%nkbz, kbz2ibz_, nk_in_star, kstar_bz_inds)
           ABI_CHECK(nk_in_star > 0, "Something wrong in get_star_from_ibz")
           do ii=1,nk_in_star
             ik_bz = kstar_bz_inds(ii)
             select_kbz_spin(ik_bz, spin) = iflag
           end do
           ABI_FREE(kstar_bz_inds)
         end select
       end do

     end do
   end do

   call xmpi_sum(select_kbz_spin, comm, ierr)

   ! Now the tricky part as we want to remove q-points that
   ! do not lead to any scattering process between two states on the FS
   ! Remember that k+q is a sub-mesh of the ebands k-mesh.

   select_qbz_spin = 0

   if (gstore%kzone == "ibz") kpts_ptr => kibz_
   if (gstore%kzone == "bz")  kpts_ptr => kbz_
   len_kpts_ptr = size(kpts_ptr, dim=2)
   ABI_MALLOC(map_kq, (6, len_kpts_ptr))

   select case (gstore%qzone)
   case ("ibz")
     do iq_ibz=1,gstore%nqibz
       if (mod(iq_ibz, all_nproc) /= my_rank) cycle ! MPI parallelism.
       qpt = gstore%qibz(:, iq_ibz)
       iq_bz = qibz2bz(iq_ibz)
       ! k + q_ibz --> k IBZ --> k BZ

       call gstore%krank%get_mapping(len_kpts_ptr, kpts_ptr, &
                                     dksqmax, cryst%gmet, map_kq, cryst%nsym, cryst%symafm, cryst%symrel, &
                                     ebands_timrev, use_symrec=.False., qpt=qpt)

       if (dksqmax > tol12) then
          ABI_ERROR("Cannot map k+q to IBZ!")
       end if

       do ii=1,len_kpts_ptr
         ikq_ibz = map_kq(1, ii)
         ikq_bz = kibz2bz(ikq_ibz)
         select_qbz_spin(iq_bz, :) = select_kbz_spin(ikq_bz, :)
       end do

     end do

   case ("bz")
     do iq_bz=1,gstore%nqbz
       if (mod(iq_bz, all_nproc) /= my_rank) cycle ! MPI parallelism.
       qpt = qbz_(:, iq_bz)
       iq_ibz = qbz2ibz_(1, iq_bz)
       ! k + q_bz --> k IBZ --> k BZ

       call gstore%krank%get_mapping(len_kpts_ptr, kpts_ptr, &
                                     dksqmax, cryst%gmet, map_kq, cryst%nsym, cryst%symafm, cryst%symrel, &
                                     ebands_timrev, use_symrec=.False., qpt=qpt)

       if (dksqmax > tol12) then
          ABI_ERROR("Cannot map k+q to IBZ!")
       end if

       do ii=1,len_kpts_ptr
         ikq_ibz = map_kq(1, ii)
         ikq_bz = kibz2bz(ikq_ibz)
         select_qbz_spin(iq_bz, :) = select_kbz_spin(ikq_bz, :)
       end do
     end do
   end select

   call xmpi_sum(select_qbz_spin, comm, ierr)

   ABI_FREE(map_kq)
   ABI_FREE(eig_ibz)
   ABI_FREE(indkk)
   call ktetra%free()

 case default
   ABI_ERROR(sjoin("Invalid kfilter:", gstore%kfilter))
 end select

 ! Total number of k/q points for each spin after filtering (if any)
 glob_nk_spin(:) = count(select_kbz_spin /= 0, dim=1)
 glob_nq_spin(:) = count(select_qbz_spin /= 0, dim=1)

 if (my_rank == master) then
   unts = [std_out, ab_out]
   call wrtout(unts, "=== Gstore parameters ===")
   call wrtout(unts, sjoin(" kzone:", gstore%kzone))
   call wrtout(unts, sjoin(" qzone:", gstore%qzone))
   call wrtout(unts, sjoin(" kfilter:", gstore%kfilter))
   call wrtout(unts, sjoin(" with_vk:", itoa(gstore%with_vk)))
   call wrtout(unts, sjoin(" nqibz, nkibz:", itoa(gstore%nqibz), itoa(gstore%nkibz)))
   call wrtout(unts, sjoin(" nqbz, nkbz:", itoa(gstore%nqbz), itoa(gstore%nkbz)))
   call wrtout(unts, sjoin(" glob_nk_spin:", ltoa(glob_nk_spin)))
   call wrtout(unts, sjoin(" glob_nq_spin:", ltoa(glob_nq_spin)))
 end if

 ! We need another table mapping the global index in the gqk matrix to the q/k index in the BZ
 ! so that one can extract the symmetry tables computed above.
 ! Again this is needed as the global sizes of the gqk matrix
 ! is not necessarily equal to the size of the BZ/IBZ if we have filtered the wave vectors.

 ABI_ICALLOC(qglob2bz_idx, (maxval(glob_nq_spin), nsppol))
 ABI_ICALLOC(kglob2bz_idx, (maxval(glob_nk_spin), nsppol))

 do spin=1,nsppol
   cnt = 0
   do iq_bz=1,gstore%nqbz
     if (select_qbz_spin(iq_bz, spin) /= 0) then
       cnt = cnt + 1; qglob2bz_idx(cnt, spin) = iq_bz
     end if
   end do

   cnt = 0
   do ik_bz=1,gstore%nkbz
     if (select_kbz_spin(ik_bz, spin) /= 0) then
       cnt = cnt + 1; kglob2bz_idx(cnt, spin) = ik_bz
     end if
   end do
 end do

 ! =============================================
 ! Initialize gqk basic dimensions and set flags
 ! =============================================

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   gqk%spin = spin
   gqk%natom3 = natom3
   gqk%cplex = 1 ! or 2 depending on dtset
   !gqk%cplex = dtset%gstore_cplex; if (present(gstore_cplex)) gqk%cplex = gstore_cplex
   ABI_CHECK_IRANGE(gqk%cplex, 1, 2, "gstore_cplex")

   ! The global shape of the q/k matrix for this spin.
   gqk%glob_nq = glob_nq_spin(spin)
   gqk%glob_nk = glob_nk_spin(spin)
 end do

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   ! Init for sequential execution.
   gqk%qpt_comm%nproc = 1
   gqk%kpt_comm%nproc = 1
   gqk%pert_comm%nproc = 1
   gqk%my_npert = natom3

   if (any(dtset%eph_np_pqbks /= 0)) then
     ! Use parameters from input file. Need to perform sanity check though.
     gqk%pert_comm%nproc = dtset%eph_np_pqbks(1)
     gqk%qpt_comm%nproc  = dtset%eph_np_pqbks(2)
     ABI_CHECK(dtset%eph_np_pqbks(3) == 1, "Band parallelism not implemented in gstore")
     gqk%kpt_comm%nproc = dtset%eph_np_pqbks(4)
     !gqk%spin_comm%nproc = dtset%eph_np_pqbks(5)
     gqk%my_npert = natom3 / gqk%pert_comm%nproc
     ABI_CHECK(gqk%my_npert > 0, "pert_comm_nproc cannot be greater than 3 * natom.")
     ABI_CHECK(mod(natom3, gqk%pert_comm%nproc) == 0, "pert_comm_nproc must divide 3 * natom.")

   else
     ! Automatic grid generation
     ! TODO: Add smart logic to find "optimal" distribution.
     ! Keep in mind that in gstore_build, the first loop is over q-points
     ! in order to reduce the number of interpolations of the DFPT potentials in q-space.
     ! So the q-point parallelism is expected to be more efficient.

     np = nproc_spin(spin)
     gqk%qpt_comm%nproc = np
     gqk%kpt_comm%nproc = 1
     gqk%pert_comm%nproc = 1

     ! Handle parallelism over perturbations first.
     ! Use MPI communicator to distribute the 3 * natom perturbations to reduce memory requirements for DFPT potentials.
     ! Ideally, perturbations are equally distributed --> total number of CPUs should be divisible by 3 * natom.
     ! or at least, divisible by one integer i for i in [2, 3 * natom - 1].

     ! Try to have 3 perts per proc first because the q-point parallelism is more efficient.
     ! The memory for W(R,r,ipert) will increase though.
     !do cnt=natom,2,-1
     !  if (mod(all_nproc, cnt) == 0 .and. mod(natom3, cnt) == 0) then
     !    pert_comm%nproc = cnt; new%my_npert = natom3 / cnt; exit
     !  end if
     !end do

     !if (pert_comm%nproc == 1) then
     !  ! Try again with more procs.
     !  do cnt=natom3,2,-1
     !    if (mod(all_nproc, cnt) == 0 .and. mod(natom3, cnt) == 0) then
     !      pert_comm%nproc = cnt; new%my_npert = natom3 / cnt; exit
     !    end if
     !  end do
     !end if

     !if (new%my_npert == natom3 .and. all_nproc > 1) then
     !  ABI_WARNING("The number of MPI procs should be divisible by 3*natom to reduce memory requirements!")
     !end if
   end if

   ! Distribute perturbations TODO
   ABI_MALLOC(gqk%my_iperts, (gqk%my_npert))
   gqk%my_iperts = [(ii, ii=1, natom3)]

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

#ifdef HAVE_MPI
 ! For each spin treated by this rank, create 3d cartesian communicator (q-points, k-points, perturbations)
 periods(:) = .False.; reorder = .False.

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   dims = [gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, gqk%pert_comm%nproc]

   ! Note comm_spin(spin)
   gqk%grid_comm = xcomm_from_mpi_int(comm_spin(spin))
   !call xmpi_comm_free(comm_spin(spin))

   call MPI_CART_CREATE(gqk%grid_comm, ndims, dims, periods, reorder, comm_cart, ierr)

   ! Find the index and coordinates of the current processor
   call MPI_COMM_RANK(comm_cart, me_cart, ierr)
   call MPI_CART_COORDS(comm_cart, me_cart, ndims, gqk%coords_qkp, ierr)

   ! Create communicator for q-points
   keepdim = .False.; keepdim(1) = .True.
   call MPI_CART_SUB(comm_cart, keepdim, gqk%qpt_comm%value, ierr); gqk%qpt_comm%me = xmpi_comm_rank(gqk%qpt_comm%value)

   ! Create communicator for k-points
   keepdim = .False.; keepdim(2) = .True.
   call MPI_CART_SUB(comm_cart, keepdim, gqk%kpt_comm%value, ierr); gqk%kpt_comm%me = xmpi_comm_rank(gqk%kpt_comm%value)

   ! Create communicator for perturbations.
   keepdim = .False.; keepdim(3) = .True.
   call MPI_CART_SUB(comm_cart, keepdim, gqk%pert_comm%value, ierr); gqk%pert_comm%me = xmpi_comm_rank(gqk%pert_comm%value)

   call xmpi_comm_free(comm_cart)

   if (my_rank == master) then
     write(std_out, "(/,a)")" === Gstore MPI distribution ==="
     write(std_out, "(a,i0)")"P Number of CPUs for parallelism over perturbations: ", gqk%pert_comm%nproc
     write(std_out, "(a,i0)")"P Number of perturbations treated by this CPU: ", gqk%my_npert
     write(std_out, "(a,i0)")"P Number of CPUs for parallelism over q-points: ", gqk%qpt_comm%nproc
     write(std_out, "(a,i0)")"P Number of CPUs for parallelism over k-points: ", gqk%kpt_comm%nproc
   end if
 end do ! my_is
#endif

 ! At this point, we have the Cartesian grid (one per spin if any)
 ! and we can finally distribute other arrays.

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   ! Compute bstart and band size for this spin.
   gqk%bstart = bands_spin(1, spin)
   gqk%nb = bands_spin(2, spin) - bands_spin(1, spin) + 1
   !print *, "gqk%bstart:", gqk%bstart; print *, "gqk%nb:", gqk%nb

   ! Split q-points and transfer symmetry tables.
   ! FIXME: Note that glob_nq and glob_nk does not necessarily correspond to the size of the BZ
   ! First of all we have to consider kzone
   ! Even if kzone == "bz" we may have filtered the wavevectors e.g. Fermi surface.
   call xmpi_split_block(gqk%glob_nq, gqk%kpt_comm%value, gqk%my_nq, myq2glob)
   ABI_CHECK(gqk%my_nq > 0, "my_nq == 0")
   gqk%my_qstart = myq2glob(1)
   ABI_FREE(myq2glob)

   ABI_MALLOC(gqk%my_q2ibz, (6, gqk%my_nq))
   do my_iq=1,gqk%my_nq
     iq_glob = my_iq + gqk%my_qstart - 1
     iq_bz = qglob2bz_idx(iq_glob, spin)
     gqk%my_q2ibz(:, my_iq) = qbz2ibz_(:, iq_bz)
   end do

   ! Split k-points and transfer symmetry tables
   call xmpi_split_block(gqk%glob_nk, gqk%kpt_comm%value, gqk%my_nk, myk2glob)
   ABI_CHECK(gqk%my_nk > 0, "my_nk == 0")
   gqk%my_kstart = myk2glob(1)
   ABI_FREE(myk2glob)

   ABI_MALLOC(gqk%my_k2ibz, (6, gqk%my_nk))
   do my_ik=1,gqk%my_nk
     ik_glob = my_ik + gqk%my_kstart - 1
     ik_bz = kglob2bz_idx(ik_glob, spin)
     gqk%my_k2ibz(:, my_ik) = kbz2ibz_(:, ik_bz)
   end do

   mem_mb = gqk%cplex * gqk%my_npert * gqk%nb ** 2 * gqk%my_nq * gqk%my_nk * eight * b2Mb
   call wrtout(std_out, sjoin(" Local memory for e-ph matrix elements:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))

   ! Allocate storage for MPI-distributed e-ph matrix elements.
   !ABI_MALLOC_OR_DIE(gqk%my_vals, (gqk%cplex, gqk%nb, gqk%my_nq, gqk%nb, gqk%my_nk, gqk%my_npert), ierr)
   ABI_MALLOC_OR_DIE(gqk%my_vals, (gqk%cplex, gqk%my_npert, gqk%nb, gqk%my_nq, gqk%nb, gqk%my_nk), ierr)
   gqk%my_vals = zero !; gqk%my_vals = huge(one)

   ! Allocate storage for MPI-distributed dH/dk matrix elements.
   select case (gstore%with_vk)
   case (1)
     mem_mb = 3 * gqk%nb * gqk%my_nk * eight * b2Mb
     call wrtout(std_out, sjoin(" Local memory for diagonal vk:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))
     ABI_MALLOC(gqk%my_vk_cart, (3, gqk%nb, gqk%my_nk))
   case (2)
     mem_mb = 3 * gqk%nb ** 2 * gqk%my_nk * eight * b2Mb
     call wrtout(std_out, sjoin(" Local memory for vkmat:", ftoa(mem_mb, fmt="f8.1"), " [Mb] <<< MEM"))
     ABI_MALLOC(gqk%my_vkmat_cart, (3, gqk%nb, gqk%nb, gqk%my_nk))
   end select

 end do ! my_is

 ABI_FREE(qbz_)
 ABI_FREE(qbz2ibz_)

 ! TODO: Use ebands%kptns
 ABI_FREE(wtk_)
 ABI_FREE(kbz_)
 ABI_FREE(kbz2ibz_)

 ABI_FREE(qibz2bz)
 ABI_FREE(kibz2bz)

 ABI_FREE(select_qbz_spin)
 ABI_FREE(select_kbz_spin)
 ABI_FREE(qglob2bz_idx)
 ABI_FREE(kglob2bz_idx)

 call cwtime_report(" gstore_new:", cpu, wall, gflops)

end function gstore_new
!!***

subroutine get_ibz2bz(nibz, nbz, bz2ibz, ibz2bz, ierr)

 integer,intent(in) :: nibz, nbz
 integer,intent(in) :: bz2ibz(6, nbz)
 integer,intent(out) :: ierr
 integer,allocatable,intent(out) :: ibz2bz(:)

!Local variables-------------------------------
!scalars
 integer :: iq_bz, iq_ibz, isym_q, trev_q, cnt, g0_q(3)
 logical :: isirr_q

!----------------------------------------------------------------------

 ABI_MALLOC(ibz2bz, (nibz))

 cnt = 0
 do iq_bz=1,nbz
   iq_ibz = bz2ibz(1, iq_bz); isym_q = bz2ibz(2, iq_bz)
   trev_q = bz2ibz(6, iq_bz); g0_q = bz2ibz(3:5,iq_bz)
   isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
   if (isirr_q) then
     cnt = cnt + 1
     ibz2bz(iq_ibz) = iq_bz
   end if
 end do

 ierr = 0
 if (cnt /= nibz) ierr = 1

end subroutine get_ibz2bz
!!***

subroutine get_star_from_ibz(ik_ibz, nkbz, bz2ibz, nk_in_star, kstar_bz_inds)

 integer,intent(in) :: ik_ibz, nkbz
 integer,intent(out) :: nk_in_star
 integer,intent(in) :: bz2ibz(6, nkbz)
 integer,allocatable,intent(out) :: kstar_bz_inds(:)

!Local variables-------------------------------
!scalars
 integer :: iq_bz

!----------------------------------------------------------------------

 nk_in_star = count(bz2ibz(1, :) == ik_ibz)

 ABI_MALLOC(kstar_bz_inds, (nk_in_star))

 nk_in_star = 0
 do iq_bz=1,nkbz
   if (bz2ibz(1, iq_bz) /= ik_ibz) continue
   nk_in_star = nk_in_star + 1
   kstar_bz_inds(nk_in_star) = iq_bz
 end do

end subroutine get_star_from_ibz
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_spin2my_is
!! NAME
!! gstore_spin2my_is
!!
!! FUNCTION
!!  Return the local spin index from the global spin index. 0 if not treated.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gstore_fill_bks_mask(gstore, mband, nkibz, nsppol, bks_mask)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 integer,intent(in) :: mband, nkibz, nsppol
 logical,intent(out) :: bks_mask(mband, nkibz, nsppol)

!Local variables-------------------------------
!scalars
 integer :: my_is, my_ik, my_iq, spin, ik_ibz, iqk_ibz
 real(dp) :: dksqmax
 type(gqk_t),pointer :: gqk
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: indkk_kq(:,:)
 real(dp) :: qpt(3)
 real(dp),allocatable :: my_kpts(:,:)

!----------------------------------------------------------------------

 bks_mask = .False.
 cryst => gstore%cryst

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   ! We need the image of this k in the IBZ.
   do my_ik=1,gqk%my_nk
     ik_ibz = gqk%my_k2ibz(1, my_ik)
     bks_mask(gqk%bstart:gqk%bstart + gqk%nb - 1, ik_ibz, spin) = .True.
   end do

   ! as well as the image of k+q in the IBZ.
   call gqk%get_all_mykpts(gstore, my_kpts)
   ABI_MALLOC(indkk_kq, (6, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     qpt = gqk%get_myqpt(my_iq, gstore)

     call gstore%krank%get_mapping(gqk%my_nk, my_kpts, dksqmax, cryst%gmet, indkk_kq, &
                                   cryst%nsym, cryst%symafm, cryst%symrel, kpts_timrev_from_kptopt(gstore%ebands%kptopt), &
                                   use_symrec=.False., qpt=qpt)

     if (dksqmax > tol12) then
       ABI_ERROR(sjoin("Cannot map k+q to IBZ with q:", ktoa(qpt), ", dkqsmax:", ftoa(dksqmax)))
     end if

     do my_ik=1,gqk%my_nk
       iqk_ibz = indkk_kq(1, my_ik)
       bks_mask(gqk%bstart:gqk%bstart + gqk%nb - 1, iqk_ibz, spin) = .True.
     end do
   end do

   ABI_FREE(my_kpts)
   ABI_FREE(indkk_kq)
 end do

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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gstore_get_mpw_gmax(gstore, ecut, mpw, gmax)

!Arguments ------------------------------------
 class(gstore_t),target,intent(in) :: gstore
 real(dp),intent(in) :: ecut
 integer,intent(out) :: mpw, gmax(3)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, spin, onpw, ii, ipw, ierr, my_mpw
 type(gqk_t),pointer :: gqk
!arrays
 integer :: my_gmax(3)
 integer,allocatable :: gtmp(:,:)
 real(dp) :: kk(3), qpt(3)

!----------------------------------------------------------------------

 mpw = 0; gmax = 0

 ! FIXME: This is an hotspot due to the double loop over k and q.
 ! Should use a geometrical approach to compute mpw and gmax.

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   do my_ik=1,gqk%my_nk
     kk = gqk%get_mykpt(my_ik, gstore)

     ! Compute G sphere, returning npw. Note istwfk == 1.
     call get_kg(kk, 1, ecut, gstore%cryst%gmet, onpw, gtmp)
     mpw = max(mpw, onpw)
     do ipw=1,onpw
       do ii=1,3
         gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
       end do
     end do
     ABI_FREE(gtmp)

     do my_iq=1,gqk%my_nq
       qpt = gqk%get_myqpt(my_iq, gstore)

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
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gstore_calc_my_phonons(gstore, store_phdispl)

!Arguments ------------------------------------
 class(gstore_t),target,intent(inout) :: gstore
 logical,intent(in) :: store_phdispl

!Local variables-------------------------------
 integer :: my_is, my_iq, ierr
 real(dp) :: qpt(3)
 type(crystal_t),pointer :: cryst
 type(gqk_t),pointer :: gqk
 real(dp),allocatable :: phfrq(:), displ_cart(:,:,:,:) !, displ_red(:,:,:,:)

!----------------------------------------------------------------------

 cryst => gstore%cryst

 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, 3 * cryst%natom))
 ABI_MALLOC(phfrq, (3 * cryst%natom))

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)

   ! Get phonon frequencies and eigenvectors for this q-point.
   ABI_CALLOC(gqk%my_wqnu, (gqk%my_npert, gqk%my_nq))
   if (store_phdispl) then
     ABI_CALLOC(gqk%my_displ_cart, (2, 3, cryst%natom, gqk%my_npert, gqk%my_nq))
   end if

   do my_iq=1,gqk%my_nq
     if (gqk%kpt_comm%skip(my_iq)) cycle ! MPI parallelism inside kpt_comm
     qpt = gqk%get_myqpt(my_iq, gstore)
     call gstore%ifc%fourq(cryst, qpt, phfrq, displ_cart) !, out_displ_red=displ_red)

     gqk%my_wqnu(:, my_iq) = phfrq
     gqk%my_wqnu(:, my_iq) = phfrq(gqk%my_iperts(:))
     if (store_phdispl) gqk%my_displ_cart(:,:,:,:,my_iq) = displ_cart(:,:,:,gqk%my_iperts(:))
   end do

   call xmpi_sum(gqk%my_wqnu, gqk%kpt_comm%value, ierr)
   if (store_phdispl) call xmpi_sum(gqk%my_displ_cart, gqk%kpt_comm%value, ierr)
   ! Note that the computation is replicated across the perturbation communicator.
   ! This should not represent a problem since the phase in the displacement is fixed in dfpt_phfrq
   ! hence results should be coherent across different procs.
 end do

 ABI_FREE(displ_cart)
 ABI_FREE(phfrq)

end subroutine gstore_calc_my_phonons
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_free
!! NAME
!! gstore_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
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
 ABI_SFREE(gstore%qibz)

 ABI_SFREE(gstore%my_spins)
 call gstore%spin_comm%free()

 call gstore%krank%free()

end subroutine gstore_free
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_get_mykpt
!! NAME
!! gqk_get_mykpt
!!
!! FUNCTION
!!  Return the reduced coordinates of the k-point from the local index my_ik
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function gqk_get_mykpt(gqk, my_ik, gstore) result (kpt)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 class(gstore_t),intent(in) :: gstore
 integer,intent(in) :: my_ik
 real(dp) :: kpt(3)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz, isym_k, trev_k, tsign, g0_k(3)
 logical :: isirr_k

!----------------------------------------------------------------------

 ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
 trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5, my_ik)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 tsign = 1; if (trev_k == 1) tsign = -1

 ! symrel^T convention for k
 kpt = tsign * matmul(transpose(gstore%cryst%symrel(:,:,isym_k)), gstore%kibz(:, ik_ibz)) + g0_k

end function gqk_get_mykpt
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_get_all_mykpts
!! NAME
!! gqk_get_all_mykpts
!!
!! FUNCTION
!!  Allocate and return array with the all reduced coordinates of the k-point
!!  treated by this MPI rank
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gqk_get_all_mykpts(gqk, gstore, my_kpts)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 class(gstore_t),intent(in) :: gstore
 real(dp),allocatable,intent(out) :: my_kpts(:,:)

!Local variables ------------------------------
!scalars
 integer :: my_ik

!----------------------------------------------------------------------

 ABI_MALLOC(my_kpts, (3, gqk%my_nk))
 do my_ik=1,gqk%my_nk
   my_kpts(:, my_ik) = gqk%get_mykpt(my_ik, gstore)
 end do

end subroutine gqk_get_all_mykpts
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_get_myqpt
!! NAME
!! gqk_get_myqpt
!!
!! FUNCTION
!!  Return the reduced coordinates of the q-point from the local index my_iq
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function gqk_get_myqpt(gqk, my_iq, gstore) result (qpt)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 class(gstore_t),intent(in) :: gstore
 integer,intent(in) :: my_iq
 real(dp) :: qpt(3)

!Local variables ------------------------------
!scalars
 integer :: iq_ibz, isym_q, trev_q, tsign, g0_q(3)
 logical :: isirr_q

!----------------------------------------------------------------------

 iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
 trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
 isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
 tsign = 1; if (trev_q == 1) tsign = -1

 ! symrec convention for q
 !qpt = tsign * matmul(transpose(gqk%cryst%symrel(:,:,isym_q)), gqk%qibz(:, iq_ibz)) + g0_q
 qpt = tsign * matmul(gstore%cryst%symrec(:,:,isym_q), gstore%qibz(:, iq_ibz)) + g0_q

end function gqk_get_myqpt
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_free
!! NAME
!! gqk_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gqk_free(gqk)

!Arguments ------------------------------------
 class(gqk_t),intent(inout) :: gqk

!----------------------------------------------------------------------

 !ABI_SFREE(gqk%kibz)
 ABI_SFREE(gqk%my_k2ibz)
 ABI_SFREE(gqk%my_q2ibz)

 ABI_SFREE(gqk%my_wqnu)
 ABI_SFREE(gqk%my_displ_cart)
 ABI_SFREE(gqk%my_vals)
 ABI_SFREE(gqk%my_iperts)
 ABI_SFREE(gqk%my_vk_cart)
 ABI_SFREE(gqk%my_vkmat_cart)

 ! Communicators
 call gqk%kpt_comm%free()
 call gqk%qpt_comm%free()
 call gqk%pert_comm%free()
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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gstore_compute(gstore, wfk0_path, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, &
                          pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 !type(datafiles_type),intent(in) :: dtfil
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
 integer,parameter :: tim_getgh1c = 1, berryopt0 = 0, ider0 = 0, idir0 = 0
 integer,parameter :: useylmgr = 0, useylmgr1 = 0, master = 0, ndat1 = 1
 integer :: my_rank,nproc,mband,nsppol,nkibz,idir,ipert,ebands_timrev !iq_ibz,
 integer :: cplex,natom,natom3,ipc,nspinor
 integer :: bstart_k,bstart_kq,nband_k,nband_kq,band_k, ib_k, ib_kq !ib1,ib2, band_kq,
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k, trev_kq
 integer :: my_ik, my_is, comm_rpt, my_npert, my_ip, my_iq
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: mpw,mnb,ierr !,cnt ii,jj,
 integer :: n1,n2,n3,n4,n5,n6,nspden
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf, nkpg,nkpg1
 real(dp) :: cpu, wall, gflops, cpu_q, wall_q, gflops_q, cpu_k, wall_k, gflops_k, cpu_all, wall_all, gflops_all
 real(dp) :: ecut, eshift, eig0nk, dksqmax
 logical :: gen_eigenpb, isirr_k, isirr_kq
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(ddkop_t) :: ddkop
 type(gqk_t),pointer :: gqk
 character(len=500) :: msg
!arrays
 integer :: g0_k(3), g0_kq(3), work_ngfft(18),gmax(3),gamma_ngqpt(3), indkk_kq(6,1)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),nband(:,:),wfd_istwfk(:)
 !integer,allocatable :: my_pinfo(:,:) !, pert_table(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3), vk(3) !, vkq(3)
 real(dp) :: phfrq(3*cryst%natom)
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: displ_cart(:,:,:,:), displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:), kinpw1(:), kpg1_k(:,:), kpg_k(:,:), dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:), ffnl1(:,:,:,:), ph3d(:,:,:), ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:), gkk_atm(:,:,:,:),gkq_nu(:,:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:), kets_k(:,:,:), h1kets_kq(:,:,:) !, cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:), vlocal(:,:,:,:), vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:), ylm_k(:,:), ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:), gvnlx1(:,:), work(:,:,:,:)
 real(dp),allocatable :: gs1c(:,:), vcar_ibz(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)

!************************************************************************

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

 ! Define q-mesh. eph_ngqpt_fine activates the Fourier interpolation of the DFPT potentials.
 gamma_ngqpt = ifc%ngqpt; if (all(dtset%eph_ngqpt_fine /= 0)) gamma_ngqpt = dtset%eph_ngqpt_fine

 !call wrtout(std_out, sjoin("q-mesh for the phonon linewidths:", ltoa(gamma_ngqpt)))
 !call wrtout(std_out, sjoin("Will compute", itoa(gams%nqibz), "q-points in the IBZ"))

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   if (gqk%pert_comm%nproc > 1) then
     ! Activate parallelism over perturbations
     ! Build table with list of perturbations treated by this CPU inside pert_comm
     ABI_ERROR("pert_comm%nproc > 1 not tested")
     !call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, gqk%pert_comm%value)
     !call dvdb%set_pert_distrib(my_npert, natom3, my_pinfo, pert_table, gqk%pert_comm%value)
     !ABI_FREE(my_pinfo)
     !ABI_FREE(pert_table)
   end if
 end do

 !call wrtout([std_out, ab_out], " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
 ! Prepare Fourier interpolation of DFPT potentials.
 comm_rpt = xmpi_comm_self
 !comm_rpt = bqs_comm%value
 call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k/k+q wavevectors
 ! treated by this MPI rank are stored.

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
 call cwtime(cpu, wall, gflops, "start")
 call wrtout(std_out, " Computing mpw. This may take some time for dense k/q meshes...")

 call gstore%get_mpw_gmax(ecut, mpw, gmax)

 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))
 call cwtime_report(" gmax and mpw", cpu, wall, gflops)

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
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, use_gpu_cuda=dtset%use_gpu_cuda)

 ! Allocate vlocal. Note nvloc
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red, (2, 3, cryst%natom, natom3))
 ABI_CALLOC(dummy_vtrial, (nfftf, nspden))

 ! Create ddkop object to compute group velocities (if needed)
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ! Loop over my spins.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   my_npert = gqk%my_npert

   ! Allocate workspace for wavefunctions. Make npw larger than expected.
   ! maxnb is the maximum number of bands crossing the FS, used to dimension arrays.
   !mnb = fs%maxnb
   mnb = gqk%nb
   ABI_MALLOC(bras_kq, (2, mpw*nspinor, mnb))
   ABI_MALLOC(kets_k, (2, mpw*nspinor, mnb))
   ABI_MALLOC(h1kets_kq, (2, mpw*nspinor, mnb))
   ABI_MALLOC(gkk_atm, (2, mnb, mnb, natom3))
   ABI_MALLOC(gkq_nu, (2, mnb, mnb, natom3))

   ! Loop over my set of q-points
   do my_iq=1,gqk%my_nq
     call cwtime(cpu_q, wall_q, gflops_q, "start")

     qpt = gqk%get_myqpt(my_iq, gstore)
     !iq_ibz = gqk%my_q2ibz(1, my_iq)

     ! Use Fourier interpolation of DFPT potentials to get my_npert potentials.
     cplex = 2
     ABI_MALLOC(v1scf, (cplex, nfft, nspden, dvdb%my_npert))
     call dvdb%ftinterp_qpt(qpt, nfftf, ngfftf, v1scf, dvdb%comm_rpt)

     ! Get phonon frequencies and eigenvectors for this q-point.
     call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

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
       call cwtime(cpu_k, wall_k, gflops_k, "start")

       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       kk = gqk%get_mykpt(my_ik, gstore)
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
       kq = kk + qpt

       call gstore%krank%get_mapping(1, kq, dksqmax, cryst%gmet, indkk_kq, cryst%nsym, cryst%symafm, cryst%symrel, &
                                    ebands_timrev, use_symrec=.False.)

       if (dksqmax > tol12) then
         write(msg, '(3a,es16.6,6a)' ) &
          "The WFK file cannot be used to compute phonon linewidths.",ch10, &
          "At least one of the k-points on the FS could not be generated from a symmetrical one. dksqmax: ", dksqmax,ch10, &
          "Q-mesh: ", trim(ltoa(gamma_ngqpt)), ", K-mesh (from kptrlatt) ", trim(ltoa(get_diag(ebands%kptrlatt))), &
          'Action: check your WFK file and the (k, q) point input variables.'
          ABI_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1, 1); isym_kq = indkk_kq(2, 1)
       trev_kq = indkk_kq(6, 1); g0_kq = indkk_kq(3:5, 1)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)

       ! If we have used the KERANGE trick, we may have k or k+q points with just one G component set to zero
       ! so we skip this transition immediately. This should happen only if fsewin > sigma_erange.
       if (wfd%npwarr(ik_ibz) == 1 .or. wfd%npwarr(ikq_ibz) == 1) cycle

       ! Number of bands crossing the Fermi level at k+q
       !bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz)
       bstart_kq = gqk%bstart; nband_kq = gqk%nb

       ABI_CHECK(nband_k <= mnb .and. nband_kq <= mnb, "wrong nband")

       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kk, kk_ibz, bstart_k, nband_k, spin, mpw, gqk%my_k2ibz(:, my_ik), cryst, &
                          work_ngfft, work, istwf_k, npw_k, kg_k, kets_k)

       ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kq, kq_ibz, bstart_kq, nband_kq, spin, mpw, indkk_kq(:,1), cryst, &
                          work_ngfft, work, istwf_kq, npw_kq, kg_kq, bras_kq)

       ! if PAW, one has to solve a generalized eigenproblem
       ! Be careful here because I will need sij_opt==-1
       gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2, npw_kq*nspinor*((sij_opt+1)/2)))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm == 1) then
       !   optder = 0; if (useylmgr == 1) optder = 1
       !   call initylmg(cryst%gprimd, kg_k, kk, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_k], dtset%nsppol, optder, cryst%rprimd, ylm_k, ylmgr)
       !   call initylmg(cryst%gprimd, kg_kq, kq, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_kq], dtset%nsppol, optder, cryst%rprimd, ylm_kq, ylmgr_kq)
       !end if

       ! Compute k+G vectors
       nkpg = 3 * dtset%nloalg(3)
       ABI_MALLOC(kpg_k, (npw_k, nkpg))
       if (nkpg > 0) call mkkpg(kg_k, kpg_k, kk, nkpg, npw_k)

       ! Compute nonlocal form factors ffnlk at (k+G)
       ABI_MALLOC(ffnlk, (npw_k, 1, psps%lmnmax, psps%ntypat))
       !call mkffnl(psps%dimekb, 1, psps%ekb, ffnlk, psps%ffspl,&
       !            cryst%gmet, cryst%gprimd, ider0, idir0, psps%indlmn, kg_k, kpg_k, kk, psps%lmnmax, &
       !            psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, psps%ntypat, &
       !            psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum, &
       !            comm=gqk%pert_comm%value)

       call mkffnl_objs(cryst, psps, 1, ffnlk, ider0, idir0, kg_k, kpg_k, kk, nkpg, npw_k, ylm_k, ylmgr_dum, &
                        comm=gqk%pert_comm%value)

       ! Compute k+q+G vectors
       nkpg1 = 3 * dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G)
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       !call mkffnl(psps%dimekb, 1, psps%ekb, ffnl1, psps%ffspl, cryst%gmet, cryst%gprimd, ider0, idir0, &
       !            psps%indlmn, kg_kq, kpg1_k, kq, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg1, &
       !            npw_kq, psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_kq, ylmgr_kq, &
       !            comm=gqk%pert_comm%value)

       call mkffnl_objs(cryst, psps, 1, ffnl1, ider0, idir0, kg_kq, kpg1_k, kq, nkpg1, npw_kq, ylm_kq, ylmgr_kq, &
                        comm=gqk%pert_comm%value)

       ! Loop over my atomic perturbations and compute gkk_atm.
       gkk_atm = zero
       do my_ip=1,my_npert
         idir = dvdb%my_pinfo(1, my_ip); ipert = dvdb%my_pinfo(2, my_ip); ipc = dvdb%my_pinfo(3, my_ip)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)

         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,my_ip), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk, kq, idir, ipert, &                   ! In
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
         ! The array eig1_k contains:
         !
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
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
       call ephtk_gkknu_from_atm(mnb, mnb, 1, natom, gkk_atm, phfrq, displ_red, gkq_nu)

       ! Save e-ph matrix elements. Remember the shape of the two arrays.
       !
       !    gkq_nu(2, mnb, mnb, natom3)
       !    my_vals(cplex, my_npert, nb, my_nq, nb, my_nk)
       !
       select case (gqk%cplex)
       case (1)
        do my_ip=1,my_npert
          ipert = gqk%my_iperts(my_ip)
          gqk%my_vals(1, my_ip, :, my_iq, :, my_ik) = gkq_nu(1, :, :, ipert) ** 2 + gkq_nu(2, :, :, ipert) ** 2
         end do
       case (2)
         do my_ip=1,my_npert
           ipert = gqk%my_iperts(my_ip)
           gqk%my_vals(:, my_ip, :, my_iq, :, my_ik) = gkq_nu(:, :, :, ipert)
         end do
       case default
         ABI_ERROR(sjoin("Invalid cplex:", itoa(gqk%cplex)))
       end select

       if (gstore%with_vk /= 0) then
         ! Compute diagonal matrix elements of velocity operator with DFPT routines
         ! Velocities are in Cartesian coordinates.
         !
         ! If k+q is not in the IBZ, we need to recostruct the value by symmetry using v(Sq) = S v(q).
         ! Use transpose(R) because we are using the tables for the wavefunctions
         ! In this case listkk has been called with symrec and use_symrec=False
         ! so q_bz = S^T q_ibz where S is the isym_kq symmetry

         call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, kg_k)

         if (gstore%with_vk == 1) then
           do ib_k=1,nband_k
             band_k = ib_k + bstart_k - 1
             !vk = vcar_ibz(:, band_k, ik_ibz, spin)
             !if (.not. isirr_k) then
             !  vk = matmul(transpose(cryst%symrel_cart(:,:,isym_k)), vk)
             !  if (trev_k /= 0) vk = -vk
             !end if
             !fs%vk(:,ib_k) = vk
             vk = ddkop%get_vdiag(ebands%eig(band_k, ik_ibz, spin), &
                                  istwf_k, npw_k, wfd%nspinor, kets_k(:,:,ib_k), cwaveprj0)

             gqk%my_vk_cart(:, ib_k, my_ik) = vk
           end do

         else if (gstore%with_vk == 2) then
           ABI_ERROR("with_vk 2")
           do ib_k=1,nband_k
             band_k = ib_k + bstart_k - 1
           end do
           !gqk%my_vk_cart(3, nb, nb, my_ik) =
         end if

         !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kq, istwf_kq, npw_kq, kg_kq)
         !do ib_kq=1,nband_kq
         !  band_kq = ib_kq + bstart_kq - 1
         !  !vkq = ddkop%get_vdiag(ebands%eig(band_kq, ikq_ibz, spin), &
         !  !                                 istwf_kq, npw_kq, wfd%nspinor, bras_kq(:,:,ib_kq), cwaveprj0)
         !  !fs%vkq(:,ib_kq) = vkq
         !end do
       end if
     end do ! my_ik

     ABI_FREE(v1scf)
     ABI_FREE(vlocal1)

     if (my_iq <= 10 .or. mod(my_iq, 100) == 0) then
       write(msg,'(2(a,i0),a)')" Computation of my q-point [", my_iq, "/", gqk%my_nq, "]"
       call cwtime_report(msg, cpu_q, wall_q, gflops_q) !, end_str=ch10)
       if (my_iq == 10) call wrtout(std_out, " >>> Will print next iteration when mod(my_iq, 100) == 0) ...")
     end if
   end do ! my_iq

   ABI_FREE(bras_kq)
   ABI_FREE(kets_k)
   ABI_FREE(h1kets_kq)
   ABI_FREE(gkk_atm)
   ABI_FREE(gkq_nu)
 end do ! my_is

 call cwtime_report(" gstore_compute", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10)

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
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_SFREE(vcar_ibz)

 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 call ddkop%free()
 call gs_hamkq%free()
 call wfd%free()

end subroutine gstore_compute
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_ncwrite_path
!! NAME
!! gstore_ncwrite_path
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gstore_ncwrite_path(gstore, path, cryst, ebands)

!Arguments ------------------------------------
!scalars
 class(gstore_t),intent(in) :: gstore
 character(len=*),intent(in) :: path
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_rank, ncid, spin, my_is, ncerr

! *************************************************************************

 my_rank = xmpi_comm_rank(gstore%comm)

 if (my_rank == master) then
   ! =======================================================
   ! Master node writes basic objects and gstore dimensions
   ! =======================================================
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))

   ! Add extra gstore stuff
   ncerr = nctk_def_dims(ncid, [ &
      nctkdim_t("nkibz", gstore%nkibz), &
      nctkdim_t("nkbz", gstore%nkbz), &
      nctkdim_t("nqibz", gstore%nqibz), &
      nctkdim_t("nqbz", gstore%nqbz), &
      nctkdim_t("nctk_slen", nctk_slen) &
     ], &
   defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "with_vk"])
   NCF_CHECK(ncerr)
   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
   !NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("qibz", "dp", "three, nqibz"), &
     nctkarr_t("kibz", "dp", "three, nkibz"), &
     nctkarr_t("kzone", "c", "nctk_slen"), &
     nctkarr_t("qzone", "c", "nctk_slen"), &
     nctkarr_t("kfilter", "c", "nctk_slen") &
   ])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, vid("with_vk"), gstore%with_vk))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qibz"), gstore%qibz))
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kibz"), gstore%kibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kzone"), gstore%kzone))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qzone"), gstore%qzone))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kfilter"), gstore%kfilter))
   NCF_CHECK(nf90_close(ncid))
 end if ! master

 call xmpi_barrier(gstore%comm)

 ! Now we write the big arrays with MPI-IO and hdf5 groups.
 ! Note the loop over spin to account as each gqk has its own dimensions.
 !
 do spin=1,gstore%nsppol
   my_is = gstore%spin2my_is(spin)
   if (my_is /= 0) call gstore%gqk(my_is)%ncwrite_path(path, gstore)
   call xmpi_barrier(gstore%comm)
 end do

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine gstore_ncwrite_path
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_ncwrite_path
!! NAME
!! gqk_ncwrite_path
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gqk_ncwrite_path(gqk, path, gstore)

!Arguments ------------------------------------
!scalars
 class(gqk_t),intent(in) :: gqk
 character(len=*),intent(in) :: path
 class(gstore_t),intent(in) :: gstore

!Local variables-------------------------------
!scalars
 integer :: root_ncid, spin_ncid, ncerr

! *************************************************************************

 !my_rank = xmpi_comm_rank(gstore%comm)
 NCF_CHECK(nctk_open_modify(root_ncid, path, gqk%grid_comm%value))

 ! Create group to store gkq datastructure for this spin.
 NCF_CHECK(nf90_def_grp(root_ncid, strcat("gqk", "_spin", itoa(gqk%spin)), spin_ncid))

 ! =====================
 ! === Write dimensions
 ! =====================
 ncerr = nctk_def_dims(spin_ncid, [ &
    nctkdim_t("gqk_cplex", gqk%cplex), &
    !nctkdim_t("nkibz", gqk%cplex), &
    nctkdim_t("nb", gqk%nb), &
    nctkdim_t("natom3", gqk%natom3), &
    nctkdim_t("glob_nk", gqk%glob_nk), &
    nctkdim_t("glob_nq", gqk%glob_nq)  &
   ], &
   defmode=.True.)
 NCF_CHECK(ncerr)

 ! ==============
 ! Define scalars
 ! ==============

 ncerr = nctk_def_iscalars(spin_ncid, [character(len=nctk_slen) :: "bstart"])
 NCF_CHECK(ncerr)
 !ncerr = nctk_def_dpscalars(spin_ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
 !NCF_CHECK(ncerr)

 ! Here I will need to gather stuff along the k/q axis before writing
 ! so it makes more sense to open file using a single row/col comm in k/q grid
 ! while keeping in mind the per_comm.

 ! =============
 ! Define arrays
 ! =============
 ncerr = nctk_def_arrays(spin_ncid, [ &
   nctkarr_t("vals", "dp", "gqk_cplex, natom3, nb, glob_nq, nb, glob_nk") &  ! Global matrix
   !nctkarr_t("q2ibz", "int", "six, glob_nq"), &
   !nctkarr_t("vals", "dp", "gqk_cplex, natom3, nb, glob_nq, nb, glob_nk") &
 ])
 NCF_CHECK(ncerr)

 if (gstore%with_vk == 1) then
   ncerr = nctk_def_arrays(spin_ncid, [ &
     nctkarr_t("vk_cart", "dp", "three, nb, glob_nk") &  ! Global matrix
   ])
   NCF_CHECK(ncerr)
 else if (gstore%with_vk == 2) then
   ncerr = nctk_def_arrays(spin_ncid, [ &
     nctkarr_t("vkmat_cart", "dp", "three, nb, nb, glob_nk") &  ! Global matrix
   ])
   NCF_CHECK(ncerr)
 end if

 ! ==========
 ! Write data
 ! ==========
 NCF_CHECK(nctk_set_datamode(spin_ncid))
 NCF_CHECK(nf90_put_var(spin_ncid, vid("bstart"), gqk%bstart))

 if (gstore%with_vk == 1) then
   !NCF_CHECK(nf90_put_var(spin_ncid, vid("vk_cart"), ??))
 else if (gstore%with_vk == 2) then
   !NCF_CHECK(nf90_put_var(spin_ncid, vid("vkmat_cart"), ??))
 end if

 NCF_CHECK(nf90_close(root_ncid))

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(spin_ncid, vname)
 end function vid

end subroutine gqk_ncwrite_path
!!***

end module m_gstore
!!***

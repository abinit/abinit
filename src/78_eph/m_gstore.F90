!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
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
 use m_tetrahedron
 use m_ifc
 use m_ebands
 use m_fstab
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
 use m_fstrings,       only : toupper, itoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, wrap2_pmhalf, simpson_int, simpson, mkherm, get_diag, isdiagmat
 use m_io_tools,       only : iomode_from_fname ! open_file
 use m_fftcore,        only : ngfft_seq, get_kg
 use m_cgtools,        only : cg_zdotc
 use m_kg,             only : getph, mkkpg
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_bz_mesh,        only : kpath_t, kpath_new
 use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk, kpts_timrev_from_kptopt
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 !use m_symtk,           only : matr3inv
 !use m_numeric_tools,   only : arth, inrange, wrap2_pmhalf
 !use m_special_funcs,   only : gaussian
 !use m_fstrings,        only : strcat, ltoa, itoa, ftoa, ktoa, sjoin
 !use m_kpts,            only : kpts_timrev_from_kptopt, listkk, kpts_ibz_from_kptrlatt
 !use m_occ,             only : occ_fd, occ_be

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_gstore/gqk_t
!! NAME
!! gqk_t
!!
!! FUNCTION
!!  Stores e-p matrix elements for a given spin as several dimensions depende on
!!  isppol when nsppol == 2.
!!
!! NOTES
!!  global and local mean that the ...
!!
!!
!! SOURCE

type, public :: gqk_t

  integer :: cplex
  ! 1 if |g|^2,
  ! 2 if complex g (global)

  integer :: nb
  ! Number of bands included in the calculation for this spin
  ! Global as this dimension is not MPI-distributed due to (m,n) pairs.

  integer :: bstart = 1
  ! The fist band starts at bstart. (global index)
  ! Used to select bands inside energy window.

  integer :: my_npert

  integer :: glob_nk, glob_nq
  ! Total number of k/q points in global matrix.
  ! Use kzone_type and qzone_type to interpret these values (IBZ or BZ)
  ! Note that k-points can be filtered by erange.

  integer :: my_nk, my_nq
  ! Number of k/q points treated by this MPI proc.
  ! Used to loop and allocate local arrays.

  integer :: my_kstart, my_qstart
  ! Index of the first k/q point in the global matrix treated by this MPI proc

  integer :: nkibz, nqibz
  ! Number of k/q points in the IBZ.

  real(dp),allocatable :: kibz(:,:)
  ! k-points in the IBZ
  ! (3, nkibz)

  real(dp),allocatable :: qibz(:,:)
  ! q-points in the IBZ
  ! (3, nqibz)

  integer,allocatable :: my_k2ibz(:,:)
  ! (6, my_nk)
  ! Mapping my_kpoints --> kibz

  integer,allocatable :: my_q2ibz(:,:)
  ! (6, my_nq)
  ! Mapping my_qpoints --> qibz

  integer,allocatable :: my_qpk2ibz_k(:,:,:)
  ! (6, my_nq, my_nk)
  ! Mapping q + k --> kibz for each (q, k) treated by this MPI proc.

  logical :: need_velocities = .False.

  ! Compute group velocities if we are in transport mode or adaptive gaussian or
  ! tetrahedron with libtetrabz returning nesting condition.

  !complex(dp),allocatable :: my_vcar_k(:,:,:,:)
  ! group velocities v_{nm,k} for the k-points treated by this MPI proc.
  ! (3, nb, nb, my_nk)

  !complex(dp),allocatable :: my_vcar_qk(:,:,:,:)
  ! group velocities v_{nm,q+k} for the k/q-points treated by this MPI proc.
  ! (3, nb, nb, my_nq, my_nk)

  ! or alternatively:
  !complex(dp),allocatable :: my_vcar_kibz(:,:,:,:)

  !real(dp) :: erange(2) = zero

  character(len=10) :: kzone_type, qzone_type
    ! Either "ibz" or "bz"

  real(dp), allocatable :: vals(:,:,:,:,:,:)
    ! (cplex, nb, my_nq, nb, my_nk, my_npert)

  integer :: coords_qkp(3)
  ! coordinates of this processor in the Cartesian grid.

  type(xcomm_t) :: kpt_comm
    ! MPI communicator for parallelism over k-points

  type(xcomm_t) :: qpt_comm
   ! MPI communicator for parallelism over q-points

  type(xcomm_t) :: pert_comm
   ! MPI communicator for parallelism over atomic perturbations.

  !type(xcomm_t) :: gkq_comm

  !integer :: natom3
  ! 3 * natom

  !integer :: kptopt, qptopt
  ! Option for k-point generation.

  !integer :: timrev
  ! 1 if the use of time-reversal is allowed; 0 otherwise

  !real(dp),allocatable :: phfrq_ibz(:,:)
  ! (nibz, natom3)
  ! Phonon frequencies in the IBZ

  !real(dp),allocatable :: eigkbs_ibz(:, :, :)
  ! (nibz, nbcount, nsppol)
  ! Electron eigenvalues in the IBZ for nbcount states
  ! (not necessarly equal to global nband, see also bstart and bcount)

  !type(crystal_t), pointer :: cryst => null()
  ! Pointer to input structure (does not own memory)

  !type(htetra_t) :: tetra_k
  ! Used to evaluate delta(w - e_{k+q} +/- phw_q) with tetrahedron method.

 contains

  procedure :: get_mykpt => gqk_get_mykpt
  ! Return the k-point from my local index my_ik

  procedure :: get_myqpt => gqk_get_myqpt
  ! Return the q-point from my local index my_iq

  !procedure :: setup_kpoint => gstore_setup_kpoint
  ! Prepare tetrahedron method for given external k-point.

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
   ! Number of spins treated by this MPI rank

  integer :: comm
   ! Global communicator, inherited by caller so do not free in gstore_free.

  integer,allocatable :: my_spins(:)
   ! my_spins(my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2.

  !kmesh, qmesh
  integer :: kptrlatt(3, 3), qptrlatt(3, 3)
   ! k-mesh and q-mesh

  real(dp) :: kshift(3, 1), qshift(3, 1)
  ! k/q-mesh shift (well, q-mesh is usually gamma-centered)

  real(dp) :: klatt(3, 3), qlatt(3, 3)
  ! Reciprocal of lattice vectors for full kpoint grid.
  ! Used by tetra routines.

  type(xcomm_t) :: spin_comm
    ! MPI communicator over spin (nsppol = 2)

  type(gqk_t), allocatable :: gqk(:)
  ! (my_nspins)

contains

  procedure :: fill_bks_mask => gstore_fill_bks_mask
  ! Fill table used to read (b, k, s) wavefunctions from file.

  procedure :: get_mpw_gmax => gstore_get_mpw_gmax

  procedure :: spin2myis => gstore_spin2myis

  procedure :: free => gstore_free
  ! Free memory

  procedure :: ncwrite_path => gstore_ncwrite_path

  !procedure :: ncread => gstore_ncread_path

  procedure :: compute => gstore_compute

end type gstore_t

public :: gstore_build
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_build
!! NAME
!! gstore_build
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

function gstore_build(cplex, dtset, kzone_type, qzone_type, cryst, ebands, ifc, comm) result (gstore)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex, comm
 character(len=*),intent(in) :: kzone_type, qzone_type
 type(dataset_type),intent(in) :: dtset
 class(crystal_t),target,intent(in) :: cryst
 class(ebands_t),intent(in) :: ebands
 class(ifc_type),intent(in) :: ifc
 type(gstore_t), target :: gstore

!Local variables-------------------------------
!scalars
 integer :: all_nproc, my_rank, ik, ierr, out_nkibz
 integer :: my_is, my_ik, my_iq, spin, ik_ibz, iqk_ibz, natom3
 real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
 type(gqk_t),pointer :: gqk
!arrays
 integer :: comm_spin(ebands%nsppol)
 integer :: glob_nk_spin(ebands%nsppol), glob_nq_spin(ebands%nsppol)
 integer,allocatable :: myq2glob(:), myk2glob(:)
 !real(dp) :: rlatt(3,3)
 !integer :: out_kptrlatt(3,3)
 !real(dp),allocatable :: out_kibz(:,:), out_wtk(:)
#ifdef HAVE_MPI
 integer :: ndims, comm_cart, me_cart
 logical :: reorder
 integer,allocatable :: dims(:)
 logical,allocatable :: periods(:), keepdim(:)
#endif

!----------------------------------------------------------------------

 !print *, "in Gstore" !; return
 call cwtime(cpu, wall, gflops, "start")
 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 natom3 = 3 * cryst%natom

 ! Set basic parameters.
 gstore%nsppol = ebands%nsppol
 gstore%comm = comm

 !gstore%kptrlatt(3, 3)
 !gstore%qptrlatt(3, 3)
 !gstore%kshift(3, 1)
 !gstore%qshift(3, 1)
 !gstore%klatt(3, 3)
 !gstore%qlatt(3, 3)
 !gstore%spin_comm

 ! Define q-mesh for integration of the self-energy.
 ! Either q-mesh from DVDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation if q not in DDB)
 !new%ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 !if (all(dtset%eph_ngqpt_fine /= 0)) then
 !  new%ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 !end if

 ! TODO: Should fix bz2ibz to use the same convenetions as krank and listkk
 ! NB: only sigmaph seems to be using this optional argument

 !! Setup IBZ, weights and BZ. Always use q --> -q symmetry for phonons even in systems without inversion
 !qptrlatt = 0; qptrlatt(1, 1) = new%ngqpt(1); qptrlatt(2, 2) = new%ngqpt(2); qptrlatt(3, 3) = new%ngqpt(3)
 !!my_shiftq(:,1) = [0.1, 0, 0]
 !call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
 !                            new%nqibz, new%qibz, new%wtq, new%nqbz, new%qbz, bz2ibz=new%ind_qbz2ibz)

 ! Distribute spins and create mapping to spin index.
 !gstore%spin_comm%nproc = 1
 !if (gstore%nsppol == 2 .and. mod(all_nproc, 2) == 0) then
 !  gstore%spin_comm%nproc = 2
 !end if
 !if (any(dtset%eph_np_pqbks /= 0)) gstore%spin_comm%nproc = dtset%eph_np_pqbks(5)

 !if (gstore%nsppol == 2) then
 !  call xmpi_split_block(gstore%nsppol, gstore%spin_comm%value, gstore%my_nspins, gstore%my_spins)
 !  msg = sjoin("nsppol (", itoa(gstore%nsppol), ") < spin_comm_nproc (", itoa(gstore%spin_comm%nproc), ")")
 !  ABI_CHECK(gstore%my_nspins > 0, msg)
 !else
   ! No nsppol parallelism DOH!
   gstore%my_nspins = 1
   ABI_MALLOC(gstore%my_spins, (gstore%my_nspins))
   gstore%my_spins = 1
 !end if

 comm_spin(:) = comm

 ABI_MALLOC(gstore%gqk, (gstore%my_nspins))

 glob_nk_spin(:) = 20
 glob_nq_spin(:) = 20

 !new%natom3 = ifc%natom * 3
 !new%bstart = bstart
 !new%kptopt = kptopt
 !new%timrev = kpts_timrev_from_kptopt(new%kptopt)
 !new%nibz = nkibz
 !new%cryst => cryst
 !call alloc_copy(kibz, new%ibz)

 !! Get full BZ (new%nbz, new%bz) and new kptrlatt for tetra.
 !call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, out_nkibz, out_kibz, out_wtk, new%nbz, new%bz, &
 !                            new_kptrlatt=out_kptrlatt)

 !new%kptrlatt = out_kptrlatt
 !rlatt = out_kptrlatt; call matr3inv(rlatt, new%klatt)

 !ABI_CHECK(size(out_kibz, dim=2) == new%nibz, "mismatch in nkibz!")
 !ABI_FREE(out_kibz)
 !ABI_FREE(out_wtk)

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! Init for sequential execution.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   gqk%qpt_comm%nproc = 1
   gqk%kpt_comm%nproc = 1
   gqk%pert_comm%nproc = 1
   gqk%my_npert = natom3

   if (any(dtset%eph_np_pqbks /= 0)) then
     ! Use parameters from input file.
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
     ! Not the most efficient distribution if large number of MPI procs.
     ! TODO: Add parallelism over perturbations although it's less efficient than the parallelism over k
     ! It starts to be interesting if we implement symmetries in the k-integration though.

     gqk%kpt_comm%nproc = all_nproc
     gqk%qpt_comm%nproc = 1
     gqk%pert_comm%nproc = 1

     ! Handle parallelism over perturbations first.
     ! Use MPI communicator to distribute the 3 * natom perturbations to reduce memory requirements for DFPT potentials.
     ! Ideally, perturbations are equally distributed --> total number of CPUs should be divisible by 3 * natom.
     ! or at least, divisible by one integer i for i in [2, 3 * natom - 1].

     !if (nsppol == 2 .and. mod(nproc, 2) == 0) then
     !  spin_comm%nproc = 2
     !  kpt_comm%nproc = nproc / 2
     !else
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

   ! Consistency check.
   !if (pert_comm%nproc * qpt_comm%nproc * band_comm%nproc * kpt_comm%nproc * spin_comm%nproc /= nproc) then
   !  write(msg, "(a,i0,3a, 6(a,1x,i0))") &
   !    "Cannot create 5d Cartesian grid with total nproc: ", nproc, ch10, &
   !    "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
   !    "pert_nproc (", pert_comm%nproc, ") x qpt_nproc (", qpt_comm%nproc, ") x bsum_nproc (", band_comm%nproc, &
   !    ") x kcalc_nproc (", kpt_comm%nproc, ") x spin_nproc (", spin_comm%nproc, ") != ", &
   !    pert_comm%nproc * qpt_comm%nproc * band_comm%nproc * kpt_comm%nproc * spin_comm%nproc
   !  ABI_ERROR(msg)
   !end if

 end do

#ifdef HAVE_MPI
 ! For each spin, create 3d cartesian communicator: q-points, k-points and perturbations
 ndims = 3
 ABI_MALLOC(dims, (ndims))
 ABI_MALLOC(periods, (ndims))
 ABI_MALLOC(keepdim, (ndims))
 periods(:) = .False.; reorder = .False.

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)

   dims = [gqk%qpt_comm%nproc, gqk%kpt_comm%nproc, gqk%pert_comm%nproc]

   ! Note comm_spin
   call MPI_CART_CREATE(comm_spin(spin), ndims, dims, periods, reorder, comm_cart, ierr)

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

   if (my_rank == 0) then
     write(std_out, "(/,a)")" === MPI parallelism ==="
     write(std_out, "(a,i0)")"P Number of CPUs for parallelism over perturbations: ", gqk%pert_comm%nproc
     write(std_out, "(a,i0)")"P Number of perturbations treated by this CPU: ", gqk%my_npert
     write(std_out, "(a,i0)")"P Number of CPUs for parallelism over q-points: ", gqk%qpt_comm%nproc
     write(std_out, "(a,i0)")"P Number of CPUs for parallelism over k-points: ", gqk%kpt_comm%nproc
   end if

 end do ! my_is

 ABI_FREE(dims)
 ABI_FREE(periods)
 ABI_FREE(keepdim)
#endif

 ! At this point, we have the Cartesian grid (one per spin if any)
 ! and we can finally distribute dimensions.

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   gqk%kzone_type = trim(kzone_type)
   gqk%qzone_type = trim(qzone_type)

   gqk%glob_nk = glob_nk_spin(spin)
   gqk%glob_nq = glob_nq_spin(spin)
   gqk%cplex = 1
   gqk%nb = 1
   gqk%bstart = 1

   gqk%nkibz = 1
   gqk%nqibz = 1

   ! Split q-points and transfer symmetry tables.
   call xmpi_split_block(gqk%glob_nq, gqk%kpt_comm%value, gqk%my_nq, myq2glob)
   ABI_CHECK(gqk%my_nq > 0, "my_nq == 0")
   gqk%my_qstart = myq2glob(1)
   ABI_MALLOC(gqk%my_q2ibz, (6, gqk%my_nq))
   do my_iq=1,gqk%my_nq
     !my2glob(my_iq)
   end do
   ABI_FREE(myq2glob)

   ! Split k-points and transfer symmetry tables
   call xmpi_split_block(gqk%glob_nk, gqk%kpt_comm%value, gqk%my_nk, myk2glob)
   ABI_CHECK(gqk%my_nk > 0, "my_nk == 0")
   gqk%my_kstart = myk2glob(1)
   ABI_MALLOC(gqk%my_k2ibz, (6, gqk%my_nk))
   do my_ik=1,gqk%my_nk
     !my2glob(my_ik)
   end do
   ABI_FREE(myk2glob)

   ! Symmetry tables for q+k
   ABI_MALLOC(gqk%my_qpk2ibz_k, (6, gqk%my_nq, gqk%my_nk))

   ! Allocate storage for e-ph matrix elements.
   ABI_MALLOC(gqk%vals, (gqk%cplex, gqk%nb, gqk%my_nq, gqk%nb, gqk%my_nk, gqk%my_npert))
   gqk%vals = zero
 end do

 call cwtime_report(" gstore_build:", cpu, wall, gflops)

end function gstore_build
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_spin2myis
!! NAME
!! gstore_spin2myis
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

integer pure function gstore_spin2myis(gstore, spin) result(my_is)

!Arguments ------------------------------------
 class(gstore_t),intent(in) :: gstore
 integer,intent(in) :: spin

!----------------------------------------------------------------------

 do my_is=1,gstore%my_nspins
   if (gstore%my_spins(my_is) == spin) return
 end do
 my_is = 0

end function gstore_spin2myis
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_fill_bks_mask
!! NAME
!! gstore_fill_bks_mask
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

subroutine gstore_fill_bks_mask(gstore, mband, nkibz, nsppol, bks_mask)

!Arguments ------------------------------------
 class(gstore_t),target,intent(in) :: gstore
 integer,intent(in) :: mband, nkibz, nsppol
 logical,intent(out) :: bks_mask(mband, nkibz, nsppol)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, spin, ik_ibz, iqk_ibz
 type(gqk_t),pointer :: gqk

!----------------------------------------------------------------------

 bks_mask = .False.

 do my_is=1,gstore%my_nspins
   gqk => gstore%gqk(my_is)
   spin = gstore%my_spins(my_is)

   do my_ik=1,gqk%my_nk
     ik_ibz = gqk%my_k2ibz(1, my_ik)
     bks_mask(gqk%bstart:gqk%bstart + gqk%nb - 1, ik_ibz, spin) = .True.
     do my_iq=1,gqk%my_nq
       iqk_ibz = gqk%my_qpk2ibz_k(1, my_iq, my_ik)
       bks_mask(gqk%bstart:gqk%bstart + gqk%nb - 1, iqk_ibz, spin) = .True.
     end do
   end do
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

subroutine gstore_get_mpw_gmax(gstore, cryst, ecut, mpw, gmax)

!Arguments ------------------------------------
 class(gstore_t),target,intent(in) :: gstore
 class(crystal_t),intent(in) :: cryst
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
     kk = gqk%get_mykpt(my_ik, cryst)

     ! Compute G sphere, returning npw. Note istwfk == 1.
     call get_kg(kk, 1, ecut, cryst%gmet, onpw, gtmp)
     mpw = max(mpw, onpw)
     do ipw=1,onpw
       do ii=1,3
         gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
       end do
     end do
     ABI_FREE(gtmp)

     do my_iq=1,gqk%my_nq
       qpt = gqk%get_myqpt(my_iq, cryst)

       ! TODO: g0 umklapp here can enter into play!
       ! fstab should contains the max of the umlapp G-vectors.
       ! gmax could not be large enough!
       call get_kg(kk + qpt, 1, ecut, cryst%gmet, onpw, gtmp)
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

 ABI_SFREE(gstore%my_spins)
 call gstore%spin_comm%free()

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

pure function gqk_get_mykpt(gqk, my_ik, cryst) result (kpt)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 integer,intent(in) :: my_ik
 class(crystal_t),intent(in) :: cryst
 real(dp) :: kpt(3)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz

!----------------------------------------------------------------------

 ik_ibz = gqk%my_k2ibz(1, my_ik)
 !kpt = matmul(gqk%kibz(:, ik_ibz))

end function gqk_get_mykpt
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

pure function gqk_get_myqpt(gqk, my_iq, cryst) result (qpt)

!Arguments ------------------------------------
 class(gqk_t),intent(in) :: gqk
 integer,intent(in) :: my_iq
 class(crystal_t),intent(in) :: cryst
 real(dp) :: qpt(3)

!Local variables ------------------------------
!scalars
 integer :: iq_ibz

!----------------------------------------------------------------------

 iq_ibz = gqk%my_q2ibz(1, my_iq)
 !qpt = matmul(gqk%qibz(:, iq_ibz))

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

 ABI_SFREE(gqk%kibz)
 ABI_SFREE(gqk%qibz)
 ABI_SFREE(gqk%my_k2ibz)
 ABI_SFREE(gqk%my_q2ibz)
 ABI_SFREE(gqk%my_qpk2ibz_k)

 ABI_SFREE(gqk%vals)

 ! Communicators
 call gqk%kpt_comm%free()
 call gqk%qpt_comm%free()
 call gqk%pert_comm%free()

end subroutine gqk_free
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_compute
!! NAME
!!  gstore_compute
!!
!! FUNCTION
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

subroutine gstore_compute(gstore, wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, &
                          pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 class(gstore_t),target,intent(inout) :: gstore
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
 integer,parameter :: tim_getgh1c = 1, berryopt0 = 0, ider0 = 0, idir0 = 0
 integer,parameter :: useylmgr = 0, useylmgr1 = 0, master = 0, ndat1 = 1
 integer :: my_rank,nproc,mband,nsppol,nkibz,idir,ipert,iq_ibz, timrev
 integer :: cplex,natom,natom3,ipc,nspinor,onpw
 integer :: bstart_k,bstart_kq,nband_k,nband_kq,band_k, band_kq, ib_k, ib_kq !ib1,ib2,
 integer :: ik_ibz,ik_bz,ikq_bz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq
 integer :: ik_fs, my_ik, my_is !, myiq
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: ii,jj,ipw,mpw,my_mpw,mnb,ierr,cnt
 integer :: n1,n2,n3,n4,n5,n6,nspden
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,kq_count,nkpg,nkpg1
 integer :: comm_rpt, nesting, my_npert, my_ip, my_iq
 real(dp) :: cpu, wall, gflops, cpu_q, wall_q, gflops_q, cpu_k, wall_k, gflops_k, cpu_all, wall_all, gflops_all
 real(dp) :: sigma, ecut, eshift, eig0nk, dksqmax
 logical :: gen_eigenpb, isirr_k, isirr_kq
 type(wfd_t) :: wfd
 type(fstab_t),pointer :: fs
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 !type(phgamma_t) :: gams
 type(ddkop_t) :: ddkop
 type(xcomm_t) :: pert_comm, qpt_comm, kpt_comm, spin_comm
 type(krank_t) :: krank
 type(gqk_t),pointer :: gqk
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0bz_kq(3),g0_kq(3)
 integer :: work_ngfft(18),gmax(3),my_gmax(3),gamma_ngqpt(3) !g0ibz_kq(3),
 integer :: indkk_kq(6,1)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),gtmp(:,:),nband(:,:),wfd_istwfk(:)
 integer,allocatable :: my_pinfo(:,:), pert_table(:,:) !, qibz_done(:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3), lf(2),rg(2),res(2), vk(3), vkq(3)
 real(dp) :: phfrq(3*cryst%natom)
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: displ_cart(:,:,:,:), displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:), kinpw1(:), kpg1_k(:,:), kpg_k(:,:), dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:), ffnl1(:,:,:,:), ph3d(:,:,:), ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:), gkk_atm(:,:,:,:),gkq_nu(:,:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:), kets_k(:,:,:), h1kets_kq(:,:,:), cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:), vlocal(:,:,:,:), vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:), ylm_k(:,:), ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:), gvnlx1(:,:), work(:,:,:,:)
 real(dp),allocatable :: gs1c(:,:), v1_work(:,:,:,:), vcar_ibz(:,:,:,:)
 real(dp),allocatable :: wt_ek(:,:), wt_ekq(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(fstab_t),target,allocatable :: fstab(:)
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
 timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ! FFT meshes
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Find Fermi surface k-points
 ! TODO: support kptopt, change setup of k-points if tetra: fist tetra weights then k-points on the Fermi surface!
 ABI_MALLOC(fstab, (nsppol))
 call fstab_init(fstab, ebands, cryst, dtset, comm)

 ! Define q-mesh. eph_ngqpt_fine activates the Fourier interpolation of the DFPT potentials.
 gamma_ngqpt = ifc%ngqpt; if (all(dtset%eph_ngqpt_fine /= 0)) gamma_ngqpt = dtset%eph_ngqpt_fine

 !call wrtout(std_out, sjoin("q-mesh for the phonon linewidths:", ltoa(gamma_ngqpt)))
 !call wrtout(std_out, sjoin("Will compute", itoa(gams%nqibz), "q-points in the IBZ"))

 !call xmpi_barrier(comm)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 ! FIXME
 gqk => gstore%gqk(1)
 my_npert = gqk%my_npert

 if (pert_comm%nproc > 1) then
   ! Activate parallelism over perturbations
   ! Build table with list of perturbations treated by this CPU inside pert_comm
   call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, pert_comm%value)
   call dvdb%set_pert_distrib(my_npert, natom3, my_pinfo, pert_table, pert_comm%value)
   ABI_FREE(my_pinfo)
   ABI_FREE(pert_table)
 end if

 call wrtout([std_out, ab_out], " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
 ! Prepare Fourier interpolation of DFPT potentials.
 comm_rpt = xmpi_comm_self
 !comm_rpt = bqs_comm%value
 call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k/k+q treated by this proc are store.

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

 call wfd%print(header="Wavefunctions for gstore calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

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

 call gstore%get_mpw_gmax(cryst, ecut, mpw, gmax)

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

 ! Create ddkop object to compute group velocities if needed.
 !
 !   1) precompute group velocities in the IBZ and the ihave_ikibz_spin file (common to all procs)
 !   2) Use symmetries to reconstruct v_kq from vcar_ibz
 !
 ! NB: All procs store in memory the same set of Bloch states inside the energy window.

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 call cwtime(cpu, wall, gflops, "start", msg=" Computing v_nk matrix elements for all states on the FS...")
 ii = huge(1); jj = -1
 do spin=1,nsppol
   ii = min(ii, fstab(spin)%bmin)
   jj = max(jj, fstab(spin)%bmax)
 end do
 ABI_CALLOC(vcar_ibz, (3, ii:jj, nkibz, nsppol))
 ABI_MALLOC(cgwork, (2, mpw * wfd%nspinor))

 cnt = 0
 do spin=1,nsppol
   fs => fstab(spin)
   do ik_ibz=1,ebands%nkpt
     kk = ebands%kptns(:, ik_ibz)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     ! NB: The two checks below are global --> all procs will cycle.
     if (all(bks_mask(:, ik_ibz, spin) .eqv. .False.)) cycle
     if (npw_k == 1) cycle
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.

     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

     do band_k=fs%bmin,fs%bmax
       if (.not. bks_mask(band_k, ik_ibz, spin)) cycle
       call wfd%copy_cg(band_k, ik_ibz, spin, cgwork)
       eig0nk = ebands%eig(band_k, ik_ibz, spin)
       vk = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cgwork, cwaveprj0)
       vcar_ibz(:, band_k, ik_ibz, spin) = vk

     end do
   end do
 end do ! spin

 call xmpi_sum(vcar_ibz, comm, ierr)
 ABI_FREE(cgwork)
 call cwtime_report(" Velocities", cpu, wall, gflops)

 ABI_FREE(bks_mask)

 ! Build krank object to find k-points
 krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

 ! Loop over my spins.
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   !spin = gams%my_spins(my_is)
   !fs => fstab(spin)

   gqk => gstore%gqk(my_is)

   ! Loop over my set of q-points
   !do my_iq=1,gams%my_nqibz
   do my_iq=1,gqk%my_nq
     call cwtime(cpu_q, wall_q, gflops_q, "start")

     !iq_ibz = gqk%my_q2ibz(1, my_iq)
     qpt = gqk%get_myqpt(my_iq, cryst)

     !iq_ibz = gams%my_iqibz(my_iq)
     !qpt = gams%qibz(:, iq_ibz)

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

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ! maxnb is the maximum number of bands crossing the FS, used to dimension arrays.
     !mnb = fs%maxnb
     mnb = gqk%nb
     ABI_MALLOC(bras_kq, (2, mpw*nspinor, mnb))
     ABI_MALLOC(kets_k, (2, mpw*nspinor, mnb))
     ABI_MALLOC(h1kets_kq, (2, mpw*nspinor, mnb))
     ABI_MALLOC(gkk_atm, (2, mnb, mnb, natom3))
     ABI_MALLOC(gkq_nu, (2, mnb, mnb, natom3))

     !do my_ik=1,gams%my_nfsk_q
     do my_ik=1,gqk%my_nk
       call cwtime(cpu_k, wall_k, gflops_k, "start")

       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       !ik_fs = gams%my_ifsk_q(my_ik)
       kk = fs%kpts(:, ik_fs)
       ik_ibz = fs%indkk_fs(1, ik_fs); isym_k = fs%indkk_fs(2, ik_fs)
       trev_k = fs%indkk_fs(6, ik_fs); g0_k = fs%indkk_fs(3:5,ik_fs)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       kk_ibz = ebands%kptns(:,ik_ibz)

       !kk = gqk%get_mykpt(my_ik, cryst)
       !gqk%my_k2ibz(:, my_ik)
       !gqk%my_qpk2ibz_k(:, my_iq, my_ik)

       ! Number of bands crossing the Fermi level at k
       bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)

       ! Find k+q in the extended zone and extract symmetry info. cycle if k+q not in FS.
       ! Be careful here because there are two umklapp vectors to be considered:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qpt; ikq_bz = fs%findkg0(kq, g0bz_kq)

       ! Skip this point if kq does not belong to the FS window.
       if (ikq_bz == -1) cycle

       call krank%get_mapping(1, kq, dksqmax, cryst%gmet, indkk_kq, cryst%nsym, cryst%symafm, cryst%symrel, timrev, &
                              use_symrec=.False.)

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
       bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz)
       ABI_CHECK(nband_k <= mnb .and. nband_kq <= mnb, "wrong nband")

       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kk, kk_ibz, bstart_k, nband_k, spin, mpw, fs%indkk_fs(:,ik_fs), cryst, &
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
       call mkffnl(psps%dimekb, 1, psps%ekb, ffnlk, psps%ffspl,&
                   cryst%gmet, cryst%gprimd, ider0, idir0, psps%indlmn, kg_k, kpg_k, kk, psps%lmnmax, &
                   psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, psps%ntypat, &
                   psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum, &
                   comm=pert_comm%value)

       ! Compute k+q+G vectors
       nkpg1 = 3 * dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G)
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       call mkffnl(psps%dimekb, 1, psps%ekb, ffnl1, psps%ffspl, cryst%gmet, cryst%gprimd, ider0, idir0, &
                   psps%indlmn, kg_kq, kpg1_k, kq, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg1, &
                   npw_kq, psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_kq, ylmgr_kq, &
                   comm=pert_comm%value)

       ! Loop over all my atomic perturbations and compute gkk_atm.
       gkk_atm = zero
       do my_ip=1,my_npert
         idir = dvdb%my_pinfo(1, my_ip); ipert = dvdb%my_pinfo(2, my_ip); ipc = dvdb%my_pinfo(3, my_ip)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)

         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,my_ip), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk, kq, idir, ipert, &                    ! In
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

       end do ! my_ip (loop over my_npert atomic perturbations)

       ABI_FREE(gs1c)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)

       ! Collect gkk_atm inside pert_comm so that all procs can operate on the data.
       if (pert_comm%nproc > 1) call xmpi_sum(gkk_atm, pert_comm%value, ierr)

       ! Get g in the phonon representation.
       call ephtk_gkknu_from_atm(mnb, mnb, 1, natom, gkk_atm, phfrq, displ_red, gkq_nu)

       ! Save results.
       !  (2, mnb, mnb, natom3))
       ! (cplex, nb, my_nq, nb, my_nk, my_npert)
       !if (gqk%cplex == 1) then
       !  gqk%vals(1, :, my_iq, :, my_ik, my_ip) =
       !else if (gqk%cplex == 2) then
       !  gqk%vals(1:2, :, my_iq, :, my_ik, my_ip) =
       !end if

       ! Compute group velocities if we are in transport mode or adaptive gaussian or
       ! tetrahedron with libtetrabz returning nesting condition.
       !need_velocities = .False.

       if (gqk%need_velocities) then
         ! Compute diagonal matrix elements of velocity operator with DFPT routines
         ! Velocities are in Cartesian coordinates.
         !
         ! If k+q is not in the IBZ, we need to recostruct the value by symmetry using v(Sq) = S v(q).
         ! Use transpose(R) because we are using the tables for the wavefunctions
         ! In this case listkk has been called with symrec and use_symrec=False
         ! so q_bz = S^T q_ibz where S is the isym_kq symmetry

         !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, kg_k)
         do ib_k=1,nband_k
           band_k = ib_k + bstart_k - 1
           vk = vcar_ibz(:, band_k, ik_ibz, spin)
           if (.not. isirr_k) then
             vk = matmul(transpose(cryst%symrel_cart(:,:,isym_k)), vk)
             if (trev_k /= 0) vk = -vk
           end if
           !vk = ddkop%get_vdiag(ebands%eig(band_k, ik_ibz, spin), &
           !                     istwf_k, npw_k, wfd%nspinor, kets_k(:,:,ib_k), cwaveprj0)
           fs%vk(:,ib_k) = vk
         end do

         !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kq, istwf_kq, npw_kq, kg_kq)
         do ib_kq=1,nband_kq
           band_kq = ib_kq + bstart_kq - 1
           vkq = vcar_ibz(:, band_kq, ikq_ibz, spin)
           if (.not. isirr_kq) then
             vkq = matmul(transpose(cryst%symrel_cart(:,:,isym_kq)), vkq)
             if (trev_kq /= 0) vkq = -vkq
           end if
           !vkq = ddkop%get_vdiag(ebands%eig(band_kq, ikq_ibz, spin), &
           !                                 istwf_kq, npw_kq, wfd%nspinor, bras_kq(:,:,ib_kq), cwaveprj0)
           fs%vkq(:,ib_kq) = vkq
         end do
       end if

       !if (my_ik < 20 .or. (fs%nkfs > 100 .and. mod(my_ik, 200) == 0)) then
       !  write(msg,'(4(a,i0),a,f8.2)')" q-point [", iq_ibz, "/", gams%nqibz, "] k-point [", my_ik, "/", gams%my_nfsk_q, "]"
       !  call cwtime_report(msg, cpu_k, wall_k, gflops_k)
       !end if
     end do ! my_ik

     ABI_FREE(bras_kq)
     ABI_FREE(kets_k)
     ABI_FREE(h1kets_kq)
     ABI_FREE(gkk_atm)
     ABI_FREE(gkq_nu)
     ABI_FREE(v1scf)
     ABI_FREE(vlocal1)

     write(msg,'(2(a,i0),a)')" Computation of my q-point [", my_iq, "/", gqk%my_nq, "]"
     call cwtime_report(msg, cpu_q, wall_q, gflops_q, end_str=ch10)
   end do ! my_iq
 end do ! my_is

 call cwtime_report(" phonon linewidths k-loop", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10)

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
 call krank%free()

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
 integer :: my_rank, ncid, spin, my_is

! *************************************************************************

 my_rank = xmpi_comm_rank(gstore%comm)

 ! Master node writes basic objects and dimensions
 ! then we are gonna write the big arrays with MPI-IO and hdf5 groups
 ! to account for spin and spin-dependent gqk dimensions.
 if (my_rank == 0) then
   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))
   NCF_CHECK(nf90_close(ncid))
 end if
 call xmpi_barrier(gstore%comm)

 do spin=1,gstore%nsppol
   my_is = gstore%spin2myis(spin)
   if (my_is /= 0) then
     ! this is the tricky business
     !call gstore%gqk(my_is)%ncwrite_path(path)
   end if
   call xmpi_barrier(gstore%comm)
 end do

 ! =====================
 ! === Write dimensions
 ! =====================
 !ncerr = nctk_def_dims(ncid, [ &
 !  nctkdim_t("max_number_of_states", ebands%mband), &
 !  nctkdim_t("number_of_spinor_components", ebands%nspinor), &
 !  nctkdim_t("nshiftk", ebands%nshiftk)], &
 !  defmode=.True.)
 !NCF_CHECK(ncerr)

 !! Define k-points
 !ncerr = nctk_def_arrays(ncid, [&
 !  nctkarr_t("reduced_coordinates_of_kpoints", "dp", "number_of_reduced_dimensions, number_of_kpoints"), &
 !  nctkarr_t("monkhorst_pack_folding", "int", "number_of_vectors") &
 !])
 !NCF_CHECK(ncerr)

 !! Define states section.
 !ncerr = nctk_def_arrays(ncid, [&
 !  nctkarr_t("number_of_states", "int", "number_of_kpoints, number_of_spins"), &
 !  nctkarr_t("smearing_scheme", "char", "character_string_length")  &
 !])
 !NCF_CHECK(ncerr)

 !ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "number_of_electrons"])
 !NCF_CHECK(ncerr)
 !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
 !NCF_CHECK(ncerr)

 ! Some variables require the specifications of units.
 !NCF_CHECK(nctk_set_atomic_units(ncid, "eigenvalues"))
 !NCF_CHECK(nctk_set_atomic_units(ncid, "fermi_energy"))

 ! Write data.
 !NCF_CHECK(nctk_set_datamode(ncid))
 !NCF_CHECK(nf90_put_var(ncid, vid("fermi_energy"), ebands%fermie))

 ! Write Abinit variables
 !NCF_CHECK(nctk_set_datamode(ncid))
 !NCF_CHECK(nf90_put_var(ncid, vid("tphysel"), ebands%tphysel))

!contains
! integer function vid(vname)
!   character(len=*),intent(in) :: vname
!   vid = nctk_idname(ncid, vname)
! end function vid

end subroutine gstore_ncwrite_path
!!***

!----------------------------------------------------------------------

end module m_gstore
!!***

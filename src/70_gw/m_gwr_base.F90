!!****m* ABINIT/m_gwr_base
!! NAME
!!  m_gwr_base
!!
!! FUNCTION
!!  Objects and procedures to implement GW in real-space and imaginary time.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2021 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

module m_gwr_base

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_slk
 use m_wfd
 use m_hdr
 use m_ebands

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use defs_abitypes,   only : mpi_type
 use m_time,          only : cwtime, cwtime_report, sec2str
 use m_io_tools,      only : iomode_from_fname !, file_exists, is_open, open_file
 use m_numeric_tools, only : get_diag, isdiagmat
 use m_fstrings,      only : sjoin, itoa
 use m_krank,         only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_crystal,       only : crystal_t
 use m_dtset,         only : dataset_type
 use m_fftcore,       only : get_kg, sphereboundary, ngfft_seq, getng, print_ngfft !, kgindex
 use m_fft,           only : fft_ug, fftbox_plan3, fftbox_plan3_init
 use m_fft_mesh,      only : times_eikr !, times_eigr, ig2gfft, get_gftt, calc_ceikr, calc_eigr rotate_fft_mesh
 use m_kpts,          only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_sort, kpts_pack_in_stars
 use m_wfk,           only : wfk_read_ebands
 use m_pawtab,        only : pawtab_type

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_gwr_base/desc_t
!! NAME
!! desc_t
!!
!! FUNCTION
!!  Parameters related to a two-point function such as
!!  gvectors, k/q point and tables used for zero padded FFTs.
!!
!! SOURCE

 type,public :: desc_t

   integer :: istwfk = 1
   ! Storage mode for this k point.

   integer :: npw = -1
   ! Number of plane-waves for this k-point.

   integer :: nfft = -1
   ! Number of FFT points

   integer :: mgfft = -1
   ! maximum size of 1D FFTs

   integer :: ngfft(18) = -1
   ! FFT mesh

   real(dp) :: point(3)
   ! k/q point in reduced coordinates.

   integer,allocatable :: gvec(:,:)
   ! gvec(3, npw)
   ! G vectors in reduced coordinates.
   ! For each k, the g-vectors are ordered according to |k+g|
   ! so that we can use the same array for Green, chi and Sigma_x
   ! Note that this array is global i.e. it is not MPI-distributed in G-space, only in k-space.

   integer,allocatable :: gbound(:,:)
   ! gbound(2 * mgfft + 8, 2)
   ! sphere boundary info for zero-padded FFT

 contains

   !procedure :: set_ngfft => desc_set_ngfft
   !procedure :: rotate => rotate

   procedure :: free => desc_free
   ! Free memory.

 end type desc_t
!!***

!type, extends(matrix_scalapack) :: twopoint_func_t
!   type(desc_t), pointer :: desc => null()
!   character(len=100) :: name
!
!   type(matrix_scalapack) :: g_g
!   type(matrix_scalapack) :: r_r
!   type(matrix_scalapack) :: g_r
!
!end type twopoint_func_t

!----------------------------------------------------------------------

!!****t* m_gwr_base/gwr_t
!! NAME
!! gwr_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: gwr_t

   integer :: nsppol = 1, nspinor = -1
   ! Number of independent spin polarizations and number of spinor components.

   integer :: my_nspins = -1
   ! Number of independent spin polarizations treated by this MPI proc

   integer :: nkbz = -1, nkibz = -1
   ! Number of k-points in the BZ/IBZ

   !integer :: nk = -1
   ! nkbz if GWr, nkibz if GWrk

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

   integer,allocatable :: my_qibz_inds(:), my_qbz_inds(:)
   ! FIXME: Do we really need my_qibz_inds and my_qibz_inds(:)

   integer :: ntau = -1
   ! Total number of imaginary time points.

   integer :: my_ntau = -1
   ! Number of imaginary time/frequency points treated by this MPI rank.

   !integer :: nomega = -1, my_nomega = -1

   logical :: use_supercell = .True.
   ! True if we are using the supercell formalism.
   ! False if we are using the mixed-space approach with convolutions in k-space.

   integer :: ngkpt(3) = -1
   ! Number of divisions in k-point mesh.

   integer,allocatable :: my_spins(:)
   ! (my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2.

   integer,allocatable :: my_itaus(:)
   ! (my_ntau)
   ! Indirect table giving the tau indices treated by this MPI rank.

   real(dp),allocatable :: tau_mesh(:)
   ! tau_mesh(ntau)
   ! Imaginary tau mesh
   ! PS: The mesh deepends on the quantity we want to compute e.g. RPA or Sigma.

   real(dp),allocatable :: iw_mesh(:)
   ! tau_mesh(ntau)
   ! Imaginary frequency mesh.

   real(dp),allocatable :: t2w_cos_wgs(:,:)
   ! (ntau, ntau)
   ! weights for cosine transform.

   real(dp),allocatable :: t2w_sin_wgs(:,:)
   ! (ntau, ntau)
   ! weights for sine transform.

   integer :: green_mpw = -1
   ! Max number of g-vectors for Green's function over k-points.

   integer :: chi_mpw = -1
   ! Max number of g-vectors for Chi over q-points.

   !integer :: sigma_mpw = -1
   ! Max number of g-vectors for Sigma over q-points.

   integer :: ngfft(18) = -1
   integer :: mgfft = -1, nfft = -1

   !integer :: chi_ngfft(18) = -1
   !integer :: sigma_ngfft(18) = -1

   type(desc_t),allocatable :: green_desc_mykf(:)
   ! (my_nkbz)
   ! Descriptor for Green's functions

   !type(desc_t),allocatable :: sigma_desc_k(:)
   ! (nk)
   ! Descriptor for self-energy

   type(desc_t),allocatable :: chi_desc_myqf(:)
   ! (my_nqbz)
   ! Descriptor for chi

   !integer,allocatable :: gvec(:,:)
   ! gvec(3, gvec_size, nk)
   ! For each k, the gvectors are ordered according to |k+g|
   ! so that we can use the same array for G, chi and Sigma_x
   ! Note that this array is global i.e. is not MPI-distributed in G-space, only in k-space.

   !integer, allocatable :: npw_green(:)
   !integer, allocatable :: npw_chi(:)
   !integer, allocatable :: npw_sigma(:)
   ! Arrays of shape (nk) giving the number of G-vectors in gvec.
   ! These quantities are computed from the values of ecut, ecuteps and ecutsigx
   ! specified in input.

   integer :: coords_gtks(4) = 0
   ! Cartesian coordinates of this processor in the Cartesian grid.

   integer :: comm
   ! Initial MPI communicator with all procs involved in the calculation

   type(xcomm_t) :: spin_comm
   ! MPI communicator for parallelism over spins.

   type(xcomm_t) :: kpt_comm
   ! MPI communicator for k/q-point distribution.

   type(xcomm_t) :: g_comm
   ! MPI communicator for g/r distribution

   type(xcomm_t) :: tau_comm
   ! MPI communicator for imag time distribution

   type(xcomm_t) :: gtau_comm

   !type(xcomm_t) :: omega_comm
   ! MPI communicator for imag frequency distribution

   type(dataset_type), pointer :: dtset => null()

   type(crystal_t), pointer  :: cryst => null()

   type(ebands_t), pointer  :: ebands => null()

   type(pseudopotential_type), pointer :: psps => null()

   type(pawtab_type), pointer :: pawtab(:) => null()

   type(mpi_type),pointer :: mpi_enreg => null()

   type(matrix_scalapack),allocatable :: gf_g_gp(:,:,:,:)
   ! Occupied/Empty Green's function.
   ! (2, my_nkbz, my_ntau, my_nspins)

   type(matrix_scalapack),allocatable :: chi_g_gp(:,:,:)
   ! (my_nqbz, my_ntau, my_nspins)

   type(matrix_scalapack),allocatable :: wc_g_gp(:,:)
   ! Correlated screened Coulomb interaction

   !type(matrix_scalapack),allocatable :: sigo_g_gp(:,:,:)

   !type(matrix_scalapack),allocatable :: sige_g_gp(:,:,:)

   !character(len=fnlen) :: wfk0_path = ABI_NOFILE
   ! Path to the WFK file with the KS wavefunctions.

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

   !integer,allocatable :: my_ibz_inds(:,:)

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

 contains

   procedure :: ggp_to_gpr => gwr_ggp_to_gpr

   procedure :: cos_transform  => gwr_cos_transform

   procedure :: free => gwr_free
   ! Free memory.

   !procedure :: print => gwr_print

   procedure :: build_gtau_from_wfk => gwr_build_gtau_from_wfk

   procedure :: build_chi => gwr_build_chi

   !procedure :: rpa_energy => gwr_rpa_energy
   !procedure :: g0w0 => gwr_g0w0

 end type gwr_t

 public :: gwr_new
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_new
!! NAME
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

function gwr_new(dtset, cryst, psps, pawtab, ebands, mpi_enreg, comm) result (gwr)

!Arguments ------------------------------------
!scalars
 type(gwr_t),target :: gwr
 type(dataset_type),target,intent(in) :: dtset
 type(crystal_t),target,intent(in) :: cryst
 type(pseudopotential_type),target,intent(in) :: psps
 type(pawtab_type),target,intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(ebands_t),target,intent(in) :: ebands
 type(mpi_type),target,intent(in) :: mpi_enreg
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: me_fft0 = 0, paral_fft0 = 0, nproc_fft1 = 1
 integer,parameter :: qptopt1 = 1, timrev1 = 1, master = 0
 integer :: my_is, my_it, my_ikf, my_iqf, ii, my_ir, my_nr, npwsp, col_bsize, ebands_timrev
 integer :: ig1, ig2, sc_nfft, sc_augsize, my_nshiftq, spin, ik_ibz, ik_bz, iq_bz, iq_ibz
 integer :: all_nproc, my_rank !, mgfft, nfft
 real(dp) :: ecut_eff
 real(dp) :: cpu, wall, gflops
 character(len=5000) :: msg
 !type(desc_t) :: desc_k, desc_kpq
 type(desc_t),pointer :: desc_k, desc_q
 type(processor_scalapack) :: slk_processor
 type(krank_t) :: qrank, krank_ibz
!arrays
 integer :: sc_ngfft(18), ngqpt(3), qptrlatt(3,3)
 integer,allocatable :: kbz2ibz(:,:), qbz2ibz(:,:)
 !integer,allocatable :: kibz2bz(:), qibz2bz(:), qglob2bz(:,:), kglob2bz(:,:)
 real(dp) :: my_shiftq(3,1), kk_ibz(3), kk_bz(3), qq_bz(3), kpq_bz(3), qq_ibz(3)
 real(dp),allocatable :: wtk(:), kibz(:,:)
 complex(dp),allocatable :: chi_scbox(:), go_scbox(:), ge_scbox(:)
 integer,allocatable :: count_ibz(:)
 integer,parameter :: ndims = 4
 integer :: comm_cart, me_cart, ierr
 logical :: reorder
 integer :: dims(ndims)
 logical :: periods(ndims), keepdim(ndims)

! *************************************************************************

 all_nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 call cwtime(cpu, wall, gflops, "start")

 gwr%dtset => dtset
 gwr%cryst => cryst
 gwr%psps => psps
 gwr%pawtab => pawtab
 gwr%ebands => ebands
 ! TODO Make sure fermie is inside the gap if semiconductor. perhaps this shoud be delegated to gwr_driver
 gwr%kibz => ebands%kptns
 gwr%wtk => ebands%wtk
 gwr%mpi_enreg => mpi_enreg

 gwr%comm = comm
 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%use_supercell = .True.

 ! =======================
 ! Setup k-mesh and q-mesh
 ! =======================

 ! Get full kBZ associated to ebands
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
                             gwr%nkibz, kibz, wtk, gwr%nkbz, gwr%kbz) !, bz2ibz=bz2ibz)
                             !new_kptrlatt=gwr%kptrlatt, new_shiftk=gwr%kshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ABI_FREE(wtk)

 ! In principle kibz should be equal to ebands%kptns
 ABI_CHECK(gwr%nkibz == ebands%nkpt, "nkibz != ebands%nkpt")
 ABI_CHECK(all(abs(ebands%kptns - kibz) < tol12), "ebands%kibz != kibz")

 if (.not. (isdiagmat(ebands%kptrlatt) .and. ebands%nshiftk == 1)) then
   ABI_ERROR("GWR requires ngkpt with one shift!")
 end if
 gwr%ngkpt = get_diag(ebands%kptrlatt)

 ! Note symrel and use_symrec=.False. in get_mapping.
 ! This means that this table can be used to symmetrize wavefunctions in cgtk_rotate.
 ! TODO This ambiguity should be removed. Change cgtk_rotate so that we can use the symrec convention.

 ABI_MALLOC(kbz2ibz, (6, gwr%nkbz))
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 krank_ibz = krank_from_kptrlatt(gwr%nkibz, kibz, ebands%kptrlatt, compute_invrank=.False.)

 if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, gwr%nkbz, gwr%kbz, kbz2ibz) /= 0) then
   ABI_ERROR("Cannot map kBZ to IBZ!")
 end if
 call krank_ibz%free()

 !call get_ibz2bz(gwr%nkibz, gwr%nkbz, kbz2ibz, kibz2bz, ierr)
 !ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for k-points")

 ABI_FREE(kibz)

 ! Setup qIBZ, weights and BZ.
 ! Always use q --> -q symmetry even in systems without inversion
 ! TODO: Might add input variable to rescale the q-mesh.
 my_nshiftq = 1; my_shiftq = zero
 qptrlatt = ebands%kptrlatt
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt1, my_nshiftq, my_shiftq, &
                             gwr%nqibz, gwr%qibz, gwr%wtq, gwr%nqbz, gwr%qbz)
                             !new_kptrlatt=gwr%qptrlatt, new_shiftk=gwr%qshift,
                             !bz2ibz=new%ind_qbz2ibz)  # FIXME

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(qbz2ibz, (6, gwr%nqbz))

 qrank = krank_from_kptrlatt(gwr%nqibz, gwr%qibz, qptrlatt, compute_invrank=.False.)

 if (kpts_map("symrec", timrev1, cryst, qrank, gwr%nqbz, gwr%qbz, qbz2ibz) /= 0) then
   ABI_ERROR("Cannot map qBZ to IBZ!")
 end if
 call qrank%free()

 ! Order qbz by stars and rearrange entries in qbz2ibz table.
 call kpts_pack_in_stars(gwr%nqbz, gwr%qbz, qbz2ibz)

 ! ====================
 ! Setup tau/omega mesh
 ! ====================
 gwr%ntau = 1
 ABI_MALLOC(gwr%tau_mesh, (gwr%ntau))
 ABI_MALLOC(gwr%iw_mesh, (gwr%ntau))
 ABI_MALLOC(gwr%t2w_cos_wgs, (gwr%ntau, gwr%ntau))
 ABI_MALLOC(gwr%t2w_sin_wgs, (gwr%ntau, gwr%ntau))

 ! ====================
 ! Planewave basis set
 ! ====================

 ! ============================
 ! Find good nproc distribution
 ! ============================
 ! TODO: for the time being all procs for tau.

 gwr%tau_comm%nproc = all_nproc

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

 call xmpi_comm_free(comm_cart)
#endif

 ! Block-distribute dimensions and allocate redirection table: local index --> global index.
 call xmpi_split_block(gwr%ntau, gwr%tau_comm%value, gwr%my_ntau, gwr%my_itaus)
 call xmpi_split_block(gwr%nsppol, gwr%spin_comm%value, gwr%my_nspins, gwr%my_spins)

 ! Distribute k-points in the full BZ, transfer symmetry tables.
 ! Finally find the number of IBZ points that should be stored in memory.

 call xmpi_split_block(gwr%nkbz, gwr%kpt_comm%value, gwr%my_nkbz, gwr%my_kbz_inds)
 ABI_MALLOC(gwr%my_kbz2ibz, (6, gwr%my_nkbz))
 gwr%my_kbz2ibz = kbz2ibz(:, gwr%my_kbz_inds(:))

 ABI_ICALLOC(count_ibz, (gwr%nkibz))
 do my_ikf=1,gwr%my_nkbz
   ik_ibz = gwr%my_kbz2ibz(1, my_ikf)
   count_ibz(ik_ibz) = count_ibz(ik_ibz) + 1
 end do

 gwr%my_nkibz = count(count_ibz > 0)
 ABI_MALLOC(gwr%my_kibz_inds, (gwr%my_nkibz))
 ii = 0
 do ik_ibz=1,gwr%nkibz
   if (count_ibz(ik_ibz) > 0) then
     ii = ii + 1; gwr%my_kibz_inds(ii) = ik_ibz
   end if
 end do

 ABI_FREE(kbz2ibz)
 ABI_FREE(count_ibz)

 ! Distribute q-points in the full BZ, transfer symmetry tables.
 ! Finally find the number of IBZ points that should be stored in memory.
 ABI_REMALLOC(count_ibz, (gwr%nqibz))
 count_ibz = 0

 call xmpi_split_block(gwr%nqbz, gwr%kpt_comm%value, gwr%my_nqbz, gwr%my_qbz_inds)
 ABI_MALLOC(gwr%my_qbz2ibz, (6, gwr%my_nqbz))
 gwr%my_qbz2ibz = qbz2ibz(:, gwr%my_qbz_inds(:))

 do my_iqf=1,gwr%my_nqbz
   iq_ibz = gwr%my_qbz2ibz(1, my_iqf)
   count_ibz(iq_ibz) = count_ibz(iq_ibz) + 1
 end do

 gwr%my_nqibz = count(count_ibz > 0)
 ABI_MALLOC(gwr%my_qibz_inds, (gwr%my_nqibz))
 ii = 0
 do iq_ibz=1,gwr%nqibz
   if (count_ibz(iq_ibz) > 0) then
     ii = ii + 1; gwr%my_qibz_inds(ii) = iq_ibz
   end if
 end do

 ABI_FREE(qbz2ibz)
 ABI_FREE(count_ibz)

 ! Build FFT descriptors for Green's functions.
 ABI_MALLOC(gwr%green_desc_mykf, (gwr%my_nkbz))

 gwr%ngfft = 0
 do my_ikf=1,gwr%my_nkbz
   ik_bz = gwr%my_kbz_inds(my_ikf)
   kk_bz = gwr%kbz(:, ik_bz)
   desc_k => gwr%green_desc_mykf(my_ikf)
   desc_k%point = kk_bz; desc_k%istwfk = 1
   !gwr%green_desc_mykf(my_ikf) = desc_new(ecut_eff, kpoint, gwr%cryst%gmet, ngfft=ngfft)
   ! Calculate G-sphere from input ecut.
   ecut_eff = dtset%ecut * dtset%dilatmx ** 2
   call get_kg(kk_bz, desc_k%istwfk, ecut_eff, gwr%cryst%gmet, desc_k%npw, desc_k%gvec)
   call getng(dtset%boxcutmin, dtset%chksymtnons, ecut_eff, cryst%gmet, &
              kk_bz, me_fft0, gwr%mgfft, gwr%nfft, gwr%ngfft, nproc_fft1, cryst%nsym, paral_fft0, &
              cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%green_mpw = max(gwr%green_mpw, desc_k%npw)
 end do

 ! Build FFT descriptors for chi.
 ABI_MALLOC(gwr%chi_desc_myqf, (gwr%my_nqbz))

 do my_iqf=1,gwr%my_nqbz
   iq_bz= gwr%my_qbz_inds(my_iqf)
   !iq_ibz= gwr%my_qibz_inds(my_iqd)
   qq_bz = gwr%qbz(:, my_iqf)
   desc_q => gwr%chi_desc_myqf(my_iqf)
   desc_k%point = qq_bz; desc_q%istwfk = 1
   !gwr%chi_desc_myqf(my_iqf) = desc_new(ecut_eff, qpoint, gwr%cryst%gmet, ngfft=chi_ngfft)
   call get_kg(qq_bz, desc_q%istwfk, dtset%ecuteps, gwr%cryst%gmet, desc_q%npw, desc_q%gvec)
   call getng(dtset%boxcutmin, dtset%chksymtnons, dtset%ecuteps, cryst%gmet, &
              qq_bz, me_fft0, gwr%mgfft, gwr%nfft, gwr%ngfft, nproc_fft1, cryst%nsym, &
              paral_fft0, cryst%symrel, cryst%tnons, unit=dev_null)

   gwr%chi_mpw = max(gwr%chi_mpw, gwr%chi_desc_myqf(my_iqf)%npw)
 end do
 call print_ngfft(gwr%ngfft, header="FFT mesh for Chi", unit=std_out)

 ! Now we know the value of ngfft. Setup tables for zero-padded FFTs.
 !nfft = product(ngfft(1:3)); mgfft = maxval(ngfft(1:3))
 do my_ikf=1,gwr%my_nkbz
   desc_k => gwr%green_desc_mykf(my_ikf)
   desc_k%nfft = gwr%nfft; desc_k%mgfft = gwr%mgfft; desc_k%ngfft = gwr%ngfft
   ABI_MALLOC(desc_k%gbound, (2 * gwr%mgfft + 8, 2))
   call sphereboundary(desc_k%gbound, desc_k%istwfk, desc_k%gvec, gwr%mgfft, desc_k%npw)
 end do

 do my_iqf=1,gwr%my_nqbz
   desc_q => gwr%chi_desc_myqf(my_iqf)
   desc_q%nfft = gwr%nfft; desc_q%mgfft = gwr%mgfft; desc_q%ngfft = gwr%ngfft
   ABI_MALLOC(desc_q%gbound, (2 * gwr%mgfft + 8, 2))
   call sphereboundary(desc_q%gbound, desc_q%istwfk, desc_q%gvec, gwr%mgfft, desc_q%npw)
 end do

 ! ======================================
 ! Allocate arrays storing Go_k and Chi_q
 ! ======================================
 ABI_MALLOC(gwr%gf_g_gp, (2, gwr%my_nkbz, gwr%my_ntau, gwr%my_nspins))
 ABI_MALLOC(gwr%chi_g_gp, (gwr%my_nqbz, gwr%my_ntau, gwr%my_nspins))

 call init_scalapack(slk_processor, gwr%g_comm%value) !# grid_shape=[1, gwr%g_comm%value]

 do my_is=1,gwr%my_nspins
   do my_it=1,gwr%my_ntau

     ! Allocate G_k(g, g')
     do my_ikf=1,gwr%my_nkbz
       npwsp = gwr%green_desc_mykf(my_ikf)%npw * gwr%nspinor
       col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
       call init_matrix_scalapack(gwr%gf_g_gp(1, my_ikf, my_it, my_is), &
          npwsp, npwsp, slk_processor, 1) ! , tbloc=[npwsp, col_bsize])
       call init_matrix_scalapack(gwr%gf_g_gp(2, my_ikf, my_it, my_is), &
          npwsp, npwsp, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
     end do

    ! Allocate chi_k(g, g')
    do my_iqf=1,gwr%my_nqbz
      npwsp = gwr%chi_desc_myqf(my_iqf)%npw * gwr%nspinor
      col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
      call init_matrix_scalapack(gwr%chi_g_gp(my_iqf, my_it, my_is), &
         npwsp, npwsp, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
    end do

   end do
 end do

 call cwtime_report(" gwr_new:", cpu, wall, gflops)

end function gwr_new
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_free
!! NAME
!! gwr_free
!!
!! FUNCTION
!!  Free memory
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

subroutine gwr_free(gwr)

!Arguments ------------------------------------
 class(gwr_t), intent(inout) :: gwr

!Local variables-------------------------------
 integer :: ii, jj, my_ikf, my_iqf

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
 ABI_SFREE(gwr%tau_mesh)
 ABI_SFREE(gwr%iw_mesh)
 ABI_SFREE(gwr%t2w_cos_wgs)
 ABI_SFREE(gwr%t2w_sin_wgs)

 ! Free descriptors
 if (allocated(gwr%green_desc_mykf)) then
   do my_ikf=1,gwr%my_nkbz
     call gwr%green_desc_mykf(my_ikf)%free()
   end do
   ABI_FREE(gwr%green_desc_mykf)
 end if
 if (allocated(gwr%chi_desc_myqf)) then
   do my_iqf=1,gwr%my_nqbz
     call gwr%chi_desc_myqf(my_iqf)%free()
   end do
   ABI_FREE(gwr%chi_desc_myqf)
 end if
 ! Free scalapack matrices
 if (allocated(gwr%gf_g_gp)) then
   call slk_free_array(gwr%gf_g_gp)
   ABI_FREE(gwr%gf_g_gp)
 end if
 if (allocated(gwr%chi_g_gp)) then
   call slk_free_array(gwr%chi_g_gp)
   ABI_FREE(gwr%chi_g_gp)
 end if
 if (allocated(gwr%wc_g_gp)) then
   call slk_free_array(gwr%wc_g_gp)
   ABI_FREE(gwr%wc_g_gp)
 end if

 ! Free MPI communicators
 call gwr%spin_comm%free()
 call gwr%g_comm%free()
 call gwr%tau_comm%free()
 call gwr%kpt_comm%free()
 call gwr%gtau_comm%free()

end subroutine gwr_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_build_gtau_from_wfk
!! NAME
!!  gwr_build_gtau_from_wfk
!!
!! FUNCTION
!!  Build the Green's function in imaginary time for k-points in the IBZ
!!  from the WFK file
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

subroutine gwr_build_gtau_from_wfk(gwr, wfk0_path)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 character(len=fnlen),intent(in) :: wfk0_path

!Local variables-------------------------------
!scalars
 integer :: mband, nkibz, nsppol, my_is, my_it, my_iki, spin, ik_ibz, my_nb, band, ib, npw_k, istwf_k
 real(dp) :: f_nk, eig_nk, tau, ef
 real(dp) :: cpu, wall, gflops
 type(wfd_t) :: wfd
 type(ebands_t) :: ebands
 type(hdr_type) :: wfk0_hdr
 type(dataset_type), pointer :: dtset
!arrays
 integer,allocatable :: nband(:,:), wfd_istwfk(:), my_bands(:)
 real(dp),allocatable :: cgwork(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 dtset => gwr%dtset

 ebands = wfk_read_ebands(wfk0_path, gwr%comm, out_hdr=wfk0_hdr)
 call wfk0_hdr%vs_dtset(gwr%dtset)

 !cryst = wfk0_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k wavevectors
 ! treated by this MPI rank are stored.

 nkibz = ebands%nkpt; nsppol = ebands%nsppol; mband = ebands%mband

 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 ! Init bks_mask: each MPI proc reads only the (k-points, spin).
 ! Distribute bands inside g_comm
 ! FIXME: Smart algo to make memory scale.
 call xmpi_split_block(mband, gwr%g_comm%value, my_nb, my_bands)

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     bks_mask(my_bands(:), ik_ibz, spin) = .True.
   end do
 end do

 ! Impose istwfk = 1 for all k-points.
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 call wfd_init(wfd, gwr%cryst, gwr%pawtab, gwr%psps, keep_ur, mband, nband, nkibz, dtset%nsppol, bks_mask, &
   dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, gwr%ngfft, &
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, gwr%comm)

 call wfd%print(header="Wavefunctions for GWR calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)
 call ebands_free(ebands)
 call wfk0_hdr%free()

 ! Read KS wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ! ==================================
 ! Build Green's functions in g-space
 ! ==================================
 ! TODO: Make sure that gvec in gwr and wfd agree with each other.
 ! NB: G_k is constructed for k in the IBZ, then we rotate the k-point to obtain G_k in the BZ.
 ef = ebands%fermie

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     !ik_bz = gwr%kibz2bz(ik_ibz)

     ABI_MALLOC(cgwork, (2, npw_k * wfd%nspinor))
     !wfd%kdata(ik_ibz)%kg_k

     do ib=1,my_nb
       band = my_bands(ib)
       f_nk = gwr%ebands%occ(band, ik_ibz, spin)
       eig_nk = gwr%ebands%eig(band, ik_ibz, spin) - ef
       call wfd%copy_cg(band, ik_ibz, spin, cgwork)

       !ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)

       !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk_ibz, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)
       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       !call wfd%sym_ug_kg(ecut, kk_bz, kk_ibz, bstart_k, nband_k, spin, mpw, gqk%my_k2ibz(:, my_ikf), cryst, &
       !                   work_ngfft, work, istwf_k, npw_k, kg_k, kets_k)

       do my_it=1,gwr%my_ntau
         tau = gwr%tau_mesh(gwr%my_itaus(my_it))
         !cgwork(ig1) x cgwork(ig2) * exp(-tau * eig_nk)
         !gwr%gf_g_gp(1, my_ikf, my_it, my_is)%buffer_cplx =
         !gwr%gf_g_gp(2, my_ikf, my_it, my_is)%buffer_cplx =
       end do ! my_it

     end do ! ib
     ABI_FREE(cgwork)

   end do ! my_ikf
 end do ! my_is

 call wfd%free()
 ABI_FREE(my_bands)

 call cwtime_report(" gwr_build_gtau_from_wfk:", cpu, wall, gflops)

end subroutine gwr_build_gtau_from_wfk

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_ggp_to_gpr
!! NAME
!!  gwr_ggp_to_gpr
!!
!! FUNCTION
!!
!! G_k(g, g') --> G_k(g', r) for each k in the BZ treated by me for given spin and tau.
!!
!!   1) FFT Transform the first index (local): G(g,g',it) --> G(r,g',it)
!!   2) MPI transposition: G(r,g',it) --> G^*(g',r,it)
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

subroutine gwr_ggp_to_gpr(gwr, my_it, my_is, gp_r)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 integer,intent(in) :: my_it, my_is
 type(matrix_scalapack),intent(inout) :: gp_r(2, gwr%my_nkbz)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1
 integer :: my_ikf, ig2, ioe
 integer :: ik_ibz, isym_k, trev_k, g0_k(3)
 real(dp) :: cpu, wall, gflops
 logical :: isirr_k
 type(matrix_scalapack), pointer :: g_gp
 type(matrix_scalapack) :: r_gp, g_r
 type(desc_t),pointer :: desc

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 do my_ikf=1,gwr%my_nkbz
   desc => gwr%green_desc_mykf(my_ikf)

   ik_ibz = gwr%my_kbz2ibz(my_ikf, 1); isym_k = gwr%my_kbz2ibz(my_ikf, 2)
   trev_k = gwr%my_kbz2ibz(my_ikf, 6); g0_k = gwr%my_kbz2ibz(my_ikf, 3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))

   ! Allocate r_gp scalapack matrix to store X(r, g')
   !call init_matrix_scalapack(r_gp, desc%nfft * gwr%nspinor, g_gp%sizeb_global(2), &
   !                           g_gp%processor, desc%istwfk) !, tbloc=tbloc)

   do ioe=1,2
     g_gp => gwr%gf_g_gp(ioe, my_ikf, my_it, my_is)
     !g_gp => gwr%gf_g1g2(my_ikf, ioe, my_it, my_is)

     ! Rotate k to get G_k in the BZ.
     if (.not. isirr_k) then
       !g_gp => gwr%gf_g_gp(ioe, my_ikf, my_it, my_is)
     end if

     ! FFT. results stored in r_gp
     do ig2=1,g_gp%sizeb_local(2)
       call fft_ug(desc%npw, desc%nfft, gwr%nspinor, ndat1, &
                   desc%mgfft, desc%ngfft, desc%istwfk, desc%gvec, desc%gbound, &
                   g_gp%buffer_cplx(:, ig2), r_gp%buffer_cplx(:, ig2))

       ! Multiply by e^{ik.r}
       !call times_eikr(desc%point, desc%ngfft, desc%nfft, ndat1, r_gp%buffer_cplx(:, ig2)
     end do ! ig2

     ! Transpose: X(r, g') -> Y(g', r) and take complex conjugate
     !call gp_r%free()
     !call init_matrix_scalapack(gp_r(ioe, my_ikf), g_gp%sizeb_global(2), &
     !                           desc%nfft * gwr%nspinor, g_gp%processor, desc%istwfk) !, tbloc=tbloc)

     !call r_gp%ptranspose("C", gp_r)
   end do ! ioe
 end do ! my_ikf

 call r_gp%free()

 call cwtime_report(" gwr_ggp_to_gpr:", cpu, wall, gflops)

end subroutine gwr_ggp_to_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_cos_transform
!! NAME
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

subroutine gwr_cos_transform(gwr, what, mode, sum_spins)

!Arguments ------------------------------------
 class(gwr_t),target,intent(in) :: gwr
 character(len=*),intent(in) :: what, mode
 logical,optional,intent(in) :: sum_spins

!Local variables-------------------------------
!scalars
 integer :: my_iqf, my_is, ig1, ig2, my_it, iw, ierr
 real(dp) :: cpu, wall, gflops
 logical :: sum_spins_
 type(desc_t),pointer :: desc
 type(matrix_scalapack), pointer :: chit(:)
!arrays
 real(dp),allocatable :: my_weights(:,:)
 complex(dp),allocatable :: cwork_t(:), cwork_wglb(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 sum_spins_ = .False.; if (present(sum_spins)) sum_spins_ = sum_spins

 ABI_MALLOC(my_weights, (gwr%ntau, gwr%my_ntau))
 ABI_MALLOC(cwork_t, (gwr%my_ntau))
 ABI_MALLOC(cwork_wglb, (gwr%my_ntau))

 ! Extract my weights from the global array according to mode.
 do my_it=1,gwr%my_ntau
   !itau = gwr%my_itaus(my_it)
   !my_weights(:, my_it) = cos(omega * tau) * global_weighs(:, itau)
   select case (mode)
   case ("w2t")
   case ("t2w")
   case default
     ABI_ERROR(sjoin("Wrong mode", mode))
   end select
 end do

 !my_nspins = gwr%my_nspins
 !if (what == "wc") my_nspins = 1

 do my_is=1,gwr%my_nspins
   do my_iqf=1,gwr%my_nqbz

     select case (what)
     case ("chi")
       desc => gwr%chi_desc_myqf(my_iqf)
       chit => gwr%chi_g_gp(my_iqf, 1:gwr%my_ntau, my_is)
     !case ("wc")
     !  desc => gwr%chi_desc_myqf(my_iqf)
     !  chit => gwr%wc_g_gp(my_iqf, 1:gwr%my_ntau)
     case default
       ABI_ERROR(sjoin("Invalid value for what:", what))
     end select

     do ig2=1,chit(1)%sizeb_local(2)
       do ig1=1,chit(1)%sizeb_local(1)
         ! Extract matrix elements as a function of tau
         ! Here we can block over ig1 and call zgemm
         ! in order to reduced the number of MPI communications.
         do my_it=1,gwr%my_ntau
           cwork_t(my_it) = chit(my_it)%buffer_cplx(ig1, ig2)
         end do
         do iw=1,gwr%ntau
           cwork_wglb(iw) = dot_product(my_weights(iw, :), cwork_t)
         end do

         call xmpi_sum(cwork_wglb, gwr%tau_comm%value, ierr)

         ! update values for this (g1, g2), now we have a slice of gg' in frequency space.
         do my_it=1,gwr%my_ntau
           !itau = gwr%my_itaus(my_it)
           !chit(my_it)%buffer_cplx(ig1, ig2) = cwork_wglb(itau)
         end do

       end do ! ig1
     end do ! ig2

   end do ! my_iqf
 end do ! spin

 if (gwr%nsppol == 2 .and. sum_spins_) then
    !call xmpi_sum(cwork_wglb, gwr%spin_comm%value, ierr)
 endif

 ABI_FREE(cwork_t)
 ABI_FREE(cwork_wglb)
 ABI_FREE(my_weights)

 call cwtime_report(" gwr_cos_transform:", cpu, wall, gflops)

end subroutine gwr_cos_transform
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_wt_from_chiw
!! NAME
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

subroutine gwr_wt_from_chiw(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer :: my_iqf, my_iw, ig1, ig2, ig1_glob, ig2_glob, ierr, my_is
 real(dp) :: cpu, wall, gflops
 logical :: is_gamma
 type(desc_t), pointer :: desc_q
 type(matrix_scalapack), pointer :: chi
 type(matrix_scalapack) :: eps

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 !ABI_MALLOC(gwr%wc_g_gp, (gwr%my_nqbz, gwr%my_ntau))

 do my_iqf=1,gwr%my_nqbz
   desc_q => gwr%chi_desc_myqf(my_iqf)
   is_gamma = all(abs(desc_q%point) < tol12)  ! TODO: Handle q --> 0 limit

   !npwsp = gwr%chi_desc_myqf(my_iqf)%npw * gwr%nspinor
   !col_bsize = npwsp / gwr%g_comm%nproc
   !if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
   !!col_bsize = 50
   !call init_matrix_scalapack(gwr%chi_g_gp(my_iqf, my_iw, my_is), &
   !   npwsp, npwsp, slk_processor, 1) ! , tbloc=[npwsp, col_bsize])

   do my_is=1,gwr%my_nspins  ! FIXME
     do my_iw=1,gwr%my_ntau
       chi => gwr%chi_g_gp(my_iqf, my_iw, my_is)

       ! Build epsilon(q, iw) =  delta_{g,g'} − v_q(g,g') chi_q(g,g, iw).
       do ig2=1,chi%sizeb_local(2)
         do ig1=1,chi%sizeb_local(1)
           call chi%loc2glob(ig1, ig2, ig1_glob, ig2_glob)
           !vc_qg1 = ??
           !vc_qg2 = ??
           !ceps = - chi%buffer_cplx(ig1, ig2) ! * vc_qg1 * vc_qg2
           !if (ig1_glob == ig2_glob) eps = one - ceps
         end do ! ig1
       end do ! ig2

       ! IMPORTANT NOTE: PZGETRF requires square block decomposition i.e., MB_A = NB_A.
       ! This means here I need to redistribute the data before calling zinvert.
       !call slk_zinvert(eps)
       !call chi%free()

       ! Build W(q, iw) = e^{-1}_q(g,g',iw) v_q(g,g')
       ! Remove bare v
       do ig2=1,eps%sizeb_local(2)
         do ig1=1,eps%sizeb_local(1)
           call eps%loc2glob(ig1, ig2, ig1_glob, ig2_glob)
           !vc_qg1 = ??
           !vc_qg2 = ??
           !eps%buffer_cplx(ig1, ig2) = eps%buffer_cplx(ig1, ig2) ! * vc_qg1 * vc_qg2
         end do ! ig1
       end do ! ig2

       !gwr%wc_g_gp, (my_iqf, my_iw, my_is) = ??
       call eps%free()
     end do  ! my_it
   end do ! my_is
 end do ! my_iqf

 ! Perform cosine transform from iw to tau to get W(tau)
 call gwr%cos_transform("wc", "w2t")

 call cwtime_report(" gwr_wt_from_chiw:", cpu, wall, gflops)

end subroutine gwr_wt_from_chiw
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/slk_alloc_transpose
!! NAME
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

type(matrix_scalapack) function slk_init_transpose(in_mat) result (out_mat)

!Arguments ------------------------------------
 class(matrix_scalapack), intent(in) :: in_mat
 !class(matrix_scalapack), intent(out) :: out_mat

!Local variables-------------------------------
 !type(processor_scalapack) :: slk_processor

! *************************************************************************

 !call init_matrix_scalapack(out_mat, in_mat%sizeb_global(2), in_mat%sizeb_global(1), slk_processor, 1, tbloc=tbloc)

end function slk_init_transpose
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/slk_ptranspose
!! NAME
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

subroutine slk_ptranspose(in_mat, trans, out_mat)

!Arguments ------------------------------------
 class(matrix_scalapack), intent(in) :: in_mat
 character(len=1),intent(in) :: trans
 class(matrix_scalapack), intent(inout) :: out_mat

! *************************************************************************

 !call out_mat%free()
 !call init_matrix_scalapack(out_mat, in_mat%sizeb_global(2), in_mat%sizeb_global(1), in_mat%slk_processor, 1, tbloc=tbloc)

 ! prototype:
 !call pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)

 if (allocated(in_mat%buffer_cplx)) then
   ! Transposes a complex distributed matrix, conjugated.
   ! sub(C):=beta*sub(C) + alpha*conjg(sub(A)'),
   select case (trans)
   case ("C")
#ifdef HAVE_LINALG_SCALAPACK
     call pztranc(in_mat%sizeb_global(2), in_mat%sizeb_global(1), cone, in_mat%buffer_cplx, 1, 1, &
                  in_mat%descript%tab, czero, out_mat%buffer_cplx, 1, 1, out_mat%descript%tab)
#else
     out_mat%buffer_cplx = transpose(conjg(in_mat%buffer_cplx))
#endif

   case ("N")
     ! sub(C):=beta*sub(C) + alpha*sub(A)',
#ifdef HAVE_LINALG_SCALAPACK
     call pztranu(in_mat%sizeb_global(2), in_mat%sizeb_global(1), cone, in_mat%buffer_cplx, 1, 1, &
                  in_mat%descript%tab, czero, out_mat%buffer_cplx, 1, 1, out_mat%descript%tab)
#else
     out_mat%buffer_cplx = transpose(in_mat%buffer_cplx)
#endif

   case default
     ABI_ERROR(sjoin("Invalid value for trans:", trans))
   end select

 else if (allocated(in_mat%buffer_real)) then
#ifdef HAVE_LINALG_SCALAPACK
   call pdtran(in_mat%sizeb_global(2), in_mat%sizeb_global(1), one, in_mat%buffer_real, 1, 1, &
               in_mat%descript%tab, zero, out_mat%buffer_real, 1, 1, out_mat%descript%tab)
#else
     out_mat%buffer_real = transpose(in_mat%buffer_real)
#endif

 else
   ABI_BUG("Neither buffer_cplx nor buffer_real are allocated!")
 end if

end subroutine slk_ptranspose
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/desc_new
!! NAME
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

type(desc_t) function desc_new(ecut, kpoint, gmet) result (new)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: ecut
!arrays
 real(dp),intent(in) :: kpoint(3), gmet(3,3)
 !integer,intent(out) :: npw_k
 !integer,allocatable,intent(out) :: kg_k(:,:)

!Local variables-------------------------------
!scalars

! *************************************************************************

 ! Calculate G-sphere from input ecut.
 !call get_kg(kpoint, istwf_1, ecut, gmet, new%npw_k, new%kg_k)

 !if (present(ngfft))
 !  !mgfft = maxval(ngfft(1:3))
 !  !ABI_MALLOC(gbound, (2*mgfft+8, 2))
 !  !call sphereboundary(gbound, istwf_1, kg_k, mgfft, npw_k)
 !end if

end function desc_new
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/desc_free
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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine desc_free(desc)

!Arguments ------------------------------------
 class(desc_t) :: desc

! *************************************************************************

 ABI_SFREE(desc%gvec)
 ABI_SFREE(desc%gbound)

end subroutine desc_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr_base/gwr_build_chi
!! NAME
!!  gwr_build_chi
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

subroutine gwr_build_chi(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer :: my_is, my_it, my_ikf, my_iqf, ii, my_ir, my_nr, npwsp, col_bsize
 integer :: ig1, ig2, sc_nfft, sc_augsize, spin, ik_ibz, ik_bz, iq_bz, iq_ibz, ierr
 !integer :: all_nproc, my_rank !, mgfft, nfft
 real(dp) :: cpu, wall, gflops
 character(len=5000) :: msg
 type(desc_t),pointer :: desc_k, desc_q
 type(processor_scalapack) :: slk_processor
 type(matrix_scalapack) :: go_k_gpr, go_kpq_gpr, ge_k_gpr, ge_kpq_gpr, chi_r_gp
!arrays
 integer :: sc_ngfft(18)
 integer,allocatable :: green_scg(:,:), chi_scg(:,:)
 real(dp) :: kk_ibz(3), kk_bz(3), qq_bz(3), kpq_bz(3), qq_ibz(3)
 complex(dp),allocatable :: chi_scbox(:), go_scbox(:), ge_scbox(:)
 type(matrix_scalapack),allocatable :: gf_gp_r(:,:), chiq_gp_r(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ! The g-vectors in the supercell for G and chi.
 ABI_MALLOC(green_scg, (3, gwr%green_mpw))
 ABI_MALLOC(chi_scg, (3, gwr%chi_mpw))

 ! TODO: Very likely this is needed only for GWR
 ! Can be allocated inside loops to render code a bit clearer.
 !if (gwr%use_supercell) then
 ABI_MALLOC(chiq_gp_r, (gwr%my_nqbz))
 do my_iqf=1,gwr%my_nqbz
   npwsp = gwr%chi_desc_myqf(my_iqf)%npw * gwr%nspinor
   col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
   call init_matrix_scalapack(chiq_gp_r(my_iqf), &
        npwsp, gwr%nfft * gwr%nspinor, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
 end do

 ! Prepare FFT in the supercell from R' -> G'
 ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
 ! Perhaps the safest approach would be to generate the plan on the fly.
 sc_ngfft = gwr%ngfft
 sc_ngfft(1:3) = gwr%ngkpt * gwr%ngfft(1:3)
 sc_ngfft(4:6) = sc_ngfft(1:3)
 sc_nfft = product(sc_ngfft(1:3))
 sc_augsize = product(sc_ngfft(4:6))

 !cal fftbox_plan3_init(sc_plan_gp2rp, ndat, dims, embed, fftalg, fftcache, isign)
 !cal fftbox_plan3_init(sc_plan_rp2gp, ndat, dims, embed, fftalg, fftcache, isign)

 ABI_MALLOC(go_scbox, (sc_nfft))
 ABI_MALLOC(ge_scbox, (sc_nfft))
 ABI_MALLOC(chi_scbox, (sc_nfft))
 !end if ! gwr%use_supercell

 ABI_MALLOC(gf_gp_r, (2, gwr%my_nkbz))

 ! Loop over my spins and my taus.
 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_it=1,gwr%my_ntau
     !tau = gwr%tau_mesh(gwr%my_itaus(my_it))

     ! G_k(g, g') --> G_k(g', r) for each k in the BZ treated by me.
     call gwr%ggp_to_gpr(my_it, my_is, gf_gp_r)

     ! ============================
     ! GWr algorithm with supercell
     ! ============================
     !if (gwr%use_supercell) then

     ! Loop over the second index r (MPI-distributed in g_comm)
     my_nr = gf_gp_r(1, 1)%sizeb_local(2)
     do my_ir=1,my_nr

       ! Insert G_k(g',r) in G'-space in the supercell FFT box.
       ! Note that we need to take the union of (k, g') for k in the BZ.
       ! TODO: Here I need a specialized version of sphere that does not initialize out to zero as I need to accumulate
       do my_ikf=1,gwr%my_nkbz
         desc_k => gwr%green_desc_mykf(my_ikf)
         do ii=1,desc_k%npw
           green_scg(:,ii) = nint(gwr%kbz(:,my_ikf) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ii)
         end do

         gf_gp_r(1, my_ikf)%buffer_cplx(:, my_ir) = zero
         !call sphere(cg,ndat,npw,go_scbox,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)

         gf_gp_r(2, my_ikf)%buffer_cplx(:, my_ir) = zero
         !call sphere(cg,ndat,npw,ge_scbox,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
       end do

       ! Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
       call xmpi_sum(go_scbox, gwr%kpt_comm%value, ierr)
       call xmpi_sum(ge_scbox, gwr%kpt_comm%value, ierr)

       ! Compute chi(R',r) for this r and go back to G'=k+g' space immediately.
       !call fftbox_execute_ip(sc_plan_gp2rp, go_scbox)
       !call fftbox_execute_ip(sc_plan_gp2rp, ge_scbox)

       chi_scbox = go_scbox * ge_scbox
       !call fftbox_execute_ip(sc_plan_rp2gp, chi_scbox)

       ! Extract chi_q(g, r) on the g-sphere from box.
       ! not necessarly equal to the one used for the Green's function.
       ! NB: For the time being I'm assuming k-mesh == q-mesh.
       do my_iqf=1,gwr%my_nqbz
         desc_q => gwr%chi_desc_myqf(my_iqf)
         do ii=1,desc_q%npw
           chi_scg(:,ii) = nint(gwr%qbz(:, my_iqf) * gwr%ngkpt) + gwr%ngkpt * desc_q%gvec(:,ii)
         end do

         !call sphere(cg,ndat,npw,cfft_o,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
         chiq_gp_r(my_iqf)%buffer_cplx(:, my_ir) = zero
       end do

     end do ! my_ir

     do my_iqf=1,gwr%my_nqbz
       desc_q => gwr%chi_desc_myqf(my_iqf)
       !chi_r_gp = chiq_gp_r(my_iq)%alloc_transpose()
       !call chiq_gp_r(my_iqf)%ptranspose("C", chi_r_gp)
       call chiq_gp_r(my_iqf)%free()
       ! FFT along the first dimension: chi_q(r, g') --> chi_q(g, g')
       do ig2=1,chi_r_gp%sizeb_local(2)
         chi_r_gp%buffer_cplx(:, ig2) = zero
         !call fft_ur(npw_k, nfft, nspinor, ndat, mgfft, ngfft, istwf_k, kg_k, gbound_k, ur, ug)
         gwr%chi_g_gp(my_iqf, my_it, my_is)%buffer_cplx(:, ig2) = zero ! ug
       end do
       call chi_r_gp%free()
     end do

     !else
     ! ====================================================
     ! Mixed-space algorithm in unit cell with convolutions
     ! ====================================================
     ! Each MPI proc has all the G_k in the IBZ thus we only need to rotate to get the BZ on the fly.

     !do iq_ibz=1,gwr%nqibz
     do my_iqf=1,gwr%my_nqibz
       ! FIXME
       qq_ibz = gwr%qibz(:, iq_ibz)
       !qq_bz = gwr%qbz(:, iq_bz)
       desc_q => gwr%chi_desc_myqf(iq_ibz)
       do my_ikf=1,gwr%my_nkbz
         ! Use symmetries to get Go_kq and Ge_k from the corresponding images in the IBZ.
         kk_bz = gwr%kbz(:, my_ikf)
         kpq_bz = kk_bz + qq_ibz
         !ik_ibz = ??
         !ikpq_ibz = ??
         !desc_k = gwr%green_desc_mykf(ik_ibz)%rotate(symtab_k)
         !desc_kpq = gwr%green_desc_mykf(ikpq_ibz)%rotate(symtab_kpq)

         !go_k_gpr = gwr%gf_g_gp(1, ik_ibz, my_it, my_is)%rotate_fftg1_and_ptran(desc_k, "C", symtab_k)
         !ge_k_gpr = gwr%ge_g_gp(2, ik_ibz, my_it, my_is)%rotate_fftg1_and_ptran(desc_k, "C", symtab_k)
         !go_kpq_gpr = gwr%gf_g_gp(1, ikpq_ibz, my_it, my_is)%rotate_fftg1_and_ptran(desc_kpq, "C", symtab_kpq)
         !ge_kpq_gpr = gwr%ge_g_gp(2, ikpq_ibz, my_it, my_is)%rotate_fftg1_and_ptran(desc_kpq, "C", symtab_kpq)

         my_nr = go_k_gpr%sizeb_local(2)
         do my_ir=1,my_nr
            go_k_gpr%buffer_cplx(:, my_ir) = zero
            ge_k_gpr%buffer_cplx(:, my_ir) = zero
            go_kpq_gpr%buffer_cplx(:, my_ir) = zero
            ge_kpq_gpr%buffer_cplx(:, my_ir) = zero

            !chi_fft = go_rpr * ge_rpr
            ! from r' to g'
            !chi_gp_r(:, my_ir) = FFT_rg(chi_fft)
         end do

         !call desc_k%free()
         !call desc_kpq%free()

         call go_kpq_gpr%free()
         call go_k_gpr%free()
         call ge_kpq_gpr%free()
         call ge_k_gpr%free()

         !call chi_gp_r%ptranspose("C", chi_r_gp)
         !call chi_gp_r%free()
         !do ig2=1, chi_rp_g%sizeb_local(2)
           !chi_r_gp%buffer_cplx(:, ig2)
           !!gwr%chi_g_gp(iq_ibz, my_it, my_is)%buffer_cplx(:, ig2) = FFT_rg
         !end do
         !call chi_r_gp%free()
       end do ! my_ikf
     end do ! iq_ibz

   end do ! my_it
 end do ! spin

 do my_ikf=1,gwr%my_nkbz
   call gf_gp_r(1, my_ikf)%free()
   call gf_gp_r(2, my_ikf)%free()
 end do
 ABI_FREE(gf_gp_r)

 if (allocated(chiq_gp_r)) then
   do my_iqf=1,gwr%my_nqbz
     call chiq_gp_r(my_iqf)%free()
   end do
   ABI_FREE(chiq_gp_r)
 end if

 ABI_FREE(go_scbox)
 ABI_FREE(ge_scbox)
 ABI_FREE(chi_scbox)
 ABI_FREE(green_scg)
 ABI_FREE(chi_scg)

 call cwtime_report(" gwr_buid_chi:", cpu, wall, gflops)

 ! Transform irreducible chi from imaginary tau to imaginary omega
 ! and sum spins in order to get total chi if collinear spins.
 call gwr%cos_transform("chi", "t2w", sum_spins=.True.)

end subroutine gwr_build_chi
!!***

end module m_gwr_base
!!***

!!****m* ABINIT/m_gwr
!! NAME
!!  m_gwr
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

module m_gwr

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_slk
 use m_hdr
 use m_ebands
 use netcdf
 use m_nctk
 use m_dtfil
 use m_yaml

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use defs_abitypes,   only : mpi_type
 use m_time,          only : cwtime, cwtime_report, sec2str
 use m_io_tools,      only : iomode_from_fname !, file_exists, is_open, open_file
 use m_numeric_tools, only : get_diag, isdiagmat, arth, print_arr
 use m_copy,          only : alloc_copy
 use m_fstrings,      only : sjoin, itoa, strcat
 use m_krank,         only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_crystal,       only : crystal_t
 use m_dtset,         only : dataset_type
 use m_fftcore,       only : get_kg, sphereboundary, ngfft_seq, getng, print_ngfft !, kgindex
 use m_fft,           only : fft_ug, fft_ur, fftbox_plan3_t
 use m_fft_mesh,      only : times_eikr !, times_eigr, ig2gfft, get_gftt, calc_ceikr, calc_eigr
 use m_kpts,          only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_sort, kpts_pack_in_stars
 use m_wfk,           only : wfk_read_ebands
 use m_wfd,           only : wfd_t, wfd_init, wave_t
 use m_pawtab,        only : pawtab_type

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
!!  gvectors and tables used for zero padded FFTs.
!!
!!  TODO: Address nspinor = 2 and scalapack distribution as MPI proc can have both spinors in memory
!!  In other words we should store the first/last index in gvec for each spinor
!!
!! SOURCE

 type,public :: desc_t

   integer :: istwfk = 1
   ! Storage mode for this k point.

   integer :: npw = -1
   ! Number of plane-waves for this k-point.

   !integer :: npwsp = -1

   integer,allocatable :: gvec(:,:)
   ! gvec(3, npw)
   ! G vectors in reduced coordinates.
   ! Note that this array is global i.e. it is not MPI-distributed in G-space, only in k-space.

   integer,allocatable :: gbound(:,:)
   ! gbound(2 * mgfft + 8, 2)
   ! sphere boundary info for zero-padded FFT

 contains

   procedure :: copy => desc_copy
   ! Copy object.

   procedure :: free => desc_free
   ! Free memory.

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

   type(desc_t),allocatable :: green_desc_kibz(:)
   ! (nkibz)
   ! Descriptor for Green's functions

   type(desc_t),allocatable :: chi_desc_qibz(:)
   ! (nqibz)
   ! Descriptor for chi

   !type(desc_t),allocatable :: sigma_desc_k(:)
   ! (nkibz)
   ! Descriptor for self-energy

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

   type(datafiles_type), pointer :: dtfil => null()

   type(crystal_t), pointer  :: cryst => null()

   type(ebands_t), pointer  :: ks_ebands => null()
   type(ebands_t) :: qp_ebands

   type(pseudopotential_type), pointer :: psps => null()

   type(pawtab_type), pointer :: pawtab(:) => null()

   type(mpi_type),pointer :: mpi_enreg => null()

   type(processor_scalapack) :: g_slkproc

   type(matrix_scalapack),allocatable :: gt_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)
   ! Occupied/Empty Green's function G_k(g, g')

   type(matrix_scalapack),allocatable :: chi_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! chi_q(g, g')

   !character(len=10) :: chi_space = "itau"

   type(matrix_scalapack),allocatable :: wc_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! Correlated screened Coulomb interaction (does not depend on the spin, though)

   !character(len=10) :: em1c_space = "itau"

   type(matrix_scalapack),allocatable :: sigc_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)

   !character(len=10) :: sigc_space = "itau"

   !character(len=fnlen) :: wfk0_path = ABI_NOFILE
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

 contains

   procedure :: init => gwr_init

   procedure :: get_green_gpr => gwr_get_green_gpr

   procedure :: cos_transform  => gwr_cos_transform

   procedure :: free => gwr_free
   ! Free memory.

   procedure :: print => gwr_print
   ! Print info on the object.

   procedure :: print_trace => gwr_print_trace

   procedure :: build_gtau_from_wfk => gwr_build_gtau_from_wfk
   ! Build the Green's function in imaginary time for k-points in the IBZ from the WFK file

   procedure :: rotate_gt => gwr_rotate_gt
   procedure :: build_chi => gwr_build_chi
   procedure :: build_wc => gwr_build_wc
   procedure :: build_sigmac => gwr_build_sigmac
   procedure :: rpa_energy => gwr_rpa_energy
   procedure :: run_g0w0 => gwr_run_g0w0

 end type gwr_t
!!***

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
!! PARENTS
!!
!! CHILDREN
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
 integer :: my_nshiftq, ik_ibz, ik_bz, iq_bz, iq_ibz, npw_, ncid ! ig1, ig2,
 integer :: comm_cart, me_cart, ierr, all_nproc, my_rank
 real(dp) :: ecut_eff
 real(dp) :: cpu, wall, gflops
 character(len=5000) :: msg
 logical :: reorder
 type(desc_t),pointer :: desc_k, desc_q
 type(krank_t) :: qrank, krank_ibz
!arrays
 integer :: qptrlatt(3,3), dims(ndims)
 integer,allocatable :: kbz2ibz(:,:), qbz2ibz(:,:), gvec_(:,:)
 !integer,allocatable :: kibz2bz(:), qibz2bz(:), qglob2bz(:,:), kglob2bz(:,:)
 integer,allocatable :: count_ibz(:)
 real(dp) :: my_shiftq(3,1), kk_ibz(3), kk_bz(3), qq_bz(3), qq_ibz(3)
 real(dp),allocatable :: wtk(:), kibz(:,:)
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
 ! TODO Make sure fermie is inside the gap if semiconductor. perhaps this shoud be delegated to gwr_driver
 gwr%kibz => ks_ebands%kptns
 gwr%wtk => ks_ebands%wtk
 gwr%mpi_enreg => mpi_enreg

 ! Initialize qp_ebands with KS values.
 call ebands_copy(ks_ebands, gwr%qp_ebands)

 gwr%comm = comm
 gwr%nspinor = dtset%nspinor
 gwr%nsppol = dtset%nsppol
 gwr%use_supercell = .True.

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
 call krank_ibz%free()

 !call get_ibz2bz(gwr%nkibz, gwr%nkbz, kbz2ibz, kibz2bz, ierr)
 !ABI_CHECK(ierr == 0, "Something wrong in symmetry tables for k-points")

 ABI_FREE(kibz)

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

 ! ====================
 ! Setup tau/omega mesh
 ! ====================
 gwr%ntau = dtset%gwr_ntau
 ABI_MALLOC(gwr%tau_mesh, (gwr%ntau))
 ABI_MALLOC(gwr%iw_mesh, (gwr%ntau))
 gwr%tau_mesh = arth(zero, one, gwr%ntau)
 gwr%iw_mesh = arth(zero, one, gwr%ntau)

 ABI_MALLOC(gwr%t2w_cos_wgs, (gwr%ntau, gwr%ntau))
 ABI_MALLOC(gwr%t2w_sin_wgs, (gwr%ntau, gwr%ntau))
 gwr%t2w_cos_wgs = one
 gwr%t2w_sin_wgs = one

 ! ==========================
 ! Setup k-points in Sigma_nk
 ! ==========================
 ! TODO:


 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================

 if (any(dtset%gwr_np_gtks /= 0)) then
   ! Use parameters from input file.
   gwr%g_comm%nproc    = dtset%gwr_np_gtks(1)
   gwr%tau_comm%nproc  = dtset%gwr_np_gtks(2)
   gwr%kpt_comm%nproc  = dtset%gwr_np_gtks(3)
   gwr%spin_comm%nproc = dtset%gwr_np_gtks(4)
   !gwr%my_npert = natom3 / gwr%tau_comm%nproc
   !ABI_CHECK(gwr%my_npert > 0, "pert_comm_nproc cannot be greater than 3 * natom.")
   !ABI_CHECK(mod(natom3, gwr%tau_comm%nproc) == 0, "pert_comm_nproc must divide 3 * natom.")
 else
   ! Automatic grid generation.
   ! TODO: for the time being use all procs for tau.
   gwr%tau_comm%nproc = all_nproc
   !gwr%g_comm%nproc = all_nproc
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

 call xmpi_comm_free(comm_cart)
#endif

 ! Block-distribute dimensions and allocate redirection table: local index --> global index.
 call xmpi_split_block(gwr%ntau, gwr%tau_comm%value, gwr%my_ntau, gwr%my_itaus)
 ABI_CHECK(gwr%my_ntau > 0, "my_ntau == 0, decrease number of procs for tau level")

 call xmpi_split_block(gwr%nsppol, gwr%spin_comm%value, gwr%my_nspins, gwr%my_spins)
 ABI_CHECK(gwr%my_nspins > 0, "my_nspins == 0, decrease number of procs for spin level")

 ! Distribute k-points in the full BZ, transfer symmetry tables.
 ! Finally, find the number of IBZ points to be stored in memory by this MPI rank.

 call xmpi_split_block(gwr%nkbz, gwr%kpt_comm%value, gwr%my_nkbz, gwr%my_kbz_inds)
 ABI_CHECK(gwr%my_nkbz > 0, "my_nkbz == 0, decrease number of procs for k-point level")


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

 ! =========================================
 ! Find FFT mesh and max number of g-vectors
 ! =========================================

 ecut_eff = dtset%ecut * dtset%dilatmx ** 2
 call ngfft_seq(gwr%ngfft, [0, 0, 0])

 do ik_bz=1,gwr%nkbz
   kk_bz = gwr%kbz(:, ik_bz)
   !ik_bz = gwr%my_kbz_inds(my_ikf)
   call get_kg(kk_bz, istwfk1, ecut_eff, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%boxcutmin, dtset%chksymtnons, ecut_eff, cryst%gmet, &
              kk_bz, me_fft0, gwr%mgfft, gwr%nfft, gwr%ngfft, nproc_fft1, cryst%nsym, paral_fft0, &
              cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%green_mpw = max(gwr%green_mpw, npw_)
 end do

 do iq_bz=1,gwr%nqbz
   qq_bz = gwr%qbz(:, iq_bz)
   call get_kg(qq_bz, istwfk1, dtset%ecuteps, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%boxcutmin, dtset%chksymtnons, dtset%ecuteps, cryst%gmet, &
              qq_bz, me_fft0, gwr%mgfft, gwr%nfft, gwr%ngfft, nproc_fft1, cryst%nsym, &
              paral_fft0, cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%chi_mpw = max(gwr%chi_mpw, npw_)
 end do

 if (my_rank == master) then
   call print_ngfft(gwr%ngfft, header="FFT mesh for GWR", unit=std_out)
   call print_ngfft(gwr%ngfft, header="FFT mesh for GWR", unit=ab_out)
 end if

 ! Now we know the value of ngfft. Setup tables for zero-padded FFTs.
 ! Build FFT descriptors for Green's functions and chi.
 ! and setup tables for zero-padded FFTs.
 ABI_MALLOC(gwr%green_desc_kibz, (gwr%nkibz))

 do my_iki=1,gwr%my_nkibz
   ik_ibz = gwr%my_kibz_inds(my_iki)
   kk_ibz = gwr%kibz(:, ik_ibz)
   desc_k => gwr%green_desc_kibz(ik_ibz)
   desc_k%istwfk = istwfk1
   call get_kg(kk_ibz, desc_k%istwfk, ecut_eff, gwr%cryst%gmet, desc_k%npw, desc_k%gvec)
   ABI_MALLOC(desc_k%gbound, (2 * gwr%mgfft + 8, 2))
   call sphereboundary(desc_k%gbound, desc_k%istwfk, desc_k%gvec, gwr%mgfft, desc_k%npw)
 end do

 ABI_MALLOC(gwr%chi_desc_qibz, (gwr%nqibz))

 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   desc_q => gwr%chi_desc_qibz(iq_ibz)
   desc_q%istwfk = istwfk1
   call get_kg(qq_ibz, desc_q%istwfk, dtset%ecuteps, gwr%cryst%gmet, desc_q%npw, desc_q%gvec)
   ABI_MALLOC(desc_q%gbound, (2 * gwr%mgfft + 8, 2))
   call sphereboundary(desc_q%gbound, desc_q%istwfk, desc_q%gvec, gwr%mgfft, desc_q%npw)
 end do

 ! =============================================================
 ! Allocate scalapack arrays for G_kibz(g,g') and Chi_qibz(g,g')
 ! =============================================================
 ABI_MALLOC(gwr%gt_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))
 ABI_MALLOC(gwr%chi_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))

 call init_scalapack(gwr%g_slkproc, gwr%g_comm%value, grid_dims=[1, gwr%g_comm%nproc])

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_it=1,gwr%my_ntau
     itau = gwr%my_itaus(my_it)

     ! Allocate G_k(g,g') for k in my IBZ, MPI distributed over g' in blocks
     do my_iki=1,gwr%my_nkibz
       ik_ibz = gwr%my_kibz_inds(my_iki)
       npwsp = gwr%green_desc_kibz(ik_ibz)%npw * gwr%nspinor
       col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
       call gwr%gt_kibz(1, ik_ibz, itau, spin)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
       call gwr%gt_kibz(2, ik_ibz, itau, spin)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
       !call gwr%gt_kibz(1, ik_ibz, itau, spin)%print(header=sjoin("gt_kibz(1) for ik_ibz:", itoa(ik_ibz)))
     end do

    ! Allocate chi_q(g,g') for q in my IBZ, MPI distributed over g' in blocks
    do my_iqi=1,gwr%my_nqibz
      iq_ibz = gwr%my_qibz_inds(my_iqi)
      npwsp = gwr%chi_desc_qibz(iq_ibz)%npw * gwr%nspinor
      col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
      call gwr%chi_qibz(iq_ibz, itau, spin)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
      !call gwr%chi_qibz(iq_ibz, itau, spin)%print(header=sjoin("chi_qibz for iq_ibz:", itoa(iq_ibz)))
    end do

   end do ! my_it
 end do ! my_is

 ! Create netcdf file to store results.
 gwr%gwrnc_path = strcat(dtfil%filnam_ds(4), "_GWR.nc")
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, gwr%gwrnc_path, xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ks_ebands, ncid))
   !NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_close(ncid))
 end if

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
!! PARENTS
!!
!! CHILDREN
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
 ABI_SFREE(gwr%tau_mesh)
 ABI_SFREE(gwr%iw_mesh)
 ABI_SFREE(gwr%t2w_cos_wgs)
 ABI_SFREE(gwr%t2w_sin_wgs)

 call ebands_free(gwr%qp_ebands)

 ! Free descriptors
 if (allocated(gwr%green_desc_kibz)) then
   call desc_array_free(gwr%green_desc_kibz)
   ABI_FREE(gwr%green_desc_kibz)
 end if
 if (allocated(gwr%chi_desc_qibz)) then
   call desc_array_free(gwr%chi_desc_qibz)
   ABI_FREE(gwr%chi_desc_qibz)
 end if

 ! Free scalapack matrices
 if (allocated(gwr%gt_kibz)) then
   call slk_array_free(gwr%gt_kibz)
   ABI_FREE(gwr%gt_kibz)
 end if
 if (allocated(gwr%chi_qibz)) then
   call slk_array_free(gwr%chi_qibz)
   ABI_FREE(gwr%chi_qibz)
 end if
 if (allocated(gwr%wc_qibz)) then
   call slk_array_free(gwr%wc_qibz)
   ABI_FREE(gwr%wc_qibz)
 end if
 if (allocated(gwr%sigc_kibz)) then
   call slk_array_free(gwr%sigc_kibz)
   ABI_FREE(gwr%sigc_kibz)
 end if

 ! Free MPI communicators
 call gwr%spin_comm%free()
 call gwr%g_comm%free()
 call gwr%tau_comm%free()
 call gwr%kpt_comm%free()
 call gwr%gtau_comm%free()

 call end_scalapack(gwr%g_slkproc)

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
 integer :: mband, nkibz, nsppol, my_is, my_it, my_iki, spin, ik_ibz
 integer :: my_nb, band, ib, npw_k, istwf_k, itau, ioe, ig1, ig2, il_g1, il_g2
 real(dp) :: f_nk, eig_nk, tau, ef, fact
 real(dp) :: cpu, wall, gflops
 character(len=5000) :: msg
 type(wfd_t),target :: wfd
 type(wave_t),pointer :: wave
 type(ebands_t) :: ks_ebands
 type(hdr_type) :: wfk0_hdr
 type(dataset_type), pointer :: dtset
 type(desc_t),pointer :: desc_k
 type(matrix_scalapack),pointer :: gt
!arrays
 integer,allocatable :: nband(:,:), wfd_istwfk(:), my_bands(:)
 real(dp),allocatable :: cgwork(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 dtset => gwr%dtset

 ks_ebands = wfk_read_ebands(wfk0_path, gwr%comm, out_hdr=wfk0_hdr)
 call wfk0_hdr%vs_dtset(dtset)
 ! TODO: More consistency checks e.g. nkibz,...

 !cryst = wfk0_hdr%get_crystal()
 !call cryst%print(header="crystal structure from WFK file")

 ! Initialize the wave function descriptor.
 ! Only wavefunctions for the symmetrical imagine of the k wavevectors
 ! treated by this MPI rank are stored.

 nkibz = ks_ebands%nkpt; nsppol = ks_ebands%nsppol; mband = ks_ebands%mband

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
   dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ks_ebands%kptns, gwr%ngfft, &
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, gwr%comm)

 call wfd%print(header="Wavefunctions for GWR calculation")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)
 call ebands_free(ks_ebands)
 call wfk0_hdr%free()

 ! Read KS wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ! ==================================
 ! Build Green's functions in g-space
 ! ==================================
 ! NB: G_k is constructed for k in the IBZ, then we rotate the k-point to obtain G_k in the BZ.
 ! TODO:
 ! 1) Make sure that gvec in gwr and wfd agree with each other.
 ! 2) May implement the following trick to accelerate convergence wrt nband:
 !     - Init nb states with random numbers and orthogonalize wrt the nband states found in the WFK file
 !     - Compute <i|H|j> in the subspace spanned by the perp states
 !     - Diagonalize H in the subspace to get variational approximation to the KS states
 !     - Add these approximated eigenstats to G.

 ef = ks_ebands%fermie

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)

   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     !ik_bz = gwr%kibz2bz(ik_ibz)
     desc_k => gwr%green_desc_kibz(ik_ibz)

     ABI_CHECK(npw_k == desc_k%npw, "npw_k != desc_k%npw")
     ABI_CHECK(all(wfd%kdata(ik_ibz)%kg_k == desc_k%gvec), "kg_k != desc_k%gvec")

     ABI_MALLOC(cgwork, (2, npw_k * wfd%nspinor))

     do ib=1,my_nb
       band = my_bands(ib)
       f_nk = gwr%ks_ebands%occ(band, ik_ibz, spin)
       eig_nk = gwr%ks_ebands%eig(band, ik_ibz, spin) - ef
       if (eig_nk < tol6) then
         ioe = 1
       else if (eig_nk > tol6) then
         ioe = 2
       else
         ABI_WARNING("Metallic system of semiconductor with Fermi level inside bands!!!!")
       end if

       call wfd%copy_cg(band, ik_ibz, spin, cgwork)
       ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)
       !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk_ibz, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

       do my_it=1,gwr%my_ntau
         itau = gwr%my_itaus(my_it)
         tau = gwr%tau_mesh(itau)
         fact = exp(-tau * eig_nk)
         gt => gwr%gt_kibz(ioe, ik_ibz, itau, spin)
         !call wave%outer_prod()
         do il_g2=1, gt%sizeb_local(2)
           ig2 = gt%loc2gcol(il_g2)
           !ig2 = gt%glob_col(il_g2)
           do il_g1=1, gt%sizeb_local(1)
             ! TODO: spinor and use BLAS2 but take into account scalapack distribution
             ig1 = gt%loc2grow(il_g1)
             !cgwork(ig1) x cgwork(ig2) *
             !gt%buffer_cplx(il_g1, il_g2) = gt%buffer_cplx(il_g1, il_g2) + &
             !  fact * wave%ug(ig1) * wave%ug(ig2)
           end do
         end do

       end do ! my_it

     end do ! ib
     ABI_FREE(cgwork)

   end do ! my_ikf
 end do ! my_is

 call wfd%free()
 ABI_FREE(my_bands)

 call cwtime_report(" gwr_build_gtau_from_wfk:", cpu, wall, gflops)

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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gwr_rotate_gt(gwr, my_ikf, my_it, my_is, desc_kbz, gt_kbz)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 integer,intent(in) :: my_ikf, my_it, my_is
 type(desc_t),intent(out) :: desc_kbz
 type(matrix_scalapack),intent(inout) :: gt_kbz(2)

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, ioe, spin, itau
 integer :: ik_ibz, isym_k, trev_k, g0_k(3)
 logical :: isirr_k

! *************************************************************************

 spin = gwr%my_spins(my_is)
 itau = gwr%my_itaus(my_it)

 ik_ibz = gwr%my_kbz2ibz(1, my_ikf); isym_k = gwr%my_kbz2ibz(2, my_ikf)
 trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 !tsign = merge(1, -1, trev_k == 0)

 ! Copy descriptor from IBZ, rotate gvec and recompute gbound.
 call gwr%green_desc_kibz(ik_ibz)%copy(desc_kbz)

 if (isirr_k) then
   do ioe=1,2
     call gwr%gt_kibz(ioe, ik_ibz, itau, spin)%copy(gt_kbz(ioe))
     call gt_kbz(ioe) %print(header=sjoin("isirr_k with gt_kibz:", itoa(ik_ibz)))
   end do
   return
 end if

 ! TODO: Handle TR
 do ig1=1,desc_kbz%npw
   desc_kbz%gvec(:,ig1) = matmul(gwr%cryst%symrec(:,:,isym_k), gwr%green_desc_kibz(ik_ibz)%gvec(:,ig1)) ! - g0_k
 end do

 call sphereboundary(desc_kbz%gbound, desc_kbz%istwfk, desc_kbz%gvec, gwr%mgfft, desc_kbz%npw)

 ! Get G_k with k in the BZ.
 do ioe=1,2
   call gwr%gt_kibz(ioe, ik_ibz, itau, spin)%copy(gt_kbz(ioe))
   do il_g2=1,gt_kbz(ioe)%sizeb_local(2)
     ig2 = gt_kbz(ioe)%loc2gcol(il_g2)
     do il_g1=1, gt_kbz(ioe)%sizeb_local(1)
       ig1 = gt_kbz(ioe)%loc2grow(il_g1)
       !call gt_kbz(ioe)%loc2glob(il_g1, il_g2, ig1, ig2)
       !gt_kbz(ioe)%buffer_cplx(il_g1, il_g2) = gt_kbz(ioe)%buffer_cplx(il_g1, il_g2) ! * ??
     end do
   end do
   !call gt_kbz(ioe)%print(header=sjoin("not isirr_k with gt_kibz:", itoa(ik_ibz)))
 end do

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
!!  for each k in the BZ treated by me for given spin and tau.
!!
!!  1) FFT Transform the first index (local): G(g,g',it) --> G(r,g',it)
!!  2) MPI transposition: G(r,g',it) --> G^*(g',r,it)
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

subroutine gwr_get_green_gpr(gwr, my_it, my_is, desc_kbz, gt_gpr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 integer,intent(in) :: my_it, my_is
 type(desc_t),target,intent(out) :: desc_kbz(gwr%my_nkbz)
 type(matrix_scalapack),intent(inout) :: gt_gpr(2, gwr%my_nkbz)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1
 integer :: my_ikf, ig2, ioe, npwsp, col_bsize !, ncol_glob,
 !integer :: ik_ibz, isym_k, trev_k, g0_k(3)
 real(dp) :: cpu, wall, gflops
 !logical :: isirr_k
 type(matrix_scalapack), pointer :: ggp
 type(matrix_scalapack) :: rgp
 type(matrix_scalapack),target :: gt_kbz(2)
 type(desc_t),pointer :: desc_k

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 !call wrtout(std_out, " in gwr_green_gpr")
 !call gwr%print(header="GWR_GREEN_GPR")
 !print *, "gwr%g_slkproc%grid%dims", gwr%g_slkproc%grid%dims

 do my_ikf=1,gwr%my_nkbz
   !ik_ibz = gwr%my_kbz2ibz(1, my_ikf); isym_k = gwr%my_kbz2ibz(2, my_ikf)
   !trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
   !isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !ik_bz = gwr%my_kbz_inds(my_ikf)

   ! Get G_k in the BZ.
   call gwr%rotate_gt(my_ikf, my_it, my_is, desc_kbz(my_ikf), gt_kbz)
   desc_k => desc_kbz(my_ikf)

   do ioe=1,2
     ggp => gt_kbz(ioe)

     ! Allocate rgp scalapack matrix to store G(r, g')
     npwsp = desc_k%npw * gwr%nspinor
     col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
     !print *, "desc_k%npw", desc_k%npw
     !print *, "npw", size(desc_k%gvec, dim=2)
     !print *, "desc_k%istwfk", desc_k%istwfk
     call rgp%init(gwr%nfft * gwr%nspinor, npwsp, gwr%g_slkproc, desc_k%istwfk, &
                   size_blocs=[gwr%nfft * gwr%nspinor, col_bsize])

     ABI_CHECK(size(ggp%buffer_cplx, dim=2) == size(rgp%buffer_cplx, dim=2), "len2")

     ! FFT. Results stored in rgp
     do ig2=1,ggp%sizeb_local(2)

       ABI_CHECK(size(ggp%buffer_cplx(:, ig2)) == desc_k%npw, "npw")
       ABI_CHECK(size(rgp%buffer_cplx(:, ig2)) == gwr%nfft * gwr%nspinor, "gwr%nfft * gwr%nspinor")

       call fft_ug(desc_k%npw, gwr%nfft, gwr%nspinor, ndat1, &
                   gwr%mgfft, gwr%ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
                   ggp%buffer_cplx(:, ig2), rgp%buffer_cplx(:, ig2))

       ! Multiply by e^{ik.r}
       !call times_eikr(desc%point, desc%ngfft, desc%nfft, ndat1, rgp%buffer_cplx(:, ig2)
     end do ! ig2

     ! Transpose: G(r, g') -> G(g', r) and take complex conjugate.
     call rgp%ptrans("C", gt_gpr(ioe, my_ikf))
     call rgp%free()
   end do ! ioe

   call slk_array_free(gt_kbz)
 end do ! my_ikf

 call cwtime_report(" gwr_get_green_gpr:", cpu, wall, gflops)

end subroutine gwr_get_green_gpr
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
 integer :: my_iqi, my_is, ig1, ig2, my_it, it, ierr, iq_ibz, itau, spin, it0, idx
 real(dp) :: cpu, wall, gflops
 logical :: sum_spins_
 type(desc_t),pointer :: desc_q
!arrays
 real(dp),allocatable :: my_weights(:,:)
 complex(dp):: cwork_t(gwr%my_ntau), cwork_wglb(gwr%ntau)
 type(matrix_scalapack), pointer :: mats(:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 sum_spins_ = .False.; if (present(sum_spins)) sum_spins_ = sum_spins

 ABI_MALLOC(my_weights, (gwr%ntau, gwr%my_ntau))

 ! Extract my weights from the global array according to mode.
 idx = -1
 if (mode == "w2t") idx = 1
 if (mode == "t2w") idx = 2
 ABI_CHECK(idx /= -1, sjoin("Wrong mode:", mode))

 do my_it=1,gwr%my_ntau
   itau = gwr%my_itaus(my_it)
   !tau = gwr%tau_mesh(itau)
   !my_weights(:, my_it) = cos(omega * tau) * gwr%global_weighs(:, itau, idx)
 end do

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     desc_q => gwr%chi_desc_qibz(iq_ibz)

     mats => null()
     if (what == "chi") mats => gwr%chi_qibz(iq_ibz, :, spin)
     if (what =="wc")   mats => gwr%wc_qibz(iq_ibz, :, spin)
     ABI_CHECK(associated(mats), sjoin("Invalid value for what:", what))

     it0 = gwr%my_itaus(1)

     do ig2=1,mats(it0)%sizeb_local(2)
       do ig1=1,mats(it0)%sizeb_local(1)
         ! Extract matrix elements as a function of tau
         ! TODO: Here we can block over ig1 and call zgemm
         ! to reduced the number of MPI communications.
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           cwork_t(my_it) = mats(itau)%buffer_cplx(ig1, ig2)
         end do

         do it=1,gwr%ntau
           cwork_wglb(it) = dot_product(my_weights(it, :), cwork_t)
         end do

         call xmpi_sum(cwork_wglb, gwr%tau_comm%value, ierr)

         ! update values for this (g1, g2), now we have a slice of gg' in frequency space.
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           mats(itau)%buffer_cplx(ig1, ig2) = cwork_wglb(itau)
         end do

       end do ! ig1
     end do ! ig2

   end do ! my_iqi
 end do ! my_is

 !gwr%chi_space = "itau"
 !gwr%chi_space = "iomega"

 if (gwr%nsppol == 2 .and. sum_spins_ .and. gwr%spin_comm%nproc > 1) then
   do my_iqi=1,gwr%my_nqibz
      iq_ibz = gwr%my_qibz_inds(my_iqi)
      desc_q => gwr%chi_desc_qibz(iq_ibz)

      do my_is=1,gwr%my_nspins
        spin = gwr%my_spins(my_is)
        mats => null()
        if (what == "chi") mats => gwr%chi_qibz(iq_ibz,:,spin)
        !if (what =="wc")   mats => gwr%wc_qibz(iq_ibz, :, spin)
        ABI_CHECK(associated(mats), sjoin("Invalid value for what:", what))

        do my_it=1,gwr%my_ntau
          itau = gwr%my_itaus(my_it)
          call xmpi_sum(mats(itau)%buffer_cplx, gwr%spin_comm%value, ierr)
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
!!  won't see the allocation done by the compiler.
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

subroutine desc_copy(in_desc, new_desc)

!Arguments ------------------------------------
 class(desc_t),intent(in) :: in_desc
 class(desc_t),intent(out) :: new_desc

! *************************************************************************

 new_desc%istwfk = in_desc%istwfk
 new_desc%npw = in_desc%npw

 call alloc_copy(in_desc%gvec, new_desc%gvec)
 call alloc_copy(in_desc%gbound, new_desc%gbound)

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
 class(desc_t),intent(inout) :: desc

! *************************************************************************

 ABI_SFREE(desc%gvec)
 ABI_SFREE(desc%gbound)

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
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gwr_print(gwr, header, unit)

!Arguments ------------------------------------
 class(gwr_t),target,intent(in) :: gwr
 character(len=*),optional,intent(in) :: header
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
!scalars
 integer :: my_is, spin, my_it, itau, my_iki, ik_ibz, unt
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
 call ydoc%add_int("chi_mpw", gwr%chi_mpw)
 call ydoc%add_int1d("ngfft", gwr%ngfft(1:8))
 call ydoc%add_int1d("P gwr_np_gtks", [gwr%g_comm%nproc, gwr%tau_comm%nproc, gwr%kpt_comm%nproc, gwr%spin_comm%nproc])
 call ydoc%write_and_free(unt)

 if (gwr%dtset%prtvol > 10) then
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       do my_iki=1,gwr%my_nkibz
         ik_ibz = gwr%my_kibz_inds(my_iki)
         call gwr%gt_kibz(1, ik_ibz, itau, spin)%print(header=sjoin("gt_kibz(1) for ik_ibz:", itoa(ik_ibz)))
         call gwr%gt_kibz(2, ik_ibz, itau, spin)%print(header=sjoin("gt_kibz(2) for ik_ibz:", itoa(ik_ibz)))
       end do
     end do
   end do
 end if

end subroutine gwr_print
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_chi
!! NAME
!!  gwr_build_chi
!!
!! FUNCTION
!!  Compute polarizability.
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
 integer,parameter :: ndat1 = 1, istwfk1 = 1
 integer :: my_is, my_it, my_ikf, ig, my_ir, my_nr, npwsp, ncol_glob, col_bsize, my_iqi !, my_iki ! my_iqf,
 integer :: sc_nfft, spin, ik_ibz, ik_bz, iq_ibz, ierr, ioe, itau, ig2 ! ig1
 real(dp) :: cpu_tau, wall_tau, gflops_tau, cpu_all, wall_all, gflops_all
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

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 !call gwr%print(header="In build_chi")

 if (gwr%use_supercell) then
   ! ============================
   ! GWr algorithm with supercell
   ! ============================
   call wrtout(std_out, " Building chi0 = G0G0 with FFTs in the supercell ...")

   ! Set FFT mesh in the supercell
   ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
   ! Perhaps the safest approach would be to generate the plan on the fly.
   sc_ngfft = gwr%ngfft
   sc_ngfft(1:3) = gwr%ngkpt * gwr%ngfft(1:3)
   sc_ngfft(4:6) = sc_ngfft(1:3)
   sc_nfft = product(sc_ngfft(1:3))
   !sc_augsize = product(sc_ngfft(4:6))

   ABI_CALLOC(gt_scbox, (sc_nfft * gwr%nspinor * ndat1, 2))
   ABI_CALLOC(chit_scbox, (sc_nfft * gwr%nspinor * ndat1))

   ! (ngfft, ndat, isign)
   call plan_gp2rp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1 * 2, +1)
   call plan_rp2gp%from_ngfft(sc_ngfft, gwr%nspinor * ndat1,     -1)

   ! The g-vectors in the supercell for G and chi.
   ABI_MALLOC(green_scg, (3, gwr%green_mpw))
   ABI_MALLOC(chi_scg, (3, gwr%chi_mpw))

   ! Allocate Scalapack arrays for chi_q(g',r) for all q in the IBZ treated by this MPI rank.
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     npwsp = gwr%chi_desc_qibz(iq_ibz)%npw * gwr%nspinor
     ncol_glob = gwr%nfft * gwr%nspinor
     col_bsize = ncol_glob / gwr%g_comm%nproc; if (mod(ncol_glob, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
     call chiq_gpr(my_iqi)%init(npwsp, gwr%nfft * gwr%nspinor, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
   end do

   ! Loop over my spins and my taus.
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
       itau = gwr%my_itaus(my_it)

       ! G_k(g,g') --> G_k(g',r) for each k in the BZ treated by me.
       call gwr%get_green_gpr(my_it, my_is, desc_kbz, gt_gpr)

       ! Loop over the second index r (MPI-distributed in g_comm).
       my_nr = gt_gpr(1, 1)%sizeb_local(2)
       do my_ir=1,my_nr

         ! Insert G_k(g',r) in G'-space in the supercell FFT box for fixed r.
         ! Note that we need to take the union of (k, g') for k in the BZ.
         gt_scbox = zero
         do my_ikf=1,gwr%my_nkbz
           ik_bz = gwr%my_kbz_inds(my_ikf)
           ik_ibz = gwr%my_kbz2ibz(1, my_ikf)

           desc_k => desc_kbz(my_ikf)
           do ig=1,desc_k%npw
             green_scg(:,ig) = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ig)
           end do

           do ioe=1,2
             call gsph2box(sc_ngfft, desc_k%npw, gwr%nspinor * ndat1, green_scg, &
                           gt_gpr(ioe, my_ikf)%buffer_cplx(:, my_ir), &  ! in
                           gt_scbox(:, ioe))                             ! inout
           end do
         end do ! my_ikf

         ! TODO: Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
         call xmpi_sum(gt_scbox, gwr%kpt_comm%value, ierr)

         ! G(G',r) --> G(R',r)
         call plan_gp2rp%execute_ip_dpc(gt_scbox)

         ! Compute chi(R',r) for this r
         chit_scbox = gt_scbox(:, 1) * gt_scbox(:, 2)

         ! Go back to G'= q + g' space immediately.
         call plan_rp2gp%execute(chit_scbox)

         ! Extract chi_q(g', r) on the g-sphere from the FFT box in the supercell.
         ! not necessarly equal to the one used for the Green's function.
         ! NB: For the time being I'm assuming k-mesh == q-mesh.
         ! Only q-points in the IBZ are considered.
         do my_iqi=1,gwr%my_nqibz
           iq_ibz = gwr%my_qibz_inds(my_iqi)
           qq_ibz = gwr%qibz(:, iq_ibz)
           desc_q => gwr%chi_desc_qibz(iq_ibz)

           do ig=1,desc_q%npw
             chi_scg(:,ig) = nint(qq_ibz * gwr%ngqpt) + gwr%ngqpt * desc_q%gvec(:,ig)
           end do

           call box2gsph(sc_ngfft, desc_q%npw, gwr%nspinor * ndat1, chi_scg, &
                         chit_scbox, &                            ! in
                         chiq_gpr(my_iqi)%buffer_cplx(:, my_ir))  ! out
         end do ! my_iqi
       end do ! my_ir

       ! Free descriptors and scalapack matrices in kBZ.
       call desc_array_free(desc_kbz)
       call slk_array_free(gt_gpr)

       ! chi_q(g',r) --> chi_q(r,g') --> chi_q(g,g') stored in gwr%chi_qibz
       ! only for the q-points in the IBZ treated by this MPI proc.
       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         desc_q => gwr%chi_desc_qibz(iq_ibz)

         call chiq_gpr(my_iqi)%ptrans("C", chi_rgp)

         !call wrtout(std_out, "DEBUG")
         !call chiq_gpr(my_iqi)%print(header="chiq_gpr(my_iqi)")
         !call gwr%chi_qibz(iq_ibz, itau, spin)%print(header="chi_qibz")
         !call chi_rgp%print(header="chi_rgp")
         !call wrtout(std_out, "END DEBUG")

         ABI_CHECK(size(gwr%chi_qibz(iq_ibz, itau, spin)%buffer_cplx, dim=2) == size(chi_rgp%buffer_cplx, dim=2), "len2")

         ! FFT r --> g along the first dimension: chi_q(r,g') --> chi_q(g,g')
         do ig2=1,chi_rgp%sizeb_local(2)

           ABI_CHECK(size(chi_rgp%buffer_cplx(:, ig2)) == gwr%nfft * gwr%nspinor, "gwr%nfft * gwr%nspinor")
           ABI_CHECK(size(gwr%chi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2)) == desc_q%npw, "npw")

           call fft_ur(desc_q%npw, gwr%nfft, gwr%nspinor, ndat1, gwr%mgfft, gwr%ngfft, &
                       istwfk1, desc_q%gvec, desc_q%gbound, &
                       chi_rgp%buffer_cplx(:, ig2),         &                  ! ur(in)
                       gwr%chi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2))   ! ug(out)
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
   call slk_array_free(chiq_gpr)

  else
    ! ===================================================================
    ! Mixed-space algorithm in the unit cell with convolutions in k-space
    ! ===================================================================
    ! Each MPI proc has all the G_k in the IBZ thus we only need to rotate to get the BZ on the fly.
    do my_is=1,gwr%my_nspins
      spin = gwr%my_spins(my_is)
      do my_it=1,gwr%my_ntau
        call cwtime(cpu_tau, wall_tau, gflops_tau, "start")
        itau = gwr%my_itaus(my_it)

        !need_kibz(gwr%nkibz)
        !got_kibz(gwr%nkibz)
        !call gwr%green_redistribute(my_it, my_is, need_kibz, got_kibz)
        !call gwr%green_free(my_it, my_is, got_kibz)

        do my_ikf=1,gwr%my_nkbz
          !ik_ibz = gwr%my_kbz2ibz(1, my_ikf); isym_k = gwr%my_kbz2ibz(2, my_ikf)
          !trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
          !isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
          ik_bz = gwr%my_kbz_inds(my_ikf)
          kk_bz = gwr%kbz(:, ik_bz)

          ! Use symmetries to get G_k from the corresponding image in the IBZ.
          ! G_k(g,g') --> G_k(r,r')
          !call gwr%gk_rrp(my_ikf, my_it, my_is, gk_rrp)

          do my_iqi=1,gwr%my_nqibz
            iq_ibz = gwr%my_qibz_inds(my_iqi)
            qq_ibz = gwr%qibz(:, iq_ibz)
            desc_q => gwr%chi_desc_qibz(iq_ibz)
            !iq_bz = gwr%my_qbz_inds(my_iqf)
            !qq_bz = gwr%qbz(:, iq_bz)

            kpq_bz = kk_bz + qq_ibz
            !ikpq_ibz = ??

            ! Use symmetries to get G_kq from the corresponding image in the IBZ.
            ! G_kq(g,g') --> G_kq(r,r')
            !call gwr%gk_rrp(my_ikqf, my_it, my_is, gkq_rrp)

            !desc_k = gwr%green_desc_kibz(ik_ibz)%rotate(symtab_k)
            !desc_kpq = gwr%green_desc_kibz(ikpq_ibz)%rotate(symtab_kpq)

            !gwr%gt_kibz(1, ik_ibz, itau, spin)%buffer_cplx(:,ig2)

            !do ig2=1,ggp%sizeb_local(2)
            !  call fft_ug(desc_k%npw, gwr%nfft, gwr%nspinor, ndat1, &
            !              gwr%mgfft, gwr%ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
            !              ggp%buffer_cplx(:, ig2),
            !              rgp%buffer_cplx(:, ig2))
            !end do ! ig2
            !call rgp%ptrans("C", gpr)
            !call rgp%free()

            !do ir1=1,ggp%sizeb_local(2)
            !  call fft_ug(desc_k%npw, gwr%nfft, gwr%nspinor, ndat1, &
            !              gwr%mgfft, gwr%ngfft, desc_k%istwfk, desc_k%gvec, desc_k%gbound, &
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
              !!gwr%chi_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) = FFT_rg
            !end do
            !call chi_r_gp%free()
          end do ! my_ikf
        end do ! iq_ibz

       write(msg,'(2(a,i0),a)')" My itau [", my_it, "/", gwr%my_ntau, "]"
       call cwtime_report(msg, cpu_tau, wall_tau, gflops_tau)
     end do ! my_it
   end do ! spin
 end if

 ! Print trace of chiq matrices for testing purposes.
 call gwr%print_trace("chi_qibz")

 ! Transform irreducible chi from imaginary tau to imaginary omega.
 ! Also sum over spins to get total chi if collinear spin.
 call gwr%cos_transform("chi", "t2w", sum_spins=.True.)

 !call gwr%print(header="GWR before leaving build_chi")
 call cwtime_report(" gwr_build_chi:", cpu_all, wall_all, gflops_all)

end subroutine gwr_build_chi
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_print_trace
!! NAME
!!  gwr_print_trace
!!
!! FUNCTION
!!  Print traces of scalapack matrices to std_out and ab_out.
!!  NB: This is a global routine that should be called by all procs inside gwr%comm.
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

subroutine gwr_print_trace(gwr, what)

!Arguments ------------------------------------
 class(gwr_t),target,intent(in) :: gwr
 character(len=*),intent(in) :: what

!Local variables-------------------------------
 integer :: my_is, spin, my_it, itau, iq_ibz, ierr, my_iqi
 complex(dp),allocatable :: ctrace(:,:,:)

! *************************************************************************

 select case (what)
 case ("chi_qibz")
   ABI_CALLOC(ctrace, (gwr%nqibz, gwr%ntau, gwr%nsppol))

   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         ctrace(iq_ibz, itau, spin) = gwr%chi_qibz(iq_ibz, itau, spin)%get_trace()
       end do
     end do
   end do

 case default
   ABI_ERROR(sjoin("Invalid value of what:", what))
 end select

 call xmpi_sum(ctrace, gwr%comm, ierr)

 if (xmpi_comm_rank(gwr%comm) == 0) then
   call wrtout(ab_out, sjoin("Trace of", what, "matrix for testing purposes:"))
   do spin=1,gwr%nsppol
     call print_arr(ctrace(:,:,spin), unit=ab_out)
   end do
 end if

 ABI_FREE(ctrace)

end subroutine gwr_print_trace
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_wc
!! NAME
!!  gwr_build_wc
!!
!! FUNCTION
!!  Compute Wc(tau) from chi(i omega)
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

subroutine gwr_build_wc(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer :: my_iqi, my_it, il_g1, il_g2, ig1_glob, ig2_glob, my_is, iq_ibz, spin, itau
 !integer :: npwsp !, col_bsize ! ierr,
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_q, wall_q, gflops_q
 logical :: is_gamma, free_chi
 character(len=500) :: msg
 type(desc_t), pointer :: desc_q
 type(matrix_scalapack),pointer :: wc
 type(matrix_scalapack) :: tmp_mat
 !arrays
 real(dp) :: qq_ibz(3)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call wrtout(std_out, " Building correlated screening Wc ...")

 ! ================================================
 ! Allocate scalapack arrays for wct_kibz(g,g')
 ! ================================================
 ! Note that we have already summed chi over spin.
 ABI_MALLOC(gwr%wc_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))
 free_chi = .False.

 do my_iqi=1,gwr%my_nqibz
   call cwtime(cpu_q, wall_q, gflops_q, "start")
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   desc_q => gwr%chi_desc_qibz(iq_ibz)
   is_gamma = all(abs(qq_ibz) < tol12) ! Handle q --> 0 limit

   do my_is=1,gwr%my_nspins  ! This loop is needed so that procs in the spin pool can operate on their matrix.
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)

       !npwsp = gwr%chi_desc_qibz(iq_ibz)%npw * gwr%nspinor
       !col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
       !call gwr%wc_qibz(iq_ibz, itau, spin)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])

       call gwr%chi_qibz(iq_ibz, itau, spin)%copy(gwr%wc_qibz(iq_ibz, itau, spin), free=free_chi)
       wc => gwr%wc_qibz(iq_ibz, itau, spin)
       !call wc%print(header="wc")

       ! Build epsilon(q,iw) = delta_{g,g'} - v_q(g,g') chi_q(g,g,iw).
       wc%buffer_cplx = zero
       do il_g2=1,wc%sizeb_local(2)
         ig2_glob = wc%loc2gcol(il_g2)
         do il_g1=1,wc%sizeb_local(1)
           ig1_glob = wc%loc2grow(il_g1)
           !vc_qg1 = ??
           !vc_qg2 = ??
           !ceps = - wc%buffer_cplx(ig1, ig2) ! * vc_qg1 * vc_qg2
           if (ig1_glob == ig2_glob) wc%buffer_cplx(il_g1, il_g2) = one ! - ceps
         end do ! il_g1
       end do ! il_g2

       ! IMPORTANT: PZGETRF requires square block decomposition i.e., MB_A = NB_A.
       ! This means I need to redistribute the data before calling zinvert.
       call wc%change_size_blocs(tmp_mat)
       call tmp_mat%zinvert()
       call wc%take_from(tmp_mat)
       call tmp_mat%free()

       ! Build W(q, iw) = e^{-1}_q(g,g',iw) v_q(g,g')
       ! Remove bare vc
       do il_g2=1,wc%sizeb_local(2)
         ig2_glob = wc%loc2gcol(il_g2)
         do il_g1=1,wc%sizeb_local(1)
           ig1_glob = wc%loc2grow(il_g1)
           !vc_qg1 = ??
           !vc_qg2 = ??
           !wc%buffer_cplx(il_g1, il_g2) = wc%buffer_cplx(il_g1, il_g2) ! * vc_qg1 * vc_qg2
         end do ! il_g1
       end do ! il_g2

     end do  ! my_it
   end do ! my_is

   write(msg,'(2(a,i0),a)')" My iqi [", my_iqi, "/", gwr%my_nqibz, "]"
   call cwtime_report(msg, cpu_q, wall_q, gflops_q)
 end do ! my_iqi

 ! Cosine transform from iomega to itau to get Wc(i tau)
 call gwr%cos_transform("wc", "w2t")

 call cwtime_report(" gwr_build_wc:", cpu_all, wall_all, gflops_all)

end subroutine gwr_build_wc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_sigmac
!! NAME
!!  gwr_build_sigmac
!!
!! FUNCTION
!!  Compute Sigma_c(i tau)
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

subroutine gwr_build_sigmac(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
 integer :: my_is, my_it, my_iki, spin, ik_ibz !, ioe my_iw,
 integer :: my_iqi, iq_ibz, itau, col_bsize, npwsp ! ierr, ig1, ig2, ig1_glob, ig2_glob,
 real(dp) :: cpu_all, wall_all, gflops_all !, cpu_k, wall_k, gflops_k
 logical :: is_gamma
 type(desc_t), pointer :: desc_q, desc_k
 !type(matrix_scalapack),pointer :: wc
!arrays
 real(dp) :: kk_ibz(3), qq_ibz(3)
 type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz)
 type(desc_t), target :: desc_kbz(gwr%my_nkbz)

! *************************************************************************

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call wrtout(std_out, " Building Sigma_c ...")

 !if (gwr%use_supercell) then

 ! ================================================
 ! Allocate scalapack arrays for sigmac_kibz(g,g')
 ! ================================================
 ABI_MALLOC(gwr%sigc_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_it=1,gwr%my_ntau
     itau = gwr%my_itaus(my_it)

     do my_iki=1,gwr%my_nkibz
       ik_ibz = gwr%my_kibz_inds(my_iki)
       npwsp = gwr%green_desc_kibz(ik_ibz)%npw * gwr%nspinor
       col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
       call gwr%sigc_kibz(1, ik_ibz, itau, spin)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
       call gwr%sigc_kibz(2, ik_ibz, itau, spin)%init(npwsp, npwsp, gwr%g_slkproc, 1, size_blocs=[npwsp, col_bsize])
       !call gwr%sigc_kibz(1, ik_ibz, itau, spin)%print(header=sjoin("sigc_kibz(1) for ik_ibz:", itoa(ik_ibz)))
     end do

   end do ! my_it
 end do ! my_is

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     kk_ibz = gwr%kibz(:, ik_ibz)
     desc_k => gwr%green_desc_kibz(ik_ibz)

     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       !tau = gwr%tau_mesh(itau)
       ! G_k(g, g') --> G_k(g', r) for each k in the BZ treated by me.
       call gwr%get_green_gpr(my_it, my_is, desc_kbz, gt_gpr)
     end do

     ! Free descriptors and scalapack matrices in kBZ.
     call desc_array_free(desc_kbz)
     call slk_array_free(gt_gpr)

   end do ! my_ikq
 end do ! my_is

 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   desc_q => gwr%chi_desc_qibz(iq_ibz)
   is_gamma = all(abs(qq_ibz) < tol12) ! Handle q --> 0 limit

   do my_is=1,gwr%my_nspins  ! FIXME
     do my_it=1,gwr%my_ntau
     end do  ! my_it
   end do ! my_is
 end do ! my_iqi

 ! Cosine transform from iw to tau to get W(tau)
 !call gwr%cos_transform("wc", "w2t")

 call cwtime_report(" gwr_build_sigmac:", cpu_all, wall_all, gflops_all)

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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gwr_rpa_energy(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
!arrays

! *************************************************************************

 call gwr%build_chi()

 ! Compute RPA total energy with different number of PWs and extrapolate.

 ! Print results to ab_out and nc file

end subroutine gwr_rpa_energy
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_run_g0w0
!! NAME
!!  gwr_rpa_energy
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

subroutine gwr_run_g0w0(gwr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr

!Local variables-------------------------------
!scalars
!arrays

! *************************************************************************

 !call gwr%get_sigx_mels()
 call gwr%build_chi()
 call gwr%build_wc()
 !call gwr%build_sigmac()

 !gwr%dtset%prtvol = 100
 !call gwr%print(header="Before returning from run_g0w0")
 !call xmpi_barrier(gwr%comm)

end subroutine gwr_run_g0w0
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gsph2box
!! NAME
!! gsph2box
!!
!! FUNCTION
!! Insert cg_k array defined on the k-centered g-sphere with npw points inside supercell FFT box.
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
!! PARENTS
!!
!! CHILDREN
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
!! Extract cg_k array defined on the k-centered g-sphere with npw points from the supercell FFT box.
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
!! PARENTS
!!
!! CHILDREN
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

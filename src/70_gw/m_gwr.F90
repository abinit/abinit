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
 use m_wfd
 use m_hdr
 use m_ebands
 use netcdf
 use m_nctk
 use m_dtfil

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use defs_abitypes,   only : mpi_type
 use m_time,          only : cwtime, cwtime_report, sec2str
 use m_io_tools,      only : iomode_from_fname !, file_exists, is_open, open_file
 use m_numeric_tools, only : get_diag, isdiagmat
 use m_copy,          only : alloc_copy
 use m_fstrings,      only : sjoin, itoa, strcat
 use m_krank,         only : krank_t, krank_new, krank_from_kptrlatt, get_ibz2bz, star_from_ibz_idx
 use m_crystal,       only : crystal_t
 use m_dtset,         only : dataset_type
 use m_fftcore,       only : get_kg, sphereboundary, ngfft_seq, getng, print_ngfft !, kgindex
 use m_fft,           only : fft_ug, fftbox_plan3, fftbox_plan3_init
 use m_fft_mesh,      only : times_eikr !, times_eigr, ig2gfft, get_gftt, calc_ceikr, calc_eigr
 use m_kpts,          only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_sort, kpts_pack_in_stars
 use m_wfk,           only : wfk_read_ebands
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
!! SOURCE

 type,public :: desc_t

   integer :: istwfk = 1
   ! Storage mode for this k point.

   integer :: npw = -1
   ! Number of plane-waves for this k-point.

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

   procedure :: copy => desc_copy
   ! Copy object.

   procedure :: free => desc_free
   ! Free memory.

 end type desc_t
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

   type(pseudopotential_type), pointer :: psps => null()

   type(pawtab_type), pointer :: pawtab(:) => null()

   type(mpi_type),pointer :: mpi_enreg => null()

   type(matrix_scalapack),allocatable :: gt_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)
   ! Occupied/Empty Green's function G_k(g, g')

   type(matrix_scalapack),allocatable :: chit_qibz(:,:,:)
   ! (nqibz, ntau, nsppol)
   ! chi_q(g, g')

   type(matrix_scalapack),allocatable :: wct_qibz(:,:)
   ! (nqibz, ntau)
   ! Correlated screened Coulomb interaction

   type(matrix_scalapack),allocatable :: sigct_kibz(:,:,:,:)
   ! (2, nkibz, ntau, nsppol)

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

   procedure :: green_gpr => gwr_green_gpr

   procedure :: cos_transform  => gwr_cos_transform

   procedure :: free => gwr_free
   ! Free memory.

   !procedure :: print => gwr_print

   procedure :: build_gtau_from_wfk => gwr_build_gtau_from_wfk
   ! Build the Green's function in imaginary time for k-points in the IBZ from the WFK file

   procedure :: rotate_gt => gwr_rotate_gt
   procedure :: build_chi => gwr_build_chi
   procedure :: build_wc => gwr_build_wc
   procedure :: build_sigmac => gwr_build_sigmac
   procedure :: rpa_energy => gwr_rpa_energy
   procedure :: run_g0w0 => gwr_run_g0w0

 end type gwr_t

 public :: gwr_new
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_new
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

function gwr_new(dtset, dtfil, cryst, psps, pawtab, ks_ebands, mpi_enreg, comm) result (gwr)

!Arguments ------------------------------------
!scalars
 type(gwr_t),target :: gwr
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
 integer,parameter :: qptopt1 = 1, timrev1 = 1, master = 0
 integer :: my_is, my_it, my_ikf, my_iqf, ii, npwsp, col_bsize, ebands_timrev, my_iki, my_iqi, itau, spin
 integer :: my_nshiftq, ik_ibz, ik_bz, iq_bz, iq_ibz, npw_, ncid ! ig1, ig2,
 integer :: all_nproc, my_rank
 real(dp) :: ecut_eff
 real(dp) :: cpu, wall, gflops
 character(len=5000) :: msg
 type(desc_t),pointer :: desc_k, desc_q
 type(processor_scalapack) :: slk_processor
 type(krank_t) :: qrank, krank_ibz
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: kbz2ibz(:,:), qbz2ibz(:,:), gvec_(:,:)
 !integer,allocatable :: kibz2bz(:), qibz2bz(:), qglob2bz(:,:), kglob2bz(:,:)
 real(dp) :: my_shiftq(3,1), kk_ibz(3), kk_bz(3), qq_bz(3), qq_ibz(3)
 real(dp),allocatable :: wtk(:), kibz(:,:)
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
 gwr%dtfil => dtfil
 gwr%cryst => cryst
 gwr%psps => psps
 gwr%pawtab => pawtab
 gwr%ks_ebands => ks_ebands
 ! TODO Make sure fermie is inside the gap if semiconductor. perhaps this shoud be delegated to gwr_driver
 gwr%kibz => ks_ebands%kptns
 gwr%wtk => ks_ebands%wtk
 gwr%mpi_enreg => mpi_enreg

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
 ! Finally, find the number of IBZ points to be stored in memory by this MPI rank.

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

 ! Find FFT mesh and compute max number of g-vectors

 ecut_eff = dtset%ecut * dtset%dilatmx ** 2
 gwr%ngfft = 0

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
 !ii = gwr%green_mpw; call xmpi_max(ii, gwr%green_mpw, gwr%comm, ierr)

 do iq_bz=1,gwr%nqbz
   qq_bz = gwr%qbz(:, iq_bz)
   call get_kg(qq_bz, istwfk1, dtset%ecuteps, gwr%cryst%gmet, npw_, gvec_)
   ABI_FREE(gvec_)
   call getng(dtset%boxcutmin, dtset%chksymtnons, dtset%ecuteps, cryst%gmet, &
              qq_bz, me_fft0, gwr%mgfft, gwr%nfft, gwr%ngfft, nproc_fft1, cryst%nsym, &
              paral_fft0, cryst%symrel, cryst%tnons, unit=dev_null)
   gwr%chi_mpw = max(gwr%chi_mpw, npw_)
 end do
 !ii = gwr%chi_mpw; call xmpi_max(ii, gwr%chi_mpw, gwr%comm, ierr)

 call print_ngfft(gwr%ngfft, header="FFT mesh for Chi", unit=std_out)

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

 ! =====================================
 ! Allocate arrays storing G_k and Chi_q
 ! =====================================
 ABI_MALLOC(gwr%gt_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))
 ABI_MALLOC(gwr%chit_qibz, (gwr%nqibz, gwr%ntau, gwr%nsppol))

 call init_scalapack(slk_processor, gwr%g_comm%value) !# grid_dims=[1, gwr%g_comm%nproc]

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_it=1,gwr%my_ntau
     itau = gwr%my_itaus(my_it)

     ! Allocate G_k(g, g')
     do my_iki=1,gwr%my_nkibz
       ik_ibz = gwr%my_kibz_inds(my_iki)
       npwsp = gwr%green_desc_kibz(ik_ibz)%npw * gwr%nspinor
       col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
       call init_matrix_scalapack(gwr%gt_kibz(1, ik_ibz, itau, spin), &
          npwsp, npwsp, slk_processor, 1) ! , tbloc=[npwsp, col_bsize])
       call init_matrix_scalapack(gwr%gt_kibz(2, ik_ibz, itau, spin), &
          npwsp, npwsp, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
     end do

    ! Allocate chi_k(g, g')
    do my_iqi=1,gwr%my_nqibz
      iq_ibz = gwr%my_qibz_inds(my_iqi)
      npwsp = gwr%chi_desc_qibz(iq_ibz)%npw * gwr%nspinor
      col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
      call init_matrix_scalapack(gwr%chit_qibz(iq_ibz, itau, spin), &
         npwsp, npwsp, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
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

 call cwtime_report(" gwr_new:", cpu, wall, gflops)

end function gwr_new
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_free
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

 ! Free descriptors
 if (allocated(gwr%green_desc_kibz)) then
   call desc_free_array(gwr%green_desc_kibz)
   ABI_FREE(gwr%green_desc_kibz)
 end if
 if (allocated(gwr%chi_desc_qibz)) then
   call desc_free_array(gwr%chi_desc_qibz)
   ABI_FREE(gwr%chi_desc_qibz)
 end if
 ! Free scalapack matrices
 if (allocated(gwr%gt_kibz)) then
   call slk_free_array(gwr%gt_kibz)
   ABI_FREE(gwr%gt_kibz)
 end if
 if (allocated(gwr%chit_qibz)) then
   call slk_free_array(gwr%chit_qibz)
   ABI_FREE(gwr%chit_qibz)
 end if
 if (allocated(gwr%wct_qibz)) then
   call slk_free_array(gwr%wct_qibz)
   ABI_FREE(gwr%wct_qibz)
 end if
 if (allocated(gwr%sigct_kibz)) then
   call slk_free_array(gwr%sigct_kibz)
   ABI_FREE(gwr%sigct_kibz)
 end if

 ! Free MPI communicators
 call gwr%spin_comm%free()
 call gwr%g_comm%free()
 call gwr%tau_comm%free()
 call gwr%kpt_comm%free()
 call gwr%gtau_comm%free()

contains
  subroutine desc_free_array(desc_array1)
    type(desc_t),intent(inout) :: desc_array1(:)
    integer :: ii
    do ii=1,size(desc_array1, dim=1)
      call desc_array1(ii)%free()
    end do
  end subroutine desc_free_array

end subroutine gwr_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_gtau_from_wfk
!! NAME
!!  gwr_build_gtau_from_wfk
!!
!! FUNCTION
!!  Build the Green's function in imaginary time for k-points in the IBZ from the WFK file
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
 type(ebands_t) :: ks_ebands
 type(hdr_type) :: wfk0_hdr
 type(dataset_type), pointer :: dtset
!arrays
 integer,allocatable :: nband(:,:), wfd_istwfk(:), my_bands(:)
 real(dp),allocatable :: cgwork(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 dtset => gwr%dtset

 ks_ebands = wfk_read_ebands(wfk0_path, gwr%comm, out_hdr=wfk0_hdr)

 call wfk0_hdr%vs_dtset(gwr%dtset)
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

     ABI_MALLOC(cgwork, (2, npw_k * wfd%nspinor))
     !wfd%kdata(ik_ibz)%kg_k

     do ib=1,my_nb
       band = my_bands(ib)
       f_nk = gwr%ks_ebands%occ(band, ik_ibz, spin)
       eig_nk = gwr%ks_ebands%eig(band, ik_ibz, spin) - ef
       call wfd%copy_cg(band, ik_ibz, spin, cgwork)

       !ABI_CHECK(wfd%get_wave_ptr(band, ik_ibz, spin, wave, msg) == 0, msg)
       !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk_ibz, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

       do my_it=1,gwr%my_ntau
         tau = gwr%tau_mesh(gwr%my_itaus(my_it))
         !cgwork(ig1) x cgwork(ig2) * exp(-tau * eig_nk)
         !gwr%gt_kibz(1, ik_ibz, itau, spin)%buffer_cplx =
         !gwr%gt_kibz(2, ik_ibz, itau, spin)%buffer_cplx =
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
 type(matrix_scalapack),intent(out) :: gt_kbz(2)

!Local variables-------------------------------
!scalars
 integer :: ig1, ig2, il_g1, il_g2, ioe, spin, itau
 integer :: ik_ibz, isym_k, trev_k, g0_k(3)
 logical :: isirr_k

! *************************************************************************

 ik_ibz = gwr%my_kbz2ibz(1, my_ikf); isym_k = gwr%my_kbz2ibz(2, my_ikf)
 trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 !tsign = merge(1, -1, trev_k == 0)

 spin = gwr%my_spins(my_is)
 itau = gwr%my_itaus(my_it)

 ! Copy descriptor from IBZ, rotate gvec and recompute gbound.
 desc_kbz = gwr%green_desc_kibz(ik_ibz)%copy()

 if (isirr_k) then
   do ioe=1,2
     gt_kbz(ioe) = gwr%gt_kibz(ioe, ik_ibz, itau, spin)%copy()
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
   gt_kbz(ioe) = gwr%gt_kibz(ioe, ik_ibz, itau, spin)%copy()
   do il_g2=1, gt_kbz(ioe)%sizeb_local(2)
     do il_g1=1, gt_kbz(ioe)%sizeb_local(1)
       ! TODO: factorize call for efficiency reasons.
       call gt_kbz(ioe)%loc2glob(il_g1, il_g2, ig1, ig2)
       gt_kbz(ioe)%buffer_cplx(il_g1, il_g2) = gt_kbz(ioe)%buffer_cplx(il_g1, il_g2) ! * ??
     end do
   end do
 end do

end subroutine gwr_rotate_gt
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_green_gpr
!! NAME
!!  gwr_green_gpr
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

subroutine gwr_green_gpr(gwr, my_it, my_is, desc_kbz, gt_gpr)

!Arguments ------------------------------------
 class(gwr_t),target,intent(inout) :: gwr
 integer,intent(in) :: my_it, my_is
 type(matrix_scalapack),intent(inout) :: gt_gpr(2, gwr%my_nkbz)
 type(desc_t),target,intent(out) :: desc_kbz(gwr%my_nkbz)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1 = 1
 integer :: my_ikf, ig2, ioe
 !integer :: ik_ibz, isym_k, trev_k, g0_k(3)
 real(dp) :: cpu, wall, gflops
 !logical :: isirr_k
 type(matrix_scalapack), pointer :: ggp
 type(matrix_scalapack) :: rgp !, gr
 type(matrix_scalapack),target :: gt_kbz(2)
 type(desc_t),pointer :: desc_k

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

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

     ! Allocate rgp scalapack matrix to store X(r, g')
     call init_matrix_scalapack(rgp, gwr%nfft * gwr%nspinor, ggp%sizeb_global(2), &
                                ggp%processor, desc_k%istwfk) !, tbloc=tbloc)

     ! FFT. Results stored in rgp
     do ig2=1,ggp%sizeb_local(2)
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

   call slk_free_array(gt_kbz)
 end do ! my_ikf

 call cwtime_report(" gwr_green_gpr:", cpu, wall, gflops)

end subroutine gwr_green_gpr
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_cos_transform
!! NAME
!!  gwr_cos_transform
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
 integer :: my_iqi, my_is, ig1, ig2, my_it, iw, ierr, iq_ibz, itau, spin, it0
 real(dp) :: cpu, wall, gflops
 logical :: sum_spins_
 type(desc_t),pointer :: desc_q
 type(matrix_scalapack), pointer :: chit(:)
!arrays
 real(dp),allocatable :: my_weights(:,:)
 complex(dp):: cwork_t(gwr%my_ntau), cwork_wglb(gwr%my_ntau)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 sum_spins_ = .False.; if (present(sum_spins)) sum_spins_ = sum_spins

 ABI_MALLOC(my_weights, (gwr%ntau, gwr%my_ntau))

 ! Extract my weights from the global array according to mode.
 do my_it=1,gwr%my_ntau
   itau = gwr%my_itaus(my_it)
   !my_weights(:, my_it) = cos(omega * tau) * global_weighs(:, itau)
   select case (mode)
   case ("w2t")
   case ("t2w")
   case default
     ABI_ERROR(sjoin("Wrong mode", mode))
   end select
 end do

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     !qq_ibz = gwr%qibz(:, iq_ibz)
     !desc_q => gwr%chi_desc_qibz(iq_ibz)
     !is_gamma = all(abs(qq_ibz) < tol12) ! Handle q --> 0 limit

     select case (what)
     case ("chi")
       desc_q => gwr%chi_desc_qibz(iq_ibz)
       chit => gwr%chit_qibz(iq_ibz,:,spin)
     !case ("wc")
     !  desc_q => gwr%chi_desc_qibz(iq_ibz)
     !  chit => gwr%wct_qibz(iq_ibz, :)
     case default
       ABI_ERROR(sjoin("Invalid value for what:", what))
     end select

     it0 = gwr%my_itaus(1)

     do ig2=1,chit(it0)%sizeb_local(2)
       do ig1=1,chit(it0)%sizeb_local(1)
         ! Extract matrix elements as a function of tau
         ! TODO: Here we can block over ig1 and call zgemm
         ! to reduced the number of MPI communications.
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           cwork_t(my_it) = chit(itau)%buffer_cplx(ig1, ig2)
         end do
         do iw=1,gwr%ntau
           cwork_wglb(iw) = dot_product(my_weights(iw, :), cwork_t)
         end do

         call xmpi_sum(cwork_wglb, gwr%tau_comm%value, ierr)

         ! update values for this (g1, g2), now we have a slice of gg' in frequency space.
         do my_it=1,gwr%my_ntau
           itau = gwr%my_itaus(my_it)
           !chit(itau)%buffer_cplx(ig1, ig2) = cwork_wglb(itau)
         end do

       end do ! ig1
     end do ! ig2

   end do ! my_iqi
 end do ! spin

 if (gwr%nsppol == 2 .and. sum_spins_) then
   !call xmpi_sum(cwork_wglb, gwr%spin_comm%value, ierr)
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

type(desc_t) function desc_copy(desc) result(new_desc)

!Arguments ------------------------------------
 class(desc_t),intent(in) :: desc

! *************************************************************************

 new_desc%istwfk = desc%istwfk
 new_desc%npw = desc%npw

 call alloc_copy(desc%gvec, new_desc%gvec)
 call alloc_copy(desc%gbound, new_desc%gbound)

end function desc_copy
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

!!****f* m_gwr/gwr_build_chi
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
 integer :: my_is, my_it, my_ikf, my_iqf, ig, my_ir, my_nr, npwsp, col_bsize, my_iqi
 integer :: ig1, ig2, sc_nfft, sc_augsize, spin, ik_ibz, ik_bz, iq_bz, iq_ibz, ierr, ioe, itau
 real(dp) :: cpu, wall, gflops
 !character(len=5000) :: msg
 type(desc_t),pointer :: desc_k, desc_q
 type(processor_scalapack) :: slk_processor
 type(matrix_scalapack) :: chi_r_gp
!arrays
 integer :: sc_ngfft(18)
 integer,allocatable :: green_scg(:,:), chi_scg(:,:)
 real(dp) :: kk_bz(3), kpq_bz(3), qq_ibz(3) ! kk_ibz(3), qq_bz(3),
 complex(dp),allocatable :: gf_scbox(:,:), chi_scbox(:)
 type(matrix_scalapack),allocatable :: chiq_gpr(:)
 type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz)
 type(desc_t), target :: desc_kbz(gwr%my_nkbz)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 if (gwr%use_supercell) then
   ! ============================
   ! GWr algorithm with supercell
   ! ============================

   ! Set FFT mesh in the supercell
   ! Be careful when using the FFT plan with ndat as ndat can change inside the loop if we start to block.
   ! Perhaps the safest approach would be to generate the plan on the fly.
   sc_ngfft = gwr%ngfft
   sc_ngfft(1:3) = gwr%ngkpt * gwr%ngfft(1:3)
   sc_ngfft(4:6) = sc_ngfft(1:3)
   sc_nfft = product(sc_ngfft(1:3))
   sc_augsize = product(sc_ngfft(4:6))

   ABI_MALLOC(gf_scbox, (sc_nfft, 2))
   ABI_MALLOC(chi_scbox, (sc_nfft))

   !cal fftbox_plan3_init(sc_plan_gp2rp, ndat, dims, embed, fftalg, fftcache, isign)
   !cal fftbox_plan3_init(sc_plan_rp2gp, ndat, dims, embed, fftalg, fftcache, isign)

   ! The g-vectors in the supercell for Green and chi.
   ABI_MALLOC(green_scg, (3, gwr%green_mpw))
   ABI_MALLOC(chi_scg, (3, gwr%chi_mpw))
   ABI_MALLOC(chiq_gpr, (gwr%my_nqibz))

   do my_iqi=1,gwr%my_nqibz
     iq_ibz = gwr%my_qibz_inds(my_iqi)
     !qq_ibz = gwr%qibz(:, iq_ibz)
     !desc_q => gwr%chi_desc_qibz(iq_ibz)
     npwsp = gwr%chi_desc_qibz(iq_ibz)%npw * gwr%nspinor
     col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
     call init_matrix_scalapack(chiq_gpr(my_iqi), &
          npwsp, gwr%nfft * gwr%nspinor, slk_processor, 1) !, tbloc=[npwsp, col_bsize])
   end do

   ! Loop over my spins and my taus.
   do my_is=1,gwr%my_nspins
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       !tau = gwr%tau_mesh(itau)

       ! G_k(g,g') --> G_k(g',r) for each k in the BZ treated by me.
       call gwr%green_gpr(my_it, my_is, desc_kbz, gt_gpr)

       ! Loop over the second index r (MPI-distributed in g_comm).
       my_nr = gt_gpr(1, 1)%sizeb_local(2)
       do my_ir=1,my_nr

         ! Insert G_k(g',r) in G'-space in the supercell FFT box.
         ! Note that we need to take the union of (k, g') for k in the BZ.
         ! TODO: Here I need a specialized version of sphere that does not initialize
         ! out array to zero as I need to accumulate.
         gf_scbox = zero
         do my_ikf=1,gwr%my_nkbz
           ik_bz = gwr%my_kbz_inds(my_ikf)
           ik_ibz = gwr%my_kbz2ibz(1, my_ikf) !; isym_k = gwr%my_kbz2ibz(2, my_ikf)
           !trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
           !isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))

           desc_k => desc_kbz(my_ikf)
           do ig=1,desc_k%npw
             green_scg(:,ig) = nint(gwr%kbz(:, ik_bz) * gwr%ngkpt) + gwr%ngkpt * desc_k%gvec(:,ig)
           end do

           do ioe=1,2
             gt_gpr(ioe, my_ikf)%buffer_cplx(:, my_ir) = zero
             !call insert_in_box(gt_gpr(ioe, my_ikf), green_scg, gf_scbox(:,ioe))
             !call sphere(cg,ndat,npw,gf_scbox,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
           end do
         end do ! my_ikf

         ! Should block using nproc in kpt_comm, scatter data and perform multiple FFTs in parallel.
         call xmpi_sum(gf_scbox, gwr%kpt_comm%value, ierr)

         ! Compute chi(R',r) for this r and go back to G'=k+g' space immediately.
         !call fftbox_execute_ip(sc_plan_gp2rp, gf_scbox)

         chi_scbox = gf_scbox(:, 1) * gf_scbox(:, 2)
         !call fftbox_execute_ip(sc_plan_rp2gp, chi_scbox)

         ! Extract chi_q(g, r) on the g-sphere from the FFT box in the supercell.
         ! not necessarly equal to the one used for the Green's function.
         ! NB: For the time being I'm assuming k-mesh == q-mesh.
         ! Only q-points in the IBZ are considered.
         do my_iqi=1,gwr%my_nqibz
           iq_ibz = gwr%my_qibz_inds(my_iqi)
           qq_ibz = gwr%qibz(:, iq_ibz)
           desc_q => gwr%chi_desc_qibz(iq_ibz)
           do ig=1,desc_q%npw
             chi_scg(:,ig) = nint(gwr%qibz(:, iq_ibz) * gwr%ngqpt) + gwr%ngqpt * desc_q%gvec(:,ig)
           end do

           !call sphere(cg,ndat,npw,cfft_o,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,me_g0,shiftg,symm,xnorm)
           chiq_gpr(my_iqi)%buffer_cplx(:, my_ir) = zero
         end do

       end do ! my_ir

       do my_iqi=1,gwr%my_nqibz
         iq_ibz = gwr%my_qibz_inds(my_iqi)
         !qq_ibz = gwr%qibz(:, iq_ibz)
         desc_q => gwr%chi_desc_qibz(iq_ibz)
         !iq_bz = gwr%my_qbz_inds(my_iqf)
         !chi_r_gp = chiq_gpr(my_iqi)%alloc_trans()
         !call chiq_gpr(my_iqi)%ptrans("C", chi_r_gp)
         call chiq_gpr(my_iqi)%free()
         ! FFT along the first dimension: chi_q(r, g') --> chi_q(g, g')
         do ig2=1,chi_r_gp%sizeb_local(2)
           chi_r_gp%buffer_cplx(:, ig2) = zero
           !call fft_ur(npw_k, nfft, nspinor, ndat, mgfft, ngfft, istwf_k, kg_k, gbound_k, ur, ug)
           gwr%chit_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) = zero ! ug
         end do
         call chi_r_gp%free()
       end do

     end do ! my_it
   end do ! spin

   ABI_FREE(gf_scbox)
   ABI_FREE(chi_scbox)
   ABI_FREE(green_scg)
   ABI_FREE(chi_scg)
   call slk_free_array(chiq_gpr)
   ABI_FREE(chiq_gpr)

  else
    ! ===================================================================
    ! Mixed-space algorithm in the unit cell with convolutions in k-space
    ! ===================================================================
    ! Each MPI proc has all the G_k in the IBZ thus we only need to rotate to get the BZ on the fly.
    do my_is=1,gwr%my_nspins
      spin = gwr%my_spins(my_is)
      do my_it=1,gwr%my_ntau
        itau = gwr%my_itaus(my_it)
        !tau = gwr%tau_mesh(itau)

        do my_iqi=1,gwr%my_nqibz
          iq_ibz = gwr%my_qibz_inds(my_iqi)
          qq_ibz = gwr%qibz(:, iq_ibz)
          desc_q => gwr%chi_desc_qibz(iq_ibz)
          !iq_bz = gwr%my_qbz_inds(my_iqf)
          !qq_bz = gwr%qbz(:, iq_bz)

          do my_ikf=1,gwr%my_nkbz
            !ik_ibz = gwr%my_kbz2ibz(1, my_ikf); isym_k = gwr%my_kbz2ibz(2, my_ikf)
            !trev_k = gwr%my_kbz2ibz(6, my_ikf); g0_k = gwr%my_kbz2ibz(3:5, my_ikf)
            !isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
            ik_bz = gwr%my_kbz_inds(my_ikf)

            ! Use symmetries to get G_kq and G_k from the corresponding images in the IBZ.
            kk_bz = gwr%kbz(:, ik_bz)
            kpq_bz = kk_bz + qq_ibz
            !ikpq_ibz = ??

            ! G_k(g, g') --> G_k(r, r')
            !call gwr%green_rrp(my_ikf, my_it, my_is, gtk_rrp)
            !call gwr%green_rrp(my_ikqf, my_it, my_is, gtkq_rrp)

            !desc_k = gwr%green_desc_kibz(ik_ibz)%rotate(symtab_k)
            !desc_kpq = gwr%green_desc_kibz(ikpq_ibz)%rotate(symtab_kpq)

            !my_nr = go_k_gpr%sizeb_local(2)
            !do my_ir=1,my_nr
            !   !chi_fft = go_rpr * ge_rpr
            !   ! from r' to g'
            !   !chi_gp_r(:, my_ir) = FFT_rg(chi_fft)
            !end do

            !call desc_k%free()
            !call desc_kpq%free()

            !call chi_gp_r%ptrans("C", chi_r_gp)
            !call chi_gp_r%free()
            !do ig2=1, chi_rp_g%sizeb_local(2)
              !chi_r_gp%buffer_cplx(:, ig2)
              !!gwr%chit_qibz(iq_ibz, itau, spin)%buffer_cplx(:, ig2) = FFT_rg
            !end do
            !call chi_r_gp%free()
          end do ! my_ikf
        end do ! iq_ibz

     end do ! my_it
   end do ! spin
 end if

 ! Free scalapack matrices
 call slk_free_array(gt_gpr)

 ! Transform irreducible chi from imaginary tau to imaginary omega.
 ! Also sum over spins to get total chi if collinear spins.
 call gwr%cos_transform("chi", "t2w", sum_spins=.True.)

 call cwtime_report(" gwr_buid_chi:", cpu, wall, gflops)

end subroutine gwr_build_chi
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_wc
!! NAME
!!  gwr_build_wc
!!
!! FUNCTION
!!  Compute W(tau) from chi(i omega)
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
 integer :: my_iqi, my_it, ig1, ig2, ig1_glob, ig2_glob, my_is, iq_ibz, spin, itau ! ierr,
 real(dp) :: cpu, wall, gflops
 logical :: is_gamma
 type(desc_t), pointer :: desc_q
 type(matrix_scalapack) :: eps
 !arrays
 real(dp) :: qq_ibz(3)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ABI_MALLOC(gwr%wct_qibz, (gwr%nqibz, gwr%ntau))

 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   desc_q => gwr%chi_desc_qibz(iq_ibz)
   is_gamma = all(abs(qq_ibz) < tol12) ! Handle q --> 0 limit

   !npwsp = gwr%chi_desc_qibz(iq_ibz)%npw * gwr%nspinor
   !col_bsize = npwsp / gwr%g_comm%nproc; if (mod(npwsp, gwr%g_comm%nproc) /= 0) col_bsize = col_bsize + 1
   !call init_matrix_scalapack(gwr%chi(my_iqf, my_it, my_is), &
   !   npwsp, npwsp, slk_processor, 1) ! , tbloc=[npwsp, col_bsize])

   do my_is=1,gwr%my_nspins  ! FIXME
     spin = gwr%my_spins(my_is)
     do my_it=1,gwr%my_ntau
       itau = gwr%my_itaus(my_it)
       !gwr%wct_qibz(iq_ibz, itau) =
       eps = gwr%chit_qibz(iq_ibz, itau, spin)%copy()

       ! Build epsilon(q,iw) = delta_{g,g'} - v_q(g,g') chi_q(g,g,iw).
       do ig2=1,eps%sizeb_local(2)
         do ig1=1,eps%sizeb_local(1)
           call eps%loc2glob(ig1, ig2, ig1_glob, ig2_glob)
           !vc_qg1 = ??
           !vc_qg2 = ??
           !ceps = - eps%buffer_cplx(ig1, ig2) ! * vc_qg1 * vc_qg2
           !if (ig1_glob == ig2_glob) eps = one - ceps
         end do ! ig1
       end do ! ig2

       ! IMPORTANT NOTE: PZGETRF requires square block decomposition i.e., MB_A = NB_A.
       ! This means here I need to redistribute the data before calling zinvert.
       !call slk_zinvert(eps)
       !call eps%free()

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

       !gwr%wct_qibz(iq_ibz, itau) = ??
       call eps%free()
     end do  ! my_it
   end do ! my_is
 end do ! my_iqi

 ! Cosine transform from iw to tau to get Wc(tau)
 call gwr%cos_transform("wc", "w2t")

 call cwtime_report(" gwr_build_wc:", cpu, wall, gflops)

end subroutine gwr_build_wc
!!***

!----------------------------------------------------------------------

!!****f* m_gwr/gwr_build_sigmac
!! NAME
!!  gwr_build_sigmac
!!
!! FUNCTION
!!  Compute Sigma(tau)
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
 integer :: my_is, my_it, my_iw, my_iki, spin, ik_ibz, ioe
 integer :: my_iqi, ig1, ig2, ig1_glob, ig2_glob, iq_ibz ! ierr,
 real(dp) :: cpu, wall, gflops
 logical :: is_gamma
 type(desc_t), pointer :: desc_q, desc_k
 !type(matrix_scalapack) :: eps
!arrays
 real(dp) :: kk_ibz(3), qq_ibz(3)
 type(matrix_scalapack) :: gt_gpr(2, gwr%my_nkbz)
 type(desc_t), target :: desc_kbz(gwr%my_nkbz)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 !if (gwr%use_supercell) then

 ABI_MALLOC(gwr%sigct_kibz, (2, gwr%nkibz, gwr%ntau, gwr%nsppol))

 do my_is=1,gwr%my_nspins
   spin = gwr%my_spins(my_is)
   do my_iki=1,gwr%my_nkibz
     ik_ibz = gwr%my_kibz_inds(my_iki)
     kk_ibz = gwr%kibz(:, ik_ibz)
     desc_k => gwr%green_desc_kibz(ik_ibz)

     do my_it=1,gwr%my_ntau
       !tau = gwr%tau_mesh(gwr%my_itaus(my_it))
       ! G_k(g, g') --> G_k(g', r) for each k in the BZ treated by me.
       call gwr%green_gpr(my_it, my_is, desc_kbz, gt_gpr)
     end do
   end do
 end do

 do my_iqi=1,gwr%my_nqibz
   iq_ibz = gwr%my_qibz_inds(my_iqi)
   qq_ibz = gwr%qibz(:, iq_ibz)
   desc_q => gwr%chi_desc_qibz(iq_ibz)
   is_gamma = all(abs(qq_ibz) < tol12) ! Handle q --> 0 limit

   !do my_is=1,gwr%my_nspins  ! FIXME
   !  do my_it=1,gwr%my_ntau
   !  end do  ! my_it
   !end do ! my_is
 end do ! my_iqi

 ! Cosine transform from iw to tau to get W(tau)
 !call gwr%cos_transform("wc", "w2t")

 call cwtime_report(" gwr_build_sigmac:", cpu, wall, gflops)

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

 call gwr%build_chi()
 call gwr%build_wc()
 call gwr%build_sigmac()

end subroutine gwr_run_g0w0
!!***

end module m_gwr
!!***

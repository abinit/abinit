!!****m* ABINIT/m_gemm_nonlop_projectors
!! NAME
!! m_gemm_nonlop
!!
!! FUNCTION
!!  This module provides functions to compute the nonlocal operator by means of the BLAS GEMM
!!  routine. By treating ndat simultaneous wavefunctions, it is able to exploit BLAS3 routines,
!!  which leads to excellent CPU efficiency and OpenMP scalability.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2026 ABINIT group (AL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

! TODO list :
! Don't allocate the full nkpt structures, only those that are treated by this proc: use same init as in m_bandfft_kpt
! support more options (forces & stresses mostly)
! Support RF/other computations (only GS right now)
! handle the case where nloalg(2) < 0, ie no precomputation of ph3d
! more systematic checking of the workflow (right now, only works if init/make/gemm/destroy, no multiple makes, etc)
! Avoid allocating the complex matrix when istwfk > 1
! Merge with chebfi's invovl


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gemm_nonlop_projectors

 use defs_basis
 use m_errors
 use m_abicore
 use m_xomp
 use m_xmpi
 use m_fstrings,    only : itoa, ftoa, sjoin
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_kg, only : mkkpg
 use m_hamiltonian, only : KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

#if defined(HAVE_GPU_CUDA)
 use m_alloc_hamilt_gpu, only : gemm_nonlop_gpu_data
#endif

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_int32_t, c_int64_t, c_float, c_double, c_size_t, c_loc, c_ptr
#endif

 implicit none

 private

 public :: init_gemm_nonlop
 public :: destroy_gemm_nonlop
 public :: reset_gemm_nonlop
 public :: set_gemm_nonlop_ikpt
 public :: prep_projectors
 public :: prep_dprojectors
 public :: refresh_projectors

!!***

!----------------------------------------------------------------------

!!****t* m_gemm_nonlop_projectors/gemm_nonlop_type
!! NAME
!! gemm_nonlop_type
!!
!! FUNCTION
!! Contains information needed to apply the nonlocal operator
!!
!! SOURCE
 type,public :: gemm_nonlop_type

   integer :: npw
   integer :: nprojs
   integer :: ngrads
   integer :: ngrads2

   integer :: nprojs_blk
   integer :: nprojs_last_blk

   real(dp), allocatable :: projs(:, :, :)
   ! (2, npw, nprojs)

   real(dp), allocatable :: projs_r(:, :, :)
   ! (1, npw, nprojs)

   real(dp), allocatable :: projs_i(:, :, :)
   ! (1, npw, nprojs)

   real(dp), allocatable :: dprojs(:, :, :)
   ! (2, npw, nprojs*ngrads)
   real(dp), allocatable :: dprojs_r(:, :, :)
   ! (1, npw, nprojs*ngrads)
   real(dp), allocatable :: dprojs_i(:, :, :)
   ! (1, npw, nprojs*ngrads)

   real(dp), allocatable :: d2projs(:, :, :)
   ! (2, npw, nprojs*ngrads)

   integer :: idir
   integer :: ikpt
   integer :: choice

 end type gemm_nonlop_type
!!***

 type(gemm_nonlop_type), save, public, target :: gemm_nonlop_kpt(2)
 !(nkpt)

 integer, save, public :: gemm_nonlop_ikpt_this_proc_being_treated
 !! This is oh so very crude, but I can't find any other way to do it without passing ikpt deep down to nonlop

 logical, save, public :: gemm_nonlop_use_gemm = .false.
 ! Public variable indicating whether we should call gemm_nonlop or fall back to the usual nonlop. Set to false
 ! in order not to interfere with non-GS calls to nonlop.

 logical, save, public :: gemm_nonlop_is_distributed = .false.
 ! Public variable indicating whether we should gemm_nonlop operated in a distributed manner. Set to false by default
 ! but might be enabled by memory constraints or forced by user through parameters.

 integer, save :: gemm_nonlop_nblocks = 1
 ! How many blocks of MPI tasks should the projs arrays be ditributed.

 integer, save, public :: gemm_nonlop_block_comm = xmpi_comm_null
 ! MPI communicator for MPI tasks processing the same gemm_nonlop block for projs array distribution

 integer, save, public :: gemm_nonlop_block_size = 0
 ! Public variable indicating size of a block (ie: number of MPI tasks in gemm_nonlop_block_comm)
 ! Default size 0 indicates no distribution at all.

 integer, save, public :: gemm_nonlop_choice = -1

 integer, save, public :: gemm_nonlop_gpu_option = ABI_GPU_DISABLED

 real(dp),save, allocatable, target :: atom_projs(:,:,:)
 real(dp),save, allocatable, target :: atom_dprojs(:,:,:,:)
 real(dp),save, allocatable, target :: atom_d2projs(:,:,:,:)
 integer,save, allocatable, target :: scal(:)
 integer,save, allocatable, target :: lmn_parity(:)
 integer, save :: mod__lmnmax, mod__npw, mod__ndprojs, mod__nd2projs
 ! Work arrays for prep_*projectors functions. Sized after mod__lmnmax, mod__npw and mod__ndprojs.

#if defined(HAVE_FC_ISO_C_BINDING) && defined(HAVE_GPU_CUDA)

 type, bind(c), public :: gemm_nonlop_gpu_type

   integer(kind=c_int32_t) :: npw
   integer(kind=c_int32_t) :: nprojs

   ! array of double on GPU, dimensions are (2, npw, nprojs)
   type(c_ptr) :: projs

   ! array of double on GPU, dimensions are (1, npw, nprojs)
   type(c_ptr) :: projs_r

   ! array of double on GPU, dimensions are (1, npw, nprojs)
   type(c_ptr) :: projs_i

 end type gemm_nonlop_gpu_type

 !! array of size nkpt of sobjects of type gemm_nonlop_gpu_type, array size is nkpt
 type(gemm_nonlop_gpu_type), save, public, target :: gemm_nonlop_kpt_gpu(2)
 !(nkpt)

#endif

!!***

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/init_gemm_nonlop
!! NAME
!! init_gemm_nonlop
!!
!! FUNCTION
!! Initalization of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine init_gemm_nonlop(gpu_option)

  integer,intent(in) :: gpu_option

! *************************************************************************

  gemm_nonlop_kpt(:)%npw = -1
  gemm_nonlop_kpt(:)%nprojs = -1
  gemm_nonlop_kpt(:)%ngrads = -1
  gemm_nonlop_kpt(:)%ngrads2 = -1
  gemm_nonlop_kpt(:)%choice = -1
  gemm_nonlop_kpt(:)%idir = -1
  gemm_nonlop_kpt(:)%ikpt = -1

  if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    gemm_nonlop_kpt_gpu(:)%npw = -1
    gemm_nonlop_kpt_gpu(:)%nprojs = -1
    gemm_nonlop_gpu_data % allocated = .false.
#endif
  end if

  gemm_nonlop_block_comm=xmpi_comm_null
  gemm_nonlop_block_size=0
  gemm_nonlop_nblocks=1
  gemm_nonlop_gpu_option=gpu_option

 end subroutine init_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/destroy_gemm_nonlop
!! NAME
!! destroy_gemm_nonlop
!!
!! FUNCTION
!! Destruction of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine destroy_gemm_nonlop(gpu_option)

  integer,intent(in) :: gpu_option

! *************************************************************************

  call free_gemm_nonlop_ikpt(1,gpu_option)
  call free_gemm_nonlop_ikpt(2,gpu_option)
  call destroy_work_arrays(gpu_option)
  if(gemm_nonlop_block_comm/=xmpi_comm_null) call xmpi_comm_free(gemm_nonlop_block_comm)

 end subroutine destroy_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/alloc_work_arrays
!! NAME
!! alloc_work_arrays
!!
!! FUNCTION
!! Allocation of work arrays
!!
!! INPUTS
!! gpu_option = which GPU code path is used
!!
!! SOURCE
 subroutine alloc_work_arrays(lmnmax,npw,ndprojs,nd2projs,gpu_option)

  integer,intent(in) :: lmnmax,npw,ndprojs,nd2projs,gpu_option

! *************************************************************************

  !FIXME Would be nice to not allocate/reallocate at each call, but it seem troublesome in practice
  !if(mod__lmnmax>=lmnmax .and. mod__npw>=npw .and. mod__ndprojs>=ndprojs .and. mod__nd2projs>=nd2projs) then
  !  return ! Nothing to do
  !end if

  call destroy_work_arrays(gpu_option)

  ABI_MALLOC(atom_projs, (2, npw, lmnmax))
  ABI_MALLOC(scal, (lmnmax))
  ABI_MALLOC(lmn_parity, (lmnmax))
#ifdef HAVE_OPENMP_OFFLOAD
  !$OMP TARGET ENTER DATA MAP(alloc:atom_projs,scal,lmn_parity) IF(gpu_option==ABI_GPU_OPENMP)
#endif

  if(ndprojs>0) then
    ABI_MALLOC(atom_dprojs, (2, npw, ndprojs, lmnmax))
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET ENTER DATA MAP(alloc:atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
#endif
  end if

  if(nd2projs>0) then
    ABI_MALLOC(atom_d2projs, (2, npw, nd2projs, lmnmax))
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET ENTER DATA MAP(alloc:atom_d2projs) IF(gpu_option==ABI_GPU_OPENMP)
#endif
  end if

  mod__lmnmax=lmnmax
  mod__npw=npw
  mod__ndprojs=ndprojs
  mod__nd2projs=nd2projs

 end subroutine alloc_work_arrays
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/destroy_work_arrays
!! NAME
!! destroy_work_arrays
!!
!! FUNCTION
!! Destruction of work arrays
!!
!! INPUTS
!! gpu_option = which GPU code path is used
!!
!! SOURCE
 subroutine destroy_work_arrays(gpu_option)

  integer,intent(in) :: gpu_option

! *************************************************************************
  ABI_UNUSED((/gpu_option/))

  if(allocated(atom_projs)) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET EXIT DATA MAP(delete:atom_projs) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    ABI_FREE(atom_projs)
  end if
  if(allocated(atom_dprojs)) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET EXIT DATA MAP(delete:atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    ABI_FREE(atom_dprojs)
  end if
  if(allocated(atom_d2projs)) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET EXIT DATA MAP(delete:atom_d2projs) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    ABI_FREE(atom_d2projs)
  end if
  if(allocated(scal)) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET EXIT DATA MAP(delete:scal) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    ABI_FREE(scal)
  end if
  if(allocated(lmn_parity)) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET EXIT DATA MAP(delete:lmn_parity) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    ABI_FREE(lmn_parity)
  end if
  mod__lmnmax=0
  mod__npw=0
  mod__ndprojs=0
  mod__nd2projs=0

 end subroutine destroy_work_arrays
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/free_gemm_nonlop_ikpt
!! NAME
!! free_destroy_gemm_nonlop_ikpt
!!
!! FUNCTION
!! Release memory for one kpt value of the gemm_nonlop_kpt array
!!
!! INPUTS
!! ikpt= index of gemm_nonlop_kptto be released
!!
!! SOURCE
 subroutine free_gemm_nonlop_ikpt(ik, gpu_option)

  integer,intent(in) :: ik, gpu_option

! *************************************************************************

 if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
   if(gemm_nonlop_kpt_gpu(ik)%nprojs /= -1) then
     ! deallocate arrays projs, projs_r and projs_i
     if (allocated(gemm_nonlop_kpt(ik)%projs)) then
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ik)%projs)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%projs_r)) then
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ik)%projs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%projs_i)) then
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ik)%projs_i)
     end if
     gemm_nonlop_kpt_gpu(ik)%nprojs = -1
     gemm_nonlop_kpt_gpu(ik)%npw = -1
   end if
#endif
 end if

 if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   call free_ompgpu_current_ikpt(ik)
#endif
 end if

 if(gemm_nonlop_kpt(ik)%nprojs /= -1) then
   if (allocated(gemm_nonlop_kpt(ik)%projs)) then
     ABI_FREE(gemm_nonlop_kpt(ik)%projs)
   end if
   if (allocated(gemm_nonlop_kpt(ik)%projs_r)) then
     ABI_FREE(gemm_nonlop_kpt(ik)%projs_r)
   end if
   if (allocated(gemm_nonlop_kpt(ik)%projs_i)) then
   ABI_FREE(gemm_nonlop_kpt(ik)%projs_i)
   end if
   gemm_nonlop_kpt(ik)%nprojs = -1
   if(gemm_nonlop_kpt(ik)%ngrads /= -1) then
     if (allocated(gemm_nonlop_kpt(ik)%dprojs)) then
       ABI_FREE(gemm_nonlop_kpt(ik)%dprojs)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%dprojs_r)) then
       ABI_FREE(gemm_nonlop_kpt(ik)%dprojs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%dprojs_i)) then
       ABI_FREE(gemm_nonlop_kpt(ik)%dprojs_i)
     end if
     gemm_nonlop_kpt(ik)%ngrads = -1
   end if
   if(gemm_nonlop_kpt(ik)%ngrads2 /= -1) then
     if (allocated(gemm_nonlop_kpt(ik)%d2projs)) then
       ABI_FREE(gemm_nonlop_kpt(ik)%d2projs)
     end if
     gemm_nonlop_kpt(ik)%ngrads2 = -1
   end if
 end if
 gemm_nonlop_kpt(ik)%choice = -1
 gemm_nonlop_kpt(ik)%idir = -1
 gemm_nonlop_kpt(ik)%ikpt = -1

 if(gemm_nonlop_is_distributed) then
   gemm_nonlop_kpt(ik)%nprojs_blk = -1
   gemm_nonlop_kpt(ik)%nprojs_last_blk = -1
 end if

 end subroutine free_gemm_nonlop_ikpt
!!***

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt(ik)

  integer,intent(in) :: ik

#ifdef HAVE_OPENMP_OFFLOAD

  call free_ompgpu_current_ikpt_projs(ik)
  call free_ompgpu_current_ikpt_dprojs(ik)

#else
  ABI_UNUSED((/ik/))
#endif
 end subroutine free_ompgpu_current_ikpt

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt_projs(ik)

  integer,intent(in) :: ik
#ifdef HAVE_OPENMP_OFFLOAD
  !NOTE: Those pointers exists to be served to OpenMP TARGET directives to hide
  !      the datastructure gemm_nonlop_kpt which is not supported in GCC, LLVM and Cray
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_projs_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_projs_r_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_projs_i_ompptr(:,:,:)


  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%projs))) then
    gemm_nonlop_kpt_projs_ompptr => gemm_nonlop_kpt(ik)%projs
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_projs_ompptr)
  end if
  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%projs_r))) then
    gemm_nonlop_kpt_projs_r_ompptr => gemm_nonlop_kpt(ik)%projs_r
    gemm_nonlop_kpt_projs_i_ompptr => gemm_nonlop_kpt(ik)%projs_i
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_projs_r_ompptr)
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_projs_i_ompptr)
  end if

#else
  ABI_UNUSED((/ik/))
#endif
 end subroutine free_ompgpu_current_ikpt_projs

!----------------------------------------------------------------------


 subroutine free_ompgpu_current_ikpt_dprojs(ik)

  integer,intent(in) :: ik
#ifdef HAVE_OPENMP_OFFLOAD
  !NOTE: Those pointers exists to be served to OpenMP TARGET directives to hide
  !      the datastructure gemm_nonlop_kpt which is not supported in GCC, LLVM and Cray
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_dprojs_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_dprojs_r_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_dprojs_i_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_d2projs_ompptr(:,:,:)

  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%dprojs))) then
    gemm_nonlop_kpt_dprojs_ompptr => gemm_nonlop_kpt(ik)%dprojs
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_dprojs_ompptr)
  end if
  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%dprojs_i))) then
    gemm_nonlop_kpt_dprojs_r_ompptr => gemm_nonlop_kpt(ik)%dprojs_r
    gemm_nonlop_kpt_dprojs_i_ompptr => gemm_nonlop_kpt(ik)%dprojs_i
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_dprojs_r_ompptr)
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_dprojs_i_ompptr)
  end if

  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%d2projs))) then
    gemm_nonlop_kpt_d2projs_ompptr => gemm_nonlop_kpt(ik)%d2projs
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt_d2projs_ompptr)
  end if

#else
  ABI_UNUSED((/ik/))
#endif
 end subroutine free_ompgpu_current_ikpt_dprojs

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/set_gemm_nonlop_ikpt
!! NAME
!! set_gemm_nonlop_ikpt
!!
!! FUNCTION
!! Set the K-point upon which projectors will be computed and
!! pre-allocate projectors buffers.
!!
!! INPUTS
!! ikpt= K-point id
!! npw= number of plane-wave
!! istwf_k=option parameter that describes the storage of wfs
!! indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!! ntypat=number of atoms types
!! nattyp(ntypat)=number of atoms of each type
!! gpu_option=which variant of GEMM nonlop is used
!!
!! SOURCE
 subroutine set_gemm_nonlop_ikpt(ikpt,npw,istwf_k,indlmn,ntypat,nattyp,gpu_option)

  integer,intent(in) :: ikpt,istwf_k,npw,ntypat,gpu_option
  integer,intent(in) :: indlmn(:,:,:), nattyp(ntypat)

  integer :: nprojs, itypat

! *************************************************************************

  gemm_nonlop_ikpt_this_proc_being_treated=ikpt

  nprojs=0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  ! Call a "dummy" refresh of projectors buffers
  ! This is mostly a work-around in GPU workloads to ensure there
  ! is a buffer allocated in GPU memory.
  call refresh_projectors(npw,istwf_k,nprojs,0,0,.false.,gpu_option)

 end subroutine set_gemm_nonlop_ikpt
!!***

!!****f* m_gemm_nonlop_projectors/reset_gemm_nonlop
!! NAME
!! reset_gemm_nonlop
!!
!! FUNCTION
!! Reset projectors to trigger their recomputation
!!
!! INPUTS
!!
!! SOURCE
 subroutine reset_gemm_nonlop()

! *************************************************************************

  gemm_nonlop_kpt(:)%ikpt   = -1
  gemm_nonlop_kpt(:)%choice = -1
  gemm_nonlop_kpt(:)%idir   = -1

 end subroutine reset_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/refresh_projectors
!! NAME
!! prep_projectors
!!
!! FUNCTION
!! Check allocation of projectors arrays for GEMM nonlop, and resize them if need be
!!
!! INPUTS
!!
!! SOURCE
 subroutine refresh_projectors(npw,istwf_k,nprojs,ndgxdt,nd2gxdt,&
 &                             is_kprime,gpu_option)
  integer,intent(in) :: npw,istwf_k,nprojs,ndgxdt,nd2gxdt,gpu_option
  logical,intent(in) :: is_kprime
  integer :: ik,rank,nprojs_blk,nprojs_my_blk,nprojs_last_blk,ierr,nprocs
  logical :: is_last_rank
#ifdef HAVE_OPENMP_OFFLOAD
  !NOTE: Those pointers exists to be served to OpenMP TARGET directives to hide
  !      the datastructure gemm_nonlop_kpt which is not supported in GCC, LLVM and Cray
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_projs_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_projs_r_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_projs_i_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_dprojs_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_dprojs_r_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_dprojs_i_ompptr(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: gemm_nonlop_kpt_d2projs_ompptr(:,:,:)
#endif

  ik=1; if(is_kprime) ik=2
  if(gemm_nonlop_kpt(ik)%ikpt/=gemm_nonlop_ikpt_this_proc_being_treated &
  &   .or. npw/=gemm_nonlop_kpt(ik)%npw .or. nprojs/=gemm_nonlop_kpt(ik)%nprojs) then
    call free_gemm_nonlop_ikpt(ik, gpu_option)
  end if

  if(gemm_nonlop_is_distributed) then
    nprocs = xmpi_comm_size(xmpi_world)
    ! If split size has changed, reset array and init MPI communicator
    if(gemm_nonlop_block_comm==xmpi_comm_null .or. gemm_nonlop_nblocks /= nprocs/gemm_nonlop_block_size) then
      call free_gemm_nonlop_ikpt(ik, gpu_option)
      if(gemm_nonlop_block_comm/=xmpi_comm_null) call xmpi_comm_free(gemm_nonlop_block_comm)
      rank = xmpi_comm_rank(xmpi_world);
      gemm_nonlop_nblocks=nprocs/gemm_nonlop_block_size
      write(std_out,'(A,I3,A,I3,A)')  "Splitting GEMM nonlop projectors on ",&
      &    gemm_nonlop_nblocks, " blocks of ", gemm_nonlop_block_size, " MPI tasks..."
      call xmpi_comm_split(xmpi_world, rank/gemm_nonlop_block_size, rank, gemm_nonlop_block_comm, ierr)
      if(ierr/=0) ABI_BUG("MPI_comm_split failed!")
    end if
  end if

  nprojs_last_blk = nprojs
  nprojs_my_blk = nprojs
  nprojs_blk = nprojs
  rank = 0; is_last_rank = .true.

  if(gemm_nonlop_block_size > 1) then
    nprojs_blk = nprojs / gemm_nonlop_block_size
    nprojs_last_blk = nprojs_blk + modulo(nprojs,nprojs_blk)

    if(gemm_nonlop_is_distributed) then
      rank = xmpi_comm_rank(gemm_nonlop_block_comm);
      is_last_rank = (rank==gemm_nonlop_block_size-1)
      if(is_last_rank) then
        nprojs_my_blk = nprojs_last_blk
      else
        nprojs_my_blk = nprojs_blk
      end if
    end if
  end if


  ! Allocation of buffers for 1st and 2nd order derivatives of projectors
  ! NOTE: those are allocated, if needed, before regular projectors buffers
  !       for optimization purposes, regarding GPU memory pool.
  if(nprojs>0) then
  if(ndgxdt>0) then
    if(npw/=gemm_nonlop_kpt(ik)%npw .or. nprojs/=gemm_nonlop_kpt(ik)%nprojs &
    &    .or. ndgxdt /= gemm_nonlop_kpt(ik)%ngrads .or. nd2gxdt /=  gemm_nonlop_kpt(ik)%ngrads2) then
      if(gpu_option == ABI_GPU_OPENMP) call free_ompgpu_current_ikpt_dprojs(ik)
      ABI_SFREE(gemm_nonlop_kpt(ik)%dprojs)
      ABI_SFREE(gemm_nonlop_kpt(ik)%dprojs_r)
      ABI_SFREE(gemm_nonlop_kpt(ik)%dprojs_i)
      ABI_SFREE(gemm_nonlop_kpt(ik)%d2projs)
      gemm_nonlop_kpt(ik)%ngrads = -1
      gemm_nonlop_kpt(ik)%ngrads2 = -1

      if(istwf_k <= 1) then
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs, (2, npw, nprojs_last_blk*ndgxdt))
#ifdef HAVE_OPENMP_OFFLOAD
        gemm_nonlop_kpt_dprojs_ompptr => gemm_nonlop_kpt(ik)%dprojs
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_dprojs_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
        if(nd2gxdt>0) then
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs, (2, npw, nprojs_last_blk*nd2gxdt))
#ifdef HAVE_OPENMP_OFFLOAD
          gemm_nonlop_kpt_d2projs_ompptr => gemm_nonlop_kpt(ik)%d2projs
          !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_d2projs_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
        end if
      else
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs_r, (1, npw, nprojs_last_blk*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs_i, (1, npw, nprojs_last_blk*ndgxdt))
#ifdef HAVE_OPENMP_OFFLOAD
        gemm_nonlop_kpt_dprojs_r_ompptr => gemm_nonlop_kpt(ik)%dprojs_r
        gemm_nonlop_kpt_dprojs_i_ompptr => gemm_nonlop_kpt(ik)%dprojs_i
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_dprojs_r_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_dprojs_i_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
      end if
    end if
  end if


  ! Allocation of projectors buffers
  if(npw/=gemm_nonlop_kpt(ik)%npw .or. nprojs/=gemm_nonlop_kpt(ik)%nprojs) then
    if(istwf_k <= 1) then
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs, (2, npw, nprojs_last_blk))
#ifdef HAVE_OPENMP_OFFLOAD
      gemm_nonlop_kpt_projs_ompptr => gemm_nonlop_kpt(ik)%projs
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_projs_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    else
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs_r, (1, npw, nprojs_last_blk))
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs_i, (1, npw, nprojs_last_blk))
#ifdef HAVE_OPENMP_OFFLOAD
      gemm_nonlop_kpt_projs_r_ompptr => gemm_nonlop_kpt(ik)%projs_r
      gemm_nonlop_kpt_projs_i_ompptr => gemm_nonlop_kpt(ik)%projs_i
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_projs_r_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt_projs_i_ompptr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
    end if
  end if
  end if


  if (nprojs>0) gemm_nonlop_kpt(ik)%nprojs = nprojs
  if (nprojs>0) gemm_nonlop_kpt(ik)%npw = npw
  if (ndgxdt>0) gemm_nonlop_kpt(ik)%ngrads = ndgxdt
  if (nd2gxdt>0) gemm_nonlop_kpt(ik)%ngrads2 = nd2gxdt
  if(gemm_nonlop_block_size > 1) then
    if(nprojs_blk>0) gemm_nonlop_kpt(ik)%nprojs_blk = nprojs_blk
    if(nprojs_last_blk>0) gemm_nonlop_kpt(ik)%nprojs_last_blk = nprojs_last_blk
  end if

  !!!!! CUDA stuff
  if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    if(gemm_nonlop_kpt_gpu(ik)%nprojs==-1) then
      gemm_nonlop_kpt_gpu(ik)%npw    = npw
      gemm_nonlop_kpt_gpu(ik)%nprojs = nprojs

#ifdef DEBUG_VERBOSE_GPU
      if(xmpi_comm_rank(xmpi_world) == 0) then
        call check_gpu_mem("refresh_projectors begin")
        call wrtout(std_out,sjoin(" npw    .......", itoa(npw)),    'COLL')
        call wrtout(std_out,sjoin(" nprojs .......", itoa(nprojs)), 'COLL')
      end if
#endif

      if(istwf_k <= 1) then
        call alloc_on_gpu(gemm_nonlop_kpt_gpu(ik)%projs, INT(2,c_size_t)*npw*nprojs*dp)
        ! TODO : gradients
      else
        call alloc_on_gpu(gemm_nonlop_kpt_gpu(ik)%projs_r, INT(1, c_size_t)*npw*nprojs*dp)
        call alloc_on_gpu(gemm_nonlop_kpt_gpu(ik)%projs_i, INT(1, c_size_t)*npw*nprojs*dp)
        ! TODO : gradients
      end if

#ifdef DEBUG_VERBOSE_GPU
      if(xmpi_comm_rank(xmpi_world) == 0) then
        call check_gpu_mem("refresh_projectors end  ")
      end if
#endif
    end if

#endif
  end if
 end subroutine refresh_projectors
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/prep_projectors
!! NAME
!! prep_projectors
!!
!! FUNCTION
!! Prepare projectors array for GEMM nonlop (choice=={0,1})
!!
!! INPUTS
!!
!! SOURCE
 subroutine prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
 &                          ucvol,ffnl,ph3d,dimffnl,matblk,&
 &                          nprojs,is_kprime,gpu_option,&
 &                          iblock)

  integer, intent(in) :: npw,lmnmax,ntypat,dimffnl,matblk
  integer, intent(in) :: istwf_k
  integer, intent(in) :: nprojs,gpu_option,iblock
  logical, intent(in) :: is_kprime
  real(dp), intent(in) :: ucvol
  ! arrays
  integer, intent(in) :: indlmn(6,lmnmax,ntypat),nattyp(ntypat)
  real(dp),intent(in),target :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(in),target :: ph3d(2,npw,matblk)

  logical :: map_ffnl,map_ph3d,is_last_rank
  integer :: il, ipw, ik, ilmn_p, nlmn_p
  integer :: itypat, ilmn, nlmn, ia, iaph3d, shift, nprojs_my_blk
  integer :: lmn_beg,ibeg,iend,shift_do,nlmn_o
  real(dp):: wt,tmp
  real(dp), ABI_CONTIGUOUS pointer ::   projs(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_r(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_i(:,:,:)

  ik=1; if(is_kprime) ik=2
  if(istwf_k <= 1) then
    projs => gemm_nonlop_kpt(ik)%projs
  else
    projs_r => gemm_nonlop_kpt(ik)%projs_r
    projs_i => gemm_nonlop_kpt(ik)%projs_i
  end if

  if(gpu_option==ABI_GPU_OPENMP) then

    if(istwf_k <= 1) then
      call gpu_set_to_zero(projs,int(2,c_size_t)*npw*nprojs)
    else
      call gpu_set_to_zero(projs_r,int(npw,c_size_t)*nprojs)
      call gpu_set_to_zero(projs_i,int(npw,c_size_t)*nprojs)
    end if

  else

    if(istwf_k <= 1) then
      projs(:,:,:) = zero
    else
      projs_r(:,:,:) = zero
      projs_i(:,:,:) = zero
    end if

 end if

  iaph3d = 1
  wt=four_pi/sqrt(ucvol)

  ! Allocate atom_projs and other work arrays if need be
  call alloc_work_arrays(lmnmax,npw,-1,-1,gpu_option)

  map_ph3d=.false.; map_ffnl=.false.
#ifdef HAVE_OPENMP_OFFLOAD
  if(.not. xomp_target_is_present(c_loc(ffnl))) map_ffnl = .true.
  if(.not. xomp_target_is_present(c_loc(ph3d))) map_ph3d = .true.
  !$OMP TARGET ENTER DATA MAP(to:ffnl) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl)
  !$OMP TARGET ENTER DATA MAP(to:ph3d) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d)
#endif

  shift = 0
  lmn_beg = 1

  if(gemm_nonlop_block_size > 1) then
    is_last_rank = (iblock==gemm_nonlop_block_size)
    nprojs_my_blk=gemm_nonlop_kpt(ik)%nprojs_blk
    if(is_last_rank) nprojs_my_blk=gemm_nonlop_kpt(ik)%nprojs_last_blk
    ibeg = (iblock-1)*gemm_nonlop_kpt(ik)%nprojs_blk+1
    iend = ibeg+nprojs_my_blk
    shift_do = 0
  end if

  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)
    nlmn_o = nlmn

    do ia = 1, nattyp(itypat)

      ! In distributed mode, loops are skipped until reach the section
      ! of "ilmn" to be stored by local rank
      if(gemm_nonlop_block_size > 1) then
        if(shift_do+nlmn < ibeg) then
          shift_do = shift_do + nlmn
          iaph3d = iaph3d + 1
          cycle
        end if

        lmn_beg = max(1,ibeg-shift_do)
        if(shift_do+nlmn > iend) nlmn = iend - shift_do - 1
      end if

      !! build atom_projs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l

      ! start from 4pi/sqrt(ucvol)*ffnl
      ! atom_projs(1, :, lmn_beg:nlmn) = four_pi/sqrt(ham%ucvol) * ham%ffnl(:, 1, lmn_beg:nlmn)
      ! TODO vectorize (DCOPY with stride)
      if(gpu_option==ABI_GPU_OPENMP) then
        call gpu_set_to_zero(atom_projs,int(2,c_size_t)*npw*lmnmax)
      else
        atom_projs(:,:,:) = zero
      end if
#ifdef HAVE_OPENMP_OFFLOAD
      !$OMP TARGET TEAMS DISTRIBUTE  PARALLEL DO COLLAPSE(2) &
      !$OMP& PRIVATE(ilmn,ipw) MAP(to:atom_projs,ffnl) &
      !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
      do ilmn=1,nlmn_o
        do ipw=1, npw
          atom_projs(1,ipw, ilmn) = wt * ffnl(ipw, 1, ilmn, itypat)
        end do
      end do

      nlmn_p=0
      ! multiply by (-i)^l
      do ilmn=1,nlmn_o
        il=mod(indlmn(1,ilmn, itypat),4);
        if(.not. (mod(il,2)==0)) then
          nlmn_p=nlmn_p+1
          lmn_parity(nlmn_p)=ilmn
        end if
        scal(ilmn)=1; if(il>1) scal(ilmn)=-1
      end do
      ! multiply by -1
      if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET UPDATE TO(scal,lmn_parity)
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !$OMP& PRIVATE(ipw,ilmn) MAP(to:atom_projs,scal) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=1,nlmn_o
          do ipw=1,npw
            atom_projs(1,ipw,ilmn) = atom_projs(1,ipw,ilmn) * scal(ilmn)
            atom_projs(2,ipw,ilmn) = atom_projs(2,ipw,ilmn) * scal(ilmn)
          end do
        end do
#endif
      else
        do ilmn=1,nlmn_o
          atom_projs(:,:,ilmn) = atom_projs(:,:,ilmn) * scal(ilmn)
        end do
      end if
      ! multiply by -i
#ifdef HAVE_OPENMP_OFFLOAD
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
      !$OMP& MAP(to:atom_projs,lmn_parity) PRIVATE(ilmn,ilmn_p,ipw,tmp) &
      !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
      do ilmn_p=1,nlmn_p
        do ipw=1,npw
          ilmn=lmn_parity(ilmn_p)
          tmp = atom_projs(2,ipw,ilmn)
          atom_projs(2,ipw,ilmn) = -atom_projs(1,ipw,ilmn)
          atom_projs(1,ipw,ilmn) =  tmp
        end do
      end do

      ! multiply by conj(ph3d)
#ifdef HAVE_OPENMP_OFFLOAD
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
      !$OMP& PRIVATE(ilmn,ipw,tmp) MAP(to:atom_projs,ph3d) &
      !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
      do ilmn=1,nlmn_o
        do ipw=1,npw
          tmp = atom_projs(1, ipw, ilmn)
          atom_projs(1, ipw, ilmn) = atom_projs(1, ipw, ilmn) * ph3d(1, ipw, iaph3d) &
          &                        + atom_projs(2, ipw, ilmn) * ph3d(2, ipw, iaph3d)
          atom_projs(2, ipw, ilmn) = atom_projs(2, ipw, ilmn) * ph3d(1, ipw, iaph3d) &
          &                        - tmp                      * ph3d(2, ipw, iaph3d)
        end do
      end do

      !! atom_projs is built, copy to projs

      if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
        if(istwf_k <= 1) then
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
          !$OMP& PRIVATE(ipw,ilmn) MAP(to:projs,atom_projs)
          do ilmn=1,nlmn-(lmn_beg-1)
            do ipw=1,npw
              projs(:, ipw, shift+ilmn) = atom_projs(:, ipw, ilmn+(lmn_beg-1))
            end do
          end do
        else ! istwf_k>1
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
          !$OMP& PRIVATE(ipw,ilmn) MAP(to:projs_r,projs_i,atom_projs)
          do ilmn=1,nlmn-(lmn_beg-1)
            do ipw=1,npw
              projs_r(1, ipw, shift+ilmn) = atom_projs(1, ipw, ilmn+(lmn_beg-1))
              projs_i(1, ipw, shift+ilmn) = atom_projs(2, ipw, ilmn+(lmn_beg-1))
            end do
          end do
        end if
#endif
      else
        if(istwf_k <= 1) then
          projs(1:2, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(:, :, lmn_beg:nlmn)
        else ! istwf_k>1
          projs_r(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(1, :, lmn_beg:nlmn)
          projs_i(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(2, :, lmn_beg:nlmn)
        end if
      end if

      iaph3d = iaph3d + 1
      if(gemm_nonlop_block_size > 1) then
        shift = shift + nlmn - (lmn_beg-1)
        shift_do = shift_do + nlmn
        if(shift_do >= iend - 1) exit
      else
        shift = shift + nlmn
      end if
    end do
    if(gemm_nonlop_block_size > 1 .and. shift_do >= iend - 1) exit
  end do


  !!!!! CUDA stuff
  if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    ! upload data to gpu memory
    if(istwf_k <= 1) then
      call copy_on_gpu(gemm_nonlop_kpt(ik)%projs, gemm_nonlop_kpt_gpu(ik)%projs, INT(2, c_size_t)*npw*nprojs*dp)
      ! TODO : gradients
    else
      call copy_on_gpu(gemm_nonlop_kpt(ik)%projs_r, gemm_nonlop_kpt_gpu(ik)%projs_r, &
        &                    INT(1, c_size_t)*npw*nprojs*dp)
      call copy_on_gpu(gemm_nonlop_kpt(ik)%projs_i, gemm_nonlop_kpt_gpu(ik)%projs_i, &
        &                    INT(1, c_size_t)*npw*nprojs*dp)
    end if
#endif
  end if

#ifdef HAVE_OPENMP_OFFLOAD
  !$OMP TARGET EXIT DATA MAP(delete:ffnl) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl)
  !$OMP TARGET EXIT DATA MAP(delete:ph3d) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d)
#endif

 end subroutine prep_projectors
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_projectors/prep_dprojectors
!! NAME
!! prep_projectors
!!
!! FUNCTION
!! Prepare projectors' derivatives for given choice (choice={2,3,4,5,51,54,55})
!!
!! INPUTS
!!
!! SOURCE
 subroutine prep_dprojectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
 &                          ucvol,ffnl,ph3d,kpg,nkpg,dimffnl,matblk,&
 &                          nprojs,ngrads,ngrads2,choice,signs,idir_pert,&
 &                          is_kprime,gpu_option,iblock)

  integer, intent(in) :: npw,lmnmax,ntypat,nkpg,dimffnl,matblk
  integer, intent(in) :: istwf_k,iblock
  integer, intent(in) :: nprojs,ngrads,ngrads2,choice,signs,idir_pert,gpu_option
  logical, intent(in) :: is_kprime
  real(dp), intent(in) :: ucvol
  ! arrays
  integer, intent(in) :: indlmn(6,lmnmax,ntypat),nattyp(ntypat)
  real(dp), intent(in), target :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp), intent(in), target :: ph3d(2,npw,matblk)
  real(dp), intent(in), target :: kpg(npw,nkpg)

  logical :: map_ffnl,map_ph3d,is_last_rank
  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
  integer :: ndprojs, nd2projs, ilmn_p, nlmn_p
  integer :: il, ipw, ik, idir, idir1, idir2, jdir1, jdir2, kdir1, kdir2, ldir1, ldir2, ldir3, ldir4, ffnl_dir
  integer :: itypat, ilmn, nlmn, ia, iaph3d, igrad, shift, shift_grad, shift_grad2
  integer :: lmn_beg,ibeg,iend,shift_do,nlmn_o,lmn_grad_beg,nprojs_my_blk
  real(dp), parameter :: two_pi2=two_pi*two_pi
  real(dp):: wt,tmp
  real(dp), ABI_CONTIGUOUS pointer :: projs  (:,:,:)  ,dprojs(:,:,:)  ,d2projs(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_r(:,:,:),dprojs_r(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_i(:,:,:),dprojs_i(:,:,:)


  ik=1; if(is_kprime) ik=2
  if(istwf_k <= 1) then
    projs => gemm_nonlop_kpt(ik)%projs
    dprojs => gemm_nonlop_kpt(ik)%dprojs
    if(ngrads2 > 0) then
      d2projs => gemm_nonlop_kpt(ik)%d2projs
    end if
  else
    projs_r => gemm_nonlop_kpt(ik)%projs_r
    projs_i => gemm_nonlop_kpt(ik)%projs_i
    dprojs_r => gemm_nonlop_kpt(ik)%dprojs_r
    dprojs_i => gemm_nonlop_kpt(ik)%dprojs_i
  end if

  if(gpu_option==ABI_GPU_OPENMP) then

    if(istwf_k <= 1) then
      call gpu_set_to_zero(dprojs,int(2,c_size_t)*npw*nprojs*ngrads)
      if(ngrads2 > 0) then
        call gpu_set_to_zero(d2projs,int(2,c_size_t)*npw*nprojs*ngrads2)
      end if
    else
      call gpu_set_to_zero(dprojs_r,int(npw,c_size_t)*nprojs*ngrads)
      call gpu_set_to_zero(dprojs_i,int(npw,c_size_t)*nprojs*ngrads)
    end if

  else

    if(istwf_k <= 1) then
      dprojs(:,:,:) = zero
      if(ngrads2 > 0) then
        d2projs(:,:,:) = zero
      end if
    else
      dprojs_r(:,:,:) = zero
      dprojs_i(:,:,:) = zero
    end if

  end if

  iaph3d = 1
  wt=four_pi/sqrt(ucvol)
  ffnl_dir=1; if(dimffnl>2) ffnl_dir=idir_pert

  ndprojs = 0
  if (signs==1 .and. (choice==3 .or. choice==23 .or. choice==54 .or. choice==55 .or. choice==6)) then
    ndprojs = 3
  else if(signs==2 .and. (choice==5 .or. choice==51 .or. choice==3)) then
    ndprojs = 1
  end if
  nd2projs = 0
  if(signs==1 .and. choice==54) then
    nd2projs = 3
  else if(signs==1 .and. choice==55) then
    nd2projs = 6
  else if(signs==1 .and. choice==6) then
    nd2projs = 10
  end if

  ! Allocate atom_projs and other work arrays if need be
  call alloc_work_arrays(lmnmax,npw,ndprojs,nd2projs,gpu_option)

  map_ph3d=.false.; map_ffnl=.false.
  if(.not. xomp_target_is_present(c_loc(ffnl))) map_ffnl = .true.
  if(.not. xomp_target_is_present(c_loc(ph3d))) map_ph3d = .true.
#ifdef HAVE_OPENMP_OFFLOAD
  !$OMP TARGET ENTER DATA MAP(to:ffnl) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl)
  !$OMP TARGET ENTER DATA MAP(to:ph3d) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d)
  !!$OMP TARGET ENTER DATA MAP(to:kpg)  IF(gpu_option==ABI_GPU_OPENMP)
#endif

  shift = 0 ; shift_grad = 0; shift_grad2 = 0
  lmn_beg = 1

  if(gemm_nonlop_block_size > 1) then
    is_last_rank = (iblock==gemm_nonlop_block_size)
    nprojs_my_blk=gemm_nonlop_kpt(ik)%nprojs_blk
    if(is_last_rank) nprojs_my_blk=gemm_nonlop_kpt(ik)%nprojs_last_blk
    ibeg = (iblock-1)*gemm_nonlop_kpt(ik)%nprojs_blk+1
    iend = ibeg+nprojs_my_blk
    shift_do = 0
    lmn_grad_beg = -1

  end if

  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)
    nlmn_o = nlmn

    do ia = 1, nattyp(itypat)

      ! In distributed mode, loops are skipped until reach the section
      ! of "ilmn" to be stored by local rank
      if(gemm_nonlop_block_size > 1) then
        if(shift_do+nlmn < ibeg) then
          shift_do = shift_do + nlmn
          iaph3d = iaph3d + 1
          cycle
        end if

        lmn_beg = max(1,ibeg-shift_do)
        if(lmn_grad_beg==-1) lmn_grad_beg = (lmn_beg-1)*ngrads
        if(shift_do+nlmn > iend) nlmn = iend - shift_do - 1
      end if

      !! build atom_dprojs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l

      ! start from 4pi/sqrt(ucvol)*ffnl
      if(gpu_option==ABI_GPU_OPENMP) then
        if (ndprojs>0) then
          call gpu_set_to_zero(atom_dprojs,int(2,c_size_t)*npw*ndprojs*lmnmax)
          if(ngrads2>0) then
            call gpu_set_to_zero(atom_d2projs,int(2,c_size_t)*npw*nd2projs*lmnmax)
          end if
        end if
      else
        if (ndprojs>0) atom_dprojs(:,:,:,:) = zero
        if (nd2projs>0) atom_d2projs(:,:,:,:) = zero
      end if
      if (signs==1 .and. (choice==3 .or. choice==23 .or. choice==54 .or. choice==55 .or. choice==6)) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_dprojs,ffnl) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ipw=1, npw
          atom_dprojs(1,ipw, 1:ndprojs, 1:nlmn_o) = wt * ffnl(ipw, 2:ndprojs+1, 1:nlmn_o, itypat)
        end do
      end if
      if (signs==2 .and. (choice==3 .or. choice==5 .or. choice==51)) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_dprojs,ffnl) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ipw=1, npw
          atom_dprojs(1,ipw, 1, 1:nlmn_o) = wt * ffnl(ipw, 1+ffnl_dir, 1:nlmn_o, itypat)
        end do
      end if
      if(signs==1 .and. choice==54) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ipw=1, npw
          atom_d2projs(1,ipw, 1:nd2projs, 1:nlmn_o) = wt * ffnl(ipw, 2:nd2projs+1, 1:nlmn_o, itypat)
        end do
      else if(signs==1 .and. choice==55) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ipw=1, npw
          atom_d2projs(1,ipw, 1:nd2projs, 1:nlmn_o) = wt * ffnl(ipw, 5:nd2projs+4, 1:nlmn_o, itypat)
        end do
      else if(signs==1 .and. choice==6) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ipw=1, npw
          atom_d2projs(1,ipw, 1:10, 1:nlmn_o) = wt * ffnl(ipw, 1:10, 1:nlmn_o, itypat)
        end do
      end if

      nlmn_p=0
      ! multiply by (-i)^l
      if (ndprojs>0) then
        do ilmn=1,nlmn_o
          il=mod(indlmn(1,ilmn, itypat),4);
          if(.not. (mod(il,2)==0)) then
            nlmn_p=nlmn_p+1
            lmn_parity(nlmn_p)=ilmn
          end if
          scal(ilmn)=1; if(il>1) scal(ilmn)=-1
        end do
        ! multiply by -1
        if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET UPDATE TO(scal,lmn_parity)
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& PRIVATE(ilmn,idir,ipw) MAP(to:atom_dprojs,scal) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do ilmn=1,nlmn_o
            do idir=1,ndprojs
              do ipw=1,npw
                atom_dprojs(1,ipw,idir,ilmn) = atom_dprojs(1,ipw,idir,ilmn) * scal(ilmn)
                atom_dprojs(2,ipw,idir,ilmn) = atom_dprojs(2,ipw,idir,ilmn) * scal(ilmn)
              end do
            end do
          end do
          if(nd2projs>0) then
            !$OMP TARGET TEAMS DISTRIBUTE &
            !$OMP& PRIVATE(ilmn) MAP(to:atom_d2projs,scal) &
            !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
            do ilmn=1,nlmn_o
              !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(idir,ipw)
              do idir=1,nd2projs
                do ipw=1,npw
                  atom_d2projs(1,ipw,idir,ilmn) = atom_d2projs(1,ipw,idir,ilmn) * scal(ilmn)
                  atom_d2projs(2,ipw,idir,ilmn) = atom_d2projs(2,ipw,idir,ilmn) * scal(ilmn)
                end do
              end do
            end do
          end if
#endif
        else
          do ilmn=1,nlmn_o
            atom_dprojs(:,:,:,ilmn) = atom_dprojs(:,:,:,ilmn) * scal(ilmn)
          end do
          if(nd2projs>0) then
            do ilmn=1,nlmn_o
              atom_d2projs(:,:,:,ilmn) = atom_d2projs(:,:,:,ilmn) * scal(ilmn)
            end do
          end if
        end if
        ! multiply by -i
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP MAP(to:atom_dprojs,lmn_parity) PRIVATE(ilmn,ilmn_p,idir,ipw,tmp) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn_p=1,nlmn_p
          do idir=1,ndprojs
            do ipw=1,npw
              ilmn=lmn_parity(ilmn_p)
              tmp = atom_dprojs(2,ipw,idir,ilmn)
              atom_dprojs(2,ipw,idir,ilmn) = -atom_dprojs(1,ipw,idir,ilmn)
              atom_dprojs(1,ipw,idir,ilmn) =  tmp
            end do
          end do
        end do
        if(ngrads2>0) then
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& MAP(to:atom_d2projs) PRIVATE(idir,ipw,ilmn) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
          do ilmn=1,nlmn_o
            do idir=1,nd2projs
              do ipw=1,npw
                atom_d2projs(1,ipw,idir,ilmn) = -atom_d2projs(1,ipw,idir,ilmn)
                atom_d2projs(2,ipw,idir,ilmn) = -atom_d2projs(2,ipw,idir,ilmn)
              end do
            end do
          end do
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& MAP(to:atom_d2projs,lmn_parity) PRIVATE(ilmn,ilmn_p,idir,ipw,tmp) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
          do ilmn_p=1,nlmn_p
            do idir=1,nd2projs
              do ipw=1,npw
                ilmn=lmn_parity(ilmn_p)
                tmp = atom_d2projs(2,ipw,idir,ilmn)
                atom_d2projs(2,ipw,idir,ilmn) = -atom_d2projs(1,ipw,idir,ilmn)
                atom_d2projs(1,ipw,idir,ilmn) =  tmp
              end do
            end do
          end do
        end if
      end if

      ! multiply by conj(ph3d)
      if (ndprojs>0) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP& PRIVATE(ilmn,tmp,ipw,idir) MAP(to:atom_dprojs,ph3d) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=1,nlmn_o
          do idir=1,ndprojs
            do ipw=1,npw
              tmp = atom_dprojs(1, ipw, idir,ilmn)
              atom_dprojs(1, ipw, idir,ilmn) = atom_dprojs(1, ipw, idir,ilmn) * ph3d(1, ipw, iaph3d) &
              &                              + atom_dprojs(2, ipw, idir,ilmn) * ph3d(2, ipw, iaph3d)
              atom_dprojs(2, ipw, idir,ilmn) = atom_dprojs(2, ipw, idir,ilmn) * ph3d(1, ipw, iaph3d) &
              &                              - tmp                 * ph3d(2, ipw, iaph3d)
            end do
          end do
        end do
      end if
      if (nd2projs>0) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE &
        !$OMP& PRIVATE(ilmn,idir,ipw,tmp) MAP(to:atom_d2projs,ph3d) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=1,nlmn_o
          !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ipw,idir,tmp)
          do idir=1,nd2projs
            do ipw=1,npw
              tmp = atom_d2projs(1, ipw, idir,ilmn)
              atom_d2projs(1, ipw, idir,ilmn) = atom_d2projs(1, ipw, idir,ilmn) * ph3d(1, ipw, iaph3d) &
              &                               + atom_d2projs(2, ipw, idir,ilmn) * ph3d(2, ipw, iaph3d)
              atom_d2projs(2, ipw, idir,ilmn) = atom_d2projs(2, ipw, idir,ilmn) * ph3d(1, ipw, iaph3d) &
              &                               - tmp                 * ph3d(2, ipw, iaph3d)
            end do
          end do
        end do
      end if

      !! Handling dprojs

      if(signs==1 .and. (choice==3 .or. choice==23 .or. choice==55 .or. choice==6)) then
        if(istwf_k <= 1) then
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& PRIVATE(ilmn,ipw,idir,idir1,idir2) MAP(to:atom_dprojs,dprojs,kpg,ipw,idir,idir1,idir2) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
          do ilmn=lmn_beg,nlmn
            do idir=1,6
              do ipw=1,npw
                idir1=alpha(idir);idir2=beta(idir)
                dprojs(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+idir) = &
                &     -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))
                dprojs(2, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+idir) = &
                &     -half*(atom_dprojs(2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(2, ipw, idir2, ilmn)*kpg(ipw,idir1))
              end do
            end do
          end do
        else ! istwf_k>1
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& PRIVATE(ilmn,ipw,idir,idir1,idir2) MAP(to:atom_dprojs,dprojs_r,dprojs_i,kpg) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
          do ilmn=lmn_beg,nlmn
            do idir=1,6
              do ipw=1,npw
                idir1=alpha(idir);idir2=beta(idir)
                dprojs_r(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+idir) = &
                &     -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))

                dprojs_i(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+idir) = &
                &     -half*(atom_dprojs(2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(2, ipw, idir2, ilmn)*kpg(ipw,idir1))
              end do
            end do
          end do
        end if
      end if


      if(signs==1 .and. (choice==2 .or. choice==23 .or. choice==4 .or. choice==54 .or. choice==6)) then
        igrad=0; if(choice==23 .or. choice==6) igrad=6
        if(istwf_k <= 1) then
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& PRIVATE(ilmn,ipw,idir) MAP(to:projs,dprojs,kpg) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
          do ilmn=lmn_beg,nlmn
            do idir=1,3
              do ipw=1,npw
                dprojs(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+igrad+idir) = &
                &     +projs(2, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir)*two_pi
                dprojs(2, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+igrad+idir) = &
                &     -projs(1, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir)*two_pi
              end do
            end do
          end do
        else ! istwf_k>1
#ifdef HAVE_OPENMP_OFFLOAD
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& PRIVATE(ilmn,ipw,idir) MAP(to:projs_r,projs_i,dprojs_r,dprojs_i,kpg) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
          do ilmn=lmn_beg,nlmn
            do idir=1,3
              do ipw=1,npw
                dprojs_r(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+igrad+idir) = &
                &     +projs_i(1, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir)*two_pi
                dprojs_i(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+igrad+idir) = &
                &     -projs_r(1, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir)*two_pi
              end do
            end do
          end do
        end if
      end if

      if(signs==1 .and. (choice==5 .or. choice==51 .or. choice==54 .or. choice==55)) then
        igrad=0; if(choice==54) igrad=3; if(choice==55) igrad=6
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP& PRIVATE(ilmn,ipw,idir) MAP(to:atom_dprojs,dprojs) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do idir=1,3
            do ipw=1,npw
              dprojs(1, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+igrad+idir) = &
              &     +atom_dprojs(1, ipw, idir, ilmn)
              dprojs(2, ipw, shift_grad+(ilmn-lmn_beg)*ngrads+igrad+idir) = &
              &     +atom_dprojs(2, ipw, idir, ilmn)
            end do
          end do
        end do
      end if

      if(signs==2 .and. (choice==2)) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !$OMP& PRIVATE(ilmn,ipw) MAP(to:projs,dprojs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+ilmn) = &
            &      projs(2, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir_pert)*two_pi
            dprojs(2, ipw, shift_grad+ilmn) = &
            &     -projs(1, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir_pert)*two_pi
          end do
        end do
      end if


      if(signs==2 .and. (choice==3)) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !$OMP& PRIVATE(ilmn,ipw) MAP(to:atom_dprojs,dprojs) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+ilmn) = &
            &     -atom_dprojs(1, ipw, 1, ilmn)
            dprojs(2, ipw, shift_grad+ilmn) = &
            &     -atom_dprojs(2, ipw, 1, ilmn)
          end do
        end do
      end if


      if(signs==2 .and. (choice==5 .or. choice==51)) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !$OMP& PRIVATE(ilmn,ipw) MAP(to:atom_dprojs,dprojs) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+ilmn) = &
            &     +atom_dprojs(1, ipw, 1, ilmn)
            dprojs(2, ipw, shift_grad+ilmn) = &
            &     +atom_dprojs(2, ipw, 1, ilmn)
          end do
        end do
      end if



      ! Handling d2projs

      if(signs==1 .and. choice==4) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
        !$OMP& PRIVATE(idir,ilmn,ipw) MAP(to:projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do idir=1,6
            do ipw=1,npw
              d2projs(1, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+idir) = &
              &     -projs(1, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir+3)*two_pi2
              d2projs(2, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+idir) = &
              &     -projs(2, ipw, shift+ilmn-lmn_beg+1)*kpg(ipw,idir+3)*two_pi2
            end do
          end do
        end do
      end if

      if(signs==1 .and. choice==54) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP& PRIVATE(idir1,idir2,ilmn,ipw) MAP(to:atom_d2projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do idir1=1,3
            do idir2=1,3
              do ipw=1,npw
                d2projs(1, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+(idir1-1)*3+idir2) = &
                &     -atom_d2projs(2, ipw, idir2, ilmn)*kpg(ipw,idir1)*two_pi
                d2projs(2, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+(idir1-1)*3+idir2) = &
                &     +atom_d2projs(1, ipw, idir2, ilmn)*kpg(ipw,idir1)*two_pi
              end do
            end do
          end do
        end do
      end if

      if(signs==1 .and. choice==55) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP& PRIVATE(idir1,idir2,ilmn,ipw) MAP(to:atom_d2projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do idir1=1,6
            do idir2=1,3
              do ipw=1,npw
                d2projs(1, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+(idir1-1)*3+idir2) = &
                &     +atom_d2projs(1, ipw, idir1, ilmn)*kpg(ipw,idir2)
                d2projs(2, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+(idir1-1)*3+idir2) = &
                &     +atom_d2projs(2, ipw, idir1, ilmn)*kpg(ipw,idir2)
              end do
            end do
          end do
        end do
      end if

      if(signs==1 .and. choice==6) then
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP& PRIVATE(idir1,idir2,jdir1,jdir2,kdir1,kdir2,ldir1,ldir2,ldir3,ldir4,ilmn,ipw) &
        !$OMP& MAP(to:atom_d2projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do idir1=1,6
            do idir2=1,6
              do ipw=1,npw
                jdir1=alpha(idir1);jdir2=beta(idir1)
                kdir1=alpha(idir2);kdir2=beta(idir2)
                ldir1=gamma(jdir1,kdir1)
                ldir2=gamma(jdir2,kdir1)
                ldir3=gamma(jdir1,kdir2)
                ldir4=gamma(jdir2,kdir2)
                d2projs(1, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+(idir1-1)*6+idir2) = &
                &     -atom_d2projs(1, ipw, 4+ldir1, ilmn)*kpg(ipw,jdir2)*kpg(ipw,kdir2)*quarter &
                &     -atom_d2projs(1, ipw, 4+ldir2, ilmn)*kpg(ipw,jdir1)*kpg(ipw,kdir2)*quarter &
                &     -atom_d2projs(1, ipw, 4+ldir3, ilmn)*kpg(ipw,jdir2)*kpg(ipw,kdir1)*quarter &
                &     -atom_d2projs(1, ipw, 4+ldir4, ilmn)*kpg(ipw,jdir1)*kpg(ipw,kdir1)*quarter
                d2projs(2, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+(idir1-1)*6+idir2) = &
                &     -atom_d2projs(2, ipw, 4+ldir1, ilmn)*kpg(ipw,jdir2)*kpg(ipw,kdir2)*quarter &
                &     -atom_d2projs(2, ipw, 4+ldir2, ilmn)*kpg(ipw,jdir1)*kpg(ipw,kdir2)*quarter &
                &     -atom_d2projs(2, ipw, 4+ldir3, ilmn)*kpg(ipw,jdir2)*kpg(ipw,kdir1)*quarter &
                &     -atom_d2projs(2, ipw, 4+ldir4, ilmn)*kpg(ipw,jdir1)*kpg(ipw,kdir1)*quarter
              end do
            end do
          end do
        end do
        igrad=36
#ifdef HAVE_OPENMP_OFFLOAD
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
        !$OMP& PRIVATE(idir1,idir2,jdir1,jdir2,ilmn,ipw) &
        !$OMP& MAP(to:atom_d2projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
        do ilmn=lmn_beg,nlmn
          do idir1=1,6 !mub
            do idir2=1,3 !mua
              do ipw=1,npw
                jdir1=alpha(idir1);jdir2=beta(idir1)
                d2projs(1, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+igrad+(idir1-1)*3+idir2) = &
                &     +kpg(ipw,idir2)*(atom_d2projs(2,ipw,1+jdir1,ilmn)*kpg(ipw,jdir2)*pi &
                &     +                atom_d2projs(2,ipw,1+jdir2,ilmn)*kpg(ipw,jdir1)*pi)

                d2projs(2, ipw, shift_grad2+(ilmn-lmn_beg)*ngrads2+igrad+(idir1-1)*3+idir2) = &
                &     -kpg(ipw,idir2)*(atom_d2projs(1,ipw,1+jdir1,ilmn)*kpg(ipw,jdir2)*pi &
                &     +                atom_d2projs(1,ipw,1+jdir2,ilmn)*kpg(ipw,jdir1)*pi)
              end do
            end do
          end do
        end do
      end if

      iaph3d = iaph3d + 1
      shift_grad2 = shift_grad2 + ngrads2*(nlmn-lmn_beg+1)
      shift_grad  = shift_grad  + ngrads*(nlmn-lmn_beg+1)

      if(gemm_nonlop_block_size > 1) then
        shift = shift + nlmn - (lmn_beg-1)
        shift_do = shift_do + nlmn
        if(shift_do >= iend - 1) exit
      else
        shift = shift + nlmn
      end if

    end do
    if(gemm_nonlop_block_size > 1 .and. shift_do >= iend - 1) exit
  end do

#ifdef HAVE_OPENMP_OFFLOAD
  !!$OMP TARGET EXIT DATA MAP(delete:kpg) IF(gpu_option==ABI_GPU_OPENMP)
  !$OMP TARGET EXIT DATA MAP(delete:ffnl) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl)
  !$OMP TARGET EXIT DATA MAP(delete:ph3d) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d)
#endif

 end subroutine prep_dprojectors
!!***

end module m_gemm_nonlop_projectors
!!***

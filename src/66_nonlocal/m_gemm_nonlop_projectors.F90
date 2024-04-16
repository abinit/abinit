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
!! Copyright (C) 2014-2024 ABINIT group (AL)
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

#ifdef HAVE_GPU_MARKERS
 use m_nvtx
#endif

 implicit none

 private

 ! Use these routines in order: first call init, then call make_gemm_nonlop for each k point,
 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: init_gemm_nonlop
 public :: destroy_gemm_nonlop
 public :: make_gemm_nonlop
 public :: prep_projectors
 public :: prep_dprojectors

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
   real(dp), allocatable :: d2projs_r(:, :, :)
   ! (1, npw, nprojs*ngrads)
   real(dp), allocatable :: d2projs_i(:, :, :)
   ! (1, npw, nprojs*ngrads)
   integer :: idir
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

 integer, save, public :: gemm_nonlop_nblocks = 1
 ! Public variable indicating in how many blocks of MPI tasks should the projs arrays be ditributed.
 ! Default size 1 indicates no distribution at all.

 integer, save, public :: gemm_nonlop_block_comm = xmpi_comm_null
 ! MPI communicator for MPI tasks processing the same gemm_nonlop block for projs array distribution

 integer, save, public :: gemm_nonlop_choice = -1

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

 integer, save, public :: current_ikpt_in_gpu=-1

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
  integer :: ikpt

! *************************************************************************

  gemm_nonlop_kpt(:)%nprojs = -1
  gemm_nonlop_kpt(:)%ngrads = -1
  gemm_nonlop_kpt(:)%ngrads2 = -1
  gemm_nonlop_kpt(:)%choice = -1
  gemm_nonlop_kpt(:)%idir = -1
  current_ikpt_in_gpu = -1

  if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    gemm_nonlop_kpt_gpu(:)%npw = -1
    gemm_nonlop_kpt_gpu(:)%nprojs = -1
    gemm_nonlop_gpu_data % allocated = .false.
#endif
  end if


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
  integer :: ikpt

! *************************************************************************

  call free_gemm_nonlop_ikpt(1,gpu_option)
  call free_gemm_nonlop_ikpt(2,gpu_option)

 end subroutine destroy_gemm_nonlop
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
     !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%projs) IF(gpu_option==ABI_GPU_OPENMP)
     ABI_FREE(gemm_nonlop_kpt(ik)%projs)
   end if
   if (allocated(gemm_nonlop_kpt(ik)%projs_r)) then
     !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%projs_r) IF(gpu_option==ABI_GPU_OPENMP)
     ABI_FREE(gemm_nonlop_kpt(ik)%projs_r)
   end if
   if (allocated(gemm_nonlop_kpt(ik)%projs_i)) then
     !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%projs_i) IF(gpu_option==ABI_GPU_OPENMP)
   ABI_FREE(gemm_nonlop_kpt(ik)%projs_i)
   end if
   gemm_nonlop_kpt(ik)%nprojs = -1
   if(gemm_nonlop_kpt(ik)%ngrads /= -1) then
     if (allocated(gemm_nonlop_kpt(ik)%dprojs)) then
       !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%dprojs) IF(gpu_option==ABI_GPU_OPENMP)
       ABI_FREE(gemm_nonlop_kpt(ik)%dprojs)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%dprojs_r)) then
       !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%dprojs_r) IF(gpu_option==ABI_GPU_OPENMP)
       ABI_FREE(gemm_nonlop_kpt(ik)%dprojs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%dprojs_i)) then
       !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%dprojs_i) IF(gpu_option==ABI_GPU_OPENMP)
       ABI_FREE(gemm_nonlop_kpt(ik)%dprojs_i)
     end if
     gemm_nonlop_kpt(ik)%ngrads = -1
   end if
   if(gemm_nonlop_kpt(ik)%ngrads2 /= -1) then
     if (allocated(gemm_nonlop_kpt(ik)%d2projs)) then
       !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%d2projs) IF(gpu_option==ABI_GPU_OPENMP)
       ABI_FREE(gemm_nonlop_kpt(ik)%d2projs)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%d2projs_r)) then
       !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%d2projs_r) IF(gpu_option==ABI_GPU_OPENMP)
       ABI_FREE(gemm_nonlop_kpt(ik)%d2projs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ik)%d2projs_i)) then
       !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%d2projs_i) IF(gpu_option==ABI_GPU_OPENMP)
       ABI_FREE(gemm_nonlop_kpt(ik)%d2projs_i)
     end if
     gemm_nonlop_kpt(ik)%ngrads2 = -1
   end if
 end if
 current_ikpt_in_gpu = -1

 end subroutine free_gemm_nonlop_ikpt
!!***

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt(ik)

  integer,intent(in) :: ik
#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  call free_ompgpu_current_ikpt_projs(ik)
  call free_ompgpu_current_ikpt_dprojs(ik)

#endif
 end subroutine free_ompgpu_current_ikpt

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt_projs(ik)

  integer,intent(in) :: ik
#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%projs))) then
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%projs)
  end if
  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%projs_r))) then
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%projs_r)
  end if

#endif
 end subroutine free_ompgpu_current_ikpt_projs

!----------------------------------------------------------------------


 subroutine free_ompgpu_current_ikpt_dprojs(ik)

  integer,intent(in) :: ik
#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%dprojs))) then
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%dprojs)
  end if
  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%dprojs_i))) then
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%dprojs_i)
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%dprojs_r)
  end if

  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%d2projs))) then
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%d2projs)
  end if
  if(xomp_target_is_present(c_loc(gemm_nonlop_kpt(ik)%d2projs_i))) then
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%d2projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:gemm_nonlop_kpt(ik)%d2projs_r)
  end if

#endif
 end subroutine free_ompgpu_current_ikpt_dprojs

!----------------------------------------------------------------------


!!****f* m_gemm_nonlop_projectors/make_gemm_nonlop
!! NAME
!! make_gemm_nonlop
!!
!! FUNCTION
!! Build the gemm_nonlop array
!!
!! INPUTS
!!
!! SOURCE
 subroutine make_gemm_nonlop(ikpt,signs,choice,npw_k,npw_kp,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol, &
&                            ffnl_k,ph3d_k,kpt_k,kg_k,kpg_k,&
&                            ffnl_kp,ph3d_kp,kpt_kp,kg_kp,kpg_kp,&
&                            select_k, &
&                            idir_pert,gpu_option) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw_k, npw_kp, lmnmax, ntypat,select_k
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:), kg_kp(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  integer, intent(in) :: signs,choice
  integer, intent(in), optional :: idir_pert,gpu_option
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:),ffnl_kp(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:),ph3d_kp(:,:,:)
  real(dp), intent(in) :: kpt_k(:),kpt_kp(:)
  real(dp), intent(in), target :: kpg_k(:,:),kpg_kp(:,:)

  integer :: nprojs,ndprojs,ndgxdt,nd2gxdt,npw

  integer :: rank, nprocs, ierr, itypat, i, ipw, ik
  integer :: nprojs_blk, nprojs_last_blk, nprojs_my_blk
  integer :: idir_beg,idir_end,idir_pert_,gpu_option_
  logical :: is_last_rank,compute_dprojs
  real(dp),allocatable :: dprojs_tmp(:,:,:),dprojs_r_tmp(:,:,:),dprojs_i_tmp(:,:,:)

! *************************************************************************

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxStartRange("make_gemm_nonlop")
#endif
  ABI_CHECK(size(ph3d_k)>0,'nloalg(2)<0 not compatible with use_gemm_nonlop=1!')
!  ABI_CHECK((.not.my_compute_grad_strain).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
!  ABI_CHECK((.not.my_compute_grad_atom).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
  idir_pert_=-1; if(present(idir_pert)) idir_pert_=idir_pert
  gpu_option_=ABI_GPU_DISABLED; if(present(gpu_option)) gpu_option_=gpu_option

  if(choice>1) then
    ABI_CHECK(signs==2.or. choice==2 .or. choice==3 .or. choice==23 .or. choice==4 .or. choice==54 .or. choice==55,'signs/=2 and idir_pert not compatible with GEMM nonlop.')
    ABI_CHECK(gpu_option_/=ABI_GPU_LEGACY,'CUDA GEMM nonlop not compatible with respfn workloads.')
    ABI_CHECK(gpu_option_/=ABI_GPU_KOKKOS,'KOKKOS GEMM nonlop not compatible with respfn workloads.')
  end if

  ! build nprojs, ndgxdt
  nprojs = 0 ; ndgxdt = 0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do

!Define dimensions of projected scalars
  ndgxdt=-1;nd2gxdt=-1
  if (choice==2) then
    if (signs==1) ndgxdt=3
    if (signs==2) ndgxdt=1
  end if
  if (choice==3) then
    if (signs==1) ndgxdt=6
    if (signs==2) ndgxdt=1
  end if
  if (choice==23) then
    if (signs==1) ndgxdt=9
  end if
  if (choice==5) then
    !if(signs==1) ndgxdt=3
    if(signs==2) ndgxdt=1
  end if
  if (choice==51) then
    !if(signs==1) ndgxdt=3
    if(signs==2) ndgxdt=1
  end if
  if (choice==4) then
    if(signs==1) ndgxdt=3
    if(signs==1) nd2gxdt=6
  end if
  if (choice==54) then
    if(signs==1) ndgxdt=6
    if(signs==1) nd2gxdt=9

    !if(signs==2) ndgxdt=1
    !if(signs==2) nd2gxdt=1
  end if
  if (choice==55) then
    if(signs==1) ndgxdt=9
    if(signs==1) nd2gxdt=18
  end if
  !if (choice==6) then
  !  if(signs==1) ndgxdt=9
  !  if(signs==1) nd2gxdt=54
  !end if




  if(current_ikpt_in_gpu/=ikpt) then
    call free_gemm_nonlop_ikpt(1, gpu_option_)
    call free_gemm_nonlop_ikpt(2, gpu_option_)
  end if

  do ik=1,2
  if(ik==1 .and. select_k==KPRIME_H_KPRIME) cycle
  if(ik==2 .and. select_k==K_H_K) cycle
  npw=npw_k; if(ik==2) npw=npw_kp
  compute_dprojs = .false.

  if(nprojs/=gemm_nonlop_kpt(ik)%nprojs) then
    call free_gemm_nonlop_ikpt(ik, gpu_option_)

    if(istwf_k <= 1) then
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs, (2, npw, nprojs))
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs_r, (1,1,1))
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs_i, (1,1,1))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%projs) IF(gpu_option_==ABI_GPU_OPENMP)
    else
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs_r, (1, npw, nprojs))
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs_i, (1, npw, nprojs))
      ABI_MALLOC(gemm_nonlop_kpt(ik)%projs, (1, 1, 1))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%projs_r) IF(gpu_option_==ABI_GPU_OPENMP)
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%projs_i) IF(gpu_option_==ABI_GPU_OPENMP)
    end if

    if(ik==1) then
    call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
    &                    ucvol,ffnl_k,ph3d_k,&
    &                    nprojs,choice,gpu_option_,&
    &                    gemm_nonlop_kpt(ik)%projs,&
    &                    gemm_nonlop_kpt(ik)%projs_r,gemm_nonlop_kpt(ik)%projs_i)
    end if
    if(ik==2) then
    call prep_projectors(npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,&
    &                    ucvol,ffnl_kp,ph3d_kp,&
    &                    nprojs,choice,gpu_option_,&
    &                    gemm_nonlop_kpt(ik)%projs,&
    &                    gemm_nonlop_kpt(ik)%projs_r,gemm_nonlop_kpt(ik)%projs_i)
    end if
  end if


  if(ndgxdt>0) then
    if(nprojs/=gemm_nonlop_kpt(ik)%nprojs .or. ndgxdt /= gemm_nonlop_kpt(ik)%ngrads &
    &    .or. nd2gxdt /=  gemm_nonlop_kpt(ik)%ngrads2) then
      if(gpu_option_ == ABI_GPU_OPENMP) call free_ompgpu_current_ikpt_dprojs(ik)
        if(allocated(gemm_nonlop_kpt(ik)%dprojs)) ABI_FREE(gemm_nonlop_kpt(ik)%dprojs)
        if(allocated(gemm_nonlop_kpt(ik)%dprojs_r)) ABI_FREE(gemm_nonlop_kpt(ik)%dprojs_r)
        if(allocated(gemm_nonlop_kpt(ik)%dprojs_i)) ABI_FREE(gemm_nonlop_kpt(ik)%dprojs_i)
        if(allocated(gemm_nonlop_kpt(ik)%d2projs)) ABI_FREE(gemm_nonlop_kpt(ik)%d2projs)
        if(allocated(gemm_nonlop_kpt(ik)%d2projs_r)) ABI_FREE(gemm_nonlop_kpt(ik)%d2projs_r)
        if(allocated(gemm_nonlop_kpt(ik)%d2projs_i)) ABI_FREE(gemm_nonlop_kpt(ik)%d2projs_i)
        gemm_nonlop_kpt(ik)%ngrads = -1
        gemm_nonlop_kpt(ik)%ngrads2 = -1
      if(istwf_k <= 1) then
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs, (2, npw, nprojs*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs_r, (1, 1, 1))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs_i, (1, 1, 1))
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%dprojs) IF(gpu_option_==ABI_GPU_OPENMP)
        if(nd2gxdt>0) then
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs, (2, npw, nprojs*nd2gxdt))
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs_r, (1, 1, 1))
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs_i, (1, 1, 1))
          !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%d2projs) IF(gpu_option_==ABI_GPU_OPENMP)
        end if
      else
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs_r, (1, npw, nprojs*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs_i, (1, npw, nprojs*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%dprojs, (1, 1, 1))
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%dprojs_r) IF(gpu_option_==ABI_GPU_OPENMP)
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%dprojs_i) IF(gpu_option_==ABI_GPU_OPENMP)
        if(nd2gxdt>0) then
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs_r, (1, npw, nprojs*nd2gxdt))
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs_i, (1, npw, nprojs*nd2gxdt))
          ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs, (1, 1, 1))
          !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%d2projs_r) IF(gpu_option_==ABI_GPU_OPENMP)
          !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%d2projs_i) IF(gpu_option_==ABI_GPU_OPENMP)
        end if
      end if
      if(nd2gxdt==0) then
        ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs, (1, 1, 1))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs_r, (1, 1, 1))
        ABI_MALLOC(gemm_nonlop_kpt(ik)%d2projs_i, (1, 1, 1))
        !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ik)%d2projs) IF(gpu_option_==ABI_GPU_OPENMP)
      end if
      compute_dprojs=.true.
    end if

    if(choice/=gemm_nonlop_kpt(ik)%choice .or. idir_pert_/=gemm_nonlop_kpt(ik)%idir &
    &  .or. compute_dprojs) then
      if(ik==1) then
      call prep_dprojectors(npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k,&
      &                    ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k,&
      &                    nprojs,ndgxdt,nd2gxdt,choice,signs,idir_pert_,gpu_option_,&
      &                    gemm_nonlop_kpt(ik)%projs,&
      &                    gemm_nonlop_kpt(ik)%projs_r,gemm_nonlop_kpt(ik)%projs_i,&
      &                    gemm_nonlop_kpt(ik)%dprojs,&
      &                    gemm_nonlop_kpt(ik)%dprojs_r,gemm_nonlop_kpt(ik)%dprojs_i,&
      &                    gemm_nonlop_kpt(ik)%d2projs,&
      &                    gemm_nonlop_kpt(ik)%d2projs_r,gemm_nonlop_kpt(ik)%d2projs_i)
      end if
      if(ik==2) then
      call prep_dprojectors(npw,lmnmax,ntypat,indlmn,kg_kp,nattyp,istwf_k,&
      &                    ucvol,ffnl_kp,ph3d_kp,kpt_kp,kpg_kp,&
      &                    nprojs,ndgxdt,nd2gxdt,choice,signs,idir_pert_,gpu_option_,&
      &                    gemm_nonlop_kpt(ik)%projs,&
      &                    gemm_nonlop_kpt(ik)%projs_r,gemm_nonlop_kpt(ik)%projs_i,&
      &                    gemm_nonlop_kpt(ik)%dprojs,&
      &                    gemm_nonlop_kpt(ik)%dprojs_r,gemm_nonlop_kpt(ik)%dprojs_i,&
      &                    gemm_nonlop_kpt(ik)%d2projs,&
      &                    gemm_nonlop_kpt(ik)%d2projs_r,gemm_nonlop_kpt(ik)%d2projs_i)
      end if
    end if
  end if

  if (nprojs>0) gemm_nonlop_kpt(ik)%nprojs = nprojs
  if (ndgxdt>0) gemm_nonlop_kpt(ik)%ngrads = ndgxdt
  if (nd2gxdt>0) gemm_nonlop_kpt(ik)%ngrads2 = nd2gxdt
  if (ndgxdt>0) gemm_nonlop_kpt(ik)%choice = choice
  if (ndgxdt>0) gemm_nonlop_kpt(ik)%idir = idir_pert_
  end do
  current_ikpt_in_gpu = ikpt


  !!!!! GPU stuff

  if(gpu_option_ == ABI_GPU_LEGACY .or. gpu_option_ == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    do ik=1,2
      if(ik==1 .and. select_k==KPRIME_H_KPRIME) cycle
      if(ik==2 .and. select_k==K_H_K) cycle
      gemm_nonlop_kpt_gpu%npw    = npw
      gemm_nonlop_kpt_gpu%nprojs = nprojs

#ifdef DEBUG_VERBOSE_GPU
      if(xmpi_comm_rank(xmpi_world) == 0) then
        call check_gpu_mem("make_gemm_nonlop begin")
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
        call check_gpu_mem("make_gemm_nonlop end  ")
      end if
#endif

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

    end do
#endif
  end if
  gemm_nonlop_choice=choice
#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxEndRange()
#endif

 end subroutine make_gemm_nonlop
!!***


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
 &                          ucvol,ffnl_k,ph3d_k,&
 &                          nprojs,choice,gpu_option,&
 &                          projs,projs_r,projs_i)

  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  integer, intent(in) :: nprojs,choice,gpu_option
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp),intent(inout), target ::   projs(:,:,:)
  real(dp),intent(inout), target :: projs_r(:,:,:)
  real(dp),intent(inout), target :: projs_i(:,:,:)

  logical :: parity
  integer :: il, ipw, idir, idir1, idir2, ffnl_dir, dimffnl
  integer :: itypat, ilmn, nlmn, ia, iaph3d, shift
  integer :: lmn_beg,ibeg,iend,shift_do,nlmn_o,lmn_grad_beg
  real(dp):: wt
  real(dp),allocatable, target :: atom_projs(:,:,:), temp(:,:)

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxStartRange("prep_projectors")
#endif

  if(gpu_option==ABI_GPU_OPENMP) then

    if(istwf_k <= 1) then
      if(.not. xomp_target_is_present(c_loc(projs)) .and. istwf_k==1) ABI_BUG("Stopii......")
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
  dimffnl = size(ffnl_k, dim=2)

  ABI_MALLOC(atom_projs, (2, npw, lmnmax))
  !$OMP TARGET ENTER DATA MAP(alloc:atom_projs) IF(gpu_option==ABI_GPU_OPENMP)

  ABI_MALLOC(temp, (npw, lmnmax))
  !$OMP TARGET ENTER DATA MAP(alloc:temp) IF(gpu_option==ABI_GPU_OPENMP)

  !!$OMP TARGET ENTER DATA MAP(to:ffnl_k,ph3d_k)

  shift = 0
  lmn_beg = 1

  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)
    nlmn_o = nlmn

    do ia = 1, nattyp(itypat)

      !! build atom_projs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l

      ! start from 4pi/sqrt(ucvol)*ffnl
      ! atom_projs(1, :, lmn_beg:nlmn) = four_pi/sqrt(ham%ucvol) * ham%ffnl_k(:, 1, lmn_beg:nlmn)
      ! TODO vectorize (DCOPY with stride)
      !if(gpu_option==ABI_GPU_OPENMP) then
      !  !$OMP TARGET DATA USE_DEVICE_PTR(atom_projs)
      !  call gpu_memset(c_loc(atom_projs), 0, int(2,c_size_t)*npw*lmnmax*dp)
      !  !$OMP END TARGET DATA
      !else
        atom_projs(:,:,:) = zero
      !end if
      !!$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(atom_projs,ffnl_k) &
      !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
      do ipw=1, npw
        atom_projs(1,ipw, 1:nlmn_o) = wt * ffnl_k(ipw, 1, 1:nlmn_o, itypat)
      end do

      !!$OMP TARGET UPDATE FROM(atom_projs)

      ! multiply by (-i)^l
      do ilmn=1,nlmn_o
        il=mod(indlmn(1,ilmn, itypat),4);
        parity=(mod(il,2)==0)
        if (il>1) then
          ! multiply by -1
          atom_projs(:,:,ilmn) = -atom_projs(:,:,ilmn)
        end if
        if(.not. parity) then
          ! multiply by -i
          temp(:,ilmn) = atom_projs(2,:,ilmn)
          atom_projs(2,:,ilmn) = -atom_projs(1,:,ilmn)
          atom_projs(1,:,ilmn) =  temp(:,ilmn)
        end if
      end do

      ! multiply by conj(ph3d)
      !!$OMP TARGET PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:temp,atom_projs,ph3d_k) &
      !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
      do ilmn=1,nlmn_o
        do ipw=1,npw
          temp(ipw,ilmn) = atom_projs(1, ipw, ilmn)
          atom_projs(1, ipw, ilmn) = atom_projs(1, ipw, ilmn) * ph3d_k(1, ipw, iaph3d) &
          &                      + atom_projs(2, ipw, ilmn) * ph3d_k(2, ipw, iaph3d)
          atom_projs(2, ipw, ilmn) = atom_projs(2, ipw, ilmn) * ph3d_k(1, ipw, iaph3d) &
          &                      - temp(ipw,ilmn)           * ph3d_k(2, ipw, iaph3d)
        end do
      end do

      !!$OMP TARGET UPDATE FROM(atom_projs) if(choice>1)
      !! atom_projs is built, copy to projs

      if(gpu_option==ABI_GPU_OPENMP) then
        if(istwf_k <= 1) then
          !!$OMP TARGET DATA USE_DEVICE_PTR(projs)
          !call copy_on_gpu(atom_projs(:, :, lmn_beg), c_loc(projs(:,:,shift+1)), int(2,c_size_t)*npw*(nlmn-(lmn_beg-1))*dp)
          !!$OMP END TARGET DATA
          !$OMP TARGET UPDATE TO(atom_projs)
          !$OMP TARGET DATA USE_DEVICE_PTR(projs,atom_projs)
          call copy_gpu_to_gpu(c_loc(projs(:,:,shift+1)), c_loc(atom_projs(:, :, lmn_beg)), int(2,c_size_t)*npw*(nlmn-(lmn_beg-1))*dp)
          !$OMP END TARGET DATA
          projs(1:2, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(:, :, lmn_beg:nlmn)
          !!$OMP TARGET PARALLEL DO COLLAPSE(2) PRIVATE(ipw,ilmn) MAP(to:projs,atom_projs)
          !do ilmn=lmn_beg,nlmn
          !  do ipw=1,npw
          !    projs(:, ipw, shift+ilmn) = atom_projs(:, ipw, ilmn)
          !  end do
          !end do
          !!$OMP TARGET UPDATE FROM(projs)
        else ! istwf_k>1
          !$OMP TARGET UPDATE TO(atom_projs)
          !$OMP TARGET MAP(to:atom_projs,projs_r,projs_i)
          projs_r(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(1, :, lmn_beg:nlmn)
          projs_i(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(2, :, lmn_beg:nlmn)
          !$OMP END TARGET
          projs_r(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(1, :, lmn_beg:nlmn)
          projs_i(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(2, :, lmn_beg:nlmn)
        end if
      else
        if(istwf_k <= 1) then
          projs(1:2, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(:, :, lmn_beg:nlmn)
        else ! istwf_k>1
          projs_r(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(1, :, lmn_beg:nlmn)
          projs_i(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(2, :, lmn_beg:nlmn)
        end if
      end if

      iaph3d = iaph3d + 1
      shift = shift + nlmn

    end do
  end do

  !$OMP TARGET EXIT DATA MAP(delete:atom_projs,temp) IF(gpu_option==ABI_GPU_OPENMP)
  ABI_FREE(atom_projs)
  ABI_FREE(temp)
  !!$OMP TARGET EXIT DATA MAP(delete:ffnl_k,ph3d_k)

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxEndRange()
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
 subroutine prep_dprojectors(npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k,&
 &                          ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k,&
 &                          nprojs,ngrads,ngrads2,choice,signs,idir_pert,gpu_option,&
 &                          projs,projs_r,projs_i,&
 &                          dprojs,dprojs_r,dprojs_i,&
 &                          d2projs,d2projs_r,d2projs_i)

  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  integer, intent(in) :: nprojs,ngrads,ngrads2,choice,signs,idir_pert,gpu_option
  real(dp), intent(in) :: ucvol
  real(dp), intent(in), target :: ffnl_k(:,:,:,:)
  real(dp), intent(in), target :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)
  real(dp),intent(inout), target :: projs  (:,:,:)  ,dprojs(:,:,:)  ,d2projs(:,:,:)
  real(dp),intent(inout), target :: projs_r(:,:,:),dprojs_r(:,:,:),d2projs_r(:,:,:)
  real(dp),intent(inout), target :: projs_i(:,:,:),dprojs_i(:,:,:),d2projs_i(:,:,:)

  logical,allocatable :: parity(:)
  logical :: map_ffnl_k,map_ph3d_k
  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer :: ndprojs, nd2projs
  integer :: il, ipw, idir, idir1, idir2, nkpg_local, ffnl_dir, dimffnl
  integer :: itypat, ilmn, nlmn, ia, iaph3d, igrad, shift, shift_grad, shift_grad2
  integer :: lmn_beg,ibeg,iend,shift_do,nlmn_o,lmn_grad_beg
  real(dp), parameter :: two_pi2=two_pi*two_pi
  real(dp):: wt,tmp
  real(dp),allocatable, target :: atom_projs(:,:,:), atom_dprojs(:,:,:,:), atom_d2projs(:,:,:,:), temp(:,:), scal(:)
  real(dp),pointer :: kpg(:,:)

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxStartRange("prep_dprojectors")
#endif

  !if(gpu_option==ABI_GPU_OPENMP) then

  !  if(istwf_k <= 1) then
  !    call gpu_set_to_zero(dprojs,2*npw*nprojs*ngrads)
  !    if(ngrads2 > 0) then
  !      call gpu_set_to_zero(d2projs,2*npw*nprojs*ngrads2)
  !    end if
  !  else
  !    call gpu_set_to_zero(dprojs_r,npw*nprojs*ngrads)
  !    call gpu_set_to_zero(dprojs_i,npw*nprojs*ngrads)
  !    if(ngrads2 > 0) then
  !      call gpu_set_to_zero(d2projs_r,npw*nprojs*ngrads2)
  !      call gpu_set_to_zero(d2projs_i,npw*nprojs*ngrads2)
  !    end if
  !  end if

  !else

    if(istwf_k <= 1) then
      dprojs(:,:,:) = zero
      if(ngrads2 > 0) then
        d2projs(:,:,:) = zero
      end if
    else
      dprojs_r(:,:,:) = zero
      dprojs_i(:,:,:) = zero
      if(ngrads2 > 0) then
        d2projs_r(:,:,:) = zero
        d2projs_i(:,:,:) = zero
      end if
    end if

  !end if

  iaph3d = 1
  wt=four_pi/sqrt(ucvol)
  dimffnl = size(ffnl_k, dim=2)
  ffnl_dir=1; if(dimffnl>2) ffnl_dir=idir_pert

  ndprojs = 0
  if (signs==1 .and. (choice==3 .or. choice==23 .or. choice==54 .or. choice==55)) then
    ndprojs = 3
  else if(signs==2 .and. (choice==5 .or. choice==51 .or. choice==3)) then
    ndprojs = 1
  end if
  nd2projs = 0
  if(signs==1 .and. choice==54) then
    nd2projs = 3
  else if(signs==1 .and. choice==55) then
    nd2projs = 6
  end if


  if (ndprojs>0) then
    ABI_MALLOC(atom_dprojs, (2, npw, ndprojs, lmnmax))
    !!$OMP TARGET ENTER DATA MAP(alloc:atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
    if(ngrads2>0) then
      ABI_MALLOC(atom_d2projs, (2, npw, nd2projs, lmnmax))
      !!$OMP TARGET ENTER DATA MAP(alloc:atom_d2projs) IF(gpu_option==ABI_GPU_OPENMP)
    end if
  end if

  ABI_MALLOC(temp, (npw, lmnmax))
  ABI_MALLOC(parity, (lmnmax))
  ABI_MALLOC(scal, (lmnmax))
  !!$OMP TARGET ENTER DATA MAP(alloc:temp,scal) IF(gpu_option==ABI_GPU_OPENMP)

  ! Compute (k+G) vectors if needed
  nkpg_local=0
  if ((choice==2 .or. choice==3 .or. choice==23 .or. choice==54 .or. choice==55).and.size(kpg_k)==0) then
    nkpg_local=3
    ABI_MALLOC(kpg,(npw,nkpg_local))
    call mkkpg(kg_k,kpg,kpt_k,nkpg_local,npw)
  else if (choice==4.and.size(kpg_k,2)<9) then
    nkpg_local=9
    ABI_MALLOC(kpg,(npw,nkpg_local))
    call mkkpg(kg_k,kpg,kpt_k,nkpg_local,npw)
  else
    kpg => kpg_k
  end if

  !!$OMP TARGET ENTER DATA MAP(to:kpg) IF(gpu_option==ABI_GPU_OPENMP)
  !map_ph3d_k=.false.; map_ffnl_k=.false.
  !if(.not. xomp_target_is_present(c_loc(ffnl_k))) map_ffnl_k = .true.
  !if(.not. xomp_target_is_present(c_loc(ph3d_k))) map_ph3d_k = .true.
  !!$OMP TARGET ENTER DATA MAP(to:ffnl_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl_k)
  !!$OMP TARGET ENTER DATA MAP(to:ph3d_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d_k)

  shift = 0 ; shift_grad = 0; shift_grad2 = 0
  lmn_beg = 1

  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)
    nlmn_o = nlmn

    do ia = 1, nattyp(itypat)

      !! build atom_dprojs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l

      ! start from 4pi/sqrt(ucvol)*ffnl
      !if(gpu_option==ABI_GPU_OPENMP) then
      !  if (ndprojs>0) then
      !    call gpu_set_to_zero(atom_dprojs,2*npw*ndprojs*lmnmax)
      !    if(ngrads2>0) then
      !      call gpu_set_to_zero(atom_d2projs,2*npw*nd2projs*lmnmax)
      !    end if
      !  end if
      !else
        if (ndprojs>0) atom_dprojs(:,:,:,:) = zero
        if (nd2projs>0) atom_d2projs(:,:,:,:) = zero
      !end if
      if (signs==1 .and. (choice==3 .or. choice==23 .or. choice==54 .or. choice==55)) then
        !!$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_dprojs,ffnl_k) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_dprojs(1,ipw, 1:ndprojs, 1:nlmn_o) = wt * ffnl_k(ipw, 2:ndprojs+1, 1:nlmn_o, itypat)
        end do
      end if
      if (signs==2 .and. (choice==3 .or. choice==5 .or. choice==51)) then
        !!$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_dprojs,ffnl_k) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_dprojs(1,ipw, 1, 1:nlmn_o) = wt * ffnl_k(ipw, 1+ffnl_dir, 1:nlmn_o, itypat)
        end do
      end if
      if(signs==1 .and. choice==54) then
        !!$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl_k) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_d2projs(1,ipw, 1:nd2projs, 1:nlmn_o) = wt * ffnl_k(ipw, 2:nd2projs+1, 1:nlmn_o, itypat)
        end do
      else if(signs==1 .and. choice==55) then
        !!$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl_k) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_d2projs(1,ipw, 1:nd2projs, 1:nlmn_o) = wt * ffnl_k(ipw, 5:nd2projs+4, 1:nlmn_o, itypat)
        end do
      end if

      ! multiply by (-i)^l
      if (ndprojs>0) then
        do ilmn=1,nlmn_o
          il=mod(indlmn(1,ilmn, itypat),4);
          parity(ilmn)=(mod(il,2)==0)
          scal(ilmn)=1; if(il>1) scal(ilmn)=-1
        end do
        ! multiply by -1
        if(gpu_option==ABI_GPU_OPENMP) then
          !!$OMP TARGET UPDATE FROM(atom_dprojs) if(gpu_option==ABI_GPU_OPENMP)
          !!$OMP TARGET UPDATE FROM(atom_d2projs) if(ngrads2>0 .and. gpu_option==ABI_GPU_OPENMP)
          !!$OMP TARGET UPDATE TO(scal)
          !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !!$OMP& PRIVATE(idir,ipw,ilmn) MAP(to:atom_dprojs,scal) &
          !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do ilmn=1,nlmn_o
            do idir=1,ndprojs
              do ipw=1,npw
                atom_dprojs(1,ipw,idir,ilmn) = atom_dprojs(1,ipw,idir,ilmn) * scal(ilmn)
                atom_dprojs(2,ipw,idir,ilmn) = atom_dprojs(2,ipw,idir,ilmn) * scal(ilmn)
              end do
            end do
          end do
          if(nd2projs>0) then
            !!$OMP TARGET UPDATE TO(scal)
            !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
            !!$OMP& PRIVATE(idir,ipw,ilmn) MAP(to:atom_d2projs,scal) &
            !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
            do ilmn=1,nlmn_o
              do idir=1,ndprojs
                do ipw=1,npw
                  atom_d2projs(1,ipw,idir,ilmn) = atom_d2projs(1,ipw,idir,ilmn) * scal(ilmn)
                  atom_d2projs(2,ipw,idir,ilmn) = atom_d2projs(2,ipw,idir,ilmn) * scal(ilmn)
                end do
              end do
            end do
          end if
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
        !!$OMP TARGET TEAMS DISTRIBUTE MAP(to:atom_dprojs,parity) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=1,nlmn_o
          if(.not. parity(ilmn)) then
            !!$OMP PARALLEL DO PRIVATE(idir,ipw,tmp) COLLAPSE(2)
            do idir=1,ndprojs
              do ipw=1,npw
                tmp = atom_dprojs(2,ipw,idir,ilmn)
                atom_dprojs(2,ipw,idir,ilmn) = -atom_dprojs(1,ipw,idir,ilmn)
                atom_dprojs(1,ipw,idir,ilmn) =  tmp
              end do
            end do
          end if
        end do
        if(ngrads2>0) then
          !!$OMP TARGET TEAMS DISTRIBUTE MAP(to:atom_d2projs,parity) &
          !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do ilmn=1,nlmn_o
            if(.not. parity(ilmn)) then
              !!$OMP PARALLEL DO PRIVATE(idir,ipw,tmp) COLLAPSE(2)
              do idir=1,nd2projs
                do ipw=1,npw
                  tmp = atom_d2projs(2,ipw,idir,ilmn)
                  atom_d2projs(2,ipw,idir,ilmn) =  atom_d2projs(1,ipw,idir,ilmn)
                  atom_d2projs(1,ipw,idir,ilmn) = -tmp
                end do
              end do
            else
              !!$OMP PARALLEL DO PRIVATE(idir,ipw,tmp) COLLAPSE(2)
              do idir=1,nd2projs
                do ipw=1,npw
                  atom_d2projs(1,ipw,idir,ilmn) = -atom_d2projs(1,ipw,idir,ilmn)
                  atom_d2projs(2,ipw,idir,ilmn) = -atom_d2projs(2,ipw,idir,ilmn)
                end do
              end do
            end if
          end do
        end if
        !!$OMP TARGET UPDATE TO(atom_dprojs) if(gpu_option==ABI_GPU_OPENMP)
        !!$OMP TARGET UPDATE TO(atom_d2projs) if(ngrads2>0 .and. gpu_option==ABI_GPU_OPENMP)
      end if

      ! multiply by conj(ph3d)
      if (ndprojs>0) then
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,idir,ipw,tmp) COLLAPSE(3) MAP(to:atom_dprojs,ph3d_k) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=1,nlmn_o
          do idir=1,ndprojs
            do ipw=1,npw
              tmp = atom_dprojs(1, ipw, idir,ilmn)
              atom_dprojs(1, ipw, idir,ilmn) = atom_dprojs(1, ipw, idir,ilmn) * ph3d_k(1, ipw, iaph3d) &
              &                              + atom_dprojs(2, ipw, idir,ilmn) * ph3d_k(2, ipw, iaph3d)
              atom_dprojs(2, ipw, idir,ilmn) = atom_dprojs(2, ipw, idir,ilmn) * ph3d_k(1, ipw, iaph3d) &
              &                              - tmp                 * ph3d_k(2, ipw, iaph3d)
            end do
          end do
        end do
      end if
      if (nd2projs>0) then
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,idir,ipw,tmp) COLLAPSE(3) MAP(to:atom_d2projs,ph3d_k) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=1,nlmn_o
          do idir=1,nd2projs
            do ipw=1,npw
              tmp = atom_d2projs(1, ipw, idir,ilmn)
              atom_d2projs(1, ipw, idir,ilmn) = atom_d2projs(1, ipw, idir,ilmn) * ph3d_k(1, ipw, iaph3d) &
              &                               + atom_d2projs(2, ipw, idir,ilmn) * ph3d_k(2, ipw, iaph3d)
              atom_d2projs(2, ipw, idir,ilmn) = atom_d2projs(2, ipw, idir,ilmn) * ph3d_k(1, ipw, iaph3d) &
              &                               - tmp                 * ph3d_k(2, ipw, iaph3d)
            end do
          end do
        end do
      end if

      !! Handling dprojs

      if(signs==1 .and. (choice==3 .or. choice==23 .or. choice==55)) then
        if(istwf_k <= 1) then
          !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw,idir,idir1,idir2) COLLAPSE(3) MAP(to:atom_dprojs,dprojs,kpg) &
          !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do ilmn=lmn_beg,nlmn
            do idir=1,6
              do ipw=1,npw
                idir1=alpha(idir);idir2=beta(idir)
                dprojs(1, ipw, shift_grad+(ilmn-1)*ngrads+idir) = &
                &     -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))
                dprojs(2, ipw, shift_grad+(ilmn-1)*ngrads+idir) = &
                &     -half*(atom_dprojs(2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(2, ipw, idir2, ilmn)*kpg(ipw,idir1))
              end do
            end do
          end do
        else ! istwf_k>1
          !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3)&
          !!$OMP& PRIVATE(ilmn,ipw,idir,idir1,idir2) MAP(to:atom_dprojs,dprojs_r,dprojs_i,kpg) &
          !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do idir=1,6
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                idir1=alpha(idir);idir2=beta(idir)
                dprojs_r(1, ipw, shift_grad+(ilmn-1)*ngrads+idir) = &
                &     -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))

                dprojs_i(1, ipw, shift_grad+(ilmn-1)*ngrads+idir) = &
                &     -half*(atom_dprojs(2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(2, ipw, idir2, ilmn)*kpg(ipw,idir1))
              end do
            end do
          end do
        end if
      end if


      if(signs==1 .and. (choice==2 .or. choice==23 .or. choice==4 .or. choice==54)) then
        igrad=0; if(choice==23) igrad=6
        if(istwf_k <= 1) then
          !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw,idir) COLLAPSE(3) MAP(to:projs,dprojs,kpg) &
          !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do ilmn=lmn_beg,nlmn
            do idir=1,3
              do ipw=1,npw
                dprojs(1, ipw, shift_grad+(ilmn-1)*ngrads+igrad+idir) = &
                &     +projs(2, ipw, shift+ilmn)*kpg(ipw,idir)*two_pi
                dprojs(2, ipw, shift_grad+(ilmn-1)*ngrads+igrad+idir) = &
                &     -projs(1, ipw, shift+ilmn)*kpg(ipw,idir)*two_pi
              end do
            end do
          end do
        else ! istwf_k>1
          !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw,idir) COLLAPSE(3) MAP(to:projs_r,projs_i,dprojs_r,dprojs_i,kpg) &
          !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do idir=1,3
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                dprojs_r(1, ipw, shift_grad+(ilmn-1)*ngrads+igrad+idir) = &
                &     +projs_i(1, ipw, shift+ilmn)*kpg(ipw,idir)*two_pi
                dprojs_i(1, ipw, shift_grad+(ilmn-1)*ngrads+igrad+idir) = &
                &     -projs_r(1, ipw, shift+ilmn)*kpg(ipw,idir)*two_pi
              end do
            end do
          end do
        end if
      end if

      if(signs==1 .and. (choice==5 .or. choice==51 .or. choice==54 .or. choice==55)) then
        igrad=0; if(choice==54) igrad=3; if(choice==55) igrad=6
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:atom_dprojs,dprojs) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=lmn_beg,nlmn
          do idir=1,3
            do ipw=1,npw
              dprojs(1, ipw, shift_grad+(ilmn-1)*ngrads+igrad+idir) = &
              &     +atom_dprojs(1, ipw, idir, ilmn)
              dprojs(2, ipw, shift_grad+(ilmn-1)*ngrads+igrad+idir) = &
              &     +atom_dprojs(2, ipw, idir, ilmn)
            end do
          end do
        end do
      end if

      if(signs==2 .and. (choice==2)) then
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:projs,dprojs,kpg) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+ilmn) = &
            &      projs(2, ipw, shift+ilmn)*kpg(ipw,idir_pert)*two_pi
            dprojs(2, ipw, shift_grad+ilmn) = &
            &     -projs(1, ipw, shift+ilmn)*kpg(ipw,idir_pert)*two_pi
          end do
        end do
      end if


      if(signs==2 .and. (choice==3)) then
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:atom_dprojs,dprojs) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:atom_dprojs,dprojs) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(idir,ilmn,ipw) COLLAPSE(3) MAP(to:projs,d2projs,kpg) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=lmn_beg,nlmn
          do idir=1,6
            do ipw=1,npw
              d2projs(1, ipw, shift_grad2+(ilmn-1)*ngrads2+idir) = &
              &     -projs(1, ipw, shift+ilmn)*kpg(ipw,idir+3)*two_pi2
              d2projs(2, ipw, shift_grad2+(ilmn-1)*ngrads2+idir) = &
              &     -projs(2, ipw, shift+ilmn)*kpg(ipw,idir+3)*two_pi2
            end do
          end do
        end do
      end if

      if(signs==1 .and. choice==54) then
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(idir1,idir2,ilmn,ipw) COLLAPSE(4) MAP(to:atom_d2projs,d2projs,kpg) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=lmn_beg,nlmn
          do idir1=1,3
            do idir2=1,3
              do ipw=1,npw
                d2projs(1, ipw, shift_grad2+(ilmn-1)*ngrads2+(idir1-1)*3+idir2) = &
                &     -atom_d2projs(2, ipw, idir2, ilmn)*kpg(ipw,idir1)*two_pi
                d2projs(2, ipw, shift_grad2+(ilmn-1)*ngrads2+(idir1-1)*3+idir2) = &
                &     +atom_d2projs(1, ipw, idir2, ilmn)*kpg(ipw,idir1)*two_pi
              end do
            end do
          end do
        end do
      end if

      if(signs==1 .and. choice==55) then
        !!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(idir1,idir2,ilmn,ipw) COLLAPSE(4) MAP(to:atom_d2projs,d2projs,kpg) &
        !!$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=lmn_beg,nlmn
          do idir1=1,6
            do idir2=1,3
              do ipw=1,npw
                d2projs(1, ipw, shift_grad2+(ilmn-1)*ngrads2+(idir1-1)*3+idir2) = &
                &     +atom_d2projs(1, ipw, idir1, ilmn)*kpg(ipw,idir2)
                d2projs(2, ipw, shift_grad2+(ilmn-1)*ngrads2+(idir1-1)*3+idir2) = &
                &     +atom_d2projs(2, ipw, idir1, ilmn)*kpg(ipw,idir2)
!                print*,d2projs(:, ipw, shift_grad2+(ilmn-1)*ngrads2+(idir-1)*3+idir2); flush(6)
              end do
            end do
          end do
        end do
      end if

      iaph3d = iaph3d + 1
      shift_grad2 = shift_grad2 + ngrads2*nlmn_o
      shift_grad  = shift_grad  + ngrads*nlmn_o
      shift       = shift       + nlmn

    end do
  end do

  if(istwf_k==1) then
    !$OMP TARGET UPDATE TO(dprojs) IF(gpu_option==ABI_GPU_OPENMP)
    !$OMP TARGET UPDATE TO(d2projs) IF(ngrads2>0 .and. gpu_option==ABI_GPU_OPENMP)
  end if
  if(signs==1) then
    if(istwf_k==2) then
      !$OMP TARGET UPDATE TO(dprojs_i,dprojs_r) IF(gpu_option==ABI_GPU_OPENMP)
    end if
  end if

  !!$OMP TARGET EXIT DATA MAP(delete:temp,scal) IF(gpu_option==ABI_GPU_OPENMP)
  ABI_FREE(temp)
  ABI_FREE(parity)
  ABI_FREE(scal)
  if (allocated(atom_dprojs)) then
    !!$OMP TARGET EXIT DATA MAP(delete:atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
    ABI_FREE(atom_dprojs)
  end if
  if (allocated(atom_d2projs)) then
    !!$OMP TARGET EXIT DATA MAP(delete:atom_d2projs) IF(gpu_option==ABI_GPU_OPENMP)
    ABI_FREE(atom_d2projs)
  end if
  !!$OMP TARGET EXIT DATA MAP(delete:kpg) IF(gpu_option==ABI_GPU_OPENMP)
  if (nkpg_local>0) then
    ABI_FREE(kpg)
  end if
  !!$OMP TARGET EXIT DATA MAP(delete:ffnl_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl_k)
  !!$OMP TARGET EXIT DATA MAP(delete:ph3d_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d_k)

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxEndRange()
#endif
 end subroutine prep_dprojectors
!!***

end module m_gemm_nonlop_projectors
!!***

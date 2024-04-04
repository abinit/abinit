!!****m* ABINIT/m_gemm_nonlop
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

module m_gemm_nonlop

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_fstrings,    only : itoa, ftoa, sjoin
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernla_gemm, only : opernla_gemm
 use m_opernlb_gemm, only : opernlb_gemm
 use m_opernld_ylm_allwf, only : opernld_ylm_allwf
 use m_opernld_ylm, only : opernld_ylm
 use m_pawcprj, only : pawcprj_type
 use m_geometry, only : strconv
 use m_kg, only : mkkpg

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
 public :: gemm_nonlop
 public :: prep_dprojectors


!!***

!----------------------------------------------------------------------

!!****t* m_gemm_nonlop/gemm_nonlop_type
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

 type(gemm_nonlop_type), save, public, allocatable, target :: gemm_nonlop_kpt(:)
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
 type(gemm_nonlop_gpu_type), save, public, target, allocatable :: gemm_nonlop_kpt_gpu(:)
 !(nkpt)

#endif

#ifdef HAVE_OPENMP_OFFLOAD
 integer, save, public :: current_ikpt_in_gpu=-1
 type(gemm_nonlop_type), pointer, save, public :: gpu_nonlop_current_ikpt
#endif

contains

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/init_gemm_nonlop
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
 subroutine init_gemm_nonlop(nkpt,gpu_option)

  integer,intent(in) :: nkpt, gpu_option
  integer :: ikpt

! *************************************************************************

  ! TODO only allocate the number of kpt treated by this proc
  ABI_MALLOC(gemm_nonlop_kpt, (nkpt))
  do ikpt=1,nkpt
    gemm_nonlop_kpt(ikpt)%nprojs = -1
    gemm_nonlop_kpt(ikpt)%ngrads = -1
    gemm_nonlop_kpt(ikpt)%ngrads2 = -1
    gemm_nonlop_kpt(ikpt)%choice = -1
    gemm_nonlop_kpt(ikpt)%idir = -1
  end do

  if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    ABI_MALLOC(gemm_nonlop_kpt_gpu, (nkpt))
    do ikpt=1,nkpt
      gemm_nonlop_kpt_gpu(ikpt)%npw = -1
      gemm_nonlop_kpt_gpu(ikpt)%nprojs = -1
    end do
    gemm_nonlop_gpu_data % allocated = .false.
#endif
  end if


 end subroutine init_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/destroy_gemm_nonlop
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
 subroutine destroy_gemm_nonlop(nkpt, gpu_option)

  integer,intent(in) :: nkpt, gpu_option
  integer :: ikpt

! *************************************************************************

! TODO add cycling if kpt parallelism
  do ikpt = 1,nkpt
    call free_gemm_nonlop_ikpt(ikpt, gpu_option)
  end do

 if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
   ABI_FREE(gemm_nonlop_kpt_gpu)
#endif
 end if
 ABI_FREE(gemm_nonlop_kpt)

 end subroutine destroy_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/free_gemm_nonlop_ikpt
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
 subroutine free_gemm_nonlop_ikpt(ikpt, gpu_option)

  integer,intent(in) :: ikpt, gpu_option

! *************************************************************************

 if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
   if(gemm_nonlop_kpt_gpu(ikpt)%nprojs /= -1) then
     ! deallocate arrays projs, projs_r and projs_i
     if (allocated(gemm_nonlop_kpt(ikpt)%projs)) then
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%projs_r)) then
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%projs_i)) then
       call dealloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_i)
     end if
     gemm_nonlop_kpt_gpu(ikpt)%nprojs = -1
   end if
#endif
 end if

 if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   if(ikpt==current_ikpt_in_gpu) call free_ompgpu_current_ikpt()
#endif
 end if

 if(gemm_nonlop_kpt(ikpt)%nprojs /= -1) then
   if (allocated(gemm_nonlop_kpt(ikpt)%projs)) then
     ABI_FREE(gemm_nonlop_kpt(ikpt)%projs)
   end if
   if (allocated(gemm_nonlop_kpt(ikpt)%projs_r)) then
     ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_r)
   end if
   if (allocated(gemm_nonlop_kpt(ikpt)%projs_i)) then
   ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_i)
   end if
   gemm_nonlop_kpt(ikpt)%nprojs = -1
   if(gemm_nonlop_kpt(ikpt)%ngrads /= -1) then
     if (allocated(gemm_nonlop_kpt(ikpt)%dprojs)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%dprojs_r)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%dprojs_i)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs_i)
     end if
     gemm_nonlop_kpt(ikpt)%ngrads = -1
   end if
   if(gemm_nonlop_kpt(ikpt)%ngrads2 /= -1) then
     if (allocated(gemm_nonlop_kpt(ikpt)%d2projs)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%d2projs)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%d2projs_r)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%d2projs_r)
     end if
     if (allocated(gemm_nonlop_kpt(ikpt)%d2projs_i)) then
       ABI_FREE(gemm_nonlop_kpt(ikpt)%d2projs_i)
     end if
     gemm_nonlop_kpt(ikpt)%ngrads2 = -1
   end if
 end if

 end subroutine free_gemm_nonlop_ikpt
!!***

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt()

#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  call free_ompgpu_current_ikpt_projs()
  call free_ompgpu_current_ikpt_dprojs()

  current_ikpt_in_gpu=-1
  nullify(gpu_nonlop_current_ikpt)
#endif
 end subroutine free_ompgpu_current_ikpt

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt_projs()

#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  if(xomp_target_is_present(c_loc(gpu_nonlop_current_ikpt%projs))) then
    !$OMP TARGET EXIT DATA MAP(delete:gpu_nonlop_current_ikpt%projs)
  end if
  if(xomp_target_is_present(c_loc(gpu_nonlop_current_ikpt%projs_r))) then
    !$OMP TARGET EXIT DATA MAP(delete:gpu_nonlop_current_ikpt%projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:gpu_nonlop_current_ikpt%projs_r)
  end if

#endif
 end subroutine free_ompgpu_current_ikpt_projs

!----------------------------------------------------------------------


 subroutine free_ompgpu_current_ikpt_dprojs()

#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  if(xomp_target_is_present(c_loc(gpu_nonlop_current_ikpt%dprojs))) then
    !$OMP TARGET EXIT DATA MAP(delete:gpu_nonlop_current_ikpt%dprojs)
  end if
  if(xomp_target_is_present(c_loc(gpu_nonlop_current_ikpt%dprojs_i))) then
    !$OMP TARGET EXIT DATA MAP(delete:gpu_nonlop_current_ikpt%dprojs_i)
    !$OMP TARGET EXIT DATA MAP(delete:gpu_nonlop_current_ikpt%dprojs_r)
  end if

#endif
 end subroutine free_ompgpu_current_ikpt_dprojs

!----------------------------------------------------------------------


!!****f* m_gemm_nonlop/make_gemm_nonlop
!! NAME
!! make_gemm_nonlop
!!
!! FUNCTION
!! Build the gemm_nonlop array
!!
!! INPUTS
!!
!! SOURCE
 subroutine make_gemm_nonlop(ikpt,signs,choice,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol, &
&                            ffnl_k,ph3d_k,kpt_k,kg_k,kpg_k, &
&                            idir_pert,gpu_option) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  integer, intent(in) :: signs,choice
  integer, intent(in), optional :: idir_pert,gpu_option
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)

  integer :: nprojs,ndprojs,ndgxdt,nd2gxdt

  integer :: rank, nprocs, ierr, itypat, i, ipw
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

  ABI_CHECK(allocated(gemm_nonlop_kpt),"init_gemm_nonlop wasn't called prior make_gemm_nonlop !")

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
  ndgxdt=0;nd2gxdt=0
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





  compute_dprojs = .false.
  if(nprojs/=gemm_nonlop_kpt(ikpt)%nprojs) then

    call free_gemm_nonlop_ikpt(ikpt,gpu_option_)
    if(gpu_option_ == ABI_GPU_OPENMP) then
      call free_ompgpu_current_ikpt_projs()
      call free_ompgpu_current_ikpt_dprojs()
    end if

    if(istwf_k <= 1) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs, (2, npw, nprojs))
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1,1,1))
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1,1,1))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%projs) IF(gpu_option_==ABI_GPU_OPENMP)
    else
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1, npw, nprojs))
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1, npw, nprojs))
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs, (1, 1, 1))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%projs_r) IF(gpu_option_==ABI_GPU_OPENMP)
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%projs_i) IF(gpu_option_==ABI_GPU_OPENMP)
    end if

    call prep_projectors(npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k,&
    &                    ucvol,ffnl_k,ph3d_k,&
    &                    nprojs,choice,gpu_option_,&
    &                    gemm_nonlop_kpt(ikpt)%projs,&
    &                    gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i)
  else
#ifdef HAVE_OPENMP_OFFLOAD
    if(gpu_option_==ABI_GPU_OPENMP .and. ikpt/=current_ikpt_in_gpu) then
      call free_ompgpu_current_ikpt_projs()
      call free_ompgpu_current_ikpt_dprojs()
      if(istwf_k <= 1) then
        !$OMP TARGET ENTER DATA MAP(to:gemm_nonlop_kpt(ikpt)%projs)
      else
        !$OMP TARGET ENTER DATA MAP(to:gemm_nonlop_kpt(ikpt)%projs_r)
        !$OMP TARGET ENTER DATA MAP(to:gemm_nonlop_kpt(ikpt)%projs_i)
      end if
    end if
#endif
  end if


  if(ndgxdt>0) then
    if(nprojs/=gemm_nonlop_kpt(ikpt)%nprojs .or. ndgxdt /= gemm_nonlop_kpt(ikpt)%ngrads) then
      if(gpu_option_ == ABI_GPU_OPENMP) call free_ompgpu_current_ikpt_dprojs()
        if(allocated(gemm_nonlop_kpt(ikpt)%dprojs)) ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs)
        if(allocated(gemm_nonlop_kpt(ikpt)%dprojs_r)) ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs_r)
        if(allocated(gemm_nonlop_kpt(ikpt)%dprojs_i)) ABI_FREE(gemm_nonlop_kpt(ikpt)%dprojs_i)
        if(allocated(gemm_nonlop_kpt(ikpt)%d2projs)) ABI_FREE(gemm_nonlop_kpt(ikpt)%d2projs)
        if(allocated(gemm_nonlop_kpt(ikpt)%d2projs_r)) ABI_FREE(gemm_nonlop_kpt(ikpt)%d2projs_r)
        if(allocated(gemm_nonlop_kpt(ikpt)%d2projs_i)) ABI_FREE(gemm_nonlop_kpt(ikpt)%d2projs_i)
        gemm_nonlop_kpt(ikpt)%ngrads = -1
        gemm_nonlop_kpt(ikpt)%ngrads2 = -1
      if(istwf_k <= 1) then
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs, (2, npw, nprojs*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_r, (1, 1, 1))
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_i, (1, 1, 1))
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs) IF(gpu_option_==ABI_GPU_OPENMP)
        if(nd2gxdt>0) then
          ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs, (2, npw, nprojs*nd2gxdt))
          ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs_r, (1, 1, 1))
          ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs_i, (1, 1, 1))
          !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs) IF(gpu_option_==ABI_GPU_OPENMP)
        end if
      else
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_r, (1, npw, nprojs*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_i, (1, npw, nprojs*ndgxdt))
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs, (1, 1, 1))
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs_r) IF(gpu_option_==ABI_GPU_OPENMP)
        !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs_i) IF(gpu_option_==ABI_GPU_OPENMP)
        if(nd2gxdt>0) then
          ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs_r, (1, npw, nprojs*nd2gxdt))
          ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs_i, (1, npw, nprojs*nd2gxdt))
          ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs, (1, 1, 1))
          !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs_r) IF(gpu_option_==ABI_GPU_OPENMP)
          !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs_i) IF(gpu_option_==ABI_GPU_OPENMP)
        end if
      end if
      if(nd2gxdt==0) then
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs, (1, 1, 1))
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs_r, (1, 1, 1))
        ABI_MALLOC(gemm_nonlop_kpt(ikpt)%d2projs_i, (1, 1, 1))
        !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs) IF(gpu_option_==ABI_GPU_OPENMP)
      end if
      compute_dprojs=.true.
    else
      if(gpu_option_==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
        if((istwf_k==1 .and. .not. xomp_target_is_present(c_loc(gemm_nonlop_kpt(ikpt)%dprojs))) .or. &
        &  (istwf_k==2 .and. .not. xomp_target_is_present(c_loc(gemm_nonlop_kpt(ikpt)%dprojs_r)))) then
          if (ikpt/=current_ikpt_in_gpu) then
            call free_ompgpu_current_ikpt_dprojs()
          end if
          if(istwf_k <= 1) then
            !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs)
            !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs) if(nd2gxdt>0)
          else
            !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs_r)
            !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs_i)
            !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs_r) if(nd2gxdt>0)
            !!$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%d2projs_i) if(nd2gxdt>0)
          end if
          compute_dprojs = .true.
        end if
#endif
      end if
    end if

    if(choice/=gemm_nonlop_kpt(ikpt)%choice .or. idir_pert_/=gemm_nonlop_kpt(ikpt)%idir &
    &  .or. compute_dprojs) then
      call prep_dprojectors(npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k,&
      &                    ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k,&
      &                    nprojs,ndgxdt,nd2gxdt,choice,signs,idir_pert_,gpu_option_,&
      &                    gemm_nonlop_kpt(ikpt)%projs,&
      &                    gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i,&
      &                    gemm_nonlop_kpt(ikpt)%dprojs,&
      &                    gemm_nonlop_kpt(ikpt)%dprojs_r,gemm_nonlop_kpt(ikpt)%dprojs_i,&
      &                    gemm_nonlop_kpt(ikpt)%d2projs,&
      &                    gemm_nonlop_kpt(ikpt)%d2projs_r,gemm_nonlop_kpt(ikpt)%d2projs_i)
    end if
  end if

  if (nprojs>0) gemm_nonlop_kpt(ikpt)%nprojs = nprojs
  if (ndgxdt>0) gemm_nonlop_kpt(ikpt)%ngrads = ndgxdt
  if (nd2gxdt>0) gemm_nonlop_kpt(ikpt)%ngrads2 = nd2gxdt
  if (ndgxdt>0) gemm_nonlop_kpt(ikpt)%choice = choice
  if (ndgxdt>0) gemm_nonlop_kpt(ikpt)%idir = idir_pert_
  !!!!! GPU stuff

  if(gpu_option_ == ABI_GPU_LEGACY .or. gpu_option_ == ABI_GPU_KOKKOS) then
#ifdef HAVE_GPU_CUDA
    gemm_nonlop_kpt_gpu(ikpt)%npw    = npw
    gemm_nonlop_kpt_gpu(ikpt)%nprojs = nprojs

#ifdef DEBUG_VERBOSE_GPU
    if(xmpi_comm_rank(xmpi_world) == 0) then
      call check_gpu_mem("make_gemm_nonlop begin")
      call wrtout(std_out,sjoin(" npw    .......", itoa(npw)),    'COLL')
      call wrtout(std_out,sjoin(" nprojs .......", itoa(nprojs)), 'COLL')
    end if
#endif

    if(istwf_k <= 1) then
      call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs, INT(2,c_size_t)*npw*nprojs*dp)
      ! TODO : gradients
    else
      call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_r, INT(1, c_size_t)*npw*nprojs*dp)
      call alloc_on_gpu(gemm_nonlop_kpt_gpu(ikpt)%projs_i, INT(1, c_size_t)*npw*nprojs*dp)
      ! TODO : gradients
    end if

#ifdef DEBUG_VERBOSE_GPU
    if(xmpi_comm_rank(xmpi_world) == 0) then
      call check_gpu_mem("make_gemm_nonlop end  ")
    end if
#endif

    ! upload data to gpu memory
    if(istwf_k <= 1) then
      call copy_on_gpu(gemm_nonlop_kpt(ikpt)%projs, gemm_nonlop_kpt_gpu(ikpt)%projs, INT(2, c_size_t)*npw*nprojs*dp)
      ! TODO : gradients
    else
      call copy_on_gpu(gemm_nonlop_kpt(ikpt)%projs_r, gemm_nonlop_kpt_gpu(ikpt)%projs_r, &
        &                    INT(1, c_size_t)*npw*nprojs*dp)
      call copy_on_gpu(gemm_nonlop_kpt(ikpt)%projs_i, gemm_nonlop_kpt_gpu(ikpt)%projs_i, &
        &                    INT(1, c_size_t)*npw*nprojs*dp)
    end if

#endif
  else if(gpu_option_ == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    current_ikpt_in_gpu = ikpt
    gpu_nonlop_current_ikpt => gemm_nonlop_kpt(ikpt)
#endif
  end if
  gemm_nonlop_choice=choice
#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxEndRange()
#endif

 end subroutine make_gemm_nonlop
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/prep_projectors
!! NAME
!! prep_projectors
!!
!! FUNCTION
!! Prepare projectors array and derivatives for given choice
!!
!! INPUTS
!!
!! SOURCE
 subroutine prep_projectors(npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k,&
 &                          ucvol,ffnl_k,ph3d_k,&
 &                          nprojs,choice,gpu_option,&
 &                          projs,projs_r,projs_i)

  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  integer, intent(in) :: nprojs,choice,gpu_option
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp),intent(inout), target :: projs(2,npw,nprojs)
  real(dp),intent(inout), target :: projs_r(1,npw,nprojs)
  real(dp),intent(inout), target :: projs_i(1,npw,nprojs)

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
          ABI_BUG("toto")
          !$OMP TARGET DATA USE_DEVICE_PTR(projs_r,projs_i,atom_projs)
          call copy_gpu_to_gpu(c_loc(projs_r(:,:,shift+1)), c_loc(atom_projs(:, :, lmn_beg)), int(1,c_size_t)*npw*(nlmn-lmn_beg)*dp)
          call copy_gpu_to_gpu(c_loc(projs_i(:,:,shift+1)), c_loc(atom_projs(:, :, lmn_beg)), int(1,c_size_t)*npw*(nlmn-lmn_beg)*dp)
          !$OMP END TARGET DATA
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

!!****f* m_gemm_nonlop/prep_dprojectors
!! NAME
!! prep_projectors
!!
!! FUNCTION
!! Prepare projectors array and derivatives for given choice
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
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)
  real(dp),intent(inout), target :: projs(2,npw,nprojs)  ,dprojs(2,npw,nprojs*ngrads)  ,d2projs(2,npw,nprojs*ngrads2)
  real(dp),intent(inout), target :: projs_r(1,npw,nprojs),dprojs_r(1,npw,nprojs*ngrads),d2projs_r(1,npw,nprojs*ngrads2)
  real(dp),intent(inout), target :: projs_i(1,npw,nprojs),dprojs_i(1,npw,nprojs*ngrads),d2projs_i(1,npw,nprojs*ngrads2)

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

  if(gpu_option==ABI_GPU_OPENMP) then

    if(istwf_k <= 1) then
      call gpu_set_to_zero(dprojs,int(2,c_size_t)*npw*nprojs*ngrads)
      if(ngrads2 > 0) then
        call gpu_set_to_zero(d2projs,int(2,c_size_t)*npw*nprojs*ngrads2)
      end if
    else
      call gpu_set_to_zero(dprojs_r,int(npw,c_size_t)*nprojs*ngrads)
      call gpu_set_to_zero(dprojs_i,int(npw,c_size_t)*nprojs*ngrads)
      if(ngrads2 > 0) then
        call gpu_set_to_zero(d2projs_r,int(npw,c_size_t)*nprojs*ngrads2)
        call gpu_set_to_zero(d2projs_i,int(npw,c_size_t)*nprojs*ngrads2)
      end if
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
      if(ngrads2 > 0) then
        d2projs_r(:,:,:) = zero
        d2projs_i(:,:,:) = zero
      end if
    end if

  end if

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
    !$OMP TARGET ENTER DATA MAP(alloc:atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
    if(ngrads2>0) then
      ABI_MALLOC(atom_d2projs, (2, npw, nd2projs, lmnmax))
      !$OMP TARGET ENTER DATA MAP(alloc:atom_d2projs) IF(gpu_option==ABI_GPU_OPENMP)
    end if
  end if

  ABI_MALLOC(temp, (npw, lmnmax))
  ABI_MALLOC(parity, (lmnmax))
  ABI_MALLOC(scal, (lmnmax))
  !$OMP TARGET ENTER DATA MAP(alloc:temp,scal) IF(gpu_option==ABI_GPU_OPENMP)

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

  !$OMP TARGET ENTER DATA MAP(to:kpg) IF(gpu_option==ABI_GPU_OPENMP)
  map_ph3d_k=.false.; map_ffnl_k=.false.
  if(.not. xomp_target_is_present(c_loc(ffnl_k))) map_ffnl_k = .true.
  if(.not. xomp_target_is_present(c_loc(ph3d_k))) map_ph3d_k = .true.
  !$OMP TARGET ENTER DATA MAP(to:ffnl_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl_k)
  !$OMP TARGET ENTER DATA MAP(to:ph3d_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d_k)

  shift = 0 ; shift_grad = 0; shift_grad2 = 0
  lmn_beg = 1

  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)
    nlmn_o = nlmn

    do ia = 1, nattyp(itypat)

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
      if (signs==1 .and. (choice==3 .or. choice==23 .or. choice==54 .or. choice==55)) then
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_dprojs,ffnl_k) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_dprojs(1,ipw, 1:ndprojs, 1:nlmn_o) = wt * ffnl_k(ipw, 2:ndprojs+1, 1:nlmn_o, itypat)
        end do
      end if
      if (signs==2 .and. (choice==3 .or. choice==5 .or. choice==51)) then
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_dprojs,ffnl_k) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_dprojs(1,ipw, 1, 1:nlmn_o) = wt * ffnl_k(ipw, 1+ffnl_dir, 1:nlmn_o, itypat)
        end do
      end if
      if(signs==1 .and. choice==54) then
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl_k) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ipw=1, npw
          atom_d2projs(1,ipw, 1:nd2projs, 1:nlmn_o) = wt * ffnl_k(ipw, 2:nd2projs+1, 1:nlmn_o, itypat)
        end do
      else if(signs==1 .and. choice==55) then
        !$OMP TARGET PARALLEL DO PRIVATE(ipw) MAP(to:atom_d2projs,ffnl_k) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
          !$OMP TARGET UPDATE TO(scal)
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
          !$OMP& PRIVATE(idir,ipw,ilmn) MAP(to:atom_dprojs,scal) &
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
            !$OMP TARGET UPDATE TO(scal)
            !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
            !$OMP& PRIVATE(idir,ipw,ilmn) MAP(to:atom_d2projs,scal) &
            !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE MAP(to:atom_dprojs,parity) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
        do ilmn=1,nlmn_o
          if(.not. parity(ilmn)) then
            !$OMP PARALLEL DO PRIVATE(idir,ipw,tmp) COLLAPSE(2)
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
          !$OMP TARGET TEAMS DISTRIBUTE MAP(to:atom_d2projs,parity) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
          do ilmn=1,nlmn_o
            if(.not. parity(ilmn)) then
              !$OMP PARALLEL DO PRIVATE(idir,ipw,tmp) COLLAPSE(2)
              do idir=1,nd2projs
                do ipw=1,npw
                  tmp = atom_d2projs(2,ipw,idir,ilmn)
                  atom_d2projs(2,ipw,idir,ilmn) =  atom_d2projs(1,ipw,idir,ilmn)
                  atom_d2projs(1,ipw,idir,ilmn) = -tmp
                end do
              end do
            else
              !$OMP PARALLEL DO PRIVATE(idir,ipw,tmp) COLLAPSE(2)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,idir,ipw,tmp) COLLAPSE(3) MAP(to:atom_dprojs,ph3d_k) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,idir,ipw,tmp) COLLAPSE(3) MAP(to:atom_d2projs,ph3d_k) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw,idir,idir1,idir2) COLLAPSE(3) MAP(to:atom_dprojs,dprojs,kpg) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
          !$OMP TARGET UPDATE FROM(atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
          do idir=1,6
            idir1=alpha(idir);idir2=beta(idir)
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                dprojs_r(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     -half*(atom_dprojs(1, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(1, ipw, idir2, ilmn)*kpg(ipw,idir1))

                dprojs_i(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
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
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw,idir) COLLAPSE(3) MAP(to:projs,dprojs,kpg) &
          !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
          !$OMP TARGET UPDATE FROM(projs_r,projs_i) IF(gpu_option==ABI_GPU_OPENMP)
          do idir=1,3
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                dprojs_r(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     +projs_i(1, ipw, shift+ilmn)*kpg(ipw,idir)*two_pi
                dprojs_i(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     -projs_r(1, ipw, shift+ilmn)*kpg(ipw,idir)*two_pi
              end do
            end do
          end do
        end if
      end if

      if(signs==1 .and. (choice==5 .or. choice==51 .or. choice==54 .or. choice==55)) then
        igrad=0; if(choice==54) igrad=3; if(choice==55) igrad=6
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:atom_dprojs,dprojs) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:projs,dprojs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:atom_dprojs,dprojs) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(ilmn,ipw) COLLAPSE(2) MAP(to:atom_dprojs,dprojs) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(idir,ilmn,ipw) COLLAPSE(2) MAP(to:projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(idir1,idir2,ilmn,ipw) COLLAPSE(2) MAP(to:atom_d2projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(idir1,idir2,ilmn,ipw) COLLAPSE(2) MAP(to:atom_d2projs,d2projs,kpg) &
        !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
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

  !if(.not. (signs==2 .and. (choice==5 .or. choice==51 .or. choice==3))) then
  if(signs==1) then
    if(istwf_k==2) then
      !$OMP TARGET UPDATE TO(dprojs_i,dprojs_r) IF(gpu_option==ABI_GPU_OPENMP)
    end if
  end if

  !$OMP TARGET EXIT DATA MAP(delete:temp,scal) IF(gpu_option==ABI_GPU_OPENMP)
  ABI_FREE(temp)
  ABI_FREE(parity)
  ABI_FREE(scal)
  if (allocated(atom_dprojs)) then
    !$OMP TARGET EXIT DATA MAP(delete:atom_dprojs) IF(gpu_option==ABI_GPU_OPENMP)
    ABI_FREE(atom_dprojs)
  end if
  if (allocated(atom_d2projs)) then
    !$OMP TARGET EXIT DATA MAP(delete:atom_d2projs) IF(gpu_option==ABI_GPU_OPENMP)
    ABI_FREE(atom_d2projs)
  end if
  !$OMP TARGET EXIT DATA MAP(delete:kpg) IF(gpu_option==ABI_GPU_OPENMP)
  if (nkpg_local>0) then
    ABI_FREE(kpg)
  end if
  !$OMP TARGET EXIT DATA MAP(delete:ffnl_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ffnl_k)
  !$OMP TARGET EXIT DATA MAP(delete:ph3d_k) IF(gpu_option==ABI_GPU_OPENMP .and. map_ph3d_k)

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxEndRange()
#endif
 end subroutine prep_dprojectors
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/gemm_nonlop
!! NAME
!! gemm_nonlop
!!
!! FUNCTION
!! Replacement of nonlop. same prototype as nonlop although not all options are implemented.
!!
!! INPUTS
!! [gpu_option] = GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!
!! SOURCE
 subroutine gemm_nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,atom_proj_shift,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm,atom_proj_shift
  integer,optional,intent(in) :: gpu_option
  real(dp),intent(in) :: lambda(ndat),ucvol
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3)
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin*useylm)
  real(dp),intent(in) :: kpgout(npwout,nkpgout*useylm),kptin(3),kptout(3)
  real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom),ph1d(2,*)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
  real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out),target :: svectout(:,:)
  real(dp),intent(inout),target :: vectout(:,:)
  real(dp),intent(inout),optional, ABI_CONTIGUOUS target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  complex(dpc), parameter :: cminusone  = (-1._dp,0._dp)
  integer :: ii, ia, idat, igrad, nprojs, ngrads, ngrads2, shift, iatom, nlmn, ierr, ibeg, iend, ikpt
  integer :: cplex, cplex_enl, cplex_fac, proj_shift, grad_shift
  integer :: projs_beg,projs_end,dprojs_beg,dprojs_end,d2projs_beg,d2projs_end
  integer :: nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(9), cplex_d2gxdt(18)
  logical :: local_vectproj
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: sij_typ(:)
  real(dp), ABI_CONTIGUOUS pointer :: projections(:,:,:)
  real(dp), allocatable :: s_projections(:,:,:), vnl_projections(:,:,:)
  real(dp), allocatable :: dprojections(:,:,:), temp_realvec(:)
  real(dp), allocatable, target :: s_dprojections(:,:,:), vnl_dprojections(:,:,:)
  real(dp), allocatable, target :: d2projections(:,:,:)
  integer :: nprojs_my_blk
  integer :: rank, nprocs
  logical :: is_last
  real(dp), allocatable :: projs_local(:,:,:)
  real(dp), allocatable :: projs_recv(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_(:,:,:),dprojs_(:,:,:),d2projs_(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_r_(:,:,:),projs_i_(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: dprojs_r_(:,:,:),dprojs_i_(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: d2projs_r_(:,:,:),d2projs_i_(:,:,:)
  integer :: ngrads_tmp,ngrads2_tmp
  real(dp), allocatable :: enlk(:),fnlk(:,:),ddkk(:,:),strnlk(:,:),gmet2(:,:)
  real(dp), allocatable :: work1(:),work2(:),work3(:,:),work4(:,:),work5(:,:,:),work6(:,:,:),work7(:,:,:)
  integer :: idbeg,idend,dshift,id2beg,id2end,d2shift,enlout_shift,ipw,i
  real(dp) :: work(6)
  integer :: mu0,ic,nu,mu,jc
  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,kpgin,kpgout/))
  ABI_UNUSED((/ph1d(1,1),ph3din,ph3dout/))
  ABI_UNUSED((/phkxredin,phkxredout,ucvol/))
  ABI_UNUSED((/mgfft,mpsang,mpssoang/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,gpu_option/))

  ! Check supported options
  if (.not.gemm_nonlop_use_gemm) then
    ABI_BUG('computation not prepared for gemm_nonlop use!')
  end if
  if ( (choice>3.and.choice/=7.and.choice/=5.and.choice/=51.and.signs==2) .or. &
&      (choice>3.and.choice/=7.and.choice/=23.and.choice/=4.and.choice/=54.and.choice/=55.and.signs==1) .or. &
&      (useylm/=1) ) then
    ABI_BUG('gemm_nonlop option not supported!')
  end if
  if (signs==1) then
    nnlout_test=0
    if (choice==1) nnlout_test=1
    if (choice==2) nnlout_test=3*natom
    if (choice==3) nnlout_test=6
    if (choice==23) nnlout_test=6+3*natom
    if (nnlout<nnlout_test) then
      ABI_BUG('wrong nnlout size!')
    end if
  end if

  ikpt=gemm_nonlop_ikpt_this_proc_being_treated
  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?

  ngrads = gemm_nonlop_kpt(ikpt)%ngrads
  ngrads2 = gemm_nonlop_kpt(ikpt)%ngrads2
  nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_blk

  ! The number of projectors used for computation may vary among
  ! nonlop calls, from computing on all atoms to a select one for
  ! some perturbations.
  ! In such cases, projs arrays must be recomputed
  nprojs=0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  projs_beg=1; projs_end=nprojs;
  dprojs_beg=1; dprojs_end=nprojs*ngrads;
  d2projs_beg=1; d2projs_end=nprojs*ngrads2;
  if((choice==2 .and. signs==2)) then
    projs_beg=atom_proj_shift+1
    projs_end=projs_beg+nprojs-1
    dprojs_beg=atom_proj_shift*ngrads+1
    dprojs_end=dprojs_beg+nprojs*ngrads-1
    d2projs_beg=atom_proj_shift*ngrads2+1
    d2projs_end=dprojs_beg+nprojs*ngrads2-1
  end if
  if(istwf_k == 1) then
    projs_  => gemm_nonlop_kpt(ikpt)%projs(:,:,projs_beg:projs_end)
    if(ngrads>0)  dprojs_ => gemm_nonlop_kpt(ikpt)%dprojs(:,:,dprojs_beg:dprojs_end)
    if(ngrads2>0) d2projs_ => gemm_nonlop_kpt(ikpt)%d2projs(:,:,d2projs_beg:d2projs_end)
  else
    projs_r_  => gemm_nonlop_kpt(ikpt)%projs_r(:,:,projs_beg:projs_end)
    projs_i_  => gemm_nonlop_kpt(ikpt)%projs_i(:,:,projs_beg:projs_end)
    if(ngrads>0)  dprojs_r_ => gemm_nonlop_kpt(ikpt)%dprojs_r(:,:,dprojs_beg:dprojs_end)
    if(ngrads>0)  dprojs_i_ => gemm_nonlop_kpt(ikpt)%dprojs_i(:,:,dprojs_beg:dprojs_end)
    if(ngrads2>0) d2projs_r_ => gemm_nonlop_kpt(ikpt)%d2projs_r(:,:,d2projs_beg:d2projs_end)
    if(ngrads2>0) d2projs_i_ => gemm_nonlop_kpt(ikpt)%d2projs_i(:,:,d2projs_beg:d2projs_end)
  end if

  if(gemm_nonlop_is_distributed) then
    rank = xmpi_comm_rank(gemm_nonlop_block_comm); nprocs = xmpi_comm_size(gemm_nonlop_block_comm)
    is_last = (rank==nprocs-1)
    if(is_last) nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_last_blk
  end if

  if(choice==1 .or. choice==0 .or. choice==7) ngrads=0
  if(signs==1) then
    ABI_CHECK(ngrads>=3.or.choice/=2 ,"Bad allocation in gemm_nonlop (2)!")
    ABI_CHECK(ngrads>=6.or.choice/=3 ,"Bad allocation in gemm_nonlop (3)!")
    ABI_CHECK(ngrads>=9.or.choice/=23,"Bad allocation in gemm_nonlop (23)!")
  else if(signs==2) then
    ABI_CHECK(ngrads==1.or.(choice==1.or.choice==0.or.choice==7) ,"Bad allocation in gemm_nonlop !")
  end if

  ! If vectproj is provided, use it for further calculations, use allocated array otherwise
  local_vectproj=.false.
  if(PRESENT(vectproj)) then
    if(size(vectproj)>1) local_vectproj=.true.
  end if
  if (local_vectproj) projections => vectproj


  if(gemm_nonlop_is_distributed) then
    ABI_MALLOC(projs_recv, (cplex, npwin, MAX(ngrads,1)*gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
    ABI_MALLOC(projs_local, (cplex, npwin, MAX(ngrads,1)*gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
  end if

  if(nprojs == 0) then
    ! TODO check if this is correct
    if(signs == 1) then
      enlout=zero
      return
    end if
    if(signs == 2) then
      vectout = zero
      if(paw_opt>0) svectout = vectin
      return
    end if
  end if

  if(signs == 1) then
    enlout=zero
    ABI_MALLOC(enlk,(ndat))
    enlk=zero
    ABI_MALLOC(fnlk,(3*natom,ndat))
    fnlk=zero
    ABI_MALLOC(ddkk,(6,ndat))
    ddkk=zero
    ABI_MALLOC(strnlk,(6,ndat))
    strnlk=zero
  end if

  ndgxdt = ngrads
  nd2gxdt = ngrads2

  ndgxdtfac = 0; nd2gxdtfac = 0
  if (choice==2) then
    if (signs==2) ndgxdtfac=1
  end if
  if (choice==22) then
    if (signs==2) ndgxdtfac=1
  end if
  if (choice==3) then
    if (signs==2) ndgxdtfac=1
  end if
  if (choice==4) then
    if(signs==1) ndgxdtfac=3
  end if
  if (choice==5) then
    if(signs==2) ndgxdtfac=1
  end if
  if (choice==51) then
    if(signs==2) ndgxdtfac=1
  end if
  if (choice==54) then
    if(signs==1) ndgxdtfac=6
    if(signs==2) ndgxdtfac=1
    if(signs==2) nd2gxdtfac=1
  end if
  if (choice==55) then
    if(signs==1) ndgxdtfac=9
  end if
  if (choice==6) then
    if(signs==1) ndgxdtfac=9
  end if
  ABI_CHECK(ndgxdtfac<=ndgxdt,"BUG: ndgxdtfac>ndgxdt!")
  optder = 0;if (ndgxdtfac>0) optder = 1
  if (nd2gxdtfac>0) optder=2
  cplex_dgxdt(:) = 1 ; cplex_d2gxdt(:) = 1
  ! When istwf_k > 1, gx derivatives can be real or pure imaginary
  ! cplex_dgxdt(i)  = 1 if dgxdt(1,i,:,:)  is real, 2 if it is pure imaginary
  ! cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
  if(ndgxdt > 0) then
   if (choice==5.or.choice==51) cplex_dgxdt(:) = 2
   if (choice==54.and.signs==1) cplex_dgxdt(4:6) = 2
   !if (choice==54.and.signs==2) cplex_dgxdt(:)   = 2
   if (choice==55.and.signs==1) cplex_dgxdt(7:9) = 2
  end if
  if(nd2gxdt > 0) then
    if (choice==54) cplex_d2gxdt(:) = 2
    if (choice==55.and.signs==1) cplex_d2gxdt(1:18)= 2
  end if

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  if(.not. local_vectproj) then
    ABI_MALLOC(projections,(cplex, nprojs,nspinor*ndat))
  end if
  ABI_MALLOC(s_projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs,nspinor*ndat))

  if(.not. local_vectproj) projections = zero
  s_projections = zero
  vnl_projections = zero

  ! Working buffers for storing derivative
  !if (ngrads>0) then
  ngrads_tmp=1; if (ngrads>0) ngrads_tmp=ngrads
    ABI_MALLOC(dprojections,(cplex, ngrads_tmp*nprojs,nspinor*ndat))
    ABI_MALLOC(s_dprojections,(cplex, ngrads_tmp*nprojs,nspinor*ndat))
    ABI_MALLOC(vnl_dprojections,(cplex_fac, ngrads_tmp*nprojs,nspinor*ndat))
    dprojections(:,:,:) = zero
    s_dprojections(:,:,:) = zero
    vnl_dprojections(:,:,:) = zero
  !end if

  ! Working buffers for storing 2nd-derivative
  !if (ngrads2>0) then
  ngrads2_tmp=1; if (ngrads>0) ngrads2_tmp=ngrads2
    ABI_MALLOC(d2projections,(cplex, ngrads2_tmp*nprojs,nspinor*ndat))
    d2projections(:,:,:) = zero
  !end if


  ! determine precisely when temp_realvec needs to be allocated
  ! to factorize allocate (resp. deallocate) at the begining (resp. at the end) of subroutine
  ! to avoid multiple allocate/deallocate that can be costly
  if (cplex /= 2) then
    if ( (cpopt < 2) .or. &
      &  (paw_opt == 3 .or. paw_opt == 4) .or. &
      &  (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)) then
       ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
    end if
  end if

  if(cpopt >= 2) then
    ! retrieve from cprjin
    if(.not. local_vectproj .and. cpopt/=3) then
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          projections(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
          shift = shift + nlmn
        end do
      end do
    end if
    if(cpopt==4.and.allocated(dprojections)) then
      ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (1)")
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn  = cprjin(iatom, idat)%nlmn
          do ilmn=1,nlmn
            do igrad=1,ngrads
              dprojections(1:cplex, shift + igrad, idat) = &
                cprjin(iatom, idat)%dcp(1:cplex,igrad,ilmn)
            end do
            shift = shift + ngrads
          end do
        end do
      end do
    end if
  end if ! cpopt

  if(cpopt<=1.or.(cpopt<=3.and.(choice==2.or.choice==3.or.choice==5.or.choice==51.or.choice==23.or.choice==54.or.choice==55.or.choice==4))) then

    call opernla_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,d2projections,dprojections,projections,&
    &       idir,istwf_k,mpi_enreg,nd2gxdt,ngrads,&
    &       npwin,nspinor,signs,ndat,rank,&
    &       cpopt,nprocs,&
    &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
    &       vectin,projs_,dprojs_,d2projs_,&
    &       projs_r_,projs_i_,&
    &       dprojs_r_,dprojs_i_,&
    &       d2projs_r_,d2projs_i_,&
    &       temp_realvec,&
    &       projs_local,projs_recv,&
    &       gpu_option,gemm_nonlop_is_distributed)

    if(cpopt >= 0) then
      ! store in cprjin
      if(.not. local_vectproj .and. cpopt/=3) then
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections(1:cplex, shift+1:shift+nlmn, idat)
            shift = shift + nlmn
          end do
        end do
      end if
      if(cpopt==1 .or. cpopt==3) then
        ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (2)")
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            do ilmn=1,nlmn
              do igrad=1,ngrads
                cprjin(iatom, idat)%dcp(1:cplex,igrad,ilmn) = &
                &                   dprojections(1:cplex, shift + igrad, idat)
              end do
              shift = shift + ngrads
            end do
          end do
        end do
      end if
    end if ! cpopt >= 0
  end if ! cpopt >= 2

  if(choice > 0) then

    if(choice /= 7) then
      ! opernlc
      iatm = 0
      shift = 0; dshift = 0; d2shift = 0
      ABI_MALLOC(sij_typ,(((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2))
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)
        if (paw_opt>=2) then
          if (cplex_enl==1) then
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij(ilmn,itypat)
            end do
          else
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij(2*ilmn-1,itypat)
            end do
          end if
        end if

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        idbeg = dshift+1
        idend = dshift+nattyp(itypat)*nlmn*ngrads_tmp

        id2beg = d2shift+1
        id2end = d2shift+nattyp(itypat)*nlmn*ngrads2_tmp

        do idat = 1,ndat
          call opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
&         dprojections(:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_dprojections(:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
&         s_dprojections(:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
&         d2projections(:, id2beg:id2end, 1+nspinor*(idat-1):nspinor*idat),&
&         d2gxdt_dum_out,d2gxdt_dum_out2,&
&         dimenl1,dimenl2,dimekbq,enl,&
&         projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         s_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         iatm,indlmn(:,:,itypat),itypat,lambda(idat),mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&         nattyp(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)
        end do

        shift = shift + nattyp(itypat)*nlmn
        dshift = dshift + nattyp(itypat)*nlmn*ngrads_tmp
        d2shift = d2shift + nattyp(itypat)*nlmn*ngrads2_tmp
        iatm = iatm+nattyp(itypat)
      end do
      ABI_FREE(sij_typ)
    else
      s_projections = projections
    end if ! choice /= 7

    ! opernlb
    if(signs==2) then
      call opernlb_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
      &       d2gxdt_dum_out,d2gxdt_dum_out,&
      &       vnl_dprojections,s_dprojections,&
      &       vnl_projections,s_projections,&
      &       idir,istwf_k,mpi_enreg,nd2gxdt,nd2gxdtfac,ndgxdt,ndgxdtfac,&
      &       npwout,nspinor,signs,ndat,rank,&
      &       cpopt,nprocs,paw_opt,&
      &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
      &       vectin,vectout,svectout,&
      &       projs_,&
      &       dprojs_,&
      &       projs_r_,projs_i_,&
      &       dprojs_r_,dprojs_i_,&
      &       temp_realvec,temp_realvec,&
      &       projs_local,projs_recv,&
      &       gpu_option,gemm_nonlop_is_distributed)
    end if

    ! opernld
    if(signs==1) then
      if(choice==2 .or. choice==3 .or. choice==23) then
        call opernld_ylm_allwf(choice,cplex,cplex_fac,&
        &       dprojections,vnl_dprojections,s_dprojections,d2projections,&
        &       enlk,enlout,projections,vnl_projections,s_projections,&
        &       ndat,nd2gxdt,ndgxdt,&
        &       ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
        &       nattyp,gpu_option)
      else
        shift=0; dshift=0; d2shift = 0; iatm=1
        do itypat=1, ntypat
          nlmn=count(indlmn(3,:,itypat)>0)

          ibeg = shift+1
          iend = shift+nattyp(itypat)*nlmn

          idbeg = dshift+1
          idend = dshift+nattyp(itypat)*nlmn*ngrads_tmp

          id2beg = d2shift+1
          id2end = d2shift+nattyp(itypat)*nlmn*ngrads2_tmp

          do idat=1,ndat
            call opernld_ylm             (choice,cplex,cplex_fac,ddkk(:,idat),&
            &       dprojections    (:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
            &       vnl_dprojections(:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
            &       s_dprojections  (:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
            &       d2projections (:, id2beg:id2end, 1+nspinor*(idat-1):nspinor*idat),&
            &       enlk(idat),enlout(nnlout*(idat-1)+1:nnlout*idat),fnlk(:,idat),&
            &       projections    (:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
            &       vnl_projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
            &       s_projections  (:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
            &       iatm,natom,1,nd2gxdt,ndgxdt,ndgxdtfac,&
            &       nattyp(itypat),nlmn,nnlout,nspinor,paw_opt,strnlk(:,idat))
          end do

          shift = shift + nattyp(itypat)*nlmn
          dshift = dshift + nattyp(itypat)*nlmn*ngrads_tmp
          d2shift = d2shift + nattyp(itypat)*nlmn*ngrads2_tmp
          iatm = iatm+nattyp(itypat)
        end do
      end if

      ! Reduction in case of parallelism
      if (mpi_enreg%paral_spinor==1) then
        if (size(enlout)>0) then
          call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
        end if
        if (choice==3.or.choice==23) then
          call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
        end if
        if (choice==55) then
          call xmpi_sum(ddkk,mpi_enreg%comm_spinor,ierr)
        end if
      end if

      !Need sometimes gmet
      if ((signs==1.and.paw_opt<=3).and. &
          & (choice==5 .or.choice==51.or.choice==52.or.choice==53.or.&
          & choice==54.or.choice==55)) then
        ABI_MALLOC(gmet2,(3,3))
        gmet2 = MATMUL(TRANSPOSE(gprimd),gprimd)
      end if

      !Coordinate transformations

      ! Derivatives wrt strain
      !  - Convert from reduced to cartesian coordinates
      !  - Substract volume contribution
      if ((choice==3.or.choice==23).and.paw_opt<=3) then
        do idat=1,ndat
          enlout_shift=(idat-1)*nnlout
          call strconv(enlout(enlout_shift+1:enlout_shift+6),gprimd,work)
          enlout(enlout_shift+1:enlout_shift+3)=(work(1:3)-enlk(idat))
          enlout(enlout_shift+4:enlout_shift+6)= work(4:6)
        end do
      end if

      !2nd derivative wrt to k wave vector and atomic position (effective charges):
      ! - convert from cartesian to reduced coordinates
      if (choice==54.and.signs==1.and.paw_opt<=3) then
        ABI_MALLOC(work1,(3))
        ABI_MALLOC(work2,(3))
        do idat=1,ndat
          mu0=0 ! Shift to be applied in enlout array
          enlout_shift=(idat-1)*nnlout
          do mu=1,3*natom
        !   First, real part
            work1(1)=enlout(enlout_shift+mu0+1);work1(2)=enlout(enlout_shift+mu0+3);work1(3)=enlout(enlout_shift+mu0+5)
            work2(:)=gmet2(:,1)*work1(1)+gmet2(:,2)*work1(2)+gmet2(:,3)*work1(3)
            enlout(enlout_shift+mu0+1)=work2(1);enlout(enlout_shift+mu0+3)=work2(2);enlout(enlout_shift+mu0+5)=work2(3)
        !   Then imaginary part
            work1(1)=enlout(enlout_shift+mu0+2);work1(2)=enlout(enlout_shift+mu0+4);work1(3)=enlout(enlout_shift+mu0+6)
            work2(:)=gmet2(:,1)*work1(1)+gmet2(:,2)*work1(2)+gmet2(:,3)*work1(3)
            enlout(enlout_shift+mu0+2)=work2(1);enlout(enlout_shift+mu0+4)=work2(2);enlout(enlout_shift+mu0+6)=work2(3)
            mu0=mu0+6
          end do
        end do !idat
        ABI_FREE(work1)
        ABI_FREE(work2)
      end if

      !2nd derivative wrt to k wave vector and strain (piezoelectric tensor):
      ! - convert from cartesian to reduced coordinates (k point)
      ! - convert from reduced to cartesian coordinates (strain)
      ! - substract volume contribution
      ! - symetrize strain components
      if (choice==55.and.signs==1.and.paw_opt<=3) then
        ABI_MALLOC(work3,(2,3))
        ABI_MALLOC(work4,(2,3))
        ABI_MALLOC(work5,(2,3,6))
        ABI_MALLOC(work7,(2,3,6))
        ABI_MALLOC(work6,(2,3,3))
        do idat=1,ndat
          enlout_shift=(idat-1)*nnlout
          do ic=1,3 ! gamma
            work5=zero
            do jc=1,3 ! nu
              do ii=1,3 ! lambda
                mu=(gamma(jc,ii)-1)*3+1
                work5(1,jc,ii)=gmet2(ic,1)*enlout(enlout_shift+2*mu-1)+gmet2(ic,2)*enlout(enlout_shift+2*mu+1) &
       &         +gmet2(ic,3)*enlout(enlout_shift+2*mu+3)
                work5(2,jc,ii)=gmet2(ic,1)*enlout(enlout_shift+2*mu  )+gmet2(ic,2)*enlout(enlout_shift+2*mu+2) &
       &         +gmet2(ic,3)*enlout(enlout_shift+2*mu+4)
              end do
            end do
            work6=zero
            do jc=1,3 ! nu
              do ii=1,3 ! beta
                work6(1:cplex,ii,jc)=gprimd(ii,1)*work5(1:cplex,jc,1)+gprimd(ii,2)*work5(1:cplex,jc,2) &
       &         +gprimd(ii,3)*work5(1:cplex,jc,3)
              end do
            end do
            do jc=1,3 ! alpha
              do ii=1,3 ! beta
                mu=gamma(jc,ii)
                work7(1:cplex,ic,mu)=gprimd(jc,1)*work6(1:cplex,ii,1)+gprimd(jc,2)*work6(1:cplex,ii,2) &
       &         +gprimd(jc,3)*work6(1:cplex,ii,3)
              end do
            end do
          end do ! gamma

          do ii=1,3 ! alpha
            work3(1,ii)=gprimd(ii,1)*ddkk(2*1-1,idat)+gprimd(ii,2)*ddkk(2*2-1,idat) &
       &     +gprimd(ii,3)*ddkk(2*3-1,idat)
            work3(2,ii)=gprimd(ii,1)*ddkk(2*1  ,idat)+gprimd(ii,2)*ddkk(2*2  ,idat) &
       &     +gprimd(ii,3)*ddkk(2*3  ,idat)
          end do
          do ii=1,3 ! gamma
            work4(1,ii)=gmet2(ii,1)*ddkk(2*1-1,idat)+gmet2(ii,2)*ddkk(2*2-1,idat) &
       &     +gmet2(ii,3)*ddkk(2*3-1,idat)
            work4(2,ii)=gmet2(ii,1)*ddkk(2*1  ,idat)+gmet2(ii,2)*ddkk(2*2  ,idat) &
       &     +gmet2(ii,3)*ddkk(2*3  ,idat)
          end do

          do mu=1,6
            ii=alpha(mu) ! alpha
            ic=beta(mu) ! beta
            do jc=1,3 ! gamma
              work7(1:cplex,jc,mu)=work7(1:cplex,jc,mu)-half &
       &       *(gprimd(ic,jc)*work3(1:cplex,ii)+gprimd(ii,jc)*work3(1:cplex,ic))
              if (ii==ic) work7(1:cplex,jc,mu)=work7(1:cplex,jc,mu)-work4(1:cplex,jc)
            end do
          end do
          do mu=1,6 ! alpha,beta
            do nu=1,3 ! gamma
              mu0=3*(mu-1)+nu
              enlout(enlout_shift+2*mu0-1)=work7(1,nu,mu)
              enlout(enlout_shift+2*mu0  )=work7(2,nu,mu)
            end do
          end do
        end do !idat
        ABI_FREE(work3)
        ABI_FREE(work4)
        ABI_FREE(work5)
        ABI_FREE(work6)
        ABI_FREE(work7)
      end if

    end if !opernld

  end if ! choice>0

! Release memory

  if (allocated(enlk)) then
    ABI_FREE(enlk)
    ABI_FREE(fnlk)
    ABI_FREE(strnlk)
    ABI_FREE(ddkk)
  end if

  if (allocated(gmet2)) then
    ABI_FREE(gmet2)
  end if

  if(.not. local_vectproj) then
    ABI_FREE(projections)
  end if
  ABI_FREE(s_projections)
  ABI_FREE(vnl_projections)
  if(gemm_nonlop_is_distributed) then
    ABI_FREE(projs_local)
    ABI_FREE(projs_recv)
  end if
  if (allocated(dprojections)) then
    ABI_FREE(dprojections)
  end if
  if (allocated(s_dprojections)) then
    ABI_FREE(s_dprojections)
  end if
  if (allocated(vnl_dprojections)) then
    ABI_FREE(vnl_dprojections)
  end if
  if (allocated(d2projections)) then
    ABI_FREE(d2projections)
  end if
  if (allocated(temp_realvec)) then
    ABI_FREE(temp_realvec)
  end if

 end subroutine gemm_nonlop
!***

end module m_gemm_nonlop
!!***

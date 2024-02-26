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
 use m_fstrings,    only : itoa, ftoa, sjoin
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernla_gemm, only : opernla_gemm
 use m_opernlb_gemm, only : opernlb_gemm
 use m_opernld_ylm_allwf_cpu, only : opernld_ylm_allwf_cpu
 use m_pawcprj, only : pawcprj_type
 use m_kg, only : mkkpg

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

#if defined(HAVE_GPU_CUDA)
 use m_alloc_hamilt_gpu, only : gemm_nonlop_gpu_data
#endif

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr, c_int32_t, c_int64_t, c_float, c_double, c_size_t, c_loc
#endif

 implicit none

 private

 ! Use these routines in order: first call init, then call make_gemm_nonlop for each k point,
 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: init_gemm_nonlop
 public :: destroy_gemm_nonlop
 public :: make_gemm_nonlop
 public :: gemm_nonlop


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
 real(dp), pointer, save, public :: current_ikpt_projs(:, :, :)
 real(dp), pointer, save, public :: current_ikpt_projs_r(:, :, :)
 real(dp), pointer, save, public :: current_ikpt_projs_i(:, :, :)
 real(dp), pointer, save, public :: current_ikpt_dprojs(:, :, :)
 real(dp), pointer, save, public :: current_ikpt_dprojs_r(:, :, :)
 real(dp), pointer, save, public :: current_ikpt_dprojs_i(:, :, :)
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
   if(ikpt==current_ikpt_in_gpu) call free_ompgpu_current_ikpt()
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
 end if

 end subroutine free_gemm_nonlop_ikpt
!!***

!----------------------------------------------------------------------

 subroutine free_ompgpu_current_ikpt()

#ifdef HAVE_OPENMP_OFFLOAD

  if(current_ikpt_in_gpu == -1) return

  if(associated(current_ikpt_projs)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs)
    nullify(current_ikpt_projs)
  end if
  if(associated(current_ikpt_projs_r)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_r)
    nullify(current_ikpt_projs_r)
    nullify(current_ikpt_projs_i)
  end if
  if(associated(current_ikpt_dprojs)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs)
    nullify(current_ikpt_dprojs)
  end if
  if(associated(current_ikpt_dprojs_i)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_i)
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_r)
    nullify(current_ikpt_dprojs_r)
    nullify(current_ikpt_dprojs_i)
  end if

  current_ikpt_in_gpu=-1
#endif
 end subroutine free_ompgpu_current_ikpt

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

  integer :: nprojs,ndprojs,ngrads

  integer :: rank, nprocs, ierr, itypat
  integer :: nprojs_blk, nprojs_last_blk, nprojs_my_blk
  integer :: idir_beg,idir_end,idir_pert_,gpu_option_
  logical :: is_last_rank
  real(dp),allocatable :: dprojs_tmp(:,:,:),dprojs_r_tmp(:,:,:),dprojs_i_tmp(:,:,:)

! *************************************************************************

  ABI_CHECK(size(ph3d_k)>0,'nloalg(2)<0 not compatible with use_gemm_nonlop=1!')
!  ABI_CHECK((.not.my_compute_grad_strain).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
!  ABI_CHECK((.not.my_compute_grad_atom).or.size(kpg_k)>0,"kpg_k should be allocated to compute gradients!")
  idir_pert_=0; if(present(idir_pert)) idir_pert_=idir_pert
  gpu_option_=ABI_GPU_DISABLED; if(present(gpu_option)) gpu_option_=gpu_option

  ABI_CHECK(allocated(gemm_nonlop_kpt),"init_gemm_nonlop wasn't called prior make_gemm_nonlop !")

  if(choice>1) then
    ABI_CHECK(signs==2.or. choice==2 .or. choice==3 .or. choice==23,'signs/=2 and idir_pert not compatible with GEMM nonlop.')
    ABI_CHECK(gpu_option_/=ABI_GPU_LEGACY,'CUDA GEMM nonlop not compatible with respfn workloads.')
    ABI_CHECK(gpu_option_/=ABI_GPU_KOKKOS,'KOKKOS GEMM nonlop not compatible with respfn workloads.')
  end if

  call free_gemm_nonlop_ikpt(ikpt,gpu_option_)

  ! build nprojs, ngrads
  nprojs = 0 ; ngrads = 0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  if(signs==2) then
    ngrads=ngrads+1
  else ! signs==1
    if (choice==3) ngrads=ngrads+6
    if (choice==2) ngrads=ngrads+3
    if (choice==23) ngrads=ngrads+9
  end if
  if (nprojs>0) gemm_nonlop_kpt(ikpt)%nprojs = nprojs
  if (ngrads>0) gemm_nonlop_kpt(ikpt)%ngrads = ngrads

  if(istwf_k <= 1) then
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs, (2, npw, nprojs))
    !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%projs)
    if(ngrads>0) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs, (2, npw, nprojs*ngrads))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs)
    end if
  else
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1, npw, nprojs))
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1, npw, nprojs))
    !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%projs_r)
    !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%projs_i)
    if(ngrads>0) then
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_r, (1, npw, nprojs*ngrads))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs_r)
      ABI_MALLOC(gemm_nonlop_kpt(ikpt)%dprojs_i, (1, npw, nprojs*ngrads))
      !$OMP TARGET ENTER DATA MAP(alloc:gemm_nonlop_kpt(ikpt)%dprojs_i)
    end if
  end if

  call prep_projectors(npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k,&
  &                    ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k,&
  &                    nprojs,ngrads,choice,signs,idir_pert_,gpu_option_,&
  &                    gemm_nonlop_kpt(ikpt)%projs,&
  &                    gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i,&
  &                    gemm_nonlop_kpt(ikpt)%dprojs,&
  &                    gemm_nonlop_kpt(ikpt)%dprojs_r,gemm_nonlop_kpt(ikpt)%dprojs_i)


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
    call free_ompgpu_current_ikpt()
    gpu_nonlop_current_ikpt => gemm_nonlop_kpt(ikpt)
    if(allocated(gpu_nonlop_current_ikpt%projs)) then
      current_ikpt_projs   => gpu_nonlop_current_ikpt%projs
      !$OMP TARGET UPDATE TO(current_ikpt_projs)
      nullify(current_ikpt_projs_r)
      nullify(current_ikpt_projs_i)
    else
      current_ikpt_projs_r => gpu_nonlop_current_ikpt%projs_r
      current_ikpt_projs_i => gpu_nonlop_current_ikpt%projs_i
      !$OMP TARGET UPDATE TO(current_ikpt_projs_r)
      !$OMP TARGET UPDATE TO(current_ikpt_projs_i)
      nullify(current_ikpt_projs)
    end if
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      if(allocated(gpu_nonlop_current_ikpt%dprojs)) then
        current_ikpt_dprojs   => gpu_nonlop_current_ikpt%dprojs
        !$OMP TARGET UPDATE TO(current_ikpt_dprojs)
        nullify(current_ikpt_dprojs_r)
        nullify(current_ikpt_dprojs_i)
      else
        current_ikpt_dprojs_r => gpu_nonlop_current_ikpt%dprojs_r
        current_ikpt_dprojs_i => gpu_nonlop_current_ikpt%dprojs_i
        !$OMP TARGET UPDATE TO(current_ikpt_dprojs_i)
        !$OMP TARGET UPDATE TO(current_ikpt_dprojs_r)
        nullify(current_ikpt_dprojs)
      end if
    end if

    current_ikpt_in_gpu=ikpt
#endif
  end if
  gemm_nonlop_choice=choice

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
 &                          ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k,&
 &                          nprojs,ngrads,choice,signs,idir_pert,gpu_option,&
 &                          projs,projs_r,projs_i,&
 &                          dprojs,dprojs_r,dprojs_i)

  integer, intent(in) :: npw, lmnmax, ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  integer, intent(in) :: nprojs,ngrads,choice,signs,idir_pert,gpu_option
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)
  real(dp),intent(out) :: projs(:,:,:),dprojs(:,:,:)
  real(dp),intent(out) :: projs_r(:,:,:),dprojs_r(:,:,:)
  real(dp),intent(out) :: projs_i(:,:,:),dprojs_i(:,:,:)

  logical :: parity
  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer :: ndprojs
  integer :: il, ipw, idir, idir1, idir2, nkpg_local, ffnl_dir, dimffnl
  integer :: itypat, ilmn, nlmn, ia, iaph3d, igrad, shift, shift_grad
  integer :: lmn_beg,ibeg,iend,shift_do,nlmn_o,lmn_grad_beg
  real(dp):: wt
  real(dp),allocatable :: atom_projs(:,:,:), atom_dprojs(:,:,:,:), temp(:)
  real(dp),pointer :: kpg(:,:)

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
       call nvtxStartRange("prep_projectors")
#endif

!  if(gpu_option==ABI_GPU_OPENMP) then
!
!    if(istwf_k <= 1) then
!      projs = zero
!      if(ngrads>0) then
!        !$OMP TARGET ENTER DATA MAP(alloc:dprojs)
!        !$OMP TARGET DATA USE_DEVICE_PTR(dprojs)
!        call gpu_memset(c_loc(dprojs), 0, int(2,c_size_t)*npw*nprojs*ngrads*dp)
!        !$OMP END TARGET DATA
!      end if
!    else
!      !$OMP TARGET ENTER DATA MAP(alloc:projs_r,projs_i)
!      !$OMP TARGET DATA USE_DEVICE_PTR(projs_r,projs_i)
!      call gpu_memset(c_loc(projs_r), 0, int(1,c_size_t)*npw*nprojs*dp)
!      call gpu_memset(c_loc(projs_i), 0, int(1,c_size_t)*npw*nprojs*dp)
!      !$OMP END TARGET DATA
!      if(ngrads>0) then
!        !$OMP TARGET ENTER DATA MAP(alloc:dprojs_r,dprojs_i)
!        !$OMP TARGET DATA USE_DEVICE_PTR(dprojs_r,dprojs_i)
!        call gpu_memset(c_loc(dprojs_r), 0, int(1,c_size_t)*npw*nprojs*ngrads*dp)
!        call gpu_memset(c_loc(dprojs_i), 0, int(1,c_size_t)*npw*nprojs*ngrads*dp)
!        !$OMP END TARGET DATA
!      end if
!    end if
!
!  else
!
    if(istwf_k <= 1) then
      projs(:,:,:) = zero
      !!$OMP TARGET UPDATE TO(projs) IF(gpu_option==ABI_GPU_OPENMP)
      if(ngrads>0) then
        dprojs(:,:,:) = zero
        !!$OMP TARGET UPDATE TO(dprojs) IF(gpu_option==ABI_GPU_OPENMP)
      end if
    else
      projs_r(:,:,:) = zero
      projs_i(:,:,:) = zero
      !!$OMP TARGET UPDATE TO(projs_r,projs_i) IF(gpu_option==ABI_GPU_OPENMP)
      if(ngrads>0) then
        dprojs_r(:,:,:) = zero
        dprojs_i(:,:,:) = zero
        !!$OMP TARGET UPDATE TO(dprojs_r,dprojs_i) IF(gpu_option==ABI_GPU_OPENMP)
      end if
    end if

! end if

  iaph3d = 1
  wt=four_pi/sqrt(ucvol)
  dimffnl = size(ffnl_k, dim=2)
  ffnl_dir=1; if(dimffnl>2) ffnl_dir=idir_pert

  ndprojs = 0
  if (signs==1 .and. (choice==3 .or. choice==23)) then
    ndprojs = 3
  else if(signs==2 .and. (choice==5 .or. choice==51 .or. choice==3)) then
    ndprojs = 1
  end if

  ABI_MALLOC(atom_projs, (2, npw, lmnmax))
  if (ndprojs>0) then
    ABI_MALLOC(atom_dprojs, (2, npw, ndprojs, lmnmax))
  end if

  ABI_MALLOC(temp, (npw))

  ! Compute (k+G) vectors if needed
  nkpg_local=0
  if ((choice==2 .or. choice==3 .or. choice==23).and.size(kpg_k)==0) then
    nkpg_local=3
    ABI_MALLOC(kpg,(npw,nkpg_local))
    call mkkpg(kg_k,kpg,kpt_k,nkpg_local,npw)
  else
    kpg => kpg_k
  end if

  shift = 0 ; shift_grad = 0
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
      atom_projs(:,:,:) = zero
      do ipw=1, npw
        atom_projs(1,ipw, 1:nlmn_o) = wt * ffnl_k(ipw, 1, 1:nlmn_o, itypat)
      end do
      if (ndprojs>0) atom_dprojs(:,:,:,:) = zero
      if (signs==1 .and. (choice==3 .or. choice==23)) then
        do ipw=1, npw
          atom_dprojs(1,ipw, 1:ndprojs, 1:nlmn_o) = wt * ffnl_k(ipw, 2:ndprojs+1, 1:nlmn_o, itypat)
        end do
      end if
      if (signs==2 .and. (choice==3 .or. choice==5 .or. choice==51)) then
        do ipw=1, npw
          atom_dprojs(1,ipw, 1, 1:nlmn_o) = wt * ffnl_k(ipw, 1+ffnl_dir, 1:nlmn_o, itypat)
        end do
      end if

      ! multiply by (-i)^l
      do ilmn=1,nlmn_o
        il=mod(indlmn(1,ilmn, itypat),4);
        parity=(mod(il,2)==0)
        if (il>1) then
          ! multiply by -1
          atom_projs(:,:,ilmn) = -atom_projs(:,:,ilmn)
          if (ndprojs>0) atom_dprojs(:,:,:,ilmn) = -atom_dprojs(:,:,:,ilmn)
        end if
        if(.not. parity) then
          ! multiply by -i
          temp = atom_projs(2,:,ilmn)
          atom_projs(2,:,ilmn) = -atom_projs(1,:,ilmn)
          atom_projs(1,:,ilmn) =  temp
          if (ndprojs>0) then
            do idir=1,ndprojs
              temp = atom_dprojs(2,:,idir,ilmn)
              atom_dprojs(2,:,idir,ilmn) = -atom_dprojs(1,:,idir,ilmn)
              atom_dprojs(1,:,idir,ilmn) =  temp
            end do
          end if
        end if
      end do

      ! multiply by conj(ph3d)
      do ilmn=1,nlmn_o
        temp = atom_projs(1, :, ilmn)
        atom_projs(1, :, ilmn) = atom_projs(1, :, ilmn) * ph3d_k(1, :, iaph3d) &
        &                      + atom_projs(2, :, ilmn) * ph3d_k(2, :, iaph3d)
        atom_projs(2, :, ilmn) = atom_projs(2, :, ilmn) * ph3d_k(1, :, iaph3d) &
        &                      - temp                   * ph3d_k(2, :, iaph3d)
      end do
      if (ndprojs>0) then
        do ilmn=1,nlmn_o
          do idir=1,ndprojs
            temp = atom_dprojs(1, :, idir,ilmn)
            atom_dprojs(1, :, idir,ilmn) = atom_dprojs(1, :, idir,ilmn) * ph3d_k(1, :, iaph3d) &
            &                            + atom_dprojs(2, :, idir,ilmn) * ph3d_k(2, :, iaph3d)
            atom_dprojs(2, :, idir,ilmn) = atom_dprojs(2, :, idir,ilmn) * ph3d_k(1, :, iaph3d) &
            &                            - temp                         * ph3d_k(2, :, iaph3d)
          end do
        end do
      end if


      !! atom_projs is built, copy to projs

      if(istwf_k <= 1) then
        projs(1:2, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(:, :, lmn_beg:nlmn)
      else ! istwf_k>1
        projs_r(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(1, :, lmn_beg:nlmn)
        projs_i(1, :, shift+1:shift+(nlmn-(lmn_beg-1))) = atom_projs(2, :, lmn_beg:nlmn)
      end if

      !! Handling dprojs

      igrad=0

      if(signs==1 .and. (choice==3 .or. choice==23)) then
        if(istwf_k <= 1) then
          do idir=1,6
            idir1=alpha(idir);idir2=beta(idir)
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                dprojs(1:2, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     -half*(atom_dprojs(1:2, ipw, idir1, ilmn)*kpg(ipw,idir2) &
                &     +atom_dprojs(1:2, ipw, idir2, ilmn)*kpg(ipw,idir1))
              end do
            end do
            igrad=igrad+1
          end do
        else ! istwf_k>1
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
            igrad=igrad+1
          end do
        end if
      end if


      if(signs==1 .and. (choice==2 .or. choice==23)) then
        if(istwf_k <= 1) then
          do idir=1,3
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     +atom_projs(2, ipw, ilmn)*kpg(ipw,idir)*two_pi
                dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     -atom_projs(1, ipw, ilmn)*kpg(ipw,idir)*two_pi
              end do
            end do
            igrad=igrad+1
          end do
        else ! istwf_k>1
          do idir=1,3
            do ilmn=lmn_beg,nlmn
              do ipw=1,npw
                dprojs_r(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     +atom_projs(2, ipw, ilmn)*kpg(ipw,idir)*two_pi
                dprojs_i(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
                &     -atom_projs(1, ipw, ilmn)*kpg(ipw,idir)*two_pi
              end do
            end do
            igrad=igrad+1
          end do
        end if
      end if


      if(signs==2 .and. (choice==2)) then
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
            &     +atom_projs(2, ipw, ilmn)*kpg(ipw,idir_pert)*two_pi
            dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
            &     -atom_projs(1, ipw, ilmn)*kpg(ipw,idir_pert)*two_pi
          end do
        end do
        igrad=igrad+1
      end if


      if(signs==2 .and. (choice==3)) then
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
            &     +atom_dprojs(1, ipw, 1, ilmn)
            dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
            &     +atom_dprojs(2, ipw, 1, ilmn)
          end do
        end do
        igrad=igrad+1
      end if


      if(signs==2 .and. (choice==5 .or. choice==51)) then
        do ilmn=lmn_beg,nlmn
          do ipw=1,npw
            dprojs(1, ipw, shift_grad+nlmn*igrad+ilmn) = &
            &     +atom_dprojs(1, ipw, 1, ilmn)
            dprojs(2, ipw, shift_grad+nlmn*igrad+ilmn) = &
            &     +atom_dprojs(2, ipw, 1, ilmn)
          end do
        end do
        igrad=igrad+1
      end if

      iaph3d = iaph3d + 1
      shift_grad = shift_grad + ngrads*nlmn_o
      shift = shift + nlmn

    end do
  end do

  ABI_FREE(atom_projs)
  ABI_FREE(temp)
  if (allocated(atom_dprojs)) then
    ABI_FREE(atom_dprojs)
  end if
  if (nkpg_local>0) then
    ABI_FREE(kpg)
  end if
 end subroutine prep_projectors
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
  integer :: ii, ia, idat, igrad, nprojs, ngrads, shift, iatom, nlmn, ierr, ibeg, iend, ikpt
  integer :: cplex, cplex_enl, cplex_fac, proj_shift, grad_shift
  integer :: projs_beg,projs_end,dprojs_beg,dprojs_end
  integer :: nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  logical :: local_vectproj
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: sij_typ(:)
  real(dp), ABI_CONTIGUOUS pointer :: projections(:,:,:)
  real(dp), allocatable :: s_projections(:,:,:), vnl_projections(:,:,:)
  real(dp), allocatable :: dprojections(:,:,:), temp_realvec(:)
  real(dp), allocatable, target :: s_dprojections(:,:,:), vnl_dprojections(:,:,:)
  integer :: nprojs_my_blk
  integer :: rank, nprocs
  logical :: is_last
  real(dp), allocatable :: projs_local(:,:,:)
  real(dp), allocatable :: projs_recv(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_(:,:,:),dprojs_(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_r_(:,:,:),projs_i_(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: dprojs_r_(:,:,:),dprojs_i_(:,:,:)
  integer :: ngrads_tmp

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
&      (choice>3.and.choice/=7.and.choice/=23.and.signs==1) .or. &
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
  if((choice==2 .and. signs==2) .or. ((choice==2 .or. choice==3 .or. choice==23) .and. signs==1)) then
    projs_beg=atom_proj_shift+1
    projs_end=projs_beg+nprojs-1
    dprojs_beg=atom_proj_shift*ngrads+1
    dprojs_end=dprojs_beg+nprojs*ngrads-1
  end if
  if(istwf_k == 1) then
    projs_  => gemm_nonlop_kpt(ikpt)%projs(:,:,projs_beg:projs_end)
    dprojs_ => gemm_nonlop_kpt(ikpt)%dprojs(:,:,dprojs_beg:dprojs_end)
  else
    projs_r_  => gemm_nonlop_kpt(ikpt)%projs_r(:,:,projs_beg:projs_end)
    projs_i_  => gemm_nonlop_kpt(ikpt)%projs_i(:,:,projs_beg:projs_end)
    dprojs_r_ => gemm_nonlop_kpt(ikpt)%dprojs_r(:,:,dprojs_beg:dprojs_end)
    dprojs_i_ => gemm_nonlop_kpt(ikpt)%dprojs_i(:,:,dprojs_beg:dprojs_end)
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

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  if(.not. local_vectproj) then
    ABI_MALLOC(projections,(cplex, nprojs,nspinor*ndat))
  end if
  ABI_MALLOC(s_projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs,nspinor*ndat))

  if(.not. local_vectproj) projections = zero
  s_projections = zero
  vnl_projections = zero
  !if (ngrads>0) then
  ngrads_tmp=1; if (ngrads>0) ngrads_tmp=ngrads
    ABI_MALLOC(dprojections,(cplex, ngrads_tmp*nprojs,nspinor*ndat))
    ABI_MALLOC(s_dprojections,(cplex, ngrads_tmp*nprojs,nspinor*ndat))
    ABI_MALLOC(vnl_dprojections,(cplex_fac, ngrads_tmp*nprojs,nspinor*ndat))
    dprojections(:,:,:) = zero
    s_dprojections(:,:,:) = zero
    vnl_dprojections(:,:,:) = zero
  !end if

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

    call opernla_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,d2gxdt_dum_in,dprojections,projections,&
    &       idir,istwf_k,mpi_enreg,nd2gxdt,ngrads,&
    &       npwin,nspinor,signs,ndat,rank,&
    &       cpopt,nprocs,&
    &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
    &       vectin,projs_,dprojs_,&
    &       projs_r_,projs_i_,&
    &       dprojs_r_,dprojs_i_,&
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
      ndgxdt = ngrads
      ndgxdtfac = ngrads
      nd2gxdt = 0
      nd2gxdtfac = 0
      optder = 0;if (ndgxdtfac>0 .and. signs == 2) optder=1
      cplex_dgxdt(:) = 1;
      if(ndgxdt > 0) then
       if (choice==5.or.choice==51) cplex_dgxdt(:) = 2
      end if

      shift = 0
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

        do idat = 1,ndat
          call opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
&         dprojections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_dprojections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         s_dprojections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         d2gxdt_dum_in,d2gxdt_dum_out,d2gxdt_dum_out2,dimenl1,dimenl2,dimekbq,enl,&
&         projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         s_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         iatm,indlmn(:,:,itypat),itypat,lambda(idat),mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&         nattyp(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)
        end do

        shift = shift + nattyp(itypat)*nlmn
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
      call opernld_ylm_allwf_cpu(choice,cplex,cplex_fac,&
      &       dprojections,vnl_dprojections,s_dprojections,d2gxdt_dum_in,&
      &       enlout,projections,vnl_projections,s_projections,&
      &       ndat,nd2gxdt,ndgxdt,&
      &       ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
      &       gprimd,nattyp,mpi_enreg)
    end if !opernld

  end if ! choice>0

! Release memory
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
  if (allocated(temp_realvec)) then
    ABI_FREE(temp_realvec)
  end if

 end subroutine gemm_nonlop
!***

end module m_gemm_nonlop
!!***

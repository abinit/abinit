!!****m* ABINIT/m_gemm_nonlop_ompgpu
!! NAME
!! m_gemm_nonlop_ompgpu
!!
!! FUNCTION
!!  This module provides functions to compute the nonlocal operator by means of the BLAS GEMM
!!  routine. By treating ndat simultaneous wavefunctions, it is able to exploit BLAS3 routines,
!!  which leads to excellent CPU efficiency and OpenMP scalability.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2024 ABINIT group (MS)
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

module m_gemm_nonlop_ompgpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_abi_linalg

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#endif

#ifdef HAVE_GPU
  use m_gpu_toolbox
#endif

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm_ompgpu, only : opernlc_ylm_ompgpu
 use m_pawcprj, only : pawcprj_type
 use m_geometry, only : strconv
 use m_kg, only : mkkpg
 use m_gemm_nonlop

#if defined HAVE_MPI2
 use mpi
#endif

 implicit none

 private

 ! Use these routines in order: first call init, then call make_gemm_nonlop for each k point,
 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: init_gemm_nonlop_ompgpu
 public :: destroy_gemm_nonlop_ompgpu
 public :: make_gemm_nonlop_ompgpu
 public :: gemm_nonlop_ompgpu

 ! Those routines are here to assess memory requirements
 public :: gemm_nonlop_ompgpu_work_mem
 public :: gemm_nonlop_ompgpu_static_mem
!!***

!----------------------------------------------------------------------

#ifdef HAVE_OPENMP_OFFLOAD
 !(nkpt)
 real(dp), pointer :: current_ikpt_projs(:, :, :)
 ! (2, npw, nprojs)
 real(dp), pointer :: current_ikpt_projs_r(:, :, :)
 ! (1, npw, nprojs)
 real(dp), pointer :: current_ikpt_projs_i(:, :, :)
 ! (1, npw, nprojs)
 real(dp), pointer :: current_ikpt_dprojs(:, :, :)
 ! (2, npw, nprojs*ngrads)
 real(dp), pointer :: current_ikpt_dprojs_r(:, :, :)
 ! (1, npw, nprojs*ngrads)
 real(dp), pointer :: current_ikpt_dprojs_i(:, :, :)
 ! (1, npw, nprojs*ngrads)
 type(gemm_nonlop_type), pointer :: gpu_nonlop_current_ikpt
 integer, save :: current_ikpt_in_gpu=-1
 integer, save :: gpu_initialised=0
 integer, save :: mod__cplex=0
 integer, save :: mod__nkpt=0
 integer, save :: mod__ndat=0
 real(dp), save, allocatable, target :: projections(:,:,:), s_projections(:,:,:), vnl_projections(:,:,:), dprojections(:,:,:)
 real(dp), save, allocatable :: temp_realvec_r(:),temp_realvec_i(:)
 real(dp), save, allocatable :: sij_typ(:,:)

#endif

contains

 function gemm_nonlop_ompgpu_work_mem(istwfk, ndat, npw, indlmn, nattyp, ntypat, lmnmax) result(req_mem)
   implicit none

   integer, intent(in) :: istwfk, ndat, npw, ntypat, lmnmax
   integer, intent(in) :: indlmn(:,:,:), nattyp(ntypat)

   integer :: nprojs, cplex, itypat
   real(dp) :: req_mem

! *************************************************************************

   cplex=2;if (istwfk>1) cplex=1
   nprojs=0
   do itypat=1,ntypat
     nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
   end do

   req_mem = 0

   if(cplex == 1) then
     req_mem = req_mem + dp * npw * ndat ! temp_realvec_r
     req_mem = req_mem + dp * npw * ndat ! temp_realvec_i
   end if

   req_mem = req_mem + dp * lmnmax * (lmnmax+1)/2 * ntypat  ! sij_typ

   req_mem = req_mem + dp * cplex * int(nprojs, c_size_t) * int(ndat, c_size_t)  ! projections
   req_mem = req_mem + dp * cplex * int(nprojs, c_size_t) * int(ndat, c_size_t)  ! s_projections
   req_mem = req_mem + dp * cplex * int(nprojs, c_size_t) * int(ndat, c_size_t)  ! vnl_projections

 end function gemm_nonlop_ompgpu_work_mem

!----------------------------------------------------------------------

 function gemm_nonlop_ompgpu_static_mem(npw, indlmn, nattyp, ntypat, nblocks,&
     compute_grad_strain,compute_grad_atom) result(req_mem)
   implicit none

   integer, intent(in) :: npw, ntypat, nblocks
   integer, intent(in) :: indlmn(:,:,:), nattyp(ntypat)
   logical, intent(in), optional :: compute_grad_strain,compute_grad_atom

   integer :: nprojs, ngrads, itypat
   logical :: my_compute_grad_strain,my_compute_grad_atom
   integer(kind=c_size_t) :: req_mem

! *************************************************************************

   my_compute_grad_strain=.false. ; if (present(compute_grad_strain)) my_compute_grad_strain=compute_grad_strain
   my_compute_grad_atom=.false. ; if (present(compute_grad_atom)) my_compute_grad_atom=compute_grad_atom

   nprojs = 0
   do itypat=1,ntypat
     nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
   end do
   nprojs = nprojs / nblocks

   ngrads = 0
   if (my_compute_grad_strain) ngrads = ngrads + 6
   if (my_compute_grad_atom)   ngrads = ngrads + 3

   req_mem = 0

   if(nblocks>1) then
     req_mem = req_mem + dp * 2 * int(npw, c_size_t) * int(nprojs, c_size_t)          !projs_recv
   end if
   ! projs or projs_r + projs_i
   req_mem = req_mem + 2 * dp * int(npw, c_size_t) * int(nprojs, c_size_t)
   if(ngrads>0) then
     ! dprojs or dprojs_r + dprojs_i
     req_mem = req_mem + 2 * dp * int(npw, c_size_t) * int(ngrads, c_size_t) * int(nprojs, c_size_t)
     if(nblocks>1) then
       req_mem = req_mem + dp * 2 * int(npw, c_size_t) * int(ngrads, c_size_t)*int(nprojs, c_size_t)   !dprojs_recv
     end if
   end if

 end function gemm_nonlop_ompgpu_static_mem

!----------------------------------------------------------------------

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! - AMD GPUs    : FC_LLVM + HIP
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

 subroutine refresh_gemm_nonlop_kpt_ompgpu(ikpt)
  integer,intent(in) :: ikpt

  if(mod__cplex == 2) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs)
    end if
  else if(mod__cplex == 1) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_r)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_i)
      !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_r)
    end if
  end if

  gpu_nonlop_current_ikpt => gemm_nonlop_kpt(ikpt)
  current_ikpt_projs   => gemm_nonlop_kpt(ikpt)%projs
  current_ikpt_projs_r => gemm_nonlop_kpt(ikpt)%projs_r
  current_ikpt_projs_i => gemm_nonlop_kpt(ikpt)%projs_i
  if(gpu_nonlop_current_ikpt%ngrads /= -1) then
    current_ikpt_dprojs   => gemm_nonlop_kpt(ikpt)%dprojs
    current_ikpt_dprojs_r => gemm_nonlop_kpt(ikpt)%dprojs_r
    current_ikpt_dprojs_i => gemm_nonlop_kpt(ikpt)%dprojs_i
  end if

  if(mod__cplex == 2) then
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs)
    end if
  else if(mod__cplex == 1) then
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs_i)
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs_r)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs_i)
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs_r)
    end if
  end if
  current_ikpt_in_gpu=ikpt

 end subroutine refresh_gemm_nonlop_kpt_ompgpu

 subroutine alloc_gemm_nonlop_kpt_ompgpu(ikpt, cplex)

  integer,intent(in) :: ikpt, cplex

  character(len=500) :: msg


! *************************************************************************

  if(current_ikpt_in_gpu == ikpt) then
    msg="GPU nonlop arrays are already initialized for K-Point %d. Redundant call"
    ABI_ERROR(msg)
  end if
  if(gemm_nonlop_kpt(ikpt)%nprojs == -1) then
    msg="Was asked to work on a K-point index for which arrays aren't built.\n" // &
&    "Requested K-point index was %d.\n" // &
&    "Please make sure that make_gemm_nonlop_ompgpu is called for this K-point before reaching this routine."
    ABI_ERROR(msg)
  end if
  if(current_ikpt_in_gpu /= ikpt) then
    call free_gemm_nonlop_kpt_ompgpu ! FIXME : Dealloc because projs arrays sizes may varry among K-Points ?
  end if
  gpu_nonlop_current_ikpt => gemm_nonlop_kpt(ikpt)
  current_ikpt_projs   => gemm_nonlop_kpt(ikpt)%projs
  current_ikpt_projs_r => gemm_nonlop_kpt(ikpt)%projs_r
  current_ikpt_projs_i => gemm_nonlop_kpt(ikpt)%projs_i
  if(gpu_nonlop_current_ikpt%ngrads /= -1) then
    current_ikpt_dprojs   => gemm_nonlop_kpt(ikpt)%dprojs
    current_ikpt_dprojs_r => gemm_nonlop_kpt(ikpt)%dprojs_r
    current_ikpt_dprojs_i => gemm_nonlop_kpt(ikpt)%dprojs_i
  end if

  if(cplex == 2) then
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs)
    end if
  else if(cplex == 1) then
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs_i)
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs_r)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs_i)
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs_r)
    end if
  end if

  mod__cplex = cplex
  current_ikpt_in_gpu=ikpt

 end subroutine alloc_gemm_nonlop_kpt_ompgpu
!!***

!----------------------------------------------------------------------

 subroutine free_gemm_nonlop_kpt_ompgpu()

  if(current_ikpt_in_gpu == -1) return

  if(mod__cplex == 2) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs)
    end if
  end if
  if(mod__cplex == 1) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_r)
    if(gpu_nonlop_current_ikpt%ngrads /= -1) then
      !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_i)
      !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_r)
    end if
  end if

  current_ikpt_in_gpu=-1

 end subroutine free_gemm_nonlop_kpt_ompgpu

!----------------------------------------------------------------------

 subroutine alloc_work_buffers(cplex, cplex_fac, ndat, nprojs, ntypat, lmnmax, npw)

  integer,intent(in) :: cplex, cplex_fac, ndat, nprojs, ntypat, lmnmax, npw

! *************************************************************************

  call free_work_buffers()

  mod__ndat=ndat

  if(cplex == 1) then
    ABI_MALLOC(temp_realvec_r,(npw*ndat))
    ABI_MALLOC(temp_realvec_i,(npw*ndat))
    !$OMP TARGET ENTER DATA MAP(alloc:temp_realvec_r,temp_realvec_i)
  end if

  ABI_MALLOC(sij_typ,(lmnmax*(lmnmax+1)/2,ntypat))
  !$OMP TARGET ENTER DATA MAP(alloc:sij_typ)

  ABI_MALLOC(projections,(cplex, nprojs, ndat))
  ABI_MALLOC(s_projections,(cplex, nprojs, ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs, ndat))

  !FIXME Smarter buffer management
#ifdef HAVE_GPU_HIP
  !Work buffer allocated once to save time in HIP (alloc costful)
  !$OMP TARGET ENTER DATA MAP(alloc:projections,s_projections,vnl_projections)
#endif

  gpu_initialised=1

 end subroutine alloc_work_buffers
!!***

!----------------------------------------------------------------------

 subroutine free_work_buffers()


  !$OMP TARGET EXIT DATA MAP(delete:sij_typ)
  if(allocated(sij_typ)) then
    ABI_FREE(sij_typ)
  end if

  if(allocated(temp_realvec_r)) then
    !$OMP TARGET EXIT DATA MAP(delete:temp_realvec_r,temp_realvec_i)
    ABI_FREE(temp_realvec_r)
    ABI_FREE(temp_realvec_i)
  end if

  !FIXME Smarter buffer management
#ifdef HAVE_GPU_HIP
  !$OMP TARGET EXIT DATA MAP(delete:projections,s_projections,vnl_projections)
#endif

  if(allocated(projections)) then
    ABI_FREE(projections)
    ABI_FREE(s_projections)
    ABI_FREE(vnl_projections)
  end if

  mod__ndat=0
  gpu_initialised=0

 end subroutine free_work_buffers

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/init_gemm_nonlop_ompgpu
!! NAME
!! init_gemm_nonlop_ompgpu
!!
!! FUNCTION
!! Initalization of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine init_gemm_nonlop_ompgpu(nkpt)

  integer,intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

  ! TODO only allocate the number of kpt treated by this proc
  ABI_MALLOC(gemm_nonlop_kpt, (nkpt))
  do ikpt=1,nkpt
    gemm_nonlop_kpt(ikpt)%nprojs = -1
    gemm_nonlop_kpt(ikpt)%ngrads = -1
  end do
  mod__nkpt=nkpt

 end subroutine init_gemm_nonlop_ompgpu
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/destroy_gpu_nonlop_ompgpu
!! NAME
!! destroy_gemm_nonlop_ompgpu
!!
!! FUNCTION
!! Destruction of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine destroy_gemm_nonlop_ompgpu()

  integer :: ikpt
  character(len=200) :: msg

! *************************************************************************

  if(mod__nkpt==0) then
    msg="GPU nonlop arrays aren't initialized. Redundant call"
    ABI_ERROR(msg)
  end if

  call free_work_buffers()
  call free_gemm_nonlop_kpt_ompgpu()
  ! TODO add cycling if kpt parallelism
  do ikpt = 1,mod__nkpt
    call free_gemm_nonlop_ompgpu_ikpt(ikpt)
  end do

  ABI_FREE(gemm_nonlop_kpt)
  mod__nkpt=0

 end subroutine destroy_gemm_nonlop_ompgpu
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/free_gemm_nonlop_ompgpu_ikpt
!! NAME
!! free_destroy_gemm_nonlop_ompgpuikpt
!!
!! FUNCTION
!! Release memory for one kpt value of the gemm_nonlop_kpt array
!!
!! INPUTS
!! ikpt= index of gemm_nonlop_kptto be released
!!
!! SOURCE
 subroutine free_gemm_nonlop_ompgpu_ikpt(ikpt)

  integer,intent(in) :: ikpt

! *************************************************************************

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

 end subroutine free_gemm_nonlop_ompgpu_ikpt
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/make_gemm_nonlop_ompgpu
!! NAME
!! make_gemm_nonlop_ompgpu
!!
!! FUNCTION
!! Build the gemm_nonlop array
!!
!! INPUTS
!!
!! SOURCE

 subroutine make_gemm_nonlop_ompgpu(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
&                            ph3d_k,kpt_k,kg_k,kpg_k, &
&                            compute_grad_strain,compute_grad_atom) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax,ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  logical, intent(in), optional :: compute_grad_strain,compute_grad_atom
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)

  integer :: cplex
! *************************************************************************

  call free_gemm_nonlop_kpt_ompgpu()
  call make_gemm_nonlop(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
&                            ph3d_k,kpt_k,kg_k,kpg_k, &
&                            compute_grad_strain=compute_grad_strain,compute_grad_atom=compute_grad_atom)
  cplex=2;if (istwf_k>1) cplex=1
  call alloc_gemm_nonlop_kpt_ompgpu(ikpt,cplex)

 end subroutine make_gemm_nonlop_ompgpu
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/gemm_nonlop_ompgpu_distributed_gemm_opernla
!! NAME
!! gemm_nonlop_ompgpu_distributed_gemm_opernla
!!
!! FUNCTION
!! Distributed version of "opernla" GEMM called in gemm_nonlop.
!!
!! INPUTS
!!
!! SOURCE
 subroutine gemm_nonlop_ompgpu_distributed_gemm_opernla(rank,nprocs,npwin,ndat,nspinor,&
 &                                               nprojs,nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex,beta,&
 &                                               projs_local,projs_recv,vectin,cprojs)
   integer,  intent(in)     :: rank,nprocs,npwin,ndat,nspinor
   integer,  intent(in)     :: nprojs,nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex
   real(dp), intent(in)     :: vectin(*)
   complex(dpc), intent(in) :: beta
   real(dp), intent(inout), target  :: projs_local(cplex,npwin,nprojs_last_blk),projs_recv(cplex,npwin,nprojs_last_blk)
   real(dp), intent(out)  :: cprojs(:,:,:)

   !Local variables
   integer :: iblock,ibeg,iend,req(2),ierr,nprojs_cur_blk,rank_prev,rank_next,idat
   real(dp), ABI_CONTIGUOUS pointer :: recv_buf(:,:,:), work_buf(:,:,:)
   real(dp), ABI_CONTIGUOUS pointer :: recv_buf_f(:), work_buf_f(:)


! *************************************************************************

   rank_next=modulo(rank + 1,nprocs)
   rank_prev=rank - 1
   if(rank_prev == -1) rank_prev = nprocs - 1
   do iblock=1,nprocs

     if(rank+iblock == nprocs) then
       nprojs_cur_blk = nprojs_last_blk
     else
       nprojs_cur_blk = nprojs_blk
     end if

     if(modulo(iblock,2)==1) then
       work_buf => projs_local(1:cplex,1:npwin,1:nprojs_last_blk)
       recv_buf => projs_recv(1:cplex,1:npwin,1:nprojs_last_blk)
     else
       work_buf => projs_recv(1:cplex,1:npwin,1:nprojs_last_blk)
       recv_buf => projs_local(1:cplex,1:npwin,1:nprojs_last_blk)
     end if

#ifndef HAVE_GPU_MPI

       ! GPU-aware MPI not available : perform MPI comms on CPU
       call xmpi_isend(work_buf,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call xmpi_irecv(recv_buf,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
       !$OMP TARGET UPDATE TO(work_buf)

#else

       ! GPU-aware MPI available : pass GPU buffers to MPI
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,recv_buf)
       call c_f_pointer(c_loc(work_buf), work_buf_f, [cplex*npwin*nprojs_last_blk])
       call c_f_pointer(c_loc(recv_buf), recv_buf_f, [cplex*npwin*nprojs_last_blk])
       call MPI_ISEND(work_buf_f,cplex*npwin*nprojs_cur_blk,MPI_DOUBLE_PRECISION,&
       &    rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call MPI_IRECV(recv_buf_f,cplex*npwin*nprojs_cur_blk,MPI_DOUBLE_PRECISION,&
       &    rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
       !$OMP END TARGET DATA

#endif


     ibeg = 1 + modulo(rank+iblock-1,nprocs)*nprojs_blk
     iend = ibeg+nprojs_cur_blk-1
     if(cplex == 2) then
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,vectin,cprojs)
       call abi_gpu_xgemm(cplex, 'C','N', &
               nprojs_cur_blk, ndat*nspinor, npwin, cone, &
               c_loc(work_buf), npwin, &
               c_loc(vectin), npwin, &
               beta, &
               c_loc(cprojs(1,ibeg,1)), nprojs)
       !$OMP END TARGET DATA
     else
       call DGEMM('T', 'N', nprojs_cur_blk, ndat*nspinor, npwin, one, &
       &          work_buf, npwin, &
       &          vectin, npwin, real(beta), cprojs(:,ibeg:iend,:), nprojs_cur_blk)
     end if

     call xmpi_waitall(req,ierr)

   end do

   if(modulo(iblock,2)==1) then
#ifndef HAVE_GPU_MPI
       call DCOPY(cplex*npwin*nprojs_cur_blk, recv_buf, 1, work_buf, 1)
#else
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,recv_buf)
       call copy_gpu_to_gpu(c_loc(work_buf), c_loc(recv_buf), INT(cplex, c_size_t)*npwin*nprojs_last_blk*dp)
       !$OMP END TARGET DATA
#endif
   end if
 end subroutine gemm_nonlop_ompgpu_distributed_gemm_opernla
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/gemm_nonlop_ompgpu_distributed_gemm_opernlb
!! NAME
!! gemm_nonlop_ompgpu_distributed_gemm_opernlb
!!
!! FUNCTION
!! Distributed version of "opernlb" GEMM called in gemm_nonlop.
!!
!! INPUTS
!!
!! SOURCE
 subroutine gemm_nonlop_ompgpu_distributed_gemm_opernlb(rank,nprocs,npwout,ndat,nspinor,&
 &                                               nprojs,nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex,&
 &                                               projs_local,projs_recv,cprojs_s_vnl,vectout)
   integer,  intent(in)     :: rank,nprocs,npwout,ndat,nspinor
   integer,  intent(in)     :: nprojs,nprojs_blk,nprojs_last_blk,nprojs_my_blk,cplex
   real(dp), intent(in)     :: cprojs_s_vnl(:,:,:)
   real(dp), intent(inout), target  :: projs_local(cplex,npwout,nprojs_last_blk),projs_recv(cplex,npwout,nprojs_last_blk)
   real(dp), intent(out)    :: vectout(*)

   !Local variables
   integer :: iblock,ibeg,iend,req(2),ierr,nprojs_cur_blk,rank_prev,rank_next
   complex(dpc) :: beta
   real(dp), ABI_CONTIGUOUS pointer :: recv_buf(:,:,:), work_buf(:,:,:)
   real(dp), ABI_CONTIGUOUS pointer :: recv_buf_f(:), work_buf_f(:)

! *************************************************************************

   rank_next=modulo(rank + 1,nprocs)
   rank_prev=rank - 1
   if(rank_prev == -1) rank_prev = nprocs - 1

   beta = czero

   do iblock=1,nprocs

     if(rank+iblock == nprocs) then
       nprojs_cur_blk = nprojs_last_blk
     else
       nprojs_cur_blk = nprojs_blk
     end if

     if(modulo(iblock,2)==1) then
       work_buf => projs_local(1:cplex,1:npwout,1:nprojs_last_blk)
       recv_buf => projs_recv(1:cplex,1:npwout,1:nprojs_last_blk)
     else
       work_buf => projs_recv(1:cplex,1:npwout,1:nprojs_last_blk)
       recv_buf => projs_local(1:cplex,1:npwout,1:nprojs_last_blk)
     end if

#ifndef HAVE_GPU_MPI

       ! GPU-aware MPI not available : perform MPI comms on CPU
       call xmpi_isend(work_buf,rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call xmpi_irecv(recv_buf,rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
       !$OMP TARGET UPDATE TO(work_buf)

#else

       ! GPU-aware MPI available : pass GPU buffers to MPI
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,recv_buf)
       call c_f_pointer(c_loc(work_buf), work_buf_f, [cplex*npwout*nprojs_last_blk])
       call c_f_pointer(c_loc(recv_buf), recv_buf_f, [cplex*npwout*nprojs_last_blk])
       call MPI_ISEND(work_buf_f,cplex*npwout*nprojs_cur_blk,MPI_DOUBLE_PRECISION,&
       &    rank_prev,iblock,gemm_nonlop_block_comm,req(1),ierr)
       call MPI_IRECV(recv_buf_f,cplex*npwout*nprojs_cur_blk,MPI_DOUBLE_PRECISION,&
       &    rank_next,iblock,gemm_nonlop_block_comm,req(2),ierr)
       !$OMP END TARGET DATA

#endif

     ibeg = 1 + modulo(rank+iblock-1,nprocs)*nprojs_blk
     iend = ibeg+nprojs_cur_blk-1
     if(cplex==2) then
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,vectout,cprojs_s_vnl)
       call abi_gpu_xgemm(cplex, 'N','N', &
               npwout, ndat*nspinor, nprojs_cur_blk, cone, &
               c_loc(work_buf), npwout, &
               c_loc(cprojs_s_vnl(1,ibeg,1)), nprojs, &
               beta, &
               c_loc(vectout), npwout)
       !$OMP END TARGET DATA
     else
       call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs_cur_blk, one, &
       &          work_buf, npwout, &
       &          cprojs_s_vnl(:,ibeg:iend,:), nprojs_cur_blk, real(beta), vectout, npwout)
     end if

     beta = cone

     call xmpi_waitall(req,ierr)

   end do

   if(modulo(iblock,2)==1) then
#ifndef HAVE_GPU_MPI
       call DCOPY(cplex*npwout*nprojs_cur_blk, recv_buf, 1, work_buf, 1)
#else
       !$OMP TARGET DATA USE_DEVICE_PTR(work_buf,recv_buf)
       call copy_gpu_to_gpu(c_loc(work_buf), c_loc(recv_buf), INT(cplex, c_size_t)*npwout*nprojs_last_blk*dp)
       !$OMP END TARGET DATA
#endif
   end if

 end subroutine gemm_nonlop_ompgpu_distributed_gemm_opernlb
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/gemm_nonlop_ompgpu
!! NAME
!! gemm_nonlop_ompgpu
!!
!! FUNCTION
!! Replacement of nonlop. same prototype as nonlop although not all options are implemented.
!!
!! INPUTS
!! [gpu_option] = GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!
!! SOURCE
 subroutine gemm_nonlop_ompgpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm
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
  integer :: ii, ia, idat, igrad, nprojs, ngrads, shift, iatom, nlmn, ierr, ibeg, iend, i1, i2, i, ikpt
  integer :: cplex, cplex_enl, cplex_fac, proj_shift, grad_shift, nattyp_i
  integer :: enlout_shift, force_shift, nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  real(dp) :: esum,esumk(9)
  real(dp) :: work(6)
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: enlk(:)
  real(dp),pointer :: projections_ptr(:,:,:)
  integer :: nprojs_my_blk, ipw, iproj, iblock, nprojs_blk
  integer :: rank, nprocs, req(2)
  logical :: is_last
  real(dp), allocatable :: projs_recv(:,:,:), dprojs_recv(:,:,:)
  logical :: local_vectproj
  logical :: transfer_vectin,transfer_vectout,transfer_svectout
  character(len=500) :: msg
  integer(C_SIZE_T) :: byte_count
  real(dp), ABI_CONTIGUOUS pointer :: vectin_(:,:),vectout_(:,:),svectout_(:,:)

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,kpgin,kpgout/))
  ABI_UNUSED((/ph1d(1,1),ph3din,ph3dout/))
  ABI_UNUSED((/phkxredin,phkxredout,ucvol/))
  ABI_UNUSED((/mgfft,mpsang,mpssoang/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,gpu_option/))

  ! Check supported options
  if ( (choice>1.and.choice/=7.and.signs==2) .or. &
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

  nprojs = gemm_nonlop_kpt(ikpt)%nprojs
  ngrads = gemm_nonlop_kpt(ikpt)%ngrads
  nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_blk

  if(gemm_nonlop_is_distributed) then
    rank = xmpi_comm_rank(gemm_nonlop_block_comm); nprocs = xmpi_comm_size(gemm_nonlop_block_comm)
    is_last = (rank==nprocs-1)
    if(is_last) nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_last_blk
  end if

  if(choice==1) ngrads=0
  ABI_CHECK(ngrads>=3.or.choice/=2 ,"Bad allocation in gemm_nonlop (2)!")
  ABI_CHECK(ngrads>=6.or.choice/=3 ,"Bad allocation in gemm_nonlop (3)!")
  ABI_CHECK(ngrads>=9.or.choice/=23,"Bad allocation in gemm_nonlop (23)!")

  ! Allocate and copy GPU buffers if user doesn't manage them
  transfer_vectin=.not. xomp_target_is_present(c_loc(vectin)) &
      .and. (cpopt < 2 .or. (choice/=7 .and.paw_opt >=3))
  transfer_vectout=.not. xomp_target_is_present(c_loc(vectout)) &
      .and. (signs==2 .and. (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4))
  transfer_svectout=.not. xomp_target_is_present(c_loc(svectout)) &
      .and. (signs==2 .and. (paw_opt == 3 .or. paw_opt == 4))
  !$OMP TARGET ENTER DATA MAP(to:vectin)      IF(transfer_vectin)
  !$OMP TARGET ENTER DATA MAP(alloc:vectout)  IF(transfer_vectout)
  !$OMP TARGET ENTER DATA MAP(alloc:svectout) IF(transfer_svectout)

  !FIXME These seemingly useless pointers are used in BLAS operations for
  ! working around a bug in AOMP LLVM misreading mapped device pointers.
  ! I chose to generalise the workaround to avoid dupplicating each
  ! BLAS call specifically for handling AOMP LLVM.
  vectin_   => vectin
  vectout_  => vectout
  svectout_ => svectout

  !$OMP TARGET ENTER DATA MAP(to:atindx1,indlmn,enl)
  if(gpu_initialised == 0 .or. mod__ndat /= ndat*nspinor) then
    call alloc_work_buffers(cplex, cplex_fac,&
&        nspinor*ndat, nprojs, ntypat, lmnmax, MAX(npwin,npwout))
    if(paw_opt>=2 .and. choice > 0 .and. choice /= 7) then
      if (cplex_enl==1) then
        do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)
          do ilmn=1,nlmn*(nlmn+1)/2
            sij_typ(ilmn,itypat)=sij(ilmn,itypat)
          end do
        end do
      else
        do itypat=1, ntypat
          nlmn=count(indlmn(3,:,itypat)>0)
          do ilmn=1,nlmn*(nlmn+1)/2
            sij_typ(ilmn,itypat)=sij(2*ilmn-1,itypat)
          end do
        end do
      end if
    end if
    !$OMP TARGET UPDATE TO(sij_typ)
  end if
  if(current_ikpt_in_gpu /= gemm_nonlop_ikpt_this_proc_being_treated) then
    call refresh_gemm_nonlop_kpt_ompgpu(gemm_nonlop_ikpt_this_proc_being_treated)
  end if

   ! If vectproj is provided, use it for further calculations, use static array otherwise
   projections_ptr => projections
   local_vectproj=.false.
   if(PRESENT(vectproj)) then
     if(size(vectproj)>1) local_vectproj=.true.
   end if
   if (local_vectproj) projections_ptr => vectproj

#ifdef HAVE_GPU_CUDA
  !Work buffers allocated at each call to save memory in CUDA
  !$OMP TARGET ENTER DATA MAP(alloc:s_projections,vnl_projections)
  if(.not. local_vectproj) then
    !$OMP TARGET ENTER DATA MAP(alloc:projections_ptr)
  end if
#endif

  if(gemm_nonlop_is_distributed) then
    ABI_MALLOC(projs_recv,   (cplex, npwin, gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
    !$OMP TARGET ENTER DATA MAP(alloc:projs_recv)
    if (signs==1.and.ngrads>0) then
      ABI_MALLOC(dprojs_recv,   (cplex, npwin, ngrads*gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
      !$OMP TARGET ENTER DATA MAP(alloc:dprojs_recv)
    end if
  end if

  ! These will store the non-local factors for vectin, svectout and vectout respectively
#if defined HAVE_GPU_CUDA
  byte_count=sizeof(projections_ptr)
  !$OMP TARGET DATA USE_DEVICE_PTR(projections_ptr,s_projections,vnl_projections)
  if(cpopt < 2) call gpu_memset(c_loc(projections_ptr),     0, byte_count)
  call gpu_memset(c_loc(s_projections),   0, byte_count)
  call gpu_memset(c_loc(vnl_projections), 0, byte_count)
  !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
  if(cpopt < 2) then
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP& MAP(to:projections_ptr) PRIVATE(i1,i2)
    do i2=1, nspinor*ndat
      do i1=1, nprojs
        projections_ptr(:,i1,i2) = zero
      end do
    end do
  end if
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
  !$OMP& MAP(to:s_projections,vnl_projections) PRIVATE(i1,i2)
  do i2=1, nspinor*ndat
    do i1=1, nprojs
      s_projections(:,i1,i2) = zero
      vnl_projections(:,i1,i2) = zero
    end do
  end do
#endif

  if (signs==1.and.ngrads>0) then
    ABI_MALLOC(dprojections,(cplex, ngrads*nprojs, nspinor*ndat))
    !$OMP TARGET ENTER DATA MAP(alloc:dprojections)
#if defined HAVE_GPU_CUDA
    byte_count=sizeof(dprojections)
    !$OMP TARGET DATA USE_DEVICE_PTR(dprojections)
    call gpu_memset(c_loc(dprojections), 0, byte_count)
    !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP& MAP(to:dprojections) PRIVATE(i1,i2)
    do i2=1, nspinor*ndat
      do i1=1, nprojs*ngrads
        dprojections(:,i1,i2) = zero
      end do
    end do
#endif
    if(choice==1.or.choice==3.or.choice==23) then
      ABI_MALLOC(enlk,(ndat))
      enlk=zero
      !$OMP TARGET ENTER DATA MAP(to:enlk,nattyp)
    end if
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
    !$OMP TARGET ENTER DATA MAP(to:enlout)
  end if

  !$OMP TASKWAIT
  if(cpopt >= 2) then
    ! retrieve from cprjin
    if(.not. local_vectproj .and. cpopt/=3) then
      !TODO This use-case is extremely painful for GEMM OpenGPU nonlop performance
      ABI_WARNING("Inefficient call of OpenGPU nonlop. Was vectproj provided with OpenMP mapping?")
      !$OMP TARGET UPDATE FROM(projections_ptr)
      !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,nlmn)
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          projections_ptr(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
          shift = shift + nlmn
        end do
      end do
      !$OMP TARGET UPDATE TO(projections_ptr)
    end if
    if(cpopt==4.and.allocated(dprojections)) then
      !$OMP TARGET UPDATE FROM(dprojections)
      ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (1)")
      !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,igrad,nlmn)
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn  = cprjin(iatom, idat)%nlmn
          do igrad=1,ngrads
            dprojections(1:cplex, shift+1:shift+nlmn, idat) = &
&                   cprjin(iatom, idat)%dcp(1:cplex,igrad,1:ilmn)
            shift = shift + nlmn
          end do
        end do
      end do
      !$OMP TARGET UPDATE TO(dprojections)
    end if
  else
    ! opernla
    if(cplex == 2) then
      if(.not. gemm_nonlop_is_distributed) then
        !$OMP TARGET DATA USE_DEVICE_PTR(projections_ptr,current_ikpt_projs,vectin_)
        call abi_gpu_xgemm(cplex, 'C', 'N', nprojs, ndat*nspinor, npwin, cone, &
&                c_loc(current_ikpt_projs), npwin,&
&                c_loc(vectin_), npwin, czero, c_loc(projections_ptr), nprojs)
        !$OMP END TARGET DATA
        if(signs==1.and.ngrads>0) then
          !$OMP TARGET DATA USE_DEVICE_PTR(dprojections,current_ikpt_dprojs,vectin_)
          call abi_gpu_xgemm(cplex, 'C', 'N', ngrads*nprojs, ndat*nspinor, npwin, cone, &
                 c_loc(current_ikpt_dprojs), npwin,&
                 c_loc(vectin_), npwin, czero, c_loc(dprojections), ngrads*nprojs)
          !$OMP END TARGET DATA
        end if
      else
        call gemm_nonlop_ompgpu_distributed_gemm_opernla(rank,nprocs,npwin,ndat,nspinor,&
        &                                         nprojs,&
        &                                         gemm_nonlop_kpt(ikpt)%nprojs_blk,&
        &                                         gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
        &                                         nprojs_my_blk,cplex,czero,&
        &                                         current_ikpt_projs,projs_recv,&
        &                                         vectin,projections_ptr)
        if(signs==1.and.ngrads>0) then
          call gemm_nonlop_ompgpu_distributed_gemm_opernla(rank,nprocs,npwin,ndat,nspinor,&
          &                                         ngrads*nprojs,&
          &                                         ngrads*gemm_nonlop_kpt(ikpt)%nprojs_blk,&
          &                                         ngrads*gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
          &                                         ngrads*nprojs_my_blk,cplex,czero,&
          &                                         current_ikpt_dprojs,dprojs_recv,&
          &                                         vectin,dprojections)
        end if
      end if
    else
      ! only compute real part of projections = P^* psi => projections_r = P_r^T psi_r + P_i^T psi_i
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !!$OMP TARGET LOOP &
      !$OMP& MAP(to:temp_realvec_r,vectin) PRIVATE(i)
      do i=1, npwin*nspinor*ndat
        temp_realvec_r(i) = vectin(1,i)
      end do
      if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
        !!$OMP TARGET LOOP &
        !$OMP& MAP(to:temp_realvec_r) PRIVATE(idat)
        do idat=1, ndat*nspinor
          temp_realvec_r(1+(idat-1)*npwin) = temp_realvec_r(1+(idat-1)*npwin)/2
        end do
      end if
      !$OMP TARGET DATA USE_DEVICE_PTR(projections_ptr,current_ikpt_projs_r,temp_realvec_r)
      call abi_gpu_xgemm(cplex, 'T', 'N', nprojs, ndat*nspinor, npwin, cone, &
&                c_loc(current_ikpt_projs_r), npwin, &
&                c_loc(temp_realvec_r), npwin, czero, c_loc(projections_ptr), nprojs)
      !$OMP END TARGET DATA
      if(signs==1.and.ngrads>0) then
        !$OMP TARGET DATA USE_DEVICE_PTR(dprojections,current_ikpt_dprojs_r,temp_realvec_r)
        call abi_gpu_xgemm(cplex, 'T', 'N', ngrads*nprojs, ndat*nspinor, npwin, cone, &
&                  c_loc(current_ikpt_dprojs_r), npwin, &
&                  c_loc(temp_realvec_r), npwin, czero, c_loc(dprojections), ngrads*nprojs)
        !$OMP END TARGET DATA
      end if

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
      !!$OMP TARGET LOOP &
      !$OMP& MAP(to:temp_realvec_i,vectin) PRIVATE(i)
      do i=1, npwin*nspinor*ndat
        temp_realvec_i(i) = vectin(2,i)
      end do
      if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
        !!$OMP TARGET LOOP &
        !$OMP& MAP(to:temp_realvec_i) PRIVATE(idat)
        do idat=1, ndat*nspinor
          temp_realvec_i(1+(idat-1)*npwin) = zero
        end do
      end if

      !$OMP TARGET DATA USE_DEVICE_PTR(projections_ptr,current_ikpt_projs_i,temp_realvec_i)
      call abi_gpu_xgemm(cplex, 'T', 'N', nprojs, ndat*nspinor, npwin, cone, &
               c_loc(current_ikpt_projs_i), npwin, &
               c_loc(temp_realvec_i), npwin, cone , c_loc(projections_ptr), nprojs)
      !$OMP END TARGET DATA

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
      !!$OMP TARGET LOOP &
      !$OMP& MAP(to:projections_ptr) PRIVATE(i1,i2)
      do i2=1, nspinor*ndat
        do i1=1, nprojs
          projections_ptr(1,i1,i2) = projections_ptr(1,i1,i2) * 2
        end do
      end do

      if(signs==1.and.ngrads>0) then
        !$OMP TARGET DATA USE_DEVICE_PTR(dprojections,current_ikpt_dprojs_i,temp_realvec_i)
        call abi_gpu_xgemm(cplex, 'T', 'N', ngrads*nprojs, ndat*nspinor, npwin, cone, &
                  c_loc(current_ikpt_dprojs_i), npwin, &
                  c_loc(temp_realvec_i), npwin, cone , c_loc(dprojections), ngrads*nprojs)
        !$OMP END TARGET DATA

        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
        !!$OMP TARGET LOOP &
        !$OMP& MAP(to:dprojections) PRIVATE(i1,i2)
        do i2=1, nspinor*ndat
          do i1=1, ngrads*nprojs
            dprojections(1,i1,i2) = dprojections(1,i1,i2) * 2
          end do
        end do
      end if
    end if

    if(cpopt >= 0) then
      ! store in cprjin
      if(.not. local_vectproj .and. cpopt/=3) then
        !TODO This use-case is extremely painful for GEMM OpenGPU nonlop performance
        !$OMP TARGET UPDATE FROM(projections_ptr)
        !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,nlmn)
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections_ptr(1:cplex, shift+1:shift+nlmn, idat)
            shift = shift + nlmn
          end do
        end do
      end if
      if(cpopt==3) then
        ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (2)")
        !$OMP TARGET UPDATE FROM(dprojections)
        !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,igrad,nlmn)
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            do igrad=1,ngrads
              cprjin(iatom, idat)%dcp(1:cplex,igrad,1:nlmn) = &
  &              dprojections(1:cplex, shift+1:shift+nlmn, idat)
              shift = shift + nlmn
            end do
          end do
        end do
      end if
    end if
  end if



  if(choice > 0) then

    if(choice /= 7) then
      ! opernlc
      iatm = 0
      ndgxdt = 0
      ndgxdtfac = 0
      nd2gxdt = 0
      nd2gxdtfac = 0
      optder = 0

      shift = 0
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        call opernlc_ylm_ompgpu(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
&         dgxdt_dum_in,dgxdt_dum_out,dgxdt_dum_out2,&
&         d2gxdt_dum_in,d2gxdt_dum_out,d2gxdt_dum_out2,dimenl1,dimenl2,dimekbq,enl,&
&         projections_ptr,&
&         vnl_projections,&
&         s_projections,&
&         iatm,indlmn,itypat,lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&         nattyp(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ,ndat,ibeg-1,iend,nprojs,ntypat)

        shift = shift + nattyp(itypat)*nlmn
        iatm = iatm+nattyp(itypat)
      end do
    else
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
      !!$OMP TARGET LOOP &
      !$OMP& MAP(to:projections_ptr,s_projections) PRIVATE(i1,i2)
      do i2=1, nspinor*ndat
        do i1=1, nprojs
          s_projections(1,i1,i2) = projections_ptr(1,i1,i2)
          s_projections(2,i1,i2) = projections_ptr(2,i1,i2)
        end do
      end do
    end if

    ! opernlb (only choice=1)
    if(signs==2) then
      if(paw_opt == 3 .or. paw_opt == 4) then
        ! Get svectout from s_projections
        if(cplex == 2) then
          if(.not. gemm_nonlop_is_distributed) then
            !$OMP TARGET DATA USE_DEVICE_PTR(s_projections,current_ikpt_projs,svectout_)
            call abi_gpu_xgemm(cplex, 'N', 'N', npwout, ndat*nspinor, nprojs, cone, &
                    c_loc(current_ikpt_projs), npwout,&
                    c_loc(s_projections), nprojs, czero, c_loc(svectout_), npwout)
            !$OMP END TARGET DATA
          else
            call gemm_nonlop_ompgpu_distributed_gemm_opernlb(rank,nprocs,npwout,ndat,nspinor,&
            &                                         nprojs,&
            &                                         gemm_nonlop_kpt(ikpt)%nprojs_blk,&
            &                                         gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
            &                                         nprojs_my_blk,cplex,&
            &                                         current_ikpt_projs,projs_recv,&
            &                                         s_projections,svectout)
          end if
        else
          !$OMP TARGET DATA USE_DEVICE_PTR(s_projections,current_ikpt_projs_r,temp_realvec_r)
          call abi_gpu_xgemm(cplex, 'N', 'N', npwout, ndat*nspinor, nprojs, cone, &
                    c_loc(current_ikpt_projs_r), npwout, &
                    c_loc(s_projections), nprojs, czero, c_loc(temp_realvec_r), npwout)
          !$OMP END TARGET DATA
          !$OMP TARGET DATA USE_DEVICE_PTR(s_projections,current_ikpt_projs_i,temp_realvec_i)
          call abi_gpu_xgemm(cplex, 'N', 'N', npwout, ndat*nspinor, nprojs, cone, &
                    c_loc(current_ikpt_projs_i), npwout,&
                    c_loc(s_projections), nprojs, czero, c_loc(temp_realvec_i), npwout)
          !$OMP END TARGET DATA

          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
          !!$OMP TARGET LOOP &
          !$OMP& MAP(to:temp_realvec_r,temp_realvec_i,svectout) PRIVATE(i)
          do i=1, npwin*nspinor*ndat
            svectout(1,i) = temp_realvec_r(i)
            svectout(2,i) = temp_realvec_i(i)
          end do
        end if
        if(choice /= 7) then
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
          !!$OMP TARGET LOOP &
          !$OMP& MAP(to:vectin,svectout) PRIVATE(i)
          do i=1, npwin*nspinor*ndat
            svectout(1,i) = svectout(1,i) + vectin(1,i)
            svectout(2,i) = svectout(2,i) + vectin(2,i)
          end do
        end if
      end if
      if(paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then
        ! Get vectout from vnl_projections
        if(cplex_fac == 2) then
          if(.not. gemm_nonlop_is_distributed) then
            !$OMP TARGET DATA USE_DEVICE_PTR(vnl_projections,current_ikpt_projs,vectout_)
            call abi_gpu_xgemm(cplex, 'N', 'N', npwout, ndat*nspinor, nprojs, cone, &
                    c_loc(current_ikpt_projs), npwout, &
                    c_loc(vnl_projections), nprojs, czero, c_loc(vectout_), npwout)
            !$OMP END TARGET DATA
          else
            call gemm_nonlop_ompgpu_distributed_gemm_opernlb(rank,nprocs,npwout,ndat,nspinor,&
            &                                         nprojs,&
            &                                         gemm_nonlop_kpt(ikpt)%nprojs_blk,&
            &                                         gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
            &                                         nprojs_my_blk,cplex,&
            &                                         current_ikpt_projs,projs_recv,&
            &                                         vnl_projections,vectout)
          end if
        else
          !$OMP TARGET DATA USE_DEVICE_PTR(vnl_projections,current_ikpt_projs_r,temp_realvec_r)
          call abi_gpu_xgemm(cplex, 'N', 'N', npwout, ndat*nspinor, nprojs, cone, &
                  c_loc(current_ikpt_projs_r), npwout, &
                  c_loc(vnl_projections), nprojs, czero, c_loc(temp_realvec_r), npwout)
          !$OMP END TARGET DATA
          !$OMP TARGET DATA USE_DEVICE_PTR(vnl_projections,current_ikpt_projs_i,temp_realvec_i)
          call abi_gpu_xgemm(cplex, 'N', 'N', npwout, ndat*nspinor, nprojs, cone, &
                  c_loc(current_ikpt_projs_i), npwout, &
                  c_loc(vnl_projections), nprojs, czero, c_loc(temp_realvec_i), npwout)
          !$OMP END TARGET DATA

          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
          !!$OMP TARGET LOOP &
          !$OMP& MAP(to:temp_realvec_r,temp_realvec_i,vectout) PRIVATE(i)
          do i=1, npwin*nspinor*ndat
            vectout(1,i) = temp_realvec_r(i)
            vectout(2,i) = temp_realvec_i(i)
          end do
        end if
      end if
    end if ! opernlb

    ! opernld
    if(signs==1) then
#ifdef HAVE_GPU_HIP
      !$OMP TARGET UPDATE FROM(vnl_projections,projections_ptr,dprojections)
#endif
      if(choice==1.or.choice==3.or.choice==23) then
        shift=0
        iatm=0
        esum=zero
        do itypat=1, ntypat
          nlmn=count(indlmn(3,:,itypat)>0)
          ibeg = shift+1
          iend = shift+nattyp(itypat)*nlmn
          nattyp_i = nattyp(itypat)
#ifndef HAVE_GPU_HIP
          !$OMP TARGET TEAMS DISTRIBUTE &
          !$OMP& MAP(to:vnl_projections,projections_ptr,enlk) &
          !$OMP& FIRSTPRIVATE(idat,itypat,nlmn,esum)
#endif
          do idat=1,ndat*nspinor
            esum=zero
#ifndef HAVE_GPU_HIP
            !$OMP PARALLEL DO COLLAPSE(3) REDUCTION(+:esum) &
            !$OMP& PRIVATE(ia,ilmn,ii)
#endif
            do ia=1,nattyp_i
              do ilmn=1,nlmn
                do ii=1,cplex
                  esum=esum +vnl_projections(ii,shift+(ia-1)*nlmn+ilmn,idat) &
&                           *projections_ptr    (ii,shift+(ia-1)*nlmn+ilmn,idat)
                end do
              end do
            end do
            enlk(idat) = enlk(idat) + esum
          end do
          shift = shift + nattyp(itypat)*nlmn
          iatm = iatm+nattyp(itypat)
        end do
        if (choice==1) then
#ifndef HAVE_GPU_HIP
          !$OMP TARGET PARALLEL DO MAP(to:enlout,enlk) PRIVATE(idat)
#endif
          do idat=1,ndat
            enlout(idat)=enlk(idat)
          end do
        end if
      end if ! choice=1/3/23
      if(choice==2.or.choice==3.or.choice==23) then
        grad_shift=merge(9,6,choice==23)
        if (choice==3.or.choice==23) then
          shift=0
          iatm=0
          do itypat=1, ntypat
            nlmn=count(indlmn(3,:,itypat)>0)
            ibeg = shift+1
            iend = shift+nattyp(itypat)*nlmn
            nattyp_i = nattyp(itypat)

#ifndef HAVE_GPU_HIP
            !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
            !$OMP& MAP(to:vnl_projections,dprojections,enlout) &
            !$OMP& FIRSTPRIVATE(idat,itypat,nlmn,esum)
#endif
            do idat=1,ndat*nspinor
              do igrad=1,6
                esum=zero
#ifndef HAVE_GPU_HIP
                !$OMP PARALLEL DO COLLAPSE(3) REDUCTION(+:esum) &
                !$OMP& PRIVATE(ia,ilmn,ii)
#endif
                do ia=1,nattyp_i
                  !Following loops are a [D][Z]DOT
                  do ilmn=1,nlmn
                    do ii=1,cplex
                      esum=esum +vnl_projections(ii,shift+(ia-1)*nlmn+ilmn,idat) &
&                               *dprojections   (ii,grad_shift*shift + (ia-1)*nlmn*grad_shift + (igrad-1)*nlmn +ilmn,idat)
                    end do
                  end do
                end do
                enlout((idat-1)*nnlout+igrad) = enlout((idat-1)*nnlout+igrad) + two*esum
              end do
            end do

            shift = shift + nattyp(itypat)*nlmn
            iatm = iatm+nattyp(itypat)
          end do
        end if

        if (choice==2.or.choice==23) then
          shift=0
          iatm=0
          force_shift=merge(6,0,choice==23)
          do itypat=1, ntypat
            nlmn=count(indlmn(3,:,itypat)>0)
            nattyp_i = nattyp(itypat)

#ifndef HAVE_GPU_HIP
            !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(3) &
            !$OMP& MAP(to:vnl_projections,dprojections,enlout) &
            !$OMP& FIRSTPRIVATE(idat,itypat,nlmn,esum)
#endif
            do idat=1,ndat*nspinor
              do ia=1,nattyp_i
                do igrad=1,3
                  !Following loops are a [D][Z]DOT
                  esum=zero
#ifndef HAVE_GPU_HIP
                  !$OMP PARALLEL DO COLLAPSE(2) REDUCTION(+:esum) &
                  !$OMP& PRIVATE(ilmn,ii)
#endif
                  do ilmn=1,nlmn
                    do ii=1,cplex
                      esum=esum +vnl_projections(ii,shift+(ia-1)*nlmn+ilmn,idat) &
&                               *dprojections(ii,grad_shift*shift+(ia-1)*nlmn*grad_shift+(igrad-1+force_shift)*nlmn +ilmn,idat)
                    end do
                  end do
                  enlout((idat-1)*nnlout + force_shift + (iatm+ia-1)*3 + igrad)= &
&                               enlout((idat-1)*nnlout + force_shift + (iatm+ia-1)*3 + igrad) + two*esum
                end do
              end do
            end do

            shift = shift + nattyp(itypat)*nlmn
            iatm = iatm+nattyp(itypat)
          end do
        end if
      end if ! choice=2, 3 or 23
#ifndef HAVE_GPU_HIP
      !$OMP TARGET UPDATE FROM(enlout,enlk)
#endif
    end if !opernld

  end if ! choice>0

  !$OMP TASKWAIT
! Reduction in case of parallelism
  if (signs==1.and.mpi_enreg%paral_spinor==1) then
    if (size(enlout)>0) then
      call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
    end if
    if (choice==3.or.choice==23) then
      call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
    end if
  end if

! Derivatives wrt strain
!  - Convert from reduced to cartesian coordinates
!  - Substract volume contribution
 if ((choice==3.or.choice==23).and.signs==1.and.paw_opt<=3) then
   do idat=1,ndat
     enlout_shift=(idat-1)*nnlout
     call strconv(enlout(enlout_shift+1:enlout_shift+6),gprimd,work)
     enlout(enlout_shift+1:enlout_shift+3)=(work(1:3)-enlk(idat))
     enlout(enlout_shift+4:enlout_shift+6)= work(4:6)
   end do
 end if

 ! Retrieve and release allocated buffers
 !$OMP TARGET EXIT DATA MAP(delete:vectin) IF(transfer_vectin)
 !$OMP TARGET EXIT DATA MAP(from:vectout)   IF(transfer_vectout)
 !$OMP TARGET EXIT DATA MAP(from:svectout)  IF(transfer_svectout)

#ifdef HAVE_GPU_CUDA
 !$OMP TARGET EXIT DATA MAP(delete:s_projections,vnl_projections)
 if(.not. local_vectproj) then
   !$OMP TARGET EXIT DATA MAP(delete:projections_ptr)
 end if
#endif

! Release memory
  if(signs == 1) then
    !$OMP TARGET EXIT DATA MAP(delete:enlout)
  end if
  if(gemm_nonlop_is_distributed) then
    !$OMP TARGET EXIT DATA MAP(delete:projs_recv)
    ABI_FREE(projs_recv)
    if (signs==1.and.ngrads>0) then
      !$OMP TARGET EXIT DATA MAP(delete:dprojs_recv)
      ABI_FREE(dprojs_recv)
    end if
  end if
  !$OMP TARGET EXIT DATA MAP(delete:atindx1,indlmn,enl)
  if (allocated(dprojections)) then
    !$OMP TARGET EXIT DATA MAP(delete:dprojections)
    ABI_FREE(dprojections)
  end if
  if (allocated(enlk)) then
    !$OMP TARGET EXIT DATA MAP(delete:enlk,nattyp)
    ABI_FREE(enlk)
  end if

 end subroutine gemm_nonlop_ompgpu

!***


#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! stubs for compiling with OpenMP GPU offload disabled.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine init_gemm_nonlop_ompgpu(nkpt)
  integer,intent(in) :: nkpt

  ABI_UNUSED((/nkpt/))
  ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
 end subroutine init_gemm_nonlop_ompgpu
!!***

 subroutine destroy_gemm_nonlop_ompgpu()
  ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")
 end subroutine destroy_gemm_nonlop_ompgpu
!!***

 subroutine make_gemm_nonlop_ompgpu(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k, &
&                            ph3d_k,kpt_k,kg_k,kpg_k, &
&                            compute_grad_strain,compute_grad_atom) ! Optional parameters

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax,ntypat
  integer, intent(in) :: indlmn(:,:,:), kg_k(:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  logical, intent(in), optional :: compute_grad_strain,compute_grad_atom
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)
  real(dp), intent(in) :: kpt_k(:)
  real(dp), intent(in), target :: kpg_k(:,:)

  ABI_UNUSED((/ikpt,npw,lmnmax,ntypat,indlmn,kg_k,nattyp,istwf_k/))
  ABI_UNUSED((/ucvol,ffnl_k,ph3d_k,kpt_k,kpg_k/))
  ABI_UNUSED((/compute_grad_strain,compute_grad_atom/))
  ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")

 end subroutine make_gemm_nonlop_ompgpu
!!***

 subroutine gemm_nonlop_ompgpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm
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
  real(dp),intent(out),target :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
  real(dp),intent(inout),target :: vectout(2,npwout*nspinor*ndat) !vz_i
  real(dp),intent(inout),optional,target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ABI_UNUSED((/choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir/))
  ABI_UNUSED((/istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin/))
  ABI_UNUSED((/nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO/))
  ABI_UNUSED((/paw_opt,signs,tim_nonlop,useylm,gpu_option/))
  ABI_UNUSED((/atindx1,indlmn,kgin,kgout,nattyp,ngfft,nloalg/))
  ABI_UNUSED((/enl,ffnlin,ffnlout,gmet,gprimd,kpgin,kpgout,kptin,kptout,phkxredin,phkxredout/))
  ABI_UNUSED((/ucvol,lambda,sij,ph1d(1,1),ph3din,ph3dout,vectin,enlout,svectout,vectout,vectproj/))
  ABI_UNUSED_A(cprjin)
  ABI_UNUSED_A(mpi_enreg)
  ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")

 end subroutine gemm_nonlop_ompgpu

#endif

end module m_gemm_nonlop_ompgpu
!!***

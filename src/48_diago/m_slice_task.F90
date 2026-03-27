!!****f* ABINIT/m_slice_task
!! NAME
!! m_slice_task
!!
!! FUNCTION
!! This module contains types and routines to implement scheduler and resource allocator 
!! for slice tasks. Implements scheduler and resource allocator for slice tasks.
!! Implements logic for asynchronous slice memory avoiding race condition in read/write.
!!
!! NOTES
!! The logic for parallel slice treatment (a task) is the following.
!! Memory and workload are distributed using a 2D cartesian grid. Let's assume 
!! for simplicity that we have four MPI processes in the spacecom communicator. 
!! Matrix X is distributed along plane-waves at the beginning:
!!
!!                    bands
!!            |-------------------|
!!            |        P0         |
!!            |                   |
!!            |-------------------|
!!            |        P1         |
!!            |                   |
!!        pw  |-------------------|
!!            |        P2         |
!!            |                   |
!!            |-------------------|
!!            |        P3         |
!!            |                   |
!!            |-------------------|
!!
!! At the start, we use xgTransposer to MPI transpose the matrix X
!! achieving a custom layout for bandpp, and we end up with:
!! 
!!                    bands
!!            |-------|---|---|---|
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!        pw  |  P0   |P1 |P2 |P3 |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |       |   |   |   |
!!            |-------|---|---|---|
!!
!! From there, we can define slices acting on subgroup of processes.
!! For example, slice one can have process 0 and slice two the remaining 
!! 1,2,3 processes. MPI transposing to Linalg representation using
!! the slice sub-communicators yields:
!!
!!                    bands
!!            |-------|-----------|
!!            |       |           |
!!            |       |    P1     |
!!            |       |           |
!!            |       |-----------|
!!            |       |           |
!!        pw  |  P0   |    P2     |
!!            |       |           |
!!            |       |-----------|
!!            |       |           |
!!            |       |    P3     |
!!            |       |           |
!!            |-------|-----------|
!!
!! At this point, parallel Rayleigh-Ritz is possible.
!
!!
!! COPYRIGHT
!! Copyright (C) 2018-2026 ABINIT group (IML)
!! This file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

module m_slice_task

    use defs_basis
    use defs_abitypes
    use m_abicore
    use m_errors
    use m_time, only : timab
    use m_sort, only: sort_dp

    use m_cgtools
    use m_xg
    use m_xgTransposer
    use m_xg_ortho_RR

    use m_xmpi
    use m_xomp
#ifdef HAVE_OPENMP
    use omp_lib
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize
#endif

#if defined(HAVE_GPU_MARKERS)
    use m_nvtx_data
#endif

    implicit none

    private

    ! Timers
    !---------------------------------------------------
    integer, parameter :: tim_swap        = 1761
    integer, parameter :: tim_RR_q        = 1759
    integer, parameter :: tim_barrier     = 1764
    integer, parameter :: tim_copy        = 1765
    integer, parameter :: tim_Bortho_X    = 1641
    integer, parameter :: tim_getAX_BX    = 1754
    integer, parameter :: tim_invovl      = 1755

    ! Public 'matrixInfo' datatype
    !-------------------------------------------------
    type, public :: matrixInfo_t

        integer :: comm_rows                    ! xmpi_comm_self ...
        integer :: comm_cols                    ! same as spacecom
        integer :: spacecom                     ! same as comm_cols
        integer :: neigenpairs                  ! total number of bands (=number of eigenpairs)
        integer :: total_spacedim               ! total number of plane-waves
        integer :: spacedim                     ! nb of plane-waves per process in linalg representation
        integer :: space                        ! real or complex eigenvectors
        integer :: gpu_kokkos_nthrd                 
        integer :: gpu_thread_limit 
        integer :: gpu_option                   ! enable GPU
        integer :: paral_kgb                    ! enable parallel (k-points, G basis, bands)
        integer :: me_g0
        integer :: me_g0_fft

    end type matrixInfo_t

    ! Public 'activeTask' datatype for active slice in use (me=current MPI process)
    !-------------------------------------------------
    type, public :: activeTask_t

        ! MPI-related information for active task
        integer :: me_g0                                    ! process contains G(0,0,0) 
        integer :: me_g0_fft                                ! process contains G(0,0,0) for fft
        integer :: me_nproc_slice                           ! number of processes reserved to slice in use
        integer :: me_comm_slice                            ! (sub-)communicator reserved to slice in use
        integer :: me_comm_rows                             ! (sub-)communicator for rows in Transposer
        integer :: me_comm_cols                             ! (sub-)communicator for cols in Transposer
        integer :: me_id_slice                              ! identifier from 1 to nslice of slice task in use
        integer :: me_neigenpairs_slice                     ! total number of eigenpairs of slice in use
        integer :: me_bandpp_slice                          ! number of distributed bands per process for slice in use
        integer :: me_lowb
        integer :: me_uppb
        integer :: me_ndeg_slice                            ! polynomial filter degree for slice in use

        ! flags for slice location
        logical :: is_lowpass
        logical :: is_last

        integer, allocatable :: me_cols_X(:)                ! columns of shared memory used in task
        integer, allocatable :: me_cols_Xext(:)             ! columns of asynchronous memory used in task
        integer, allocatable :: me_mask_Xext(:)             ! converged columns

        ! MPI column and row distribution for active task
        integer, allocatable :: me_ncolsColsRows_slice(:)   ! ncol of colsrows representation for my slice
                                                            !   (column capacity per process in colsrows repr)
        integer, allocatable :: me_nrowsLinalg_slice(:)     ! nrow of linalg representation for my slice
                                                            !   (row capacity per process in linalg repr)
        
        ! Memory address used for reads/writes
        type(xgBlock_t) :: me_Xext                          ! eigenvector memory in use by active slice
        type(xgBlock_t) :: me_resid
        type(xgBlock_t) :: me_eigen

        ! Memory space used for computation
        !type(chebfi_t) :: chebfi ! should not use chebfi objects

    end type activeTask_t

    ! Public 'taskScheduler' datatype for asynchronous slice treatment
    ! it is basically for asynchronous treatment of the 'async' memory
    !-------------------------------------------------
    type, public :: taskScheduler_t

        integer :: next_task_id
        integer :: ntasks
        integer :: nresources

        integer, allocatable :: ncols_per_task(:)       ! number of bands per slice
        integer, allocatable :: ncols_per_proc(:)       ! number of bands per MPI rank
        integer, allocatable :: nresources_per_task(:)  ! number of available processes per slice
        integer, allocatable :: lookup_proc(:)          ! which slice each MPI rank serves

    end type taskScheduler_t

    ! Public 'asyncMemory' datatype for slice input/output without race condition
    ! is basically created by blocks of number of asynchronous tasks
    !-------------------------------------------------
    type, public :: asyncMemory_t

        integer :: neigenpairs_ext                          ! total number of extended columns
        integer :: paral_kgb 
        logical :: has_transposer

        ! Memory buffers
        type(xg_t) :: X_ext                                 ! eigenvector memory used by all slices for I/O
        type(xg_t) :: eigen_ext                             ! eigenvalue memory used by all slices for I/O
        type(xg_t) :: resid_ext                             ! residual memory used by all slices for I/O
        type(xgTransposer_t) :: xgTransposerXext            ! transposer datastructure for eigenvectors
        
        ! Pointers
        type(xgBlock_t) :: XextLinalg

        ! Arrays related to MPI (all-slice phase)
        integer, allocatable :: ncolsColsRows(:)            ! ncol of colsrows representation for global X

        ! Arrays related to data distribution (all-slice phase)
        integer, allocatable :: lookup_cols_Xext(:)         ! which slice each column serves in extented memory

        ! Flags for MPI distribution state
        logical :: use_linalg = .false.                     ! use linalg representation
        logical :: use_colsrows = .false.                   ! use colsrows representation

    end type asyncMemory_t
 
    ! Public methods
    !-------------------------------------------------
    public :: init_matrixInfo                   ! wrapper for various xgBlock parameters
    public :: slice_task_allocateAsyncMemory    ! allocates async memory buffer
    public :: slice_task_initAsyncMemory        ! fills async memory colwise
    public :: slice_task_initSchedule           ! compute 'process-to-slices' distribution
    public :: free_extended_memory
    public :: init_active_memory
    public :: mark_active_task
    public :: allocate_active_task
    public :: free_active_task
    public :: execute_active_task
    public :: mask_active_task
    public :: compress_extended_memory

    CONTAINS  
!=====================================================================
!!***

!!****f* m_slice_task/init_matrixInfo
!! NAME
!! init_matrixInfo
!! 
!! SOURCE

  subroutine init_matrixInfo(matrixInfo, comm_rows, comm_cols, spacecom, neigenpairs, total_spacedim, &
          spacedim, space, gpu_kokkos_nthrd, gpu_thread_limit, gpu_option, paral_kgb, me_g0, me_g0_fft)

      implicit none

      type(matrixInfo_t), intent(inout) :: matrixInfo
      integer, intent(in) :: comm_rows, comm_cols, spacecom, neigenpairs, total_spacedim, spacedim, space
      integer, intent(in) :: gpu_kokkos_nthrd, gpu_thread_limit, gpu_option, paral_kgb, me_g0, me_g0_fft

      matrixInfo%comm_rows      = comm_rows
      matrixInfo%comm_cols      = comm_cols
      matrixInfo%spacecom       = spacecom
      matrixInfo%neigenpairs    = neigenpairs
      matrixInfo%total_spacedim = total_spacedim
      matrixInfo%spacedim       = spacedim
      matrixInfo%space          = space
      matrixInfo%gpu_kokkos_nthrd = gpu_kokkos_nthrd
      matrixInfo%gpu_thread_limit = gpu_thread_limit
      matrixInfo%gpu_option     = gpu_option
      matrixInfo%paral_kgb      = paral_kgb
      matrixInfo%me_g0          = me_g0
      matrixInfo%me_g0_fft      = me_g0_fft

  end subroutine init_matrixInfo
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/slice_task_initSchedule
!! NAME
!! slice_task_initSchedule
!!
!! FUNCTION
!! Initialization of scheduler object from 'npro'c available resources
!! using execution options given by 'paral_kgb' and 'paral_task'.
!! Essentially allows to compute and apply (via wrapper to xgTransposer) 
!! an intermediate level of MPI distribution along tasks, on top of paral_kgb level. 
!! - If enable_paral then we should first divide processes to slices 
!! - If disable_paral then we should use all available processes for every slice (no division)
!! 
!! SOURCE

  subroutine slice_task_initSchedule(scheduler, nproc, paral_kgb, paral_task)
      implicit none
      type(taskScheduler_t), intent(inout) :: scheduler
      integer, intent(in) :: nproc, paral_kgb, paral_task
      scheduler%ntasks = nslice
      scheduler%nresources = nproc
      scheduler%next_task_id = 0
      call alloc_schedule(scheduler)
      if (paral_kgb==0 .or. paral_task==0) then
          ! do not divide available resources directly say that all processes are
          ! available per slice
      else 
          ! compute
      end if
  end subroutine slice_task_initSchedule
!!***

!! Constructor for scheduler object
subroutine alloc_schedule(scheduler)
    implicit none
    type(taskScheduler_t), intent(inout) :: scheduler
    call free_schedule(scheduler)
    ABI_MALLOC_IFNOT(scheduler%neigenpairs_per_slice, (scheduler%ntasks))
    ABI_MALLOC_IFNOT(scheduler%nproc_per_slice, (scheduler%ntasks))
    ABI_MALLOC_IFNOT(scheduler%lookup_proc, (scheduler%nresources))
end subroutine alloc_schedule

!! Destructor for scheduler object
subroutine free_schedule(scheduler)
    implicit none
    type(taskScheduler_t), intent(inout) :: scheduler
    ABI_SFREE(scheduler%neigenpairs_per_slice)
    ABI_SFREE(scheduler%nproc_per_slice)
    ABI_SFREE(scheduler%lookup_proc)
end subroutine free_schedule

!! Logic for parallel execution of tasks (a slice has SOME MPI processes)
subroutine schedule_parallel_tasks(scheduler)
    implicit none
    type(taskScheduler_t), intent(inout) :: scheduler
    integer :: ntasks, nproc
    integer, allocatable :: weights(:)
    nproc = scheduler%nresources
    ntasks = scheduler%ntasks 
    !ABI_MALLOC_IFNOT(weights, (slice%ntasks))
    weights = 1 
    ! weights = slice%poly_degrees
    !call fair_allocation(slice%ntasks, slice%neigenpairs_per_slice, weights, nproc, slice%nproc_per_slice)
    !call assign_tasks_to_processes(slice%nproc_per_slice, slice%lookup_proc)
    ABI_SFREE(weights)
end subroutine schedule_parallel_tasks

!! Logic for sequential execution of tasks (a slice has ALL MPI processes)
subroutine schedule_next_task(scheduler, load_size)
    implicit none
    type(taskScheduler_t), intent(inout) :: scheduler
    integer, intent(in) :: load_size
    scheduler%next_task_id = scheduler%next_task_id + 1
    scheduler%neigenpairs_per_slice = load_size
    scheduler%nproc_per_slice = scheduler%nresources !! use all MPI, can be 1
    scheduler%lookup_proc = scheduler%next_task_id
end subroutine schedule_next_task

!----------------------------------------------------------------------

!!****f* m_slice_task/slice_task_allocateAsyncMemory
!! NAME
!! slice_task_allocateAsyncMemory
!! 
!! FUNCTION
!! Allocate async memory buffers and distribute them according
!! to slice logic. Essentially allocates memory work%Xext
!! favoring data overlap over communication overlap.
!! 
!! SOURCE

subroutine slice_task_allocateAsyncMemory(work, minfo, ncol_ext) 

    implicit none
    type(asyncMemory_t), intent(inout) :: work
    type(matrixInfo_t), intent(inout) :: minfo
    integer, intent(in) :: ncol_ext
    integer :: gpu_option

    work%paral_kgb = minfo%paral_kgb
    work%neigenpairs_ext = ncol_ext
    work%use_linalg = .true.  ! for sanity
    work%use_colsrows = .false. ! for sanity
    work%has_transposer = .false.
    gpu_option = minfo%gpu_option

    ABI_MALLOC_IFNOT(work%lookup_cols_Xext, (ncol_ext))
    
    ! Allocate extended space in linalg representation
    call xg_init(work%X_ext, minfo%space, minfo%spacedim, work%neigenpairs_ext, &
        minfo%spacecom, me_g0=minfo%me_g0, gpu_option=gpu_option)

    call xg_init(work%resid_ext, SPACE_R, work%neigenpairs_ext, 1, gpu_option=gpu_option)
    call xg_init(work%eigen_ext, SPACE_R, work%neigenpairs_ext, 1, gpu_option=gpu_option)

end subroutine slice_task_allocateAsyncMemory
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/slice_task_initAsyncMemory
!! NAME
!! slice_task_initAsyncMemory
!! 
!! FUNCTION
!! Initialize async memory (for all tasks)
!! IML dev note: 
!! current version initializes using sketching of wanted size. Might also need
!! to test if choosing directly random vectors in better.
!! 
!! SOURCE

subroutine slice_task_initAsyncMemory(work, X0, minfo, mapper, ncols_per_task)

    implicit none
    
    !Arguments ------------------------------------
    type(asyncMemory_t), intent(inout) :: work
    type(matrixInfo_t), intent(inout) :: minfo
    type(xgBlock_t), intent(in) :: X0
    logical, pointer, intent(in) :: mapper(:,:)
    integer, intent(in) :: ncols_per_task(:)

    !Local variables-------------------------------
    integer :: nslice, islice, fcol, fcol_ext_prev, fcol_ext
    integer :: nrows, ncols, k_sketch, me_ncompl, nselect, fcol_compl
    integer :: space, gpu_option, spacecom, fcol_ext_sketch
    integer, allocatable :: nrand_per_task(:)
    integer, allocatable :: ncompl_per_task(:)
    ! typed
    type(xg_t) :: X0_compl, X_sketch
    type(xgBlock_t) :: col_in, col_out, Xext_last
    
    ! *********************************************************************

    ! Sanity check
    if ((.not. work%use_linalg) .or. work%use_colsrows) then
        ABI_ERROR("not in linalg representation")
    end if

    nrows = rows(X0)
    ncols = cols(X0)
    space = minfo%space
    gpu_option = minfo%gpu_option
    spacecom = minfo%spacecom

    work%XextLinalg = work%X_ext%self
    nslice = size(mapper, dim=2)
    ABI_MALLOC(nrand_per_task, (nslice))
    ABI_MALLOC(ncompl_per_task, (nslice))

    ! Copy X to XextLinalg to achieve contiguous column blocks
    fcol_ext_prev = 1
    fcol_ext = 0
    do islice=1, nslice
        nselect = 0
        do fcol=1, ncols
            if (.not. mapper(fcol, islice)) then
                exit
            end if
            nselect = nselect + 1
            fcol_ext = fcol_ext + nselect
            call xgBlock_setBlock(X0, col_in, nrows, 1, fcol=fcol)
            call xgBlock_setBlock(work%XextLinalg, col_out, nrows, 1, fcol=fcol_ext)

            ! Reminder: xgBlock_copy is always on CPU expect if both blocks are on GPU
            call xgBlock_copy(col_in, col_out)
        end do
        nrand_per_task(islice) = ncols_per_task(islice) - nselect
        ncompl_per_task(islice) = ncols - nselect
        write(std_out,*) 'selected', nselect, 'out of', ncols_per_task(islice), 'then rand is', &
            nrand_per_task(islice)
        fcol_ext = ncols_per_task(islice) ! jump to end of slice
        work%lookup_cols_Xext(fcol_ext_prev:fcol_ext) = islice
        fcol_ext_prev = fcol_ext
    end do

    ! recover the remaining columns not mapped to the slice and sketch them
    do islice=1, nslice
        k_sketch = nrand_per_task(islice) ! size of sketch
        me_ncompl = ncompl_per_task(islice) ! size of complement in X0
        if (k_sketch==0 .or. me_ncompl==0) then
            exit
        end if
        call xg_init(X0_compl, space, nrows, me_ncompl, spacecom, gpu_option=gpu_option)
        fcol_compl = 0
        fcol_ext_sketch = 1
        do fcol=1, ncols
            if (.not. mapper(fcol, islice)) then ! for remaining dimensions not in slice
                fcol_compl = fcol_compl + 1 
                call xgBlock_setBlock(X0, col_in, nrows, 1, fcol=fcol)
                call xgBlock_setBlock(X0_compl%self, col_out, nrows, 1, fcol=fcol_compl)
                call xgBlock_copy(col_in, col_out)
            end if
        end do
        ! Y = X * Omega where Omega sketch matrix to capture all directions at once (linalg distribution)
        call xg_init(X_sketch, space, nrows, k_sketch, spacecom, gpu_option=gpu_option)
        call xgBlock_randomSketching(X0_compl%self, X_sketch%self, k_sketch)
        call xgBlock_setBlock(work%XextLinalg, Xext_last, nrows, k_sketch, fcol=fcol_ext_sketch)
        call xgBlock_copy(X_sketch%self, Xext_last)
        call xg_free(X0_compl)
        call xg_free(X_sketch)
        fcol_ext_sketch = ncols_per_task(islice) ! jump to end of slice
    end do

    ABI_FREE(nrand_per_task)
    ABI_FREE(ncompl_per_task)

end subroutine slice_task_initAsyncMemory
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/free_extended_memory
!! NAME
!! free_extended_memory
!! 
!! SOURCE

subroutine free_extended_memory(work) 

    implicit none
    type(asyncMemory_t), intent(inout) :: work

    call xg_free(work%X_ext)
    call xg_free(work%eigen_ext)
    call xg_free(work%resid_ext)
    ABI_SFREE(work%lookup_cols_Xext)
    if (work%has_transposer) then
        call xgTransposer_free(work%xgTransposerXext)
    end if

end subroutine free_extended_memory
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/init_active_memory
!! NAME
!! init_active_memory
!!
!! FUNCTION
!! apply offset, in linalg MPI distribution
!! 
!! SOURCE

subroutine init_active_memory(work, task, X0, p)

    implicit none

    type(asyncMemory_t), intent(inout) :: work
    type(activeTask_t), intent(inout) :: task
    type(xgBlock_t), intent(in) :: X0
    integer, intent(in) :: p
    integer :: k_sketch, m, m_wanted
    type(xg_t) :: X_sketch
    type(xgBlock_t) :: Xext_last
    type(xg_t) :: X_compl
    integer :: space_, nrows, ncols_large, spacecom, gpu_option_

    space_ = space(task%me_Xext)
    nrows = rows(task%me_Xext)
    ncols_large = cols(X0)
    spacecom = comm(task%me_Xext)
    gpu_option_ = gpu_option(task%me_Xext)

    m = size(task%me_cols_X)
    m_wanted = task%me_neigenpairs_slice
    k_sketch = m_wanted - m + p

    call xgBlock_setBlock(task%me_Xext, Xext_last, nrows, k_sketch)
    
    !call xg_init(X_compl, ..)
    ! todo gather complement of selected indices...

    ! TODO 
    ! Est-ce que c'est mieux de prendre un melange aléatoire des autres directions ou 
    ! de prendre simplement des vecteurs aléatoires?
    ! V1 essayer avec des vecteurs aléatoires, c'est plus facile d'implémenter

    ! for remaining dimensions not in p
    ! Y = X * Omega where Omega sketch matrix to capture all directions at once (linalg distribution)
    call xg_init(X_sketch, space_, nrows, ncols_large, spacecom, gpu_option=gpu_option_)
    call xgBlock_randomSketching(X0, X_sketch%self, k_sketch)
    call xgBlock_copy(X_sketch%self, Xext_last)
    call xg_free(X_sketch)

end subroutine init_active_memory
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/mark_active_task
!! NAME
!! mark_active_task
!!
!! FUNCTION
!! Mark allocated portion as actively used by setting all me_* variables.
!! My process only marks resources its assigned slice has reserved.
!!
!! SOURCE

subroutine mark_active_task(task, scheduler, spacecom, paral_slice)

    implicit none

    ! Arguments
    type(activeTask_t), intent(inout) :: task
    type(taskScheduler_t), intent(inout) :: scheduler
    integer, intent(in) :: spacecom
    integer, intent(in) :: paral_slice

    ! Local variables
    integer :: my_rank, my_id, my_rank_sub, ierr

    ! *********************************************************************
 
    my_rank = xmpi_comm_rank(spacecom) 
    my_id = scheduler%lookup_proc(my_rank + 1) + 1
    task%me_id_slice = my_id

    task%me_neigenpairs_slice = scheduler%neigenpairs_per_slice(my_id)
    task%me_nproc_slice = scheduler%nproc_per_slice(my_id)
    !task%me_ndeg_slice = scheduler%poly_degrees(my_id)

    if (task%me_nproc_slice<scheduler%nresources) then
        ! can split communicator
    else
        ! use spacecom entire one
    end if

    ABI_MALLOC_IFNOT(task%me_ncolsColsRows_slice, (task%me_nproc_slice))
    ABI_MALLOC_IFNOT(task%me_nrowsLinalg_slice, (task%me_nproc_slice))

!    if (slice%paral_kgb==1) then
!        ! Compute column distribution across active resources
!        call distribute_vectors(slice%me_neigenpairs_slice, slice%me_nproc_slice, task%me_ncolsColsRows_slice)
!    
!        ! Compute row distribution across active resources
!        call distribute_vectors(slice%total_spacedim, slice%me_nproc_slice, task%me_nrowsLinalg_slice)
!    end if
!
!    ! Split global comm into disjoint sub-comms, only procs with the same color communicate
!    if (slice%paral_kgb==0) then ! todo also implement the case of sequential slices?
!        slice%me_comm_slice = slice%spacecom
!        slice%me_comm_rows = slice%comm_rows
!        slice%me_comm_cols = slice%comm_cols
!    else
!        call xmpi_comm_split(slice%spacecom, slice%me_id_slice, my_rank, slice%me_comm_slice, ierr)        
!        if ( ierr /= xmpi_success ) then
!            ABI_ERROR("Error while creating slice spacecom subcommunicator")
!        end if
!        call xmpi_comm_split(slice%comm_rows, slice%me_id_slice, my_rank, slice%me_comm_rows, ierr)          
!        if ( ierr /= xmpi_success ) then
!            ABI_ERROR("Error while creating slice row subcommunicator")
!        end if      
!        call xmpi_comm_split(slice%comm_cols, slice%me_id_slice, my_rank, slice%me_comm_cols, ierr)         
!        if ( ierr /= xmpi_success ) then
!            ABI_ERROR("Error while creating slice col subcommunicator")
!        end if 
!    end if
!    
!    ! process waits for others to create their subcommunicators before using its own
!    call xmpi_barrier(slice%spacecom) 
!
!    ! Concatenate slice%me_ncolsColsRows_slice into collective slice%ncolsColsRows
!    my_rank_sub = xmpi_comm_rank(slice%me_comm_slice)
!    slice%me_bandpp_slice = slice%bandpp
!    if (slice%paral_kgb==1) then
!        slice%me_bandpp_slice = slice%me_ncolsColsRows_slice(my_rank_sub + 1)
!    end if
!    call xmpi_allgather(slice%me_bandpp_slice, slice%ncolsColsRows, slice%spacecom, ierr)
!    if ( ierr /= xmpi_success ) then
!        ABI_ERROR("Error while gathering number of columns in colsrows for all slices")
!    end if

end subroutine mark_active_task
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/allocate_active_task
!! NAME
!! allocate_active_task
!! 
!! FUNCTION
!! Distributes extended columns across **all** MPI processes
!! After the transposition each process contains the correct
!! bandpp corresponding to the slice so that no additional communication
!! has to be performed in order to bring band slices to processes.
!! Input is Xext (linalg state) distributed across global processes.
!! Attention chebfi%X is distributed across slice processes =/= Xext per process.
!! Extracting chebfi%X in the slice distribution from Xext would require comms.
!!
!! OUTPUT
!! eigen=converged eigenvalue array of size (neigenpairs,1), first rows written only
!! residu=slice residual array of size (neigenpairs,1), first rows written only
!! work%XextLinalg= guess/converged eigenvectors for all slices
!! 
!! SOURCE

subroutine allocate_active_task(work, task)

    implicit none
    !type(slice_t), intent(inout) :: slice ! should not use slice objects at all!!!
    type(activeTask_t), intent(inout) :: task
    type(asyncMemory_t), intent(inout) :: work
    !type(chebfi_t), intent(inout) :: chebfi ! should not use chebfi objects at all !!!
    
    integer :: nbdbuf, oracle, num_proc
    real(dp) :: oracle_factor, oracle_min_occ

!    ! ========================== Transpose ===================================
!    !! Function
!    ! This transposition allows for a slice to not see others. 
!    ! It serves as a transition from global communicator to slice communicator.
!    ! 
!    
!    ! what are all those things comm_rows etc must use the ones from task TODO
!
!
!    ncolsColsRows_ptr => slice%ncolsColsRows ! todo this is not necessary but ok
!
!    if (slice%paral_kgb==1) then
!        ! Allocate slice%me_Xext according to the target MPI distribution for slices
!        call xgTransposer_constructor(slice%xgTransposerXext, work%XextLinalg, slice%me_Xext,&
!            nspinor, STATE_LINALG, TRANS_ALL2ALL, slice%comm_rows, slice%comm_cols, 0, 0, slice%me_g0_fft,&
!            gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit,&
!            custom_ncolsColsRows=.true., ncolsColsRows_sub=ncolsColsRows_ptr)
!   
!        slice%xgTransposerXext%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
!    
!        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
!        call xgTransposer_transpose(slice%xgTransposerXext, STATE_COLSROWS)
!        ABI_NVTX_END_RANGE()
!
!        slice%use_colsrows = .true.
!        slice%use_linalg = .false.
!
!        ! Unitary test
!        if ( cols(slice%me_Xext) /= slice%ncolsColsRows(xmpi_comm_rank(slice%spacecom)+1) ) then
!            ABI_ERROR('wrong colsrows representation')
!        end if
!        write(std_out,'(a,i6,i6,i6)') '# proc has # cols of Xext ', xmpi_comm_rank(slice%spacecom),cols(slice%me_Xext)
!        work%has_transposer = .true.
!    else
!        call xgBlock_setBlock(work%XextLinalg, slice%me_Xext, rows(work%XextLinalg), cols(work%XextLinalg))
!    end if
!
!    ! Get parameters of active task
!    neigenpairs = task%me_neigenpairs_slice
!    bandpp = task%me_bandpp_slice
!    ndeg_filter = task%me_ndeg_slice
!    comm = task%me_comm_slice
!    comm_rows = task%me_comm_rows
!    comm_cols = task%me_comm_cols
!    if (slice%me_id_slice==1) then   
!        is_lowpass = .true. ! [lambda_minus,lambda_plus) to be diminished
!        lambda_minus = slice%poly_upp_bounds(task%me_id_slice)
!        lambda_plus = slice%maxeig_global
!    else
!        is_lowpass = .false. ! [lambda_minus,lambda_plus) to be amplified
!        lambda_minus = slice%poly_low_bounds(task%me_id_slice)
!        lambda_plus = slice%poly_upp_bounds(task%me_id_slice)
!    end if
!    ! deactivate chebfi oracle
!    oracle = 0
!    nbdbuf = 0
!    oracle_factor = 1.d0
!    oracle_min_occ = 0.d0
!
!    num_proc = xmpi_comm_size(comm)
!    ABI_MALLOC_IFNOT(nrowsLinalg,(num_proc))
!    nrowsLinalg_ptr => nrowsLinalg
!    nrowsLinalg = task%me_nrowsLinalg_slice
!
!    write(std_out,*) "Allocating slice space..", slice%total_spacedim, neigenpairs
!    write(std_out,*) "bands per process=", bandpp
!    flush(std_out)
!
!    ! Initialize chebfi object in MPI Colsrows distribution
!    call chebfi_init(chebfi,neigenpairs,slice%total_spacedim,slice%tolerance,slice%ecut,slice%paral_kgb,bandpp,&
!        ndeg_filter,nbdbuf,slice%space,1,comm,task%me_g0,task%me_g0_fft,slice%paw,comm_rows,comm_cols,&
!        oracle,oracle_factor,oracle_min_occ,slice%gpu_option,gpu_kokkos_nthrd=slice%gpu_kokkos_nthrd,&
!        gpu_thread_limit=slice%gpu_thread_limit,from_linalg=.false.)

end subroutine allocate_active_task
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/free_active_task
!! NAME
!! free_active_task
!! 
!! SOURCE

subroutine free_active_task(task)

    implicit none
    type(activeTask_t), intent(inout) :: task
    
    ! *********************************************************************

    ABI_SFREE(task%me_ncolsColsRows_slice)   
    ABI_SFREE(task%me_nrowsLinalg_slice)

end subroutine free_active_task
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/execute_active_task
!! NAME
!! execute_active_task
!! 
!! SOURCE

subroutine execute_active_task(task)

    type(activeTask_t), intent(inout) :: task

    type(xgBlock_t) :: X0_active
    type(xgBlock_t) :: eigen_active
    type(xgBlock_t) :: residu_active
    
    ! *********************************************************************

!    ! Recover actively used array
!    X0_active = task%me_Xext
!    
!    ! fixme move this inside chebfi_runSI?
!    task%chebfi%xXColsRows = X0_active
!        
!    !call chebfi_runSlice(chebfi, X0_active, getAX_BX, getBm1X, eigen_active, residu_active, nspinor,&
!    !    slice%mineig_global, slice%maxeig_global, lambda_minus, lambda_plus, is_lowpass, slice%neigenpairs,&
!    !    nrowsLinalg_ptr)
!
!    ! todo give k=m+p where p is oversample
!    ! residual will be converged for m values. Give m as input
!    k_conv = slice%neigenpairs - 20 ! hardcoded assuming offset 20 fixme 
!
!    !call chebfi_runSubspaceIteration(chebfi, X0_active, getAX_BX, getBm1X, eigen_active, residu_active, &
!    !    nspinor, slice%mineig_global, slice%maxeig_global, lambda_minus, lambda_plus, is_lowpass, &
!    !    k_conv, nrowsLinalg_ptr)
!
!    call chebfi_runSubspaceIterationDummy(task%chebfi, X0_active, getAX_BX, getBm1X, eigen_active, residu_active, &
!        nspinor, slice%mineig_global, slice%maxeig_global, lambda_minus, lambda_plus, is_lowpass, &
!        k_conv, nrowsLinalg_ptr)
!
!    ! why?
!    if (slice%gpu_option == ABI_GPU_OPENMP) then
!        call xgBlock_copy_to_gpu(eigen_active)
!        call xgBlock_copy_to_gpu(residu_active)
!    end if


end subroutine execute_active_task
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/mask_active_task
!! NAME
!! mask_active_task
!! 
!! FUNCTION
!! Mask converged solutions in active task (asynchronous).
!! Create mask for converged solutions in active task
!! that will be used to combine all active tasks to extended memory
!!
!! SOURCE

subroutine mask_active_task(task, tol)

    implicit none

    ! Arguments ------------------------------------
    type(activeTask_t), intent(inout) :: task
    real(dp), intent(inout) :: tol

    ! Local variables-------------------------------
    integer :: n_active, iband
    logical :: selected
    real(dp) :: theta, res
    real(dp), pointer :: thetas_conv(:,:) => null()
    real(dp), pointer :: residu_conv(:,:) => null()

    ! *********************************************************************

    n_active = rows(task%me_eigen)
    call xgBlock_reverseMap(task%me_eigen, thetas_conv, rows=n_active, cols=1)
    call xgBlock_reverseMap(task%me_resid, residu_conv, rows=n_active, cols=1)
   
    ! Hard acceptance criterion so that slices do not overlap
    do iband=1, n_active
        res = residu_conv(iband, 1)
        theta = thetas_conv(iband, 1)
        selected = .false.
        if (task%is_lowpass) then
            selected = (res < tol .and. theta < task%me_uppb)
        else if (task%is_last) then
            selected = (res < tol .and. theta > task%me_lowb)
        else
            selected = (res < tol .and. theta < task%me_uppb .and. theta > task%me_lowb)
        end if
        if (selected) then
            task%me_mask_Xext(iband) = 1
        end if
    end do

end subroutine mask_active_task
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/compress_extended_memory
!! NAME
!! compress_extended_memory 
!!
!! FUNCTION
!! Copy data from asyncMemory to spectrum I/O memory (shared)
!! Copy masked async memory to spectrum memory
!!
!! SOURCE

subroutine compress_extended_memory(task, work, X0, eigen, resid)

    implicit none
    
    ! Arguments ------------------------------------
    type(activeTask_t), intent(inout) :: task
    type(asyncMemory_t), intent(inout) :: work
    type(xgBlock_t), intent(inout) :: X0
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: resid

    ! Local variables-------------------------------
    integer :: nrows, ncols, fcol
    
    ! *********************************************************************

!    if (work%paral_kgb==1) then
!        
!        ! Sanity check
!        if ((.not. work%use_linalg) .or. work%use_colsrows) ) then
!            ABI_ERROR("not in linalg")
!        end if
!
!        ! Recover pointer task%me_Xext into memory work%XextLinalg
!        call xmpi_barrier(slice%spacecom)
!        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
!        call xgTransposer_transpose(work%xgTransposerXext, STATE_LINALG)
!        ABI_NVTX_END_RANGE()
!    else
!
!        nrows = 
!        ncols = 
!        fcol = 
!        xgBlock_setBlock(work%XextLinalg, task%me_Xext, rows(work%XextLinalg), cols(
!
!    end if
!
!    ! Detect missing or extra eigenvalues
!    if (tot_ncols_kept < slice%neigenpairs) then
!        ABI_WARNING("Not enough converged eigenvalues in slice")
!    else if (tot_ncols_kept > slice%neigenpairs) then
!        ABI_WARNING("Too many converged eigenvalues kept. Decrease tolfilter or nstep_mixed.")
!    end if    
!
!    ! Copy from async memory to regular memory
!    do islice=1,slice%nslice
!        fcol = slice%fcol_in_X(islice)
!        fcol_ext = slice%fcol_in_Xext(islice)
!        neigenpairs_slice = slice%neigenpairs_per_slice(islice)
!        write(std_out,*) 'block copy from fcol, ncols=', fcol_ext, neigenpairs_slice
!        write(std_out,*) 'block copy to fcol, ncols=', fcol, neigenpairs_slice
!        ! Blocks to copy from
!        call xgBlock_setBlock(slice%XextLinalg, X_kept, rows=slice%spacedim, cols=neigenpairs_slice, fcol=fcol_ext)
!        call xgBlock_setBlock(eigen_ext%self, eigen_kept, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
!        call xgBlock_setBlock(resid_ext%self, resid_kept, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
!        ! Blocks to copy to
!        call xgBlock_setBlock(X0, X0_out, rows=slice%spacedim, cols=neigenpairs_slice, fcol=fcol)
!        call xgBlock_setBlock(eigen, eigen_out, rows=1, cols=neigenpairs_slice, fcol=fcol)
!        call xgBlock_setBlock(resid, resid_out, rows=1, cols=neigenpairs_slice, fcol=fcol)
!        ! copy
!        call xgBlock_copy(X_kept, X0_out)
!        call xgBlock_copy(eigen_kept, eigen_out)
!        call xgBlock_copy(resid_kept, resid_out)
!    end do
!
!    ! Recover dimensions
!    call xgBlock_reshape(eigen, slice%neigenpairs, 1) 
!    call xgBlock_reshape(resid, slice%neigenpairs, 1)
!
!    ! Free memory
!    call xg_free(eigen_ext)
!    call xg_free(resid_ext)
!
end subroutine compress_extended_memory
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/assign_tasks_to_processes
!! NAME
!! assign_tasks_to_processes
!! 
!! FUNCTION
!! Perform the inverse of the allocation operation, assigning 
!! processes to slices based on the allocation array.
!! 
!! SOURCE

subroutine assign_tasks_to_processes(allocations, processes)

    implicit none
    integer, intent(in) :: allocations(:)
    integer, intent(out) :: processes(:)
    integer :: i, j
    
    ! *********************************************************************

    j = 1
    do i = 1, size(allocations)
        processes(j:j + allocations(i) - 1) = i - 1  
        j = j + allocations(i)
    end do

end subroutine assign_tasks_to_processes
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/distribute_vectors
!! NAME
!! distribute_vectors
!!
!! FUNCTION
!! Distribute m vectors as uniformly as possible across n processes.
!! Assumptions:
!! -The remainder should be distributed evenly to the first n-1 processes.
!! -The last process should always get fewer vectors.
!! -The sum of the allocations should be exactly m.
!! 
!! SOURCE

subroutine distribute_vectors(m, n, allocation)

    implicit none
    integer, intent(in) :: m, n
    integer, intent(out) :: allocation(n)
    integer :: i, base, remainder
    
    ! *********************************************************************

    base = m / n
    remainder = m - base * n
    allocation = base
    if (remainder /= 0) then
        do i = 1, n-1
            if (remainder > 0) then
                allocation(i) = allocation(i) + 1
                remainder = remainder - 1
            end if
        end do
    end if

end subroutine distribute_vectors
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/fair_allocation
!! NAME
!! fair_allocation
!! 
!! FUNCTION
!! Solve integer optimization problem under constraint: 
!! 
!!     min_{x_1,..,x_s} max_{1,..,s} f_i(x_i)
!!     subject to:   x_1 + .. + x_s = p
!!                   x_i integers
!! 
!! with objective cost function f_i(x)=m_i*w_i/x.
!! The solution x_i is the amount of resource allocated to the i-th task.
!! The algorithm uses binary search for integer rounding. 
!! Note that this is better than greedy but not optimal. 
!! Exhaustive search is too expensive (=(p+1)^n combinations).
!! Feature: ensures the total allocation is exactly equal to p while 
!! minimizing the allocation imbalance.
!! 
!! INPUTS
!! arrays m and w (length n), and integer p
!! m can be the group size, w can be another measure or need (weight)
!! p in the number of total resources
!! 
!! OUTPUT
!! integer array x of size s such that sum(x) = p and max(m_i*n_i/x_i) is minimized
!!
!! SOURCE

subroutine fair_allocation(n, m, w, p, x)

    implicit none

    ! Arguments
    integer, intent(in) :: n            ! Number of groups (slices)
    integer, intent(in) :: m(n)         ! Array: size of each group
    integer, intent(in) :: p            ! Total resources to allocate
    integer, intent(in) :: w(n)         ! Weight of each group
    integer, intent(out) :: x(n)        ! Array: allocated resources per group

    ! Local variables
    integer :: i, total_allocated
    real(dp) :: total_work, lower, upper, mid, multiplier
    integer :: allocation(n)

    ! *********************************************************************

    total_work = dot_product(m, w)

    ! Binary search for the optimal multiplier: divide search interval in half
    lower = 0.d0
    upper = real(p)
    do while (upper - lower > 1.d0)
        mid = (lower + upper) / 2.d0
        do i = 1, n
            allocation(i) = int((real(w(i)) * real(m(i)) / total_work) * mid + 0.d5)
        end do
        total_allocated = sum(allocation(1:n))

        ! Adjust binary search bounds
        if (total_allocated > p) then
            upper = mid
        else
            lower = mid
        end if
    end do

    ! Fair rounding: tasks with heavier workload get extra resource units

    ! Final allocation after binary search converges
    multiplier = (lower + upper) / 2.d0
    do i = 1, n
        allocation(i) = int((real(w(i)) * real(m(i)) / total_work) * multiplier + 0.d5)
    end do

    ! Adjust total allocation to exactly match p
    total_allocated = sum(allocation(1:n))

    if (total_allocated < p) then
        do while (total_allocated < p)
            ! Add one resource to the group closest to its ideal allocation
            call adjust_allocation(n, m, w, allocation, total_work, p, total_allocated)
            total_allocated = sum(allocation(1:n))
        end do
    else if (total_allocated > p) then ! FIXME error infinite loop
        do while (total_allocated > p)
            ! Remove one resource from the over-allocated group
            call reduce_allocation(n, m, w, allocation, total_work, p, total_allocated)
            total_allocated = sum(allocation(1:n))
        end do
    end if

    ! Assign the final allocation to the output variable
    x = allocation

    !write(std_out,'(a)') 'Memory allocation info:'
    !do i=1,n
    !    if (x(i) .ne. 0) then
    !        write(std_out,'(a,i4,i6,i5)') '#task #workload #allocated resources', i, w(i)*m(i)/x(i), x(i)
    !    else
    !        write(std_out,'(a,i4,i4,i4,i4)') 'Allocation error: x(i)= w(i)= m(i)= for i=', x(i), w(i), m(i), i
    !    end if
    !end do

end subroutine fair_allocation
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/adjust_allocation
!! NAME
!! adjust_allocation
!! 
!! FUNCTION
!! Adjust allocation by adding resources to the group closest to its ideal allocation

subroutine adjust_allocation(n, m, w, allocation, total_weight, p, total_allocated)

    implicit none

    integer, intent(in) :: n, m(n)
    integer, intent(in) :: w(n), p
    real(dp), intent(in) :: total_weight
    integer, intent(inout) :: allocation(n)
    integer, intent(inout) :: total_allocated

    integer :: i, closest_group
    real(dp) :: max_diff, diff, ideal

    ! *********************************************************************

    ! Find the group with the largest difference between current and ideal allocation
    max_diff = -1.0_dp
    closest_group = -1

    do i = 1, n
        ideal = (real(w(i), dp) * real(m(i), dp) * real(p, dp)) / total_weight
        diff = abs(real(allocation(i), dp) - ideal)
        if (diff > max_diff) then
            max_diff = diff
            closest_group = i
        end if
    end do

    ! Only allocate if total does not exceed limit
    if (total_allocated < p .and. closest_group > 0) then
        allocation(closest_group) = allocation(closest_group) + 1
        total_allocated = total_allocated + 1
    end if

end subroutine adjust_allocation
!!***

!----------------------------------------------------------------------

!!****f* m_slice_task/reduce_allocation
!! NAME
!! reduce_allocation
!! 
!! FUNCTION
!! Reduce allocation by removing resources from the over-allocated group

subroutine reduce_allocation(n, m, w, allocation, total_weight, p, total_allocated)
  
    implicit none

    integer, intent(in) :: n, m(n)
    integer, intent(in) :: w(n), p
    real(dp), intent(in) :: total_weight
    integer, intent(inout) :: allocation(n)
    integer, intent(inout) :: total_allocated

    integer :: i, target_group
    real(dp) :: max_diff, expected, diff
    integer :: reduce_by

    ! *********************************************************************

    ! Find the group with the largest *positive* over-allocation
    max_diff = -1.0_dp
    target_group = -1

    do i = 1, n
        expected = (real(w(i), dp) * real(m(i), dp) * real(p, dp)) / total_weight
        diff = real(allocation(i), dp) - expected
        if (diff > max_diff .and. diff > 0.0_dp .and. allocation(i) > 0) then
            max_diff = diff
            target_group = i
        end if
    end do

    ! If we found an over-allocated group, reduce its allocation
    if (target_group > 0 .and. max_diff > 0.0_dp) then
        reduce_by = min(1, allocation(target_group))  ! Only reduce if it's > 0
        reduce_by = min(1, allocation(target_group))
        allocation(target_group) = allocation(target_group) - reduce_by
        total_allocated = total_allocated - reduce_by
    end if

end subroutine reduce_allocation
!!***

end module m_slice_task
!!***

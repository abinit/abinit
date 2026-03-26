!!****f* ABINIT/m_slice
!! NAME
!! m_slice
!!
!! FUNCTION
!! This module contains the types and routines used to apply the Spectrum Slicing 
!! method. It mainly defines 'slice' datatypes and associated methods. 
!!
!! Main features:
!! - uses 'xgTools' implementation for matrix data structure.
!! - adopts 'chebfi' functionalities for most matrix calculations.
!! - uses polynomial filtering, Chebyshev lowpass, next Chebyshev-Jackson.
!! - uses scheduler and resource allocator for slice tasks.
!! - uses logic for asynchronous slice memory avoiding race condition in read/write.
!! - applies Rayleigh-Ritz for individual slices in parallel or sequentially.
!!
!! NOTES
!! Dependence with other modules in this directory:
!!   m_slice uses:
!!      |- m_polynomial_filter
!!      |- m_task_scheduler
!!      |- m_trace_estimation
!!      |- m_chebfi2
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

module m_slice

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
    use m_chebfi2

    use m_trace_estimation
    use m_task_scheduler
    use m_polynomial_filter

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

    ! Spectral interval decomposition strategy (spectral_cut input variable)
    ! todo delete this
    !---------------------------------------------------
    integer, parameter :: DIVIDE_INTERVAL_WIDTH    = 0
    integer, parameter :: DIVIDE_NUMBER_OF_VECTORS = 1
    integer, parameter :: DIVIDE_SPECTRAL_GAPS     = 2

    ! Load balance criteria for fair resource allocation (paral_slice input variable)
    ! todo should be either 0 or 1
    !---------------------------------------------------
    integer, parameter :: FAIR_BANDPP      = 1              ! bandpp (unweigted)
    integer, parameter :: FAIR_BANDPP_WDEG = 2              ! bandpp weighted by degree
    integer, parameter :: EVEN_SLICES      = 3              ! all slices have the same num of procs
    integer, parameter :: SEQUENTIAL_SLICE = 4              ! sequential slice treatment

    integer, parameter :: ENABLE_PARAL  = 0
    integer, parameter :: DISABLE_PARAL = 1

    ! Timers
    !---------------------------------------------------
    integer, parameter :: tim_swap        = 1761
    integer, parameter :: tim_RR_q        = 1759
    integer, parameter :: tim_barrier     = 1764
    integer, parameter :: tim_copy        = 1765
    integer, parameter :: tim_Bortho_X    = 1641
    integer, parameter :: tim_getAX_BX    = 1754
    integer, parameter :: tim_invovl      = 1755
    integer, parameter :: tim_slice_sched = 2161
    integer, parameter :: tim_slice_1     = 2162
    integer, parameter :: tim_slice_2     = 2163
    integer, parameter :: tim_slice_3     = 2164
    integer, parameter :: tim_slice_X     = 2165
    integer, parameter :: tim_lanczos     = 2167
    integer, parameter :: tim_trace       = 2168
      
    ! Public 'slice' datatype
    !-------------------------------------------------
    type, public :: slice_t

        ! MPI-related common to all-slice phase
        integer :: nproc                                    ! number of available processes
        integer :: comm_rows                                ! xmpi_comm_self ...
        integer :: comm_cols                                ! same as spacecom
        integer :: spacecom                                 ! same as comm_cols
        integer :: bandpp                                   ! nb of bands per process in colsrows representation 
        integer :: me_g0                                    ! 1 if this processors treats G=0, 0 otherwise
        integer :: me_g0_fft                                ! 1 if this processors treats G=0 in FFT, 0 otherwise

        ! Eigenpair parameters
        integer :: neigenpairs                              ! total number of bands (=number of eigenpairs)
        integer :: total_spacedim                           ! total number of plane-waves
        integer :: spacedim                                 ! nb of plane-waves per process in linalg representation
        integer :: nslice                                   ! number of spectral slices
        integer :: space                                    ! real or complex eigenvectors
        integer :: nbdbuf                                   ! fixme for the moment not clear if useful
        
        ! GPU-related
        integer :: gpu_kokkos_nthrd                 
        integer :: gpu_thread_limit 
        
        ! Various options
        integer :: gpu_option                               ! enable GPU
        integer :: paral_kgb                                ! enable parallel (k-points, G basis, bands)
        integer :: paral_slice                              ! how to allocate resources for parallel slices
        integer :: spectral_cut                             ! how to decompose spectrum into slices
        
        ! Various model parameters
        logical :: paw                                      ! use PAW or not 
        integer :: ndeg_filter                              ! lowpass degree of polynomial filter
        real(dp) :: tolerance                               ! tolerance on the residu to stop the minimization
        real(dp) :: ramp                                    ! bandpass filter tolerance
        real(dp) :: ecut                                    ! Ecut Fermi level
        real(dp) :: mineig_global                           ! guaranteed lower spectral bound
        real(dp) :: maxeig_global                           ! guaranteed upper spectral bound

        ! Arrays related to polynomial filtering
        integer, allocatable :: neigenpairs_per_slice(:)    ! number of total eigenpairs per slice
        integer, allocatable :: poly_degrees(:)            ! polynomial filter degrees
        real(dp), allocatable :: low_bounds(:)             ! lower bounds in spectral partition (disjoint)
        real(dp), allocatable :: upp_bounds(:)             ! upper bounds in spectral partition (disjoint) 

        ! Main object of per-slice phase 
        type(activeSlice_t) :: activeSlice

        ! Main objects of all-slice phase
        type(taskScheduler_t) :: taskScheduler
        type(extendedMemory_t) :: extendedMemory

    end type slice_t

    ! Public methods associated to 'slice' datatype
    !-------------------------------------------------
    public :: slice_init                                    ! initialize slice datatype object
    public :: slice_free                                    ! free slice datatype object
    public :: slice_run                                     ! run Spectrum Slicing for active slice task

    CONTAINS  
!=====================================================================
!!***

!!****f* m_slice/slice_init
!! NAME
!! slice_init
!!
!! FUNCTION
!! Initialize a 'slice' datastructure.
!!
!! INPUTS
!!  nslice= number of spectral slices
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  spacedim= space dimension for one vector
!!  tolerance= tolerance criterion on the residu to stop the minimization
!!  ecut= plane-wave cut-off energy
!!  paral_kgb= flag controlling (k,g,bands) parallelization
!!  bandpp= number of 'bands' handled by a processor
!!  ndeg_filter= polynomial degree of the polynomial filter (.i.e. number of H applications)
!!  mineig_global= lower spectral bound, to be scaled to -1
!!  maxeig_global= upper spectral bound, to be scaled to 1
!!  lambda_minus= lower bound of interval to amplify
!!  lambda_plus= upper bound of interval to amplify
!!  is_lowpass= flag. True: use lowpass Chebyshev, false: use bandpass Heaviside expanded on Chebyshev
!!  space= defines in which space we are (columns, rows, etc.)
!!  eigenProblem= type of eigenpb: 1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
!!  spacecom= MPI communicator
!!  me_g0= 1 if this processors treats G=0, 0 otherwise
!!  me_g0_ftt= 1 if this processors treats G=0 in FFT, 0 otherwise
!!  paw= flag. TRUE if current calculation ses the PAW approach
!!  comm_rows= "rows" communicator
!!  comm_cols= "cols" communicator
!!  gpu_option= flag. Enable GPU if true
!!  gpu_kokkos_nthrd= number of OpenMP offloaded threads used
!!  gpu_thread_limit= maximum number of OpenMP offloaded threads
!!
!! SOURCE

subroutine slice_init(slice,nslice,neigenpairs,spacedim,tolerance,paral_kgb,&
        paral_slice,ndeg_filter,nbdbuf,ramp,ecut,bandpp,space,spacecom,me_g0,me_g0_fft,&
        paw,comm_rows,comm_cols,spectral_cut,gpu_option,gpu_kokkos_nthrd,gpu_thread_limit)

    implicit none

    ! Arguments ------------------------------------
    integer      , intent(in   ) :: bandpp
    integer      , intent(in   ) :: nslice
    integer      , intent(in   ) :: me_g0
    integer      , intent(in   ) :: me_g0_fft
    integer      , intent(in   ) :: neigenpairs
    integer      , intent(in   ) :: comm_rows
    integer      , intent(in   ) :: comm_cols
    integer      , intent(in   ) :: paral_kgb
    integer      , intent(in   ) :: paral_slice
    integer      , intent(in   ) :: space
    integer      , intent(in   ) :: spacecom
    integer      , intent(in   ) :: spacedim
    integer      , intent(in   ) :: ndeg_filter
    integer      , intent(in   ) :: nbdbuf
    integer      , intent(in   ) :: spectral_cut
    integer      , intent(in   ) :: gpu_option
    logical      , intent(in   ) :: paw
    real(dp)     , intent(in   ) :: ramp
    real(dp)     , intent(in   ) :: tolerance
    real(dp)     , intent(in   ) :: ecut
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ), optional :: gpu_kokkos_nthrd
    integer      , intent(in   ), optional :: gpu_thread_limit
    
    ! Local variables --------------------------------
    integer :: total_spacedim, ierr
    logical :: on_host, on_device

    ! *********************************************************************

    if (paral_slice==0) then
        ABI_ERROR("Sequential slices not implemented.")
    end if

    slice%bandpp        = bandpp
    slice%nslice        = nslice
    slice%me_g0         = me_g0
    slice%me_g0_fft     = me_g0_fft
    slice%neigenpairs   = neigenpairs
    slice%comm_rows     = comm_rows
    slice%comm_cols     = comm_cols
    slice%paral_kgb     = paral_kgb
    slice%paral_slice   = paral_slice
    slice%space         = space
    slice%spacecom      = spacecom
    slice%spacedim      = spacedim
    slice%ndeg_filter   = ndeg_filter
    slice%nbdbuf        = nbdbuf
    slice%spectral_cut  = spectral_cut
    slice%gpu_option    = gpu_option
    slice%paw           = paw
    slice%ramp          = ramp
    slice%tolerance     = tolerance
    slice%ecut          = ecut

    slice%gpu_kokkos_nthrd = 1
    if (present(gpu_kokkos_nthrd)) slice%gpu_kokkos_nthrd = gpu_kokkos_nthrd
    slice%gpu_thread_limit = 0
    if (present(gpu_thread_limit)) slice%gpu_thread_limit = gpu_thread_limit

    ! Total number of rows (used in colsrows representation)
    if (paral_kgb==0) then
        slice%total_spacedim = spacedim
    else
        total_spacedim = spacedim
        call xmpi_sum(total_spacedim, spacecom, ierr)
        slice%total_spacedim = total_spacedim
    end if

    ! Total number of processes
    slice%nproc = xmpi_comm_size(spacecom)

    ! Arrays
    call slice_allocateAll(slice)

end subroutine slice_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allocateAll
!! name
!! slice_allocateAll

subroutine slice_allocateAll(slice)

    implicit none
 
    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice

    ! *********************************************************************

    call slice_free(slice)

    ABI_MALLOC_IFNOT(slice%neigenpairs_per_slice, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%poly_degrees, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%low_bounds, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%upp_bounds, (slice%nslice))

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_free
!! name
!! slice_free

subroutine slice_free(slice)

    implicit none
    type(slice_t), intent(inout) :: slice

    ! *********************************************************************

    ABI_SFREE(slice%neigenpairs_per_slice)
    ABI_SFREE(slice%poly_degrees)
    ABI_SFREE(slice%low_bounds)
    ABI_SFREE(slice%upp_bounds)

end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_run
!! NAME
!! slice_run
!!
!! FUNCTION
!! Run Spectrum slicing on a given set of active vectors. 
!! Computation uses marked resources for current task **only**.
!! 
!! INPUTS
!! X=            size (spacedim, neigenpairs)
!! eigen,residu= size (neigenpairs, 1) (column vectors)
!! getAX_BX= pointer to the function giving A|X> and B|X>
!!           A is typically the Hamiltonian H, and B the overlap operator S
!! getBm1X= pointer to the function giving B^-1|X>
!!          B is typically the overlap operator S
!! nspinor= number of spinorial components of the wavefunctions
!! 
!! SIDE EFFECTS
!! slice <type(slice_t)>= memory workspace used for Spectrum slicing
!! 
!! SOURCE

subroutine slice_run(slice, X, getAX_BX, getBm1X, eigen, residu, nspinor)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X      ! size (spacedim, neigenpairs)
    type(xgBlock_t), intent(inout) :: eigen  ! size (neigenpairs,1)
    type(xgBlock_t), intent(inout) :: residu ! size (neigenpairs,1)
    integer, intent(in) :: nspinor
    interface
        subroutine getAX_BX(X,AX,BX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
            type(xgBlock_t), intent(inout) :: BX
        end subroutine getAX_BX
    end interface
    interface
        subroutine getBm1X(X,Bm1X)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: Bm1X
        end subroutine getBm1X
    end interface

    ! Local scalars
    type(chebfi_t) :: chebfi
    integer :: tim_slice_me
    integer :: i, iband, nrows
    integer :: me_nbdbuf
    integer :: neigenpairs, bandpp, ndeg_filter
    integer :: comm, comm_rows, comm_cols
    integer :: num_restart
    integer :: num_kept
    integer :: k_conv
    integer :: ierr
    logical :: has_converged
    real(dp) :: lambda_minus, lambda_plus
    real(dp) :: theta
    real(dp) :: safe, tol ! for slice selection window
    real(dp) :: a_part, b_part, max_resid_kept
    logical :: is_lowpass, on_host, on_device
    ! todo use slice%..
    type(activeSlice_t) :: task
    type(taskScheduler_t) :: scheduler
    type(extendedMemory_t) :: extendedMemory
    ! Arrays
    real(dp) :: tsec(2)
    integer, allocatable :: mapper(:,:)
    integer, allocatable, target :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null() 
    integer, pointer :: ncolsColsRows_ptr(:) => null()
    
    ! *********************************************************************
    
    if (task%me_id_slice==1) then
        tim_slice_me = tim_slice_1
    else if (task%me_id_slice==2) then
        tim_slice_me = tim_slice_2
    else if (task%me_id_slice==3) then
        tim_slice_me = tim_slice_3
    else
        tim_slice_me = tim_slice_X
    end if
    
    nrows = slice%spacedim
    if (slice%paral_kgb==1) then
        nrows = slice%total_spacedim
    end if

    ! ================================== Prepare spectral slices ===========================================

!    call timab(tim_slice_sched,1,tsec)
!    ABI_NVTX_START_RANGE(NVTX_SLICE_SCHEDULE)
!
!    ! Split spectrum and query X0
!    ABI_NVTX_START_RANGE(NVTX_SLICE_RRQ)
!    call slice_prepareSpectrum(slice, X0, eigen, residu, getAX_BX, getBm1X, nspinor)
!    ABI_NVTX_END_RANGE()
!
!    ABI_NVTX_END_RANGE()
!    call timab(tim_slice_sched,2,tsec)
!
!    ! Attribute column vectors of X0 to slices using query results
!    ABI_MALLOC(mapper, (slice%neigenpairs, slice%nslice))
!    safe = 2
!    tol = 1.d0
!    call slice_applySelectionWindow(slice, mapper, eigen, residu, tol, safe)
!
!    write(std_out,*) 'created the following map slices to cols='
!    write(std_out,*) mapper(1,:)
!    flush(std_out)

    !! IML debug start
    !! debug up to here is independent of parallelism

    ! what to do next
    ! sketch remaining offset using random vectors
    !k_sketch = k - m_conv

    ! IML here start the complicated part that depends on parallelism
 
    ! ============================== Initialize guess for subspace iteration ===============================

    ! Resource management system 
!    if (slice%paral_kgb==1) then
!        call init_schedule(scheduler)
!        ! Divide resources into slice tasks
!        select case(slice%paral_slice)
!        case(DISABLE_PARAL)
!            nband_ext = max(k+p)
!            call schedule_next_task(scheduler)
!        case(ENABLE_PARAL)
!            nband_ext = xmpi_sum(k+p)
!            call schedule_parallel_tasks(scheduler)
!            call xmpi_barrier(slice%spacecom)
!            do iproc = 1, slice%nproc
!                write(std_out,'(a,i5,a,i5)') "Process ", iproc-1, " allocated to task ", slice%lookup_proc(iproc)
!            end do
!        end select
!        ! Mark my slice task and resources as actively in use
!        call mark_active_task(scheduler, task) ! called once for parallel slices and in loop for sequential
!
!        ! Allocate and fill extended memory buffer
!        nband_ext = sum(slice%neigenpairs_per_slice)
!        call allocate_extended_memory(extendedMemory, nband_ext)
!
!        ! k is the number of columbs from X0
!        ! p is the offset
!        call init_extended_memory(extendedMemory, task, X0, k, p)
!
!    else
!        ! configure task without any communicators
!    end if
!
!    call allocate_active_task(slice, scheduler, task)
!
!  
!    ! will use mapper to copy columns of X0 to extended memory
!    !! assumes linalg representation of both ext and spectrum mem
!    call init_extended_memory(extendedMemory, X0, mapper)
!    
!    ! ============================ Active task execution =======================================
!
!    call timab(tim_slice_me,1,tsec)
!    
!    if (slice%paral_kgb==0 .or. slice%paral_slice==DISABLE_PARAL) then
!        ! execute active tasks sequentially
!
!        do islice=1, nslice
!            call schedule_next_task(scheduler, neigenpairs)
!            call init_active_memory(extendedMemory, task, X0, p)
!            call execute_active_task(task)
!            call mask_active_task(task, tol) ! mask extendedMem
!        end do
!       
!    else
!        ! execute active tasks in parallel
!
!        call init_active_memory(extendedMemory, task, X0, p)
!        call execute_active_task(task)
!        call mask_active_task(task, tol) ! mask extendedMem
!
!    end if
!   
!    ! Copy to spectrum memory only when active tasks have finished
!    call compress_extended_memory(work) ! extendendMem -> spectrumMem
!
!
!    call free_schedule(scheduler)
!
!
!
!    ! Free temporary memory
!    call chebfi_free(chebfi)
!    ABI_SFREE(nrowsLinalg)
!
!    ! Timer is BEFORE the barrier !!
!    call timab(tim_slice_me,2,tsec)
!
!    call free_extended_memory(extendedMemory)
!    ABI_FREE(mapper)

end subroutine slice_run
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_prepareSpectrum
!! NAME
!! slice_prepareSpectrum
!!
!! FUNCTION
!! Compute spectral intervals, cut them in half etc.
!! 
!! INPUTS
!! X        =eigenvector guess
!! getAX_BX =Hamiltonian application
!! 
!! SOURCE

subroutine slice_prepareSpectrum(slice, X, eigen, resid, getAX_BX, getBm1X, nspinor)

    implicit none

    ! Arguments
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: resid
    integer, intent(in) :: nspinor
    interface
        subroutine getAX_BX(X,AX,BX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
            type(xgBlock_t), intent(inout) :: BX
        end subroutine getAX_BX
    end interface
    interface
        subroutine getBm1X(X,Bm1X)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: Bm1X
        end subroutine getBm1X
    end interface

    ! Local variables
    ! Scalars
    integer :: spacedim, tot_spacedim
    integer :: ndeg_filter_max, neigenpairs
    integer :: ierr
    integer :: ideg
    integer :: nstep_bisect
    integer :: kmax
    integer :: m_probe
    integer :: my_shift, my_rank
    integer :: space, spacecom, gpu_option
    integer :: k_sketch
    real(dp) :: lambda_min, res_norm
    real(dp) :: center, radius
    real(dp) :: lowb, mineig, maxeig
    real(dp) :: mineig_global, maxeig_global
    real(dp) :: lanczos_lowb, lanczos_lowb_global
    ! Derived types
#ifdef HAVE_OPENMP_OFFLOAD
    integer :: me_g0
    type(xg_t) :: W_dummy
    integer :: work_size
#endif
    type(xg_t) :: BX
    type(xgBlock_t) :: xXColsRows
    type(xgBlock_t) :: eigen_me, resid_me
    type(matrixInfo_t) :: matrixInfo
    type(xgTransposer_t) :: xgTransposerX
    ! Arrays
    complex(dp), allocatable :: cheby_moments(:,:)
    real(dp) :: tsec(2)

    ! *********************************************************************

    neigenpairs = slice%neigenpairs
    tot_spacedim = slice%total_spacedim
    spacedim = slice%spacedim
    gpu_option = slice%gpu_option
    spacecom = slice%spacecom
    space = slice%space
#ifdef HAVE_OPENMP_OFFLOAD
    me_g0 = slice%me_g0
    if (slice%paral_kgb==1) then
        me_g0 = slice%me_g0_fft
    end if
#endif

    ! Y = X * Omega where Omega sketch matrix to capture all directions at once
    !k_sketch = neigenpairs
    !call xg_init(X_sketch, space, spacedim, neigenpairs, spacecom, gpu_option=gpu_option)
    !call randomSketching(slice, X, X_sketch%self, k_sketch)
    !call xgBlock_copy(X_sketch%self, X)
    !call xg_free(X_sketch)
    ! fixme this has been moved elsewhere

    ! ============== Transpose ==============
    if (slice%paral_kgb==1) then

        ! Allocate memory for X in colsrows representation
        call xgTransposer_constructor(xgTransposerX, X, xXColsRows, nspinor, STATE_LINALG,&
            TRANS_ALL2ALL, slice%comm_rows, slice%comm_cols, 0, 0, slice%me_g0_fft,&
            gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit)
         
        xgTransposerX%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
        
        !call xmpi_barrier(slice%spacecom)
        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
        call xgTransposer_transpose(xgTransposerX, STATE_COLSROWS)
        ABI_NVTX_END_RANGE()

    else

        ! Use colsrows notion instead of X notion
        call xgBlock_setBlock(X, xXColsRows, spacedim, neigenpairs)

    end if

    ! Prevent invovl error
    ! dummy invovl calculation to allocate buffers of full size
    ! Workaround: getBm1X error related to invovl allocated buffers. 
    ! With this setup we allocate large then use 1 in Lanczos.
    ! The inverse is not possible with current implementation. This order allows to avoid
    ! make_invovl for 1 vector then make_invovl for nband vectors.
#ifdef HAVE_OPENMP_OFFLOAD
    if (slice%paw) then
        work_size = maxval(slice%neigenpairs_per_slice)
        call xg_init(W_dummy, slice%space, tot_spacedim, work_size, xmpi_comm_null, &
            me_g0=me_g0, gpu_option=slice%gpu_option)
        call timab(tim_invovl, 1, tsec)
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_BM1X)
        call getBm1X(W_dummy%self, W_dummy%self)
        ABI_NVTX_END_RANGE()
        call timab(tim_invovl, 2, tsec)
        call xg_free(W_dummy)
    end if
#endif

    call init_matrixInfo(matrixInfo, slice%comm_rows, slice%comm_cols, slice%spacecom, slice%neigenpairs,&
        slice%total_spacedim, slice%spacedim, slice%space, slice%gpu_kokkos_nthrd, slice%gpu_thread_limit,&
        slice%gpu_option, slice%paral_kgb, slice%me_g0, slice%me_g0_fft)

    kmax = 30
    call computeBLanczos(matrixInfo, slice%paw, getAX_BX, getBm1X, kmax, lambda_min, res_norm)

    lanczos_lowb = lambda_min - res_norm
    call xmpi_min(lanczos_lowb,lanczos_lowb_global,slice%spacecom,ierr)
    lowb = lanczos_lowb_global

    write(std_out,*) 'Lanczos lambda_min=', lambda_min
    write(std_out,*) 'Lanczos res_norm  =', res_norm
    write(std_out,*) 'Lanczos guarantee =', lanczos_lowb
    write(std_out,*) 'Lanczos guarantee(global) =', lanczos_lowb_global
    flush(std_out)

    ! Perform sensitivity study of trace estimation for these parameters
    ! Keep m_probe small allows to reduce noise
    ndeg_filter_max = 25
    m_probe = 5

    write(std_out,*) 'Here I compute Stochastic Trace Estimation'
    write(std_out,*) 'm_probe=    ', m_probe
    write(std_out,*) 'ndeg_filter=', ndeg_filter_max
    flush(std_out)

    call computeTraceEstimation(matrixInfo, slice%nslice, slice%ecut, slice%paw, &
        slice%tolerance, getAX_BX, getBm1X, ndeg_filter_max, m_probe, lanczos_lowb_global)
    
    write(std_out,*) 'STE exited'
    flush(std_out)
 
    ! Compute Rayleigh quotients (colsrows distribution)
    !ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RRQ)
    !call timab(tim_RR_q, 1, tsec)

    if (slice%paral_kgb==1) then
        my_shift = my_rank * slice%bandpp
        call xgBlock_reshape(eigen, 1, slice%neigenpairs)
        call xgBlock_reshape(resid, 1, slice%neigenpairs)
        call xgBlock_setBlock(eigen, eigen_me, 1, slice%bandpp, fcol=my_shift)
        call xgBlock_setBlock(resid, resid_me, 1, slice%bandpp, fcol=my_shift)
        call xgBlock_reshape(eigen, slice%neigenpairs, 1)
        call xgBlock_reshape(resid, slice%neigenpairs, 1)
    else
        call xgBlock_setBlock(eigen, eigen_me, rows(eigen), 1)
        call xgBlock_setBlock(resid, resid_me, rows(resid), 1)
    end if

    call slice_queryCandidates(slice, xXColsRows, eigen_me, resid_me, getAX_BX)
    !call timab(tim_RR_q, 2, tsec)
    !ABI_NVTX_END_RANGE()

    write(std_out,*) 'write residual0='
    call xgBlock_print(resid, std_out)
    flush(std_out)

    write(std_out,*) 'write eigen0='
    call xgBlock_print(eigen, std_out)
    flush(std_out)

    ! ============== Transpose ==============
    if (slice%paral_kgb == 1) then
        call xmpi_barrier(slice%spacecom)
        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
        call xgTransposer_transpose(xgTransposerX, STATE_LINALG)
        ABI_NVTX_END_RANGE()
        
        ! reset buffers to right address
        if (xmpi_comm_size(slice%spacecom) == 1) then
            call xgBlock_setBlock(xXColsRows, X, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(xXColsRows, X, spacedim, neigenpairs)
    end if

    if (slice%paral_kgb == 1) then
        call xgTransposer_free(xgTransposerX)
    end if

end subroutine slice_prepareSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_queryCandidates
!! NAME
!! slice_queryCandidates
!! 
!! FUNCTION
!! Compute Rayleigh-Ritz quotients and residuals in colsrows MPI representation.
!! Return bandpp indicators per MPI process.
!! These two indicators are useful to rank candidates per slices.
!! 
!! SOURCE

subroutine slice_queryCandidates(slice, xXColsRows, eigen, resid, getAX_BX)

    implicit none

    !Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: xXColsRows
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: resid
    interface
        subroutine getAX_BX(X,AX,BX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
            type(xgBlock_t), intent(inout) :: BX
        end subroutine getAX_BX
    end interface

    !Local variables-------------------------------
    type(xg_t) :: XAB
    type(xg_t) :: Results1
    type(xg_t) :: Results2
    type(xg_t) :: norml
    type(xgBlock_t) :: xAXColsRows
    type(xgBlock_t) :: xBXColsRows
    type(xgBlock_t) :: xRXColsRows
    integer :: space_res
    integer :: me_g0, nrows, ncols, gpu_option
    real(dp) :: tsec(2)

    ! *********************************************************************
    
    if (slice%space==SPACE_C) then
        space_res = SPACE_C
    else if (slice%space==SPACE_CR) then
        space_res = SPACE_R
    else
        ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
    end if

    ncols = slice%neigenpairs
    nrows = slice%spacedim
    me_g0 = slice%me_g0
    gpu_option = slice%gpu_option
    if (slice%paral_kgb==1) then
        ncols = slice%bandpp
        nrows = slice%total_spacedim
        me_g0 = slice%me_g0_fft
    end if
    call xg_init(Results1, space_res, ncols, 1, gpu_option=slice%gpu_option)
    call xg_init(Results2, space_res, ncols, 1, gpu_option=slice%gpu_option)
    call xg_init(norml, SPACE_R, slice%bandpp, 1)

    call xg_init(XAB,slice%space,nrows,3*ncols,slice%spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_setBlock(XAB, xAXColsRows, nrows, ncols)
    call xg_setBlock(XAB, xBXColsRows, nrows, ncols, fcol=ncols+1)
    call xg_setBlock(XAB, xRXColsRows, nrows, ncols, fcol=2*ncols+1)
    
    ! Compute A*Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
    call getAX_BX(xXColsRows, xAXColsRows, xBXColsRows)
    call xgBlock_zero_im_g0(xAXColsRows)
    call xgBlock_zero_im_g0(xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ! <Psi|H|Psi>
    call xgBlock_colwiseDotProduct(xXColsRows, xAXColsRows, Results1%self, comm_loc=xmpi_comm_null)

    ! <Psi|S|Psi>
    call xgBlock_colwiseDotProduct(xXColsRows, xBXColsRows, Results2%self, comm_loc=xmpi_comm_null)

    ! eig = <Psi|H|Psi> / <Psi|S|Psi>
    call xgBlock_colwiseDivision(Results1%self, Results2%self, eigen)

    ! xRXColsRows = S|Psi>
    call xgBlock_copy(xBXColsRows,xRXColsRows)
    ! xRXColsRows = - eigen * S|Psi>
    call xgBlock_ymax(xRXColsRows,eigen,0,1)
    ! norml = |eigen*S|Psi>|^2
    call xgBlock_colwiseNorm2(xBXColsRows, norml%self, comm_loc=xmpi_comm_null)
    ! xRXColsRows = H|Psi> - eigen * S|Psi>
    call xgBlock_add(xRXColsRows,xAXColsRows)
    ! resid = |xRXColsRows|^2
    call xgBlock_colwiseNorm2(xRXColsRows, resid, comm_loc=xmpi_comm_null)

    ! compute relative resid=resid/(|eigen|*|BX|)^2
    ! scales with the eigenvalue
    call xgBlock_colwiseDivision(resid, norml%self, resid)

    call xg_free(Results1)
    call xg_free(Results2)
    call xg_free(XAB)
    call xg_free(norml)

end subroutine slice_queryCandidates
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_unitTest
!! NAME
!! slice_unitTest
!! 
!! FUNCTION
!! Given a set of vectors X computes the columnwise norm squared, on:
!! - CPU sequential,
!! - CPU MPI,
!! - GPU sequential (if GPU available),
!! - GPU MPI (if GPU available).
!! Checks that all four are the same.
!! 
!! SOURCE

function slice_unitTest(X) result(ierr)

    implicit none

    !Arguments ------------------------------------
    type(xgBlock_t), intent(in) :: X
    integer                     :: ierr         
    
    !Local variables-------------------------------
    integer :: nrows,ncols
    integer :: comm
    integer :: gpu_option
    real(dp) :: id_cpu, id_from_gpu, id_to_gpu

    ! *********************************************************************

    ! Get nrows, nrcols, gpu_option of X
    call xgBlock_getSize(X,nrows,ncols)
    call xgBlock_get_gpu_option(X,gpu_option)
    call xgBlock_get_communicator(X,comm)
    ! above this is just comm(X)
    write(std_out,*) '<--*--> Start unitary test for matrix'
    write(std_out,*) 'nrows       =', nrows
    write(std_out,*) 'ncols       =', ncols
    write(std_out,*) 'gpu_option  =', gpu_option
    write(std_out,*) 'communicator=', comm

    id_cpu = -1.0d0
    id_from_gpu = -1.0d0
    id_to_gpu = -1.0d0

    if (gpu_option==ABI_GPU_DISABLED) then
        id_cpu = xgBlock_getid(X, comm)
        write(std_out,*) 'id_cpu     =', id_cpu
    end if        

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    if (gpu_option/=ABI_GPU_DISABLED) then
        call xgBlock_copy_from_gpu(X)
        id_from_gpu = xgBlock_getid(X, comm)
        call xgBlock_copy_to_gpu(X)
        write(std_out,*) 'id_from_gpu=', id_from_gpu
    else
        call xgBlock_copy_to_gpu(X)
        call xgBlock_copy_from_gpu(X)
        id_to_gpu = xgBlock_getid(X, comm)
        call xgBlock_copy_to_gpu(X)
        call xgBlock_copy_from_gpu(X)
        write(std_out,*) 'id_to_gpu  =', id_to_gpu
    end if
#endif
    
    ! TODO getid on xmpi_comm_null
    ! then sum across comm
    ! should give the same as
    ! getid on comm
    ! this is useful for MPI

    ierr = 0
    if (id_cpu < 0) ierr = -1
    if (id_from_gpu < 0) ierr = -1
    if (id_to_gpu < 0) ierr = -1

end function slice_unitTest
!!***

!----------------------------------------------------------------------

!!****f* m_slice/splitSpectrum
!! NAME
!! splitSpectrum
!! 
!! FUNCTION
!! Split working spectrum [a,b) mapped to [-1,1) with center and radius.
!! 
!! NOTES
!! Working spectrum is initialized as [a,b) and assumes that b may have an error
!! so adapts it. A way to incorporate case of large error on b is to increase b 
!! by a step until all the mass equal to nband is included. 
!! While missing mass then extend
!! 
!! SOURCE

subroutine splitSpectrum(nband_tot, nstep_bisect, center, radius, a, b, c_split, &
        cheby_moments, comm)

        implicit none

        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: c_split
        real(dp), intent(in) :: center, radius
        integer, intent(in) :: nstep_bisect
        integer, intent(in) :: nband_tot
        complex(dp), intent(in) :: cheby_moments(:,:)
        integer, intent(in), optional :: comm

        integer :: neigenpairs
        integer :: ndeg_filter
        integer :: ishift
        integer :: iext_step
        integer :: comm_, ierr
        real(dp) :: mass_diff, mass_diff_prev
        real(dp) :: tot_mass, mass_left, mass_right
        real(dp) :: b_ext, c, width, width_ext
        real(dp), allocatable :: spectral_mass_left(:)
        real(dp), allocatable :: spectral_mass_right(:)
        real(dp), allocatable :: energy_per_band(:)
        real(dp), allocatable :: cja(:)

        ! *********************************************************************

        ! todo might be a good idea to put here what there is in 
        ! computeTraceEstimation

end subroutine splitSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_applySelectionWindow
!! NAME
!! slice_applySelectionWindow
!!
!! FUNCTION
!! Input: k=m+p ou m la masse et p oversampling parameter
!! step 1) calculer les valeurs de Ritz
!! step 2) calculer les résidus
!! step 3) mettre les valeurs de Ritz dans part_low_bounds, part_upp_bounds sous condition que res<tol
!! step 4) sinon faire un sketch de taille K du reste
!!
!! INPUT
!! safe is between 2 and 10
!!
!! OUTPUT
!! mapper contains the column indices per slice (allowing repetitions)
!! 
!! SOURCE

  subroutine slice_applySelectionWindow(slice, mapper, eigen, resid, tol, safe)

    implicit none

    type(slice_t), intent(inout) :: slice
    integer, pointer, intent(inout) :: mapper(:,:)
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: resid
    real(dp), intent(in) :: tol
    integer, intent(in) :: safe

    integer :: nband, iband, islice, m, k_kept
    real(dp) :: uppb, lowb, uppb_plus, lowb_minus
    real(dp) :: theta, relaxf, relres
    real(dp), pointer :: resid_vals(:,:) => null()
    real(dp), pointer :: thetas(:,:) => null()

  ! *********************************************************************

    nband = slice%neigenpairs

    call xgBlock_reverseMap(eigen, thetas, rows=nband, cols=1)
    call xgBlock_reverseMap(resid, resid_vals, rows=nband, cols=1)

    ! apply selection window per slice
    do islice=1,slice%nslice
        m = slice%neigenpairs_per_slice(islice)
        lowb = slice%low_bounds(islice)
        uppb = slice%upp_bounds(islice)

        ! keep candidates in expanded interval
        k_kept = 0
        do iband=1, nband
            relres = resid_vals(iband, 1)
            if (relres < tol) then
                ! Soft selection:
                ! relaxed bounds for converged eigenpairs near boundaries
                relaxf = safe * relres
                lowb_minus = lowb - relaxf
                uppb_plus = uppb + relaxf
                theta = thetas(iband, 1)
            
                ! accept if lambda in relaxed interval
                if (theta < uppb_plus .and. theta > lowb_minus) then
                    mapper(k_kept, islice) = iband
                    k_kept = k_kept + 1
                end if
            end if
        end do

        write(std_out,*) 'selection window kept', k_kept, 'for eigendimension', m
        flush(std_out)
    end do 
   
  end subroutine slice_applySelectionWindow
!!***

end module m_slice
!!***

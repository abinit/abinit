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
!! Dependence with other modules and hierarchy in this directory:
!!   m_slice uses:
!!      |- m_polynomial_filter | various computational routines
!!      |- m_trace_estimation  | various computational routines
!!      |- m_chebfi2           ! the Chebyshev recursion
!!      |- m_slice_task        | Low-level
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
    use m_slice_task
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
        real(dp) :: mineig                                  ! slightly larger than mineig_global
        real(dp) :: mineig_global                           ! guaranteed lower spectral bound
        real(dp) :: maxeig_global                           ! guaranteed upper spectral bound

        ! Arrays related to polynomial filtering
        integer, allocatable :: neigenpairs_per_slice(:)    ! number of total eigenpairs per slice
        integer, allocatable :: poly_degrees(:)             ! polynomial filter degrees
        real(dp), allocatable :: low_bounds(:)              ! lower bounds in spectral partition (disjoint)
        real(dp), allocatable :: upp_bounds(:)              ! upper bounds in spectral partition (disjoint) 

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
    if (present(gpu_kokkos_nthrd)) then
        slice%gpu_kokkos_nthrd = gpu_kokkos_nthrd
    end if
    slice%gpu_thread_limit = 0
    if (present(gpu_thread_limit)) then
        slice%gpu_thread_limit = gpu_thread_limit
    end if

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
!! X0=           size (spacedim, neigenpairs)
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

subroutine slice_run(slice, X0, getAX_BX, getBm1X, eigen, residu, nspinor)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X0     ! size (spacedim, neigenpairs)
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
    integer :: i, iband, nrows, ncols
    integer :: space_res
    integer :: neigenpairs, bandpp, ndeg_filter
    integer :: neigenpairs_ext
    integer :: comm, comm_rows, comm_cols
    integer :: num_kept
    integer :: k_conv, ndeg_filter_max
    integer :: ierr
    real(dp) :: theta
    real(dp) :: safe, tol ! for slice selection window
    real(dp) :: a_part, b_part, max_resid_kept
    logical :: is_lowpass, on_host, on_device
    ! todo use slice%..
    ! IML ----> variables for logic of slice solver
    type(matrixInfo_t) :: matrixInfo
    type(xg_t) :: DivResults
    type(activeTask_t) :: task
    type(taskScheduler_t) :: scheduler
    type(asyncMemory_t) :: asyncMemory
    ! <----- IML
    ! Arrays
    real(dp) :: tsec(2)
    real(dp), allocatable :: moments(:)
    logical, allocatable, target :: mapper(:,:)
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
    ncols = slice%neigenpairs
    if (slice%paral_kgb==1) then ! colsrows distribution
        nrows = slice%total_spacedim
        ncols = slice%bandpp
    end if

    if (slice%space==SPACE_C) then
        space_res = SPACE_C
    else if (slice%space==SPACE_CR) then
        space_res = SPACE_R
    else
        ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
    end if

    call xg_init(DivResults, space_res, slice%neigenpairs, 1, gpu_option=slice%gpu_option)

    call init_matrixInfo(matrixInfo, slice%comm_rows, slice%comm_cols, slice%spacecom, slice%neigenpairs,&
        slice%total_spacedim, slice%spacedim, slice%space, slice%gpu_kokkos_nthrd, slice%gpu_thread_limit,&
        slice%gpu_option, slice%paral_kgb, slice%me_g0, slice%me_g0_fft)

    write(std_out,*) 'in slice_run'; flush(std_out)

    ! ================================== Prepare spectral slices ===========================================

    call timab(tim_slice_sched,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_SCHEDULE)

    ndeg_filter_max = 10 ! todo IML autotune based on the presence of oscillations
    ABI_MALLOC(moments, (ndeg_filter_max+1))

    ! Split spectrum and query X0
    ABI_NVTX_START_RANGE(NVTX_SLICE_RRQ)
    call slice_prepareSpectrum(slice, X0, matrixInfo, moments, DivResults%self, residu, getAX_BX, getBm1X, nspinor)
    ABI_NVTX_END_RANGE()

    ABI_NVTX_END_RANGE()
    call timab(tim_slice_sched,2,tsec)

    call slice_splitSpectrum(slice, moments)
    ABI_FREE(moments)

    ! Attribute column vectors of X0 to slices using query results
    ABI_MALLOC(mapper, (slice%neigenpairs, slice%nslice))
    mapper = .false.
    write(std_out,*) 'tolerance in residual window=', slice%tolerance
    flush(std_out)

    safe = 2.d0
    tol = 1e-2 ! slice%tolerance
    call slice_applySelectionWindow(slice, mapper, DivResults%self, residu, tol, safe)

    write(std_out,*) 'created the following map slices to cols='
    write(std_out,*) mapper(1,:)
    flush(std_out)

    ! Include m+p where p offset 10%
    slice%neigenpairs_per_slice = ceiling(slice%neigenpairs_per_slice*1.1d0)
    neigenpairs_ext = sum(slice%neigenpairs_per_slice)

    write(std_out,*) 'allocating async memory of size', neigenpairs_ext
    flush(std_out)

    ! Allocate extended buffer in Linalg representation. Notice spacecom communicator (global)
    call slice_task_allocateAsyncMemory(asyncMemory, matrixInfo, neigenpairs_ext)

    ! Initilize extended buffer with vectors from X or random vectors. This part also uses global comm.
    call slice_task_initAsyncMemory(asyncMemory, X0, matrixInfo, mapper, slice%neigenpairs_per_slice)

    write(std_out,*) 'getid after init', xgBlock_getid(asyncMemory%XextLinalg)
    flush(std_out)

    ! MPI phase I: Compute 'process-to-slices' distribution according to paral options
    call slice_task_initSchedule(scheduler, slice%paral_kgb, slice%paral_slice)

    ! MPI phase II: Compute 'slice columns-to-process distribution

    !! todo add in function
    !! Here we compute the distribution of columns in colsrows representation.

    ! Use scheduler to distribute the extended memory to processes



    ! At the end we should compute the distribution of columns in colsrows representation.
    ! ** This phase depends on the execution of slices **
    ! - If enable_paral then we should first divide processes to slices then divide slice bands to processes.
    ! - If disable_paral then we should use all available processes for every slice then divice slice bands to processes.
    ! Notice that the last step where we divide slice bands to number of processes is the same.
    ! Only the number of available processes to use per slice changes.

    ! This should allow to mark every process with the color of the slice.

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
!        ! Allocate and fill async memory buffer
!        nband_ext = sum(slice%neigenpairs_per_slice)
!        call allocate_extended_memory(asyncMemory, nband_ext)
!
!        ! k is the number of columbs from X0
!        ! p is the offset
!        call init_extended_memory(asyncMemory, task, X0, k, p)
!
!    else
!        ! configure task without any communicators
!    end if
!
!    call allocate_active_task(slice, scheduler, task)
!
!  
!    ! will use mapper to copy columns of X0 to async memory
!    !! assumes linalg representation of both ext and spectrum mem
!    call init_extended_memory(asyncMemory, X0, mapper)
!    
!    ! ============================ Active task execution =======================================
!
!    call timab(tim_slice_me,1,tsec)
!    
!    if (slice%paral_kgb==0 .or. slice%paral_slice==DISABLE_PARAL) then
!        ! execute active tasks sequentially
!
!        do islice=1, nslice
!            if (task%active(islice)) then                   ! <--- normally it should not be that different.. 
!            call schedule_next_task(scheduler, neigenpairs)
!            call init_active_memory(asyncMemory, task, X0, p)
!            call execute_active_task(task)
!            call mask_active_task(task, tol) ! mask extendedMem
!        end do
!       
!    else
!        ! execute active tasks in parallel
!
!        call init_active_memory(asyncMemory, task, X0, p)
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
!    call free_extended_memory(asyncMemory)
!    ABI_FREE(mapper)
    call xg_free(DivResults)

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

subroutine slice_prepareSpectrum(slice, X, matrixInfo, moments, eigen, resid, getAX_BX, getBm1X, nspinor)

    implicit none

    ! Arguments
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X
    type(matrixInfo_t), intent(inout) :: matrixInfo
    real(dp), intent(inout) :: moments(:)
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
    type(xgTransposer_t) :: xgTransposerX
    ! Arrays
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
        work_size = slice%neigenpairs ! fixme will complain
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
    m_probe = 5
    ndeg_filter_max = size(moments) - 1
    write(std_out,*) 'Here I compute Stochastic Trace Estimation'
    write(std_out,*) 'm_probe=    ', m_probe
    write(std_out,*) 'ndeg_filter=', ndeg_filter_max
    flush(std_out)

    call computeTraceEstimation(matrixInfo, slice%nslice, slice%ecut, slice%paw, &
        slice%tolerance, getAX_BX, getBm1X, ndeg_filter_max, m_probe, lanczos_lowb_global, moments)
   
    ! Define working spectrum to be splitted to slices
    slice%mineig = lanczos_lowb 
    slice%mineig_global = lanczos_lowb_global
    slice%maxeig_global = slice%ecut

    write(std_out,*) 'STE exited'
    flush(std_out)
 
    ! Compute Rayleigh quotients (colsrows distribution)
    !ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RRQ)
    !call timab(tim_RR_q, 1, tsec)

    if (slice%paral_kgb==1) then
        my_rank = xmpi_comm_rank(slice%spacecom)
        my_shift = my_rank * slice%bandpp
        call xgBlock_reshape(eigen, 1, slice%neigenpairs)
        call xgBlock_reshape(resid, 1, slice%neigenpairs)
        call xgBlock_setBlock(eigen, eigen_me, 1, slice%bandpp, fcol=my_shift+1)
        call xgBlock_setBlock(resid, resid_me, 1, slice%bandpp, fcol=my_shift+1)
        call xgBlock_reshape(eigen, slice%neigenpairs, 1)
        call xgBlock_reshape(resid, slice%neigenpairs, 1)
        call xgBlock_reshape(eigen_me, slice%bandpp, 1)
        call xgBlock_reshape(resid_me, slice%bandpp, 1)
    else
        call xgBlock_setBlock(eigen, eigen_me, rows(eigen), 1)
        call xgBlock_setBlock(resid, resid_me, rows(resid), 1)
    end if

    call slice_queryCandidates(slice, xXColsRows, eigen_me, resid_me, getAX_BX)
    !call timab(tim_RR_q, 2, tsec)
    !ABI_NVTX_END_RANGE()

    call xgBlock_mpi_sum(eigen, comm=slice%spacecom)
    call xgBlock_mpi_sum(resid, comm=slice%spacecom)

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

    ncols = slice%bandpp
    nrows = rows(xXColsRows)
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

!!****f* m_slice/slice_splitSpectrum
!! NAME
!! slice_splitSpectrum
!! 
!! FUNCTION
!! Split working spectrum [a,b) mapped to [-1,1) with center and radius.
!!
!! INPUT
!! a= slice%mineig_global
!! b= slice%maxeig_global 
!! 
!! SIDE EFFECTS
!! slice%low_bounds
!! slice%upp_bounds
!! slice%neigenpairs_per_slice
!! 
!! SOURCE

  subroutine slice_splitSpectrum(slice, moments)

      implicit none

      ! Arguments ------------------------------------
      type(slice_t), intent(inout) :: slice
      real(dp), intent(inout) :: moments(:)

      ! Local variables --------------------------------
      integer :: uppb_loc, ib
      integer :: ngrid_fine, ngrid_coarse
      integer :: half_index, i_left, i_right
      integer :: i_split, num_moments
      logical :: gap_left, gap_right
      real(dp) :: a, b
      real(dp) :: center, radius
      real(dp) :: step_fine, step_coarse
      real(dp) :: partial_mass, half_mass
      real(dp) :: b_scaled
      real(dp) :: mass_left, mass_right
      real(dp) :: lambda_plus_wanted
      real(dp), allocatable :: cumm_eigen_count(:)
      real(dp), allocatable :: bgrid_fine(:)
      real(dp), allocatable :: bgrid_coarse(:)
      real(dp), allocatable :: work(:)
      character(len=500) :: msg

    ! *********************************************************************

      if (slice%nslice/=2) then
          ABI_ERROR("spectral bidirectional split not yet implemented for 3 slices")
      end if

      ngrid_coarse = 10 ! coarse, just to find uppb, hardcoded
      ngrid_fine = 30 ! used for cumulative eigenvalue count, hardcoded
    
      a = slice%mineig_global
      b = slice%maxeig_global

      num_moments = size(moments) 
      ABI_MALLOC(work, (num_moments))
      ABI_MALLOC(bgrid_coarse, (ngrid_coarse))
      ABI_MALLOC(bgrid_fine, (ngrid_fine))
      ABI_MALLOC(cumm_eigen_count, (ngrid_fine))

      center = (a + b) / 2.d0
      radius = (b - a) / 2.d0 

      ! #########################################
      ! ########## Coarse resolution ############
      ! #########################################

      ! scan with lowpass
      ! first pass uppb is actually unknown
      step_coarse = (b - a) / (ngrid_coarse - 1)
      bgrid_coarse = (/ ( a + (ib-1)*step_coarse, ib=1,ngrid_coarse ) /)
      uppb_loc = -1
      do ib=1, ngrid_coarse
          b_scaled = (bgrid_coarse(ib) - center) / radius
          partial_mass = get_eigenvalue_count(b_scaled, moments, work) 

          write(std_out,*) ib, 'scan: <=', bgrid_coarse(ib), 'mass=', partial_mass 
          flush(std_out)

          if (partial_mass > slice%neigenpairs) then
              uppb_loc = ib
              exit
          end if

      end do

      lambda_plus_wanted = (bgrid_coarse(uppb_loc-1) + bgrid_coarse(uppb_loc)) / 2.d0
      write(std_out,*) 'found upp bound in', lambda_plus_wanted
      write(std_out,*) 'estimated mass=', get_eigenvalue_count((lambda_plus_wanted-center)/radius, moments, work)
      write(std_out,*) 'starting adaptive refinement ..'
      flush(std_out)

      ! #########################################
      ! ########### Fine resolution #############
      ! #########################################

      ! Now compute eigenvalue count
      step_fine = (lambda_plus_wanted - a) / (ngrid_fine - 1)
      bgrid_fine = (/ (a + (ib-1)*step_fine, ib=1,ngrid_fine) /) 
      do ib=1, ngrid_fine
          b_scaled = (bgrid_fine(ib) - center) / radius
          partial_mass = get_eigenvalue_count(b_scaled, moments, work)
          cumm_eigen_count(ib) = partial_mass
          if (partial_mass < 0.d0 .and. abs(partial_mass) > 1.d0) then
              write(msg,'(a)') "fine resolution for eigenvalue count failed due to oscillations. ",&
                  "Reduce ndeg_filter_max (hard-coded) as a solution"
              ABI_ERROR(msg)
          end if
          write(std_out,*) ib, 'scan: <=', bgrid_fine(ib), 'mass=', partial_mass
      end do
      !! todo if not found then add points in the fine grid...

      ! Prepare: detect gap existence in the interior of slice
      ! define value of smallest_gap
      ! todo 

      ! step 1 cut in half balanced mass
      half_mass = cumm_eigen_count(ngrid_fine)/2.d0 
      half_index = minloc(abs(half_mass - cumm_eigen_count), dim=1)
      ! step 2 adjust so that it is on constant mass (predicts gap)
      ! bidirectional search left and right
      write(std_out,*) 'bidirectional search from index=', half_index
      write(std_out,*) 'half mass=', half_mass
      flush(std_out)
      i_left = half_index
      i_right = half_index
      gap_left = .false.
      gap_right = .false.
      do while((.not.gap_right .and. .not.gap_left) .and. (i_left >= 2 .and. i_right <=ngrid_fine-1))
          ! this is if gap exists. If it does not exist.. must minimize using smallest_gap 
          gap_left = abs(cumm_eigen_count(i_left) - cumm_eigen_count(i_left-1)) < 1e-4
          gap_right = abs(cumm_eigen_count(i_right) - cumm_eigen_count(i_right+1)) < 1e-4
          i_left = i_left - 1
          i_right = i_right + 1
      end do
      if (gap_right) then
          i_split = i_right
      end if
      if (gap_left) then
          i_split = i_left
      end if
      if (.not.gap_right .and. .not.gap_left) then
          ! todo treat this case
          ABI_ERROR("spectrum has no gap..")
      end if
      mass_left = get_eigenvalue_count((bgrid_fine(i_split)-center)/radius,moments,work)
      mass_right = slice%neigenpairs - mass_left
      write(std_out,*) 'i_split val=', i_split, bgrid_fine(i_split)
      write(std_out,*) 'mass left=', mass_left
      write(std_out,*) 'mass right=', mass_right
      flush(std_out)

      ! Store results into the slice
      slice%neigenpairs_per_slice(1) = ceiling(mass_left)
      slice%neigenpairs_per_slice(2) = ceiling(mass_right)

      slice%low_bounds(1) = slice%mineig
      slice%low_bounds(2) = bgrid_fine(i_split)

      slice%upp_bounds(1) = bgrid_fine(i_split)
      slice%upp_bounds(2) = bgrid_fine(ngrid_fine) 
    
      ABI_FREE(bgrid_fine)
      ABI_FREE(bgrid_coarse)
      ABI_FREE(cumm_eigen_count)
      ABI_FREE(work)

  end subroutine slice_splitSpectrum
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
    logical, intent(inout) :: mapper(:,:)
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: resid
    real(dp), intent(in) :: tol
    real(dp), intent(in) :: safe

    integer :: nband, iband, islice, m, k_kept, ierr
    real(dp) :: uppb, lowb, uppb_plus, lowb_minus
    real(dp) :: theta, relaxf, relres
    real(dp), pointer :: resid_vals(:,:) => null()
    complex(dp), pointer :: thetas(:,:) => null()

  ! *********************************************************************

    nband = slice%neigenpairs

    call xgBlock_reverseMap(eigen, thetas, rows=nband, cols=1)
    call xgBlock_reverseMap(resid, resid_vals, rows=nband, cols=1)
   
    ! apply selection window per slice
    do islice=1,slice%nslice
        m = slice%neigenpairs_per_slice(islice)
        lowb = slice%low_bounds(islice)
        uppb = slice%upp_bounds(islice)
        write(std_out,*) islice, 'slice bounds', lowb, uppb

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
                theta = real(thetas(iband, 1))
                write(std_out,*) iband, 'theta=', theta, 'resid=', relres
                flush(std_out)
            
                ! accept if lambda in relaxed interval
                if (theta < uppb_plus .and. theta > lowb_minus) then
                    k_kept = k_kept + 1
                    mapper(iband, islice) = .true.
                end if
            end if
        end do

        ! todo m not set must modify computeTraceEstimation and separate the split spectrum from that
        ! because split spectrum will store info on slice object that is invisible from trace
        write(std_out,*) 'selection window kept', k_kept, 'for eigendimension', m
        flush(std_out)
    end do 
   
  end subroutine slice_applySelectionWindow
!!***

end module m_slice
!!***

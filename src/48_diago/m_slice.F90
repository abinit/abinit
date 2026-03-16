!!****f* ABINIT/m_slice
!! NAME
!! m_slice
!!
!! FUNCTION
!! This module contains the types and routines used to apply the Spectrum Slicing 
!! method. It mainly defines 'slice' datatypes and associated methods. 
!!
!! NOTES
!! Main features:
!! - uses 'xgTools' implementation for matrix data structure.
!! - adopts 'chebfi' functionalities for most matrix calculations.
!! - uses polynomial filtering, Chebyshev lowpass, next Chebyshev-Jackson.
!! - implements scheduler and resource allocator for slice tasks.
!! - applies Rayleigh-Ritz for individual slices in parallel.
!! 
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
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (IL)
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
    use m_slice_cprj, only: smallestTridiagEigenpair, buildChebyshevJacksonCoeffs

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

    ! Spectral interval decomposition strategy (spectral_cut option)
    !---------------------------------------------------
    integer, parameter :: DIVIDE_INTERVAL_WIDTH    = 0
    integer, parameter :: DIVIDE_NUMBER_OF_VECTORS = 1
    integer, parameter :: DIVIDE_SPECTRAL_GAPS     = 2

    ! Load balance criteria for fair resource allocation (paral_slice option)
    !---------------------------------------------------
    integer, parameter :: FAIR_BANDPP      = 1              ! bandpp (unweigted)
    integer, parameter :: FAIR_BANDPP_WDEG = 2              ! bandpp weighted by degree
    integer, parameter :: EVEN_SLICES      = 3              ! all slices have the same num of procs 

    ! Timers
    !---------------------------------------------------
    integer, parameter :: tim_swap      = 1761
    integer, parameter :: tim_RR_q      = 1759
    integer, parameter :: tim_barrier   = 1764
    integer, parameter :: tim_copy      = 1765
    integer, parameter :: tim_Bortho_X  = 1641
    integer, parameter :: tim_getAX_BX  = 1754
    integer, parameter :: tim_invovl    = 1755
    integer, parameter :: tim_lanczos   = 2167
    integer, parameter :: tim_trace     = 2168
    
    ! Public 'slice' datatype
    !-------------------------------------------------
    type, public :: slice_t

        ! MPI-related (me=current MPI process)
        integer :: nproc                                    ! number of available processes
        integer :: comm_rows                                ! xmpi_comm_self ...
        integer :: comm_cols                                ! same as spacecom
        integer :: spacecom                                 ! same as comm_cols
        integer :: me_g0                                    ! process contains G(0,0,0) 
        integer :: me_g0_fft                                ! process contains G(0,0,0) for fft
        integer :: me_nproc_slice                           ! number of processes reserved to slice in use
        integer :: me_comm_slice                            ! sub-communicator reserved to slice in use
        integer :: me_comm_rows                             ! sub-communicator for rows in Transposer
        integer :: me_comm_cols                             ! sub-communicator for cols in Transposer
        integer :: me_id_slice                              ! identifier in 1,..,nslice of slice task in use
        integer :: me_neigenpairs_slice                     ! total number of eigenpairs of slice in use
        integer :: me_bandpp_slice                          ! number of distributed bands per process for slice in use
        integer :: me_ndeg_slice                            ! polynomial filter degree for slice in use
        
        ! Eigenpair parameters
        integer :: neigenpairs                              ! total number of bands (=number of eigenpairs)
        integer :: neigenpairs_ext                          ! total numner of extended columns
        integer :: total_spacedim                           ! total number of plane-waves
        integer :: bandpp                                   ! nb of bands per process in colsrows representation 
        integer :: spacedim                                 ! nb of plane-waves per process in linalg representation
        integer :: nslice                                   ! number of spectral slices
        integer :: space                                    ! real or complex eigenvectors
        integer :: nbdbuf
        
        ! GPU-related
        integer :: gpu_kokkos_nthrd                 
        integer :: gpu_thread_limit
        
        ! Flags
        logical :: on_host = .false.                        ! running on CPU
        logical :: on_device = .false.                      ! running on GPU
        logical :: use_linalg = .false.                     ! use linalg representation
        logical :: use_colsrows = .false.                   ! use colsrows representation
        
        ! Options
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

        ! Memory buffers
        type(xg_t) :: X_ext                                 ! eigenvector memory used by all slices
        type(xgTransposer_t) :: xgTransposerXext            ! transposer datastructure 
        
        ! Pointers
        type(xgBlock_t) :: me_Xext_active                   ! eigenvector memory in use by active slice
        type(xgBlock_t) :: XextLinalg

        ! Arrays for my slice only
        integer, allocatable :: me_ncolsColsRows_slice(:)   ! ncol of colsrows representation for my slice
        integer, allocatable :: me_nrowsLinalg_slice(:)     ! nrow of linalg representation for my slice

        ! Arrays related to MPI (all slices)
        integer, allocatable :: neigenpairs_per_slice(:)    ! number of total eigenpairs per slice
        integer, allocatable :: nproc_per_slice(:)          ! number of processes per slice
        integer, allocatable :: lookup_proc(:)              ! which slice each process serves
        integer, allocatable :: ncolsColsRows(:)            ! ncol of colsrows representation for global X

        ! Arrays related to data distribution (all slices)
        integer, allocatable :: fcol_in_X(:)                ! first band of slice in spectrum memory
        integer, allocatable :: fcol_in_Xext(:)             ! first band of slice in extended memory

        ! Arrays related to polynomial filtering (all slices)
        integer, allocatable :: poly_degrees(:)            ! polynomial filter degrees
        real(dp), allocatable :: part_low_bounds(:)         ! lower bounds in spectral partition (disjoint)
        real(dp), allocatable :: part_upp_bounds(:)         ! upper bounds in spectral partition (disjoint) 
        real(dp), allocatable :: poly_low_bounds(:)         ! lower bounds used to define polynomials (overlap)
        real(dp), allocatable :: poly_upp_bounds(:)         ! upper bounds used to define polynomials (overlap)

    end type slice_t

    ! Public methods associated to 'slice' datatype
    !-------------------------------------------------
    public :: slice_init                                    ! initialize slice datatype object
    public :: slice_free                                    ! free slice datatype object
    public :: slice_allschedule                             ! build slice distributed workspace
    public :: slice_run                                     ! run Spectrum Slicing for active slice task
    public :: slice_allmerge                                ! merge slice result to distributed workspace
    public :: slice_unitTest                                ! used for debugging GPU device

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

    ! Set flags to initial state (before Transpose)
    slice%use_linalg = .true.
    slice%use_colsrows = .false.
    call slice_queryHostDevice(slice,on_host,on_device)

    write(std_out,*) 'At slice_init:'
    write(std_out,*) 'slice%spacecom id=', slice%spacecom 
    write(std_out,*) 'slice%spacecom size=', xmpi_comm_size(slice%spacecom)
    write(std_out,*) 'on_host=, on_device', on_host, on_device
    write(std_out,*) 'slice%spacedim=', slice%spacedim
    write(std_out,*) 'slice%total_spacedim=', slice%total_spacedim
    write(std_out,*) 'slice%bandpp=', slice%bandpp

    !if (xmpi_comm_size(slice%spacecom)==1) then
    !    ABI_ERROR("Slicing with 1 MPI process not implemented")
    !end if

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
    ABI_MALLOC_IFNOT(slice%nproc_per_slice, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%lookup_proc, (slice%nproc))

    ABI_MALLOC_IFNOT(slice%fcol_in_X, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%fcol_in_Xext, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%ncolsColsRows, (slice%nproc))

    ABI_MALLOC_IFNOT(slice%poly_degrees, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%part_low_bounds, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%part_upp_bounds, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%poly_low_bounds, (slice%nslice))
    ABI_MALLOC_IFNOT(slice%poly_upp_bounds, (slice%nslice))

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_free
!! name
!! slice_free

subroutine slice_free(slice)

    implicit none
    
    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice

    ! *********************************************************************

    call xg_free(slice%X_ext)
    call xgTransposer_free(slice%xgTransposerXext)

    ABI_SFREE(slice%me_ncolsColsRows_slice)   
    ABI_SFREE(slice%me_nrowsLinalg_slice)

    ABI_SFREE(slice%neigenpairs_per_slice)
    ABI_SFREE(slice%nproc_per_slice)
    ABI_SFREE(slice%lookup_proc)
    ABI_SFREE(slice%ncolsColsRows)   

    ABI_SFREE(slice%fcol_in_X)
    ABI_SFREE(slice%fcol_in_Xext)

    ABI_SFREE(slice%poly_degrees)
    ABI_SFREE(slice%part_low_bounds)
    ABI_SFREE(slice%part_upp_bounds)
    ABI_SFREE(slice%poly_low_bounds)
    ABI_SFREE(slice%poly_upp_bounds)

    !if (slice%me_comm_slice /= slice%spacecom) then
    !    call xmpi_comm_free(slice%me_comm_slice)
    !end if
    !if (slice%me_comm_rows /= slice%comm_rows) then
    !    call xmpi_comm_free(slice%me_comm_rows)
    !end if
    !if (slice%me_comm_cols /= slice%comm_cols) then
    !    call xmpi_comm_free(slice%me_comm_cols)
    !end if

end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allschedule
!! NAME
!! slice_allschedule
!! 
!! FUNCTION
!! Resource allocation for individual slice tasks.
!! Create extended memory buffers and distribute them according
!! to slice logic. Essentially allocates memory slice%me_Xext_active 
!! favoring data overlap over communication overlap.
!!
!! INPUT 
!! X0                =eigenvector guess used to split spectrum
!!        not suitable for parallel Rayleigh-Ritz calculations
!! 
!! OUTPUT
!! eigen= Rayleigh quotients associated to X0 of size (1,neigenpairs)
!! 
!! SIDE EFFECTS
!! slice%me_Xext_active  = adapted version distributed correctly
!!        and free of data overlap. Suitable for parallel RR.
!! 
!! SOURCE

subroutine slice_allschedule(slice, X0, getAX_BX, getBm1X, eigen, nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), target, intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X0
    type(xgBlock_t), intent(inout) :: eigen
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

    ! Local variables --------------------------------
    integer, parameter :: tim_slice_sched = 2161
    integer :: neigenpairs, iband, min_loc, islice
    integer :: fcol, fcol_ext, ncols, itest, iparal
    integer :: npband_test
    integer :: wanted_mass
    logical :: on_host, on_device
    real(dp) :: lowb, uppb, c_split
    real(dp) :: lambda_minus, lambda_plus
    real(dp) :: tol12 = 1.0e-12
    ! Derived types
    type(xg_t) :: resid0
    type(xgBlock_t) :: eigen_sorted
    type(xgBlock_t) :: slicecols_in
    type(xgBlock_t) :: slicecols_ext_out
    ! Arrays
    integer :: npband_list(4)
    real(dp) :: tsec(2)
    integer, allocatable :: weights(:)
    integer, allocatable :: npband_per_slice(:)
    integer, allocatable :: nband_per_slice(:)
    integer, allocatable :: bands_left(:)
    integer, allocatable :: bands_right(:)
    
    ! *********************************************************************

    call timab(tim_slice_sched,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_SCHEDULE)
    
    ! Sanity check for X0 (on GPU): verify we are in target enter data map
    ! TODO debug not working
    !if (slice%gpu_option==ABI_GPU_OPENMP) then
    !    call slice_queryHostDevice(slice, on_host, on_device)
    !    if (on_host .or. (.not. on_device)) then
    !        ABI_ERROR("not in target data map region")
    !    end if
    !end if

    neigenpairs = slice%neigenpairs
  
    ! ===================== Compute Rayleigh quotients and residuals ===================================
    
    wanted_mass = ceiling(neigenpairs * 0.6d0) ! plus 10% extra vectors
    
    ABI_MALLOC(bands_left, (wanted_mass))
    ABI_MALLOC(bands_right, (wanted_mass)) 

    write(std_out,*) 'X0 (init)=', xgBlock_getid(X0); flush(std_out)
    
    ABI_NVTX_START_RANGE(NVTX_SLICE_RRQ)
    call slice_prepareSpectrum(slice, X0, lowb, uppb, c_split, bands_left, bands_right, &
        getAX_BX, getBm1X, nspinor)
    ABI_NVTX_END_RANGE()
    
    write(std_out,*) 'X0 (sketched)=', xgBlock_getid(X0); flush(std_out)
    
    !write(std_out,*) 'wanted mass=', wanted_mass
    !write(std_out,*) 'bands_left=', bands_left(:)
    !write(std_out,*) 'bands_right=', bands_right(:)
    !flush(std_out)

    ABI_FREE(bands_left)
    ABI_FREE(bands_right)

    !! TODO actually some bands are never assigned to a slice..
   
    ! Output:
    ! - c_split                 : between [a,b)
    ! - nvec_left, nvec_right   : trial vectors per slice
    ! - index_left, index_right : column indices per slice
    ! 

    ! ===================== Compute guaranteed spectral bounds ======================================== 

    slice%mineig_global = lowb
    slice%maxeig_global = slice%ecut

    ! ===================== Split interval [lambda_minus,lambda_plus) into slices ======================
   
    ! todo simplify fix polynomial degree and give it here

    ! Slice left
    !slice%neigenpairs_per_slice(1) = 180 ! number of TRUE eigenvalues in (poly_low, poly_upp)
    !slice%poly_degrees(1) = 10           ! filter degree
    !slice%part_low_bounds(1) = -0.14     ! first slice is lowpass so interval to suppress
    !slice%part_upp_bounds(1) = 2.1       ! used for convergence <----
    !slice%poly_low_bounds(1) = -0.14     ! with overlap
    !slice%poly_upp_bounds(1) = 2.5       ! with overlap

    ! Slice right
    !slice%neigenpairs_per_slice(2) = 180 ! number of TRUE eigenvalues in (poly_low, poly_upp)
    !slice%poly_degrees(2) = 50           ! slice%ndeg_filter
    !slice%part_low_bounds(2) = 2.1       ! used for convergence <-----
    !slice%part_upp_bounds(2) = 5.0       ! used for convergence <-----
    !slice%poly_low_bounds(2) = 1.9       ! with overlap
    !slice%poly_upp_bounds(2) = 5.2       ! with overlap

    ! todo deduce from vector pruning
    ! Indices of test vectors
    ! todo big modif col_in_X should be replaced by a simple index set
    slice%fcol_in_X(1) = 1
    slice%fcol_in_Xext(1) = 1
    
    slice%fcol_in_X(2) = 1
    slice%fcol_in_Xext(2) = slice%neigenpairs_per_slice(1) + 1

    if (slice%nslice==3) then
        slice%fcol_in_X(3) = 1
        slice%fcol_in_Xext(3) = slice%fcol_in_Xext(2) + slice%neigenpairs_per_slice(2) 
    end if

    ! Slice three
    !slice%neigenpairs_per_slice(3) = 96 ! wanted_mass
    !slice%fcol_in_X(3) = 1
    !slice%fcol_in_Xext(3) = 97
    !slice%poly_degrees(3) = 50! slice%ndeg_filter
    !slice%part_low_bounds(3) = 2.3
    !slice%part_upp_bounds(3) = 5.0
    !slice%poly_low_bounds(3) = 2.2        ! with overlap
    !slice%poly_upp_bounds(3) = 5.2        ! with overlap

    ! TODO analyze the lambda of Cheby converged in order to find out how many vectors to put per slice
    ! then fix true slices and see what happens
    ! You can try this with the small system alu as well. It has the mass splitted into two parts
    ! then we can compare the total time spent to filter etc

    ! ===================== Resource management system ================================================= 

    ! TODO half of the processes are assigned per slice, not true actually
    ! under the uniform mass splitting that simplifies things and we no longer need fair allocation
    ! Process: 
    ! 
    
    ! Divide resources into slice tasks
    call slice_allocateResources(slice)

    ! Run on all ranks of spacecom: Mark my slice task and resources as actively in use
    call slice_markActiveTask(slice)

    ! todo ok so normally here when tasks have been allocated we should sketch using G
    ! of size nxk where k is the estimated rank of the eigenspace. Then we also sketch Omega nxp
    ! with p oversampling. When applying the filter we should also restart while a condition is
    ! satisfied. This condition can be: if the Ritz value is outside the slice for some
    ! offset then we consider the slice is full and we can stop. One possibility is 
    ! r = F(A)q - q. In the idea that F(A)q ~ q is q is already is the subspace.
    ! the thing is that I don't store q because I apply F(A) in place and overwrite q.
    ! Idea is to measure energy change which is ||F(A)Q||^2. If the change between two iterations is
    ! small then it means that applying the filter does not capture more eigenvectors therefore
    ! it is enough.

    ! ===================== Allocate and fill extended memory buffer ================================== 
  
    ! Sanity check
    if ((.not. slice%use_linalg) .or. slice%use_colsrows) then
        ABI_ERROR("not in linalg representation")
    end if

    slice%neigenpairs_ext = sum(slice%neigenpairs_per_slice)

    if (slice%paral_kgb==0) then
        slice%XextLinalg = X0
    else

        write(std_out,*) 'allocating extended space of size', slice%neigenpairs_ext
        flush(std_out)

        ! Allocate extended space in linalg representation
        call xg_init(slice%X_ext, slice%space, slice%spacedim, slice%neigenpairs_ext, &
            slice%spacecom, me_g0=slice%me_g0, gpu_option=slice%gpu_option)
        
        slice%XextLinalg = slice%X_ext%self

        ! Copy X to XextLinalg by column blocks
        do islice=1,slice%nslice
            ncols = slice%neigenpairs_per_slice(islice)
            fcol = slice%fcol_in_X(islice)
            fcol_ext = slice%fcol_in_Xext(islice)
            call xgBlock_setBlock(X0, slicecols_in, slice%spacedim, ncols, fcol=fcol)
            call xgBlock_setBlock(slice%XextLinalg, slicecols_ext_out, slice%spacedim, ncols, fcol=fcol_ext)
            call xgBlock_copy(slicecols_in, slicecols_ext_out)
            ! Remember xgBlock_copy is always on CPU expect if both blocks are on GPU
        end do
    end if

    ! Unitary test
    if (slice%paral_kgb==1 .and. cols(slice%XextLinalg) /= slice%neigenpairs_ext) then
        ABI_ERROR('wrong linalg representation')
    end if
    write(std_out,'(a,i6,i6)') '# proc has # cols of Xext_linalg ', xmpi_comm_rank(slice%spacecom), cols(slice%XextLinalg)

    ABI_NVTX_END_RANGE()
    call timab(tim_slice_sched,2,tsec)

end subroutine slice_allschedule
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_run
!! NAME
!! slice_run
!!
!! FUNCTION
!! Run Spectrum slicing on a given set of active vectors. 
!! Computation uses marked resources for current process **only**.
!! 
!! INPUTS
!! getAX_BX= pointer to the function giving A|X> and B|X>
!!           A is typically the Hamiltonian H, and B the overlap operator S
!! getBm1X= pointer to the function giving B^-1|X>
!!          B is typically the overlap operator S
!! nspinor= number of spinorial components of the wavefunctions
!!
!! OUTPUT
!! eigen=converged eigenvalue array of size (neigenpairs,1), first rows written only
!! residu=slice residual array of size (neigenpairs,1), first rows written only
!! slice%XextLinalg= guess/converged eigenvectors for all slices
!! 
!! SIDE EFFECTS
!! slice <type(slice_t)>= memory workspace used for Spectrum slicing
!! 
!! SOURCE

subroutine slice_run(slice, getAX_BX, getBm1X, eigen, residu, nspinor)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: residu
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

    ! Timers
    integer, parameter :: tim_slice_1 = 2162
    integer, parameter :: tim_slice_2 = 2163
    integer, parameter :: tim_slice_3 = 2164
    integer, parameter :: tim_slice_X = 2165
    integer :: tim_slice_me
    ! Variables
    type(chebfi_t) :: chebfi
    type(xgBlock_t) :: X0_active
    type(xgBlock_t) :: eigen_active
    type(xgBlock_t) :: residu_active
    integer :: nbdbuf, oracle, num_proc
    integer :: i, iband, nrows
    integer :: me_nbdbuf
    integer :: neigenpairs, bandpp, ndeg_filter
    integer :: comm, comm_rows, comm_cols
    integer :: num_restart
    integer :: num_kept
    integer :: ierr
    logical :: has_converged
    real(dp) :: lambda_minus, lambda_plus
    real(dp) :: theta
    real(dp) :: oracle_factor, oracle_min_occ
    real(dp) :: a_part, b_part, max_resid_kept
    logical :: is_lowpass, on_host, on_device
    ! Arrays
    real(dp) :: tsec(2)
    integer, allocatable, target :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null() 
    integer, pointer :: ncolsColsRows_ptr(:) => null()
    real(dp), pointer :: thetas_conv(:,:) => null()
    real(dp), pointer :: residu_conv(:,:) => null()
    
    ! *********************************************************************
    
    if (slice%me_id_slice==1) then
        tim_slice_me = tim_slice_1
    else if (slice%me_id_slice==2) then
        tim_slice_me = tim_slice_2
    else if (slice%me_id_slice==3) then
        tim_slice_me = tim_slice_3
    else
        tim_slice_me = tim_slice_X
    end if
    
    call timab(tim_slice_me,1,tsec)

    nrows = slice%spacedim
    if (slice%paral_kgb==1) then
        nrows = slice%total_spacedim
    end if
    
    ! Sanity check
    if ( (.not. slice%use_linalg) .and. slice%use_colsrows) then
        ABI_ERROR("should be in linalg representation")
    end if
    !if (slice%gpu_option==ABI_GPU_OPENMP) then
    !    call slice_queryHostDevice(slice, on_host, on_device)
    !    ABI_CHECK(on_device,"GPU not used when it should be!")
    !end if
   
    ! ========================== Transpose ===================================
    !! Function
    ! This transposition allows for a slice to not see others. 
    ! It serves as a transition from global communicator to slice communicator.
    ! 
    !! Notes
    ! Distributes extended columns across **all** MPI processes
    ! After the transposition each process contains the correct
    ! bandpp corresponding to the slice so that no additional communication
    ! has to be performed in order to bring band slices to processes.
    ! Input is Xext (linalg state) distributed across global processes.
    ! Attention chebfi%X is distributed across slice processes =/= Xext per process.
    ! Extracting chebfi%X in the slice distribution from Xext would require comms.
    ! 
    ncolsColsRows_ptr => slice%ncolsColsRows

    ! Allocate slice%me_Xext_active according to the target MPI distribution for slices
    call xgTransposer_constructor(slice%xgTransposerXext, slice%XextLinalg, slice%me_Xext_active,&
        nspinor, STATE_LINALG, TRANS_ALL2ALL, slice%comm_rows, slice%comm_cols, 0, 0, slice%me_g0_fft,&
        gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit,&
        custom_ncolsColsRows=.true., ncolsColsRows_sub=ncolsColsRows_ptr)
   
    slice%xgTransposerXext%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
    
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(slice%xgTransposerXext, STATE_COLSROWS)
    ABI_NVTX_END_RANGE()

    slice%use_colsrows = .true.
    slice%use_linalg = .false.

    ! Unitary test
    if ( slice%paral_kgb == 1 .and.&
        cols(slice%me_Xext_active) /= slice%ncolsColsRows(xmpi_comm_rank(slice%spacecom)+1) ) then
        ABI_ERROR('wrong colsrows representation')
    end if
    write(std_out,'(a,i6,i6,i6)') '# proc has # cols of Xext ', xmpi_comm_rank(slice%spacecom),&
        cols(slice%me_Xext_active)

    ! Get parameters of active task
    neigenpairs = slice%me_neigenpairs_slice
    bandpp = slice%me_bandpp_slice
    ndeg_filter = slice%me_ndeg_slice
    comm = slice%me_comm_slice
    comm_rows = slice%me_comm_rows
    comm_cols = slice%me_comm_cols
    if (slice%me_id_slice==1) then   
        is_lowpass = .true. ! [lambda_minus,lambda_plus) to be diminished
        lambda_minus = slice%poly_upp_bounds(slice%me_id_slice)
        lambda_plus = slice%maxeig_global
    else
        is_lowpass = .false. ! [lambda_minus,lambda_plus) to be amplified
        lambda_minus = slice%poly_low_bounds(slice%me_id_slice)
        lambda_plus = slice%poly_upp_bounds(slice%me_id_slice)
    end if
    ! deactivate chebfi oracle
    oracle = 0
    nbdbuf = 0
    oracle_factor = 1.d0
    oracle_min_occ = 0.d0

    num_proc = xmpi_comm_size(comm)
    ABI_MALLOC_IFNOT(nrowsLinalg,(num_proc))
    nrowsLinalg_ptr => nrowsLinalg
    nrowsLinalg = slice%me_nrowsLinalg_slice

    write(std_out,*) "Allocating slice space..", slice%total_spacedim, neigenpairs
    write(std_out,*) "bands per process=", bandpp
    flush(std_out)

    ! Initialize chebfi object in MPI Colsrows distribution
    call chebfi_init(chebfi,neigenpairs,slice%total_spacedim,slice%tolerance,slice%ecut,slice%paral_kgb,bandpp,&
        ndeg_filter,nbdbuf,slice%space,1,comm,slice%me_g0,slice%me_g0_fft,slice%paw,comm_rows,comm_cols,&
        oracle,oracle_factor,oracle_min_occ,slice%gpu_option,gpu_kokkos_nthrd=slice%gpu_kokkos_nthrd,&
        gpu_thread_limit=slice%gpu_thread_limit,from_linalg=.false.)
 
    ! Define pointers to actively used arrays
    !call xmpi_barrier(slice%spacecom)
    X0_active = slice%me_Xext_active
    call xgBlock_reshape(eigen, 1, slice%neigenpairs)
    call xgBlock_setBlock(eigen, eigen_active, rows=1, cols=neigenpairs, fcol=slice%fcol_in_X(slice%me_id_slice))
    call xgBlock_reshape(eigen_active, neigenpairs, 1)
    call xgBlock_reshape(eigen, slice%neigenpairs, 1)
    call xgBlock_setBlock(residu, residu_active, rows=neigenpairs, cols=1)
   
    write(std_out,*) 'calling runSlice from rank and subrank', xmpi_comm_rank(slice%spacecom), xmpi_comm_rank(comm)
 
    a_part = slice%part_low_bounds(slice%me_id_slice)
    b_part = slice%part_upp_bounds(slice%me_id_slice)
    if (is_lowpass) then
        write(std_out,*) 'spectral offset left=', lambda_minus - b_part
    else
        write(std_out,*) 'spectral offset left=', a_part - lambda_minus
        write(std_out,*) 'spectral offset right=', lambda_plus - b_part
    end if
    
    !write(std_out,*) 'eigen_active='
    !call xgBlock_print(eigen_active, std_out)
   
    ! TODO 
    ! 2) make for 1 MPI

    !num_restart = 200
    num_restart = 1

    i = 0
    max_resid_kept = 1e10
    do while ( (max_resid_kept > slice%ramp) .and. (i < num_restart) )
       
        i = i + 1
        chebfi%xXColsRows = X0_active
        
        write(std_out,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(std_out,*) 'Restart', i, 'with dimensions' 
        write(std_out,*) 'rows=', rows(X0_active), rows(chebfi%xXColsRows), rows(chebfi%X_next)
        write(std_out,*) '     ', rows(chebfi%xAXColsRows), rows(chebfi%xBXColsRows)
        write(std_out,*) 'cols=', cols(X0_active), cols(chebfi%xXColsRows), cols(chebfi%X_next)
        write(std_out,*) '     ', cols(chebfi%xAXColsRows), cols(chebfi%xBXColsRows)
        write(std_out,*) 'target tol=', slice%ramp
        write(std_out,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        flush(std_out)
    
        !call chebfi_runSlice(chebfi, X0_active, getAX_BX, getBm1X, eigen_active, residu_active, nspinor,&
        !    slice%mineig_global, slice%maxeig_global, lambda_minus, lambda_plus, is_lowpass, slice%neigenpairs,&
        !    nrowsLinalg_ptr)

        call chebfi_runSubspaceIteration(chebfi, X0_active, getAX_BX, getBm1X, eigen_active, residu_active, &
            nspinor, slice%mineig_global, slice%maxeig_global, lambda_minus, lambda_plus, is_lowpass, &
            slice%neigenpairs, nrowsLinalg_ptr)


        ! compute residuals of eigenvalues in slice
        ! reverseMap for active eigen and residuals
        write(std_out,*) 'eigen_active dims=', rows(eigen_active), cols(eigen_active)
        write(std_out,*) 'reisd_active dim=', rows(residu_active), cols(residu_active)
        write(std_out,*) 'eigen_active space=', space(eigen_active)
        write(std_out,*) 'reisd_active space=', space(residu_active)
        write(std_out,*) 'slice%neigenpairs=', slice%neigenpairs
        write(std_out,*) 'neigenpairs=', neigenpairs
        flush(std_out)

        if (slice%gpu_option == ABI_GPU_OPENMP) then
            call xgBlock_copy_to_gpu(eigen_active)
            call xgBlock_copy_to_gpu(residu_active)
        end if
       
        !me_nbdbuf = neigenpairs + 100
        !if (slice%me_id_slice==slice%nslice) then
        !    if (nbdbuf > bandpp) then
        !        ABI_WARNING("nbdbuf too large")
        !    end if
        !    if (xmpi_comm_rank(comm)==xmpi_comm_size(comm)-1) then
        !        me_nbdbuf = neigenpairs - slice%nbdbuf
        !        write(std_out,*) 'me_nbdbuf=', me_nbdbuf, " neigenpairs=", neigenpairs
        !    end if
        !end if
        !flush(std_out)

        call xgBlock_reverseMap(eigen_active , thetas_conv, rows=neigenpairs, cols=1)
        call xgBlock_reverseMap(residu_active, residu_conv, rows=neigenpairs, cols=1)
        max_resid_kept = -1e10 ! reset
        num_kept = 0
        do iband=1, neigenpairs
            !if (iband>me_nbdbuf) then
            !    write(std_out,*) 'excluding ', iband, ' in nbdbuf ', me_nbdbuf
            !    flush(std_out)
            !else
                theta = thetas_conv(iband, 1)
                has_converged = .false.
                if (is_lowpass) then
                    has_converged = ( theta < b_part )
                else 
                    has_converged = ( (a_part < theta) .and. (theta < b_part) )
                end if
                if (has_converged) then
                    max_resid_kept = max(max_resid_kept, residu_conv(iband, 1))
                    num_kept = num_kept + 1
                end if
            !end if
        end do
        !write(std_out,*) 'resid debug=', max_resid_kept
        !write(std_out,*) residu_conv(:,1)
        flush(std_out)
        call xmpi_sum(num_kept, comm, ierr)
        call xmpi_max(max_resid_kept, comm, ierr) ! entire slice
        write(std_out,*) '################################################# '
        write(std_out,'(a,i5)') ' Convergence of inner iteration=', i
        write(std_out,*) 'partition             =', a_part, b_part
        write(std_out,*) 'eigenspace dimension  =', num_kept
        write(std_out,*) 'max resid(excl nbdbuf)=', max_resid_kept
        write(std_out,*) '################################################# '
        flush(std_out)

        write(std_out,*) 'residuals='
        do iband=1, neigenpairs
            write(std_out,*) residu_conv(iband, 1)
            flush(std_out)
        end do

        ! todo diagnostic
        ! count how may eigenvalues converged in slice and outside slice but in overlap
        ! compare with expected count

        ! Prepare next iteration
        ! reinitialize pointers to workspaces ... otherwise invovl complains
        call xg_setBlock(chebfi%X_NP, chebfi%X_next, nrows, bandpp)
        call xg_setBlock(chebfi%X_NP, chebfi%X_prev, nrows, bandpp, fcol=bandpp+1)

    end do

    write(std_out,*) 'converged at ninner=', i
    flush(std_out)

    !write(std_out,*) 'getid after runSlice X0_active', xgBlock_getId(X0_active) 

    write(std_out,*) 'chebfi%eigenvalues converged='
    call xgBlock_print(chebfi%eigenvalues,std_out)

    !write(std_out,*) 'residuals='
    !call xgBlock_print(residu_active,std_out)
    !flush(std_out)

    ! Free temporary memory
    call chebfi_free(chebfi)
    ABI_SFREE(nrowsLinalg)

    ! Sanity check
    if ( (.not. slice%use_colsrows) .and. slice%use_linalg) then
        ABI_ERROR("should be in colsrows representation")
    end if
    !if (slice%gpu_option==ABI_GPU_OPENMP) then
    !    call slice_queryHostDevice(slice, on_host, on_device)
    !    ABI_CHECK(on_device,"GPU not used when it should be!")
    !end if

    ! Timer is BEFORE the barrier !!
    call timab(tim_slice_me,2,tsec)

    ! Actually do the transposition to linalg
    call xmpi_barrier(slice%spacecom)
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(slice%xgTransposerXext, STATE_LINALG)
    ABI_NVTX_END_RANGE()
    ! Note: At this point slice%me_Xext_active is recovered into slice%XextLinalg

    slice%use_colsrows = .false.
    slice%use_linalg = .true.

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

subroutine slice_prepareSpectrum(slice, X, lowb, uppb, c_split, bands_left, bands_right, &
        getAX_BX, getBm1X, nspinor)

    implicit none

    ! Arguments
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X
    real(dp), intent(out) :: lowb, uppb, c_split
    integer, intent(out) :: bands_left(:)
    integer, intent(out) :: bands_right(:)
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
    integer :: space, spacecom, gpu_option
    integer :: k_sketch
    real(dp) :: lambda_min, res_norm
    real(dp) :: center, radius
    real(dp) :: mineig, maxeig
    real(dp) :: mineig_global, maxeig_global
    real(dp) :: lanczos_lowb, lanczos_lowb_global
    ! Derived types
#ifdef HAVE_OPENMP_OFFLOAD
    integer :: me_g0
    type(xg_t) :: W_dummy
    integer :: work_size
#endif
    type(xg_t) :: BX
    type(xg_t) :: X_sketch
    type(xgBlock_t) :: xXColsRows
    type(xgTransposer_t) :: xgTransposerX
    type(mpi_type) :: mpi_enreg_old
    type(mpi_type) :: mpi_enreg_new
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
    k_sketch = neigenpairs
    call xg_init(X_sketch, space, spacedim, neigenpairs, spacecom, gpu_option=gpu_option)
    call randomSketching(slice, X, X_sketch%self, k_sketch)
    call xgBlock_copy(X_sketch%self, X)
    call xg_free(X_sketch)

    ! ============== Transpose ==============
    if (slice%paral_kgb==1) then

        ! Allocate memory for X in colsrows representation
        call xgTransposer_constructor(xgTransposerX, X, xXColsRows, nspinor, STATE_LINALG,&
            TRANS_ALL2ALL, slice%comm_rows, slice%comm_cols, 0, 0, slice%me_g0_fft,&
            gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit)
         
        xgTransposerX%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
        
        ! Sanity check
        if ((.not. slice%use_linalg) .or. slice%use_colsrows) then
            ABI_ERROR("not in linalg")
        end if
    
        !call xmpi_barrier(slice%spacecom)
        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
        call xgTransposer_transpose(xgTransposerX, STATE_COLSROWS)
        ABI_NVTX_END_RANGE()

        slice%use_linalg = .false.
        slice%use_colsrows = .true.

    else

        ! Use colsrows notion instead of X notion
        call xgBlock_setBlock(X, xXColsRows, spacedim, neigenpairs)

    end if

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

    kmax = 30
    call computeBLanczos(slice, getAX_BX, getBm1X, kmax, lambda_min, res_norm)

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

    call computeTraceEstimation(slice, getAX_BX, getBm1X, ndeg_filter_max, m_probe,& 
        lanczos_lowb_global)
    
    write(std_out,*) 'STE exited'
    flush(std_out)

    ! TODO 
    ! it would be nice to support nslice=2 and nslice=3
    ! how to split: sketching the restart technique.
    ! Step 1. find mass flip c_1; a < c_1 < b.
    ! Step 2. Cut to this mass flip and restrict to [c_1,b)
    ! Step 3. Repeat step 1 to find mass flip c_2; c_1 < c_2 < b.
    !> other technique
    ! capable of detecting gaps
    ! detects steps where the mass stays constant. This is the
    ! criterion of the constant mass.
    

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

    slice%use_linalg = .true.
    slice%use_colsrows = .false.

    if (slice%paral_kgb == 1) then
        call xgTransposer_free(xgTransposerX)
    end if

end subroutine slice_prepareSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_cutSpectrum
!! NAME
!! slice_cutSpectrum
!! 
!! FUNCTION
!! Split spectrum theta \in [lambda_minus, lambda_plus) into overlapping slices. 
!! Store parameters into arrays of size nslice into the datatype 'slice'.
!! 
!! INPUTS
!! lambda_minus= lower bound of interval to split
!! lambda_plus= upper bound of interval to split
!! theta= array of **sorted** Rayleigh quotients of size (neigenpairs)
!! plot_filter= (option) true if print x,f(x) 
!!
!! SIDE EFFECTS
!! slice%neigenpairs_per_slice
!! slice%fcol_in_X
!! slice%fcol_in_Xext
!! slice%poly_degrees
!! slice%part_low_bounds
!! slice%part_upp_bounds
!! slice%poly_low_bounds
!! slice%poly_upp_bounds
!! 
!! SOURCE

subroutine slice_cutSpectrum(slice, lambda_minus, lambda_plus, theta, plot_filter)

    implicit none

    !Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    real(dp), intent(in) :: lambda_minus
    real(dp), intent(in) :: lambda_plus
    real(dp), pointer, intent(in) :: theta(:)
    logical, optional, intent(in) :: plot_filter
    
    !Local variables-------------------------------
    integer :: iband, islice, jmax, nvec, ndeg, nslice
    integer :: neigenpairs, n_frac
    integer :: first_col, first_col_ext, last_col, last_col_ext
    integer :: ndeg_max = 200
    real(dp) :: ramp, width, wovlp, center, radius
    real(dp) :: poly_low, poly_upp, part_low, part_upp
    real(dp) :: l, u,lw,uw,f_l,f_lw,f_u,f_uw
    real(dp) :: tol12 = 1.0e-12
    logical :: plot_filter_
    ! arrays
    integer, allocatable :: jperm(:)
    real(dp), allocatable :: consdiff(:)
    real(dp), allocatable :: spectral_partition(:)

    ! *********************************************************************

    ramp = slice%ramp
    nslice = slice%nslice
    neigenpairs = slice%neigenpairs
    plot_filter_ = .false.
    if (present(plot_filter)) plot_filter_ = plot_filter

    ABI_MALLOC_IFNOT(spectral_partition,(nslice+1)) 
   
    ! Compute spectral partition using spectral cut strategy
    spectral_partition = 0.d0
    spectral_partition(1) = lambda_minus
    spectral_partition(nslice+1) = lambda_plus
    select case(slice%spectral_cut)
    case(DIVIDE_INTERVAL_WIDTH)

        width = (lambda_plus - lambda_minus) / nslice
        spectral_partition(2:nslice) = (/ (lambda_minus + width*islice, islice=1,nslice-1) /)

    case(DIVIDE_NUMBER_OF_VECTORS)

        n_frac = neigenpairs / nslice
        spectral_partition(2:nslice) = (/ (theta(n_frac * islice), islice=1,nslice-1) /)

    case(DIVIDE_SPECTRAL_GAPS)

        ABI_MALLOC_IFNOT(jperm, (neigenpairs-1))
        ABI_MALLOC_IFNOT(consdiff, (neigenpairs-1))
        
        jperm = (/ (iband, iband=1,neigenpairs-1) /)
        ! weighted by position in spectrum: extremal has large weight
        consdiff = (/ ((theta(iband + 1) - theta(iband))*&
            1.d0/max(theta(iband)-lambda_minus,lambda_plus-theta(iband+1)), iband=1,neigenpairs-1) /)
        
        ! Sort consecutive differences (=gaps) by increasing order
        call sort_dp(neigenpairs-1, consdiff, jperm, tol12)

        ! Take median of largest gaps
        do islice=1,nslice-1
            jmax = jperm(neigenpairs - islice)
            spectral_partition(islice + 1) = (theta(jmax) + theta(jmax + 1)) / 2.d0
        end do

        ABI_SFREE(jperm)
        ABI_SFREE(consdiff)

    end select

    ! Center and radius of the entire spectrum mapped to -1,1
    center = (slice%mineig_global + slice%maxeig_global) / 2.d0
    radius = (slice%maxeig_global - slice%mineig_global) / 2.d0

    ! Define spectral subintervals and optimize degrees for individual slices
    first_col_ext = 1
    do islice=1, nslice 

        part_low = spectral_partition(islice)
        part_upp = spectral_partition(islice+1)
        wovlp = (part_upp - part_low)/8.d0 ! FIXME add abi parameter to tune this
        poly_low = part_low - wovlp
        poly_upp = part_upp + wovlp

        if (islice==1) then
            ! Amplification for first slice is f(u+w)/f(u)
            ndeg = slice%ndeg_filter - 1
            f_uw = 1.d0; f_u = 1.d0
            do while( f_uw/f_u > ramp .and. ndeg < ndeg_max )
                ndeg = ndeg + 1
                f_u = cheb_poly(part_upp,ndeg,poly_upp,slice%maxeig_global)
                f_uw = cheb_poly(poly_upp,ndeg,poly_upp,slice%maxeig_global)
            end do
            write(std_out,*) 'slice 1 amplif factor f(out)/f(in)=', f_uw / f_u 
        else
            ! Amplification is f(l-w)/f(l) (left) and f(u+w)/f(u) (upper)
            lw = (poly_low - center)/radius ! scaled point outside slice
            uw = (poly_upp - center)/radius ! scaled point outside slice
            l = (part_low - center)/radius ! scaled point inside slice
            u = (part_upp - center)/radius ! scaled point inside slice
            ndeg = 8
            f_lw = 1.d0; f_uw = 1.d0; f_l = 1.d0; f_u = 1.d0
            do while ( (f_lw/f_l > ramp .or. f_uw/f_u > ramp) .and. ndeg < ndeg_max )
                ndeg = ndeg + 1
                f_l  = bandpassIndicator_sca(l ,lw,uw,ndeg)
                f_lw = bandpassIndicator_sca(lw,lw,uw,ndeg)
                f_u  = bandpassIndicator_sca(u ,lw,uw,ndeg)
                f_uw = bandpassIndicator_sca(uw,lw,uw,ndeg)
            end do
            write(std_out,*) 'left/right amplif factor f(out)/f(in)=', f_lw/f_l, f_uw/f_u
        end if

        ! Plot filter in interval [glb, ub) (set manually because depends on the case)
        if (plot_filter_) then
            call print_scalar_filter(slice%mineig_global,poly_upp,poly_low,&
                poly_upp,slice%mineig_global,slice%maxeig_global,ndeg,(islice==1))
        end if

        ! Count theta eigenvalues in current spectral partition with overlap 
        if (nslice==1) then
            first_col = 1
            last_col = neigenpairs
        else
            first_col = maxloc(theta, dim=1, mask=(theta < poly_low)) + 1
            last_col = maxloc(theta, dim=1, mask=(theta < poly_upp))       
            if (islice == 1) first_col = 1
            if (islice == nslice) last_col = neigenpairs
        end if
        nvec = last_col - first_col + 1

        ! Print slice interval info
        write(std_out,'(a,i2,a,i6,a,i6,a)') '======= Slice ', islice, ' nvec=', nvec, ' ndeg=', ndeg,&
&               '      lb, ub, width'
        write(std_out,*) part_low, part_upp, part_upp - part_low
        write(std_out,*) poly_low, poly_upp, poly_upp - poly_low

        ! Compute last index in extended memory (without ovlp)
        last_col_ext = first_col_ext + nvec - 1

        slice%neigenpairs_per_slice(islice) = nvec
        slice%fcol_in_X(islice) = first_col
        slice%fcol_in_Xext(islice) = first_col_ext
        slice%poly_degrees(islice) = ndeg
        slice%part_low_bounds(islice) = part_low
        slice%part_upp_bounds(islice) = part_upp
        slice%poly_low_bounds(islice) = poly_low
        slice%poly_upp_bounds(islice) = poly_upp
 
        ! Update starting index of next slice in extended
        first_col_ext = last_col_ext + 1
           
    end do

    ABI_SFREE(spectral_partition)
   
end subroutine slice_cutSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allocateResources
!! NAME
!! slice_allocateResources
!! 
!! FUNCTION
!! Assign available resources to smaller slice tasks
!! in an efficient way to achieve load balancing.
!! Every process knows the resource allocation.
!!
!! INPUTS
!! slice%nproc                 =number of resources to be splitted
!! slice%nslice                =number of tasks
!! slice%neigenpairs_per_slice =number of eigenpairs per slice
!! slice%poly_degree           =polynomial filter degree per slice
!! 
!! SIDE EFFECTS
!! slice%nproc_per_slice       =number of allocated resources per task
!! slice%lookup_proc           =identifier of assigned task per resource
!! 
!! SOURCE

subroutine slice_allocateResources(slice)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice

    ! Local variables
    integer :: iproc
    ! Arrays
    integer, allocatable :: weights(:)
    
    ! *********************************************************************

    if (slice%paral_slice==EVEN_SLICES) then
       
        slice%nproc_per_slice = ceiling(real(slice%nproc) / real(slice%nslice))

        ! remove from first slice if sum exceeds total nproc, allocation must sum to nproc
        if (modulo(slice%nproc,slice%nslice)/=0) then
            slice%nproc_per_slice(1) = 0
            slice%nproc_per_slice(1) = slice%nproc - sum(slice%nproc_per_slice)
        end if
        
        call assign_tasks_to_processes(slice%nproc_per_slice, slice%lookup_proc)
    else

        ABI_MALLOC_IFNOT(weights, (slice%nslice))

        ! Apply weighted fair allocation with a load balance criterion
        select case(slice%paral_slice)
        case(FAIR_BANDPP)
            weights = 1 
        case(FAIR_BANDPP_WDEG)
            weights = slice%poly_degrees
        end select
    
        ! Solve allocation problem to find the amount of resource allocated to each slice
        call fair_allocation(slice%nslice, slice%neigenpairs_per_slice, weights, slice%nproc, &
                slice%nproc_per_slice)
    
        ! Call the subroutine to assign tasks (=slices) to processes
        call assign_tasks_to_processes(slice%nproc_per_slice, slice%lookup_proc)

        ! Free temporary memory
        ABI_SFREE(weights) 

    end if
   
    call xmpi_barrier(slice%spacecom)
   
    do iproc = 1, slice%nproc
        write(std_out,'(a,i5,a,i5)') "Process ", iproc-1, " allocated to task ", slice%lookup_proc(iproc)
    end do

end subroutine slice_allocateResources
!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_markActiveTask
!! NAME
!! slice_markActiveTask
!!
!! FUNCTION
!! Mark allocated portion as actively used. Also mark all me_ variables.
!! My process only marks resources its assigned slice has reserved. 
!!
!! SIDE EFFECTS
!! slice%me_comm_slice= sub-communicator of active slice
!! slice%me_comm_rows= used in Transposer
!! slice%me_comm_cols= used in Transposer
!! slice%me_ncolsColsRows_slice= column capacity per process in colsrows repr
!! slice%me_nrowsLinalg_slice= row capacity per process in linalg repr
!! slice%me_id_slice= slice index starting from 1,..,nslice
!! slice%me_nproc_slice= number of processes reserved for slice
!! slice%me_ndeg_slice= polynomial degree of slice filter
!! 
!! SOURCE

subroutine slice_markActiveTask(slice)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice

    ! Local variables
    integer :: my_rank, my_rank_sub, ierr

    ! *********************************************************************
 
    my_rank = xmpi_comm_rank(slice%spacecom) 
    slice%me_id_slice = slice%lookup_proc(my_rank + 1) + 1
    slice%me_neigenpairs_slice = slice%neigenpairs_per_slice(slice%me_id_slice)
    slice%me_nproc_slice = slice%nproc_per_slice(slice%me_id_slice)
    slice%me_ndeg_slice = slice%poly_degrees(slice%me_id_slice)

    ABI_MALLOC_IFNOT(slice%me_ncolsColsRows_slice, (slice%me_nproc_slice))
    ABI_MALLOC_IFNOT(slice%me_nrowsLinalg_slice, (slice%me_nproc_slice))

    if (slice%paral_kgb==1) then
        ! Compute column distribution across active resources
        call distribute_vectors(slice%me_neigenpairs_slice, slice%me_nproc_slice, slice%me_ncolsColsRows_slice)
    
        ! Compute row distribution across active resources
        call distribute_vectors(slice%total_spacedim, slice%me_nproc_slice, slice%me_nrowsLinalg_slice)
    end if

    ! Split global comm into disjoint sub-comms, only procs with the same color communicate
    if (slice%paral_kgb==0) then
        slice%me_comm_slice = slice%spacecom
        slice%me_comm_rows = slice%comm_rows
        slice%me_comm_cols = slice%comm_cols
    else
        call xmpi_comm_split(slice%spacecom, slice%me_id_slice, my_rank, slice%me_comm_slice, ierr)        
        if ( ierr /= xmpi_success ) then
            ABI_ERROR("Error while creating slice spacecom subcommunicator")
        end if
        call xmpi_comm_split(slice%comm_rows, slice%me_id_slice, my_rank, slice%me_comm_rows, ierr)          
        if ( ierr /= xmpi_success ) then
            ABI_ERROR("Error while creating slice row subcommunicator")
        end if      
        call xmpi_comm_split(slice%comm_cols, slice%me_id_slice, my_rank, slice%me_comm_cols, ierr)         
        if ( ierr /= xmpi_success ) then
            ABI_ERROR("Error while creating slice col subcommunicator")
        end if 
    end if
    
    ! process waits for others to create their subcommunicators before using its own
    call xmpi_barrier(slice%spacecom) 

    ! Concatenate slice%me_ncolsColsRows_slice into collective slice%ncolsColsRows
    my_rank_sub = xmpi_comm_rank(slice%me_comm_slice)
    slice%me_bandpp_slice = slice%bandpp
    if (slice%paral_kgb==1) then
        slice%me_bandpp_slice = slice%me_ncolsColsRows_slice(my_rank_sub + 1)
    end if
    call xmpi_allgather(slice%me_bandpp_slice, slice%ncolsColsRows, slice%spacecom, ierr)
    if ( ierr /= xmpi_success ) then
        ABI_ERROR("Error while gathering number of columns in colsrows for all slices")
    end if

end subroutine slice_markActiveTask
!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allmerge
!! NAME
!! slice_allmerge
!! 
!! FUNCTION
!! Filter converged eigenvalues of each slice using a criterion
!! based on spectral partition interval bounds.
!!
!! INPUTS
!! eigen= slice eigenvalues written in first 1,..,neigenpairs_slice rows of (neigenpairs,1) 
!! resid= slice residuals written in first 1,..,neigenpairs_slice rows of (neigenpairs,1)
!! slice%me_Xext_active= converged slice eigenvectors of size (total_spacedim,bandpp)
!! 
!! OUTPUT
!! X0= converged eigenvectors in Linalg representation of size (spacedim, neigenpairs)
!! eigen= all converged eigenvalues of size (neigenpairs, 1)
!! resid= all residuals of size (neigenpairs,1)
!! 
!! SIDE EFFECTS
!! slice <type(slice_t)>= memory workspace used for Spectrum slicing
!! slice%XextLinalg= converged eigenvectors in extended space of size (spacedim,neigenpairs_ext)
!! slice%neigenpairs_per_slice= number of kept eigenpairs after merging
!! slice%fcol_in_Xext= first index to copy from extended memory in 1,..,neigenpairs_ext
!! slice%fcol_in_X= first index to copy to regular memory in 1,..,neigenpairs
!!
!! SOURCE

subroutine slice_allmerge(slice, X0, eigen, resid)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X0
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: resid

    ! Local variables-------------------------------
    integer :: my_rank, my_slice, neigenpairs_slice, nkept
    integer :: fcol_ext, fcol, tot_ncols_kept, lcol_ext
    integer :: islice, fcol_in_slice, lcol_in_slice, rem
    real(dp) :: part_low_bound, part_upp_bound
    logical :: on_host, on_device
    ! Derived types
    type(xg_t) :: eigen_ext
    type(xg_t) :: resid_ext
    type(xgBlock_t) :: eigen_conv
    type(xgBlock_t) :: resid_conv
    type(xgBlock_t) :: eigen_ext_slice
    type(xgBlock_t) :: resid_ext_slice
    type(xgBlock_t) :: X_kept, eigen_kept, resid_kept
    type(xgBlock_t) :: X0_out, eigen_out, resid_out
    ! arrays
    real(dp), pointer :: theta_ext(:,:) => null()
    real(dp), allocatable :: theta_reshaped(:)
 
    ! *********************************************************************

    ! Two ways to get slice eigenvalues to filter
    ! 1) (implemented)
    !    each slice writes its converged eigenvalues into a range of eigen_ext
    !    without overlap. Then we recover this range and we filter it.
    !    This used an extra memory space called eigen_ext.
    ! 2) We recover the slice converged eigenvalues directly from eigen.
    !    The convention is that we wrote to the first neigenpairs_slice
    !    rows. This also uses a pointer pointing to the part of eigen
    !    we want to discard. This part will be set to zero. Finally we 
    !    call xgBlock_mpi_sum(eigen, spacecom) in order to sum non-discarded parts.
    !    This is complicated because we also have to shift parts.

    ! Allocate extended space for all slice eigenvalues and residuals
    call xg_init(eigen_ext, SPACE_R, rows=1, cols=slice%neigenpairs_ext, gpu_option=slice%gpu_option)
    call xg_init(resid_ext, SPACE_R, rows=1, cols=slice%neigenpairs_ext, gpu_option=slice%gpu_option)

    write(std_out,*) 'db xgBlock_zero call'
    call xgBlock_zero(eigen_ext%self)
    call xgBlock_zero(resid_ext%self)

    ! Bring to columns to be able to select column range
    call xgBlock_reshape(eigen, 1, slice%neigenpairs)
    call xgBlock_reshape(resid, 1, slice%neigenpairs)

    ! MPI communication to gather slice eigen/resid to eigen_ext/resid_ext
    if (slice%paral_kgb==1) then

        if (xmpi_comm_size(slice%spacecom) > 1) then
        
            ! Copy only once, eg for first process in slice subcomm 
            if (xmpi_comm_rank(slice%me_comm_slice)==0) then
                my_rank = xmpi_comm_rank(slice%spacecom)
                my_slice = slice%lookup_proc(my_rank + 1)
                neigenpairs_slice = slice%neigenpairs_per_slice(my_slice + 1)
                fcol_ext = slice%fcol_in_Xext(my_slice + 1)
                fcol = slice%fcol_in_X(my_slice + 1)
                ! set blocks copy from 
                call xgBlock_setBlock(eigen, eigen_conv, rows=1, cols=neigenpairs_slice, fcol=fcol)
                call xgBlock_setBlock(resid, resid_conv, rows=1, cols=neigenpairs_slice, fcol=fcol)
                ! set blocks copy to
                call xgBlock_setBlock(eigen_ext%self, eigen_ext_slice, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
                call xgBlock_setBlock(resid_ext%self, resid_ext_slice, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
                ! perform copy
                call xgBlock_copy(eigen_conv, eigen_ext_slice)
                call xgBlock_copy(resid_conv, resid_ext_slice)
            end if
    
            ! All processes wait to finish copying before summing 
            call xmpi_barrier(slice%spacecom)

            call xgBlock_mpi_sum(eigen_ext%self, comm=slice%spacecom)
            call xgBlock_mpi_sum(resid_ext%self, comm=slice%spacecom)

        else
            ABI_BUG("Not implemented!") 
            call xgBlock_copy(eigen, eigen_ext%self)
            call xgBlock_copy(resid, resid_ext%self)
        end if
    else
        ABI_BUG("Not implemented!") 
        call xgBlock_copy(eigen, eigen_ext%self)
        call xgBlock_copy(resid, resid_ext%self)
    end if 

    ! Copy is on CPU so update CPU data from GPU
    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call xgBlock_copy_from_gpu(eigen_ext%self)
        call xgBlock_copy_from_gpu(resid_ext%self)
    end if

    !write(std_out,*) 'residuals='
    !call xgBlock_print(resid_ext%self,std_out)

    ! Results could be complex, so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(eigen_ext%self, theta_ext, rows=1, cols=slice%neigenpairs_ext)

    ! Filter eigenvalues in extended space using spectral partition
    tot_ncols_kept = 0
    do islice=1,slice%nslice

        ! Before merge: get slice eigenvalues to filter
        neigenpairs_slice = slice%neigenpairs_per_slice(islice)
        ABI_MALLOC_IFNOT(theta_reshaped, (neigenpairs_slice)) 
        fcol_ext = slice%fcol_in_Xext(islice)
        lcol_ext = fcol_ext + neigenpairs_slice - 1
        if (islice == slice%nslice) then
            lcol_ext = slice%neigenpairs_ext
        end if
        theta_reshaped(1:neigenpairs_slice) = theta_ext(1,fcol_ext:lcol_ext)

        ! Apply filter criterion to find kept first and last column in slice
        part_low_bound = slice%part_low_bounds(islice)
        part_upp_bound = slice%part_upp_bounds(islice)
        fcol_in_slice = maxloc(theta_reshaped, dim=1, mask=(theta_reshaped < part_low_bound)) + 1
        lcol_in_slice = maxloc(theta_reshaped, dim=1, mask=(theta_reshaped < part_upp_bound))            
        if (islice == 1) then
            fcol_in_slice = 1
        else if (islice == slice%nslice) then
            rem = slice%neigenpairs - tot_ncols_kept
            if (rem < 0) then
                ABI_ERROR("Not enough eigenvalues in last slice. Decrease tolfilter or nstep_mixed.")
            else
                lcol_in_slice = min(neigenpairs_slice, fcol_in_slice + rem - 1)
                write(std_out,*) 'rem= lcol_in_slice=', rem, lcol_in_slice 
            end if
        end if

        nkept = lcol_in_slice - fcol_in_slice + 1

!        write(std_out,*) 'Filter in ', part_low_bound, part_upp_bound
!        write(std_out,*) 'kept indices', fcol_in_slice, lcol_in_slice, nkept
!        !write(std_out,*) 'filtered eigenvalues=', theta_reshaped
!        write(std_out,*) 'filtered eigval(first,last)=', theta_reshaped(1), theta_reshaped(neigenpairs_slice)
!        !write(std_out,*) 'kept eigenvalues=', theta_reshaped(fcol_in_slice:lcol_in_slice)
!        write(std_out,*) 'kept eigval(first,last)=', theta_reshaped(fcol_in_slice), theta_reshaped(lcol_in_slice)
!        write(std_out,*) 'tot_ncols_kept(prev)=', tot_ncols_kept 

        ! After merge: Update first columns to copy from Xext to X
        slice%fcol_in_X(islice)= tot_ncols_kept + 1
        slice%fcol_in_Xext(islice) = fcol_ext + fcol_in_slice - 1
        slice%neigenpairs_per_slice(islice) = nkept
        tot_ncols_kept = tot_ncols_kept + slice%neigenpairs_per_slice(islice)
        
        write(std_out,*) 'tot_ncols_kept(next)=', tot_ncols_kept 

        ABI_SFREE(theta_reshaped)

    end do
    
    ! Detect missing or extra eigenvalues
    if (tot_ncols_kept < slice%neigenpairs) then
        ABI_WARNING("Not enough converged eigenvalues in slice")
    else if (tot_ncols_kept > slice%neigenpairs) then
        ABI_WARNING("Too many converged eigenvalues kept. Decrease tolfilter or nstep_mixed.")
    end if

    ! Copy from extended memory to regular memory
    do islice=1,slice%nslice
        fcol = slice%fcol_in_X(islice)
        fcol_ext = slice%fcol_in_Xext(islice)
        neigenpairs_slice = slice%neigenpairs_per_slice(islice)
        write(std_out,*) 'block copy from fcol, ncols=', fcol_ext, neigenpairs_slice
        write(std_out,*) 'block copy to fcol, ncols=', fcol, neigenpairs_slice
        ! Blocks to copy from
        call xgBlock_setBlock(slice%XextLinalg, X_kept, rows=slice%spacedim, cols=neigenpairs_slice, fcol=fcol_ext)
        call xgBlock_setBlock(eigen_ext%self, eigen_kept, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
        call xgBlock_setBlock(resid_ext%self, resid_kept, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
        ! Blocks to copy to
        call xgBlock_setBlock(X0, X0_out, rows=slice%spacedim, cols=neigenpairs_slice, fcol=fcol)
        call xgBlock_setBlock(eigen, eigen_out, rows=1, cols=neigenpairs_slice, fcol=fcol)
        call xgBlock_setBlock(resid, resid_out, rows=1, cols=neigenpairs_slice, fcol=fcol)
        ! copy
        call xgBlock_copy(X_kept, X0_out)
        call xgBlock_copy(eigen_kept, eigen_out)
        call xgBlock_copy(resid_kept, resid_out)
    end do

    ! Recover dimensions
    call xgBlock_reshape(eigen, slice%neigenpairs, 1) 
    call xgBlock_reshape(resid, slice%neigenpairs, 1)

    !write(std_out,*) 'residuals after merge=', xgBlock_getid(resid)
    !flush(std_out)
 
    !write(std_out,*) 'residuals squared after all merge='
    !call xgBlock_print(resid,std_out)

    !write(std_out,*) 'KEPT slice eigs='
    !call xgBlock_print(eigen, std_out)

    ! Free memory
    call xg_free(eigen_ext)
    call xg_free(resid_ext)

end subroutine slice_allmerge
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_queryHostDevice
!! NAME
!! slice_queryHostDevice
!! 
!! FUNCTION
!! Query OpenMP offloading API to update/get(optional) GPU/CPU flags
!! 
!! SOURCE

subroutine slice_queryHostDevice(slice,on_host,on_device)

    implicit none
    type(slice_t), intent(inout) :: slice
    logical, optional, intent(inout) :: on_host
    logical, optional, intent(inout) :: on_device
    integer :: device_id
    
    ! *********************************************************************

    device_id = xomp_get_device_num()
    slice%on_host = (device_id < 1) ! outside target (-1), inside target host (0)
    slice%on_device = (device_id > 0) ! inside target not host (device number)

    if (present(on_host)) on_host = slice%on_host
    if (present(on_device)) on_device = slice%on_device

end subroutine slice_queryHostDevice
!***

!----------------------------------------------------------------------

!!****f* m_slice/assign_tasks_to_processes
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
!***

!----------------------------------------------------------------------

!!****f* m_slice/distribute_vectors
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

!!****f* m_slice/fair_allocation
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

!!****f* m_slice/adjust_allocation
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

!!****f* m_slice/reduce_allocation
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

!!****f* m_slice/cheb_poly
!! NAME
!! cheb_poly
!!
!! FUNCTION
!! Compute Chebyshev polynomial
!!
!! INPUTS
!!  xx= input variable
!!  aa= left bound of the interval
!!  bb= right bound of the interval
!!  nn=
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

function cheb_poly(xx,nn,aa,bb) result(yy)

    implicit none

    ! Arguments ------------------------------------
    integer,  intent(in) :: nn
    real(dp), intent(in) :: xx, aa, bb
    real(dp)             :: yy

    ! Local variables-------------------------------
    integer  :: ii
    real(dp) :: xred,yim1,temp

    ! *************************************************************************

    xred = (xx-(aa+bb)/2)/(bb-aa)*2
    yy = xred
    yim1 = 1
    do ii= 2, nn
        temp = yy
        yy = 2*xred*yy - yim1
        yim1 = temp
    end do

end function cheb_poly
!!***

!----------------------------------------------------------------------

!!****f* m_slice/bandpass_sca
!! NAME
!! bandpass_sca
!!
!! FUNCTION
!! Computes delta-Dirac polynomial filter f(x) approximated by a Chebyshev
!! expansion of order deg centered at gamma, evaluated at point x=t
!!
!! INPUTS
!!  t=      scalar to evaluate filter on
!!  deg=    order of Chebyshev expansion
!!  gam=    center of Chebyshev expansion
!!
!! OUTPUT
!!  res
!!
!! SOURCE

function bandpass_sca(t, deg, gam) result(f_t)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in ) :: t, gam
    integer , intent(in ) :: deg

    real(dp) :: f_t
    
    !Local variables-------------------------------
    real(dp) :: yt0, yt, yg0, yg, yt_swap, yg_swap
    real(dp) :: mu, damp, rho, rhog, theta
    integer  :: i
    
    ! *********************************************************************

    ! init cheby of deg=0,1 eval at t,gamma
    yt0 = 1.d0
    yt = t
    
    yg0 = 1.d0
    yg = gam

    ! init delta-Dirac filters of deg=0
    theta = Pi/(deg + 1)
    damp = SIN(theta) / theta
    rho = 0.5d0 + gam * damp * yt
    rhog = 0.5d0 + gam * damp * yg

    do i=2,deg 

        ! Update Chebyshev polynomials
        yt_swap = yt
        yt = 2 * t * yt - yt0
        yt0 = yt_swap
        
        yg_swap = yg
        yg = 2 * gam * yg - yg0
        yg0 = yg_swap

        ! Update delta-Dirac filters
        mu = COS(i * ACOS(gam))
        damp = SIN(i * theta) / (i * theta)
        rho = rho + mu * damp * yt
        rhog = rhog + mu * damp * yg
        
    end do

    f_t = rho / rhog

end function bandpass_sca
!!***

!----------------------------------------------------------------------

!!****f* m_slice/print_scalar_filter
!! NAME
!! print_scalar_filter
!! 
!! FUNCTION
!! Print x,f(x) for every x in (a,b).
!! 
!! INPUTS
!! lb,ub= interval to amplify/vanish
!! glb,gub= guaranteed spectral bounds used to scale to [-1,1]
!! ndeg= polynomial degree of filter
!! is_lowpass= flag. If true then Chebyshev if false Chebyshev-Jackson

subroutine print_scalar_filter(a, b, lb, ub, glb, gub, ndeg, is_lowpass)

    implicit none

    real(dp), intent(in) :: a,b,lb,ub,glb,gub
    integer, intent(in) :: ndeg
    logical, intent(in) :: is_lowpass

    integer :: npt, ipt
    real(dp) :: c, r, pt, fun_pt

    npt = 100
    c = (gub + glb)/2.d0
    r = (gub - glb)/2.d0
    write(std_out,*) ' '
    write(std_out,*) 'Plot filter ==== x | f(x)'
    if (is_lowpass) then
        do ipt=1,npt
            pt = a + (ipt-1)*(b-a)/npt ! unscaled!
            fun_pt = cheb_poly(pt,ndeg,ub,gub)
            write(std_out,*) pt, fun_pt
        end do
    else
        do ipt=1,npt
            pt = (a + (ipt-1)*(b-a)/npt - c)/r ! scaled!
            fun_pt = bandpassIndicator_sca(pt,(a-c)/r,(b-c)/r,ndeg)
            write(std_out,*) pt, fun_pt
        end do
    end if
    write(std_out,*) ' '

end subroutine print_scalar_filter
!!***

!----------------------------------------------------------------------

!!****f* m_slice/computeBLanczos
!! NAME
!! computeBLanczos
!! 
!! FUNCTION
!! B-Lanczos three-term recurrence (using B-inner product).
!! Performs k Lanczos iterations on a column vector.
!!
!! SOURCE
  
  subroutine computeBLanczos(slice, getAX_BX, getBm1X, k, lambda_min, res_norm)

    implicit none

    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: k
    real(dp), intent(out) :: lambda_min, res_norm
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
    
    type(xg_t) :: W_vcol
    type(xg_t) :: W_dot
    type(xgBlock_t) :: q, v, Bv, Bm1v, qprev
    type(xgBlock_t) :: dot_qTBv, dot_qTv, dot_vTBv
    real(dp) :: alpha(k), beta(k-1)
    real(dp) :: beta_prev
    real(dp) :: norml_q
    real(dp), pointer :: dot_qTBv_layout(:,:) => null()
    real(dp), pointer :: dot_qTv_layout(:,:) => null()
    real(dp), pointer :: dot_vTBv_layout(:,:) => null()
    real(dp), allocatable :: v_min(:)

    integer :: i, j
    integer :: rank
    integer :: space
    integer :: tot_spacedim
    integer :: me_g0
    integer :: gpu_option
    real(dp) :: tsec(2)

    ! *********************************************************************

    call timab(tim_lanczos,1,tsec)
    
    space = slice%space
    tot_spacedim = slice%total_spacedim
    gpu_option = slice%gpu_option
    me_g0 = slice%me_g0
    if (slice%paral_kgb==1) then
        me_g0 = slice%me_g0_fft
    end if
    rank = xmpi_comm_rank(slice%spacecom)
    beta_prev = 0.0_dp

    ABI_MALLOC(v_min, (tot_spacedim))
    
    write(std_out,*) 'Lanczos in rank=', rank; flush(std_out)

    ! workspace size (npw,5)
    call xg_init(W_vcol, space, tot_spacedim, 5, xmpi_comm_null, me_g0=me_g0, gpu_option=gpu_option)
    call xgBlock_setBlock(W_vcol%self,     q, tot_spacedim, 1)         ! q
    call xgBlock_setBlock(W_vcol%self,     v, tot_spacedim, 1, fcol=2) ! Aq
    call xgBlock_setBlock(W_vcol%self,    Bv, tot_spacedim, 1, fcol=3) ! Bq
    call xgBlock_setBlock(W_vcol%self,  Bm1v, tot_spacedim, 1, fcol=4) ! Bm1 v
    call xgBlock_setBlock(W_vcol%self, qprev, tot_spacedim, 1, fcol=5) ! q_prev

    call xg_init(W_dot, space, 1, 3, xmpi_comm_null, me_g0=me_g0, gpu_option=gpu_option)
    call xgBlock_setBlock(W_dot%self, dot_qTBv, 1, 1)
    call xgBlock_setBlock(W_dot%self,  dot_qTv, 1, 1, fcol=2)
    call xgBlock_setBlock(W_dot%self, dot_vTBv, 1, 1, fcol=3)

    ! q = random column vector
    call xgBlock_colwiseRandom(q, rank, 1)
    !write(std_out,*) 'Random id=', xgBlock_getid(q) 
    !flush(std_out)

    ! Bv = B * q / norml_q
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
    call getAX_BX(q, v, Bv)
    call xgBlock_zero_im_g0(v) ! v stores Aq
    call xgBlock_zero_im_g0(Bv) ! Bv stores Bq
    ABI_NVTX_END_RANGE() 
    call timab(tim_getAX_BX,2,tsec)

    call xgBlock_colwiseDotProduct(q, Bv, dot_qTBv)
    call xgBlock_reverseMap(dot_qTBv,dot_qTBv_layout,rows=1,cols=1)    
    
    norml_q = 1.d0 / sqrt(dot_qTBv_layout(1,1))
    call xgBlock_scale(q, norml_q, 1)
    call xgBlock_scale(v, norml_q, 1)
    call xgBlock_scale(Bv, norml_q, 1)

    alpha = 0.d0
    beta = 0.d0

    do j = 1, k

        if (j>1) then ! rewrite v
            ! v = A * q
            call timab(tim_getAX_BX,1,tsec)
            ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
            call getAX_BX(q, v, Bv)
            call xgBlock_zero_im_g0(v) ! v stores Aq
            call xgBlock_zero_im_g0(Bv) ! Bv stores Bq
            ABI_NVTX_END_RANGE()
            call timab(tim_getAX_BX,2,tsec)
        end if

        ! alpha_j = q^T * Aq
        call xgBlock_colwiseDotProduct(q, v, dot_qTv)
        call xgBlock_reverseMap(dot_qTv,dot_qTv_layout,rows=1,cols=1)
        alpha(j) = dot_qTv_layout(1,1)

        ! v = B^{-1} * A * q
        if (slice%paw) then
            call timab(tim_invovl, 1, tsec)
            ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_BM1X)
            call getBm1X(v, Bm1v)
            ABI_NVTX_END_RANGE()
            call timab(tim_invovl, 2, tsec)
            call timab(tim_copy, 1, tsec)
            call xgBlock_copy(Bm1v, v)
            call timab(tim_copy, 2, tsec)
        end if

        ! v = B^{-1} A q - alpha q - beta_prev q_prev
        call xgBlock_saxpy(v, -1.d0 * alpha(j), q)
        if (j > 1) then
            call xgBlock_saxpy(v, -beta_prev, qprev)
        end if

        ! Compute beta_j if j<k
        if (j < k) then
            
            ! Bv = B * v
            call timab(tim_getAX_BX,1,tsec)
            ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
            call getAX_BX(v, Bm1v, Bv)
            call xgBlock_zero_im_g0(Bm1v) ! Bm1v dummy workspace
            call xgBlock_zero_im_g0(Bv) ! Bv 
            ABI_NVTX_END_RANGE()
            call timab(tim_getAX_BX,2,tsec)
 
            call xgBlock_colwiseDotProduct(v, Bv, dot_vTBv)
            call xgBlock_reverseMap(dot_vTBv,dot_vTBv_layout,rows=1,cols=1)
            beta(j) = sqrt(dot_vTBv_layout(1,1))

            ! Update q_prev, q=v/beta_j, beta_prev
            call xgBlock_copy(q, qprev)
            call xgBlock_scale(v, 1.d0/beta(j), 1)
            call xgBlock_copy(v, q)
            beta_prev = beta(j)
        end if
    end do

    ! Diagonalize T (always on CPU)
    call smallestTridiagEigenpair(k, alpha, beta, lambda_min, v_min)

    ! residual norm using Lanczos shortcut
    res_norm = abs(beta(k-1)*v_min(k))

    call xg_free(W_vcol)
    call xg_free(W_dot)
    ABI_FREE(v_min)
    
    call timab(tim_lanczos,2,tsec)

  end subroutine computeBLanczos
!!***

!----------------------------------------------------------------------

!!****f* m_slice/randomSketching
!! NAME
!! randomSketching
!! 
!! SOURCE
  
  subroutine randomSketching(slice, X, X_sketch, k_sketch)
    
    implicit none

    type(slice_t), intent(in) :: slice
    type(xgBlock_t), intent(in) :: X
    type(xgBlock_t), intent(inout) :: X_sketch
    integer, intent(in) :: k_sketch

    integer :: k
    integer :: rank
    integer :: spacecom, space
    integer :: nband, tot_spacedim
    integer :: gpu_option
    type(xg_t) :: Omega
    type(xgBlock_t) :: q

    ! *********************************************************************

    space = slice%space
    nband = slice%neigenpairs
    tot_spacedim = slice%total_spacedim
    gpu_option = slice%gpu_option
    spacecom = slice%spacecom 

    if (k_sketch > nband) then
        ABI_ERROR("sketching dimension cannot be more than initial one")
    end if

    ! Each MPI has the same sketch matrix
    call xg_init(Omega, space, nband, k_sketch, xmpi_comm_null, gpu_option=gpu_option) 
    
    rank = xmpi_comm_rank(spacecom)

    do k = 1, k_sketch
        ! seed depends on column index
        ! rank * offset + k, with offset > nband to avoid overlap between columns across ranks
        call xgBlock_colwiseRandomGaussian(Omega%self, rank*(k_sketch+10)+k, k)
        
        ! test
        ! q = random column vector
        !call xgBlock_setBlock(Omega%self, q, nband, 1, fcol=k)
        !write(std_out,*) 'Random id=', xgBlock_getid(q) 
        !flush(std_out)
    end do

    ! Compute X * Omega
    call xgBlock_gemm('n','n',1.0d0,X,Omega%self,0.d0,X_sketch,comm=xmpi_comm_null)
    !call xgBlock_copy(Omega%self, X_sketch)

    call xg_free(Omega)

  end subroutine randomSketching
!!***

!----------------------------------------------------------------------

!!****f* m_slice/computeChebyshevMoments
!! NAME
!! computeChebyshevMoments
!! 
!! FUNCTION
!! Compute Chebyshev moments up to maximal degree all centered in [A,B)
!! Upper bound is computed as Rayleigh quotient
!!
!! OUTPUT
!! M_n = <X, f_n(B^{-1}AX) X> for n=1,..,ndeg_filter_max
!! 
!! SOURCE

subroutine computeChebyshevMoments(slice, X0, getAX_BX, getBm1X, &
        lambda_minus, lambda_plus, ndeg_filter, cheby_moments, maxeig_global)

    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X0
    integer, intent(in) :: ndeg_filter
    real(dp), intent(in) :: lambda_minus, lambda_plus
    real(dp), intent(out), optional :: maxeig_global
    complex(dp), intent(out) :: cheby_moments(:,:)
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

    integer :: neigenpairs, nband
    integer :: ideg
    integer :: iband
    integer :: tot_spacedim
    integer :: space_res, me_g0
    integer :: space, spacecom, gpu_option
    integer :: ierr
    logical :: compute_QR
    real(dp) :: maxeig, mineig
    real(dp) :: center, radius
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: tsec(2)
    type(xg_t) :: Moments
    type(xg_t) :: DivResults
    type(xg_t) :: X0_backup
    type(xgBlock_t) :: Moment_ideg
    type(chebfi_t) :: chebfi
    complex(dp), pointer :: momvals(:,:) => null()

    ! *********************************************************************

    ! todo add timer
    ! tim_cheby_moments

    space = slice%space
    spacecom = slice%spacecom
    neigenpairs = slice%neigenpairs
    tot_spacedim = slice%total_spacedim
    gpu_option = slice%gpu_option
    me_g0 = slice%me_g0
    if (slice%paral_kgb==1) then
        me_g0 = slice%me_g0_fft
    end if
    nband = cols(X0)
 
    compute_QR = .false.
    if (present(maxeig_global)) then
        compute_QR = .true.
    end if
    if (compute_QR) then
        if (slice%space==SPACE_C) then
            space_res = SPACE_C
        else if (slice%space==SPACE_CR) then
            space_res = SPACE_R
        else
            ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
        end if
        ! Workspace for Rayleigh quotients
        call xg_init(DivResults, space_res, nband, 1, gpu_option=gpu_option)
    end if
    
    ! Moment workspace size (1, ndeg+1)
    call xg_init(Moments, space, 1, ndeg_filter+1, gpu_option=gpu_option) ! M_n=<X0,f_n(A)X0>
    call xg_init(X0_backup, space, tot_spacedim, nband, spacecom, me_g0=me_g0, gpu_option=gpu_option) ! X0

    ! Initialize chebfi object in MPI Colsrows distribution
    call chebfi_init(chebfi,nband,tot_spacedim,slice%tolerance,slice%ecut,slice%paral_kgb,&
        nband,ndeg_filter,0,space,1,xmpi_comm_null,slice%me_g0,slice%me_g0_fft,&
        slice%paw,slice%comm_rows,slice%comm_cols,0,1.d0,0.d0,gpu_option,&
        gpu_kokkos_nthrd=slice%gpu_kokkos_nthrd,gpu_thread_limit=slice%gpu_thread_limit,&
        from_linalg=.false.)

    ! Initialize Chebyshev recursion
    call xgBlock_copy(X0, X0_backup%self)
    chebfi%xXColsRows = X0    

    ! Compute A*Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows, chebfi%xAXColsRows, chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ! Initialize Chebyshev moment at k=1
    call xgBlock_setBlock(Moments%self, Moment_ideg, 1, 1, fcol=1) 
    call xgBlock_dot(X0_backup%self, chebfi%xXColsRows, Moment_ideg)

    if (compute_QR) then
        ! Compute upper bound of interval as Rayleigh quotient
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RRQ)
        call timab(tim_RR_q, 1, tsec)
        call chebfi_rayleighRitzQuotients(chebfi, maxeig, mineig, DivResults%self)
        call timab(tim_RR_q, 2, tsec)
        ABI_NVTX_END_RANGE()
    
        call xmpi_max(maxeig, maxeig_global, spacecom, ierr)

        write(std_out,*) 'maxeig_global=', maxeig_global
        !write(std_out,*) 'divresults=', xgBlock_getid(DivResults%self)
        !write(std_out,*) 'X0=', xgBlock_getid(X0)
        !write(std_out,*) 'xX=', xgBlock_getid(chebfi%xXColsRows)
        !call xgBlock_print(DivResults%self, std_out)
        flush(std_out)
    end if

    ! Spectral interval to be amplified scaled to [-1,1)
    center = (lambda_plus + lambda_minus)/2.d0
    radius = (lambda_plus - lambda_minus)/2.d0 
    one_over_r = 1.d0/radius
    two_over_r = 2.d0/radius

    do ideg = 0, ndeg_filter - 1
     
        !chebfi%paw = .false.
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()
        !chebfi%paw = .true.

        ! chebfi%xXColsRows = f_ideg X0
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
        call timab(tim_swap,1,tsec)
        call chebfi_swapInnerBuffers(chebfi, tot_spacedim, nband)
        call timab(tim_swap,2,tsec)
        ABI_NVTX_END_RANGE()
        
        !A * Psi    
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows, chebfi%xAXColsRows, chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)

        ! M_ideg = < X0, f_ideg X0 >_B
        call xgBlock_setBlock(Moments%self, Moment_ideg, 1, 1, fcol=ideg+2) 
        call xgBlock_dot(X0_backup%self, chebfi%xXColsRows, Moment_ideg)

    end do 
    
    call xgBlock_reverseMap(Moments%self, momvals, 1, ndeg_filter+1)
    cheby_moments(:,:) = momvals(:,:)

    ! Free memory
    call xg_free(Moments)
    if (compute_QR) then
        call xg_free(DivResults)
    end if
    call chebfi_free(chebfi)
    call xg_free(X0_backup)
    
end subroutine computeChebyshevMoments
!!***

!----------------------------------------------------------------------

!!****f* m_slice/computeTraceEstimation
!! NAME
!! computeTraceEstimation
!! 
!! FUNCTION
!! Compute Girard-Hutchinson trace estimator
!! ndeg_filter -> number of Chebyshev moments
!! m_probe -> number of stochastic probes
!! moments(1:ndeg_filter) -> Chebyshev moments computed with probes
!! and B-inner product
!!
!! SOURCE

subroutine computeTraceEstimation(slice, getAX_BX, getBm1X, ndeg_filter, m_probe,&
        min_low_bound, xXColsRows)

    implicit none

    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: ndeg_filter
    integer, intent(in) :: m_probe
    real(dp), intent(in) :: min_low_bound
    type(xgBlock_t), optional, intent(inout) :: xXColsRows
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

    integer :: neigenpairs, nband
    integer :: ideg
    integer :: iband
    integer :: spacedim, tot_spacedim
    integer :: space_res
    integer :: space, spacecom, gpu_option
    integer :: ierr
    integer :: me_g0
    integer :: jcol
    integer :: my_rank
    integer :: seed
    integer :: m_probe_tot
    integer :: i, k, ib
    integer :: num_moments
    integer :: uppb_loc
    integer :: ngrid_fine, ngrid_coarse
    real(dp) :: ncount_ovlp
    real(dp) :: lambda_plus_wanted
    real(dp) :: b_init, b_ext, ecut
    real(dp) :: maxeig, mineig
    real(dp) :: center, radius
    real(dp) :: partial_mass
    real(dp) :: step_fine, step_coarse
    real(dp) :: b_scaled
    real(dp) :: deriv_ib
    type(xg_t) :: X_probe
    complex(dp), allocatable :: cheby_moments(:,:)
    real(dp), allocatable :: bgrid_fine(:)
    real(dp), allocatable :: bgrid_coarse(:)
    real(dp), allocatable :: moments(:)
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: cumm_eigen_count(:)
    real(dp) :: tsec(2)
    
    ! *********************************************************************
    
    call timab(tim_trace,1,tsec)

    ecut = slice%ecut
    spacecom = slice%spacecom
    spacedim = slice%spacedim
    neigenpairs = slice%neigenpairs
    tot_spacedim = slice%total_spacedim
    space = slice%space
    my_rank = xmpi_comm_rank(slice%spacecom)
    gpu_option = slice%gpu_option
    me_g0 = slice%me_g0
    if (slice%paral_kgb==1) then
        me_g0 = slice%me_g0_fft
    end if

    num_moments = ndeg_filter + 1
    ABI_MALLOC(moments, (num_moments))
    ABI_MALLOC(work, (num_moments))
    ABI_MALLOC(cheby_moments, (1, num_moments) )
    
    ngrid_coarse = 10 ! coarse, just to find uppb
    ngrid_fine = 30 ! used for cumulative eigenvalue count
    ABI_MALLOC(bgrid_coarse, (ngrid_coarse))
    ABI_MALLOC(bgrid_fine, (ngrid_fine))
    ABI_MALLOC(cumm_eigen_count, (ngrid_fine))

    ! total number of probes is m_probes * number of MPI processes
    call xg_init(X_probe, slice%space, tot_spacedim, m_probe, spacecom, &
        me_g0=me_g0, gpu_option=gpu_option)

    ! Define random isotropic probes (unit variance NOT unit norm!)
    do jcol = 1, m_probe
        seed = my_rank*(m_probe+10)+jcol ! seed depends on column index
        call xgBlock_colwiseRandomRademacher(X_probe%self, seed, jcol)
    end do

    write(std_out,*) 'moments are scaled in'
    write(std_out,*) min_low_bound, ecut
    flush(std_out)

    ! Lowpass scan - coarse resolution with low degree
    ! Chebyshev moments in maximal [a,ecut)
    ! todo <X_probe, f(A) X_probe> (trace) and <X0, f(A) X_probe> (principal angles)
    ! maybe <xX, f(A)X_probe> is useful for principal angles and column selection
    ! in that case incorporate it in the loop
    if (present(xXColsRows)) then
        call computeChebyshevMoments(slice, xXColsRows, getAX_BX, getBm1X, &
            min_low_bound, ecut, ndeg_filter, cheby_moments, b_init) ! debug
    else
        call computeChebyshevMoments(slice, X_probe%self, getAX_BX, getBm1X, &
            min_low_bound, ecut, ndeg_filter, cheby_moments, b_init)
    end if

    write(std_out,*) 'b_init=', b_init
    flush(std_out)

    ! Sum real part of moments and divide by number of probes
    ! moments(k) = 1/Nv * Sum_{i=1}^Nv v_i^T T_k(A)v_i
    m_probe_tot = m_probe
    call xmpi_sum(m_probe_tot, spacecom, ierr)
    !moments(1:num_moments) = (/ (real(cheby_moments(1,k)), k=1,num_moments) /)
    moments = real(cheby_moments(1,:))
    call xmpi_sum(moments, spacecom, ierr)
    moments(1:num_moments) = moments(1:num_moments)/m_probe_tot
    
    write(std_out,*) 'moments k=0=', real(cheby_moments(1,1))
    write(std_out,*) 'moments k=1=', real(cheby_moments(1,2))
    write(std_out,*) 'moments k=3=', real(cheby_moments(1,3))
    flush(std_out)

    center = (ecut + min_low_bound) / 2.d0
    radius = (ecut - min_low_bound) / 2.d0
   
    ! #########################################
    ! ########## Coarse resolution ############
    ! #########################################

    ! scan with lowpass
    ! first pass uppb is actually unknown
    step_coarse = (ecut - min_low_bound) / (ngrid_coarse - 1)
    bgrid_coarse = (/ ( min_low_bound + (ib-1)*step_coarse, ib=1,ngrid_coarse ) /)
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

    ! Fine grid resolution to find upper bound
    ! todo refinement if needed

    ! #########################################
    ! ########### Fine resolution #############
    ! #########################################

    ! Now compute eigenvalue count
    step_fine = (lambda_plus_wanted - min_low_bound) / (ngrid_fine - 1)
    bgrid_fine = (/ (min_low_bound + (ib-1)*step_fine, ib=1,ngrid_fine) /) 
    deriv_ib = 0.d0
    do ib=1, ngrid_fine
        b_scaled = (bgrid_fine(ib) - center) / radius
        partial_mass = get_eigenvalue_count(b_scaled, moments, work)
        cumm_eigen_count(ib) = partial_mass
        if (ib>1) then
            deriv_ib = partial_mass - cumm_eigen_count(ib-1)
        end if 
        write(std_out,*) ib, 'scan: <=', bgrid_fine(ib), 'mass=', partial_mass, 'deriv=', deriv_ib
    end do

    ! #########################################
    ! ########### Final decision  #############
    ! #########################################

    ! todo integrate this procedure in the refinement
    ! like refine until target is reached then find maximal overlap staying in the gap etc

    ! Bound placement
    ! Constraints: 
    ! 1) number of bands per slice balanced
    write(std_out,*) 'target bands per slice=', slice%neigenpairs/slice%nslice
    flush(std_out)

    ! 2) gap: bound is located at a region where deriv is zero 

    ! 3) no cluster is present after the bound, like the next bound does not contain a lot eigs

    slice%neigenpairs_per_slice(1) = get_eigenvalue_count((slice%poly_upp_bounds(1)-center)/radius, moments, work)

    ncount_ovlp = get_eigenvalue_count((slice%poly_low_bounds(2)-center)/radius, moments, work)

    slice%neigenpairs_per_slice(2) = neigenpairs - 2*slice%neigenpairs_per_slice(1) + ncount_ovlp

    write(std_out,*) 'testing ncount 1=', slice%neigenpairs_per_slice(1)
    write(std_out,*) 'testing ncount 2=', slice%neigenpairs_per_slice(2)
    flush(std_out)

    !do k=1, slice%nslice
        ! find the closest zero

        ! interval limits without overlap
        !slice%part_low_bounds(k) = 
        !slice%part_upp_bounds(k) = 

        ! interval limits with overlap
        !slice%poly_low_bounds(k) = 
        !slice%poly_upp_bounds(k) = 

        ! count with overlap
        !slice%neigenpairs_per_slice(k) =

    !end do

    ! Free memory
    ABI_FREE(cheby_moments)
    ABI_FREE(moments)
    ABI_FREE(bgrid_fine)
    ABI_FREE(bgrid_coarse)
    ABI_FREE(work)
    ABI_FREE(cumm_eigen_count)

    call xg_free(X_probe)
    
    call timab(tim_trace,2,tsec)
    
end subroutine computeTraceEstimation
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

        comm_ = xmpi_comm_null
        if (present(comm)) then
            comm_ = comm
        end if
        neigenpairs = size(cheby_moments, 1)
        ndeg_filter = size(cheby_moments, 2)-1

        ABI_MALLOC(energy_per_band, (neigenpairs))
        ABI_MALLOC(cja, (ndeg_filter + 1))
        ABI_MALLOC(spectral_mass_left, (nstep_bisect))
        ABI_MALLOC(spectral_mass_right, (nstep_bisect))

        ! Step 1: increase b until all mass is included            
        b_ext = b
        call buildChebyshevJacksonCoeffs((a-center)/radius, (b_ext-center)/radius, ndeg_filter, cja)
        call computeFilterEnergy(cja, cheby_moments, energy_per_band)
        tot_mass = norm2(energy_per_band)**2
        call xmpi_sum(tot_mass, comm, ierr)
        iext_step = 0
        width_ext = (b - a)/12.d0
        write(std_out,*) 'width_ext=', width_ext
        flush(std_out)
        do while (tot_mass < nband_tot .and. iext_step < 100)
            b_ext = b_ext + width_ext
            call buildChebyshevJacksonCoeffs((a-center)/radius, (b_ext-center)/radius, ndeg_filter, cja)
            call computeFilterEnergy(cja, cheby_moments, energy_per_band)
            tot_mass = sum(energy_per_band)
            call xmpi_sum(tot_mass, comm, ierr)
            iext_step = iext_step + 1
            write(std_out,*) 'extended interval to=', b_ext, 'nvec=', tot_mass
            flush(std_out)
        end do

        width = (b_ext - a) / (nstep_bisect + 1)
        mass_diff_prev = 0.d0
        c_split = a
        ! TODO verify that it converges towards a limit that is upper bound of b_true

        do ishift = 1, nstep_bisect
    
            c = a + ishift * width
           
            write(std_out,*) '======================================'
            write(std_out,*) 'ishift=', ishift, 'c=', c
            flush(std_out)

            ! Slice Left [a,c)
            call buildChebyshevJacksonCoeffs((a-center)/radius, (c-center)/radius, ndeg_filter, cja)
            call computeFilterEnergy(cja, cheby_moments, energy_per_band)
            mass_left = sum(energy_per_band)
            write(std_out,*) 'energy_per_band=', energy_per_band; flush(std_out)
            call xmpi_sum(mass_left, comm, ierr)
            spectral_mass_left(ishift) = mass_left

            write(std_out,*) '                      Nvec left=', mass_left
            flush(std_out)

            ! Slice Right [c,b)
            call buildChebyshevJacksonCoeffs((c-center)/radius, (b_ext-center)/radius, ndeg_filter, cja)
            call computeFilterEnergy(cja, cheby_moments, energy_per_band)
            write(std_out,*) 'energy_per_band=', energy_per_band; flush(std_out)
            mass_right = sum(energy_per_band)
            call xmpi_sum(mass_right, comm_, ierr)
            spectral_mass_right(ishift) = mass_right
           
            write(std_out,*) '                      Nvec right=', mass_right
            flush(std_out)
            
            ! Detect spectral mass flip
            mass_diff = mass_left - mass_right
            if (ishift > 1 .and. mass_diff * mass_diff_prev < 0) then
                c_split = (2 * a + (2*ishift - 1) * width) / 2.d0
            end if
            mass_diff_prev = mass_diff

        end do

        ABI_FREE(energy_per_band)
        ABI_FREE(spectral_mass_left)
        ABI_FREE(spectral_mass_right)
        ABI_FREE(cja)

end subroutine splitSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/spectralPruning
!! NAME
!! spectralPruning
!! 
!! FUNCTION
!! Selects max probes in [a,b) mapped to [-1,1) with center and radius.
!! Returns indices out of tot_nband to keep of size wanted_mass.
!! 
!! SOURCE

subroutine spectralPruning(a, b, center, radius, cheby_moments, tot_nband, idx, comm)

        implicit none

        real(dp), intent(in) :: a, b
        real(dp), intent(in) :: center, radius
        integer, intent(in) :: tot_nband
        integer, intent(out) :: idx(:)
        complex(dp), intent(in) :: cheby_moments(:,:)
        integer, intent(in), optional :: comm

        integer :: neigenpairs
        integer :: ndeg_filter
        integer :: ishift
        integer :: wanted_mass
        integer :: comm_, my_rank, ierr
        integer :: iband
        real(dp) :: tol12 = 1.0e-12
        integer, allocatable :: jperm(:)
        real(dp), allocatable :: energy_per_band(:)
        real(dp), allocatable :: energy_per_band_global(:)
        real(dp), allocatable :: cja(:)

        ! *********************************************************************

        comm_ = xmpi_comm_null
        if (present(comm)) then
            comm_ = comm
        end if
        neigenpairs = size(cheby_moments, 1)
        ndeg_filter = size(cheby_moments, 2)-1
        wanted_mass = size(idx)

        ABI_MALLOC(jperm, (tot_nband))
        ABI_MALLOC(energy_per_band, (neigenpairs))
        ABI_MALLOC(energy_per_band_global, (tot_nband))
        ABI_MALLOC(cja, (ndeg_filter + 1))

        call buildChebyshevJacksonCoeffs((a-center)/radius, (b-center)/radius, ndeg_filter, cja)
        call computeFilterEnergy(cja, cheby_moments, energy_per_band)
      
        energy_per_band_global = 0.d0
        if (xmpi_comm_size(comm_) > 1) then
            ! sum contributions across procs
            my_rank = xmpi_comm_rank(comm_)
            ishift = my_rank * neigenpairs
            energy_per_band_global(ishift+1:ishift+neigenpairs) = energy_per_band(:)
            call xmpi_sum(energy_per_band_global, comm_, ierr)
        else
            energy_per_band_global(:) = energy_per_band(:)
        end if
    
        jperm = (/ (iband, iband=1, tot_nband) /)
        call sort_dp(tot_nband, energy_per_band_global, jperm, tol12)
        idx(1:wanted_mass) = jperm(tot_nband - wanted_mass + 1:tot_nband) 

        ABI_FREE(jperm)
        ABI_FREE(energy_per_band)
        ABI_FREE(energy_per_band_global)
        ABI_FREE(cja)

end subroutine spectralPruning
!!***

!----------------------------------------------------------------------

!!****f* m_slice/computeFilterEnergy
!! NAME
!! computeFilterEnergy
!! 
!! SOURCE

subroutine computeFilterEnergy(cja, cheby_moments, energy_per_band)

      implicit none

      real(dp), intent(in) :: cja(:)
      complex(dp), intent(in) :: cheby_moments(:,:)
      real(dp), intent(out) :: energy_per_band(:)
 
      complex(dp), allocatable :: energy(:)
      integer :: neigenpairs, ndeg_filter
      integer :: ideg, j

      ! *********************************************************************

      neigenpairs = size(cheby_moments,1)
      ndeg_filter = size(cheby_moments,2)-1
      ABI_MALLOC(energy, (neigenpairs))
      energy = dcmplx(0.0d0,0.0d0)
      !$omp parallel do private(ideg)
      do j = 1, neigenpairs
        do ideg = 1, ndeg_filter+1
            energy(j) = energy(j) + cja(ideg) * cheby_moments(j, ideg)
        end do
      end do
      !$omp end parallel do
      !$omp parallel do
      do j = 1, neigenpairs
          energy_per_band(j) = real(energy(j))
      end do
      !$omp end parallel do
      ABI_FREE(energy)

end subroutine computeFilterEnergy
!!***

!----------------------------------------------------------------------

!!****f* m_slice/jackson_step_coeffs
!! NAME
!! jackson_step_coeffs
!! 
!! FUNCTION
!! Jackson damped coefficients
!! 
!! SOURCE

  subroutine jackson_step_coeffs(a,b,lmin,lmax,deg,ctilde)
      
      implicit none

      integer, intent(in) :: deg
      real(dp), intent(in) :: a,b,lmin,lmax
      real(dp), intent(out) :: ctilde(1:deg+1)
      integer :: k
      real(dp) :: c, r, a_, b_, ck, gk, alpha, cotv

      c  = (lmin + lmax)/2.d0
      r  = (lmax - lmin)/2.d0
      a_ = (a - c)/r
      b_ = (b - c)/r
      alpha = pi/(deg+2)
      cotv = cos(alpha)/sin(alpha)

      do k = 0, deg
         ! Chebyshev coefficient for step function
         if (k == 0) then
            ck = (acos(a_) - acos(b_))/pi
         else
            ck = 2.0d0/pi*(sin(k*acos(a_)) - sin(k*acos(b_)))/k
         end if

         ! Jackson damping factor
         gk = ((deg - k + 1)*cos(pi*k/(deg+1)) + sin(pi*k/(deg+1))*cotv) / (deg+1)

         ctilde(k+1) = ck * gk
      end do

  end subroutine jackson_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_slice/erf_step_coeffs
!! NAME
!! erf_step_coeffs
!! 
!! FUNCTION
!! Erf damped coefficients
!! 
!! SOURCE

  function erf_step_coeffs(b, ndeg, sigma, Ngrid) result(coeffs)

      implicit none
      real(dp), intent(in) :: b, sigma
      integer, intent(in) :: ndeg, Ngrid

      real(dp) :: coeffs(ndeg+1)
      integer :: i, k
      real(dp) :: x(Ngrid+1)
      real(dp) :: f(Ngrid+1)

      x = (/ ( cos(Pi*(i-2)/(Ngrid-1) ) , i=1,Ngrid+1) /)
      f = (/ ( 0.5 * (1.0 - fast_erf((x(i) - b)/sigma)), i=1,Ngrid+1) /)
      
      do k=0, ndeg
        coeffs(k+1) = (2.0/Ngrid) * sum( (/ (f(i)*cos(k*acos(x(i))), i=1,Ngrid+1 ) /) )
      end do
      coeffs(1) = coeffs(1) / 2.d0 

  end function erf_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_slice/fast_erf
!! NAME
!! fast_erf
!! 
!! FUNCTION
!! Fast approximation of erf(x) using Abramowitz & Stegun 7.1.26

  function fast_erf(x) result(erf_val)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: erf_val
    real(dp) :: t, tau, ax
    real(dp), parameter :: p  = 0.3275911_dp
    real(dp), parameter :: a1 = 0.254829592_dp
    real(dp), parameter :: a2 = -0.284496736_dp
    real(dp), parameter :: a3 = 1.421413741_dp
    real(dp), parameter :: a4 = -1.453152027_dp
    real(dp), parameter :: a5 = 1.061405429_dp

    ax = abs(x)
    t = 1.0_dp / (1.0_dp + p * ax)
    tau = (((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t) * exp(-ax*ax)
    erf_val = 1.0_dp - tau
    if (x < 0.0_dp) erf_val = -erf_val
  end function fast_erf
!!***

!----------------------------------------------------------------------

!!****f* m_slice/smooth_step_coeffs
!! NAME
!! smooth_step_coeffs
!! 
!! FUNCTION
!! Smooth damped coefficients
!! 
!! SOURCE

  function smooth_step_coeffs(b, ndeg, alpha, Ngrid) result(coeffs)

      implicit none
      real(dp), intent(in) :: b, alpha
      integer, intent(in) :: ndeg, Ngrid

      real(dp) :: coeffs(ndeg+1)
      integer :: i, k
      real(dp) :: x(Ngrid+1)
      real(dp) :: f(Ngrid+1)

      x = (/ ( cos(Pi*(i-2)/(Ngrid-1) ) , i=1,Ngrid+1) /)
      f = (/ ( 0.5 * (1.0 - tanh( alpha*(x(i) - b) )) , i=1,Ngrid+1) /)
      
      do k=0, ndeg
        coeffs(k+1) = (2.0/Ngrid) * sum( (/ (f(i)*cos(k*acos(x(i))), i=1,Ngrid+1 ) /) )
      end do
      coeffs(1) = coeffs(1) / 2.d0 

  end function smooth_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_slice/lanczos_step_coeffs
!! NAME
!! lanczos_step_coeffs
!! 
!! FUNCTION
!! Chebyshev coefficients with Lanczos damping
!! 
!! SOURCE

function lanczos_step_coeffs(b, ndeg) result(coeffs)

    implicit none
    real(dp), intent(in) :: b
    integer, intent(in) :: ndeg

    real(dp) :: coeffs(ndeg+1)
    integer :: j
    real(dp) :: x(ndeg+1)
    real(dp) :: f(ndeg+1)
    real(dp) :: k_array(ndeg+1)

    x = (/ ( cos(Pi*(2*j-1)/(2*(ndeg+1)) ) , j=1,ndeg+1) /)
    f = merge(1.0, 0.0, x < b) ! f=1 if x<b else 0
    k_array = (/ (j, j=1,ndeg+1) /)
    do j = 1, ndeg+1
        coeffs(j) = 2.0d0 / (ndeg+1) * &
            sum( f(:) * cos( pi*(j-1)*(2.0d0*k_array(:)-1.0d0) / (2.0d0*(ndeg+1)) ) )
    end do
    coeffs(1) = coeffs(1) / 2.d0 
    do j = 2, ndeg+1
        coeffs(j) = coeffs(j) * sin(pi*(j-1)/(ndeg+1)) / (pi*(j-1)/(ndeg+1))
    end do 

end function lanczos_step_coeffs
!!***

!----------------------------------------------------------------------

!!****f* m_slice/get_eigenvalue_count
!! NAME
!! get_eigenvalue_count
!! 
!! SOURCE

function get_eigenvalue_count(b, moments, work) result(mass)

    implicit none

    real(dp), intent(in) :: b
    real(dp), intent(in) :: moments(:)
    real(dp), intent(inout) :: work(:)
    real(dp) :: mass
    integer :: ndeg_filter
    real(dp) :: sigma, alpha
    integer :: Ngrid

    ndeg_filter = size(moments)-1

    ! Erf damping coefficients
    !sigma = 4.d0 / ndeg_filter
    !Ngrid = 500
    !work = erf_step_coeffs(b, ndeg_filter, sigma, Ngrid)

    !alpha = 20
    !Ngrid = 500
    !work = smooth_step_coeffs(b, ndeg_filter, alpha, Ngrid)
    
    ! Steep Lanczos
    work = lanczos_step_coeffs(b, ndeg_filter)

    mass = dot_product(work, moments)

end function get_eigenvalue_count
!!***

end module m_slice
!!***

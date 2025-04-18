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
    use m_io_tools, only : flush_unit

    use m_cgtools
    use m_xg
    use m_xgTransposer
    use m_xg_ortho_RR
    use m_chebfi2

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
        integer :: me_id_slice                              ! identifier of slice task in use
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
        paral_slice,ndeg_filter,ramp,ecut,bandpp,space,spacecom,me_g0,me_g0_fft,&
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
    call slice_queryHostDevice(slice)

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


    if (slice%nslice /= 1 .and. slice%me_comm_slice /= slice%spacecom) then
        call xmpi_comm_free(slice%me_comm_slice)
    end if

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

subroutine slice_allschedule(slice, X0, getAX_BX, eigen, nspinor)

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

    ! Local variables --------------------------------
    integer :: neigenpairs, iband, min_loc, islice
    integer :: fcol, fcol_ext, ncols
    logical :: on_host, on_device
    real(dp) :: lambda_minus, lambda_plus
    real(dp) :: tol12 = 1.0e-12
    ! Derived types
    type(xg_t) :: resid0
    type(xgBlock_t) :: eigen_sorted
    type(xgBlock_t) :: slicecols_in
    type(xgBlock_t) :: slicecols_ext_out
    ! Arrays
    integer, allocatable, target :: permute_cols(:)
    real(dp), allocatable, target :: theta_reshaped(:)
    integer, pointer :: permute_cols_ptr(:) => null()
    real(dp), pointer :: theta_(:,:) => null()
    real(dp), pointer :: resid_(:,:) => null()
    real(dp), pointer :: theta_reshaped_ptr(:) => null()
    
    ! *********************************************************************

    ABI_NVTX_START_RANGE(NVTX_SLICE_SCHEDULE)
    
    ! Sanity check for X0 (on GPU): verify we are in target enter data map
    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call slice_queryHostDevice(slice, on_host, on_device)
        if (on_host .or. (.not. on_device)) then
            ABI_ERROR("not in target data map region")
        end if
    end if

    neigenpairs = slice%neigenpairs
    
    ! Space for computed residuals (not distributed)
    call xg_init(resid0, SPACE_R, rows=1, cols=neigenpairs, gpu_option=slice%gpu_option)

    ABI_MALLOC_IFNOT(theta_reshaped, (neigenpairs))
    ABI_MALLOC_IFNOT(permute_cols, (neigenpairs))

    call xgBlock_reshape(eigen, 1, neigenpairs)
    call xgBlock_zero(eigen)
    call xgBlock_zero(resid0%self)
    
    ! ===================== Compute Rayleigh quotients and residuals ===================================
    
    ABI_NVTX_START_RANGE(NVTX_SLICE_RRQ)
    call slice_computeSpectrum(slice, X0, getAX_BX, eigen, resid0%self, nspinor)
    ABI_NVTX_END_RANGE()

    ! ===================== Compute guaranteed spectral bounds ======================================== 

    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call xgBlock_copy_from_gpu(eigen)
        call xgBlock_copy_from_gpu(resid0%self)
    end if

    ! Results could be complex, so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(eigen, theta_, rows=1, cols=neigenpairs)
    call xgBlock_reverseMap(resid0%self, resid_, rows=1, cols=neigenpairs)
 
    ! Sort thetas in increasing order and store result to eigen
    theta_reshaped(1:neigenpairs) = theta_(1,1:neigenpairs)
    permute_cols_ptr => permute_cols
    permute_cols(1:neigenpairs) = (/ (iband, iband=1,neigenpairs) /)
    call sort_dp(neigenpairs, theta_reshaped, permute_cols_ptr, tol12)
    theta_(1,1:neigenpairs) = theta_reshaped(1:neigenpairs)
    call xgBlock_map(eigen_sorted, theta_, SPACE_R, rows=1, cols=neigenpairs, gpu_option=slice%gpu_option)
    call xgBlock_copy(eigen_sorted, eigen)

    ! Minimum and maximum quotients
    lambda_minus = theta_reshaped(1)
    lambda_plus = theta_reshaped(neigenpairs)

    ! Guaranteed spectral bounds
    min_loc = permute_cols(1)
    slice%mineig_global = lambda_minus - sqrt(resid_(1, min_loc))
    slice%maxeig_global = slice%ecut

    ! ===================== Split interval [lambda_minus,lambda_plus) into slices ======================
    
    if (slice%nslice==1) then
        slice%neigenpairs_per_slice = slice%neigenpairs
        slice%fcol_in_X = 1
        slice%fcol_in_Xext = 1
        slice%poly_degrees = slice%ndeg_filter
        slice%part_low_bounds = slice%mineig_global
        slice%part_upp_bounds = slice%maxeig_global
        slice%poly_low_bounds = lambda_minus
        slice%poly_upp_bounds = lambda_plus
    else
        theta_reshaped_ptr => theta_reshaped
        call slice_cutSpectrum(slice, lambda_minus, lambda_plus, theta_reshaped_ptr, plot_filter=.false.)
    end if

    ! ===================== Resource management system =================================================
    
    ! Divide resources into slice tasks
    if (slice%nslice==1) then
        slice%nproc_per_slice = xmpi_comm_size(slice%spacecom)
        slice%lookup_proc = 0
    else
        call slice_allocateResources(slice)
    end if

    ! Run on all ranks of spacecom: Mark my slice task and resources as actively in use
    call slice_markActiveTask(slice)

    ! ===================== Allocate and fill extended memory buffer ================================== 
  
    ! Sanity check
    if ((.not. slice%use_linalg) .or. slice%use_colsrows) then
        ABI_ERROR("not in linalg representation")
    end if

    ! Permute column vectors in Rayleigh quotient increasing order
    if (slice%nslice>1) then
        call xgBlock_permuteCols(X0, slice%total_spacedim, neigenpairs, permute_cols_ptr)
    end if

    slice%neigenpairs_ext = sum(slice%neigenpairs_per_slice)

    if (slice%nslice==1) then
        slice%XextLinalg = X0
    else 
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
    if (cols(slice%XextLinalg) /= slice%neigenpairs_ext) then
        ABI_ERROR('wrong linalg representation')
    end if
    write(std_out,'(a,i6,i6)') '# proc has # cols of Xext_linalg ', xmpi_comm_rank(slice%spacecom), cols(slice%XextLinalg)

    write(std_out,*) 'X0', xgBlock_getid(X0)
    write(std_out,*) 'slice%XextLinalg', xgBlock_getid(slice%XextLinalg)
    
    ! Recover dimensions
    call xgBlock_reshape(eigen, neigenpairs, 1)

    ! Free temporary memory
    call xg_free(resid0) 
    ABI_SFREE(permute_cols)
    ABI_SFREE(theta_reshaped)

    ABI_NVTX_END_RANGE()

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

    ! Variables
    type(chebfi_t) :: chebfi
    type(xgBlock_t) :: X0_active
    type(xgBlock_t) :: eigen_active
    type(xgBlock_t) :: residu_active
    integer :: nbdbuf, oracle, num_proc
    integer :: neigenpairs, bandpp, ndeg_filter, comm
    real(dp) :: lambda_minus, lambda_plus
    real(dp) :: oracle_factor, oracle_min_occ
    logical :: is_lowpass, on_host, on_device
    ! Arrays
    integer, allocatable, target :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null() 
    integer, pointer :: ncolsColsRows_ptr(:) => null()

    ! *********************************************************************
    
    ! Sanity check
    if ( (.not. slice%use_linalg) .and. slice%use_colsrows) then
        ABI_ERROR("should be in linalg representation")
    end if
    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call slice_queryHostDevice(slice, on_host, on_device)
        ABI_CHECK(on_device,"GPU not used when it should be!")
    end if
    
    write(std_out,*) 'slice%XextLinalg', xgBlock_getid(slice%XextLinalg)

    ! Distribute extended columns across **all** MPI processes
    ! After the transposition each process contains the correct
    ! bandpp corresponding to the slice so that no additional communication
    ! has to be performed in order to bring band slices to processes.
    ncolsColsRows_ptr => slice%ncolsColsRows

    write(std_out,*) ncolsColsRows_ptr

    ! Allocate slice%me_Xext_active according to the target MPI distribution for slices
    call xgTransposer_constructor(slice%xgTransposerXext, slice%XextLinalg, slice%me_Xext_active,&
        nspinor, STATE_LINALG, TRANS_ALL2ALL, slice%comm_rows, slice%comm_cols, 0, 0, slice%me_g0_fft,&
        gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit,&
        custom_ncolsColsRows=.true., ncolsColsRows_sub=ncolsColsRows_ptr)
   
    slice%xgTransposerXext%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
    
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(slice%xgTransposerXext, STATE_COLSROWS)
    ABI_NVTX_END_RANGE()
     
    write(std_out,*) 'slice%me_Xext_active', xgBlock_getid(slice%me_Xext_active)
    ! to be compared with chebfi_run

    slice%use_colsrows = .true.
    slice%use_linalg = .false.

    ! Unitary test
    if ( cols(slice%me_Xext_active) /= slice%ncolsColsRows(xmpi_comm_rank(slice%spacecom)+1) ) then
        ABI_ERROR('wrong colsrows representation')
    end if
    write(std_out,'(a,i6,i6,i6)') '# proc has # cols of Xext ', xmpi_comm_rank(slice%spacecom),&
        cols(slice%me_Xext_active)

    ! Get parameters of active task
    neigenpairs = slice%me_neigenpairs_slice
    bandpp = slice%me_bandpp_slice
    ndeg_filter = slice%me_ndeg_slice
    comm = slice%me_comm_slice 
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

    ! Initialize chebfi object in MPI Colsrows distribution
    call chebfi_init(chebfi,neigenpairs,slice%total_spacedim,slice%tolerance,slice%ecut,slice%paral_kgb,bandpp,&
        ndeg_filter,nbdbuf,slice%space,1,comm,slice%me_g0,slice%me_g0_fft,slice%paw,xmpi_comm_self,comm,&
        oracle,oracle_factor,oracle_min_occ,slice%gpu_option,gpu_kokkos_nthrd=slice%gpu_kokkos_nthrd,&
        gpu_thread_limit=slice%gpu_thread_limit,from_linalg=.false.)
 
    ! Define pointers to actively used arrays
    call xmpi_barrier(slice%spacecom)
    X0_active = slice%me_Xext_active
    call xgBlock_setBlock(eigen, eigen_active, rows=neigenpairs, cols=1)
    call xgBlock_setBlock(residu, residu_active, rows=neigenpairs, cols=1)
    
    write(std_out,*) 'X0_active', xgBlock_getid(X0_active)
    
    ! Restrict to sub-communicator
    call xgBlock_setComm(X0_active, comm)
    call xgBlock_setComm(slice%me_Xext_active, comm)
    !call xgBlock_setComm(eigen_active, comm)
    !call xgBlock_setComm(residu_active, comm)
    
    write(std_out,*) 'X0_active', xgBlock_getid(X0_active)

    call chebfi_runSlice(chebfi, X0_active, getAX_BX, getBm1X, eigen_active, residu_active, nspinor,&
        slice%mineig_global, slice%maxeig_global, lambda_minus, lambda_plus, is_lowpass, nrowsLinalg_ptr)

    ! Restore global comm
    !call xgBlock_setComm(slice%me_Xext_active, slice%spacecom)
    !call xgBlock_setComm(eigen, slice%spacecom)
    !call xgBlock_setComm(residu, slice%spacecom)

    ! Free temporary memory
    call chebfi_free(chebfi)
    ABI_SFREE(nrowsLinalg)

    ! Sanity check
    if ( (.not. slice%use_colsrows) .and. slice%use_linalg) then
        ABI_ERROR("should be in colsrows representation")
    end if
    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call slice_queryHostDevice(slice, on_host, on_device)
        ABI_CHECK(on_device,"GPU not used when it should be!")
    end if

    ! Actually do the transposition to linalg
    call xmpi_barrier(slice%spacecom)
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(slice%xgTransposerXext, STATE_LINALG)
    ABI_NVTX_END_RANGE()
    ! Note: At this point slice%me_Xext_active is recovered into slice%XextLinalg

    write(std_out,*) 'slice%XextLinalg', xgBlock_getid(slice%XextLinalg)

    slice%use_colsrows = .false.
    slice%use_linalg = .true.

end subroutine slice_run
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_computeSpectrum
!! NAME
!! slice_computeSpectrum
!!
!! FUNCTION
!! Compute Rayleigh quotients and residuals
!! 
!! INPUTS
!! X        =eigenvector guess
!! getAX_BX =Hamiltonian application
!! 
!! OUTPUT
!! eigen  =Rayleigh Quotients from X, size (neigenpairs,1)
!! residu =residuals from X, size (neigenpairs,1)
!! 
!! SOURCE

subroutine slice_computeSpectrum(slice, X, getAX_BX, eigen, resid, nspinor)

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

    ! Local variables
    ! Scalars
    integer :: my_rank, shift, bandpp, space_res
    real(dp) :: mineig, maxeig
    ! Derived types
    type(xg_t) :: Results1, Results2
    type(xg_t) :: eigen_mpi, resid_mpi
    type(xg_t) :: X_NAB ! vector memory for X_next, AX, BX
    type(xgBlock_t) :: eigen_block, resid_block
    type(xgBlock_t) :: xAXColsRows
    type(xgBlock_t) :: xBXColsRows
    type(xgBlock_t) :: X_next
    type(xgBlock_t) :: xXColsRows
    type(xgBlock_t) :: eigen_mpi_reshaped
    type(xgTransposer_t) :: xgTransposerX
    ! Arrays
    real(dp), allocatable, target :: theta_mpi_reshaped(:)
    real(dp), pointer :: theta_mpi_reshaped_ptr(:) => null()
    real(dp), pointer :: theta_mpi(:,:) => null()
    integer :: maxeig_pos(2)
    integer :: mineig_pos(2)

    ! *********************************************************************

    if (slice%paral_kgb==0) then
        ABI_ERROR("Sequential bands not implemented")
    end if

    ! Space of eigenvalues
    if (slice%space==SPACE_C) then
        space_res = SPACE_C
    else if (slice%space==SPACE_CR) then
        space_res = SPACE_R
    end if
    bandpp = slice%bandpp

    ! Allocate temporary memory (distributed in colsrows representation)
    call xg_init(X_NAB, slice%space, slice%total_spacedim, 3*bandpp, slice%spacecom, &
        me_g0=slice%me_g0, gpu_option=slice%gpu_option)

    call xg_setBlock(X_NAB, X_next, slice%total_spacedim, bandpp)                           ! X_next
    call xg_setBlock(X_NAB, xAXColsRows, slice%total_spacedim, bandpp, fcol=bandpp + 1)     ! xAXColsRows
    call xg_setBlock(X_NAB, xBXColsRows, slice%total_spacedim, bandpp, fcol=2*bandpp + 1)   ! xBXColsRows
    
    ! Allocate one-dimensional memory (distributed)
    ! using space_res
    call xg_init(Results1, space_res, rows=bandpp, cols=1, gpu_option=slice%gpu_option)
    call xg_init(Results2, space_res, rows=bandpp, cols=1, gpu_option=slice%gpu_option)
    call xg_init(eigen_mpi, space_res, rows=bandpp, cols=1, comm=slice%spacecom, gpu_option=slice%gpu_option)
    ! using SPACE_R
    call xg_init(resid_mpi, SPACE_R, rows=bandpp, cols=1, comm=slice%spacecom, gpu_option=slice%gpu_option)
 
    ! Allocate memory for X in colsrows representation
    call xgTransposer_constructor(xgTransposerX, X, xXColsRows, nspinor, STATE_LINALG,&
        TRANS_ALL2ALL, xmpi_comm_self, slice%spacecom, 0, 0, slice%me_g0_fft,&
        gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit)
        
    xgTransposerX%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
        
    ! Sanity check
    if ((.not. slice%use_linalg) .or. slice%use_colsrows) then
        ABI_ERROR("not in linalg")
    end if

    ! ============== Transpose ==============
    call xmpi_barrier(slice%spacecom)
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(xgTransposerX, STATE_COLSROWS)
    ABI_NVTX_END_RANGE()

    slice%use_linalg = .false.
    slice%use_colsrows = .true.

    ! Now apply A and B to X (requires colsrows representation) to create AX and BX in colsrows
    ! Remember that this function will copy X to BX if paw
    ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
    call getAX_BX(xXColsRows, xAXColsRows, xBXColsRows)
    call xgBlock_zero_im_g0(xAXColsRows)
    call xgBlock_zero_im_g0(xBXColsRows)
    ABI_NVTX_END_RANGE()

    ! Compute Rayleigh quotients
    ! <Psi|H|Psi>
    call xgBlock_colwiseDotProduct(xXColsRows, xAXColsRows, Results1%self, comm_loc=xmpi_comm_null)
    ! <Psi|S|Psi>
    call xgBlock_colwiseDotProduct(xXColsRows, xBXColsRows, Results2%self, comm_loc=xmpi_comm_null)
    ! eigen = <Psi|H|Psi> / <Psi|S|Psi>
    call xgBlock_colwiseDivision(Results1%self, Results2%self, eigen_mpi%self, &
        maxeig, maxeig_pos, mineig, mineig_pos)

    ! In order to avoid transposing AX,BX, we compute residuals in colsrows representation
    ! TODO IL 20/01/2025 ymax has not been tested on GPU
    call xgBlock_copy(xBXColsRows, X_next)                             ! X_next = S|Psi>
    call xgBlock_ymax(X_next, eigen_mpi%self, 0, 1)                    ! X_next = - eig * S|Psi>
    call xgBlock_add(X_next, xAXColsRows)                              ! X_next = H|Psi> - eig * S|Psi>
    call xgBlock_colwiseNorm2(X_next, resid_mpi%self, comm_loc=xmpi_comm_null)  ! resid = |X_next|^2
   
    ! MPI communication for gathering all thetas (using summation strategy on columns)
    my_rank = xmpi_comm_rank(slice%spacecom)
    shift = my_rank * bandpp
    ! for eigen
    call xgBlock_setBlock(eigen, eigen_block, rows=1, cols=bandpp, fcol=1+shift)
    call xgBlock_reshape(eigen_mpi%self, 1, bandpp)
    if (space_res==SPACE_R) then     
        call xgBlock_copy(eigen_mpi%self, eigen_block)
        call xgBlock_mpi_sum(eigen, comm=slice%spacecom)
    else
        ! workaround to copy from SPACE_C to SPACE_R
        ABI_MALLOC_IFNOT(theta_mpi_reshaped,(bandpp))
        theta_mpi_reshaped_ptr => theta_mpi_reshaped
        call xgBlock_reverseMap(eigen_mpi%self, theta_mpi, rows=1, cols=bandpp)
        theta_mpi_reshaped(1:bandpp) = theta_mpi(1,1:bandpp)
        call xgBlock_map_1d(eigen_mpi_reshaped, theta_mpi_reshaped_ptr, SPACE_R, bandpp, gpu_option=slice%gpu_option)
        call xgBlock_copy(eigen_mpi_reshaped, eigen_block)
        call xgBlock_mpi_sum(eigen, comm=slice%spacecom)
        ABI_SFREE(theta_mpi_reshaped)
    end if
    ! for resid (SPACE_R)
    call xgBlock_setBlock(resid, resid_block, rows=1, cols=bandpp, fcol=1+shift)
    call xgBlock_reshape(resid_mpi%self, 1, bandpp)     
    call xgBlock_copy(resid_mpi%self, resid_block)
    call xgBlock_mpi_sum(resid, comm=slice%spacecom)

    ! ============== Transpose ==============
    call xmpi_barrier(slice%spacecom)
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(xgTransposerX, STATE_LINALG)
    ABI_NVTX_END_RANGE()

    slice%use_linalg = .true.
    slice%use_colsrows = .false.

    ! Free temporary memory
    call xg_free(Results1)
    call xg_free(Results2)
    call xg_free(eigen_mpi)
    call xg_free(resid_mpi)
    call xg_free(X_NAB)
    call xgTransposer_free(xgTransposerX)

end subroutine slice_computeSpectrum
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
!! theta= sorted Rayleigh quotients of size (neigenpairs)
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
        consdiff = (/ (theta(iband + 1) - theta(iband), iband=1,neigenpairs-1) /)
        
        ! Sort consecutive differences (=gaps) by increasing order
        call sort_dp(neigenpairs-1, consdiff, jperm, tol12)

        ! Take median of largest gaps
        do islice=1,nslice
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
            ndeg = slice%ndeg_filter
            !do while(1.d0/cheb_poly(part_low,ndeg,poly_upp,slice%maxeig_global)<ramp) 
            !    ndeg = ndeg + 1
            !end do
            if (1.d0/cheb_poly(part_low,ndeg,poly_upp,slice%maxeig_global)<ramp) then
                ! FIXME understand why this is always true?
                write(std_out,*) 'Warn: Chebyshev polynomial degree does not amplify enough'
                write(std_out,*) 1.d0/cheb_poly(part_low,ndeg,poly_upp,slice%maxeig_global)
            end if
        else
            ! ********* optimize amplification ratio ********
            ! The convergence ratio r0/rN is approximated by amplification ratios 
            ! f(l)/f(l-w) and f(u)/f(u+w). 
            lw = (poly_low - center)/radius ! scaled point outside slice
            uw = (poly_upp - center)/radius ! scaled point outside slice
            l = (part_low - center)/radius ! scaled point inside slice
            u = (part_upp - center)/radius ! scaled point inside slice
            ndeg = 4
            f_l = 0.d0; f_u = 0.d0; f_lw = 1.d0; f_uw = 1.d0
            do while ( (f_l/f_lw < ramp) .and. (f_u/f_uw < ramp) .and. (ndeg < ndeg_max) )
                f_l  = bandpassIndicator_sca(l ,lw,uw,ndeg)
                f_lw = bandpassIndicator_sca(lw,lw,uw,ndeg)
                f_u  = bandpassIndicator_sca(u ,lw,uw,ndeg)
                f_uw = bandpassIndicator_sca(uw,lw,uw,ndeg)
                ndeg = ndeg + 1
            end do
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
        write(std_out,'(a,i2)') '======= Slice ', islice
        write(std_out,*) '   Partition, width=', part_low, part_upp, part_upp - part_low
        write(std_out,*) 'With overlap, width=', poly_low, poly_upp, poly_upp - poly_low
        write(std_out,'(a,i6,a,i6)') 'nvec= ', nvec, ' ndeg= ', ndeg
        write(std_out,*) ' '

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
    integer, allocatable, target :: weights(:)
    integer, pointer :: weights_ptr(:) => null()
    integer, pointer :: task_sizes(:) => null()
    integer, pointer :: task_nprocs(:) => null()
    integer, pointer :: assigned_task(:) => null()
    
    ! *********************************************************************

    ABI_MALLOC_IFNOT(weights, (slice%nslice))

    weights_ptr => weights
    task_sizes => slice%neigenpairs_per_slice
    task_nprocs => slice%nproc_per_slice
    assigned_task => slice%lookup_proc

    ! Apply weighted fair allocation with a load balance criterion
    select case(slice%paral_slice)
    case(FAIR_BANDPP)
        weights = 1 
    case(FAIR_BANDPP_WDEG)
        weights = slice%poly_degrees
    end select
   
    ! Solve allocation problem to find the amount of resource allocated to each slice
    call fair_allocation(slice%nslice, task_sizes, weights_ptr, slice%nproc, task_nprocs)
    
    ! Call the subroutine to assign tasks (=slices) to processes
    call assign_tasks_to_processes(task_nprocs, assigned_task)
    do iproc = 1, slice%nproc
        write(std_out,'(a,i5,a,i5)') "Process ", iproc-1, " allocated to task ", slice%lookup_proc(iproc)
    end do

    call xmpi_barrier(slice%spacecom)

    ! Free temporary memory
    ABI_SFREE(weights) 

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
!! slice%me_comm_slice          = sub-communicator of active slice
!! slice%me_ncolsColsRows_slice = column capacity per process in colsrows repr
!! slice%me_nrowsLinalg_slice   = row capacity per process in linalg repr
!! slice%me_id_slice
!! slice%me_nproc_slice
!! slice%me_ndeg_slice
!! 
!! SOURCE

subroutine slice_markActiveTask(slice)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice

    ! Local variables
    integer :: my_rank, my_rank_sub, ierr
    integer, pointer :: me_colsrows_ptr(:) => null()
    integer, pointer :: me_linalg_ptr(:) => null()
    integer, pointer :: all_colsrows_ptr(:) => null()

    ! *********************************************************************
 
    my_rank = xmpi_comm_rank(slice%spacecom) 
    slice%me_id_slice = slice%lookup_proc(my_rank + 1) + 1
    slice%me_neigenpairs_slice = slice%neigenpairs_per_slice(slice%me_id_slice)
    slice%me_nproc_slice = slice%nproc_per_slice(slice%me_id_slice)
    slice%me_ndeg_slice = slice%poly_degrees(slice%me_id_slice)

    ABI_MALLOC_IFNOT(slice%me_ncolsColsRows_slice, (slice%me_nproc_slice))
    ABI_MALLOC_IFNOT(slice%me_nrowsLinalg_slice, (slice%me_nproc_slice))

    all_colsrows_ptr => slice%ncolsColsRows
    me_colsrows_ptr => slice%me_ncolsColsRows_slice
    me_linalg_ptr => slice%me_nrowsLinalg_slice

    ! Compute column distribution across active resources
    call distribute_vectors(slice%me_neigenpairs_slice, slice%me_nproc_slice, me_colsrows_ptr)
    
    ! Compute row distribution across active resources
    call distribute_vectors(slice%total_spacedim, slice%me_nproc_slice, me_linalg_ptr)

    ! Split global comm into disjoint sub-comms, only procs with the same color communicate
    if (slice%nslice==1) then
        slice%me_comm_slice = slice%spacecom
    else
        call xmpi_comm_split(slice%spacecom, slice%me_id_slice, my_rank, slice%me_comm_slice, ierr)        
        if ( ierr /= xmpi_success ) then
            ABI_ERROR("Error while creating slice subcommunicators")
        end if
    end if
    
    ! process waits for others to create their subcommunicators before using its own
    call xmpi_barrier(slice%spacecom) 

    ! Concatenate slice%me_ncolsColsRows_slice into collective slice%ncolsColsRows
    my_rank_sub = xmpi_comm_rank(slice%me_comm_slice)
    slice%me_bandpp_slice = slice%me_ncolsColsRows_slice(my_rank_sub + 1)
    call xmpi_allgather(slice%me_bandpp_slice, all_colsrows_ptr, slice%spacecom, ierr)
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
    integer :: my_rank, my_slice, neigenpairs_slice
    integer :: fcol_ext, fcol, tot_ncols_kept, lcol_ext
    integer :: islice, fcol_in_slice, lcol_in_slice
    real(dp) :: part_low_bound, part_upp_bound
    logical :: on_host, on_device
    ! Derived types
    type(xg_t) :: eigen_ext
    type(xg_t) :: resid_ext
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

    call xgBlock_zero(eigen_ext%self)
    call xgBlock_zero(resid_ext%self)

    ! Bring to columns to be able to select column range
    call xgBlock_reshape(eigen, 1, slice%neigenpairs)
    call xgBlock_reshape(resid, 1, slice%neigenpairs)

    ! MPI communication to gather slice eigen/resid to eigen_ext/resid_ext
    if (slice%paral_kgb==1) then
        if (xmpi_comm_size(slice%spacecom) > 1) then
        
            my_rank = xmpi_comm_rank(slice%spacecom)
            my_slice = slice%lookup_proc(my_rank + 1)
            neigenpairs_slice = slice%neigenpairs_per_slice(my_slice + 1)
            fcol_ext = slice%fcol_in_Xext(my_slice + 1)

            call xgBlock_setBlock(eigen_ext%self, eigen_ext_slice, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
            call xgBlock_setBlock(resid_ext%self, resid_ext_slice, rows=1, cols=neigenpairs_slice, fcol=fcol_ext)
            call xgBlock_copy(eigen, eigen_ext_slice)
            call xgBlock_copy(resid, resid_ext_slice)

            ! All processes wait before summing 
            call xmpi_barrier(slice%spacecom)

            call xgBlock_mpi_sum(eigen_ext%self, comm=slice%spacecom)
            call xgBlock_mpi_sum(resid_ext%self, comm=slice%spacecom)

        end if
    else
        call xgBlock_copy(eigen, eigen_ext%self)
        call xgBlock_copy(resid, resid_ext%self)
    end if 

    if (slice%gpu_option==1) then
        call xgBlock_copy_from_gpu(eigen_ext%self)
        call xgBlock_copy_from_gpu(resid_ext%self)
    end if

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
        theta_reshaped(1:neigenpairs_slice) = theta_ext(1, fcol_ext:lcol_ext)

        ! Apply filter criterion to find kept first and last column in slice
        part_low_bound = slice%part_low_bounds(islice)
        part_upp_bound = slice%part_upp_bounds(islice)
        fcol_in_slice = maxloc(theta_reshaped, dim=1, mask=(theta_reshaped < part_low_bound)) + 1
        lcol_in_slice = maxloc(theta_reshaped, dim=1, mask=(theta_reshaped < part_upp_bound))            
        if (islice == 1     ) fcol_in_slice = 1
        if (islice == slice%nslice) lcol_in_slice = neigenpairs_slice

        !write(std_out,*) 'Filter in ', part_low_bound, part_upp_bound
        !write(std_out,*) 'kept indices', fcol_in_slice, lcol_in_slice
        !write(std_out,*) theta_reshaped

        ! After merge: Update first columns to copy from Xext to X
        slice%fcol_in_X(islice)= tot_ncols_kept + 1
        slice%fcol_in_Xext(islice) = fcol_ext + fcol_in_slice - 1
        slice%neigenpairs_per_slice(islice) = lcol_in_slice - fcol_in_slice + 1
        tot_ncols_kept = tot_ncols_kept + slice%neigenpairs_per_slice(islice) 

        ABI_SFREE(theta_reshaped)

    end do
    
    ! Detect missing or extra eigenvalues
    if (tot_ncols_kept < slice%neigenpairs) then
        ABI_ERROR("Too few converged eigenvalues kept")
    else if (tot_ncols_kept > slice%neigenpairs) then
        ABI_ERROR("Too many converged eigenvalues kept")
    end if

    ! Copy from extended memory to regular memory
    do islice=1,slice%nslice
        fcol = slice%fcol_in_X(islice)
        fcol_ext = slice%fcol_in_Xext(islice)
        neigenpairs_slice = slice%neigenpairs_per_slice(islice)
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

    do i=1,n
        write(std_out,*) '# proc # workload # allocated resources', i, real(w(i)*m(i))/real(x(i)), x(i)
    end do

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

!!****f* m_slice/bandpassIndicator_sca
!! NAME
!! bandpassIndicator_sca
!!
!! FUNCTION
!! Scalar Chebyshev-Jackson polynomial filter f(x) approximating an 
!! indicator function, using degree deg evaluated at point x=t
!!
!! INPUTS
!!  a,b=    interval to amplify included in -1,1
!!  t=      scalar to evaluate filter on
!!  deg=    order of Chebyshev expansion
!!
!! OUTPUT
!!  res
!!
!! SOURCE

function bandpassIndicator_sca(t,a,b,deg) result(f_t)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in ) :: t,a,b
    integer , intent(in ) :: deg

    real(dp) :: f_t
    
    !Local variables-------------------------------
    real(dp) :: yt0,yt,yt_swap,ck,mu,damp
    integer  :: i
    
    ! *********************************************************************

    ! init cheby of deg=0,1 eval at t
    yt0 = 1.d0
    yt = t

    ! init filter for deg=0
    ck = Pi/(deg+2)
    mu = 1/Pi*(ACOS(a)-ACOS(b))
    damp = 1.d0
    f_t = mu * damp * yt0

    do i=1,deg 
        
        ! Update damping and expansion coefficient
        mu = 2/Pi * (SIN(i*ACOS(a)) - SIN(i*ACOS(b)))/i
        damp = ((1 - i/(deg+2))*SIN(ck)*COS(i*ck) + 1/(deg+2)*COS(ck)*SIN(i*ck))/SIN(ck)

        ! Sum terms
        f_t = f_t + mu * damp * yt

        ! Update Chebyshev polynomial
        yt_swap = yt
        yt = 2 * t * yt - yt0
        yt0 = yt_swap
        
    end do

end function bandpassIndicator_sca
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

end module m_slice
!!***

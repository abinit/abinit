!!****f* ABINIT/m_slice
!! NAME
!! m_slice
!!
!! FUNCTION
!! This module contains the types and routines used to apply 
!! the Spectrum Slicing method. It mainly defines 'slice' 
!! datatypes and associated methods. Features:
!! 
!! - uses 'xgTools' implementation for matrix data structure.
!! - adopts 'chebfi' functionalities for most matrix calculations.
!! - uses polynomial filtering, Chebyshev for first slice and
!!   Chebyshev-Jackson expansion of Heaviside otherwise.
!! - implements scheduler and resource allocator for slice tasks.
!! - applies Rayleigh-Ritz for individual slices in parallel.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (IL, LB)
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
    use m_polyfi

    use m_xmpi
    use m_xomp
#ifdef HAVE_OPENMP
    use omp_lib
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize
#endif

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
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
    integer, parameter :: FAIR_BANDPP      = 1 ! bandpp (unweigted)
    integer, parameter :: FAIR_BANDPP_WDEG = 2 ! bandpp weighted by degree 

    ! Public 'slice' datatype
    !-------------------------------------------------
    type, public :: slice_t

        ! MPI-related (me=current MPI process)
        integer :: nproc                ! number of available processes
        integer :: comm_rows            ! xmpi_comm_self ...
        integer :: comm_cols            ! same as spacecom
        integer :: spacecom             ! same as comm_cols
        integer :: me_g0                ! process contains G(0,0,0) 
        integer :: me_g0_fft            ! process contains G(0,0,0) for fft
        integer :: me_nproc_slice       ! number of processes reserved to slice task in use
        integer :: me_comm_slice        ! sub-communicator reserved to slice task in use
        integer :: me_id_slice          ! identifier of slice task in use
        integer :: me_neigenpairs_slice ! total number of eigenpairs of slice task in use
        integer :: me_bandpp_slice      ! number of distributed bands per process for task in use

        ! Eigenpair parameters
        integer :: neigenpairs       ! total number of bands (=number of eigenpairs)
        integer :: neigenpairs_ext   ! total numner of extended columns
        integer :: total_spacedim    ! total number of plane-waves
        integer :: bandpp            ! nb of bands per process in colsrows representation 
        integer :: spacedim          ! nb of plane-waves per process in linalg representation
        integer :: nslice            ! number of spectral slices
        integer :: space             ! real or complex eigenvectors
        integer :: space_res         ! real or complex eigenvalues
        
        ! GPU-related
        integer :: gpu_kokkos_nthrd
        integer :: gpu_thread_limit

        ! Flags
        logical :: on_host = .false.        ! running on CPU
        logical :: on_device = .false.      ! running on GPU
        logical :: use_linalg = .false.     ! use linalg representation
        logical :: use_colsrows = .false.   ! use colsrows representation

        ! Options
        integer :: gpu_option    ! enable GPU
        integer :: paral_kgb     ! enable parallel (k-points, G basis, bands)
        integer :: paral_slice   ! how to allocate resources for parallel slices
        integer :: spectral_cut  ! how to decompose spectrum into slices

        ! Various model parameters
        logical :: paw              ! use PAW or not 
        integer :: ndeg_filter      ! lowpass degree of polynomial filter
        real(dp) :: ramp            ! bandpass filter tolerance
        real(dp) :: ecut            ! Ecut Fermi level
        real(dp) :: mineig_global   ! guaranteed lower spectral bound
        real(dp) :: maxeig_global   ! guaranteed upper spectral bound

        ! Memory buffers
        type(xg_t) :: DivResults
        type(xg_t) :: X_ext         ! eigenvector memory used by all slices
        type(xgTransposer_t) :: xgTransposerXext
        
        ! Pointers
        type(xgBlock_t) :: me_Xext_active   ! eigenvector memory in use by active slice
        type(xgBlock_t) :: XextLinalg

        ! Arrays for my slice only
        integer, allocatable :: me_ncolsColsRows_slice(:)   ! ncol of colsrows representation for my slice
        integer, allocatable :: me_nrowsLinalg_slice(:)     ! nrow of linalg representation for my slice

        ! Arrays related to MPI (all slices)
        integer, allocatable :: neigenpairs_per_slice(:)   ! number of total eigenpairs per slice
        integer, allocatable :: nproc_per_slice(:)         ! number of processes per slice
        integer, allocatable :: lookup_proc(:)             ! which slice each process serves
        integer, allocatable :: ncolsColsRows(:)           ! ncol of colsrows representation for global X

        ! Arrays related to data distribution (all slices)
        integer, allocatable :: fcol_in_X(:)           ! first band of slice in spectrum memory
        integer, allocatable :: fcol_in_Xext(:)        ! first band of slice in extended memory

        ! Arrays related to polynomial filtering (all slices)
        integer, allocatable :: poly_degrees()         ! polynomial filter degrees
        real(dp), allocatable :: part_low_bounds(:)    ! lower bounds in spectral partition (disjoint)
        real(dp), allocatable :: part_upp_bounds(:)    ! upper bounds in spectral partition (disjoint) 
        real(dp), allocatable :: poly_low_bounds(:)    ! lower bounds used to define polynomials (overlap)
        real(dp), allocatable :: poly_upp_bounds(:)    ! upper bounds used to define polynomials (overlap)

    end type slice_t

    ! Public methods associated to 'slice' datatype
    !-------------------------------------------------
    public :: slice_init                ! initialize slice datatype object
    public :: slice_free                ! free slice datatype object
    public :: slice_schedule            ! build slice distributed workspace
    public :: slice_run                 ! run Spectrum Slicing for active slice task
    public :: slice_merge               ! merge slice result to distributed workspace
    public :: slice_unitTest            ! used for debugging GPU device

    CONTAINS  
!=====================================================================
!!***

!!****f* m_slice/slice_init
!! NAME
!! slice_init
!!
!! FUNCTION
!! Initialize a 'slice' datastructure. Memory needed to compute
!! Rayleigh quotients for neigenpairs and create nslice slices.
!!
!! SOURCE

subroutine slice_init(slice,nslice,neigenpairs,spacedim,paral_kgb,paral_slice,ndeg_filter,&
        ramp,ecut,bandpp,space,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
        spectral_cut,gpu_option,gpu_kokkos_nthrd,gpu_thread_limit)

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
    real(dp)     , intent(in   ) :: ecut
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ), optional :: gpu_kokkos_nthrd
    integer      , intent(in   ), optional :: gpu_thread_limit
    
    ! Local variables --------------------------------
    integer :: ierr

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
    slice%ecut          = ecut

    slice%gpu_kokkos_nthrd = 1
    if (present(gpu_kokkos_nthrd)) slice%gpu_kokkos_nthrd = gpu_kokkos_nthrd
    slice%gpu_thread_limit = 0
    if (present(gpu_thread_limit)) slice%gpu_thread_limit = gpu_thread_limit

    ! Space of eigenvalues
    if (space==SPACE_C) then
        slice%space_res = SPACE_C
    else if (space==SPACE_CR) then
        slice%space_res = SPACE_R
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

    if (slice%paral_kgb == 0) then
        call xg_init(slice%DivResults, slice%space_res, rows=slice%neigenpairs, cols=1, &
            gpu_option=slice%%gpu_option)
    else
        call xg_init(slice%DivResults, slice%space_res, rows=slice%bandpp, cols=1, &
            gpu_option=slice%gpu_option)
    end if

    if(.not.allocated(slice%neigenpairs_per_slice)) ABI_MALLOC(slice%neigenpairs_per_slice, (slice%nslice))
    if(.not.allocated(slice%nproc_per_slice)) ABI_MALLOC(slice%nproc_per_slice, (slice%nslice))
    if(.not.allocated(slice%lookup_proc)) ABI_MALLOC(slice%lookup_proc, (slice%nproc))

    if(.not.allocated(slice%fcol_in_X)) ABI_MALLOC(slice%fcol_in_X, (slice%nslice))
    if(.not.allocated(slice%fcol_in_Xext)) ABI_MALLOC(slice%fcol_in_Xext, (slice%nslice))
    if (.not.allocated(slice%ncolsColsRows)) ABI_MALLOC(slice%ncolsColsRows, (slice%nproc))

    if(.not.allocated(slice%poly_degrees)) ABI_MALLOC(slice%poly_degrees, (slice%nslice))
    if(.not.allocated(slice%part_low_bounds)) ABI_MALLOC(slice%part_low_bounds, (slice%nslice))
    if(.not.allocated(slice%part_upp_bounds)) ABI_MALLOC(slice%part_upp_bounds, (slice%nslice))
    if(.not.allocated(slice%poly_low_bounds)) ABI_MALLOC(slice%poly_low_bounds, (slice%nslice))
    if(.not.allocated(slice%poly_upp_bounds)) ABI_MALLOC(slice%poly_upp_bounds, (slice%nslice))

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

    call xg_free(slice%DivResults)
    call xg_free(slice%X_ext)
    call xgTransposer_free(slice%xgTransposerX)

    if(allocated(slice%me_ncolsColsRows_slice)) ABI_FREE(slice%me_ncolsColsRows_slice)   
    if(allocated(slice%me_nrowsLinalg_slice)) ABI_FREE(slice%me_nrowsLinalg_slice)

    if(allocated(slice%neigenpairs_per_slice)) ABI_FREE(slice%neigenpairs_per_slice)
    if(allocated(slice%nproc_per_slice)) ABI_FREE(slice%nproc_per_slice)
    if(allocated(slice%lookup_proc)) ABI_FREE(slice%lookup_proc)
    if(allocated(slice%ncolsColsRows)) ABI_FREE(slice%ncolsColsRows)   

    if(allocated(slice%fcol_in_X)) ABI_FREE(slice%fcol_in_X)
    if(allocated(slice%fcol_in_Xext)) ABI_FREE(slice%fcol_in_Xext)

    if(allocated(slice%poly_degrees)) ABI_FREE(slice%poly_degrees)
    if(allocated(slice%part_low_bounds)) ABI_FREE(slice%part_low_bounds)
    if(allocated(slice%part_upp_bounds)) ABI_FREE(slice%part_upp_bounds)
    if(allocated(slice%poly_low_bounds)) ABI_FREE(slice%poly_low_bounds)
    if(allocated(slice%poly_upp_bounds)) ABI_FREE(slice%poly_upp_bounds)

end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_schedule
!! NAME
!! slice_schedule
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
!! SIDE EFFECTS
!! slice%me_Xext_active  = adapted version distributed correctly
!!        and free of data overlap. Suitable for parallel RR.
!! 
!! SOURCE

subroutine slice_schedule(slice, X0, getAX_BX, nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X0
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
    integer :: neigenpairs, iband, min_loc
    integer :: ierr
    logical :: on_host, on_device
    type(xg_t) :: eigen0
    type(xg_t) :: resid0
    real(dp) :: lambda_minus, lambda_plus
    real(dp) :: tol12 = 1.0e-12
    ! Derived types
    type(xgBlock_t) :: slicecols_in
    type(xgBlock_t) :: slicecols_ext_out
    ! Arrays
    integer, allocatable, target :: permute_cols(:)
    integer, pointer :: permute_cols(:) => null()
    real(dp), pointer :: theta_(:,:)
    real(dp), pointer :: resid_(:,:)
    real(dp), allocatable, target :: theta_reshaped(:)
    real(dp), pointer, theta_reshaped_ptr => null()
    
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
    
    ! Space for computed eigenvalues and residuals (not distributed)
    call xg_init(eigen0, slice%space_res, rows=1, cols=neigenpairs, gpu_option=slice%gpu_option)
    call xg_init(resid0, SPACE_R        , rows=1, cols=neigenpairs, gpu_option=slice%gpu_option)
    call xgBlock_zero(eigen0%self)
    call xgBlock_zero(resid0%self)

    if(.not.allocated(theta_reshaped)) ABI_MALLOC(theta_reshaped, (neigenpairs))
    if(.not.allocated(permute_cols)) ABI_MALLOC(permute_cols, (neigenpairs))
    
    ! ===================== Compute Rayleigh quotients and residuals =====================
    
    ABI_NVTX_START_RANGE(NVTX_SLICE_RRQ)
    call slice_computeSpectrum(slice, X0, getAX_BX, eigen0%self, resid0%self, nspinor)
    ABI_NVTX_END_RANGE()

    ! ===================== Compute guaranteed spectral bounds ===================== 

    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call xgBlock_copy_from_gpu(eigen0%self)
        call xgBlock_copy_from_gpu(resid0%self)
    end if

    ! Results could be complex, so neigenpairs has to be in cols, not rows
    call xgBlock_reshape(eigen0%self, neigenpairs, 1)     
    call xgBlock_reshape(resid0%self, neigenpairs, 1)
    call xgBlock_reverseMap(eigen0%self, theta_, rows=1, cols=neigenpairs)
    call xgBlock_reverseMap(resid0%self, resid_, rows=1, cols=neigenpairs)
 
    ! Sort thetas in increasing order and store permutation
    theta_reshaped(1:neigenpairs) = theta_(1,1:neigenpairs)
    permute_cols_ptr => permute_cols
    permute_cols(1:neigenpairs) = (/ (iband, iband=1,neigenpairs) /)
    call sort_dp(neigenpairs, theta_reshaped, permute_cols_ptr, tol12)

    ! Minimum and maximum quotients
    lambda_minus = theta_reshaped(1)
    lambda_plus = theta_reshaped(neigenpairs)

    ! Guaranteed spectral bounds
    min_loc = permute_cols(1)
    slice%mineig_global = lambda_minus - sqrt(resid_(1, min_loc))
    slice%maxeig_global = slice%ecut

    ! ===================== Decompose interval [lambda_minus,lambda_plus) to slices =====================
    
    theta_reshaped_ptr => theta_reshaped
    call slice_cutSpectrum(slice, lambda_minus, lambda_plus, theta_reshaped_ptr, plot_filter=.false.)

    ! ======== Resource management system ====================
    
    ! Divide resources into slice tasks
    call slice_allocateResources(slice)

    ! Run on all ranks of spacecom: Mark my slice resources as actively in use
    call slice_markActiveResources(slice)

    ! ======== Allocate and distribute extended memory buffer ====================
  
    ! Sanity check
    if ((.not. slice%use_linalg) .or. slice%use_colsrows) then
        ABI_ERROR("not in linalg representation")
    end if

    ! Permute column vectors in Rayleigh quotient increasing order (on CPU)
    if (slice%gpu_option==ABI_GPU_OPENMP) then
        call xgBlock_copy_from_gpu(X0)
    end if
    call xgBlock_permuteCols(X0, slice%total_spacedim, neigenpairs, permute_cols_ptr)


    slice%neigenpairs_ext = sum(slice%neigenpairs_per_slice)
    
    if (slice%nslice==1) then
        slice%XextLinalg = X0
    else
 
        ! Allocate extended space in linalg representation (on CPU)
        call xg_init(slice%X_ext, slice%space, slice%spacedim, slice%neigenpairs_ext, &
            slice%spacecom, me_g0=slice%me_g0, gpu_option=ABI_GPU_DISABLED)
        
        slice%XextLinalg = slice%X_ext%self

        ! Copy X to XextLinalg by column blocks
        do islice=1, nslice
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
    ABI_CHECK(cols(slice%XextLinalg)==slice%neigenpairs_ext,'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # cols of Xext_linalg ', xmpi_comm_rank(slice%spacecom), cols(slice%XextLinalg)

    ! Distribute extended columns across **all** MPI processes
    ! After the transposition each process contains the correct
    ! bandpp corresponding to the slice so that no additional communication
    ! has to be performed in order to bring band slices to processes.
    ncolsColsRows_ptr => slice%ncolsColsRows

    ! Allocate slice%me_Xext_active according to the target MPI distribution for slices
    call xgTransposer_constructor(slice%xgTransposerX, slice%XextLinalg, slice%me_Xext_active,&
        nspinor, STATE_LINALG, TRANS_ALL2ALL, slice%comm_rows, slice%comm_cols, 0, 0, slice%me_g0_fft,&
        gpu_option=slice%gpu_option, gpu_thread_limit=slice%gpu_thread_limit,&
        custom_ncolsColsRows=.true., ncolsColsRows_sub=ncolsColsRows_ptr)
   
    slice%xgTransposerX%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
   
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(slice%xgTransposerX, STATE_COLSROWS)
    ABI_NVTX_END_RANGE()

    ! Unitary test
    if ( cols(slice%me_Xext_active)==slice%ncolsColsRows(xmpi_comm_rank(slice%spacecom)) ) then
        ABI_ERROR('wrong colsrows representation')
    end if
    write(*,'(a,i6,i6)') '# proc has # cols of Xext ', xmpi_comm_rank(slice%spacecom), cols(slice%me_Xext_active)

    ! Free memory
    call xg_free(eigen0)
    call xg_free(resid0) 
    if (allocated(theta)) ABI_FREE(theta)
    if (allocated(resid)) ABI_FREE(resid)

    ABI_NVTX_END_RANGE()

end subroutine slice_schedule
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
!! slice%me_Xext_active  =eigenvector guess in colsrows representation 
!!
!! OUTPUT
!! eigen                 =converged eigenvalue array of size neigenpairs
!! residu                =slice residual array of size neigenpairs
!! 
!! SIDE EFFECTS
!! slice%me_Xext_active  =converged eigenvectors on slice
!! 
!! SOURCE

subroutine slice_run(slice, getAX_BX, getBm1X, eigen, residu, nspinor)

    implicit none

    ! Arguments
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X
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

    ! Variables
    type(polyfi_t) :: polyfi

    ! *********************************************************************
 
    call slice_getActiveTask(slice,nband_sub,spacecom_sub,mineig_global,maxeig_global,&
         lambda_minus,lambda_plus,nrowsLinalg_ptr)

    call polyfi_init(polyfi,nband_sub,dtset%tolwfr_diago,dtset%ecut,&
        dtset%paral_kgb,space,1,spacecom_sub,&
        me_g0,me_g0_fft,l_paw,l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
        mineig_global,maxeig_global,lambda_minus,lambda_plus,nrowsLinalg,
        l_gs_hamk%gpu_option,gpu_kokkos_nthrd=dtset%gpu_kokkos_nthrd,&
        gpu_thread_limit=dtset%gpu_thread_limit)

    call polyfi_run(polyfi,slice%Xext,getghc_gsc1,getBm1X,xgeigenslice,xgresiduslice,nspinor)

    call polyfi_free(polyfi)

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
!! eigen  =Rayleigh Quotients from X
!! residu =residuals from X
!! 
!! SIDE EFFECTS
!! slice%DivResults= contains eigen_mpi or eigen depending on paral_kgb
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
    integer :: my_rank, shift, bandpp
    real(dp) :: mineig, maxeig
    ! Derived types
    type(xg_t) :: Results1
    type(xg_t) :: Results2
    type(xg_t) :: Results3
    type(xg_t) :: eigen_mpi
    type(xg_t) :: resid_mpi
    type(xg_t) :: X_NAB ! vector memory for X_next, AX, BX
    type(xgBlock_t) :: xAXColsRows
    type(xgBlock_t) :: xBXColsRows
    type(xgBlock_t) :: X_next
    type(xgTransposer_t) :: xgTransposerX
    type(xgBlock_t) :: xXColsRows
    ! Arrays
    integer :: maxeig_pos(2)
    integer :: mineig_pos(2)

    ! *********************************************************************

    if (slice%paral_kgb==0) then
        ABI_ERROR("Sequential bands not implemented")
    end if

    total_spacedim = slice%total_spacedim
    bandpp = slice%bandpp

    ! Allocate temporary memory (distributed in colsrows representation)
    call xg_init(X_NAB, slice%space, total_spacedim, 3*bandpp, slice%spacecom, &
        me_g0=slice%me_g0, gpu_option=slice%gpu_option)

    call xg_setBlock(X_NAB, X_next, total_spacedim, bandpp)                               ! X_next
    call xg_setBlock(X_NAB, xAXColsRows, total_spacedim, bandpp, fcol=bandpp + 1)         ! xAXColsRows
    call xg_setBlock(X_NAB, xBXColsRows, total_spacedim, bandpp, fcol=2*bandpp + 1)       ! xBXColsRows
    
    ! Allocate one-dimensional memory (distributed)
    ! for eigenvalues
    call xg_init(eigen_mpi, slice%space_res, rows=bandpp, cols=1, comm=slice%spacecom, gpu_option=slice%gpu_option)
    call xg_init(Results1, slice%space_res, rows=bandpp, cols=1, gpu_option=slice%gpu_option)
    call xg_init(Results2, slice%space_res, rows=bandpp, cols=1, gpu_option=slice%gpu_option)
    ! for residuals
    call xg_init(resid_mpi, SPACE_R, rows=bandpp, cols=1, comm=slice%spacecom, gpu_option=slice%gpu_option)
    call xg_init(Results3, SPACE_R, rows=bandpp, cols=1, gpu_option=slice%gpu_option)
 
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
    call xgBlock_colwiseNorm2(X_next, resid, comm_loc=xmpi_comm_null)  ! resid = |X_next|^2
   
    ! MPI communication for gathering all thetas (using summation strategy on columns)
    call xgBlock_reshape(eigen_mpi%self, 1, bandpp) 
    call xgBlock_reshape(resid_mpi%self, 1, bandpp)     
    if (xmpi_comm_size(slice%spacecom) > 1) then
        
        my_rank = xmpi_comm_rank(slice%spacecom)
        shift = my_rank * bandpp

        call xgBlock_setBlock(eigen, Results1%self, nrows=1, ncols=bandpp, fcol=1+shift)
        call xgBlock_copy(eigen_mpi%self, Results1%self)
        call xgBlock_mpi_sum(eigen, comm=slice%spacecom)

        call xgBlock_setBlock(resid, Results3%self, nrows=1, ncols=bandpp, fcol=1+shift)
        call xgBlock_copy(resid_mpi%self, Results3%self)
        call xgBlock_mpi_sum(resid, comm=slice%spacecom)
    else
        call xgBlock_copy(eigen_mpi%self, eigen)
        call xgBlock_copy(resid_mpi%self, resid)
    end if

    ! Save the eigenvalue result to slice%DivResults depending on the parallelization
    if (slice%paral_kgb==0) then
        call xgBlock_reshape(eigen_mpi%self, bandpp, 1) 
        call xgBlock_copy(eigen_mpi%self, slice%DivResults%self)
    else if (slice%paral_kgb==1) then
        call xgBlock_copy(eigen, slice%DivResults%self)
    end if

    ! ============== Transpose ==============
    call xmpi_barrier(slice%spacecom)
    ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
    call xgTransposer_transpose(xgTransposerX, STATE_LINALG)
    ABI_NVTX_END_RANGE()

    slice%use_linalg = .true.
    slice%use_colsrows = .false.

    ! Free memory
    call xg_free(Results1)
    call xg_free(Results2)
    call xg_free(Results3)
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
!! (Optional) Print filter as {x,f(x)} in std_out 
!! 
!! INPUTS
!! theta =sorted Rayleigh quotients
!! 
!! SOURCE

subroutine slice_cutSpectrum(slice, lambda_minus, lambda_plus, theta, plot_filter)

    implicit none

    !Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    real(dp), pointer, intent(in) :: theta
    logical, optional, intent(in) :: plot_filter
    
    !Local variables-------------------------------
    integer :: j,k,jmax,spos,nvec,nvec,k1,k2
    integer :: nv_pad,ndeg,npband,comm_cols
    integer :: nslice,neigenpairs,ndeg_filter, n_frac
    integer :: ndeg_max = 200
    integer :: ipt,npt,iptL,iptR
    logical :: plot_filter_
    real(dp) :: ramp, width
    real(dp) :: tol12 = 1.0e-12
    real(dp) :: ecut,low,upp,glb,gub,c,r
    real(dp) :: lj,uj,wj,finL,finR,foutL,foutR
    real(dp) :: f_l, f_lw, f_u, f_uw
    real(dp) :: wlj,wuj
    real(dp) :: fun_pt,pt
    real(dp) :: a_,b_ ! target interval scaled in -1,1
    ! arrays
    integer :: jperm(nband-1)
    real(dp) :: consdiff(nband-1)
    real(dp), allocatable :: spectral_partition(:)

    ! *********************************************************************

    ramp = slice%ramp
    nslice = slice%nslice
    ndeg_filter = slice%ndeg_filter
    plot_filter_ = .false.
    if (present(plot_filter)) plot_filter_ = plot_filter

    if(.not.allocated(spectral_partition)) ABI_MALLOC(spectral_partition,(nslice+1)) 
   
    ! Compute spectral partition
    spectral_partition(:) = 0.d0
    spectral_partition(1) = theta(1)                  ! lambda_minus
    spectral_partition(nslice+1) = theta(neigenpairs) ! lambda_plus

    ! Je suis bête on a besoin de theta ici
    select case(slice%spectral_cut)
    case(DIVIDE_INTERVAL_WIDTH)

        width = (lambda_plus - lambda_minus) / nslice
        spectral_partition(2:nslice) = (/ (lambda_minus + width*islice, islice=1,nslice-1) /)

    case(DIVIDE_NUMBER_OF_VECTORS)

        n_frac = neigenpairs / nslice
        spectral_partition(2:nslice) = (/ (theta(n_frac * islice), islice=1,nslice-1) /)

    case(DIVIDE_SPECTRAL_GAPS)

        jperm = (/ (iband, iband=1,neigenpairs-1) /)
        consdiff = (/ (theta(iband + 1) - theta(iband), iband=1,neigenpairs-1) /)
        
        ! Sort consecutive differences (=gaps) by increasing order
        call sort_dp(neigenpairs-1, consdiff, jperm, tol12)

        ! Take median of largest gaps
        do islice=1,nslice
            jmax = jperm(neigenpairs -i)
            spectral_partition(islice + 1) = (theta(jmax) + theta(jmax + 1)) / 2.d0
        end do

    end select

    ! Center and radius of the entire spectrum mapped to -1,1
    center = (slice%mineig_global + slice%maxeig_global) / 2.d0
    radius = (slice%maxeig_global - slice%mineig_global) / 2.d0

    ! Define spectral subintervals and optimize degrees for individual slices
    first_col_ext = 1
    do islice=1, nslice 

        part_low = spectral_partition(islice)
        part_upp = spectral_partition(islice+1)
        wovlp = (part_upp - part_low)/10.d0 ! FIXME allow tuning from abi param
        l = (part_low - center)/radius
        u = (part_upp - center)/radius

        ! TODO Priority IL 14/4
        ! Very much attention to this. It is messed up.
        ! chebfi_run normally takes lambda_minus = maxeig_global and lambda_plus = ecut.
        ! So lambda_minus should really be the largest wanted eigenvalue. We must
        ! thus use the extended interval including overlap in there (poly_upp).

        if (islice==1) then
            ! Spectral interval to amplify is [-oo, lambda_minus)
            poly_low = part_upp + wovlp
            poly_upp = slice%ecut
            ndeg_filter = slice%ndeg_filter
        else 
            poly_low = part_low - wovlp
            poly_upp = part_upp + wovlp
        end if
        lw = (poly_low - center)/radius
        uw = (poly_upp - center)/radius


        if (islice==1) then
            ! uj,gub is the interval mapped to -1,1
            ! in this interval Chebyshev poly is bounded by 1
            ! remember uj,gub is the interval to ignore
            ndeg = 4
            ! FIXME is the scaling in [-1,1] OK?
            do while(1.d0/cheb_poly(poly_low,ndeg,poly_upp,ecut)<ramp) 
                ndeg = ndeg + 1
            end do
        else
            ! ********* optimize amplification ratio ********
            ! The convergence ratio r0/rN is approximated by amplification ratios 
            ! f(l)/f(l-w) and f(u)/f(u+w). 
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

        ! Plot filter
        if (plot_filter_) then
            call print_scalar_filter(low,upp,center,radius,npt,is_lowpass=(islice==1))
        end if

        ! Print slice interval info
        write(std_out,*) 'Without overlap=', low, upp
        write(std_out,*) '          width=', upp-low
        write(std_out,*) '      scaled to=', (low-c)/r,(upp-c)/r
        write(std_out,*) 'With    overlap=', poly_low, poly_upp
        write(std_out,*) '          width=', poly_upp-poly_low
        write(std_out,*) '      scaled to=', a,b
        write(std_out,*) '           nvec=', nvec
        write(std_out,*) '           nvec=', nvec
        write(std_out,*) '           ndeg=', ndeg
        write(std_out,*) ' '
        ! end print

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

    if (allocated(spectral_partition)) ABI_FREE(spectral_partition)
   
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

    if(.not.allocated(weights)) ABI_MALLOC(weights, (slice%nslice))

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
    do iproc = 1, slice%nslice
        write(*,'(a,i5,a,i5)') "Process ", iproc, " is in task ", slice%lookup_proc(iproc)
    end do

    ! Free temporary memory
    if (allocated(weights)) ABI_FREE(weights) 

end subroutine slice_allocateResources
!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_markActiveResources
!! NAME
!! slice_markActiveResources
!!
!! FUNCTION
!! Mark allocated portion as actively used. 
!! My process only marks resources its assigned slice has reserved. 
!!
!! SIDE EFFECTS
!! slice%me_comm_slice          = sub-communicator of active slice
!! slice%me_ncolsColsRows_slice = column capacity per process in colsrows repr
!! slice%me_nrowsLinalg_slice   = row capacity per process in linalg repr
!! 
!! SOURCE

subroutine slice_markActiveResources(slice)

    implicit none

    ! Arguments
    type(slice_t), target, intent(inout) :: slice

    ! Local variables
    integer :: color, my_rank, my_rank_sub, ierr
    integer, pointer :: me_colsrows_ptr(:) => null()
    integer, pointer :: me_linalg_ptr(:) => null()
    integer, pointer :: all_colsrows_ptr(:) => null()

    ! *********************************************************************
 
    my_rank = xmpi_comm_rank(slice%spacecom) 
    slice%me_id_slice = slice%lookup_proc(my_rank + 1)
    slice%me_neigenpairs_slice = slice%neigenpairs_per_slice(slice%me_id_slice)
    slice%me_nproc_slice = slice%nproc_per_slice(slice%me_id_slice)

    if (.not.allocated(slice%me_ncolsColsRows_slice)) then
        ABI_MALLOC(slice%me_ncolsColsRows_slice, (slice%me_nproc_slice))
    end if
    if (.not.allocated(slice%me_nrowsLinalg_slice)) then
        ABI_MALLOC(slice%me_nrowsLinalg_slice, (slice%me_nproc_slice))
    end if

    all_colsrows_ptr => slice%ncolsColsRows
    me_colsrows_ptr => slice%me_ncolsColsRows_slice
    me_linalg_ptr => slice%_me_nrowsLinalg_slice

    ! Compute column distribution across active resources
    call distribute_vectors(slice%me_nproc_slice, slice%me_neigenpairs_slice, me_colsrows_ptr)
    
    ! Compute row distribution across active resources
    call distribute_vectors(slice%me_nproc_slice, slice%total_spacedim, me_linalg_ptr)

    ! Split global comm into disjoint sub-comms, only procs with the same color communicate
    call xmpi_comm_split(slice%spacecom, slice%me_id_slice, my_rank, slice%me_comm_slice, ierr)        
    if ( ierr /= xmpi_success ) then
        ABI_ERROR("Error while creating slice subcommunicators")
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

end subroutine slice_markActiveResources
!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_getActiveTask
!! NAME
!! slice_getActiveTask
!!
!! FUNCTION
!! Get slice parameters for active slice based on actived resources.
!!
!! SOURCE

subroutine slice_getActiveTask(slice,glb,gub,lb,ub)
    
    implicit none

    ! Arguments
    type(slice_t), target, intent(in) :: slice

    ! Variables

    ! *********************************************************************

    ! TODO should give comm_slice, mineig_global, maxeig_global,
    ! lambda_minus, lambda_plus, ndeg_filter
    ! also bandpp, nband
    ! all information used to define chebfi on spectral slice
    bandpp = slice%ncolsColsRows_slice(my_rank+1)
    nband = slice%neigenpairs_slice
    is_lowpass = (slice%lookup_proc(my_rank+1) == 1)

    glb = slice%mineig_global
    gub = slice%maxeig_global
    lambda_min = slice%poly_low_bounds(islice)
    ub = slice%poly_upp_bounds(islice)

end subroutine slice_getActiveTask
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_merge
!! NAME
!! slice_merge
!! 
!! FUNCTION
!! Filter converged eigenvalues of each slice using a criterion
!! based on spectral partition interval bounds.
!! Assumes linalg representation, so after all transpositions.
!!
!! INPUTS
!! slice%XextLinalg= converged slice eigenvectors in extended column space
!!                   of size (spacedim, neigenpairs_ext) in linalg representation (after transpose)
!! eigen = (1,neigenpairs_per_slice)=slice eigenvalues, (1,neigenpairs_per_slice+1:neigenpairs)=0
!! resid = (1,neigenpairs_per_slice)=slice residuals, (1,neigenpairs_per_slice+1:neigenpairs)=0
!! 
!! OUTPUT
!! X0     = converged eigenvector array of size (spacedim, neigenpairs)
!! eigen  = converged eigenvalues array of size (1, neigenpairs) 
!! resid  = converged residual array of size (1, neigenpairs)
!! 
!! SIDE EFFECTS
!! slice%neigenpairs_per_slice= number of kept eigenpairs after merging
!! slice%fcol_in_Xext         = first index to copy from extended memory (*,neigenpairs_ext)
!! slice%fcol_in_X            = first index to copy to regular memory (*,neigenpairs)
!!
!! SOURCE

subroutine slice_merge(slice, X0, eigen, resid)

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
    integer :: ncols_kept
    real(dp) :: part_low_bound, part_upp_bound
    ! Derived types
    type(xg_t) :: eigen_ext
    type(xg_t) :: resid_ext
    type(xgBlock_t) :: eigen_ext_slice
    type(xgBlock_t) :: resid_ext_slice
    type(xgBlock_t) :: X_kept, eigen_kept, resid_kept
    type(xgBlock_t) :: X0_out, eigen_out, resid_out
    ! arrays
    real(dp), pointer :: theta_ext(:,:)
    real(dp), allocatable :: theta_reshaped(:)
    real(dp), allocatable :: theta_ext_reshaped(:)
 
    ! *********************************************************************

    ! Sanity check
    if ( (.not. slice%use_linalg) .and. slice%use_colsrows) then
        ABI_ERROR("should be in linalg representation")
    end if
    ! TODO also check that we are within a GPU enter data map
    ! perform necessary host to device transfers

    if (.not.allocated(theta_ext_reshaped)) ABI_MALLOC(theta_ext_reshaped, (slice%neigenpairs_ext))

    ! Allocate extended space for all slice eigenvalues and residuals
    call xg_init(eigen_ext, slice%space_res, rows=1, cols=slice%neigenpairs_ext, gpu_option=slice%gpu_option)
    call xg_init(resid_ext, slice%space_res, rows=1, cols=slice%neigenpairs_ext, gpu_option=slice%gpu_option)

    call xgBlock_zero(eigen_ext%self)
    call xgBlock_zero(resid_ext%self)

    ! reshape not sure to verify TODO
    call xgBlock_reshape(eigen, 1, slice%neigenpairs)
    call xgBlock_reshape(resid, 1, slice%neigenpairs)

    ! MPI communication to gather slice eigen/resid to eigen_ext/resid_ext
    if (slice%paral_kgb==1) then
        if (xmpi_comm_size(slice%spacecom) > 1) then
        
            my_rank = xmpi_comm_rank(slice%spacecom)
            my_slice = slice%lookup_proc(my_rank + 1)
            neigen_slice = slice%neigenpairs_per_slice(my_slice)
            fcol_ext = slice%fcol_in_Xext(my_slice)

            call xgBlock_setBlock(eigen_ext, eigen_ext_slice, nrows=1, ncols=neigenpairs_slice, fcol=fcol_ext)
            call xgBlock_setBlock(resid_ext, resid_ext_slice, nrows=1, ncols=neigenpairs_slice, fcol=fcol_ext)
            call xgBlock_copy(eigen, eigen_ext_slice)
            call xgBlock_copy(resid, resid_ext_slice)

            ! All processes wait before summing 
            call xmpi_comm_barier(slice%spacecom)

            call xgBlock_mpi_sum(eigen_ext, comm=slice%spacecom)
            call xgBlock_mpi_sum(resid_ext, comm=slice%spacecom)

        end if
    else
        call xgBlock_copy(eigen, eigen_ext)
        call xgBlock_copy(resid, resid_ext)
    end if 

    if (slice%gpu_option==1) then
        call xgBlock_copy_from_gpu(eigen_ext)
        call xgBlock_copy_from_gpu(resid_ext)
    end if

    ! Results could be complex, so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(eigen_ext, theta_ext, rows=1, cols=slice%neigenpairs_ext)

    ! Filter eigenvalues in extended space using spectral partition
    tot_ncols_kept = 0
    do islice=1,slice%nslice

        ! Before merge: get slice eigenvalues to filter
        neigenpairs_slice = slice%neigenpairs_per_slice(islice)
        fcol_ext = slice%fcol_in_Xext(islice)
        lcol_ext = fcol_ext + neigenpairs_slice - 1
        if (.not.allocated(theta_reshaped)) ABI_MALLOC(theta_reshaped, (neigenpairs_slice)) 
        theta_reshaped(1:neigenpairs_slice) = theta_ext(1, fcol_ext:lcol_ext)

        ! Apply filter criterion to find kept first and last column in slice
        part_low_bound = slice%part_low_bounds(islice)
        part_upp_bound = slice%part_upp_bounds(islice)
        fcol_in_slice = maxloc(theta_reshaped, dim=1, mask=(theta_reshaped < part_low_bound)) + 1
        lcol_in_slice = maxloc(theta_reshaped, dim=1, mask=(theta_reshaped < part_upp_bound))            
        if (islice == 1     ) fcol_in_slice = 1
        if (islice == nslice) lcol_in_slice = neigenpairs_slice

        ! After merge: Update first columns to copy from Xext to X
        slice%fcol_in_X(islice)= tot_ncols_kept + 1
        slice%fcol_in_Xext(islice) = fcol_ext + fcol_in_slice - 1
        slice%neigenpairs_per_slice(islice) = lcol_in_slice - fcol_in_slice + 1
        tot_ncols_kept = tot_ncols_kept + slice%neigenpairs_per_slice(islice) 

        if (allocated(theta_reshaped)) ABI_FREE(theta_reshaped)

    end do
    
    ! Detect missing or extra eigenvalues
    if (tot_ncols_kept < slice%neigenpairs) then
        ABI_ERROR("Not enough converged eigenvalues")
    else if (tot_ncols_kept > slice%neigenpairs) then
        ABI_ERROR("Too many converged eigenvalues")
    end if

    ! Copy from extended memory to regular memory
    do islice=1,slice%nslice
        fcol = slice%fcol_in_X(islice)
        fcol_ext = slice%fcol_in_Xext(islice)
        neigenpairs_slice = slice%neigenpairs_per_slice(islice)
        ! Blocks to copy from
        call xgBlock_setBlock(slice%XextLinalg, X_kept, nrows=slice%spacedim, ncols=neigenpairs_slice, fcol=fcol_ext)
        call xgBlock_setBlock(eigen_ext, eigen_kept, nrows=1, ncols=neigenpairs_slice, fcol=fcol_ext)
        call xgBlock_setBlock(resid_ext, resid_kept, nrows=1, ncols=neigenpairs_slice, fcol=fcol_ext)
        ! Blocks to copy to
        call xgBlock_setBlock(X0, X0_out, nrows=slice%spacedim, ncols=neigenpairs_slice, fcol=fcol)
        call xgBlock_setBlock(eigen, eigen_out, nrows=1, ncols=neigenpairs_slice, fcol=fcol)
        call xgBlock_setBlock(resid, resid_out, nrows=1, ncols=neigenpairs_slice, fcol=fcol)
        ! copy
        call xgBlock_copy(X_kept, X0_out)
        call xgBlock_copy(eigen_kept, eigen_out)
        call xgBlock_copy(resid_kept, resid_out)
    end do

    ! reshape not sure to verify TODO
    call xgBlock_reshape(eigen, slice%neigenpairs, 1) 
    call xgBlock_reshape(resid, slice%neigenpairs, 1) 

    ! Free memory
    call ABI_FREE(eigen_ext)
    call ABI_FREE(resid_ext)
    if (allocated(theta_ext_reshaped)) ABI_FREE(theta_ext_reshaped)
 
end subroutine slice_merge
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
    integer :: i
    real(dp) :: total_weight, lower, upper, mid, total_allocated, multiplier
    integer :: allocation(n)

    ! *********************************************************************

    ! Calculate total weight
    total_weight = real(sum(w))

    ! Binary search for the optimal multiplier
    lower = 0.0
    upper = real(p)
    do while (upper - lower > 1.0)
        mid = (lower + upper) / 2.0
        total_allocated = 0.0
        do i = 1, n
            allocation(i) = int((real(w(i)) * real(m(i)) / total_weight) * mid + 0.5)
            total_allocated = total_allocated + allocation(i)
        end do

        ! Adjust binary search bounds
        if (total_allocated > real(p)) then
            upper = mid
        else
            lower = mid
        end if
    end do

    ! Final allocation after binary search converges
    multiplier = (lower + upper) / 2.0
    do i = 1, n
        allocation(i) = int((real(w(i)) * real(m(i)) / total_weight) * multiplier + 0.5)
    end do

    ! Adjust total allocation to exactly match p
    total_allocated = 0
    do i = 1, n
        total_allocated = total_allocated + allocation(i)
    end do

    if (total_allocated < p) then
        do while (total_allocated < p)
            ! Add one resource to the group closest to its ideal allocation
            call adjust_allocation(n, m, w, allocation, total_weight, p, total_allocated)
            total_allocated = 0
            do i = 1, n
                total_allocated = total_allocated + allocation(i)
            end do
        end do
    else if (total_allocated > p) then
        do while (total_allocated > p)
            ! Remove one resource from the over-allocated group
            call reduce_allocation(n, m, w, allocation, total_weight, p, total_allocated)
            total_allocated = 0
            do i = 1, n
                total_allocated = total_allocated + allocation(i)
            end do
        end do
    end if

    ! Assign the final allocation to the output variable
    x = allocation

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
    integer, intent(in) :: w(n), total_weight, p
    integer, intent(inout) :: allocation(n)
    integer, intent(inout) :: total_allocated

    integer :: i, closest_group
    real(dp) :: max_diff, diff
    
    ! *********************************************************************

    ! Find the group with the largest difference between allocation and ideal allocation
    max_diff = -1.0
    closest_group = 1
    do i = 1, n
        diff = abs(real(allocation(i)) - (real(w(i)) * real(m(i)) * real(p)) / real(total_weight))
        if (diff > max_diff) then
            max_diff = diff
            closest_group = i
        end if
    end do

    ! Add one resource to the group with the largest difference
    allocation(closest_group) = allocation(closest_group) + 1
    total_allocated = total_allocated + 1

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
    integer, intent(in) :: w(n), total_weight, p
    integer, intent(inout) :: allocation(n)
    integer, intent(inout) :: total_allocated

    integer :: i, closest_group
    real(dp) :: max_diff, diff
    
    ! *********************************************************************

    ! Find the group with the smallest over-allocation
    max_diff = -1.0
    closest_group = 1
    do i = 1, n
        diff = abs(real(allocation(i)) - (real(w(i)) * real(m(i)) * real(p)) / real(total_weight))
        if (diff < max_diff) then
            max_diff = diff
            closest_group = i
        end if
    end do

    ! Remove one resource from the group with the smallest over-allocation
    allocation(closest_group) = allocation(closest_group) - 1
    total_allocated = total_allocated - 1

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

function print_scalar_filter(upp, low, center, radius, npt, is_lowpass)

    implicit none

    real(dp), intent(in) :: upp
    real(dp), intent(in) :: low
    real(dp), intent(in) :: center
    real(dp), intent(in) :: radius
    integer, intent(in) :: npt
    logical, intent(in) :: is_lowpass

            write(std_out,*) ' '
            write(std_out,*) 'Plot filter ==== x | f(x)'
            npt = 100
            if (islice==1) then
                do ipt=1,npt
                    pt = low + (ipt-1)*(upp-low)/npt
                    fun_pt = cheb_poly(pt,ndeg,poly_upp,gub)
                    write(std_out,*) pt, fun_pt
                end do
            else
                do ipt=1,npt
                    pt = (low + (ipt-1)*(upp-low)/npt - c)/r
                    fun_pt = bandpassIndicator_sca(pt,a,b,ndeg)
                    write(std_out,*) pt, fun_pt
                end do
            end if
            write(std_out,*) ' '




end function print_scalar_filter
!!***

end module m_slice
!!***

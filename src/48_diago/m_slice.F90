!!****f* ABINIT/m_slice
!! NAME
!! m_slice
!!
!! FUNCTION
!! This module contains the types and routines used to apply 
!! the Spectrum Slicing method. It mainly defines 'slice' 
!! datatypes and associated methods. Features:
!! - uses xgTools implementation as matrix data structure.
!! - based on 'chebfi' data structure for most vector routines.
!! - polynomial filters are Chebyshev for first slice and
!!   Chebyshev-Jackson expansion of indicator otherwise.
!! - implements new parallel level between slices.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (IML, LB)
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

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
    use m_nvtx_data
#endif

    implicit none

    private

    ! Load balance criteria for fair resource allocation
    !---------------------------------------------------
    integer, parameter :: FAIR_BANDPP      = 0 ! bandpp (unweigted)
    integer, parameter :: FAIR_BANDPP_WDEG = 1 ! bandpp weighted by degree 

    ! Private 'sliceTasks' datatype
    ! Parameters of individual slice tasks executed in parallel (scheduler)
    !-------------------------------------------------
    type, private :: sliceTasks_t

        ! Slice parameter object
        type(slice_t), allocatable :: slice

        ! MPI-related
        integer :: ntasks        ! number of tasks (=slices)
        integer :: nprocs_tot    ! total number of available resources (process)
        integer :: task_me       ! which task current process serves
        integer :: ncols_me      ! number of columns of current task
        integer :: nprocs_me     ! number of resources used by current task
        integer :: comm_global   ! global communicator (inter-task)

        ! Data buffer
        integer :: nrows_tot  ! total number of rows (planewaves)

        ! Flags
        logical :: on_host = .false.
        logical :: on_device = .false.
        logical :: use_blockrows = .false.
        logical :: use_blockcols = .false.

        ! Options
        integer :: resource_allocation

        ! Arrays
        integer, allocatable :: tasks_ncols(:)   ! number of columns per task
        integer, allocatable :: tasks_nprocs(:)  ! number of processes per task
        integer, allocatable :: assigned_task(:) ! which task each process serves
        integer, allocatable :: blockcols_me(:)  ! ncol of colsrows-blocks in my task
        integer, allocatable :: blockrows_me(:)  ! nrow of linalg-blocks in my task

    end type sliceTasks_t

    ! Public 'cgSliced' datatype
    ! cg does not have a distribution/structure suitable for parallel spectrum slicing
    ! cgSliced is a version of cg, distributed correctly and free of data overlap
    !-------------------------------------------------
    type, public :: cgSliced_t

        ! MPI-related
        integer :: spacecom
        type(sliceTasks_t), private :: schedule

        ! Pointers to existing memory
        type(xgBlock_t) :: Xblockcols     ! if use_blockcols=.true.
        type(xgBlock_t) :: Xblockrows     ! if use_blockrows=.true.
        
        ! Proper memory (extension of cg)
        type(xg_t) :: xgXblockrows
        
        ! Transposer
        type(xgTransposer_t) :: xgTransposerX
    
    end type cgSliced_t

    ! Public 'slice' datatype
    ! [IML 9/4 work in progress] I don't think this is necessary
    ! Parameters specific to individual slices
    !-------------------------------------------------
    type, public :: sliceWorker_t

        integer :: islice                       ! slice index in 1,..,nslice
        integer :: nband                        ! number of bands in slice
        integer :: i1                           ! first band index in global
        integer :: i2                           ! last band index in global
        integer :: degree                       ! filter degree
        integer :: bandpp                       ! number of slice bands per process
        real(dp) :: ramp                        ! filter amplification factor
        real(dp) :: low                         ! filter support lower bound
        real(dp) :: upp                         ! filter support upper bound
        real(dp) :: glb                         ! spectrum lower bound
        real(dp) :: gub                         ! spectrum upper bound

        ! Memory allocated for a slice
        type(xg_t) :: XW                        ! input/ output eigenvectors
        type(xg_t) :: EW                        ! output eigenvalues
        type(xg_t) :: RW                        ! output residuals
        type(xg_t) :: OCCW                      ! unused memory space in Slicing
        type(chebfi_t) :: chebfi                ! workspace for slice diago
   
        ! Pointers to slice memory
        type(xgBlock_t) :: xgx0
        type(xgBlock_t) :: xgeigen
        type(xgBlock_t) :: xgresidu
        type(xgBlock_t) :: xgocc ! not used at all

    end type sliceWorker_t

    ! Public 'slice' datatype
    ! Parameters common to all slices
    !-------------------------------------------------
    type, public :: sliceAll_t

        integer :: nslice                        ! Total number of slices
        integer :: space
        integer :: space_res
        integer :: spacedim                      ! Space dimension for one vector
        integer :: total_spacedim                ! Maybe not needed
        integer :: neigenpairs                   ! Number of eigen values/vectors we want
        integer :: mdeg_filter                   ! Degree of the polynomial filter
        integer :: spacecom                      ! Communicator for MPI
        real(dp) :: tolerance                    ! Tolerance on the residu to stop the minimization
        real(dp) :: ecut                         ! Ecut used for Polynomial filtering
       
        ! Variables deactivated for chebfi
        integer :: oracle = 0                        
        integer :: nbdbuf = 0                       
        real(dp) :: oracle_factor = 1.d0                
        real(dp) :: oracle_min_occ = 0.d0              

        integer :: paral_kgb                     ! allow MPI distribution of k-points, bands or PW
        integer :: paral_slice                   ! distribute slices over MPI processes
        integer :: bandpp                        ! nbands per MPI process
        integer :: comm_cols                     ! MPI column communicator
        integer :: comm_rows                     ! MPI row communicator
        integer :: me_g0           
        integer :: me_g0_fft

        ! General slicing params
        integer :: nband_ovlp                        ! number of columns in overlap-free space
        integer :: npband                            ! number of MPI processes
        integer :: spectral_cut                         ! spectral partition technique
        real(dp) :: glb                              ! global lower spectral bound
        real(dp) :: gub                              ! global upper spectral bound
        real(dp) :: ramp                             ! amplification factor for slices
        
        logical :: paw
        integer :: eigenProblem   ! 1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)

        ! when GPU is enabled, currently OpenMP is not fully supported, abinit is launched
        ! with OMP_NUM_THREADS=1, but we may locally increase the number of OpenMP threads
        ! wherever it is safe to do; in that case we use gpu_kokkos_nthrd to specify
        ! the number of OpenMP threads. This value is controlled by dtset variable
        ! dtset%gpu_kokkos_nthrd
        integer :: gpu_option
        integer :: gpu_kokkos_nthrd = 1 ! only used if gpu is enabled, number of OpenMP threads used
        integer :: gpu_thread_limit = 1 ! only used if GPU is enabled, max number of OpenMP threads used in sensitive areas

        ! DOS (all eigenvalues+residuals) before Slicing. Not distributed across MPI proc. 
        ! Useful for convergence test
        type(xg_t) :: Eig0                           ! Rayleigh quotients before Slicing
        type(xg_t) :: Res0                           ! residual norms before Slicing

        ! Initial memory space (in/out)
        type(xgBlock_t) :: X
        type(xgBlock_t) :: eigen
        type(xgBlock_t) :: residu

        ! Extended memory space
        type(xg_t) :: Xext_                        ! input/ output eigenvectors
        type(xg_t) :: Eext_                        ! output eigenvalues
        type(xg_t) :: Rext_                        ! output residuals

        ! Pointers to extended memory (linalg representation)
        type(xgBlock_t) :: Xext_linalg
        type(xgBlock_t) :: Eext_linalg
        type(xgBlock_t) :: Rext_linalg
        ! colsrows representation
        type(xgBlock_t) :: Xext

        ! Transposer
        type(xgBlock_t) :: Transposer_Xext

        integer, allocatable :: permute_cols(:)        ! eigenvector permutation (TODO move out)
        integer, allocatable :: slice_fcol(:)          ! first col of slice in spectrum memory
        integer, allocatable :: slice_ncols(:)         ! number of cols per slice
        integer, allocatable :: slice_fcol_ext(:)      ! first col of slice in extended memory
        integer, allocatable :: conv_fcol_slice(:)     ! first converged in slice memory
        integer, allocatable :: poly_degrees()         ! polynomial filter degrees
        real(dp), allocatable :: part_low_bounds(:)    ! lower bounds in spectral partition
        real(dp), allocatable :: part_upp_bounds(:)    ! upper bounds in spectral partition
        real(dp), allocatable :: poly_low_bounds(:)    ! lower bounds used to define polynomials
        real(dp), allocatable :: poly_upp_bounds(:)    ! upper bounds used to define polynomials

    end type slice_t

    ! Public methods associated to 'slice' datatype
    !-------------------------------------------------
    public :: cgSliced_init                     ! initialize cgSliced data type object
    public :: cgSliced_buildFrom                ! build cgSliced from cg
    public :: cgSliced_mergeTo                  ! merge cgSliced to cg
    public :: cgSliced_free                     ! free cgSliced data type object
    public :: slice_init                        ! initiate slice data type object
    public :: slice_run                         ! diagonalize individual slices
    public :: slice_free                        ! free slice data type object
    public :: slice_mergeConverged              ! merge converged slices by removing duplicates
    public :: slice_unitTest                    ! used for debugging

    CONTAINS  
!=====================================================================
!!***

!!****f* m_slice/cgSliced_init
!! NAME
!! cgSliced_init
!!
!! FUNCTION
!! Initialize a 'cgSliced' datastructure. See chebfi_init().
!!
!! SOURCE

subroutine cgSliced_init(cgSliced,neigenpairs,spacedim,tolerance,ecut,&
        paral_kgb,paral_slice,bandpp,space,spacecom,me_g0,me_g0_fft,paw,&
        comm_rows,comm_cols,nslice,ramp,spectral_cut,gpu_option,&
        gpu_kokkos_nthrd,gpu_thread_limit)

    implicit none

    ! Arguments ------------------------------------
    integer      , intent(in   ) :: bandpp
    integer      , intent(in   ) :: npband
    integer      , intent(in   ) :: nslice
    integer      , intent(in   ) :: eigenProblem
    integer      , intent(in   ) :: me_g0
    integer      , intent(in   ) :: me_g0_fft
    integer      , intent(in   ) :: neigenpairs
    integer      , intent(in   ) :: mdeg_filter
    integer      , intent(in   ) :: comm_cols
    integer      , intent(in   ) :: comm_rows
    integer      , intent(in   ) :: paral_kgb
    integer      , intent(in   ) :: paral_slice
    integer      , intent(in   ) :: space
    integer      , intent(in   ) :: spacecom
    integer      , intent(in   ) :: spacedim
    integer      , intent(in   ) :: spectral_cut
    integer      , intent(in   ) :: gpu_option
    logical      , intent(in   ) :: paw
    real(dp)     , intent(in   ) :: ramp
    real(dp)     , intent(in   ) :: ecut
    real(dp)     , intent(in   ) :: tolerance
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ), optional :: gpu_kokkos_nthrd
    integer      , intent(in   ), optional :: gpu_thread_limit

    ! *********************************************************************

    slice%space        = space
    slice%neigenpairs  = neigenpairs
    slice%spacedim     = spacedim
    slice%tolerance    = tolerance
    slice%ecut         = ecut
    slice%paral_kgb    = paral_kgb
    slice%paral_slice  = paral_slice
    slice%comm_cols    = comm_cols
    slice%bandpp       = bandpp
    slice%comm_rows    = comm_rows
    slice%mdeg_filter  = mdeg_filter
    slice%spacecom     = spacecom
    slice%eigenProblem = eigenProblem
    slice%me_g0        = me_g0
    slice%me_g0_fft    = me_g0_fft
    slice%paw          = paw
    slice%gpu_option   = gpu_option
    slice%nslice       = nslice
    slice%npband       = xmpi_comm_size(comm_cols)
    slice%ramp         = ramp
    slice%spectral_cut = spectral_cut
    slice%nband_ovlp   = neigenpairs ! see initExtended
    slice%nbdbuf       = nbdbuf

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

    ! Dimensions of linalg representation
    cgslice%nrows_blockrows = spacedim
    cgslice%ncols_blockrows = neigenpairs
    ! Dimensions of colsrows representation
    if (paral_kgb == 0) then
        ! No distribution
        cgslice%nrows_blockcols = spacedim
        cgslice%ncols_blockcols = neigenpairs
    else if (paral_kgb == 1) then
        total_spacedim = spacedim
        call xmpi_sum(total_spacedim,cgslice%spacecom,ierr)
        cgslice%nrows_blockcols = total_spacedim
        cgslice%ncols_blockcols = bandpp
    end if

end subroutine cgSliced_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/cgSliced_buildFrom
!! NAME
!! cgslice_buildFrom
!! 
!! FUNCTION
!! Split then distribute vectors X across processes.
!! Then create extended memory buffers and distribute them.
!! 
!! SOURCE

subroutine cgSlice_divide(cgslice,X,getAX_BX,nspinor)

    implicit none

    ! Arguments    
    type(cgSlice_t), intent(inout) :: cgslice
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

    ! Local variables
    integer :: resource_allocation
    integer :: nrows_total
    integer :: spacecom
    type(xg_t) :: eigen0
    type(xg_t) :: resid0
    ! Arrays
    integer, pointer :: slice_ncols(:) => null()
    integer, pointer :: slice_degrees(:) => null()
    
    ! *********************************************************************

    space_res = cgslice%space_res
    ncols_tot = cgslice%ncols_blockrows
    gpu_option = cgslice%gpu_option

    ! There are the ENTIRE values used for later. Every process contains entire array
    call xg_init(eigen0, space_res, rows=1, cols=ncols_tot, gpu_option=gpu_option)
    call xg_init(resid0, space_res, rows=1, cols=ncols_tot, gpu_option=gpu_option)
    call xgBlock_zero(eigen0%self)
    call xgBlock_zero(resid0%self)

    ! Compute Rayleigh quotients and residuals
    call cgslice_computeSpectrum(cgslice,X,getAX_BX,eigen0,resid0,nspinor)
    ! now eigen0%self holds thetas
    ! en vrai pour cutSpectrum on n'a pas besoin de theta.
    ! juste du maxeig, mineig et mineig_pos pour faire resid(mineig_pos)
    ! après il faut aussi faire la permutation du X SUR CPU en LINALG!!!!

    ! Map to fortran arrays
    call xgBlock_reshape(eigen0%self, (/ncols_tot,1/))     
    call xgBlock_reshape(resid0%self, (/ncols_tot,1/))

    ! normally here we should move some data to CPU

    ! Compute nband per slice
    call slice_cutSpectrum(slice,pband_ptr,spectral_cut)

    ! TODO actually store these parameters into the sliceTasks object

    nrows_total = slice%total_spacedim
    spacecom = slice%spacecom

    if (paral_slice==0) then
        ABI_ERROR("Sequential slices not implemented")
    else if (paral_slice==1) then
        resource_allocation = FAIR_BANDPP
    else if (paral_slice==2) then
        resource_allocation = FAIR_BANDPP_WDEG
    end if

    ! Compute parameters of slice tasks
    ! TODO associate schedule into slice in some way?
    ! TODO free scheduler at some point
    slice_ncols => slice%slice_ncols
    slice_degrees => slice%slice_degrees
    call sliceTasks_init(schedule,slice_ncols,slice_degrees,nrows_total,&
        spacecom,resource_allocation)

    ! Create the extended workspaces on CPU
    call slice_initExtended(slice,pband_ptr)

    ! Distribute extended space across MPI processes
    call slice_distributeExtended(slice)



end subroutine cgslice_divide
!!***

!! FUNCTION
!! Compute ordered Rayleigh quotients and residuals
!! TODO IL 11/4 restore chebfi from develop
!! this is a new version without chebfi

subroutine cgslice_computeSpectrum(cgslice,X0,getAX_BX,eigen,resid,nspinor)

    implicit none

    ! Arguments
    type(cgSlice_t), intent(inout) :: cgslice
    type(xgBlock_t), intent(inout) :: X0
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
    integer :: space, space_res
    integer :: nrows_blockcols, ncols_blockcols
    integer :: paral_kgb, comm
    integer :: me_g0, gpu_option
    integer :: my_rank, shift, bandpp
    real(dp) :: mineig, maxeig
    ! Derived types
    type(xg_t) :: Results1
    type(xg_t) :: Results2
    type(xg_t) :: Results3
    type(xg_t) :: eigen_mpi
    type(xg_t) :: resid_mpi
    type(xg_t) :: X_NAB ! vector memory for X_next, AX, BX
    type(xgBlock_t) :: xXColsRows
    type(xgBlock_t) :: xAXColsRows
    type(xgBlock_t) :: xBXColsRows
    type(xgBlock_t) :: X_next
    type(xgTransposer_t) :: xgTransposerX
    ! Arrays
    integer :: maxeig_pos(2)
    integer :: mineig_pos(2)

! *********************************************************************

    ! Priority IL 11/4
    ! TODO rename cgslice to xgXslice or something. The structure is xg not cg!!

    space = cgslice%space
    space_res = cgslice%space_res
    paral_kgb = cgslice%paral_kgb
    nrows_blockcols = cgslice%nrows_blockcols
    ncols_blockcols = cgslice%ncols_blockcols
    comm = cgslice%comm
    me_g0 = cgslice%me_g0
    gpu_option = cgslice%gpu_option

    ! Allocate temporary memory space
    call xg_init(X_NAB, space, nrows_blockcols, 3*ncols_blockcols, comm, me_g0=me_g0, gpu_option=gpu_option)
    call xg_setBlock(X_NAB, X_next, nrows_blockcols, ncols_blockcols)
    call xg_setBlock(X_NAB, xAXColsRows, nrows_blockcols, ncols_blockcols, fcol=ncols_blockcols + 1)
    call xg_setBlock(X_NAB, xBXColsRows, nrows_blockcols, ncols_blockcols, fcol=2*ncols_blockcols + 1)
    
    ! one-dimensional memory spaces
    call xg_init(eigen_mpi, space_res, rows=ncols_blockcols, cols=1, comm=comm, gpu_option=gpu_option)
    call xg_init(resid_mpi, SPACE_R, rows=ncols_blockcols, cols=1, comm=comm, gpu_option=gpu_option)
    call xg_init(Results1, space_res, rows=ncols_blockcols, cols=1, gpu_option=gpu_option)
    call xg_init(Results2, space_res, rows=ncols_blockcols, cols=1, gpu_option=gpu_option)
    call xg_init(Results3, SPACE_R, rows=ncols_blockcols, cols=1, gpu_option=gpu_option)

    ! First transpose existing cg
    if (paral_kgb == 1) then
        call xmpi_barrier(comm)
        call xgTransposer_constructor(xgTransposerX,X0,xXColsRows,nspinor,&
            STATE_LINALG,TRANS_ALL2ALL,xmpi_comm_self,comm,0,0,cgslice%me_g0_fft,&
            gpu_option=gpu_option,gpu_thread_limit=cgslice%gpu_thread_limit)
        xgTransposerX%gpu_kokkos_nthrd  = cgslice%gpu_kokkos_nthrd
        call xgTransposer_transpose(xgTransposerX,STATE_COLSROWS)
    end if

    ! Now apply AX BX (colsrows representation)
    ! Remember that this function will copy X to BX if paw
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
    call getAX_BX(xXColsRows,xAXColsRows,xBXColsRows)
    call xgBlock_zero_im_g0(xAXColsRows)
    call xgBlock_zero_im_g0(xBXColsRows)
    ABI_NVTX_END_RANGE()

    ! Compute Rayleigh quotients
    ! <Psi|H|Psi>
    call xgBlock_colwiseDotProduct(xXColsRows, xAXColsRows, Results1%self, comm_loc=xmpi_comm_null)
    ! <Psi|S|Psi>
    call xgBlock_colwiseDotProduct(xXColsRows, xBXColsRows, Results2%self, comm_loc=xmpi_comm_null)
    ! eigen = <Psi|H|Psi> / <Psi|S|Psi>
    call xgBlock_colwiseDivision(Results1%self, Results2%self, eigen_mpi, maxeig, maxeig_pos, mineig, mineig_pos)

    ! In order to avoid transposing AX,BX, we compute residuals in colsrows representation
    ! TODO IL 20/01/2025 ymax has not been tested on GPU
    call xgBlock_copy(xBXColsRows,X_next)                             ! X_next = S|Psi>
    call xgBlock_ymax(X_next,eigen_mpi,0,1)                           ! X_next = - eig * S|Psi>
    call xgBlock_add(X_next,xAXColsRows)                              ! X_next = H|Psi> - eig * S|Psi>
    call xgBlock_colwiseNorm2(%X_next,residu,comm_loc=xmpi_comm_null) ! resid = |X_next|^2
   
    ! MPI communication for gathering all thetas (using summation strategy)
    if (xmpi_comm_size(comm)>1) then
        my_rank = xmpi_comm_rank(comm)
        bandpp = ncols_blockcols
        shift = my_rank * bandpp

        call xgBlock_setBlock(eigen%self, Results1%self, nrows=1, ncols=bandpp, fcol=1+shift)
        call xgBlock_reshape(eigen_mpi, (/1,bandpp/)) 
        call xgBlock_copy(eigen_mpi, Results1%self)
        call xgBlock_mpi_sum(eigen%self,comm=comm)

        call xgBlock_setBlock(resid%self, Results3%self, nrows=1, ncols=bandpp, fcol=1+shift)
        call xgBlock_reshape(resid_mpi, (/1,bandpp/))     
        call xgBlock_copy(resid_mpi, Results3%self)
        call xgBlock_mpi_sum(resid%self,comm=comm)
    else
        call xgBlock_copy(eigen_mpi,eigen%self)
        call xgBlock_copy(resid_mpi,resid%self)
    end if

    ! Restore cg to linalg representation
    if (paral_kgb == 1) then
        call xmpi_barrier(comm)
        call xgTransposer_transpose(xgTransposerX, STATE_LINALG)
    end if

    ! Free memory
    call xg_free(Results1)
    call xg_free(Results2)
    call xg_free(Results3)
    call xg_free(eigen_mpi)
    call xg_free(resid_mpi)
    call xg_free(X_NAB)

end subroutine cgslice_computeSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_init
!! NAME
!! slice_init
!!
!! FUNCTION
!! Initialize a 'slice' datastructure. See chebfi_init().
!!
!! SOURCE

subroutine slice_init(slice,neigenpairs,spacedim,tolerance,ecut,paral_kgb,paral_slice,&
        bandpp,mdeg_filter,space,eigenProblem,spacecom,me_g0,me_g0_fft,&
        paw,comm_rows,comm_cols,nslice,ramp,balance,gpu_option,&
        gpu_kokkos_nthrd,gpu_thread_limit)

    implicit none

    ! Arguments ------------------------------------
    integer      , intent(in   ) :: bandpp
    integer      , intent(in   ) :: npband
    integer      , intent(in   ) :: nslice
    integer      , intent(in   ) :: eigenProblem
    integer      , intent(in   ) :: me_g0
    integer      , intent(in   ) :: me_g0_fft
    integer      , intent(in   ) :: neigenpairs
    integer      , intent(in   ) :: mdeg_filter
    integer      , intent(in   ) :: comm_cols
    integer      , intent(in   ) :: comm_rows
    integer      , intent(in   ) :: paral_kgb
    integer      , intent(in   ) :: paral_slice
    integer      , intent(in   ) :: space
    integer      , intent(in   ) :: spacecom
    integer      , intent(in   ) :: spacedim
    integer      , intent(in   ) :: balance
    integer      , intent(in   ) :: gpu_option
    logical      , intent(in   ) :: paw
    real(dp)     , intent(in   ) :: ramp
    real(dp)     , intent(in   ) :: ecut
    real(dp)     , intent(in   ) :: tolerance
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ), optional :: gpu_kokkos_nthrd
    integer      , intent(in   ), optional :: gpu_thread_limit

    ! *********************************************************************

    slice%space        = space
    slice%neigenpairs  = neigenpairs
    slice%spacedim     = spacedim
    slice%tolerance    = tolerance
    slice%ecut         = ecut
    slice%paral_kgb    = paral_kgb
    slice%paral_slice  = paral_slice
    slice%comm_cols    = comm_cols
    slice%bandpp       = bandpp
    slice%comm_rows    = comm_rows
    slice%mdeg_filter  = mdeg_filter
    slice%spacecom     = spacecom
    slice%eigenProblem = eigenProblem
    slice%me_g0        = me_g0
    slice%me_g0_fft    = me_g0_fft
    slice%paw          = paw
    slice%gpu_option   = gpu_option
    slice%nslice       = nslice
    slice%npband       = xmpi_comm_size(comm_cols)
    slice%ramp         = ramp
    slice%spectral_cut = balance
    slice%nband_ovlp   = neigenpairs ! see initExtended
    slice%nbdbuf       = nbdbuf

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

    call slice_allocateAll(slice)
    
end subroutine slice_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allocateAll
!! NAME
!! slice_allocateAll
!! 
 
subroutine slice_allocateAll(slice)

    implicit none
    
    type(slice_t), intent(inout) :: slice
    integer :: space_res
    integer :: neigenpairs
    integer :: comm_cols
    integer :: gpu_option

    call slice_free(slice)

    space_res = slice%space_res
    neigenpairs = slice%neigenpairs
    comm_cols = slice%comm_cols
    !gpu_option = slice%gpu_option
    ! FIXME forced CPU
    gpu_option = ABI_GPU_DISABLED

    ! Eigenvalues and residuals before slicing
    ! Every MPI process contains this array
    ! IML 9/4 FIXME this is useful if we compute convergence rates aposteriori
    ! otherwise, I am not sure how to use this information
    ! IML 11/4 TODO maybe we don't need to store this information as member. 
    ! Try to localize as much as possible to use in required routines only
    call xg_init(slice%Eig0,space_res,rows=1,cols=neigenpairs,comm=comm_cols,gpu_option=gpu_option)
    call xg_init(slice%Res0,SPACE_R,rows=1,cols=neigenpairs,comm=comm_cols,gpu_option=gpu_option)

    if(.not.allocated(slice%slice_fcol)) ABI_MALLOC(slice%slice_fcol,(nslice))
    if(.not.allocated(slice%slice_ncols)) ABI_MALLOC(slice%slice_ncols,(nslice))
    if(.not.allocated(slice%slice_fcol_ext)) ABI_MALLOC(slice%slice_fcol_ext,(nslice))
    if(.not.allocated(slice%poly_degrees)) ABI_MALLOC(slice%poly_degrees,(nslice))
    if(.not.allocated(slice%part_low_bounds)) ABI_MALLOC(slice%part_low_bounds,(nslice))
    if(.not.allocated(slice%part_upp_bounds)) ABI_MALLOC(slice%part_upp_bounds,(nslice))
    if(.not.allocated(slice%poly_low_bounds)) ABI_MALLOC(slice%poly_low_bounds,(nslice))
    if(.not.allocated(slice%poly_upp_bounds)) ABI_MALLOC(slice%poly_upp_bounds,(nslice))

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_free
!! NAME
!! slice_free
!! 
 
subroutine slice_free(slice)

    implicit none
    
    type(slice_t), intent(inout) :: slice

    call slice_freeOverlapFree(slice)
    call xg_free(slice%Eig0)
    call xg_free(slice%Res0)
 
    ! Free slice parameters
    if(allocated(slice%slice_fcol)) ABI_FREE(slice%slice_fcol)
    if(allocated(slice%slice_ncols)) ABI_FREE(slice%slice_ncols)
    if(allocated(slice%slice_fcol_ext)) ABI_FREE(slice%slice_fcol_ext)
    if(allocated(slice%poly_degrees)) ABI_FREE(slice%poly_degrees)
    if(allocated(slice%part_low_bounds)) ABI_FREE(slice%part_low_bounds)
    if(allocated(slice%part_upp_bounds)) ABI_FREE(slice%part_upp_bounds)
    if(allocated(slice%poly_low_bounds)) ABI_FREE(slice%poly_low_bounds)
    if(allocated(slice%poly_upp_bounds)) ABI_FREE(slice%poly_upp_bounds)
    
end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_cutSpectrum
!! NAME
!! slice_cutSpectrum
!! 
!! FUNCTION
!! Split spectrum into overlapping slices. 
!!
!! INPUT
!! sliceAll=         parameters common to all slices, such as
!!                   balance_option=balance number of vectors(1)
!!                                  balance filter degrees(2)
!!                   ramp=filter convergence threshold, greater than 1
!! idxAll=           start and end index per slice
!! ndegAll=          filter degree per slice
!! sboundAll=        various spectral bounds
!! pband=            eigenvector permutation
!!
!! SOURCE

subroutine slice_cutSpectrum(slice,nband,idxAll,ndegAll,sboundAll,pband,npbandSlice,plot_filter)

    implicit none

    !Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: nband
    integer, pointer, intent(inout) :: idxAll(:,:)
    integer, pointer, intent(inout) :: ndegAll(:)
    integer, pointer, intent(inout) :: pband(:)
    integer, pointer, intent(inout) :: npbandSlice(:)
    real(dp), pointer, intent(inout) :: sboundAll(:,:)
    logical, optional, intent(in) :: plot_filter
    
    !Local variables-------------------------------
    integer :: j,k,jmax,spos,nvec_ovlp,nvec,k1,k2
    integer :: nv_pad,ndeg,npband,comm_cols
    integer :: paral_slice
    integer :: nslice,neigenpairs,nline
    integer :: ndeg_max = 200
    integer :: balance_option
    integer :: ipt,npt,iptL,iptR
    logical :: plot_filter_
    real(dp) :: ramp
    real(dp) :: tol12 = 1.0e-12
    real(dp) :: ecut,low,upp,glb,gub,c,r
    real(dp) :: lj,uj,wj,finL,finR,foutL,foutR
    real(dp) :: wlj,wuj
    real(dp) :: fun_pt,pt
    real(dp) :: a_,b_ ! target interval scaled in -1,1
    type(xgBlock_t) :: Eig0_all, Res0_all
    ! arrays
    integer :: jperm(nband-1)
    real(dp) :: consdiff(nband-1)
    real(dp), allocatable :: slice_cut(:)
    real(dp), allocatable :: resid_cut(:)
    real(dp), pointer :: resid_(:,:)
    real(dp), pointer :: theta_(:,:)
    real(dp), allocatable, target :: resid(:)
    real(dp), allocatable, target :: theta(:)

    ! *********************************************************************

    plot_filter_ = .false.
    if (present(plot_filter)) plot_filter_ = plot_filter

    ! Interval slicing parameters
    npband         = slice%npband ! number of MPI processes
    bandpp         = slice%bandpp ! fixed number of bands per MPI process
    comm_cols      = slice%comm_cols
    neigenpairs    = slice%neigenpairs
    nslice         = slice%nslice
    ecut           = slice%ecut
    nline          = slice%mdeg_filter
    paral_slice    = slice%paral_slice
    balance_option = slice%spectral_cut
    ramp           = slice%ramp

    ! Set pointers to eigenvalue and residual memory
    Eig0_all = slice%Eig0%self
    Res0_all = slice%Res0%self

    ! Results could be complex, so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(Eig0_all,theta_,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(Res0_all,resid_,rows=1,cols=neigenpairs)

    ! Save theta_,resid_ first row in theta,resid
    ABI_MALLOC(theta,(neigenpairs))
    ABI_MALLOC(resid,(neigenpairs))
    theta(1:neigenpairs) = theta_(1,1:neigenpairs)
    resid(1:neigenpairs) = resid_(1,1:neigenpairs)

    ! Sort thetas in increasing order and store permutation
    call sort_dp(neigenpairs,theta,pband,tol12)
    
    !write(std_out,*) 'Rayleigh Values after sort:'
    !write(std_out,*) theta(:)

    ! Spectrum bounds
    low = theta(1)
    upp = theta(neigenpairs)
    ! Bounds that contain entire spectrum (guaranteed)
    glb = low - sqrt(resid(pband(1)))
    gub = ecut 
    ! ecut too large resulting in fine slices, narrow down:
    !gub = upp + sqrt(resid(pband(neigenpairs))) + 2.d0
    ! this results in eigenvalues between gub and ecut

    ! Center and radius of the entire spectrum mapped to -1,1
    c = (glb + gub) / 2.d0
    r = (gub - glb) / 2.d0
    
    write(std_out,*) ' '
    write(std_out,*) '------------ A priori spectrum '
    write(std_out,*) 'Spectral bounds=', glb,gub
    write(std_out,*) ' Minmax eigvals=', low,upp
    write(std_out,*) '           ecut=', ecut
    write(std_out,*) '         center=', c
    write(std_out,*) '         radius=', r
    write(std_out,*) ' '

    ! Store to slicing object
    slice%gub = gub
    slice%glb = glb
   
    ! Initialize spectral cuts
    ABI_MALLOC(slice_cut,(nslice+1))
    slice_cut(:) = 0.d0
    slice_cut(1) = low
    slice_cut(nslice+1) = upp
    
    ! Compute spectral cuts on interior slices
    select case(balance_option)
    case(1)
        ! Balance interval widths
        slice_cut(2:nslice) = (/ (low+(upp-low)/nslice*j, j=1,nslice-1) /)
    case(2)
        ! Balance number of vectors
        slice_cut(2:nslice) = (/ (theta(neigenpairs/nslice*j), j=1,nslice-1) /)
    case(3)
        ! Cut on spectral gaps=where eigenvalues are less concentrated
        ! Sort consecutive differences by increasing order
        jperm = (/ (k, k=1,nband-1) /)
        consdiff = (/ (theta(k+1) - theta(k), k=1,nband-1) /)
        call sort_dp(nband-1,consdiff,jperm,tol12)
        ! Take median of largest gaps
        do i=1,nslice
            jmax = jperm(nband-i)
            slice_cut(i+1) = (theta(jmax) + theta(jmax+1)) / 2.d0
        end do
    end select

    ! Define spectral subintervals and optimize degrees for individual slices
    j1 = 1 ! band index
    do islice=1,nslice
            
        ! Slice position, first (1), interior (2), last (3)
        spos = 2
        if (islice==1) spos = 1
        if (islice==nslice) spos = 3

        ! Spectral slice of interest is [l,u)
        low = slice_cut(islice)
        upp = slice_cut(islice+1)
        
        ! Try different tuning???
        !wovlp = (upp - low) / 8.d0
        wovlp = (upp - low) / 10.d0
       
        poly_low = low - wovlp
        poly_upp = upp + wovlp
       
        ! Count theta eigenvalues in slice cut plus overlap 
        call count_values(poly_low,poly_upp,theta,spos,neigenpairs,k1,k2,nvec)
        nvec_ovlp = nvec
        write(std_out,*) '    after overlap', k1,k2

        ! How to compute optimal degree in slice with overlap:
        ! ********* control amplification ratio ********
        ! The convergence ratio r0/rN depends on the amplification
        ! ratios f(l)/f(l-w) and f(u)/f(u+w). Increasing w should
        ! improve the convergence ratio.

        ! Filter support [l-w,u+w) scaled to [-1,1)
        a = (poly_low-c)/r
        b = (poly_upp-c)/r

        if (islice==1) then
            ! uj,gub is the interval mapped to -1,1
            ! in this interval Chebyshev poly is bounded by 1
            ! remember uj,gub is the interval to ignore
            ndeg = 4
            ! FIXME is the scaling in [-1,1] OK?
            do while(1.d0/cheb_poly(low,ndeg,poly_upp,ecut)<ramp) 
                ndeg = ndeg + 1
            end do
        else
            ndeg = 4
            finL = 0.d0; finR = 0.d0; foutL = 1.d0; foutR = 1.d0
            do while ( (finL/foutL < ramp) .and. (finR/foutR < ramp) .and. (ndeg<ndeg_max) )
                finL  = bandpassIndicator_sca((low-c)/r,a,b,ndeg)
                foutL = bandpassIndicator_sca(a      ,a,b,ndeg)
                finR  = bandpassIndicator_sca((upp-c)/r,a,b,ndeg)
                foutR = bandpassIndicator_sca(b      ,a,b,ndeg)
                ndeg = ndeg + 1
            end do
        end if

        ! Plot filter
        if (plot_filter_) then
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
        end if

        ! Print slice interval info
        write(std_out,*) 'Without overlap=', low, upp
        write(std_out,*) '          width=', upp-low
        write(std_out,*) '      scaled to=', (low-c)/r,(upp-c)/r
        write(std_out,*) 'With    overlap=', poly_low, poly_upp
        write(std_out,*) '          width=', poly_upp-poly_low
        write(std_out,*) '      scaled to=', a,b
        write(std_out,*) '      nvec_ovlp=', nvec_ovlp
        write(std_out,*) '           nvec=', nvec
        write(std_out,*) '           ndeg=', ndeg
        write(std_out,*) ' '
        ! end print

        ! Compute last index in extended memory (without ovlp)
        j2 = j1 + nvec_ovlp - 1

        ! TODO fix notation according to this
        slice%slice_fcol(islice) = k1
        slice%slice_ncols(islice) = k2-k1+1 ! nvec_ovlp
        slice%slice_fcol_ext(islice) = j1
        slice%poly_degrees(islice) = ndeg
        slice%part_low_bounds(islice) = low
        slice%part_upp_bounds(islice) = upp
        slice%poly_low_bounds(islice) = poly_low
        slice%poly_upp_bounds(islice) = poly_upp
 
        ! Update starting index of next slice in extended
        j1 = j2 + 1
           
    end do

    ! Free temporary memory space
    if (allocated(slice_cut)) ABI_FREE(slice_cut)
    if (allocated(resid_cut)) ABI_FREE(resid_cut)
    if (allocated(theta)) ABI_FREE(theta)
    if (allocated(resid)) ABI_FREE(resid)
   
end subroutine slice_cutSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_mergeConverged
!! NAME
!! slice_mergeConverged
!! 
!! FUNCTION
!! Find converged eigenvalues in each slice using a criterion.
!! [^Not true. For the moment we brutally merge using interval limits]
!! Return range of first and last index to merge per slice, both 
!! computed from a residual criterion on eigenvalues.
!! slice_mergeConverged is done locally on slice
!! while on linalg representation
!! so it is after Rayleigh-Ritz and before the last transposition
!! 
!! SOURCE

subroutine slice_mergeConverged(sliceAll,idx_ovlp,idx_merge,merge_option)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t) , intent(inout) :: sliceAll
    integer          , intent(in   ) :: merge_option
    integer, pointer , intent(inout) :: idx_merge(:,:)
    integer, pointer , intent(in   ) :: idx_ovlp(:,:)

    ! Local variables-------------------------------    
    integer :: neigenpairs,nslice,islice,i1,i2,j1,j2
    integer :: iband,nband_slice,nband_ovlp
    integer :: k,k1,k2,nvec_count,spos
    integer :: nband_conv,ndeg,k_out
    real(dp) :: ramp,slow,supp,flow,fupp
    real(dp) :: glb,gub
    real(dp) :: c,r ! center, radius
    type(xgBlock_t) :: Eig0_all, Res0_all
    ! arrays
    integer, pointer :: pband(:)
    real(dp), pointer :: theta0(:,:), resid0(:,:)
    real(dp), pointer :: thetaN(:,:), residN(:,:)
    !real(dp), pointer :: theta0_(:), resid0_(:)
    !real(dp), pointer :: thetaN_(:), residN_(:)
    real(dp), allocatable :: resid0_slice(:)
    real(dp), allocatable :: residN_slice(:)
    real(dp), allocatable :: theta_slice(:)
 
    ! *********************************************************************

    ! Various variables
    nslice = sliceAll%nslice
    neigenpairs = sliceAll%neigenpairs
    nband_ovlp = sliceAll%nband_ovlp
    glb = sliceAll%glb
    gub = sliceAll%gub
    ramp = sliceAll%ramp
    r = (gub - glb)/2.d0
    c = (gub + glb)/2.d0
    
    ! Memory pointers
    Eig0_all = sliceAll%Eig0%self
    Res0_all = sliceAll%Res0%self
    pband => sliceAll%pband ! sorting theta0 in increasing order 
    
    ! Results could be complex (with null imaginary part), so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(Eig0_all,theta0,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(Res0_all,resid0,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(sliceAll%xgeigen_ovlp,thetaN,rows=1,cols=nband_ovlp)
    call xgBlock_reverseMap(sliceAll%xgresidu_ovlp,residN,rows=1,cols=nband_ovlp)

    ! Use index maps to recover bands associated to a slice
    nband_conv = 0
    do islice=1,nslice
        slow = sliceAll%sbound(islice,1) ! slice low
        supp = sliceAll%sbound(islice,2) ! slice upp
        flow = sliceAll%sbound(islice,3) ! filter low
        fupp = sliceAll%sbound(islice,4) ! filter upp
        ndeg = sliceAll%ndeg(islice)

        ! Indices in overlapping objects theta0,resid0
        i1 = sliceAll%idx(islice,1)
        i2 = sliceAll%idx(islice,2)
        nband_slice = i2 - i1 + 1

        ! Indices in overlap-free objects thetaN,residN
        j1 = idx_ovlp(islice,1)
        j2 = idx_ovlp(islice,2)

        ABI_MALLOC(theta_slice,(nband_slice))
        
        ABI_MALLOC(residN_slice,(nband_slice))
        ABI_MALLOC(resid0_slice,(nband_slice))

        ! Assumes row distribution (process has all nbands_slice in memory)
        theta_slice(1:nband_slice) = thetaN(1, j1:j2) ! after slicing
        residN_slice(1:nband_slice) = sqrt(residN(1,j1:j2))
        resid0_slice(1:nband_slice) = sqrt(resid0(1,pband(i1:i2)))

        ! Slice position, first (1), interior (2), last (3)
        spos = 2
        if (islice==1) spos = 1
        if (islice==nslice) spos = 3

        write(std_out,*) ''
        write(std_out,*) '-------------- /Mark/ Slice', islice
        write(std_out,*) '      nband_slice=', nband_slice
        write(std_out,*) ''

        ! Number of eigenvalues in slice interval
        ! This is the ones we keep
        call count_values(slow,supp,theta_slice,spos,nband_slice,k1,k2,nvec_count)
        write(std_out,*) 'Partition:', slow, supp
        write(std_out,*) 'number of eigenvalues in Partition =', nvec_count
        write(std_out,*) 'max residual rN       in Partition =', maxval(sqrt(residN_slice(k1:k2)))
        ! Print boundaries of this kept range
        write(std_out,*) 'min eigenvalue in Partition (kept) =', minval(theta_slice(k1:k2))
        write(std_out,*) 'max eigenvalue in Partition (kept) =', maxval(theta_slice(k1:k2))
 
        ! Store limited indices in overlap-free memory, attention shift k1,k2 by j1
        if (nband_conv+k2-k1+1>neigenpairs) then
            ! excess of eigenvalues, ignore last ones
            write(std_out,'(a,i0)') 'last one is ', k2
            k2 = k2 - (nband_conv+k2-k1+1-neigenpairs)
            write(std_out,'(a,i0)') 'due to excess, shift last one to ', k2
            idx_merge(islice,1) = j1 + k1 - 1
            idx_merge(islice,2) = j1 + k2 - 1
            nband_conv = nband_conv + k2 - k1 + 1
        else
            idx_merge(islice,1) = j1 + k1 - 1
            idx_merge(islice,2) = j1 + k2 - 1
            nband_conv = nband_conv + k2 - k1 + 1
        end if 
        write(std_out,*) ''

        ! Number of eigenvalues in filter support
        call count_values(flow,fupp,theta_slice,spos,nband_slice,k1,k2,nvec_count)
        write(std_out,*) 'Support:', flow, fupp
        write(std_out,*) 'number of eigenvalues in Support      =', nvec_count
        write(std_out,*) 'max residual rN       in Support      =', maxval(sqrt(residN_slice(k1:k2)))
        ! Print boundaries of this kept range
        write(std_out,*) 'min eigenvalue in Support (converged) =', minval(theta_slice(k1:k2))
        write(std_out,*) 'max eigenvalue in Support (converged) =', maxval(theta_slice(k1:k2))

        if (allocated(theta_slice)) ABI_FREE(theta_slice)
        if (allocated(residN_slice)) ABI_FREE(residN_slice)
        if (allocated(resid0_slice)) ABI_FREE(resid0_slice)

    end do
    
    ! Report missing or extra eigenvalues
    if (nband_conv<neigenpairs) then
        ABI_ERROR("Not enough converged eigenvalues")
    else if (nband_conv>neigenpairs) then
        ABI_ERROR("Too many converged eigenvalues")
    end if

    ! TODO 
    ! * count how many thetaN_ are in low,upp for every slice
    ! * count how many thetaN_ are outside current slice, and if they converged
    ! * count how many are in overlap region
    ! This will help diagnostic convergence "slice full"

    ! FIXME actually do the copy from extended to io
    ! replace blockCopy by xgBlock_copy
    call xgBlock_reshape(xgeigen, (/1,nband/))
    call xgBlock_reshape(xgresidu, (/1,nband/))
    j1 = 1                      ! start copy to overlapping mem (cg,eig,resid)
    do islice=1,nslice
        i1 = idx_merge(islice,1) ! start read from overlap-free mem
        i2 = idx_merge(islice,2)
        nband_merge = i2 - i1 + 1
        j2 = j1 + nband_merge - 1
        call slice_blockCopy(sliceAll%xgx0_ovlp,xgx0,i1,j1,i2,j2)
        call slice_blockCopy(sliceAll%xgeigen_ovlp,xgeigen,i1,j1,i2,j2)
        call slice_blockCopy(sliceAll%xgresidu_ovlp,xgresidu,i1,j1,i2,j2)
        j1 = j2 + 1
    end do
    ! Write clean as this:
    !do islice=1,nslice
    !    ncols = spsl%nband_slice(islice)
    !    fcol = spsl%fcol_slice_merge(islice)
    !    fcol_buf = spsl%fcol_buf_merge(islice)
    !    call xgBlock_setBlock(X0,xgcols_out,spacedim,ncols,fcol=fcol)
    !    call xgBlock_setBlock(spsl%Bufr,xgcols_in,spacedim,ncols,fcol=fcol_buf)
    !    call xgBlock_copy(xgcols_in,xgcols_out)
    !end do
    call xgBlock_reshape(xgeigen, (/nband,1/))
    call xgBlock_reshape(xgresidu, (/nband,1/))
 
end subroutine slice_mergeConverged
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_init
!! NAME
!! slice_init
!! 
!! SOURCE

subroutine slice_init(sliceAll,slice,islice)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t), intent(in   ) :: sliceAll
    type(slice_t   ), intent(inout) :: slice
    integer         , intent(in   ) :: islice

    ! Local variables-------------------------------    
    integer :: bandpp_slice,i1,i2,npband_slice,nband_slice
    integer :: tim_slice_init
    real(dp) :: tsec(2)

    ! *********************************************************************

    tim_slice_init = tim_slice2_init + islice - 1
    call timab(tim_slice_init,1,tsec)
    
    ! Input variables
    slice%islice = islice
    
    ! Variables inherited from sliceAll
    slice%glb = sliceAll%glb
    slice%gub = sliceAll%gub
    slice%low = sliceAll%sbound(islice,3)
    slice%upp = sliceAll%sbound(islice,4)
    slice%degree = sliceAll%ndeg(islice)
    npband_slice = sliceAll%npbandSlice(islice)
    i1 = sliceAll%idx(islice,1)
    i2 = sliceAll%idx(islice,2)

    ! Number of bands in slice
    nband_slice = i2 - i1 + 1
    slice%nband = nband_slice

    ! Number of bands per MPI process
    bandpp_slice = nband_slice
    if (sliceAll%paral_kgb==1) bandpp_slice=nband_slice/npband_slice
    slice%bandpp = bandpp_slice

    call slice_allocateAll(sliceAll,slice)

    ! Make dimension compatible with xgeigen,xgresidu of slicewf
    call xgBlock_reshape(slice%EW%self, (/nband_slice,1/))
    call xgBlock_reshape(slice%RW%self, (/nband_slice,1/))

    ! Attention T_n<0 for n odd
    ! TODO add a warning for fist slice

    ! Set pointers to slice memory
    slice%xgx0 = slice%XW%self
    slice%xgeigen = slice%EW%self
    slice%xgresidu = slice%RW%self
    slice%xgocc = slice%OCCW%self
    
    call timab(tim_slice_init,2,tsec)

end subroutine slice_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allocateAll
!! NAME
!! slice_allocateAll
!! 
!! FUNCTION
!! Allocate slice memory spaces on CPU/GPU
!! 
!! SOURCE

subroutine slice_allocateAll(sliceAll,slice)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t), intent(in   ) :: sliceAll
    type(slice_t)   , intent(inout) :: slice

    ! Local variables-------------------------------    
    integer :: spacedim,paral_kgb,nband,bandpp
    integer :: ndeg,space,eigenProblem,spacecom
    integer :: me_g0,me_g0_fft,comm_rows,comm_cols
    integer :: gpu_option,gpu_kokkos_nthrd,gpu_thread_limit
    integer :: nbdbuf, oracle
    integer  :: total_spacedim, ierr
    logical :: paw
    real(dp) :: tolerance,ecut
    real(dp) :: oracle_factor, oracle_min_occ

    ! *********************************************************************
    
    ! Various parameters common to all slices
    space            = sliceAll%space
    spacedim         = sliceAll%spacedim
    spacecom         = sliceAll%spacecom
    me_g0            = sliceAll%me_g0
    me_g0_fft        = sliceAll%me_g0_fft
    tolerance        = sliceAll%tolerance
    ecut             = sliceAll%ecut
    paral_kgb        = sliceAll%paral_kgb
    comm_rows        = sliceAll%comm_rows
    comm_cols        = sliceAll%comm_cols
    paw              = sliceAll%paw
    gpu_option       = sliceAll%gpu_option
    gpu_kokkos_nthrd = sliceAll%gpu_kokkos_nthrd
    gpu_thread_limit = sliceAll%gpu_thread_limit

    ! not used but passed to chebfi with deactivated values
    oracle           = 0
    nbdbuf           = 0
    oracle_factor    = 1.d0
    oracle_min_occ   = 0.d0

    ! Various parameters of current slice
    nband  = slice%nband
    bandpp = slice%bandpp
    ndeg   = slice%degree
  
    call slice_free(slice)

    ! With current def each xg is distributed along spacecom=plane-wave MPI distr
    ! TODO Every slice MPI proc has entire space
    ! FIXME EW,RW devrait être alloués exactement comme xgx0 et xgresidu in chebfiwf 
    !! FIXME this is not necessary because extended spaces does the pointer
    write(std_out,'(a,i0)') '-----> Allocating slice result memory nband_slice=', nband
    call xg_init(slice%XW,space,spacedim,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
        
    ! TODO each column MPI has all rows and bandpp columns
    ! The subcomm_band should be defined outside this function. So in 
    ! slice_init 
    call xg_init(slice%XW,space,total_spacedim,bandpp,subcomm_band,&
&                me_g0=me_g0,gpu_option=gpu_option)
    
    call xg_init(slice%EW,SPACE_R,1,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%RW,SPACE_R,1,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%OCCW,SPACE_R,nband,1,gpu_option=gpu_option)! fill with zero

    ! Get total number of rows that are distributed along MPI processes
    total_spacedim = spacedim
    if (paral_kgb == 1) then
        call xmpi_sum(total_spacedim,spacecom,ierr)
    end if
    ! TODO paral_slice values should be encoded
    if (paral_slice == 0) then
        ! the communicator of chebfi is the large one
    else if (paral_slice == 1) then
        ! the communicator of chebfi is the subcomm one
    end if

    ! Every slice MPI proc has a part (bandpp) of this space
    write(std_out,'(a,i0)') '-----> Allocating slice working memory bandpp=', bandpp
    call chebfi_init(slice%chebfi,nband,spacedim,tolerance,ecut,paral_kgb,bandpp,&
&                    ndeg,nbdbuf,space,1,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
&                    oracle,oracle_factor,oracle_min_occ,gpu_option,&
&                    gpu_kokkos_nthrd=gpu_kokkos_nthrd,gpu_thread_limit=gpu_thread_limit)

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_free
!! NAME
!! slice_free
!! 
!! SOURCE

subroutine slice_free(slice)

    implicit none

    type(slice_t), intent(inout) :: slice
    integer :: tim_slice_free
    real(dp) :: tsec(2)

    tim_slice_free = tim_slice2_free + slice%islice - 1
    call timab(tim_slice_free,1,tsec)
    
    call xg_free(slice%XW)
    call xg_free(slice%EW)
    call xg_free(slice%RW)
    call xg_free(slice%OCCW)
    call chebfi_free(slice%chebfi)
    
    call timab(tim_slice_free,2,tsec)

end subroutine slice_free
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

!!****f* m_slice/count_values
!! NAME
!! count_values
!! 
!! FUNCTION
!! Return number of consecutive theta values in interval [a,b)
!! and their index range (first and last indices).
!! Slice position (spos) takes into account extremal indices.
!! 
!! SOURCE

subroutine count_values(low,upp,theta,spos,nband,i1,i2,nvec)

    implicit none

    !Arguments ------------------------------------
    integer, intent(in) :: spos,nband
    integer, intent(inout) :: i1,i2,nvec
    real(dp), intent(in) :: low,upp
    real(dp), intent(in) :: theta(1:nband)
    !Local variables-------------------------------

! *********************************************************************

    i1 = maxloc(theta, dim=1, mask=(theta < low)) + 1
    i2 = maxloc(theta, dim=1, mask=(theta < upp))            
    if (spos == 1) i1 = 1
    if (spos == 3) i2 = nband
    nvec = i2 - i1 + 1
    nvec = i2 - i1 + 1

end subroutine count_values
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_blockCopy
!! NAME
!! slice_blockCopy
!! 
!! FUNCTION
!! Read from A deep copy to B, given column ranges.
!! Handles GPU / CPU data location by calling xgBlock_copy as well as
!! enter data map OpenMP ranges for execution on CPU or GPU.
!! 
!! SOURCE

subroutine slice_blockCopy(A_in,B_out,a1,b1,a2,b2)

    implicit none

    !Arguments ------------------------------------
    type(xgBlock_t)  , intent(in   ) :: A_in
    type(xgBlock_t)  , intent(inout) :: B_out
    integer          , intent(in   ) :: a1,b1 ! first columns
    integer, optional, intent(in   ) :: a2,b2 ! last columns
    
    !Local variables-------------------------------
    integer :: nrowsA,ncolsA
    integer :: nrowsB,ncolsB
    integer :: a2_,b2_
    logical :: on_host 
    integer :: gpua,gpub ! gpu_options
    integer :: ncolsA_block, ncolsB_block
    type(xgBlock_t) :: A_block, B_block
    real(dp) :: tsec(2)

! *********************************************************************

    call timab(tim_slice_Acopy,1,tsec)
    
    call xgBlock_getSize(A_in ,nrowsA,ncolsA)
    call xgBlock_getSize(B_out,nrowsB,ncolsB)

    if (nrowsA/=nrowsB) ABI_ERROR("A and B should have the same number of rows")

    ! DEBUG print gpu_option of in/out objects
    ! OMP query: are we on CPU or not?
    call xgBlock_get_gpu_option(A_in ,gpua)
    call xgBlock_get_gpu_option(B_out,gpub)
    write(std_out,'(a,i0,a,i0,a,i0)') 'Memcopy: gpu_optionA ',gpua,' gpu_optionB ',gpub
    if (gpua==ABI_GPU_OPENMP .or. gpub==ABI_GPU_OPENMP) then 
        on_host = xomp_is_initial_device()
        write(std_out,*) 'Memcopy: on_host', on_host
    end if
    !! IML 14/03 TODO decide whether we perform actions of copy from to gpu in here ..
    ! first detect cases of incompatiblity
    ! print on_host has no point because only the CPU prints anyway .. always true
    ! ..

    ! Get last index in block
    a2_ = ncolsA
    b2_ = ncolsB
    if (present(a2)) a2_ = a2
    if (present(b2)) b2_ = b2
    
    ! Number of block columns in range
    ncolsA_block = a2 - a1 + 1
    ncolsB_block = b2 - b1 + 1

    ! Deep copy from X to Y
    ! ============================
    ! xgBlock_copy will do FIXME also depends on on_host!! 
    ! if A CPU and B GPU then: copy from B GPU to B CPU, copy from A CPU to B CPU
    ! so the copy is performed on CPU if one of A and B is on CPU.
    call xgBlock_setBlock(A_in ,A_block,rows=nrowsA,cols=ncolsA_block,fcol=a1)
    call xgBlock_setBlock(B_out,B_block,rows=nrowsB,cols=ncolsB_block,fcol=b1)
    call xgBlock_copy(A_block,B_block)
    
    call timab(tim_slice_Acopy,2,tsec)

end subroutine slice_blockCopy
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_initExtended
!! NAME
!! slice_initExtended
!! 
!! FUNCTION
!! Allocate extended buffer by replicating overlapping slice data.
!! This creates a buffer without data overlap on which we can safely
!! read and write data avoiding concurrent memory access. The 
!! overlapping data has two independent copies for adjacent slices.
!! The extended memory space is bigger than original one.
!! 
!! SOURCE

subroutine slice_initExtended(slice,pband_ptr)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, pointer, intent(in) :: pband_ptr
    
    ! Local variables-------------------------------    
    integer :: me_g0,nband,comm
    integer :: nrows,ncols
    integer :: islice,nslice
    integer :: fcol,fcol_ext
    type(xgBlock_t) :: xgcols_in, xgcols_out 
    
    ! *********************************************************************

    nslice = slice%nslice
    nrows = slice%spacedim
    comm = slice%spacecom
    nband = slice%nband
    me_g0 = slice%me_g0
    ncol = sum(slice%nband_slice)

    ! Permute columns of X in linalg representation
    call xgBlock_permuteCols(slice%X,nrow,nband,pband_ptr)

    if (nslice==1) then
        slice%Xext_linalg = slice%X
        slice%Eext_linalg = slice%eigen
        slice%Rext_linalg = slice%residu
    else
        ! Allocate extended spaces in linalg representation (on CPU)
        call xg_init(slice%Xext_,space,nrow,ncol,comm,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
        call xg_init(slice%Eext_,SPACE_R,1,ncol,comm,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
        call xg_init(slice%Rext_,SPACE_R,1,ncol,comm,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
    
        ! Define pointers
        slice%Xext_linalg = slice%Xext_%self
        slice%Eext_linalg = slice%Eext_%self
        slice%Rext_linalg = slice%Rext_%self

        ! Copy X to extended by blocks
        ! Reminder: xgBlock_copy is always on CPU expect if both blocks are on GPU
        do islice=1,nslice
            ncols = slice%slice_ncols(islice)
            fcol = slice%slice_fcol(islice)
            fcol_ext = slice%slice_fcol_ext(islice)
            call xgBlock_setBlock(slice%X,xgcols_in,nrows,ncols,fcol=fcol)
            call xgBlock_setBlock(slice%Xext_linalg,xgcols_out,nrows,ncols,fcol=fcol_ext)
            call xgBlock_copy(xgcols_in,xgcols_out)
        end do
    end if

    ! Unitary test
    ABI_CHECK(cols(slice%Xext_linalg)==ncol,'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # cols of Xext_linalg ', xmpi_comm_rank(comm), cols(slice%Xext_linalg)

end subroutine slice_initExtended
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_freeExtended
!! NAME
!! slice_freeExtended
!! 
!! SOURCE

subroutine slice_freeExtended(slice)

    implicit none

    type(slice_t), intent(inout) :: slice
   
    if (slice%nslice>1) then
        call xg_free(slice%Xext_)
        call xg_free(slice%Eext_)
        call xg_free(slice%Rext_)
    end if

end subroutine slice_freeExtended
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_distributeExtended
!! NAME
!! slice_distributeExtended
!!
!! FUNCTION
!! Distribute Xext_linalg across *all* MPI processes and allocate 
!! its distributed version Xext on individual processes.
!! This routine applied transposition across *all* MPI processes.
!! After the transposition each process contains the correct
!! bandpp corresponding to the slice so that no additional communication
!! has to be performed in order to bring band slices to processes.
!!
!! SOURCE

subroutine slice_distributeExtended(slice,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: nspinor

    ! Local variables -------------------------------
    integer :: specedim,neigenpairs
    integer, pointer :: ncolsColsRows_ptr(:) => null()

    ! *********************************************************************
 
    spacedim = slice%spacedim
    neigenpairs = slice%neigenpairs

    if (chebfi%paral_kgb == 1) then

        nprocs = xmpi_comm_size(comm(X0))

        ! Rule for number of bands per process
        ncolsColsRows_ptr => slice%mpiData%ncolsColsRows

        ! Allocate slice%Xext according to the target MPI distribution for slices
        call xgTransposer_constructor(slice%xgTransposerX,slice%Xext_linalg,slice%Xext,nspinor,&
            STATE_LINALG,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0_fft,&
            gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
            custom_ncolsColsRows=.true.,ncolsColsRows_sub=ncolsColsRows_ptr)
   
        slice%xgTransposerX%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
   
        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE_XEXPAND)
        call xgTransposer_transpose(slice%xgTransposerX,STATE_COLSROWS)
        ABI_NVTX_END_RANGE()

    else
        call xgBlock_setBlock(slice%Xext, slice%Xext_linalg, spacedim, neigenpairs)
    end if

    ! Unitary test
    ABI_CHECK(cols(slice%Xext)==ncolsColsRows(xmpi_comm_rank(spacecom)),'wrong colsrows representation')
    write(*,'(a,i6,i6)') '# proc has # cols of Xext ', xmpi_comm_rank(spacecom), cols(slice%Xext)
    ! cols(slice%X_expand) should be equal to ncolsColsRows(i) where i is the rank of MPI process
    
end subroutine slice_distributeExtended
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_run
!! NAME
!! slice_run
!! 
!! FUNCTION
!! Diagonalize slices in parallel. Input/output is the 
!! extended buffer in colsrows representation.
!! Notice that initial objects xgx0,eigen,residu are not input
!!
!! IML TODO add count
!! In early SCF iterations, count the number of Ritz values that fall within the
!! perturbed spectral interval of each slice, where the size of the perturbation
!! is related to the residual norm of each Ritz pair. We could therefore terminate
!! the subspace iterations when the counts no longer change.
!!
!! SOURCE

subroutine slice_run(slice,getAX_BX,getBm1X,nspinor)

    implicit none

    !Arguments ------------------------------------    
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ) :: nspinor
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

    !Local variables-------------------------------
    integer :: tot_nrows,ncols_slice,spacecom
    integer :: color,key,slice_comm,ierr
    integer, target, allocatable :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null()
    type(xgBlock_t) :: X0

    ! *********************************************************************

    call xgBlock_getSize(slice%X_expand,tot_nrows,ncols_slice)
    spacecom = comm(xgBlock_colsrows) ! global communicator
    ! spacedim is smaller than tot_nrows

    ! slice colsrows workspace points to column range of distributed X_expand
    call xgBlock_setBlock(slice%X_expand,X0,rows=tot_nrows,cols=ncols_slice)
    slice%X = X0

    ! if gpu: This is important! Because slice%Xext is on CPU
    ! call xgBlock_copy_to_gpu(slice%X)
    ! call xgBlock_set_gpu_option(slice%X)
    
    ! Synchronize before creating a subcomm
    call xmpi_barrier(schedule%comm_global)
    
    ! Split global communicator so that only procs with the same color communicate
    key = sliceTasks_queryRank(env)
    color = sliceTasks_queryTask(env)
    call xmpi_comm_split(schedule%comm_global,color,key,slice_comm,ierr)

    ! Restrict all communications to slice subcommunicator
    call xgBlock_setComm(slice%X,slice_comm) ! colsrows representation
    call xgBlock_setComm(slice%X_linalg,slice_comm) ! linalg representation

    ! Compute distribution of rows (plane waves) in the subcommunicator
    call mpiSlice_distributeRows(mpi_slice,slice_comm)
    nrowsLinalg_ptr => mpi_slice%lookup
    nproc = xmpi_comm_size(slice_comm)
    ABI_MALLOC(nrowsLinalg, (nproc)) ! FIXME size known before slice construction. Move out
    nrowsLinalg_ptr => nrowsLinalg
    call compute_uniform_distribution(nrowsLinalg_ptr,nproc,tot_nrows)

    ! Create slice transposer using subcommunicator and allocate slice%X_linalg
    call xgTransposer_constructor(slice%xgTransposerX,slice%X_linalg,slice%X,nspinor,&
        STATE_COLSROWS,TRANS_ALL2ALL,chebfi%comm_rows,slice_comm,0,0,chebfi%me_g0_fft,&
        gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
        custom_ncolsColsRows=.true.,nrowsLinalg_sub=nrowsLinalg_ptr)
        ! true to allow different bandpp (avoid pad)

    ncols_slice = cols(slice%X)

    ! Do the same for AX and BX ..
    call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,&
        chebfi%AX%self,chebfi%xAXColsRows,STATE_COLSROWS)
    call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,&
        chebfi%BX%self,chebfi%xBXColsRows,STATE_COLSROWS)

    chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
    chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
    chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

    ! Body of computation
    ! ===============
    
    ! Apply polynomial filtering (requires colsrows state)
    if (islice==0) then
        !call slice_applyLowpassFilter(slice,getAX_BX,getBm1X,nspinor)
        call slice_applyLowpassFilter(chebfi,getAX_BX,getBm1X,nspinor)
    else
        !call slice_applyBandpassFilter(slice,getAX_BX,getBm1X,nspinor)
        call slice_applyBandpassFilter(chebfi,getAX_BX,getBm1X,nspinor)
    end if

    ! Transpose to linalg state
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then

        ! All MPI columns wait to finish
        call xmpi_barrier(chebfi%spacecom)

        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_LINALG)
        call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG)
        call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_LINALG)

        !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
        if (xmpi_comm_size(chebfi%spacecom) == 1) then 
            call xgBlock_setBlock(chebfi%xXColsRows,  chebfi%X,       spacedim, neigenpairs)
            call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, spacedim, neigenpairs)
            call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(chebfi%xXColsRows,  chebfi%X,       spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, spacedim, neigenpairs)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_transpose,2,tsec)

    ! Unitary test
    ABI_CHECK(rows(slice%X_linalg)==nrowsLinalg(xmpi_comm_rank(spacecom)),'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # rows of slice X ', xmpi_comm_rank(spacecom), rows(slice%X_linalg)

    ! Perform Rayleigh-Ritz and compute residuals (requires linalg state)
    call slice_RayleighRitz(slice,eigen,residu)
    
    ! Transpose to colsrows state (X only)
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then

        ! All MPI rows wait to finish
        call xmpi_barrier(chebfi%spacecom)

        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_COLSROWS)

        !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
        if (xmpi_comm_size(chebfi%spacecom) == 1) then 
            call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_transpose,2,tsec)

    ! TODO Deal with xgeigen and xgresidu
    !call xgBlock_reshape(slice%xgeigen, (/1,nband_slice/))
    !call xgBlock_reshape(slice%xgresidu, (/1,nband_slice/))
    !call xgBlock_copy_from_gpu(slice%xgeigen)
    !call xgBlock_copy_from_gpu(slice%xgresidu)
    !call slice_blockCopy(slice%xgeigen,sliceAll%xgeigen_ovlp,1,j1,nband_slice,j2)
    !call slice_blockCopy(slice%xgresidu,sliceAll%xgresidu_ovlp,1,j1,nband_slice,j2)

    ! Unitary test
    ABI_CHECK(cols(slice%X)==ncols_slice,'wrong colsrows representation')
    write(*,'(a,i6,i6)') '# proc has # cols of slice X ', xmpi_comm_rank(spacecom), cols(slice%X)

    ! Copy slice solution to the extended buffer (requires colsrows state)
    call xgBlock_copy(slice%X,X0)
    ! FIXME same for eigen, residu?

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    if (gpu_option==ABI_GPU_KOKKOS) then
        call gpu_device_synchronize()
    end if
#endif

    ! Free transposer objects
    if (chebfi%paral_kgb == 1) then
        call xgTransposer_free(chebfi%xgTransposerX)
        call xgTransposer_free(chebfi%xgTransposerAX)
        call xgTransposer_free(chebfi%xgTransposerBX)
    end if
    
    ! Free memory
    if (allocated(nrowsLinalg)) ABI_FREE(nrowsLinalg)
    call xgTransposer_free(slice%xgTransposerX)

end subroutine slice_run
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_applyLowpassFilter
!! NAME
!! slice_applyLowpassFilter
!!
!! FUNCTION
!! Apply Lowpass filter using Chebyshev polynomial on a set of vectors.
!!
!! INPUTS
!!  slice   = spectral slice parameters
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm 
!!  on a single spectral slice
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine slice_applyLowpassFilter(slice,getAX_BX,getBm1X,nspinors)

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
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
    !Local variables-------------------------------
    type(chebfi_t) :: chebfi

    ! *********************************************************************

    chebfi = slice%chebfi

    ! When entering DivResult%self should contain all eigenvalues in slice
    ! in the form (bandpp,1) as a ROW vector
    comm_slice = slice%mpiData%comm_sub

    ! This is the maximum and the minimum eigenvalue in the slice
    if (slice%paral_kgb == 1) then
        call xmpi_max(slice%maxeig,maxeig_global,comm_slice,ierr)
        call xmpi_min(slice%mineig,mineig_global,comm_slice,ierr)
    else
        call xmpi_max(slice%maxeig,maxeig_global,comm,ierr)
        call xmpi_min(slice%mineig,mineig_global,comm,ierr)
    end if

    eigenvalues = DivResults%self !! ....
    ! DivResults%self is already constructed from previous routines
    ! this routine should also set the maxeig_global, mineig_global then

    if (chebfi%paral_kgb == 0) then
        ABI_MALLOC(ndeg_filter_bands,(neigenpairs))
    else
        ABI_MALLOC(ndeg_filter_bands,(bandpp))
    end if
    ndeg_filter_bands(:) = ndeg_filter
    
    ! Spectral interval to amplify is [-oo, lambda_minus)
    lambda_minus = maxeig_global
    lambda_plus = slice%%ecut

    center = (lambda_plus + lambda_minus)*0.5
    radius = (lambda_plus - lambda_minus)*0.5

    one_over_r = 1/radius
    two_over_r = 2/radius

    !A * Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
    do ideg = 0, ndeg_filter - 1

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, bandpp)
        end if
        ABI_NVTX_END_RANGE()

        !A * Psi
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)

    end do ! ideg
    ABI_NVTX_END_RANGE()

    ! Scale X,AX,BX by amplification factor to reduce large values
    call chebfi_ampfactor(chebfi, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)

    call xg_free(DivResults) ! en fait ne pas faire ça ici
    ! car on garde DivResults à l'extérieur des slices aussi pour
    ! comparer les convergences
    ABI_FREE(ndeg_filter_bands)

end subroutine slice_applyLowpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_applyBandpassFilter
!! NAME
!! slice_applyBandpassFilter
!!
!! FUNCTION
!! Apply Bandpass filter using Chebyshev-Jackson polynomial on a set of vectors.
!!
!! INPUTS
!!  slice   = spectral slice parameters
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm 
!!  on a single spectral slice
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine slice_applyBandpassFilter(slice,getAX_BX,getBm1X,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
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

    ! Local variables-------------------------------
    type(chebfi_t) :: chebfi
    integer :: space
    integer :: spacedim
    integer :: neigenpairs
    integer :: nline
    integer :: gpu_option
    integer :: nrows, ncols
    integer :: iline, ilinep1, iband, ierr
    real(dp) :: tolerance
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: center, radius
    real(dp) :: ck, mu, damp, tau  ! bandpass filter parameters
    real(dp) :: alow,bupp,low,upp  ! slice interval 
    type(xg_t) :: ChebyExpansion   ! Chebyshev expansion for vectors

    ! *********************************************************************

    chebfi = slice%chebfi

    ! Initialize solutions for slice using chebfi 
    space = chebfi%space
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs 
    tolerance = chebfi%tolerance
    gpu_option = chebfi%gpu_option
    chebfi%eigenvalues = eigen
    chebfi%X = X0
    nrows = spacedim
    ncols = neigenpairs
    if (chebfi%paral_kgb==1) then
        nrows = chebfi%total_spacedim
        ncols = chebfi%bandpp
    end if
    
    ! Global spectral interval
    radius = (slice%gub - slice%glb)/2.d0   ! entire spectrum radius
    center = (slice%gub + slice%glb)/2.d0   ! entire spectrum center
    one_over_r = 1/radius
    two_over_r = 2/radius
    ! Target local to amplify scaled to [-1,1)
    nline = slice%degree                    ! polynomial filter degree
    low = slice%low                         ! filter support low bound
    upp = slice%upp                         ! filter support upper bound
    alow = (low - center) / radius          ! scaled filter support low
    bupp = (upp - center) / radius          ! scaled filter support upp

    ! AX_next=A*X -> 1 Hamiltonian application
    call chebfi_getAX_BX(chebfi, getAX_BX)

    ! B-orthonormalize X, BX and AX
    !call xg_Borthonormalize(chebfi%xXColsRows,chebfi%xBxColsRows,ierr,1,gpu_option,AX=chebfi%xAXColsRows)
    ! IL TODO Deflate vectors 10/03/2025

    ! Why do this here?
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    write(std_out,*) 'TRACE initialize Chebyshev expansion (hopefuly on GPU)'
    ! Compute Chebyshev polynomial expansion on X iteratively on iline=0,nline
    ! Initialize Xsum = 0 (bands are distributed)
    ABI_NVTX_START_RANGE(NVTX_SLICE_EXPANSION)
    call xg_init(ChebyExpansion, chebfi%space, nrows, ncols, chebfi%spacecom, gpu_option=gpu_option)
    call xgBlock_zero(ChebyExpansion%self)
    ! X_next=X -> iline=0 Hamiltonian applications
    ck = Pi/(nline+2)
    mu = 1/Pi*(ACOS(alow)-ACOS(bupp))
    damp = 1.d0
    !Xsum = mu(0)*damp(0)*X_next + Xsum
    call xgBlock_saxpy(ChebyExpansion%self, mu*damp, chebfi%xXColsRows)
    ABI_NVTX_END_RANGE()

    write(std_out,*) 'TRACE start Slice core'
    ABI_NVTX_START_RANGE(NVTX_BANDPASS_CORE)
    do iline = 0, nline - 1  

        ! X_next=2/r*(AX_next-c*X_next)-X_prev, -> iline+1 Hamiltonian applications
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ! xXColsRows=X_next
        call chebfi_swapInnerBuffers(chebfi, nrows, ncols)

        ! Add new term to the Chebyshev expansion
        !Xsum = damp(i+1)*mu(i+1)*X_next + Xsum
        ABI_NVTX_START_RANGE(NVTX_CHEBFI_EXPANSION)
        ilinep1 = iline + 1
        mu = 2/Pi * (SIN(ilinep1*ACOS(alow)) - SIN(ilinep1*ACOS(bupp)))/ilinep1
        damp = ((1 - ilinep1/(nline+2))*SIN(ck)*COS(ilinep1*ck) + 1/(nline+2)*COS(ck)*SIN(ilinep1*ck))/SIN(ck)
        call xgBlock_saxpy(ChebyExpansion%self, mu*damp, chebfi%xXColsRows)
   
        ! Store term before exit
        ! AX_next=A*X_next -> iline+2 Hamiltonian applications
        if (iline==nline-1) then
            ! X_next=Xsum (copy Xsum to X_next)
            call xgBlock_copy(ChebyExpansion%self, chebfi%xXColsRows)
        end if
        ABI_NVTX_END_RANGE()

        ! Apply A and B to X
        call chebfi_getAX_BX(chebfi, getAX_BX)
    
    end do ! end iline
    ABI_NVTX_END_RANGE()

    ! All slice processes wait for filter done
    ! FIXME Is this necessary? No communication happens actually
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Free Chebyshev expansion workspace
    call xg_free(ChebyExpansion)

end subroutine slice_applyBandpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_RayleighRitz
!! NAME
!! slice_RayleighRitz
!!
!! FUNCTION
!! Wrapper for xg_RayleighRitz + compute residual.
!! Assumes (X,AX,BX) have MPI row distribution.
!!
!! SOURCE

subroutine slice_RayleighRitz(slice,eigen,residu)

    ! Arguments ***
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: residu
    ! Local variables ***
    type(chebfi_t) :: chebfi
    integer :: ierr

    chebfi = slice%chebfi
    chebfi%eigenvalues = eigen

    ! Apply Rayleigh-Ritz for each MPI row
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
    call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,chebfi%eigenvalues,ierr,0,tim_RR,&
        chebfi%gpu_option,solve_ax_bx=.true.)
    ABI_NVTX_END_RANGE()
    
    if ( ierr /= 0 ) then
        ABI_WARNING("RayleighRitz did not work, but continue anyway.")
    end if

    ! Compute residual for each MPI row and store it to AX
    if (chebfi%paw) then
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
    else
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
    end if

    ! Wait until all MPI rows have computed their residual
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Communicate MPI rows to compute residual norm squared
    call xgBlock_colwiseNorm2(chebfi%AX%self,residu)

end subroutine slice_RayleighRitz
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceTasks_init
!! NAME
!! sliceTasks_init
!! 
!! FUNCTION
!! Set parameters of slice parallelisation environment.
!! Assumes that spectral slices have already been splitted
!! (uses degree). Note that schedule is target to allow pointers
!! targeting its member variables (smart!).
!! 
!! SOURCE

subroutine sliceTasks_init(schedule,slice_sizes,slice_degrees,nrows_tot,spacecom,&
    resource_allocation)

    implicit none

    ! Arguments
    type(sliceTasks_t), target, intent(inout) :: schedule
    integer, pointer, intent(in) :: slice_sizes(:)
    integer, pointer, intent(in) :: slice_degrees(:)
    integer, intent(in) :: nrows_tot
    integer, intent(in) :: spacecom
    integer, intent(in) :: resource_allocation

    ! Local variables
    integer :: ntasks, nprocs, iproc, task_me
    ! Arrays
    integer, allocatable, target :: weights(:)
    integer, pointer :: weights_ptr(:) => null()
    integer, pointer :: task_ncols_ptr(:) => null()
    integer, pointer :: taks_nprocs_ptr(:) => null() 
    integer, pointer :: assigned_task_ptr(:) => null()
    integer, pointer :: blockcols_me_ptr(:) => null()
    integer, pointer :: blocrows_me_ptr(:) => null()
    
! *********************************************************************

    call sliceTasks_free(env)

    ntasks = size(slice_sizes)
    nprocs = xmpi_comm_size(spacecom)

    schedule%ntasks = ntasks
    schedule%nprocs = nprocs
    schedule%comm_global = spacecom
    schedule%nrows_tot = nrows_tot
    schedule%resource_allocation = resource_allocation

    ! Initialize MPI distribution flags to linalg (no transpose yet)
    schedule%use_blockrows = .true.
    schedule%use_blockcols = .false.

    ! Set CPU/GPU flags
    call sliceTasks_queryTarget(schedule)

    if(.not.allocated(schedule%task_ncols)) ABI_MALLOC(schedule%task_ncols, (ntasks))
    if(.not.allocated(schedule%task_nprocs)) ABI_MALLOC(schedule%task_nprocs, (ntasks))
    if(.not.allocated(schedule%assigned_task)) ABI_MALLOC(schedule%assigned_task, (nprocs))
    if(.not.allocated(weights)) ABI_MALLOC(weights, (ntasks))

    schedule%task_ncols(:) = slice_sizes(:)
    task_ncols_ptr => schedule%task_ncols
    task_nprocs_ptr => schedule%task_nprocs
    assigned_task_ptr => schedule%assigned_task
    weights_ptr => weights

    ! Apply weighted fair allocation with various balance criteria for load balance
    select case(resource_allocation)
    case(FAIR_BANDPP)
        weights(:) = 1 
    case(FAIR_BANDPP_WDEG)
        weights(:) = slice_degrees(:)
    end select
    
    ! Solve allocation problem to find the amount of resource allocated to each slice
    call fair_allocation(ntasks, slice_sizes_ptr, weights_ptr, nprocs, task_nprocs_ptr)
    
    ! Call the subroutine to assign tasks (=slices) to processes
    call assign_tasks_to_processes(task_nprocs_ptr, assigned_task_ptr)
    do iproc = 1, ntasks
        write(*,'(a,i5,a,i5)') "Process ", itask, " is in task ", schedule%assigned_task(i)
    end do

    ! Deduce band capacity per process (=bandpp) on individual slice
    ! a process only has the slice that corresponds to it
    task_me = sliceTasks_queryTask(schedule)
    ncols_me = schedule%task_ncols(task_me)
    nprocs_me = schedule%task_nprocs(task_me)

    schedule%task_me = task_me
    schedule%ncols_me = ncols_me
    schedule%nprocs_me = nprocs_me

    ! Series of allocations corresponding to _me
    if (.not.allocated(schedule%blockcols_me)) ABI_MALLOC(schedule%blockcols_me,(nprocs_me))
    if (.not.allocated(schedule%blockrows_me)) ABI_MALLOC(schedule%blockrows_me,(nprocs_me))

    blockcols_me_ptr => schedule%blockcols_me
    blockrows_me_ptr => schedule%blockrows_me

    ! Compute size of column-blocks in MPI col distribution
    call distribute_vectors(nprocme,ncols_me,blockcols_me_ptr)
    
    ! Compute size of row-blocks in MPI row distribution 
    call distribute_vectors(nprocs_me,nrows_tot,blockrows_me_ptr)

    ! Free temporary memory
    if (allocated(weights)) ABI_FREE(weights) 

end subroutine sliceTasks_init
!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceTasks_free
!! NAME
!! sliceTasks_free

subroutine sliceTasks_free(schedule)

    implicit none 

    type(sliceTasks_t), intent(inout) :: schedule

    if(allocated(schedule%task_ncols)) ABI_FREE(schedule%task_ncols)
    if(allocated(schedule%task_nprocs)) ABI_FREE(schedule%task_nprocs)
    if(allocated(schedule%assigned_task)) ABI_FREE(schedule%assigned_task)
    if(allocated(schedule%blockcols_me)) ABI_FREE(schedule%blockcols_me)   
    if(allocated(schedule%blocrows_me)) ABI_FREE(schedule%blockrows_me)

end subroutine sliceTasks_free
!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceTasks_queryRank
!! NAME
!! sliceTasks_queryRank

integer function sliceTasks_queryRank(schedule) result(rank_me)

    implicit none
    type(sliceTasks_t), intent(in) :: schedule

    rank_me = xmpi_comm_rank(schedule%comm_global)

end function sliceTasks_queryRank
!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceTasks_queryTask
!! NAME
!! sliceTasks_queryTask

integer function sliceTasks_queryTask(schedule) result(task_me)

    type(sliceTasks_t), intent(in) :: schedule
    integer :: rank_me
    
    rank_me = xmpi_comm_rank(schedule%comm_global)
    task_me = schedule%assigned_task(rank_me+1)

end function sliceTasks_queryTask
!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceTasks_queryTarget
!! NAME
!! sliceTasks_queryTarget
!! 
!! FUNCTION
!! Query OpenMP offloading API to set/get GPU/CPU flags
!! 
!! SOURCE

subroutine sliceTasks_queryTarget(schedule,on_host,on_device)

    implicit none
    type(sliceTasks_t), intent(inout) :: schedule
    logical, optional, intent(inout) :: on_host
    logical, optional, intent(inout) :: on_device
    logical :: on_host_
    logical :: on_device_
    integer :: device_id

    device_id = xomp_get_device_num()
    on_host_ = (device_id < 1) ! outside target (-1), inside target host (0)
    on_device_ = (device_id > 0) ! inside target not host (device number)

    schedule%on_host = on_host_
    schedule%on_device = on_device_

    if (present(on_host)) on_host = on_host_
    if (present(on_device)) on_device = on_device_

end subroutine sliceTasks_queryTarget
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

end module m_slice
!!***



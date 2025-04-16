!!****f* ABINIT/m_polyfi
!! NAME
!! m_polyfi
!!
!! FUNCTION
!! This module contains the types and routines used to apply the Polynomial Filtering method, 
!! extending Chebyshev filtering ('chebfi') to arbitrary spectral interval to be amplified and 
!! new types of bandpass polynomials. It mainly defines 'polyfi' datatypes and associated methods. 
!!
!! COPYRIGHT
!! Copyright (C) 2025- ABINIT group (IL)
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

module m_polyfi

    use defs_basis
    use defs_abitypes
    use m_abicore
    use m_errors
    use m_time, only : timab
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

    !Several (private) parameters
    !-------------------------------------------------

    integer, parameter :: tim_transpose    = 1758
    integer, parameter :: tim_RR           = 1757

    ! Public 'polyfi' datatype
    !-------------------------------------------------
    type, public :: polyfi_t 

        ! Memory workspace
        type(chebfi_t) :: chebfi  

        ! Spectral interval bounds
        real(dp) :: mineig_global        ! lower spectral bound (will be scaled to -1)
        real(dp) :: maxeig_global        ! upper spectral bound (will be scaled to 1)
        real(dp) :: lambda_minus         ! lower bound of interval to amplify
        real(dp) :: lambda_plus          ! upper bound of interval to amplify

        ! Flags
        logical :: is_lowpass = .true.   ! use lowpass polynomial (Chebyshev polynomial)
        logical :: is_bandpass = .false. ! use bandpass polynomial (Heaviside expanded on Chebyshev)

    end type polyfi_t

    ! Public methods associated to 'polyfi' datatype
    !-------------------------------------------------
    public :: polyfi_init
    public :: polyfi_run
    public :: polyfi_free

    CONTAINS  
!=====================================================================
!!***

!!****f* m_polyfi/polyfi_init
!! NAME
!! polyfi_init
!! 
!! FUNCTION
!! Initialize a 'polyfi' datastructure.
!!
!! INPUTS
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
!! OUTPUT
!!
!! SIDE EFFECTS
!!  polyfi <type(polyfi_t)>=all data used to apply Polynomial Filtering algorithm
!!
!! SOURCE

 subroutine polyfi_init(polyfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp,&
         ndeg_filter,mineig_global,maxeig_global,lambda_minus,lambda_plus,is_lowpass,&
         space,eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
         gpu_option,gpu_kokkos_nthrd,gpu_thread_limit)

     implicit none

    ! Arguments ------------------------------------
    integer       , intent(in   ) :: bandpp
    integer       , intent(in   ) :: eigenProblem
    integer       , intent(in   ) :: me_g0
    integer       , intent(in   ) :: me_g0_fft
    integer       , intent(in   ) :: neigenpairs
    integer       , intent(in   ) :: ndeg_filter
    integer       , intent(in   ) :: comm_cols
    integer       , intent(in   ) :: comm_rows
    integer       , intent(in   ) :: paral_kgb
    integer       , intent(in   ) :: space
    integer       , intent(in   ) :: spacecom
    integer       , intent(in   ) :: spacedim
    integer       , intent(in   ) :: gpu_option
    logical       , intent(in   ) :: paw
    logical       , intent(in   ) :: is_lowpass
    real(dp)      , intent(in   ) :: mineig_global
    real(dp)      , intent(in   ) :: maxeig_global
    real(dp)      , intent(in   ) :: lambda_minus
    real(dp)      , intent(in   ) :: lambda_plus
    real(dp)      , intent(in   ) :: ecut
    real(dp)      , intent(in   ) :: tolerance
    type(polyfi_t), intent(inout) :: polyfi
    integer       , intent(in   ), optional :: gpu_kokkos_nthrd
    integer       , intent(in   ), optional :: gpu_thread_limit

    ! Local variables-------------------------------
    integer :: nbdbuf
    integer :: oracle
    real(dp) :: oracle_factor
    real(dp) :: oracle_min_occ

    ! *********************************************************************

    polyfi%mineig_global = mineig_global
    polyfi%maxeig_global = maxeig_global
    polyfi%lambda_minus = lambda_minus
    polyfi%lambda_plus = lambda_plus
    polyfi%is_lowpass = is_lowpass
    polyfi%is_bandpass = (.not. is_lowpass)
    
    ! Oracle not used but passed to chebfi with deactivated values
    oracle = 0
    nbdbuf = 0
    oracle_factor = 1.d0
    oracle_min_occ = 0.d0

    call chebfi_free(polyfi%chebfi)

    ! Define chebfi object from Colsrows representation
    call chebfi_init(polyfi%chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp,&
        ndeg_filter,nbdbuf,space,1,comm,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
        oracle,oracle_factor,oracle_min_occ,gpu_option,gpu_kokkos_nthrd=gpu_kokkos_nthrd,&
        gpu_thread_limit=gpu_thread_limit,from_linalg=.false.)

 end subroutine polyfi_init
 !!***
 
!----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_free
!! NAME
!! polyfi_free
!! 
!! SOURCE

subroutine polyfi_free(polyfi)

    implicit none

    type(polyfi_t), intent(inout) :: polyfi

    call chebfi_free(polyfi%chebfi)

end subroutine polyfi_free
!!***

!----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_run
!! NAME
!! polyfi_run
!! 
!! FUNCTION
!! Apply Polynomial filtering to amplify spectral interval [lambda_minus,lambda_plus)
!! from a given set of vectors, using computational resources reserved by spacecom only.
!!
!! INPUT
!! polyfi <type(polyfi_t)>= all data used to apply Polynomial Filtering algorithm
!! X0= eigenvector guess distributed in ColsRows representation
!! getAX_BX= pointer to the function giving A|X> and B|X>
!!           A is typically the Hamiltonian H, and B the overlap operator S
!! getBm1X= pointer to the function giving B^-1|X>
!!          B is typically the overlap operator S
!! eigen= Rayleigh quotients associated to X0
!! residu= empty array
!! nspinor= number of spinorial components of the wavefunctions
!! 
!! SIDE EFFECTS
!! polyfi= workspaces used
!! X0= full set of vectors distributed in ColsRows representation 
!! eigen= full eigenvalues
!! residu= residuals, i.e. norm of (A-lambdaB)|X>
!! 
!! SOURCE


!! Genralization of chebfi_run() to be able to amplify any given spectral interval
!! using a lowpass OR a bandpass filter.

subroutine chebfi_runSlice(chebfi,X0,getAX_BX,getBm1X,eigen,residu,mineig_global,maxeig_global,&
        lambda_minus,lambda_plus,is_lowpass,nspinor)

subroutine polyfi_run(polyfi, X0, getAX_BX, getBm1X, eigen, residu, nspinor)

    implicit none

    !Arguments ------------------------------------    
    type(polyfi_t), intent(inout) :: polyfi
    type(xgBlock_t), intent(inout) :: X0
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: residu
    integer       , intent(in   ) :: nspinor
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
    integer :: color,my_rank,ierr
    integer, target, allocatable :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null()
    type(xgBlock_t) :: X0
    type(chebfi_t) :: chebfi

    ! *********************************************************************
 
    ! Memory workspace Is Chebfi
    chebfi = polyfi%chebfi

    ! Parameters
    total_spacedim = chebfi%total_spacedim
    neigenpairs = chebfi%neigenpairs
    spacecom = chebfi%spacecom
    num_proc = xmpi_comm_size(spacecom)

    ! 
    chebfi%xXColsRows = X0
    chebfi%eigenvalues = eigen

    ! thing is
    !! eigen%rows == AX%cols 
    ! this eigen has neigenpairs rows!!!!

    ! Arrays
    if(.not.allocated(nrowsLinalg)) ABI_MALLOC(nrowsLinalg,(num_proc))
    nrowsLinalg_ptr => nrowsLinalg
    nrowsLinalg(:) = polyfi%nrows_blockrows

    ! Compute row distribution across active resources
    ! TODO decide either we do this here OR we read it from polyfi members
    ! call distribute_vectors(num_proc, chebfi%total_spacedim, nrowsLinalg_ptr)

    ! Allocate chebfi%X in Linalg representation
    call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nspinor,&
        STATE_COLSROWS,TRANS_ALL2ALL,chebfi%comm_rows,spacecom,0,0,chebfi%me_g0_fft,&
        gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
        custom_ncolsColsRows=.true.,nrowsLinalg_sub=nrowsLinalg_ptr)
        ! true to allow different bandpp (avoid pad)

    call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,&
        chebfi%AX%self,chebfi%xAXColsRows,STATE_COLSROWS)
    call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,&
        chebfi%BX%self,chebfi%xBXColsRows,STATE_COLSROWS)

    chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
    chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
    chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

    ! Body of computation
    ! ===============
    
    ! Apply polynomial filtering (uses ColsRows representation)
    if (polyfi%is_lowpass) then

        ! Note: mineig_global and maxeig_global are not used more in chebfi
        !lambda_minus = maxeig_global
        !lambda_plus = chebfi%ecut

        call chebfi_lowpassFilter(polyfi%chebfi,eigen,lambda_minus,lambda_plus,getAX_BX,getBm1X)

    else if (polyfi%is_bandpass) then

        !lambda_minus = slice%lb(islice) = polyfi%lambda_minus
        !lambda_plus = slice%ub(islice) = polyfi%lambda_plus
        !mineig_global = guaranteed_lb = polyfi%mineig_global
        !maxeig_global = chebfi%ecut

        call chebfi_bandpassFilter(polyfi%chebfi,lambda_minus,lambda_plus,mineig_global,&
            maxeig_global,getAX_BX,getBm1X)

    end if

    ! Transpose to linalg state
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then

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
    ABI_CHECK(rows(chebfi%X)==nrowsLinalg(xmpi_comm_rank(spacecom)),'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # rows of slice X ', xmpi_comm_rank(spacecom), rows(chebfi%X)

    ! Apply Rayleigh-Ritz for each MPI row
    ABI_NVTX_START_RANGE(NVTX_POLYFI_RR)
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

    ! Communicate MPI rows to compute residual norm squared
    call xgBlock_colwiseNorm2(chebfi%AX%self,residu)
 
    ! Copy in Linalg representation (see chebfi_run, kept for reference)
    ! call xgBlock_copy(chebfi%X,X0)
    
    ! Transpose to colsrows state (X only)
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_COLSROWS)
        if (xmpi_comm_size(chebfi%spacecom) == 1) then 
            call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_transpose,2,tsec)
 
    ! Copy in ColsRows representation
    call xgBlock_copy(chebfi%xXColsRows, X0)

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
    
    ! Free temporary memory
    if (allocated(nrowsLinalg)) ABI_FREE(nrowsLinalg)

end subroutine polyfi_run
!!***

end module m_polyfi
!!***

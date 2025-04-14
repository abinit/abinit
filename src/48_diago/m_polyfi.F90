!!****f* ABINIT/m_polyfi
!! NAME
!! m_polyfi
!!
!! FUNCTION
!! This module contains the types and routines used to apply 
!! Polynomial Filtering method. It mainly defines 'polyfi' 
!! datatypes and associated methods. It generalizes m_chebfi2
!! for other filtering polynomials and arbitrary spectral bounds.
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

    !Several (private) parameters
    !-------------------------------------------------

    integer, parameter :: tim_getAX_BX     = 1754
    integer, parameter :: tim_transpose    = 1758

    ! Public 'polyfi' datatype
    !-------------------------------------------------
    type, public :: polyfi_t 

        type(chebfi_t) :: chebfi
        type(xg_t) :: DivResults
        real(dp) :: mineig_global
        real(dp) :: maxeig_global
        real(dp) :: lambda_minus
        real(dp) :: lambda_plus
        logical :: is_lowpass = .false.
        logical :: is_bandpass = .true.

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

 subroutine polyfi_init(polyfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp,&
         mineig_global,maxeig_global,lambda_minus,lambda_plus,ndeg_filter,space,&
         eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
         is_lowpass,is_bandpass,gpu_option,gpu_kokkos_nthrd,gpu_thread_limit)

     implicit none

    ! Arguments ------------------------------------
    integer       , intent(in   ) :: neigenpairs
    integer       , intent(in   ) :: bandpp
    integer       , intent(in   ) :: eigenProblem
    integer       , intent(in   ) :: me_g0
    integer       , intent(in   ) :: me_g0_fft
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
    logical       , intent(in   ) :: is_bandpass
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
    polyfi%is_bandpass = is_bandpass
    
    ! Oracle not used but passed to chebfi with deactivated values
    oracle = 0
    nbdbuf = 0
    oracle_factor = 1.d0
    oracle_min_occ = 0.d0

    ! Define chebfi object from Colsrows representation
    call chebfi_init(polyfi%chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp,&
        ndeg_filter,nbdbuf,space,1,comm,me_g0,me_g0_fft,paw,comm_rows,comm,&
        oracle,oracle_factor,oracle_min_occ,gpu_option,&
        gpu_kokkos_nthrd=gpu_kokkos_nthrd,gpu_thread_limit=gpu_thread_limit,&
        from_linalg=.false.)

    call polyfi_allocateAll(polyfi)

 end subroutine polyfi_init
 !!***
 
 !----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_allocateAll
!! NAME
!! polyfi_allocateAll

 subroutine polyfi_allocateAll(polyfi)

    implicit none

    type(polyfi_t), intent(inout) :: polyfi
    integer :: space_res
   
    if (polyfi%chebfi%space==SPACE_C) then
        space_res = SPACE_C
    else if (polyfi%chebfi%space==SPACE_CR) then
        space_res = SPACE_R
    else
        ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
    end if
    call xg_init(polyfi%DivResults, space_res, polyfi%chebfi%bandpp, 1, gpu_option=polyfi%chebfi%gpu_option)

 end subroutine
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
    call xg_free(polyfi%DivResults)

end subroutine polyfi_free
!!***

!----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_run
!! NAME
!! polyfi_run
!! 
!! FUNCTION
!! Diagonalize individual slice in parallel. Remember this is in distributed 
!! parallel region and there is no inter-slice communication.
!!
!! INPUT
!! X0     =eigenvector guess distributed for slice in colsrows representation
!! chebfi =data structure used to compute Chebyshev polynomial 
!! 
!! SOURCE

subroutine polyfi_run(polyfi,X0,getAX_BX,getBm1X,eigen,residu,nspinor)

    implicit none

    !Arguments ------------------------------------    
    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X0
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
    comm_size = xmpi_comm_size(spacecom)

    ! Arrays
    if(.not.allocated(nrowsLinalg)) ABI_MALLOC(nrowsLinalg,(comm_size))
    nrowsLinalg_ptr => nrowsLinalg
    nrowsLinalg(:) = polyfi%nrows_blockrows

    ! Note: we want to use the memory space of X0 but not the
    ! same pointer because it is common for all slices. For this
    ! reason we create a new xgBlock independent of X0 for slice.
    ! Do not do that: chebfi%xXColsRows = X0!!
    call xgBlock_setBlock(X0, chebfi%xXColsRows, total_spacedim, neigenpairs)

    if (polyfi%gpu_option==ABI_GPU_OFFLOAD) then
        ! Because X0 is on CPU but chebfi%xXColsRows on GPU
        ! FIXME either after or before setBlock
        ! call xgBlock_copy_from_gpu(chebfi%xXColsRows)
        ! call xgBlock_set_gpu_option(chebfi%xXColsRows,ABI_GPU_OFFLOAD)
    end if

    ! Restrict all communications to current subcommunicator
    call xgBlock_setComm(chebfi%xXColsRows,spacecom)
    call xgBlock_setComm(chebfi%X,spacecom)

    ! Create slice transposer using subcommunicator and allocate chebfi%X
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
    
    ! Apply polynomial filtering (requires colsrows state)
    if (polyfi%is_lowpass) then

        call polyfi_lowpass(polyfi,getAX_BX,getBm1X,nspinor)

    else if (polyfi%is_bandpass) then

        call polyfi_bandpass(polyfi,getAX_BX,getBm1X,nspinor)

    end if

    ! Wait filtering to finish before communication
    if (chebfi%paral_kgb==1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Transpose to linalg state
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then

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

    ! Wait until all MPI rows have computed their residual
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Communicate MPI rows to compute residual norm squared
    call xgBlock_colwiseNorm2(chebfi%AX%self,residu)

    
    ! Transpose to colsrows state (X only)
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_TRANSPOSE)
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

    if (polyfi%gpu_option==ABI_GPU_OFFLOAD) then
        ! TODO Deal with xgeigen and xgresidu
        !call xgBlock_reshape(chebfi%xgeigen, (/1,nband_slice/))
        !call xgBlock_reshape(chebfi%xgresidu, (/1,nband_slice/))
        !call xgBlock_copy_from_gpu(chebfi%xgeigen)
        !call xgBlock_copy_from_gpu(chebfi%xgresidu)
        !call slice_blockCopy(chebfi%xgeigen,sliceAll%xgeigen_ovlp,1,j1,nband_slice,j2)
        !call slice_blockCopy(chebfi%xgresidu,sliceAll%xgresidu_ovlp,1,j1,nband_slice,j2)
    end if

    ! Unitary test
    ABI_CHECK(cols(chebfi%xXColsRows)==ncols_slice,'wrong colsrows representation')
    write(*,'(a,i6,i6)') '# proc has # cols of slice X ', xmpi_comm_rank(spacecom), cols(chebfi%xXColsRows)

    ! if gpu: This is important! Because X0 is on CPU
    ! call xgBlock_copy_from_gpu(chebfi%xXColsRows)
    ! call xgBlock_set_gpu_option(chebfi%xXColsRows,ABI_GPU_DISABLED)

    ! Copy slice solution to the extended buffer (requires colsrows state)
    call xgBlock_copy(chebfi%xXColsRows,X0)
    ! FIXME same for eigen, residu?

    ! TODO we can also compute the merge on individual slices
    ! otherwise we do it outside?

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

end subroutine polyfi_run
!!***

!----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_lowpass
!! NAME
!! polyfi_lowpas
!!
!! FUNCTION
!! Apply Lowpass filter using Chebyshev polynomial on a set of vectors.
!!
!! INPUTS
!!  polyfi  = polynomial filtering datastructure
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  polyfi <type(polyfi_t)>=all data used to apply Polynomial Filtering algorithm 
!!  polyfi%chebfi%xXColsRows= Filtered vectors
!!
!! SOURCE

subroutine polyfi_lowpass(polyfi,getAX_BX,getBm1X,nspinors)

    implicit none

    ! Arguments ------------------------------------
    type(polyfi_t), intent(inout) :: polyfi
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
    integer :: ndeg_filter, ideg, ierr
    integer, allocatable :: ndeg_filter_bands(:)
    real(dp) :: center, radius, one_over_r, two_over_r

    ! *********************************************************************

    chebfi = polyfi%chebfi
    ndeg_filter = chebfi%ndeg_filter 

    if (chebfi%paral_kgb == 0) then
        ABI_MALLOC(ndeg_filter_bands,(chebfi%neigenpairs))
    else
        ABI_MALLOC(ndeg_filter_bands,(chebfi%bandpp))
    end if
    ndeg_filter_bands(:) = ndeg_filter
    
    ! Filter parameters
    center = (polyfi%lambda_plus + polyfi%lambda_minus)*0.5
    radius = (polyfi%lambda_plus - polyfi%lambda_minus)*0.5

    one_over_r = 1/radius
    two_over_r = 2/radius

    !A * Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ABI_NVTX_START_RANGE(NVTX_POLYFI_CORE)
    do ideg = 0, ndeg_filter - 1

        ABI_NVTX_START_RANGE(NVTX_POLYFI_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_POLYFI_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, chebfi%spacedim, chebfi%neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
        end if
        ABI_NVTX_END_RANGE()

        !A * Psi
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)

    end do ! ideg
    ABI_NVTX_END_RANGE()

    ! Scale X,AX,BX by amplification factor to reduce large values
    ! FIXME DivResults obscure
    ! il est juste utilisé pour multiplier le chebfi%xX,xAX,xBX par les valeurs de DivResults
    ! ce serait mieux de renommer DivResults à RRQ
    call chebfi_ampfactor(chebfi, polyfi%DivResults%self, polyfi%lambda_minus, &
        polyfi%lambda_plus, ndeg_filter_bands)

    if (allocated(ndeg_filter_bands)) then
        ABI_FREE(ndeg_filter_bands)
    end if

end subroutine polyfi_lowpass
!!***

!----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_bandpass
!! NAME
!! polyfi_bandpass
!!
!! FUNCTION
!! Apply Bandpass filter using Chebyshev-Jackson polynomial on a set of vectors.
!!
!! INPUTS
!!  polyfi  = polynomial filtering datastructure
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  polyfi <type(polyfi_t)>=all data used to apply Polynomial Filtering algorithm 
!!  polyfi%chebfi%xXColsRows= Filtered vectors
!!
!! SOURCE

subroutine polyfi_bandpass(polyfi,getAX_BX,getBm1X,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(polyfi_t), intent(inout) :: polyfi
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
    integer :: ndeg, n, ierr
    real(dp) :: center, radius, one_over_r, two_over_r
    real(dp) :: ls, us, cdeg, mu, damp
    type(xg_t) :: PolySum

    ! *********************************************************************

    chebfi = polyfi%chebfi
    ndeg = chebfi%ndeg_filter 
  
    ! Allocate memory for Chebyshev expansion of Heaviside step function
    call xg_init(PolySum, chebfi%space, chebfi%total_spacedim, chebfi%bandpp, chebfi%spacecom, &
        gpu_option=chebfi%gpu_option)
    call xgBlock_zero(PolySum%self)

    ! Filter parameters
    radius = (polyfi%maxeig_global - polyfi%mineig_global) / 2.d0
    center = (polyfi%maxeig_global + polyfi%mineig_global) / 2.d0

    one_over_r = 1/radius
    two_over_r = 2/radius

    ! Scaled slice bounds
    ls = (polyfi%lambda_minus - center) / radius
    us = (polyfi%lambda_plus - center) / radius
    
    ! A*Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ! TODO IL 10/3/2025 
    ! Deflate vectors ==========
    !call xg_Borthonormalize(chebfi%xXColsRows,chebfi%xBxColsRows,ierr,1,chebfi%gpu_option,AX=chebfi%xAXColsRows)

    ! Compute Chebyshev polynomial expansion on X iteratively on iline=0,nline
    ! X_next=X -> iline=0 Hamiltonian applications
    cdeg = Pi/(ndeg+2)
    mu = 1/Pi*(ACOS(ls)-ACOS(us))
    damp = 1.d0
    !Xsum = mu(0)*damp(0)*X_next + Xsum
    call xgBlock_saxpy(PolySum%self, mu*damp, chebfi%xXColsRows)

    ABI_NVTX_START_RANGE(NVTX_POLYFI_CORE)
    do n = 0, ndeg - 1  

        ! X_next=2/r*(AX_next-c*X_next)-X_prev
        ABI_NVTX_START_RANGE(NVTX_POLYFI_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, n, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ! xXColsRows=X_next
        call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)

        !Xsum = damp(i+1)*mu(i+1)*X_next + Xsum
        mu = 2/Pi * (SIN(*ACOS(ls)) - SIN((n+1)*ACOS(us)))/(n+1)
        damp = ((1 - (n+1)/(ndeg+2))*SIN(cdeg)*COS((n+1)*cdeg) + 1/(ndeg+2)*COS(cdeg)*SIN((n+1)*cdeg))/SIN(cdeg)
        call xgBlock_saxpy(PolySum%self, mu*damp, chebfi%xXColsRows)
   
        ! Store term before exit
        ! AX_next=A*X_next -> iline+2 Hamiltonian applications
        if (n==ndeg-1) then
            ! X_next=Xsum (copy Xsum to X_next)
            call xgBlock_copy(PolySum%self, chebfi%xXColsRows)
        end if

        ! Apply A and B to X
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)
    
    end do ! end iline
    ABI_NVTX_END_RANGE()

    ! Free Chebyshev expansion workspace
    call xg_free(PolySum)

end subroutine polyfi_bandpass
!!***

end module m_polyfi
!!***

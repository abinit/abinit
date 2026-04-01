!!****f* ABINIT/m_trace_estimation
!! NAME
!! m_trace_estimation
!!
!! FUNCTION
!! This module contains routines used to compute the stochastic Lanczos trace estimation.
!! It also computes cummulative eigenvalue counts using differences of trace estimation
!! on consecutive intervals. 
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

module m_trace_estimation

    use defs_basis
    use defs_abitypes
    use m_abicore
    use m_errors
    use m_time, only : timab
    use m_sort, only: sort_dp

    use m_cgtools
    use m_xg
    use m_xgTransposer
    
    use m_chebfi2
    use m_polynomial_filter
    use m_slice_task, only: matrixInfo_t

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
    integer, parameter :: tim_getAX_BX    = 1754
    integer, parameter :: tim_invovl      = 1755
    integer, parameter :: tim_lanczos     = 2167
    integer, parameter :: tim_trace       = 2168

    ! Public methods
    !-------------------------------------------------
    public :: computeBLanczos           ! for lower bound estimation
    public :: get_eigenvalue_count      ! count eigenvalues for lowpass intervals 
    public :: computeTraceEstimation    ! compute moments etc 
    public :: smallestTridiagEigenpair  ! used in slice_cprj (experimental)

    CONTAINS  
!=====================================================================
!!***

!!****f* m_trace_estimation/computeBLanczos
!! NAME
!! computeBLanczos
!! 
!! FUNCTION
!! B-Lanczos three-term recurrence (using B-inner product).
!! Performs k Lanczos iterations on a column vector.
!!
!! SOURCE
  
  subroutine computeBLanczos(minfo, paw, getAX_BX, getBm1X, k, lambda_min, res_norm)

    implicit none

    type(matrixInfo_t), intent(in) :: minfo
    logical, intent(in) :: paw
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

    integer :: j
    integer :: rank
    integer :: me_g0, space, tot_spacedim, gpu_option
    real(dp) :: tsec(2)

  ! *********************************************************************

    call timab(tim_lanczos,1,tsec)
    
    tot_spacedim = minfo%total_spacedim
    gpu_option = minfo%gpu_option
    space = minfo%space
    me_g0 = minfo%me_g0
    if (minfo%paral_kgb==1) then
        me_g0 = minfo%me_g0_fft
    end if
    rank = xmpi_comm_rank(minfo%spacecom)
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
        if (paw) then
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

!!****f* m_trace_estimation/computeChebyshevMoments
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

subroutine computeChebyshevMoments(minfo, tolerance, ecut, paw, X0, getAX_BX, getBm1X, &
        lambda_minus, lambda_plus, ndeg_filter, cheby_moments)

    implicit none

    type(matrixInfo_t), intent(in) :: minfo
    real(dp), intent(in) :: ecut, tolerance
    logical, intent(in) :: paw
    type(xgBlock_t), intent(inout) :: X0
    integer, intent(in) :: ndeg_filter
    real(dp), intent(in) :: lambda_minus, lambda_plus
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

    integer :: ideg
    integer :: nband, tot_spacedim, space, me_g0, gpu_option
    real(dp) :: center, radius
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: tsec(2)
    type(xg_t) :: Moments
    type(xg_t) :: X0_backup
    type(xgBlock_t) :: Moment_ideg
    type(chebfi_t) :: chebfi
    complex(dp), pointer :: momvals(:,:) => null()

    ! *********************************************************************

    ! todo add timer
    ! tim_cheby_moments

    gpu_option = minfo%gpu_option
    tot_spacedim = minfo%total_spacedim
    space = minfo%space
    me_g0 = minfo%me_g0
    if (minfo%paral_kgb==1) then
        me_g0 = minfo%me_g0_fft
    end if
    nband = cols(X0)
 
    ! Moment workspace size (1, ndeg+1)
    call xg_init(Moments, space, 1, ndeg_filter+1, gpu_option=gpu_option) ! M_n=<X0,f_n(A)X0>
    call xg_init(X0_backup, space, tot_spacedim, nband, minfo%spacecom, me_g0=me_g0, gpu_option=gpu_option) ! X0

    ! Initialize chebfi object in MPI Colsrows distribution
    call chebfi_init(chebfi,nband,tot_spacedim,tolerance,ecut,minfo%paral_kgb,&
        nband,ndeg_filter,0,minfo%space,1,xmpi_comm_null,minfo%me_g0,minfo%me_g0_fft,&
        paw,minfo%comm_rows,minfo%comm_cols,0,1.d0,0.d0,minfo%gpu_option,&
        gpu_kokkos_nthrd=minfo%gpu_kokkos_nthrd,gpu_thread_limit=minfo%gpu_thread_limit,&
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
    call chebfi_free(chebfi)
    call xg_free(X0_backup)
    
end subroutine computeChebyshevMoments
!!***

!----------------------------------------------------------------------

!!****f* m_trace_estimation/computeTraceEstimation
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

subroutine computeTraceEstimation(minfo, ecut, paw, tolerance, getAX_BX, getBm1X, ndeg_filter, &
        m_probe, min_low_bound, moments)

    implicit none

    type(matrixInfo_t), intent(in) :: minfo
    integer, intent(in) :: ndeg_filter, m_probe
    logical, intent(in) :: paw
    real(dp), intent(in) :: ecut, min_low_bound, tolerance
    real(dp), intent(inout) :: moments(:)
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

    !integer :: nband
    integer :: tot_spacedim, gpu_option
    integer :: ierr
    integer :: me_g0
    integer :: jcol
    integer :: my_rank
    integer :: seed
    integer :: m_probe_tot
    integer :: num_moments
    type(xg_t) :: X_probe
    complex(dp), allocatable :: cheby_moments(:,:)
    real(dp) :: tsec(2)
    
    ! *********************************************************************
    
    call timab(tim_trace,1,tsec)

    tot_spacedim = minfo%total_spacedim
    gpu_option = minfo%gpu_option
    my_rank = xmpi_comm_rank(minfo%spacecom)
    me_g0 = minfo%me_g0
    if (minfo%paral_kgb==1) then
        me_g0 = minfo%me_g0_fft
    end if

    num_moments = ndeg_filter + 1
    ABI_MALLOC(cheby_moments, (1, num_moments) ) 

    ! total number of probes is m_probes * number of MPI processes
    call xg_init(X_probe, minfo%space, tot_spacedim, m_probe, minfo%spacecom, &
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
    call computeChebyshevMoments(minfo, tolerance, ecut, paw, X_probe%self, &
        getAX_BX, getBm1X, min_low_bound, ecut, ndeg_filter, cheby_moments)

    ! Sum real part of moments and divide by number of probes
    ! moments(k) = 1/Nv * Sum_{i=1}^Nv v_i^T T_k(A)v_i
    m_probe_tot = m_probe
    call xmpi_sum(m_probe_tot, minfo%spacecom, ierr)
    !moments(1:num_moments) = (/ (real(cheby_moments(1,k)), k=1,num_moments) /)
    moments = real(cheby_moments(1,:))
    call xmpi_sum(moments, minfo%spacecom, ierr)
    moments(1:num_moments) = moments(1:num_moments)/m_probe_tot
    
    !write(std_out,*) 'moments k=0=', real(cheby_moments(1,1))
    !write(std_out,*) 'moments k=1=', real(cheby_moments(1,2))
    !write(std_out,*) 'moments k=3=', real(cheby_moments(1,3))
    !flush(std_out)

    ! Free memory
    ABI_FREE(cheby_moments)
    call xg_free(X_probe)
    
    call timab(tim_trace,2,tsec)
    
end subroutine computeTraceEstimation
!!***

!----------------------------------------------------------------------

!!****f* m_trace_estimation/get_eigenvalue_count
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
    !real(dp) :: sigma, alpha
    !integer :: Ngrid

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

!----------------------------------------------------------------------

!!****f* m_trace_estimation/smallestTridiagEigenpair
!! NAME
!! smallestTridiagEigenpair
!! 
!! FUNCTION
!! Smallest eigenvalue of symmetric tridiagonal
!!
!! SOURCE

  subroutine smallestTridiagEigenpair(n, d, e, lambda, v)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in)  :: d(n), e(n-1)
    real(dp), intent(out) :: lambda
    real(dp), intent(out) :: v(n)

    ! Local copies (DSTEVX overwrites input)
    real(dp) :: dloc(n), eloc(n-1)
    real(dp), allocatable :: z(:,:), work(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info, m
    
    ! *********************************************************************

    dloc = d
    eloc = e

    ABI_MALLOC(z, (n,1))
    ABI_MALLOC(work, (5*n))
    ABI_MALLOC(iwork, (5*n))
    ABI_MALLOC(ifail, (n))

    ! DSTEVX computes selected eigenpairs (here smallest: index 1)
    call dstevx('V', 'I', n, dloc, eloc, 0.0d0, 0.0d0, 1, 1, 1.0d-12, m, dloc, z, &
        n, work, iwork, ifail, info)

    if (info /= 0) then
       ABI_ERROR('DSTEVX failed')
    end if

    lambda = dloc(1)
    v      = z(:,1)

    ABI_FREE(z)
    ABI_FREE(work)
    ABI_FREE(iwork)
    ABI_FREE(ifail)

  end subroutine smallestTridiagEigenpair
!!***

end module m_trace_estimation
!!***

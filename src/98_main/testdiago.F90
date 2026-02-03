!!****p* ABINIT/testdiago
!! NAME
!! testdiago
!!
!! FUNCTION
!! Testing algorithms related to diagonalisation of Hermitian matrices.
!!
!! COPYRIGHT
!! Copyright (C) 2017-2026 ABINIT group (IML)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

program testdiago

    use m_xg
    use m_xgTransposer
    use m_xmpi
    use m_time
    use defs_basis
    use m_profiling_abi
    use m_slice_cprj
    use m_errors

    implicit none

  integer, parameter :: n = 500
  integer, parameter :: kmax = 20

  real(dp) :: A(n,n), B(n,n)
  real(dp) :: evals(n)
  real(dp) :: alpha(kmax), beta(kmax+1)
  real(dp) :: v0(n)
  real(dp) :: wB(n)
  real(dp) :: evalsB(n)
  real(dp) :: Bwork(n,n)
  real(dp) :: v_lanczos(n)
  real(dp) :: eig_lanczos
  real(dp) :: shift
  real(dp) :: res_norm
  real(dp) :: low_bound
  real(dp), allocatable :: workB(:)
  integer, allocatable :: seed(:)

  integer :: i, info, lwork
  integer :: nseed
  logical :: success

  call random_seed(size=nseed)
  allocate(seed(nseed))
  ! Case n=500 with kmax=10 F and kmax=20 T
  !seed = (/ 1239937284, -897626825, 1137443096, -1034675534, &
  !          124992374, 488443746, -644767705, -1072627500 /)
  ! Case n=500 with kmax=10 F and kmax=20 T
  !seed = (/ 1555480729, 1032432162, 113869316, 1593505229, &
  !          461320453, 1082979611, -238973894, -1158470550 /)
  !call random_seed(put=seed)
  call random_seed(get=seed)
  print *, "Seed=", seed
  deallocate(seed)

  ! -------------------------------------------------
  ! 1. Random symmetric matrix A
  ! -------------------------------------------------
  call random_number(A)
  A = 0.5_dp * (A + transpose(A))

  ! -------------------------------------------------
  ! 2. Random symmetric positive definite matrix B
  !    B = C^T C + n I
  ! -------------------------------------------------
  call random_number(B)
  B = matmul(transpose(B), B)
  Bwork = B
  lwork = max(1, 4*n)
  allocate(workB(lwork))
  call dgeev('N','N',n,Bwork,n,evalsB,wB,0,1,0,1,workB,lwork,info)
  if (info /= 0) stop "DSYEV failed"
  shift = minval(evalsB) - 0.1d0
  do i = 1, n
    B(i,i) = B(i,i) - shift
  end do
  deallocate(workB)

  ! -------------------------------------------------
  ! 3. Exact generalized eigenvalue solve
  !    A x = lambda B x
  ! -------------------------------------------------
  call exact_lowest_gep(n, A, B, evals)

  print *, "Exact lowest generalized eigenvalue:"
  print *, evals(1)

  ! -------------------------------------------------
  ! 4. B-Lanczos approximation
  ! -------------------------------------------------
  call random_number(v0)
  call b_lanczos(A, B, n, kmax, eig_lanczos, res_norm)

  print *, "B-Lanczos lowest eigenvalue approximation (upper bound):"
  print *, eig_lanczos

  print *, "residual error norm"
  print *, res_norm

  low_bound = eig_lanczos - res_norm
  print *, "Final estimation (guaranteed lower bound)"
  print *, low_bound

  success = evals(1) > low_bound
  print *, "Status"
  print *, success

  if (.not. success) then
      print *, "must increase kmax"
  end if

  ! --- Call block B-Lanczos ---
  call block_b_lanczos(A, B, n, kmax, eig_lanczos, res_norm)

  print *, "Lowest Ritz value from block B-Lanczos:", eig_lanczos
  print *, "Residual norm:", res_norm

contains

! =====================================================
! B-inner product
! =====================================================
real(dp) function bdot(n, x, y, B)
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n), y(n), B(n,n)
  real(dp) :: tmp(n)

  tmp = matmul(B, x)
  bdot = dot_product(y, tmp)
end function bdot

  !-----------------------------
  ! Factor B (Cholesky)
  !-----------------------------
  subroutine factor_B(B, n, info)
    integer, intent(in) :: n
    integer, intent(out) :: info
    real(dp), intent(inout) :: B(n,n)

    call dpotrf('U', n, B, n, info)
    if (info /= 0) then
       print *, "B is not SPD!"
       stop
    end if
  end subroutine factor_B

  !-----------------------------
  ! Binv_mul: v = B^{-1} u
  !-----------------------------
  subroutine Binv_mul(B, n, u, v)
    integer, intent(in) :: n
    real(dp), intent(in) :: B(n,n), u(n)
    real(dp), intent(out) :: v(n)
    integer :: nrhs, info

    v = u
    nrhs = 1
    call dpotrs('U', n, nrhs, B, n, v, n, info)
    if (info /= 0) then
       print *, "DSYTRS failed!"
       stop
    end if
  end subroutine Binv_mul

! =====================================================
! B-Lanczos three-term recurrence
! =====================================================
subroutine b_lanczos(A, B, n, k, lambda_min, res_norm)
    integer, intent(in) :: n, k
    real(dp), intent(in) :: A(n,n), B(n,n)
    real(dp), intent(out) :: lambda_min, res_norm

    ! Arrays
    real(dp), allocatable :: q(:), q_prev(:), v(:), Bv(:)
    real(dp), allocatable :: alpha(:), beta(:)
    real(dp) :: beta_prev, normB, theta
    integer :: j, info
    real(dp), allocatable :: v_min(:)
    integer :: i

    ! Copy B to factor
    real(dp), allocatable :: Bfact(:,:)
    allocate(Bfact(n,n))
    Bfact = B
    call factor_B(Bfact, n, info)

    ! Allocate
    allocate(q(n), q_prev(n), v(n), Bv(n))
    allocate(v_min(k))
    allocate(alpha(k), beta(k-1))

    ! Initialize q randomly and B-normalize
    call random_number(q)
    Bv = matmul(B, q)
    normB = sqrt(dot_product(q, Bv))
    q = q / normB

    q_prev = 0.0_dp
    beta_prev = 0.0_dp

    do j = 1, k
        ! Aq = A * q
        v = matmul(A, q)

        ! α_j = q^T * Aq
        alpha(j) = dot_product(q, v)

        ! v = B^{-1} A q - α q - β_prev q_prev
        call Binv_mul(Bfact, n, v, v)
        v = v - alpha(j)*q
        if (j > 1) v = v - beta_prev*q_prev

        ! Compute β_j if j<k
        if (j < k) then
            Bv = matmul(B, v)
            beta(j) = sqrt(dot_product(v, Bv))

            ! Update q_prev, q, beta_prev
            q_prev = q
            q = v / beta(j)
            beta_prev = beta(j)
        end if
    end do

    ! Diagonalize T
    call smallestTridiagEigenpair(k, alpha, beta, lambda_min, v_min)

    ! residual norm using Lanczos shortcut
    res_norm = abs(beta(k-1)*v_min(k))

    ! Clean up
    deallocate(q, q_prev, v, Bv, v_min, alpha, beta, Bfact)

end subroutine b_lanczos

! =====================================================
! Exact lowest generalized eigenvalue via LAPACK
! =====================================================
subroutine exact_lowest_gep(n, A, B, evals)
  integer, intent(in) :: n
  real(dp), intent(in) :: A(:,:), B(:,:)
  real(dp), intent(out) :: evals(n)
  
  real(dp), allocatable :: Awork(:,:), Bwork(:,:)
  real(dp), allocatable :: rwork(:)
  integer :: info
  integer :: lwork
  integer :: lrwork

  allocate(Awork(n,n), Bwork(n,n)) ! force contiguous memory
  Awork = A
  Bwork = B

  lwork  = 3*n
  lrwork = 2*n*n + 6*n + 1
  allocate(rwork(lrwork))

  call dsygv(1,'N','U', n, Awork, n, Bwork, n, evals, rwork, lrwork, info)
  if (info /= 0) stop "DSYGV failed"

  deallocate(Awork, Bwork, rwork)
end subroutine exact_lowest_gep


  !-----------------------------
  ! Binv multiplication for each column of block
  !-----------------------------
  subroutine Binv_mul_block(Bfact, n, U, V)
    integer, intent(in) :: n
    real(dp), intent(in) :: Bfact(n,n), U(n,n)
    real(dp), intent(out) :: V(n, size(U,2))
    integer :: i, nrhs, info

    do i = 1, size(U,2)
       V(:,i) = U(:,i)
       nrhs = 1
       call dpotrs('U', n, nrhs, Bfact, n, V(:,i), n, info)
       if (info /= 0) stop "DSYTRS failed in Binv_mul_block"
    end do
  end subroutine Binv_mul_block

  !-----------------------------
  ! B-orthonormalize a block Q (n x k)
  ! Q^T B Q = I
  !-----------------------------
  subroutine B_orthonormalize(Q, B, n, k)
    integer, intent(in) :: n, k
    real(dp), intent(in) :: B(n,n)
    real(dp), intent(inout) :: Q(n,k)
    real(dp), allocatable :: M(:,:), L(:,:)
    integer :: info

    allocate(M(k,k), L(k,k))

    ! Compute M = Q^T B Q
    M = matmul(transpose(Q), matmul(B, Q))

    ! Cholesky: M = L^T L
    L = M
    call dpotrf('U', k, L, k, info)
    if (info /= 0) stop "Cholesky failed in B_orthonormalize"

    ! Q <- Q * inv(L^T)
    call dtrsm('R', 'U', 'T', 'N', n, k, 1.0_dp, L, k, Q, n)

    deallocate(M,L)
  end subroutine B_orthonormalize

  ! Simple Euclidean norm
   real(dp) function norm2(x)
      real(dp), intent(in) :: x(:)
      norm2 = sqrt(dot_product(x,x))
    end function norm2

  !-----------------------------
  ! Generate a B-orthonormal independent random block
  ! Q(:,1:k) will satisfy Q^T B Q = I
  !-----------------------------

subroutine generate_independent_block(Q, B, n, k)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: n, k
  real(dp), intent(in) :: B(n,n)
  real(dp), intent(out) :: Q(n,k)

  integer :: i, j
  real(dp) :: normB
  real(dp), allocatable :: v(:), Bv(:)

  allocate(v(n), Bv(n))

  do j = 1, k
     ! Generate a random vector
     call random_number(v)

     ! Make it B-orthogonal to previous columns
     do i = 1, j-1
        Bv = matmul(B, Q(:,i))
        v = v - (dot_product(v, Bv)) * Q(:,i)
     end do

     ! Normalize w.r.t. B-inner product
     Bv = matmul(B, v)
     normB = sqrt(dot_product(v, Bv))
     Q(:,j) = v / normB
  end do

  deallocate(v, Bv)
end subroutine generate_independent_block

  !-----------------------------
  ! Single-block B-Lanczos
  !-----------------------------
  subroutine block_b_lanczos(A, B, n, k, lambda_min, res_norm)
    integer, intent(in) :: n, k
    real(dp), intent(in) :: A(n,n), B(n,n)
    real(dp), intent(out) :: lambda_min, res_norm

    real(dp), allocatable :: Q(:,:), V(:,:), T(:,:), evals(:), eigvecs(:,:)
    real(dp), allocatable :: Bfact(:,:)
    real(dp), allocatable :: work(:)
    integer :: info, lwork

    ! Allocate
    allocate(Q(n,k), V(n,k), T(k,k), evals(k), eigvecs(k,k), Bfact(n,n))

    ! Factor B once
    Bfact = B
    call dpotrf('U', n, Bfact, n, info)
    if (info /= 0) stop "B not SPD"

    ! Initialize random block Q
    !call random_number(Q)
    !call B_orthonormalize(Q, B, n, k)

    ! --- Generate independent random block Q ---
    call generate_independent_block(Q, B, n, k)

    ! Apply A once to the block
    V = matmul(A, Q)

    ! Projected T = Q^T A Q
    T = matmul(transpose(Q), V)

    ! Diagonalize T
    lwork = 2*k*k + 6*k + 1
    allocate(work(lwork))
    call dsyev('V','U', k, T, k, evals, work, lwork, info)
    if (info /= 0) stop "DSYEV failed"

    lambda_min = evals(1)

    ! Compute residual norm for the lowest Ritz value
    ! r = A*q_min - lambda_min * B*q_min
    ! Use first column of eigenvectors of T
    eigvecs = T
    call dsyev('V','U', k, eigvecs, k, evals, work, lwork, info)
    res_norm = norm2(matmul(V, eigvecs(:,1)) - lambda_min*matmul(B, Q(:,1)))

    ! Deallocate
    deallocate(Q,V,T,evals,eigvecs,Bfact,work)
  end subroutine block_b_lanczos

!**********************************************************************

 end program testdiago
!!***


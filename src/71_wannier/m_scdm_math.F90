!!****m*ABINIT/m_scdm_math
!! NAME
!!  m_scdm_math
!!
!! FUNCTION
!!  Math functions used by SCDM
!! (select columns of density matrix method) for getting wannier functions.
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2022 ABINIT group (hexu)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif



#include "abi_common.h"

module m_scdm_math
  use defs_basis
  use m_abicore
  use m_errors
  implicit none
  real(dp), parameter, public:: tpi = 2*PI
  complex(dp), parameter, public:: tpi_im = cmplx(0.0_dp, tpi, kind = dp)

  type, public:: eigensolver
     integer:: ndim = -1, lwork = 0
     complex(dp), allocatable  :: work(:)
     real(dp), allocatable:: rwork(:)
   contains
     procedure:: run => eigensolver_run
     procedure:: finalize => eigensolver_finalize
  end type eigensolver


  public:: complex_QRCP_piv_only
  public:: real_svd
  public:: complex_svd
  public:: gaussian
  public:: fermi
  public:: insertion_sort_double
  public:: build_Rgrid
  private

contains

  subroutine complex_QRCP_Piv_only(A, Piv)
    complex(dp), intent(in):: A(:, :)
    integer, intent(inout):: Piv(:)
    real(dp):: rwork(size(A, 2)*2)
    complex(dp):: tau(min(size(A, 1), size(A, 2)))
    integer:: m, n
    complex(dp), allocatable:: work(:)
    complex(dp):: tmp_work(2)
    integer:: lwork
    integer:: info
    EXTERNAL ZGEQP3

    m = size(A, 1)
    n = size(A, 2)
    rwork(:)=0.0_dp
    call ZGEQP3(m, n, A, m, piv, tau, tmp_work, -1, rwork, info)
    if(info /= 0) then
       ABI_ERROR("Error in doing QRCP")
    endif
    Piv(:) = 0
    rwork(:) = 0.0_DP
    lwork = INT(AINT(REAL(tmp_work(1))))

    ABI_MALLOC(work, (lwork))
    work(:) = (0.0_DP, 0.0_DP)

    CALL ZGEQP3(m, n, A, m, piv, tau, work, lwork, rwork, info)
    if(info /= 0) then
       ABI_ERROR("Error in doing QRCP")
    endif
    ABI_SFREE(work)
  end subroutine complex_QRCP_Piv_only



  subroutine real_svd(A,  U, S, VT)
    real(dp), intent(in):: A(:, :)
    real(dp), intent(inout):: U(:, :), S(:), VT(:,:)
    integer:: LWMAX
    real(dp):: tmp(5)
    real(dp), allocatable::  WORK(:)
    integer  :: M, N
    integer::         LDA, LDU, LDVT
    integer::         INFO, LWORK
    EXTERNAL         DGESVD
    M = size(A, 1)
    N = size(A, 2)
    LDA = M
    LDU = M
    LDVT = N
    LWORK = -1
    LWMAX = max(size(A, 1), size(A, 2))*10
    CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, &
         & tmp, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( tmp( 1 ) ) )

    ABI_MALLOC(work, (lwork))
    CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, &
         WORK, LWORK, INFO )
    IF( INFO .GT. 0 ) THEN
       ABI_ERROR('The algorithm computing SVD failed to converge.')
       STOP
    END IF
    ABI_SFREE(work)
  end subroutine real_svd

  subroutine complex_svd(A,  U, S, VT, mode)
    complex(dp), intent(in):: A(:, :)
    complex(dp), intent(inout):: U(:, :),  VT(:,:)
    real(dp), intent(inout):: S(:)
    character, intent(in):: mode  ! A, or S, or N/O
    integer:: LWMAX
    complex(dp):: tmp(2)
    complex(dp), allocatable::  WORK(:)
    real(dp), allocatable:: rwork(:)
    integer  :: M, N
    integer::         LDA, LDU, LDVT
    integer::         INFO, LWORK
    EXTERNAL         ZGESVD

    M = size(A, 1)
    N = size(A, 2)
    LWMAX = max(size(A, 1), size(A, 2))*10
    LDA = M
    LDU = M
    LDVT = N
    LWORK = -1
    ABI_MALLOC(rwork, (min(M, N)*6))
    CALL ZGESVD( mode, mode, M, N, A, LDA, S, U, LDU, VT, LDVT, &
         & tmp, LWORK, rwork, INFO )
    LWORK = MIN( LWMAX, INT( tmp( 1 ) ) )
    ABI_MALLOC(work, (lwork))
    CALL ZGESVD( mode, mode, M, N, A, LDA, S, U, LDU, VT, LDVT, &
         WORK, LWORK, rwork, INFO )
    IF( INFO .GT. 0 ) THEN
       ABI_ERROR('The algorithm computing SVD failed to converge.')
       STOP
    END IF
    ABI_SFREE(work)
    ABI_SFREE(rwork)
  end subroutine complex_svd


  function gaussian(x, mu, sigma) result(y)
    real(dp), intent(in):: x, mu, sigma
    real(dp):: y
    y = exp(-1.0 * (x-mu)**2/sigma**2)
  end function gaussian

  !===============================================================
  ! Complementary error function.
  !===============================================================
  pure function erfcd(x) result(y)
    real(dp), intent(in):: x
    real(dp):: y
    real(dp):: t, z
    z = abs(x)
    t = 1.0 / ( 1.0+0.5*z )

    y = t*exp( -z*z - 1.26551223+t *   &
         ( 1.00002368+t * ( 0.37409196+t * &
         ( 0.09678418+t * (-0.18628806+t * &
         ( 0.27886807+t * (-1.13520398+t * &
         ( 1.48851587+t * (-0.82215223+t * 0.17087277 )))))))))
    if ( x .lt. 0.0 ) y = 2.0-y

  end function erfcd

  function fermi(x, mu, sigma) result(y)
    real(dp), intent(in):: x, mu, sigma
    real(dp):: y
    y = 0.5*erfc((x-mu) / sigma)
  end function fermi

  subroutine eigensolver_run(self, evals, evecs)
    class(eigensolver), intent(inout):: self
    real(dp), intent(inout):: evals(:)
    complex(dp), intent(inout):: evecs(:,:)
    integer ::  info
    external ZHEEV
    if (self%ndim == -1) then
       self%ndim = size(evecs, 1)
       self%lwork = -1
       ABI_MALLOC(self%work, (1))
       ABI_MALLOC(self%rwork, (3*size(evecs, 1)-2))
       call ZHEEV('V', 'U', self%ndim, evecs, self%ndim, evals, self%work, self%lwork, self%rwork, info)
       self%lwork = INT(self%work(1))
       ABI_SFREE(self%work)
       ABI_MALLOC(self%work, (self%lwork))
    else if (self%ndim /= size(evecs, 1)) then
       ABI_ERROR("Eigensovler: The size of the evecs is not the same as previous one.")
    end if
    call ZHEEV('V', 'U', self%ndim, evecs, self%ndim, evals, self%work, self%lwork, self%rwork, info)
  end subroutine eigensolver_run

  subroutine eigensolver_finalize(self)
    class(eigensolver), intent(inout):: self
    self%ndim = -1
    self%lwork = -1
    ABI_SFREE(self%work)
    ABI_SFREE(self%rwork)
  end subroutine eigensolver_finalize

  !----------------------------------------------------------------------
  !> @brief insertion_sort_int: sort a array using insertion sort algorithm
  !>  it is a memory safe method but is generally slow.
  !> @param[inout]  a: the array to be sorted. and will output inplace
  !> @param[inout] order (optional) the sorted index, it can be used to sort
  !>  other arrays so that the order in consistent.
  !----------------------------------------------------------------------
  subroutine insertion_sort_double(a, order)
    real(dp), intent(inout):: a(:)
    integer, optional, intent(inout):: order(size(a))
    integer:: n, i, j
    real(dp):: v
    n = size(a)
    if (present(order)) then
       do i = 1, n
          order(i)=i
       end do
    end if
    do i = 2, n
       v = a(i)
       j = i-1
       do while(j >= 1 )
          if (a(j)<=v) exit
          a(j+1)=a(j)
          if(present(order)) order(j+1)=order(j)
          j = j-1
       end do
       a(j+1)=v
       if(present(order)) order(j+1)=i
    end do

  end subroutine insertion_sort_double


  ! -4/2-> -1, 4/2->2. Then (-1, 0, 1, 2)
  pure function div2(x) result(y)
    integer, intent(in):: x
    integer:: y
    if(mod(x, 2)==0 .and. x < 0) then 
       y = x/2+1
    else
       y = x/2
    endif
  end function div2

  subroutine build_Rgrid(kmesh, Rlist)
    integer, intent(in):: kmesh(3)
    integer, allocatable, intent(inout):: Rlist(:,:)
    integer:: n, i1, i2, i3, i

    n = kmesh(1) * kmesh(2) *kmesh(3)
    ABI_MALLOC(Rlist, (3, n))
    i = 0
    ! Note that C/Fortran integer division is "truncate towards 0" division, 
    ! whereas Python one is "floor" division.
    ! For C/Fortran, the behavior for even and odd numbers is 
    !  not consistent and need special treatment in div2.
    do i3 = div2(-kmesh(3)), div2(kmesh(3))
       do i2 = div2(-kmesh(2)), div2(kmesh(2))
          do i1 = div2(-kmesh(1)), div2(kmesh(1))
             i = i+1
             Rlist(:, i) = [i1, i2, i3]
          end do
       end do
    end do
  end subroutine build_Rgrid



end module m_scdm_math

!!***

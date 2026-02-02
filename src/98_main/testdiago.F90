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

program testdiago

    use iso_fortran_env, only: output_unit
    implicit none

    !Local variables-------------------------------
    !scalars
    real :: xTx
    real(8) :: lambda
    integer :: n
    !arrays
    real :: x(10)
    real(8), allocatable :: d(:), e(:), v(:)

!**********************************************************************

    call random_number(x) 
    xTx = dot_product(x, x)
    write(output_unit,*) 'testdiago OK', xTx

      ! Matrix size
  n = 10
  allocate(d(n), e(n-1), v(n))

  ! Generate random tridiagonal matrix
  call generate_random_tridiagonal(n, d, e)

  ! Compute smallest eigenpair
  call smallest_tridiag_eigenpair(n, d, e, lambda, v)

  ! Print results
  write(output_unit,*) 'Smallest eigenvalue =', lambda
  write(output_unit,*) 'Corresponding eigenvector:'
  write(output_unit,*) v

contains

  !------------------------------------------
  subroutine generate_random_tridiagonal(n, d, e)
    implicit none
    integer, intent(in) :: n
    real(8), intent(out) :: d(n), e(n-1)
    integer :: i

    call random_seed()

    do i = 1, n
       call random_number(d(i))
    end do

    do i = 1, n-1
       call random_number(e(i))
    end do
  end subroutine generate_random_tridiagonal

  !------------------------------------------
  subroutine smallest_tridiag_eigenpair(n, d, e, lambda, v)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in)  :: d(n), e(n-1)
    real(8), intent(out) :: lambda
    real(8), intent(out) :: v(n)

    ! Local copies (DSTEVX overwrites input)
    real(8) :: dloc(n), eloc(n-1)
    real(8), allocatable :: z(:,:), work(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info, m

    dloc = d
    eloc = e

    allocate(z(n,1))
    allocate(work(5*n))
    allocate(iwork(5*n), ifail(n))

    ! DSTEVX computes selected eigenpairs (here smallest: index 1)
    call dstevx('V', 'I', n, dloc, eloc, 0.0d0, 0.0d0, 1, 1, 1.0d-12, m, dloc, z, n, work, iwork, ifail, info)

    if (info /= 0) then
       write(output_unit,*) 'DSTEVX failed, INFO =', info
       stop
    end if

    lambda = dloc(1)
    v      = z(:,1)

    deallocate(z, work, iwork, ifail)
  end subroutine smallest_tridiag_eigenpair

!**********************************************************************

 end program testdiago
!!***


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

    !Local variables-------------------------------
    !scalars
    real :: xTx
    real(dp) :: lambda
    integer :: n
    !arrays
    real :: x(10)
    real(dp), allocatable :: d(:), e(:), v(:)

!**********************************************************************

    call random_number(x) 
    xTx = dot_product(x, x)
    write(std_out,*) 'testdiago OK', xTx

    ! Matrix size
    n = 10
    allocate(d(n), e(n-1), v(n))

    ! Generate random tridiagonal matrix
    call generate_random_tridiagonal(n, d, e)

    ! Compute smallest eigenpair
    call smallestTridiagEigenpair(n, d, e, lambda, v)

    ! Print results
    write(std_out,*) 'Smallest eigenvalue =', lambda
    write(std_out,*) 'Corresponding eigenvector Yeyyy:'
    write(std_out,*) v

contains

  !------------------------------------------
  subroutine generate_random_tridiagonal(n, d, e)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: d(n), e(n-1)
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
  ! new routine ... 

!**********************************************************************

 end program testdiago
!!***


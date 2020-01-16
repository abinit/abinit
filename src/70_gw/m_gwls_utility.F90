!!****m* ABINIT/m_gwls_utility
!! NAME
!! m_gwls_utility
!!
!! FUNCTION
!!  .
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (JLJ, BR, MC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


!---------------------------------------------------------------------
!  Utility modules, which are not directly related to the
!  problem to be solved.
!---------------------------------------------------------------------

module m_gwls_utility
!----------------------------------------------------------------------------------------------------
! This module contains useful functions and subroutines, which are not necessarily numerical
! in nature, for timing, opening files, etc...
!----------------------------------------------------------------------------------------------------

! abinit modules
use defs_basis
use m_abicore
use m_xmpi

use m_io_tools, only : get_unit

implicit none

private

complex(dpc), public, parameter :: cmplx_i = (0.0_dp,1.0_dp)
complex(dpc), public, parameter :: cmplx_1 = (1.0_dp,0.0_dp)
complex(dpc), public, parameter :: cmplx_0 = (0.0_dp,0.0_dp)

logical, public  :: master_debug
character(len=100), public :: files_status_new='new'
character(len=100), public :: files_status_old='unknown'

public :: driver_invert_positive_definite_hermitian_matrix
public :: ritz_analysis_general, orthogonalize
public :: complex_vector_product
!!***
contains

!!****f* m_gwls_utility/complex_vector_product
!! NAME
!!  complex_vector_product
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

complex(dpc) function complex_vector_product(v1,v2,l)
!--------------------------------------------------------------------------
! This function computes the vector product of two complex vectors.
!--------------------------------------------------------------------------
implicit none

integer,     intent(in)  :: l
complex(dpc),intent(in)  :: v1(l), v2(l)

! *************************************************************************

complex_vector_product = sum(conjg(v1(:))*v2(:))

end function complex_vector_product
!!***

!!****f* m_gwls_utility/orthogonalize
!! NAME
!!  orthogonalize
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_GWlanczos
!!
!! CHILDREN
!!      matrix_function,xmpi_sum,zgemm
!!
!! SOURCE

subroutine orthogonalize(mpi_communicator, Hsize,Qsize,Xsize,Q,X)
!--------------------------------------------------------------------------
! This function wraps the relevant LAPACK routines to perform
!
!                X = (1 - Q.Q^H ) . X
!
! This essentially projects out  the subspace in Q from X.
! The array dimensions are:
!
!                                Q [ Hsize, Qsize]
!                                X [ Hsize, Xsize]
!                                Y [ Hsize, Xsize]
!
!  Hsize means "dimension of the Hilbert space", so typically the number
!  of plane waves...
!--------------------------------------------------------------------------
implicit none

integer,     intent(in)  :: mpi_communicator
integer,     intent(in)  :: Hsize, Qsize, Xsize
complex(dpc),intent(in)  :: Q(Hsize,Qsize)

complex(dpc),intent(inout)  :: X(Hsize,Xsize)

complex(dpc),allocatable :: C(:,:)

integer :: ierr

! *************************************************************************


ABI_ALLOCATE(C,(Qsize,Xsize))

! Compute Q^dagger . X
call ZGEMM(            'C',   & ! Hermitian conjugate the first array
'N',   & ! Leave second array as is
Qsize,   & ! the number of rows of the  matrix op( A )
Xsize,   & ! the number of columns of the  matrix op( B )
Hsize,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
cmplx_1,   & ! alpha constant
Q,   & ! matrix A
Hsize,   & ! LDA
X,   & ! matrix B
Hsize,   & ! LDB
cmplx_0,   & ! beta constant
C,   & ! matrix C
Qsize)     ! LDC


call xmpi_sum(C,mpi_communicator,ierr) ! sum on all processors working on FFT!

! Compute X - Q.(Q^dagger . X)
call ZGEMM(            'N',   & ! Leave first array as is
'N',   & ! Leave second array as is
Hsize,   & ! the number of rows of the  matrix op( A )
Xsize,   & ! the number of columns of the  matrix op( B )
Qsize,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
-cmplx_1,   & ! alpha constant
Q,   & ! matrix A
Hsize,   & ! LDA
C,   & ! matrix B
Qsize,   & ! LDB
cmplx_1,   & ! beta constant
X,   & ! matrix C
Hsize)     ! LDC

ABI_DEALLOCATE(C)

end subroutine orthogonalize
!!***

!!****f* m_gwls_utility/driver_invert_positive_definite_hermitian_matrix
!! NAME
!!  driver_invert_positive_definite_hermitian_matrix
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_GenerateEpsilon
!!
!! CHILDREN
!!      matrix_function,xmpi_sum,zgemm
!!
!! SOURCE

subroutine driver_invert_positive_definite_hermitian_matrix(matrix,ldim)
!----------------------------------------------------------------------------------------------------
!        This is a utility-type subroutine, which encapsulates the many Lapack steps necessary
!        to invert a positive definite hermitian matrix, and returns the full inverse to avoid
!        errors!
!
!        The subroutine overwrites the input.
!----------------------------------------------------------------------------------------------------
implicit none

integer     , intent(in)    :: ldim
complex(dpc), intent(inout) :: matrix(ldim,ldim)

! local variables
integer      :: i, j


integer      :: info

! *************************************************************************



! First, peform a decomposition
call zpotrf( 'U', ldim,matrix, ldim, info )

! Second, inverse the matrix in the new format
call zpotri( 'U', ldim,matrix, ldim, info )


! Finally,  properly symmetrise the matrix so that it is hermitian!
! The upper triangular part of the matrix is correct; the lower triangular must be built
do j=1, ldim
do i=j+1, ldim
matrix(i,j)= conjg(matrix(j,i))
end do
end do


end subroutine driver_invert_positive_definite_hermitian_matrix
!!***

!!****f* m_gwls_utility/ritz_analysis_general
!! NAME
!!  ritz_analysis_general
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_ComputePoles,gwls_GenerateEpsilon
!!
!! CHILDREN
!!      matrix_function,xmpi_sum,zgemm
!!
!! SOURCE

subroutine ritz_analysis_general(mpi_communicator,matrix_function,lmax,Hsize,Lbasis,eigenvalues)
!----------------------------------------------------------------------
!
! This subroutine is mostly for testing purposes.
!
! Given a matrix A_{NxN} which undergoes Lanczos analysis to yield
! a trigonal matrix T_{k x k}, for k lanczos steps, it is useful to
! test how well the eivenalues of T reproduce the eigenvalues of A.
!
! Following Matrix computations by Golub and Van Loan, define
!
!                Q = [ |  |  ...  | ],        T = S^H.D.S,  D diagonal, S unitary
!                    [ q1 q2 ...  qk]
!                    [ |  |  ...  | ]
!
!                Y = Q.S
!
! If tridiagonalisation was taken all the way to k = N, we would expect
! Y to contain the eigenvectors of A. It is useful to ask if, for a finite
! k, Y already contains good approximations to eigenvectors, by computing
! the Ritz residual,
!
!        Ri = || A.yi - di yi ||.
!
!        INPUTS:
!                     matrix_function    : the function which yields the action
!                                        of the implicit matrix on a vector
!                        lmax            : the total number of Lanczos steps
!                        Hsize           : the dimension of the Hilbert space
!                        Lbasis          : the Y matrix
!                       eigenvalues      : the computed approximate eigenvalues
!                                        of the matrix
!----------------------------------------------------------------------
implicit none

interface
  subroutine matrix_function(v_out,v_in,l)
  use defs_basis

  integer,     intent(in)  :: l
  complex(dp), intent(out) :: v_out(l)
  complex(dp), intent(in)  :: v_in(l)

  end subroutine matrix_function
end interface


integer,      intent(in)    :: Hsize, lmax , mpi_communicator
complex(dpc), intent(in)    :: Lbasis(Hsize,lmax)
real(dp),     intent(in)    :: eigenvalues(lmax)


! local variables
complex(dpc),allocatable :: check_matrix(:,:)
complex(dpc),allocatable :: yl(:), rl(:), Ayl(:)

real(dp)     :: lambda_l
real(dp)     :: check_norm
integer      :: l, i

real(dp)     :: norm_ritz, norm_ritz_squared

character(128) :: filename
logical        :: file_exists
integer        :: io_unit
integer        :: ierr
integer        :: mpi_rank

logical        :: head_node

! *************************************************************************

ABI_ALLOCATE(check_matrix,(lmax,lmax))
ABI_ALLOCATE(yl,(Hsize))
ABI_ALLOCATE(rl,(Hsize))
ABI_ALLOCATE(Ayl,(Hsize))

mpi_rank  = xmpi_comm_rank(mpi_communicator)

head_node = mpi_rank == 0

if (head_node) then
  io_unit = get_unit()

  i = 0
  file_exists = .true.
  do while (file_exists)
  i = i+1
  write(filename,'(A,I0.4,A)') "General_Ritz_Analisis_",i,".log"
  inquire(file=filename,exist=file_exists)
  end do


  open(io_unit,file=filename,status=files_status_new)

  write(io_unit,10) ''
  write(io_unit,10) '#===================================================================================================='
  write(io_unit,10) '# Entering ritz_analisis'
  write(io_unit,10) '# '
  write(io_unit,10) '#  parameters'
  write(io_unit,16) '#          Dimension of Hilbert space    : ',Hsize
  write(io_unit,16) '#          total number of Lanczos steps : ',lmax
  write(io_unit,10) '# '
  write(io_unit,10) '#===================================================================================================='
  write(io_unit,10) ''
  flush(io_unit)
end if

! NEVER use MATMUL! It stores temp arrays on stack, which kills executables compiled with intel!
! check_matrix(:,:) = matmul(transpose(conjg(Lbasis)),Lbasis)
call ZGEMM(            'C',   & ! Hermitian conjugate the first array
'N',   & ! Leave second array as is
lmax,   & ! the number of rows of the  matrix op( A )
lmax,   & ! the number of columns of the  matrix op( B )
Hsize,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
cmplx_1,   & ! alpha constant
Lbasis,   & ! matrix A
Hsize,   & ! LDA
Lbasis,   & ! matrix B
Hsize,   & ! LDB
cmplx_0,   & ! beta constant
check_matrix,   & ! matrix C
lmax)     ! LDC


call xmpi_sum(check_matrix,mpi_communicator,ierr) ! sum on all processors working on FFT!

do l = 1, lmax
check_matrix(l,l) =  check_matrix(l,l) - cmplx_1
end do
check_norm = sqrt(sum(abs(check_matrix(:,:))**2))


if (head_node) then
  write(io_unit,10) '#'
  write(io_unit,11) '#  Is the basis orthonormal?   || I - Y^H . Y || = ',check_norm
  write(io_unit,10) '#'
  flush(io_unit)


  ! Compute Ritz norms
  write(io_unit,10) ''
  write(io_unit,10) '#  Ritz analysis Lanczos Basis'
  write(io_unit,10) "#   l               lambda_l                       || R_l ||                                         "
  write(io_unit,10) '#===================================================================================================='
  flush(io_unit)
end if

do l = 1, lmax

lambda_l = eigenvalues(l)

yl       = Lbasis(:,l)
call matrix_function(Ayl,yl,Hsize)
rl       = Ayl - lambda_l*yl

norm_ritz_squared = sum(abs(rl(:))**2)
call xmpi_sum(norm_ritz_squared ,mpi_communicator,ierr)

norm_ritz = sqrt(norm_ritz_squared )


if (head_node) write(io_unit,14) l, lambda_l, norm_ritz

end do

if (head_node) then
  flush(io_unit)
  close(io_unit)
end if


ABI_DEALLOCATE(check_matrix)
ABI_DEALLOCATE(yl)
ABI_DEALLOCATE(rl)
ABI_DEALLOCATE(Ayl)


10 format(A)
11 format(A,ES24.16)
14 format(I5,2(5X,F24.12))
16 format(A,I5)


end subroutine ritz_analysis_general
!!***

end module m_gwls_utility

!!***

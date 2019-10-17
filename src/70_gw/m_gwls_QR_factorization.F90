!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gwls_QR_factorization
!! NAME
!! m_gwls_QR_factorization
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



module m_gwls_QR_factorization
!----------------------------------------------------------------------------------------------------
! This module implements the QR factorization using various algorithms, for the specific
! data distribution corresponding to FFT parallelism.
!
! There are standard routines which do this (lapack, scalapack, etc), however it is complicated
! to get scalapack to run properly in parallel. Implementing ourselves is the shortest path to
! a working solution.
!----------------------------------------------------------------------------------------------------
!local modules
use m_gwls_utility
use m_gwls_TimingLog
use m_gwls_wf
use m_gwls_hamiltonian

!abinit modules
use defs_basis
use defs_abitypes
use defs_wvltypes
use m_abicore
use m_xmpi
use m_errors

use defs_abitypes, only : MPI_type
use m_io_tools,  only : get_unit
use m_time,      only : timab


implicit none
save
private
!!***

logical, private ::  debug = .false.
!!***

public :: extract_QR, extract_SVD
!!***

contains

!!****f* m_hamiltonian/extract_QR
!! NAME
!!  extract_QR
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_GWlanczos,gwls_QR_factorization
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_sum
!!
!! SOURCE

subroutine extract_QR(mpi_communicator,Hsize,Xsize,Xmatrix,Rmatrix)
!--------------------------------------------------------------------------
! This function computes the QR factorization:
!
!                X = Q . R
!
! in order to extract the matrix of orthonormal vectors Q and
! the R matrix.
!
! On output, the matrix X is replaced by Q.
!
! If the code is running with only one processor, this routine
! simply invokes extract_QR_serial, which wraps standard Lapack routines.
! If we are running in MPI parallel, the serial Lapack routines no
! longer work (and understanding scalapack is too complicated right now).
! Thus, in that case, this routine implements some old school Gram-Schmidt
! algorithm.
!--------------------------------------------------------------------------
implicit none

integer,        intent(in) :: Hsize, Xsize, mpi_communicator
complex(dpc),intent(inout) :: Xmatrix(Hsize,Xsize)

complex(dpc),  intent(out),optional :: Rmatrix(Xsize,Xsize)

! local variables

real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB

! *************************************************************************



GWLS_TIMAB   = 1519
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


!--------------------------------------------------------------------------------
! Implement Gram-Schmidt.
!--------------------------------------------------------------------------------
call extract_QR_Householder(mpi_communicator,Hsize,Xsize,Xmatrix,Rmatrix)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

end subroutine extract_QR
!!***


!!****f* m_hamiltonian/extract_SVD
!! NAME
!!  extract_SVD
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_sum
!!
!! SOURCE

subroutine extract_SVD(mpi_communicator, Hsize,lsolutions_max,svd_matrix,svd_values)
!--------------------------------------------------------------------------
! This function computes the singular value decomposition
!
!                X = U . SIGMA . V^dagger
!
! More specifically,  the matrix U of orthonormal vectors and SIGMA
! the eigenvalues are returned.
!
! different algorithms are used, depending on parallelisation scheme.
!--------------------------------------------------------------------------
implicit none

integer,      intent(in)    :: mpi_communicator
integer,      intent(in)    :: Hsize, lsolutions_max
complex(dpc), intent(inout) :: svd_matrix(Hsize,lsolutions_max)
real(dp),     intent(out)   :: svd_values(lsolutions_max)


complex(dpc), allocatable   :: Rmatrix(:,:)
complex(dpc), allocatable   :: svd_tmp(:,:)

real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB

! *************************************************************************



GWLS_TIMAB   = 1520
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


!if ( mpi_enreg%nproc_fft ==1 ) then
if ( .false. ) then

  call extract_SVD_lapack(Hsize,lsolutions_max,svd_matrix,svd_values)

else

  ABI_ALLOCATE(Rmatrix,(lsolutions_max,lsolutions_max))

  ! perform QR first
  call extract_QR(mpi_communicator, Hsize,lsolutions_max,svd_matrix,Rmatrix)

  ! perform SVD on the much smaller Rmatrix!
  call extract_SVD_lapack(lsolutions_max,lsolutions_max,Rmatrix,svd_values)

  ABI_ALLOCATE(svd_tmp,(Hsize,lsolutions_max))

  ! Rmatrix is overwritten with U matrix from SVD. Update the svd_matrix
  call ZGEMM(            'N',   & ! Leave first array as is
  'N',   & ! Leave second array as is
  Hsize,   & ! the number of rows of the  matrix op( A )
  lsolutions_max,   & ! the number of columns of the  matrix op( B )
  lsolutions_max,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
  cmplx_1,   & ! alpha constant
  svd_matrix,   & ! matrix A
  Hsize,   & ! LDA
  Rmatrix,   & ! matrix B
  lsolutions_max,   & ! LDB
  cmplx_0,   & ! beta constant
  svd_tmp,   & ! matrix C
  Hsize)     ! LDC

  svd_matrix(:,:) = svd_tmp(:,:)

  ABI_DEALLOCATE(svd_tmp)
  ABI_DEALLOCATE(Rmatrix)
end if
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


end subroutine extract_SVD
!!***

!!****f* m_hamiltonian/extract_SVD_lapack
!! NAME
!!  extract_SVD_lapack
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_QR_factorization
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_sum
!!
!! SOURCE

subroutine extract_SVD_lapack(Hsize,lsolutions_max,svd_matrix,svd_values)
!--------------------------------------------------------------------------
! This function computes the singular value decomposition
! using lapack routines. This is not appropriate in MPI parallel!
!
!
!--------------------------------------------------------------------------
implicit none


integer,      intent(in)    :: Hsize, lsolutions_max
complex(dpc), intent(inout) :: svd_matrix(Hsize,lsolutions_max)
real(dp),     intent(out)   :: svd_values(lsolutions_max)



integer                   :: info_zgesvd
integer                   :: lwork_svd
complex(dpc), allocatable :: work_svd(:)
complex(dpc), allocatable :: svd_U(:,:), svd_V(:,:)
real   (dp ), allocatable :: rwork_svd(:)

integer        :: debug_unit
character(50)  :: debug_filename

! *************************************************************************




! allocate arrays for the svd
ABI_ALLOCATE(svd_U                  ,(1,1))
ABI_ALLOCATE(svd_V                  ,(1,1))


! DIMENSION QUERRY for singluar decomposition problem

ABI_ALLOCATE(rwork_svd    ,(5*min(Hsize,lsolutions_max)))
ABI_ALLOCATE(work_svd,(1))
lwork_svd = -1

call zgesvd('O',            & ! The first min(m,n) columns of U (the left singular vectors) are overwritten on the array A;
'N',            & ! no column vectors of V are computed
Hsize,          & ! number of rows of the matrix
lsolutions_max, & ! number of columns of the matrix
svd_matrix,     & ! matrix to be decomposed
Hsize,          & ! LDA
svd_values,     & ! singular values
svd_U,          & ! dummy U; not referenced
1,              & ! size of U
svd_V,          & ! dummy V; not referenced
1,              & ! size of V
work_svd,       & ! work array
lwork_svd,      & ! size of work array
rwork_svd,      & ! work array
info_zgesvd )

if ( info_zgesvd /= 0) then
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info_zgesvd,' in ZGESVD(1), gwls_QR_factorization'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if






lwork_svd = nint(dble(work_svd(1)))

ABI_DEALLOCATE(work_svd)

ABI_ALLOCATE(work_svd,(lwork_svd))

! computation run

call zgesvd('O',            & ! The first min(m,n) columns of U (the left singular vectors) are overwritten on the array A;
'N',            & ! no column vectors of V are computed
Hsize,          & ! number of rows of the matrix
lsolutions_max, & ! number of columns of the matrix
svd_matrix,     & ! matrix to be decomposed
Hsize,          & ! LDA
svd_values,     & ! singular values
svd_U,          & ! dummy U; not referenced
1,              & ! size of U
svd_V,          & ! dummy V; not referenced
1,              & ! size of V
work_svd,       & ! work array
lwork_svd,      & ! size of work array
rwork_svd,      & ! work array
info_zgesvd )

if ( info_zgesvd /= 0) then
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info_zgesvd,' in ZGESVD(2), gwls_QR_factorization'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if




ABI_DEALLOCATE(work_svd)
ABI_DEALLOCATE(rwork_svd)
ABI_DEALLOCATE(svd_U)
ABI_DEALLOCATE(svd_V)

end subroutine extract_SVD_lapack
!!***





!!****f* m_hamiltonian/extract_QR_Householder
!! NAME
!!  extract_QR_Householder
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_QR_factorization
!!
!! CHILDREN
!!      xmpi_allgather,xmpi_sum
!!
!! SOURCE

subroutine extract_QR_Householder(mpi_communicator,Hsize,Xsize,Xmatrix,Rmatrix)
!--------------------------------------------------------------------------
! This function computes the QR factorization:
!
!                X = Q . R
!
! in order to extract the matrix of orthonormal vectors Q and
! the R matrix.
!
! On output, the matrix X is replaced by Q.
!
! This routine uses Householder operations to generate Q and R.
! Special attention is given to the fact that the matrix may be
! distributed across processors and MPI communication is necessary.
!
!--------------------------------------------------------------------------
implicit none

integer,        intent(in) :: Hsize, Xsize, mpi_communicator
complex(dpc),intent(inout) :: Xmatrix(Hsize,Xsize)

complex(dpc),  intent(out),optional :: Rmatrix(Xsize,Xsize)

! local variables
integer        :: numbrer_of_plane_waves

integer        :: io_unit
integer        :: ierr
character(50)  :: filename
logical        :: file_exists

integer,save :: counter = 0

integer :: i, j, l_local
integer :: l1, l2

integer, allocatable      :: nproc_array(:)

complex(dpc), allocatable :: Qinternal(:,:)
complex(dpc), allocatable :: Rinternal(:,:)
complex(dpc), allocatable :: vj(:)

complex(dpc), allocatable :: A_matrix(:,:)
complex(dpc), allocatable :: V_matrix(:,:)

complex(dpc), allocatable :: list_beta(:)

complex(dpc) :: cmplx_value
real   (dp ) :: real_value

complex(dpc) :: norm_x
complex(dpc) :: phase

complex(dpc), allocatable :: error(:,:)
complex(dpc), allocatable :: coeff(:)


integer :: mpi_rank
integer :: mpi_nproc
logical :: head_node

! *************************************************************************


!--------------------------------------------------------------------------------
! Implement Householder algorithm, in parallel
!--------------------------------------------------------------------------------
mpi_nproc        = xmpi_comm_size(mpi_communicator)

! extract the rank of this processor in the communicator
mpi_rank         = xmpi_comm_rank(mpi_communicator)

! only head node will write!
head_node        = mpi_rank == 0


!--------------------------------------------------------------------------------
! Open a log file for the output of extract_QR
!--------------------------------------------------------------------------------
if (debug .and.  head_node ) then

  io_unit  = get_unit()
  write(filename,'(A,I0.4,A)') "extract_QR_",mpi_rank,".log"
  inquire(file=trim(filename),exist=file_exists)

  if (file_exists) then
    open(io_unit,file=trim(filename),position='append',status=files_status_old)
  else
    open(io_unit,file=trim(filename),status=files_status_new)
    write(io_unit,10) "#======================================================================================="
    write(io_unit,10) "#                                                                                       "
    write(io_unit,10) "#   This file contains information regarding the QR factorization, from extract_QR      "
    write(io_unit,25) "#   The algorithm is running in MPI parallel with ",mpi_nproc," processors"
    write(io_unit,10) "#                                                                                       "
    write(io_unit,10) "#======================================================================================="
  end if

  counter = counter + 1

  write(io_unit,10) "#                                                                                       "
  write(io_unit,11) "#   Call # ", counter
  write(io_unit,10) "#                                                                                       "
  write(io_unit,11) "#                        Hsize    = ",Hsize
  write(io_unit,11) "#                        Xsize    = ",Xsize
  write(io_unit,13) "#                Rmatrix present? = ",present(Rmatrix)


end if

!--------------------------------------------------------------------------------
! Get the number of plane waves on every processor
!--------------------------------------------------------------------------------
ABI_ALLOCATE(nproc_array,(mpi_nproc))

nproc_array = 0

numbrer_of_plane_waves = Hsize ! do this to avoid "intent" problems
call xmpi_allgather(numbrer_of_plane_waves, nproc_array, mpi_communicator, ierr)


!--------------------------------------------------------------------------------
! Get the offset for each processor
!
! The global index is then given by
!                   I_{global} = nproc_array(1+rank)+i_{local}
!
! similarly, the local index is given by
!                  i_{local} = I_{global} - nproc_array(1+rank)
! which is only meaningful if 1 <= i_{local} <= Hsize
!--------------------------------------------------------------------------------

do j = mpi_nproc, 1, -1
nproc_array(j) = 0
do i = 1, j-1
nproc_array(j) = nproc_array(j)+nproc_array(i)
end do
end do


!--------------------------------------------------------------------------------
! Act on the A matrix, following the book by Golub  (more or less ;) )
!
!--------------------------------------------------------------------------------

ABI_ALLOCATE(A_matrix, (Hsize,Xsize))
ABI_ALLOCATE(V_matrix, (Hsize,Xsize))
ABI_ALLOCATE(list_beta, (Xsize))
ABI_ALLOCATE(coeff   , (Xsize))

A_matrix(:,:) = Xmatrix(:,:)
V_matrix(:,:) = cmplx_0
list_beta(:)  = cmplx_0


ABI_ALLOCATE(vj, (Hsize))

do j = 1, Xsize

! Store xj in vj, for now
vj(:) = A_matrix(:,j)


if (j > 1) then
  !------------------------------------------
  ! set the array to zero all the way to j-1
  !------------------------------------------
  l_local = j-1-nproc_array(1+mpi_rank)

  if ( l_local > Hsize) then
    vj(:) = cmplx_0
  else if ( l_local <= Hsize  .and. l_local >= 1) then
    vj(1:l_local) = cmplx_0
  end if

end if

! compute the norm of x
norm_x = sum(conjg(vj(:))*vj(:))
call xmpi_sum(norm_x,mpi_communicator,ierr) ! sum on all processors in communicator
norm_x = sqrt(norm_x)


if (abs(norm_x) > tol14) then
  ! if |x| ~ 0, there is nothing to do! the column in A is full of zeros.

  ! find the j^th component of x
  l_local = j-nproc_array(1+mpi_rank)

  ! update vj, on the right processor!
  if ( l_local <= Hsize  .and. l_local >= 1) then

    phase = vj(l_local)

    if (abs(phase) < tol14) then
      phase = cmplx_1
    else
      phase = phase/abs(phase)
    end if

    vj(l_local) = vj(l_local) + phase*norm_x

  end if

  !compute beta
  cmplx_value = sum(conjg(vj(:))*vj(:))
  call xmpi_sum(cmplx_value,mpi_communicator,ierr) ! sum on all processors

  list_beta(j) = 2.0_dp/cmplx_value

  ! store v for later use; this is less efficient than storing it in the null part of A,
  ! but it is less of a pain to implement. Feel free to clean this up.
  V_matrix(:,j) = vj(:)


  ! Update the A matrix

  ! Compute v^dagger . A
  !call ZGEMM(            'C',   & ! Hermitian conjugate the first array
  !                       'N',   & ! Leave second array as is
  !                         1,   & ! the number of rows of the  matrix op( A )
  !                 Xsize-j+1,   & ! the number of columns of the  matrix op( B )
  !                     Hsize,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
  !                   cmplx_1,   & ! alpha constant
  !                        vj,   & ! matrix A
  !                     Hsize,   & ! LDA
  !       A_matrix(:,j:Xsize),   & ! matrix B
  !                     Hsize,   & ! LDB
  !                   cmplx_0,   & ! beta constant
  !          coeff(:,j:Xsize),   & ! matrix C
  !                         1)     ! LDC
  do i = j, Xsize
  coeff(i) = sum(conjg(vj)*A_matrix(:,i))
  end do
  call xmpi_sum(coeff,mpi_communicator,ierr) ! sum on all processors in the communicator

  ! update A
  do i = j, Xsize
  A_matrix(:,i) = A_matrix(:,i) - list_beta(j)*coeff(i)*vj(:)
  end do

end if

end do

!--------------------------------------------------------------------------------
! Extract the R matrix
!
!--------------------------------------------------------------------------------

ABI_ALLOCATE(Rinternal,(Xsize,Xsize))

Rinternal = cmplx_0


do j = 1, Xsize
do i = 1, j
l_local = i-nproc_array(1+mpi_rank)

if ( l_local <= Hsize  .and. l_local >= 1) then
  Rinternal(i,j) = A_matrix(l_local,j)
end if

end do ! i
end do ! j

call xmpi_sum(Rinternal,mpi_communicator,ierr) ! sum on all processors



!--------------------------------------------------------------------------------
! Extract the Q matrix
!
!--------------------------------------------------------------------------------

ABI_ALLOCATE( Qinternal, (Hsize,Xsize))

! initialize Q to the identity in the top corner
Qinternal = cmplx_0

do j = 1, Xsize
l_local = j-nproc_array(1+mpi_rank)

if ( l_local <= Hsize  .and. l_local >= 1 ) then
  Qinternal(l_local,j) = cmplx_1
end if

end do ! j


! Build Q interatively
do j = Xsize,1, -1

vj(:) = V_matrix(:,j)


! Update the A matrix

! Compute v^dagger . A
!call ZGEMM(            'C',   & ! Hermitian conjugate the first array
!                       'N',   & ! Leave second array as is
!                         1,   & ! the number of rows of the  matrix op( A )
!                     Xsize,   & ! the number of columns of the  matrix op( B )
!                     Hsize,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
!                   cmplx_1,   & ! alpha constant
!                        vj,   & ! matrix A
!                     Hsize,   & ! LDA
!                 Qinternal,   & ! matrix B
!                     Hsize,   & ! LDB
!                   cmplx_0,   & ! beta constant
!                     coeff,   & ! matrix C
!                         1)     ! LDC

do i = 1, Xsize
coeff(i) = sum(conjg(vj)*Qinternal(:,i))
end do
call xmpi_sum(coeff,mpi_communicator,ierr) ! sum on all processors in communicator


! update Q
do i = 1, Xsize
Qinternal(:,i) = Qinternal(:,i) - list_beta(j)*coeff(i)*vj(:)
end do

end do ! j



! clean up
ABI_DEALLOCATE(V_matrix)
ABI_DEALLOCATE(coeff)
ABI_DEALLOCATE(vj)

!--------------------------------------------------------------------------------
! Do some debug, if requested
!
!--------------------------------------------------------------------------------

if (debug ) then

  if ( head_node ) then

    write(io_unit,20) "#    nproc_array   = ",nproc_array
    flush(io_unit)

    write(io_unit,40) "#    list_beta     = ",real(list_beta)
    flush(io_unit)
  end if


  ABI_ALLOCATE(error,(Xsize,Xsize))

  error = cmplx_0

  do l2=1,Xsize
  error(l2,l2) = error(l2,l2) - cmplx_1
  do l1=1,Xsize

  cmplx_value = complex_vector_product(Qinternal(:,l1),Qinternal(:,l2),Hsize)
  call xmpi_sum(cmplx_value,mpi_communicator,ierr) ! sum on all processors working on FFT!

  error(l1,l2) = error(l1,l2)+cmplx_value

  end do
  end do

  if ( head_node ) then
    write(io_unit,12) "#               || Q^t.Q - I ||   = ",sqrt(sum(abs(error(:,:))**2))
    flush(io_unit)
  end if


  ABI_DEALLOCATE(error)

  ABI_ALLOCATE(error,(Hsize,Xsize))

  error = Xmatrix

  do l2=1,Xsize
  do l1=1,Xsize
  error(:,l2) = error(:,l2) - Qinternal(:,l1)*Rinternal(l1,l2)
  end do
  end do

  real_value = zero
  do l2=1,Xsize
  do l1=1,Xsize
  cmplx_value = complex_vector_product(error(:,l1),error(:,l2),Hsize)

  real_value  = real_value + abs(cmplx_value)**2
  end do
  end do

  call xmpi_sum(real_value,mpi_communicator,ierr) ! sum on all processors


  real_value = sqrt(real_value)

  if ( head_node) then
    write(io_unit,12) "#               || Xin - Q.R ||   = ",real_value

    if ( real_value > 1.0D-10 ) write(io_unit,10) "#               ERROR!              "

    write(io_unit,10) '#'
    write(io_unit,10) '# R matrix'
    write(io_unit,10) '#'
    do l1=1, Xsize
    write(io_unit,30) Rinternal(l1,:)
    end do
    write(io_unit,10) ''
    write(io_unit,10) ''


    write(io_unit,10) '#'
    write(io_unit,10) '# top of A matrix'
    write(io_unit,10) '#'
    do l1=1, 2*Xsize
    write(io_unit,30) A_matrix(l1,:)
    end do
    write(io_unit,10) ''
    write(io_unit,10) ''

    close(io_unit)

  end if

  ABI_DEALLOCATE(error)

end if

!--------------------------------------------------------------------------------
! Final assignment
!
!--------------------------------------------------------------------------------

Xmatrix = Qinternal

if (present(Rmatrix)) then
  Rmatrix = Rinternal
end if

ABI_DEALLOCATE(Qinternal)
ABI_DEALLOCATE(Rinternal)
ABI_DEALLOCATE(nproc_array)
ABI_DEALLOCATE(A_matrix)
ABI_DEALLOCATE(list_beta)


10 format(A)
11 format(A,I8)
12 format(A,E24.16)
13 format(A,2X,L10)
20 format(A,1000I10)
25 format(A,I5,A)
30 format(1000(F22.16,2X,F22.16,5X))
40 format(A,1000(F22.16,5X))

end subroutine extract_QR_Householder
!!***




end module m_gwls_QR_factorization
!!***

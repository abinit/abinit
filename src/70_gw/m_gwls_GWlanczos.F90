!!****m* ABINIT/m_gwls_GWlanczos
!! NAME
!! m_gwls_GWlanczos
!!
!! FUNCTION
!!  .
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (JLJ, BR, MC)
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


module m_gwls_GWlanczos
!----------------------------------------------------------------------------------------------------
! This module implements the Lanczos scheme to band diagonalize an implicit operator.
!----------------------------------------------------------------------------------------------------
!local modules
use m_gwls_utility
use m_gwls_TimingLog
use m_gwls_wf
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_polarisability
use m_gwls_QR_factorization

!abinit modules
use defs_basis
use defs_wvltypes
use m_abicore
use m_xmpi
use m_pawang
use m_errors

use m_io_tools,         only : get_unit
use m_paw_dmft,         only : paw_dmft_type
use m_ebands,           only : ebands_init, ebands_free

implicit none
save
private
!!***

!!***
public :: block_lanczos_algorithm

public :: diagonalize_lanczos_banded
public :: get_seeds
!!***
contains


!!****f* m_gwls_GWlanczos/get_seeds
!! NAME
!!  get_seeds
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputePoles,m_gwls_GenerateEpsilon
!!
!! CHILDREN
!!      zgemm,zhbev
!!
!! SOURCE

subroutine get_seeds(first_seed, nseeds, seeds)
!----------------------------------------------------------------------------------------------------
!
! This subroutine compute the seeds using the eigenstates of the Hamiltonian
!
!----------------------------------------------------------------------------------------------------
implicit none

integer,      intent(in)  :: first_seed, nseeds
complex(dpc), intent(out) :: seeds(npw_k,nseeds)

real(dp)    , allocatable :: psik_out(:,:)
real(dp)    , allocatable :: psikb_e(:,:)
real(dp)    , allocatable :: psig_e(:,:)
real(dp)    , allocatable :: psikb_s(:,:)
real(dp)    , allocatable :: psig_s(:,:)



! local variables
integer  :: n
integer  :: i, j, nsblk

! *************************************************************************

! Generate the seeds for the Lanczos algorithm
ABI_MALLOC(psik_out,(2,npw_k))
ABI_MALLOC(psikb_e,(2,npw_kb))
ABI_MALLOC(psig_e,(2,npw_g))
ABI_MALLOC(psikb_s,(2,npw_kb))
ABI_MALLOC(psig_s,(2,npw_g))

nsblk = ceiling(1.0*nseeds/blocksize)

do i=1,nsblk
do j=1,blocksize
psikb_e(:,(j-1)*npw_k+1:j*npw_k) = cg(:,(e-1)*npw_k+1:e*npw_k)
end do

psig_e = zero
call wf_block_distribute(psikb_e,  psig_e,1) ! LA -> FFT

do j=1,blocksize
n = (i-1)*blocksize + j-1 + first_seed
if ((i-1)*blocksize + j <= nseeds) then
  psikb_s(:,(j-1)*npw_k+1:j*npw_k) = cg(:,(n-1)*npw_k+1:n*npw_k)
else
  psikb_s(:,(j-1)*npw_k+1:j*npw_k) = zero
end if
end do

psig_s = zero
call wf_block_distribute(psikb_s,  psig_s,1) ! LA -> FFT

! Fourier transform valence wavefunction, to real space
call g_to_r(psir1,psig_s)

psir1(2,:,:,:) = -psir1(2,:,:,:)

call gr_to_g(psig_s,psir1,psig_e)

! return to LA configuration, in order to apply Coulomb potential
call wf_block_distribute(psikb_s,  psig_s,2) ! FFT -> LA

do j=1, blocksize
n = (i-1)*blocksize + j
if(n<=nseeds) then
  psik_out = psikb_s(:,(j-1)*npw_k+1:j*npw_k)
  call sqrt_vc_k(psik_out)
  seeds(:,n) = cmplx_1*psik_out(1,:) + cmplx_i*psik_out(2,:)
end if
end do
end do

ABI_FREE(psik_out)
ABI_FREE(psikb_e)
ABI_FREE(psig_e)
ABI_FREE(psikb_s)
ABI_FREE(psig_s)

end subroutine get_seeds
!!***

!!****f* m_gwls_GWlanczos/block_lanczos_algorithm
!! NAME
!!  block_lanczos_algorithm
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputePoles,m_gwls_GenerateEpsilon,m_gwls_LanczosResolvents
!!
!! CHILDREN
!!      zgemm,zhbev
!!
!! SOURCE

subroutine block_lanczos_algorithm(mpi_communicator,matrix_function,kmax,nseeds,Hsize,seeds,alpha,beta,Lbasis,X0,beta0,Qk)
!----------------------------------------------------------------------------------------------------
!
! This subroutine implements the Block Lanczos algorithm for an arbitrary, user-supplied function which
! returns the action of the implicit matrix, which is assumed to be Hermitian.
!
!
!----------------------------------------------------------------------------------------------------
implicit none

!-----------------------------------------
! interface with implicit matrix function
!-----------------------------------------
!********************************************************
!***        NOTE:                                           ***
!***                                                        ***
!***        There *appears* to be a bug in makemake;        ***
!***        when the name of the function in the interface  ***
!***        below contains the word "implicit", makemake    ***
!***        yields an error. I suppose makemake simply      ***
!***        parses for the word "implicit" without regards  ***
!***        for the fact that it may simply be part of a    ***
!***        naive developper's chosen name for the routine. ***
!***        Correspondingly, I change the name to           ***
!***        "matrix_function" to avoid problems.            ***
!***                                                        ***
!***                                     Bruno Rousseau     ***
!***                                     08/13/2012         ***
!********************************************************
interface
  subroutine matrix_function(v_out,v_in,l)

  use defs_basis

  integer,     intent(in)   :: l
  complex(dpc), intent(out) :: v_out(l)
  complex(dpc), intent(in)  :: v_in(l)

  end subroutine matrix_function
end interface

!------------------------------
! input/output variables
!------------------------------

integer, intent(in) :: mpi_communicator
integer, intent(in) :: kmax        ! number of Lanczos blocks
integer, intent(in) :: nseeds      ! size of each blocks
integer, intent(in) :: Hsize       ! size of the Hilbert space in which the matrix lives

complex(dpc), intent(inout):: seeds(Hsize,nseeds) ! seed vectors for the algorithm
! overwritten by X_{k+1} on output

!logical,      intent(in) :: ortho           ! should the Lanczos vector be orthogonalized?

complex(dpc), intent(out) :: alpha(nseeds,nseeds,kmax)  ! the alpha array from the Lanczos algorithm
complex(dpc), intent(out) :: beta(nseeds,nseeds,kmax)   ! the  beta array from the Lanczos algorithm
complex(dpc), intent(out) :: Lbasis(Hsize,nseeds*kmax)  ! array containing the Lanczos basis


complex(dpc), intent(in),optional :: X0(Hsize,nseeds)
complex(dpc), intent(in),optional :: beta0(nseeds,nseeds)
complex(dpc), intent(in),optional :: Qk(:,:)  ! array containing vectors to which

! the basis must be orthonormalized



!------------------------------
! local variables
!------------------------------
integer     :: k, seed1
integer     :: dum(2), lk

complex(dpc), allocatable :: xk(:,:), xkm1(:,:), rk(:,:)

integer     :: ntime, itime
real(dp)    :: total_time1, total_time2
real(dp)    :: time1, time2
integer     :: ierr
real(dp),allocatable :: list_time(:)

! *************************************************************************

call cpu_time(total_time1)


ntime = 7
ABI_MALLOC(list_time,(ntime))
list_time(:) = zero

if(present(Qk)) then
  dum   = shape(Qk)
  lk    = dum(2)
end if

ABI_MALLOC( xk,  (Hsize,nseeds))
ABI_MALLOC( xkm1,(Hsize,nseeds))
ABI_MALLOC( rk  ,(Hsize,nseeds))


alpha = cmplx_0
beta  = cmplx_0

!------------------------------------------------
! Orthonormalize the seeds
!------------------------------------------------
! initialize the xk array with the seeds
xk(:,:) = seeds(:,:)


! orthonormalize the block using the QR algorithm
! xk is overwritten by Q, the array of orthonormal vectors

call extract_QR(mpi_communicator, Hsize,nseeds,xk)

!------------------------------------------------
! Loop on all blocks
!------------------------------------------------


do k = 1, kmax

itime = 0

! tabulate basis, computed at previous step
Lbasis(:,nseeds*(k-1)+1:nseeds*k) = xk(:,:)


itime = itime+1
call cpu_time(time1)

! Initialize the residual array
do seed1 = 1, nseeds
! If we are constructing the $\hat \epsilon(i\omega = 0)$ matrix (and the Lanczos basis at the same time),
! note the index in which the Sternheimer solutions will be stored (for use in the projected Sternheimer section).
if(write_solution) index_solution = (k-1)*nseeds + seed1

call matrix_function(rk(:,seed1),xk(:,seed1),Hsize)
end do

call cpu_time(time2)
list_time(itime) = list_time(itime) + time2-time1

itime = itime+1
call cpu_time(time1)
! compute the alpha array, alpha = X^d.A.X

call ZGEMM(              'C',   & ! take Hermitian conjugate of first array
'N',   & ! leave second array as is
nseeds,   & ! the number  of rows of the  matrix op( A )
nseeds,   & ! the number  of columns of the  matrix op( B )
Hsize,   & ! the number  of columns of the  matrix op( A ) == rows of matrix op( B )
cmplx_1,   & ! alpha constant
xk,   & ! matrix A
Hsize,   & ! LDA
rk,   & ! matrix B
Hsize,   & ! LDB
cmplx_0,   & ! beta constant
alpha(:,:,k),   & ! matrix C
nseeds)     ! LDC
call xmpi_sum(alpha(:,:,k),mpi_communicator,ierr) ! sum on all processors



call cpu_time(time2)
list_time(itime) = list_time(itime) + time2-time1

itime = itime+1
call cpu_time(time1)
! update the residual array, rk = rk-X.alpha
call ZGEMM(              'N',   & ! leave first array as is
'N',   & ! leave second array as is
Hsize,   & ! the number  of rows of the  matrix op( A )
nseeds,   & ! the number  of columns of the  matrix op( B )
nseeds,   & ! the number  of columns of the  matrix op( A ) == rows of matrix op( B )
-cmplx_1,   & ! alpha constant
xk,   & ! matrix A
Hsize,   & ! LDA
alpha(:,:,k),   & ! matrix B
nseeds,   & ! LDB
cmplx_1,   & ! beta constant
rk,   & ! matrix C
Hsize)     ! LDC

call cpu_time(time2)
list_time(itime) = list_time(itime) + time2-time1

if (k .eq. 1 .and. present(X0) .and. present(beta0)) then
  ! if k == 1, and X0,beta0 are present,
  !  update the residual array, r1 = r1-X_{0}.beta^d_{0}
  call ZGEMM(                'N',   & ! leave first array as is
  'C',   & ! Hermitian conjugate the second array
  Hsize,   & ! the number  of rows of the  matrix op( A )
  nseeds,   & ! the number  of columns of the  matrix op( B )
  nseeds,   & ! the number  of columns of the  matrix op( A ) == rows of matrix op( B )
  -cmplx_1,   & ! alpha constant
  X0,   & ! matrix A
  Hsize,   & ! LDA
  beta0(:,:),   & ! matrix B
  nseeds,   & ! LDB
  cmplx_1,   & ! beta constant
  rk,   & ! matrix C
  Hsize)     ! LDC
end if


itime = itime+1
call cpu_time(time1)
if (k .gt. 1) then

  ! if k > 1, update the residual array, rk = rk-X_{k-1}.beta^d_{k-1}
  call ZGEMM(                'N',   & ! leave first array as is
  'C',   & ! Hermitian conjugate the second array
  Hsize,   & ! the number  of rows of the  matrix op( A )
  nseeds,   & ! the number  of columns of the  matrix op( B )
  nseeds,   & ! the number  of columns of the  matrix op( A ) == rows of matrix op( B )
  -cmplx_1,   & ! alpha constant
  xkm1,   & ! matrix A
  Hsize,   & ! LDA
  beta(:,:,k-1),   & ! matrix B
  nseeds,   & ! LDB
  cmplx_1,   & ! beta constant
  rk,   & ! matrix C
  Hsize)     ! LDC

end if
call cpu_time(time2)
list_time(itime) = list_time(itime) + time2-time1

! store xk for next iteration
xkm1(:,:) = xk(:,:)


! Orthonormalize THE RESIDUAL to all previously calculated directions

itime = itime+1
call cpu_time(time1)
!if ( ortho .and. (dtset%gwcalctyp/=1) ) then !This is a test to obtain the CPU time taken by the orthogonalizations.

if(present(Qk)) then
  ! Orthonormalize to all previously calculated directions, if
  ! this is a restarted Lanczos step
  call orthogonalize(mpi_communicator, Hsize,lk,nseeds,Qk,rk)
end if

call orthogonalize(mpi_communicator, Hsize,k*nseeds,nseeds,Lbasis(:,1:k*nseeds),rk)

!end if
call cpu_time(time2)
list_time(itime) = list_time(itime) + time2-time1


itime = itime+1
call cpu_time(time1)

! perform QR decomposition to extract X_{k+1} and beta_{k}
call extract_QR(mpi_communicator, Hsize,nseeds,rk,beta(:,:,k))

call cpu_time(time2)
list_time(itime) = list_time(itime) + time2-time1
! copy the Q matrix (written on rk) in xk, which becomes X_{k+1}
xk(:,:) = rk(:,:)


end do !end loop on k

! overwrite the seeds with the last vector block.
seeds(:,:) = xk(:,:)

ABI_FREE( xk  )
ABI_FREE( xkm1)
ABI_FREE( rk  )
call cpu_time(total_time2)

list_time(7) = total_time2-total_time1

call write_block_lanczos_timing_log(list_time,ntime)

ABI_FREE(list_time)

end subroutine block_lanczos_algorithm
!!***

!!****f* m_gwls_GWlanczos/diagonalize_lanczos_banded
!! NAME
!!  diagonalize_lanczos_banded
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputePoles,m_gwls_GenerateEpsilon,m_gwls_LanczosResolvents
!!
!! CHILDREN
!!      zgemm,zhbev
!!
!! SOURCE

subroutine diagonalize_lanczos_banded(kmax,nseeds,Hsize,alpha,beta,Lbasis,eigenvalues,debug)
!-----------------------------------------------------------------------------------
! Given the result of the Lanczos algorithm, this subroutine diagonalize the banded
! matrix as well as updates the basis.
!-----------------------------------------------------------------------------------
implicit none

integer, intent(in)  :: kmax        ! number of Lanczos blocks
integer, intent(in)  :: nseeds      ! size of each blocks
integer, intent(in)  :: Hsize       ! size of the Hilbert space in which the matrix lives
logical, intent(in)  :: debug

complex(dpc), intent(in) :: alpha(nseeds,nseeds,kmax)  ! the alpha array from the Lanczos algorithm
complex(dpc), intent(in) :: beta (nseeds,nseeds,kmax)  ! the  beta array from the Lanczos algorithm

complex(dpc), intent(inout) :: Lbasis(Hsize,nseeds*kmax)  ! array containing the Lanczos basis


real(dp), intent(out) :: eigenvalues(nseeds*kmax)


! local variables

integer :: kd   ! number of superdiagonal above the diagonal in banded storage
integer :: ldab ! dimension of banded storage matrix

complex(dpc), allocatable :: band_storage_matrix(:,:)
complex(dpc), allocatable :: saved_band_storage_matrix(:,:)

complex(dpc), allocatable :: eigenvectors(:,:)

complex(dpc), allocatable :: Lbasis_tmp(:,:)

integer :: i, j
integer :: k
integer :: s1, s2
integer :: info


complex(dpc), allocatable :: work(:)
real(dp),     allocatable :: rwork(:)

integer        :: io_unit
character(128) :: filename
logical        :: file_exists

integer        :: debug_unit
character(50)  :: debug_filename

! *************************************************************************



! number of superdiagonals
kd   = nseeds
ldab = kd + 1

ABI_MALLOC(      band_storage_matrix, (ldab,nseeds*kmax))
ABI_MALLOC(saved_band_storage_matrix, (ldab,nseeds*kmax))
!---------------------------------------------------------
! Store banded matrix in banded format
!---------------------------------------------------------
! for UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).

band_storage_matrix(:,:) = cmplx_0

!-----------------------------------------
! Store the alpha and beta matrices
!-----------------------------------------

! loop on all blocks
do k=1,kmax

! alpha blocks
do s2 = 1, nseeds
j = (k-1)*nseeds+s2

do s1 = s2, nseeds
i = j+s1-s2
band_storage_matrix(1+i-j,j) = alpha(s1,s2,k)
end do
end do

! exit when k = kmax, as this beta block does not contribute.
if  (k .eq. kmax) exit

! beta blocks
do s2 = 1, nseeds
j = (k-1)*nseeds+s2

do s1 = 1, s2
i = j+s1-s2+nseeds
band_storage_matrix(1+i-j,j) = beta(s1,s2,k)
end do

end do
end do

saved_band_storage_matrix(:,:) = band_storage_matrix(:,:)

!-----------------------------------------
! Diagonalize the banded matrix
!-----------------------------------------

ABI_MALLOC(eigenvectors, (nseeds*kmax,nseeds*kmax))

ABI_MALLOC(work,(nseeds*kmax))
ABI_MALLOC(rwork,(3*nseeds*kmax-2))

call ZHBEV(                     'V',      & ! compute eigenvalues and eigenvectors
'L',      & ! lower triangular part of matrix is stored in banded_matrix
nseeds*kmax,      & ! dimension of matrix
kd,      & ! number of superdiagonals in banded matrix
band_storage_matrix,      & ! matrix in banded storage
ldab,      & ! leading dimension of banded_matrix
eigenvalues,      & ! eigenvalues of matrix
eigenvectors,      & ! eigenvectors of matrix
nseeds*kmax,      & ! dimension of eigenvector matrix
work, rwork, info )  ! work arrays and info


if ( info /= 0) then
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info,' in ZHBEV (1), gwls_GWlanczos'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if



!----------------------------------------------------------------------------------
! update the Lanczos basis to reflect diagonalization of T matrix
!
!        Note that by definition
!
!                        Q^H . A . Q = T  ==>   A = Q . T . Q^H
!
!        where Q (Lbasis) contains the Lanczos basis.
!
!        Diagonalizing T, ie  T = U . LAMBDA . U^H leads to
!
!                         A  = [ Q.U] . LAMBDA . [Q.U]^H
!
!     The updated basis is thus Q.U == Lbasis . eigenvectors
!----------------------------------------------------------------------------------

! NEVER use matmul!!! It sends temporary arrays to the stack, which can be much smaller
! than needed; this leads to mysterious segfaults!

! Lbasis = matmul(Lbasis,eigenvectors)

! use temporary array, which is PROPERLY ALLOCATED, to perform matrix multiplication
ABI_MALLOC(Lbasis_tmp, (Hsize,nseeds*kmax))

! Compute C = A * B, where A = Lbasis, B = eigenvectors, and C = Lbasis_tmp
call ZGEMM(     'N',     & ! leave array A as is
'N',     & ! leave array B as is
Hsize,     & ! number of rows of A
nseeds*kmax,     & ! number of columns of B
nseeds*kmax,     & ! number of columns of A /rows of B
cmplx_1,     & ! constant alpha
Lbasis,     & ! matrix A
Hsize,     & ! LDA
eigenvectors,     & ! matrix B
nseeds*kmax,     & ! LDB
cmplx_0,     & ! constant beta
Lbasis_tmp,     & ! matrix C
Hsize)       ! LDC

! overwrite initial array
Lbasis(:,:) = Lbasis_tmp(:,:)


ABI_FREE(Lbasis_tmp)

if ( debug .and. mpi_enreg%me == 0)  then
  !----------------------------------------------------------------------------------
  ! For the purpose of debugging, print relevant results to file to check
  ! data consistency. This may not be necessary once the code has been shown to
  ! work properly.
  !----------------------------------------------------------------------------------

  io_unit  = get_unit()
  i = 0
  file_exists = .true.
  do while (file_exists)
  i = i+1
  write(filename,'(A,I0.4,A)') "diagonalize_banded_matrix_",i,".log"
  inquire(file=filename,exist=file_exists)
  end do


  open(io_unit,file=filename,status=files_status_new)
  write(io_unit,10) "#======================================================================================="
  write(io_unit,10) "#                                                                                       "
  write(io_unit,10) "#   This file contains information pertaining to the diagonalization of a banded        "
  write(io_unit,10) "#   matrix, expressed in terms of the Lanczos alpha and beta block arrays.              "
  write(io_unit,10) "#                                                                                       "
  write(io_unit,10) "#======================================================================================="
  write(io_unit,10) "#                                                                                       "
  write(io_unit,12) "#   diagonalization info : ",info,"                                                     "
  write(io_unit,10) "#                                                                                       "
  write(io_unit,10) "#   l                 lambda_l                                                          "
  write(io_unit,10) "#======================================================================================="

  do i = 1, nseeds*kmax
  write(io_unit,13) i, eigenvalues(i)
  end do

  write(io_unit,10) "                                                                                        "
  write(io_unit,10) "#                                                                                       "
  write(io_unit,10) "#    alpha and beta blocks                                                              "
  write(io_unit,10) "#                                                                                       "
  write(io_unit,10) "#======================================================================================="

  ! loop on all blocks
  do k=1,kmax
  write(io_unit,10) "#                        "
  write(io_unit,12) "#   block k = ",k,"      "
  write(io_unit,10) "#                        "
  write(io_unit,15) "#                  alpha:   ||alpha^H-alpha|| = ",  &
  sqrt(sum(abs(alpha(:,:,k)-transpose(conjg(alpha(:,:,k))))**2))
  do s1 = 1, nseeds
  write(io_unit,14) alpha(s1,:,k)
  end do

  write(io_unit,10) "#                        "
  write(io_unit,10) "#                  beta  "
  do s1 = 1, nseeds
  write(io_unit,14) beta(s1,:,k)
  end do


  end do




  write(io_unit,10) "#                                                                                       "
  write(io_unit,10) "#   band storage matrix:                                                                "
  write(io_unit,10) "#======================================================================================="

  do i = 1, ldab
  write(io_unit,14) saved_band_storage_matrix(i,:)
  end do

  close(io_unit)
end if


! clean up memory
ABI_FREE(eigenvectors)
ABI_FREE( work)
ABI_FREE(rwork)
ABI_FREE(band_storage_matrix)
ABI_FREE(saved_band_storage_matrix)


10 format(A)
12 format(A,I5,A)
13 format(I5,5X,F24.12)
14 format(4X,1000(F12.8,SP,F12.8,1X,'i',2X))
15 format(A,ES24.8)

end subroutine diagonalize_lanczos_banded
!!***

end module m_gwls_GWlanczos
!!***

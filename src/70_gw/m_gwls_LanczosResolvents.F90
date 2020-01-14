!!****m* ABINIT/m_gwls_LanczosResolvents
!! NAME
!! m_gwls_LanczosResolvents
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


module m_gwls_LanczosResolvents
!----------------------------------------------------------------------------------------------------
! This module contains routines to compute the matrix elements of the resolent, namely
!
!                       < phi_l | 1/ (H-z) | psi_l' >
!
! where H is the Hamiltonian, and z should be such that we are far from its poles.
!
!----------------------------------------------------------------------------------------------------
! local modules

use m_gwls_utility
use m_gwls_wf
use m_gwls_TimingLog
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_GWlanczos

use defs_basis
use m_abicore
use m_xmpi
use m_io_tools,  only : get_unit

implicit none
save
private
!!***


integer :: LR_kmax, LR_nseeds

complex(dpc), allocatable :: precondition_C(:)       ! this operator is diagonal!
complex(dpc), allocatable :: precondition_one_on_C(:)  ! this operator is diagonal!


complex(dpc), allocatable, public :: Hamiltonian_Qk(:,:) ! Lanczos basis of the Hamiltonian

complex(dpc), allocatable :: LR_alpha(:,:,:)
complex(dpc), allocatable :: LR_beta (:,:,:)
complex(dpc), allocatable :: LR_seeds(:,:)
real(dp),     allocatable :: LR_Hamiltonian_eigenvalues(:)

complex(dpc), allocatable, public :: LR_M_matrix(:,:)
!!***

public :: setup_LanczosResolvents
public :: cleanup_LanczosResolvents
public :: compute_resolvent_column_shift_lanczos
public :: compute_resolvent_column_shift_lanczos_right_vectors
public :: build_preconditioned_Hamiltonian_Lanczos_basis
public :: invert_general_matrix
!!***

contains


!!****f* m_hamiltonian/setup_LanczosResolvents
!! NAME
!!  setup_LanczosResolvents
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_Projected_AT,gwls_Projected_BT
!!
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine setup_LanczosResolvents(kmax, prec)
!----------------------------------------------------------------------------------------------------
!
! This subroutine prepares global work arrays in order to use the Lanczos scheme to 
! compute matrix elements of inverse operators.
!
!
!----------------------------------------------------------------------------------------------------
implicit none

!------------------------------
! input/output variables
!------------------------------

integer, intent(in) :: kmax  ! number of Lanczos steps
logical, intent(in) :: prec  ! should there be preconditioning?

! *************************************************************************


LR_kmax   = kmax
LR_nseeds = 1      ! only one seed

! Prepare the algorithm

ABI_ALLOCATE(precondition_C,        (npw_g))
ABI_ALLOCATE(precondition_one_on_C, (npw_g))


if ( prec ) then
  call set_precondition()

  precondition_one_on_C(:) = sqrt(pcon(:))
  !precondition_one_on_C(:) = pcon(:) ! yields very bad results!
  precondition_C(:)        = cmplx_1/precondition_one_on_C(:)

else 
  precondition_C(:)        = cmplx_1
  precondition_one_on_C(:) = cmplx_1
end if

ABI_ALLOCATE(LR_alpha,       (LR_nseeds,LR_nseeds,LR_kmax))
ABI_ALLOCATE(LR_beta ,       (LR_nseeds,LR_nseeds,LR_kmax))
ABI_ALLOCATE(LR_seeds,       (npw_g,LR_nseeds))
ABI_ALLOCATE(Hamiltonian_Qk, (npw_g,LR_kmax))
ABI_ALLOCATE(LR_Hamiltonian_eigenvalues, (LR_kmax))
ABI_ALLOCATE(LR_M_matrix, (LR_kmax,LR_kmax))

end subroutine setup_LanczosResolvents
!!***

!!****f* m_hamiltonian/cleanup_LanczosResolvents
!! NAME
!!  cleanup_LanczosResolvents
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_Projected_AT,gwls_Projected_BT
!!
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine cleanup_LanczosResolvents
!----------------------------------------------------------------------------------------------------
!
! This subroutine cleans up global arrays.
!
!
!----------------------------------------------------------------------------------------------------
implicit none

! *************************************************************************

ABI_DEALLOCATE(precondition_C)
ABI_DEALLOCATE(precondition_one_on_C)
ABI_DEALLOCATE(LR_alpha)
ABI_DEALLOCATE(LR_beta )
ABI_DEALLOCATE(LR_seeds)
ABI_DEALLOCATE(Hamiltonian_Qk)
ABI_DEALLOCATE(LR_Hamiltonian_eigenvalues)
ABI_DEALLOCATE(LR_M_matrix)

end subroutine cleanup_LanczosResolvents
!!***

!!****f* m_hamiltonian/matrix_function_preconditioned_Hamiltonian
!! NAME
!!  matrix_function_preconditioned_Hamiltonian
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
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine matrix_function_preconditioned_Hamiltonian(vector_out,vector_in,Hsize)
!----------------------------------------------------------------------------------------------------
! This subroutine is a simple wrapper around the preconditioned Hamiltonian operator which acts as
!
!                               C^{-1} . H . C^{-1}
!
! where C^{-2} ~ 1 / T, where T is the kinetic energy.
!
!----------------------------------------------------------------------------------------------------
implicit none
integer,      intent(in)  :: Hsize
complex(dpc), intent(out) :: vector_out(Hsize)
complex(dpc), intent(in)  :: vector_in(Hsize)

! local variables
real(dp)  :: psikg(2,npw_g)

complex(dpc)  :: tmp_vector(Hsize)

! *************************************************************************

! convert from one format to the other

tmp_vector(:) = precondition_one_on_C(:)*vector_in(:)

psikg(1,:) = dble (tmp_vector(:))
psikg(2,:) = dimag(tmp_vector(:))

call Hpsik(psikg)

tmp_vector(:) = cmplx_1*psikg(1,:)+cmplx_i*psikg(2,:)

vector_out(:) = precondition_one_on_C(:)*tmp_vector(:)


end subroutine matrix_function_preconditioned_Hamiltonian
!!***

!!****f* m_hamiltonian/build_preconditioned_Hamiltonian_Lanczos_basis
!! NAME
!!  build_preconditioned_Hamiltonian_Lanczos_basis
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_LanczosResolvents
!!
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine build_preconditioned_Hamiltonian_Lanczos_basis(seed_vector)
!----------------------------------------------------------------------------------------------------
! This function Computes the Lanczos basis of the preconditioned Hamiltonian.
!----------------------------------------------------------------------------------------------------
implicit none
complex(dpc), intent(in) :: seed_vector(npw_g)

integer :: l
integer :: ierr
integer :: mpi_communicator

! *************************************************************************



mpi_communicator = mpi_enreg%comm_fft

! Compute the Lanczos basis as well as the eigenvalues of the T matrix
! Seed must also be preconditioned!
LR_seeds(:,1) = precondition_one_on_C(:)*seed_vector(:)

call block_lanczos_algorithm(mpi_communicator, matrix_function_preconditioned_Hamiltonian, LR_kmax, LR_nseeds, npw_g, LR_seeds, &
&                            LR_alpha, LR_beta, Hamiltonian_Qk)


call diagonalize_lanczos_banded(LR_kmax,LR_nseeds,npw_g, LR_alpha,LR_beta,  & 
Hamiltonian_Qk, LR_Hamiltonian_eigenvalues,.false.)



! Update the Lanczos vectors with the preconditioning matrix
do l = 1, LR_kmax
Hamiltonian_Qk(:,l) =  precondition_C(:)*Hamiltonian_Qk(:,l)
end do

! compute the M matrix

! Compute Q^dagger . Q
call ZGEMM('C','N',LR_kmax,LR_kmax,npw_g,cmplx_1,Hamiltonian_Qk,npw_g,Hamiltonian_Qk,npw_g,cmplx_0,LR_M_matrix,LR_kmax)
call xmpi_sum(LR_M_matrix,mpi_communicator,ierr) ! sum on all processors working on FFT!

do l = 1, LR_kmax
LR_M_matrix(l,:) = LR_M_matrix(l,:)*LR_Hamiltonian_eigenvalues(l)
end do

end subroutine build_preconditioned_Hamiltonian_Lanczos_basis
!!***


!!****f* m_hamiltonian/compute_resolvent_column_shift_lanczos
!! NAME
!!  compute_resolvent_column_shift_lanczos
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_Projected_AT
!!
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine compute_resolvent_column_shift_lanczos(nz, list_z, nvec, list_left_vectors, seed_vector, matrix_elements_resolvent)
!----------------------------------------------------------------------------------------------------
! This function Computes one column of the resolvent using the shift lanczos algorithm,
! namely
!
!
!               I_l(z) = < left_vector_l | xi(z) > 
!
! where |xi(z) > is obtained from the shift lanczos scheme. 
!----------------------------------------------------------------------------------------------------
implicit none

integer     , intent(in) :: nz
complex(dpc), intent(in) :: list_z(nz)

integer     , intent(in) :: nvec
complex(dpc), intent(in) :: list_left_vectors(npw_g,nvec)

complex(dpc), intent(in) :: seed_vector(npw_g)

complex(dpc), intent(out) :: matrix_elements_resolvent(nz,nvec)

! local variables
integer        :: iz, l
complex(dpc)   :: z
integer        :: ierr
integer        :: mpi_communicator

complex(dpc)   :: right_vec(LR_kmax)
complex(dpc)   :: left_vecs(LR_kmax,nvec)
complex(dpc)   :: work_vec(LR_kmax)
complex(dpc)   :: shift_lanczos_matrix(LR_kmax,LR_kmax)

complex(dpc)   :: work_array(npw_g)

! *************************************************************************


! processes will communicate along FFT rows
mpi_communicator = mpi_enreg%comm_fft

!----------------------------------------------------------------------------------------------------
! First, build the Lanczos basis for the preconditioned Hamiltonian
!----------------------------------------------------------------------------------------------------
call build_preconditioned_Hamiltonian_Lanczos_basis(seed_vector)


!----------------------------------------------------------------------------------------------------
! Next, generate arrays which do not depend on z, the external shift.
!----------------------------------------------------------------------------------------------------

! Q^dagger . C^{-2} | seed > 
work_array(:) = precondition_one_on_C(:)**2*seed_vector 

call ZGEMV('C',npw_g,LR_kmax,cmplx_1,Hamiltonian_Qk,npw_g,work_array,1,cmplx_0,right_vec,1)
call xmpi_sum(right_vec,mpi_communicator,ierr) ! sum on all processors working on FFT!

!  Q^dagger | left_vectors > 
call ZGEMM('C','N',LR_kmax,nvec,npw_g,cmplx_1,Hamiltonian_Qk,npw_g,list_left_vectors,npw_g,cmplx_0,left_vecs,LR_kmax)
call xmpi_sum(left_vecs,mpi_communicator,ierr) ! sum on all processors working on FFT!

!----------------------------------------------------------------------------------------------------
! Use shift Lanczos to compute all matrix elements!
!----------------------------------------------------------------------------------------------------

do iz = 1, nz

z = list_z(iz)

! Generate the matrix to be inverted
shift_lanczos_matrix(:,:) = LR_M_matrix(:,:)

do l = 1, LR_kmax
shift_lanczos_matrix(l,l) = shift_lanczos_matrix(l,l)-z
end do 

! since z could be complex, the matrix is not necessarily hermitian. Invert using general Lapack scheme
call invert_general_matrix(LR_kmax,shift_lanczos_matrix)
! the matrix now contains the inverse!


! | work_vec > = M^{-1} . | right_vec >
call ZGEMV('N',LR_kmax,LR_kmax,cmplx_1,shift_lanczos_matrix,LR_kmax,right_vec,1,cmplx_0,work_vec,1)


! matrix_elements = < right_vecs | work_vec > 
call ZGEMV('C',LR_kmax,nvec,cmplx_1,left_vecs,LR_kmax,work_vec,1,cmplx_0,matrix_elements_resolvent(iz,:),1)

end do


end subroutine compute_resolvent_column_shift_lanczos
!!***


!!****f* m_hamiltonian/compute_resolvent_column_shift_lanczos_right_vectors
!! NAME
!!  compute_resolvent_column_shift_lanczos_right_vectors
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_Projected_BT
!!
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine compute_resolvent_column_shift_lanczos_right_vectors(seed_vector, right_vec)
!----------------------------------------------------------------------------------------------------
! The goal of this routine is to participate in the computation of a column of the resolvent 
! using the shift lanczos algorithm. This routine should be called FIRST.
!
! This routine:
!
!       1) builds the preconditioned Hamiltonian Lanczos basis, using the seed_vector, and 
!          stores result in the Module variables.
!
!       2) prepares and returns the Lanczos basis-projected right vectors, ready for further processing.
!
!
!----------------------------------------------------------------------------------------------------
implicit none

complex(dpc), intent(in) :: seed_vector(npw_g)
complex(dpc), intent(out):: right_vec(LR_kmax)

! local variables
integer        :: mpi_communicator
integer        :: ierr

complex(dpc)   :: work_array(npw_g)

! *************************************************************************


! processes will communicate along FFT rows
mpi_communicator = mpi_enreg%comm_fft

!----------------------------------------------------------------------------------------------------
! First, build the Lanczos basis for the preconditioned Hamiltonian
!----------------------------------------------------------------------------------------------------
call build_preconditioned_Hamiltonian_Lanczos_basis(seed_vector)

!----------------------------------------------------------------------------------------------------
! Next, generate arrays which do not depend on z, the external shift.
!----------------------------------------------------------------------------------------------------

! Q^dagger . C^{-2} | seed > 
work_array(:) = precondition_one_on_C(:)**2*seed_vector(:)

call ZGEMV('C', npw_g, LR_kmax, cmplx_1, Hamiltonian_Qk, npw_g, work_array, 1, cmplx_0, right_vec, 1)
call xmpi_sum(right_vec,mpi_communicator,ierr) ! sum on all processors working on FFT!


end subroutine compute_resolvent_column_shift_lanczos_right_vectors
!!***



!!****f* m_hamiltonian/invert_general_matrix
!! NAME
!!  invert_general_matrix
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_LanczosResolvents,gwls_Projected_BT
!!
!! CHILDREN
!!      zgetrf,zgetri
!!
!! SOURCE

subroutine invert_general_matrix(n,matrix)
!----------------------------------------------------------------------------------------------------
! Simple wrapper around lapack routines to invert a general matrix
!----------------------------------------------------------------------------------------------------
implicit none

integer,      intent(in)    :: n
complex(dpc), intent(inout) :: matrix(n,n)

integer :: info
integer :: ipiv(n)

integer        :: debug_unit
character(50)  :: debug_filename




complex(dpc) :: work(n)

! *************************************************************************

call ZGETRF( n, n, matrix, n, ipiv, info )
if ( info /= 0) then        
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info,' in ZGETRF(1), gwls_LanczosResolvents'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if


call ZGETRI( n, matrix, n, ipiv, work, n, info )

if ( info /= 0) then        
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info,' in ZGETRI(1), gwls_LanczosResolvents'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if


end subroutine invert_general_matrix
!!***


end module m_gwls_LanczosResolvents
!!***

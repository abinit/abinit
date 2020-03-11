!!****m* ABINIT/m_gwls_Projected_BT
!! NAME
!! m_gwls_Projected_BT
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


module m_gwls_Projected_BT
!----------------------------------------------------------------------------------------------------
! This module contains routines to compute the projections of the Sternheimer equations for the 
! BT operator, which appears in the Numerical integral.
!----------------------------------------------------------------------------------------------------
! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_TimingLog
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_GWlanczos
use m_gwls_LanczosBasis 
use m_gwls_DielectricArray
use m_gwls_LanczosResolvents
! abinit modules
use defs_basis
use m_abicore
use m_xmpi

implicit none
save
private
!!***

! Module arrays, to be set and used in this module 
integer, public  :: BT_lsternheimer 
complex(dp), public, allocatable :: projected_BT_A_matrix(:,:)
complex(dp), public, allocatable :: projected_BT_L_matrix(:,:)
complex(dp), public, allocatable :: projected_BT_BETA_matrix(:,:)

integer, public  :: BT_lsternheimer_Lanczos 
! We use pointers so that arrays may be allocated in a subroutine
complex(dp), public, pointer:: projected_BT_A_matrix_Lanczos(:,:)
complex(dp), public, pointer:: projected_BT_L_matrix_Lanczos(:,:)
complex(dp), public, pointer:: projected_BT_BETA_matrix_Lanczos(:,:)


integer, public  :: BT_lsternheimer_model_Lanczos 
! We use pointers so that arrays may be allocated in a subroutine
complex(dp), public, pointer :: projected_BT_A_matrix_model_Lanczos(:,:)
complex(dp), public, pointer :: projected_BT_L_matrix_model_Lanczos(:,:)
complex(dp), public, pointer :: projected_BT_BETA_matrix_model_Lanczos(:,:)
!!***


! Public methods
public :: compute_projected_BT_shift_Lanczos
public :: compute_projected_BT_shift_Lanczos_DISTRIBUTED
!!***
contains

!!****f* m_hamiltonian/compute_projected_BT_shift_Lanczos
!! NAME
!!  compute_projected_BT_shift_Lanczos
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_ComputeCorrelationEnergy
!!
!! CHILDREN
!!      cleanup_lanczosresolvents
!!      compute_resolvent_column_shift_lanczos_right_vectors
!!      invert_general_matrix,setup_lanczosresolvents,wf_block_distribute
!!      xmpi_allgather,xmpi_sum,zgemm,zgemv
!!
!! SOURCE

subroutine compute_projected_BT_shift_Lanczos(nfreq, list_external_omega, lmax, modified_Lbasis,         &
kmax_numeric, npt_gauss, dielectric_array, array_integrand )
!----------------------------------------------------------------------------------------------------
!
! This function returns the integrand
!
!         I(w', w) =  sum_{l1,l2} DielectricArray_{l1,l2}(w') * B_{l1,l2}(w',w)
!
! where BT is obtained by Shift Lanczos for all frequencies.  It is assumed that modified_Lbasis
! already contains the properly modified basis vectors.
!
!----------------------------------------------------------------------------------------------------
implicit none

integer,      intent(in) :: nfreq
real(dp),     intent(in) :: list_external_omega(nfreq)
integer,      intent(in) :: lmax
integer,      intent(in) :: kmax_numeric
integer,      intent(in) :: npt_gauss
complex(dpc), intent(in) :: modified_Lbasis(npw_k,lmax)
complex(dpc), intent(in) :: dielectric_array(lmax,lmax,npt_gauss+1)
complex(dpc), intent(out):: array_integrand(npt_gauss+1,nfreq)


! local variables

logical   :: prec

integer   :: l, l1, iw_ext, iw_prime, iw, mb, iblk, nbdblock_lanczos
integer   :: ierr

integer      :: mpi_band_rank 

integer      :: k
integer      :: iz, nz
complex(dpc) :: z

complex(dpc), allocatable :: list_z(:)

real(dp)     :: external_omega, omega_prime


complex(dpc), allocatable :: matrix_elements_resolvent(:,:)

real(dp),     allocatable :: psik_wrk(:,:)
real(dp),     allocatable :: psikb_wrk(:,:)
real(dp),     allocatable :: psikg_wrk(:,:)

complex(dpc), allocatable :: seed_vector(:)

complex(dpc), allocatable :: right_vec_FFT(:)

complex(dpc), allocatable :: right_vec_LA(:,:)
complex(dpc), allocatable :: LR_M_matrix_LA(:,:,:)
complex(dpc), allocatable :: Hamiltonian_Qk_LA(:,:,:) 
complex(dpc), allocatable :: left_vecs_LA(:,:) 
complex(dpc), allocatable :: shift_lanczos_matrix(:,:) 

complex(dpc), allocatable :: work_vec(:) 


real(dp),     allocatable :: real_wrk_vec(:),   imag_wrk_vec(:)
real(dp),     allocatable :: real_wrk_mat(:,:), imag_wrk_mat(:,:)

! *************************************************************************

! prepare the complex frequency array
nz = 2*nfreq*npt_gauss
ABI_ALLOCATE(list_z,(nz))

ABI_ALLOCATE(matrix_elements_resolvent, (lmax,nz))

ABI_ALLOCATE(psik_wrk,  (2,npw_k))
ABI_ALLOCATE(psikb_wrk, (2,npw_kb))
ABI_ALLOCATE(psikg_wrk, (2,npw_g))
ABI_ALLOCATE(seed_vector, (npw_g))


ABI_ALLOCATE(real_wrk_vec, (kmax_numeric*blocksize))
ABI_ALLOCATE(imag_wrk_vec, (kmax_numeric*blocksize))
ABI_ALLOCATE(real_wrk_mat, (kmax_numeric,kmax_numeric*blocksize))
ABI_ALLOCATE(imag_wrk_mat, (kmax_numeric,kmax_numeric*blocksize))

ABI_ALLOCATE(shift_lanczos_matrix, (kmax_numeric,kmax_numeric))

ABI_ALLOCATE(work_vec, (kmax_numeric))

ABI_ALLOCATE(right_vec_FFT,(kmax_numeric) )
ABI_ALLOCATE(right_vec_LA,(kmax_numeric, blocksize) )
ABI_ALLOCATE(LR_M_matrix_LA,(kmax_numeric,kmax_numeric, blocksize) )
ABI_ALLOCATE(Hamiltonian_Qk_LA,(npw_k, kmax_numeric, blocksize) )


ABI_ALLOCATE(left_vecs_LA,(kmax_numeric, lmax) )

iw = 0
do iw_ext = 1, nfreq

external_omega = list_external_omega(iw_ext)

do iw_prime = 1, npt_gauss

! Remember! the first element of the list_omega array, which contain the imaginary 
! integration frequencies, is zero (which need not be computed explicitly)

omega_prime = list_omega(iw_prime+1)

iw = iw+1
list_z(iw) = cmplx_1*external_omega-cmplx_i*omega_prime

iw = iw+1
list_z(iw) = cmplx_1*external_omega+cmplx_i*omega_prime

end do

end do 


! initialize the array to zero
array_integrand(:,:) = cmplx_0


! prepare the shift lanczos scheme
prec = .false.  ! let's not precondition for now
call setup_LanczosResolvents(kmax_numeric,prec)

! Number of blocks of lanczos vectors
nbdblock_lanczos = lmax/blocksize
if (modulo(lmax,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1

mpi_band_rank = mpi_enreg%me_band



!-------------------------------------------------------------------
!
! The shift lanczos scheme will be implemented explicitly in the
! loop below instead of being wrapped in routines in the module
! gwls_LanczosResolvents. The task to be performed is subtle 
! because of the different data distributions and the need to 
! project on all Lanczos vectors.
!
!-------------------------------------------------------------------

! loop on all blocks of Lanczos vectors
do iblk = 1, nbdblock_lanczos

! Convert a block of Lanczos vectors to the FFT configuration

! Change the configuration of the data
do mb =1, blocksize
l = (iblk-1)*blocksize+mb
if (l <= lmax) then 
  psik_wrk(1,:)  = dble ( modified_Lbasis(:,l) )
  psik_wrk(2,:)  = dimag( modified_Lbasis(:,l) )
else
  psik_wrk(:,:)  = zero
end if

psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:) 

end do ! mb 

! change configuration of the data, from LA to FFT
call wf_block_distribute(psikb_wrk,  psikg_wrk, 1) ! LA -> FFT

! compute the seed vector for this block
seed_vector(:) = cmplx_1*psikg_wrk(1,:) + cmplx_i*psikg_wrk(2,:)


! Compute the Lanczos basis for the Hamiltonian, in FFT configuration,
! given the seed vector. Project the right vector onto the thus produced Lanczos basis.
call compute_resolvent_column_shift_lanczos_right_vectors(seed_vector, right_vec_FFT)

! Distribute the data from the FFT configuration back to the LA configuration
! We will use some public data from the LanczosResolvent module.
!
! PATTERN: xmpi_allgather(xval,nelem,recvbuf,spaceComm,ier) 
! unfortunately, there are only interfaces for real arguments, not complex.

call xmpi_allgather( dble(right_vec_FFT), kmax_numeric, real_wrk_vec, mpi_enreg%comm_band, ierr) 
call xmpi_allgather(dimag(right_vec_FFT), kmax_numeric, imag_wrk_vec, mpi_enreg%comm_band, ierr) 



call xmpi_allgather( dble(LR_M_matrix), kmax_numeric**2, real_wrk_mat, mpi_enreg%comm_band, ierr) 
call xmpi_allgather(dimag(LR_M_matrix), kmax_numeric**2, imag_wrk_mat, mpi_enreg%comm_band, ierr) 


do mb =1, blocksize
right_vec_LA    (:,mb) = cmplx_1*real_wrk_vec((mb-1)*kmax_numeric+1:mb*kmax_numeric) + &
&                                    cmplx_i*imag_wrk_vec((mb-1)*kmax_numeric+1:mb*kmax_numeric)


LR_M_matrix_LA(:,:,mb) = cmplx_1*real_wrk_mat(:,(mb-1)*kmax_numeric+1:mb*kmax_numeric) + &
&                                    cmplx_i*imag_wrk_mat(:,(mb-1)*kmax_numeric+1:mb*kmax_numeric)
end do ! mb

do k = 1, kmax_numeric

psikg_wrk(1,:) = dble ( Hamiltonian_Qk(:,k) )
psikg_wrk(2,:) = dimag( Hamiltonian_Qk(:,k) )

! change configuration of the data, from FFT to LA
call wf_block_distribute(psikb_wrk,  psikg_wrk, 2) ! FFT -> LA 

do mb=1, blocksize

Hamiltonian_Qk_LA(:, k, mb) =  cmplx_1*psikb_wrk(1,(mb-1)*npw_k+1:mb*npw_k) &
&                                              +cmplx_i*psikb_wrk(2,(mb-1)*npw_k+1:mb*npw_k) 

end do ! mb

end do ! k


! Data is now available in LA configuration, perfect for linear algebra!


do mb = 1, blocksize
! loop on all vectors in this block

l1 = (iblk-1)*blocksize+mb

if (l1 > lmax) cycle

!  Compute Q^dagger | left_vectors > 

! computes C = alpha *  op(A).op(B) + beta * C
call ZGEMM(                         'C', &  ! First array is hermitian conjugated 
'N', &  ! second array is taken as is
kmax_numeric, &  ! number of rows of matrix op(A)
lmax, &  ! number of columns of matrix op(B)
npw_k, &  ! number of columns of op(A) == number of rows of matrix op(B)
cmplx_1, &  ! alpha
Hamiltonian_Qk_LA(:, :, mb), &  ! A matrix
npw_k, &  ! LDA
modified_Lbasis, &  ! B matrix
npw_k, &  ! LDB
cmplx_0, &  ! beta
left_vecs_LA, &  ! C matrix
kmax_numeric)  ! LDC

call xmpi_sum(left_vecs_LA,mpi_enreg%comm_bandfft,ierr) ! sum on all processors working on LA

! Use shift Lanczos to compute all matrix elements! 

! THE FOLLOWING COULD BE DONE IN PARALLEL INSTEAD OF HAVING EVERY PROCESSOR DUMBLY DO THE SAME THING
do iz = 1, nz

z = list_z(iz)

! Generate the matrix to be inverted
shift_lanczos_matrix(:,:) = LR_M_matrix_LA(:,:,mb)

do k = 1, kmax_numeric
shift_lanczos_matrix(k,k) = shift_lanczos_matrix(k,k)-z
end do 

! since z could be complex, the matrix is not necessarily hermitian. Invert using general Lapack scheme
call invert_general_matrix(kmax_numeric,shift_lanczos_matrix)
! the matrix now contains the inverse!


! | work_vec > = M^{-1} . | right_vec >

! compute y = alpha op(A).x + beta y

call ZGEMV(                       'N', &! A matrix is as is
kmax_numeric, &! number of rows of matrix A
kmax_numeric, &! number of columns of matrix A
cmplx_1, &! alpha
shift_lanczos_matrix, &! matrix A
kmax_numeric, &! LDA
right_vec_LA(:,mb), &! array X
1, &! INC X
cmplx_0, &! beta
work_vec, &! Y array
1)  ! INC Y


! matrix_elements = < right_vecs | work_vec > 

call ZGEMV(                       'C', &! A matrix is hermitan conjugate
kmax_numeric, &! number of rows of matrix A
lmax, &! number of columns of matrix A
cmplx_1, &! alpha
left_vecs_LA, &! matrix A
kmax_numeric, &! LDA
work_vec, &! array X
1, &! INC X
cmplx_0, &! beta
matrix_elements_resolvent(:,iz), &! Y array
1)  ! INC Y



end do ! iz

! update integrand
iw = 0
do iw_ext = 1, nfreq
do iw_prime = 1, npt_gauss

iw = iw+1

! this expression will be normalized by 1/2pi at the end
array_integrand(iw_prime+1,iw_ext) = array_integrand(iw_prime+1,iw_ext)  +    &
sum(dielectric_array(:,l1,iw_prime+1)*  &
(matrix_elements_resolvent(:,iw)+matrix_elements_resolvent(:,iw+1)))

iw = iw+1

end do !iw_prime 
end do !iw_ext 



end do !mb
end do ! iblk

! normalize !
array_integrand(:,:)  = array_integrand(:,:)/(2.0_dp*pi)



call cleanup_LanczosResolvents


ABI_DEALLOCATE(psik_wrk)
ABI_DEALLOCATE(psikb_wrk)
ABI_DEALLOCATE(psikg_wrk)

ABI_DEALLOCATE(seed_vector)
ABI_DEALLOCATE(work_vec)


ABI_DEALLOCATE(right_vec_FFT)
ABI_DEALLOCATE(right_vec_LA)

ABI_DEALLOCATE(LR_M_matrix_LA)

ABI_DEALLOCATE(Hamiltonian_Qk_LA)

ABI_DEALLOCATE(real_wrk_vec)
ABI_DEALLOCATE(imag_wrk_vec)
ABI_DEALLOCATE(real_wrk_mat)
ABI_DEALLOCATE(imag_wrk_mat)


ABI_DEALLOCATE(shift_lanczos_matrix)
ABI_DEALLOCATE(left_vecs_LA)


ABI_DEALLOCATE(list_z)
ABI_DEALLOCATE(matrix_elements_resolvent)

end subroutine compute_projected_BT_shift_Lanczos
!!***



!!****f* m_hamiltonian/compute_projected_BT_shift_Lanczos_DISTRIBUTED
!! NAME
!!  compute_projected_BT_shift_Lanczos_DISTRIBUTED
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_ComputeCorrelationEnergy
!!
!! CHILDREN
!!      cleanup_lanczosresolvents
!!      compute_resolvent_column_shift_lanczos_right_vectors
!!      invert_general_matrix,setup_lanczosresolvents,wf_block_distribute
!!      xmpi_allgather,xmpi_sum,zgemm,zgemv
!!
!! SOURCE

subroutine compute_projected_BT_shift_Lanczos_DISTRIBUTED(nfreq, list_external_omega, lmax,blocksize_eps,       &
model_lanczos_vector_belongs_to_this_node, model_lanczos_vector_index,  &
modified_Lbasis, kmax_numeric, npt_gauss, dielectric_array, array_integrand )
!----------------------------------------------------------------------------------------------------
!
! This function returns the integrand
!
!         I(w', w) =  sum_{l1,l2} DielectricArray_{l1,l2}(w') * B_{l1,l2}(w',w)
!
! where BT is obtained by shift Lanczos for all frequencies.  It is assumed that modified_Lbasis
! already contains the properly modified basis vectors.
!
! This routine performs the same function as compute_projected_BT_shift_Lanczos, but 
! takes into account that the MODEL dielectric array is distributed over the processors.
!
!----------------------------------------------------------------------------------------------------
implicit none

integer,      intent(in) :: nfreq
real(dp),     intent(in) :: list_external_omega(nfreq)
integer,      intent(in) :: lmax, blocksize_eps
integer,      intent(in) :: kmax_numeric
integer,      intent(in) :: npt_gauss
complex(dpc), intent(in) :: modified_Lbasis(npw_k,lmax)
complex(dpc), intent(in) :: dielectric_array(lmax, blocksize_eps, npt_gauss+1)

logical, intent(in) :: model_lanczos_vector_belongs_to_this_node(lmax)
integer, intent(in) :: model_lanczos_vector_index(lmax)

complex(dpc), intent(out):: array_integrand(npt_gauss+1,nfreq)



! local variables

logical   :: prec

integer   :: l, l1, lb, iw_ext, iw_prime, iw, mb, iblk, nbdblock_lanczos
integer   :: ierr

integer      :: mpi_band_rank 

integer      :: k
integer      :: iz, nz
complex(dpc) :: z

complex(dpc), allocatable :: list_z(:)

real(dp)     :: external_omega, omega_prime


complex(dpc), allocatable :: matrix_elements_resolvent(:,:)

real(dp),     allocatable :: psik_wrk(:,:)
real(dp),     allocatable :: psikb_wrk(:,:)
real(dp),     allocatable :: psikg_wrk(:,:)

complex(dpc), allocatable :: seed_vector(:)

complex(dpc), allocatable :: right_vec_FFT(:)

complex(dpc), allocatable :: right_vec_LA(:,:)
complex(dpc), allocatable :: LR_M_matrix_LA(:,:,:)
complex(dpc), allocatable :: Hamiltonian_Qk_LA(:,:,:) 
complex(dpc), allocatable :: left_vecs_LA(:,:) 
complex(dpc), allocatable :: shift_lanczos_matrix(:,:) 

complex(dpc), allocatable :: work_vec(:) 


real(dp),     allocatable :: real_wrk_vec(:),   imag_wrk_vec(:)
real(dp),     allocatable :: real_wrk_mat(:,:), imag_wrk_mat(:,:)

! *************************************************************************

! prepare the complex frequency array
nz = 2*nfreq*npt_gauss
ABI_ALLOCATE(list_z,(nz))

ABI_ALLOCATE(matrix_elements_resolvent, (lmax,nz))

ABI_ALLOCATE(psik_wrk,  (2,npw_k))
ABI_ALLOCATE(psikb_wrk, (2,npw_kb))
ABI_ALLOCATE(psikg_wrk, (2,npw_g))
ABI_ALLOCATE(seed_vector, (npw_g))


ABI_ALLOCATE(real_wrk_vec, (kmax_numeric*blocksize))
ABI_ALLOCATE(imag_wrk_vec, (kmax_numeric*blocksize))
ABI_ALLOCATE(real_wrk_mat, (kmax_numeric,kmax_numeric*blocksize))
ABI_ALLOCATE(imag_wrk_mat, (kmax_numeric,kmax_numeric*blocksize))

ABI_ALLOCATE(shift_lanczos_matrix, (kmax_numeric,kmax_numeric))

ABI_ALLOCATE(work_vec, (kmax_numeric))

ABI_ALLOCATE(right_vec_FFT,(kmax_numeric) )
ABI_ALLOCATE(right_vec_LA,(kmax_numeric, blocksize) )
ABI_ALLOCATE(LR_M_matrix_LA,(kmax_numeric,kmax_numeric, blocksize) )
ABI_ALLOCATE(Hamiltonian_Qk_LA,(npw_k, kmax_numeric, blocksize) )


ABI_ALLOCATE(left_vecs_LA,(kmax_numeric, lmax) )

iw = 0
do iw_ext = 1, nfreq

external_omega = list_external_omega(iw_ext)

do iw_prime = 1, npt_gauss

! Remember! the first element of the list_omega array, which contain the imaginary 
! integration frequencies, is zero (which need not be computed explicitly)

omega_prime = list_omega(iw_prime+1)

iw = iw+1
list_z(iw) = cmplx_1*external_omega-cmplx_i*omega_prime

iw = iw+1
list_z(iw) = cmplx_1*external_omega+cmplx_i*omega_prime

end do

end do 


! initialize the array to zero
array_integrand(:,:) = cmplx_0


! prepare the shift lanczos scheme
prec = .false.  ! let's not precondition for now
call setup_LanczosResolvents(kmax_numeric,prec)

! Number of blocks of lanczos vectors
nbdblock_lanczos = lmax/blocksize
if (modulo(lmax,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1

mpi_band_rank = mpi_enreg%me_band



!-------------------------------------------------------------------
!
! The shift lanczos scheme will be implemented explicitly in the
! loop below instead of being wrapped in routines in the module
! gwls_LanczosResolvents. The task to be performed is subtle 
! because of the different data distributions and the need to 
! project on all Lanczos vectors.
!
!-------------------------------------------------------------------

! loop on all blocks of Lanczos vectors
do iblk = 1, nbdblock_lanczos

! Convert a block of Lanczos vectors to the FFT configuration

! Change the configuration of the data
do mb =1, blocksize
l = (iblk-1)*blocksize+mb
if (l <= lmax) then 
  psik_wrk(1,:)  = dble ( modified_Lbasis(:,l) )
  psik_wrk(2,:)  = dimag( modified_Lbasis(:,l) )
else
  psik_wrk(:,:)  = zero
end if

psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:) 

end do ! mb 

! change configuration of the data, from LA to FFT
call wf_block_distribute(psikb_wrk,  psikg_wrk, 1) ! LA -> FFT

! compute the seed vector for this block
seed_vector(:) = cmplx_1*psikg_wrk(1,:) + cmplx_i*psikg_wrk(2,:)


! Compute the Lanczos basis for the Hamiltonian, in FFT configuration,
! given the seed vector. Project the right vector onto the thus produced Lanczos basis.
call compute_resolvent_column_shift_lanczos_right_vectors(seed_vector, right_vec_FFT)

! Distribute the data from the FFT configuration back to the LA configuration
! We will use some public data from the LanczosResolvent module.
!
! PATTERN: xmpi_allgather(xval,nelem,recvbuf,spaceComm,ier) 
! unfortunately, there are only interfaces for real arguments, not complex.

call xmpi_allgather( dble(right_vec_FFT), kmax_numeric, real_wrk_vec, mpi_enreg%comm_band, ierr) 
call xmpi_allgather(dimag(right_vec_FFT), kmax_numeric, imag_wrk_vec, mpi_enreg%comm_band, ierr) 



call xmpi_allgather( dble(LR_M_matrix), kmax_numeric**2, real_wrk_mat, mpi_enreg%comm_band, ierr) 
call xmpi_allgather(dimag(LR_M_matrix), kmax_numeric**2, imag_wrk_mat, mpi_enreg%comm_band, ierr) 


do mb =1, blocksize
right_vec_LA    (:,mb) = cmplx_1*real_wrk_vec((mb-1)*kmax_numeric+1:mb*kmax_numeric) + &
&                                    cmplx_i*imag_wrk_vec((mb-1)*kmax_numeric+1:mb*kmax_numeric)


LR_M_matrix_LA(:,:,mb) = cmplx_1*real_wrk_mat(:,(mb-1)*kmax_numeric+1:mb*kmax_numeric) + &
&                                    cmplx_i*imag_wrk_mat(:,(mb-1)*kmax_numeric+1:mb*kmax_numeric)
end do ! mb

do k = 1, kmax_numeric

psikg_wrk(1,:) = dble ( Hamiltonian_Qk(:,k) )
psikg_wrk(2,:) = dimag( Hamiltonian_Qk(:,k) )

! change configuration of the data, from FFT to LA
call wf_block_distribute(psikb_wrk,  psikg_wrk, 2) ! FFT -> LA 

do mb=1, blocksize

Hamiltonian_Qk_LA(:, k, mb) =  cmplx_1*psikb_wrk(1,(mb-1)*npw_k+1:mb*npw_k) &
&                                             +cmplx_i*psikb_wrk(2,(mb-1)*npw_k+1:mb*npw_k) 

end do ! mb

end do ! k


! Data is now available in LA configuration, perfect for linear algebra!


do mb = 1, blocksize
! loop on all vectors in this block

l1 = (iblk-1)*blocksize+mb

if (l1 > lmax) cycle

!  Compute Q^dagger | left_vectors > 

! computes C = alpha *  op(A).op(B) + beta * C
call ZGEMM(                         'C', &  ! First array is hermitian conjugated 
'N', &  ! second array is taken as is
kmax_numeric, &  ! number of rows of matrix op(A)
lmax, &  ! number of columns of matrix op(B)
npw_k, &  ! number of columns of op(A) == number of rows of matrix op(B)
cmplx_1, &  ! alpha
Hamiltonian_Qk_LA(:, :, mb), &  ! A matrix
npw_k, &  ! LDA
modified_Lbasis, &  ! B matrix
npw_k, &  ! LDB
cmplx_0, &  ! beta
left_vecs_LA, &  ! C matrix
kmax_numeric)  ! LDC

call xmpi_sum(left_vecs_LA,mpi_enreg%comm_bandfft,ierr) ! sum on all processors working on LA

! Use shift Lanczos to compute all matrix elements! 

! THE FOLLOWING COULD BE DONE IN PARALLEL INSTEAD OF HAVING EVERY PROCESSOR DUMBLY DO THE SAME THING
do iz = 1, nz

z = list_z(iz)

! Generate the matrix to be inverted
shift_lanczos_matrix(:,:) = LR_M_matrix_LA(:,:,mb)

do k = 1, kmax_numeric
shift_lanczos_matrix(k,k) = shift_lanczos_matrix(k,k)-z
end do 

! since z could be complex, the matrix is not necessarily hermitian. Invert using general Lapack scheme
call invert_general_matrix(kmax_numeric,shift_lanczos_matrix)
! the matrix now contains the inverse!


! | work_vec > = M^{-1} . | right_vec >

! compute y = alpha op(A).x + beta y

call ZGEMV(                       'N', &! A matrix is as is
kmax_numeric, &! number of rows of matrix A
kmax_numeric, &! number of columns of matrix A
cmplx_1, &! alpha
shift_lanczos_matrix, &! matrix A
kmax_numeric, &! LDA
right_vec_LA(:,mb), &! array X
1, &! INC X
cmplx_0, &! beta
work_vec, &! Y array
1)  ! INC Y


! matrix_elements = < right_vecs | work_vec > 

call ZGEMV(                       'C', &! A matrix is hermitan conjugate
kmax_numeric, &! number of rows of matrix A
lmax, &! number of columns of matrix A
cmplx_1, &! alpha
left_vecs_LA, &! matrix A
kmax_numeric, &! LDA
work_vec, &! array X
1, &! INC X
cmplx_0, &! beta
matrix_elements_resolvent(:,iz), &! Y array
1)  ! INC Y



end do ! iz

if ( model_lanczos_vector_belongs_to_this_node(l1) ) then

  lb = model_lanczos_vector_index(l1)



  ! update integrand
  iw = 0
  do iw_ext = 1, nfreq
  do iw_prime = 1, npt_gauss

  iw = iw+1

  ! this expression will be normalized by 1/2pi at the end
  array_integrand(iw_prime+1,iw_ext) = array_integrand(iw_prime+1,iw_ext)  +    &
  sum(dielectric_array(:,lb,iw_prime+1)*  &
  (matrix_elements_resolvent(:,iw)+matrix_elements_resolvent(:,iw+1)))

  iw = iw+1

  end do !iw_prime 
  end do !iw_ext 

end if

end do !mb
end do ! iblk

! sum on processors !

call xmpi_sum(array_integrand, mpi_enreg%comm_bandfft, ierr) ! sum on all processors for LA configuration

! normalize !
array_integrand(:,:)  = array_integrand(:,:)/(2.0_dp*pi)



call cleanup_LanczosResolvents


ABI_DEALLOCATE(psik_wrk)
ABI_DEALLOCATE(psikb_wrk)
ABI_DEALLOCATE(psikg_wrk)

ABI_DEALLOCATE(seed_vector)
ABI_DEALLOCATE(work_vec)


ABI_DEALLOCATE(right_vec_FFT)
ABI_DEALLOCATE(right_vec_LA)

ABI_DEALLOCATE(LR_M_matrix_LA)

ABI_DEALLOCATE(Hamiltonian_Qk_LA)

ABI_DEALLOCATE(real_wrk_vec)
ABI_DEALLOCATE(imag_wrk_vec)
ABI_DEALLOCATE(real_wrk_mat)
ABI_DEALLOCATE(imag_wrk_mat)


ABI_DEALLOCATE(shift_lanczos_matrix)
ABI_DEALLOCATE(left_vecs_LA)


ABI_DEALLOCATE(list_z)
ABI_DEALLOCATE(matrix_elements_resolvent)

end subroutine compute_projected_BT_shift_Lanczos_DISTRIBUTED
!!***

end module m_gwls_Projected_BT
!!***

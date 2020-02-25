!!****m* ABINIT/m_gwls_DielectricArray
!! NAME
!! m_gwls_DielectricArray
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


module m_gwls_DielectricArray
!----------------------------------------------------------------------------------------------------
! This module generates and stores the arrays
!
!        { eps^{-1}(iw)-eps_model^{-1}(iw) }  in Lanczos basis
!        { eps_model^{-1}(iw) - 1 }           in model Lanczos basis
!
! It makes sense to build these only once, as they do not depend on the external frequency.
!
!----------------------------------------------------------------------------------------------------

! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_valenceWavefunctions
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_polarisability
use m_gwls_model_polarisability
use m_gwls_GenerateEpsilon
use m_gwls_TimingLog
use m_gwls_QR_factorization
use m_gwls_LanczosBasis

! abinit modules
use defs_basis
use m_abicore
use m_xmpi
use m_cgtools

use m_time,                only : timab
use m_io_tools,            only: get_unit
use m_gaussian_quadrature, only: get_frequencies_and_weights_legendre


implicit none
save

private
!!***

! Frequencies and weights for Legendre integration
real(dp), public, allocatable :: list_omega(:)
real(dp), public, allocatable :: list_weights(:)

! Arrays to store the combinations of dielectric operators
complex(dp), public, allocatable :: model_dielectric_Lanczos_basis(:,:,:)
complex(dp), public, allocatable :: projected_dielectric_Lanczos_basis(:,:,:)
complex(dp), public, allocatable :: eps_m1_minus_eps_model_m1(:,:,:)

complex(dp), public, allocatable :: eps_model_m1_minus_one(:,:,:)
complex(dpc),public, allocatable :: eps_model_m1_minus_one_DISTR(:,:,:)


! dimensions of blocks in the model dielectric matrix
integer,public :: nbdblock_epsilon 
integer,public :: blocksize_epsilon 
logical,public, allocatable :: model_lanczos_vector_belongs_to_this_node(:)
integer,public, allocatable :: model_lanczos_vector_index(:)


! Arrays necessary to project the Sternheimer equation within
! the computation of the dielectric matrix.
complex(dp), public, allocatable :: projected_epsilon_M_matrix(:,:,:)
complex(dp), public, allocatable :: projected_epsilon_B_matrix(:,:,:)
complex(dp), public, allocatable :: projected_epsilon_G_matrix(:,:,:)

integer,public, allocatable  :: list_lsolutions_EpsilonProjected(:)
!!***

public :: generate_frequencies_and_weights
public :: cleanup_projected_Sternheimer_epsilon
public :: compute_eps_m1_minus_eps_model_m1
public :: compute_eps_m1_minus_one
public :: compute_eps_model_m1_minus_one
public :: ProjectedSternheimerEpsilon
!!***

contains

!!****f* m_hamiltonian/generate_frequencies_and_weights
!! NAME
!!  generate_frequencies_and_weights
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
!!      cpu_time,extract_svd,g_to_r,gr_to_g,hpsik,pc_k_valence_kernel
!!      setup_pk_model,sqmr,sqrt_vc_k,wf_block_distribute
!!      write_text_block_in_timing_log,write_timing_log,xmpi_sum,zgemm,zgesv
!!
!! SOURCE

subroutine generate_frequencies_and_weights(npt_gauss)
!--------------------------------------------------------------------------------
!
! This subroutine computes the frequencies and weights necessary for Gauss-Legendre
! quadrature, and stores the results in module arrays.
!
!--------------------------------------------------------------------------------
implicit none

integer, intent(in)  :: npt_gauss


real(dp), allocatable ::   list_omega_tmp(:)
real(dp), allocatable :: list_weights_tmp(:)

integer     :: i

! *************************************************************************

ABI_ALLOCATE(list_omega_tmp,   (npt_gauss))
ABI_ALLOCATE(list_weights_tmp, (npt_gauss))

call get_frequencies_and_weights_legendre(npt_gauss,list_omega_tmp,list_weights_tmp)


ABI_ALLOCATE(list_omega,   (npt_gauss+1))
ABI_ALLOCATE(list_weights, (npt_gauss+1))

! make sure the first frequency in zero!
list_omega(1)   = zero
list_weights(1) = zero

do i = 1,npt_gauss

! inverse the order of the frequency points, as they come out
! in reverse order from the generating subroutine
list_omega  (i+1) = list_omega_tmp  (npt_gauss+1-i)
list_weights(i+1) = list_weights_tmp(npt_gauss+1-i)

end do


ABI_DEALLOCATE(list_omega_tmp)
ABI_DEALLOCATE(list_weights_tmp)


end subroutine generate_frequencies_and_weights
!!***

!!****f* m_hamiltonian/compute_eps_m1_minus_eps_model_m1
!! NAME
!!  compute_eps_m1_minus_eps_model_m1
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
!!      cpu_time,extract_svd,g_to_r,gr_to_g,hpsik,pc_k_valence_kernel
!!      setup_pk_model,sqmr,sqrt_vc_k,wf_block_distribute
!!      write_text_block_in_timing_log,write_timing_log,xmpi_sum,zgemm,zgesv
!!
!! SOURCE

subroutine compute_eps_m1_minus_eps_model_m1(lmax, npt_gauss)
!----------------------------------------------------------------------------------------------------
!
! This subroutine computes the array 
!
!                eps^{-1}(iw) - eps_model^{-1}(iw), 
!
! for all relevant frequencies in the Lanczos basis.
!----------------------------------------------------------------------------------------------------
implicit none

integer ,     intent(in)  :: lmax, npt_gauss

character(256) :: timing_string
real(dp)       :: time1, time2
real(dp)       :: time

integer        :: iw, l
complex(dpc),allocatable  :: dummy_matrix(:,:)
complex(dpc),allocatable  :: iden(:,:)

! *************************************************************************


timing_string = "#"
call write_text_block_in_Timing_log(timing_string)
timing_string = "#        Computing eps^{-1}(iw) - eps_model^{-1}(iw) "
call write_text_block_in_Timing_log(timing_string)
timing_string = "#"
call write_text_block_in_Timing_log(timing_string)


call cpu_time(time1)
! Allocate the module array
ABI_ALLOCATE(eps_m1_minus_eps_model_m1, (lmax,lmax,npt_gauss+1))
ABI_ALLOCATE(dummy_matrix, (lmax,lmax))
ABI_ALLOCATE(iden, (lmax,lmax))


iden = cmplx_0

do l = 1, lmax
iden(l,l) = cmplx_1
end do

do iw = 1, npt_gauss + 1


dummy_matrix(:,:) = projected_dielectric_Lanczos_basis(:,:,iw)

call driver_invert_positive_definite_hermitian_matrix(dummy_matrix,lmax)


eps_m1_minus_eps_model_m1(:,:,iw) = dummy_matrix(:,:)

dummy_matrix(:,:) = model_dielectric_Lanczos_basis(:,:,iw)



call driver_invert_positive_definite_hermitian_matrix(dummy_matrix,lmax)

eps_m1_minus_eps_model_m1(:,:,iw) = eps_m1_minus_eps_model_m1(:,:,iw) - dummy_matrix(:,:)


end do




! Deallocate the arrays which are no longer needed
ABI_DEALLOCATE(model_dielectric_Lanczos_basis)
ABI_DEALLOCATE(projected_dielectric_Lanczos_basis)



ABI_DEALLOCATE(dummy_matrix)
ABI_DEALLOCATE(iden)



call cpu_time(time2)
time = time2-time1

timing_string = "#        Total time   :   "
call write_timing_log(timing_string,time)




end subroutine compute_eps_m1_minus_eps_model_m1
!!***

!!****f* m_hamiltonian/compute_eps_m1_minus_one
!! NAME
!!  compute_eps_m1_minus_one
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
!!      cpu_time,extract_svd,g_to_r,gr_to_g,hpsik,pc_k_valence_kernel
!!      setup_pk_model,sqmr,sqrt_vc_k,wf_block_distribute
!!      write_text_block_in_timing_log,write_timing_log,xmpi_sum,zgemm,zgesv
!!
!! SOURCE

subroutine compute_eps_m1_minus_one(lmax, npt_gauss)
!----------------------------------------------------------------------------------------------------
!
! This subroutine computes the array 
!
!                eps^{-1}(iw) - I
!
! for all relevant frequencies in the Lanczos basis.
!----------------------------------------------------------------------------------------------------
implicit none

integer ,     intent(in)  :: lmax, npt_gauss

character(256) :: timing_string
real(dp)       :: time1, time2
real(dp)       :: time

integer        :: iw, l
complex(dpc),allocatable  :: dummy_matrix(:,:)
complex(dpc),allocatable  :: iden(:,:)

! *************************************************************************



timing_string = "#"
call write_text_block_in_Timing_log(timing_string)
timing_string = "#        Computing eps^{-1}(iw) - I "
call write_text_block_in_Timing_log(timing_string)
timing_string = "#"
call write_text_block_in_Timing_log(timing_string)

call cpu_time(time1)
! Allocate the module array

! The array eps_m1_minus_eps_model_m1 will be used to store
! eps^{-1}-1; we can think of eps_model = I in this case.
ABI_ALLOCATE(eps_m1_minus_eps_model_m1, (lmax,lmax,npt_gauss+1))

ABI_ALLOCATE(dummy_matrix, (lmax,lmax))
ABI_ALLOCATE(iden, (lmax,lmax))

iden = cmplx_0

do l = 1, lmax
iden(l,l) = cmplx_1
end do


do iw = 1, npt_gauss + 1

dummy_matrix(:,:) = projected_dielectric_Lanczos_basis(:,:,iw)

call driver_invert_positive_definite_hermitian_matrix(dummy_matrix,lmax)

eps_m1_minus_eps_model_m1(:,:,iw) = dummy_matrix(:,:)-iden(:,:)

end do


! Deallocate the arrays which are no longer needed
ABI_DEALLOCATE(projected_dielectric_Lanczos_basis)

ABI_DEALLOCATE(dummy_matrix)
ABI_DEALLOCATE(iden)

call cpu_time(time2)
time = time2-time1

timing_string = "#        Total time   :   "
call write_timing_log(timing_string,time)

end subroutine compute_eps_m1_minus_one
!!***

!!****f* m_hamiltonian/compute_eps_model_m1_minus_one
!! NAME
!!  compute_eps_model_m1_minus_one
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
!!      cpu_time,extract_svd,g_to_r,gr_to_g,hpsik,pc_k_valence_kernel
!!      setup_pk_model,sqmr,sqrt_vc_k,wf_block_distribute
!!      write_text_block_in_timing_log,write_timing_log,xmpi_sum,zgemm,zgesv
!!
!! SOURCE

subroutine compute_eps_model_m1_minus_one(lmax_model, npt_gauss, second_model_parameter, epsilon_model_eigenvalues_0)
!----------------------------------------------------------------------------------------------------
!
! This subroutine computes the array 
!
!                eps_model^{-1}(iw) - 1
!
! for all relevant frequencies in the model Lanczos basis.
!
! This array can potentially get very large, as the complementary basis gets big to achieve 
! convergence. Correspondingly, it makes sense to DISTRIBUTE this array on all processors.
!
! The algorithm will go as follows:
!
!               eps_m = 1 - V . P . V, V = sqrt{vc}
!
!       I ) compute VPV, storing blocks on different processors:
!
!               VPV = [ ------- --------      ]  = VPV[lc, nB, nW]
!                  |  [| block |  block |     ]
!                 lc  [|   1   |    2   | ... ]
!                  |  [|       |        |     ]
!                  |  [|       |        |     ]
!                  |  [ ------- ---------     ]
!                        bsize
!
!               This construction *does* involve a fair bit of communications, but it takes
!               a lot less RAM!
!
!       II) once VPV is constructed, do, one frequency at a time:
!               - import all blocks to the HEAD processor
!               - compute eps_m = 1- VPV
!               - invert eps_m^{-1}
!               - subtract identity eps_m^{-1} - 1
!               - redistribute, block by block
!        
!               Doing this frequency by frequency will reduce the RAM weight on the head node.
!
!
!----------------------------------------------------------------------------------------------------
implicit none

integer,  intent(in) :: lmax_model, npt_gauss
real(dp), intent(in) :: second_model_parameter
real(dp), intent(in) :: epsilon_model_eigenvalues_0(lmax_model)

integer  :: l, l1, l2
integer  :: iw 
integer  :: v
integer  :: ierr

real(dp),    allocatable :: psikg_valence(:,:)
real(dp),    allocatable :: psir_valence(:,:,:,:)


real(dp),    allocatable :: psik_wrk(:,:)
real(dp),    allocatable :: psikb_wrk(:,:)
real(dp),    allocatable :: psikg_wrk(:,:)
real(dp),    allocatable :: psikg_tmp(:,:)

complex(dpc),allocatable :: local_Lbasis_conjugated(:,:)


complex(dpc),allocatable :: VPV(:,:,:)
real(dp),    allocatable :: re_buffer(:,:), im_buffer(:,:) 

complex(dpc),allocatable :: vpv_row(:)


complex(dpc),allocatable :: epsilon_head(:,:)

real(dp), allocatable :: re_BUFFER_head(:,:)
real(dp), allocatable :: im_BUFFER_head(:,:)



complex(dpc),allocatable :: YL(:)

character(256) :: timing_string
real(dp)       :: time1, time2, time
real(dp)       :: fft_time1, fft_time2, fft_time
real(dp)       :: prod_time1, prod_time2, prod_time

complex(dpc)   :: z


integer   :: iblk_lanczos, nbdblock_lanczos
integer   :: mb
integer   :: lb
integer   :: ik

integer   :: mpi_communicator
integer   :: mpi_nproc
integer   :: mpi_rank
integer   :: mpi_head_rank


integer,allocatable   :: sendcounts(:), displs(:)

integer   :: sendcount, recvcount


logical   :: head


! timing
real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB

! *************************************************************************


GWLS_TIMAB   = 1541
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



timing_string = "#"
call write_text_block_in_Timing_log(timing_string)
timing_string = "#        computing eps_model_m1_minus_one"
call write_text_block_in_Timing_log(timing_string)
timing_string = "#"
call write_text_block_in_Timing_log(timing_string)


! Number of blocks of lanczos vectors (blocksize = npband)
nbdblock_lanczos = lmax_model/blocksize
if (modulo(lmax_model,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1


! communicator 
mpi_communicator = mpi_enreg%comm_bandfft

! total number of processors in the communicator
mpi_nproc        = xmpi_comm_size(mpi_communicator )

! what is the rank of this processor?
mpi_rank         = xmpi_comm_rank(mpi_communicator )

! rank of the "head" processor
mpi_head_rank    = 0


! number of blocks in the model dielectric matrix, which is equal to the number of processors
nbdblock_epsilon = mpi_nproc


! number of vectors in every block
blocksize_epsilon =  lmax_model/mpi_nproc
if (modulo(lmax_model,mpi_nproc) /= 0) blocksize_epsilon = blocksize_epsilon + 1

! attribute blocks to every nodes, and tabulate ownership in logical array
! This is not *the most efficient* implementation possible, but it is convenient
ABI_ALLOCATE( model_lanczos_vector_belongs_to_this_node, (lmax_model))
ABI_ALLOCATE( model_lanczos_vector_index, (lmax_model))


model_lanczos_vector_index = 0
model_lanczos_vector_belongs_to_this_node = .false.

do l =1, lmax_model
if (mpi_rank == (l-1)/blocksize_epsilon) then
  model_lanczos_vector_belongs_to_this_node(l) = .true.
  model_lanczos_vector_index(l) = l-mpi_rank*blocksize_epsilon
end if
end do

!write(100+mpi_rank,*) "model_lanczos_vector_belongs_to_this_node = ",model_lanczos_vector_belongs_to_this_node(:)
!write(100+mpi_rank,*) "model_lanczos_vector_index                = ",model_lanczos_vector_index(:)                
!flush(100+mpi_rank)

! Prepare the array that will contain the matrix elements of the model operator
ABI_ALLOCATE(VPV, (lmax_model,blocksize_epsilon,npt_gauss+1))

VPV(:,:,:) = cmplx_0

ABI_ALLOCATE(vpv_row, (lmax_model))


! various working arrays
GWLS_TIMAB   = 1542
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

ABI_ALLOCATE(psikg_valence,(2,npw_g))
ABI_ALLOCATE(psir_valence ,(2,n4,n5,n6))


ABI_ALLOCATE(psik_wrk,  (2,npw_k))
ABI_ALLOCATE(psikb_wrk, (2,npw_kb))
ABI_ALLOCATE(psikg_wrk, (2,npw_g))
ABI_ALLOCATE(psikg_tmp, (2,npw_g))

ABI_ALLOCATE(local_Lbasis_conjugated,(npw_k,lmax_model))
ABI_ALLOCATE(YL,(npw_k))

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



fft_time  = zero
prod_time = zero


call cpu_time(time1)
! loop on all valence bands

do v = 1, nbandv


! copy pre-calculated valence state in this covenient local array
psikg_valence(:,:) = kernel_wavefunctions_FFT(:,:,v)

! compute fourier transform of valence state, and conjugate
call g_to_r(psir_valence,psikg_valence) 
psir_valence(2,:,:,:) = -psir_valence(2,:,:,:) 

!--------------------------------------------------------------------------
!
! Step 1: build the modified basis Pc . [ (V^{1/2} l) . psi_v^*]. 
!
!--------------------------------------------------------------------------
GWLS_TIMAB   = 1543
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

! loop on all blocks of lanczos vectors
do iblk_lanczos = 1, nbdblock_lanczos
! loop on all states within this block
do mb = 1, blocksize

! Determine the index of the Lanczos vector
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax_model) then
  psik_wrk(1,:) = dble (Lbasis_model_lanczos(:,l))
  psik_wrk(2,:) = dimag(Lbasis_model_lanczos(:,l))
else
  psik_wrk(:,:) = zero
end if

! apply Coulomb potential
call sqrt_vc_k(psik_wrk)

! Store in array of blocks of wavefunctions
psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)

end do ! mb

call cpu_time(fft_time1)

! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_wrk,1) ! LA -> FFT 


!  generate the vector  Pc [ (sqrt_V_c.l) psi_v^*]

! Compute the real space product, and return to k space, in FFT configuration
call gr_to_g(psikg_tmp,psir_valence, psikg_wrk)


call cpu_time(fft_time2)
fft_time = fft_time + fft_time2-fft_time1

! project
call pc_k_valence_kernel(psikg_tmp)

! Return to LA configuration

! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_tmp, 2) ! FFT -> LA

do mb = 1, blocksize

! Determine the index of the Lanczos vector
l = (iblk_lanczos-1)*blocksize + mb


if ( l <= lmax_model) then
  psik_wrk(:,:) = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) 
  local_Lbasis_conjugated(:,l) = cmplx_1*psik_wrk(1,:)+cmplx_i*psik_wrk(2,:)
end if


end do ! mb

end do !iblk_lanczos 
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


!--------------------------------------------------------------------------
! Step 2: Now that we have the modified basis, compute the matrix
!             elements of the model dielectric operator
!
!
!--------------------------------------------------------------------------
GWLS_TIMAB   = 1544
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

do iw = 2, npt_gauss+1

call setup_Pk_model(list_omega(iw),second_model_parameter)

call cpu_time(prod_time1)
do l1 = 1, lmax_model

! Apply core function Y to left-vector; conjugate
do ik = 1, npw_k         
YL(ik) = model_Y_LA(ik)*conjg(local_Lbasis_conjugated(ik,l1))
end do

! Only compute lower diagonal part of matrix; epsilon is hermitian conjugate!
vpv_row = cmplx_0
do l2 = 1, l1
do ik = 1, npw_k         
vpv_row(l2) = vpv_row(l2) + YL(ik)*local_Lbasis_conjugated(ik,l2)
end do
end do

! the code below is twice as long!
!call ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
!call ZGEMV ( 'T', npw_k, lmax_model, cmplx_1, local_Lbasis_conjugated, npw_k, YL, 1, cmplx_0, vpv_row, 1)

!do l2 = 1, l1
!        eps_model_m1_minus_one(l1,l2,iw) = eps_model_m1_minus_one(l1,l2,iw)     &
!                       -complex_vector_product(YL, local_Lbasis_conjugated(:,l2),npw_k)
!end do

! Sum on all processors, making sure all processors have the total vpv_row
call xmpi_sum(vpv_row, mpi_communicator, ierr) ! sum on all processors for LA configuration

! Each processor takes its slice!
do l2 =1, l1
if ( model_lanczos_vector_belongs_to_this_node(l2) ) then
  lb = model_lanczos_vector_index(l2)
  VPV(l1,lb,iw) = VPV(l1,lb,iw) + vpv_row(l2)
end if
end do

end do
call cpu_time(prod_time2)

prod_time = prod_time + prod_time2-prod_time1

end do ! iw
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

end do ! v

call cpu_time(time2)

time = time2-time1
timing_string = "#        computing VPV                          : "
call write_timing_log(timing_string,time)

timing_string = "#                --- of which is FFT transforms : "
call write_timing_log(timing_string,fft_time)

timing_string = "#                --- of which is products       : "
call write_timing_log(timing_string,prod_time)





GWLS_TIMAB   = 1542
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)
ABI_DEALLOCATE(local_Lbasis_conjugated)
ABI_DEALLOCATE(YL)
ABI_DEALLOCATE(psikg_valence)
ABI_DEALLOCATE(psir_valence)
ABI_DEALLOCATE(psik_wrk)
ABI_DEALLOCATE(psikb_wrk)
ABI_DEALLOCATE(psikg_wrk)
ABI_DEALLOCATE(psikg_tmp)
ABI_DEALLOCATE(vpv_row)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



call cpu_time(time1)
GWLS_TIMAB   = 1547
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


!--------------------------------------------------------------------------------
!
!
! Gather dielectric matrix on head node, invert, and re-distribute. This
! Saves a lot of RAM, without needing the full machinery of ScaLAPACK.
!
!--------------------------------------------------------------------------------


ABI_ALLOCATE(eps_model_m1_minus_one_DISTR, (lmax_model,blocksize_epsilon,npt_gauss+1))
ABI_ALLOCATE(re_buffer, (lmax_model,blocksize_epsilon))
ABI_ALLOCATE(im_buffer, (lmax_model,blocksize_epsilon))

eps_model_m1_minus_one_DISTR(:,:,:) = cmplx_0

! Define the head node, which will invert the dielectric matrices
head = mpi_rank == mpi_head_rank

! Amount of data received and sent
sendcount = lmax_model*blocksize_epsilon
recvcount = lmax_model*blocksize_epsilon

ABI_ALLOCATE(sendcounts,(mpi_nproc))
ABI_ALLOCATE(displs    ,(mpi_nproc))
sendcounts(:) = sendcount 

do l =1, mpi_nproc
displs(l) = (l-1)*sendcount
end do

if (head) then
  ! build and invert the dielectric array
  ! careful!  lmax_model not necessarily equal to blocksize_epsilon*nbdblock_epsilon
  ABI_ALLOCATE(re_BUFFER_head, (lmax_model, blocksize_epsilon*nbdblock_epsilon))
  ABI_ALLOCATE(im_BUFFER_head, (lmax_model, blocksize_epsilon*nbdblock_epsilon))
  ABI_ALLOCATE(epsilon_head, (lmax_model, lmax_model))
else
  !This looks superfluous and it is on large number of systems, but sending these 
  !unallocated in xmpi_scatterv caused 'cannot allocate memory' cryptic errors on 
  !several parallel test farm computers (cronos_gcc46_paral, petrus_nag, inca_gcc44_sdebug)
  ABI_ALLOCATE(re_BUFFER_head, (1,1))
  ABI_ALLOCATE(im_BUFFER_head, (1,1))
  ABI_ALLOCATE(epsilon_head, (1,1))
end if

! Do  one frequency at a time, to avoid overflowing the RAM
do iw = 1, npt_gauss+1
! Gather, except for static case
if ( iw /=1 ) then
  ! Gather VPV on head node, for this frequency
  call xmpi_gather(dble(VPV(:,:,iw)),  sendcount , re_BUFFER_head, recvcount, mpi_head_rank, mpi_communicator,ierr)
  call xmpi_gather(dimag(VPV(:,:,iw)), sendcount , im_BUFFER_head, recvcount, mpi_head_rank, mpi_communicator,ierr)
end if

if ( head ) then

  ! fill the dielectric matrix

  epsilon_head(:,:) = cmplx_0
  if (iw ==1) then
    ! STATIC CASE, diagonal matrix
    do l= 1, lmax_model
    epsilon_head(l,l) = cmplx_1/epsilon_model_eigenvalues_0(l)-cmplx_1
    end do

  else
    ! DYNAMIC CASE, compute
    do l1 =1, lmax_model
    do l2 =1, l1
    z = -cmplx_1*re_BUFFER_head(l1,l2)-cmplx_i*im_BUFFER_head(l1,l2)
    epsilon_head(l1,l2) = z
    epsilon_head(l2,l1) = conjg(z)
    end do
    epsilon_head(l1,l1) = epsilon_head(l1,l1) + cmplx_1
    end do

    ! invert the matrix
    call driver_invert_positive_definite_hermitian_matrix(epsilon_head,lmax_model)

    ! subtract identity
    do l =1, lmax_model
    epsilon_head(l,l) = epsilon_head(l,l) - cmplx_1
    end do 
  end if

  ! copy in head buffer
  re_BUFFER_head(:,:) = zero
  im_BUFFER_head(:,:) = zero
  do l1 =1, lmax_model
  do l2 =1, lmax_model
  z = epsilon_head(l1,l2)
  re_BUFFER_head(l1,l2) = dble(z) 
  im_BUFFER_head(l1,l2) = dimag(z) 
  end do 
  end do 

end if 

! Scatter back the data on the head to all processors
call xmpi_scatterv(re_BUFFER_head, sendcounts, displs, re_buffer, recvcount, mpi_head_rank, mpi_communicator, ierr)
call xmpi_scatterv(im_BUFFER_head, sendcounts, displs, im_buffer, recvcount, mpi_head_rank, mpi_communicator, ierr)

eps_model_m1_minus_one_DISTR(:,:,iw) = cmplx_1*re_buffer(:,:) + cmplx_i*im_buffer(:,:)

end do

ABI_DEALLOCATE(re_BUFFER_head)
ABI_DEALLOCATE(im_BUFFER_head)
ABI_DEALLOCATE(epsilon_head)



!================================================================================
!
! For debugging purposes, store distributed dielectric matrix back in the 
! complete local copies, to insure the rest of the code works.
!
!================================================================================

if (.false.) then
  ! Prepare the array that will contain the matrix elements of the model operator
  ! THIS IS ONLY FOR THE REST OF THE CODE TO WORK; WE WILL REMOVE THIS 
  ! TO SAVE RAM LATER
  ABI_ALLOCATE(eps_model_m1_minus_one, (lmax_model,lmax_model,npt_gauss+1))

  ! initialize the array with zeros
  eps_model_m1_minus_one = cmplx_0

  ! Amount of data received and sent
  sendcount = lmax_model*blocksize_epsilon
  recvcount = lmax_model*blocksize_epsilon*nbdblock_epsilon


  ABI_ALLOCATE(re_BUFFER_head, (lmax_model,blocksize_epsilon*nbdblock_epsilon))
  ABI_ALLOCATE(im_BUFFER_head, (lmax_model,blocksize_epsilon*nbdblock_epsilon))

  do iw = 1, npt_gauss+1

  re_buffer(:,:) = dble( eps_model_m1_minus_one_DISTR(:,:,iw))
  im_buffer(:,:) = dimag(eps_model_m1_minus_one_DISTR(:,:,iw))

  call xmpi_allgather(re_buffer,sendcount,re_BUFFER_head,mpi_communicator,ierr)
  call xmpi_allgather(im_buffer,sendcount,im_BUFFER_head,mpi_communicator,ierr)

  do l = 1, lmax_model
  eps_model_m1_minus_one(:,l,iw) =  cmplx_1*re_BUFFER_head(:,l)+ cmplx_i*im_BUFFER_head(:,l)
  end do
  end do

  ABI_DEALLOCATE(re_BUFFER_head)
  ABI_DEALLOCATE(im_BUFFER_head)
end if

call cpu_time(time2)
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

!================================================================================
!
!================================================================================


time = time2-time1
timing_string = "#        inverting / distributing               :  "

call write_timing_log(timing_string,time)



ABI_DEALLOCATE(re_buffer )
ABI_DEALLOCATE(im_buffer )
ABI_DEALLOCATE(sendcounts)
ABI_DEALLOCATE(displs    )
ABI_DEALLOCATE(VPV       )



GWLS_TIMAB   = 1541
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


end subroutine compute_eps_model_m1_minus_one
!!***

!!****f* m_hamiltonian/cleanup_projected_Sternheimer_epsilon
!! NAME
!!  cleanup_projected_Sternheimer_epsilon
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
!!      cpu_time,extract_svd,g_to_r,gr_to_g,hpsik,pc_k_valence_kernel
!!      setup_pk_model,sqmr,sqrt_vc_k,wf_block_distribute
!!      write_text_block_in_timing_log,write_timing_log,xmpi_sum,zgemm,zgesv
!!
!! SOURCE

subroutine cleanup_projected_Sternheimer_epsilon

implicit none

! *************************************************************************

if(allocated(projected_epsilon_M_matrix)) then
  ABI_DEALLOCATE(projected_epsilon_M_matrix)
end if
if(allocated(projected_epsilon_B_matrix)) then
  ABI_DEALLOCATE(projected_epsilon_B_matrix)
end if
if(allocated(projected_epsilon_G_matrix)) then
  ABI_DEALLOCATE(projected_epsilon_G_matrix)
end if
if(allocated(eps_m1_minus_eps_model_m1)) then
  ABI_DEALLOCATE(eps_m1_minus_eps_model_m1)
end if
if(allocated(list_omega)) then
  ABI_DEALLOCATE(list_omega)
end if
if(allocated(list_weights)) then
  ABI_DEALLOCATE(list_weights)
end if
if(allocated(eps_model_m1_minus_one_DISTR)) then
  ABI_DEALLOCATE(eps_model_m1_minus_one_DISTR)
end if
if(allocated(model_lanczos_vector_belongs_to_this_node)) then
  ABI_DEALLOCATE(model_lanczos_vector_belongs_to_this_node)
end if
if(allocated(model_lanczos_vector_index)) then
  ABI_DEALLOCATE(model_lanczos_vector_index)
end if


end subroutine cleanup_projected_Sternheimer_epsilon
!!***


!!****f* m_hamiltonian/ProjectedSternheimerEpsilon
!! NAME
!!  ProjectedSternheimerEpsilon
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
!!      cpu_time,extract_svd,g_to_r,gr_to_g,hpsik,pc_k_valence_kernel
!!      setup_pk_model,sqmr,sqrt_vc_k,wf_block_distribute
!!      write_text_block_in_timing_log,write_timing_log,xmpi_sum,zgemm,zgesv
!!
!! SOURCE

subroutine ProjectedSternheimerEpsilon(lmax, npt_gauss, second_model_parameter, &
list_projection_frequencies,nfrequencies,&
epsilon_eigenvalues_0,debug,use_model)
!----------------------------------------------------------------------------------------------------
! This subroutine combines in a single subprogram the jobs of previous routines
!
!               - setup_projected_Sternheimer_epsilon
!               - compute_projected_Sternheimer_epsilon
!
! The purpose of this combination is to avoid independent loops on nbandv, requiring the 
! arrays
!               projected_epsilon_M_matrix
!               projected_epsilon_B_matrix
!               projected_epsilon_G_matrix
!
! from scaling like N^3, which grows very large with problem size.
! 
! Thus, this routine:
!
!       -  Computes the frequency-dependent dielectric matrix in the Lanczos basis, using
!          the projected Sternheimer equations.
!
!       -  Computes the frequency-dependent MODEL dielectric matrix in the complementary Lanczos basis.
!
!
! This routine will be verbose and write log files; indeed, large jobs crash in here, it will
! be important to know where/why!
!
! The subroutine also computes the matrix elements on epsilon_model(iw) in the Lanczos basis;
! this is done here to avoid preforming direct products with the valence states again later.
!----------------------------------------------------------------------------------------------------
implicit none

real(dp), parameter     :: svd_tolerance = 1.0e-16_dp

integer,     intent(in) :: lmax, npt_gauss
integer,     intent(in) :: nfrequencies
real(dp),    intent(in) :: list_projection_frequencies(nfrequencies)
logical,     intent(in) :: debug
real(dp),    intent(in) :: epsilon_eigenvalues_0(lmax)

logical,optional,intent(in) :: use_model

real(dp),    intent(in) :: second_model_parameter


integer :: l, l1, l2
integer :: i, iw, v
integer :: recy_i
integer :: lsolutions_max, lsolutions, ls
integer :: projection

complex(dpc), allocatable :: sternheimer_A0(:,:)
complex(dpc), allocatable :: sternheimer_A(:,:)
complex(dpc), allocatable :: sternheimer_B(:,:)
complex(dpc), allocatable :: sternheimer_X(:,:)
complex(dpc), allocatable :: sternheimer_G(:,:)


complex(dpc), allocatable :: dummy_tmp_1(:,:)
complex(dpc), allocatable :: dummy_tmp_2(:,:)

integer, allocatable      :: ipiv(:)



complex(dpc),allocatable :: local_Lbasis(:,:)
complex(dpc),allocatable :: local_Lbasis_conjugated(:,:)

complex(dpc),allocatable :: YL(:)

real(dp), allocatable :: psikg_in(:,:), psikg_out(:,:)

real(dp), allocatable :: psik_wrk(:,:), psikg_wrk(:,:), psikb_wrk(:,:)
real(dp), allocatable :: psi_gamma_l1(:,:), psi_gamma_l2(:,:)

real(dp), allocatable :: psikg_valence(:,:)
real(dp), allocatable :: psir_valence(:,:,:,:)

real(dp), allocatable :: psi_rhs(:,:,:)

real(dp), allocatable :: psikg_VL(:,:)


complex(dpc), allocatable :: check_matrix(:,:), check_matrix2(:,:)
complex(dpc), allocatable :: c_sternheimer_solutions(:,:)
complex(dpc), allocatable :: QR_orthonormal_basis(:,:)

complex(dpc), allocatable :: svd_matrix(:,:)
real   (dp ), allocatable :: svd_values(:)


integer                   :: iblk_lanczos, nbdblock_lanczos
integer                   :: iblk_solutions, nbdblock_solutions
integer                   :: mb

character(128) :: filename
logical        :: file_exists
integer        :: io_unit

character(128) :: filename_log
integer        :: io_unit_log


real(dp)      :: omega

character(256) :: timing_string
real(dp)       :: time1, time2
real(dp)       :: time_exact


integer  :: info
integer  :: ierr

real(dp)  ::  z(2)


logical        :: omega_is_imaginary
real(dp)       :: omega0

logical        :: model
logical        :: write_debug


integer        :: mpi_communicator, mpi_rank, mpi_group

! *************************************************************************


!================================================================================
! Prepare MPI information
!================================================================================

! for LA configuration ,The processors communicate over band+FFT
mpi_communicator = mpi_enreg%comm_bandfft

! what is the rank of this processor, within its group?
mpi_rank  = mpi_enreg%me_fft

! Which group does this processor belong to, given the communicator?
mpi_group = mpi_enreg%me_band




!================================================================================
! Setup a log file, to keep track of the algorithm
!================================================================================


write(filename_log,'(A,I4.4,A)') 'ProjectedSternheimerEpsilon_PROC=',mpi_enreg%me,'.log'

io_unit_log = get_unit()
open(io_unit_log,file=filename_log,status=files_status_new)
write(io_unit_log,10) ''
write(io_unit_log,10) '#===================================================================================================='
write(io_unit_log,10) "#                     ProjectedSternheimerEpsilon: log file                                          "
write(io_unit_log,10) "#                     -------------------------------------------------------------------            "
write(io_unit_log,10) "#                                                                                                    "
write(io_unit_log,10) "# This file tracks the algorithm in the routine ProjectedSternheimerEpsilon. The goal is to          " 
write(io_unit_log,10) "# establish where the algorithm crashes if it does, and/or  to track are far along the code is.      "
write(io_unit_log,10) '#'
write(io_unit_log,10) '#  MPI data for this process:'
write(io_unit_log,10) '#'
write(io_unit_log,22) '#    mpi_rank :',mpi_rank,'  (rank of this processor in its band group)'
write(io_unit_log,22) '#    mpi_group:',mpi_group,' (band group to which this processor belongs)'
write(io_unit_log,10) '#===================================================================================================='
flush(io_unit_log)


!================================================================================
! Setup timing; prepare arrays
!================================================================================

write(io_unit_log,10) " - Preparing and allocating arrays ...."
flush(io_unit_log)



timing_string = "#"
call write_text_block_in_Timing_log(timing_string)
timing_string = "#        ProjectedSternheimerEpsilon "
call write_text_block_in_Timing_log(timing_string)
timing_string = "#"
call write_text_block_in_Timing_log(timing_string)

! Allocate the module array
ABI_ALLOCATE(projected_dielectric_Lanczos_basis, (lmax,lmax,npt_gauss+1))

projected_dielectric_Lanczos_basis(:,:,:) = cmplx_0

! initialize zero frequency with exact solution
do l = 1, lmax
projected_dielectric_Lanczos_basis(l,l,1) = cmplx_1*epsilon_eigenvalues_0(l)
end do


! initialize other frequencies with the identity
do iw = 2, npt_gauss + 1
do l = 1, lmax
projected_dielectric_Lanczos_basis(l,l,iw) = cmplx_1
end do
end do

time_exact = zero


!================================================================================
! Parallelisation of the code is subtle; Hamiltonian must act over 
! FFT rows, we must be careful with memory, etc...
!
! We will parallelise in block of Lanczos vectors, not over bands.
!================================================================================

! Number of blocks of lanczos vectors
nbdblock_lanczos = lmax/blocksize
if (modulo(lmax,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1


if (present(use_model)) then
  model = use_model
else
  model = .true.
end if

if (model) then
  ! Prepare the array that will contain the matrix elements of the model operator
  ABI_ALLOCATE(model_dielectric_Lanczos_basis, (lmax,lmax,npt_gauss+1))

  ! initialize with zero. NOT with the identity, in order to avoid extra communications
  ! (see below)
  model_dielectric_Lanczos_basis(:,:,:) = cmplx_0

end if

! various working arrays

ABI_ALLOCATE(psikg_valence    ,(2,npw_g))
ABI_ALLOCATE(psir_valence     ,(2,n4,n5,n6))


ABI_ALLOCATE(psikg_VL ,(2,npw_g))


ABI_ALLOCATE(psik_wrk         ,(2,npw_k))
ABI_ALLOCATE(psikb_wrk        ,(2,npw_kb))
ABI_ALLOCATE(psikg_wrk        ,(2,npw_g))


ABI_ALLOCATE(psi_gamma_l1     ,(2,npw_k))
ABI_ALLOCATE(psi_gamma_l2     ,(2,npw_k))

ABI_ALLOCATE(psi_rhs          ,(2,npw_k,lmax))


ABI_ALLOCATE(psikg_in   ,(2,npw_g))
ABI_ALLOCATE(psikg_out  ,(2,npw_g))


! maximal possible dimension of the solution space
! +1 because the solutions at $\omega=\infty$ are free.
! +1 if recycling is activated, because the solutions at $\omega=0$ are then available.
i=1
if(dtset%gwls_recycle == 1 .or. dtset%gwls_recycle == 2) then
  i=2
end if
lsolutions_max = lmax*(nfrequencies+i)


ABI_ALLOCATE(local_Lbasis,           (npw_k,lmax))
ABI_ALLOCATE(local_Lbasis_conjugated,(npw_k,lmax))
ABI_ALLOCATE(YL,(npw_k))

ABI_ALLOCATE(c_sternheimer_solutions,(npw_k,lsolutions_max))
ABI_ALLOCATE(QR_orthonormal_basis   ,(npw_k,lsolutions_max))


omega_is_imaginary = .true.


ABI_ALLOCATE(svd_matrix,(npw_k,lsolutions_max))
ABI_ALLOCATE(svd_values,(lsolutions_max))


! Prepare files for writing
write_debug = debug .and. mpi_enreg%me == 0

if ( write_debug ) then

  write(filename,'(A)') "ProjectedSternheimerEpsilon.log"
  inquire(file=filename,exist=file_exists)

  i = 0
  do while (file_exists)
  i = i+1
  write (filename,'(A,I0.4,A)') "ProjectedSternheimerEpsilon_",i,".log"
  inquire(file=filename,exist=file_exists)
  end do


  io_unit = get_unit()
  open(io_unit,file=filename,status=files_status_new)
  write(io_unit,10) ''
  write(io_unit,10) '#===================================================================================================='
  write(io_unit,10) "#                     Building the dielectic matrix using projected Sternheimer equation             "
  write(io_unit,10) "#                     -------------------------------------------------------------------            "
  write(io_unit,10) "#                                                                                                    "
  write(io_unit,10) '# This file contains some tests to check if the various elements entering the projected  '
  write(io_unit,10) '# dielectric matrix have the right properties. At this point, this is mostly for debugging.'
  write(io_unit,10) '# The wavefunctions and other related arrays are stored in reciprocal space.'
  write(io_unit,10) '#'
  write(io_unit,10) '#===================================================================================================='
  write(io_unit,10) ''
  flush(io_unit)
end if

if (debug) then
  ABI_ALLOCATE(check_matrix ,(lsolutions_max,lsolutions_max))
  ABI_ALLOCATE(check_matrix2,(lsolutions_max,lsolutions_max))
end if

!================================================================================
!  Loop on all valence bands
!================================================================================




do v = 1, nbandv
write(io_unit_log,10) '#===================================================================================================='
write(io_unit_log,20) '#  valence band index:', v
write(io_unit_log,10) '#===================================================================================================='
flush(io_unit_log)




if ( write_debug ) then
  write(io_unit,10) '#===================================================================================================='
  write(io_unit,20) '#  valence band index:', v
  write(io_unit,10) '#===================================================================================================='
  flush(io_unit)
end if

write(io_unit_log,10) '   - Fourier transform valence state ...'
flush(io_unit_log)


! copy pre-calculated valence state in this covenient local array
psikg_valence(:,:) = kernel_wavefunctions_FFT(:,:,v)

! compute fourier transform of valence state, and conjugate
call g_to_r(psir_valence,psikg_valence) 
psir_valence(2,:,:,:) = -psir_valence(2,:,:,:) 

! loop on all blocks of lanczos vectors
write(io_unit_log,10) '   - Loop on all lanczos blocks to generate modified basis and Sternheimer RHS:'
flush(io_unit_log)
do iblk_lanczos = 1, nbdblock_lanczos
!--------------------------------------------------------------------------
! Below, we build the modified basis, [ (V^{1/2} l)^* . psi_v], 
! as well as conjugated, projected form, Pc . [ (V^{1/2} l) . psi_v^*]. 
!
! It is very irritating to have to do it this way, but I don't
! see an alternative; see discussion below.
!
!--------------------------------------------------------------------------

write(io_unit_log,23) '       iblk_lanczos  = ',iblk_lanczos," / ",nbdblock_lanczos


write(io_unit_log,10) '          -- Prepare modified basis computation...'
flush(io_unit_log)



! loop on all states within this block
do mb = 1, blocksize

! Determine the index of the Lanczos vector
l = (iblk_lanczos-1)*blocksize + mb

! take a single lanczos vector

if ( l <= lmax) then
  psik_wrk(1,:) = dble (Lbasis_lanczos(:,l))
  psik_wrk(2,:) = dimag(Lbasis_lanczos(:,l))
else
  psik_wrk(:,:) =  zero
end if


! Apply coulomb potential
call sqrt_vc_k(psik_wrk)

! Store in array of blocks of wavefunctions
psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)

end do ! mb

! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_VL,1) ! LA -> FFT 

! psikg_VL now contains | V^1/2 . l >, in FFT configuration


!----------------------------------------------------------
! i) Compute the modified basis
!----------------------------------------------------------
write(io_unit_log,10) '          -- compute modified basis ...'
flush(io_unit_log)



! Fourier transform to real space, and conjugate (psir1 is a global work array)
call g_to_r(psir1,psikg_VL) 
psir1(2,:,:,:) = -psir1(2,:,:,:) ! IS THIS STACK-DANGEROUS?

! Compute the real space product, and return to k space, in FFT configuration
call gr_to_g(psikg_wrk,psir1,psikg_valence)

! psikg_wrk contains | (V^1/2 . l)^* phi_v >, in FFT configuration

! return to LA representation
call wf_block_distribute(psikb_wrk,  psikg_wrk,2) ! FFT -> LA 

! store data, in LA representation
do mb = 1, blocksize
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax) then
  ! local_Lbasis
  psik_wrk(:,:)     = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)
  local_Lbasis(:,l) = cmplx_1*psik_wrk(1,:)+cmplx_i*psik_wrk(2,:)
end if

end do !mb 

!--------------------------------------------------------------------------
! ii) Build the Sternheimer coefficients, at various frequencies,
!     to define the Sternheimer basis.
!--------------------------------------------------------------------------
write(io_unit_log,20) '          -- Compute Sternheimer RHS...'
flush(io_unit_log)


! psikg_wrk still contains | (V^1/2 . l)^* phi_v >, in FFT configuration

!  Create right-hand-side of Sternheimer equation, in FFT configuration
call pc_k_valence_kernel(psikg_wrk)
call Hpsik(psikg_in,psikg_wrk,eig(v))
call pc_k_valence_kernel(psikg_in)
psikg_in(:,:) = -psikg_in(:,:) ! IS THIS STACK-DANGEROUS?

! return RHS  to LA representation, for explicit storage 
call wf_block_distribute(psikb_wrk,  psikg_in,2) ! FFT -> LA 

! store data, in LA representation
do mb = 1, blocksize
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax) then
  ! psi_rhs
  psi_rhs(:,:,l)    = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)
end if

end do !mb 

!----------------------------------------------------------
! iii) extract solutions for all projection frequencies
!----------------------------------------------------------
write(io_unit_log,20) '          -- Extract solutions for all projection frequencies...'
flush(io_unit_log)



do iw = 1, nfrequencies

omega0 = list_projection_frequencies(iw)
! Solve Sternheimer equation

projection = 0
if(omega0 < 1d-12) projection=1

! solve A x = b, over the whole lanczos block
call sqmr(psikg_in, psikg_out, eig(v), projection, omega0, omega_is_imaginary)

! return LA representation, for explicit storage 
call wf_block_distribute(psikb_wrk,  psikg_out, 2) ! FFT -> LA 

do mb = 1, blocksize
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax) then
  ls = (l-1)*nfrequencies+iw

  psik_wrk(:,:)     = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)

  c_sternheimer_solutions(:,ls)= cmplx_1*psik_wrk(1,:)+cmplx_i*psik_wrk(2,:)
end if

end do ! mb

end do ! iw

!----------------------------------------------------------
! iv) Compute the conjugated, projected modified basis
!----------------------------------------------------------

write(io_unit_log,20) '          -- compute the conjugated, projected modified basis...'
flush(io_unit_log)




! Compute the real space product, | (V^1/2. l) phi_v^* > and return to k space, in FFT configuration
call gr_to_g(psikg_wrk, psir_valence, psikg_VL)

! project on conduction states
call pc_k_valence_kernel(psikg_wrk)

! return to LA representation
call wf_block_distribute(psikb_wrk,  psikg_wrk,2) ! FFT -> LA 

! store back, in LA configuration
do mb = 1, blocksize
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax) then
  psik_wrk(:,:) = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)
  local_Lbasis_conjugated(:,l) = cmplx_1*psik_wrk(1,:)+cmplx_i*psik_wrk(2,:)
end if

end do !mb 

end do ! iblk_lanczos



write(io_unit_log,10) '   - Read in solutions at w=0 and/or w = infinity, if appropriate...'
flush(io_unit_log)

! Begin with the storage of the solutions at $\omega = 0.$, which are free.
if(dtset%gwls_recycle == 1) then 
  c_sternheimer_solutions(:,lsolutions_max-2*lmax+1:lsolutions_max-lmax) = cmplx_1*Sternheimer_solutions_zero(1,:,:,v) + &
  &                                                                             cmplx_i*Sternheimer_solutions_zero(2,:,:,v)
end if
if(dtset%gwls_recycle == 2) then
  do i=1,lmax
  recy_i = (i-1)*nbandv + v       
  !BUG : On petrus, NAG 5.3.1 + OpenMPI 1.6.2 cause read(...,rec=i) to read the data written by write(...,rec=i+1).
  read(recy_unit,rec=recy_i) psik_wrk
  c_sternheimer_solutions(:,lsolutions_max-2*lmax+i) = cmplx_1*psik_wrk(1,:) + cmplx_i*psik_wrk(2,:)
  end do
end if

! and then continue with the storage of the vectors on which the Sternheimer solutions will be projected.
c_sternheimer_solutions(:,lsolutions_max-lmax+1:lsolutions_max) = cmplx_1*psi_rhs(1,:,:) + cmplx_i*psi_rhs(2,:,:)
! Previously was = local_Lbasis; but analysis in the Lanczos article reveals psi_rhs should be better.
! Furthermore, tests show that, with psi_rhs, silane@1Ha has Sigma_c 0.01mHa away from the result with 
! gwls_list_proj_freq 0.0 1.0, in contrast with local_Lbasis, which has Sigma_c 0.3mHa away from the same result.

if ( model ) then

  write(io_unit_log,10) '   - USE MODEL: model = .true., hence compute model model dielectric matrix...'
  flush(io_unit_log)



  !--------------------------------------------------------------------------
  ! Now that we have the modified basis, compute the matrix
  ! elements of the model dielectric operator
  !
  ! CAREFUL!
  !
  !        The model is given by 
  !
  !                P_model(iw) = sum_{v} phi_v(r) P_c.Y(iw).P_c phi_v^*(r')
  !
  !    such that        
  !
  !    <l1 | eps_model(iw) | l2 >  = delta_{l1,l2} 
  !    - sum_{v} < (V^{1/2}.l1).phi_v^*| Pc . Y . Pc | (V^{1/2}.l2).phi_v^* >
  !
  ! But local_Lbasis defined above corresponds to
  !                                Pc | (V^{1/2} .l )^* phi_v >.
  !
  ! This is why we must define local_Lbasis_conjugated, of the form
  !                                Pc | (V^{1/2} .l ) phi_v^* >.
  !
  !--------------------------------------------------------------------------

  do iw = 1, npt_gauss+1

  call setup_Pk_model(list_omega(iw),second_model_parameter)

  ! Only build the lower triangular part; the upper triangular part is obtained from the Hermitian conjugate
  do l1 = 1, lmax

  YL(:) = model_Y_LA(:)*local_Lbasis_conjugated(:,l1)
  do l2 = 1, l1
  model_dielectric_Lanczos_basis(l1,l2,iw) = model_dielectric_Lanczos_basis(l1,l2,iw)  &
  -complex_vector_product(YL, local_Lbasis_conjugated(:,l2),npw_k)

  end do
  end do

  end do ! iw


end if


!--------------------------------------------------------------------------
! Check explicitly that solutions satisfy the Sternheimer equations
!--------------------------------------------------------------------------

if ( debug ) then

  if (write_debug) then
    write(io_unit,10) "#--------------------------------------------------------------------------------"
    write(io_unit,10) "# Check explicitly that solutions satisfy the Sternheimer equation.              "
    write(io_unit,10) "#                                                                                "
    write(io_unit,10) "# Define:                                                                        "
    write(io_unit,10) "# E_l = || (omega^2+[H-Ev]^2) |phi_l> + Pc.[H-Ev].Pc |(V^1/2.q_l)^*.phi_v >  ||  "
    write(io_unit,10) "#--------------------------------------------------------------------------------"
    write(io_unit,10) '#  l                Im[omega] (Ha)                  E_l'
    write(io_unit,10) "#--------------------------------------------------------------------------------"
    flush(io_unit)
  end if


  do iblk_lanczos = 1, nbdblock_lanczos

  ! loop on all states within this block
  do mb = 1, blocksize

  l = (iblk_lanczos-1)*blocksize + mb


  if ( l <= lmax) then
    ! psik_wrk = | (V^{1/2} .l )^* phi_v >
    psik_wrk(1,:) = dble (local_Lbasis(:,l))
    psik_wrk(2,:) = dimag(local_Lbasis(:,l))
  else
    psik_wrk(:,:) = zero
  end if

  ! Store in array of blocks of wavefunctions
  psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)
  end do ! mb

  ! Transform to FFT representation
  call wf_block_distribute(psikb_wrk,  psikg_wrk, 1) ! LA -> FFT 

  !  Create right-hand-side of Sternheimer equation
  call pc_k_valence_kernel(psikg_wrk)
  call Hpsik(psikg_in,psikg_wrk,eig(v))
  call pc_k_valence_kernel(psikg_in)
  psikg_in(:,:)  = -psikg_in(:,:) ! IS THIS STACK-DANGEROUS?

  do iw = 1, nfrequencies
  ! loop on all states within this block
  do mb = 1, blocksize

  l = (iblk_lanczos-1)*blocksize + mb

  ls = (l-1)*nfrequencies+iw

  if ( l <= lmax) then
    psik_wrk(1,:) = dble (c_sternheimer_solutions(:,ls))
    psik_wrk(2,:) = dimag(c_sternheimer_solutions(:,ls))
  else
    psik_wrk(:,:) = zero
  end if

  ! Store in array of blocks of wavefunctions
  psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)
  end do

  ! Transform to FFT representation
  call wf_block_distribute(psikb_wrk,  psikg_wrk, 1) ! LA -> FFT 


  omega0 = list_projection_frequencies(iw)

  psikg_out(:,:) = omega0**2*psikg_wrk(:,:)

  call Hpsik(psikg_wrk,cte=eig(v))
  call Hpsik(psikg_wrk,cte=eig(v))


  psikg_out(:,:) = psikg_out(:,:) + psikg_wrk(:,:)-psikg_in(:,:)
  ! psikg_out now contains [ w0^2 + [H-epsilon_v]^2 ] | x > - |RHS>, in FFT configuration.

  ! bring it back to LA configuration

  ! Transform to FFT representation
  call wf_block_distribute(psikb_wrk,  psikg_out, 2) ! FFT -> LA 

  do mb = 1, blocksize
  l = (iblk_lanczos-1)*blocksize + mb

  if ( l <= lmax) then
    psik_wrk(:,:) = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)

    z(:) = cg_zdotc(npw_k ,psik_wrk, psik_wrk)

    call xmpi_sum(z, mpi_communicator,ierr) ! sum on all processors working on FFT!
    if (write_debug) write(io_unit,21)  l, omega0, sqrt(z(1))
  end if
  end do ! mb

  end do ! iw
  end do ! iblk_lanczos 

end if

!--------------------------------------------------------------------------
! Step 5: Perform a singular value decomposition to extract a
!         linearly independent basis for the solution space.
!
!--------------------------------------------------------------------------
write(io_unit_log,10) '   - Perform SVD to extract linearly independent basis to Sternheimer equation...'
flush(io_unit_log)




svd_matrix(:,:) =  c_sternheimer_solutions(:,:)

call extract_SVD(mpi_communicator, npw_k,lsolutions_max,svd_matrix,svd_values)

if ( write_debug ) then
  write(io_unit,10) "#--------------------------------------------------------------------------------"
  write(io_unit,10) "# Check the singular value decomposition of the arrays"
  write(io_unit,10) '#  l                          svd'
  write(io_unit,10) "#--------------------------------------------------------------------------------"
  flush(io_unit)
end if

lsolutions = 0
QR_orthonormal_basis(:,:) = cmplx_0
do l=1, lsolutions_max
if (svd_values(l) > svd_tolerance ) then
  lsolutions = lsolutions + 1

  if ( write_debug ) then
    write(io_unit,14)   l,svd_values(l)
    flush(io_unit)
  end if
  QR_orthonormal_basis(:,l) = svd_matrix(:,l)

else
  if ( write_debug ) then
    write(io_unit,15)   l,svd_values(l),' SVD value too small! Vector to be discarded!'
    flush(io_unit)
  end if
end if
end do

!--------------------------------------------------------------------------
! Step 6: project all relevant arrays onto the newly defined orthonormal
!         basis.
!--------------------------------------------------------------------------
write(io_unit_log,10) '   - Compute the B matrix...'
flush(io_unit_log)

if (debug) then
  check_matrix(:,:) = cmplx_0
  do l = 1, lsolutions
  check_matrix(l,l) = -cmplx_1
  end do
end if

ABI_ALLOCATE(ipiv                ,(lsolutions))
ABI_ALLOCATE(sternheimer_A0      ,(lsolutions,lsolutions))
ABI_ALLOCATE(sternheimer_B       ,(lsolutions,lmax))
ABI_ALLOCATE(sternheimer_G       ,(lmax,lsolutions))

! Compute the X matrix and the check_matrix
do l1 = 1, lsolutions

psi_gamma_l1(1,:) = real (QR_orthonormal_basis(:,l1))
psi_gamma_l1(2,:) = dimag(QR_orthonormal_basis(:,l1))

do l2 = 1, lsolutions

if (l2 <= lmax) then
  z(:) = cg_zdotc(npw_k, psi_gamma_l1,  psi_rhs(:,:,l2))
  call xmpi_sum(z, mpi_communicator,ierr) ! sum on all processors for LA configuration

  sternheimer_B(l1,l2) = cmplx_1*z(1)+cmplx_i*z(2)
end if

if (debug)  then
  psi_gamma_l2(1,:) = dble (QR_orthonormal_basis(:,l2))
  psi_gamma_l2(2,:) = dimag(QR_orthonormal_basis(:,l2))

  z(:) = cg_zdotc(npw_k, psi_gamma_l1,  psi_gamma_l2)  
  call xmpi_sum(z, mpi_communicator,ierr) ! sum on all processors 

  check_matrix(l1,l2) = check_matrix(l1,l2) + cmplx_1*z(1)+cmplx_i*z(2)
end if

end do ! l2
end do ! l1


! Number of blocks of solution vectors
nbdblock_solutions = lsolutions/blocksize

if (modulo(lsolutions,blocksize) /= 0) nbdblock_solutions = nbdblock_solutions + 1


write(io_unit_log,10) '   - Compute the A0 matrix...'
flush(io_unit_log)

! Compute the A matrix
do iblk_solutions =1, nbdblock_solutions 

do mb = 1, blocksize
l2 = (iblk_solutions-1)*blocksize + mb

if ( l2 <= lsolutions) then
  psik_wrk(1,:) = dble (QR_orthonormal_basis(:,l2))
  psik_wrk(2,:) = dimag(QR_orthonormal_basis(:,l2))
else
  psik_wrk(:,:) =  zero
end if

! Store in array of blocks of wavefunctions
psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)
end do

! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_wrk,1) ! LA -> FFT 

! act twice with the Hamiltonian operator
call Hpsik(psikg_out,psikg_wrk,eig(v))
call Hpsik(psikg_out,cte=eig(v))

! return to LA representation
call wf_block_distribute(psikb_wrk,  psikg_out,2) ! FFT -> LA

do mb = 1, blocksize
l2 = (iblk_solutions-1)*blocksize + mb

if ( l2 <= lsolutions) then
  psik_wrk(:,:)     = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)
else
  psik_wrk(:,:)     = zero
end if

do l1 = 1, lsolutions

psi_gamma_l1(1,:) = real (QR_orthonormal_basis(:,l1))
psi_gamma_l1(2,:) = dimag(QR_orthonormal_basis(:,l1))

z(:) = cg_zdotc(npw_k, psi_gamma_l1,  psik_wrk)
call xmpi_sum(z,mpi_communicator,ierr) ! sum on all processors working on FFT!


if ( l2 <= lsolutions ) then
  sternheimer_A0(l1,l2) = cmplx_1*z(1)+cmplx_i*z(2)
end if


end do ! l1
end do ! mb
end do ! iblk_solutions 


if (debug) then
  ! HERE we use a dummy variable to avoid operations which might blow the stack!
  ! Stack overflows lead to hard-to-find bugs; let's avoid putting them in here.
  ABI_ALLOCATE(dummy_tmp_1,(lsolutions,lsolutions))
  ABI_ALLOCATE(dummy_tmp_2,(lsolutions,lsolutions))


  dummy_tmp_1(:,:)= transpose(sternheimer_A0(:,:))
  dummy_tmp_2(:,:)= conjg(dummy_tmp_1(:,:))

  check_matrix2(:,:) = sternheimer_A0(:,:)-dummy_tmp_2(:,:)

  ABI_DEALLOCATE(dummy_tmp_1)
  ABI_DEALLOCATE(dummy_tmp_2)
end if

write(io_unit_log,10) '   - Compute the GAMMA matrix...'
flush(io_unit_log)


! Compute the GAMMA matrices
do l1 = 1, lmax

! psik_wrk = | (V^{1/2} . l )^* phi_v >
psik_wrk(1,:) = dble (local_Lbasis(:,l1) )
psik_wrk(2,:) = dimag(local_Lbasis(:,l1))


do l2 = 1, lsolutions

psi_gamma_l2(1,:) = real (QR_orthonormal_basis(:,l2))
psi_gamma_l2(2,:) = dimag(QR_orthonormal_basis(:,l2))


! Note that G_{lJ} = < l | Vc^{1/2}.(gamma_J^*.phi_v)>
!                  = < gamma_J | (Vc^{1/2}.l^*).phi_v >

z(:) = cg_zdotc(npw_k,psi_gamma_l2, psik_wrk)
call xmpi_sum(z,mpi_communicator,ierr) ! sum on all processors working on FFT!
sternheimer_G(l1,l2) = cmplx_1*z(1)+cmplx_i*z(2)

end do ! l1
end do ! l2



if ( write_debug ) then
  write(io_unit,19)   '<gamma| gamma>', sqrt(sum(abs(check_matrix(:,:))**2))
  write(io_unit,19)   '  M hermitian ',sqrt(sum(abs(check_matrix2(:,:))**2))
  write(io_unit,10)   " "
  write(io_unit,10)   "# GAMMA Matrix:"
  write(io_unit,10)   " "

  do l1 = 1, lmax
  write(io_unit,30) sternheimer_G(l1,:)
  end do
  flush(io_unit)
end if

!--------------------------------------------------------------------------
! Step 7: Compute the solutions
!--------------------------------------------------------------------------

write(io_unit_log,10) '   - Compute the Projected Sternheimer solutions and build approximate dielectric operator...'
flush(io_unit_log)



ABI_ALLOCATE(sternheimer_A ,(lsolutions,lsolutions))
ABI_ALLOCATE(sternheimer_X       ,(lsolutions,lmax))
do iw = 2, npt_gauss + 1

write(io_unit_log,23) '        -- iw = ',iw,' / ',npt_gauss+1
flush(io_unit_log)


omega = list_omega(iw)

sternheimer_A(:,:) = sternheimer_A0(:,:) 

do l = 1, lsolutions
sternheimer_A(l,l) = sternheimer_A(l,l) + omega**2
end do

sternheimer_X(:,:) = sternheimer_B(:,:) 

!--------------------------------------------------------------------------
! Step 2:  solve A*X = B, a projected form of the Sternheimer equation
!--------------------------------------------------------------------------
write(io_unit_log,10) '        -- Solve A*X = B'
flush(io_unit_log)



call cpu_time(time1)
call zgesv(lsolutions,      & ! number of rows of A matrix
lmax,            & ! number of columns of B matrix
sternheimer_A,   & ! The A matrix on input, the LU factorization on output
lsolutions,      & ! leading dimension of A
ipiv,            & ! array of pivots
sternheimer_X,   & ! B matrix on input, solution X on output
lsolutions,      & ! leading dimension of B
info )
call cpu_time(time2)

time_exact = time_exact + time2-time1

!--------------------------------------------------------------------------
! Step 3: Add contribution to projected epsilon
!--------------------------------------------------------------------------
! perform E = E  -4 sum_l Gamma_v*X^*_v

write(io_unit_log,10) '        -- add contribution to dielectric matrix at this frequency'
flush(io_unit_log)



ABI_ALLOCATE(dummy_tmp_1,(lsolutions,lmax))
dummy_tmp_1(:,:) = conjg(sternheimer_X) ! DO THIS to avoid potential stack problems

call cpu_time(time1)
call zgemm(     'N',    & ! A matrix is in normal order
'N',    & ! B matrix is in normal order
lmax,    & ! number of rows of A
lmax,    & ! number of columns of B
lsolutions,    & ! number of columns of A
-4.0_dp*cmplx_1,    & ! premultiply A*B by this scalar
sternheimer_G,    & ! GAMMA matrix
lmax,    & ! leading dimension of A
dummy_tmp_1,    & ! B matrix
lsolutions,    & ! leading dimension of B
cmplx_1,    & ! beta  is one
projected_dielectric_Lanczos_basis(:,:,iw),    & ! C matrix
lmax)      ! leading dimension of C

ABI_DEALLOCATE(dummy_tmp_1)

call cpu_time(time2)
time_exact = time_exact + time2-time1

end do ! iw

write(io_unit_log,10) '   - Deallocate tmp arrays...'
flush(io_unit_log)




ABI_DEALLOCATE(ipiv           )
ABI_DEALLOCATE(sternheimer_A  )
ABI_DEALLOCATE(sternheimer_A0 )
ABI_DEALLOCATE(sternheimer_B  )
ABI_DEALLOCATE(sternheimer_X  )
ABI_DEALLOCATE(sternheimer_G  )

end do ! v

!--------------------------------------------------------------------------
! Finalize, post v loop
!--------------------------------------------------------------------------
write(io_unit_log,10) " - Finalize, after band iterations...."
flush(io_unit_log)





timing_string = "#        Exact Sector :   "
call write_timing_log(timing_string,time_exact)


ABI_ALLOCATE(dummy_tmp_1,(lmax,lmax))
ABI_ALLOCATE(dummy_tmp_2,(lmax,lmax))

do iw = 2, npt_gauss+1
! finally, make sure matrix is hermitian

! play this little game to avoid STACK problems
dummy_tmp_1(:,:) = 0.5_dp*transpose(projected_dielectric_Lanczos_basis(:,:,iw))
dummy_tmp_2(:,:) = conjg(dummy_tmp_1(:,:))
dummy_tmp_1(:,:) = dummy_tmp_2(:,:) +0.5_dp*projected_dielectric_Lanczos_basis(:,:,iw)

projected_dielectric_Lanczos_basis(:,:,iw) =    dummy_tmp_1(:,:)
end do

ABI_DEALLOCATE(dummy_tmp_1)
ABI_DEALLOCATE(dummy_tmp_2)



if ( write_debug ) then
  !write some results to a file

  write(io_unit,10) '#===================================================================================================='
  write(io_unit,10) "#                     Projected dielectric matrices                                                  "
  write(io_unit,10) "#                     -------------------------------------------------------                        "
  write(io_unit,10) "# This file contains the various projected dielectric matrices as a function of frequency.           "
  write(io_unit,10) '#===================================================================================================='
  write(io_unit,10) ''
  flush(io_unit)

  do iw = 1, npt_gauss+1
  write(io_unit,10) "#"
  write(io_unit,12) "# omega = ",list_omega(iw), " i Ha"
  write(io_unit,10) "#"

  do l =1, lmax
  write(io_unit,30)  projected_dielectric_Lanczos_basis(l,:,iw)
  end do 
  end do 

end if


if ( model ) then


  ! Add all components on the processors
  call xmpi_sum(model_dielectric_Lanczos_basis,mpi_communicator,ierr) ! sum on all processors

  ! add the identity
  do iw = 1, npt_gauss+1
  do l= 1, lmax
  model_dielectric_Lanczos_basis(l,l,iw) = model_dielectric_Lanczos_basis(l,l,iw) + cmplx_1
  end do
  end do 

  ! hermitian the operator
  do iw = 1, npt_gauss+1

  do l1 = 1, lmax
  do l2 = 1, l1
  ! operator is hermitian
  model_dielectric_Lanczos_basis(l2,l1,iw) = conjg(model_dielectric_Lanczos_basis(l1,l2,iw))
  end do
  end do

  end do ! iw

end if



if ( write_debug ) then
  close(io_unit)
end if

write(io_unit_log,10) " - Deallocate and exit...."
flush(io_unit_log)



if (debug) then
  ABI_DEALLOCATE(check_matrix)
  ABI_DEALLOCATE(check_matrix2)
end if


ABI_DEALLOCATE(psikg_valence)
ABI_DEALLOCATE(psir_valence)
ABI_DEALLOCATE(psikg_VL)


ABI_DEALLOCATE(YL)

ABI_DEALLOCATE(local_Lbasis_conjugated)
ABI_DEALLOCATE(local_Lbasis)


ABI_DEALLOCATE(svd_matrix)
ABI_DEALLOCATE(svd_values)
ABI_DEALLOCATE(psi_gamma_l1)
ABI_DEALLOCATE(psi_gamma_l2)

ABI_DEALLOCATE(psikg_in)
ABI_DEALLOCATE(psikg_out)

ABI_DEALLOCATE(psi_rhs)

ABI_DEALLOCATE(c_sternheimer_solutions)
ABI_DEALLOCATE(QR_orthonormal_basis)

ABI_DEALLOCATE(psik_wrk)
ABI_DEALLOCATE(psikb_wrk)
ABI_DEALLOCATE(psikg_wrk)


close(io_unit_log)


10 format(A)
12 format(A,F12.8,A)
14 format(I5,10X,ES24.12)
15 format(I5,10X,ES24.12,10X,A)

19 format(20X,A,15X,E24.16)

20 format(A,I5)
21 format(I5,10X,F8.4,15X,ES24.12)
22 format(A,I5,A)
23 format(A,I5,A,I5)
30 format(2X,1000(ES12.4,2X,ES12.4,5X))

end subroutine ProjectedSternheimerEpsilon
!!***


end module m_gwls_DielectricArray
!!***

!!****m* ABINIT/m_gwls_LanczosBasis
!! NAME
!! m_gwls_LanczosBasis
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


module m_gwls_LanczosBasis
!----------------------------------------------------------------------------------------------------
! This module contains the static Lanczos basis, which should be computed once and for all,
! and then be made available to other modules.
!----------------------------------------------------------------------------------------------------
! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_hamiltonian
use m_gwls_GenerateEpsilon
use m_dtset

! Abinit modules
use defs_basis
use m_abicore

implicit none
save
private
!!***

! Global arrays

! basis which diagonalizes the static dielectric matrix
complex(dpc), public, allocatable :: Lbasis_lanczos(:,:)  ! complex array which contains the Lanczos basis

! basis which diagonalizes the model static dielectric matrix
complex(dpc), public, allocatable :: Lbasis_model_lanczos(:,:)  ! complex array which contains the Lanczos basis


!------------------------------------------------------------
!
! OBSOLETE STRUCTURES, KEPT AROUND TO NOT BREAK THE CODE
!    Some cleaning will eventually have to be done!
!
!------------------------------------------------------------
! Lanczos basis which puts the static dielectric matrix in bands (OBSOLETE)
real(dp), public, allocatable :: lanczos_basis_0(:,:,:,:)

! modified basis, of the form (V^{1/2}. l^*) psie (OBSOLETE)
real(dp), public, allocatable :: basis_0(:,:,:,:)


! modified basis, of the form (V^{1/2}. l^*) psie
complex(dpc), public, allocatable :: Lbasis_modified(:,:)


integer, public  :: lanczos_basis_size
!!***

public :: setup_Lanczos_basis
public :: cleanup_Lanczos_basis
public :: modify_Lbasis_Coulomb
!!***
contains

!!****f* m_hamiltonian/setup_Lanczos_basis
!! NAME
!!  setup_Lanczos_basis
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputeCorrelationEnergy
!!
!! CHILDREN
!!      g_to_r,gr_to_g,sqrt_vc_k,wf_block_distribute
!!
!! SOURCE

subroutine setup_Lanczos_basis(lmax,lmax_model)
!----------------------------------------------------------------------------------------------------
! Set up the lanczos basis
!----------------------------------------------------------------------------------------------------

integer, intent(in) :: lmax, lmax_model

! *************************************************************************

ABI_ALLOCATE(Lbasis_lanczos,(npw_k,lmax))

if (lmax_model > 0) then
  ABI_ALLOCATE(Lbasis_model_lanczos,(npw_k,lmax_model))
end if

end subroutine setup_Lanczos_basis
!!***

!!****f* m_hamiltonian/cleanup_Lanczos_basis
!! NAME
!!  cleanup_Lanczos_basis
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputeCorrelationEnergy
!!
!! CHILDREN
!!      g_to_r,gr_to_g,sqrt_vc_k,wf_block_distribute
!!
!! SOURCE

subroutine cleanup_Lanczos_basis()


! *************************************************************************

! if(allocated()) ABI_DEALLOCATE() can cause a segfault if used in this form, without the THEN.
! This is because ABI_ALLOCATE is expanded at compilation time in many statements; and the if() can only prevent the execution of
! the first.
! So, the deallocation takes place even if allocated() returns .false. without the THEN.
if (allocated(Lbasis_lanczos))  then
  ABI_DEALLOCATE(Lbasis_lanczos)
end if
if (allocated(Lbasis_model_lanczos))  then
  ABI_DEALLOCATE(Lbasis_model_lanczos)
end if

end subroutine cleanup_Lanczos_basis
!!***

!!****f* m_hamiltonian/modify_Lbasis_Coulomb
!! NAME
!!  modify_Lbasis_Coulomb
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_ComputeCorrelationEnergy
!!
!! CHILDREN
!!      g_to_r,gr_to_g,sqrt_vc_k,wf_block_distribute
!!
!! SOURCE

subroutine modify_Lbasis_Coulomb(psie_k, lmax, lmax_model)
!----------------------------------------------------------------------------------------------------
! This subroutine computes, once and for all, the vectors (V^1/2.L)^*.psie
!----------------------------------------------------------------------------------------------------
real(dp),intent(in)  :: psie_k(2,npw_k)
integer ,intent(in)  :: lmax, lmax_model


real(dp), allocatable  :: psik_wrk(:,:), psikb_wrk(:,:), psikg_wrk(:,:)
real(dp), allocatable  :: psikg_e(:,:)

integer  :: iblk_lanczos, nbdblock_lanczos

integer   :: l, mb

! *************************************************************************

ABI_ALLOCATE(psik_wrk  ,(2,npw_k))
ABI_ALLOCATE(psikb_wrk ,(2,npw_kb))
ABI_ALLOCATE(psikg_wrk ,(2,npw_g))
ABI_ALLOCATE(psikg_e ,(2,npw_g))

! copy the valence state on every row of FFT processors
do mb = 1, blocksize
psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psie_k(:,:)
end do ! mb
! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_e,1) ! LA -> FFT


! Number of blocks of lanczos vectors
nbdblock_lanczos = lmax/blocksize
if (modulo(lmax,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1

! loop on all blocks of lanczos vectors
do iblk_lanczos = 1, nbdblock_lanczos

! loop on all states within this block
do mb = 1, blocksize
! Determine the index of the Lanczos vector
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax) then
  psik_wrk(1,:) = dble (Lbasis_lanczos(:,l))
  psik_wrk(2,:) = dimag(Lbasis_lanczos(:,l))
else
  psik_wrk(:,:) = zero
end if
! Apply coulomb potential
call sqrt_vc_k(psik_wrk)

! Store in array of blocks of wavefunctions
psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)

end do ! mb

! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_wrk,1) ! LA -> FFT

! fourier transform, conjugate
call g_to_r(psir1,psikg_wrk)
psir1(2,:,:,:) = -psir1(2,:,:,:)

! compute the product with the state "e"
call gr_to_g(psikg_wrk, psir1, psikg_e)

! return to LA configuration

! Transform to LA representation
call wf_block_distribute(psikb_wrk,  psikg_wrk, 2) ! FFT -> LA

do mb = 1, blocksize
! Determine the index of the Lanczos vector
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax) then
  psik_wrk(:,:) = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)
  Lbasis_lanczos(:,l) = cmplx_1*psik_wrk(1,:)+cmplx_i*psik_wrk(2,:)
end if

end do ! mb

end do ! iblk_lanczos

! repeat, for model basis!
! Number of blocks of lanczos vectors
nbdblock_lanczos = lmax_model/blocksize
if (modulo(lmax_model,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1

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
! Apply coulomb potential
call sqrt_vc_k(psik_wrk)

! Store in array of blocks of wavefunctions
psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:)

end do ! mb

! Transform to FFT representation
call wf_block_distribute(psikb_wrk,  psikg_wrk,1) ! LA -> FFT

! fourier transform, conjugate
call g_to_r(psir1,psikg_wrk)
psir1(2,:,:,:) = -psir1(2,:,:,:)

! compute the product with the state "e"
call gr_to_g(psikg_wrk, psir1, psikg_e)

! return to LA configuration
! Transform to LA representation
call wf_block_distribute(psikb_wrk,  psikg_wrk, 2) ! FFT -> LA

do mb = 1, blocksize
! Determine the index of the Lanczos vector
l = (iblk_lanczos-1)*blocksize + mb

if ( l <= lmax_model) then
  psik_wrk(:,:) = psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k)
  Lbasis_model_lanczos(:,l) = cmplx_1*psik_wrk(1,:)+cmplx_i*psik_wrk(2,:)
end if

end do ! mb

end do

ABI_DEALLOCATE(psik_wrk  )
ABI_DEALLOCATE(psikb_wrk )
ABI_DEALLOCATE(psikg_wrk )
ABI_DEALLOCATE(psikg_e )


end subroutine modify_Lbasis_Coulomb
!!***

end module m_gwls_LanczosBasis
!!***

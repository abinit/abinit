!!****m* ABINIT/m_gwls_model_polarisability
!! NAME
!! m_gwls_model_polarisability
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


module m_gwls_model_polarisability

! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_valenceWavefunctions
use m_gwls_hamiltonian
use m_gwls_lineqsolver

! abinit modules
use defs_basis
use m_abicore
use m_bandfft_kpt
use m_errors

use m_time,      only : timab

implicit none
save
private
!!***

real(dp), public :: model_polarizability_epsilon_0  ! model parameter
!integer          :: model_polarizability_model_type ! how is epsilon_0 used to model?


real(dp), allocatable, private :: psir_model(:,:,:,:), psir_ext_model(:,:,:,:)


integer, public  :: dielectric_model_type = 1

real(dp),public, allocatable  :: model_Y(:) ! model susceptibility, distributed according to FFT configuration
real(dp),public, allocatable  :: model_Y_LA(:) ! model susceptibility, distributed according to LA configuration

real(dp),public, allocatable  :: sqrt_density(:,:,:,:) ! average valence wave function...
!!***

public :: Pk_model


public :: epsilon_k_model
public :: setup_Pk_model
public :: cleanup_Pk_model

public :: matrix_function_epsilon_model_operator
!!***

contains

!!****f* m_hamiltonian/epsilon_k_model
!! NAME
!!  epsilon_k_model
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_model_polarisability
!!
!! CHILDREN
!!      epsilon_k_model
!!
!! SOURCE

subroutine epsilon_k_model(psi_out,psi_in)

implicit none

real(dp), intent(out) :: psi_out(2,npw_k)
real(dp), intent(in)  :: psi_in(2,npw_k)

real(dp) :: psik(2,npw_k)

! *************************************************************************

psik = psi_in
call sqrt_vc_k(psik)
call Pk_model(psi_out ,psik)
call sqrt_vc_k(psi_out)
psi_out = psi_in - psi_out

end subroutine epsilon_k_model
!!***

!!****f* m_hamiltonian/setup_Pk_model
!! NAME
!!  setup_Pk_model
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_ComputeCorrelationEnergy,gwls_DielectricArray,gwls_GenerateEpsilon
!!
!! CHILDREN
!!      epsilon_k_model
!!
!! SOURCE

subroutine setup_Pk_model(omega,epsilon_0)
!---------------------------------------------------------------
!
! This subroutine prepares a global array in order to 
! act with the model susceptibility, given by
!
! Pk_model(r,r',i omega) = sum_v phi_v(r) Y(r-r',i omega) phi^*_v(r')
!
!
! This subroutine computes the Fourier transform of Y(r,i omega),
! Y(G,omega), for a given omega and epsilon_0.
!
! It is assumed that omega and epsilon_0 >= 0. Also, the model
! describes an IMAGINARY frequency i omega.
!---------------------------------------------------------------
real(dp), intent(in) :: epsilon_0, omega

real(dp) ::  theta, R_omega
real(dp) ::  x, y


integer  :: ig

real(dp),allocatable ::  G_array(:)

! *************************************************************************

if (.not. allocated(psir_model)) then
  ABI_ALLOCATE(psir_model, (2,n4,n5,n6))
end if

if (.not. allocated(psir_ext_model)) then
  ABI_ALLOCATE(psir_ext_model, (2,n4,n5,n6))
end if

R_omega = 2.0_dp*sqrt(epsilon_0**2+omega**2)
if(abs(omega) > tol16 .or. abs(epsilon_0) > tol16) then 
  theta   = atan2(omega,epsilon_0)
else
  theta   = zero
end if

x = sqrt(R_omega)*cos(0.5_dp*theta)
y = sqrt(R_omega)*sin(0.5_dp*theta)


if (.not. allocated(model_Y)) then
  ABI_ALLOCATE(model_Y, (npw_g))
end if

if (.not. allocated(model_Y_LA)) then
  ABI_ALLOCATE(model_Y_LA, (npw_k))
end if

!================================================================================
! Compute model_Y, in FFT configuration
!================================================================================

ABI_ALLOCATE(G_array,(npw_g))
G_array(:) = sqrt(2.0_dp*kinpw_gather(:))

model_Y(:) = zero
do ig = 1, npw_g


if (G_array(ig) > tol12) then
  ! G != 0.
  model_Y(ig) = -4.0_dp/G_array(ig)*                             &
  ((G_array(ig)+y)/((G_array(ig)+y)**2+x**2)     &
  + (G_array(ig)-y)/((G_array(ig)-y)**2+x**2))


else
  if ( abs(epsilon_0) < tol12 ) then
    model_Y(ig) = zero
  else 
    model_Y(ig) = -4.0_dp*epsilon_0/(epsilon_0**2+omega**2)
  end if
end if
end do ! ig

ABI_DEALLOCATE(G_array)

!================================================================================
! Compute model_Y_LA, in LA configuration
!================================================================================

ABI_ALLOCATE(G_array,(npw_k))
G_array(:) = sqrt(2.0_dp*kinpw(:))

model_Y_LA(:) = zero
do ig = 1, npw_k

if (G_array(ig) > tol12) then
  ! G != 0.
  model_Y_LA(ig) = -4.0_dp/G_array(ig)*                          &
  ((G_array(ig)+y)/((G_array(ig)+y)**2+x**2)     &
  + (G_array(ig)-y)/((G_array(ig)-y)**2+x**2))


else
  if ( abs(epsilon_0) < tol12 ) then
    model_Y_LA(ig) = zero
  else 
    model_Y_LA(ig) = -4.0_dp*epsilon_0/(epsilon_0**2+omega**2)
  end if
end if
end do ! ig

ABI_DEALLOCATE(G_array)



if (dielectric_model_type == 2) then

  MSG_BUG('dielectric_model_type == 2 not properly implemented. Review code or input!')

  !ABI_ALLOCATE(sqrt_density,(2,n4,n5,n6))
  !sqrt_density(:,:,:,:) = zero
  !do v= 1, nbandv
  !        sqrt_density(1,:,:,:) = sqrt_density(1,:,:,:) + valence_wfr(1,:,:,:,v)**2+valence_wfr(2,:,:,:,v)**2
  !end do 
  !sqrt_density(1,:,:,:) = sqrt(sqrt_density(1,:,:,:))

end if 


end subroutine setup_Pk_model
!!***

!!****f* m_hamiltonian/cleanup_Pk_model
!! NAME
!!  cleanup_Pk_model
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
!!      epsilon_k_model
!!
!! SOURCE

subroutine cleanup_Pk_model()

implicit none

! *************************************************************************

if (allocated(model_Y)) then
  ABI_DEALLOCATE(model_Y)
end if

if (allocated(model_Y_LA)) then
  ABI_DEALLOCATE(model_Y_LA)
end if



if (allocated(sqrt_density)) then
  ABI_DEALLOCATE(sqrt_density)
end if

if ( allocated(psir_model)) then
  ABI_DEALLOCATE(psir_model)
end if

if (allocated(psir_ext_model)) then
  ABI_DEALLOCATE(psir_ext_model)
end if

end subroutine cleanup_Pk_model
!!***

!!****f* m_hamiltonian/Pk_model
!! NAME
!!  Pk_model
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_model_polarisability
!!
!! CHILDREN
!!      epsilon_k_model
!!
!! SOURCE

subroutine Pk_model(psi_out,psi_in)
!------------------------------------------------------------------------------------------------------------------------
! Returns the action of a frequency-dependent model susceptibility
!------------------------------------------------------------------------------------------------------------------------
implicit none

real(dp), intent(out) :: psi_out(2,npw_k)
real(dp), intent(in)  :: psi_in(2,npw_k)

! *************************************************************************

if ( dielectric_model_type == 1) then
  call Pk_model_implementation_1(psi_out ,psi_in)
else if ( dielectric_model_type == 2) then
!  call Pk_model_implementation_2(psi_out ,psi_in)
end if

end subroutine Pk_model
!!***

!!****f* m_hamiltonian/Pk_model_implementation_1
!! NAME
!!  Pk_model_implementation_1
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_model_polarisability
!!
!! CHILDREN
!!      epsilon_k_model
!!
!! SOURCE

subroutine Pk_model_implementation_1(psi_out,psi_in)
!------------------------------------------------------------------------------------------------------------------------
! Returns the action of a frequency-dependent model susceptibility
!------------------------------------------------------------------------------------------------------------------------
implicit none

real(dp), intent(out) :: psi_out(2,npw_k)
real(dp), intent(in)  :: psi_in(2,npw_k)

integer :: v

real(dp), allocatable ::   psik(:,:), psik_g(:,:)

integer, save  :: icounter = 0

integer  :: mb, iblk

integer  :: mpi_band_rank    

real(dp) :: time1, time2
real(dp) :: total_time1, total_time2

real(dp), save :: fft_time        = zero
real(dp), save :: projection_time = zero
real(dp), save :: Y_time          = zero
real(dp), save :: total_time      = zero


real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB

! *************************************************************************

GWLS_TIMAB   = 1534
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



call cpu_time(total_time1)
icounter = icounter + 1

GWLS_TIMAB   = 1535
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


ABI_ALLOCATE(psik,           (2,npw_kb))
ABI_ALLOCATE(psik_g,         (2,npw_g))

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


! initialize the output to zero
psi_out = zero

! MPI information
mpi_band_rank = mpi_enreg%me_band


!-----------------------------------------------------------------
! Put a copy of the external state on every row of FFT processors.
!
! The, copy conjugate of initial wavefunction in local array,
! and set inout array to zero.
!-----------------------------------------------------------------
call cpu_time(time1)
GWLS_TIMAB   = 1536
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

! fill the array psik_ext with copies of the external state
do mb = 1, blocksize
psik(:,(mb-1)*npw_k+1:mb*npw_k)   = psi_in(:,:)
end do

! change configuration of the data, from LA to FFT
call wf_block_distribute(psik,  psik_g,1) ! LA -> FFT
! Now every row of FFT processors has a copy of the external state.


! Store external state in real-space format, inside module-defined work array psir_ext_model
! Each row of FFT processors will have a copy!
call g_to_r(psir_ext_model,psik_g)

! Conjugate the external wavefunction; result will be conjugated again later,
! insuring we are in fact acting on psi, and not psi^*. This conjugation
! is only done for algorithmic convenience.
psir_ext_model(2,:,:,:) = -psir_ext_model(2,:,:,:)


OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
fft_time = fft_time+time2-time1


! Loop on all blocks of eigenstates
do iblk = 1, nbdblock

v = (iblk-1)*blocksize + mpi_band_rank + 1 ! CAREFUL! This is a guess. Revisit this if code doesn't work as intended. 

call cpu_time(time1)
GWLS_TIMAB   = 1537
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


! Multiply valence state by external state, yielding  |psik_g> = | phi_v x psi_in^* > 
call gr_to_g(psik_g, psir_ext_model, valence_wavefunctions_FFT(:,:,iblk))

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
fft_time = fft_time+time2-time1

! Project out to conduction space
call cpu_time(time1)
GWLS_TIMAB   = 1538
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call pc_k_valence_kernel(psik_g)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
projection_time = projection_time+time2-time1


! act with model susceptibility
call cpu_time(time1)
GWLS_TIMAB   = 1539
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


psik_g(1,:)  = psik_g(1,:)*model_Y(:)
psik_g(2,:)  = psik_g(2,:)*model_Y(:)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
Y_time = Y_time+time2-time1

! Project out to conduction space, again!
call cpu_time(time1)
GWLS_TIMAB   = 1538
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call pc_k_valence_kernel(psik_g)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)
call cpu_time(time2)

projection_time= projection_time+time2-time1

! Express result in real space, in module-defined work array psir_model
call cpu_time(time1)
GWLS_TIMAB   = 1536
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call g_to_r(psir_model,psik_g)
! conjugate the result, cancelling the initial conjugation described earlier.
psir_model(2,:,:,:) = -psir_model(2,:,:,:)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
fft_time = fft_time+time2-time1

!  Multiply by valence state in real space 
call cpu_time(time1)
GWLS_TIMAB   = 1537
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call gr_to_g(psik_g,psir_model, valence_wavefunctions_FFT(:,:,iblk))

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)
call cpu_time(time2)
fft_time = fft_time+time2-time1



! Return to linear algebra format, and add condtribution
GWLS_TIMAB   = 1540
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call wf_block_distribute(psik,  psik_g,2) ! FFT -> LA

do mb = 1, blocksize

v = (iblk-1)*blocksize + mb
if ( v <= nbandv) then
  psi_out(:,:) = psi_out(:,:) + psik(:,(mb-1)*npw_k+1:mb*npw_k)
end if
end do

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


end do ! iblk

call cpu_time(total_time2)
total_time = total_time + total_time2-total_time1

GWLS_TIMAB   = 1535
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

ABI_DEALLOCATE(psik)
ABI_DEALLOCATE(psik_g)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


GWLS_TIMAB   = 1534
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


end subroutine Pk_model_implementation_1 
!!***

!!****f* m_hamiltonian/matrix_function_epsilon_model_operator
!! NAME
!!  matrix_function_epsilon_model_operator
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_GenerateEpsilon
!!
!! CHILDREN
!!      epsilon_k_model
!!
!! SOURCE

subroutine matrix_function_epsilon_model_operator(vector_out,vector_in,Hsize)
!----------------------------------------------------------------------------------------------------
! This function returns the action of the operator epsilon_model on a given vector.
! It is assumed that the frequency has been set in the module using setup_Pk_model.
!    
!
!----------------------------------------------------------------------------------------------------
implicit none
integer,      intent(in)  :: Hsize
complex(dpc), intent(out) :: vector_out(Hsize)
complex(dpc), intent(in)  :: vector_in(Hsize)

! local variables
real(dp)     :: psik (2,Hsize)
real(dp)     :: psik2(2,Hsize)

! *************************************************************************

! convert from one format to the other
psik(1,:) = dble (vector_in(:))
psik(2,:) = dimag(vector_in(:))

call epsilon_k_model(psik2 ,psik)

! Act with  epsilon_model
vector_out = cmplx_1*psik2(1,:)+cmplx_i*psik2(2,:)

end subroutine matrix_function_epsilon_model_operator
!!***


end module m_gwls_model_polarisability
!!***

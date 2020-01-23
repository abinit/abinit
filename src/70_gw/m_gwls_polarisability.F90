!!****m* ABINIT/m_gwls_polarisability
!! NAME
!! m_gwls_polarisability
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


module m_gwls_polarisability
! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_valenceWavefunctions
use m_gwls_hamiltonian
use m_gwls_lineqsolver

! abinit modules
use defs_basis
use m_errors
use m_abicore
use m_bandfft_kpt

use m_time,      only : timab

implicit none
save
private
!!***

real(dp),public :: matrix_function_omega(2)


! Some timing variables
integer,  public :: counter_fft, counter_sqmr, counter_rprod, counter_proj, counter_H

real(dp), public :: time1, time2, time_fft, time_sqmr, time_rprod,time_proj,time_H

real(dp), allocatable, public :: Sternheimer_solutions_zero(:,:,:,:)
integer, public :: index_solution=0
integer, public :: recy_unit
logical, public :: write_solution=.false.

!integer          :: io_unit
!!***

public :: Pk, epsilon_k
public :: matrix_function_epsilon_k
public :: set_dielectric_function_frequency
!!***

contains

!!****f* m_hamiltonian/Pk
!! NAME
!!  Pk
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_polarisability
!!
!! CHILDREN
!!      epsilon_k
!!
!! SOURCE

subroutine Pk(psi_inout,omega)
!===============================================================================
! 
! This routine applies the polarizability operator to an arbitrary state psi_inout.
!
!
! A note about parallelism: 
! -------------------------
!  The input/output is in "linear algebra" configuration, which is to say
!  that ALL processors have a fraction of the G-vectors for this state. Internally,
!  this routine will parallelise over bands and FFT, thus using the "FFT" 
!  configuration of the data. For more information, see the lobpcgwf.F90, and 
!  the article by F. Bottin et al.
!===============================================================================
implicit none

real(dp), intent(inout) :: psi_inout(2,npw_k)
real(dp), intent(in)  :: omega(2)

logical :: omega_imaginary

integer :: v, mb, iblk

real(dp):: norm_omega
real(dp), allocatable ::   psik(:,:)
real(dp), allocatable ::   psik_alltoall(:,:),    psik_wrk_alltoall(:,:)
real(dp), allocatable ::   psik_in_alltoall(:,:), psik_tmp_alltoall(:,:)

real(dp), allocatable ::   psik_ext(:,:), psik_ext_alltoall(:,:)


real(dp), allocatable ::   psir(:,:,:,:), psir_ext(:,:,:,:) 

integer         :: cplex
integer :: recy_i


integer      :: mpi_band_rank

real(dp) ::  list_SQMR_frequencies(2)
real(dp) ::  list_QMR_frequencies(2,2)


integer,  save ::  icounter = 0
real(dp), save ::  total_time1, total_time2, total_time


integer :: num_op_v, i_op_v, case_op_v
real(dp):: factor_op_v 
real(dp):: lbda
real(dp):: zz(2)

character(len=500) :: message

real(dp) :: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB

! *************************************************************************

GWLS_TIMAB   = 1524
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



icounter = icounter + 1
call cpu_time(total_time1)


!========================================
! Allocate work arrays and define 
! important parameters
!========================================
GWLS_TIMAB   = 1525
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


cplex  = 2 ! complex potential

mpi_band_rank    = mpi_enreg%me_band

ABI_ALLOCATE(psik,                (2,npw_kb))
ABI_ALLOCATE(psik_alltoall,       (2,npw_g))
ABI_ALLOCATE(psik_wrk_alltoall,   (2,npw_g))
ABI_ALLOCATE(psik_tmp_alltoall,   (2,npw_g))
ABI_ALLOCATE(psik_in_alltoall,    (2,npw_g))

ABI_ALLOCATE(psir,    (2,n4,n5,n6))
ABI_ALLOCATE(psir_ext,(2,n4,n5,n6))


OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


!--------------------------------------------------------------------------------
! 
! The polarizability acting on a state | PSI > is given by
!
!     Pk | PSI >  = 2 sum_{v} | h_v^* phi_v >    ( factor of 2 comes from sum on spin)
!
!        | h_v >  = Pc . OPERATOR_v . Pc | PSI^* phi_v >
!
!     There are multiple cases to consider:
!
!     CASE:
!                 1)  |omega| = 0
!                                 OPERATOR_v =  1/[H-ev]
!                                 prefactor  = -4
!
!                 2)   omega  =  i lbda (imaginary)
!                                 OPERATOR_v =   [H-ev]/(lbda^2+[H-ev]^2)
!                                 prefactor  = -4
!
!                 3)   omega  =    lbda (real) (USE SQMR)
!                                 OPERATOR_v =  { 1/[H-ev+lbda]+ 1/[H-ev-lbda] }
!                                 prefactor  = -2
!
!                 4)   omega  =    lbda (real) (USE QMR)
!                                 OPERATOR_v =  { 1/[H-ev+lbda]+ 1/[H-ev-lbda] }
!                                 prefactor  = -2
!
! It simplifies the code below to systematize the algorithm.
!
!--------------------------------------------------------------------------------

! Check which part of omega is non-zero. Default is that omega is real.
if (abs(omega(2)) < 1.0d-12) then
  omega_imaginary=.false.
  norm_omega     = omega(1)

  elseif (abs(omega(1)) < 1.0d-12 .and. abs(omega(2)) > 1.0d-12) then
  omega_imaginary = .true.
  norm_omega      = omega(2)

else
  write(message,"(a,es16.8,3a)")&
    "omega=",omega,",",ch10,&
    "but either it's real or imaginary part need to be 0 for the polarisability routine to work."
  MSG_ERROR(message)
end if


!-----------------------------------------------------------------
! I) Prepare global values depending on the CASE being considered
!-----------------------------------------------------------------

if (norm_omega < 1.0D-12 ) then
  case_op_v = 1
  num_op_v  = 1
  factor_op_v  = -4.0_dp

else if (omega_imaginary) then
  case_op_v = 2
  num_op_v  = 1
  factor_op_v  = -4.0_dp

else if( .not. activate_inf_shift_poles) then
  case_op_v = 3
  num_op_v  = 2
  factor_op_v  = -2.0_dp

else 
  case_op_v = 4
  num_op_v  = 2
  factor_op_v  = -2.0_dp

  inf_shift_poles = dtset%zcut


  write(message,*) " inf_shift_poles = ",inf_shift_poles
  call wrtout(std_out,message,'COLL')
end if


write(message,10)" "
call wrtout(std_out,message,'COLL')
write(message,10) "     Pk: applying the polarizability on states"
call wrtout(std_out,message,'COLL')
write(message,10) "     =============================================================="
call wrtout(std_out,message,'COLL')
write(message,12) "     CASE : ",case_op_v
call wrtout(std_out,message,'COLL')



!-----------------------------------------------------------------
! II) copy conjugate of initial wavefunction in local array,
!     and set inout array to zero. Each FFT row of processors
!     must have a copy of the initial wavefunction in FFT 
!     configuration!
!-----------------------------------------------------------------
call cpu_time(time1)

GWLS_TIMAB   = 1526
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


ABI_ALLOCATE(psik_ext,(2,npw_kb))
ABI_ALLOCATE(psik_ext_alltoall,(2,npw_g))

! fill the array psik_ext with copies of the external state
do mb = 1, blocksize
psik_ext(:,(mb-1)*npw_k+1:mb*npw_k)   = psi_inout(:,:)
end do

! change configuration of the data, from LA to FFT
call wf_block_distribute(psik_ext,  psik_ext_alltoall,1) ! LA -> FFT
! Now every row of FFT processors has a copy of the external state.

! Copy the external state to the real space format, appropriate for real space products to be 
! used later.

call g_to_r(psir_ext,psik_ext_alltoall)
psir_ext(2,:,:,:) = -psir_ext(2,:,:,:)


! Don't need these arrays anymore...
ABI_DEALLOCATE(psik_ext)
ABI_DEALLOCATE(psik_ext_alltoall)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time_fft    = time_fft + time2-time1
counter_fft = counter_fft+1 

! set external state to zero, ready to start cumulating the answer
psi_inout = zero

!-----------------------------------------------------------------
! III) Iterate on all valence bands
!-----------------------------------------------------------------

! Loop on all blocks of eigenstates
do iblk = 1, nbdblock


! What is the valence band index for this block and this row of FFT processors? It is not clear from the
! code in lobpcgwf.F90; I'm going to *guess*.

v = (iblk-1)*blocksize + mpi_band_rank + 1 ! CAREFUL! This is a guess. Revisit this if code doesn't work as intended. 


!Solving of Sternheiner equation
write(message,12) "          band            :", v
call wrtout(std_out,message,'COLL')
write(message,14) "          eigenvalue (Ha) :  ",eig(v)
call wrtout(std_out,message,'COLL')
write(message,14) "          Re[omega]  (Ha) :  ",omega(1)
call wrtout(std_out,message,'COLL')
write(message,14) "          Im[omega]  (Ha) :  ",omega(2)
call wrtout(std_out,message,'COLL')

!-----------------------------------------------------------------
! IV) prepare some arrays, if they are needed
!-----------------------------------------------------------------
if      (case_op_v == 3) then

  list_SQMR_frequencies(1) = eig(v) - norm_omega
  list_SQMR_frequencies(2) = eig(v) + norm_omega

else if (case_op_v == 4) then
  list_QMR_frequencies(:,1) = (/eig(v)-norm_omega,-inf_shift_poles/)
  list_QMR_frequencies(:,2) = (/eig(v)+norm_omega, inf_shift_poles/)
end if


!-----------------------------------------------------------------
! V) Compute the real-space product of the input wavefunction 
!    with the valence wavefunction.
!-----------------------------------------------------------------
! Unfortunately, the input wavefunctions will be FFT-transformed k-> r nbandv times
! because fourwf cannot multiply a real-space potential with a real-space wavefunction!
call cpu_time(time1)

GWLS_TIMAB   = 1527
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call gr_to_g(psik_in_alltoall,psir_ext,valence_wavefunctions_FFT(:,:,iblk))

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
time_fft    = time_fft + time2-time1
counter_fft = counter_fft+1 

!-----------------------------------------------------------------
! VI) Project out to conduction space
!-----------------------------------------------------------------
call cpu_time(time1)

GWLS_TIMAB   = 1528
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

!call pc_k(psik_in)
call pc_k_valence_kernel(psik_in_alltoall)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time_proj   = time_proj + time2-time1
counter_proj= counter_proj+1 


!-----------------------------------------------------------------
! VII) Loop on potential valence operators,
!-----------------------------------------------------------------
do i_op_v = 1, num_op_v

if      (case_op_v == 1) then

  ! frequency is zero, operator to apply is 1/[H-ev]
  call cpu_time(time1)
  GWLS_TIMAB   = 1529
  OPTION_TIMAB = 1
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

  call sqmr(psik_in_alltoall,psik_tmp_alltoall,eig(v),1)

  OPTION_TIMAB = 2
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


  call cpu_time(time2)
  time_sqmr   = time_sqmr+ time2-time1
  counter_sqmr= counter_sqmr+1 

  ! If we are constructing the $\hat \epsilon(i\omega = 0)$ matrix (and the Lanczos basis at the same time),
  ! keep the Sternheimer solutions for use in the projected Sternheimer section (in LA configuration).
  if(write_solution .and. ((iblk-1)*blocksize < nbandv)) then
    call wf_block_distribute(psik, psik_tmp_alltoall, 2) ! FFT -> LA
    do mb=1,blocksize
    v = (iblk-1)*blocksize+mb
    if(v <= nbandv) then
      if(dtset%gwls_recycle == 1) then
        Sternheimer_solutions_zero(:,:,index_solution,v) = psik(:,(mb-1)*npw_k+1:mb*npw_k)
      end if
      if(dtset%gwls_recycle == 2) then
        recy_i = (index_solution-1)*nbandv + v
        !BUG : On petrus, NAG 5.3.1 + OpenMPI 1.6.2 cause read(...,rec=i) to read the data written by write(...,rec=i+1). 
        !Workaround compatible only with nag : write(recy_unit,rec=recy_i+1).
        write(recy_unit,rec=recy_i) psik(:,(mb-1)*npw_k+1:mb*npw_k)
      end if
    end if
    end do
  end if

else if (case_op_v == 2) then

  ! frequency purely imaginary, operator to apply is 
  ! [H-ev]/(lbda^2+[H-ev]^2)
  call cpu_time(time1)

  GWLS_TIMAB   = 1533
  OPTION_TIMAB = 1
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

  psik_wrk_alltoall = psik_in_alltoall
  call Hpsik(psik_wrk_alltoall,eig(v))

  OPTION_TIMAB = 2
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


  call cpu_time(time2)
  time_H    = time_H + time2-time1
  counter_H = counter_H+1 

  call cpu_time(time1)

  GWLS_TIMAB   = 1530
  OPTION_TIMAB = 1
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

  call sqmr(psik_wrk_alltoall,psik_tmp_alltoall,eig(v),0,norm_omega,omega_imaginary)

  OPTION_TIMAB = 2
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


  call cpu_time(time2)
  time_sqmr   = time_sqmr+ time2-time1
  counter_sqmr= counter_sqmr+1 

else if (case_op_v == 3) then

  ! frequency purely real, operator to apply is 
  !        1/[H-ev +/- lbda]
  !        TREATED WITH SQMR


  lbda = list_SQMR_frequencies(i_op_v)
  call cpu_time(time1)

  GWLS_TIMAB   = 1531
  OPTION_TIMAB = 1
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

  call sqmr(psik_in_alltoall,psik_tmp_alltoall,lbda,1)

  OPTION_TIMAB = 2
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)



  call cpu_time(time2)
  time_sqmr   = time_sqmr+ time2-time1
  counter_sqmr= counter_sqmr+1 

else if (case_op_v == 4) then
  ! frequency purely real, operator to apply is 
  !        1/[H-ev +/- lbda]
  !        TREATED WITH QMR

  zz(:) = list_QMR_frequencies(:,i_op_v)

  call cpu_time(time1)
  GWLS_TIMAB   = 1532
  OPTION_TIMAB = 1
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

  call qmr(psik_in_alltoall,psik_tmp_alltoall,zz) !,0)

  OPTION_TIMAB = 2
  call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


  call cpu_time(time2)
  time_sqmr   = time_sqmr+ time2-time1
  counter_sqmr= counter_sqmr+1 

end if
!-----------------------------------------------------------------
! VIII) Project on conduction states
!-----------------------------------------------------------------
call cpu_time(time1)
GWLS_TIMAB   = 1528
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call pc_k_valence_kernel(psik_tmp_alltoall)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
time_proj   = time_proj + time2-time1
counter_proj= counter_proj+1 

!-----------------------------------------------------------------
! IX) Conjugate result, and express in denpot format
!-----------------------------------------------------------------

call cpu_time(time1)

GWLS_TIMAB   = 1526
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call g_to_r(psir,psik_tmp_alltoall)
psir(2,:,:,:) = -psir(2,:,:,:)


OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time_fft    = time_fft + time2-time1
counter_fft = counter_fft+1 

!-----------------------------------------------------------------
! X) Multiply by valence state in real space 
!-----------------------------------------------------------------
call cpu_time(time1)
GWLS_TIMAB   = 1527
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call gr_to_g(psik_alltoall,psir,valence_wavefunctions_FFT(:,:,iblk))


OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
time_fft    = time_fft + time2-time1
counter_fft = counter_fft+1 

!-----------------------------------------------------------------
! XI) Return to LA configuration, and cumulate the sum
!-----------------------------------------------------------------
call wf_block_distribute(psik,  psik_alltoall,2) ! FFT -> LA

do mb = 1, blocksize

v = (iblk-1)*blocksize + mb

if (v <= nbandv) then
  ! only add contributions from valence
  psi_inout(:,:) = psi_inout(:,:) + psik(:,(mb-1)*npw_k+1:mb*npw_k)
end if

end do


end do ! i_op_v

end do ! iblk


!-----------------------------------------------------------------
! XII) account for prefactor
!-----------------------------------------------------------------
psi_inout = factor_op_v*psi_inout


GWLS_TIMAB   = 1525
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


ABI_DEALLOCATE(psik)
ABI_DEALLOCATE(psik_alltoall)
ABI_DEALLOCATE(psik_wrk_alltoall)
ABI_DEALLOCATE(psik_tmp_alltoall)
ABI_DEALLOCATE(psik_in_alltoall)

ABI_DEALLOCATE(psir)
ABI_DEALLOCATE(psir_ext)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(total_time2)

total_time = total_time + total_time2 - total_time1


GWLS_TIMAB   = 1524
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

10 format(A)
12 format(A,I5)
14 format(A,F24.16)

end subroutine Pk
!!***

!!****f* m_hamiltonian/epsilon_k
!! NAME
!!  epsilon_k
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_polarisability
!!
!! CHILDREN
!!      epsilon_k
!!
!! SOURCE

subroutine epsilon_k(psi_out,psi_in,omega)

implicit none

real(dp), intent(out) :: psi_out(2,npw_k)
real(dp), intent(in)  :: psi_in(2,npw_k), omega(2)

! *************************************************************************

psi_out = psi_in
call sqrt_vc_k(psi_out)
call Pk(psi_out,omega)
call sqrt_vc_k(psi_out)

psi_out = psi_in - psi_out
end subroutine epsilon_k
!!***

!!****f* m_hamiltonian/set_dielectric_function_frequency
!! NAME
!!  set_dielectric_function_frequency
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_ComputeCorrelationEnergy,gwls_ComputePoles,gwls_GenerateEpsilon
!!
!! CHILDREN
!!      epsilon_k
!!
!! SOURCE

subroutine set_dielectric_function_frequency(omega)
!----------------------------------------------------------------------------------------------------
! This routine sets the value of the module's frequency.
!----------------------------------------------------------------------------------------------------
implicit none
real(dp), intent(in) :: omega(2)

! *************************************************************************

matrix_function_omega(:) = omega(:)

end subroutine set_dielectric_function_frequency
!!***

!!****f* m_hamiltonian/matrix_function_epsilon_k
!! NAME
!!  matrix_function_epsilon_k
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
!!      epsilon_k
!!
!! SOURCE

subroutine matrix_function_epsilon_k(vector_out,vector_in,Hsize)
!----------------------------------------------------------------------------------------------------
! This function is a simple wrapper around epsilon_k to be fed to the Lanczos
! algorithm.
!----------------------------------------------------------------------------------------------------
implicit none
integer,      intent(in)  :: Hsize
complex(dpc), intent(out) :: vector_out(Hsize)
complex(dpc), intent(in)  :: vector_in(Hsize)


! local variables
real(dp)  :: psik (2,npw_k)
real(dp)  :: psik2(2,npw_k)

! *************************************************************************

! convert from one format to the other
psik(1,:) = dble (vector_in(:))
psik(2,:) = dimag(vector_in(:))

! act on vector


call epsilon_k(psik2 ,psik, matrix_function_omega)

! convert back
vector_out = cmplx_1*psik2(1,:)+cmplx_i*psik2(2,:)

end subroutine matrix_function_epsilon_k
!!***

end module m_gwls_polarisability
!!***

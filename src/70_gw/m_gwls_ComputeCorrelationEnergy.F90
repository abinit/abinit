!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gwls_ComputeCorrelationEnergy
!! NAME
!! m_gwls_ComputeCorrelationEnergy
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

module m_gwls_ComputeCorrelationEnergy

! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_polarisability
use m_gwls_model_polarisability
use m_gwls_DielectricArray
use m_gwls_ComputePoles
use m_gwls_Projected_AT
use m_gwls_Projected_BT
use m_gwls_GWlanczos
use m_gwls_GenerateEpsilon
use m_gwls_GWanalyticPart
use m_gwls_TimingLog
use m_gwls_LanczosBasis

! abinit modules
use defs_basis
use defs_wvltypes
use m_abicore
use m_xmpi
use m_pawang
use m_errors
use m_dtset

use m_time,             only : timab
use m_io_tools,         only : get_unit,open_file
use m_paw_dmft,         only : paw_dmft_type
use m_ebands,           only : ebands_init, ebands_free
use m_gaussian_quadrature, only: gaussian_quadrature_gegenbauer, gaussian_quadrature_legendre


implicit none
save
private
!!***

!!***
public :: compute_correlations_shift_lanczos
public :: compute_correlations_no_model_shift_lanczos
!!***
contains

!!****f* m_gwls_ComputeCorrelationEnergy/compute_correlations_shift_lanczos
!! NAME
!!  compute_correlations_shift_lanczos
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_sternheimer
!!
!! CHILDREN
!!
!! SOURCE

subroutine compute_correlations_shift_lanczos(dtset, Sigma_x,Vxc_energy,debug)
!----------------------------------------------------------------------------------------------------
!
!
! This function computes the correlation energy Sigma_c(w) for the state we wish to correct.
! The shift lanczos algorithm is used.
!
!----------------------------------------------------------------------------------------------------

type(dataset_type),intent(in) :: dtset

real(dp),intent(in)  :: Sigma_x,Vxc_energy
logical, intent(in)  :: debug


!Local variables

integer :: n_ext_freq

integer :: npt_gauss
integer :: print_debug


integer :: kmax_poles
integer :: kmax_model
integer :: lmax_model


integer :: kmax_analytic
integer :: kmax_numeric

real(dp) :: omega_static

real(dp) :: lorentzian
integer  :: iw_ext
integer  :: iw

real(dp) , allocatable :: psie_k(:,:)


real(dp),  allocatable :: epsilon_eigenvalues_0(:)
real(dp),  allocatable :: epsilon_model_eigenvalues_0(:)


complex(dpc), allocatable   :: AT_Lanczos(:,:)
complex(dpc), allocatable   :: AT_model_Lanczos(:,:)


! To use lanczos instead of sqmr
complex(dpc), allocatable   :: Lbasis_diagonalize_dielectric_terms(:,:)
complex(dpc), allocatable   :: hermitian_static_eps_m1_minus_eps_model_m1(:,:)
real(dp), allocatable       :: eigenvalues_static_eps_m1_minus_eps_model_m1(:)
real(dp), allocatable       :: eigenvalues_static_eps_model_m1_minus_one(:)
complex(dpc), allocatable   :: work(:)
real(dp), allocatable       :: rwork(:)
integer                     :: lwork
integer                     :: info
integer                     :: l


integer        :: debug_unit
character(50)  :: debug_filename

real(dp)       :: time1, time2, time
real(dp)       :: total_time1,total_time2,total_time
real(dp)       :: setup_time1, setup_time2, setup_time
real(dp)       :: freq_time1, freq_time2, freq_time

integer              :: nfrequencies
real(dp),allocatable :: list_projection_frequencies(:)

complex(dpc), allocatable   :: array_integrand_exact_sector(:,:)
complex(dpc), allocatable   :: array_integrand_model_sector(:,:)


complex(dpc), allocatable :: tmp_dielectric_array(:,:,:)

real(dp)        :: external_omega

character(256)  :: timing_string

integer :: recy_line_size
character(128) :: recy_name
logical :: local_tmp_exist
logical :: use_model

! Energy contributions

real(dp)       :: pole_energy

real(dp)       :: sigma_A_Lanczos
real(dp)       :: sigma_A_model_Lanczos
real(dp)       :: sigma_B_Lanczos
real(dp)       :: sigma_B_model_Lanczos

real(dp)       :: correlations
real(dp)       :: renormalized_energy

real(dp)       :: second_model_parameter

real(dp):: tsec(2)
integer :: GWLS_TIMAB, OPTION_TIMAB
character(500) :: msg

! *************************************************************************

!--------------------------------------------------------------------------------
!
! Set up variables and allocate arrays
!
!--------------------------------------------------------------------------------

call cpu_time(total_time1)
!Variable allocation and initialization
model_number        = dtset%gwls_diel_model
model_parameter     = dtset%gwls_model_parameter
npt_gauss           = dtset%gwls_npt_gauss_quad
print_debug         = dtset%gwls_print_debug

first_seed          = dtset%gwls_first_seed
e                   = dtset%gwls_band_index


!second_model_parameter  = dtset%gwls_second_model_parameter
second_model_parameter  = zero



! set variables from gwls_GenerateEpsilon module
kmax   = dtset%gwls_stern_kmax
nseeds = dtset%gwls_nseeds

kmax_model    = dtset%gwls_kmax_complement
kmax_poles    = dtset%gwls_kmax_poles
kmax_analytic = dtset%gwls_kmax_analytic
kmax_numeric  = dtset%gwls_kmax_numeric

n_ext_freq    = dtset%gw_customnfreqsp

use_model     = .True.

call cpu_time(setup_time1)

!--------------------------------------------------------------------------------
!
! Extract the frequencies at which the integrand will be evaluated
! add the value zero in the set.
!
!--------------------------------------------------------------------------------

call generate_frequencies_and_weights(npt_gauss)

!--------------------------------------------------------------------------------
!
!
! Compute the static bases for the exact and model
! dielectric operator
!
!
!--------------------------------------------------------------------------------

omega_static = zero
! define dimensions
lmax         = nseeds*kmax
lmax_model   = nseeds*kmax_model

! Allocate arrays which will contain basis
call setup_Lanczos_basis(lmax,lmax_model)

! allocate eigenvalues array
ABI_ALLOCATE(epsilon_eigenvalues_0, (lmax))
ABI_ALLOCATE(epsilon_model_eigenvalues_0, (lmax_model))

! set omega=0 for exact dielectric operator
call set_dielectric_function_frequency([0.0_dp,omega_static])

! and make note that the Sternheimer solutions must be kept (for use in the projected Sternheimer section).
if(dtset%gwls_recycle == 1) then
  ABI_ALLOCATE(Sternheimer_solutions_zero,(2,npw_k,lmax,nbandv))
  Sternheimer_solutions_zero = zero
  write_solution = .true.
end if
if(dtset%gwls_recycle == 2) then
  write(recy_name,'(A,I0.4,A)') "Recycling_",mpi_enreg%me,".dat"

  inquire(iolength=recy_line_size) cg(:,1:npw_k)

  inquire(file='local_tmp', exist=local_tmp_exist)
  if(local_tmp_exist) recy_name = 'local_tmp/' // recy_name(1:118)

  if (open_file(file=recy_name,iomsg=msg,newunit=recy_unit,access='direct',form='unformatted',&
&               status='replace',recl=recy_line_size)/=0) then
    MSG_ERROR(msg)
  end if

  write_solution = .true.
end if

call cpu_time(time1)
! Compute the Lanczos basis using Block Lanczos; store
! basis in Lbasis_lanczos
GWLS_TIMAB   = 1504
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call driver_generate_dielectric_matrix( matrix_function_epsilon_k, &
nseeds, kmax, &
epsilon_eigenvalues_0, &
Lbasis_lanczos, debug)
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time = time2-time1

write(timing_string,'(A)')  "Time to compute the EXACT Static Dielectric Matrix  :   "
call write_timing_log(timing_string,time)


call output_epsilon_eigenvalues(lmax,epsilon_eigenvalues_0,1)


! The Sternheimer solutions at $\omega = 0$ have been stored.
write_solution = .false.


! Prepare the model dielectric operator
call setup_Pk_model(omega_static,second_model_parameter)

call cpu_time(time1)
! Compute the Lanczos basis of the model operator using Block Lanczos; store
! basis in Lbasis_model_lanczos
GWLS_TIMAB   = 1505
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call driver_generate_dielectric_matrix(matrix_function_epsilon_model_operator, &
nseeds, kmax_model, &
epsilon_model_eigenvalues_0, &
Lbasis_model_lanczos, debug)
OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time = time2-time1

write(timing_string,'(A)')  "Time to compute the MODEL Static Dielectric Matrix  :   "
call write_timing_log(timing_string,time)


call output_epsilon_eigenvalues(lmax_model,epsilon_model_eigenvalues_0,2)


ABI_ALLOCATE(eigenvalues_static_eps_model_m1_minus_one, (lmax_model))

do l = 1, lmax_model
eigenvalues_static_eps_model_m1_minus_one(l) = one/epsilon_model_eigenvalues_0(l)-one
end do



!--------------------------------------------------------------------------------
!
!
! Prepare and compute the projection of the dielectric Sternheimer equations
!
!
!--------------------------------------------------------------------------------

!  Setup various arrays necessary for the Sternheimer projection scheme
nfrequencies = dtset%gwls_n_proj_freq
ABI_ALLOCATE(list_projection_frequencies,(nfrequencies))

list_projection_frequencies = dtset%gwls_list_proj_freq

call cpu_time(time1)
GWLS_TIMAB   = 1507
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

!call setup_projected_Sternheimer_epsilon(lmax, npt_gauss, second_model_parameter, &
!                                        list_projection_frequencies,nfrequencies,debug)


call ProjectedSternheimerEpsilon(lmax, npt_gauss, second_model_parameter, &
list_projection_frequencies,nfrequencies,&
epsilon_eigenvalues_0,debug,use_model)




! The Sternheimer solutions at $\omega = 0$ have been used to make the basis for the projected Sternheimer equations.
if(dtset%gwls_recycle == 1) then
  ABI_DEALLOCATE(Sternheimer_solutions_zero)
end if
if(dtset%gwls_recycle == 2) then
  close(recy_unit,status='delete')
end if

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to setup and compute the projected Sternheimer epsilon     :   "
call write_timing_log(timing_string,time)


call cpu_time(time1)

GWLS_TIMAB   = 1508
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call compute_eps_m1_minus_eps_model_m1(lmax, npt_gauss)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to compute eps^{-1}-eps_model^{-1}             :   "
call write_timing_log(timing_string,time)


!--------------------------------------------------------------------------------
!
!
! Compute the model dielectric array
!
!
!--------------------------------------------------------------------------------

call cpu_time(time1)

GWLS_TIMAB   = 1509
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call compute_eps_model_m1_minus_one(lmax_model, npt_gauss, second_model_parameter, epsilon_model_eigenvalues_0)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)
ABI_DEALLOCATE(epsilon_model_eigenvalues_0)


call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to compute eps_model^{-1}-1                    :   "
call write_timing_log(timing_string,time)


!--------------------------------------------------------------------------------
!
!
! We no longer need the Lanczos basis in its current form!
! Modify the basis so that it now contains (V^1/2.L)^*.psie
!
!
!--------------------------------------------------------------------------------

call cpu_time(time1)
ABI_ALLOCATE(psie_k, (2,npw_k))

GWLS_TIMAB   = 1510
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


psie_k = cg(:,(e-1)*npw_k+1:e*npw_k)

call modify_Lbasis_Coulomb(psie_k, lmax, lmax_model)

ABI_DEALLOCATE(psie_k)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to modify the Lanczos basis                    :   "
call write_timing_log(timing_string,time)

!--------------------------------------------------------------------------------
!
!
!
! Diagonalize the static array eps^{-1} - eps^{-1}_model in order to
! be able to apply diagonal shift Lanczos.
!
!--------------------------------------------------------------------------------

call cpu_time(time1)
! diagonalize the static eps^{-1} - eps^{-1}_model array, so as the use the diagonal Lanczos procedure
! for the analytical term
GWLS_TIMAB   = 1511
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

ABI_ALLOCATE(Lbasis_diagonalize_dielectric_terms, (npw_k,lmax))
ABI_ALLOCATE(hermitian_static_eps_m1_minus_eps_model_m1, (lmax,lmax))

ABI_ALLOCATE(eigenvalues_static_eps_m1_minus_eps_model_m1, (lmax))

ABI_ALLOCATE(rwork, (3*lmax-2))


hermitian_static_eps_m1_minus_eps_model_m1(:,:) = eps_m1_minus_eps_model_m1(:,:,1)



! WORK QUERRY
lwork = -1
ABI_ALLOCATE(work, (1))
call ZHEEV( 'V',        & ! Compute eigenvectors and eigenvalues
'U',        & ! use Upper triangular part
lmax,        & ! order of matrix
hermitian_static_eps_m1_minus_eps_model_m1,  & ! initial matrix on input; eigenvectors on output
lmax,        & ! LDA
eigenvalues_static_eps_m1_minus_eps_model_m1,& ! eigenvalues
work, lwork, rwork, info )  ! work stuff

if ( info /= 0) then
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info,' in ZHEEV (1), gwls_ComputeCorrelationEnergy'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if

! COMPUTATION
lwork = nint(dble(work(1)))
ABI_DEALLOCATE(work)
ABI_ALLOCATE(work, (lwork))
call ZHEEV( 'V',        & ! Compute eigenvectors and eigenvalues
'U',        & ! use Upper triangular part
lmax,        & ! order of matrix
hermitian_static_eps_m1_minus_eps_model_m1,  & ! initial matrix on input; eigenvectors on output
lmax,        & ! LDA
eigenvalues_static_eps_m1_minus_eps_model_m1,& ! eigenvalues
work, lwork, rwork, info )  ! work stuff

if ( info /= 0) then
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*********************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info,' in ZHEEV (2), gwls_ComputeCorrelationEnergy'
  write(debug_unit,'(A)')      '*********************************************************************************************'

  close(debug_unit)

end if



ABI_DEALLOCATE(work)
ABI_DEALLOCATE(rwork)

!--------------------------------------------------------------------------------
!
! update basis: L' = L . Q
!
!
! CAREFUL!!! We must multiply by conjg(hermitian_static_eps_m1_minus_eps_model_m1),
!            which is the COMPLEX CONJUGATE of the eigenvectors of the matrix
!            eps^{-1}-eps_m^{-1}, because we have MODIFIED the basis Lbasis_lanczos
!            to contain the complex conjugate of the eigenvectors of eps.
!            This is somewhat subtle, but forgetting to do this leads to small errors
!            in the results...
!--------------------------------------------------------------------------------
hermitian_static_eps_m1_minus_eps_model_m1 = conjg(hermitian_static_eps_m1_minus_eps_model_m1)

call ZGEMM('N','N',npw_k,lmax,lmax,cmplx_1,Lbasis_lanczos,npw_k,  &
hermitian_static_eps_m1_minus_eps_model_m1, &
lmax,cmplx_0,Lbasis_diagonalize_dielectric_terms,npw_k)


OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to diagonalize eps^{-1}(0)-eps^{-1}(0)_model   :   "
call write_timing_log(timing_string,time)

ABI_DEALLOCATE(hermitian_static_eps_m1_minus_eps_model_m1)


call cpu_time(setup_time2)
setup_time = setup_time2 - setup_time1

write(timing_string,'(A)')  "       TOTAL DIELECTRIC SETUP TIME                  :   "
call write_timing_log(timing_string,setup_time)

!--------------------------------------------------------------------------------
!
! Compute the Analytic energy using Shift Lanczos
!
!--------------------------------------------------------------------------------

call cpu_time(time1)

GWLS_TIMAB   = 1512
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


ABI_ALLOCATE(AT_Lanczos,(n_ext_freq,lmax))
call compute_AT_shift_Lanczos(n_ext_freq,dtset%gw_freqsp, model_parameter, lmax, Lbasis_diagonalize_dielectric_terms,&
&                             kmax_analytic, AT_Lanczos)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


ABI_DEALLOCATE(Lbasis_diagonalize_dielectric_terms)

call cpu_time(time2)
time = time2 - time1

write(timing_string,'(A)')  "Time to compute analytical term by SHIFT LANCZOS    :   "
call write_timing_log(timing_string,time)


call cpu_time(time1)

GWLS_TIMAB   = 1513
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

ABI_ALLOCATE(AT_model_Lanczos,(n_ext_freq,lmax_model))
call compute_AT_shift_Lanczos(n_ext_freq,dtset%gw_freqsp, model_parameter, lmax_model,  &
Lbasis_model_lanczos, kmax_analytic, AT_model_Lanczos)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)
time = time2 - time1

write(timing_string,'(A)')  "Time to compute analytical MODEL by SHIFT LANCZOS   :  "
call write_timing_log(timing_string,time)


!--------------------------------------------------------------------------------
!
! Compute the Numeric energy using Shift Lanczos
!
!--------------------------------------------------------------------------------

call cpu_time(time1)

ABI_ALLOCATE(array_integrand_exact_sector,(npt_gauss+1,n_ext_freq))

ABI_ALLOCATE( tmp_dielectric_array, (lmax,lmax,npt_gauss+1))

do iw = 1, npt_gauss + 1

lorentzian      = model_parameter**2/(list_omega(iw)**2+model_parameter**2)
tmp_dielectric_array(:,:,iw) = eps_m1_minus_eps_model_m1(:,:,iw)-lorentzian*eps_m1_minus_eps_model_m1(:,:,1)

end do
GWLS_TIMAB   = 1514
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call compute_projected_BT_shift_Lanczos(n_ext_freq, dtset%gw_freqsp, lmax, Lbasis_lanczos,         &
kmax_numeric, npt_gauss, tmp_dielectric_array, array_integrand_exact_sector )

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

ABI_DEALLOCATE( tmp_dielectric_array)
call cpu_time(time2)
time = time2 - time1
write(timing_string,'(A)')  "Time to compute numerical term by SHIFT LANCZOS     :   "
call write_timing_log(timing_string,time)



call cpu_time(time1)

ABI_ALLOCATE(array_integrand_model_sector,(npt_gauss+1,n_ext_freq))

!ABI_ALLOCATE( tmp_dielectric_array, (lmax_model,lmax_model,npt_gauss+1))
ABI_ALLOCATE( tmp_dielectric_array, (lmax_model,blocksize_epsilon,npt_gauss+1))


do iw = 1, npt_gauss + 1

lorentzian      = model_parameter**2/(list_omega(iw)**2+model_parameter**2)
!tmp_dielectric_array(:,:,iw) = eps_model_m1_minus_one(:,:,iw)-lorentzian*eps_model_m1_minus_one(:,:,1)
tmp_dielectric_array(:,:,iw) = eps_model_m1_minus_one_DISTR(:,:,iw)-lorentzian*eps_model_m1_minus_one_DISTR(:,:,1)

end do

GWLS_TIMAB   = 1515
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

!call compute_projected_BT_shift_Lanczos(n_ext_freq , dtset%gw_freqsp, lmax_model, Lbasis_model_lanczos,         &
!                                        kmax_numeric, npt_gauss, tmp_dielectric_array, array_integrand_model_sector )


call compute_projected_BT_shift_Lanczos_DISTRIBUTED(n_ext_freq, dtset%gw_freqsp, lmax_model, blocksize_epsilon, &
model_lanczos_vector_belongs_to_this_node, model_lanczos_vector_index,  &
Lbasis_model_lanczos, kmax_numeric, npt_gauss, tmp_dielectric_array, &
array_integrand_model_sector )

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

ABI_DEALLOCATE(tmp_dielectric_array)
call cpu_time(time2)
time = time2 - time1
write(timing_string,'(A)')  "Time to compute numerical model   SHIFT LANCZOS     :   "
call write_timing_log(timing_string,time)


!--------------------------------------------------------------------------------
!
! set up arrays for poles
!
!--------------------------------------------------------------------------------



!call generate_degeneracy_table_for_poles(debug) ! so we can compute Poles contributions
call generate_degeneracy_table_for_poles(.true.) ! so we can compute Poles contributions

!--------------------------------------------------------------------------------
!
! Print contributions to sigma_A_Lanczos  to a file
!
!--------------------------------------------------------------------------------
call output_Sigma_A_by_eigenvalues(n_ext_freq,lmax,dtset%gw_freqsp,AT_Lanczos,eigenvalues_static_eps_m1_minus_eps_model_m1,2)
call output_Sigma_A_by_eigenvalues(n_ext_freq,lmax_model,dtset%gw_freqsp,AT_model_Lanczos,&
&                                  eigenvalues_static_eps_model_m1_minus_one,3)

!epsilon_eigenvalues_0

ABI_DEALLOCATE(epsilon_eigenvalues_0)
!--------------------------------------------------------------------------------
!
! Iterate on external frequencies
!
!--------------------------------------------------------------------------------


do iw_ext = 1, dtset%gw_customnfreqsp

call cpu_time(freq_time1)
external_omega = dtset%gw_freqsp(iw_ext)

write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A,I4,A,F8.4,A)')  "#  Frequency # ",iw_ext," omega = ",external_omega," Ha"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)


!--------------------------------------------------------------------------------
!
! compute the pole term
!                        CAREFUL! The real valence states must still be allocated
!                        for the dielectric operator to work properly
!
!--------------------------------------------------------------------------------

call cpu_time(time1)
GWLS_TIMAB   = 1516
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

pole_energy    = compute_Poles(external_omega,kmax_poles,debug)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call cpu_time(time2)

time = time2-time1
write(timing_string,'(A)')  "Time to compute the Poles contribution              :   "
call write_timing_log(timing_string,time)



!================================================================================
! Compute the contributions from the analytic term
!================================================================================
GWLS_TIMAB   = 1517
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)


call cpu_time(time1)

sigma_A_Lanczos      = dble(sum(AT_Lanczos(iw_ext,:)*eigenvalues_static_eps_m1_minus_eps_model_m1(:)))

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time  Tr[(eps^{-1}-eps_model^{-1}).AT] AFTER SHIFT  :   "
call write_timing_log(timing_string,time)


call cpu_time(time1)

sigma_A_model_Lanczos= dble(sum(AT_model_Lanczos(iw_ext,:)*eigenvalues_static_eps_model_m1_minus_one(:)))

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time for Tr[ (eps_model^{-1}-1) . AT ] AFTER SHIFT  :   "
call write_timing_log(timing_string,time)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)
!--------------------------------------------------------------------------------
!
! compute integrand
!
!--------------------------------------------------------------------------------
GWLS_TIMAB   = 1518
OPTION_TIMAB = 1
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

call compute_integrands_shift_lanczos(iw_ext, n_ext_freq, npt_gauss, array_integrand_exact_sector, &
array_integrand_model_sector, sigma_B_Lanczos, sigma_B_model_Lanczos)

OPTION_TIMAB = 2
call timab(GWLS_TIMAB,OPTION_TIMAB,tsec)

!--------------------------------------------------------------------------------
!
! Output results
!
!--------------------------------------------------------------------------------


call output_results(iw_ext,npt_gauss, lmax,lmax_model, model_parameter, second_model_parameter,  &
external_omega, Sigma_x,Vxc_energy,pole_energy,                          &
sigma_A_Lanczos,sigma_A_model_Lanczos,sigma_B_Lanczos,sigma_B_model_Lanczos)


call cpu_time(freq_time2)
freq_time = freq_time2-freq_time1

write(timing_string,'(A)')  "               TOTAL FREQUENCY TIME                 :   "
call write_timing_log(timing_string,freq_time)



correlations        = pole_energy+sigma_A_Lanczos+sigma_A_model_Lanczos+sigma_B_Lanczos+sigma_B_model_Lanczos

renormalized_energy = eig(e) + Sigma_x-Vxc_energy +correlations
write(std_out,10) '                               '
write(std_out,14) ' For omega                   : ',external_omega     ,' Ha = ',external_omega     *Ha_eV,' eV'
write(std_out,14) '   <psi_e | Sigma_c  | psi_e>: ',correlations       ,' Ha = ',correlations       *Ha_eV,' eV'
write(std_out,14) '   eps_e + <Sigma_xc - V_xc> : ',renormalized_energy,' Ha = ',renormalized_energy*Ha_eV,' eV'

write(ab_out,10) '                               '
write(ab_out,14) ' For omega                   : ',external_omega     ,' Ha = ',external_omega     *Ha_eV,' eV'
write(ab_out,14) '   <psi_e | Sigma_c  | psi_e>: ',correlations       ,' Ha = ',correlations       *Ha_eV,' eV'
write(ab_out,14) '   eps_e + <Sigma_xc - V_xc> : ',renormalized_energy,' Ha = ',renormalized_energy*Ha_eV,' eV'


end do

ABI_DEALLOCATE(AT_Lanczos)
ABI_DEALLOCATE(AT_model_Lanczos)
ABI_DEALLOCATE(array_integrand_exact_sector)
ABI_DEALLOCATE(array_integrand_model_sector)
ABI_DEALLOCATE(eigenvalues_static_eps_model_m1_minus_one)
ABI_DEALLOCATE(eigenvalues_static_eps_m1_minus_eps_model_m1)
call clean_degeneracy_table_for_poles()
call cleanup_Pk_model()
call cleanup_Lanczos_basis()
call cleanup_projected_Sternheimer_epsilon()

call cpu_time(total_time2)
total_time = total_time2-total_time1
write(timing_string,'(A)')  "               TOTAL TIME                           :   "
call write_timing_log(timing_string,total_time)


10 format(A)
14 format(A,ES24.16,A,F16.8,A)

end subroutine compute_correlations_shift_lanczos
!!***

!!****f* m_gwls_ComputeCorrelationEnergy/compute_correlations_no_model_shift_lanczos
!! NAME
!!  compute_correlations_no_model_shift_lanczos
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_sternheimer
!!
!! CHILDREN
!!
!! SOURCE

subroutine compute_correlations_no_model_shift_lanczos(dtset, Sigma_x,Vxc_energy,debug)
!----------------------------------------------------------------------------------------------------
!
!
! This function computes the correlation energy Sigma_c(w) for the state we wish to correct.
!
! this subroutine does not rely on the use of a model dielectric operator. Thus
!
!        Sigma^A(w) = Tr[ (eps^{-1}(0)-1). A^T(w)]
!        Sigma^N(w) = int dw' Tr[{(eps^{-1}(w')-1)-f(w')(eps^{-1}(0)-1)}B^T(w';w)]
!
! Shift lanczos is used for the resolvents.
!----------------------------------------------------------------------------------------------------

type(dataset_type),intent(in) :: dtset

real(dp),intent(in)  :: Sigma_x, Vxc_energy
logical, intent(in)  :: debug


!Local variables

real(dp):: Sigma_x_Lanczos_projected
integer :: npt_gauss
integer :: print_debug


integer :: kmax_poles

integer :: lmax_model

integer :: kmax_analytic
integer :: kmax_numeric
integer :: n_ext_freq


real(dp) :: omega_static

real(dp) :: lorentzian
integer  :: iw_ext
integer  :: iw



real(dp) , allocatable :: psie_k(:,:)

real(dp),  allocatable :: epsilon_eigenvalues_0(:)


complex(dpc), allocatable   :: AT_Lanczos(:,:)



real(dp)       :: time1, time2, time
real(dp)       :: total_time1,total_time2,total_time
real(dp)       :: setup_time1, setup_time2, setup_time
real(dp)       :: freq_time1, freq_time2, freq_time

integer              :: nfrequencies
real(dp),allocatable :: list_projection_frequencies(:)



complex(dpc), allocatable   :: array_integrand_exact_sector(:,:)
complex(dpc), allocatable   :: array_integrand_model_sector(:,:)
complex(dpc), allocatable   :: tmp_dielectric_array(:,:,:)

real(dp)        :: external_omega

character(256)  :: timing_string

integer :: recy_line_size
character(128) :: recy_name
logical :: local_tmp_exist
character(500) :: msg

! Energy contributions

real(dp)       :: pole_energy

real(dp)       :: sigma_A_Lanczos
real(dp)       :: sigma_A_model_Lanczos
real(dp)       :: sigma_B_Lanczos
real(dp)       :: sigma_B_model_Lanczos

real(dp)       :: correlations
real(dp)       :: renormalized_energy

logical :: use_model

! *************************************************************************


!--------------------------------------------------------------------------------
!
! Set up variables and allocate arrays
!
!--------------------------------------------------------------------------------

!Variable allocation and initialization
model_number        = dtset%gwls_diel_model
model_parameter     = dtset%gwls_model_parameter
npt_gauss           = dtset%gwls_npt_gauss_quad
print_debug         = dtset%gwls_print_debug

first_seed          = dtset%gwls_first_seed
e                   = dtset%gwls_band_index


! set variables from gwls_GenerateEpsilon module
kmax   = dtset%gwls_stern_kmax
nseeds = dtset%gwls_nseeds

kmax_poles  = dtset%gwls_kmax_poles

kmax_poles    = dtset%gwls_kmax_poles
kmax_analytic = dtset%gwls_kmax_analytic
kmax_numeric  = dtset%gwls_kmax_numeric

n_ext_freq    = dtset%gw_customnfreqsp


use_model     = .False.


call cpu_time(setup_time1)
!--------------------------------------------------------------------------------
!
! Extract the frequencies at which the integrand will be evaluated
! add the value zero in the set.
!
!--------------------------------------------------------------------------------

call generate_frequencies_and_weights(npt_gauss)

!--------------------------------------------------------------------------------
!
!
! Compute the static bases for the exact dielectric operator
!
!
!--------------------------------------------------------------------------------


omega_static = zero
! define dimensions
lmax         = nseeds*kmax

! Allocate arrays which will contain basis
! the 0 indicates we will not use arrays for the model dielectric operator
call setup_Lanczos_basis(lmax,0)

! allocate eigenvalues array
ABI_ALLOCATE(epsilon_eigenvalues_0, (lmax))

! set omega=0 for exact dielectric operator
call set_dielectric_function_frequency([0.0_dp,omega_static])

! and make note that the Sternheimer solutions must be kept (for use in the projected Sternheimer section).
if(dtset%gwls_recycle == 1) then
  ABI_ALLOCATE(Sternheimer_solutions_zero,(2,npw_k,lmax,nbandv))
  Sternheimer_solutions_zero = zero
  write_solution = .true.
end if
if(dtset%gwls_recycle == 2) then
  write(recy_name,'(A,I0.4,A)') "Recycling_",mpi_enreg%me,".dat"

  inquire(iolength=recy_line_size) cg(:,1:npw_k)

  inquire(file='local_tmp', exist=local_tmp_exist)
  if(local_tmp_exist) recy_name = 'local_tmp/' // recy_name(1:118)

  if (open_file(file=recy_name,iomsg=msg,newunit=recy_unit,access='direct',form='unformatted',&
&               status='replace',recl=recy_line_size)/=0) then
    MSG_ERROR(msg)
  end if

  write_solution = .true.
end if


call cpu_time(time1)
! Compute the Lanczos basis using Block Lanczos; store
! basis in Lbasis_lanczos
call driver_generate_dielectric_matrix( matrix_function_epsilon_k, &
nseeds, kmax, &
epsilon_eigenvalues_0, &
Lbasis_lanczos, debug)
call cpu_time(time2)
time = time2-time1

write(timing_string,'(A)')  "Time to compute the EXACT Static Dielectric Matrix  :   "
call write_timing_log(timing_string,time)

! The Sternheimer solutions at $\omega = 0$ have been stored.
write_solution = .false.


call output_epsilon_eigenvalues(lmax,epsilon_eigenvalues_0,1)


! compute the Exchange energy, when it is projected on the Lanczos basis
! CAREFUL! This must be done BEFORE we modify the lanczos basis
Sigma_x_Lanczos_projected =  exchange(e, Lbasis_lanczos)



!--------------------------------------------------------------------------------
!
!
! Prepare and compute the projection of the dielectric Sternheimer equations
!
!
!--------------------------------------------------------------------------------

!  Setup various arrays necessary for the Sternheimer projection scheme
nfrequencies = dtset%gwls_n_proj_freq
ABI_ALLOCATE(list_projection_frequencies,(nfrequencies))

list_projection_frequencies = dtset%gwls_list_proj_freq

call cpu_time(time1)
! The explicit "false" as the last argument is for the optional
! variable "use_model"; we are not using a model here!

!call setup_projected_Sternheimer_epsilon(lmax, npt_gauss, zero, &
!                list_projection_frequencies,nfrequencies,debug,.false.)


call ProjectedSternheimerEpsilon(lmax, npt_gauss, zero, &
list_projection_frequencies,nfrequencies,&
epsilon_eigenvalues_0,debug,use_model)



!call cpu_time(time2)
!time = time2-time1
!write(timing_string,'(A)')  "Time to setup the projected Sternheimer epsilon     :   "
!call write_timing_log(timing_string,time)


! The Sternheimer solutions at $\omega = 0$ have been used to make the basis for the projected Sternheimer equations.
if(dtset%gwls_recycle == 1) then
  ABI_DEALLOCATE(Sternheimer_solutions_zero)
end if
if(dtset%gwls_recycle == 2) then
  close(recy_unit,status='delete')
end if


!call cpu_time(time1)
!call compute_projected_Sternheimer_epsilon(lmax, npt_gauss, epsilon_eigenvalues_0,debug)



call cpu_time(time2)

time = time2-time1
write(timing_string,'(A)')  "Time to compute the projected Sternheimer epsilon   :   "
call write_timing_log(timing_string,time)



call cpu_time(time1)
call compute_eps_m1_minus_one(lmax, npt_gauss)
call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to compute eps^{-1}-I                          :   "
call write_timing_log(timing_string,time)


!--------------------------------------------------------------------------------
!
!
! We no longer need the Lanczos basis in its current form!
! Modify the basis so that it now contains (V^{1/2}.l)
!
!
!--------------------------------------------------------------------------------

call cpu_time(time1)
ABI_ALLOCATE(psie_k, (2,npw_k))
psie_k = cg(:,(e-1)*npw_k+1:e*npw_k)

lmax_model = 0
call modify_Lbasis_Coulomb(psie_k, lmax, lmax_model) ! lmax_model is set to zero, such that
! the model lanczos basis (which doesn't exist
! in this case) will not be modified

ABI_DEALLOCATE(psie_k)

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to modify the Lanczos basis                    :   "
call write_timing_log(timing_string,time)


call cpu_time(setup_time2)
setup_time = setup_time2 - setup_time1

write(timing_string,'(A)')  "       TOTAL DIELECTRIC SETUP TIME                  :   "
call write_timing_log(timing_string,setup_time)

!--------------------------------------------------------------------------------
!
! Compute the Analytic energy using Shift Lanczos
!
!--------------------------------------------------------------------------------

call cpu_time(time1)

ABI_ALLOCATE(AT_Lanczos,(n_ext_freq,lmax))
! Note that the array eps^{-1}(0) - 1 is diagonal in the lanczos basis already! No need to diagonalize, so we can
! use Lbasis_lanczos directly...
call compute_AT_shift_Lanczos(n_ext_freq,dtset%gw_freqsp, model_parameter, lmax, Lbasis_lanczos, kmax_analytic, AT_Lanczos)

call cpu_time(time2)
time = time2 - time1

write(timing_string,'(A)')  "Time to compute analytical term by SHIFT LANCZOS    :   "
call write_timing_log(timing_string,time)


!--------------------------------------------------------------------------------
!
! Compute the Numeric energy using Shift Lanczos
!
!--------------------------------------------------------------------------------

call cpu_time(time1)

ABI_ALLOCATE(array_integrand_exact_sector,(npt_gauss+1,n_ext_freq))
ABI_ALLOCATE(array_integrand_model_sector,(npt_gauss+1,n_ext_freq))


ABI_ALLOCATE( tmp_dielectric_array, (lmax,lmax,npt_gauss+1))

do iw = 1, npt_gauss + 1

lorentzian      = model_parameter**2/(list_omega(iw)**2+model_parameter**2)
tmp_dielectric_array(:,:,iw) = eps_m1_minus_eps_model_m1(:,:,iw)-lorentzian*eps_m1_minus_eps_model_m1(:,:,1)

end do

call compute_projected_BT_shift_Lanczos(n_ext_freq, dtset%gw_freqsp, lmax, Lbasis_lanczos,         &
kmax_numeric, npt_gauss, tmp_dielectric_array, array_integrand_exact_sector )

array_integrand_model_sector = zero ! just a dummy array in this case

ABI_DEALLOCATE( tmp_dielectric_array)
call cpu_time(time2)
time = time2 - time1
write(timing_string,'(A)')  "Time to compute numerical term by SHIFT LANCZOS     :   "
call write_timing_log(timing_string,time)



!--------------------------------------------------------------------------------
!
! set up arrays for poles
!
!--------------------------------------------------------------------------------

!call generate_degeneracy_table_for_poles(debug) ! so we can compute Poles contributions
call generate_degeneracy_table_for_poles(.true.) ! so we can compute Poles contributions

!--------------------------------------------------------------------------------
!
! Print contributions to sigma_A_Lanczos  to a file
!
!--------------------------------------------------------------------------------
call output_Sigma_A_by_eigenvalues(n_ext_freq,lmax,dtset%gw_freqsp,AT_Lanczos,one/epsilon_eigenvalues_0-one,1)


!--------------------------------------------------------------------------------
!
! Iterate on external frequencies
!
!--------------------------------------------------------------------------------


do iw_ext = 1, dtset%gw_customnfreqsp

call cpu_time(freq_time1)
external_omega = dtset%gw_freqsp(iw_ext)

write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A,I4,A,F8.4,A)')  "#  Frequency # ",iw_ext," omega = ",external_omega," Ha"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)
write(timing_string,'(A)')  "#"
call write_text_block_in_Timing_log(timing_string)


!--------------------------------------------------------------------------------
!
! compute the pole term
!                        CAREFUL! The real valence states must still be allocated
!                        for the dielectric operator to work properly
!
!--------------------------------------------------------------------------------

call cpu_time(time1)

pole_energy    = compute_Poles(external_omega,kmax_poles,debug)

call cpu_time(time2)

time = time2-time1
write(timing_string,'(A)')  "Time to compute the Poles contribution              :   "
call write_timing_log(timing_string,time)

!================================================================================
! Compute the contributions from the analytic term
!================================================================================

call cpu_time(time1)


sigma_A_Lanczos      = dble(sum(AT_Lanczos(iw_ext,:)*(one/epsilon_eigenvalues_0(:)-one)))

ABI_DEALLOCATE(epsilon_eigenvalues_0)

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time for Tr[ (eps_model^{-1}-1) . AT ] AFTER SHIFT  :   "
call write_timing_log(timing_string,time)

!--------------------------------------------------------------------------------
!
! compute integrand
!
!--------------------------------------------------------------------------------


call compute_integrands_shift_lanczos(iw_ext, n_ext_freq, npt_gauss, array_integrand_exact_sector, &
array_integrand_model_sector, sigma_B_Lanczos, sigma_B_model_Lanczos)



!--------------------------------------------------------------------------------
!
! Output results
!
!--------------------------------------------------------------------------------

! just dummy variables so we can use the output_results routine
sigma_A_model_Lanczos = zero
sigma_B_model_Lanczos = zero
call output_results(iw_ext,npt_gauss, lmax, lmax_model, model_parameter, zero,  &
external_omega, Sigma_x,Vxc_energy,pole_energy,        &
sigma_A_Lanczos,sigma_A_model_Lanczos,sigma_B_Lanczos,sigma_B_model_Lanczos, &
Sigma_x_Lanczos_projected )


call cpu_time(freq_time2)
freq_time = freq_time2-freq_time1

write(timing_string,'(A)')  "               TOTAL FREQUENCY TIME                 :   "
call write_timing_log(timing_string,freq_time)



correlations        = pole_energy+sigma_A_Lanczos+sigma_B_Lanczos

renormalized_energy = eig(e) + Sigma_x-Vxc_energy +correlations
write(std_out,10) '                               '
write(std_out,14) ' For omega                   : ',external_omega     ,' Ha = ',external_omega     *Ha_eV,' eV'
write(std_out,14) '  <psi_e | Sigma_c  | psi_e> : ',correlations       ,' Ha = ',correlations       *Ha_eV,' eV'
write(std_out,14) '  eps_e + <Sigma_xc - V_xc>  : ',renormalized_energy,' Ha = ',renormalized_energy*Ha_eV,' eV'

write(ab_out,10) '                               '
write(ab_out,14) ' For omega                   : ',external_omega     ,' Ha = ',external_omega     *Ha_eV,' eV'
write(ab_out,14) '   <psi_e | Sigma_c  | psi_e>: ',correlations       ,' Ha = ',correlations       *Ha_eV,' eV'
write(ab_out,14) '   eps_e + <Sigma_xc - V_xc> : ',renormalized_energy,' Ha = ',renormalized_energy*Ha_eV,' eV'



end do

ABI_DEALLOCATE(AT_Lanczos)
ABI_DEALLOCATE(array_integrand_exact_sector)
ABI_DEALLOCATE(array_integrand_model_sector)
call clean_degeneracy_table_for_poles()
call cleanup_Pk_model()
call cleanup_Lanczos_basis()
call cleanup_projected_Sternheimer_epsilon()

call cpu_time(total_time2)
total_time = total_time2-total_time1
write(timing_string,'(A)')  "               TOTAL TIME                           :   "
call write_timing_log(timing_string,total_time)



10 format(A)
14 format(A,ES24.16,A,F16.8,A)

end subroutine compute_correlations_no_model_shift_lanczos
!!***

!!****f* m_gwls_ComputeCorrelationEnergy/compute_integrands_shift_lanczos
!! NAME
!!  compute_integrands_shift_lanczos
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
!!
!! SOURCE

subroutine compute_integrands_shift_lanczos(iw_ext,n_ext_freq,npt_gauss, array_integrand_exact_sector, &
array_integrand_model_sector, sigma_B_Lanczos, sigma_B_model_Lanczos)
!----------------------------------------------------------------------------------------------------
!
! This subroutine computes the integrands, assuming data was generated by shift lanczos.
!----------------------------------------------------------------------------------------------------

integer,  intent(in)   :: iw_ext, npt_gauss, n_ext_freq

complex(dpc),  intent(in)   :: array_integrand_exact_sector(npt_gauss+1,n_ext_freq)
complex(dpc),  intent(in)   :: array_integrand_model_sector(npt_gauss+1,n_ext_freq)

real(dp), intent(out)  :: sigma_B_Lanczos, sigma_B_model_Lanczos

real(dp) :: integrand_Lanczos , integrand_model_Lanczos
real(dp) :: omega_prime

integer         :: iw
integer         :: io_unit
character(256)  :: title_string
character(256)  :: timing_string


real(dp)       :: time1, time2, time

! *************************************************************************


if (mpi_enreg%me == 0) then
  io_unit = get_unit()

  if (iw_ext < 10) then
    write(title_string,'(A,I1,A)')  'APPROXIMATE_INTEGRANDS_',iw_ext,'.dat'
  else if (iw_ext < 100) then
    write(title_string,'(A,I2,A)')  'APPROXIMATE_INTEGRANDS_',iw_ext,'.dat'
  else
    write(title_string,'(A,I3,A)')  'APPROXIMATE_INTEGRANDS_',iw_ext,'.dat'
  end if

  open(file=title_string,status=files_status_new,unit=io_unit)
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'
  write(io_unit,10) '#'
  write(io_unit,10) '#               Approximate Integrands as a function of frequency'
  write(io_unit,10) '#'
  write(io_unit,10) '#        I1 = Tr[ (eps^{-1}(iw) - eps_model^{-1}(iw)) -                      '
  write(io_unit,10) '#                          f(w)(eps^{-1}(0) - eps_model^{-1}(0)) BT(w) ]'
  write(io_unit,10) '#'
  write(io_unit,10) '#        I2 = Tr[ (eps_model^{-1}(iw) - 1) - f(w)(eps_model^{-1}(0)-1) BT(w)]'
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#         DIAG[I] will represent the contribution coming from taking only the diagonal elements of'
  write(io_unit,10) '#         the arrays in the trace.'
  write(io_unit,10) '#'
  write(io_unit,10) '#      omega (Ha)                  I1                      I2                 gaussian weight'
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'
  flush(io_unit)
end if

call cpu_time(time1)

sigma_B_Lanczos              = zero
sigma_B_model_Lanczos        = zero

do iw = 1,npt_gauss+1

omega_prime     = list_omega(iw)

integrand_Lanczos         =  dble(array_integrand_exact_sector(iw,iw_ext))
integrand_model_Lanczos   =  dble(array_integrand_model_sector(iw,iw_ext))

sigma_B_Lanczos       = sigma_B_Lanczos       + integrand_Lanczos*list_weights(iw)
sigma_B_model_Lanczos = sigma_B_model_Lanczos + integrand_model_Lanczos*list_weights(iw)



if (mpi_enreg%me == 0) write(io_unit,8) omega_prime, integrand_Lanczos , integrand_model_Lanczos, list_weights(iw)



end do
call cpu_time(time2)

if (mpi_enreg%me == 0) then
  write(io_unit,10) ''
  write(io_unit,14) '# Value of the I1 integral: ',sigma_B_Lanczos ,' Ha'
  write(io_unit,14) '# Value of the I2 integral: ',sigma_B_model_Lanczos ,' Ha'
  write(io_unit,10) ''
  write(io_unit,10) ''
  close(io_unit)
end if

time = time2-time1
write(timing_string,'(A)')  "Time compute the integrands and integrals           :   "
call write_timing_log(timing_string,time)

8  format(4ES24.16)
10 format(A)
14 format(A,ES24.16,A)


end subroutine compute_integrands_shift_lanczos
!!***

!!****f* m_gwls_ComputeCorrelationEnergy/output_results
!! NAME
!!  output_results
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
!!
!! SOURCE

subroutine output_results(iw_ext,npt_gauss, lmax,lmax_model, model_parameter, second_model_parameter,  external_omega, &
Sigma_x,Vxc_energy,pole_energy,sigma_A_Lanczos,sigma_A_model_Lanczos,sigma_B_Lanczos,sigma_B_model_Lanczos,&
Sigma_x_Lanczos_projected )
!----------------------------------------------------------------------------------------------------
!
! This subroutine computes the integrands
!----------------------------------------------------------------------------------------------------

integer,  intent(in)   :: iw_ext,lmax,lmax_model,npt_gauss
real(dp), intent(in)   :: model_parameter, second_model_parameter,  external_omega
real(dp), intent(in)   :: Sigma_x,Vxc_energy,pole_energy,sigma_A_Lanczos
real(dp), intent(in)   :: sigma_A_model_Lanczos,sigma_B_Lanczos,sigma_B_model_Lanczos

real(dp), optional, intent(in)   :: Sigma_x_Lanczos_projected

integer         :: io_unit

character(128) :: filename

real(dp)       :: Sigma_c

! *************************************************************************


if (mpi_enreg%me == 0) then
  io_unit = get_unit()
  if (iw_ext < 10) then
    write(filename,'(A,I1,A)')  'ALL_ENERGY_',iw_ext,'.dat'
  else if (iw_ext < 100) then
    write(filename,'(A,I2,A)')  'ALL_ENERGY_',iw_ext,'.dat'
  else
    write(filename,'(A,I3,A)')  'ALL_ENERGY_',iw_ext,'.dat'
  end if


  open(file=filename,status=files_status_new,unit=io_unit)
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'
  write(io_unit,10) '#'
  write(io_unit,10) '#               This file contains the results of the Correlation energy calculation.'
  write(io_unit,10) '# '
  write(io_unit,10) '# Definitions:'
  write(io_unit,10) '# '
  write(io_unit,10) '#                eps_e     = Bare DFT energy of the state                                           '
  write(io_unit,10) '# '
  write(io_unit,10) '#                Sigma_A_1 = Tr[(eps^{-1}(0) - eps_model^{-1}(0)) AT(W) ]                           '
  write(io_unit,10) '# '
  write(io_unit,10) '#                Sigma_A_2 = Tr[(eps_model^{-1}(0) - 1) AT(W) ]                           '
  write(io_unit,10) '# '
  write(io_unit,10) '#                Sigma_B_1 = Int dw I1(w)                                                '
  write(io_unit,10) '# '
  write(io_unit,10) '#                Sigma_B_2 = Int dw I2(w)                                                '
  write(io_unit,10) '# '
  write(io_unit,10) '#        I1(w) = Tr[ (eps^{-1}(iw) - eps_model^{-1}(iw)) -                      '
  write(io_unit,10) '#                          f(w)(eps^{-1}(0) - eps_model^{-1}(0)) BT(w,W) ]'
  write(io_unit,10) '#'
  write(io_unit,10) '#        I2(w) = Tr[ (eps_model^{-1}(iw) - 1) - f(w)(eps_model^{-1}(0)-1) BT(w,W)]'
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#  Parameters:                                                               '
  write(io_unit,10) '#                                                                            '
  write(io_unit,25) '#                 lmax          = ', lmax
  write(io_unit,25) '#                 lmax_model    = ', lmax_model
  write(io_unit,25) '#                 npt_gauss     = ', npt_gauss
  write(io_unit,12) '#                 omega0        = ', model_parameter,'  Ha'
  write(io_unit,12) '#                 epsilon0      = ', second_model_parameter,'  Ha'
  write(io_unit,14) '#                 omega_ext (W) = ', external_omega,'  Ha'
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#   NOTE: if lmax_model = 0, then eps_model = I, the identity.               '
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'
  write(io_unit,10) '                                                                             '
  write(io_unit,30) '       eps_e    (Ha)      :     ', eig(e)
  write(io_unit,10) '                                                                             '
  write(io_unit,30) '      Sigma_x   (Ha)      :     ', Sigma_x

  if (present(Sigma_x_Lanczos_projected) ) then
    write(io_unit,30) '   Sigma_x_PROJECTED (Ha) :     ', Sigma_x_Lanczos_projected
  end if



  write(io_unit,30) '    < V_xc >_e  (Ha)      :     ', Vxc_energy
  write(io_unit,10) '                                                                             '
  write(io_unit,30) '       poles    (Ha)      :     ', pole_energy
  write(io_unit,30) '     Sigma_A_1  (Ha)      :     ', sigma_A_Lanczos
  write(io_unit,30) '     Sigma_A_2  (Ha)      :     ', sigma_A_model_Lanczos
  write(io_unit,30) '     Sigma_B_1  (Ha)      :     ', sigma_B_Lanczos
  write(io_unit,30) '     Sigma_B_2  (Ha)      :     ', sigma_B_model_Lanczos
  write(io_unit,10) '                                                                             '

  Sigma_c = pole_energy+sigma_A_Lanczos+sigma_A_model_Lanczos+sigma_B_Lanczos+sigma_B_model_Lanczos

  write(io_unit,30) '      Sigma_c   (Ha)      :     ',Sigma_c
  write(io_unit,10) '                                                                             '
  write(io_unit,30) '       E_e      (Ha)      :     ', eig(e)+Sigma_x-Vxc_energy+Sigma_c


  close(io_unit)
end if



10 format(A)
12 format(A,ES10.3,A)
14 format(A,ES24.16,A)
25 format(A,I5)
30 format(A,ES24.16)

end subroutine output_results
!!***

!!****f* m_gwls_ComputeCorrelationEnergy/output_epsilon_eigenvalues
!! NAME
!!  output_epsilon_eigenvalues
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
!!
!! SOURCE

subroutine output_epsilon_eigenvalues(lmax,eigenvalues,which_case)
!----------------------------------------------------------------------------------------------------
!       This routine outputs the eigenvalues of the static dielectric matrix
!
!       There are two cases to consider:
!                               1 )   the exact dielectric matrix
!                               2 )   the model dielectric matrix
!----------------------------------------------------------------------------------------------------

integer,      intent(in) ::  lmax, which_case
real(dp),     intent(in) :: eigenvalues(lmax)


integer         :: io_unit
character(128) :: filename

integer         :: l

! *************************************************************************

if (mpi_enreg%me == 0) then
  io_unit = get_unit()

  if (which_case == 1) then
    write(filename,'(A)') "EPSILON_EIGENVALUES.dat"
  else if (which_case == 2) then
    write(filename,'(A)') "MODEL_EPSILON_EIGENVALUES.dat"
  end if


  open(file=filename,status=files_status_new,unit=io_unit)
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'
  write(io_unit,10) '#'
  write(io_unit,10) '#  This file contains the computed eigenvalues of the static dielectric operator'
  write(io_unit,10) '#  either exact or model, as indicated by the name of this file.'
  write(io_unit,10) '#'
  write(io_unit,10) '#     l                           epsilon_l                                                     '
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'

  do l = 1, lmax

  write(io_unit,20) l, eigenvalues(l)
  end do


  close(io_unit)


end if

10 format(A)
20 format(I7,20X,ES24.16)


end subroutine output_epsilon_eigenvalues
!!***


!!****f* m_gwls_ComputeCorrelationEnergy/output_Sigma_A_by_eigenvalues
!! NAME
!!  output_Sigma_A_by_eigenvalues
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
!!
!! SOURCE

subroutine output_Sigma_A_by_eigenvalues(n_ext_freq,lmax,external_frequencies,AT_Lanczos,eigenvalues_array,which_case)
!----------------------------------------------------------------------------------------------------
!       This routine outputs the eigenvalues of the static dielectric matrix, as well as the
!       contributions to Sigma_A, decomposed by eigenvalues.
!
!       There are three cases to consider:
!                               1 )   no model is being used; we are printing A
!                               2 )   a model is being used; we are printing A1
!                               2 )   a model is being used; we are printing A2
!----------------------------------------------------------------------------------------------------

integer,      intent(in) :: n_ext_freq, lmax, which_case
complex(dpc), intent(in) :: AT_Lanczos(n_ext_freq,lmax)
real(dp),     intent(in) :: eigenvalues_array(lmax)
real(dp),     intent(in) :: external_frequencies(n_ext_freq)

integer   :: iw_ext, l
real(dp)  :: external_omega


complex(dpc)  :: matrix, eig

integer         :: io_unit
character(128) :: filename

! *************************************************************************


if (mpi_enreg%me == 0) then
  io_unit = get_unit()


  do iw_ext = 1, n_ext_freq

  if (which_case == 1) then
    write(filename,'(A,I0.4,A)') "SIGMA_A_BY_EIGENVALUES_",iw_ext,".dat"
  else if (which_case == 2) then
    write(filename,'(A,I0.4,A)') "SIGMA_A1_BY_EIGENVALUES_",iw_ext,".dat"
  else if (which_case == 3) then
    write(filename,'(A,I0.4,A)') "SIGMA_A2_BY_EIGENVALUES_",iw_ext,".dat"
  end if


  external_omega = external_frequencies(iw_ext)

  open(file=filename,status=files_status_new,unit=io_unit)
  write(io_unit,10) '#-----------------------------------------------------------------------------------------------'
  write(io_unit,10) '#'
  write(io_unit,10) '#  This file contains the contributions to Sigma_A, the analytical self-energy term,'
  write(io_unit,10) '#  as a function of the eigenvalue of the dielectric operator.'
  write(io_unit,10) '#'
  write(io_unit,10) '# Definitions:'
  write(io_unit,10) '# '
  write(io_unit,10) '# '

  if (which_case == 1) then
    write(io_unit,10) '#            Sigma_A = Tr[(eps^{-1}(0) -I ) AT(W) ]                           '
    write(io_unit,10) '#                    = sum_{l} ( 1/eps_l -1 ) < V_l | AT(W) | V_l >           '
  else if (which_case == 2) then
    write(io_unit,10) '#            Sigma_A1= Tr[(eps^{-1}(0) - eps_model^{-1}(0) ) A1T(W) ]         '
    write(io_unit,10) '#                    = sum_{l} ( LBDA_l ) < V_l | A1T(W) | V_l >          '
    write(io_unit,10) '#                    where LBDA_l are the eigenvalues of eps^{-1}(0)-eps_model^{-1}(0) '
  else if (which_case == 3) then
    write(io_unit,10) '#            Sigma_A2= Tr[(eps_model^{-1}(0)- I ) A2T(W) ]         '
    write(io_unit,10) '#                    = sum_{l} (1/eps^{model}_l -1 ) < V_l | A2T(W) | V_l >          '
  end if


  write(io_unit,10) '# '
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#  Parameters:                                                               '
  write(io_unit,10) '#                                                                            '

  if (which_case == 1 .or. which_case == 2) then
    write(io_unit,25) '#                 lmax          = ', lmax
  else if (which_case == 3) then
    write(io_unit,25) '#                 lmax_model    = ', lmax
  end if

  write(io_unit,25) '#                 n_ext_freq    = ', n_ext_freq
  write(io_unit,25) '#                 iw_ext        = ', iw_ext
  write(io_unit,14) '#                 omega_ext (W) = ', external_omega,'  Ha'
  write(io_unit,10) '#                                                                            '
  write(io_unit,10) '#                                                                            '
  if (which_case == 1) then
    write(io_unit,10) '#     (1/eps_l -1)                  < V_l | AT(W) | V_l >  (Ha)                ( 1/eps_l -1 ) '&
    &                     //'< V_l | AT(W) | V_l > (Ha)'
  else if (which_case == 2) then
    write(io_unit,10) '#        LBDA_1                     < V_l | A1T(W) | V_l >  (Ha)                   LBDA_l     '&
    &                     //'< V_l | A1T(W) | V_l > (Ha)'
  else if (which_case == 3) then
    write(io_unit,10) '#     (1/eps_m_l -1)                < V_l | A2T(W) | V_l >  (Ha)               ( 1/eps_m_l -1 ) '&
    &                     //'< V_l | A2T(W) | V_l > (Ha)'
  end if

  write(io_unit,10) '#                                 real                 imaginary                    real                  '&
  &                   //'imaginary'
  write(io_unit,10) '#---------------------------------------------------------------------------------------------------------'&
  &                   //'----------------'

  do l = 1, lmax

  matrix = AT_Lanczos(iw_ext,l)
  eig    = eigenvalues_array(l)


  write(io_unit,20) real(eig), matrix, eig*matrix
  end do


  close(io_unit)

  end do

end if

10 format(A)
14 format(A,ES24.16,A)
20 format(ES24.16,2ES24.16,2X,2ES24.16)
25 format(A,I5)



end subroutine output_Sigma_A_by_eigenvalues
!!***


end module m_gwls_ComputeCorrelationEnergy
!!***

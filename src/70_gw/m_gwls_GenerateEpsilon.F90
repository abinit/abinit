!!****m* ABINIT/m_gwls_GenerateEpsilon
!! NAME
!! m_gwls_GenerateEpsilon
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


module m_gwls_GenerateEpsilon
!----------------------------------------------------------------------------------------------------
! This module contains routines to compute and store the dielectric matrix, which plays a
! central role in the self energy computations. In particular, global arrays are used to store
! the static dielectric matrix.
!----------------------------------------------------------------------------------------------------
! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_TimingLog
use m_gwls_polarisability
use m_gwls_model_polarisability
use m_gwls_GWlanczos
! Abinit modules
use m_abicore
use defs_basis
use m_dtset

use m_io_tools,    only : get_unit

implicit none
save
private
!!***

! Global arrays

real(dp), public, allocatable  :: epsilon_eigenvalues_0(:)     ! eigenvalues of the static dielectric matrix

complex(dpc), public, allocatable  :: epsilon_inverse_0(:,:)   ! eps^{-1}-1 in diagonal basis

integer, public ::  kmax, nseeds, lmax
integer, public ::  first_seed
!!***

public :: driver_generate_dielectric_matrix
public :: GeneratePrintDielectricEigenvalues
public :: Driver_GeneratePrintDielectricEigenvalues
!!***

contains

!!****f* m_gwls_GenerateEpsilon/driver_generate_dielectric_matrix
!! NAME
!!  driver_generate_dielectric_matrix
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
!!      cpu_time,diagonalize_lanczos_banded
!!      driver_invert_positive_definite_hermitian_matrix
!!      generateprintdielectriceigenvalues
!!      matrix_function_epsilon_model_operator
!!      set_dielectric_function_frequency,setup_pk_model,write_timing_log,zgemm
!!      zheevd
!!
!! SOURCE

subroutine driver_generate_dielectric_matrix(epsilon_matrix_function,nseeds,kmax,&
epsilon_eigenvalues,Lbasis,debug)
!----------------------------------------------------------------------
! This routine computes the Lanczos approximate representation of the
! implicit dielectic operator and then diagonalizes the banded
! Lanczos matrix.
!----------------------------------------------------------------------
interface
  subroutine epsilon_matrix_function(v_out,v_in,l)
  use defs_basis

  integer,     intent(in)  :: l
  complex(dp), intent(out) :: v_out(l)
  complex(dp), intent(in)  :: v_in(l)

  end subroutine epsilon_matrix_function
end interface

integer,       intent(in) :: nseeds, kmax
logical,       intent(in) :: debug

real   (dp),  intent(out) :: epsilon_eigenvalues(nseeds*kmax)
complex(dpc), intent(out) :: Lbasis(npw_k,nseeds*kmax)  ! array containing the Lanczos basis


! local variables

complex(dpc), allocatable :: seeds(:,:)

complex(dpc),allocatable :: alpha(:,:,:)
complex(dpc),allocatable :: beta (:,:,:)

integer :: mpi_communicator

! *************************************************************************


! The epsilon operator will act in LA mode.
mpi_communicator = mpi_enreg%comm_bandfft


!Create seeds
ABI_ALLOCATE(seeds,(npw_k,nseeds))
call get_seeds(first_seed, nseeds, seeds)

! compute the Lanczos basis
ABI_ALLOCATE(alpha,(nseeds,nseeds,kmax))
ABI_ALLOCATE(beta ,(nseeds,nseeds,kmax))

call block_lanczos_algorithm(mpi_communicator,epsilon_matrix_function,kmax,nseeds,npw_k,        &
seeds,alpha,beta,Lbasis)

! Diagonalize the epsilon matrix, which is banded
call diagonalize_lanczos_banded(kmax,nseeds,npw_k,alpha,beta,Lbasis,epsilon_eigenvalues,debug)

if (debug) then
  call ritz_analysis_general(mpi_communicator, epsilon_matrix_function,nseeds*kmax,npw_k,Lbasis,epsilon_eigenvalues)
end if

ABI_DEALLOCATE(seeds)
ABI_DEALLOCATE(alpha)
ABI_DEALLOCATE(beta)

end subroutine driver_generate_dielectric_matrix
!!***

!!****f* m_gwls_GenerateEpsilon/GeneratePrintDielectricEigenvalues
!! NAME
!!  GeneratePrintDielectricEigenvalues
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_GenerateEpsilon
!!
!! CHILDREN
!!      cpu_time,diagonalize_lanczos_banded
!!      driver_invert_positive_definite_hermitian_matrix
!!      generateprintdielectriceigenvalues
!!      matrix_function_epsilon_model_operator
!!      set_dielectric_function_frequency,setup_pk_model,write_timing_log,zgemm
!!      zheevd
!!
!! SOURCE

subroutine GeneratePrintDielectricEigenvalues(epsilon_matrix_function,kmax,output_filename,Lbasis,alpha,beta)
!----------------------------------------------------------------------
! This routine computes the Lanczos approximate representation of the
! implicit dielectic operator and then diagonalizes the banded
! Lanczos matrix.
!----------------------------------------------------------------------
interface
  subroutine epsilon_matrix_function(v_out,v_in,l)
  use defs_basis

  integer,     intent(in)  :: l
  complex(dp), intent(out) :: v_out(l)
  complex(dp), intent(in)  :: v_in(l)

  end subroutine epsilon_matrix_function
end interface

integer,       intent(in) :: kmax

character(*),  intent(in) :: output_filename


complex(dpc), intent(out) :: Lbasis(npw_k,nseeds*kmax)
complex(dpc), intent(out) :: alpha(nseeds,nseeds,kmax)
complex(dpc), intent(out) :: beta (nseeds,nseeds,kmax)


! local variables


complex(dpc),allocatable :: seeds(:,:)

complex(dpc),allocatable :: Lbasis_diag(:,:)


real(dp),    allocatable :: psik(:,:)
real(dp),    allocatable :: psir(:,:,:,:)

real(dp),    allocatable :: epsilon_eigenvalues(:)


integer :: mpi_communicator
integer :: io_unit
integer :: lmax
integer :: l
integer :: ir1, ir2, ir3
integer :: n1, n2, n3

real(dp) :: R, G
real(dp) :: sigma_R, sigma_G
real(dp) :: x, y, z

real(dp),allocatable :: G_array(:)
real(dp),allocatable :: R_array(:,:,:)

logical :: debug

! *************************************************************************


debug = .false.
lmax = kmax*nseeds
mpi_communicator = mpi_enreg%comm_bandfft
!Create seeds
ABI_ALLOCATE(seeds,(npw_k,nseeds))
call get_seeds(first_seed, nseeds, seeds)

! compute the Lanczos basis
ABI_ALLOCATE(Lbasis_diag,(npw_k,lmax))
ABI_ALLOCATE(epsilon_eigenvalues,(lmax))

ABI_ALLOCATE(psik,(2,npw_k))
ABI_ALLOCATE(psir,(2,n4,n5,n6))
ABI_ALLOCATE(G_array,(npw_k))
ABI_ALLOCATE(R_array,(n4,n5,n6))

psir = zero
R_array = zero

n1 = n4-1
n2 = n5-1
n3 = n6

! Generate the Lanczos basis and banded eigenvalue representation
call block_lanczos_algorithm(mpi_communicator, epsilon_matrix_function,kmax,nseeds,npw_k,        seeds,alpha,beta,Lbasis)

Lbasis_diag = Lbasis

! Diagonalize the epsilon matrix, which is banded
call diagonalize_lanczos_banded(kmax,nseeds,npw_k,alpha,beta,Lbasis_diag,epsilon_eigenvalues,debug)

call ritz_analysis_general(mpi_communicator, epsilon_matrix_function,lmax,npw_k,Lbasis_diag,epsilon_eigenvalues)

io_unit = get_unit()
open(file=output_filename,status=files_status_new,unit=io_unit)
write(io_unit,10) '#----------------------------------------------------------------------------'
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                       Partial eigenvalues                                  '
write(io_unit,10) '#           ==========================================                       '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#  Tabulate the eigenvalues of the dielectic matrix, as well as some         '
write(io_unit,10) '#  information regarding the eigenstates.                                    '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#  definitions:                                                              '
write(io_unit,10) '#                l      index of the eigenvalue                              '
write(io_unit,10) '#                eig    eigenvalue                                           '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                (the following vectors are expressed in crystal units)      '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                R       =   < l | |r| | l >                                 '
write(io_unit,10) '#                sigma_R = sqrt{< l | (|r|-R)^2 | l > }                      '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                G       =   < l | |G| | l >                                 '
write(io_unit,10) '#                sigma_G = sqrt{< l | (|G|-G)^2 | l > }                      '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#  l            eig                r         sigma_r       G         sigma_G '
write(io_unit,10) '#----------------------------------------------------------------------------'
flush(io_unit)


G_array(:) = kg_k(1,:)**2+ kg_k(2,:)**2+ kg_k(3,:)**2

G_array(:) = sqrt(G_array(:))

R_array  = zero

do ir3=1,n3

if (ir3 <= n3/2 ) then
  z = (one*ir3)/(one*n3)
else
  z = (one*ir3)/(one*n3)-one
end if

do ir2=1,n2

if (ir2 <= n2/2 ) then
  y = (one*ir2)/(one*n2)
else
  y = (one*ir2)/(one*n2)-one
end if

do ir1=1,n1

if (ir1 <= n1/2 ) then
  x = (one*ir1)/(one*n1)
else
  x = (one*ir1)/(one*n1)-one
end if


R_array(ir1,ir2,ir3)  = sqrt(x**2+y**2+z**2)
end do
end do
end do



do l=1, lmax

psik(1,:) = dble (Lbasis_diag(:,l))
psik(2,:) = dimag(Lbasis_diag(:,l))

call g_to_r(psir ,psik)


G       = sum(G_array(:)*(psik(1,:)**2+psik(2,:)**2))
R       = sum(R_array(:,:,:)*(psir(1,:,:,:)**2+psir(2,:,:,:)**2) )*ucvol/nfft

sigma_G = sqrt(sum((G_array(:)    -G)**2*(psik(1,:)**2    +psik(2,:)**2)))
sigma_R = sqrt(sum((R_array(:,:,:)-R)**2*(psir(1,:,:,:)**2+psir(2,:,:,:)**2))*ucvol/nfft)



write(io_unit,20) l, epsilon_eigenvalues(l), R,sigma_R, G,sigma_G

end do

close(io_unit)

ABI_DEALLOCATE(seeds)
ABI_DEALLOCATE(Lbasis_diag)
ABI_DEALLOCATE(psik)
ABI_DEALLOCATE(psir)
ABI_DEALLOCATE(G_array)
ABI_DEALLOCATE(R_array)

ABI_DEALLOCATE(epsilon_eigenvalues)


10 format(A)
20 format(I5,ES24.16,4F12.8)

end subroutine GeneratePrintDielectricEigenvalues
!!***

!!****f* m_gwls_GenerateEpsilon/Driver_GeneratePrintDielectricEigenvalues
!! NAME
!!  Driver_GeneratePrintDielectricEigenvalues
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
!!      cpu_time,diagonalize_lanczos_banded
!!      driver_invert_positive_definite_hermitian_matrix
!!      generateprintdielectriceigenvalues
!!      matrix_function_epsilon_model_operator
!!      set_dielectric_function_frequency,setup_pk_model,write_timing_log,zgemm
!!      zheevd
!!
!! SOURCE

subroutine Driver_GeneratePrintDielectricEigenvalues(dtset)
!----------------------------------------------------------------------
! Compute the eigenvalues of the various dielectric operators
!----------------------------------------------------------------------
type(dataset_type),intent(in) :: dtset

integer  ::  kmax_exact, kmax_model, kmax
real(dp) :: second_model_parameter


integer  ::  lm, k, lmax, l1, l2
integer  ::  io_unit
integer  ::  io_unit2

real(dp) :: time1, time2, time
! local variables
character(128)  :: output_filename
character(256)  :: timing_string


complex(dpc), allocatable :: Lbasis_exact(:,:)
complex(dpc), allocatable :: Lbasis_model(:,:)

complex(dpc), allocatable :: sub_Lbasis_exact(:,:)
complex(dpc), allocatable :: sub_Lbasis_model(:,:)

complex(dpc), allocatable :: dummy(:,:)
complex(dpc), allocatable :: dummy2(:,:)
complex(dpc), allocatable :: dummy3(:,:)

complex(dpc), allocatable :: alpha_exact(:,:,:)
complex(dpc), allocatable :: beta_exact (:,:,:)

complex(dpc), allocatable :: alpha_model(:,:,:)
complex(dpc), allocatable :: beta_model (:,:,:)

real(dp), allocatable :: eig_exact(:)
real(dp), allocatable :: eig_model(:)

complex(dpc), allocatable :: model_epsilon_matrix(:,:)
complex(dpc), allocatable :: vector(:)

real(dpc) :: tr_eps_1, tr_eps_2, tr_eps_3


integer   ::  lwork, lrwork, liwork, info
complex(dpc), allocatable :: work(:)
real(dp)    , allocatable :: rwork(:)
integer     , allocatable :: iwork(:)

integer        :: debug_unit
character(50)  :: debug_filename

! *************************************************************************




kmax_exact   = dtset%gwls_stern_kmax
kmax_model   = dtset%gwls_kmax_complement

!second_model_parameter  = dtset%gwls_second_model_parameter
second_model_parameter  = zero


! global stuff
nseeds       = dtset%gwls_nseeds
first_seed   = dtset%gwls_first_seed
e            = dtset%gwls_band_index



ABI_ALLOCATE(Lbasis_exact,(npw_k,kmax_exact*nseeds))
ABI_ALLOCATE(Lbasis_model,(npw_k,kmax_model*nseeds))

ABI_ALLOCATE(alpha_exact, (nseeds,nseeds,kmax_exact))
ABI_ALLOCATE(beta_exact , (nseeds,nseeds,kmax_exact))
ABI_ALLOCATE(alpha_model, (nseeds,nseeds,kmax_model))
ABI_ALLOCATE(beta_model , (nseeds,nseeds,kmax_model))


! set omega=0 for exact dielectric operator
call set_dielectric_function_frequency([zero,zero])



call cpu_time(time1)
output_filename = 'EIGENVALUES_EXACT.dat'
call GeneratePrintDielectricEigenvalues(matrix_function_epsilon_k, kmax_exact, &
output_filename, Lbasis_exact, alpha_exact, beta_exact)



call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to compute the EXACT Static Dielectric Matrix  :   "
call write_timing_log(timing_string,time)



call cpu_time(time1)
call setup_Pk_model(zero,second_model_parameter)
output_filename = 'EIGENVALUES_MODEL.dat'
call GeneratePrintDielectricEigenvalues(matrix_function_epsilon_model_operator, kmax_model, output_filename,&
Lbasis_model, alpha_model, beta_model)
call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to compute the MODEL Static Dielectric Matrix  :   "
call write_timing_log(timing_string,time)


call cpu_time(time1)
if (kmax_exact <= kmax_model) then
  kmax  = kmax_exact
else
  kmax  = kmax_model
end if
lmax = nseeds*kmax

! Build model operator matrix elements in the exact basis
ABI_ALLOCATE(model_epsilon_matrix, (lmax,lmax))
ABI_ALLOCATE(vector, (npw_k))

model_epsilon_matrix = cmplx_0

do l2 =1 , lmax
call matrix_function_epsilon_model_operator(vector ,Lbasis_exact(:,l2),npw_k)


do l1 =1, lmax

model_epsilon_matrix(l1, l2) = complex_vector_product(Lbasis_exact(:,l1),vector,npw_k)

end do


end do

ABI_DEALLOCATE(vector)

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Compute MODEL matrix elements in EXACT basis        :   "
call write_timing_log(timing_string,time)




! Compare the traces


io_unit = get_unit()
open(file='DIELECTRIC_TRACE.dat',status=files_status_new,unit=io_unit)
write(io_unit,10) '#----------------------------------------------------------------------------'
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                       Partial traces                                       '
write(io_unit,10) '#           ==========================================                       '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#  Tabulate the trace of various operators, as function of the number        '
write(io_unit,10) '#  of lanczos steps performed.                                               '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#  NOTES:                                                                    '
write(io_unit,10) '#         Tr[1-eps^{-1}] is evaluated in the Lanczos basis of eps            '
write(io_unit,10) '#         Tr[1-eps_m^{-1}] is evaluated in the Lanczos basis of eps_m        '
write(io_unit,10) '#         Tr[eps_m^{-1}-eps^{-1}] is evaluated in the Lanczos basis of eps   '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#                                                                            '
write(io_unit,10) '#  k       Tr[1-eps^{-1}]        Tr[1-eps_m^{-1}]     Tr[eps_m^{-1}-eps^{-1}]'
write(io_unit,10) '#----------------------------------------------------------------------------'
flush(io_unit)


io_unit2 = get_unit()
open(file='RPA_ENERGY.dat',status=files_status_new,unit=io_unit2)
write(io_unit2,10) '#----------------------------------------------------------------------------'
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#                       RPA TOTAL ENERGY                                     '
write(io_unit2,10) '#           ==========================================                       '
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#  It can be shown that the correlation energy, within the RPA, is given     '
write(io_unit2,10) '#  by:                                                                       '
write(io_unit2,10) '#       E_c = int_0^{infty} dw/(2pi) Tr[ ln(eps{iw)}+1-eps(iw)]              '
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#  As a gauge of what can be expected as far as convergence is concerned,    '
write(io_unit2,10) '#  the following will be printed.                                            '
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#  I_1 = Tr[  ln(eps)  + 1 - eps   ]                                         '
write(io_unit2,10) '#  I_2 = Tr[ ln(eps_m) + 1 - eps_m ]                                         '
write(io_unit2,10) '#  I_3 = Tr[ ln(eps^{-1}_m . eps ) + eps_m - eps ]                           '
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#                                                                            '
write(io_unit2,10) '#  k       I_1                   I_2                  I_3                    '
write(io_unit2,10) '#----------------------------------------------------------------------------'
flush(io_unit2)


! Iterate every 10 values of k max, or else the linear algebra gets too expensive...
do k = 4, kmax, 4

ABI_ALLOCATE(sub_Lbasis_exact,(npw_k,k*nseeds))
ABI_ALLOCATE(sub_Lbasis_model,(npw_k,k*nseeds))


ABI_ALLOCATE(eig_exact,(k*nseeds))
ABI_ALLOCATE(eig_model,(k*nseeds))
ABI_ALLOCATE(dummy,(k*nseeds,k*nseeds))
ABI_ALLOCATE(dummy2,(k*nseeds,k*nseeds))

sub_Lbasis_exact(:,:) = Lbasis_exact(:,1:k*nseeds)
sub_Lbasis_model(:,:) = Lbasis_model(:,1:k*nseeds)

! Diagonalize the epsilon matrix, which is banded
call diagonalize_lanczos_banded(k,nseeds,npw_k,alpha_exact(:,:,1:k),beta_exact(:,:,1:k),sub_Lbasis_exact,eig_exact,.false.)
call diagonalize_lanczos_banded(k,nseeds,npw_k,alpha_model(:,:,1:k),beta_model(:,:,1:k),sub_Lbasis_model,eig_model,.false.)

tr_eps_1 = sum(one-one/eig_exact(:))
tr_eps_2 = sum(one-one/eig_model(:))

tr_eps_3 = -sum(one/eig_exact(:))

dummy(:,:) = model_epsilon_matrix(1:k*nseeds, 1:k*nseeds)

call driver_invert_positive_definite_hermitian_matrix(dummy,k*nseeds)

do lm = 1, k*nseeds

tr_eps_3 = tr_eps_3 +  dble(dummy(lm,lm))
end do


write(io_unit,20) k, tr_eps_1, tr_eps_2, tr_eps_3
flush(io_unit)



tr_eps_1 = sum(log(eig_exact(:))+one-eig_exact(:))
tr_eps_2 = sum(log(eig_model(:))+one-eig_model(:))

dummy2(:,:) = zero
tr_eps_3    = zero
do lm = 1, k*nseeds
dummy2(lm,lm) = eig_exact(lm)
tr_eps_3 = tr_eps_3 + dble(model_epsilon_matrix(lm,lm)) - dble(eig_exact(lm))
end do

ABI_ALLOCATE(dummy3,(k*nseeds,k*nseeds))
call ZGEMM(      'N',   & ! Hermitian conjugate the first array
'N',   & ! Leave second array as is
k*nseeds,   & ! the number of rows of the  matrix op( A )
k*nseeds,   & ! the number of columns of the  matrix op( B )
k*nseeds,   & ! the number of columns of the  matrix op( A ) == rows of matrix op( B )
cmplx_1,   & ! alpha constant
dummy2,   & ! matrix A
k*nseeds,   & ! LDA
dummy,   & ! matrix B
k*nseeds,   & ! LDB
cmplx_0,   & ! beta constant
dummy3,   & ! matrix C
k*nseeds)     ! LDC

dummy2(:,:) = dummy3(:,:)
ABI_DEALLOCATE(dummy3)

! find eigenvalues
!call heevd(dummy2, eig_exact)

lwork  = k*nseeds+1
lrwork = k*nseeds
liwork = 1

ABI_ALLOCATE(work,(lwork))
ABI_ALLOCATE(rwork,(lrwork))
ABI_ALLOCATE(iwork,(liwork))

call zheevd('N', 'U',k*nseeds, dummy2, k*nseeds, eig_exact, work, lwork, rwork, lrwork, iwork, liwork, info)
if ( info /= 0) then
  debug_unit = get_unit()
  write(debug_filename,'(A,I4.4,A)') 'LAPACK_DEBUG_PROC=',mpi_enreg%me,'.log'

  open(debug_unit,file=trim(debug_filename),status='unknown')

  write(debug_unit,'(A)')      '*************************************************************************************'
  write(debug_unit,'(A,I4,A)') '*      ERROR: info = ',info,' in ZHEEVD(1), gwls_GenerateEpsilon'
  write(debug_unit,'(A)')      '*************************************************************************************'

  close(debug_unit)

end if





tr_eps_3 =  tr_eps_3 + sum(log(eig_exact))

write(io_unit2,20) k, tr_eps_1, tr_eps_2, tr_eps_3
flush(io_unit2)



ABI_DEALLOCATE(work)
ABI_DEALLOCATE(rwork)
ABI_DEALLOCATE(iwork)



ABI_DEALLOCATE(sub_Lbasis_exact)
ABI_DEALLOCATE(sub_Lbasis_model)

ABI_DEALLOCATE(eig_exact)
ABI_DEALLOCATE(eig_model)
ABI_DEALLOCATE(dummy)
ABI_DEALLOCATE(dummy2)
end do

close(io_unit)
close(io_unit2)

call cpu_time(time2)
time = time2-time1
write(timing_string,'(A)')  "Time to compute the TRACES of the Dielectric Matrices:   "
call write_timing_log(timing_string,time)

ABI_DEALLOCATE(Lbasis_exact)
ABI_DEALLOCATE(Lbasis_model)
ABI_DEALLOCATE(alpha_exact)
ABI_DEALLOCATE(beta_exact )
ABI_DEALLOCATE(alpha_model)
ABI_DEALLOCATE(beta_model )
ABI_DEALLOCATE(model_epsilon_matrix)



10 format(A)
20 format(I5,3ES24.16)

end subroutine Driver_GeneratePrintDielectricEigenvalues
!!***

end module m_gwls_GenerateEpsilon
!!***

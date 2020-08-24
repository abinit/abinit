!!****m* ABINIT/m_gwls_ComputePoles
!! NAME
!! m_gwls_ComputePoles
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

module m_gwls_ComputePoles

use m_gwls_utility
use m_gwls_wf
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_polarisability
use m_gwls_GWlanczos
use m_gwls_GenerateEpsilon
use m_gwls_GWanalyticPart
use m_gwls_TimingLog
use m_gwls_LanczosBasis

use defs_basis
use defs_wvltypes
use m_abicore
use m_xmpi
use m_errors

use m_io_tools,         only : get_unit


implicit none
save
private

integer  :: number_of_denerate_sets
integer  :: largest_degeneracy

integer, allocatable :: degeneracy_table(:,:)
integer, allocatable :: number_of_degenerate_states(:)

real(dp) :: En_m_omega_2

public :: compute_Poles
public :: generate_degeneracy_table_for_poles
public :: clean_degeneracy_table_for_poles

CONTAINS
!!***

!!****f* m_gwls_ComputePoles/generate_degeneracy_table_for_poles
!! NAME
!!  generate_degeneracy_table_for_poles
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
!!      block_lanczos_algorithm,diagonalize_lanczos_banded
!!      ritz_analysis_general,xmpi_sum
!!
!! SOURCE

subroutine generate_degeneracy_table_for_poles(debug)
!----------------------------------------------------------------------
! This subroutine groups, once and for all, the indices of
! degenerate eigenstates. This will be useful to compute the poles
!
!----------------------------------------------------------------------

logical, intent(in) :: debug

real(dp) :: degeneracy_tolerance

integer  :: nbands
integer  :: n, i
integer  :: i_set, j_deg
integer  :: n_degeneracy

real(dp) :: energy

integer        :: io_unit
character(128) :: filename
logical        :: file_exists

! *************************************************************************


degeneracy_tolerance = 1.0D-8
!--------------------------------------------------------------------------------
!
! First, find the largest degeneracy in the eigenvalue spectrum
!
!--------------------------------------------------------------------------------

nbands = size(eig)

! initialize
largest_degeneracy      = 1
number_of_denerate_sets = 1

n_degeneracy = 1
energy       = eig(1)

!============================================================
! Notes for the code block below:
!
! The electronic eigenvalues can be grouped in degenerate
! sets. The code below counts how many such sets there are,
! and how many eigenvalues belong to each set.
!
! It is important to have a robust algorithm to do this
! properly. In particular, the edge case where the LAST SET
! is degenerate must be treated with care .

! The algorithm.I'm looking at at the time of this writing
! has a bug in it and cannot handle a degenerate last set...
! Let's fix that!
!============================================================

do n = 2, nbands
  if (abs(eig(n) - energy) < degeneracy_tolerance) then
    ! A degenerate state! add one
    n_degeneracy  = n_degeneracy + 1
  else
    ! We are no longer degenerate. Update
    if ( n_degeneracy > largest_degeneracy ) largest_degeneracy = n_degeneracy

    n_degeneracy = 1
    energy       = eig(n)
    number_of_denerate_sets = number_of_denerate_sets + 1
  end if

  ! If this is the last index, update the largest_degeneracy if necessary
  if ( n == nbands .and. n_degeneracy > largest_degeneracy ) largest_degeneracy = n_degeneracy

end do

!--------------------------------------------------------------------------------
!
! Allocate the array which will contain the indices of the degenerate
! states, and populate it.
!--------------------------------------------------------------------------------
ABI_ALLOCATE(degeneracy_table,            (number_of_denerate_sets,largest_degeneracy))
ABI_ALLOCATE(number_of_degenerate_states, (number_of_denerate_sets))

degeneracy_table(:,:) = 0

i_set = 1
j_deg = 1

! initialize
energy = eig(1)

degeneracy_table(i_set,j_deg)       = 1
number_of_degenerate_states(i_set)  = 1

do n = 2, nbands

  if (abs(eig(n) - energy) < degeneracy_tolerance) then
    ! A degenerate state! add one
    j_deg = j_deg + 1

  else
    ! We are no longer degenerate. Update
    j_deg = 1
    i_set = i_set+1
    energy = eig(n)

  end if


  number_of_degenerate_states(i_set) = j_deg
  degeneracy_table(i_set,j_deg)      = n

end do

if (debug .and. mpi_enreg%me == 0) then
  io_unit  = get_unit()
  filename = "degeneracy_table.log"

  i = 0
  inquire(file=filename,exist=file_exists)
  do while (file_exists)
  i = i+1
  write (filename,'(A,I0,A)') "degeneracy_table_",i,".log"
  inquire(file=filename,exist=file_exists)
  end do

  io_unit = get_unit()

  open(io_unit,file=filename,status=files_status_new)

  write(io_unit,10) " "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,10) "#                     Degeneracy table : tabulate the degenerate states                         "
  write(io_unit,10) "#                     -------------------------------------------------------                   "
  write(io_unit,10) "#                                                                                               "
  write(io_unit,14) "#     number_of_denerate_sets = ", number_of_denerate_sets
  write(io_unit,10) "#                                                                                               "
  write(io_unit,14) "#      largest_degeneracy     = ", largest_degeneracy
  write(io_unit,10) "#                                                                                               "
  write(io_unit,10) "# Eigenvalues (Ha)                                                                              "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,16) eig(:)



  write(io_unit,10) "#==============================================================================================="
  write(io_unit,10) "#  i_set    number of states           States                                                   "
  write(io_unit,10) "#==============================================================================================="
  flush(io_unit)

  do i_set = 1, number_of_denerate_sets
  write(io_unit,12) i_set, number_of_degenerate_states(i_set), degeneracy_table(i_set,:)
  end do

  flush(io_unit)

  close(io_unit)
end if

10 format(A)
12 format(I5,10X,I5,15X,1000I5)
14 format(A,I5)
16 format(1000F12.8,2X)

end subroutine generate_degeneracy_table_for_poles
!!***

!!****f* m_gwls_ComputePoles/clean_degeneracy_table_for_poles
!! NAME
!!  clean_degeneracy_table_for_poles
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
!!      block_lanczos_algorithm,diagonalize_lanczos_banded
!!      ritz_analysis_general,xmpi_sum
!!
!! SOURCE

subroutine clean_degeneracy_table_for_poles()

! *************************************************************************

if(allocated(degeneracy_table)) then
  ABI_DEALLOCATE(degeneracy_table)
end if
if(allocated(number_of_degenerate_states)) then
  ABI_DEALLOCATE(number_of_degenerate_states)
end if

end subroutine clean_degeneracy_table_for_poles
!!***

!!****f* m_gwls_ComputePoles/compute_Poles
!! NAME
!!  compute_Poles
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
!!
!! CHILDREN
!!
!!
!! SOURCE

function compute_Poles(external_omega,kmax_poles,debug)
!----------------------------------------------------------------------
! This function extract the Pole contributions to the correlation
! energy, as a function of the external frequency.
!
! The algorithm builds a Lanczos chain for each contributing subspace;
! the number of steps is controlled by kmax_poles.
!
! This function will take in explicit arguments, as it is simpler
! to do this than to define global arrays.
!----------------------------------------------------------------------
real(dp) :: compute_Poles

real(dp),     intent(in) :: external_omega
integer,      intent(in) :: kmax_poles
logical,      intent(in) :: debug

real(dp) :: energy_tolerance
real(dp) :: pole_contribution



integer :: number_of_seeds
integer :: i_set
integer :: n

real(dp)     :: prefactor
real(dp)     :: En_m_omega

logical      :: pole_is_valence
logical      :: pole_is_conduction
logical      :: pole_is_in_gap

integer        :: io_unit, i
character(128) :: filename
logical        :: file_exists



complex(dpc), allocatable :: seeds(:,:)

! *************************************************************************

compute_Poles    =  zero

energy_tolerance = 1.0D-8


if (debug .and. mpi_enreg%me == 0) then
  io_unit = get_unit()

  i = 0

  file_exists = .true.
  do while (file_exists)
  i = i+1
  write (filename,'(A,I0.4,A)') "ComputePoles_",i,".log"
  inquire(file=filename,exist=file_exists)
  end do


  open(io_unit,file=filename,status=files_status_new)

  write(io_unit,10) " "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,10) "#                     ComputePoles: debug information for the pole computation                  "
  write(io_unit,10) "#                     --------------------------------------------------------                  "
  write(io_unit,10) "#                                                                                               "
  write(io_unit,10) "# This file contains data describing the computation of the pole contribution to the            "
  write(io_unit,10) "# correlation self energy.                                                                      "
  write(io_unit,10) "#                                                                                               "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,10) " "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,10) "#                                                                                               "
  write(io_unit,10) "#  parameters:                                                                                  "
  write(io_unit,10) "#                                                                                               "
  write(io_unit,11) "#          external_omega  = ",external_omega," Ha                                              "
  write(io_unit,12) "#               kmax_poles = ",kmax_poles
  write(io_unit,12) "#               nbandv     = ",nbandv
  write(io_unit,10) "#                                                                                               "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,10) "#                                                                                               "
  write(io_unit,10) "#   DFT Eigenvalues (Ha)                                                                        "
  write(io_unit,10) "#==============================================================================================="
  write(io_unit,13) eig(:)
  flush(io_unit)

end if



!--------------------------------------------------------------------------------
!
! Determine if external frequency corresponds to the valence or the conduction
! manifold.
!
!--------------------------------------------------------------------------------

compute_Poles = zero

pole_is_conduction = .false.
pole_is_valence    = .false.
pole_is_in_gap     = .false.



! Careful here! there may be only nbandv states in memory; nbandv+1 causes segfaults!
if ( external_omega <= eig(nbandv)) then
  pole_is_valence    = .true.
else if ( external_omega > eig(nbandv)) then
  pole_is_conduction = .true.
else
  pole_is_in_gap = .true.
end if




if (debug .and. mpi_enreg%me == 0 ) then
  write(io_unit,10) "#===================================================================================================="
  write(io_unit,10) "#                                                                                                    "
  write(io_unit,10) "#  Determine where the external energy is:                                                           "
  write(io_unit,10) "#                                                                                                    "
  write(io_unit,14) "#             pole_is_valence    = ",pole_is_valence
  write(io_unit,14) "#             pole_is_conduction = ",pole_is_conduction
  write(io_unit,14) "#             pole_is_in_gap     = ",pole_is_in_gap
  write(io_unit,10) "#                                                                                                    "
  write(io_unit,10) "#===================================================================================================="
  flush(io_unit)

end if

if ( pole_is_in_gap) return

!--------------------------------------------------------------------------------
!
! Loop on all degenerate sets
!
!--------------------------------------------------------------------------------

if (debug .and. mpi_enreg%me == 0 ) then
  write(io_unit,10) "#===================================================================================================="
  write(io_unit,10) "#                                                                                                    "
  write(io_unit,10) "#  Iterating over all degenerate sets of eigenvalues:                                                "
  write(io_unit,10) "#                                                                                                    "
  write(io_unit,10) "#===================================================================================================="
  flush(io_unit)

end if


do i_set =1, number_of_denerate_sets

n = degeneracy_table(i_set,1)

En_m_omega = eig(n)-external_omega

if (debug .and. mpi_enreg%me == 0) then
  write(io_unit,12) "# i_set = ", i_set
  write(io_unit,12) "#                           n = ",n
  write(io_unit,16) "#                eig(n)-omega = ",En_m_omega," Ha"
  flush(io_unit)
end if

!------------------------------------------
! Test if we need to exit the loop
!------------------------------------------
if (pole_is_valence ) then

  ! If the pole is valence, get out when we enter conduction states
  if (  n > nbandv  ) then
    if (debug.and. mpi_enreg%me == 0) then
      write(io_unit,10) "#"
      write(io_unit,10) "#                n > nbandv : exit loop!"
      flush(io_unit)
    end if

    exit
  end if

  ! if the valence energy is smaller than the external frequency,
  ! then there is no contribution

  if (En_m_omega < zero  .and. abs(En_m_omega) > energy_tolerance) then
    ! careful close to zero!
    if (debug .and. mpi_enreg%me == 0) then
      write(io_unit,10) "# "
      write(io_unit,10) "#                 eig(n) < omega : cycle!"
      flush(io_unit)
    end if
    cycle
  end if


  ! if we are still here, there is a valence contribution
  prefactor = -one

else if ( pole_is_conduction ) then

  ! If the pole is conduction, get out when the conduction state is
  ! larger than the frequency (careful close to zero!)
  if ( En_m_omega > energy_tolerance ) then
    if (debug .and. mpi_enreg%me == 0) then
      write(io_unit,10) "#"
      write(io_unit,10) "#                eig(n) > omega : exit!"
      flush(io_unit)
    end if
    exit
  end if

  ! If the pole is conduction, there is no contribution while
  ! we are in the valence states
  if (  n <= nbandv  ) then
    if (debug .and. mpi_enreg%me == 0) then
      write(io_unit,10) "#"
      write(io_unit,10) "#                n <= nbandv : cycle!"
      flush(io_unit)
    end if

    cycle
  end if

  ! if we are still here, there is a conduction contribution
  prefactor = one
end if


!-------------------------------------------------
! If we made it this far, we have a contribution!
!-------------------------------------------------

if (abs(En_m_omega) < energy_tolerance ) then

  if (debug .and. mpi_enreg%me == 0) then
    write(io_unit,10) "# "
    write(io_unit,10) "#                En - omega ~ 0: pole at the origin, multiply by 1/2!"
    flush(io_unit)
  end if

  ! The factor of 1/2 accounts for the fact that
  ! the pole is at the origin!
  prefactor = 0.5_dp*prefactor
end if



number_of_seeds = number_of_degenerate_states(i_set)

ABI_ALLOCATE(seeds, (npw_k,number_of_seeds))

call get_seeds(n, number_of_seeds, seeds) !Missing wrappers

call set_dielectric_function_frequency([En_m_omega,zero])
if (debug .and. mpi_enreg%me == 0) then
  write(io_unit,10) "#                Compute pole contribution:"
  write(io_unit,12) "#                        number of seeds = ",number_of_seeds
  write(io_unit,16) "#                        ||   seeds   || = ",sqrt(sum(abs(seeds(:,:))**2))  !Missing xmpi_sum
  write(io_unit,16) "#                        eig(n)-omega    = ",En_m_omega, " Ha"
  write(io_unit,17) "#                        prefactor       = ",prefactor
  flush(io_unit)
end if
if(dtset%zcut > tol12) activate_inf_shift_poles = .true.
En_m_omega_2 = En_m_omega
pole_contribution =                                             &
compute_pole_contribution(matrix_function_epsilon_k,    &
number_of_seeds, kmax_poles,    &
seeds,debug)
if(dtset%zcut > tol12) activate_inf_shift_poles = .false.
if (debug .and. mpi_enreg%me == 0) then
  write(io_unit,16) "#                      pole contribution = ",prefactor*pole_contribution, " Ha"
  flush(io_unit)
end if


compute_Poles = compute_Poles + prefactor*pole_contribution
ABI_DEALLOCATE(seeds)
end do

if (debug .and. mpi_enreg%me == 0) then
  close(io_unit)
end if


10 format(A)
11 format(A,F8.4,A)
12 format(A,I5)
13 format(1000F16.8)
14 format(A,L10)
16 format(A,ES12.4,A)
17 format(A,F8.4)

end function compute_Poles
!!***

!!****f* m_gwls_ComputePoles/compute_pole_contribution
!! NAME
!! compute_pole_contribution
!!
!! FUNCTION
!! .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

function compute_pole_contribution(epsilon_matrix_function,nseeds,kmax,seeds,debug)
!----------------------------------------------------------------------
! This routine computes the contribution to the  poles energy
! coming from the states in the seeds.
!----------------------------------------------------------------------
interface
  subroutine epsilon_matrix_function(v_out,v_in,l)

  use defs_basis

  integer,     intent(in)  :: l
  complex(dp), intent(out) :: v_out(l)
  complex(dp), intent(in)  :: v_in(l)

  end subroutine epsilon_matrix_function
end interface

real(dp) :: compute_pole_contribution

integer,       intent(in) :: nseeds, kmax
complex(dpc),  intent(in) :: seeds(npw_k,nseeds)
logical,       intent(in) :: debug

! local variables

integer  :: mpi_communicator

complex(dpc),allocatable :: local_seeds(:,:)

complex(dpc),allocatable :: Lbasis(:,:)  ! array containing the Lanczos basis
complex(dpc),allocatable :: alpha(:,:,:)
complex(dpc),allocatable :: beta (:,:,:)
real(dp),    allocatable :: epsilon_eigenvalues(:)

real(dp):: matrix_elements

complex(dpc) :: cmplx_value
integer :: l, s
integer :: ierr

! *************************************************************************

! compute the Lanczos basis
ABI_ALLOCATE(alpha,(nseeds,nseeds,kmax))
ABI_ALLOCATE(beta ,(nseeds,nseeds,kmax))
ABI_ALLOCATE(Lbasis,(npw_k,nseeds*kmax))
ABI_ALLOCATE(local_seeds,(npw_k,nseeds))
ABI_ALLOCATE(epsilon_eigenvalues, (nseeds*kmax))


mpi_communicator = mpi_enreg%comm_bandfft !Missing maybe something for easy access of LA and FFT comms?

local_seeds(:,:) = seeds(:,:)

call block_lanczos_algorithm(mpi_communicator, epsilon_matrix_function,kmax,nseeds,npw_k,        &
&                                local_seeds,alpha,beta,Lbasis)

write(std_out,*) "alpha:"
do l=1,kmax
do s=1,nseeds
write(std_out,*) alpha(:,s,l)
end do
write(std_out,*) " "
end do

write(std_out,*) "beta:"
do l=1,kmax
do s=1,nseeds
write(std_out,*) beta(:,s,l)
end do
write(std_out,*) " "
end do

ABI_DEALLOCATE(local_seeds)

! Diagonalize the epsilon matrix, which is banded
call diagonalize_lanczos_banded(kmax,nseeds,npw_k,alpha,beta,Lbasis,epsilon_eigenvalues,debug)

if (debug) then
  call ritz_analysis_general(mpi_communicator ,epsilon_matrix_function,nseeds*kmax,npw_k,Lbasis,epsilon_eigenvalues)
end if

compute_pole_contribution = zero

do l = 1, nseeds*kmax

matrix_elements = zero

do s = 1, nseeds

cmplx_value = complex_vector_product(seeds(:,s),Lbasis(:,l),npw_k)

call xmpi_sum(cmplx_value,mpi_communicator,ierr) ! sum on all processors working on FFT!

matrix_elements = matrix_elements + abs(cmplx_value)**2


end do


compute_pole_contribution = compute_pole_contribution  + &
matrix_elements *(one/epsilon_eigenvalues(l)-one)
end do



ABI_DEALLOCATE(alpha)
ABI_DEALLOCATE(beta)
ABI_DEALLOCATE(Lbasis)
ABI_DEALLOCATE(epsilon_eigenvalues)

end function compute_pole_contribution

end module m_gwls_ComputePoles
!!***

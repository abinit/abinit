!!****m* ABINIT/m_gwls_GWanalyticPart
!! NAME
!! m_gwls_GWanalyticPart
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


module m_gwls_GWanalyticPart
!----------------------------------------------------------------------------------------------------
! This module contains routines to compute the contribution to the GW correlation energy coming
! from the analytic integral of the trial frequency function.
!----------------------------------------------------------------------------------------------------
! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_GWlanczos
use m_gwls_GenerateEpsilon
use m_gwls_LanczosBasis
use m_gwls_model_polarisability
use m_gwls_polarisability
!use m_gwls_ComplementSpacePolarizability
use m_gwls_GWlanczos
use m_gwls_TimingLog

! abinit modules
use defs_basis
use m_abicore

implicit none
save
private
!!***

!real(dp),allocatable :: epsilon_eigenvalues_complement(:)
!real(dp),allocatable :: lanczos_basis_complement(:,:,:)


complex(dpc),public, allocatable :: A_array(:,:)

integer, public  :: model_number
real(dp),public  :: model_parameter

real(dp), public :: G0_model_epsilon_0 ! model parameter
!!***

public :: get_projection_band_indices
!!***

contains

!!****f* m_hamiltonian/get_projection_band_indices
!! NAME
!!  get_projection_band_indices
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
!!
!! SOURCE

subroutine get_projection_band_indices(omega,band_index_below, band_index_above)
!----------------------------------------------------------------------------------------------------
! This subroutine computes the band indices necessary for properly projecting the Sternheimer equations
! where
!                Pe : projection on states such that epsilon_n < omega
!                Qe : projection on states such that epsilon_n > omega
!----------------------------------------------------------------------------------------------------
implicit none

real(dp),intent(in)  :: omega 
integer, intent(out) :: band_index_below, band_index_above

! *************************************************************************

! First, find the indices for the projections

band_index_below = 1
! it is assumed that the eigenvalues are sorted. As soon as the current eigenvalue
! is equal or larger than epsilon_e,  exit!
do
if ( omega - eig(band_index_below) <= 1.0e-8 ) exit
band_index_below = band_index_below + 1

if (band_index_below > size(eig)) then

  write(std_out,*) '************************************************************'
  write(std_out,*) '***    ERROR IN ROUTINE get_projection_band_indices      ***'
  write(std_out,*) '***                                                      ***'
  write(std_out,*) '***    The index of the DFT eigenvalue larger than       ***'
  write(std_out,*) '***    the target frequency is larger than the number    ***'
  write(std_out,*) '***    of explicitly calculated DFT eigenvalues.         ***'
  write(std_out,*) '***    The computation cannot go on to produce a         ***'
  write(std_out,*) '***    meaningful result: review your input!             ***'
  write(std_out,*) '***                                                      ***'
  write(std_out,*) '***               program stops.                         ***'
  write(std_out,*) '************************************************************'
  stop
end if
end do
! We have overshooted in order to exit the do loop , so decrement by 1
band_index_below = band_index_below - 1

! band_index_below  is now the index of the highest eigenvalue BELOW omega
! NOTE THAT band_index_below = 0 if omega is smaller than all eigenvalues!
! A test should be performed on the indices obtained from this routine

band_index_above = band_index_below+1
do
if ( eig(band_index_above)-omega > 1.0e-8 ) exit
band_index_above = band_index_above + 1
end do
! band_index_above is now the index of the lowest eigenvalue ABOVE eig(e)
band_index_above = band_index_above - 1
! band_index_above is now the highest index with eigenvalue equal to eig(e).
! This is the correct index to remove the states with energy <= eig(e)

end subroutine get_projection_band_indices
!!***

end module m_gwls_GWanalyticPart
!!***

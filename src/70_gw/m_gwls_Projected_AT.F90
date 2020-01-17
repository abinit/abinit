!!****m* ABINIT/m_gwls_Projected_AT
!! NAME
!! m_gwls_Projected_AT
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


module m_gwls_Projected_AT
!----------------------------------------------------------------------------------------------------
! This module contains routines to compute the matrix elements of the so-called A operator,
! which accounts for the Static correlation energy term.
!----------------------------------------------------------------------------------------------------
! local modules
use m_gwls_utility
use m_gwls_wf
use m_gwls_TimingLog
use m_gwls_hamiltonian
use m_gwls_lineqsolver
use m_gwls_GWlanczos
use m_gwls_LanczosBasis 
use m_gwls_LanczosResolvents

use m_gwls_GWanalyticPart, only : get_projection_band_indices
! abinit modules
use defs_basis
use m_abicore
use m_xmpi
use m_io_tools,  only : get_unit



implicit none
save
private
!!***

!!***

! Public methods
public :: compute_AT_shift_Lanczos
!!***

contains

!!****f* m_hamiltonian/compute_AT_shift_Lanczos
!! NAME
!!  compute_AT_shift_Lanczos
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
!!      cleanup_lanczosresolvents,compute_resolvent_column_shift_lanczos
!!      get_projection_band_indices,pc_k_valence_kernel,setup_lanczosresolvents
!!      wf_block_distribute,xmpi_sum
!!
!! SOURCE

subroutine compute_AT_shift_Lanczos(nfreq,list_external_omega,model_parameter,lmax, modified_Lbasis,kmax_analytic,list_AT_Lanczos)
!----------------------------------------------------------------------------------------------------
! This function returns the diagonal matrix elements of the so-called A^T operator, which is pertinent to the
! computation of the analytic energy term.
!
! The operator is given by
!
!                 Am(W) = w0/2 Um^dagger  . [ PW/(H-W-w0)+QW/(H-W+w0)] . Um
!
!  where W is the external frequency of the self energy, and w0 is the lorentzian parameter.
!  The operator PW projects on states of energy lower than W, and QW on states or energy higher than W.
!
!
!  For every l, We seek to compute  AT_l = < l |  A^T | l > = < l^* | A | l^* > 
!
!  It will be assumed that the array modified_Lbasis already contains the basis vectors U_m | l^* >.
!
!  This function does not use SQMR, but rather shift Lanczos to extract the values of the matrix elements for 
!  all external frequencies.
!----------------------------------------------------------------------------------------------------
implicit none

integer,      intent(in) :: nfreq
real(dp),     intent(in) :: list_external_omega(nfreq)
real(dp),     intent(in) :: model_parameter
integer,      intent(in) :: lmax
integer,      intent(in) :: kmax_analytic
complex(dpc), intent(in) :: modified_Lbasis(npw_k,lmax)
complex(dpc), intent(out):: list_AT_Lanczos(nfreq,lmax)


real(dp),     allocatable :: psik_wrk(:,:)
real(dp),     allocatable :: psikb_wrk(:,:)
real(dp),     allocatable :: psikg_wrk(:,:)

real(dp),     allocatable :: psikg(:,:)

complex(dpc), allocatable :: seed_vector(:)
complex(dpc), allocatable :: list_left_vectors(:,:)
complex(dpc), allocatable :: matrix_elements_resolvent(:,:)

complex(dpc), allocatable :: list_z_P(:), list_z_Q(:)
integer,      allocatable :: frequency_indices_array(:,:)


real(dp):: external_omega 
logical :: prec
integer :: nvec
integer :: band_index_below, band_index_above, bib0, bia0, bib, bia

integer :: l, lloc, iw_ext, iw_ext_min, iw_ext_max
integer :: ierr

integer :: number_of_frequency_blocks, ifreq_block

integer :: iblk, nbdblock_lanczos
integer :: mb 

integer :: nz

integer :: io_unit
integer :: mpi_band_rank

! *************************************************************************

mpi_band_rank    = mpi_enreg%me_band

!=================================================
!
! Find the number of frequency blocks, where
! the projection operators are constant within a
! block.
!=================================================

if (mpi_enreg%me == 0 ) then
  io_unit  = get_unit()
  open(io_unit,file='Frequency_blocks_AT.log',position='append')

  write(io_unit,10) "#================================================================================"
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  This log file documents how the algorithm in compute_AT_shift_Lanczos         "
  write(io_unit,10) "#  separates the external frequencies in blocks. It is quite easy to put         "
  write(io_unit,10) "#  bugs in this algorithm, so monitoring is critical.                            "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#================================================================================"


  write(io_unit,10) "                                                                                 "
  write(io_unit,10) "#================================================================================"
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  input parameters:                                                             "
  write(io_unit,10) "#                                                                                "
  write(io_unit,15) "#                      nfreq    : ",nfreq                                                     
  write(io_unit,17) "#          list_external_omega  : ",list_external_omega  
  write(io_unit,17) "#              model_parameter  : ",model_parameter
  write(io_unit,15) "#                       lmax    : ",lmax
  write(io_unit,15) "#                kmax_analytic  : ",kmax_analytic
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#================================================================================"
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  Building the blocks                                                           "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#         bib : band index below                                                 "
  write(io_unit,10) "#         bia : band index above                                                 "
  write(io_unit,10) "#         nfb : number_of_frequency_blocks                                       "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  iw_ext      external_omega (Ha)     bib  bia  nfb                             "
  write(io_unit,10) "#================================================================================"
end if
external_omega = list_external_omega(1)


! bib == band_index_below
! bia == band_index_above
call get_projection_band_indices(external_omega,bib0, bia0)

number_of_frequency_blocks = 1

! loop on all external frequencies
do iw_ext = 1 , nfreq

external_omega = list_external_omega(iw_ext)
! Define the energy of the state to be corrected

! Find the indices for the projections PW and QW
call get_projection_band_indices(external_omega,bib, bia)


if (mpi_enreg%me == 0 ) write(io_unit,20) iw_ext, external_omega, bib, bia, number_of_frequency_blocks 

if (bib /= bib0 .or. bia /= bia0) then

  if (mpi_enreg%me == 0 )   write(io_unit,10) "*************    new block!    **********************"
  bib0 = bib
  bia0 = bia

  number_of_frequency_blocks = number_of_frequency_blocks + 1

end if
end do

!=================================================
!
! fill the frequency indices array
!
!=================================================

if (mpi_enreg%me == 0 ) then
  write(io_unit,10) "                                                                                 "
  write(io_unit,10) "#================================================================================"
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  Filling the frequency indices array                                           "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  iw_ext      iw_ext_min     iw_ext_max  ifreq_block                            "
  write(io_unit,10) "#================================================================================"
end if


ABI_ALLOCATE(frequency_indices_array, (2,number_of_frequency_blocks))
frequency_indices_array = 0

iw_ext_min = 1
iw_ext_max = 1

! loop on all external frequencies
external_omega = list_external_omega(1)
call get_projection_band_indices(external_omega,bib0, bia0)
ifreq_block = 1

do iw_ext = 1 , nfreq

if (mpi_enreg%me == 0 ) write(io_unit,30) iw_ext, iw_ext_min, iw_ext, ifreq_block

external_omega = list_external_omega(iw_ext)
! Define the energy of the state to be corrected

! Find the indices for the projections PW and QW
call get_projection_band_indices(external_omega,bib, bia)

if (bib /= bib0 .or. bia /= bia0) then

  if (mpi_enreg%me == 0 )   write(io_unit,10) "*************    new block!    **********************"
  bib0 = bib
  bia0 = bia

  ! write previous block
  frequency_indices_array(1,ifreq_block) = iw_ext_min
  frequency_indices_array(2,ifreq_block) = iw_ext-1 ! we went one too far

  ifreq_block  = ifreq_block  + 1
  iw_ext_min = iw_ext
end if 

end do

! write last block!
frequency_indices_array(1,ifreq_block) = iw_ext_min
frequency_indices_array(2,ifreq_block) = nfreq

if (mpi_enreg%me == 0 ) then
  write(io_unit,10) "                                                                                 "
  write(io_unit,10) "#================================================================================"
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  Final frequency_indices_array :                                               "
  write(io_unit,10) "#                                                                                "
  write(io_unit,10) "#  ifreq_block iw_ext_min     iw_ext_max                                         "
  write(io_unit,10) "#================================================================================"

  do ifreq_block = 1 , number_of_frequency_blocks 

  write(io_unit,40) ifreq_block, frequency_indices_array(1,ifreq_block), frequency_indices_array(2,ifreq_block)

  end do
end if 




!=================================================
!
! Set up the LanczosResolvents algorithm
!
!=================================================

prec = .false.  ! let's not precondition for now
call setup_LanczosResolvents(kmax_analytic, prec)

!=================================================
!
! iterate on frequency blocks!
!
!=================================================
nvec = 1

ABI_ALLOCATE(list_left_vectors, (npw_g,nvec))
ABI_ALLOCATE(seed_vector, (npw_g))

ABI_ALLOCATE(psik_wrk,  (2,npw_k))
ABI_ALLOCATE(psikb_wrk, (2,npw_kb))

ABI_ALLOCATE(psikg_wrk, (2,npw_g))
ABI_ALLOCATE(psikg,     (2,npw_g))

list_AT_Lanczos(:,:) = cmplx_0

! Number of blocks of lanczos vectors
nbdblock_lanczos = lmax/blocksize
if (modulo(lmax,blocksize) /= 0) nbdblock_lanczos = nbdblock_lanczos + 1


do ifreq_block = 1, number_of_frequency_blocks

!-------------------------------
! Build the frequency arrays
!-------------------------------
iw_ext_min = frequency_indices_array(1,ifreq_block) 
iw_ext_max = frequency_indices_array(2,ifreq_block)


nz = iw_ext_max -iw_ext_min + 1

ABI_ALLOCATE( list_z_P, (nz))
ABI_ALLOCATE( list_z_Q, (nz))
ABI_ALLOCATE(matrix_elements_resolvent, (nz,nvec))


external_omega = list_external_omega(iw_ext_min)
call get_projection_band_indices(external_omega,band_index_below, band_index_above)


do iw_ext = iw_ext_min, iw_ext_max

list_z_P(iw_ext-iw_ext_min+1) = list_external_omega(iw_ext)+ model_parameter
list_z_Q(iw_ext-iw_ext_min+1) = list_external_omega(iw_ext)- model_parameter

end do

!-----------------------------------------
! Compute with shift Lanczos, using blocks
!-----------------------------------------
do iblk = 1, nbdblock_lanczos
! which l is on this row of FFT processors?
l = (iblk-1)*blocksize + mpi_band_rank + 1 

! Change the configuration of the data
do mb =1, blocksize
lloc = (iblk-1)*blocksize+mb
if (lloc <= lmax) then 
  psik_wrk(1,:)  = dble ( modified_Lbasis(:,lloc) )
  psik_wrk(2,:)  = dimag( modified_Lbasis(:,lloc) )
else
  psik_wrk(:,:)  = zero
end if

psikb_wrk(:,(mb-1)*npw_k+1:mb*npw_k) = psik_wrk(:,:) 

end do ! mb 

! change configuration of the data, from LA to FFT
call wf_block_distribute(psikb_wrk,  psikg_wrk, 1) ! LA -> FFT


if (band_index_below == 0 .and. l <= lmax) then
  ! There are no DFT states with energy below W!
  ! PW is thus zero, and QW = I.

  seed_vector(:)         = cmplx_1*psikg_wrk(1,:) + cmplx_i*psikg_wrk(2,:)
  list_left_vectors(:,1) = seed_vector(:)

  call compute_resolvent_column_shift_lanczos(nz, list_z_Q, nvec, list_left_vectors, seed_vector, &
  &                                                                   matrix_elements_resolvent)

  list_AT_Lanczos(iw_ext_min:iw_ext_max, l) = 0.5_dp*model_parameter*matrix_elements_resolvent(:,1)


else if ( l <= lmax ) then

  !----------------------------------------------------------------------------------
  !
  ! CAREFUL HERE! PW + QW != I. If the external frequency is exactly equal to a 
  !               DFT eigenvalue, then the projections PW and QW both project OUT
  !               of the subspace corresponding to this energy! 
  !
  !----------------------------------------------------------------------------------

  !-------------------
  ! Treat first PW
  !-------------------
  psikg(:,:) = psikg_wrk(:,:)
  call pc_k_valence_kernel(psikg,band_index_below)

  seed_vector = cmplx_1*psikg(1,:)+cmplx_i*psikg(2,:)

  list_left_vectors(:,1) = seed_vector(:)

  call compute_resolvent_column_shift_lanczos(nz, list_z_P, nvec, list_left_vectors, seed_vector, &
  &                                                                   matrix_elements_resolvent)

  list_AT_Lanczos(iw_ext_min:iw_ext_max, l) = 0.5_dp*model_parameter*matrix_elements_resolvent(:,1)


  !-------------------
  ! Treat second QW
  !-------------------

  psikg(:,:) = psikg_wrk(:,:)

  call pc_k_valence_kernel(psikg,band_index_above)

  psikg(:,:) = psikg_wrk(:,:)-psikg(:,:)

  seed_vector = cmplx_1*psikg(1,:)+cmplx_i*psikg(2,:)

  list_left_vectors(:,1) = seed_vector(:)

  call compute_resolvent_column_shift_lanczos(nz, list_z_Q, nvec, list_left_vectors, seed_vector, &
  &                                                                   matrix_elements_resolvent)

  list_AT_Lanczos(iw_ext_min:iw_ext_max, l) = list_AT_Lanczos(iw_ext_min:iw_ext_max, l) +   &
  0.5_dp*model_parameter*matrix_elements_resolvent(:,1)

end if

end do

ABI_DEALLOCATE( list_z_P)
ABI_DEALLOCATE( list_z_Q)
ABI_DEALLOCATE(matrix_elements_resolvent)

end do


! sum all results
call xmpi_sum(list_AT_Lanczos, mpi_enreg%comm_band,ierr) ! sum on all processors working on bands


if (mpi_enreg%me == 0 ) then
  close(io_unit)
end if 

ABI_DEALLOCATE(list_left_vectors)
ABI_DEALLOCATE(seed_vector)
ABI_DEALLOCATE(frequency_indices_array)


ABI_DEALLOCATE(psik_wrk)
ABI_DEALLOCATE(psikb_wrk)

ABI_DEALLOCATE(psikg_wrk)
ABI_DEALLOCATE(psikg)




call cleanup_LanczosResolvents

10 format(A)
15 format(A,I5)
17 format(A,1000F8.4)
20 format(I5,7X,ES24.16,3I5)
30 format(4(I5,10X))
40 format(3(I5,10x))

end subroutine compute_AT_shift_Lanczos
!!***


end module m_gwls_Projected_AT
!!***

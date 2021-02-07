!!****m* ABINIT/m_gwls_valenceWavefunctions
!! NAME
!! m_gwls_valenceWavefunctions
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


module m_gwls_valenceWavefunctions

! local modules
use m_gwls_utility
use m_gwls_hamiltonian

! abinit modules
use defs_basis
use m_abicore
use m_xmpi
implicit none
save
private
!!***

real(dp), public, allocatable :: valence_wfr(:,:,:,:,:)
real(dp), public, allocatable :: valence_wfr_fftpac(:,:,:)
!!***

public :: prepareValenceWavefunctions
public :: cleanupValenceWavefunctions
public :: compute_Exchange_and_Correlation_energies
!!***

contains

!!****f* m_hamiltonian/prepareValenceWavefunctions
!! NAME
!!  prepareValenceWavefunctions
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

subroutine prepareValenceWavefunctions()
!--------------------------------------------------------------------------------
!
! This subroutine allocates and fills an array valence_wfr which will contain
! the valence wavefunctions in real space. We do this once and for all, instead
! of constantly performing FFTs throughout the program.
!
! NOTE: fftpac is not reversible while using MPI!  At first, we tried storing
!       valence wavefunctions in real space using fftpac, but the unpacking routine
!       doesn't return the original (n4,n5,n6) array !
!
!      Valence wavefunctions will be stored in Fourier space instead.
!
!--------------------------------------------------------------------------------
implicit none

!integer  :: v, kmin, kmax
!
!
!real(dp), allocatable :: psir(:,:,:,:)

! *************************************************************************


! Routine is left blank for now


! old code, broken because fftpac isn't invertible !
! ABI_MALLOC(valence_wfr_fftpac,(2,nfft,nbandv)) ! nfft is a public variable from the gwls_hamiltonian module

! ABI_MALLOC(psir,(2,n4,n5,n6))
! psir = zero
! do v=1,nbandv
!        kmin = 1+(v-1)*npw_k
!        kmax =    v   *npw_k
!
!       ! transform from k to r, storing wavefunction in real space in work array psir1
!       call g_to_r(psir,cg(:,kmin:kmax))
!
!       ! pack the real-space wavefunction in the purpose-made array
!       call sg_to_dg(valence_wfr_fftpac(:,:,v), psir)
! end do
! ABI_FREE(psir)

! old code!
!ABI_MALLOC(valence_wfr,(2,n4,n5,n6,nbandv))
!valence_wfr = zero

!do v=1,nbandv
!       kmin = 1+(v-1)*npw_k
!       kmax =    v   *npw_k
!       call g_to_r(valence_wfr(:,:,:,:,v),cg(:,kmin:kmax))
!end do

end subroutine prepareValenceWavefunctions
!!***

!!****f* m_hamiltonian/cleanupValenceWavefunctions
!! NAME
!!  cleanupValenceWavefunctions
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

subroutine cleanupValenceWavefunctions()
!--------------------------------------------------------------------------------
!
! This subroutine deallocates the array valence_wfr once it is no longer needed.
!
!--------------------------------------------------------------------------------
implicit none

! *************************************************************************

!if (allocated(valence_wfr)) ABI_FREE(valence_wfr)
!if (allocated(valence_wfr_fftpac)) ABI_FREE(valence_wfr_fftpac)

end subroutine cleanupValenceWavefunctions
!!***

!!****f* m_hamiltonian/compute_Exchange_and_Correlation_energies
!! NAME
!!  compute_Exchange_and_Correlation_energies
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

subroutine compute_Exchange_and_Correlation_energies(e_index, exchange_energy, Vxc_energy)

!--------------------------------------------------------------------------------
!
! This subroutine computes the exchange and correlation energies.
!
!--------------------------------------------------------------------------------
implicit none

integer, intent(in)    :: e_index
real(dp), intent(out)  :: exchange_energy
real(dp), intent(out)  :: vxc_energy

! *************************************************************************

vxc_energy         = dft_xc_energy(e_index)

exchange_energy    = exchange(e_index)

end subroutine compute_Exchange_and_Correlation_energies
!!***

end module m_gwls_valenceWavefunctions
!!***

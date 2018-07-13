!{\src2tex{textfont=tt}}
!!****f* ABINIT/diag_occ
!! NAME
!! diag_occ
!!
!! FUNCTION
!! Use for DMFT in KGB parallelisation. Diagonalise the occupation matrix 
!! and return diagonalised occupations and associated eigen vectors sorted
!! with descending occupation
!! 
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (TCavignac)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   occ_nd_cpx(nband, nband) = matrix of non diagonal occupations for DMFT
!!   nband = number of band to be processed
!!   
!!
!! OUTPUT
!!   occ_diag(2, nband) = diagonal occupation in the new band space 
!!   cwavef_rot(2, npw, nband) = fourier coefficient of wave functions for all bands rotated in the band space
!!
!! PARENTS
!!      mkrho
!!
!! CHILDREN
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !
!! TODO add the possibility of using ScaLAPACK to do computation in parallel

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine diag_occ(occ_nd_cpx, nband, occ_diag) 

  use defs_basis
  use defs_abitypes
  use m_profiling_abi
  use m_xmpi
  use m_errors

! use m_paw_dmft,     only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diag_occ'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: nband
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: paw_dmft
!no_abirules
  complex(kind=dpc), intent(inout) :: occ_nd_cpx(nband, nband)
  real(kind=dp), intent(out) :: occ_diag(nband)

!Local variables-------------------------------

!scalars
  integer :: info, lwork
  character(len=500) :: message

!arrays
  real(kind=dp) :: rwork(3*nband-1)
  complex(kind=dpc), allocatable :: work(:)

! *************************************************************************

  DBG_ENTER("COLL")

!! Initialisation
  rwork = zero
  info = 0

!! use the opposite to have zheev orders the eigenvalues descending (once the
!! opposite have been taken again)
  occ_nd_cpx = -occ_nd_cpx

!! Get diagonal occupations and associeted base

! Compute the optimal working array size
  ABI_ALLOCATE(work,(1))
  work = czero
  call zheev('V', 'U', nband, occ_nd_cpx, nband, occ_diag, work, -1, rwork, info)
  lwork = work(1)
  ABI_DEALLOCATE(work)

! Compute the eigenvalues (occ_diag) and vectors
  ABI_ALLOCATE(work,(lwork))
  work = czero

  call zheev('V', 'U', nband, occ_nd_cpx, nband, occ_diag, work, lwork, rwork, info)

!! obtain the true eigen values of occupation matrix in descending order
  occ_diag = -occ_diag

  ABI_DEALLOCATE(work)

  if (info > 0) then
    message=""
    write(message, "(a,i5)") " something wrong happened with the diagonalisation of the occupation matrix (did't converge), info=",info
    MSG_ERROR(message)
  else if (info < 0) then
    message=""
    write(message, "(a,i5)") " something wrong happened with the diagonalisation of the occupation matrix (bad input argument), info=",info
    MSG_ERROR(message)
  end if
  
  DBG_EXIT("COLL")

end subroutine diag_occ
!!***

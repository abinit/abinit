!{\src2tex{textfont=tt}}
!!****f* ABINIT/diag_occ
!! NAME
!! m_rot_cg
!!
!! FUNCTION
!!  Rotate the cg coefficient with the rotation matrix obtained from the
!!  diagonalisation of the non-diagonal occupation matrix produced by DMFT.
!! 
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (TCavignac)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mkrho
!!
!! CHILDREN
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rot_cg

  use defs_basis
  use defs_abitypes
  use m_profiling_abi
  use m_xmpi
  use m_errors

  implicit none

  private

  public :: rot_cg

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
!!
!! PARENTS
!!      rot_cg
!!
!! CHILDREN
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !
subroutine diag_occ(occ_nd_cpx, nband, occ_diag) 

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
  complex(kind=dpc) :: rwork(3*nband-1)
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
!{\src2tex{textfont=tt}}
!!****f* ABINIT/rot_cg
!! NAME
!! rot_cg
!!
!! FUNCTION
!! Use for DMFT in KGB parallelisation. Diagonalise the occupation matrix 
!! and use the resulting base to represent the wave functions.
!! 
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (TCavignac)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   occ_nd(2, nband, nband) = matrix of non diagonal occupations for DMFT
!!   cwavef(2, npw, nband) = fourier coefficient of wave functions for all bands
!!   npw = number of G vectors computed in this iteration
!!   nband = number of band to be processed
!!   blocksize = size of the block for th LO.. algorithm
!!               still have to be equal to nband
!!   nspinor = number of spinor components
!!   
!!
!! OUTPUT
!!   occ_diag(2, nband) = diagonal occupation in the new band space 
!!
!! SIDE EFFECT
!!   cwavef is rotated with the unitary matrix obtained from the diagonalisation
!!   of occupations (occ_nd)
!!
!! PARENTS
!!      mkrho 
!!
!! CHILDREN
!!      diag_occ
!!
!! SOURCE
!!
!! TODO /!\ No parallel computing yet !

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rot_cg(occ_nd, cwavef, npw, nband, blocksize, nspinor, occ_diag)

  use defs_basis
  use defs_abitypes
  use m_profiling_abi
  use m_xmpi
  use m_errors

! use m_paw_dmft,     only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rot_cg'
 use interfaces_67_common, except_this_one => rot_cg
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: npw, nband, blocksize, nspinor
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: paw_dmft
!no_abirules
  real(kind=dp), intent(in) :: occ_nd(2, blocksize, blocksize)
  real(kind=dp), intent(inout) :: cwavef(2, npw, blocksize, nspinor)
  real(kind=dp), intent(out) :: occ_diag(blocksize)

!Local variables-------------------------------

!scalars
  integer :: n, np, ig
  character(len=500) :: message

!arrays
  complex(kind=dpc) :: occ_nd_cpx(blocksize, blocksize)
  complex(kind=dpc) :: cwavef_rot_g(blocksize, nspinor)

! *************************************************************************

  DBG_ENTER("COLL")

  if(nband /= blocksize) then
    message = " DMFT in KGB cannot be used with multiple blocks yet. Make sure that bandpp*npband = nband."
    MSG_ERROR(message)
  end if

!! Initialisation

  do n=1,blocksize
    do np=1,blocksize
      occ_nd_cpx(n,np) = cmplx(occ_nd(1,n,np), occ_nd(2,n,np), kind=dpc)
    end do
  end do

!! Get diagonal occupations and associeted base

  call diag_occ(occ_nd_cpx, blocksize, occ_diag)

!! Compute the corresponding wave functions if nothing wrong happened
  ! $c^{rot}_{n,k}(g) =  \sum_{n'} [\bar{f_{n',n}} * c_{n',k}(g)]$
  do ig=1,npw
    cwavef_rot_g(:,:) = czero
    do n=1,blocksize
      do np=1,blocksize
        cwavef_rot_g(n,:) = cwavef_rot_g(n,:) + dconjg(occ_nd_cpx(np, n)) * &
&                           cmplx(cwavef(1,ig,np,:), cwavef(2,ig,np,:), kind=dpc)
      end do
    end do
    cwavef(1,ig,:,:) = dreal(cwavef_rot_g)
    cwavef(2,ig,:,:) = dimag(cwavef_rot_g)
  end do
  
  DBG_EXIT("COLL")

end subroutine rot_cg
!!***

end module m_rot_cg

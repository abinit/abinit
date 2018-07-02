!{\src2tex{textfont=tt}}
!!****f* ABINIT/diag_occ_rot_cg
!! NAME
!! diag_occ_rot_cg
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
!!   cwavef_toberot(2, npw, nband) = fourier coefficient of wave functions for all bands
!!   npw = number of G vectors computed in this iteration
!!   nband = number of band to be processed
!!   blocksize = size of the block for th LO.. algorithm
!!               still have to be equal to nband
!!   nspinor = number of spinor components
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
!! TODO Make the computation of the new wf parallel

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine diag_occ_rot_cg(occ_nd, cwavef_toberot, npw, nband, blocksize, nspinor, occ_diag, &
&                          cwavef_rot) 

  use defs_basis
  use defs_abitypes
  use m_profiling_abi
  use m_xmpi
  use m_errors

! use m_paw_dmft,     only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diag_occ_rot_cg'
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
  real(kind=dp), intent(in) :: cwavef_toberot(2, npw, blocksize, nspinor)
  real(kind=dp), intent(out) :: occ_diag(blocksize)
  real(kind=dp), intent(out) :: cwavef_rot(2, npw, blocksize, nspinor)

!Local variables-------------------------------

!scalars
  integer :: info, lwork, n, np
  character(len=500) :: message

!arrays
  complex(kind=dpc) :: occ_nd_cpx(blocksize, blocksize), rwork(3*blocksize-1), fc_tmp(npw,nspinor)
  complex(kind=dpc), allocatable :: work(:)

! *************************************************************************

  DBG_ENTER("COLL")

  if(nband /= blocksize) then
    message = " DMFT in KGB cannot be used with multiple blocks yet. Make sure that bandpp*npband = nband."
    MSG_ERROR(message)
  end if

!! Initialisation
  rwork = zero
  info = 0
  cwavef_rot(:,:,:,:) = zero

  do n=1,blocksize
    do np=1,blocksize
      occ_nd_cpx(n,np) = -cmplx(occ_nd(1,n,np), occ_nd(2,n,np), kind=dpc)
    end do
  end do

!! Get diagonal occupations and associeted base

! Compute the optimal working array size
  ABI_ALLOCATE(work,(1))
  work = czero
  call zheev('V', 'U', blocksize, occ_nd_cpx, blocksize, occ_diag, work, -1, rwork, info)
  lwork = work(1)
  ABI_DEALLOCATE(work)

! Compute the eigenvalues (occ_diag) and vectors
  ABI_ALLOCATE(work,(lwork))
  work = czero

  call zheev('V', 'U', blocksize, occ_nd_cpx, blocksize, occ_diag, work, lwork, rwork, info)
! occ_nd_cpx is now the eigen vectors (P) of the occ_nd matrix

  occ_diag = -occ_diag
  write (*,*) "))) Valeures propres des occupations ", occ_diag
  ABI_DEALLOCATE(work)

!! Compute the corresponding wave functions if nothing wrong happened
  if(info == 0) then
    ! $c^{rot}_{n,k}(g) =  \sum_{n'} [\bar{f_{n',n}} * c_{n',k}(g)]$
    fc_tmp = zero
    do n=1,blocksize
      do np=1,blocksize
        fc_tmp = dconjg(occ_nd_cpx(np, n)) * cmplx(cwavef_toberot(1,:,np,:), cwavef_toberot(2,:,np,:), kind=dpc)
        cwavef_rot(1,:,n,:) = cwavef_rot(1,:,n,:) + dreal(fc_tmp)
        cwavef_rot(2,:,n,:) = cwavef_rot(2,:,n,:) + dimag(fc_tmp)
!        cwavef_rot(1,:,n,:) = cwavef_rot(1,:,n,:) + dreal(occ_nd_cpx(np, n)) * cwavef_toberot(1,:,np,:) &
!&                                                 + dimag(occ_nd_cpx(np, n)) * cwavef_toberot(2,:,np,:)
!        cwavef_rot(2,:,n,:) = cwavef_rot(2,:,n,:) + dreal(occ_nd_cpx(np, n)) * cwavef_toberot(2,:,np,:) &
!&                                                 - dimag(occ_nd_cpx(np, n)) * cwavef_toberot(1,:,np,:)
      end do
    end do
  else if (info > 0) then
    message=""
    write(message, "(a,i5)") " something wrong happened with the diagonalisation of the occupation matrix (did't converge), info=",info
    MSG_ERROR(message)
  else
    message=""
    write(message, "(a,i5)") " something wrong happened with the diagonalisation of the occupation matrix (bad input argument), info=",info
    MSG_ERROR(message)
  end if
  
  DBG_EXIT("COLL")

end subroutine diag_occ_rot_cg
!!***

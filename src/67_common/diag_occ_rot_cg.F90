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
!!   spinpol = number of spin polarisation
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine diag_occ_rot_cg(occ_nd, cwavef_toberot,npw, nband, blocksize, spinpol, occ_diag, &
&                          cwavef_rot) 

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_paw_dmft,     only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diag_occ_rot_cg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw, nband, blocksize, spinpol
!! type(MPI_type),intent(inout) :: mpi_enreg
!! type(dataset_type),intent(in) :: dtset
!! type(paw_dmft_type), intent(in)  :: paw_dmft
!no_abirules
 real(dp), intent(in) :: occ_nd(2, nband, nband)
 real(dp), intent(in) :: cwavef_toberot(2, npw, nband)
 real(dp), intent(out) :: occ_diag(2, nband)
 real(dp), intent(out) :: cwavef_rot(2, npw, nband)

!Local variables-------------------------------

!scalars

!arrays

! *************************************************************************

 DBG_ENTER("COLL")


 DBG_EXIT("COLL")

end subroutine diag_occ_rot_cg
!!***

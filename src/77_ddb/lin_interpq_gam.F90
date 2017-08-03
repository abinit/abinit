!{\src2tex{textfont=tt}}
!!****f* ABINIT/lin_interpq_gam
!!
!! NAME
!! lin_interpq_gam
!!
!! FUNCTION
!! This routine interpolates the electron phonon coupling matrix
!! to a give _qpoint_ by linear methods
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   nbranch=Number of phonon branches 3*natom
!!   nqbz=Number of q-points in full BZ.
!!   gamma_qpt(2,nbranch**2,nsppol,nqbz)=Gamma matrix in full BZ.
!!   nsppol=Number of independent spin polarizations.
!!   qpt = q-point in reciprocal space for interpolation
!!
!! OUTPUT
!!   gam_now = interpolated gamma matrix at qpt
!!
!! PARENTS
!!
!! CHILDREN
!!      wrap2_zero_one
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine lin_interpq_gam(gamma_qpt,nbranch,nqbz,nsppol,gam_now,isppol,kptrlatt,qpt)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,  only : wrap2_zero_one, interpol3d

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lin_interpq_gam'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nbranch,isppol,nsppol,nqbz
!arrays
 integer, intent(in) :: kptrlatt(3,3)
 real(dp), intent(in) :: qpt(3)
 real(dp),intent(in) :: gamma_qpt(2,nbranch**2,nsppol,nqbz)
 real(dp), intent(out) :: gam_now(2,nbranch**2)

!Local variables-------------------------------
!scalars
 integer :: ibranch
 real(dp) :: res
 character(len=500) :: msg
!arrays
 real(dp) :: qpt_incube(3)

! *************************************************************************

!the array gamma_qpt has dimensions: (2,nbranch**2,nsppol,elph_ds%k_fine%nkpt)
 if (kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3) /= nqbz) then
   write(msg,'(a,2i0)')"Wrong dimensions in gamma_qpt ",kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3),nqbz
   MSG_ERROR(msg)
 end if

 call wrap2_zero_one(qpt(1),qpt_incube(1),res)
 call wrap2_zero_one(qpt(2),qpt_incube(2),res)
 call wrap2_zero_one(qpt(3),qpt_incube(3),res)
 
 do ibranch=1,nbranch**2
   gam_now(1,ibranch) = interpol3d(qpt_incube,kptrlatt(1,1),kptrlatt(2,2),kptrlatt(3,3),gamma_qpt(1,ibranch,isppol,:))
   gam_now(2,ibranch) = interpol3d(qpt_incube,kptrlatt(1,1),kptrlatt(2,2),kptrlatt(3,3),gamma_qpt(2,ibranch,isppol,:))
 end do

end subroutine lin_interpq_gam
!!***

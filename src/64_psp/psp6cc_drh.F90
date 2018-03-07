!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp6cc_drh
!! NAME
!! psp6cc_drh
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 6 of pseudopotentials.
!! Version modified by DHamann, with consistent treatment
!! of the derivatives in this routine and the remaining of the code.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (AF,DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! PARENTS
!!      psp6in
!!
!! CHILDREN
!!      cc_derivatives
!!
!! NOTES
!! Test version by DRH - requires very smooth model core charge
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp6cc_drh'
 use interfaces_64_psp, except_this_one => psp6cc_drh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: irad
 character(len=500) :: errmsg
!arrays
 real(dp),allocatable :: ff(:),ff1(:),ff2(:),rad(:)

!**********************************************************************

 ABI_ALLOCATE(ff,(mmax))
 ABI_ALLOCATE(ff1,(mmax))
 ABI_ALLOCATE(ff2,(mmax))
 ABI_ALLOCATE(rad,(mmax))

!
!read from pp file the model core charge (ff) and first (ff1) and
!second (ff2) derivative on logarithmic mesh mmax; rad is the radial grid
!the input functions contain the 4pi factor, it must be rescaled.

!***drh test
 write(std_out,'(a,2i6)') 'drh:psp6cc_drh - mmax,n1xccc',mmax,n1xccc
!***end drh test
 do irad=1,mmax
   read(tmp_unit,*,err=10,iomsg=errmsg) rad(irad),ff(irad),ff1(irad),ff2(irad)
   ff(irad)=ff(irad)/4.d0/pi
   ff1(irad)=ff1(irad)/4.d0/pi
   ff2(irad)=ff2(irad)/4.d0/pi
 end do
 rad(1)=0.d0

 call cc_derivatives(rad,ff,ff1,ff2,mmax,n1xccc,rchrg,xccc1d)
 
 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(ff1)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(rad)

 return

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine psp6cc_drh
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/radsmear
!! NAME
!! radsmear
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2012-2018 ABINIT group (ILuk)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
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

function radsmear(r, rsph, rsm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'radsmear'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: radsmear
 real(dp), intent(in) :: r, rsph, rsm

!Local variables ------------------------------
!scalars
 real(dp) :: xx

!******************************************************************

 radsmear = zero
 if (r < rsph - rsm - tol12) then
   radsmear = one
 else if (r < rsph - tol12) then
   xx = (rsph - r) / rsm
   radsmear = xx**2*(3+xx*(1+xx*(-6+3*xx)))
 end if

end function radsmear
!!***

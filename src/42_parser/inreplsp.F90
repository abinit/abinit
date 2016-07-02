!{\src2tex{textfont=tt}}
!!****f* ABINIT/inreplsp
!! NAME
!! inreplsp
!!
!! FUNCTION
!! Replace all occurrences of characters lexically less than SP (blank)
!! by SP in the input string, returning modified string of same length.
!! Also replace a '=' by a SP.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string=character string to be modified
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  string=same character string with ASCII (decimal) 0-31 replaced by 32.
!!
!! PARENTS
!!      incomprs
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inreplsp(string)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inreplsp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
!scalars
 integer :: ilenth,length

! *************************************************************************

!Get length of string
 length=len(string)

!Proceed only if string has nonzero length
 if (length>0) then

!  Do replacement by going through input
!  character string one character at a time
   do ilenth=1,length
     if (llt(string(ilenth:ilenth),' ')) then
       string(ilenth:ilenth)=' '
     end if
     if(string(ilenth:ilenth)=='=')then
       string(ilenth:ilenth)=' '
     end if
   end do

 end if

end subroutine inreplsp
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/consist
!! NAME
!! consist
!!
!! FUNCTION
!! Checking of the consistency between the values of input variables
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  aim_dtset= the structured entity containing all input variables
!!  tstngr= information about the test on the ngrid input variable
!!  tstvpt= information about the test on the vpts input variable
!!
!! OUTPUT
!!  (only checking : print error message and stop if there is a problem)
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      adini
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine consist(aim_dtset,tstngr,tstvpt)

 use defs_basis
 use defs_aimprom
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'consist'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tstngr,tstvpt
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------

! *********************************************************************

!write(std_out,*) tstngr, tstvpt

 if (((aim_dtset%denout/=0).or.(aim_dtset%lapout/=0)).and.((tstngr < 1).or.(tstvpt < 2))) then
   MSG_ERROR('in input1 - I cannot do the output !')
 end if
 if ((aim_dtset%denout > 0).and.(aim_dtset%lapout>0)) then
   if (aim_dtset%denout/=aim_dtset%lapout) then
     write(std_out,*) 'ERROR in input - when both denout and lapout are positive non-zero,'
     write(std_out,*) 'they must be equal.'
     MSG_ERROR("Aborting now")
   end if
   if ((tstvpt < aim_dtset%denout+1).or.(tstngr < aim_dtset%denout)) then
     write(std_out,*) 'ERROR in input2 - I cannot do the output !'
     MSG_ERROR("Aborting now")
   end if
 elseif (aim_dtset%denout > 0) then
   if ((tstvpt < aim_dtset%denout+1).or.(tstngr < aim_dtset%denout)) then
     write(std_out,*) 'ERROR in input - I cannot do the output !'
     MSG_ERROR("Aborting now")
   end if
 elseif (aim_dtset%lapout > 0) then
   if ((tstvpt < aim_dtset%lapout+1).or.(tstngr < aim_dtset%lapout)) then
     write(std_out,*) 'ERROR in input - I cannot do the output !'
     MSG_ERROR("Aborting now")
   end if
 end if

 if ((aim_dtset%isurf==1).and.(aim_dtset%crit==0)) then
   write(std_out,*) 'ERROR in input - must have crit/=0 for isurf==1'
   MSG_ERROR("Aborting now")
 end if

 if (((aim_dtset%ivol/=0).or.(aim_dtset%irho/=0)).and.(aim_dtset%isurf==0)) then
   MSG_ERROR('in input - I cannot integrate without surface !')
 end if

end subroutine consist
!!***

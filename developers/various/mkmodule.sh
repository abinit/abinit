# Copyright (C) 2010-2020 ABINIT group (Yann Pouillon)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
# 
# The purpose of this script is to create a new module
# ${module_name}.F90, from the embedded template, 
# where "${module_name}" is the argument of the script.
# Suppose that one is in a source directory, and that
# Utility/template.F90 is accessible as ../Utilities/template.F90

if [ $# -lt 1 ]
then
	echo "Usage: mkmodule module_name [routine_name]"
	exit 1
fi

module_name="${1}"
routine_name="${2}"
test "${routine_name}" = "" && routine_name="dummy_routine"
this_year=`date '+%Y'`

echo -n "mkmodule: creating ${module_name}.F90..."

cat > ${module_name}.F90 <<EOF
!{\src2tex{textfont=tt}}
!!****m* ABINIT/${module_name}
!! NAME
!!  ${module_name}
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) ${this_year} ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module ${module_name}
    
 use defs_basis
 use m_profiling_abi
 use m_errors

 implicit none

 private

! *************************************************************************

contains 
!!***

!!****f* ABINIT/${module_name}/${routine_name}
!! NAME
!!  ${routine_name}
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ${routine_name}(argin,argout,option,sizein,sizeout)
    
 use defs_basis

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: option,sizein,sizeout
 integer , intent(in)  :: argin(sizein)
 integer , intent(out) :: argout(sizeout)
 real(dp), intent(out) ::                        

!Local variables-------------------------------
 integer ::                                      
 real(dp) ::                                     
!character(len=500) :: msg                   
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!  write(msg,'(3a,i0)')&
!&  'The argument option should be 1 or 2,',ch10,&
!&  'however, option=',option
!  MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!  write(msg,'(3a,i0)')&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!  MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

end subroutine ${routine_name}
!!***

end module ${module_name}
!!***
EOF

echo "done."


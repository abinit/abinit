!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddb_chkname
!! NAME ddb_chkname
!! ddb_chkname
!!
!!
!! FUNCTION
!! This small subroutine check the identity of its argument,
!! who are a6 names, and eventually send a message and stop
!! if they are found unequal
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nmfond= name which has to be checked
!! nmxpct= name expected for nmfond
!! nmxpct2= eventual second optional name (backward compatibility)
!!
!! OUTPUT
!!
!! TODO
!! Describe the inputs
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ddb_chkname(nmfond,nmxpct,nmxpct2)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_chkname'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: nmfond,nmxpct
 character(len=*),intent(in),optional :: nmxpct2

!Local variables-------------------------------
!scalars
 logical :: found
 character(len=500) :: nmfond_,nmxpct_,nmxpct2_
 character(len=500) :: message

! *********************************************************************

 nmxpct_ = trim(adjustl(nmxpct))
 nmfond_ = trim(adjustl(nmfond))

 found = (nmxpct_ == nmfond_)

 if (present(nmxpct2) .and. .not. found) then
   nmxpct2_ = trim(adjustl(nmxpct2))
   found = (nmxpct2_==nmfond_)
 end if

 if (.not. found) then
   write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&   '  Reading DDB, expected name was "',trim(nmxpct_),'"',ch10,&
&   '               and name found is "',trim(nmfond_),'"',ch10,&
&   '  Likely your DDB is incorrect.',ch10,&
&   '  Action : correct your DDB, or contact the ABINIT group.'
   MSG_ERROR(message)
 end if

end subroutine ddb_chkname
!!***

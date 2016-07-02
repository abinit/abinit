!{\src2tex{textfont=tt}}
!!****f* ABINIT/time_accu
!! NAME
!!  time_accu
!!
!! FUNCTION
!!  Timing subroutine.  Calls machine-dependent "timein" which
!!  returns elapsed cpu and wall clock times in sec.
!!  Also return the number of times the counter has been called
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed);
!!
!! OUTPUT
!!  tottim(2)=accumulated time for accumulator nn
!!  totftimes(2)=accumulated time for accumulator nn evaluated by papi
!!  totffops =accumulated number of flops for accumulator nn evaluated by papi
!!  return_ncount gives the number of times that the accumulator has been incremented
!!
!! PARENTS
!!      timana
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine time_accu(nn,return_ncount,tottim,totflops,totftimes)

 use defs_basis
 use defs_time
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_accu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 integer,intent(out) :: return_ncount
 real(dp),intent(out) :: totflops
!arrays
 real(dp),intent(out) :: totftimes(2),tottim(2)

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

!Check that nn lies in sensible bounds
 if (nn<0.or.nn>mtim) then
   write(message,'(a,i6,a,i8,a)')' dim mtim=',mtim,' but input nn=',nn,'.'
   MSG_BUG(message)
 end if

!return accumulated time for nn
 tottim(1)=acctim(1,nn)
 tottim(2)=acctim(2,nn)

!return accumulated number flops for nn
 totflops = papi_accflops(nn) 

!return accumulated time for nn evaluated by papi
 totftimes(1) = papi_acctim(1,nn) 
 totftimes(2) = papi_acctim(2,nn) 
 return_ncount=ncount(nn)

end subroutine time_accu
!!***

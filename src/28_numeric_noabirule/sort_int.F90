!{\src2tex{textfont=tt}}
!!****f* ABINIT/sort_int
!! NAME
!!  sort_int
!!
!! FUNCTION
!!   Sort integer array list(n) into ascending numerical order using Heapsort
!!   algorithm, while making corresponding rearrangement of the integer
!!   array iperm. 
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2016 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  n        intent(in)    dimension of the list
!!  list(n)  intent(inout) list of double precision numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation given the right ascending order
!!
!! PARENTS
!!      getng,getngrec,initmpi_img,invars2,irrzg,m_dvdb,m_hdr,m_nesting,m_wfk
!!      mkfskgrid,shellstruct
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sort_int(n,list,iperm)

 use defs_basis, only : std_out
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sort_int'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n 
 integer,intent(inout) :: list(n),iperm(n)

!Local variables-------------------------------
!scalars
 integer :: l,ir,i,j,ip,ipp
! *************************************************************************

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(std_out,1000) n
  1000  format(/,' sort_int has been called with array length n=',i12,/, &
&  ' having a value less than 1.  This is not allowed.')
  MSG_ERROR("Aborting now")

 else ! n>1

! Conduct the usual sort

  l=n/2+1
  ir=n

  do   ! Infinite do-loop
 
   if (l>1) then

    l=l-1
    ip=list(l)
    ipp=iperm(l)

   else

    ip=list(ir)
    ipp=iperm(ir)
    list(ir)=list(1)
    iperm(ir)=iperm(1)
    ir=ir-1

    if (ir==1) then
     list(1)=ip
     iperm(1)=ipp
     exit   ! This is the end of this algorithm
    end if

   end if ! l>1

   i=l
   j=l+l

   do while (j<=ir)
    if (j<ir) then
     if (list(j).lt.list(j+1)) j=j+1
    end if
    if (ip.lt.list(j)) then
     list(i)=list(j)
     iperm(i)=iperm(j)
     i=j
     j=j+j
    else
     j=ir+1
    end if
   enddo

   list(i)=ip
   iperm(i)=ipp

  enddo ! End infinite do-loop

 end if ! n>1

end subroutine sort_int
!!***

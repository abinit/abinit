!!****m* ABINIT/m_pretty_rec
!! NAME
!!  m_pretty_rec
!!
!! FUNCTION
!!  This module  provides some minor functions applied in the
!!  recursion
!!
!! COPYRIGHT
!! Copyright (C) 2002-2020 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


MODULE m_pretty_rec

 use defs_basis
 use m_abicore

 implicit none

 private

 public ::  prtwork                ! calculate the work done during recursion

 interface prtwork
  module procedure prtworksim
  module procedure prtworksiRe
  module procedure prtworkadv
 end interface prtwork

CONTAINS  !===========================================================
!!***

!!****f* m_pretty_rec/prtworksim
!! NAME
!! prtworksim
!!
!! FUNCTION
!! Calculate the work done during recursion
!!
!! INPUTS
!! counter
!!
!! SIDE EFFECTS
!!
!! OUTPUT
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prtworksim(work_now)

 implicit none
!Arguments ------------------------------------
! scalars
 integer,intent(in) :: work_now
! arrays
!Local ---------------------------
 integer,save :: work_done=0
 integer :: ii,pt
 character(500) :: msg
! *********************************************************************
 if(mod(work_done-work_now,5)==0 .and. work_now/=work_done .and. work_now <105) then
  work_done = work_now
  pt = work_done/5
  write(msg,'(a,23a,i3,a)')'work done:',('-',ii=1,pt),(' ',ii=1,23-pt),work_done,'%'
  call wrtout(std_out,msg,'COLL')
 endif

end subroutine prtworksim
!!***

!!****f* m_pretty_rec/prtworksiRe
!! NAME
!! prtworksiRe
!!
!! FUNCTION
!! Calculate the work done during recursion
!!
!! INPUTS
!! counter
!!
!! SIDE EFFECTS
!!
!! OUTPUT
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prtworksiRe(work_now)

 implicit none
!Arguments ------------------------------------
! scalars
 real(dp),intent(in) :: work_now
! arrays
!Local ---------------------------
 integer,save :: work_done=0
 integer :: ii,pt,work_now_int
 character(500) :: msg
! *********************************************************************
 work_now_int = int(work_now)
 if(mod(work_done-work_now_int,5)==0 .and. work_now_int/=work_done .and. work_now_int <105) then
  work_done = work_now_int
  pt = work_done/5
  write(msg,'(a,23a,i3,a)')'work done:',('-',ii=1,pt),(' ',ii=1,23-pt),work_done,'%'
  call wrtout(std_out,msg,'COLL')
 endif

end subroutine prtworksiRe
!!***

!!****f* m_pretty_rec/prtworkadv
!! NAME
!! prtworkadv
!!
!! FUNCTION
!! Calculate the work done during recursion
!!
!! INPUTS
!! counter
!!
!! SIDE EFFECTS
!!
!! OUTPUT
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prtworkadv(work_now,adv)

 implicit none

!Arguments ------------------------------------
! scalars
 integer,intent(in) :: work_now
 character(1),intent(in) ::adv
! arrays
!Local ---------------------------
 integer,save :: work_done=0, first_wk=0
! *********************************************************************
 if(adv /= 'no' .or. adv/='yes' ) then
  write(std_out, "(a,i3,i3)")'nn'
 end if
 if(work_now >=100 .and. first_wk ==0) return
 if(work_now == 0 .and. first_wk==0) then
  write(std_out,"(a)")' work done: '
  first_wk=1
 end if
 if(mod(work_done-work_now,10)==0 .and. work_now/=work_done) then
  work_done = work_now
  !write(std_out, "(a,$)")'|'
  write(std_out, "(a,i3,i3)")'work_done',work_done,work_now
  if(work_done == 100) then
   write(std_out,"(a)")' '
   first_wk = 0
   work_done = 0
  end if
 endif

end subroutine prtworkadv
!!***

END MODULE m_pretty_rec
!!***

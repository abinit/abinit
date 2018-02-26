!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpar
!! NAME
!! inpar
!!
!! FUNCTION
!! Parser for the aim utility (shorter than the one of ABINIT)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  This routine uses data from the defs_aimprom module
!!
!! OUTPUT
!!  instr=string of character containing the input data
!!  lenstr=actual length of the character string
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inpar(instr,lenstr)

 use defs_basis
 use defs_aimprom
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inpar'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: lenstr
 character(len=*),intent(out) :: instr

!Local variables ------------------------------
 character(len=1),parameter :: space=' '
 character(len=26),parameter :: uplett='ABCDEFGHIJKLMNOPQRSTUVWXYZ', lolett='abcdefghijklmnopqrstuvwxyz'
!scalars
 integer,parameter :: nline=100
 integer :: ii,inxh,inxl,ios,jj,kk,ll
 character(len=fnlen) :: line

! *********************************************************************

 lenstr=0

 do ii=1,26
   inxh=index(lolett,uplett(ii:ii))
   if (inxh > 0) then
     write(std_out,*) 'ERROR The ', uplett(ii:ii) ,' is considered come lowcase !'
     MSG_ERROR("Aborting now")
   end if
 end do
 rewind(unt0)
 do ii=1,nline
   read(unt0,'(A)',iostat=ios) line(1:fnlen)
   if (ios/=0) exit
   inxh=index(line,'#')
   if (inxh == 1) then
     cycle
   elseif (inxh > 0) then
     inxl=inxh-1
     line(inxh:inxh)=space
   else
     inxl=len_trim(line)
     if (inxl==0) cycle
   end if
   inxh=index(line(1:inxl),char(9))
   if (inxh/=0) line(inxh:inxh)=space
   do ll=1,inxl
     if (iachar(line(ll:ll)) < 32) line(ll:ll)=space
   end do
   inxh=index(line(1:inxl),'- ')
   if (inxh/=0) then
     write(std_out,*) 'ERROR sign minus with white space in input file'
     MSG_ERROR("Aborting now")
   end if
   line(1:inxl)=adjustl(line(1:inxl))
   inxl=len_trim(line(1:inxl))+1
   jj=2;kk=0
   line(1:inxl)=adjustl(line(1:inxl))
   kk=len_trim(line(1:inxl))+1
   do ll=1,inxl
     inxh=index(line(jj:kk),space)
     if ((inxh==0).or.((jj+inxh-1)==kk)) exit
     line(inxh+jj:kk)=adjustl(line(inxh+jj:kk))
     kk=len_trim(line(1:inxl))
     if (kk == inxl) then
       exit
     end if
     jj=jj+inxh
   end do
   inxl=len_trim(line(1:inxl))+1
   do ll=1,inxl-1
     inxh=index(lolett,line(ll:ll))
     if (inxh/=0) line(ll:ll)=uplett(inxh:inxh)
   end do
   if ((lenstr+inxl) > strlen ) then
     write(std_out,*) 'ERROR Too large input !'
     MSG_ERROR("Aborting now")
   else
     instr(lenstr+1:lenstr+inxl)=line(1:inxl)
     lenstr=lenstr+inxl
   end if
 end do
end subroutine inpar
!!***

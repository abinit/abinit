!!****m* ABINIT/m_io_redirect
!! NAME
!!  m_io_redirect
!!
!! FUNCTION
!!  management of output and log files when parallelisation on cells is activated
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (FJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_io_redirect

 use defs_basis
 use m_abicore
 use m_errors

 use m_xmpi,        only : xmpi_comm_rank, xmpi_barrier
 use m_libpaw_tools,only : libpaw_write_comm_set
 use m_io_tools,    only : delete_file, flush_unit, open_file, get_unit
 use m_fstrings,    only : int2char4

 implicit none

 private

 character(len=fnlen),allocatable,save :: fillog(:),filout(:)

!Procedures ------------------------------------
 public :: localfilnam    ! Define new appended name
 public :: localwrfile    ! write local files
 public :: localrdfile    ! read local files
 public :: localredirect  ! redirect communicators for local files

contains
!!***

!!****f* m_io_redirect/localfilnam
!! NAME
!! localfilnam
!!
!! FUNCTION
!! Define new appended name
!!
!! INPUTS
!!   commspace= communicator between cells
!!   commworld= communicator of the world
!!   filnam(5)=character strings giving file names
!!   nam= name of the extension of the files
!!   nfil= number of files to be named
!!
!! PARENTS
!!      m_dfpt_looppert,m_gstateimg
!!
!! CHILDREN
!!      abi_io_redirect,libpaw_write_comm_set
!!
!! SOURCE
!!
 subroutine localfilnam(commspace,commspace1,commworld,filnam,nam,nfil)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: commspace,commspace1,commworld,nfil
 character(len=4) :: nam
!arrays
 character(len=fnlen),intent(in) :: filnam(:)

!Local variables-------------------------------
!scalars
 integer :: ii,ios,me,me_loc
 character(len=10) :: appen,tag

! *********************************************************************
 
 if (nfil<=1) return

 ABI_ALLOCATE(filout,(nfil))
 ABI_ALLOCATE(fillog,(nfil))
 me=xmpi_comm_rank(commworld)
 me_loc=xmpi_comm_rank(commspace1)
 call int2char4(me,tag)
 ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
 do ii=1,nfil
   call appdig(ii,'',appen)
   filout(ii)=trim(filnam(2))//nam//trim(appen)
   if (me_loc==0) then
     fillog(ii)=trim(filnam(5))//nam//trim(appen)//"_LOG"
   else
     fillog(ii)=trim(filnam(5))//nam//trim(appen)//"_LOG_"//trim(tag)
   end if
 end do

 if (me==0) then
   do ii=1,nfil
     call delete_file(filout(ii),ios)
   end do
 end if

 call xmpi_barrier(commspace)

 end subroutine localfilnam
!!***

!!****f* m_io_redirect/localwrfile
!! NAME
!!   localwrfile
!!
!! FUNCTION
!!   Write local (for each cell) output and log files
!!
!! INPUTS
!!   commspace= communicator inside cells
!!   ii=number of the iteration on the cells
!!   nfil= number of files to be named
!!   paral= flag to activate parallelisation
!!   prtvol= flag to activate printing
!!
!! PARENTS
!!      m_dfpt_looppert,m_gstateimg
!!
!! CHILDREN
!!      abi_io_redirect,libpaw_write_comm_set
!!
!! SOURCE
!!
 subroutine localwrfile(commspace,ii,nfil,paral,prtvol)

!Arguments ------------------------------------
 integer, intent(in) :: commspace,ii,nfil,paral,prtvol

!Local variables ------------------------------
 integer :: me
 character(len=500) :: msg

! *********************************************************************

 if (nfil<=1) return

 me=xmpi_comm_rank(commspace)
 if (paral==1) then
   call abi_io_redirect(new_io_comm=commspace)
   call libpaw_write_comm_set(commspace)
 end if
 if (prtvol>0) then
   if (do_write_log) then
     call abi_io_redirect(new_ab_out=get_unit())
     if (me==0) then 
       if (open_file(NULL_FILE,msg,unit=ab_out,status='unknown') /= 0) then
         MSG_ERROR(msg)
       end if
     end if
   else
     call abi_io_redirect(new_ab_out=std_out)
   end if
 else if (paral==1) then
   call abi_io_redirect(new_ab_out=get_unit())
   if (me==0) then 
     if (open_file(filout(ii),msg,unit=ab_out,status='unknown') /= 0) then
       MSG_ERROR(msg)
     end if
   end if
 end if
 if (paral==1.and.me==0.and.do_write_log) then
   call abi_io_redirect(new_std_out=get_unit())
   if (open_file(fillog(ii),msg,unit=std_out,status='unknown') /= 0) then
     MSG_ERROR(msg)
   end if
 end if

 end subroutine localwrfile
!!***

!!****f* m_io_redirect/localrdfile
!! NAME
!!   localrdfile
!!
!! FUNCTION
!!   read local (for each cell) output and log files
!!   to gather in one output or log file
!!
!! INPUTS
!!   commspace= communicator between cells
!!   commworld= communicator of the world
!!   compute_all= flag to activate the reading on all cells
!!   [dyn(nfil)]= --optional-- indicate if this cell is to be taken into account
!!   nfil= number of files to be named
!!   paral= flag to activate parallelisation
!!   prtvol= flag to activate printing
!!
!! PARENTS
!!      m_dfpt_looppert,m_gstateimg
!!
!! CHILDREN
!!      abi_io_redirect,libpaw_write_comm_set
!!
!! SOURCE
!!
 subroutine localrdfile(commspace,commworld,compute_all,nfil,paral,prtvol,dyn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: commspace,commworld,nfil,paral,prtvol
 logical,intent(in) :: compute_all
!arrays
 integer,intent(in),target,optional :: dyn(nfil)

!Local variables ------------------------------
!scalars
 integer :: ii,ios,me,temp_unit
 logical :: eof
 character(len=500) :: msg
 character(len=fnlen) :: line
!arrays
 integer,pointer :: my_dyn(:)

! *********************************************************************

 if (nfil<=1) return
 
 me=xmpi_comm_rank(commworld)
 call xmpi_barrier(commspace) !waiting for all files to be written and close

 if (prtvol<=0.and.(paral==1)) then
   if (me==0) then
     if (present(dyn)) then
       my_dyn => dyn
     else
       ABI_ALLOCATE(my_dyn,(nfil))
       my_dyn(:)=1
     end if
     do ii=1,nfil
       if (open_file(filout(ii),msg,newunit=temp_unit,status='old') == 0) then
         if (my_dyn(ii)==1.or.compute_all) then
           eof=.false.
           do while (.not.eof)
             read(temp_unit,fmt='(a)',err=111,end=111,iostat=ios) line
             if (ios==0) write(ab_out,'(a)') trim(line)
             goto 112
             111              eof=.true.
             112              continue
           end do
         end if
         close(unit=temp_unit,status='delete')
       end if
     end do
     call flush_unit(ab_out)
     if (.not.present(dyn)) then
       ABI_DEALLOCATE(my_dyn)
     end if
   end if
   call xmpi_barrier(commspace)
 end if

 if (paral==1.and.do_write_log) then
   if (me==0) then
     do ii=1,nfil
       if (open_file(fillog(ii),msg,newunit=temp_unit,status='old') == 0) then
         eof=.false.
         do while (.not.eof)
           read(temp_unit,fmt='(a)',err=113,end=113,iostat=ios) line
           if (ios==0) write(std_out,'(a)') trim(line)
           goto 114
           113            eof=.true.
           114            continue
         end do
         close(unit=temp_unit,status='delete')
       end if
     end do
   end if
 end if

 call xmpi_barrier(commspace)

 ABI_DEALLOCATE(filout)
 ABI_DEALLOCATE(fillog)

 end subroutine localrdfile
!!***

!!****f* m_io_redirect/localredirect
!! NAME
!!   localredirect
!!
!! FUNCTION
!!   Close output units ; restore defaults
!!
!! INPUTS
!!   commspace= communicator between cells
!!   commworld= communicator of the world
!!   nfil= number of files to be named
!!   paral= flag to activate parallelisation
!!   prtvol= flag to activate printing  
!!
!! PARENTS
!!      m_dfpt_looppert,m_gstateimg
!!
!! CHILDREN
!!      abi_io_redirect,libpaw_write_comm_set
!!
!! SOURCE
!!
 subroutine localredirect(commspace,commworld,nfil,paral,prtvol)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: commspace,commworld,nfil,paral,prtvol

!Local variables ------------------------------
 integer :: me

! *********************************************************************

 if (nfil<=1) return

 me=xmpi_comm_rank(commspace)
 if (prtvol>0) then
   if (do_write_log.and.me==0) close(unit=ab_out)
 else if (paral==1) then
   if (me==0) close(unit=ab_out)
 end if
 if (paral==1.and.me==0.and.do_write_log) close(unit=std_out)
 call abi_io_redirect(new_ab_out=ab_out_default,new_std_out=std_out_default,new_io_comm=commworld)
 call libpaw_write_comm_set(commworld)

 end subroutine localredirect
!!***

end module m_io_redirect
!!***

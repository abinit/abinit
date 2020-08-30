!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_io_tools
!! NAME
!!  m_io_tools
!!
!! FUNCTION
!!  This module contains basic tools to deal with Fortran IO
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_io_tools

 use defs_basis

 implicit none

 private

 public :: get_unit           ! Get a free unit if no argument is specified or report the unit associated to a file name
 public :: file_exists        ! Return .TRUE. if file exists.
 public :: delete_file        ! Delete a file if present.
 public :: is_open            ! .TRUE. if file is open
 public :: is_connected       ! .TRUE. if file is connected to a logical unit number
 public :: prompt             ! Simple prompt
 public :: read_string        ! Read string from unit ignoring blank lines and deleting comments beginning with ! or #
 public :: flush_unit         ! Wrapper to the intrinsic flush routine, not implemented by every compiler
 public :: pick_aname         ! Returns the name of a non-existent file to be used for temporary storage.
 public :: isncfile           ! .TRUE. if we have a NETCDF file.
 public :: iomode_from_fname  ! Automatic selection of the IO mode based on the file extension.
 public :: iomode2str         ! Convert iomode to string
 public :: enforce_fortran_io ! Set the value of enforce_fortran_io__
 public :: mvrecord           ! Moves forward or backward in a Fortran binary file by nn records.
 public :: open_file          ! Helper function to open a file in sequential mode with improved error handling.
 public :: close_unit         ! Helper function to close a Fortran unit with improved error handling.
 public :: write_lines        ! split a string in lines and output the text to the specified unit
 public :: lock_and_write     ! Write a string to a file with locking mechanism.
 public :: num_opened_units   ! Return the number of opened units.
 public :: show_units         ! Print info on the logical units.

 interface get_unit
   module procedure get_free_unit
   module procedure get_unit_from_fname
 end interface

 interface is_open
   module procedure is_open_unit
   module procedure is_open_fname
 end interface

 interface prompt
   module procedure prompt_int0D
   module procedure prompt_rdp0D
   module procedure prompt_string
   module procedure prompt_int1D
   module procedure prompt_int2D
   module procedure prompt_rdp1D
   module procedure prompt_rdp2D
 end interface

  integer,parameter :: MIN_UNIT_NUMBER=10  ! Fortran does not define the range for logical unit numbers (they not be negative)
#ifdef FC_NAG
  integer,parameter :: MAX_UNIT_NUMBER=64    ! There's a serious problem in Nag6.0. In principle
                                             ! Maximum unit number: 2147483647
#else
  integer,parameter :: MAX_UNIT_NUMBER=1024  ! The following values should be safe
#endif
  integer,parameter :: IO_MAX_LEN=500
  character(len=1),parameter :: BLANK=' '

  ! For interactive sessions
  integer,parameter :: IO_EOT=-1           ! End of transmission i.e CTRL+D
  !character(len=4),parameter :: PS1='>>> '
  ! Prepend prompt with `-` to bypass bug in intel18-19 so that flddiff.py will ignore the line
  character(len=4),parameter :: PS1='->> '
  character(len=4),parameter :: PS2='??? '

  integer,parameter :: IO_NO_AVAILABLE_UNIT  =-1   ! No units are available for Fortran I/O
  integer,parameter :: IO_FILE_NOT_ASSOCIATED=-2   ! File is not associated with any unit

  ! Enforce IO_MODE_FORTRAN in iomode_from_fname
  logical,save,protected :: enforce_fortran_io__ = .False.

CONTAINS  !===========================================================
!!***

!!****f* m_io_tools/get_unit
!! NAME
!!  get_unit
!!
!! FUNCTION
!!  Obtain a logical Fortran unit.
!!  A free unit is reported if no argument is specified.
!!  If the file name is supplied, the function reports the unit number
!!  associated to the file
!!  Note that GET_UNIT assumes that units 0, 5, 6 (stderr, stdin, std_out)
!!  are special, and will never return those values.
!!
!! TODO
!!   One should define an abinit-specific function with a list of reserved units!
!!
!! OUTPUT
!!  The unit number (free unit or unit associated to the file)
!!  Raises:
!!   IO_NO_AVAILABLE_UNIT if no logical unit is free (!)
!!   IO_FILE_NOT_ASSOCIATED if the file is not linked to a logical unit
!!
!! PARENTS
!!
!! SOURCE

integer function get_free_unit()

!Local variables-------------------------------
 integer :: iunt
 logical :: isopen
! *********************************************************************

 do iunt=MAX_UNIT_NUMBER,MIN_UNIT_NUMBER,-1
   if (any(iunt == [std_err, std_in, std_out])) cycle
   inquire(unit=iunt, opened=isopen)
   if (.not.isopen) then
      get_free_unit = iunt; return
   end if
 end do
 get_free_unit = IO_NO_AVAILABLE_UNIT

end function get_free_unit
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/get_unit_from_fname
!! NAME
!! get_unit_from_fname
!!
!! FUNCTION
!!  Returns the unit number associated to an open file whose name is fname.
!!  If the file is not connected to an unit number, returns IO_FILE_NOT_ASSOCIATED
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

integer function get_unit_from_fname(fname)

!Arguments ------------------------------------
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: unit
! *********************************************************************

 inquire(file=fname,number=unit)

 get_unit_from_fname=unit
 if (unit==-1) get_unit_from_fname=IO_FILE_NOT_ASSOCIATED

end function get_unit_from_fname
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/file_exists
!! NAME
!!  file_exists
!!
!! FUNCTION
!!  Return .TRUE. if file existent (function version of inquire).
!!
!! INPUTS
!!  fname=The name of the file.
!!
!! PARENTS
!!
!! SOURCE

logical function file_exists(fname)

 character(len=*),intent(in) :: fname

! *********************************************************************

 inquire(file=fname, exist=file_exists)

end function file_exists
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/delete_file
!! NAME
!!  delete_file
!!
!! FUNCTION
!!  Delete a file if present.
!!
!! INPUTS
!!  fname=The name of the file.
!!
!! OUTPUT
!!  ierr=Non-zero value indicates that a problem occured.
!!   111 = To signal that the file does not exist.
!!   112 = File exist, is open but no associated unit is found!
!!   Other values are system-dependent as the value is returned by a open or close
!!   instruction.
!!
!! SIDE EFFECTS
!!  The specified file is deleted.
!!
!! PARENTS
!!      ioprof,m_dvdb,m_io_redirect,m_mlwfovlp,m_nctk,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine delete_file(fname,ierr)

 integer,intent(out) :: ierr
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: tmp_unt
 logical :: exists
! *********************************************************************

 ierr = 0

 inquire(file=fname, exist=exists)

 if (.not.exists) then
   ierr = 111
   write(std_out,*)" Asked to delete not existent file: ",TRIM(fname)
   return
 end if

 if (is_open_fname(fname)) then
   tmp_unt = get_unit_from_fname(fname)
   if (tmp_unt == IO_FILE_NOT_ASSOCIATED) then
    write(std_out,*) "File is opened but no associated unit found!"
    ierr = 112; return
   end if
   close(tmp_unt)
 else
   tmp_unt = get_unit()
 end if

 ! Now close the file.
 open(unit=tmp_unt, file=trim(fname), status="OLD", iostat=ierr)
 if (ierr==0) close(unit=tmp_unt, status="DELETE", iostat=ierr)

end subroutine delete_file
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/is_connected
!! NAME
!!  is_connected
!!
!! FUNCTION
!!  Returns .TRUE. if unit is connected to fname.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

logical function is_connected(unit, fname)

 integer,intent(in) :: unit
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: unt_found
 logical :: isopen
! *********************************************************************

 inquire(file=fname, number=unt_found, opened=isopen)
 is_connected=(isopen .and. (unt_found == unit))

end function is_connected
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/is_open
!! NAME
!!  is_open
!!
!! FUNCTION
!!  Returns .TRUE. if unit is associated to an open file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

logical function is_open_unit(unit)

 integer,intent(in) :: unit
! *********************************************************************

 inquire(unit=unit, opened=is_open_unit)

end function is_open_unit
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/is_open_fname
!! NAME
!!  is_open_fname
!!
!! FUNCTION
!!  Returns .TRUE. if the file name fname is open.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

logical function is_open_fname(fname)

 character(len=*),intent(in) :: fname
! *********************************************************************

 inquire(file=fname,opened=is_open_fname)

end function is_open_fname
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_int0D
!! NAME
!!  prompt_int0D
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
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

subroutine prompt_int0D(msg,ivalue)

 character(len=*),intent(in) :: msg
 integer,intent(out) :: ivalue

!Local variables-------------------------------
 integer :: ios
 character(len=4) :: PS
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
  write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
  call flush_unit(std_out)
  read(std_in,*,IOSTAT=ios)ivalue
  if (ios==IO_EOT) call prompt_exit()
  PS=PS2
 end do
 write(std_out,*)

end subroutine prompt_int0D
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_rdp0d
!! NAME
!!  prompt_rdp0d
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
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

subroutine prompt_rdp0D(msg,rvalue)

 character(len=*),intent(in) :: msg
 real(dp),intent(out) :: rvalue

!Local variables-------------------------------
 integer :: ios
 character(len=4) :: PS
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
  write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
  call flush_unit(std_out)
  read(std_in,*,IOSTAT=ios)rvalue
  if (ios==IO_EOT) call prompt_exit()
  PS=PS2
 end do
 write(std_out,*)

end subroutine prompt_rdp0D
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_string
!! NAME
!!  prompt_string
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
!!  If strip_comment is True (default), all the characters after "#" or "!" are ignored.
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

subroutine prompt_string(msg,string,strip_comment)

 character(len=*),intent(in) :: msg
 logical,optional,intent(in) :: strip_comment
 character(len=*),intent(out) :: string

!Local variables-------------------------------
 integer :: ios,ic
 logical :: do_strip
 character(len=4) :: PS
 !character(len=len(string)) :: tmps
! *********************************************************************

 do_strip = .True.; if (present(strip_comment)) do_strip = strip_comment

 ios=-1 ; PS=PS1
 do while (ios/=0)
   write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
   call flush_unit(std_out)
   read(std_in,'(a)',IOSTAT=ios)string
   if (ios==IO_EOT) call prompt_exit()
   PS=PS2
 end do
 write(std_out,*)

 if (do_strip) then
   ic = INDEX(string, "#"); if (ic /= 0) string(:) = string(:ic-1)
   ic = INDEX(string, "!"); if (ic /= 0) string(:) = string(:ic-1)
 end if

end subroutine prompt_string
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_int1D
!! NAME
!!  prompt_int1D
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
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

subroutine prompt_int1D(msg,ivect)

!Arguments ------------------------------------
 character(len=*),intent(in) :: msg
 integer,intent(out) :: ivect(:)

!Local variables-------------------------------
 integer :: ios
 character(len=4) :: PS
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
   write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
   call flush_unit(std_out)
   read(std_in,*,IOSTAT=ios)ivect(:)
   if (ios==IO_EOT) call prompt_exit()
   PS=PS2
 end do
 write(std_out,*)

end subroutine prompt_int1D
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_int2D
!! NAME
!!  prompt_int2d
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
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

subroutine prompt_int2D(msg,iarr)

 character(len=*),intent(in) :: msg
 integer,intent(out) :: iarr(:,:)

!Local variables-------------------------------
 integer :: ios
 character(len=4) :: PS
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
   write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
   call flush_unit(std_out)
   read(std_in,*,IOSTAT=ios)iarr(:,:)
   if (ios==IO_EOT) call prompt_exit()
   PS=PS2
 end do
 write(std_out,*)

end subroutine prompt_int2D
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_rdp1D
!! NAME
!!  prompt_rdp1D
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
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

subroutine prompt_rdp1D(msg,rvect)

!Arguments ------------------------------------
 character(len=*),intent(in) :: msg
 real(dp),intent(out) :: rvect(:)
 character(len=4) :: PS
!Local variables-------------------------------
 integer :: ios
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
   write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
   call flush_unit(std_out)
   read(std_in,*,IOSTAT=ios)rvect(:)
   if (ios==IO_EOT) call prompt_exit()
   PS=PS2
 end do
 write(std_out,*)

end subroutine prompt_rdp1D
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_rdp2D
!! NAME
!!  prompt_rdp2D
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
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

subroutine prompt_rdp2D(msg,rarr)

 character(len=*),intent(in) :: msg
 real(dp),intent(out) :: rarr(:,:)

!Local variables-------------------------------
 integer :: ios
 character(len=4) :: PS
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
   write(std_out,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
   call flush_unit(std_out)
   read(std_in,*,IOSTAT=ios)rarr(:,:)
   if (ios==IO_EOT) call prompt_exit()
   PS=PS2
 end do
 write(std_out,*)

end subroutine prompt_rdp2D
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/prompt_exit
!! NAME
!!  prompt_exit
!!
!! FUNCTION
!!  A primitive prompt. Writes msg on std_out and reads the value entered by the user.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_io_tools
!!
!! CHILDREN
!!
!! SOURCE

subroutine prompt_exit()

 integer,parameter :: NASK=5
 integer :: ios,iask
 character(len=IO_MAX_LEN) :: ans
! *********************************************************************

 write(std_out,*)
 ios=-1 ; iask=0
 do while (ios/=0.or.(ans/='y'.or.ans/='n'))
   iask=iask+1
   write(std_out,'(a)')' Do you really want to exit (y/n)? '
   call flush_unit(std_out)
   read(std_in,*,IOSTAT=ios)ans
   if (ans=='y'.or.iask>NASK) STOP
   if (ans=='n') RETURN
 end do

end subroutine prompt_exit
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/read_string
!! NAME
!!  read_string
!!
!! FUNCTION
!!  Reads string from unit=std_in_ or unit if specified, ignoring blank lines
!!  and deleting comments beginning with `!`. Return exit code.
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

integer function read_string(string, unit) result(ios)

 character(len=*),intent(out):: string
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
 integer :: ipos,unt
! *********************************************************************

 unt=std_in; if (present(unit)) unt=unit

 read(unt,'(a)', iostat=ios) string  ! read input line
 if (ios/=0) return
 string = ADJUSTL(string)

 ! Ignore portion after comments
 ipos = INDEX(string, "!")
 if (ipos /= 0) string=string(:ipos-1)
 ipos = INDEX(string, "#")
 if (ipos /= 0) string=string(:ipos-1)

end function read_string
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/flush_unit
!! NAME
!! flush_unit
!!
!! FUNCTION
!! Wrapper for the standard flush_unit routine
!!
!! INPUTS
!!  unit=Fortran logical Unit number
!!
!! OUTPUT
!!
!! NOTES
!!  Available only if the compiler implements this intrinsic procedure.
!!
!! PARENTS
!!      abinit,anaddb,cut3d,fftprof,m_chi0,m_common,m_dmft,m_errors,m_forctqmc
!!      m_hdr,m_io_redirect,m_io_tools,m_matlu,m_mpi_setup,m_paw_mkaewf
!!      m_prep_calc_ucrpa,m_specialmsg,m_tdep_shell,m_vtorho,m_xc_vdw,mrggkk
!!      mrgscr,multibinit,optic,vdw_kernelgen
!!
!! CHILDREN
!!
!! SOURCE

subroutine flush_unit(unit)

 integer,intent(in) :: unit

!Local variables-------------------------------
 logical :: isopen

!************************************************************************

 if (unit == dev_null) return

 inquire(unit=unit,opened=isopen)

!FLUSH on unconnected unit is illegal: F95 std., 9.3.5.
#if defined HAVE_FC_FLUSH
 if (isopen) call flush(unit)
#elif defined HAVE_FC_FLUSH_
 if (isopen) call flush_(unit)
#endif

end subroutine flush_unit
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/pick_aname
!! NAME
!!  pick_aname
!!
!! FUNCTION
!!  Returns the name of a non-existent file to be used for temporary storage.
!!
!! PARENTS
!!
!! SOURCE

function pick_aname() result(aname)

 character(len=fnlen) :: aname

!Local variables-------------------------------
 integer :: ii,spt,ept
 real(dp) :: xrand(fnlen)
!************************************************************************

 aname="__TMP_FILE__"

 spt=LEN(aname); ept=spt

 do while (file_exists(aname))
   call RANDOM_NUMBER(xrand(spt:ept))
   xrand(spt:ept) = 64+xrand(spt:ept)*26
   do ii=spt,ept
     aname(ii:ii) = ACHAR(NINT(xrand(ii)))
   end do
   ept = MIN(ept+1,fnlen)
 end do

end function pick_aname
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/isncfile
!! NAME
!! isncfile
!!
!! FUNCTION
!!  Return .TRUE. if fname is a NETCDF file.
!!
!! INPUTS
!!  fname(len=*)=The name of the file to be tested.
!!
!! NOTES
!!  The idea is extremely simple: a NETCDF file terminates with ".nc".
!!  Obviously this approach is not bulletproof but it will work
!!  provided that we continue to append the ".nc" string to any NETCDF
!!  file produced by abinit.
!!
!! PARENTS
!!
!! SOURCE

pure logical function isncfile(fname)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
!scalars
 integer :: ic,nch_trim

! *************************************************************************

 nch_trim=LEN_TRIM(fname)
 ic = INDEX (TRIM(fname), ".", back=.TRUE.)

 isncfile=.FALSE.
 if (ic >= 1 .and. ic <= nch_trim-1) then ! there is stuff after the .
   isncfile = (fname(ic+1:nch_trim)=="nc")
 end if

end function isncfile
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/iomode_from_fname
!! NAME
!! iomode_from_fname
!!
!! FUNCTION
!!  Automatic selection of the IO mode based on the file extension.
!!
!! INPUTS
!!  fname = Name of the file.
!!
!! NOTES
!!  if fname has extension '.nc', IO_MODE_ETSF is used
!!  else:
!!    IO_MODE_MPI if available
!!    IO_MODE_FORTRAN if HAVE_MPI_IO is not defined.
!!
!! PARENTS
!!
!! SOURCE

pure integer function iomode_from_fname(fname) result(iomode)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname

! *************************************************************************

 if (isncfile(fname)) then
   iomode = IO_MODE_ETSF
 else
#ifdef HAVE_MPI_IO
   iomode = IO_MODE_MPI
#else
   iomode = IO_MODE_FORTRAN
#endif

   if (enforce_fortran_io__) iomode = IO_MODE_FORTRAN
 end if

end function iomode_from_fname
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/enforce_fortran_io
!! NAME
!! enforce_fortran_io
!!
!! FUNCTION
!!  Set the value of the enforce_fortran__ global variable.
!!
!! PARENTS
!!
!! SOURCE

subroutine enforce_fortran_io(bool)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: bool

! *************************************************************************

 enforce_fortran_io__ = bool

end subroutine enforce_fortran_io
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/iomode2str
!! NAME
!! iomode2str
!!
!! FUNCTION
!!  Convert iomode to string
!!
!! PARENTS
!!
!! SOURCE

pure function iomode2str(iomode)

!Arguments ------------------------------------
!scalars
 character(len=48) :: iomode2str
 integer,intent(in) :: iomode

! *************************************************************************

 select case (iomode)
 case (IO_MODE_FORTRAN_MASTER)
   iomode2str = "IO_MODE_FORTRAN_MASTER"
 case (IO_MODE_FORTRAN)
   iomode2str = "IO_MODE_FORTRAN"
 case (IO_MODE_MPI)
   iomode2str = "IO_MODE_MPI"
 case (IO_MODE_NETCDF)
   iomode2str = "IO_MODE_NETCDF"
 case (IO_MODE_ETSF)
   iomode2str = "IO_MODE_ETSF"
 case default
   iomode2str = "Unknown!"
 end select

end function iomode2str
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/mvrecord
!! NAME
!! mvrecord
!!
!! FUNCTION
!! This subroutine moves forward or backward in a Fortran binary file by nn records.
!!
!! INPUTS
!! funt= Fortran file unit number
!! nrec=number of records
!!
!! OUTPUT
!! ierr=error code
!!
!! TODO
!! One should treat the possible errors of backspace
!!
!! PARENTS
!!      m_wffile,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine mvrecord(funt,nrec,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: funt,nrec
 integer,intent(out) :: ierr

!Local variables-------------------------------
!scalars
 integer :: irec

! *************************************************************************

 ierr = 0
 if (nrec > 0) then ! Move forward nrec records
   do irec=1,nrec
     read(funt,iostat=ierr)
     if (ierr /= 0) EXIT
   end do
 else if (nrec < 0) then ! Move backward nrec records
   do irec=1,-nrec
     backspace (unit=funt,iostat=ierr)
     if (ierr /= 0) EXIT
   end do
 end if

end subroutine mvrecord
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/open_file
!! NAME
!! open_file
!!
!! FUNCTION
!!  Open a file in sequential mode and associate it to the unit number number.
!!  The main differences wrt the intrinsic open:
!!
!!    * Function statement that returns the value of iostat
!!    * Emulate iomsg (F2003)
!!    * Accepts either unit (user-specified unit number, input) or
!!      newunit (free unit not associated to any file, output).
!!      The two options are mutually exclusive.
!!
!!  See Fortran intrinsic for a more detailed description of the variables
!!
!! OUTPUT
!!  iostat=Exit status
!!  iomsg=Error message
!!
!! PARENTS
!!
!! SOURCE

function open_file(file,iomsg,unit,newunit,access,form,status,action,recl) result(iostat)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: file
 character(len=*),optional,intent(in) :: access,form,status,action
 character(len=*),intent(out) :: iomsg
 integer,optional,intent(in) :: recl,unit
 integer,optional,intent(out) :: newunit
 integer :: iostat

!Local variables-------------------------------
!scalars
 character(len=500) :: my_access,my_form,my_status,my_action,msg

! *************************************************************************

 my_access = "sequential"; if (present(access)) my_access = access
 my_form = "formatted"; if (present(form)) my_form = form
 my_status = "unknown"; if (present(status)) my_status = status
 my_action = "readwrite"; if (present(action)) my_action = action ! default is system dependent. Enforce RW mode

 iomsg = ""  ! iomsg is not changed if open succeeds

 if (present(unit)) then
   if (present(recl)) then
     open(file=trim(file),unit=unit,form=my_form,status=my_status,access=my_access,iostat=iostat,recl=recl, iomsg=iomsg)
   else
     open(file=trim(file),unit=unit,form=my_form,status=my_status,access=my_access,iostat=iostat, iomsg=iomsg)
   end if
   if (present(newunit)) iostat = -666 ! wrong call

 else if (present(newunit)) then
   ! Get free unit (emulate newunit of F2008)
   newunit = get_unit()
   if (present(recl)) then
     open(file=trim(file),unit=newunit,form=my_form,status=my_status,access=my_access,iostat=iostat,recl=recl, iomsg=iomsg)
   else
     open(file=trim(file),unit=newunit,form=my_form,status=my_status,access=my_access,iostat=iostat, iomsg=iomsg)
   end if
   if (present(unit)) iostat = -666  ! wrong call

 else
   iomsg = "Either unit or newunit must be specified"
   iostat = -1
 end if

 if (iostat /= 0) then
   write(msg, "(a,i0,2a)")"Fortran open returned iostat: ",iostat," while opening file: "//trim(file)
   iomsg = trim(msg)//ch10//"Runtime error message: "//trim(iomsg)
 end if

end function open_file
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/close_unit
!! NAME
!! close_unit
!!
!! FUNCTION
!!  close a Fortran unit
!!  The main differences wrt the intrinsic close:
!!
!!    * Function statement that returns the value of iostat
!!    * Emulate iomsg (F2003)
!!
!!  See Fortran intrinsic for a more detailed description of the variables
!!
!! OUTPUT
!!  iostat=Exit status
!!  iomsg=Error message
!!
!! PARENTS
!!
!! SOURCE

function close_unit(unit,iomsg,status) result(iostat)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: unit
 character(len=*),optional,intent(in) :: status
 character(len=*),intent(out) :: iomsg
 integer :: iostat

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 iomsg = "" ! iomsg is not changed if close succeeds

 if (.not.present(status)) then ! Use Fortran default e.g delete for scratch files.
#ifdef HAVE_FC_IOMSG
   close(unit=unit,iostat=iostat,iomsg=iomsg)
#else
   close(unit=unit,iostat=iostat)
#endif
 else
#ifdef HAVE_FC_IOMSG
   close(unit=unit,iostat=iostat,status=status,iomsg=iomsg)
#else
   close(unit=unit,iostat=iostat,status=status)
#endif
 end if

 ! TODO: Add more info for example the filename.
 if (iostat /= 0) then
   write(msg,'(2(a,i0),a)')"Fortran close returned iostat ",iostat," while closing unit: ",unit,ch10
   iomsg = trim(msg)//ch10//"IOMSG: "//trim(msg)
 end if

end function close_unit
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/write_lines
!! NAME
!!  write_lines
!!
!! FUNCTION
!!  This routine receives a string, split the message in lines according to the
!!  ch10 character and output the text to the specified unit
!!
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!  [toflush]=flag to activate immediate flush of the I/O buffer (default=FALSE)
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_io_tools,m_specialmsg
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_lines(unit,message,toflush)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 logical,intent(in),optional :: toflush
 character(len=*),intent(in) :: message

!Local variables-------------------------------
!scalars
 integer :: msg_size,ii,jj,rtnpos
 logical :: toflush_

!******************************************************************

 msg_size = len_trim(message)
 toflush_=.false.;if (present(toflush)) toflush_=toflush

 if (msg_size == 0) then
   write(unit,*)
   return
 end if

 ! Here, split the message, according to the char(10) characters (carriage return).
 ! This technique is portable accross different OS.
 rtnpos = index(message,ch10)

 if (rtnpos == 0) then
   write(unit,"(a)")message(1:msg_size)
   if (toflush_) call flush_unit(unit)
   return
 end if

 ii = 1; jj = rtnpos
 do
   if (ii == jj) then
     write(unit,*)
     if (toflush_) call flush_unit(unit)
   else
     write(unit, '(a)' ) message(ii:jj-1)
     if (toflush_) call flush_unit(unit)
   end if
   ii = jj + 1
   if (ii > msg_size) exit
   jj = index(message(ii:msg_size),ch10)
   if (jj == 0) then
     ! Will write the last line at the next iteration and exit .
     jj = msg_size + 1
   else
     jj = jj + ii - 1
   end if
   !write(*,*)"ii, jj, msg_size",ii, jj, msg_size
 end do

 ! This is needed to preserve the od behaviour: a ch10 at the
 ! end of the string was causing an extra newline!
 if (message(msg_size:msg_size) == ch10) write(unit,*)

end subroutine write_lines
!!***

!!****f* m_io_tools/lock_and_write
!! NAME
!!  lock_and_write
!!
!! FUNCTION
!!  Writes a string to filename with locking mechanism.
!!
!! INPUTS
!!  filename: Name of the file.
!!  string: Input string.
!!  ierr: Exit status, 0 is string has been written to filename.
!!
!! PARENTS
!!      m_errors
!!
!! CHILDREN
!!
!! SOURCE

subroutine lock_and_write(filename, string, ierr)

 integer,intent(out) :: ierr
 character(len=*),intent(in) :: filename,string

!Local variables-------------------------------
 integer :: lock_unit,file_unit
 character(len=len(filename) + 5) :: lock
 !character(len=500) :: msg

! *********************************************************************

 ierr = 0

 ! Try to acquire the lock.
 lock = trim(filename)//".lock"
 lock_unit = get_unit()
 open(unit=lock_unit, file=trim(lock), status='new', err=99)

 file_unit = get_unit()
 open(unit=file_unit, file=trim(filename), form="formatted")
 call write_lines(file_unit, string, toflush=.true.)
 close(lock_unit, status="delete")
 close(file_unit)
 return

99 ierr = 1

end subroutine lock_and_write
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/num_opened_units
!! NAME
!!  num_opened_units
!!
!! FUNCTION
!!  Return the number of opened units.
!!  Unit numbers listed in the optional argument `ignore` are not considered.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function num_opened_units(ignore) result(nn)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: ignore(:)

!Local variables-------------------------------
 integer :: ii,iostat
 logical  :: opened

! *********************************************************************

 nn = 0
 do ii=0, max_unit_number
   if (present(ignore)) then
     if (any(ii == ignore)) cycle
   end if
   inquire(ii, opened=opened, iostat=iostat)
   if (iostat == 0 .and. opened) nn = nn + 1
 end do

end function num_opened_units
!!***

!----------------------------------------------------------------------

!!****f* m_io_tools/show_units
!! NAME
!!  show_units
!!
!! FUNCTION
!!  Print info on the logical units
!!
!! PARENTS
!!      m_errors
!!
!! CHILDREN
!!
!! SOURCE

subroutine show_units(ount)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ount

!Local variables-------------------------------
 integer :: ii,iostat
 logical  :: named, opened
 character(len=fnlen) :: filename,form

! *********************************************************************

 write(ount,'(a)') '******** Fortran Logical Units ********'

 do ii=0,max_unit_number
   inquire(ii, opened=opened, named=named, name=filename, form=form, iostat=iostat)
   if (iostat == 0) then
      if (opened) then
         if (named) then
            write(ount,*)"unit: ", ii, "form: ", trim(form), ", filename: ", trim(filename)
         else
            write(ount,*)"unit: ", ii, "form: ",form, ', No name available'
         endif
      else
        !write(ount,*)"unit: ", ii, " is not opened"
      endif
   else
      write(ount,*)" unit: ", ii, ' Iostat error'
   endif
 end do

end subroutine show_units
!!***

!----------------------------------------------------------------------

END MODULE m_io_tools
!!***

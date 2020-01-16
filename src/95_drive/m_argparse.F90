!!****m* ABINIT/m_argparse
!! NAME
!! m_argparse
!!
!! FUNCTION
!!   Simple argument parser used in main programs
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_argparse

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use m_xieee
 use m_abi_linalg
 use m_fft
 use m_exit
 use m_clib

 use m_build_info,      only : dump_config, abinit_version
 use m_io_tools,        only : open_file, file_exists
 use m_cppopts_dumper,  only : dump_cpp_options
 use m_optim_dumper,    only : dump_optim
 use m_fstrings,        only : atoi, atof, itoa, firstchar, startswith, sjoin
 use m_time,            only : str2sec
 use m_libpaw_tools,    only : libpaw_log_flag_set

 implicit none

 private

 public :: get_arg         !  Parse scalar argument from command line. Return exit code.

 interface get_arg
   module procedure get_arg_int
   module procedure get_arg_dp
   module procedure get_arg_str
 end interface get_arg

 public :: get_arg_list    ! Parse array argument from command line. Return exit code.

 interface get_arg_list
   module procedure get_arg_list_int
   module procedure get_arg_list_dp
 end interface get_arg_list

 public :: parse_kargs    !  Parse command line arguments, return options related to k-point sampling
!!***

!!****t* m_argparse/args_t
!! NAME
!! args_t
!!
!! FUNCTION
!! Stores the command line options
!!
!! SOURCE

 type,public :: args_t

   integer :: exit = 0
     ! /=0 to exit after having parsed the command line options.

   integer :: abimem_level = 0

   integer :: dry_run = 0
     ! /= 0 to exit after the validation of the input file.

   real(dp) :: abimem_limit_mb = 20_dp
     ! Optional memory limit in Mb. used when abime_level == 3

   character(len=500) :: cmdline = ""
     ! The entire command line

   character(len=fnlen) :: input_path = ""

   !! Below are for multibinit
   integer :: multibinit_F03_mode = 0
   !1: legacy mode
   !0: use full F03 implementation mode
   ! TODO: It will be deprecated when everything is ready in
   ! and the new mode will be default.

 end type args_t

 public :: args_parser   ! Parse command line options.
!!***

contains

!----------------------------------------------------------------------

!!****f* m_argparse/args_parser
!! NAME
!!  args_parser
!!
!! FUNCTION
!!  Simple command line argument parser for abinit.
!!
!! PARENTS
!!      abinit
!!
!! SOURCE

type(args_t) function args_parser() result(args)

!Local variables-------------------------------
 integer :: ii,ierr
 logical :: iam_master,verbose
 real(dp) :: timelimit
 character(len=500) :: arg !,msg

! *************************************************************************

 args%exit = 0; ierr=0; verbose = .False.

#ifndef HAVE_FC_COMMAND_ARGUMENT
 call wrtout(std_out,"get_command_argument is not supported by FC. Ignoring command lines options!")
 return ! Nothing to do
#else

 if (command_argument_count() == 0) return

 iam_master = (xmpi_comm_rank(xmpi_world) == 0)

 ! Store full command line for future reference.
 call get_command(args%cmdline)

  do ii=1,command_argument_count()
    call get_command_argument(ii, arg)
    !write(std_out,*)"arg", TRIM(arg)

    if (ii == 1 .and. .not. firstchar(arg, "-")) then
      ! `abinit path` syntax reads input from path and deactivates files file mode.
      args%input_path = trim(arg)
      if (iam_master) then
         ABI_CHECK(file_exists(args%input_path), sjoin("Cannot find input file:", args%input_path))
      end if
      cycle
    end if

    if (arg == "-v" .or. arg == "--version") then
      call wrtout(std_out,TRIM(abinit_version),"COLL")
      args%exit = args%exit + 1

    else if (arg == "-b" .or. arg == "--build") then
      call print_kinds(unit=std_out)
      call xmpi_show_info(unit=std_out)
      call dump_cpp_options(std_out)
      call dump_config(std_out)
      call dump_optim(std_out)

      args%exit = args%exit + 1

    else if (arg == "-d" .or. arg == "--dry-run") then
      args%dry_run = 1

    else if (arg == "--abimem-level") then
      call get_command_argument(ii + 1, arg)
      args%abimem_level = atoi(arg)

    else if (arg == "--abimem-limit-mb") then
      call get_command_argument(ii + 1, arg)
      args%abimem_limit_mb = atof(arg)

    else if (arg == "-j" .or. arg == "--omp-num-threads") then
      call get_command_argument(ii + 1, arg)
      call xomp_set_num_threads(atoi(arg))

    ! timelimit handler.
    else if (arg == "-t" .or. arg == "--timelimit") then
      call get_command_argument(ii + 1, arg)
      timelimit = str2sec(arg)
      if (timelimit < zero) then
        write(std_out,*)"Wrong timelimit argument:",trim(arg)
        args%exit = args%exit + 1
      else
        call exit_init(timelimit)
      end if

    ! IEEE exceptions.
    else if (arg == "--ieee-halt") then
      call xieee_halt_ifexc(.True.)
    else if (arg == "--ieee-signal") then
      call xieee_signal_ifexc(.True.)

    ! Enable/disable non-blocking ialltoall in MPI-FFT
    else if (begins_with(arg, "--fft-ialltoall")) then
      call fft_allow_ialltoall(parse_yesno(arg, "--fft-ialltoall"))

    ! Enable/disable [Z,C]GEMM3
    else if (begins_with(arg, "--xgemm3m")) then
      call linalg_allow_gemm3m(parse_yesno(arg, "--xgemm3m"))

    ! Enable/disable PLASMA
    else if (begins_with(arg, "--plasma")) then
      call linalg_allow_plasma(parse_yesno(arg, "--plasma"))

    else if (arg == "--gnu-mtrace") then
      if (iam_master) then
        call clib_mtrace(ierr)
        ABI_CHECK(ierr == 0, sjoin("clib_mtrace returned ierr:", itoa(ierr)))
      end if

    else if (arg == "--log") then
       ! Enable logging
       call abi_log_status_state(new_do_write_log=.True., new_do_write_status=.True.)
       call libpaw_log_flag_set(.True.)

    !else if (arg == "-i" .or. arg == "--input") then
    !  call get_command_argument(ii + 1, arg)
    !  args%input_path = trim(arg)
    !  ! Redirect stdin to file.
    !  if (iam_master) then
    !     ABI_CHECK(file_exists(args%input_path), sjoin("Cannot find input file:", args%input_path))
    !     !close(std_in)
    !     !if (open_file(arg, msg, unit=std_in, form='formatted', status='old', action="read") /= 0) then
    !     !  MSG_ERROR(msg)
    !     !end if
    !  end if

    !else if (arg == "-o" .or. arg == "--output") then
    !  ! Redirect stdout to file.
    !  call get_command_argument(ii + 1, arg)
    !  if (iam_master) then
    !  close(std_out)
    !  if (open_file(arg, msg, unit=std_out, form='formatted', status='new', action="write") /= 0) then
    !    MSG_ERROR(message)
    !  end if

    ! For multibinit only
    else if (arg == "--F03") then
       args%multibinit_F03_mode = 1

    else if (arg == "-h" .or. arg == "--help") then
      if (iam_master) then
        ! Document the options.
        write(std_out,*)"-v, --version              Show version number and exit."
        write(std_out,*)"-b, --build                Show build parameters and exit."
        write(std_out,*)"-d, --dry-run              Validate input file and exit."
        write(std_out,*)"-j, --omp-num-threads      Set the number of OpenMp threads."
        write(std_out,*)"--abimem-level NUM         Set memory profiling level. Requires HAVE_MEM_PROFILING"
        write(std_out,*)"--abimem-limit-mb NUM      Log malloc/free only if size > limit in Megabytes. Requires abimem-level 3"
        write(std_out,*)"--ieee-halt                Halt the code if one of the *usual* IEEE exceptions is raised."
        write(std_out,*)"--ieee-signal              Signal the occurrence of the *usual* IEEE exceptions."
        write(std_out,*)"--fft-ialltoall[=bool]     Use non-blocking ialltoall in MPI-FFT (used only if ndat>1 and MPI3)."
        write(std_out,*)"--xgemm3m[=bool]           Use [Z,C]GEMM3M]"
        write(std_out,*)"--gnu-mtrace               Enable mtrace (requires GNU and clib)."
        write(std_out,*)"--log                      Enable log files and status files in parallel execution."
        write(std_out,*)"-i FILE                    Read input data from FILE instead of stdin."
        write(std_out,*)"                           Useful if MPI library does not handle `abinit < files` redirection"
        write(std_out,*)"-t, --timelimit            Set the timelimit for the run. Accepts time in Slurm notation:"
        write(std_out,*)"                               days-hours"
        write(std_out,*)"                               days-hours:minutes"
        write(std_out,*)"                               days-hours:minutes:seconds"
        write(std_out,*)"                               minutes"
        write(std_out,*)"                               minutes:seconds"
        write(std_out,*)"                               hours:minutes:seconds"
        write(std_out,*)"--verbose                  Verbose mode"
        write(std_out,*)"-h, --help                 Show this help and exit."

        ! Multibinit
        write(std_out,*)"--F03                      Run F03 mode (for Multibinit only)."
      end if
      args%exit = args%exit + 1

    else if (arg == "--verbose") then
      verbose = .True.

    else
      if (firstchar(arg, "-")) then
        MSG_WARNING("Unsupported option: "//trim(arg))
        args%exit = args%exit + 1
      else
        continue
      end if
    end if
  end do

  if (verbose) call args_print(args)
#endif

end function args_parser
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/args_print
!! NAME
!!  args_print
!!
!! FUNCTION
!!  Print object.
!!
!! PARENTS
!!      m_argparse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine args_print(args)

!Arguments ------------------------------------
 type(args_t),intent(in) :: args

! *************************************************************************

 call wrtout(std_out, sjoin("Command line:", args%cmdline))
 call wrtout(std_out, sjoin("exit:", itoa(args%abimem_level)))
 call wrtout(std_out, sjoin("abimem_level:", itoa(args%abimem_level)))
 call wrtout(std_out, sjoin("dry_run:", itoa(args%abimem_level)))

end subroutine args_print
!!***

!!****f* m_argparse/begins_with
!! NAME
!!  begins_with
!!
!! FUNCTION
!!  Returns true if argument arg begins with string
!!
!! SOURCE

pure logical function begins_with(arg, string) result(bool)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arg,string

 bool = .False.; if (len(arg) >= len(string)) bool = (arg(1:len(string)) == string)

end function begins_with
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/parse_yesno
!! NAME
!!  parse_yesno
!!
!! FUNCTION
!!  This function receives an argument, arg of the form --foo[=bool_value]
!!  that begins with optname (i.e. --foo) and returns the value of bool_value
!!  If bool_value is not present, returns default (.True. if not specified)
!!
!! SOURCE

logical function parse_yesno(arg, optname, default) result(bool)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arg,optname
 logical,optional,intent(in) :: default

 bool = .True.; if (present(default)) bool = default

 ! Assume default if value is not given
 if (len_trim(optname) == len_trim(arg)) return

 select case (arg(len(optname)+1:))
 case ("=yes", "=y")
   bool = .True.
 case ("=no", "=n")
   bool = .False.
 case default
   write(std_out,*)"Wrong option ",trim(arg),". Will default to ",bool
   MSG_ERROR("Aborting now")
 end select

end function parse_yesno
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/get_arg_int
!! NAME
!!  get_arg_int
!!
!! FUNCTION
!!  Parse scalar argument from command line. Return exit code.
!!
!! INPUT
!!  argname= Argument name
!!  [default]= Default value.
!!  [exclude]= argname and exclude are mutually exclusive.
!!
!! OUTPUT
!!  argval= Value of argname
!!  msg= Error message
!!
!! FUNCTION
!!
!! SOURCE

integer function get_arg_int(argname, argval, msg, default, exclude) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 integer,intent(out) :: argval
 character(len=*),intent(out) :: msg
 integer,optional,intent(in) :: default
 character(len=*),optional,intent(in) :: exclude

!Local variables-------------------------------
 integer :: ii, istat
 logical :: found_argname, found_excl
 character(len=500) :: arg, iomsg

! *************************************************************************

 ierr = 0; msg = ""; if (present(default)) argval = default
 found_argname = .False.; found_excl = .False.

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (present(exclude)) then
     if (arg == "--" // trim(exclude)) found_excl = .True.
   end if
   if (arg == "--" // trim(argname)) then
     found_argname = .True.
     call get_command_argument(ii + 1, arg, status=istat)
     if (istat == 0) then
       read(arg, *, iostat=istat, iomsg=iomsg) argval
       if (istat /= 0) then
         ierr = ierr + 1; msg = sjoin(msg, ch10, iomsg)
       end if
     else
       ierr = ierr + 1; msg = sjoin(msg, ch10, "Error in get_command_argument")
     end if
   end if
 end do

 if (ierr /= 0) msg = sjoin("Error while reading argument: ", argname, ch10, msg)
 if (found_argname .and. found_excl) then
   ierr = ierr + 1; msg = sjoin("Variables", argname, "and", exclude, "are mutually exclusive", ch10, msg)
 end if

end function get_arg_int
!!***

!!****f* m_argparse/get_arg_dp
!! NAME
!!  get_arg_dp
!!
!! FUNCTION
!!  Parse scalar argument from command line. Return exit code.
!!
!! INPUTS
!!  argname= Argument name
!!  [default]= Default value
!!  [exclude]= argname and exclude are mutually exclusive.
!!
!! OUTPUT
!!   argval= Value of argname
!!   msg= Error message
!!
!! SOURCE

integer function get_arg_dp(argname, argval, msg, default, exclude) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 real(dp),intent(out) :: argval
 character(len=*),intent(out) :: msg
 real(dp),optional,intent(in) :: default
 character(len=*),optional,intent(in) :: exclude

!Local variables-------------------------------
 integer :: ii, istat
 logical :: found_argname, found_excl
 character(len=500) :: arg, iomsg

! *************************************************************************

 ierr = 0; msg = ""; if (present(default)) argval = default
 found_argname = .False.; found_excl = .False.

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (present(exclude)) then
     if (arg == "--" // trim(exclude)) found_excl = .True.
   end if
   if (arg == "--" // trim(argname)) then
     found_argname = .True.
     call get_command_argument(ii + 1, arg, status=istat)
     if (istat == 0) then
       read(arg, *, iostat=istat, iomsg=iomsg) argval
       if (istat /= 0) then
         ierr = ierr + 1; msg = sjoin(msg, ch10, iomsg)
       end if
     else
       ierr = ierr + 1; msg = sjoin(msg, ch10, "Error in get_command_argument")
     end if
   end if
 end do

 if (ierr /= 0) msg = sjoin("Error while reading argument: ", argname, ch10, msg)
 if (found_argname .and. found_excl) then
   ierr = ierr + 1; msg = sjoin("Variables", argname, "and", exclude, "are mutually exclusive", ch10, msg)
 end if

end function get_arg_dp
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/get_arg_str
!! NAME
!!  get_arg_str
!!
!! FUNCTION
!!  Parse scalar argument from command line. Return exit code.
!!
!! INPUTS
!!  argname= Argument name
!!  [default]= Default value
!!  [exclude]= argname and exclude are mutually exclusive.
!!
!! OUTPUT
!!   argval= Value of argname
!!   msg= Error message
!!
!! SOURCE

integer function get_arg_str(argname, argval, msg, default, exclude) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 character(len=*),intent(out) :: argval
 character(len=*),intent(out) :: msg
 character(len=*),optional,intent(in) :: default
 character(len=*),optional,intent(in) :: exclude

!Local variables-------------------------------
 integer :: ii, istat
 logical :: found_argname, found_excl
 character(len=500) :: arg

! *************************************************************************

 ierr = 0; msg = ""; if (present(default)) argval = default
 found_argname = .False.; found_excl = .False.

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (present(exclude)) then
     if (arg == "--" // trim(exclude)) found_excl = .True.
   end if
   if (arg == "--" // trim(argname)) then
     found_argname = .True.
     call get_command_argument(ii + 1, argval, status=istat)
     if (istat /= 0) then
       ierr = ierr + 1; msg = sjoin(msg, ch10, "Error in get_command_argument")
     end if
   end if
 end do

 if (ierr /= 0) msg = sjoin("Error while reading argument: ", argname, ch10, msg)
 if (found_argname .and. found_excl) then
   ierr = ierr + 1; msg = sjoin("Variables", argname, "and", exclude, "are mutually exclusive", ch10, msg)
 end if

end function get_arg_str
!!***

!!****f* m_argparse/get_arg_list_int
!! NAME
!!  get_arg_list_int
!!
!! FUNCTION
!!  Parse array argument from command line. Return exit code.
!!
!! INPUT
!!  argname= Argument name
!!  [default]= Default value (scalar)
!!  [default_list]= Default value (vector)
!!  [exclude]= argname and exclude are mutually exclusive.
!!  [want_len]= Require want_len items in CLI.
!!
!! OUTPUT
!!  argval= Value of argname
!!  msg= Error message
!!
!! SOURCE

integer function get_arg_list_int(argname, argval, lenr, msg, default, default_list, exclude, want_len) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 integer,intent(out) :: argval(:)
 integer,intent(out) :: lenr
 character(len=*),intent(out) :: msg
 character(len=*),optional,intent(in) :: exclude
 integer,optional,intent(in) :: default
 integer,optional,intent(in) :: default_list(:)
 integer,optional,intent(in) :: want_len

!Local variables-------------------------------
 integer :: ii, istat, iarg, maxlen
 logical :: found_argname, found_excl
 character(len=500) :: arg, iomsg

! *************************************************************************

 ierr = 0; msg = ""; lenr = 0
 found_argname = .False.; found_excl = .False.

 maxlen = size(argval);
 if (maxlen == 0) then
   ierr = ierr + 1; msg = "zero-sized argval!"; return
 end if

 if (present(default)) argval = default
 if (present(default_list)) argval = default_list

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (present(exclude)) then
     if (arg == "--" // trim(exclude)) found_excl = .True.
   end if
   if (arg == "--" // trim(argname)) then
     ! Read list of values
     found_argname = .True.
     do iarg=1,maxlen
       call get_command_argument(ii + iarg, arg, status=istat)
       if (istat == 0) then
         !write(std_out, *)"arg:", trim(arg)
         if (startswith(arg, "--")) exit
         read(arg,*, iostat=istat, iomsg=iomsg) argval(iarg)
         if (istat == 0) then
           lenr = lenr + 1
         else
           ierr = ierr + 1; msg = sjoin(msg, ch10, iomsg)
         end if
       else
         ! If there are less than NUMBER arguments specified at the command line, VALUE will be filled with blanks.
         if (arg == "") exit
         ierr = ierr + 1; msg = sjoin(msg, ch10, "Error in get_command_argument")
       end if
     end do
   end if
 end do

 if (ierr /= 0) msg = sjoin("Error while reading argument: ", argname, ch10, msg)
 if (found_argname .and. found_excl) then
   ierr = ierr + 1; msg = sjoin("Variables", argname, "and", exclude, "are mutually exclusive", ch10, msg)
 end if

 if (present(want_len)) then
   if (found_argname) then
     if (want_len /= lenr) then
       ierr = ierr + 1
       msg = sjoin(argname, "requires", itoa(want_len), " tokens while found ", itoa(lenr), ch10, msg)
     end if
   else
     ierr = ierr + 1
     msg = sjoin("Cannot find --", argname, "option in CLI and want_len:", itoa(want_len), ch10, msg)
   end if
 end if

end function get_arg_list_int
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/get_arg_list_dp
!! NAME
!!  get_arg_list_dp
!!
!! FUNCTION

!! INPUT
!!  argname
!!  [default]
!!  [default_list]
!!  [exclude]
!!  [want_len]
!!
!! OUTPUT
!!  argval
!!  msg
!!
!! SOURCE

integer function get_arg_list_dp(argname, argval, lenr, msg, default, default_list, exclude, want_len) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 real(dp),intent(out) :: argval(:)
 integer,intent(out) :: lenr
 character(len=*),intent(out) :: msg
 character(len=*),optional,intent(in) :: exclude
 real(dp),optional,intent(in) :: default
 real(dp),optional,intent(in) :: default_list(:)
 integer,optional,intent(in) :: want_len

!Local variables-------------------------------
 integer :: ii, istat, iarg, maxlen
 logical :: found_argname, found_excl
 character(len=500) :: arg, iomsg

! *************************************************************************

 ierr = 0; msg = ""; lenr = 0
 found_argname = .False.; found_excl = .False.

 maxlen = size(argval);
 if (maxlen == 0) then
   ierr = ierr + 1; msg = "zero-sized argval!"; return
 end if

 if (present(default)) argval = default
 if (present(default_list)) argval = default_list

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (present(exclude)) then
     if (arg == "--" // trim(exclude)) found_excl = .True.
   end if
   if (arg == "--" // trim(argname)) then
     ! Read list of values
     found_argname = .True.
     do iarg=1,maxlen
       call get_command_argument(ii + iarg, arg, status=istat)
       if (istat == 0) then
         !write(std_out, *)"arg:", trim(arg)
         if (startswith(arg, "--")) exit
         read(arg,*, iostat=istat, iomsg=iomsg) argval(iarg)
         if (istat == 0) then
           lenr = lenr + 1
         else
           ierr = ierr + 1; msg = sjoin(msg, ch10, iomsg)
         end if
       else
         ! If there are less than NUMBER arguments specified at the command line, VALUE will be filled with blanks.
         if (arg == "") exit
         ierr = ierr + 1; msg = sjoin(msg, ch10, "Error in get_command_argument")
       end if
     end do
   end if
 end do

 if (ierr /= 0) msg = sjoin("Error while reading argument: ", argname, ch10, msg)
 if (found_argname .and. found_excl) then
   ierr = ierr + 1; msg = sjoin("Variables", argname, "and", exclude, "are mutually exclusive", ch10, msg)
 end if

 if (present(want_len)) then
   if (found_argname) then
     if (want_len /= lenr) then
       ierr = ierr + 1
       msg = sjoin(argname, "requires", itoa(want_len), " tokens while found ", itoa(lenr), ch10, msg)
     end if
   else
     ierr = ierr + 1
     msg = sjoin("Cannot find --", argname, " option in CLI and want_len:", itoa(want_len), ch10, msg)
   end if
 end if

end function get_arg_list_dp
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/parse_kargs
!! NAME
!! parse_kargs
!!
!! FUNCTION
!!  Parse command line arguments, return options related to k-point sampling
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

subroutine parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)

!Arguments ------------------------------------
 integer,intent(out) :: kptopt, nshiftk, chksymbreak
 integer,intent(out) :: kptrlatt(3,3)
 real(dp),allocatable,intent(out) :: shiftk(:,:)

!Local variables-------------------------------
 integer :: ii, lenr
 character(len=500) :: msg
 integer :: ivec9(9), ngkpt(3)
 real(dp) :: my_shiftk(3 * MAX_NSHIFTK)

! *************************************************************************

 ABI_CHECK(get_arg("kptopt", kptopt, msg, default=1) == 0, msg)
 ABI_CHECK(get_arg("chksymbreak", chksymbreak, msg, default=1) == 0, msg)
 ABI_CHECK(get_arg_list("ngkpt", ngkpt, lenr, msg, exclude="kptrlatt", want_len=3) == 0, msg)
 if (lenr == 3) then
   kptrlatt = 0
   do ii=1,3
     kptrlatt(ii, ii) = ngkpt(ii)
   end do
 end if
 ABI_CHECK(get_arg_list("kptrlatt", ivec9, lenr, msg, exclude="ngkpt", want_len=9) == 0, msg)
 if (lenr == 9) kptrlatt = transpose(reshape(ivec9, [3, 3]))

 ! Init default
 ABI_CHECK(get_arg_list("shiftk", my_shiftk, lenr, msg) == 0, msg)
 if (lenr /= 0) then
   ABI_CHECK(mod(lenr, 3) == 0, "Expecting 3 * nshift array")
   nshiftk = lenr / 3
   ABI_MALLOC(shiftk, (3, nshiftk))
   shiftk = reshape(my_shiftk(1:lenr), [3, nshiftk])
 else
   nshiftk = 1
   ABI_MALLOC(shiftk, (3, nshiftk))
   shiftk(:, 1) = [half, half, half]
 end if
 !write(std_out, *)"kptopt = ", kptopt, ", chksymbreak = ", chksymbreak, ", nshiftk = ", nshiftk, ", kptrlatt = ", kptrlatt

end subroutine parse_kargs
!!***

end module m_argparse
!!***

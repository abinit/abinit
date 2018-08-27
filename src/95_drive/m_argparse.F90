!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_argparse
!! NAME
!! m_argparse
!!
!! FUNCTION
!!   Simple argument parser used in main programs
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
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

 use m_build_info,      only : dump_config, abinit_version
 use m_io_tools,        only : open_file
 use m_cppopts_dumper,  only : dump_cpp_options
 use m_optim_dumper,    only : dump_optim
 use m_fstrings,        only : atoi, itoa, firstchar, sjoin
 use m_time,            only : str2sec
 use m_libpaw_tools,    only : libpaw_log_flag_set

 implicit none

 private
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

   integer :: exit=0
     ! /=0 to exit after having parsed the command line options.

   integer :: abimem_level=0

   integer :: dry_run=0
     ! /= 0 to exit after the validation of the input file.

   character(len=500) :: cmdline=""
     ! The entire command line

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'args_parser'
!End of the abilint section

 implicit none

!Local variables-------------------------------
 integer :: ii,ierr
 logical :: iam_master,verbose
 real(dp) :: timelimit
 character(len=500) :: arg,msg

! *************************************************************************

 args%exit = 0; ierr=0; verbose = .False.

#ifndef HAVE_FC_COMMAND_ARGUMENT
 call wrtout(std_out,"get_command_argument is not supported by FC. Ignoring command lines options!")
 return ! Nothing to do
#else

 if (command_argument_count()==0) return

 iam_master = (xmpi_comm_rank(xmpi_world) == 0)

 ! Store full command line for future reference.
 call get_command(args%cmdline)

  do ii=1,command_argument_count()
    call get_command_argument(ii, arg)
    !write(std_out,*)"arg", TRIM(arg)

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
      call get_command_argument(ii+1, arg)
      args%abimem_level = atoi(arg)

    else if (arg == "-j" .or. arg == "--omp-num-threads") then
      call get_command_argument(ii+1, arg)
      call xomp_set_num_threads(atoi(arg))

    ! timelimit handler.
    else if (arg == "-t" .or. arg == "--timelimit") then
      call get_command_argument(ii+1, arg)
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
        if (ierr/=0 .and. iam_master) write(std_out,"(a,i0)")"clib_mtrace returned ierr: ",ierr
    end if

    else if (arg == "--log") then
       ! Enable logging
       call abi_log_status_state(new_do_write_log=.True.,new_do_write_status=.True.)
       call libpaw_log_flag_set(.True.)

    else if (arg == "-i" .or. arg == "--input") then
      ! Redirect stdin to file.
      call get_command_argument(ii+1, arg)
      if (iam_master) then
        close(std_in)
        if (open_file(arg, msg, unit=std_in, form='formatted', status='old', action="read") /= 0) then
          MSG_ERROR(msg)
        end if
      end if

    !else if (arg == "-o" .or. arg == "--output") then
    !  ! Redirect stdout to file.
    !  call get_command_argument(ii+1, arg)
    !  if (iam_master) then
    !  close(std_out)
    !  if (open_file(arg, msg, unit=std_out, form='formatted', status='new', action="write") /= 0) then
    !    MSG_ERROR(message)
    !  end if
    !  end if

    else if (arg == "-h" .or. arg == "--help") then
      if (iam_master) then
        ! Document the options.
        write(std_out,*)"-v, --version              Show version number and exit."
        write(std_out,*)"-b, --build                Show build parameters and exit."
        write(std_out,*)"-d, --dry-run              Validate input file and exit."
        write(std_out,*)"-j, --omp-num-threads      Set the number of OpenMp threads."
        write(std_out,*)"--abimem-level NUM         Set memory profiling level. Requires HAVE_MEM_PROFILING"
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'args_print'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'begins_with'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arg,string

 bool = .False.; if (len(arg) >= len(string)) bool = (arg(1:len(string)) == string)

end function begins_with
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/parse_yesno
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'parse_yesno'
!End of the abilint section

 implicit none

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

end module m_argparse
!!***

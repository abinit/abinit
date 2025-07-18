!!****m* ABINIT/m_argparse
!! NAME
!! m_argparse
!!
!! FUNCTION
!!   Simple argument parser used in main programs
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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
 use m_nctk

 use m_build_info,      only : dump_config, abinit_version
 use m_io_tools,        only : open_file, file_exists, enforce_fortran_io
 use m_cppopts_dumper,  only : dump_cpp_options
 use m_optim_dumper,    only : dump_optim
 use m_fstrings,        only : atoi, atof, itoa, firstchar, startswith, endswith, sjoin, find_and_select
 use m_time,            only : str2sec
 use m_libpaw_tools,    only : libpaw_log_flag_set
 use m_ipi,             only : ipi_setup

 implicit none

 private

 public :: get_arg         !  Parse scalar argument from command line. Return exit code.

 interface get_arg
   module procedure get_arg_int
   module procedure get_arg_dp
   module procedure get_arg_str
   module procedure get_arg_bool
 end interface get_arg

 public :: get_arg_list    ! Parse array argument from command line. Return exit code.

 interface get_arg_list
   module procedure get_arg_list_int
   module procedure get_arg_list_dp
 end interface get_arg_list

 public :: get_start_step_num    ! Parse string from command line in the format "start:step:num"
                                 ! defining an arithmetic progression.
 public :: parse_kargs           !  Parse command line arguments, return options related to k-point sampling
!!***

!!****t* m_argparse/args_t
!! NAME
!! args_t
!!
!! FUNCTION
!! Stores command line options
!!
!! SOURCE

 type,public :: args_t

   integer :: exit = 0
    ! /=0 to exit after having parsed the command line options.

   integer :: abimem_level = 0
    ! Options for memory profiling. See m_profiling_abi

   integer :: dry_run = 0
    ! /= 0 to exit after the validation of the input file.

   real(dp) :: abimem_limit_mb = 20.0_dp
    ! Optional memory limit in Mb. used when abimem_level == 3

   character(len=500) :: cmdline = ""
    ! The entire command line

   character(len=fnlen) :: input_path = ""

   !! Below are for multibinit
   integer :: multibinit_F03_mode = 0
   !1: legacy mode
   !0: use full F03 implementation mode
   ! TODO: It will be deprecated when everything is ready in and the new mode will be default.

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
!!  Simple command line argument parser for abinit and other main programs.
!!
!! SOURCE

type(args_t) function args_parser() result(args)

!Local variables-------------------------------
 integer :: ii, ierr, ntasks_per_node = -1
 logical :: iam_master, verbose
 real(dp) :: timelimit, memb_per_node = -one, memb_per_cpu = -one
 character(len=500) :: arg !,msg

! *************************************************************************

 args%exit = 0; ierr=0; verbose = .False.

#ifndef HAVE_FC_COMMAND_ARGUMENT
 call wrtout(std_out,"get_command_argument is not supported by FC. Ignoring command lines options!")
 return ! Nothing to do
#else

 if (command_argument_count() == 0) return

 iam_master = xmpi_comm_rank(xmpi_world) == 0

 ! Store full command line for future reference.
 call get_command(args%cmdline)

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   !write(std_out,*)"arg", trim(arg)

   if (ii == 1 .and. .not. firstchar(arg, "-")) then
     ! `abinit path` syntax reads input from path and deactivates files file mode.
     args%input_path = trim(arg)
     if (iam_master) then
        ABI_CHECK(file_exists(args%input_path), sjoin("Cannot find input file:", args%input_path))
     end if
     cycle
   end if

   if (arg == "-v" .or. arg == "--version") then
     call wrtout(std_out, trim(abinit_version))
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

   else if (arg == "-t" .or. arg == "--timelimit") then
     ! timelimit handler.
     call get_command_argument(ii + 1, arg)
     timelimit = str2sec(arg)
     if (timelimit < zero) then
       write(std_out,*)"Wrong timelimit argument: ",trim(arg)
       args%exit = args%exit + 1
     else
       call exit_init(timelimit)
     end if

   else if (arg == "--ieee-halt") then
     ! IEEE exceptions.
     call xieee_halt_ifexc(.True.)

   else if (arg == "--ieee-signal") then
     call xieee_signal_ifexc(.True.)

   else if (begins_with(arg, "--fft-ialltoall")) then
     ! Enable/disable non-blocking ialltoall in MPI-FFT
     call fft_allow_ialltoall(parse_yesno(arg, "--fft-ialltoall"))

   else if (begins_with(arg, "--ipi")) then
     call get_command_argument(ii + 1, arg)
     call ipi_setup(arg, xmpi_world)

   else if (begins_with(arg, "--use-xgemm3m")) then
     ! Enable/disable [Z,C]GEMM3
     call linalg_allow_gemm3m(parse_yesno(arg, "--use-xgemm3m"), write_msg=iam_master)

   else if (begins_with(arg, "--use-mpi-in-place")) then
     ! Enable/disable usage of MPI_IN_PLACE.
     call xmpi_set_inplace_operations(parse_yesno(arg, "--use-mpi-in-place"))

   else if (begins_with(arg, "--plasma")) then
     ! Enable/disable PLASMA
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

   else if (arg == "--netcdf-classic") then
     ! Use netcdf classic mode for new files when only sequential-IO needs to be performed
     call nctk_use_classic_for_seq()

   else if (arg == "--enforce-fortran-io") then
     call enforce_fortran_io(.True.)

   else if (begins_with(arg, "--mem-per-cpu=")) then
     memb_per_cpu = parse_slurm_mem(arg, "--mem-per-cpu=")
     call set_mem_per_cpu_mb(memb_per_cpu)

   else if (begins_with(arg, "--mem=")) then
     memb_per_node = parse_slurm_mem(arg, "--mem=")

   else if (begins_with(arg, "--ntasks-per-node=")) then
     call get_command_argument(ii + 1, arg)
     ntasks_per_node = atoi(arg)

   else if (arg == "--F03") then
     ! For multibinit only
     args%multibinit_F03_mode = 1

   else if (arg == "-h" .or. arg == "--help") then
     if (iam_master) then
       ! Document the options.
       write(std_out,*)"-v, --version              Show version number and exit."
       write(std_out,*)"-b, --build                Show build parameters and exit."
       write(std_out,*)"-d, --dry-run              Validate input file and exit."
       write(std_out,*)"-j, --omp-num-threads      Set the number of OpenMp threads."
       write(std_out,*)"--use-xgemm3m[=yesno]      Use ZGEMM3M routines instead of ZGEMM. Default: no "
       write(std_out,*)"--use-mpi-in-place[=yesno] Enable/disable usage of MPI_IN_PLACE in e.g. xmpi_sum. Default: no"
       write(std_out,*)"                           Note that some MPI libs e.g. intel-mpi may not implement this feature"
       write(std_out,*)"                           correctly so it is adviced to test this option with e.g. structural"
       write(std_out,*)"                           relaxations before running production calculations."
       write(std_out,*)"--ipi                      Activate socket-driven calculation using i-pi protocol."
       write(std_out,*)"                           For UNIX socket, use: --ipi {unixsocket}:UNIX"
       write(std_out,*)"                           For INET socket, use  --ipi {host}:{port}. Usage example:"
       write(std_out,*)"                           `abinit run.abi --ipi {unixsocket}:UNIX > run.log`"
       write(std_out,*)"                           NB: Requires ionmov 28 and some tuning of input variables. See:"
       write(std_out,*)"                           https://wiki.fysik.dtu.dk/ase/dev/ase/calculators/socketio/socketio.html"
       write(std_out,*)"--log                      Enable log files and status files in parallel execution."
       write(std_out,*)"--netcdf-classic           Use netcdf classic mode for new files if parallel-IO is not needed."
       write(std_out,*)"                           Default is netcdf4/hdf5"
       write(std_out,*)"--enforce-fortran-io       Use Fortran-IO instead of MPI-IO when operating on Fortran files"
       write(std_out,*)"                           Useful to read files when the MPI-IO library is not efficient."
       write(std_out,*)"                           DON'T USE this option when the code needs to write large files e.g. WFK"
       write(std_out,*)"-t, --timelimit            Set the timelimit for the run. Accepts time in Slurm syntax:"
       write(std_out,*)"                               days-hours"
       write(std_out,*)"                               days-hours:minutes"
       write(std_out,*)"                               days-hours:minutes:seconds"
       write(std_out,*)"                               minutes"
       write(std_out,*)"                               minutes:seconds"
       write(std_out,*)"                               hours:minutes:seconds"
       write(std_out,*)"                           At present only GS, relaxations and MD runs support this option"
       write(std_out,*)"--mem-per-cpu=<size>[units] Set memory per cpu using Slurm syntax. Default units are megabytes."
       write(std_out,*)"                           Different units can be specified using the suffix [K|M|G|T]."
       write(std_out,*)"--mem=<size>[units]        Set memory per node using Slurm syntax. Default units are megabytes."
       write(std_out,*)"                           Requires `ntasks-per-node`. Not compatibile with `-mem-per-cpu`."
       write(std_out,*)"--ntasks-per-node=INT      Set number of tasks per node. Used in conjunction with --mem`"
       write(std_out,*)"--verbose                  Enable verbose mode in argparse"
       write(std_out,*)"-h, --help                 Show this help and exit."

       write(std_out,*)""
       write(std_out,*)""
       write(std_out,*)"=============================="
       write(std_out,*)"=== Options for developers ==="
       write(std_out,*)"=============================="
       write(std_out,*)"--abimem-level NUM         Set memory profiling level. Requires HAVE_MEM_PROFILING in config.h"
       write(std_out,*)"       0 -> no file abimem.mocc is created, only memory allocation counters running."
       write(std_out,*)"       1 -> light version. Only memory peaks are written."
       write(std_out,*)"       2 -> file abimem.mocc is created with full information inside."
       write(std_out,*)"       3 -> Write info only if allocation/deallocation is larger or smaller than limit_mb"
       write(std_out,*)"            depending on of the sign of abimem-limit-mb."
       write(std_out,*)"    NOTE: By default, only master node writes, use negative values to make all MPI procs write info to disk."
       write(std_out,*)"--abimem-limit-mb NUM      Log malloc/free only if size > limit in Megabytes. Requires abimem-level 3"
       write(std_out,*)"--fft-ialltoall[=yesno]    Use non-blocking ialltoall in MPI-FFT (used only if ndat > 1 and MPI2+)."
       write(std_out,*)"--gnu-mtrace               Enable mtrace (requires GNU and clib)."
       write(std_out,*)"--ieee-halt                Halt the code if one of the *usual* IEEE exceptions is raised."
       write(std_out,*)"--ieee-signal              Signal the occurrence of the *usual* IEEE exceptions."
       ! Multibinit
       write(std_out,*)"--F03                      Run F03 mode (Multibinit only)."
     end if
     args%exit = args%exit + 1

   else if (arg == "--verbose") then
     verbose = .True.

   else
     if (firstchar(arg, "-")) then
       ABI_WARNING("Unsupported option: "//trim(arg))
       args%exit = args%exit + 1
     else
       continue
     end if
   end if
 end do

 if (ntasks_per_node /= -1 .or. memb_per_node /= -one) then
   ! Set mem_per_cpu from node info.
   ABI_CHECK(ntasks_per_node /= -1, "`mem-per-node` requires `ntasks-per-node`")
   ABI_CHECK(memb_per_node /= -one, "`ntasks-per-node` requires `mem-per-node`")
   ABI_CHECK(memb_per_cpu == -one, "`mem-per-cpu` and `mem-per-node` are mutually exclusive!")
   call set_mem_per_cpu_mb(memb_per_node / ntasks_per_node)
 end if

#endif

end function args_parser
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
 character(len=*),intent(in) :: arg,string
! *************************************************************************

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
 character(len=*),intent(in) :: arg,optname
 logical,optional,intent(in) :: default
! *************************************************************************

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
   ABI_ERROR("Aborting now")
 end select

end function parse_yesno
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/parse_slurm_mem
!! NAME
!!  parse_slurm_mem
!!
!! FUNCTION
!!  Parse `arg` string with memory given in Slurm syntax. Return value in Mb.
!!  From https://slurm.schedmd.com/sbatch.html
!!
!!  --mem=<size>[units]
!!
!!  Default units are megabytes. Different units can be specified using the suffix [K|M|G|T].
!!
!!  For a list of slurm env variables that can be used to pass options to Abinit via the submission script, see:
!!  https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/SLURM-environment-variables.html
!!
!! SOURCE

real(dp) function parse_slurm_mem(arg, optname) result(mem_mb)

!Arguments ------------------------------------
 character(len=*),intent(in) :: arg,optname

!Local variables-------------------------------
 integer :: istop, istat
 real(dp) :: fact
 character(len=500) :: iomsg
! *************************************************************************

 fact = one
 istop = find_and_select(arg, &
                         ["K", "M", "G", "T"], &
                         [one/1024._dp, one, 1024._dp, 1024._dp ** 2], fact, iomsg, default=one)

 ABI_CHECK(istop /= -1, iomsg)
 istop = merge(len_trim(arg), istop - 1, istop == 0)

 read(arg(len(optname) + 1: istop), *, iostat=istat, iomsg=iomsg) mem_mb
 ABI_CHECK(istat == 0, sjoin("Invalid syntax for memory string:", arg, ch10, "iomsg", iomsg))
 ABI_CHECK(mem_mb > zero, "mem_mb must be positive!")
 mem_mb = mem_mb * fact

end function parse_slurm_mem
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
!!  Parse scalar string argument from command line. Return exit code.
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
 character(len=*),intent(out) :: argval, msg
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

!----------------------------------------------------------------------

!!****f* m_argparse/get_arg_bool
!! NAME
!!  get_arg_bool
!!
!! FUNCTION
!!  Parse scalar boolean argument from command line. Return exit code.
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

integer function get_arg_bool(argname, argval, msg, default, exclude) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 logical,intent(out) :: argval
 character(len=*),intent(out) :: msg
 logical,optional,intent(in) :: default
 character(len=*),optional,intent(in) :: exclude

!Local variables-------------------------------
 integer :: ii
 logical :: found_argname, found_excl
 character(len=500) :: arg

! *************************************************************************

 ierr = 0; msg = ""; if (present(default)) argval = default
 found_argname = .False.; found_excl = .False.
 argval = .False.

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (present(exclude)) then
     if (arg == "--" // trim(exclude)) found_excl = .True.
   end if
   if (begins_with(arg, "--" // trim(argname))) then
     argval = parse_yesno(arg, "--" // trim(argname), default=.True.)
     found_argname = .True.
   end if
 end do

 if (ierr /= 0) msg = sjoin("Error while reading argument: ", argname, ch10, msg)
 if (found_argname .and. found_excl) then
   ierr = ierr + 1; msg = sjoin("Variables", argname, "and", exclude, "are mutually exclusive", ch10, msg)
 end if

end function get_arg_bool
!!***

!----------------------------------------------------------------------

!!****f* m_argparse/get_start_step_num
!! NAME
!!  get_start_step_num
!!
!! FUNCTION
!!  Parse string from command line in the format "start:step:num" defining an arithmetic progression.
!!  Return exit code.
!!
!! INPUTS
!!  argname= Argument name
!!  [default]= Default value
!!  [exclude]= argname and exclude are mutually exclusive.
!!
!! OUTPUT
!!   ilist= [start, step, num]
!!   msg= Error message
!!
!! SOURCE

integer function get_start_step_num(argname, ilist, msg, default, exclude) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: argname
 integer,intent(out) :: ilist(3)
 character(len=*),intent(out) :: msg
 integer,optional,intent(in) :: default(3)
 character(len=*),optional,intent(in) :: exclude

!Local variables-------------------------------
 integer :: ii, jj
 character(len=500) :: str

! *************************************************************************

 if (present(exclude)) then
   ierr = get_arg_str(argname, str, msg, default="", exclude=exclude)
 else
   ierr = get_arg_str(argname, str, msg, default="")
 end if
 if (ierr /= 0) return

 if (len_trim(str) == 0) then
   if (present(default)) then
     ilist = default
   else
     ierr = ierr + 1; msg = sjoin("Variables", argname, "is not found and default is not given")
   end if
   return
 end if

 ! We got a non-empty string. Let's parse it.
 ii = index(str, ":")
 if (ii <= 1) then
   msg = sjoin("Cannot find first `:` in string:", str)
   ierr = ierr + 1; return
 end if
 ilist(1) = atoi(str(1:ii-1))

 jj = index(str(ii+1:), ":")
 if (jj == 0) then
   msg = sjoin("Cannot find second `:` in string:", str)
   ierr = ierr + 1; return
 end if

 ilist(2) = atoi(str(ii+1: jj+ii-1))
 ilist(3) = atoi(str(jj+ii+1:))
 !print *, "ilist:", ilist

end function get_start_step_num
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
!!
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
!! SOURCE

subroutine parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)

!Arguments ------------------------------------
 integer,intent(out) :: kptopt, nshiftk, chksymbreak
 integer,intent(out) :: kptrlatt(3,3)
 real(dp),allocatable,intent(out) :: shiftk(:,:)

!Local variables-------------------------------
 integer :: ii, lenr, ierr
 character(len=500) :: msg
 integer :: ivec9(9), ngkpt(3)
 real(dp) :: my_shiftk(3 * MAX_NSHIFTK)

! *************************************************************************

 ABI_CHECK(get_arg("kptopt", kptopt, msg, default=1) == 0, msg)
 ABI_CHECK(get_arg("chksymbreak", chksymbreak, msg, default=1) == 0, msg)

 ierr = get_arg_list("ngkpt", ngkpt, lenr, msg, exclude="kptrlatt", want_len=3)
 if (ierr == 0) then
 !if (lenr == 3) then
   kptrlatt = 0
   do ii=1,3
     kptrlatt(ii, ii) = ngkpt(ii)
   end do
 else
   ABI_CHECK(get_arg_list("kptrlatt", ivec9, lenr, msg, exclude="ngkpt", want_len=9) == 0, msg)
   ABI_CHECK(lenr == 9, "Expecting 9 values for kptrlatt")
   kptrlatt = transpose(reshape(ivec9, [3, 3]))
 end if

 ! Init default
 ABI_CHECK(get_arg_list("shiftk", my_shiftk, lenr, msg) == 0, msg)
 if (lenr /= 0) then
   ABI_CHECK(mod(lenr, 3) == 0, "Expecting 3 * nshift array")
   nshiftk = lenr / 3
   ABI_MALLOC(shiftk, (3, nshiftk))
   shiftk = reshape(my_shiftk(1:lenr), [3, nshiftk])
 else
   nshiftk = 1
   ABI_CALLOC(shiftk, (3, nshiftk))
   !shiftk(:, 1) = [half, half, half]
 end if
 !write(std_out, *)"kptopt = ", kptopt, ", chksymbreak = ", chksymbreak, ", nshiftk = ", nshiftk, ", kptrlatt = ", kptrlatt

end subroutine parse_kargs
!!***

end module m_argparse
!!***

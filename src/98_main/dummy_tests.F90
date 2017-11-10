!!****p* ABINIT/dummy_tests
!! NAME
!! dummy_test
!!
!! FUNCTION
!! This code is supposed to be compiled, and generate warnings when it does not fulfills the abirules.
!! This is a way to check that the testing capabilities are not lost when the test farm is modified ...
!!
!! COPYRIGHT
!! Copyright (C) 2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! The heading of the constituted database is read,
!! then the heading of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!! This process can be iterated.
!! The whole database will be stored in
!! central memory. One could introduce a third mode in which
!! only the temporary DDB is in central memory, while the
!! input DDB is read twice : first to make a table of blocks,
!! counting the final number of blocks, and second to merge
!! the two DDBs. This would save memory.
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,ddb_hdr_free
!!      ddb_hdr_open_read,get_command_argument,herald,mblktyp1,mblktyp5,timein
!!      wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program dummy_tests

 use defs_basis
 use m_build_info
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_dummy_tests

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dummy_tests'
 use interfaces_10_dumpinfo
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Local variables-------------------------------
!scalars
 integer :: comm,dummy_out,my_rank
 integer :: unused_arg,unused_variable,used_arg,used_variable
 real(dp) :: tcpu,tcpui,twall,twalli
 character(len=24) :: codename
!arrays
 real(dp) :: tsec(2)
 character(len=10) :: dummy_string

!******************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

 codename='DUMMY_TESTS'//repeat(' ',13)
 call herald(codename,abinit_version,std_out)

!**********************************************************************

 used_variable=1
 used_arg=used_variable
 call test_unused_arg(used_arg,unused_arg)
 write(std_out,'(a,i4)')' dummy_tests : used_arg is ',used_arg

 call test_same_actual_arg(dummy_out,dummy_out,used_arg)

 dummy_string="This is too long !"

 call test_dummy(dummy_out,used_arg)

!**********************************************************************

 call timein(tcpu,twall)

 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 write(std_out, '(a,a,a,f13.1,a,f13.1)' ) &
& '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 call wrtout(std_out,'+dummy_tests : the run completed successfully ','COLL', do_flush=.True.)

 call abinit_doctor("__dummy_tests")

 100 call xmpi_end()

 end program dummy_tests
!!***

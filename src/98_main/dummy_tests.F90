!!****p* ABINIT/dummy_tests
!! NAME
!! dummy_tests
!!
!! FUNCTION
!! This code is supposed to be compiled, and generate warnings when it does not fulfills the abirules.
!! This is a way to check that the testing capabilities are not lost when the test farm is modified ...
!!
!! COPYRIGHT
!! Copyright (C) 2017-2019 ABINIT group (XG)
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
!!      abi_io_redirect,abimem_init,test_dummy,test_same_actual_arg
!!      test_unused_arg,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program dummy_tests

 use defs_basis
 use m_build_info
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dummy_tests
 implicit none

!Local variables-------------------------------
!scalars
 integer :: comm,dummy_out,my_rank
 integer :: unused_arg,unused_variable,used_arg,used_variable
!arrays
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

!**********************************************************************

 used_variable=1
 used_arg=used_variable
 call test_unused_arg(used_arg,unused_arg)
 write(std_out,'(a,i4)')' dummy_tests : used_arg is ',used_arg

 call test_same_actual_arg(dummy_out,dummy_out,used_arg)

 dummy_string="This is too long !"

 call test_dummy(dummy_out,used_arg)

!**********************************************************************

 100 call xmpi_end()

 end program dummy_tests
!!***

!!****p* ABINIT/mrgddb
!! NAME
!! mrgddb
!!
!! FUNCTION
!! This code merges the derivative databases.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, SP)
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
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,flush_unit,herald,mrgddb_init,initmpi_seq
!!      inprep8,mblktyp1,mblktyp5,timein,wrtout,xmpi_end,xmpi_init
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,ddb_getdims,destroy_mpi_enreg
!!      flush_unit,get_command_argument,herald,initmpi_seq,mblktyp1,mblktyp5
!!      mrgddb_init,timein,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program mrgddb

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_io_tools, only : flush_unit
 use m_mpinfo,   only : destroy_mpi_enreg
 use m_ddb,      only : ddb_getdims, DDB_VERSION

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrgddb'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Local variables-------------------------------
! Set array dimensions
!  mddb=maximum number of databases (cannot be made dynamic)
 integer,parameter :: mddb=5000,ddbun=2
 integer :: chkopt,dummy,dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7
 integer :: iddb,ii,mblktyp,mblktyptmp,nddb
!  msym=maximum number of symmetry elements in space group
 integer :: msym,comm,my_rank
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 character(len=24) :: arg
 character(len=24) :: codename
 character(len=fnlen) :: dscrpt
 character(len=fnlen) :: filnam(mddb+1)
 type(MPI_type) :: mpi_enreg

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()
 comm = xmpi_world
 my_rank = xmpi_comm_rank(comm)

 ABI_CHECK(xmpi_comm_size(comm)==1, "mrgddb not programmed for parallel execution")

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

!Default for sequential use
!Other values of mpi_enreg are dataset dependent, and should NOT be initialized inside mrgddb.F90.
 call initmpi_seq(mpi_enreg)

 codename='MRGDDB'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

!Initialise the code : write heading,
!read names of files, operating mode (also check its value),
!and short description of new database.
 call mrgddb_init(dscrpt,filnam,mddb,nddb)

!Detect whether to run mrgddb in a strict way or not (making consistency checks
!or not) from command line
 chkopt = 1
 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "--nostrict") then
     chkopt = 0
   end if
 end do

!Evaluate the mblktyp of the databases
 mblktyptmp=1
 do iddb=1,nddb
   call ddb_getdims(dummy,filnam(iddb+1),dummy1,dummy2,mblktyp,&
&   msym,dummy3,dummy4,dummy5,dummy6,ddbun,dummy7,DDB_VERSION,comm)

   if(mblktyp > mblktyptmp) mblktyptmp = mblktyp
 end do

 mblktyp = mblktyptmp
!write(std_out,*),'mblktyp',mblktyp

 if (mblktyp==5) then
!  Memory optimized routine
   call mblktyp5(chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,DDB_VERSION)
 else
!  Speed optimized routine
   call mblktyp1(chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,DDB_VERSION)
 end if

!**********************************************************************

 call timein(tcpu,twall)

 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 write(std_out, '(a,a,a,f13.1,a,f13.1)' ) &
& '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 call wrtout(std_out,'+mrgddb : the run completed successfully ','COLL')

 call destroy_mpi_enreg(mpi_enreg)

 call flush_unit(std_out)

 call abinit_doctor("__mrgddb")

 call xmpi_end()

 end program mrgddb
!!***


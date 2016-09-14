!{\src2tex{textfont=tt}}
!!****p* ABINIT/epigene
!! NAME
!! epigene
!!
!! FUNCTION
!! Main routine the Effective PotentIel GENErator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2015 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirec,flush_unit,herald,int2char4,
!!      init10,timein,xmpi_bcast,wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program epigene

 use defs_basis
 use defs_abitypes
 use m_abimover
 use m_build_info
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_effective_potential
 use m_epigene_dataset
 use m_effective_potential_file
 use m_libxml

 use m_io_tools,       only : get_unit, flush_unit
 use m_fstrings,       only : int2char4
 use m_time ,          only : asctime

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'epigene'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_78_effpot
 use interfaces_95_drive
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
! Set array dimensions
 integer,parameter :: ddbun=2,master=0 ! FIXME: these should not be reserved unit numbers!
 integer :: comm,ii,ierr,lenstr
 integer :: natom,nph1l,nrpt,ntypat,nproc,my_rank
 integer :: option
 logical :: iam_master
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 character(len=24) :: codename,start_datetime
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(15),tmpfilename,name
 character(len=500) :: message
 type(epigene_dataset_type) :: inp
 type(effective_potential_type),target :: reference_effective_potential
!TEST_AM
 real(dp),allocatable :: dynmat(:,:,:,:,:)
!******************************************************************
!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()
 comm = xmpi_world

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Initialisation of the timing
 call timein(tcpui,twalli)

 if (iam_master) then
   codename='EPIGENE'//repeat(' ',17)
   call herald(codename,abinit_version,std_out)
 end if

 start_datetime = asctime()

!Initialise the code : write heading, and read names of files.
 if (iam_master) then
   call init10(filnam)
 end if
 call xmpi_bcast (filnam, master, comm, ierr)

!make log file for non-master procs
! if (.not. iam_master.and.prt_log) then
!   std_out = get_unit()
!   call int2char4(my_rank, procstr)
!   ABI_CHECK((procstr(1:1)/='#'),'Bug: string length too short!')
!   tmpfilename = trim(filnam(2)) // "_LOG_P" // trim(procstr)
!   open (unit=std_out, file=tmpfilename)
! end if

!******************************************************************

 call timein(tcpu,twall)

 write(message, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

! Open output files and ab_out (might change its name if needed)
! MJV 1/2010 : now output file is open, but filnam(2) continues unmodified
! so the other output files are overwritten instead of accumulating.
 if (iam_master) then
   tmpfilename = filnam(2)
   call isfile(tmpfilename,'new')
   open (unit=ab_out,file=tmpfilename,form='formatted',status='new')
   rewind (unit=ab_out)
   call herald(codename,abinit_version,ab_out)
 else
   ab_out = dev_null
 end if

 write(message, '(a,(80a),a,a)' ) ch10,&
&  ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!To automate a maximum calculation, epigine reads the number of atoms 
!in the file (ddb or xml). If DDB file is present in input, the ifc calculation
!will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10
 write(message, '(5a)' )' Read the information in the reference structure in ',ch10,&
&    trim(filnam(3)),ch10,' to initialize the epigene input'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 call effective_potential_file_getDim(filnam(3),natom,ntypat,nph1l,nrpt,comm)

!Read the input file, and store the information in a long string of characters
!strlen from defs_basis module
 option=1
 if (iam_master) then
   call instrng (filnam(1),lenstr,option,strlen,string)
   !To make case-insensitive, map characters to upper case:
   call inupper(string(1:lenstr))
 end if

 call xmpi_bcast(string,master, comm, ierr)
 call xmpi_bcast(lenstr,master, comm, ierr)

!Read the input file
 call invars10(inp,lenstr,natom,string)

 if (iam_master) then
!  Echo the inputs to console and main output file
   call outvars_epigene(inp,std_out)
   call outvars_epigene(inp,ab_out)
 end if

 if (iam_master) then

! First step: Treat the reference structure 
!**********************************************************************
   call effective_potential_file_read(filnam(3),reference_effective_potential,inp,comm)
!**********************************************************************

!**********************************************************************
!Second step: Compute the third order derivative with finite differences
   if (inp%prt_3rd > 0) then 
     call strain_phonon_coupling(reference_effective_potential,filnam,inp,comm)
   end if
!**********************************************************************

!**********************************************************************
!  Print the effective potential
   if(inp%prt_effpot<=-1.or.inp%prt_effpot>=3) then
     select case(inp%prt_effpot)
     case (-1)  
       name = "system.xml"
       call effective_potential_writeXML(reference_effective_potential,1,filename=name)
     case(-2)
       name = "system.nc"
       call effective_potential_writeNETCDF(reference_effective_potential,1,filename=name)
     case (3)  
       name = "system.xml"
       call effective_potential_writeXML(reference_effective_potential,1,filename=name)
     end select
   end if

!*********************************************************************
!TEST_AM SECTION
! Print the Phonon dos/spectrum
!   if(inp%prt_phfrq > 0) then
!     call effective_potential_printPDOS(reference_effective_potential,filnam(2),&
!&           inp%n_cell,inp%nph1l,inp%prt_phfrq,inp%qph1l)
!   end if

! Just for TEST 
!Intialisation of the effective potential type
!   call effective_potential_file_read(filnam(4),reference_effective_potential,inp,comm)
!   name = "test.xml"
!   call effective_potential_writeXML(reference_effective_potential,1,filename=name)
! just for TEST

!TEST_AM
   if(inp%prt_phfrq > 0) then
     ABI_ALLOCATE(dynmat,(2,3,reference_effective_potential%supercell%natom_supercell,3,reference_effective_potential%supercell%natom_supercell))

     call effective_potential_effpot2dynmat(dynmat,inp%delta_df,reference_effective_potential,&
&                                           reference_effective_potential%supercell%natom_supercell,&
&                                           int(reference_effective_potential%supercell%qphon),3)

     ABI_DEALLOCATE(dynmat)
   end if
!TEST_AM
!**********************************************************************


!TEST_AM
   if(inp%monte_carlo>=1) then   
!    Compute the monte carlo, molecular dynamics of compute specific energy 
     call mover_effpot(inp,filnam,reference_effective_potential,comm)
   end if
!TEST_AM


!********************************************************************** 
!Free the effective_potential 
   call effective_potential_free(reference_effective_potential)
!**********************************************************************

   write(message,'(a,a,a,(80a))') ch10,('=',ii=1,80),ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   
 end if

 call timein(tcpu,twall)
 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 write(message, '(a,i4,a,f13.1,a,f13.1)' )' Proc.',my_rank,' individual time (sec): cpu=',&
&                 tsec(1),'  wall=',tsec(2)
 call wrtout(std_out,message,"COLL")

 if (iam_master) then
   write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' )'-',ch10,&
&   '- Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 end if

 call xmpi_sum(tsec,comm,ierr)

 write(message, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
& ('=',ii=1,80),ch10,ch10,&
& '+Total cpu time',tsec(1),&
& '  and wall time',tsec(2),' sec',ch10,ch10,&
& ' epigene : the run completed succesfully.'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
  
 call flush_unit(ab_out)
 call flush_unit(std_out)

 if (iam_master) close(ab_out)

 call xmpi_end()
 
end program epigene
!!***

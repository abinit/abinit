!{\src2tex{textfont=tt}}
!!****p* ABINIT/multibinit
!! NAME
!! multibinit
!!
!! FUNCTION
!! Main routine MULTIBINIT.
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
!!      abi_io_redirect,abimem_init,abinit_doctor,compute_anharmonics,
!!      effective_potential_free,effective_potential_file_getDimSystem,effective_potential_file_read,
!!      effective_potential_writeNETCDF,effective_potential_writeXML,flush_unit,herald
!!      init10,instrng,invars10,inupper, isfile,mover_effpot,multibinit_dtset_fre
!!      outvars_multibinit,timein,xmpi_bcast,xmpi_end,xmpi_init,xmpi_sum,wrtout
!!      abi_io_redirec,flush_unit,herald,int2char4,
!!      init10,timein,xmpi_bcast,wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program multibinit

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_xomp
 use m_profiling_abi
 use m_errors
 use m_effective_potential
 use m_multibinit_dataset
 use m_effective_potential_file

 use m_io_tools,   only : get_unit, flush_unit,open_file
 use m_fstrings,   only : int2char4
 use m_time ,      only : asctime

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multibinit'
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
 integer :: comm,filetype,ii,ierr,lenstr
 integer :: natom,nph1l,nrpt,ntypat,nproc,my_rank
 integer :: option
 logical :: iam_master
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 character(len=24) :: codename,start_datetime
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(15),tmpfilename,name
 character(len=500) :: message
 type(multibinit_dataset_type) :: inp
 type(effective_potential_type),target :: reference_effective_potential
!TEST_AM
! integer :: natom_sp
! real(dp),allocatable :: dynmat(:,:,:,:,:)
!TEST_AM
!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

!MPI variables
 comm = xmpi_world
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Initialize memory profiling if it is activated !if a full abimem.mocc report is desired, 
!set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

!Initialisation of the timing
 call timein(tcpui,twalli)

 if (iam_master) then
   codename='MULTIBINIT'//repeat(' ',17)
   call herald(codename,abinit_version,std_out)
 end if

 start_datetime = asctime()

!Print the number of cpus in log file
 write(message,'(a,i5,a)') '-  nproc =',nproc,ch10
 call wrtout(std_out,message,'COLL')

!Initialise the code : write heading, and read names of files.
 call init10(filnam,comm)

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
   if (open_file(tmpfilename,message,unit=ab_out,form="formatted",status="new",&
     &              action="write") /= 0) then
     MSG_ERROR(message)
   end if
!   call open_file(unit=ab_out,file=tmpfilename,form='formatted',status='new')
   rewind (unit=ab_out)
   call herald(codename,abinit_version,ab_out)
!  Print the number of cpus in output
   write(message,'(a,i5,a)') '-  nproc =',nproc
   call wrtout(ab_out,message,'COLL')
 else
   ab_out = dev_null
 end if

 write(message, '(a,(80a),a)' ) ch10,&
&  ('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')


!To automate a maximum calculation, multibinit reads the number of atoms 
!in the file (ddb or xml). If DDB file is present in input, the ifc calculation
!will be initilaze array to the maximum of atoms (natifc=natom,atifc=1,natom...) in invars10
 write(message, '(6a)' )' Read the information in the reference structure in ',ch10,&
&    '-',trim(filnam(3)),ch10,' to initialize the multibinit input'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 call effective_potential_file_getDimSystem(filnam(3),natom,ntypat,nph1l,nrpt,comm)
 
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
   call outvars_multibinit(inp,std_out)
   call outvars_multibinit(inp,ab_out)
 end if

! Read and treat the reference structure 
!****************************************************************************************
!Read the harmonics parts
 call effective_potential_file_read(filnam(3),reference_effective_potential,inp,comm)
!Read the coefficient from fit
 if(filnam(4)/='')then
   call effective_potential_file_getType(filnam(4),filetype)
   if(filetype==3) then
     call effective_potential_file_read(filnam(4),reference_effective_potential,inp,comm)
   end if
 else
   write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
&                       'There is no file for the coefficients from polynomial fitting'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
!****************************************************************************************

! Compute the third order derivative with finite differences
!****************************************************************************************
 if (inp%prt_3rd > 0) then 
   call compute_anharmonics(reference_effective_potential,filnam,inp,comm)
 end if
!****************************************************************************************

!If needed, fit the anharmonic part
!****************************************************************************************
 if (inp%ncoeff > 0) then
!   call fit_polynomial_coeff_init
!   call fit_polynomial_coeff_init(reference_effective_potential%,filnam,inp,comm)
 end if
!****************************************************************************************

!****************************************************************************************
!Print the effective potential system + coefficients (only master CPU)
 if(iam_master.and.(inp%prt_effpot<=-1.or.inp%prt_effpot>=3)) then
   select case(inp%prt_effpot)
   case (-1)  
     name = "system.xml"
     call effective_potential_writeXML(reference_effective_potential,-1,filename=name)
   case(-2)
     name = "system.nc"
     call effective_potential_writeNETCDF(reference_effective_potential,1,filename=name)
   case (3)  
     name = "system.xml"
     call effective_potential_writeXML(reference_effective_potential,1,filename=name)
   end select
 end if
!****************************************************************************************

!TEST_AM SECTION
! Print the Phonon dos/spectrum
!   if(inp%prt_phfrq > 0) then
!     call effective_potential_printPDOS(reference_effective_potential,filnam(2),&
!&           inp%n_cell,inp%nph1l,inp%prt_phfrq,inp%qph1l)
!   end if

!Intialisation of the effective potential type
!   call effective_potential_file_read(filnam(4),reference_effective_potential,inp,comm)
!   name = "test.xml"
!   call effective_potential_writeXML(reference_effective_potential,1,filename=name)
! just for TEST
!   if(inp%prt_phfrq > 0) then
!     natom_sp = reference_effective_potential%supercell%natom_supercell
!     ABI_ALLOCATE(dynmat,(2,3,natom_sp,3,natom_sp))
!     call effective_potential_effpot2dynmat(dynmat,inp%delta_df,reference_effective_potential,&
!&                                           reference_effective_potential%supercell%natom_supercell,&
!&                                           int(reference_effective_potential%supercell%qphon),3)
!
!     ABI_DEALLOCATE(dynmat)
!   end if
!TEST_AM SECTION


! Compute the monte carlo, molecular dynamics of compute specific energy 
!****************************************************************************************
 if(inp%dynamics>=1) then
   call mover_effpot(inp,filnam,reference_effective_potential,comm)
 end if
!****************************************************************************************    


!Free the effective_potential and dataset
!**************************************************************************************** 
 call effective_potential_free(reference_effective_potential)
 call multibinit_dtset_free(inp)
!****************************************************************************************

 write(message,'(a,a,a,(80a))') ch10,('=',ii=1,80),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

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
& ' multibinit : the run completed succesfully.'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 if (iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
   write(std_out,"(a)")"--- !FinalSummary"
   write(std_out,"(a)")"program: multibinit"
   write(std_out,"(2a)")"version: ",trim(abinit_version)
   write(std_out,"(2a)")"start_datetime: ",start_datetime
   write(std_out,"(2a)")"end_datetime: ",asctime()
   write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
   write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
   write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
   write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out,"(a,i0)")"num_warnings: ",nwarning
   !write(std_out,"(a,i0)")"num_comments: ",ncomment
   write(std_out,"(a)")"..."
   call flush_unit(std_out)
 end if

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__multibinit")
  
 call flush_unit(ab_out)
 call flush_unit(std_out)

 if (iam_master) close(ab_out)

 call xmpi_end()
 
end program multibinit
!!***

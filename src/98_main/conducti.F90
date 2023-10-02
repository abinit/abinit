!!****p* ABINIT/conducti
!! NAME
!! conducti
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2022 ABINIT group (FJ,SMazevet)
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
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program conducti

 use defs_basis
 use m_xmpi
 use m_errors
 use m_abicore
#if defined HAVE_MPI2
 use mpi
#endif
 use m_conducti

 use m_io_tools,  only : open_file
 use m_time,      only : timein
 use m_fstrings,  only : sjoin, itoa
 use m_nctk,      only : nctk_test_mpiio
 use m_paw_optics,only : linear_optics_paw

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: incpaw,nproc,comm,inunt,my_rank,mpierr
 real(dp) :: tcpu,tcpui,twall,twalli
 character(len=fnlen) :: filnam,filnam_out
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)

! *********************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Start, MPI init, ...
 call xmpi_init()
 call timein(tcpui,twalli)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

!Test if the netcdf library supports MPI-IO
 call nctk_test_mpiio(print_warning=.false.)

!Extract MPI parallelization data
 comm = xmpi_world
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)


!Read some input data
 if (my_rank==master) then

!  Read data file name
   write(std_out,'(a)')' Please, give the name of the data file ...'
   read(std_in, '(a)') filnam
   write(std_out,'(2a)')' The name of the data file is: ',trim(filnam)

!  Read type of calculation
   if (open_file(filnam,msg,newunit=inunt,form='formatted')==0) then
     rewind(inunt)
     read(inunt,*) incpaw
     close(inunt)
   else
     ABI_ERROR('Error opening data file!')
   end if

!  Read output file name
   write(std_out,'(a)')' Give the name of the output file ...'
   read(std_in, '(a)') filnam_out
   write(std_out,'(2a)')' The name of the output file is: ',filnam_out

 end if

!Broadcast input data
 call xmpi_bcast(incpaw,master,comm,mpierr)
 call xmpi_bcast(filnam,master,comm,mpierr)
 call xmpi_bcast(filnam_out,master,comm,mpierr)

!Call main routine
 if (incpaw==1) then
   if (my_rank==master) then
     call conducti_nc(filnam,filnam_out)
   end if
 elseif (incpaw==2) then
   call conducti_paw(filnam,filnam_out)
 elseif (incpaw==3) then
   if (my_rank==master) then
     call linear_optics_paw(filnam,filnam_out)
   end if
 elseif (incpaw==4) then
   call conducti_paw(filnam,filnam_out)
   call conducti_paw_core(filnam,filnam_out,with_absorption=.true.,with_emissivity=.true.)
 elseif (incpaw==5) then
   call conducti_paw_core(filnam,filnam_out,with_absorption=.true.)
 elseif (incpaw==6) then
   call conducti_paw_core(filnam,filnam_out,with_absorption=.true.,with_emissivity=.true.)
 elseif (incpaw==42) then
   call conducti_paw(filnam,filnam_out,varocc=1)
 else
   ABI_ERROR(sjoin("Wrong incpaw:", itoa(incpaw)))
 end if

!End, memory cleaning
 call timein(tcpu,twall)
 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli
 if (my_rank==0) then
   write(std_out, '(a,a,a,f13.1,a,f13.1)' ) &
&   '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 end if

 call abinit_doctor("__conducti")

 call xmpi_end()

 end program conducti
!!***

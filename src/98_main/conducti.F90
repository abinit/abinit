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
!! Copyright (C) 2006-2019 ABINIT group (FJ,SMazevet)
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
!!      abi_io_redirect,abimem_init,abinit_doctor,conducti_nc,conducti_paw
!!      conducti_paw_core,destroy_mpi_enreg,emispec_paw,linear_optics_paw
!!      timein,xmpi_end,xmpi_init
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

 use defs_abitypes, only : MPI_type
 use m_io_tools,  only : open_file
 use m_time,      only : timein
 use m_fstrings,  only : sjoin, itoa
 use m_mpinfo,    only : destroy_mpi_enreg
 use m_paw_optics,only : linear_optics_paw

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 integer :: incpaw,nproc,comm,inunt,my_rank
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 character(len=fnlen) :: filnam,filnam_out
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg_seq ! needed for calling rwwf

! *********************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 call xmpi_init()

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

 comm = xmpi_world
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 if ( nproc > 1) then
   MSG_ERROR("conducti is not parallelized: run with one processor.")
 end if

!Read data file name
 write(std_out,'(a)')' Please, give the name of the data file ...'
 read(std_in, '(a)')filnam
 write(std_out,'(a)')' The name of the data file is :',trim(filnam)

 if (open_file(filnam,msg,newunit=inunt,form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if
 rewind(inunt)
 read(inunt,*) incpaw
 close(inunt)

 write(std_out,'(a)')' Give the name of the output file ...'
 read(std_in, '(a)') filnam_out
 write(std_out,'(a)')' The name of the output file is :',filnam_out

 if (incpaw==1) then
   call conducti_nc(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==2) then
   call conducti_paw(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==3) then
   call linear_optics_paw(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==4) then
   call conducti_paw(filnam,filnam_out,mpi_enreg_seq)
   call conducti_paw_core(filnam,filnam_out,mpi_enreg_seq)
   call emispec_paw(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==5) then
   call conducti_paw_core(filnam,filnam_out,mpi_enreg_seq)
 elseif (incpaw==6) then
   call emispec_paw(filnam,filnam_out,mpi_enreg_seq)
 else
   MSG_ERROR(sjoin("Wrong incpaw:", itoa(incpaw)))
 end if

 call destroy_mpi_enreg(mpi_enreg_seq)

 call timein(tcpu,twall)

 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli
 write(std_out, '(a,a,a,f13.1,a,f13.1)' ) &
& '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 call abinit_doctor("__conducti")

 call xmpi_end()

 end program conducti
!!***

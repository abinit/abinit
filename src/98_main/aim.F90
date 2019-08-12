!{\src2tex{textfont=tt}}
!!****p* ABINIT/aim
!! NAME
!! aim
!!
!! FUNCTION
!! Main routine for Bader Atom-In-Molecule analysis.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2019 ABINIT group (PCasek,FF,XG)
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
!! WARNING
!! ABINIT rules are not yet followed in the present routine.
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,adini,aim_shutdown,defad,drvaim,herald
!!      inpar,int2char4,timein,xmpi_bcast,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program aim

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_xmpi
 use m_build_info
 use m_errors
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,     only : timein
 use m_io_tools, only : open_file, file_exists
 use m_specialmsg,  only : specialmsg_getcount, herald
 use m_fstrings, only : int2char4
 use m_bader !,    only : adini, drvaim, inpar, defad, aim_shutdown

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: fin,ii,ios,iunt,ierr,nfcfile
! Allow for maximum of 100 fc files
 integer,parameter :: natm=500
 integer :: lenstr,me,nproc,comm
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 character(len=fnlen) :: dnfile,hname,infile,ofile,tmpfilename
 character(len=strlen) :: instr
 character(len=10) :: procstr
 character(len=24) :: codename
 type(aim_dataset_type) :: aim_dtset
 character(len=500) :: msg
 character(len=fnlen) :: fcfile(natm)

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

 comm = xmpi_world

 me=xmpi_comm_rank(comm); nproc=xmpi_comm_size(comm)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

!Initialize the code, master only : write heading, and read names of files.
 if (me == master) then
   read(*,'(a)') infile
   infile = trim(infile)
   read(*,'(a)') dnfile
   dnfile = trim(dnfile)
   read(*,'(a)') ofile
   ofile = trim(ofile)
   do ii=1,natm
     iunt=unt+ii
     read(*,'(a)',iostat=ios) fcfile(ii)
     if (ios /=0) exit
   end do
   nfcfile=ii-1
 end if

!Transfer file names to other procs
 call xmpi_bcast (infile, master, comm, ierr)
 call xmpi_bcast (dnfile, master, comm, ierr)
 call xmpi_bcast (ofile, master, comm, ierr)
 call xmpi_bcast (nfcfile, master, comm, ierr)
 do ii=1,nfcfile
   call xmpi_bcast (fcfile(ii), master, comm, ierr)
 end do

!Prepare initialization of the log and output files
 fin=len_trim(ofile)
 hname(1:fin)=ofile(1:fin)
 hname(fin+1:fin+1)='.'
 untout=14
 codename='AIM   '//repeat(' ',18)

!Open main output file and main log file, then print herald at top of files
 if(me==master)then

   hname(fin+2:fin+4)='out'
   if (open_file(hname(1:fin+4),msg,unit=untout,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if
   rewind (unit=untout)
   call herald(codename,abinit_version,untout)
   call herald(codename,abinit_version,std_out)
 end if

!Open log file for non-master procs, then print herald at top of files
 if (me /= master) then
   call int2char4(me, procstr)
   tmpfilename = hname(1:fin) // "_LOG_P" // trim(procstr)
   !close(std_out)
   if (open_file(tmpfilename, msg, newunit=std_out, status='unknown',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if
   rewind (unit=std_out)
   call herald(codename,abinit_version,std_out)
 end if

 unt=21 ! WARNING : this number is used to define other unit numbers, in init.f
 unt0=9

 unto=6  ! XG020629 use standard IO unit => easier testing . So, should replace everywhere unto by std_out

 untc=11
 unts=12
 untad=19
 untd=17
 untl=18
 unta=15
 untp=13
 untg=20

!OPENING OF THE INPUT FILES

 aim_iomode = IO_MODE_FORTRAN
 if(me==master)then
   if (open_file(infile,msg,unit=unt0,status='old',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   if (file_exists(dnfile)) then
     if (open_file(dnfile,msg,unit=untad,status='old',form='unformatted') /=0) then
       MSG_ERROR(msg)
     end if
   else
     if (file_exists(nctk_ncify(dnfile))) then
       ! Use netcdf-io
       write(std_out,"(3a)")"- File: ",trim(dnfile)," does not exist but found netcdf file with similar name."
       dnfile = nctk_ncify(dnfile)
       aim_iomode = IO_MODE_ETSF
#ifdef HAVE_NETCDF
       NCF_CHECK(nctk_open_read(untad, dnfile, xmpi_comm_self))
#else
       MSG_ERROR("Cannot read netcdf file because netcdf support in Abinit is missing.")
#endif
     else
       MSG_ERROR('Missing data file: '//TRIM(dnfile))
     end if
   end if

   do ii=1,nfcfile
     iunt=unt+ii
     if (open_file(fcfile(ii),msg,unit=iunt,status='old',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
   end do
 end if
 call xmpi_bcast(aim_iomode, master, comm, ierr)

!call dump_config(std_out)

!READING OF THE MAIN INPUT FILE

!Setting the input variables to their default values
 call defad(aim_dtset)

!Reading of the input file -> one string called instr
 if(me==master)then
   call inpar(instr,lenstr)
 end if
 call xmpi_bcast (lenstr, master, comm, ierr)
 call xmpi_bcast (instr(1:lenstr), master, comm, ierr)

!Analysis of the input string, setting of input variables in aim_dtset
 call adini(aim_dtset,instr,lenstr)

!OPENING OF THE OUTPUT FILES
 if(me==master)then

   if (aim_dtset%isurf/=0) hname(fin+2:fin+5)='surf'

   if (aim_dtset%isurf==1) then
     ierr = open_file(hname(1:fin+5),msg,unit=unts,status='unknown',form='formatted')
   elseif (aim_dtset%isurf==-1) then
     ierr = open_file(hname(1:fin+5),msg,unit=unts,status='old',action='read',form='formatted')
   end if
   if (ierr /= 0) then
     MSG_ERROR(msg)
   end if

   if (aim_dtset%crit/=0) hname(fin+2:fin+5)='crit'

   if (aim_dtset%crit>0) then
     ierr = open_file(hname(1:fin+5),msg,unit=untc,status='unknown',form='formatted')
   else if (aim_dtset%crit==-1) then
     ierr = open_file(hname(1:fin+5),msg,unit=untc,status='old',action='read',form='formatted')
   end if
   if (ierr /= 0) then
     MSG_ERROR(msg)
   end if

   if (aim_dtset%denout==1) then
     hname(fin+2:fin+4)='dn1'
     ierr = open_file(hname(1:fin+4),msg,unit=untd,status='unknown',form='formatted')
   elseif (aim_dtset%denout==2) then
     hname(fin+2:fin+4)='dn2'
     ierr = open_file(hname(1:fin+4),msg,unit=untd,status='unknown',form='formatted')
   elseif (aim_dtset%denout==3) then
     hname(fin+2:fin+4)='dn3'
     ierr = open_file(hname(1:fin+4),msg,unit=untd,status='unknown',form='formatted')
   elseif (aim_dtset%denout==-1) then
     hname(fin+2:fin+4)='dna'
     ierr = open_file(hname(1:fin+4),msg,unit=untd,status='unknown',form='unformatted')
   end if

   if (ierr /= 0) then
     MSG_ERROR(msg)
   end if

   if (aim_dtset%lapout==1) then
     hname(fin+2:fin+4)='lp1'
     ierr = open_file(hname(1:fin+4),msg,unit=untl,status='unknown',form='formatted')
   elseif (aim_dtset%lapout==2) then
     hname(fin+2:fin+4)='lp2'
     ierr = open_file(hname(1:fin+4),msg,unit=untl,status='unknown',form='formatted')
   elseif (aim_dtset%lapout==3) then
     hname(fin+2:fin+4)='lp3'
     ierr = open_file(hname(1:fin+4),msg,unit=untl,status='unknown',form='formatted')
   elseif (aim_dtset%lapout==-1) then
     hname(fin+2:fin+4)='lpa'
     ierr = open_file(hname(1:fin+4),msg,unit=untl,status='unknown',form='unformatted')
   end if

   if (ierr /= 0) then
     MSG_ERROR(msg)
   end if

   if (aim_dtset%gpsurf==1) then
     hname(fin+2:fin+3)='gp'
     if (open_file(hname(1:fin+3),msg,unit=untg,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
   end if

   if (aim_dtset%plden==1) then
     hname(fin+2:fin+4)='pld'
     if (open_file(hname(1:fin+5),msg,unit=untp,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
   end if

 end if

 call timein(tcpu,twall)
 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli
 write(std_out, '(5a,f13.1,a,f13.1)' ) &
& '-',ch10,'- After reading the input file and opening the output files ',ch10,&
& '- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

!MAIN DRIVER OF THE ANALYSIS

 write(std_out,'(a,a,a,a,i4)' )char(10),&
& ' aim : read density file ',trim(dnfile),' from unit number ',untad

 call drvaim(aim_dtset,tcpui,twalli)

!THE TOTAL TIME NEEDED
 call timein(tcpu,twall)
 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli
 write(std_out, '(a,a,a,f13.1,a,f13.1)' ) &
& '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 write(std_out,'(/," Time needed (seconds) - total, CP analyse, SURF determination:",/,/,"-         ",3F16.8)')&
& tsec(2),ttcp,ttsrf


 if(me==0)then
   write(untout,*)
   write(untout,*) "TIME ANALYSIS"
   write(untout,*) "============"
   write(untout,'(/," Time needed (seconds) - total, CP analyse, SURF determination:",/,/,"-         ",3F16.8)') &
&   tsec(1),ttcp,ttsrf
   write(untout,'(a,a,f11.3,a,f11.3,a,a,a,a)') char(10),&
&   '+Total cpu time',tsec(1),&
&   '  and wall time',tsec(2),' sec',char(10),char(10),&
&   ' aim : the run completed succesfully.'
 end if

!CLOSING OF THE FILES
 if(me==0)then
   close(unt0)
   close(untout)
   if (aim_dtset%isurf/=0) close(unts)
   if (aim_dtset%crit/=0) close(untc)
   if (aim_dtset%denout/=0) close(untd)
   if (aim_dtset%lapout/=0) close(untl)
   if (aim_dtset%gpsurf/=0) close(untg)
   if (aim_dtset%plden/=0) close(untp)
 end if

 call aim_shutdown()

 !call abinit_doctor("__aim")

!Eventual cleaning of MPI run
 call xmpi_end()

 end program aim
!!***

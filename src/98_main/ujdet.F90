!{\src2tex{textfont=tt}}
!!****p* ABINIT/ujdet
!! NAME
!! ujdet
!!
!! FUNCTION
!!  Main routine. Determines U from inputfile ujdet.in (reduced inputfile containing mainly the atomic occupancies,
!!  atomic positions, and potential shifts.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! INPUTS
!! OUTPUT
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,destroy_mpi_enreg,herald
!!      initmpi_seq,intagm,isfile,parsefile,pawuj_det,pawuj_free,pawuj_ini
!!      wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program ujdet

 use defs_basis
 use defs_datatypes
 use m_xmpi
 use m_abicore
 use m_build_info
 use m_errors

 use defs_abitypes, only : MPI_type, macro_uj_type
 use m_specialmsg,  only : specialmsg_getcount, herald
 use m_io_tools,    only : open_file
 use m_parser,      only : intagm, parsefile
 use m_mpinfo,      only : destroy_mpi_enreg, initmpi_seq
 use m_dtfil,       only : isfile
 use m_paw_uj,      only : pawuj_ini,pawuj_free,pawuj_det

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter     :: ndtpawuj=4,nwfchr=6,master=0
 integer               :: ii,jdtset,lenstr,marr,ndtset,tread,comm
 integer               :: nat,nproc,my_rank,nspden,macro_uj,pawujat,pawprtvol
 logical               :: iam_master
 character(len=24)     :: codename
 character(len=30)     :: token
 character(len=fnlen)  :: filnam(5)
 character(len=strlen) :: string
 character(len=500)    :: message
 type(MPI_type)        :: mpi_enreg
 real(dp)              :: ures
!arrays
 integer,allocatable   :: idumar1(:),intarr(:),jdtset_(:)
 real(dp),allocatable  :: dpdumar(:),dprarr(:)
 type(macro_uj_type),allocatable   :: dtpawuj(:)

! *********************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI (not used but necessary)
 call xmpi_init()
 comm = xmpi_world
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

!Default for sequential use
 call initmpi_seq(mpi_enreg)

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 if(xmpi_paral==0)then
   write(message,'(3a)')&
&   'In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   'Action : recompile your code with different CPP flags.'
   MSG_ERROR(message)
 end if
#endif

!----------------------------------------------------------------------

!Check that old input file exists (18seqpar/iofn1.F90)
 filnam(1)='ujdet.in'
 call isfile(filnam(1),'old')

 codename='UJDET'//repeat(' ',18)
 filnam(2)='ujdet.out'

!Check that new output file does NOT exist (18seqpar/iofn1.F90)
 if (iam_master) then
   call isfile(filnam(2),'new')
   if (open_file(filnam(2),message,unit=ab_out,form="formatted",status="new") /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unit=ab_out)
!  Print header
   call herald(codename,abinit_version,ab_out)
 end if

!Print header
 call herald(codename,abinit_version,std_out)

!Read the file, stringify it and return the number of datasets. (main/abinit.F90)
 call parsefile(filnam(1), lenstr, ndtset, string, comm)

!DEBUG
!write(std_out,*)'ujdet: trim(string):',trim(adjustl(string))
!write(std_out,*)'ujdet: lenstr: ',lenstr
!write(std_out,*)'ujdet: ndtset: ',ndtset
!ENDDEBUG

 ABI_DATATYPE_ALLOCATE(dtpawuj,(0:ndtpawuj))

 call pawuj_ini(dtpawuj,ndtset)

 dtpawuj(1:ndtset)%ndtset=ndtset

 marr=ndtset
 ABI_ALLOCATE(idumar1,(marr))
 ABI_ALLOCATE(dpdumar,(marr))
 ABI_ALLOCATE(jdtset_,(ndtset))
 jdtset_=(/ ( ii,ii=1,ndtset )/)

!Read integers (main dimensions)
 do jdtset=1,ndtset
!  Integer
   token = 'pawprtvol'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%pawprtvol=idumar1(1)

   token = 'nat'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%nat=idumar1(1)

   token = 'nspden'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%nspden=idumar1(1)

   token = 'macro_uj'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%macro_uj=idumar1(1)

   token = 'pawujat'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%pawujat=idumar1(1)

   token = 'dmatpuopt'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%dmatpuopt=idumar1(1)

!  DEBUG
!  write(std_out,*)'ujdet: read pawujat ; jdtset, idumar1(1)',jdtset, idumar1(1)
!  write(std_out,*)'ujdet: read pawujat ; dtpawuj(1:ndtset)%pawujat ', dtpawuj(1:ndtset)%pawujat
!  END DEBUG

   token = 'pawujopt'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%option=idumar1(1)

!  Real

   token = 'diemix'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%diemix=idumar1(1)

   token = 'mdist'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%mdist=idumar1(1)

   token = 'pawujga'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%pawujga=dpdumar(1)

   token = 'ph0phiint'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%ph0phiint=dpdumar(1)

   token = 'pawrad'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1) dtpawuj(1:ndtset)%pawrad=dpdumar(1)

   token = 'pawujrad'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1) dtpawuj(1:ndtset)%pawujrad=dpdumar(1)

 end do !ndtset

 nat=dtpawuj(1)%nat
 nspden=dtpawuj(1)%nspden
 macro_uj=dtpawuj(1)%macro_uj
 pawujat=dtpawuj(1)%pawujat
 pawprtvol=dtpawuj(1)%pawprtvol

!Minimal tests
 if (.not.all(dtpawuj(1:ndtset)%nat==nat).and.all(dtpawuj(1:ndtset)%pawujat==pawujat)&
& .and.all(dtpawuj(1:ndtset)%macro_uj==macro_uj)) then
   write(message,'(a)')'Problems with nat, pawujat or nspden. Values for datasets differ.'
   call wrtout(std_out,message,'COLL')
   write(message,*) ' ujdet: nat ', all(dtpawuj(1:ndtset)%nat==nat)
   call wrtout(std_out,message,'COLL')
   write(message,*) ' ujdet: pawujat ', all(dtpawuj(1:ndtset)%pawujat==pawujat)
   call wrtout(std_out,message,'COLL')
   write(message,*) ' ujdet: macro_uj', all(dtpawuj(1:ndtset)%macro_uj==macro_uj)
   call wrtout(std_out,message,'COLL')
 end if

 if (pawprtvol>1) then
   write(message,fmt='(a,150i3)')' ujdet: dtpawuj(0:ndtset)%macro_uj ',dtpawuj(0:ndtset)%macro_uj
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150i3)')' ujdet: dtpawuj(0:ndtset)%pawujat ',dtpawuj(0:ndtset)%pawujat
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150i3)')' ujdet: dtpawuj(0:ndtset)%nspden ',dtpawuj(0:ndtset)%nspden
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150i3)')' ujdet: nat*nspden ',nat*nspden
   call wrtout(std_out,message,'COLL')
 end if

!Read arrays

 marr=maxval((/nspden*nat*nat, 3*3 ,nat*3/))
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 do jdtset=0,ndtset
!  DEBUG
!  write(std_out,*)'ujdet 2a, jdtset ',jdtset
!  END DEBUG
   ABI_ALLOCATE(dtpawuj(jdtset)%vsh,(nspden,nat))
   ABI_ALLOCATE(dtpawuj(jdtset)%occ,(nspden,nat))
   ABI_ALLOCATE(dtpawuj(jdtset)%xred,(3,nat))

   dtpawuj(jdtset)%iuj=jdtset
   dtpawuj(jdtset)%vsh=zero
   dtpawuj(jdtset)%occ=zero
   dtpawuj(jdtset)%xred=zero
 end do

 do jdtset=1,ndtset
   dprarr=-1_dp
   intarr=-1

!  Integer arrays

!  scdim, wfchr and rprimd allocated in pawuj_ini
   token = 'scdim'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(jdtset)%scdim(1:3)=intarr(1:3)

   token = 'wfchr'
   call intagm(dprarr,intarr,jdtset,marr,nwfchr,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%wfchr(1:nwfchr)=dprarr(1:nwfchr)

   write(std_out,*)'ujdet: wfchr ',dtpawuj(jdtset)%wfchr

   token = 'wfchr'
   call intagm(dprarr,intarr,jdtset,marr,nwfchr,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%wfchr(1:nwfchr)=dprarr(1:nwfchr)

   write(std_out,*)'ujdet: wfchr ',dtpawuj(jdtset)%wfchr

   token = 'rprimd'
   call intagm(dprarr,intarr,jdtset,marr,3*3,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%rprimd(1:3,1:3)=reshape(dprarr(1:3*3),(/ 3,3/))

   token = 'occ'
   call intagm(dprarr,intarr,jdtset,marr,nspden*nat,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%occ(1:nspden,1:nat)=reshape(dprarr(1:nspden*nat),(/ nspden,nat/))

   token = 'vsh'
   call intagm(dprarr,intarr,jdtset,marr,nat*nspden,string(1:lenstr),token,tread,'ENE')
   if(tread==1) dtpawuj(jdtset)%vsh(1:nspden,1:nat)=reshape(dprarr(1:nspden*nat),(/ nspden,nat/))

   token = 'xred'
   call intagm(dprarr,intarr,jdtset,marr,nat*3,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%xred(1:3,1:nat)=reshape(dprarr(1:3*nat),(/ 3, nat/))

!  DEBUG
!  write(std_out,*)'dtpawuj(',jdtset,')%occ ',dtpawuj(jdtset)%occ,ch10
!  write(std_out,*)'dtpawuj(',jdtset,')%rprimd ',dtpawuj(jdtset)%rprimd,ch10
!  write(std_out,*)'dtpawuj(',jdtset,')%vsh ',dtpawuj(jdtset)%vsh,ch10
!  write(std_out,*)'dtpawuj(',jdtset,')%xred ',dtpawuj(jdtset)%xred,ch10,ch10
!  END DEBUG

 end do ! jdtset

!Not yet really parallelized ...
 if (iam_master) then
   call pawuj_det(dtpawuj,ndtset, filnam(2)//"_UJDET.nc",ures)
 end if

!Close files
 if (iam_master) then
   close(ab_out)
 end if

!Deallocations
 do jdtset=0,ndtset
   call pawuj_free(dtpawuj(jdtset))
 end do
 ABI_DATATYPE_DEALLOCATE(dtpawuj)

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(idumar1)
 ABI_DEALLOCATE(dpdumar)
 ABI_DEALLOCATE(jdtset_)

 call destroy_mpi_enreg(mpi_enreg)

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__ujdet")

 call xmpi_end()

 end program ujdet
!!***

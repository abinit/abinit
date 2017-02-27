!{\src2tex{textfont=tt}}
!!****f* ABINIT/distrb2
!! NAME
!!  distrb2
!!
!! FUNCTION
!!  This routine creates the tabs of repartition of processors
!!  for sharing the jobs on k-points, spins and bands.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2016 ABINIT group (AR,XG,MB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mband = maximum number of bands
!!  nband(nkpt*nsppol) = number of bands per k point, for each spin
!!  nkpt = number of k-points
!!  nproc= number of processors available for this distribution
!!  nsppol = 1 for unpolarized, 2 for polarized
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!   mpi_enreg%proc_distrb(nkpt,mband,nsppol)=number of the processor
!!       that will treat each band in each k point.
!!   mpi_enreg%nproc_kpt is set
!!
!!
!! NOTES
!!  For the time being, the band parallelisation works only
!!  when the number of bands is identical for spin up and spin down
!!  at the same k point. The problem is the most clearly seen
!!  in the kpgio routine, where a different parallel repartition
!!  of k points for spin up and spin down would conflict with the
!!  present computation of k+G sphere, independent of the spin.
!!
!! PARENTS
!!      dfpt_looppert,eig2stern,eig2tot,mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine distrb2(mband,nband,nkpt,nproc,nsppol,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_io_tools,    only : file_exists, open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'distrb2'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mband,nkpt,nproc,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: inb,inb1,ind,ind0,nband_k,proc_max,proc_min
 integer :: iiband,iikpt,iisppol,ikpt_this_proc,nbsteps,nproc_kpt,temp_unit
 integer :: kpt_distrb(nkpt)
 logical,save :: first=.true.,has_file
 character(len=500) :: message

!******************************************************************

 nproc_kpt=mpi_enreg%nproc_kpt
 if (mpi_enreg%paral_pert==1) nproc_kpt=nproc

!Initialization of proc_distrb
 do iisppol=1,nsppol
   do iiband=1,mband
     do iikpt=1,nkpt
       mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=nproc_kpt-1
     end do
   end do
 end do
!That s all for an empty communication space
 if (nproc==0) return

!Some checks
 if (mpi_enreg%paralbd==0) then
!  Check if nkpt and nproc_kpt match
   if(nproc_kpt>nkpt*nsppol) then
!    Too much proc. with respect to nkpt
     write(message,'(a,i0,a,i0,a,i0,2a)')&
&     'nproc_kpt=',nproc_kpt,' >= nkpt=',nkpt,'* nsppol=',nsppol,ch10,&
&     'The number of processors is larger than nkpt*nsppol. This is a waste.'
     MSG_WARNING(message)
   else if(mod(nkpt*nsppol,nproc_kpt)/=0) then
!    nkpt not a multiple of nproc_kpt
     write(message,'(a,i0,a,i0,3a)')&
&     'nkpt*nsppol (', nkpt*nsppol, ') is not a multiple of nproc_kpt (',nproc_kpt, ')', ch10,&
&     'The k-point parallelisation is not efficient.'
     MSG_WARNING(message)
   end if
 end if

!Inquire whether there exist a file containing the processor distribution
 if (first) then
!  Case first time : test file to do
!  Open the file containing the k-point distribution
   first=.false.; has_file = file_exists("kpt_distrb")
 end if

!Initialize the processor distribution, either from a file, or from an algorithm
 if (has_file) then
   if (open_file('kpt_distrb',message,newunit=temp_unit,form='formatted',status='old') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unit=temp_unit)
   if (mpi_enreg%paralbd == 1) then
!    -> read bands distribution
     read(temp_unit,*) mpi_enreg%proc_distrb
   else
     read(temp_unit,*) kpt_distrb
   end if
   close(temp_unit)
   proc_max=0
   proc_min=nproc_kpt
!  -> determine the range of proc. requested
   if (mpi_enreg%paralbd == 1) then
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nband_k = nband(iikpt+(iisppol-1)*nkpt)
         proc_max=maxval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
         proc_min=minval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
       end do
     end do
   else
     proc_max=maxval(kpt_distrb(1:nkpt))
     proc_min=minval(kpt_distrb(1:nkpt))
!    -> fill the tab proc_distrb with kpt_distrb
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nband_k = nband(iikpt+(iisppol-1)*nkpt)
         do iiband=1,nband_k
           mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=kpt_distrb(iikpt)
         end do
       end do
     end do
   end if ! mpi_enreg%paralbd

   if(proc_max>(nproc_kpt-1)) then
!    Too much proc. requested
     write(message, '(a,a,a,i0,a,a,a)' )&
&     'The number of processors mentioned in the kpt_distrb file',ch10,&
&     'must be lower or equal to the actual number of processors =',nproc_kpt-1,ch10,&
&     'Action: change the kpt_distrb file, or increase the','  number of processors.'
     MSG_ERROR(message)
   end if

   if(proc_max/=(nproc_kpt-1)) then
!    Too few proc. used
     write(message, '(a,i0,a,a,a,i0,a,a,a)' )&
&     'Only ',proc_max+1,' processors are used (from kpt_distrb file),',ch10,&
&     'when',nproc_kpt,' processors are available.',ch10,&
&     'Action : adjust number of processors and kpt_distrb file.'
     MSG_ERROR(message)
   end if

   if(proc_min<0) then
     write(message, '(a,a,a)' )&
&     'The number of processors must be bigger than 0 in kpt_distrb file.',ch10,&
&     'Action: modify kpt_distrb file.'
     MSG_ERROR(message)
   end if

 else
!  'kpt_distrb' file does not exist

   if (mpi_enreg%paralbd==1) then

!    No possible band parallelization
     if (nproc<(nkpt*nsppol)) then

!      Does not allow a processor to treat different spins
       ind0=0
       inb1=(nkpt*nsppol)/nproc;if (mod((nkpt*nsppol),nproc)/=0) inb1=inb1+1
       do iikpt=1,nkpt
         nband_k=nband(iikpt)
         ind=ind0/inb1
         do iiband=1,nband_k
           mpi_enreg%proc_distrb(iikpt,iiband,1)=ind
           if (nsppol==2) mpi_enreg%proc_distrb(iikpt,iiband,2)=nproc-ind-1
         end do
         ind0=ind0+1
       end do

!      MT130831 : OLD CODING
!      do iisppol=1,nsppol;do iikpt=1,nkpt
!      ind=(iikpt+(iisppol-1)*nkpt-1)/inb1
!      nband_k=nband(iikpt+(iisppol-1)*nkpt)
!      do iiband=1,nband_k
!      mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind
!      end do;end do;end do
!      MT130831 : END OF OLD CODING

!    Possible band parallelization
     else
!      Does not allow a processor to treat different spins
       ind0=0;inb=nproc/(nkpt*nsppol)
       do iikpt=1,nkpt
         nband_k=nband(iikpt)
         inb1=nband_k/inb;if (mod(nband_k,inb)/=0) inb1=inb1+1
         do iiband=1,nband_k
           ind=(iiband-1)/inb1+ind0
           mpi_enreg%proc_distrb(iikpt,iiband,1)=ind
           if (nsppol==2) mpi_enreg%proc_distrb(iikpt,iiband,2)=nproc-ind-1
         end do
         ind0=ind+1
       end do

!      MT130831 : OLD CODING
!      ind0=0;inb=nproc/(nkpt*nsppol)
!      do iisppol=1,nsppol;do iikpt=1,nkpt
!      nband_k=nband(iikpt+(iisppol-1)*nkpt)
!      inb1=nband_k/inb;if (mod(nband_k,inb)/=0) inb1=inb1+1
!      do iiband=1,nband_k
!      ind=(iiband-1)/inb1+ind0
!      mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind
!      end do
!      ind0=ind+1
!      end do;end do
!      MT130831 : END OF OLD CODING

     end if

!    XG060807 : OLD CODING
!    ind=0
!    do iisppol=1,nsppol;do iikpt=1,nkpt
!    nband_k=nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind/nbsteps
!    ind=ind+1
!    end do;end do;end do
!    XG060807 : END OF OLD CODING

   elseif (mpi_enreg%paralbd==0) then

!    Does not allow a processor to treat different spins
     ind0=0
     nbsteps=(nsppol*nkpt)/nproc_kpt
     if (mod((nsppol*nkpt),nproc_kpt)/=0) nbsteps=nbsteps+1
     do iikpt=1,nkpt
       nband_k=nband(iikpt)
       ind=ind0/nbsteps
       do iiband=1,nband_k
         mpi_enreg%proc_distrb(iikpt,iiband,1)=ind
         if (nsppol==2) mpi_enreg%proc_distrb(iikpt,iiband,2)=nproc_kpt-ind-1
       end do
       ind0=ind0+1
     end do

!    XG060807 : OLD CODING
!    ind=0
!    do iisppol=1,nsppol;do iikpt=1,nkpt
!    nband_k = nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    Distribute k-points homogeneously
!    proc_distrb(iikpt,iiband,iisppol)=mod(iikpt-1,nproc_kpt)
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind/nbsteps
!    end do
!    ind=ind + 1
!    end do;end do
!    XG060807 : END OF OLD CODING

   end if ! mpi_enreg%paralbd

 end if ! has_file

 mpi_enreg%my_kpttab(:)=0
 mpi_enreg%my_isppoltab(:)=0
 do iisppol=1,nsppol
   ikpt_this_proc=0
   do iikpt=1,nkpt
     nband_k=nband(iikpt+(iisppol-1)*nkpt)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,iikpt,1,nband_k,iisppol,mpi_enreg%me_kpt)) cycle
     ikpt_this_proc=ikpt_this_proc+1
!    This test should be done when dataset are read and slipt of work do between processor
!    If this test is not good for one proc then other procs fall in deadlock->so PERS and MPI_ABORT
!    if (ikpt_this_proc > mkmem) then
!    message = ' this bandfft tab cannot be allocated !'
!    MSG_BUG(message)
!    end if
     mpi_enreg%my_kpttab(iikpt)=ikpt_this_proc
     mpi_enreg%my_isppoltab(iisppol)=1
   end do
 end do

end subroutine distrb2
!!***

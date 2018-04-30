!{\src2tex{textfont=tt}}
!!****f* ABINIT/distrb2_hf
!! NAME
!!  distrb2_hf
!!
!! FUNCTION
!!  This routine creates the tabs of repartition of processors
!!  for sharing the jobs on occupied states (labeled by k-points, 
!!  bands and spin indices) for an Hartree-Fock calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2018 ABINIT group (CMartins)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nbandhf = maximum number of occupied bands
!!  nkpthf = number of k-points in full BZ
!!  nproc= number of processors available for this distribution
!!  nsppol = 1 for unpolarized, 2 for polarized
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!   mpi_enreg%proc_distrb(nkpthf,nbandhf,nsppol)=number of the processor
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
!!      mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine distrb2_hf(nbandhf,nkpthf, nproc, nsppol, mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'distrb2_hf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nbandhf,nkpthf,nproc,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: ind,iiband,iikpt,iistep,nproc_hf
 character(len=500) :: message

!******************************************************************

 nproc_hf=mpi_enreg%nproc_hf

!* Initialize distrb_hf (the array exists necessarily)
 do iiband=1,nbandhf
   do iikpt=1,nkpthf
     mpi_enreg%distrb_hf(iikpt,iiband,1)=nproc_hf-1
   end do
 end do

!* End of the routine for an empty communication space
 if (nproc==0) return

!*** Testing section ***

 if (nsppol==2) then
!* Check that the distribution over (spin,k point) allow to consider spin up and spin dn independently.
   if (mpi_enreg%nproc_kpt/=1.and.mod(mpi_enreg%nproc_kpt,2)/=0) then
     MSG_ERROR( 'The variable nproc_kpt is not even but nssppol= 2')
!* In this case, one processor will carry both spin. (this will create pbms for the calculation)
   end if
!* Check that the number of band is the same for each spin, at each k-point. (It should be)
!*   do iikpt=1,nkpthf
!*     if (nband(iikpt)/=nband(iikpt+nkpthf)) then
!*     message = ' WARNING - the number of bands is different for spin up or spin down. '
!*     MSG_ERROR(message)
!*     end if
!*    end do
!* If one of this test is not good for one proc then other procs fall in deadlock, according to distrb2.
!* What does it mean ???
 end if


!* Check if nkpthf and nproc_hf match
 if (nproc_hf>nkpthf*nbandhf) then
!* There are too many processors with respect to nkpthf*nbandhf
   write(message, '(a,a,i4,a,i4,a,i4,a,a)' ) ch10,&
&   'nproc_hf=',nproc_hf,' >= nkpthf=',nkpthf,'* nbandhf=',nbandhf,ch10,&
&   'The number of processors is larger than nkpthf*nbandhf. This is a waste.'
   MSG_WARNING(message)

 else if(mod(nkpthf*nbandhf,nproc_hf)/=0) then
!* nkpthf*nbandhf is not a multiple of nproc_hf
   write(message, '(2a,i5,a,i5,3a)' ) ch10,&
&   'nkpthf*nbandhf (', nkpthf*nbandhf, ') is not a multiple of nproc_hf (',nproc_hf, ')', ch10,&
&   'The parallelisation may not be efficient.'
   MSG_WARNING(message)
 end if

!*** End of testing section ***

!*** Initialize the processor distribution from a simple algorithm ***

 if (nproc_hf<nkpthf) then
!* In this case, a parallelization over kpts only.
   iistep=nkpthf/nproc_hf
   if (mod(nkpthf,nproc_hf) /=0) iistep=iistep+1
   ind=0
   do iikpt=1,nkpthf
!*** Only the first "nbandhf" bands are considered (they are assumed to be the only occupied ones)
     do iiband=1,nbandhf
       mpi_enreg%distrb_hf(iikpt,iiband,1)=ind/iistep
     end do
     ind=ind+1
   end do

 else 
!* In this case, a parallelization over all the occupied states is possible.
   if (nproc_hf < nbandhf*nkpthf) then
     iistep=(nbandhf*nkpthf)/nproc_hf;
     if (mod((nbandhf*nkpthf),nproc_hf) /=0) iistep=iistep+1
   else 
     iistep=1
   end if
   ind=0
   do iikpt=1,nkpthf
!*** Only the first "nbandhf" bands are considered (they are assumed to be the only occupied ones)
     do iiband=1,nbandhf
       mpi_enreg%distrb_hf(iikpt,iiband,1)=ind/iistep
       ind=ind+1
     end do
   end do
 end if

!*** Initialization of processor distribution from a file (simple copy from distrb2, not yet implemented) ***

! !* Inquire whether there exist a file containing the processor distribution
!  if (first) then
! !  Case first time : test file to do
! !  Open the file containing the (k-points,bands) distribution
!    open(unit=temp_unit,file='kpt_distrb_hf',form='formatted',status='old',iostat=ios)
!    if(ios==0) then
! !    'kpt_distrb_hf' file exists
!      file_exist=1
!      close(temp_unit)
!    else
!      file_exist=0
!    end if
!    first=.false.
!  end if
! 
! !* Initialize the processor distribution, either from a file, or from an algorithm
!  if (file_exist == 1) then
! !* Read (k-points,bands) distribution out of the file
!    open(unit=temp_unit,file='kpt_distrb_hf',form='formatted',status='old',iostat=ios)
!    rewind(unit=temp_unit)
!    read(temp_unit,*) mpi_enreg%distrb_hf
!    close(temp_unit)
! !* Determine the range of processors requested
!    proc_max=0
!    proc_min=nproc_hf
!    do iikpt=1,nkpthf
!      mband_occ_k = mband_occ(iikpt+(iisppol-1)*nkpthf)
!      proc_max=maxval(mpi_enreg%distrb_hf(iikpt,1:mband_occ_k,1))
!      proc_min=minval(mpi_enreg%distrb_hf(iikpt,1:mband_occ_k,1))
!    end do
! 
!    if(proc_max>(nproc_hf-1)) then
! !*    Too much proc. requested
!      write(message, '(a,a,a,i4,a,a,a)' )&
! &     '  The number of processors mentioned in the kpt_distrb file',ch10,&
! &     '  must be lower or equal to the actual number of processors =',&
! &     nproc_hf-1,ch10,&
! &     '  Action : change the kpt_distrb file, or increase the',&
! &     '  number of processors.'
!      MSG_ERROR(message)
!    end if
! 
!    if(proc_max/=(nproc_hf-1)) then
! !*    Too few proc. used
!      write(message, '(a,i4,a,a,a,i4,a,a,a)' )&
! &     '  Only ',proc_max+1,' processors are used (from kpt_distrb file),',ch10,&
! &     '  when',nproc_hf,' processors are available.',ch10,&
! &     '  Action : adjust number of processors and kpt_distrb file.'
!      MSG_ERROR(message)
!    end if
! 
!    if(proc_min<0) then
!      write(message, '(a,a,a)' )&
! &     '  The number of processors must be bigger than 0 in kpt_distrb file.',ch10,&
! &     ' Action : modify kpt_distrb file.'
!      MSG_ERROR(message)
!    end if
!  else
! !* The file does not exist...
!  end if ! file_exist

!DEBUG
!write(std_out,*)' distrb2_hf: exit '
!ENDDEBUG

end subroutine distrb2_hf
!!***

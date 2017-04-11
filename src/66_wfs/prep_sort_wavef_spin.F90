!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_sort_wavef_spin
!!
!! NAME
!! prep_sort_wavef_spin
!!
!! FUNCTION
!! Compute index used to sort a spinorial wave-function by spin
!! Sort to have all nspinor=1 fisrt, then all nspinor=2
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MD)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nproc_band=size of "band" communicator
!!  nspinor=number of spinorial components of the wavefunction
!!  ndatarecv=total number of values on all processors
!!  recvcounts(nproc_band)= number of received values by the processor
!!  rdispls(nproc_band)= offsets of the received values by the processor
!!
!! OUTPUT
!!  index_wavef(:)=array containing the sorted indexes (pointer, allocated in this routine)
!!
!! PARENTS
!!      prep_getghc,prep_nonlop
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef)


 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_sort_wavef_spin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,allocatable,intent(out) :: index_wavef(:)

!Local variables-------------------------------
!scalars
 integer :: isft,isft1,iproc,iindex
!arrays
 integer,allocatable :: recvcountsloc(:),rdisplsloc(:)

! *********************************************************************

 ABI_ALLOCATE(index_wavef,(ndatarecv*nspinor))

 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc,(nproc_band))
 recvcountsloc(:)=recvcounts(:)*2*nspinor
 rdisplsloc(:)=rdispls(:)*2*nspinor

!---------------------------------------------
!Loops on bandpp and processors band
!---------------------------------------------
 isft=0
 do iproc=1,nproc_band

!  ===== Spin up
   if (iproc==1) then
     isft= 0
   else
     isft= sum(recvcounts(1: (iproc-1)))
   end if
   isft1 = 0.5*rdisplsloc(iproc)

   index_wavef(1+isft:isft+recvcounts(iproc))= &
&   (/(iindex,iindex=isft1+1,isft1+recvcounts(iproc))/)

!  =====Spin down
   if (iproc==1)then
     isft=sum(recvcounts(1:nproc_band))
   else
     isft=sum(recvcounts(1:nproc_band)) &
&     +sum(recvcounts(1:iproc-1))
   end if
   isft1 = 0.5 * rdisplsloc(iproc) + recvcounts(iproc)

   index_wavef(1+isft:isft+recvcounts(iproc))= &
&   (/(iindex,iindex=isft1+1,isft1+ recvcounts(iproc))/)

 end do

 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)

end subroutine prep_sort_wavef_spin
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/ep_setupqpt
!!
!! NAME
!! ep_setupqpt
!!
!! FUNCTION
!!  set up qpoint grid for elphon.
!!  2 modes, either uniform grid from anaddb input nqpt
!!  or take qpt from anaddb input (explicitly listed)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   crystal>crystal_t>=data type gathering info on the crystalline structure.
!!   anaddb_dtset=dataset with input variables
!!     %qgrid_type gives type of q grid 1=uniform 2=take from input
!!     %ep_nqpt    number of auxiliary qpoints 
!!     %ep_qptlist list of qpoints, 
!!
!! OUTPUT
!!   
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      getkgrid,smpbz,symkpt,wrap2_pmhalf,wrtout
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ep_setupqpt (elph_ds,crystal,anaddb_dtset,qptrlatt,timrev)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,   only : wrap2_pmhalf
 use m_anaddb_dataset,  only : anaddb_dataset_type
 use m_crystal,         only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ep_setupqpt'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer, intent(in) :: timrev
 type(crystal_t),intent(in) :: crystal
 type(anaddb_dataset_type), intent(in) :: anaddb_dtset
 type(elph_type), intent(inout) :: elph_ds
!arrays
 integer, intent(out) :: qptrlatt(3,3) 

!Local variables -------------------------
!scalars
 integer :: nqshft,option,iqpt, nqpt1
 integer :: iscf,mqpt,iout,berryopt,nqpt_computed
 real(dp) :: qptrlen, res
 character(len=500) :: message
!arrays
 integer :: vacuum(3)
 integer,allocatable :: indqpt1(:)
 real(dp) :: kpt(3)
 real(dp),allocatable :: wtq_folded(:)
 real(dp), allocatable :: wtq(:),qpt_full(:,:),tmpshifts(:,:)

! *********************************************************************

!default is to expect a uniform grid
 elph_ds%tuniformgrid = 1

!if we use the normal grid way of generating the qpoints:
 if (anaddb_dtset%qgrid_type==1) then
!  qpoint lattice vectors (inverse, like kptrlatt)
   qptrlatt(:,:)=0
   qptrlatt(1,1)=anaddb_dtset%ngqpt(1)
   qptrlatt(2,2)=anaddb_dtset%ngqpt(2)
   qptrlatt(3,3)=anaddb_dtset%ngqpt(3)
   
   if (anaddb_dtset%nqshft /= 1) then
!    try to reduce the qpoint grid to a single qshift, otherwise stop
!    dummy args for call to getkgrid
     vacuum(:) = 0
     iscf = 3
     
     mqpt = anaddb_dtset%ngqpt(1)*anaddb_dtset%ngqpt(2)*anaddb_dtset%ngqpt(3)*anaddb_dtset%nqshft
     ABI_ALLOCATE(qpt_full,(3,mqpt))
     ABI_ALLOCATE(wtq,(mqpt))
     ABI_ALLOCATE(tmpshifts,(3,210))

     wtq(:) = one

     tmpshifts(:,:) = zero
     tmpshifts(:,1:4) = anaddb_dtset%q1shft(:,:)

     iout=6

     berryopt = 1

!    just call with identity, to get full set of kpts in qpt_full, but
!    reduce qshfts
     
     nqshft=anaddb_dtset%nqshft
     call getkgrid(0,0,iscf,qpt_full,3,qptrlatt,qptrlen, &
&     1,mqpt,nqpt_computed,nqshft,1,crystal%rprimd,tmpshifts,crystal%symafm, &
&     crystal%symrel,vacuum,wtq)
     ABI_DEALLOCATE(qpt_full)
     ABI_DEALLOCATE(wtq)
     ABI_DEALLOCATE(tmpshifts)

     if (anaddb_dtset%nqshft /= 1) then
       write (message,'(a,i0)')&
&       ' multiple qpt shifts not treated yet (should be possible), nqshft= ', anaddb_dtset%nqshft
       MSG_ERROR(message)
     end if
   end if  ! end multiple shifted qgrid


   write(message,'(a,9(i0,1x))')' elphon : enter smpbz with  qptrlatt = ',qptrlatt 
   call wrtout(std_out,message,'COLL')

   option=1
!  mqpt=anaddb_dtset%ngqpt(1)*anaddb_dtset%ngqpt(2)*anaddb_dtset%ngqpt(3)*anaddb_dtset%nqshft
   mqpt= qptrlatt(1,1)*qptrlatt(2,2)*qptrlatt(3,3) &
&   +qptrlatt(1,2)*qptrlatt(2,3)*qptrlatt(3,1) &
&   +qptrlatt(1,3)*qptrlatt(2,1)*qptrlatt(3,2) &
&   -qptrlatt(1,2)*qptrlatt(2,1)*qptrlatt(3,3) &
&   -qptrlatt(1,3)*qptrlatt(2,2)*qptrlatt(3,1) &
&   -qptrlatt(1,1)*qptrlatt(2,3)*qptrlatt(3,2)

   ABI_ALLOCATE(qpt_full,(3,mqpt))
   iout = 6
   call smpbz(anaddb_dtset%brav,iout,qptrlatt,mqpt,elph_ds%nqpt_full,anaddb_dtset%nqshft,option,anaddb_dtset%q1shft,qpt_full)


!  save the q-grid for future reference
   ABI_ALLOCATE(elph_ds%qpt_full,(3,elph_ds%nqpt_full))

!  reduce qpt_full to correct zone
   do iqpt=1,elph_ds%nqpt_full
     call wrap2_pmhalf(qpt_full(1,iqpt),kpt(1),res)
     call wrap2_pmhalf(qpt_full(2,iqpt),kpt(2),res)
     call wrap2_pmhalf(qpt_full(3,iqpt),kpt(3),res)
     qpt_full(:,iqpt) = kpt
     elph_ds%qpt_full(:,iqpt)=kpt
   end do
   ABI_DEALLOCATE(qpt_full)

 else if (anaddb_dtset%qgrid_type==2) then ! use explicit list of qpoints from anaddb input
   qptrlatt(:,:)=0
   qptrlatt(1,1)=1
   qptrlatt(2,2)=1
   qptrlatt(3,3)=1

   elph_ds%nqpt_full=anaddb_dtset%ep_nqpt
   ABI_ALLOCATE(elph_ds%qpt_full,(3,elph_ds%nqpt_full))

   elph_ds%qpt_full = anaddb_dtset%ep_qptlist

   elph_ds%tuniformgrid = 0
 end if ! type of qgrid for elphon

!=================================================================
!Calculate weights, needed to estimate lambda using the weighted
!sum of the uninterpolated e-ph matrix elements
!=================================================================
 call wrtout(std_out,' setqgrid : calling symkpt to find irred q points',"COLL")

 ABI_ALLOCATE(indqpt1,(elph_ds%nqpt_full))
 ABI_ALLOCATE(wtq_folded,(elph_ds%nqpt_full))
 ABI_ALLOCATE(wtq,(elph_ds%nqpt_full))

 wtq(:) = one/dble(elph_ds%nqpt_full) !weights normalized to unity

!
!NOTE: this reduction of irred qpt may not be identical to that in GKK file
!which would be more practical to use.
!
 iout=0 !do not write to ab_out
!should we save indqpt1 for use inside elph_ds?
 call symkpt(0,crystal%gmet,indqpt1,iout,elph_ds%qpt_full,elph_ds%nqpt_full,nqpt1,crystal%nsym,crystal%symrec,&
& timrev,wtq,wtq_folded)

 write (message,'(2a,i0)')ch10,' Number of irreducible q-points = ',nqpt1
 call wrtout(std_out,message,'COLL')
 elph_ds%nqptirred=nqpt1

 call wrtout(std_out,' === Irreducible q points with weights ==== ','COLL')

 do iqpt=1,elph_ds%nqpt_full
   if (wtq_folded(iqpt) /= zero) then
     write (message,'(1x,i4,a2,4es16.8)')iqpt,') ',elph_ds%qpt_full(:,iqpt),wtq_folded(iqpt)
     call wrtout(std_out,message,'COLL')
   end if
 end do

 call wrtout(std_out,ch10,'COLL')

 ABI_ALLOCATE(elph_ds%wtq,(elph_ds%nqpt_full))

 elph_ds%wtq(:)=wtq_folded(:)
!MEMO indqpt could be useful to test the qgrid read by abinit
 ABI_DEALLOCATE(indqpt1)
 ABI_DEALLOCATE(wtq_folded)
 ABI_DEALLOCATE(wtq)

end subroutine ep_setupqpt
!!***

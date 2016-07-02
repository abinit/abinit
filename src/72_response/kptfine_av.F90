!!****f* ABINIT/kptfine_av
!! NAME
!! kptfine_av
!!
!! FUNCTION
!! This routine compute the k-points of a fine grid that are around a 
!! define k-point of a coarse mesh.  
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (SP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  center(3) = the point of the coarse mesh around which you want know which
!!              k-points of the fine mesh belong to.
!!  qptrlatt(3,3) = qptrlatt of the considered calculation (this is obtained
!!              from the input variable ngqpt and shiftq.
!!  kpt_fine(3,nkpt_fine) = this table contain all the k-points of the fine grid
!!              in the full BZ (no sym op. allowed) and is read from the header
!!              of the dense WF file.
!!  nkpt_fine = number of k-points of the fine grid read from the header of the
!!              dense WF file. 
!!
!! OUTPUT
!!  kpt_fine_sub(nkpt_sub) = k-points of the fine grid that are around center(3)
!!  nkpt_sub = number of k-points of the fine grid that are around center(3)
!!  wgt_sub(nkpt_sub) = weight of the k-points of the fine grid that are around
!!              center(3).  
!!
!! PARENTS
!!      eig2stern,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine kptfine_av(center,qptrlatt,kpt_fine,nkpt_fine,kpt_fine_sub,&
&                nkpt_sub,wgt_sub)

 use defs_basis
 use m_profiling_abi
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kptfine_av'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: nkpt_fine
 integer,intent(out)  :: nkpt_sub
!arrays
 integer,intent(in)   :: qptrlatt(3,3)
 real(dp),intent(in)  :: kpt_fine(3,nkpt_fine)
 real(dp),intent(in)  :: center(3)
 integer,pointer      :: kpt_fine_sub(:)
 real(dp),pointer     :: wgt_sub(:)

!Local variables-------------------------------
!scalars
 integer :: ikpt,aa,bb,cc
 integer :: ii,jj
!arrays
 real(dp) :: center_ref(3)
 real(dp) :: kpt_fine_ref(3)
 real(dp) :: kpt_tmp(3),kpt_tmp2(3)
 integer,allocatable  :: kpt_fine_sub_tmp(:)
 real(dp),allocatable :: wgt_sub_tmp(:)
 logical :: found(3)

! *************************************************************************

 ABI_ALLOCATE(kpt_fine_sub_tmp,(nkpt_fine))
 ABI_ALLOCATE(wgt_sub_tmp,(nkpt_fine))

!It is easier to work in real space using the qptrlatt matrices because in this
!referential any k-points sampling will be cast into an orthorhombic shape.
!In that space we can simply take all k-points of the fine grid that between
!center_ref-0.5 and center_ref+0.5

 center_ref = MATMUL(qptrlatt,center)

!When considering points center(3) that lying close or on a BZ edge we need to
!take the k-points of the fine grid taking into account unklamp vectors. This
!is done with the aa, bb and cc loops. 

 ii = 1
 do ikpt=1,nkpt_fine
   kpt_tmp = kpt_fine(:,ikpt)
   do aa=-1,1
     kpt_tmp2(1) = kpt_tmp(1)+aa
     do bb=-1,1
       kpt_tmp2(2) = kpt_tmp(2)+bb
       do cc=-1,1
         kpt_tmp2(3) = kpt_tmp(3)+cc
         kpt_fine_ref = MATMUL(qptrlatt,kpt_tmp2)
         if((kpt_fine_ref(1)>=center_ref(1)-0.5-tol8).and.&
&         (kpt_fine_ref(1)<=center_ref(1)+0.5+tol8)) then
           if((kpt_fine_ref(2)>=center_ref(2)-0.5-tol8).and.&
&           (kpt_fine_ref(2)<=center_ref(2)+0.5+tol8)) then
             if((kpt_fine_ref(3)>=center_ref(3)-0.5-tol8).and.&
&             (kpt_fine_ref(3)<=center_ref(3)+0.5+tol8)) then
               kpt_fine_sub_tmp(ii) = ikpt
               ii = ii +1
             end if
           end if
         end if
       end do
     end do
   end do
 end do

 nkpt_sub = ii-1
 ABI_ALLOCATE(kpt_fine_sub,(nkpt_sub))
 ABI_ALLOCATE(wgt_sub,(nkpt_sub))

 do jj=1,nkpt_sub
   kpt_fine_sub(jj) = kpt_fine_sub_tmp(jj)
 end do

!We then compute a weight function. This weight function is simply a
!rectangular weight function that take the value 1 for k-points of the fine
!grid inside the cube, 0.5 for k-points that are lying on one face of the cube,
!0.25 for k-points that are lying on an edge of the cube and 0.125 for k-points
!that are lying on a peak of the cube.  

 wgt_sub(:) = 1.0

 do ikpt=1,nkpt_sub
   found(:) = .True.
   kpt_tmp = kpt_fine(:,kpt_fine_sub(ikpt))
   do aa=-1,1
     kpt_tmp2(1) = kpt_tmp(1)+aa
     do bb=-1,1
       kpt_tmp2(2) = kpt_tmp(2)+bb
       do cc=-1,1
         kpt_tmp2(3) = kpt_tmp(3)+cc
         kpt_fine_ref = MATMUL(qptrlatt,kpt_tmp2)
         if((ABS(kpt_fine_ref(1)-center_ref(1)-0.5)< tol8) .or.&
&         (ABS(kpt_fine_ref(1)-center_ref(1)+0.5) < tol8)) then
           if(found(1)) then
             wgt_sub(ikpt) = wgt_sub(ikpt)*0.5
             found(1) = .False.
           end if
         end if
         if((ABS(kpt_fine_ref(2)-center_ref(2)-0.5) < tol8) .or.&
&         (ABS(kpt_fine_ref(2)-center_ref(2)+0.5) < tol8)) then
           if(found(2)) then
             wgt_sub(ikpt) = wgt_sub(ikpt)*0.5
             found(2) = .False.
           end if
         end if
         if((ABS(kpt_fine_ref(3)-center_ref(3)-0.5)< tol8) .or.&
&         (ABS(kpt_fine_ref(3)-center_ref(3)+0.5) < tol8)) then
           if(found(3)) then
             wgt_sub(ikpt) = wgt_sub(ikpt)*0.5
             found(3) = .False.
           end if
         end if
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(kpt_fine_sub_tmp)
 ABI_DEALLOCATE(wgt_sub_tmp)

end subroutine kptfine_av
!!***

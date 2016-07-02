!{\src2tex{textfont=tt}}
!!****f* ABINIT/invacuum
!!
!! NAME
!! invacuum
!!
!! FUNCTION
!! Determine whether there is vacuum along some of the primitive directions
!! in real space.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! natom=number of atoms
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! string*(*)=character string containing all the input data.
!!  Initialized previously in instrng.
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!! vacuum(3)= for each direction, 0 if no vacuum, 1 if vacuum
!!
!! PARENTS
!!      invars1,invars2
!!
!! CHILDREN
!!      intagm,metric,sort_dp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine invacuum(jdtset,lenstr,natom,rprimd,string,vacuum,xred)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invacuum'
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: jdtset,lenstr,natom
 character(len=*),intent(in) :: string
!arrays
 integer,intent(out) :: vacuum(3)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ii,marr,tread
 real(dp) :: max_diff_xred,ucvol,vacwidth,vacxred
!arrays
 integer,allocatable :: list(:)
 integer,allocatable :: intarr(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: xred_sorted(:)
 real(dp),allocatable :: dprarr(:)

! *************************************************************************

!Compute the maximum size of arrays intarr and dprarr
 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!Get metric quantities
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Read vacwidth, or set the default
 vacwidth=10.0_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'vacwidth',tread,'LEN')
 if(tread==1) vacwidth=dprarr(1)

!Read vacuum, or compute it using the atomic coordinates and vacwidth.
 vacuum(1:3)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'vacuum',tread,'INT')

 if(tread==1)then
   vacuum(1:3)=intarr(1:3)
 else
!  For each direction, determine whether a vacuum space exists
   ABI_ALLOCATE(list,(natom))
   ABI_ALLOCATE(xred_sorted,(natom))
   do ii=1,3
!    This is the minimum xred difference needed to have vacwidth
     vacxred=vacwidth*sqrt(sum(gprimd(:,ii)**2))
!    Project the reduced coordinate in the [0.0_dp,1.0_dp[ interval
     xred_sorted(:)=mod(xred(ii,:),1.0_dp)
!    list is dummy
     list(:)=0
!    Sort xred_sorted
     call sort_dp(natom,xred_sorted,list,tol14)
     if(natom==1)then
       max_diff_xred=1.0_dp
     else
!      Compute the difference between each pair of atom in the sorted order
       max_diff_xred=0.0_dp
       do ia=1,natom-1
         max_diff_xred=max(max_diff_xred,xred_sorted(ia+1)-xred_sorted(ia))
       end do
!      Do not forget the image of the first atom in the next cell
       max_diff_xred=max(max_diff_xred,1.0_dp+xred_sorted(1)-xred_sorted(ia))
     end if
     if(vacxred<max_diff_xred+tol10)vacuum(ii)=1
   end do
   ABI_DEALLOCATE(list)
   ABI_DEALLOCATE(xred_sorted)
 end if

!DEBUG
!write(std_out,*)' invacuum : vacuum=',vacuum(1:3)
!ENDDEBUG

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine invacuum
!!***

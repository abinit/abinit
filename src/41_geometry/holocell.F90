!{\src2tex{textfont=tt}}
!!****f* ABINIT/holocell
!! NAME
!! holocell
!!
!! FUNCTION
!! Examine whether the trial conventional cell described by cell_base
!! is coherent with the required holohedral group.
!! Possibly enforce the holohedry and modify the basis vectors.
!! Note : for iholohedry=4, the tetragonal axis is not required to be
!! along the C axis.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  enforce= if 0, only check; if =1, enforce exactly the holohedry
!!  iholohedry=required holohegral group
!!  iholohedry=1   triclinic      1bar
!!  iholohedry=2   monoclinic     2/m
!!  iholohedry=3   orthorhombic   mmm
!!  iholohedry=4   tetragonal     4/mmm
!!  iholohedry=5   trigonal       3bar m
!!  iholohedry=6   hexagonal      6/mmm
!!  iholohedry=7   cubic          m3bar m
!!  tolsym=tolerance for the symmetry operations
!!
!! OUTPUT
!!  foundc=1 if the basis vectors supports the required holohedry ; =0 otherwise
!!
!! SIDE EFFECTS
!!  cell_base(3,3)=basis vectors of the conventional cell  (changed if enforce==1, otherwise unchanged)
!!
!! PARENTS
!!      symlatt,symmetrize_rprimd
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine holocell(cell_base,enforce,foundc,iholohedry,tolsym)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'holocell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enforce,iholohedry
 integer,intent(out) :: foundc
 real(dp),intent(in) :: tolsym
!arrays
 real(dp),intent(inout) :: cell_base(3,3)

!Local variables ------------------------------
!scalars
 integer :: allequal,ii,jj,orth
 real(dp):: aa,reldiff,scprod1
 character(len=500) :: msg
!arrays
 integer :: ang90(3),equal(3)
 real(dp) :: length(3),metric(3,3),norm(3),rbasis(3,3),rconv(3,3),rconv_new(3,3)
 real(dp) :: rnormalized(3,3),symmetrized_length(3)

!**************************************************************************

 do ii=1,3
   metric(:,ii)=cell_base(1,:)*cell_base(1,ii)+&
&   cell_base(2,:)*cell_base(2,ii)+&
&   cell_base(3,:)*cell_base(3,ii)
 end do

!Examine the angles and vector lengths
 ang90(:)=0
 if(metric(1,2)**2<tolsym**2*metric(1,1)*metric(2,2))ang90(3)=1
 if(metric(1,3)**2<tolsym**2*metric(1,1)*metric(3,3))ang90(2)=1
 if(metric(2,3)**2<tolsym**2*metric(2,2)*metric(3,3))ang90(1)=1
 orth=0
 if(ang90(1)==1 .and. ang90(2)==1 .and. ang90(3)==1) orth=1
 equal(:)=0
 if(abs(metric(1,1)-metric(2,2))<tolsym*half*(metric(1,1)+metric(2,2)))equal(3)=1
 if(abs(metric(1,1)-metric(3,3))<tolsym*half*(metric(1,1)+metric(3,3)))equal(2)=1
 if(abs(metric(2,2)-metric(3,3))<tolsym*half*(metric(2,2)+metric(3,3)))equal(1)=1
 allequal=0
 if(equal(1)==1 .and. equal(2)==1 .and. equal(3)==1) allequal=1

 foundc=0
 if(iholohedry==1)                                      foundc=1
 if(iholohedry==2 .and. ang90(1)+ang90(3)==2 )          foundc=1
 if(iholohedry==3 .and. orth==1)                        foundc=1
 if(iholohedry==4 .and. orth==1 .and.          &
& (equal(3)==1 .or. equal(2)==1 .or. equal(1)==1) ) foundc=1
 if(iholohedry==5 .and. allequal==1 .and. &
& (abs(metric(1,2)-metric(2,3))<tolsym*metric(2,2)) .and. &
& (abs(metric(1,2)-metric(1,3))<tolsym*metric(1,1))         )      foundc=1
 if(iholohedry==6 .and. equal(3)==1 .and. &
& ang90(1)==1 .and. ang90(2)==1 .and. &
& (2*metric(1,2)-metric(1,1))<tolsym*metric(1,1) )      foundc=1
 if(iholohedry==7 .and. orth==1 .and. allequal==1)      foundc=1

!-------------------------------------------------------------------------------------
!Possibly enforce the holohedry (if it is to be enforced !)

 if(foundc==1.and.enforce==1.and.iholohedry/=1)then

!  Copy the cell_base vectors, and possibly fix the tetragonal axis to be the c-axis
   if(iholohedry==4.and.equal(1)==1)then
     rconv(:,3)=cell_base(:,1) ; rconv(:,1)=cell_base(:,2) ; rconv(:,2)=cell_base(:,3)
   else if (iholohedry==4.and.equal(2)==1)then
     rconv(:,3)=cell_base(:,2) ; rconv(:,2)=cell_base(:,1) ; rconv(:,1)=cell_base(:,3)
   else
     rconv(:,:)=cell_base(:,:)
   end if

!  Compute the length of the three conventional vectors
   length(1)=sqrt(sum(rconv(:,1)**2))
   length(2)=sqrt(sum(rconv(:,2)**2))
   length(3)=sqrt(sum(rconv(:,3)**2))

!  Take care of the first conventional vector aligned with rbasis(:,3) (or aligned with the trigonal axis if rhombohedral)
!  and choice of the first normalized direction 
   if(iholohedry==5)then
     rbasis(:,3)=third*(rconv(:,1)+rconv(:,2)+rconv(:,3))
   else
     rbasis(:,3)=rconv(:,3)
   end if
   norm(3)=sqrt(sum(rbasis(:,3)**2))
   rnormalized(:,3)=rbasis(:,3)/norm(3)

!  Projection of the first conventional vector perpendicular to rbasis(:,3)
!  and choice of the first normalized direction 
   scprod1=sum(rnormalized(:,3)*cell_base(:,1))
   rbasis(:,1)=rconv(:,1)-rnormalized(:,3)*scprod1
   norm(1)=sqrt(sum(rbasis(:,1)**2))
   rnormalized(:,1)=rbasis(:,1)/norm(1)
   
!  Generation of the second vector, perpendicular to the third and first
   rnormalized(1,2)=rnormalized(2,3)*rnormalized(3,1)-rnormalized(3,3)*rnormalized(2,1)
   rnormalized(2,2)=rnormalized(3,3)*rnormalized(1,1)-rnormalized(1,3)*rnormalized(3,1)
   rnormalized(3,2)=rnormalized(1,3)*rnormalized(2,1)-rnormalized(2,3)*rnormalized(1,1)

!  Compute the vectors of the conventional cell, on the basis of iholohedry
   if(iholohedry==2)then
     rconv_new(:,3)=rconv(:,3)
     rconv_new(:,1)=rconv(:,1)
     rconv_new(:,2)=rnormalized(:,2)*length(2) ! Now, the y axis is perpendicular to the two others, that have not been changed
   else if(iholohedry==3.or.iholohedry==4.or.iholohedry==7)then
     if(iholohedry==7)then
       symmetrized_length(1:3)=sum(length(:))*third
     else if(iholohedry==4)then
       symmetrized_length(3)=length(3)
       symmetrized_length(1:2)=half*(length(1)+length(2))
     else if(iholohedry==3)then
       symmetrized_length(:)=length(:)
     end if
     do ii=1,3
       rconv_new(:,ii)=rnormalized(:,ii)*symmetrized_length(ii)
     end do
   else if(iholohedry==5)then
!    In the normalized basis, they have coordinates (a,0,c), and (-a/2,+-sqrt(3)/2*a,c)
!    c is known, but a is computed from the knowledge of the average length of the initial vectors
     aa=sqrt(sum(length(:)**2)*third-norm(3)**2)
     rconv_new(:,1)=aa*rnormalized(:,1)+rbasis(:,3)
     rconv_new(:,2)=aa*half*(-rnormalized(:,1)+sqrt(three)*rnormalized(:,2))+rbasis(:,3)
     rconv_new(:,3)=aa*half*(-rnormalized(:,1)-sqrt(three)*rnormalized(:,2))+rbasis(:,3)
   else if(iholohedry==6)then

!    In the normalized basis, they have coordinates (a,0,0), (-a/2,+-sqrt(3)/2*a,0), and (0,0,c)
!    c is known, but a is computed from the knowledge of the average length of the initial vectors
     aa=half*(length(1)+length(2))
     rconv_new(:,1)=aa*rnormalized(:,1)
     rconv_new(:,2)=aa*half*(-rnormalized(:,1)+sqrt(three)*rnormalized(:,2))
     rconv_new(:,3)=rconv(:,3)
   end if

!  Check whether the modification make sense 
   do ii=1,3
     do jj=1,3
       reldiff=(rconv_new(ii,jj)-rconv(ii,jj))/length(jj)
!      Allow for twice tolsym
       if(abs(reldiff)>two*tolsym)then
         write(msg,'(a,6(2a,3es14.6))')&
&         'Failed rectification of lattice vectors to comply with Bravais lattice identification, modifs are too large',ch10,&
&         '  rconv    =',rconv(:,1),ch10,&
&         '            ',rconv(:,2),ch10,&
&         '            ',rconv(:,3),ch10,&
&         '  rconv_new=',rconv_new(:,1),ch10,&
&         '            ',rconv_new(:,2),ch10,&
&         '            ',rconv_new(:,3)
         MSG_ERROR_CLASS(msg, "TolSymError")
       end if
     end do
   end do

!  Copy back the cell_base vectors
   if(iholohedry==4.and.equal(1)==1)then
     cell_base(:,3)=rconv_new(:,2) ; cell_base(:,2)=rconv_new(:,1) ; cell_base(:,1)=rconv_new(:,3)
   else if (iholohedry==4.and.equal(2)==1)then
     cell_base(:,3)=rconv_new(:,1) ; cell_base(:,1)=rconv_new(:,2) ; cell_base(:,2)=rconv_new(:,3)
   else
     cell_base(:,:)=rconv_new(:,:)
   end if

 end if

end subroutine holocell
!!***

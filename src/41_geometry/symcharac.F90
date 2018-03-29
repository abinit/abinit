!{\src2tex{textfont=tt}}
!!****f* ABINIT/symcharac
!! NAME
!! symcharac
!!
!! FUNCTION
!! Get the type of axis for the symmetry.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT group (RC, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! center=bravais(2)
!! determinant=the value of the determinant of sym
!! iholohedry=bravais(1)
!! isym=number of the symmetry operation that is currently analyzed
!! order=the order of the symmetry
!! symrel(3,3)= the symmetry matrix
!! tnons(3)=nonsymmorphic translations
!!
!! OUTPUT
!! label=a human readable text for the characteristic of the symmetry
!! type_axis=an identifier for the type of symmetry
!!
!! PARENTS
!!      m_ab7_symmetry,symspgr
!!
!! CHILDREN
!!      symaxes,symplanes,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symcharac(center, determinant, iholohedry, isym, label, symrel, tnons, type_axis)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symcharac'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => symcharac
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: determinant, center, iholohedry, isym
 integer, intent(out) :: type_axis
 character(len = 128), intent(out) :: label
 !arrays
 integer,intent(in) :: symrel(3,3)
 real(dp),intent(in) :: tnons(3)

 !Local variables-------------------------------
 !scalars
 logical,parameter :: verbose=.FALSE.
 integer :: tnons_order, identified, ii, order, iorder
 character(len=500) :: message
 !arrays
 integer :: identity(3,3),matrix(3,3),trial(3,3)
 real(dp) :: reduced(3),trialt(3)

 !**************************************************************************

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
 trial(:,:)=identity(:,:)
 matrix(:,:)=symrel(:,:)

 order=0
 do iorder=1,6
   trial=matmul(matrix,trial)
   if(sum((trial-identity)**2)==0)then
     order=iorder
     exit
   end if
   if(sum((trial+identity)**2)==0)then
     order=iorder
     exit
   end if
 end do

 if(order==0)then
   type_axis = -2
   return
 end if

!Determination of the characteristics of proper symmetries (rotations)
 if(determinant==1)then

!  Determine the translation vector associated to the rotations
!  and its order : apply the symmetry operation
!  then analyse the resulting vector.
   identified=0
   trialt(:)=zero
   do ii=1,order
     trialt(:)=matmul(symrel(:,:),trialt(:))+tnons(:)
   end do
!  Gives the associated translation, with components in the
!  interval [-0.5,0.5] .
   reduced(:)=trialt(:)-nint(trialt(:)-tol6)

   if(sum(abs(reduced(:)))<tol6)identified=1
   if( (center==1 .or. center==-3) .and. &
&   sum(abs(reduced(:)-(/zero,half,half/)))<tol6 )identified=2
   if( (center==2 .or. center==-3) .and. &
&   sum(abs(reduced(:)-(/half,zero,half/)))<tol6 )identified=3
   if( (center==3 .or. center==-3) .and. &
&   sum(abs(reduced(:)-(/half,half,zero/)))<tol6 )identified=4
   if(center==-1.and. sum(abs(reduced(:)-(/half,half,half/)))<tol6 )identified=5

!  If the symmetry operation has not been identified, there is a problem ...
   if(identified==0) then
     type_axis = -1
     return
   end if

!  Compute the translation vector associated with one rotation
   trialt(:)=trialt(:)/order
   trialt(:)=trialt(:)-nint(trialt(:)-tol6)

!  Analyse the resulting vector.
   identified=0
   do ii=1,order
     reduced(:)=ii*trialt(:)-nint(ii*trialt(:)-tol6)
     if(sum(abs(reduced(:)))<tol6)identified=1
     if( (center==1 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/zero,half,half/)))<tol6 )identified=2
     if( (center==2 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/half,zero,half/)))<tol6 )identified=3
     if( (center==3 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/half,half,zero/)))<tol6 )identified=4
     if(center==-1.and. sum(abs(reduced(:)-(/half,half,half/)))<tol6 )identified=5

     if(identified/=0)then
       tnons_order=ii
       exit
     end if
   end do ! ii

!  Determinant (here=+1, as we are dealing with proper symmetry operations),
!  order, tnons_order and identified are enough to
!  determine the kind of symmetry operation

   select case(order)
   case(1)                       ! point symmetry 1
     if(identified==1) then
       type_axis=8                 ! 1
       write(label,'(a)') 'the identity'
     else
       type_axis=7                 ! t
       write(label,'(a)') 'a pure translation '
     end if

     if (verbose) then
       write(message,'(a,i3,2a)')' symspgr : the symmetry operation no. ',isym,' is ',trim(label)
       call wrtout(std_out,message,'COLL')
     end if

   case(2,3,4,6)                 ! point symmetry 2,3,4,6 - rotations
     call symaxes(center,iholohedry,isym,symrel,label,order,tnons_order,trialt,type_axis)
   end select

 else if (determinant==-1)then

!  Now, take care of the improper symmetry operations.
!  Their treatment is relatively easy, except for the mirror planes
   select case(order)
   case(1)                       ! point symmetry 1
     type_axis=5                  ! -1
     write(label,'(a)') 'an inversion'
   case(2)                       ! point symmetry 2 - planes
     call symplanes(center,iholohedry,isym,symrel,tnons,label,type_axis)
   case(3)                       ! point symmetry 3
     type_axis=3                  ! -3
     write(label,'(a)') 'a -3 axis '
   case(4)                       ! point symmetry 1
     type_axis=2                  ! -4
     write(label,'(a)') 'a -4 axis '
   case(6)                       ! point symmetry 1
     type_axis=1                  ! -6
     write(label,'(a)') 'a -6 axis '
   end select

   if (order /= 2 .and. verbose) then
     write(message,'(a,i3,2a)')' symspgr : the symmetry operation no. ',isym,' is ',trim(label)
     call wrtout(std_out,message,'COLL')
   end if

 end if ! determinant==1 or -1

end subroutine symcharac
!!***

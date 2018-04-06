!{\src2tex{textfont=tt}}
!!****f* ABINIT/symplanes
!! NAME
!! symplanes
!!
!! FUNCTION
!! Determines the type of symmetry mirror planes: m,a,b,c,d,n,g.
!! This is used (see symlist.f) to identify the space group.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT group (RC, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! center=type of bravais lattice centering
!!   center=0        no centering
!!   center=-1       body-centered
!!   center=-3       face-centered
!!   center=1        A-face centered
!!   center=2        B-face centered
!!   center=3        C-face centered
!! iholohedry=type of holohedry
!!   iholohedry=1   triclinic      1bar
!!   iholohedry=2   monoclinic     2/m
!!   iholohedry=3   orthorhombic   mmm
!!   iholohedry=4   tetragonal     4/mmm
!!   iholohedry=5   trigonal       3bar m
!!   iholohedry=6   hexagonal      6/mmm
!!   iholohedry=7   cubic          m3bar m
!! isym=number of the symmetry operation that is currently analyzed
!! isymrelconv=symrel matrix for the particular operation, in conv. coord.
!! itnonsconv=tnons vector for the particular operation, in conv. coord
!!
!! OUTPUT
!! label=user friendly label of the plane
!! type_axis=type of the symmetry operation
!!
!! NOTES
!! One follows the
!! conventions explained in table 1.3 of the international tables for
!! crystallography. In the case of the rhombohedral system,
!! one takes into account the first footnote of this table 1.3 .
!! In general, we will assign the different symmetries to
!! the following numbers :  m -> 15 , (a, b or c) -> 16,
!!  d -> 17, n -> 18 , g -> 19
!! However, there is the same problem as for binary axes,
!! namely, for parallel mirror planes, one can find different
!! translation vectors, and these might be found at random,
!! depending on the input tnons.
!! (1) In the tP case, one will distinguish tertiary
!!  mirror plane, for which it is important to know whether they are
!!  m or c (for tertiary planes in tP, g is equivalent to m and n is equivalent to c).
!!  On the other hand, it is important to distinguish among
!!  primary and secondary mirror planes, those that are m,(a or b),c, or n.
!!  To summarize, the number of the symmetry will be :
!!  m (primary, secondary or tertiary) -> 15 ,
!!  secondary (a or b) -> 16, secondary c -> 17,
!!  primary or secondary n -> 18 , tertiary c -> 19
!! (2) In the tI case, one will distinguish tertiary
!!  mirror plane, for which it is important to know whether they are
!!  m or d (for tertiary planes in tI, c is equivalent to m.
!!  On the other hand, it is important to distinguish among
!!  primary and secondary mirror planes, those that are m (equivalent to n),
!!  or a,b or c.
!!  To summarize, the number of the symmetry will be :
!!  m (primary, secondary, tertiary) -> 15 ,
!!  a,b or c (primary or secondary) -> 16, tertiary d -> 17
!! (3) For hP and hR, a m plane is always coupled to a a or b plane,
!!  while a c plane is always coupled to an n plane. On the other
!!  hand, it is important to distinguish between primary or secondary
!!  mirror planes, and tertiary mirror planes. So we will keep the
!!  following sets : m non-tertiary (that includes a or b non-tertiary) -> 15,
!!  c non-tertiary (that includes n non-tertiary) -> 16,
!!  m tertiary (that includes a or b non-tertiary) -> 17,
!!  c tertiary (that includes n non-tertiary) -> 18.
!!  For hR, all mirror planes are secondary.
!! (4) For the cP lattice, in the same spirit, one can see that
!!  the tertiary m and g mirror planes are to be classified as "m" -> 15,
!!  while n, a and c are to be classified as "n" -> 18. There is no need
!!  to distinguish between primary, secondary or tertiary axes.
!!
!! PARENTS
!!      symcharac
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symplanes(center,iholohedry,isym,isymrelconv,itnonsconv,label,type_axis)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symplanes'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: center,iholohedry,isym
 integer,intent(out) :: type_axis
 character(len = 128), intent(out) :: label
!arrays
 integer,intent(in) :: isymrelconv(3,3)
 real(dp),intent(in) :: itnonsconv(3)

!Local variables-------------------------------
!scalars
 logical,parameter :: verbose=.FALSE.
 character(len=500) :: message
 integer :: directiontype,sum_elements
 real(dp),parameter :: nzero=1.0d-6
!arrays
 integer :: identity(3,3),mirrormxy(3,3),mirrormyz(3,3),mirrormzx(3,3)
 integer :: mirrorx(3,3),mirrorxy(3,3),mirrory(3,3),mirroryz(3,3),mirrorz(3,3)
 integer :: mirrorzx(3,3)
 real(dp) :: trialt(3)
! real(dp) :: itnonsconv2(3),trialt2(3)

!**************************************************************************

!DEBUG
!write(std_out,*)' symplanes : enter'
!write(std_out,*)' center,iholohedry,isym,isymrelconv,itnonsconv=',center,iholohedry,isym,isymrelconv,itnonsconv
!stop
!ENDDEBUG

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1


!Will be a mirror plane, but one must characterize
!(1) the type of plane (primary, secondary or tertiary)
!(2) the gliding vector. One now defines a few matrices.
 mirrorx(:,:)=identity(:,:) ; mirrorx(1,1)=-1
 mirrory(:,:)=identity(:,:) ; mirrory(2,2)=-1
 mirrorz(:,:)=identity(:,:) ; mirrorz(3,3)=-1
 mirrorxy(:,:)=0 ; mirrorxy(1,2)=1 ; mirrorxy(2,1)=1 ; mirrorxy(3,3)=1
 mirrorzx(:,:)=0 ; mirrorzx(1,3)=1 ; mirrorzx(3,1)=1 ; mirrorzx(2,2)=1
 mirroryz(:,:)=0 ; mirroryz(2,3)=1 ; mirroryz(3,2)=1 ; mirroryz(1,1)=1
 mirrormxy(:,:)=0 ; mirrormxy(1,2)=-1 ; mirrormxy(2,1)=-1 ; mirrormxy(3,3)=1
 mirrormzx(:,:)=0 ; mirrormzx(1,3)=-1 ; mirrormzx(3,1)=-1 ; mirrormzx(2,2)=1
 mirrormyz(:,:)=0 ; mirrormyz(2,3)=-1 ; mirrormyz(3,2)=-1 ; mirrormyz(1,1)=1

!Determine the type of plane. At the end,
!directiontype=1 will correspond to a primary axis (or equivalent
!axes for orthorhombic)
!directiontype=2 will correspond to a secondary axis
!directiontype=3 will correspond to a tertiary axis
!See table 2.4.1, 11.2 and 11.3 of the international tables for crystallography
 directiontype=0
!The sum of elements of the matrices allow to characterize them
 sum_elements=sum(isymrelconv(:,:))

 if(sum_elements==1)then
!  The mirror plane perpendicular to the c axis is always primary
   if( sum(abs(isymrelconv(:,:)-mirrorz(:,:)))==0 )then
     directiontype=1
!    All the other planes with a symrel matrix whose sum of elements is 1
!    are a or b planes. They are primary or
!    secondary planes, depending the holohedry.
   else if(sum(isymrelconv(:,:))==1)then
     if( iholohedry==2 .or. iholohedry==3 .or. iholohedry==7 )then
       directiontype=1
     else if(iholohedry==4 .or. iholohedry==6)then
       directiontype=2
     end if
   end if
 end if

!All the planes with a symrel matrix whose sum of elements
!is 2 are secondary planes (table 11.3).
 if( sum_elements==2 ) directiontype=2

!The planes with a symrel matrix whose sum of elements
!is 3 or 0 are tertiary planes
 if( sum_elements==3 .or. sum_elements==0 )directiontype=3

!One is left with sum_elements=-1, tertiary for tetragonal
!or cubic, secondary for hexagonal
 if( sum_elements==-1)then
   if(iholohedry==4 .or. iholohedry==7)directiontype=3
   if(iholohedry==6)directiontype=2
 end if


!Now, determine the gliding vector
!First, apply the symmetry operation
!to itnonsconv, in order to get the translation vector
!under the application of twice the symmetry operation
 trialt(:)=matmul(isymrelconv(:,:),itnonsconv(:)) +itnonsconv(:)
!Get the translation associated with one application,
!and force its components to be in the interval ]-0.5,0.5] .
 trialt(:)=trialt(:)*half
 trialt(:)=trialt(:)-nint(trialt(:)-nzero)

!If there is a glide vector for the initial choice of itnonsconv,
!it might be that it disappears if itnonsconv is translated by a 
!lattice vector of the conventional cell
!if(trialt(1)**2+trialt(2)**2+trialt(3)**2>tol5)then
!do ii=1,3  
!itnonsconv2(:)=itnonsconv(:)
!itnonsconv2(ii)=itnonsconv(ii)+one
!trialt2(:)=matmul(isymrelconv(:,:),itnonsconv2(:)) +itnonsconv2(:)
!trialt2(:)=trialt2(:)*half
!trialt2(:)=trialt2(:)-nint(trialt2(:)-nzero)
!if(trialt2(1)**2+trialt2(2)**2+trialt2(3)**2<tol5)then
!trialt(:)=trialt2(:)
!endif
!enddo
!endif

 write(message,'(a)') ' symplanes...'

!Must use the convention of table 1.3 of the international
!tables for crystallography, see also pp 788 and 789.
!Often, one needs to specialize the selection according
!to the Bravais lattice or the system.

 if(sum(abs(trialt(:)))<nzero .and. iholohedry/=6)then
   type_axis=15  ! m
   write(label,'(a)') 'a mirror plane'
 else if(iholohedry==4 .and. center==0)then    ! primitive tetragonal

   if(directiontype==1)then
     type_axis=18  ! primary n
     write(label,'(a)') 'a primary n plane'
   else if(directiontype==2)then
     if(sum(abs(trialt(:)-(/half,zero,zero/)))<nzero .or. &
&     sum(abs(trialt(:)-(/zero,half,zero/)))<nzero       )then
       type_axis=16  ! secondary a or b
       write(label,'(a)') 'a secondary a or b plane'
     else if(sum(abs(trialt(:)-(/zero,zero,half/)))<nzero)then
       type_axis=17    ! secondary c
       write(label,'(a)') 'a secondary c plane'
     else
       type_axis=18    ! secondary n
       write(label,'(a)') 'a secondary n plane'
     end if ! directiontype==2
   else if(directiontype==3)then
     if( abs(trialt(3))<nzero )then
       type_axis=15    ! tertiary m
       write(label,'(a)') 'a tertiary m plane'
     else if( abs(trialt(3)-half)<nzero )then
       type_axis=19    ! tertiary c
       write(label,'(a)') 'a tertiary c plane'
     end if
   end if

 else if(iholohedry==4 .and. center==-1)then    ! inner tetragonal

   if(directiontype==1 .or. directiontype==2)then
     if(sum(abs(trialt(:)-(/half,zero,zero/)))<nzero .or. &
&     sum(abs(trialt(:)-(/zero,half,zero/)))<nzero .or. &
&     sum(abs(trialt(:)-(/zero,zero,half/)))<nzero      )then
       type_axis=16    ! a, b, or c
       write(label,'(a)') 'an a, b or c plane'
     else if(sum(abs(trialt(:)-(/half,half,zero/)))<nzero .or. &
&       sum(abs(trialt(:)-(/zero,half,half/)))<nzero .or. &
&       sum(abs(trialt(:)-(/half,zero,half/)))<nzero       )then
       type_axis=15    ! n plane, equivalent to m
       write(label,'(a)') 'a m plane'
     end if ! directiontype==1 or 2
   else if(directiontype==3)then
     if( abs(trialt(3))<nzero .or. abs(trialt(3)-half)<nzero )then
       type_axis=15    ! tertiary c, equivalent to m
       write(label,'(a)') 'a tertiary m plane'
     else
       type_axis=17    ! tertiary d
       write(label,'(a)') 'a tertiary d plane'
     end if
   end if

 else if(iholohedry==5)then    ! hR

   if( abs(sum(abs(trialt(:)))-one) < nzero) then
     type_axis=15    ! secondary m
     write(label,'(a)') 'a secondary m plane'
   else if( abs(sum(abs(trialt(:)))-half) < nzero .or. &
&     abs(sum(abs(trialt(:)))-three*half) < nzero )then
     type_axis=16    ! secondary c
     write(label,'(a)') 'a secondary c plane'
   end if

 else if(iholohedry==6)then    ! hP

   if(directiontype==1)then
     if( abs(trialt(3)) <nzero )then
       type_axis=15    ! primary m
       write(label,'(a)') 'a primary m plane'
     end if
   else if(directiontype==2)then
     if( abs(trialt(3)) <nzero )then
       type_axis=15    ! secondary m
       write(label,'(a)') 'a secondary m plane'
     else if( abs(trialt(3)-half) < nzero ) then
       type_axis=16    ! secondary c
       write(label,'(a)') 'a secondary c plane'
     end if
   else if(directiontype==3)then
     if( abs(trialt(3)) <nzero )then
       type_axis=17    ! tertiary m
       write(label,'(a)') 'a tertiary m plane'
     else if( abs(trialt(3)-half) < nzero ) then
       type_axis=18    ! tertiary c
       write(label,'(a)') 'a tertiary c plane'
     end if
   end if ! directiontype

!  else if(iholohedry==7 .and. center==0)then    ! cP
 else if(iholohedry==7)then    ! cP

   if(directiontype==1)then
     if((sum(abs(isymrelconv(:,:)-mirrorx(:,:)))==0 .and.  &
&     sum(abs(two*abs(trialt(:))-(/zero,half,half/)))<nzero   ).or. &
&     (sum(abs(isymrelconv(:,:)-mirrory(:,:)))==0 .and.  &
&     sum(abs(two*abs(trialt(:))-(/half,zero,half/)))<nzero   ).or. &
&     (sum(abs(isymrelconv(:,:)-mirrorz(:,:)))==0 .and.  &
&     sum(abs(two*abs(trialt(:))-(/half,half,zero/)))<nzero   )    ) then
       type_axis=17     ! d
       write(label,'(a)') 'a d plane'
     else
       type_axis=18    ! primary n
       write(label,'(a)') 'a primary n plane'
     end if
   else if(directiontype==3)then
     if(sum(abs(two*abs(trialt(:))-(/half,half,half/)))<nzero       )then
       type_axis=17     ! d
       write(label,'(a)') 'a d plane'
     else if( abs(sum(abs(trialt(:)))-half) < nzero .or. &
&       abs(sum(abs(trialt(:)))-three*half) < nzero ) then
       type_axis=18    ! tertiary n
       write(label,'(a)') 'a tertiary n plane'
     else if( abs(sum(abs(trialt(:)))-one) < nzero )then
       type_axis=15    ! tertiary m
       write(label,'(a)') 'a tertiary m plane'
     end if
   end if 

!  Now, treat all other cases (including other centered Bravais lattices)
 else if( sum(abs(trialt(:)-(/half,zero,zero/)))<nzero .or. &
&   sum(abs(trialt(:)-(/zero,half,zero/)))<nzero .or. &
&   sum(abs(trialt(:)-(/zero,zero,half/)))<nzero       )then
   type_axis=16     ! a, b or c
   write(label,'(a)') 'an a,b, or c plane'
 else if( (directiontype==1 .or. directiontype==2) .and. &
&   (sum(abs(trialt(:)-(/half,half,zero/)))<nzero .or. &
&   sum(abs(trialt(:)-(/zero,half,half/)))<nzero .or. &
&   sum(abs(trialt(:)-(/half,zero,half/)))<nzero     ) )then
   type_axis=18     ! n
   write(label,'(a)') 'an n plane'
 else if( directiontype==3 .and. &
&   sum(abs(trialt(:)-(/half,half,half/)))<nzero )then
   type_axis=18     ! n
   write(label,'(a)') 'an n plane'
 else if((sum(abs(isymrelconv(:,:)-mirrorx(:,:)))==0 .and.  &
&   sum(abs(two*abs(trialt(:))-(/zero,half,half/)))<nzero   ).or. &
&   (sum(abs(isymrelconv(:,:)-mirrory(:,:)))==0 .and.  &
&   sum(abs(two*abs(trialt(:))-(/half,zero,half/)))<nzero   ).or. &
&   (sum(abs(isymrelconv(:,:)-mirrorz(:,:)))==0 .and.  &
&   sum(abs(two*abs(trialt(:))-(/half,half,zero/)))<nzero   )    ) then
   type_axis=17     ! d
   write(label,'(a)') 'a d plane'
 else if( directiontype==3 .and. &
&   sum(abs(two*abs(trialt(:))-(/half,half,half/)))<nzero       )then
   type_axis=17     ! d
   write(label,'(a)') 'a d plane'
 else
   type_axis=19     ! g (all other planes with
!  unconventional glide vector)
   write(label,'(a)') 'a g plane'
 end if

 if (verbose) then
   write(message,'(a,i3,a,a)')' symplanes : the symmetry operation no. ',isym,' is ', trim(label)
   call wrtout(std_out,message,'COLL')
 end if

!call symlist(brvltt,nsymconv,n_axes,problem)

end subroutine symplanes
!!***

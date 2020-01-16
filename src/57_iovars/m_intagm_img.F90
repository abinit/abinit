!!****m* ABINIT/m_intagm_img
!!
!! NAME
!! m_intagm_img
!!
!! FUNCTION
!! This module provides a generic interface that allows to
!! initialize some of the geometry variables in the case of "images".
!! Set up: acell, scalecart, rprim, angdeg, xred, xangst, xcart, vel
!! These variables can be defined for a set of images of the cell.
!! They also can be be defined along a path (in the configuration space).
!! The path must be defined with its first and last points, but also
!! with intermediate points.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2019 ABINIT group (XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iimage=index of the current image
!!  jdtset=number of the dataset looked for
!!  lenstr=actual length of the input string
!!  nimage=number of images
!!  size1,size2, ...: size of array to be read (dp_data)
!!  string=character string containing 'tags' and data.
!!  token=character string for tagging the data to be read in input string
!!  typevarphys= variable type (for dimensionality purposes)
!!
!! SIDE EFFECTS
!!  dp_data(size1,size2,...)=data to be read (double precision)
!!  tread_ok=flag to be set to 1 if the data have been found in input string
!!
!! NOTES
!! The routine is a generic interface calling subroutine according to the
!! number of arguments of the variable to be read
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_intagm_img

 use defs_basis
 use m_abicore
 use m_errors

 use m_parser,    only : intagm

 implicit none

 private
!!***

 public :: intagm_img   !  Read input file variables according to images path definition (1D array)

interface intagm_img
  module procedure intagm_img_1D
  module procedure intagm_img_2D
end interface intagm_img

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_intagm_img/ingeo_img_1D
!! NAME
!!  intagm_img_1D
!!
!! FUNCTION
!!  Read input file variables according to images path definition (1D array)
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

subroutine intagm_img_1D(dp_data,iimage,jdtset,lenstr,nimage,size1,string,token,tread_ok,typevarphys)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iimage,jdtset,lenstr,nimage,size1
 integer,intent(inout) :: tread_ok
 real(dp),intent(inout) :: dp_data(size1)
 character(len=*),intent(in) :: typevarphys
 character(len=*),intent(in) :: token
 character(len=*),intent(in) :: string
!arrays

!Local variables-------------------------------
!scalars
 integer :: iimage_after,iimage_before,marr,tread_after,tread_before,tread_current
 real(dp) :: alpha
 character(len=10) :: stringimage
 character(len=3*len(token)+10) :: token_img
!arrays
 integer, allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),dp_data_after(:),dp_data_before(:)

! *************************************************************************

!Nothing to do in case of a single image
 if (nimage<=1) return

 marr=size1
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!First, try to read data for current image
 tread_current=0
 write(stringimage,'(i10)') iimage
 token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
 call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&            token_img,tread_current,typevarphys)
 if (tread_current==1)then
   dp_data(1:size1)=dprarr(1:size1)
   tread_ok=1
 end if
 if (tread_current==0.and.iimage==nimage) then
!  If the image is the last one, try to read data for last image (_lastimg)
   token_img=trim(token)//'_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&              token_img,tread_current,typevarphys)
   if (tread_current==1)then
     dp_data(1:size1)=dprarr(1:size1)
     tread_ok=1
   end if
 end if

 if (tread_current==0) then

!  The current image is not directly defined in the input string
   ABI_ALLOCATE(dp_data_before,(size1))
   ABI_ALLOCATE(dp_data_after,(size1))

!  Find the nearest previous defined image
   tread_before=0;iimage_before=iimage
   do while (iimage_before>1.and.tread_before/=1)
     iimage_before=iimage_before-1
     write(stringimage,'(i10)') iimage_before
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&                token_img,tread_before,typevarphys)
     if (tread_before==1) dp_data_before(1:size1)=dprarr(1:size1)
   end do
   if (tread_before==0) then
     iimage_before=1
     dp_data_before(1:size1)=dp_data(1:size1)
   end if

!  Find the nearest following defined image
   tread_after=0;iimage_after=iimage
   do while (iimage_after<nimage.and.tread_after/=1)
     iimage_after=iimage_after+1
     write(stringimage,'(i10)') iimage_after
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&                token_img,tread_after,typevarphys)
     if (tread_after==1) dp_data_after(1:size1)=dprarr(1:size1)
     if (tread_after==0.and.iimage_after==nimage) then
       token_img=trim(token)//'_lastimg'
       call intagm(dprarr,intarr,jdtset,marr,size1,string(1:lenstr),&
&                  token_img,tread_after,typevarphys)
       if (tread_after==1) dp_data_after(1:size1)=dprarr(1:size1)
     end if
   end do
   if (tread_after==0) then
     iimage_after=nimage
     dp_data_after(1:size1)=dp_data(1:size1)
   end if

!  Interpolate image data
   if (tread_before==1.or.tread_after==1) then
     alpha=real(iimage-iimage_before,dp)/real(iimage_after-iimage_before,dp)
     dp_data(1:size1)=dp_data_before(1:size1) &
&                    +alpha*(dp_data_after(1:size1)-dp_data_before(1:size1))
     tread_ok=1
   end if

   ABI_DEALLOCATE(dp_data_before)
   ABI_DEALLOCATE(dp_data_after)

 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine intagm_img_1D
!!***


!----------------------------------------------------------------------


!!****f* m_intagm_img/ingeo_img_2D
!! NAME
!!  intagm_img_2D
!!
!! FUNCTION
!!  Read input file variables according to images path definition (2D array)
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

subroutine intagm_img_2D(dp_data,iimage,jdtset,lenstr,nimage,size1,size2,string,token,tread_ok,typevarphys)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iimage,jdtset,lenstr,nimage,size1,size2
 integer,intent(inout) :: tread_ok
 real(dp),intent(inout) :: dp_data(size1,size2)
 character(len=*),intent(in) :: typevarphys
 character(len=*),intent(in) :: token
 character(len=*),intent(in) :: string
!arrays

!Local variables-------------------------------
!scalars
 integer :: iimage_after,iimage_before,marr,tread_after,tread_before,tread_current
 real(dp) :: alpha
 character(len=10) :: stringimage
 character(len=3*len(token)+10) :: token_img
!arrays
 integer, allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),dp_data_after(:,:),dp_data_before(:,:)

! *************************************************************************

!Nothing to do in case of a single image
 if (nimage<=1) return

 marr=size1*size2
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!First, try to read data for current image
 tread_current=0
 write(stringimage,'(i10)') iimage
 token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
 call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&            token_img,tread_current,typevarphys)
 if (tread_current==1)then
   dp_data(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
   tread_ok=1
 end if
 if (tread_current==0.and.iimage==nimage) then
!  In the image is the last one, try to read data for last image (_lastimg)
   token_img=trim(token)//'_lastimg'
   call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&              token_img,tread_current,typevarphys)
   if (tread_current==1)then
     dp_data(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
     tread_ok=1
   end if
 end if

 if (tread_current==0) then

!  The current image is not directly defined in the input string
   ABI_ALLOCATE(dp_data_before,(size1,size2))
   ABI_ALLOCATE(dp_data_after,(size1,size2))

!  Find the nearest previous defined image
   tread_before=0;iimage_before=iimage
   do while (iimage_before>1.and.tread_before/=1)
     iimage_before=iimage_before-1
     write(stringimage,'(i10)') iimage_before
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&                token_img,tread_before,typevarphys)
     if (tread_before==1) &
&      dp_data_before(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
   end do
   if (tread_before==0) then
     iimage_before=1
     dp_data_before(1:size1,1:size2)=dp_data(1:size1,1:size2)
   end if

!  Find the nearest following defined image
   tread_after=0;iimage_after=iimage
   do while (iimage_after<nimage.and.tread_after/=1)
     iimage_after=iimage_after+1
     write(stringimage,'(i10)') iimage_after
     token_img=trim(token)//'_'//trim(adjustl(stringimage))//'img'
     call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&                token_img,tread_after,typevarphys)
     if (tread_after==1) &
&      dp_data_after(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
     if (tread_after==0.and.iimage_after==nimage) then
       token_img=trim(token)//'_lastimg'
       call intagm(dprarr,intarr,jdtset,marr,size1*size2,string(1:lenstr),&
&                  token_img,tread_after,typevarphys)
       if (tread_after==1) &
&        dp_data_after(1:size1,1:size2)=reshape( dprarr(1:size1*size2),(/size1,size2/) )
     end if
   end do
   if (tread_after==0) then
     iimage_after=nimage
     dp_data_after(1:size1,1:size2)=dp_data(1:size1,1:size2)
   end if

!  Interpolate image data
   if (tread_before==1.or.tread_after==1) then
     alpha=real(iimage-iimage_before,dp)/real(iimage_after-iimage_before,dp)
     dp_data(1:size1,1:size2)=dp_data_before(1:size1,1:size2) &
&       +alpha*(dp_data_after(1:size1,1:size2)-dp_data_before(1:size1,1:size2))
     tread_ok=1
   end if

   ABI_DEALLOCATE(dp_data_before)
   ABI_DEALLOCATE(dp_data_after)

 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine intagm_img_2D
!!***

!------------------------------------------------------------------------------------

END MODULE m_intagm_img
!!***

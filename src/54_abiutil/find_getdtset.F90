!{\src2tex{textfont=tt}}
!!****f* ABINIT/find_getdtset
!! NAME
!! find_getdtset
!!
!! FUNCTION
!! Find the number of the dataset (iget) for a given value of a "get" variable (getvalue)
!! of name getname, given the number of the current dataset (idtset).
!! Also find the coefficients of mixing of the images of the old dataset, to initialize the new dataset images
!!     (use a linear interpolation)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtsets(0:ndtset_alloc)=<type datasets_type>contains all input variables
!! getvalue=value of the get variable
!! getname=name of the get variable
!! idtset=number of the current dataset
!! mxnimage=dimension of miximage
!! ndtset_alloc=dimension of dtsets
!!
!! OUTPUT
!! iget=number of the dataset from which the value must be get, 0 if the data should not be got from another dataset
!! miximage(mxnimage,mxnimage)=coefficients of mixing of the images of the old dataset, to initialize the new dataset images
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine find_getdtset(dtsets,getvalue,getname,idtset,iget,miximage,mxnimage,ndtset_alloc)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_getdtset'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: getvalue,idtset,mxnimage,ndtset_alloc
 integer, intent(out) :: iget
 real(dp), intent(out) :: miximage(mxnimage,mxnimage)
 character(len=*),intent(in) :: getname
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
 integer :: iimage
 real(dp) :: newimage_get,ratio
 character(len=500) :: message

! *************************************************************************

 iget=0
 if(getvalue>0 .or. (getvalue<0 .and. idtset+getvalue>0) )then
!  In case getvalue is a negative number (so must add to idtset)
   if(getvalue<0 .and. idtset+getvalue>0) iget=idtset+getvalue
   if(getvalue>0)then
     do iget=1,idtset
       if( dtsets(iget)%jdtset==getvalue )exit
     end do
     if(iget==idtset)then
!      The index of the dataset, from which the data ought to be taken,
!      does not correspond to a previous dataset.
       write(message, '(a,i3,4a,i3,7a)' )&
&       'The component number',idtset,' of the input variable ',trim(getname),',',&
&       ' equal to',getvalue,',',ch10,&
&       'does not correspond to an existing index.',ch10,&
&       'Action : correct ',trim(getname),' or jdtset in your input file.'
       MSG_ERROR(message)
     end if
   end if
   write(message, '(3a,i3,2a)' )&
&   ' find_getdtset : ',trim(getname),'/=0, take data from output of dataset with index',dtsets(iget)%jdtset,'.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!For the time being, uses a simple interpolation when the images do not match. If only one image, take the first get image.
 miximage(:,:)=zero
 if(dtsets(idtset)%nimage==1)then
   miximage(1,1)=one
 else
   do iimage=1,dtsets(idtset)%nimage
     ratio=(iimage-one)/real(dtsets(idtset)%nimage-one)
     newimage_get=one+ratio*(dtsets(iget)%nimage-one)
     if(abs(newimage_get-nint(newimage_get))<tol8)then
       miximage(iimage,nint(newimage_get))=one
     else
       miximage(iimage,floor(newimage_get))=one-(newimage_get-floor(newimage_get))
       miximage(iimage,ceiling(newimage_get))=one-miximage(iimage,floor(newimage_get))
     end if
   end do
 end if

end subroutine find_getdtset
!!***

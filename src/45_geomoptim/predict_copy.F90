!{\src2tex{textfont=tt}}
!!****f* ABINIT/predict_copy
!! NAME
!! predict_copy
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images.
!! Here, simple copy of the previous image.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
!! ndynimage=number of dynamical images
!! nimage=number of images
!! ntimimage=dimension of several arrays
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage,nimage)=datastructure that hold all the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images, up to itimimage
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images, up to itimimage
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images, up to itimimage
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images, up to itimimage
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images, up to itimimage
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images, up to itimimage
!!    at output, the predicted values of xred for all images
!!
!! PARENTS
!!      predictimg
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine predict_copy(itimimage,list_dynimage,ndynimage,nimage,ntimimage,results_img)

 use defs_basis
 use m_profiling_abi

 use m_results_img, only : results_img_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predict_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,ndynimage,nimage,ntimimage
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type),intent(inout) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer :: idynimage,iimage
!arrays

! *************************************************************************

 do idynimage=1,ndynimage

   iimage=list_dynimage(idynimage)

   results_img(iimage,itimimage+1)%acell(:)     =results_img(iimage,itimimage)%acell(:)
   results_img(iimage,itimimage+1)%rprim(:,:)   =results_img(iimage,itimimage)%rprim(:,:)
   results_img(iimage,itimimage+1)%vel(:,:)     =results_img(iimage,itimimage)%vel(:,:)
   results_img(iimage,itimimage+1)%vel_cell(:,:)=results_img(iimage,itimimage)%vel_cell(:,:)
   results_img(iimage,itimimage+1)%xred(:,:)    =results_img(iimage,itimimage)%xred(:,:)

 end do  ! idynimage

end subroutine predict_copy
!!***

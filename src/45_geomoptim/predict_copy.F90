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
!! Copyright (C) 2009-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
!! ndynimage=number of dynamical images
!! nimage=number of images
!! ntimimage_stored=number of time steps stored in the history
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage_stored,nimage)=datastructure that holds the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images
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


subroutine predict_copy(itimimage_eff,list_dynimage,ndynimage,nimage,&
&                       ntimimage_stored,results_img)

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
 integer,intent(in) :: itimimage_eff,ndynimage,nimage,ntimimage_stored
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type),intent(inout) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: idynimage,iimage,next_itimimage
!arrays

! *************************************************************************

 next_itimimage=itimimage_eff+1
 if (next_itimimage>ntimimage_stored) next_itimimage=1

 do idynimage=1,ndynimage

   iimage=list_dynimage(idynimage)

   results_img(iimage,next_itimimage)%acell(:)     =results_img(iimage,itimimage_eff)%acell(:)
   results_img(iimage,next_itimimage)%rprim(:,:)   =results_img(iimage,itimimage_eff)%rprim(:,:)
   results_img(iimage,next_itimimage)%vel(:,:)     =results_img(iimage,itimimage_eff)%vel(:,:)
   results_img(iimage,next_itimimage)%vel_cell(:,:)=results_img(iimage,itimimage_eff)%vel_cell(:,:)
   results_img(iimage,next_itimimage)%xred(:,:)    =results_img(iimage,itimimage_eff)%xred(:,:)

 end do  ! idynimage

end subroutine predict_copy
!!***

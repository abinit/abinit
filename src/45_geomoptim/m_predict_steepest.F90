!!****m* ABINIT/m_predict_steepest
!! NAME
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_predict_steepest

 use defs_basis
 use m_abicore
 use m_mep

 use m_results_img, only : results_img_type,get_geometry_img

 implicit none

 private
!!***

 public :: predict_steepest
!!***

contains
!!***

!!****f* ABINIT/predict_steepest
!! NAME
!! predict_steepest
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images.
!! Here, simple steepest descent algorithm, based on the value of the forces on the current timimage step.
!! No change of acell, rprim and vel at present.
!!
!! INPUTS
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
!! mep_param=several parameters for Minimal Energy Path (MEP) search
!! natom=dimension of vel_timimage and xred_timimage
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
!!      get_geometry_img,mep_steepest
!!
!! SOURCE

subroutine predict_steepest(itimimage,itimimage_eff,list_dynimage,mep_param,natom,&
&                           ndynimage,nimage,ntimimage_stored,results_img)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,itimimage_eff,natom,ndynimage
 integer,intent(in) :: nimage,ntimimage_stored
 type(mep_type),intent(inout) :: mep_param
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: iimage,next_itimimage
!arrays
 real(dp),allocatable :: etotal(:),fcart(:,:,:),rprimd(:,:,:),xcart(:,:,:),xred(:,:,:)

! *************************************************************************

!Retrieve positions and forces
 ABI_ALLOCATE(etotal,(nimage))
 ABI_ALLOCATE(xred,(3,natom,nimage))
 ABI_ALLOCATE(xcart,(3,natom,nimage))
 ABI_ALLOCATE(fcart,(3,natom,nimage))
 ABI_ALLOCATE(rprimd,(3,3,nimage))
 call get_geometry_img(etotal,natom,nimage,results_img(:,itimimage_eff),&
& fcart,rprimd,xcart,xred)

!Compute new atomic positions in each cell
 call mep_steepest(fcart,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)

!Store acell, rprim, xred and vel for the new iteration
 next_itimimage=itimimage+1
 if (next_itimimage>ntimimage_stored) next_itimimage=1
 do iimage=1,nimage
   results_img(iimage,next_itimimage)%xred(:,:)    =xred(:,:,iimage)
   results_img(iimage,next_itimimage)%acell(:)     =results_img(iimage,itimimage_eff)%acell(:)
   results_img(iimage,next_itimimage)%rprim(:,:)   =results_img(iimage,itimimage_eff)%rprim(:,:)
   results_img(iimage,next_itimimage)%vel(:,:)     =results_img(iimage,itimimage_eff)%vel(:,:)
   results_img(iimage,next_itimimage)%vel_cell(:,:)=results_img(iimage,itimimage_eff)%vel_cell(:,:)
 end do

 ABI_DEALLOCATE(etotal)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(rprimd)

end subroutine predict_steepest
!!***

end module m_predict_steepest
!!***

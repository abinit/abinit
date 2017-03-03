!{\src2tex{textfont=tt}}
!!****f* ABINIT/predict_steepest
!! NAME
!! predict_steepest
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images.
!! Here, simple steepest descent algorithm, based on the value of the forces on the current timimage step.
!! No change of acell, rprim and vel at present.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
!! mep_param=several parameters for Minimal Energy Path (MEP) search
!! natom=dimension of vel_timimage and xred_timimage
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
!!      get_geometry_img,mep_steepest
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine predict_steepest(itimimage,list_dynimage,mep_param,natom,ndynimage,nimage,&
&                           ntimimage,results_img)

 use m_profiling_abi

 use defs_basis
 use m_results_img, only : results_img_type,get_geometry_img
 use m_mep

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predict_steepest'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,ndynimage,nimage,ntimimage
 type(mep_type),intent(inout) :: mep_param
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer :: iimage
!arrays
 real(dp),allocatable :: etotal(:),fcart(:,:,:),rprimd(:,:,:),xcart(:,:,:),xred(:,:,:)

! *************************************************************************

!Retrieve positions and forces
 ABI_ALLOCATE(etotal,(nimage))
 ABI_ALLOCATE(xred,(3,natom,nimage))
 ABI_ALLOCATE(xcart,(3,natom,nimage))
 ABI_ALLOCATE(fcart,(3,natom,nimage))
 ABI_ALLOCATE(rprimd,(3,3,nimage))
 call get_geometry_img(etotal,natom,nimage,results_img(:,itimimage),fcart,rprimd,xcart,xred)

!Compute new atomic positions in each cell
 call mep_steepest(fcart,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)

!Store acell, rprim, xred and vel for the new iteration
 do iimage=1,nimage
   results_img(iimage,itimimage+1)%xred(:,:)    =xred(:,:,iimage)
   results_img(iimage,itimimage+1)%acell(:)     =results_img(iimage,itimimage)%acell(:)
   results_img(iimage,itimimage+1)%rprim(:,:)   =results_img(iimage,itimimage)%rprim(:,:)
   results_img(iimage,itimimage+1)%vel(:,:)     =results_img(iimage,itimimage)%vel(:,:)
   results_img(iimage,itimimage+1)%vel_cell(:,:)=results_img(iimage,itimimage)%vel_cell(:,:)
 end do

 ABI_DEALLOCATE(etotal)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(rprimd)

end subroutine predict_steepest
!!***

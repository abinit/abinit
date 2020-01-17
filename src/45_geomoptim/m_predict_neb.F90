!!****m* ABINIT/m_predict_neb
!! NAME
!!  m_predict_neb
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2020 ABINIT group (MT)
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

module m_predict_neb

 use defs_basis
 use m_splines
 use m_mep
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_results_img, only : results_img_type, gather_array_img, scatter_array_img, get_geometry_img

 implicit none

 private
!!***

 public :: predict_neb
!!***

contains
!!***

!!****f* ABINIT/predict_neb
!! NAME
!! predict_neb
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images using Nudged Elastic Band method.
!! No change of acell, rprim and vel at present.
!!
!! INPUTS
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example: in the NEB of string method, one expect the two end images to be fixed.
!! This is quite useful when ground states of the A and B states is known
!! mpi_enreg=MPI-parallelisation information
!! natom=dimension of vel_timimage and xred_timimage
!! ndynimage=number of dynamical images
!! nimage=number of images (on current proc)
!! nimage_tot=total number of images
!! ntimimage_stored=number of time steps stored in the history
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! mep_param=several parameters for Minimal Energy Path (MEP) search
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
!!      gather_array_img,get_geometry_img,mep_gbfgs,mep_lbfgs,mep_qmin
!!      mep_steepest,xmpi_bcast
!!
!! SOURCE

subroutine predict_neb(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&                      ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,itimimage_eff,natom,ndynimage
 integer,intent(in) :: nimage,nimage_tot,ntimimage_stored
 type(mep_type),intent(inout) :: mep_param
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in)     :: list_dynimage(ndynimage)
 type(results_img_type),intent(inout) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: ierr,ii,iimage,iimage_min,iimage_max,next_itimimage
 real(dp) :: dvmax,dvmin,ecur,emax,emin,eref,f_para1,f_para2
 logical :: test_minus_one,test_plus_one
!arrays
 real(dp),allocatable :: buffer(:,:),buffer_all(:,:)
 real(dp),allocatable :: coordif(:,:,:),dimage(:),rprimd(:,:,:),spring(:)
 real(dp),allocatable :: tangent(:,:,:),vect(:,:)
 real(dp),allocatable,target :: etotal(:),fcart(:,:,:),neb_forces(:,:,:)
 real(dp),allocatable,target :: xcart(:,:,:),xred(:,:,:)
 real(dp),pointer :: etotal_all(:),fcart_all(:,:,:),neb_forces_all(:,:,:)
 real(dp),pointer :: xcart_all(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(xred,(3,natom,nimage))

!Parallelism over images: only one process per image of the cell
 if (mpi_enreg%me_cell==0) then

!  Retrieve positions and forces
   ABI_ALLOCATE(etotal,(nimage))
   ABI_ALLOCATE(xcart,(3,natom,nimage))
   ABI_ALLOCATE(fcart,(3,natom,nimage))
   ABI_ALLOCATE(rprimd,(3,3,nimage))
   call get_geometry_img(etotal,natom,nimage,results_img(:,itimimage_eff),&
&   fcart,rprimd,xcart,xred)

!  Array containing effective NEB forces (F_ortho+F_spring)
   ABI_ALLOCATE(neb_forces,(3,natom,nimage))

!  Parallelism: gather data of all images
   if (mpi_enreg%paral_img==1) then
     ABI_ALLOCATE(buffer,(6*natom+1,nimage))
     ABI_ALLOCATE(buffer_all,(6*natom+1,nimage_tot))
     ABI_ALLOCATE(xcart_all,(3,natom,nimage_tot))
     ABI_ALLOCATE(fcart_all,(3,natom,nimage_tot))
     ABI_ALLOCATE(etotal_all,(nimage_tot))
     ABI_ALLOCATE(neb_forces_all,(3,natom,nimage_tot))
     buffer=zero;ii=0
     buffer(ii+1:ii+3*natom,1:nimage)=reshape(xcart,(/3*natom,nimage/));ii=3*natom
     buffer(ii+1:ii+3*natom,1:nimage)=reshape(fcart,(/3*natom,nimage/));ii=6*natom
     buffer(ii+1           ,1:nimage)=etotal(1:nimage);ii=0
     call gather_array_img(buffer,buffer_all,mpi_enreg,allgather=.true.)
     xcart_all(:,:,:)=reshape(buffer_all(ii+1:ii+3*natom,1:nimage_tot),(/3,natom,nimage_tot/));ii=3*natom
     fcart_all(:,:,:)=reshape(buffer_all(ii+1:ii+3*natom,1:nimage_tot),(/3,natom,nimage_tot/));ii=6*natom
     etotal_all(:)=buffer_all(ii+1,1:nimage_tot);ii=0
     ABI_DEALLOCATE(buffer)
     ABI_DEALLOCATE(buffer_all)
   else
     xcart_all      => xcart
     fcart_all      => fcart
     etotal_all     => etotal
     neb_forces_all => neb_forces
   end if

!  coordif is the vector between two images
!  dimage is the distance between two images
!  tangent is the tangent at image i
   ABI_ALLOCATE(coordif,(3,natom,nimage_tot))
   ABI_ALLOCATE(tangent,(3,natom,nimage_tot))
   ABI_ALLOCATE(dimage,(nimage_tot))

!  Compute distances between images
   coordif(:,:,1)=zero;dimage(1)=zero
   do iimage=2,nimage_tot
     coordif(:,:,iimage)=xcart_all(:,:,iimage)-xcart_all(:,:,iimage-1)
     dimage(iimage)=mep_img_norm(coordif(:,:,iimage))
   end do

!  Compute tangent (not normalized)
   tangent(:,:,1)=zero
   tangent(:,:,nimage_tot)=zero
!  === Original definition of tangent
   if (mep_param%neb_algo==0) then
     do iimage=2,nimage_tot-1
!      tangent(:,:,iimage)=xcart_all(:,:,iimage+1)-xcart_all(:,:,iimage-1)
       tangent(:,:,iimage)=coordif(:,:,iimage  )/dimage(iimage) &
&                         +coordif(:,:,iimage+1)/dimage(iimage+1)
     end do
!    === Improved tangent (J. Chem. Phys. 113, 9978 (2000) [[cite:Henkelman2000]])
   else
     do iimage=2,nimage_tot-1
       test_minus_one=(etotal_all(iimage-1)<etotal_all(iimage))
       test_plus_one =(etotal_all(iimage  )<etotal_all(iimage+1))
       if (test_minus_one.and.test_plus_one)then       ! V_i-1 < V_i < V_i+1
         tangent(:,:,iimage)=coordif(:,:,iimage+1)
       else if ((.not.test_minus_one).and.(.not.test_plus_one))then     ! V_i-1 > V_i > V_i+1
         tangent(:,:,iimage)=coordif(:,:,iimage)
       else                             ! V_i-1 < V_i > V_i+1  OR  V_i-1 > V_i < V_i+1
         dvmax=max(abs(etotal_all(iimage+1)-etotal_all(iimage)),&
&         abs(etotal_all(iimage-1)-etotal_all(iimage)))
         dvmin=min(abs(etotal_all(iimage+1)-etotal_all(iimage)),&
&         abs(etotal_all(iimage-1)-etotal_all(iimage)))
         if (etotal_all(iimage+1)>etotal_all(iimage-1)) then   ! V_i+1 > V_i-1
           tangent(:,:,iimage)=coordif(:,:,iimage+1)*dvmax &
&                             +coordif(:,:,iimage  )*dvmin
         else                                                  ! V_i+1 < V_i-1
           tangent(:,:,iimage)=coordif(:,:,iimage+1)*dvmin &
&                             +coordif(:,:,iimage  )*dvmax
         end if
       end if
     end do
   end if

!  Normalize tangent
   do iimage=2,nimage_tot-1
     tangent(:,:,iimage)=tangent(:,:,iimage)/mep_img_norm(tangent(:,:,iimage))
   end do

!  Compute spring constant(s)
   ABI_ALLOCATE(spring,(nimage_tot))
   if (abs(mep_param%neb_spring(2)-mep_param%neb_spring(1))<tol8) then ! Constant spring
     spring(:)=mep_param%neb_spring(1)
   else                                                                ! Variable spring
     emax=maxval(etotal_all(:))
     eref=max(etotal_all(1),etotal_all(nimage_tot))
     spring(:)=mep_param%neb_spring(1)
     do iimage=2,nimage_tot-1
       ecur=max(etotal_all(iimage-1),etotal_all(iimage))
       if (ecur<eref) then
         spring(iimage)=mep_param%neb_spring(1)
       else
         spring(iimage)=mep_param%neb_spring(2) &
&         -(mep_param%neb_spring(2)-mep_param%neb_spring(1)) &
&         *(emax-ecur)/(emax-eref)
       end if
     end do
   end if

!  CI-NEB: determine image(s) with maximal energy
   iimage_min=-1;iimage_max=-1
   if (mep_param%neb_algo==2.and.itimimage>=mep_param%cineb_start) then
     emin=min(etotal_all(1),etotal_all(nimage_tot))
     emax=max(etotal_all(1),etotal_all(nimage_tot))
     do iimage=2,nimage_tot-1
       if (etotal_all(iimage)<emin) then
         iimage_min=iimage;emin=etotal_all(iimage)
       end if
       if (etotal_all(iimage)>emax) then
         iimage_max=iimage;emax=etotal_all(iimage)
       end if
     end do
   end if

!  Compute NEB forces
   neb_forces_all(:,:,1)=fcart_all(:,:,1)
   neb_forces_all(:,:,nimage_tot)=fcart_all(:,:,nimage_tot)
   do iimage=2,nimage_tot-1
!    === Standard NEB
     if (iimage/=iimage_max) then
!    if (iimage/=iimage_min.and.iimage/=iimage_max) then
       f_para1=mep_img_dotp(fcart_all(:,:,iimage),tangent(:,:,iimage))
       if (mep_param%neb_algo==0) then   ! Original NEB algo
         ABI_ALLOCATE(vect,(3,natom))
         vect(:,:)=spring(iimage+1)*coordif(:,:,iimage+1)-spring(iimage)*coordif(:,:,iimage)
         f_para2=mep_img_dotp(vect,tangent(:,:,iimage))
         ABI_DEALLOCATE(vect)
       else       ! Modification from J. Chem. Phys. 113, 9978 (2000) [[cite:Henkelman2000]]
         f_para2=spring(iimage+1)*dimage(iimage+1)-spring(iimage)*dimage(iimage)
       end if
       neb_forces_all(:,:,iimage)=fcart_all(:,:,iimage) &       ! F_ortho
&      -f_para1*tangent(:,:,iimage) &  ! F_ortho
&      +f_para2*tangent(:,:,iimage)    ! F_spring
!      === CI-NEB for the image(s) with maximal energy
     else
       f_para1=mep_img_dotp(fcart_all(:,:,iimage),tangent(:,:,iimage))
       neb_forces_all(:,:,iimage)=fcart_all(:,:,iimage)-two*f_para1*tangent(:,:,iimage)
     end if
   end do

!  Parallelism: distribute NEB forces
   if (mpi_enreg%paral_img==1) then
     do iimage=1,nimage
       ii=mpi_enreg%my_imgtab(iimage)
       neb_forces(:,:,iimage)=neb_forces_all(:,:,ii)
     end do
   end if

!  Compute new atomic positions in each cell
   if      (mep_param%mep_solver==0) then ! Steepest-descent
     call mep_steepest(neb_forces,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)
   else if (mep_param%mep_solver==1) then ! Quick-min
     call mep_qmin(neb_forces,itimimage,list_dynimage,mep_param,natom,ndynimage, &
&     nimage,rprimd,xcart,xred)
   else if (mep_param%mep_solver==2) then ! Local BFGS
     call mep_lbfgs(neb_forces,itimimage,list_dynimage,mep_param,natom,ndynimage,&
&     nimage,rprimd,xcart,xred)
   else if (mep_param%mep_solver==3) then ! Global BFGS
     call mep_gbfgs(neb_forces,itimimage,list_dynimage,mep_param,mpi_enreg,natom,ndynimage,&
&     nimage,nimage_tot,rprimd,xcart,xred)
   else
     MSG_BUG("Inconsistent solver !")
   end if

!  Free memory
   ABI_DEALLOCATE(spring)
   ABI_DEALLOCATE(coordif)
   ABI_DEALLOCATE(tangent)
   ABI_DEALLOCATE(dimage)
   if (mpi_enreg%paral_img==1)  then
     ABI_DEALLOCATE(xcart_all)
     ABI_DEALLOCATE(fcart_all)
     ABI_DEALLOCATE(etotal_all)
     ABI_DEALLOCATE(neb_forces_all)
   end if

   ABI_DEALLOCATE(neb_forces)
   ABI_DEALLOCATE(etotal)
   ABI_DEALLOCATE(xcart)
   ABI_DEALLOCATE(fcart)
   ABI_DEALLOCATE(rprimd)
 end if ! mpi_enreg%me_cell==0

!Store acell, rprim, xred and vel for the new iteration
 call xmpi_bcast(xred,0,mpi_enreg%comm_cell,ierr)
 next_itimimage=itimimage_eff+1
 if (next_itimimage>ntimimage_stored) next_itimimage=1
 do iimage=1,nimage
   results_img(iimage,next_itimimage)%xred(:,:)    =xred(:,:,iimage)
   results_img(iimage,next_itimimage)%acell(:)     =results_img(iimage,itimimage_eff)%acell(:)
   results_img(iimage,next_itimimage)%rprim(:,:)   =results_img(iimage,itimimage_eff)%rprim(:,:)
   results_img(iimage,next_itimimage)%vel(:,:)     =results_img(iimage,itimimage_eff)%vel(:,:)
   results_img(iimage,next_itimimage)%vel_cell(:,:)=results_img(iimage,itimimage_eff)%vel_cell(:,:)
 end do
 ABI_DEALLOCATE(xred)

end subroutine predict_neb
!!***

end module m_predict_neb
!!***

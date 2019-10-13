!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_predict_string
!! NAME
!!  m_predict_string
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2019 ABINIT group (XG,ARom,MT)
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

module m_predict_string

 use defs_basis
 use m_abicore
 use m_splines
 use m_mep
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_results_img, only : results_img_type, gather_array_img, get_geometry_img

 implicit none

 private
!!***

 public :: predict_string
!!***

contains
!!***

!!****f* ABINIT/predict_string
!! NAME
!! predict_string
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images using String Method.
!! The changes on the geometry and others are predicted by rescaling the path
!! No change of acell, rprim and vel at present.
!!
!! INPUTS
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB of string method, one expect the two end images to be fixed.
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
!!      gather_array_img,get_geometry_img,mep_rk4,mep_steepest,spline,splint
!!      xmpi_bcast
!!
!! SOURCE

subroutine predict_string(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&                         ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,itimimage_eff,natom,ndynimage
 integer,intent(in) :: nimage,nimage_tot,ntimimage_stored
 type(mep_type),intent(inout) :: mep_param
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in)     :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: idynimage,ierr,ii,iimage,iatom,next_itimimage
 real(dp) :: emax,emin,step
!arrays
 real(dp),allocatable :: buffer(:,:),buffer_all(:,:)
 real(dp),allocatable :: darc(:),dimage(:),fcart(:,:,:),rprimd(:,:,:),vect(:,:),wimage(:)
 real(dp),allocatable :: x(:),y(:),z(:),x2(:),y2(:),z2(:)
 real(dp),allocatable :: xout(:),yout(:),zout(:)
 real(dp),allocatable,target :: etotal(:),xcart(:,:,:),xred(:,:,:)
 real(dp),pointer :: etotal_all(:),xcart_all(:,:,:),xred_all(:,:,:)

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

!  EVOLUTION STEP
!  ===============================================

!  Compute new atomic positions in each cell
   if      (mep_param%mep_solver==0) then ! Steepest-descent
     call mep_steepest(fcart,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)
   else if (mep_param%mep_solver==4) then ! 4th-order Runge-Kutta
     call mep_rk4(fcart,itimimage,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)
   else
     MSG_BUG("Inconsistent solver !")
   end if

!  REPARAMETRIZATION STEP
!  ===============================================

!  No reparametrization step in case of Runge-Kutta and mod(istep,4)>0
   if (mep_param%mep_solver/=4.or.mod(itimimage,4)==0) then

!    Parallelism: gather data of all images
     if (mpi_enreg%paral_img==1) then
       ABI_ALLOCATE(buffer,(6*natom+1,nimage))
       ABI_ALLOCATE(buffer_all,(6*natom+1,nimage_tot))
       ABI_ALLOCATE(xred_all,(3,natom,nimage_tot))
       ABI_ALLOCATE(xcart_all,(3,natom,nimage_tot))
       ABI_ALLOCATE(etotal_all,(nimage_tot))
       buffer=zero;ii=0
       buffer(ii+1:ii+3*natom,1:nimage)=reshape(xred ,(/3*natom,nimage/));ii=3*natom
       buffer(ii+1:ii+3*natom,1:nimage)=reshape(xcart,(/3*natom,nimage/));ii=6*natom
       buffer(ii+1           ,1:nimage)=etotal(1:nimage);ii=0
       call gather_array_img(buffer,buffer_all,mpi_enreg,allgather=.true.)
       xred_all(:,:,:) =reshape(buffer_all(ii+1:ii+3*natom,1:nimage_tot),(/3,natom,nimage_tot/));ii=3*natom
       xcart_all(:,:,:)=reshape(buffer_all(ii+1:ii+3*natom,1:nimage_tot),(/3,natom,nimage_tot/));ii=6*natom
       etotal_all(:)=buffer_all(ii+1,1:nimage_tot);ii=0
       ABI_DEALLOCATE(buffer)
       ABI_DEALLOCATE(buffer_all)
     else
       xred_all   => xred
       xcart_all  => xcart
       etotal_all => etotal
     end if

!    dimage is the distance between two images
!    darc is the parametrization on the string
!    wimage is the weight
     ABI_ALLOCATE(darc,(nimage_tot))
     ABI_ALLOCATE(dimage,(nimage_tot))
     ABI_ALLOCATE(wimage,(nimage_tot))

!    === Weights for equal arc length
     if (mep_param%string_algo/=2) then
       wimage(:)=one

!      === Weights for energy-weight arc length
     else
       emin=min(etotal_all(1),etotal_all(nimage_tot))
       emax=maxval(etotal_all(:))
       wimage(1)=one
       do iimage=2,nimage_tot
         wimage(iimage)=exp((half*(etotal_all(iimage)+etotal_all(iimage-1))-emin)/(emax-emin))
       end do
     end if

!    The distance between images is calculated
!    and normalized to a string length of 1.0
     dimage=zero
     ABI_ALLOCATE(vect,(3,natom))
     do iimage=2,nimage_tot
       dimage(iimage)=dimage(iimage-1)
!      MT april 2012: distance must be computed with cartesian coordinates
!      vect(:,:)=xred_all(:,:,iimage)-xred_all(:,:,iimage-1)
       vect(:,:)=xcart_all(:,:,iimage)-xcart_all(:,:,iimage-1)
       dimage(iimage)=dimage(iimage)+wimage(iimage)*mep_img_norm(vect)
     end do
     dimage(:)=dimage(:)/dimage(nimage_tot)
     ABI_DEALLOCATE(vect)

!    Arc lengths
     darc(1)=zero
     step=one/dble(nimage_tot-1)
     do iimage=2,nimage_tot
       darc(iimage)=darc(iimage-1)+step
     end do

!    New image coordinates are calculated and such that now the mesh is uniform
     ABI_ALLOCATE(x,(nimage_tot))
     ABI_ALLOCATE(y,(nimage_tot))
     ABI_ALLOCATE(z,(nimage_tot))
     ABI_ALLOCATE(x2,(nimage_tot))
     ABI_ALLOCATE(y2,(nimage_tot))
     ABI_ALLOCATE(z2,(nimage_tot))
     ABI_ALLOCATE(xout,(nimage_tot))
     ABI_ALLOCATE(yout,(nimage_tot))
     ABI_ALLOCATE(zout,(nimage_tot))
     do iatom=1,natom
       do iimage=1,nimage_tot
         x(iimage)=xred_all(1,iatom,iimage)
         y(iimage)=xred_all(2,iatom,iimage)
         z(iimage)=xred_all(3,iatom,iimage)
       end do
       call spline(dimage,x,nimage_tot,greatest_real,greatest_real,x2)
       call spline(dimage,y,nimage_tot,greatest_real,greatest_real,y2)
       call spline(dimage,z,nimage_tot,greatest_real,greatest_real,z2)
       call splint(nimage_tot,dimage,x,x2,nimage_tot,darc,xout)
       call splint(nimage_tot,dimage,y,y2,nimage_tot,darc,yout)
       call splint(nimage_tot,dimage,z,z2,nimage_tot,darc,zout)
!      After a spline, new image coordinate for that particular
!      atom are generated only if they are dynamical
       do idynimage=1,ndynimage
         iimage=list_dynimage(idynimage)
         ii=mpi_enreg%my_imgtab(iimage)
         xred(1,iatom,iimage)=xout(ii)
         xred(2,iatom,iimage)=yout(ii)
         xred(3,iatom,iimage)=zout(ii)
       end do
     end do  ! iatom

!    Free memory
     ABI_DEALLOCATE(x)
     ABI_DEALLOCATE(y)
     ABI_DEALLOCATE(z)
     ABI_DEALLOCATE(x2)
     ABI_DEALLOCATE(y2)
     ABI_DEALLOCATE(z2)
     ABI_DEALLOCATE(xout)
     ABI_DEALLOCATE(yout)
     ABI_DEALLOCATE(zout)
     ABI_DEALLOCATE(darc)
     ABI_DEALLOCATE(dimage)
     ABI_DEALLOCATE(wimage)
     if (mpi_enreg%paral_img==1)  then
       ABI_DEALLOCATE(xred_all)
       ABI_DEALLOCATE(xcart_all)
       ABI_DEALLOCATE(etotal_all)
     end if

   end if ! Reparametrization

!  ===============================================

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

end subroutine predict_string
!!***

end module m_predict_string
!!***

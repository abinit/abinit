!!****m* ABINIT/m_predict_pimd
!! NAME
!!  m_predict_pimd
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (GG)
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

module m_predict_pimd

 use defs_basis
 use m_abicore
 use m_pimd
 use m_xmpi
 use m_results_img

 use defs_abitypes,    only : MPI_type
 use m_geometry,       only : mkradim, mkrdim
 use m_pimd_langevin,  only : pimd_langevin_npt, pimd_langevin_nvt
 use m_pimd_nosehoover, only : pimd_nosehoover_npt, pimd_nosehoover_nvt

 implicit none

 private
!!***

 public :: predict_pimd
!!***

contains
!!***

!!****f* ABINIT/predict_pimd
!! NAME
!! predict_pimd
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics
!! Given the positions at time t and t-dtion, an estimation of the velocities at time t,
!! the forces and an estimation of the stress at time t, and an estimation of the cell at time t,
!! computes in the Path Integral Molecular Dynamics framework the new positions at time t+dtion,
!! computes self-consistently the velocities, the stress and the cell at time t and produces
!! an estimation of the velocities, stress and new cell at time t+dtion
!! No change of acell and rprim at present.
!!
!! INPUTS
!! imgmov=gives the algorithm to be used for prediction of new set of images
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! mpi_enreg=MPI-parallelisation information
!! natom=dimension of vel_timimage and xred_timimage
!! nimage=number of images (treated by current proc)
!! nimage_tot=total number of images
!! ntimimage_stored=number of time steps stored in the history
!! results_gs_timimage(ntimimage,nimage)=datastructure that hold all the history of previous computations.
!! pimd_param=datastructure that contains all the parameters necessary to Path-Integral MD
!! prtvolimg=printing volume
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
!!      m_gstateimg
!!
!! CHILDREN
!!      gather_array_img,mkradim,mkrdim,pimd_langevin_npt,pimd_langevin_nvt
!!      pimd_nosehoover_npt,pimd_nosehoover_nvt,scatter_array_img,xmpi_bcast
!!      xmpi_gather
!!
!! SOURCE

subroutine predict_pimd(imgmov,itimimage,itimimage_eff,mpi_enreg,natom,nimage,nimage_tot,&
&                       ntimimage_stored,pimd_param,prtvolimg,results_img)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imgmov,itimimage,itimimage_eff,natom
 integer,intent(in) :: nimage,nimage_tot,ntimimage_stored,prtvolimg
 type(MPI_type),intent(in) :: mpi_enreg
 type(pimd_type),intent(inout) :: pimd_param
!arrays
 type(results_img_type) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer :: ierr,ii,itime,itime_next,itime_prev
 real(dp) :: volume
!arrays
 real(dp) :: rprimd(3,3),rprimd_next(3,3),rprimd_prev(3,3),vel_cell(3,3)
 real(dp),allocatable :: mpibuf(:),mpibuffer(:,:,:),mpibuffer_all(:,:,:)
 real(dp),allocatable :: etotal(:),forces(:,:,:),stressin(:,:,:),vel(:,:,:)
 real(dp),allocatable :: xred(:,:,:),xred_next(:,:,:),xred_prev(:,:,:)

! *************************************************************************

!############# Parallelism stuff 1 #######################

!Parallelism over image: only one process per image of the cell
 if (mpi_enreg%me_cell==0) then

   itime=itimimage_eff
   itime_prev=itime-1;if (itime_prev<1) itime_prev=ntimimage_stored

   if (mpi_enreg%paral_img==0.or.mpi_enreg%me_img==0) then
     ABI_ALLOCATE(xred,(3,natom,nimage_tot))
     ABI_ALLOCATE(xred_prev,(3,natom,nimage_tot))
     ABI_ALLOCATE(xred_next,(3,natom,nimage_tot))
     ABI_ALLOCATE(etotal,(nimage_tot))
     ABI_ALLOCATE(forces,(3,natom,nimage_tot))
     ABI_ALLOCATE(stressin,(3,3,nimage_tot))
     ABI_ALLOCATE(vel,(3,natom,nimage_tot))
   end if

!  Parallelism: Gather positions/forces/velocities/stresses/energy from all images
   if (mpi_enreg%paral_img==1) then
     ABI_ALLOCATE(mpibuffer,(12,natom+1,nimage))
     do ii=1,nimage
       mpibuffer(1:3  ,1:natom,ii)=results_img(ii,itime)%xred(1:3,1:natom)
       mpibuffer(4:6  ,1:natom,ii)=results_img(ii,itime_prev)%xred(1:3,1:natom)
       mpibuffer(7:9  ,1:natom,ii)=results_img(ii,itime)%results_gs%fcart(1:3,1:natom)
       mpibuffer(10:12,1:natom,ii)=results_img(ii,itime)%vel(1:3,1:natom)
       mpibuffer(1:6  ,natom+1,ii)=results_img(ii,itime)%results_gs%strten(1:6)
       mpibuffer(7:12 ,natom+1,ii)=zero
     end do
     if (mpi_enreg%me_img==0)  then
       ABI_ALLOCATE(mpibuffer_all,(12,natom+1,nimage_tot))
     end if
     call gather_array_img(mpibuffer,mpibuffer_all,mpi_enreg,only_one_per_img=.true.,allgather=.false.)
     ABI_DEALLOCATE(mpibuffer)
     if (mpi_enreg%me_img==0) then
       do ii=1,nimage_tot
         xred     (1:3,1:natom,ii)=mpibuffer_all(1:3  ,1:natom,ii)
         xred_prev(1:3,1:natom,ii)=mpibuffer_all(4:6  ,1:natom,ii)
         forces   (1:3,1:natom,ii)=mpibuffer_all(7:9  ,1:natom,ii)
         vel      (1:3,1:natom,ii)=mpibuffer_all(10:12,1:natom,ii)
         stressin (1,1,ii)        =mpibuffer_all(1,natom+1,ii)
         stressin (2,2,ii)        =mpibuffer_all(2,natom+1,ii)
         stressin (3,3,ii)        =mpibuffer_all(3,natom+1,ii)
         stressin (3,2,ii)        =mpibuffer_all(4,natom+1,ii)
         stressin (3,1,ii)        =mpibuffer_all(5,natom+1,ii)
         stressin (2,1,ii)        =mpibuffer_all(6,natom+1,ii)
         stressin (2,3,ii)=stressin (3,2,ii)
         stressin (1,3,ii)=stressin (3,1,ii)
         stressin (1,2,ii)=stressin (2,1,ii)
       end do
       ABI_DEALLOCATE(mpibuffer_all)
     end if
     ABI_ALLOCATE(mpibuf,(nimage))
     if (mpi_enreg%me_img/=0) then
       ABI_ALLOCATE(etotal,(0))
     end if
     do ii=1,nimage
       mpibuf(ii)=results_img(ii,itime)%results_gs%etotal
     end do
     call xmpi_gather(mpibuf,nimage,etotal,nimage,0,mpi_enreg%comm_img,ierr)
     ABI_DEALLOCATE(mpibuf)
     if (mpi_enreg%me_img/=0) then
       ABI_DEALLOCATE(etotal)
     end if

!    No parallelism: simply copy positions/forces/velocities/stresses/energy
   else
     do ii=1,nimage
       xred     (:,:,ii)=results_img(ii,itime)%xred(:,:)
       xred_prev(:,:,ii)=results_img(ii,itime_prev)%xred(:,:)
       forces   (:,:,ii)=results_img(ii,itime)%results_gs%fcart(:,:)
       vel      (:,:,ii)=results_img(ii,itime)%vel(:,:)
       etotal   (    ii)=results_img(ii,itime)%results_gs%etotal
       stressin (1,1,ii)=results_img(ii,itime)%results_gs%strten(1)
       stressin (2,2,ii)=results_img(ii,itime)%results_gs%strten(2)
       stressin (3,3,ii)=results_img(ii,itime)%results_gs%strten(3)
       stressin (3,2,ii)=results_img(ii,itime)%results_gs%strten(4)
       stressin (3,1,ii)=results_img(ii,itime)%results_gs%strten(5)
       stressin (2,1,ii)=results_img(ii,itime)%results_gs%strten(6)
       stressin (2,3,ii)=stressin (3,2,ii)
       stressin (1,3,ii)=stressin (3,1,ii)
       stressin (1,2,ii)=stressin (2,1,ii)
     end do
   end if

!  Parallelism over image: only one process does the job
   if (mpi_enreg%paral_img==0.or.mpi_enreg%me_img==0) then

!    ############# PIMD MD algorithm #########################

!    Some useful quantities about the cells (common to all images)
!    Take acell and rprim from 1st image
     call mkrdim(results_img(1,itime)%acell,results_img(1,itime)%rprim,rprimd)
     call mkrdim(results_img(1,itime_prev)%acell,results_img(1,itime_prev)%rprim,rprimd_prev)
     rprimd_next(:,:)=rprimd_prev(:,:)
     vel_cell(:,:)=results_img(1,1)%vel_cell(:,:)

!    Compute the volume of the supercell
     volume=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
&     rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
&     rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))
     volume=abs(volume)

     select case(imgmov)

     case(9,10)  !Langevin

       select case(pimd_param%optcell)
       case(0)  !NVT
         call pimd_langevin_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&         rprimd,stressin,nimage_tot,vel,volume,xred,xred_next,xred_prev)
       case(2)  !NPT
         call pimd_langevin_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&         rprimd,rprimd_next,rprimd_prev,stressin,nimage_tot,vel,vel_cell,&
&         volume,xred,xred_next,xred_prev)
       end select

     case(13)  !Nose Hoover chains

       select case(pimd_param%optcell)
       case(0)  !NVT
         call pimd_nosehoover_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&         rprimd,stressin,nimage_tot,vel,volume,xred,xred_next,xred_prev)
       case(2)  !NPT
         call pimd_nosehoover_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&         rprimd,rprimd_next,rprimd_prev,stressin,nimage_tot,vel,vel_cell,&
&         volume,xred,xred_next,xred_prev)
       end select

     end select

!    ############# Parallelism stuff 2 ########################

   end if ! mpi_enreg%me_img==0

!  Parallelism: dispatch results
!  The trick: use (9,natom) to store xred,xred_next,vel for all atoms
!  use (9      ) to store rprimd_next
   ABI_ALLOCATE(mpibuffer,(9,natom+2,nimage))
   if (mpi_enreg%paral_img==1) then
     if (mpi_enreg%me_img==0) then
       ABI_ALLOCATE(mpibuffer_all,(9,natom+2,nimage_tot))
       do ii=1,nimage_tot
         mpibuffer_all(1:3,1:natom,ii)=xred_next(1:3,1:natom,ii)
         mpibuffer_all(4:6,1:natom,ii)=xred(1:3,1:natom,ii)
         mpibuffer_all(7:9,1:natom,ii)=vel(1:3,1:natom,ii)
         mpibuffer_all(1:3,natom+1,ii)=rprimd_next(1:3,1)
         mpibuffer_all(4:6,natom+1,ii)=rprimd_next(1:3,2)
         mpibuffer_all(7:9,natom+1,ii)=rprimd_next(1:3,3)
         mpibuffer_all(1:3,natom+2,ii)=vel_cell(1:3,1)
         mpibuffer_all(4:6,natom+2,ii)=vel_cell(1:3,2)
         mpibuffer_all(7:9,natom+2,ii)=vel_cell(1:3,3)
       end do
     end if
     call scatter_array_img(mpibuffer,mpibuffer_all,mpi_enreg)
     if (mpi_enreg%me_img==0)  then
       ABI_DEALLOCATE(mpibuffer_all)
     end if
   else
     do ii=1,nimage
       mpibuffer(1:3,1:natom,ii)=xred_next(1:3,1:natom,ii)
       mpibuffer(4:6,1:natom,ii)=xred(1:3,1:natom,ii)
       mpibuffer(7:9,1:natom,ii)=vel(1:3,1:natom,ii)
       mpibuffer(1:3,natom+1,ii)=rprimd_next(1:3,1)
       mpibuffer(4:6,natom+1,ii)=rprimd_next(1:3,2)
       mpibuffer(7:9,natom+1,ii)=rprimd_next(1:3,3)
       mpibuffer(1:3,natom+2,ii)=vel_cell(1:3,1)
       mpibuffer(4:6,natom+2,ii)=vel_cell(1:3,2)
       mpibuffer(7:9,natom+2,ii)=vel_cell(1:3,3)
     end do
   end if

   if (mpi_enreg%paral_img==0.or.mpi_enreg%me_img==0) then
     ABI_DEALLOCATE(xred)
     ABI_DEALLOCATE(xred_prev)
     ABI_DEALLOCATE(xred_next)
     ABI_DEALLOCATE(etotal)
     ABI_DEALLOCATE(forces)
     ABI_DEALLOCATE(stressin)
     ABI_DEALLOCATE(vel)
   end if

 else
   ABI_ALLOCATE(mpibuffer,(9,natom+2,nimage))

 end if ! mpi_enreg%me_cell==0

!Send results to all procs treating the same image
 call xmpi_bcast(mpibuffer,0,mpi_enreg%comm_cell,ierr)

!Store results in final place
 itime=itimimage_eff
 itime_prev=itime-1;if (itime_prev<1) itime_prev=ntimimage_stored
 itime_next=itime+1;if (itime_next>ntimimage_stored) itime_next=1
 do ii=1,nimage
   results_img(ii,itime_next)%xred(1:3,1:natom)=mpibuffer(1:3,1:natom,ii)
   results_img(ii,itime)%xred(1:3,1:natom)=mpibuffer(4:6,1:natom,ii)
   results_img(ii,itime_next)%vel(1:3,1:natom)=mpibuffer(7:9,1:natom,ii)
 end do
 if (pimd_param%optcell/=0) then
   do ii=1,nimage
     results_img(ii,itime)%acell(1:3)=results_img(ii,itime_prev)%acell(1:3)
     results_img(ii,itime)%rprim(1:3,1:3)=results_img(ii,itime_prev)%rprim(1:3,1:3)
     rprimd(1:3,1)=mpibuffer(1:3,natom+1,ii)
     rprimd(1:3,2)=mpibuffer(4:6,natom+1,ii)
     rprimd(1:3,3)=mpibuffer(7:9,natom+1,ii)
     results_img(ii,itime_next)%vel_cell(1:3,1)=mpibuffer(1:3,natom+2,ii)
     results_img(ii,itime_next)%vel_cell(1:3,2)=mpibuffer(4:6,natom+2,ii)
     results_img(ii,itime_next)%vel_cell(1:3,3)=mpibuffer(7:9,natom+2,ii)
     call mkradim(results_img(ii,itime_next)%acell,results_img(ii,itime_next)%rprim,rprimd)
   end do
 end if
 ABI_DEALLOCATE(mpibuffer)

end subroutine predict_pimd
!!***

end module m_predict_pimd
!!***

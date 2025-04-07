!!****m* ABINIT/m_predict_neb
!! NAME
!!  m_predict_neb
!!
!! FUNCTION
!! This module implement the Nudged Elastic Band method (several variants)
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2025 ABINIT group (MT,QD)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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
 use m_results_img, only : results_img_type, gather_array_img, scatter_array_img, &
&                          get_geometry_img
 use m_geometry, only : mkradim, metric, fcart2gred
 use m_matrix, only : matr3inv

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
!! SOURCE

subroutine predict_neb(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&                      ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

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
 integer :: ierr,ii,jj,iimage,iimage_min,iimage_max,natom_eff,next_itimimage
 logical :: test_minus_one,test_plus_one,use_reduced_coord
 real(dp) :: dvmax,dvmin,ecur,emax,emin,eref,f_para1,f_para2
 character(len=100) :: msg
!arrays
 integer,parameter :: voigt(3,3)=reshape([1,6,3,6,2,4,3,4,3],[3,3])
 real(dp),parameter :: identity_real(3,3)=reshape([one,zero,zero,zero,one,zero,zero,zero,one],[3,3])
 real(dp) :: mat3_1(3,3),mat3_2(3,3)
 real(dp),allocatable :: acell(:,:),buffer(:,:),buffer_all(:,:)
 real(dp),allocatable :: coordif(:,:,:),dimage(:),spring(:)
 real(dp),allocatable :: rprim(:,:,:),rprimd_start(:,:,:),rprimd_start_inv(:,:,:)
 real(dp),allocatable :: fcart(:,:,:),xcart(:,:,:),xred(:,:,:),strainfact_jj(:)
 real(dp),allocatable :: strten(:,:),strten_mat(:,:,:),rmet(:,:,:),ucvol(:),pressure(:)
 real(dp),allocatable :: tangent(:,:,:),vect(:,:)
 real(dp),allocatable,target :: etotal(:),neb_forces(:,:,:),rprimd(:,:,:)
 real(dp),allocatable,target :: fcart_eff(:,:,:),xcart_eff(:,:,:),xred_eff(:,:,:)
 real(dp),pointer :: coord_(:,:,:),etotal_all(:),fcart_eff_all(:,:,:), neb_forces_all(:,:,:)
 real(dp),pointer :: rprimd_all(:,:,:),xcart_eff_all(:,:,:),xred_eff_all(:,:,:)

! *************************************************************************

!Check options
 if (mep_param%neb_cell_algo/=NEB_CELL_ALGO_NONE.and.&
&    mep_param%mep_solver/=MEP_SOLVER_STEEPEST) then
   msg='Variable cell NEB only allowed with steepest descent algo!'
   ABI_ERROR(msg)
 end if

!In  case of variable-cell, 3 additional "fake" atoms are added to unit cell vectors
 natom_eff=natom ; if (mep_param%neb_cell_algo/=NEB_CELL_ALGO_NONE) natom_eff=natom_eff+3

!VC-NEB uses reduced coordinates
 use_reduced_coord=.false.
 if (mep_param%neb_cell_algo==NEB_CELL_ALGO_VCNEB) use_reduced_coord=.true.

 ABI_MALLOC(acell,(3,nimage))
 ABI_MALLOC(rprim,(3,3,nimage))
 ABI_MALLOC(xred_eff,(3,natom_eff,nimage))

!Parallelism over images: only one process per image of the cell
 if (mpi_enreg%me_cell==0) then

!  Retrieve positions and forces
   ABI_MALLOC(etotal,(nimage))
   ABI_MALLOC(xcart_eff,(3,natom_eff,nimage))
   ABI_MALLOC(fcart_eff,(3,natom_eff,nimage))
   ABI_MALLOC(rprimd,(3,3,nimage))
   ABI_MALLOC(strainfact_jj,(nimage))
   strainfact_jj(:)=one
   ABI_MALLOC(rprimd_start,(3,3,nimage))
   do iimage=1,nimage
     ii=mpi_enreg%my_imgtab(iimage)
     rprimd_start(:,:,iimage)=mep_param%rprimd_start(:,:,ii)
   end do
   ABI_MALLOC(xred,(3,natom,nimage))
   ABI_MALLOC(xcart,(3,natom,nimage))
   ABI_MALLOC(fcart,(3,natom,nimage))
   ABI_MALLOC(strten,(6,nimage))
   call get_geometry_img(results_img(:,itimimage_eff),etotal,natom,nimage,fcart,rprimd,strten,xcart,xred)
   xred_eff(1:3,1:natom,1:nimage)=xred(1:3,1:natom,1:nimage)
   xcart_eff(1:3,1:natom,1:nimage)=xcart(1:3,1:natom,1:nimage)
   fcart_eff(1:3,1:natom,1:nimage)=fcart(1:3,1:natom,1:nimage)

!  Retrieve unit cell vectors and derivatives (forces) in case of variable-cell NEB
!   => Retrieve pressure, volume, rmet, strten, ...
   if (mep_param%neb_cell_algo/=NEB_CELL_ALGO_NONE) then
     ABI_MALLOC(rmet,(3,3,nimage))
     ABI_MALLOC(ucvol,(nimage))
     ABI_MALLOC(pressure,(nimage))
     ABI_MALLOC(strten_mat,(3,3,nimage))
     ABI_MALLOC(rprimd_start_inv,(3,3,nimage))
     pressure(1:nimage)=-(strten(1,1:nimage)+strten(2,1:nimage)+strten(3,1:nimage))*third
     do iimage=1,nimage
       xcart_eff(1:3,natom+1:natom+3,iimage)=zero
       xred_eff(1:3,natom+1:natom+3,iimage)=zero
       call metric(mat3_1,mat3_2,-1,rmet(:,:,iimage),rprimd(:,:,iimage),ucvol(iimage))
       do jj=1,3; do ii=1,3
         strten_mat(ii,jj,iimage)=strten(voigt(ii,jj),iimage)
       end do; end do
       call matr3inv(rprimd_start(:,:,iimage),mat3_1(:,:))
       rprimd_start_inv(:,:,iimage) = transpose(mat3_1(:,:))
     end do

!    ==== Generalized Solid-State NEB (GSS-NEB)
!         See: Sheppard, Xiao, Chemelewski, Johnson, Henkelman, J. Chem. Phys. 136, 074103 (2012)
     if (mep_param%neb_cell_algo==NEB_CELL_ALGO_GSSNEB) then
       strainfact_jj(1:nimage)=(ucvol(1:nimage)**third)*(natom**sixth)
       do iimage=1,nimage
         mat3_1(1:3,1:3)=rprimd(1:3,1:3,iimage)-rprimd_start(1:3,1:3,iimage)
         xcart_eff(1:3,natom+1:natom+3,iimage)=strainfact_jj(iimage)*matmul(rprimd_start_inv(:,:,iimage),mat3_1)       
         strten_mat(1:3,1:3,iimage)=strten_mat(1:3,1:3,iimage)+pressure(iimage)*identity_real(1:3,1:3)
         fcart_eff(1:3,natom+1:natom+3,iimage)= &
&                 -(ucvol(iimage)/strainfact_jj(iimage))*strten_mat(1:3,1:3,iimage)
       end do
     end if

!    ==== Variable-cell NEB (VC-NEB, in reduced coordinates) - At present, doesnt work
!         See: Qian, Dong, Zhou, Tian, Oganov, Wang, Comp. Phys. Comm. 184, 2111 (2013)
     if (mep_param%neb_cell_algo==NEB_CELL_ALGO_VCNEB) then
       do iimage=1,nimage
         xred_eff(1:3,natom+1:natom+3,iimage)= &
&           matmul(rprimd(1:3,1:3,iimage),rprimd_start_inv(1:3,1:3,iimage))-identity_real(1:3,1:3)
         call fcart2gred(fcart(1:3,1:natom,iimage),fcart_eff(1:3,1:natom,iimage),rprimd(1:3,1:3,iimage),natom)
         fcart_eff(1:3,1:natom,iimage)=matmul(rmet(1:3,1:3,iimage),fcart_eff(1:3,1:natom,iimage))
         mat3_1(1:3,1:3)=identity_real(1:3,1:3)+transpose(xred_eff(1:3,natom+1:natom+3,iimage))
         call matr3inv(mat3_1,mat3_2)
         mat3_1=transpose(mat3_2)
         mat3_2(1:3,1:3)=(strten_mat(1:3,1:3,iimage)+pressure(iimage))*identity_real(1:3,1:3)*ucvol(iimage)
         fcart_eff(1:3,natom+1:natom+3,iimage)=-matmul(mat3_2,mat3_1)
       end do
     end if

     ABI_FREE(rmet)
     ABI_FREE(ucvol)
     ABI_FREE(pressure)
     ABI_FREE(strten_mat)
     ABI_FREE(rprimd_start_inv)
   end if

   ABI_FREE(xred)
   ABI_FREE(xcart)
   ABI_FREE(fcart)
   ABI_FREE(strten)

!  Array containing effective NEB forces (F_ortho+F_spring)
   ABI_MALLOC(neb_forces,(3,natom_eff,nimage))

!  Parallelism: gather data of all images
   if (mpi_enreg%paral_img==1) then
     ABI_MALLOC(buffer,(9*natom_eff+10,nimage))
     ABI_MALLOC(buffer_all,(9*natom_eff+10,nimage_tot))
     ABI_MALLOC(etotal_all,(nimage_tot))
     ABI_MALLOC(xcart_eff_all,(3,natom_eff,nimage_tot))
     ABI_MALLOC(fcart_eff_all,(3,natom_eff,nimage_tot))
     ABI_MALLOC(xred_eff_all,(3,natom_eff,nimage_tot))
     ABI_MALLOC(rprimd_all,(3,3,nimage_tot))
     ABI_MALLOC(neb_forces_all,(3,natom_eff,nimage_tot))
     buffer=zero;ii=0
     buffer(ii+1               ,1:nimage)=etotal(1:nimage);ii=ii+1
     buffer(ii+1:ii+9          ,1:nimage)=reshape(rprimd,(/9,nimage/));ii=ii+9
     buffer(ii+1:ii+3*natom_eff,1:nimage)=reshape(xcart_eff,(/3*natom_eff,nimage/));ii=ii+3*natom_eff
     buffer(ii+1:ii+3*natom_eff,1:nimage)=reshape(xred_eff,(/3*natom_eff,nimage/));ii=ii+3*natom_eff
     buffer(ii+1:ii+3*natom_eff,1:nimage)=reshape(fcart_eff,(/3*natom_eff,nimage/));ii=ii+3*natom_eff
     call gather_array_img(buffer,buffer_all,mpi_enreg,allgather=.true.)
     ii=0
     etotal_all(:)=buffer_all(ii+1,1:nimage_tot);ii=ii+1
     rprimd_all(:,:,:)=reshape(buffer_all(ii+1:ii+9,1:nimage_tot),(/3,3,nimage_tot/));ii=ii+9
     xcart_eff_all(:,:,:)=reshape(buffer_all(ii+1:ii+3*natom_eff,1:nimage_tot),(/3,natom_eff,nimage_tot/));ii=ii+3*natom_eff
     xred_eff_all(:,:,:)=reshape(buffer_all(ii+1:ii+3*natom_eff,1:nimage_tot),(/3,natom_eff,nimage_tot/));ii=ii+3*natom_eff
     fcart_eff_all(:,:,:)=reshape(buffer_all(ii+1:ii+3*natom_eff,1:nimage_tot),(/3,natom_eff,nimage_tot/));ii=ii+3*natom_eff
     ABI_FREE(buffer)
     ABI_FREE(buffer_all)
   else
     etotal_all        => etotal
     rprimd_all        => rprimd
     xcart_eff_all     => xcart_eff
     xred_eff_all      => xred_eff
     fcart_eff_all     => fcart_eff
     neb_forces_all    => neb_forces
   end if

!  coordif is the vector between two images
!  dimage is the distance between two images
!  tangent is the tangent at image i
   ABI_MALLOC(coordif,(3,natom_eff,nimage_tot))
   ABI_MALLOC(tangent,(3,natom_eff,nimage_tot))
   ABI_MALLOC(dimage,(nimage_tot))

!  Compute distances between images
   coordif(:,:,1)=zero;dimage(1)=zero
   coord_ => xcart_eff_all ; if (use_reduced_coord) coord_ => xred_eff_all
   do iimage=2,nimage_tot
     coordif(:,:,iimage)=coord_(:,:,iimage)-coord_(:,:,iimage-1)
   end do

   do iimage=2,nimage_tot
     dimage(iimage)=mep_img_norm(coordif(:,:,iimage))
   end do

!  Compute tangent (not normalized)
   tangent(:,:,1)=zero
   tangent(:,:,nimage_tot)=zero
!  === Original definition of tangent
   if (mep_param%neb_algo==NEB_ALGO_STANDARD) then
     do iimage=2,nimage_tot-1
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
   ABI_MALLOC(spring,(nimage_tot))
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
   if (mep_param%neb_algo==NEB_ALGO_CINEB.and.itimimage>=mep_param%cineb_start) then
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
   neb_forces_all(:,:,1)=fcart_eff_all(:,:,1)
   neb_forces_all(:,:,nimage_tot)=fcart_eff_all(:,:,nimage_tot)
   do iimage=2,nimage_tot-1
!    === Standard NEB
     if (iimage/=iimage_max) then
!    if (iimage/=iimage_min.and.iimage/=iimage_max) then
       f_para1=mep_img_dotp(fcart_eff_all(:,:,iimage),tangent(:,:,iimage))
       if (mep_param%neb_algo==NEB_ALGO_STANDARD) then   ! Original NEB algo
         ABI_MALLOC(vect,(3,natom_eff))
         vect(:,:)=spring(iimage+1)*coordif(:,:,iimage+1)-spring(iimage)*coordif(:,:,iimage)
         f_para2=mep_img_dotp(vect(:,:),tangent(:,:,iimage))
         ABI_FREE(vect)
       else       ! Modification from J. Chem. Phys. 113, 9978 (2000) [[cite:Henkelman2000]]
         f_para2=spring(iimage+1)*dimage(iimage+1)-spring(iimage)*dimage(iimage)
       end if
       neb_forces_all(:,:,iimage)=fcart_eff_all(:,:,iimage) &       ! F_ortho
&      -f_para1*tangent(:,:,iimage) &  ! F_ortho
&      +f_para2*tangent(:,:,iimage)    ! F_spring
!      === CI-NEB for the image(s) with maximal energy
     else
       f_para1=mep_img_dotp(fcart_eff_all(:,:,iimage),tangent(:,:,iimage))
       neb_forces_all(:,:,iimage)=fcart_eff_all(:,:,iimage)-two*f_para1*tangent(:,:,iimage)
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
   if (mep_param%mep_solver==MEP_SOLVER_STEEPEST) then ! Steepest-descent
     call mep_steepest(neb_forces,list_dynimage,mep_param,natom,natom_eff,&
&                      ndynimage,nimage,rprimd,xcart_eff,xred_eff,&
&                      strainfact=strainfact_jj,rprimd_start=rprimd_start,&
&                      use_reduced_coord=use_reduced_coord)
   else if (mep_param%mep_solver==MEP_SOLVER_QUICKMIN) then ! Quick-min
     call mep_qmin(neb_forces,itimimage,list_dynimage,mep_param,natom,ndynimage, &
&     nimage,rprimd,xcart_eff,xred_eff)
   else if (mep_param%mep_solver==MEP_SOLVER_LBFGS) then ! Local BFGS
     call mep_lbfgs(neb_forces,itimimage,list_dynimage,mep_param,natom,ndynimage,&
&     nimage,rprimd,xcart_eff,xred_eff)
   else if (mep_param%mep_solver==MEP_SOLVER_GBFGS) then ! Global BFGS
     call mep_gbfgs(neb_forces,itimimage,list_dynimage,mep_param,mpi_enreg,natom,ndynimage,&
&     nimage,nimage_tot,rprimd,xcart_eff,xred_eff)
   else
     ABI_BUG("Inconsistent solver !")
   end if

!  Compute new acell and rprim from new rprimd
   !if (mep_param%neb_cell_algo/=NEB_CELL_ALGO_NONE) then
     do iimage=1,nimage
       call mkradim(acell(:,iimage),rprim(:,:,iimage),rprimd(:,:,iimage))
     end do
   !end if

!  Free memory
   ABI_FREE(spring)
   ABI_FREE(coordif)
   ABI_FREE(tangent)
   ABI_FREE(dimage)
   if (mpi_enreg%paral_img==1)  then
     ABI_FREE(etotal_all)
     ABI_FREE(rprimd_all)
     ABI_FREE(xred_eff_all)
     ABI_FREE(xcart_eff_all)
     ABI_FREE(fcart_eff_all)
     ABI_FREE(neb_forces_all)
   end if

   ABI_FREE(neb_forces)
   ABI_FREE(etotal)
   ABI_FREE(rprimd)
   ABI_FREE(xcart_eff)
   ABI_FREE(fcart_eff)
   ABI_FREE(strainfact_jj)
   ABI_FREE(rprimd_start)

 end if ! mpi_enreg%me_cell==0

!Store acell, rprim, xred and vel for the new iteration
 call xmpi_bcast(xred_eff,0,mpi_enreg%comm_cell,ierr)
 call xmpi_bcast(acell,0,mpi_enreg%comm_cell,ierr)
 call xmpi_bcast(rprim,0,mpi_enreg%comm_cell,ierr)
 next_itimimage=itimimage_eff+1
 if (next_itimimage>ntimimage_stored) next_itimimage=1
 do iimage=1,nimage
   results_img(iimage,next_itimimage)%xred(:,1:natom)=xred_eff(:,1:natom,iimage)
   if (mep_param%neb_cell_algo==NEB_CELL_ALGO_NONE) then
     results_img(iimage,next_itimimage)%acell(:)   =results_img(iimage,itimimage_eff)%acell(:)
     results_img(iimage,next_itimimage)%rprim(:,:) =results_img(iimage,itimimage_eff)%rprim(:,:)
   else
     results_img(iimage,next_itimimage)%acell(:)   =acell(:,iimage)
     results_img(iimage,next_itimimage)%rprim(:,:) =rprim(:,:,iimage)
   end if
   results_img(iimage,next_itimimage)%vel(:,:)     =results_img(iimage,itimimage_eff)%vel(:,:)
   results_img(iimage,next_itimimage)%vel_cell(:,:)=results_img(iimage,itimimage_eff)%vel_cell(:,:)
 end do

 ABI_FREE(xred_eff)
 ABI_FREE(acell)
 ABI_FREE(rprim)

end subroutine predict_neb
!!***

end module m_predict_neb
!!***

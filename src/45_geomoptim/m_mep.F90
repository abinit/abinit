!!****m* ABINIT/m_mep
!! NAME
!!  m_mep
!!
!! FUNCTION
!!  This module provides several routines and datatypes for the
!!  Minimal Energy Path (MEP) search implementation.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2020 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_mep

 use defs_basis
 use m_abicore
 use m_errors
 use m_dtset
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_geometry,    only : fred2fcart, fcart2fred, xcart2xred, xred2xcart, metric
 use m_bfgs,        only : hessupdt
 use m_results_img, only : results_img_type, gather_array_img

 implicit none

 private

!public procedures
 public :: mep_init
 public :: mep_destroy
 public :: mep_steepest
 public :: mep_qmin
 public :: mep_lbfgs
 public :: mep_gbfgs
 public :: mep_rk4
 public :: mep_img_dotp
 public :: mep_img_norm
 public :: mep_img_dotp_red
 public :: mep_img_norm_red
!!***

!!****t* m_mep/mep_type
!! NAME
!! mep_type
!!
!! FUNCTION
!! Datatype with the variables required to perform MEP search
!!
!! SOURCE

 type,public :: mep_type
! Scalars
  integer  :: cineb_start  ! Starting iteration for the CI-NEB
  integer  :: mep_solver   ! Selection of a solver for the ODE
  integer  :: neb_algo     ! Selection of the variant of the NEB method
  integer  :: string_algo  ! Selection of the variant of the String Method
  real(dp) :: fxcartfactor ! Time step for steepest descent or RK4
  real(dp) :: mep_mxstep   ! Selection of a max. step size for the ODE
! Arrays
  integer,pointer      :: iatfix(:,:)=>null() ! Atoms to fix (this pointer is associated with dtset%iatfix)
  real(dp)             :: neb_spring(2)       ! Spring constants for the NEB method
  real(dp),allocatable :: bfgs_xprev(:,:,:)   ! BFGS storage (prev positions)
  real(dp),allocatable :: bfgs_fprev(:,:,:)   ! BFGS storage (prev forces)
  real(dp),allocatable :: gbfgs_hess(:,:)     ! global-BFGS storage (Hessian matrix)
  real(dp),allocatable :: lbfgs_hess(:,:,:)   ! local-BFGS storage (Hessian matrix)
  real(dp),allocatable :: qmin_vel(:,:,:)     ! Quick-min algo storage (velocities)
  real(dp),allocatable :: rk4_xcart1(:,:,:)   ! 4th-order Runge-Kutta storage
  real(dp),allocatable :: rk4_fcart1(:,:,:)   ! 4th-order Runge-Kutta storage
  real(dp),allocatable :: rk4_fcart2(:,:,:)   ! 4th-order Runge-Kutta storage
  real(dp),allocatable :: rk4_fcart3(:,:,:)   ! 4th-order Runge-Kutta storage
 end type mep_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_mep/mep_init
!! NAME
!!  mep_init
!!
!! FUNCTION
!!  Initialize a datastructure of type mep_type.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in current dataset
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mep_param=datastructure of type mep_type.
!!            several parameters for Minimal Energy Path (MEP) search.
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_init(dtset,mep_param)

!Arguments ------------------------------------
!scalars
 type(dataset_type),target,intent(in) :: dtset
 type(mep_type),intent(inout) :: mep_param

!************************************************************************

 if((dtset%imgmov==1).or.(dtset%imgmov==2).or.(dtset%imgmov==5))then
   mep_param%cineb_start  = dtset%cineb_start
   mep_param%mep_solver   = dtset%mep_solver
   mep_param%neb_algo     = dtset%neb_algo
   mep_param%string_algo  = dtset%string_algo
   mep_param%fxcartfactor = dtset%fxcartfactor
   mep_param%mep_mxstep   = dtset%mep_mxstep
   mep_param%neb_spring   = dtset%neb_spring
   mep_param%iatfix       =>dtset%iatfix
 else
   mep_param%cineb_start  = -1
   mep_param%mep_solver   = -1
   mep_param%neb_algo     = -1
   mep_param%string_algo  = -1
   mep_param%fxcartfactor = zero
   mep_param%mep_mxstep   = 100._dp
   mep_param%neb_spring   = zero
   nullify(mep_param%iatfix)
 end if

end subroutine mep_init
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_destroy
!! NAME
!!  mep_destroy
!!
!! FUNCTION
!!  Destroy the content of a datastructure of type mep_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mep_param=datastructure of type mep_type.
!!            several parameters for Minimal Energy Path (MEP) search.
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_destroy(mep_param)

!Arguments ------------------------------------
!scalars
 type(mep_type),intent(inout) :: mep_param

!************************************************************************

 ABI_SFREE(mep_param%bfgs_xprev)
 ABI_SFREE(mep_param%gbfgs_hess)
 ABI_SFREE(mep_param%bfgs_fprev)
 ABI_SFREE(mep_param%lbfgs_hess)
 ABI_SFREE(mep_param%qmin_vel)
 ABI_SFREE(mep_param%rk4_xcart1)
 ABI_SFREE(mep_param%rk4_fcart1)
 ABI_SFREE(mep_param%rk4_fcart2)
 ABI_SFREE(mep_param%rk4_fcart3)

 nullify(mep_param%iatfix)

end subroutine mep_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_steepest
!! NAME
!!  mep_steepest
!!
!! FUNCTION
!!  Make a path (string of images) evolve according to a steepest descent algorithm
!!
!! INPUTS
!!  fcart(3,natom,nimage)=cartesian forces in each image along the path
!!  list_dynimage(nimage)=list of dynamical images.
!!  mep_param=datastructure of type mep_type.
!!            several parameters for Minimal Energy Path (MEP) search.
!!  natom=number of atoms
!!  ndynimage=number of dynamical images along the path
!!  nimage=number of images (including static ones)
!!  results_img(nimage)=datastructure that hold data for each image
!!                      (positions, forces, energy, ...)
!!  rprimd(3,3,nimage)=dimensional primitive translations for each image along the path
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xcart(3,natom,nimage)=cartesian coordinates of atoms in each image along the path
!!                        before and after time evolution
!!  xred(3,natom,nimage)=reduced coordinates of atoms in each image along the path
!!                       before and after time evolution
!!
!! PARENTS
!!      m_predict_neb,m_predict_steepest,m_predict_string
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_steepest(fcart,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndynimage,nimage
 type(mep_type),intent(in) :: mep_param
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 real(dp),intent(in) :: fcart(3,natom,nimage),rprimd(3,3,nimage)
 real(dp),intent(inout) :: xcart(3,natom,nimage),xred(3,natom,nimage)
!Local variables-------------------------------
!scalars
 integer :: iatom,idynimage,iimage
 real(dp) :: stepsize
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: xred_old(:,:),xstep(:,:)

!************************************************************************

 ABI_MALLOC(xred_old,(3,natom))
 ABI_MALLOC(xstep,(3,natom))

 do idynimage=1,ndynimage
   iimage=list_dynimage(idynimage)
   xred_old(:,:)=xred(:,:,iimage)

!  Compute image step
!  Note that one uses fcart, for which the sum of forces on all atoms vanish
   xstep(:,:)=mep_param%fxcartfactor*fcart(:,:,iimage)
   stepsize=mep_img_norm(xstep)
   if (stepsize>=mep_param%mep_mxstep) then
     xstep=xstep*mep_param%mep_mxstep/stepsize
     write(msg,'(a,i3,a)') " Restricting step size of image ",iimage,"."
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out ,msg,'COLL')
   end if

!  Update positions
   xcart(:,:,iimage)=xcart(:,:,iimage)+xstep(:,:)
   call xcart2xred(natom,rprimd(:,:,iimage),xcart(:,:,iimage),xred(:,:,iimage))

!  In case atom is fixed, we restore its previous value
   do iatom=1,natom
     if (any(mep_param%iatfix(:,iatom)==1)) then
       where(mep_param%iatfix(:,iatom)==1)
         xred(:,iatom,iimage)=xred_old(:,iatom)
       end where
       call xred2xcart(1,rprimd(:,:,iimage),xcart(:,iatom,iimage),xred(:,iatom,iimage))
     end if
   end do

 end do

 ABI_FREE(xred_old)
 ABI_FREE(xstep)

end subroutine mep_steepest
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_qmin
!! NAME
!!  mep_qmin
!!
!! FUNCTION
!!  Make a path (string of images) evolve according to a quick-minimizer algorithm
!!
!! INPUTS
!!  fcart(3,natom,nimage)=cartesian forces in each image along the path
!!  itime=time step
!!  list_dynimage(nimage)=list of dynamical images.
!!  mep_param=datastructure of type mep_type.
!!            several parameters for Minimal Energy Path (MEP) search.
!!  natom=number of atoms
!!  ndynimage=number of dynamical images along the path
!!  nimage=number of images (including static ones)
!!  results_img(nimage)=datastructure that hold data for each image
!!                      (positions, forces, energy, ...)
!!  rprimd(3,3,nimage)=dimensional primitive translations for each image along the path
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xcart(3,natom,nimage)=cartesian coordinates of atoms in each image along the path
!!                        before and after time evolution
!!  xred(3,natom,nimage)=reduced coordinates of atoms in each image along the path
!!                       before and after time evolution
!!
!! PARENTS
!!      m_predict_neb
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_qmin(fcart,itime,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,natom,ndynimage,nimage
 type(mep_type),intent(inout) :: mep_param
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 real(dp),intent(in) :: fcart(3,natom,nimage),rprimd(3,3,nimage)
 real(dp),intent(inout) :: xcart(3,natom,nimage),xred(3,natom,nimage)
!Local variables-------------------------------
!scalars
 integer :: iatom,idynimage,iimage
 real(dp) :: stepsize,vdotf
 character(len=500) :: msg
!arrays
 real(dp) :: vel_red(3)
 real(dp),allocatable :: xred_old(:,:),xstep(:,:)

!***********************************************************************

!Allocate history array (at first time step)
 if (itime==1) then
   if (allocated(mep_param%qmin_vel)) then
     ABI_FREE(mep_param%qmin_vel)
   end if
   ABI_MALLOC(mep_param%qmin_vel,(3,natom,ndynimage))
   mep_param%qmin_vel=zero
 end if

 ABI_MALLOC(xred_old,(3,natom))
 ABI_MALLOC(xstep,(3,natom))

 do idynimage=1,ndynimage
   iimage=list_dynimage(idynimage)
   xred_old(:,:)=xred(:,:,iimage)

!  Compute velocities
   vdotf=mep_img_dotp(mep_param%qmin_vel(:,:,idynimage),fcart(:,:,iimage))
   if (vdotf>=zero) then
     mep_param%qmin_vel(:,:,idynimage)=vdotf*fcart(:,:,iimage) &
&                              /mep_img_norm(fcart(:,:,iimage))
   else
     mep_param%qmin_vel(:,:,idynimage)=zero
     write(msg,'(a,i3,a)') " Setting velocities of image ",iimage," to zero."
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out ,msg,'COLL')
   end if
   mep_param%qmin_vel(:,:,idynimage)=mep_param%qmin_vel(:,:,idynimage) &
&                   +mep_param%fxcartfactor*fcart(:,:,iimage)

!  Compute image step
   xstep(:,:)=mep_param%fxcartfactor*mep_param%qmin_vel(:,:,idynimage)
   stepsize=mep_img_norm(xstep)
   if (stepsize>=mep_param%mep_mxstep) then
     xstep=xstep*mep_param%mep_mxstep/stepsize
     write(msg,'(a,i3,a)') " Restricting step size of image ",iimage,"."
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out ,msg,'COLL')
   end if

!  Update positions
   xcart(:,:,iimage)=xcart(:,:,iimage)+xstep(:,:)
   call xcart2xred(natom,rprimd(:,:,iimage),xcart(:,:,iimage),xred(:,:,iimage))

!  In case atom is fixed, we restore its previous value
   do iatom=1,natom
     if (any(mep_param%iatfix(:,iatom)==1)) then
       call xcart2xred(1,rprimd(:,:,iimage),mep_param%qmin_vel(:,iatom,idynimage),vel_red)
       where(mep_param%iatfix(:,iatom)==1)
         xred(:,iatom,iimage)=xred_old(:,iatom)
         vel_red(:)=zero
       end where
       call xred2xcart(1,rprimd(:,:,iimage),xcart(:,iatom,iimage),xred(:,iatom,iimage))
       call xred2xcart(1,rprimd(:,:,iimage),mep_param%qmin_vel(:,iatom,idynimage),vel_red)
     end if
   end do

 end do

 ABI_FREE(xred_old)
 ABI_FREE(xstep)

end subroutine mep_qmin
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_lbfgs
!! NAME
!!  mep_lbfgs
!!
!! FUNCTION
!!  Make a path (string of images) evolve according to a
!!  local Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
!!
!! INPUTS
!!  itime=time step
!!  list_dynimage(nimage)=list of dynamical images.
!!  mep_param=datastructure of type mep_type.
!!            several parameters for Minimal Energy Path (MEP) search.
!!  natom=number of atoms
!!  ndynimage=number of dynamical images along the path
!!  nimage=number of images (including static ones)
!!  rprimd(3,3,nimage)=dimensional primitive translations for each image along the path
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mep_param=datastructure of type mep_type.
!!            History for Runge-Kutta algorithm is filled up
!!  xcart(3,natom,nimage)=cartesian coordinates of atoms in each image along the path
!!                        before and after time evolution
!!  xred(3,natom,nimage)=reduced coordinates of atoms in each image along the path
!!                       before and after time evolution
!!
!! NOTES
!!  Could see Numerical Recipes (Fortran), 1986, page 307.
!!
!! PARENTS
!!      m_predict_neb
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_lbfgs(fcart,itime,list_dynimage,mep_param,natom,ndynimage,&
&                    nimage,rprimd,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,natom,ndynimage,nimage
 type(mep_type),intent(inout) :: mep_param
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 real(dp),intent(in) :: rprimd(3,3,nimage)
 real(dp),intent(in) :: fcart(3,natom,nimage)
 real(dp),intent(inout) :: xcart(3,natom,nimage),xred(3,natom,nimage)
!Local variables-------------------------------
!scalars
 integer :: iatom,idynimage,ii,iimage,indi,indj,jatom,jj
 logical :: reset
 real(dp),parameter :: initial_Hessian=1._dp ! in Bohr^2/Hartree
 real(dp) :: dot1,dot2,stepsize,ucvol
 character(len=500) :: msg
!arrays
 real(dp) :: gprimd(3,3),gmet(3,3),rmet(3,3)
 real(dp),allocatable :: fred(:,:),xstep(:,:)

!************************************************************************

!Allocate history array (at first time step)
 if (itime==1) then
   if (allocated(mep_param%bfgs_xprev)) then
     ABI_FREE(mep_param%bfgs_xprev)
   end if
   if (allocated(mep_param%bfgs_fprev)) then
     ABI_FREE(mep_param%bfgs_fprev)
   end if
   if (allocated(mep_param%lbfgs_hess)) then
     ABI_FREE(mep_param%lbfgs_hess)
   end if
   ABI_MALLOC(mep_param%bfgs_xprev,(3,natom,ndynimage))
   ABI_MALLOC(mep_param%bfgs_fprev,(3,natom,ndynimage))
   ABI_MALLOC(mep_param%lbfgs_hess,(3*natom,3*natom,ndynimage))
   mep_param%bfgs_xprev=zero
   mep_param%bfgs_fprev=zero
 end if

!Temporary storage
 ABI_MALLOC(fred,(3,natom))
 ABI_MALLOC(xstep,(3,natom))

!Loop over images
 do idynimage=1,ndynimage
   iimage=list_dynimage(idynimage)
   call metric(gmet,gprimd,-1,rmet,rprimd(:,:,iimage),ucvol)
   call fcart2fred(fcart(:,:,iimage),fred,rprimd(:,:,iimage),natom)

!  Test if a reset is needed
   reset=.false.
   if (itime>1) then
     dot1=mep_img_dotp(mep_param%bfgs_fprev(:,:,idynimage),fred)
     dot2=mep_img_dotp(mep_param%bfgs_fprev(:,:,idynimage), &
&                      mep_param%bfgs_fprev(:,:,idynimage))
!     dot1=mep_img_dotp_red(rmet,mep_param%bfgs_fprev(:,:,idynimage),fred)
!     dot2=mep_img_dotp_red(rmet,mep_param%bfgs_fprev(:,:,idynimage), &
!&                               mep_param%bfgs_fprev(:,:,idynimage))
     reset=((dot2<two*abs(dot1)).or.abs(dot2)<tol8)
     if (reset) then
       write(msg,'(a,i3,a)') " Resetting Hessian matrix for image ",iimage,"."
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out ,msg,'COLL')
     end if
   end if

!  ===> First step or reset: initialize the Hessian matrix (in reduced coordinates)
   if (itime==1.or.reset) then
     mep_param%lbfgs_hess(:,:,idynimage)=zero
     do iatom=1,natom
       indi=3*(iatom-1)
       do ii=1,3
         do jj=1,3
           if (mep_param%iatfix(ii,iatom)==0.and. &
&              mep_param%iatfix(jj,iatom)==0) then
             mep_param%lbfgs_hess(indi+ii,indi+jj,idynimage)=gmet(ii,jj)*initial_Hessian
           end if
         end do
       end do
     end do

!  ===> Other steps: update the Hessian matrix
   else
     call hessupdt(mep_param%lbfgs_hess(:,:,idynimage),&
&                  mep_param%iatfix,natom,3*natom, &
                   xred(:,:,iimage),mep_param%bfgs_xprev(:,:,idynimage),&
                   fred(:,:),mep_param%bfgs_fprev(:,:,idynimage))
   end if

!  Update history
   mep_param%bfgs_xprev(:,:,idynimage)=xred(:,:,iimage)
   mep_param%bfgs_fprev(:,:,idynimage)=fred(:,:)

!  Compute image step
   xstep=zero
   do iatom=1,natom
     indi=3*(iatom-1)
     do ii=1,3
       do jatom=1,natom
         indj=3*(jatom-1)
         do jj=1,3
           xstep(ii,iatom)=xstep(ii,iatom) &
&             -mep_param%lbfgs_hess(indi+ii,indj+jj,idynimage)*fred(jj,jatom)
         end do
       end do
     end do
   end do

!  Restrict image step size
   stepsize=mep_img_norm_red(rmet,xstep)
   if (stepsize>=mep_param%mep_mxstep) then
     xstep=xstep*mep_param%mep_mxstep/stepsize
     write(msg,'(a,i3,a)') " Restricting BFGS step size of image ",iimage,"."
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out ,msg,'COLL')
   end if

!  Update positions
   xred(:,:,iimage)=xred(:,:,iimage)+xstep(:,:)

!  In case atom is fixed, we restore its previous value
   where(mep_param%iatfix(:,:)==1)
     xred(:,:,iimage)=mep_param%bfgs_xprev(:,:,idynimage)
   end where

   call xred2xcart(natom,rprimd(:,:,iimage),xcart(:,:,iimage),xred(:,:,iimage))

!End loop over images
 end do

 ABI_FREE(fred)
 ABI_FREE(xstep)

end subroutine mep_lbfgs
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_gbfgs
!! NAME
!!  mep_gbfgs
!!
!! FUNCTION
!!  Make a path (string of images) evolve according to a
!!  global Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
!!
!! INPUTS
!!  itime=time step
!!  list_dynimage(nimage)=list of dynamical images.
!!  mep_param=datastructure of type mep_type.
!!  mpi_enreg=MPI-parallelisation information
!!  mep_param=several parameters for Minimal Energy Path (MEP) search.
!!  natom=number of atoms
!!  ndynimage=number of dynamical images along the path
!!  nimage=number of images (including static ones)
!!  nimage_tot=total number of images
!!  rprimd(3,3,nimage)=dimensional primitive translations for each image along the path
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mep_param=datastructure of type mep_type.
!!            History for Runge-Kutta algorithm is filled up
!!  xcart(3,natom,nimage)=cartesian coordinates of atoms in each image along the path
!!                        before and after time evolution
!!  xred(3,natom,nimage)=reduced coordinates of atoms in each image along the path
!!                       before and after time evolution
!!
!! NOTES
!!  Could see Numerical Recipes (Fortran), 1986, page 307.
!!  Has to work in cartesian coordinates
!!
!! PARENTS
!!      m_predict_neb
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_gbfgs(fcart,itime,list_dynimage,mep_param,mpi_enreg,natom,&
&                    ndynimage,nimage,nimage_tot,rprimd,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,natom,ndynimage,nimage,nimage_tot
 type(mep_type),intent(inout) :: mep_param
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 real(dp),intent(in) :: fcart(3,natom,nimage),rprimd(3,3,nimage)
 real(dp),intent(inout) :: xcart(3,natom,nimage),xred(3,natom,nimage)
!Local variables-------------------------------
!scalars
 integer :: iatom,idynimage,ii,iimage,iimage_tot,indi,indj,ierr
 integer :: jdynimage,jatom,jj,mu,ndynimage_tot,nu
 logical :: reset
 real(dp),parameter :: initial_Hessian=1._dp ! in Bohr^2/Hartree
 real(dp) :: dot1,dot2,stepsize,ucvol
 character(len=500) :: msg
!arrays
 integer,allocatable :: dynimage_tot(:),iatfix_fake(:,:),ind_dynimage_tot(:)
 integer,allocatable :: list_dynimage_tot(:)
 real(dp) :: favg(3),gprimd(3,3),gmet(3,3),rmet(3,3)
 real(dp),allocatable :: buffer(:,:),buffer_all(:,:),fred(:,:)
 real(dp),allocatable :: fcart_all(:,:,:),fcartp_all(:,:,:)
 real(dp),allocatable :: gmet_all(:,:,:),gprimd_all(:,:,:),rprimd_all(:,:,:)
 real(dp),allocatable :: xcart_all(:,:,:),xcartp_all(:,:,:),xred_old(:,:),xstep_all(:,:,:)

!************************************************************************

!Retrieve indexes of all dynamical images
 ABI_MALLOC(ind_dynimage_tot,(nimage_tot))
 if (mpi_enreg%paral_img==1) then
   ABI_MALLOC(dynimage_tot,(nimage_tot))
   dynimage_tot=0
   do idynimage=1,ndynimage
     iimage=list_dynimage(idynimage)
     iimage_tot=mpi_enreg%my_imgtab(iimage)
     dynimage_tot(iimage_tot)=1
   end do
   call xmpi_sum(dynimage_tot,mpi_enreg%comm_img,ierr)
   ndynimage_tot=count(dynimage_tot(:)>0)
   ABI_MALLOC(list_dynimage_tot,(ndynimage_tot))
   idynimage=0;ind_dynimage_tot(:)=-1
   do iimage_tot=1,nimage_tot
     if (dynimage_tot(iimage_tot)>0) then
       idynimage=idynimage+1
       ind_dynimage_tot(iimage_tot)=idynimage
       list_dynimage_tot(idynimage)=iimage_tot
     end if
   end do
   ABI_FREE(dynimage_tot)
 else
   ndynimage_tot=ndynimage
   ABI_MALLOC(list_dynimage_tot,(ndynimage))
   ind_dynimage_tot(:)=-1
   do idynimage=1,ndynimage
     ind_dynimage_tot(list_dynimage(idynimage))=idynimage
     list_dynimage_tot(idynimage)=list_dynimage(idynimage)
   end do
 end if

!Allocate history array (at first time step)
 if (itime==1) then
   if (allocated(mep_param%bfgs_xprev)) then
     ABI_FREE(mep_param%bfgs_xprev)
   end if
   if (allocated(mep_param%bfgs_fprev)) then
     ABI_FREE(mep_param%bfgs_fprev)
   end if
   if (allocated(mep_param%gbfgs_hess)) then
     ABI_FREE(mep_param%gbfgs_hess)
   end if
   ABI_MALLOC(mep_param%bfgs_xprev,(3,natom,ndynimage))
   ABI_MALLOC(mep_param%bfgs_fprev,(3,natom,ndynimage))
   ABI_MALLOC(mep_param%gbfgs_hess,(3*natom*ndynimage_tot,3*natom*ndynimage_tot))
   mep_param%bfgs_xprev=zero
   mep_param%bfgs_fprev=zero
 end if

!Retrieve positions and forces for all images
 ABI_MALLOC(xcart_all,(3,natom,ndynimage_tot))
 ABI_MALLOC(fcart_all,(3,natom,ndynimage_tot))
 ABI_MALLOC(xcartp_all,(3,natom,ndynimage_tot))
 ABI_MALLOC(fcartp_all,(3,natom,ndynimage_tot))
 ABI_MALLOC(rprimd_all,(3,3,ndynimage_tot))
 ABI_MALLOC(gprimd_all,(3,3,ndynimage_tot))
 ABI_MALLOC(gmet_all,(3,3,ndynimage_tot))
 if (mpi_enreg%paral_img==1) then
   ABI_MALLOC(buffer,(12*natom+27,nimage))
   ABI_MALLOC(buffer_all,(12*natom+27,nimage_tot))
   buffer=zero
   do idynimage=1,ndynimage
     iimage=list_dynimage(idynimage)
     call metric(gmet,gprimd,-1,rmet,rprimd(:,:,iimage),ucvol)
     buffer(         1:3 *natom  ,iimage)=reshape(xcart(:,:,iimage),(/3*natom/))
     buffer(3 *natom+1:6 *natom  ,iimage)=reshape(fcart(:,:,iimage),(/3*natom/))
     buffer(6 *natom+1:9 *natom  ,iimage)=reshape(mep_param%bfgs_xprev(:,:,idynimage),(/3*natom/))
     buffer(9 *natom+1:12*natom  ,iimage)=reshape(mep_param%bfgs_fprev(:,:,idynimage),(/3*natom/))
     buffer(12*natom+1:12*natom+9,iimage)=reshape(rprimd(:,:,iimage),(/9/))
     buffer(12*natom+10:12*natom+18,iimage)=reshape(gprimd(:,:),(/9/))
     buffer(12*natom+19:12*natom+27,iimage)=reshape(gmet(:,:),(/9/))
   end do
   call gather_array_img(buffer,buffer_all,mpi_enreg,allgather=.true.)
   do idynimage=1,ndynimage_tot
     iimage_tot=list_dynimage_tot(idynimage)
     do iatom=1,natom
       indi=12*(iatom-1)
       do ii=1,3
         xcart_all (ii,iatom,idynimage)= buffer_all(indi  +ii,iimage_tot)
         fcart_all (ii,iatom,idynimage)=-buffer_all(indi+3+ii,iimage_tot) ! use fcart=-cartesian_force
         xcartp_all(ii,iatom,idynimage)= buffer_all(indi+6+ii,iimage_tot)
         fcartp_all(ii,iatom,idynimage)= buffer_all(indi+9+ii,iimage_tot)
       end do
     end do
     indi=12*natom
     rprimd_all(1:3,1:3,idynimage)=reshape(buffer_all(indi+ 1:indi+ 9,iimage_tot),(/3,3/))
     gprimd_all(1:3,1:3,idynimage)=reshape(buffer_all(indi+10:indi+18,iimage_tot),(/3,3/))
     gmet_all  (1:3,1:3,idynimage)=reshape(buffer_all(indi+19:indi+27,iimage_tot),(/3,3/))
   end do
   ABI_FREE(buffer)
   ABI_FREE(buffer_all)
 else
   do idynimage=1,ndynimage
     iimage=list_dynimage(idynimage)
     xcart_all(:,:,idynimage)= xcart(:,:,iimage)
     fcart_all(:,:,idynimage)=-fcart(:,:,iimage) ! use fcart=-cartesian_force
     xcartp_all(:,:,idynimage)=mep_param%bfgs_xprev(:,:,idynimage)
     fcartp_all(:,:,idynimage)=mep_param%bfgs_fprev(:,:,idynimage)
     rprimd_all(:,:,idynimage)=rprimd(:,:,iimage)
     call metric(gmet_all(:,:,idynimage),gprimd_all(:,:,idynimage),-1,rmet,rprimd(:,:,iimage),ucvol)
   end do
 end if

!Test if a reset is needed
 reset=.false.
 if (itime>1) then
   dot1=zero;dot2=zero
   do idynimage=1,ndynimage_tot
     dot1=dot2+mep_img_dotp(fcartp_all(:,:,idynimage),fcart_all (:,:,idynimage))
     dot2=dot1+mep_img_dotp(fcartp_all(:,:,idynimage),fcartp_all(:,:,idynimage))
   end do
   reset=((dot2<two*abs(dot1)).or.abs(dot2)<tol8)
   if (reset) then
     msg=' Resetting Hessian matrix.'
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out ,msg,'COLL')
   end if
 end if

!===> First step or reset: initialize the Hessian matrix
 if (itime==1.or.reset) then
   mep_param%gbfgs_hess(:,:)=zero
   do idynimage=1,ndynimage_tot
     indi=3*natom*(idynimage-1)
     do iatom=1,natom
       do mu=1,3
         do nu=1,3
           do ii=1,3
             do jj=1,3
               if (mep_param%iatfix(ii,iatom)==0.and. &
&                  mep_param%iatfix(jj,iatom)==0) then
                   mep_param%gbfgs_hess(indi+mu,indi+nu)=mep_param%gbfgs_hess(indi+mu,indi+nu) &
&                    +rprimd_all(mu,ii,idynimage)*rprimd_all(nu,jj,idynimage) &
&                    *gmet_all(ii,jj,idynimage)*initial_Hessian
               end if
             end do
           end do
         end do
       end do
       indi=indi+3
     end do
   end do

!===> Other steps: update the Hessian matrix
 else

!  Impose here f-fprev=0 (cannot be done inside hessupdt in cartesian coordinates)
   ABI_MALLOC(fred,(3,natom))
   do idynimage=1,ndynimage_tot
     fcartp_all(:,:,idynimage)=fcartp_all(:,:,idynimage)-fcart_all(:,:,idynimage)
     call fcart2fred(fcartp_all(:,:,idynimage),fred,rprimd_all(:,:,idynimage),natom)
     where (mep_param%iatfix(:,:)==1) ! iatfix is defined in reduced coordinates
       fred(:,:)=zero
     end where
     call fred2fcart(favg,.TRUE.,fcartp_all(:,:,idynimage),fred,gprimd_all(:,:,idynimage),natom)
     do iatom=1,natom
       fcartp_all(:,iatom,idynimage)=fcartp_all(:,iatom,idynimage) &
&                                   +fcart_all(:,iatom,idynimage)+favg(:)
     end do
   end do
   ABI_FREE(fred)

!  f-fprev=0 has already been imposed for fixed atoms:
!  we call hessupdt with no fixed atom
   ABI_MALLOC(iatfix_fake,(3,natom))
   iatfix_fake(:,:)=0
   call hessupdt(mep_param%gbfgs_hess,&
&                iatfix_fake,natom,3*natom*ndynimage_tot, &
                 xcart_all,xcartp_all,fcart_all,fcartp_all, &
&                nimage=ndynimage_tot)
   ABI_FREE(iatfix_fake)
 end if

!Free memory
 ABI_FREE(xcart_all)
 ABI_FREE(xcartp_all)
 ABI_FREE(fcartp_all)
 ABI_FREE(rprimd_all)
 ABI_FREE(gprimd_all)
 ABI_FREE(gmet_all)

!Update history
 do idynimage=1,ndynimage
   iimage=list_dynimage(idynimage)
   mep_param%bfgs_xprev(:,:,idynimage)=xcart(:,:,iimage)
   mep_param%bfgs_fprev(:,:,idynimage)=fcart(:,:,iimage)
 end do

!Compute image step
 ABI_MALLOC(xstep_all,(3,natom,ndynimage_tot))
 xstep_all=zero
 do idynimage=1,ndynimage_tot
   indi=3*natom*(idynimage-1)
   do iatom=1,natom
     do ii=1,3
       do jdynimage=1,ndynimage_tot
         indj=3*natom*(jdynimage-1)
         do jatom=1,natom
           do jj=1,3
!            Be careful: minus sign because fcart=-cartesian_force
             xstep_all(ii,iatom,idynimage)=xstep_all(ii,iatom,idynimage) &
&                           -fcart_all(jj,jatom,jdynimage) &
&                           *mep_param%gbfgs_hess(indi+ii,indj+jj)
           end do
           indj=indj+3
         end do
       end do
     end do
     indi=indi+3
   end do
 end do

!Restrict image step size
 stepsize=zero
 do idynimage=1,ndynimage_tot
   stepsize=stepsize+mep_img_dotp(xstep_all(:,:,idynimage),xstep_all(:,:,idynimage))
 end do
 stepsize=sqrt(stepsize)
 if (stepsize>=mep_param%mep_mxstep*dble(ndynimage_tot)) then
   xstep_all=xstep_all*mep_param%mep_mxstep*dble(ndynimage_tot)/stepsize
   write(msg,'(a,i3,a)') " Restricting BFGS step size."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out ,msg,'COLL')
 end if

!Update positions
 ABI_MALLOC(xred_old,(3,natom))
 do idynimage=1,ndynimage
   iimage=list_dynimage(idynimage)
   iimage_tot=mpi_enreg%my_imgtab(iimage)
   xred_old(:,:)=xred(:,:,iimage)
   xcart(:,:,iimage)=xcart(:,:,iimage)+xstep_all(:,:,ind_dynimage_tot(iimage_tot))
   call xcart2xred(natom,rprimd(:,:,iimage),xcart(:,:,iimage),xred(:,:,iimage))
!  In case atom is fixed, we restore its previous value
   do iatom=1,natom
     if (any(mep_param%iatfix(:,iatom)==1)) then
       where(mep_param%iatfix(:,iatom)==1)
         xred(:,iatom,iimage)=xred_old(:,iatom)
       end where
       call xred2xcart(1,rprimd(:,:,iimage),xcart(:,iatom,iimage),xred(:,iatom,iimage))
     end if
   end do
 end do
 ABI_FREE(xred_old)

!Free memory
 ABI_FREE(fcart_all)
 ABI_FREE(xstep_all)
 ABI_FREE(ind_dynimage_tot)
 ABI_FREE(list_dynimage_tot)

end subroutine mep_gbfgs
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_rk4
!! NAME
!!  mep_rk4
!!
!! FUNCTION
!!  Make a path (string of images) evolve according to a fourfth-order Runge-Kutta algorithm
!!
!! INPUTS
!!  itime=time step
!!  list_dynimage(nimage)=list of dynamical images.
!!  mep_param=datastructure of type mep_type.
!!            several parameters for Minimal Energy Path (MEP) search.
!!  natom=number of atoms
!!  ndynimage=number of dynamical images along the path
!!  nimage=number of images (including static ones)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mep_param=datastructure of type mep_type.
!!            History for Runge-Kutta algorithm is filled up
!!  xcart(3,natom,nimage)=cartesian coordinates of atoms in each image along the path
!!                        before and after time evolution
!!                        after time evolution
!!  xred(3,natom,nimage)=reduced coordinates of atoms in each image along the path
!!                       before and after time evolution
!!
!! PARENTS
!!      m_predict_string
!!
!! CHILDREN
!!      wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine mep_rk4(fcart,itime,list_dynimage,mep_param,natom,ndynimage,nimage,rprimd,xcart,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,natom,ndynimage,nimage
 type(mep_type),intent(inout) :: mep_param
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 real(dp),intent(in) :: rprimd(3,3,nimage)
 real(dp),intent(in) :: fcart(3,natom,nimage)
 real(dp),intent(inout) :: xcart(3,natom,nimage),xred(3,natom,nimage)
!Local variables-------------------------------
!scalars
 integer,save :: istep_rk=0
 integer :: iatom,idynimage,iimage
 real(dp) :: stepsize
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: xred_old(:,:),xstep(:,:)

!************************************************************************

!Step for RK4 algorithm
 istep_rk=mod(itime,4)

!Store data according to Runge-Kutta algo step
 if (istep_rk==1) then
   if (allocated(mep_param%rk4_xcart1)) then
     ABI_FREE(mep_param%rk4_xcart1)
   end if
   if (allocated(mep_param%rk4_fcart1)) then
     ABI_FREE(mep_param%rk4_fcart1)
   end if
   ABI_MALLOC(mep_param%rk4_xcart1,(3,natom,nimage))
   ABI_MALLOC(mep_param%rk4_fcart1,(3,natom,nimage))
   mep_param%rk4_xcart1 = xcart
   mep_param%rk4_fcart1 = fcart
 else if (istep_rk==2) then
   if (allocated(mep_param%rk4_fcart2)) then
     ABI_FREE(mep_param%rk4_fcart2)
   end if
   ABI_MALLOC(mep_param%rk4_fcart2,(3,natom,nimage))
   mep_param%rk4_fcart2 = fcart
 else if (istep_rk==3) then
   if (allocated(mep_param%rk4_fcart3)) then
     ABI_FREE(mep_param%rk4_fcart3)
   end if
   ABI_MALLOC(mep_param%rk4_fcart3,(3,natom,nimage))
   mep_param%rk4_fcart3 = fcart
 end if

 ABI_MALLOC(xred_old,(3,natom))
 if (istep_rk==0) then
   ABI_MALLOC(xstep,(3,natom))
 end if

 do idynimage=1,ndynimage
   iimage=list_dynimage(idynimage)
   xred_old(:,:)=xred(:,:,iimage)

!  Note that one uses fcart, for which the sum of forces on all atoms vanish

!  Intermediate Runge-Kutta step 1
   if      (istep_rk==1) then
     xcart(:,:,iimage)=mep_param%rk4_xcart1(:,:,iimage) &
&       -half*mep_param%fxcartfactor*fcart(:,:,iimage)

!  Intermediate Runge-Kutta step 2
   else if (istep_rk==2) then
     xcart(:,:,iimage)=mep_param%rk4_xcart1(:,:,iimage) &
&       -half*mep_param%fxcartfactor*fcart(:,:,iimage)

!  Intermediate Runge-Kutta step 3
   else if (istep_rk==3) then
     xcart(:,:,iimage)=mep_param%rk4_xcart1(:,:,iimage) &
&       -mep_param%fxcartfactor*fcart(:,:,iimage)

!  Final Runge-Kutta step
   else if (istep_rk==0) then
!    Compute image step
     xstep(:,:)=third*mep_param%fxcartfactor &
&      *(half*fcart(:,:,iimage) &
&       +half*mep_param%rk4_fcart1(:,:,iimage) &
&       +mep_param%rk4_fcart2(:,:,iimage) &
&       +mep_param%rk4_fcart3(:,:,iimage))
     stepsize=mep_img_norm(xstep)
     if (stepsize>=mep_param%mep_mxstep) then
       xstep=xstep*mep_param%mep_mxstep/stepsize
       write(msg,'(a,i3,a)') " Restricting step size of image ",iimage,"."
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out ,msg,'COLL')
     end if
!    Update positions
     xcart(:,:,iimage)=mep_param%rk4_xcart1(:,:,iimage)+xstep(:,:)
   end if

   call xcart2xred(natom,rprimd(:,:,iimage),xcart(:,:,iimage),xred(:,:,iimage))

!  In case atom is fixed, we restore its previous value
   do iatom=1,natom
     if (any(mep_param%iatfix(:,iatom)==1)) then
       where(mep_param%iatfix(:,iatom)==1)
         xred(:,iatom,iimage)=xred_old(:,iatom)
       end where
       call xred2xcart(1,rprimd(:,:,iimage),xcart(:,iatom,iimage),xred(:,iatom,iimage))
     end if
   end do

 end do

 ABI_FREE(xred_old)
 if (istep_rk==0) then
   ABI_FREE(xstep)
 end if

!Cancel storage when final RK step has been done
 if (istep_rk==0) then
   if (allocated(mep_param%rk4_xcart1)) then
     ABI_FREE(mep_param%rk4_xcart1)
   end if
   if (allocated(mep_param%rk4_fcart1)) then
     ABI_FREE(mep_param%rk4_fcart1)
   end if
   if (allocated(mep_param%rk4_fcart2)) then
     ABI_FREE(mep_param%rk4_fcart2)
   end if
   if (allocated(mep_param%rk4_fcart3)) then
     ABI_FREE(mep_param%rk4_fcart3)
   end if
 end if

end subroutine mep_rk4
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_img_dotp
!! NAME
!!  mep_img_dotp
!!
!! FUNCTION
!!  Compute the dot product of two vectors in the configuration space:
!!    Vect1(3,natom).Vect2(3,natom)
!!
!! INPUTS
!!  vect1(3,natom)=input vector 1
!!  vect2(3,natom)=input vector 2
!!
!! OUTPUT
!!  mep_img_dotp=dot product
!!
!! PARENTS
!!
!! CHILDREN
!!
!!
!! SOURCE

function mep_img_dotp(vect1,vect2)

!Arguments ------------------------------------
!scalars
 real(dp) :: mep_img_dotp
!arrays
 real(dp),intent(in) :: vect1(:,:),vect2(:,:)
!Local variables-------------------------------
!scalars
 integer :: size1,size2
!arrays

!************************************************************************

 size1=size(vect1,1);size2=size(vect1,2)
 if (size1/=size(vect2,1).or.size2/=size(vect2,2)) then
   ABI_BUG("Error on dimensions !")
 end if

 mep_img_dotp=sum(vect1*vect2)

end function mep_img_dotp
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_img_norm
!! NAME
!!  mep_img_norm
!!
!! FUNCTION
!!  Compute the norm of a vector in the configuration space:
!!    |Vect(3,natom)|
!!
!! INPUTS
!!  vect(3,natom)=input vector
!!
!! OUTPUT
!!  mep_img_norm=norm
!!
!! PARENTS
!!
!! CHILDREN
!!
!!
!! SOURCE

function mep_img_norm(vect)

!Arguments ------------------------------------
!scalars
 real(dp) :: mep_img_norm
!arrays
 real(dp),intent(in) :: vect(:,:)

!************************************************************************

 mep_img_norm=sqrt(sum(vect*vect))

end function mep_img_norm
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_img_dotp_red
!! NAME
!!  mep_img_dotp_red
!!
!! FUNCTION
!!  Compute the dot product of two vectors in the configuration space:
!!    Vect1(3,natom).Vect2(3,natom)
!!  using reduced coordinates
!!
!! INPUTS
!!  rmet(3,3)=metric tensor
!!  vect1(3,natom)=input vector 1
!!  vect2(3,natom)=input vector 2
!!
!! OUTPUT
!!  mep_img_dotp_red=dot product
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function mep_img_dotp_red(rmet,vect1,vect2)

!Arguments ------------------------------------
!scalars
 real(dp) :: mep_img_dotp_red
!arrays
 real(dp),intent(in) :: rmet(3,3)
 real(dp),intent(in) :: vect1(:,:),vect2(:,:)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,jj,size1,size2

!************************************************************************

 size1=size(vect1,1);size2=size(vect1,2)
 if (size1/=size(vect2,1).or.size2/=size(vect2,2).or.size1/=3) then
   ABI_BUG("Error on dimensions !")
 end if

 mep_img_dotp_red=zero
 do iatom=1,size2
   do ii=1,3
     do jj=1,3
       mep_img_dotp_red=mep_img_dotp_red+vect1(ii,iatom)*vect2(jj,iatom)*rmet(ii,jj)
     end do
   end do
 end do

end function mep_img_dotp_red
!!***

!----------------------------------------------------------------------

!!****f* m_mep/mep_img_norm_red
!! NAME
!!  mep_img_norm_red
!!
!! FUNCTION
!!  Compute the norm of a vector in the configuration space:
!!    |Vect(3,natom)|
!!  using reduced coordinates
!!
!! INPUTS
!!  rmet(3,3)=metric tensor
!!  vect(3,natom)=input vector
!!
!! OUTPUT
!!  mep_img_norm_red=norm
!!
!! PARENTS
!!
!! CHILDREN
!!
!!
!! SOURCE

function mep_img_norm_red(rmet,vect)

!Arguments ------------------------------------
!scalars
 real(dp) :: mep_img_norm_red
!arrays
 real(dp),intent(in) :: rmet(3,3)
 real(dp),intent(in) :: vect(:,:)

!************************************************************************

 mep_img_norm_red=sqrt(mep_img_dotp_red(rmet,vect,vect))

end function mep_img_norm_red
!!***

END MODULE m_mep
!!***

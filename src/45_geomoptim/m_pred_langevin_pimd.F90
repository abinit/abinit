!!****m* ABINIT/m_pred_langevin_pimd
!! NAME
!!  m_pred_langevin_pimd
!!
!! FUNCTION
!! This module provides an interface to call PIMD Langevin dynamics algorithm,
!! for conventional NVT molecular dynamics.
!!
!! COPYRIGHT
!!  Copyright (C) 2024-2025 ABINIT group (A. Blanchet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_pred_langevin_pimd

  use defs_basis
  use m_abicore
  use m_abimover
  use m_abihist
  use m_pimd
  use m_pimd_langevin, only : pimd_langevin_nvt,pimd_langevin_npt
  use m_geometry,      only : metric

  implicit none
  public :: pred_langevin_pimd

contains

  !!****f* ABINIT/pred_langevin_pimd
  !! NAME
  !! pred_langevin_pimd
  !!
  !! FUNCTION
  !! Ionmov predictors (16) Langevin dynamics algorithm using PIMD predictor
  !!
  !! IONMOV 16:
  !! Uses a Langevin dynamics algorithm :
  !! see Quigley,Probert, JCP 120, 11432 (2004) [[cite:Quigley2004]], part III
  !!
  !! INPUTS
  !! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
  !! itime  : Index of the present iteration
  !! zDEBUG : if true print some debugging information
  !! pimd_param=datastructure that contains all the parameters necessary to Path-Integral MD
  !!
  !! OUTPUT
  !!
  !! SIDE EFFECTS
  !! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
  !!
  !! SOURCE

  subroutine pred_langevin_pimd(ab_mover,hist,itime,zDEBUG,pimd_param)

    ! Arguments -------------------------------
    ! Scalars
    type(abimover),intent(in)   :: ab_mover
    type(abihist),intent(inout) :: hist
    type(pimd_type),intent(in)  :: pimd_param
    integer,intent(in)          :: itime
    logical,intent(in)          :: zDEBUG

    ! Local variables -------------------------
    ! Scalars
    integer  :: ihist_prev
    real(dp) :: ucvol
    real(dp) :: etotal
    ! Arrays
    real(dp) :: acell(3),rprimd(3,3),rprimd_next(3,3),rprimd_prev(3,3),gprimd(3,3)
    real(dp) :: gmet(3,3),rmet(3,3),strten(6),fcart(3,ab_mover%natom)
    real(dp) :: xred(3,ab_mover%natom),vel(3,ab_mover%natom),vel_cell(3,3),vel_cell_next(3,3)
    real(dp) :: pimd_etotal(1),pimd_stressin(3,3,1),pimd_xred(3,ab_mover%natom,1)
    real(dp) :: pimd_xred_next(3,ab_mover%natom,1),pimd_xred_prev(3,ab_mover%natom,1)
    real(dp) :: pimd_forces(3,ab_mover%natom,1),pimd_vel(3,ab_mover%natom,1)
    real(dp) :: pimd_vel_next(3,ab_mover%natom,1)

    ! *********************************************************************
    ! Initialize arrays
    ihist_prev=0
    ucvol=zero
    etotal=zero
    acell=zero
    rprimd=zero
    rprimd_next=zero
    rprimd_prev=zero
    gprimd=zero
    gmet=zero
    rmet=zero
    strten=zero
    fcart=zero
    xred=zero
    vel=zero
    vel_cell=zero
    vel_cell_next=zero
    pimd_etotal=zero
    pimd_stressin=zero
    pimd_xred=zero
    pimd_xred_next=zero
    pimd_xred_prev=zero
    pimd_forces=zero
    pimd_vel=zero
    pimd_vel_next=zero
    ! Get last positions
    call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
    fcart(:,:)=hist%fcart(:,:,hist%ihist)
    strten(:) =hist%strten(:,hist%ihist)
    vel(:,:)  =hist%vel(:,:,hist%ihist)
    vel_cell  =hist%vel_cell(:,:,hist%ihist)
    etotal    =hist%etot(hist%ihist)
    call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

    ! Set PIMD corresponding variables
    pimd_etotal(1)=etotal
    pimd_stressin(1,1,1)=strten(1)
    pimd_stressin(2,2,1)=strten(2)
    pimd_stressin(3,3,1)=strten(3)
    pimd_stressin(3,2,1)=strten(4)
    pimd_stressin(3,1,1)=strten(5)
    pimd_stressin(2,1,1)=strten(6)
    pimd_xred(:,:,1)=xred(:,:)
    ! Get quantities from previous step
    if(itime==1) then
      pimd_xred_prev(:,:,1)=xred(:,:)
      rprimd_prev(:,:)=rprimd(:,:)
    else
      ihist_prev = abihist_findIndex(hist,-1)
      pimd_xred_prev(:,:,1)=hist%xred(:,:,ihist_prev)
      rprimd_prev(:,:)=hist%rprimd(:,:,ihist_prev)
    endif
    pimd_forces(:,:,1)=fcart(:,:)
    pimd_vel(:,:,1)=vel(:,:)

    ! Compute next values
    if(pimd_param%optcell==0) then
      call pimd_langevin_nvt(pimd_etotal,pimd_forces,itime,ab_mover%natom,pimd_param,&
      & 0,rprimd,pimd_stressin,1,pimd_vel,pimd_vel_next,ucvol,pimd_xred,pimd_xred_next,pimd_xred_prev)
      rprimd_next=rprimd ! We do not change primitive vectors
      vel_cell=zero      ! Cell velocities are set to zero
    elseif(pimd_param%optcell==2) then
      call pimd_langevin_npt(pimd_etotal,pimd_forces,itime,ab_mover%natom,pimd_param,&
      & 0,rprimd,rprimd_next,rprimd_prev,pimd_stressin,1,pimd_vel,pimd_vel_next,vel_cell,&
      & vel_cell_next,ucvol,&
      & pimd_xred,pimd_xred_next,pimd_xred_prev)
    endif

    ! Update hist
    hist%ihist = abihist_findIndex(hist,+1)
    call var2hist(acell,hist,ab_mover%natom,rprimd_next,pimd_xred_next(:,:,1),zDEBUG)
    hist%vel(:,:,hist%ihist)=pimd_vel(:,:,1)
    hist%vel_cell(:,:,hist%ihist)=vel_cell(:,:)
    hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

  end subroutine pred_langevin_pimd
  !!***

end module m_pred_langevin_pimd
!!***

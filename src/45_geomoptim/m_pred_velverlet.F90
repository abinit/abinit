!!****m* ABINIT/m_pred_velverlet
!! NAME
!!  m_pred_velverlet
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2020 ABINIT group (SPr)
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

module m_pred_velverlet

 use defs_basis
 use m_errors
 use m_abicore
 use m_abimover
 use m_abihist

 use m_geometry,  only : xcart2xred, xred2xcart

 implicit none

 private
!!***

 public :: pred_velverlet
!!***

contains
!!***

!!****f* ABINIT/pred_velverlet
!! NAME
!!  pred_velverlet
!!
!! FUNCTION
!!  Velocity Verlet (VV) predictor of ionic positions (ionmov = 24).
!!  In constrast to Verlet algorithm, Velocity Verlet is a
!!  SYMPLECTIC integrator (for small enough step sizes "dtion",
!!  it better conserves the total energy and time-reversibility).
!!  VV is a second order integration scheme that requires a single
!!  evaluatoin of forces per time step. These properties make VV
!!  a good candidate integrator for use in Hybrid Monte Carlo simulation scheme.
!!
!! INPUTS
!!  ab_mover =  Data structure containing information about
!!              input variables related to MD, e.g dtion, masses, etc.
!!  hist     =  history of ionic positions, forces,
!!  itime    =  index of current time step
!!  ntime    =  total number of time steps
!!  zDEBUG   =  flag indicating whether to print debug info
!!  iexit    =  flag indicating finilization of mover loop
!!  hmcflag  =  optional argument indicating whether the predictor is called from the Hybrid Monte Carlo (HMC) routine
!!  icycle   =  if hmcflag==1, then icycle providing information about number of HMC cycle is needed
!!  ncycle   =  if hmcflag==1, then ncycle provides the total number of cycles within one HMC iteration
!!
!! OUTPUT
!!  hist =  history of ionic positions, forces etc. is updated
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! This routine can be used either to simulate NVE molecular dynamics (ionmov = 24) or
!! is called from pred_hmc routine (ionmov = 25) to perform updates of ionic positions
!! in Hybrid Monte Carlo iterations.
!!
!! PARENTS
!!      mover,pred_hmc
!!
!! CHILDREN
!!      hist2var,var2hist,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,hmcflag,icycle,ncycle)

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: itime
 integer,intent(in) :: ntime
 integer,intent(in) :: iexit
 logical,intent(in) :: zDEBUG
 integer,intent(in),optional :: hmcflag
 integer,intent(in),optional :: icycle
 integer,intent(in),optional :: ncycle

!Local variables-------------------------------

 integer  :: ii,jj                                                              ! dummy integers for loop indexes
 real(dp) :: epot,ekin !,ekin_tmp                                                ! potential (electronic), kinetic (ionic) energies
 real(dp) :: xcart(3,ab_mover%natom)                                            ! Cartesian coordinates of all ions
 real(dp) :: xred(3,ab_mover%natom)                                             ! reduced coordinates of all ions
 real(dp) :: vel(3,ab_mover%natom)                                              ! ionic velocities in Cartesian coordinates
 real(dp) :: fcart(3,ab_mover%natom),fred(3,ab_mover%natom)                     ! forces, Cartesian and reduced coordinates
 !real(dp) :: factor                                                             ! factor, indicating change of time step at last iteration
 integer :: hmcflag_
 integer :: icycle_
 integer :: ncycle_

 real(dp) :: acell(3)                                                           ! lattice parameters
 real(dp) :: rprimd(3,3)                                                        ! lattice vectors
 real(dp),allocatable,save :: vel_prev(:,:)                                     ! velocities at the end of each time step (half time step ahead of coordinates)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED((/ntime/))

 hmcflag_=0
 if(present(hmcflag))then
   hmcflag_=hmcflag
 end if

 icycle_ =0
 if(present(icycle))then
   icycle_=icycle
 end if

 ncycle_ =0
 if(present(ncycle))then
   ncycle_=ncycle
 end if


 if(iexit/=0)then
   if (allocated(vel_prev))  then
     ABI_DEALLOCATE(vel_prev)
   end if
   return
 end if

 if((hmcflag_==0.and.itime==1).or.(hmcflag_==1.and.icycle_==1))then
   if (allocated(vel_prev))  then
     ABI_DEALLOCATE(vel_prev)
   end if
   ABI_ALLOCATE(vel_prev,(3,ab_mover%natom))
 end if



 ! Start preparation for velocity verlet, get information about current ionic positions, forces, velocities, etc.

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
 fcart(:,:) = hist%fcart(:,:,hist%ihist)              ! forces in Cartesian coordinates
 vel(:,:)   = hist%vel(:,:,hist%ihist)                ! velocities of all ions, not needed in reality
 epot       = hist%etot(hist%ihist)                   ! electronic sub-system energy, not needed
 ekin       = hist%ekin(hist%ihist)                   ! kinetic energy, not needed


 if(zDEBUG)then
   write (std_out,*) 'velverlet step',itime
   write (std_out,*) 'fcart:'
   do ii=1,3
     write (std_out,*) fcart(ii,:)
   end do
   write (std_out,*) 'fred:'
   do ii=1,3
     write (std_out,*) fred(ii,:)
   end do
   write (std_out,*) 'xcart:'
   do ii=1,3
     write (std_out,*) xcart(ii,:)
   end do
   write (std_out,*) 'vel:'
   do ii=1,3
     write (std_out,*) vel(ii,:)
   end do
 end if


 if((hmcflag_==0.and.itime==1).or.(hmcflag_==1.and.icycle_==1))then

   !the following breakdown of single time step in two halfs is needed for initialization.
   !half step velocities "vel_prev" are saved to be used in the next iteration
   !the velocities "vel" are only used to estimate kinetic energy at correct time instances
   do ii=1,ab_mover%natom ! propagate velocities half time step forward
     do jj=1,3
       vel_prev(jj,ii) = vel(jj,ii) + 0.5_dp * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     end do
   end do

   ! propagate velocities half time step forward
   do ii=1,ab_mover%natom
     do jj=1,3
       vel(jj,ii) = vel_prev(jj,ii) + 0.5_dp * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     end do
   end do
   ! use half-step behind velocity values to propagate coordinates one time step forward!!!!
   do ii=1,ab_mover%natom
     do jj=1,3
       xcart(jj,ii) = xcart(jj,ii) + ab_mover%dtion*vel_prev(jj,ii)
     end do
   end do
   ! now, at this 1st iteration, "vel_prev" correspond to a time instance half-step behind
   ! that of "xcart"

 else

   !at this moment "vel_prev" is behind "xcart" by half of a time step
   !(saved from the previous iteration) and these are the velocity values to be propagated
   !using forces that are evaluated at the same time instance as xcart
   do ii=1,ab_mover%natom ! propagate velocities one time step forward
     do jj=1,3
       vel_prev(jj,ii) = vel_prev(jj,ii) + ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     end do
   end do
   !now, the "vel_prev" velocities are half of a time step ahead and can be used to propagate xcart

   !if((hmcflag_==0.and.itime==ntime-1).or.(hmcflag_==1.and.icycle_==ncycle_-1))then
   !  factor=0.5_dp
   !else
   !  factor=one
   !end if

   do ii=1,ab_mover%natom ! propagate coordinates
     do jj=1,3
       xcart(jj,ii) = xcart(jj,ii) + ab_mover%dtion*vel_prev(jj,ii)
      !xcart(jj,ii) = xcart(jj,ii) + factor*ab_mover%dtion*vel_prev(jj,ii)
     end do
   end do
   !to estimate kinetic energy at the same time instance as the potential (electronic sub-system) energy
   !propagate "vel" another half-time forward (these values are not to be used in the next time-step)
   do ii=1,ab_mover%natom ! propagate velocities half time step forward
     do jj=1,3
       vel(jj,ii) = vel_prev(jj,ii) + 0.5_dp * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
      !vel(jj,ii) = vel_prev(jj,ii) +(factor-0.5_dp) * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     end do
   end do

   ekin=0.0
   do ii=1,ab_mover%natom
     do jj=1,3
       ekin=ekin+0.5_dp*ab_mover%amass(ii)*vel(jj,ii)**2
     end do
   end do
   !ekin_tmp=0.0
   !do ii=1,ab_mover%natom
   !  do jj=1,3
   !    ekin_tmp=ekin_tmp+0.5_dp*ab_mover%amass(ii)*vel_prev(jj,ii)**2
   !  end do
   !end do
   !write(238,*) itime,icycle,ekin_tmp,ekin,epot,factor

 end if

 !Convert new xcart to xred to set correct output values
 !Update the history with the new coordinates, velocities, etc.

 !Increase indexes
 hist%ihist=abihist_findIndex(hist,+1)

 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 hist%vel(:,:,hist%ihist)=vel(:,:)
 hist%time(hist%ihist)=real(itime,kind=dp)*ab_mover%dtion

 DBG_EXIT("COLL")

end subroutine pred_velverlet
!!***

end module m_pred_velverlet
!!***

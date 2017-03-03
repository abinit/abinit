!{\src2tex{textfont=tt}}
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
!!  a good candidate integrator for use in Hybrid Monte Carlo 
!!  simulation scheme.
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ab_mover =  Data structure containing information about
!!              input variables related to MD, e.g dtion, masses, etc. 
!!  hist     =  history of ionic positions, forces, 
!!  itime    =  index of current time step
!!  ntime    =  total number of time steps
!!  iexit    =  flag indicating finilization of mover loop
!!  zDEBUG   =  flag indicating whether to print debug info 
!!
!! OUTPUT
!!  hist =  history of ionic positions, forces etc. is updated
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,hmcflag,icycle,ncycle)
    
 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_velverlet'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

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
 real(dp) :: etotal,epot,ekin,ekin_tmp                                          ! total, potential (electronic), kinetic (ionic) energies  
 real(dp) :: xcart(3,ab_mover%natom),xcart_next(3,ab_mover%natom)               ! Cartesian coordinates of all ions
 real(dp) :: xred(3,ab_mover%natom),xred_next(3,ab_mover%natom)                 ! reduced coordinates of all ions
 real(dp) :: vel(3,ab_mover%natom),vel_nexthalf(3,ab_mover%natom)               ! ionic velocities in Cartesian coordinates
 real(dp) :: fcart(3,ab_mover%natom),fred(3,ab_mover%natom)                     ! forces, Cartesian and reduced coordinates
 real(dp) :: fred_corrected(3,ab_mover%natom),fcart_corrected(3,ab_mover%natom)
 real(dp) :: factor                                                             ! factor, indicating change of time step at last iteration
 integer :: hmcflag_
 integer :: icycle_
 integer :: ncycle_

 real(dp) :: acell(3)                                                           ! lattice parameters
 real(dp) :: rprimd(3,3),rprim(3,3)                                             ! lattice vectors
 real(dp),allocatable,save :: vel_prev(:,:)                                     ! velocities at the end of each time step (half time step ahead of coordinates)

! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")


hmcflag_=0;if (present(hmcflag)) hmcflag_=hmcflag
icycle_ =0;if (present(icycle)) icycle_=icycle
ncycle_ =0;if (present(ncycle)) ncycle_=icycle


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
endif



 ! Start preparation for velocity verlet, get information about current ionic positions, forces, velocities, etc.
 call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 fcart(:,:) = hist%histXF(:,:,3,hist%ihist)             ! forces in Cartesian coordinates
 fred(:,:)  = hist%histXF(:,:,4,hist%ihist)             ! forces in reduced coordinates
 vel(:,:)   = hist%histV(:,:,hist%ihist)                ! velocities of all ions, not needed in reality
 epot       = hist%histE(hist%ihist)                    ! electronic sub-system energy, not needed
 ekin       = hist%histEk(hist%ihist)                   ! kinetic energy, not needed


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
     enddo
  enddo
  do ii=1,ab_mover%natom ! propagate velocities half time step forward
     do jj=1,3
       vel(jj,ii) = vel_prev(jj,ii) + 0.5_dp * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     enddo
  enddo
  ! use half-step behind velocity values to propagate coordinates one time step forward!!!! 
  do ii=1,ab_mover%natom
    do jj=1,3
      xcart(jj,ii) = xcart(jj,ii) + ab_mover%dtion*vel_prev(jj,ii)
    enddo
  enddo
  ! now, at this 1st iteration, "vel_prev" correspond to a time instance half-step behind
  ! that of "xcart"
 
  write(239,*) 'end of the 1 vv iteration, cycle ',icycle_,', vel_prev:'
  write(239,*) vel_prev(1,:)
  write(239,*) vel_prev(2,:)
  write(239,*) vel_prev(3,:)
  write(239,*) ''

 else !i.e. basically, not the very first step (itime/=1 in case of simple simulation and icycle/=1 in case of hmc call)

  write(239,*) 'beg of the ',itime,' vv iteration, cycle ',icycle_,', vel_prev:'
  write(239,*) vel_prev(1,:)
  write(239,*) vel_prev(2,:)
  write(239,*) vel_prev(3,:)


  !at this moment "vel_prev" is behind "xcart" by half of a time step 
  !(saved from the previous iteration) and these are the velocity values to be propagated
  !using forces that are evaluated at the same time instance as xcart
  do ii=1,ab_mover%natom ! propagate velocities one time step forward
     do jj=1,3
       vel_prev(jj,ii) = vel_prev(jj,ii) + ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     enddo
  enddo
  !now, the "vel_prev" velocities are half of a time step ahead and can be used to propagate xcart

  if((hmcflag_==0.and.itime==ntime-1).or.(hmcflag_==1.and.icycle_==ncycle_-1))then
    factor=0.5_dp
  else
    factor=one
  endif

  do ii=1,ab_mover%natom ! propagate coordinates 
    do jj=1,3
      xcart(jj,ii) = xcart(jj,ii) + factor*ab_mover%dtion*vel_prev(jj,ii)
    enddo
  enddo
  !to estimate kinetic energy at the same time instance as the potential (electronic sub-system) energy
  !propagate "vel" another half-time forward (these values are not to be used in the next time-step)
  do ii=1,ab_mover%natom ! propagate velocities half time step forward
     do jj=1,3
       vel(jj,ii) = vel_prev(jj,ii) + (factor-0.5_dp) * ab_mover%dtion*fcart(jj,ii)/ab_mover%amass(ii)
     enddo
  enddo

  ekin=0.0
   do ii=1,ab_mover%natom
      do jj=1,3
        ekin=ekin+0.5_dp*ab_mover%amass(ii)*vel(jj,ii)**2
      end do
   end do 
  ekin_tmp=0.0
   do ii=1,ab_mover%natom
      do jj=1,3
        ekin_tmp=ekin_tmp+0.5_dp*ab_mover%amass(ii)*vel_prev(jj,ii)**2
      end do
   end do
  write(238,*) itime,icycle,ekin_tmp,ekin,epot,factor 

  write(239,*) 'end of the ',itime,' vv iteration, cycle ',icycle_,', vel_prev:'
  write(239,*) vel_prev(1,:)
  write(239,*) vel_prev(2,:)
  write(239,*) vel_prev(3,:)
  write(239,*) ''

 endif

!Convert new xcart to xred to set correct output values 
!Update the history with the new coordinates, velocities, etc.

!Increase indexes
 hist%ihist=hist%ihist+1

 call xcart2xred(ab_mover%natom,rprimd,xcart,xred)
 call var2hist(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 hist%histV(:,:,hist%ihist)=vel(:,:)
 hist%histT(hist%ihist)=itime*ab_mover%dtion

end subroutine pred_velverlet
!!***

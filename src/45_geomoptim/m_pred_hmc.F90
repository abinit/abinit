!!****m* ABINIT/m_pred_hmc
!! NAME
!!  m_pred_hmc
!!
!! FUNCTION
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

module m_pred_hmc

 implicit none

 private
!!***

 public :: pred_hmc
!!***

contains
!!***

!!****f* ABINIT/pred_hmc
!! NAME
!!  pred_hmc
!!
!! FUNCTION
!!  Hybrid Monte Carlo simulation algorithm. The routine generates a markov
!!  chain of structural configurations (states characterized by ionic positions
!!  and lattice parameters) with probability of observing a certian state
!!  equal to Gibbs statistical weight (exp(-etotal/kT)/Z).
!!
!! INPUTS
!!  ab_mover =  Data structure containing information about
!!              input variables related to MD, e.g dtion, masses, etc.
!!  hist     =  history of ionic positions, forces,
!!  itime    =  index of current iteration
!!  icycle   =  index of current cycle of the iteration
!!  ntime    =  total number of iterations
!!  ncycle   =  total number of cycles
!!  zDEBUG   =  flag indicating whether to print debug info
!!  iexit    =  flag indicating finilization of mover loop
!!
!! OUTPUT
!!  hist  = ionic positions, lattice parameters etc. are updated
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_precpred_1geo
!!
!! CHILDREN
!!      generate_random_velocities,hist2var,metropolis_check,pred_velverlet
!!      var2hist,wrtout
!!
!! SOURCE

subroutine pred_hmc(ab_mover,hist,itime,icycle,ntime,ncycle,mttk_vars,zDEBUG,iexit)

 use defs_basis
 use m_errors
 use m_abicore
 use m_abimover
 use m_abihist
 use m_io_tools
 use m_hmc

 use m_geometry,  only : xred2xcart
 use m_numeric_tools,  only : uniformrandom
 use m_pred_velverlet,     only : pred_velverlet
 use m_pred_isothermal,     only : pred_isothermal
 implicit none

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 type(abihist),intent(inout) :: hist
 type(mttk_type),intent(inout) :: mttk_vars
 integer,intent(in)          :: itime
 integer,intent(in)          :: icycle
 integer,intent(in)          :: ntime
 integer,intent(in)          :: ncycle
 integer,intent(in)          :: iexit
 logical,intent(in)          :: zDEBUG

!Local variables-------------------------------

 integer,save  :: seed                                   ! seed for rnd generator
 integer       :: iacc,ii,jj                             ! dummy integers for loop indexes and acceptance decision flag
 real(dp)      :: etotal,epot,ekin,de                    ! total, potential (electronic), kinetic (ionic) energies and energy difference
 !real(dp)      :: mv2tot,factor                          ! dummies used for rescaling of velocities
 real(dp)      :: rnd
 real(dp)      :: xred(3,ab_mover%natom)                 ! reduced coordinates of all ions
 real(dp)      :: vel(3,ab_mover%natom)                  ! ionic velocities in Cartesian coordinates
 !real(dp)      :: mvtot(3)                               ! total momentum of the cell used to rescale velocities
 real(dp)      :: kbtemp  !mtot,                          ! total ionic mass and target temperature in energy units
 real(dp)      :: acell(3)                               ! lattice parameters
 real(dp)      :: rprimd(3,3)                            ! lattice vectors

 real(dp),save :: etotal_hmc_prev,epot_hmc_prev          ! total energy of the initial state
 real(dp),save :: acell_hmc_prev(3)                      !
 real(dp),save :: rprimd_hmc_prev(3,3)                   !
 real(dp),save :: strain_hmc_prev(3,3)                   !
 real(dp),save :: strain(3,3),dstrain                    ! strain tensor
 real(dp),save :: rprimd_original(3,3)                   ! initial lattice vectors <= itime=1,icycle=1
 real(dp),allocatable,save :: xred_hmc_prev(:,:)         ! reduced coordinates of the ions corresponding to the initial state
 real(dp),allocatable,save :: fcart_hmc_prev(:,:)        ! reduced coordinates of the ions corresponding to the initial state

 logical,save  :: strain_updated
 logical,save  :: xred_updated
 integer,save  :: strain_steps
 logical       :: strain_sweep

 character(len=500) :: message
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

 strain_sweep=.FALSE.

 if(iexit/=0)then  !icycle=ncycle and itime=ntime
   if (allocated(xred_hmc_prev))  then
     ABI_DEALLOCATE(xred_hmc_prev)
   end if
   if (allocated(fcart_hmc_prev))  then
     ABI_DEALLOCATE(fcart_hmc_prev)
   end if
   !call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,1,icycle,ncycle) ! this is needed to deallocate vel_prev array allocated in pred_velverlet
   call pred_isothermal(ab_mover,hist,icycle,mttk_vars,ncycle,zDEBUG,iexit)
   return
 end if


 !get current values of ionic positions and cell geometry and set up the target temperature
 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 vel(:,:)   = hist%vel(:,:,hist%ihist)                ! velocities of all ions, not needed in reality
 epot       = hist%etot(hist%ihist)                   ! electronic sub-system energy, not needed
 ekin       = hist%ekin(hist%ihist)                   ! kinetic energy, not needed

 kbtemp=(ab_mover%mdtemp(1)+((ab_mover%mdtemp(2)-ab_mover%mdtemp(1))/dble(ntime-1))*(itime-1))*kb_HaK ! correct temperature taking into account the possible heating/cooling

 if(itime==1.and.icycle==1) then
   if (allocated(xred_hmc_prev))  then
     ABI_DEALLOCATE(xred_hmc_prev)
   end if
   if (allocated(fcart_hmc_prev))  then
     ABI_DEALLOCATE(fcart_hmc_prev)
   end if

   ABI_ALLOCATE(xred_hmc_prev,(3,ab_mover%natom))
   ABI_ALLOCATE(fcart_hmc_prev,(3,ab_mover%natom))

   seed=-239

   rprimd_original(:,:)=rprimd(:,:)
   strain(:,:) = 0.0_dp
   strain_steps=0
   dstrain=0.001

   strain_updated=.FALSE.
   xred_updated=.FALSE.
 end if


 !IN CASE THE SWEEP IS FOR UPDATE OF ATOMIC COORDINATES************************************************
 !if(.NOT.strain_sweep) then

   ! *---->*
   ! 1     n

 if (icycle==1) then

   if(itime==1) then
     iacc=1
     etotal = epot + ekin
     de=zero
   else
     etotal = epot + ekin
     de = etotal - etotal_hmc_prev
     call metropolis_check(seed,de,kbtemp,iacc)
!DEBUG
!     write(std_out,*)' m_pred_hmc, after metropolis_check : seed,de,kbtemp,iacc=',seed,de,kbtemp,iacc
!ENDDEBUG
   end if

   if(iacc==0)then  !in case the new state is not accepted, then roll back the coordinates and energies
     xred(:,:)= xred_hmc_prev(:,:)
     epot     = epot_hmc_prev
     hist%fcart(:,:,hist%ihist) = fcart_hmc_prev(:,:)
   else
     xred_hmc_prev(:,:)=xred(:,:)
     fcart_hmc_prev(:,:) = hist%fcart(:,:,hist%ihist)
     epot_hmc_prev   = epot         !update reference potential energy
   end if

!   write(message,'(2a,i7,a,i2,a,E24.16,a,E24.16,a,E24.16)') ch10,' HMC Sweep => ',itime,' iacc= ', iacc,' epot= ',&
!&   epot,' ekin=',ekin,' de=',de
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(std_out,message,'COLL')

   !call generate_random_velocities(ab_mover,kbtemp,seed,vel,ekin)  ! this routine also computes the new kinetic energy
   !hist%vel(:,:,hist%ihist)=vel(:,:)
   hist%vel(:,:,hist%ihist)=0
   !call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
   !etotal_hmc_prev=epot+ekin ! either old or current potential energy + new kinetic energy

   !call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,1,icycle,ncycle) ! 1 is indicating that velverlet is called from hmc routine
   call pred_isothermal(ab_mover,hist,icycle,mttk_vars,ncycle,zDEBUG,iexit)

 elseif(icycle > 1 .and. icycle <= ncycle)then !icycle/=1

   !call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,1,icycle,ncycle) ! 1 is indicating that velverlet is called from hmc routine
   call pred_isothermal(ab_mover,hist,icycle,mttk_vars,ncycle,zDEBUG,iexit)

 !end if
 !END OF ATOMIC COORDINATES SWEEP************************************************

! else if(icycle>ncycle) then ! strain update
!   strain_updated = .TRUE.
!   strain_steps   = strain_steps + 1
!! Metropolis update of lattice vectors and parameters in case optcell/=0
!   if(icycle==ncycle+1.and.xred_updated) then
!     !save rprimd_hmc_prev and total electronic energy etotal_hmc_prev
!     call hist2var(acell_hmc_prev,hist,ab_mover%natom,rprimd_hmc_prev,xred,zDEBUG)
!     etotal_hmc_prev = hist%etot(hist%ihist)
!     strain_hmc_prev(:,:) = strain(:,:)
!
!     select case (ab_mover%optcell)
!     case (1) !volume optimization only
!       acell(:)=acell(:)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case (2,3,7,8,9) !full geometry optimization
!       !suggest new strain tensor values
!       do ii=1,3
!         do jj=ii,3
!           strain(ii,jj) = strain(ii,jj)+ 2.0_dp*dstrain*(uniformrandom(seed)-0.5_dp)
!           strain(jj,ii) = strain(ii,jj)
!         enddo
!       enddo
!       if(ab_mover%optcell==3) then !eliminate volume change if optcell==3
!         do ii=1,3
!           strain(ii,ii) = strain(ii,ii) -(strain(1,1)+strain(2,2)+strain(3,3))
!         enddo
!       endif
!       do jj=1,3    ! sum over three lattice vectors
!         do ii=1,3  ! sum over Cart components
!           rprimd(ii,jj)=rprimd_original(ii,jj)+&
!&                        rprimd_original(1,jj)*strain(ii,1)+&
!&                        rprimd_original(2,jj)*strain(ii,2)+&
!&                        rprimd_original(3,jj)*strain(ii,3)
!         enddo
!       enddo
!       if(ab_mover%optcell==7) then
!         rprimd(:,1)=rprimd_original(:,1)
!       else if (ab_mover%optcell==8) then
!         rprimd(:,2)=rprimd_original(:,2)
!       else if (ab_mover%optcell==9) then
!         rprimd(:,3)=rprimd_original(:,3)
!       endif
!     case(4)
!       acell(1)=acell(1)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case(5)
!       acell(2)=acell(2)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case(6)
!       acell(3)=acell(3)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     !case default
!     !  write(message,"(a,i0)") "Wrong value of optcell: ",ab_mover%optcell
!     !  MSG_ERROR(message)
!     end select
!
!     !update the new suggested rprimd and or acell in the history record
!     hist%ihist=abihist_findIndex(hist,+1)
!     call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
!   else
!
!     etotal = hist%etot(hist%ihist)
!     de = etotal - etotal_hmc_prev
!
!     iacc=0
!     rnd=uniformrandom(seed)
!     if(de<0)then
!       iacc=1
!     else
!       if(exp(-de/kbtemp)>rnd)then
!         iacc=1
!       end if
!     end if
!
!     if(iacc==0) then
!      strain(:,:)=strain_hmc_prev(:,:)
!      acell(:)=acell_hmc_prev(:)
!     else
!      call hist2var(acell_hmc_prev,hist,ab_mover%natom,rprimd_hmc_prev,xred,zDEBUG)
!      strain_hmc_prev(:,:) = strain(:,:)
!      etotal_hmc_prev=etotal
!     endif
!
!    !suggest new acell/rprimd values depending on the optcell value
!     select case (ab_mover%optcell)
!     case (1) !volume optimization only
!       acell(:)=acell(:)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case (2,3,7,8,9) !full geometry optimization
!       !suggest new strain tensor values
!       do ii=1,3
!         do jj=ii,3
!           strain(ii,jj) = strain(ii,jj)+ 2.0_dp*dstrain*(uniformrandom(seed)-0.5_dp)
!           strain(jj,ii) = strain(ii,jj)
!         enddo
!       enddo
!       if(ab_mover%optcell==3) then !eliminate volume change if optcell==3
!         do ii=1,3
!           strain(ii,ii) = strain(ii,ii) -(strain(1,1)+strain(2,2)+strain(3,3))
!         enddo
!       endif
!       do jj=1,3    ! sum over three lattice vectors
!         do ii=1,3  ! sum over Cart components
!           rprimd(ii,jj)=rprimd_original(ii,jj)+&
!&                        rprimd_original(1,jj)*strain(ii,1)+&
!&                        rprimd_original(2,jj)*strain(ii,2)+&
!&                        rprimd_original(3,jj)*strain(ii,3)
!         enddo
!       enddo
!       if(ab_mover%optcell==7) then
!         rprimd(:,1)=rprimd_original(:,1)
!       else if (ab_mover%optcell==8) then
!         rprimd(:,2)=rprimd_original(:,2)
!       else if (ab_mover%optcell==9) then
!         rprimd(:,3)=rprimd_original(:,3)
!       endif
!     case(4)
!       acell(1)=acell(1)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case(5)
!       acell(2)=acell(2)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case(6)
!       acell(3)=acell(3)*(1.0_dp+dstrain*2.0_dp*(uniformrandom(seed)-0.5_dp))
!     case default
!     !  write(message,"(a,i0)") "Wrong value of optcell: ",ab_mover%optcell
!     !  MSG_ERROR(message)
!     end select
!
!     !update the new suggested rprimd/acell in the history record
!     hist%ihist=abihist_findIndex(hist,+1)
!     call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
!
!   endif

 end if

end subroutine pred_hmc
!!***

end module m_pred_hmc
!!***

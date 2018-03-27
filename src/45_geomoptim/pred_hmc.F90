!{\src2tex{textfont=tt}}
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
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2018 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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
!!  hist = ionic positions, lattice parameters etc. are updated
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      hist2var,pred_velverlet,var2hist,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pred_hmc(ab_mover,hist,itime,icycle,ntime,ncycle,zDEBUG,iexit)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist

 use m_geometry,  only : xred2xcart
 use m_numeric_tools,  only : uniformrandom

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_hmc'
 use interfaces_45_geomoptim, except_this_one => pred_hmc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in)          :: itime
 integer,intent(in)          :: icycle
 integer,intent(in)          :: ntime
 integer,intent(in)          :: ncycle
 integer,intent(in)          :: iexit
 logical,intent(in)          :: zDEBUG

!Local variables-------------------------------

 integer       :: seed                                                         ! seed for rnd generator
 integer       :: ii,jj,iacc                                                   ! dummy integers for loop indexes and acceptance decision flag
 real(dp)      :: etotal,epot,ekin,de                                          ! total, potential (electronic), kinetic (ionic) energies and energy difference between initial and proposed states
 real(dp)      :: mv2tot,factor                                                ! dummies used for rescaling of velocities
 real(dp)      :: rnd
 real(dp)      :: xred(3,ab_mover%natom)                                       ! reduced coordinates of all ions
 real(dp)      :: vel(3,ab_mover%natom)                                        ! ionic velocities in Cartesian coordinates
 real(dp)      :: mvtot(3)                                                     ! total momentum of the cell used to rescale velocities
 real(dp)      :: mtot,kbtemp                                                  ! total ionic mass and target temperature in energy units
 real(dp)      :: acell(3)                                                     ! lattice parameters
 real(dp)      :: rprimd(3,3)                                                  ! lattice vectors

 real(dp),save :: etotal_hmc_prev                                              ! total energy of the initial state
 real(dp),save :: acell_hmc_prev(3)                                            !
 real(dp),save :: rprimd_hmc_prev(3,3)                                         !
 real(dp),allocatable,save :: xcart_hmc_prev(:,:)                              ! Cart. coordinates of the ions corresponding to the initial state
 real(dp),allocatable,save :: xred_hmc_prev(:,:)                               ! reduced coordinates of the ions corresponding to the initial state


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


 if(iexit/=0)then
   if (allocated(xcart_hmc_prev))  then
     ABI_DEALLOCATE(xcart_hmc_prev)
   end if
   if (allocated(xred_hmc_prev))  then
     ABI_DEALLOCATE(xred_hmc_prev)
   end if
   return
 end if


 !get current values of ionic positions and cell geometry and set up the target temperature
 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 vel(:,:)   = hist%vel(:,:,hist%ihist)                ! velocities of all ions, not needed in reality
 epot       = hist%etot(hist%ihist)                   ! electronic sub-system energy, not needed
 ekin       = hist%ekin(hist%ihist)                   ! kinetic energy, not needed


 kbtemp=(ab_mover%mdtemp(1)*kb_HaK)  ! get target temperature in energy units

 if(icycle==1)then

  !write(239117,*) ' entering first cycle of HMC iteration ',itime

   if(itime==1)then
     seed=-239

     if (allocated(xcart_hmc_prev))  then
       ABI_DEALLOCATE(xcart_hmc_prev)
     end if
     if (allocated(xred_hmc_prev))  then
       ABI_DEALLOCATE(xred_hmc_prev)
     end if
     ABI_ALLOCATE(xcart_hmc_prev,(3,ab_mover%natom))
     ABI_ALLOCATE(xred_hmc_prev,(3,ab_mover%natom))
   end if

 !generate new set of velocities and get rid of the possible overall momentum
   mtot=sum(ab_mover%amass(:))         ! total mass to eventually get rid of total center of mass (CoM) momentum

  ! generate velocities from normal distribution with zero mean and correct standard deviation
   do ii=1,ab_mover%natom
     do jj=1,3
       vel(jj,ii)=sqrt(kbtemp/ab_mover%amass(ii))*cos(two_pi*uniformrandom(seed))
       vel(jj,ii)=vel(jj,ii)*sqrt(-2.0*log(uniformrandom(seed)))
     end do
   end do
  !since number of atoms is most probably not big enough to obtain overall zero CoM momentum, shift the velocities
  !and then renormalize
   mvtot=zero ! total momentum
   do ii=1,ab_mover%natom
     do jj=1,3
       mvtot(jj)=mvtot(jj)+ab_mover%amass(ii)*vel(jj,ii)
     end do
   end do
   do ii=1,ab_mover%natom
     do jj=1,3
       vel(jj,ii)=vel(jj,ii)-(mvtot(jj)/mtot)
     end do
   end do
  !now the total cell momentum is zero
   mv2tot=0.0
   do ii=1,ab_mover%natom
     do jj=1,3
       mv2tot=mv2tot+ab_mover%amass(ii)*vel(jj,ii)**2
     end do
   end do
   factor = mv2tot/(dble(3*ab_mover%natom))
   factor = sqrt(kbtemp/factor)
   vel(:,:)=vel(:,:)*factor
   hist%vel(:,:,hist%ihist)=vel(:,:)



  !save the starting values of ionic positions (before an update attempt is performed)
   call hist2var(acell_hmc_prev,hist,ab_mover%natom,rprimd_hmc_prev,xred_hmc_prev,zDEBUG)
   call xred2xcart(ab_mover%natom,rprimd_hmc_prev,xcart_hmc_prev,xred_hmc_prev)

  !also save the initial total energy
   ekin=0.0
   do ii=1,ab_mover%natom
     do jj=1,3
       ekin=ekin+half*ab_mover%amass(ii)*vel(jj,ii)**2
     end do
   end do
   epot = hist%etot(hist%ihist)                     ! electronic sub-system energy
   etotal_hmc_prev = epot + ekin                    ! total energy before an attempt of variables update is performed


   call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,1,icycle,ncycle)

  !call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)
  !write(std_out,*) '  xcart_after_update:'
  !write(std_out,*) '  ',xcart(1,:)
  !write(std_out,*) '  ',xcart(2,:)
  !write(std_out,*) '  ',xcart(3,:)



 else if(icycle<ncycle)then

   call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit,1,icycle,ncycle)
  !call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)
  !vel(:,:)   = hist%vel(:,:,hist%ihist)                ! velocities of all ions, not needed in reality

 else !icycle==ncycle

  !the only thing left to do, is to compute the difference of the total energies and decide whether to accept the new state
   vel(:,:)=hist%vel(:,:,hist%ihist)
   ekin=0.0
   do ii=1,ab_mover%natom
     do jj=1,3
       ekin=ekin+half*ab_mover%amass(ii)*vel(jj,ii)**2
     end do
   end do
   epot = hist%etot(hist%ihist)      ! electronic sub-system energy
   etotal = epot + ekin
   de = etotal - etotal_hmc_prev

   iacc=0
   rnd=uniformrandom(seed)
   if(de<0)then
     iacc=1
   else
     if(exp(-de/kbtemp)>rnd)then
       iacc=1
     end if
   end if
  !write(238,*) '  random number: ',rnd,' -de/kbtemp:',-de/kbtemp,' acceptance decision: ',iacc
  !write(239,*) '  de: ',de,' estart: ',etotal_hmc_prev,' efin:', etotal
  !write(239,*) 'ekin= ',ekin,'  epot= ',epot,'  e_start= ',etotal_hmc_prev,'  e_end= ',etotal,'  iacc=',iacc, '  dekT=',-de/kbtemp

   call hist2var(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)
  !call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

   if(iacc==0)then  !in case the new state is not accepted, then roll back the coordinates
    !write(std_out,*) '  the proposed state was not accepted, roll back the configuration'
    !xcart(:,:)=xcart_hmc_prev(:,:)
    !no need to roll back xcart any longer, everything is stored in xred now (v8.3.1)
     xred(:,:)=xred_hmc_prev(:,:)
   end if

   hist%ihist=abihist_findIndex(hist,+1)
   call var2hist(acell,hist,ab_mover%natom,rprimd,xred,zDEBUG)

 end if

! note: add update of lattice vectors and parameters in case optcell/=0

end subroutine pred_hmc
!!***

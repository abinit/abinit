!{\src2tex{textfont=tt}}
!!****f* ABINIT/pred_hmc
!! NAME
!!  pred_hmc
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
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

subroutine pred_hmc(ab_mover,hist,itime,icycle,ntime,ncycle,zDEBUG,iexit)
    
 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_hmc'
 use interfaces_28_numeric_noabirule
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
 real(dp)      :: xcart(3,ab_mover%natom)                                      ! Cartesian coordinates of all ions
 real(dp)      :: xred(3,ab_mover%natom)                                       ! reduced coordinates of all ions
 real(dp)      :: vel(3,ab_mover%natom)                                        ! ionic velocities in Cartesian coordinates
 real(dp)      :: mvtot(3)                                                     ! total momentum of the cell used to rescale velocities
 real(dp)      :: mtot,kbtemp                                                  ! total ionic mass and target temperature in energy units
 real(dp)      :: acell(3)                                                     ! lattice parameters
 real(dp)      :: rprimd(3,3),rprim(3,3)                                       ! lattice vectors
 
 real(dp),save :: etotal_hmc_prev                                              ! total energy of the initial state
 real(dp),save :: acell_hmc_prev(3)                                            !
 real(dp),save :: rprimd_hmc_prev(3,3),rprim_hmc_prev(3,3)                     !
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


 write(std_out,*) ' entering HMC '


if(iexit/=0)then
   if (allocated(xcart_hmc_prev))  then
     ABI_DEALLOCATE(xcart_hmc_prev)
   end if
   if (allocated(xred_hmc_prev))  then
     ABI_DEALLOCATE(xred_hmc_prev)
   end if
   return
end if


 kbtemp=(ab_mover%mdtemp(1)*kb_HaK)  ! temperature in energy units

 write(std_out,*) ' HMC temperature: ', ab_mover%mdtemp(1),' K'


 if(icycle==1)then

  write(239117,*) ' entering first cycle of HMC iteration ',itime
  write(239117,*) '  target temperature: ',ab_mover%mdtemp(1)

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
    write(239117,*) '  allocated xcart_hmc_prev and xred_hmc_prev'
  endif

  write(239117,*) '  rnd seed: ', seed


 !generate new set of velocities and get rid of the possible overall momentum
  mtot=sum(ab_mover%amass(:))         ! total mass to eventually get rid of total center of mass (CoM) momentum
  write(239117,*) '  generating initial velocities'

  ! generate velocities from normal distribution with zero mean and correct standard deviation
  do ii=1,ab_mover%natom
    do jj=1,3
      vel(jj,ii)=sqrt(kbtemp/ab_mover%amass(ii))*cos(two_pi*uniformrandom(seed))
      vel(jj,ii)=vel(jj,ii)*sqrt(-2.0*log(uniformrandom(seed)))
    end do
  end do
  write(239117,*) '  vel:'
  write(239117,*) '  ',vel(1,:)
  write(239117,*) '  ',vel(2,:)
  write(239117,*) '  ',vel(3,:)
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
  hist%histV(:,:,hist%ihist)=vel(:,:)
  write(239117,*) '  factor: ',factor
  write(239117,*) '  vel_renorm:'
  write(239117,*) '  ',vel(1,:)
  write(239117,*) '  ',vel(2,:)
  write(239117,*) '  ',vel(3,:)

  write(239117,*) '  xcart:'
  write(239117,*) '  ',xcart(1,:)
  write(239117,*) '  ',xcart(2,:)
  write(239117,*) '  ',xcart(3,:)


  !now the starting values of ionic positions (before an update attempt is done)
  call hist2var(acell_hmc_prev,hist,ab_mover%natom,rprim_hmc_prev,rprimd_hmc_prev,xcart_hmc_prev,xred_hmc_prev,zDEBUG)

  !also save the initial total energy
  ekin=0.0
  do ii=1,ab_mover%natom
     do jj=1,3
       ekin=ekin+half*ab_mover%amass(ii)*vel(jj,ii)**2
     end do
  end do
  epot = hist%histE(hist%ihist)                    ! electronic sub-system energy
  etotal_hmc_prev = epot + ekin                    ! total energy before an attempt of variables update is performed

  write(239117,*) '  ekin,temp: ', ekin,kbtemp


  write(239117,*) '  calling velverlet predictor'
  !now everything is ready to make a prediction
  call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit)
 
 else if(icycle<ncycle)then

  write(239117,*) ' entering ',icycle,'nt cycle of HMC iteration ',itime
  call pred_velverlet(ab_mover,hist,itime,ntime,zDEBUG,iexit)

 else !presumably icycle==ncycle 

  write(239117,*) ' entering last cycle of HMC iteration ',itime
  !the only thing left to do, is to compute the difference of the total energies and decide whether to accept the new state
  vel(:,:)=hist%histV(:,:,hist%ihist)
  ekin=0.0
  do ii=1,ab_mover%natom
     do jj=1,3
       ekin=ekin+half*ab_mover%amass(ii)*vel(jj,ii)**2
     end do
  end do
  epot = hist%histE(hist%ihist)      ! electronic sub-system energy
  etotal = epot + ekin
  de = etotal - etotal_hmc_prev

  iacc=0
  rnd=uniformrandom(seed)
  if(de<0)then
    iacc=1
  else
    if(exp(-de/kbtemp)>rnd)then
     iacc=1
    endif
  endif
  write(239117,*) '  random number: ',rnd,' acceptance decision: ',iacc

  call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)
  if(iacc==0)then  !in case the new state is not accepted, then roll back the coordinates
    write(239117,*) '  the proposed state was not accepted, roll back the configuration'
    xcart(:,:)=xcart_hmc_prev(:,:)
    xred(:,:)=xred_hmc_prev(:,:)
  endif
  write(239117,*) '  xcart:'
  write(239117,*) '  ',xcart(1,:)
  write(239117,*) '  ',xcart(2,:)
  write(239117,*) '  ',xcart(3,:)

  hist%ihist=hist%ihist+1
  call var2hist(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,zDEBUG)

 endif

 !save current values of energies, as well as xred,xcart,etc.

! elseif(icycle==12)then

 !get the new values of energies

 !check if the new state is to be accepted

 !if the state is accepted, update the history and coordinates, otherwise retract

! elseif(icycle>12)then

! endif

end subroutine pred_hmc
!!***

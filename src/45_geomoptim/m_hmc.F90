!!****m* ABINIT/m_hmc
!! NAME
!!  m_hmc
!!
!! FUNCTION
!!  Auxiliary hmc functions
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2025 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_hmc

 use defs_basis
 use m_abicore
 use m_errors
 use m_abimover
 use m_io_tools

!use m_geometry,       only : xred2xcart
 use m_numeric_tools,  only : uniformrandom

 implicit none

 private

! *************************************************************************
 public :: compute_kinetic_energy
 public :: generate_random_velocities
 public :: metropolis_check

contains
!!***

!!****f* ABINIT/m_hmc/compute_kinetic_energy
!! NAME
!!  comute_kintic_energy
!!
!! FUNCTION
!!  Computes kintetic energy
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine compute_kinetic_energy(ab_mover,vel,ekin)

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 real(dp),      intent(in)   :: vel(3,ab_mover%natom) ! velocities
 real(dp),      intent(out)  :: ekin                  ! output kinetic energy

!Local variables-------------------------------
 integer :: ii
!character(len=500) :: msg

! *************************************************************************

 ekin=0.0
do ii = 1, ab_mover%natom
  ekin = ekin + half * ab_mover%amass(ii) * DOT_PRODUCT(vel(:, ii), vel(:, ii))
end do

end subroutine compute_kinetic_energy
!!***






!!****f* ABINIT/m_hmc/generate_random_velocities
!! NAME
!!  generate_random_velocities
!!
!! FUNCTION
!!  Generate normally distributed random velocities
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine generate_random_velocities(ab_mover,kbtemp,seed,vel,ekin)

!Arguments ------------------------------------
 type(abimover),intent(in)   :: ab_mover
 integer,       intent(inout):: seed
 real(dp),      intent(in)   :: kbtemp
 real(dp),      intent(inout):: vel(3,ab_mover%natom) ! velocities
 real(dp),      intent(out)  :: ekin                  ! output kinetic energy
!Local variables-------------------------------
 integer :: ii,jj,natom
 real(dp):: mtot,mvtot(3),mv2tot,factor
!character(len=500) :: msg

! *************************************************************************


 natom = ab_mover%natom
 mtot=sum(ab_mover%amass(:))         ! total mass to eventually get rid of total center of mass (CoM) momentum
 !generate velocities from normal distribution with zero mean and correct standard deviation
 do ii=1,ab_mover%natom
   do jj=1,3
     vel(jj,ii)=sqrt(kbtemp/ab_mover%amass(ii))*cos(two_pi*uniformrandom(seed))
     vel(jj,ii)=vel(jj,ii)*sqrt(-2.0*log(uniformrandom(seed)))
   end do
 end do
 !since number of atoms is most probably not big enough to obtain overall zero CoM momentum, shift the velocities
 !and then renormalize
 ! mvtot -> total momentum
 mvtot(:) = MATMUL(vel(1:3,1:natom), ab_mover%amass(1:natom))
 do ii=1,ab_mover%natom
   vel(:,ii)=vel(1:3,ii)-(mvtot(1:3)/mtot)
 end do
 !now the total cell momentum is zero
 mv2tot=0.0
 do ii=1,ab_mover%natom
   mv2tot=mv2tot+ab_mover%amass(ii)* DOT_PRODUCT(vel(:,ii), vel(:,ii))
 end do
 factor = mv2tot/(dble(3*ab_mover%natom))
 factor = sqrt(kbtemp/factor)
 vel(:,:)=vel(:,:)*factor

 call compute_kinetic_energy(ab_mover,vel,ekin)

end subroutine generate_random_velocities
!!***




!!****f* ABINIT/m_hmc/metropolis_check
!! NAME
!!  metropolis_check
!!
!! FUNCTION
!!  Make an acceptance decision based on the energy differences
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

subroutine metropolis_check(seed,de,kbtemp,iacc)

!Arguments ------------------------------------
 integer,       intent(inout):: seed
 real(dp),      intent(in)   :: de
 real(dp),      intent(in)   :: kbtemp
 integer,       intent(inout):: iacc

!Local variables-------------------------------
 real(dp)   :: rnd
!character(len=500) :: msg

! *************************************************************************

 iacc=0
 rnd=uniformrandom(seed)
 if(de<0)then
   iacc=1
 else
   if(exp(-de/kbtemp)>rnd)then
      iacc=1
   end if
 end if

end subroutine metropolis_check
!!***




end module m_hmc
!!***

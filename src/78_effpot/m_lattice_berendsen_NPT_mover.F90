!!****m* ABINIT/m_lattice_berendsen_NPT_mover
!! TODO: This is not yet implemented.
!! NAME
!! m_lattice_berendsen_NPT_mover
!!
!! FUNCTION
!! This module contains the berendsen  (NPT) lattice mover.
!! The method is described in  
!! H.J.C. Berendsen, J.P.M. Postma, A. DiNola, and J.R. Haak,
!! "Molecular dynamics with coupling to an external bath,"
!!  J. Chem. Phys., 81 3684-3690 (1984)
!! NOTE: that this method does NOT generate properly the thermostated
!! ensemble. It does not have the correct distribution of the kinetic energy.
!! However, it approches the target temperature exponentially without oscillation, 
!! for which the steps can be easily controlled.
!!
!! Datatypes:
!!
!! * lattice_berendsen_NPT_mover_t: defines the lattice movers
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE



#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lattice_berendsen_NPT_mover
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_mover, only: lattice_mover_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_random_xoroshiro128plus, only:  rng_t
  use m_hashtable_strval, only: hash_table_t
!!***

  implicit none

  private

  type, public, extends(lattice_mover_t) :: lattice_berendsen_NPT_mover_t
     real(dp) :: taut ! the characteristic time of the relaxation of velocity.
     ! it is usually larger than the time step.
     real(dp) :: taup !  the characteristic time of the relaxation of pressure.
     real(dp) :: compressibility
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_one_step
     procedure :: scale_velocities
  end type lattice_berendsen_NPT_mover_t

contains



  subroutine initialize(self,params, supercell, rng)
    class(lattice_berendsen_NPT_mover_t), intent(inout) :: self
    type(multibinit_dtset_type), target, intent(in):: params
    type(mbsupercell_t), target, intent(in) :: supercell
    type(rng_t), target, intent(in) :: rng
    self%taut = params%latt_taut
    self%taup = params%latt_taup
    self%compressibility =params%latt_compressibility
    call self%lattice_mover_t%initialize(params, supercell, rng)
    MSG_ERROR("The Berendsen NPT mover has not yet been implemented")
    !TODO: Implement
  end subroutine initialize


  subroutine finalize(self)
    class(lattice_berendsen_NPT_mover_t), intent(inout) :: self
    call self%lattice_mover_t%finalize()
  end subroutine finalize


  !-------------------------------------------------------------------!
  ! scale_velocities:
  !   scale the velocities so that they get close to the required temperture
  !
  !-------------------------------------------------------------------!
  subroutine scale_velocities(self)
    class(lattice_berendsen_NPT_mover_t), intent(inout) :: self
    real(dp) :: tautscl, old_temperature, scale_temperature, tmp
    tautscl = self%dt / self%taut
    old_temperature=self%T_ob
    tmp=1.0 +(self%temperature / old_temperature - 1.0) *    tautscl
    if(tmp< 0.0) then
       MSG_ERROR("The time scale for the Berendsen Algorithm should be at least larger than dtion.")
    else
       scale_temperature=sqrt(tmp)
    end if
    ! Limit the velocity scaling to reasonable values
    if( scale_temperature > 1.1) then
       scale_temperature = 1.1
    elseif (scale_temperature < 0.9) then
       scale_temperature = 0.9
    endif
    self%current_vcart(:,:) = self%current_vcart(:,:) * scale_temperature
  end subroutine scale_velocities
 

  !-------------------------------------------------------------------!
  ! run_one_step.
  ! The algorithm is almost the same as the velocity verlet algorithm,
  ! except at the begining, the velocities are scaled so that the temperature
  ! is getting closer to the required temperature.
  !-------------------------------------------------------------------!
  subroutine run_one_step(self, effpot,displacement, strain, spin, lwf, energy_table)
    class(lattice_berendsen_NPT_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    integer :: i
    character(len=40) :: key

    ABI_UNUSED(displacement)
    ABI_UNUSED(strain)

    ! scale the velocity.
    call self%scale_velocities()

    self%energy=0.0
    self%forces(:,:) =0.0
    call effpot%calculate( displacement=self%displacement, strain=self%strain, &
         & spin=spin, lwf=lwf, force=self%forces, stress=self%stress, &
         & energy=self%energy, energy_table=energy_table)
    do i=1, self%natom
       self%current_vcart(:,i) = self%current_vcart(:,i) + &
            & (0.5_dp * self%dt) * self%forces(:,i)/self%masses(i)
    end do
    call self%force_stationary()
    self%displacement(:,:) = self%displacement(:,:)+self%current_vcart(:,:) * self%dt


    ! second half of velocity update.
    ! v(t+dt) = v(t + 1/2 dt) + F/m * 1/2 dt
    ! NOTE: energy and forces should be initialized before every calculation!
    self%energy=0.0
    self%forces(:,:)=0.0
    call effpot%calculate( displacement=self%displacement, &
         & strain=self%strain, spin=spin, lwf=lwf, force=self%forces, &
         & stress=self%stress,  energy=self%energy, energy_table=energy_table)
    do i=1, self%natom
       self%current_vcart(:,i) = self%current_vcart(:,i) &
            & + (0.5_dp * self%dt) * self%forces(:,i)/self%masses(i)
    end do
    call self%force_stationary()

    call self%get_T_and_Ek()
    if (present(energy_table)) then
      key = 'Lattice kinetic energy'
      call energy_table%put(key, self%Ek)
    end if


  end subroutine run_one_step

end module m_lattice_berendsen_NPT_mover


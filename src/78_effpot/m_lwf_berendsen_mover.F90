!!****m* ABINIT/m_lwf_berendsen_mover
!! NAME
!! m_lwf_berendsen_mover
!!
!! FUNCTION
!! This module contains the lwf berensen NVT mover 
!!
!!
!! Datatypes:
!!
!! * lwf_berendsen_mover_t
!!
!! Subroutines:
!!
!! * TODO: update this when F2003 documentation format decided.
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
module m_lwf_berendsen_mover
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_nctk
#define HAVE_NETCDF 1
#if defined HAVE_NETCDF
  use netcdf
#endif
  use m_mpi_scheduler, only: mpi_scheduler_t, init_mpi_info
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_hashtable_strval, only: hash_table_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_lwf_hist, only: lwf_hist_t
  use m_lwf_observables, only: lwf_observables_t
  use m_lwf_ncfile, only: lwf_ncfile_t
  use m_lwf_mover, only: lwf_mover_t

  implicit none
  private
  !!***

  type, public, extends(lwf_mover_t) :: lwf_berendsen_mover_t
     real(dp) :: taut ! the characteristic time of the relaxation of velocity.
   contains
     procedure :: set_params
     procedure :: scale_velocities
     procedure :: run_one_step
  end type lwf_berendsen_mover_t

  contains

    subroutine set_params(self, params)
      class(lwf_berendsen_mover_t), intent(inout) :: self
      type(multibinit_dtset_type) :: params
      call self%lwf_mover_t%set_params(params)
      self%taut=params%lwf_taut
    end subroutine set_params


    subroutine run_one_step(self, effpot, displacement, strain, spin, lwf, energy_table)
      class(lwf_berendsen_mover_t), intent(inout) :: self
      real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
      class(abstract_potential_t), intent(inout) :: effpot
      type(hash_table_t),optional, intent(inout) :: energy_table
      integer :: i
      character(len=40) :: key
      ABI_UNUSED_A(lwf)
      ! scale the velocity.
      call self%scale_velocities()
      self%energy=0.0
      self%lwf_force(:) =0.0
      call effpot%calculate( displacement=displacement, strain=strain, &
           & spin=spin, lwf=self%lwf, lwf_force=self%lwf_force, &
           & energy=self%energy, energy_table=energy_table)
      !print *, "lwf_force", self%lwf_force
      !print *, "lwf_masses", self%lwf_masses
      do i=1, self%nlwf
         self%vcart(i) = self%vcart(i) + &
              & (0.5_dp * self%dt) * self%lwf_force(i)/self%lwf_masses(i)
      end do
      !print *, "vcart", self%vcart
      !call self%force_stationary()
      self%lwf= self%lwf+self%vcart * self%dt
      !print *, "lwf", self%lwf


      self%energy=0.0
      self%lwf_force(:)=0.0
      call effpot%calculate( displacement=displacement, strain=strain, &
           & spin=spin, lwf=self%lwf, lwf_force=self%lwf_force, &
           & energy=self%energy, energy_table=energy_table)
      call effpot%calculate( displacement=displacement, strain=strain, &
           & spin=spin, lwf=self%lwf, lwf_force=self%lwf_force, &
           & energy=self%energy, energy_table=energy_table)
      do i=1, self%nlwf
         self%vcart(i) = self%vcart(i) + &
              & (0.5_dp * self%dt) * self%lwf_force(i)/self%lwf_masses(i)
      end do
      !call self%force_stationary()
      self%lwf= self%lwf+self%vcart * self%dt
      call self%get_T_and_Ek()

      if (present(energy_table)) then
         key = 'Lwf kinetic energy'
         call energy_table%put(key, self%Ek)
      end if

    end subroutine run_one_step


  !-------------------------------------------------------------------!
  ! scale_velocities:
  !   scale the velocities so that they get close to the required temperture
  !
  !-------------------------------------------------------------------!
  subroutine scale_velocities(self)
    class(lwf_berendsen_mover_t), intent(inout) :: self
    real(dp) :: tautscl, old_temperature, scale_temperature, tmp
    tautscl = self%dt / self%taut
    old_temperature=self%T_ob
    if (old_temperature< 1e-19) then
       old_temperature=1e-19
    end if
    tmp=1.0 +(self%temperature / old_temperature - 1.0) *    tautscl
    if(tmp< 0.0) then
       MSG_ERROR("The time scale for the Berendsen algorithm should be at least larger than dtion")
    else
       scale_temperature=sqrt(tmp)
    end if
    ! Limit the velocity scaling to reasonable values
    if( scale_temperature > 1.1) then
       scale_temperature = 1.1
    elseif (scale_temperature < 0.9) then
       scale_temperature = 0.9
    endif
    self%vcart = self%vcart * scale_temperature
  end subroutine scale_velocities


end module m_lwf_berendsen_mover



!!****m* ABINIT/m_lwf_mc_mover
!! NAME
!! m_lwf_mc_mover
!!
!! FUNCTION
!! This module contains the lwf Markov chain Monte Carlo functions for lwf mover .
!!
!!
!! Datatypes:
!!
!! * lwf_mc_t : MCMC. It defines how to move lwfs in one step,
!! attempt function: whether to accept move
!! accecpt/reject method which define what to do if move is
!! accepted or rejected!! . 
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
    module m_lwf_mc_mover
    use defs_basis
    use m_abicore
    use m_errors
    use m_abstract_potential, only: abstract_potential_t
    use m_random_xoroshiro128plus, only: rng_t
    use m_hashtable_strval, only: hash_table_t
    use m_multibinit_dataset, only: multibinit_dtset_type
    use m_lwf_mover, only: lwf_mover_t
    use m_multibinit_cell, only: mbcell_t, mbsupercell_t
    implicit none
!!***
    private

    !----------------------------------------------------------------------
    !> @brief An helper type to run lwf dynamics
    ! The Metropolis-Hasting algorithm is used.
    !----------------------------------------------------------------------
    type,public, extends(lwf_mover_t) :: lwf_mc_t
       real(dp) :: avg_amp! an angle to rotate by average
       real(dp) ::  deltaE ! energy and the change of energy when change one lwf
       real(dp) :: lwf_new, lwf_old! energy and the change of energy when change one lwf
       real(dp) :: beta  ! 1/(kb T)
       integer :: nstep  ! number of steps
       integer :: imove ! index of lwf to be moved
       integer :: naccept  ! number of accepted steps
     contains
       procedure :: initialize
       procedure :: finalize
       procedure, private :: attempt
       procedure, private :: accept
       procedure, private :: reject
       procedure  :: run_one_step
       procedure :: set_temperature
       procedure, private :: run_one_mc_step
    end type lwf_mc_t


  contains

    subroutine initialize(self, params, supercell, rng)
      class(lwf_mc_t), intent(inout) :: self
      type(multibinit_dtset_type),target, intent(in) :: params
      type(mbsupercell_t),target, intent(in) :: supercell
      type(rng_t), target, intent(in) :: rng
      call self%lwf_mover_t%initialize(params, supercell, rng)
      self%nstep=self%nlwf
      self%avg_amp=params%lwf_mc_avg_amp
      self%temperature=params%lwf_temperature
      self%beta=1.0/self%temperature ! Kb in a.u. is 1.
      self%lwf_new=0.0_dp
      self%lwf_old=0.0_dp
    end subroutine initialize

    !----------------------------------------------------------------------
    !> @brief finalize
    !----------------------------------------------------------------------
    subroutine finalize(self)
      class(lwf_mc_t), intent(inout) :: self
      call self%lwf_mover_t%finalize()
      self%nstep=0
      self%avg_amp=0.0
      self%temperature=0.0
      self%beta=0.0
    end subroutine finalize

    subroutine set_temperature(self, temperature)
      class(lwf_mc_t), intent(inout) :: self
      real(dp), intent(in) :: temperature
      call self%lwf_mover_t%set_temperature(temperature)
      self%temperature=temperature
      self%beta=1.0/self%temperature ! Kb in a.u. is 1.
    end subroutine set_temperature


    !----------------------------------------------------------------------
    !> @brief run one monte carlo step
    !> @param[in]   rngL rundom number generator
    !> @param[in] effpot: effective lwf potential
    !----------------------------------------------------------------------
   subroutine run_one_mc_step(self,  effpot)
     class(lwf_mc_t) :: self
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp) :: r

     ! try to change lwf
     r=self%attempt(self%rng, effpot)
     ! metropolis-hastings
     if(self%rng%rand_unif_01()< min(1.0_dp, r) .and. abs(self%lwf_new)<0.5 ) then
        self%naccept=self%naccept+1
        call self%accept()
     else
        call self%reject()
     end if
   end subroutine run_one_mc_step


  !-------------------------------------------------------------------!
  !> @brief: Run_one_step
  !>   effpot: the potential (which do the calculation of E and dE/dvar)
  !> param[in]: effpot
  !> param[in]: (optional) displacement
  !> param[in]: (optional) strain
  !> param[in]: (optional) spin
  !> param[in]: (optional) lwf
  ! NOTE: No need to pass the variable already saved in the mover.
  !     e.g. For spin mover, do NOT pass the spin to it.
  !    The other variables are only required if there is coupling with
  !    the mover variable.
  !-------------------------------------------------------------------!
  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf,  energy_table)
    ! run one step. (For MC also?)
    class(lwf_mc_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    class(abstract_potential_t), intent(inout) :: effpot
    integer :: i
    character(len=40) :: key

    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)

    !self%lwf(:)=lwf(:)
    self%lwf_force(:) = 0.0_dp
    ! calculate energy at first step
    self%energy=0.0


    call effpot%calculate( displacement=displacement, strain=strain, &
         & spin=spin, lwf=self%lwf, lwf_force=self%lwf_force, &
         & energy=self%energy, energy_table=energy_table)
    !print *, "Calculated energy", self%energy/self%supercell%ncell
    do i = 1, self%nstep
       call self%run_one_mc_step(effpot)
    end do
    lwf(:)=self%lwf(:)
    if (present(energy_table)) then
       key = 'LWF energy'
       call energy_table%put(key, self%energy)
    end if

  end subroutine run_one_step


   !----------------------------------------------------------------------
   !> @brief accept the trail step, which update the lwf and energy
   !----------------------------------------------------------------------
   subroutine accept(self)
     class(lwf_mc_t), intent(inout) :: self
     self%lwf(self%imove)=self%lwf_new
     self%energy=self%energy+self%deltaE
     !print *, "E:", self%energy/self%supercell%ncell
   end subroutine accept

   !----------------------------------------------------------------------
   !> @brief reject the trail step, changes nothing.
   !----------------------------------------------------------------------
   subroutine reject(self)
     class(lwf_mc_t), intent(inout) :: self
     ! do nothing.
     ABI_UNUSED_A(self)
   end subroutine reject

   !----------------------------------------------------------------------
   !> @brief define a trail step  using Hinzke_nowak method and calculate energy difference
   !----------------------------------------------------------------------
   function attempt(self,rng, effpot) result(r)
     class(lwf_mc_t) :: self
     class(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp) :: r
     ! choose one site
     self%imove = rng%rand_choice(self%nlwf)
     self%lwf_old= self%lwf(self%imove)
     self%deltaE=0.0
     call move(rng, self%lwf_old, self%lwf_new, self%avg_amp)
     call effpot%get_delta_E_lwf( self%lwf, self%imove, self%lwf_new, self%deltaE)
     r=exp(-self%deltaE *self%beta)
   end function attempt

   !----------------------------------------------------------------------
   !> @brief  rotate the  lwf by the average of angle (normal distribution)
   !----------------------------------------------------------------------
   subroutine move(rng, lwf_old, lwf_new, avg_amp)
     type(rng_t) :: rng
     real(dp), intent(in) :: lwf_old, avg_amp
     real(dp), intent(out) :: lwf_new
     real(dp):: dlwf
     dlwf=rng%rand_normal()
     lwf_new=lwf_old + dlwf*avg_amp
   end subroutine move

end module m_lwf_mc_mover

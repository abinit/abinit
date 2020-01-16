!!****m* ABINIT/m_spin_mc_mover
!! NAME
!! m_spin_mc_mover
!!
!! FUNCTION
!! This module contains the spin Markov chain Monte Carlo functions for spin mover .
!!
!!
!! Datatypes:
!!
!! * spin_mc_t : MCMC. It defines how to move spins in one step,
!! attempt function: whether to accept move
!! accecpt/reject method which define what to do if move is
!! accepted or rejected!! . 
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
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
    module m_spin_mc_mover
    use defs_basis
    use m_abicore
    use m_errors
    use m_abstract_potential, only: abstract_potential_t
    use m_random_xoroshiro128plus, only: rng_t
    use m_hashtable_strval, only: hash_table_t
    implicit none
!!***
    private

    !----------------------------------------------------------------------
    !> @brief An helper type to run spin dynamics
    ! The Metropolis-Hasting algorithm is used.
    !----------------------------------------------------------------------
    type,public :: spin_mc_t
       real(dp), allocatable :: S(:,:) ! the spin for the whole structure
       real(dp) ::  Sold(3), Snew(3) ! old and new S for one spin
       real(dp) :: angle   ! an angle to rotate by average
       real(dp) :: energy, deltaE ! energy and the change of energy when change one spin
       real(dp) :: temperature  
       real(dp) :: beta  ! 1/(kb T)
       integer :: nspin   ! number of spins
       integer :: nstep  ! number of steps
       integer :: imove ! index of spin to be moved
       integer :: naccept  ! number of accepted steps
     contains
       procedure :: initialize
       procedure :: finalize
       procedure, private :: attempt
       procedure, private :: accept
       procedure, private :: reject
       procedure, private :: run_one_step
       procedure :: run_MC
    end type spin_mc_t


  contains
    !----------------------------------------------------------------------
    !> @brief initialize mc helper class
    !>
    !> @param[in]  nspin: number of spins
    !> @param[in]  angle: a angle to rotate
    !> @param[in]  temperature: temperature
    !----------------------------------------------------------------------
    subroutine initialize(self, nspin, angle, temperature)
      class(spin_mc_t), intent(inout) :: self
      integer, intent(in) :: nspin
      real(dp), intent(in) :: angle, temperature
      self%nspin=nspin
      self%nstep=self%nspin
      ABI_ALLOCATE(self%S, (3, self%nspin))
      self%angle=angle
      self%temperature=temperature
      self%beta=1.0/temperature ! Kb in a.u. is 1.
      self%Sold(:)=0.0_dp
      self%Snew(:)=0.0_dp
    end subroutine initialize

    !----------------------------------------------------------------------
    !> @brief finalize
    !----------------------------------------------------------------------
    subroutine finalize(self)
      class(spin_mc_t), intent(inout) :: self
      if (allocated(self%S)) then
         ABI_DEALLOCATE(self%S)
      end if
      self%Sold=zero
      self%Snew=zero
      self%nspin=0
      self%nstep=0
    end subroutine finalize


    !----------------------------------------------------------------------
    !> @brief run one monte carlo step
    !> @param[in]   rngL rundom number generator
    !> @param[in] effpot: effective spin potential
    !----------------------------------------------------------------------
   subroutine run_one_step(self, rng, effpot)
     class(spin_mc_t) :: self
     class(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp) :: r

     ! try to change spin
     r=self%attempt(rng, effpot)
     ! metropolis-hastings
     if(rng%rand_unif_01()< min(1.0_dp, r) ) then
        self%naccept=self%naccept+1
        call self%accept()
     else
        call self%reject()
     end if
   end subroutine run_one_step

   !----------------------------------------------------------------------
   !> @brief run a number of MC steps. Since one step only changes
   !> too little things, a few steps are bunched as one. Then things like
   !> output or calculation of observables are done after the big step.
   !>
   !> @param[in]  rng: random number generator
   !> @param[in]  effpot:  the spin potential
   !> @param[in]  S_in:  the intial spin state
   !> @param[out]  etot:  the final total energy
   !----------------------------------------------------------------------
   subroutine run_MC(self, rng, effpot, S_in, etot)
     class(spin_mc_t), intent(inout) :: self
     type(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp), intent(inout) :: S_in(3,self%nspin)
     real(dp), intent(out) ::  etot
     integer :: i
     self%S(:,:)=S_in(:,:)
     call effpot%calculate(spin=S_in, energy=self%energy)
     do i = 1, self%nstep
        call self%run_one_step(rng, effpot)
     end do
     S_in(:, :)=self%S(:,:)
     etot=self%energy
   end subroutine run_MC

   !----------------------------------------------------------------------
   !> @brief accept the trail step, which update the spin and energy
   !----------------------------------------------------------------------
   subroutine accept(self)
     class(spin_mc_t), intent(inout) :: self
     self%S(:,self%imove)=self%Snew(:)
     self%energy=self%energy+self%deltaE
   end subroutine accept

   !----------------------------------------------------------------------
   !> @brief reject the trail step, changes nothing.
   !----------------------------------------------------------------------
   subroutine reject(self)
     class(spin_mc_t), intent(inout) :: self
     ! do nothing.
     ABI_UNUSED_A(self)
   end subroutine reject

   !----------------------------------------------------------------------
   !> @brief define a trail step  using Hinzke_nowak method and calculate energy difference
   !----------------------------------------------------------------------
   function attempt(self,rng, effpot) result(r)
     class(spin_mc_t) :: self
     class(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp) :: r
     ! choose one site
     self%imove = rng%rand_choice(self%nspin)
     self%Sold(:)= self%S(:,self%imove)
     self%deltaE=0.0
     call move_hinzke_nowak(rng, self%Sold, self%Snew, self%angle)
     call effpot%get_delta_E( self%S, self%imove, self%Snew, self%deltaE)
     r=exp(-self%deltaE *self%beta)
   end function attempt

   !----------------------------------------------------------------------
   !> @brief  rotate the  spin by the average of angle (normal distribution)
   !----------------------------------------------------------------------
   subroutine move_angle(rng, Sold, Snew, angle)
     type(rng_t) :: rng
     real(dp), intent(in) :: Sold(3), angle
     real(dp), intent(out) :: Snew(3)
     call rng%rand_normal_array(Snew, 3)
     Snew(:)=Sold(:) + Snew(:)*angle
     Snew(:)=Snew(:)/sqrt(Snew(1)*Snew(1)+Snew(2)*Snew(2)+Snew(3)*Snew(3))
   end subroutine move_angle

   !----------------------------------------------------------------------
   !> @brief  flip one spin
   !----------------------------------------------------------------------
   subroutine move_flip(Sold, Snew)
     real(dp), intent(in) :: Sold(3)
     real(dp), intent(out) :: Snew(3)
     Snew(:)=-Sold(:)
   end subroutine move_flip

   !----------------------------------------------------------------------
   !> @brief  set spin to random orientation
   !----------------------------------------------------------------------
   subroutine move_uniform(rng, Snew)
     type(rng_t), intent(inout) :: rng
     real(dp), intent(out) :: Snew(3)
     call rng%rand_normal_array(Snew, 3)
     Snew(:)=Snew(:)/sqrt(Snew(1)*Snew(1)+Snew(2)*Snew(2)+Snew(3)*Snew(3))
   end subroutine move_uniform

   !----------------------------------------------------------------------
   !> @brief combine rotate, flip and random set.
   !----------------------------------------------------------------------
   subroutine move_hinzke_nowak(rng, Sold, Snew, angle)
     type(rng_t), intent(inout) :: rng
     real(dp), intent(in) :: Sold(3), angle
     real(dp), intent(out) :: Snew(3)
     integer :: move
     move=rng%rand_choice(3)
     select case (move)
     case (1)
        call move_angle(rng, Sold, Snew, angle)
     case(2)
        call move_flip(Sold, Snew)
     case(3)
        call move_uniform(rng, Snew)
     case default
        call move_angle(rng, Sold, Snew, angle)
     end select
   end subroutine move_hinzke_nowak

end module m_spin_mc_mover

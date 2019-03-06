!{\src2tex{textfont=tt}}
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
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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
    use m_abstract_potential, only: abstract_potential_t
    use m_random_xoroshiro128plus, only: rng_t
    implicit none
!!***
    private

    type,public :: spin_mc_t
       real(dp), allocatable :: S(:,:)
       real(dp) ::  Sold(3), Snew(3)
       real(dp) :: angle
       real(dp) :: energy, deltaE
       real(dp) :: temperature
       real(dp) :: beta  ! 1/(kb T)
       integer :: nspin
       integer :: nstep
       integer :: imove ! index of spin to be moved
       integer :: naccept
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


   subroutine run_one_step(self, rng, effpot)
     class(spin_mc_t) :: self
     class(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp) :: S_out(3,self%nspin), etot, r
     integer :: i, j
     r=self%attempt(rng, effpot)
     !print *, "r", r
     if(rng%rand_unif_01()< min(1.0, r) ) then
        self%naccept=self%naccept+1
        call self%accept()
     else
        call self%reject()
     end if
   end subroutine run_one_step

   subroutine run_MC(self, rng, effpot, S_in, S_out, etot)
     class(spin_mc_t), intent(inout) :: self
     type(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp), intent(inout) :: S_in(3,self%nspin)
     real(dp), intent(out) :: S_out(3,self%nspin), etot
     integer :: i
     self%S(:,:)=S_in(:,:)
     !print*, "S_in", S_in
     call effpot%calculate(spin=S_in, energy=self%energy)
     !print *, self%energy
     do i = 1, self%nstep
        call self%run_one_step(rng, effpot)
     end do
     S_out(:, :)=self%S(:,:)
     etot=self%energy
   end subroutine run_MC

   subroutine accept(self)
     class(spin_mc_t), intent(inout) :: self
     self%S(:,self%imove)=self%Snew(:)
     self%energy=self%energy+self%deltaE
   end subroutine accept

   subroutine reject(self)
     class(spin_mc_t), intent(inout) :: self
     ! do nothing.
   end subroutine reject

   function attempt(self,rng, effpot) result(r)
     class(spin_mc_t) :: self
     class(rng_t) :: rng
     class(abstract_potential_t), intent(inout) :: effpot
     real(dp) :: r
     ! choose one site
     self%imove = rng%rand_choice(self%nspin)
     self%Sold(:)= self%S(:,self%imove)
     call move_hinzke_nowak(rng, self%Sold, self%Snew, self%angle)
     call effpot%get_delta_E( self%S, self%imove, self%Snew, self%deltaE)
     r=exp(-self%deltaE *self%beta)
   end function attempt

   subroutine move_angle(rng, Sold, Snew, angle)
     type(rng_t) :: rng
     real(dp), intent(in) :: Sold(3), angle
     real(dp), intent(out) :: Snew(3)
     real(dp) :: m
     call rng%rand_normal_array(Snew, 3)
     Snew(:)=Sold(:) + Snew(:)*angle
     Snew(:)=Snew(:)/sqrt(Snew(1)*Snew(1)+Snew(2)*Snew(2)+Snew(3)*Snew(3))
   end subroutine move_angle

   subroutine move_flip(Sold, Snew)
     real(dp), intent(in) :: Sold(3)
     real(dp), intent(out) :: Snew(3)
     Snew(:)=-Sold(:)
   end subroutine move_flip

   subroutine move_uniform(rng, Snew)
     type(rng_t), intent(inout) :: rng
     real(dp), intent(out) :: Snew(3)
     call rng%rand_normal_array(Snew, 3)
     Snew(:)=Snew(:)/sqrt(Snew(1)*Snew(1)+Snew(2)*Snew(2)+Snew(3)*Snew(3))
   end subroutine move_uniform

   subroutine move_hinzke_nowak(rng, Sold, Snew, angle)
     type(rng_t), intent(inout) :: rng
     real(dp), intent(in) :: Sold(3), angle
     real(dp), intent(out) :: Snew(3)
     integer :: move
     move=rng%rand_choice(3)
     !print *, "choice of move", move
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

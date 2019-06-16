!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_mover
!! NAME
!! m_lattice_mover
!!
!! FUNCTION
!! This module contains the lattice mover.
!!
!!
!! Datatypes:
!!
!! * lattice_mover_t: defines the lattice movers
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

module m_lattice_mover
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_random_xoroshiro128plus, only:  rng_t

!!***

  implicit none
  private

  type ,public, extends(abstract_mover_t) :: lattice_mover_t
     ! This is the abstract class of mover
     ! It do the following things:
     ! calculate d(var)/dt and integrate new var. 
     ! call functions to calculate observables.
     ! interact with hist file.
     type(multibinit_dtset_type), pointer :: params=>null()
     real(dp), allocatable :: masses(:)
     integer :: natom
     real(dp) :: stress(3,3), strain(3,3)
     real(dp), allocatable :: current_xcart(:,:), current_vcart(:,:), &
          & forces(:,:), displacement(:,:)
     real(dp) :: energy
     !type(lattice_hist_t) :: hist
   contains
     procedure:: initialize       ! perhaps each effpot type should have own 
     procedure :: finalize
     procedure :: set_params
     procedure :: set_initial_state ! initial state
     procedure:: run_one_step 
     procedure :: reset            ! reset the mover
     procedure :: calc_observables ! call functions to calculate observables
     procedure :: write_hist       ! write hist file
  end type lattice_mover_t


contains

  subroutine initialize(self, params, supercell)
    class(lattice_mover_t), intent(inout) :: self
    type(multibinit_dtset_type),target, intent(in) :: params
    type(mbsupercell_t),target, intent(in) :: supercell
    self%params=>params
    self%supercell=>supercell
    self%label="Lattice Mover"
    !TODO: self%natom = supercell%natom
    ABI_ALLOCATE(self%masses, (self%natom))
    ABI_ALLOCATE(self%displacement, (3, self%natom))
    ABI_ALLOCATE(self%current_xcart, (3, self%natom))
    ABI_ALLOCATE(self%current_vcart, (3, self%natom))
    ABI_ALLOCATE(self%forces, (3,self%natom))

  end subroutine initialize

  subroutine finalize(self)
    class(lattice_mover_t), intent(inout) :: self
    ABI_DEALLOCATE(self%masses)
    nullify(self%supercell)
    nullify(self%params)
    self%label="Destroyed lattice mover"
    ABI_DEALLOCATE(self%masses)
    ABI_DEALLOCATE(self%current_xcart)
    ABI_DEALLOCATE(self%current_vcart)
    ABI_DEALLOCATE(self%forces)
  end subroutine finalize

  subroutine set_params(self, params)
    ! set parameters from input file. (something else, like temperature for MvT calculation?)
    class(lattice_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(params)

  end subroutine set_params

  subroutine set_initial_state(self, mode)
    ! set initial positions, spin, etc
    class(lattice_mover_t), intent(inout) :: self
    integer, optional, intent(in) :: mode
    real(dp) :: xi(3, self%natom)
    integer :: i

    if(mode==1) then ! using a boltzmann distribution. 
       call self%rng%rand_normal_array(xi, 3*self%natom)
       do i=1, self%natom
          self%current_vcart(:,i) = xi(:, i) *sqrt(3.0*self%temperature/self%masses(i))
       end do
       self%current_xcart(:, :) = self%supercell%lattice%xcart(:,:)
    else 
       

    end if

  end subroutine set_initial_state

  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf)
    ! run one step. (For MC also?)
    class(lattice_mover_t), intent(inout) :: self
    ! array of effective potentials so that there can be multiple of them.
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    if(present(displacement) .or. present(strain)) then
       MSG_ERROR("displacement and strain should not be input for lattice mover")
    end if
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)

  end subroutine run_one_step

  subroutine reset(self)
    ! reset the state of mover (e.g. counter->0)
    ! so it can be reused.
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine reset

  subroutine calc_observables(self)
    ! call functions to calculate observables.
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine calc_observables

  subroutine write_hist(self)
    ! write to hist file
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine write_hist

  subroutine get_state(self, displacement, strain, spin, lwf, ihist)
    ! get the state of the ihist(th) step. ihist can be 0 (current), -1 (last), ... -maxhist..
    class(lattice_mover_t), intent(in):: self
    real(dp), optional, intent(inout) :: displacement, strain, spin, lwf
    integer, optional, intent(in):: ihist
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(ihist)
  end subroutine get_state

end module m_lattice_mover


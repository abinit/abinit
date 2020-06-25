!!****m* ABINIT/m_lattice_harmonic_potential
!! NAME
!! m_lattice_harmonic_potential
!!
!! FUNCTION
!! This module contains an harmonic lattice potential. 
!!
!! Datatypes:
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


module m_lattice_harmonic_potential
  use defs_basis
  use m_errors
  use m_abicore
  use m_abstract_potential, only: abstract_potential_t
  use m_spmat_coo, only: COO_mat_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_hashtable_strval, only: hash_table_t
  implicit none
!!***

  private

  type, public, extends(abstract_potential_t) :: lattice_harmonic_potential_t
     integer :: natom           ! number of atoms
     real(dp) :: ref_energy=0.0      ! reference energy
     !real(dp):: ref_cell(3,3)   ! reference structure cell parameters.Not needed.
     !real(dp), allocatable :: ref_xcart(:,:) ! reference xcart. not needed
     type(COO_mat_t) :: coeff  ! coefficient. A COO sparse matrix (3N*3N).
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_supercell
     procedure :: set_ref_energy
     procedure :: calculate
     procedure :: add_term
  end type lattice_harmonic_potential_t

contains


  !-------------------------------------------------------------------!
  ! initialize
  ! Input:
  !  natom: number of atoms
  !-------------------------------------------------------------------!
  subroutine initialize(self, natom)
    class(lattice_harmonic_potential_t), intent(inout) :: self
    integer, intent(in) :: natom
    self%has_displacement = .True.
    self%is_null = .False.
    self%label="Lattice_harmonic_potential"
    self%natom=natom
    call self%coeff%initialize(mshape= [3*self%natom, 3*self%natom])
    !ABI_MALLOC(self%ref_xcart, (3, self%natom))
  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(lattice_harmonic_potential_t), intent(inout) :: self
    self%has_displacement=.False.
    call self%coeff%finalize()
    call self%abstract_potential_t%finalize()
    self%is_null=.True.
    self%natom=0
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Add a term to the potential
  !-------------------------------------------------------------------!
  subroutine add_term(self, i,j, val)
    class(lattice_harmonic_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val
    call self%coeff%add_entry([i,j], val)
  end subroutine add_term


  !-------------------------------------------------------------------!
  ! Set the reference energy.
  !-------------------------------------------------------------------!
  subroutine set_ref_energy(self, ref_energy)
    class(lattice_harmonic_potential_t), intent(inout) :: self
    real(dp), intent(in) :: ref_energy
    self%ref_energy=ref_energy
  end subroutine set_ref_energy

  !-------------------------------------------------------------------!
  ! set_supercell
  !  link the supercell with potential.
  ! Inputs:
  !   supercell: mbsupercell_t
  !-------------------------------------------------------------------!
  subroutine set_supercell(self, supercell)
    class(lattice_harmonic_potential_t), intent(inout) :: self
    type(mbsupercell_t), target, intent(inout) :: supercell
    self%supercell => supercell
  end subroutine set_supercell


  !-------------------------------------------------------------------!
  ! calculate force and energy from harmonic potential
  ! F= - IFC .matmul. displacement
  ! E = 1/2 (-F) .dot. displacement = 1/2<disp|IFC|disp>
  ! Input:
  !   displacement: required.
  !-------------------------------------------------------------------!
  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, energy, energy_table)
    ! This function calculate the energy and its first derivative
    ! the inputs and outputs are optional so that each effpot can adapt to its
    ! own.
    ! In principle, the 1st derivatives are only calculated if asked to (present). However, they can be computed if it is simply convinient to do.
    class(lattice_harmonic_potential_t), intent(inout) :: self  ! the effpot may save the states.

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    type(hash_table_t),optional, intent(inout) :: energy_table
    real(dp) :: etmp
    real(dp) :: f(3,self%natom)
    ! if present in input
    ! calculate if required
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)
    ABI_UNUSED_A(lwf_force)

    etmp=0.0_dp

    call self%coeff%mv(displacement, f)
    if (present(force)) then
       force(:,:) = force(:,:) - f
    endif

    !energy =energy + 0.5_dp * sum(f*displacement)
    etmp=0.5_dp * sum(f*displacement)
    if (present(energy)) then
       energy=energy+etmp
    endif
    if(present(energy_table)) then
       call energy_table%put(self%label, etmp)
    endif
  end subroutine calculate

end module m_lattice_harmonic_potential

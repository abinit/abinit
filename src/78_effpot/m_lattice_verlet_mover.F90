!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_verlet_mover
!! NAME
!! m_lattice_verlet_mover
!!
!! FUNCTION
!! This module contains the verlet  (NVE) lattice mover.
!! 
!!
!!
!! Datatypes:
!!
!! * lattice_verlet_mover_t: defines the lattice movers
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

module m_lattice_verlet_mover
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

  type, public, extends(lattice_mover_t) :: lattice_verlet_mover_t
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_one_step
  end type lattice_verlet_mover_t

contains


  subroutine initialize(self,params, supercell, rng)
    class(lattice_verlet_mover_t), intent(inout) :: self
    type(multibinit_dtset_type), target, intent(in):: params
    type(mbsupercell_t), target, intent(in) :: supercell
    type(rng_t), target, intent(in) :: rng
    call self%lattice_mover_t%initialize(params, supercell, rng)
    self%label = "Velocity Verlet lattice mover"
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_verlet_mover_t), intent(inout) :: self
    call self%lattice_mover_t%finalize()
  end subroutine finalize


  !===================== run_one_step===============================!
  ! run one md step
  ! effpot: effective potential
  ! displacement: Should NOT be given, because it is stored in the mover already.
  ! strain: Should Not be given. Because 1) it is stored in the mover,
  !          and 2) this is a constant volume mover.
  ! spin: spin of atoms. Useful with spin-lattice coupling.
  ! lwf: lattice wannier function. Useful with lattice-lwf coupling (perhaps useless.)
  ! energy_table: energy table
  subroutine run_one_step(self, effpot,displacement, strain, spin, lwf, energy_table)
    class(lattice_verlet_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    integer :: i

    ! first half of velocity update. And full displacement update.
    ! v(t+1/2 dt) = v(t) + F/m * 1/2 dt
    ! x(t+dt) = x(t) + v(t+1/2 dt) * dt
    self%forces(:, :) =0.0
    self%energy = 0.0
    call effpot%calculate( displacement=self%displacement, strain=self%strain, &
         & spin=spin, lwf=lwf, force=self%forces, stress=self%stress,  energy=self%energy, energy_table=energy_table)
    do i=1, self%natom
       self%current_vcart(:,i) = self%current_vcart(:,i) + &
            & (0.5_dp * self%dt) * self%forces(:,i)/self%masses(i)
    end do
    !call self%force_stationary()

    self%displacement(:,:) = self%displacement(:,:)+self%current_vcart(:,:) * self%dt
    ! No need to update xcart.
    !self%current_xcart(:,:) = self%supercell%lattice%xcart(:,:) + self%displacement(:,:)


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
    !call self%force_stationary()
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(displacement)
  end subroutine run_one_step


end module m_lattice_verlet_mover


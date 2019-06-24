!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_langevin_mover
!! NAME
!! m_lattice_langevin_mover
!!
!! FUNCTION
!! This module contains the langevin  (NVT) lattice mover.
!! 
!!
!!
!! Datatypes:
!!
!! * lattice_langevin_mover_t: defines the lattice movers
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

module m_lattice_langevin_mover
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_mover, only: lattice_mover_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t

!!***

  implicit none

  private

  type, public, extends(lattice_mover_t) :: lattice_langevin_mover_t
     real(dp) :: c1, c2
     real(dp), allocatable :: c3(:), c4(:), c5(:)
     real(dp) :: fr   ! friction
     real(dp), allocatable :: xi(:), eta(:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_one_step
     procedure :: update_vars
  end type lattice_langevin_mover_t

contains


  subroutine initialize(self,params, supercell)
    class(lattice_langevin_mover_t), intent(inout) :: self
    type(multibinit_dtset_type), target, intent(in):: params
    type(mbsupercell_t), target, intent(in) :: supercell
    call self%lattice_mover_t%initialize(params, supercell)
    ABI_ALLOCATE(self%c3, (self%natom))
    ABI_ALLOCATE(self%c4, (self%natom))
    ABI_ALLOCATE(self%c5, (self%natom))
    ABI_ALLOCATE(self%xi, (3,self%natom))
    ABI_ALLOCATE(self%eta, (3,self%natom))

    call self%update_vars()
  end subroutine initialize


  subroutine finalize(self)
    class(lattice_langevin_mover_t), intent(inout) :: self
    ABI_DEALLOCATE(self%c3)
    ABI_DEALLOCATE(self%c4)
    ABI_DEALLOCATE(self%c5)
    ABI_DEALLOCATE(self%xi)
    ABI_DEALLOCATE(self%eta)
    call self%lattice_mover_t%finalize()
  end subroutine finalize
 

  subroutine update_vars(self):
    class(lattice_langevin_mover_t), intent(inout) :: self
    real(dp) :: dt, T, fr, sigma(self%natom)

    dt=self%dt
    T= self%temperature
    fr=self%fr
    sigma(:) = sqrt(2.0*T*fr/self%masses(:))
    self%c1 = dt / 2.0 - dt * dt * fr / 8.0
    self%c2 = dt * fr / 2.0 - dt * dt * fr * fr / 8.0
    self%c3 = sqrt(dt) * sigma / 2.0 - dt**1.5 * fr * sigma / 8.0
    self%c5 = dt**1.5 * sigma / (2.0 * sqrt(3.0))
    self%c4 = fr / 2. * self%c5
  end subroutine update_vars


  subroutine run_one_step(self, effpot,displacement, strain, spin, lwf )
    class(lattice_langevin_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    integer :: i

    call effpot%calculate( displacement=self%displacement, strain=self%strain, &
         & spin=spin, lwf=lwf, force=self%forces, stress=self%stress,  energy=self%energy)
    call self%rng%rand_normal_array(xi, 3*self%natom)
    call self%rng%rand_normal_array(eta, 3*self%natom)

    ! First half of velocity update
    do i =1, self%natom
       self%current_vcart(:, i) = self%current_vcart(:,i) + &
            & (self%c1 * self%forces(:,i) / self.masses(i) - &
            & self.c2 * self.current_vcart(:,i) + &
            & self.c3(i) * self.xi(:, i) - self.c4(i) * self.eta(:,i) )

       self%displacement(:, i) = self%displacement(:, i) &
            & + self%dt * self%current_vcart(:, i) * self%c5( i) *self%eta(:,i)
    end do

    ! second half, update the velocity but not the displacement.
    call effpot%calculate( displacement=self%displacement, strain=self%strain, &
         & spin=spin, lwf=lwf, force=self%forces, stress=self%stress,  energy=self%energy)
    do i =1, self%natom
       self%current_vcart(:, i) = self%current_vcart(:,i) + &
            & (self%c1 * self%forces(:,i) / self.masses(i) - &
            & self.c2 * self.current_vcart(:,i) + &
            & self.c3(i) * self.xi(:, i) - self.c4(i) * self.eta(:,i) )
    end do

  end subroutine run_one_step

end module m_lattice_langevin_mover


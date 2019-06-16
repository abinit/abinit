!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_Langevin_mover
!! NAME
!! m_lattice_Langevin_mover
!!
!! FUNCTION
!! This module contains the Langevin  (NVT) lattice mover.
!! 
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

module m_lattice_Langevin_mover
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_mover, only: lattice_movert_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t

!!***

  implicit none

  private

  type, public, extends(lattice_mover_t) :: lattice_Langevin_mover_t
     real(dp) :: fr, c1, c2, c3, c4, c5
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_one_step
  end type lattice_Langevin_mover_t

contains


  subroutine initialize(self)
    class(lattice_langevin_mover_t), intent(inout) :: self
    
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_langevin_mover_t), intent(inout) :: self
    call self%lattice_mover_t%finalize()
  end subroutine finalize

  subroutine update_vars(self)
    class(lattice_langevin_mover_t), intent(inout) :: self
    real(dp) :: dt, T, fr, sigma(self%natom)
    dt=self%dt
    T=self%temperature
    fr=self%fr
    sigma = sqrt(2.0*T*fr/self%masses)
    self%c1 = dt/2.0_dp - dt*dt * fr/ 8.0_dp
    self%c2 = dt * fr / 2 - dt * dt * fr * fr / 8._dp
    self%c3 = sqrt(dt) * sigma / 2.0_dp - dt**1.5 * fr * sigma / 8.
    self%c5 = dt**1.5 * sigma / (2.0_dp * sqrt(3.0_dp))
    self%c4 = fr / 2.0_dp * self.c5
  end subroutine update_vars
  
  subroutine run_one_step(self, effpot,displacement, strain, spin, lwf )
    class(lattice_langevin_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
  end subroutine run_one_step

end module m_lattice_Langevin_mover


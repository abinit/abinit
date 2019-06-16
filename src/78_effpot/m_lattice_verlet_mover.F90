!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_verlet_mover
!! NAME
!! m_lattice_Verlet_mover
!!
!! FUNCTION
!! This module contains the Verlet(NVE) lattice mover.
!! 
!!
!!
!! Datatypes:
!!
!! * lattice_verlet_mover_t: defines the lattice mover using verlet method
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

module m_lattice_Verlet_mover
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

  type, public, extends(lattice_mover_t) :: lattice_Verlet_mover_t
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: run_one_step
  end type lattice_Verlet_mover_t

contains


  subroutine initialize(self)
    class(lattice_Verlet_mover_t), intent(inout) :: self
    
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_Verlet_mover_t), intent(inout) :: self
    call self%lattice_mover_t%finalize()
  end subroutine finalize

  
  subroutine run_one_step(self, effpot,displacement, strain, spin, lwf )
    class(lattice_Verlet_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
  end subroutine run_one_step

end module m_lattice_verlet_mover


!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_effpot
!! NAME
!! m_lattice_effpot
!!
!! FUNCTION
!! place holder for lattice potential
!!
!!
!! Datatypes:
!!
!! * lattice_potential_t
!!
!! Subroutines:
!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (TO, hexu)
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

module m_lattice_effpot
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only : abstract_potential_t
  use m_hashtable_strval, only: hash_table_t
!!***

  implicit none
  private
  type ,public, extends(abstract_potential_t) :: lattice_effpot_t
   contains
     !procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     !procedure :: set_deformation
     !procedure :: get_force
     !procedure :: get_stress
     procedure :: calculate
  end type lattice_effpot_t

contains

  subroutine finalize(self)
    class(lattice_effpot_t), intent(inout) :: self
    ABI_UNUSED_A(self)
    MSG_ERROR("finalize for lattice_effpot not yet implemented")
  end subroutine finalize

  subroutine set_params(self, params)
    class(lattice_effpot_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(params)
    MSG_ERROR("set_params for lattice_effpot not yet implemented")
  end subroutine set_params


  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, &
          & energy, energy_table)
    class(lattice_effpot_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    type(hash_table_t), optional, intent(inout) :: energy_table
    if(present(force)) then
       force(:,:)=zero
    end if
    if (present(stress)) then
       stress(:,:)=zero
    end if
    if (present(energy)) then
       energy=zero
    end if
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(force)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)
    ABI_UNUSED_A(lwf_force)
    ABI_UNUSED_A(energy)
    ABI_UNUSED_A(energy_table)

    MSG_ERROR("calculate for lattice_effpot not yet implemented.")
  end subroutine calculate

end module m_lattice_effpot

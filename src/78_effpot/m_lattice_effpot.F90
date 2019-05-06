!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_potential
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

  implicit none
  private
  type ,public, extends(abstract_potential_t) :: lattice_effpot_t
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     !procedure :: set_deformation
     !procedure :: get_force
     !procedure :: get_stress
     procedure :: calculate
  end type lattice_effpot_t

contains

  subroutine initialize(self, params, fnames)
    class(lattice_effpot_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(inout) :: params
    character(*), intent(in) :: fnames(:)
    MSG_ERROR("initialize for lattice_effpot not yet implemented")
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_effpot_t), intent(inout) :: self
    MSG_ERROR("finalize for lattice_effpot not yet implemented")
  end subroutine finalize

  subroutine set_params(self, params)
    class(lattice_effpot_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    MSG_ERROR("set_params for lattice_effpot not yet implemented")
  end subroutine set_params

  ! subroutine set_deformation(self, displacements, strain )
  !   class(lattice_effpot_t), intent(inout) :: self
  !   real(dp), intent(in) :: displacements(:,:), strain(:,:)
  ! end subroutine set_deformation

  ! subroutine get_force(self, force)
  !   class(lattice_effpot_t), intent(inout) :: self
  !   real(dp), intent(out) :: force(:,:)
  ! end subroutine get_force

  ! subroutine get_stress(self, stress )
  !   class(lattice_effpot_t), intent(inout) :: self
  !   real(dp), intent(out) :: stress(:,:)
  ! end subroutine get_stress

  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, energy)
    class(lattice_effpot_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    if(present(force)) then
       force(:,:)=zero
    end if
    if (present(stress)) then
       stress(:,:)=zero
    end if
    if (present(energy)) then
       energy=zero
    end if
    MSG_ERROR("calculate for lattice_effpot not yet implemented.")
  end subroutine calculate

end module m_lattice_effpot

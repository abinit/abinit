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
  use m_effpot_api, only : effpot_t

  implicit none
  private
  type ,public, extends(effpot_t) :: lattice_effpot_t
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
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_effpot_t), intent(inout) :: self

  end subroutine finalize

  subroutine set_params(self, params)
    class(lattice_effpot_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
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

   subroutine calculate(self, displacement, strain, spin, force, stress, bfield, energy)
    class(lattice_effpot_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), optional, intent(in) :: displacement(:,:), strain(:,:), spin(:,:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
    if(present(force)) then
       force(:,:)=zero
    end if
    if (present(stress)) then
       stress(:,:)=zero
    end if
    if (present(energy)) then
       energy=zero
    end if
  end subroutine calculate
 
end module m_lattice_effpot

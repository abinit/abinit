#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lattice_api
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_effpot_api, only : effpot_t

  implicit none
  private
  type ,public, extends(effpot_t) :: lattice_api_t
   contains
     procedure :: set_params
     procedure :: get_force
     procedure :: get_stress
  end type lattice_api_t

contains

  subroutine set_params(self, params)
    class(lattice_api_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
  end subroutine set_params

  subroutine set_deformation(self, displacements, strain, )
    class(lattice_api_t), intent(inout) :: self
    real(dp), intent(in) :: displacements(:,:), strain(:,:)
  end subroutine set_deformation

  subroutine get_force(self, force)
    class(lattice_api_t), intent(inout) :: self
    real(dp), intent(out) :: force(:,:)
  end subroutine get_force

  subroutine get_stress(self, stress )
    class(lattice_api_t), intent(inout) :: self
    real(dp), intent(out) :: stress(:,:)
  end subroutine get_stress

end module m_lattice_api

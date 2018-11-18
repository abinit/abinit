#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_effpot_api
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_supercell, only: supercell_type

  implicit none
  private
  type ,public :: effpot_t
     ! This is the abstract class of effective potential.
     ! It do the following things:
     !  - read from file (which has the effective potential in primitive cell, corresponding to xml file)
     !  - build supercell.
     !  - calculate 0th, 1st,... derivative of energy related, e.g. Force for lattice model
     !type(prim_model_t) :: prim_model
     !type(sc_model_t) :: sc_model_t

     ! labels for variables
     logical :: has_displacement=.False.
     logical :: has_strain=.False.
     logical :: has_spin=.False.
     logical :: is_null=.False.   ! if is_null, this term does not exist.
     integer :: natom=0, nspin=0
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params      ! parameters from input file
     procedure :: read_potential  ! read effpot from file (primtive cell) (e.g. DDB, xml)
     procedure :: make_supercell  ! build supercell potential

     !procedure :: set_variables
     !procedure :: get_1st_deriv

     procedure :: set_distortion
     procedure :: set_spin
     procedure :: get_energy           ! energy
     procedure :: get_force            ! force
     procedure :: get_stress           ! stress
     procedure :: get_effective_Bfield ! effective Bfield
  end type effpot_t

contains

  subroutine initialize(self, params, fnames)
    class(effpot_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params  ! read from input file
    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).
    ! hexu comment  It's better to get rid of it and put everything into input file!?)

  end subroutine initialize

  subroutine read_potential(self, fnames)
    class(effpot_t), intent(inout) :: self
    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).
  end subroutine read_potential

  subroutine finalize(self)
    class(effpot_t), intent(inout) :: self
  end subroutine finalize

  subroutine set_params(self, params)
    class(effpot_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
  end subroutine set_params


  subroutine make_supercell(self, supercell)
    class(effpot_t), intent(inout) :: self
    type(supercell_type), intent(in) :: supercell
  end subroutine make_supercell

  ! hexu comment : which one is better, more general variables,
  !  or one function for each type of var?
  !subroutine set_variables(self, displacements, strain, spin)
  !  class(lattice_api_t), intent(inout) :: self
  !  real(dp), optional, intent(in) :: displacements(:,:), strain(:,:), spin(:,:)
  !end subroutine set_variables

  !subroutine get_1st_deriv(self, force, stress, bfield)
  !  class(lattice_api_t), intent(inout) :: self
  !  real(dp), optional, intent(out) :: force(:,:), stress(:,:), bfield(:,:)
  !end subroutine get_1st_deriv

  subroutine set_distortion(self, displacements, strain)
    class(effpot_t), intent(inout) :: self
    real(dp), optional, intent(in) :: displacements(:,:), strain(:,:)
  end subroutine set_distortion

  subroutine get_energy(self, energy)
    class(effpot_t), intent(inout) :: self
    real(dp) , intent(inout) :: energy
  end subroutine get_energy


  subroutine set_spin(self, spin)
    class(effpot_t), intent(inout) :: self
    real(dp), optional, intent(in) :: spin
  end subroutine set_spin

  subroutine get_force(self, force)
    class(effpot_t), intent(inout) :: self
    real(dp), intent(out) :: force(:,:)
  end subroutine get_force

  subroutine get_stress(self, stress)
    class(effpot_t), intent(inout) :: self
    real(dp), intent(out) :: stress(:,:)
  end subroutine get_stress

  subroutine get_effective_Bfield(self, bfield)
    class(effpot_t), intent(in) :: self
    real(dp) :: bfield(:,:)
  end subroutine get_effective_Bfield

end module m_effpot_api

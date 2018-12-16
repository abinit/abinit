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
     integer :: natoms=0, nspins=0
     real(dp), allocatable :: ms(:)
   contains
     procedure :: set_params      ! parameters from input file
     procedure :: read_potential  ! read effpot from file (primtive cell) (e.g. DDB, xml)
     procedure :: make_supercell  ! build supercell potential
     procedure :: calculate       ! get energy and 1st derivative from input state
     procedure :: get_delta_E=> effpot_t_get_delta_E ! calculate energy diffence if one component is changed for Monte carlo algorithm

     !procedure :: set_variables
     !procedure :: get_1st_deriv

     procedure :: set_distortion
     procedure :: set_spin
     !procedure :: get_energy           ! energy
     !procedure :: get_force            ! force
     !procedure :: get_stress           ! stress
     !procedure :: get_effective_Bfield ! effective Bfield
  end type effpot_t

  type, public:: effpot_pointer_t ! pointer to effpot
     class(effpot_t) , pointer :: obj
  end type effpot_pointer_t

  type, public, extends(effpot_t) :: effpot_list_t
     integer :: size=0
     integer :: MAXSIZE=8
     type(effpot_pointer_t) :: effpots(8)  ! Maximum number of 8, should we make it more dynamic?
     real(dp) :: stress_tmp(3,3)
     real(dp), allocatable :: force_tmp(:,:), bfield_tmp(:,:)
   contains
     procedure :: initialize => effpot_list_t_initialize
     procedure :: finalize => effpot_list_t_finalize
     procedure :: append => effpot_list_t_append
     procedure :: calculate => effpot_list_t_calculate
  end type effpot_list_t

contains

  subroutine read_potential(self, fnames)
    class(effpot_t), intent(inout) :: self
    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).
    ! TODO
  end subroutine read_potential

  subroutine set_params(self, params)
    class(effpot_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    ! TODO
  end subroutine set_params


  subroutine make_supercell(self, supercell)
    class(effpot_t), intent(inout) :: self
    type(supercell_type), intent(in) :: supercell
    ! TODO 
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

  subroutine set_distortion(self, displacement, strain)
    class(effpot_t), intent(inout) :: self
    real(dp), optional, intent(in) :: displacement(:,:), strain(:,:)
  end subroutine set_distortion

  subroutine set_spin(self, spin)
    class(effpot_t), intent(inout) :: self
    real(dp), optional, intent(in) :: spin
  end subroutine set_spin

  subroutine calculate(self, displacement, strain, spin, force, stress, bfield, energy)
    class(effpot_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), optional, intent(in) :: displacement(:,:), strain(:,:), spin(:,:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
    ! if present in input
    ! calculate if required
  end subroutine calculate

  subroutine effpot_t_get_delta_E(self, S, ispin, Snew, deltaE)
    ! for spin monte carlo
    ! calculate energy difference if one spin is moved.
    class(effpot_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), intent(in) :: S(:,:),  Snew(:)
    integer, intent(in) :: ispin
    real(dp), intent(out) :: deltaE
  end subroutine effpot_t_get_delta_E

!   subroutine get_energy(self, energy)
!     class(effpot_t), intent(inout) :: self
!     real(dp) , intent(inout) :: energy
!   end subroutine get_energy


!   subroutine get_force(self, force)
!     class(effpot_t), intent(inout) :: self
!     real(dp), intent(out) :: force(:,:)
!   end subroutine get_force

!   subroutine get_stress(self, stress)
!     class(effpot_t), intent(inout) :: self
!     real(dp), intent(out) :: stress(:,:)
!   end subroutine get_stress

!   subroutine get_effective_Bfield(self, spin,bfield)
!     class(effpot_t), intent(in) :: self
!     real(dp), intent(in) :: spin(:,:)
!     real(dp), intent(inout) :: bfield(:,:)
!   end subroutine get_effective_Bfield


  !======================== Effpot_list_t======================
  ! TODO : is this an overkill?

  subroutine effpot_list_t_initialize(self, natoms, nspins)
    class (effpot_list_t), intent(inout) :: self
    integer, intent(in) :: natoms, nspins
    self%natoms=natoms
    self%nspins=nspins
    ABI_ALLOCATE(self%force_tmp, (3, natoms))
    ABI_ALLOCATE(self%bfield_tmp, (3, nspins))
  end subroutine effpot_list_t_initialize

  subroutine effpot_list_t_finalize(self)
    class (effpot_list_t) ,intent(inout) :: self
    integer :: i
    do i=1, self%size
       nullify(self%effpots(i)%obj)
    end do
    self%size=0
    self%is_null=.True.
    self%has_displacement=.False.
    self%has_strain=.False.
    self%has_spin=.False.

    ABI_DEALLOCATE(self%force_tmp)
    ABI_DEALLOCATE(self%bfield_tmp)
  end subroutine effpot_list_t_finalize


  subroutine effpot_list_t_append(self, effpot)
    class (effpot_list_t) :: self
    class (effpot_t), target :: effpot
    self%size=self%size + 1
    if(self%size>self%MAXSIZE) then
       write(std_err, *) "Number of effpot larger than the maximum 8"
    end if
    self%effpots(self%size)%obj => effpot
    self%is_null= (self%is_null .or. effpot%is_null)
    self%has_spin= (self%has_spin .or. effpot%has_spin)
    self%has_displacement= (self%has_displacement .or. effpot%has_displacement)
    self%has_strain= (self%has_strain.or. effpot%has_strain)
  end subroutine effpot_list_t_append

  subroutine effpot_list_t_calculate(self, displacement, strain, spin, force, stress, bfield, energy)
    class(effpot_list_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), optional, intent(in) :: displacement(:,:), strain(:,:), spin(:,:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
    integer :: i
    real(dp) :: e
    ! if present in input
    ! calculate if required
    if( present(force) .and. present(stress)) then
       force(:,:)=0.0
       stress(:,:)=0.0
    endif

    if(present(bfield)) then
       bfield(:,:)=0.0
    end if

    energy=0.0

    ! calculate force and strain if asked to
    if (present(force) .and. present(strain)) then
       do i=1, self%size
          ! only for these effpot has force term
          if(self%effpots(i)%obj%has_displacement .and. self%effpots(i)%obj%has_strain) then
             call self%effpots(i)%obj%calculate(displacement=displacement, strain=strain, &
                     & spin=spin, force=self%force_tmp, stress=self%stress_tmp, energy=e)
          endif
          force(:,:) = force(:,:)+self%force_tmp(:,:)
          stress(:,:) = stress(:,:) + self%stress_tmp(:,:)
          energy=energy+e
       end do
    else if(present(bfield)) then
       do i=1, self%size
          ! only for these effpot has force term
          if(self%effpots(i)%obj%has_spin) then
             ! uncomment if F2008 style passing optional is not allowed.
             !if (present(displacement) .and. present(strain)) then
             call self%effpots(i)%obj%calculate(displacement=displacement, strain=strain, &
                     & spin=spin, bfield=self%bfield_tmp, energy=e)
             !else
             !   call self%effpots(i)%obj%calculate( spin=spin,  &
             !        & bfield=bfield, energy=e)
             !endif
          endif
          bfield(:,:)=bfield(:,:)+self%bfield_tmp
       end do
    endif

  end subroutine effpot_list_t_calculate

end module m_effpot_api

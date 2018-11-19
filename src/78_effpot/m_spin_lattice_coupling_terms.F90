#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_spin_lattice_coupling_terms
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_effpot_api, only: effpot_t

  implicit none
  private

  type, public, extends(effpot_t) :: abstract_spin_lattice_coupling_term_t
   contains
     ! things all spin_lattice_couplgin terms share
  end type abstract_spin_lattice_coupling_term_t

  type, public, extends(abstract_spin_lattice_coupling_term_t) :: slc_Oiju_t
   contains
     procedure, public :: initialize => Oiju_initialize
     procedure, public :: finalize => Oiju_finalize
     !procedure, public :: get_force => Oiju_get_force
     !procedure, public :: get_bfield => Oiju_get_bfield
     procedure, public :: calculate => Oiju_calculate
  end type slc_Oiju_t

  type, public, extends(abstract_spin_lattice_coupling_term_t) :: slc_Tijuv_t
   contains
     procedure, public :: initialize => Tijuv_initialize
     procedure, public :: finalize => Tijuv_finalize
     !procedure, public :: get_force => Tijuv_get_force
     !procedure, public :: get_bfield=> Tijuv_get_bfield
     procedure, public :: calculate => Tijuv_calculate
  end type slc_Tijuv_t


  type, public, extends(abstract_spin_lattice_coupling_term_t) :: spin_lattice_coupling_term_t
     type(slc_Oiju_t):: Oiju
     type(slc_Tijuv_t) :: Tijuv
   contains
     procedure, public :: initialize => slc_initialize
     procedure, public :: finalize => slc_finalize
     !procedure, public :: get_force => slc_get_force
     !procedure, public :: get_bfield => slc_get_bfield
     procedure, public :: calculate=> slc_calculate
  end type spin_lattice_coupling_term_t


contains
!!! --------------- Oiju -----------------

 subroutine Oiju_initialize(self, params, fnames)
   class(slc_Oiju_t), intent(inout) :: self
   type(multibinit_dtset_type) :: params  ! read from input file
   character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).

 end subroutine Oiju_initialize

 subroutine Oiju_finalize(self)
   class(slc_Oiju_t) , intent(inout):: self
 end subroutine Oiju_finalize

 subroutine Oiju_get_force(self, force)
   class(slc_Oiju_t) , intent(inout):: self
   real(dp), intent(out) :: force(:,:)
 end subroutine Oiju_get_force

 subroutine Oiju_get_force(self, force)
   class(slc_Oiju_t) , intent(inout):: self
   real(dp), intent(out) :: force(:,:)
 end subroutine Oiju_get_force



!!! -------------- Tijuv ---------------
 subroutine Tijuv_initialize(self, params, fnames)
   class(slc_Tijuv_t), intent(inout) :: self
   type(multibinit_dtset_type) :: params  ! read from input file
   character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).

 end subroutine Tijuv_initialize

 subroutine Tijuv_finalize(self)
   class(slc_Tijuv_t), intent(inout) :: self
 end subroutine Tijuv_finalize

 subroutine Tijuv_get_force(self, force)
   class(slc_Tijuv_t) , intent(inout):: self
   real(dp), intent(out) :: force(:,:)
 end subroutine Tijuv_get_force

 subroutine slc_calculate(self, displacement, strain, spin, force, stress, bfield, energy)
   class(spin_lattice_coupling_term_t), intent(inout) :: self  ! the effpot may save the states.
   real(dp), optional, intent(in) :: displacement(:,:), strain(:,:), spin(:,:)
   real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
 end subroutine slc_calculate



!!! ------------- spin_lattice_coupling_term--------
 subroutine slc_initialize(self, params, fnames)
   class(spin_lattice_coupling_term_t), intent(inout) :: self
   type(multibinit_dtset_type) :: params  ! read from input file
   character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).

   call self%Oiju%initialize(params, fnames)
   call self%Tijuv%initialize(params, fnames)
 end subroutine slc_initialize

 subroutine slc_make_supercell(self, supercell)
   class(spin_lattice_coupling_term_t), intent(inout) :: self
   type(supercell_type), intent(in) :: supercell
 end subroutine slc_make_supercell


 subroutine slc_finalize(self)
   class(spin_lattice_coupling_term_t), intent(inout) :: self
   if(.not. self%Oiju%is_null) then
      call self%Oiju%finalize()
   end if

   if(.not. self%Tijuv%is_null) then
      call self%Tijuv%finalize()
   end if

 end subroutine slc_finalize

 subroutine slc_calculate(self, displacement, strain, spin, force, stress, bfield, energy)
   class(spin_lattice_coupling_term_t), intent(inout) :: self  ! the effpot may save the states.
   real(dp), optional, intent(in) :: displacement(:,:), strain(:,:), spin(:,:)
   real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
 end subroutine slc_calculate


end module m_spin_lattice_coupling_terms

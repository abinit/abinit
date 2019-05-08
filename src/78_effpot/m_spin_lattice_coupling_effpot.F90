!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_lattice_coupling_effpot
!! NAME
!! m_spin_lattice_coupling_effpot
!!
!! FUNCTION
!! This module define the spin-lattice coupling terms.
!! The terms are derived from the abstract_potential_t type.
!! It can calculate the energy, first derivative to spin or lattice.
!!
!! Datatypes:
!!  base_spin_lattice_coupling_effpot_t
!!  slc_Oiju_t: Oiju (one lattice, two spin)
!!  slc_Tijuv_t: Tijuv term (two lattice, two spin)
!!
!! Subroutines:
!! TODO: add this when F2003 documentation system is determined
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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

module m_spin_lattice_coupling_effpot
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_polynomial_potential, only: polynomial_potential_t
  use m_supercell, only : supercell_type

  implicit none
!!***
  private

  type, public, extends(polynomial_potential_t) :: slc_Oiju_t
   !contains
     !procedure, public :: initialize => Oiju_initialize
     !procedure, public :: finalize => Oiju_finalize
     !procedure, public :: calculate => Oiju_calculate
  end type slc_Oiju_t

  type, public, extends(polynomial_potential_t) :: slc_Tijuv_t
   !contains
     !procedure, public :: initialize => Tijuv_initialize
     !procedure, public :: finalize => Tijuv_finalize
     !procedure, public :: calculate => Tijuv_calculate
  end type slc_Tijuv_t


  type, public, extends(abstract_potential_t) :: spin_lattice_coupling_effpot_t
     type(slc_Oiju_t):: Oiju
     type(slc_Tijuv_t) :: Tijuv
   !contains
     !procedure, public :: initialize => slc_initialize
     !procedure, public :: finalize => slc_finalize
     !procedure, public :: get_force => slc_get_force
     !procedure, public :: get_bfield => slc_get_bfield
     !procedure, public :: calculate=> slc_calculate
  end type spin_lattice_coupling_effpot_t


contains
!!! --------------- Oiju -----------------

!  subroutine Oiju_initialize(self, params, fnames)
!    class(slc_Oiju_t), intent(inout) :: self
!    type(multibinit_dtset_type) :: params  ! read from input file
!    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).

!  end subroutine Oiju_initialize

!  subroutine Oiju_finalize(self)
!    class(slc_Oiju_t) , intent(inout):: self
!  end subroutine Oiju_finalize

!  subroutine Oiju_get_force(self, force)
!    class(slc_Oiju_t) , intent(inout):: self
!    real(dp), intent(out) :: force(:,:)
!  end subroutine Oiju_get_force

!  subroutine Oiju_get_stress(self, force)
!    class(slc_Oiju_t) , intent(inout):: self
!    real(dp), intent(out) :: force(:,:)
!  end subroutine Oiju_get_stress

!  subroutine Oiju_calculate(self, displacement, strain, spin, lwf, force, stress, bfield,lwf_force, energy)
!      class(slc_Oiju_t), intent(inout) :: self  ! the effpot may save the states.
!      real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
!      real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
!      ! if present in input
!      ! calculate if required
!    end subroutine Oiju_calculate




! !!! -------------- Tijuv ---------------
!  subroutine Tijuv_initialize(self, params, fnames)
!    class(slc_Tijuv_t), intent(inout) :: self
!    type(multibinit_dtset_type) :: params  ! read from input file
!    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).

!  end subroutine Tijuv_initialize

!  subroutine Tijuv_finalize(self)
!    class(slc_Tijuv_t), intent(inout) :: self
!  end subroutine Tijuv_finalize

!  subroutine Tijuv_get_force(self, force)
!    class(slc_Tijuv_t) , intent(inout):: self
!    real(dp), intent(out) :: force(:,:)
!  end subroutine Tijuv_get_force

!  subroutine Tijuv_calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, energy)
!    class(slc_Tijuv_t), intent(inout) :: self  ! the effpot may save the states.
!    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
!    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
!  end subroutine Tijuv_calculate



! !!! ------------- spin_lattice_coupling_effpot--------
!  subroutine slc_initialize(self, params, fnames)
!    class(spin_lattice_coupling_effpot_t), intent(inout) :: self
!    type(multibinit_dtset_type) :: params  ! read from input file
!    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).

!    call self%Oiju%initialize(params, fnames)
!    call self%Tijuv%initialize(params, fnames)
!  end subroutine slc_initialize

!  subroutine slc_make_supercell(self, supercell)
!    class(spin_lattice_coupling_effpot_t), intent(inout) :: self
!    type(supercell_type), intent(in) :: supercell
!  end subroutine slc_make_supercell


!  subroutine slc_finalize(self)
!    class(spin_lattice_coupling_effpot_t), intent(inout) :: self
!    if(.not. self%Oiju%is_null) then
!       call self%Oiju%finalize()
!    end if

!    if(.not. self%Tijuv%is_null) then
!       call self%Tijuv%finalize()
!    end if

!  end subroutine slc_finalize

!  subroutine slc_calculate(self, displacement, strain, spin, lwf, force, stress, bfield,lwf_force, energy)
!    class(spin_lattice_coupling_effpot_t), intent(inout) :: self  ! the effpot may save the states.
!    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
!    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
!  end subroutine slc_calculate


end module m_spin_lattice_coupling_effpot

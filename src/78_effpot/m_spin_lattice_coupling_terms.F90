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
     integer :: nspin, natom
   contains
     ! things all spin_lattice_couplgin terms share
  end type abstract_spin_lattice_coupling_term_t

  type, public, extends(abstract_spin_lattice_coupling_term_t) :: slc_Oiju_t
   contains
     procedure, public :: initialize => Oiju_initialize
     procedure, public :: finalize => Oiju_finalize
     procedure, public :: get_force => Oiju_get_force
  end type slc_Oiju_t

  type, public, extends(abstract_spin_lattice_coupling_term_t) :: slc_Tijuv_t
   contains
     procedure, public :: initialize => Tijuv_initialize
     procedure, public :: finalize => Tijuv_finalize
     procedure, public :: get_force => Tijuv_get_force
  end type slc_Tijuv_t


  type, public, extends(abstract_spin_lattice_coupling_term_t) :: spin_lattice_coupling_term_t
     type(slc_Oiju_t):: Oiju
     type(slc_Tijuv_t) :: Tijuv
   contains
     procedure, public :: initialize => slc_initialize
     procedure, public :: finalize => slc_finalize
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



!!! ------------- spin_lattice_coupling_term--------
 subroutine slc_initialize(self, params, fnames)
   class(spin_lattice_coupling_term_t), intent(inout) :: self
   type(multibinit_dtset_type) :: params  ! read from input file
   character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).


   call self%Oiju%initialize(params, fnames)
   call self%Tijuv%initialize(params, fnames)
 end subroutine slc_initialize

 subroutine slc_finalize(self)
   class(spin_lattice_coupling_term_t), intent(inout) :: self
   if(.not. self%Oiju%is_null) then
      call self%Oiju%finalize()
   end if

   if(.not. self%Tijuv%is_null) then
      call self%Tijuv%finalize()
   end if

 end subroutine slc_finalize


end module m_spin_lattice_coupling_terms

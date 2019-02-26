!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_primitive_potential_list
!! NAME
!! m_primitive_potential_list
!!
!! FUNCTION
!! This module define the primitive potential list type, which is a list of primitive potentials
!! 
!! Datatypes:
!!  primitive_potential_list_t
!!
!! Subroutines:
!! 
!!  * fill_supercell: use translation symmetry to fill the supercell.
!!  * load_from_file: load potential from file.
!!  * save_to_file: save to file.
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

module m_primitive_potential_list
  use defs_basis
  use m_abicore
  use m_errors
  use m_multibinit_dataset , only: multibinit_dtset_type
  use m_unitcell, only: unitcell_t
  use m_supercell_maker, only: supercell_maker_t
  use m_abstract_potential, only: abstract_potential_t, effpot_list_t
  use m_primitive_potential, only: primitive_potential_t
  implicit none


  type, public:: primitive_potential_pointer_t ! pointer to effpot
     class(primitive_potential_t), pointer :: obj
  end type primitive_potential_pointer_t


  type, public, extends(primitive_potential_t):: primitive_potential_list_t
     type(primitive_potential_pointer_t), allocatable :: pots(:)
     type(unitcell_t), pointer :: primcell
     integer :: size=0
     integer :: capacity=0
   contains
     procedure :: initialize => initialize
     procedure :: finalize => finalize
     procedure :: append => append
     procedure :: load_from_files
     procedure :: fill_supercell
     procedure :: fill_supercell_list
  end type primitive_potential_list_t
contains

subroutine fill_supercell(self, sc_maker, sc_pot)
    class(primitive_potential_list_t), intent(inout) :: self
    class(supercell_maker_t), intent(in) :: sc_maker
    class(abstract_potential_t), pointer, intent(inout) :: sc_pot
    ! Note that sc_pot is a pointer

    ! use a pointer to the specific potential which will be filled
    ! e.g. type(spin_potential_t), pointer :: tmp
    type(abstract_potential_t), pointer :: tmp
    allocate(abstract_potential_t:: tmp)
    ! call tmp%initialize(....)
    ! set tmp
    sc_pot=>tmp
    nullify(tmp)
  end subroutine fill_supercell

  subroutine load_from_files(self, params,  fnames)
    class(primitive_potential_list_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
  end subroutine load_from_files


  subroutine save_to_file(self, fname)
    class(primitive_potential_list_t), intent(inout) :: self
    character(*), intent(in) :: fname
  end subroutine save_to_file


  subroutine initialize(self)
    class(primitive_potential_list_t), intent(inout):: self
    self%size=0
  end subroutine initialize

  subroutine finalize(self)
    class(primitive_potential_list_t), intent(inout):: self
    self%size=0
  end subroutine finalize

  subroutine append(self, pot)
    class(primitive_potential_list_t), intent(inout):: self
    class(primitive_potential_t), intent(inout) :: pot
  end subroutine append

  subroutine fill_supercell_ptr(self, sc_maker, sc_pot)
    class(primitive_potential_list_t), intent(inout) :: self
    class(supercell_maker_t), intent(in) :: sc_maker
    class(abstract_potential_t), pointer, intent(inout) :: sc_pot
    ! Note that sc_pot is a pointer
    ! use a pointer to the specific potential which will be filled
    type(effpot_list_t), pointer :: tmp
    allocate(tmp)
    call self%fill_supercell_list( sc_maker, tmp)
    ! call tmp%initialize(....)
    ! set tmp
    sc_pot=>tmp
    nullify(tmp)
  end subroutine fill_supercell_ptr

  subroutine fill_supercell_list(self, sc_maker, sc_pots)
    class(primitive_potential_list_t), intent(inout) :: self
    class(supercell_maker_t), intent(in) :: sc_maker
    class(effpot_list_t), intent(inout) :: sc_pots
    ! Note that sc_pot is a pointer
    ! use a pointer to the specific potential which will be filled
    class(abstract_potential_t), pointer :: tmp
    integer i
    do i =1, self%size
      call self%pots(i)%obj%fill_supercell(sc_maker, tmp)
    end do
    call sc_pots%append(tmp)
    nullify(tmp)
  end subroutine fill_supercell_list

end module m_primitive_potential_list

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_primitive_potential
!! NAME
!! m_primitive_potential
!!
!! FUNCTION
!! This module define the primitive potential type, which is a direct map to the xml file (or other format).
!! 
!! Datatypes:
!!  primitive_potential_t
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

module m_primitive_potential
  use defs_basis
  use m_abicore
  use m_errors
  use m_supercell
  use m_un
  use m_supercell_maker, only: supercell_maker_t
  use m_abstract_potential, only: abstract_potential_t
  implicit none

  type ,public :: primitive_potential_t
     ! This is the abstract class of primitive potential
     ! It do the following things:
     type(unitcell_t)
   contains
     !procedure:: initialize       ! perhaps each effpot type should have own 
     !procedure :: finalize
     procedure :: fill_supercell
     procedure :: load_from_file
  end type primitive_potential_t


  type, public:: primitive_potential_pointer_t ! pointer to effpot
     class(primitive_potential_t), pointer :: obj
  end type primitive_potential_pointer_t


  type, public, extends(primitive_potential_t):: primtive_potential_list_t
     type(primtive_potential_pointer_t), allocatable :: pots(:)
     type(primitive_cell_t), pointer :: primcell
     integer :: size
   contains
     procedure :: initialize => list_initialize
     procedure :: finalize => list_finalize
     procedure :: append => list_append
     procedure :: fill_supercell => list_fill_supercell
  end type primtive_potential_list_t
contains

function fill_supercell(self, sc_maker) result(sc_pot)
    class(primitive_potential_t), intent(inout) :: self
    class(supercell_maker_t), intent(in) :: sc_maker
    class(abstract_potential_t), pointer :: sc_pot
    ! Note that sc_pot is a pointer

    ! use a pointer to the specific potential which will be filled
    type(abstract_potential_t), pointer :: tmp
    allocate(abstract_potential_t:: tmp)

    ! call tmp%initialize(....)
    ! set tmp
    sc_pot=>tmp
    nullify(tmp)
  end subroutine fill_supercell

  subroutine load_from_file(self, fname)
    class(primitive_potential_t), intent(inout) :: self
    character(*), intent(in) :: fname
  end subroutine load_from_file

  
  subroutine save_to_file(self, fname)
    class(primitive_potential_t), intent(inout) :: self
    character(*), intent(in) :: fname
  end subroutine save_to_file


  subroutine list_initialize(self)
    class(primitive_potential_list_t), intent(inout):: self
    self%size=0
  end subroutine list_initialize

end module m_primitive_potential

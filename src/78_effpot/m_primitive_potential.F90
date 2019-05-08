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
  use m_multibinit_dataset , only: multibinit_dtset_type
  use m_unitcell, only: unitcell_t
  use m_supercell_maker, only: supercell_maker_t
  use m_abstract_potential, only: abstract_potential_t
  use m_potential_list, only: potential_list_t
  use m_mpi_scheduler, only: init_mpi_info
  implicit none
  private
!!***


  type ,public :: primitive_potential_t
     ! This is the abstract class of primitive potential
     ! It do the following things:
     type(unitcell_t), pointer :: primcell=>null()
     character(len=200) :: label="Abstract primitive potential"
     logical::has_displacement=.False.
     logical::has_strain=.False.
     logical::has_spin=.False.
     logical::has_lwf=.False.
   contains
     !procedure:: initialize       ! perhaps each effpot type should have own 
     procedure :: finalize
     procedure :: fill_supercell
     procedure :: load_from_files
  end type primitive_potential_t

contains

subroutine finalize(self)
  class(primitive_potential_t), intent(inout) :: self
  nullify(self%primcell)
  self%label="Destroyed primitive potential"
  self%has_displacement=.False.
  self%has_strain=.False.
  self%has_spin=.False.
  self%has_lwf=.False.
end subroutine finalize

subroutine fill_supercell(self, scmaker, scpot)
    class(primitive_potential_t), intent(inout) :: self
    type(supercell_maker_t), intent(inout) :: scmaker
    class(abstract_potential_t), pointer, intent(inout) :: scpot
    ! Note that sc_pot is a pointer

    ! use a pointer to the specific potential which will be filled
    ! e.g. type(spin_potential_t), pointer :: tmp
    ! call tmp%initialize(....)
    ! set tmp
  end subroutine fill_supercell

  subroutine load_from_files(self, params,  fnames)
    class(primitive_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    write(*, *) "load_from_files for abstract primitive potential is not implemented."
  end subroutine load_from_files


  subroutine save_to_file(self, fname)
    class(primitive_potential_t), intent(inout) :: self
    character(*), intent(in) :: fname
  end subroutine save_to_file

end module m_primitive_potential

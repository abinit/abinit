!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_potential_list
!! NAME
!! m_potential_list
!!
!! FUNCTION
!! This module contains the potential_list potential. It is made of list of pointers to potentials,
!! and has all the functionality of a potential
!! It is to represent the equation H=\sum_i H_i
!! with dH/dx = \sum_i dH_i/dx
!!
!! Datatypes:
!!
!! * potential_list_t: list of abstract_potential_t, which is essentially a list of pointer to abstract_potential_t
!!    itself is also a effpot type, and its energy, 1st derivative to energy are the sum of all items.
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
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

module m_potential_list
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_mpi_scheduler, only: init_mpi_info
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t


  implicit none
  !!***

  !!****t* defs_abitypes/effpot_pointer_t
  !! NAME
  !! effpot_pointer_t
  !!
  !! FUNCTION
  !! class pointer to abstract potentials. Only here because
  !! Fortran does not allow array of pointers! So we must put
  !! pointer to an type and make an array of this type.
  !! It should only be used in potential_list_t.
  !!
  !! SOURCE
  type, private:: effpot_pointer_t ! pointer to effpot
     class(abstract_potential_t) , pointer :: ptr=>null()
  end type effpot_pointer_t
  !!***

  !!****t* defs_abitypes/potential_list_t
  !! NAME
  !! potential_list_t
  !!
  !! FUNCTION
  !! A potential which is made of a list of potentials.
  !! H=\sum_i H_i
  !!
  !! SOURCE
  type, public, extends(abstract_potential_t) :: potential_list_t
     integer :: size=0  ! size is the USED number of potentials in the list
     integer :: capacity=0  ! capacity is the number of places allocated for potentials
     type(effpot_pointer_t), allocatable :: list(:) ! A list of pointer type to potentials
   contains
     procedure :: initialize ! make an empty list
     procedure :: ipot       ! return the i'th potential
     procedure :: set_supercell ! pointer to supercell and set_supercell for all pot in list
     procedure :: finalize
     procedure :: append  ! add a potential to the list
     procedure :: calculate ! each potential in list do calculate and then sum.
     procedure :: get_delta_E ! currently only used in spin potential, 
                              ! to calculate energy difference when one spin changes.
  end type potential_list_t
  !!***

contains

  !----------------------------------------------------------------------
  !> @brief initialize 
  !>
  !----------------------------------------------------------------------
  subroutine initialize(self)
    class (potential_list_t), intent(inout) :: self
    self%size=0
    self%capacity=0
    self%label="ListPotential"
  end subroutine initialize

  !----------------------------------------------------------------------
  !> @brief    return the i'th potential in the list
  !>
  !> @param[in]  i: the index of the spin potentails
  ! Unfortunately, fortran does not allow things  call potlist%ipot(i)%method...
  !----------------------------------------------------------------------
  function ipot(self, i)
    class (potential_list_t), target, intent(in) :: self
    integer, intent(in) :: i
    class(abstract_potential_t), pointer :: ipot
    ipot=>self%list(i)%ptr
  end function ipot

  !----------------------------------------------------------------------
  !> @brief  set_supercell for every member of the potential list
  !> @param[in]  supercell
  !----------------------------------------------------------------------
  subroutine set_supercell(self, supercell)
    class (potential_list_t), intent(inout) :: self
    type (mbsupercell_t), target, intent(inout) :: supercell
    integer :: i
    self%supercell => supercell
    do i=1, self%size
       call self%list(i)%ptr%set_supercell(supercell)
    end do
  end subroutine set_supercell

  !----------------------------------------------------------------------
  !> @brief finalize
  !----------------------------------------------------------------------
  subroutine finalize(self)
    class (potential_list_t) ,intent(inout) :: self
    integer :: i
    do i=1, self%size
       call self%list(i)%ptr%finalize()
       ! Intel compiler complains
       if(associated(self%list(i)%ptr)) then
          ABI_DATATYPE_DEALLOCATE_SCALAR(self%list(i)%ptr)
       endif 
       
       nullify(self%list(i)%ptr)
    end do
    if (allocated(self%list)) then
       ABI_DEALLOCATE(self%list)
    end if
    self%size=0
    self%capacity=0
    self%is_null=.True.
    self%has_displacement=.False.
    self%has_strain=.False.
    self%has_spin=.False.
  end subroutine finalize


  !----------------------------------------------------------------------
  !> @brief append a potential to the list
  !> It also update the has_* according to the added effpot.
  !> @param[in]  effpot: the effective potential to be added.
  !----------------------------------------------------------------------
  subroutine append(self, effpot)
    ! Add a pointer to an effpot term to list.
    class (potential_list_t) :: self
    class (abstract_potential_t), target :: effpot
    type(effpot_pointer_t), allocatable :: temp(:)
    self%size=self%size + 1
    if(self%size==1) then
       self%capacity=8
       ABI_MALLOC(self%list, (self%capacity))
    else if ( self%size>self%capacity ) then
       ! fancy increasing equation to allocate new array.
       self%capacity = self%size + self%size / 4 + 8
       ABI_MALLOC(temp, (self%capacity))
       temp(:self%size-1) = self%list(:)
       ABI_MOVE_ALLOC(temp, self%list) !temp gets deallocated
    end if
    self%list(self%size)%ptr=>effpot
    self%is_null= (self%is_null .and. effpot%is_null)
    self%has_spin= (self%has_spin .or. effpot%has_spin)
    self%has_displacement= (self%has_displacement .or. effpot%has_displacement)
    self%has_strain= (self%has_strain.or. effpot%has_strain)
    self%has_lwf =(self%has_lwf.or. effpot%has_lwf)
  end subroutine append

  !----------------------------------------------------------------------
  !> @brief calculate energy and 1st derivatives
  !>  The list will add the sum of all the results from its components.
  !>  If one optional variable is passed to a subroutine, its "presence" will be kept.
  !> @param[in] (optional) displacement, strain, spin, lwf
  !> @param[out](optional) force, stress, bfield, lwf_force, energy
  !----------------------------------------------------------------------
  subroutine calculate(self, displacement, strain, spin, lwf,  force, stress, bfield, lwf_force, energy)
    ! calculate energy and its first derivatives.
    class(potential_list_t), intent(inout) :: self  ! the effpot may save the states.
    ! inputs
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    ! outputs
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    integer :: i
    ! Note
    ! calculate force and strain if asked to
    if(present(force)) force(:,:)=0.0d0
    if(present(stress)) stress(:,:)=0.0d0
    if(present(bfield)) bfield(:,:)=0.0d0
    if(present(lwf_force)) lwf_force(:)=0.0d0
    if(present(energy)) energy =0.0
    do i=1, self%size
       call self%list(i)%ptr%calculate(displacement=displacement, strain=strain, &
            & spin=spin, lwf=lwf, force=force, stress=stress, bfield=bfield, &
            lwf_force=lwf_force, energy=energy)
    end do
  end subroutine calculate

  !----------------------------------------------------------------------
  !> @brief get_delta_E: calculate the energy difference when a given spin
  !> is changed. This is to be used for spin Monte Carlo. Currently the
  !> only supported is the spin model.
  !>
  !> @param[in]  S: spin of full structure. array of (3, nspin)
  !> @param[in]  ispin: the index of spin changed. integer
  !> @param[in]  snew: the new value of the changed spin.
  !> @param[out] deltaE: the energy difference
  !----------------------------------------------------------------------
  subroutine get_delta_E(self, S, ispin, Snew, deltaE)
    class(potential_list_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), intent(inout) :: S(:,:),  Snew(:)
    integer, intent(in) :: ispin
    real(dp), intent(inout) :: deltaE
    integer :: i
    do i=1, self%size
       call self%list(i)%ptr%get_delta_E(S=S, ispin=ispin, Snew=Snew, deltaE=deltaE)
    end do
    end subroutine get_delta_E

end module m_potential_list

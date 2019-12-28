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
  use m_xmpi

  use m_mpi_scheduler, only: init_mpi_info
  use m_multibinit_dataset , only: multibinit_dtset_type
  !use m_multibinit_cell, only: mbcell_t
  use m_supercell_maker, only: supercell_maker_t
  use m_abstract_potential, only: abstract_potential_t
  use m_potential_list, only: potential_list_t
  use m_primitive_potential, only: primitive_potential_t
  implicit none
  private
!!***

  !-------------------------------------------------------------------!
  ! primitve_potential_pointer_t
  !-------------------------------------------------------------------!
  type, public:: primitive_potential_pointer_t ! pointer to effpot
     class(primitive_potential_t), pointer :: obj=>null()
  end type primitive_potential_pointer_t


  !-------------------------------------------------------------------!
  ! primitive_potential_list_t
  ! abstract type for primitive potential list
  !-------------------------------------------------------------------!
  type, public, extends(primitive_potential_t):: primitive_potential_list_t
     type(primitive_potential_pointer_t), allocatable :: data(:) ! list of pointer type
     integer :: size=0   ! number of components.
     integer :: capacity=0  ! number of slots allocated for saving the pointers. 
   contains
     procedure :: initialize
     procedure :: append
     procedure :: load_from_files
     procedure :: fill_supercell
     procedure :: fill_supercell_list
     procedure:: finalize
  end type primitive_potential_list_t
contains

  !-------------------------------------------------------------------!
  ! fill supercell: 
  !-------------------------------------------------------------------!
  subroutine fill_supercell(self, scmaker, params, scpot)
    class(primitive_potential_list_t), intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: scmaker
    type(multibinit_dtset_type),       intent(inout) :: params
    class(abstract_potential_t), pointer, intent(inout) :: scpot
    ! Note that sc_pot is a pointer
    ! use a pointer to the specific potential which will be filled
    ! e.g. type(spin_potential_t), pointer :: tmp
    type(abstract_potential_t), pointer :: tmp
    ABI_DATATYPE_ALLOCATE_SCALAR(abstract_potential_t, tmp)
    !call tmp%initialize(....)
    ! set tmp
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(scmaker)
    ABI_UNUSED_A(params)
    ABI_UNUSED_A(scpot)

    nullify(tmp)
  end subroutine fill_supercell

  !-------------------------------------------------------------------!
  ! load primitive potential from file
  !-------------------------------------------------------------------!
  subroutine load_from_files(self, params,  fnames)
    class(primitive_potential_list_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(params)
    ABI_UNUSED_A(fnames)


  end subroutine load_from_files

  !-------------------------------------------------------------------!
  ! save primitive potential to file
  !-------------------------------------------------------------------!
  subroutine save_to_file(self, fname)
    class(primitive_potential_list_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: fname
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(fname)

  end subroutine save_to_file


  !-------------------------------------------------------------------!
  ! Initialize
  !-------------------------------------------------------------------!
  subroutine initialize(self)
    class(primitive_potential_list_t), intent(inout):: self
    self%label="ListPotential"
    self%size=0
    self%capacity=0
  end subroutine initialize

  !----------------------------------------------------------------------
  !> @brief finalize, all the pots in the list will also be finalized.
  !> and the pointers will be nullified.
  !----------------------------------------------------------------------
  subroutine finalize(self)
    class(primitive_potential_list_t), intent(inout):: self
    integer :: i

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    call xmpi_bcast(self%size, master, comm, ierr)
    do i=1, self%size
       call self%data(i)%obj%finalize()
       if(associated(self%data(i)%obj)) then
           ABI_DATATYPE_DEALLOCATE_SCALAR(self%data(i)%obj)
       endif
       nullify(self%data(i)%obj)
    end do
    if(allocated(self%data)) then
       ABI_DEALLOCATE(self%data)
    end if
    nullify(self%primcell)
    self%size=0
    self%capacity=0
    self%has_displacement=.False.
    self%has_strain=.False.
    self%has_spin=.False.
    self%has_lwf=.False.
  end subroutine finalize

  !----------------------------------------------------------------------
  !> @brief append a potential to the list
  !>    The meta data will be updated accordingly.
  !> @param[in]  input
  !> @param[out] output
  !----------------------------------------------------------------------
  subroutine append(self, pot)
    class(primitive_potential_list_t), intent(inout):: self
    class(primitive_potential_t), target, intent(inout) :: pot
    type(primitive_potential_pointer_t), allocatable :: temp(:)
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    self%size=self%size + 1
    if(self%size==1) then
       self%capacity=8
       ABI_ALLOCATE(self%data, (self%capacity))
    else if ( self%size>self%capacity ) then
       self%capacity = self%size + self%size / 4 + 8
       ABI_MALLOC(temp, (self%capacity))
       temp(1:self%size-1) = self%data(:)
       ABI_MOVE_ALLOC(temp, self%data) !temp gets deallocated
    end if

    call xmpi_barrier(comm)
    self%data(self%size)%obj=>pot
    self%has_spin= (self%has_spin .or. pot%has_spin)
    self%has_displacement= (self%has_displacement .or. pot%has_displacement)
    self%has_strain= (self%has_strain.or. pot%has_strain)
    self%has_lwf= (self%has_lwf.or. pot%has_lwf)
    call xmpi_bcast(self%size, master, comm, ierr)
    call xmpi_bcast(self%capacity, master, comm, ierr)
    call xmpi_bcast(self%has_spin, master, comm, ierr)
    call xmpi_bcast(self%has_displacement, master, comm, ierr)
    call xmpi_bcast(self%has_strain, master, comm, ierr)
    call xmpi_bcast(self%has_lwf, master, comm, ierr)

  end subroutine append

  
  !----------------------------------------------------------------------
  !> @brief build supercell potential for every component in the list
  !> Here sc_pot is an pointer.
  !> @param[in]  sc_maker: the helper class for uilder supercell
  !> @param[out] sc_pots: the potential list of supercell pots.
  !----------------------------------------------------------------------
  subroutine fill_supercell_ptr(self, sc_maker, params, sc_pot)
    class(primitive_potential_list_t), intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: sc_maker
    type(multibinit_dtset_type),       intent(inout) :: params
    class(abstract_potential_t), pointer, intent(inout) :: sc_pot

    ! Note that sc_pot is a pointer
    ! use a pointer to the specific potential which will be filled
    type(potential_list_t), pointer :: tmp
    ABI_MALLOC_SCALAR(tmp)
    call self%fill_supercell_list( sc_maker, params, tmp)
    sc_pot=>tmp
    nullify(tmp)
  end subroutine fill_supercell_ptr

  !----------------------------------------------------------------------
  !> @brief build supercell potential for every component in the list
  !>
  !> @param[in]  sc_maker: the helper class for uilder supercell
  !> @param[out] sc_pots: the potential list of supercell pots.
  !----------------------------------------------------------------------
  subroutine fill_supercell_list(self, sc_maker, params, sc_pots)
    class(primitive_potential_list_t), intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: sc_maker
    type(multibinit_dtset_type),       intent(inout) :: params
    type(potential_list_t),            intent(inout) :: sc_pots

    ! Note that sc_pot is a pointer
    ! use a pointer to the specific potential which will be filled
    class(abstract_potential_t), pointer :: tmp
    integer :: i
    tmp=>null()
    do i =1, self%size
      call self%data(i)%obj%fill_supercell(sc_maker, params, tmp)
      call sc_pots%append(tmp)
    end do

    ABI_UNUSED_A(params)

  end subroutine fill_supercell_list

end module m_primitive_potential_list


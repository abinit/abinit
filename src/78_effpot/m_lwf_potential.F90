!!****m* ABINIT/m_lwf_potential
!! NAME
!! m_lwf_potential
!!
!! FUNCTION
!! This module contains an LWF potential. 
!!
!! Datatypes:
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
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


module m_lwf_potential
  use defs_basis
  use m_errors
  use m_abicore
  use m_abstract_potential, only: abstract_potential_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_spmat_coo, only: COO_mat_t
  use m_spmat_csr, only: CSR_mat_t
  use m_spmat_convert, only: spmat_convert
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_hashtable_strval, only: hash_table_t
  use m_twobody_interaction, only: get_twobody_dEdx, get_twobody_delta_E
  implicit none
!!***

  private

  type, public, extends(abstract_potential_t) :: lwf_potential_t
     integer :: nlwf ! number of lwf
     real(dp) :: ref_energy=0.0      ! reference energy
     logical :: csr_mat_ready=.False.
     type(COO_mat_t) :: coeff_coo  ! coefficient. A COO sparse matrix (3N*3N).
     type(CSR_mat_t) :: coeff
     type(NDCOO_mat_t) :: coeff2

     logical :: has_self_bound_term = .False.
     integer :: self_bound_order=0
     real(dp) :: self_bound_coeff=0.0_dp
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_supercell
     procedure :: set_ref_energy
     procedure :: calculate
     procedure :: add_term
     procedure :: convert_coeff_to_csr
     procedure :: get_delta_E_lwf
     procedure :: add_self_bound_term
  end type lwf_potential_t

contains


  !-------------------------------------------------------------------!
  ! initialize
  ! Input:
  !  natom: number of atoms
  !-------------------------------------------------------------------!
  subroutine initialize(self, nlwf)
    class(lwf_potential_t), intent(inout) :: self
    integer, intent(in) :: nlwf
    self%has_lwf= .False.
    self%is_null = .False.
    self%label="lwf_potential"
    self%nlwf=nlwf
    self%ref_energy=0.0_dp
    call self%coeff_coo%initialize(mshape= [self%nlwf, self%nlwf])
    !call self%coeff%initialize(mshape= [self%nlwf, self%nlwf])
    call self%coeff2%initialize(mshape= [self%nlwf, self%nlwf])
    self%csr_mat_ready=.False.
  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(lwf_potential_t), intent(inout) :: self
    self%has_displacement=.False.
    call self%coeff%finalize()
    call self%coeff2%finalize()
    call self%abstract_potential_t%finalize()
    self%is_null=.True.
    self%nlwf=0
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Add a term to the potential
  !-------------------------------------------------------------------!
  subroutine add_term(self, i,j, val)
    class(lwf_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val
    call self%coeff_coo%add_entry([i,j], val)
  end subroutine add_term

  subroutine add_higher_order_term(self, i, j, orderi, orderj, val)
    class(lwf_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j, orderi, orderj
    real(dp), intent(in) :: val
    call self%coeff2%add_entry([i,j, orderi, orderj], val)
  end subroutine add_higher_order_term

  subroutine convert_coeff_to_csr(self)
    class(lwf_potential_t), intent(inout) :: self

    if (.not. self%csr_mat_ready) then
       !call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
       !if(iam_master) then
       call spmat_convert(self%coeff_coo, self%coeff)
       call self%coeff_coo%finalize()
       !endif
       !call self%bilinear_csr_mat%sync(master=master, comm=comm, nblock=1)
       self%csr_mat_ready=.True.
       !call xmpi_bcast(self%csr_mat_ready, master, comm, ierr)
    endif

  end subroutine convert_coeff_to_csr


  !-------------------------------------------------------------------!
  ! Set the reference energy.
  !-------------------------------------------------------------------!
  subroutine set_ref_energy(self, ref_energy)
    class(lwf_potential_t), intent(inout) :: self
    real(dp), intent(in) :: ref_energy
    self%ref_energy=ref_energy
  end subroutine set_ref_energy

  !-------------------------------------------------------------------!
  ! set_supercell
  !  link the supercell with potential.
  ! Inputs:
  !   supercell: mbsupercell_t
  !-------------------------------------------------------------------!
  subroutine set_supercell(self, supercell)
    class(lwf_potential_t), intent(inout) :: self
    type(mbsupercell_t), target, intent(inout) :: supercell
    self%supercell => supercell
  end subroutine set_supercell


  !-------------------------------------------------------------------!
  ! calculate force and energy from harmonic potential
  ! F= - IFC .matmul. displacement
  ! E = 1/2 (-F) .dot. displacement = 1/2<disp|IFC|disp>
  ! Input:
  !   displacement: required.
  !-------------------------------------------------------------------!
  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, energy, energy_table)
    ! This function calculate the energy and its first derivative
    ! the inputs and outputs are optional so that each effpot can adapt to its
    ! own.
    ! In principle, the 1st derivatives are only calculated if asked to (present). However, they can be computed if it is simply convinient to do.
    class(lwf_potential_t), intent(inout) :: self  ! the effpot may save the states.

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    type(hash_table_t),optional, intent(inout) :: energy_table
    real(dp) :: etmp
    real(dp) :: f(self%nlwf)
    ! if present in input
    ! calculate if required
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(force)

    if (.not. present(lwf)) then
       MSG_BUG("lwf should be present in lwf calculate function")
    end if
    if (.not. present(lwf_force)) then
       MSG_BUG("lwf_force should be present in lwf calculate function")
    end if

    etmp=0.0_dp
    call self%convert_coeff_to_csr()
    call self%coeff%mv(lwf, f)

    if (present(lwf_force)) then
       lwf_force(:) = lwf_force(:) - f(:)
       if (self%has_self_bound_term) then
          lwf_force(:) = lwf_force(:) - &
               & self%self_bound_coeff*self%self_bound_order* lwf**(self%self_bound_order-1)
       endif
    endif

    !energy =energy + 0.5_dp * sum(f*displacement)
    etmp=0.5_dp * sum(f*lwf) 
    if (self%has_self_bound_term) then
       etmp = etmp + &
            & self%self_bound_coeff*sum(lwf**(self%self_bound_order))
    end if

    if (present(energy)) then
       energy=energy+etmp
    endif
    if(present(energy_table)) then
       call energy_table%put(self%label, etmp)
    endif
  end subroutine calculate

  !----------------------------------------------------------------------
  !> @brief get_delta_E_lwf: calculate the energy difference when a given lwf
  !> is changed. This is to be used for spin Monte Carlo. Currently the
  !> only supported is the spin model. 
  !>
  !> @param[in]  lwf: lwf of full structure. array of (nlwf)
  !> @param[in]  ilwf: the index of spin changed. integer
  !> @param[in]  lwf_new: the new value of the changed spin. 
  !> @param[out] deltaE: the energy difference
  !----------------------------------------------------------------------
  subroutine get_delta_E_lwf(self, lwf, ilwf, lwf_new, deltaE)
    class(lwf_potential_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp),  intent(inout) :: lwf(:),  lwf_new
    integer,  intent(in) :: ilwf
    real(dp), intent(inout) :: deltaE
    real(dp) :: tmp, dlwf, lold

    lold=lwf(ilwf)
    dlwf=lwf_new-lold
    lwf(ilwf) = lwf_new

    call self%convert_coeff_to_csr()
    call self%coeff%mv_one_row(ilwf, lwf, tmp)
    deltaE=deltaE+tmp*dlwf

    if (self%has_self_bound_term) then
       deltaE= deltaE+ &
            & self%self_bound_coeff*(lwf_new**(self%self_bound_order) &
            & - lold**(self%self_bound_order))
    end if

    lwf(ilwf)=lwf(ilwf)-dlwf
  end subroutine get_delta_E_lwf

  subroutine add_self_bound_term(self, order, coeff)
    class(lwf_potential_t), intent(inout) :: self 
    integer, intent(in) :: order
    real(dp), intent(in) :: coeff
    if (order /= 0) then
       self%has_self_bound_term=.True.
       self%self_bound_order=order
       self%self_bound_coeff=coeff
    end if
  end subroutine add_self_bound_term

end module m_lwf_potential

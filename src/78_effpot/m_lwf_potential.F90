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
  use m_spmat_spvec, only: sp_real_vec
  use m_spmat_convert, only: spmat_convert
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_hashtable_strval, only: hash_table_t
  use m_twobody_interaction, only: get_twobody_dEdx, get_twobody_delta_E
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_lattice_lwf_map, only: lwf_force_to_lattice, lattice_to_lwf_projection
  implicit none
!!***

  private

  type, public, extends(abstract_potential_t) :: lwf_potential_t
     integer :: nlwf ! number of lwf
     real(dp) :: ref_energy=0.0      ! reference energy
     logical :: csr_mat_ready=.False.
     type(COO_mat_t) :: coeff_coo  ! coefficient. A COO sparse matrix (3N*3N).
     type(CSR_mat_t) :: coeff
     type(NDCOO_mat_t) :: onebody_coeff ! onebody anharmonic
     type(NDCOO_mat_t) :: coeff2 ! twobody anharmonic
     type(sp_real_vec), allocatable :: lwf_latt_coeffs(:) ! ( nlwf)

     logical :: use_harmonic = .True.
     logical :: has_self_bound_term = .False.
     integer :: self_bound_order=0
     real(dp) :: self_bound_coeff=0.0_dp

     real(dp) :: beta
     real(dp), allocatable :: coeff_diag(:)

     logical :: as_lattice_anharmonic=.False.
     ! tmp arrays for forces and lwf
     real(dp), allocatable :: lwf_force(:), lwf_amp(:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: use_as_lattice_anharmonic
     procedure :: set_supercell
     procedure :: set_ref_energy
     procedure :: set_params
     procedure :: calculate
     procedure :: add_term
     procedure :: convert_coeff_to_csr
     procedure :: get_delta_E_lwf
     procedure :: add_onebody_term
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
    self%has_lwf= .True.
    self%is_null = .False.
    self%label="lwf_potential"
    self%nlwf=nlwf
    self%ref_energy=0.0_dp
    call self%coeff_coo%initialize(mshape= [self%nlwf, self%nlwf])
    call self%onebody_coeff%initialize(mshape= [self%nlwf, -1])
    call self%coeff2%initialize(mshape= [self%nlwf, self%nlwf])
    self%csr_mat_ready=.False.
    ABI_MALLOC(self%coeff_diag, (self%nlwf))
    ABI_MALLOC(self%lwf_latt_coeffs, (self%nlwf))

    ABI_MALLOC(self%lwf_force, (self%nlwf))
    ABI_MALLOC(self%lwf_amp, (self%nlwf))

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(lwf_potential_t), intent(inout) :: self
    integer :: ilwf
    self%has_displacement=.False.
    call self%onebody_coeff%finalize()
    call self%coeff2%finalize()
    call self%abstract_potential_t%finalize()

    if (.not. self%csr_mat_ready) then
       call self%coeff_coo%finalize()
    else
       call self%coeff%finalize()
    end if
    ABI_SFREE(self%coeff_diag)
    do ilwf=1, self%nlwf
       call self%lwf_latt_coeffs(ilwf)%finalize()
    end do
    ABI_SFREE(self%lwf_latt_coeffs)
    self%nlwf=0
    self%is_null=.True.
    self%as_lattice_anharmonic=.False.
    ABI_SFREE(self%lwf_force)
    ABI_SFREE(self%lwf_amp)
    self%csr_mat_ready=.False.
  end subroutine finalize

  subroutine use_as_lattice_anharmonic(self)
    class(lwf_potential_t), intent(inout) :: self
    self%as_lattice_anharmonic=.True.
    self%use_harmonic = .False.
    self%has_displacement = .True.
    self%has_lwf= .False.
    self%is_null = .False.
    self%label="latt_lwf_potential"
  end subroutine use_as_lattice_anharmonic



  !===============================================================
  !
  !> @
  !===============================================================
  subroutine add_lattice_coeffs(self, ilwf, ilatt, val)
    class(lwf_potential_t), intent(inout) :: self
    integer , intent(in) :: ilwf, ilatt
    real(dp) , intent(in) :: val
    call self%lwf_latt_coeffs(ilwf)%push(ilatt, val)
  end subroutine add_lattice_coeffs

  !-------------------------------------------------------------------!
  ! Add a term to the potential
  !-------------------------------------------------------------------!
  subroutine add_term(self, i,j, val)
    class(lwf_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val
    call self%coeff_coo%add_entry([i,j], val)
  end subroutine add_term

  subroutine add_onebody_term(self, i, order, val)
    class(lwf_potential_t), intent(inout) :: self
    integer, intent(in) :: i, order
    real(dp), intent(in) :: val
    call self%onebody_coeff%add_entry([i, order], val)
  end subroutine add_onebody_term

  subroutine add_higher_order_term(self, i, j, orderi, orderj, val)
    class(lwf_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j, orderi, orderj
    real(dp), intent(in) :: val
    call self%coeff2%add_entry([i,j, orderi, orderj], val)
  end subroutine add_higher_order_term

  subroutine convert_coeff_to_csr(self)
    class(lwf_potential_t), intent(inout) :: self
    integer :: i

    if (.not. self%csr_mat_ready) then
       !call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
       !if(iam_master) then
       call spmat_convert(self%coeff_coo, self%coeff)
       call self%coeff_coo%diag(self%coeff_diag)
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
    real(dp) :: etmp, val
    integer :: i,ilwf, order
    ! if present in input
    ! calculate if required
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)

    !if(.not. present(lwf)) then
    !    MSG_BUG("lwf not exist")
    !end if


    if(self%as_lattice_anharmonic) then
       call lattice_to_lwf_projection(self%lwf_latt_coeffs, displacement, lwf)
    end if


    self%lwf_force(:) =0.0_dp
    etmp=0.0_dp
    ! Harmonic term
    if (self%use_harmonic) then
       call self%convert_coeff_to_csr()
       call self%coeff%mv(lwf, self%lwf_force)
       etmp=etmp+0.5_dp * sum(self%lwf_force*lwf)
       self%lwf_force(:) = -self%lwf_force(:)
    end if

    ! self_bound_term as from the input parameters
    if (self%has_self_bound_term) then
       self%lwf_force(:) = self%lwf_force(:) - &
            & self%self_bound_coeff*self%self_bound_order* lwf**(self%self_bound_order-1)
       etmp = etmp + &
            & self%self_bound_coeff*sum(lwf**(self%self_bound_order))
    endif

    if (self%onebody_coeff%nnz/= 0) then
       do i =1, self%onebody_coeff%nnz
          ilwf= self%onebody_coeff%ind%data(1, i)
          order= self%onebody_coeff%ind%data(2, i)
          val=self%onebody_coeff%val%data(i)
          etmp= etmp + val*lwf(ilwf)**order
          self%lwf_force(ilwf) =self%lwf_force(ilwf) - val*order*lwf(ilwf)**(order-1)
       end do
    end if

    !TODO: remove. For testing only
    !etmp = etmp+ self%beta*sum(lwf(::2)**2 * lwf(1::2)**2)

    if(self%as_lattice_anharmonic) then
       call lwf_force_to_lattice(self%lwf_latt_coeffs, self%lwf_force, force)
    else
       lwf_force=lwf_force+self%lwf_force
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
    real(dp) :: tmp, dlwf, lold, val

    lold=lwf(ilwf)
    dlwf=lwf_new-lold
    lwf(ilwf) = lwf_new

    call self%convert_coeff_to_csr()

    deltaE=0.0_dp
    tmp=0.0_dp
    call self%coeff%mv_one_row(ilwf, lwf, tmp)
    deltaE=deltaE+ tmp*dlwf - 0.5* self%coeff_diag(ilwf)*dlwf*dlwf

    ! bound term
    if (self%has_self_bound_term) then
       deltaE= deltaE+ &
            & self%self_bound_coeff*(lwf_new**(self%self_bound_order) &
            & - lold**(self%self_bound_order))
    ! Adding x^6 and x^8 for VO2
        deltaE=deltaE -1.1344*(lwf_new**6-lold**6) + 0.438*(lwf_new**8-lold**8)
    end if

    ! (Q1 Q2)^2 term
   ! if(modulo(ilwf, 2)==0) then
   !    deltaE=deltaE+ self%beta*((lwf_new**2- lold**2)*lwf(ilwf-1)**2)
   ! else
   !    deltaE=deltaE+ self%beta*((lwf_new**2- lold**2)*lwf(ilwf+1)**2)
   ! endif

    ! if (self%onebody_coeff%nnz/= 0) then
    !    do i =1, self%onebody_coeff%nnz
    !       ilwf= self%onebody_coeff%ind%data(i, 1)
    !       order= self%onebody_coeff%ind%data(i, 2)
    !       val=self%onebody_coeff%val%data(i)
    !       deltaE= deltaE + val*self%lwf_amp(ilwf)**order
    !       self%lwf_force(ilwf) =self%lwf_force(ilwf) - val*order*self%lwf_amp(ilwf)
    !    end do
    ! end if

    if(modulo(ilwf, 2)==0) then
       deltaE=deltaE+ self%beta*((lwf_new**2- lold**2)*lwf(ilwf-1)**2)
    else
       deltaE=deltaE+ self%beta*((lwf_new**2- lold**2)*lwf(ilwf+1)**2)
    endif

    deltaE=deltaE*340.0/1200.0

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

  !----------------------------------------------------------------------
  !> @brief set_params: set the parameters from input file parameters
  !>
  !> @param[in]  params: multibinit_dtset_type: from input file
  !----------------------------------------------------------------------
  subroutine set_params(self, params)
    class(lwf_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(inout) :: params
    self%beta=params%spin_damping
  end subroutine set_params


end module m_lwf_potential

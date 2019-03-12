!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_potential
!! NAME
!! m_spin_potential
!!
!! FUNCTION
!! This module contains the spin hamiltonian, and the methods for
!! calculating effective magnetic field (torque), dS/dt, and total_energy
!!
!!
!! Datatypes:
!!
!! * spin_potential_t
!!
!! Subroutines:
!!
!! * spin_potential_t_initialize
!! * spin_potential_t_finalize
!! * spin_potential_t_total_Heff : calculate total Heff (no Langevin term)
!! * spin_potential_t_Heff_to_dSdt: 
!!  * spin_potential_t_get_dSdt : dSdt, Langevin term is an input.
!!  * spin_potential_t_get_Langevin_Heff

!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (TO, hexu)
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
module  m_spin_potential
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_multibinit_global
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_supercell, only: mb_supercell_t
  use m_spmat_coo, only: coo_mat_t
  use m_spmat_csr, only : CSR_mat_t
  use m_spmat_convert, only : coo_to_csr
  use m_abstract_potential, only : abstract_potential_t
  implicit none
  !!***
  private

  type, public, extends(abstract_potential_t) :: spin_potential_t
     integer :: nspin=0
     logical :: has_external_hfield, has_dipdip, has_bilinear


     ! Array or scalar?
     real(dp), allocatable :: external_hfield(:,:)

     ! Exchange/DMI/dipdip stored like COO sparse matrix form.
     type(coo_mat_t) :: coeff_coo
     logical :: csr_mat_ready= .False.
     type(CSR_mat_t) :: bilinear_csr_mat
     ! 3, 3, ninit

     ! vector for calculating effective field
     real(dp), allocatable :: Htmp(:, :)
     real(dp), allocatable :: ms(:)
   CONTAINS
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     procedure :: set_supercell
     procedure :: get_Heff => spin_potential_t_total_Heff
     procedure :: set_external_hfield
     procedure :: calc_bilinear_term_Heff
     procedure :: calc_external_Heff
     procedure :: calculate => spin_potential_t_calculate
     procedure :: get_energy => spin_potential_t_get_energy
     procedure :: get_delta_E => spin_potential_t_get_delta_E
     procedure :: add_bilinear_term
     procedure :: add_bilinear_term_spin_block
     procedure :: set_bilinear_term
     procedure, private :: prepare_csr_matrix
     procedure :: set_terms
  end type spin_potential_t

contains

  subroutine initialize(self, nspin)
    !Arguments ------------------------------------
    !scalars
    class(spin_potential_t), intent(inout) :: self
    integer :: nspin
    self%label="SpinPotential"
    self%has_spin=.True.
    self%has_displacement=.False.
    self%has_strain=.False.
    self%is_null=.False.
    self%nspin=nspin
    call xmpi_bcast(self%nspin, master, comm, ierr)

    ABI_ALLOCATE( self%ms, (self%nspin))
    if(iam_master) then
       call self%coeff_coo%initialize([self%nspin*3, self%nspin*3])
    endif
    ABI_ALLOCATE(self%external_hfield, (3,self%nspin))
    self%external_hfield=0.0_dp
    ABI_ALLOCATE( self%Htmp, (3, self%nspin))
    self%Htmp(:,:)=0.0_dp

    self%has_external_hfield=.False.
    self%has_dipdip=.False.
    self%has_bilinear=.False.
  end subroutine initialize

  subroutine set_supercell(self, supercell)
    class(spin_potential_t), intent(inout) :: self
    type(mb_supercell_t), target, intent(inout) :: supercell
    self%supercell=>supercell
    self%ms(:)=supercell%ms(:)
    call xmpi_bcast(self%ms, master, comm, ierr)
  end subroutine set_supercell

  subroutine set_terms(self, &
       &     external_hfield, &
       & bilinear_i, bilinear_j, bilinear_val)

    implicit none
    !Arguments ------------------------------------
    !scalars
    !arrays
    class(spin_potential_t), intent(inout) :: self

    ! Terms.
    real(dp), optional, intent(in) :: external_hfield(:,:)
    integer, optional, intent(in) :: bilinear_i(:), bilinear_j(:)
    real(dp), optional,intent(in) :: bilinear_val(:)

    !Local variables-------------------------------
    ! *************************************************************************


    if(present(external_hfield)) then
       call self%set_external_hfield( external_hfield)
    end if

    if ( present(bilinear_i) .and. present( bilinear_j) .and. present(bilinear_val) ) then
       call self%set_bilinear_term(bilinear_i, bilinear_j, bilinear_val)
    endif

  end subroutine set_terms

  !-------------------------------------------------------------------!
  !set_params :
  !-------------------------------------------------------------------!
  subroutine set_params(self, params)
    class(spin_potential_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    real(dp) :: tmp(3, self%nspin)
    integer :: i
    if(iam_master) then
       do i=1, self%nspin
          tmp(:, i)=params%spin_mag_field(:)
       end do
    endif
    call self%set_external_hfield(tmp)
  end subroutine set_params


  subroutine set_external_hfield(self, external_hfield)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(in) :: external_hfield(:,:)
    self%has_external_hfield = .true.
    self%external_hfield = external_hfield
    call xmpi_bcast(self%has_external_hfield, master, comm, ierr)
    call xmpi_bcast(self%external_hfield, master, comm, ierr)
  end subroutine set_external_hfield

  subroutine calc_external_Heff(self, Heff)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(out) :: Heff(:,:)
    Heff(:,:)= Heff(:,:) +self%external_hfield(:,:)
  end subroutine calc_external_Heff

  subroutine add_bilinear_term(self, i,j, val)
    class(spin_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val
    if(iam_master) then
       call self%coeff_coo%add_entry(ind=[i,j],val=val)
       self%has_bilinear=.True.
    endif
  end subroutine add_bilinear_term

  subroutine set_bilinear_term(self, ilist,jlist, vallist)
    class(spin_potential_t), intent(inout) :: self
    integer, intent(in) :: ilist(:), jlist(:)
    real(dp), intent(in) :: vallist(:)
    integer :: i
    if (iam_master) then
       do i = 1, size(ilist)
          call self%coeff_coo%add_entry(ind=[ilist(i),jlist(i)],val=vallist(i))
       end do
    endif
    self%has_bilinear=.True.
  end subroutine set_bilinear_term


  subroutine add_bilinear_term_spin_block(self, ispin, jspin, val)

    class(spin_potential_t), intent(inout) :: self
    integer, intent(in) :: ispin, jspin
    real(dp), intent(in) :: val(:,:)
    integer :: ia, ib
    self%has_bilinear=.True.
    if(iam_master) then
       do ia = 1, 3, 1
          do ib=1, 3, 1
             call self%coeff_coo%add_entry([(ispin-1)*3+ia, &
                  (jspin-1)*3+ib], val=val(ia,ib))
          end do
       end do
    endif
  end subroutine add_bilinear_term_spin_block

  subroutine prepare_csr_matrix(self)
    class(spin_potential_t), intent(inout) :: self
    if (.not. self%csr_mat_ready) then
       if(iam_master) then
          call coo_to_csr(self%coeff_coo, self%bilinear_csr_mat)
          call self%coeff_coo%finalize()
       endif
       call self%bilinear_csr_mat%sync()
       self%csr_mat_ready=.True.
    endif
  end subroutine prepare_csr_matrix


  subroutine calc_bilinear_term_Heff(self, S, Heff)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(inout) :: S(:,:)
    real(dp), intent(out) :: Heff(3,self%nspin)
    integer :: i
    !call self%bilinear_csr_mat%mv(S ,Heff)
    call self%prepare_csr_matrix()
    call self%bilinear_csr_mat%mv_mpi(S ,Heff)
    do i =1, self%nspin
       Heff(:, i)=Heff(:,i)/self%ms(i)*2.0_dp
    end do
  end subroutine calc_bilinear_term_Heff

  subroutine spin_potential_t_calculate(self, displacement, strain, spin, lwf, &
       force, stress, bfield, lwf_force, energy)
    class(spin_potential_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    ! if present in input
    ! calculate if required
    if (present(bfield) .and. present(energy)) then
       call self%get_Heff(spin, bfield, energy)
    else
       call self%get_energy(spin, energy)
    end if
  end subroutine spin_potential_t_calculate


  subroutine spin_potential_t_total_Heff(self,S, Heff, energy)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(inout):: S(3,self%nspin)
    real(dp), intent(inout):: Heff(3,self%nspin)
    real(dp), intent(inout) :: energy
    integer :: i, j

    if(self%has_bilinear) then
       self%Htmp(:,:)=0.0_dp
       call self%calc_bilinear_term_Heff(S,self%Htmp)
       Heff=Heff+self%Htmp
    endif
    if(iam_master) then
       if (self%has_dipdip) then
          continue
          ! TODO implement dipdip and add it
       endif
       ! calculate energy from bilinear terms (all the above ones)
       do i=1, self%nspin
          do j=1, 3
             energy=energy-(Heff(j, i)*S(j,i)*self%ms(i))*0.5_dp
          end do
       end do
       if (self%has_external_hfield) then
          call self%calc_external_Heff(self%Htmp)
          Heff = Heff+self%Htmp
          do i=1, self%nspin
             do j=1, 3
                energy= energy- self%Htmp(j, i)*S(j, i)*self%ms(i)
             end do
          end do
       endif
    endif

  end subroutine spin_potential_t_total_Heff

  subroutine spin_potential_t_get_energy(self, S, energy)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(inout):: S(3,self%nspin)
    real(dp), intent(inout) :: energy
    integer :: i, j
    if(self%has_bilinear) then
       if (.not. self%csr_mat_ready) then
          call coo_to_csr(self%coeff_coo, self%bilinear_csr_mat)
          call self%bilinear_csr_mat%sync()
          self%csr_mat_ready=.True.
       endif
       call self%bilinear_csr_mat%mv_mpi(S ,self%Htmp)
       if(iam_master) then
          energy=energy - sum(sum(self%Htmp* S, dim=1))
       endif
    end if
    if(iam_master) then
       if (self%has_external_hfield) then
          call self%calc_external_Heff(self%Htmp)
          do i=1, self%nspin
             do j=1, 3
                energy= energy- self%Htmp(j, i)*S(j, i)*self%ms(i)
             end do
          end do
       endif
    end if
  end subroutine spin_potential_t_get_energy


  subroutine spin_potential_t_get_delta_E(self, S, ispin, Snew, deltaE)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(inout):: S(:,:), Snew(:)
    integer, intent(in) :: ispin
    real(dp), intent(inout) ::deltaE
    !real(dp) ::  Eold, Enew
    real(dp) :: tmp(3), dS(3)
    ! naive implementation, for test only
    !call self%get_Heff(S, self%Htmp, Eold)
    !call self%bilinear_csr_mat%mv(S, self%Htmp)
    !Eold=-sum(sum(self%Htmp(:,:)*S(:,:), dim=1))

    !Stmp(:,:)=S(:,:)
    !Stmp(:, ispin)= Snew(:)
    !call self%get_Heff(Stmp, self%Htmp, Enew)
    !deltaE=Enew-Eold

    ! test
    dS(:)=Snew(:)-S(:, ispin)
    S(:, ispin)= S(:, ispin)+ dS

    if(self%has_bilinear) then
       if (.not. self%csr_mat_ready) then
          call coo_to_csr(self%coeff_coo, self%bilinear_csr_mat)
          self%csr_mat_ready=.True.
       endif
       call self%bilinear_csr_mat%mv_select_row(3, [3*ispin-2, 3*ispin-1, 3*ispin], S, tmp)
       deltaE=deltaE-dot_product(tmp, dS ) *2.0
    end if

    if(self%has_external_hfield) then
       deltaE=deltaE - dot_product(self%external_hfield(:, ispin), dS)*self%ms(ispin)
    end if
    S(:, ispin)=S(:,ispin)-dS
    !print*, deltaE
  end subroutine spin_potential_t_get_delta_E



  subroutine finalize(self)
    class(spin_potential_t), intent(inout):: self
    if (allocated(self%ms)) then
       ABI_DEALLOCATE(self%ms)
    end if

    if (allocated(self%Htmp)) then
       ABI_DEALLOCATE(self%Htmp)
    endif

    if (allocated(self%external_hfield)) then
       ABI_DEALLOCATE(self%external_hfield)
    endif

    self%has_bilinear=.False.
    ! destroy LIL an CSR
    call self%bilinear_csr_mat%finalize()
    if (.not. self%csr_mat_ready) then
       call self%coeff_coo%finalize()
    end if

  end subroutine finalize

end module m_spin_potential

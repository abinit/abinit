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
!! Copyright (C) 2001-2020 ABINIT group (TO, hexu)
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
  use m_mpi_scheduler, only: mb_mpi_info_t, init_mpi_info, mpi_scheduler_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_spmat_coo, only: coo_mat_t
  use m_spmat_lil, only: lil_mat_t
  use m_spmat_csr, only : CSR_mat_t
  use m_spmat_convert, only : spmat_convert
  use m_abstract_potential, only : abstract_potential_t
  use m_hashtable_strval, only: hash_table_t
  implicit none
  !!***
  private

  type, public, extends(abstract_potential_t) :: spin_potential_t
     integer :: nspin=0
     logical :: has_external_hfield, has_dipdip


     ! Array or scalar?
     real(dp), allocatable :: external_hfield(:,:)

     ! Exchange/DMI/dipdip stored like COO sparse matrix form.
     type(lil_mat_t) :: coeff_coo
     logical :: csr_mat_ready= .False.
     type(CSR_mat_t) :: bilinear_csr_mat
     ! 3, 3, ninit

     ! vector for calculating effective field
     real(dp), allocatable :: Htmp(:, :)
     real(dp), allocatable :: ms(:)
     type(mb_mpi_info_t) :: mpiinfo
     type(mpi_scheduler_t) :: mps
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
     !procedure :: get_energy => spin_potential_t_get_energy
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

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc)
    call self%mps%initialize(ntasks=nspin, master=master, comm=comm)
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
  end subroutine initialize

  subroutine set_supercell(self, supercell)
    class(spin_potential_t), intent(inout) :: self
    type(mbsupercell_t), target, intent(inout) :: supercell
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
 
    self%supercell=>supercell
    self%ms(:)=supercell%spin%ms(:)
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
    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
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
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    self%has_external_hfield = .true.
    self%external_hfield = external_hfield
    call xmpi_bcast(self%has_external_hfield, master, comm, ierr)
    call xmpi_bcast(self%external_hfield, master, comm, ierr)
  end subroutine set_external_hfield

  subroutine calc_external_Heff(self, Heff)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(out) :: Heff(:,:)
    integer :: i
     do i= self%mps%istart, self%mps%iend
       Heff(:, i)=Heff(:,i) + self%external_hfield(:, i)
    end do
  end subroutine calc_external_Heff

  subroutine add_bilinear_term(self, i,j, val)
    class(spin_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val
    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if(iam_master) then
       call self%coeff_coo%add_entry(ind=[i,j],val=val)
    endif
  end subroutine add_bilinear_term

  subroutine set_bilinear_term(self, ilist,jlist, vallist)
    class(spin_potential_t), intent(inout) :: self
    integer, intent(in) :: ilist(:), jlist(:)
    real(dp), intent(in) :: vallist(:)
    integer :: i
    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if (iam_master) then
       do i = 1, size(ilist)
          call self%coeff_coo%add_entry(ind=[ilist(i),jlist(i)],val=vallist(i))
       end do
    endif
  end subroutine set_bilinear_term


  subroutine add_bilinear_term_spin_block(self, ispin, jspin, val)

    class(spin_potential_t), intent(inout) :: self
    integer, intent(in) :: ispin, jspin
    real(dp), intent(in) :: val(:,:)
    integer :: ia, ib

    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

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
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master

    if (.not. self%csr_mat_ready) then
       call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
       if(iam_master) then
          call spmat_convert(self%coeff_coo, self%bilinear_csr_mat)
          call self%coeff_coo%finalize()
       endif
       call self%bilinear_csr_mat%sync(master=master, comm=comm, nblock=1)
       self%csr_mat_ready=.True.
       call xmpi_bcast(self%csr_mat_ready, master, comm, ierr)
    endif
  end subroutine prepare_csr_matrix


  subroutine calc_bilinear_term_Heff(self, S, Heff)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(inout) :: S(:,:)
    real(dp), intent(out) :: Heff(3,self%nspin)
    integer :: i
    !call self%bilinear_csr_mat%mv(S ,Heff)
    call self%prepare_csr_matrix()
    call self%bilinear_csr_mat%mv_mpi(x=S ,b=Heff, bcastx=.False., syncb=.False.)
    do i= self%mps%istart, self%mps%iend
       Heff(:, i)=Heff(:,i)/self%ms(i)*2.0_dp
    end do
  end subroutine calc_bilinear_term_Heff

  subroutine spin_potential_t_calculate(self, displacement, strain, spin, lwf, &
       force, stress, bfield, lwf_force, energy, energy_table)
    class(spin_potential_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    type(hash_table_t),optional, intent(inout) :: energy_table
    real(dp) :: etmp
    ! if present in input
    ! calculate if required

    if (present(bfield) ) then
       call self%get_Heff(spin, bfield, etmp)

       ! only update energy when bfield is asked for.
       if ( present(energy)) then
          energy=energy+etmp
       end if
       if (present(energy_table)) then
          call energy_table%put(self%label, etmp)
       end if
    end if


    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(force)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)
    ABI_UNUSED_A(lwf_force)
    ABI_UNUSED_A(energy)
  end subroutine spin_potential_t_calculate


  subroutine spin_potential_t_total_Heff(self,S, Heff, energy)
    class(spin_potential_t), intent(inout) :: self
    real(dp), intent(inout):: S(3,self%nspin)
    real(dp), intent(inout):: Heff(3,self%nspin)
    real(dp), intent(inout) :: energy
    real(dp) :: etmp
    integer :: i, j
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    self%Htmp(:,:)=0.0_dp
    etmp=0.0_dp


    ! calculate energy from bilinear terms (all the above ones)
    call self%calc_bilinear_term_Heff(S,self%Htmp)
    do i= self%mps%istart, self%mps%iend
       Heff(:,i)=Heff(:,i)+self%Htmp(:,i)
       do j=1, 3
          etmp=etmp-(self%Htmp(j, i)*S(j,i)*self%ms(i))*0.5_dp
       end do
    enddo
    if(iam_master) then
       if (self%has_dipdip) then
          continue
          ! TODO implement dipdip and add it
       endif
    endif

    ! linear terms
    if (self%has_external_hfield) then
       call self%calc_external_Heff(self%Htmp)
        do i= self%mps%istart, self%mps%iend
           Heff(:,i)=Heff(:,i)+self%Htmp(:,i)
           do j=1, 3
              etmp=etmp-(self%Htmp(j, i)*S(j,i)*self%ms(i))
           end do
        enddo
    endif

    call xmpi_sum_master(etmp, 0, xmpi_world, ierr )
    ! NOTE: here energy is not added to input energy.
    energy=etmp

  end subroutine spin_potential_t_total_Heff

!  subroutine spin_potential_t_get_energy(self, S, energy)
!    class(spin_potential_t), intent(inout) :: self
!    real(dp), intent(inout):: S(3,self%nspin)
!    real(dp), intent(inout) :: energy
!    integer :: i, j
!    integer :: master, my_rank, comm, nproc, ierr
!    logical :: iam_master
!    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
!
!
!    if (.not. self%csr_mat_ready) then
!       call spmat_convert(self%coeff_coo, self%bilinear_csr_mat)
!       call self%bilinear_csr_mat%sync()
!       self%csr_mat_ready=.True.
!    endif
!    call self%bilinear_csr_mat%mv_mpi(S ,self%Htmp)
!    if(iam_master) then
!       energy=energy - sum(sum(self%Htmp* S, dim=1))
!    endif
!    if(iam_master) then
!       if (self%has_external_hfield) then
!          call self%calc_external_Heff(self%Htmp)
!          do i=1, self%nspin
!             do j=1, 3
!                energy= energy- self%Htmp(j, i)*S(j, i)*self%ms(i)
!             end do
!          end do
!       endif
!    end if
!  end subroutine spin_potential_t_get_energy


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

       if (.not. self%csr_mat_ready) then
          call spmat_convert(self%coeff_coo, self%bilinear_csr_mat)
          self%csr_mat_ready=.True.
       endif
       call self%bilinear_csr_mat%mv_select_row(3, [3*ispin-2, 3*ispin-1, 3*ispin], S, tmp)
       deltaE=deltaE-dot_product(tmp, dS ) *2.0

    if(self%has_external_hfield) then
       deltaE=deltaE - dot_product(self%external_hfield(:, ispin), dS)*self%ms(ispin)
    end if
    S(:, ispin)=S(:,ispin)-dS
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

    call self%bilinear_csr_mat%finalize()
    if (.not. self%csr_mat_ready) then
       call self%coeff_coo%finalize()
    end if

    call self%mps%finalize()

  end subroutine finalize

end module m_spin_potential

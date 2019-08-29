!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_slc_potential
!! NAME
!! m_slc_potential
!!
!! FUNCTION
!! This module contains the spin-lattice coupling, and the methods for
!! calculating effective magnetic field (torque), force, and total_energy
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (TO, hexu, NH)
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
module  m_slc_potential
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi

  use m_abstract_potential, only : abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_mpi_scheduler, only: mb_mpi_info_t, init_mpi_info, mpi_scheduler_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_spmat_convert, only : spmat_convert
  use m_spmat_coo, only: coo_mat_t
  use m_spmat_csr, only : CSR_mat_t
  use m_spmat_lil, only: lil_mat_t
  use m_spmat_ndcoo, only: ndcoo_mat_t

  implicit none
  private

  type, public, extends(abstract_potential_t) :: slc_potential_t
     integer :: nspin=0
     integer :: natom=0
     logical :: has_bilin=.False.   ! bilinear coupling term, i.e. liu        
     logical :: has_linquad=.False. ! spin first then lattice, i.e. niuv
     logical :: has_quadlin=.False. ! spin first then lattice, i.e. oiju
     logical :: has_biquad=.False.  ! biquadratic coupling term, i.e. tijuv

     type(ndcoo_mat_t) :: liu_sc           ! parameter values bilin term
     type(ndcoo_mat_t) :: niuv_sc          ! parameter values linquad term
     type(ndcoo_mat_t) :: oiju_sc          ! parameter values quadlin term
     type(ndcoo_mat_t) :: tijuv_sc         ! parameter values biquad term

     ! magnetic moments
     real(dp), allocatable :: ms(:)

     ! mpi
     type(mb_mpi_info_t) :: mpiinfo
     type(mpi_scheduler_t) :: mpsspin
     type(mpi_scheduler_t) :: mpslatt

   CONTAINS
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     procedure :: set_supercell
     procedure :: add_oiju_term
     procedure :: add_tijuv_term
     procedure :: calculate
  end type slc_potential_t
 
contains

  subroutine initialize(self, nspin, natom)
    !Arguments ------------------------------------
    !scalars
    class(slc_potential_t), intent(inout) :: self
    integer,                intent(in)    :: nspin
    integer,                intent(in)    :: natom

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc)
    call self%mpsspin%initialize(ntasks=nspin, master=master, comm=comm)
    call self%mpslatt%initialize(ntasks=natom, master=master, comm=comm)
    self%label="SLCPotential"
    self%has_spin=.True.
    self%has_displacement=.True.
    self%has_strain=.False.
    self%is_null=.False.
    self%nspin=nspin
    self%natom=natom
    
    call xmpi_bcast(self%nspin, master, comm, ierr)
    call xmpi_bcast(self%natom, master, comm, ierr)
    ABI_ALLOCATE(self%ms, (self%nspin))
  end subroutine initialize

  subroutine finalize(self)

    class(slc_potential_t), intent(inout):: self

    if(self%has_bilin) call self%liu_sc%finalize()
    if(self%has_quadlin) call self%oiju_sc%finalize()
    !if(self%has_linquad) call self%niuv%finalize() 
    !if(self%has_biquad) call self%tijuv%finalize()

    call self%mpsspin%finalize()
    call self%mpslatt%finalize()

  end subroutine finalize

  !-------------------------------------------------------------------!
  !set_params: which coupling terms are used
  !-------------------------------------------------------------------!
  subroutine set_params(self, params)
    class(slc_potential_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params

    integer :: master, my_rank, comm, nproc, ierr, coupling
    logical :: iam_master

    coupling = params%slc_coupling

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      if(coupling .ge. 1000) then
        self%has_biquad=.True.
        call xmpi_bcast(self%has_biquad, master, comm, ierr)
        coupling=coupling - 1000
      endif

      if(coupling .ge. 100) then
        self%has_linquad=.True.
        call xmpi_bcast(self%has_linquad, master, comm, ierr)
        coupling=coupling - 100
      endif

      if(coupling .ge. 10) then
        self%has_quadlin=.True.
        call xmpi_bcast(self%has_quadlin, master, comm, ierr)
        coupling=coupling - 10
      endif

      if(coupling .ge. 1) then   
        self%has_bilin=.True.
        call xmpi_bcast(self%has_bilin, master, comm, ierr)
      endif
    endif
  end subroutine set_params

  !-------------------------------------------------------------------!
  !set_supercell: use the same supercell for all terms
  !               copy magnetic moments
  !-------------------------------------------------------------------!
  subroutine set_supercell(self, supercell)
    class(slc_potential_t),      intent(inout) :: self
    type(mbsupercell_t), target, intent(inout) :: supercell
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    self%supercell=>supercell
    self%ms(:)=supercell%spin%ms(:)

    call xmpi_bcast(self%ms, master, comm, ierr)
  end subroutine set_supercell


  subroutine add_oiju_term(self, i,j,u, val)
    class(slc_potential_t), intent(inout) :: self
    integer,                intent(in)    :: i, j, u
    real(dp),               intent(in)    :: val

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if(iam_master) then
       call self%oiju_sc%add_entry(ind=[i,j,u],val=val)
    endif
  end subroutine add_oiju_term

  subroutine add_tijuv_term(self, i,j,u,v, val)
    class(slc_potential_t), intent(inout) :: self
    integer,                intent(in)    :: i, j, u, v
    real(dp),               intent(in)    :: val

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if(iam_master) then
       call self%tijuv_sc%add_entry(ind=[i,j,u,v],val=val)
    endif
  end subroutine add_tijuv_term


  subroutine calculate(self, displacement, strain, spin, lwf, &
       force, stress, bfield, lwf_force, energy)
    class(slc_potential_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
 
    integer :: ii
    real(dp) :: f1(1:3*self%natom), disp(1:3*self%natom), b1(1:3*self%nspin), sp(1:3*self%nspin)
    real(dp) :: btmp(3, self%nspin)

    sp(:) = reshape(spin, (/ 3*self%nspin /))
    disp(:) = reshape(displacement, (/ 3*self%natom /))

    write(*,*) 'Calculating coupling terms'

    ! oiju contribution
    if(present(bfield)) then
      b1(:) = 0.0d0
      call self%oiju_sc%vec_product(1, sp, 3, disp, 2, b1)
      btmp = reshape(b1, (/ 3, self%nspin /))
      do ii = 1, self%nspin
        btmp(:, ii) = btmp(:,ii)/self%ms(ii)
      enddo
      bfield(:,:) = bfield(:,:) + btmp(:,:)
      write(*,*) 'Magnetic fields are'
      do ii = 1, self%nspin
        !if(dot_product(bfield(:,ii), bfield(:,ii)).gt.1d-16) then
        !  write(*,*) ii, bfield(:,ii)
        !endif
      enddo
    endif

    if(present(force)) then
      f1(:) = 0.0d0
      call self%oiju_sc%vec_product(1, sp, 2, sp, 3, f1)
      force(:,:) = force(:,:) + reshape(f1, (/ 3, self%natom /))
      write(*,*) 'Forces are'
      do ii = 1, self%natom
        !if(dot_product(force(:,ii), force(:,ii)).gt.1d-16) then
        !  write(*,*) ii, force(:,ii)
        !endif
      enddo
    endif

    if(present(energy)) then
      if(present(bfield)) then
        b1=reshape(btmp, (/ 3*self%nspin /))
        energy = energy + dot_product(b1, sp)
      else 
        if(present(force)) then
          energy = energy - 0.5_dp*dot_product(f1, disp)
        else
          f1(:) = 0.0d0
          call self%oiju_sc%vec_product(1, sp, 2, sp, 3, f1)
          energy = energy + dot_product(f1, disp)
        endif
      endif   
    endif
 
  end subroutine calculate


end module m_slc_potential

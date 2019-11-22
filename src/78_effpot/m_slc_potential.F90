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

  use m_hashtable_strval, only: hash_table_t
  use m_abstract_potential, only : abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_mpi_scheduler, only: mb_mpi_info_t, init_mpi_info, mpi_scheduler_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_multibinit_dataset, only: multibinit_dtset_type
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

     ! precalculated things for reference structure
     real(dp), allocatable :: fref(:)        ! force from liu and oiju for reference spin structure
     type(ndcoo_mat_t) :: matrixref    ! matrix in u and v from niuv and tijuv for reference spin structure

     ! mpi
     type(mb_mpi_info_t) :: mpiinfo
     type(mpi_scheduler_t) :: mpsspin
     type(mpi_scheduler_t) :: mpslatt

   CONTAINS
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     procedure :: set_supercell
     procedure :: calculate_ref
     procedure :: add_liu_term
     procedure :: add_niuv_term
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
    if(self%has_linquad) call self%niuv_sc%finalize() 
    if(self%has_biquad) call self%tijuv_sc%finalize()

    call self%mpsspin%finalize()
    call self%mpslatt%finalize()
    ABI_SFREE(self%fref)
    if(self%has_linquad .or. self%has_biquad) call self%matrixref%finalize()

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

  !-------------------------------------------------------------------!
  !calculate ref: calculates terms necessary for the force and energy
  !               of the reference structure terms
  !-------------------------------------------------------------------!

  subroutine calculate_ref(self)
    class(slc_potential_t), intent(inout) :: self

    integer :: ii
    real(dp) :: spref(1:3*self%nspin), beta
    real(dp), allocatable :: force(:)
    type(ndcoo_mat_t) :: m2dim

    integer :: master, my_rank, comm, nproc
    logical :: iam_master


    spref(:) = reshape(self%supercell%spin%Sref, (/ 3*self%nspin/))

    beta = 0.5_dp

    if(self%has_bilin .or. self%has_quadlin) then
      ABI_ALLOCATE(self%fref, (3*self%natom))
      ABI_ALLOCATE(force, (3*self%natom))
      self%fref=0.0d0
      if(self%has_bilin) then
        force = 0.0d0
        call self%liu_sc%vec_product2d(1, spref, 2, force)
        self%fref(:) = self%fref(:) - force(:)
      endif
      if(self%has_quadlin) then
        force = 0.0d0
        call self%oiju_sc%vec_product(1, spref, 2, spref, 3, force)
        self%fref(:) = self%fref(:) - beta*force(:)
      endif
      ABI_SFREE(force)
    endif

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(self%has_linquad .or. self%has_biquad) then
      if(iam_master) then
        call self%matrixref%initialize(mshape=[self%natom*3, self%natom*3])
      endif
  
      if(self%has_linquad .and. self%has_biquad) then
        call self%niuv_sc%mv1vec(spref, 1, self%matrixref)
        if(iam_master) then
          call m2dim%initialize(mshape=[self%natom*3, self%natom*3])
        endif
        call self%tijuv_sc%mv2vec(0.5*spref, spref, 1, 2, m2dim)
        do ii = 1, m2dim%nnz
          call self%matrixref%add_entry(m2dim%ind%data(:,ii), m2dim%val%data(ii))
        enddo
        call self%matrixref%sum_duplicates()
        if(iam_master) then
          call m2dim%finalize()
        endif
      else
        if(self%has_linquad) call self%niuv_sc%mv1vec(spref, 1, self%matrixref)
        if(self%has_biquad)  call self%tijuv_sc%mv2vec(0.5*spref, spref, 1, 2, self%matrixref)
      endif
    endif

  end subroutine calculate_ref

  !-----------------------------------------------------
  ! add different coupling terms to supercell potential
  ! TODO: test liu, niuv, tijuv
  !-----------------------------------------------------
  subroutine add_liu_term(self, i, u, val)
    class(slc_potential_t), intent(inout) :: self
    integer,                intent(in)    :: i, u
    real(dp),               intent(in)    :: val

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if(iam_master) then
       call self%liu_sc%add_entry(ind=[i,u],val=val)
    endif
  end subroutine add_liu_term


  subroutine add_niuv_term(self, i, u, v,  val)
    class(slc_potential_t), intent(inout) :: self
    integer,                intent(in)    :: i, u, v
    real(dp),               intent(in)    :: val

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if(iam_master) then
       call self%niuv_sc%add_entry(ind=[i,u,v],val=val)
    endif
  end subroutine add_niuv_term



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


  !-----------------------------------------------------------------
  ! Calculate forces, magnetic fields and energy for coupling terms
  ! TODO: precalculate terms containing Sref?
  !-----------------------------------------------------------------

  subroutine calculate(self, displacement, strain, spin, lwf, &
       force, stress, bfield, lwf_force, energy, energy_table)
    class(slc_potential_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    type(hash_table_t),optional, intent(inout) :: energy_table
 
    integer :: ii
    real(dp) :: eslc, beta
    real(dp) :: disp(1:3*self%natom), sp(1:3*self%nspin), spref(1:3*self%nspin) 
    real(dp) :: f1(1:3*self%natom), b1(1:3*self%nspin)
    real(dp) :: btmp(3, self%nspin), bslc(1:3*self%nspin), fslc(1:3*self%natom)

    ABI_UNUSED(lwf)
    ABI_UNUSED(strain)
    ABI_UNUSED(lwf_force)
    ABI_UNUSED(stress)

    sp(:) = reshape(spin, (/ 3*self%nspin /))
    spref(:) = reshape(self%supercell%spin%Sref, (/ 3*self%nspin/))
    disp(:) = reshape(displacement, (/ 3*self%natom /))

    beta = 0.5_dp

    ! Magnetic field
    if(present(bfield)) then
      bslc(:) = 0.0d0
      if(self%has_bilin) then
        b1(:) = 0.0d0
        call self%liu_sc%vec_product2d(2, disp, 1, b1)
        bslc(:) = bslc(:) + b1(:)
      endif      
      if(self%has_linquad) then
        b1(:) = 0.0d0
        call self%niuv_sc%vec_product(2, disp, 3, disp, 1, b1)
        bslc(:) = bslc(:) + beta*b1(:)
      endif      
      if(self%has_quadlin) then
        b1(:) = 0.0d0
        call self%oiju_sc%vec_product(1, sp, 3, disp, 2, b1)
        bslc(:) = bslc(:) + 2.0d0*beta*b1(:)
      endif
      if(self%has_biquad) then
        b1(:) = 0.0d0
        call self%tijuv_sc%vec_product4d(1, sp, 3, disp, 4, disp, 2, b1)
        bslc(:) = bslc(:) + beta*b1(:)
      endif
      btmp = reshape(bslc, (/ 3, self%nspin /))
      do ii = 1, self%nspin
        btmp(:,ii) = btmp(:,ii)/self%ms(ii)
      enddo
      bfield(:,:) = bfield(:,:) + btmp(:,:)

      ! TESTING: write magnetic fields to a file
      write(201,*) 'Magnetic fields are'
      do ii = 1, self%nspin
        !if(dot_product(bfield(:,ii), bfield(:,ii)).gt.1d-16) then
          write(201,*) ii, bfield(:,ii)
        !endif
      enddo
    endif

    ! Force and energy
    if(present(force) .or. present(energy)) then
      fslc(:) = 0.0d0
      eslc = 0.0d0
      if(self%has_bilin) then
        f1(:) = 0.0d0
        call self%liu_sc%vec_product2d(1, sp, 2, f1)
        fslc(:) = fslc(:) + f1(:)
        eslc = eslc - dot_product(f1, disp)
      endif      
      if(self%has_linquad) then
        f1(:) = 0.0d0
        call self%niuv_sc%vec_product(1, sp, 2, disp, 3, f1)
        fslc(:) = fslc(:) + 2.0d0*beta*f1(:)
        eslc = eslc - beta*dot_product(f1, disp)
      endif      
      if(self%has_quadlin) then
        f1(:) = 0.0d0
        call self%oiju_sc%vec_product(1, sp, 2, sp, 3, f1)
        fslc(:) = fslc(:) + beta*f1(:)
        eslc = eslc - beta*dot_product(f1, disp)
      endif
      if(self%has_biquad) then
        f1(:) = 0.0d0
        call self%tijuv_sc%vec_product4d(1, sp, 2, sp, 3, disp, 4, f1)
        fslc(:) = fslc(:) + beta*f1(:)
        eslc = eslc - 0.5_dp*beta*dot_product(f1, disp)
      endif

      ! add forces and energy for reference spin structure terms
      if(self%has_bilin .or. self%has_quadlin) then
        fslc(:) = fslc(:) + self%fref(:)
        eslc = eslc - dot_product(self%fref, disp)
      endif
      if(self%has_linquad .or. self%has_biquad) then
        f1(:) = 0.0d0
        call self%matrixref%vec_product2d(1, disp, 2, f1)
        fslc(:)= fslc+2.0*beta*f1(:)
        eslc = eslc + beta*dot_product(f1, disp)
      endif
    endif !energy or force

    if(present(force)) then
      force(:,:) = force(:,:) + reshape(fslc, (/3, self%natom /))
      !TESTING write forces to file
      write(200,*) 'Forces are'
      do ii = 1, self%natom
        !if(dot_product(force(:,ii), force(:,ii)).gt.1d-16) then
          write(200,*) ii, force(:,ii)
        !endif
      enddo
    endif

    !if(present(energy)) then
    !   energy= eslc
    !end if

    if(present(energy_table)) then
       call energy_table%put(self%label, eslc)
       ! TODO: use terms instead of total, e.g.
       !call energy_table%put("SLC Oiju term", eOiju)
       !call energy_table%put("SLC Tijuv term", eTijuv)
    end if

    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(lwf_force)
    
  end subroutine calculate

  !!***
end module m_slc_potential

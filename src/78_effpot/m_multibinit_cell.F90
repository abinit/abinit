!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_cell
!! NAME
!! m_multibinit_cell
!!
!! FUNCTION
!! This module define the m_unitcell type, which defines an atomic structure, with spin, lwf, electron
!!
!! Datatypes:
!!  
!!
!! Subroutines:
!!
!! COPYRIGHT !! Copyright (C) 2001-2019 ABINIT group (hexu) !! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"


module m_multibinit_cell
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_mpi_scheduler, only: init_mpi_info
  !use m_multibinit_supercell, only: mb_supercell_t
  use m_supercell_maker , only: supercell_maker_t
  use m_multibinit_dataset, only : multibinit_dtset_type
  use m_crystal, only : crystal_t
  implicit none

!!***
  private
  ! TODO: use crystal_t
  !type, public :: mbcell_lattice_t
  !type, public, extends(crystal_t) :: mbcell_lattice_t
  type, public ::mbcell_lattice_t
     integer:: natom
     integer, allocatable  :: zion(:)
     real(dp), allocatable :: masses(:), xcart(:,:)
     real(dp) :: cell(3,3)
   contains
     procedure :: initialize => latt_initialize
     procedure :: finalize => latt_finalize
     procedure :: fill_supercell=>latt_fill_supercell
     procedure :: from_unitcell => latt_from_unitcell
  end type mbcell_lattice_t

  type, public :: mbcell_spin_t
     integer :: nspin
     real(dp) :: rprimd(3,3)
     real(dp), allocatable :: ms(:)
     real(dp), allocatable :: gyro_ratio(:)
     real(dp), allocatable :: gilbert_damping(:)
     real(dp), allocatable :: spin_positions(:,:)
     integer, allocatable ::  ispin_prim(:), rvec(:,:)
   contains
     procedure :: initialize => spin_initialize
     procedure :: set => spin_set
     procedure :: finalize => spin_finalize
     procedure :: fill_supercell=>spin_fill_supercell
     procedure :: from_unitcell=>spin_from_unitcell
  end type mbcell_spin_t

  type, public :: mbcell_lwf_t
     integer :: nlwf
   contains
     procedure :: initialize => lwf_initialize
     procedure :: finalize => lwf_finalize
     procedure :: fill_supercell => lwf_fill_supercell
     procedure :: from_unitcell => lwf_from_unitcell
  end type mbcell_lwf_t

  type ,public :: mbcell_t
     logical :: has_lattice=.False.
     logical :: has_spin=.False.
     logical :: has_lwf=.False.
     logical :: has_electron=.False.
     type(mbcell_lattice_t) :: lattice
     type(mbcell_spin_t) :: spin
     type(mbcell_lwf_t) :: lwf
   contains
     procedure:: initialize
     procedure :: finalize
     procedure :: set_lattice
     procedure :: set_spin
     procedure :: fill_supercell
     !procedure :: from_unitcell
  end type mbcell_t

  type, public, extends(mbcell_t):: mbsupercell_t
     integer :: sc_matrix(3,3)
     integer :: ncell
     type(supercell_maker_t), pointer :: supercell_maker
     type(mbcell_t), pointer :: unitcell
   contains
     procedure :: from_unitcell
  end type mbsupercell_t

contains

  subroutine initialize(self)
    class(mbcell_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine initialize

  subroutine set_lattice(self, natom, cell, xcart, masses, zion)
    class(mbcell_t), intent(inout) :: self
    integer, intent(in) :: natom, zion(:)
    real(dp), intent(in) :: cell(3,3), xcart(:,:), masses(:)
    self%has_lattice=.True.
    call self%lattice%initialize(natom=natom, cell=cell, &
         &xcart=xcart, masses=masses, zion=zion)
  end subroutine set_lattice


  subroutine set_spin(self,nspin, ms, spin_positions, gyro_ratio, gilbert_damping, rvec,  ispin_prim)
    class(mbcell_t) , intent(inout):: self
    integer, intent(in) :: nspin
    real(dp), intent(in) :: ms(nspin), spin_positions(3, nspin), gyro_ratio(nspin), gilbert_damping(nspin)
    integer, optional, intent(in) :: rvec(3, nspin), ispin_prim(nspin)
    self%has_spin=.True.
    call self%spin%initialize(nspin)
    call self%spin%set(nspin, ms, spin_positions, gyro_ratio, gilbert_damping, rvec, ispin_prim)
  end subroutine set_spin

  subroutine set_lwf(self)
    class(mbcell_t) , intent(inout):: self
    self%has_lwf=.True.
    call self%lwf%initialize()
  end subroutine set_lwf

  subroutine read_from_file(self, params, fnames)
    class(mbcell_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(17)

    ABI_UNUSED_A(self)
    ABI_UNUSED_A(params)
    ABI_UNUSED_A(fnames)
  end subroutine read_from_file


  subroutine finalize(self)
    class(mbcell_t), intent(inout) :: self
    call self%lattice%finalize()
    call self%spin%finalize()
    call self%lwf%finalize()
  end subroutine finalize

  subroutine fill_supercell(self, sc_maker, supercell)
    class(mbcell_t), target, intent(inout) :: self
    type(supercell_maker_t), target, intent(inout) :: sc_maker
    class(mbsupercell_t), intent(inout) :: supercell
    call supercell%initialize()
    supercell%ncell=sc_maker%ncells
    supercell%sc_matrix=sc_maker%scmat
    supercell%has_lattice = self%has_lattice
    supercell%has_spin = self%has_spin
    supercell%has_lwf =self%has_lwf
    supercell%supercell_maker => sc_maker
    supercell%unitcell => self
    if (self%has_lattice) then
       call self%lattice%fill_supercell(sc_maker, supercell%lattice)
    endif

    if (self%has_spin) then
       call self%spin%fill_supercell(sc_maker, supercell%spin)
    endif

    if (self%has_lwf) then
       call self%lwf%fill_supercell(sc_maker, supercell%lwf)
    endif

  end subroutine fill_supercell


!=================== Supercell =====================================
  subroutine supercell_initialize(self)
    class(mbsupercell_t), intent(inout) :: self
    call self%mbcell_t%initialize()
  end subroutine supercell_initialize

  subroutine from_unitcell(self, sc_maker, unitcell)
    class(mbsupercell_t), intent(inout) :: self
    type(supercell_maker_t),  intent(inout) :: sc_maker
    class(mbcell_t), intent(inout) :: unitcell
    call unitcell%fill_supercell(sc_maker, self)
  end subroutine from_unitcell

  subroutine supercell_finalize(self)
    class(mbsupercell_t), intent(inout) :: self
    call self%mbcell_t%finalize()
    nullify(self%supercell_maker)
    nullify(self%unitcell)
  end subroutine supercell_finalize

!================================Lattice====================================

  subroutine latt_initialize(self, natom, cell, xcart, masses, zion)
    class(mbcell_lattice_t), intent(inout) :: self
    integer, intent(in) :: natom, zion(:)
    real(dp), intent(in) :: cell(3,3), xcart(:,:), masses(:)
    self%natom=natom
    ABI_ALLOCATE(self%zion, (natom))
    ABI_ALLOCATE(self%xcart, (3, natom))
    ABI_ALLOCATE(self%masses, (natom))

    self%zion(:) = zion(:)
    self%cell(:,:) =cell(:,:)
    self%xcart(:,:) = xcart(:,:)
    self%masses(:) = masses(:)
  end subroutine latt_initialize

  subroutine latt_finalize(self)
    class(mbcell_lattice_t) :: self
    self%natom=0
    if (allocated(self%xcart)) then
       ABI_DEALLOCATE(self%xcart)
    endif
    if (allocated(self%masses)) then
       ABI_DEALLOCATE(self%masses)
    endif
    if (allocated(self%zion)) then
       ABI_DEALLOCATE(self%zion)
    endif
  end subroutine latt_finalize


  !-------------------------------------------------------------------!
  ! latt_fill_supercell: make a lattice supercell from primitive cell
  ! Inputs:
  !> sc_maker: supercell maker
  ! Output:
  !> supercell: supercell
  !-------------------------------------------------------------------!
  subroutine latt_fill_supercell(self, sc_maker, supercell)
    class(mbcell_lattice_t), intent(inout):: self
    type(supercell_maker_t), intent(inout):: sc_maker
    type(mbcell_lattice_t), intent(inout):: supercell

    real(dp) :: sc_cell(3,3)
    real(dp), allocatable :: sc_xcart(:,:)
    real(dp), allocatable :: sc_masses(:)
    integer, allocatable :: sc_zion(:)

    sc_cell(:,:) = sc_maker%sc_cell(self%cell)
    call sc_maker%trans_xcart(self%cell, self%xcart, sc_xcart)
    call sc_maker%repeat(self%masses, sc_masses)
    call sc_maker%repeat(self%zion, sc_zion)
    call supercell%initialize(natom=self%natom*sc_maker%ncells, cell=sc_cell, xcart= sc_xcart, masses=sc_masses, zion=sc_zion)

    ABI_SFREE(sc_xcart)
    ABI_SFREE(sc_masses)
    ABI_SFREE(sc_zion)
  end subroutine latt_fill_supercell


  !-------------------------------------------------------------------!
  !latt_from_unitcell: build lattice supercell from primitive cell
  ! same as above, only the order of argument differs.
  !-------------------------------------------------------------------!
  subroutine latt_from_unitcell(self, sc_maker, unitcell)
    class(mbcell_lattice_t), intent(inout):: self
    type(supercell_maker_t), intent(inout):: sc_maker
    type(mbcell_lattice_t), intent(inout):: unitcell
    call unitcell%fill_supercell(sc_maker, self)
  end subroutine latt_from_unitcell

  !========================= SPIN =================================


  Subroutine spin_initialize(self, nspin)
    class(mbcell_spin_t) , intent(inout):: self
    integer, intent(in) :: nspin
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    self%nspin=nspin
    call xmpi_bcast(self%nspin, master, comm, ierr)
    ABI_ALLOCATE(self%spin_positions, (3, self%nspin))
    ABI_ALLOCATE(self%ms, (self%nspin))
    ABI_ALLOCATE(self%gyro_ratio, (self%nspin))
    ABI_ALLOCATE(self%gilbert_damping, (self%nspin))
    ABI_ALLOCATE(self%rvec,(3, self%nspin) )
    ABI_ALLOCATE(self%ispin_prim,(self%nspin) )
  end subroutine spin_initialize


  Subroutine spin_set(self, nspin, ms, spin_positions, gyro_ratio, gilbert_damping, rvec, ispin_prim)
    class(mbcell_spin_t) , intent(inout):: self
    integer, intent(in) :: nspin
    real(dp), intent(in) :: ms(nspin), spin_positions(3, nspin), gyro_ratio(nspin), gilbert_damping(nspin)
    integer, optional, intent(in) :: rvec(3, nspin), ispin_prim(nspin)
    integer :: i
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if (iam_master) then
       self%ms(:) = ms(:)
       self%spin_positions(:,:)=spin_positions(:,:)
       self%gyro_ratio(:)=gyro_ratio(:)
       self%gilbert_damping(:)=gilbert_damping(:)
       if (present(rvec)) then
          self%rvec(:,:)=rvec(:,:)
       else
          self%rvec(:,:)=0
       end if
       if (present(ispin_prim)) then
          self%ispin_prim(:)=ispin_prim
       else
          do i =1 , nspin
             self%ispin_prim(i)=i
          end do
       end if
    endif
    call xmpi_bcast(self%spin_positions, master, comm, ierr)
    call xmpi_bcast(self%ms, master, comm, ierr)
    call xmpi_bcast(self%gyro_ratio, master, comm, ierr)
    call xmpi_bcast(self%gilbert_damping, master, comm, ierr)
    call xmpi_bcast(self%rvec, master, comm, ierr)
    call xmpi_bcast(self%ispin_prim, master, comm, ierr)
  end Subroutine spin_set


  subroutine spin_finalize(self)
    class(mbcell_spin_t) :: self
    self%nspin=0
    if (allocated(self%ms)) then
       ABI_DEALLOCATE(self%ms)
    end if
    if (allocated(self%spin_positions)) then
       ABI_DEALLOCATE(self%spin_positions)
    end if
    if (allocated(self%gyro_ratio)) then
       ABI_DEALLOCATE(self%gyro_ratio)
    end if
    if (allocated(self%gilbert_damping)) then
       ABI_DEALLOCATE(self%gilbert_damping)
    end if
    if (allocated(self%ispin_prim)) then
       ABI_DEALLOCATE(self%ispin_prim)
    end if
    if (allocated(self%rvec)) then
       ABI_DEALLOCATE(self%rvec)
    end if
  end subroutine spin_finalize

  subroutine spin_fill_supercell(self, sc_maker, supercell)
    class(mbcell_spin_t),intent(inout) :: self
    type(supercell_maker_t), intent(inout):: sc_maker
    type(mbcell_spin_t), intent(inout) :: supercell
    integer :: i, nspin
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    nspin=sc_maker%ncells*self%nspin
    call supercell%initialize(nspin=nspin)

    if(iam_master) then
       call sc_maker%repeat(self%ms, supercell%ms)
       call sc_maker%repeat(self%gyro_ratio, supercell%gyro_ratio)
       call sc_maker%repeat(self%gilbert_damping, supercell%gilbert_damping)
       call sc_maker%repeat([(i ,i=1, self%nspin)], supercell%ispin_prim)
       supercell%rprimd(:,:)=sc_maker%sc_cell(self%rprimd)
       call sc_maker%trans_xcart(self%rprimd, self%spin_positions, supercell%spin_positions)
       call sc_maker%rvec_for_each(self%nspin, supercell%rvec)
    end if
    call xmpi_bcast(supercell%nspin, master, comm, ierr)
    call xmpi_bcast(supercell%ms, master, comm, ierr)
    call xmpi_bcast(supercell%gyro_ratio, master, comm, ierr)
    call xmpi_bcast(supercell%gilbert_damping, master, comm, ierr)
    call xmpi_bcast(supercell%ispin_prim, master, comm, ierr)
    call xmpi_bcast(supercell%spin_positions, master, comm, ierr)
    call xmpi_bcast(supercell%rvec, master, comm, ierr)
    call xmpi_bcast(supercell%rprimd, master, comm, ierr)
  end subroutine spin_fill_supercell

  subroutine spin_from_unitcell(self, sc_maker, unitcell)
    class(mbcell_spin_t),intent(inout) :: self
    type(supercell_maker_t), intent(inout):: sc_maker
    type(mbcell_spin_t), intent(inout) :: unitcell
    call unitcell%fill_supercell(sc_maker, self)
  end subroutine spin_from_unitcell


  !========================= LWF =================================
  Subroutine lwf_initialize(self)
    class(mbcell_lwf_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine lwf_initialize

  subroutine lwf_finalize(self)
    class(mbcell_lwf_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine lwf_finalize

  subroutine lwf_fill_supercell(self, sc_maker,supercell)
    class(mbcell_lwf_t) :: self
    type(supercell_maker_t):: sc_maker
    type(mbcell_lwf_t) :: supercell
    call supercell%from_unitcell(sc_maker, self)
  end subroutine lwf_fill_supercell

  subroutine lwf_from_unitcell(self, sc_maker, unitcell)
    class(mbcell_lwf_t) :: self
    type(supercell_maker_t):: sc_maker
    type(mbcell_lwf_t) :: unitcell
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(sc_maker)
    ABI_UNUSED_A(unitcell)
  end subroutine lwf_from_unitcell

end module m_multibinit_cell

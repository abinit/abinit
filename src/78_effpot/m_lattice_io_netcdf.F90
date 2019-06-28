!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lattice_harmonic_primitive_potential
!! NAME
!! m_lattice_harmonic_primitive_potential
!!
!! FUNCTION
!! This module contains the example lattice primitive potential. It is not the
!! one really used. The purpose of this is to show how to extend multibinit
!! to new degree of freedom
!!
!! Datatypes:
!!  lattice_harmonic_primitive_potential_t
!!
!! Subroutines:
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

module m_lattice_harmonic_primitive_potential
  use iso_c_binding
  !use m_dynamic_array, only: int_array_type, real_array_type, int2d_array_type
  !use m_mathfuncs
  use defs_basis
  use m_abicore
  use m_errors
  use m_nctk
  !#if defined HAVE_NETCDF
  use netcdf
  !#endif
  use m_xmpi
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t
  use m_primitive_potential, only: primitive_potential_t
  use m_abstract_potential, only: abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_supercell_maker, only: supercell_maker_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_lattice_harmonic_potential, only: lattice_harmonic_potential_t
  implicit none
  private
  !!*** 
  type, public, extends(primitive_potential_t) :: lattice_harmonic_primitive_potential_t
     integer :: natom
     type(ndcoo_mat_t) :: coeff  ! 
     integer, allocatable :: Rlist(:,:) ! 
   contains
     procedure:: initialize
     procedure:: finalize
     procedure :: load_from_netcdf
     procedure:: fill_supercell
  end type lattice_harmonic_primitive_potential_t

contains

  subroutine initialize(self, primcell)
    class(lattice_harmonic_primitive_potential_t), intent(inout) :: self
    type(mbcell_t), target, intent(inout) :: primcell
    !integer, intent(in) :: nspin
    self%primcell=>primcell
    self%label="lattice_harmonic_primitive_potential"
    self%has_spin=.False.
    self%has_displacement=.True.
    self%has_strain=.False.
    self%has_lwf=.False.
  end subroutine initialize


  subroutine finalize(self)
    class(lattice_harmonic_primitive_potential_t), intent(inout) :: self
    call self%coeff%finalize()
    ABI_DEALLOCATE(self%Rlist)
    nullify(self%primcell)
    self%natom=0
    self%label="Destroyed lattice_harmonic_primitive_potential"
    call self%primitive_potential_t%finalize()
  end subroutine finalize


  !-------------------------------------------------------------------!
  ! Load from files:
  !> params: parameters
  !> fnames: file names from files file
  !-------------------------------------------------------------------!
  subroutine load_from_files(self, params, fnames)
    class(lattice_harmonic_primitive_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    call self%load_from_netcdf( fnames(1))
    ABI_UNUSED_A(params)
  end subroutine load_from_files

  subroutine load_from_netcdf(self, fname)
    class(lattice_harmonic_primitive_potential_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: fname
    integer :: ncid, ierr
    integer :: iR, R(3), nR, natom, natom3
    real(dp) :: cell(3,3)
    real(dp), allocatable :: xcart(:,:), masses(:)
    integer, allocatable :: zion(:)
    real(dp), allocatable :: ifc_vallist(:,:,:)
    integer :: varid, i, j, d(3)


    ierr=nf90_open(fname, NF90_NOWRITE, ncid)


    ierr=nctk_get_dim(ncid, "ifc_nR" , nR)
    ierr=nctk_get_dim(ncid, "natom", natom)
    ierr=nctk_get_dim(ncid, "natom3", natom3)

    ABI_ALLOCATE(masses, (natom))
    ABI_ALLOCATE(xcart, (3, natom))
    ABI_ALLOCATE(zion,(nR))

    ABI_ALLOCATE(ifc_vallist, (nR, natom3, natom3))
    ABI_ALLOCATE(self%Rlist,(3, nR))

    call self%coeff%initialize([ nR, natom3, natom3 ])
    call self%primcell%set_lattice(natom, cell, xcart, masses, zion)

    ierr =nf90_inq_varid(ncid, "ref_masses", varid)
    ierr = nf90_get_var(ncid, varid, masses)

    ierr =nf90_inq_varid(ncid, "ref_xcart", varid)
    ierr = nf90_get_var(ncid, varid, xcart)

    ierr =nf90_inq_varid(ncid, "ref_cell", varid)
    ierr = nf90_get_var(ncid, varid, cell)

    ierr =nf90_inq_varid(ncid, "ifc_Rlist", varid)
    ierr = nf90_get_var(ncid, varid, self%Rlist)

    ierr =nf90_inq_varid(ncid, "ifc_varlist", varid)
    ierr = nf90_get_var(ncid, varid, ifc_vallist)

    do iR =1, nR
       do i=1 , natom3
          do j=1, natom3
             call self%coeff%add_entry([iR, i, j], ifc_vallist(iR, i, j))
          end do
       end do
    end do

    ierr=nf90_close(ncid)

    ABI_DEALLOCATE(masses)
    ABI_DEALLOCATE(xcart)
    ABI_DEALLOCATE(zion)
    ABI_DEALLOCATE(ifc_vallist)
  end subroutine load_from_netcdf


  subroutine fill_supercell(self, scmaker, scpot)
    class(lattice_harmonic_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t), intent(inout):: scmaker
    class(abstract_potential_t), pointer, intent(inout) :: scpot
    integer :: natom, sc_natom
    integer, allocatable :: ilist_sc(:), jlist_sc(:), Rlist_sc(:,:)
    real, allocatable :: vallist_sc(:)

    natom=self%natom
    sc_natom= natom* scmaker%ncells

    !! NOTE: the code below can be used as a pattern to build supercell from primitivecell.
    ! Step 1: allocate the scpot as a corresponding supercell potential
    ABI_DATATYPE_ALLOCATE_SCALAR(lattice_harmonic_potential_t, scpot)
    ! Fortran does not know the type of the generic pointer
    ! unless select type is used:
    select type(scpot)
    type is (lattice_harmonic_potential_t)
       call scpot%initialize(sc_natom)
       ! IFC is an COO_mat_t, which has the index of R1, R2, R3, i, j and the value of val
       ! list of index R: coeff%ind%data(1, 1:coeff%nnz)
       ! list of i: coeff%ind%data(2, 1:coeff%nnz)
       ! list of j: coeff%ind%data(3, 1:coeff%nnz)
       ! IFC: (ind) = val
       ! TODO
       ! trans_ilist: for each i, it will find one i' for each cell in supercell (total: ncell)
       ! nbasis: number of degrees of freedom in primitive cell, here is 3natom.
    end select
    ABI_DEALLOCATE(ilist_sc)
    ABI_DEALLOCATE(jlist_sc)
    ABI_DEALLOCATE(Rlist_sc)
    ABI_DEALLOCATE(vallist_sc)
  end subroutine fill_supercell




end module m_lattice_harmonic_primitive_potential

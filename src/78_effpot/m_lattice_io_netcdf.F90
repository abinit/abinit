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
  use defs_basis
  use m_abicore
  use m_errors
  use m_nctk
!#if defined HAVE_NETCDF
  use netcdf
!#endif
  use m_xmpi
  use m_mathfuncs, only: eigensh
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


  !-------------------------------------------------------------------!
  ! An harmonic potential which has only the IFC (with no dipole-dipole)
  !
  ! IFC is written in a coefficient matrix M(R, i, j)=val
  !-------------------------------------------------------------------!
  type, public, extends(primitive_potential_t) :: lattice_harmonic_primitive_potential_t
     integer :: natom    ! number of atoms
     type(ndcoo_mat_t) :: coeff  !  A N-dimensional COO matrix.
                                 !  The indices are (ind_R, i, j). Note that xyz is included in i and j.
     integer, allocatable :: Rlist(:,:) ! The list of R points (3, number of R-points)
                                        !. Rlist(:, ind_R) is a R-vector.
   contains
     procedure:: initialize
     procedure:: finalize
     procedure :: load_from_files   ! load potential from files listed in files file
     procedure :: load_from_netcdf  ! load potential from a netcdf file
     procedure:: fill_supercell     ! fill a supercell potential
     procedure :: get_hamk          ! generate hamiltonian for one k point.
     procedure :: get_eigen         ! eigen values and eigen vectors
  end type lattice_harmonic_primitive_potential_t

contains

  !-------------------------------------------------------------------!
  ! Initialize
  !  Input:
  !    primcell: the reference primitive cell.
  !-------------------------------------------------------------------!
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


  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
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
    call self%load_from_netcdf( fnames(3))
    ABI_UNUSED_A(params)
  end subroutine load_from_files

  !-------------------------------------------------------------------!
  ! load potential from netcdf file
  ! Note that the lattic part of the primitive cell is also loaded.
  ! Input:
  !  fname: filename
  !-------------------------------------------------------------------!
  subroutine load_from_netcdf(self, fname)
    class(lattice_harmonic_primitive_potential_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: fname
    integer :: ncid, ierr
    integer :: iR, nR, natom, natom3
    real(dp) :: cell(3,3)
    real(dp), allocatable :: xcart(:,:), masses(:)
    integer, allocatable :: zion(:)
    real(dp), allocatable :: ifc_vallist(:,:,:)
    integer :: varid, i, j
!#if defined HAVE_NETCDF
    ierr=nf90_open(trim(fname), NF90_NOWRITE, ncid)
    call nc_handle_err(ierr)

    ierr=nctk_get_dim(ncid, "ifc_nR" , nR)
    ierr=nctk_get_dim(ncid, "natom", natom)
    ierr=nctk_get_dim(ncid, "natom3", natom3)


    ABI_ALLOCATE(masses, (natom))
    ABI_ALLOCATE(xcart, (3, natom))
    ABI_ALLOCATE(zion,(natom))

    ABI_ALLOCATE(ifc_vallist, (natom3, natom3, nR))
    ABI_ALLOCATE(self%Rlist,(3, nR))


    ierr =nf90_inq_varid(ncid, "ref_masses", varid)
    call nc_handle_err(ierr, "ref_masses")
    ierr = nf90_get_var(ncid, varid, masses)
    call nc_handle_err(ierr, "ref_masses")

    ierr =nf90_inq_varid(ncid, "ref_xcart", varid)
    call nc_handle_err(ierr, "ref_xcart")
    ierr = nf90_get_var(ncid, varid, xcart)
    call nc_handle_err(ierr, "ref_xcart")

    ierr =nf90_inq_varid(ncid, "ref_cell", varid)
    call nc_handle_err(ierr, "rec_cell")
    ierr = nf90_get_var(ncid, varid, cell)
    call nc_handle_err(ierr, "ref_cell")

    ierr =nf90_inq_varid(ncid, "ref_zion", varid)
    call nc_handle_err(ierr, "ref_zion")
    ierr = nf90_get_var(ncid, varid, zion)
    call nc_handle_err(ierr, "ref_zion")

    masses(:)  = masses(:) * amu_emass
    cell(:,:) = cell(:,:) / Bohr_Ang
    xcart(:,:) = xcart(:,:) /Bohr_Ang

    self%natom=natom
    call self%primcell%set_lattice(natom, cell, xcart, masses, zion)

    call self%coeff%initialize([ nR, natom3, natom3 ])
    ierr =nf90_inq_varid(ncid, "ifc_Rlist", varid)
    call nc_handle_err(ierr, "ifc_Rlist")
    ierr = nf90_get_var(ncid, varid, self%Rlist)
    call nc_handle_err(ierr, "ifc_Rlist")

    ierr =nf90_inq_varid(ncid, "ifc_vallist", varid)
    call nc_handle_err(ierr, "ifc_vallist")
    ierr = nf90_get_var(ncid, varid, ifc_vallist)
    call nc_handle_err(ierr, "ifc_vallist")

    ifc_vallist(:,:,:) = ifc_vallist(:,:,:) * eV_Ha * (Bohr_Ang * Bohr_Ang)

    do iR =1, nR
       do i=1 , natom3
          do j=1, natom3
             if (abs(ifc_vallist(j, i, iR))>1e-9) then
                ! NOTE: in fortran the order of index in reversed when reading netcdf array.
                call self%coeff%add_entry([iR, i, j], ifc_vallist(j, i, iR))
             end if
          end do
       end do
    end do

    ierr=nf90_close(ncid)
    call nc_handle_err(ierr)

    ABI_DEALLOCATE(masses)
    ABI_DEALLOCATE(xcart)
    ABI_DEALLOCATE(zion)
    ABI_DEALLOCATE(ifc_vallist)
!# else
!   MSG_ERROR("Error: Netcdf should be installed with Multibinit to read the netcdf potential file.")
!#endif
  end subroutine load_from_netcdf


  !-------------------------------------------------------------------!
  !Fill supercell
  ! Inputs:
  !  scmaker : supercell_maker_t
  ! Output:
  !  scpot: a class pointer to an ABSTRACT potential.
  !         Not a harmonic potential because this is inherited from
  !         an abstract_primitive_potential_t, which doesn't know
  !         the type of the supercell potential.
  !
  !-------------------------------------------------------------------!
  subroutine fill_supercell(self, scmaker, scpot)
    use m_spmat_convert, only: COO_to_dense
    class(lattice_harmonic_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t), intent(inout):: scmaker
    class(abstract_potential_t), pointer, intent(inout) :: scpot
    integer :: natom, sc_natom
    integer :: inz, iR, R(3), i, j, icell
    integer, allocatable :: ilist_sc(:), jlist_sc(:), Rlist_sc(:,:)
    real(dp):: val


    natom=self%natom
    sc_natom= natom* scmaker%ncells

    !! NOTE: the code below can be used as a pattern to build supercell from primitivecell.
    ! Step 1: allocate the scpot as a corresponding supercell potential
    ABI_DATATYPE_ALLOCATE_SCALAR(lattice_harmonic_potential_t, scpot)
    ! Fortran does not know the functions specific to the derived class pointer.
    ! Only the ones inheritated from abstract class,
    ! unless select type is used:
    select type(scpot)
    type is (lattice_harmonic_potential_t)
       call scpot%initialize(sc_natom)
       ! IFC is an COO_mat_t, which has the index of R1, R2, R3, i, j and the value of val
       ! list of index R: coeff%ind%data(1, 1:coeff%nnz)
       ! list of i: coeff%ind%data(2, 1:coeff%nnz)
       ! list of j: coeff%ind%data(3, 1:coeff%nnz)
       ! IFC: (ind) = val
       do inz =1 , self%coeff%nnz
          ! For each non-zero entry in the coeff matrix
          ! get the R, i, and j, val
          iR=self%coeff%ind%data(1,inz)
          R(:) = self%Rlist(:, iR)
          i=self%coeff%ind%data(2, inz)
          j=self%coeff%ind%data(3, inz)
          val=self%coeff%val%data(inz)
          ! translate i to i in supercell.
          ! No need to allocate, it is done by trans_i . but remember to deallocate!
          ! nbasis is the number in one primitive cell.
          ! e.g. there are 3*natom possible i (3: x, y, z) in each primitive cell.
          call scmaker%trans_i(nbasis=self%natom*3, i=i, i_sc=ilist_sc )
          ! translate j, Rj to supercell.
          call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=j, Rj=R, j_sc=jlist_sc, Rj_sc=Rlist_sc)
          ! values are repeated in cells
          do icell=1, scmaker%ncells
             call scpot%add_term(ilist_sc(icell), jlist_sc(icell), val )
          end do
       end do

       ! Test the phonon energy
       !call COO_to_dense(scpot%coeff, real_sc_evecs)
       !sc_evecs(:,:) = real_sc_evecs(:,:)
       !call eigensh(sc_evals, sc_evecs)

    end select


    ABI_SFREE(ilist_sc)
    ABI_SFREE(jlist_sc)
    ABI_SFREE(Rlist_sc)

  end subroutine fill_supercell


  !-------------------------------------------------------------------!
  ! calculate hamiltonian at k point (or more precisely, q point)
  ! H_ij(k) = sum_R H_ij(R) exp( i2pi k.dot.R)
  !-------------------------------------------------------------------!
  subroutine get_hamk(self, kpoint, hamk)
    class(lattice_harmonic_primitive_potential_t), intent(in) :: self
    real(dp), intent(in) :: kpoint(3)
    complex(dp), intent(inout) :: hamk(:,:)
    integer :: inz, iR, R(3), i, j
    real(dp) :: val
    hamk(:,:) = cmplx(0.0,0.0)
    do inz =1, self%coeff%nnz
       iR=self%coeff%ind%data(1,inz)
       R(:) = self%Rlist(:, iR)
       i=self%coeff%ind%data(2, inz)
       j=self%coeff%ind%data(3, inz)
       val=self%coeff%val%data(inz)
       hamk(i, j) =hamk(i, j) + val*exp(cmplx(0.0,2.0) *pi * dot_product(kpoint, R))
    end do
  end subroutine get_hamk

  !-------------------------------------------------------------------!
  !calculate eigenvalue and eigen vector
  !-------------------------------------------------------------------!
  subroutine get_eigen(self, kpoint, evals, evecs)
    class(lattice_harmonic_primitive_potential_t), intent(in) :: self
    real(dp), intent(in) :: kpoint(3)
    real(dp), intent(inout) :: evals(:)
    complex(dp), intent(inout) :: evecs(:,:)
    call self%get_hamk(kpoint, evecs)
    call eigensh(evals, evecs)
  end subroutine get_eigen

  !-------------------------------------------------------------------!
  ! nc_handle_err:
  ! Print the error message of an netcdf from an ierr number
  ! Inputs:
  !  ierr: the error number
  !  name: a tag printed with the message so it is easier to debug.
  !-------------------------------------------------------------------!
  subroutine nc_handle_err(ierr, name)
    integer, intent ( in) ::ierr
    character(*), optional, intent(in) :: name
    if(ierr/= nf90_noerr) then
       if (present(name)) then
          write(std_out, *)  trim(nf90_strerror(ierr)), "when trying to read ", name
       else
          write(std_out, *)  trim(nf90_strerror(ierr))
       end if
       stop "Stopped"
    end if
  end subroutine nc_handle_err

end module m_lattice_harmonic_primitive_potential

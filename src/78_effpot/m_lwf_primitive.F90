!!****m* ABINIT/m_lwf_primitive_potential
!! NAME
!! m_lwf_primitive_potential
!!
!! FUNCTION
!! This module contains the example lwf primitive potential. It is not the
!! one really used. The purpose of this is to show how to extend multibinit
!! to new degree of freedom
!!
!! Datatypes:
!!  LWF_primitive_potential_t
!!
!! Subroutines:
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

#define HAVE_NETCDF 1

module m_lwf_primitive_potential
  use iso_c_binding
  use defs_basis
  use m_abicore
  use m_errors
  use m_nctk
#if defined HAVE_NETCDF
  use netcdf
#endif
  use m_xmpi
  use m_mathfuncs, only: eigensh
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t
  use m_primitive_potential, only: primitive_potential_t
  use m_abstract_potential, only: abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_supercell_maker, only: supercell_maker_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_lwf_potential, only: lwf_potential_t
  implicit none
  private
  !!***

  !-------------------------------------------------------------------!
  ! An harmonic potential which has only the IFC (with no dipole-dipole)
  !
  ! IFC is written in a coefficient matrix M(R, i, j)=val
  !-------------------------------------------------------------------!
  type, public, extends(primitive_potential_t) :: lwf_primitive_potential_t
     integer :: nlwf=0, natom=0, nR=0 ! number of LWF
     type(ndcoo_mat_t) :: coeff  !  A N-dimensional COO matrix.
     !  The indices are (ind_R, i, j). Note that xyz is included in i and j.
     real(dp), allocatable :: lwf_masses(:)

     integer :: onebody_nterm=0

     real(dp), allocatable :: lattice_coeffs(:,:, :) ! (natom*3, nlwf, nR)

     integer, allocatable :: onebody_i(:), onebody_order(:)
     real(dp), allocatable :: onebody_val(:)

     type(ndcoo_mat_t) :: coeff_twobody ! Two body interaction parameters of higher order.
     ! The indices are (indR, i, j, orderi, orderj), H(R, i, j)= val A_i^(orderi) A_j^(orderj)
     integer, allocatable :: Rlist(:,:) ! The list of R points (3, number of R-points)

     ! For testing the bounding.
     logical :: has_self_bound_term= .False.
     integer :: self_bound_order=0
     real(dp) :: self_bound_coeff=0.0_dp

     !. Rlist(:, ind_R) is a R-vector.
     real(dp) :: ref_energy=0.0                  ! reference energy

     logical :: as_lattice_anharmonic=.False.

   contains
     procedure:: initialize
     procedure:: finalize
     procedure :: load_from_files   ! load potential from files listed in files file
     procedure :: load_from_netcdf  ! load potential from a netcdf file
     procedure:: fill_supercell     ! fill a supercell potential
     procedure :: get_hamk          ! generate hamiltonian for one k point.
     procedure :: get_eigen         ! eigen values and eigen vectors
     procedure :: add_self_bound_term ! add a self bound term
  end type lwf_primitive_potential_t


contains

  !-------------------------------------------------------------------!
  ! Initialize
  !  Input:
  !    nlwf: number of Wannier functions
  !    primcell: the reference primitive cell.
  !-------------------------------------------------------------------!
  subroutine initialize(self, primcell)
    class(lwf_primitive_potential_t), intent(inout) :: self
    type(mbcell_t), target, intent(inout) :: primcell
    !integer, intent(in) :: nlwf
    self%primcell=>primcell
    self%label="lwf_primitive_potential"
    self%has_spin=.False.
    self%has_displacement=.False.
    self%has_strain=.False.
    self%has_lwf=.True.
  end subroutine initialize


  !-------------------------------------------------------------------!
  ! Finalize
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(lwf_primitive_potential_t), intent(inout) :: self
    call self%coeff%finalize()
    call self%coeff_twobody%finalize()
    ABI_SFREE(self%Rlist)
    nullify(self%primcell)
    self%nlwf=0
    self%label="Destroyed lwf_primitive_potential"
    call self%primitive_potential_t%finalize()

    ABI_SFREE(self%lattice_coeffs)
    ABI_SFREE(self%lwf_masses)
    self%natom=0
    self%nR=0

    if (self%onebody_nterm /= 0) then
       ABI_SFREE(self%onebody_i)
       ABI_SFREE(self%onebody_val)
       ABI_SFREE(self%onebody_order)
    end if
    self%onebody_nterm=0
  end subroutine finalize


  !-------------------------------------------------------------------!
  ! Load from files:
  !> params: parameters
  !> fnames: file names from files file
  !-------------------------------------------------------------------!
  subroutine load_from_files(self, params, fnames)
    class(lwf_primitive_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    call self%load_from_netcdf(fnames(1))
    call self%add_self_bound_term(params%lwf_self_bound_order, &
         & params%lwf_self_bound_coeff)
    self%as_lattice_anharmonic= (params%latt_lwf_anharmonic==1)
  end subroutine load_from_files

  !-------------------------------------------------------------------!
  ! load potential from netcdf file
  ! Note that the lattic part of the primitive cell is also loaded.
  ! Input:
  !  fname: filename
  !-------------------------------------------------------------------!
  subroutine load_from_netcdf(self, fname)
    class(lwf_primitive_potential_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: fname
    integer :: ncid, ierr
    integer :: iR,  nR, nlwf, natom, twobody_nterm, onebody_nterm
    real(dp) :: cell(3,3)
    real(dp), allocatable :: ifc(:, :, :), xcart(:,:), masses(:), zion(:)
    real(dp), allocatable ::  twobody_val(:)
    integer, allocatable ::   twobody_iR(:), twobody_i(:), twobody_j(:), twobody_orderi(:), twobody_orderj(:)
    integer :: varid, i, j
    !#if defined HAVE_NETCDF
    ierr=nf90_open(trim(fname), NF90_NOWRITE, ncid)
    NCF_CHECK_MSG(ierr, "Open netcdf file")

    ! read primcell info
    ierr=nctk_get_dim(ncid, "wann_natom", natom)
    NCF_CHECK_MSG(ierr, "getting natom in lwf potential file")

    ierr=nctk_get_dim(ncid, "nR" , nR)
    NCF_CHECK_MSG(ierr, "getting nR in lwf potential file")
    self%nR=nR

    ierr=nctk_get_dim(ncid, "wann_nwann", nlwf)
    NCF_CHECK_MSG(ierr, "getting wann_nwann in lwf potential file")

    call self%primcell%set_lwf(nlwf)
    self%nlwf=nlwf
    call self%coeff%initialize([ nR, nlwf, nlwf])
    call self%coeff_twobody%initialize([nR, nlwf, nlwf, -1, -1])

    ABI_ALLOCATE(ifc, (nlwf, nlwf, nR))
    ABI_ALLOCATE(self%Rlist,(3, nR))
    ABI_ALLOCATE(xcart, (3,natom))
    ABI_ALLOCATE(masses, (natom))
    ABI_ALLOCATE(zion, (natom))

    ABI_ALLOCATE(self%lwf_masses, (nlwf))

    ierr =nf90_inq_varid(ncid, "wann_lwf_masses", varid)
    NCF_CHECK_MSG(ierr, "lwf_masses")
    ierr = nf90_get_var(ncid, varid, self%lwf_masses)
    NCF_CHECK_MSG(ierr, "lwf_masses")
    self%primcell%lwf%lwf_masses=self%lwf_masses

    ierr=nf90_inq_dimid(ncid, "wann_onebody_nterm", onebody_nterm)
    if (ierr/=nf90_noerr) then
       onebody_nterm=0
    else
       ierr=nctk_get_dim(ncid, "wann_onebody_nterm", onebody_nterm)
       NCF_CHECK_MSG(ierr, "getting wann_onebody_nterm in lwf potential file")
    end if
    self%onebody_nterm=onebody_nterm


    ierr=nf90_inq_dimid(ncid, "wann_twobody_nterm", twobody_nterm)
    if (ierr/=nf90_noerr) then
       twobody_nterm=0
    else
       ierr=nctk_get_dim(ncid, "wann_twobody_nterm", twobody_nterm)
       NCF_CHECK_MSG(ierr, "getting wann_twobody_nterm in lwf potential file")
    end if

    self%natom=natom

    ABI_UNUSED(cell)

    !ierr =nf90_inq_varid(ncid, "ref_cell", varid)
    !NCF_CHECK_MSG(ierr, "ref_cell")
    !ierr = nf90_get_var(ncid, varid, cell)
    !NCF_CHECK_MSG(ierr, "ref_cell")
    !cell(:,:)=cell(:,:)/ Bohr_Ang


    !ierr =nf90_inq_varid(ncid, "ref_xcart", varid)
    !NCF_CHECK_MSG(ierr, "ref_xcart")
    !ierr = nf90_get_var(ncid, varid, xcart)
    !NCF_CHECK_MSG(ierr, "ref_xcart")

    !xcart(:,:)=xcart(:,:)/ Bohr_Ang

    ABI_MALLOC(self%lattice_coeffs, (nlwf, natom*3, nR))

    ierr =nf90_inq_varid(ncid, "wann_wannier_function_real", varid)
    NCF_CHECK_MSG(ierr, "wann_wannier_function_real")
    ierr = nf90_get_var(ncid, varid, self%lattice_coeffs)
    NCF_CHECK_MSG(ierr, "wann_wannier_function_real")

    ierr =nf90_inq_varid(ncid, "wann_Rlist", varid)
    NCF_CHECK_MSG(ierr, "wann_Rlist")
    ierr = nf90_get_var(ncid, varid, self%Rlist)
    NCF_CHECK_MSG(ierr, "wann_Rlist")

    ierr =nf90_inq_varid(ncid, "wann_ifc_real", varid)
    NCF_CHECK_MSG(ierr, "wann_ifc_real")
    ierr = nf90_get_var(ncid, varid, ifc)
    NCF_CHECK_MSG(ierr, "wann_ifc_real")


    ifc(:,:,:) = ifc(:,:,:) * eV_Ha * (Bohr_Ang * Bohr_Ang)

    do iR =1, nR
       do i=1 , nlwf
          do j=1, nlwf
             if (abs(ifc(j, i, iR))>1e-4) then
                ! NOTE: in fortran the order of index in reversed when reading netcdf array.
                call self%coeff%add_entry([iR, i, j], ifc(j, i, iR))
             end if
          end do
       end do
    end do


    if(self%onebody_nterm /= 0) then
       ABI_MALLOC(self%onebody_i, (self%onebody_nterm))
       ABI_MALLOC(self%onebody_order, (self%onebody_nterm))
       ABI_MALLOC(self%onebody_val, (self%onebody_nterm))

       ierr =nf90_inq_varid(ncid, "wann_onebody_i", varid)
       NCF_CHECK_MSG(ierr, "wann_onebody_i")
       ierr = nf90_get_var(ncid, varid, self%onebody_i)
       NCF_CHECK_MSG(ierr, "wann_onebody_i")

       self%onebody_i=self%onebody_i +1

       ierr =nf90_inq_varid(ncid, "wann_onebody_val", varid)
       NCF_CHECK_MSG(ierr, "wann_onebody_val")
       ierr = nf90_get_var(ncid, varid, self%onebody_val)
       NCF_CHECK_MSG(ierr, "wann_onebody_val")

       ierr =nf90_inq_varid(ncid, "wann_onebody_order", varid)
       NCF_CHECK_MSG(ierr, "wann_onebody_order")
       ierr = nf90_get_var(ncid, varid, self%onebody_order)
       NCF_CHECK_MSG(ierr, "wann_onebody_order")

       do i =1, self%onebody_nterm
          self%onebody_val(i)=self%onebody_val(i) * eV_Ha * (Bohr_Ang ** self%onebody_order(i))
       end do
    end if

    if(twobody_nterm /= 0) then
       ABI_MALLOC(twobody_iR, (twobody_nterm))
       ABI_MALLOC(twobody_i, (twobody_nterm))
       ABI_MALLOC(twobody_j, (twobody_nterm))
       ABI_MALLOC(twobody_orderi, (twobody_nterm))
       ABI_MALLOC(twobody_orderj, (twobody_nterm))
       ABI_MALLOC(twobody_val, (twobody_nterm))

       ierr =nf90_inq_varid(ncid, "wann_twobody_iR", varid)
       NCF_CHECK_MSG(ierr, "wann_twobody_iR")
       ierr = nf90_get_var(ncid, varid, twobody_iR)
       NCF_CHECK_MSG(ierr, "wann_twobody_iR")

       ierr =nf90_inq_varid(ncid, "wann_twobody_i", varid)
       NCF_CHECK_MSG(ierr, "wann_twobody_i")
       ierr = nf90_get_var(ncid, varid, twobody_i)
       NCF_CHECK_MSG(ierr, "wann_twobody_i")

       ierr =nf90_inq_varid(ncid, "wann_twobody_j", varid)
       NCF_CHECK_MSG(ierr, "wann_twobody_j")
       ierr = nf90_get_var(ncid, varid, twobody_j)
       NCF_CHECK_MSG(ierr, "wann_twobody_j")

       ierr =nf90_inq_varid(ncid, "wann_twobody_orderi", varid)
       NCF_CHECK_MSG(ierr, "wann_twobody_orderi")
       ierr = nf90_get_var(ncid, varid, twobody_orderi)
       NCF_CHECK_MSG(ierr, "wann_twobody_i")

       ierr =nf90_inq_varid(ncid, "wann_twobody_orderj", varid)
       NCF_CHECK_MSG(ierr, "wann_twobody_orderj")
       ierr = nf90_get_var(ncid, varid, twobody_orderj)
       NCF_CHECK_MSG(ierr, "wann_twobody_orderj")

       ierr =nf90_inq_varid(ncid, "wann_twobody_val", varid)
       NCF_CHECK_MSG(ierr, "wann_twobody_val")
       ierr = nf90_get_var(ncid, varid, twobody_val)
       NCF_CHECK_MSG(ierr, "wann_twobody_val")

       ierr=nf90_close(ncid)
       NCF_CHECK_MSG(ierr, "Close netcdf file")

       do i=1, twobody_nterm
          call self%coeff_twobody%add_entry([twobody_iR(i), twobody_i(i), &
               & twobody_j(i), twobody_orderi(i), twobody_orderj(i)], twobody_val(i))
          !TODO change to a.u.
       end do
    end if

    ABI_SFREE(masses)
    ABI_SFREE(xcart)
    ABI_SFREE(zion)
    ABI_SFREE(ifc)

    if(twobody_nterm /=0) then
       ABI_SFREE(twobody_i)
       ABI_SFREE(twobody_j)
       ABI_SFREE(twobody_orderi)
       ABI_SFREE(twobody_orderj)
       ABI_SFREE(twobody_val)
    endif
    !#else
    !NETCDF_NOTENABLED_ERROR()
    !#endif

  end subroutine load_from_netcdf

  !-------------------------------------------------------------------!
  !Fill supercell
  ! Inputs:
  !  scmaker : supercell_maker_t
  ! Output:
  !  scpot: a class pointer to an ABSTRACT potential.
  !         Not a  potential because this is inherited from
  !         an abstract_primitive_potential_t, which doesn't know
  !         the type of the supercell potential.
  !
  !-------------------------------------------------------------------!
  subroutine fill_supercell(self, scmaker, params, scpot)
    use m_spmat_convert, only: COO_to_dense

    class(lwf_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t),                        intent(inout) :: scmaker
    type(multibinit_dtset_type),                    intent(inout) :: params
    class(abstract_potential_t), pointer,           intent(inout) :: scpot

    integer :: nlwf, sc_nlwf
    integer :: inz, iR, R(3), i, j, icell
    integer, allocatable :: ilist_sc(:), jlist_sc(:), Rlist_sc(:,:)
    real(dp):: val

    ABI_UNUSED_A(params)

    nlwf=self%nlwf
    sc_nlwf= nlwf* scmaker%ncells

    ! Step 1: allocate the scpot as a corresponding supercell potential
    ABI_DATATYPE_ALLOCATE_SCALAR(lwf_potential_t, scpot)
    ! Fortran does not know the functions specific to the derived class pointer.
    ! Only the ones inheritated from abstract class,
    ! unless select type is used:
    select type(scpot)
    type is (lwf_potential_t)
       call scpot%initialize(sc_nlwf)
       ! IFC is an COO_mat_t, which has the index of R1, R2, R3, i, j and the value of val
       ! list of index R: coeff%ind%data(1, 1:coeff%nnz)
       ! list of i: coeff%ind%data(2, 1:coeff%nnz)
       ! list of j: coeff%ind%data(3, 1:coeff%nnz)
       ! IFC: (ind) = val
       if(self%as_lattice_anharmonic) then
          call scpot%use_as_lattice_anharmonic()
       else
          ! harmonic terms
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
             call scmaker%trans_i(nbasis=self%nlwf, i=i, i_sc=ilist_sc )
             ! translate j, Rj to supercell.
             call scmaker%trans_j_and_Rj(nbasis=self%nlwf, j=j, Rj=R, j_sc=jlist_sc, Rj_sc=Rlist_sc)
             ! values are repeated in cells
             do icell=1, scmaker%ncells
                call scpot%add_term(ilist_sc(icell), jlist_sc(icell), val )
             end do
             ABI_SFREE(ilist_sc)
             ABI_SFREE(jlist_sc)
             ABI_SFREE(Rlist_sc)
          end do
          end if

       ! coefficients of atomic displacements in supercell
       do i=1, self%nlwf
          call scmaker%trans_i(nbasis=nlwf, i=i, i_sc=ilist_sc)
          do icell =1, scmaker%ncells
              call scpot%lwf_latt_coeffs(ilist_sc(icell))%initialize(scmaker%ncells*self%natom*3)
          end do
          do iR=1, self%nR
             R(:) = self%Rlist(:, iR)
             do j =1, self%natom*3
                val=self%lattice_coeffs(i, j, iR) ! ilwf, iatom3, iR
                if (abs(val) > 1e-4) then
                   call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=j, Rj=R, j_sc=jlist_sc, Rj_sc=Rlist_sc)
                   do icell =1, scmaker%ncells
                      call scpot%lwf_latt_coeffs(ilist_sc(icell))%push(jlist_sc(icell), val)
                   end do
                   ABI_SFREE(jlist_sc)
                   ABI_SFREE(Rlist_sc)
                end if
             end do
          end do
          ABI_SFREE(ilist_sc)
       end do

       ! anharmonic terms
       !call scpot%set_ref_energy(self%ref_energy * scmaker%ncells)
       if (self%has_self_bound_term) then
          call scpot%add_self_bound_term(self%self_bound_order, self%self_bound_coeff)
       end if

       if (self%onebody_nterm/=0) then
          do i=1, self%onebody_nterm
             call scmaker%trans_i(nbasis=self%nlwf, i=self%onebody_i(i), i_sc=ilist_sc)
             do icell=1, scmaker%ncells
                call scpot%add_onebody_term(i=ilist_sc(icell), order=self%onebody_order(i), val=self%onebody_val(i))
             end do
             ABI_SFREE(ilist_sc)
          end do
       endif
       ! Test the phonon energy
       !call COO_to_dense(scpot%coeff, real_sc_evecs)
       !sc_evecs(:,:) = real_sc_evecs(:,:)
       !call eigensh(sc_evals, sc_evecs)

    end select


  end subroutine fill_supercell

  !-------------------------------------------------------------------!
  ! calculate hamiltonian at k point (or more precisely, q point)
  ! H_ij(k) = sum_R H_ij(R) exp( i2pi k.dot.R)
  !-------------------------------------------------------------------!
  subroutine get_hamk(self, kpoint, hamk)
    class(lwf_primitive_potential_t), intent(in) :: self
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
    class(lwf_primitive_potential_t), intent(in) :: self
    real(dp), intent(in) :: kpoint(3)
    real(dp), intent(inout) :: evals(:)
    complex(dp), intent(inout) :: evecs(:,:)
    call self%get_hamk(kpoint, evecs)
    ! The evecs array is reused both as the matrix and eigenvectors.
    call eigensh(evals, evecs)
  end subroutine get_eigen

  subroutine add_self_bound_term(self, order, coeff)
    class(lwf_primitive_potential_t), intent(inout) :: self 
    integer, intent(in) :: order
    real(dp), intent(in) :: coeff
    if (order /= 0) then
       self%has_self_bound_term=.True.
       self%self_bound_order = order
       self%self_bound_coeff= coeff
    end if

  end subroutine add_self_bound_term


end module m_lwf_primitive_potential

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spmat_ndcoo
!! NAME
!! m_spmat_ndcoo
!!
!! FUNCTION
!! This module contains the a NDCOO (n-dimensional coordinate) format of sparse matrix.
!! Datatypes:
!!  NDCOO_mat_t: ND COO matrix
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
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

module m_spmat_NDCOO
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_spmat_base
  use m_dynamic_array, only: int2d_array_type, real_array_type, int_array_type
  implicit none
  !!***
  private

  type, public :: ndcoo_mat_t
     integer :: ndim=0
     integer :: nnz=0
     integer, allocatable :: mshape(:)
     type(int2d_array_type) :: ind
     type(real_array_type) :: val
     logical :: is_sorted = .False.
     logical :: is_unique = .False.
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: add_entry
     procedure :: remove_zeros
     procedure :: sort_indices
     procedure :: sum_duplicates
     procedure :: get_val_inz
     procedure :: get_ind_inz
     procedure :: get_ind
     procedure :: group_by_1dim
     procedure :: vec_product
     procedure :: mv1vec
     !procedure :: print
  end type ndcoo_mat_t

  public:: test_ndcoo
contains

  !-------------------------------------------------------------------!
  ! ndcoo_mat_t initializer:
  ! Input:
  !  mshape: the shape of the N-dimension matrix. array(ndim)
  !-------------------------------------------------------------------!
  subroutine initialize(self, mshape)
    class(ndcoo_mat_t), intent(inout) :: self
    integer, intent(in) :: mshape(:)
    self%ndim=size(mshape)
    ABI_ALLOCATE(self%mshape, (self%ndim))
    self%mshape=mshape
    self%nnz=0
    self%is_sorted=.False.
    self%is_unique=.False.
  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Finalizer of ndcoo_mat_t
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(ndcoo_mat_t), intent(inout) :: self
    self%ndim=0
    self%nnz=0
    self%is_sorted=.False.
    self%is_unique=.False.
    if (allocated(self%mshape)) then
       ABI_DEALLOCATE(self%mshape)
    endif
    call self%ind%finalize()
    call self%val%finalize()
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Add one entry to the ndcoo_mat_t
  ! Inputs:
  !   ind: indices of the matrix.
  !   val: value of matrix.
  ! Example:
  !  call m%add_entry([1,2,3], 0.5)
  !-------------------------------------------------------------------!
  subroutine add_entry(self, ind, val)
    class(ndcoo_mat_t), intent(inout) :: self
    integer, intent(in) :: ind(self%ndim)
    real(dp), intent(in) :: val
    self%nnz=self%nnz+1
    call self%ind%push(ind)
    call self%val%push(val)
    self%is_sorted=.False.
    self%is_unique=.False.
  end subroutine add_entry


  !-------------------------------------------------------------------!
  ! sort the entries by indices. (left to right)
  !-------------------------------------------------------------------!
  subroutine sort_indices(self)
    class(ndcoo_mat_t), intent(inout) :: self
    real(dp) :: tmp(self%nnz)
    integer :: reorder(self%nnz)
    if(self%is_sorted .or. self%nnz==0) return
    call self%ind%sort(order=reorder)
    tmp(:)=self%val%data(1:self%nnz)
    self%val%data(1:self%nnz)=tmp(reorder)
    self%is_sorted=.True.
  end subroutine sort_indices


  !-------------------------------------------------------------------!
  ! Remove zero entries in coo matrix.
  !-------------------------------------------------------------------!
  subroutine remove_zeros(self, eps)
    class(ndcoo_mat_t), intent(inout) :: self
    real(dp), optional, intent(in) :: eps
    real(dp) :: eps1
    integer :: i, counter
    if (present(eps)) then
       eps1=eps
    else
       eps1=epsilon(1.0_dp)
    end if
    counter=0
    do i=1, self%nnz
       if (abs(self%val%data(i))> epsilon(1.0)) then
          counter=counter+1
          self%ind%data(:,counter) =self%ind%data(:, i)
          self%val%data(counter) = self%val%data(i)
       end if
    end do
    self%nnz=counter
    self%ind%size=counter
    self%val%size=counter
  end subroutine remove_zeros

  !-------------------------------------------------------------------!
  ! sum duplicate entries (also sort by indices)
  !-------------------------------------------------------------------!
  subroutine sum_duplicates(self)
    class(ndcoo_mat_t), intent(inout) :: self
    integer :: new_ind(self%ndim, self%nnz), i, counter
    real(dp) :: new_val(self%nnz)
    if (self%nnz==0) then
       self%is_unique=.True.
       return
    end if
    call self%remove_zeros()
    if (.not. self%is_sorted) then
       call self%sort_indices()
    end if
    counter=1
    new_ind(:, counter)= self%ind%data(:, 1)
    new_val(counter)=self%val%data(1)
    do i=2, self%nnz
       if (all(self%ind%data(:, i)==self%ind%data(:, i-1))) then
          new_val(counter)=new_val(counter)+self%val%data(i)
       else
          counter=counter+1
          new_ind(:, counter)= self%ind%data(:, i)
          new_val(counter)=self%val%data(i)
       end if
    end do
    self%nnz=counter
    self%ind%data(:,1:counter)=new_ind(:,1:counter)
    self%val%data(1:counter)=new_val(1:counter)
    self%ind%size=self%nnz
    self%val%size=self%nnz
    self%is_unique=.True.
  end subroutine sum_duplicates

  !-------------------------------------------------------------------!
  ! Get the i'th value of the matrix.
  !-------------------------------------------------------------------!
  function get_val_inz(self, i) result(v)
    class(ndcoo_mat_t), intent(inout) :: self
    integer, intent(in) :: i
    real(dp) :: v
    v= self%val%data(i)
  end function get_val_inz


  !-------------------------------------------------------------------!
  ! get all the indices for the ith entry
  ! Input:
  !  i: ith entry
  ! Return:
  !  a integer array of indices.
  !-------------------------------------------------------------------!
  function get_ind_inz(self, i) result(ind)
    class(ndcoo_mat_t), intent(inout) :: self
    integer, intent(in) :: i
    integer :: ind(self%ndim)
    ind(:)=self%ind%data(:,i)
  end function get_ind_inz

  !-------------------------------------------------------------------!
  ! Group the sparse matrix by the first dimension
  !> Output:
  !> ngroup: number of groups
  !> i1_list: list of 1st indices (array(ngroup))
  !> istartend: start and end of each group (array(ngroup+1))
  !>           The starts will be istartend(1:ngroup)
  !>           The ends will be istartend(2: ngroup+1)-1
  !-------------------------------------------------------------------!
  subroutine group_by_1dim(self, ngroup, i1_list, istartend)
    class(ndcoo_mat_t), intent(inout) :: self
    integer, intent(inout) :: ngroup
    integer, allocatable, intent(inout) :: i1_list(:), istartend(:)
    integer :: i, ii
    type(int_array_type) :: j1, jstartend
    if (.not. (self%is_unique))  then
       call self%sum_duplicates()
    end if
    if (self%nnz<1) then
       ngroup=0
    else if (self%nnz==1) then
       i=1
       ii=self%ind%data(1,i)
       call j1%push(ii)
       call jstartend%push(1)
       call jstartend%push(2)
    else
       i=1
       ii=self%ind%data(1,i)
       call j1%push(ii)
       call jstartend%push(i)
       do i=2, self%nnz
          ii=self%ind%data(1,i)
          if(ii == self%ind%data(1, i-1)) then
             cycle
          else
             call j1%push(ii)
             call jstartend%push(i)
          end if
       end do
       call jstartend%push(self%nnz+1)
    end if
    ngroup=j1%size
    if(ngroup>0) then
       ABI_ALLOCATE(i1_list, (ngroup))
       ABI_ALLOCATE(istartend, (ngroup+1))
       i1_list(:)=j1%data(1: j1%size)
       istartend(:)=jstartend%data(1: jstartend%size)
    end if
    call j1%finalize()
    call jstartend%finalize()
  end subroutine group_by_1dim

  !-------------------------------------------------------------------!
  ! Get the indices of the dim'th dimension
  ! Input:
  !   dim: dimension
  ! Returns:
  !   a integer array(nnz)
  !-------------------------------------------------------------------!
  function get_ind(self, dim) result(ilist)
    class(ndcoo_mat_t), intent(inout) :: self
    integer, intent(in) :: dim
    integer :: ilist(self%nnz)
    ilist(:)=self%ind%data(dim, 1:self%nnz)
  end function get_ind


  ! matrix vector product. 
  subroutine mv1vec(self, vec, i, res)
    class(ndcoo_mat_t), intent(inout) :: self
    real(dp), intent(in) :: vec(:)
    integer ,intent(in) :: i               !
    class(ndcoo_mat_t), intent(inout) :: res ! result
    !integer :: iind
    !TODO: to be implemented
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(vec)
    ABI_UNUSED_A(i)
    ABI_UNUSED_A(res)
  end subroutine mv1vec

  ! matrix vector vector  product. matrix should be dim3.
  ! n(vector)=ndim-1
  ! which returns a vecor
  ! res_r = \sum_ij M_{ijr} V_i V_j
  ! i, j, r can be in any order.
  subroutine vec_product(self, iv, veci, jv, vecj, rv, res)
    class(ndcoo_mat_t), intent(inout) :: self
    real(dp), intent(in) :: veci(:), vecj(:)
    integer ,intent(in) :: iv, jv, rv               !
    real(dp), intent(inout) :: res(:)
    integer :: iind, iiv, ijv, irv
    do iind =1 , self%nnz
      iiv=self%ind%data(iv, iind)
      ijv=self%ind%data(jv, iind)
      irv=self%ind%data(rv, iind)
      res(irv) = res(irv) + self%val%data(iind) * veci(iiv)*vecj(ijv)
    end do
  end subroutine vec_product


  subroutine test_ndcoo()
    type(ndcoo_mat_t) :: m
    integer :: ngroup
    integer, allocatable :: i1list(:), ise(:)
    call m%initialize(mshape=[3,3,3])
    call m%add_entry(ind=[3, 2,1], val=0.3d0)
    call m%add_entry(ind=[1, 2,1], val=0.3d0)
    call m%add_entry(ind=[1, 2,1], val=0.4d0)
    call m%add_entry(ind=[3, 2,1], val=0.5d0)
    call m%add_entry(ind=[1, 1,2], val=0.5d0)
    call m%add_entry(ind=[2,5,1], val=0.0d0)
    !call m%print()
    call m%sort_indices()
    call m%sum_duplicates()
    !print *, "After sum"
    !call m%print()
    !print *, "Grouping"
    call m%group_by_1dim(ngroup, i1list, ise)
    !print *,  "ngroup: ", ngroup
    !print *, "i1list: ", i1list
    !print *, "ise: ", ise
    if(allocated(i1list)) ABI_DEALLOCATE(i1list)
    if(allocated(ise)) ABI_DEALLOCATE(ise)
  end subroutine test_ndcoo

end module m_spmat_NDCOO


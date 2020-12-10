!!****m* ABINIT/m_spmat_lco
!! NAME
!! m_spmat_lco
!!
!! FUNCTION
!! This module contains the a LCO (lcordinate) format of sparse matrix.
!! The efficiency of mat vec multiplication is fine but not as good as CSR
!! Datatypes:
!!  LCO_mat_t: LCO matrix
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

module m_spmat_lco
  use defs_basis  
  use m_xmpi
  use m_errors
  use m_abicore
  use m_dynamic_array, only: int_array_type, real_array_type
  use m_spmat_base
  implicit none
  private
!!***

  !!----------- LCO ------------------------
  ! LCO sparse matrix.
  ! i, j, val are the row index, col index and value of each entry.
  ! nnz: number of non-zeros.
  type, public, extends(base_mat2d_t) :: LCO_mat_t
     integer :: nnz
     type(int_array_type),allocatable :: icol(:)
     type(real_array_type), allocatable :: val(:)
     logical :: is_unique=.False. , is_sorted=.False.
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: add_entry
     procedure :: sort_indices
     procedure :: sum_duplicates

  end type LCO_mat_t


contains

  subroutine initialize(self, mshape)
    class(lco_mat_t), intent(inout) :: self
    integer, intent(in) :: mshape(:)
    if (size(mshape)/=2) ABI_ERROR("LCO_mat should be of dimension 2.")
    call self%base_mat2d_t%initialize(mshape)
    ABI_MALLOC(self%icol, (self%nrow))
    ABI_MALLOC(self%val, (self%nrow))
    self%is_sorted=.False.
    self%is_unique=.False.


  end subroutine initialize

  subroutine finalize(self)
    class(lco_mat_t), intent(inout) :: self
    integer :: i
    do i=1, self%nrow
       call self%icol(i)%finalize()
       call self%val(i)%finalize()
    end do
    if(allocated(self%icol)) ABI_FREE(self%icol)
    if(allocated(self%val)) ABI_FREE(self%val)
    self%nnz=0
    call self%base_mat2d_t%finalize()
  end subroutine finalize

  subroutine add_entry(self, ind, val)
    class(LCO_mat_t), intent(inout) :: self
    integer, intent(in) :: ind(self%ndim)
    real(dp), intent(in) :: val
    integer :: irow, icol
    irow=ind(1)
    icol=ind(2)
    self%nnz=self%nnz+1
    call self%icol(irow)%push(icol)
    call self%val(irow)%push(val)
    self%is_sorted=.False.
    self%is_unique=.False.
  end subroutine add_entry

  subroutine sum_duplicates(self, eps)
    class(LCO_mat_t), intent(inout) :: self
    real(dp), optional, intent(in) :: eps
    real(dp) :: eps1
    integer :: i, counter, irow,nnz
    if (present(eps)) then
       eps1=eps
    else
       eps1=epsilon(1.0)
    end if
    nnz=0
    do irow=1, self%nrow
       associate(indcol=>self%icol(irow)%data)
         associate(val=>self%val(irow)%data)
           counter=0
           do i=2, self%icol(irow)%size
              if (abs(val(i))> epsilon(1.0)) then
                 if (indcol(i)==indcol(i-1)) then
                    val(counter) = val(counter)+ val(i)
                 else
                    counter=counter+1
                    indcol(counter) =indcol(i)
                    val(counter) = val(i)
                 end if
              end if
           end do
           self%icol(irow)%size=counter
           self%val(irow)%size=counter
           nnz=nnz+counter
         end associate
       end associate
    end do
    self%nnz=nnz
  end subroutine sum_duplicates


  subroutine sort_indices(self)
    class(LCO_mat_t), intent(inout) :: self
    integer ::  irow
    integer, allocatable :: order(:)
    do irow=1, self%nrow
       associate(indcol=>self%icol(irow))
         associate(val=>self%val(irow))
           ABI_MALLOC(order, (indcol%size))
           call indcol%sort(order=order)
           val%data(1:val%size) =val%data(order)
           ABI_FREE(order)
         end associate
       end associate
    end do
  end subroutine sort_indices


  ! LCO sparse matrix-vector multiplication. naive implementation.
  subroutine mv(self, x, b)
    class(LCO_mat_t), intent(in) :: self
    real(dp), intent(in) :: x(self%mshape(1))
    real(dp), intent(out) :: b(self%mshape(2))
    integer:: i, ind_i, ind_j
    b(:)=0.0D0
    do ind_i = 1, self%nrow
       do i=1, self%icol(ind_i)%size
          ind_j=self%icol(ind_i)%data(i)
          b(ind_i)=b(ind_i)+self%val(ind_i)%data(i)*x(ind_j)
       end do
    end do
  end subroutine  mv



end module m_spmat_LCO

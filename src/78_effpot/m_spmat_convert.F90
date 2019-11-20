!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spmat_convert
!! NAME
!! m_spmat_convert
!!
!! FUNCTION
!! This module contains the functions to convert between sparse matrix formats
!!
!!
!! Subroutines:
!! - dense_to_LIL
!! - LIL_to_dense
!! - LIL_to_COO
!! - LIL_to_CSR
!! - Dense_to_CSR
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


module m_spmat_convert
  use defs_basis
  use m_xmpi
  use m_abicore
  use m_errors
  use m_spmat_dense
  use m_linked_list
  use m_spmat_lil
  use m_spmat_csr
  use m_spmat_coo
  use m_spmat_lco
  implicit none
!!****m*
  public
  !-----------------------------------------------------------------------
  !> @brief convert one type of matrix Amat to the other Bmat
  !> Note that the dense matrix is not the dense_mat_t because a
  !> simple 2d array is more generally used.
  !> @param [in] Amat
  !> @param [out] Bmat
  !-----------------------------------------------------------------------
  interface spmat_convert
     !procedure  dense_to_LIL
     procedure  LIL_to_dense
     procedure  LIL_to_COO
     procedure  LIL_to_CSR
     procedure dense_to_CSR
     procedure dense_to_COO
     procedure  COO_to_CSR
     !procedure LCO_to_CSR
  end interface spmat_convert

  public :: COO_to_dense
contains

  !-----------------------------------------------------------------------
  !> @brief dense matrix to LIL matrix 
  !> @param [in] mat: 2D array
  !> @param [out] ll : LIL matrix
  !-----------------------------------------------------------------------
  subroutine dense_to_LIL(mat, ll)
    ! check shape
    real(dp), intent(in):: mat(:, :)
    type(LIL_mat_t), intent(inout) :: ll
    integer :: s(2), irow, icol
    s=shape(mat)
    call ll%initialize(s)
    do icol=1, ll%ncol, 1
       do irow=1, ll%nrow, 1
          call ll%insert(irow,icol,mat(irow, icol),0)
       enddo
    enddo
  end subroutine dense_to_LIL

  !-----------------------------------------------------------------------
  !> @brief LIL matrix to dense matrix
  !> @param [in] ll : LIL matrix
  !> @param [out] mat: 2D array
  !-----------------------------------------------------------------------
  subroutine LIL_to_dense(ll, mat)
    class(LIL_mat_t) , intent(inout):: ll
    real(dp), intent(out):: mat(ll%nrow,ll%ncol)
    integer:: irow
    mat(:,:)=0.0d0
    do irow=1, ll%nrow
       call llist_iter_restart(ll%rows(irow))
       do while(associated(ll%rows(irow)%iter))
          mat(irow, ll%rows(irow)%iter%i)=ll%rows(irow)%iter%val
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
    enddo
  end subroutine LIL_to_dense

  !-----------------------------------------------------------------------
  !> @brief LIL matrix to COO matrix
  !> @param [in] ll : LIL matrix
  !> @param [out] csrmat: COO matrix
  !-----------------------------------------------------------------------
  subroutine LIL_to_COO(ll, COO)
    class(LIL_mat_t) , intent(inout):: ll
    class(COO_mat_t), intent(out):: COO
    integer ::  irow
    call  COO%initialize(mshape=ll%mshape)
    do irow=1, ll%nrow
       call llist_iter_restart(ll%rows(irow))
       do while(associated(ll%rows(irow)%iter))
          call coo%add_entry([irow, ll%rows(irow)%iter%i], ll%rows(irow)%iter%val)
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
    end do

  end subroutine LIL_to_COO

  !-----------------------------------------------------------------------
  !> @brief LIL matrix to CSR matrix
  !> @param [in] ll : LIL matrix
  !> @param [out] csrmat: CSR matrix
  !-----------------------------------------------------------------------
  subroutine LIL_to_CSR(ll, csrmat)
    type(LIL_mat_t) , intent(inout):: ll
    type(CSR_mat_t), intent(inout):: csrmat
    integer:: irow,  i, nzrow
    call  csrmat%initialize([ll%nrow, ll%ncol])
    call csrmat%set(nnz=ll%get_nnz())
    i=0
    csrmat%row_shift(1)=1
    do irow=1, ll%nrow
       nzrow=0
       call llist_iter_restart(ll%rows(irow))
       do while(associated(ll%rows(irow)%iter))
          nzrow=nzrow+1
          i=i+1
          csrmat%icol(i)=ll%rows(irow)%iter%i
          csrmat%val(i)=ll%rows(irow)%iter%val
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
       csrmat%row_shift(irow+1)=csrmat%row_shift(irow)+nzrow
    enddo
  end subroutine LIL_to_CSR

  !-----------------------------------------------------------------------
  !> @brief dense matrix to CSR matrix
  !> @param [in] mat: dense matrix
  !> @param [out] csrmat: CSR matrix
  !-----------------------------------------------------------------------
  subroutine dense_to_CSR(mat, csrmat)
    real(dp), intent(in):: mat(:, :)
    type(csr_mat_t), intent(out) :: csrmat
    !type(lil_mat_t):: lilmat
    integer :: i, j, nnz, nrow, ncol, inz, nzrow
    !call dense_to_lil(mat, lilmat)
    !call lil_to_csr(lilmat, csrmat)
    !call lilmat%finalize()
    nnz=count(mat /= 0.0_dp)
    nrow=size(mat, 1)
    ncol=size(mat, 2)
    call  csrmat%initialize([nrow, ncol])
    call csrmat%set(nnz=nnz)
    csrmat%row_shift(1)=1
    inz=0
    do i=1, nrow
       nzrow=0
       do j=1, ncol
          if (mat(i, j)/=0.0_dp) then
             inz=inz+1
             nzrow=nzrow+1
             csrmat%icol(inz)=j
             csrmat%val(inz)=mat(i,j)
          end if
       end do
       csrmat%row_shift(i+1)=csrmat%row_shift(i)+nzrow
    end do
  end subroutine dense_to_CSR

  !-------------------------------------------------------------------!
  ! dense_to_coo
  ! Dense to COO matrix convertion
  !-------------------------------------------------------------------!
  subroutine dense_to_coo(mat, coo)
    real(dp), intent(in) :: mat(:,:)
    type(coo_mat_t), intent(inout) :: coo
    integer:: i, j
    call coo%initialize(shape(mat))
    do j=1, size(mat, dim=2)
       do i=1, size(mat, dim=1)
          if (abs(mat(i, j)) >tiny(1.0d0) ) then
             call coo%add_entry(ind=[i,j], val=mat(i,j))
          end if
       end do
    end do
  end subroutine dense_to_coo

  !-------------------------------------------------------------------!
  ! COO_to_CSR:
  !  translate COO matrix to CSR matrix
  !  NOTE: This can be quite slow when it is large due to the sort algorithm 
  ! Input:
  !  COO matrix
  ! Output:
  !   CSR matrix
  !-------------------------------------------------------------------!
  subroutine COO_to_CSR(coo, csr)
    type(COO_mat_t), intent(inout) :: coo
    type(CSR_mat_t), intent(inout) :: csr
    integer :: ngroup
    integer, allocatable :: i1_list(:), istartend(:)
    integer :: row_nz(coo%mshape(1))
    integer :: i, irow
    call coo%group_by_1dim(ngroup, i1_list, istartend)
    call csr%initialize(coo%mshape)
    call csr%set(nnz=coo%nnz)
    csr%icol(:)=coo%ind%data(2,1:coo%nnz)
    csr%val(:) = coo%val%data(1:coo%nnz)
    csr%row_shift(:)=0
    row_nz(:)=0
    do i=1, ngroup
       irow=i1_list(i)
       row_nz(irow) = istartend(i+1)-istartend(i)
    end do

    csr%row_shift(1)=1
    do i=2, csr%nrow+1
       csr%row_shift(i)= csr%row_shift(i-1)+row_nz(i-1)
    end do
    if(allocated(i1_list)) then
       ABI_DEALLOCATE(i1_list)
    endif
    if(allocated(istartend)) then
       ABI_DEALLOCATE(istartend)
    endif
  end subroutine COO_to_CSR


  !-------------------------------------------------------------------!
  ! COO_to_dense
  !   COO matrix to dense matrix
  ! Input:
  !   COO: coo matrix
  ! Output:
  !   dense: dense matrix
  !-------------------------------------------------------------------!
  subroutine COO_to_dense(coo, dense)
    type(COO_mat_t), intent(inout) :: coo
    real(dp), intent(inout) :: dense(:,:)
    integer :: i, j, inz
    real(dp) :: val
    dense(:,:) =0.0
    do inz =1, coo%nnz
       i=coo%ind%data(1, inz)
       j=coo%ind%data(2, inz)
       val= coo%val%data(inz)
       dense(i,j) = dense(i,j) + val
    end do
  end subroutine COO_to_dense


!  subroutine LCO_to_CSR(lco, csr)
!    type(LCO_mat_t), intent(inout) :: lco
!    type(CSR_mat_t), intent(inout) :: csr
!    integer, allocatable :: i1_list(:), istartend(:)
!    integer :: row_nz(lco%mshape(1))
!    integer :: i, irow
!
!  end subroutine LCO_to_CSR



  subroutine spmat_convert_unittest()
    real(dp) ::mat(4,4), x(4), b(4)
    type(coo_mat_t) :: coo
    type(csr_mat_t) :: csr, csr2
    integer :: ngroup
    integer, allocatable :: i1list(:), ise(:)
    mat=reshape([1,2,0,0, 8, 3, 0,0, 0,2,0,0, 0,0,0,9], [4,4])*1.0d0
    call dense_to_coo(mat, coo)
    call coo%sum_duplicates()
    call dense_to_csr(mat, csr)
    call coo%group_by_1dim(ngroup, i1list, ise)
    call coo%sum_duplicates()
    if(allocated(i1list)) ABI_DEALLOCATE(i1list)
    if(allocated(ise)) ABI_DEALLOCATE(ise)
    call coo_to_csr(coo, csr2)

    x=[1,2,3,4]
    b=matmul(mat, x)
    b=0
    call csr%mv(x, b)
    b=0
    call csr2%mv(x, b)

    b=0
    call csr2%mv(x, b)



    call coo%finalize()
    call csr%finalize()
    call csr2%finalize()
  end subroutine spmat_convert_unittest


end module m_spmat_convert

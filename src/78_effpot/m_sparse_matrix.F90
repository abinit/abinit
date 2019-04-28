!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sparse_matrix
!!
!! NAME
!! m_sparse_matrix
!!
!! FUNCTION
!!  This module implement several sparse matrix data structure and matrix-vector product
!!  including LIL_mat: for constructing CSR matrix
!!  CSR_mat : for m-v product.
!!  dense matrix and COO_mat: for reference. Maybe dense matrix will be used if dipdip is considered?
!! Datatypes:
!!  * lnode: a node in linked list
!!  * llist: linked list
!!  * LIL_mat: linked list sparse matrix
!!  * dense_mat: dense matrix (to have compatible syntax with sparse matrices so it's easy to test.)
!!  * COO_mat : COO matrix
!!  * CSR_mat : CSR matrix
!! Subroutines:
!!
!!  * dense_mat_initialize
!!  * dense_mat_finalize
!!  * dense_mat_mv
!!  * dense_mat_spmv 
!!  * dense_mat_insert
!!
!!  * LIL_mat_initialize
!!  * LIL_mat_finalize
!!  * LIL_mat_insert : insert value to LIL matrix
!!  * mat_to_LIL_mat: dense->LIL
!!  * LIL_mat_to_mat: LIL->dense
!!  * LIL_mat_to_COO: LIL->COO
!!  * LIL_mat_to_CSR: LIL->CSR
!!
!!  * CSR_mat_initialize: 
!!  * CSR_mat_finalize:
!!  * CSR_mat_mv : matrix-vector product
!!
!!  * COO_mat_initialize
!!  * COO_mat_finalize TODO hexu: not yet!
!!  * COO_mat_mv
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


! This file implement several sparse matrix format and matrix vector product.
! Dense matrix is also here for test.
! Dense_mat
! LIL_mat
! COO_mat
! CSR_mat
! LIL mat is use only for constructing matrix, therefore there is no mat-vec production.
! No good to use COO in this implementation. (perhaps)
! CSR mat is better used for matrix-vector product.
! Dense_mat is better if matrix is dense or small.
! Dense<->LIL->CSR and LIL->COO translations are possible.

! A distributed counterpart of the module is m_sparse_matrix_distributed. (TODO: To be implemented)
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module m_sparse_matrix

  use defs_basis
  use m_errors
  use m_abicore

  implicit none

  !!***

  ! node of linked list, which will be one non-zero entry in LIL matrix
  type lnode
     integer :: i
     real(dp):: val
     type(lnode), pointer :: next
  end type lnode

  ! linked list of (i, val), it can be a column or a row of LIL sparse matrix
  type llist
     type(lnode), pointer :: first=>null()
     type(lnode), pointer :: last=>null()
     type(lnode), pointer :: iter=>null()
     integer :: length =0
  end type llist

  ! linked list type sparse matrix
  ! a array of rows, each row is a linked list.
  ! used for constructing sparse matrix. Not to do calculations.
  type LIL_mat
     integer :: nrow=0, ncol=0
     type(llist), allocatable :: rows(:)
   !contains
   !  procedure :: print => LIL_mat_print
  end type LIL_mat

  ! dense matrix
  type dense_mat
     integer ::  nrow, ncol
     real(dp), allocatable :: mat(:,:)
  end type dense_mat

  ! COO sparse matrix.
  ! i, j, val are the row index, col index and value of each entry.
  ! nnz: number of non-zeros.
  type COO_mat
     integer :: nnz, nrow, ncol
     integer, allocatable :: i(:), j(:)
     real(dp), allocatable:: val(:)
  end type COO_mat

  ! CSR sparse matrix
  ! nnz: number of non-zeros.
  ! icol : column index of entries .size:nnz
  ! row_shift: row_shift(irow) to row_shift(irow+1)-1  are the index
  !            of column and values for entries. size : nrow+1
  ! val: values of non-zero entries. size(nnz)
  type CSR_mat
     integer :: nnz,  nrow, ncol
     integer, allocatable :: icol(:), row_shift(:)
     real(dp), allocatable:: val(:)
   !contains
   !  procedure :: print=> CSR_mat_print
  end type CSR_mat

  ! linked list format sparse matrix
contains
  subroutine dense_mat_initialize(self,  nrow, ncol)

    type(dense_mat), intent(inout) :: self
    integer, intent(in):: nrow, ncol
    self%nrow=nrow
    self%ncol=ncol
    ABI_MALLOC(self%mat, (nrow, ncol))
    self%mat(:,:)=0.0d0
  end subroutine dense_mat_initialize

  subroutine dense_mat_finalize(self)

    type(dense_mat), intent(inout) :: self
    ABI_SFREE(self%mat)
    self%ncol=0
    self%nrow=0
  end subroutine dense_mat_finalize

  subroutine dense_mat_insert(self, irow, icol, val)

    type(dense_mat), intent(inout) :: self
    integer, intent(inout) :: irow, icol
    real(dp), intent(in) :: val
    self%mat(irow, icol)=val
  end subroutine dense_mat_insert

  ! dense matrix-vector multiplication, using blas DGEMV
  subroutine dense_mat_mv(self, x, b)

    type(dense_mat) :: self
    real(dp), intent(in) :: x
    real(dp), intent(inout) :: b
    call dgemv("N", self%nrow, self%ncol, 1.0d0,self%mat , 2,  x, 1, 0.0d0,  b, 1)
  end subroutine dense_mat_mv

  ! dense symmetric matrix-vector multiplication, using blas DSPMV
  subroutine dense_mat_spmv(self, x, b)

    type(dense_mat) :: self
    real(dp), intent(in) :: x
    real(dp), intent(inout) :: b
    !  dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    call dspmv('U', self%nrow, 1.0d0,self%mat , x, 1, 0.0d0, b, 1)
  end subroutine dense_mat_spmv

  ! COO matrix
  subroutine COO_mat_initialize(A, nrow, ncol, nnz, i, j, val)

    type(COO_mat), intent(inout) :: A
    integer, intent(in) :: nnz, nrow, ncol
    integer , intent(in), optional :: i(:), j(:)
    real(dp), intent(in), optional:: val(:)
    ABI_MALLOC(A%i, (nnz))
    ABI_MALLOC(A%j, (nnz))
    ABI_MALLOC(A%val, (nnz))
    A%nrow=nrow
    A%ncol=ncol
    A%nnz=nnz
    if(present(i)) then
    A%i(:)=i(:)
end if

    if(present(j)) then
    A%j(:)=j(:)
end if
    if(present(val)) then
    A%val(:)=val(:)
endif
  end subroutine COO_mat_initialize

  ! COO sparse matrix-vector multiplication. naive implementation.
  subroutine COO_mat_mv(A, x, b)

    type(COO_mat), intent(in) :: A
    real(dp), intent(in):: x(:)
    real(dp), intent(inout):: b(:)
    integer:: ind, ind_i, ind_j
    b(:)=0.0D0
    do ind = 1, A%nnz, 1
       ind_i=A%i(ind)
       ind_j=A%j(ind)
       b(ind_i)=b(ind_i)+A%val(ind)*x(ind_j)
    end do
  end subroutine  COO_mat_mv

  subroutine LIL_mat_initialize(self, nrow, ncol)

    type(LIL_mat) , intent(inout):: self
    integer, intent(in):: nrow, ncol
    self%nrow=nrow
    self%ncol=ncol
    ABI_ALLOCATE(self%rows, (nrow))
  end subroutine LIL_mat_initialize

  subroutine LIL_mat_finalize(self)

    type(LIL_mat) , intent(inout):: self
    integer :: i

    if (allocated(self%rows)) then
      do i=1, self%nrow, 1
        call llist_finalize(self%rows(i))
      end do
      ABI_FREE(self%rows)
    endif
    self%ncol=0
    self%nrow=0
  end subroutine LIL_mat_finalize

  subroutine LIL_mat_insert(self, irow, icol, val, mode)

    type(LIL_mat) , intent(inout):: self
    integer, intent(in):: irow, icol, mode
    real(dp), intent(in):: val
    if(abs(val)>tiny(0.0d0)) then
       call llist_sorted_insert(self%rows(irow), icol, val, mode)
    end if
  end subroutine LIL_mat_insert

  subroutine mat_to_LIL_mat(mat, ll)
    ! check shape
    real(dp), intent(in):: mat(:, :)
    type(LIL_mat), intent(inout) :: ll
    integer :: s(2), irow, icol
    s=shape(mat)
    !if ( s(1) /= ll%nrow .or. s(2)/=ll%ncol ) then
       !print *, "mat_to_LIL_mat: nrow or ncol not equal."
    !endif
    do icol=1, ll%ncol, 1
       do irow=1, ll%nrow, 1
          call LIL_mat_insert(ll,irow,icol,mat(irow, icol),0)
       enddo
    enddo
  end subroutine mat_to_LIL_mat

  subroutine LIL_mat_print(self, mat)

    class(LIL_mat) , intent(inout)::self 
    real(dp), intent(out):: mat(self%nrow,self%ncol)
    integer:: irow !, icol
    !real(dp):: val
    mat(:,:)=0.0d0
    do irow=1, self%nrow
       call llist_iter_restart(self%rows(irow))
       do while(associated(self%rows(irow)%iter))
          !write(std_out,*) "Irow: " ,irow, "Icol: ", self%rows(irow)%iter%i, "  val: ", self%rows(irow)%iter%val
          !TODO print with wrtout
          self%rows(irow)%iter=>self%rows(irow)%iter%next
       enddo
    enddo
  end subroutine LIL_mat_print

  
  subroutine LIL_mat_to_mat(ll, mat)

    type(LIL_mat) , intent(inout):: ll
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
  end subroutine LIL_mat_to_mat


  function LIL_mat_get_nnz(ll) result(nnz)

    type(LIL_mat), intent(in)::ll
    integer ::nnz, irow
    nnz=0
    do irow=1, ll%nrow
       nnz=nnz+ll%rows(irow)%length
    enddo
  end function LIL_mat_get_nnz

  subroutine LIL_mat_to_COO(ll, COO)

    type(LIL_mat) , intent(inout):: ll
    type(COO_mat), intent(out):: COO
    integer :: nnz, counter, irow
    
    nnz=LIL_mat_get_nnz(ll)
    
    call  COO_mat_initialize(COO, ll%nrow, ll%ncol, LIL_mat_get_nnz(ll))
    counter=0
    do irow=1, ll%nrow
       call llist_iter_restart(ll%rows(irow))
    do while(associated(ll%rows(irow)%iter))
        counter=counter+1
          !write(std_out,*) "I: ", self%iter%i, "  val: ", self%iter%val
          !mat(irow, ll%rows(irow)%iter%i)=ll%rows(irow)%iter%val
          !ll%rows(irow)%iter=>ll%rows(irow)%iter%next
          !COO_mat_insert(COO, irow, ll%rows(irow)%iter%i, ll%rows(irow)%iter%val )
          COO%i(counter)=irow
          COO%j(counter)=ll%rows(irow)%iter%i
          COO%val(counter)= ll%rows(irow)%iter%val 
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
       end do

  end subroutine LIL_mat_to_COO

  subroutine LIL_mat_to_CSR(ll, csrmat)

    type(LIL_mat) , intent(inout):: ll
    type(CSR_mat), intent(out):: csrmat
    integer:: irow, i, nzrow !, nnz
    !CSR_mat_initialize(A,nrow,ncol,nnz,i,j,val)
    call  CSR_mat_initialize(csrmat, ll%nrow, ll%ncol, LIL_mat_get_nnz(ll))
    i=0
    csrmat%row_shift(1)=1
    do irow=1, ll%nrow
       nzrow=0
       call llist_iter_restart(ll%rows(irow))
       do while(associated(ll%rows(irow)%iter))
          nzrow=nzrow+1
          i=i+1
          !mat(irow, ll%rows(irow)%iter%i)=ll%rows(irow)%iter%val
          !ll%rows(irow)%iter=>ll%rows(irow)%iter%next
          csrmat%icol(i)=ll%rows(irow)%iter%i
          csrmat%val(i)=ll%rows(irow)%iter%val
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
       csrmat%row_shift(irow+1)=csrmat%row_shift(irow)+nzrow
    enddo

  end subroutine LIL_mat_to_CSR

  subroutine CSR_mat_print(csr)

    class(CSR_mat), intent(in) :: csr
    ABI_UNUSED(csr%val)
    !print *, "Val: ", csr%val
    !print *, "iCol: ", csr%icol
    !print *, "row_shift:", csr%row_shift
  end subroutine CSR_mat_print

  recursive subroutine llist_finalize(self)

    type(llist), intent(inout) ::self
    type(lnode), pointer :: iter, tmp
    iter=>self%first
    do while(associated(iter))
       tmp=>iter
       iter=>iter%next
       if( associated(tmp)) then
         ABI_FREE_SCALAR(tmp)
       endif
    enddo
    nullify(self%first)
    nullify(self%last)
    self%length=0
  end subroutine llist_finalize

  subroutine llist_append(self, i, val)
    ! append a element at the end of list
    type(llist), intent(inout) ::self
    integer, intent(in)::i
    real(dp), intent(in)::val
    if(.not. associated(self%last)) then
       ABI_MALLOC_SCALAR(self%first)
       self%last=>self%first
    else
       ABI_MALLOC_SCALAR(self%last%next)
       self%last=>self%last%next
    endif
    self%last%i=i
    self%last%val=val
    self%last%next=>null()
    self%length = self%length+1
    !print *, 'length', self%length
  end subroutine llist_append

  subroutine llist_iter_restart(self)

    type(llist):: self
    self%iter=>self%first
  end subroutine llist_iter_restart

  subroutine llist_insert_after(self, ptr, i, val)
    !insert a element so i is sorted.
    ! if mode=0: if i already exist, substitute i, val
    ! if mode=1: val+=val
    type(llist):: self
    integer, intent(in) :: i
    real(dp), intent(in):: val
    type(lnode), pointer, intent(in):: ptr
    type(lnode), pointer:: tmp=>null()
    if(.not.associated(ptr%next)) then
       call llist_append(self,i,val)
    else
       ABI_MALLOC_SCALAR(tmp)
       tmp%i=i
       tmp%val=val
       tmp%next=>ptr%next
       ptr%next=>tmp
       self%length=self%length+1
    endif
  end subroutine llist_insert_after


  subroutine llist_insert_head(self, i, val)

    type(llist):: self
    integer, intent(in) :: i
    real(dp), intent(in):: val
    type(lnode), pointer:: tmp=>null()
    ABI_MALLOC_SCALAR(tmp)
    tmp%i=i
    tmp%val=val
    tmp%next=>self%first
    self%first=>tmp
    if (self%length==0) then
       self%last=>tmp
    endif
    self%length=self%length+1
  end subroutine llist_insert_head

  subroutine llist_sorted_insert(self, i, val, mode)
    !insert a element so i is sorted.
    ! if mode=0: if i already exist, substitute i, val
    ! if mode=1: val+=val
    type(llist):: self
    integer, intent(in) :: i, mode
    real(dp), intent(in):: val
    !print *, "debug insert"
    !print *, "first i", self%first%i
    !print *, "last i", self%last%i
    !print *, i
    call llist_iter_restart(self)
    if(.not.associated(self%last)) then
       !print *, "append"
       ! no element in list
       call llist_append(self,i,val)
       !print *, self%last%i
       !print *, self%last%val
    else if (i<self%first%i) then
       !print *, "insert head"
       call llist_insert_head(self, i, val)
    else
       !print *, "insert middle"
       do while(associated(self%iter))
          ! at the begining i<i0
          ! before the end,
          if (i>self%iter%i .and. associated(self%iter%next) .and. i<self%iter%next%i) then
             call llist_insert_after(self,self%iter,i,val)
             return
          else if(i==self%iter%i) then
             ! i<i0 or i>i
             if(mode==0) then
                self%iter%val=val
             else if(mode==1) then
                self%iter%val=self%iter%val+val
             endif
             return
          endif
          self%iter=>self%iter%next
       enddo
       ! i>last i
       if(i>self%last%i) then
          call llist_append(self,i,val)
          return
       else
          write(std_out, *) "cannot find proper place to insert"
       endif
    endif

    !allocate(self%last%next)
    !self%last=>self%last%next

  end subroutine llist_sorted_insert

  subroutine llist_print_all(self)

    type(llist), intent(inout)::self
    call llist_iter_restart(self)
    write(std_out,*) "linkedlist of length ", self%length
    do while(associated(self%iter))
       write(std_out,*) "I: ", self%iter%i, "  val: ", self%iter%val
       self%iter=>self%iter%next
    enddo
  end subroutine llist_print_all


  subroutine llist_get_data(self, ilist, vallist)

    type(llist), intent(inout)::self
    integer, allocatable, intent(inout)::ilist(:)
    real(dp),allocatable, intent(inout)::vallist(:)
    integer::ind=1
    ABI_MALLOC(ilist, (self%length))
    ABI_MALLOC(vallist, (self%length))
    call llist_iter_restart(self)
    do while(associated(self%iter))
       ilist(ind)=self%iter%i
       vallist(ind)=self%iter%val
       self%iter=>self%iter%next
       ind=ind+1
    enddo
  end subroutine llist_get_data

  subroutine CSR_mat_finalize(A)

    type(CSR_mat), intent(inout) :: A
    A%ncol=0
    A%nrow=0
    A%nnz=0
    ABI_SFREE(A%icol)
    ABI_SFREE(A%row_shift)
    ABI_SFREE(A%val)
  end subroutine CSR_mat_finalize


  ! COO matrix
  subroutine CSR_mat_initialize(A, nrow, ncol, nnz, icol, row_shift, val)

    type(CSR_mat), intent(inout) :: A
    integer, intent(in) :: nnz, nrow, ncol
    ! i: col number of each entry
    ! j: first 0, n1, n1+n2, ...
    ! val(irow, irow+1) are the values of entries in row irow.
    integer , intent(in), optional :: icol(:), row_shift(:)
    real(dp), intent(in), optional:: val(:)

    if(.not. allocated(A%icol)) then
       ABI_ALLOCATE(A%icol, (nnz))
    endif

    if(.not. allocated(A%row_shift)) then
       ABI_ALLOCATE(A%row_shift, (nrow+1))
    endif

    if(.not. allocated(A%val)) then
       ABI_ALLOCATE(A%val, (nnz))
    endif

    A%nrow=nrow
    A%ncol=ncol
    A%nnz=nnz
    if(present(icol)) then
       A%icol(:)=icol(:)
    endif
    if(present(row_shift)) then
       A%row_shift(:)=row_shift(:)
    endif
    if(present(val)) then
       A%val(:)=val(:)
    endif
  end subroutine CSR_mat_initialize

  subroutine CSR_mat_mv(A, x, y)

    type(CSR_mat), intent(in):: A
    real(dp), intent(in) :: x(A%ncol)
    real(dp), intent(out) :: y(A%nrow)
    integer::irow, i1, i2, i
    !real(dp)::ddot
    !external ddot
    y(:)=0.0d0
    !$OMP PARALLEL DO private(i, i1, i2)
    do irow=1, A%nrow

        ! benchmark: this version is the second fastest
        !y(irow)=dot_product(A%val(A%row_shift(irow):A%row_shift(irow+1)-1),  x(A%icol(A%row_shift(irow):A%row_shift(irow+1)-1)))

        ! third fastest
        !y(irow)=ddot(A%row_shift(irow+1)-A%row_shift(irow), A%val(A%row_shift(irow):A%row_shift(irow+1)-1),1, &
        !x(A%icol(A%row_shift(irow):A%row_shift(irow+1)-1)), 1)

        ! Slowest
        ! Do not use this one. cost twice time
        !i1=A%row_shift(irow)
        !i2=A%row_shift(irow+1)-1
        !y(irow)=dot_product(A%val(i1:i2), x(A%icol(i1:i2)))


        ! Comment hexu: The Champion of Speed, however, it relies on the optimization of compilers. need at least -O2 with both ifort &
        ! gfortran. It may also depend on SIMD instructions the CPU support. Tested on a cpu with AVX2.
        i1=A%row_shift(irow)
        i2=A%row_shift(irow+1)-1
        do i=i1, i2
            y(irow)=y(irow)+ A%val(i)*x(A%icol(i))
        end do

        ! Same speed as previous 
        !do i=A%row_shift(irow), A%row_shift(irow+1)-1
        !    y(irow)=y(irow)+ A%val(i)*x(A%icol(i))
        !end do

    enddo
    !$OMP END PARALLEL DO
  end subroutine CSR_mat_mv

!#ifdef DMKL
! wrapper to mkl CSR matrix mv mkl_dcsrgemv. For test only.
!  subroutine CSR_mat_mv_mkl(A, x, y)
!    type(CSR_mat), intent(in):: A
!    real(dp), intent(in) :: x(A%ncol)
!    real(dp), intent(out) :: y(A%nrow)
!    !call mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
!    call mkl_dcsrgemv('N', A%nrow, A%val, A%row_shift, A%icol, x, y)
!  end subroutine CSR_mat_mv_mkl
!#endif

end module m_sparse_matrix

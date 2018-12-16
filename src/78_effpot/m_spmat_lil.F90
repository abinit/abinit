#include "abi_common.h"

module m_spmat_lil
  use defs_basis
  use m_spmat_base
  use m_linked_list
  implicit none
  private

  ! linked list type sparse matrix
  ! a array of rows, each row is a linked list.
  ! used for constructing sparse matrix. Not to do calculations.
  type,public, extends(base_mat_t):: LIL_mat_t
     type(llist), allocatable :: rows(:)
   contains
     procedure :: initialize => lil_mat_t_initialize
     procedure :: finalize => lil_mat_t_finalize
     procedure :: insert => lil_mat_t_insert
     procedure :: get_nnz => lil_mat_t_get_nnz
     procedure :: print => LIL_mat_t_print
    end type LIL_mat_t

  contains
  !------------------ LIL ---------------------
  subroutine LIL_mat_t_initialize(self, nrow, ncol)
    class(LIL_mat_t) , intent(inout):: self
    integer, intent(in):: nrow, ncol
    self%nrow=nrow
    self%ncol=ncol
    ABI_ALLOCATE(self%rows, (nrow))
  end subroutine LIL_mat_t_initialize

  subroutine LIL_mat_t_finalize(self)
    class(LIL_mat_t) , intent(inout):: self
    integer :: i
    if (allocated(self%rows)) then
    do i=1, self%nrow, 1
          call llist_finalize(self%rows(i))
    end do
    ABI_DEALLOCATE(self%rows)
    endif
    self%ncol=0
    self%nrow=0
  end subroutine LIL_mat_t_finalize

  subroutine LIL_mat_t_insert(self, irow, icol, val, mode)
    class(LIL_mat_t) , intent(inout):: self
    integer, intent(in):: irow, icol, mode
    real(dp), intent(in):: val
    if(abs(val)>tiny(0.0d0)) then
       call llist_sorted_insert(self%rows(irow), icol, val, mode)
    end if
  end subroutine LIL_mat_t_insert


  subroutine LIL_mat_t_print(self, mat)

    class(LIL_mat_t) , intent(inout)::self 
    real(dp), intent(out):: mat(self%nrow,self%ncol)
    integer:: irow, icol
    real(dp):: val
    mat(:,:)=0.0d0
    do irow=1, self%nrow
       call llist_iter_restart(self%rows(irow))
       do while(associated(self%rows(irow)%iter))
          !print*, "Irow: " ,irow, "Icol: ", self%rows(irow)%iter%i, "  val: ", self%rows(irow)%iter%val
          !TODO print with wrtout
          self%rows(irow)%iter=>self%rows(irow)%iter%next
       enddo
    enddo
  end subroutine LIL_mat_t_print

 
  function LIL_mat_t_get_nnz(ll) result(nnz)

    class(LIL_mat_t), intent(in)::ll
    integer ::nnz, irow
    nnz=0
    do irow=1, ll%nrow
       nnz=nnz+ll%rows(irow)%length
    enddo
  end function LIL_mat_t_get_nnz


end module m_spmat_lil

#include "abi_common.h"
module m_spmat_convert
  use defs_basis
  use m_spmat_dense
  use m_linked_list
  use m_spmat_lil
  use m_spmat_csr
  use m_spmat_coo
  implicit none
  public
contains
  subroutine dense_to_LIL(mat, ll)
    ! check shape
    real(dp), intent(in):: mat(:, :)
    type(LIL_mat_t), intent(inout) :: ll
    integer :: s(2), irow, icol
    s=shape(mat)
    call ll%initialize(nrow=s(1), ncol=s(2))
    do icol=1, ll%ncol, 1
       do irow=1, ll%nrow, 1
          call ll%insert(irow,icol,mat(irow, icol),0)
       enddo
    enddo
  end subroutine dense_to_LIL

  subroutine LIL_to_dense(ll, mat)

    class(LIL_mat_t) , intent(inout):: ll
    real(dp), intent(out):: mat(ll%nrow,ll%ncol)
    integer:: irow, icol
    real(dp):: val
    mat(:,:)=0.0d0
    do irow=1, ll%nrow
       call llist_iter_restart(ll%rows(irow))
       do while(associated(ll%rows(irow)%iter))
          mat(irow, ll%rows(irow)%iter%i)=ll%rows(irow)%iter%val
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
    enddo
  end subroutine LIL_to_dense


  subroutine LIL_to_COO(ll, COO)
    class(LIL_mat_t) , intent(inout):: ll
    class(COO_mat_t), intent(out):: COO
    integer ::  counter, irow
    call  COO%initialize(ll%nrow, ll%ncol, ll%get_nnz())
    counter=0
    do irow=1, ll%nrow
       call llist_iter_restart(ll%rows(irow))
       do while(associated(ll%rows(irow)%iter))
          counter=counter+1
          COO%i(counter)=irow
          COO%j(counter)=ll%rows(irow)%iter%i
          COO%val(counter)= ll%rows(irow)%iter%val 
          ll%rows(irow)%iter=>ll%rows(irow)%iter%next
       enddo
    end do

  end subroutine LIL_to_COO

  subroutine LIL_to_CSR(ll, csrmat)
    type(LIL_mat_t) , intent(inout):: ll
    type(CSR_mat_t), intent(inout):: csrmat
    integer:: irow, icol, i, nzrow, nnz
    real(dp):: val
    call  csrmat%initialize(ll%nrow, ll%ncol, ll%get_nnz())
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
  end subroutine LIL_to_CSR

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
    call  csrmat%initialize(nrow, ncol, nnz)

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

end module m_spmat_convert

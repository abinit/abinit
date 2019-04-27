
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module m_mathfuncs

  use defs_basis, only: dp, PI, std_out
  use m_errors
  use m_abicore
  !use ziggurat
  use m_random_xoroshiro128plus
  implicit none

contains

  ! vector cross production
  function cross(a, b) result(c)

    real(dp), intent(in) :: a(3), b(3)
    real(dp)  :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    return
  end function cross


  ! defines outer product of two vectors.
  function outer_product(a,b) ! result (c)

    real(dp), intent(in) :: a(:), b(:)
    real(dp)  :: outer_product(size(b, dim=1), size(a, dim=1))
    integer:: i, j
      do i=1, size(a, dim=1)
        do j=1, size(b, dim=1)
          outer_product(j,i) = a(i)*b(j)
      enddo
    enddo
    return
  end function outer_product

  subroutine set_random_seed(seed)

      integer , intent(in) :: seed(:)
      write(std_out,*) "Warning! Currently I'm not sure about how this function  &
      &(set_random_seed,which calls RANDOM_SEED) works. Do test it!"
      call RANDOM_SEED(put=seed(:))
  end subroutine set_random_seed

  ! Random number generator; Normal (Gaussian) dist.
  ! a is a array.
  subroutine rand_normal_builtin(a)

    real(dp), intent(out)::a(:,:)
    real(dp), allocatable :: b(:,:)
    ABI_ALLOCATE(b, (size(a,dim=1), size(a, dim=2)))
    call random_number(a)
    b(:,:) = sqrt(-2*dlog(1.0-a(:,:)))
    call random_number(a)
    a(:,:)=b(:,:)*cos(PI*a(:,:))
    ABI_DEALLOCATE(b)
    !a(:,:)=1.0
    !a(:,1)=0.2
  end subroutine rand_normal_builtin

  ! wrapper of  normal random number generator using ziggurat method.
  ! rnor implemented in ziggruat.F90
  ! TODO hexu: implement a random number generator module which contains
  ! all the related subroutines. & Move this one there.
  !  subroutine rand_normal_ziggurat(a)


!    real(dp), intent(out)::a(:,:)
!    integer :: i, j
!    do i=1, size(a, dim=2)
!      do j=1, size(a, dim=1)
!        a(j, i)=rnor()
!      end do
!    end do
!    end subroutine rand_normal_ziggurat

end module m_mathfuncs

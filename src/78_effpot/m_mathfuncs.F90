#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module m_mathfuncs

  use defs_basis, only: dp, PI
  use m_errors
  use m_abicore
  !use ziggurat
  use m_random_xoroshiro128plus
  implicit none

  interface mat33det
     procedure  real_mat33det
     procedure  int_mat33det
  end interface mat33det


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
      print*, "Warning! Currently I'm not sure about how this function  &
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


  function real_mat33det(A) result(det)
    real(dp), intent(in) :: A(3,3)
    real(dp) :: det
    DET =   A(1,1)*A(2,2)*A(3,3)  &
         - A(1,1)*A(2,3)*A(3,2)  &
         - A(1,2)*A(2,1)*A(3,3)  &
         + A(1,2)*A(2,3)*A(3,1)  &
         + A(1,3)*A(2,1)*A(3,2)  &
         - A(1,3)*A(2,2)*A(3,1)
  end function real_mat33det

  function int_mat33det(A) result(det)
    integer, intent(in) :: A(3,3)
    integer :: det
    DET =   A(1,1)*A(2,2)*A(3,3)  &
         - A(1,1)*A(2,3)*A(3,2)  &
         - A(1,2)*A(2,1)*A(3,3)  &
         + A(1,2)*A(2,3)*A(3,1)  &
         + A(1,3)*A(2,1)*A(3,2)  &
         - A(1,3)*A(2,2)*A(3,1)
  end function int_mat33det


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

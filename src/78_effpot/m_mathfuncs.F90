#include "abi_common.h"
module m_mathfuncs
  use defs_basis, only: dp, PI
  use ziggurat
  implicit none
  !real(dp), parameter :: PI=4.D0*DATAN(1.D0)
contains

  ! vector cross production
  function cross(a, b) result(c)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cross'
!End of the abilint section

    real(dp), intent(in) :: a(3), b(3)
    real(dp)  :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    return
  end function cross


  ! defines outer product of two vectors.
  function outer_product(a,b) ! result (c)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outer_product'
!End of the abilint section

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
      ! FIXME Not sure how does this work. One number for each thread?

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'set_random_seed'
!End of the abilint section

      integer , intent(in) :: seed(:)
      print*, "Warning! Currently I'm not sure about how this function  &
      &(set_random_seed,which calls RANDOM_SEED) works. Do test it!"
      call RANDOM_SEED(put=seed(:))
  end subroutine set_random_seed

  ! Random number generator; Normal (Gaussian) dist.
  ! a is a array.
  subroutine rand_normal(a)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rand_normal'
!End of the abilint section

    real(dp), intent(out)::a(:,:)
    real(dp), allocatable, save :: b(:,:)
    if (.not. allocated(b)) then
    ABI_ALLOCATE(b, (size(a,dim=1), size(a, dim=2)))
    ! TODO: this will not be deallocated. should implement a RNG module,
    ! and this will be reimplemented.
    end if
    call random_number(a)
    b(:,:) = sqrt(-2*dlog(1.0-a(:,:)))
    call random_number(a)
    a(:,:)=b(:,:)*cos(PI*a(:,:))
    !ABI_DEALLOCATE(b)
    !a(:,:)=1.0
    !a(:,1)=0.2
  end subroutine rand_normal

  ! wrapper of  normal random number generator using ziggurat method.
  ! rnor implemented in ziggruat.F90
  ! TODO hexu: implement a random number generator module which contains 
  ! all the related subroutines. & Move this one there.
  subroutine rand_normal_ziggurat(a)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rand_normal_ziggurat'
!End of the abilint section

    real(dp), intent(out)::a(:,:)
    integer :: i, j
    do i=1, size(a, dim=2)
      do j=1, size(a, dim=1)
        a(j, i)=rnor()
      end do
    end do
    end subroutine rand_normal_ziggurat
  
end module

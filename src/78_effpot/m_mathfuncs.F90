
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module m_mathfuncs

  use defs_basis, only: dp, PI, std_out
  use m_errors
  use m_abicore
  use m_random_xoroshiro128plus
  implicit none


  ! the determinant of a 3*3 matrix
  interface mat33det
     procedure  real_mat33det
     procedure  int_mat33det
  end interface mat33det

  ! integer/real
  ! change a diagonal (a array) to a matrix
  ! or get the diagonal of a 2D matrix
  ! similar to matlab diag function
  interface diag
     procedure diag_mat_int
     procedure diag_mat_real
     procedure diag_array_int
     procedure diag_array_real
  end interface diag

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

  !----------------------------------------------------------------------
  !> @brief compare two arrays a<b
  !>
  !> @param[in]  a: integer array
  !> @param[in]  b: integer array
  !> @param[in]  N: length of a and b
  !> @return bool.
  !----------------------------------------------------------------------
  function array_lessthan(a, b, N) result (y)
    integer, intent(in)  :: a(:), b(:), N
    logical :: y
    integer :: i
    y=.False.
    do i =1, N
       if (a(i).lt.b(i)) then
          y=.True.
          exit
       elseif (a(i).gt.b(i)) then
          y=.False.
          exit
       end if
    end do
  end function array_lessthan


  !----------------------------------------------------------------------
  !> @brief compare two arrays a>b
  !>
  !> @param[in]  a: integer array
  !> @param[in]  b: integer array
  !> @param[in]  N: length of a and b
  !> @return bool.
  !----------------------------------------------------------------------
  function array_morethan(a, b, N) result (y)
    integer, intent(in)  :: a(:), b(:), N
    logical :: y
    integer :: i
    y=.False.
    do i =1, N
       if (a(i).gt.b(i)) then
          y=.True.
          exit
       elseif (a(i).lt.b(i)) then
          y=.False.
          exit
       end if
    end do
  end function array_morethan


  !----------------------------------------------------------------------
  !> @brief find an integer from a array
  !>
  !> @param[in]  a: the array to find from
  !> @param[in] x: the value to be find
  !> @param[out] ix: the index found.  0 if not found.
  !----------------------------------------------------------------------
  function find_int(a, x) result(ix)
    integer, intent(in):: a(:), x
    integer :: ix, i
    ix=0
    do i=1, size(a)
       if( a(i)==x ) then
          ix=i
       endif
    end do
  end function find_int

  !----------------------------------------------------------------------
  !> @brief find an integer from a SORTED array using binary search
  !>
  !> @param[in]  a: the array to find from
  !> @param[in] x: the value to be find
  !> @param[out] ix: the index found.  0 if not found.
  !----------------------------------------------------------------------
  function binsearch_left_integer(a, x) result(ix)
    integer, intent(in):: a(:), x
    integer :: n,ix, ub, lb
    integer , save :: i=1
    n=size(a)
    if (i<0 .or. i>n) i=(size(a)/2+1)
    if (a(i)==x) then
       ix=i
    else
       ub=n
       lb=1
       do while (lb<ub)
          i=floor((lb+ub)/2.0)
          if (a(i)< x) then
             lb=i+1
          else
             ub=i
          end if
       end do
       if(a(lb)==x) then
          ix=lb
          i=ix
       else
          ix=0
       end if
    endif
  end function binsearch_left_integer

  !-------------------------------------------------------------------!
  !Binaray search in a interger list.
  !  Once it find one match, the index is returned
  ! Input:
  !  a: the list.
  !  x: the element to search for
  ! Output:
  !  ix: the index of x. If x is not in a, ix =0. 
  !-------------------------------------------------------------------!
  function binsearch_left_integerlist(a, x) result(ix)
    integer, intent(in):: a(:,:), x(:)
    integer :: n,ix, ub, lb, nx
    integer , save :: i=1
    nx=size(x)
    n=size(a, dim=2)
    if (i<0 .or. i>n) i=(size(a, dim=2)/2+1)
    if (all(a(:,i)==x(:))) then
       ix=i
    else
       ub=n
       lb=1
       do while (lb<ub)
          i=floor((lb+ub)/2.0)
          if (array_lessthan(a(:, i), x, nx)) then
             lb=i+1
          else
             ub=i
          end if
       end do
       if(all(a(:, lb)==x)) then
          ix=lb
          i=ix
       else
          ix=0
       end if
    endif
  end function binsearch_left_integerlist

  subroutine set_random_seed(seed)

      integer , intent(in) :: seed(:)
      write(std_out,*) "Warning! Currently I'm not sure about how this function  &
      &(set_random_seed,which calls RANDOM_SEED) works. Do test it!"
      call RANDOM_SEED(put=seed(:))
  end subroutine set_random_seed

  !----------------------------------------------------------------------
  !> @brief get the diagonal of a matrix
  !>
  !> @param[in]  mat: matrix
  !> @param[out] ret: the diagonal, a 1d array
  !----------------------------------------------------------------------
  pure function diag_mat_real(mat) result (ret)
    real(dp), intent(in) :: mat(:, :)
    real(dp):: ret(size(mat, dim=1))
    integer :: n, i
    n=size(mat, dim=1)
    do i=1, n
       ret(i)=mat(i,i)
    end do
  end function diag_mat_real

  !----------------------------------------------------------------------
  !> @brief get the diagonal of a integer matrix
  !>
  !> @param[in]  mat: matrix
  !> @param[out] ret: the diagonal, a 1d array
  !----------------------------------------------------------------------
  pure function diag_mat_int(mat) result (ret)
    integer, intent(in) :: mat(:, :)
    integer:: ret(size(mat, dim=1))
    integer :: n, i
    n=size(mat, dim=1)
    do i=1, n
       ret(i)=mat(i,i)
    end do
  end function diag_mat_int

  !----------------------------------------------------------------------
  !> @brief build a matrix from the diagonal (int)
  !>
  !> @param[in]  a: the diagonal array
  !> @param[out] the matrix
  !----------------------------------------------------------------------
  pure function diag_array_int(a) result (ret)
    integer, intent(in) :: a(:)
    integer:: ret(size(a), size(a))
    integer :: i
    ret(:,:)=0
    do i=1, size(a)
       ret(i, i)=a(i)
    end do
  end function diag_array_int

  !----------------------------------------------------------------------
  !> @brief build a matrix from the diagonal (real)
  !>
  !> @param[in]  a: the diagonal array
  !> @param[out] the matrix
  !----------------------------------------------------------------------
  pure function diag_array_real(a) result (ret)
    real(dp), intent(in) :: a(:)
    real(dp):: ret(size(a), size(a))
    integer :: i
    ret(:,:)=0.0_dp
    do i=1, size(a)
       ret(i, i)=a(i)
    end do
  end function diag_array_real



  ! Random number generator; Normal (Gaussian) dist.
  ! This should NOT be used in practice.
  ! There is no standart builtin random number in fortran compiler,
  ! which could be a very bad one.
  ! also the algorithm used here is not efficient.
  ! Only for test purpose
  subroutine rand_normal_builtin(a)

    real(dp), intent(out)::a(:,:)
    real(dp), allocatable :: b(:,:)
    ABI_ALLOCATE(b, (size(a,dim=1), size(a, dim=2)))
    call random_number(a)
    b(:,:) = sqrt(-2*dlog(1.0-a(:,:)))
    call random_number(a)
    a(:,:)=b(:,:)*cos(PI*a(:,:))
    ABI_DEALLOCATE(b)
  end subroutine rand_normal_builtin

  !-------------------------------------------------------------------!
  ! real_mat33det
  ! 3*3 real(dp) matrix determinant.
  ! Input:
  !  A: 3*3 real matrix
  ! Output:
  !  det: the determinant
  !-------------------------------------------------------------------!
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

  !-------------------------------------------------------------------!
  ! int_mat33det
  ! 3*3 integer matrix determinant.
  ! Input:
  !  A: 3*3 interger matrix
  ! Output:
  !  det: the determinant
  !-------------------------------------------------------------------!
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


  !-------------------------------------------------------------------!
  ! Shortcut for Hermitian matrix eigen value and eigen vectors.
  ! Input:
  !   evecs
  ! Output:
  !  evals: eigen values
  !  evecs: eigen vectors
  !-------------------------------------------------------------------!
  subroutine eigensh(evals, evecs)
    real(dp), intent(inout) :: evals(:)
    complex(dp), intent(inout) :: evecs(:,:)
    complex(dp), allocatable  :: work(:)
    integer ::  info, lwork
    real(dp) ::  Rwork(3*size(evecs,1)-2)
    integer :: ndim
    external ZHEEV
    ndim = size(evecs, 1)
    lwork = -1
    ABI_ALLOCATE(work , (1))
    call ZHEEV('V', 'U', ndim, evecs, ndim, evals, work, lwork, rwork, info)
    lwork= INT(work(1))
    ABI_SFREE(work)

    ABI_ALLOCATE(work , (lwork))
    call ZHEEV('V', 'U', ndim, evecs, ndim, evals, work, lwork, rwork, info)
    ABI_SFREE(work)
    IF( INFO.gt.0 ) THEN
       MSG_ERROR('The zheev algorithm failed to compute eigenvalues.')
       STOP
    END IF
  end subroutine eigensh

  !-----------------------------------------------------------------------
  !> @brief rotate vector around axis by angle
  !>  Using the quaternion rotation algorithm
  !>   \vec{v}_new = \vec{v} + 2 \vec{r} \cross (\vec{r}\cross\vec{v}
  !>                  + w \vec{v})
  !> see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
  !> @param [in] angle : angle
  !> @param [in] axis  : axis vector. The norm of axis does not matter.
  !> @param [in] vec  : vector to roate
  !> @param [out] vec2 : result vector
  !-----------------------------------------------------------------------
  function rotate_by_angle_around_axis(angle, axis, vec) result(vec2)
    real(dp), intent(in) :: angle, axis(3), vec(3)
    real(dp) :: vec2(3)
    real(dp) :: half_angle, r(3), w, norm
    half_angle=angle/2.0_dp
    norm=sqrt(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3))
    r(:)=axis(:)/norm * sin(half_angle)
    w=cos(half_angle)
    ! (w, r) is the quaternion
    vec2(:) = vec2(:) + 2.0 * cross(r, (cross(r, vec) + w*vec))
  end function rotate_by_angle_around_axis

end module m_mathfuncs

!!****m* ABINIT/m_polynomial
!! NAME
!! m_polynomial
!!
!! FUNCTION
!! This module contains the datastructure of a polynomial
!!
!! Datatypes:
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


module m_polynomial
  use defs_basis
  use m_errors
  use m_abicore

  implicit none
!!***

  type, public :: polynomial_element_t
     real(dp), allocatable :: coeff
     integer :: n
     integer, allocatable :: ids(:)
     integer, allocatable :: orders(:)
   contains
     procedure :: alloc => elem_alloc
     procedure :: initialize => elem_initialize
     procedure :: finalize => elem_finalize
     procedure :: evaluate => elem_evaluate
     procedure :: add_to_first_derivative => elem_add_to_first_derivative
  end type polynomial_element_t

  type, public :: polynomial_t
     integer:: size=0, capacity=0
     type(polynomial_element_t) , allocatable :: terms(:)
   contains
     procedure :: push => polynomial_push
     procedure :: finalize => polynomial_finalize
     procedure :: evaluate => polynomial_evaluate
     procedure :: add_to_first_derivative=> polynomial_add_to_first_derivative
  end type polynomial_t

  public :: test_polynomial

contains

  subroutine elem_alloc(self, n)
    class(polynomial_element_t) ,intent(inout) :: self
    integer, intent(in) :: n
    ABI_ALLOCATE(self%ids, (n))
    ABI_ALLOCATE(self%orders, (n))
    self%n=n
  end subroutine elem_alloc

  subroutine elem_initialize(self, n,  coeff, ids, orders)
    class(polynomial_element_t) ,intent(inout) :: self
    real(dp), intent(in) :: coeff
    integer, intent(in) :: n, ids(n), orders(n)
    call self%alloc(n)
    self%coeff = coeff
    self%ids = ids
    self%orders = orders
  end subroutine elem_initialize

  subroutine elem_finalize(self)
    class(polynomial_element_t) ,intent(inout) :: self
    self%n=0
    ABI_SFREE(self%ids)
    ABI_SFREE(self%orders)
  end subroutine elem_finalize

  function elem_evaluate(self, vals) result (ret)
    class(polynomial_element_t) ,intent(in) :: self
    real(dp), intent(in) :: vals(self%n)
    real(dp):: ret
    integer :: i
    ret = self%coeff
    do i = 1, self%n
       ret =ret * vals(self%ids(i))**self%orders(i)
    end do
  end function elem_evaluate


  subroutine elem_add_to_first_derivative(self, vals, d1)
    class(polynomial_element_t) ,intent(in) :: self
    real(dp), intent(in) :: vals(self%n)
    real(dp), intent(inout):: d1(:)
    integer :: i, j
    real(dp) :: t1
    t1 = self%coeff
    do i = 1, self%n
       t1 =t1* vals(self%ids(i))**self%orders(i)
    end do
    if (t1==0.0_dp) then
       return
    end if
    do i=1, self%n
       j=self%ids(i)
       d1(j) = d1(j) + t1 * self%orders(i) / vals(j)
    end do
  end subroutine elem_add_to_first_derivative

  subroutine polynomial_push(self, elem)
    class(polynomial_t), intent(inout):: self
    class(polynomial_element_t), intent(in) :: elem
    type(polynomial_element_t), allocatable :: temp(:)
    self%size=self%size+1
    if(self%size==1) then
       self%capacity=8
       ABI_MALLOC(self%terms, (self%capacity))
    else if ( self%size>self%capacity ) then
       self%capacity = self%size + self%size / 4 + 16
       ABI_MALLOC(temp, (self%capacity))
       temp(:self%size-1) = self%terms
       !temp gets deallocated
       ABI_MOVE_ALLOC(temp, self%terms)
    end if
    self%terms(self%size)=elem
  end subroutine polynomial_push

  subroutine polynomial_finalize(self)
    class(polynomial_t), intent(inout):: self
    integer :: i
    do i =1, self%size
       call self%terms(i)%finalize()
    end do
    ABI_SFREE(self%terms)
    self%size=0
    self%capacity=0
  end subroutine polynomial_finalize

  function polynomial_evaluate(self, vals) result (ret)
    class(polynomial_t), intent(inout):: self
    real(dp), intent(in) :: vals(:)
    real(dp) :: ret
    integer :: i
    ret=0.0_dp
    do i =1, self%size
       ret = ret + self%terms(i)%evaluate(vals)
    end do
  end function polynomial_evaluate

  subroutine polynomial_add_to_first_derivative(self, vals, d1)
    class(polynomial_t), intent(inout):: self
    real(dp), intent(in) :: vals(:)
    real(dp), intent(inout) :: d1(:)
    integer :: i
    do i=1, self%size
       call self%terms(i)%add_to_first_derivative(vals, d1)
    end do
  end subroutine polynomial_add_to_first_derivative

  subroutine test_polynomial
    real(dp) :: vals(5)=[1,0,3,4,5]
    real(dp) :: d1(5)
    type(polynomial_element_t) :: p1
    type(polynomial_t) :: p
    real(dp) :: e
    call p1%initialize(2, 0.1_dp, [1, 2,4], [2,1] )
    e=p1%evaluate(vals) 
    print *, "e: ", e
    d1=0.0
    call p1%add_to_first_derivative(vals, d1)
    !print *, "d1", d1
    call p1%finalize()


    call p1%initialize(2, 0.1_dp, [1, 2], [2,1] )
    call p%push(p1)
    call p1%finalize()

    call p1%initialize(2, 0.1_dp, [3, 4], [2,3] )
    call p%push(p1)
    call p1%finalize()

    e=p%evaluate(vals)
    !print *, "e:", e

    d1=0.0
    call p%add_to_first_derivative(vals, d1)
    !print *, "d1:", d1

  end subroutine test_polynomial


end module m_polynomial


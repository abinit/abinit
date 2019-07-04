!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dynarray
!!
!! NAME
!! m_dynarray
!!
!! FUNCTION
!! Module for int and real(dp) array which allocate memory dynamically
!! real_array_type for real(dp) and int_array_type for integer.
!! they have push (but no pop) and finalize methods.
!! TODO hexu: Is this already implemented somewhere in abinit.
!! If not, should this file be moved to the place to make it more general usable?
!!
!! MG: Yes, this module should be moved to a lower level directory so that one can reuse it in other 
!! parts of the code.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module m_dynamic_array

  use defs_basis
  use m_abicore
  use m_errors
  use m_mathfuncs, only: array_morethan, binsearch_left_integerlist,binsearch_left_integer

  implicit none
  private
!!***

!!****t* defs_abitypes/real_array_type
!! NAME
!! real_array_type
!!
!! FUNCTION
!! datatype of real(dp) array which can be dynamically allocated
!!
!! SOURCE
  type, public:: real_array_type
    integer:: size=0, capacity=0
    real(dp), allocatable :: data(:)
  CONTAINS
    procedure :: push => real_array_type_push
    procedure :: finalize => real_array_type_finalize
  end type real_array_type


!!***

!!****t* defs_abitypes/int_array_type
!! NAME
!! int_array_type
!!
!! FUNCTION
!! datatype of real(dp) array which can be dynamically allocated
!!
!! SOURCE
  type , public::int_array_type
    integer:: size=0, capacity=0
    integer, allocatable :: data(:)
  CONTAINS
    procedure :: push => int_array_type_push
    procedure :: finalize => int_array_type_finalize
    procedure :: sort => int_array_type_sort
  end type int_array_type
!!***

  !!****t* defs_abitypes/int2d_array_type
  !! NAME
  !! int2d_array_type
  !!
  !! FUNCTION
  !! datatype of integer array which the dim=2 can be dynamically allocated
  !!
  !! SOURCE
  type , public::int2d_array_type
     integer:: size=0, capacity=0
     logical :: sorted=.False.
     integer, allocatable :: data(:,:)
   CONTAINS
     procedure :: push => int2d_array_type_push
     procedure :: push_unique => int2d_array_type_push_unique
     procedure :: sort => int2d_array_type_sort
     procedure :: binsearch =>int2d_array_type_binsearch
     procedure :: finalize => int2d_array_type_finalize
  end type int2d_array_type
  !!***

  public::  dynamic_array_unittest
CONTAINS


!****f* m_dynarray/real_array_type_push
!!
!! NAME
!! real_array_type_push
!!
!! FUNCTION
!! push data to a real_array_type
!!
!! INPUTS
!! self = real_array_type object
!! val= data to be pushed
!! OUTPUT
!! real_array<type(real_array_type)()> = real_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine real_array_type_push(self, val)

    class(real_array_type), intent(inout):: self
    real(dp) :: val
    real(dp), allocatable :: temp(:)
    self%size=self%size+1
    if(self%size==1) then
      self%capacity=8
      ABI_MALLOC(self%data, (self%capacity))
    else if ( self%size>self%capacity ) then
      self%capacity = self%size + self%size / 4 + 8
      ABI_MALLOC(temp, (self%capacity))
      temp(:self%size-1) = self%data
      !temp gets deallocated
      ABI_MOVE_ALLOC(temp, self%data) 
    end if
    self%data(self%size)=val

end subroutine real_array_type_push
!!***

!****f* m_disarray/real_array_type_finalize
!!
!! NAME
!! real_array_type_finalize
!!
!! FUNCTION
!! destroy real_array_type
!!
!! INPUTS
!! self= real_array_type object
!! OUTPUT
!! real_array<type(real_array_type)()> = real_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine real_array_type_finalize(self)

  class(real_array_type), intent(inout):: self
  ABI_SFREE(self%data)
  self%size=0
  self%capacity=0

end subroutine real_array_type_finalize
!!***

!****f* m_dynarray/int_array_type_push
!!
!! NAME
!! int_array_type_push
!!
!! FUNCTION
!! push data to a int_array_type
!!
!! INPUTS
!! self = int_array_type object
!! val= data to be pushed
!! OUTPUT
!! int_array<type(real_array_type)()> = int_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine int_array_type_push(self, val)

    class(int_array_type), intent(inout):: self
    integer :: val
    integer, allocatable :: temp(:)
    self%size=self%size+1
    if(self%size==1) then
      self%capacity=8
      ABI_MALLOC(self%data, (self%capacity))
    else if ( self%size>self%capacity ) then
      self%capacity = self%size + self%size / 4 + 8
      ABI_MALLOC(temp, (self%capacity))
      temp(:self%size-1) = self%data
      !temp gets deallocated
      ABI_MOVE_ALLOC(temp, self%data) 
    end if
    self%data(self%size)=val

end subroutine int_array_type_push
!!***


!----------------------------------------------------------------------
!> @brief insertion_sort_int: sort a array using insertion sort algorithm
!>  it is a memory safe method but is generally slow.
!> @param[inout]  a: the array to be sorted. and will output inplace
!> @param[inout] order (optional) the sorted index, it can be used to sort
!>  other arrays so that the order in consistent.
!----------------------------------------------------------------------
subroutine insertion_sort_int(a, order)
  integer, intent(inout) :: a(:)
  integer, optional, intent(inout):: order(size(a))
  integer :: n,i,j, v
  n=size(a)
  if (present(order)) then
     do i = 1,n
        order(i)=i
     end do
  end if
  do i = 2,n
     v=a(i)
     j=i-1
     do while(j>=1 )
        if (a(j)<=v) exit
        a(j+1)=a(j)
        if(present(order)) order(j+1)=order(j)
        j=j-1
     end do
     a(j+1)=v
     if(present(order)) order(j+1)=i
  end do

end subroutine insertion_sort_int

!----------------------------------------------------------------------
!> @brief insertion_sort_int: sort a DYNAMIC INT array using insertion sort algorithm
!>  it is a memory safe method but is generally slow.
!> @param[inout]  a: an dynamic array. the array to be sorted. and will output inplace
!> @param[inout] order (optional) the sorted index, it can be used to sort
!>  other arrays so that the order in consistent.
!----------------------------------------------------------------------
subroutine int_array_type_sort(self, order)
  class(int_array_type), intent(inout):: self
  integer, optional, intent(inout):: order(self%size)
  integer :: i,j, v
  if (present(order)) then
     do i = 1, self%size
        order(i)=i
     end do
  end if
  do i = 2, self%size
     v=self%data(i)
     j=i-1
     do while(j>=1 )
        if(.not. self%data(j)>v) exit
        self%data(j+1)=self%data(j)
        if(present(order)) order(j+1)=order(j)
        j=j-1
     end do
     self%data(j+1)=v
     if(present(order)) order(j+1)=i
  end do
end subroutine int_array_type_sort


!****f* m_dynarray/int_array_type_finalize
!!
!! NAME
!! int_array_type_finalize
!!
!! FUNCTION
!! destroy int_array_type
!!
!! INPUTS
!! self= int_array_type object
!! OUTPUT
!! int_array<type(int_array_type)()> = int_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine int_array_type_finalize(self)

  class(int_array_type), intent(inout):: self
  ABI_SFREE(self%data)
  self%size=0
  self%capacity=0

end subroutine int_array_type_finalize


!==================================================================

!****f* m_dynarray/int2d_array_type_push
!!
!! NAME
!! int2d_array_type_push
!!
!! FUNCTION
!! push data to a int2d_array_type
!!
!! INPUTS
!! self = int2d_array_type object
!! val= data to be pushed
!! OUTPUT
!! int_array<type(real_array_type)()> = int2d_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine int2d_array_type_push(self, val)

    class(int2d_array_type), intent(inout):: self
    integer :: val(:)
    integer, allocatable :: temp(:,:)
    self%size=self%size+1
    if(self%size==1) then
      self%capacity=8
      ABI_MALLOC(self%data, (size(val), self%capacity))
    else if ( self%size>self%capacity ) then
      self%capacity = self%size + self%size / 4 + 8
      ABI_MALLOC(temp, (size(val), self%capacity))
      temp(:,:self%size-1) = self%data
      ABI_MOVE_ALLOC(temp, self%data) !temp gets deallocated
    end if
    self%data(:,self%size)=val
end subroutine int2d_array_type_push
!!***



!****f* m_dynarray/int2d_array_type_push
!!
!! NAME
!! int2d_array_type_push
!!
!! FUNCTION
!! push data to a int2d_array_type
!!
!! INPUTS
!! self = int2d_array_type object
!! val= data to be pushed
!! OUTPUT
!! int_array<type(real_array_type)()> = int2d_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine int2d_array_type_push_unique(self, val, position)

    class(int2d_array_type), intent(inout):: self
    integer, intent(in) :: val(:)
    integer, optional, intent(out) :: position
    integer :: i
    logical :: inside
    inside=.False.
    do i=1, self%size
       if(all(self%data(:,i)==val)) then
          inside=.True.
          if (present(position)) position=i
          exit
       endif
    enddo
    if(.not. inside) then
       call self%push(val)
       if (present(position)) position=self%size
    end if
  end subroutine int2d_array_type_push_unique
!!***



!****f* m_dynarray/int2d_array_type_finalize
!!
!! NAME
!! int2d_array_type_finalize
!!
!! FUNCTION
!! destroy int2d_array_type
!!
!! INPUTS
!! self= int2d_array_type object
!! OUTPUT
!! int_array<type(int2d_array_type)()> = int2d_array_type data
!! PARENTS
!!      m_dynarray
!!
!! CHILDREN
!!
!! SOURCE
subroutine int2d_array_type_finalize(self)

  class(int2d_array_type), intent(inout):: self
  if ( allocated(self%data) ) then
      ABI_SFREE(self%data)
  end if
  self%size=0
  self%capacity=0

end subroutine int2d_array_type_finalize

!----------------------------------------------------------------------
!> @brief sort a 2D DYNAMIC INT array using insertion sort algorithm
!>  it is a memory safe method but is generally slow.
!> it compares the elements in first dimension i.e. A(:, i) and sort the second dim.
!>  The comparing is from left to right.
!> @param[inout]  a: an dynamic array. the array to be sorted. and will output inplace
!> @param[inout] order (optional) the sorted index, it can be used to sort
!>  other arrays so that the order in consistent.
!----------------------------------------------------------------------

subroutine int2d_array_type_sort(self, order)
  class(int2d_array_type), intent(inout):: self
  integer, optional, intent(inout):: order(self%size)
  integer :: i,j, v(size(self%data, dim=1))
  if (present(order)) then
     do i = 1, self%size
        order(i)=i
     end do
  end if
  do i = 2, self%size
     v(:)=self%data(:,i)
     j=i-1
     do while(j>=1)
        if (.not. (array_morethan(self%data(:,j),v, size(self%data, dim=1)))) exit
        self%data(:,j+1)=self%data(:,j)
        if(present(order)) order(j+1)=order(j)
        j=j-1
     end do
     self%data(:,j+1)=v(:)
     if(present(order)) order(j+1)=i
  end do
  self%sorted=.True.
end subroutine int2d_array_type_sort

!----------------------------------------------------------------------
!> @brief binary search 
!>
!> @param[in] self: the 2D array to be searched from
!> @param[in] val: the value to be searched
!> @param[out] i: the index of the first one found. returns 0 if not found.
!----------------------------------------------------------------------

function int2d_array_type_binsearch(self, val) result(i)
  class(int2d_array_type), intent(inout):: self
  integer, intent(inout) :: val(:)
  integer :: i
  i=binsearch_left_integerlist(self%data(:,1:self%size), val)
end function int2d_array_type_binsearch



!====================== Unit tests======================

subroutine binsearch_test()
  integer :: a(4)=[1,2,3,4], b(3,3)=reshape([0,0,0,0,1,0,1,0,0], [3,3])
  integer :: i
  i=binsearch_left_integer(a, 5)
  i=binsearch_left_integerlist(b, [0,0,0])
end subroutine binsearch_test

subroutine insertion_sort_int_test()
  integer :: a(4), order(4), b(4), a2(8)
  a=[1,5,3,4]
  b=a
  call insertion_sort_int(a, order)
  a2=[3,6,2,4, 3, 5, 0, 9]
  call insertion_sort_int(a2)
end subroutine insertion_sort_int_test

subroutine int2d_array_test()
  type(int2d_array_type) :: t
  call t%push_unique([1,1,2])
  call t%push_unique([1,2,2])
  call t%push_unique([1,1,2])
  call t%push_unique([1,1,1])
  call t%push_unique([-1, 3, 3])
  call t%push_unique([2,1, 4])
  call t%sort()
end subroutine int2d_array_test

subroutine dynamic_array_unittest()
  call binsearch_test()
  call insertion_sort_int_test()
  call int2d_array_test()
end subroutine dynamic_array_unittest

end module m_dynamic_array
!!***

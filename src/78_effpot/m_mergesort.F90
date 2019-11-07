!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mergesort
!!
!! NAME
!! m_mergesort
!!
!! FUNCTION
!! Module for sorting integer arrays using merge sorting algorithm
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


module m_mergesort
!!***
  use m_mathfuncs, only: array_lessthan, array_morethan, array_le
  implicit none
  private
  public :: MergeSort
  public :: MergeSort2D
contains

!======== Merge sort algorithm for 1d array ===============
! Modified from https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
! The original code is public domain.

  subroutine merge(A, B, C)
    implicit none
    ! The targe attribute is necessary, because A .or. B might overlap with C.
    integer, target, intent(in) :: A(:), B(:)
    integer, target, intent(inout) :: C(:)
    integer :: i, j, k

    if (size(A) + size(B) > size(C)) then
       stop (1)
    end if

    i = 1; j = 1
    do k = 1, size(C)
       if (i <= size(A) .and. j <= size(B)) then
          if (A(i) <= B(j)) then
             C(k) = A(i)
             i = i + 1
          else
             C(k) = B(j)
             j = j + 1
          end if
       else if (i <= size(A)) then
          C(k) = A(i)
          i = i + 1
       else if (j <= size(B)) then
          C(k) = B(j)
          j = j + 1
       end if
    end do
  end subroutine merge

  subroutine merge_with_order(A, B, C, orderA, orderB, orderC)
    implicit none
    ! The targe attribute is necessary, because A .or. B might overlap with C.
    integer, target, intent(in) :: A(:), B(:), orderA(:), orderB(:)
    integer, target, intent(inout) :: C(:), orderC(:)
    integer :: i, j, k

    if (size(A) + size(B) > size(C)) then
       stop (1)
    end if

    i = 1; j = 1
    do k = 1, size(C)
       if (i <= size(A) .and. j <= size(B)) then
          if (A(i) <= B(j)) then
             C(k) = A(i)
             orderC(k) = orderA(i)
             i = i + 1
          else
             C(k) = B(j)
             orderC(k) =orderB(j)
             j = j + 1
          end if
       else if (i <= size(A)) then
          C(k) = A(i)
          orderC(k) = orderA(i)
          i = i + 1
       else if (j <= size(B)) then
          C(k) = B(j)
          orderC(k) = orderB(j)
          j = j + 1
       end if
    end do
  end subroutine merge_with_order


  subroutine swap(x, y)
    implicit none
    integer, intent(inout) :: x, y
    integer :: tmp
    tmp = x; x = y; y = tmp
  end subroutine swap

  recursive subroutine MergeSort_no_init_order(A, work, order, worder)
    implicit none
    integer, intent(inout) :: A(:)
    integer, intent(inout) :: work(:)
    integer, optional, intent(inout):: order(size(A)), worder(size(work))
    integer :: half, n
    logical :: ordered
    n=size(A)
    ordered=.False.
    if (present(order)) then
       ordered=.True.
    end if

    half = (size(A) + 1) / 2
    if (size(A) < 2) then
       continue
    else if (size(A) == 2) then
       if (A(1) > A(2)) then
          call swap(A(1), A(2))
          if(ordered) call swap(order(1), order(2))
       end if
    else
       if(ordered)  then
          call MergeSort_no_init_order(A( : half), work, order(:half), worder)
          call MergeSort_no_init_order(A(half + 1 :), work, order(half+1:), worder)
          if (A(half) > A(half + 1)) then
             work(1 : half) = A(1 : half)
             worder(1:half) = order(1: half)
             call merge_with_order(work(1 : half), A(half + 1:), A &
                  &, worder(1 : half), order(half + 1:), order)
          endif
       else
          call MergeSort_no_init_order(A( : half), work)
          call MergeSort_no_init_order(A(half + 1 :), work)
          if (A(half) > A(half + 1)) then
             work(1 : half) = A(1 : half)
             call merge(work(1 : half), A(half + 1:), A)
          endif
       endif
    end if
  end subroutine MergeSort_No_Init_Order

  subroutine MergeSort(A, work, order, worder)
    implicit none
    integer, intent(inout) :: A(:)
    integer, intent(inout) :: work(:)
    integer, optional, intent(inout):: order(size(A)), worder(size(work))
    integer :: i
    if (present(order)) then
       do i =1 , size(A)
          order(i) = i
       end do
    end if
    call MergeSort_no_init_order(A, work, order, worder)
  end subroutine MergeSort

!== Merge sort for 2D integer array=================================


  subroutine merge2D(A, B, C)
    implicit none
    ! The targe attribute is necessary, because A .or. B might overlap with C.
    integer, target, intent(in) :: A(:, :), B(:, :)
    integer, target, intent(inout) :: C(:,:)
    integer :: i, j, k

    if (size(A, 2) + size(B, 2) > size(C, 2)) then
       stop (1)
    end if

    i = 1; j = 1
    do k = 1, size(C, 2)
       if (i <= size(A, 2) .and. j <= size(B,2)) then
          if (array_le (A(:, i) ,  B(:,j), size(A, 1))) then
             C(:,k) = A(:,i)
             i = i + 1
          else
             C(:,k) = B(:, j)
             j = j + 1
          end if
       else if (i <= size(A, 2)) then
          C(:,k) = A(:, i)
          i = i + 1
       else if (j <= size(B, 2)) then
          C(:,k) = B(:,j)
          j = j + 1
       end if
    end do
  end subroutine merge2D

  subroutine merge2D_with_order(A, B, C, orderA, orderB, orderC)
    implicit none
    ! The targe attribute is necessary, because A .or. B might overlap with C.
    integer, target, intent(in) :: A(:, :), B(:,:), orderA(:), orderB(:)
    integer, target, intent(inout) :: C(:,:), orderC(:)
    integer :: i, j, k

    if (size(A, 2) + size(B,2) > size(C,2)) then
       stop (1)
    end if

    i = 1; j = 1
    do k = 1, size(C, 2)
       if (i <= size(A, 2) .and. j <= size(B, 2)) then
          if (array_le (A(:, i) ,  B(:,j), size(A, 1))) then
             C(:,k) = A(:,i)
             orderC(k) = orderA(i)
             i = i + 1
          else
             C(:,k) = B(:,j)
             orderC(k) =orderB(j)
             j = j + 1
          end if
       else if (i <= size(A,2)) then
          C(:,k) = A(:,i)
          orderC(k) = orderA(i)
          i = i + 1
       else if (j <= size(B, 2)) then
          C(:,k) = B(:,j)
          orderC(k) = orderB(j)
          j = j + 1
       end if
    end do
  end subroutine merge2D_with_order


  subroutine swap2D(x, y)
    implicit none
    integer, intent(inout) :: x(:), y(:)
    integer :: tmp, i
    do i =1, size(x)
       tmp = x(i); x(i) = y(i); y(i) = tmp
    end do
  end subroutine swap2D

  recursive subroutine MergeSort2D_no_init_order(A, work, order, worder)
    implicit none
    integer, intent(inout) :: A(:, :)
    integer, intent(inout) :: work(:, :)
    integer, optional, intent(inout):: order(size(A, 2)), worder(size(work, 2))
    integer :: half, n
    logical :: ordered
    n=size(A,2)
    ordered=.False.
    if (present(order)) then
       ordered=.True.
    end if

    half = (n + 1) / 2
    if (n < 2) then
       continue
    else if (n == 2) then
       if (array_morethan(A(:, 1) , A(:,2), size(A, 1))) then
          call swap2D(A(:,1), A(:,2))
          if(ordered) call swap(order(1), order(2))
       end if
    else
       if(ordered)  then
          call MergeSort2D_no_init_order(A(:, : half), work, order(:half), worder)
          call MergeSort2D_no_init_order(A(:, half + 1 :), work, order(half+1:), worder)
          if (array_morethan(A(:,half) , A(:,half + 1), size(A, 1))) then
             work(:, 1 : half) = A(:, 1 : half)
             worder(1:half) = order(1: half)
             call merge2D_with_order(work(:, 1 : half), A(:, half + 1:), A &
                  &, worder(1 : half), order(half + 1:), order)
          endif
       else
          call MergeSort2D_no_init_order(A(:, : half), work)
          call MergeSort2D_no_init_order(A(:, half + 1 :), work)
          if (array_morethan(A(:,half) , A(:,half + 1), size(A, 1))) then
             work(:, 1 : half) = A(:,1 : half)
             call merge2D(work(:, 1 : half), A(:, half + 1:), A)
          endif
       endif
    end if
  end subroutine MergeSort2D_No_Init_Order

  subroutine MergeSort2D(A, work, order, worder)
    implicit none
    integer, intent(inout) :: A(:,:)
    integer, intent(inout) :: work(:, :)
    integer, optional, intent(inout):: order(size(A, 2)), worder(size(work, 2))

    integer :: i
    if(present(order)) then
       do i =1 , size(A, 2)
          order(i) = i
       end do
    end if
    call MergeSort2D_no_init_order(A, work, order, worder)
  end subroutine MergeSort2D



end module m_mergesort

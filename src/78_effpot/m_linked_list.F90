!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_linked_list
!! NAME
!! m_linked_list
!!
!! FUNCTION
!! This module contains the a linked list for real(dp).
!! It is used to build the LIL (list of linked list) format of sparse matrix.
!!
!! Datatypes:
!!
!! * lnode: a node in the linked list
!! * llist: linked list
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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

module m_linked_list
  use defs_basis
  implicit none
!!***

  ! node of linked list, which will be one non-zero entry in LIL matrix
  type, public:: lnode
     integer :: i
     real(dp):: val
     type(lnode), pointer :: next
  end type lnode

  ! linked list of (i, val), it can be a column or a row of LIL sparse matrix
  type, public:: llist
     type(lnode), pointer :: first=>null()
     type(lnode), pointer :: last=>null()
     type(lnode), pointer :: iter=>null()
     integer :: length =0
     contains
         procedure :: finalize=>llist_finalize
         procedure :: append => llist_append
         procedure :: iter_restart => llist_iter_restart
         procedure :: insert_after => llist_insert_after
         procedure :: insert_head => llist_insert_head
         procedure :: sorted_insert => llist_sorted_insert
         procedure :: get_data => llist_get_data
  end type llist

  contains

  recursive subroutine llist_finalize(self)
    class(llist), intent(inout) ::self
    type(lnode), pointer :: iter, tmp
    iter=>self%first
    do while(associated(iter))
       tmp=>iter
       iter=>iter%next
       if( associated(tmp)) then
          deallocate(tmp)
       endif
    enddo
    nullify(self%first)
    nullify(self%last)
    self%length=0
  end subroutine llist_finalize

  subroutine llist_append(self, i, val)
    ! append a element at the end of list
    class(llist), intent(inout) ::self
    integer, intent(in)::i
    real(dp), intent(in)::val
    if(.not. associated(self%last)) then
       allocate(self%first)
       self%last=>self%first
    else
       allocate(self%last%next)
       self%last=>self%last%next
    endif
    self%last%i=i
    self%last%val=val
    self%last%next=>null()
    self%length = self%length+1
    !print *, 'length', self%length
  end subroutine llist_append

  subroutine llist_iter_restart(self)

    class(llist):: self
    self%iter=>self%first
  end subroutine llist_iter_restart

  subroutine llist_insert_after(self, ptr, i, val)
    !insert a element so i is sorted.
    ! if mode=0: if i already exist, substitute i, val
    ! if mode=1: val+=val
    class(llist):: self
    integer, intent(in) :: i
    real(dp), intent(in):: val
    type(lnode), pointer, intent(in):: ptr
    type(lnode), pointer:: tmp=>null()
    if(.not.associated(ptr%next)) then
       call llist_append(self,i,val)
    else
       allocate(tmp)
       tmp%i=i
       tmp%val=val
       tmp%next=>ptr%next
       ptr%next=>tmp
       self%length=self%length+1
    endif
  end subroutine llist_insert_after


  subroutine llist_insert_head(self, i, val)
    class(llist):: self
    integer, intent(in) :: i
    real(dp), intent(in):: val
    type(lnode), pointer:: tmp=>null()
    allocate(tmp)
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
    class(llist):: self
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
          print*, "cannot find proper place to insert"
       endif
    endif

    !allocate(self%last%next)
    !self%last=>self%last%next

  end subroutine llist_sorted_insert

  subroutine llist_print_all(self)

    class(llist), intent(inout)::self
    call llist_iter_restart(self)
    print*, "linkedlist of length ", self%length
    do while(associated(self%iter))
       print*, "I: ", self%iter%i, "  val: ", self%iter%val
       self%iter=>self%iter%next
    enddo
  end subroutine llist_print_all


  subroutine llist_get_data(self, ilist, vallist)

    class(llist), intent(inout)::self
    integer, allocatable, intent(inout)::ilist(:)
    real(dp),allocatable, intent(inout)::vallist(:)
    integer::ind=1
    allocate(ilist(self%length))
    allocate(vallist(self%length))
    call llist_iter_restart(self)
    do while(associated(self%iter))
       ilist(ind)=self%iter%i
       vallist(ind)=self%iter%val
       self%iter=>self%iter%next
       ind=ind+1
    enddo
  end subroutine llist_get_data

end module m_linked_list

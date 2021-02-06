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
!! * lnode: a node in the linked list, in each node, there is an integer and a real(dp)
!! * llist: linked list of lnode
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

module m_linked_list
  use defs_basis
  use m_errors
  use m_xmpi
  use m_profiling_abi
  implicit none
!!***

  ! node of linked list, which will be used as one non-zero entry in LIL matrix
  ! each node has an integer (i) and a real(dp) (val).
  ! That is because in a LIL sparse matrix, it saves the information of a
  ! column/row indices (int) and values (real)
  ! and the pointer points to the next node
  type, public:: lnode
     integer :: i
     real(dp):: val
     type(lnode), pointer :: next=>null()
  end type lnode

  ! linked list of (i, val), it can be a column or a row of LIL sparse matrix
  type, public:: llist
     type(lnode), pointer :: first=>null() ! pointer to first node
     type(lnode), pointer :: last=>null()  ! pointer to last node
     type(lnode), pointer :: iter=>null()  ! pointer to the current node (e.g. in a loop)
     integer :: length =0                  ! number of nodes
   contains
     procedure :: finalize=>llist_finalize      ! free memory
     procedure :: append => llist_append        ! add a new entry
     procedure :: iter_restart => llist_iter_restart
     procedure :: insert_after => llist_insert_after
     procedure :: insert_head => llist_insert_head
     procedure :: sorted_insert => llist_sorted_insert
     procedure :: get_data => llist_get_data
  end type llist

  contains

    !----------------------------------------------------------------------
    !> @brief finalize: free memory
    !----------------------------------------------------------------------
  recursive subroutine llist_finalize(self)
    class(llist), intent(inout) ::self
    type(lnode), pointer :: iter, tmp
    iter=>self%first
    do while(associated(iter))
       tmp=>iter
       iter=>iter%next
       if( associated(tmp)) then
          ABI_FREE_SCALAR(tmp)
       endif
    enddo
    nullify(self%first)
    nullify(self%last)
    self%length=0
  end subroutine llist_finalize

  !----------------------------------------------------------------------
  !> @brief append to the end of the list
  !> @param[in]  i: int 
  !> @param[in]  val: real value
  !----------------------------------------------------------------------
  subroutine llist_append(self, i, val)
    ! append a element at the end of list
    class(llist), intent(inout) ::self
    integer, intent(in)::i
    real(dp), intent(in)::val
    if(.not. associated(self%last)) then
       ABI_MALLOC_SCALAR(self%first)
       self%last=>self%first
    else
       ABI_MALLOC_SCALAR(self%last%next)
       self%last=>self%last%next
    endif
    self%last%i=i
    self%last%val=val
    self%last%next=>null()
    self%length = self%length+1
  end subroutine llist_append

  !----------------------------------------------------------------------
  !> @brief restart the iteration. set iter to first node
  !----------------------------------------------------------------------
  subroutine llist_iter_restart(self)
    class(llist):: self
    self%iter=>self%first
  end subroutine llist_iter_restart

  !----------------------------------------------------------------------
  !> @brief insert to the node after one node ptr
  !>
  !> @param[in] ptr: a pointer to a node in the list
  !> @param[in] i: the integer 
  !> @param[in] val: the real value
  !----------------------------------------------------------------------
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
       ABI_MALLOC_SCALAR(tmp)
       tmp%i=i
       tmp%val=val
       tmp%next=>ptr%next
       ptr%next=>tmp
       self%length=self%length+1
    endif
  end subroutine llist_insert_after

!----------------------------------------------------------------------
  !> @brief inset to the head of the list
  !>
  !> @param[in]  i
  !> @param[in] val
  !----------------------------------------------------------------------
  subroutine llist_insert_head(self, i, val)
    class(llist):: self
    integer, intent(in) :: i
    real(dp), intent(in):: val
    type(lnode), pointer:: tmp=>null()
    ABI_MALLOC_SCALAR(tmp)
    tmp%i=i
    tmp%val=val
    tmp%next=>self%first
    self%first=>tmp
    if (self%length==0) then
       self%last=>tmp
    endif
    self%length=self%length+1
  end subroutine llist_insert_head

  !----------------------------------------------------------------------
  !> @brief insert to a node so that the list is sorted by i
  !>
  !> @param[in]  i :  integer value
  !> @param[in]  val : real value
  !> @param[in]  mode : 
  !> if mode=0: if i already exist, substitute i, val
  !> if mode=1: val+=val
  !----------------------------------------------------------------------
  subroutine llist_sorted_insert(self, i, val, mode)
    !insert a element so i is sorted.
    class(llist):: self
    integer, intent(in) :: i, mode
    real(dp), intent(in):: val
    call llist_iter_restart(self)
    if(.not.associated(self%last)) then
       call llist_append(self,i,val)
    else if (i<self%first%i) then
       call llist_insert_head(self, i, val)
    else
       do while(associated(self%iter))
          ! at the begining i<i0
          ! before the end,
          if (i>self%iter%i) then
             if (associated(self%iter%next)) then
                 if (i<self%iter%next%i) then
                    call llist_insert_after(self,self%iter,i,val)
                    return
                 end if
             end if
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
          ABI_BUG("m_linked_list cannot find proper place to insert")
       endif
    endif


  end subroutine llist_sorted_insert


  !----------------------------------------------------------------------
  !> @brief get all the data to arrays of i and val
  !> @param[out]  ilist: the array of i
  !> @param[out]  vallist: the array of val
  !----------------------------------------------------------------------
  subroutine llist_get_data(self, ilist, vallist)

    class(llist), intent(inout)::self
    integer, allocatable, intent(inout)::ilist(:)
    real(dp),allocatable, intent(inout)::vallist(:)
    integer::ind=1
    ABI_MALLOC(ilist,(self%length))
    ABI_MALLOC(vallist, (self%length))
    call llist_iter_restart(self)
    do while(associated(self%iter))
       ilist(ind)=self%iter%i
       vallist(ind)=self%iter%val
       self%iter=>self%iter%next
       ind=ind+1
    enddo
  end subroutine llist_get_data

end module m_linked_list

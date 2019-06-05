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

module m_dynmaic_array

  use defs_basis
  use m_abicore
  use m_errors

  implicit none
!!***

!!****t* defs_abitypes/real_array_type
!! NAME
!! real_array_type
!!
!! FUNCTION
!! datatype of real(dp) array which can be dynamically allocated
!!
!! SOURCE
  type real_array_type
    integer:: size=0, capacity=0
    real(dp), allocatable :: data(:)
  !CONTAINS
  !  procedure :: push => real_array_type_push
  !  procedure :: finalize => real_array_type_finalize
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
  type int_array_type
    integer:: size=0, capacity=0
    integer, allocatable :: data(:)
  !CONTAINS
  !  procedure :: push => int_array_type_push
  !  procedure :: finalize => int_array_type_finalize
  end type int_array_type
!!***
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

!****f* m_dynarray/real_array_type_finalize
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

end module m_dynmaic_array
!!***

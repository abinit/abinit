! This module defines arrays which can dynamically allocate memory.
! real_array_type for real(dp) and int_array_type for integer.
! they have push (but no pop) and finalize methods.

! TODO hexu: Is this already implemented somewhere in abinit. 
! If not, should this file  be moved to the place to make it more general usable?

#include "abi_common.h"
module m_dynmaic_array
  use defs_basis
  implicit none

  type real_array_type
    integer:: size=0, capacity=0
    real(dp), allocatable :: data(:)
  !CONTAINS
  !  procedure :: push => real_array_type_push
  !  procedure :: finalize => real_array_type_finalize
  end type real_array_type

  type int_array_type
    integer:: size=0, capacity=0
    integer, allocatable :: data(:)
  !CONTAINS
  !  procedure :: push => int_array_type_push
  !  procedure :: finalize => int_array_type_finalize
  end type int_array_type

CONTAINS
subroutine real_array_type_push(self, val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'real_array_type_push'
!End of the abilint section

    class(real_array_type), intent(inout):: self
    real(dp) :: val
    real(dp), allocatable :: temp(:)
    self%size=self%size+1
    if(self%size==1) then
      self%capacity=8
      ABI_ALLOCATE(self%data, (self%capacity))
    else if ( self%size>self%capacity ) then
      self%capacity = self%size + self%size / 4 + 8
      ABI_ALLOCATE(temp,(self%capacity))
      temp(:self%size-1) = self%data
      call move_alloc(temp, self%data) !temp gets deallocated
    end if
    self%data(self%size)=val
end subroutine

subroutine real_array_type_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'real_array_type_finalize'
!End of the abilint section

  class(real_array_type), intent(inout):: self
  if ( allocated(self%data) ) then
      ABI_DEALLOCATE(self%data)
  end if
  self%size=0
  self%capacity=0
end subroutine

subroutine int_array_type_push(self, val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int_array_type_push'
!End of the abilint section

    class(int_array_type), intent(inout):: self
    integer :: val
    integer, allocatable :: temp(:)
    self%size=self%size+1
    if(self%size==1) then
      self%capacity=8
      ABI_ALLOCATE(self%data, (self%capacity))
    else if ( self%size>self%capacity ) then
      self%capacity = self%size + self%size / 4 + 8
      ABI_ALLOCATE(temp,(self%capacity))
      temp(:self%size-1) = self%data
      call move_alloc(temp, self%data) !temp gets deallocated
    end if
    self%data(self%size)=val
end subroutine int_array_type_push

subroutine int_array_type_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'int_array_type_finalize'
!End of the abilint section

  class(int_array_type), intent(inout):: self
  if ( allocated(self%data) ) then
      ABI_DEALLOCATE(self%data)
  end if
  self%size=0
  self%capacity=0
end subroutine

end module m_dynmaic_array

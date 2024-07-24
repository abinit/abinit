!
! yakl simple demo
!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program Fortran_Gator

  use defs_basis
  use gator_mod

#if defined HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_ptr,c_loc,c_intptr_t
#endif

  implicit none
  real(8), pointer, contiguous :: a(:) => null()
  integer(kind=c_intptr_t), parameter :: ptr_int = 2

  call gator_init()

  write(std_out,*) "===== before alloc ===== : c_loc(a)=",transfer(c_loc(a), ptr_int)

  if (associated(a)) then
    write(std_out,*) "a is associated"
  else
    write(std_out,*) "a is not associated"
  end if

  call gator_allocate( a , (/1024/) )

  write(std_out,*) "===== after  alloc ===== : c_loc(a)=",transfer(c_loc(a), ptr_int)

  if (associated(a)) then
    write(std_out,*) "a is associated"
  else
    write(std_out,*) "a is not associated"
  end if

  a(1) = 42.0
  write(std_out,*) "a(1)=",a(1)

  call gator_deallocate( a )

  if (associated(a)) then
    write(std_out,*) "a is associated"
  else
    write(std_out,*) "a is not associated"
  end if

  call gator_finalize()

end program Fortran_Gator

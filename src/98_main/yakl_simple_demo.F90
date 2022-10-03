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

  implicit none
  real(8), pointer, contiguous :: a(:) => null()

  call gator_init()

  write(std_out,*) "===== before alloc ===== : c_loc(a)=",c_loc(a)

  if (associated(a)) then
    write(std_out,*) "a is associated"
  else
    write(std_out,*) "a is not associated"
  end if

  call gator_allocate( a , (/1024/) )

  write(std_out,*) "===== after  alloc ===== : c_loc(a)=",c_loc(a)

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

!
!
!
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

#ifdef HAVE_KOKKOS

module m_kokkos_utils

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  implicit none

  private

  public :: &
    & abinit_kokkos_print_config, &
    & kokkos_initialize, &
    & kokkos_initialize_without_args, &
    & kokkos_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Kokkos interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface
    subroutine f_kokkos_print_config() &
      & bind(c, name="c_abinit_kokkos_print_config")
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine f_kokkos_print_config
  end interface

  interface
    subroutine f_kokkos_initialize(argc, argv) &
      & bind(c, name="c_kokkos_initialize")
      use, intrinsic :: iso_c_binding, only : c_int, c_ptr
      integer(c_int), intent(inout) :: argc
      type(c_ptr), value :: argv
    end subroutine f_kokkos_initialize
  end interface

  interface
    subroutine f_kokkos_initialize_without_args() &
      bind(c, name='c_kokkos_initialize_without_args')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine f_kokkos_initialize_without_args
  end interface

  interface
    subroutine f_kokkos_finalize() &
      & bind(c, name="c_kokkos_finalize")
    end subroutine f_kokkos_finalize
  end interface

  interface
    function f_kokkos_is_initialized() result(is_init) &
      & bind(c, name='c_kokkos_is_initialized')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(c_bool) :: is_init
    end function f_kokkos_is_initialized
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine abinit_kokkos_print_config()
    use, intrinsic :: iso_c_binding
    implicit none
    call f_kokkos_print_config()
  end subroutine abinit_kokkos_print_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kokkos_initialize()
    use, intrinsic :: iso_c_binding
    implicit none
    integer :: arg_count, max_length = 0, str_length, n, cli_count
    character(kind=c_char, len=:), allocatable :: str
    character(kind=c_char), allocatable, target :: strs_array(:,:)
    type(c_ptr), allocatable, target :: c_strs(:)

    arg_count = command_argument_count()
    ! include command name
    do n = 0, arg_count
      call get_command_argument(n, length=str_length)
      max_length = max(max_length, str_length)
    end do

    ABI_MALLOC(strs_array, (1:max_length + 1, 0:arg_count))
    ABI_MALLOC(c_strs, (0:arg_count))
    ABI_MALLOC_TYPE_SCALAR(character(max_length),str)

    do n = 0, arg_count
      call get_command_argument(n, value=str, length=str_length)
      strs_array(1:str_length, n) = transfer(str, ' ', size=str_length)
      strs_array(str_length+1, n) = c_null_char
      c_strs(n) = c_loc(strs_array(:,n))
    end do

    cli_count = arg_count + 1
    if( cli_count .le. 0 ) then
       ! must call the without_args variant or get array-oob error
       call f_kokkos_initialize_without_args()
    else
       call f_kokkos_initialize(cli_count, c_loc(c_strs(0)))
    endif

    ABI_FREE(strs_array)
    ABI_FREE(c_strs)
    ABI_FREE(str)

  end subroutine kokkos_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kokkos_initialize_without_args()
    use, intrinsic :: iso_c_binding
    implicit none
    call f_kokkos_initialize_without_args()
  end subroutine kokkos_initialize_without_args

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kokkos_finalize()
    use, intrinsic :: iso_c_binding
    implicit none
    call f_kokkos_finalize
  end subroutine kokkos_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function kokkos_is_initialized() result(is_init)
    use, intrinsic :: iso_c_binding
    implicit none
    logical :: is_init
    logical(c_bool) :: c_is_init
    c_is_init = f_kokkos_is_initialized()
    is_init = logical(c_is_init)
  end function kokkos_is_initialized

end module m_kokkos_utils

#endif

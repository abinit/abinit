
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_neat

  use defs_basis
  use m_abicore
  use m_pair_list
  use m_stream_string
  use m_yaml_out

  implicit none

  private
  !public :: neat_start_iter
  public :: stream_string

  contains

!!***f* m_neat/neat_start_iter
!!
!! NAME
!! neat_start_iter
!!
!! FUNCTION
!! Mark the beginning of an iteration
!!
!! INPUTS
!!     n integer= index of the iteration
!!     unit integer= file descriptor to write in
!!     name character(len=*)= name of the iteration (without i or n prefix)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!     Output a YAML document to start the iteration
!!
!! PARENTS
!!
!! CHILDREN
!!     yaml_iterstart
!!
!! SOURCE
!  subroutine neat_start_iter(n, name, unit)
!    integer,intent(in) :: n, unit
!    character(len=*),intent(in) :: name
!    type(stream_string) :: stream
!
!    call yaml_iterstart(trim(name), n, stream=stream)
!    call stream%flush(unit)
!  end subroutine neat_start_iter
!!***

end module m_neat
!!***

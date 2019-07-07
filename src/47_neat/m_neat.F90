
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
  public :: neat_energies, neat_start_iter
  public :: neat_etot_add_line, neat_open_etot, neat_finish_etot
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
  subroutine neat_start_iter(n, name, unit)
    integer,intent(in) :: n, unit
    character(len=*),intent(in) :: name
    type(stream_string) :: stream

    call yaml_iterstart(trim(name), n, stream=stream)
    call stream%dump(unit)
  end subroutine neat_start_iter
!!***

!!****f* m_neat/neat_energies
!!
!! NAME
!! neat_energies
!!
!! FUNCTION
!! Write components of total energies in a structured way
!!
!! INPUTS
!!  energies <type(pair_list)>=values of parts of total energy
!!  unit= unit of output file
!!  tag= optional tag to distinct the main Etot document from secondary Etot(DC) for example
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!    Print a YAML document to output file
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  subroutine neat_energies(energies, unit, tag)
    type(pair_list),intent(inout) :: energies
    integer,intent(in) :: unit
    character(len=*),intent(in),optional :: tag
!Local variables-------------------------------
    type(stream_string) :: stream

    if(present(tag)) then
      call yaml_single_dict(tag, '', energies, 35, 500, width=20, stream=stream, real_fmt='(ES25.18)')
    else
      call yaml_single_dict('EnergyTerms', '', energies, 35, 500, width=20, stream=stream, real_fmt='(ES25.18)')
    end if

    call stream%dump(unit)
  end subroutine neat_energies
!!***

!!****f* m_neat/neat_open_etot
!!
!! NAME
!! neat_open_etot
!!
!! FUNCTION
!! Open a document for ETOT
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!    Write the beginning of the document to stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  subroutine neat_open_etot(stream, comment, header)
    type(stream_string),intent(inout) :: stream
    character(len=*),intent(in) :: header, comment

    call yaml_open_doc('EtotSteps', comment, stream=stream)

    call yaml_open_tabular('data', stream=stream, tag='EtotIters')
    call yaml_add_tabular_line(header, stream=stream)

  end subroutine neat_open_etot
!!***

!!***f* m_neat/neat_etot_add_line
!!
!! NAME
!! neat_etot_add_line
!!
!! FUNCTION
!! Add a line to the ETOT table
!!
!! INPUTS
!!     stream <type(stream_string)>= stream to accumulate the document
!!     line <character(len=*)>= line to add
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  subroutine neat_etot_add_line(stream, line)
    type(stream_string),intent(inout) :: stream
    character(len=*),intent(in) :: line

    call yaml_add_tabular_line('  '//line(6:), stream=stream)
  end subroutine neat_etot_add_line
!!***

!!****f* m_neat/neat_finish_etot
!!
!! NAME
!! neat_finish_etot
!!
!! FUNCTION
!! Close the document and write the document to unit
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  subroutine neat_finish_etot(stream, unit)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: unit

    if(stream%length > 0) then
      call yaml_close_doc(stream=stream)
      call stream%dump(unit)
    end if
  end subroutine neat_finish_etot
!!***
end module m_neat
!!***

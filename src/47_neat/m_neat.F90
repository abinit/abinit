module m_neat
  use m_abicore
  use m_pair_list
  use m_stream_string
  use m_yaml_out

  implicit none

  private
  public :: neat_energies
! ! Expose pair list api
!   public :: pair_list_init, pair_list_set, pair_list_get, pair_list_free
!   public :: pair_list_next, pair_list_look, pair_list_iter, pair_list_restart
!   public :: pair_list
!   public :: TC_EMPTY, TC_NOTFOUND, TC_INT, TC_REAL, TC_STRING
! ! Expose stream string api
!   public :: stream_write, stream_get_chunk, stream_free, stream_string
!   public :: stream_to_string, stream_copy, stream_transfer, stream_to_file
!   public :: stream_debug
  contains
  
!!****f* m_neat/neat_energies
!!
!! NAME
!! neat_write_energies
!!
!! FUNCTION
!! Write components of total energies in a structured way
!!
!! INPUTS
!!  energies <type(energies_type)>=values of parts of total energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryphase
!!   | kptopt
!!   | occopt
!!   | positron=option for electron-positron calculation
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  iout= unit of output file
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
  subroutine neat_energies(energies, iout)
    type(pair_list),intent(inout) :: energies
    integer,intent(in) :: iout
!Local variables-------------------------------
    type(stream_string) :: stream
    character(len=:),allocatable :: s

    call yaml_single_dict('Etot', energies, 35, 500, key_width=20, stream=stream, real_fmt='(ES25.18)')

    allocate(character(stream%length) :: s)

    call stream_to_string(stream, s)
    call wrtout(iout, s, 'COLL')

    deallocate(s)
  end subroutine neat_energies
end module m_neat

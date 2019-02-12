module m_neat
  use defs_basis
  use m_abicore
  use m_pair_list
  use m_stream_string
  use m_yaml_out
  use m_results_gs, only: results_gs_type
  use m_crystal, only: crystal_t

  implicit none

  private
  public :: neat_energies, neat_results_gs, neat_crystal
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

  subroutine wrtout_stream(stream, iout)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: iout

    character(len=stream%length) :: s

    call stream_to_string(stream, s)
    call wrtout(iout, s, 'COLL')
  end subroutine wrtout_stream
  
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

    call yaml_single_dict('Etot', '', energies, 35, 500, tag='ETOT', width=20, stream=stream, real_fmt='(ES25.18)')

    call wrtout_stream(stream, iout)
  end subroutine neat_energies

!!****f* m_neat/neat_results_gs
!!
!! NAME
!! neat_results_gs
!!
!! FUNCTION
!! Write neat_results_gs
!!
!! INPUTS
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
  subroutine neat_results_gs(results, iout, ecut, pawecutdg, comment)
    type(results_gs_type),intent(in) :: results
    integer,intent(in) :: iout
    real(dp),intent(in) :: ecut, pawecutdg
    character(len=*),intent(in),optional :: comment

    type(stream_string) :: stream
    type(pair_list) :: dict
    real(dp) :: strten(2,3)

    if(present(comment)) then
      call yaml_open_doc('results_gs', comment, width=10, stream=stream)
    else
      call yaml_open_doc('results_gs', '', width=10, stream=stream)
    end if

    call yaml_add_intfield('natom', results%natom, width=10, stream=stream)
    call yaml_add_intfield('nsppol', results%nsppol, width=10, stream=stream)

    call pair_list_set(dict, 'ecut', r=ecut)
    call pair_list_set(dict, 'pawecutdg', r=pawecutdg)
    call yaml_add_dict('cut', dict, width=10, stream=stream)
    call pair_list_free(dict)

    call pair_list_set(dict, 'deltae', r=results%deltae)
    call pair_list_set(dict, 'res2', r=results%res2)
    call pair_list_set(dict, 'residm', r=results%residm)
    call pair_list_set(dict, 'diffor', r=results%diffor)
    call yaml_add_dict('convergence', dict, width=10, multiline_trig=2, stream=stream)
    call pair_list_free(dict)

    call yaml_add_realfield('etotal', results%etotal, width=10, stream=stream)
    call yaml_add_realfield('entropy', results%entropy, width=10, stream=stream)
    call yaml_add_realfield('fermie', results%fermie, width=10, stream=stream)

    strten(1,:) = results%strten(1:3)
    strten(2,:) = results%strten(4:6)
    call yaml_add_real2d('stress tensor', 2, 3, strten, width=10, stream=stream, tag='Tensor32')    
    call yaml_add_real2d('cartesian forces', 3, results%natom, results%fcart, width=10, stream=stream, tag='CartForces')    

    call yaml_close_doc(stream=stream)

    call wrtout_stream(stream, iout)
  end subroutine neat_results_gs

!!****f* m_neat/neat_crystal
!!
!! NAME
!! neat_crystal
!!
!! FUNCTION
!! Write neat_crystal
!!
!! INPUTS
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
  subroutine neat_crystal(crystal, iout, comment)
    type(crystal_t),intent(in) :: crystal
    integer,intent(in) :: iout
    character(len=*),intent(in),optional :: comment

    
  end subroutine neat_crystal
end module m_neat

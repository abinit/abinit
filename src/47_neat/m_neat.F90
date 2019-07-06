
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
  use m_results_gs, only: results_gs_type

  implicit none

  private
  public :: neat_energies, neat_results_gs, neat_crystal, neat_start_iter
  public :: neat_open_gw_sigma_pert, neat_gw_sigma_pert_add_line, neat_finish_gw_sigma_pert
  public :: neat_etot_add_line, neat_open_etot, neat_finish_etot
  public :: stream_string, enable_yaml

  logical :: enable = .false., switch_lock = .false.

  contains

  ! Can only be called once, then it lock. This is to prevent unconsistent
  ! behaviour with activating and disabling on demand (the parser on the python
  ! side won't like it)
  subroutine enable_yaml(yes)
    logical, intent(in) :: yes
    if (.not. switch_lock) then
      enable = yes
      switch_lock = .true.
    end if
  end subroutine enable_yaml

  subroutine wrtout_stream(stream, iout)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: iout

    character(len=stream%length) :: s

    if (enable) then
      call stream%to_string(s)
      call wrtout(iout, s, 'COLL')
    endif

    call stream%free()
  end subroutine wrtout_stream

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
!!     iout integer= file descriptor to write in
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
!!     yaml_iterstart, wrtout_stream
!!
!! SOURCE
  subroutine neat_start_iter(n, name, iout)
    integer,intent(in) :: n, iout
    character(len=*),intent(in) :: name
    type(stream_string) :: stream

    call yaml_iterstart(trim(name), n, stream=stream)
    call wrtout_stream(stream, iout)
  end subroutine neat_start_iter
!!*** m_neat/neat_start_iter

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
!!  iout= unit of output file
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
  subroutine neat_energies(energies, iout, tag)
    type(pair_list),intent(inout) :: energies
    integer,intent(in) :: iout
    character(len=*),intent(in),optional :: tag
!Local variables-------------------------------
    type(stream_string) :: stream

    if(present(tag)) then
      call yaml_single_dict(tag, '', energies, 35, 500, width=20, stream=stream, real_fmt='(ES25.18)')
    else
      call yaml_single_dict('EnergyTerms', '', energies, 35, 500, width=20, stream=stream, real_fmt='(ES25.18)')
    end if

    call wrtout_stream(stream, iout)
  end subroutine neat_energies
!!*** m_neat/neat_energies

!!****f* m_neat/neat_results_gs
!!
!! NAME
!! neat_results_gs
!!
!! FUNCTION
!! Write neat_results_gs
!!
!! INPUTS
!!  results <type(results_gs_type)>=miscellaneous informations about the system after ground state computation
!!  iout= unit of output file
!!  ecut= cut energy
!!  pawecutdg= PAW cut energy
!!  comment= optional comment for the final docuemtn
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
    real(dp) :: strten(3,3)
    real(dp) :: forces(results%natom, 3)
    integer :: j

    if(present(comment)) then
      call yaml_open_doc('ResultsGS', comment, width=10, stream=stream)
    else
      call yaml_open_doc('ResultsGS', '', width=10, stream=stream)
    end if

    call yaml_add_intfield('natom', results%natom, width=10, stream=stream)
    call yaml_add_intfield('nsppol', results%nsppol, width=10, stream=stream)

    call dict%set('ecut', r=ecut)
    call dict%set('pawecutdg', r=pawecutdg)
    call yaml_add_dict('cutoff_energies', dict, width=10, stream=stream)
    call dict%free()

    call dict%set('deltae', r=results%deltae)
    call dict%set('res2', r=results%res2)
    call dict%set('residm', r=results%residm)
    call dict%set('diffor', r=results%diffor)
    call yaml_add_dict('convergence', dict, width=10, multiline_trig=2, stream=stream)
    call dict%free()

    call yaml_add_realfield('etotal', results%etotal, width=10, stream=stream)
    call yaml_add_realfield('entropy', results%entropy, width=10, stream=stream)
    call yaml_add_realfield('fermie', results%fermie, width=10, stream=stream)

    strten(1,1) = results%strten(1)
    strten(2,2) = results%strten(2)
    strten(3,3) = results%strten(3)

    strten(2,3) = results%strten(4)
    strten(3,2) = results%strten(4)
    strten(1,3) = results%strten(5)
    strten(3,1) = results%strten(5)
    strten(1,2) = results%strten(6)
    strten(2,1) = results%strten(6)
    call yaml_add_real2d('stress tensor', 3, 3, strten, width=10, stream=stream, tag='Tensor')
    call stream%write(ch10)

    do j=1,3
      forces(:,j) = results%fcart(j,:)
    end do
    call yaml_add_real2d('cartesian forces', results%natom, 3, forces, width=10, stream=stream, tag='CartForces')

    call yaml_close_doc(stream=stream)

    call wrtout_stream(stream, iout)
  end subroutine neat_results_gs
!!***


!!****f* m_neat/neat_open_gw_sigma_pert
!!
!! NAME
!! neat_open_gw_sigma_pert
!!
!! FUNCTION
!! Open a document for GW Sigma_perturbative
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
  subroutine neat_open_gw_sigma_pert(stream, comment, k, e0, egw, degw, header)
    type(stream_string),intent(inout) :: stream
    real(kind=dp),intent(in) :: k(3), e0, egw, degw
    character(len=*),intent(in) :: header, comment

    call yaml_open_doc('GwSigmaPerturbative', comment, stream=stream)

    call yaml_add_real1d('k point', 3, k, stream=stream, real_fmt='(3f8.3)')
    call yaml_add_realfield('E^0_gap', e0, stream=stream)
    call yaml_add_realfield('E^GW_gap', egw, stream=stream)
    call yaml_add_realfield('DeltaE^GW_gap', degw, stream=stream)

    call yaml_open_tabular('data', stream=stream, tag='GwSigmaData')
    call yaml_add_tabular_line(header, stream=stream)

  end subroutine neat_open_gw_sigma_pert
!!***

!!****f* m_neat/neat_gw_sigma_pert_add_line
!!
!! NAME
!! neat_gw_sigma_pert_add_line
!!
!! FUNCTION
!! Add a line to GW Sigma_perturbative document
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
  subroutine neat_gw_sigma_pert_add_line(stream, line)
    type(stream_string),intent(inout) :: stream
    character(len=*),intent(in) :: line

    call yaml_add_tabular_line(line, stream=stream)
  end subroutine neat_gw_sigma_pert_add_line
!!***

!!****f* m_neat/neat_finish_gw_sigma_pert
!!
!! NAME
!! neat_finish_gw_sigma_pert
!!
!! FUNCTION
!! Close the document and write it to iout
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
  subroutine neat_finish_gw_sigma_pert(stream, iout)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: iout

    call yaml_close_doc(stream=stream)

    call wrtout_stream(stream, iout)
  end subroutine neat_finish_gw_sigma_pert
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
!!*** m_neat/neat_open_etot

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
!! SIDE EFFECTS
!!
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
!!*** m_neat/neat_etot_add_line

!!****f* m_neat/neat_finish_etot
!!
!! NAME
!! neat_finish_etot
!!
!! FUNCTION
!! Close the document and write the document to iout
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
  subroutine neat_finish_etot(stream, iout)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: iout

    if(stream%length > 0) then
      call yaml_close_doc(stream=stream)

      call wrtout_stream(stream, iout)
    end if
  end subroutine neat_finish_etot
!!*** m_neat/neat_finish_etot
end module m_neat
!!***

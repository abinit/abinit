!!****m* ABINIT/m_slc_dynamics
!! NAME
!! m_slc_potential
!!
!! FUNCTION
!! This module contains the dynamics for coupled spin-lattice dynamics
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu, NH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module  m_slc_dynamics

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_io_tools, only : get_unit, open_file, close_unit
  use m_mpi_scheduler, only: mpi_scheduler_t, init_mpi_info


  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_spin_mover, only: spin_mover_t
  use m_lattice_mover, only: lattice_mover_t
  use m_hashtable_strval, only: hash_table_t
  use m_slc_potential, only: slc_potential_t

  private

  type, public, extends(abstract_mover_t) :: slc_mover_t
    class(spin_mover_t),    pointer :: spin_mover
    class(lattice_mover_t), pointer :: lattice_mover    

  CONTAINS
    procedure :: initialize
    procedure :: finalize
    procedure :: run_one_step
    procedure :: run_time

  end type slc_mover_t

  contains

  subroutine initialize(self, spin_mover, lattice_mover)

    class(slc_mover_t) :: self
    class(spin_mover_t), target :: spin_mover
    class(lattice_mover_t), target :: lattice_mover

    self%spin_mover => spin_mover
    self%lattice_mover => lattice_mover

  end subroutine initialize

  subroutine finalize(self)

    class(slc_mover_t) :: self
  
    nullify(self%spin_mover)
    nullify(self%lattice_mover)

  end subroutine finalize

  subroutine run_time(self, calculator, displacement, strain, spin, lwf, energy_table)

    class(slc_mover_t),          intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: calculator

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:), spin(:,:)
    type(hash_table_t),optional, intent(inout) :: energy_table

    real(dp):: t
    integer :: counter
    character(len=80) :: msg, msg_empty
    real(dp):: etotal, espin, elatt, eslc

    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 


   !if(self%spin_mover%total_time .ne. self%lattice_mover%total_time) then
   !  MSG_ERROR("Total time for spin and lattice dynamics are different, check your input file.")
   !endif

   !if(self%spin_mover%dt .ne. self%lattice_mover%dt) then
   !  MSG_ERROR("Different time steps for spin and lattice dynamics not yet implemented, check your input file.")
   !endif
     

   t=0.0
   counter=0

   if(iam_master) then
     msg_empty=ch10

     msg=repeat("=", 80)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')
     write(msg, '(A20)') "Coupled spin-lattice dynamic steps:"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')
     msg=repeat("=", 80)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

     write(msg, "(5X, A9, 6X, A10, 5X, A10, 5X, A10, 5X, A10)")  "Iteration", "E_spin(Ha)", "E_latt(Ha)", "E_slc(Ha)", "E_tot(Ha)"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

     msg=repeat("-", 80)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')
    end if


    do while(t<self%spin_mover%total_time)
      counter=counter+1
      !one step in coupled spin-lattice dynamics
      call self%run_one_step(effpot=calculator, spin=spin, displacement=displacement, strain=strain, lwf=lwf, &
           & energy_table=energy_table)

      if (iam_master) then
        call self%spin_mover%hist%set_vars(time=t,  inc=.True.)
        call self%spin_mover%spin_ob%get_observables(self%spin_mover%hist%S(:,:, self%spin_mover%hist%ihist_prev), &
             self%spin_mover%hist%Snorm(:,self%spin_mover%hist%ihist_prev), &
             self%spin_mover%hist%etot(self%spin_mover%hist%ihist_prev))
        if(modulo(counter, self%spin_mover%hist%spin_nctime)==0) then
          call self%spin_mover%spin_ncfile%write_one_step(self%spin_mover%hist)
          etotal = energy_table%sum_val()
          espin = energy_table%sum_val(label='SpinPotential')
          elatt = energy_table%sum_val(label='Lattice_harmonic_potential')
          elatt = elatt + energy_table%sum_val(label='Lattice kinetic energy')
          eslc = energy_table%sum_val(prefix='SLCPotential')
          write(msg, "(A1, 1X, I13, 2X, ES13.5, 2X, ES13.5, 2X, ES13.5, 2X, ES13.5)") "-", counter, espin, elatt, eslc, etotal
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
        endif
      end if

      t=t+self%spin_mover%dt
    enddo

  end subroutine run_time


  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf, energy_table)

    class(slc_mover_t),          intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:), spin(:,:)
    type(hash_table_t), optional, intent(inout) :: energy_table


    call self%spin_mover%run_one_step(effpot=effpot, displacement=displacement, strain=strain, &
        &    lwf=lwf, energy_table=energy_table)

    call self%lattice_mover%run_one_step(effpot=effpot, spin=spin, lwf=lwf, energy_table=energy_table)

  end subroutine run_one_step

  !!***
end module m_slc_dynamics


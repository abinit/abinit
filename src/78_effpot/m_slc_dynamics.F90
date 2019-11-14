!{\src2tex{textfont=tt}}
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
  use m_spin_mover, only: spin_mover_t
  use m_lattice_mover, only: lattice_mover_t

  contains

  subroutine slc_run_time(spin_mover, lattice_mover, calculator, displacement, strain, spin, lwf)

    class(spin_mover_t),    intent(inout) :: spin_mover
    class(lattice_mover_t), intent(inout) :: lattice_mover
    class(abstract_potential_t), intent(inout) :: calculator

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:), spin(:,:)

    real(dp):: t
    integer :: counter
    character(len=80) :: msg, msg_empty

    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 


   if(spin_mover%total_time .ne. lattice_mover%total_time) then
     MSG_ERROR("Total time for spin and lattice dynamics are different, check your input file.")
   endif

   if(spin_mover%dt .ne. lattice_mover%dt) then
     MSG_ERROR("Different time steps for spin and lattice dynamics not yet implemented, check your input file.")
   endif
     

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

     write(msg, "(A13, 4X, A13, 6X, A13, 4X, A13)")  "Iteration", "time(s)", "Avg_Mst/Ms", "ETOT(Ha/uc)"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

     msg=repeat("-", 80)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')
    end if


    do while(t<spin_mover%total_time)
      counter=counter+1
      !one step in spin dynamics
      call spin_mover%run_one_step(effpot=calculator, displacement=displacement, strain=strain, lwf=lwf)
      if (iam_master) then
        call spin_mover%hist%set_vars(time=t,  inc=.True.)
        call spin_mover%spin_ob%get_observables(spin_mover%hist%S(:,:, spin_mover%hist%ihist_prev), &
             spin_mover%hist%Snorm(:,spin_mover%hist%ihist_prev), spin_mover%hist%etot(spin_mover%hist%ihist_prev))
        if(modulo(counter, spin_mover%hist%spin_nctime)==0) then
          call spin_mover%spin_ncfile%write_one_step(spin_mover%hist)
          write(msg, "(A1, 1X, I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") "-", counter, t*Time_Sec, &
                & spin_mover%spin_ob%Mst_norm_total/spin_mover%spin_ob%Snorm_total, &
                & spin_mover%hist%etot(spin_mover%hist%ihist_prev)/spin_mover%spin_ob%nscell
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
        endif
      end if
      !one step in lattice dynamics
      call lattice_mover%run_one_step(effpot=calculator, spin=spin, lwf=lwf)
      ! TODO: Adjust output
      !write(msg, "(A13, 4X, A15, 4X, A15, 4X, A15, 4X, A15)") &
      !      &  "Iteration", "temperature(K)", "Ekin(Ha/uc)", &
      !      & "Epot(Ha/uc)", "ETOT(Ha/uc)"
      !call wrtout(std_out,msg,'COLL')
      !call wrtout(ab_out, msg, 'COLL')

      t=t+spin_mover%dt
    enddo


  end subroutine slc_run_time

end module m_slc_dynamics


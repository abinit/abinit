!!****m* ABINIT/m_lattice_lwf_mover
!! NAME
!! m_lattice_lwf_mover
!!
!! FUNCTION
!! This module contains the mover for coupled lattice lwf dynamics 
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
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

module  m_lattice_lwf_mover

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_io_tools, only : get_unit, open_file, close_unit
  use m_mpi_scheduler, only: mpi_scheduler_t, init_mpi_info


  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_lattice_mover, only: lattice_mover_t
  use m_lwf_mover, only: lwf_mover_t
  use m_hashtable_strval, only: hash_table_t

  private

  type, public, extends(abstract_mover_t) :: lattice_lwf_mover_t
    class(lattice_mover_t), pointer :: lattice_mover => null()
    class(lwf_mover_t),    pointer :: lwf_mover => null()
  CONTAINS
    procedure :: initialize
    procedure :: finalize
    procedure :: run_time

 end type lattice_lwf_mover_t

contains
  subroutine initialize(self, lattice_mover, lwf_mover)
    class(lattice_lwf_mover_t),intent(inout) :: self
    class(lattice_mover_t), target:: lattice_mover
    class(lwf_mover_t), target:: lwf_mover
    self%lattice_mover=>lattice_mover
    self%lwf_mover=>lwf_mover
    self%dt=self%lattice_mover%dt
    self%thermal_time=self%lattice_mover%thermal_time
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_lwf_mover_t),intent(inout) :: self
    nullify(self%lattice_mover)
    nullify(self%lwf_mover)
  end subroutine finalize

  subroutine run_time(self, effpot, displacement, strain, spin, lwf, energy_table)
    class(lattice_lwf_mover_t),intent(inout) :: self
    ! array of effective potentials so that there can be multiple of them.
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    integer :: i, nstep
    character(len=90) :: msg
    if(present(displacement) .or. present(strain)) then
       MSG_ERROR("displacement and strain should not be input for lattice mover")
    end if
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(energy_table)


    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')
    write(msg, '(A22)') "Lattice dynamic steps:"
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')
    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')

    write(msg, "(A13, 4X, A15, 4X, A15, 4X, A15, 4X, A15)") &
            &  "Iteration", "temperature(K)", "Ekin(Ha/uc)", &
            & "Epot(Ha/uc)", "ETOT(Ha/uc)"
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')
    
    nstep=floor(self%lattice_mover%thermal_time/self%lattice_mover%dt)
    do i =1, nstep
       call self%lattice_mover%run_one_step(effpot=effpot, displacement=self%lattice_mover%displacement, &
            & strain=strain, spin=spin, lwf=self%lwf_mover%lwf, energy_table=energy_table)
       !call self%lwf_mover%run_one_step(effpot=effpot, displacement=self%lattice_mover%displacement, strain=strain, spin=spin, &
       !     & lwf=lwf, energy_table=energy_table)
       call self%lwf_mover%hist%set_hist(lwf=self%lwf_mover%lwf, energy=self%lattice_mover%energy )
    end do

    nstep=floor(self%lattice_mover%total_time/self%lattice_mover%dt)
    do i =1, nstep
       call self%lattice_mover%run_one_step(effpot=effpot, displacement=self%lattice_mover%displacement, &
                                & lwf=self%lwf_mover%lwf, energy_table=energy_table)
       call self%lwf_mover%hist%set_hist(lwf=self%lwf_mover%lwf, energy=self%lattice_mover%energy )

       if(modulo(i, self%lattice_mover%params%nctime)==0) then
          write(msg, "(I13, 4X, F15.5, 4X, ES15.5, 4X, ES15.5, 4X, ES15.5)")  i, self%lattice_mover%T_ob*Ha_K, &
               & self%lattice_mover%Ek/self%lattice_mover%supercell%ncell, &
               & self%lattice_mover%energy/self%lattice_mover%supercell%ncell, &
               & (self%lattice_mover%Ek+self%lattice_mover%energy)/self%lattice_mover%supercell%ncell
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')

          self%lattice_mover%current_xcart(:, :) = self%lattice_mover%supercell%lattice%xcart(:,:)+self%lattice_mover%displacement
          call self%lattice_mover%ncfile%write_one_step(self%lattice_mover%current_xcart, &
               & self%lattice_mover%current_vcart, self%lattice_mover%energy, self%lattice_mover%Ek)

          call self%lwf_mover%ncfile%write_one_step(self%lwf_mover%hist)

       end if
    end do

    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')
  end subroutine run_time

end module m_lattice_lwf_mover

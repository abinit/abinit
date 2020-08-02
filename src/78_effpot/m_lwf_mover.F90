!!****m* ABINIT/m_lwf_mover
!! NAME
!! m_lwf_mover
!!
!! FUNCTION
!! This module contains the lwf mover, which controls how the lattice wannier function move.
!!
!!
!! Datatypes:
!!
!! * lwf_mover_t
!!
!! Subroutines:
!!
!! * lwf_mover_t_initialize
!! * lwf_mover_t_run_one_step
!! * lwf_mover_t_run_time
!! * TODO: update this when F2003 documentation format decided.
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
module m_lwf_mover
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_mpi_scheduler, only: mpi_scheduler_t, init_mpi_info
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_hashtable_strval, only: hash_table_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_lwf_hist, only: lwf_hist_t
  use m_lwf_observables, only: lwf_observables_t
  use m_lwf_ncfile, only: lwf_ncfile_t


  implicit none
  private
  !!***

  type, public, extends(abstract_mover_t) :: lwf_mover_t
     type(multibinit_dtset_type), pointer :: params
     real(dp) :: lwf_temperature, energy
     integer :: nlwf
     real(dp), allocatable :: lwf(:), lwf_force(:)
     type(lwf_ncfile_t) :: ncfile
     type(lwf_hist_t) :: hist
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     procedure :: set_initial_state
     procedure :: run_one_step
     procedure :: run_time
     procedure :: prepare_ncfile
     procedure :: set_ncfile_name
  end type lwf_mover_t

contains

  subroutine initialize(self, params, supercell, rng)
    class(lwf_mover_t), intent(inout) :: self
    type(multibinit_dtset_type),target, intent(in) :: params
    type(mbsupercell_t),target, intent(in) :: supercell
    type(rng_t), target, intent(in) :: rng
    self%params=>params
    self%supercell=>supercell
    self%label="LWF Mover"
    call self%set_params(params)
    call self%set_rng(rng)
    self%nlwf=self%supercell%lwf%nlwf
    ABI_ALLOCATE(self%lwf, (self%nlwf))
    self%lwf(:) = 0.0_dp
    ABI_ALLOCATE(self%lwf_force, (self%nlwf))
    self%lwf_force(:) = 0.0_dp
    self%energy=0.0_dp
    call self%hist%initialize(nlwf=self%nlwf, mxhist=1)
  end subroutine initialize


  subroutine finalize(self)
    class(lwf_mover_t), intent(inout) :: self
    nullify(self%supercell)
    nullify(self%params)
    ABI_SFREE(self%lwf)
    ABI_SFREE(self%lwf_force)
    call self%hist%finalize()
    !call self%ncfile%finalize()
  end subroutine finalize

  subroutine set_params(self, params)
    class(lwf_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    self%dt=params%lwf_dt
    self%total_time=params%lwf_ntime*params%lwf_dt
    self%lwf_temperature=params%lwf_temperature
  end subroutine set_params

  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf,  energy_table)
    ! run one step. (For MC also?)
    class(lwf_mover_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    class(abstract_potential_t), intent(inout) :: effpot
    type(hash_table_t),optional, intent(inout) :: energy_table
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(energy_table)
    MSG_ERROR("run_one_step not implemented for this mover")
  end subroutine run_one_step

  !-------------------------------------------------------------------!
  ! run from begining to end.
  !-------------------------------------------------------------------!
  subroutine run_time(self, effpot, displacement, strain, spin, lwf, energy_table)
    ! run one step. (For MC also?)
    class(lwf_mover_t), intent(inout) :: self
    ! array of effective potentials so that there can be multiple of them.
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    integer :: i, nstep
    character(len=90) :: msg
    if(present(lwf)) then
       MSG_ERROR("lwf should not be input for lwf mover")
    end if
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
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

    !nstep=floor(self%total_time/self%dt)
    !do i =1, nstep
    !   call self%run_one_step(effpot=effpot, spin=spin, lwf=lwf, energy_table=energy_table)
    !end do

    nstep=floor(self%total_time/self%dt)
    do i =1, nstep
       !print *, "Step: ", i,  "    T: ", self%T_ob*Ha_K, "    Ek:", self%Ek, "Ev", self%energy, "Etot", self%energy+self%Ek
       call self%run_one_step(effpot=effpot, spin=spin, lwf=self%lwf, energy_table=energy_table)

       call self%hist%set_hist(lwf=self%lwf, energy=self%energy )
       if(modulo(i, self%params%lwf_nctime)==0) then
          call self%ncfile%write_one_step(self%hist)
       end if
       !print *, "Step: ", i,   "Ev", self%energy, "Etot"

!       write(msg, "(I13, 4X, F15.5, 4X, ES15.5, 4X, ES15.5, 4X, ES15.5)")  i, self%T_ob*Ha_K, &
!            & self%Ek/self%supercell%ncell, self%energy/self%supercell%ncell, &
!            & (self%Ek+self%energy)/self%supercell%ncell
!       call wrtout(std_out,msg,'COLL')
!       call wrtout(ab_out, msg, 'COLL')
       write(msg, "(I13, 4X,  ES15.5)")  i, self%energy
       !            & (self%Ek+self%energy)/self%supercell%ncell
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       !TODO: output, observables
    end do

    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')

  end subroutine run_time


    !----------------------------------------------------------------------
    !> @brief set initial state.
    !>
    !> @param[in]  mode: a integer to define the kind of initial state.
    !----------------------------------------------------------------------
    subroutine set_initial_state(self, mode, restart_hist_fname)
      ! set initial positions, spin, etc
      class(lwf_mover_t), intent(inout) :: self
      integer, optional, intent(in) :: mode
      character(len=*), optional, intent(in) :: restart_hist_fname

      select case(mode)
      case(0)
         self%lwf(:)=0.0
      end select

      call self%hist%set_hist(lwf=self%lwf, energy=0.0_dp)

      ABI_UNUSED(mode)
      ABI_UNUSED(restart_hist_fname)
    end subroutine set_initial_state

    subroutine prepare_ncfile(self, params, fname)
      class(lwf_mover_t), intent(inout) :: self
      type(multibinit_dtset_type) :: params
      character(len=*), intent(in) :: fname
      integer :: master, my_rank, comm, nproc
      logical :: iam_master
      ABI_UNUSED_A(params)
      call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
      if(iam_master) then
         call self%ncfile%initialize( trim(fname), 1)
         call self%ncfile%write_cell(self%supercell)
         call self%ncfile%def_lwf_var(self%hist)
      end if
    end subroutine prepare_ncfile

    !-------------------------------------------------------------------!
    !set_ncfile_name :
    !-------------------------------------------------------------------!
    subroutine set_ncfile_name(self, params, fname)
      class(lwf_mover_t), intent(inout) :: self
      type(multibinit_dtset_type) :: params
      character(len=fnlen), intent(in) :: fname
      integer :: master, my_rank, comm, nproc
      logical :: iam_master
      call init_mpi_info(master, iam_master, my_rank, comm, nproc)
      if (iam_master) then
         call self%prepare_ncfile(params, trim(fname)//'_lwfhist.nc')
         call self%ncfile%write_one_step(self%hist)
      endif
    end subroutine set_ncfile_name



end module m_lwf_mover


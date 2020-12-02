!!****m* ABINIT/m_lattice_mover
!! NAME
!! m_lattice_mover
!!
!! FUNCTION
!! This module contains the lattice mover.
!!
!!
!! Datatypes:
!!
!! * lattice_mover_t: defines the lattice movers
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributorsi see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE



#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lattice_mover
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_random_xoroshiro128plus, only:  rng_t
  use m_hashtable_strval, only: hash_table_t
  use m_lattice_ncfile, only: lattice_ncfile_t
  use m_mpi_scheduler, only: init_mpi_info
!!***

  implicit none
  private

  type ,public, extends(abstract_mover_t) :: lattice_mover_t
     !> This is the abstract lattice mover

     type(multibinit_dtset_type), pointer :: params=>null() ! input parameters
     integer :: natom=0     ! number of atoms
     real(dp) :: stress(3,3), strain(3,3)  ! stress and strain
     real(dp), allocatable :: masses(:)  ! masses
     integer :: latt_dynamics=0 ! type of lattice dynamics

     ! hexu: is xcart needed?
     real(dp), allocatable :: current_xcart(:,:) ! xcart of current step
     real(dp), allocatable :: current_vcart(:,:) ! vcart of current step
     real(dp), allocatable :: forces(:,:)        ! forces
     real(dp), allocatable :: displacement(:,:)  ! displacement
     real(dp) :: energy  ! total energy
     real(dp) :: Ek      ! kinetic energy
     real(dp) :: T_ob    ! observed temperature
     logical :: is_null = .True.
     !> TODO: hist
     !type(lattice_hist_t) :: hist

     real(dp) :: mass_total 
     type(lattice_ncfile_t) :: ncfile
   contains
     procedure:: initialize       ! perhaps each effpot type should have own 
     procedure :: finalize
     procedure :: set_params
     procedure :: prepare_ncfile
     procedure :: set_ncfile_name
     procedure :: set_initial_state ! initial state
     procedure :: set_temperature
     procedure :: force_stationary
     procedure :: run_one_step
     procedure :: run_time
     procedure :: run_varT
     procedure :: reset            ! reset the mover
     procedure :: get_T_and_Ek     ! calculate temperature and kinetic energy.
     procedure :: calc_observables ! call functions to calculate observables
     procedure :: write_hist       ! write hist file
  end type lattice_mover_t


contains

  !-------------------------------------------------------------------!
  ! Initialize:
  !
  ! Inputs:
  !> params: input parameters
  !> supercell: supercell.
  !> rng: random number generator
  !-------------------------------------------------------------------!
  subroutine initialize(self, params, supercell, rng)
    class(lattice_mover_t), intent(inout) :: self
    type(multibinit_dtset_type),target, intent(in) :: params
    type(mbsupercell_t),target, intent(in) :: supercell
    type(rng_t), target, intent(in) :: rng
    self%params=>params
    self%supercell=>supercell
    self%label="Lattice Mover"
    self%natom = supercell%lattice%natom
    ABI_ALLOCATE(self%masses, (self%natom))
    ABI_ALLOCATE(self%displacement, (3, self%natom))
    ABI_ALLOCATE(self%current_xcart, (3, self%natom))
    ABI_ALLOCATE(self%current_vcart, (3, self%natom))
    ABI_ALLOCATE(self%forces, (3,self%natom))
    self%is_null=.False.
    self%strain(:,:) = 0.0
    self%stress(:,:) = 0.0
    self%forces(:,:) = 0.0
    self%displacement(:,:) = 0.0
    self%current_vcart(:,:) = 0.0
    call self%set_params(params)
    call self%set_rng(rng)
  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Finalize:
  !-------------------------------------------------------------------!
  subroutine finalize(self)
    class(lattice_mover_t), intent(inout) :: self
    nullify(self%supercell)
    nullify(self%params)
    self%label="Destroyed lattice mover"
    if (.not.self%is_null) then
       ABI_DEALLOCATE(self%masses)
       ABI_DEALLOCATE(self%current_xcart)
       ABI_DEALLOCATE(self%current_vcart)
       ABI_DEALLOCATE(self%forces)
       ABI_DEALLOCATE(self%displacement)
    endif
    self%is_null=.True.
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Set the mover using the input parameters
  !-------------------------------------------------------------------!
  subroutine set_params(self, params)
    ! set parameters from input file. (something else, like temperature for MvT calculation?)
    class(lattice_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    self%temperature = params%temperature !TODO: to Hartree ??
    self%dt =params%dtion 
    self%masses(:)=self%supercell%lattice%masses(:)
    self%mass_total = sum(self%masses)
    self%total_time = self%dt * params%ntime 
    self%latt_dynamics = params%dynamics
  end subroutine set_params

  !-------------------------------------------------------------------!
  ! Set the mover temperature
  !-------------------------------------------------------------------!
  subroutine set_temperature(self, temperature)
    class(lattice_mover_t), intent(inout) :: self
    real(dp),intent(in) :: temperature
    self%temperature = temperature !TODO: to Hartree ??
  end subroutine set_temperature


  !-------------------------------------------------------------------!
  ! set initial state:
  !  Inputs:
  !   mode: integer
  !    if mode=1, use a Boltzman distribution to init the velocities.
  !    if mode=2, ...
  !-------------------------------------------------------------------!
  subroutine set_initial_state(self, mode, restart_hist_fname)
    ! set initial positions, spin, etc
    class(lattice_mover_t), intent(inout) :: self
    integer, optional, intent(in) :: mode
    character(len=*), optional, intent(in) :: restart_hist_fname


    real(dp) :: xi(3, self%natom)
    integer :: i

    ABI_UNUSED(restart_hist_fname)

    if(mode==1) then ! using a boltzmann distribution. 
       ! Should only be used for a constant Temperature mover
       ! which includes:
       !   102:    Langevin
       !   103:    Brendesen
       if (.not.( &
          self%latt_dynamics==101 .or.  &  ! TODO remove
          self%latt_dynamics==102 .or. self%latt_dynamics==103 ) ) then
          MSG_ERROR("Only set lattice initial state with a Boltzmann distribution in a constant T mover.")
       end if
       call self%rng%rand_normal_array(xi, 3*self%natom)
       do i=1, self%natom
          self%current_vcart(:,i) = xi(:, i) *sqrt(self%temperature/self%masses(i))
       end do
       call self%force_stationary()
       call self%get_T_and_Ek()
    else if(mode==2) then ! Use reference structure and 0 velocity.
       ! other modes.
       if(self%latt_dynamics==102 .or. self%latt_dynamics==103 ) then
           MSG_ERROR("Displacement and velocity set to zero in a NVT mover.")
       end if
       do i=1, self%natom
          self%current_vcart(:,i) = 0.0
       end do
       self%current_xcart(:, :) = self%supercell%lattice%xcart(:,:)
       call self%get_T_and_Ek()
    end if

    ABI_UNUSED(restart_hist_fname)

  end subroutine set_initial_state

  subroutine prepare_ncfile(self, params, fname)
    class(lattice_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    character(len=*), intent(in) :: fname
    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    ABI_UNUSED_A(params)
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    if(iam_master) then
       call self%ncfile%initialize( trim(fname), 1)
       call self%ncfile%write_cell(self%supercell)
       call self%ncfile%def_lattice_var()
    end if
  end subroutine prepare_ncfile


  !-------------------------------------------------------------------!
  !set_ncfile_name :
  !-------------------------------------------------------------------!
  subroutine set_ncfile_name(self, params, fname)
    class(lattice_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    character(len=fnlen), intent(in) :: fname
    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc)
    if (iam_master) then
       call self%prepare_ncfile(params, trim(fname)//'_latthist.nc')
       call self%ncfile%write_one_step(self%current_xcart, self%current_vcart, self%energy, self%Ek )
    endif
  end subroutine set_ncfile_name



  !-------------------------------------------------------------------!
  ! Make sure the mass center does not move.
  !-------------------------------------------------------------------!
  subroutine force_stationary(self)
    class(lattice_mover_t), intent(inout) :: self
    integer :: i
    real(dp) :: p(3), pavg(3)
    p(:)=0.0
    do i = 1, self%natom
       p(:)=p(:)+self%current_vcart(:,i) * self%masses(i)
    end do
    pavg=p/self%mass_total
    do i = 1, self%natom
       self%current_vcart(:, i) = self%current_vcart(:, i) - pavg(:)
    end do
  end subroutine force_stationary


  !-------------------------------------------------------------------!
  ! Force the temperature strictly.
  ! Since the boltzman distribution has some fluctuation
  !-------------------------------------------------------------------!
  subroutine force_temperature(self)
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine force_temperature


  !-------------------------------------------------------------------!
  ! run_one_step
  !  run one step of dynamics.
  !  Should be overrided.
  !  Inputs:
  !> effpot: effective potential
  !> displacement: should NOT be provided, since it is already stored.
  !> strain: Also should NOT be provided.
  !> spin: should be provided only if there is spin-lattice coupling
  !> lwf : should be provided only if there is lattice-lwf coupling (unlikely)
  !> energy_table: energy_table.
  !-------------------------------------------------------------------!
  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf, energy_table)
    ! run one step. (For MC also?)
    class(lattice_mover_t),      intent(inout) :: self    ! array of effective potentials so that there can be multiple of them.
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional,          intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table

    character(len=40) :: key

    if(present(displacement) .or. present(strain)) then
       MSG_ERROR("displacement and strain should not be input for lattice mover")
    end if

    MSG_BUG("The abstract lattice mover is used, which should be a bug.")

    ABI_UNUSED_A(self)
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(energy_table)

    call self%get_T_and_Ek()
    if (present(energy_table)) then
      key = 'Lattice kinetic energy'
      call energy_table%put(key, self%Ek)
    end if
  end subroutine run_one_step


  !-------------------------------------------------------------------!
  !get_temperature_and_kinetic_energy
  ! Ek = 1/2 \sum m_i vi^2
  ! T = 2/3 Ek/natom  (in a.u.)
  !-------------------------------------------------------------------!
  subroutine get_T_and_Ek(self)
    class(lattice_mover_t), intent(inout) :: self
    integer :: i
    self%Ek=0.0
    do i =1, self%natom
       self%Ek = self%Ek+ 0.5* self%masses(i) *  &
            &sum(self%current_vcart(:,i)*self%current_vcart(:,i))
    end do
    ! temperature
    self%T_ob = 2.0*self%Ek/(3*self%natom)
  end subroutine get_T_and_Ek


  !-------------------------------------------------------------------!
  ! run from begining to end.
  !-------------------------------------------------------------------!
  subroutine run_time(self, effpot, displacement, strain, spin, lwf, energy_table)
    class(lattice_mover_t), intent(inout) :: self
    ! array of effective potentials so that there can be multiple of them.
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    integer :: i, nstep
    character(len=90) :: msg
    if(present(displacement) .or. present(strain)) then
       MSG_ERROR("displacement and strain should not be input for lattice mover")
    end if
    ABI_UNUSED_A(self)
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
    
    nstep=floor(self%thermal_time/self%dt)
    do i =1, nstep
       call self%run_one_step(effpot=effpot, spin=spin, lwf=lwf, energy_table=energy_table)
    end do

    nstep=floor(self%total_time/self%dt)
    do i =1, nstep
       !print *, "Step: ", i,  "    T: ", self%T_ob*Ha_K, "    Ek:", self%Ek, "Ev", self%energy, "Etot", self%energy+self%Ek
       call self%run_one_step(effpot=effpot, spin=spin, lwf=lwf, energy_table=energy_table)
       if(modulo(i, self%params%nctime)==0) then
          write(msg, "(I13, 4X, F15.5, 4X, ES15.5, 4X, ES15.5, 4X, ES15.5)")  i, self%T_ob*Ha_K, &
               & self%Ek/self%supercell%ncell, self%energy/self%supercell%ncell, &
               & (self%Ek+self%energy)/self%supercell%ncell
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')

          self%current_xcart(:, :) = self%supercell%lattice%xcart(:,:)+self%displacement
          call self%ncfile%write_one_step(self%current_xcart, self%current_vcart, self%energy, self%Ek)
       end if
       !TODO: output, observables
    end do

    msg=repeat("=", 90)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')

  end subroutine run_time


  !-------------------------------------------------------------------!
  ! Reset:
  ! It reset the counter of steps, but does not set the initial state
  ! again.
  !-------------------------------------------------------------------!
  subroutine reset(self)
    ! reset the state of mover (e.g. counter->0)
    ! so it can be reused.
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine reset

  !-------------------------------------------------------------------!
  !Calc_observables
  ! 
  !-------------------------------------------------------------------!
  subroutine calc_observables(self)
    ! call functions to calculate observables.
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
  end subroutine calc_observables

  !-------------------------------------------------------------------!
  ! Write_hist: write to hist file
  !-------------------------------------------------------------------!
  subroutine write_hist(self)
    ! write to hist file
    class(lattice_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)

  end subroutine write_hist

  !-------------------------------------------------------------------!
  ! Get_state: get the current state
  !-------------------------------------------------------------------!
  subroutine get_state(self, displacement, strain, spin, lwf, ihist)
    ! get the state of the ihist(th) step. ihist can be 0 (current), -1 (last), ... -maxhist..
    class(lattice_mover_t), intent(in):: self
    real(dp), optional, intent(inout) :: displacement, strain, spin, lwf
    integer, optional, intent(in):: ihist
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(ihist)
  end subroutine get_state


  !!****f* m_lwf_mover/run_varT
  !!
  !! NAME
  !! run_varT
  !!
  !! FUNCTION
  !! run M vs Temperature
  !!
  !! INPUTS
  !! pot: potential
  !! T_start, Tend, T_nstep
  !u
  !! OUTPUT
  !!
  !! PARENTS
!!
  !! CHILDREN
!!      self%hist%finalize,self%mps%finalize,self%spin_mc%finalize
!!      self%spin_ob%finalize
!!
  !! SOURCE
  subroutine  run_varT(self, pot, ncfile_prefix, displacement, strain, spin, lwf, energy_table)
    class(lattice_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: pot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:), spin(:,:)
    character(fnlen), intent(inout) :: ncfile_prefix
    type(hash_table_t), optional, intent(inout) :: energy_table
    real(dp) :: T_start, T_end
    integer :: T_nstep
    !type(lwf_ncfile_t) :: lwf_ncfile
    character(len=4) :: post_fname
    real(dp) :: T, T_step
    integer :: i
    !integer :: Tfile, iostat
    character(len=90) :: msg
    !character(len=4200) :: Tmsg ! to write to var T file
    !character(len=150) :: iomsg
    !character(fnlen) :: Tfname ! file name for output various T calculation
    !real(dp), allocatable :: Tlist(:), chi_list(:), Cv_list(:), binderU4_list(:)
    !real(dp), allocatable :: Mst_sub_norm_list(:, :)
    !real(dp), allocatable ::  Mst_norm_total_list(:)

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if (iam_master) then
       T_start=self%params%latt_temperature_start
       T_end=self%params%latt_temperature_end
       T_nstep=self%params%latt_temperature_nstep
       !Tfile=get_unit()
       !Tfname = trim(ncfile_prefix)//'.varT'
       !iostat=open_file(file=Tfname, unit=Tfile, iomsg=iomsg )
       if (T_nstep<=1) then
          T_step=0.0
       else
          T_step=(T_end-T_start)/(T_nstep-1)
       endif
       write(msg, "(A52, ES13.5, A11, ES13.5, A1)") & 
            & "Starting temperature dependent calculations. T from ", &
            & T_start*Ha_K, "K to ", T_end*Ha_K, " K."
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")
    end if

    call xmpi_bcast(T_nstep, 0, comm, ierr)
    do i=1, T_nstep
       if(iam_master) then
          T=T_start+(i-1)*T_step
          msg=repeat("=", 79)
          call wrtout(std_out, msg, "COLL")
          call wrtout(ab_out, msg, "COLL")

          write(msg, "(A13, 5X, ES13.5, A3)") "Temperature: ", T*Ha_K, " K."
          call wrtout(std_out, msg, "COLL")
          call wrtout(ab_out,  msg, "COLL")

          ! set temperature
          ! TODO make this into a subroutine set_params
       endif
       call self%set_temperature(temperature=T)
       if(iam_master) then
          if(i==1) then
             call self%set_initial_state(mode=1)
          endif

          write(post_fname, "(I4.4)") i
          call self%prepare_ncfile( self%params, &
               & trim(ncfile_prefix)//'_T'//post_fname//'_latthist.nc')
          call self%ncfile%write_one_step(self%current_xcart, self%current_vcart, self%energy, self%Ek )
       endif

       call self%run_time(pot, spin=spin, &
            & lwf=lwf, energy_table=energy_table)

       if(iam_master) then
          call self%ncfile%finalize()
       endif
    end do

  end subroutine run_varT
  !!***



end module m_lattice_mover


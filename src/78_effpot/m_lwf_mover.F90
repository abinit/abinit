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
  use m_nctk
#define HAVE_NETCDF 1
#if defined HAVE_NETCDF
  use netcdf
#endif
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
     real(dp) ::  energy
     integer :: nlwf
     real(dp), allocatable :: lwf(:), lwf_force(:), vcart(:)
     type(lwf_ncfile_t) :: ncfile
     type(lwf_hist_t) :: hist
     real(dp), pointer :: lwf_masses(:) => null()
     real(dp) :: Ek=0.0_dp     ! kinetic energy
     real(dp) :: T_ob=0.0_dp    ! observed temperature
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_temperature
     procedure :: set_params
     procedure :: set_initial_state
     procedure :: get_T_and_Ek
     procedure :: run_one_step
     procedure :: run_time
     procedure :: run_varT
     procedure :: prepare_ncfile
     procedure :: set_ncfile_name
     procedure :: read_hist_lwf_state
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
    ABI_ALLOCATE(self%vcart, (self%nlwf))
    ABI_ALLOCATE(self%lwf_force, (self%nlwf))
    self%lwf(:) = 0.0_dp
    self%lwf_force(:) = 0.0_dp
    self%vcart(:) = 0.0_dp
    self%energy=0.0_dp
    self%lwf_masses=>self%supercell%lwf%lwf_masses
    call self%hist%initialize(nlwf=self%nlwf, mxhist=1)
  end subroutine initialize


  subroutine finalize(self)
    class(lwf_mover_t), intent(inout) :: self
    nullify(self%supercell)
    nullify(self%params)
    ABI_SFREE(self%lwf)
    ABI_SFREE(self%vcart)
    ABI_SFREE(self%lwf_force)
    nullify(self%lwf_masses)
    call self%hist%finalize()
    !call self%ncfile%finalize()
  end subroutine finalize

  subroutine set_params(self, params)
    class(lwf_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    self%dt=params%lwf_dt
    self%total_time=params%lwf_ntime*params%lwf_dt
    self%temperature=params%lwf_temperature
  end subroutine set_params

  subroutine set_temperature(self, temperature)
    class(lwf_mover_t), intent(inout) :: self
    real(dp), intent(in) :: temperature
    self%temperature=temperature
  end subroutine set_temperature

  !-------------------------------------------------------------------!
  !get_temperature_and_kinetic_energy
  ! Ek = 1/2 \sum m_i vi^2
  ! T = 2 Ek/nlwf (in a.u.)
  !-------------------------------------------------------------------!
  subroutine get_T_and_Ek(self)
    class(lwf_mover_t), intent(inout) :: self
    integer :: i
    self%Ek= sum(self%lwf_masses * (self%vcart * self%vcart))
    self%T_ob = 2.0*self%Ek/self%nlwf
  end subroutine get_T_and_Ek



  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf,  energy_table)
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
    write(msg, '(A22)') "LWF dynamic steps:"
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

       call self%hist%set_hist(lwf=self%lwf, vcart=self%vcart, energy=self%energy )
       if(modulo(i, self%params%lwf_nctime)==0) then
          call self%ncfile%write_one_step(self%hist)
       !print *, "Step: ", i,   "Ev", self%energy, "Etot"

       write(msg, "(I13, 4X, F15.5, 4X, ES15.5, 4X, ES15.5, 4X, ES15.5)")  i, self%T_ob*Ha_K, &
            & self%Ek/self%supercell%ncell, self%energy/self%supercell%ncell, &
            & (self%Ek+self%energy)/self%supercell%ncell
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')
       !write(msg, "(I13, 4X,  ES15.5)")  i, self%energy/self%supercell%ncell
       !            & (self%Ek+self%energy)/self%supercell%ncell

       end if
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
      integer :: i
      real(dp) :: tmp
      real(dp) :: kpoint(3)

      self%lwf(:)=0.0
      select case(mode)
      case(0)
         kpoint(:)=[0.5_dp, 0.0_dp, 0.5_dp]
         do i=1, self%supercell%ncell
           tmp=0.2*real(exp(cmplx(0.0,two_pi, kind=dp) * &
                               &dot_product(kpoint, self%supercell%supercell_maker%rvecs(:, i))), kind=dp)
           self%lwf(i*2-1)=tmp
           self%lwf(i*2)=tmp
         enddo
      ! random
      case(1)
         call self%rng%rand_unif_01_array(self%lwf, self%nlwf)
         self%lwf=(self%lwf-0.5)*0.1
      ! zero
      case(2)
         self%lwf(:)=0.0
      ! read from lwf hist file
      case(4)
         !print*, "Reading from lwf hist file: ", trim(restart_hist_fname)
         call self%read_hist_lwf_state(restart_hist_fname)
      end select

      call self%rng%rand_normal_array(self%vcart(:), self%nlwf)
      do i=1, self%nlwf
         self%vcart(i) = self%vcart(i) *sqrt(self%temperature/self%lwf_masses(i))
      end do

      call self%hist%set_hist(lwf=self%lwf, vcart=self%vcart, energy=0.0_dp)

    end subroutine set_initial_state

  !-------------------------------------------------------------------!
  ! read_hist_lwf_state
  !  read the last step of spin from hist file.
  !-------------------------------------------------------------------!
  subroutine read_hist_lwf_state(self, fname)
    class(lwf_mover_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: fname
    integer :: ierr, ncid, varid
    integer :: nlwf, ntime
    character(len=118) :: msg
    ! open file

#if defined HAVE_NETCDF
    ierr=nf90_open(trim(fname), NF90_NOWRITE, ncid)
    NCF_CHECK_MSG(ierr, "The lwf_init_state is set to 4. But opening netcdf file "//trim(fname)//" Failed. ")

    ! sanity check. If the hist file is consistent with the current calculation
    ierr=nctk_get_dim(ncid, "nlwf" , nlwf)
    NCF_CHECK_MSG(ierr, "when reading nlwf")

    msg="The number of lwfs in histfile is not equal & & to the present calculation." // &
         & " Please check if the file is consistent."
    if (nlwf/= self%nlwf) then
       MSG_ERROR(msg)
    end if

    ierr=nctk_get_dim(ncid, "ntime", ntime)
    NCF_CHECK_MSG(ierr, "when reading ntime")

    ! read lwf and set as initial state
    ierr =nf90_inq_varid(ncid, "lwf", varid)
    NCF_CHECK_MSG(ierr, "when reading lwf.")

    ierr = nf90_get_var(ncid=ncid, varid=varid, values=self%lwf(:), &
         & start=(/ 1, ntime/), count=(/nlwf,1/))
    NCF_CHECK_MSG(ierr, "when reading lwf from lwf hist file")

    ! close file
    ierr=nf90_close(ncid)
    NCF_CHECK_MSG(ierr, "Close netcdf file")
#else
    MSG_ERROR("lwf_init_state set to 4 but abinit is not compiled with netcdf.")
#endif

  end subroutine read_hist_lwf_state



    subroutine prepare_ncfile(self, params, fname)
      class(lwf_mover_t), intent(inout) :: self
      type(multibinit_dtset_type) :: params
      character(len=*), intent(in) :: fname
      integer :: master, my_rank, comm, nproc
      logical :: iam_master
      ABI_UNUSED_A(params)
      call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
      if(iam_master) then
         call self%ncfile%initialize(fname, 1)
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
    class(lwf_mover_t), intent(inout) :: self
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
       T_start=self%params%lwf_temperature_start
       T_end=self%params%lwf_temperature_end
       T_nstep=self%params%lwf_temperature_nstep
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

          call self%hist%reset(array_to_zero=.False.)
          ! set temperature
          ! TODO make this into a subroutine set_params
       endif
       call self%set_temperature(temperature=T)
       if(iam_master) then
          if(i==1) then
             call self%set_initial_state(mode=self%params%lwf_init_state)
          endif

          write(post_fname, "(I4.4)") i
          call self%prepare_ncfile( self%params, &
               & trim(ncfile_prefix)//'_T'//post_fname//'_lwfhist.nc')
          call self%ncfile%write_one_step(self%hist)
       endif

       call self%run_time(pot, displacement=displacement, strain=strain, spin=spin, &
            & lwf=lwf, energy_table=energy_table)

       if(iam_master) then
          call self%ncfile%finalize()
       endif
    end do

  end subroutine run_varT
  !!***


end module m_lwf_mover


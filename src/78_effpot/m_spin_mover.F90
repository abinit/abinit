!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_mover
!! NAME
!! m_spin_mover
!!
!! FUNCTION
!! This module contains the spin mover, which controls how the spin 
!!
!!
!! Datatypes:
!!
!! * spin_mover_t
!!
!! Subroutines:
!!
!! * spin_mover_t_initialize
!! * spin_mover_t_run_one_step
!! * spin_mover_t_run_time
!! * TODO: update this when F2003 documentation format decided.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
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

module m_spin_mover

  use defs_basis
  use m_nctk
#if defined HAVE_NETCDF
  use netcdf
#endif
  use m_errors
  use m_abicore
  use m_xmpi
  use m_io_tools, only : get_unit, open_file, close_unit
  use m_mpi_scheduler, only: mpi_scheduler_t, init_mpi_info
  use m_mathfuncs, only : cross
  use m_spin_observables , only : spin_observable_t
  use m_spin_potential, only:  spin_potential_t
  use m_spin_hist, only: spin_hist_t
  use m_spin_ncfile, only: spin_ncfile_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_spin_mc_mover, only : spin_mc_t
  use m_hashtable_strval, only: hash_table_t
  implicit none
  private
  !!***


  !!****t* m_spin_mover/spin_mover_t
  !! NAME
  !! spin_mover_t
  !!
  !! FUNCTION
  !! this type contains the parameters for the spin mover.
  !!
  !! It contains:
  !! dt: time step
  !! total_time
  !! temperature.
  !! nspin number of magnetic atoms
  !! SOURCE


  type, public, extends(abstract_mover_t) :: spin_mover_t
     integer :: nspin, method
     real(dp), allocatable :: gyro_ratio(:), damping(:), gamma_L(:), H_lang_coeff(:), ms(:), Stmp(:,:), Stmp2(:,:)
     real(dp), allocatable :: Heff_tmp(:,:), Htmp(:,:), Hrotate(:,:), H_lang(:,:), buffer(:,:)
     real(dp) :: init_qpoint(3), init_rotate_axis(3) ! qpoint and rotation axis to set up initial spin configuration
     real(dp) :: init_orientation(3) ! spin orientation in primitive cell which is then rotated
     type(spin_hist_t) :: hist
     logical :: gamma_l_calculated
     type(spin_mc_t) :: spin_mc
     type(mpi_scheduler_t) :: mps
     type(spin_observable_t) :: spin_ob
     type(spin_ncfile_t) :: spin_ncfile
     type(multibinit_dtset_type), pointer :: params
   CONTAINS
     procedure :: initialize
     procedure :: finalize
     procedure :: set_initial_state
     procedure :: read_hist_spin_state
     procedure, private :: run_one_step_DM => spin_mover_t_run_one_step_DM
     procedure, private :: run_one_step_HeunP => spin_mover_t_run_one_step_HeunP
     procedure, private :: run_one_step_dummy=> spin_mover_t_run_one_step_dummy
     procedure, private :: run_one_step_MC=> spin_mover_t_run_one_step_MC
     procedure :: run_one_step => spin_mover_t_run_one_step
     procedure :: run_time => spin_mover_t_run_time
     procedure :: run_MvT
     procedure :: set_temperature
     procedure, private :: prepare_ncfile
     procedure, private ::get_Langevin_Heff
     procedure :: current_spin
     procedure :: set_ncfile_name
  end type spin_mover_t
  !!***

contains

  !!****f* m_spin_mover/initialize
  !!
  !! NAME
  !!  initialize
  !!
  !! FUNCTION
  !!  initialize the spin mover
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine initialize(self, params, supercell, rng, restart_hist_fname)
    class(spin_mover_t), intent(inout) :: self
    type(multibinit_dtset_type), target :: params
    type(mbsupercell_t), target :: supercell
    type(rng_t), target, intent(in) :: rng
    character(len=fnlen), optional, intent(in) :: restart_hist_fname
    integer ::  nspin

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    self%params=>params
    self%supercell=>supercell
    self%rng => rng
    if (iam_master) then
       nspin=supercell%spin%nspin
       self%nspin=nspin
       self%dt= params%spin_dt
       self%thermal_time= params%spin_ntime_pre * self%dt
       self%total_time= params%spin_ntime * self%dt
       self%temperature=params%spin_temperature
       if(params%spin_dynamics>=0) then
          self%method=params%spin_dynamics
       endif
       if(params%spin_init_state==3) then
         self%init_qpoint = params%spin_init_qpoint
         self%init_rotate_axis = params%spin_init_rotate_axis
         self%init_orientation = params%spin_init_orientation
       endif
    end if
    if(params%spin_dynamics==3) then ! Monte carlo
       call self%spin_mc%initialize(nspin=nspin, angle=1.0_dp, temperature=params%spin_temperature)
    end if
    call xmpi_bcast(self%nspin, master, comm, ierr)
    call xmpi_bcast(self%dt, master, comm, ierr)
    call xmpi_bcast(self%thermal_time, master, comm, ierr)
    call xmpi_bcast(self%total_time, master, comm, ierr)
    call xmpi_bcast(self%temperature, master, comm, ierr)
    call xmpi_bcast(self%method, master, comm, ierr)

    ABI_ALLOCATE(self%ms, (self%nspin) )
    ABI_ALLOCATE(self%gyro_ratio, (self%nspin) )
    ABI_ALLOCATE(self%damping, (self%nspin) )
    ABI_ALLOCATE(self%gamma_l, (self%nspin) )
    ABI_ALLOCATE(self%H_lang_coeff, (self%nspin) )

    ABI_ALLOCATE(self%Heff_tmp, (3,self%nspin) )
    ABI_ALLOCATE(self%Htmp, (3,self%nspin) )
    ABI_ALLOCATE(self%Hrotate, (3,self%nspin) )
    ABI_ALLOCATE(self%Stmp, (3,self%nspin) )
    ABI_ALLOCATE(self%Stmp2, (3,self%nspin) )
    ABI_ALLOCATE(self%buffer, (3,self%nspin) )
    ABI_ALLOCATE(self%H_lang, (3,self%nspin) )

    self%gamma_l_calculated=.False.
    call self%mps%initialize(ntasks=nspin,master=master, comm=comm)


    call xmpi_bcast(params%spin_damping, master, comm, ierr)

    if (iam_master) then
       if (params%spin_damping >=0) then
          self%damping(:)= params%spin_damping
       else
          self%damping(:)=supercell%spin%gilbert_damping(:)
       end if

       self%gyro_ratio(:)=supercell%spin%gyro_ratio(:)
       self%ms(:)=supercell%spin%ms(:)
    endif

    call xmpi_bcast(self%damping, master, comm, ierr)
    call xmpi_bcast(self%gyro_ratio, master, comm, ierr)
    call xmpi_bcast(self%ms, master, comm, ierr)
    call self%set_temperature(temperature=params%spin_temperature)

    ! Hist and set initial spin state
    if(iam_master) then
       call self%hist%initialize(nspin=self%nspin, &
            &   mxhist=3, has_latt=.False.)
       call self%hist%set_params(spin_nctime=params%spin_nctime, &
            &     spin_temperature=params%spin_temperature)
    endif

    if(present(restart_hist_fname)) then 
      call self%set_initial_state(mode=params%spin_init_state, restart_hist_fname=restart_hist_fname)
    else
      call self%set_initial_state(mode=params%spin_init_state)
    endif
      
    ! observable
    if(iam_master) then
       call self%spin_ob%initialize(self%supercell, params)
    endif

  end subroutine initialize
  !!***

  !-------------------------------------------------------------------!
  !set_ncfile_name :
  !-------------------------------------------------------------------!
  subroutine set_ncfile_name(self, params, fname)
    class(spin_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    character(len=fnlen), intent(in) :: fname
    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc)
    if (iam_master) then
       call self%prepare_ncfile(params, trim(fname)//'_spinhist.nc')
       call self%spin_ncfile%write_one_step(self%hist)
    endif
  end subroutine set_ncfile_name


  !-------------------------------------------------------------------!
  ! read_hist_spin_state
  !  read the last step of spin from hist file.
  !  and save if to self%Stmp
  !-------------------------------------------------------------------!
  subroutine read_hist_spin_state(self, fname)
    class(spin_mover_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: fname
    integer :: ierr, ncid, varid
    integer :: nspin, ntime
    ! open file

#if defined HAVE_NETCDF
    ierr=nf90_open(trim(fname), NF90_NOWRITE, ncid)
    NCF_CHECK_MSG(ierr, "The spin_init_mode is set to 4. But opening netcdf file "//trim(fname)//" Failed. ")

    ! sanity check. If the hist file is consistent with the current calculation
    ierr=nctk_get_dim(ncid, "nspin" , nspin)
    NCF_CHECK_MSG(ierr, "when reading nspin")

    if (nspin /= self%nspin) then
       MSG_ERROR("The number of spins in histfile is not equal & 
           & to the present calculation. &
           & Please check if the file is consistent.")
    end if


    ierr=nctk_get_dim(ncid, "ntime", ntime)
    NCF_CHECK_MSG(ierr, "when reading ntime")


    ! TODO: more check ???

    ! read Spin and set as initial state
    ierr =nf90_inq_varid(ncid, "S", varid)
    NCF_CHECK_MSG(ierr, "when reading S. Try using spin_init_state=3 option instead (specify spin_init_qpoint,&
      &  spin_init_rotate_axis and spin_init_orientation as needed).")
    
    ierr = nf90_get_var(ncid=ncid, varid=varid, values=self%Stmp(:,:), start=(/1, 1, ntime/), count=(/3, nspin,1/))
    NCF_CHECK_MSG(ierr, "when reading S from spin hist file")

    ! close file
    ierr=nf90_close(ncid)
    NCF_CHECK_MSG(ierr, "Close netcdf file")
#else
    MSG_ERROR("spin_init_state set to 4 but abinit is not compiled with netcdf.")
#endif 

  end subroutine read_hist_spin_state

  !----------------------------------------------------------------------------!
  !set_initial_state:
  ! mode: which configuration to use
  !   1. Random
  !   2. reference state from potential file
  !   3. spin configuration using qpoint and rotation axis (e.g. for FM or AFM)
  !   4. Restart from last entry of hist netcdf file
  !----------------------------------------------------------------------------!
  subroutine set_initial_state(self, mode, restart_hist_fname)
    class(spin_mover_t),            intent(inout) :: self
    integer,              optional, intent(in)    :: mode
    character(len=*), optional, intent(in)    :: restart_hist_fname

    integer :: i, init_mode
    character(len=500) :: msg

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    real(dp), allocatable :: Sprim(:,:)

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
       if (present(mode)) then
          init_mode=mode
       else
          init_mode=1
       end if

       select case (init_mode)
         case (1)
           ! randomize S using uniform random number
           write(msg,*) "Initial spins set to random values."
           call wrtout(ab_out,msg,'COLL')
           call wrtout(std_out,msg,'COLL')
           call random_number(self%Stmp)
           self%Stmp=self%Stmp-0.5
           do i=1, self%nspin
             self%Stmp(:,i)=self%Stmp(:,i)/sqrt(sum(self%Stmp(:, i)**2))
           end do

         case (2)
           ! set spin to reference state using the reference qpoint and rotation axis from potential file
           write(msg,*) "Initial spins set to reference configuration."
           call wrtout(ab_out,msg,'COLL')
           call wrtout(std_out,msg,'COLL')

           do i=1, self%nspin
             self%Stmp(:,:) = self%supercell%spin%Sref(:,:)
           end do

         case (3)
           write(msg,*) "Initial spins set according to spin_init_* variables."
           call wrtout(ab_out,msg,'COLL')
           call wrtout(std_out,msg,'COLL')

           ABI_ALLOCATE(Sprim, (3,self%supercell%unitcell%spin%nspin) )

           ! set inital spin state using the input variables
           ! set spin to ferromagnetic along init_orientation then rotate
           do i=1, self%supercell%unitcell%spin%nspin
             Sprim(:,i)=self%init_orientation(:)
           enddo
           self%Stmp(:,:) = 0.0d0

           call self%supercell%supercell_maker%generate_spin_wave_vectorlist(A=Sprim, &
             & kpoint=self%init_qpoint, axis=self%init_rotate_axis, A_sc=self%Stmp)

           ABI_SFREE(SPrim)

         case (4)
          ! read from last step of hist file
          write(msg,'(a,a,a)') "Initial spins set to input spin hist file ",&
             &  trim(restart_hist_fname), '.'  
          call wrtout(ab_out,msg,'COLL')
          call wrtout(std_out,msg,'COLL')
          if (.not. present(restart_hist_fname)) then
             MSG_ERROR("Spin initialize mode set to 4, but restart_hist_fname is not used.")
           end if
           call self%read_hist_spin_state(fname=restart_hist_fname)

       end select

       call self%hist%set_vars(S=self%Stmp, Snorm=self%supercell%spin%ms, &
            &  time=0.0_dp, ihist_latt=0, inc=.True.)

    endif
    call xmpi_bcast(self%Stmp, 0, comm, ierr)
  end subroutine set_initial_state


  !-------------------------------------------------------------------!
  ! prepare_ncfile:
  !-------------------------------------------------------------------!
  subroutine prepare_ncfile(self, params, fname)
    class(spin_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    character(len=*), intent(in) :: fname

    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
       call self%spin_ncfile%initialize( trim(fname), params%spin_write_traj)
       call self%spin_ncfile%def_spindynamics_var(self%hist)
       call self%spin_ncfile%def_observable_var(self%spin_ob)
       call self%spin_ncfile%write_primitive_cell(self%supercell%unitcell)
       call self%spin_ncfile%write_supercell(self%supercell)
       call self%spin_ncfile%write_parameters(params)
    endif
  end subroutine prepare_ncfile



  subroutine set_temperature(self, temperature)
    class(spin_mover_t), intent(inout) :: self
    real(dp), optional, intent(in) ::  temperature

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(present(temperature)) self%temperature=temperature
    call xmpi_bcast(self%temperature, master, comm, ierr)
    if(self%method==3) then
       if(iam_master) self%spin_mc%temperature = temperature
       if(iam_master) self%spin_mc%beta=1.0_dp/temperature
       call xmpi_bcast(self%spin_mc%temperature, master, comm, ierr)
       call xmpi_bcast(self%spin_mc%beta, master, comm, ierr)
    end if
    self%gamma_l(:)= self%gyro_ratio(:)/(1.0_dp+ self%damping(:)**2)
    self%gamma_l_calculated=.True.
    self%H_lang_coeff(:)=sqrt(2.0*self%damping(:)* self%temperature &
         &  /(self%gyro_ratio(:)* self%dt *self%ms(:)))
  end subroutine set_temperature


  subroutine get_Langevin_Heff(self, H_lang)
    class(spin_mover_t), intent(inout) :: self
    real(dp), intent(inout):: H_lang(3,self%nspin)
    integer :: i
    if ( self%temperature .gt. 1d-7) then
       call rand_normal_array(self%rng, H_lang(:, self%mps%istart:self%mps%iend), 3*self%mps%ntask)
       do i = self%mps%istart, self%mps%iend
          H_lang(:,i)= H_lang(:,i) * self%H_lang_coeff(i)
       end do
    else
       H_lang(:,:)=0.0_dp
    end if
  end subroutine get_Langevin_Heff


  !!****f* m_spin_mover/spin_mover_t_run_one_step_HeunP
  !!
  !! NAME
  !!  spin_mover_t_run_one_step_HeunP
  !!
  !! FUNCTION
  !! run one spin step using HeunP method
  !!
  !! INPUTS
  !! effpot: abstract_potential_t type.
  !! S_in : input spin. (3*nspin)
  !!
  !! OUTPUT
  !! etot: energy (scalar)
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_mover_t_run_one_step_HeunP(self, effpot, S_in, &
       & etot, displacement, strain, lwf, energy_table)
    !class (spin_mover_t), intent(inout):: self
    class(spin_mover_t), intent(inout):: self
    class(abstract_potential_t), intent(inout) :: effpot

    real(dp), optional, intent(inout):: displacement(:,:), &
         strain(:,:), lwf(:)
    real(dp), intent(inout) :: S_in(3,self%nspin)
    real(dp), intent(out) ::  etot
    integer :: i
    real(dp) :: dSdt(3), Htmp(3), Ri(3)
    type(hash_table_t),optional, intent(inout) :: energy_table

    !integer :: master, my_rank, comm, nproc, ierr
    !logical :: iam_master
    !call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    ! predict
    etot=0.0
    self%Heff_tmp(:,:)=0.0
    call effpot%calculate(displacement=displacement, strain=strain, lwf=lwf, spin=S_in, &
         & bfield=self%Heff_tmp, energy=etot, energy_table=energy_table)
    call self%get_Langevin_Heff(self%H_lang)
    do i=self%mps%istart, self%mps%iend
       Htmp=self%Heff_tmp(:,i)+self%H_lang(:,i)
       Ri = cross(S_in(:,i),Htmp)
       dSdt = -self%gamma_L(i)*(Ri+self%damping(i)* cross(S_in(:,i), Ri))
       Ri=S_in(:,i)+dSdt*self%dt
       Ri=Ri/sqrt(Ri(1)*Ri(1)+Ri(2)*Ri(2)+Ri(3)*Ri(3))
       self%Stmp2(:,i)=Ri
    end do
    call self%mps%allgatherv_dp2d(self%Stmp2, 3, buffer=self%buffer)

    ! correction
    self%Htmp(:,:)=0.0
    etot=0.0
    call effpot%calculate(displacement=displacement, strain=strain, lwf=lwf,spin=self%Stmp2, &
         & bfield=self%Htmp, energy=etot, energy_table=energy_table)
    do i=self%mps%istart, self%mps%iend
       Htmp=(self%Heff_tmp(:,i)+self%Htmp(:,i))*0.5_dp+self%H_lang(:,i)
       Ri = cross(S_in(:,i),Htmp)
       dSdt = -self%gamma_L(i)*(Ri+self%damping(i)* cross(S_in(:,i), Ri))
       Ri=S_in(:,i)+dSdt*self%dt
       Ri=Ri/sqrt(Ri(1)*Ri(1)+Ri(2)*Ri(2)+Ri(3)*Ri(3))
       self%Stmp(:,i)=Ri
    end do
    call self%mps%allgatherv_dp2d(self%Stmp, 3, buffer=self%buffer)
  end subroutine spin_mover_t_run_one_step_HeunP
  !!***



  !!****f* m_spin_mover/spin_mover_t_run_one_step_dummy
  !!
  !! NAME
  !!  spin_mover_t_run_one_step_dummy
  !!
  !! FUNCTION
  !! run one spin step using dummy method
  !!
  !! INPUTS
  !! effpot: abstract_potential_t type.
  !! S_in : input spin. (3*nspin)
  !!
  !! OUTPUT
  !! etot: energy (scalar)
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_mover_t_run_one_step_dummy(self, effpot, S_in, etot, &
       & displacement, strain, lwf, energy_table)
    !class (spin_mover_t), intent(inout):: self
    class(spin_mover_t), intent(inout):: self
    class(abstract_potential_t), intent(inout) :: effpot

    real(dp), optional, intent(inout):: displacement(:,:), &
         strain(:,:), lwf(:)
    real(dp), intent(inout) :: S_in(3,self%nspin)
    real(dp), intent(out) ::  etot
    type(hash_table_t),optional, intent(inout) :: energy_table
    integer :: i
    real(dp) ::  Htmp(3), Ri(3)

    ! predict
    etot=0.0
    self%Heff_tmp(:,:)=0.0
    call effpot%calculate(displacement=displacement, strain=strain, lwf=lwf, spin=S_in, &
         & bfield=self%Heff_tmp, energy=etot, energy_table=energy_table)
    call self%get_Langevin_Heff(self%H_lang)
    do i=self%mps%istart, self%mps%iend
       Htmp=self%Heff_tmp(:,i)+self%H_lang(:,i)
       !Ri = cross(S_in(:,i),Htmp)
       !dSdt = -self%gamma_L(i)*(Ri+self%damping(i)* cross(S_in(:,i), Ri))
       Ri=S_in(:,i)!+dSdt*self%dt
       Ri=Ri/sqrt(Ri(1)*Ri(1)+Ri(2)*Ri(2)+Ri(3)*Ri(3))
       self%Stmp(:,i)=Ri
    end do
    call self%mps%allgatherv_dp2d(self%Stmp2, 3, buffer=self%buffer)

  end subroutine spin_mover_t_run_one_step_dummy
  !!***



  !----------------------------------------------------------------------
  !> @brief rotate spin with a rotation matrix
  !>
  !> @param[in]  S_in: input spin array(3)
  !> @param[in]  Heff: effective field
  !> @param[in]  dt: time step
  !> @param[out]  S_out: output spin
  !----------------------------------------------------------------------

  pure function rotate_S_DM(S_in, Heff, dt) result(S_out)
    ! Depondt & Mertens method to rotate S_in
    real(dp), intent(in) :: S_in(3), Heff(3), dt
    real(dp) :: S_out(3)
    real(dp) :: B(3) , w, u, Bnorm, R(3,3), cosw, sinw
    Bnorm=sqrt(sum(Heff*Heff)) 
    B(:)=Heff(:)/Bnorm   ! axis of rotation
    w=Bnorm*dt             ! amplitude of rotation
    sinw=sin(w)
    cosw=cos(w)
    u=1.0d0-cosw

    ! R is rotation matrix
    R(1,1)=B(1)*B(1)*u+cosw
    R(2,1)=B(1)*B(2)*u+B(3)*sinw
    R(3,1)=B(1)*B(3)*u-B(2)*sinw

    R(1,2)=B(1)*B(2)*u-B(3)*sinw
    R(2,2)=B(2)*B(2)*u+cosw
    R(3,2)=B(2)*B(3)*u+B(1)*sinw

    R(1,3)=B(1)*B(3)*u+B(2)*sinw
    R(2,3)=B(2)*B(3)*u-B(1)*sinw
    R(3,3)=B(3)*B(3)*u+cosw
    ! rotate
    S_out=matmul(R, S_in)
  end function rotate_S_DM

  subroutine spin_mover_t_run_one_step_DM(self, effpot, S_in, etot, displacement, strain,&
       & lwf, energy_table)
    ! Depondt & Mertens (2009) method, using a rotation matrix so length doesn't change.
    !class (spin_mover_t), intent(inout):: self
    class(spin_mover_t), intent(inout):: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:)
    real(dp), intent(inout) :: S_in(3,self%nspin)
    real(dp), intent(out) ::  etot
    type(hash_table_t),optional, intent(inout) :: energy_table
    real(dp) :: Htmp(3)
    integer :: i

    ! predict
    self%Stmp2(:,:)=0.0_dp
    etot=0.0
    self%Heff_tmp(:,:)=0.0_dp
    call effpot%calculate(displacement=displacement, strain=strain, lwf=lwf,spin=S_in, &
         bfield=self%Heff_tmp, energy=etot, energy_table=energy_table)
    call self%get_Langevin_Heff(self%H_lang)
    do i=self%mps%istart, self%mps%iend
       Htmp=self%Heff_tmp(:,i)+self%H_lang(:,i)
       ! Note that there is no - , because dsdt =-cross (S, Hrotate) 
       self%Hrotate(:,i) = self%gamma_L(i) * (Htmp + self%damping(i)* cross(S_in(:,i), Htmp))
       self%Stmp2(:,i)= rotate_S_DM(S_in(:,i), self%Hrotate(:,i), self%dt)
    end do
    call self%mps%allgatherv_dp2d(self%Stmp2, 3, self%buffer)

    ! correction
    self%Htmp(:,:)=0.0_dp
    etot=0.0
    call effpot%calculate(displacement=displacement, strain=strain, lwf=lwf, spin=self%Stmp2, &
         bfield=self%Htmp, energy=etot, energy_table=energy_table)

    do i=self%mps%istart, self%mps%iend
       Htmp=(self%Heff_tmp(:,i)+self%Htmp(:,i))*0.5_dp + self%H_lang(:,i)
       self%Hrotate(:,i) = self%gamma_L(i) * (Htmp + self%damping(i)* cross(S_in(:,i), Htmp))
       self%Stmp(:, i)= rotate_S_DM(S_in(:,i), self%Hrotate(:,i), self%dt)
    end do
    call self%mps%allgatherv_dp2d(self%Stmp, 3, self%buffer)
  end subroutine spin_mover_t_run_one_step_DM

  subroutine spin_mover_t_run_one_step_MC(self, effpot, S_in,  etot, displacement, strain,  lwf, energy_table)
    class(spin_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:, :), strain(:,:), lwf(:)
    real(dp), intent(inout) :: S_in(3,self%nspin)
    real(dp), intent(out) ::  etot
    type(hash_table_t),optional, intent(inout) :: energy_table
    if(present(displacement) .or. present(lwf) .or. present(strain)) then
       MSG_BUG("Monte Carlo only implemented for spin.")
       call self%spin_mc%run_MC(self%rng, effpot, S_in, etot)
    end if
    call energy_table%put(self%label, etot)
  end subroutine spin_mover_t_run_one_step_MC

  subroutine spin_mover_t_run_one_step(self, effpot, displacement, strain, spin, lwf, energy_table)
    class(spin_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: effpot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp) ::  etot
    type(hash_table_t),optional, intent(inout) :: energy_table

    if(present(spin)) MSG_ERROR("spin should not be input for spin mover.")
    if(self%method==1) then
       call self%run_one_step_HeunP(effpot=effpot, S_in=self%Stmp, etot=etot, &
            displacement=displacement, strain=strain, lwf=lwf, energy_table=energy_table)
    else if (self%method==2) then
       call self%run_one_step_DM(effpot=effpot, S_in=self%Stmp, etot=etot,&
            displacement=displacement, strain=strain, lwf=lwf, energy_table=energy_table)
    else if (self%method==3) then
       if(present(displacement) .or. present(strain) .or. present(lwf)) then
          MSG_ERROR("Monte carlo not implemented for lattice and lwf yet.")
       endif
       call self%run_one_step_MC(effpot, self%Stmp, etot, energy_table=energy_table)
    else if (self%method==20) then
       call self%run_one_step_dummy(effpot=effpot, S_in=self%Stmp, etot=etot, &
            displacement=displacement, strain=strain, lwf=lwf, energy_table=energy_table)
    end if

    ! do not inc until time is set to hist.
    ! run one step does not know about time. So it will be done in the outer loop.
    if(self%mps%irank==0) then
       call self%hist%set_vars(S=self%Stmp, Snorm=effpot%supercell%spin%ms, etot=etot, inc=.False.)
    end if
  end subroutine spin_mover_t_run_one_step

  !!****f* m_spin_mover/spin_mover_t_run_time
  !!
  !! NAME
  !!  spin_mover_t_run_time
  !!
  !! FUNCTION
  !! run all spin step
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_mover_t_run_time(self, calculator, displacement, strain, spin, lwf, energy_table)

    class(spin_mover_t), intent(inout):: self
    class(abstract_potential_t), intent(inout) :: calculator

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:), spin(:,:)

    type(hash_table_t),optional, intent(inout) :: energy_table
    !type(spin_hist_t), intent(inout) :: hist
    !type(spin_ncfile_t), intent(inout) :: ncfile
    !type(spin_observable_t), intent(inout) :: ob
    !real(dp) ::  S(3, self%nspin)
    real(dp):: t
    integer :: counter, i, ii
    character(len=80) :: msg, msg_empty

    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    t=0.0
    counter=0
    if(iam_master) then
       msg_empty=ch10

       msg=repeat("=", 80)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')
       write(msg, '(A20)') "Spin dynamic steps:"
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

    if (abs(self%thermal_time) > 1e-30) then
       if (iam_master) then
          msg="Thermalization run:"
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
       end if

       do while(t<self%thermal_time)
          counter=counter+1
          call self%run_one_step(effpot=calculator, displacement=displacement, strain=strain, &
               & lwf=lwf, energy_table=energy_table)
          if (iam_master) then
             call self%hist%set_vars( time=t,  inc=.True.)
             if(mod(counter, self%hist%spin_nctime)==0) then
                call self%spin_ob%get_observables( self%hist%S(:,:, self%hist%ihist_prev), &
                     self%hist%Snorm(:,self%hist%ihist_prev),self%hist%etot(self%hist%ihist_prev))
                write(msg, "(A1, 1X, I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") "-", counter, t*Time_Sec, &
                     & self%spin_ob%Mst_norm_total/self%spin_ob%Snorm_total, &
                     & self%hist%etot(self%hist%ihist_prev)/self%spin_ob%nscell
                ! total : 13+4+...= 64 
                call wrtout(std_out,msg,'COLL')
                call wrtout(ab_out, msg, 'COLL')
             endif
          end if
          t=t+self%dt
       end do

       t=0.0
       counter=0
       if (iam_master) then
          call self%hist%reset(array_to_zero=.False.)
          msg="Measurement run:"
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
       end if
    endif
    if(iam_master) then
       call self%spin_ob%reset()
    endif

    do while(t<self%total_time)
       counter=counter+1
       call self%run_one_step(effpot=calculator, displacement=displacement, strain=strain, &
            & spin=spin, lwf=lwf, energy_table=energy_table)
       if (iam_master) then
          call self%hist%set_vars(time=t,  inc=.True.)
          call self%spin_ob%get_observables(self%hist%S(:,:, self%hist%ihist_prev), &
               self%hist%Snorm(:,self%hist%ihist_prev), self%hist%etot(self%hist%ihist_prev))
          if(modulo(counter, self%hist%spin_nctime)==0) then
             call self%spin_ncfile%write_one_step(self%hist)
             write(msg, "(A1, 1X, I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") "-", counter, t*Time_Sec, &
                  & self%spin_ob%Mst_norm_total/self%spin_ob%Snorm_total, &
                  & self%hist%etot(self%hist%ihist_prev)/self%spin_ob%nscell
             call wrtout(std_out,msg,'COLL')
             call wrtout(ab_out, msg, 'COLL')
          endif
       end if
       t=t+self%dt
    enddo

    if (iam_master) then
       msg=repeat("-", 80)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       write(msg, "(A27)") "Summary of spin dynamics:"
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       write(msg, "(A65)") "At the end of the run, the average spin at each sublattice is"
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       write(msg, "(6X, A10, 5X, 3A10, A11)")  'Sublattice', '<M_i>(x)', '<M_i>(y)', '<M_i>(z)', '||<M_i>||'
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       do i =1, self%spin_ob%nsublatt
          write(msg, "(A1, 5X, 2X, I5.4, 8X, 4F10.5)") '-', i, &
               (self%spin_ob%Mst_sub(ii,i)/self%spin_ob%nspin_sub(i)/mu_B , ii=1, 3), &
               sqrt(sum((self%spin_ob%Mst_sub(:, i)/self%spin_ob%nspin_sub(i)/mu_B)**2))
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
       end do

       call wrtout(std_out,msg_empty,'COLL')
       call wrtout(ab_out, msg_empty, 'COLL')

       write(msg, "(A1, 1X, A11, 3X, A13, 3X, A13, 3X, A13, 3X, A13 )" ) &
            "#", "Temperature", "Cv", "chi",  "BinderU4", "Mst"
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")
       write(msg, "(2X, F11.5, 3X, ES13.5, 3X, ES13.5, 3X, E13.5, 3X, ES13.5, 3X )" ) &
            self%temperature*Ha_K , self%spin_ob%Cv, self%spin_ob%chi, &
            self%spin_ob%binderU4, self%spin_ob%Avg_Mst_norm_total/self%spin_ob%snorm_total
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out,  msg, "COLL")

       msg=repeat("=", 80)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')
    end if
  end subroutine spin_mover_t_run_time
  !!***




  !!****f* m_spin_mover/run_MvT
  !!
  !! NAME
  !! run_MvT
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
  !!
  !! CHILDREN
  !!
  !!
  !! SOURCE
  subroutine  run_MvT(self, pot, ncfile_prefix, displacement, strain, spin, lwf, energy_table)
    class(spin_mover_t), intent(inout) :: self
    class(abstract_potential_t), intent(inout) :: pot
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), lwf(:), spin(:,:)
    character(fnlen), intent(inout) :: ncfile_prefix
    type(hash_table_t), optional, intent(inout) :: energy_table
    real(dp) :: T_start, T_end
    integer :: T_nstep
    type(spin_ncfile_t) :: spin_ncfile
    character(len=4) :: post_fname
    real(dp) :: T, T_step
    integer :: i, ii, Tfile, iostat
    character(len=90) :: msg
    character(len=4200) :: Tmsg ! to write to var T file
    character(len=150) :: iomsg
    character(fnlen) :: Tfname ! file name for output various T calculation
    real(dp), allocatable :: Tlist(:), chi_list(:), Cv_list(:), binderU4_list(:)
    real(dp), allocatable :: Mst_sub_norm_list(:, :)
    real(dp), allocatable ::  Mst_norm_total_list(:)

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if (iam_master) then
       T_start=self%params%spin_temperature_start
       T_end=self%params%spin_temperature_end
       T_nstep=self%params%spin_temperature_nstep
       Tfile=get_unit()
       Tfname = trim(ncfile_prefix)//'.varT'
       iostat=open_file(file=Tfname, unit=Tfile, iomsg=iomsg )
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

       ABI_ALLOCATE(Tlist, (T_nstep))
       ABI_ALLOCATE(chi_list, (T_nstep))
       ABI_ALLOCATE(Cv_list, (T_nstep))
       ABI_ALLOCATE(binderU4_list, (T_nstep))
       ABI_ALLOCATE(Mst_sub_norm_list, (self%spin_ob%nsublatt, T_nstep))
       ABI_ALLOCATE( Mst_norm_total_list, (T_nstep))
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
          self%params%spin_temperature=T
       endif
       call self%set_temperature(temperature=T)
       if(iam_master) then
          call self%hist%set_params(spin_nctime=self%params%spin_nctime, &
               &     spin_temperature=T)
          call self%spin_ob%reset(self%params)
          ! uncomment if then to use spin initializer at every temperature. otherwise use last temperature
          if(i==0) then
             call self%set_initial_state()
          else
             call self%hist%inc1()
          endif

          write(post_fname, "(I4.4)") i
          call self%prepare_ncfile( self%params, &
               & trim(ncfile_prefix)//'_T'//post_fname//'_spinhist.nc')
          call spin_ncfile%write_one_step(self%hist)
       endif

       ! run in parallel
       call self%run_time(pot, displacement=displacement, strain=strain, spin=spin, &
            & lwf=lwf, energy_table=energy_table)

       if(iam_master) then
          call spin_ncfile%close()
          ! save observables
          Tlist(i)=T
          chi_list(i)=self%spin_ob%chi
          Cv_list(i)=self%spin_ob%Cv
          binderU4_list(i)=self%spin_ob%binderU4
          !Mst_sub_list(:,:,i)=self%spin_ob%Mst_sub(:,:)  ! not useful
          Mst_sub_norm_list(:,i)=self%spin_ob%Avg_Mst_sub_norm(:)
          Mst_norm_total_list(i)=self%spin_ob%Avg_Mst_norm_total
       endif
    end do


    if(iam_master) then
       ! write summary of MvT run
       msg=repeat("=", 79)
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")

       write(msg, *) "Summary of various T run: "
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")

       write(msg, "(A1, 1X, A11, 3X, A13, 3X, A13, 3X, A13, 3X, A13)" ) &
            "#", "Temperature", "Cv", "chi",  "BinderU4", "Mst"
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out,  msg, "COLL")

       do i = 1, T_nstep
          write(msg, "(2X, F11.5, 3X, ES13.5, 3X, ES13.5, 3X, E13.5, 3X, ES13.5 )" ) &
               Tlist(i)*Ha_K, Cv_list(i), chi_list(i),  binderU4_list(i), Mst_norm_total_list(i)/self%spin_ob%snorm_total
          call wrtout(std_out, msg, "COLL")
          call wrtout(ab_out, msg, "COLL")
       end do

       msg=repeat("=", 79)
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")


       ! write to .varT file
       write(Tmsg, "(A1, 1X, A11, 3X, A13, 3X, A13, 3X, A13, 3X, A13, 3X, *(I13, 3X) )" ) &
            "#", "Temperature (K)", "Cv (1)", "chi (1)",  "BinderU4 (1)", "Mst/Ms(1)", (ii, ii=1, self%spin_ob%nsublatt)
       call wrtout(Tfile, Tmsg, "COLL")

       do i = 1, T_nstep
          write(Tmsg, "(2X, F11.5, 3X, ES13.5, 3X, ES13.5, 3X, E13.5, 3X, ES13.5, 3X, *(ES13.5, 3X) )" ) &
               Tlist(i)*Ha_K, Cv_list(i), chi_list(i),  binderU4_list(i), Mst_norm_total_list(i)/self%spin_ob%snorm_total,&
               & (Mst_sub_norm_list(ii,i)/mu_B, ii=1, self%spin_ob%nsublatt)
          call wrtout(Tfile, Tmsg, "COLL")
       end do
       iostat= close_unit(unit=Tfile, iomsg=iomsg)

       ABI_DEALLOCATE(Tlist)
       ABI_DEALLOCATE(chi_list)
       ABI_DEALLOCATE(Cv_list)
       ABI_DEALLOCATE(binderU4_list)
       ABI_DEALLOCATE(Mst_sub_norm_list)
       ABI_DEALLOCATE( Mst_norm_total_list)

    endif
  end subroutine run_MvT
  !!***

  !!****f* m_spin_mover/current_spin
  !!
  !! NAME
  !! current_spin
  !!
  !! FUNCTION
  !! return the current spin state
  !!
  !! INPUTS
  !! 
  !!
  !! OUTPUT
  !!
  !! PARENTS
  !!
  !!
  !! CHILDREN
  !!
  !!
  !! SOURCE
  function current_spin(self) result(ret)
    class(spin_mover_t), target, intent(inout) :: self
    real(dp), pointer :: ret(:,:)
    integer :: i
    i=self%hist%findIndex(step=0)
    ret => self%hist%S(:,:,i)
  end function current_spin
  !!***




  !!****f* m_spin_mover/finalize
  !!
  !! NAME
  !! finalize
  !!
  !! FUNCTION
  !! finalize spin mover.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!   does nothing. But it's better to preserve initialize-finalize symmetry.
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine finalize(self)

    class(spin_mover_t), intent(inout):: self
    if(allocated(self%gyro_ratio) ) then
       ABI_DEALLOCATE(self%gyro_ratio)
    end if

    if(allocated(self%damping) ) then
       ABI_DEALLOCATE(self%damping)
    end if

    if(allocated(self%gamma_l) ) then
       ABI_DEALLOCATE(self%gamma_l)
    end if

    if(allocated(self%H_lang_coeff) ) then
       ABI_DEALLOCATE(self%H_lang_coeff)
    end if

    if(allocated(self%ms) ) then
       ABI_DEALLOCATE(self%ms)
    end if



    if(self%method==3) then
       call self%spin_mc%finalize()
    end if

    if(allocated(self%Stmp)) then
       ABI_DEALLOCATE(self%Stmp)
    end if

    if(allocated(self%Stmp2)) then
       ABI_DEALLOCATE(self%Stmp2)
    end if


    if(allocated(self%Heff_tmp)) then
       ABI_DEALLOCATE(self%Heff_tmp)
    end if

    if(allocated(self%Htmp)) then
       ABI_DEALLOCATE(self%Htmp)
    end if

    if(allocated(self%Hrotate)) then
       ABI_DEALLOCATE(self%Hrotate)
    end if

    if(allocated(self%H_lang)) then
       ABI_DEALLOCATE(self%H_lang)
    end if

    if(allocated(self%buffer)) then
       ABI_DEALLOCATE(self%buffer)
    end if


    nullify(self%supercell)
    nullify(self%params)
    nullify(self%rng)
    call self%mps%finalize()
    call self%hist%finalize()
    call self%spin_ob%finalize()
    !call self%spin_ncfile%close()
    call self%spin_ob%finalize()
  end subroutine finalize
  !!***

end module m_spin_mover

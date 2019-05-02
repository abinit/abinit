!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_model
!! NAME
!! m_spin_model
!!
!! FUNCTION
!! This module contains the manager of spin dynamics. It call other
!! related modules to run spin dynamics. Itself does nothing in detail.
!!
!!
!! Datatypes:
!!
!! * spin_model_t
!!
!! Subroutines:
!!
!! * spin_model_t
!! * spin_model_t_run
!! * spin_model_t_initialize
!! * spin_model_t_set_initial_spin
!! * spin_model_t_finalize
!! * spin_model_t_read_xml
!! * spin_model_t_make_supercell
!! * spin_model_t_run_one_step
!! * spin_model_t_run_time
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

module m_spin_model

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_io_tools, only : get_unit, open_file, close_unit

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_spin_terms, only: spin_terms_t, spin_terms_t_set_params,  spin_terms_t_finalize, &
       & spin_terms_t_set_external_hfield, spin_terms_t_add_SIA
  use m_spin_model_primitive, only: spin_model_primitive_t, &
       & spin_model_primitive_t_initialize, &
       & spin_model_primitive_t_print_terms, &
       & spin_model_primitive_t_finalize, &
       & spin_model_primitive_t_read_xml, &
       & spin_model_primitive_t_make_supercell
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_init, spin_hist_t_get_s, spin_hist_t_free, &
       & spin_hist_t_set_params, spin_hist_t_reset, spin_hist_t_inc
  use m_spin_mover, only: spin_mover_t, spin_mover_t_initialize, spin_mover_t_finalize, &
       & spin_mover_t_run_time, spin_mover_t_run_one_step
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_init, spin_ncfile_t_close, spin_ncfile_t_def_sd, &
       & spin_ncfile_t_write_primitive_cell, spin_ncfile_t_write_supercell, spin_ncfile_t_write_parameters, &
       & spin_ncfile_t_write_one_step, spin_ncfile_t_def_ob
  use m_spin_observables, only: spin_observable_t, ob_reset, ob_initialize, ob_finalize, ob_calc_observables

  implicit none
  !!***

  !!****t* m_spin_model/spin_model_t
  !! NAME
  !! spin_model_t
  !!
  !! FUNCTION
  !! This type contains all the data for spin dynamics.
  !!
  !! It contains:
  !! spin_primitive = information of the primitive cell, which is a map to the xml file.
  !! spin_calculator = the calculator for spin dynamics, include information of supercell.
  !! spin_mover = include information for how to run spin dynamics
  !! params = parameters from input file
  !! spin_ncfile = wrapper for spin hist netcdf file.
  !! nspins= number of magnetic atoms
  !! SOURCE

  type spin_model_t
     type(spin_model_primitive_t) :: spin_primitive
     type(spin_terms_t) :: spin_calculator
     type(spin_hist_t):: spin_hist
     type(spin_mover_t):: spin_mover
     type(multibinit_dtset_type) :: params
     type(spin_ncfile_t) :: spin_ncfile
     type(spin_observable_t) :: spin_ob
     integer :: nspins
     character(len=fnlen) :: in_fname, out_fname, xml_fname
     !  CONTAINS
     !    procedure :: run => spin_model_t_run
     !    procedure :: initialize=>spin_model_t_initialize
     !    procedure :: set_initial_spin => spin_model_t_set_initial_spin
     !    procedure :: finalize => spin_model_t_finalize
     !    procedure :: read_xml => spin_model_t_read_xml
     !    procedure :: make_supercell => spin_model_t_make_supercell
     !    procedure :: run_one_step => spin_model_t_run_one_step
     !    procedure :: run_time => spin_model_t_run_time
     !    procedure :: run_MvT => spin_model_t_run_MvT
  end type spin_model_t
  !!***
contains

  !!****f* m_spin_model/spin_model_t_run
  !!
  !! NAME
  !!  spin_model_t_run
  !!
  !! FUNCTION
  !!  run the whole spin dynamics.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  currently it just call run_time.
  !!  It will be able to call e.g. run_MvT to do other things.
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run(self)

    class(spin_model_t), intent(inout) :: self
    if(self%params%spin_var_temperature==1) then
       call spin_model_t_run_various_T(self, self%params%spin_temperature_start, &
            & self%params%spin_temperature_end, self%params%spin_temperature_nstep)
    else
       call spin_model_t_run_time(self)
    end if
  end subroutine spin_model_t_run
  !!***

  !!****f* m_spin_model/spin_model_t_initialize
  !!
  !! NAME
  !!  spin_model_t_initialize
  !!
  !! FUNCTION
  !! initialize spin dynamics from input
  !!
  !! INPUTS
  !!  xml_fname= file of xml file
  !! params = multibinit dataset from input file
  !! OUTPUT
  !!  self <spin_model_t>
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_initialize(self, filenames,  params)

    class(spin_model_t), intent(inout) :: self
    character(len=fnlen), intent(in) :: filenames(:)
    integer :: sc_matrix(3,3)
    type(multibinit_dtset_type), intent(in) :: params

    self%in_fname=filenames(1)
    self%out_fname=filenames(2)
    self%xml_fname=filenames(3)

    ! read input
    !call self%spin_primitive%initialize()
    call spin_model_primitive_t_initialize(self%spin_primitive)
    self%params=params
    ! TODO: remove this and use a.u. everywhere.
    call spin_model_t_unit_conversion(self)
    !call self%read_xml(xml_fname)
    call spin_model_t_read_xml(self, trim(self%xml_fname)//char(0))
    !call self%spin_primitive%print_terms()
    call spin_model_primitive_t_print_terms(self%spin_primitive)

    ! make supercell
    sc_matrix(:,:)=0.0_dp
    sc_matrix(1,1)=self%params%ncell(1)
    sc_matrix(2,2)=self%params%ncell(2)
    sc_matrix(3,3)=self%params%ncell(3)
    !call self%make_supercell(sc_matrix)
    call spin_model_t_make_supercell(self, sc_matrix)

    ! set parameters to hamiltonian and mover
    self%nspins= self%spin_calculator%nspins

    call spin_model_t_set_params(self)

    !TODO hexu: mxhist, has_latt, natoms should be input with their true values when lattice part also added
    call spin_hist_t_init(hist=self%spin_hist, nspins=self%nspins, &
         &   mxhist=3, has_latt=.False.)

    call spin_hist_t_set_params(self%spin_hist, spin_nctime=self%params%spin_nctime, &
         &     spin_temperature=self%params%spin_temperature)

    !call self%set_initial_spin(mode=1)
    call spin_model_t_set_initial_spin(self)

    !call self%spin_mover%initialize(self%nspins, dt=params%dtspin, & 
    !    & total_time=params%dtspin*params%ntime_spin, temperature=self%params%self)
    call spin_mover_t_initialize(self%spin_mover, self%nspins, dt=self%params%spin_dt, &
         &  total_time=self%params%spin_dt*self%params%spin_ntime, temperature=self%params%spin_temperature, &
         & pre_time=self%params%spin_dt*self%params%spin_ntime_pre, method=self%params%spin_dynamics)

    call ob_initialize(self%spin_ob, self%spin_calculator, self%params)

    call spin_model_t_prepare_ncfile(self, self%spin_ncfile, trim(self%out_fname)//'_spinhist.nc')

    call spin_ncfile_t_write_one_step(self%spin_ncfile, self%spin_hist)
  end subroutine spin_model_t_initialize
  !!***



  !!****f* m_spin_model/spin_model_t_finalize
  !!
  !! NAME
  !!  spin_model_t_finalize
  !!
  !! FUNCTION
  !! clean memory for spin dynamics
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
  subroutine spin_model_t_finalize(self)

    class(spin_model_t), intent(inout) :: self
    !call self%spin_primitive%finalize()
    call spin_model_primitive_t_finalize(self%spin_primitive)
    call spin_terms_t_finalize(self%spin_calculator)
    call spin_hist_t_free(self%spin_hist)
    call spin_mover_t_finalize(self%spin_mover)
    call ob_finalize(self%spin_ob)
    call spin_ncfile_t_close(self%spin_ncfile)
    ! TODO: finalize others
  end subroutine spin_model_t_finalize
  !!***

  !!****f* m_spin_model/spin_model_t_set_params
  !!
  !! NAME
  !!  spin_model_t_set_params
  !!
  !! FUNCTION
  !! set parameters for spin dynamics from input.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !! Most are already done in initialize. (should be moved to here or remove this function)
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_set_params(self)

    class(spin_model_t), intent(inout) :: self
    real(dp):: mfield(3, self%nspins), damping(self%nspins)
    integer ::  i
    ! params -> mover
    ! params -> calculator


    if (self%params%spin_damping >=0) then
       damping(:)= self%params%spin_damping
       call spin_terms_t_set_params(self%spin_calculator, dt=self%params%spin_dt, &
            & temperature=self%params%spin_temperature, gilbert_damping=damping)
    else
       call spin_terms_t_set_params(self%spin_calculator, dt=self%params%spin_dt, &
            & temperature=self%params%spin_temperature)
    end if

    do i=1, self%nspins
       mfield(:, i)=self%params%spin_mag_field(:)
    end do
    call spin_terms_t_set_external_hfield(self%spin_calculator, mfield)

    if (self%params%spin_sia_add /= 0 ) then
       call spin_terms_t_add_SIA(self%spin_calculator, self%params%spin_sia_add, &
            &  self%params%spin_sia_k1amp, self%params%spin_sia_k1dir)
    end if

    ! params -> hist

  end subroutine spin_model_t_set_params
  !!***

  !!****f* m_spin_model/spin_model_t_read_xml
  !!
  !! NAME
  !!  spin_model_t_read_xml
  !!
  !! FUNCTION
  !! read xml file and set primitive cell parameters.
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
  subroutine spin_model_t_read_xml(self, xml_fname)

    class(spin_model_t), intent(inout) :: self
    character(len=*), intent(in) :: xml_fname
    logical:: use_sia, use_exchange, use_dmi, use_bi
    !call self%spin_primitive%read_xml(xml_fname)
    use_exchange=.True.
    use_sia=.True.
    use_dmi=.True.
    use_bi=.True.

    ! Do not use sia term in xml if spin_sia_add is set to 1.
    if(self%params%spin_sia_add == 1) use_sia=.False.

    call spin_model_primitive_t_read_xml(self%spin_primitive, xml_fname, &
         & use_exchange=use_exchange,  use_sia=use_sia, use_dmi=use_dmi, use_bi=use_bi)
  end subroutine spin_model_t_read_xml
  !!***


  !!****f* m_spin_model/spin_model_t_make_supercell
  !!
  !! NAME
  !!  spin_model_t_make_supercell
  !!
  !! FUNCTION
  !!  make supercell and prepare calculator
  !!
  !! INPUTS
  !!  sc_mat(3,3) = supercell matrix
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_make_supercell(self, sc_mat)

    class(spin_model_t), intent(inout) :: self
    integer , intent(in):: sc_mat(3, 3)
    !call self%spin_primitive%make_supercell(sc_mat, self%spin_calculator)
    call spin_model_primitive_t_make_supercell(self%spin_primitive, sc_mat, self%spin_calculator)
  end subroutine spin_model_t_make_supercell
  !!***

  !!****f* m_spin_model/spin_model_t_prepare_ncfile
  !!
  !! FUNCTION
  !!  set initial spin state
  !!
  !! INPUTS
  !!  mode= 0 : all along z. 1: random
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_prepare_ncfile(self, spin_ncfile, fname)

    class(spin_model_t), intent(inout) :: self
    type(spin_ncfile_t), intent(out) :: spin_ncfile
    character(len=*) :: fname
    call spin_ncfile_t_init(spin_ncfile, trim(fname), self%params%spin_write_traj)
    call spin_ncfile_t_def_sd(spin_ncfile, self%spin_hist )
    call spin_ncfile_t_def_ob(spin_ncfile, self%spin_ob)
    !call spin_ncfile_t_write_primitive_cell(self%spin_ncfile, self%spin_primitive)
    call spin_ncfile_t_write_supercell(spin_ncfile, self%spin_calculator)
    call spin_ncfile_t_write_parameters(spin_ncfile, self%params)
  end subroutine spin_model_t_prepare_ncfile
  !!***


  !!****f* m_spin_model/spin_model_t_set_initial_spin
  !!
  !! NAME
  !!  spin_model_t_set_initial_spin
  !!
  !! FUNCTION
  !!  set initial spin state
  !!
  !! INPUTS
  !!  mode= 0 : all along z. 1: random
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_set_initial_spin(self)

    class(spin_model_t), intent(inout) :: self
    integer :: mode
    integer :: i
    real(dp) :: S(3, self%nspins)
    character(len=500) :: msg
    mode=self%params%spin_init_state
    if(mode==2) then
       ! set all spin to z direction.
       S(1,:)=0.0d0
       S(2,:)=0.0d0
       S(3,:)=1.0d0
    else if (mode==1) then
       ! randomize S using uniform random number
       ! print *, "Initial spin set to random value"
       write(msg,*) "Initial spin set to random value."
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
       call random_number(S)
       S=S-0.5
       do i=1, self%nspins
          S(:,i)=S(:,i)/sqrt(sum(S(:, i)**2))
       end do
    else
       write(msg,*) "Error: Set initial spin: mode should be 2(FM) or 1 (random). Others are not yet implemented."
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')

    end if

    call spin_hist_t_set_vars(self%spin_hist, S=S, Snorm=self%spin_calculator%ms, &
         &  time=0.0_dp, ihist_latt=0, inc=.True.)
  end subroutine spin_model_t_set_initial_spin
  !!***


  !!****f* m_spin_model/spin_model_t_run_one_step
  !!
  !! NAME
  !!  spin_model_t_run_one_step
  !!
  !! FUNCTION
  !!  run one spin dynamics step
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  This is not called by run_time. insted run_time call mover%run_time
  !!  (in the spirit that this module should do nothing in detail!)
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run_one_step(self)

    class(spin_model_t), intent(inout) :: self
    call spin_mover_t_run_one_step(self%spin_mover, self%spin_calculator, &
         self%spin_hist)
  end subroutine spin_model_t_run_one_step
  !!***

  !!****f* m_spin_model/spin_model_t_run_time
  !!
  !! NAME
  !!  spin_model_t_run_time
  !!
  !! FUNCTION
  !!  run all spin time step
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
  subroutine spin_model_t_run_time(self)

    class(spin_model_t), intent(inout) :: self
    call spin_mover_t_run_time(self%spin_mover,self%spin_calculator,self%spin_hist, self%spin_ncfile, self%spin_ob)
  end subroutine spin_model_t_run_time
  !!***


  subroutine spin_model_t_run_various_T(self, T_start, T_end, T_nstep)

    class(spin_model_t), intent(inout) :: self
    real(dp), intent(in) :: T_start, T_end
    integer, intent(in) :: T_nstep
    type(spin_ncfile_t) :: spin_ncfile
    character(len=4) :: post_fname
    real(dp) :: T, T_step
    integer :: i, ii, Tfile, iostat !, sublatt(self%spin_ob%nsublatt)
    character(len=90) :: msg
    character(len=4200) :: Tmsg ! to write to var T file
    character(len=150) :: iomsg
    character(fnlen) :: Tfname ! file name for output various T calculation


    real(dp) :: Tlist(T_nstep), chi_list(T_nstep), Cv_list(T_nstep), binderU4_list(T_nstep)
    !real(dp) :: Mst_sub_list(3, self%spin_ob%nsublatt, T_nstep)
    real(dp) :: Mst_sub_norm_list(self%spin_ob%nsublatt, T_nstep)
    real(dp) ::  Mst_norm_total_list(T_nstep)

    Tfile=get_unit()
    Tfname = trim(self%out_fname)//'.varT'
    iostat=open_file(file=Tfname, unit=Tfile, iomsg=iomsg )

    !do i=1, self%spin_ob%nsublatt
    !   sublatt(i)=i
    !end do
    T_step=(T_end-T_start)/(T_nstep-1)

    write(msg, "(A52, ES13.5, A11, ES13.5, A1)") & 
         & "Starting temperature dependent calculations. T from ", &
         & T_start, "K to ", T_end, " K."
    call wrtout(std_out, msg, "COLL")
    call wrtout(ab_out, msg, "COLL")

    do i=1, T_nstep
       T=T_start+(i-1)*T_step

       msg=repeat("=", 79)
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")

       write(msg, "(A13, 5X, ES13.5, A3)") "Temperature: ", T, " K."
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out,  msg, "COLL")

       call spin_hist_t_reset(self%spin_hist, array_to_zero=.False.)
       ! set temperature
       ! TODO make this into a subroutine set_params
       self%params%spin_temperature=T
       call spin_terms_t_set_params(self%spin_calculator, temperature=T)
       self%spin_mover%temperature=T
       call spin_hist_t_set_params(self%spin_hist, spin_nctime=self%params%spin_nctime, &
            &     spin_temperature=T)
       call ob_reset(self%spin_ob, self%params)
       ! uncomment if then to use spin initializer at every temperature. otherwise use last temperature
       if(i==0) then
          call spin_model_t_set_initial_spin(self)
       else
          call spin_hist_t_inc(self%spin_hist)
       endif

       write(post_fname, "(I4.4)") i+1
       call spin_model_t_prepare_ncfile(self, spin_ncfile, & 
            & trim(self%out_fname)//'_T'//post_fname//'_spinhist.nc')
       call spin_ncfile_t_write_one_step(spin_ncfile, self%spin_hist)
       call spin_mover_t_run_time(self%spin_mover, self%spin_calculator, & 
            & self%spin_hist, ncfile=spin_ncfile, ob=self%spin_ob)

       call spin_ncfile_t_close(spin_ncfile)
       ! save observables
       Tlist(i)=T
       chi_list(i)=self%spin_ob%chi
       Cv_list(i)=self%spin_ob%Cv
       binderU4_list(i)=self%spin_ob%binderU4
       !Mst_sub_list(:,:,i)=self%spin_ob%Mst_sub(:,:)  ! not useful
       Mst_sub_norm_list(:,i)=self%spin_ob%Avg_Mst_sub_norm(:)
       Mst_norm_total_list(i)=self%spin_ob%Avg_Mst_norm_total
    end do


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
            Tlist(i), Cv_list(i), chi_list(i),  binderU4_list(i), Mst_norm_total_list(i)/self%spin_ob%snorm_total
       call wrtout(std_out, msg, "COLL")
       call wrtout(ab_out, msg, "COLL")
    end do

    msg=repeat("=", 79)
    call wrtout(std_out, msg, "COLL")
    call wrtout(ab_out, msg, "COLL")
 

    ! write to .varT file
    write(Tmsg, "(A1, 1X, A11, 3X, A13, 3X, A13, 3X, A13, 3X, A13, 3X, *(I13, 3X) )" ) &
         "#", "Temperature", "Cv", "chi",  "BinderU4", "Mst", (ii, ii=1, self%spin_ob%nsublatt)
    call wrtout(Tfile, Tmsg, "COLL")

    do i = 1, T_nstep
       write(Tmsg, "(2X, F11.5, 3X, ES13.5, 3X, ES13.5, 3X, E13.5, 3X, ES13.5, 3X, *(ES13.5, 3X) )" ) &
            Tlist(i), Cv_list(i), chi_list(i),  binderU4_list(i), Mst_norm_total_list(i)/self%spin_ob%snorm_total,&
            & (Mst_sub_norm_list(ii,i)/mu_B_SI, ii=1, self%spin_ob%nsublatt)
       call wrtout(Tfile, Tmsg, "COLL")
    end do
    iostat= close_unit(unit=Tfile, iomsg=iomsg)

  end subroutine spin_model_t_run_various_T

  !!****f* m_spin_model/spin_model_t_run_MvH
  !!
  !! NAME
  !!  spin_model_t_run_MvT
  !!
  !! FUNCTION
  !!  run spin vs external magnetic field
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!  Not yet implemented.
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine spin_model_t_run_MvH(self)

    class(spin_model_t), intent(inout) :: self
    ABI_UNUSED(self%params%spin_dt)
    write(std_out, *) "MvH is not yet implemented. "
  end subroutine spin_model_t_run_MvH
  !!***

  !! convert unit of input variables into S.I. unit
  !! TODO This is temporary and should be removed
  !!         a.u. should be used internally. 
  subroutine spin_model_t_unit_conversion(self)

    class(spin_model_t), intent(inout) :: self
    self%params%spin_dt = self%params%spin_dt*Time_Sec
  end subroutine spin_model_t_unit_conversion

end module m_spin_model

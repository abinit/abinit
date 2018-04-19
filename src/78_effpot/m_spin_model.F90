

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_spin_model
  use defs_basis
  use m_profiling_abi
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_spin_terms, only: spin_terms_t, spin_terms_t_finalize
  use m_spin_model_primitive, only: spin_model_primitive_t, &
       & spin_model_primitive_t_initialize, &
       & spin_model_primitive_t_print_terms, &
       & spin_model_primitive_t_finalize, &
       & spin_model_primitive_t_read_xml, &
       & spin_model_primitive_t_make_supercell
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_init, spin_hist_t_get_s, spin_hist_t_free, &
       &   spin_hist_t_set_params
  use m_spin_mover, only: spin_mover_t, spin_mover_t_initialize, spin_mover_t_finalize, &
       &  spin_mover_t_run_time, spin_mover_t_run_one_step
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_init, spin_ncfile_t_close, spin_ncfile_t_def_sd, &
       & spin_ncfile_t_write_primitive_cell, spin_ncfile_t_write_supercell, spin_ncfile_t_write_parameters, &
       & spin_ncfile_t_write_one_step
  implicit none
  type spin_model_t
     type(spin_model_primitive_t) :: spin_primitive
     type(spin_terms_t) :: spin_calculator
     type(spin_hist_t):: spin_hist
     type(spin_mover_t):: spin_mover
     type(multibinit_dtset_type) :: params
     type(spin_ncfile_t) :: spin_ncfile
     integer :: nmatoms
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
     !    procedure :: run_MvH => spin_model_t_run_MvH
  end type spin_model_t
contains
  ! This is the main routine of all spin model.
  subroutine spin_model_t_run(self)
    !TODO add input as parameters.

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    integer :: i
    !call spin_mover_t_run_time(self%spin_mover, self%spin_calculator, self%spin_hist)
    call spin_model_t_run_time(self)
  end subroutine spin_model_t_run

  ! initialize
  subroutine spin_model_t_initialize(self, xml_fname,  params)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_initialize'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    character(len=40), intent(in) :: xml_fname
    integer :: sc_matrix(3,3)
    type(multibinit_dtset_type), intent(in) :: params

    ! read input
    !call self%spin_primitive%initialize()
    call spin_model_primitive_t_initialize(self%spin_primitive)

    self%params=params
    !call self%read_xml(xml_fname)
    call spin_model_t_read_xml(self, xml_fname)
    !call self%spin_primitive%print_terms()
    call spin_model_primitive_t_print_terms(self%spin_primitive)

    ! make supercell
    sc_matrix(:,:)=0.0_dp
    sc_matrix(1,1)=params%n_cell(1)
    sc_matrix(2,2)=params%n_cell(2)
    sc_matrix(3,3)=params%n_cell(3)
    !call self%make_supercell(sc_matrix)
    call spin_model_t_make_supercell(self, sc_matrix)

    ! set parameters to hamiltonian and mover
    self%nmatoms= self%spin_calculator%nmatoms
    ! TODO hexu: max_save, step_save should be defined in input file and params

    !TODO hexu: mxhist, has_latt, natom should be input with their true values when lattice part also added
    call spin_hist_t_init(hist=self%spin_hist, nmatom=self%nmatoms, mxhist=3, has_latt=.False.)
    call spin_hist_t_set_params(self%spin_hist, spin_nctime=self%params%spin_nctime, &
            &     spin_temperature=self%params%spin_temperature)
    !TODO
    !call spin_hist_t_set_atomic_structure(self%spin_hist, acell, rprimd, xred, spin_index, ntypat,  typat, znucl)

    !call self%set_initial_spin(mode=1)
    call spin_model_t_set_initial_spin(self, mode=0)

    !call self%spin_mover%initialize(self%nmatoms, dt=params%dtspin, total_time=params%dtspin*params%ntime_spin, temperature=self%params%self)
    call spin_mover_t_initialize(self%spin_mover, self%nmatoms, dt=params%spin_dt, &
         &  total_time=params%spin_dt*params%spin_ntime, temperature=self%params%spin_temperature)

    call spin_ncfile_t_init(self%spin_ncfile, 'spinhist.nc')
    !call spin_ncfile_t_write_primitive_cell(self%spin_ncfile, self%spin_primitive)
    !call spin_ncfile_t_write_supercell(self%spin_ncfile, self%spin_calculator)
    call spin_ncfile_t_def_sd(self%spin_ncfile, self%spin_hist )
    call spin_ncfile_t_write_one_step(self%spin_ncfile, self%spin_hist)
  end subroutine spin_model_t_initialize

  ! finalize
  subroutine spin_model_t_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_finalize'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    !call self%spin_primitive%finalize()
    call spin_model_primitive_t_finalize(self%spin_primitive)
    call spin_terms_t_finalize(self%spin_calculator)
    call spin_hist_t_free(self%spin_hist)
    call spin_mover_t_finalize(self%spin_mover)
    call spin_ncfile_t_close(self%spin_ncfile)
    ! TODO: finalize others
  end subroutine spin_model_t_finalize

  ! read xml file and set primitive cell parameters.
  subroutine spin_model_t_read_xml(self, xml_fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_read_xml'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    character(len=40), intent(in) :: xml_fname
    !call self%spin_primitive%read_xml(xml_fname)
    call spin_model_primitive_t_read_xml(self%spin_primitive, xml_fname)
  end subroutine spin_model_t_read_xml

  ! from primitive cell parameter make supercell and prepare calculator
  subroutine spin_model_t_make_supercell(self, sc_mat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_make_supercell'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    integer , intent(in):: sc_mat(3, 3)
    !call self%spin_primitive%make_supercell(sc_mat, self%spin_calculator)
    call spin_model_primitive_t_make_supercell(self%spin_primitive, sc_mat, self%spin_calculator)
  end subroutine spin_model_t_make_supercell

  subroutine spin_model_t_set_initial_spin(self, mode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_set_initial_spin'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    integer, intent(in) :: mode
    integer :: i
    real(dp) :: S(3, self%nmatoms)
    if(mode==0) then
       ! set all spin to z direction.
       S(1,:)=0.0d0
       S(2,:)=0.0d0
       S(3,:)=1.0d0
    else if (mode==1) then
       ! randomize S using uniform random number
       print *, "Initial spin set to random value"
       call random_number(S)
       do i=1, self%nmatoms
          S(:,i)=S(:,i)/sqrt(sum(S(:, i)**2))
       end do
    else
       print *, "Error: Set initial spin: mode should be 0 (FM) or 1 (random)"
    end if
    !call self%spin_hist%insert(S)
    !call spin_hist_t_insert(self%spin_hist, S)

    ! spin_hist_t_set_vars(hist, S, Snorm, dSdt, Heff, etot, entropy, time, ihist_latt, inc)
    !TODO initialize lattice structure, spin_type, Snorm
    call spin_hist_t_set_vars(self%spin_hist, S=S, time=0.0_dp, ihist_latt=0, inc=.True.)
    !print *, "initial spin", self%spin_hist%current_S
  end subroutine spin_model_t_set_initial_spin

  ! run one time step
  subroutine spin_model_t_run_one_step(self)
!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_one_step'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    real(dp) :: S_tmp(3,self%nmatoms)
    call spin_mover_t_run_one_step(self%spin_mover, self%spin_calculator, &
         spin_hist_t_get_S(self%spin_hist),S_tmp)
    call spin_hist_t_set_vars(self%spin_hist, S=S_tmp, inc=.False.)
  end subroutine spin_model_t_run_one_step

  ! run all time step
  subroutine spin_model_t_run_time(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_time'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    real(dp) :: t
    integer :: i
    t=0.0d0
    i=0

    call spin_mover_t_run_time(self%spin_mover,self%spin_calculator,self%spin_hist, self%spin_ncfile)
    do while(i<self%params%spin_ntime)
       print *, "==================", i
       i=i+1
       t=t+self%params%spin_dt
       !print *, "t= ", t
       !call self%run_one_step()
       !call spin_model_t_run_one_step(self)
       !call spin_hist_t_set_vars(self%spin_hist, time=t, inc=.True.)
       !call spin_ncfile_t_write_one_step(self%spin_ncfile, self%spin_hist)
    end do
  end subroutine spin_model_t_run_time

  ! run spin .vs. temperature
  subroutine spin_model_t_run_MvT(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_MvT'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    !TODO
    write(std_err, *) "MvH is not yet implemented. "
  end subroutine spin_model_t_run_MvT

  ! run spin .vs. external magnetic field.
  subroutine spin_model_t_run_MvH(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_t_run_MvH'
!End of the abilint section

    class(spin_model_t), intent(inout) :: self
    !TODO
    write(std_err, *) "MvH is not yet implemented. "
  end subroutine spin_model_t_run_MvH

end module m_spin_model

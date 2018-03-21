module m_spin_model
  use m_spin_terms
  use m_spin_model_primitive, only: spin_model_primitive_t, &
         & spin_model_primitive_t_initialize, &
         & spin_model_primitive_t_print_terms, &
         & spin_model_primitive_t_finalize, &
         & spin_model_primitive_t_read_xml, & 
         & spin_model_primitive_t_make_supercell

  use m_spin_mover
  use m_spin_hist
  use m_multibinit_dataset
  implicit none
  type spin_model_t
     type(spin_model_primitive_t) :: spin_primitive
     type(spin_terms_t) :: spin_calculator
     type(spin_hist_t):: spin_hist
     type(spin_mover_t):: spin_mover
     type(multibinit_dtset_type) :: params
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
    do i=1, self%params%ntime_spin
        call spin_model_t_run_one_step(self)
    enddo
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
    !call self%spin_hist%initialize(self%nmatoms, max_save=40000, step_save=1000, dt=self%params%dtspin)
    call spin_hist_t_initialize(self%spin_hist, self%nmatoms, max_save=40000, step_save=100, dt=self%params%dtspin)

    !call self%set_initial_spin(mode=1)
    call spin_model_t_set_initial_spin(self, mode=0)

    !call self%spin_mover%initialize(self%nmatoms, dt=params%dtspin, total_time=params%dtspin*params%ntime_spin, temperature=self%params%self)
    call spin_mover_t_initialize(self%spin_mover, self%nmatoms, dt=params%dtspin, &
         &  total_time=params%dtspin*params%ntime_spin, temperature=self%params%temperature_spin)
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
    call spin_hist_t_insert(self%spin_hist, S)
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
    !call self%spin_mover%run_one_step(self%spin_calculator, &
    !     self%spin_hist%current_S,S_tmp)
    call spin_mover_t_run_one_step(self%spin_mover, self%spin_calculator, &
         self%spin_hist%current_S,S_tmp)

    !call self%spin_hist%insert(S_tmp)
    call spin_hist_t_insert(self%spin_hist, S_tmp)
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
    !print *, "params: dtspin", self%params%dtspin
    do while(i<self%params%ntime_spin)
       i=i+1
       t=t+self%params%dtspin
       !print *, "t= ", t
       !call self%run_one_step()
       call spin_model_t_run_one_step(self)
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
  end subroutine spin_model_t_run_MvH

end module m_spin_model

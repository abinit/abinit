module m_spin_mover
  use defs_basis
  use m_spin_terms
  use m_mathfuncs
  use m_spin_hist
  implicit none

  type spin_mover_t
     integer :: nmatoms
     real(dp) :: dt, total_time, temperature
     !CONTAINS
     !   procedure :: initialize => spin_mover_t_initialize
     !   procedure :: run_one_step => spin_mover_t_run_one_step
     !   procedure :: run_time => spin_mover_t_run_time
  end type

contains
  subroutine spin_mover_t_initialize(self, nmatoms, dt, total_time, temperature)
    !class (spin_mover_t):: self

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_initialize'
!End of the abilint section

    type(spin_mover_t), intent(inout) :: self
    real(dp), intent(in) :: dt, total_time, temperature
    integer, intent(in) :: nmatoms
    self%nmatoms=nmatoms
    self%dt=dt
    self%total_time=total_time
    self%temperature=temperature
  end subroutine spin_mover_t_initialize

  ! Heun's integration Method
  subroutine spin_mover_t_run_one_step(self, calculator, S_in, S_out)
    !class (spin_mover_t), intent(inout):: self

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_one_step'
!End of the abilint section

    type(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nmatoms)
    real(dp), intent(out) :: S_out(3,self%nmatoms)
    integer :: i
    real(dp) ::  dSdt(3, self%nmatoms), dSdt2(3, self%nmatoms), &
         & Heff(3, self%nmatoms), Heff2(3, self%nmatoms), &
         & H_lang(3, self%nmatoms)
    ! predict

    !call calculator%get_Langevin_Heff(self%dt, self%temperature, H_lang)
    call spin_terms_t_get_Langevin_Heff(calculator, self%dt, self%temperature, H_lang)

    !call calculator%get_dSdt(S_in, H_lang, dSdt)
    call spin_terms_t_get_dSdt(calculator, S_in, H_lang, dSdt)
    !$OMP PARALLEL DO
    do i =1, self%nmatoms
       S_out(:,i)=  S_in(:,i) +dSdt(:,i) * self%dt
    end do
    !$OMP END PARALLEL DO

    !do i=1, self%nmatoms
    !   S_out(:,i)=S_out(:,i)/sqrt(sum(S_out(:,i)**2))
    !end do

    ! correction

    !call calculator%get_dSdt(S_out, H_lang, dSdt2)
    call spin_terms_t_get_dSdt(calculator, S_out, H_lang, dSdt2)
    !$OMP PARALLEL DO
    do i =1, self%nmatoms
       S_out(:,i)=  S_in(:,i) +(dSdt(:,i)+dSdt2(:,i)) * (0.5_dp*self%dt)
    end do
    !$OMP END PARALLEL DO
    !print *, "S before norm", S_out
    do i=1, self%nmatoms
       S_out(:,i)=S_out(:,i)/sqrt(sum(S_out(:,i)**2))
    end do
    !print *, "dt: ",self%dt
    !print *, "S:", S_out
  end subroutine spin_mover_t_run_one_step

  subroutine spin_mover_t_run_time(self, calculator, hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_time'
!End of the abilint section

    class(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    type(spin_hist_t), intent(inout) :: hist
    real(dp) ::  S(3, self%nmatoms)
    real(dp):: t
    t=0.0
    do while(t<self%total_time)
       !call self%run_one_step(calculator, hist%current_S, S)
       call spin_mover_t_run_one_step(self, calculator, hist%current_S, S)
       !call hist%insert(S)
       call spin_hist_t_insert(hist, S)
       t=t+self%dt
    enddo
  end subroutine spin_mover_t_run_time

end module m_spin_mover

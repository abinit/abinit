#include "abi_common.h"
module m_spin_hist
  use defs_basis
  implicit none
  type spin_hist_t
     integer:: nmatoms, max_save, step_save, current_step, n_saved_step
     integer, allocatable :: steps(:)
     real(dp), allocatable :: current_S(:, :), last_S(:,:), &
          & saved_S(:,:,:), saved_time(:)
     real(dp), allocatable :: average_S(:,:), average_S_per_site(:,:,:)
     real(dp) :: current_time, last_time, dt, current_average_S
  !  CONTAINS
  !    procedure :: initialize => spin_hist_t_initialize
  !    procedure :: finalize => spin_hist_t_finalize
  !    procedure :: insert => spin_hist_t_insert
  !    procedure :: save_file => spin_hist_t_save_file
  !    procedure :: calc_observables => spin_hist_t_calc_observables
  end type spin_hist_t
CONTAINS
  subroutine spin_hist_t_initialize(self, nmatoms, max_save, step_save, dt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_hist_t_initialize'
!End of the abilint section

    class(spin_hist_t), intent(inout) :: self
    integer, intent(in) :: nmatoms, max_save, step_save
    real(dp), intent(in) :: dt

    self%nmatoms=nmatoms
    self%max_save=max_save
    self%step_save=step_save
    self%dt=dt

    self%current_time=0.0d0
    self%current_step=0
    self%n_saved_step=0

    if(.not. allocated(self%last_S)) then
       ABI_ALLOCATE(self%last_S, (3, self%nmatoms))
    endif
    self%last_S(:,:)=0.0d0

    if(.not. allocated(self%current_S)) then
       ABI_ALLOCATE(self%current_S, (3, self%nmatoms))
    endif
    self%current_S(:,:)=0.0d0

    if(.not. allocated(self%saved_S)) then
       ABI_ALLOCATE(self%saved_S, (3, self%nmatoms, self%max_save))
    endif

    if(.not. allocated(self%saved_time)) then
       ABI_ALLOCATE(self%saved_time, (self%max_save))
    endif

    if(.not. allocated(self%average_S)) then
       ABI_ALLOCATE(self%average_S, (3, self%max_save))
    endif

    if(.not. allocated(self%average_S_per_site)) then
       ABI_ALLOCATE(self%average_S_per_site, (3, self%nmatoms, self%max_save))
    endif

  end subroutine spin_hist_t_initialize

  subroutine spin_hist_t_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_hist_t_finalize'
!End of the abilint section

    class(spin_hist_t), intent(inout) :: self
    self%nmatoms=0
    self%max_save=0
    self%step_save=0
    self%dt=0

    self%current_time=0.0d0
    self%current_step=0
    self%n_saved_step=0

    if( allocated(self%last_S)) then
       ABI_DEALLOCATE(self%last_S)
    endif

    if( allocated(self%current_S)) then
       ABI_DEALLOCATE(self%current_S)
    endif

    if( allocated(self%saved_S)) then
       ABI_DEALLOCATE(self%saved_S)
    endif

    if( allocated(self%saved_time)) then
       ABI_DEALLOCATE(self%saved_time)
    endif

    if( allocated(self%average_S)) then
       ABI_DEALLOCATE(self%average_S)
    endif

    if(allocated(self%average_S_per_site)) then
       ABI_DEALLOCATE(self%average_S_per_site)
    endif

  end subroutine spin_hist_t_finalize

  subroutine spin_hist_t_insert(self, S)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_hist_t_insert'
!End of the abilint section

    class(spin_hist_t), intent(inout) :: self
    real(dp), intent(in) :: S(3, self%nmatoms)
    real(dp), save:: avg_S, total_S
    !print *, "Spin hist: insert new S. current_step= ", self%current_step, mod( self%current_step, self%step_save)
    self%last_S(:,:)=self%current_S(:,:)
    self%current_S(:,:)=S(:,:)
    self%current_step=self%current_step+1
    self%current_time=self%current_time+self%dt
    avg_S=0
    if ( mod( self%current_step-1, self%step_save)==0  ) then
       if ( self%n_saved_step<self%max_save) then
          !print *, "spin saved to hist."
          self%n_saved_step=self%n_saved_step+1
          self%saved_S(:,:,self%n_saved_step)=S(:,:)
          self%saved_time(self%n_saved_step)=self%current_time
          !call self%calc_observables()
          call spin_hist_t_calc_observables(self)
          if(self%n_saved_step>100) then
             total_S=total_S+self%current_average_S
             print *, "average: ", total_S/(self%n_saved_step-100)
          endif
          !self%average_S(:,self%n_saved_step)=
       else
          print *, "Warning: the steps to save exceeded the maximum number of spin in the hist."
       endif
    end if
  end subroutine spin_hist_t_insert

  subroutine spin_hist_t_calc_observables(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_hist_t_calc_observables'
!End of the abilint section

    class(spin_hist_t), intent(inout) :: self
    call spin_hist_t_calc_average_S(self)
  end subroutine spin_hist_t_calc_observables

  subroutine spin_hist_t_calc_average_S(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_hist_t_calc_average_S'
!End of the abilint section

    class(spin_hist_t), intent(inout) :: self
    real(dp) :: a(3), anorm
    a(:)=sum(self%current_S(:,:), dim=2)/self%nmatoms
    anorm=sqrt(sum(a(:)**2))
    self%current_average_S=anorm
    print *, "average S_total, Sx, Sy, Sz: ", anorm,  a(:)

  end subroutine spin_hist_t_calc_average_S


  subroutine spin_hist_t_save_file(self, filename)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_hist_t_save_file'
!End of the abilint section

    class(spin_hist_t), intent(inout) :: self
    character (len=40), intent(in) :: filename
    !TODO
  end subroutine spin_hist_t_save_file

end module m_spin_hist

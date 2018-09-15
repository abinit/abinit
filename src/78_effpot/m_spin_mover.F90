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
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE



#if defined HAVE_CONFIG_H
#include "config.h"
#endif
module m_spin_mover
  use defs_basis
  use m_errors
  use m_abicore
  use m_spin_terms, only: spin_terms_t_get_dSdt, spin_terms_t_get_Langevin_Heff, spin_terms_t_get_gamma_l, spin_terms_t
  !use m_mathfuncs
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_get_s
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_write_one_step
  implicit none
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
  !! nmatoms number of magnetic atoms
  !! SOURCE
  
  
  type spin_mover_t
     integer :: nmatoms
     real(dp) :: dt, total_time, temperature
     !CONTAINS
     !   procedure :: initialize => spin_mover_t_initialize
     !   procedure :: run_one_step => spin_mover_t_run_one_step
     !   procedure :: run_time => spin_mover_t_run_time
  end type spin_mover_t
  !!***

contains

  !!****f* m_spin_mover/spin_mover_t_initialize
  !!
  !! NAME
  !!  spin_mover_t_initialize
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
!!***


  !!****f* m_spin_mover/spin_mover_t_run_one_step
  !!
  !! NAME
  !!  spin_mover_t_run_one_step
  !!
  !! FUNCTION
  !! run one spin step
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! NOTES
  !!
  !! Currently Heun's  (HeunP) integration Method is implemented
  !! should be able to call multi method when implemented
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  
  subroutine spin_mover_t_run_one_step(self, calculator, S_in, S_out, etot)
    !class (spin_mover_t), intent(inout):: self

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_one_step'
!End of the abilint section

    type(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nmatoms)
    real(dp), intent(out) :: S_out(3,self%nmatoms), etot
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

    ! correction
    !call calculator%get_dSdt(S_out, H_lang, dSdt2)
    call spin_terms_t_get_dSdt(calculator, S_out, H_lang, dSdt2)
    etot=calculator%etot
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
  !!***


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
  subroutine spin_mover_t_run_time(self, calculator, hist, ncfile)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_time'
!End of the abilint section

    class(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    type(spin_hist_t), intent(inout) :: hist
    type(spin_ncfile_t), intent(inout) :: ncfile
    real(dp) ::  S(3, self%nmatoms), etot
    real(dp):: t
    integer :: counter
    character(len=100) :: msg
    t=0.0
    counter=0
    write(msg, *) " Begining spin dynamic steps :"
    write(std_out,*) msg
    write(ab_out, *) msg
    write(msg, *)  "==================================================================================" 
    write(std_out,*) msg
    write(ab_out, *) msg

    write(msg, "(A13, 4X, A13, 4X, A13, 4X, A13)")  "Iteration", "time(s)", "average M", "Energy"
    write(std_out,*) msg
    write(ab_out, *) msg
    write(msg, *)  "-----------------------------------------------------------------------------------" 
    write(std_out,*) msg
    write(ab_out, *) msg

    do while(t<self%total_time)
       counter=counter+1
       !call self%run_one_step(calculator, hist%current_S, S)
       call spin_mover_t_run_one_step(self, calculator, spin_hist_t_get_S(hist), S, etot)
       call spin_hist_t_set_vars(hist=hist, S=S, time=t,etot=etot, inc=.True.)
       if(mod(counter, hist%spin_nctime)==0) then
          call spin_ncfile_t_write_one_step(ncfile, hist)
          write(msg, "(I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") counter, t, sqrt(sum((sum(S, dim=2)/self%nmatoms)**2)), &
               & hist%etot(hist%ihist_prev)
          write(std_out,*) msg
          write(ab_out, *) msg
       endif
       !call wrtout_myproc(std_out,msg)
          !print "(I8, 4X, 4ES13.5)", self%current_step, anorm, (a(i) , i=1, 3)
       t=t+self%dt
    enddo
    write(msg, *)  "=====================================================================================" 
    write(std_out,*) msg
    write(ab_out, *) msg

  end subroutine spin_mover_t_run_time
!!***

  
  !!****f* m_spin_mover/spin_mover_t_finalize
  !!
  !! NAME
  !!  spin_mover_t_finalize
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
  subroutine spin_mover_t_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_finalize'
!End of the abilint section

    class(spin_mover_t), intent(inout):: self
  end subroutine spin_mover_t_finalize
  !!***

end module m_spin_mover

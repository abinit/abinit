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
  use m_spin_observables , only : spin_observable_t, ob_calc_observables
  use m_spin_terms, only: spin_terms_t_get_dSdt, spin_terms_t_get_Langevin_Heff, spin_terms_t_get_gamma_l, spin_terms_t
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_get_s, spin_hist_t_reset
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
  !! nspins number of magnetic atoms
  !! SOURCE
  
  
  type spin_mover_t
     integer :: nspins
     real(dp) :: dt, total_time, temperature, pre_time
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
  subroutine spin_mover_t_initialize(self, nspins, dt, total_time, temperature, pre_time)
    !class (spin_mover_t):: self

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_initialize'
!End of the abilint section

    type(spin_mover_t), intent(inout) :: self
    real(dp), intent(in) :: dt, total_time, pre_time,temperature
    integer, intent(in) :: nspins
    self%nspins=nspins
    self%dt=dt
    self%pre_time=pre_time
    self%total_time=total_time
    self%temperature=temperature
  end subroutine spin_mover_t_initialize
!!***



  !!****f* m_spin_mover/spin_mover_t_run_one_step_HeunP
  !!
  !! NAME
  !!  spin_mover_t_run_one_step_HeunP
  !!
  !! FUNCTION
  !! run one spin step using HeunP method
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
  
  subroutine spin_mover_t_run_one_step_HeunP(self, calculator, S_in, S_out, etot)

    !class (spin_mover_t), intent(inout):: self

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_one_step_HeunP'
!End of the abilint section

    type(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nspins)
    real(dp), intent(out) :: S_out(3,self%nspins), etot
    integer :: i
    real(dp) ::  dSdt(3, self%nspins), dSdt2(3, self%nspins), &
         & H_lang(3, self%nspins)
    ! predict
    !call calculator%get_Langevin_Heff(self%dt, self%temperature, H_lang)
    call spin_terms_t_get_Langevin_Heff(calculator, self%dt, self%temperature, H_lang)

    !call calculator%get_dSdt(S_in, H_lang, dSdt)
    call spin_terms_t_get_dSdt(calculator, S_in, H_lang, dSdt)
    !$OMP PARALLEL DO
    do i =1, self%nspins
       S_out(:,i)=  S_in(:,i) +dSdt(:,i) * self%dt
    end do
    !$OMP END PARALLEL DO

    ! correction
    !call calculator%get_dSdt(S_out, H_lang, dSdt2)
    call spin_terms_t_get_dSdt(calculator, S_out, H_lang, dSdt2)
    etot=calculator%etot
    !$OMP PARALLEL DO private(i)
    do i =1, self%nspins
       S_out(:,i)=  S_in(:,i) +(dSdt(:,i)+dSdt2(:,i)) * (0.5_dp*self%dt)
       S_out(:,i)=S_out(:,i)/sqrt(sum(S_out(:,i)**2))
    end do
  end subroutine spin_mover_t_run_one_step_HeunP
  !!***
 
  subroutine spin_mover_t_run_one_step_DM(self, calculator, S_in, S_out, etot)

    ! TODO: implement Depondt & Mertens (2009) method. 
    !class (spin_mover_t), intent(inout):: self

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_one_step_DM'
!End of the abilint section

    type(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nspins)
    real(dp), intent(out) :: S_out(3,self%nspins), etot
    integer :: i
    real(dp) ::  dSdt(3, self%nspins), dSdt2(3, self%nspins), &
         & H_lang(3, self%nspins)
    ! predict
    !call calculator%get_Langevin_Heff(self%dt, self%temperature, H_lang)
    call spin_terms_t_get_Langevin_Heff(calculator, self%dt, self%temperature, H_lang)

    !call calculator%get_dSdt(S_in, H_lang, dSdt)
    call spin_terms_t_get_dSdt(calculator, S_in, H_lang, dSdt)
    !$OMP PARALLEL DO
    do i =1, self%nspins
       S_out(:,i)=  S_in(:,i) +dSdt(:,i) * self%dt
    end do
    !$OMP END PARALLEL DO

    ! correction
    !call calculator%get_dSdt(S_out, H_lang, dSdt2)
    call spin_terms_t_get_dSdt(calculator, S_out, H_lang, dSdt2)
    etot=calculator%etot
    !$OMP PARALLEL DO private(i)
    do i =1, self%nspins
       S_out(:,i)=  S_in(:,i) +(dSdt(:,i)+dSdt2(:,i)) * (0.5_dp*self%dt)
       S_out(:,i)=S_out(:,i)/sqrt(sum(S_out(:,i)**2))
    end do
  end subroutine spin_mover_t_run_one_step_DM
  !!***



  subroutine spin_mover_t_run_one_step(self, calculator, hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_one_step'
!End of the abilint section

    type(spin_mover_t), intent(inout) :: self
    type(spin_terms_t), intent(inout) :: calculator
    type(spin_hist_t),intent(inout) :: hist
    real(dp) :: S_out(3,self%nspins), etot
    call spin_mover_t_run_one_step_HeunP(self, calculator, spin_hist_t_get_S(hist), S_out, etot)
    ! do not inc until time is set to hist.
    call spin_hist_t_set_vars(hist=hist, S=S_out, Snorm=calculator%ms, etot=etot, inc=.False.)
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
  subroutine spin_mover_t_run_time(self, calculator, hist, ncfile, ob)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_mover_t_run_time'
!End of the abilint section

    class(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    type(spin_hist_t), intent(inout) :: hist
    type(spin_ncfile_t), intent(inout) :: ncfile
    type(spin_observable_t), intent(inout) :: ob
    !real(dp) ::  S(3, self%nspins)
    real(dp):: t
    integer :: counter, i, ii
    character(len=80) :: msg
    t=0.0
    counter=0
    write(msg, *) " Begining spin dynamic steps :"

     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')
    msg=repeat("=", 65)

     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

    write(msg, "(A13, 4X, A13, 4X, A13, 4X, A13)")  "Iteration", "time(s)", "Avg_Mst/Ms", "Energy (Ha)"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

    msg=repeat("-", 65)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

    if (abs(self%pre_time) > 1e-30) then
       msg="Pre-run:"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

 
       do while(t<self%pre_time)
          counter=counter+1
          call spin_mover_t_run_one_step(self, calculator, hist)
          call spin_hist_t_set_vars(hist=hist, time=t,  inc=.True.)
          if(mod(counter, hist%spin_nctime)==0) then
             call ob_calc_observables(ob, hist%S(:,:, hist%ihist_prev), hist%Snorm(:,hist%ihist_prev))
             write(msg, "(I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") counter, t, &
                  & ob%Mst_norm_total/ob%Snorm_total, &
                  & hist%etot(hist%ihist_prev)/Ha_J
          ! total : 13+4+...= 64 
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

 
          endif
          t=t+self%dt
       end do
       t=0.0
       counter=0
       call spin_hist_t_reset(hist,array_to_zero=.False.)
       msg="Formal run:"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')


    endif

    do while(t<self%total_time)
       counter=counter+1
       call spin_mover_t_run_one_step(self, calculator, hist)
       call spin_hist_t_set_vars(hist=hist, time=t,  inc=.True.)
       if(mod(counter, hist%spin_nctime)==0) then
          call ob_calc_observables(ob, hist%S(:,:, hist%ihist_prev), hist%Snorm(:,hist%ihist_prev))
          call spin_ncfile_t_write_one_step(ncfile, hist)
          write(msg, "(I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") counter, t, &
               & ob%Mst_norm_total/ob%Snorm_total, &
               & hist%etot(hist%ihist_prev)/Ha_J
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')


          write(std_out,*) msg
          write(ab_out, *) msg
       endif
       t=t+self%dt
    enddo

    msg=repeat("-", 65)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')



    write(msg, "(A30)") "Summary of spin dynamics:"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')

    write(msg, "(A65)") "At the end of the run, the average spin at each sublattice is"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')


    write(msg, "(8X, 10X, A5, A2, 4X, 4A10)")  'ID', ": ", '<M_i>(x)', '<M_i>(y)', '<M_i>(z)', '||<M_i>||'
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')


    do i =1, ob%nsublatt
       write(msg, "(8X, A10, I5.4, A2, 4X, 4F10.5)") "Sublattice", i, ": ", (ob%Mst_sub(ii,i)/ob%nspins_sub(i)/mu_B_SI , ii=1, 3), &
            sqrt(sum((ob%Mst_sub(:, i)/ob%nspins_sub(i)/mu_B_SI)**2))
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')
    end do

    msg=repeat("=", 65)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out, msg, 'COLL')



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

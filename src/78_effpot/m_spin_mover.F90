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
  use m_errors
  use m_abicore
  use m_spin_observables , only : spin_observable_t, ob_calc_observables, ob_reset
  use m_spin_terms, only: spin_terms_t_get_dSdt, spin_terms_t_get_Langevin_Heff, &
       & spin_terms_t_get_gamma_l, spin_terms_t, spin_terms_t_total_Heff, spin_terms_t_get_etot, &
       & spin_terms_t_Hrotate
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
     integer :: nspins, method
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
  subroutine spin_mover_t_initialize(self, nspins, dt, total_time, temperature, pre_time, method)
    !class (spin_mover_t):: self
    type(spin_mover_t), intent(inout) :: self
    real(dp), intent(in) :: dt, total_time, pre_time,temperature
    integer, intent(in) :: nspins, method
    self%nspins=nspins
    self%dt=dt
    self%pre_time=pre_time
    self%total_time=total_time
    self%temperature=temperature
    self%method=method
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

  subroutine rotate_S_DM(S_in, Heff, dt, S_out)
    ! Depondt & Mertens method to roate S_in
    real(dp), intent(in) :: S_in(3), Heff(3), dt
    real(dp), intent(inout) :: S_out(3)
    real(dp) :: B(3) , w, u, Bnorm, R(3,3), cosw, sinw
    Bnorm=sqrt(sum(Heff*Heff))
    B(:)=Heff(:)/Bnorm
    w=Bnorm*dt
    !print *, w
    cosw=cos(w)
    sinw=1.0d0-cosw*cosw
    u=1.0d0-cosw
    R(1,1)=B(1)*B(1)*u+cosw
    R(2,1)=B(1)*B(2)*u+B(3)*sinw
    R(3,1)=B(1)*B(3)*u-B(2)*sinw

    R(1,2)=B(1)*B(2)*u-B(3)*sinw
    R(2,2)=B(2)*B(2)*u+cosw
    R(3,2)=B(2)*B(3)*u+B(1)*sinw

    R(1,3)=B(1)*B(3)*u+B(2)*sinw
    R(2,3)=B(2)*B(3)*u-B(1)*sinw
    R(3,3)=B(3)*B(3)*u+cosw
    S_out=matmul(R, S_in)
  end subroutine rotate_S_DM

  subroutine spin_mover_t_run_one_step_DM(self, calculator, S_in, S_out, etot)
    ! Depondt & Mertens (2009) method, using a rotation matrix so length doesn't change.
    !class (spin_mover_t), intent(inout):: self
    type(spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nspins)
    real(dp), intent(out) :: S_out(3,self%nspins), etot
    integer :: i
    real(dp) :: H_lang(3, self%nspins),  Heff(3, self%nspins), Htmp(3, self%nspins), Hrotate(3, self%nspins)
    ! predict
    call spin_terms_t_total_Heff(calculator, S=S_in, Heff=Heff)

    call spin_terms_t_get_Langevin_Heff(calculator, self%dt, self%temperature, H_lang)
    Htmp(:,:)=Heff(:,:)+H_lang(:,:)

    call spin_terms_t_Hrotate(calculator, Htmp, S_in, Hrotate)
    do i=1, self%nspins
     call rotate_S_DM(S_in(:,i), Htmp(:,i), self%dt, S_out(:,i))
    end do

    ! correction
    call spin_terms_t_total_Heff(self=calculator, S=S_in, Heff=Htmp)
    Heff(:,:)=(Heff(:,:)+Htmp(:,:))*0.5_dp
    Htmp(:,:)=Heff(:,:)+H_lang(:,:)
    call spin_terms_t_Hrotate(calculator, Htmp, S_in, Hrotate)

    do i=1, self%nspins
       call rotate_S_DM(S_in(:,i), Hrotate(:,i), self%dt, S_out(:,i))
    end do

    call spin_terms_t_get_etot(calculator, S_out, Heff, etot)
  end subroutine spin_mover_t_run_one_step_DM


  subroutine spin_mover_t_run_one_step(self, calculator, hist)

    type(spin_mover_t), intent(inout) :: self
    type(spin_terms_t), intent(inout) :: calculator
    type(spin_hist_t),intent(inout) :: hist
    real(dp) :: S_out(3,self%nspins), etot
    if(self%method==1) then
       call spin_mover_t_run_one_step_HeunP(self, calculator, spin_hist_t_get_S(hist), S_out, etot)
    else if (self%method==2) then
       call spin_mover_t_run_one_step_DM(self, calculator, spin_hist_t_get_S(hist), S_out, etot)
    end if

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
    write(msg, *) " Beginning spin dynamic steps :"

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
       msg="Thermolization run:"
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')

       do while(t<self%pre_time)
          counter=counter+1
          call spin_mover_t_run_one_step(self, calculator, hist)
          call spin_hist_t_set_vars(hist=hist, time=t,  inc=.True.)
          if(mod(counter, hist%spin_nctime)==0) then
             call ob_calc_observables(ob, hist%S(:,:, hist%ihist_prev), &
                  hist%Snorm(:,hist%ihist_prev), hist%etot(hist%ihist_prev))
             write(msg, "(A1, 1X, I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") "-", counter, t, &
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
       msg="Measurement run:"
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')


    endif
    call ob_reset(ob)

    do while(t<self%total_time)
       counter=counter+1
       call spin_mover_t_run_one_step(self, calculator, hist)
       call spin_hist_t_set_vars(hist=hist, time=t,  inc=.True.)
       call ob_calc_observables(ob, hist%S(:,:, hist%ihist_prev), &
            hist%Snorm(:,hist%ihist_prev), hist%etot(hist%ihist_prev))
       if(mod(counter, hist%spin_nctime)==0) then
          call spin_ncfile_t_write_one_step(ncfile, hist)
          write(msg, "(A1, 1X, I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") "-", counter, t, &
               & ob%Mst_norm_total/ob%Snorm_total, &
               & hist%etot(hist%ihist_prev)/Ha_J
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
       endif
       t=t+self%dt
    enddo

    msg=repeat("-", 65)
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

    do i =1, ob%nsublatt
       write(msg, "(A1, 5X, 2X, I5.4, 8X, 4F10.5)") '-', i, (ob%Mst_sub(ii,i)/ob%nspins_sub(i)/mu_B_SI , ii=1, 3), &
            sqrt(sum((ob%Mst_sub(:, i)/ob%nspins_sub(i)/mu_B_SI)**2))
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out, msg, 'COLL')
    end do
    msg=repeat("-", 65)
    call wrtout(std_out,msg,'COLL')
    call wrtout(ab_out, msg, 'COLL')


    write(msg, "(A1, 1X, A11, 3X, A13, 3X, A13, 3X, A13, 3X, A13 )" ) &
         "#", "Temperature", "Cv", "chi",  "BinderU4", "Mst"
    call wrtout(std_out, msg, "COLL")
    call wrtout(ab_out, msg, "COLL")
    write(msg, "(2X, F11.5, 3X, ES13.5, 3X, ES13.5, 3X, E13.5, 3X, ES13.5, 3X )" ) &
         self%temperature , ob%Cv, ob%chi,  ob%binderU4, ob%Avg_Mst_norm_total/ob%snorm_total
    call wrtout(std_out, msg, "COLL")
    call wrtout(ab_out,  msg, "COLL")

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

    class(spin_mover_t), intent(inout):: self
    ABI_UNUSED((/self%nspins/))
  end subroutine spin_mover_t_finalize
  !!***

end module m_spin_mover

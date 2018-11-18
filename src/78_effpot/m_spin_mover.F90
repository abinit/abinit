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
#include "abi_common.h"
module m_spin_mover
  use defs_basis
  use m_errors
  use m_abicore
  use m_mathfuncs, only : cross
  use m_spin_observables , only : spin_observable_t, ob_calc_observables, ob_reset
  use m_spin_terms, only:  spin_terms_t, spin_terms_t_total_Heff, spin_terms_t_get_etot 
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_get_s, spin_hist_t_reset
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_write_one_step
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
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
     real(dp), allocatable :: gyro_ratio(:), damping(:), gamma_L(:), H_lang_coeff(:), ms(:)
     type(rng_t) :: rng
     logical :: gamma_l_calculated
     CONTAINS
        procedure :: initialize => spin_mover_t_initialize
        procedure :: finalize => spin_mover_t_finalize
        procedure :: run_one_step => spin_mover_t_run_one_step
        procedure :: run_time => spin_mover_t_run_time
        procedure :: set_Langevin_params
        procedure :: get_Langevin_Heff
        procedure :: get_dSdt
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
  subroutine spin_mover_t_initialize(self, params, nspins)
    class(spin_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    integer, intent(in) :: nspins
    self%nspins=nspins
    self%dt= params%spin_dt
    self%pre_time= params%spin_ntime_pre * self%dt
    self%total_time= params%spin_ntime * self%dt
    self%temperature=params%spin_temperature
    if(params%spin_dynamics>=0) then
       self%method=params%spin_dynamics
    endif
    call set_seed(self%rng, [111111_dp, 2_dp])


    ABI_ALLOCATE(self%ms, (nspins) )
    ABI_ALLOCATE(self%gyro_ratio, (nspins) )
    ABI_ALLOCATE(self%damping, (nspins) )
    ABI_ALLOCATE(self%gamma_l, (nspins) )
    ABI_ALLOCATE(self%H_lang_coeff, (nspins) )

    self%gamma_l_calculated=.False.
  end subroutine spin_mover_t_initialize
  !!***


  subroutine set_Langevin_params(self, gyro_ratio, damping, temperature, ms)
    class(spin_mover_t), intent(inout) :: self
    real(dp), optional, intent(in) :: damping(self%nspins), temperature, &
         &    gyro_ratio(self%nspins), ms(self%nspins)

    if(present(damping)) then
       self%damping(:)=damping(:)
    endif

    if(present(gyro_ratio)) then
       self%gyro_ratio(:)=gyro_ratio(:)
    endif
    if(present(ms)) then
       self%ms(:)=ms(:)
    endif


    if(present(temperature)) then
       self%temperature=temperature
    end if


    self%gamma_l(:)= self%gyro_ratio(:)/(1.0_dp+ self%damping(:)**2)
    self%gamma_l_calculated=.True.

    self%H_lang_coeff(:)=sqrt(2.0*self%damping(:)* self%temperature &
         &  /(self%gyro_ratio(:)* self%dt *self%ms(:)))
  end subroutine set_Langevin_params


  subroutine get_Langevin_Heff(self, H_lang)
    class(spin_mover_t), intent(inout) :: self
    real(dp), intent(inout):: H_lang(3,self%nspins)
    integer :: i
      if ( self%temperature .gt. 1d-7) then
         call rand_normal_array(self%rng, H_lang, 3*self%nspins)
         do i = 1, self%nspins
            H_lang(:,i)= H_lang(:,i) * self%H_lang_coeff(i)
         end do
      else
         H_lang(:,:)=0.0_dp
      end if
    end subroutine get_Langevin_Heff

  subroutine get_dSdt(self, S, Heff, dSdt)
    class(spin_mover_t), intent(inout) :: self
    real(dp), intent(in) :: Heff(3,self%nspins), S(3,self%nspins)
    real(dp), intent(out) :: dSdt(3, self%nspins)
    integer :: i
    real(dp) :: Ri(3)
    !$OMP PARALLEL DO private(Ri, i)
    do i=1,self%nspins
       Ri = cross(S(:,i),Heff(:,i))
       dSdt(:,i) = -self%gamma_L(i)*(Ri+self%damping(i)* cross(S(:,i), Ri))
    end do
    !$OMP END PARALLEL DO
  end subroutine get_dSdt



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
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE

  subroutine spin_mover_t_run_one_step_HeunP(self, calculator, S_in, S_out, etot)

    class (spin_mover_t), intent(inout):: self
    type(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nspins)
    real(dp), intent(out) :: S_out(3,self%nspins), etot
    integer :: i
    real(dp) ::  dSdt(3, self%nspins), dSdt2(3, self%nspins), &
         & H_lang(3, self%nspins), Htmp(3, self%nspins), Heff(3, self%nspins) 
    ! predict
    call self%get_Langevin_Heff(H_lang)
    call spin_terms_t_total_Heff(calculator, S=S_in, Heff=Heff)
    Htmp(:,:)=Heff(:,:)+H_lang(:,:)

    call self%get_dSdt(S_in, Heff, dSdt)

    !$OMP PARALLEL DO private(i)
    do i =1, self%nspins
       S_out(:,i)=  S_in(:,i) +dSdt(:,i) * self%dt
    end do
    !$OMP END PARALLEL DO

    ! correction
    !call calculator%get_dSdt(S_out, H_lang, dSdt2)
    call spin_terms_t_total_Heff(calculator, S=S_in, Heff=Heff)
    etot=calculator%etot
    !$OMP PARALLEL DO private(i)
    do i =1, self%nspins
       S_out(:,i)=  S_in(:,i) +(dSdt(:,i)+dSdt2(:,i)) * (0.5_dp*self%dt)
       S_out(:,i)=S_out(:,i)/sqrt(sum(S_out(:,i)**2))
    end do
    !$OMP END PARALLEL DO
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
    sinw=sin(w)
    cosw=sqrt(1.0_dp- sinw*sinw)
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
    class(spin_mover_t), intent(inout):: self
    class(spin_terms_t), intent(inout) :: calculator
    real(dp), intent(in) :: S_in(3,self%nspins)
    real(dp), intent(out) :: S_out(3,self%nspins), etot
    integer :: i
    real(dp) :: H_lang(3, self%nspins),  Heff(3, self%nspins), Htmp(3, self%nspins), Hrotate(3, self%nspins)
    ! predict
    call spin_terms_t_total_Heff(calculator, S=S_in, Heff=Heff)

    call self%get_Langevin_Heff(H_lang)
    Htmp(:,:)=Heff(:,:)+H_lang(:,:)

!$OMP PARALLEL DO private(i)
    do i=1,self%nspins
       ! Note that there is no - , because dsdt =-cross (S, Hrotate) 
       Hrotate(:,i) = self%gamma_L(i) * ( Htmp(:,i) + self%damping(i)* cross(S_in(:,i), Htmp(:,i)))
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO 
    do i=1, self%nspins
     call rotate_S_DM(S_in(:,i), Hrotate(:,i), self%dt, S_out(:,i))
    end do
!$OMP END PARALLEL DO

    ! correction
    call spin_terms_t_total_Heff(self=calculator, S=S_in, Heff=Htmp)
    Heff(:,:)=(Heff(:,:)+Htmp(:,:))*0.5_dp
    Htmp(:,:)=Heff(:,:)+H_lang(:,:)
    !call spin_terms_t_Hrotate(calculator, Htmp, S_in, Hrotate)
    !$OMP PARALLEL DO private(i)
    do i=1,self%nspins
       Hrotate(:,i) = self%gamma_L(i) * ( Htmp(:,i) + self%damping(i)* cross(S_in(:,i), Htmp(:,i)))
    end do
    !$OMP END PARALLEL DO

    do i=1, self%nspins
       call rotate_S_DM(S_in(:,i), Hrotate(:,i), self%dt, S_out(:,i))
    end do

    call spin_terms_t_get_etot(calculator, S_out, Heff, etot)
  end subroutine spin_mover_t_run_one_step_DM


  subroutine spin_mover_t_run_one_step(self, calculator, hist)
    class(spin_mover_t), intent(inout) :: self
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
    character(len=80) :: msg, msg_empty
    t=0.0
    counter=0

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
                  & hist%etot(hist%ihist_prev)/ob%nscell
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
       call self%run_one_step(calculator, hist)
       call spin_hist_t_set_vars(hist=hist, time=t,  inc=.True.)
       call ob_calc_observables(ob, hist%S(:,:, hist%ihist_prev), &
            hist%Snorm(:,hist%ihist_prev), hist%etot(hist%ihist_prev))
       if(mod(counter, hist%spin_nctime)==0) then
          call spin_ncfile_t_write_one_step(ncfile, hist)
          write(msg, "(A1, 1X, I13, 4X, ES13.5, 4X, ES13.5, 4X, ES13.5)") "-", counter, t*Time_Sec, &
               & ob%Mst_norm_total/ob%Snorm_total, &
               & hist%etot(hist%ihist_prev)/ob%nscell
          call wrtout(std_out,msg,'COLL')
          call wrtout(ab_out, msg, 'COLL')
       endif
       t=t+self%dt
    enddo

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

    do i =1, ob%nsublatt
       write(msg, "(A1, 5X, 2X, I5.4, 8X, 4F10.5)") '-', i, (ob%Mst_sub(ii,i)/ob%nspins_sub(i)/mu_B , ii=1, 3), &
            sqrt(sum((ob%Mst_sub(:, i)/ob%nspins_sub(i)/mu_B)**2))
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
         self%temperature*Ha_K , ob%Cv, ob%chi,  ob%binderU4, ob%Avg_Mst_norm_total/ob%snorm_total
    call wrtout(std_out, msg, "COLL")
    call wrtout(ab_out,  msg, "COLL")

    msg=repeat("=", 80)
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




  end subroutine spin_mover_t_finalize
  !!***

end module m_spin_mover

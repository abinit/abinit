#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"
    module m_spin_mc_mover
    use defs_basis
    use m_errors
    use m_abicore
    use m_mathfuncs, only : cross
    use m_spin_observables , only : spin_observable_t, ob_calc_observables, ob_reset
    use m_spin_terms, only:  spin_terms_t
    use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_get_s, spin_hist_t_reset
    use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_write_one_step
    use m_multibinit_dataset, only: multibinit_dtset_type
    use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
    use m_effpot_api, only: effpot_t
    use m_mover_api, only: abstract_mover_t
    implicit none
    private

    type, extends(spin_mover_t):: spin_mc_mover_t
       type(rng_t) :: rng
     contains
    end type spin_mc_mover_t

  contains
    subroutine initialize(self)
    end subroutine initialize

   subroutine run_one_step(self, effpot)
     class(spin_mc_mover_t) :: self
     class(effpot_t), intent(inout) :: effpot
     real(dp) :: S_out(3,self%nspins), etot
     integer :: i, j
     if(self%rng%rand_unif_01< min(1.0, self%attempt())) then
        self%naccept=self%naccept+1
        call self%accept()
     else
        call self%reject()
     end if
   end subroutine run_one_step

   function attempt(self) result(y)
     class(spin_mc_mover_t) :: self
     real(dp) :: y, Sold(3), Snew(3)
     integer :: ispin
     ! choose one site
     ispin = self%rand_choice(self%nspins)
     Sold(3)= self%S(:,ispin)
     call move_hinzke_nowak(Sold, Snew, self%angle)
     call self%effpot%calc_delta_energy(ispin, S, Sold, Snew)
   end function attempt

   subroutine move_angle(Sold, Snew, angle)
     real(dp), intent(in) :: Sold(3), angle
     real(dp), intent(out) :: Snew(3)
     real(dp) :: m
     call rng%rand_normal_array(Snew, 3)
     Snew(:)=Sold(:) + Snew(:)*angle
     Snew(:)=Snew(:)/sqrt(Snew(1)*Snew(1)+Snew(2)*Snew(2)+Snew(3)*Snew(3))
   end subroutine move_angle

   subroutine move_flip(Sold, Snew)
     real(dp), intent(in) :: Sold(3)
     real(dp), intent(out) :: Snew(3)
     Snew(:)=-Sold(:)
   end subroutine move_flip

   subroutine move_uniform( Snew)
     real(dp), intent(out) :: Snew(3)
     call rng%rand_normal_array(Snew, 3)
     Snew(:)=Snew(:)/sqrt(Snew(1)*Snew(1)+Snew(2)*Snew(2)+Snew(3)*Snew(3))
   end subroutine move_uniform

   subroutine move_hinzke_nowak(Sold, Snew, angle)
     real(dp), intent(in) :: Sold(3), angle
     real(dp), intent(out) :: Snew(3)
     integer :: move
     move=rng%rand_choice(3)
     select case (move)
     case (1)
        call move_angle(Sold, Snew, angle)
     case(2)
        call move_flip(Sold, Snew)
     case(3)
        call move_uniform(Snew)
     case default
        call move_angle(Sold, Snew, angle)
     end select
   end subroutine move_hinzke_nowak

   subroutine accept(self)
   end subroutine accept

   subroutine reject(self)
   end subroutine reject
end module m_spin_mc_mover

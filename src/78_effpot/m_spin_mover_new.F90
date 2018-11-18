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
     CONTAINS
        procedure :: initialize => spin_mover_t_initialize
        procedure :: run_one_step => spin_mover_t_run_one_step
        procedure :: run_time => spin_mover_t_run_time

        procedure :: get_Langevin_Heff
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
  subroutine spin_mover_t_initialize(self, params)
    class(spin_mover_t), intent(inout) :: self
    type(m_multibinit_dataset), target, intent(in):: params
    self%params => params

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



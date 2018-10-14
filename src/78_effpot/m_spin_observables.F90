!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_model
!! NAME
!! m_spin_observables
!!
!! FUNCTION
!! This module contains the subroutines to calculate the observables of spin dynamics
!!
!!
!! Datatypes:
!!
!! Subroutines:
!!
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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

module m_spin_observables

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_spin_model_primitive, only: spin_model_primitive_t
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_init, spin_hist_t_get_s, spin_hist_t_free, &
       & spin_hist_t_set_params
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_write_one_step
  implicit none
  !!***
  public :: spin_hist_t_ob_calc_staggered_M
  public :: spin_hist_t_ob_calc_thermo_obs
  public :: spin_hist_t_ob_calc_correlation_obs
  public :: spin_hist_t_ob_calc_traj_obs

contains

  subroutine spin_hist_t_ob_calc_staggered_M(self)
    class(spin_hist_t) :: self
  end subroutine spin_hist_t_ob_calc_staggered_M

  subroutine spin_hist_t_ob_calc_traj_obs(self)
      class(spin_hist_t) :: self
  end subroutine spin_hist_t_ob_calc_traj_obs

  subroutine spin_hist_t_ob_calc_thermo_obs(self)
    class(spin_hist_t) :: self
  end subroutine spin_hist_t_ob_calc_thermo_obs

  subroutine spin_hist_t_ob_calc_correlation_obs(self)
    class(spin_hist_t) :: self
  end subroutine spin_hist_t_ob_calc_correlation_obs


end module m_spin_observables

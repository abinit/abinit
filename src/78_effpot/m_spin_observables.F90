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

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_spin_terms, only: spin_terms_t, spin_terms_t_finalize, spin_terms_t_set_external_hfield
  use m_spin_model_primitive, only: spin_model_primitive_t
  use m_spin_hist, only: spin_hist_t, spin_hist_t_set_vars, spin_hist_t_init, spin_hist_t_get_s, spin_hist_t_free, &
       & spin_hist_t_set_params
  use m_spin_mover, only: spin_mover_t, spin_mover_t_initialize, spin_mover_t_finalize, &
       & spin_mover_t_run_time, spin_mover_t_run_one_step
  use m_spin_ncfile, only: spin_ncfile_t, spin_ncfile_t_init, spin_ncfile_t_close, spin_ncfile_t_def_sd, &
       & spin_ncfile_t_write_primitive_cell, spin_ncfile_t_write_supercell, spin_ncfile_t_write_parameters, &
       & spin_ncfile_t_write_one_step
  implicit none
!!***


end module m_spin_observables

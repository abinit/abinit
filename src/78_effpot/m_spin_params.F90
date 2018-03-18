#include "abi_common.h"
module m_spin_params
  use defs_basis
  implicit none
  type spin_params_t
     integer :: nmatoms=0
     integer :: hist_max_save=1000
     integer :: hist_step_save=100
     real(dp) :: dtspin=1d-14
     real(dp) :: damping_factor=1.0
     real(dp) :: total_time=1e-11
  end type spin_params_t
end module m_spin_params

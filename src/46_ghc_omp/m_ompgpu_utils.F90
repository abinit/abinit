!!****m* ABINIT/m_ompgpu_utils
!! NAME
!!  m_ompgpu_utils
!!
!! FUNCTION
!!  Collection of routines useful for leveraging OpenMP GPU offload capabilities.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2025 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ompgpu_utils

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 implicit none

 private


#ifdef HAVE_OPENMP_OFFLOAD
 ! pointer to GPU buffers
 integer, pointer, save :: cur_kg_k(:,:) => null()
 integer, pointer, save :: cur_kg_kp(:,:) => null()
 real(dp), pointer, save :: cur_ffnl_k(:,:,:,:) => null()
 real(dp), pointer, save :: cur_ph3d_k(:,:,:) => null()
 integer, pointer, save :: cur_kg_k_gather(:,:) => null()
#endif

 public :: ompgpu_load_hamilt_buffers
 public :: ompgpu_free_hamilt_buffers

contains

 subroutine ompgpu_load_hamilt_buffers(kg_k,kg_kp,ffnl_k,ph3d_k,kg_k_gather)
   integer,intent(in),target :: kg_k(:,:), kg_kp(:,:)
   real(dp),intent(in),target :: ffnl_k(:,:,:,:),ph3d_k(:,:,:)
   integer,intent(in),target,optional :: kg_k_gather(:,:)

#ifdef HAVE_OPENMP_OFFLOAD
   if(associated(cur_kg_k)) call ompgpu_free_hamilt_buffers()

   cur_kg_k            => kg_k
   cur_kg_kp           => kg_kp
   cur_ffnl_k          => ffnl_k
   cur_ph3d_k          => ph3d_k
   if(present(kg_k_gather))  cur_kg_k_gather     => kg_k_gather
   !$OMP TARGET ENTER DATA MAP(to:cur_kg_k,cur_kg_kp,cur_ffnl_k,cur_ph3d_k,cur_kg_k_gather)
#else
   ABI_UNUSED((/kg_k,kg_kp,kg_k_gather/))
   ABI_UNUSED((/ffnl_k,ph3d_k/))
   ABI_BUG("ABINIT wasn't compiled with OpenMP GPU offloading, aborting.")
#endif
 end subroutine ompgpu_load_hamilt_buffers

 subroutine ompgpu_free_hamilt_buffers()
#ifdef HAVE_OPENMP_OFFLOAD

   !$OMP TARGET EXIT DATA MAP(release:cur_kg_k,cur_kg_kp,cur_ffnl_k,cur_ph3d_k,cur_kg_k_gather)
   cur_kg_k            => null()
   cur_kg_kp           => null()
   cur_kg_k_gather     => null()
#else
   ABI_BUG("ABINIT wasn't compiled with OpenMP GPU offloading, aborting.")
#endif
 end subroutine ompgpu_free_hamilt_buffers

end module m_ompgpu_utils
!!***


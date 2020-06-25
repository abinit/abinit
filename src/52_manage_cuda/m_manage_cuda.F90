!!****m* ABINIT/m_manage_cuda
!! NAME
!!  m_manage_cuda
!!
!! FUNCTION
!!  Fake module to dupe the build system and allow it to include cuda files
!!   in the chain of dependencies.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2020 ABINIT group (MT)
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

module m_manage_cuda

 use m_alloc_hamilt_gpu

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#endif

 implicit none

!Interfaces for C bindings --- To be completed
#ifdef HAVE_FC_ISO_C_BINDING
#if defined HAVE_GPU_CUDA
!interface
!  integer(C_INT) function cuda_func() bind(C)
!    use iso_c_binding, only : C_INT,C_PTR
!    type(C_PTR) :: ptr
!  end function cuda_func
#endif
#endif

contains
!!***

end module m_manage_cuda
!!***

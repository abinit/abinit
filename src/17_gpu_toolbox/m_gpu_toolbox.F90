!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gpu_toolbox
!! NAME
!!  m_gpu_toolbox
!!
!! FUNCTION
!!  Fake module to dupe the build system and allow it to include cuda files
!!   in the chain of dependencies.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2019 ABINIT group (MT)
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

module m_gpu_toolbox

 use m_initcuda
 use m_gpu_detect

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

end module m_gpu_toolbox
!!***

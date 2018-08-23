!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abicore
!! NAME
!!  m_abicore
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

module m_abicore

 !public

 !use defs_basis
 !use m_build_info
 use m_profiling_abi
 use m_specialmsg,  only : herald, specialmsg_setcount, specialmsg_getcount, specialmsg_mpisum, wrtout
 !use m_errors

 !implicit none
 !private
!!***

end module m_abicore
!!***

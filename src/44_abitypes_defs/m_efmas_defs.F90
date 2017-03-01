!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_efmas_defs
!! NAME
!! m_efmas_defs
!!
!! FUNCTION
!! This module contains datatypes for efmas functionalities.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_efmas_defs

 use defs_basis

 implicit none
 public
 !Put eventual module variables here...
!!***

!!****t* m_efmas_defs/efmasfr_type
!! NAME
!! efmasfr_type
!!
!! FUNCTION
!! The efmasfr_type structured datatype
!!
!! SOURCE

 type efmasfr_type
   
   !For k-point
   complex(dpc),allocatable :: ch2c(:,:,:,:)

 end type efmasfr_type
!!***

!!****t* m_efmas_defs/efmasdeg_type
!! NAME
!! efmasdeg_type
!!
!! FUNCTION
!! The efmasdeg_type structured datatype
!!
!! SOURCE

 type efmasdeg_type
   
   !For k-point
   integer :: ndegs
   integer, allocatable :: degs_bounds(:,:)
   !For band
   logical,allocatable :: degenerate(:), treated(:)
   integer :: band_range(2), deg_range(2)
   integer,allocatable :: deg_dim(:), degl(:), ideg(:)

 end type efmasdeg_type
!!***

end module m_efmas_defs
!!***

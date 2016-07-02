!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_parameters
!! NAME
!! defs_parameters
!!
!! FUNCTION
!! Default parameters.
!! Typical use : parameters that are used in different routines, of which
!! the value is fixed here (still, they are not as fundamental as the
!! definitions in defs_basis).
!! Also, some input variables, of integer type, might be associated
!! with a name for selected values. One might associated here these
!! these particular values with a variable name.
!! Please, make sure that these parameters are easy to trace, in the
!! different routines they are used. Hence, give them a characteristic,
!! informative, name.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG,YMN)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_parameters

 use defs_basis

 implicit none

!- Set of parameters for the aim utility -----------------------------------
 real(dp), parameter :: aim_rhocormin=1.d-10  ! the minimal core density
 real(dp), parameter :: aim_epstep=0.5
 real(dp), parameter :: aim_rhomin=1.d-5,aim_dgmin=1.d-9,aim_dmaxcrit=5.d-2
 real(dp), parameter :: aim_dmin=1.d-3,aim_hmax=2.d7,aim_fac0=2.1_dp,aim_facmin=1.d-3
 real(dp), parameter :: aim_hmult=15._dp,aim_tiny=1.d-4,aim_snull=1.d-6
 real(dp), parameter :: aim_deltarmin=1.d-7
!the minimal length of one step following the gradient line
 real(dp), parameter :: aim_fac=1.2_dp,aim_drmin=1.d-5
 real(dp), parameter :: aim_dlimit=1.d-4,aim_dmaxcs=3.d-1
 real(dp), parameter :: aim_dpc0=1.d-2
 integer, parameter :: aim_maxstep=100
 real(dp), parameter :: aim_xymin=1.d-10
 integer, parameter :: aim_npmaxin=17
 real(dp), parameter :: aim_stmax=0.05
 real(dp), parameter :: aim_dmaxc1=1.d-1, aim_dmaxcl=5.d-2

!- Particular values of some input variables -----------------------------------

end module defs_parameters

!!***

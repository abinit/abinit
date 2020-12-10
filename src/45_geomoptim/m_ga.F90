!!****m* ABINIT/m_ga
!! NAME
!!  m_ga
!!
!! FUNCTION
!!  This module provides several routines and datatypes for the
!!  Genetic algorithm stochastic search implementation.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2020 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ga

 use defs_basis
 use m_errors
 use m_abicore
 use m_dtset

 implicit none

 private

!public procedures
 public :: ga_init
 public :: ga_destroy
 public :: ga_nullify
!!***

!!****t* m_ga/ga_type
!! NAME
!! ga_type
!!
!! FUNCTION
!! Datatype with the variables required to perform GA random search
!!
!! NOTES
!!
!! SOURCE

 type,public :: ga_type
! Scalars
  integer  :: ga_start     ! Starting iteration for the Genetic algorithm
  integer  :: ga_n_rules   ! Number of genetic rules chosen by the user
  integer  :: ga_fitness   ! Fitness function chosen by the user
  integer  :: ga_algor     ! Genetic algorithm
  real(dp) :: ga_opt_percent   ! Percentage of the population that passes from one iteration to the next
! Arrays same dimension as defined for the whole data input data structure
  integer  :: ga_rules(30)         ! Genetic rules chosen by the user
  integer,allocatable  :: ga_iatfix(:,:)  ! defined fixed atoms and directions
 end type ga_type
!!***

CONTAINS !===========================================================
!!***

!!****f* m_ga/ga_init
!! NAME
!!  ga_init
!!
!! FUNCTION
!!  Initialize a datastructure of type ga_type.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in current dataset
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ga_param=datastructure of type ga_type.
!!            several parameters for Genetic Algorithm random search.
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine ga_init(dtset,ga_param)

!Arguments ------------------------------------
!scalars
 type(dataset_type),target,intent(in) :: dtset
 type(ga_type),intent(inout) :: ga_param

!************************************************************************

 if(dtset%imgmov==4)then
   ga_param%ga_n_rules           = dtset%ga_n_rules
   ga_param%ga_rules(:)          = dtset%ga_rules
   ga_param%ga_fitness           = dtset%ga_fitness
   ga_param%ga_opt_percent       = dtset%ga_opt_percent
   ga_param%ga_algor             = dtset%ga_algor
   ABI_MALLOC(ga_param%ga_iatfix,(3,dtset%natom))
   ga_param%ga_iatfix            = dtset%iatfix
 end if

end subroutine ga_init
!!***

!----------------------------------------------------------------------

!!****f* m_ga/ga_destroy
!! NAME
!!  ga_destroy
!!
!! FUNCTION
!!  Destroy the content of a datastructure of type ga_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ga_param=datastructure of type ga_type.
!!            parameters for Genetic algorithm search.
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine ga_destroy(ga_param)

!Arguments ------------------------------------
!scalars
 type(ga_type),intent(inout) :: ga_param

!************************************************************************

 if(allocated(ga_param%ga_iatfix))then
   ABI_FREE(ga_param%ga_iatfix)
 endif
 call ga_nullify(ga_param)

end subroutine ga_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_ga/ga_nullify
!! NAME
!!  ga_nullify
!!
!! FUNCTION
!!  Nullify the content of a datastructure of type ga_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ga_param=datastructure of type ga_type.
!!            several parameters for Genetic algorithm search.
!!
!! PARENTS
!!      m_ga
!!
!! CHILDREN
!!
!! SOURCE

subroutine ga_nullify(ga_param)

!Arguments ------------------------------------
!scalars
 type(ga_type),intent(inout) :: ga_param

!************************************************************************

 ga_param%ga_start  = -1
 ga_param%ga_n_rules   = -1
! ga_param%ga_algor     = -1
 ga_param%ga_opt_percent  = zero
 ga_param%ga_rules(:) = -1

end subroutine ga_nullify
!!***

END MODULE m_ga
!!***

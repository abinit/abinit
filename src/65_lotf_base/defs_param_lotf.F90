!!****m* ABINIT/defs_param_lotf
!! NAME
!! defs_param_lotf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module defs_param_lotf
 use defs_basis
 implicit none

 private

 type lotf_control_type
   integer :: natom !--number of atoms


   !--Control variables
   integer :: n0
   integer :: classic, version, me, nproc

   integer :: nitex  !--number of LOTF iterations

   !--Atomflags variables
   integer :: nneigx !--max number of neighbours

 end type lotf_control_type

 !--Contains all control variables
 type(lotf_control_type) :: lotfvar


 public ::             &
   lotfvar,            &
   lotfparam_init

CONTAINS !===========================================================
!!***


!!****f* defs_param_lotf/lotfparam_init
!! NAME
!! lotfparam_init
!!
!! FUNCTION
!!  set some internal variable of lotf
!! INPUTS
!!  natom=number of atoms
!!  version=set type of MD algo
!!  nstart=initial step
!!  nitex=number of LOTF steps
!!  nneigx=roughly decide the number of neighbours
!!  upd=....
!!  me,nproc =disabled parallel LOTF
!!  classic=stick with the adaptable Glue model (rough version)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

 subroutine lotfparam_init(natom,version,nstart,nitex,neeigx,&
   &                    classic,me,nproc)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: natom,version,nstart,neeigx
  integer,intent(in) :: classic,me,nproc
  integer,intent(in) :: nitex

  !--switch number of atom, LOTF notation
  lotfvar%natom = natom

  !--set type of MD algo
  lotfvar%version = version

  !--probably useful in upd_lis0
  lotfvar%n0 = nstart

  !--number of LOTF steps :
  lotfvar%nitex = nitex

  !--roughly decide the number of neighbours :
  lotfvar%nneigx = neeigx

  !--disable LOTF parallel version :
  lotfvar%me  = me
  lotfvar%nproc = nproc

  !--stick with the adaptable Glue model (rough version):
  lotfvar%classic = classic

 end subroutine lotfparam_init

end module defs_param_lotf
!!***

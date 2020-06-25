!!****m* ABINIT/tools_lotf
!! NAME
!! tools_lotf
!!
!! FUNCTION
!!  Contains simple functions for LOTF
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module tools_lotf

 use defs_basis
 use m_errors

 implicit none

 private

 public ::             &
   pinterp,            &
   pinterp_nolinear,   &
   dlvsum,             &
   icf


contains
!!***

 !!****f* tools_lotf/pinterp
 !! NAME
 !! pinterp
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!
 !! SOURCE
 subroutine  pinterp(a0,a1,ainterpoled,ndim,nitex,n)
  ! interpolator : a0 for n=0
  !                a1 for n=nx
  !                a1 for nx = 0
  implicit none

  !Arguments ------------------------
  integer,intent(in) :: ndim,nitex,n 
  real(dp),intent(in) :: a0(ndim),a1(ndim)
  real(dp),intent(out) :: ainterpoled(ndim)
  !Local --------------------------- 
  character(len=500) :: message

! *************************************************************************

  !--Control
  if(nitex < 0) then
    write(message,'(a,i4,a)') '  LOTF: pinterp nitex =',nitex,'smaller than 0'
    MSG_ERROR(message)

  elseif(nitex >= 1) then
    ainterpoled = a0 + (n/real(nitex,dp))*(a1-a0) 
  elseif(nitex == 0) then
    ainterpoled = a1 
  endif
 end subroutine pinterp
 !!***


 !!****f* tools_lotf/pinterp_nolinear
 !! NAME
 !! pinterp_nolinear
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!
 !! CHILDREN
!!
 !! SOURCE
 subroutine  pinterp_nolinear(a0,a1,ainterpoled,ndim,nitex,n)
  ! interpolator : a0 for n=0
  !                a1 for n=nx
  !                a1 for nx = 0
  implicit none

  !Arguments ------------------------
  integer,intent(in) :: ndim,nitex,n 
  real(dp),intent(in) :: a0(ndim),a1(ndim)
  real(dp),intent(out) :: ainterpoled(ndim)
  !Local --------------------------- 
  real(dp) :: lambda

! *************************************************************************

  !--Control
  lambda = n/real(nitex,dp)
  ainterpoled = a0 - six*(a1-a0)*(lambda**3/three-lambda**2/two) 
 end subroutine pinterp_nolinear
 !!***


 !!****f* tools_lotf/dlvsum
 !! NAME
 !! dlvsum
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!
 !! SOURCE
 subroutine dlvsum(local_pe,npes,sum_local_pe,ndimve)
  !------------------------------------------------------------
  ! (ADV) GLOBAL VECTOR SUM, REAL 28/5/94
  implicit none

  !Arguments ------------------------
  integer :: local_pe,npes,ndimve
  real(dp) :: sum_local_pe(ndimve)

! *************************************************************************

  ! do nothing here !!! 
  ! better remove dlvsum  everywhere in this version.
  return
 end subroutine dlvsum
 !!***


 !!****f* tools_lotf/icf
 !! NAME
 !! icf
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
 !!
 !! CHILDREN
 !!
 !! SOURCE
 FUNCTION icf(ix,iy,iz,ic1,ic2,ic3) 

  implicit none
  !Arguments ------------------------
  integer :: ix,iy,iz,ic1,ic2,ic3,icf  
  icf = 1 &
    + mod( ix - 1 + ic1, ic1 ) &
    + mod( iy - 1 + ic2, ic2 ) * ic1 &
    + mod( iz - 1 + ic3, ic3 ) * ic1 * ic2
 end FUNCTION icf

end module tools_lotf
!!***

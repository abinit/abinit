!!****m* ABINIT/pbc_lotf
!! NAME
!! pbc_lotf
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

module pbc_lotf

 use defs_basis

 implicit none
 private

 real(dp),private :: aa(3,3),bb(3,3)
 real(dp),public :: rd(3)  ! these are pbc output
 real(dp),public  :: r2  ! these are pbc output

 public ::           &
   dist_pbc,         &
   pbc_init,         &
   pbc_aa_contract,  &
   pbc_bb_contract,  &
   pbc_bb_proj

 private ::          &
   vecp


 interface dist_pbc
   module procedure dist_pbc_ext
   module procedure dist_pbc_int
 end interface dist_pbc

contains
!!***

!!****f* pbc_lotf/vecp
!! NAME
!! vecp
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_pbc_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine vecp(a,b,c)

  real(dp),intent(in) :: a(3),b(3)
  real(dp),intent(out) :: c(3)

! *************************************************************************

  c(1) = a(2) * b(3) - b(2) * a(3)
  c(2) = a(3) * b(1) - b(3) * a(1)
  c(3) = a(1) * b(2) - b(1) * a(2)
 end subroutine vecp
 !!***

!!****f* pbc_lotf/pbc_init
!! NAME
!! pbc_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine pbc_init(rprimd)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: rprimd(3,3)
  !Local ---------------------------
  real(dp) :: avol

  aa(:,:) =  rprimd(:,:)
  call vecp(aa(1,2),aa(1,3),bb(1,1))     !--b^c
  call vecp(aa(1,3),aa(1,1),bb(1,2))     !--c^a
  call vecp(aa(1,1),aa(1,2),bb(1,3))     !--a^b
  avol = dot_product(aa(:,3),bb(:,3))    !--c.(a^b)

  bb   = (one/avol)*bb
 end subroutine pbc_init
 !!***

!!****f* pbc_lotf/pbc_aa_contract
!! NAME
!! pbc_aa_contract
!!
!! FUNCTION
!!  Compute the contraction of aa
!!
!! INPUTS
!!
!! CHILDREN
!!
!! SOURCE

 function pbc_aa_contract()

  implicit none

  !Arguments ------------------------
  real(dp) :: pbc_aa_contract(3)
  pbc_aa_contract = sqrt(sum(aa(:,:)**two,dim=1))
  return
 end function pbc_aa_contract
 !!***

!!****f* pbc_lotf/pbc_bb_contract
!! NAME
!! pbc_bb_contract
!!
!! FUNCTION
!!  Compute the contraction of bb
!! INPUTS
!!
!! CHILDREN
!!
!! SOURCE

 function pbc_bb_contract()

  implicit none

  !Arguments ------------------------
  real(dp) :: pbc_bb_contract(3)

! *************************************************************************

  pbc_bb_contract = sqrt(sum(bb(:,:)**two,dim=1))
  return
 end function pbc_bb_contract
 !!***


!!****f* pbc_lotf/pbc_bb_proj
!! NAME
!! pbc_bb_proj
!!
!! FUNCTION
!!  Compute the application of a vector on bb
!! INPUTS
!!  vi(3)=real vector
!! CHILDREN
!!
!! SOURCE

 function pbc_bb_proj(vi)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: vi(3)
  real(dp) :: pbc_bb_proj(3)

! *************************************************************************

  pbc_bb_proj = matmul(vi,bb)
  return
 end function pbc_bb_proj
 !!***

!!****f* pbc/dist_pbc_ext
!! NAME
!! dist_pbc_ext
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine dist_pbc_ext(RI,RJ,r2,RD)
  ! ONLY aa AND bb MATRICES ARE USED IN THIS VERSION
  implicit none

  !Arguments ------------------------
  real(dp),intent(out) :: r2
  real(dp),intent(in),dimension(3) :: RI, RJ
  real(dp),intent(out),dimension(3) :: RD
  !Local ---------------------------
  integer :: ii
  real(dp),dimension(3) :: rad, radi

! *************************************************************************

  !--at These are cartesian:
  rad = RI - RJ

  !--at These are monoclinic:
  radi(:) = matmul(rad,bb)

  do ii=1,3
    if(radi(ii) < -half) then
      rad = rad + aa(:,ii)
    elseif(radi(ii) > half) then
      rad = rad - aa(:,ii)
    endif
  enddo

  r2 = dot_product(rad,rad)
  RD = rad
 end subroutine dist_pbc_ext
 !!***


!!****f* pbc/dist_pbc_int
!! NAME
!! dist_pbc_int
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine dist_pbc_int(RI,RJ)
  ! ONLY aa AND bb MATRICES ARE USED IN THIS VERSION
  implicit none

  !Arguments ------------------------
  real(dp),intent(in),dimension(3) :: RI, RJ
  !Local ---------------------------
  integer :: ii
  real(dp),dimension(3) :: rad, radi

! *************************************************************************

  !--at These are cartesian:
  rad = RI - RJ

  !--at These are monoclinic:
  radi(:) = matmul(rad,bb)

  do ii=1,3
    if(radi(ii) < -half) then
      rad = rad + aa(:,ii)
    elseif(radi(ii) > half) then
      rad = rad - aa(:,ii)
    endif
  enddo

  r2 = dot_product(rad,rad)
  RD = rad
 end subroutine dist_pbc_int

end module pbc_lotf
!!***

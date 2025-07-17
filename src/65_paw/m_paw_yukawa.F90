!!****m* m_paw_yukawa/m_paw_yukawa
!! NAME
!!  m_paw_yukawa
!!
!! FUNCTION
!!  This module contains several routines related to the Yukawa parametrization
!!  of Coulomb interactions in the PAW approach.
!!
!! COPYRIGHT
!! Copyright (C) 2025-2025 ABINIT group
!! These routines are inspired by K. Haule routines in embedded DMFT.
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_yukawa

 use defs_basis
 use m_abicore
 use m_errors
 use m_pawrad, only : pawrad_type

 implicit none

 private

 public :: compute_slater
 public :: get_lambda

CONTAINS  !========================================================================================
!!***

!!****f* m_paw_yukawa/compute_slater
!! NAME
!! compute_slater
!!
!! FUNCTION
!!
!! Compute Slater integrals for a screened Yukawa potential v(r,r') = exp(-lambda*(r-r'))/(epsilon*(r-r'))
!! This is eq. 36 in the supplementary of Physical review letters, Haule, K. (2015), 115(19), 196403
!!
!! INPUTS
!!  lpawu = angular momentum
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  proj2 = u(r)**2 where u(r) is the atomic orbital multiplied by r
!!  meshsz = size of the radial mesh
!!  lambda, eps = parameters of the Yukawa potential
!!
!! OUTPUT
!!  fk(lpawu+1)= Slater integrals
!!
!! SOURCE

 subroutine compute_slater(lpawu,pawrad,proj2,meshsz,lambda,eps,fk)

 use m_pawrad, only : pawrad_type,simp_gen
 use m_bessel2, only : bessel_iv,bessel_kv

!Arguments ------------------------------------
 integer, intent(in) :: lpawu,meshsz
 real(dp), intent(in) :: lambda,eps
 real(dp), intent(in) :: proj2(meshsz)
 real(dp), intent(inout) :: fk(lpawu+1)
 type(pawrad_type), intent(in) :: pawrad
!Local variables ------------------------------
 integer :: ir,k,mesh_type
 real(dp) :: dum,r_for_intg,y0
 real(dp), allocatable :: r_k(:),u_inside(:),u_outside(:),y1(:),y2(:)
 !************************************************************************

 mesh_type  = pawrad%mesh_type
 r_for_intg = pawrad%rad(meshsz)
 ABI_MALLOC(r_k,(meshsz))
 ABI_MALLOC(u_inside,(meshsz))
 ABI_MALLOC(u_outside,(meshsz))

 if (lambda == zero) then

   do k=0,2*lpawu+1,2
     u_inside(1) = zero
     r_k(:) = pawrad%rad(1:meshsz)**k
     do ir=2,meshsz
       if (ir == 2) then
         ! Use a trapezoidal rule
         u_inside(ir) = half * (proj2(1)*r_k(1)+proj2(2)*r_k(2)) * &
                 & (pawrad%rad(2)-pawrad%rad(1))
       else if (ir == 3 .and. mesh_type == 3) then
         ! simp_gen doesn't handle this case, so we use a trapezoidal rule instead
         u_inside(ir) = u_inside(2) + half*(proj2(2)*r_k(2)+proj2(3)*r_k(3))* &
                 & (pawrad%rad(3)-pawrad%rad(2))
       else
         ! Use Simpson rule when enough points are available
         call simp_gen(u_inside(ir),proj2(1:ir)*r_k(1:ir),pawrad,r_for_intg=pawrad%rad(ir))
       end if ! ir=2
     end do ! ir
     u_outside(1) = zero
     u_outside(2:meshsz) = two * u_inside(2:meshsz) * proj2(2:meshsz) / (pawrad%rad(2:meshsz)*r_k(2:meshsz))
     call simp_gen(fk(k/2+1),u_outside(:),pawrad,r_for_intg=r_for_intg)
   end do ! k

 else

   ABI_MALLOC(y1,(meshsz))
   ABI_MALLOC(y2,(meshsz))

   r_k(:) = sqrt(pawrad%rad(1:meshsz))

   do k=0,2*lpawu+1,2

     y0 = zero
     if (k == 0) y0 = lambda * sqrt(two/pi)
     y1(1) = y0
     u_inside(1) = zero

     do ir=2,meshsz

       call bessel_iv(half+dble(k),lambda*pawrad%rad(ir),zero,y1(ir),dum)
       y1(ir) = y1(ir) / r_k(ir)

       if (ir == 2) then
         ! Use a trapezoidal rule
         u_inside(ir) = half * (proj2(1)*y1(1)+proj2(2)*y1(2)) * &
                 & (pawrad%rad(2)-pawrad%rad(1))
       else if (ir == 3 .and. mesh_type == 3) then
         ! simp_gen doesn't handle this case, so we use a trapezoidal rule instead
         u_inside(ir) = u_inside(2) + half*(proj2(2)*y1(2)+proj2(3)*y1(3)) * &
                 & (pawrad%rad(3)-pawrad%rad(2))
       else
         ! Use Simpson rule when enough points are available
         call simp_gen(u_inside(ir),proj2(1:ir)*y1(1:ir),pawrad,r_for_intg=pawrad%rad(ir))
       end if ! ir

       call bessel_kv(half+dble(k),lambda*pawrad%rad(ir),zero,y2(ir),dum)
       y2(ir) = y2(ir) / r_k(ir)
     end do ! ir

     u_outside(1) = zero
     u_outside(2:meshsz) = two * (two*dble(k)+one)*u_inside(2:meshsz)*proj2(2:meshsz)*y2(2:meshsz)
     call simp_gen(fk(k/2+1),u_outside(:),pawrad,r_for_intg=r_for_intg)

   end do ! k

   ABI_FREE(y1)
   ABI_FREE(y2)

 end if ! lambda

 fk(1:lpawu+1) = fk(1:lpawu+1) / eps

 ABI_FREE(r_k)
 ABI_FREE(u_inside)
 ABI_FREE(u_outside)

 end subroutine compute_slater
!!***

!----------------------------------------------------------------------

!!****f* m_paw_yukawa/get_lambda
!! NAME
!! get_lambda
!!
!! FUNCTION
!!
!! Conversion from U,J,f4/f2,f6/f2 parametrization of Slater integrals
!! to lambda,epsilon parametrization, with lambda and epsilon the parameters
!! of the screened Yukawa potential v(r,r') = exp(-lambda*(r-r'))/(epsilon*(r-r')).
!!
!! CAREFUL: this routine does not handle custom f4/f2 and f6/f2 values, and set them
!!          to their default values f4/f2=0.625 for l=2 and f4/f2=0.6681, f6/f2=0.4943 for l=3.
!!
!! CAREFUL: For l>=2 we lose information since we have more input parameters than output
!!          parameters. In this case, a compromise has to be made, and you will no longer
!!          have the exact same Slater integrals as before.
!!
!! INPUTS
!!  lpawu = angular momentum
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  proj2 = u(r)**2 where u(r) is the atomic orbital multiplied by r
!!  meshsz = size of the radial mesh
!!  upawu,jpawu = parameters for Slater integrals
!!  yukawa_param = if set to 1, search for lambda and epsilon yielding the values closest to u and j
!!                 if set to 2, search for lambda yielding u, and set epsilon to 1
!!
!! OUTPUT
!!  lambda,epsilon = parameters of the corresponding Yukawa potential
!!
!! SOURCE

 subroutine get_lambda(lpawu,pawrad,proj2,meshsz,upawu,jpawu,lambda,eps,yukawa_param)

 use m_brentq, only : brentq
 use m_hybrd, only : hybrd

!Arguments ------------------------------------
 integer, intent(in) :: lpawu,meshsz,yukawa_param
 real(dp), intent(in) :: upawu,jpawu
 real(dp), intent(out) :: lambda,eps
 real(dp), intent(in) :: proj2(meshsz)
 type(pawrad_type), intent(in) :: pawrad
!Local variables ------------------------------
 integer  :: i,ierr,info,ldfjac,lr,maxfev,ml,mode,mu,n,nfev,nprint
 real(dp) :: epsfcn,fac,lmb_temp,upbound,xtol
 real(dp) :: diag(2),fjac(2,2),fkk(lpawu+1),fvec(2),lmb_eps(2)
 real(dp) :: r(3),qtf(2),wa1(2),wa2(2),wa3(2),wa4(2)
 character(len=500) :: message
 !************************************************************************

 ! Find suitable upper bound for brentq routine
 upbound = five
 do i=1,10

   call compute_slater(lpawu,pawrad,proj2(:),meshsz,upbound,one,fkk(:))
   if (fkk(1) < upawu) exit
   upbound = two * upbound

 end do ! i

 write(message,'(4a)') "An error occurred when trying to find a suitable lambda ", &
                     & "for your input values of upawu.", ch10, &
                     & "Either try a different value of upawu or use dmft_lambda_yukawa."

 if (fkk(1) > upawu) ABI_ERROR(message)

 ! First, set epsilon to 1, and find lambda which yields the correct F0=upawu, to have a good starting point
 call brentq(get_coulomb_u,zero,upbound,two*tol12,four*epsilon(one),100,lmb_temp,ierr)

 if (ierr == 0) ABI_ERROR(message)

 ! Initial values for lambda and epsilon
 lmb_eps(1) = lmb_temp
 lmb_eps(2) = one

 lambda = lmb_temp
 eps    = one

 if (yukawa_param == 1) return

 ! This part finds epsilon (not used anymore)
 if (lpawu > 0) then

   ! Default values from scipy
   epsfcn = epsilon(one) ; fac = dble(100.) ; n = 2
   ldfjac = n ; lr = n * (n+1) / 2
   maxfev = 200 * (n+1) ; ml = n - 1 ; mode = 1
   mu = n - 1 ; nprint = 0 ; xtol = dble(1.49012e-8)

   ! Now find lambda and epsilon
   call hybrd(get_coulomb_uj,2,lmb_eps(:),fvec(:),xtol,maxfev,ml,mu,epsfcn,diag(:),mode, &
            & fac,nprint,info,nfev,fjac(:,:),ldfjac,r(:),lr,qtf(:),wa1(:),wa2(:),wa3(:),wa4(:))

   if (info /= 1) ABI_ERROR(message)

 end if ! lpawu > 0

 lambda = lmb_eps(1)
 eps    = lmb_eps(2)

 contains

 subroutine get_coulomb_u(lmb,uu)

!Arguments ------------------------------------
 real(dp), intent(in) :: lmb
 real(dp), intent(out) :: uu
!Local variables ------------------------------
 real(dp) :: fk(lpawu+1)
!************************************************************************

 call compute_slater(lpawu,pawrad,proj2(:),meshsz,lmb,one,fk(:))
 uu = fk(1) - upawu

 end subroutine get_coulomb_u

 subroutine get_coulomb_uj(n,lmb_eps,uj,iflag)

!Arguments ------------------------------------
 integer, intent(in) :: iflag,n
 real(dp), intent(in) :: lmb_eps(n)
 real(dp), intent(inout) :: uj(n)
 !Local variables ------------------------------
 real(dp) :: eps,f4of2,f6of2,factor,j2,j4,j6,jh,lmb
 real(dp) :: fk(lpawu+1)
 character(len=500) :: message
!************************************************************************

 ABI_UNUSED(iflag)

 lmb = lmb_eps(1)
 eps = lmb_eps(2)
 call compute_slater(lpawu,pawrad,proj2(:),meshsz,lmb,eps,fk(:))
 uj(1) = fk(1) - upawu

 if (lpawu == 1) then
   j2 = fk(2) * fifth
   jh = j2
 else if (lpawu == 2) then
   f4of2  = dble(0.625)
   factor = (one+f4of2) / dble(14)
   j2 = fk(2)
   j4 = fk(3) / f4of2
   jh = (j2+j4) * factor * half
 else if (lpawu == 3) then
   f4of2  = dble(0.6681)
   f6of2  = dble(0.4943)
   factor = (dble(286.)+dble(195.)*f4of2+dble(250.)*f6of2) / dble(6435.)
   j2 = fk(2)
   j4 = fk(3) / f4of2
   j6 = fk(4) / f6of2
   jh = (j2+j4+j6) * factor * third
 else
   write(message,'(a,i0,2a)') ' lpawu=',lpawu,ch10,' lpawu not equal to 0, 1, 2 or 3 is not allowed'
   ABI_ERROR(message)
 end if ! lpawu

 uj(2) = jh - jpawu

 end subroutine get_coulomb_uj

 end subroutine get_lambda
!!***

!----------------------------------------------------------------------

END MODULE m_paw_yukawa
!!***

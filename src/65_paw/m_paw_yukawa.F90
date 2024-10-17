!!****m* m_paw_yukawa/m_paw_yukawa
!! NAME
!!  m_paw_yukawa
!!
!! FUNCTION
!!  This module contains several routines related to the Yukawa parametrization
!!  of Coulomb interactions in the PAW approach.
!!
!! COPYRIGHT
!! Copyright (C) 2024 K. Haule 
!! These routines were translated from embedded DMFT.
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
 use m_errors
 use m_pawrad, only : pawrad_type

 implicit none

 private

!public procedures.
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
!!  Fk(lpawu+1)= Slater integrals
!!
!! SOURCE

 subroutine compute_slater(lpawu,pawrad,proj2,meshsz,lambda,eps,Fk)
 
 use m_pawrad, only : pawrad_type,simp_gen
 use m_bessel2, only : bessel_iv,bessel_kv

!Arguments ------------------------------------
 integer, intent(in) :: lpawu,meshsz
 type(pawrad_type), intent(in) :: pawrad
 real(dp), intent(in) :: proj2(:)
 real(dp), intent(in) :: lambda,eps
 real(dp), intent(out) :: Fk(:)
!Local variables ------------------------------
 integer :: ir,k,mesh_type
 real(dp) :: dum1,dum2,r_for_intg,y0
 real(dp), allocatable :: U_inside(:),U_outside(:),y1(:),y2(:)
 !************************************************************************

 mesh_type  = pawrad%mesh_type 
 r_for_intg = pawrad%rad(meshsz)
 ABI_MALLOC(U_inside,(meshsz)) 
 ABI_MALLOC(U_outside,(meshsz)) 

 if (lambda == zero) then
   do k=0,2*lpawu+1,2
     U_inside(1) = zero
     do ir=2,meshsz
       if (ir == 2) then
         ! Use a trapezoidal rule 
         U_inside(ir) = half * (proj2(1)*(pawrad%rad(1)**k)+proj2(2)*(pawrad%rad(2)**k)) * &
                 & (pawrad%rad(2)-pawrad%rad(1))
       else if (ir == 3 .and. mesh_type == 3) then 
         ! simp_gen doesn't handle this case, so we use a trapezoidal rule instead
         U_inside(ir) = U_inside(2) + half*(proj2(2)*(pawrad%rad(2)**k)+proj2(3)*(pawrad%rad(3)**k))* &
                 & (pawrad%rad(3)-pawrad%rad(2))
       else 
         ! Use Simpson rule when enough points are available 
         call simp_gen(U_inside(ir),proj2(1:ir)*(pawrad%rad(1:ir)**k),pawrad,r_for_intg=pawrad%rad(ir))
       end if ! ir=2 
     end do ! ir
     U_outside(1) = zero
     U_outside(2:meshsz) = two * U_inside(2:meshsz) * proj2(2:meshsz) / (pawrad%rad(2:meshsz)**(k+1))
     call simp_gen(Fk(k/2+1),U_outside(:),pawrad,r_for_intg=r_for_intg) 
   end do ! k

 else
 
   ABI_MALLOC(y1,(meshsz)) 
   ABI_MALLOC(y2,(meshsz))
 
   do k=0,2*lpawu+1,2
    
     y0 = zero
     if (k == 0) y0 = lambda * sqrt(two/pi)
     y1(1) = y0 
     U_inside(1) = zero
     
     do ir=2,meshsz
     
       dum1 = zero ! very important
       
       call bessel_iv(half+dble(k),lambda*pawrad%rad(ir),dum1,y1(ir),dum2)
       y1(ir) = y1(ir) / sqrt(pawrad%rad(ir))
       
       if (ir == 2) then
         ! Use a trapezoidal rule
         U_inside(ir) = half * (proj2(1)*y1(1)+proj2(2)*y1(2)) * &
                 & (pawrad%rad(2)-pawrad%rad(1))
       else if (ir == 3 .and. mesh_type == 3) then
         ! simp_gen doesn't handle this case, so we use a trapezoidal rule instead
         U_inside(ir) = U_inside(2) + half*(proj2(2)*y1(2)+proj2(3)*y1(3)) * &
                 & (pawrad%rad(3)-pawrad%rad(2))
       else
         ! Use Simpson rule when enough points are available
         call simp_gen(U_inside(ir),proj2(1:ir)*y1(1:ir),pawrad,r_for_intg=pawrad%rad(ir))
       end if ! ir
       
       dum1 = zero ! very important
       call bessel_kv(half+dble(k),lambda*pawrad%rad(ir),dum1,y2(ir),dum2)
       y2(ir) = y2(ir) / sqrt(pawrad%rad(ir))
     end do ! ir
     
     U_outside(1) = zero
     U_outside(2:meshsz) = two * (two*dble(k)+one)*U_inside(2:meshsz)*proj2(2:meshsz)*y2(2:meshsz)
     call simp_gen(Fk(k/2+1),U_outside(:),pawrad,r_for_intg=r_for_intg)
     
   end do ! k
   
   ABI_FREE(y1)
   ABI_FREE(y2)
   
 end if ! lambda

 Fk(1:lpawu+1) = Fk(1:lpawu+1) / eps

 ABI_FREE(U_inside) 
 ABI_FREE(U_outside) 

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
!!
!! OUTPUT
!!  lambda,epsilon = parameters of the corresponding Yukawa potential
!!
!! SOURCE

 subroutine get_lambda(lpawu,pawrad,proj2,meshsz,upawu,jpawu,lambda,eps)
 
 use m_brentq, only : brentq
 use m_hybrd, only : hybrd

!Arguments ------------------------------------
 integer, intent(in) :: lpawu,meshsz
 type(pawrad_type), intent(in) :: pawrad
 real(dp), intent(in) :: proj2(:)
 real(dp), intent(in) :: upawu,jpawu
 real(dp), intent(out) :: lambda,eps
!Local variables ------------------------------
 integer  :: i,ierr,info,ldfjac,lr,maxfev,ml,mode,mu,n,nfev,nprint
 real(dp) :: epsfcn,fac,lmb_temp,upbound,xtol
 real(dp) :: diag(2),fjac(2,2),Fkk(lpawu+1),fvec(2),lmb_eps(2)
 real(dp) :: r(3),qtf(2),wa1(2),wa2(2),wa3(2),wa4(2)
 character(len=500) :: message
 !************************************************************************
 
 ! Find suitable upper bound for brentq routine
 upbound = five
 do i=1,10

   call compute_slater(lpawu,pawrad,proj2(:),meshsz,upbound,one,Fkk(:))
   if (Fkk(1) < upawu) exit
   upbound = two * upbound
 
 end do ! i

 if (Fkk(1) > upawu) ABI_ERROR("Could not find a suitable lambda in get_lambda subroutine")

 ! First, set epsilon to 1, and find lambda which yields the correct F0=upawu, to have a good starting point
 call brentq(get_coulomb_u,zero,upbound,two*tol12,four*epsilon(one),100,lmb_temp,ierr)

 if (ierr == 0) ABI_ERROR("The brentq method did not converge in the get_lambda subroutine")

 ! Initial values for lambda and epsilon
 lmb_eps(1) = lmb_temp 
 lmb_eps(2) = one

 if (lpawu > 0) then
 
   ! Default values from scipy
   n = 2 
   xtol = dble(1.49012e-8) 
   maxfev = 200 * (n+1) 
   ml = n - 1 
   mu = n - 1
   epsfcn = epsilon(one) 
   mode = 1 
   fac = dble(100.) 
   nprint = 0 
   ldfjac = n 
   lr = n * (n+1) / 2
 
   ! Now find lambda and epsilon
   call hybrd(get_coulomb_uj,2,lmb_eps(:),fvec(:),xtol,maxfev,ml,mu,epsfcn,diag(:),mode, &
            & fac,nprint,info,nfev,fjac(:,:),ldfjac,r(:),lr,qtf(:),wa1(:),wa2(:),wa3(:),wa4(:))
  
   if (info /= 1) then
     write(message,*) "Error in hybrd, info is equal to",info
     ABI_ERROR(message)
   end if
 
 end if ! lpawu > 0

 lambda = lmb_eps(1) 
 eps    = lmb_eps(2)

 contains

 subroutine get_coulomb_u(lmb,uu)
 
!Arguments ------------------------------------
 real(dp), intent(in) :: lmb
 real(dp), intent(out) :: uu
!Local variables ------------------------------
 real(dp) :: Fk(lpawu+1)
!************************************************************************

 call compute_slater(lpawu,pawrad,proj2(:),meshsz,lmb,one,Fk(:))
 uu = Fk(1) - upawu
 
 end subroutine get_coulomb_u

 subroutine get_coulomb_uj(n,lmb_eps,uj,iflag)

!Arguments ------------------------------------
 integer, intent(in) :: iflag,n
 real(dp), intent(in) :: lmb_eps(n)
 real(dp), intent(inout) :: uj(n)
 !Local variables ------------------------------
 real(dp) :: eps,f4of2,f6of2,factor,J2,J4,J6,Jh,lmb
 real(dp) :: Fk(lpawu+1)
 character(len=500) :: message
!************************************************************************

 ABI_UNUSED(iflag)

 lmb = lmb_eps(1) 
 eps = lmb_eps(2)
 call compute_slater(lpawu,pawrad,proj2(:),meshsz,lmb,eps,Fk(:))
 uj(1) = Fk(1) - upawu

 if (lpawu == 1) then
   J2 = Fk(2) * fifth 
   Jh = J2
 else if (lpawu == 2) then
   f4of2  = dble(0.625) 
   factor = (one+f4of2) / dble(14)
   J2 = Fk(2) 
   J4 = Fk(3) / f4of2
   Jh = (J2+J4) * factor * half
 else if (lpawu == 3) then
   f4of2  = dble(0.6681) 
   f6of2  = dble(0.4943) 
   factor = (dble(286.)+dble(195.)*f4of2+dble(250.)*f6of2) / dble(6435.)
   J2 = Fk(2)  
   J4 = Fk(3) / f4of2 
   J6 = Fk(4) / f6of2 
   Jh = (J2+J4+J6) * factor * third
 else
   write(message,'(a,i0,2a)') ' lpawu=',lpawu,ch10,' lpawu not equal to 0, 1, 2 or 3 is not allowed'
   ABI_ERROR(message)
 end if ! lpawu
         
 uj(2) = Jh - jpawu

 end subroutine get_coulomb_uj

 end subroutine get_lambda
!!***

!----------------------------------------------------------------------

END MODULE m_paw_yukawa
!!***

!!****m* m_paw_exactDC/m_paw_exactDC
!! NAME
!!  m_paw_exactDC
!!
!! FUNCTION
!!  This module contains several routines related to the exact formula for the double counting.
!!
!! COPYRIGHT
!! Copyright (C) 2025-2025 ABINIT group
!! These routines were inspired by K. Haule routines in embedded DMFT.
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_exactDC

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private

 public :: compute_exactDC

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_exactDC/compute_exactDC
!! NAME
!! compute_exactDC
!!
!! FUNCTION
!!
!! Compute the exact formula for the double counting.
!! See Physical review letters, Haule, K. (2015), 115(19), 196403 for formula.
!!
!! INPUTS
!!  lpawu = angular momentum for correlated species
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  occ(2*lpawu+1,2*lpawu+1) = occupation matrix (summed over spins) in
!!    the real spherical harmonics basis, with the convention
!!    occ(i,j) = <c_j^dagger c_i>
!!  ixc = index of the XC functional
!!
!! OUTPUT
!!  vdc(2*lpawu+1,2*lpawu+1) = double counting potential
!!  edc = double counting energy
!!  edcdc = integral of Vdc(r)*rho_loc(r)
!!
!! SOURCE

 subroutine compute_exactDC(lpawu,pawtab,pawrad,occ,vdc,edc,edcdc,ixc)

 use m_pawtab, only : pawtab_type
 use m_pawrad, only : pawrad_type,simp_gen
 use m_paw_sphharm, only : ylmc,ylmcd
 use m_splines, only : spline2

!Arguments ------------------------------------
 integer, intent(in) :: ixc,lpawu
 type(pawtab_type), intent(in) :: pawtab
 type(pawrad_type), intent(in) :: pawrad
 complex(dpc), intent(in) :: occ(2*lpawu+1,2*lpawu+1)
 complex(dpc), intent(inout) :: vdc(2*lpawu+1,2*lpawu+1)
 real(dp), intent(out) :: edc,edcdc
!Local variables ------------------------------
 integer :: i,ir,j,k,l,m,m1,mm,minusm,meshsz,ndim
 logical :: need_gradient
 real(dp) :: cphi,ctheta,ecc,excint,exx,grad,gradphi,gradr,gradth
 real(dp) :: onemsqrt2,phi,rad,rho,rhor,sexc,sphi,stheta,svxc
 real(dp) :: theta,vcc,vxcf,vxx
 real(dp) :: kcart(3)
 complex(dpc) :: dphi,dth,ylm_c
 integer, parameter :: ln = 13
 real(dp), allocatable :: exc(:,:,:),exci(:),phis(:),proj_dr(:),rho_angle(:,:)
 real(dp), allocatable :: rho_angle_dphi(:,:),rho_angle_dtheta(:,:)
 real(dp), allocatable :: slm_dphi(:,:,:),slm_dtheta(:,:,:),thetas(:),tweights(:)
 real(dp), allocatable :: vxc(:,:,:),vxci(:,:,:),ylm(:,:,:)
 !************************************************************************

 ! VERY IMPORTANT: This routine assumes that we work in the real spherical harmonic basis.
 ! If you want to generalize that to the complex case for some reason, the formulas
 ! below are not valid, so you would need to add some complex conjugates.
 ! Also, we work with the convention occ(i,j) = <c_j^dagger c_i>.

 if (ixc /= 7 .and. ixc /= 11 .and. ixc /= -1012 .and. ixc /= -101130) &
    & ABI_ERROR("Only PW92 and PBE are handled!")

 need_gradient = (ixc == 11 .or. ixc == -101130)

 meshsz = size(pawtab%proj2(:)) ! Size of radial mesh
 ndim = 2*lpawu + 1

 ABI_MALLOC(exc,(ln+1,2*ln+1,meshsz))
 ABI_MALLOC(exci,(meshsz))
 ABI_MALLOC(phis,(2*ln+1))
 ABI_MALLOC(rho_angle,(ln+1,2*ln+1))
 ABI_MALLOC(thetas,(ln+1))
 ABI_MALLOC(tweights,(ln+1))
 ABI_MALLOC(vxc,(ln+1,2*ln+1,meshsz))
 ABI_MALLOC(vxci,(meshsz,ndim,ndim))
 ABI_MALLOC(ylm,(ln+1,2*ln+1,ndim))

 if (need_gradient) then
   ABI_MALLOC(proj_dr,(meshsz))
   ABI_MALLOC(rho_angle_dphi,(ln+1,2*ln+1))
   ABI_MALLOC(rho_angle_dtheta,(ln+1,2*ln+1))
   ABI_MALLOC(slm_dphi,(ln+1,2*ln+1,ndim))
   ABI_MALLOC(slm_dtheta,(ln+1,2*ln+1,ndim))
   call spline2(pawrad%rad(1:meshsz),pawtab%proj2(:)/max(epsilon(one),pawrad%rad(1:meshsz)**2), &
              & meshsz,proj_dr(:),zero,zero,3,3)
   rho_angle_dphi(:,:)   = zero
   rho_angle_dtheta(:,:) = zero
 end if ! gradient

 vdc(:,:) = czero
 edc = zero

 ! Hartree contribution

 do l=1,ndim
   do j=1,ndim
     do k=1,ndim
       do i=1,ndim
         ! In the complex case, you should use conjg(pawtab%vee(i,k,j,l))
         vdc(i,j) = vdc(i,j) + pawtab%vee(i,k,j,l)*occ(k,l)
         edc = edc + pawtab%vee(i,k,j,l)*dble(occ(i,j)*occ(k,l))
       end do ! i
     end do ! k
   end do ! j
 end do ! l

 edc = half * edc

 ! XC contribution

 ! Prepare angular mesh
 call angular_mesh(thetas(:),phis(:),tweights(:),ln)

 ! Compute the real spherical harmonics on the angular mesh
 do m=0,lpawu
   mm = m + lpawu + 1 ; minusm = - m + lpawu + 1
   onemsqrt2 = (-one)**m * sqrt2
   do j=1,2*ln+1
     phi = phis(j) ; cphi = cos(phi) ; sphi = sin(phi)
     do i=1,ln+1
       theta = thetas(i) ; ctheta = cos(theta) ; stheta = sin(theta)
       kcart(1) = cphi * stheta ; kcart(2) = sphi * stheta ; kcart(3) = ctheta
       ylm_c = ylmc(lpawu,m,kcart(:))
       ! Convert complex harmonics to real harmonics
       if (m == 0) then
         ylm(i,j,mm) = dble(ylm_c)
       else
         ylm(i,j,mm)     = onemsqrt2 * dble(ylm_c)
         ylm(i,j,minusm) = onemsqrt2 * aimag(ylm_c)
       end if ! m=0
       if (need_gradient) then
         call ylmcd(lpawu,m,kcart(:),dth,dphi)
         if (m == 0) then
             slm_dphi(i,j,mm)   = dble(dphi)
             slm_dtheta(i,j,mm) = dble(dth)
         else
             slm_dphi(i,j,mm)       = onemsqrt2 * dble(dphi)
             slm_dtheta(i,j,mm)     = onemsqrt2 * dble(dth)
             slm_dphi(i,j,minusm)   = onemsqrt2 * aimag(dphi)
             slm_dtheta(i,j,minusm) = onemsqrt2 * aimag(dth)
         end if ! m=0
       end if ! gradient
     end do ! i
   end do ! j
 end do ! m

 rho_angle(:,:) = zero

 do m1=1,ndim
   do m=1,ndim
     ! In the complex case, you should use ylm(:,:,m)*conjg(ylm(:,:,m1))
     rho_angle(:,:) = rho_angle(:,:) + dble(occ(m,m1))*ylm(:,:,m)*ylm(:,:,m1)
     if (need_gradient) then
       ! In the complex case, there should be a conjg each time there is a m1
       rho_angle_dphi(:,:)   = rho_angle_dphi(:,:) + dble(occ(m,m1))*(slm_dphi(:,:,m)*ylm(:,:,m1)+ylm(:,:,m)*slm_dphi(:,:,m1))
       rho_angle_dtheta(:,:) = rho_angle_dtheta(:,:) + dble(occ(m,m1))*(slm_dtheta(:,:,m)*ylm(:,:,m1)+ylm(:,:,m)*slm_dtheta(:,:,m1))
     end if ! gradient
   end do ! m
 end do ! m1

 if (maxval(rho_angle(:,:)-abs(rho_angle(:,:))) > tol10) ABI_WARNING("WARNING: the density is negative !")
 rho_angle(:,:) = abs(rho_angle(:,:))

 do ir=1,meshsz
   rad  = pawrad%rad(ir)
   rhor = pawtab%proj2(ir) / max(rad**2,epsilon(one))
   do j=1,2*ln+1
     do i=1,ln+1

       rho  = rho_angle(i,j) * rhor
       grad = zero

       if (need_gradient) then
         stheta  = sin(thetas(i))
         gradr   = rho_angle(i,j) * proj_dr(ir)
         gradth  = (rho_angle_dtheta(i,j)*rhor) / max(rad,epsilon(one))
         gradphi = (rho_angle_dphi(i,j)*rhor) / (max(rad,epsilon(one))*stheta)
         grad = gradr*gradr + gradth*gradth + gradphi*gradphi
       end if ! gradient

       call exchange_yukawa(exx,vxx,rho,pawtab%lambda,pawtab%eps,grad,ixc)
       exc(i,j,ir) = exx ; vxc(i,j,ir) = vxx

       call correlation_yukawa(ecc,vcc,rho,pawtab%lambda,pawtab%eps,grad,ixc)
       exc(i,j,ir) = exc(i,j,ir) + ecc ; vxc(i,j,ir) = vxc(i,j,ir) + vcc

     end do ! i
   end do ! j
 end do ! ir

 exci(:) = zero

 ! Integrals over theta and phi

 do m1=1,ndim
   do m=1,ndim
     do ir=1,meshsz
       svxc = zero
       sexc = zero
       do j=1,2*ln+1
         ! In the complex case, you should use ylm(:,j,m)*conjg(ylm(:,j,m1))
         svxc = svxc + sum(ylm(:,j,m)*ylm(:,j,m1)*tweights(:)*vxc(:,j,ir))
         sexc = sexc + sum(ylm(:,j,m)*ylm(:,j,m1)*tweights(:)*exc(:,j,ir))
       end do ! j
       vxci(ir,m,m1) = two_pi * svxc / dble(2*ln+1)
       exci(ir) = exci(ir) + two_pi * sexc * dble(occ(m,m1)) / dble(2*ln+1)
     end do ! ir
   end do ! m
 end do ! m1

 ! Integrals over r

 call simp_gen(excint,exci(:)*pawtab%proj2(:),pawrad,r_for_intg=pawrad%rad(meshsz))
 edc = edc + excint
 edcdc = zero

 ! Careful, vdc(i,j) as defined here is the i,j-th matrix element of the TRANSPOSE of vdc !

 do m1=1,ndim
   do m=1,ndim
     call simp_gen(vxcf,vxci(:,m,m1)*pawtab%proj2(:),pawrad,r_for_intg=pawrad%rad(meshsz))
     vdc(m,m1) = vdc(m,m1) + cmplx(vxcf,zero,kind=dp)
     edcdc = edcdc + dble(vdc(m,m1)*occ(m,m1))
   end do ! m
 end do ! m1

 ABI_FREE(exc)
 ABI_FREE(exci)
 ABI_FREE(phis)
 ABI_FREE(rho_angle)
 ABI_FREE(thetas)
 ABI_FREE(tweights)
 ABI_FREE(vxc)
 ABI_FREE(vxci)
 ABI_FREE(ylm)
 ABI_SFREE(proj_dr)
 ABI_SFREE(rho_angle_dphi)
 ABI_SFREE(rho_angle_dtheta)
 ABI_SFREE(slm_dphi)
 ABI_SFREE(slm_dtheta)

 end subroutine compute_exactDC
!!***

!----------------------------------------------------------------------

!!****f* m_paw_exactDC/angular_mesh
!! NAME
!! angular_mesh
!!
!! FUNCTION
!!
!! Prepare the angular mesh for the integration over thetas and phis.
!!
!! INPUTS
!!  ln = controls the number of integration points
!!
!! OUTPUT
!!  thetas(ln+1) = Gauss-Legendre integration points for theta grid
!!  phis(2*ln+1) = uniform points for phi grid
!!  tweights(ln+1) = Gauss-Legendre weights for theta grid
!!
!! SOURCE

subroutine angular_mesh(thetas,phis,tweights,ln)

use m_numeric_tools, only : coeffs_gausslegint

!Arguments ------------------------------------
 integer, intent(in) :: ln
 real(dp), intent(inout) :: phis(2*ln+1),thetas(ln+1),tweights(ln+1)
!Local variables ------------------------------
 integer :: i,lg
 real(dp) :: dsum,phi
 real(dp), parameter :: fake_shift = exp(-4.0_dp) ! small shift such that we do not start at phi=0
!************************************************************************

 lg = ln + 1
 call coeffs_gausslegint(-1.0_dp,1.0_dp,thetas(:),tweights(:),lg)
 do i=1,lg
   thetas(i) = acos(thetas(i))
 end do

 do i=0,2*lg-2       ! 2*lg-1 points in phi direction
   phi = pi * (two*dble(i)/dble(2*lg-1)+fake_shift)
   phis(i+1) = phi
 end do ! i

 dsum = sum(tweights(:))
 tweights(:) = tweights(:) * 2.0_dp / dsum

end subroutine angular_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_paw_exactDC/exchange_yukawa
!! NAME
!! exchange_yukawa
!!
!! FUNCTION
!!
!! Compute the exchange contribution with a Yukawa potential
!!
!! INPUTS
!!  rho = density at current point
!!  lambda = parameter for Yukawa potential (inverse screening length)
!!  eps = parameter for Yukawa potential (dielectric constant)
!!  grad = square of the gradient of the density
!!  ixc = index of the XC functional
!!
!! OUTPUT
!!  ex = exchange energy per particle
!!  vx = exchange potential
!!
!! SOURCE

subroutine exchange_yukawa(ex,vx,rho,lambda,eps,grad,ixc)

!Arguments ------------------------------------
 integer, intent(in) :: ixc
 real(dp), intent(in) :: grad,lambda,eps,rho
 real(dp), intent(out) :: ex,vx
!Local variables ------------------------------
 logical :: islambda
 real(dp) :: dfdkappa,dfdmu,dfdrho,dfdss,dfx,div,div2,dkappadx,dlogf,dmudx,dssdrho
 real(dp) :: dxdrho,fx,fx_pbe,kappa,mu,rhothird,rhotwothird,rsinv,ss,x,x2
 real(dp), parameter :: c0 = (9.0_dp/(4.0_dp*(pi**2)))**third * three_quarters
 real(dp), parameter :: c1 = 1.804_dp
 real(dp), parameter :: kf_fac = (3.0_dp*(pi**2))**third
 real(dp), parameter :: mu0 = 0.2195149727645171_dp
 real(dp), parameter :: rsinv_fac = (4.0_dp*pi/3.0_dp)**third
 real(dp), parameter :: twotwothird = (2.0_dp)**(2.0_dp*third)
 real(dp), parameter :: x_fac = (9.0_dp*pi/4.0_dp)**third
!************************************************************************

 rhothird = rho**third ; rsinv = rsinv_fac * rhothird
 islambda = (abs(lambda) > tol10)

 ! LDA exchange
 if (islambda) then
   x = x_fac * rsinv / lambda
   call fexchange(x,fx,dfx)
 else
   fx  = 1.0_dp
   dfx = 0.0_dp
 end if

 ex = - c0 * fx * rsinv / eps
 vx = 4.0_dp*ex/3.0_dp - c0*x*dfx*rsinv/(3.0_dp*eps)

 if (ixc == 7 .or. ixc == -1012) return

 ! PBE exchange
 if (abs(rho) < tol30) then
   ex = zero ; vx = zero
   return
 end if

 if (islambda) then
   x2 = x * x
   kappa = x2 / (x2 + twotwothird)
 else
   kappa = 1.0_dp
 end if

 kappa = c1*kappa/fx - 1.0_dp
 rhotwothird = rhothird * rhothird
 ss = grad / (4.0_dp*(rho**2)*rhotwothird*(kf_fac**2))
 mu = mu0 / fx
 div = 1.0_dp+mu*ss/kappa
 div2 = div * div
 fx_pbe = 1.0_dp + kappa*(1.0_dp-1.0_dp/div)
 dssdrho = -8.0_dp * third * ss / rho
 dfdss = mu / div2
 dfdrho = dfdss * dssdrho
 if (islambda) then
   dlogf = dfx / fx
   dxdrho = third * x / rho
   dkappadx = (kappa+1.0_dp) * (2.0_dp*twotwothird/(x*(x2+twotwothird))-dlogf)
   dfdkappa = 1.0_dp - (1.0_dp+2.0_dp*mu*ss/kappa)/div2
   dmudx = -mu * dlogf
   dfdmu = ss / div2
   dfdrho = dfdrho + (dfdkappa*dkappadx+dfdmu*dmudx)*dxdrho
 end if
 vx = vx*fx_pbe + rho*ex*dfdrho
 ex = ex * fx_pbe

end subroutine exchange_yukawa
!!***

!----------------------------------------------------------------------

!!****f* m_paw_exactDC/fexchange
!! NAME
!! fexchange
!!
!! FUNCTION
!!
!! Compute the function in the HEG exchange energy with Yukawa potential
!!
!! INPUTS
!!  x = input of the function (=(9*pi/4)**(1/3) / (lambda*rs)
!!      with lambda the Yukawa parameter and rs the Wigner-Seitz radius)
!!
!! OUTPUT
!!  fx  = value of the function at x
!!  dfx = derivative of the function at x
!!
!! SOURCE

subroutine fexchange(x,fx,dfx)

!Arguments ------------------------------------
 real(dp), intent(in) :: x
 real(dp), intent(out) :: fx,dfx
!Local variables ------------------------------
 real(dp) :: at2,lg2,x2,x3,x4,x5
!************************************************************************

 x2 =  x * x ; x3 = x2 * x ; x4 = x3 * x ; x5 = x4 * x
 if (x < tol2) then ! Taylor expansion
   fx  = 4.0_dp * x2 * (1.0_dp/9.0_dp-2.0_dp*x2/15.0_dp+8.0_dp*x4/35.0_dp)
   dfx = 8.0_dp * x * (1.0_dp/9.0_dp-4.0_dp*x2/15.0_dp+24.0_dp*x4/35.0_dp)
 else
   at2 = atan(2.0_dp*x) ; lg2 = log(1.0_dp+4.0_dp*x2)
   fx  = 1.0_dp - 1.0_dp/(6.0_dp*x2) - 4.0_dp*at2/(3.0_dp*x) + &
      & (1.0_dp+12.0_dp*x2)*lg2/(24.0_dp*x4)
   dfx = (8.0_dp*x3*at2+4.0_dp*x2-(1.0_dp+6.0_dp*x2)*lg2)/(6.0_dp*x5)
 end if

end subroutine fexchange
!!***

!!****f* m_paw_exactDC/correlation_yukawa
!! NAME
!! correlation_yukawa
!!
!! FUNCTION
!!
!! Compute the correlation contribution with a Yukawa potential
!!
!! INPUTS
!!  rho = density at current point
!!  lambda = parameter for screened potential (inverse screening length)
!!  eps = parameter for screened potential (dielectric constant)
!!  grad = square of the gradient of the density
!!  ixc = index of the XC functional
!!
!! OUTPUT
!!  ec = correlation energy per particle
!!  vc = correlation potential
!!
!! SOURCE

subroutine correlation_yukawa(ec,vc,rho,lambda,eps,grad,ixc)

!Arguments ------------------------------------
 integer, intent(in) :: ixc
 real(dp), intent(in) :: eps,grad,lambda,rho
 real(dp), intent(out) :: ec,vc
!Local variables ------------------------------
 logical :: islambda
 real(dp) :: aa_pbe,alb,alb_pow,arg_log,daadec,decdrho,den,df,dhdaa
 real(dp) :: dhdrho,dhdtt,div,div2,dttdrho,eps2,eps3,exp_pbe,fx,h_pbe,lg,pade,pow,q0
 real(dp) :: q1,q1p,rhothird,rs,sqr_rs,tt,xx
 real(dp), parameter :: aa = 0.031091_dp,a1 = 0.21370_dp
 real(dp), parameter :: b1 = 7.5957_dp,b2 = 3.5876_dp
 real(dp), parameter :: b3 = 1.6382_dp,b4 = 0.49294_dp
 real(dp), parameter :: a = 0.47808102_dp,b = 0.84449703_dp
 real(dp), parameter :: c = 1.30089155_dp,d = 0.02949437_dp,beta = 1.34835105_dp
 real(dp), parameter :: kf_fac = (3.0_dp*(pi**2))**third
 real(dp), parameter :: rs_fac = (3.0_dp/(4.0_dp*pi))**third
 real(dp), parameter :: beta_pbe = 0.066725_dp,gamma_pbe = (1.0_dp-log(2.0_dp))/(pi**2)
!************************************************************************

 if (abs(rho) < tol30) then
   ec = 0.0_dp
   vc = 0.0_dp
   return
 end if

 rhothird = rho**third
 rs = rs_fac / (rhothird*eps) ! Scaling law for r_s
 sqr_rs = sqrt(rs)

 eps2 = eps * eps

 ! LDA correlation for rescaled Coulomb potential (PW91 parametrization)
 q0  = -2.0_dp * aa * (1.0_dp+a1*rs)
 q1  = 2.0_dp * aa * (b1*sqr_rs+b2*rs+b3*rs*sqr_rs+b4*rs*rs)
 q1p = aa * (b1/sqr_rs+2._dp*b2+3._dp*b3*sqr_rs+4._dp*b4*rs)
 den = 1.0_dp / (q1*q1+q1)
 lg  = -log(q1*q1*den)
 ec  = q0 * lg
 vc  = -2.0_dp*aa*a1*lg - q0*q1p*den
 vc  = ec - rs*vc/3.0_dp

 islambda = (abs(lambda) > tol10)

 ! Correction factor for Yukawa potential
 ! Parametrization based on data from Savin, "Beyond the Kohn-Sham Determinant"
 pow = c + d*log(1.0_dp+rs)
 alb = a*lambda*eps*(rs**b) ; alb_pow = alb**pow  ! Scaling law for lambda
 fx  = (1.0_dp+alb_pow)**(-beta-1.0_dp)
 if (islambda) then
   df = -beta * fx * alb_pow * (pow*b/rs+d*log(alb)/(1.0_dp+rs))
 else
   df = zero
 end if
 fx = fx * (1.0_dp+alb_pow)

 vc = (vc*fx-ec*rs*df/3.0_dp) / eps2 ! Scaling law
 ec = ec * fx / eps2

 if (ixc == 7 .or. ixc == -1012) return

 eps3 = eps2 * eps
 tt = grad * pi / ((rho**2)*16.0_dp*kf_fac*rhothird)
 exp_pbe = exp(-ec*eps2/gamma_pbe)
 if (abs(exp_pbe-1.0_dp) < tol30) return

 aa_pbe = eps * beta_pbe / (gamma_pbe*(exp_pbe-1.0_dp))
 daadec = eps3 * beta_pbe * exp_pbe / (gamma_pbe*(exp_pbe-1.0_dp))**2
 decdrho = (vc-ec) / rho
 xx = aa_pbe * tt
 div = 1.0_dp + xx + xx**2
 div2 = div * div
 pade = (1.0_dp+xx) / div
 arg_log = 1.0_dp + eps*beta_pbe*tt*pade/gamma_pbe
 h_pbe = gamma_pbe * log(arg_log) / eps2
 dhdaa = -beta_pbe * (tt**2) * xx * (2.0_dp+xx) / (eps*div2*arg_log)
 dttdrho = -7.0_dp * tt / (3.0_dp*rho)
 dhdtt = beta_pbe * (pade-(xx**2)*(2.0_dp+xx)/div2) / (arg_log*eps)
 dhdrho = dhdtt*dttdrho + dhdaa*daadec*decdrho
 vc = vc + h_pbe + rho*dhdrho
 ec = ec + h_pbe

end subroutine correlation_yukawa
!!***

END MODULE m_paw_exactDC
!!***


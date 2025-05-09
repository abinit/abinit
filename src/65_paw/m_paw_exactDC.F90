!!****m* m_paw_exactDC/m_paw_exactDC
!! NAME
!!  m_paw_exactDC
!!
!! FUNCTION
!!  This module contains several routines related to the "exact" double counting.
!!
!! COPYRIGHT
!! Copyright (C) 2020 Kristjan Haule
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

MODULE m_paw_exactDC

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private

 public :: compute_exactDC
 public :: grule

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_exactDC/compute_exactDC
!! NAME
!! compute_exactDC
!!
!! FUNCTION
!!
!! Compute the "exact" double counting.
!! See Physical review letters, Haule, K. (2015), 115(19), 196403  for formula.
!!
!! INPUTS
!!  lpawu = angular momentum for correlated species
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  occ(2*lpawu+1,2*lpawu+1) = occupation matrix (summed over spins) in
!!    the real spherical harmonics basis, with the convention
!!    occ(i,j) = <c_j^dagger c_i>
!!
!! OUTPUT
!!  vdc(2*lpawu+1,2*lpawu+1) = double counting potential
!!  edc = double counting energy
!!  edcdc = integral of Vdc(r)*rho_loc(r)
!!
!! SOURCE

 subroutine compute_exactDC(lpawu,pawtab,pawrad,occ,vdc,edc,edcdc)

 use m_pawtab, only : pawtab_type
 use m_pawrad, only : pawrad_type,simp_gen
 use m_paw_sphharm, only : ylmc

!Arguments ------------------------------------
 integer, intent(in) :: lpawu
 type(pawtab_type), intent(in) :: pawtab
 type(pawrad_type), intent(in) :: pawrad
 complex(dpc), intent(in) :: occ(2*lpawu+1,2*lpawu+1)
 complex(dpc), intent(inout) :: vdc(2*lpawu+1,2*lpawu+1)
 real(dp), intent(out) :: edc,edcdc
!Local variables ------------------------------
 integer :: i,ir,j,k,l,ln,m,m1,mm,minusm,meshsz,ndim
 real(dp) :: excint,phi,rho,rs1,sexc,svxc,theta,vxcf
 real(dp) :: ecc(1),exx(1),kcart(3),rs1s(1),vcc(1),vxx(1)
 complex(dpc) :: ylm_c
 real(dp), allocatable :: exc(:,:,:),exci(:),phis(:),rho_angle(:,:)
 real(dp), allocatable :: thetas(:),tweights(:),vxc(:,:,:),vxci(:,:,:),ylm(:,:,:)
 !************************************************************************

 ! VERY IMPORTANT: This routine assumes that we work in the real spherical harmonic basis.
 ! If you want to generalize that to the complex case for some reason, the formulas
 ! below are not valid, so you would need to add some complex conjugates.
 ! Also, we work with the convention occ(i,j) = <c_j^dagger c_i>.

 ln = 13  ! Size of angular mesh, same as eDMFT
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
 call tfpoint(thetas(:),phis(:),tweights(:),ln)

 ! Compute the real spherical harmonics on the angular mesh
 do m=0,lpawu
   mm = m + lpawu + 1
   minusm = - m + lpawu + 1
   do j=1,2*ln+1
     phi = phis(j)
     do i=1,ln+1
       theta = thetas(i)
       ! The definition of theta and phi is swapped in ylmc
       kcart(1) = cos(phi) * sin(theta)
       kcart(2) = sin(phi) * sin(theta)
       kcart(3) = cos(theta)
       ylm_c = ylmc(lpawu,m,kcart(:))
       ! Convert complex harmonics to real harmonics
       if (m == 0) then
         ylm(i,j,mm) = dble(ylm_c)
       else
         ylm(i,j,mm) = (-one)**m * sqrt2 * dble(ylm_c)
         ylm(i,j,minusm) = (-one)**m * sqrt2 * aimag(ylm_c)
       end if ! m=0
     end do ! i
   end do ! j
 end do ! m

 rho_angle(:,:) = zero

 do m1=1,ndim
   do m=1,ndim
     ! In the complex case, you should use ylm(:,:,m)*conjg(ylm(:,:,m1))
     rho_angle(:,:) = rho_angle(:,:) + dble(occ(m,m1))*ylm(:,:,m)*ylm(:,:,m1)
   end do ! m
 end do ! m1

 do ir=1,meshsz
   do j=1,2*ln+1
     do i=1,ln+1
       rho = rho_angle(i,j) * pawtab%proj2(ir) / (max(pawrad%rad(ir),epsilon(one))**2)
       rs1 = (four_pi * rho / three)**(third)
       rs1s(1) = rs1
       call ExchangeLDA(exx(:),vxx(:),rs1s(:),pawtab%lambda,1)
       exc(i,j,ir) = exx(1) * half / pawtab%eps ! convert from Rydberg to Hartree
       vxc(i,j,ir) = vxx(1) * half / pawtab%eps
       rs1s(1) = one / (max(rs1,epsilon(one)))
       call CorrLDA_2(ecc(:),vcc(:),rs1s(:),pawtab%lambda,pawtab%eps,1)
       exc(i,j,ir) = exc(i,j,ir) + ecc(1)*half
       vxc(i,j,ir) = vxc(i,j,ir) + vcc(1)*half
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
         svxc = svxc + two_pi*sum(ylm(:,j,m)*ylm(:,j,m1)*tweights(:)* &
              & vxc(:,j,ir))/dble(2*ln+1)
         sexc = sexc + two_pi*sum(ylm(:,j,m)*ylm(:,j,m1)*tweights(:)* &
              & exc(:,j,ir))*dble(occ(m,m1))/dble(2*ln+1)
       end do ! j
       vxci(ir,m,m1) = svxc
       exci(ir) = exci(ir) + sexc
     end do ! ir
   end do ! m
 end do ! m1

 ! Integrals over r

 call simp_gen(excint,exci(:)*pawtab%proj2(:),pawrad,r_for_intg=pawrad%rad(meshsz))
 edc = edc + excint
 edcdc = zero

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

 end subroutine compute_exactDC
!!***

SUBROUTINE tfpoint(thetas,phis,tweights,ln)
  !************************************************!
  !     GENERATES POINTS ON A UNIT SPHERE          !
  !************************************************!
  IMPLICIT NONE
  !
  INTEGER, intent(in)  :: ln
  REAL*8, intent(out)  :: thetas(ln+1), phis(2*ln+1), tweights(ln+1)! spt(2,(ln+1)*(2*ln+1)),weight((ln+1)*(2*ln+1))
  !
  REAL*8 :: dsum
  REAL*8 :: fake_shift, phi, pi
  INTEGER :: i, j, lg
  !
  PI=4.D0*ATAN(1.D0)
  lg = ln+1
  !
  CALL GRULE(lg,thetas,tweights) ! Gauss points in the interval [1...0]
  !
  do i=1,(lg+1)/2                ! adding symmetrical points to [0,...-1]
     j = lg+1-i
     thetas(j)=-thetas(i)
     tweights(j)=tweights(i)
  enddo
  do i=1,lg
     thetas(i) = dacos(thetas(i))
  enddo
  !
  fake_shift=EXP(-4.D0)  ! small shift such that we do not start at phi=0
  DO i=0,2*lg-2       ! 2*l-1 points in phi direction
     PHI=PI*(2*dble(i)/dble(2*lg-1) + fake_shift)
     phis(i+1)=PHI
  ENDDO
  !
  dsum = sum(tweights)
  tweights(:) = tweights(:)*2./dsum
  !write(6,*)'Gauss-Legendre grid of ',lg,'x',2*lg-1
  return
end SUBROUTINE tfpoint

SUBROUTINE GRULE(N,X,W)
  IMPLICIT NONE
  !
  !     DETERMINES THE (N+1)/2 NONNEGATIVE POINTS X(I) AND
  !     THE CORRESPONDING WEIGHTS W(I) OF THE N-POINT
  !     GAUSS-LEGENDRE INTEGRATION RULE, NORMALIZED TO THE
  !     INTERVAL \-1,1\. THE X(I) APPEAR IN DESCENDING ORDER.
  !
  !     THIS ROUTINE IS FROM 'METHODS OF NUMERICAL INTEGRATION',
  !     P.J. DAVIS AND P. RABINOWITZ, PAGE 369.
  !
  INTEGER, intent(in) :: N
  REAL*8, intent(out) :: X(*), W(*)
  !
  REAL*8  :: D1, d2pn, d3pn, d4pn, den, dp, dpn, e1, fx, h, p, pi, pk, pkm1, pkp1, t, t1, u, v, x0
  INTEGER :: i, it, k, m
  !
  PI=4.D0*ATAN(1.D0)
  M=(N+1)/2
  E1=N*(N+1)
  DO I=1,M         ! 1
     T=(4*I-1)*PI/(4*N+2)
     X0=(1.D0-(1.D0-1.D0/N)/(8.D0*N*N))*COS(T)
     !--->    ITERATE ON THE VALUE  (M.W. JAN. 1982)
     DO IT=1,3     ! 2
        PKM1=1.D0
        PK=X0
        DO K=2,N   ! 3
           T1=X0*PK
           PKP1=T1-PKM1-(T1-PKM1)/K+T1
           PKM1=PK
           PK=PKP1
        ENDDO      ! 3
        DEN=1.D0-X0*X0
        D1=N*(PKM1-X0*PK)
        DPN=D1/DEN
        D2PN=(2.D0*X0*DPN-E1*PK)/DEN
        D3PN=(4.D0*X0*D2PN+(2.D0-E1)*DPN)/DEN
        D4PN=(6.D0*X0*D3PN+(6.D0-E1)*D2PN)/DEN
        U=PK/DPN
        V=D2PN/DPN
        H=-U*(1.D0+.5D0*U*(V+U*(V*V-U*D3PN/(3.D0*DPN))))
        P=PK+H*(DPN+.5D0*H*(D2PN+H/3.D0*(D3PN+.25D0*H*D4PN)))
        DP=DPN+H*(D2PN+.5D0*H*(D3PN+H*D4PN/3.D0))
        H=H-P/DP
        X0=X0+H
     ENDDO
     X(I)=X0
     FX=D1-H*E1*(PK+.5D0*H*(DPN+H/3.D0*(D2PN+.25D0*H*(D3PN+.2D0*H*D4PN))))
     W(I)=2.D0*(1.D0-X(I)*X(I))/(FX*FX)
  ENDDO
  IF(M+M.GT.N) X(M)=0.D0
  RETURN
END SUBROUTINE GRULE

subroutine ExchangeLDA(Ex, Vx, rs_1, lambda, N)
  IMPLICIT NONE
  REAL*8, intent(in) :: lambda
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs_1(N)
  REAL*8, intent(out):: Ex(N), Vx(N)
  ! locals
  REAL*8 :: kf_rs, c0, pi
  REAL*8 :: xs(N), dEx(N)
  !
  pi = 4.*atan(1.)
  !
  kf_rs = (9*pi/4.)**(1./3.)
  c0 = (3./(2.*pi))*(9.*pi/4.)**(1./3.)
  if (abs(lambda).gt.1e-10) then
     !print *, 'lambda too small'
     xs(:) = rs_1(:)*(kf_rs/lambda)       ! kf/(rs*lambda)
     CALL fexchange(xs,Ex,dEx,N)
  else
     !print *, 'lambda zero', lambda, abs(lambda).gt.1e-10
     xs(:) = rs_1(:)*kf_rs                ! kf/(rs*lambda)
     Ex(:)=1.0
     dEx(:)=0.0
  endif
  Ex(:) = -(c0*rs_1)*Ex(:)
  Vx(:) = 4./3.*Ex(:) - (c0/3.)*rs_1(:)*xs(:)*dEx(:)
end subroutine ExchangeLDA

SUBROUTINE fexchange(xs,ex,dex,N)
  !* Evaluating a function (and its derivative) for exchange energy, which has a formula
  !*  f(x) = 1-1/(6*x^2)-4/(3*x)*atan(2*x)+(1+1/(12*x^2))/(2*x^2)*log(1+4*x^2)
  !* df/dx = 2/(3*x^3) + 4/(3*x^2)*atan(2*x) - (1+6*x^2)/(6*x^5)*log(1+4*x^2)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: xs(N)
  REAL*8, intent(out):: ex(N), dex(N)
  ! locals
  INTEGER :: i
  integer, volatile :: iCopy
  REAL*8 :: x, x2, x3, at2, lg2
  do i=1,N
     x = xs(i)
     x2 = x*x
     x3 = x2*x
     iCopy = i  ! A fix for intel compiler bug, so that it does not optimize out the loop
     if (x<0.01) then
        ex(i) = 4.*x2/9.*(1-6*x2/5.)
        dex(i) = 8.*x/9.*(1-12*x2/5.)
     else
        at2 = atan(2*x)*4./(3.*x)
        lg2 = log(1+4.*x2)/(2.*x2)
        ex(i) = 1-1/(6.*x2) - at2+ (1.+1/(12.*x2))*lg2
        dex(i) = 2./(3.*x3) + at2/x - (1.+6*x2)/(3.*x3)*lg2
     end if
  end do
END SUBROUTINE fexchange

subroutine CorrLDA_2(Ec, Vc, rs, lambda, eps, N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: lambda, eps
  REAL*8, intent(out):: Ec(N), Vc(N)
  !locals
  INTEGER :: i
  REAL*8 :: fn_l(N), qn_l(N), fn_e(N), qn_e(N), fn(N), qn(N)

  do i=1,N
     CALL CorLDA(Ec(i),Vc(i),rs(i))
  enddo
  CALL EcVc_reduce_yw_2(rs,lambda,fn_l,qn_l,N)
  !print *, 'l=', lambda, 'fn_l=', fn_l
  !print *, 'rs=', rs, 'qn_l=', qn_l
  CALL EcVc_reduce_di_2(rs,eps,fn_e,qn_e,N)
  !print *, 'e=', eps, 'fn_e=', fn_e
  !print *, 'rs=', rs, 'qn_e=', qn_e

  fn(:) = fn_l(:)*fn_e(:)
  qn(:) = qn_l(:)*fn_e(:)+fn_l(:)*qn_e(:)

  Vc(:) = Vc(:)*fn(:) + Ec(:)*qn(:)
  Ec(:) = Ec(:)*fn(:)
end subroutine CorrLDA_2

subroutine CorLDA(Ec,Vc,rs)
  IMPLICIT NONE
  !  UNITS OF Rydberg
  !  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
  !  INPUT: SEITZ RADIUS (rs)
  !  OUTPUT: CORRELATION ENERGY PER ELECTRON (Ec) and POTENTIALS (Vc)
  !
  REAL*8, intent(out):: Ec, Vc
  REAL*8, intent(in) :: rs
  ! locals
  REAL*8 :: H2Ry, Ecrs
  H2Ry = 2.0    ! Conversion from Hartree to Rydberg
  CALL Gcor(Ec, Ecrs,    rs, 0.0310907D0,  0.21370D0,  7.5957D0, 3.5876D0, 1.6382D0, 0.49294D0, 1)
  Vc = Ec - rs*Ecrs/3.
  Vc = Vc*H2Ry  ! From Hartree to Rydbergs
  Ec = Ec*H2Ry  ! From Hartree to Rydbergs
  !if (ISNAN(GG)) then
  !print *, 'rs=', rs, 'RS12=', RS12, 'RS32=', RS32, 'RSP=', RSP, 'Q1=', Q1, 'Q2=', Q2, 'GG=', GG, 'GGRS=', GGRS
  !endif
  !print *, 'rs=', rs, Ec, Vc
end subroutine CorLDA

SUBROUTINE  Gcor(GG,GGrs,rs, A,A1,B1,B2,B3,B4,P)
  IMPLICIT NONE
  REAL*8, intent(in)  :: rs, A, A1, B1, B2, B3, B4
  INTEGER, intent(in) :: P
  REAL*8, intent(out) :: GG, GGrs
  ! locals
  REAL*8 :: P1, Q0, RS12, RS32, RSP, Q1, Q2, Q3
  P1 = P + 1.
  Q0 = -2.*A*(1.+A1*rs)
  RS12 = sqrt(rs)
  RS32 = RS12**3
  RSP = rs**P
  Q1 = 2.*A*(B1*RS12+B2*rs+B3*RS32+B4*rs*RSP)
  Q2 = log(1.+1./Q1)
  GG = Q0*Q2
  Q3 = A*(B1/RS12+2.*B2+3.*B3*RS12+2.*B4*P1*RSP)
  GGRS = -2.*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
END SUBROUTINE Gcor

subroutine EcVc_reduce_yw_2(rs,lambda,fn,qn,N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: lambda
  REAL*8, intent(out):: fn(N), qn(N)
  ! locals
  INTEGER:: i, m
  REAL*8 :: lmrs, dlms, te, das
  REAL*8 :: C(7,6)
  REAL*8 :: an(6)
  C = Reshape((/0.15805009, -0.77391602, 1.23971169, -1.04865383, 0.47809619, -0.11057964, 0.01016968,&
               -0.306851, -0.77296572, 0.8791705, -0.69185034, 0.33779654, -0.08858483, 0.00935635,&
               0.13215843, -0.2776552, 0.45727548, -0.31469164, 0.10787374, -0.01661214, 0.0007591,&
               -0.03086548, 0.0549528, -0.07252823, 0.04177618, -0.01084882, 0.00062192, 0.0001177,&
               0.00273230889, -0.00357007233, 0.00425309814, -0.00198811211, 0.000233761378, 0.000106803015, -2.50612307e-05,&
               -9.28530649e-05, 8.09009085e-05, -9.43747991e-05, 3.89520548e-05, -3.10149723e-07, &
               -4.23041605e-06, 8.02291467e-07/),(/7,6/))
  do m=1,6
     an(m) = lambda * (c(1,m) + lambda * (c(2,m) + lambda * (c(3,m) + lambda * (c(4,m) + lambda * (c(5,m) + lambda * (c(6,m) &
             + lambda*c(7,m)))))))
  enddo
  !print *, an

  do i=1,N
     te = exp(an(1)+rs(i)*(an(2)+rs(i)*(an(3)+rs(i)*(an(4)+rs(i)*(an(5)+rs(i)*an(6))))))
     lmrs = 0.008 - 0.00112 * rs(i)**2
     dlms = -2*0.00112*rs(i)
     fn(i) = te*(1-lmrs) + lmrs
     das = rs(i)*(an(2) + rs(i)*(2*an(3) + rs(i)*(3*an(4) + rs(i)*(4*an(5) + rs(i)*5*an(6)))))
     qn(i) = -1./3.*te*(1-lmrs)*das - 1./3.*rs(i)*(1-te)*dlms
  enddo
  ! Ec = Ev * fn
  ! Vc = Vc * fn + Ec * qn
end subroutine EcVc_reduce_yw_2

subroutine EcVc_reduce_di_2(rs,eps,fn,qn,N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: eps
  REAL*8, intent(out):: fn(N), qn(N)
  ! locals
  INTEGER:: i
  REAL*8 :: a, b, d, eps_a, eps_d, eps_abd, b2r_9, r_da_dr, r_db_dr, r_dd_dr
  REAL*8 :: ca(3), cb(3), cd(4)
  ! The coefficients a, b, d depend on rs, and here is their more precise fit:
  ca = (/1.74596971, -0.0892907,   0.00658866/)
  cb = (/ 1.63289109,  1.15291480, 0.149402/)
  cd = (/3.64370598, 0.03636027, -0.03886317, 0.00693599/)

  do i=1,N
     if (rs(i) .gt. 40) then
        qn(i)=0.0
        CYCLE
     endif
     a = ca(1) + rs(i) * ca(2) + rs(i)**2 * ca(3)
     b = 0.001*(cb(1)*sqrt(rs(i))+cb(2)*rs(i))/(1+ (cb(3)*rs(i))**9)
     d = cd(1) + rs(i) * cd(2) + rs(i)**2 * cd(3) + rs(i)**3 * cd(4)
     eps_a = eps**a
     eps_d = eps**d
     eps_abd = eps_a+b*eps_d
     b2r_9 = (cb(3)*rs(i))**9
     fn(i) = (1.+b)/eps_abd
     r_da_dr = rs(i) * ca(2) + 2*rs(i)**2 * ca(3)
     r_dd_dr = rs(i) * cd(2) + 2*rs(i)**2 * cd(3) + 3*rs(i)**3 * cd(4)
     r_db_dr = 0.001*(cb(1)*sqrt(rs(i))*(0.5 - 17./2.*b2r_9)+cb(2)*rs(i)*(1-b2r_9*8))/(1+ b2r_9)**2
     qn(i) = -1./3.*(r_db_dr*(eps_a-eps_d)-(eps_a*r_da_dr + b*eps_d*r_dd_dr)*(1+b)*log(eps))/eps_abd**2
  enddo
  ! Ec = Ev * fn
  ! Vc = Vc * fn + Ec * qn
end subroutine EcVc_reduce_di_2

END MODULE m_paw_exactDC
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/nfourier
!! NAME
!!  nfourier
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_lib_four

 use defs_basis
 use m_errors
 use m_abicore

contains

! This routine contains direct and inverse fourier transformation
! It is a modification of a routine of the GNU GPL
! code available on http://dmft.rutgers.edu/ and
! described in the RMP 2006 paper written by
! G.Kotliar, S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti
!=======+=========+=========+=========+=========+=========+=========+=$
!       TYPE   : SUBROUTINE
!       PROGRAM: nfourier
!       PURPOSE: fourier-transform the natural-spline interpolation
!                of function Green(tau)
!                calculate function Green(omega)
!       I/O    :
!       VERSION: 2-16-92
!                29-Nov-95 removal of minimal bug concerning
!                          DIMENSION of rindata
!       COMMENT: cf J. Stoer R. Bulirsch, Introduction to numerical
!                analysis (Springer, New York, 1980)
!=======+=========+=========+=========+=========+=========+=========+=$
!
      SUBROUTINE nfourier(rindata,coutdata,iflag,Iwmax,L,Beta)

!include 'param.dat'
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER(I-N)
       integer  Iwmax,L,iflag
       DIMENSION rindata(L)
       DIMENSION rincopy(L+1),a(L),b(L),c(L),d(L),u(L+1), q(L+1),XM(L+1)
       complex*16 :: coutdata(Iwmax+1)
       complex*16 cdummy,explus,ex
       xpi = ACOS(-One)
       delta = Beta/L
       DO i = 1,L
          rincopy(i) = rindata(i)
       ENDDO
       if(iflag==1) then
         rincopy(L+1) = -1-rindata(1)
       else
         rincopy(L+1) = -rindata(1)
       endif
!       Three = Two+One
!       six = Two*Three

!c
!c     spline interpolation:  the spline is given by
!c     G(tau) = a(i) + b(i) (tau-tau_i) + c(i) ( )^2 + d(i) ( )^3
!c     The following formulas are taken directly from  Stoer and
!c     Bulirsch p. 102
!c
       q(1) = Zero
       u(1) = Zero
       DO k = 2,L
          p = q(k-1)/Two+Two
          q(k)=-One/Two/p
          u(k)=Three/delta**2*(rincopy(k+1)+rincopy(k-1)-Two*rincopy(k))
          u(k)=(u(k)-u(k-1)/Two)/p
       ENDDO
       XM(L+1) = 0
       DO k = L,1,-1
          XM(k) = q(k)*XM(k+1)+u(k)
       ENDDO
!c
!c     The following formulas are taken directly from  Stoer and
!c     Bulirsch p. 98
!c
       DO j = 1, L
          a(j) = rincopy(j)
          c(j) = XM(j)/Two
          b(j) = (rincopy(j+1)-rincopy(j))/delta - &
     &       (Two*XM(j)+XM(j+1))*delta/6.
          d(j) = (XM(j+1)-XM(j))/(6.*delta)
       ENDDO

!c
!c     The Spline multiplied by the exponential can now be exlicitely
!c     integrated. The following formulas were obtained using
!c     MATHEMATICA
!c
        DO i = 0,Iwmax
           om = (Two*(i)+One)*xpi/Beta
           coutdata(i+1) = czero
           DO j = 1,L
              cdummy = j_dpc*om*delta*j
              explus = exp(cdummy)
              cdummy = j_dpc*om*delta*(j-1)
              ex = exp(cdummy)
              coutdata(i+1) = coutdata(i+1) + explus*(&
     &         ( -six* d(j) )/om**4 + &
     &         ( Two*j_dpc*c(j) + six*delta*j_dpc*d(j)  )/om**3 +&
     &         ( b(j)+ Two*delta*c(j)+ three*delta**2*d(j) )/om**2 +&
     &         (- j_dpc*a(j) - delta*j_dpc*b(j) - delta**2*j_dpc*c(j) -&
     &         delta**3*j_dpc*d(j))/om)

              coutdata(i+1) = coutdata(i+1) + ex*(&
     &        six*d(j)/om**4 - Two*j_dpc*c(j)/om**3 &
     &        -b(j)/om**2 + j_dpc*a(j)/om)
           ENDDO
        ENDDO
        end subroutine nfourier
!!***

! This routine contains direct and inverse fourier transformation
! It is a modification of a routine of the GNU GPL
! code available on http://dmft.rutgers.edu/ and
! described in the RMP 2006 paper written by
! G.Kotliar, S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti
!=======+=========+=========+=========+=========+=========+=========+=$
!       TYPE   : SUBROUTINE
!       PROGRAM: nfourier
!       PURPOSE: fourier-transform the natural-spline interpolation
!                of function Green(tau)
!                calculate function Green(omega)
!       I/O    :
!       VERSION: 2-16-92
!                29-Nov-95 removal of minimal bug concerning
!                          DIMENSION of rindata
!       COMMENT: cf J. Stoer R. Bulirsch, Introduction to numerical
!                analysis (Springer, New York, 1980)
!=======+=========+=========+=========+=========+=========+=========+=$
!
      SUBROUTINE nfourier2(rindata,coutdata,iflag,om,L,Beta)
!       include 'param.dat'
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER(I-N)
       integer  L,iflag
       DIMENSION rindata(L)
       DIMENSION rincopy(L+1),a(L),b(L),c(L),d(L),u(L+1), q(L+1),XM(L+1)
       complex*16 :: coutdata
       real*8 :: om
       complex*16 cdummy,explus,ex
       xpi = ACOS(-One)
       delta = Beta/L
       DO i = 1,L
          rincopy(i) = rindata(i)
       ENDDO
       if(iflag==1 .and. L.ge.1) then
         rincopy(L+1) = -1-rindata(1)
       elseif(iflag==0 .and. L.ge.1) then
         rincopy(L+1) = -rindata(1)
       else
         write(std_out,*) "Warning : Check nfourier2"
       endif
!       Three = Two+One
!       six = Two*Three

!c
!c     spline interpolation:  the spline is given by
!c     G(tau) = a(i) + b(i) (tau-tau_i) + c(i) ( )^2 + d(i) ( )^3
!c     The following formulas are taken directly from  Stoer and
!c     Bulirsch p. 102
!c
       q(1) = Zero
       u(1) = Zero
       DO k = 2,L
          p = q(k-1)/Two+Two
          q(k)=-One/Two/p
          u(k)=Three/delta**2*(rincopy(k+1)+rincopy(k-1)-Two*rincopy(k))
          u(k)=(u(k)-u(k-1)/Two)/p
       ENDDO
       XM(L+1) = 0
       DO k = L,1,-1
          XM(k) = q(k)*XM(k+1)+u(k)
       ENDDO
!c
!c     The following formulas are taken directly from  Stoer and
!c     Bulirsch p. 98
!c
       DO j = 1, L
          a(j) = rincopy(j)
          c(j) = XM(j)/Two
          b(j) = (rincopy(j+1)-rincopy(j))/delta - &
     &       (Two*XM(j)+XM(j+1))*delta/6.
          d(j) = (XM(j+1)-XM(j))/(6.*delta)
       ENDDO

!c
!c     The Spline multiplied by the exponential can now be exlicitely
!c     integrated. The following formulas were obtained using
!c     MATHEMATICA
!c
        coutdata = czero
        DO j = 1,L
           cdummy = j_dpc*om*delta*j
           explus = exp(cdummy)
           cdummy = j_dpc*om*delta*(j-1)
           ex = exp(cdummy)
           coutdata = coutdata + explus*(&
     &      ( -six* d(j) )/om**4 + &
     &      ( Two*j_dpc*c(j) + six*delta*j_dpc*d(j)  )/om**3 +&
     &      ( b(j)+ Two*delta*c(j)+ three*delta**2*d(j) )/om**2 +&
     &      (- j_dpc*a(j) - delta*j_dpc*b(j) - delta**2*j_dpc*c(j) -&
     &         delta**3*j_dpc*d(j))/om)

           coutdata = coutdata + ex*(&
     &     six*d(j)/om**4 - Two*j_dpc*c(j)/om**3 &
     &     -b(j)/om**2 + j_dpc*a(j)/om)
        ENDDO
        end subroutine nfourier2
!=======+=========+=========+=========+=========+=========+=========+=$
!       TYPE   : SUBROUTINE
!       PROGRAM: invfourier
!       PURPOSE: inverse fourier transform
!                Greent, Greenw use physical definition
!                Greent(i) = G((i-1)*deltau) for i = 1,...,L
!                Greenw(n) = G(i w_n), for n = 0,L/2-1
!                       w_n = (2*n+1)pi/beta
!                Symmetry property:
!                G(iw_(-n) = G(iw_(n-1))*
!                coupled to the impurity
!       I/O    :
!       VERSION: 6-16-92
!       COMMENT:
!=======+=========+=========+=========+=========+=========+=========+=$
!
       SUBROUTINE invfourier(cindata,routdata,Iwmax,L,iflag,beta)

!       include 'param.dat'
       implicit none
       integer, intent(in) :: Iwmax
       complex*16, intent(in) :: cindata(1:Iwmax)    !vz_d
       integer, intent(in) :: L
       complex*16, intent(inout) :: routdata(1:L)    !vz_d
       integer, intent(in) :: iflag
       double precision, intent(in) :: beta

       double precision :: xpi
       double precision :: tau
       double precision :: om
       complex*16 :: cdummy,dummy
       integer :: i,j

       xpi = ACOS(-One)
       DO 1 i = 1,L
       routdata(i) = Zero
          tau = (i-1)*beta/real(L)
          DO 2 j = 1,Iwmax
!               om = mod((2*(j)+One)*xpi/Beta*tau,2*xpi)
                om =    ((2*(j)-One)*xpi/Beta*tau)
               cdummy = CMPLX(Zero,om)
               dummy = cindata(j)*exp(-cdummy)
           routdata(i) = routdata(i)+Two/beta*dummy
2         CONTINUE
!           write(std_out,*) "FT",i,routdata(i)
1      CONTINUE
!c
!c     special treatment for tau = 0
!c
       if(iflag==1 .and. L.ge.1 ) then
         routdata(1) = -One/Two+routdata(1)
       endif
       END SUBROUTINE invfourier

END MODULE m_lib_four
!!***

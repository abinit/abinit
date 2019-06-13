!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_simtet
!! NAME
!! m_simtet
!!
!! FUNCTION
!!
!! COPYRIGHT
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

module m_simtet

  use defs_basis
  use m_errors, only: unused_var

 implicit none

 private

 public :: sim0onei
 public :: sim0twoi
!!***

 contains

!C file: sim0onei.f
!C date: 2011-04-01
!C who:  S.Kaprzyk
!C what: Complex linear form integral over standard terahedron
!C ------------------------------------------------------------
!C -         *                                                -
!C -        *                         1                       -
!C - SIM0=6 * dt1*dt2*dt3 ----------------------------------  -
!C -        *             verm(4)+(verm(1)-verm(4))*t1+ ..... -
!C -       *                                                  -
!C -   0<t1+t2+t3<1                                           -
!C ------------------------------------------------------------
SUBROUTINE SIM0ONEI(SIM0, SIM0I, VERM)

      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX SIM0, SIM0I
      DOUBLE COMPLEX VERM(4)
!C
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(4),LEPS(3)
!c
      DOUBLE COMPLEX   VERM2D(3)
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3), AL(4)
      INTEGER          I, I1, I2, I3, K, N
      DOUBLE PRECISION ZERO, EPS, SMALL
      DATA             EPS/1.0D-6/, SMALL/1.0D-5/
      DATA             ZERO/0.0D0/
!C
      DO 100 I = 1, 3
       IF (DIMAG(VERM(4))*DIMAG(VERM(I)).LT.ZERO) THEN
       WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,4)
 9010  FORMAT (' ***sim0onei: not signed ImgVERM()=',4(D13.6,1X))
!c      STOP ' ***sim0onei: '
       END IF
 100  CONTINUE
      DO 200 I = 1, 4
       AL(I) = CDABS(VERM(I))
 200  CONTINUE
      DO 300 N = 1, 4
       IF (AL(N).LT.SMALL) THEN
        I1 = MOD(N+0,4) + 1
        I2 = MOD(N+1,4) + 1
        I3 = MOD(N+2,4) + 1
        VERM2D(1) = VERM(I1)
        VERM2D(2) = VERM(I2)
        VERM2D(3) = VERM(I3)
        CALL S2D0ONEI(SIM,SIMI,VERM2D)
        SIM0 =  3*SIM/2
        SIM0I = SIMI + 1.75D0
        RETURN
       END IF
 300  CONTINUE
!c
      N = 4  !   N = 3, 2, 1
      DO 201 I = 1, 3
       IF (CDABS(VERM(I)).GT.CDABS(VERM(N))) N = I
 201  CONTINUE
!C
!C Here are 4 cases a), b), c), d); see, Fig.3.1
!C LONE(1..4) [w1<w2<w3<ONE]; [w1<w2<ONE<w3];
!C            [w1<ONE<w2<w3]; [ONE<w1<w2<w3]
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e]; [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
!*-ASIS
      CALL SIM0LEPS(VERM, N, U, AS, LONE, EPS, LEPS, iuerr)
!c U(1) = W(1); U(2) = W(2); U(3) = W(3)
      IF (LONE(1)) THEN
       CALL SIM1ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!c U(1) = W(1); U(2) = W(2); U(3) = CONE/W(3)
      IF (LONE(2)) THEN
       CALL SIM2ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!c U(1) = W(1); U(2) = CONE/W(2); U(3) = CONE/W(3)
      IF (LONE(3)) THEN
       CALL SIM3ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!c U(1) = CONE/W(1); U(2) = CONE/W(2); U(3) = CONE/W(3)
      IF (LONE(4)) THEN
       CALL SIM4ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
      SIM0 = 0.75D0*SIM/VERM(N)
      SIM0I = log(VERM(N)) + 0.25D0*SIMI
      SIM0I = SIM0I + (0.25D0,0.0D0)
      RETURN
      END SUBROUTINE SIM0ONEI
!C
!C -----------------------------------------------------------------
!C file: sim0twoi.f
!C date: 2011-04-02
!C who:  S.Kaprzyk
!C what: Complex linear form integrals over standard tetrahedron
!C ------------------(i=1,2,3)--------------------------------------
!C            *                                                    -
!C           *            t_i                                      -
!C  VERL(i)=6*dt1dt2dt3 ------------------------------------------ -
!C           *          VERM(4)+[VERM(1)-VERM(4)]*t1+...           -
!C          *                                                      -
!C     0<t1+t2+t3<1                                                -
!C -------------------(i=4)-----------------------------------------
!C            *                                                    -
!C           *           1- t1 -t2 - t3                            -
!C  VERL(4)=6*dt1dt2dt3 ------------------------------------------ -
!C           *          VERM(4)+[VERM(1)-VERM(4)]*t1+...           -
!C          *                                                      -
!C     0<t1+t2+t3<1                                                -
!C -----------------------------------------------------------------
SUBROUTINE   SIM0TWOI(VERL, VERLI, VERM)

      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX VERL(4), VERLI(4), VERM(4)
!C
      DOUBLE COMPLEX   SIM, SIMI
      DOUBLE COMPLEX   VERL2D(3), VERL2DI(3), VERM2D(3)
      LOGICAL          LONE(4), LEPS(3)
!c
      INTEGER          I, I1, I2, I3, K, N
      DOUBLE PRECISION AS(3), AL(4)
      DOUBLE COMPLEX   U(3)
      DOUBLE COMPLEX   CZERO
      DOUBLE PRECISION EPS, SMALL
      DOUBLE PRECISION ZERO
      DATA             EPS/1.0D-6/, SMALL/1.0D-5/
      DATA             ZERO/0.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      DO 50 I = 1, 3
        IF (DIMAG(VERM(4))*DIMAG(VERM(I)).LT.ZERO) THEN
        WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,4)
 9010   FORMAT (' ***sim0twoi: not signed ImgVERM()=',4(D13.6,1X))
!c        STOP ' ***sim0twoi: '
        END IF
 50    CONTINUE
       DO 200 I = 1, 4
       AL(I) = CDABS(VERM(I))
 200  CONTINUE
!c
      DO 100 N = 1, 4
       VERL(N) = CZERO
       VERLI(N) = CZERO
       I1 = MOD(N+0,4) + 1
       I2 = MOD(N+1,4) + 1
       I3 = MOD(N+2,4) + 1
       IF (AL(I1).LT.SMALL) THEN
        VERM2D(1) = VERM(I2)
        VERM2D(2) = VERM(I3)
        VERM2D(3) = VERM(N)
        CALL S2D0TWOI(VERL2D,VERL2DI,VERM2D)
        VERL(N) = VERL2D(3)
        VERLI(N) = VERL2DI(3)*3/4
        GO TO 100
       END IF
       IF (AL(I2).LT.SMALL) THEN
        VERM2D(1) = VERM(I3)
        VERM2D(2) = VERM(N)
        VERM2D(3) = VERM(I1)
        CALL S2D0TWOI(VERL2D,VERL2DI,VERM2D)
        VERL(N) = VERL2D(2)
        VERLI(N) = VERL2DI(2)*3/4
        GO TO 100
       END IF
       IF (AL(I3).LT.SMALL) THEN
        VERM2D(1) = VERM(N)
        VERM2D(2) = VERM(I1)
        VERM2D(3) = VERM(I2)
        CALL S2D0TWOI(VERL2D,VERL2DI,VERM2D)
        VERL(N) = VERL2D(1)
        VERLI(N) = VERL2DI(1)*3/4
        GO TO 100
       END IF
       IF (AL(N).LT.SMALL) THEN
        VERM2D(1) = VERM(I1)
        VERM2D(2) = VERM(I2)
        VERM2D(3) = VERM(I3)
        CALL S2D0ONEI(SIM,SIMI,VERM2D)
        VERL(N) = SIM/2
        VERLI(N) = SIMI/4 +1.75D0
        GO TO 100
       END IF
!C Here are 4 cases a), b), c), d); see, Fig.3.1
!C LONE(1..4) [w1<w2<w3<ONE]; [w1<w2<ONE<w3];
!C            [w1<ONE<w2<w3]; [ONE<w1<w2<w3]
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e]; [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
!*-ASIS
       CALL SIM0LEPS(VERM, N, U, AS, LONE, EPS, LEPS, iuerr)
!c U(1) = W(1); U(2) = W(2); U(3) = W(3)
       IF (LONE(1)) THEN
        CALL SIM1TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c U(1) = W(1); U(2) = W(2); U(3) = CONE/W(3)
       IF (LONE(2)) THEN
        CALL SIM2TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c U(1) = W(1); U(2) = CONE/W(2); U(3) = CONE/W(3)
       IF (LONE(3)) THEN
        CALL SIM3TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c U(1) = CONE/W(1); U(2) = CONE/W(2); U(3) = CONE/W(3)
       IF (LONE(4)) THEN
        CALL SIM4TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c
       VERL(N) = SIM/VERM(N)/8
       VERLI(N) = log(VERM(N))/4 + SIMI/32
       VERLI(N) = VERLI(N) + (0.25D0,0.0D0)
 100  CONTINUE ! end DO ... N = 1, 4
      RETURN
END SUBROUTINE SIM0TWOI
!C
!C ---------------------------------------------------------------------
!C file: sim0ur0.f
!C date: 2011-03-16
!C who:  S.Kaprzyk
!C what: R0(U)=Arth(w)/w=(ln(1+U)-ln(1-U))/(2*U)
!C what: R0(U)= 1 + U^2/3 + U^4/5 + ...
      DOUBLE COMPLEX FUNCTION SIM0UR0(U)

      IMPLICIT       NONE
      DOUBLE COMPLEX U
!c
!c      DOUBLE PRECISION DLAMCH
!c      EXTERNAL         DLAMCH
!c
      INTEGER          K
      DOUBLE COMPLEX   AV, BV, CV, UU
      DOUBLE COMPLEX   CONE
      DOUBLE PRECISION DEV, TOL
      DATA             TOL /1.0D-13/
!c
!c      TOL = DLAMCH('e')*10
      CONE = (1.0D0,0.0D0)
      UU = U*U
      IF (CDABS(U).GT.0.27D0) THEN
       SIM0UR0 = log((CONE+U)/(CONE-U))/(2*U)
!c       SIM0UR0 = (log(CONE+U)-log(CONE-U))/(2*U)
      ELSE
       CV = UU
       AV = CONE
       DO 50 K = 3, 13, 2
        AV = AV + CV/DBLE(K)
        CV = CV*UU
 50    CONTINUE
       DO 100 K = 15, 49, 2
        BV = CV/DBLE(K)
        AV = AV + BV
        DEV = CDABS(BV)
        IF (DEV.LE.TOL*CDABS(AV)) GO TO 150
        CV = CV*UU
 100   CONTINUE
       WRITE (std_out,9010) U, DEV, TOL
 9010  FORMAT ('***sim0ur0: U=',2(D15.8,1X),' DEV=',D14.7,' TOL=',D14.7)
       STOP '***sim0ur0: accuracy not reached'
 150   CONTINUE
       SIM0UR0 = AV
      END IF
!C
      RETURN
END FUNCTION SIM0UR0
!C
!C ----------------------------------------------------------------------
!C file: sim0ux0.f
!C date: 2010-08-17
!C who:  S.Kaprzyk
!C what: Complex version of X0(U)=(R0(U)/U-1-U^2/3-U^4/5)/U^6
!C       R0(U)=Arth(w)/w=Log((1+U)/(1-U))/(2*U)
      DOUBLE COMPLEX FUNCTION SIM0UX0(U)

      IMPLICIT       NONE
      DOUBLE COMPLEX U
!c
!c      DOUBLE PRECISION DLAMCH
!c      EXTERNAL         DLAMCH
!c
      INTEGER          K
      DOUBLE COMPLEX   AV, CV, UU
      DOUBLE COMPLEX   CONE
      DOUBLE PRECISION DEV, TOL
      DATA TOL /1.0D-13/
!c      TOL = DLAMCH('e')*10
      CONE = (1.0D0,0.0D0)
      UU = U*U
      IF (CDABS(U).GT.0.27D0) THEN
        SIM0UX0 =(log((CONE+U)/(CONE-U))/(2*U) - (CONE+UU*(CONE/3.0D0+UU*CONE/5.0D0)))/(UU*UU*UU)
      ELSE
       CV = UU
       AV = CONE/7.0D0
       DO 50 K = 9, 13, 2
        AV = AV + CV/DBLE(K)
        CV = CV*UU
 50    CONTINUE
       DO 100 K = 15, 49, 2
        AV = AV + CV/DBLE(K)
        DEV = CDABS(CV/AV)/DBLE(K)
        CV = CV*UU
        IF (DEV.LE.TOL) GO TO 150
 100   CONTINUE
       WRITE (std_out,9010) U, DEV, TOL
 9010  FORMAT ('***sim0ux0: U=',2(D15.8,1X), ' DEV=',D14.7,' TOL=',D14.7)
       STOP '***sim0ux0: accuracy not reached'
 150   CONTINUE
       SIM0UX0 = AV
      END IF
!C
      RETURN
      END FUNCTION SIM0UX0
!C
!C ----------------------------------------------------------------------
!C file: sim0udx0.f
!C date: 2010-08-17
!C who:  S.Kaprzyk
!C what: Complex version of DX0(u)=dX0/du
      SUBROUTINE SIM0UDX0(U, X0, DX0)

      IMPLICIT       NONE
      DOUBLE COMPLEX U, X0, DX0(4)
!c
      INTEGER          I
      DOUBLE COMPLEX   AV, BV, CV, DV
      DOUBLE COMPLEX   CZERO, CONE
!C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
!c
      AV = U
      BV = AV*AV
      IF (CDABS(BV).GE.0.37D0) THEN
       BV = CONE/(CONE-BV)
       CV = BV*BV
       DX0(1) = (BV-7.0D0*X0)/AV
       DX0(2) = 2.0D0*(CV-4.0D0*DX0(1)/AV)
       DX0(3) = 8.0D0*AV*CV*BV + (2.0D0*CV-9.0D0*DX0(2))/AV
       DX0(4) = 24.0D0*CV*BV*(CONE+2.0D0*AV*AV*BV) - 10.0D0*DX0(3)/AV
      ELSE
       DX0(3) = CZERO
       DX0(2) = (2.0D0,0.0D0)/9.0D0
       DX0(1) = DX0(2)*AV
       DX0(4) = (24.0D0,0.0D0)/11.0D0
       CV = AV
       BV = AV*AV
       DO 50 I = 2, 48, 2
        DV = DBLE(2+I)/DBLE(9+I)*CV
        DX0(4) = DX0(4) + DBLE((4+I)*(3+I)*(2+I)*(1+I))/DBLE(11+I)*CV*AV
        DX0(3) = DX0(3) + DBLE((1+I)*I)*DV
        DX0(2) = DX0(2) + DBLE(1+I)*DV*AV
        DX0(1) = DX0(1) + DV*BV
        CV = CV*BV
 50    CONTINUE
      END IF
!C
      RETURN
END SUBROUTINE SIM0UDX0
!C
!C ----------------------------------------------------------------------
!C file: sim1onei.f
!C date: 2010-09-22
!C who:  S.Kaprzyk
!C what: CASE a): [w1<w2<w3<ONE]
!C       U(1)=W(1); U(2)=W(2); U(3)=W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE SIM1ONEI(U, AS, LEPS, SIM1, SIM1I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM1, SIM1I
!C
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(3), QX0(3), DX0(4)
      DOUBLE COMPLEX   R0(3), QR0(3)
      DOUBLE COMPLEX   S0(3), QS0(3)
      DOUBLE COMPLEX   AV, BV, CV, DU, G(5)
      INTEGER          I, K !, N
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(as)

      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
!C ------------------------------------------------------------------
!C      R0(U) = ARTH(U)/U
!C      R0(U) = 1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   The DX0(I) contains  derivatives of the function X0(u)      -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
      DO 200 I = 1, 3, 2
       DU = U(I) - U(2)
       IF (.NOT.LEPS(I)) THEN
        QX0(I) = (X0(I)-X0(2))/DU
       ELSE
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*DU) *DU)*DU
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
       QS0(I) = (-S0(2)-R0(2)+(CONE-U(I))*QR0(I))/(CONE+U(I))
 200  CONTINUE ! end DO ... I = 1, 3, 2
      IF (LEPS(2)) THEN
!*-ASIS
       QX0(2) = DX0(2)/2.0D0 + DX0(3)*(U(1)+U(3)-2.0D0*U(2))/6.0D0+ &
     &          DX0(4)*((U(3)-U(2))*(U(3)-U(2))+(U(1)-U(2))*(U(1)-U(2))+ &
     &         (U(3)-U(2))*(U(1)-U(2)))/24.0D0
      ELSE
       QX0(2) = (QX0(3)-QX0(1))/(U(3)-U(1))
      END IF
      AV = U(1)
      G(1) = U(1) + U(3)
      DO 300 K = 2, 5
       AV = AV*U(1)
       G(K) = G(K-1)*U(3) + AV
 300  CONTINUE
!*-ASIS
      QR0(2) = CONE/3.0D0 + (G(2)+U(2)*(U(2)+G(1)))/5.0D0 + &
     &   (G(4)+U(2)*(G(3)+U(2)*(G(2)+U(2)*(G(1)+U(2)))))*X0(2)+ &
     &    G(5)*QX0(3) + AV*U(1)*QX0(2)
      QS0(2) = (-QS0(3)-QR0(3)+(CONE-U(1))*QR0(2))/(CONE+U(1))
!c
      AV = (CONE+U(1))*(CONE+U(2))*(CONE+U(3))
      BV = U(3) + U(1) - (2.0D0,0.0D0)
      CV = (CONE-U(1))*(CONE-U(1))
      SIM1 = AV*(R0(2)+QR0(3)*BV+QR0(2)*CV)
      SIM1I = AV*(S0(2)+QS0(3)*BV+QS0(2)*CV)
      RETURN
END SUBROUTINE SIM1ONEI
!C
!C ----------------------------------------------------------------------
!C file: sim1twoi.f
!C date: 2010-09-22
!C who:  S.Kaprzyk
!C what: CASE a):  [w1<w2<w3<ONE]
!C       U(1)=W(1); U(2)=W(2); U(3)=W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE SIM1TWOI(U, AS, LEPS, SIM, SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX R4(3), X0(3), DX0(4)
      DOUBLE COMPLEX QR4(3), QX0(3)
      DOUBLE COMPLEX S4(3), QS4(3)
      DOUBLE COMPLEX AV, BV, CV, UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(as)
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
!C -----------------------------------------------------------------
!C -                                                               -
!C -        R4(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
!C -    R4(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C -                                                               -
!C -----------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R4(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
       S4(I) = (CONE-U(I))/(CONE+U(I))*R4(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   The DX0(I) Contains  Derivatives Of The Function X0(W)      -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!C
      DO 200 I = 1, 3, 2
       UV = U(I) - U(2)
       IF (LEPS(I)) THEN
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UV) *UV)*UV
       ELSE
        QX0(I) = (X0(I)-X0(2))/UV
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR4(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 + G(3) &
     &  /5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
       QS4(I) = (-S4(2)-R4(2)+(CONE-U(I))*QR4(I))/(CONE+U(I))
 200  CONTINUE  ! DO 750 I = 1, 3, 2
      IF (LEPS(2)) THEN
!*-ASIS
        QX0(2) = DX0(2)/2.0D0 + DX0(3)*(U(1)+U(3)-2.0D0*U(2))/6.0D0 + &
     &  DX0(4)*((U(3)-U(2))*(U(3)-U(2))+(U(1)-U(2))*(U(1)-U(2))+ &
     &  (U(3)-U(2))*(U(1)-U(2)))/24.0D0
      ELSE
       QX0(2) = (QX0(3)-QX0(1))/(U(3)-U(1))
      END IF
      AV = U(1)
      G(1) = U(1) + U(3)
      DO 300 K = 2, 5
       AV = AV*U(1)
       G(K) = G(K-1)*U(3) + AV
 300  CONTINUE
!*-ASIS
      QR4(2) = CONE/3.0D0- (U(2)+G(1))/5.0D0 + &
     &  (G(2)+U(2)*(U(2)+G(1)))/5.0D0 + &
     &  (G(4)+(U(2)-CONE)*(G(3) + G(2)*U(2) + G(1)*U(2)*U(2) + &
     &   U(2)*U(2)*U(2)))*X0(2) + (G(5)-G(4))*QX0(3) + &
     &  (U(1)-CONE)*AV*QX0(2)
      QS4(2) = (-QS4(3)-QR4(3)+(CONE-U(1))*QR4(2))/(CONE+U(1))
!c
      AV = (CONE+U(1))*(CONE+U(2))*(CONE+U(3))
      BV = U(3) + U(1) - (2.0D0,0.0D0)
      CV = (CONE-U(1))*(CONE-U(1))
      SIM = AV*(R4(2)+QR4(3)*BV+QR4(2)*CV)
      SIMI = AV*(S4(2)+QS4(3)*BV+QS4(2)*CV)
      RETURN
END SUBROUTINE SIM1TWOI
!C
!C ----------------------------------------------------------------------
!C file: sim2onei.f
!C date: 2010-08-19
!C who:  S.Kaprzyk
!C what: CASE b) : [w1<w2<ONE<w3]
!C       U(1)=W(1); U(2)=W(2); U(3)=CONE/W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
      SUBROUTINE  SIM2ONEI(U,AS,LEPS, SIM2, SIM2I)

      IMPLICIT         NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM2, SIM2I
!c
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(3), QX0(3), DX0(4)
      DOUBLE COMPLEX   R0(3), QR0(3)
      DOUBLE COMPLEX   S0(3), QS0(3)
      DOUBLE COMPLEX   AV, BV, CV, DV, G(5)
      DOUBLE COMPLEX   UD
      INTEGER          I, K !, N
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
!C ------------------------------------------------------------------
!C      R0(U) = ARTH(U)/U
!C      R0(U) = 1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   THE DX0(I) CONTAINS  DERIVATIVES OF THE FUNCTION X0(W)         -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!c      DO 400 I = 1, 3, 2
      I = 1
      UD = U(I) - U(2)
      IF (.NOT.LEPS(I)) THEN
       QX0(I) = (X0(I)-X0(2))/UD
      ELSE
       QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UD) *UD)*UD
      END IF
      AV = U(I)
      G(1) = U(I) + U(2)
      DO 200 K = 2, 5
       AV = AV*U(I)
       G(K) = G(K-1)*U(2) + AV
 200  CONTINUE
      QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
      QS0(I) = (-S0(2)-R0(2)+(CONE-U(I))*QR0(I))/(CONE+U(I))
!c 400  CONTINUE ! end DO ... I = 1, 3, 2
      AV = (CONE+U(1))*(CONE+U(2))*(U(3)+CONE) /((U(3)*U(1)-CONE)*(U(3)*U(2)-CONE))
      BV = CTWO - (U(1)+U(2)) - U(3)*(CONE-U(1)*U(2))
      CV = (U(3)*U(2)-CONE)*(CONE-U(1))**2
      DV = (U(3)-CONE)**2
      SIM2 = AV*(BV*R0(2)+CV*QR0(1)+DV*(DCMPLX(ZERO,AS(3))+U(3)*R0(3)))
      DV = DV*(U(3)-CONE)/(U(3)+CONE)
      SIM2I = AV*(BV*S0(2)+CV*QS0(1)+DV*(DCMPLX(ZERO,AS(3))+U(3)*R0(3)))
!c
      RETURN
END SUBROUTINE SIM2ONEI
!C
!C ----------------------------------------------------------------------
!C file: sim2twoi.f
!C date: 2010-08-19
!C who:  S.Kaprzyk
!C what: CASE b) : [w1<w2<ONE<w3]
!C       U(1)=W(1); U(2)=W(2); U(3)=CONE/W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE SIM2TWOI(U,AS,LEPS, SIM, SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX R4(3), QR4(3)
      DOUBLE COMPLEX S4(3), QS4(3)
      DOUBLE COMPLEX X0(3), QX0(3), DX0(4), G(5)
      DOUBLE COMPLEX AV, BV, CV, DV, UV
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
!C   -----------------------------------------------------------------
!C   -                                                               -
!C   -        R4(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
!C   -    R4(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C   -                                                               -
!C   -----------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R4(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
       S4(I) = (CONE-U(I))/(CONE+U(I))*R4(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   The DX0(I) Contains  Derivatives Of The Function X0(W)      -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!c      DO 200 I = 1, 3, 2
      I = 1
      UV = U(I) - U(2)
      IF (LEPS(I)) THEN
       QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UV) *UV)*UV
      ELSE
       QX0(I) = (X0(I)-X0(2))/UV
      END IF
      AV = U(I)
      G(1) = U(I) + U(2)
      DO 200 K = 2, 5
       AV = AV*U(I)
       G(K) = G(K-1)*U(2) + AV
 200  CONTINUE
!*-ASIS
       QR4(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 + &
     &   G(3) /5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
      QS4(I) = (-S4(2)-R4(2)+(CONE-U(I))*QR4(I))/(CONE+U(I))
!c  200  CONTINUE  ! DO ... I = 1, 3, 2
      AV = (CONE+U(1))*(CONE+U(2))*(U(3)+CONE) &
     & /((U(3)*U(1)-CONE)*(U(3)*U(2)-CONE))
      BV = CTWO - (U(1)+U(2)) - U(3)*(CONE-U(1)*U(2))
      CV = (U(3)*U(2)-CONE)*(CONE-U(1))**2
      DV = (U(3)-CONE)**2
      SIM = AV*(BV*R4(2)+CV*QR4(1) &
     & +DV*((CONE-U(3))*DCMPLX(ZERO,AS(3))+CONE+U(3)-U(3)*U(3)*R4(3)))
      DV = DV*(U(3)-CONE)/(U(3)+CONE)
      SIMI = AV*(BV*S4(2)+CV*QS4(1) &
     & +DV*((CONE-U(3))*DCMPLX(ZERO,AS(3))+CONE+U(3)-U(3)*U(3)*R4(3)))
!C
      RETURN
END SUBROUTINE SIM2TWOI
!C
!C ----------------------------------------------------------------------
!C file: sim3onei.f
!C date: 2010-08-19
!C who:  S.Kaprzyk
!C what: CASE c) : [w1<ONE<w2<w3]
!C       U(1)=W(1); U(2)=CONE/W(2); U(3)=CONE/W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
      SUBROUTINE SIM3ONEI(U,AS,LEPS,SIM3,SIM3I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM3, SIM3I
!C
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(3), QX0(3), DX0(4)
      DOUBLE COMPLEX   R0(3), QR0(3)
      DOUBLE COMPLEX   S0(3), QS0(3)
      DOUBLE COMPLEX   AV, BV, CV, DV, UD, G(5)
      INTEGER          I,  K !, N
      DOUBLE COMPLEX   CZERO, CONE, CTWO, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
!C ------------------------------------------------------------------
!C      R0(U) = ARTH(U)/U
!C      R0(U) = 1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
!c       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   THE DX0(I) CONTAINS  DERIVATIVES OF THE FUNCTION X0(U)         -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!c      DO 400 I = 1, 3, 2
      I = 3
      UD = U(I) - U(2)
      IF (.NOT.LEPS(I)) THEN
       QX0(I) = (X0(I)-X0(2))/UD
      ELSE
       QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UD) *UD)*UD
      END IF
      AV = U(I)
      G(1) = U(I) + U(2)
      DO 200 K = 2, 5
       AV = AV*U(I)
       G(K) = G(K-1)*U(2) + AV
 200  CONTINUE
      QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
!c 400  CONTINUE ! end DO ... I = 1, 3, 2
      QR0(3) = R0(2) + U(3)*QR0(3)
      IF (.NOT.LEPS(3)) THEN
       QR0(3) = QR0(3) + DCMPLX(ZERO,(AS(3)-AS(2)))/(U(3)-U(2))
      END IF
      R0(2) = DCMPLX(ZERO,AS(2)) + U(2)*R0(2)
      AV = (ONE+U(1))*(U(2)+CONE)*(U(3)+CONE) /((U(1)*U(2)-CONE)*(U(1)*U(3)-CONE))
      BV = (CONE-U(1))**2
      CV = CTWO - (U(2)+U(3)) - U(1)*(CONE-U(2)*U(3))
      DV = (U(1)*U(2)-CONE)*(U(3)-CONE)**2
      SIM3 = AV*(BV*R0(1)+CV*R0(2)+DV*QR0(3))
!c
      S0(1) = (CONE-U(1))/(CONE+U(1))*R0(1)
      S0(2) = (U(2)-CONE)/(U(2)+CONE)*R0(2)
      QS0(3) = (-S0(2)+R0(2)+(U(3)-CONE)*QR0(3))/(U(3)+CONE)
      SIM3I = AV*(BV*S0(1)+CV*S0(2)+DV*QS0(3))
!c
      RETURN
END SUBROUTINE SIM3ONEI
!C
!C ----------------------------------------------------------------------
!C file: sim3twoi.f
!C date: 2010-08-19
!C who:  S.Kaprzyk
!C what: CASE c) : [w1<ONE<w2<w3]
!C       U(1)=W(1); U(2)=CONE/W(2); U(3)=CONE/W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE SIM3TWOI(U,AS,LEPS,SIM,SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!C
      DOUBLE COMPLEX X0(3), QX0(3), DX0(4)
      DOUBLE COMPLEX R4(3), QR4(3)
      DOUBLE COMPLEX S4(3), QS4(3)
      DOUBLE COMPLEX AV, BV, CV, DV, UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CIMAG, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
!C   -----------------------------------------------------------------
!C   -                                                               -
!C   -        R4(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
!C   -    R4(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C   -                                                               -
!C   -----------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R4(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
 100  CONTINUE  ! end DO 350 I = 1, 3
!C -----------------------------------------------------------------
!C -   The DX0(I) Contains  Derivatives Of The Function X0(W)      -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!c      DO 750 I = 1, 3, 2
      I = 3
      UV = U(I) - U(2)
      IF (.NOT.LEPS(I)) THEN
       QX0(I) = (X0(I)-X0(2))/UV
      ELSE
       QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UV) *UV)*UV
      END IF
      AV = U(I)
      G(1) = U(I) + U(2)
      DO 200 K = 2, 5
       AV = AV*U(I)
       G(K) = G(K-1)*U(2) + AV
 200  CONTINUE
!*-ASIS
      QR4(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 + G(3)/5.0D0 + &
     &  (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
!C 750   CONTINUE ! DO .. I = 1, 3, 2
      QR4(3) = -DCMPLX(ZERO,AS(2)) + CONE - (U(3)+U(2))*R4(2) - U(3)*U(3)*QR4(3)
      IF (.NOT.LEPS(3)) THEN
       QR4(3) = QR4(3) + (CONE-U(3))*DCMPLX(ZERO,(AS(3)-AS(2))) /(U(3)-U(2))
      END IF
      R4(2) = (CONE-U(2))*DCMPLX(ZERO,AS(2)) + CONE + U(2) - U(2)*U(2)*R4(2)
!c
      AV = (ONE+U(1))*(U(2)+CONE)*(U(3)+CONE)/((U(1)*U(2)-CONE)*(U(1)*U(3)-CONE))
      BV = (CONE-U(1))**2
      CV = CTWO - (U(2)+U(3)) - U(1)*(CONE-U(2)*U(3))
      DV = (U(1)*U(2)-CONE)*(U(3)-CONE)**2
      SIM = AV*(BV*R4(1)+CV*R4(2)+DV*QR4(3))
!c
      S4(1) = (CONE-U(1))/(CONE+U(1))*R4(1)
      S4(2) = (U(2)-CONE)/(U(2)+CONE)*R4(2)
      QS4(3) = (-S4(2)+R4(2)+(U(3)-CONE)*QR4(3))/(U(3)+CONE)
      SIMI = AV*(BV*S4(1)+CV*S4(2)+DV*QS4(3))
!c
      RETURN
      END SUBROUTINE SIM3TWOI
!C
!C ----------------------------------------------------------------------
!C file: sim4onei.f
!C date: 2010-09-22
!C who:  S.Kaprzyk
!C what: CASE d): [ ONE<w1<w2<w3 ]
!C       U(1)=CONE/W(1); U(2)=CONE/W(2); U(3)=CONE/W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE  SIM4ONEI(U,AS,LEPS,SIM4,SIM4I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM4, SIM4I
!c
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(3), QX0(3), DX0(4)
      DOUBLE COMPLEX   R0(3), QR0(3)
      DOUBLE COMPLEX   S0(3), QS0(3)
      DOUBLE COMPLEX   QAS(3)
      DOUBLE COMPLEX   AV, BV, CV, UD,  G(5)
      INTEGER          I, K !, N
      DOUBLE COMPLEX   CZERO, CONE, CTWO, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
!C ------------------------------------------------------------------
!C       R0(U)=ARTH(U)/U
!C       R0(U)=1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
!c       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   THE DX0(I) CONTAINS  DERIVATIVES OF THE FUNCTION X0(U)
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
      DO 200 I = 1, 3, 2
       UD = U(I) - U(2)
       IF (.NOT.LEPS(I)) THEN
        QX0(I) = (X0(I)-X0(2))/UD
        QAS(I) = (AS(I)-AS(2))/UD
       ELSE
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UD) *UD)*UD
        QAS(I) = CZERO
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
 200  CONTINUE ! end DO ... I = 1, 3, 2
      IF (LEPS(2)) THEN
       QAS(2) = CZERO
!*-ASIS
       QX0(2) = DX0(2)/2.0D0 + DX0(3)*(U(1)+U(3)-2.0D0*U(2))/6.0D0 + &
     &         DX0(4)*((U(3)-U(2))*(U(3)-U(2))+(U(1)-U(2))*(U(1)-U(2)) + &
     &        (U(3)-U(2))*(U(1)-U(2)))/24.0D0
      ELSE
       QAS(2) = (QAS(3)-QAS(1))/(U(3)-U(1))
       QX0(2) = (QX0(3)-QX0(1))/(U(3)-U(1))
      END IF
      AV = U(1)
      G(1) = U(1) + U(3)
      DO 300 K = 2, 5
       AV = AV*U(1)
       G(K) = G(K-1)*U(3) + AV
 300  CONTINUE
!*-ASIS
       QR0(2) = CONE/3.0D0 + (G(2)+U(2)*(U(2)+G(1)))/5.0D0 + &
     &   (G(4)+U(2)*(G(3)+U(2)*(G(2)+U(2)*(G(1)+U(2)))))*X0(2)+ &
     &    G(5)*QX0(3) + AV*U(1)*QX0(2)
!c
      QR0(2) = CIMAG*QAS(2) + QR0(3) + U(1)*QR0(2)
      QR0(3) = CIMAG*QAS(3) + R0(2) + U(3)*QR0(3)
      R0(2) = CIMAG*AS(2) + U(2)*R0(2)
!c
      S0(2) = (U(2)-CONE)/(U(2)+CONE)*R0(2)
      QS0(3) = (-S0(2)+R0(2)+(U(3)-CONE)*QR0(3))/(U(3)+CONE)
      QS0(2) = (-QS0(3)+QR0(3)+(U(1)-CONE)*QR0(2))/(U(1)+CONE)
!c
      AV = (CONE+U(1))*(CONE+U(2))*(CONE+U(3))
      BV = U(3) + U(1) - CTWO
      CV = (CONE-U(1))**2
      SIM4 = AV*(R0(2)+BV*QR0(3)+CV*QR0(2))
      SIM4I = AV*(S0(2)+BV*QS0(3)+CV*QS0(2))
!c
      RETURN
END SUBROUTINE SIM4ONEI
!C
!C ----------------------------------------------------------------------
!C file: sim4twoi.f
!C date: 2010-09-22
!C who:  S.Kaprzyk
!C what: CASE d): [ ONE<w1<w2<w3 ]
!C       U(1)=CONE/W(1); U(2)=CONE/W(2); U(3)=CONE/W(3)
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE SIM4TWOI(U,AS,LEPS,SIM,SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(3)
      DOUBLE PRECISION AS(3)
      LOGICAL          LEPS(3)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX X0(3), QX0(3), DX0(4)
      DOUBLE COMPLEX R4(3), QR4(3)
      DOUBLE COMPLEX S4(3), QS4(3)
!c      DOUBLE COMPLEX RU(3), TU(3)
      DOUBLE COMPLEX QAS(3), AV, BV, CV, UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      CTWO = (2.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
!C   -----------------------------------------------------------------
!C   -                                                               -
!C   -        R4(U)=1+(U-1)*(ARTH(U)/U-1)/U                          -
!C   -    R4(U)=1-U/3+U**2/3-U**3/5+U**4/5+(U-1)*U**5*X0(U)          -
!C   -                                                               -
!C   -----------------------------------------------------------------
      DO 100 I = 1, 3
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R4(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &  AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
 100  CONTINUE  ! end DO ... I = 1, 3
!C -----------------------------------------------------------------
!C -   The DX0(I) Contains  Derivatives Of The Function X0(U)      -
!C -----------------------------------------------------------------
      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!C
      DO 200 I = 1, 3, 2
       UV = U(I) - U(2)
       IF (LEPS(I)) THEN
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UV) *UV)*UV
        QAS(I) = CZERO
       ELSE
        QX0(I) = (X0(I)-X0(2))/UV
        QAS(I) = (AS(I)-AS(2))/UV
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR4(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 + G(3) &
     &  /5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
 200  CONTINUE  ! end DO ... I = 1, 3, 2
!C
      IF (LEPS(2)) THEN
       QAS(2) = CZERO
!*-ASIS
        QX0(2) = DX0(2)/2.0D0 + DX0(3)*(U(1)+U(3)-2.0D0*U(2))/6.0D0 + &
     &  DX0(4)*((U(3)-U(2))*(U(3)-U(2))+(U(1)-U(2))*(U(1)-U(2)) + &
     &  (U(3)-U(2))*(U(1)-U(2)))/24.0D0
      ELSE
       QAS(2) = (QAS(3)-QAS(1))/(U(3)-U(1))
       QX0(2) = (QX0(3)-QX0(1))/(U(3)-U(1))
      END IF
      AV = U(1)
      G(1) = U(1) + U(3)
      DO 300 K = 2, 5
       AV = AV*U(1)
       G(K) = G(K-1)*U(3) + AV
 300  CONTINUE
!*-ASIS
      QR4(2) = CONE/3.0D0- (U(2)+G(1))/5.0D0 + &
     &  (G(2)+U(2)*(U(2)+G(1)))/5.0D0 + &
     &  (G(4)+(U(2)-CONE)*(G(3) + G(2)*U(2) + G(1)*U(2)*U(2) + &
     &   U(2)*U(2)*U(2)))*X0(2) + (G(5)-G(4))*QX0(3) + &
     &   (U(1)-CONE)*AV*QX0(2)
!C
      QR4(2) = CIMAG*(-QAS(3)+(CONE-U(1))*QAS(2)) &
     & - (R4(2)+(U(3)+U(1))*QR4(3)+U(1)*U(1)*QR4(2))
      QR4(3) = CIMAG*(-AS(2)+(CONE-U(3))*QAS(3)) + CONE - (U(3)+U(2)) &
     & *R4(2) - U(3)*U(3)*QR4(3)
      R4(2) = CIMAG*(CONE-U(2))*AS(2) + CONE + U(2) - U(2)*U(2)*R4(2)
!c
      S4(2) = (U(2)-CONE)/(U(2)+CONE)*R4(2)
      QS4(3) = (-S4(2)+R4(2)+(U(3)-CONE)*QR4(3))/(U(3)+CONE)
      QS4(2) = (-QS4(3)+QR4(3)+(U(1)-CONE)*QR4(2))/(U(1)+CONE)
!c
      AV = (CONE+U(1))*(CONE+U(2))*(CONE+U(3))
      BV = U(3) + U(1) - CTWO
      CV = (CONE-U(1))**2
      SIM = AV*(R4(2)+BV*QR4(3)+CV*QR4(2))
      SIMI = AV*(S4(2)+BV*QS4(3)+CV*QS4(2))
!c
      RETURN
END SUBROUTINE SIM4TWOI
!C
!C ----------------------------------------------------------------------
!C file: sim0leps.f
!C date: 2011-03-12
!C who:  S.Kaprzyk
!C what: Return a proper case a), b), c) or d): see, Fig.3.1
!C INPUT:
!C VERM(1..4) - four complex numbers
!C N - integer number for W(I)=(VERM(I)-VERM(N))/(VERM(I)+VERM(N))
!C EPS  - critical distance |w1-w2| etc.
!C OUTPUT:
!C LONE(1..4) [w1<w2<w3<ONE]; [w1<w2<ONE<w3];
!C            [w1<ONE<w2<w3]; [ONE<w1<w2<w3]
!C LEPS(1..3) [|w1-w2|<e,|w3-w2|>e];
!C            [|w1-w2|<e,|w3-w2|<e];
!C            [|w1-w2|>e,|w3-w2|<e]
SUBROUTINE  SIM0LEPS(VERM, N, W, AS, LONE, EPS, LEPS, iuerr)

      IMPLICIT       NONE
      DOUBLE COMPLEX   VERM(4)
      INTEGER          N
      DOUBLE COMPLEX   W(3)
      DOUBLE PRECISION AS(3)
      DOUBLE PRECISION EPS
      LOGICAL          LONE(4), LEPS(3)
      INTEGER          iuerr
!C
!c      DOUBLE PRECISION DLAMCH
!c      EXTERNAl         DLAMCH
!C
      LOGICAL          LWONE(3), LV
      DOUBLE COMPLEX   AV
      DOUBLE PRECISION AL(3), AA, AIW, OZERO
      INTEGER          I, II, K
      DOUBLE PRECISION ZERO, ONE, PI
      DATA             ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/
      DATA             OZERO/ 1.0D-13/
!C
!c      OZERO = DLAMCH('e')*10
      IF ((N.GT.4).OR.(N.LT.1)) THEN
       WRITE(iuerr,9003) N
 9003  FORMAT(' ***sim0leps: N =',I6,' must be 1,2,3 or 4' )
       STOP ' ***sim0leps: '
      ENDIF
!c      DO 100 I = 1, 4
!c       IF (I.EQ.N) GO TO 100
!c       IF (DIMAG(VERM(N))*DIMAG(VERM(I)).LT.ZERO) THEN
!c        WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,4)
!c 9010   FORMAT (' ***sim0leps: not signed ImgVERM()=',4(D13.6,1X))
!c        STOP ' ***sim0leps: '
!c       END IF
!c 100  CONTINUE
      II = 0
      DO 200 I = 1, 4
       IF (I.EQ.N) GO TO 200
       II = II + 1
       IF (CDABS(VERM(N)-VERM(I)).LT.CDABS(VERM(N)+VERM(I))) THEN
        LWONE(II) = .TRUE.
        W(II) = (VERM(N)-VERM(I))/(VERM(N)+VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = 0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       ELSE
        LWONE(II) = .FALSE.
        W(II) = (VERM(N)+VERM(I))/(VERM(N)-VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = -0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       END IF
       AL(II) = CDABS(W(II))
 200  CONTINUE
      DO 300 I = 1, 3
       DO 250 K = I, 3
        IF (LWONE(I).AND.LWONE(K).AND.(AL(I).LT.AL(K))) GO TO 250
        IF (.NOT.LWONE(I).AND..NOT.LWONE(K).AND.(AL(K).LT.AL(I))) GO TO 250
        IF (.NOT.LWONE(I).AND.LWONE(K).AND.(ONE.LT.AL(I)*AL(K))) GO TO 250
        IF (LWONE(I).AND..NOT.LWONE(K).AND.(AL(I)*AL(K).LT.ONE)) GO TO 250
        AA = AL(K)
        AL(K) = AL(I)
        AL(I) = AA
        AA = AS(K)
        AS(K) = AS(I)
        AS(I) = AA
        AV = W(K)
        W(K) = W(I)
        W(I) = AV
        LV = LWONE(K)
        LWONE(K) = LWONE(I)
        LWONE(I) = LV
 250   CONTINUE
 300  CONTINUE
!C
      LONE(1) = .FALSE.
      LONE(2) = .FALSE.
      LONE(3) = .FALSE.
      LONE(4) = .FALSE.
      IF (LWONE(1).AND.LWONE(2).AND.LWONE(3)) THEN
       LONE(1) = .TRUE.
      END IF
      IF (LWONE(1).AND.LWONE(2).AND.(.NOT.LWONE(3))) THEN
       LONE(2) = .TRUE.
      END IF
      IF (LWONE(1).AND.(.NOT.LWONE(2)).AND.(.NOT.LWONE(3))) THEN
       LONE(3) = .TRUE.
      END IF
      IF ((.NOT.LWONE(1)).AND.(.NOT.LWONE(2)).AND.(.NOT.LWONE(3))) THEN
       LONE(4) = .TRUE.
      END IF
!C here are 3-cases, how close are w()
      LEPS(1) = .FALSE.
      LEPS(2) = .FALSE.
      LEPS(3) = .FALSE.
      IF ( LONE(1).OR.LONE(4) ) THEN
       IF (ABS(W(1)-W(2)).LE.EPS) THEN
        LEPS(1) = .TRUE.
       END IF
       IF (ABS(W(3)-W(2)).LE.EPS) THEN
        LEPS(3) = .TRUE.
       END IF
       IF ((ABS(W(1)-W(2)).LE.EPS).AND.(ABS(W(3)-W(2)).LE.EPS)) THEN
        LEPS(2) = .TRUE.
       END IF
      END IF
!c
      IF ( LONE(2)) THEN
       IF (ABS(W(1)-W(2)).LE.EPS) THEN
        LEPS(1) = .TRUE.
       END IF
      END IF
!c
      IF ( LONE(3)) THEN
       IF (ABS(W(2)-W(3)).LE.EPS) THEN
         LEPS(3) = .TRUE.
       END IF
      END IF
!C
      IF ((.NOT.LONE(1)).AND.(.NOT.LONE(2)).AND.(.NOT.LONE(3)).AND. (.NOT.LONE(4))) THEN
       WRITE (iuerr,9020) (LONE(I),I=1,4)
 9020  FORMAT ('***sim0leps: LONE()=',4(L1,1X),'not in order')
       STOP '***sim0leps: '
      END IF
      RETURN
END SUBROUTINE SIM0LEPS
!C
!C file: s2d0onei.f
!C date: 2011-04-01
!C who:  S.Kaprzyk
!C what: Complex linear form integral over standard 2-d simplex
!C ----------------------------------------------------
!C -         *                                        -
!C -        *                         1               -
!C - SIM0= 2* dt1dt2 ----------------------------------
!C -        *     verm(3)+(verm(1)-verm(3))*t1+..     -
!C -       *                                          -
!C -    0<t1<1                                        -
!C ----------------------------------------------------
SUBROUTINE S2D0ONEI(SIM0, SIM0I, VERM)

      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX SIM0, SIM0I
      DOUBLE COMPLEX VERM (3)
!C
      DOUBLE COMPLEX   VERM1D(2)
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(3),LEPS(1)
!c
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2), AL(3)
      INTEGER          I, I1, I2, K, N
      DOUBLE PRECISION ZERO, EPS, SMALL
      DATA             EPS/1.0D-6/, SMALL/1.0D-5/
      DATA             ZERO/0.0D0/
!C
      DO 100 I = 1, 2
       IF (DIMAG(VERM(3))*DIMAG(VERM(I)).LT.ZERO) THEN
        WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,3)
 9010   FORMAT (' ***s2d0onei: not signed ImgVERM()=',3(D13.6,1X))
!c        STOP ' ***s2d0onei: '
       END IF
 100  CONTINUE
!C
      DO 200 I = 1, 3
       AL(I) = CDABS(VERM(I))
 200  CONTINUE
      DO 300 N = 1, 3
       IF (AL(N).LT.SMALL) THEN
        I1 = MOD(N,3) + 1
        I2 = MOD(N+1,3) + 1
        VERM1D(1) = VERM(I1)
        VERM1D(2) = VERM(I2)
        CALL S1D0ONEI(SIM,SIMI,VERM1D)
        SIM0  = 2*SIM
        SIM0I = SIMI - 0.5D0
        RETURN
       END IF
 300  CONTINUE
!C
      N = 3
      DO 400 I = 1, 2
       IF (AL(I).GT.AL(N)) N = I
 400  CONTINUE
!C Here are 3 cases, how w() are placed on complex-plane
!C LONE(1..3) [w1<w2<ONE];[w1<ONE<w2];[ONE<w1<w2]
!C LEPS(1) [|w1-w2|<e
!*-ASIS
      CALL S2D0LEPS(VERM, N, U, AS, LONE, EPS, LEPS, iuerr)
!c U(1) = W(1); U(2) = W(2);
      IF (LONE(1)) THEN
       CALL S2D1ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!c U(1) = W(1); U(2) = 1/W(2);
      IF (LONE(2)) THEN
       CALL S2D2ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!c U(1) = 1/W(1); U(2) = 1/W(2);
      IF (LONE(3)) THEN
       CALL S2D3ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!C
      SIM0 = SIM/VERM(N)
      SIM0I = log(VERM(N))+0.5D0*SIMI - 1.5D0
      RETURN
END SUBROUTINE S2D0ONEI
!C
!C ----------------------------------------------------------------------
!C file: s2d0twoi.f
!C date: 2011-04-01
!C who:  S.Kaprzyk
!C what: Complex linear form integrals over standard 2-d simplex
!C ------------------(i=1,2)----------------------------------------
!C            *                                                    -
!C           *            t_i                                      -
!C  VERL(i)=2*dt1dt2 --------------------------------------------- -
!C           *          VERM(3)+[VERM(1)-VERM(3)]*t1+...           -
!C          *                                                      -
!C     0<t1+t2<1                                                   -
!C -------------------(i=3)-----------------------------------------
!C            *                                                    -
!C           *           1- t1 -t2
!C  VERL(3)=2*dt1dt2 -----------------------------------------------
!C           *          VERM(3)+[VERM(1)-VERM(3)]*t1+...           -
!C          *                                                      -
!C     0<t1+t2<1                                                   -
!C -----------------------------------------------------------------
SUBROUTINE   S2D0TWOI(VERL, VERLI, VERM)

      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX VERL(3), VERLI(3), VERM(3)
!C
      DOUBLE COMPLEX  SIM, SIMI
      DOUBLE COMPLEX  VERL1D(2), VERL1DI(2), VERM1D(2)
      LOGICAL         LONE(3),LEPS(2)
!c
      INTEGER          I, I1, I2, K, N
      DOUBLE PRECISION AS(2), AL(3)
      DOUBLE COMPLEX   U(2)
      DOUBLE COMPLEX   CZERO
      DOUBLE PRECISION EPS, SMALL
      DOUBLE PRECISION ZERO
      DATA             EPS/1.0D-6/, SMALL/1.0D-5/
      DATA             ZERO/0.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      DO 100 I = 1, 2
       IF (DIMAG(VERM(3))*DIMAG(VERM(I)).LT.ZERO) THEN
        WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,3)
 9010   FORMAT (' ***s2d0twoi: not signed ImgVERM()=',3(D13.6,1X))
        STOP ' ***s2d0twoi: '
       END IF
 100  CONTINUE
      DO 200 I = 1, 3
       AL(I) = CDABS(VERM(I))
 200  CONTINUE
!C
      DO 300 N = 1, 3
       VERL(N) = CZERO
       VERLI(N) = CZERO
       I1 = MOD(N,3) + 1
       I2 = MOD(N+1,3) + 1
       IF (AL(I1).LT.SMALL) THEN
        VERM1D(1) = VERM(I2)
        VERM1D(2) = VERM(N)
        CALL S1D0TWOI(VERL1D,VERL1DI,VERM1D)
        VERL(N) = VERL1D(2)
        VERLI(N) = VERL1DI(2)*2/3
        GO TO 300
       END IF
       IF (AL(I2).LT.SMALL) THEN
        VERM1D(1) = VERM(I1)
        VERM1D(2) = VERM(N)
        CALL S1D0TWOI(VERL1D,VERL1DI,VERM1D)
        VERL(N) = VERL1D(2)
        VERLI(N) = VERL1DI(2)*2/3
        GO TO 300
       END IF
       IF (AL(N).LT.SMALL) THEN
        VERM1D(1) = VERM(I1)
        VERM1D(2) = VERM(I2)
        CALL S1D0ONEI(SIM,SIMI,VERM1D)
        VERL(N) = SIM
        VERLI(N) = SIMI/3 - 0.5D0
        GO TO 300
       END IF
!C Here are 3 cases, how w() are placed on complex-plane
!C LONE(1..3) [w1<w2<ONE];[w1<ONE<w2];[ONE<w1<w2]
!C LEPS(1) [|w1-w2|<e
!*-ASIS
       CALL S2D0LEPS(VERM, N, U, AS, LONE, EPS, LEPS, iuerr)
!c U(1) = W(1); U(2) = W(2)
       IF (LONE(1)) THEN
        CALL S2D1TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c U(1) = W(1); U(2) = CONE/W(2)
       IF (LONE(2)) THEN
        CALL S2D2TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c U(1) = CONE/W(1); U(2) = CONE/W(2)
       IF (LONE(3)) THEN
        CALL S2D3TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c
       VERL(N) = SIM/VERM(N)/4
       VERLI(N) = (log(VERM(N))+SIMI/4)/3 - 5.0D0/18.0D0
 300  CONTINUE ! end DO ... N = 1, 3
      RETURN
END SUBROUTINE S2D0TWOI
!C
!C ----------------------------------------------------------------------
!C file: s2d1onei.f
!C date: 2011-03-23
!C who:  S.Kaprzyk
!C what: CASE a) : [w1<w2<ONE]
!C       U(1)=W(1); U(2)=W(2)
!C LEPS(1)  [|w1-w2|<e]
SUBROUTINE S2D1ONEI(U,AS,LEPS, SIM1, SIM1I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM1, SIM1I
!C
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(2), DX0(4), R0(2), S0(2)
      DOUBLE COMPLEX   QX0(1), QR0(1), QS0(1)
      DOUBLE COMPLEX   AV, BV, DU, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/

      ABI_UNUSED((/as(1)/))
!C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
!C ------------------------------------------------------------------
!C      R0(U) = ARTH(U)/U
!C      R0(U) = 1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   The DX0(I) contains  derivatives of the function X0(u)      -
!C -----------------------------------------------------------------
!c      IF (LEPS(1).OR.LEPS(2).OR.LEPS(3)) THEN
      IF (LEPS(1)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!c      DO 200 I = 1, 3, 2
      I = 1
       DU = U(I) - U(2)
       IF (.NOT.LEPS(I)) THEN
        QX0(I) = (X0(I)-X0(2))/DU
       ELSE
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+ DX0(4)/24.0D0*DU)*DU)*DU
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
       QS0(I) = (-S0(2)-R0(2)+(CONE-U(I))*QR0(I))/(CONE+U(I))
!C
      AV = (CONE+U(1))*(CONE+U(2))
      BV = (U(2)-CONE)
      SIM1  = AV*(R0(1) + QR0(1)*BV)
      SIM1I = AV*(S0(1) + QS0(1)*BV)
      RETURN
END SUBROUTINE S2D1ONEI
!C
!C ----------------------------------------------------------------------
!C file: s2d1twoi.f
!C date: 2010-09-22
!C who:  S.Kaprzyk
!C what: CASE a) : [w1<w2<ONE]
!C       U(1)=W(1); U(2)=W(2)
!C LEPS(1)  [|w1-w2|<e]
SUBROUTINE S2D1TWOI(U,AS,LEPS, SIM, SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(2)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX X0(2), DX0(4), R3(2), S3(2)
      DOUBLE COMPLEX QR3(1), QX0(1)
      DOUBLE COMPLEX QS3(1)
      DOUBLE COMPLEX AV, BV, UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C

      ABI_UNUSED((/as/))
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CTWO  = (2.0D0,0.0D0)
!C -----------------------------------------------------------------
!C -                                                               -
!C -        R3(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
!C -    R3(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C -                                                               -
!C -----------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R3(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
       S3(I) = (CONE-U(I))/(CONE+U(I))*R3(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   The DX0(I) Contains  Derivatives Of The Function X0(W)      -
!C -----------------------------------------------------------------
      IF (LEPS(1)) THEN
        CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!C
!c      DO 200 I = 1, 3, 2
      I = 1
       UV = U(I) - U(2)
       IF (LEPS(I)) THEN
       QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+ DX0(4)/24.0D0*UV)*UV)*UV
       ELSE
        QX0(I) = (X0(I)-X0(2))/UV
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR3(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 + &
     &   G(3)/5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
       QS3(I) = (-S3(2)-R3(2)+(CONE-U(I))*QR3(I))/(CONE+U(I))
!c 200  CONTINUE  ! DO 750 I = 1, 3, 2
!c
      AV = (CONE+U(1))*(CONE+U(2))
      BV = U(2)-CONE
      SIM  = AV*(R3(1) + QR3(1)*BV)
      SIMI = AV*(S3(1) + QS3(1)*BV)
      RETURN
END SUBROUTINE S2D1TWOI
!C
!C ----------------------------------------------------------------------
!C file: s2d2onei.f
!C date: 2011-03-24
!C who:  S.Kaprzyk
!C what: CASE b) : [w1<ONE<w2]
!C       U(1)=W(1); U(2)= CONE/W(2)
!C LEPS(1) [|w1-w2|<e]
SUBROUTINE  S2D2ONEI(U,AS,LEPS, SIM2, SIM2I)

      IMPLICIT         NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM2, SIM2I
!c
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
      DOUBLE COMPLEX   X0(2), R0(2), S0(2) !, DX0(4)
!c      DOUBLE COMPLEX   QX0(1), QR0(1), QS0(1)
      DOUBLE COMPLEX   AV, BV, CV !, G(5)
      !DOUBLE COMPLEX   UD
      INTEGER          I !, K, N
      DOUBLE COMPLEX   CZERO, CONE, CTWO
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(leps)
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CTWO  = (2.0D0,0.0D0)
!C ------------------------------------------------------------------
!C      R0(U) = ARTH(U)/U
!C      R0(U) = 1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       S0(I) = (CONE-U(I))/(CONE+U(I))*R0(I)
 100  CONTINUE
!c
      AV = (CONE+U(1))*(U(2)+CONE)/(CONE-U(2)*U(1))
      BV = -(U(1)-CONE)
      CV = (CONE-U(2))
      SIM2 = AV*(BV*R0(1)+CV*(DCMPLX(ZERO,AS(2))+U(2)*R0(2)))
      CV = CV*(U(2)-CONE)/(U(2)+CONE)
      SIM2I= AV*(BV*S0(1)+CV*(DCMPLX(ZERO,AS(2))+U(2)*R0(2)))
!c
      RETURN
END SUBROUTINE S2D2ONEI
!C
!C ----------------------------------------------------------------------
!C file: s2d2twoi.f
!C date: 2011-03-26
!C who:  S.Kaprzyk
!C what: CASE b) : [w1<ONE<w2]
!C       U(1)=W(1); U(2)= CONE/W(2)
!C LEPS(1) [|w1-w2|<e]
      SUBROUTINE S2D2TWOI(U,AS,LEPS, SIM, SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX X0(2), R3(2), S3(2)
      DOUBLE COMPLEX AV, BV
      INTEGER          I
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C

      ABI_UNUSED((/leps/))
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
!C   -----------------------------------------------------------------
!C   -                                                               -
!C   -        R3(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
!C   -    R3(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C   -                                                               -
!C   -----------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R3(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
 100  CONTINUE
      S3(1) = (CONE-U(1))/(CONE+U(1))*R3(1)
      AV = (CONE+U(1))*(U(2)+CONE)/(U(1)*U(2)-CONE)
      BV = (CONE-U(2))*DCMPLX(ZERO,AS(2))+CONE+U(2)-U(2)*U(2)*R3(2)
      SIM = AV*((U(1)-CONE)*R3(1) - (CONE-U(2))*BV)
      BV = (U(2)-CONE)/(U(2)+CONE)*BV
      SIMI= AV*((U(1)-CONE)*S3(1) - (CONE-U(2))*BV)
!C
      RETURN
END SUBROUTINE S2D2TWOI
!C
!C ----------------------------------------------------------------------
!C file: s2d3onei.f
!C date: 2011-04-01
!C who:  S.Kaprzyk
!C what: CASE c) : [ ONE<w1<w2]
!C       U(1)=CONE/W(1); U(2)=CONE/W(2)
!C LEPS(1) [|w1-w2|<e]
SUBROUTINE  S2D3ONEI(U,AS,LEPS,SIM3,SIM3I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM3, SIM3I
!c
      !DOUBLE COMPLEX   SIM0UX0
      !EXTERNAL         SIM0UX0
!c
      DOUBLE COMPLEX   X0(2), QX0(1), DX0(4)
      DOUBLE COMPLEX   R0(2), UR0(2), QR0(1), QUR0(1)
      DOUBLE COMPLEX   US0(2), QUS0(1), QAS(1)
      DOUBLE COMPLEX   AV, UD,  G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
!C ------------------------------------------------------------------
!C       R0(U)=ARTH(U)/U
!C       R0(U)=1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)*U(I)
       X0(I) = SIM0UX0(U(I))
       R0(I) = CONE + AV*(CONE/3.0D0+AV*(CONE/5.0D0+AV*X0(I)))
       UR0(I)= DCMPLX(ZERO,AS(I)) + U(I)*R0(I)
 100  CONTINUE
!C -----------------------------------------------------------------
!C -   THE DX0(I) CONTAINS  DERIVATIVES OF THE FUNCTION X0(U)
!C -----------------------------------------------------------------
      IF (LEPS(1)) THEN
       CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!c
      I = 1
       UD = U(I) - U(2)
       IF (.NOT.LEPS(I)) THEN
        QX0(I) = (X0(I)-X0(2))/UD
        QAS(I) = (AS(I)-AS(2))/UD
       ELSE
         QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0 + DX0(4)/24.0D0*UD)*UD)*UD
        QAS(I) = CZERO
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR0(I) = G(1)/3.0D0 + G(3)/5.0D0 + G(5)*X0(2) + AV*U(I)*QX0(I)
!c
       QUR0(1) = CIMAG*QAS(1) + R0(2) +U(1)*QR0(1)
       US0(2) = (U(2)-CONE)/(U(2)+CONE)*UR0(2)
       QUS0(1) = ( -US0(2)+UR0(2)+(U(1)-CONE)*QUR0(1))/(U(1)+CONE)
!c
      AV = (U(1)+CONE)*(U(2)+CONE)
      SIM3  = AV*(UR0(2) + (U(1)-CONE)*QUR0(1))
      SIM3I = AV*(US0(2) + (U(1)-CONE)*QUS0(1))
      RETURN
END SUBROUTINE S2D3ONEI
!C
!C ----------------------------------------------------------------------
!C file: s2d3twoi.f
!C date: 2011-03-28
!C who:  S.Kaprzyk
!C what: CASE c) : [ ONE<w1<w2]
!C       U(1)=CONE/W(1); U(2)=CONE/W(2)
!C LEPS(1) [|w1-w2|<e]
SUBROUTINE S2D3TWOI(U,AS,LEPS,SIM,SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(2)
      DOUBLE PRECISION AS(2)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX X0(2), QX0(1), DX0(4)
      DOUBLE COMPLEX R3(2), UR3(2), QR3(1), QUR3(1)
      DOUBLE COMPLEX US3(2), QUS3(1)
      DOUBLE COMPLEX QAS(1), AV,  UV, G(5)
      INTEGER          I, K
      DOUBLE COMPLEX   CZERO, CONE, CTWO, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CTWO  = (2.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
!C   -----------------------------------------------------------------
!C   -                                                               -
!C   -        R3(U)=1+(U-1)*(ARTH(U)/U-1)/U                          -
!C   -    R3(U)=1-U/3+U**2/3-U**3/5+U**4/5+(U-1)*U**5*X0(U)          -
!C   -                                                               -
!C   -----------------------------------------------------------------
      DO 100 I = 1, 2
       AV = U(I)
       X0(I) = SIM0UX0(AV)
!*-ASIS
       R3(I) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(I)))))
       UR3(I)=(CONE-U(I))*DCMPLX(ZERO,AS(I))+CONE+U(I)-U(I)*U(I)*R3(I)
 100  CONTINUE  ! end DO ... I = 1, 2
!C -----------------------------------------------------------------
!C -   The DX0(I) Contains  Derivatives Of The Function X0(U)      -
!C -----------------------------------------------------------------
      IF (LEPS(1)) THEN
        CALL SIM0UDX0(U(2),X0(2),DX0)
      END IF
!C
       I = 1
       UV = U(I) - U(2)
       IF (LEPS(I)) THEN
        QX0(I) = DX0(1) + (DX0(2)/2.0D0+(DX0(3)/6.0D0+DX0(4)/24.0D0*UV)*UV)*UV
        QAS(I) = CZERO
       ELSE
        QX0(I) = (X0(I)-X0(2))/UV
        QAS(I) = (AS(I)-AS(2))/UV
       END IF
       AV = U(I)
       G(1) = U(I) + U(2)
       DO 150 K = 2, 5
        AV = AV*U(I)
        G(K) = G(K-1)*U(2) + AV
 150   CONTINUE
       QR3(I) = -CONE/3.0D0 + G(1)/3.0D0 - G(2)/5.0D0 + &
     &   G(3)/5.0D0 + (G(5)-G(4))*X0(2) + (U(I)-CONE)*AV*QX0(I)
       QUR3(1) = -DCMPLX(ZERO,AS(2))+CIMAG*(CONE-U(1))*QAS(1)+ CONE - &
     &            (U(1)+U(2))*R3(2)-U(1)*U(1)*QR3(1)
!C
       US3(2) = (U(2)-CONE)/(U(2)+CONE)*UR3(2)
       QUS3(1)=( -US3(2)+UR3(2)+(U(1)-CONE)*QUR3(1))/(U(1)+CONE)
!c
      AV = (U(1)+CONE)*(U(2)+CONE)
      SIM = AV*(UR3(2) + (U(1)-CONE)*QUR3(1))
      SIMI= AV*(US3(2) + (U(1)-CONE)*QUS3(1))
!c
      RETURN
END SUBROUTINE S2D3TWOI
!C
!C ----------------------------------------------------------------------
!C file: s2d0leps.f
!C date: 2011-04-02
!C who:  S.Kaprzyk
!C what: Return a proper CASE:
!C       a) [w1<w2<ONE]; b) [w1<ONE<w2]; c) [ONE<w1<w2]
!C INPUT:
!C VERM(3) - three complex numbers
!C N - integer number for W(I)=(VERM(I)-VERM(N))/(VERM(I)+VERM(N))
!C EPS  - critical distance |w1-w2| etc.
!C OUTPUT:
!C LONE(3) [w1<w2<ONE]; [w1<ONE<w2]; [ONE<w1<w2]
!C LEPS(1) [|w1-w2|<eps]
      SUBROUTINE  S2D0LEPS(VERM, N, W, AS, LONE, EPS, LEPS, iuerr)

      IMPLICIT       NONE
      DOUBLE COMPLEX   VERM(3)
      INTEGER          N
      DOUBLE COMPLEX   W(2)
      DOUBLE PRECISION AS(2)
      DOUBLE PRECISION EPS
      LOGICAL          LONE(3), LEPS(1)
      INTEGER          iuerr
!C
!c      DOUBLE PRECISION DLAMCH
!c      EXTERNAl         DLAMCH
!C
      LOGICAL          LWONE(2), LV
      DOUBLE COMPLEX   AV
      DOUBLE PRECISION AL(2), AA, AIW, OZERO
      INTEGER          I, II, K
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE, PI
      DATA             ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/
      DATA             OZERO/ 1.0D-13/
!c      OZERO = DLAMCH('e')*10
      CZERO = (0.0D0,0.0D0)
      CONE = (1.0D0,0.0D0)
      IF ((N.GT.3).OR.(N.LT.1)) THEN
       WRITE (iuerr,9010) N
 9010  FORMAT (' ***s2d0leps: N =',I6,' must be 1, or 2')
       STOP ' ***s2d0leps: '
      END IF
      II = 0
      DO 200 I = 1, 3
       IF (I.EQ.N) GO TO 200
       II = II + 1
       IF (CDABS(VERM(N)-VERM(I)).LT.CDABS(VERM(N)+VERM(I))) THEN
        LWONE(II) = .TRUE.
        W(II) = (VERM(N)-VERM(I))/(VERM(N)+VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = 0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       ELSE
        LWONE(II) = .FALSE.
        W(II) = (VERM(N)+VERM(I))/(VERM(N)-VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = -0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       END IF
       AL(II) = CDABS(W(II))
 200  CONTINUE
      DO 300 I = 1, 2
       DO 250 K = I, 2
        IF (LWONE(I).AND.LWONE(K).AND.(AL(I).LT.AL(K))) GO TO 250
        IF (.NOT.LWONE(I).AND..NOT.LWONE(K).AND.(AL(K).LT.AL(I))) GO TO 250
        IF (.NOT.LWONE(I).AND.LWONE(K).AND.(ONE.LT.AL(I)*AL(K))) GO TO 250
        IF (LWONE(I).AND..NOT.LWONE(K).AND.(AL(I)*AL(K).LT.ONE)) GO TO 250
        AA = AL(K)
        AL(K) = AL(I)
        AL(I) = AA
        AA = AS(K)
        AS(K) = AS(I)
        AS(I) = AA
        AV = W(K)
        W(K) = W(I)
        W(I) = AV
        LV = LWONE(K)
        LWONE(K) = LWONE(I)
        LWONE(I) = LV
 250   CONTINUE
 300  CONTINUE
!C
      LONE(1) = .FALSE.
      LONE(2) = .FALSE.
      LONE(3) = .FALSE.
      IF (LWONE(1).AND.LWONE(2)) THEN
       LONE(1) = .TRUE.
      END IF
      IF (LWONE(1).AND..NOT.LWONE(2)) THEN
       LONE(2) = .TRUE.
      END IF
      IF (.NOT.LWONE(1).AND..NOT.LWONE(2)) THEN
       LONE(3) = .TRUE.
      END IF
!C here are cases, how close are w()
      LEPS(1) = .FALSE.
      IF ( (LONE(1).OR.LONE(3)).AND. (CDABS(W(1)-W(2)).LT.EPS) ) LEPS(1)=.TRUE.
!c
      IF ((.NOT.LONE(1)).AND.(.NOT.LONE(2)).AND.(.NOT.LONE(3))) THEN
        WRITE (iuerr,9030) (LONE(I),I=1,3)
 9030   FORMAT ('***s2d0leps: LONE()=',3(L1,1X),'not in order')
        STOP '***s2d0leps: '
      END IF
!C
      RETURN
END SUBROUTINE S2D0LEPS
!C s2d0leps()
!C file: s1d0onei.f
!C date: 2011-03-20
!C who:  S.Kaprzyk
!C what: Complex linear form integral over standard 1-d simplex
!C ----------------------------------------------------
!C -         *                                        -
!C -        *                         1               -
!C - SIM0=  * dt1 ----------------------------------  -
!C -        *     verm(2)+(verm(1)-verm(2))*t1        -
!C -       *                                          -
!C -    0<t1<1                                        -
!C ----------------------------------------------------
SUBROUTINE S1D0ONEI(SIM0, SIM0I, VERM)

      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX SIM0, SIM0I
      DOUBLE COMPLEX VERM(*)
!C
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(2),LEPS(1)
!c
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      INTEGER          I, K, N
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, EPS
      DATA             EPS/1.0D-6/
      DATA             ZERO/0.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      N = 2
      DO 100 I = 1, 1
       IF (DIMAG(VERM(N))*DIMAG(VERM(I)).LT.ZERO) THEN
         WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,2)
 9010    FORMAT (' ***s1d0onei: not signed ImgVERM()=',4(D13.6,1X))
!c        STOP ' ***s1d0onei: '
       END IF
 100  CONTINUE
      DO 101 I = 1, 1
      IF (CDABS(VERM(I)).GT.CDABS(VERM(N))) N = I
 101  CONTINUE
!C Here are 2 cases, how w() are placed on complex-plane
!C LONE(1..2) [w1<ONE]; [ONE<w1]
!C LEPS(1) null
!*-ASIS
      CALL S1D0LEPS(VERM, N, U, AS, LONE, EPS, LEPS, iuerr)
!c U(1) = W(1);
      IF (LONE(1)) THEN
        CALL S1D1ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!c U(1) = CONE/W(1)
      IF (LONE(2)) THEN
        CALL S1D2ONEI(U,AS,LEPS,SIM,SIMI)
      END IF
!C
      SIM0 =  SIM/VERM(N)
      SIM0I= log(VERM(N)) + SIMI - CONE
      RETURN
END SUBROUTINE S1D0ONEI
!C
!C ----------------------------------------------------------------------
!C file: s1d0twoi.f
!C date: 2011-03-24
!C who:  S.Kaprzyk
!C what: Complex linear form integrals over standard 1d-simplex
!C ------------------(i=1)--------------------------------------
!C            *                                                -
!C           *            t_i                                  -
!C  VERL(i)= *dt1 ------------------------------------------   -
!C           *        VERM(2)+[VERM(1)-VERM(2)]*t1             -
!C          *                                                  -
!C     0<t1<1                                                  -
!C -------------------(i=2)-------------------------------------
!C            *                                                -
!C           *           1- t1                                 -
!C  VERL(2)= *dt1 ------------------------------------------   -
!C           *          VERM(2)+[VERM(1)-VERM(2)]*t1           -
!C          *                                                  -
!C     0<t1<1                                                  -
!C -------------------------------------------------------------
SUBROUTINE   S1D0TWOI(VERL, VERLI, VERM)

      IMPLICIT       NONE
      INTEGER        iuerr
      PARAMETER     (iuerr=6)
      DOUBLE COMPLEX VERL(*), VERLI(*), VERM(*)
!C
      DOUBLE COMPLEX   SIM, SIMI
      LOGICAL          LONE(2),LEPS(1)
!c
      INTEGER          I, K, N
      DOUBLE PRECISION AS(1)
      DOUBLE COMPLEX   U(1)
      DOUBLE COMPLEX   CZERO
      DOUBLE PRECISION EPS
      DOUBLE PRECISION ZERO
      DATA             EPS/1.0D-6/
      DATA             ZERO/0.0D0/
!C
      CZERO = (0.0D0,0.0D0)
      DO 50 I = 1, 1
        IF (DIMAG(VERM(2))*DIMAG(VERM(I)).LT.ZERO) THEN
         WRITE (iuerr,9010) (DIMAG(VERM(K)),K=1,2)
 9010    FORMAT (' ***s1d0twoi: not signed ImgVERM()=',2(D13.6,1X))
         STOP ' ***s1dt0woi: '
        END IF
 50    CONTINUE
      DO 100 N = 1, 2
      VERL(N) = CZERO
!C Here are 2 cases, how w() are placed on complex-plane
!C LONE(1..2) [w1<ONE]; [ONE<w1]
!C LEPS(1) null
!*-ASIS
       CALL S1D0LEPS(VERM, N, U, AS, LONE, EPS, LEPS, iuerr)
!c U(1) = W(1);
       IF (LONE(1)) THEN
        CALL S1D1TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c U(1) = CONE/W(1);
       IF (LONE(2)) THEN
        CALL S1D2TWOI(U,AS,LEPS,SIM,SIMI)
       END IF
!c
       VERL(N)  = 0.50D0*SIM/VERM(N)
       VERLI(N) = 0.50D0*log(VERM(N)) + 0.25D0*SIMI -0.25D0
 100  CONTINUE ! end DO .. N = 1, 2
      RETURN
END SUBROUTINE S1D0TWOI
!C
!C ----------------------------------------------------------------------
!C file: s1d1onei.f
!C date: 2010-09-22
!C who:  S.Kaprzyk
!C what: CASE a) : [w1<ONE]
!C       U(1)=W(1)
!C LEPS() null
      SUBROUTINE S1D1ONEI(U,AS,LEPS, SIM1, SIM1I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM1, SIM1I
!C
      !DOUBLE COMPLEX   SIM0UR0
      !EXTERNAL         SIM0UR0
      DOUBLE COMPLEX   R0(1)
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(as)
      ABI_UNUSED(leps)
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
!C ------------------------------------------------------------------
!C      R0(U) = ARTH(U)/U
!C      R0(U) = 1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      R0(1) = SIM0UR0(U(1))
      SIM1  = (CONE+U(1))*R0(1)
      SIM1I = (CONE-U(1))*R0(1)
      RETURN
      END SUBROUTINE S1D1ONEI
!C
!C ----------------------------------------------------------------------
!C file: s1d1twoi.f
!C date: 2011-03-25
!C who:  S.Kaprzyk
!C what: CASE a) : [w1<ONE]
!C       U(1)=W(1)
!C LEPS() null
SUBROUTINE S1D1TWOI(U,AS,LEPS, SIM, SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX X0(1), R2(1), S2(1)
      DOUBLE COMPLEX AV
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(as)
      ABI_UNUSED(leps)
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
!C -----------------------------------------------------------------
!C -                                                               -
!C -    R2(W)=1+(W-1)*(ARTH(W)/W-1)/W=[1+(W-1)*R0(W)]/W            -
!C -    R2(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C -                                                               -
!C -----------------------------------------------------------------
       AV = U(1)
       X0(1) = SIM0UX0(AV)
!*-ASIS
       R2(1) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &   AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(1)))))
       S2(1) = (CONE-U(1))/(CONE+U(1))*R2(1)
      SIM  = (CONE+U(1))*R2(1)
      SIMI = (CONE+U(1))*S2(1)
      RETURN
END SUBROUTINE S1D1TWOI
!C
!C ----------------------------------------------------------------------
!C file: s1d2onei.f
!C date: 2011-03-20
!C who:  S.Kaprzyk
!C what: CASE b) : [ ONE<w1]
!C       U(1) = CONE/W(1)
!C LEPS(1) null
SUBROUTINE  S1D2ONEI(U,AS,LEPS,SIM2,SIM2I)

      IMPLICIT       NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM2, SIM2I
!c
      !DOUBLE COMPLEX   SIM0UR0
      !EXTERNAL         SIM0UR0
      DOUBLE COMPLEX   R0(1)
      DOUBLE COMPLEX   CZERO, CONE, CIMAG
      DOUBLE PRECISION ZERO, ONE
      DATA             ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(leps)
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      CIMAG = (0.0D0,1.0D0)
!C ------------------------------------------------------------------
!C       R0(U)=ARTH(U)/U
!C       R0(U)=1+U**2/3+U**6*X0(U)
!C ------------------------------------------------------------------
      R0(1) = SIM0UR0(U(1))
      SIM2 = (U(1)+CONE)*(CIMAG*AS(1)+U(1)*R0(1))
      SIM2I= (U(1)-CONE)*(CIMAG*AS(1)+U(1)*R0(1))
!c
      RETURN
END SUBROUTINE S1D2ONEI
!C
!C ----------------------------------------------------------------------
!C s1d2onei()
!C file: s1d2twoi.f
!C date: 2011-03-25
!C who:  S.Kaprzyk
!C what: CASE b) : [ ONE<w1]
!C       U(1) = CONE/W(1)
!C LEPS(1) null
SUBROUTINE S1D2TWOI(U,AS,LEPS, SIM, SIMI)

      IMPLICIT   NONE
      DOUBLE COMPLEX   U(1)
      DOUBLE PRECISION AS(1)
      LOGICAL          LEPS(1)
      DOUBLE COMPLEX   SIM, SIMI
!C
      !DOUBLE COMPLEX SIM0UX0
      !EXTERNAL       SIM0UX0
!c
      DOUBLE COMPLEX X0(2), R2(2), S2(2)
      DOUBLE COMPLEX AV
      !INTEGER        I
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE
      DATA  ZERO/0.0D0/, ONE/1.0D0/
!C
      ABI_UNUSED(leps)
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
!C   -----------------------------------------------------------------
!C   -                                                               -
!C   -        R2(W)=1+(W-1)*(ARTH(W)/W-1)/W                          -
!C   -    R2(W)=1-W/3+W**2/3-W**3/5+W**4/5+(W-1)*W**5*X0(W)          -
!C   -                                                               -
!C   -----------------------------------------------------------------
       AV = U(1)
       X0(1) = SIM0UX0(AV)
!*-ASIS
       R2(1) = CONE + AV*(-CONE/3.0D0 + AV*(CONE/3.0D0 + &
     &    AV*((-0.2D0,0.0D0)+AV*((0.2D0,0.0D0)+AV*(AV-CONE)*X0(1)))))
!C
      R2(1)= (CONE-U(1))*DCMPLX(ZERO,AS(1))+CONE+U(1)-U(1)*U(1)*R2(1)
      S2(1) = (U(1)-CONE)/(U(1)+CONE)*R2(1)
      SIM  = (U(1)+CONE)*R2(1)
      SIMI = (U(1)+CONE)*S2(1)
!C
      RETURN
END SUBROUTINE S1D2TWOI
!C
!C ----------------------------------------------------------------------
!C file: s1d0leps.f
!C date: 2011-04-02
!C who:  S.Kaprzyk
!C what: Return a proper case:
!C       a) [w1<ONE]; b) [ONE<w1]
!C INPUT:
!C VERM(1,2) - two complex numbers
!C N - integer number for W(I)=(VERM(I)-VERM(N))/(VERM(I)+VERM(N))
!C EPS  - critical distance
!C OUTPUT:
!C LONE(2) [w1<ONE]; [ONE<w1]
!C LEPS(1) null
      SUBROUTINE  S1D0LEPS(VERM, N, W, AS, LONE, EPS, LEPS, iuerr)

      IMPLICIT       NONE
      DOUBLE COMPLEX   VERM(2)
      INTEGER          N
      DOUBLE COMPLEX   W(1)
      DOUBLE PRECISION AS(1)
      DOUBLE PRECISION EPS
      LOGICAL          LONE(2), LEPS(1)
      INTEGER          iuerr
!C
!c      DOUBLE PRECISION DLAMCH
!c      EXTERNAl         DLAMCH
!C
      LOGICAL          LWONE(1)
      DOUBLE PRECISION AL(1), AIW, OZERO
      INTEGER          I, II !, K
      DOUBLE COMPLEX   CZERO, CONE
      DOUBLE PRECISION ZERO, ONE, PI
      DATA             ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/
      DATA             OZERO/ 1.0D-13/
      ABI_UNUSED(eps)
!c      OZERO = DLAMCH('e')*10
      CZERO = (0.0D0,0.0D0)
      CONE  = (1.0D0,0.0D0)
      IF ((N.GT.2).OR.(N.LT.1)) THEN
       WRITE (iuerr,9010) N
 9010  FORMAT (' ***s1d0leps: N =',I6,' must be 1, or 2')
       STOP ' ***s1d0leps: '
      END IF
!C
      II = 0
      DO 200 I = 1, 2
       IF (I.EQ.N) GO TO 200
       II = II + 1
       IF (CDABS(VERM(N)-VERM(I)).LT.CDABS(VERM(N)+VERM(I))) THEN
        LWONE(II) = .TRUE.
        W(II) = (VERM(N)-VERM(I))/(VERM(N)+VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = 0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       ELSE
        LWONE(II) = .FALSE.
        W(II) = (VERM(N)+VERM(I))/(VERM(N)-VERM(I))
        AIW = DIMAG(W(II))
        IF (DABS(AIW).GE.OZERO) THEN
        AS(II) = -0.5D0*PI*DSIGN(ONE,AIW)
        ELSE
        AS(II) = 0.5D0*PI*DSIGN(ONE,DREAL(VERM(I)-VERM(N)))
        END IF
       END IF
       AL(II) = CDABS(W(II))
 200  CONTINUE
!C
      LONE(1) = .FALSE.
      LONE(2) = .FALSE.
      IF (LWONE(1)) THEN
       LONE(1) = .TRUE.
      END IF
      IF (.NOT.LWONE(1)) THEN
       LONE(2) = .TRUE.
      END IF
!C here are cases, how close are w()
      LEPS(1) = .FALSE.
!c
      IF ((.NOT.LONE(1)).AND.(.NOT.LONE(2))) THEN
       WRITE (iuerr,9030) (LONE(I),I=1,2)
 9030  FORMAT ('***s1d0leps: LONE()=',2(L1,1X),'not in order')
       STOP '***s1d0leps: '
      END IF
!c
      RETURN
END SUBROUTINE S1D0LEPS
!C

end module m_simtet
!!***

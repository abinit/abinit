#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

SUBROUTINE CALJY0(ARG,RESULT,JINT)

!---------------------------------------------------------------------
!
! This packet computes zero-order Bessel functions of the first and
!   second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!   for Y0, and |X| <= XMAX for J0.  It contains two function-type
!   subprograms,  BESJ0  and  BESY0,  and one subroutine-type
!   subprogram,  CALJY0.  The calling statements for the primary
!   entries are:
!
!           Y = BESJ0(X)
!   and
!           Y = BESY0(X),
!
!   where the entry points correspond to the functions J0(X) and Y0(X),
!   respectively.  The routine  CALJY0  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY0  with
!   the statement
!           CALL CALJY0(ARG,RESULT,JINT),
!   where the parameter usage is as follows:
!
!      Function                  Parameters for CALJY0
!       call              ARG             RESULT          JINT
!
!     BESJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
!     BESY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1
!
!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the 
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ0(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1. 
!   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).
!
!     Approximate values for some important machines are
!
!                          eps      XMAX     XSMALL      XINF  
!
!  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
!
!*******************************************************************
!*******************************************************************
!
! Error Returns
!
!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY0 is called with a negative or zero argument.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, COS, LOG, SIN, SQRT
!
!
!  Latest modification: June 2, 1989
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division 
!          Argonne National Laboratory
!          Argonne, IL 60439
! 
!  Taken from http://www.netlib.org/specfun/j0y0
!
!--------------------------------------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'CALJY0'
!End of the abilint section

      IMPLICIT NONE
      INTEGER :: I,JINT
!CS    REAL
       DOUBLE PRECISION :: &
&            ARG,AX,CONS,DOWN,EIGHT,FIVE5,FOUR,ONE,ONEOV8,PI2,PJ0,&
&            PJ1,PLG,PROD,PY0,PY1,PY2,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,&
&            QY2,Q0,Q1,RESJ,RESULT,R0,R1,SIXTY4,THREE,TWOPI,TWOPI1,&
&            TWOPI2,TWO56,UP,W,WSQ,XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,&
&            XJ1,XJ01,XJ02,XJ11,XJ12,XY,XY0,XY01,XY02,XY1,XY11,XY12,&
&            XY2,XY21,XY22,Z,ZERO,ZSQ
      DIMENSION :: PJ0(7),PJ1(8),PLG(4),PY0(6),PY1(7),PY2(8),P0(6),P1(6),&
&               QJ0(5),QJ1(7),QLG(4),QY0(5),QY1(6),QY2(7),Q0(5),Q1(5)
!-------------------------------------------------------------------
!  Mathematical constants
!    CONS = ln(.5) + Euler's gamma
!-------------------------------------------------------------------
!CS    DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0E0,1.0E0,3.0E0,4.0E0,8.0E0/,
!CS   1     FIVE5,SIXTY4,ONEOV8,P17/5.5E0,64.0E0,0.125E0,1.716E-1/,
!CS   2     TWO56,CONS/256.0E0,-1.1593151565841244881E-1/,
!CS   3     PI2,TWOPI/6.3661977236758134308E-1,6.2831853071795864769E0/,
!CS   4     TWOPI1,TWOPI2/6.28125E0,1.9353071795864769253E-3/
     DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0D0,1.0D0,3.0D0,4.0D0,8.0D0/,&
&        FIVE5,SIXTY4,ONEOV8,P17/5.5D0,64.0D0,0.125D0,1.716D-1/,&
&        TWO56,CONS/256.0D0,-1.1593151565841244881D-1/,&
&        PI2,TWOPI/6.3661977236758134308D-1,6.2831853071795864769D0/,&
&        TWOPI1,TWOPI2/6.28125D0,1.9353071795864769253D-3/
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
!CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-09/,XINF/1.7E+38/
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
!CS    DATA XJ0/2.4048255576957727686E+0/,XJ1/5.5200781102863106496E+0/,
!CS   1     XY0/8.9357696627916752158E-1/,XY1/3.9576784193148578684E+0/,
!CS   2     XY2/7.0860510603017726976E+0/,
!CS   3     XJ01/ 616.0E+0/, XJ02/-1.4244423042272313784E-03/,
!CS   4     XJ11/1413.0E+0/, XJ12/ 5.4686028631064959660E-04/,
!CS   5     XY01/ 228.0E+0/, XY02/ 2.9519662791675215849E-03/,
!CS   6     XY11/1013.0E+0/, XY12/ 6.4716931485786837568E-04/,
!CS   7     XY21/1814.0E+0/, XY22/ 1.1356030177269762362E-04/
      DATA XJ0/2.4048255576957727686D+0/,XJ1/5.5200781102863106496D+0/,&
&          XY0/8.9357696627916752158D-1/,XY1/3.9576784193148578684D+0/,&
&          XY2/7.0860510603017726976D+0/,&
&          XJ01/ 616.0D+0/, XJ02/-1.4244423042272313784D-03/,&
&          XJ11/1413.0D+0/, XJ12/ 5.4686028631064959660D-04/,&
&          XY01/ 228.0D+0/, XY02/ 2.9519662791675215849D-03/,&
&          XY11/1013.0D+0/, XY12/ 6.4716931485786837568D-04/,&
&          XY21/1814.0D+0/, XY22/ 1.1356030177269762362D-04/
!C-------------------------------------------------------------------
!C  Coefficients for rational approximation to ln(x/a)
!C--------------------------------------------------------------------
!CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
!CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
!CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
!CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
       DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,&
&               -5.4989956895857911039D+02,3.5687548468071500413D+02/
       DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,&
&               -3.3442903192607538956D+02,1.7843774234035750207D+02/
!C-------------------------------------------------------------------
!C  Coefficients for rational approximation of
!C  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
!C--------------------------------------------------------------------
!CS    DATA PJ0/6.6302997904833794242E+06,-6.2140700423540120665E+08,
!CS   1         2.7282507878605942706E+10,-4.1298668500990866786E+11,
!CS   2        -1.2117036164593528341E-01, 1.0344222815443188943E+02,
!CS   3        -3.6629814655107086448E+04/
!CS    DATA QJ0/4.5612696224219938200E+05, 1.3985097372263433271E+08,
!CS   1         2.6328198300859648632E+10, 2.3883787996332290397E+12,
!CS   2         9.3614022392337710626E+02/
    DATA PJ0/6.6302997904833794242D+06,-6.2140700423540120665D+08,&
&            2.7282507878605942706D+10,-4.1298668500990866786D+11,&
&           -1.2117036164593528341D-01, 1.0344222815443188943D+02,&
&           -3.6629814655107086448D+04/
    DATA QJ0/4.5612696224219938200D+05, 1.3985097372263433271D+08,&
&            2.6328198300859648632D+10, 2.3883787996332290397D+12,&
&            9.3614022392337710626D+02/
!C-------------------------------------------------------------------
!C  Coefficients for rational approximation of
!C  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
!C-------------------------------------------------------------------
!CS    DATA PJ1/4.4176707025325087628E+03, 1.1725046279757103576E+04,
!CS   1         1.0341910641583726701E+04,-7.2879702464464618998E+03,
!CS   2        -1.2254078161378989535E+04,-1.8319397969392084011E+03,
!CS   3         4.8591703355916499363E+01, 7.4321196680624245801E+02/
!CS    DATA QJ1/3.3307310774649071172E+02,-2.9458766545509337327E+03,
!CS   1         1.8680990008359188352E+04,-8.4055062591169562211E+04,
!CS   2         2.4599102262586308984E+05,-3.5783478026152301072E+05,
!CS   3        -2.5258076240801555057E+01/
    DATA PJ1/4.4176707025325087628D+03, 1.1725046279757103576D+04,&
&            1.0341910641583726701D+04,-7.2879702464464618998D+03,&
&           -1.2254078161378989535D+04,-1.8319397969392084011D+03,&
&            4.8591703355916499363D+01, 7.4321196680624245801D+02/
    DATA QJ1/3.3307310774649071172D+02,-2.9458766545509337327D+03,&
&            1.8680990008359188352D+04,-8.4055062591169562211D+04,&
&            2.4599102262586308984D+05,-3.5783478026152301072D+05,&
&           -2.5258076240801555057D+01/
!C-------------------------------------------------------------------
!C  Coefficients for rational approximation of
!C    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
!C        XSMALL  <  |X|  <=  3.0
!C--------------------------------------------------------------------
!CS    DATA PY0/1.0102532948020907590E+04,-2.1287548474401797963E+06,
!CS   1         2.0422274357376619816E+08,-8.3716255451260504098E+09,
!CS   2         1.0723538782003176831E+11,-1.8402381979244993524E+01/
!CS    DATA QY0/6.6475986689240190091E+02, 2.3889393209447253406E+05,
!CS   1         5.5662956624278251596E+07, 8.1617187777290363573E+09,
!CS   2         5.8873865738997033405E+11/
     DATA PY0/1.0102532948020907590D+04,-2.1287548474401797963D+06,&
&             2.0422274357376619816D+08,-8.3716255451260504098D+09,&
&             1.0723538782003176831D+11,-1.8402381979244993524D+01/
     DATA QY0/6.6475986689240190091D+02, 2.3889393209447253406D+05,&
&             5.5662956624278251596D+07, 8.1617187777290363573D+09,&
&             5.8873865738997033405D+11/
!C-------------------------------------------------------------------
!C  Coefficients for rational approximation of
!C    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
!C        3.0  <  |X|  <=  5.5
!C--------------------------------------------------------------------
!CS    DATA PY1/-1.4566865832663635920E+04, 4.6905288611678631510E+06,
!CS   1         -6.9590439394619619534E+08, 4.3600098638603061642E+10,
!CS   2         -5.5107435206722644429E+11,-2.2213976967566192242E+13,
!CS   3          1.7427031242901594547E+01/
!CS    DATA QY1/ 8.3030857612070288823E+02, 4.0669982352539552018E+05,
!CS   1          1.3960202770986831075E+08, 3.4015103849971240096E+10,
!CS   2          5.4266824419412347550E+12, 4.3386146580707264428E+14/
      DATA PY1/-1.4566865832663635920D+04, 4.6905288611678631510D+06,&
&              -6.9590439394619619534D+08, 4.3600098638603061642D+10,&
&              -5.5107435206722644429D+11,-2.2213976967566192242D+13,&
&               1.7427031242901594547D+01/
      DATA QY1/ 8.3030857612070288823D+02, 4.0669982352539552018D+05,&
&               1.3960202770986831075D+08, 3.4015103849971240096D+10,&
&               5.4266824419412347550D+12, 4.3386146580707264428D+14/
!C-------------------------------------------------------------------
!C  Coefficients for rational approximation of
!C    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
!C        5.5  <  |X|  <=  8.0
!C--------------------------------------------------------------------
!CS    DATA PY2/ 2.1363534169313901632E+04,-1.0085539923498211426E+07,
!CS   1          2.1958827170518100757E+09,-1.9363051266772083678E+11,
!CS   2         -1.2829912364088687306E+11, 6.7016641869173237784E+14,
!CS   3         -8.0728726905150210443E+15,-1.7439661319197499338E+01/
!CS    DATA QY2/ 8.7903362168128450017E+02, 5.3924739209768057030E+05,
!CS   1          2.4727219475672302327E+08, 8.6926121104209825246E+10,
!CS   2          2.2598377924042897629E+13, 3.9272425569640309819E+15,
!CS   3          3.4563724628846457519E+17/
    DATA PY2/ 2.1363534169313901632D+04,-1.0085539923498211426D+07,&
&             2.1958827170518100757D+09,-1.9363051266772083678D+11,&
&            -1.2829912364088687306D+11, 6.7016641869173237784D+14,&
&            -8.0728726905150210443D+15,-1.7439661319197499338D+01/
    DATA QY2/ 8.7903362168128450017D+02, 5.3924739209768057030D+05,&
&             2.4727219475672302327D+08, 8.6926121104209825246D+10,&
&             2.2598377924042897629D+13, 3.9272425569640309819D+15,&
&             3.4563724628846457519D+17/
!C-------------------------------------------------------------------
!C  Coefficients for Hart,s approximation,  |X| > 8.0
!C-------------------------------------------------------------------
!CS    DATA P0/3.4806486443249270347E+03, 2.1170523380864944322E+04,
!CS   1        4.1345386639580765797E+04, 2.2779090197304684302E+04,
!CS   2        8.8961548424210455236E-01, 1.5376201909008354296E+02/
!CS    DATA Q0/3.5028735138235608207E+03, 2.1215350561880115730E+04,
!CS   1        4.1370412495510416640E+04, 2.2779090197304684318E+04,
!CS   2        1.5711159858080893649E+02/
!CS    DATA P1/-2.2300261666214198472E+01,-1.1183429920482737611E+02,
!CS   1        -1.8591953644342993800E+02,-8.9226600200800094098E+01,
!CS   2        -8.8033303048680751817E-03,-1.2441026745835638459E+00/
!CS    DATA Q1/1.4887231232283756582E+03, 7.2642780169211018836E+03,
!CS   1        1.1951131543434613647E+04, 5.7105024128512061905E+03,
!CS   2        9.0593769594993125859E+01/
    DATA P0/3.4806486443249270347D+03, 2.1170523380864944322D+04,&
&           4.1345386639580765797D+04, 2.2779090197304684302D+04,&
&           8.8961548424210455236D-01, 1.5376201909008354296D+02/
    DATA Q0/3.5028735138235608207D+03, 2.1215350561880115730D+04,&
&           4.1370412495510416640D+04, 2.2779090197304684318D+04,&
&           1.5711159858080893649D+02/
    DATA P1/-2.2300261666214198472D+01,-1.1183429920482737611D+02,&
&           -1.8591953644342993800D+02,-8.9226600200800094098D+01,&
&           -8.8033303048680751817D-03,-1.2441026745835638459D+00/
    DATA Q1/1.4887231232283756582D+03, 7.2642780169211018836D+03,&
&           1.1951131543434613647D+04, 5.7105024128512061905D+03,&
&           9.0593769594993125859D+01/
!C-------------------------------------------------------------------
!C  Check for error conditions
!C-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. (ARG .LE. ZERO)) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) GO TO 800
      IF (AX .LE. XSMALL) THEN
         IF (JINT .EQ. 0) THEN
               RESULT = ONE
            ELSE
               RESULT = PI2 * (LOG(AX) + CONS)
         END IF
         GO TO 2000
      END IF
!C-------------------------------------------------------------------
!C  Calculate J0 for appropriate interval, preserving
!C     accuracy near the zero of J0
!C-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(5) * ZSQ + PJ0(6)) * ZSQ + PJ0(7)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            WSQ = ONE - ZSQ / SIXTY4
            XNUM = PJ1(7) * WSQ + PJ1(8)
            XDEN = WSQ + QJ1(7)
            DO 220 I = 1, 6
               XNUM = XNUM * WSQ + PJ1(I)
               XDEN = XDEN * WSQ + QJ1(I)
  220       CONTINUE
            PROD = (AX + XJ1) * ((AX - XJ11/TWO56) - XJ12)
      END IF
      RESULT = PROD * XNUM / XDEN
      IF (JINT .EQ. 0) GO TO 2000
!C-------------------------------------------------------------------
!C  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!C    where xn is a zero of Y0
!C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE IF (AX .LE. FIVE5) THEN
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
         ELSE
            UP = (AX-XY21/TWO56)-XY22
            XY = XY2
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
!C-------------------------------------------------------------------
!C  Now calculate Y0 for appropriate interval, preserving
!C     accuracy near the zero of Y0
!C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            XNUM = PY0(6) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 5
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE IF (AX .LE. FIVE5) THEN
            XNUM = PY1(7) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 6
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
         ELSE
            XNUM = PY2(8) * ZSQ + PY2(1)
            XDEN = ZSQ + QY2(1)
            DO 380 I = 2, 7
               XNUM = XNUM * ZSQ + PY2(I)
               XDEN = XDEN * ZSQ + QY2(I)
  380       CONTINUE
      END IF
      RESULT = RESJ + UP * DOWN * XNUM / XDEN
      GO TO 2000
!C-------------------------------------------------------------------
!C  Calculate J0 or Y0 for |ARG|  >  8.0
!C-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AX / TWOPI
      W = AINT(W) + ONEOV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(5) * ZSQ + P0(6)
      XDEN = ZSQ + Q0(5)
      UP = P1(5) * ZSQ + P1(6)
      DOWN = ZSQ + Q1(5)
      DO 850 I = 1, 4
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = SQRT(PI2/AX) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = SQRT(PI2/AX) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
 2000 RETURN
!C---------- Last line of CALJY0 ----------
      END
      DOUBLE PRECISION FUNCTION BESJ0(X)
!CS    REAL FUNCTION BESJ0(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for Bessel functions
!   of the first kind of order zero for arguments  |X| <= XMAX
!   (see comments heading CALJY0).
!
!--------------------------------------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'BESJ0'
 use interfaces_28_numeric_noabirule, except_this_one => BESJ0
!End of the abilint section

      IMPLICIT NONE
      INTEGER :: JINT
!S    REAL  X, RESULT
      DOUBLE PRECISION :: X, RESULT
!--------------------------------------------------------------------
      JINT=0
      CALL CALJY0(X,RESULT,JINT)
      BESJ0 = RESULT
      RETURN
!---------- Last line of BESJ0 ----------
      END
      DOUBLE PRECISION FUNCTION BESY0(X)
!CS    REAL FUNCTION BESY0(X)
!C--------------------------------------------------------------------
!C
!C This subprogram computes approximate values for Bessel functions
!C   of the second kind of order zero for arguments 0 < X <= XMAX
!C   (see comments heading CALJY0).
!C
!C--------------------------------------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'BESY0'
 use interfaces_28_numeric_noabirule, except_this_one => BESY0
!End of the abilint section

      IMPLICIT NONE
      INTEGER JINT
!CS    REAL  X, RESULT
      DOUBLE PRECISION :: X, RESULT
!--------------------------------------------------------------------
      JINT=1
      CALL CALJY0(X,RESULT,JINT)
      BESY0 = RESULT
      RETURN
!---------- Last line of BESY0 ----------
      END

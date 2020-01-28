!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_bessel
!! NAME
!! m_bessel
!!
!! FUNCTION
!! Bessel functions
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_bessel

 implicit none

 private

 public :: CALJY0
 public :: CALJY1
 public :: CALCK0
 public :: CALCK1

contains
!!***

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
      END subroutine caljy0

      DOUBLE PRECISION FUNCTION BESJ0(X)
!CS    REAL FUNCTION BESJ0(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for Bessel functions
!   of the first kind of order zero for arguments  |X| <= XMAX
!   (see comments heading CALJY0).
!
!--------------------------------------------------------------------
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
      END function besj0
      DOUBLE PRECISION FUNCTION BESY0(X)
!CS    REAL FUNCTION BESY0(X)
!C--------------------------------------------------------------------
!C
!C This subprogram computes approximate values for Bessel functions
!C   of the second kind of order zero for arguments 0 < X <= XMAX
!C   (see comments heading CALJY0).
!C
!C--------------------------------------------------------------------
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
      END function besy0

SUBROUTINE CALJY1(ARG,RESULT,JINT)

!---------------------------------------------------------------------
!
! This packet computes first-order Bessel functions of the first and
!   second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!   for Y1, and |X| <= XMAX for J1.  It contains two function-type
!   subprograms,  BESJ1  and  BESY1,  and one subroutine-type
!   subprogram,  CALJY1.  The calling statements for the primary
!   entries are:
!
!           Y = BESJ1(X)
!   and
!           Y = BESY1(X),
!
!   where the entry points correspond to the functions J1(X) and Y1(X),
!   respectively.  The routine  CALJY1  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY1  with
!   the statement
!           CALL CALJY1(ARG,RESULT,JINT),
!   where the parameter usage is as follows:
!
!      Function                  Parameters for CALJY1
!       call              ARG             RESULT          JINT
!
!     BESJ1(ARG)     |ARG| .LE. XMAX       J1(ARG)          0
!     BESY1(ARG)   0 .LT. ARG .LE. XMAX    Y1(ARG)          1
!
!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons,
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ1(X), but cannot match that version for
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
!   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
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
!    -XINF when BESLY1 is called with a negative or zero argument.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, COS, LOG, SIN, SQRT
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: November 10, 1987
!
!  Taken from http://www.netlib.org/specfun/j1y1
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,JINT
      DIMENSION :: PJ0(7),PJ1(8),PLG(4),PY0(7),PY1(9),P0(6),P1(6),&
&                  QJ0(5),QJ1(7),QLG(4),QY0(6),QY1(8),Q0(6),Q1(6)
!CS    REAL
       DOUBLE PRECISION :: &
&        ARG,AX,DOWN,EIGHT,FOUR,HALF,PI2,PJ0,PJ1,PLG,PROD,PY0,&
&        PY1,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,Q0,Q1,RESJ,RESULT,&
&        RTPI2,R0,R1,THROV8,TWOPI,TWOPI1,TWOPI2,TWO56,UP,W,WSQ,&
&        XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,XJ1,XJ01,XJ02,XJ11,XJ12,&
&        XY,XY0,XY01,XY02,XY1,XY11,XY12,Z,ZERO,ZSQ
!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
!CS    DATA EIGHT/8.0E0/,
!CS   1     FOUR/4.0E0/,HALF/0.5E0/,THROV8/0.375E0/,
!CS   2     PI2/6.3661977236758134308E-1/,P17/1.716E-1/
!CS   3     TWOPI/6.2831853071795864769E+0/,ZERO/0.0E0/,
!CS   4     TWOPI1/6.28125E0/,TWOPI2/1.9353071795864769253E-03/
!CS   5     TWO56/256.0E+0/,RTPI2/7.9788456080286535588E-1/
      DATA EIGHT/8.0D0/,&
&          FOUR/4.0D0/,HALF/0.5D0/,THROV8/0.375D0/,&
&          PI2/6.3661977236758134308D-1/,P17/1.716D-1/&
&          TWOPI/6.2831853071795864769D+0/,ZERO/0.0D0/,&
&          TWOPI1/6.28125D0/,TWOPI2/1.9353071795864769253D-03/&
&          TWO56/256.0D+0/,RTPI2/7.9788456080286535588D-1/
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
!CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-09/,XINF/1.7E+38/
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
!CS    DATA XJ0/3.8317059702075123156E+0/,XJ1/7.0155866698156187535E+0/,
!CS   1     XY0/2.1971413260310170351E+0/,XY1/5.4296810407941351328E+0/,
!CS   2     XJ01/ 981.0E+0/, XJ02/-3.2527979248768438556E-04/,
!CS   3     XJ11/1796.0E+0/, XJ12/-3.8330184381246462950E-05/,
!CS   4     XY01/ 562.0E+0/, XY02/ 1.8288260310170351490E-03/,
!CS   5     XY11/1390.0E+0/, XY12/-6.4592058648672279948E-06/
      DATA XJ0/3.8317059702075123156D+0/,XJ1/7.0155866698156187535D+0/,&
&          XY0/2.1971413260310170351D+0/,XY1/5.4296810407941351328D+0/,&
&          XJ01/ 981.0D+0/, XJ02/-3.2527979248768438556D-04/,&
&          XJ11/1796.0D+0/, XJ12/-3.8330184381246462950D-05/,&
&          XY01/ 562.0D+0/, XY02/ 1.8288260310170351490D-03/,&
&          XY11/1390.0D+0/, XY12/-6.4592058648672279948D-06/
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
!CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
!CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
!CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
!CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,&
&              -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,&
&              -3.3442903192607538956D+02,1.7843774234035750207D+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
!CS    DATA PJ0/9.8062904098958257677E+05,-1.1548696764841276794E+08,
!CS   1       6.6781041261492395835E+09,-1.4258509801366645672E+11,
!CS   2      -4.4615792982775076130E+03, 1.0650724020080236441E+01,
!CS   3      -1.0767857011487300348E-02/
!CS    DATA QJ0/5.9117614494174794095E+05, 2.0228375140097033958E+08,
!CS   1       4.2091902282580133541E+10, 4.1868604460820175290E+12,
!CS   2       1.0742272239517380498E+03/
      DATA PJ0/9.8062904098958257677D+05,-1.1548696764841276794D+08,&
&              6.6781041261492395835D+09,-1.4258509801366645672D+11,&
&             -4.4615792982775076130D+03, 1.0650724020080236441D+01,&
&             -1.0767857011487300348D-02/
      DATA QJ0/5.9117614494174794095D+05, 2.0228375140097033958D+08,&
&              4.2091902282580133541D+10, 4.1868604460820175290D+12,&
&              1.0742272239517380498D+03/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
!CS    DATA PJ1/4.6179191852758252280E+00,-7.1329006872560947377E+03,
!CS   1       4.5039658105749078904E+06,-1.4437717718363239107E+09,
!CS   2       2.3569285397217157313E+11,-1.6324168293282543629E+13,
!CS   3       1.1357022719979468624E+14, 1.0051899717115285432E+15/
!CS    DATA QJ1/1.1267125065029138050E+06, 6.4872502899596389593E+08,
!CS   1       2.7622777286244082666E+11, 8.4899346165481429307E+13,
!CS   2       1.7128800897135812012E+16, 1.7253905888447681194E+18,
!CS   3       1.3886978985861357615E+03/
     DATA PJ1/4.6179191852758252280D+00,-7.1329006872560947377D+03,&
&             4.5039658105749078904D+06,-1.4437717718363239107D+09,&
&             2.3569285397217157313D+11,-1.6324168293282543629D+13,&
&             1.1357022719979468624D+14, 1.0051899717115285432D+15/
     DATA QJ1/1.1267125065029138050D+06, 6.4872502899596389593D+08,&
&             2.7622777286244082666D+11, 8.4899346165481429307D+13,&
&             1.7128800897135812012D+16, 1.7253905888447681194D+18,&
&             1.3886978985861357615D+03/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
!CS    DATA PY0/2.2157953222280260820E+05,-5.9157479997408395984E+07,
!CS   1         7.2144548214502560419E+09,-3.7595974497819597599E+11,
!CS   2         5.4708611716525426053E+12, 4.0535726612579544093E+13,
!CS   3        -3.1714424660046133456E+02/
!CS    DATA QY0/8.2079908168393867438E+02, 3.8136470753052572164E+05,
!CS   1         1.2250435122182963220E+08, 2.7800352738690585613E+10,
!CS   2         4.1272286200406461981E+12, 3.0737873921079286084E+14/
      DATA PY0/2.2157953222280260820D+05,-5.9157479997408395984D+07,&
&              7.2144548214502560419D+09,-3.7595974497819597599D+11,&
&              5.4708611716525426053D+12, 4.0535726612579544093D+13,&
&             -3.1714424660046133456D+02/
      DATA QY0/8.2079908168393867438D+02, 3.8136470753052572164D+05,&
&              1.2250435122182963220D+08, 2.7800352738690585613D+10,&
&              4.1272286200406461981D+12, 3.0737873921079286084D+14/
!--------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
!        4.0  <  |X|  <=  8.0
!--------------------------------------------------------------------
!CS    DATA PY1/ 1.9153806858264202986E+06,-1.1957961912070617006E+09,
!CS   1          3.7453673962438488783E+11,-5.9530713129741981618E+13,
!CS   2          4.0686275289804744814E+15,-2.3638408497043134724E+16,
!CS   3         -5.6808094574724204577E+18, 1.1514276357909013326E+19,
!CS   4         -1.2337180442012953128E+03/
!CS    DATA QY1/ 1.2855164849321609336E+03, 1.0453748201934079734E+06,
!CS   1          6.3550318087088919566E+08, 3.0221766852960403645E+11,
!CS   2          1.1187010065856971027E+14, 3.0837179548112881950E+16,
!CS   3          5.6968198822857178911E+18, 5.3321844313316185697E+20/
       DATA PY1/ 1.9153806858264202986D+06,-1.1957961912070617006D+09,&
&                3.7453673962438488783D+11,-5.9530713129741981618D+13,&
&                4.0686275289804744814D+15,-2.3638408497043134724D+16,&
&               -5.6808094574724204577D+18, 1.1514276357909013326D+19,&
&               -1.2337180442012953128D+03/
      DATA QY1/ 1.2855164849321609336D+03, 1.0453748201934079734D+06,&
&               6.3550318087088919566D+08, 3.0221766852960403645D+11,&
&               1.1187010065856971027D+14, 3.0837179548112881950D+16,&
&               5.6968198822857178911D+18, 5.3321844313316185697D+20/
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
!CS    DATA P0/-1.0982405543459346727E+05,-1.5235293511811373833E+06,
!CS   1         -6.6033732483649391093E+06,-9.9422465050776411957E+06,
!CS   2         -4.4357578167941278571E+06,-1.6116166443246101165E+03/
!CS    DATA Q0/-1.0726385991103820119E+05,-1.5118095066341608816E+06,
!CS   1         -6.5853394797230870728E+06,-9.9341243899345856590E+06,
!CS   2         -4.4357578167941278568E+06,-1.4550094401904961825E+03/
!CS    DATA P1/ 1.7063754290207680021E+03, 1.8494262873223866797E+04,
!CS   1          6.6178836581270835179E+04, 8.5145160675335701966E+04,
!CS   2          3.3220913409857223519E+04, 3.5265133846636032186E+01/
!CS    DATA Q1/ 3.7890229745772202641E+04, 4.0029443582266975117E+05,
!CS   1          1.4194606696037208929E+06, 1.8194580422439972989E+06,
!CS   2          7.0871281941028743574E+05, 8.6383677696049909675E+02/
      DATA P0/-1.0982405543459346727D+05,-1.5235293511811373833D+06,&
&             -6.6033732483649391093D+06,-9.9422465050776411957D+06,&
&             -4.4357578167941278571D+06,-1.6116166443246101165D+03/
      DATA Q0/-1.0726385991103820119D+05,-1.5118095066341608816D+06,&
&             -6.5853394797230870728D+06,-9.9341243899345856590D+06,&
&             -4.4357578167941278568D+06,-1.4550094401904961825D+03/
      DATA P1/ 1.7063754290207680021D+03, 1.8494262873223866797D+04,&
&              6.6178836581270835179D+04, 8.5145160675335701966D+04,&
&              3.3220913409857223519D+04, 3.5265133846636032186D+01/
      DATA Q1/ 3.7890229745772202641D+04, 4.0029443582266975117D+05,&
&              1.4194606696037208929D+06, 1.8194580422439972989D+06,&
&              7.0871281941028743574D+05, 8.6383677696049909675D+02/
!-------------------------------------------------------------------
!  Check for error conditions
!-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. ((ARG .LE. ZERO) .OR.&
&        ((ARG .LT. HALF) .AND. (AX*XINF .LT. PI2)))) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) THEN
            GO TO 800
         ELSE IF (AX .LE. XSMALL) THEN
            IF (JINT .EQ. 0) THEN
                  RESULT = ARG * HALF
               ELSE
                  RESULT = -PI2 / AX
            END IF
            GO TO 2000
      END IF
!-------------------------------------------------------------------
!  Calculate J1 for appropriate interval, preserving
!     accuracy near the zero of J1
!-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(7) * ZSQ + PJ0(6)) * ZSQ + PJ0(5)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ARG * ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            XNUM = PJ1(1)
            XDEN = (ZSQ + QJ1(7)) * ZSQ + QJ1(1)
            DO 220 I = 2, 6
               XNUM = XNUM * ZSQ + PJ1(I)
               XDEN = XDEN * ZSQ + QJ1(I)
  220       CONTINUE
            XNUM = XNUM * (AX - EIGHT) * (AX + EIGHT) + PJ1(7)
            XNUM = XNUM * (AX - FOUR) * (AX + FOUR) + PJ1(8)
            PROD = ARG * ((AX - XJ11/TWO56) - XJ12) * (AX + XJ1)
      END IF
      RESULT = PROD * (XNUM / XDEN)
      IF (JINT .EQ. 0) GO TO 2000
!-------------------------------------------------------------------
!  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
!    where xn is a zero of Y1
!-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
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
!-------------------------------------------------------------------
!  Now calculate Y1 for appropriate interval, preserving
!     accuracy near the zero of Y1
!-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            XNUM = PY0(7) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 6
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE
            XNUM = PY1(9) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 8
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
      END IF
      RESULT = RESJ + (UP*DOWN/AX) * XNUM / XDEN
      GO TO 2000
!-------------------------------------------------------------------
!  Calculate J1 or Y1 for |ARG|  >  8.0
!-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AINT(AX/TWOPI) + THROV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(6)
      XDEN = ZSQ + Q0(6)
      UP = P1(6)
      DOWN = ZSQ + Q1(6)
      DO 850 I = 1, 5
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = (RTPI2/SQRT(AX)) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = (RTPI2/SQRT(AX)) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
      IF ((JINT .EQ. 0) .AND. (ARG .LT. ZERO)) RESULT = -RESULT
 2000 RETURN
!---------- Last card of CALJY1 ----------
      END subroutine caljy1

      DOUBLE PRECISION FUNCTION BESJ1(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for Bessel functions
!   of the first kind of order zero for arguments  |X| <= XMAX
!   (see comments heading CALJY1).
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: JINT
!CS    REAL
     DOUBLE PRECISION :: &
&        RESULT,X
!--------------------------------------------------------------------
      JINT=0
      CALL CALJY1(X,RESULT,JINT)
      BESJ1 = RESULT
      RETURN
!---------- Last card of BESJ1 ----------
      END function besj1


      DOUBLE PRECISION FUNCTION BESY1(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for Bessel functions
!   of the second kind of order zero for arguments 0 < X <= XMAX
!   (see comments heading CALJY1).
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: JINT
!CS    REAL
     DOUBLE PRECISION :: &
&        RESULT,X
!--------------------------------------------------------------------
      JINT=1
      CALL CALJY1(X,RESULT,JINT)
      BESY1 = RESULT
      RETURN
!---------- Last card of BESY1 ----------
      END function besy1
!!***

!!****f* ABINIT/CALCK0
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

SUBROUTINE CALCK0(ARG,RESULT,JINT)

!--------------------------------------------------------------------
!
! This packet computes modified Bessel functions of the second kind
!   and order zero, K0(X) and EXP(X)*K0(X), for real
!   arguments X.  It contains two function type subprograms, BESK0
!   and BESEK0, and one subroutine type subprogram, CALCK0.
!   the calling statements for the primary entries are
!
!                   Y=BESK0(X)
!   and
!                   Y=BESEK0(X)
!
!   where the entry points correspond to the functions K0(X) and
!   EXP(X)*K0(X), respectively.  The routine CALCK0 is
!   intended for internal packet use only, all computations within
!   the packet being concentrated in this routine.  The function
!   subprograms invoke CALCK0 with the statement
!          CALL CALCK0(ARG,RESULT,JINT)
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALCK0
!       Call              ARG                  RESULT          JINT
!
!     BESK0(ARG)   0 .LT. ARG .LE. XMAX       K0(ARG)           1
!     BESEK0(ARG)     0 .LT. ARG           EXP(ARG)*K0(ARG)     2
!
!   The main computation evaluates slightly modified forms of near
!   minimax rational approximations generated by Russon and Blair,
!   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!   1969.  This transportable program is patterned after the
!   machine-dependent FUNPACK packet NATSK0, but cannot match that
!   version for efficiency or accuracy.  This version uses rational
!   functions that theoretically approximate K-SUB-0(X) to at
!   least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!   XSMALL = Argument below which BESK0 and BESEK0 may
!            each be represented by a constant and a log.
!            largest X such that  1.0 + X = 1.0  to machine
!            precision.
!   XINF   = Largest positive machine number; approximately
!            beta**maxexp
!   XMAX   = Largest argument acceptable to BESK0;  Solution to
!            equation:
!               W(X) * (1-1/8X+9/128X**2) = beta**minexp
!            where  W(X) = EXP(-X)*SQRT(PI/2X)
!
!
!     Approximate values for some important machines are:
!
!
!                           beta       minexp       maxexp
!
!  CRAY-1        (S.P.)       2        -8193         8191
!  Cyber 180/185
!    under NOS   (S.P.)       2         -975         1070
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)       2         -126          128
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)       2        -1022         1024
!  IBM 3033      (D.P.)      16          -65           63
!  VAX D-Format  (D.P.)       2         -128          127
!  VAX G-Format  (D.P.)       2        -1024         1023
!
!
!                          XSMALL       XINF         XMAX
!
! CRAY-1        (S.P.)    3.55E-15   5.45E+2465    5674.858
! Cyber 180/855
!   under NOS   (S.P.)    1.77E-15   1.26E+322      672.788
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)    5.95E-8    3.40E+38        85.337
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)    1.11D-16   1.79D+308      705.342
! IBM 3033      (D.P.)    1.11D-16   7.23D+75       177.852
! VAX D-Format  (D.P.)    6.95D-18   1.70D+38        86.715
! VAX G-Format  (D.P.)    5.55D-17   8.98D+307      706.728
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns the value XINF for ARG .LE. 0.0, and the
!  BESK0 entry returns the value 0.0 for ARG .GT. XMAX.
!
!
!  Intrinsic functions required are:
!
!     EXP, LOG, SQRT
!
!  Latest modification: March 19, 1990
!
!  Authors: W. J. Cody and Laura Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!  Original subroutine from netlib http://www.netlib.org/specfun/k0
!  Slightly modified by MG to follow f90 rules and double precision arithmetic
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,JINT
!CS    REAL
      DOUBLE PRECISION :: ARG,RESULT,SUMF,SUMG,SUMP,SUMQ,TEMP
!CS    REAL
      DOUBLE PRECISION :: X,XX
!CS    REAL
      DOUBLE PRECISION :: P(6),Q(2),PP(10),QQ(10),F(4),G(3)
!C--------------------------------------------------------------------
!C  Mathematical constants
!C--------------------------------------------------------------------
!CS    REAL, PARAMETER  ::             ONE=1.0E0,ZERO=0.0E0
      DOUBLE PRECISION,PARAMETER ::  ONE=1.0D0,ZERO=0.0D0
!C--------------------------------------------------------------------
!C  Machine-dependent constants
!C--------------------------------------------------------------------
!CS    REAL.PARAMETER ::             XSMALL=5.95E-8, XINF=3.40E+38 ,XMAX=85.337E0
      DOUBLE PRECISION,PARAMETER :: XSMALL=1.11D-16,XINF=1.79D+308,XMAX=705.342D0
!--------------------------------------------------------------------
!
!     Coefficients for XSMALL .LE.  ARG  .LE. 1.0
!
!--------------------------------------------------------------------
!S    DATA   P/ 5.8599221412826100000E-04, 1.3166052564989571850E-01,
!S   1          1.1999463724910714109E+01, 4.6850901201934832188E+02,
!S   2          5.9169059852270512312E+03, 2.4708152720399552679E+03/
!S    DATA   Q/-2.4994418972832303646E+02, 2.1312714303849120380E+04/
!S    DATA   F/-1.6414452837299064100E+00,-2.9601657892958843866E+02,
!S   1         -1.7733784684952985886E+04,-4.0320340761145482298E+05/
!S    DATA   G/-2.5064972445877992730E+02, 2.9865713163054025489E+04,
!S   1         -1.6128136304458193998E+06/
     DATA    P/5.8599221412826100000D-04,1.3166052564989571850D-01,&
&              1.1999463724910714109D+01,4.6850901201934832188D+02,&
&              5.9169059852270512312D+03,2.4708152720399552679D+03/
     DATA    Q/-2.4994418972832303646D+02, 2.1312714303849120380D+04/
     DATA    F/-1.6414452837299064100D+00,-2.9601657892958843866D+02,&
&              -1.7733784684952985886D+04,-4.0320340761145482298D+05/
     DATA    G/-2.5064972445877992730D+02, 2.9865713163054025489D+04,&
&              -1.6128136304458193998D+06/
!--------------------------------------------------------------------
!
!     Coefficients for  1.0 .LT. ARG
!
!--------------------------------------------------------------------
!S    DATA  PP/ 1.1394980557384778174E+02, 3.6832589957340267940E+03,
!S   1          3.1075408980684392399E+04, 1.0577068948034021957E+05,
!S   2          1.7398867902565686251E+05, 1.5097646353289914539E+05,
!S   3          7.1557062783764037541E+04, 1.8321525870183537725E+04,
!S   4          2.3444738764199315021E+03, 1.1600249425076035558E+02/
!S    DATA  QQ/ 2.0013443064949242491E+02, 4.4329628889746408858E+03,
!S   1          3.1474655750295278825E+04, 9.7418829762268075784E+04,
!S   2          1.5144644673520157801E+05, 1.2689839587977598727E+05,
!S   3          5.8824616785857027752E+04, 1.4847228371802360957E+04,
!S   4          1.8821890840982713696E+03, 9.2556599177304839811E+01/
     DATA  PP/  1.1394980557384778174D+02, 3.6832589957340267940D+03,&
&               3.1075408980684392399D+04, 1.0577068948034021957D+05,&
&               1.7398867902565686251D+05, 1.5097646353289914539D+05,&
&               7.1557062783764037541D+04, 1.8321525870183537725D+04,&
&               2.3444738764199315021D+03, 1.1600249425076035558D+02/
     DATA  QQ/  2.0013443064949242491D+02, 4.4329628889746408858D+03, &
&                3.1474655750295278825D+04, 9.7418829762268075784D+04,&
&                1.5144644673520157801D+05, 1.2689839587977598727D+05,&
&                5.8824616785857027752D+04, 1.4847228371802360957D+04,&
&                1.8821890840982713696D+03, 9.2556599177304839811D+01/
!--------------------------------------------------------------------
      X = ARG
      IF (X .GT. ZERO) THEN
            IF (X .LE. ONE) THEN
!--------------------------------------------------------------------
!     0.0 .LT.  ARG  .LE. 1.0
!--------------------------------------------------------------------
                  TEMP = LOG(X)
                  IF (X .LT. XSMALL) THEN
!--------------------------------------------------------------------
!     Return for small ARG
!--------------------------------------------------------------------
                        RESULT = P(6)/Q(2) - TEMP
                     ELSE
                        XX = X * X
                        SUMP = ((((P(1)*XX + P(2))*XX + P(3))*XX +&
                               P(4))*XX + P(5))*XX + P(6)
                        SUMQ = (XX + Q(1))*XX + Q(2)
                        SUMF = ((F(1)*XX + F(2))*XX + F(3))*XX + F(4)
                        SUMG = ((XX + G(1))*XX + G(2))*XX + G(3)
                        RESULT = SUMP/SUMQ - XX*SUMF*TEMP/SUMG - TEMP
                        IF (JINT .EQ. 2) RESULT = RESULT * EXP(X)
                  END IF
               ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
!--------------------------------------------------------------------
!     Error return for ARG .GT. XMAX
!--------------------------------------------------------------------
                  RESULT = ZERO
               ELSE
!--------------------------------------------------------------------
!     1.0 .LT. ARG
!--------------------------------------------------------------------
                  XX = ONE / X
                  SUMP = PP(1)
                  DO 120 I = 2, 10
                     SUMP = SUMP*XX + PP(I)
  120             CONTINUE
                  SUMQ = XX
                  DO 140 I = 1, 9
                     SUMQ = (SUMQ + QQ(I))*XX
  140             CONTINUE
                  SUMQ = SUMQ + QQ(10)
                  RESULT = SUMP / SUMQ / SQRT(X)
                  IF (JINT .EQ. 1) RESULT = RESULT * EXP(-X)
            END IF
         ELSE
!--------------------------------------------------------------------
!     Error return for ARG .LE. 0.0
!--------------------------------------------------------------------
            RESULT = XINF
      END IF
!--------------------------------------------------------------------
!     Update error counts, etc.
!--------------------------------------------------------------------
      RETURN
!---------- Last line of CALCK0 ----------
      END subroutine calck0
!!***

!S    REAL
      DOUBLE PRECISION FUNCTION BESK0(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   modified Bessel function of the second kind of order zero
!   for arguments 0.0 .LT. ARG .LE. XMAX (see comments heading
!   CALCK0).
!
!  Authors: W. J. Cody and Laura Stoltz
!
!  Latest Modification: January 19, 1988
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: JINT
!S    REAL
      DOUBLE PRECISION :: X, RESULT
!--------------------------------------------------------------------
      JINT = 1
      CALL CALCK0(X,RESULT,JINT)
      BESK0 = RESULT
      RETURN
!---------- Last line of BESK0 ----------
      END function besk0
!!***

!S    REAL
      DOUBLE PRECISION FUNCTION BESEK0(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   modified Bessel function of the second kind of order zero
!   multiplied by the Exponential function, for arguments
!   0.0 .LT. ARG.
!
!  Authors: W. J. Cody and Laura Stoltz
!
!  Latest Modification: January 19, 1988
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER JINT
!S    REAL
      DOUBLE PRECISION :: X,RESULT
!--------------------------------------------------------------------
      JINT = 2
      CALL CALCK0(X,RESULT,JINT)
      BESEK0 = RESULT
      RETURN
!---------- Last line of BESEK0 ----------
      END function BESEK0
!!***

!!****f* ABINIT/CALCK1
!! NAME
!!  CALCK1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

SUBROUTINE CALCK1(ARG,RESULT,JINT)

!--------------------------------------------------------------------
!
! This packet computes modified Bessel functions of the second kind
!   and order one,  K1(X)  and  EXP(X)*K1(X), for real arguments X.
!   It contains two function type subprograms, BESK1  and  BESEK1,
!   and one subroutine type subprogram, CALCK1.  The calling
!   statements for the primary entries are
!
!                   Y=BESK1(X)
!   and
!                   Y=BESEK1(X)
!
!   where the entry points correspond to the functions K1(X) and
!   EXP(X)*K1(X), respectively.  The routine CALCK1 is intended
!   for internal packet use only, all computations within the
!   packet being concentrated in this routine.  The function
!   subprograms invoke CALCK1 with the statement
!          CALL CALCK1(ARG,RESULT,JINT)
!   where the parameter usage is as follows
!
!      Function                      Parameters for CALCK1
!        Call             ARG                  RESULT          JINT
!
!     BESK1(ARG)  XLEAST .LT. ARG .LT. XMAX    K1(ARG)          1
!     BESEK1(ARG)     XLEAST .LT. ARG       EXP(ARG)*K1(ARG)    2
!
!   The main computation evaluates slightly modified forms of near
!   minimax rational approximations generated by Russon and Blair,
!   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!   1969.  This transportable program is patterned after the
!   machine-dependent FUNPACK packet NATSK1, but cannot match that
!   version for efficiency or accuracy.  This version uses rational
!   functions that theoretically approximate K-SUB-1(X) to at
!   least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!   XLEAST = Smallest acceptable argument, i.e., smallest machine
!            number X such that 1/X is machine representable.
!   XSMALL = Argument below which BESK1(X) and BESEK1(X) may
!            each be represented by 1/X.  A safe value is the
!            largest X such that  1.0 + X = 1.0  to machine
!            precision.
!   XINF   = Largest positive machine number; approximately
!            beta**maxexp
!   XMAX   = Largest argument acceptable to BESK1;  Solution to
!            equation:
!               W(X) * (1+3/8X-15/128X**2) = beta**minexp
!            where  W(X) = EXP(-X)*SQRT(PI/2X)
!
!
!     Approximate values for some important machines are:
!
!                           beta       minexp       maxexp
!
!  CRAY-1        (S.P.)       2        -8193         8191
!  Cyber 180/185
!    under NOS   (S.P.)       2         -975         1070
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)       2         -126          128
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)       2        -1022         1024
!  IBM 3033      (D.P.)      16          -65           63
!  VAX D-Format  (D.P.)       2         -128          127
!  VAX G-Format  (D.P.)       2        -1024         1023
!
!
!                         XLEAST     XSMALL      XINF       XMAX
!
! CRAY-1                1.84E-2466  3.55E-15  5.45E+2465  5674.858
! Cyber 180/855
!   under NOS   (S.P.)  3.14E-294   1.77E-15  1.26E+322    672.789
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)  1.18E-38    5.95E-8   3.40E+38      85.343
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)  2.23D-308   1.11D-16  1.79D+308    705.343
! IBM 3033      (D.P.)  1.39D-76    1.11D-16  7.23D+75     177.855
! VAX D-Format  (D.P.)  5.88D-39    6.95D-18  1.70D+38      86.721
! VAX G-Format  (D.P.)  1.12D-308   5.55D-17  8.98D+307    706.728
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns the value XINF for ARG .LE. 0.0 and the
!   BESK1 entry returns the value 0.0 for ARG .GT. XMAX.
!
!
!  Intrinsic functions required are:
!
!     LOG, SQRT, EXP
!
!
!  Authors: W. J. Cody and Laura Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!  Latest modification: January 28, 1988
!  Taken from http://www.netlib.org/specfun/k1
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,JINT
!CS    REAL
      DOUBLE PRECISION :: &
&         ARG,F,G,ONE,P,PP,Q,QQ,RESULT,SUMF,SUMG,&
&         SUMP,SUMQ,X,XINF,XMAX,XLEAST,XSMALL,XX,ZERO
      DIMENSION P(5),Q(3),PP(11),QQ(9),F(5),G(3)
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
!CS    DATA ONE/1.0E0/,ZERO/0.0E0/
       DATA ONE/1.0D0/,ZERO/0.0D0/
!--------------------------------------------------------------------
!  Machine-dependent constants
!--------------------------------------------------------------------
!CS    DATA XLEAST/1.18E-38/,XSMALL/5.95E-8/,XINF/3.40E+38/,
!CS   1     XMAX/85.343E+0/
     DATA XLEAST/2.23D-308/,XSMALL/1.11D-16/,XINF/1.79D+308/,&
&         XMAX/705.343D+0/
!--------------------------------------------------------------------
!  Coefficients for  XLEAST .LE.  ARG  .LE. 1.0
!--------------------------------------------------------------------
!CS    DATA   P/ 4.8127070456878442310E-1, 9.9991373567429309922E+1,
!CS   1          7.1885382604084798576E+3, 1.7733324035147015630E+5,
!CS   2          7.1938920065420586101E+5/
!CS    DATA   Q/-2.8143915754538725829E+2, 3.7264298672067697862E+4,
!CS   1         -2.2149374878243304548E+6/
!CS    DATA   F/-2.2795590826955002390E-1,-5.3103913335180275253E+1,
!CS   1         -4.5051623763436087023E+3,-1.4758069205414222471E+5,
!CS   2         -1.3531161492785421328E+6/
!CS    DATA   G/-3.0507151578787595807E+2, 4.3117653211351080007E+4,
!CS   2         -2.7062322985570842656E+6/
      DATA   P/ 4.8127070456878442310D-1, 9.9991373567429309922D+1,&
&             7.1885382604084798576D+3, 1.7733324035147015630D+5,&
&             7.1938920065420586101D+5/
    DATA   Q/-2.8143915754538725829D+2, 3.7264298672067697862D+4,&
&            -2.2149374878243304548D+6/
    DATA   F/-2.2795590826955002390D-1,-5.3103913335180275253D+1,&
&            -4.5051623763436087023D+3,-1.4758069205414222471D+5,&
&            -1.3531161492785421328D+6/
    DATA   G/-3.0507151578787595807D+2, 4.3117653211351080007D+4,&
&            -2.7062322985570842656D+6/
!--------------------------------------------------------------------
!  Coefficients for  1.0 .LT.  ARG
!--------------------------------------------------------------------
!CS    DATA  PP/ 6.4257745859173138767E-2, 7.5584584631176030810E+0,
!CS   1          1.3182609918569941308E+2, 8.1094256146537402173E+2,
!CS   2          2.3123742209168871550E+3, 3.4540675585544584407E+3,
!CS   3          2.8590657697910288226E+3, 1.3319486433183221990E+3,
!CS   4          3.4122953486801312910E+2, 4.4137176114230414036E+1,
!CS   5          2.2196792496874548962E+0/
!CS    DATA  QQ/ 3.6001069306861518855E+1, 3.3031020088765390854E+2,
!CS   1          1.2082692316002348638E+3, 2.1181000487171943810E+3,
!CS   2          1.9448440788918006154E+3, 9.6929165726802648634E+2,
!CS   3          2.5951223655579051357E+2, 3.4552228452758912848E+1,
!CS   4          1.7710478032601086579E+0/
    DATA  PP/ 6.4257745859173138767D-2, 7.5584584631176030810D+0,&
&             1.3182609918569941308D+2, 8.1094256146537402173D+2,&
&             2.3123742209168871550D+3, 3.4540675585544584407D+3,&
&             2.8590657697910288226D+3, 1.3319486433183221990D+3,&
&             3.4122953486801312910D+2, 4.4137176114230414036D+1,&
&             2.2196792496874548962D+0/
    DATA  QQ/ 3.6001069306861518855D+1, 3.3031020088765390854D+2,&
&             1.2082692316002348638D+3, 2.1181000487171943810D+3,&
&             1.9448440788918006154D+3, 9.6929165726802648634D+2,&
&             2.5951223655579051357D+2, 3.4552228452758912848D+1,&
&             1.7710478032601086579D+0/
!--------------------------------------------------------------------
      X = ARG
      IF (X .LT. XLEAST) THEN
!--------------------------------------------------------------------
!  Error return for  ARG  .LT. XLEAST
!--------------------------------------------------------------------
            RESULT = XINF
         ELSE IF (X .LE. ONE) THEN
!--------------------------------------------------------------------
!  XLEAST .LE.  ARG  .LE. 1.0
!--------------------------------------------------------------------
            IF (X .LT. XSMALL) THEN
!--------------------------------------------------------------------
!  Return for small ARG
!--------------------------------------------------------------------
                  RESULT = ONE / X
               ELSE
                  XX = X * X
                  SUMP = ((((P(1)*XX + P(2))*XX + P(3))*XX + P(4))*XX &
&                        + P(5))*XX + Q(3)
                  SUMQ = ((XX + Q(1))*XX + Q(2))*XX + Q(3)
                  SUMF = (((F(1)*XX + F(2))*XX + F(3))*XX + F(4))*XX &
&                        + F(5)
                  SUMG = ((XX + G(1))*XX + G(2))*XX + G(3)
                  RESULT = (XX * LOG(X) * SUMF/SUMG + SUMP/SUMQ) / X
                  IF (JINT .EQ. 2) RESULT = RESULT * EXP(X)
            END IF
         ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
!--------------------------------------------------------------------
!  Error return for  ARG  .GT. XMAX
!--------------------------------------------------------------------
            RESULT = ZERO
         ELSE
!--------------------------------------------------------------------
!  1.0 .LT.  ARG
!--------------------------------------------------------------------
            XX = ONE / X
            SUMP = PP(1)
            DO 120 I = 2, 11
               SUMP = SUMP * XX + PP(I)
  120       CONTINUE
            SUMQ = XX
            DO 140 I = 1, 8
               SUMQ = (SUMQ + QQ(I)) * XX
  140       CONTINUE
            SUMQ = SUMQ + QQ(9)
            RESULT = SUMP / SUMQ / SQRT(X)
            IF (JINT .EQ. 1) RESULT = RESULT * EXP(-X)
      END IF
      RETURN
!---------- Last line of CALCK1 ----------
      END subroutine calck1
!!***

!CS    REAL
      DOUBLE PRECISION FUNCTION BESK1(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   modified Bessel function of the second kind of order one
!   for arguments  XLEAST .LE. ARG .LE. XMAX.
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: JINT
!CS    REAL
      DOUBLE PRECISION :: &
&         X, RESULT
!--------------------------------------------------------------------
      JINT = 1
      CALL CALCK1(X,RESULT,JINT)
      BESK1 = RESULT
      RETURN
!---------- Last line of BESK1 ----------
      END function besk1
!!***

!CS    REAL
      DOUBLE PRECISION  FUNCTION BESEK1(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   modified Bessel function of the second kind of order one
!   multiplied by the exponential function, for arguments
!   XLEAST .LE. ARG .LE. XMAX.
!
!--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER JINT
!CS    REAL
      DOUBLE PRECISION :: &
&         X, RESULT
!--------------------------------------------------------------------
      JINT = 2
      CALL CALCK1(X,RESULT,JINT)
      BESEK1 = RESULT
      RETURN
!---------- Last line of BESEK1 ----------
      END function besek1
!!***

end module m_bessel
!!***

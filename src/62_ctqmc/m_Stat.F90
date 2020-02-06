
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_Stat
!! NAME
!!  m_Stat
!! 
!! FUNCTION 
!!  FIXME: add description. 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#include "defs.h"
MODULE m_Stat
USE m_global
IMPLICIT NONE
  
PRIVATE

PUBLIC :: Stat_average
PUBLIC :: Stat_variance
PUBLIC :: Stat_coVariance
PUBLIC :: Stat_deviation
PUBLIC :: Stat_linearReg
PUBLIC :: Stat_powerReg

CONTAINS
!!***

DOUBLE PRECISION FUNCTION Stat_average(tab)

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab
  !INTEGER :: sizet

  !sizet = SIZE(tab)
  Stat_average = 0.d0
  !IF ( sizet .GT. 0 ) &
  Stat_average = SUM(tab)/DBLE(SIZE(tab))
END FUNCTION Stat_average
!!***

DOUBLE PRECISION FUNCTION Stat_variance(tab)

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab
  INTEGER                                    :: sizet
  INTEGER                                    :: i
  DOUBLE PRECISION                           :: average
  DOUBLE PRECISION                           :: tab2
  
  sizet = SIZE(tab)
  Stat_variance = 0.d0

  !IF ( sizet .GT. 0 ) THEN
    average = Stat_average(tab)
    DO i = 1, sizet
      tab2 = tab(i)-average
      tab2 = tab2 * tab2
      Stat_variance = Stat_variance + tab2
    END DO
    Stat_variance = Stat_variance / DBLE(sizet)
  !ELSE
  !  Stat_variance = 0.d-16
  !END IF
END FUNCTION Stat_variance
!!***

DOUBLE PRECISION FUNCTION Stat_coVariance(tab1, tab2)

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab1
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab2
  INTEGER                                    :: size1
  INTEGER                                    :: size2
  INTEGER                                    :: i
  DOUBLE PRECISION                           :: average1
  DOUBLE PRECISION                           :: average2
  DOUBLE PRECISION                           :: tmp

  size1    = SIZE(tab1)
  size2    = SIZE(tab2)

  IF ( size1 .NE. size2 ) &
    CALL ERROR("Stat_coVariance : Array sizes mismatch            ")

  average1 = Stat_average(tab1)
  average2 = Stat_average(tab2)
  Stat_coVariance = 0.d0

  !IF ( size1 .GT. 0 ) THEN
    DO i = 1, size1
      tmp = (tab1(i)-average1)*(tab2(i)-average2)
      Stat_coVariance = Stat_coVariance + tmp
    END DO
    Stat_coVariance = Stat_coVariance / DBLE(size1)
  !ELSE
  !  Stat_coVariance = 0.d-16
  !END IF
END FUNCTION Stat_coVariance
!!***

DOUBLE PRECISION FUNCTION Stat_deviation(tab1)

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab1

  Stat_deviation = SQRT(Stat_variance(tab1))
END FUNCTION Stat_deviation
!!***

!!****f* ABINIT/m_Stat/Stat_linearReg
!! NAME
!!  Stat_linearReg
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE Stat_linearReg(tabX, tabY, a, b, R)
!Arguments ------------------------------------
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN ) :: tabX
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN ) :: tabY
  DOUBLE PRECISION              , INTENT(OUT) :: a
  DOUBLE PRECISION              , INTENT(OUT) :: b
  DOUBLE PRECISION              , INTENT(OUT) :: R
  DOUBLE PRECISION                            :: coVar
  DOUBLE PRECISION                            :: Var

  coVar = Stat_coVariance(tabX, tabY)
  Var   = Stat_variance(tabX)
  a = coVar / var
  b = Stat_average(tabY) - a* Stat_average(tabX)
  R = ABS(a * SQRT(var) / Stat_deviation(tabY))
END SUBROUTINE Stat_linearReg
!!***
  
!!****f* ABINIT/m_Stat/Stat_powerReg
!! NAME
!!  Stat_powerReg
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE Stat_powerReg(tabX, tabY, a, b, R)
!Arguments ------------------------------------
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN ) :: tabX
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN ) :: tabY
  DOUBLE PRECISION              , INTENT(OUT) :: a
  DOUBLE PRECISION              , INTENT(OUT) :: b
  DOUBLE PRECISION              , INTENT(OUT) :: R
  INTEGER                                     :: size1
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tab1
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tab2
  INTEGER :: i

  size1 = SIZE(tabX)
  IF ( size1 .NE. SIZE(tabY) ) &
    CALL ERROR("Stat_powerReg : Array sizes mismatch              ")

  FREEIF(tab1)
  FREEIF(tab2)
  MALLOC(tab1,(1:size1))
  MALLOC(tab2,(1:size1))


!  IF ( ISNAN(a) .OR. ISNAN(b) ) THEN
    DO i = 1, size1 
      tab1(i) = LOG(tabX(i))
      IF ( tabY(i) .NE. 0.d0 ) THEN
        tab2(i) = LOG(tabY(i))
      ELSE
        tab2(i) = LOG(1.d-16)
      END IF
      !WRITE(92,'(4E22.14)') tabX(i), tab1(i), tabY(i), tab2(i)
    END DO
    !WRITE(92,*)
  !END IF
  CALL Stat_linearReg(tab1, tab2, b, a, R)
  a = EXP(a)
  FREE(tab1)
  FREE(tab2)
END SUBROUTINE Stat_powerReg
!!***

DOUBLE PRECISION FUNCTION Stat_simpson(func, a, b, N)

  DOUBLE PRECISION, DIMENSION(:), POINTER :: func     !vz_i
  DOUBLE PRECISION                       , INTENT(IN) :: a
  DOUBLE PRECISION                       , INTENT(IN) :: b
  INTEGER                                , INTENT(IN) :: N
  INTEGER :: i
  INTEGER :: x0
  INTEGER :: x1
  INTEGER :: x2 
  DOUBLE PRECISION :: dtau
  DOUBLE PRECISION :: J

  dtau = (b-a)/DBLE(N-1)
  J=0.d0
  x2 = 1
  x1 = 0
  DO i = 1, N-2, 2
    x0 = x2
    x1 = x1 + 2    
    x2 = x0 + 2   
    J = J + func(x0) + 4.d0*func(x1) + func(x2)
  END DO
  J = J * dtau / 3.d0
  Stat_simpson = J
END FUNCTION Stat_simpson
!!***

END  MODULE m_Stat

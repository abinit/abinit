#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_Statoffdiag
!! NAME
!!  m_Statoffdiag
!! 
!! FUNCTION 
!!  FIXME: add description. 
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
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
MODULE m_Statoffdiag
USE m_global
IMPLICIT NONE
  
CONTAINS
!!***

DOUBLE PRECISION FUNCTION Statoffdiag_average(tab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_average'
!End of the abilint section

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab
  !INTEGER :: sizet

  !sizet = SIZE(tab)
  Statoffdiag_average = 0.d0
  !IF ( sizet .GT. 0 ) &
  Statoffdiag_average = SUM(tab)/DBLE(SIZE(tab))
END FUNCTION Statoffdiag_average
!!***

DOUBLE PRECISION FUNCTION Statoffdiag_variance(tab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_variance'
!End of the abilint section

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab
  INTEGER                                    :: sizet
  INTEGER                                    :: i
  DOUBLE PRECISION                           :: average
  DOUBLE PRECISION                           :: tab2
  
  sizet = SIZE(tab)
  Statoffdiag_variance = 0.d0

  !IF ( sizet .GT. 0 ) THEN
    average = Statoffdiag_average(tab)
    DO i = 1, sizet
      tab2 = tab(i)-average
      tab2 = tab2 * tab2
      Statoffdiag_variance = Statoffdiag_variance + tab2
    END DO
    Statoffdiag_variance = Statoffdiag_variance / DBLE(sizet)
  !ELSE
  !  Statoffdiag_variance = 0.d-16
  !END IF
END FUNCTION Statoffdiag_variance
!!***

DOUBLE PRECISION FUNCTION Statoffdiag_coVariance(tab1, tab2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_coVariance'
!End of the abilint section

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
    CALL ERROR("Statoffdiag_coVariance : Array sizes mismatch            ")

  average1 = Statoffdiag_average(tab1)
  average2 = Statoffdiag_average(tab2)
  Statoffdiag_coVariance = 0.d0

  !IF ( size1 .GT. 0 ) THEN
    DO i = 1, size1
      tmp = (tab1(i)-average1)*(tab2(i)-average2)
      Statoffdiag_coVariance = Statoffdiag_coVariance + tmp
    END DO
    Statoffdiag_coVariance = Statoffdiag_coVariance / DBLE(size1)
  !ELSE
  !  Statoffdiag_coVariance = 0.d-16
  !END IF
END FUNCTION Statoffdiag_coVariance
!!***

DOUBLE PRECISION FUNCTION Statoffdiag_deviation(tab1)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_deviation'
!End of the abilint section

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: tab1

  Statoffdiag_deviation = SQRT(Statoffdiag_variance(tab1))
END FUNCTION Statoffdiag_deviation
!!***

!!****f* ABINIT/m_Statoffdiag/Statoffdiag_linearReg
!! NAME
!!  Statoffdiag_linearReg
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
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

SUBROUTINE Statoffdiag_linearReg(tabX, tabY, a, b, R)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_linearReg'
!End of the abilint section

  DOUBLE PRECISION, DIMENSION(:), INTENT(IN ) :: tabX
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN ) :: tabY
  DOUBLE PRECISION              , INTENT(OUT) :: a
  DOUBLE PRECISION              , INTENT(OUT) :: b
  DOUBLE PRECISION              , INTENT(OUT) :: R
  DOUBLE PRECISION                            :: coVar
  DOUBLE PRECISION                            :: Var

  coVar = Statoffdiag_coVariance(tabX, tabY)
  Var   = Statoffdiag_variance(tabX)
  a = coVar / var
  b = Statoffdiag_average(tabY) - a* Statoffdiag_average(tabX)
  R = ABS(a * SQRT(var) / Statoffdiag_deviation(tabY))
END SUBROUTINE Statoffdiag_linearReg
!!***
  
!!****f* ABINIT/m_Statoffdiag/Statoffdiag_powerReg
!! NAME
!!  Statoffdiag_powerReg
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
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

SUBROUTINE Statoffdiag_powerReg(tabX, tabY, a, b, R)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_powerReg'
!End of the abilint section

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
    CALL ERROR("Statoffdiag_powerReg : Array sizes mismatch              ")

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
  CALL Statoffdiag_linearReg(tab1, tab2, b, a, R)
  a = EXP(a)
  FREE(tab1)
  FREE(tab2)
END SUBROUTINE Statoffdiag_powerReg
!!***

DOUBLE PRECISION FUNCTION Statoffdiag_simpson(func, a, b, N)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Statoffdiag_simpson'
!End of the abilint section

  DOUBLE PRECISION, DIMENSION(:), POINTER, INTENT(IN) :: func
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
  Statoffdiag_simpson = J
END FUNCTION Statoffdiag_simpson
!!***

END  MODULE m_Statoffdiag

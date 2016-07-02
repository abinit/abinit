 
! Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

      real(kind=kind(0.0d0)) FUNCTION SURPLA( C )
 
!  CALCULATES THE SRFACE OF THE UNIT CELL NORMAL TO THE INTERFACE

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'SURPLA'
!End of the abilint section

      IMPLICIT NONE

      real(kind=kind(0.0d0)) C(3,3)
      SURPLA = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) **2 +&
     &         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) **2 +&
     &         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) **2
      SURPLA = SQRT( ABS( SURPLA ) )
      END FUNCTION SURPLA

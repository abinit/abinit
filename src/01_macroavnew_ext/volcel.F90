#if defined HAVE_CONFIG_H
#include "config.h"
#endif

      DOUBLE PRECISION FUNCTION VOLCEL( C )
 
! Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .
 
!  CALCULATES THE VOLUME OF THE UNIT CELL

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'VOLCEL'
!End of the abilint section

      IMPLICIT NONE

      DOUBLE PRECISION C(3,3)
      VOLCEL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) +&
     &         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) +&
     &         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOLCEL = ABS( VOLCEL )
      END FUNCTION VOLCEL

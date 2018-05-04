#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_Vectoroffdiag
!! NAME
!!  m_Vectoroffdiag
!! 
!! FUNCTION 
!!  Manage a double precision vector
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
MODULE m_Vectoroffdiag
USE m_Global
IMPLICIT NONE

!!***

!!****t* m_Vectoroffdiag/Vectoroffdiag
!! NAME
!!  Vectoroffdiag
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE Vectoroffdiag
  INTEGER :: size
  INTEGER :: tail
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vec 
END TYPE Vectoroffdiag
!!***

CONTAINS
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_init
!! NAME
!!  Vectoroffdiag_init
!!
!! FUNCTION
!!  initialize
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=ListCdagCoffdiag
!!  size=size of initialization
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_init(vector_1, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_init'
!End of the abilint section

  TYPE(Vectoroffdiag)     , INTENT(INOUT) :: vector_1
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  vector_1%size = size_val
  FREEIF(vector_1%vec)
  MALLOC(vector_1%vec,(1:size_val))
  vector_1%tail     = 0 
  vector_1%vec = 0.d0
END SUBROUTINE Vectoroffdiag_init
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_setSize
!! NAME
!!  Vectoroffdiag_setSize
!!
!! FUNCTION
!!  impose size
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=vector
!!  new_tail=new_size
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_setSize(vector_1,new_tail)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_setSize'
!End of the abilint section

  TYPE(Vectoroffdiag), INTENT(INOUT) :: vector_1
  INTEGER     , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                     :: size

  IF ( .NOT. ALLOCATED(vector_1%vec) ) THEN
    CALL Vectoroffdiag_init(vector_1,new_tail)
  ELSE
    size = vector_1%size
    IF( new_tail .GT. size ) CALL Vectoroffdiag_enlarge(vector_1,MAX(Global_SIZE,new_tail-size))
  END IF
  vector_1%tail = new_tail
END SUBROUTINE Vectoroffdiag_setSize  
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_enlarge
!! NAME
!!  Vectoroffdiag_enlarge
!!
!! FUNCTION
!!  enlarge memory size
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=vector
!!  size=memory size to add
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_enlarge(vector_1, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_enlarge'
!End of the abilint section

  TYPE(Vectoroffdiag)     , INTENT(INOUT)        :: vector_1
  INTEGER          , INTENT(IN   )        :: size
!Local variables ------------------------------
  INTEGER                                 :: width
  INTEGER                                 :: tail
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vector_temp 
  INTEGER                                 :: size_val

  IF ( ALLOCATED(vector_1%vec) ) THEN
    FREEIF(vector_temp)
    width = vector_1%size
    tail  = vector_1%tail
    size_val = size 
    MALLOC(vector_temp,(1:tail))
    vector_temp(1:tail) = vector_1%vec(1:tail)
    FREE(vector_1%vec)
    vector_1%size = width + size_val
    MALLOC(vector_1%vec,(1:vector_1%size))
    vector_1%vec(1:tail) = vector_temp(1:tail)
    FREE(vector_temp)
  ELSE
    CALL Vectoroffdiag_init(vector_1, Global_SIZE)
  END IF
END SUBROUTINE Vectoroffdiag_enlarge
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_pushBack
!! NAME
!!  Vectoroffdiag_pushBack
!!
!! FUNCTION
!!  push an element at the end
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=vector
!!  value=value to add
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_pushBack(vector_1, value)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_pushBack'
!End of the abilint section

  TYPE(Vectoroffdiag)    , INTENT(INOUT) :: vector_1
  DOUBLE PRECISION, INTENT(IN   ) :: value
!Local variables ------------------------------
  INTEGER                         :: tail

  IF ( vector_1%size .EQ. 0 ) CALL Vectoroffdiag_init(vector_1, Global_SIZE)
  tail = vector_1%tail
  tail = tail + 1
  IF ( tail .GT. vector_1%size ) CALL Vectoroffdiag_enlarge(vector_1,Global_SIZE)
  vector_1%vec(tail) = value
  vector_1%tail      = tail
END SUBROUTINE Vectoroffdiag_pushBack
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_clear
!! NAME
!!  Vectoroffdiag_clear
!!
!! FUNCTION
!!  Clear vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=vector
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_clear(vector_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_clear'
!End of the abilint section

  TYPE(Vectoroffdiag), INTENT(INOUT) :: vector_1
  vector_1%tail = 0 
END SUBROUTINE Vectoroffdiag_clear
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_print
!! NAME
!!  Vectoroffdiag_print
!!
!! FUNCTION
!!  print vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=vector
!!  ostream=file stream
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_print(vector_1,ostream)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_print'
!End of the abilint section

  TYPE(Vectoroffdiag), INTENT(IN) :: vector_1
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: it1
  CHARACTER(LEN=4 )             :: size
  CHARACTER(LEN=15)             :: string

  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  WRITE(size,'(I4)') vector_1%tail
  WRITE(ostream_val,'(A)') "("
  string ='(1x,1ES10.2)'
  DO it1 = 1, vector_1%tail
    WRITE(ostream_val,string) vector_1%vec(it1)
  END DO
  WRITE(ostream_val,'(A)') ")"
END SUBROUTINE Vectoroffdiag_print
!!***

!!****f* ABINIT/m_Vectoroffdiag/Vectoroffdiag_destroy
!! NAME
!!  Vectoroffdiag_destroy
!!
!! FUNCTION
!!  Destroy vector 
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  vector_1=vector
!!
!! OUTPUT
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

SUBROUTINE Vectoroffdiag_destroy(vector_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Vectoroffdiag_destroy'
!End of the abilint section

  TYPE(Vectoroffdiag), INTENT(INOUT) :: vector_1

  FREEIF(vector_1%vec)

  vector_1%tail     = 0
  vector_1%size     = 0
END SUBROUTINE Vectoroffdiag_destroy
!!***

END MODULE m_Vectoroffdiag
!!***


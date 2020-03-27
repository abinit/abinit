
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_Vector
!! NAME
!!  m_Vector
!! 
!! FUNCTION 
!!  Manage a double precision vector
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
MODULE m_Vector
USE m_Global
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_Vector/Vector
!! NAME
!!  Vector
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE, PUBLIC :: Vector
  INTEGER         :: size
  INTEGER         :: tail
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)         :: vec 
END TYPE Vector
!!***

PUBLIC :: Vector_init
PUBLIC :: Vector_setSize
PUBLIC :: Vector_enlarge
PUBLIC :: Vector_pushBack
PUBLIC :: Vector_clear
PUBLIC :: Vector_print
PUBLIC :: Vector_destroy

CONTAINS
!!***

!!****f* ABINIT/m_Vector/Vector_init
!! NAME
!!  Vector_init
!!
!! FUNCTION
!!  initialize
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_init(this, size)

!Arguments ------------------------------------
  TYPE(Vector)     , INTENT(INOUT) :: this
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  this%size = size_val
  FREEIF(this%vec)
  MALLOC(this%vec,(1:size_val))
  this%tail     = 0 
  this%vec = 0.d0
END SUBROUTINE Vector_init
!!***

!!****f* ABINIT/m_Vector/Vector_setSize
!! NAME
!!  Vector_setSize
!!
!! FUNCTION
!!  impose size
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_setSize(this,new_tail)

!Arguments ------------------------------------
  TYPE(Vector), INTENT(INOUT) :: this
  INTEGER     , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                     :: size

  IF ( .NOT. ALLOCATED(this%vec) ) THEN
    CALL Vector_init(this,new_tail)
  ELSE
    size = this%size
    IF( new_tail .GT. size ) THEN
      CALL Vector_enlarge(this,MAX(Global_SIZE,new_tail-size))
    END IF
  END IF
  this%tail = new_tail
END SUBROUTINE Vector_setSize  
!!***

!!****f* ABINIT/m_Vector/Vector_enlarge
!! NAME
!!  Vector_enlarge
!!
!! FUNCTION
!!  enlarge memory size
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_enlarge(this, size)

!Arguments ------------------------------------
  TYPE(Vector)     , INTENT(INOUT)        :: this
  INTEGER          , INTENT(IN   )        :: size
!Local variables ------------------------------
  INTEGER                                 :: width
  INTEGER                                 :: tail
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: thistemp 
  INTEGER                                 :: size_val

  IF ( ALLOCATED(this%vec) ) THEN
    FREEIF(thistemp)
    width = this%size
    tail  = this%tail
    size_val = size 
    MALLOC(thistemp,(1:tail))
    thistemp(1:tail) = this%vec(1:tail)
    FREE(this%vec)
    this%size = width + size_val
    MALLOC(this%vec,(1:this%size))
    this%vec(1:tail) = thistemp(1:tail)
    FREE(thistemp)
  ELSE
    CALL Vector_init(this, Global_SIZE)
  END IF
END SUBROUTINE Vector_enlarge
!!***

!!****f* ABINIT/m_Vector/Vector_pushBack
!! NAME
!!  Vector_pushBack
!!
!! FUNCTION
!!  push an element at the end
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_pushBack(this, value)

!Arguments ------------------------------------
  TYPE(Vector)    , INTENT(INOUT) :: this
  DOUBLE PRECISION, INTENT(IN   ) :: value
!Local variables ------------------------------
  INTEGER                         :: tail

  IF ( this%size .EQ. 0 ) THEN
    CALL Vector_init(this, Global_SIZE)
  END IF
  tail = this%tail
  tail = tail + 1
  IF ( tail .GT. this%size ) THEN
    CALL Vector_enlarge(this,Global_SIZE)
  END IF
  this%vec(tail) = value
  this%tail      = tail
END SUBROUTINE Vector_pushBack
!!***

!!****f* ABINIT/m_Vector/Vector_clear
!! NAME
!!  Vector_clear
!!
!! FUNCTION
!!  Clear vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_clear(this)

!Arguments ------------------------------------
  TYPE(Vector), INTENT(INOUT) :: this
  this%tail = 0 
END SUBROUTINE Vector_clear
!!***

!!****f* ABINIT/m_Vector/Vector_print
!! NAME
!!  Vector_print
!!
!! FUNCTION
!!  print vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_print(this,ostream)

!Arguments ------------------------------------
  TYPE(Vector), INTENT(IN) :: this
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: it1
  CHARACTER(LEN=4 )             :: size
  CHARACTER(LEN=15)             :: string

  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  WRITE(size,'(I4)') this%tail
  WRITE(ostream_val,'(A)') "("
  string ='(1x,1ES10.2)'
  DO it1 = 1, this%tail
    WRITE(ostream_val,string) this%vec(it1)
  END DO
  WRITE(ostream_val,'(A)') ")"
END SUBROUTINE Vector_print
!!***

!!****f* ABINIT/m_Vector/Vector_destroy
!! NAME
!!  Vector_destroy
!!
!! FUNCTION
!!  Destroy vector 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=vector
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

SUBROUTINE Vector_destroy(this)

!Arguments ------------------------------------
  TYPE(Vector), INTENT(INOUT) :: this

  FREEIF(this%vec)

  this%tail     = 0
  this%size     = 0
END SUBROUTINE Vector_destroy
!!***

END MODULE m_Vector
!!***


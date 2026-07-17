
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_VectorComplex
!! NAME
!!  m_VectorComplex
!!
!! FUNCTION
!!  Manage a double precision vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#include "defs.h"
MODULE m_VectorComplex
USE m_Global
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_VectorComplex/VectorComplex
!! NAME
!!  VectorComplex
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE, PUBLIC :: VectorComplex
  INTEGER         :: size
  INTEGER         :: tail
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:)         :: vec
END TYPE VectorComplex
!!***

PUBLIC :: VectorComplex_init
PUBLIC :: VectorComplex_setSize
PUBLIC :: VectorComplex_enlarge
PUBLIC :: VectorComplex_pushBack
PUBLIC :: VectorComplex_clear
PUBLIC :: VectorComplex_print
PUBLIC :: VectorComplex_destroy

CONTAINS
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_init
!! NAME
!!  VectorComplex_init
!!
!! FUNCTION
!!  initialize
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_init(this, size)

!Arguments ------------------------------------
  TYPE(VectorComplex)     , INTENT(INOUT) :: this
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  this%size = size_val
  FREEIF(this%vec)
  MALLOC(this%vec,(1:size_val))
  this%tail     = 0
  this%vec = cmplx(0.d0,0.d0,kind=8)
END SUBROUTINE VectorComplex_init
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_setSize
!! NAME
!!  VectorComplex_setSize
!!
!! FUNCTION
!!  impose size
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_setSize(this,new_tail)

!Arguments ------------------------------------
  TYPE(VectorComplex), INTENT(INOUT) :: this
  INTEGER     , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                     :: size

  IF ( .NOT. ALLOCATED(this%vec) ) THEN
    CALL VectorComplex_init(this,new_tail)
  ELSE
    size = this%size
    IF( new_tail .GT. size ) THEN
      CALL VectorComplex_enlarge(this,MAX(Global_SIZE,new_tail-size))
    END IF
  END IF
  this%tail = new_tail
END SUBROUTINE VectorComplex_setSize
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_enlarge
!! NAME
!!  VectorComplex_enlarge
!!
!! FUNCTION
!!  enlarge memory size
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_enlarge(this, size)

!Arguments ------------------------------------
  TYPE(VectorComplex)     , INTENT(INOUT)        :: this
  INTEGER          , INTENT(IN   )        :: size
!Local variables ------------------------------
  INTEGER                                 :: width
  INTEGER                                 :: tail
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: thistemp
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
    CALL VectorComplex_init(this, Global_SIZE)
  END IF
END SUBROUTINE VectorComplex_enlarge
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_pushBack
!! NAME
!!  VectorComplex_pushBack
!!
!! FUNCTION
!!  push an element at the end
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_pushBack(this, value)

!Arguments ------------------------------------
  TYPE(VectorComplex)    , INTENT(INOUT) :: this
  COMPLEX(KIND=8), INTENT(IN   ) :: value
!Local variables ------------------------------
  INTEGER                         :: tail

  IF ( this%size .EQ. 0 ) THEN
    CALL VectorComplex_init(this, Global_SIZE)
  END IF
  tail = this%tail
  tail = tail + 1
  IF ( tail .GT. this%size ) THEN
    CALL VectorComplex_enlarge(this,Global_SIZE)
  END IF
  this%vec(tail) = value
  this%tail      = tail
END SUBROUTINE VectorComplex_pushBack
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_clear
!! NAME
!!  VectorComplex_clear
!!
!! FUNCTION
!!  Clear vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_clear(this)

!Arguments ------------------------------------
  TYPE(VectorComplex), INTENT(INOUT) :: this
  this%tail = 0
END SUBROUTINE VectorComplex_clear
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_print
!! NAME
!!  VectorComplex_print
!!
!! FUNCTION
!!  print vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_print(this,ostream)

!Arguments ------------------------------------
  TYPE(VectorComplex), INTENT(IN) :: this
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
END SUBROUTINE VectorComplex_print
!!***

!!****f* ABINIT/m_VectorComplex/VectorComplex_destroy
!! NAME
!!  VectorComplex_destroy
!!
!! FUNCTION
!!  Destroy vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
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
!! SOURCE

SUBROUTINE VectorComplex_destroy(this)

!Arguments ------------------------------------
  TYPE(VectorComplex), INTENT(INOUT) :: this

  FREEIF(this%vec)

  this%tail     = 0
  this%size     = 0
END SUBROUTINE VectorComplex_destroy
!!***

END MODULE m_VectorComplex
!!***


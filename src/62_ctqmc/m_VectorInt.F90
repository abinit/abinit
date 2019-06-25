
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_VectorInt
!! NAME
!!  m_VectorInt
!! 
!! FUNCTION 
!!  Manage an integer vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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
MODULE m_VectorInt
USE m_Global
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_VectorInt/VectorInt
!! NAME
!!  VectorInt
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE, PUBLIC :: VectorInt
  INTEGER _PRIVATE :: size
  INTEGER          :: tail
  INTEGER, ALLOCATABLE, DIMENSION(:)         :: vec
END TYPE VectorInt
!!***

PUBLIC  :: VectorInt_init
PUBLIC  :: VectorInt_setSize
PRIVATE :: VectorInt_enlarge
PUBLIC  :: VectorInt_pushBack
PUBLIC  :: VectorInt_clear
PUBLIC  :: VectorInt_print
PUBLIC  :: VectorInt_destroy

CONTAINS
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_init
!! NAME
!!  VectorInt_init
!!
!! FUNCTION
!!  initialize
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_init(this, size)

!Arguments ------------------------------------
  TYPE(VectorInt)     , INTENT(INOUT) :: this
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  this%size = size_val
  FREEIF(this%vec)
  MALLOC(this%vec,(1:size_val))
  this%tail = 0 
  this%vec  = 0
END SUBROUTINE VectorInt_init
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_setSize
!! NAME
!!  VectorInt_setSize
!!
!! FUNCTION
!!  impose size
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_setSize(this,new_tail)

!Arguments ------------------------------------
  TYPE(VectorInt), INTENT(INOUT) :: this
  INTEGER     , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                     :: size

  IF ( .NOT. ALLOCATED(this%vec) ) THEN
    CALL VectorInt_init(this,new_tail)
  ELSE
    size = this%size
    IF( new_tail .GT. size ) THEN
      CALL VectorInt_enlarge(this,MAX(new_tail-size,Global_SIZE))
    END IF
  END IF
  this%tail = new_tail
END SUBROUTINE VectorInt_setSize  
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_enlarge
!! NAME
!!  VectorInt_enlarge
!!
!! FUNCTION
!!  enlarge memory size
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_enlarge(this, size)

!Arguments ------------------------------------
  TYPE(VectorInt)     , INTENT(INOUT)        :: this
  INTEGER             , INTENT(IN   )        :: size
!Local variables ------------------------------
  INTEGER                                 :: width
  INTEGER                                 :: tail
  INTEGER, ALLOCATABLE, DIMENSION(:) :: thistemp 
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
    CALL VectorInt_init(this, Global_SIZE)
  END IF
END SUBROUTINE VectorInt_enlarge
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_pushBack
!! NAME
!!  VectorInt_pushBack
!!
!! FUNCTION
!!  push an element at the end
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_pushBack(this, value)

!Arguments ------------------------------------
  TYPE(VectorInt)    , INTENT(INOUT) :: this
  INTEGER, INTENT(IN   ) :: value
!Local variables ------------------------------
  INTEGER                         :: tail

  IF ( this%size .EQ. 0 ) THEN
    CALL VectorInt_init(this, Global_SIZE)
  END IF
  tail = this%tail
  tail = tail + 1
  IF ( tail .GT. this%size ) THEN
    CALL VectorInt_enlarge(this,Global_SIZE)
  END IF
  this%vec(tail) = value
  this%tail      = tail
END SUBROUTINE VectorInt_pushBack
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_clear
!! NAME
!!  VectorInt_clear
!!
!! FUNCTION
!!  Clear vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_clear(this)

!Arguments ------------------------------------
  TYPE(VectorInt), INTENT(INOUT) :: this
  this%tail = 0 
END SUBROUTINE VectorInt_clear
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_print
!! NAME
!!  VectorInt_print
!!
!! FUNCTION
!!  print vector
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_print(this,ostream)

!Arguments ------------------------------------
  TYPE(VectorInt), INTENT(IN) :: this
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
END SUBROUTINE VectorInt_print
!!***

!!****f* ABINIT/m_VectorInt/VectorInt_destroy
!! NAME
!!  VectorInt_destroy
!!
!! FUNCTION
!!  Destroy vector 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

SUBROUTINE VectorInt_destroy(this)

!Arguments ------------------------------------
  TYPE(VectorInt), INTENT(INOUT) :: this

  FREEIF(this%vec)

  this%tail     = 0
  this%size     = 0
END SUBROUTINE VectorInt_destroy
!!***

END MODULE m_VectorInt
!!***


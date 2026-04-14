
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_MapHybComplex
!! NAME
!!  m_MapHybComplex
!!
!! FUNCTION
!!  map template integer/double
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
MODULE m_MapHybComplex
USE m_Global
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_MapHybComplex/MapHybComplex
!! NAME
!!  MapHybComplex
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

TYPE, PUBLIC :: MapHybComplex
  INTEGER _PRIVATE :: size
  INTEGER          :: tail
  INTEGER         , ALLOCATABLE, DIMENSION(:) :: listINT
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: listDBLE
END TYPE MapHybComplex
!!***

INTERFACE MapHybComplex_sort
  MODULE PROCEDURE MapHybComplex_quickSort, MapHybComplex_sort
END INTERFACE

!PUBLIC INTERFACE ASSIGNMENT(=)
!  MODULE PROCEDURE MapHybComplex_assign
!END INTERFACE

PUBLIC  :: MapHybComplex_init
PUBLIC  :: MapHybComplex_setSize
PRIVATE :: MapHybComplex_enlarge
PUBLIC  :: MapHybComplex_assign
PUBLIC  :: MapHybComplex_sort
PUBLIC  :: MapHybComplex_quickSort
PUBLIC  :: MapHybComplex_print
PUBLIC  :: MapHybComplex_clear
PUBLIC  :: MapHybComplex_destroy

CONTAINS
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_init
!! NAME
!!  MapHybComplex_init
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
!!  this=Map
!!  size=memory size for initialization
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_init(this, size)

!Arguments ------------------------------------
  TYPE(MapHybComplex)     , INTENT(INOUT) :: this
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  this%size = size_val
  FREEIF(this%listINT)
  MALLOC(this%listINT,(1:size_val))
  FREEIF(this%listDBLE)
  MALLOC(this%listDBLE,(1:size_val))
  this%tail     = 0
END SUBROUTINE MapHybComplex_init
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_setSize
!! NAME
!!  MapHybComplex_setSize
!!
!! FUNCTION
!!  impose size of the this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!  new_tail=new size
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_setSize(this,new_tail)

!Arguments ------------------------------------
  TYPE(MapHybComplex), INTENT(INOUT) :: this
  INTEGER     , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                     :: size

  IF ( .NOT. ALLOCATED(this%listINT) ) THEN
    CALL MapHybComplex_init(this)
  END IF
  size = this%size
  IF( new_tail .GT. size ) THEN
    CALL MapHybComplex_enlarge(this, MAX(new_tail-size,Global_SIZE))
  END IF
  this%tail = new_tail
END SUBROUTINE MapHybComplex_setSize
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_enlarge
!! NAME
!!  MapHybComplex_enlarge
!!
!! FUNCTION
!!  enlarge memory space
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!  size=new memory size
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_enlarge(this, size)

!Arguments ------------------------------------
  TYPE(MapHybComplex)     , INTENT(INOUT)       :: this
  INTEGER, OPTIONAL, INTENT(IN   )       :: size
!Local variables ------------------------------
  INTEGER                                :: width
  INTEGER                                :: tail
  INTEGER                                :: size_val
  INTEGER         , ALLOCATABLE, DIMENSION(:) :: listINT_temp
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: listDBLE_temp

  IF ( ALLOCATED(this%listINT) ) THEN
    FREEIF(listINT_temp)
    width = this%size
    tail  = this%tail
    size_val = width
    IF ( PRESENT(size) ) size_val = size
    ! listINT enlarge
    MALLOC(listINT_temp,(1:tail))
    listINT_temp(1:tail) = this%listINT(1:tail)
    FREE(this%listINT)
    this%size = width + size_val
    MALLOC(this%listINT,(1:this%size))
    this%listINT(1:tail) = listINT_temp(1:tail)
    FREE(listINT_temp)
    ! listDBLE enlarge
    MALLOC(listDBLE_temp,(1:tail))
    listDBLE_temp(1:tail) = this%listDBLE(1:tail)
    FREE(this%listDBLE)
    MALLOC(this%listDBLE,(1:this%size))
    this%listDBLE(1:tail) = listDBLE_temp(1:tail)
    FREE(listDBLE_temp)
  ELSE
    CALL MapHybComplex_init(this, Global_SIZE)
  END IF
END SUBROUTINE MapHybComplex_enlarge
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_assign
!! NAME
!!  MapHybComplex_assign
!!
!! FUNCTION
!!  assign this=map
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=this
!!  this=Map
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_assign(this, map)

!Arguments ------------------------------------
  TYPE(MapHybComplex), INTENT(INOUT) :: this
  TYPE(MapHybComplex), INTENT(IN   ) :: map
!Local variables ------------------------------
  INTEGER                     :: tail

  tail = map%tail
  CALL MapHybComplex_setSize(this, tail)
  this%listINT(1:tail)  = map%listINT(1:tail)
  this%listDBLE(1:tail) = map%listDBLE(1:tail)

END SUBROUTINE MapHybComplex_assign
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_sort
!! NAME
!!  MapHybComplex_sort
!!
!! FUNCTION
!!  sort the this with respect to the integer array
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_sort(this)

!Arguments ------------------------------------
  TYPE(MapHybComplex), INTENT(INOUT) :: this

  IF ( this%tail .EQ. 1 ) RETURN
  CALL MapHybComplex_quickSort(this, 1, this%tail)
END SUBROUTINE MapHybComplex_sort
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_quickSort
!! NAME
!!  MapHybComplex_quickSort
!!
!! FUNCTION
!!  sort the this with respect to the integer array
!!  with the quickSort algo
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!  begin=first element to consider
!!  end=last element to consider
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

RECURSIVE SUBROUTINE MapHybComplex_quickSort(this, begin, end)

!Arguments ------------------------------------
  TYPE(MapHybComplex), INTENT(INOUT) :: this
  INTEGER     , INTENT(IN   ) :: begin
  INTEGER     , INTENT(IN   ) :: end
!Local variables ------------------------------
  INTEGER                     :: it1
  INTEGER                     :: it2
  INTEGER                     :: pivot
  INTEGER                     :: Iswap
  COMPLEX(KIND=8)             :: Dswap

  pivot = this%listINT((end-begin)/2 + begin) ! not the betterchoice.... FIXME
  it1 = begin
  it2 = end
  DO WHILE (it1 .LE. it2)
    DO WHILE ( this%listINT(it1) .LT. pivot )
      it1 = it1 + 1
    END DO
    DO WHILE ( this%listINT(it2) .GT. pivot )
      it2 = it2 - 1
    END DO
    IF ( it1 .LE. it2) THEN
      Iswap = this%listINT(it1)
      Dswap = this%listDBLE(it1)
      this%listINT(it1)  = this%listINT(it2)
      this%listDBLE(it1) = this%listDBLE(it2)
      this%listINT(it2)  = Iswap
      this%listDBLE(it2) = Dswap
      it1 = it1 + 1
      it2 = it2 - 1
    END IF
  END DO
  IF ( begin < it2 ) THEN
    CALL MapHybComplex_quickSort(this,begin,it2)
  END IF
  !!it2= it1+1
  IF ( it1 < end ) THEN
    CALL MapHybComplex_quickSort(this,it1,end)
  END IF

END SUBROUTINE MapHybComplex_quickSort
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_print
!! NAME
!!  MapHybComplex_print
!!
!! FUNCTION
!!  print the this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!  ostream=file stream
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_print(this,ostream)

!Arguments ------------------------------------
  TYPE(MapHybComplex)     , INTENT(IN) :: this
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: it

  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  WRITE(ostream_val,'(A,2x,A5,2x,A5)') "#","Index", "Value"
  DO it = 1, this%tail
    WRITE(ostream_val,'(3x,I5,2x,ES22.14)') this%listINT(it), this%listDBLE(it)
  END DO
END SUBROUTINE MapHybComplex_print
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_clear
!! NAME
!!  MapHybComplex_clear
!!
!! FUNCTION
!!  Clear the this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_clear(this)

!Arguments ------------------------------------
  TYPE(MapHybComplex), INTENT(INOUT) :: this
  this%tail = 0
END SUBROUTINE MapHybComplex_clear
!!***

!!****f* ABINIT/m_MapHybComplex/MapHybComplex_destroy
!! NAME
!!  MapHybComplex_destroy
!!
!! FUNCTION
!!  destroy and deallocate the this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Map
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MapHybComplex_destroy(this)

!Arguments ------------------------------------
  TYPE(MapHybComplex), INTENT(INOUT) :: this

  FREEIF(this%listINT)
  FREEIF(this%listDBLE)

  this%tail     = 0
  this%size     = 0
END SUBROUTINE MapHybComplex_destroy
!!***

END MODULE m_MapHybComplex
!!***



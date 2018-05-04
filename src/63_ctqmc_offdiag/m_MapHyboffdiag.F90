#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_MapHyboffdiag
!! NAME
!!  m_MapHyboffdiag
!! 
!! FUNCTION 
!!  map template integer/double
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
MODULE m_MapHyboffdiag
USE m_Global
IMPLICIT NONE


!!***

!!****t* m_MapHyboffdiag/MapHyboffdiag
!! NAME
!!  MapHyboffdiag
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

TYPE MapHyboffdiag
  INTEGER :: size
  INTEGER :: tail
  INTEGER         , ALLOCATABLE, DIMENSION(:) :: listINT
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: listDBLE
END TYPE MapHyboffdiag
!!***

INTERFACE MapHyboffdiag_sort
  MODULE PROCEDURE MapHyboffdiag_quickSort, MapHyboffdiag_sort
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE MapHyboffdiag_assign
END INTERFACE

CONTAINS
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_init
!! NAME
!!  MapHyboffdiag_init
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
!!  map=map
!!  size=memory size for initialization
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

SUBROUTINE MapHyboffdiag_init(map, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_init'
!End of the abilint section

  TYPE(MapHyboffdiag)     , INTENT(INOUT) :: map
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  map%size = size_val
  FREEIF(map%listINT)
  MALLOC(map%listINT,(1:size_val))
  map%listINT=0
  FREEIF(map%listDBLE)
  MALLOC(map%listDBLE,(1:size_val))
  map%listDBLE=0.d0
  map%tail     = 0
END SUBROUTINE MapHyboffdiag_init
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_setSize
!! NAME
!!  MapHyboffdiag_setSize
!!
!! FUNCTION
!!  impose size of the map
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
!!  new_tail=new size
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

SUBROUTINE MapHyboffdiag_setSize(map,new_tail)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_setSize'
!End of the abilint section

  TYPE(MapHyboffdiag), INTENT(INOUT) :: map
  INTEGER     , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                     :: size

  IF ( .NOT. ALLOCATED(map%listINT) ) CALL MapHyboffdiag_init(map)
  size = map%size
  IF( new_tail .GT. size ) CALL MapHyboffdiag_enlarge(map, MAX(new_tail-size,Global_SIZE))
  map%tail = new_tail
END SUBROUTINE MapHyboffdiag_setSize  
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_enlarge
!! NAME
!!  MapHyboffdiag_enlarge
!!
!! FUNCTION
!!  enlarge memory space
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
!!  size=new memory size
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

SUBROUTINE MapHyboffdiag_enlarge(map, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_enlarge'
!End of the abilint section

  TYPE(MapHyboffdiag)     , INTENT(INOUT)       :: map
  INTEGER, OPTIONAL, INTENT(IN   )       :: size
!Local variables ------------------------------
  INTEGER                                :: width
  INTEGER                                :: tail
  INTEGER                                :: size_val
  INTEGER         , ALLOCATABLE, DIMENSION(:) :: listINT_temp 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: listDBLE_temp 

  IF ( ALLOCATED(map%listINT) ) THEN
    FREEIF(listINT_temp)
    width = map%size
    tail  = map%tail
    size_val = width
    IF ( PRESENT(size) ) size_val = size 
    ! listINT enlarge
    MALLOC(listINT_temp,(1:tail))
    listINT_temp(1:tail) = map%listINT(1:tail)
    FREE(map%listINT)
    map%size = width + size_val
    MALLOC(map%listINT,(1:map%size))
    map%listINT(1:tail) = listINT_temp(1:tail)
    FREE(listINT_temp)
    ! listDBLE enlarge
    MALLOC(listDBLE_temp,(1:tail))
    listDBLE_temp(1:tail) = map%listDBLE(1:tail)
    FREE(map%listDBLE)
    MALLOC(map%listDBLE,(1:map%size))
    map%listDBLE(1:tail) = listDBLE_temp(1:tail)
    FREE(listDBLE_temp)
  ELSE
    CALL MapHyboffdiag_init(map, Global_SIZE)
  END IF
END SUBROUTINE MapHyboffdiag_enlarge
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_assign
!! NAME
!!  MapHyboffdiag_assign
!!
!! FUNCTION
!!  assign
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map_1=map
!!  map_2=map
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

SUBROUTINE MapHyboffdiag_assign(map_1, map_2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_assign'
!End of the abilint section

  TYPE(MapHyboffdiag), INTENT(INOUT) :: map_1
  TYPE(MapHyboffdiag), INTENT(IN   ) :: map_2
!Local variables ------------------------------
  INTEGER                     :: tail

  tail = map_2%tail
  CALL MapHyboffdiag_setSize(map_1, tail)
  map_1%listINT(1:tail) = map_2%listINT(1:tail)
  map_1%listDBLE(1:tail) = map_2%listDBLE(1:tail)

END SUBROUTINE MapHyboffdiag_assign
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_sort
!! NAME
!!  MapHyboffdiag_sort
!!
!! FUNCTION
!!  sort the map with respect to the integer array
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
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

SUBROUTINE MapHyboffdiag_sort(map)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_sort'
!End of the abilint section

  TYPE(MapHyboffdiag), INTENT(INOUT) :: map
 
  IF ( map%tail .EQ. 1 ) RETURN
  CALL MapHyboffdiag_quickSort(map, 1, map%tail)
END SUBROUTINE MapHyboffdiag_sort
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_quickSort
!! NAME
!!  MapHyboffdiag_quickSort
!!
!! FUNCTION
!!  sort the map with respect to the integer array
!!  with the quickSort algo
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
!!  begin=first element to consider
!!  end=last element to consider
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

RECURSIVE SUBROUTINE MapHyboffdiag_quickSort(map, begin, end)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_quickSort'
!End of the abilint section

  TYPE(MapHyboffdiag), INTENT(INOUT) :: map
  INTEGER     , INTENT(IN   ) :: begin
  INTEGER     , INTENT(IN   ) :: end
!Local variables ------------------------------
  INTEGER                     :: it1
  INTEGER                     :: it2
  INTEGER                     :: pivot
  INTEGER                     :: Iswap
  DOUBLE PRECISION            :: Dswap

  pivot = map%listINT((end-begin)/2 + begin) ! not the betterchoice.... FIXME
  it1 = begin
  it2 = end
  DO WHILE (it1 .LE. it2)
    DO WHILE ( map%listINT(it1) .LT. pivot )
      it1 = it1 + 1
    END DO
    DO WHILE ( map%listINT(it2) .GT. pivot )
      it2 = it2 - 1
    END DO
    IF ( it1 .LE. it2) THEN
      Iswap = map%listINT(it1)
      Dswap = map%listDBLE(it1)
      map%listINT(it1)  = map%listINT(it2)
      map%listDBLE(it1) = map%listDBLE(it2)
      map%listINT(it2)  = Iswap
      map%listDBLE(it2) = Dswap
      it1 = it1 + 1
      it2 = it2 - 1
    END IF
  END DO
  IF ( begin < it2 ) CALL MapHyboffdiag_quickSort(map,begin,it2)
  !!it2= it1+1
  IF ( it1 < end ) CALL MapHyboffdiag_quickSort(map,it1,end)

END SUBROUTINE MapHyboffdiag_quickSort
!!***
 
!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_print
!! NAME
!!  MapHyboffdiag_print
!!
!! FUNCTION
!!  print the map
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
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

SUBROUTINE MapHyboffdiag_print(map,ostream)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_print'
!End of the abilint section

  TYPE(MapHyboffdiag)     , INTENT(IN) :: map
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: it

  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  WRITE(ostream_val,'(A,2x,A5,2x,A5)') "#","Index", "Value"
  DO it = 1, map%tail
    WRITE(ostream_val,'(3x,I5,2x,ES22.14)') map%listINT(it), map%listDBLE(it) 
  END DO
END SUBROUTINE MapHyboffdiag_print
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_clear
!! NAME
!!  MapHyboffdiag_clear
!!
!! FUNCTION
!!  Clear the map
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
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

SUBROUTINE MapHyboffdiag_clear(map)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_clear'
!End of the abilint section

  TYPE(MapHyboffdiag), INTENT(INOUT) :: map
  map%tail = 0 
END SUBROUTINE MapHyboffdiag_clear
!!***

!!****f* ABINIT/m_MapHyboffdiag/MapHyboffdiag_destroy
!! NAME
!!  MapHyboffdiag_destroy
!!
!! FUNCTION
!!  destroy and deallocate the map
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  map=Map
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

SUBROUTINE MapHyboffdiag_destroy(map)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MapHyboffdiag_destroy'
!End of the abilint section

  TYPE(MapHyboffdiag), INTENT(INOUT) :: map

  FREEIF(map%listINT)
  FREEIF(map%listDBLE)

  map%tail     = 0
  map%size     = 0
END SUBROUTINE MapHyboffdiag_destroy
!!***

END MODULE m_MapHyboffdiag
!!***



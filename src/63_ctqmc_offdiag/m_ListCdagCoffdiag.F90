#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ListCdagCoffdiag
!! NAME
!!  m_ListCdagCoffdiag
!! 
!! FUNCTION 
!!  Manage a 2D vector to store couple of c+c
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
MODULE m_ListCdagCoffdiag
!USE m_CdagC
USE m_Global
IMPLICIT NONE


!!***

!!****t* m_ListCdagCoffdiag/ListCdagCoffdiag
!! NAME
!!  ListCdagCoffdiag
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

TYPE ListCdagCoffdiag
  INTEGER :: size = 0
!  max size of matrix list

  INTEGER :: tail = 0 
!  the size of matrix list that contains physical data (ie number of
!  segment)

  !DOUBLE PRECISION :: inv_dt = 0.d0
!  TYPE(CdagC), ALLOCATABLE, DIMENSION(:) :: list => NULL()
  !INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: ind
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: list
!  for all elements i below itail, list(i,1:2) are times for creation
!  and destruction of particles.

END TYPE ListCdagCoffdiag
!!***

INTERFACE ListCdagCoffdiag_firstHigher
  MODULE PROCEDURE ListCdagCoffdiag_firstHigherThanReal
END INTERFACE

INTERFACE ListCdagCoffdiag_sort
  MODULE PROCEDURE ListCdagCoffdiag_quickSort, ListCdagCoffdiag_sort
END INTERFACE

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE ListCdagCoffdiag_assign
END INTERFACE

CONTAINS
!!***

!SUBROUTINE ListCdagCoffdiag_init(list_1, inv_dt, size)
!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_init
!! NAME
!!  ListCdagCoffdiag_init
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
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_init(list_1, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_init'
!End of the abilint section

  TYPE(ListCdagCoffdiag)  , INTENT(INOUT) :: list_1
  !DOUBLE PRECISION , INTENT(IN   ) :: inv_dt
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  !list_1%inv_dt = inv_dt
  IF ( PRESENT(size) ) size_val = size
  list_1%size = size_val
  FREEIF(list_1%list)
  MALLOC(list_1%list,(0:size_val,1:2))
  !FAKEFREEIF(list_1%ind)
  !FAKEMALLOC(list_1%ind,(0:size_val,1:2))
  list_1%tail     = 0
END SUBROUTINE ListCdagCoffdiag_init
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_setSize
!! NAME
!!  ListCdagCoffdiag_setSize
!!
!! FUNCTION
!!  Impose size of the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_setSize(list_1,new_tail)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_setSize'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  INTEGER        , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                        :: size

  !IF ( .NOT. ALLOCATED(list_1%list) ) CALL ListCdagCoffdiag_init(list_1, list_1%inv_dt)
  IF ( .NOT. ALLOCATED(list_1%list) ) CALL ListCdagCoffdiag_init(list_1)
  size = list_1%size
  IF( new_tail .GT. size ) CALL ListCdagCoffdiag_enlarge(list_1,MAX(Global_SIZE, new_tail-size))
  list_1%tail = new_tail
END SUBROUTINE ListCdagCoffdiag_setSize  
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_enlarge
!! NAME
!!  ListCdagCoffdiag_enlarge
!!
!! FUNCTION
!!  Enlarge memory space of the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_enlarge(list_1, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_enlarge'
!End of the abilint section

  TYPE(ListCdagCoffdiag),   INTENT(INOUT)       :: list_1
  INTEGER, OPTIONAL, INTENT(IN   )       :: size
!Local variables ------------------------------
  INTEGER                                :: width
  INTEGER                                :: tail
  INTEGER                                :: size_val
  !INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: ind_temp 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: list_temp 

  IF ( ALLOCATED(list_1%list) ) THEN
    FREEIF(list_temp)
    !FAKEFREEIF(ind_temp )
    width = list_1%size
    tail  = list_1%tail
    size_val = width
    IF ( PRESENT(size) ) size_val = size 
    MALLOC(list_temp,(0:tail,1:2))
    !MALLOC( ind_temp,(0:width,1:2))
    list_temp(0:tail,:) = list_1%list(0:tail,:)
    !ind_temp  = list_1%ind
    FREE(list_1%list)
    !FREE(list_1%ind )
    list_1%size = width + size_val
    MALLOC(list_1%list,(0:list_1%size,1:2))
    !MALLOC(list_1%ind ,(0:list_1%size,1:2))
    list_1%list(0:tail,1:2) = list_temp(0:tail,1:2)
    !list_1%ind (0:width,1:2) = ind_temp (0:width,1:2)
    FREE(list_temp)
  ELSE
    !CALL ListCdagCoffdiag_init(list_1, list_1%inv_dt, Global_SIZE)
    CALL ListCdagCoffdiag_init(list_1, Global_SIZE)
  END IF
END SUBROUTINE ListCdagCoffdiag_enlarge
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/listCdagCoffdiag_assign
!! NAME
!!  listCdagCoffdiag_assign
!!
!! FUNCTION
!!  assign routine
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  list_2=ListCdagCoffdiag
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

SUBROUTINE listCdagCoffdiag_assign(list_1, list_2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'listCdagCoffdiag_assign'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  TYPE(ListCdagCoffdiag), INTENT(IN   ) :: list_2
!Local variables ------------------------------
  INTEGER                        :: tail

  tail = list_2%tail
  CALL ListCdagCoffdiag_setSize(list_1, tail)
  list_1%list(0:tail,1:2) = list_2%list(0:tail,1:2)
  !list_1%ind (0:tail,1:2) = list_2%ind (0:tail,1:2)

END SUBROUTINE ListCdagCoffdiag_assign
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_swap
!! NAME
!!  ListCdagCoffdiag_swap
!!
!! FUNCTION
!!  Swap two lists
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  list_2=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_swap(list_1,list_2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_swap'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_2
!Local variables ------------------------------
  INTEGER :: tail1
  INTEGER :: tail2
  INTEGER :: i
  INTEGER :: j
  !INTEGER         , DIMENSION(1:2) :: ind_tmp
  DOUBLE PRECISION, DIMENSION(1:2) :: CdagC_tmp

  tail1 = list_1%tail
  tail2 = list_2%tail

  i = MAX(tail1,tail2)
  IF ( list_1%size .LT. i ) CALL ListCdagCoffdiag_enlarge(list_1,i)
  IF ( list_2%size .LT. i ) CALL ListCdagCoffdiag_enlarge(list_2,i)

  DO j = 0, i
    CdagC_tmp(1:2) = list_1%list(j,1:2)
    !ind_tmp  (1:2) = list_1%ind (j,1:2)
    list_1%list(j,1:2) = list_2%list(j,1:2)
    !list_1%ind (j,1:2) = list_2%ind (j,1:2)
    list_2%list(j,1:2) = CdagC_tmp(1:2)
    !list_2%ind (j,1:2) = ind_tmp  (1:2)
  END DO
  list_2%tail = tail1
  list_1%tail = tail2
END SUBROUTINE ListCdagCoffdiag_swap
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_pushBack
!! NAME
!!  ListCdagCoffdiag_pushBack
!!
!! FUNCTION
!!  push at the end of the list a couple
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  CdagC_1=couple
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

SUBROUTINE ListCdagCoffdiag_pushBack(list_1, CdagC_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_pushBack'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT)       :: list_1
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1
!Local variables ------------------------------
  INTEGER                              :: tail

  !IF ( list_1%size .EQ. 0 ) CALL ListCdagCoffdiag_init(list_1, list_1%inv_dt, Global_SIZE)
  IF ( list_1%size .EQ. 0 ) CALL ListCdagCoffdiag_init(list_1, Global_SIZE)
  tail = list_1%tail
  tail = tail + 1
  IF ( tail .GT. list_1%size ) CALL ListCdagCoffdiag_enlarge(list_1)
  list_1%list(tail,1:2) = CdagC_1
  !list_1%ind (tail,Cdag_) = INT(CdagC_1(Cdag_) * list_1%inv_dt + 0.5d0)
  !list_1%ind (tail,C_   ) = INT(CdagC_1(C_   ) * list_1%inv_dt + 0.5d0)
  list_1%tail       = tail
END SUBROUTINE ListCdagCoffdiag_pushBack
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_insert
!! NAME
!!  ListCdagCoffdiag_insert
!!
!! FUNCTION
!!  insert somewhere a couple
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  CdagC_1=couple
!!  position=where to insert
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

SUBROUTINE ListCdagCoffdiag_insert(list_1, CdagC_1, position)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_insert'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1 
  INTEGER        , INTENT(IN   ) :: position
!Local variables ------------------------------
  INTEGER                        :: new_position
  INTEGER                        :: tail
  
  tail         = list_1%tail + 1
  new_position = position 
  IF ( tail .GT. list_1%size ) CALL ListCdagCoffdiag_enlarge(list_1)
  IF ( position .EQ. -1 ) THEN
    new_position = tail
  ELSE IF ( position .LE. tail ) THEN
  ! new_position = position 
    list_1%list(tail:position+1:-1,1:2) = list_1%list(list_1%tail:position:-1,1:2)
  ELSE 
    CALL ERROR("ListCdagCoffdiag_insert : position > tail                ")
  END IF
  
  list_1%list(new_position,1:2) = CdagC_1
  !list_1%ind (new_position,Cdag_) = INT(CdagC_1(Cdag_) * list_1%inv_dt + 0.5d0)
  !list_1%ind (new_position,C_   ) = INT(CdagC_1(C_   ) * list_1%inv_dt + 0.5d0)
  list_1%tail = tail
END SUBROUTINE ListCdagCoffdiag_insert
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_popBack
!! NAME
!!  ListCdagCoffdiag_popBack
!!
!! FUNCTION
!!  Remove the last element
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_popBack(list_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_popBack'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
!Local variables ------------------------------
  INTEGER                        :: tail

  tail = list_1%tail
  IF ( tail .EQ. 0 ) RETURN
  list_1%tail = tail - 1
END SUBROUTINE ListCdagCoffdiag_popBack
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_erase
!! NAME
!!  ListCdagCoffdiag_erase
!!
!! FUNCTION
!!  Erase a couple at a given position
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  position=position of the element to remove
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

SUBROUTINE ListCdagCoffdiag_erase(list_1,position)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_erase'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  INTEGER,         INTENT(IN   ) :: position
!Local variables ------------------------------
  INTEGER                        :: tail
  INTEGER                        :: new_tail
  INTEGER                        :: continueing
  
  tail = list_1%tail
  IF ( position .GT. tail ) &
    CALL ERROR("ListCdagCoffdiag_erase : position > tail                 ")
  new_tail    = tail - 1
  continueing = position + 1
  list_1%list(new_tail:position:-1,1:2) = list_1%list(tail:continueing:-1,1:2)
  list_1%tail = new_tail
END SUBROUTINE ListCdagCoffdiag_erase
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_firstHigherThanReal
!! NAME
!!  ListCdagCoffdiag_firstHigherThanReal
!!
!! FUNCTION
!!  search for the first element higher than the real time
!!  assume the list is already sorted
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  time=reference
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

INTEGER FUNCTION ListCdagCoffdiag_firstHigherThanReal(list_1, time)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_firstHigherThanReal'
!End of the abilint section

  TYPE(ListCdagCoffdiag),  INTENT(IN) :: list_1
  DOUBLE PRECISION, INTENT(IN) :: time
#include "ListCdagCoffdiag_firstHigher.h"
  ! Dichotomic research
#include "ListCdagCoffdiag_firstHigher"
!  unefficient function for long list  
!  it = 1
!  DO WHILE ( it .LE. list_1%tail .AND. list_1%list(it) .LE. value )
!    it = it + 1
!  END DO
!  IF ( it .GT. list_1%tail ) it = -1
!  ListCdagCoffdiag_firstHigherThanReal = it
  ListCdagCoffdiag_firstHigherThanReal = firstHigher
END FUNCTION ListCdagCoffdiag_firstHigherThanReal
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_sort
!! NAME
!!  ListCdagCoffdiag_sort
!!
!! FUNCTION
!!  sort the list by c+ increasing
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_sort(list_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_sort'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
 
  IF ( list_1%tail .EQ. 1 ) RETURN
  CALL ListCdagCoffdiag_quickSort(list_1, 1, list_1%tail)
END SUBROUTINE ListCdagCoffdiag_sort
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_quickSort
!! NAME
!!  ListCdagCoffdiag_quickSort
!!
!! FUNCTION
!!  sort the list by c+ increasing with the quick sort algo
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  begin=from element to consider
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

RECURSIVE SUBROUTINE ListCdagCoffdiag_quickSort(list_1, begin, end)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_quickSort'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  INTEGER,         INTENT(IN   ) :: begin
  INTEGER,         INTENT(IN   ) :: end
!Local variables k-----------------------------
  INTEGER                        :: it1
  INTEGER                        :: it2
  DOUBLE PRECISION               :: pivot
  DOUBLE PRECISION, DIMENSION(1:2):: CdagC_swap
  !DOUBLE PRECISION, DIMENSION(1:2):: ind_swap

  pivot = list_1%list((end-begin)/2 + begin,Cdag_) ! not the betterchoice.... FIXME
  it1 = begin
  it2 = end
  DO WHILE (it1 .LE. it2)
    DO WHILE ( list_1%list(it1,Cdag_) .LT. pivot )
      it1 = it1 + 1
    END DO
    DO WHILE ( list_1%list(it2,Cdag_) .GT. pivot )
      it2 = it2 - 1
    END DO
    IF ( it1 .LE. it2) THEN
      CdagC_swap = list_1%list(it1,1:2)
      !ind_swap   = list_1%ind (it1,1:2)
      list_1%list(it1,1:2) = list_1%list(it2,1:2)
      !list_1%ind (it1,1:2) = list_1%ind (it2,1:2)
      list_1%list(it2,1:2) = CdagC_swap
      !list_1%ind (it2,1:2) = ind_swap
      it1 = it1 + 1
      it2 = it2 - 1
    END IF
  END DO
  IF ( begin < it2 ) CALL ListCdagCoffdiag_quickSort(list_1,begin,it2)
  !!it2= it1+1
  IF ( it1 < end ) CALL ListCdagCoffdiag_quickSort(list_1,it1,end)

END SUBROUTINE ListCdagCoffdiag_quickSort
!!***
 
!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_print
!! NAME
!!  ListCdagCoffdiag_print
!!
!! FUNCTION
!!  print the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
!!  ostrean=file stream
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

SUBROUTINE ListCdagCoffdiag_print(list_1,ostream)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_print'
!End of the abilint section

  TYPE(ListCdagCoffdiag)  , INTENT(IN) :: list_1
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: it

  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  WRITE(ostream_val,'(A,2x,A4,22x,A)') "#","Cdag", "C"
  DO it = 1, list_1%tail
    WRITE(ostream_val,*) list_1%list(it,Cdag_), list_1%list(it,C_) 
  END DO
END SUBROUTINE ListCdagCoffdiag_print
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_clear
!! NAME
!!  ListCdagCoffdiag_clear
!!
!! FUNCTION
!!  Clear the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_clear(list_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_clear'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1
  list_1%tail = 0 
END SUBROUTINE ListCdagCoffdiag_clear
!!***

!!****f* ABINIT/m_ListCdagCoffdiag/ListCdagCoffdiag_destroy
!! NAME
!!  ListCdagCoffdiag_destroy
!!
!! FUNCTION
!!  destroy the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagCoffdiag
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

SUBROUTINE ListCdagCoffdiag_destroy(list_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ListCdagCoffdiag_destroy'
!End of the abilint section

  TYPE(ListCdagCoffdiag), INTENT(INOUT) :: list_1

  FREEIF(list_1%list)
  !FAKEFREEIF(list_1%ind )

  list_1%tail     = 0
  list_1%size     = 0
  !list_1%inv_dt   = 0.d0
END SUBROUTINE ListCdagCoffdiag_destroy
!!***

END MODULE m_ListCdagCoffdiag
!!***


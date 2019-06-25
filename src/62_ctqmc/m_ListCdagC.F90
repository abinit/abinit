
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ListCdagC
!! NAME
!!  m_ListCdagC
!! 
!! FUNCTION 
!!  Manage a 2D vector to store couple of c+c
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
MODULE m_ListCdagC
USE m_Global

IMPLICIT NONE

!!***

PRIVATE

!!****t* m_ListCdagC/ListCdagC
!! NAME
!!  ListCdagC
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

TYPE, PUBLIC :: ListCdagC
  INTEGER _PRIVATE :: size = 0
!  max size of matrix list

  INTEGER          :: tail = 0 
!  the size of matrix list that contains physical data (ie number of
!  segment)
  !DOUBLE PRECISION :: inv_dt = 0.d0
!  TYPE(CdagC), ALLOCATABLE, DIMENSION(:) :: list => NULL()
  !INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: ind
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: list
!  for all elements i below itail, list(i,1:2) are times for creation
!  and destruction of particles.
END TYPE ListcdagC
!!***

INTERFACE ListCdagC_firstHigher
  MODULE PROCEDURE ListCdagC_firstHigherThanReal
END INTERFACE

INTERFACE ListCdagC_sort
  MODULE PROCEDURE ListCdagC_quickSort, ListCdagC_sort
END INTERFACE

!INTERFACE ASSIGNMENT(=)
!  MODULE PROCEDURE ListCdagC_assign
!END INTERFACE

PUBLIC  :: ListCdagC_init
PUBLIC  :: ListCdagC_setSize
PRIVATE :: ListCdagC_enlarge
PUBLIC  :: listCdagC_assign
PUBLIC  :: ListCdagC_swap
PUBLIC  :: ListCdagC_pushBack
PUBLIC  :: ListCdagC_insert
PUBLIC  :: ListCdagC_popBack
PUBLIC  :: ListCdagC_erase
PUBLIC  :: ListCdagC_firstHigher
PUBLIC  :: ListCdagC_sort
PUBLIC  :: ListCdagC_quickSort
PUBLIC  :: ListCdagC_print
PUBLIC  :: ListCdagC_clear
PUBLIC  :: ListCdagC_destroy

CONTAINS
!!***

!SUBROUTINE ListCdagC_init(this, inv_dt, size)
!!****f* ABINIT/m_ListCdagC/ListCdagC_init
!! NAME
!!  ListCdagC_init
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
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_init(this, size)

!Arguments ------------------------------------
  TYPE(ListCdagC)  , INTENT(INOUT) :: this
  !DOUBLE PRECISION , INTENT(IN   ) :: inv_dt
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  !this%inv_dt = inv_dt
  IF ( PRESENT(size) ) size_val = size
  this%size = size_val
  FREEIF(this%list)
  MALLOC(this%list,(0:size_val,1:2))
  !FAKEFREEIF(this%ind)
  !FAKEMALLOC(this%ind,(0:size_val,1:2))
  this%tail     = 0
END SUBROUTINE ListCdagC_init
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_setSize
!! NAME
!!  ListCdagC_setSize
!!
!! FUNCTION
!!  Impose size of the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_setSize(this,new_tail)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  INTEGER        , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                        :: size

  !IF ( .NOT. ALLOCATED(this%list) ) THEN
  !  CALL ListCdagC_init(this, this%inv_dt)
  !END IF
  IF ( .NOT. ALLOCATED(this%list) ) THEN
    CALL ListCdagC_init(this)
  END IF
  size = this%size
  IF( new_tail .GT. size ) THEN
    CALL ListCdagC_enlarge(this,MAX(Global_SIZE, new_tail-size))
  END IF
  this%tail = new_tail
END SUBROUTINE ListCdagC_setSize  
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_enlarge
!! NAME
!!  ListCdagC_enlarge
!!
!! FUNCTION
!!  Enlarge memory space of the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_enlarge(this, size)

!Arguments ------------------------------------
  TYPE(ListCdagC),   INTENT(INOUT)       :: this
  INTEGER, OPTIONAL, INTENT(IN   )       :: size
!Local variables ------------------------------
  INTEGER                                :: width
  INTEGER                                :: tail
  INTEGER                                :: size_val
  !INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: ind_temp 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: list_temp 

  IF ( ALLOCATED(this%list) ) THEN
    FREEIF(list_temp)
    !FAKEFREEIF(ind_temp )
    width = this%size
    tail  = this%tail
    size_val = width
    IF ( PRESENT(size) ) size_val = size 
    MALLOC(list_temp,(0:tail,1:2))
    !MALLOC( ind_temp,(0:width,1:2))
    list_temp(0:tail,:) = this%list(0:tail,:)
    !ind_temp  = this%ind
    FREE(this%list)
    !FREE(this%ind )
    this%size = width + size_val
    MALLOC(this%list,(0:this%size,1:2))
    !MALLOC(this%ind ,(0:this%size,1:2))
    this%list(0:tail,1:2) = list_temp(0:tail,1:2)
    !this%ind (0:width,1:2) = ind_temp (0:width,1:2)
    FREE(list_temp)
  ELSE
    !CALL ListCdagC_init(this, this%inv_dt, Global_SIZE)
    CALL ListCdagC_init(this, Global_SIZE)
  END IF
END SUBROUTINE ListCdagC_enlarge
!!***

!!****f* ABINIT/m_ListCdagC/listCdagC_assign
!! NAME
!!  listCdagC_assign
!!
!! FUNCTION
!!  assign routine
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
!!  list_2=ListCdagC
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

SUBROUTINE listCdagC_assign(this, list_2)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  TYPE(ListCdagC), INTENT(IN   ) :: list_2
!Local variables ------------------------------
  INTEGER                        :: tail

  tail = list_2%tail
  CALL ListCdagC_setSize(this, tail)
  this%list(0:tail,1:2) = list_2%list(0:tail,1:2)
  !this%ind (0:tail,1:2) = list_2%ind (0:tail,1:2)

END SUBROUTINE ListCdagC_assign
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_swap
!! NAME
!!  ListCdagC_swap
!!
!! FUNCTION
!!  Swap two lists
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
!!  list_2=ListCdagC
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

SUBROUTINE ListCdagC_swap(this,list_2)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  TYPE(ListCdagC), INTENT(INOUT) :: list_2
!Local variables ------------------------------
  INTEGER :: tail1
  INTEGER :: tail2
  INTEGER :: i
  INTEGER :: j
  !INTEGER         , DIMENSION(1:2) :: ind_tmp
  DOUBLE PRECISION, DIMENSION(1:2) :: CdagC_tmp

  tail1 = this%tail
  tail2 = list_2%tail

  i = MAX(tail1,tail2)
  IF ( this%size .LT. i ) THEN
    CALL ListCdagC_enlarge(this,i)
  END IF
  IF ( list_2%size .LT. i ) THEN
    CALL ListCdagC_enlarge(list_2,i)
  END IF

  DO j = 0, i
    CdagC_tmp(1:2) = this%list(j,1:2)
    !ind_tmp  (1:2) = this%ind (j,1:2)
    this%list(j,1:2) = list_2%list(j,1:2)
    !this%ind (j,1:2) = list_2%ind (j,1:2)
    list_2%list(j,1:2) = CdagC_tmp(1:2)
    !list_2%ind (j,1:2) = ind_tmp  (1:2)
  END DO
  list_2%tail = tail1
  this%tail = tail2
END SUBROUTINE ListCdagC_swap
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_pushBack
!! NAME
!!  ListCdagC_pushBack
!!
!! FUNCTION
!!  push at the end of the list a couple
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_pushBack(this, CdagC_1)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT)       :: this
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1
!Local variables ------------------------------
  INTEGER                              :: tail

  !IF ( this%size .EQ. 0 ) THEN
  !  CALL ListCdagC_init(this, this%inv_dt, Global_SIZE)
  !ENDIF
  IF ( this%size .EQ. 0 ) THEN
    CALL ListCdagC_init(this, Global_SIZE)
  END IF
  tail = this%tail
  tail = tail + 1
  IF ( tail .GT. this%size ) THEN
    CALL ListCdagC_enlarge(this)
  END IF
  this%list(tail,1:2) = CdagC_1
  !this%ind (tail,Cdag_) = INT(CdagC_1(Cdag_) * this%inv_dt + 0.5d0)
  !this%ind (tail,C_   ) = INT(CdagC_1(C_   ) * this%inv_dt + 0.5d0)
  this%tail       = tail
END SUBROUTINE ListCdagC_pushBack
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_insert
!! NAME
!!  ListCdagC_insert
!!
!! FUNCTION
!!  insert somewhere a couple
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_insert(this, CdagC_1, position)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1 
  INTEGER        , INTENT(IN   ) :: position
!Local variables ------------------------------
  INTEGER                        :: new_position
  INTEGER                        :: tail
  
  tail         = this%tail + 1
  new_position = position 
  IF ( tail .GT. this%size ) THEN
    CALL ListCdagC_enlarge(this)
  END IF
  IF ( position .EQ. -1 ) THEN
    new_position = tail
  ELSE IF ( position .LE. tail ) THEN
  ! new_position = position 
    this%list(tail:position+1:-1,1:2) = this%list(this%tail:position:-1,1:2)
  ELSE 
    CALL ERROR("ListCdagC_insert : position > tail                ")
  END IF
  
  this%list(new_position,1:2) = CdagC_1
  !this%ind (new_position,Cdag_) = INT(CdagC_1(Cdag_) * this%inv_dt + 0.5d0)
  !this%ind (new_position,C_   ) = INT(CdagC_1(C_   ) * this%inv_dt + 0.5d0)
  this%tail = tail
END SUBROUTINE ListCdagC_insert
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_popBack
!! NAME
!!  ListCdagC_popBack
!!
!! FUNCTION
!!  Remove the last element
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_popBack(this)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                        :: tail

  tail = this%tail
  IF ( tail .EQ. 0 ) RETURN
  this%tail = tail - 1
END SUBROUTINE ListCdagC_popBack
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_erase
!! NAME
!!  ListCdagC_erase
!!
!! FUNCTION
!!  Erase a couple at a given position
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_erase(this,position)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  INTEGER,         INTENT(IN   ) :: position
!Local variables ------------------------------
  INTEGER                        :: tail
  INTEGER                        :: new_tail
  INTEGER                        :: continueing
  
  tail = this%tail
  IF ( position .GT. tail ) &
    CALL ERROR("ListCdagC_erase : position > tail                 ")
  new_tail    = tail - 1
  continueing = position + 1
  this%list(new_tail:position:-1,1:2) = this%list(tail:continueing:-1,1:2)
  this%tail = new_tail
END SUBROUTINE ListCdagC_erase
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_firstHigherThanReal
!! NAME
!!  ListCdagC_firstHigherThanReal
!!
!! FUNCTION
!!  search for the first element higher than the real time
!!  assume the list is already sorted
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

INTEGER FUNCTION ListCdagC_firstHigherThanReal(this, time)

!Arguments ------------------------------------
  TYPE(ListCdagC),  INTENT(IN) :: this
  DOUBLE PRECISION, INTENT(IN) :: time
#include "ListCdagC_firstHigher.h"
  ! Dichotomic research
#define list_1 this
#include "ListCdagC_firstHigher"
#undef list_1
!  unefficient function for long list  
!  it = 1
!  DO WHILE ( it .LE. this%tail .AND. this%list(it) .LE. value )
!    it = it + 1
!  END DO
!  IF ( it .GT. this%tail ) it = -1
!  ListCdagC_firstHigherThanReal = it
  ListCdagC_firstHigherThanReal = firstHigher
END FUNCTION ListCdagC_firstHigherThanReal
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_sort
!! NAME
!!  ListCdagC_sort
!!
!! FUNCTION
!!  sort the list by c+ increasing
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_sort(this)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
 
  IF ( this%tail .EQ. 1 ) RETURN
  CALL ListCdagC_quickSort(this, 1, this%tail)
END SUBROUTINE ListCdagC_sort
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_quickSort
!! NAME
!!  ListCdagC_quickSort
!!
!! FUNCTION
!!  sort the list by c+ increasing with the quick sort algo
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

RECURSIVE SUBROUTINE ListCdagC_quickSort(this, begin, end)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  INTEGER,         INTENT(IN   ) :: begin
  INTEGER,         INTENT(IN   ) :: end
!Local variables k-----------------------------
  INTEGER                        :: it1
  INTEGER                        :: it2
  DOUBLE PRECISION               :: pivot
  DOUBLE PRECISION, DIMENSION(1:2):: CdagC_swap
  !DOUBLE PRECISION, DIMENSION(1:2):: ind_swap

  pivot = this%list((end-begin)/2 + begin,Cdag_) ! not the betterchoice.... FIXME
  it1 = begin
  it2 = end
  DO WHILE (it1 .LE. it2)
    DO WHILE ( this%list(it1,Cdag_) .LT. pivot )
      it1 = it1 + 1
    END DO
    DO WHILE ( this%list(it2,Cdag_) .GT. pivot )
      it2 = it2 - 1
    END DO
    IF ( it1 .LE. it2) THEN
      CdagC_swap = this%list(it1,1:2)
      !ind_swap   = this%ind (it1,1:2)
      this%list(it1,1:2) = this%list(it2,1:2)
      !this%ind (it1,1:2) = this%ind (it2,1:2)
      this%list(it2,1:2) = CdagC_swap
      !this%ind (it2,1:2) = ind_swap
      it1 = it1 + 1
      it2 = it2 - 1
    END IF
  END DO
  IF ( begin < it2 ) THEN
    CALL ListCdagC_quickSort(this,begin,it2)
  END IF
  !!it2= it1+1
  IF ( it1 < end ) THEN
    CALL ListCdagC_quickSort(this,it1,end)
  END IF

END SUBROUTINE ListCdagC_quickSort
!!***
 
!!****f* ABINIT/m_ListCdagC/ListCdagC_print
!! NAME
!!  ListCdagC_print
!!
!! FUNCTION
!!  print the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_print(this,ostream)

!Arguments ------------------------------------
  TYPE(ListCdagC)  , INTENT(IN) :: this
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: it

  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  WRITE(ostream_val,'(A,2x,A4,22x,A)') "#","Cdag", "C"
  DO it = 1, this%tail
    WRITE(ostream_val,*) this%list(it,Cdag_), this%list(it,C_) 
  END DO
END SUBROUTINE ListCdagC_print
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_clear
!! NAME
!!  ListCdagC_clear
!!
!! FUNCTION
!!  Clear the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_clear(this)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this
  this%tail = 0 
END SUBROUTINE ListCdagC_clear
!!***

!!****f* ABINIT/m_ListCdagC/ListCdagC_destroy
!! NAME
!!  ListCdagC_destroy
!!
!! FUNCTION
!!  destroy the list
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  list_1=ListCdagC
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

SUBROUTINE ListCdagC_destroy(this)

!Arguments ------------------------------------
  TYPE(ListCdagC), INTENT(INOUT) :: this

  FREEIF(this%list)
  !FAKEFREEIF(this%ind )

  this%tail     = 0
  this%size     = 0
  !this%inv_dt   = 0.d0
END SUBROUTINE ListCdagC_destroy
!!***

END MODULE m_ListCdagC
!!***


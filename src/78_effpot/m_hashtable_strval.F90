!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hastable
!! NAME
!! m_hashtable
!!
!! FUNCTION
!! This module provide a string: value pair hash table
!! COPYRIGHT
!! Taken from http://fortranwiki.org/fortran/show/hash+table+example
!! The code is originally written by Izaak Beekman under the LGPL license.
!! Adapted for usage in Abinit by hexu
!! Below is the original Copyright.
!!=======================================
!! Copyright (c) Izaak Beekman 2010
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_hashtable_strval

!!***
  use defs_basis
  use m_errors
  use m_abicore

  IMPLICIT NONE ! Use strong typing
  INTEGER, PARAMETER :: tbl_size = 50

  TYPE sllist
     TYPE(sllist), POINTER :: child => NULL()
     CHARACTER(len=:), ALLOCATABLE :: key
     real(dp) :: val
   CONTAINS
     PROCEDURE :: put  => put_sll
     PROCEDURE :: get  => get_sll
     PROCEDURE :: free => free_sll
  END TYPE sllist

  TYPE hash_table_t
     TYPE(sllist), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
   CONTAINS
     PROCEDURE :: init => init_hash_table_t
     PROCEDURE :: put  => put_hash_table_t
     PROCEDURE :: get  => get_hash_table_t
     PROCEDURE :: free => free_hash_table_t
  END TYPE hash_table_t

  PUBLIC :: hash_table_t
CONTAINS

  RECURSIVE SUBROUTINE put_sll(list,key,val)
    CLASS(sllist),    INTENT(inout) :: list
    CHARACTER(len=*), INTENT(in)    :: key
    real(dp), intent(in)  :: val
    INTEGER                         :: keylen

    keylen = LEN(key)
    IF (ALLOCATED(list%key)) THEN
       IF (list%key /= key) THEN
          IF ( .NOT. ASSOCIATED(list%child)) then
             !ABI_ALLOC_SCALAR(list%child)
             allocate(list%child)
             !call abimem_record(0, QUOTE(list%child), _LOC(list%child), "A", storage_size(scalar, kind=8),  __FILE__, __LINE__)
          end IF

          CALL put_sll(list%child,key,val)
       END IF
    ELSE
       IF (.NOT. ALLOCATED(list%key)) &
            ABI_DATATYPE_ALLOCATE_SCALAR(CHARACTER(len=keylen), list%key)
       list%key = key
       list%val = val
    END IF
  END SUBROUTINE put_sll


  RECURSIVE SUBROUTINE get_sll(list,key,val)
    CLASS(sllist),                 INTENT(in)    :: list
    CHARACTER(len=*),              INTENT(in)    :: key
    real(dp),                      INTENT(out)   :: val
    INTEGER                                      :: vallen

    vallen = 0
    IF (ALLOCATED(list%key) .AND. (list%key == key)) THEN
       val = list%val
    ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
       CALL get_sll(list%child,key,val)
    ELSE ! At the end of the list, no key found
       return
    END IF
  END SUBROUTINE get_sll


  RECURSIVE SUBROUTINE free_sll(list)
    CLASS(sllist), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll(list%child)
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    IF (ALLOCATED(list%key)) DEALLOCATE(list%key)
  END SUBROUTINE free_sll


  SUBROUTINE init_hash_table_t(tbl,tbl_len)
    CLASS(hash_table_t),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_hash_table_t

  ! The first part of the hashing procedure using the string
  ! collating sequence
  ELEMENTAL FUNCTION sum_string(str) RESULT(sig)
    CHARACTER(len=*), INTENT(in)   :: str
    INTEGER                        :: sig
    CHARACTER, DIMENSION(LEN(str)) :: tmp
    INTEGER :: i

    FORALL (i=1:LEN(str))
       tmp(i) = str(i:i)
    END FORALL
    sig = SUM(ICHAR(tmp))
  END FUNCTION sum_string


  SUBROUTINE put_hash_table_t(tbl,key,val)
    CLASS(hash_table_t), INTENT(inout) :: tbl
    CHARACTER(len=*),    INTENT(in)    :: key
    real(dp),            INTENT(in)    :: val
    INTEGER                            :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%put(key=key,val=val)
  END SUBROUTINE put_hash_table_t


  SUBROUTINE get_hash_table_t(tbl,key,val)
    CLASS(hash_table_t),           INTENT(in)    :: tbl
    CHARACTER(len=*),              INTENT(in)    :: key
    real(dp),                      INTENT(out)   :: val
    INTEGER                                      :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%get(key=key,val=val)
  END SUBROUTINE get_hash_table_t


  SUBROUTINE free_hash_table_t(tbl)
    CLASS(hash_table_t), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (ALLOCATED(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free()
       END DO
       DEALLOCATE(tbl%vec)
    END IF
    tbl%is_init = .FALSE.
  END SUBROUTINE free_hash_table_t


end module m_hashtable_strval

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_MatrixHyboffdiag
!! NAME
!!  m_MatrixHyboffdiag
!! 
!! FUNCTION 
!!  Module to deals with a matrix (mainly used for M matrix in m_BathOperatoroffdiag).
!!  Perform varius operation on matrices.
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

MODULE m_MatrixHyboffdiag

USE m_Global
IMPLICIT NONE

!!***

!!****t* m_MatrixHyboffdiag/MatrixHyboffdiag
!! NAME
!!  MatrixHyboffdiag
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

TYPE MatrixHyboffdiag
 
  INTEGER :: size
  ! size of the matrix

  INTEGER :: tail
  ! the size of the matrix that is actually used.

  INTEGER :: iTech = -1
  ! precise if time or frequency values will be used.

  INTEGER :: Wmax
  ! size if the frequency grid for mat_omega

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: mat
  ! matrix of size "size"

  INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: mat_tau
  ! matrix of size "size"

  COMPLEX(KIND=8) , ALLOCATABLE, DIMENSION(:,:,:) :: mat_omega
  ! array of size (Wmax, size, size)

END TYPE MatrixHyboffdiag
!!***

INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE MatrixHyboffdiag_assign
END INTERFACE

CONTAINS
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_init
!! NAME
!!  MatrixHyboffdiag_init
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
!!  matrix_1=matrix
!!  iTech=DO NOT USE => BUG
!!  size=memory size for initialization
!!  Wmax=maximum frequency number
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

SUBROUTINE MatrixHyboffdiag_init(matrix_1, iTech, size, Wmax)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_init'
!End of the abilint section

  TYPE(MatrixHyboffdiag)     , INTENT(INOUT) :: matrix_1
  INTEGER             , INTENT(IN   ) :: iTech
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
  INTEGER, OPTIONAL, INTENT(IN   ) :: Wmax
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  matrix_1%size = size_val
  FREEIF(matrix_1%mat)
  MALLOC(matrix_1%mat,(1:size_val,1:size_val))
  matrix_1%tail  = 0 
  matrix_1%mat   = 0.d0
  matrix_1%iTech = iTech
  SELECT CASE(matrix_1%iTech)
  CASE (GREENHYB_TAU)
    FREEIF(matrix_1%mat_tau)
    MALLOC(matrix_1%mat_tau,(1:size_val,1:size_val))
    matrix_1%mat_tau=0
  CASE (GREENHYB_OMEGA)
    IF ( PRESENT(Wmax) .AND. Wmax .GT. 0 ) THEN
      FREEIF(matrix_1%mat_omega)
      MALLOC(matrix_1%mat_omega,(1:Wmax,1:size_val,1:size_val))
      matrix_1%Wmax=Wmax
    ELSE
      CALL ERROR("MatrixHyboffdiag_init : Missing argument Wmax for measurement in omega")
    END IF
  CASE DEFAULT
      CALL WARNALL("MatrixHyboffdiag_init : Wrong input argument iTech")
  END SELECT
END SUBROUTINE MatrixHyboffdiag_init
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_setSize
!! NAME
!!  MatrixHyboffdiag_setSize
!!
!! FUNCTION
!!  impose size of the matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=matrix
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

SUBROUTINE MatrixHyboffdiag_setSize(matrix_1,new_tail)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_setSize'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(INOUT) :: matrix_1
  INTEGER        , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                        :: size

  IF ( .NOT. ALLOCATED(matrix_1%mat) ) &
    CALL MatrixHyboffdiag_init(matrix_1,matrix_1%iTech,Wmax=matrix_1%Wmax)
  size = matrix_1%size
  IF( new_tail .GT. size ) CALL MatrixHyboffdiag_enlarge(matrix_1, MAX(Global_SIZE,new_tail-size))
  matrix_1%tail = new_tail
END SUBROUTINE MatrixHyboffdiag_setSize  
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_enlarge
!! NAME
!!  MatrixHyboffdiag_enlarge
!!
!! FUNCTION
!!  This subroutine enlarges memory space
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=matrix
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

SUBROUTINE MatrixHyboffdiag_enlarge(matrix_1, size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_enlarge'
!End of the abilint section

  TYPE(MatrixHyboffdiag)     , INTENT(INOUT)          :: matrix_1
  INTEGER, OPTIONAL, INTENT(IN   )          :: size
!Local variables ------------------------------
  INTEGER                                   :: width
  INTEGER                                   :: tail
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: matrix_temp 
  INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: matrix_temp_tau
  COMPLEX(KIND=8) , ALLOCATABLE, DIMENSION(:,:,:) :: matrix_temp_omega
  INTEGER                                   :: size_val

  IF ( ALLOCATED(matrix_1%mat) ) THEN
    FREEIF(matrix_temp)
    FREEIF(matrix_temp_tau)
    width = matrix_1%size
    tail  = matrix_1%tail 
    size_val = width
    IF ( PRESENT(size) ) size_val = size 

!   change size of mat
    MALLOC(matrix_temp,(1:tail,1:tail))
    matrix_temp(1:tail,1:tail) = matrix_1%mat(1:tail,1:tail)
    FREE(matrix_1%mat)
    matrix_1%size = width + size_val
    MALLOC(matrix_1%mat,(1:matrix_1%size,1:matrix_1%size))
    matrix_1%mat(1:tail,1:tail) = matrix_temp(1:tail,1:tail)
    FREE(matrix_temp)

    SELECT CASE(matrix_1%iTech)
    CASE (GREENHYB_TAU)

!   change size of mat_tau
      MALLOC(matrix_temp_tau,(1:tail,1:tail))
      matrix_temp_tau(1:tail,1:tail) = matrix_1%mat_tau(1:tail,1:tail)
      FREE(matrix_1%mat_tau)
      MALLOC(matrix_1%mat_tau,(1:matrix_1%size,1:matrix_1%size))
      matrix_1%mat_tau(1:tail,1:tail) = matrix_temp_tau(1:tail,1:tail)
      FREE(matrix_temp_tau)

    CASE (GREENHYB_OMEGA)

!   change size of mat_omega
      MALLOC(matrix_temp_omega,(1:matrix_1%Wmax,1:tail,1:tail))
      matrix_temp_omega(1:matrix_1%Wmax,1:tail,1:tail) = matrix_1%mat_omega(1:matrix_1%Wmax,1:tail,1:tail)
      FREE(matrix_1%mat_omega)
      MALLOC(matrix_1%mat_omega,(1:matrix_1%Wmax,1:matrix_1%size,1:matrix_1%size))
      matrix_1%mat_omega(1:matrix_1%Wmax,1:tail,1:tail) = matrix_temp_omega(1:matrix_1%Wmax,1:tail,1:tail)
      FREE(matrix_temp_omega)

    END SELECT

  ELSE
    CALL MatrixHyboffdiag_init(matrix_1, matrix_1%iTech, size=Global_SIZE, Wmax=matrix_1%Wmax)
  END IF
END SUBROUTINE MatrixHyboffdiag_enlarge
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_clear
!! NAME
!!  MatrixHyboffdiag_clear
!!
!! FUNCTION
!!  Clear matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=matrix
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

SUBROUTINE MatrixHyboffdiag_clear(matrix_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_clear'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(INOUT) :: matrix_1
  matrix_1%tail = 0 
END SUBROUTINE MatrixHyboffdiag_clear
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_assign
!! NAME
!!  MatrixHyboffdiag_assign
!!
!! FUNCTION
!!  assign: copy matrix_2 in matrix_1
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=Matrix 1
!!  matrix_2=Matrix 2
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

SUBROUTINE MatrixHyboffdiag_assign(matrix_1, matrix_2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_assign'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(INOUT) :: matrix_1
  TYPE(MatrixHyboffdiag), INTENT(IN   ) :: matrix_2
!Local variables ------------------------------
  INTEGER                        :: tail

  tail = matrix_2%tail
  CALL MatrixHyboffdiag_setSize(matrix_1, tail)
  matrix_1%mat(1:tail,1:tail) = matrix_2%mat(1:tail,1:tail)
  IF ( matrix_1%iTech .NE. matrix_2%iTech ) & 
    CALL ERROR("MatrixHyboffdiag_assign : not compatible matrices")
  SELECT CASE(matrix_1%iTech)
  CASE (GREENHYB_TAU)
    matrix_1%mat_tau(1:tail,1:tail) = matrix_2%mat_tau(1:tail,1:tail)
  CASE (GREENHYB_OMEGA)
    IF ( matrix_1%Wmax .NE. matrix_2%Wmax ) &
      CALL ERROR("MatrixHyboffdiag_assig : different numbers of omega")
    matrix_1%mat_omega(1:matrix_1%Wmax,1:tail,1:tail) = &
    matrix_2%mat_omega(1:matrix_1%Wmax,1:tail,1:tail)
  END SELECT

END SUBROUTINE MatrixHyboffdiag_assign
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_inverse
!! NAME
!!  MatrixHyboffdiag_inverse
!!
!! FUNCTION
!!  inverse the  matrix and compute the determinant
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=matrix
!!
!! OUTPUT
!!  determinant=determinant of the matrix before inversion
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

SUBROUTINE MatrixHyboffdiag_inverse(matrix_1,determinant)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_inverse'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(INOUT) :: matrix_1
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: determinant
!Local variables ------------------------------
  !DOUBLE PRECISION :: reste
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: invMatrix
  INTEGER :: ligne
  INTEGER :: ligne_virtuelle
  INTEGER :: colonne
  INTEGER :: sys
  INTEGER :: sys_vir
  INTEGER :: tail
  INTEGER, DIMENSION(:), ALLOCATABLE :: pivot 
  DOUBLE PRECISION :: det

  tail = matrix_1%tail
  IF ( tail .EQ. 0 ) THEN
    IF ( PRESENT(determinant) ) determinant = 1.d0
    RETURN
  END IF

  MALLOC(invMatrix,(1:tail,1:tail))

  FREEIF(pivot)
  MALLOC(pivot,(1:tail))

  CALL MatrixHyboffdiag_LU(matrix_1, pivot=pivot, determinant=det)

  !det = 1.d0
  DO sys = 1, tail
    sys_vir = pivot(sys)
!!    DO ligne=1,tail
!!      ligne_virtuelle = pivot(ligne)
!!      reste = 0.d0
!!      DO colonne = 1, ligne-1
!!        reste = reste + matrix_1%mat(ligne_virtuelle,colonne)*invMatrix%mat(colonne,sys_vir)
!!      END DO
!!      IF ( ligne .EQ. sys ) THEN
!!        invMatrix%mat(ligne,sys_vir) = 1.d0 - reste
!!      ELSE
!!        invMatrix%mat(ligne,sys_vir) = 0.d0 - reste
!!      END IF
!!    END DO
!!    
!!    DO ligne = tail, 1, -1
!!      ligne_virtuelle=pivot(ligne)
!!      reste = 0.d0
!!        DO colonne = ligne+1, tail
!!          reste= reste + matrix_1%mat(ligne_virtuelle,colonne)*invMatrix%mat(colonne,sys_vir)
!!        END DO
!!        invMatrix%mat(ligne, sys_vir) = (invMatrix%mat(ligne, sys_vir)-reste)/matrix_1%mat(ligne_virtuelle,ligne)
!!    END DO
    !det = det*matrix_1%mat(sys_vir,sys)
    invMatrix(:,sys_vir) = 0.d0
    invMatrix(sys,sys_vir) = 1.d0
    DO colonne = 1, tail-1
      DO ligne=colonne+1,tail
        ligne_virtuelle = pivot(ligne)
        invMatrix(ligne, sys_vir) = invMatrix(ligne,sys_vir) &
                                                - matrix_1%mat(ligne_virtuelle,colonne)  &
                                                * invMatrix(colonne,sys_vir)
      END DO
    END DO

    DO colonne = tail, 1, -1
      invMatrix(colonne, sys_vir) = invMatrix(colonne, sys_vir) / matrix_1%mat(pivot(colonne),colonne)
      DO ligne = 1, colonne-1
        ligne_virtuelle = pivot(ligne)
        invMatrix(ligne, sys_vir) = invMatrix(ligne,sys_vir) &
                                                 - matrix_1%mat(ligne_virtuelle,colonne) &
                                                 * invMatrix(colonne, sys_vir)
      END DO
    END DO

  END DO
  FREE(pivot)
  matrix_1%mat(1:tail,1:tail) = invMatrix(1:tail,1:tail)
  IF ( PRESENT(determinant) ) THEN
    determinant = det
  END IF
  FREE(invMatrix)
END SUBROUTINE MatrixHyboffdiag_inverse
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_LU
!! NAME
!!  MatrixHyboffdiag_LU
!!
!! FUNCTION
!!  LU decomposition
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mat_a=matrix
!!
!! OUTPUT
!!  pivot=gauss pivot
!!  determinant=determinant of the matrix
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

SUBROUTINE MatrixHyboffdiag_LU(mat_a,pivot,determinant)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_LU'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(INOUT) :: mat_a
  INTEGER, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: pivot
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: determinant
!Local variables ------------------------------
  INTEGER :: ligne
  INTEGER :: colonne 
  INTEGER :: colonne_de_1_ligne
  INTEGER :: max_ind_lig
  INTEGER :: ligne_virtuelle
  INTEGER :: tail
  INTEGER, DIMENSION(:), ALLOCATABLE :: pivot_tmp
  DOUBLE PRECISION :: max_col
  DOUBLE PRECISION :: inverse_pivot
  DOUBLE PRECISION :: coef
  DOUBLE PRECISION :: det
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mat_tmp

  tail = mat_a%tail
  det = 1.d0
  FREEIF(pivot_tmp)
  MALLOC(pivot_tmp,(1:tail))

  DO ligne=1, tail
    pivot_tmp(ligne)=ligne
  END DO

  MALLOC(mat_tmp,(1:tail,1:tail))
  mat_tmp(1:tail,1:tail) = mat_a%mat(1:tail,1:tail)
  DO colonne = 1, tail-1 
    max_col = ABS(mat_tmp(pivot_tmp(colonne),colonne))
    max_ind_lig = colonne

    DO ligne = colonne+1, tail 
      ligne_virtuelle = pivot_tmp(ligne)
      IF ( ABS( mat_tmp(ligne_virtuelle,colonne)).GT.max_col) then
        max_col = ABS(mat_tmp(ligne_virtuelle,colonne))
        max_ind_lig = ligne
      ENDIF
    END DO

    ligne              = pivot_tmp(colonne) 
    pivot_tmp(colonne)     = pivot_tmp(max_ind_lig)
    pivot_tmp(max_ind_lig) = ligne
    IF ( pivot_tmp(colonne) .NE. pivot_tmp(max_ind_lig) ) det = det * (-1.d0)
    inverse_pivot=1.d0 / mat_tmp(pivot_tmp(colonne),colonne)
    det = det * mat_tmp(pivot_tmp(colonne),colonne)
    DO ligne = colonne+1, tail 
      ligne_virtuelle = pivot_tmp(ligne)
      coef = mat_tmp(ligne_virtuelle,colonne)*inverse_pivot
      mat_tmp(ligne_virtuelle,colonne) = coef
      DO colonne_de_1_ligne = colonne+1, tail 
        mat_tmp(ligne_virtuelle,colonne_de_1_ligne)= mat_tmp(ligne_virtuelle,colonne_de_1_ligne)&
                                         -coef * mat_tmp(pivot_tmp(colonne) ,colonne_de_1_ligne)
      END DO
    END DO
  END DO
  det = det * mat_tmp(pivot_tmp(tail),tail)
  IF ( PRESENT(determinant) ) &
    determinant = det
  IF ( PRESENT(pivot) ) THEN
    mat_a%mat(1:tail,1:tail) = mat_tmp(1:tail,1:tail)
    IF ( ALLOCATED(pivot) .AND. SIZE(pivot) .NE. tail ) THEN
      FREE(pivot)
      MALLOC(pivot,(1:tail))
    ELSE IF ( .NOT. ALLOCATED(pivot) ) THEN
      MALLOC(pivot,(1:tail))
    END IF
    pivot(1:tail)=pivot_tmp(1:tail)
  END IF
  FREE(mat_tmp)
  FREE(pivot_tmp)
END SUBROUTINE MatrixHyboffdiag_LU
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_getDet
!! NAME
!!  MatrixHyboffdiag_getDet
!!
!! FUNCTION
!!  Just get the determinant 
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_a=matrix
!!
!! OUTPUT
!!  det=determinant
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

SUBROUTINE MatrixHyboffdiag_getDet(matrix_a,det)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_getDet'
!End of the abilint section

  TYPE(MatrixHyboffdiag) , INTENT(INOUT) :: matrix_a
  DOUBLE PRECISION, INTENT(  OUT) :: det

  IF ( matrix_a%tail .EQ. 0 ) THEN
    det = 1.d0
    RETURN
  END IF
  CALL MatrixHyboffdiag_LU(matrix_a, determinant=det)
END SUBROUTINE MatrixHyboffdiag_getDet
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_print
!! NAME
!!  MatrixHyboffdiag_print
!!
!! FUNCTION
!!  print Matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=Matrix
!!  ostream=file stream
!!  opt_print=0 mat
!!            1 mat_tau
!!            2 both
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

SUBROUTINE MatrixHyboffdiag_print(matrix_1,ostream,opt_print)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_print'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(IN) :: matrix_1
  INTEGER, OPTIONAL, INTENT(IN) :: ostream
  INTEGER, OPTIONAL, INTENT(IN) :: opt_print
!Local variables ------------------------------
  INTEGER                       :: ostream_val
  INTEGER                       :: opt_val
  INTEGER                       :: it1
  INTEGER                       :: it2
  CHARACTER(LEN=4 )             :: size
  CHARACTER(LEN=22)             :: string

  ostream_val = 6
  opt_val=0
  IF ( PRESENT(ostream) ) ostream_val = ostream
  IF ( PRESENT(opt_print) ) opt_val = opt_print
  WRITE(size,'(I4)') matrix_1%tail
  IF ( MOD(opt_val,2) .EQ. 0 ) THEN
    string ='(1x,A,1x,'//TRIM(ADJUSTL(size))//'(ES10.2,1x),A)'
    WRITE(ostream_val,'(A)') "["
    DO it1 = 1, matrix_1%tail
      WRITE(ostream_val,string) "[",(/ (matrix_1%mat(it1,it2),it2=1,matrix_1%tail)  /)," ]"
    END DO
    WRITE(ostream_val,'(A)') "]"
  END IF  
  IF ( opt_val .GE.1 ) THEN
    string ='(1x,A,1x,'//TRIM(ADJUSTL(size))//'(I4,1x),A)'
    WRITE(ostream_val,'(A)') "["
    DO it1 = 1, matrix_1%tail
      WRITE(ostream_val,string) "[",(/ (matrix_1%mat_tau(it1,it2),it2=1,matrix_1%tail)  /),"]"
    END DO
    WRITE(ostream_val,'(A)') "]"
  END IF
END SUBROUTINE MatrixHyboffdiag_print
!!***

!!****f* ABINIT/m_MatrixHyboffdiag/MatrixHyboffdiag_destroy
!! NAME
!!  MatrixHyboffdiag_destroy
!!
!! FUNCTION
!!  Destroy
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matrix_1=Matrix
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

SUBROUTINE MatrixHyboffdiag_destroy(matrix_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatrixHyboffdiag_destroy'
!End of the abilint section

  TYPE(MatrixHyboffdiag), INTENT(INOUT) :: matrix_1

  FREEIF(matrix_1%mat)
  SELECT CASE(matrix_1%iTech)
  CASE (GREENHYB_TAU)
    FREEIF(matrix_1%mat_tau)
  CASE (GREENHYB_OMEGA)
    FREEIF(matrix_1%mat_omega)
  END SELECT

  matrix_1%tail     = 0
  matrix_1%size     = 0
  matrix_1%iTech    = -1
END SUBROUTINE MatrixHyboffdiag_destroy
!!***

END MODULE m_MatrixHyboffdiag
!!***


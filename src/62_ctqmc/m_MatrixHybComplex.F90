
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_MatrixHybComplex
!! NAME
!!  m_MatrixHybComplex
!! 
!! FUNCTION 
!!  Module to deals with a matrix (mainly used for M matrix in m_BathOperatoroffdiagComplex).
!!  Perform varius operation on matrices.
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
MODULE m_MatrixHybComplex
USE m_Global
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_MatrixHybComplex/MatrixHybComplex
!! NAME
!!  MatrixHybComplex
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

TYPE, PUBLIC :: MatrixHybComplex
  INTEGER _PRIVATE :: size
  ! size of the matrix

  INTEGER          :: tail
  ! the size of the matrix that is actually used.

  INTEGER _PRIVATE :: iTech = -1
  ! precise if time or frequency values will be used.

  INTEGER _PRIVATE :: Wmax
  ! size if the frequency grid for mat_omega

  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:)           :: mat
  ! matrix of size "size"

  INTEGER         , ALLOCATABLE, DIMENSION(:,:)           :: mat_tau
  ! matrix of size "size"

  COMPLEX(KIND=8) , ALLOCATABLE, DIMENSION(:,:,:)         :: mat_omega
  ! array of size (Wmax, size, size)

END TYPE MatrixHybComplex
!!***

!PUBLIC INTERFACE ASSIGNMENT (=)
!  MODULE PROCEDURE MatrixHybComplex_assign
!END INTERFACE

PUBLIC  :: MatrixHybComplex_init
PUBLIC  :: MatrixHybComplex_setSize
PRIVATE :: MatrixHybComplex_enlarge
PUBLIC  :: MatrixHybComplex_clear
PUBLIC  :: MatrixHybComplex_assign
PUBLIC  :: MatrixHybComplex_inverse
PRIVATE :: MatrixHybComplex_LU
PUBLIC  :: MatrixHybComplex_getDet
PUBLIC  :: MatrixHybComplex_print
PUBLIC  :: MatrixHybComplex_destroy

CONTAINS
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_init
!! NAME
!!  MatrixHybComplex_init
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
!!  this=matrix
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
!! SOURCE

SUBROUTINE MatrixHybComplex_init(this, iTech, size, Wmax)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex)     , INTENT(INOUT) :: this
  INTEGER             , INTENT(IN   ) :: iTech
  INTEGER, OPTIONAL, INTENT(IN   ) :: size
  INTEGER, OPTIONAL, INTENT(IN   ) :: Wmax
!Local variables ------------------------------
  INTEGER                          :: size_val

  size_val = Global_SIZE
  IF ( PRESENT(size) ) size_val = size
  this%size = size_val
  FREEIF(this%mat)
  MALLOC(this%mat,(1:size_val,1:size_val))
  this%tail  = 0 
  this%mat   = 0.d0
  this%iTech = iTech
  SELECT CASE(this%iTech)
  CASE (GREENHYB_TAU)
    FREEIF(this%mat_tau)
    MALLOC(this%mat_tau,(1:size_val,1:size_val))
  CASE (GREENHYB_OMEGA)
    IF ( PRESENT(Wmax) .AND. Wmax .GT. 0 ) THEN
      FREEIF(this%mat_omega)
      MALLOC(this%mat_omega,(1:Wmax,1:size_val,1:size_val))
      this%Wmax=Wmax
    ELSE
      CALL ERROR("MatrixHybComplex_init : Missing argument Wmax for measurement in omega")
    END IF
  CASE DEFAULT
      CALL WARNALL("MatrixHybComplex_init : Wrong input argument iTech")
  END SELECT
END SUBROUTINE MatrixHybComplex_init
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_setSize
!! NAME
!!  MatrixHybComplex_setSize
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
!!  this=matrix
!!  new_tail=new size
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_setSize(this,new_tail)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(INOUT) :: this
  INTEGER        , INTENT(IN   ) :: new_tail
!Local variables ------------------------------
  INTEGER                        :: size

  IF ( .NOT. ALLOCATED(this%mat) ) &
    CALL MatrixHybComplex_init(this,this%iTech,Wmax=this%Wmax)
  size = this%size
  IF( new_tail .GT. size ) THEN
    CALL MatrixHybComplex_enlarge(this, MAX(Global_SIZE,new_tail-size))
  END IF
  this%tail = new_tail
END SUBROUTINE MatrixHybComplex_setSize  
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_enlarge
!! NAME
!!  MatrixHybComplex_enlarge
!!
!! FUNCTION
!!  This subroutine enlarges memory space
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=matrix
!!  size=new memory size
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_enlarge(this, size)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex)     , INTENT(INOUT)          :: this
  INTEGER, OPTIONAL, INTENT(IN   )          :: size
!Local variables ------------------------------
  INTEGER                                   :: width
  INTEGER                                   :: tail
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: this_temp 
  INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: this_temp_tau
  COMPLEX(KIND=8) , ALLOCATABLE, DIMENSION(:,:,:) :: this_temp_omega
  INTEGER                                   :: size_val

  IF ( ALLOCATED(this%mat) ) THEN
    FREEIF(this_temp)
    FREEIF(this_temp_tau)
    width = this%size
    tail  = this%tail 
    size_val = width
    IF ( PRESENT(size) ) size_val = size 

!   change size of mat
    MALLOC(this_temp,(1:tail,1:tail))
    this_temp(1:tail,1:tail) = this%mat(1:tail,1:tail)
    FREE(this%mat)
    this%size = width + size_val
    MALLOC(this%mat,(1:this%size,1:this%size))
    this%mat(1:tail,1:tail) = this_temp(1:tail,1:tail)
    FREE(this_temp)
    SELECT CASE(this%iTech)
    CASE (GREENHYB_TAU)

!   change size of mat_tau
      MALLOC(this_temp_tau,(1:tail,1:tail))
      this_temp_tau(1:tail,1:tail) = this%mat_tau(1:tail,1:tail)
      FREE(this%mat_tau)
      MALLOC(this%mat_tau,(1:this%size,1:this%size))
      this%mat_tau(1:tail,1:tail) = this_temp_tau(1:tail,1:tail)
      FREE(this_temp_tau)
    CASE (GREENHYB_OMEGA)

!   change size of mat_omega
      MALLOC(this_temp_omega,(1:this%Wmax,1:tail,1:tail))
      this_temp_omega(1:this%Wmax,1:tail,1:tail) = this%mat_omega(1:this%Wmax,1:tail,1:tail)
      FREE(this%mat_omega)
      MALLOC(this%mat_omega,(1:this%Wmax,1:this%size,1:this%size))
      this%mat_omega(1:this%Wmax,1:tail,1:tail) = this_temp_omega(1:this%Wmax,1:tail,1:tail)
      FREE(this_temp_omega)
    END SELECT
  ELSE
    CALL MatrixHybComplex_init(this, this%iTech, size=Global_SIZE, Wmax=this%Wmax)
  END IF
END SUBROUTINE MatrixHybComplex_enlarge
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_clear
!! NAME
!!  MatrixHybComplex_clear
!!
!! FUNCTION
!!  Clear this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=matrix
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_clear(this)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(INOUT) :: this
  this%tail = 0 
END SUBROUTINE MatrixHybComplex_clear
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_assign
!! NAME
!!  MatrixHybComplex_assign
!!
!! FUNCTION
!!  assign this=matrix2
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Matrix 1
!!  Matrix2=Matrix 2
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_assign(this, matrix)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(INOUT) :: this
  TYPE(MatrixHybComplex), INTENT(IN   ) :: matrix
!Local variables ------------------------------
  INTEGER                        :: tail

  tail = matrix%tail
  CALL MatrixHybComplex_setSize(this, tail)
  this%mat(1:tail,1:tail) = matrix%mat(1:tail,1:tail)
  IF ( this%iTech .NE. matrix%iTech ) & 
    CALL ERROR("MatrixHybComplex_assign : not compatible matrices")
  SELECT CASE(this%iTech)
  CASE (GREENHYB_TAU)
    this%mat_tau(1:tail,1:tail) = matrix%mat_tau(1:tail,1:tail)
  CASE (GREENHYB_OMEGA)
    IF ( this%Wmax .NE. matrix%Wmax ) &
      CALL ERROR("MatrixHybComplex_assig : different numbers of omega")
    this%mat_omega(1:this%Wmax,1:tail,1:tail) = &
    matrix%mat_omega(1:this%Wmax,1:tail,1:tail)
  END SELECT

END SUBROUTINE MatrixHybComplex_assign
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_inverse
!! NAME
!!  MatrixHybComplex_inverse
!!
!! FUNCTION
!!  inverse the matrix and compute the determinant
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=this
!!
!! OUTPUT
!!  determinant=determinant of the this before inversion
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_inverse(this,determinant)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(INOUT) :: this
  COMPLEX(KIND=8), OPTIONAL, INTENT(OUT) :: determinant
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: invMatrix
!Local variables ------------------------------
  !DOUBLE PRECISION :: reste
  INTEGER :: ligne
  INTEGER :: ligne_virtuelle
  INTEGER :: colonne
  INTEGER :: sys
  INTEGER :: sys_vir
  INTEGER :: tail
  INTEGER, DIMENSION(:), ALLOCATABLE :: pivot 
  COMPLEX(KIND=8) :: det

  tail = this%tail
  IF ( tail .EQ. 0 ) THEN
    IF ( PRESENT(determinant) ) determinant = cmplx(1.d0,0.d0,kind=8)
    RETURN
  END IF

  MALLOC(invMatrix,(1:tail,1:tail))

  FREEIF(pivot)
  MALLOC(pivot,(1:tail))

  CALL MatrixHybComplex_LU(this, pivot=pivot, determinant=det)

  !det = 1.d0
  DO sys = 1, tail
    sys_vir = pivot(sys)
!!    DO ligne=1,tail
!!      ligne_virtuelle = pivot(ligne)
!!      reste = 0.d0
!!      DO colonne = 1, ligne-1
!!        reste = reste + this%mat(ligne_virtuelle,colonne)*invMatrix%mat(colonne,sys_vir)
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
!!          reste= reste + this%mat(ligne_virtuelle,colonne)*invMatrix%mat(colonne,sys_vir)
!!        END DO
!!        invMatrix%mat(ligne, sys_vir) = (invMatrix%mat(ligne, sys_vir)-reste)/this%mat(ligne_virtuelle,ligne)
!!    END DO
    !det = det*this%mat(sys_vir,sys)
    invMatrix(:,sys_vir) = cmplx(0.d0,0.d0,kind=8)
    invMatrix(sys,sys_vir) = cmplx(1.d0,0.d0,kind=8)
    DO colonne = 1, tail-1
      DO ligne=colonne+1,tail
        ligne_virtuelle = pivot(ligne)
        invMatrix(ligne, sys_vir) = invMatrix(ligne,sys_vir) &
                                                - this%mat(ligne_virtuelle,colonne)  &
                                                * invMatrix(colonne,sys_vir)
      END DO
    END DO

    DO colonne = tail, 1, -1
      invMatrix(colonne, sys_vir) = invMatrix(colonne, sys_vir) / this%mat(pivot(colonne),colonne)
      DO ligne = 1, colonne-1
        ligne_virtuelle = pivot(ligne)
        invMatrix(ligne, sys_vir) = invMatrix(ligne,sys_vir) &
                                                 - this%mat(ligne_virtuelle,colonne) &
                                                 * invMatrix(colonne, sys_vir)
      END DO
    END DO

  END DO
  FREE(pivot)
  this%mat(1:tail,1:tail) = invMatrix(1:tail,1:tail)
  IF ( PRESENT(determinant) ) THEN
    determinant = det
  END IF
  FREE(invMatrix)
END SUBROUTINE MatrixHybComplex_inverse
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_LU
!! NAME
!!  MatrixHybComplex_LU
!!
!! FUNCTION
!!  LU decomposition
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=matrix
!!
!! OUTPUT
!!  pivot=gauss pivot
!!  determinant=determinant of the this
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_LU(this,pivot,determinant)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(INOUT) :: this
  INTEGER, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: pivot
  COMPLEX(KIND=8), OPTIONAL, INTENT(OUT) :: determinant
!Local variables ------------------------------
  INTEGER :: ligne
  INTEGER :: colonne 
  INTEGER :: colonne_de_1_ligne
  INTEGER :: max_ind_lig
  INTEGER :: ligne_virtuelle
  INTEGER :: tail
  INTEGER, DIMENSION(:), ALLOCATABLE :: pivot_tmp
  DOUBLE PRECISION :: max_col
  COMPLEX(KIND=8) :: inverse_pivot
  COMPLEX(KIND=8) :: coef
  COMPLEX(KIND=8) :: det
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: mat_tmp

  tail = this%tail
  det = cmplx(1.d0,0.d0,kind=8)
  FREEIF(pivot_tmp)
  MALLOC(pivot_tmp,(1:tail))

  DO ligne=1, tail
    pivot_tmp(ligne)=ligne
  END DO

  MALLOC(mat_tmp,(1:tail,1:tail))
  mat_tmp(1:tail,1:tail) = this%mat(1:tail,1:tail)
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
    this%mat(1:tail,1:tail) = mat_tmp(1:tail,1:tail)
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
END SUBROUTINE MatrixHybComplex_LU
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_getDet
!! NAME
!!  MatrixHybComplex_getDet
!!
!! FUNCTION
!!  Just get the determinant 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=matrix
!!
!! OUTPUT
!!  det=determinant
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_getDet(this,det)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex) , INTENT(INOUT) :: this
  COMPLEX(KIND=8), INTENT(  OUT) :: det

  IF ( this%tail .EQ. 0 ) THEN
    det = cmplx(1.d0,0.d0,kind=8)
    RETURN
  END IF
  CALL MatrixHybComplex_LU(this, determinant=det)
END SUBROUTINE MatrixHybComplex_getDet
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_print
!! NAME
!!  MatrixHybComplex_print
!!
!! FUNCTION
!!  print Matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Matrix
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
!! SOURCE

SUBROUTINE MatrixHybComplex_print(this,ostream,opt_print)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(IN) :: this
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
  WRITE(size,'(I4)') this%tail
  IF ( MOD(opt_val,2) .EQ. 0 ) THEN
    string ='(1x,A,1x,'//TRIM(ADJUSTL(size))//'(ES10.2,1x),A)'
    WRITE(ostream_val,'(A)') "["
    DO it1 = 1, this%tail
      WRITE(ostream_val,string) "[",(/ (this%mat(it1,it2),it2=1,this%tail)  /)," ]"
    END DO
    WRITE(ostream_val,'(A)') "]"
  END IF  
  IF ( opt_val .GE.1 ) THEN
    string ='(1x,A,1x,'//TRIM(ADJUSTL(size))//'(I4,1x),A)'
    WRITE(ostream_val,'(A)') "["
    DO it1 = 1, this%tail
      WRITE(ostream_val,string) "[",(/ (this%mat_tau(it1,it2),it2=1,this%tail)  /),"]"
    END DO
    WRITE(ostream_val,'(A)') "]"
  END IF
END SUBROUTINE MatrixHybComplex_print
!!***

!!****f* ABINIT/m_MatrixHybComplex/MatrixHybComplex_destroy
!! NAME
!!  MatrixHybComplex_destroy
!!
!! FUNCTION
!!  Destroy
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Matrix
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

SUBROUTINE MatrixHybComplex_destroy(this)

!Arguments ------------------------------------
  TYPE(MatrixHybComplex), INTENT(INOUT) :: this

  FREEIF(this%mat)
  SELECT CASE(this%iTech)
  CASE (GREENHYB_TAU)
    FREEIF(this%mat_tau)
  CASE (GREENHYB_OMEGA)
    FREEIF(this%mat_omega)
  END SELECT

  this%tail     = 0
  this%size     = 0
  this%iTech    = -1
END SUBROUTINE MatrixHybComplex_destroy
!!***

END MODULE m_MatrixHybComplex
!!***


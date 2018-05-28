#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ImpurityOperatoroffdiag
!! NAME
!!  m_ImpurityOperatoroffdiag
!! 
!! FUNCTION 
!!  manage all related to Impurity
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
MODULE m_ImpurityOperatoroffdiag
USE m_ListCdagCoffdiag
USE m_Global
IMPLICIT NONE

!!***

!!****t* m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag
!! NAME
!!  ImpurityOperatoroffdiag
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

TYPE ImpurityOperatoroffdiag
  LOGICAL                   :: doCheck = .FALSE.
  INTEGER                   :: flavors
   !  Number of flavors
  INTEGER                   :: activeFlavor
   !  Flavor considered e.g when a segment is added

!  TYPE(ListCdagCoffdiag) , POINTER :: activeParticle => NULL()
  DOUBLE PRECISION          :: beta
   !  Inverse of temperature.

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: mat_U 
   !  for iflavor1 and iflavor2, mat_U(iflavor1,iflavor2) is the
   !  coulomb interaction between iflavor1 and iflavor2.

 !DOUBLE PRECISION, POINTER, DIMENSION(:  ) :: mu    => NULL()
 !DOUBLE PRECISION                          :: shift_mu   

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: overlaps   
   !  for iflavor1 and iflavor2 overlaps(iflavor1,iflavor2) is the total
   !  overlap between segments of iflavor1 and segments of iflavor2.

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:  ) :: updates    ! new_(anti)seg 
   !  For a given flavor (activeflavor), gives for each other flavors, the
   !  supplementary overlaps, called updates(otherflavor).

  TYPE(ListCdagCoffdiag) , ALLOCATABLE, DIMENSION(:  ) :: particles 
   !  for each flavor, particles(iflavor)%list(2,maxnbofsegment) 
   !  gives the beginning and end of each segment.


!#ifdef CTQMC_CHECK
  DOUBLE PRECISION                          :: checkNumber
  DOUBLE PRECISION                          :: tolerance
  DOUBLE PRECISION                          :: meanError
!#endif
END TYPE ImpurityOperatoroffdiag
!!***

INTERFACE ImpurityOperatoroffdiag_setU
  MODULE PROCEDURE ImpurityOperatoroffdiag_computeU, ImpurityOperatoroffdiag_setUmat
END INTERFACE

CONTAINS
!!***

!SUBROUTINE ImpurityOperatoroffdiag_init(op, flavors, beta, N)
!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_init
!! NAME
!!  ImpurityOperatoroffdiag_init
!!
!! FUNCTION
!!  Initialize and allocate
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurtiyOperator
!!  flavors=number of flavors
!!  beta=inverse temperature
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

SUBROUTINE ImpurityOperatoroffdiag_init(op, flavors, beta)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_init'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER               , INTENT(IN   ) :: flavors
  !DOUBLE PRECISION      , INTENT(IN   ) :: U
  !DOUBLE PRECISION      , INTENT(IN   ) :: J
  DOUBLE PRECISION      , INTENT(IN   ) :: beta
  !INTEGER               , INTENT(IN   ) :: N
!Local variables ------------------------------
  INTEGER                               :: IT
  
  op%flavors      = flavors
  op%activeFlavor = 0
  op%beta         = beta

  IF ( MOD(flavors,2) .NE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_init : flavors is not even        ")

!#ifdef CTQMC_CHECK
  op%meanError    = 0.d0
  op%checkNumber  = 0.d0
  op%tolerance    = 0.d0
  op%doCheck      = .FALSE.
!#endif
 ! DT_FREEIF(op%particles)
  IF(ALLOCATED(op%particles)) THEN
    DO IT = 1,flavors
      CALL ListCdagCoffdiag_destroy(op%particles(IT)) 
    END DO
    ABI_DATATYPE_DEALLOCATE(op%particles) 
  END IF
  DT_MALLOC(op%particles,(1:flavors))
  FREEIF(op%mat_U)
  MALLOC(op%mat_U,(1:flavors,1:flavors))
  FREEIF(op%overlaps)
  MALLOC(op%overlaps,(1:flavors,1:flavors))
  op%overlaps = 0.d0
 ! if(.NOT.ALLOCATED(op%updates)) write(6,*) "NOT ALLOCATED"
  !FREEIF(op%updates)
  IF(ALLOCATED(op%updates)) THEN 
    write(6,*) "ALLOCATED",size(op%updates)
    ABI_DEALLOCATE(op%updates) 
    write(6,*) "ALLOCATED",size(op%updates)
  END IF
  MALLOC(op%updates,(1:flavors))
  op%updates = 0.d0
  !CALL ImpurityOperatoroffdiag_computeU(op, U, J)
  !op%mat_U = U
  !IF ( ASSOCIATED(op%mu) ) FREE(op%mu)
  !MALLOC(op%mu,(1:flavors))
 
  !op%shift_mu = SUM(op%mat_U(:,1)) * .5d0 
  DO IT = 1,flavors
    !CALL ListCdagCoffdiag_init(op%particles(IT), DBLE(N)/beta,100) !FIXME size of the List

!  ----- For Each flavor, construct op%particles(IT)%list(100,2)
    CALL ListCdagCoffdiag_init(op%particles(IT),100) !FIXME size of the List
    op%particles(IT)%list(0,C_   ) = beta ! Empty orbital 
    op%particles(IT)%list(0,Cdag_) = 0.d0

   ! op%particles(IT)%list(0,Cdag_) = beta ! Full orbital 
   ! op%particles(IT)%list(0,C_)    = 0.d0
    op%particles(IT)%tail=0
  END DO
  op%activeFlavor = 0
END SUBROUTINE ImpurityOperatoroffdiag_init 
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_reset
!! NAME
!!  ImpurityOperatoroffdiag_reset
!!
!! FUNCTION
!!  reset operator
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
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

SUBROUTINE ImpurityOperatoroffdiag_reset(op)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_reset'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER                               :: IT

  op%activeFlavor = 0
  op%overlaps = 0.d0
  op%updates = 0.d0
!#ifdef CTQMC_CHECK
  op%meanError    = 0.d0
  op%checkNumber  = 0.d0
  op%tolerance    = 0.d0
  op%doCheck      = .FALSE.
!#endif
  DO IT = 1,op%flavors
    CALL ListCdagCoffdiag_clear(op%particles(IT)) 
    op%particles(IT)%list(0,C_   )    = op%beta ! Empty orbital 
    op%particles(IT)%list(0,Cdag_) = 0.d0
   ! op%particles(IT)%list(0,Cdag_) = op%beta ! Full orbital 
   ! op%particles(IT)%list(0,C_)    = 0.d0
  END DO

END SUBROUTINE ImpurityOperatoroffdiag_reset
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_computeU
!! NAME
!!  ImpurityOperatoroffdiag_computeU
!!
!! FUNCTION
!!  Compute an interaction matrix for t2g like interaction
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  U=Coulomb scrren interaction
!!  J=Hund couplage
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

SUBROUTINE ImpurityOperatoroffdiag_computeU(op, U, J)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_computeU'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  DOUBLE PRECISION      , INTENT(IN   ) :: U
  DOUBLE PRECISION      , INTENT(IN   ) :: J
!Local variables ------------------------------
  INTEGER                               :: flavor11
  INTEGER                               :: flavor12
  INTEGER                               :: flavor21
  INTEGER                               :: flavor22
  INTEGER                               :: flavors
  INTEGER                               :: flavors_2
  DOUBLE PRECISION                      :: Uprime

  Uprime = U - 2.d0 * J
  flavors = op%flavors
  flavors_2 = flavors / 2
  DO flavor11 = 1, flavors_2
    flavor12 = flavors - flavor11 + 1
    op%mat_U(flavor11, flavor11) = 0.d0
    op%mat_U(flavor12, flavor12) = 0.d0
    op%mat_U(flavor11+flavors_2, flavor11) = U
    op%mat_U(flavor12-flavors_2, flavor12) = U
    DO flavor21 = flavor11+1, flavors_2
      flavor22 = flavors - flavor21 + 1
      op%mat_U(flavor21, flavor11) = Uprime
      op%mat_U(flavor22-flavors_2, flavor12-flavors_2) = Uprime
      op%mat_U(flavor21+flavors_2, flavor11+flavors_2) = Uprime
      op%mat_U(flavor22, flavor12) = Uprime
    END DO
    DO flavor21 = flavor11+flavors_2+1, flavors
      flavor22 = flavors - flavor21 + 1
      op%mat_U(flavor21, flavor11) = Uprime - J
      op%mat_U(flavor22+flavors_2, flavor12-flavors_2) = Uprime - J
      op%mat_U(flavor21-flavors_2, flavor11+flavors_2) = Uprime - J
      op%mat_U(flavor22, flavor12) = Uprime - J
    END DO
  END DO
END SUBROUTINE ImpurityOperatoroffdiag_computeU 
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_setUmat
!! NAME
!!  ImpurityOperatoroffdiag_setUmat
!!
!! FUNCTION
!!  Set directly the U interaction matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  matU=interaction matrix
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

SUBROUTINE ImpurityOperatoroffdiag_setUmat(op, matU)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_setUmat'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN   ) :: matU

  IF ( SIZE(matU) .NE. op%flavors*op%flavors ) &
    CALL ERROR("ImpurityOperatoroffdiag_setUmat : Wrong interaction matrix")

  op%mat_U = matU
END SUBROUTINE ImpurityOperatoroffdiag_setUmat
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_activateParticle
!! NAME
!!  ImpurityOperatoroffdiag_activateParticle
!!
!! FUNCTION
!!  active a flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  flavor=the flavor
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

SUBROUTINE ImpurityOperatoroffdiag_activateParticle(op,flavor)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_activateParticle'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER               , INTENT(IN   ) :: flavor

  IF ( flavor .GT. op%flavors ) &
    CALL ERROR("ImpurityOperatoroffdiag_activateParticle : out of range  ")
  IF ( ALLOCATED(op%particles) ) THEN 
    op%activeFlavor   =  flavor
  ELSE
    CALL ERROR("ImpurityOperatoroffdiag_activateParticle : not allocated  ")
  END IF
END SUBROUTINE ImpurityOperatoroffdiag_activateParticle
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getAvailableTime
!! NAME
!!  ImpurityOperatoroffdiag_getAvailableTime
!!
!! FUNCTION
!!  get the time available and the position of the segment to consider
!!  negative if on a segment
!!  positive if outside a segment
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  time=time to look for
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_getAvailableTime=Time available
!!  position=position of the next segment
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_getAvailableTime(op, time, position)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getAvailableTime'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN   ) :: op
  DOUBLE PRECISION      , INTENT(IN   ) :: time
  INTEGER               , INTENT(OUT  ) :: position
!Local variables ------------------------------
  DOUBLE PRECISION                      :: t_avail
  INTEGER                               :: position_dwn
  INTEGER                               :: aF
#include "ListCdagCoffdiag_firstHigher.h"
  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_getAvailableTime : no active flav")
  
  IF ( op%particles(aF)%tail .EQ. 0 ) THEN
    t_avail = op%particles(aF)%list(0,C_) - op%particles(aF)%list(0,Cdag_)
    position = SIGN(1,INT(t_avail))
  ELSE
!    position = ListCdagCoffdiag_firstHigher( op%particles(aF), time ) 
#define list_1 op%particles(aF) 
#include "ListCdagCoffdiag_firstHigher"
#undef list_1
! FirstHigher is the index of the segment just higher.
! FirstHigher is positive, except if the next segment is found after
! winding around beta, in this case FirstHigher=-1
! ======================================================================
! if tail==1, firsthigher can be 1 or -1 and thus position_dwn is 0 or -2=>1
    position = firstHigher
    position_dwn = position - 1
    IF ( position_dwn .LE. 0) position_dwn = op%particles(aF)%tail
  
!    t_avail = (time - op%particles(aF)%list(position_dwn)) .MOD. op%beta
    t_avail = time - op%particles(aF)%list(position_dwn,C_)
    IF ( op%particles(aF)%list(position_dwn,Cdag_) .GT. time ) &
      t_avail = t_avail + op%beta 
  
    IF ( t_avail .GT. 0.d0 ) THEN  !! We are outside the position_dwn segment
!      t_avail = (op%particles(aF)%list(ABS(position)) - time ) .MOD. op%beta 
      t_avail = op%particles(aF)%list(ABS(position),Cdag_) - time 
      IF ( op%particles(aF)%list(ABS(position),Cdag_) .LT. time ) &
        t_avail = t_avail + op%beta
      ! ABS is used to prevent position to be -1 which is HERE the same as 1
    ELSE
      position = - position_dwn
    END IF
  END IF
  
    ImpurityOperatoroffdiag_getAvailableTime = t_avail

END FUNCTION ImpurityOperatoroffdiag_getAvailableTime
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getAvailedTime
!! NAME
!!  ImpurityOperatoroffdiag_getAvailedTime
!!
!! FUNCTION
!!  get the time available without the segment "position"
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  position=position of the segment
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_getAvailedTime=time available before ...
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_getAvailedTime(op, position)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getAvailedTime'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN   ) :: op
  INTEGER               , INTENT(IN   ) :: position
  DOUBLE PRECISION                      :: T_avail
  INTEGER                               :: Pup
  INTEGER                               :: ABSp
  INTEGER                               :: tail
  INTEGER                               :: aF

  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_getAvailedTime : no active flavor")
  ABSp = ABS(position)
!  position_up = (ABSposition+1).MOD.op%particles(aF)%tail
 tail = op%particles(aF)%tail
  MODCYCLE(ABSp+1,tail,Pup)
  IF ( position .GT. 0 ) THEN
    t_avail = op%particles(aF)%list(Pup, Cdag_) &
            - op%particles(aF)%list(ABSp,Cdag_)
  ELSE
    t_avail = op%particles(aF)%list(Pup ,C_) &
            - op%particles(aF)%list(ABSp,C_)
  END IF
  IF ( t_avail .LE. 0.d0 ) t_avail = t_avail + op%beta
  ImpurityOperatoroffdiag_getAvailedTime = t_avail
END FUNCTION ImpurityOperatoroffdiag_getAvailedTime
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_add
!! NAME
!!  ImpurityOperatoroffdiag_add
!!
!! FUNCTION
!!  add a segment to the active flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  CdagC_1=couple of times
!!  position_val=position of the CdagC_1 couple in the list
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  op=ImpurityOperatoroffdiag
!!   op%particles(aF)%list is updated
!!   op%overlaps  is updated
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE ImpurityOperatoroffdiag_add(op, CdagC_1, position_val)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_add'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1
  INTEGER               , INTENT(IN   ) :: position_val
!Local variables ------------------------------
  INTEGER                               :: position
  INTEGER                               :: aF
  INTEGER                               :: i
  DOUBLE PRECISION, DIMENSION(1:2)      :: C2modify
  DOUBLE PRECISION, DIMENSION(1:2)      :: C2add
  DOUBLE PRECISION                      :: TCdag
  DOUBLE PRECISION                      :: TC

  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_add : no active flavor           ")
  
  position = position_val

  IF ( CdagC_1(C_) .GT. CdagC_1(Cdag_) ) THEN ! Ajout d'un segment
    C2add = CdagC_1
  ELSE                                        ! Ajout d'un antisegment
    IF ( (op%particles(aF)%tail .EQ. 0) .AND. (ABS(op%particles(aF)%list(0,C_)) .LE. 0.0001d0)) THEN ! should be full orbital
      IF ( CdagC_1(Cdag_) .GT. op%beta ) THEN
!        CALL CdagC_init(C2add,CdagC_1%Cdag-op%beta,CdagC_1%C)
        ! From the IF condition and the creation of CdagC in TryAddRemove, we have
        ! CdagC_1(Cdag_) > beta
        ! CdagC_1(C_)    < beta
        C2add(Cdag_) = CdagC_1(Cdag_)-op%beta
        C2add(C_   ) = CdagC_1(C_)
        ! Now C2add(Cdag_) < beta
        ! still C2add(C_)  < beta
      ELSE
!        CALL CdagC_init(C2add,CdagC_1%Cdag,CdagC_1%C+op%beta)
        ! CdagC_1(Cdag_) < beta
        ! CdagC_1(C_)    < beta
        C2add(Cdag_) = CdagC_1(Cdag_)
        C2add(C_   ) = CdagC_1(C_)+op%beta
        ! C2add(Cdag_) < beta
        ! C2ass(C_)    > beta
      END IF
      position = 0
      ! See impurityoperator_init to understand this. This is due to the
      ! convention for the full orbital case.
      op%particles(aF)%list(0,C_   ) = op%beta
      op%particles(aF)%list(0,Cdag_) = 0.d0
    ELSE IF ( op%particles(aF)%tail .GT. 0 ) THEN
      position = ABS(position)
      TCdag = op%particles(aF)%list(position,Cdag_)
      TC    = CdagC_1(C_)
      IF ( TCdag .GT. TC ) TC = TC + op%beta
!      CALL CdagC_init(C2modify,TCdag,TC)
      C2modify(Cdag_) = TCdag
      C2modify(C_   ) = TC
  
!      TCdag    = CdagC_1%Cdag.MOD.op%beta
      MODCYCLE(CdagC_1(Cdag_),op%beta,TCdag)
      TC       = op%particles(aF)%list(position,C_)
!      CALL CdagC_init(C2add,TCdag,TC)
      C2add(Cdag_) = TCdag
      C2add(C_   ) = TC
  
      op%particles(aF)%list(position,:) = C2modify
      IF ( C2modify(Cdag_) .GT. C2add(Cdag_) ) THEN
        position = 0
!        C2add%C = C2add%C.MOD.op%beta
        MODCYCLE(C2add(C_),op%beta,C2add(C_))
      END IF
    ELSE
      CALL ERROR("ImpurityOperatoroffdiag_add : try to add an antisegment to an empty orbital")
    END IF
    position = position + 1
  END IF
  CALL ListCdagCoffdiag_insert(op%particles(aF), C2add, position)
  DO i = 1, op%flavors
    op%overlaps(i,aF) = op%overlaps(i,aF) + op%updates(i)
    op%overlaps(aF,i) = op%overlaps(i,aF)
  END DO

END SUBROUTINE ImpurityOperatoroffdiag_add
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getSegment
!! NAME
!!  ImpurityOperatoroffdiag_getSegment
!!
!! FUNCTION
!!  Return the segment at position_val
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  position_val=position of the asked segment
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_getSegment(2)=the couple of time
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

FUNCTION ImpurityOperatoroffdiag_getSegment(op,position_val)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getSegment'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER               , INTENT(IN   ) :: position_val
!Local variables ------------------------------
  INTEGER                               :: position
  INTEGER                               :: tail
  INTEGER                               :: aF
  DOUBLE PRECISION                      :: beta
  DOUBLE PRECISION                      :: ImpurityOperatoroffdiag_getSegment(1:2)

  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_getSegment : no active flavor    ")

  IF ( position_val .GT. 0 ) THEN
    ImpurityOperatoroffdiag_getSegment = op%particles(aF)%list(position_val,1:2)
  ELSE
    position = ABS(position_val)
    tail = op%particles(aF)%tail
    beta = op%beta
    ImpurityOperatoroffdiag_getSegment(C_)  = op%particles(aF)%list(position,C_)
    position = position + 1
    IF ( position .GT. tail ) THEN
      IF ( ImpurityOperatoroffdiag_getSegment(C_) .LT. beta ) THEN
        ImpurityOperatoroffdiag_getSegment(Cdag_) = op%particles(aF)%list(1,Cdag_) + beta
      ELSE
        ImpurityOperatoroffdiag_getSegment(Cdag_) = op%particles(aF)%list(1,Cdag_)
        ImpurityOperatoroffdiag_getSegment(C_)    = ImpurityOperatoroffdiag_getSegment(C_) -beta
      END IF
    ELSE
      ImpurityOperatoroffdiag_getSegment(Cdag_) = op%particles(aF)%list(position,Cdag_)
    END IF

  END IF
END FUNCTION ImpurityOperatoroffdiag_getSegment
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_remove
!! NAME
!!  ImpurityOperatoroffdiag_remove
!!
!! FUNCTION
!!  Remove a segment for the active flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  ieme=segment to remove
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

SUBROUTINE ImpurityOperatoroffdiag_remove(op,ieme)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_remove'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER               , INTENT(IN   ) :: ieme
!Local variables ------------------------------
  DOUBLE PRECISION, DIMENSION(1:2)      :: CdagC_1
  INTEGER                               :: position
  INTEGER                               :: position_dwn
  INTEGER                               :: i
  INTEGER                               :: tail
  INTEGER                               :: aF
!  DOUBLE PRECISION                      :: toRemove

  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_removeIeme : no active flavor    ")
  position = ABS(ieme)
  IF ( position .GT. op%particles(aF)%tail ) &
    CALL ERROR("ImpurityOperatoroffdiag_removeIeme : out of range        ")

  IF ( (ieme .LT. 0)  .AND. (op%particles(aF)%tail .GT. 1) ) THEN 
    position_dwn = position
!    position = (position+1).MOD.op%particles(aF)%tail
    tail = op%particles(aF)%tail
    MODCYCLE((position+1),tail,position)
    CdagC_1(Cdag_) = op%particles(aF)%list(position_dwn,Cdag_)
    CdagC_1(C_   ) = op%particles(aF)%list(position,C_)
    IF (position_dwn .GT. position) CdagC_1(C_) = CdagC_1(C_) + op%beta
!    toRemove  = op%particles(aF)%list(position)%C - (CdagC_1%C.MOD.op%beta)
!    CdagC_1%C = CdagC_1%C + toRemove  
    op%particles(aF)%list(position_dwn,:) = CdagC_1
  END IF

  IF ( position .EQ. 1 ) THEN
    SELECT CASE (ieme)
      CASE (1) 
        op%particles(aF)%list(0,C_   ) = op%beta
        op%particles(aF)%list(0,Cdag_) = 0.d0 
      CASE (-1)
        op%particles(aF)%list(0,C_   ) = 0.d0 
        op%particles(aF)%list(0,Cdag_) = op%beta
    END SELECT
  END IF
  CALL ListCdagCoffdiag_erase(op%particles(aF),position)
  DO i = 1, op%flavors
    op%overlaps(i,aF) = op%overlaps(i,aF) - op%updates(i)
    op%overlaps(aF,i) = op%overlaps(i,aF)
  END DO
END SUBROUTINE ImpurityOperatoroffdiag_remove
!!***



!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getsign
!! NAME
!!  ImpurityOperatoroffdiag_getsign
!!
!! FUNCTION
!!  Get the sign of the ratio of impurity traces
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2014 ABINIT group (B. Amadon)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op     = ImpurityOperatoroffdiag
!!  time2    = for segment/antisegment addition, end of segment
!!  position = for segment/antisegment removal, position  of segment/antisegment removed
!!  action = > 0.5 addition 
!!           < 0.5 removal
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_getsign = sign of ratio of impurity traces
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Ctqmcoffdiag_tryAddRemove
!!
!! CHILDREN
!!
!! SOURCE

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_getsign(op, time2, i, action, position)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getsign'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
  DOUBLE PRECISION, INTENT(IN) :: time2, action
  INTEGER ,  INTENT(IN) :: i,position
!Local variables ------------------------------
  INTEGER                            :: tailint
  DOUBLE PRECISION                   :: sign_imp
! ************************************************************************
  tailint=op%particles(op%activeflavor)%tail
  if(action < 0.5d0) then
    if(tailint>=1) then
      if ( op%particles(op%activeFlavor)%list(tailint,2)>op%beta ) then ! segment winds around
        if (i==1) then ! add segment do not change winding
           sign_imp = 1
        else if (i==2) then ! antisegment
           if(time2>op%beta) then ! suppress winding around
             sign_imp = -1
           else   ! winding around still here
             sign_imp = 1
           endif
        endif
      else ! segment do not wind around
        if (i==1) then ! segment
          if(time2>op%beta) then ! create winding
            sign_imp = -1
          else   ! do not create winding
            sign_imp = 1
          endif
        else if (i==2) then ! no winding in any case
          sign_imp = 1
        endif
      endif
    else if (tailint==0) then
      if (i==1) then ! segment
        if(time2>op%beta) then ! create winding
           sign_imp = -1
        else   ! do not create winding
           sign_imp = 1
        endif
      else if (i==2) then ! antisegment
        if(time2>op%beta) then ! do not create winding
          sign_imp = 1
        else   ! create winding
          sign_imp = -1
        endif
      endif
    endif
  else
    if ( op%particles(op%activeFlavor)%list(tailint,2)>op%beta ) then ! segment winds around
      if (i==1) then ! remove segment
        if(position==tailint) then ! suppress winding around
          sign_imp = -1
        else  ! winding around still here
          sign_imp = 1
        endif
      else if (i==2) then ! remove antisegment
        if(tailint==1) then ! if tailint=1, create full orbital
          sign_imp = -1
        else  ! if tailint >1 preserve winding
          sign_imp = 1
        endif
      endif
    else ! segments do not wind around
      if (i==1) then ! suppress segment do not change winding
        sign_imp = 1
      else if (i==2) then ! antisegment 
        if(abs(position)==tailint) then  ! create winding around only tailint >=1
          if(tailint==1)  then 
            sign_imp = 1
          else 
            sign_imp = -1
          endif
        else  !do not create winding around
          sign_imp = 1
        endif
      endif
    endif
  endif

  ImpurityOperatoroffdiag_getsign=sign_imp


END FUNCTION ImpurityOperatoroffdiag_getsign
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getNewOverlap
!! NAME
!!  ImpurityOperatoroffdiag_getNewOverlap
!!
!! FUNCTION
!!  Get the overlap induced by CdagC_1 in the current configuration
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  CdagC_1=the segment
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_getNewOverlap=overlap..
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_getNewOverlap(op, CdagC_1)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getNewOverlap'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN) :: CdagC_1
!Local variables ------------------------------
  DOUBLE PRECISION, DIMENSION(1:2)   :: CdagC_2
  DOUBLE PRECISION                   :: overlap
  DOUBLE PRECISION                   :: totalOverlap
  DOUBLE PRECISION                   :: sign
  INTEGER                            :: flavor
  INTEGER                            :: otherFlavor

  flavor = op%activeFlavor
  IF ( flavor .LE. 0 ) &
    CALL ERROR("ImpurityOperatoroffdiag_getNewOverlap : no active flavor ")
  IF ( CdagC_1(Cdag_) .LT. CdagC_1(C_) ) THEN ! segment C*C
    CdagC_2 = CdagC_1
    sign = -1.d0
  ELSE
    CdagC_2(C_) = CdagC_1(Cdag_)
    CdagC_2(Cdag_) = CdagC_1(C_)
    sign = 1.d0
  END IF

  totalOverlap = 0.d0

  DO otherFlavor = 1, op%flavors
    IF ( otherFlavor .EQ. flavor ) CYCLE
    overlap = ImpurityOperatoroffdiag_overlapSegFlav(op,CdagC_2(1:2),otherflavor)
    totalOverlap = totalOverlap &
                 + overlap * op%mat_U(otherFlavor,flavor)
                 !!write(6,*) "overlap 1",otherFlavor,overlap,totaloverlap
    op%updates(otherFlavor) = -sign * overlap
  END DO

  totalOverlap = totalOverlap * sign
  ImpurityOperatoroffdiag_getNewOverlap = totalOverlap

END FUNCTION ImpurityOperatoroffdiag_getNewOverlap
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_overlapSegFlav
!! NAME
!!  ImpurityOperatoroffdiag_overlapSegFlav
!!
!! FUNCTION
!!  Compute the overlap of a segment with a flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  CdagC_1=segment
!!  flavor=flavor to use
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_overlapSegFlav=overlap between CdagC_1 and flavor
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_overlapSegFlav(op,CdagC_1,flavor)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_overlapSegFlav'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
!  TYPE(CdagC)           , INTENT(IN) :: CdagC_1
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN) :: CdagC_1
  INTEGER               , INTENT(IN) :: flavor
!Local variables ------------------------------
!  TYPE(CdagC), DIMENSION(:), POINTER :: list => NULL()
  DOUBLE PRECISION                   :: totalCdag
  DOUBLE PRECISION                   :: totalC
  DOUBLE PRECISION                   :: beta
  DOUBLE PRECISION                   :: Time
  DOUBLE PRECISION                   :: Tmin
  DOUBLE PRECISION                   :: Tmax
  DOUBLE PRECISION                   :: TmaxBeta
  DOUBLE PRECISION                   :: TscanMin
  DOUBLE PRECISION                   :: TscanMax
!  DOUBLE PRECISION                   :: sign
  DOUBLE PRECISION                   :: C
  DOUBLE PRECISION                   :: Cdag
  DOUBLE PRECISION                   :: loop
  DOUBLE PRECISION                   :: itmin
  DOUBLE PRECISION                   :: itmax
  INTEGER                            :: tail
  INTEGER                            :: tp1
  INTEGER                            :: scanning 
  INTEGER                            :: imin 
  INTEGER                            :: imax 
  INTEGER                            :: loops
  INTEGER                            :: iloop
#include "ListCdagCoffdiag_firstHigher.h"

  beta = op%beta
  Tmin = CdagC_1(Cdag_)
  Tmax = CdagC_1(C_)
  itmin= 0.d0

!  TmaxBeta     = Tmax.MOD.beta
  MODCYCLE(Tmax,beta,TmaxBeta)

  tail = op%particles(flavor)%tail

  totalC = 0.d0
  totalCdag = 0.d0
  IF ( tail .NE. 0 ) THEN
    tp1  = tail + 1
    loop = 0.d0
!    imin = ListCdagCoffdiag_firstHigher( op%particles(flavor), Tmin ) - 1
    Time = Tmin
#define list_1 op%particles(flavor) 
#include "ListCdagCoffdiag_firstHigher"
    imin = firstHigher - 1

    SELECT CASE ( imin ) 
      CASE(0)
        scanning = tail
        loop = -1.d0
      CASE(-2)
        scanning = tail
      CASE DEFAULT
        scanning = imin
    END SELECT
!    imax = ListCdagCoffdiag_firstHigher( op%particles(flavor), TmaxBeta ) !- 1 Jamais atteint
    Time = TmaxBeta
#include "ListCdagCoffdiag_firstHigher"
#undef list_1
    imax = firstHigher

    TscanMin = Tmin
    TscanMax = Tmax

    ! Regarder avant 
    IF ( (imin .EQ. 0) ) THEN
      C = op%particles(flavor)%list(scanning,C_) +loop*  beta
      Cdag = op%particles(flavor)%list(scanning,Cdag_) +loop* beta
      itmax = MAX(TscanMin, Cdag)
      itmin = MIN(TscanMax, C   )

      IF ( itmin .GT. itmax ) THEN ! si egal alors overlap de 0
        totalC = totalC + itmin
        totalCdag = totalCdag + itmax
      END IF
      scanning = scanning+1
      IF ( scanning .EQ. tp1 ) THEN
        scanning = 1
      END IF
    END IF

    loops = imax - scanning
    IF ( TmaxBeta .NE. Tmax ) THEN
      loops = tail - loops
    ELSE IF ( imax .EQ. -1 ) THEN
      loops = tail - imin
    END IF

    !Comparer betement 2 segments
    DO iloop =0, loops
      C = op%particles(flavor)%list(scanning,C_)
      Cdag = op%particles(flavor)%list(scanning,Cdag_)
      itmax = MAX(TscanMin, Cdag)
      itmin = MIN(TscanMax,C)

      IF ( itmin .GT. itmax ) THEN ! si egal alors overla de 0
        totalC = totalC + itmin
        totalCdag = totalCdag + itmax
      END IF
      scanning = scanning + 1
      IF ( scanning .EQ. tp1 ) THEN
        scanning = 1
        IF ( itmin .EQ. TScanMax ) EXIT
        TscanMin = TscanMin - beta
        TscanMax = TscanMax - beta
      END IF
    END DO

    ! Regarder apres le segment
    IF ( (itmin .NE. TscanMax) ) THEN
      C = op%particles(flavor)%list(scanning,C_)
      Cdag = op%particles(flavor)%list(scanning,Cdag_) 
      itmax = MAX(TscanMin, Cdag)
      itmin = MIN(TscanMax,C)

      IF ( itmin .GT. itmax ) THEN ! si egal alors overla de 0
        totalC = totalC + itmin
        totalCdag = totalCdag + itmax
      END IF
    END IF
  ELSE IF ( op%particles(flavor)%list(0,C_) .EQ. 0.d0 ) THEN ! full orbital
      totalC    = Tmax
      totalCdag = Tmin
  END IF
!#ifdef CTQMC_CHECK
  IF ( op%doCheck .EQV. .TRUE. ) &
    CALL ImpurityOperatoroffdiag_checkOverlap(op, Tmin, Tmax,totalC-totalCdag,flavor)
!#endif
  ImpurityOperatoroffdiag_overlapSegFlav = totalC - totalCdag 
                 !!write(6,*) "overlap 2",totalC, totalCdag

END FUNCTION ImpurityOperatoroffdiag_overlapSegFlav
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_overlapFlavor
!! NAME
!!  ImpurityOperatoroffdiag_overlapFlavor
!!
!! FUNCTION
!!  Returns the overlap of flavor with the others
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  flavor=the one we want
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_overlapFlavor=result
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_overlapFlavor(op,flavor)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_overlapFlavor'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
  INTEGER,      OPTIONAL, INTENT(IN) :: flavor
!Local variables ------------------------------
  INTEGER                            :: otherFlavor
  DOUBLE PRECISION                   :: overlap
  DOUBLE PRECISION                   :: totalOverlap

  totalOverlap = 0.d0
  DO otherFlavor = 1, op%flavors
    IF ( otherFlavor .EQ. flavor ) CYCLE
    overlap = op%overlaps(otherFlavor,flavor)
    totalOverlap = totalOverlap &
                 + overlap * op%mat_U(otherFlavor,flavor)
  END DO

  ImpurityOperatoroffdiag_overlapFlavor = totalOverlap

END FUNCTION ImpurityOperatoroffdiag_overlapflavor
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_overlapSwap
!! NAME
!!  ImpurityOperatoroffdiag_overlapSwap
!!
!! FUNCTION
!!  compute the interaction energy using the segment configuration of
!! flavor2 and the interaction strengh of flavor1
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  flavor1=interaction value
!!  flavor2=configuration
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_overlapSwap=new overlap
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_overlapSwap(op,flavor1,flavor2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_overlapSwap'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
  INTEGER               , INTENT(IN) :: flavor1
  INTEGER               , INTENT(IN) :: flavor2
!Local variables ------------------------------
  INTEGER                            :: otherFlavor
  DOUBLE PRECISION                   :: overlap
  DOUBLE PRECISION                   :: totalOverlap

  totalOverlap = 0.d0
! Calcul l'overlap de flavor1 en utilisant la configuration de flavor2
  DO otherFlavor = 1, op%flavors
    IF ( otherFlavor .EQ. flavor2 ) THEN
      CYCLE
    ELSE IF ( otherFlavor .EQ. flavor1 ) THEN
      overlap = op%overlaps(otherFlavor,flavor2)
      ! here  <n_1 n_2>*U_12
      totalOverlap = totalOverlap &
                   + overlap * op%mat_U(otherFlavor,flavor2)
    ELSE
      overlap = op%overlaps(otherFlavor,flavor2)
      ! here  <n_3 n_2>*U_31 + <n_4 n_2>*U_41
      totalOverlap = totalOverlap &
                   + overlap * op%mat_U(otherFlavor,flavor1)
    END IF
  END DO

  ImpurityOperatoroffdiag_overlapSwap = totalOverlap

END FUNCTION ImpurityOperatoroffdiag_overlapSwap
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_swap
!! NAME
!!  ImpurityOperatoroffdiag_swap
!!
!! FUNCTION
!!  Swap to flavors
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurtiyOperator
!!  flavor1=to swap
!!  flavor2=to swap
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

SUBROUTINE ImpurityOperatoroffdiag_swap(op,flavor1, flavor2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_swap'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER               , INTENT(IN   ) :: flavor1
  INTEGER               , INTENT(IN   ) :: flavor2
!Local variables ------------------------------
  INTEGER                               :: iflavor
  DOUBLE PRECISION                      :: overlap_tmp
  TYPE(ListCdagCoffdiag)                       :: list_tmp

  DO iflavor = 1, op%flavors
    IF ( iflavor .NE. flavor1  .AND. iflavor .NE. flavor2) THEN
      overlap_tmp = op%overlaps(iflavor,flavor1)
      op%overlaps(iflavor,flavor1) = op%overlaps(iflavor,flavor2)
      op%overlaps(flavor1,iflavor) = op%overlaps(iflavor,flavor2)
      op%overlaps(iflavor,flavor2) = overlap_tmp
      op%overlaps(flavor2,iflavor) = overlap_tmp
    END IF
  END DO

  !CALL ListCdagCoffdiag_print(op%particles(flavor1),233)
  !CALL ListCdagCoffdiag_print(op%particles(flavor2),233)
  list_tmp = op%particles(flavor1)
  op%particles(flavor1) = op%particles(flavor2)
  op%particles(flavor2) = list_tmp
  CALL ListCdagCoffdiag_destroy(list_tmp)
  !CALL ListCdagCoffdiag_swap(op%particles(flavor1),op%particles(flavor2))
  !CALL ListCdagCoffdiag_print(op%particles(flavor1),233)
  !CALL ListCdagCoffdiag_print(op%particles(flavor2),233)

END SUBROUTINE ImpurityOperatoroffdiag_swap
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_overlapIJ
!! NAME
!!  ImpurityOperatoroffdiag_overlapIJ
!!
!! FUNCTION
!!  Compute overlap between two flavors
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  i=first flavor
!!  j=second flavor
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_overlapIJ=result
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_overlapIJ(op,i,j)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_overlapIJ'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER               , INTENT(IN) :: i
  INTEGER               , INTENT(IN) :: j
!Local variables ------------------------------
!  TYPE(ListCdagCoffdiag)       , POINTER    :: particle1 => NULL()
!  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: list1 => NULL()
  INTEGER                            :: tail1
  DOUBLE PRECISION, DIMENSION(1:2)   :: CdagC_1
  INTEGER                            :: isegment

!  particle1 => op%particles(i) 
!  list1     => particle1%list
  tail1 = op%particles(i)%tail

  ImpurityOperatoroffdiag_overlapIJ = 0.d0
  IF ( tail1 .EQ. 0 .AND. op%particles(i)%list(0,C_) .EQ. 0.d0 ) THEN ! FULL
!    CALL CdagC_init(CdagC_1,0.d0,op%beta)
    CdagC_1(Cdag_) = 0.d0
    CdagC_1(C_   ) = op%beta

    ImpurityOperatoroffdiag_overlapIJ = ImpurityOperatoroffdiag_overlapSegFlav(op,CdagC_1,j)
  ELSE IF ( tail1 .GT. 0) THEN
    op%activeFlavor = i
    DO isegment = 1, tail1
        CdagC_1(:) = op%particles(i)%list(isegment,1:2)
        ImpurityOperatoroffdiag_overlapIJ = ImpurityOperatoroffdiag_overlapIJ &
                   + ImpurityOperatoroffdiag_overlapSegFlav(op,CdagC_1,j)
    END DO
  END IF

END FUNCTION ImpurityOperatoroffdiag_overlapIJ
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_measDE
!! NAME
!!  ImpurityOperatoroffdiag_measDE
!!
!! FUNCTION
!!  measure double occupancy and interaction energy
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!
!! OUTPUT
!!  DE=array accumulating duoble occupancy and energy
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

SUBROUTINE ImpurityOperatoroffdiag_measDE(op,DE)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_measDE'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: DE
!Local variables ------------------------------
  DOUBLE PRECISION                              :: localD
  DOUBLE PRECISION                              :: totalE
  INTEGER                                       :: iflavor1
  INTEGER                                       :: iflavor2
  INTEGER                                       :: flavors

  IF ( .NOT. ALLOCATED(op%particles) ) &
    CALL ERROR("ImpurityOperatoroffdiag_measD : no particle set   ")

  totalE = 0.d0
  flavors = op%flavors
  DO iflavor1 = 1, flavors
    DO iflavor2 = iflavor1+1, flavors
      !localD = ImpurityOperatoroffdiag_overlapIJ(op,iflavor1,iflavor2) 
      localD = op%overlaps(iflavor2,iflavor1)
      DE(iflavor2,iflavor1) = DE(iflavor2,iflavor1) + localD  
      totalE = totalE + localD * op%mat_U(iflavor1,iflavor2)
    END DO
  END DO

  DE(1,1) = DE(1,1) + totalE

END SUBROUTINE ImpurityOperatoroffdiag_measDE
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_cleanOverlaps
!! NAME
!!  ImpurityOperatoroffdiag_cleanOverlaps
!!
!! FUNCTION
!!  Compute from scratch all overlaps
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
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

SUBROUTINE ImpurityOperatoroffdiag_cleanOverlaps(op)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_cleanOverlaps'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER                                       :: iflavor1
  INTEGER                                       :: iflavor2
  INTEGER                                       :: flavors

  IF ( .NOT. ALLOCATED(op%particles) ) &
    CALL ERROR("ImpurityOperatoroffdiag_cleanOverlap : no particle set   ")

  flavors = op%flavors
  DO iflavor1 = 1, flavors
    DO iflavor2 = iflavor1+1, flavors
      op%overlaps(iflavor2,iflavor1) = ImpurityOperatoroffdiag_overlapIJ(op,iflavor1,iflavor2) 
    END DO
  END DO

END SUBROUTINE ImpurityOperatoroffdiag_cleanOverlaps
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_measN
!! NAME
!!  ImpurityOperatoroffdiag_measN
!!
!! FUNCTION
!!  measure the number of electrons on flavor flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  flavor=the flavor
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_measN=number of electrons
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_measN(op,flavor)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_measN'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
  INTEGER,      OPTIONAL, INTENT(IN) :: flavor
!Local variables ------------------------------
  DOUBLE PRECISION                   :: totalCdag
  DOUBLE PRECISION                   :: totalC
  INTEGER                            :: scanning
  INTEGER                            :: aF

  IF ( PRESENT(flavor) ) THEN
    aF = flavor
  ELSE
    aF = op%activeFlavor
  END IF

  IF ( aF .LE. 0 ) & 
    CALL ERROR("ImpurityOperatoroffdiag_measN : no active flavor     ")

  totalC    = (op%particles(aF)%list(0,Cdag_) - op%particles(aF)%list(0,C_) + op%beta) * .5d0 
  totalCdag = 0.d0

  DO scanning = 1, op%particles(aF)%tail
    totalCdag = totalCdag + op%particles(aF)%list(scanning,Cdag_)
    totalC    = totalC    + op%particles(aF)%list(scanning,C_   )
  END DO

  ImpurityOperatoroffdiag_measN = totalC - totalCdag

END FUNCTION ImpurityOperatoroffdiag_measN
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_destroy
!! NAME
!!  ImpurityOperatoroffdiag_destroy
!!
!! FUNCTION
!!  destroy and deallocate
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
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

SUBROUTINE ImpurityOperatoroffdiag_destroy(op)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_destroy'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER                               :: IT

  IF ( ALLOCATED(op%particles) ) THEN
    DO IT = 1, op%flavors
      CALL ListCdagCoffdiag_destroy(op%particles(IT))
    END DO
    DT_FREE(op%particles)
  ENDIF
  FREEIF(op%mat_U)
  FREEIF(op%overlaps)
  FREEIF(op%updates)
  op%activeFlavor = 0
  op%flavors      = 0
  op%beta         = 0.d0
END SUBROUTINE ImpurityOperatoroffdiag_destroy
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getErrorOverlap
!! NAME
!!  ImpurityOperatoroffdiag_getErrorOverlap
!!
!! FUNCTION
!!  compute error on the overlap (numerical accumulation)
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!
!! OUTPUT
!!  DE=save the error
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

SUBROUTINE ImpurityOperatoroffdiag_getErrorOverlap(op,DE)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getErrorOverlap'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: DE
!Local variables ------------------------------
  DOUBLE PRECISION                              :: localD1
  DOUBLE PRECISION                              :: localD2
  DOUBLE PRECISION                              :: totalE1
  DOUBLE PRECISION                              :: totalE2
  INTEGER                                       :: iflavor1
  INTEGER                                       :: iflavor2
  INTEGER                                       :: flavors

  IF ( .NOT. ALLOCATED(op%particles) ) &
    CALL ERROR("ImpurityOperatoroffdiag_getErrorOverlap : no particle set ")

  totalE1 = 0.d0
  totalE2 = 0.d0
  flavors = op%flavors
  DO iflavor1 = 1, flavors
    DO iflavor2 = iflavor1+1, flavors
      localD1 = ImpurityOperatoroffdiag_overlapIJ(op,iflavor1,iflavor2) 
      localD2 = op%overlaps(iflavor2,iflavor1)
      totalE1 = totalE1 + localD1 * op%mat_U(iflavor1,iflavor2)
      totalE2 = totalE2 + localD2 * op%mat_U(iflavor1,iflavor2)
    END DO
  END DO

  DE(2,2) = ABS(totalE1 - totalE2)

END SUBROUTINE ImpurityOperatoroffdiag_getErrorOverlap
!!***
!#ifdef CTQMC_CHECK
!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_doCheck
!! NAME
!!  ImpurityOperatoroffdiag_doCheck
!!
!! FUNCTION
!!  set the check mechanism
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  opt_check=1||3 do check
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

SUBROUTINE ImpurityOperatoroffdiag_doCheck(op,opt_check)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_doCheck'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag) , INTENT(INOUT) :: op
  INTEGER                , INTENT(IN   ) :: opt_check

  IF ( opt_check .EQ. 1 .OR. opt_check .EQ. 3 ) &
    op%doCheck = .TRUE.
END SUBROUTINE ImpurityOperatoroffdiag_doCheck
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_checkOverlap
!! NAME
!!  ImpurityOperatoroffdiag_checkOverlap
!!
!! FUNCTION
!!  check the calculation of the overlap (very very slow routine)
!!  between Tmin and Tmax (c+ and c)
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  Tmin=c+
!!  Tmax=c
!!  iOverlap=input overlap (fast calculation)
!!  iflavor=active flavor
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

SUBROUTINE ImpurityOperatoroffdiag_checkOverlap(op, Tmin, Tmax, iOverlap, iflavor)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_checkOverlap'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(INOUT)  :: op
  DOUBLE PRECISION      , INTENT(IN   )  :: Tmin
  DOUBLE PRECISION      , INTENT(IN   )  :: Tmax
  DOUBLE PRECISION      , INTENT(IN   )  :: iOverlap
  INTEGER               , INTENT(IN   )  :: iflavor 
!Local variables ------------------------------
  INTEGER, PARAMETER                     :: size=10000000
  INTEGER                                :: imin
  INTEGER                                :: imax
  INTEGER                                :: imaxbeta
  INTEGER                                :: isegment
  INTEGER                                :: tail
  INTEGER(1), DIMENSION(1:size,1:2)      :: checktab 
  CHARACTER(LEN=4)                       :: a
  DOUBLE PRECISION                       :: dt
  DOUBLE PRECISION                       :: inv_dt
  DOUBLE PRECISION                       :: overlap
  DOUBLE PRECISION                       :: erreur
  DOUBLE PRECISION                       :: weight
  INTEGER :: try

  checktab = INT(0,1)
  overlap = 0.d0

  dt = op%beta / DBLE((size-1))
  inv_dt = 1.d0 / dt
  imin = INT(Tmin / dt + 0.5d0) + 1
  imax = INT(Tmax / dt + 0.5d0) + 1
  MODCYCLE(imax, size, imaxbeta)

  tail = op%particles(iflavor)%tail

  DO try = imin, MIN(imax,size)
    checktab(try,1)=INT(1,1)!IBSET(checktab(try,1),0)
  END DO

  IF ( imax .NE. imaxbeta ) THEN 
    DO try = 1, imaxbeta
      checktab(try,1)=INT(1,1)!IBSET(checktab(try,1),0)
    END DO
  END IF

  IF ( tail .NE. 0 ) THEN
    DO isegment=1, tail
      imin = INT(op%particles(iflavor)%list(isegment,Cdag_)* inv_dt + 0.5d0) + 1
      imax = INT(op%particles(iflavor)%list(isegment,C_   )* inv_dt + 0.5d0) + 1
      MODCYCLE(imax, size, imaxbeta)
      DO try = imin, MIN(imax,size)
        checktab(try,2)=INT(1,1)!IBSET(checktab(try,2),0)
      END DO
      IF ( imax .NE. imaxbeta ) THEN
        DO try = 1, imaxbeta
          checktab(try,2)=INT(1,1)!IBSET(checktab(try,2),0)
        END DO
      END IF
    END DO
  ELSE IF ( op%particles(iflavor)%list(0,C_) .EQ. 0.d0 ) THEN
    DO try = 1, size
      checktab(try,2)=INT(1,1)!IBSET(checktab(try,2),0)
    END DO
  END IF

  DO isegment = 1, size
    IF ( IAND(checktab(isegment,1),checktab(isegment,2)) .EQ. INT(1,1) ) &
      overlap = overlap + 1.d0
  END DO

  overlap = overlap * dt

  IF ( iOverlap .EQ. 0.d0 ) THEN
    erreur = ABS(overlap)
  ELSE
    erreur = ABS(overlap                - iOverlap)
  END IF
  weight = ABS(2.d0 * DBLE(tail) * dt - iOverlap)
  IF ( erreur .GT. weight  ) THEN 
    WRITE(a,'(I4)') INT(erreur*100.d0)
    CALL WARN("ImpurityOperatoroffdiag_checkOverlap : "//a//"%              ") 
  END IF
  IF ( iOverlap .LE. (2.d0 * DBLE(tail) * dt) ) &
    op%meanError = op%meanError + 1.d0
  op%checkNumber = op%checkNumber + 1.d0 !weight 

END SUBROUTINE ImpurityOperatoroffdiag_checkOverlap
!!***

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_getError
!! NAME
!!  ImpurityOperatoroffdiag_getErro
!!
!! FUNCTION
!!  get error on computing the overlap
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!
!! OUTPUT
!!  ImpurityOperatoroffdiag_getError=percentage error
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

DOUBLE PRECISION FUNCTION ImpurityOperatoroffdiag_getError(op)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_getError'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
!Local variables ------------------------------
!  DOUBLE PRECISION :: tolerance
  DOUBLE PRECISION :: error

  IF ( op%doCheck .EQV. .TRUE. ) THEN
    error     = ABS(op%meanError/op%checkNumber) 
!  tolerance = ABS(op%tolerance/op%checkNumber) 
    ImpurityOperatoroffdiag_getError = error 
  ELSE
    ImpurityOperatoroffdiag_getError = 0.d0
  END IF
END FUNCTION ImpurityOperatoroffdiag_getError
!!***
!#endif

!!****f* ABINIT/m_ImpurityOperatoroffdiag/ImpurityOperatoroffdiag_printLatex
!! NAME
!!  ImpurityOperatoroffdiag_printLatex
!!
!! FUNCTION
!!  print in a latex format all the configuration
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ImpurityOperatoroffdiag
!!  ostream=file stream
!!  isweep=current sweep number
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

SUBROUTINE ImpurityOperatoroffdiag_printLatex(op, ostream, isweep)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ImpurityOperatoroffdiag_printLatex'
!End of the abilint section

  TYPE(ImpurityOperatoroffdiag), INTENT(IN) :: op
  INTEGER               , INTENT(IN) :: ostream
  INTEGER               , INTENT(IN) :: isweep
!Local variables ------------------------------
  INTEGER                            :: flavors
  INTEGER                            :: iflavor
  INTEGER                            :: tail
  INTEGER                            :: it
  DOUBLE PRECISION                   :: C
  DOUBLE PRECISION                   :: Cdag
  INTEGER                            :: ordo
  INTEGER                            :: y
  INTEGER                            :: lines
  INTEGER                            :: letters
  DOUBLE PRECISION                   :: length

  flavors = op%flavors

  WRITE(ostream,'(A13)')    "\begin{frame}"
  WRITE(ostream,'(2x,A14)') "\begin{figure}"
  WRITE(ostream,'(4x,A28)') "\setlength{\unitlength}{1mm}"
  WRITE(ostream,'(4x,A23)') "\begin{picture}(104,90)"
  WRITE(ostream,'(6x,A29,I6,A2)') "\put(52,00){\makebox(0,0)[c]{",isweep,"}}"
  y = INT(90.d0/DBLE(flavors+1))
  DO iflavor = 1, flavors
    tail =  op%particles(iflavor)%tail
    ordo = iflavor * y
    lines = ordo - 1
    letters = ordo - 5
    WRITE(ostream,'(6x,A6,I2)') "%ligne", iflavor
    WRITE(ostream,'(6x,A41,I2,A16)') "\linethickness{0.5pt}\color{black}\put(2,",lines,"){\line(0,1){2}}"
    WRITE(ostream,'(6x,A7,I2,A24)')  "\put(2,",letters,"){\makebox(0,0)[c]{$0$}}"
    WRITE(ostream,'(6x,A7,I2,A18)')  "\put(2,",ordo,"){\line(1,0){100}}"
    WRITE(ostream,'(6x,A9,I2,A16)')  "\put(102,",lines,"){\line(0,1){2}}"
    WRITE(ostream,'(6x,A9,I2,A28)')  "\put(102,",letters,"){\makebox(0,0)[c]{$\beta$}}"
    DO it = 1, tail 
      Cdag = 2.d0+(op%particles(iflavor)%list(it,Cdag_)/op%beta*100.d0)
      C    = 2.d0+(op%particles(iflavor)%list(it,C_   )/op%beta*100.d0)
      length = C - Cdag
      IF ( op%particles(iflavor)%list(it,C_) .LE. op%beta ) THEN
        WRITE(ostream,'(8x,A9,I2)')             "%segments", it
        WRITE(ostream,'(8x,A37,F5.1,A1,I2,A13,F5.1,A2)') &
        "\linethickness{2pt}\color{black}\put(",Cdag,",",ordo,"){\line(1,0){",length,"}}"
        WRITE(ostream,'(8x,A5)')                "%Cdag"
        WRITE(ostream,'(8x,A12)')               "\color{blue}"
        WRITE(ostream,'(8x,A5,F5.1,A1,I2,A14)') "\put(",Cdag,",",ordo,"){\circle*{1}}"
        WRITE(ostream,'(8x,A2)')                "%C"
        WRITE(ostream,'(8x,A11)')               "\color{red}"
        WRITE(ostream,'(8x,A5,F5.1,A1,I2,A14)') "\put(",C,",",ordo,"){\circle*{1}}"
      ELSE
        WRITE(ostream,'(8x,A9,I2)')             "%segments", it
        WRITE(ostream,'(8x,A37,F5.1,A1,I2,A13,F5.1,A2)') &
        "\linethickness{2pt}\color{black}\put(",Cdag,",",ordo,"){\line(1,0){",102.d0-Cdag,"}}"
        WRITE(ostream,'(8x,A7,I2,A13,F5.1,A2)') "\put(2,",ordo,"){\line(1,0){",C-102.d0,"}}"
        WRITE(ostream,'(8x,A5)')                "%Cdag"
        WRITE(ostream,'(8x,A12)')               "\color{blue}"
        WRITE(ostream,'(8x,A5,F5.1,A1,I2,A14)') "\put(",Cdag,",",ordo,"){\circle*{1}}"
        WRITE(ostream,'(8x,A2)')                "%C"
        WRITE(ostream,'(8x,A11)')               "\color{red}"
        WRITE(ostream,'(8x,A5,F5.1,A1,I2,A14)') "\put(",C-100.d0,",",ordo,"){\circle*{1}}"
      END IF
    END DO
    IF ( tail .EQ. 0 .AND. op%particles(iflavor)%list(0,C_) .EQ. 0.d0 ) THEN 
      WRITE(ostream,'(8x,A9,I2)')      "%segments", it
      WRITE(ostream,'(8x,A39,I2,A18)') "\linethickness{2pt}\color{black}\put(2,",ordo,"){\line(1,0){100}}"
    END IF
  END DO
  WRITE(ostream,'(4x,A13)') "\end{picture}"
  WRITE(ostream,'(2x,A12)') "\end{figure}"
  WRITE(ostream,'(2x,A17)') "\transduration{0}"
  WRITE(ostream,'(A11)')    "\end{frame}"
END SUBROUTINE ImpurityOperatoroffdiag_printLatex
!!***

END MODULE m_ImpurityOperatoroffdiag
!!***


#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ImpurityOperator
!! NAME
!!  m_ImpurityOperator
!! 
!! FUNCTION 
!!  manage all related to Impurity
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
MODULE m_ImpurityOperator
USE m_ListCdagC
USE m_Global
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_ImpurityOperator/ImpurityOperator
!! NAME
!!  ImpurityOperator
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

TYPE, PUBLIC :: ImpurityOperator
  LOGICAL          _PRIVATE :: doCheck = .FALSE.
  INTEGER          _PRIVATE :: flavors
   !  Number of flavors
  INTEGER                   :: activeFlavor
   !  Flavor considered e.g when a segment is added


  DOUBLE PRECISION _PRIVATE          :: beta
   !  Inverse of temperature.

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)          :: mat_U 
   !  for iflavor1 and iflavor2, mat_U(iflavor1,iflavor2) is the
   !  coulomb interaction between iflavor1 and iflavor2.


  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) _PRIVATE :: overlaps   ! total overlaps
   !  for iflavor1 and iflavor2 overlaps(iflavor1,iflavor2) is the total
   !  overlap between segments of iflavor1 and segments of iflavor2.

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:  ) _PRIVATE :: updates    ! new_(anti)seg 
   !  For a given flavor (activeflavor), gives for each other flavors, the
   !  supplementary overlaps, called updates(otherflavor).

  TYPE(ListCdagC)                               _PRIVATE :: list_swap
  TYPE(ListCdagC) , ALLOCATABLE, DIMENSION(:  )          :: particles 
   !  for each flavor, particles(iflavor)%list(2,maxnbofsegment) 
   !  gives the beginning and end of each segment.

  DOUBLE PRECISION _PRIVATE :: checkNumber
  DOUBLE PRECISION _PRIVATE :: tolerance
  DOUBLE PRECISION _PRIVATE :: meanError
END TYPE ImpurityOperator
!!***

PUBLIC  :: ImpurityOperator_init
PUBLIC  :: ImpurityOperator_reset
PUBLIC  :: ImpurityOperator_computeU
PUBLIC  :: ImpurityOperator_setUmat
PUBLIC  :: ImpurityOperator_setMu
PUBLIC  :: ImpurityOperator_activateParticle
PUBLIC  :: ImpurityOperator_getAvailableTime
PUBLIC  :: ImpurityOperator_getAvailedTime
PUBLIC  :: ImpurityOperator_add
PUBLIC :: ImpurityOperator_getSegment
PUBLIC :: ImpurityOperator_getsign
PUBLIC  :: ImpurityOperator_remove
PUBLIC :: ImpurityOperator_getNewOverlap
PUBLIC  :: ImpurityOperator_getTraceAdd
PUBLIC  :: ImpurityOperator_getTraceRemove
PRIVATE :: ImpurityOperator_overlapSegFlav
PUBLIC  :: ImpurityOperator_overlapflavor
PUBLIC  :: ImpurityOperator_overlapSwap
PUBLIC  :: ImpurityOperator_swap
PRIVATE :: ImpurityOperator_overlapIJ
PUBLIC  :: ImpurityOperator_measDE
PUBLIC  :: ImpurityOperator_cleanOverlaps
PUBLIC  :: ImpurityOperator_measN
PUBLIC  :: ImpurityOperator_destroy
PUBLIC  :: ImpurityOperator_getErrorOverlap
PUBLIC  :: ImpurityOperator_doCheck
PRIVATE :: ImpurityOperator_checkOverlap
PUBLIC  :: ImpurityOperator_getError
PUBLIC  :: ImpurityOperator_printLatex

CONTAINS
!!***

!SUBROUTINE ImpurityOperator_init(this, flavors, beta, N)
!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_init
!! NAME
!!  ImpurityOperator_init
!!
!! FUNCTION
!!  Initialize and allocate
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurtiyOperator
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

SUBROUTINE ImpurityOperator_init(this, flavors, beta)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN   ) :: flavors
  !DOUBLE PRECISION      , INTENT(IN   ) :: U
  !DOUBLE PRECISION      , INTENT(IN   ) :: J
  DOUBLE PRECISION      , INTENT(IN   ) :: beta
  !INTEGER               , INTENT(IN   ) :: N
!Local variables ------------------------------
  INTEGER                               :: IT
  
  this%flavors      = flavors
  this%activeFlavor = 0
  this%beta         = beta

  IF ( MOD(flavors,2) .NE. 0 ) &
    CALL ERROR("ImpurityOperator_init : flavors is not even        ")

!#ifdef CTQMC_CHECK
  this%meanError    = 0.d0
  this%checkNumber  = 0.d0
  this%tolerance    = 0.d0
  this%doCheck      = .FALSE.
!#endif
  DT_FREEIF(this%particles)
  DT_MALLOC(this%particles,(1:flavors))
  FREEIF(this%mat_U)
  MALLOC(this%mat_U,(1:flavors,1:flavors))
  FREEIF(this%overlaps)
  MALLOC(this%overlaps,(1:flavors,1:flavors))
  this%overlaps = 0.d0
  FREEIF(this%updates)
  MALLOC(this%updates,(1:flavors))
  this%updates = 0.d0
  !CALL ImpurityOperator_computeU(this, U, J)
  !this%mat_U = U
  !IF ( ASSOCIATED(this%mu) ) FREE(this%mu)
  !MALLOC(this%mu,(1:flavors))
 
  !this%shift_mu = SUM(this%mat_U(:,1)) * .5d0 
  DO IT = 1,flavors
    !CALL ListCdagC_init(this%particles(IT), DBLE(N)/beta,100) !FIXME size of the List
    CALL ListCdagC_init(this%particles(IT),100) !FIXME size of the List
    this%particles(IT)%list(0,C_   ) = beta ! Empty orbital 
    this%particles(IT)%list(0,Cdag_) = 0.d0
!    this%particles(IT)%list(0)%Cdag = beta ! Full orbital 
!    this%particles(IT)%list(0)%C    = 0.d0
  END DO
  this%activeFlavor = 0
END SUBROUTINE ImpurityOperator_init 
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_reset
!! NAME
!!  ImpurityOperator_reset
!!
!! FUNCTION
!!  reset operator
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurtiyOperator
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

SUBROUTINE ImpurityOperator_reset(this)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                               :: IT

  this%activeFlavor = 0
  this%overlaps = 0.d0
  this%updates = 0.d0
!#ifdef CTQMC_CHECK
  this%meanError    = 0.d0
  this%checkNumber  = 0.d0
  this%tolerance    = 0.d0
  this%doCheck      = .FALSE.
!#endif
  DO IT = 1,this%flavors
    CALL ListCdagC_clear(this%particles(IT)) 
    this%particles(IT)%list(0,C_   )    = this%beta ! Empty orbital 
    this%particles(IT)%list(0,Cdag_) = 0.d0
!    this%particles(IT)%list(0)%Cdag = beta ! Full orbital 
!    this%particles(IT)%list(0)%C    = 0.d0
  END DO

END SUBROUTINE ImpurityOperator_reset
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_computeU
!! NAME
!!  ImpurityOperator_computeU
!!
!! FUNCTION
!!  Compute an interaction this for t2g like interaction
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_computeU(this, U, J)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
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
  flavors = this%flavors
  flavors_2 = flavors / 2
  DO flavor11 = 1, flavors_2
    flavor12 = flavors - flavor11 + 1
    this%mat_U(flavor11, flavor11) = 0.d0
    this%mat_U(flavor12, flavor12) = 0.d0
    this%mat_U(flavor11+flavors_2, flavor11) = U
    this%mat_U(flavor12-flavors_2, flavor12) = U
    DO flavor21 = flavor11+1, flavors_2
      flavor22 = flavors - flavor21 + 1
      this%mat_U(flavor21, flavor11) = Uprime
      this%mat_U(flavor22-flavors_2, flavor12-flavors_2) = Uprime
      this%mat_U(flavor21+flavors_2, flavor11+flavors_2) = Uprime
      this%mat_U(flavor22, flavor12) = Uprime
    END DO
    DO flavor21 = flavor11+flavors_2+1, flavors
      flavor22 = flavors - flavor21 + 1
      this%mat_U(flavor21, flavor11) = Uprime - J
      this%mat_U(flavor22+flavors_2, flavor12-flavors_2) = Uprime - J
      this%mat_U(flavor21-flavors_2, flavor11+flavors_2) = Uprime - J
      this%mat_U(flavor22, flavor12) = Uprime - J
    END DO
  END DO
END SUBROUTINE ImpurityOperator_computeU 
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_setUmat
!! NAME
!!  ImpurityOperator_setUmat
!!
!! FUNCTION
!!  Set directly the U interaction this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurtityOperator
!!  matU=interaction this
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

SUBROUTINE ImpurityOperator_setUmat(this, matU)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN   ) :: matU
  INTEGER :: iflavor1
  INTEGER :: iflavor2

  IF ( SIZE(matU) .NE. this%flavors*this%flavors ) &
    CALL ERROR("ImpurityOperator_setUmat : Wrong interaction this")

  DO iflavor1 = 1, this%flavors
    DO iflavor2 = iflavor1+1, this%flavors
      this%mat_U(iflavor1,iflavor2) = matU(iflavor1,iflavor2)
      this%mat_U(iflavor2,iflavor1) = matU(iflavor2,iflavor1)
    END DO
  END DO
END SUBROUTINE ImpurityOperator_setUmat
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_setMu
!! NAME
!!  ImpurityOperator_setMu
!!
!! FUNCTION
!!  Set directly the chemical potential
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurtityOperator
!!  mu=chimical potential
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

SUBROUTINE ImpurityOperator_setMu(this, mu)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN   ) :: mu
  INTEGER :: iflavor

  IF ( SIZE(mu) .NE. this%flavors ) &
    CALL ERROR("ImpurityOperator_setMu : Wrong chimical potentials")

  DO iflavor = 1, this%flavors
    this%mat_U(iflavor,iflavor) = mu(iflavor)
  END DO
END SUBROUTINE ImpurityOperator_setMu
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_activateParticle
!! NAME
!!  ImpurityOperator_activateParticle
!!
!! FUNCTION
!!  active a flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_activateParticle(this,flavor)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN   ) :: flavor

  IF ( flavor .GT. this%flavors ) &
    CALL ERROR("ImpurityOperator_activateParticle : out of range  ")
  IF ( ALLOCATED(this%particles) ) THEN 
    this%activeFlavor   =  flavor
  ELSE
    CALL ERROR("ImpurityOperator_activateParticle : not allocated  ")
  END IF
END SUBROUTINE ImpurityOperator_activateParticle
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getAvailableTime
!! NAME
!!  ImpurityOperator_getAvailableTime
!!
!! FUNCTION
!!  get the time available and the position of the segment to consider
!!  negative if on a segment
!!  positive if outside a segment
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  time=time to look for
!!
!! OUTPUT
!!  ImpurityOperator_getAvailableTime=Time available
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

DOUBLE PRECISION FUNCTION ImpurityOperator_getAvailableTime(this, time, position)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN   ) :: this
  DOUBLE PRECISION      , INTENT(IN   ) :: time
  INTEGER               , INTENT(OUT  ) :: position
!Local variables ------------------------------
  DOUBLE PRECISION                      :: t_avail
  INTEGER                               :: position_dwn
  INTEGER                               :: aF
#include "ListCdagC_firstHigher.h"
  aF = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperator_getAvailableTime : no active flav")
  
  IF ( this%particles(aF)%tail .EQ. 0 ) THEN
    t_avail = this%particles(aF)%list(0,C_) - this%particles(aF)%list(0,Cdag_)
    position = SIGN(1,INT(t_avail))
  ELSE
!    position = ListCdagC_firstHigher( this%particles(aF), time ) 
#define list_1 this%particles(aF) 
#include "ListCdagC_firstHigher"
#undef list_1
    position = firstHigher
    position_dwn = position - 1
    IF ( position_dwn .LE. 0) position_dwn = this%particles(aF)%tail
  
!    t_avail = (time - this%particles(aF)%list(position_dwn)) .MOD. this%beta
    t_avail = time - this%particles(aF)%list(position_dwn,C_)
    IF ( this%particles(aF)%list(position_dwn,Cdag_) .GT. time ) &
      t_avail = t_avail + this%beta 
  
    IF ( t_avail .GT. 0.d0 ) THEN  !! We are outside the position_dwn segment
!      t_avail = (this%particles(aF)%list(ABS(position)) - time ) .MOD. this%beta 
      t_avail = this%particles(aF)%list(ABS(position),Cdag_) - time 
      IF ( this%particles(aF)%list(ABS(position),Cdag_) .LT. time ) &
        t_avail = t_avail + this%beta
      ! ABS is used to prevent position to be -1 which is HERE the same as 1
    ELSE
      position = - position_dwn
    END IF
  END IF
  
    ImpurityOperator_getAvailableTime = t_avail

END FUNCTION ImpurityOperator_getAvailableTime
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getAvailedTime
!! NAME
!!  ImpurityOperator_getAvailedTime
!!
!! FUNCTION
!!  get the time available without the segment "position"
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  position=position of the segment
!!
!! OUTPUT
!!  ImpurityOperator_getAvailedTime=time available before ...
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

DOUBLE PRECISION FUNCTION ImpurityOperator_getAvailedTime(this, position)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN   ) :: this
  INTEGER               , INTENT(IN   ) :: position
  DOUBLE PRECISION                      :: T_avail
  INTEGER                               :: Pup
  INTEGER                               :: ABSp
  INTEGER                               :: tail
  INTEGER                               :: aF

  aF = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperator_getAvailedTime : no active flavor")
  ABSp = ABS(position)
!  position_up = (ABSposition+1).MOD.this%particles(aF)%tail
 tail = this%particles(aF)%tail
  MODCYCLE(ABSp+1,tail,Pup)
  IF ( position .GT. 0 ) THEN
    t_avail = this%particles(aF)%list(Pup, Cdag_) &
            - this%particles(aF)%list(ABSp,Cdag_)
  ELSE
    t_avail = this%particles(aF)%list(Pup ,C_) &
            - this%particles(aF)%list(ABSp,C_)
  END IF
  IF ( t_avail .LE. 0.d0 ) t_avail = t_avail + this%beta
  ImpurityOperator_getAvailedTime = t_avail
END FUNCTION ImpurityOperator_getAvailedTime
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_add
!! NAME
!!  ImpurityOperator_add
!!
!! FUNCTION
!!  add a segment to the active flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  CdagC_1=couple of times
!!  position_val=position of the CdagC_1 couple in the list
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  this=ImpurityOperatoroffdiag
!!   this%particles(aF)%list is updated
!!   this%overlaps  is updated
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE ImpurityOperator_add(this, CdagC_1, position_val)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
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

  aF = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperator_add : no active flavor           ")
  
  position = position_val

  IF ( CdagC_1(C_) .GT. CdagC_1(Cdag_) ) THEN ! Ajout d'un segment
    C2add = CdagC_1
  ELSE                                        ! Ajout d'un antisegment
    IF ( (this%particles(aF)%tail .EQ. 0) .AND. (this%particles(aF)%list(0,C_) .EQ. 0d0)) THEN ! should be full orbital
      IF ( CdagC_1(Cdag_) .GT. this%beta ) THEN
!        CALL CdagC_init(C2add,CdagC_1%Cdag-this%beta,CdagC_1%C)
        ! From the IF condition and the creation of CdagC in TryAddRemove, we have
        ! CdagC_1(Cdag_) > beta
        ! CdagC_1(C_)    < beta
        C2add(Cdag_) = CdagC_1(Cdag_)-this%beta
        C2add(C_   ) = CdagC_1(C_)
        ! Now C2add(Cdag_) < beta
        ! still C2add(C_)  < beta
      ELSE
!        CALL CdagC_init(C2add,CdagC_1%Cdag,CdagC_1%C+this%beta)
        ! CdagC_1(Cdag_) < beta
        ! CdagC_1(C_)    < beta
        C2add(Cdag_) = CdagC_1(Cdag_)
        C2add(C_   ) = CdagC_1(C_)+this%beta
        ! C2add(Cdag_) < beta
        ! C2ass(C_)    > beta
      END IF
      position = 0
      ! See impurityoperator_init to understand this. This is due to the
      ! convention for the full orbital case.
      this%particles(aF)%list(0,C_   ) = this%beta
      this%particles(aF)%list(0,Cdag_) = 0.d0
    ELSE IF ( this%particles(aF)%tail .GT. 0 ) THEN
      position = ABS(position)
      TCdag = this%particles(aF)%list(position,Cdag_)
      TC    = CdagC_1(C_)
      IF ( TCdag .GT. TC ) TC = TC + this%beta
!      CALL CdagC_init(C2modify,TCdag,TC)
      C2modify(Cdag_) = TCdag
      C2modify(C_   ) = TC
  
!      TCdag    = CdagC_1%Cdag.MOD.this%beta
      MODCYCLE(CdagC_1(Cdag_),this%beta,TCdag)
      TC       = this%particles(aF)%list(position,C_)
!      CALL CdagC_init(C2add,TCdag,TC)
      C2add(Cdag_) = TCdag
      C2add(C_   ) = TC
  
      this%particles(aF)%list(position,:) = C2modify
      IF ( C2modify(Cdag_) .GT. C2add(Cdag_) ) THEN
        position = 0
!        C2add%C = C2add%C.MOD.this%beta
        MODCYCLE(C2add(C_),this%beta,C2add(C_))
      END IF
    ELSE
      CALL ERROR("ImpurityOperator_add : try to add an antisegment to an empty orbital")
    END IF
    position = position + 1
  END IF
  CALL ListCdagC_insert(this%particles(aF), c2add, position)
  DO i = 1, this%flavors
    this%overlaps(i,aF) = this%overlaps(i,aF) + this%updates(i)
    this%overlaps(aF,i) = this%overlaps(i,aF)
  END DO

END SUBROUTINE ImpurityOperator_add
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getSegment
!! NAME
!!  ImpurityOperator_getSegment
!!
!! FUNCTION
!!  Return the segment at position_val
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  position_val=position of the asked segment
!!
!! OUTPUT
!!  ImpurityOperator_getSegment(2)=the couple of time
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

FUNCTION ImpurityOperator_getSegment(this,position_val)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN   ) :: position_val
!Local variables ------------------------------
  INTEGER                               :: position
  INTEGER                               :: tail
  INTEGER                               :: aF
  DOUBLE PRECISION                      :: beta
  DOUBLE PRECISION                      :: ImpurityOperator_getSegment(1:2)

  aF = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperator_getSegment : no active flavor    ")

  IF ( position_val .GT. 0 ) THEN
    ImpurityOperator_getSegment = this%particles(aF)%list(position_val,1:2)
  ELSE
    position = ABS(position_val)
    tail = this%particles(aF)%tail
    beta = this%beta
    ImpurityOperator_getSegment(C_)  = this%particles(aF)%list(position,C_)
    position = position + 1
    IF ( position .GT. tail ) THEN
      IF ( ImpurityOperator_getSegment(C_) .LT. beta ) THEN
        ImpurityOperator_getSegment(Cdag_) = this%particles(aF)%list(1,Cdag_) + beta
      ELSE
        ImpurityOperator_getSegment(Cdag_) = this%particles(aF)%list(1,Cdag_)
        ImpurityOperator_getSegment(C_)    = ImpurityOperator_getSegment(C_) -beta
      END IF
    ELSE
      ImpurityOperator_getSegment(Cdag_) = this%particles(aF)%list(position,Cdag_)
    END IF

  END IF
END FUNCTION ImpurityOperator_getSegment
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_remove
!! NAME
!!  ImpurityOperator_remove
!!
!! FUNCTION
!!  Remove a segment for the active flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_remove(this,ieme)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN   ) :: ieme
!Local variables ------------------------------
  DOUBLE PRECISION, DIMENSION(1:2)      :: CdagC_1
  INTEGER                               :: position
  INTEGER                               :: position_dwn
  INTEGER                               :: i
  INTEGER                               :: tail
  INTEGER                               :: aF
!  DOUBLE PRECISION                      :: toRemove

  aF = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("ImpurityOperator_removeIeme : no active flavor    ")
  position = ABS(ieme)
  IF ( position .GT. this%particles(aF)%tail ) &
    CALL ERROR("ImpurityOperator_removeIeme : out of range        ")

  IF ( (ieme .LT. 0)  .AND. (this%particles(aF)%tail .GT. 1) ) THEN 
    position_dwn = position
!    position = (position+1).MOD.this%particles(aF)%tail
    tail = this%particles(aF)%tail
    MODCYCLE((position+1),tail,position)
    CdagC_1(Cdag_) = this%particles(aF)%list(position_dwn,Cdag_)
    CdagC_1(C_   ) = this%particles(aF)%list(position,C_)
    IF (position_dwn .GT. position) CdagC_1(C_) = CdagC_1(C_) + this%beta
!    toRemove  = this%particles(aF)%list(position)%C - (CdagC_1%C.MOD.this%beta)
!    CdagC_1%C = CdagC_1%C + toRemove  
    this%particles(aF)%list(position_dwn,:) = CdagC_1
  END IF

  IF ( position .EQ. 1 ) THEN
    SELECT CASE (ieme)
      CASE (1) 
        this%particles(aF)%list(0,C_   ) = this%beta
        this%particles(aF)%list(0,Cdag_) = 0.d0 
      CASE (-1)
        this%particles(aF)%list(0,C_   ) = 0.d0 
        this%particles(aF)%list(0,Cdag_) = this%beta
    END SELECT
  END IF
  CALL ListCdagC_erase(this%particles(aF),position)
  DO i = 1, this%flavors
    this%overlaps(i,aF) = this%overlaps(i,aF) - this%updates(i)
    this%overlaps(aF,i) = this%overlaps(i,aF)
  END DO
END SUBROUTINE ImpurityOperator_remove
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getNewOverlap
!! NAME
!!  ImpurityOperator_getNewOverlap
!!
!! FUNCTION
!!  Get the overlap induced by CdagC_1 in the current configuration
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  CdagC_1=the segment
!!
!! OUTPUT
!!  ImpurityOperator_getNewOverlap=overlap..
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

DOUBLE PRECISION FUNCTION ImpurityOperator_getNewOverlap(this, CdagC_1)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN) :: CdagC_1
!Local variables ------------------------------
  DOUBLE PRECISION, DIMENSION(1:2)   :: CdagC_2
  DOUBLE PRECISION                   :: overlap
  DOUBLE PRECISION                   :: totalOverlap
  DOUBLE PRECISION                   :: sign
  INTEGER                            :: flavor
  INTEGER                            :: otherFlavor

  flavor = this%activeFlavor
  IF ( flavor .LE. 0 ) &
    CALL ERROR("ImpurityOperator_getNewOverlap : no active flavor ")
  IF ( CdagC_1(Cdag_) .LT. CdagC_1(C_) ) THEN ! segment C*C
    CdagC_2 = CdagC_1
    sign = -1.d0
  ELSE
    CdagC_2(C_) = CdagC_1(Cdag_)
    CdagC_2(Cdag_) = CdagC_1(C_)
    sign = 1.d0
  END IF

  totalOverlap = 0.d0

  DO otherFlavor = 1, this%flavors
    IF ( otherFlavor .EQ. flavor ) CYCLE
    overlap = ImpurityOperator_overlapSegFlav(this,CdagC_2(1:2),otherflavor)
    totalOverlap = totalOverlap &
                 + overlap * this%mat_U(otherFlavor,flavor)
    this%updates(otherFlavor) = -sign * overlap
  END DO

  totalOverlap = totalOverlap * sign
  ImpurityOperator_getNewOverlap = totalOverlap

END FUNCTION ImpurityOperator_getNewOverlap
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getsign
!! NAME
!!  ImpurityOperator_getsign
!!
!! FUNCTION
!!  Get the sign of the ratio of impurity traces
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (B. Amadon)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this     = ImpurityOperator
!!  time2    = for segment/antisegment addition, end of segment
!!  position = for segment/antisegment removal, position  of segment/antisegment removed
!!  action = > 0.5 addition 
!!           < 0.5 removal
!!
!! OUTPUT
!!  ImpurityOperator_getsign = sign of ratio of impurity traces
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Ctqmc_tryAddRemove
!!
!! CHILDREN
!!
!! SOURCE

DOUBLE PRECISION FUNCTION ImpurityOperator_getsign(this, time2, i, action, position)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
  DOUBLE PRECISION, INTENT(IN) :: time2, action
  INTEGER ,  INTENT(IN) :: i,position
!Local variables ------------------------------
  INTEGER                            :: tailint
  DOUBLE PRECISION                   :: sign_imp
! ************************************************************************
  tailint=this%particles(this%activeflavor)%tail
  if(action < 0.5d0) then
    if(tailint>=1) then
      if ( this%particles(this%activeFlavor)%list(tailint,2)>this%beta ) then ! segment winds around
        if (i==1) then ! add segment do not change winding
           sign_imp = 1
        else if (i==2) then ! antisegment
           if(time2>this%beta) then ! suppress winding around
             sign_imp = -1
           else   ! winding around still here
             sign_imp = 1
           endif
        endif
      else ! segment do not wind around
        if (i==1) then ! segment
          if(time2>this%beta) then ! create winding
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
        if(time2>this%beta) then ! create winding
           sign_imp = -1
        else   ! do not create winding
           sign_imp = 1
        endif
      else if (i==2) then ! antisegment
        if(time2>this%beta) then ! do not create winding
          sign_imp = 1
        else   ! create winding
          sign_imp = -1
        endif
      endif
    endif
  else
    if ( this%particles(this%activeFlavor)%list(tailint,2)>this%beta ) then ! segment winds around
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

  ImpurityOperator_getsign=sign_imp


END FUNCTION ImpurityOperator_getsign
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getTraceAdd
!! NAME
!!  ImpurityOperator_getTraceAdd
!!
!! FUNCTION
!!  Get the ratio of the traces of the impurity hamiltonien with and without the
!!  new (anti-)segment.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  CdagC_1=the segment
!!
!! OUTPUT
!!  ImpurityOperator_getTraceAdd = Tr[exp(-beta !H_impurity)c(t1)cd(t1)c(t2)cd(t2)...]/Tr[..]
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

FUNCTION ImpurityOperator_getTraceAdd(this, CdagC_1) RESULT(trace)

  TYPE(ImpurityOperator)          , INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1
  LOGICAL          :: antiseg
  DOUBLE PRECISION :: trace
  DOUBLE PRECISION :: overlap
  DOUBLE PRECISION :: length
  DOUBLE PRECISION :: antisym_sign
  DOUBLE PRECISION :: beta

  beta = this%beta
  antisym_sign = 1.0d0
  overlap   = ImpurityOperator_getNewOverlap(this,CdagC_1)
  length    = CdagC_1(C_   ) - CdagC_1(Cdag_)
  antiseg    = length .LT. 0.d0
  ! length > 0 if segment; < 0 if antisegment
  if ( this%particles(this%activeFlavor)%tail .GT. 0  .AND. &
       ( ( (.NOT. antiseg) .AND. CdagC_1(C_) .GT. beta ) .OR. &! for seg only
         ( antiseg .AND. CdagC_1(C_) .LT. beta .AND. CdagC_1(Cdag_) .GT. beta ) & ! SIGN > 0 for antiseg only
       ) &
     ) THEN
    antisym_sign = -1.d0
  ELSE IF ( this%particles(this%activeFlavor)%tail .EQ. 0 .AND. &
            ( ( (.NOT. antiseg) .AND. CdagC_1(C_) .GT. beta ) .OR. & ! >beta only possible for seg
              ( antiseg .AND. CdagC_1(Cdag_) .LT. beta ) & ! antiseg cdag < beta
            ) & 
          ) THEN
    antisym_sign = -1.d0
  END IF

  trace = antisym_sign * DEXP(this%mat_U(this%activeFlavor,this%activeFlavor)*length + overlap) 

END FUNCTION ImpurityOperator_getTraceAdd
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getTraceRemove
!! NAME
!!  ImpurityOperator_getTraceRemove
!!
!! FUNCTION
!!  Get the ratio of the traces of the impurity hamiltonien without and with the
!!  (anti-)segment.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  position=position of the segment
!!
!! OUTPUT
!!  ImpurityOperator_getTraceRemove = Tr[exp(-beta !H_impurity)c(t1)cd(t1)c(t2)cd(t2)...]/Tr[..]
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

FUNCTION ImpurityOperator_getTraceRemove(this, position) RESULT(trace)

  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN   ) :: position
  INTEGER          :: tail
  DOUBLE PRECISION :: trace
  DOUBLE PRECISION :: overlap
  DOUBLE PRECISION :: length
  DOUBLE PRECISION :: antisym_sign
  DOUBLE PRECISION :: last_C
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION, DIMENSION(1:2) :: CdagC_1

  beta = this%beta
  antisym_sign = 1.0d0

  CdagC_1    = ImpurityOperator_getSegment(this,position)
  length     = CdagC_1(C_) - CdagC_1(Cdag_)
  ! length > 0 if segment; < 0 if antisegment
  overlap    = ImpurityOperator_getNewOverlap(this,CdagC_1)

  tail = this%particles(this%activeFlavor)%tail
  last_C = this%particles(this%activeFlavor)%list(tail,C_)
  IF ( last_C .GT. beta ) THEN ! tail > 0 since if tail == 0 {0,beta}
    IF ( ( position .EQ. tail ) .OR. & ! only possible for segment (<0 if antiseg)
         ( length .LT. 0.d0 .AND. tail .EQ. 1 ) ) THEN
      antisym_sign = -1.d0
    END IF
  ELSE 
    IF ( tail .GT. 1 .AND. position .EQ. -tail ) & !tail>1 and last antisegment
    antisym_sign = -1.d0
  END IF

  trace = antisym_sign * DEXP(-this%mat_U(this%activeFlavor,this%activeFlavor)*length-overlap)

END FUNCTION ImpurityOperator_getTraceRemove
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_overlapSegFlav
!! NAME
!!  ImpurityOperator_overlapSegFlav
!!
!! FUNCTION
!!  Compute the overlap of a segment with a flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  CdagC_1=segment
!!  flavor=flavor to use
!!
!! OUTPUT
!!  ImpurityOperator_overlapSegFlav=overlap between CdagC_1 and flavor
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

DOUBLE PRECISION FUNCTION ImpurityOperator_overlapSegFlav(this,CdagC_1,flavor)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
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
#include "ListCdagC_firstHigher.h"

  beta = this%beta
  Tmin = CdagC_1(Cdag_)
  Tmax = CdagC_1(C_)
  itmin= 0.d0

!  TmaxBeta     = Tmax.MOD.beta
  MODCYCLE(Tmax,beta,TmaxBeta)

  tail = this%particles(flavor)%tail

  totalC = 0.d0
  totalCdag = 0.d0
  IF ( tail .NE. 0 ) THEN
    tp1  = tail + 1
    loop = 0.d0
!    imin = ListCdagC_firstHigher( this%particles(flavor), Tmin ) - 1
    Time = Tmin
#define list_1 this%particles(flavor) 
#include "ListCdagC_firstHigher"
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
!    imax = ListCdagC_firstHigher( this%particles(flavor), TmaxBeta ) !- 1 Jamais atteint
    Time = TmaxBeta
#include "ListCdagC_firstHigher"
#undef list_1
    imax = firstHigher

    TscanMin = Tmin
    TscanMax = Tmax

    ! Regarder avant 
    IF ( (imin .EQ. 0) ) THEN
      C = this%particles(flavor)%list(scanning,C_) +loop*  beta
      Cdag = this%particles(flavor)%list(scanning,Cdag_) +loop* beta
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
      C = this%particles(flavor)%list(scanning,C_)
      Cdag = this%particles(flavor)%list(scanning,Cdag_)
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
      C = this%particles(flavor)%list(scanning,C_)
      Cdag = this%particles(flavor)%list(scanning,Cdag_) 
      itmax = MAX(TscanMin, Cdag)
      itmin = MIN(TscanMax,C)

      IF ( itmin .GT. itmax ) THEN ! si egal alors overla de 0
        totalC = totalC + itmin
        totalCdag = totalCdag + itmax
      END IF
    END IF
  ELSE IF ( this%particles(flavor)%list(0,C_) .EQ. 0.d0 ) THEN ! full orbital
      totalC    = Tmax
      totalCdag = Tmin
  END IF
!#ifdef CTQMC_CHECK
  IF ( this%doCheck .EQV. .TRUE. ) &
    CALL ImpurityOperator_checkOverlap(this, Tmin, Tmax,totalC-totalCdag,flavor)
!#endif
  ImpurityOperator_overlapSegFlav = totalC - totalCdag 

END FUNCTION ImpurityOperator_overlapSegFlav
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_overlapFlavor
!! NAME
!!  ImpurityOperator_overlapFlavor
!!
!! FUNCTION
!!  Returns the overlap of flavor with the others
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  flavor=the one we want
!!
!! OUTPUT
!!  ImpurityOperator_overlapFlavor=result
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

DOUBLE PRECISION FUNCTION ImpurityOperator_overlapFlavor(this,flavor)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
  INTEGER,      OPTIONAL, INTENT(IN) :: flavor
!Local variables ------------------------------
  INTEGER                            :: otherFlavor
  DOUBLE PRECISION                   :: overlap
  DOUBLE PRECISION                   :: totalOverlap

  totalOverlap = 0.d0
  DO otherFlavor = 1, this%flavors
    IF ( otherFlavor .EQ. flavor ) CYCLE
    overlap = this%overlaps(otherFlavor,flavor)
    totalOverlap = totalOverlap &
                 + overlap * this%mat_U(otherFlavor,flavor)
  END DO

  ImpurityOperator_overlapFlavor = totalOverlap

END FUNCTION ImpurityOperator_overlapflavor
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_overlapSwap
!! NAME
!!  ImpurityOperator_overlapSwap
!!
!! FUNCTION
!!  compute the overlap of flavor1 with the configuration of flavor2
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  flavor1=interaction value
!!  flavor2=configuration
!!
!! OUTPUT
!!  ImpurityOperator_overlapSwap=new overlap
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

DOUBLE PRECISION FUNCTION ImpurityOperator_overlapSwap(this,flavor1,flavor2)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
  INTEGER               , INTENT(IN) :: flavor1
  INTEGER               , INTENT(IN) :: flavor2
!Local variables ------------------------------
  INTEGER                            :: otherFlavor
  DOUBLE PRECISION                   :: overlap
  DOUBLE PRECISION                   :: totalOverlap

  totalOverlap = 0.d0
! Calcul l'overlap de flavor1 en utilisant la configuration de flavor2
  DO otherFlavor = 1, this%flavors
    IF ( otherFlavor .EQ. flavor2 ) THEN
      CYCLE
    ELSE IF ( otherFlavor .EQ. flavor1 ) THEN
      overlap = this%overlaps(otherFlavor,flavor2)
      totalOverlap = totalOverlap &
                   + overlap * this%mat_U(otherFlavor,flavor2)
    ELSE
      overlap = this%overlaps(otherFlavor,flavor2)
      totalOverlap = totalOverlap &
                   + overlap * this%mat_U(otherFlavor,flavor1)
    END IF
  END DO

  ImpurityOperator_overlapSwap = totalOverlap

END FUNCTION ImpurityOperator_overlapSwap
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_swap
!! NAME
!!  ImpurityOperator_swap
!!
!! FUNCTION
!!  Swap to flavors
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurtiyOperator
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

SUBROUTINE ImpurityOperator_swap(this,flavor1, flavor2)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN   ) :: flavor1
  INTEGER               , INTENT(IN   ) :: flavor2
!Local variables ------------------------------
  INTEGER                               :: iflavor
  DOUBLE PRECISION                      :: overlap_tmp

  DO iflavor = 1, this%flavors
    IF ( iflavor .NE. flavor1  .AND. iflavor .NE. flavor2) THEN
      overlap_tmp = this%overlaps(iflavor,flavor1)
      this%overlaps(iflavor,flavor1) = this%overlaps(iflavor,flavor2)
      this%overlaps(flavor1,iflavor) = this%overlaps(iflavor,flavor2)
      this%overlaps(iflavor,flavor2) = overlap_tmp
      this%overlaps(flavor2,iflavor) = overlap_tmp
    END IF
  END DO

  !CALL ListCdagC_print(this%particles(flavor1),233)
  !CALL ListCdagC_print(this%particles(flavor2),233)
  CALL ListCdagC_assign(this%list_swap, this%particles(flavor1)) !list_swap = particle
  this%particles(flavor1) = this%particles(flavor2)
  this%particles(flavor2) = this%list_swap
  !CALL ListCdagC_swap(this%particles(flavor1),this%particles(flavor2))
  !CALL ListCdagC_print(this%particles(flavor1),233)
  !CALL ListCdagC_print(this%particles(flavor2),233)

END SUBROUTINE ImpurityOperator_swap
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_overlapIJ
!! NAME
!!  ImpurityOperator_overlapIJ
!!
!! FUNCTION
!!  Compute overlap between two flavors
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  i=first flavor
!!  j=second flavor
!!
!! OUTPUT
!!  ImpurityOperator_overlapIJ=result
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

DOUBLE PRECISION FUNCTION ImpurityOperator_overlapIJ(this,i,j)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  INTEGER               , INTENT(IN) :: i
  INTEGER               , INTENT(IN) :: j
!Local variables ------------------------------
!  TYPE(ListCdagC)       , POINTER    :: particle1 => NULL()
!  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: list1 => NULL()
  INTEGER                            :: tail1
  DOUBLE PRECISION, DIMENSION(1:2)   :: CdagC_1
  INTEGER                            :: isegment

!  particle1 => this%particles(i) 
!  list1     => particle1%list
  tail1 = this%particles(i)%tail

  ImpurityOperator_overlapIJ = 0.d0
  IF ( tail1 .EQ. 0 .AND. this%particles(i)%list(0,C_) .EQ. 0.d0 ) THEN ! FULL
!    CALL CdagC_init(CdagC_1,0.d0,this%beta)
    CdagC_1(Cdag_) = 0.d0
    CdagC_1(C_   ) = this%beta

    ImpurityOperator_overlapIJ = ImpurityOperator_overlapSegFlav(this,CdagC_1,j)
  ELSE IF ( tail1 .GT. 0) THEN
    this%activeFlavor = i
    DO isegment = 1, tail1
        CdagC_1(:) = this%particles(i)%list(isegment,1:2)
        ImpurityOperator_overlapIJ = ImpurityOperator_overlapIJ &
                   + ImpurityOperator_overlapSegFlav(this,CdagC_1,j)
    END DO
  END IF

END FUNCTION ImpurityOperator_overlapIJ
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_measDE
!! NAME
!!  ImpurityOperator_measDE
!!
!! FUNCTION
!!  measure double occupancy and interaction energy
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_measDE(this,DE)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: DE
!Local variables ------------------------------
  DOUBLE PRECISION                              :: localD
  DOUBLE PRECISION                              :: totalE
  INTEGER                                       :: iflavor1
  INTEGER                                       :: iflavor2
  INTEGER                                       :: flavors

  IF ( .NOT. ALLOCATED(this%particles) ) &
    CALL ERROR("ImpurityOperator_measD : no particle set   ")

  totalE = 0.d0
  flavors = this%flavors
  DO iflavor1 = 1, flavors
    DO iflavor2 = iflavor1+1, flavors
      !localD = ImpurityOperator_overlapIJ(this,iflavor1,iflavor2) 
      localD = this%overlaps(iflavor2,iflavor1)
      DE(iflavor2,iflavor1) = DE(iflavor2,iflavor1) + localD  
      totalE = totalE + localD * this%mat_U(iflavor1,iflavor2)
    END DO
  END DO

  DE(1,1) = DE(1,1) + totalE

END SUBROUTINE ImpurityOperator_measDE
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_cleanOverlaps
!! NAME
!!  ImpurityOperator_cleanOverlaps
!!
!! FUNCTION
!!  Compute from scratch all overlaps
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_cleanOverlaps(this)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                                       :: iflavor1
  INTEGER                                       :: iflavor2
  INTEGER                                       :: flavors

  IF ( .NOT. ALLOCATED(this%particles) ) &
    CALL ERROR("ImpurityOperator_cleanOverlap : no particle set   ")

  flavors = this%flavors
  DO iflavor1 = 1, flavors
    DO iflavor2 = iflavor1+1, flavors
      this%overlaps(iflavor2,iflavor1) = ImpurityOperator_overlapIJ(this,iflavor1,iflavor2) 
    END DO
  END DO

END SUBROUTINE ImpurityOperator_cleanOverlaps
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_measN
!! NAME
!!  ImpurityOperator_measN
!!
!! FUNCTION
!!  measure the number of electrons on flavor flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!  flavor=the flavor
!!
!! OUTPUT
!!  ImpurityOperator_measN=number of electrons
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

DOUBLE PRECISION FUNCTION ImpurityOperator_measN(this,flavor)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
  INTEGER,      OPTIONAL, INTENT(IN) :: flavor
!Local variables ------------------------------
  DOUBLE PRECISION                   :: totalCdag
  DOUBLE PRECISION                   :: totalC
  INTEGER                            :: scanning
  INTEGER                            :: aF

  IF ( PRESENT(flavor) ) THEN
    aF = flavor
  ELSE
    aF = this%activeFlavor
  END IF

  IF ( aF .LE. 0 ) & 
    CALL ERROR("ImpurityOperator_measN : no active flavor     ")

  totalC    = (this%particles(aF)%list(0,Cdag_) - this%particles(aF)%list(0,C_) + this%beta) * .5d0 
  totalCdag = 0.d0

  DO scanning = 1, this%particles(aF)%tail
    totalCdag = totalCdag + this%particles(aF)%list(scanning,Cdag_)
    totalC    = totalC    + this%particles(aF)%list(scanning,C_   )
  END DO

  ImpurityOperator_measN = totalC - totalCdag

END FUNCTION ImpurityOperator_measN
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_destroy
!! NAME
!!  ImpurityOperator_destroy
!!
!! FUNCTION
!!  destroy and deallocate
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_destroy(this)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                               :: IT

  IF ( ALLOCATED(this%particles) ) THEN
    DO IT = 1, this%flavors
      CALL ListCdagC_destroy(this%particles(IT))
    END DO
    DT_FREE(this%particles)
  ENDIF
  CALL ListCdagC_destroy(this%list_swap)
  FREEIF(this%mat_U)
  FREEIF(this%overlaps)
  FREEIF(this%updates)
  this%activeFlavor = 0
  this%flavors      = 0
  this%beta         = 0.d0
END SUBROUTINE ImpurityOperator_destroy
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getErrorOverlap
!! NAME
!!  ImpurityOperator_getErrorOverlap
!!
!! FUNCTION
!!  compute error on the overlap (numerical accumulation)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_getErrorOverlap(this,DE)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: DE
!Local variables ------------------------------
  DOUBLE PRECISION                              :: localD1
  DOUBLE PRECISION                              :: localD2
  DOUBLE PRECISION                              :: totalE1
  DOUBLE PRECISION                              :: totalE2
  INTEGER                                       :: iflavor1
  INTEGER                                       :: iflavor2
  INTEGER                                       :: flavors

  IF ( .NOT. ALLOCATED(this%particles) ) &
    CALL ERROR("ImpurityOperator_getErrorOverlap : no particle set ")

  totalE1 = 0.d0
  totalE2 = 0.d0
  flavors = this%flavors
  DO iflavor1 = 1, flavors
    DO iflavor2 = iflavor1+1, flavors
      localD1 = ImpurityOperator_overlapIJ(this,iflavor1,iflavor2) 
      localD2 = this%overlaps(iflavor2,iflavor1)
      totalE1 = totalE1 + localD1 * this%mat_U(iflavor1,iflavor2)
      totalE2 = totalE2 + localD2 * this%mat_U(iflavor1,iflavor2)
    END DO
  END DO

  DE(2,2) = ABS(totalE1 - totalE2)

END SUBROUTINE ImpurityOperator_getErrorOverlap
!!***
!#ifdef CTQMC_CHECK
!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_doCheck
!! NAME
!!  ImpurityOperator_doCheck
!!
!! FUNCTION
!!  set the check mechanism
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_doCheck(this,opt_check)

!Arguments ------------------------------------
  TYPE(ImpurityOperator) , INTENT(INOUT) :: this
  INTEGER                , INTENT(IN   ) :: opt_check

  IF ( opt_check .EQ. 1 .OR. opt_check .EQ. 3 ) &
    this%doCheck = .TRUE.
END SUBROUTINE ImpurityOperator_doCheck
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_checkOverlap
!! NAME
!!  ImpurityOperator_checkOverlap
!!
!! FUNCTION
!!  check the calculation of the overlap (very very slow routine)
!!  between Tmin and Tmax (c+ and c)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_checkOverlap(this, Tmin, Tmax, iOverlap, iflavor)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(INOUT)  :: this
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

  dt = this%beta / DBLE((size-1))
  inv_dt = 1.d0 / dt
  imin = INT(Tmin / dt + 0.5d0) + 1
  imax = INT(Tmax / dt + 0.5d0) + 1
  MODCYCLE(imax, size, imaxbeta)

  tail = this%particles(iflavor)%tail

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
      imin = INT(this%particles(iflavor)%list(isegment,Cdag_)* inv_dt + 0.5d0) + 1
      imax = INT(this%particles(iflavor)%list(isegment,C_   )* inv_dt + 0.5d0) + 1
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
  ELSE IF ( this%particles(iflavor)%list(0,C_) .EQ. 0.d0 ) THEN
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
    CALL WARN("ImpurityOperator_checkOverlap : "//a//"%              ") 
  END IF
  IF ( iOverlap .LE. (2.d0 * DBLE(tail) * dt) ) &
    this%meanError = this%meanError + 1.d0
  this%checkNumber = this%checkNumber + 1.d0 !weight 

END SUBROUTINE ImpurityOperator_checkOverlap
!!***

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_getError
!! NAME
!!  ImpurityOperator_getErro
!!
!! FUNCTION
!!  get error on computing the overlap
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
!!
!! OUTPUT
!!  ImpurityOperator_getError=percentage error
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

DOUBLE PRECISION FUNCTION ImpurityOperator_getError(this)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
!Local variables ------------------------------
!  DOUBLE PRECISION :: tolerance
  DOUBLE PRECISION :: error

  IF ( this%doCheck .EQV. .TRUE. ) THEN
    error     = ABS(this%meanError/this%checkNumber) 
!  tolerance = ABS(this%tolerance/this%checkNumber) 
    ImpurityOperator_getError = error 
  ELSE
    ImpurityOperator_getError = 0.d0
  END IF
END FUNCTION ImpurityOperator_getError
!!***
!#endif

!!****f* ABINIT/m_ImpurityOperator/ImpurityOperator_printLatex
!! NAME
!!  ImpurityOperator_printLatex
!!
!! FUNCTION
!!  print in a latex format all the configuration
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ImpurityOperator
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

SUBROUTINE ImpurityOperator_printLatex(this, ostream, isweep)

!Arguments ------------------------------------
  TYPE(ImpurityOperator), INTENT(IN) :: this
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

  flavors = this%flavors

  WRITE(ostream,'(A13)')    "\begin{frame}"
  WRITE(ostream,'(2x,A14)') "\begin{figure}"
  WRITE(ostream,'(4x,A28)') "\setlength{\unitlength}{1mm}"
  WRITE(ostream,'(4x,A23)') "\begin{picture}(104,90)"
  WRITE(ostream,'(6x,A29,I6,A2)') "\put(52,00){\makebox(0,0)[c]{",isweep,"}}"
  y = INT(90.d0/DBLE(flavors+1))
  DO iflavor = 1, flavors
    tail =  this%particles(iflavor)%tail
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
      Cdag = 2.d0+(this%particles(iflavor)%list(it,Cdag_)/this%beta*100.d0)
      C    = 2.d0+(this%particles(iflavor)%list(it,C_   )/this%beta*100.d0)
      length = C - Cdag
      IF ( this%particles(iflavor)%list(it,C_) .LE. this%beta ) THEN
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
    IF ( tail .EQ. 0 .AND. this%particles(iflavor)%list(0,C_) .EQ. 0.d0 ) THEN 
      WRITE(ostream,'(8x,A9,I2)')      "%segments", it
      WRITE(ostream,'(8x,A39,I2,A18)') "\linethickness{2pt}\color{black}\put(2,",ordo,"){\line(1,0){100}}"
    END IF
  END DO
  WRITE(ostream,'(4x,A13)') "\end{picture}"
  WRITE(ostream,'(2x,A12)') "\end{figure}"
  WRITE(ostream,'(2x,A17)') "\transduration{0}"
  WRITE(ostream,'(A11)')    "\end{frame}"
END SUBROUTINE ImpurityOperator_printLatex
!!***

END MODULE m_ImpurityOperator
!!***

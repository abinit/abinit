
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_BathOperator
!! NAME
!!  m_BathOperator
!! 
!! FUNCTION 
!!  Manage all stuff related to the bath for the 
!!  simgle Anderson Impurity Model
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
MODULE m_BathOperator
USE m_MatrixHyb
USE m_Vector
USE m_VectorInt
USE m_Global
USE m_ListCdagC

IMPLICIT NONE

!!***

PRIVATE

!!****t* m_BathOperator/BathOperator
!! NAME
!!  BathOperator
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

TYPE, PUBLIC :: BathOperator
  LOGICAL _PRIVATE :: set         = .FALSE.
  LOGICAL          :: MAddFlag    = .FALSE. ! Set to true if we can compute a new M (see updateDetXX)
  LOGICAL          :: MRemoveFlag = .FALSE. ! Set to true if we can compute a new M (see updateDetXX)
  LOGICAL _PRIVATE :: antiShift   = .FALSE. ! shift when M is updated with antiseg
  LOGICAL _PRIVATE :: doCheck     = .FALSE.
  INTEGER _PRIVATE :: flavors
  INTEGER          :: activeFlavor
  INTEGER _PRIVATE :: samples
  INTEGER _PRIVATE :: sizeHybrid
  INTEGER _PRIVATE :: updatePosRow
  INTEGER _PRIVATE :: updatePosCol
  INTEGER _PRIVATE :: iTech
  INTEGER _PRIVATE :: checkNumber
  DOUBLE PRECISION _PRIVATE                   :: beta
  DOUBLE PRECISION _PRIVATE                   :: dt
  DOUBLE PRECISION _PRIVATE                   :: inv_dt
  DOUBLE PRECISION _PRIVATE                   :: meanError
  DOUBLE PRECISION _PRIVATE                   :: S
  DOUBLE PRECISION _PRIVATE                   :: Stau
  DOUBLE PRECISION _PRIVATE                   :: Stilde
  TYPE(Vector)     _PRIVATE                   :: R 
  TYPE(Vector)     _PRIVATE                   :: Q 
  TYPE(Vector)     _PRIVATE                   :: Rtau
  TYPE(Vector)     _PRIVATE                   :: Qtau
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) _PRIVATE :: F ! sample,Flavors
  TYPE(MatrixHyb) , ALLOCATABLE, DIMENSION(:)            :: M  ! Flavors
  TYPE(MatrixHyb) , ALLOCATABLE, DIMENSION(:)   _PRIVATE :: M_update  ! Flavors
END TYPE BathOperator
!!***

PUBLIC  :: BathOperator_init
PUBLIC  :: BathOperator_reset
PUBLIC  :: BathOperator_activateParticle
PRIVATE :: BathOperator_hybrid
PUBLIC  :: BathOperator_getDetAdd
PUBLIC  :: BathOperator_getDetRemove
PUBLIC  :: BathOperator_getDetF
PUBLIC  :: BathOperator_setMAdd
PUBLIC  :: BathOperator_setMRemove
PUBLIC  :: BathOperator_swap
PUBLIC  :: BathOperator_initF
PUBLIC  :: BathOperator_setF
PUBLIC  :: BathOperator_printF
PUBLIC  :: BathOperator_printM
PUBLIC  :: BathOperator_destroy
PUBLIC  :: BathOperator_doCheck
PRIVATE :: BathOperator_checkM
PUBLIC  :: BathOperator_getError

CONTAINS
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_init
!! NAME
!!  BathOperator_init
!!
!! FUNCTION
!!  Initialize and allocate data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath object
!!  flavors=numbers of flavors we have (including spin)
!!  samples=Time slices in the input file
!!  beta=inverse temperature
!!  iTech=imaginary time or frequencies
!!  It is imposes to imaginary time
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

SUBROUTINE BathOperator_init(this, flavors, samples, beta, iTech)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
  INTEGER           , INTENT(IN   ) :: flavors
  INTEGER           , INTENT(IN   ) :: samples
  DOUBLE PRECISION  , INTENT(IN   ) :: beta
!Local variables ------------------------------
  INTEGER           , INTENT(IN   ) :: iTech
  INTEGER                           :: it

  this%MAddFlag     = .FALSE.
  this%MRemoveFlag  = .FALSE.
  this%flavors      = flavors
  this%beta         = beta
  this%samples      = samples
  this%sizeHybrid   = samples + 1
  this%dt      = beta / DBLE(samples)
  this%inv_dt  = DBLE(samples) / beta
  this%activeFlavor= 0 
  this%updatePosRow = 0
  this%updatePosCol = 0
  this%iTech        = iTech
!#ifdef CTQMC_CHECK
  this%checkNumber  = 0
  this%meanError    = 0.d0
  this%doCheck = .FALSE.
!#endif

  FREEIF(this%F)
  MALLOC(this%F,(1:this%sizeHybrid+1,1:flavors))
  DT_FREEIF(this%M)
  DT_MALLOC(this%M,(1:flavors))
  DT_FREEIF(this%M_update)
  DT_MALLOC(this%M_update,(1:flavors))
  
  CALL Vector_init(this%R,100)
  CALL Vector_init(this%Q,100)
  CALL Vector_init(this%Rtau,100)
  CALL Vector_init(this%Qtau,100)

  DO it = 1, flavors
    CALL MatrixHyb_init(this%M(it),this%iTech,size=Global_SIZE,Wmax=samples) !FIXME Should be consistent with ListCagC
    CALL MatrixHyb_init(this%M_update(it),this%iTech,size=Global_SIZE,Wmax=samples) !FIXME Should be consistent with ListCagC
  END DO
  this%F       = 0.d0
  this%set     = .TRUE.
  
END SUBROUTINE BathOperator_init
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_reset
!! NAME
!!  BathOperator_reset
!!
!! FUNCTION
!!  Reset all internal variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator to reset
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

SUBROUTINE BathOperator_reset(this)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                           :: it
  this%MAddFlag     = .FALSE.
  this%MRemoveFlag  = .FALSE.
  this%activeFlavor = 0 
  this%updatePosRow = 0
  this%updatePosCol = 0
!#ifdef CTQMC_CHECK
  this%checkNumber  = 0
  this%meanError    = 0.d0
!#endif
  this%doCheck = .FALSE.
  CALL Vector_clear(this%R)
  CALL Vector_clear(this%Q)
  CALL Vector_clear(this%Rtau)
  CALL Vector_clear(this%Qtau)

  DO it = 1, this%flavors
    CALL MatrixHyb_clear(this%M(it)) !FIXME Should be consistent with ListCagC
  END DO
  this%F       = 0.d0

END SUBROUTINE BathOperator_reset
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_activateParticle
!! NAME
!!  BathOperator_activateParticle
!!
!! FUNCTION
!!  Just save on wicht flavor we are working
!!  It is better to use the macro defined in defs.h
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  flavor=the flavor to activate
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

SUBROUTINE BathOperator_activateParticle(this,flavor)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER           , INTENT(IN   ) :: flavor

  IF ( flavor .GT. this%flavors ) &
    CALL ERROR("BathOperator_activateParticle : out of range      ")
  IF ( this%set .EQV. .TRUE. .AND. ALLOCATED(this%M) ) THEN 
    this%activeFlavor =  flavor
    this%MAddFlag     = .FALSE.
    this%MRemoveFlag  = .FALSE.
  ELSE
    CALL ERROR("BathOperator_activateParticle : not allocated      ")
  END IF
END SUBROUTINE BathOperator_activateParticle
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_hybrid
!! NAME
!!  BathOperator_hybrid
!!
!! FUNCTION
!!  Compute the hybridization for the active flavor
!!  at time time
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  time=time  F(time)
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

DOUBLE PRECISION FUNCTION BathOperator_hybrid(this,time)

  TYPE(BathOperator), INTENT(IN) :: this
  DOUBLE PRECISION  , INTENT(IN) :: time
#include "BathOperator_hybrid.h"

  IF ( this%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperator_hybrid : no active hybrid func        ")
#include "BathOperator_hybrid"
  BathOperator_hybrid = hybrid

END FUNCTION BathOperator_hybrid
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_getDetAdd
!! NAME
!!  BathOperator_getDetAdd
!!
!! FUNCTION
!!  Compute the determinant ratio when a (anti)segment
!!  is trying to be added and store some array for setMadd
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  CdagC_1=segment to be added
!!  position=ordered position of the Cdag time
!!  particle=full list of CdagC for activeFlavor
!!
!! OUTPUT
!!  BathOperator_getDetAdd=the det 
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
DOUBLE PRECISION  FUNCTION BathOperator_getDetAdd(this,CdagC_1, position, particle)

!Arguments ------------------------------------
  TYPE(BathOperator)      , INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1
  INTEGER                 , INTENT(IN   ) :: position  
  TYPE(ListCdagC), INTENT(IN   ) :: particle   
!Local variables-------------------------------
  INTEGER                                 :: it1
  INTEGER                                 :: it2
  INTEGER                                 :: it3
  INTEGER                                 :: tail
  INTEGER                                 :: new_tail
  DOUBLE PRECISION                        :: C
  DOUBLE PRECISION                        :: Cbeta
  DOUBLE PRECISION                        :: Cibeta
  DOUBLE PRECISION                        :: Cdag
  DOUBLE PRECISION                        :: Cdagbeta
  DOUBLE PRECISION                        :: beta
  DOUBLE PRECISION                        :: ratio
  DOUBLE PRECISION                        :: time
!  TYPE(CdagC)    , POINTER, DIMENSION(:)  :: list => NULL()
#include "BathOperator_hybrid.h"

  this%antiShift = .FALSE.
  beta     = this%beta
  C        =  CdagC_1(C_)
!  Cbeta    = C.MOD.beta
  MODCYCLE(C,beta,Cbeta)
  Cdag     =  CdagC_1(Cdag_)
!  cdagbeta = Cdag.MOD.beta
  MODCYCLE(Cdag,beta,CdagBeta)
!  IF ( Cdag .GE. beta ) &
!    CALL ERROR("BathOperator_getDetAdd : bad case ...              ")
  IF ( this%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperator_getDetAdd : no active hybrid function ")

  tail =  particle%tail
  new_tail = tail+1
!  list => particle%list
  
  IF ( ((C .GT. Cdag) .AND. (position .EQ. -1)) &
       .OR. ((C .LT. Cdag) .AND. (tail .EQ. 0))) THEN ! Possible only if it is a segment
    this%updatePosRow = tail + 1
    this%updatePosCol = tail + 1
  ELSE
    this%updatePosRow  = ABS(position)
    this%updatePosCol  = ABS(position)
  END IF
  
  ! If antisegment, the det ratio has to be by -1 ( sign of the signature of one
  ! permutation line in the this
  IF ( C .LT. Cdag .AND. tail .GT. 0) THEN ! if antiseg
  !  ratio = -ratio 
    this%updatePosRow  = (this%updatePosRow + 1) !position in [1;tail]
    IF ( CdagBeta .LT. particle%list(this%updatePosCol,Cdag_) ) this%antiShift = .TRUE.
  END IF

!  CALL Vector_setSize(this%R,tail)
!  CALL Vector_setSize(this%Q,tail)
  Vector_QuickResize(this%R,new_tail)
  Vector_QuickResize(this%Q,new_tail)
  Vector_QuickResize(this%Rtau,new_tail)
  Vector_QuickResize(this%Qtau,new_tail)

  DO it1 = 1, tail
    it2 = it1 + ( 1+SIGN(1,it1-this%updatePosRow) )/2
    it3 = it1 + ( 1+SIGN(1,it1-this%updatePoscol) )/2

    this%Rtau%vec(it2)= C - particle%list(it1,Cdag_)
    !this%Rtau%vec(it1)= C - particle%list(it1,Cdag_)
    time = Cbeta - particle%list(it1,Cdag_)
#include "BathOperator_hybrid"
    this%R%vec(it1) = hybrid
!    this%R%vec(it) = BathOperator_hybrid(this, Cbeta - list(it)%Cdag)
!    Cibeta = list(it)%C.MOD.beta
    MODCYCLE(particle%list(it1,C_),beta,Cibeta)
    time = Cibeta - Cdagbeta
    this%Qtau%vec(it3)= time
    !this%Qtau%vec(it1)= time
#include "BathOperator_hybrid"
    this%Q%vec(it1) = hybrid
    !this%Q%vec(it3) = hybrid
!    Q(it) = BathOperator_hybrid(this, Cibeta - Cdagbeta)
  END DO
  ! Compute S
  this%Stau = C - Cdagbeta 
  this%Rtau%vec(this%updatePosRow) = this%Stau
  this%Qtau%vec(this%updatePosCol) = this%Rtau%vec(this%updatePosRow)

  time = Cbeta-Cdagbeta
#include "BathOperator_hybrid"
  this%S = hybrid

  !ratio = this%S - DOT_PRODUCT(MATMUL(this%R%vec(1:tail),this%M(this%activeFlavor)%mat(1:tail,1:tail)),this%Q%vec(1:tail))
  ratio = 0.d0
  DO it1 = 1, tail
    time = 0.d0
    DO it2 = 1, tail
      time = time + this%R%vec(it2) * this%M(this%activeFlavor)%mat(it2,it1)
    END DO
    ratio = ratio + this%Q%vec(it1) * time
  END DO
  ratio = this%S - ratio

  this%Stilde = 1.d0 / ratio

  ! This IF is the LAST "NON CORRECTION" in my opinion this should not appears.
!  IF ( MAX(C,Cdag) .GT. this%beta ) THEN
!    WRITE(*,*) this%Stilde
!    this%Stilde = - ABS(this%Stilde)
!  END IF

  ! If antisegment, the det ratio has to be by -1 ( sign of the signature of one
  ! permutation line in the this)
  IF ( C .LT. Cdag .AND. tail .GT. 0) THEN ! if antiseg
    ratio = -ratio 
  ENDIF

  BathOperator_getDetAdd = ratio
  this%MAddFlag   = .TRUE.
!#ifdef CTQMC_CHECK
!  this%ListCdagC = particle
!!write(*,*) this%Stilde
!!write(*,*) this%antishift
!!write(*,*)    this%updatePosRow 
!!write(*,*)    this%updatePosCol 
!#endif

END FUNCTION BathOperator_getDetAdd
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_getDetRemove
!! NAME
!!  BathOperator_getDetRemove
!!
!! FUNCTION
!!  Compute the determinant ratio when a (anti)segment
!!  is trying to be removed 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  position=position of segment to be removed
!!
!! OUTPUT
!!  BathOperator_getDetRemove=the det 
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

DOUBLE PRECISION FUNCTION BathOperator_getDetRemove(this,position)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
!Local arguments-------------------------------
  INTEGER           , INTENT(IN   ) :: position  
  INTEGER                           :: ABSposition  
  INTEGER                           :: tail

  IF ( this%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperator_getDetRemove : no active hybrid fun  ")

  this%antiShift = .FALSE.
  tail         = this%M(this%activeFlavor)%tail
  ABSposition  = ABS(position)
  IF ( ABSposition .GT. tail ) &
    CALL ERROR("BathOperator_getDetRemove : position > M size     ")
  this%updatePosCol = ABSposition
  this%antiShift    = .FALSE.
  IF ( position .GT. 0 ) THEN
    this%updatePosRow = ABSposition
  ELSE
    this%updatePosRow = ABSposition+1
    IF ( ABSposition .EQ. tail ) THEN 
      this%antiShift = .TRUE.
      this%updatePosRow = 1 !ABSposition - 1
!      this%updatePosRow = ABSposition    
!      IF ( this%updatePosCol .EQ. 0) this%updatePosCol = tail
    END IF
  ENDIF
  this%Stilde                 = this%M(this%activeflavor)%mat(this%updatePosRow,this%updatePosCol) 
  this%MRemoveFlag            = .TRUE.
  BathOperator_getDetRemove = this%Stilde

  ! If remove an antiseg , the det ratio has to be multiplied by -1
  IF ( position .LT. 0 .AND. tail .GT. 1 ) &
    BathOperator_getDetRemove = - BathOperator_getDetRemove
!#ifdef CTQMC_CHECK
!  this%ListCdagC = particle
!!write(*,*) this%updatePosRow, this%updatePosCol, position
!!CALL ListCdagC_print(particle)
!#endif

END FUNCTION BathOperator_getDetRemove
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_getDetF
!! NAME
!!  BathOperator_getDetF
!!
!! FUNCTION
!!  Compute the determinant of the F this
!!  using the hybridization of flavor and the 
!!  segments of particle
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  flavor=hybridization function to take
!!  particles=segments to use
!!
!! OUTPUT
!!  BathOperator_getDetF=the det 
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

DOUBLE PRECISION FUNCTION BathOperator_getDetF(this,flavor,particle)

!Arguments ------------------------------------
  TYPE(BathOperator)       , INTENT(INOUT)      :: this
  INTEGER                  , INTENT(IN   )  :: flavor
  TYPE(ListCdagC), OPTIONAL, INTENT(IN   )  :: particle
!Local arguments-------------------------------
  INTEGER :: iCdag
  INTEGER :: iC
  INTEGER :: tail
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: tC
  DOUBLE PRECISION :: tCdag
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: mbeta_two
  DOUBLE PRECISION :: signe
  DOUBLE PRECISION :: inv_dt
#include "BathOperator_hybrid.h"

  BathOperator_getDetF = 1.d0 ! pour eviter des divisions par 0
  IF ( PRESENT( particle ) ) THEN
    tail = particle%tail
    activeF = flavor
    beta = this%beta
    mbeta_two = -beta*0.5d0
    inv_dt =  this%inv_dt
    CALL MatrixHyb_setSize(this%M_update(flavor),tail)
    DO iCdag = 1, tail
      tCdag  = particle%list(iCdag,Cdag_)
      DO iC  = 1, tail
        !tC   = particle%list(C_,iC).MOD.beta
        MODCYCLE(particle%list(iC,C_),beta,tC)
        time = tC - tCdag
#include "BathOperator_hybrid"
        this%M_update(flavor)%mat(iC,iCdag) = hybrid 
      END DO
    END DO
    ! mat_tau needs to be transpose of ordered time mat (way of measuring
    ! G(tau))
    DO iC  = 1, tail
      tC   = particle%list(iC,C_)
      DO iCdag = 1, tail
        tCdag  = particle%list(iCdag,Cdag_)
        time = tC - tCdag
        signe = SIGN(1.d0,time)
        time = time + (signe-1.d0)*mbeta_two
        this%M_update(flavor)%mat_tau(iCdag,iC) = INT( ( time * inv_dt ) + 1.5d0 )
      END DO
    END DO
    CALL MatrixHyb_inverse(this%M_update(flavor),BathOperator_getDetF) ! calcul le det de la matrice et l'inverse
  ELSE
    CALL MatrixHyb_getDet(this%M(flavor),BathOperator_getDetF) ! det M = 1/detF !
    BathOperator_getDetF = 1.d0 / BathOperator_getDetF
  ENDIF
END FUNCTION BathOperator_getDetF
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_setMAdd
!! NAME
!!  BathOperator_setMAdd
!!
!! FUNCTION
!!  Update de M this inserting a row and a column
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  particle=segments of active flavor
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

SUBROUTINE BathOperator_setMAdd(this,particle) 

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
  TYPE(ListCdagC)   , INTENT(IN   ) :: particle
!Local variables ------------------------------
  INTEGER                           :: tail
  INTEGER                           :: new_tail
  INTEGER                           :: col
  INTEGER                           :: col_move
  INTEGER                           :: row_move
  INTEGER                           :: row
  INTEGER                           :: positionRow
  INTEGER                           :: positionCol
  INTEGER                           :: aF
  DOUBLE PRECISION                  :: Stilde
  DOUBLE PRECISION                  :: time
  DOUBLE PRECISION                  :: mbeta_two
  DOUBLE PRECISION                  :: inv_dt
  TYPE(Vector) :: vec_tmp
  TYPE(VectorInt) :: vecI_tmp
  INTEGER :: m
  INTEGER :: count
  INTEGER :: i
  INTEGER :: j
  INTEGER :: p

  IF ( this%MAddFlag .EQV. .FALSE. ) &
    CALL ERROR("BathOperator_setMAdd : MAddFlag turn off           ")
  af = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("BathOperator_setMAdd : no active hybrid function   ")
  tail     =  this%M(aF)%tail
  new_tail =  tail + 1
!CALL this_print(M)

  positionRow =  this%updatePosRow
  positionCol =  this%updatePosCol
  Stilde      =  this%Stilde
!  write(6,*) "before", positionRow, positionCol
  !CALL MatrixHyb_print(this%M(aF),opt_print=1)
  CALL MatrixHyb_setSize(this%M(aF),new_tail)

  ! Compute Qtilde with Q
  !this%Q%vec(1:tail) = (-1.d0) * MATMUL(this%M(aF)%mat(1:tail,1:tail),this%Q%vec(1:tail)) * Stilde
  this%Q%vec(1:tail) = MATMUL(this%M(aF)%mat(1:tail,1:tail),this%Q%vec(1:tail))
  !this%Q%vec(PositionRow:new_tail) = EOSHIFT(this%Q%vec(PositionRow:new_tail), SHIFT=-1, BOUNDARY=-1.d0, DIM=1)
!  this%Qtau%vec(PositionCol:new_tail) = EOSHIFT(this%Qtau%vec(PositionCol:new_tail), SHIFT=-1, BOUNDARY=1.d0, DIM=1)
!  this%Qtau%vec(PositionCol) = this%Stau

  !Compute Rtilde with R and without multiplying by Stilde
  !this%R%vec(1:tail) = (-1.d0) * MATMUL(this%R%vec(1:tail),this%M(aF)%mat(1:tail,1:tail))
  this%R%vec(1:tail) = MATMUL(this%R%vec(1:tail),this%M(aF)%mat(1:tail,1:tail))
  !this%R%vec(PositionCol:new_tail) = EOSHIFT(this%R%vec(PositionCol:new_tail), SHIFT=-1, BOUNDARY=-1.d0, DIM=1)
!  this%Rtau%vec(PositionRow:new_tail) = EOSHIFT(this%Rtau%vec(PositionRow:new_tail), SHIFT=-1, BOUNDARY=1.d0, DIM=1)
!  this%Rtau%vec(PositionRow) = this%Stau

  !Compute the new M this
  !this%M(aF)%mat(PositionRow:new_tail,1:new_tail) = &
  !                   EOSHIFT(this%M(aF)%mat(PositionRow:new_tail,1:new_tail),SHIFT=-1, BOUNDARY=0.d0, DIM=1)
  !this%M(aF)%mat(1:new_tail,PositionCol:new_tail) = &
  !                   EOSHIFT(this%M(aF)%mat(1:new_tail,PositionCol:new_tail),SHIFT=-1, BOUNDARY=0.d0, DIM=2)
! ! this%M(aF)%mat(1:new_tail,1:new_tail) =  this%M(aF)%mat(1:new_tail,1:new_tail) + &
! ! Stilde * MATMUL(RESHAPE(this%Q%vec(1:new_tail),(/ new_tail,1 /)),RESHAPE(this%R%vec(1:new_tail),(/ 1,new_tail /)))

  !this%M(aF)%mat_tau(PositionRow:new_tail,1:new_tail) = &
  !                   EOSHIFT(this%M(aF)%mat_tau(PositionRow:new_tail,1:new_tail),SHIFT=-1, BOUNDARY=0, DIM=1)
  !this%M(aF)%mat_tau(1:new_tail,PositionCol:new_tail) = &
  !                   EOSHIFT(this%M(aF)%mat_tau(1:new_tail,PositionCol:new_tail),SHIFT=-1, BOUNDARY=0, DIM=2)

  mbeta_two = -this%beta*0.5d0
  inv_dt = this%inv_dt
  !Shift mat_tau
  !update old m
  DO col=tail,1,-1
    col_move = col +  ( 1+SIGN(1,col-PositionCol) )/2
    DO row=tail,1,-1
      row_move = row +  ( 1+SIGN(1,row-PositionRow) )/2
      this%M(aF)%mat_tau(row_move,col_move) = this%M(aF)%mat_tau(row,col)
      this%M(aF)%mat(row_move,col_move) = this%M(aF)%mat(row,col) + this%Q%vec(row)*this%R%vec(col) * Stilde
    END DO
  END DO
  ! Add new stuff for new row
  DO row = 1, tail
    row_move = row +  ( 1+SIGN(1,row-PositionRow) )/2
    this%M(aF)%mat(row_move,PositionCol) = -this%Q%vec(row)*Stilde
    time = this%Rtau%vec(row)
    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
    this%M(aF)%mat_tau(row,PositionCol) = INT ( (time*inv_dt) +1.5d0 )
  END DO
  ! Add last time missing in the loops
  time = this%Rtau%vec(new_tail)
  time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  this%M(aF)%mat_tau(new_tail,PositionCol) = INT ( (time*inv_dt) +1.5d0 )
  ! Add new stuff for new col
  DO col = 1, tail 
    col_move = col +  ( 1+SIGN(1,col-PositionCol) )/2
    this%M(aF)%mat(PositionRow,col_move) = -this%R%vec(col)*Stilde
    time = this%Qtau%vec(col)
    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
    this%M(aF)%mat_tau(PositionRow,col) = INT ( (time*inv_dt) +1.5d0 )
  END DO
  ! Add last time missing in the loops
  time = this%Qtau%vec(new_tail)
  time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  this%M(aF)%mat_tau(PositionRow,new_tail) = INT ( (time*inv_dt) +1.5d0 )

  this%M(aF)%mat(PositionRow,PositionCol) = Stilde

  !CALL MatrixHyb_print(this%M(aF),opt_print=1)

!  DO col = 1, new_tail
!    time = this%Rtau%vec(col)
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!    this%M(aF)%mat_tau(col,PositionCol) = INT ( (time*inv_dt) +1.5d0 )
!    time = this%Qtau%vec(col)
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!    this%M(aF)%mat_tau(PositionRow,Col) = INT ( (time*inv_dt) +1.5d0 )
!    time = this%R%vec(col)*Stilde
!    DO row = 1, new_tail
!      this%M(aF)%mat(row,col) = this%M(aF)%mat(row,col) + this%Q%vec(row)*time
!    END DO
!  END DO

  !col_move = new_tail
  !col      = tail
  !DO col_move = new_tail, 1, -1
  !  IF ( col_move .EQ. positionCol ) THEN
  !    ! on calcule rajoute Q tilde
  !    !row_move = new_tail
  !    row      = tail 
  !    DO row_move = new_tail, 1, -1
  !      ! calcul itau
  !      IF ( row_move .EQ. positionRow ) THEN
  !        this%M(aF)%mat(row_move,col_move) = Stilde
  !        !time = this%Stau
  !      ELSE
  !        this%M(aF)%mat(row_move,col_move) = -this%Q%vec(row)*Stilde
  !        !time = this%Rtau%vec(row_move)
  !        row      = row      - 1 
  !      END IF
  !      !time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  !      !this%M(aF)%mat_tau(row_move,col_move) = INT ( (time*inv_dt) +1.5d0 )
  !    END DO
  !    ! realignement des indices
  !  ELSE
  !    ! on calcule Ptilde
  !    !row_move = new_tail
  !    row      = tail 
  !    DO row_move = new_tail, 1, -1
  !      IF ( row_move .EQ. positionRow ) THEN
  !        this%M(aF)%mat(row_move,col_move) = -this%R%vec(col) * Stilde
  !        ! calcul itau
  !        !time = this%Qtau%vec(col_move)
  !        !time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  !        !this%M(aF)%mat_tau(row_move,col_move) = INT ( (time*inv_dt) +1.5d0 )
  !      ELSE
  !        this%M(aF)%mat(row_move,col_move) = this%M(aF)%mat(row,col) + this%Q%vec(row)*this%R%vec(col)*Stilde
  !        ! copy itau
  !        !this%M(aF)%mat_tau(row_move,col_move) = this%M(aF)%mat_tau(row,col)
  !        row      = row      - 1 
  !      END IF
  !    END DO
  !    col      = col      - 1
  !  END IF
  !END DO
!  write(6,*) "after"
!  CALL MatrixHyb_print(this%M(aF),opt_print=1)
!CALL this_inverse(M)
!CALL MatrixHyb_print(M)
!CALL this_inverse(M)

  IF ( this%antiShift .EQV. .TRUE. ) THEN ! antisegment
    CALL Vector_init(vec_tmp,new_tail)
    CALL VectorInt_init(vecI_tmp,new_tail)
  ! Shift if necessary according to this%antishift
  ! shift DIM=2 (col)
    p = new_tail - 1
    m = 1
    count = 0
    DO WHILE ( count .NE. new_tail )
      vec_tmp%vec(1:new_tail) = this%M(aF)%mat(1:new_tail,m)
      vecI_tmp%vec(1:new_tail) = this%M(aF)%mat_tau(1:new_tail,m)
      i = m
      !j = m+p
      MODCYCLE(m+p, new_tail, j)
      DO WHILE (j .NE. m)
        this%M(aF)%mat(1:new_tail,i) = this%M(aF)%mat(1:new_tail,j)
        this%M(aF)%mat_tau(1:new_tail,i) = this%M(aF)%mat_tau(1:new_tail,j)
        i = j
        MODCYCLE(j+p, new_tail, j)
        count = count+1
      END DO
      this%M(aF)%mat(1:new_tail,i) = vec_tmp%vec(1:new_tail)
      this%M(aF)%mat_tau(1:new_tail,i) = vecI_tmp%vec(1:new_tail)
      count = count+1
      m = m+1
    END DO
    ! shift DIM=1 (row)
    p = new_tail - 1
    m = 1
    count = 0
    DO WHILE ( count .NE. new_tail)
      vec_tmp%vec(1:new_tail) = this%M(aF)%mat(m,1:new_tail)
      vecI_tmp%vec(1:new_tail) = this%M(aF)%mat_tau(m,1:new_tail)
      i = m
      !j = m+p
      MODCYCLE(m+p, new_tail, j)
      DO WHILE ( j .NE. m )
        this%M(aF)%mat(i,1:new_tail) = this%M(aF)%mat(j,1:new_tail)
        this%M(aF)%mat_tau(i,1:new_tail) = this%M(aF)%mat_tau(j,1:new_tail)
        i = j
        MODCYCLE(j+p, new_tail, j)
        count = count+1
      END DO
      this%M(aF)%mat(i,1:new_tail) = vec_tmp%vec(1:new_tail)
      this%M(aF)%mat_tau(i,1:new_tail) = vecI_tmp%vec(1:new_tail)
      count = count+1
      m = m+1
    END DO
    CALL Vector_destroy(vec_tmp)
    CALL VectorInt_destroy(vecI_tmp)
    !this%M(aF)%mat(1:new_tail,1:new_tail) = CSHIFT(this%M(aF)%mat(1:new_tail,1:new_tail), SHIFT=-1, DIM=1) ! Shift to the bottom
    !this%M(aF)%mat(1:new_tail,1:new_tail) = CSHIFT(this%M(aF)%mat(1:new_tail,1:new_tail), SHIFT=-1, DIM=2) ! Shift to the right
    !this%M(aF)%mat_tau(1:new_tail,1:new_tail) = CSHIFT(this%M(aF)%mat_tau(1:new_tail,1:new_tail), SHIFT=-1, DIM=1) ! Shift to the bottom
    !this%M(aF)%mat_tau(1:new_tail,1:new_tail) = CSHIFT(this%M(aF)%mat_tau(1:new_tail,1:new_tail), SHIFT=-1, DIM=2) ! Shift to the right
!CALL this_print(M)
  END IF

  IF ( this%doCheck .EQV. .TRUE.) THEN
!#ifdef CTQMC_CHECK
  CALL BathOperator_checkM(this,particle)
!#endif
  END IF

  this%MAddFlag = .FALSE.

END SUBROUTINE BathOperator_setMAdd
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_setMRemove
!! NAME
!!  BathOperator_setMRemove
!!
!! FUNCTION
!!  delete one row and one column of the M this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  particle=segments of the active flavor
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

SUBROUTINE BathOperator_setMRemove(this,particle) 

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT)  :: this
  TYPE(ListCdagC)   , INTENT(IN   )  :: particle
!Local variables ------------------------------
  INTEGER                            :: tail
  INTEGER                            :: new_tail
  INTEGER                            :: col
  INTEGER                            :: col_move
  INTEGER                            :: row_move
  INTEGER                            :: row
  INTEGER                            :: positionCol
  INTEGER                            :: positionRow
  INTEGER                            :: aF
  INTEGER                              :: m
  INTEGER                              :: count
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: p
  DOUBLE PRECISION                   :: invStilde
  DOUBLE PRECISION                   :: invStilde2
  TYPE(VectorInt) :: vecI_tmp
  TYPE(Vector)    :: vec_tmp

  IF ( this%MRemoveFlag .EQV. .FALSE. ) &
    CALL ERROR("BathOperator_setMRemove : MRemoveFlag turn off     ")
  af = this%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("BathOperator_setMRemove : no active hybrid func    ")
  tail        =  this%M(aF)%tail
  new_tail    =  tail - 1
  positionCol =  this%updatePosCol
  positionRow =  this%updatePosRow
  invStilde   = 1.d0 / this%Stilde

!  write(6,*) "before", positionRow, positionCol
!  CALL MatrixHyb_print(this%M(aF),opt_print=1)

!  IF ( new_tail .EQ. 0 ) THEN
!!    IF ( this%antiShift .EQV. .TRUE.  ) THEN
!!      this%M(aF)%mat(1,1) = 1.d0/BathOperator_Hybrid(this, this%beta)
!!      this%MRemoveFlag = .FALSE.
!!      RETURN
!!    END IF
!    CALL MatrixHyb_clear(this%M(aF))
!    this%MRemoveFlag = .FALSE.
!    RETURN
!  END IF

!  CALL Vector_setSize(this%Q,new_tail)
!  CALL Vector_setSize(this%R,new_tail)
  Vector_QuickResize(this%Q,new_tail)
  Vector_QuickResize(this%R,new_tail)

!  We use R and Q as this%R%vec and this%Q%vec
!  this%R%vec => this%R
!  this%Q%vec => this%Q

  !row      = 1
  !row_move = 1
  !col      = 1
  !col_move = 1
  DO row_move = 1, new_tail
    !IF ( row .EQ. positionRow ) row = row + 1
    !IF ( col .EQ. positionCol ) col = col + 1
    col = row_move + (1+SIGN(1,row_move-positionCol))/2
    row = row_move + (1+SIGN(1,row_move-positionRow))/2
    this%R%vec(row_move) = this%M(aF)%mat(positionRow,col)
    this%Q%vec(row_move) = this%M(aF)%mat(row,positionCol)
    !row      = row + 1 
    !col      = col + 1
  END DO
!!    this%R%vec(1:positionCol-1) = this%M(aF)%mat(positionRow,1:positionCol-1)
!!    this%R%vec(positionCol:new_tail) = this%M(aF)%mat(positionRow,positionCol+1:tail)
!!    this%Q%vec(1:positionRow-1) = this%M(aF)%mat(1:positionRow-1,positionCol)
!!    this%Q%vec(positionRow:new_tail) = this%M(aF)%mat(positionRow+1:tail,positionCol)
!write(*,*) positionRow, positionCol
!CALL MatrixHyb_print(M)
!CALL Vector_print(this%R)
!CALL Vector_print(this%Q)
!CALL ListCdagC_print(this%ListCdagC)

  !col      = 1
  DO col_move = 1, new_tail 
    !IF ( col_move .EQ. positionCol ) col = col + 1
    col = col_move + (1+SIGN(1,col_move-positionCol))/2
    !row      = 1
    invStilde2 = invStilde * this%R%vec(col_move)
    DO row_move = 1, new_tail
      !IF ( row_move .EQ. positionRow ) row = row + 1
      row = row_move + (1+SIGN(1,row_move-positionRow))/2
      this%M(aF)%mat(row_move,col_move) = this%M(aF)%mat(row,col) &
                                      - this%Q%vec(row_move)*invStilde2
      this%M(aF)%mat_tau(row_move,col_move) = this%M(aF)%mat_tau(row,col)
      !row      = row      + 1
    END DO
    !col      = col      + 1 
  END DO
  CALL MatrixHyb_setSize(this%M(aF),new_tail)

  IF ( this%antiShift .EQV. .TRUE. ) THEN ! antisegment
    ! Shift if necessary according to this%antishift
    ! shift DIM=2 (col)
    CALL Vector_init(vec_tmp,new_tail)
    CALL VectorInt_init(vecI_tmp,new_tail)
    p = 1
    m = 1
    count = 0
    DO WHILE ( count .NE. new_tail )
      vec_tmp%vec(1:new_tail) = this%M(aF)%mat(1:new_tail,m)
      vecI_tmp%vec(1:new_tail) = this%M(aF)%mat_tau(1:new_tail,m)
      i = m
      !j = m+p
      MODCYCLE(m+p, new_tail, j)
      DO WHILE (j .NE. m)
        this%M(aF)%mat(1:new_tail,i) = this%M(aF)%mat(1:new_tail,j)
        this%M(aF)%mat_tau(1:new_tail,i) = this%M(aF)%mat_tau(1:new_tail,j)
        i = j
        MODCYCLE(j+p, new_tail, j)
        count = count+1
      END DO
      this%M(aF)%mat(1:new_tail,i) = vec_tmp%vec(1:new_tail)
      this%M(aF)%mat_tau(1:new_tail,i) = vecI_tmp%vec(1:new_tail)
      count = count+1
      m = m+1
    END DO
    CALL Vector_destroy(vec_tmp)
    CALL VectorInt_destroy(vecI_tmp)
    !this%M(aF)%mat(1:new_tail,1:new_tail) = &
    !           CSHIFT(this%M(aF)%mat(1:new_tail,1:new_tail), SHIFT=1, DIM=2) ! Shift to the top
    !this%M(aF)%mat_tau(1:new_tail,1:new_tail) = &
    !           CSHIFT(this%M(aF)%mat_tau(1:new_tail,1:new_tail), SHIFT=1, DIM=2) ! Shift to the top
  END IF
!  write(6,*) "after "
!  CALL MatrixHyb_print(this%M(aF),opt_print=1)

  IF ( this%doCheck .EQV. .TRUE. ) THEN
!#ifdef CTQMC_CHECK
  CALL BathOperator_checkM(this,particle)
!#endif
  END IF

  this%MRemoveFlag = .FALSE.

END SUBROUTINE BathOperator_setMRemove
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_swap
!! NAME
!!  BathOperator_swap
!!
!! FUNCTION
!!  Recompute 2 M this swaping the segments
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  iflavor1=flavor to swap with the next one
!!  iflavor2=favor to swap with the previous one
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

SUBROUTINE BathOperator_swap(this, flavor1, flavor2)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
  INTEGER           , INTENT(IN   ) :: flavor1
  INTEGER           , INTENT(IN   ) :: flavor2

  !CALL MatrixHyb_print(this%M(flavor1),234)
  this%M(flavor1) = this%M_update(flavor1)
  !CALL MatrixHyb_print(this%M(flavor1),234)
  !CALL MatrixHyb_print(this%M(flavor2),234)
  this%M(flavor2) = this%M_update(flavor2)
  !CALL MatrixHyb_print(this%M(flavor2),234)

END SUBROUTINE BathOperator_swap
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_initF
!! NAME
!!  BathOperator_initF
!!
!! FUNCTION
!!  Copy input hybridization functions from a file
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  ifstream=file stream to read F
!!
!! OUTPUT
!!  argout(sizeout)=description
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

SUBROUTINE BathOperator_initF(this,ifstream)

!Arguments ----------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
  INTEGER           , INTENT(IN   ) :: ifstream
!Local variables ----------------
  INTEGER                           :: flavor
  INTEGER                           :: sample

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("BathOperator_initF : BathOperator not set         ")

  DO flavor=1,this%flavors
    DO sample = 1, this%sizeHybrid
      READ(ifstream,*) this%F(sample,flavor)
    END DO
  END DO
END SUBROUTINE BathOperator_initF
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_setF
!! NAME
!!  BathOperator_setF
!!
!! FUNCTION
!!  Copy F from input array
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  F=array of the hybridization function
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

SUBROUTINE BathOperator_setF(this,F)

!Arguments ------------------------------------
  TYPE(BathOperator)               , INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:,:) , INTENT(IN   ) :: F
!Arguments ------------------------------------
  INTEGER                                          :: flavor
  INTEGER                                          :: sample
  INTEGER                                          :: length

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("BathOperator_setF : BathOperator not set          ")

 length  = SIZE(F)
  IF ( length .NE. (this%flavors * this%sizeHybrid) ) &
    CALL ERROR("BathOperator_setF : wrong input F                 ")

  DO flavor=1,this%flavors
    DO sample = 1, this%sizeHybrid
    this%F(sample,flavor) = F(sample,flavor)
    END DO
  END DO
END SUBROUTINE BathOperator_setF
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_printF
!! NAME
!!  BathOperator_printF
!!
!! FUNCTION
!!  print F function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  ostream=file stream to write in
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

SUBROUTINE BathOperator_printF(this,ostream)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(INOUT) :: this
  INTEGER,OPTIONAL  , INTENT(IN   ) :: ostream
!Local variables ------------------------------
  CHARACTER(LEN=4)                  :: aflavor
  CHARACTER(LEN=50)                  :: string
  INTEGER                           :: flavor
  INTEGER                           :: sample
  INTEGER                           :: ostream_val

  IF ( PRESENT(ostream) ) THEN 
    ostream_val = ostream
  ELSE  
    ostream_val = 65
    OPEN(UNIT=ostream_val, FILE="F.dat")
  END IF

  WRITE(aflavor,'(I4)') this%flavors+1
  string = '(1x,'//TRIM(ADJUSTL(aflavor))//'E22.14)'
  DO sample = 1, this%sizeHybrid
    WRITE(ostream_val,string) (sample-1)*this%dt, (this%F(sample,flavor), flavor=1,this%flavors)
  END DO
  !CALL FLUSH(ostream_val)

  IF ( .NOT. PRESENT(ostream) ) &
    CLOSE(ostream_val)

END SUBROUTINE BathOperator_printF
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_printM
!! NAME
!!  BathOperator_printM
!!
!! FUNCTION
!!  print M =F^{-1} this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  ostream=file stream to write in
!!
!! OUTPUT
!!  argout(sizeout)=description
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

SUBROUTINE BathOperator_printM(this,ostream)

!Arguments ------------------------------------
  TYPE(BathOperator), INTENT(IN) :: this
  INTEGER, OPTIONAL , INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                        :: ostream_val

  IF ( this%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperator_printM : no active hybrid function    ")
  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  CALL MatrixHyb_print(this%M(this%activeFlavor),ostream_val)
END SUBROUTINE BathOperator_printM
!!***

!!****f* ABINIT/m_BathOperator/ BathOperator_destroy
!! NAME
!!   BathOperator_destroy
!!
!! FUNCTION
!!  Deallocate and reset every thing
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
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

SUBROUTINE  BathOperator_destroy(this)

  TYPE(BathOperator), INTENT(INOUT) :: this
  INTEGER  :: it

  DO it = 1, this%flavors
    CALL MatrixHyb_destroy(this%M(it))
    CALL MatrixHyb_destroy(this%M_update(it))
  END DO

  CALL Vector_destroy(this%R)
  CALL Vector_destroy(this%Q)
  CALL Vector_destroy(this%Rtau)
  CALL Vector_destroy(this%Qtau)
  FREEIF(this%F)
  DT_FREEIF(this%M)
  DT_FREEIF(this%M_update)

  this%MAddFlag     = .FALSE.
  this%MRemoveFlag  = .FALSE.
  this%flavors      = 0 
  this%beta         = 0.d0
  this%dt      = 0.d0
  this%inv_dt  = 0.d0
  this%samples      = 0
  this%sizeHybrid   = 0
  this%activeFlavor = 0 
  this%updatePosRow = 0
  this%updatePosCol = 0

END SUBROUTINE BathOperator_destroy
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_doCheck
!! NAME
!!  BathOperator_doCheck
!!
!! FUNCTION
!!  Just store if we perfom check for updates of M
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  opt_check=second bit should be one
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

SUBROUTINE BathOperator_doCheck(this,opt_check)

!Arguments ------------------------------------
  TYPE(BathOperator) , INTENT(INOUT) :: this
  INTEGER            , INTENT(IN   ) :: opt_check
  
  IF ( opt_check .GE. 2 ) &
    this%doCheck = .TRUE.
END SUBROUTINE BathOperator_doCheck
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_checkM
!! NAME
!!  BathOperator_checkM
!!
!! FUNCTION
!!  compute from scratch the M this and compar it
!!  with the already computed M this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!  particle=list of all segments of the active flavor
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

SUBROUTINE BathOperator_checkM(this,particle)

!Arguments ------------------------------------
  TYPE(BathOperator) , INTENT(INOUT) :: this
  TYPE(ListCdagC)    , INTENT(IN   ) :: particle
!Local variables ------------------------------
!  TYPE(MatrixHyb)                    :: checkMatrix
  LOGICAL :: checkTau
  INTEGER :: tail
  INTEGER :: iC
  INTEGER :: iCdag
  INTEGER :: aF
  CHARACTER(LEN=4) :: a
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: mbeta_two
  DOUBLE PRECISION :: erreur
  DOUBLE PRECISION :: tc
  DOUBLE PRECISION :: tCdag
  DOUBLE PRECISION :: sumMmat
  DOUBLE PRECISION :: sumCheck
#include "BathOperator_hybrid.h"

  aF = this%activeFlavor
  !Construction de la this
  tail = particle%tail
!  CALL MatrixHyb_init(checkMatrix,this%iTech,size=tail,Wmax=this%samples)
!  CALL MatrixHyb_setSize(checkMatrix,tail)
  CALL MatrixHyb_setSize(this%M_update(aF),tail)
  beta   =  this%beta
  mbeta_two = -beta*0.5d0
  this%checkNumber = this%checkNumber + 1
  IF ( tail .NE. this%M(aF)%tail ) THEN
    CALL WARN("BathOperator_checkM : tails are different          ")
    RETURN
  END IF

!CALL ListCdagC_print(particle)
  DO iCdag = 1, tail
    tCdag  = particle%list(iCdag,Cdag_)
    DO iC  = 1, tail
      !tC   = particle%list(C_,iC).MOD.beta
      MODCYCLE(particle%list(iC,C_),beta,tC)
      time = tC - tCdag
#include "BathOperator_hybrid"
      this%M_update(aF)%mat(iC,iCdag) = hybrid

      time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
      this%M_update(aF)%mat_tau(iCdag,iC) = INT ( (time*this%inv_dt) +1.5d0 ) 
    END DO
  END DO

!    CALL MatrixHyb_Print(checkMatrix)
  !Inversion de la this
  CALL MatrixHyb_inverse(this%M_update(aF))
!    CALL MatrixHyb_Print(checkMatrix)

  !Comparaison
  sumMmat =0.d0
  sumCheck=0.d0
  erreur = 0.d0
  checkTau = .FALSE.
  DO iCdag = 1, tail
    Do iC =1, tail
      this%M_update(aF)%mat(iC,iCdag) = ABS((this%M_update(aF)%mat(iC, iCdag) - this%M(aF)%mat(iC,iCdag))/this%M(aF)%mat(iC,iCdag))
      IF ( this%M_update(aF)%mat(iC,iCdag) .GT. erreur ) erreur = this%M_update(aF)%mat(ic,iCdag)
      IF ( this%M_update(aF)%mat_tau(iC,iCdag) .NE. this%M(aF)%mat_tau(iC,iCdag) ) checkTau = .TRUE.
    END DO
  END DO

  IF ( checkTau .EQV. .TRUE. ) THEN
    CALL WARN("BathOperator_checkM : mat_tau differs should be")
    CALL MatrixHyb_print(this%M_update(aF),opt_print=1)
    CALL WARN("BathOperator_checkM : whereas it is")
    CALL MatrixHyb_print(this%M(aF),opt_print=1)
  END IF
  this%meanError = this%meanError + erreur
  IF ( erreur .GT. 1.d0 ) THEN 
    WRITE(a,'(I4)') INT(erreur*100.d0)
!    CALL MatrixHyb_Print(this%M(aF)
    CALL WARN("BathOperator_checkM : "//a//"%                        ") 
  END IF
!  CALL MatrixHyb_destroy(checkMatrix)
END SUBROUTINE BathOperator_checkM
!!***

!!****f* ABINIT/m_BathOperator/BathOperator_getError
!! NAME
!!  BathOperator_getError
!!
!! FUNCTION
!!  compute a percentage error / checkM
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=bath operator
!!
!! OUTPUT
!!  BathOperator_getError=Error in percent
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

DOUBLE PRECISION FUNCTION BathOperator_getError(this)

  TYPE(BathOperator), INTENT(IN) :: this

  IF ( this%doCheck .EQV. .TRUE. ) THEN
    BathOperator_getError = this%meanError / DBLE(this%checkNumber)
  ELSE
    BathOperator_getError = 0.d0
  END IF
END FUNCTION BathOperator_getError
!!***
!#endif

END MODULE m_BathOperator
!!***

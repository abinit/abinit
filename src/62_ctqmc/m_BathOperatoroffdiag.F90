
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!****m* ABINIT/m_BathOperatoroffdiag
!! NAME
!!  m_BathOperatoroffdiag
!! 
!! FUNCTION 
!!  Manage all stuff related to the bath for the 
!!  simgle Anderson Impurity Model
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder, J. Denier, B. Amadon)
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
MODULE m_BathOperatoroffdiag
USE m_MatrixHyb
USE m_Vector
USE m_VectorInt
USE m_Global
USE m_ListCdagC
IMPLICIT NONE

! subroutines
 public :: BathOperatoroffdiag_init
 public :: BathOperatoroffdiag_reset
 public :: BathOperatoroffdiag_activateParticle
 public :: BathOperatoroffdiag_setMAdd
 public :: BathOperatoroffdiag_setMRemove
 public :: BathOperatoroffdiag_swap
 public :: BathOperatoroffdiag_initF
 public :: BathOperatoroffdiag_setF
 public :: BathOperatoroffdiag_printF
 public :: BathOperatoroffdiag_printM
 public :: BathOperatoroffdiag_destroy
 public :: BathOperatoroffdiag_doCheck
 public :: BathOperatoroffdiag_checkM

! functions
! public :: BathOperatoroffdiag_hybrid
 public :: BathOperatoroffdiag_getDetAdd
 public :: BathOperatoroffdiag_getDetRemove
 public :: BathOperatoroffdiag_getDetF
 public :: BathOperatoroffdiag_getError
!!***

!!****t* m_BathOperatoroffdiag/BathOperatoroffdiag
!! NAME
!!  BathOperatoroffdiag
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE BathOperatoroffdiag
  LOGICAL :: set         = .FALSE.
  ! True if the BathOperatoroffdiag is initialized in BathOperatoroffdiag_init

  LOGICAL :: MAddFlag    = .FALSE.
  ! Set to true if we can compute a new M (see updateDetXX) (ie in
  ! BathOperatoroffdiag_getDetAdd)

  LOGICAL :: MRemoveFlag = .FALSE. 
  ! Set to true if we can compute a new M (see updateDetXX) (ie in
  ! BathOperatoroffdiag_getDetRemove)

  LOGICAL :: antiShift   = .FALSE. 
  ! shift when M is updated with antiseg

  LOGICAL :: doCheck     = .FALSE.
  ! TRUE is checks are activated

  INTEGER :: opt_nondiag = 0
! if opt_nondiag = 1 F is non diagonal.

  INTEGER :: flavors
! number of flavors
! if opt_nondiag = 0 , flavors= number of flavor
! if opt_nondiag = 1 , flavors= 1

  INTEGER :: activeFlavor
  ! Active flavor on which a segment is added/suppressed...

  INTEGER :: samples
  ! Number of time slices (given in the input file)

  INTEGER :: sizeHybrid
  ! Number of time slices (given in the input file) + 1 (=qmc_l+1)

  INTEGER :: updatePosRow
  ! Gives the position of new Row to add
  ! Modified in  BathOperatoroffdiag_getDetAdd and  BathOperatoroffdiag_getDetRemove
  ! could be the Row in the Full matrix for the non diag implementation

  INTEGER :: updatePosCol
  ! Gives the position of new Col to add
  ! Modified in  BathOperatoroffdiag_getDetAdd and  BathOperatoroffdiag_getDetRemove

  INTEGER :: iTech
  ! iTech is an integer which precise the technics used to compute the
  ! Green's function (in time or frequency)

  INTEGER :: sumtails
  !  size of the full F matrix (sums of tails(iflavor) over iflavor)

  INTEGER,          ALLOCATABLE, DIMENSION(:) :: tails
  ! tails(iflavor) is the current number of segments for the flavor iflavor

  INTEGER,          ALLOCATABLE, DIMENSION(:) :: Fshift
  ! Fshift(iflavor) is the sum of number of segments for all flavors iflavor2
  ! such that iflavor< iflavor
  ! It is thus the shift in the F matrix to have the first segment of the flavor
  ! iflavor
  ! Fshift(nflavor+1) is the total nb of tails (=sumtails) 

  DOUBLE PRECISION                            :: beta
  ! Inverse of Temperature
  ! 

  DOUBLE PRECISION                            :: dt
  ! dt=beta/samples

  DOUBLE PRECISION                            :: inv_dt
  ! inv_dt=1/dt

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: F ! qmc_l+2,Flavors
  ! Hybridization function F(1:op%sizeHybrid+1,1:flavors,1:flavors)

  DOUBLE PRECISION                            :: S
  ! Sherman Morrison notations 

  DOUBLE PRECISION                            :: Stau
  ! Sherman Morrison notations 

  DOUBLE PRECISION                            :: Stilde
  ! Sherman Morrison notations 

  TYPE(Vector)                                :: R 
  ! Sherman Morrison notations R%vec(size).
  ! computed for each flavor (As matrices are made of Blocks for each
  ! flavor because the code is restricted to diagonal F matrices)

  TYPE(Vector)                                :: Q 
  ! Sherman Morrison notations 
  ! computed for each flavor (As matrices are made of Blocks for each
  ! flavor because the code is restricted to diagonal F matrices)

  TYPE(Vector)                                :: Rtau
  ! Sherman Morrison notations 
  ! Rtau gives the time length for each elements of R
  ! computed for each flavor (As matrices are made of Blocks for each
  ! flavor because the code is restricted to diagonal F matrices)

  TYPE(Vector)                                :: Qtau
  ! Sherman Morrison notations 
  ! Qtau gives the time length for each elements of Q
  ! computed for each flavor (As matrices are made of Blocks for each
  ! flavor because the code is restricted to diagonal F matrices)

  TYPE(MatrixHyb)                             :: M  ! Flavors
  ! inverse of  Hybridization matrix  M%mat(global_size,global_size)
  ! contains the value of the hybridization for all flavor and segments times, 
  ! the times (mat_tau), and possibly the
  ! frequency

  TYPE(MatrixHyb)                             :: M_update  ! Flavors
  !  used in BathOperatoroffdiag_getdetF and in BathOperatoroffdiag_checkM
  ! for checks 

!#ifdef CTQMC_CHECK
  INTEGER                                     :: checkNumber
  DOUBLE PRECISION                            :: meanError
!  TYPE(ListCdagC)                             :: ListCdagC 
!#endif
END TYPE BathOperatoroffdiag
!!***

CONTAINS
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_init
!! NAME
!!  BathOperatoroffdiag_init
!!
!! FUNCTION
!!  Initialize and allocate data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath object
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

SUBROUTINE BathOperatoroffdiag_init(op, flavors, samples, beta, iTech,opt_nondiag)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER           , INTENT(IN   ) :: flavors
  INTEGER           , INTENT(IN   ) :: samples
  INTEGER           , INTENT(IN   ) :: opt_nondiag
  DOUBLE PRECISION  , INTENT(IN   ) :: beta
!Local variables ------------------------------
  INTEGER           , INTENT(IN   ) :: iTech
  !INTEGER                           :: it

  op%MAddFlag     = .FALSE.
  op%MRemoveFlag  = .FALSE.
  op%flavors      = flavors
  op%opt_nondiag  = opt_nondiag
  op%beta         = beta
  op%samples      = samples
  op%sizeHybrid   = samples + 1
  op%dt      = beta / DBLE(samples)
  op%inv_dt  = DBLE(samples) / beta
  op%activeFlavor= 0 
  op%updatePosRow = 0
  op%updatePosCol = 0
  op%iTech        = iTech
!#ifdef CTQMC_CHECK
  op%checkNumber  = 0
  op%meanError    = 0.d0
  op%doCheck = .FALSE.
!#endif

  FREEIF(op%F)
  MALLOC(op%F,(1:op%sizeHybrid+1,1:flavors,1:flavors))
  DT_FREEIF(op%tails)
  DT_MALLOC(op%tails,(1:op%flavors))
  op%tails=0
  DT_FREEIF(op%Fshift)
  DT_MALLOC(op%Fshift,(1:op%flavors+1))
  op%Fshift=0
  
  CALL Vector_init(op%R,100*op%flavors)
  CALL Vector_init(op%Q,100*op%flavors)
  CALL Vector_init(op%Rtau,100*op%flavors)
  CALL Vector_init(op%Qtau,100*op%flavors)

  CALL MatrixHyb_init(op%M,op%iTech,size=Global_SIZE*op%flavors,Wmax=samples) !FIXME Should be consistent with ListCagC
  CALL MatrixHyb_init(op%M_update,op%iTech,size=Global_SIZE*op%flavors,Wmax=samples) !FIXME Should be consistent with ListCagC
  op%F       = 0.d0
  op%set     = .TRUE.
  
END SUBROUTINE BathOperatoroffdiag_init
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_reset
!! NAME
!!  BathOperatoroffdiag_reset
!!
!! FUNCTION
!!  Reset all internal variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator to reset
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

SUBROUTINE BathOperatoroffdiag_reset(op)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER                           :: iflavor
  op%MAddFlag     = .FALSE.
  op%MRemoveFlag  = .FALSE.
  op%activeFlavor = 0 
  op%updatePosRow = 0
  op%updatePosCol = 0
!#ifdef CTQMC_CHECK
  op%checkNumber  = 0
  op%meanError    = 0.d0
  op%sumtails    = 0
!#endif
  op%doCheck = .FALSE.
  CALL Vector_clear(op%R)
  CALL Vector_clear(op%Q)
  CALL Vector_clear(op%Rtau)
  CALL Vector_clear(op%Qtau)

  CALL MatrixHyb_clear(op%M) !FIXME Should be consistent with ListCagC
  op%F       = 0.d0
  do iflavor=1,op%flavors
    op%tails(iflavor)=0
    op%Fshift(iflavor)=0
  enddo

END SUBROUTINE BathOperatoroffdiag_reset
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_activateParticle
!! NAME
!!  BathOperatoroffdiag_activateParticle
!!
!! FUNCTION
!!  Just save on wicht flavor we are working
!!  It is better to use the macro defined in defs.h
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_activateParticle(op,flavor)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER           , INTENT(IN   ) :: flavor

  IF ( flavor .GT. op%flavors ) &
    CALL ERROR("BathOperatoroffdiag_activateParticle : out of range      ")
  IF ( op%set .EQV. .TRUE. ) THEN 
    op%activeFlavor =  flavor
    op%MAddFlag     = .FALSE.
    op%MRemoveFlag  = .FALSE.
  ELSE
    CALL ERROR("BathOperatoroffdiag_activateParticle : not allocated      ")
  END IF
END SUBROUTINE BathOperatoroffdiag_activateParticle
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_getDetAdd
!! NAME
!!  BathOperatoroffdiag_getDetAdd
!!
!! FUNCTION
!!  Compute the determinant ratio when a (anti)segment
!!  is trying to be added and store some array for setMAdd
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
!!  CdagC_1=segment to be added
!!  position=ordered position of the Cdag time
!!  particle=full list of CdagC for activeFlavor
!!
!! OUTPUT
!!  BathOperatoroffdiag_getDetAdd=the det 
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
DOUBLE PRECISION  FUNCTION BathOperatoroffdiag_getDetAdd(op,CdagC_1, position, particle)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag)      , INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(1:2), INTENT(IN   ) :: CdagC_1
  INTEGER                 , INTENT(IN   ) :: position  
  TYPE(ListCdagC), INTENT(IN   ) :: particle(:)
!Local variables-------------------------------
  INTEGER                                 :: it1
  INTEGER                                 :: it2
  INTEGER                                 :: it3,iflavor,iflavora,iflavorb
  INTEGER                                 :: iflavorbegin,iflavorend
  INTEGER                                 :: tail,tailbegin,tailend
  INTEGER                                 :: tcheck
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
#include "BathOperatoroffdiag_hybrid.h"

  op%antiShift = .FALSE.
  beta     = op%beta
  C        =  CdagC_1(C_)
!  Cbeta    = C.MOD.beta
  MODCYCLE(C,beta,Cbeta)
  Cdag     =  CdagC_1(Cdag_)
!  cdagbeta = Cdag.MOD.beta
  MODCYCLE(Cdag,beta,Cdagbeta)
!  IF ( Cdag .GE. beta ) &
!    CALL ERROR("BathOperatoroffdiag_getDetAdd : bad case ...              ")
  IF ( op%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperatoroffdiag_getDetAdd : no active hybrid function ")

  IF ( size(particle)/=op%flavors ) &
    CALL ERROR("BathOperatoroffdiag_getDetAdd : size of particle is erroneous ")
 
 ! tail is now the complete size of the F matrix Fshift(nflavors+1)
  tail =  op%sumtails
  new_tail = tail+1
!  list => particle%list

  if(op%opt_nondiag==1) then
    iflavorbegin = 1
    iflavorend   = op%flavors
    tailbegin    = 1
    tailend      = tail
  else
  !sui!write(6,*) "Bathoperator opt_nondiag=0"
    iflavorbegin = op%activeflavor
    iflavorend   = op%activeflavor
    tailbegin    = op%Fshift(op%activeflavor)+1
    tailend      = op%Fshift(op%activeflavor)+op%tails(op%activeflavor)
  endif
  
  IF ( ((C .GT. Cdag) .AND. (position .EQ. -1)) &  ! Segment added at the end of the segment
       .OR. ((C .LT. Cdag) .AND. (tail .EQ. 0))) THEN ! empty orbital case: only adding a segment is possible
   ! If ones add a segment to an empty orbital or a segment at the end
   ! of a segment, then:
    op%updatePosRow = op%tails(op%activeFlavor) + 1
    op%updatePosCol = op%tails(op%activeFlavor) + 1
  ELSE
   ! For all the other cases, ABS(position) is the true position.
    op%updatePosRow  = ABS(position)
    op%updatePosCol  = ABS(position)
  END IF
    !write(6,*) "       BathOperatoroffdiag_getDetAdd : op%updatePosRow",op%updatePosRow
    !write(6,*) "       BathOperatoroffdiag_getDetAdd : op%updatePosCol",op%updatePosCol
    !write(6,*) "       BathOperatoroffdiag_getDetAdd : C,Cdag",C,Cdag
  
  IF ( C .LT. Cdag .AND. op%tails(op%activeFlavor) .GT. 0) THEN ! only if an antisegment is added
  !  ratio = -ratio 
    op%updatePosRow  = (op%updatePosRow + 1) !position in [1;tail]
  ! If the antisegment created is such that a segment with tcdagger> tc
  ! is suppressed
    !write(6,*) "       BathOperatoroffdiag_getDetAdd : op%updatePosRow",op%updatePosRow
    !write(6,*) "       BathOperatoroffdiag_getDetAdd : op%updatePosCol",op%updatePosCol
    IF ( Cdagbeta .LT. particle(op%activeFlavor)%list(op%updatePosCol,Cdag_) ) op%antiShift = .TRUE.
  END IF

!  CALL Vector_setSize(op%R,tail)
!  CALL Vector_setSize(op%Q,tail)
  Vector_QuickResize(op%R,new_tail)
  Vector_QuickResize(op%Q,new_tail)
  Vector_QuickResize(op%Rtau,new_tail)
  Vector_QuickResize(op%Qtau,new_tail)

!  This loop compute all Row and Col except op%updatePosRow
  tcheck=0
  DO iflavor = iflavorbegin,iflavorend
    !write(6,*) "       BathOperatoroffdiag_getDetAdd : tails(iflavor)",iflavor,op%tails(iflavor)
  DO it1 = 1, op%tails(iflavor)
    tcheck=tcheck+1
    it2 = it1
    it3 = it1
    IF ( iflavor .GE. op%activeFlavor ) THEN
      it2 = it1 + 1
      it3 = it1 + 1
      IF ( iflavor .EQ. op%activeFlavor .AND. it1 .LT. op%updatePosRow ) it2 = it1
      IF ( iflavor .EQ. op%activeFlavor .AND. it1 .LT. op%updatePosCol ) it3 = it1
    !it2 = it1 + ( 1+SIGN(1,it1-op%updatePosRow) )/2
    !it3 = it1 + ( 1+SIGN(1,it1-op%updatePoscol) )/2
    ! if it1>=op%updatePosRow and iflavor> activeflavor, then it2=it1+1
    ! if it1< op%updatePosRow and iflavor> activeflavor, then it2=it1
    END IF

    !!write(6,*) size(op%Rtau%vec)
    !!write(6,*) size(particle(iflavor)%list,1)
    !!write(6,*) size(particle(iflavor)%list,2)
    !!write(6,*) size(op%Fshift)
    !!write(6,*) it1,Cdag_,op%Fshift(iflavor)+it2
    op%Rtau%vec(op%Fshift(iflavor)+it2)= C - particle(iflavor)%list(it1,Cdag_)
    !   the following line happend only for nondiag case
    IF(op%Rtau%vec(op%Fshift(iflavor)+it2) .GT. beta) op%Rtau%vec(op%Fshift(iflavor)+it2)=op%Rtau%vec(op%Fshift(iflavor)+it2)-beta
    !op%Rtau%vec(it1)= C - particle%list(it1,Cdag_)
    time = Cbeta - particle(iflavor)%list(it1,Cdag_)
    if(op%Rtau%vec(op%Fshift(iflavor)+it2)>beta) then
    !write(6,*) "Rtau sup beta",op%Rtau%vec(op%Fshift(iflavor)+it2),C,particle(iflavor)%list(it1,Cdag_)
    !write(6,*) time
    stop
    endif

! "BathOperatoroffdiag_hybrid" interpolates between known values of F for the
!  selected time.
    iflavora=iflavor
    iflavorb=op%activeFlavor
#include "BathOperatoroffdiag_hybrid"

    op%R%vec(op%Fshift(iflavor)+it1) = hybrid
!    op%R%vec(it) = BathOperatoroffdiag_hybrid(op, Cbeta - list(it)%Cdag)
!    Cibeta = list(it)%C.MOD.beta
    MODCYCLE(particle(iflavor)%list(it1,C_),beta,Cibeta)
    time = Cibeta - Cdagbeta
    op%Qtau%vec(op%Fshift(iflavor)+it3)= time
    !op%Qtau%vec(it1)= time

    iflavora=op%activeFlavor
    iflavorb=iflavor
#include "BathOperatoroffdiag_hybrid"
    op%Q%vec(op%Fshift(iflavor)+it1) = hybrid

    !op%Q%vec(it3) = hybrid
!    Q(it) = BathOperatoroffdiag_hybrid(op, Cibeta - Cdagbeta)
  END DO
  END DO
  if(tcheck.ne.tail) then
    !write(6,*) " PRB in the loop tail tcheck",tail,tcheck
    stop
  endif

  ! Compute S
  op%Stau = C - Cdagbeta 
  op%Rtau%vec(op%Fshift(op%activeFlavor)+op%updatePosRow) = op%Stau
    if(op%Rtau%vec(op%Fshift(op%activeFlavor)+op%updatePosRow)>beta) then
    !write(6,*) "Rtau sup beta", op%Stau,C,Cdagbeta
    stop
    endif
  op%Qtau%vec(op%Fshift(op%activeFlavor)+op%updatePosCol) = op%Rtau%vec(op%Fshift(op%activeFlavor)+op%updatePosRow)
  !write(6,*) "              getdetAdd op%Stau",op%Stau

  time = Cbeta-Cdagbeta
  !write(6,*) "              getdetAdd time",time
  iflavora=op%activeFlavor
  iflavorb=op%activeFlavor
  !write(6,*) "time",time
#include "BathOperatoroffdiag_hybrid"
  !write(6,*) "hybrid",hybrid
  op%S = hybrid
  !write(6,*) "              getdetAdd hybrid=op%S",hybrid

  !ratio = op%S - DOT_PRODUCT(MATMUL(op%R%vec(1:tail),op%M(op%activeFlavor)%mat(1:tail,1:tail)),op%Q%vec(1:tail))

  ! product of matrix R and M(k) is computed now:
  ratio = 0.d0
  DO it1 = tailbegin, tailend
    time = 0.d0
    DO it2 = tailbegin, tailend
      time = time + op%R%vec(it2) * op%M%mat(it2,it1)
    END DO
    ratio = ratio + op%Q%vec(it1) * time
  END DO
  !sui!write(6,*) "        = R Matrix",tail
  !sui!write(6,*) "        R      ",(op%R%vec(it1),it1=1,tail)
  !sui!write(6,*) "        = Q Matrix",tail
  !sui!do it1=1,tail
  !sui!write(6,*) "        Q      ",op%Q%vec(it1)
  !sui!enddo
  !sui!write(6,*) "        = M Matrix",tail
  !sui!do it2=1,tail
  !sui!write(6,*) "        M      ",(op%M%mat(it2,it1),it1=1,tail)
  !sui!enddo
  !sui!write(6,*) "        RMQ    =", ratio
  !sui!write(6,*) "         S     =", op%S
 ratio = op%S - ratio
 !sui!write(6,*) "         S-RMQ =", ratio
 !sui!write(6,*) "              getdetAdd ratio",ratio

  op%Stilde = 1.d0 / ratio
  ! If antisegment, the det ratio has to be multiplied by -1 ( sign of the signature of one
  ! permutation line in the matrix)
  IF ( C .LT. Cdag .AND. op%tails(op%activeFlavor) .GT. 0) THEN ! only if an antisegment is added
    ratio=-ratio
  ENDIF

  ! This IF is the LAST "NON CORRECTION" in my opinion this should not appears.
!  IF ( MAX(C,Cdag) .GT. op%beta ) THEN
!    WRITE(*,*) op%Stilde
!    op%Stilde = - ABS(op%Stilde)
!  END IF
  BathOperatoroffdiag_getDetAdd = ratio
      !write(6,*) " getdetAdd",ratio,BathOperatoroffdiag_getDetAdd
  op%MAddFlag   = .TRUE.
!#ifdef CTQMC_CHECK
!  op%ListCdagC = particle
!!write(*,*) op%Stilde
!!write(*,*) op%antishift
!!write(*,*)    op%updatePosRow 
!!write(*,*)    op%updatePosCol 
!#endif

END FUNCTION BathOperatoroffdiag_getDetAdd
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_getDetRemove
!! NAME
!!  BathOperatoroffdiag_getDetRemove
!!
!! FUNCTION
!!  Compute the determinant ratio when a (anti)segment
!!  is trying to be removed 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
!!  position=position of segment to be removed
!!
!! OUTPUT
!!  BathOperatoroffdiag_getDetRemove=the det 
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

DOUBLE PRECISION FUNCTION BathOperatoroffdiag_getDetRemove(op,position)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
!Local arguments-------------------------------
  INTEGER           , INTENT(IN   ) :: position  
  INTEGER                           :: ABSposition  
  INTEGER                           :: tail !,it,it1

  IF ( op%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperatoroffdiag_getDetRemove : no active hybrid fun  ")

  op%antiShift = .FALSE.
  tail         = op%sumtails
  ABSposition  = ABS(position)
  IF ( ABSposition .GT. op%tails(op%activeFlavor) ) &
    CALL ERROR("BathOperatoroffdiag_getDetRemove : position > M size     ")
  op%updatePosCol = ABSposition
  op%antiShift    = .FALSE.
  IF ( position .GT. 0 ) THEN
    op%updatePosRow = ABSposition
  ELSE
    op%updatePosRow = ABSposition+1
    IF ( ABSposition .EQ. op%tails(op%activeFlavor) ) THEN 
      op%antiShift = .TRUE.
      op%updatePosRow = 1 !ABSposition - 1
!      op%updatePosRow = ABSposition    
!      IF ( op%updatePosCol .EQ. 0) op%updatePosCol = tail
    END IF
  ENDIF
  op%Stilde                 = op%M%mat(op%Fshift(op%activeFlavor)+&
&                     op%updatePosRow,op%Fshift(op%activeFlavor)+op%updatePosCol) 
!sui!write(6,*) "Fshift",op%Fshift(op%activeFlavor)
!sui!write(6,*) "updatepos",op%updatePosRow,op%updatePosCol
  
 
  op%MRemoveFlag            = .TRUE.
       !write(6,*) "        getdetRemove",op%Stilde
  BathOperatoroffdiag_getDetRemove = op%Stilde
  if(position<0.and.op%tails(op%activeFlavor)>1) then
    BathOperatoroffdiag_getDetRemove = -op%Stilde
  endif
  !do it=1,op%sumtails
  !!sui!write(6,*) "        getdetRemove M",(op%M%mat(it,it1),it1=1,op%sumtails)
  !enddo
!#ifdef CTQMC_CHECK
!  op%ListCdagC = particle
!!write(*,*) op%updatePosRow, op%updatePosCol, position
!!CALL ListCdagC_print(particle)
!#endif

END FUNCTION BathOperatoroffdiag_getDetRemove
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_getDetF
!! NAME
!!  BathOperatoroffdiag_getDetF
!!
!! FUNCTION
!!  Compute the determinant of the F matrix
!!  using the hybridization of flavor and the 
!!  segments of particle
!!  used for Gloval moves only
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
!!  flavor=hybridization function to take
!!  particles=segments to use
!!
!! OUTPUT
!!  BathOperatoroffdiag_getDetF=the det 
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

DOUBLE PRECISION FUNCTION BathOperatoroffdiag_getDetF(op,particle,option)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag)       , INTENT(INOUT)      :: op
  TYPE(ListCdagC), OPTIONAL, INTENT(IN   )  :: particle(:)
  INTEGER , optional :: option
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
  INTEGER :: iflavor,iflavora
  INTEGER :: iflavordag,iflavorb
#include "BathOperatoroffdiag_hybrid.h"

  BathOperatoroffdiag_getDetF = 1.d0 ! pour eviter des divisions par 0
  IF ( PRESENT( particle ) ) THEN
    tail = op%sumtails
    beta = op%beta
    mbeta_two = -beta*0.5d0
    inv_dt =  op%inv_dt
    CALL MatrixHyb_setSize(op%M_update,tail)
    DO iflavordag=1,op%flavors
    DO iCdag = 1, op%tails(iflavordag)
      tCdag  = particle(iflavordag)%list(iCdag,Cdag_)
      DO iflavor=1,op%flavors
      DO iC  = 1, op%tails(iflavor)
        !tC   = particle%list(C_,iC).MOD.beta
        MODCYCLE(particle(iflavor)%list(iC,C_),beta,tC)
        time = tC - tCdag
        iflavora=iflavordag
        iflavorb=iflavor
#include "BathOperatoroffdiag_hybrid"
        op%M_update%mat(op%Fshift(iflavor)+iC,op%Fshift(iflavordag)+iCdag) = hybrid 
      END DO
      END DO
    END DO
    END DO
    ! mat_tau needs to be transpose of ordered time mat (way of measuring
    ! G(tau))
    DO iflavor=1,op%flavors
    DO iC  = 1, tail
      tC   = particle(iflavor)%list(iC,C_)
      DO iflavordag=1,op%flavors
      DO iCdag = 1, tail
    !sui!write(6,*) iCdag,Cdag_,size(particle(iflavordag)%list,1) 
      !stop 
        tCdag  = particle(iflavordag)%list(iCdag,Cdag_)
        time = tC - tCdag
        signe = SIGN(1.d0,time)
        time = time + (signe-1.d0)*mbeta_two
        op%M_update%mat_tau(op%Fshift(iflavordag)+iCdag,op%Fshift(iflavor)+iC) = INT( ( time * inv_dt ) + 1.5d0 )
      END DO
      END DO
    END DO
    END DO
    CALL MatrixHyb_inverse(op%M_update,BathOperatoroffdiag_getDetF) ! calcul le det de la matrice et l'inverse
  ELSE
    if(present(option)) then
      CALL MatrixHyb_getDet(op%M_update,BathOperatoroffdiag_getDetF) ! det M = 1/detF !
    else
      CALL MatrixHyb_getDet(op%M,BathOperatoroffdiag_getDetF) ! det M = 1/detF !
    endif
    BathOperatoroffdiag_getDetF = 1.d0 / BathOperatoroffdiag_getDetF
  ENDIF
END FUNCTION BathOperatoroffdiag_getDetF
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_setMAdd
!! NAME
!!  BathOperatoroffdiag_setMAdd
!!
!! FUNCTION
!!  Update de M matrix inserting a row and a column
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_setMAdd(op,particle) 

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
  TYPE(ListCdagC)   , INTENT(IN   ) :: particle(:)
!Local variables ------------------------------
  INTEGER                           :: tail
  INTEGER                           :: new_tail
  INTEGER                           :: col
  INTEGER                           :: col_move
  INTEGER                           :: row_move
  INTEGER                           :: row
  INTEGER                           :: positionRow
  INTEGER                           :: positionCol
  INTEGER                           :: aF,indice
  INTEGER                           :: tailb,taile
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
  INTEGER :: p,it !,it1



! ---  op%MAddFlag is put to .TRUE. in BathOperatoroffdiag_getDetAdd.
  IF ( op%MAddFlag .EQV. .FALSE. ) &
    CALL ERROR("BathOperatoroffdiag_setMAdd : MAddFlag turn off           ")

! ---  op%activeFlavor is put in ctqmc_loop
  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("BathOperatoroffdiag_setMAdd : no active hybrid function   ")

!!  do it=1,op%sumtails
    !write(6,*) "        setMAdd begin M",(op%M%mat(it,it1),it1=1,op%sumtails)
!!  enddo
!!  do it=1,op%sumtails
    !write(6,*) "        setMAdd begin M%mat_tau",(op%M%mat_tau(it,it1),it1=1,op%sumtails)
!!  enddo
!  old tail
  !write(6,*) "       BathOperatoroffdiag_setMAdd op%sumtails",op%sumtails
  tail = op%sumtails
  new_tail =  tail + 1
  op%tails(aF)= op%tails(aF) + 1
  DO indice = aF +1, op%flavors+1
    op%Fshift(indice) = op%Fshift(indice) + 1
  END DO
  op%sumtails = op%Fshift(op%flavors) + op%tails(op%flavors) !last slot of Fshift is the tail of full matrix
  !write(6,*) "       BathOperatoroffdiag_setMAdd op%sumtails",op%sumtails
  !write(6,*) "        setMAdd actualized Fshift",(op%Fshift(it),it=1,op%flavors+1)
  !write(6,*) "        setMAdd actualized tails",(op%tails(it),it=1,op%flavors)
  !CALL matrix_print(M)

  if(op%opt_nondiag==1) then
    tailb        = 1
    taile        = tail
  else
  !sui!write(6,*) "Bathoperator a opt_nondiag=0"
    tailb        = op%Fshift(aF)+1
    taile        = op%Fshift(aF)+op%tails(aF)
  endif

! ---  data obtained from BathOperatoroffdiag_getDetAdd
  PositionRow =  op%updatePosRow + op%Fshift(aF) ! position in the full matrix
  PositionCol =  op%updatePosCol + op%Fshift(aF) ! position in the full matrix
  Stilde      =  op%Stilde

!  !write(6,*) "before", positionRow, positionCol
  !CALL MatrixHyb_print(op%M(aF),opt_print=1)
! ---  MatrixHyb_setSize
  !write(6,*) "       BathOperatoroffdiag_setMAdd before setsize",size(op%M%mat,1)
  CALL MatrixHyb_setSize(op%M,new_tail)
  !write(6,*) "       BathOperatoroffdiag_setMAdd after setsize",size(op%M%mat,1)

  ! Compute Qtilde with Q
  !op%Q%vec(1:tail) = (-1.d0) * MATMUL(op%M(aF)%mat(1:tail,1:tail),op%Q%vec(1:tail)) * Stilde

! ---  M*Q => Q
  op%Q%vec(tailb:taile) = MATMUL(op%M%mat(tailb:taile,tailb:taile),op%Q%vec(tailb:taile))

  !op%Q%vec(PositionRow:new_tail) = EOSHIFT(op%Q%vec(PositionRow:new_tail), SHIFT=-1, BOUNDARY=-1.d0, DIM=1)
!  op%Qtau%vec(PositionCol:new_tail) = EOSHIFT(op%Qtau%vec(PositionCol:new_tail), SHIFT=-1, BOUNDARY=1.d0, DIM=1)
!  op%Qtau%vec(PositionCol) = op%Stau

  !Compute Rtilde with R and without multiplying by Stilde
  !op%R%vec(1:tail) = (-1.d0) * MATMUL(op%R%vec(1:tail),op%M(aF)%mat(1:tail,1:tail))

! ---  R*M => R
  op%R%vec(tailb:taile) = MATMUL(op%R%vec(tailb:taile),op%M%mat(tailb:taile,tailb:taile))

  !op%R%vec(PositionCol:new_tail) = EOSHIFT(op%R%vec(PositionCol:new_tail), SHIFT=-1, BOUNDARY=-1.d0, DIM=1)
!  op%Rtau%vec(PositionRow:new_tail) = EOSHIFT(op%Rtau%vec(PositionRow:new_tail), SHIFT=-1, BOUNDARY=1.d0, DIM=1)
!  op%Rtau%vec(PositionRow) = op%Stau

  !Compute the new M matrix
  !op%M(aF)%mat(PositionRow:new_tail,1:new_tail) = &
  !                   EOSHIFT(op%M(aF)%mat(PositionRow:new_tail,1:new_tail),SHIFT=-1, BOUNDARY=0.d0, DIM=1)
  !op%M(aF)%mat(1:n12 characters (ABI_ALLOCATE) instead of 4 (FREE)ew_tail,PositionCol:new_tail) = &
  !                   EOSHIFT(op%M(aF)%mat(1:new_tail,PositionCol:new_tail),SHIFT=-1, BOUNDARY=0.d0, DIM=2)
! ! op%M(aF)%mat(1:new_tail,1:new_tail) =  op%M(aF)%mat(1:new_tail,1:new_tail) + &
! ! Stilde * MATMUL(RESHAPE(op%Q%vec(1:new_tail),(/ new_tail,1 /)),RESHAPE(op%R%vec(1:new_tail),(/ 1,new_tail /)))

  !op%M(aF)%mat_tau(PositionRow:new_tail,1:new_tail) = &
  !                   EOSHIFT(op%M(aF)%mat_tau(PositionRow:new_tail,1:new_tail),SHIFT=-1, BOUNDARY=0, DIM=1)
  !op%M(aF)%mat_tau(1:new_tail,PositionCol:new_tail) = &
  !                   EOSHIFT(op%M(aF)%mat_tau(1:new_tail,PositionCol:new_tail),SHIFT=-1, BOUNDARY=0, DIM=2)

  mbeta_two = -op%beta*0.5d0
  inv_dt = op%inv_dt

! ------ Shift mat_tau and update old M=Ptilde
  DO col=tail,1,-1  ! decreasing order to avoid overwrite of data
    col_move = col +  ( 1+SIGN(1,col-PositionCol) )/2
    ! if col>= PositionCol col_move=col+1
    ! if col<  PositionCol col_move=col
    DO row=tail,1,-1
      row_move = row +  ( 1+SIGN(1,row-PositionRow) )/2
! ---  times for Ptilde are kept unchanged. But we have to copy it at the right place
      op%M%mat_tau(row_move,col_move) =  &
      op%M%mat_tau(row,col)
! ---  Update Ptilde with the same indices as mat_tau
! ---  M + M*Q Stilde R*M => Ptilde => M
      !if(row>=tailb.and.row<=taile.and.col>=tailb.and.col<=taile) then
        op%M%mat(row_move,col_move) =  &
        op%M%mat(row,col) + op%Q%vec(row)*op%R%vec(col) * Stilde
      !else
      !  op%M%mat(row_move,col_move) = op%M%mat(row,col) 
      !endif
    END DO
  END DO

! ------ Add new stuff for new row
  DO row = 1, tail
    row_move = row +  ( 1+SIGN(1,row-PositionRow) )/2
! ---  M*Q Stilde => Qtilde => M with the good indices
    !if(row>=tailb.and.row<=taile) then
      op%M%mat(row_move,PositionCol) = -op%Q%vec(row)*Stilde
    !else
    !  op%M%mat(row_move,PositionCol) = op%M%mat(row,PositionCol) 
    !endif

    time = op%Rtau%vec(row) !  pourquoi Rtau et pas Qtau ici ?
    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
    ! if time>=0 time=time
    ! if time< 0 time=time + beta
! ---  mat_tau=int(time*L/beta+1.5)
    op%M%mat_tau(row,PositionCol) = INT ( (time*inv_dt) +1.5d0 )
    !write(6,*) "     setMadd new row", op%Rtau%vec(row),op%M%mat_tau(row,PositionCol)
!    if(op%M%mat_tau(row,PositionCol)>301) then
!      !write(6,*) ">301 a", time,inv_dt, op%M%mat_tau(row,PositionCol)
!    time = op%Rtau%vec(row) !  pourquoi Rtau et pas Qtau ici ?
!      !write(6,*) time,mbeta_two
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!      !write(6,*) time
!      !write(6,*) INT ( (time*inv_dt) +1.5d0 )
!      stop
!    endif
  END DO
  ! Add last time missing in the loops
  time = op%Rtau%vec(new_tail)
  time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  op%M%mat_tau(new_tail,PositionCol) = INT ( (time*inv_dt) +1.5d0 )
    !write(6,*) "     setMadd last time", op%Rtau%vec(new_tail),op%M%mat_tau(new_tail,PositionCol)
!    if(op%M%mat_tau(new_tail,PositionCol)>301) then
!      !write(6,*) ">301 b", time,inv_dt, op%M%mat_tau(new_tail,PositionCol)
!    time = op%Rtau%vec(new_tail) !  pourquoi Rtau et pas Qtau ici ?
!      !write(6,*) time,mbeta_two
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!      !write(6,*) time
!      !write(6,*) INT ( (time*inv_dt) +1.5d0 )
!      stop
!    endif

  ! Add new stuff for new col
  DO col = 1, tail 
    col_move = col +  ( 1+SIGN(1,col-PositionCol) )/2
! ---   Stilde RN => Rtilde => M
    !if(col>=tailb.and.col<=taile) then
      op%M%mat(PositionRow,col_move) = -op%R%vec(col)*Stilde
    !else
    !  op%M%mat(PositionRow,col_move) = op%M%mat(PositionRow,col) 
    !endif
    time = op%Qtau%vec(col)
    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
    op%M%mat_tau(PositionRow,col) = INT ( (time*inv_dt) +1.5d0 )
    !write(6,*) "     setMadd new col", op%Qtau%vec(col),op%M%mat_tau(PositionRow,col)
!    if(op%M%mat_tau(PositionRow,col)>301) then
!      !write(6,*) ">301 c", time,inv_dt, op%M%mat_tau(PositionRow,col)
!    time = op%Qtau%vec(col) !  pourquoi Rtau et pas Qtau ici ?
!      !write(6,*) time,mbeta_two
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!      !write(6,*) time
!      !write(6,*) INT ( (time*inv_dt) +1.5d0 )
!      stop
!    endif
  END DO
  ! Add last time missing in the loops
  time = op%Qtau%vec(new_tail)
  time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  op%M%mat_tau(PositionRow,new_tail) = INT ( (time*inv_dt) +1.5d0 )
    !write(6,*) "     setMadd last time", op%Qtau%vec(new_tail),op%M%mat_tau(PositionRow,new_tail) 
!    if(op%M%mat_tau(PositionRow,new_tail)>301) then
!      !write(6,*) ">301 d", time,inv_dt, op%M%mat_tau(PositionRow,new_tail)
!    time = op%Qtau%vec(new_tail) !  pourquoi Rtau et pas Qtau ici ?
!      !write(6,*) time,mbeta_two
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!      !write(6,*) time
!      !write(6,*) INT ( (time*inv_dt) +1.5d0 )
!      stop
!    endif

  op%M%mat(PositionRow,PositionCol) = Stilde

  !CALL MatrixHyb_print(op%M,opt_print=1)

!  DO col = 1, new_tail
!    time = op%Rtau%vec(col)
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!    op%M(aF)%mat_tau(col,PositionCol) = INT ( (time*inv_dt) +1.5d0 )
!    time = op%Qtau%vec(col)
!    time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
!    op%M(aF)%mat_tau(PositionRow,Col) = INT ( (time*inv_dt) +1.5d0 )
!    time = op%R%vec(col)*Stilde
!    DO row = 1, new_tail
!      op%M(aF)%mat(row,col) = op%M(aF)%mat(row,col) + op%Q%vec(row)*time
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
  !        op%M(aF)%mat(row_move,col_move) = Stilde
  !        !time = op%Stau
  !      ELSE
  !        op%M(aF)%mat(row_move,col_move) = -op%Q%vec(row)*Stilde
  !        !time = op%Rtau%vec(row_move)
  !        row      = row      - 1 
  !      END IF
  !      !time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  !      !op%M(aF)%mat_tau(row_move,col_move) = INT ( (time*inv_dt) +1.5d0 )
  !    END DO
  !    ! realignement des indices
  !  ELSE
  !    ! on calcule Ptilde
  !    !row_move = new_tail
  !    row      = tail 
  !    DO row_move = new_tail, 1, -1
  !      IF ( row_move .EQ. positionRow ) THEN
  !        op%M(aF)%mat(row_move,col_move) = -op%R%vec(col) * Stilde
  !        ! calcul itau
  !        !time = op%Qtau%vec(col_move)
  !        !time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
  !        !op%M(aF)%mat_tau(row_move,col_move) = INT ( (time*inv_dt) +1.5d0 )
  !      ELSE
  !        op%M(aF)%mat(row_move,col_move) = op%M(aF)%mat(row,col) + op%Q%vec(row)*op%R%vec(col)*Stilde
  !        ! copy itau
  !        !op%M(aF)%mat_tau(row_move,col_move) = op%M(aF)%mat_tau(row,col)
  !        row      = row      - 1 
  !      END IF
  !    END DO
  !    col      = col      - 1
  !  END IF
  !END DO
!  !write(6,*) "after"
!  CALL MatrixHyb_print(op%M(aF),opt_print=1)
!CALL matrix_inverse(M)
!CALL MatrixHyb_print(M)
!CALL matrix_inverse(M)

  IF ( op%antiShift .EQV. .TRUE. ) THEN ! antisegment
  if(3==4) then
    CALL Vector_init(vec_tmp,new_tail)
    CALL VectorInt_init(vecI_tmp,new_tail)
  ! Shift if necessary according to op%antishift
  ! shift DIM=2 (col)

! For new_tail=4, the following lines transform
! M=(a,b,c,d) vith a,b,c,d column vectors into
! M=(d,a,b,c)
    p = new_tail - 1  ! = tail
    m = 1
!   count increases in the loop from 0 to new_tail-1
    count = 0
    DO WHILE ( count .NE. new_tail )
      ! put column b in vec_tmp
      vec_tmp%vec(1:new_tail) = op%M%mat(1:new_tail,m)
      vecI_tmp%vec(1:new_tail) = op%M%mat_tau(1:new_tail,m)
      i = m
      !j = m+p
      MODCYCLE(m+p, new_tail, j)   ! j=m+p modulo new_tail
      DO WHILE (j .NE. m)
        op%M%mat(1:new_tail,i) = op%M%mat(1:new_tail,j)
        op%M%mat_tau(1:new_tail,i) = op%M%mat_tau(1:new_tail,j)
        i = j
        MODCYCLE(j+p, new_tail, j)
        count = count+1
      END DO
      op%M%mat(1:new_tail,i) = vec_tmp%vec(1:new_tail)
      op%M%mat_tau(1:new_tail,i) = vecI_tmp%vec(1:new_tail)
      count = count+1
      m = m+1
    END DO
    ! shift DIM=1 (row)

!   below is similar to above but for rows instead of columns.
    p = new_tail - 1
    m = 1
    count = 0
    DO WHILE ( count .NE. new_tail)
      vec_tmp%vec(1:new_tail) = op%M%mat(m,1:new_tail)
      vecI_tmp%vec(1:new_tail) = op%M%mat_tau(m,1:new_tail)
      i = m
      !j = m+p
      MODCYCLE(m+p, new_tail, j)
      DO WHILE ( j .NE. m )
        op%M%mat(i,1:new_tail) = op%M%mat(j,1:new_tail)
        op%M%mat_tau(i,1:new_tail) = op%M%mat_tau(j,1:new_tail)
        i = j
        MODCYCLE(j+p, new_tail, j)
        count = count+1
      END DO
      op%M%mat(i,1:new_tail) = vec_tmp%vec(1:new_tail)
      op%M%mat_tau(i,1:new_tail) = vecI_tmp%vec(1:new_tail)
      count = count+1
      m = m+1
    END DO
    CALL Vector_destroy(vec_tmp)
    CALL VectorInt_destroy(vecI_tmp)
  endif
    !op%M(aF)%mat(1:new_tail,1:new_tail) = CSHIFT(op%M(aF)%mat(1:new_tail,1:new_tail), SHIFT=-1, DIM=1) ! Shift to the bottom
    !op%M(aF)%mat(1:new_tail,1:new_tail) = CSHIFT(op%M(aF)%mat(1:new_tail,1:new_tail), SHIFT=-1, DIM=2) ! Shift to the right
    !op%M(aF)%mat_tau(1:new_tail,1:new_tail) = CSHIFT(op%M(aF)%mat_tau(1:new_tail,1:new_tail), SHIFT=-1, DIM=1) ! Shift to the bottom
    !op%M(aF)%mat_tau(1:new_tail,1:new_tail) = CSHIFT(op%M(aF)%mat_tau(1:new_tail,1:new_tail), SHIFT=-1, DIM=2) ! Shift to the right
  !write(6,*) "        setMAdd size M%mat",size(op%M%mat,1),size(op%M%mat,2),new_tail
  !write(6,*) "        setMAdd arguement M%mat",aF,op%Fshift(aF)
  !write(6,*) "        setMAdd arguement M%mat",aF,op%Fshift(aF),op%Fshift(aF+1)
  do it=1,op%sumtails
    !write(6,*) "        setMAdd before antishift M%mat_tau",(op%M%mat_tau(it,it1),it1=1,op%sumtails)
  enddo
    op%M%mat(op%Fshift(aF)+1:op%Fshift(aF+1) , 1:new_tail) = & 
        CSHIFT( op%M%mat(op%Fshift(aF)+1:op%Fshift(aF+1) , 1:new_tail) , SHIFT=-1 , DIM=1) ! Shift to the bottom

    op%M%mat(1:new_tail , op%Fshift(aF)+1:op%Fshift(aF+1)) = &
        CSHIFT( op%M%mat(1:new_tail , op%Fshift(aF)+1:op%Fshift(aF+1)) , SHIFT=-1 , DIM=2) ! Shift to the right

    op%M%mat_tau(op%Fshift(aF)+1:op%Fshift(aF+1) , 1:new_tail) = &
        CSHIFT( op%M%mat_tau(op%Fshift(aF)+1:op%Fshift(aF+1) , 1:new_tail) , SHIFT=-1 , DIM=1) ! Shift to the bottom

    op%M%mat_tau(1:new_tail , op%Fshift(aF)+1:op%Fshift(aF+1)) = &
        CSHIFT( op%M%mat_tau(1:new_tail , op%Fshift(aF)+1:op%Fshift(aF+1)) , SHIFT=-1 , DIM=2) ! Shift to the right
  !CALL matrix_print(M)
  END IF

!!  do it=1,op%sumtails
    !write(6,*) "        setMAdd end M",(op%M%mat(it,it1),it1=1,op%sumtails)
!!  enddo
 !! do it=1,op%sumtails
!!    !write(6,*) "        setMAdd end M%mat_tau",(op%M%mat_tau(it,it1),it1=1,op%sumtails)
 !! enddo
  IF ( op%doCheck .EQV. .TRUE.) THEN
!#ifdef CTQMC_CHECK
    CALL BathOperatoroffdiag_checkM(op,particle)
!#endif
  END IF

  op%MAddFlag = .FALSE.

END SUBROUTINE BathOperatoroffdiag_setMAdd
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_setMRemove
!! NAME
!!  BathOperatoroffdiag_setMRemove
!!
!! FUNCTION
!!  delete one row and one column of the M matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_setMRemove(op,particle) 

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT)  :: op
  TYPE(ListCdagC)   , INTENT(IN   )  :: particle(:)
!Local variables ------------------------------
  INTEGER                            :: tail,tailb,taile
  INTEGER                            :: new_tail
  INTEGER                            :: col
  INTEGER                            :: col_move
  INTEGER                            :: row_move
  INTEGER                            :: row
  INTEGER                            :: positionCol
  INTEGER                            :: positionRow
  INTEGER                            :: aF,iaf
  INTEGER                              :: m
  INTEGER                              :: count
  INTEGER                              :: i
  INTEGER                              :: j,it !,it1
  INTEGER                              :: p
  DOUBLE PRECISION                   :: invStilde
  DOUBLE PRECISION                   :: invStilde2
  TYPE(VectorInt) :: vecI_tmp
  TYPE(Vector)    :: vec_tmp

  IF ( op%MRemoveFlag .EQV. .FALSE. ) &
    CALL ERROR("BathOperatoroffdiag_setMRemove : MRemoveFlag turn off     ")
  aF = op%activeFlavor
  IF ( aF .LE. 0 ) &
    CALL ERROR("BathOperatoroffdiag_setMRemove : no active hybrid func    ")
  do it=1,op%sumtails
    !write(6,*) "        setMRemove begin M",(op%M%mat(it,it1),it1=1,op%sumtails)
  enddo
  tail        =  op%sumtails
  new_tail    =  tail - 1
  op%tails(af)= op%tails(af) - 1
  DO iaf=af+1 , op%flavors+1
    op%Fshift(iaf) = op%Fshift(iaf) - 1
  END DO
  op%sumtails = op%Fshift(op%flavors) + op%tails(op%flavors)
  positionCol =  op%updatePosCol + op%Fshift(af)
  positionRow =  op%updatePosRow + op%Fshift(af)
  invStilde   = 1.d0 / op%Stilde
  if(op%opt_nondiag==1) then
    tailb        = 1
    taile        = new_tail
  else
  !sui!write(6,*) "Bathoperator c opt_nondiag=0"
    tailb        = op%Fshift(aF)+1
    taile        = op%Fshift(aF)+op%tails(aF)
  endif

!  !write(6,*) "before", positionRow, positionCol
!  CALL MatrixHyb_print(op%M(aF),opt_print=1)

!  IF ( new_tail .EQ. 0 ) THEN
!!    IF ( op%antiShift .EQV. .TRUE.  ) THEN
!!      op%M(aF)%mat(1,1) = 1.d0/BathOperatoroffdiag_Hybrid(op, op%beta)
!!      op%MRemoveFlag = .FALSE.
!!      RETURN
!!    END IF
!    CALL MatrixHyb_clear(op%M(aF))
!    op%MRemoveFlag = .FALSE.
!    RETURN
!  END IF

!  CALL Vector_setSize(op%Q,new_tail)
!  CALL Vector_setSize(op%R,new_tail)
  Vector_QuickResize(op%Q,new_tail)
  Vector_QuickResize(op%R,new_tail)

!  We use R and Q as op%R%vec and op%Q%vec
!  op%R%vec => op%R
!  op%Q%vec => op%Q

  row      = 1
  !row_move = 1
  col      = 1
  !col_move = 1
  DO row_move = 1, new_tail
    IF ( row .EQ. positionRow ) row = row + 1
    IF ( col .EQ. positionCol ) col = col + 1
    !col = row_move + (1+SIGN(1,row_move-positionCol))/2
    !row = row_move + (1+SIGN(1,row_move-positionRow))/2
    op%R%vec(row_move) = op%M%mat(positionRow,col)
    op%Q%vec(row_move) = op%M%mat(row,positionCol)
    row      = row + 1 
    col      = col + 1
  END DO
!!    op%R%vec(1:positionCol-1) = op%M(aF)%mat(positionRow,1:positionCol-1)
!!    op%R%vec(positionCol:new_tail) = op%M(aF)%mat(positionRow,positionCol+1:tail)
!!    op%Q%vec(1:positionRow-1) = op%M(aF)%mat(1:positionRow-1,positionCol)
!!    op%Q%vec(positionRow:new_tail) = op%M(aF)%mat(positionRow+1:tail,positionCol)
!write(*,*) positionRow, positionCol
!CALL MatrixHyb_print(M)
!CALL Vector_print(op%R)
!CALL Vector_print(op%Q)
!CALL ListCdagC_print(op%ListCdagC)

  col      = 1
  DO col_move = 1, new_tail 
    IF ( col_move .EQ. positionCol ) col = col + 1
    !col = col_move + (1+SIGN(1,col_move-positionCol))/2
    row      = 1
    invStilde2 = invStilde * op%R%vec(col_move)
    DO row_move = 1, new_tail
      IF ( row_move .EQ. positionRow ) row = row + 1
      !row = row_move + (1+SIGN(1,row_move-positionRow))/2
!    Compute for all rows and cols M <= M - Q 1/S R
      !if(row_move>=tailb.and.row_move<=taile.and.col_move>=tailb.and.col_move<=taile) then
        op%M%mat(row_move,col_move) = op%M%mat(row,col) &
                                        - op%Q%vec(row_move)*invStilde2
      !else
      !  op%M%mat(row_move,col_move) = op%M%mat(row,col) 
      !endif
      op%M%mat_tau(row_move,col_move) = op%M%mat_tau(row,col)
      row      = row      + 1
    END DO
    col      = col      + 1 
  END DO
  CALL MatrixHyb_setSize(op%M,new_tail)

  IF ( op%antiShift .EQV. .TRUE. ) THEN ! antisegment
   if(3==4) then
    ! Shift if necessary according to op%antishift
    ! shift DIM=2 (col)
    CALL Vector_init(vec_tmp,new_tail)
    CALL VectorInt_init(vecI_tmp,new_tail)
    p = 1
    m = 1
    count = 0
    DO WHILE ( count .NE. new_tail )
      vec_tmp%vec(1:new_tail) = op%M%mat(1:new_tail,m)
      vecI_tmp%vec(1:new_tail) = op%M%mat_tau(1:new_tail,m)
      i = m
      !j = m+p
      MODCYCLE(m+p, new_tail, j)
      DO WHILE (j .NE. m)
        op%M%mat(1:new_tail,i) = op%M%mat(1:new_tail,j)
        op%M%mat_tau(1:new_tail,i) = op%M%mat_tau(1:new_tail,j)
        i = j
        MODCYCLE(j+p, new_tail, j)
        count = count+1
      END DO
      op%M%mat(1:new_tail,i) = vec_tmp%vec(1:new_tail)
      op%M%mat_tau(1:new_tail,i) = vecI_tmp%vec(1:new_tail)
      count = count+1
      m = m+1
    END DO
    CALL Vector_destroy(vec_tmp)
    CALL VectorInt_destroy(vecI_tmp)
    !op%M(aF)%mat(1:new_tail,1:new_tail) = &
    !           CSHIFT(op%M(aF)%mat(1:new_tail,1:new_tail), SHIFT=1, DIM=2) ! Shift to the top
    !op%M(aF)%mat_tau(1:new_tail,1:new_tail) = &
    !           CSHIFT(op%M(aF)%mat_tau(1:new_tail,1:new_tail), SHIFT=1, DIM=2) ! Shift to the top
   endif
    op%M%mat(1:new_tail,op%Fshift(af)+1:op%Fshift(af+1)) = &
               CSHIFT(op%M%mat(1:new_tail,op%Fshift(af)+1:op%Fshift(af+1)), SHIFT=1, DIM=2) ! Shift to the top
    op%M%mat_tau(1:new_tail,op%Fshift(af)+1:op%Fshift(af+1)) = &
               CSHIFT(op%M%mat_tau(1:new_tail,op%Fshift(af)+1:op%Fshift(af+1)), SHIFT=1, DIM=2) ! Shift to the top
  END IF
!  !write(6,*) "after "
!  CALL MatrixHyb_print(op%M(aF),opt_print=1)

  IF ( op%doCheck .EQV. .TRUE. ) THEN
!#ifdef CTQMC_CHECK
  CALL BathOperatoroffdiag_checkM(op,particle)
!#endif
  END IF
  do it=1,op%sumtails
    !write(6,*) "        setMRemove end M",(op%M%mat(it,it1),it1=1,op%sumtails)
  enddo

  op%MRemoveFlag = .FALSE.

END SUBROUTINE BathOperatoroffdiag_setMRemove
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_swap
!! NAME
!!  BathOperatoroffdiag_swap
!!
!! FUNCTION
!!  Recompute 2 M matrix swaping the segments (used for Global moves)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_swap(op, flavor1, flavor2)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER           , INTENT(IN   ) :: flavor1
  INTEGER           , INTENT(IN   ) :: flavor2
  INTEGER            :: ii,iflavort,itmptail,flavora,flavorb
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: mat_temp
  INTEGER         , ALLOCATABLE, DIMENSION(:,:) :: mat_tau_temp

  if(flavor1>flavor2) then
    flavora=flavor2
    flavorb=flavor1
  else
    flavora=flavor1
    flavorb=flavor2
  endif
  MALLOC(mat_temp,(1:op%sumtails,1:op%sumtails))
  MALLOC(mat_tau_temp,(1:op%sumtails,1:op%sumtails))
  !mat_temp= op%M%mat
  !mat_tau_temp= op%M%mat_tau
  !it1=0
  !do iflav1=1,op%flavors
  !  do ii1=1,op%tails(iflav1)
  !    it1=it1+1
  !    it2=0
  !    do iflav2=1,op%flavors
  !      do ii2=1,op%tails(iflav1)
  !        it2=it2+1
  !        op%M%mat(it1,it2)=
  !      enddo
  !    enddo
  !  enddo
  !enddo
  if(3==3) then
    op%M=op%M_update
!     shift block flavorb at the place of flavora (column)
    do ii=1, op%tails(flavorb)
      op%M%mat(op%Fshift(flavora)+1:op%Fshift(flavorb+1) , 1:op%sumtails) = & 
          CSHIFT( op%M%mat(op%Fshift(flavora)+1:op%Fshift(flavorb+1) , 1:op%sumtails) , SHIFT=-1 , DIM=1) 
      op%M%mat_tau(op%Fshift(flavora)+1:op%Fshift(flavorb+1) , 1:op%sumtails) = &
          CSHIFT( op%M%mat_tau(op%Fshift(flavora)+1:op%Fshift(flavorb+1) , 1:op%sumtails) , SHIFT=-1 , DIM=1) 
    enddo

!     shift block flavora at the place of flavorb (column)
    do ii=1, op%tails(flavora)
      op%M%mat(op%Fshift(flavora)+op%tails(flavorb)+&
&      1:op%Fshift(flavorb)+op%tails(flavorb) , 1:op%sumtails) = & 
          CSHIFT( op%M%mat( op%Fshift(flavora)+op%tails(flavorb)&
&          +1:op%Fshift(flavorb)+op%tails(flavorb) , 1:op%sumtails) , SHIFT=1 , DIM=1) 
      op%M%mat_tau(op%Fshift(flavora)+op%tails(flavorb)+1:op%Fshift(flavorb)+&
&      op%tails(flavorb) , 1:op%sumtails) = & 
          CSHIFT( op%M%mat_tau( op%Fshift(flavora)+op%tails(flavorb)+&
&          1:op%Fshift(flavorb)+op%tails(flavorb) , 1:op%sumtails) , SHIFT=1 , DIM=1) 
    enddo

!     shift block flavorb at the place of flavora (row)
    do ii=1, op%tails(flavorb)
      op%M%mat(1:op%sumtails , op%Fshift(flavora)+1:op%Fshift(flavorb+1)) = &
          CSHIFT( op%M%mat(1:op%sumtails , op%Fshift(flavora)+1:op%Fshift(flavorb+1)) , SHIFT=-1 , DIM=2) 
      op%M%mat_tau(1:op%sumtails , op%Fshift(flavora)+1:op%Fshift(flavorb+1)) = &
          CSHIFT( op%M%mat_tau(1:op%sumtails , op%Fshift(flavora)+1:op%Fshift(flavorb+1)) , SHIFT=-1 , DIM=2) 
    enddo

!     shift block flavora at the place of flavorb (row)
    do ii=1, op%tails(flavora)
      op%M%mat(1:op%sumtails , op%Fshift(flavora)+op%tails(flavorb)+1:op%Fshift(flavorb)+op%tails(flavorb)) = &
          CSHIFT( op%M%mat(1:op%sumtails ,op%Fshift(flavora)+op%tails(flavorb)&
&          +1:op%Fshift(flavorb)+op%tails(flavorb)) , SHIFT=1 , DIM=2) 
      op%M%mat_tau(1:op%sumtails ,op%Fshift(flavora)+op%tails(flavorb)+&
&      1:op%Fshift(flavorb)+op%tails(flavorb) ) = &
          CSHIFT( op%M%mat_tau(1:op%sumtails ,op%Fshift(flavora)+&
&          op%tails(flavorb)+1:op%Fshift(flavorb)+op%tails(flavorb) ) , SHIFT=1 , DIM=2) 
    enddo
  endif
  if(3==4) then
    op%M=op%M_update
  endif


  do iflavort=flavora+1,flavorb
    op%Fshift(iflavort)=op%Fshift(iflavort)+op%tails(flavorb)-op%tails(flavora)
  enddo

  itmptail=op%tails(flavora)
  op%tails(flavora)=op%tails(flavorb)
  op%tails(flavorb)=itmptail
  FREE(mat_temp)
  FREE(mat_tau_temp)

END SUBROUTINE BathOperatoroffdiag_swap
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_initF
!! NAME
!!  BathOperatoroffdiag_initF
!!
!! FUNCTION
!!  Copy input hybridization functions from a file
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_initF(op,ifstream)

!Arguments ----------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER           , INTENT(IN   ) :: ifstream
!Local variables ----------------
  INTEGER                           :: iflavor1
  INTEGER                           :: iflavor2                  
  INTEGER                           :: sample

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("BathOperatoroffdiag_initF : BathOperatoroffdiag not set         ")

  DO iflavor1=1,op%flavors
    DO iflavor2=2,op%flavors
      DO sample = 1, op%sizeHybrid
        READ(ifstream,*) op%F(sample,iflavor1,iflavor2)
      END DO
    END DO
  END DO
END SUBROUTINE BathOperatoroffdiag_initF
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_setF
!! NAME
!!  BathOperatoroffdiag_setF
!!
!! FUNCTION
!!  Copy F from input array
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_setF(op,F)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag)               , INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(:,:,:) , INTENT(IN   ) :: F
!Arguments ------------------------------------
  INTEGER                                          :: iflavor1
  INTEGER                                          :: iflavor2
  INTEGER                                          :: sample
  INTEGER                                          :: length

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("BathOperatoroffdiag_setF : BathOperatoroffdiag not set          ")

 length  = SIZE(F)
  IF ( length .NE. (op%flavors * op%flavors * op%sizeHybrid) ) &
    CALL ERROR("BathOperatoroffdiag_setF : wrong input F                 ")

  DO iflavor1=1,op%flavors
    DO iflavor2=1,op%flavors
      DO sample = 1, op%sizeHybrid
      op%F(sample,iflavor1,iflavor2) = F(sample,iflavor1,iflavor2)
      END DO
    END DO
  END DO
END SUBROUTINE BathOperatoroffdiag_setF
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_printF
!! NAME
!!  BathOperatoroffdiag_printF
!!
!! FUNCTION
!!  print F function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_printF(op,ostream)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op
  INTEGER,OPTIONAL  , INTENT(IN   ) :: ostream
!Local variables ------------------------------
  CHARACTER(LEN=4)                  :: aflavor
  CHARACTER(LEN=50)                  :: string
  INTEGER                           :: iflavor1
  INTEGER                           :: iflavor2
  INTEGER                           :: sample
  INTEGER                           :: ostream_val

  IF ( PRESENT(ostream) ) THEN 
    ostream_val = ostream
  ELSE  
    ostream_val = 65
    OPEN(UNIT=ostream_val, FILE="F.dat")
  END IF

  WRITE(aflavor,'(I4)') (op%flavors*op%flavors+1)
  string = '(1x,'//TRIM(ADJUSTL(aflavor))//'E22.14)'
  DO sample = 1, op%sizeHybrid
    WRITE(ostream_val,string) (sample-1)*op%dt, ((op%F(sample,iflavor1,iflavor2),&
                                                 iflavor1=1,op%flavors),iflavor2=1,op%flavors)
  END DO
  !CALL FLUSH(ostream_val)

  IF ( .NOT. PRESENT(ostream) ) &
    CLOSE(ostream_val)

END SUBROUTINE BathOperatoroffdiag_printF
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_printM
!! NAME
!!  BathOperatoroffdiag_printM
!!
!! FUNCTION
!!  print M =F^{-1} matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_printM(op,ostream)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(IN) :: op
  INTEGER, OPTIONAL , INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                        :: ostream_val

  IF ( op%activeFlavor .LE. 0 ) &
    CALL ERROR("BathOperatoroffdiag_printM : no active hybrid function    ")
  ostream_val = 6
  IF ( PRESENT(ostream) ) ostream_val = ostream
  CALL MatrixHyb_print(op%M,ostream_val)
END SUBROUTINE BathOperatoroffdiag_printM
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_printM_matrix
!! NAME
!!  BathOperatoroffdiag_printM_matrix
!!
!! FUNCTION
!!  print M =F^{-1} matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_printM_matrix(op,ostream)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag), INTENT(IN) :: op
  INTEGER, OPTIONAL , INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                        :: iflavor1
  INTEGER                        :: i1,it1,it2
  CHARACTER(LEN=22)              :: string
  CHARACTER(LEN=22)              :: string2
  CHARACTER(LEN=4 )              :: size

  ABI_UNUSED(ostream)

  WRITE(size,'(I4)') op%sumtails
  string ='(i2,x,i3,a,'//TRIM(ADJUSTL(size))//'(E5.2,1x))'
  string2 ='(6x,'//TRIM(ADJUSTL(size))//'(i6))'
  open(unit=222, file="M_matrix.dat")
  open(unit=223, file="M_matrix_tau.dat")
  it1=0
  write(222,string2) ((i1,i1=1,op%tails(iflavor1)),iflavor1=1,op%flavors)
  do iflavor1=1, op%flavors
    do i1=1, op%tails(iflavor1)
      it1=it1+1
      write(222,string) iflavor1,i1,'|',(op%M%mat(it1,it2),it2=1,op%sumtails)
    enddo
  enddo


END SUBROUTINE BathOperatoroffdiag_printM_matrix
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/ BathOperatoroffdiag_destroy
!! NAME
!!   BathOperatoroffdiag_destroy
!!
!! FUNCTION
!!  Deallocate and reset every thing
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE  BathOperatoroffdiag_destroy(op)

  TYPE(BathOperatoroffdiag), INTENT(INOUT) :: op

  CALL MatrixHyb_destroy(op%M)
  CALL MatrixHyb_destroy(op%M_update)

  CALL Vector_destroy(op%R)
  CALL Vector_destroy(op%Q)
  CALL Vector_destroy(op%Rtau)
  CALL Vector_destroy(op%Qtau)
  FREEIF(op%F)
  FREEIF(op%Fshift)
  FREEIF(op%tails)

  op%MAddFlag     = .FALSE.
  op%MRemoveFlag  = .FALSE.
  op%flavors      = 0 
  op%beta         = 0.d0
  op%dt      = 0.d0
  op%inv_dt  = 0.d0
  op%samples      = 0
  op%sizeHybrid   = 0
  op%activeFlavor = 0 
  op%updatePosRow = 0
  op%updatePosCol = 0

END SUBROUTINE BathOperatoroffdiag_destroy
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_doCheck
!! NAME
!!  BathOperatoroffdiag_doCheck
!!
!! FUNCTION
!!  Just store if we perfom check for updates of M
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_doCheck(op,opt_check)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag) , INTENT(INOUT) :: op
  INTEGER            , INTENT(IN   ) :: opt_check
  
  IF ( opt_check .GE. 2 ) &
    op%doCheck = .TRUE.
END SUBROUTINE BathOperatoroffdiag_doCheck
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_checkM
!! NAME
!!  BathOperatoroffdiag_checkM
!!
!! FUNCTION
!!  compute from scratch the M matrix and compar it
!!  with the already computed M matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_checkM(op,particle)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag) , INTENT(INOUT) :: op
  TYPE(ListCdagC)    , INTENT(IN   ) :: particle(:)
!Local variables ------------------------------
!  TYPE(MatrixHyb)                    :: checkMatrix
  LOGICAL :: checkTau
  INTEGER :: tail
  INTEGER :: iC
  INTEGER :: iCdag
  INTEGER :: aF
  INTEGER :: iflavora
  INTEGER :: iflavorb,it !,it1
  CHARACTER(LEN=6) :: a
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: mbeta_two
  DOUBLE PRECISION :: errorabs
  DOUBLE PRECISION :: errormax
  DOUBLE PRECISION :: error1
  DOUBLE PRECISION :: errorrel
  DOUBLE PRECISION :: tc
  DOUBLE PRECISION :: tCdag
  DOUBLE PRECISION :: sumMmat
  DOUBLE PRECISION :: sumCheck
#include "BathOperatoroffdiag_hybrid.h"

  aF = op%activeFlavor
  !Construction de la matrix
  tail = op%sumtails
!  CALL MatrixHyb_init(checkMatrix,op%iTech,size=tail,Wmax=op%samples)
!  CALL MatrixHyb_setSize(checkMatrix,tail)

  ! --- set size of the matrix
  CALL MatrixHyb_setSize(op%M_update,tail)

  ! --- compute useful quantities
  beta   =  op%beta
  mbeta_two = -beta*0.5d0
  op%checkNumber = op%checkNumber + 1
  IF ( tail .NE. op%M%tail ) THEN
    CALL WARN("BathOperatoroffdiag_checkM : tails are different          ")
    RETURN
  END IF

  do it=1,op%sumtails
    !write(6,*) "        checkM begin M_update%mat_tau",(op%M_update%mat_tau(it,it1),it1=1,op%sumtails)
  enddo
  ! --- build matrix
!CALL ListCdagC_print(particle)
  DO iflavora = 1, op%flavors
  DO iCdag = 1, op%tails(iflavora)
    tCdag  = particle(iflavora)%list(iCdag,Cdag_)
      !write(6,*) "         checkM a",iflavora,tCdag
    DO iflavorb = 1, op%flavors
    DO iC  = 1, op%tails(iflavorb)
      !tC   = particle%list(C_,iC).MOD.beta
      MODCYCLE(particle(iflavorb)%list(iC,C_),beta,tC) ! tC is tC, or Tc-Beta if tc>beta
      !write(6,*) "         checkM b",iflavorb,tC
      time = tC - tCdag  ! time is positive or negative but lower than beta
      !write(6,*) "         checkM time",time

#include "BathOperatoroffdiag_hybrid"

      op%M_update%mat(op%Fshift(iflavorb)+iC,op%Fshift(iflavora)+iCdag) = hybrid

      time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
      op%M_update%mat_tau(op%Fshift(iflavora)+iCdag,op%Fshift(iflavorb)+iC) = INT ( (time*op%inv_dt) +1.5d0 ) 
      !write(6,*) "         checkM mat_tau",INT ( (time*op%inv_dt) +1.5d0 )
      !write(6,*) "         checkM shifts",op%Fshift(iflavorb),iCdag,op%Fshift(iflavora),iC
    END DO ! iC
    END DO ! iflavorb
  END DO ! iCdag
  END DO ! iflavora

!    CALL MatrixHyb_Print(checkMatrix)
  ! --- Inverse matrix
  CALL MatrixHyb_inverse(op%M_update)

!    CALL MatrixHyb_Print(checkMatrix)
  do it=1,op%sumtails
    !write(6,*) "        checkM end M_update%mat_tau",(op%M_update%mat_tau(it,it1),it1=1,op%sumtails)
  enddo
  do it=1,op%sumtails
    !write(6,*) "        checkM end M_update",(op%M%mat(it,it1),it1=1,op%sumtails)
  enddo

  ! --- Compare M_update and M to check if calculation of M is correct
  sumMmat =0.d0
  sumCheck=0.d0
  error1 = 0.d0
  errormax = 0.d0
  checkTau = .FALSE.
  DO iCdag = 1, tail
    Do iC =1, tail
        errorrel= ABS((op%M_update%mat(iC, iCdag) - & 
                  op%M%mat(iC,iCdag))/op%M_update%mat(iC,iCdag))
        errorabs= ABS(op%M_update%mat(iC, iCdag) - & 
                  op%M%mat(iC,iCdag))
        IF ( errorrel .gt. errormax .and. errorabs .gt. 0.001d0 ) errormax = errorrel
                 ! write(6,*) "     checkM ", errorrel,errorabs
        IF ( op%M_update%mat_tau(iC,iCdag) .NE. op%M%mat_tau(iC,iCdag) ) then
                checkTau = .TRUE.
                !write(6,*) "op%M_update%mat_tau(iC,iCdag), op%M%mat_tau(iC,iCdag)",op%M_update%mat_tau(iC,iCdag), op%M%mat_tau(iC,iCdag)
                !call flush(6)
         CALL ERROR("BathOperatoroffdiag_checkM : "//a//"%                        ") 
        ENDIF
  
    END DO
  END DO

  IF ( checkTau .EQV. .TRUE. ) THEN
    CALL WARN("BathOperatoroffdiag_checkM : mat_tau differs should be")
    CALL MatrixHyb_print(op%M_update,opt_print=1)
    CALL WARN("BathOperatoroffdiag_checkM : whereas it is")
    CALL MatrixHyb_print(op%M,opt_print=1)
  END IF
  op%meanError = op%meanError + errormax
  IF ( errormax .GT. 1.d0 ) THEN 
    WRITE(a,'(I4)') INT(error1*100.d0)
    !write(6,'(I4)') INT(error1*100.d0)
!    CALL MatrixHyb_Print(op%M)
    CALL WARN("BathOperatoroffdiag_checkM") 
  END IF
!  CALL MatrixHyb_destroy(checkMatrix)
END SUBROUTINE BathOperatoroffdiag_checkM
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_recomputeM
!! NAME
!!  BathOperatoroffdiag_recomputeM
!!
!! FUNCTION
!!  compute from scratch the M matrix 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (B. Amadon)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
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

SUBROUTINE BathOperatoroffdiag_recomputeM(op,particle,flav_i,flav_j)

!Arguments ------------------------------------
  TYPE(BathOperatoroffdiag) , INTENT(INOUT) :: op
  TYPE(ListCdagC)    , INTENT(IN   ) :: particle(:)
  INTEGER :: flav_i,flav_j
!Local variables ------------------------------
!  TYPE(MatrixHyb)                    :: checkMatrix
  INTEGER :: tail
  INTEGER :: iC
  INTEGER :: iCdag
  INTEGER :: aF
  INTEGER :: iflavora
  INTEGER :: iflavorb,it !,it1
  INTEGER :: iflavora_imp
  INTEGER :: iflavorb_imp
  !CHARACTER(LEN=6) :: a
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: mbeta_two
  DOUBLE PRECISION :: tc
  DOUBLE PRECISION :: tCdag
  !DOUBLE PRECISION :: sumMmat
  !DOUBLE PRECISION :: sumCheck
#include "BathOperatoroffdiag_hybrid.h"

  aF = op%activeFlavor
  !Construction de la matrix
  tail = op%sumtails
!  CALL MatrixHyb_init(checkMatrix,op%iTech,size=tail,Wmax=op%samples)
!  CALL MatrixHyb_setSize(checkMatrix,tail)

  ! --- set size of the matrix
  CALL MatrixHyb_setSize(op%M_update,tail)

  ! --- compute useful quantities
  beta   =  op%beta
  mbeta_two = -beta*0.5d0
  op%checkNumber = op%checkNumber + 1
  IF ( tail .NE. op%M%tail ) THEN
    CALL WARN("BathOperatoroffdiag_checkM : tails are different          ")
    RETURN
  END IF

  do it=1,op%sumtails
    !write(6,*) "        checkM begin M_update%mat_tau",(op%M_update%mat_tau(it,it1),it1=1,op%sumtails)
  enddo
  ! --- build matrix
!CALL ListCdagC_print(particle)
  DO iflavora = 1, op%flavors
    iflavora_imp=iflavora
    if(iflavora==flav_i) iflavora_imp=flav_j
    if(iflavora==flav_j) iflavora_imp=flav_i
    DO iCdag = 1, op%tails(iflavora_imp)
      tCdag  = particle(iflavora_imp)%list(iCdag,Cdag_)
        !write(6,*) "         checkM a",iflavora,tCdag
      DO iflavorb = 1, op%flavors
        iflavorb_imp=iflavorb
        if(iflavorb==flav_j) iflavorb_imp=flav_i
        if(iflavorb==flav_i) iflavorb_imp=flav_j
        DO iC  = 1, op%tails(iflavorb_imp)
          !tC   = particle%list(C_,iC).MOD.beta
          MODCYCLE(particle(iflavorb_imp)%list(iC,C_),beta,tC) ! tC is tC, or Tc-Beta if tc>beta
          !write(6,*) "         checkM b",iflavorb,tC
          time = tC - tCdag  ! time is positive or negative but lower than beta
          !write(6,*) "         checkM time",time

#include "BathOperatoroffdiag_hybrid"

          op%M_update%mat(op%Fshift(iflavorb_imp)+iC,op%Fshift(iflavora_imp)+iCdag) = hybrid

          time = time + ( SIGN(1.d0,time) - 1.d0 )*mbeta_two
          op%M_update%mat_tau(op%Fshift(iflavora_imp)+iCdag,op%Fshift(iflavorb_imp)+iC) = INT ( (time*op%inv_dt) +1.5d0 ) 
          !write(6,*) "         checkM mat_tau",INT ( (time*op%inv_dt) +1.5d0 )
          !write(6,*) "         checkM shifts",op%Fshift(iflavorb),iCdag,op%Fshift(iflavora),iC
        END DO ! iC
      END DO ! iflavorb
    END DO ! iCdag
  END DO ! iflavora

!    CALL MatrixHyb_Print(checkMatrix)
  ! --- Inverse matrix
  CALL MatrixHyb_inverse(op%M_update)

!    CALL MatrixHyb_Print(checkMatrix)
  do it=1,op%sumtails
    !write(6,*) "        checkM end M_update%mat_tau",(op%M_update%mat_tau(it,it1),it1=1,op%sumtails)
  enddo
  do it=1,op%sumtails
    !write(6,*) "        checkM end M_update",(op%M%mat(it,it1),it1=1,op%sumtails)
  enddo

  ! --- Compare M_update and M to check if calculation of M is correct
END SUBROUTINE BathOperatoroffdiag_recomputeM
!!***

!!****f* ABINIT/m_BathOperatoroffdiag/BathOperatoroffdiag_getError
!! NAME
!!  BathOperatoroffdiag_getError
!!
!! FUNCTION
!!  compute a percentage error / checkM
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=bath operator
!!
!! OUTPUT
!!  BathOperatoroffdiag_getError=Error in percent
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

DOUBLE PRECISION FUNCTION BathOperatoroffdiag_getError(op)

  TYPE(BathOperatoroffdiag), INTENT(IN) :: op

  IF ( op%doCheck .EQV. .TRUE. ) THEN
    BathOperatoroffdiag_getError = op%meanError / DBLE(op%checkNumber)
  ELSE
    BathOperatoroffdiag_getError = 0.d0
  END IF
END FUNCTION BathOperatoroffdiag_getError
!!***
!#endif

END MODULE m_BathOperatoroffdiag
!!***

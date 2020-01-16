
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!****m* ABINIT/m_Ctqmc
!! NAME
!!  m_Ctqmc
!! 
!! FUNCTION 
!!  Manage and drive all the CTQMC
!!  Should not be used if you don't know what you do
!!  Please use CtqmcInterface
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

MODULE m_Ctqmc

USE m_Global
USE m_GreenHyb
USE m_BathOperator
USE m_ImpurityOperator
USE m_Stat
USE m_FFTHyb
USE m_OurRng
USE m_Vector
#ifdef HAVE_MPI2
USE mpi
#endif

IMPLICIT NONE

!!***

PRIVATE

INTEGER, PARAMETER :: CTQMC_SLICE1 = 100
! Coupe Sweeps en 100
INTEGER, PARAMETER :: CTQMC_SLICE2 = 100
! Coupe modNoise1 en 100
INTEGER, PARAMETER :: CTQMC_SEGME =  1
INTEGER, PARAMETER :: CTQMC_ANTIS = -2
INTEGER, PARAMETER :: CTQMC_ADDED =  3  
INTEGER, PARAMETER :: CTQMC_REMOV =  4
INTEGER, PARAMETER :: CTQMC_DETSI =  5


!!****t* m_Ctqmc/Ctqmc
!! NAME
!!  Ctqmc
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

TYPE, PUBLIC :: Ctqmc

  LOGICAL _PRIVATE :: init = .FALSE.
! Flag: is MC initialized

  LOGICAL _PRIVATE :: set  = .FALSE.
! Flag: ??

  LOGICAL _PRIVATE :: setU = .FALSE.
! Flag: is U Set ?

  LOGICAL _PRIVATE :: inF  = .FALSE.
! Flag: is hybridization fct in input ?

  LOGICAL _PRIVATE :: done = .FALSE.
! Flag: is MC terminated ?

  LOGICAL _PRIVATE :: para = .FALSE.
! Flag:  do we have parameters in input

  LOGICAL _PRIVATE :: have_MPI = .FALSE.
! Flag: 

  INTEGER _PRIVATE :: opt_movie = 0
!

  INTEGER _PRIVATE :: opt_analysis = 0
! correlations 

  INTEGER _PRIVATE :: opt_check = 0
! various check 0
! various check 1 impurity
! various check 2 bath
! various check 3 both

  INTEGER _PRIVATE :: opt_order = 0
! nb of segments max for analysis

  INTEGER _PRIVATE :: opt_noise = 0
! compute noise

  INTEGER _PRIVATE :: opt_spectra = 0
! markov chain FT (correlation time)

  INTEGER _PRIVATE :: opt_levels = 0
! do we have energy levels

  INTEGER _PRIVATE :: flavors
!

  INTEGER _PRIVATE :: measurements
! nb of measure in the MC

  INTEGER _PRIVATE :: samples
! nb of L points

  INTEGER(8) _PRIVATE :: seed
!

  INTEGER _PRIVATE :: sweeps
!

  INTEGER _PRIVATE :: thermalization
!

  INTEGER _PRIVATE :: ostream
! output file

  INTEGER _PRIVATE :: istream
! input file

  INTEGER _PRIVATE :: modNoise1
! measure the noise each modNoise1

  INTEGER _PRIVATE :: modNoise2
! measure the noise each modNoise2

  INTEGER _PRIVATE :: activeFlavor
! orbital on which one do sth now

  INTEGER, DIMENSION(1:2) _PRIVATE :: modGlobalMove
! 1: gloabl move each modglobalmove(1)
! 2: we have done modglobalmove(2) for two different orbitals.

  INTEGER _PRIVATE :: Wmax
! Max freq for FT

  DOUBLE PRECISION, DIMENSION(1:6) _PRIVATE :: stats
! to now how many negative determinant, antisegments,seeme.e.twfs...j

  DOUBLE PRECISION _PRIVATE :: swap
! nb of successfull GM

  INTEGER _PRIVATE :: MY_COMM
! 

  INTEGER _PRIVATE :: rank
!

  INTEGER _PRIVATE :: size
! size of MY_COMM

  DOUBLE PRECISION _PRIVATE :: runTime ! time for the run routine
!  

  DOUBLE PRECISION _PRIVATE :: beta
!

  DOUBLE PRECISION _PRIVATE :: U

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) _PRIVATE :: mu
! levels

  TYPE(GreenHyb)  , ALLOCATABLE, DIMENSION(:    ) _PRIVATE :: Greens 

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:  ) _PRIVATE :: measN 
! measure of occupations (3or4,flavor) 

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:  ) _PRIVATE :: measDE
!  (flavor,flavor) double occupancies
!  (1,1): total energy of correlation.

  DOUBLE PRECISION _PRIVATE :: a_Noise
! Noise a exp (-bx) for the  noise

  DOUBLE PRECISION _PRIVATE :: b_Noise
! Noise a exp (-bx) for the  noise

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) _PRIVATE :: abNoiseG   !(ab,tau,flavor)
! Noise but for G

  TYPE(Vector)             , DIMENSION(1:2) _PRIVATE :: measNoise 
  TYPE(Vector), ALLOCATABLE, DIMENSION(:,:,:) _PRIVATE :: measNoiseG       !(tau,flavor,mod) 
! accumulate each value relataed to measurenoise 1 2

  DOUBLE PRECISION _PRIVATE                            :: inv_dt
! 1/(beta/L)

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:  ) _PRIVATE :: measPerturbation 
! opt_order,nflavor

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) _PRIVATE :: measCorrelation 
! segment,antisegment,nflavor,nflavor

  DOUBLE PRECISION _PRIVATE :: errorImpurity
! check 

  DOUBLE PRECISION _PRIVATE :: errorBath
! for check

  TYPE(BathOperator) _PRIVATE              :: Bath

  TYPE(ImpurityOperator) _PRIVATE          :: Impurity

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) _PRIVATE :: density

END TYPE Ctqmc
!!***

PUBLIC  :: Ctqmc_init
PRIVATE :: Ctqmc_setParameters
PUBLIC  :: Ctqmc_setSweeps
PUBLIC  :: Ctqmc_setSeed
PRIVATE :: Ctqmc_allocateAll
PRIVATE :: Ctqmc_allocateOpt
PUBLIC  :: Ctqmc_setG0wFile
PUBLIC  :: Ctqmc_setG0wTab
PUBLIC  :: Ctqmc_setU
PRIVATE :: Ctqmc_clear
PUBLIC  :: Ctqmc_reset
PUBLIC  :: Ctqmc_setMu
PRIVATE :: Ctqmc_computeF
PUBLIC  :: Ctqmc_run
PRIVATE :: Ctqmc_loop
PRIVATE :: Ctqmc_tryAddRemove
PRIVATE :: Ctqmc_trySwap
PRIVATE :: Ctqmc_measN
PRIVATE :: Ctqmc_measCorrelation
PRIVATE :: Ctqmc_measPerturbation
PUBLIC  :: Ctqmc_getResult
PUBLIC  :: Ctqmc_symmetrizeGreen
PUBLIC  :: Ctqmc_getGreen
PUBLIC  :: Ctqmc_getD
PUBLIC  :: Ctqmc_getE
PUBLIC  :: Ctqmc_printAll
PUBLIC  :: Ctqmc_printQMC
PUBLIC  :: Ctqmc_printGreen
PUBLIC  :: Ctqmc_printD
PUBLIC  :: Ctqmc_printE
PUBLIC  :: Ctqmc_printPerturbation
PUBLIC  :: Ctqmc_printCorrelation
PUBLIC  :: Ctqmc_printSpectra
PUBLIC  :: Ctqmc_destroy

CONTAINS
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_init
!! NAME
!!  Ctqmc_init
!!
!! FUNCTION
!!  Initialize the type Ctqmc
!!  Allocate all the non optional variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  ostream=where to write
!!  istream=where to read the input parameters if so
!!  bFile=do we read in istream ?
!!  MY_COMM=mpi communicator for the CTQMC
!!  iBuffer=input parameters if bFile is false
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

SUBROUTINE Ctqmc_init(this, ostream, istream, bFile, MY_COMM, iBuffer)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT)                      :: this
  INTEGER  , INTENT(IN   )                      :: ostream
  INTEGER  , INTENT(IN   )                      :: istream
  LOGICAL  , INTENT(IN   )                      :: bFile
  DOUBLE PRECISION, DIMENSION(1:9), OPTIONAL, INTENT(IN) :: iBuffer
  INTEGER  , OPTIONAL, INTENT(IN   )                      :: MY_COMM
!Local variables ------------------------------
#ifdef HAVE_MPI
  INTEGER                                       :: ierr
#endif
  INTEGER                                       :: iflavor
#ifdef __GFORTRAN__
!  INTEGER                                       :: pid
!  CHARACTER(LEN=5)                              :: Cpid
!
#endif
  DOUBLE PRECISION, DIMENSION(1:9)             :: buffer

  this%ostream = ostream
  this%istream = istream
  
! --- RENICE ---
!#ifdef __GFORTRAN__
!  pid = GetPid()
!  WRITE(Cpid,'(I5)') pid
!  CALL SYSTEM('renice +19 '//TRIM(ADJUSTL(Cpid))//' > /dev/null')
!#endif
!! --- RENICE ---

  IF ( PRESENT(MY_COMM)) THEN
#ifdef HAVE_MPI
    this%have_MPI = .TRUE.
    this%MY_COMM = MY_COMM
    CALL MPI_Comm_rank(this%MY_COMM, this%rank, ierr)
    CALL MPI_Comm_size(this%MY_COMM, this%size, ierr)
#else
    CALL WARN("Ctqmc_init : MPI is not used                                    ")
    this%have_MPI = .FALSE.
    this%MY_COMM = -1
    this%rank = 0
    this%size = 1
#endif
  ELSE
    this%have_MPI = .FALSE.
    this%MY_COMM = -1
    this%rank = 0
    this%size = 1
  END IF

  !IF ( this%rank .EQ. 0 ) THEN
!  WRITE(ostream,'(A20)') 'Job reniced with +19'
    !CALL FLUSH(ostream)
  !END IF

  IF ( bFile .EQV. .TRUE. ) THEN
    IF ( this%rank .EQ. 0 ) THEN

      READ(istream,*) buffer(1) !iseed
      READ(istream,*) buffer(2) !this%sweeps
      READ(istream,*) buffer(3) !this%thermalization
      READ(istream,*) buffer(4) !this%measurements
      READ(istream,*) buffer(5) !this%flavors
      READ(istream,*) buffer(6) !this%samples
      READ(istream,*) buffer(7) !this%beta
      READ(istream,*) buffer(8) !U
      READ(istream,*) buffer(9) !iTech
      !READ(istream,*) buffer(9) !Wmax
!#ifdef CTCtqmc_ANALYSIS
      !READ(istream,*) buffer(10) !order
!#endif
    END IF

#ifdef HAVE_MPI
    IF ( this%have_MPI .EQV. .TRUE. ) &
      CALL MPI_Bcast(buffer, 9, MPI_DOUBLE_PRECISION, 0,    &
                   this%MY_COMM, ierr)
#endif
  ELSE IF ( PRESENT(iBuffer) ) THEN
    buffer(1:9) = iBuffer(1:9)
  ELSE
    CALL ERROR("Ctqmc_init : No input parameters                    ")
  END IF

  CALL Ctqmc_setParameters(this, buffer)

  CALL Ctqmc_allocateAll(this)

  DO iflavor = 1, this%flavors
      CALL GreenHyb_init(this%Greens(iflavor),this%samples, this%beta, iTech=INT(buffer(9)),MY_COMM=this%MY_COMM)
  END DO


!  this%seg_added    = 0.d0
!  this%anti_added   = 0.d0
!  this%seg_removed  = 0.d0
!  this%anti_removed = 0.d0
!  this%seg_sign     = 0.d0
!  this%anti_sign    = 0.d0
  this%stats(:)     = 0.d0
  this%swap         = 0.d0
  this%runTime      = 0.d0

  CALL Vector_init(this%measNoise(1),this%sweeps/this%modNoise1)
  CALL Vector_init(this%measNoise(2),(this%sweeps/this%modNoise1+1)*CTQMC_SLICE2)
  !CALL Vector_init(this%measNoise(3),101)
  !CALL Vector_init(this%measNoise(4),101)

  this%set  = this%para .AND. this%inF
  this%done = .FALSE.
  this%init = .TRUE.

!#ifdef CTCtqmc_CHECK
  this%errorImpurity = 0.d0
  this%errorBath     = 0.d0
!#endif
END SUBROUTINE Ctqmc_init
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_setParameters
!! NAME
!!  Ctqmc_setParameters
!!
!! FUNCTION
!!  set all parameters and operators
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  buffer=input parameters
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

SUBROUTINE Ctqmc_setParameters(this,buffer)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT)                         :: this
  DOUBLE PRECISION, DIMENSION(1:9), INTENT(IN   ) :: buffer


  this%thermalization = INT(buffer(3)) !this%thermalization
  CALL Ctqmc_setSeed(this,INT(buffer(1)))
  CALL Ctqmc_setSweeps(this,buffer(2))

  this%measurements   = INT(buffer(4)) !this%measurements
  this%flavors        = INT(buffer(5)) !this%flavors
  this%samples        = INT(buffer(6)) !this%samples
  this%beta           = buffer(7)      !this%beta
  this%U              = buffer(8)      !U
!  this%mu             = buffer(9)      !this%mu
  !this%Wmax           = INT(buffer(9)) !Freq
!#ifdef CTCtqmc_ANALYSIS
!  this%order          = INT(buffer(10)) ! order
  this%inv_dt         = this%samples / this%beta
!#endif

  !CALL ImpurityOperator_init(this%Impurity,this%flavors,this%beta, this%samples)
  CALL ImpurityOperator_init(this%Impurity,this%flavors,this%beta)
  IF ( this%U .GE. 0.d0 ) THEN
    CALL ImpurityOperator_computeU(this%Impurity,this%U,0.d0)
    this%setU = .TRUE.
  END IF
!  this%mu = this%mu + this%Impurity%shift_mu

  CALL BathOperator_init(this%Bath, this%flavors, this%samples, this%beta, INT(buffer(9)))

  this%para = .TRUE.

END SUBROUTINE Ctqmc_setParameters
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_setSweeps
!! NAME
!!  Ctqmc_setSweeps
!!
!! FUNCTION
!!  set the number of sweeps
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  sweeps=asked sweeps
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

SUBROUTINE Ctqmc_setSweeps(this,sweeps)

!Arguments ------------------------------------
  TYPE(Ctqmc)         , INTENT(INOUT) :: this
  DOUBLE PRECISION  , INTENT(IN   ) :: sweeps

  this%sweeps = NINT(sweeps / DBLE(this%size))
!  write(6,*) this%sweeps,NINT(sweeps / DBLE(this%size)),ANINT(sweeps/DBLE(this%size))
  IF ( DBLE(this%sweeps) .NE. ANINT(sweeps/DBLE(this%size)) ) &
    CALL ERROR("Ctqmc_setSweeps : sweeps is negative or too big     ")
  IF ( this%sweeps .LT. 2*CTQMC_SLICE1 ) THEN  !202
    CALL WARNALL("Ctqmc_setSweeps : # sweeps automtically changed     ")
    this%sweeps = 2*CTQMC_SLICE1
!  ELSE IF ( this%sweeps .LT. this%thermalization ) THEN
!    CALL WARNALL("Ctqmc_setSweeps : Thermalization > sweeps / cpu -> auto fix")
!    this%sweeps = this%thermalization
  END IF
  IF ( DBLE(NINT(DBLE(this%sweeps)*DBLE(this%size)/DBLE(CTQMC_SLICE1))) .NE.  &
  ANINT(DBLE(this%sweeps)*DBLE(this%size)/DBLE(CTQMC_SLICE1)) ) THEN
    this%modNoise1 = this%sweeps
  ELSE
    this%modNoise1    = MIN(this%sweeps,INT(DBLE(this%sweeps)*DBLE(this%size) / DBLE(CTQMC_SLICE1))) !101
  END IF
  this%modNoise2    = MAX(this%modNoise1 / CTQMC_SLICE2, 1)   ! 100
!  this%modGlobalMove(1) = this%thermalization / 10 + 1
!  this%modGlobalMove(2) = 0

END SUBROUTINE Ctqmc_setSweeps
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_setSeed
!! NAME
!!  Ctqmc_setSeed
!!
!! FUNCTION
!!  initialize random number generator
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  iseed=seed from imput
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

SUBROUTINE Ctqmc_setSeed(this,iseed)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT)           :: this
  INTEGER  , INTENT(IN   )           :: iseed
!Local variables ------------------------------
  !INTEGER                            :: n
  !INTEGER                            :: i
  !INTEGER, DIMENSION(:), ALLOCATABLE :: seed


  !CALL RANDOM_SEED(size = n)
  !MALLOC(seed,(n))
  !seed =  iseed + (/ (i - 1, i = 1, n) /)

  !CALL RANDOM_SEED(PUT = seed+this%rank)

  !FREE(seed)

  this%seed=INT(iseed+this%rank,8)

END SUBROUTINE Ctqmc_setSeed
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_allocateAll
!! NAME
!!  Ctqmc_allocateAll
!!
!! FUNCTION
!!  Allocate all non option varibales
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_allocateAll(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                  :: flavors

  IF ( .NOT. this%para ) &
    CALL ERROR("Ctqmc_allocateAll : Ctqmc_setParameters never called  ")

  flavors = this%flavors

  DT_FREEIF(this%Greens)
  DT_MALLOC(this%Greens,(1:flavors))

  FREEIF(this%measN)
  MALLOC(this%measN,(1:4,1:flavors))
  this%measN = 0.d0

  FREEIF(this%measDE)
  MALLOC(this%measDE,(1:flavors,1:flavors) )
  this%measDE = 0.d0

  FREEIF(this%mu)
  MALLOC(this%mu,(1:flavors) )
  this%mu = 0.d0
END SUBROUTINE Ctqmc_allocateAll
!!***

!#ifdef CTCtqmc_ANALYSIS
!!****f* ABINIT/m_Ctqmc/Ctqmc_allocateOpt
!! NAME
!!  Ctqmc_allocateOpt
!!
!! FUNCTION
!!  allocate all option variables 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_allocateOpt(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: j
  INTEGER :: k

  IF ( .NOT. this%para ) &
    CALL ERROR("Ctqmc_allocateOpt : Ctqmc_setParameters never called  ")

  IF ( this%opt_analysis .EQ. 1 ) THEN
    FREEIF(this%measCorrelation)
    MALLOC(this%measCorrelation,(1:this%samples+1,1:3,1:this%flavors))
    this%measCorrelation = 0.d0
  END IF

  IF ( this%opt_order .GT. 0 ) THEN
    FREEIF(this%measPerturbation)
    MALLOC(this%measPerturbation,(1:this%opt_order,1:this%flavors))
    this%measPerturbation = 0.d0
  END IF

  IF ( this%opt_noise .EQ. 1 ) THEN
    IF ( ALLOCATED(this%measNoiseG) ) THEN
      DO i=1,2
        DO j = 1, this%flavors
          DO k= 1, this%samples+1
            CALL Vector_destroy(this%measNoiseG(k,j,i))
          END DO
        END DO
      END DO
      DT_FREE(this%measNoiseG)
    END IF
    DT_MALLOC(this%measNoiseG,(1:this%samples+1,1:this%flavors,1:2))
    !DO i=1,2
      DO j = 1, this%flavors
        DO k= 1, this%samples+1
          CALL Vector_init(this%measNoiseG(k,j,1),CTQMC_SLICE1)
        END DO
      END DO
      DO j = 1, this%flavors
        DO k= 1, this%samples+1
          CALL Vector_init(this%measNoiseG(k,j,2),CTQMC_SLICE1*CTQMC_SLICE2+1) ! +1 pour etre remplacer ceil
        END DO
      END DO
    !END DO
    FREEIF(this%abNoiseG)
    MALLOC(this%aBNoiseG,(1:2,1:this%samples+1,this%flavors))
    this%abNoiseG = 0.d0
  END IF

  IF (this%opt_spectra .GE. 1 ) THEN
    FREEIF(this%density)
    !MALLOC(this%density,(1:this%thermalization,1:this%flavors))
    i = CEILING(DBLE(this%thermalization+this%sweeps)/DBLE(this%measurements*this%opt_spectra))
    MALLOC(this%density,(1:this%flavors+1,1:i))
    this%density = 0.d0
  END IF
!#endif
END SUBROUTINE Ctqmc_allocateOpt
!!***

SUBROUTINE Ctqmc_setG0wFile(this,istream,opt_fk)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
  TYPE(Ctqmc)                    , INTENT(INOUT) :: this
  INTEGER                        , INTENT(IN   ) :: istream
  INTEGER                        , INTENT(IN   ) :: opt_fk
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: Gomega
  INTEGER                                        :: flavors
  INTEGER                                        :: iflavor
  INTEGER                                        :: iomega
#ifdef HAVE_MPI
  INTEGER                                        :: ierr
#endif

  IF ( .NOT. this%para ) &
    CALL ERROR("Ctqmc_setG0wFile : Ctqmc_setParameters never called   ") 

  flavors = this%flavors

  this%Wmax = this%samples
  MALLOC(Gomega,(1:this%Wmax,1:flavors))

  IF ( this%rank .EQ. 0 ) THEN
    DO iomega=1, this%Wmax
      READ(istream,*) (Gomega(iomega,iflavor),iflavor=1,flavors)
    END DO
  END IF

#ifdef HAVE_MPI
  CALL MPI_Bcast(Gomega, this%Wmax*flavors, MPI_DOUBLE_COMPLEX, 0,    &
                 this%MY_COMM, ierr)
#endif

  CALL Ctqmc_setG0wTab(this,Gomega,opt_fk)
  FREE(Gomega)

END SUBROUTINE Ctqmc_setG0wFile
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_setG0wTab
!! NAME
!!  Ctqmc_setG0wTab
!!
!! FUNCTION
!!  Set Gow from input array
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  Gomega=G0w
!!  opt_fk=F is already inversed with out iwn
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

SUBROUTINE Ctqmc_setG0wTab(this,Gomega,opt_fk)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT)                      :: this
  COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN ) :: Gomega
  INTEGER                         , INTENT(IN ) :: opt_fk
!Local variable -------------------------------
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F

  IF ( .NOT. this%para ) &
    CALL ERROR("Ctqmc_setG0wTab : Ctqmc_setParameters never called    ") 

  MALLOC(F,(1:this%samples+1,1:this%flavors))
  CALL Ctqmc_computeF(this,Gomega, F, opt_fk)  ! mu is changed
  CALL BathOperator_setF(this%Bath, F)
  !CALL BathOperator_printF(this%Bath)
  FREE(F)

  IF ( this%opt_levels .NE. 1 ) THEN ! We compute the mu by hand in computeF
    CALL ImpurityOperator_setMu(this%Impurity,this%mu)
  END IF

  this%inF = .TRUE.
  this%set = .TRUE. 

END SUBROUTINE Ctqmc_setG0wTab
!!***

!SUBROUTINE Ctqmc_setFwK(this,Gomega)
!  COMPLEX*16      , DIMENSION(:,:), INTENT(IN ) :: Gomega
!  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
!
!  IF ( .NOT. this%para ) &
!    CALL ERROR("Ctqmc_setG0wTab : Ctqmc_setParameters never called    ") 
!
!  MALLOC(F,(1:this%samples+1,1:this%flavors))
!  CALL Ctqmc_computeFK(this,Gomega, this%Wmax, F)  ! mu is changed
!  CALL BathOperator_setF(this%Bath, F)
!  CALL BathOperator_printF(this%Bath)
!  FREE(F)
!
!  this%inF = .TRUE.
!  this%set = .TRUE. 
!
!END SUBROUTINE Ctqmc_setFwK
!!***

!SUBROUTINE Ctqmc_setBand(this, mu, U)
!  DOUBLE PRECISION, INTENT(IN   ) :: mu
!  DOUBLE PRECISION, INTENT(IN   ) :: U
!
!  IF ( .NOT. this%para ) &
!    CALL ERROR("Ctqmc_setBand : Ctqmc_setParameters never called      ")
!
!  CALL ImpurityOperator_setU(this%Impurity, U, 0.d0)
!  !this%mu = mu + this%Impurity%shift_mu
!END SUBROUTINE Ctqmc_setBand
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_setU
!! NAME
!!  Ctqmc_setU
!!
!! FUNCTION
!!  set the interaction this
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_setU(this,matU)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: matU

  IF ( SIZE(matU) .NE. this%flavors*this%flavors ) &
    CALL ERROR("Ctqmc_setU : Wrong interaction this (size)        ")

  CALL ImpurityOperator_setUmat(this%Impurity, matU)
  this%setU = .TRUE.
END SUBROUTINE Ctqmc_setU
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_clear
!! NAME
!!  Ctqmc_clear
!!
!! FUNCTION
!!  clear a ctqmc run
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_clear(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: j
  INTEGER :: k

  this%measN(1,:) = 0.d0
  this%measN(2,:) = 0.d0
  !Do not set measN(3,:) to 0 to avoid erasing N between therm and ctqmc
  this%measN(4,:) = 0.d0
  this%measDE = 0.d0
!  this%seg_added    = 0.d0
!  this%anti_added   = 0.d0
!  this%seg_removed  = 0.d0
!  this%anti_removed = 0.d0
!  this%seg_sign     = 0.d0
!  this%anti_sign    = 0.d0
  this%stats(:)     = 0.d0
  this%swap         = 0.d0
  this%runTime      = 0.d0
  this%modGlobalMove(2) = 0 
  CALL Vector_clear(this%measNoise(1))
  CALL Vector_clear(this%measNoise(2))
!#ifdef CTCtqmc_CHECK
  this%errorImpurity = 0.d0
  this%errorBath     = 0.d0
!#endif
  DO j = 1, this%flavors
    CALL GreenHyb_clear(this%Greens(j))
  END DO
!#ifdef CTCtqmc_ANALYSIS
  IF ( this%opt_analysis .EQ. 1 .AND. ALLOCATED(this%measCorrelation) ) &    
    this%measCorrelation = 0.d0 
  IF ( this%opt_order .GT. 0 .AND. ALLOCATED(this%measPerturbation) ) &
    this%measPerturbation = 0.d0
  IF ( this%opt_noise .EQ. 1 .AND. ALLOCATED(this%measNoiseG) ) THEN
    DO i=1,2
      DO j = 1, this%flavors
        DO k= 1, this%samples+1
          CALL Vector_clear(this%measNoiseG(k,j,i))
        END DO
      END DO
    END DO
  END IF
!#endif
END SUBROUTINE Ctqmc_clear
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_reset
!! NAME
!!  Ctqmc_reset
!!
!! FUNCTION
!!  reset a ctqmc simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_reset(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                  :: iflavor
  DOUBLE PRECISION         :: sweeps

  DO iflavor = 1, this%flavors
    CALL GreenHyb_reset(this%Greens(iflavor))
  END DO
  CALL Ctqmc_clear(this)
  CALL ImpurityOperator_reset(this%Impurity)
  CALL BathOperator_reset    (this%Bath)
  this%measN(3,:) = 0.d0
  !complete restart -> measN=0
  this%done = .FALSE.
  this%set  = .FALSE.
  this%inF  = .FALSE.
  this%opt_movie = 0
  this%opt_analysis = 0
  this%opt_order = 0
  this%opt_check = 0
  this%opt_noise = 0
  this%opt_spectra = 0
  this%opt_levels = 0
  sweeps = DBLE(this%sweeps)*DBLE(this%size)
  CALL Ctqmc_setSweeps(this, sweeps)
!#ifdef HAVE_MPI
!  CALL MPI_BARRIER(this%MY_COMM,iflavor)
!  IF ( this%rank .EQ. 0 ) &
!#endif
!  WRITE(this%ostream,'(A9)') "QMC reset"
!  CALL FLUSH(this%ostream)
END SUBROUTINE Ctqmc_reset
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_setMu
!! NAME
!!  Ctqmc_setMu
!!
!! FUNCTION
!!  impose energy levels
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  levels=energy levels vector
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

SUBROUTINE Ctqmc_setMu(this, levels)

!Arguments ------------------------------------
  TYPE(Ctqmc)                   , INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN   ) :: levels

  IF ( this%flavors .NE. SIZE(levels,1) ) &
    CALL WARNALL("Ctqmc_setMu : Taking energy levels from weiss G(iw)")

  this%mu(:)=-levels  ! levels = \epsilon_j - \mu
  !this%mu =\tilde{\mu} = \mu -\epsilon_j
  CALL ImpurityOperator_setMu(this%Impurity,this%mu)
  this%opt_levels = 1
END SUBROUTINE Ctqmc_setMu
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_computeF
!! NAME
!!  Ctqmc_computeF
!!
!! FUNCTION
!!  Compute the hybridization function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  Gomega=G0 to compute F
!!  opt_fk=What is Gomega
!!
!! OUTPUT
!!  F=hybridization function
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

SUBROUTINE Ctqmc_computeF(this, Gomega, F, opt_fk)

!Arguments ------------------------------------
  TYPE(Ctqmc)                       , INTENT(INOUT) :: this
  COMPLEX(KIND=8), DIMENSION(:,:), INTENT(IN   ) :: Gomega
  !INTEGER                         , INTENT(IN   ) :: Wmax
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: F
  INTEGER                         , INTENT(IN   ) :: opt_fk
!Local variables ------------------------------
  INTEGER                                         :: flavors
  INTEGER                                         :: samples
  INTEGER                                         :: iflavor
  INTEGER                                         :: iomega
  INTEGER                                         :: itau
  DOUBLE PRECISION                                :: pi_invBeta
  DOUBLE PRECISION                                :: K
  DOUBLE PRECISION                                :: re
  DOUBLE PRECISION                                :: im
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: F_omega
  TYPE(GreenHyb)                                     :: F_tmp

  flavors    = this%flavors

  samples    = this%samples
  pi_invBeta = ACOS(-1.d0) / this%beta
  this%Wmax=SIZE(Gomega,1)

  IF ( this%have_MPI .EQV. .TRUE. ) THEN
    CALL GreenHyb_init(F_tmp,samples,this%beta,MY_COMM=this%MY_COMM)
  ELSE
    CALL GreenHyb_init(F_tmp,samples,this%beta)
  END IF
!  K = this%mu

  MALLOC(F_omega,(1:this%Wmax,1:flavors))

  !IF ( this%rank .EQ. 0 ) &
    !OPEN(UNIT=9876,FILE="K.dat",POSITION="APPEND")
  IF ( opt_fk .EQ. 0 ) THEN
    DO iflavor = 1, flavors
      DO iomega=1,this%Wmax
        re = REAL(Gomega(iomega,iflavor))
        im = AIMAG(Gomega(iomega,iflavor))
        F_omega(iomega,iflavor) = CMPLX(-re/(re*re+im*im),im/(re*re+im*im),8)
      END DO
    END DO
    !F_omega = CMPLX(-1.d0,0,8)/Gomega
  ELSE
    F_omega = Gomega
  END IF

  DO iflavor = 1, flavors
    IF ( this%opt_levels .EQ. 1 ) THEN
      K = this%mu(iflavor)
    ELSE
      K = -REAL(F_omega(this%Wmax, iflavor))
!    this%mu = K
      this%mu(iflavor) = K 
    END IF
    !IF ( this%rank .EQ. 0 ) &
    !WRITE(9876,'(I4,2E22.14)') iflavor, K, REAL(-F_omega(this%Wmax, iflavor))
    !IF(this%rank .EQ.0) &
    !WRITE(this%ostream,*) "CTQMC K, this%mu = ",K,this%mu(iflavor)
    !WRITE(this%ostream,*) "CTQMC beta     = ",this%beta
    IF ( opt_fk .EQ. 0 ) THEN
      DO iomega = 1, this%Wmax
        re = REAL(F_omega(iomega,iflavor))
        im = AIMAG(F_omega(iomega,iflavor))
        F_omega(iomega,iflavor) = CMPLX(re + K, im + (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, 8)
        !if(iflavor==1.and.this%rank==0) then
          !write(224,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(F_omega(iomega,iflavor)),imag(F_omega(iomega,iflavor))
          !write(225,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(Gomega(iomega, iflavor)),imag(Gomega(iomega, iflavor))
        !end if 
      END DO
    ELSE
      DO iomega = 1, this%Wmax
        F_omega(iomega,iflavor) = F_omega(iomega,iflavor) &
                    + CMPLX(K, 0.d0, 8)
        !if(iflavor==1.and.this%rank==0) then
          !write(224,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(F_omega(iomega,iflavor)),imag(F_omega(iomega,iflavor))
          !write(225,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(Gomega(iomega, iflavor)),imag(Gomega(iomega, iflavor))
        !end if 
      END DO
    END IF
    K = REAL(CMPLX(0,(2.d0*DBLE(this%Wmax)-1.d0)*pi_invBeta,8)*F_omega(this%Wmax,iflavor))
    CALL GreenHyb_setMuD1(this%Greens(iflavor),this%mu(iflavor),K)
    CALL GreenHyb_setOperW(F_tmp,F_omega(:,iflavor))
    !CALL GreenHyb_backFourier(F_tmp,F_omega(:,iflavor))
    CALL GreenHyb_backFourier(F_tmp)
    F(1:samples+1,iflavor) = (/ (-F_tmp%oper(samples+1-itau),itau=0,samples) /)
  END DO
  FREE(F_omega)
  CALL GreenHyb_destroy(F_tmp)
END SUBROUTINE Ctqmc_computeF
!!***

!SUBROUTINE Ctqmc_computeFK(this, Gomega, Wmax, F)
!  COMPLEX*16      , DIMENSION(:,:), INTENT(IN   ) :: Gomega
!  INTEGER                         , INTENT(IN   ) :: Wmax
!  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: F
!  INTEGER                                         :: flavors
!  INTEGER                                         :: samples
!  INTEGER                                         :: iflavor
!  INTEGER                                         :: iomega
!  INTEGER                                         :: itau
!  DOUBLE PRECISION                                :: pi_invBeta
!  DOUBLE PRECISION                                :: K
!  COMPLEX*16      , DIMENSION(:,:), ALLOCATABLE   :: F_omega
!  TYPE(GreenHyb)                                     :: F_tmp
!
!  flavors    = this%flavors
!
!  samples    = this%samples
!  pi_invBeta = ACOS(-1.d0) / this%beta
!
!  IF ( this%have_MPI .EQV. .TRUE. ) THEN
!    CALL GreenHyb_init(F_tmp,samples,this%beta,this%MY_COMM)
!  ELSE
!    CALL GreenHyb_init(F_tmp,samples,this%beta)
!  END IF
!!  K = this%mu
!
!  MALLOC(F_omega,(1:Wmax,1:flavors))
!
!  DO iflavor = 1, flavors
!    K = REAL(Gomega(Wmax, iflavor))
!    WRITE(this%ostream,*) "CTQMC K, this%mu = ",K,this%mu
!    WRITE(this%ostream,*) "CTQMC beta     = ",this%beta
!    this%mu(iflavor) = K 
!    DO iomega = 1, Wmax
!      F_omega(iomega,iflavor) = Gomega(iomega,iflavor) &
!                  - CMPLX(K, 0.d0, 8)
!      !if(iflavor==1.and.this%rank==0) then
!        !write(224,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(F_omega(iomega,iflavor)),imag(F_omega(iomega,iflavor))
!        !write(225,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(Gomega(iomega, iflavor)),imag(Gomega(iomega, iflavor))
!      !end if 
!    END DO
!    CALL GreenHyb_backFourier(F_tmp,F_omega(:,iflavor))
!    F(1:samples+1,iflavor) = (/ (-F_tmp%oper(samples+1-itau),itau=0,samples) /)
!  END DO
!  FREE(F_omega)
!  CALL GreenHyb_destroy(F_tmp)
!END SUBROUTINE Ctqmc_computeFK
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_run
!! NAME
!!  Ctqmc_run
!!
!! FUNCTION
!!  set all options and run a simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  opt_order=maximal perturbation order to scope
!!  opt_movie=draw a movie of the simulation
!!  opt_analysis=compute correlation functions
!!  opt_check=check fast calculations
!!  opt_noise=compute noise for green function
!!  opt_spectra=fourier transform of the time evolution of the number of electrons
!!  opt_gMove=steps without global move
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

SUBROUTINE Ctqmc_run(this,opt_order,opt_movie,opt_analysis,opt_check,opt_noise,opt_spectra,opt_gMove)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT)           :: this
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_order
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_movie
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_analysis
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_check
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_noise
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_spectra
  INTEGER, OPTIONAL, INTENT(IN   )  :: opt_gMove
!Local variables ------------------------------
#ifdef HAVE_MPI
  INTEGER                            :: ierr
#endif
!#ifdef CTCtqmc_MOVIE
  INTEGER                            :: ilatex
  CHARACTER(LEN=4)                   :: Cchar
!#endif
  DOUBLE PRECISION                   :: estimatedTime

  IF ( .NOT. this%set  ) &
    CALL ERROR("Ctqmc_run : QMC not set up                          ")
  IF ( .NOT. this%setU ) &
    CALL ERROR("Ctqmc_run : QMC does not have a U this            ")


! OPTIONS of the run
  IF ( PRESENT( opt_check ) ) THEN
    this%opt_check = opt_check
    CALL ImpurityOperator_doCheck(this%Impurity,opt_check)
    CALL BathOperator_doCheck(this%Bath,opt_check)
  END IF
  IF ( PRESENT( opt_movie ) ) &
    this%opt_movie = opt_movie
  IF ( PRESENT( opt_analysis ) ) &
    this%opt_analysis = opt_analysis
  IF ( PRESENT ( opt_order ) ) &
    this%opt_order = opt_order 
  IF ( PRESENT ( opt_noise ) ) THEN
    this%opt_noise = opt_noise 
  END IF
  IF ( PRESENT ( opt_spectra ) ) &
    this%opt_spectra = opt_spectra

  this%modGlobalMove(1) = this%sweeps+1 ! No Global Move
  this%modGlobalMove(2) = 0
  IF ( PRESENT ( opt_gMove ) ) THEN
    IF ( opt_gMove .LE. 0 .OR. opt_gMove .GT. this%sweeps ) THEN
      this%modGlobalMove(1) = this%sweeps+1
      CALL WARNALL("Ctqmc_run : global moves option is <= 0 or > sweeps/cpu -> No global Moves")
    ELSE 
      this%modGlobalMove(1) = opt_gMove 
    END IF
  END IF

  CALL Ctqmc_allocateOpt(this)
  
!#ifdef CTCtqmc_MOVIE  
  ilatex = 0
  IF ( this%opt_movie .EQ. 1 ) THEN
    Cchar ="0000"
    WRITE(Cchar,'(I4)') this%rank 
    ilatex = 87+this%rank
    OPEN(UNIT=ilatex, FILE="Movie_"//TRIM(ADJUSTL(Cchar))//".tex")
    WRITE(ilatex,'(A)') "\documentclass{beamer}"
    WRITE(ilatex,'(A)') "\usepackage{color}"
    WRITE(ilatex,'(A)') "\setbeamersize{sidebar width left=0pt}"
    WRITE(ilatex,'(A)') "\setbeamersize{sidebar width right=0pt}"
    WRITE(ilatex,'(A)') "\setbeamersize{text width left=0pt}"
    WRITE(ilatex,'(A)') "\setbeamersize{text width right=0pt}"
    WRITE(ilatex,*) 
    WRITE(ilatex,'(A)') "\begin{document}"
    WRITE(ilatex,*) 
  END IF
!#endif

  IF ( this%rank .EQ. 0 ) THEN
    WRITE(this%ostream,'(A29)') "Starting QMC (Thermalization)"
  END IF
  
  !=================================
  ! STARTING THERMALIZATION 
  !=================================
  CALL Ctqmc_loop(this,this%thermalization,ilatex)
  !=================================
  ! ENDING   THERMALIZATION 
  !=================================

  estimatedTime = this%runTime
#ifdef HAVE_MPI
  CALL MPI_REDUCE(this%runTime, estimatedTime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             0, this%MY_COMM, ierr)
#endif

  IF ( this%rank .EQ. 0 ) THEN
    WRITE(this%ostream,'(A26,I6,A11)') "Thermalization done in    ", CEILING(estimatedTime), "    seconds"
    WRITE(this%ostream,'(A25,I7,A15,I5,A5)') "The QMC should run in    ", &
           CEILING(estimatedTime*DBLE(this%sweeps)/DBLE(this%thermalization)),&
                        "    seconds on ", this%size, " CPUs"
  END IF

  !=================================
  ! CLEANING CTQMC          
  !=================================
  CALL Ctqmc_clear(this)

  !=================================
  ! STARTING CTQMC          
  !=================================
  CALL Ctqmc_loop(this,this%sweeps,ilatex)
  !=================================
  ! ENDING   CTQMC          
  !=================================

  IF ( this%opt_movie .EQ. 1 ) THEN
    WRITE(ilatex,*) ""
    WRITE(ilatex,'(A14)') "\end{document}"
    CLOSE(ilatex)
  END IF

  this%done     = .TRUE.

END SUBROUTINE Ctqmc_run
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_loop
!! NAME
!!  Ctqmc_loop
!!
!! FUNCTION
!!  Definition the main loop of the CT-QMC
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  itotal=number of sweeps to perform : thermalization or sweeps
!!  ilatex=unit of file to write movie if so
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

SUBROUTINE Ctqmc_loop(this,itotal,ilatex)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT)         :: this
  INTEGER    , INTENT(IN   )         :: itotal
  INTEGER    , INTENT(IN   )         :: ilatex
!Local variables ------------------------------
  LOGICAL                            :: updated 
  LOGICAL                            :: updated_seg
  LOGICAL, DIMENSION(:), ALLOCATABLE :: updated_swap

  INTEGER                            :: flavors
  INTEGER                            :: measurements
  INTEGER                            :: modNoise1
  INTEGER                            :: modNoise2
  INTEGER                            :: modGlobalMove
  INTEGER                            :: sp1
  INTEGER                            :: itau   
  INTEGER                            :: ind
  INTEGER                            :: endDensity
  INTEGER                            :: indDensity
  INTEGER                            :: swapUpdate1
  INTEGER                            :: swapUpdate2
  INTEGER                            :: old_percent
  INTEGER                            :: new_percent
  INTEGER                            :: ipercent
  INTEGER                            :: iflavor
  INTEGER                            :: isweep

  DOUBLE PRECISION                   :: cpu_time1
  DOUBLE PRECISION                   :: cpu_time2
  DOUBLE PRECISION                   :: NRJ_old1
  DOUBLE PRECISION                   :: NRJ_old2
  DOUBLE PRECISION                   :: NRJ_new
  DOUBLE PRECISION                   :: total
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: gtmp_old1
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: gtmp_old2
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: gtmp_new

  CALL CPU_TIME(cpu_time1)

  flavors        = this%flavors
  measurements   = this%measurements
  modNoise1      = this%modNoise1
  modNoise2      = this%modNoise2
  modGlobalMove  = this%modGlobalMove(1)
  sp1            = this%samples+1

  old_percent    = 0

  MALLOC(updated_swap,(1:flavors))
  updated_swap(:) = .FALSE.

  NRJ_old1  = 0.d0
  NRJ_old2  = 0.d0
  NRJ_new   = 0.d0

  MALLOC(gtmp_new,(1,1))
  gtmp_new  = 0.d0
  MALLOC(gtmp_old1,(1,1))
  gtmp_old1 = 0.d0
  MALLOC(gtmp_old2,(1,1))
  gtmp_old2 = 0.d0

  endDensity = SIZE(this%density,2)

  IF ( this%opt_noise .GT. 0 ) THEN
    FREEIF(gtmp_new)
    MALLOC(gtmp_new,(1:sp1,1:flavors))
    FREEIF(gtmp_old1)
    MALLOC(gtmp_old1,(1:sp1,1:flavors))
    FREEIF(gtmp_old2)
    MALLOC(gtmp_old2,(1:sp1,1:flavors))
  END IF

  IF ( this%rank .EQ. 0 ) THEN
    WRITE(this%ostream, '(1x,103A)') &
    "|----------------------------------------------------------------------------------------------------|"
    WRITE(this%ostream,'(1x,A)', ADVANCE="NO") "|"
  END IF

  total = DBLE(itotal)

  indDensity = 1
  DO isweep = 1, itotal
    DO iflavor = 1, flavors
      ImpurityOperator_QuickActivation(this%Impurity,iflavor)
      BathOperator_QuickActivation(this%Bath,iflavor)

      CALL Ctqmc_tryAddRemove(this,updated_seg)
      updated = updated_seg .OR.  updated_swap(iflavor)
      updated_swap(iflavor) = .FALSE.

      CALL GreenHyb_measHybrid(this%Greens(iflavor), this%Bath%M(iflavor), this%Impurity%Particles(iflavor), updated)
      CALL Ctqmc_measN        (this, iflavor, updated)
      IF ( this%opt_analysis .EQ. 1 ) &
        CALL Ctqmc_measCorrelation (this, iflavor)
      IF ( this%opt_order .GT. 0 ) &
        CALL Ctqmc_measPerturbation(this, iflavor)
    END DO

    IF ( MOD(isweep,modGlobalMove) .EQ. 0 ) THEN
      CALL Ctqmc_trySwap(this,swapUpdate1, swapUpdate2)
      IF ( swapUpdate1 .NE. 0 .AND. swapUpdate2 .NE. 0 ) THEN
        updated_swap(swapUpdate1) = .TRUE.
        updated_swap(swapUpdate2) = .TRUE.
      END IF
    END IF
    
    IF ( MOD(isweep,measurements) .EQ. 0 ) THEN
      CALL ImpurityOperator_measDE(this%Impurity,this%measDE)
      IF ( this%opt_spectra .GE. 1 .AND. MOD(isweep,measurements*this%opt_spectra) .EQ. 0 ) THEN
        this%density(1:flavors,indDensity) = this%measN(3,1:flavors)
        indDensity = indDensity+1
      END IF
    END IF

    IF ( MOD(isweep, modNoise1) .EQ. 0 ) THEN
      !modNext = isweep + modNoise2
      NRJ_new = (SUM(this%measDE(:,:))-this%measDE(1,1))*0.5d0 ! double occupation, avoid stat with 0 for U=J=0
      CALL Vector_pushBack(this%measNoise(1),NRJ_new - NRJ_old1)
      NRJ_old1 = NRJ_new

      !! Try to limit accumulation error
      CALL ImpurityOperator_cleanOverlaps(this%Impurity)

      IF ( this%opt_noise .EQ. 1 ) THEN
        DO iflavor = 1, flavors
          DO ind = 1, this%Greens(iflavor)%this%tail
            itau = this%Greens(iflavor)%this%listINT(ind)
            gtmp_new(itau,iflavor) = this%Greens(iflavor)%oper(itau) & 
                        +this%Greens(iflavor)%this%listDBLE(ind)*DBLE(this%Greens(iflavor)%factor)
          END DO
          DO itau = 1, sp1
            CALL Vector_pushBack(this%measNoiseG(itau,iflavor,1), gtmp_new(itau,iflavor) - gtmp_old1(itau,iflavor))
            gtmp_old1(itau,iflavor) = gtmp_new(itau,iflavor)
          END DO
        END DO
      END IF
    END IF

    IF ( MOD(isweep,modNoise2) .EQ. 0 ) THEN
      NRJ_new = (SUM(this%measDE(:,:))-this%measDE(1,1))*0.5d0 ! double occupation, avoid stat with 0 for U=J=0
      CALL Vector_pushBack(this%measNoise(2),NRJ_new - NRJ_old2)
      NRJ_old2 = NRJ_new
      IF ( this%opt_noise .EQ. 1 ) THEN
        DO iflavor = 1, flavors
          DO ind = 1, this%Greens(iflavor)%this%tail
            itau = this%Greens(iflavor)%this%listINT(ind)
            gtmp_new(itau,iflavor) = this%Greens(iflavor)%oper(itau) & 
                        +this%Greens(iflavor)%this%listDBLE(ind)*this%Greens(iflavor)%factor
          END DO
          DO itau = 1, sp1
            CALL Vector_pushBack(this%measNoiseG(itau,iflavor,2), gtmp_new(itau,iflavor) - gtmp_old2(itau,iflavor))
            gtmp_old2(itau,iflavor) = gtmp_new(itau,iflavor)
          END DO
        END DO 
      END IF

      IF ( this%rank .EQ. 0 ) THEN 
        new_percent = CEILING(DBLE(isweep)*100.d0/DBLE(itotal))
        DO ipercent = old_percent+1, new_percent 
          WRITE(this%ostream,'(A)',ADVANCE="NO") "-"
        END DO
        old_percent = new_percent
      END IF
    END IF

    IF ( this%opt_movie .EQ. 1 ) THEN
      WRITE(ilatex,'(A11,I9)') "%iteration ", isweep
      CALL ImpurityOperator_printLatex(this%Impurity,ilatex,isweep)
    END IF

  END DO

  IF ( this%rank .EQ. 0 ) THEN
    DO ipercent = old_percent+1, 100
      WRITE(this%ostream,'(A)',ADVANCE="NO") "-"
    END DO
    WRITE(this%ostream,'(A)') "|"
  END IF
 
  FREE(gtmp_new)
  FREE(gtmp_old1)
  FREE(gtmp_old2)
  FREE(updated_swap)

  IF ( this%opt_spectra .GE. 1 .AND. itotal .EQ. this%sweeps ) THEN
    IF ( endDensity .NE. indDensity-1 ) THEN
      this%density(:,endDensity) = -1.d0
    END IF
  END IF

  CALL CPU_TIME(cpu_time2)

  this%runTime = (cpu_time2 - cpu_time1)*1.05d0 ! facteur arbitraire de correction
END SUBROUTINE Ctqmc_loop
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_tryAddRemove
!! NAME
!!  Ctqmc_tryAddRemove
!!
!! FUNCTION
!!  Try to add or remove a segment and an anti-segment
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!
!! OUTPUT
!!  updated=something changed
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

SUBROUTINE Ctqmc_tryAddRemove(this,updated)

!Arguments ------------------------------------
  TYPE(Ctqmc)             , INTENT(INOUT) :: this
!  TYPE(BathOperator)    , INTENT(INOUT) :: Bath 
!  TYPE(ImpurityOperator), INTENT(INOUT) :: Impurity 
  LOGICAL               , INTENT(  OUT) :: updated
!Local variables ------------------------------
  INTEGER                               :: position
  INTEGER         , DIMENSION(1:2)     :: nature ! -2 for antiseg and 1 for seg
  INTEGER                               :: i! -2 for antiseg and 1 for seg
  DOUBLE PRECISION                      :: action
  DOUBLE PRECISION                      :: beta
  DOUBLE PRECISION                      :: time1
  DOUBLE PRECISION                      :: time2
  DOUBLE PRECISION                      :: time_avail
  DOUBLE PRECISION                      :: det_ratio
  DOUBLE PRECISION                      :: imp_trace
  DOUBLE PRECISION                      :: signe
  DOUBLE PRECISION                      :: tail
  DOUBLE PRECISION, DIMENSION(1:2)      :: CdagC_1

  IF ( .NOT. this%set ) &
    CALL ERROR("Ctqmc_trySegment : QMC not set                       ")

  nature(1) = CTQMC_SEGME
  nature(2) = CTQMC_ANTIS
  beta      = this%beta

  updated = .FALSE.
  tail  = DBLE(this%Impurity%particles(this%Impurity%activeFlavor)%tail)


  DO i = 1, 2
    signe = SIGN(1.d0,DBLE(nature(i))) 

    !CALL RANDOM_NUMBER(action)
    CALL OurRng(this%seed,action)

    IF ( action .LT. .5d0 ) THEN ! Ajout de segment
      !CALL RANDOM_NUMBER(time1)
      CALL OurRng(this%seed,time1)
      time1 = time1 * beta
      time_avail = ImpurityOperator_getAvailableTime(this%Impurity,time1,position) * signe
      IF ( time_avail .GT. 0.d0 ) THEN
        !CALL RANDOM_NUMBER(time2)
        CALL OurRng(this%seed,time2)
        IF ( time2 .EQ. 0.d0 ) THEN
          CALL OurRng(this%seed,time2) ! Prevent null segment
        END IF
        time2     = time1 + time2 * time_avail
        CdagC_1(Cdag_) = ((1.d0+signe)*time1+(1.d0-signe)*time2)*0.5d0
        CdagC_1(C_   ) = ((1.d0+signe)*time2+(1.d0-signe)*time1)*0.5d0
        det_ratio = BathOperator_getDetAdd(this%Bath,CdagC_1,position,this%Impurity%particles(this%Impurity%activeFlavor))
        imp_trace = ImpurityOperator_getTraceAdd(this%Impurity,CdagC_1)
        !CALL RANDOM_NUMBER(time1)
        CALL OurRng(this%seed,time1)
        IF ( det_ratio*imp_trace .LT. 0.d0 ) THEN
          this%stats(nature(i)+CTQMC_DETSI) = this%stats(nature(i)+CTQMC_DETSI) + 1.d0
        END IF
        IF ( (time1 * (tail + 1.d0 )) &
             .LT. (beta * time_avail * det_ratio * imp_trace ) ) THEN
          CALL ImpurityOperator_add(this%Impurity,CdagC_1,position)
          CALL BathOperator_setMAdd(this%bath,this%Impurity%particles(this%Impurity%activeFlavor))
          this%stats(nature(i)+CTQMC_ADDED) = this%stats(nature(i)+CTQMC_ADDED)  + 1.d0
          updated = .TRUE. .OR. updated
          tail = tail + 1.d0
        END IF 
      END IF 

    ELSE ! Supprimer un segment
      IF ( tail .GT. 0.d0 ) THEN
        !CALL RANDOM_NUMBER(time1)
        CALL OurRng(this%seed,time1)
        position = INT(((time1 * tail) + 1.d0) * signe )
        time_avail = ImpurityOperator_getAvailedTime(this%Impurity,position)
        det_ratio  = BathOperator_getDetRemove(this%Bath,position)
        imp_trace  = ImpurityOperator_getTraceRemove(this%Impurity,position)
        !CALL RANDOM_NUMBER(time1)
        CALL OurRng(this%seed,time1)
        IF ( det_ratio * imp_trace .LT. 0.d0 ) THEN
          this%stats(nature(i)+CTQMC_DETSI) = this%stats(nature(i)+CTQMC_DETSI) + 1.d0
        END IF
        IF ( (time1 * beta * time_avail ) &
             .LT. (tail * det_ratio * imp_trace) ) THEN
          CALL ImpurityOperator_remove(this%Impurity,position)
          CALL BathOperator_setMRemove(this%Bath,this%Impurity%particles(this%Impurity%activeFlavor))
          this%stats(nature(i)+CTQMC_REMOV) = this%stats(nature(i)+CTQMC_REMOV)  + 1.d0
          updated = .TRUE. .OR. updated
          tail = tail -1.d0
        END IF
      END IF
    END IF
  END DO
END SUBROUTINE Ctqmc_tryAddRemove
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_trySwap
!! NAME
!!  Ctqmc_trySwap
!!
!! FUNCTION
!!  try a global move (swap to flavors)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!
!! OUTPUT
!!  flav_i=first flavor swaped
!!  flav_j=second flavor swaped
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

SUBROUTINE Ctqmc_trySwap(this,flav_i,flav_j)

!Arguments ------------------------------------
  TYPE(Ctqmc)           , INTENT(INOUT) :: this
!  TYPE(BathOperator)    , INTENT(INOUT) :: Bath 
!  TYPE(ImpurityOperator), INTENT(INOUT) :: Impurity 
  INTEGER               , INTENT(  OUT) :: flav_i
  INTEGER               , INTENT(  OUT) :: flav_j
!Local variables ------------------------------
  INTEGER :: flavor_i
  INTEGER :: flavor_j
  DOUBLE PRECISION :: rnd
  DOUBLE PRECISION :: lengthi
  DOUBLE PRECISION :: lengthj
  DOUBLE PRECISION :: overlapic1
  DOUBLE PRECISION :: overlapjc1
  DOUBLE PRECISION :: overlapic2
  DOUBLE PRECISION :: overlapjc2
  DOUBLE PRECISION :: detic1
  DOUBLE PRECISION :: detjc1
  DOUBLE PRECISION :: detic2
  DOUBLE PRECISION :: detjc2
  DOUBLE PRECISION :: det_ratio
  DOUBLE PRECISION :: local_ratio

  !CALL RANDOM_NUMBER(rnd)
  CALL OurRng(this%seed,rnd)
  flavor_i = NINT(rnd*DBLE(this%flavors-1.d0))+1
  !CALL RANDOM_NUMBER(rnd)
  CALL OurRng(this%seed,rnd)
  flavor_j = NINT(rnd*DBLE(this%flavors-1.d0))+1
  
  flav_i = 0
  flav_j = 0

  IF ( flavor_i .NE. flavor_j ) THEN
    ! On tente d'intervertir i et j
    ! Configuration actuelle :
    this%modGlobalMove(2) = this%modGlobalMove(2)+1
    detic1     = BathOperator_getDetF(this%Bath,flavor_i)
    detjc1     = BathOperator_getDetF(this%Bath,flavor_j)
    lengthi    = ImpurityOperator_measN(this%Impurity,flavor_i)
    lengthj    = ImpurityOperator_measN(this%Impurity,flavor_j)
    overlapic1 = ImpurityOperator_overlapFlavor(this%Impurity,flavor_i)
    overlapjc1 = ImpurityOperator_overlapFlavor(this%Impurity,flavor_j)
    ! Configuration nouvelle :
    detic2     = BathOperator_getDetF(this%Bath,flavor_i,this%Impurity%particles(flavor_j))
    detjc2     = BathOperator_getDetF(this%Bath,flavor_j,this%Impurity%particles(flavor_i))
    ! lengths unchanged
    overlapic2 = ImpurityOperator_overlapSwap(this%Impurity,flavor_i,flavor_j)
    overlapjc2 = ImpurityOperator_overlapSwap(this%Impurity,flavor_j,flavor_i)

!    IF ( detic1*detjc1 .EQ. detic2*detjc2 ) THEN
!      det_ratio = 1.d0
!    ELSE IF ( detic1*detjc1 .EQ. 0.d0 ) THEN
!      det_ratio = detic2*detjc2 ! evite de diviser par 0 si pas de segment
!    ELSE
      det_ratio = detic2*detjc2/(detic1*detjc1)
!    END IF
    local_ratio = DEXP(-overlapic2*overlapjc2+overlapic1*overlapjc1 &
                      +(lengthj-lengthi)*(this%mu(flavor_i)-this%mu(flavor_j)))
    ! Wloc = exp(muN-Uo)
    !CALL RANDOM_NUMBER(rnd)
    CALL OurRng(this%seed,rnd)
    IF ( rnd .LT. local_ratio*det_ratio ) THEN ! swap accepted
      CALL ImpurityOperator_swap(this%Impurity, flavor_i,flavor_j)
      CALL BathOperator_swap    (this%Bath    , flavor_i,flavor_j)
      this%swap = this%swap + 1.d0
      flav_i = flavor_i
      flav_j = flavor_j
!    ELSE
!      CALL WARN("Swap refused")
!      WRITE(this%ostream,'(6E24.14)') local_ratio, det_ratio, detic1, detjc1, detic2, detjc2
    END IF
  END IF

END SUBROUTINE Ctqmc_trySwap
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_measN
!! NAME
!!  Ctqmc_measN
!!
!! FUNCTION
!!  measure the number of electron
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  iflavor=which flavor to measure
!!  updated=something has changed since last time
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

SUBROUTINE Ctqmc_measN(this, iflavor, updated)

!Arguments ------------------------------------
  TYPE(Ctqmc)             , INTENT(INOUT)     :: this
  !TYPE(ImpurityOperator), INTENT(IN   )     :: impurity
  INTEGER               , INTENT(IN   )     :: iflavor
  LOGICAL               , INTENT(IN   )     :: updated

!  IF ( .NOT. this%set ) &
!    CALL ERROR("Ctqmc_measN : QMC not set                           ")

  
  IF ( updated .EQV. .TRUE. ) THEN
    this%measN(1,iflavor) = this%measN(1,iflavor) + this%measN(3,iflavor)*this%measN(4,iflavor)
    this%measN(2,iflavor) = this%measN(2,iflavor) + this%measN(4,iflavor)
    this%measN(3,iflavor) = ImpurityOperator_measN(this%impurity)
    this%measN(4,iflavor) = 1.d0
  ELSE
    this%measN(4,iflavor) = this%measN(4,iflavor) + 1.d0
  END IF
END SUBROUTINE Ctqmc_measN
!!***

!#ifdef CTCtqmc_ANALYSIS
!!****f* ABINIT/m_Ctqmc/Ctqmc_measCorrelation
!! NAME
!!  Ctqmc_measCorrelation
!!
!! FUNCTION
!!  measure all correlations in times for a flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  iflavor=the flavor to measure
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

SUBROUTINE Ctqmc_measCorrelation(this, iflavor)

!Arguments ------------------------------------
  TYPE(Ctqmc)             , INTENT(INOUT)       :: this
  !TYPE(ImpurityOperator), INTENT(IN   )       :: impurity
  INTEGER               , INTENT(IN   )       :: iflavor
!Local variables ------------------------------
  INTEGER                                     :: iCdag
  INTEGER                                     :: iCdagBeta
  INTEGER                                     :: iC
  INTEGER                                     :: index
  INTEGER                                     :: size
  DOUBLE PRECISION                            :: tC
  DOUBLE PRECISION                            :: tCdag
  !DOUBLE PRECISION                            :: time
  DOUBLE PRECISION                            :: inv_dt
  DOUBLE PRECISION                            :: beta

  IF ( .NOT. this%set ) &
    CALL ERROR("Ctqmc_measCorrelation : QMC not set                 ")

  size = this%impurity%particles(this%impurity%activeFlavor)%tail
  beta = this%beta

  IF ( size .EQ. 0 ) RETURN
  
  inv_dt = this%inv_dt

  DO iCdag = 1, size ! first segments
    tCdag  = this%impurity%particles(this%impurity%activeFlavor)%list(iCdag,Cdag_)
    tC     = this%impurity%particles(this%impurity%activeFlavor)%list(iCdag,C_   )
    index = INT( ( (tC - tCdag)  * inv_dt ) + .5d0 ) + 1
    this%measCorrelation(index,1,iflavor) = this%measCorrelation(index,1,iflavor) + 1.d0
    MODCYCLE(iCdag+1,size,iCdagBeta)
    index = INT( ( ( &
                    this%impurity%particles(this%impurity%activeFlavor)%list(iCdagBeta,Cdag_) - tC &
                    + AINT(DBLE(iCdag)/DBLE(size))*beta &
                   )  * inv_dt ) + .5d0 ) + 1
    IF ( index .LT. 1 .OR. index .GT. this%samples+1 ) THEN
      CALL WARN("Ctqmc_measCorrelation : bad index line 1095         ")
    ELSE
      this%measCorrelation(index,2,iflavor) = this%measCorrelation(index,2,iflavor) + 1.d0
    END IF
!    DO iC = 1, size
!      tC = impurity%particles(impurity%activeFlavor)%list(C_,iC)
!      time = tC - tCdag
!      IF ( time .LT. 0.d0 ) time = time + beta
!      index = INT( ( time * inv_dt ) + .5d0 ) + 1
!      this%measCorrelation(index,3,iflavor) = this%measCorrelation(index,3,iflavor) + 1.d0
!    END DO
    DO iC = 1, size!  this%Greens(iflavor)%index_old%tail 
        this%measCorrelation(this%Greens(iflavor)%this%listINT(iC+(iCdag-1)*size),3,iflavor) = &
        this%measCorrelation(this%Greens(iflavor)%this%listINT(iC+(iCdag-1)*size),3,iflavor) + 1.d0
    END DO
  END DO

END SUBROUTINE Ctqmc_measCorrelation
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_measPerturbation
!! NAME
!!  Ctqmc_measPerturbation
!!
!! FUNCTION
!!  measure perturbation order
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  iflavor=the flavor to measure
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

SUBROUTINE Ctqmc_measPerturbation(this, iflavor)

!Arguments ------------------------------------
  TYPE(Ctqmc)             , INTENT(INOUT)     :: this
  !TYPE(ImpurityOperator), INTENT(IN   )     :: impurity
  INTEGER               , INTENT(IN   )     :: iflavor
!Local variables ------------------------------
  INTEGER                                   :: index

  IF ( .NOT. this%set ) &
    CALL ERROR("Ctqmc_measiPerturbation : QMC not set               ")

  index = this%impurity%particles(this%impurity%activeFlavor)%tail + 1
  IF ( index .LE. this%opt_order ) &
    this%measPerturbation(index,iflavor) = this%measPerturbation(index,iflavor) + 1.d0

END SUBROUTINE Ctqmc_measPerturbation
!!***
!#endif

!!****f* ABINIT/m_Ctqmc/Ctqmc_getResult
!! NAME
!!  Ctqmc_getResult
!!
!! FUNCTION
!!  reduce everything to get the result of the simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_getResult(this)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(Ctqmc)  , INTENT(INOUT)                    :: this
!Local variables ------------------------------
  INTEGER                                       :: iflavor
  INTEGER                                       :: flavors
  INTEGER                                       :: itau
  INTEGER                                       :: endDensity
  DOUBLE PRECISION                              :: inv_flavors
  DOUBLE PRECISION                              :: a
  DOUBLE PRECISION                              :: b
  DOUBLE PRECISION                              :: r
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: alpha
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: beta
  DOUBLE PRECISION,              DIMENSION(1:2) :: TabX
  DOUBLE PRECISION,              DIMENSION(1:2) :: TabY
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: freqs
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: counts
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: displs
  INTEGER                                       :: sp1
  INTEGER                                       :: spAll
  INTEGER                                       :: last
  INTEGER                                       :: n1
  INTEGER                                       :: n2
  INTEGER                                       :: debut
!  INTEGER                                       :: fin
#ifdef HAVE_MPI
  INTEGER                                       :: ierr
#endif
  DOUBLE PRECISION                              :: inv_size
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: buffer 
  TYPE(FFTHyb) :: FFTmrka

  IF ( .NOT. this%done ) &
    CALL ERROR("Ctqmc_getResult : Simulation not run                ")

  flavors     =  this%flavors
  inv_flavors = 1.d0 / DBLE(flavors)


  inv_size = 1.d0 / DBLE(this%size)
  sp1 = 0
  spAll = 0

!#ifdef CTCtqmc_CHECK
  IF ( this%opt_check .GT. 0 ) THEN
    this%errorImpurity = ImpurityOperator_getError(this%Impurity) * inv_flavors 
    this%errorBath     = BathOperator_getError    (this%Bath    ) * inv_flavors 
  END IF
!#endif

  MALLOC(alpha,(1,1))
  MALLOC(beta,(1,1))
  MALLOC(buffer,(1,1))
  IF ( this%opt_noise .EQ. 1) THEN
    FREEIF(alpha)
    MALLOC(alpha,(1:this%samples+1,1:flavors))
    FREEIF(beta)
    MALLOC(beta,(1:this%samples+1,1:flavors))
  END IF

  IF ( this%have_MPI .EQV. .TRUE.) THEN 
    sp1   = this%samples+1
    spALL = sp1 + flavors + 6 

!#ifdef CTCtqmc_ANALYSIS
    IF ( this%opt_analysis .EQ. 1 ) &
      spAll = spAll + 3*sp1 
    IF ( this%opt_order .GT. 0 ) &
      spAll = spAll + this%opt_order 
    IF ( this%opt_noise .EQ. 1 ) &
      spAll = spAll + 2*(this%samples + 1)
!#endif

    FREEIF(buffer)
    MALLOC(buffer,(1:spAll,1:MAX(2,flavors)))
  END IF

!  this%seg_added    = this%seg_added    * inv_flavors 
!  this%seg_removed  = this%seg_removed  * inv_flavors
!  this%seg_sign     = this%seg_sign     * inv_flavors
!  this%anti_added   = this%anti_added   * inv_flavors
!  this%anti_removed = this%anti_removed * inv_flavors
!  this%anti_sign    = this%anti_sign    * inv_flavors
  this%stats(:) = this%stats(:) * inv_flavors

  DO iflavor = 1, flavors
    CALL GreenHyb_measHybrid(this%Greens(iflavor), this%Bath%M(iflavor), this%Impurity%Particles(iflavor), .TRUE.)
    CALL GreenHyb_getHybrid(this%Greens(iflavor))
    ! Accumule les dernieres mesure de N
    this%measN(1,iflavor) = this%measN(1,iflavor) + this%measN(3,iflavor)*this%measN(4,iflavor)
    this%measN(2,iflavor) = this%measN(2,iflavor) + this%measN(4,iflavor)
    ! Reduction
    this%measN(1,iflavor)  = this%measN(1,iflavor) / ( this%measN(2,iflavor) * this%beta )
    ! Correction
    CALL GreenHyb_setN(this%Greens(iflavor), this%measN(1,iflavor))
!#ifdef CTCtqmc_ANALYSIS
    IF ( this%opt_order .GT. 0 ) &
      this%measPerturbation(:   ,iflavor) = this%measPerturbation(:,iflavor) &
                                    / SUM(this%measPerturbation(:,iflavor))
    IF ( this%opt_analysis .EQ. 1 ) THEN
      this%measCorrelation (:,1,iflavor) = this%measCorrelation  (:,1,iflavor) &
                                    / SUM(this%measCorrelation (:,1,iflavor)) &
                                    * this%inv_dt 
      this%measCorrelation (:,2,iflavor) = this%measCorrelation  (:,2,iflavor) &
                                    / SUM(this%measCorrelation (:,2,iflavor)) &
                                    * this%inv_dt 
      this%measCorrelation (:,3,iflavor) = this%measCorrelation  (:,3,iflavor) &
                                    / SUM(this%measCorrelation (:,3,iflavor)) &
                                    * this%inv_dt 
    END IF
!#endif
    IF ( this%opt_noise .EQ. 1 ) THEN
      TabX(1) = DBLE(this%modNoise2)
      TabX(2) = DBLE(this%modNoise1)
      DO itau = 1, this%samples+1
        this%measNoiseG(itau,iflavor,2)%vec = -this%measNoiseG(itau,iflavor,2)%vec*this%inv_dt &  
                                           /(this%beta*DBLE(this%modNoise2))
        this%measNoiseG(itau,iflavor,1)%vec = -this%measNoiseG(itau,iflavor,1)%vec*this%inv_dt &  
                                           /(this%beta*DBLE(this%modNoise1))
        n2 = this%measNoiseG(itau,iflavor,2)%tail
        TabY(1) = Stat_deviation(this%measNoiseG(itau,iflavor,2)%vec(1:n2))!*SQRT(n2/(n2-1))
        n1 = this%measNoiseG(itau,iflavor,1)%tail
        TabY(2) = Stat_deviation(this%measNoiseG(itau,iflavor,1)%vec(1:n1))!*SQRT(n1/(n1-1))
        CALL Stat_powerReg(TabX,SQRT(2.d0*LOG(2.d0))*TabY,alpha(itau,iflavor),beta(itau,iflavor),r)
        ! ecart type -> 60%
        ! largeur a mi-hauteur d'une gaussienne -> sqrt(2*ln(2))*sigma
      END DO
    END IF

    IF ( this%have_MPI .EQV. .TRUE. ) THEN 
      buffer(1:sp1, iflavor) = this%Greens(iflavor)%oper(1:sp1)
    END IF
  END DO
  last = sp1

  this%measDE(:,:) = this%measDE(:,:) * DBLE(this%measurements) /(DBLE(this%sweeps)*this%beta)

  n1 = this%measNoise(1)%tail
  n2 = this%measNoise(2)%tail

  ! On utilise freqs comme tableau de regroupement
  ! Gather de Noise1
  IF ( this%have_MPI .EQV. .TRUE. ) THEN
    MALLOC(counts,(1:this%size))
    MALLOC(displs,(1:this%size))
    FREEIF(freqs)
    MALLOC(freqs,(1:this%size*n1))
    freqs = 0.d0
    freqs(n1*this%rank+1:n1*(this%rank+1)) = this%measNoise(1)%vec(1:n1) 
    counts(:) = n1
    displs(:) = (/ ( iflavor*n1, iflavor=0, this%size-1 ) /)
#ifdef HAVE_MPI
    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DOUBLE_PRECISION, &
                        freqs, counts, displs, &
                        MPI_DOUBLE_PRECISION, this%MY_COMM, ierr)
#endif
    n1 = this%size*n1
    CALL Vector_setSize(this%measNoise(1),n1)
    this%measNoise(1)%vec(1:n1) = freqs(:)
    ! Gather de Noise2
    FREE(freqs)
    MALLOC(freqs,(1:this%size*n2))
    freqs = 0.d0
    freqs(n2*this%rank+1:n2*(this%rank+1)) = this%measNoise(2)%vec(1:n2) 
    counts(:) = n2
    displs(:) = (/ ( iflavor*n2, iflavor=0, this%size-1 ) /)
#ifdef HAVE_MPI
    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        freqs, counts, displs, &
                        MPI_DOUBLE_PRECISION, this%MY_COMM, ierr)
#endif
    n2 = this%size*n2
    CALL Vector_setSize(this%measNoise(2),n2)
    this%measNoise(2)%vec(1:n2) = freqs(:)
    FREE(counts)
    FREE(displs)
    FREE(freqs)
  END IF
  !n1 = this%measNoise(1)%tail
  !n2 = this%measNoise(2)%tail

  ! Transformation des paquets pour que ca fit a CTQMC_SLICE(1|2)
  IF ( n1 .GT. CTQMC_SLICE1 ) THEN
    itau = n1/CTQMC_SLICE1
    MALLOC(freqs,(1:n1/itau))
    DO debut=1, n1/itau
      freqs(debut)=SUM(this%measNoise(1)%vec((debut-1)*itau+1:itau*debut))
    END DO
    freqs(:) = freqs(:)/DBLE(itau)
    this%modNoise1 = this%modNoise1*itau
    n1 = n1/itau
    CALL Vector_setSize(this%measNoise(1),n1)
    this%measNoise(1)%vec(1:n1) = freqs(:)
    FREE(freqs)
  END IF
  IF ( n2 .GT. CTQMC_SLICE1*CTQMC_SLICE2 ) THEN
    itau = n2/(CTQMC_SLICE1*CTQMC_SLICE2)
    MALLOC(freqs,(1:n2/itau))
    DO debut=1, n2/itau
      freqs(debut)=SUM(this%measNoise(2)%vec((debut-1)*itau+1:itau*debut))
    END DO
    freqs(:) = freqs(:)/DBLE(itau)
    this%modNoise2 = this%modNoise2*itau
    n2 = n2/itau
    CALL Vector_setSize(this%measNoise(2),n2)
    this%measNoise(2)%vec(1:n2) = freqs(:)
    FREE(freqs)
  END IF
  ! On peut s'amuser avec nos valeur d'energies
  !MALLOC(TabX,(1:20))
  !MALLOC(TabY,(1:20))

  TabX(1) = DBLE(this%modNoise2)
  TabX(2) = DBLE(this%modNoise1)

  ! Il faut calculer pour chaque modulo 10 ecarts type sur les donnes acquises
  this%measNoise(1)%vec(1:n1) = this%measNoise(1)%vec(1:n1)/(this%beta*DBLE(this%modNoise1))*DBLE(this%measurements)
  this%measNoise(2)%vec(1:n2) = this%measNoise(2)%vec(1:n2)/(this%beta*DBLE(this%modNoise2))*DBLE(this%measurements)
!  CALL Vector_print(this%measNoise(1),this%rank+70)
!  CALL Vector_print(this%measNoise(2),this%rank+50)
!  DO iflavor=1,10
!    debut = (iflavor-1)*n2/10+1
!    fin   = iflavor*n2/10
!    TabY(iflavor) = Stat_deviation(this%measNoise(2)%vec(debut:fin))
!    debut = (iflavor-1)*n1/10+1
!    fin   = iflavor*n1/10
!    TabY(10+iflavor) = Stat_deviation(this%measNoise(1)%vec(debut:fin))
!  END DO
!!  TabY(1:n) = (this%measNoise(2)%vec(1:n)   &
!!              )
!!             !/(this%beta*DBLE(this%modNoise2))*DBLE(this%measurements) &
!!             !- this%measDE(1,1))
!!  TabY(this%measNoise(2)%tail+1:n+this%measNoise(2)%tail) = (this%measNoise(1)%vec(1:n)   &
!!               )
!!             ! /(this%beta*DBLE(this%modNoise1))*DBLE(this%measurements) &
!!             ! - this%measDE(1,1))
!  IF ( this%rank .EQ. 0 ) THEN
!    DO iflavor=1,20
!      write(45,*) TabX(iflavor), TabY(iflavor)
!    END DO
!  END IF
!


  TabY(1) = Stat_deviation(this%measNoise(2)%vec(1:n2))!*SQRT(n2/(n2-1))
!!  write(this%rank+10,*) TabX(2)
!!  write(this%rank+40,*) TabX(1)
!!  CALL Vector_print(this%measNoise(1),this%rank+10)
!!  CALL Vector_print(this%measNoise(2),this%rank+40)
!!  CLOSE(this%rank+10)
!!  CLOSE(this%rank+40)
  TabY(2) = Stat_deviation(this%measNoise(1)%vec(1:n1))!*SQRT(n1/(n1-1))
!!  ! Ecart carre moyen ~ ecart type mais non biaise. Serait moins precis. Aucun
  ! impact sur la pente, juste sur l'ordonnee a l'origine.

  CALL Stat_powerReg(TabX,SQRT(2.d0*LOG(2.d0))*TabY,a,b,r)
!  FREE(TabX)
!  FREE(TabY)
  ! ecart type -> 60%
  ! largeur a mi-hauteur d'une gaussienne -> sqrt(2*ln(2))*sigma

  !this%measDE(1,1) = SUM(this%measNoise(1)%vec(1:this%measNoise(1)%tail))/(DBLE(this%measNoise(1)%tail*this%modNoise1)*this%beta)
  !this%measDE(2:flavors,1:flavors) = this%measDE(2:flavors,1:flavors) /(DBLE(this%sweeps)*this%beta)
  CALL ImpurityOperator_getErrorOverlap(this%Impurity,this%measDE)
  ! Add the difference between true calculation and quick calculation of the
  ! last sweep overlap to measDE(2,2)
  !this%measDE = this%measDE * DBLE(this%measurements) 
  IF ( this%have_MPI .EQV. .TRUE. ) THEN 
    IF ( this%opt_analysis .EQ. 1 ) THEN
      buffer(last+1:last+sp1,:) = this%measCorrelation(:,1,:)
      last = last + sp1
      buffer(last+1:last+sp1,:) = this%measCorrelation(:,2,:)
      last = last + sp1
      buffer(last+1:last+sp1,:) = this%measCorrelation(:,3,:)
      last = last + sp1
    END IF
    IF ( this%opt_order .GT. 0 ) THEN
      buffer(last+1:last+this%opt_order, :) = this%measPerturbation(:,:)
      last = last + this%opt_order
    END IF
    IF ( this%opt_noise .EQ. 1 ) THEN
      buffer(last+1:last+this%samples+1,:) = alpha(:,:)
      last = last + this%samples + 1
      buffer(last+1:last+this%samples+1,:) = beta(:,:)
      last = last + this%samples + 1
    END IF
!  this%measDE(2,2) = a*EXP(b*LOG(DBLE(this%sweeps*this%size)))
    buffer(spall-(flavors+5):spAll-6,:) = this%measDE(:,:)
!    buffer(spAll  ,1) = this%seg_added   
!    buffer(spAll-1,1) = this%seg_removed 
!    buffer(spAll-2,1) = this%seg_sign    
!    buffer(spAll  ,2) = this%anti_added  
!    buffer(spAll-1,2) = this%anti_removed
!    buffer(spAll-2,2) = this%anti_sign   
    buffer(spAll  ,1) = this%stats(1)
    buffer(spAll-1,1) = this%stats(2)
    buffer(spAll-2,1) = this%stats(3)
    buffer(spAll  ,2) = this%stats(4)
    buffer(spAll-1,2) = this%stats(5)
    buffer(spAll-2,2) = this%stats(6)
    buffer(spAll-3,1) = this%swap
    buffer(spAll-3,2) = DBLE(this%modGlobalMove(2))
    buffer(spAll-4,1) = a
    buffer(spAll-4,2) = b
!#ifdef CTCtqmc_CHECK
    buffer(spAll-5,1) = this%errorImpurity
    buffer(spAll-5,2) = this%errorBath 
!#endif

#ifdef HAVE_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, buffer, spAll*flavors, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, this%MY_COMM, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%runTime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             this%MY_COMM, ierr)
#endif

  
    buffer          = buffer * inv_size
    this%measDE(:,:)  = buffer(spall-(flavors+5):spAll-6,:)
!    this%seg_added    = buffer(spAll  ,1)
!    this%seg_removed  = buffer(spAll-1,1)
!    this%seg_sign     = buffer(spAll-2,1)
!    this%anti_added   = buffer(spAll  ,2)
!    this%anti_removed = buffer(spAll-1,2)
!    this%anti_sign    = buffer(spAll-2,2)
    this%stats(1)    = buffer(spAll  ,1)
    this%stats(2)    = buffer(spAll-1,1)
    this%stats(3)    = buffer(spAll-2,1)
    this%stats(4)    = buffer(spAll  ,2)
    this%stats(5)    = buffer(spAll-1,2)
    this%stats(6)    = buffer(spAll-2,2)
    this%swap         = buffer(spAll-3,1)
    this%modGlobalMove(2) = NINT(buffer(spAll-3,2))
    a               = buffer(spAll-4,1) 
    b               = buffer(spAll-4,2)
!!#ifdef CTCtqmc_CHECK
    this%errorImpurity= buffer(spAll-5,1) 
    this%errorBath    = buffer(spAll-5,2)   
!#endif

    DO iflavor = 1, flavors
      this%Greens(iflavor)%oper          = buffer(1:sp1          , iflavor)
    END DO
    last = sp1
    IF ( this%opt_analysis .EQ. 1 ) THEN
      this%measCorrelation(:,1,:) = buffer(last+1:last+sp1,:) 
      last = last + sp1
      this%measCorrelation(:,2,:) = buffer(last+1:last+sp1,:) 
      last = last + sp1
      this%measCorrelation(:,3,:) = buffer(last+1:last+sp1,:) 
      last = last + sp1
    END IF
    IF ( this%opt_order .GT. 0 ) THEN
      this%measPerturbation(:,:) = buffer(last+1:last+this%opt_order, :)
      last = last + this%opt_order
    END IF
    IF ( this%opt_noise .EQ. 1 ) THEN
      alpha(:,:) = buffer(last+1:last+this%samples+1,:)
      last = last + this%samples + 1
      beta(:,:) = buffer(last+1:last+this%samples+1,:)
      last = last + this%samples + 1
    END IF
  END IF
  DO iflavor = 1, flavors
    ! complete DE this
    this%measDE(iflavor, iflavor+1:flavors) = this%measDE(iflavor+1:flavors,iflavor)
  END DO
  FREE(buffer)

  IF ( this%opt_spectra .GE. 1 ) THEN
    endDensity = SIZE(this%density,2)
    IF ( this%density(1,endDensity) .EQ. -1.d0 ) &
      endDensity = endDensity - 1
    CALL FFTHyb_init(FFTmrka,endDensity,DBLE(this%thermalization)/DBLE(this%measurements*this%opt_spectra))
    ! Not very Beauty 
    MALLOC(freqs,(1:FFTmrka%size/2))
    DO iflavor = 1, flavors
      ! mean value is removed to supress the continue composent 
      CALL FFTHyb_setData(FFTmrka,this%density(iflavor,1:endDensity)/this%beta+this%Greens(iflavor)%oper(this%samples+1))
      CALL FFTHyb_run(FFTmrka,1)
      CALL FFTHyb_getData(FFTmrka,endDensity,this%density(iflavor,:),freqs)
    END DO
    this%density(flavors+1,:) = -1.d0
    this%density(flavors+1,1:FFTmrka%size/2) = freqs
    CALL FFTHyb_destroy(FFTmrka)
    FREE(freqs)
  END IF

  this%a_Noise = a
  this%b_Noise = b
  IF ( this%opt_noise .EQ. 1 ) THEN
    this%abNoiseG(1,:,:) = alpha
    this%abNoiseG(2,:,:) = beta
  END IF
  FREE(alpha)
  FREE(beta)

END SUBROUTINE Ctqmc_getResult
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_symmetrizeGreen
!! NAME
!!  Ctqmc_symmetrizeGreen
!!
!! FUNCTION
!!  optionnaly symmetrize the green functions
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  syms=weight factors
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

SUBROUTINE Ctqmc_symmetrizeGreen(this, syms)

!Arguments ------------------------------------
  TYPE(Ctqmc)                     , INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN   ) :: syms
!Local variables ------------------------------
  INTEGER :: iflavor1
  INTEGER :: iflavor2
  INTEGER :: flavors
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: green_tmp
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:  ) :: n_tmp

  flavors = this%flavors
  IF ( SIZE(syms,1) .NE. flavors .OR. SIZE(syms,2) .NE. flavors ) THEN
    CALL WARNALL("Ctqmc_symmetrizeGreen : wrong opt_sym -> not symmetrizing")
    RETURN
  END IF
 
  MALLOC(green_tmp,(1:this%samples+1,flavors))
  green_tmp(:,:) = 0.d0
  MALLOC(n_tmp,(1:flavors))
  n_tmp(:) = 0.d0
  DO iflavor1=1, flavors
    DO iflavor2=1,flavors
      green_tmp(:,iflavor1) = green_tmp(:,iflavor1) &
                             + syms(iflavor2,iflavor1) * this%Greens(iflavor2)%oper(:)
      n_tmp(iflavor1) = n_tmp(iflavor1) &
                             + syms(iflavor2,iflavor1) * this%measN(1,iflavor2)
    END DO
  END DO
  DO iflavor1=1, flavors
    this%Greens(iflavor1)%oper(:) = green_tmp(:,iflavor1)
    this%measN(1,iflavor1)          = n_tmp(iflavor1)
  END DO
  FREE(green_tmp)
  FREE(n_tmp)
END SUBROUTINE Ctqmc_symmetrizeGreen
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_getGreen
!! NAME
!!  Ctqmc_getGreen
!!
!! FUNCTION
!!  Get the full green functions in time and/or frequency
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!
!! OUTPUT
!!  Gtau=green function in time
!!  Gw=green function in frequency
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

SUBROUTINE Ctqmc_getGreen(this, Gtau, Gw)

!Arguments ------------------------------------
  TYPE(Ctqmc)          , INTENT(INOUT)    :: this
  DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: Gtau
  COMPLEX(KIND=8), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: Gw
!Local variables ------------------------------
  !INTEGER                            :: itime
  INTEGER                            :: iflavor1
  INTEGER                            :: iflavor2
  INTEGER                            :: iflavor3
  INTEGER                            :: flavors
  DOUBLE PRECISION :: u1 
  DOUBLE PRECISION :: u2
  DOUBLE PRECISION :: Un
  DOUBLE PRECISION :: UUnn

  flavors = this%flavors
  DO iflavor1 = 1, flavors
    u1 = 0.d0
    u2 = 0.d0
    DO iflavor2 = 1, flavors
      IF ( iflavor2 .EQ. iflavor1 ) CYCLE
      Un = this%Impurity%mat_U(iflavor2,iflavor1) * this%measN(1,iflavor2)
      u1 = u1 + Un 
      u2 = u2 + Un*this%Impurity%mat_U(iflavor2,iflavor1) 
      DO iflavor3 = 1, flavors
        IF ( iflavor3 .EQ. iflavor2 .OR. iflavor3 .EQ. iflavor1 ) CYCLE
        UUnn = (this%Impurity%mat_U(iflavor2,iflavor1)*this%Impurity%mat_U(iflavor3,iflavor1)) * this%measDE(iflavor2,iflavor3) 
        u2 = u2 + UUnn 
      END DO
    END DO  

    CALL GreenHyb_setMoments(this%Greens(iflavor1),u1,u2)
    IF ( PRESENT( Gtau ) ) THEN
      Gtau(1:this%samples,iflavor1) = this%Greens(iflavor1)%oper(1:this%samples)
    END IF
       !write(6,*) "present gw", present(gw)
    IF ( PRESENT( Gw ) ) THEN
       !write(6,*) "size gw",SIZE(Gw,DIM=2) ,flavors+1 
      IF ( SIZE(Gw,DIM=2) .EQ. flavors+1 ) THEN
        CALL GreenHyb_forFourier(this%Greens(iflavor1), Gomega=Gw(:,iflavor1), omega=Gw(:,this%flavors+1))
        !IF ( this%rank .EQ. 0 ) write(20,*) Gw(:,iflavor1)
      ELSE IF ( SIZE(Gw,DIM=2) .EQ. flavors ) THEN  
        CALL GreenHyb_forFourier(this%Greens(iflavor1),Gomega=Gw(:,iflavor1))
      ELSE
        CALL WARNALL("Ctqmc_getGreen : Gw is not valid                    ")
        CALL GreenHyb_forFourier(this%Greens(iflavor1),Wmax=this%Wmax)
      END IF
    ELSE
      CALL GreenHyb_forFourier(this%Greens(iflavor1),Wmax=this%Wmax)
    END IF
  END DO
END SUBROUTINE Ctqmc_getGreen
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_getD
!! NAME
!!  Ctqmc_getD
!!
!! FUNCTION
!!  get double occupation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!
!! OUTPUT
!!  D=full double occupation
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

SUBROUTINE Ctqmc_getD(this, D)

!Arguments ------------------------------------
  TYPE(Ctqmc)       , INTENT(IN ) :: this
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: D
!Local variables ------------------------------
  INTEGER                       :: iflavor1
  INTEGER                       :: iflavor2

  IF ( SIZE(D,1) .LT. this%flavors .OR. SIZE(D,2) .LT. this%flavors ) &
    CALL ERROR("Ctqmc_getD : Dimensions of array D are too small")

  D = 0.d0

  DO iflavor1 = 1, this%flavors
    DO iflavor2 = 1, this%flavors
      D(iflavor2,iflavor1) =  this%measDE(iflavor2,iflavor1)
    END DO
    D(iflavor1,iflavor1) = 0.d0
  END DO

END SUBROUTINE Ctqmc_getD
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_getE
!! NAME
!!  Ctqmc_getE
!!
!! FUNCTION
!!  get interaction energy and noise on it
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!
!! OUTPUT
!!  E=interaction energy
!!  noise=noise on this value
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

SUBROUTINE Ctqmc_getE(this,E,noise)

!Arguments ------------------------------------
  TYPE(Ctqmc)       , INTENT(IN ) :: this
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: E
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: Noise

  IF ( PRESENT(E) ) &
    E = this%measDE(1,1)  
  IF ( PRESENT(Noise) ) &
    Noise = SUM(this%Impurity%mat_U)/(this%flavors*(this%flavors-1)) &
            * this%a_Noise*(DBLE(this%sweeps)*DBLE(this%size))**this%b_Noise
END SUBROUTINE Ctqmc_getE
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_printAll
!! NAME
!!  Ctqmc_printAll
!!
!! FUNCTION
!!  print different functions computed during the simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_printAll(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this

  IF ( .NOT. this%done ) &
    CALL WARNALL("Ctqmc_printAll : Simulation not run                 ")

  CALL Ctqmc_printQMC(this)

  CALL Ctqmc_printGreen(this)

  CALL Ctqmc_printD(this)

!  CALL Ctqmc_printE(this)

!#ifdef CTCtqmc_ANALYSIS
  CALL Ctqmc_printPerturbation(this)

  CALL Ctqmc_printCorrelation(this)
!#endif

  CALL Ctqmc_printSpectra(this)

END SUBROUTINE Ctqmc_printAll
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_printQMC
!! NAME
!!  Ctqmc_printQMC
!!
!! FUNCTION
!!  print ctqmc statistics
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_printQMC(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                  :: ostream
  INTEGER                  :: iflavor
  DOUBLE PRECISION         :: sweeps
  DOUBLE PRECISION         :: invSweeps
  CHARACTER(LEN=2)         :: a
  CHARACTER(LEN=15)        :: string

  !IF ( this%rank .NE. 0) RETURN
  IF ( this%rank .NE. MOD(this%size,this%size)) RETURN

  ostream   = this%ostream
  sweeps    = DBLE(this%sweeps)
  invSweeps = 1.d0/sweeps

  WRITE(ostream,'(1x,F13.0,A11,F10.2,A12,I5,A5)') sweeps*DBLE(this%size), " sweeps in ", this%runTime, &
                 " seconds on ", this%size, " CPUs"
  WRITE(ostream,'(A28,F6.2)') "Segments added        [%] : ", this%stats(CTQMC_SEGME+CTQMC_ADDED)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Segments removed      [%] : ", this%stats(CTQMC_SEGME+CTQMC_REMOV)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Segments sign         [%] : ", this%stats(CTQMC_SEGME+CTQMC_DETSI)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Anti-segments added   [%] : ", this%stats(CTQMC_ANTIS+CTQMC_ADDED)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Anti-segments removed [%] : ", this%stats(CTQMC_ANTIS+CTQMC_REMOV)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Anti-segments sign    [%] : ", this%stats(CTQMC_ANTIS+CTQMC_DETSI)*invSweeps*100.d0
  IF ( this%modGlobalMove(1) .LT. this%sweeps + 1 ) THEN
    WRITE(ostream,'(A28,F6.2)') "Global Move           [%] : ", this%swap         *invSweeps*100.d0*this%modGlobalMove(1)
    WRITE(ostream,'(A28,F6.2)') "Global Move Reduced   [%] : ", this%swap         / DBLE(this%modGlobalMove(2))*100.d0
  END IF
!#ifdef CTCtqmc_CHECK
  IF ( this%opt_check .EQ. 1 .OR. this%opt_check .EQ. 3 ) &
    WRITE(ostream,'(A28,E22.14)') "Impurity test         [%] : ", this%errorImpurity*100.d0
  IF ( this%opt_check .GE. 2 ) &
      WRITE(ostream,'(A28,E22.14)') "Bath     test         [%] : ", this%errorBath    *100.d0
!#endif
  WRITE(ostream,'(A28,ES22.14,A5,ES21.14)') "<Epot>                [U] : ", this%measDE(1,1), " +/- ",&
!#ifdef HAVE_MPI
         SUM(this%Impurity%mat_U)/(this%flavors*(this%flavors-1)) * this%a_Noise*(sweeps*DBLE(this%size))**this%b_Noise
!#else
!                                                              this%a_Noise*(sweeps)**this%b_Noise
!#endif
  WRITE(ostream,'(A28,F8.4,A3,F7.4)') "Noise                [/U] : ", this%a_Noise, " x^", this%b_Noise
  WRITE(ostream,'(A28,E10.2)')  "Niquist puls.     [/beta] : ", ACOS(-1.d0)*this%inv_dt
  WRITE(ostream,'(A28,E22.14)') "Max Acc. Epot Error   [U] : ", this%measDE(2,2)/(this%beta*this%modNoise1*2.d0)*sweeps
  
  !WRITE(ostream,'(A28,F7.4,A3,F7.4,A4,E20.14)') "Noise            [G(tau)] : ", this%a_Noise(2), "x^", this%b_Noise(2), " -> ", &
                                                              !this%a_Noise(2)*(sweeps*DBLE(this%size))**this%b_Noise(2)
  IF ( this%opt_order .GT. 0 ) THEN 
    WRITE(a,'(I2)') this%flavors
    string = '(A28,'//TRIM(ADJUSTL(a))//'(1x,I3))'
    WRITE(ostream,string) "Perturbation orders       : ", &
      (/ (MAXLOC(this%measPerturbation(:, iflavor))-1, iflavor=1, this%flavors) /)
  END IF
  !CALL FLUSH(this%ostream)
  IF ( ABS(((this%stats(CTQMC_SEGME+CTQMC_ADDED) *invSweeps*100.d0) / &
            (this%stats(CTQMC_SEGME+CTQMC_REMOV) *invSweeps*100.d0) - 1.d0)) .GE. 0.02d0 &
   .OR. ABS(((this%stats(CTQMC_ANTIS+CTQMC_ADDED)*invSweeps*100.d0) / &
             (this%stats(CTQMC_ANTIS+CTQMC_REMOV)*invSweeps*100.d0) - 1.d0)) .GE. 0.02d0 ) &
    THEN 
    CALL WARNALL("Ctqmc_printQMC : bad statistic according to moves. Increase sweeps")
  END IF
  ! Check sign problem for diagonal hybridization.
  IF ( (this%stats(CTQMC_SEGME+CTQMC_DETSI) + this%stats(CTQMC_ANTIS+CTQMC_DETSI)) .GT. 1.d-10 ) THEN
    CALL WARNALL("Ctqmc_printQMC : at least one negative sign occured. There might be a bug in the CT-QMC")
  END IF

  IF ( ABS(this%b_Noise+0.5)/0.5d0 .GE. 0.05d0 ) &
    CALL WARNALL("Ctqmc_printQMC : bad statistic according to Noise. Increase sweeps")
!  IF ( ISNAN(this%a_Noise) .OR. ISNAN(this%a_Noise) ) &
!    CALL WARNALL("Ctqmc_printQMC : NaN appeared. Increase sweeps    ")


END SUBROUTINE Ctqmc_printQMC
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_printGreen
!! NAME
!!  Ctqmc_printGreen
!!
!! FUNCTION
!!  print green functions
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  oFileIn=file stream
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

SUBROUTINE Ctqmc_printGreen(this, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmc)        , INTENT(IN)    :: this
  INTEGER  , OPTIONAL, INTENT(IN)    :: oFileIn
!Local variables ------------------------------
  INTEGER                            :: oFile
  INTEGER                            :: itime
  INTEGER                            :: sp1
  INTEGER                            :: iflavor
  INTEGER                            :: flavors
  CHARACTER(LEN=4)                   :: cflavors
  CHARACTER(LEN=50)                  :: string
  DOUBLE PRECISION                   :: dt
  DOUBLE PRECISION                   :: sweeps

  !IF ( this%rank .NE. MOD(1,this%size)) RETURN
  IF ( this%rank .NE. MOD(this%size+1,this%size)) RETURN

  oFile = 40
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Gtau.dat")
  END IF

  sp1     =  this%samples
  dt      =  this%beta / DBLE(sp1)
  sp1     =  sp1 + 1
  flavors =  this%flavors
  sweeps = DBLE(this%sweeps)*DBLE(this%size)

  IF ( this%opt_noise .EQ. 1) THEN
    WRITE(cflavors,'(I4)') 2*flavors+1
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'ES22.14)'
    DO itime = 1, sp1
      WRITE(oFile,string) DBLE(itime-1)*dt, &
      (/ (this%Greens(iflavor)%oper(itime), iflavor=1, flavors) /), &
      (/ (this%abNoiseG(1,itime,iflavor)*(sweeps)**this%abNoiseG(2,itime,iflavor), iflavor=1, flavors) /)
    END DO
  ELSE
    WRITE(cflavors,'(I4)') flavors+1
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'ES22.14)'
    DO itime = 1, sp1
      WRITE(oFile,string) DBLE(itime-1)*dt, &
      (/ (this%Greens(iflavor)%oper(itime), iflavor=1, flavors) /)
    END DO
  END IF

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmc_printGreen
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_printD
!! NAME
!!  Ctqmc_printD
!!
!! FUNCTION
!!  print individual double occupancy
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  oFileIn=file stream
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

SUBROUTINE Ctqmc_printD(this,oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmc)          , INTENT(IN)    :: this
  INTEGER  , OPTIONAL, INTENT(IN)    :: oFileIn
!Local variables ------------------------------
  INTEGER                            :: oFile
  INTEGER                            :: iflavor1
  INTEGER                            :: iflavor2

  !IF ( this%rank .NE. MOD(2,this%size)) RETURN
  IF ( this%rank .NE. MOD(this%size+2,this%size)) RETURN

  oFile = 41
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="D.dat")
  END IF

  DO iflavor1 = 1, this%flavors
    DO iflavor2 = iflavor1+1, this%flavors
      WRITE(oFile,'(1x,A8,I4,A1,I4,A3,ES21.14)') "Orbitals", iflavor1, "-", iflavor2, " : ", this%measDE(iflavor2,iflavor1)
    END DO
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmc_printD
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_printE
!! NAME
!!  Ctqmc_printE
!!
!! FUNCTION
!!  print energy and noise 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  oFileIn=file stream
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

SUBROUTINE Ctqmc_printE(this,oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmc)          , INTENT(IN)    :: this
  INTEGER  , OPTIONAL, INTENT(IN)    :: oFileIn
!Local variables ------------------------------
  INTEGER                            :: oFile
  DOUBLE PRECISION                   :: E
  DOUBLE PRECISION                   :: Noise

  !IF ( this%rank .NE. MOD(3,this%size)) RETURN
  IF ( this%rank .NE. MOD(this%size+3,this%size)) RETURN

  oFile = 42
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="BetaENoise.dat")
  END IF

  CALL Ctqmc_getE(this,E,Noise)

  WRITE(oFile,'(1x,F5.2,A2,ES21.14,A2,ES21.14)') this%beta, "  ", E, "  ",  Noise

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmc_printE
!!***

!#ifdef CTCtqmc_ANALYSIS
!!****f* ABINIT/m_Ctqmc/Ctqmc_printPerturbation
!! NAME
!!  Ctqmc_printPerturbation
!!
!! FUNCTION
!!  print perturbation order
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  oFileIn=file stream
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

SUBROUTINE Ctqmc_printPerturbation(this, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmc)          , INTENT(IN)           :: this
  INTEGER  , OPTIONAL,  INTENT(IN)          :: oFileIn
!Local variables-------------------------------
  INTEGER                                   :: oFile
  INTEGER                                   :: iorder
  INTEGER                                   :: order
  INTEGER                                   :: iflavor
  INTEGER                                   :: flavors
  CHARACTER(LEN=2)                          :: a
  CHARACTER(LEN=50)                         :: string

  !IF ( this%rank .NE. MOD(4,this%size)) RETURN
  IF ( this%rank .NE. MOD(this%size+4,this%size)) RETURN
  IF ( this%opt_order .LE. 0 ) RETURN

  oFile = 43
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Perturbation.dat")
  END IF
    
  order        =  this%opt_order
  flavors      =  this%flavors

  WRITE(a,'(I2)') flavors
  string = '(I5,'//TRIM(ADJUSTL(a))//'F19.15)'
  DO iorder = 1, order
    WRITE(oFile,string) iorder-1, &
                (/ (this%measPerturbation(iorder, iflavor), iflavor=1, flavors) /)
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)
END SUBROUTINE Ctqmc_printPerturbation
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_printCorrelation
!! NAME
!!  Ctqmc_printCorrelation
!!
!! FUNCTION
!!  print correlation fonctions
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  oFileIn=file stream
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

SUBROUTINE Ctqmc_printCorrelation(this, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmc)          , INTENT(IN)             :: this
  INTEGER  , OPTIONAL, INTENT(IN)             :: oFileIn
!Local variables ------------------------------
  INTEGER                                     :: oFile
  INTEGER                                     :: itime
  INTEGER                                     :: sp1
  INTEGER                                     :: iflavor
  INTEGER                                     :: i
  INTEGER                                     :: flavors
  CHARACTER(LEN=2)                            :: a
  CHARACTER(LEN=50)                           :: string
  DOUBLE PRECISION                            :: dt

  !IF ( this%rank .NE. MOD(5,this%size)) RETURN
  IF ( this%rank .NE. MOD(this%size+5,this%size)) RETURN
  IF ( this%opt_analysis .NE. 1 ) RETURN

  oFile = 44
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Correlation.dat")
  END IF

  sp1         =  this%samples
  dt          =  this%beta / sp1
  sp1         =  sp1 + 1
  flavors     =  this%flavors

  i = 3*flavors + 1
  WRITE(a,'(I2)') i
  WRITE(oFile,*) "# time  (/ (segement, antiseg, correl), i=1, flavor/)"
  string = '(1x,'//TRIM(ADJUSTL(a))//'F19.15)'
  DO itime = 1, sp1
    WRITE(oFile,string) DBLE(itime-1)*dt, &
                   (/ ( &
                   (/ ( this%measCorrelation(itime, i, iflavor), i=1,3) /) &
                   , iflavor=1, flavors) /)
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmc_printCorrelation
!!***
!#endif

!!****f* ABINIT/m_Ctqmc/Ctqmc_printSpectra
!! NAME
!!  Ctqmc_printSpectra
!!
!! FUNCTION
!!  print fourier transform of time evolution of number of electrons
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
!!  oFileIn=file stream
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

SUBROUTINE Ctqmc_printSpectra(this, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmc)          , INTENT(IN)             :: this
  INTEGER  , OPTIONAL, INTENT(IN)             :: oFileIn
!Local variables ------------------------------
  INTEGER                                     :: oFile
  INTEGER                                     :: flavors
  INTEGER                                     :: indDensity
  INTEGER                                     :: endDensity
  CHARACTER(LEN=4)                            :: a
  CHARACTER(LEN=16)                           :: formatSpectra

  !IF ( this%rank .NE. MOD(6,this%size)) RETURN
  IF ( this%opt_spectra .LT. 1 ) RETURN

  oFile = 45+this%rank
  a ="0000"
  WRITE(a,'(I4)') this%rank
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Markov_"//TRIM(ADJUSTL(a))//".dat")
  END IF

  flavors     =  this%flavors
  WRITE(a,'(I4)') flavors+1
  formatSpectra ='(1x,'//TRIM(ADJUSTL(a))//'ES22.14)'
  WRITE(oFile,*) "# freq[/hermalization] FFT"

  endDensity = SIZE(this%density,2)
  DO WHILE ( this%density(flavors+1,endDensity) .EQ. -1 )
    endDensity = endDensity -1
  END DO

  DO indDensity = 1, endDensity
    WRITE(oFile,formatSpectra) this%density(flavors+1,indDensity), this%density(1:flavors,indDensity)
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmc_printSpectra
!!***

!!****f* ABINIT/m_Ctqmc/Ctqmc_destroy
!! NAME
!!  Ctqmc_destroy
!!
!! FUNCTION
!!  destroy and deallocate all variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmc
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

SUBROUTINE Ctqmc_destroy(this)

!Arguments ------------------------------------
  TYPE(Ctqmc), INTENT(INOUT) :: this
!Local variables ------------------------------
  INTEGER                  :: iflavor
  INTEGER                  :: flavors
  INTEGER                  :: i
  INTEGER                  :: j
  INTEGER                  :: k

  if ( this%init .EQV. .FALSE. ) RETURN

  flavors = this%flavors

  CALL ImpurityOperator_destroy(this%Impurity)
  CALL BathOperator_destroy(this%Bath)
  CALL Vector_destroy(this%measNoise(1))
  CALL Vector_destroy(this%measNoise(2))

  IF ( ALLOCATED(this%Greens) ) THEN
    DO iflavor = 1, flavors
     CALL GreenHyb_destroy(this%Greens(iflavor))
    END DO
    DT_FREE( this%Greens )
  END IF
!#ifdef CTCtqmc_ANALYSIS
  FREEIF(this%measCorrelation)
  FREEIF(this%measPerturbation)
  FREEIF(this%measN)
  FREEIF(this%measDE)
  FREEIF(this%mu)
  FREEIF(this%abNoiseG)
  IF ( ALLOCATED(this%measNoiseG) ) THEN
    DO i=1,2
      DO j = 1, this%flavors
        DO k= 1, this%samples+1
          CALL Vector_destroy(this%measNoiseG(k,j,i))
        END DO
      END DO
    END DO
    DT_FREE(this%measNoiseG)
  END IF
  FREEIF(this%density)
!#endif
  this%ostream        = 0
  this%istream        = 0
 
  this%sweeps         = 0
  this%thermalization = 0
  this%flavors        = 0
  this%samples        = 0
  this%beta           = 0.d0
!  this%seg_added      = 0.d0
!  this%anti_added     = 0.d0
!  this%seg_removed    = 0.d0
!  this%anti_removed   = 0.d0
!  this%seg_sign       = 0.d0
!  this%anti_sign      = 0.d0
  this%stats          = 0.d0
  this%swap           = 0.d0


  this%set  = .FALSE.
  this%done = .FALSE.
  this%init = .FALSE.
END SUBROUTINE Ctqmc_destroy
!!***

END MODULE m_Ctqmc
!!***

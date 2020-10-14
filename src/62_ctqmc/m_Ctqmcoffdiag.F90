
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_Ctqmcoffdiag
!! NAME
!!  m_Ctqmcoffdiag
!! 
!! FUNCTION 
!!  Manage and drive all the CTQMC
!!  Should not be used if you don't know what you do
!!  Please use CtqmcoffdiagInterface
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder, B. Amadon, J. Denier)
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
#define CTQMC_SLICE1 100
! Coupe Sweeps en 100
#define CTQMC_SLICE2 100
! Coupe modNoise1 en 2000
#define CTQMC_SEGME  1
#define CTQMC_ANTIS -2
#define CTQMC_ADDED  3  
#define CTQMC_REMOV  4
#define CTQMC_DETSI  5
MODULE m_Ctqmcoffdiag

 USE m_Global
 USE m_GreenHyboffdiag
 USE m_BathOperatoroffdiag
 USE m_ImpurityOperator
 USE m_Stat
 USE m_FFTHyb
 USE m_OurRng
#ifdef HAVE_MPI2
 USE mpi
#endif
 IMPLICIT NONE

 public :: Ctqmcoffdiag_init
 public :: Ctqmcoffdiag_setParameters
 public :: Ctqmcoffdiag_setSweeps
 public :: Ctqmcoffdiag_setSeed
 public :: Ctqmcoffdiag_allocateAll
 public :: Ctqmcoffdiag_allocateOpt
 public :: Ctqmcoffdiag_setG0wTab
 public :: Ctqmcoffdiag_setU
 public :: Ctqmcoffdiag_clear
 public :: Ctqmcoffdiag_reset
 public :: Ctqmcoffdiag_setMu
 public :: Ctqmcoffdiag_computeF
 public :: Ctqmcoffdiag_run
 public :: Ctqmcoffdiag_tryAddRemove
 public :: Ctqmcoffdiag_trySwap
 public :: Ctqmcoffdiag_measN
 public :: Ctqmcoffdiag_measCorrelation
 public :: Ctqmcoffdiag_measPerturbation
 public :: Ctqmcoffdiag_getResult
 public :: Ctqmcoffdiag_symmetrizeGreen
 public :: Ctqmcoffdiag_getGreen
 public :: Ctqmcoffdiag_getD
 public :: Ctqmcoffdiag_getE
 public :: Ctqmcoffdiag_printAll
 public :: Ctqmcoffdiag_printQMC
 public :: Ctqmcoffdiag_printGreen
 public :: Ctqmcoffdiag_printD
 public :: Ctqmcoffdiag_printE
 public :: Ctqmcoffdiag_printPerturbation
 public :: Ctqmcoffdiag_printCorrelation
 public :: Ctqmcoffdiag_printSpectra
 public :: Ctqmcoffdiag_destroy

!!***

!!****t* m_Ctqmcoffdiag/Ctqmcoffdiag
!! NAME
!!  Ctqmcoffdiag
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

TYPE Ctqmcoffdiag

  LOGICAL :: init = .FALSE.
! Flag: is MC initialized

  LOGICAL :: set  = .FALSE.
! Flag: ??

  LOGICAL :: setU = .FALSE.
! Flag: is U Set ?

  LOGICAL :: inF  = .FALSE.
! Flag: is hybridization fct in input ?

  LOGICAL :: done = .FALSE.
! Flag: is MC terminated ?

  LOGICAL :: para = .FALSE.
! Flag:  do we have parameters in input

  LOGICAL :: have_MPI = .FALSE.
! Flag: 

  INTEGER :: opt_movie = 0
!

  INTEGER :: opt_analysis = 0
! correlations 

  INTEGER :: opt_check = 0
! various check 0
! various check 1 impurity
! various check 2 bath
! various check 3 both

  INTEGER :: opt_order = 0
! nb of segments max for analysis

  INTEGER :: opt_noise = 0
! compute noise

  INTEGER :: opt_spectra = 0
! markov chain FT (correlation time)

  INTEGER :: opt_levels = 0
! do we have energy levels

  INTEGER :: opt_hybri_limit = 0
! do we have limit of hybridization (yes=1)

  INTEGER :: opt_nondiag = 0
! if opt_nondiag = 1 F is non diagonal.

  INTEGER :: prtopt = 1
! printing

  INTEGER :: flavors
! number of flavors 

  INTEGER :: measurements
!  The modulo used to measure the interaction energy and the number of electrons. Example : 2 means the measure is perform every two sweeps. 

  INTEGER :: samples
! nb of L points (dmftqmc_l)

  INTEGER(8) :: seed
!

  INTEGER :: sweeps
!

  INTEGER :: thermalization
!

  INTEGER :: ostream
! output file

  INTEGER :: istream
! input file

  INTEGER :: modNoise1
! measure the noise each modNoise1

  INTEGER :: modNoise2
! measure the noise each modNoise2

  INTEGER :: activeFlavor
! orbital on which one do sth now

  INTEGER, DIMENSION(1:2) :: modGlobalMove
! 1: global move each modglobalmove(1)
! 2: we have done modglobalmove(2) for two different orbitals.

  INTEGER :: Wmax
! Max freq for FT

  DOUBLE PRECISION, DIMENSION(1:6) :: stats
! to now how many negative determinant, antisegments,seeme.e.twfs...j

  DOUBLE PRECISION :: swap
! nb of successfull GM

  DOUBLE PRECISION :: signvalue

  INTEGER :: MY_COMM
! 

  INTEGER :: rank
!

  INTEGER :: size
! size of MY_COMM

  DOUBLE PRECISION :: runTime ! time for the run routine
!  

  DOUBLE PRECISION :: beta
!

  DOUBLE PRECISION :: U

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: mu
! levels

  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: hybri_limit
! coeff A such that F=-A/(iwn)

  TYPE(GreenHyboffdiag)                        :: Greens 
! Green's function

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: measN 
! measure of occupations (3or4,flavor) 

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:  ) :: measDE
!  (flavor,flavor) double occupancies
!  (1,1): total energy of correlation.

  DOUBLE PRECISION :: a_Noise
! Noise a exp (-bx) for the  noise

  DOUBLE PRECISION :: b_Noise
! Noise a exp (-bx) for the  noise

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: abNoiseG   !(ab,tau,flavor)
! Noise but for G

  TYPE(Vector)             , DIMENSION(1:2) :: measNoise 
  TYPE(Vector), ALLOCATABLE, DIMENSION(:,:,:) :: measNoiseG       !(tau,flavor,mod) 
! accumulate each value relataed to measurenoise 1 2

!#ifdef CTCtqmcoffdiag_ANALYSIS
!  INTEGER                                     :: order
  DOUBLE PRECISION                            :: inv_dt
! 1/(beta/L)

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:  ) :: measPerturbation 
! opt_order,nflavor

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:  ) :: meas_fullemptylines
! opt_order,nflavor

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: measCorrelation 
! segment,antisegment,nflavor,nflavor

!#endif
!#ifdef CTCtqmcoffdiag_CHECK
  DOUBLE PRECISION :: errorImpurity
! check 

  DOUBLE PRECISION :: errorBath
! for check

!#endif
  TYPE(BathOperatoroffdiag)              :: Bath


  TYPE(ImpurityOperator)          :: Impurity

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: density

END TYPE Ctqmcoffdiag
!!***

!INTERFACE Ctqmcoffdiag_setG0w
!  MODULE PROCEDURE Ctqmcoffdiag_setG0wFile, Ctqmcoffdiag_setG0wTab
!END INTERFACE

CONTAINS
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_init
!! NAME
!!  Ctqmcoffdiag_init
!!
!! FUNCTION
!!  Initialize the type Ctqmcoffdiag
!!  Allocate all the non optional variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
!!  ostream=where to write
!!  istream=where to read the input parameters if so
!!  bFile=logical argument True if input is read from istream
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

SUBROUTINE Ctqmcoffdiag_init(op, ostream, istream, bFile, MY_COMM, iBuffer)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT)                      :: op
  INTEGER  , INTENT(IN   )                      :: ostream
  INTEGER  , INTENT(IN   )                      :: istream
  LOGICAL  , INTENT(IN   )                      :: bFile
  DOUBLE PRECISION, DIMENSION(1:10), OPTIONAL, INTENT(IN) :: iBuffer
  INTEGER  , OPTIONAL, INTENT(IN   )                      :: MY_COMM
!Local variables ------------------------------
#ifdef HAVE_MPI
  INTEGER                                       :: ierr
#endif
  !INTEGER                                       :: iflavor
#ifdef __GFORTRAN__
!  INTEGER                                       :: pid
!  CHARACTER(LEN=5)                              :: Cpid
!
#endif
  DOUBLE PRECISION, DIMENSION(1:10)             :: buffer

  op%ostream = ostream
  op%istream = istream
  
! --- RENICE ---
!#ifdef __GFORTRAN__
!  pid = GetPid()
!  WRITE(Cpid,'(I5)') pid
!  CALL SYSTEM('renice +19 '//TRIM(ADJUSTL(Cpid))//' > /dev/null')
!#endif
!! --- RENICE ---

  IF ( PRESENT(MY_COMM)) THEN
#ifdef HAVE_MPI
    op%have_MPI = .TRUE.
    op%MY_COMM = MY_COMM
    CALL MPI_Comm_rank(op%MY_COMM, op%rank, ierr)
    CALL MPI_Comm_size(op%MY_COMM, op%size, ierr)
#else
    CALL WARN("MPI is not used                                    ")
    op%have_MPI = .FALSE.
    op%MY_COMM = -1
    op%rank = 0
    op%size = 1
#endif
  ELSE
    op%have_MPI = .FALSE.
    op%MY_COMM = -1
    op%rank = 0
    op%size = 1
  END IF

  !IF ( op%rank .EQ. 0 ) THEN
!  WRITE(ostream,'(A20)') 'Job reniced with +19'
    !CALL FLUSH(ostream)
  !END IF

  IF ( bFile .EQV. .TRUE. ) THEN
    IF ( op%rank .EQ. 0 ) THEN

      READ(istream,*) buffer(1) !iseed
      READ(istream,*) buffer(2) !op%sweeps
      READ(istream,*) buffer(3) !op%thermalization
      READ(istream,*) buffer(4) !op%measurements
      READ(istream,*) buffer(5) !op%flavors
      READ(istream,*) buffer(6) !op%samples
      READ(istream,*) buffer(7) !op%beta
      READ(istream,*) buffer(8) !U
      READ(istream,*) buffer(9) !iTech
      !READ(istream,*) buffer(9) !Wmax
!#ifdef CTCtqmcoffdiag_ANALYSIS
      !READ(istream,*) buffer(10) !order
!#endif
    END IF

#ifdef HAVE_MPI
    IF ( op%have_MPI .EQV. .TRUE. ) &
      CALL MPI_Bcast(buffer, 10, MPI_DOUBLE_PRECISION, 0,    &
                   op%MY_COMM, ierr)
#endif
  ELSE IF ( PRESENT(iBuffer) ) THEN
    buffer(1:10) = iBuffer(1:10)
  ELSE
    CALL ERROR("Ctqmcoffdiag_init : No input parameters                    ")
  END IF

  CALL Ctqmcoffdiag_setParameters(op, buffer)

  CALL Ctqmcoffdiag_allocateAll(op)

  CALL GreenHyboffdiag_init(op%Greens,op%samples, op%beta,INT(buffer(5)), &
  iTech=INT(buffer(9)),MY_COMM=op%MY_COMM)


!  op%seg_added    = 0.d0
!  op%anti_added   = 0.d0
!  op%seg_removed  = 0.d0
!  op%anti_removed = 0.d0
!  op%seg_sign     = 0.d0
!  op%anti_sign    = 0.d0
  op%stats(:)     = 0.d0
 ! write(std_out,*) "op%stats",op%stats
  op%signvalue    = 1.d0
!  op%signvaluecurrent    = 0.d0
!  op%signvaluemeas = 0.d0
  op%swap         = 0.d0
  op%runTime      = 0.d0

  CALL Vector_init(op%measNoise(1),op%sweeps/op%modNoise1)
  CALL Vector_init(op%measNoise(2),(op%sweeps/op%modNoise1+1)*CTQMC_SLICE2)
  !CALL Vector_init(op%measNoise(3),101)
  !CALL Vector_init(op%measNoise(4),101)

  op%set  = op%para .AND. op%inF
  op%done = .FALSE.
  op%init = .TRUE.

!#ifdef CTCtqmcoffdiag_CHECK
  op%errorImpurity = 0.d0
  op%errorBath     = 0.d0
!#endif
END SUBROUTINE Ctqmcoffdiag_init
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_setParameters
!! NAME
!!  Ctqmcoffdiag_setParameters
!!
!! FUNCTION
!!  set all parameters and operators
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_setParameters(op,buffer)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT)                         :: op
  DOUBLE PRECISION, DIMENSION(1:10), INTENT(IN   ) :: buffer


  op%thermalization = INT(buffer(3)) !op%thermalization
  CALL Ctqmcoffdiag_setSeed(op,INT(buffer(1)))
  CALL Ctqmcoffdiag_setSweeps(op,buffer(2))

  op%measurements   = INT(buffer(4)) !op%measurements
  op%flavors        = INT(buffer(5))
  op%samples        = INT(buffer(6)) !op%samples
  op%beta           = buffer(7)      !op%beta
  op%U              = buffer(8)      !U
  op%opt_nondiag    = INT(buffer(10))
!  op%mu             = buffer(9)      !op%mu
  !op%Wmax           = INT(buffer(9)) !Freq
!#ifdef CTCtqmcoffdiag_ANALYSIS
!  op%order          = INT(buffer(10)) ! order
  op%inv_dt         = op%samples / op%beta
!#endif

  !CALL ImpurityOperator_init(op%Impurity,op%flavors,op%beta, op%samples)
  CALL ImpurityOperator_init(op%Impurity,op%flavors,op%beta)
  IF ( op%U .GE. 0.d0 ) THEN
    CALL ImpurityOperator_computeU(op%Impurity,op%U,0.d0)
    op%setU = .TRUE.
  END IF
!  op%mu = op%mu + op%Impurity%shift_mu
!sui!write(std_out,*) "op%opt_nondiag",op%opt_nondiag
  CALL BathOperatoroffdiag_init(op%Bath, op%flavors, op%samples, op%beta, INT(buffer(9)), op%opt_nondiag)

  op%para = .TRUE.

END SUBROUTINE Ctqmcoffdiag_setParameters
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_setSweeps
!! NAME
!!  Ctqmcoffdiag_setSweeps
!!
!! FUNCTION
!!  set the number of sweeps
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_setSweeps(op,sweeps)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)         , INTENT(INOUT) :: op
  DOUBLE PRECISION  , INTENT(IN   ) :: sweeps

  op%sweeps = NINT(sweeps / DBLE(op%size))
!  !write(std_out,*) op%sweeps,NINT(sweeps / DBLE(op%size)),ANINT(sweeps/DBLE(op%size))
  IF ( DBLE(op%sweeps) .NE. ANINT(sweeps/DBLE(op%size)) ) &
    CALL ERROR("Ctqmcoffdiag_setSweeps : sweeps is negative or too big     ")
  IF ( op%sweeps .LT. 2*CTQMC_SLICE1 ) THEN  !202
    CALL WARNALL("Ctqmcoffdiag_setSweeps : # sweeps automtically changed     ")
    op%sweeps = 2*CTQMC_SLICE1
!  ELSE IF ( op%sweeps .LT. op%thermalization ) THEN
!    CALL WARNALL("Ctqmcoffdiag_setSweeps : Thermalization > sweeps / cpu -> auto fix")
!    op%sweeps = op%thermalization
  END IF
  IF ( DBLE(NINT(DBLE(op%sweeps)*DBLE(op%size)/DBLE(CTQMC_SLICE1))) .NE.  &
  ANINT(DBLE(op%sweeps)*DBLE(op%size)/DBLE(CTQMC_SLICE1)) ) THEN
    op%modNoise1 = op%sweeps
  ELSE
    op%modNoise1    = MIN(op%sweeps,INT(DBLE(op%sweeps)*DBLE(op%size) / DBLE(CTQMC_SLICE1))) !101
  END IF
  op%modNoise2    = MAX(op%modNoise1 / CTQMC_SLICE2, 1)   ! 100
!  op%modGlobalMove(1) = op%thermalization / 10 + 1
!  op%modGlobalMove(2) = 0

END SUBROUTINE Ctqmcoffdiag_setSweeps
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_setSeed
!! NAME
!!  Ctqmcoffdiag_setSeed
!!
!! FUNCTION
!!  initialize random number generator
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_setSeed(op,iseed)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT)           :: op
  INTEGER  , INTENT(IN   )           :: iseed
!Local variables ------------------------------
  !INTEGER                            :: n
  !INTEGER                            :: i
  !INTEGER, DIMENSION(:), ALLOCATABLE :: seed


  !CALL RANDOM_SEED(size = n)
  !MALLOC(seed,(n))
  !seed =  iseed + (/ (i - 1, i = 1, n) /)

  !CALL RANDOM_SEED(PUT = seed+op%rank)

  !FREE(seed)

  op%seed=INT(iseed+op%rank,8)

END SUBROUTINE Ctqmcoffdiag_setSeed
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_allocateAll
!! NAME
!!  Ctqmcoffdiag_allocateAll
!!
!! FUNCTION
!!  Allocate all non option variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_allocateAll(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER                  :: flavors

  IF ( .NOT. op%para ) &
    CALL ERROR("Ctqmcoffdiag_allocateAll : Ctqmcoffdiag_setParameters never called  ")

  flavors = op%flavors


!  number of electrons
  FREEIF(op%measN)
  MALLOC(op%measN,(1:4,1:flavors))
  op%measN = 0.d0

!  double occupancies 
  FREEIF(op%measDE)
  MALLOC(op%measDE,(1:flavors,1:flavors) )
  op%measDE = 0.d0

  FREEIF(op%mu)
  MALLOC(op%mu,(1:flavors) )
  op%mu = 0.d0
  FREEIF(op%hybri_limit)
  MALLOC(op%hybri_limit,(flavors,flavors) )
  op%hybri_limit = czero
END SUBROUTINE Ctqmcoffdiag_allocateAll
!!***

!#ifdef CTCtqmcoffdiag_ANALYSIS
!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_allocateOpt
!! NAME
!!  Ctqmcoffdiag_allocateOpt
!!
!! FUNCTION
!!  allocate all option variables 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_allocateOpt(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: j
  INTEGER :: k

  IF ( .NOT. op%para ) &
    CALL ERROR("Ctqmcoffdiag_allocateOpt : Ctqmcoffdiag_setParameters never called  ")

  IF ( op%opt_analysis .EQ. 1 ) THEN
    FREEIF(op%measCorrelation)
    MALLOC(op%measCorrelation,(1:op%samples+1,1:3,1:op%flavors))
    op%measCorrelation = 0.d0
  END IF

  IF ( op%opt_order .GT. 0 ) THEN
    FREEIF(op%measPerturbation)
    MALLOC(op%measPerturbation,(1:op%opt_order,1:op%flavors))
    op%measPerturbation = 0.d0
    FREEIF(op%meas_fullemptylines)
    MALLOC(op%meas_fullemptylines,(2,1:op%flavors))
    op%meas_fullemptylines = 0.d0
  END IF

  IF ( op%opt_noise .EQ. 1 ) THEN
    IF ( ALLOCATED(op%measNoiseG) ) THEN
      DO i=1,2
        DO j = 1, op%flavors
          DO k= 1, op%samples+1
            CALL Vector_destroy(op%measNoiseG(k,j,i))
          END DO
        END DO
      END DO
      DT_FREE(op%measNoiseG)
    END IF
    DT_MALLOC(op%measNoiseG,(1:op%samples+1,1:op%flavors,1:2))
    !DO i=1,2
      DO j = 1, op%flavors
        DO k= 1, op%samples+1
          CALL Vector_init(op%measNoiseG(k,j,1),CTQMC_SLICE1)
        END DO
      END DO
      DO j = 1, op%flavors
        DO k= 1, op%samples+1
          CALL Vector_init(op%measNoiseG(k,j,2),CTQMC_SLICE1*CTQMC_SLICE2+1) ! +1 pour etre remplacer ceil
        END DO
      END DO
    !END DO
    FREEIF(op%abNoiseG)
    MALLOC(op%aBNoiseG,(1:2,1:op%samples+1,op%flavors))
    op%abNoiseG = 0.d0
  END IF

  IF (op%opt_spectra .GE. 1 ) THEN
    FREEIF(op%density)
    !MALLOC(op%density,(1:op%thermalization,1:op%flavors))
    i = CEILING(DBLE(op%thermalization+op%sweeps)/DBLE(op%measurements*op%opt_spectra))
    MALLOC(op%density,(1:op%flavors+1,1:i))
    op%density = 0.d0
  END IF
!#endif
END SUBROUTINE Ctqmcoffdiag_allocateOpt
!!***

!SUBROUTINE Ctqmcoffdiag_setG0wFile(op,istream,opt_fk)
!include 'mpif.h'
!#endif
!  TYPE(Ctqmcoffdiag), INTENT(INOUT)                      :: op
!  INTEGER  , INTENT(IN   )                      :: istream
!  INTEGER                         , INTENT(IN ) :: opt_fk
!!  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
!  COMPLEX*16      , DIMENSION(:,:), ALLOCATABLE :: Gomega
!  INTEGER                                       :: flavors
!  INTEGER                                       :: iflavor
!  INTEGER                                       :: iomega
!  TYPE(GreenHyboffdiag) :: temp
!#ifdef HAVE_MPI
!  INTEGER                                       :: ierr
!#endif
!
!  IF ( .NOT. op%para ) &
!    CALL ERROR("Ctqmcoffdiag_setG0wFile : Ctqmcoffdiag_setParameters never called   ") 
!
!  flavors = op%flavors
!
!!  MALLOC(F,(1:op%samples+1,1:flavors))
!  MALLOC(Gomega,(1:op%Wmax,1:flavors))
!  IF ( op%have_MPI .EQV. .TRUE. ) THEN
!    CALL GreenHyboffdiag_init(temp, op%samples, op%beta, op%MY_COMM)
!  ELSE
!    CALL GreenHyboffdiag_init(temp, op%samples, op%beta)
!  END IF
!
!  DO iflavor = 1, flavors
!    IF ( op%rank .EQ. 0 ) THEN
!      DO iomega=1, op%Wmax
!        READ(istream,*) Gomega(iomega,iflavor)
!      END DO
!    END IF
!  END DO
!
!#ifdef HAVE_MPI
!  CALL MPI_Bcast(Gomega, op%Wmax*flavors, MPI_DOUBLE_COMPLEX, 0,    &
!                 op%MY_COMM, ierr)
!#endif
!
!!  CALL GreenHyboffdiag_backFourier(temp,Gomega(:,1))
!!  WRITE(10,*) Gomega(:,1)
!!  CALL GreenHyboffdiag_print(temp)
!!  CALL GreenHyboffdiag_forFourier (temp,Gomega(:,1))
!!  WRITE(11,*) Gomega(:,1)
!  CALL Ctqmcoffdiag_setG0wTab(op,Gomega,opt_fk)
!!  CALL Ctqmcoffdiag_computeF(op,Gomega, op%Wmax, F) ! mu is changed
!  FREE(Gomega)
!!  CALL BathOperatoroffdiag_setF(op%Bath, F)
!!  CALL BathOperatoroffdiag_printF(op%Bath)
!!  FREE(F)
!
!!  stop
!!  op%inF = .TRUE.
!!  op%set = .TRUE. 
!
!END SUBROUTINE Ctqmcoffdiag_setG0wFile
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_setG0wTab
!! NAME
!!  Ctqmcoffdiag_setG0wTab
!!
!! FUNCTION
!!  Set Gow from input array
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_setG0wTab(op,Gomega,opt_fk)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT)                      :: op
  COMPLEX(KIND=8), DIMENSION(:,:,:), INTENT(IN ) :: Gomega
  INTEGER                         , INTENT(IN ) :: opt_fk
!Local variable -------------------------------
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: F

  IF ( .NOT. op%para ) &
    CALL ERROR("Ctqmcoffdiag_setG0wTab : Ctqmcoffdiag_setParameters never called    ") 

  MALLOC(F,(1:op%samples+1,1:op%flavors,1:op%flavors))
  CALL Ctqmcoffdiag_computeF(op,Gomega, F, opt_fk)  ! mu is changed
  CALL BathOperatoroffdiag_setF(op%Bath, F)
 ! CALL BathOperatoroffdiag_printF(op%Bath,333)
  FREE(F)

  op%inF = .TRUE.
  op%set = .TRUE. 

END SUBROUTINE Ctqmcoffdiag_setG0wTab
!!***

!SUBROUTINE Ctqmcoffdiag_setFwK(op,Gomega)
!  COMPLEX*16      , DIMENSION(:,:), INTENT(IN ) :: Gomega
!  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: F
!
!  IF ( .NOT. op%para ) &
!    CALL ERROR("Ctqmcoffdiag_setG0wTab : Ctqmcoffdiag_setParameters never called    ") 
!
!  MALLOC(F,(1:op%samples+1,1:op%flavors))
!  CALL Ctqmcoffdiag_computeFK(op,Gomega, op%Wmax, F)  ! mu is changed
!  CALL BathOperatoroffdiag_setF(op%Bath, F)
!  CALL BathOperatoroffdiag_printF(op%Bath)
!  FREE(F)
!
!  op%inF = .TRUE.
!  op%set = .TRUE. 
!
!END SUBROUTINE Ctqmcoffdiag_setFwK
!!***

!SUBROUTINE Ctqmcoffdiag_setBand(op, mu, U)
!  DOUBLE PRECISION, INTENT(IN   ) :: mu
!  DOUBLE PRECISION, INTENT(IN   ) :: U
!
!  IF ( .NOT. op%para ) &
!    CALL ERROR("Ctqmcoffdiag_setBand : Ctqmcoffdiag_setParameters never called      ")
!
!  CALL ImpurityOperator_setU(op%Impurity, U, 0.d0)
!  !op%mu = mu + op%Impurity%shift_mu
!END SUBROUTINE Ctqmcoffdiag_setBand
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_setU
!! NAME
!!  Ctqmcoffdiag_setU
!!
!! FUNCTION
!!  set the interaction matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_setU(op,matU)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) ::op
!Local variables ------------------------------
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: matU

  IF ( SIZE(matU) .NE. op%flavors*op%flavors ) &
    CALL ERROR("Ctqmcoffdiag_setU : Wrong interaction matrix (size)        ")

  CALL ImpurityOperator_setUmat(op%Impurity, matU)
  op%setU = .TRUE.
END SUBROUTINE Ctqmcoffdiag_setU
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_clear
!! NAME
!!  Ctqmcoffdiag_clear
!!
!! FUNCTION
!!  clear a ctqmc run
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_clear(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: j
  INTEGER :: k

  op%measN(1,:) = 0.d0
  op%measN(2,:) = 0.d0
  !Do not set measN(3,:) to 0 to avoid erasing N between therm and ctqmc
  op%measN(4,:) = 0.d0
  op%measDE = 0.d0
!  op%seg_added    = 0.d0
!  op%anti_added   = 0.d0
!  op%seg_removed  = 0.d0
!  op%anti_removed = 0.d0
!  op%seg_sign     = 0.d0
!  op%anti_sign    = 0.d0
  op%stats(:)     = 0.d0
!  op%signvaluecurrent    = 0.d0
!  op%signvaluemeas    = 0.d0
  op%swap         = 0.d0
  op%runTime      = 0.d0
  op%modGlobalMove(2) = 0 
  CALL Vector_clear(op%measNoise(1))
  CALL Vector_clear(op%measNoise(2))
!#ifdef CTCtqmcoffdiag_CHECK
  op%errorImpurity = 0.d0
  op%errorBath     = 0.d0
!#endif
  CALL GreenHyboffdiag_clear(op%Greens)
!#ifdef CTCtqmcoffdiag_ANALYSIS
  IF ( op%opt_analysis .EQ. 1 .AND. ALLOCATED(op%measCorrelation) ) &    
    op%measCorrelation = 0.d0 
  IF ( op%opt_order .GT. 0 .AND. ALLOCATED(op%measPerturbation) ) &
    op%measPerturbation = 0.d0
  IF ( op%opt_order .GT. 0 .AND. ALLOCATED(op%meas_fullemptylines) ) &
    op%meas_fullemptylines = 0.d0
  IF ( op%opt_noise .EQ. 1 .AND. ALLOCATED(op%measNoiseG) ) THEN
    DO i=1,2
      DO j = 1, op%flavors
        DO k= 1, op%samples+1
          CALL Vector_clear(op%measNoiseG(k,j,i))
        END DO
      END DO
    END DO
    !DO j = 1, op%flavors
    !  CALL GreenHyboffdiag_clear(op%Greens(j))
    !END DO
  END IF
!#endif
END SUBROUTINE Ctqmcoffdiag_clear
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_reset
!! NAME
!!  Ctqmcoffdiag_reset
!!
!! FUNCTION
!!  reset a ctqmc simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_reset(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  !INTEGER                  :: iflavor
  DOUBLE PRECISION         :: sweeps

  CALL GreenHyboffdiag_reset(op%Greens)
  CALL Ctqmcoffdiag_clear(op)
  CALL ImpurityOperator_reset(op%Impurity)
  CALL BathOperatoroffdiag_reset    (op%Bath)
  op%measN(3,:) = 0.d0
  !complete restart -> measN=0
  op%done = .FALSE.
  op%set  = .FALSE.
  op%inF  = .FALSE.
  op%opt_movie = 0
  op%opt_analysis = 0
  op%opt_order = 0
  op%opt_check = 0
  op%opt_noise = 0
  op%opt_spectra = 0
  op%opt_levels = 0
  sweeps = DBLE(op%sweeps)*DBLE(op%size)
  CALL Ctqmcoffdiag_setSweeps(op, sweeps)
!#ifdef HAVE_MPI
!  CALL MPI_BARRIER(op%MY_COMM,iflavor)
!  IF ( op%rank .EQ. 0 ) &
!#endif
!  WRITE(op%ostream,'(A9)') "QMC reset"
!  CALL FLUSH(op%ostream)
END SUBROUTINE Ctqmcoffdiag_reset
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_setMu
!! NAME
!!  Ctqmcoffdiag_setMu
!!
!! FUNCTION
!!  impose energy levels
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_setMu(op, levels)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)                     , INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN   ) :: levels

  IF ( op%flavors .NE. SIZE(levels,1) ) &
    CALL WARNALL("Ctqmcoffdiag_setMu : Taking energy levels from weiss G(iw)")

  op%mu(:)=-levels(:)  ! levels = \epsilon_j - \mu
  !op%mu =\tilde{\mu} = \mu -\epsilon_j
  op%opt_levels = 1
END SUBROUTINE Ctqmcoffdiag_setMu
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_sethybri_limit
!! NAME
!!  Ctqmcoffdiag_sethybri_limit
!!
!! FUNCTION
!!  use coefficient A such that F=-A/(iwn) given by DMFT code.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  hybri_limit(nflavor,nflavor)=contains the limit for each couple of flavors
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!  op(Ctqmcoffdiag_type) = is the ctqmc main variable
!!           op&limit is now filled
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE Ctqmcoffdiag_sethybri_limit(op, hybri_limit)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)                     , INTENT(INOUT) :: op
  COMPLEX(KIND=8) , DIMENSION(:,:),  INTENT(IN ) :: hybri_limit

  IF ( op%flavors .NE. SIZE(hybri_limit,1) ) &
    CALL ERROR("Error in sethybri_limit")

  op%hybri_limit(:,:)=hybri_limit(:,:)  
  op%opt_hybri_limit = 1
END SUBROUTINE Ctqmcoffdiag_sethybri_limit
!!***
!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_computeF
!! NAME
!!  Ctqmcoffdiag_computeF
!!
!! FUNCTION
!!  Compute the hybridization function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_computeF(op, Gomega, F, opt_fk)

 use m_hide_lapack,  only : xginv
!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)                       , INTENT(INOUT) :: op
  COMPLEX(KIND=8), DIMENSION(:,:,:), INTENT(IN   ) :: Gomega
  !INTEGER                         , INTENT(IN   ) :: Wmax
  DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: F
  INTEGER                         , INTENT(IN   ) :: opt_fk
!Local variables ------------------------------
  INTEGER                                         :: flavors
  INTEGER                                         :: samples
  INTEGER                                         :: iflavor,ifl
  INTEGER                                         :: iflavor2
  INTEGER                                         :: iomega
  INTEGER                                         :: itau
  DOUBLE PRECISION                                :: pi_invBeta
  DOUBLE PRECISION                                :: K
  !DOUBLE PRECISION                                :: re
  !DOUBLE PRECISION                                :: im
  !DOUBLE PRECISION                                :: det
  COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE   :: F_omega
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: F_omega_inv
  COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE   :: Gomega_tmp
  TYPE(GreenHyboffdiag)                                     :: F_tmp
  !character(len=4) :: tag_proc
  !character(len=30) :: tmpfil
  !INTEGER :: unitnb

  ABI_UNUSED((/opt_fk/))

  flavors    = op%flavors

  samples    = op%samples
  pi_invBeta = ACOS(-1.d0) / op%beta
  op%Wmax=SIZE(Gomega,1)
!sui!write(std_out,*) "op%Wmax",op%Wmax
  !=================================
  ! --- Initialize F_tmp 
  !=================================
  IF ( op%have_MPI .EQV. .TRUE. ) THEN
    CALL GreenHyboffdiag_init(F_tmp,samples,op%beta,flavors,MY_COMM=op%MY_COMM)
  ELSE
    CALL GreenHyboffdiag_init(F_tmp,samples,op%beta,flavors)
  END IF
!  K = op%mu

  !=================================
  ! --- Allocate F_omega
  !=================================
  MALLOC(F_omega,(1:op%Wmax,1:flavors,1:flavors))
  MALLOC(F_omega_inv,(1:flavors,1:flavors))
  MALLOC(Gomega_tmp,(1:op%Wmax,1:flavors,1:flavors))
  !op%hybri_limit(2,2)=op%hybri_limit(1,1)
  !op%mu(1)=op%mu(1)/10
  !op%mu(2)=op%mu(1)
  DO iomega=1,op%Wmax
    do iflavor=1,flavors
      do iflavor2=1,flavors
       ! Gomega_tmp(iomega,iflavor,iflavor2)=op%hybri_limit(iflavor,iflavor2)/(cmplx(0.d0,(2.d0*DBLE(iomega)-1.d0) * pi_invBeta))/3.d0
      enddo
    enddo
  END DO
  Gomega_tmp=Gomega

  !IF ( op%rank .EQ. 0 ) &
    !OPEN(UNIT=9876,FILE="K.dat",POSITION="APPEND")
  
  !=============================================================================================
  ! --- Compute Bath Green's function from Hybridization function in imaginary time
  !=============================================================================================
  !IF ( opt_fk .EQ. 0 ) THEN
   IF ( op%rank .EQ. 0 ) THEN
   !  DO iflavor = 1, flavors
   !    DO iflavor2 = 1, flavors
   !        write(330,*) "#",iflavor,iflavor2
   !        write(331,*) "#",iflavor,iflavor2
   !      do  iomega=1,op%Wmax
   !        write(330,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(Gomega_tmp(iomega,iflavor,iflavor2))
   !        write(331,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(Gomega_tmp(iomega,iflavor,iflavor2))
   !      enddo
   !        write(330,*) 
   !        write(331,*) 
   !    END DO
   !  END DO
   ENDIF
     DO iomega=1,op%Wmax
     !  be careful...here 
     ! Gomega in input is Fomega and 
     ! F_omega is   Gomega.
     ! COMPUTE G0 FROM F
       do iflavor=1,flavors
         do iflavor2=1,flavors
           if (iflavor==iflavor2) then
             F_omega_inv(iflavor,iflavor2)= (cmplx(0.d0,(2.d0*DBLE(iomega)-1.d0) * pi_invBeta,kind=8) &
&             + op%mu(iflavor)- Gomega_tmp(iomega,iflavor,iflavor2))
           else
             F_omega_inv(iflavor,iflavor2)= (- Gomega_tmp(iomega,iflavor,iflavor2))
           endif
         enddo
       enddo
  !   END DO
  ! IF ( op%rank .EQ. 0 ) THEN
  !   DO iflavor = 1, flavors
  !     DO iflavor2 = 1, flavors
  !         write(334,*) "#",iflavor,iflavor2
  !         write(335,*) "#",iflavor,iflavor2
  !       do  iomega=1,op%Wmax
  !         write(334,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_omega(iomega,iflavor,iflavor2))
  !         write(335,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(F_omega(iomega,iflavor,iflavor2))
  !       enddo
  !         write(334,*) 
  !         write(335,*) 
  !     END DO
  !   END DO
  ! ENDIF

  !   DO iomega=1,op%Wmax
       call xginv(F_omega_inv,flavors)
       do iflavor=1,flavors
         do iflavor2=1,flavors
           F_omega(iomega,iflavor,iflavor2) = F_omega_inv(iflavor,iflavor2)
         enddo
       enddo
     END DO

   !IF ( op%rank .EQ. 0 ) THEN
   !  DO iflavor = 1, flavors
   !    DO iflavor2 = 1, flavors
   !        write(332,*) "#",iflavor,iflavor2
   !        write(333,*) "#",iflavor,iflavor2
   !      do  iomega=1,op%Wmax
   !        write(332,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_omega(iomega,iflavor,iflavor2))
   !        write(333,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(F_omega(iomega,iflavor,iflavor2))
   !      enddo
   !        write(332,*) 
   !        write(333,*) 
   !    END DO
   !  END DO
   !ENDIF
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  for test: Fourier of G0(iwn)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !sui!write(std_out,*) "opt_fk=0"
     CALL GreenHyboffdiag_setOperW(F_tmp,F_omega)
  ! IF ( op%rank .EQ. 0 ) THEN
  !    DO iflavor = 1, flavors
  !      DO iflavor2 = 1, flavors
  !          write(336,*) "#",iflavor,iflavor2
  !          write(337,*) "#",iflavor,iflavor2
  !        do  iomega=1,op%Wmax
  !          write(336,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_tmp%oper_w(iomega,iflavor,iflavor2))
  !          write(337,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(F_tmp%oper_w(iomega,iflavor,iflavor2))
  !        enddo
  !          write(336,*) 
  !          write(337,*) 
  !      END DO
  !    END DO
  !  ENDIF
     CALL GreenHyboffdiag_backFourier(F_tmp,func="green")
     ! --- Put the result in F
     DO iflavor = 1, flavors
       DO iflavor2 = 1, flavors
         DO itau=1,samples+1
         F(itau,iflavor,iflavor2) = F_tmp%oper(itau,iflavor,iflavor2)
         END DO
       END DO
     END DO
     !IF ( op%rank .EQ. 0 ) THEN
     !  DO iflavor = 1, flavors
     !    DO iflavor2 = 1, flavors
     !        write(346,*) "#",iflavor,iflavor2
     !      do  itau=1,op%samples+1
     !        write(346,*) (itau-1)*op%beta/(op%samples),real(F(itau,iflavor,iflavor2))
     !      enddo
     !        write(346,*) 
     !    END DO
     !  END DO
     !ENDIF
     DO iflavor = 1, flavors
       DO iflavor2 = 1, flavors
         DO itau=1,samples+1
!         This symetrization is general and valid even with SOC
!         Without SOC, it leads to zero.
         F(itau,iflavor,iflavor2) = (F_tmp%oper(itau,iflavor,iflavor2)+F_tmp%oper(itau,iflavor2,iflavor))/2.d0
         END DO
       END DO
     END DO
     open (unit=4367,file='G0tau_fromF',status='unknown',form='formatted')
     rewind(4367)
     IF ( op%rank .EQ. 0 ) THEN
       DO iflavor = 1, flavors
         DO iflavor2 = 1, flavors
             write(4367,*) "#",iflavor,iflavor2
           do  itau=1,op%samples+1
             write(4367,*) (itau-1)*op%beta/(op%samples),real(F(itau,iflavor,iflavor2))
           enddo
             write(4367,*) 
         END DO
       !sui!write(std_out,'(5x,14(2f9.5,2x))') (F(op%samples+1,iflavor,iflavor2),iflavor2=1,flavors)
       END DO
     ENDIF
     !call flush(436)
     !call flush(437)
     close(4367)
     !call flush(6)
     
     call xmpi_barrier(op%MY_COMM)
     !CALL ERROR("END OF CALCULATION")
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  END OF TEST
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !DO iomega=1,op%Wmax
     !  call xginv(F_omega(iomega,:,:),flavors)
     !END DO
    !F_omega = CMPLX(-1.d0,0,8)/Gomega_tmp
  !ELSE
  !=============================================================================================
  ! --- Restore Hybridization in F_omega
  !=============================================================================================

!   Restore Hybridization in F_omega for the following operations
  F_omega = Gomega_tmp
  !END IF

  !==================================================================
  ! --- Full double loop on flavors to compute F (remove levels)
  !==================================================================
  DO iflavor = 1, flavors
    DO iflavor2 = 1, flavors

  ! --- Compute or use the levels for the diagonal hybridization (else K=0)
      IF(iflavor==iflavor2) THEN
        IF ( op%opt_levels .EQ. 1 ) THEN
          K = op%mu(iflavor)
        ELSE
          K = -REAL(F_omega(op%Wmax, iflavor,iflavor))
!        op%mu = K
          op%mu(iflavor) = K 
        END IF
      ELSE
        K=0.d0
      ENDIF
      !IF ( op%rank .EQ. 0 ) &
      !WRITE(9876,'(I4,2E22.14)') iflavor, K, REAL(-F_omega(op%Wmax, iflavor))
     ! IF(op%rank .EQ.0) &
     ! WRITE(op%ostream,*) "CTQMC K, op%mu = ",K,op%mu(iflavor)
      !WRITE(op%ostream,*) "CTQMC beta     = ",op%beta

  ! --- Compute F (by removing the levels) if opt_fk==0
    !  IF ( opt_fk .EQ. 0 ) THEN
    !   ! DO iomega = 1, op%Wmax
    !   !   re = REAL(F_omega(iomega,iflavor,iflavor2))
    !   !   im = AIMAG(F_omega(iomega,iflavor,iflavor2))
    !   !   if (iflavor==iflavor2) then
    !   !     F_omega(iomega,iflavor,iflavor) = CMPLX(re + K, im + (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, 8)
    !   !   else
    !   !     F_omega(iomega,iflavor,iflavor2) = CMPLX(re , im  , 8)
    !   !   endif
    !   !   !if(iflavor==1.and.op%rank==0) then
    !   !     !write(224,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(F_omega(iomega,iflavor)),imag(F_omega(iomega,iflavor))
    !   !     !write(225,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(Gomega_tmp(iomega, iflavor)),imag(Gomega_tmp(iomega, iflavor))
    !   !   !end if 
    !   ! END DO
    !  ELSE
    !    DO iomega = 1, op%Wmax
    !      !F_omega(iomega,iflavor,iflavor2) = F_omega(iomega,iflavor,iflavor2) + CMPLX(K, 0.d0, 8)


    !      !if(iflavor==1.and.op%rank==0) then
    !        !write(224,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(F_omega(iomega,iflavor)),imag(F_omega(iomega,iflavor))
    !        !write(225,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(Gomega_tmp(iomega, iflavor)),imag(Gomega_tmp(iomega, iflavor))
    !      !end if 
    !    END DO
    !  END IF
  ! --- compute residual K (?)
      K = REAL(CMPLX(0,(2.d0*DBLE(op%Wmax)-1.d0)*pi_invBeta,8)*F_omega(op%Wmax,iflavor,iflavor2))
      CALL GreenHyboffdiag_setMuD1(op%Greens,iflavor,iflavor2,op%mu(iflavor),K)
    END DO
  END DO

  do  iomega=1,op%Wmax
   ! write(336,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_omega(iomega,1,1)),imag(F_omega(iomega,1,1))
  enddo

  ! --- Creates F_tmp%oper_w
  CALL GreenHyboffdiag_setOperW(F_tmp,F_omega)
 ! do  iflavor=1, flavors ; do  iflavor2=1, flavors ; write(337,*) "#",iflavor,iflavor2 ; do  iomega=1,op%Wmax
 !   write(337,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_tmp%oper_w(iomega,iflavor,iflavor2)),&
 !&   imag(F_tmp%oper_w(iomega,iflavor,iflavor2))
 ! enddo ; write(337,*) ; enddo ; enddo
!  IF ( op%rank .EQ. 0 ) THEN
!    DO iflavor = 1, flavors
!      DO iflavor2 = 1, flavors
!        write(336,*) "#",iflavor,iflavor2
!        write(337,*) "#",iflavor,iflavor2
!        do  iomega=1,op%Wmax
!          write(336,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_tmp%oper_w(iomega,iflavor,iflavor2))
!          write(337,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(F_tmp%oper_w(iomega,iflavor,iflavor2))
!        enddo
!        write(336,*) 
!        write(337,*) 
!        write(136,*) "#",iflavor,iflavor2
!        write(137,*) "#",iflavor,iflavor2
!        do  iomega=1,op%Wmax
!        write(136,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_tmp%oper_w(iomega,iflavor,iflavor2)-op%hybri_limit(iflavor,iflavor2)/(cmplx(0.d0,(2.d0*DBLE(iomega)-1.d0) * pi_invBeta,kind=8)))
!        write(137,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(F_tmp%oper_w(iomega,iflavor,iflavor2)-op%hybri_limit(iflavor,iflavor2)/(cmplx(0.d0,(2.d0*DBLE(iomega)-1.d0) * pi_invBeta,kind=8)))
!        enddo
!        write(136,*) 
!        write(137,*) 
!        write(836,*) "#",iflavor,iflavor2
!        write(837,*) "#",iflavor,iflavor2
!        do  iomega=1,op%Wmax
!        write(836,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%hybri_limit(iflavor,iflavor2)/(cmplx(0.d0,(2.d0*DBLE(iomega)-1.d0) * pi_invBeta)))
!        write(837,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(op%hybri_limit(iflavor,iflavor2)/(cmplx(0.d0,(2.d0*DBLE(iomega)-1.d0) * pi_invBeta)))
!        enddo
!        write(836,*) 
!        write(837,*) 
!      END DO
!    END DO
!  ENDIF
  !CALL GreenHyboffdiag_backFourier(F_tmp,F_omega(:,iflavor))
   ! DO iflavor = 1, flavors
   !   DO iflavor2 = 1, flavors
   !   unitnb=80000+F_tmp%rank
   !   call int2char4(F_tmp%rank,tag_proc)
   !   tmpfil = 'oper_wavantFOURIER'//tag_proc
   !   open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
   !   write(unitnb,*) "#",iflavor,iflavor2
   !   ! C_omega et oper_w differents Domega identique. Est ce du a des
   !   ! diago differentes   pour chaque procs dans qmc_prep_ctqmc
   !   do  iomega=1,F_tmp%Wmax
   !   write(unitnb,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(F_tmp%oper_w(iomega,iflavor,iflavor2))
   !   enddo
   !   write(unitnb,*) 
   !   END DO
   ! END DO

  ! --- For all iflavor and iflavor2, do the Fourier transformation to
  ! --- have (F(\tau))
  !CALL GreenHyboffdiag_backFourier(F_tmp,hybri_limit=op%hybri_limit,opt_hybri_limit=op%opt_hybri_limit)
   write(std_out,*) "WARNING opt_hybri_limit==0"
  CALL GreenHyboffdiag_backFourier(F_tmp,hybri_limit=op%hybri_limit,opt_hybri_limit=0)
!  CALL GreenHyboffdiag_backFourier(F_tmp,hybri_limit=op%hybri_limit,opt_hybri_limit=1)
 ! CALL GreenHyboffdiag_backFourier(F_tmp)

  ! --- Put the result in F
  DO iflavor = 1, flavors
    DO iflavor2 = 1, flavors
      DO itau=1,samples+1
      F(itau,iflavor,iflavor2) = -F_tmp%oper(samples+2-itau,iflavor,iflavor2)
      END DO
    END DO
  END DO
! IF ( op%rank .EQ. 0 ) THEN
!   ifl=0
!   DO iflavor = 1, flavors
!     DO iflavor2 = 1, flavors
!       ifl=ifl+1
!       write(346,*) "#",iflavor,iflavor2,ifl
!       do  itau=1,op%samples+1
!         write(346,*) itau,real(F(itau,iflavor,iflavor2))
!       enddo
!       write(346,*) 
!     END DO
!   END DO
! ENDIF
! close(346)
  DO iflavor = 1, flavors
    DO iflavor2 = 1, flavors
      DO itau=1,samples+1
!      This symetrization is general and valid even with SOC
!      Without SOC, it leads to zero.
      F(itau,iflavor,iflavor2) = -(F_tmp%oper(samples+2-itau,iflavor,iflavor2)+F_tmp%oper(samples+2-itau,iflavor2,iflavor))/2.d0
      END DO
    END DO
  END DO
  !DO iflavor = 1, flavors
  !  DO iflavor2 = 1, flavors
  !    DO itau=1,samples+1
  !    F(itau,iflavor,iflavor2) = F(samples/2,iflavor,iflavor2)
  !    END DO
  !  END DO
  !END DO

 !  SOME TRY TO ADJUST F
  !DO iflavor = 1, flavors
  !  DO iflavor2 = 1, flavors
  !    do  itau=1,op%samples+1
  !    !if(iflavor/=iflavor2) F(itau,iflavor,iflavor2)=F((op%samples+1)/2,iflavor,iflavor2)
  !    !if(iflavor==iflavor2) F(itau,iflavor,iflavor2)=F((op%samples+1)/2,iflavor,iflavor2)
  !    enddo
  !  END DO
  !END DO

  IF ( op%rank .EQ. 0 ) THEN
    open (unit=436,file='Hybridization.dat',status='unknown',form='formatted')
    rewind(436)
    ifl=0
    DO iflavor = 1, flavors
      DO iflavor2 = 1, flavors
        ifl=ifl+1
          write(436,*) "#",iflavor,iflavor2,ifl,op%hybri_limit(iflavor,iflavor2)
        do  itau=1,op%samples+1
          write(436,*) itau,real(F(itau,iflavor,iflavor2))
        enddo
          write(436,*) 
      END DO
    END DO
    close(436)
  ENDIF
  FREE(Gomega_tmp)
  FREE(F_omega)
  FREE(F_omega_inv)
  CALL GreenHyboffdiag_destroy(F_tmp)


END SUBROUTINE Ctqmcoffdiag_computeF
!!***

!SUBROUTINE Ctqmcoffdiag_computeFK(op, Gomega, Wmax, F)
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
!  TYPE(GreenHyboffdiag)                                     :: F_tmp
!
!  flavors    = op%flavors
!
!  samples    = op%samples
!  pi_invBeta = ACOS(-1.d0) / op%beta
!
!  IF ( op%have_MPI .EQV. .TRUE. ) THEN
!    CALL GreenHyboffdiag_init(F_tmp,samples,op%beta,op%MY_COMM)
!  ELSE
!    CALL GreenHyboffdiag_init(F_tmp,samples,op%beta)
!  END IF
!!  K = op%mu
!
!  MALLOC(F_omega,(1:Wmax,1:flavors))
!
!  DO iflavor = 1, flavors
!    K = REAL(Gomega(Wmax, iflavor))
!    WRITE(op%ostream,*) "CTQMC K, op%mu = ",K,op%mu
!    WRITE(op%ostream,*) "CTQMC beta     = ",op%beta
!    op%mu(iflavor) = K 
!    DO iomega = 1, Wmax
!      F_omega(iomega,iflavor) = Gomega(iomega,iflavor) &
!                  - CMPLX(K, 0.d0, 8)
!      !if(iflavor==1.and.op%rank==0) then
!        !write(224,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(F_omega(iomega,iflavor)),imag(F_omega(iomega,iflavor))
!        !write(225,*) (2.d0*DBLE(iomega)-1.d0) * pi_invBeta, real(Gomega(iomega, iflavor)),imag(Gomega(iomega, iflavor))
!      !end if 
!    END DO
!    CALL GreenHyboffdiag_backFourier(F_tmp,F_omega(:,iflavor))
!    F(1:samples+1,iflavor) = (/ (-F_tmp%oper(samples+1-itau),itau=0,samples) /)
!  END DO
!  FREE(F_omega)
!  CALL GreenHyboffdiag_destroy(F_tmp)
!END SUBROUTINE Ctqmcoffdiag_computeFK
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_run
!! NAME
!!  Ctqmcoffdiag_run
!!
!! FUNCTION
!!  set all options and run a simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_run(op,opt_order,opt_movie,opt_analysis,opt_check,opt_noise,opt_spectra,opt_gMove)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT)           :: op
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
!#ifdef CTCtqmcoffdiag_MOVIE
  INTEGER                            :: ilatex
  CHARACTER(LEN=4)                   :: Cchar
!#endif
  DOUBLE PRECISION                   :: estimatedTime

  IF ( .NOT. op%set  ) &
    CALL ERROR("Ctqmcoffdiag_run : QMC not set up                          ")
  IF ( .NOT. op%setU ) &
    CALL ERROR("Ctqmcoffdiag_run : QMC does not have a U matrix            ")


! OPTIONS of the run
  IF ( PRESENT( opt_check ) ) THEN
    op%opt_check = opt_check
    CALL ImpurityOperator_doCheck(op%Impurity,opt_check)
    CALL BathOperatoroffdiag_doCheck(op%Bath,opt_check)
  END IF
  IF ( PRESENT( opt_movie ) ) &
    op%opt_movie = opt_movie
  IF ( PRESENT( opt_analysis ) ) &
    op%opt_analysis = opt_analysis
  IF ( PRESENT ( opt_order ) ) &
    op%opt_order = opt_order 
  IF ( PRESENT ( opt_noise ) ) THEN
    op%opt_noise = opt_noise 
  END IF
  IF ( PRESENT ( opt_spectra ) ) &
    op%opt_spectra = opt_spectra

  op%modGlobalMove(1) = max(op%sweeps,op%thermalization)+1 ! No Global Move
!!sui!write(std_out,*) "op%sweeps",op%thermalization,op%sweeps,opt_gMove
  op%modGlobalMove(2) = 0
  IF ( PRESENT ( opt_gMove ) ) THEN
    IF ( opt_gMove .LE. 0 .OR. opt_gMove .GT. op%sweeps ) THEN
     ! op%modGlobalMove(1) = op%sweeps+1
      op%modGlobalMove(1) = max(op%sweeps,op%thermalization)+1 ! No Global Move
      !write(std_out,*) "op%sweeps",op%sweeps, op%modGlobalMove(1)
      CALL WARNALL("Ctqmcoffdiag_run : global moves option is <= 0 or > sweeps/cpu -> No global Moves")
    ELSE 
      op%modGlobalMove(1) = opt_gMove 
    END IF
  END IF
!sui!write(std_out,*) "op%sweeps",op%thermalization,op%sweeps

  CALL Ctqmcoffdiag_allocateOpt(op)
  
!#ifdef CTCtqmcoffdiag_MOVIE  
  ilatex = 0
  IF ( op%opt_movie .EQ. 1 ) THEN
    Cchar ="0000"
    WRITE(Cchar,'(I4)') op%rank 
    ilatex = 87+op%rank
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

  IF ( op%rank .EQ. 0 ) THEN
    WRITE(op%ostream,'(A29)') "Starting QMC (Thermalization)"
  END IF
  
  !=================================
  ! STARTING THERMALIZATION 
  !=================================
  !write(std_out,*) "sweeps before thermalization",op%sweeps
  !write(std_out,*) "op%stats",op%stats
  CALL Ctqmcoffdiag_loop(op,op%thermalization,ilatex)
  !=================================
  ! ENDING   THERMALIZATION 
  !=================================

  estimatedTime = op%runTime
#ifdef HAVE_MPI
  CALL MPI_REDUCE(op%runTime, estimatedTime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             0, op%MY_COMM, ierr)
#endif

  IF ( op%rank .EQ. 0 ) THEN
    WRITE(op%ostream,'(A26,I6,A11)') "Thermalization done in    ", CEILING(estimatedTime), "    seconds"
    WRITE(op%ostream,'(A25,I7,A15,I5,A5)') "The QMC should run in    ", &
           CEILING(estimatedTime*DBLE(op%sweeps)/DBLE(op%thermalization)),&
                        "    seconds on ", op%size, " CPUs"
  END IF

  !=================================
  ! CLEANING CTQMC          
  !=================================
  CALL Ctqmcoffdiag_clear(op)

  !=================================
  ! STARTING CTQMC          
  !=================================
  !write(std_out,*) "sweeps before loop",op%sweeps
  !write(std_out,*) "op%stats",op%stats
  CALL Ctqmcoffdiag_loop(op,op%sweeps,ilatex)
  !=================================
  ! ENDING   CTQMC          
  !=================================

  IF ( op%opt_movie .EQ. 1 ) THEN
    WRITE(ilatex,*) ""
    WRITE(ilatex,'(A14)') "\end{document}"
    CLOSE(ilatex)
  END IF

  op%done     = .TRUE.
!sui!write(std_out,*) "op%stats en of ctqmc_run",op%stats

END SUBROUTINE Ctqmcoffdiag_run
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_loop
!! NAME
!!  Ctqmcoffdiag_loop
!!
!! FUNCTION
!!  Definition the main loop of the CT-QMC
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_loop(op,itotal,ilatex)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT)         :: op
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
  INTEGER                            :: ipercent !,ii
  INTEGER                            :: iflavor,ifl1,iflavor_d
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

  flavors        = op%flavors
  measurements   = op%measurements
  modNoise1      = op%modNoise1
  modNoise2      = op%modNoise2
  modGlobalMove  = op%modGlobalMove(1)
  sp1            = op%samples+1

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

  endDensity = SIZE(op%density,2)

  IF ( op%opt_noise .GT. 0 ) THEN
    FREEIF(gtmp_new)
    MALLOC(gtmp_new,(1:sp1,1:flavors))
    FREEIF(gtmp_old1)
    MALLOC(gtmp_old1,(1:sp1,1:flavors))
    FREEIF(gtmp_old2)
    MALLOC(gtmp_old2,(1:sp1,1:flavors))
  END IF

  IF ( op%rank .EQ. 0 ) THEN
    WRITE(op%ostream, '(1x,103A)') &
    "|----------------------------------------------------------------------------------------------------|"
    WRITE(op%ostream,'(1x,A)', ADVANCE="NO") "|"
  END IF

  total = DBLE(itotal)
  !write(std_out,*) "itotal",itotal
  indDensity = 1
  !write(std_out,*) "op%stats",op%stats
  DO isweep = 1, itotal
  !ii if(op%prtopt==1) write(std_out,*) "======== Isweep = ",isweep
    !updated_seg=.FALSE.
    DO iflavor = 1, flavors
     ! if(isweep==itotal) write(std_out,*) "    Iflavor = ",iflavor,op%Impurity%Particles(iflavor)%tail
   !ii if(op%prtopt==1)  write(std_out,*) "      ===Iflavor = ",iflavor
      op%Impurity%activeFlavor=iflavor
      op%Bath%activeFlavor=iflavor ; op%Bath%MAddFlag= .FALSE. ; op%Bath%MRemoveFlag = .FALSE.

      !write(std_out,*) "before tryaddremove"

      ! For iflavor, Try a move
      !==========================
      CALL Ctqmcoffdiag_tryAddRemove(op,updated_seg)
    !sui!write(std_out,*) "after tryaddremove",updated_seg

      updated = updated_seg .OR.  updated_swap(iflavor).OR.(isweep==1)
      updated_swap(iflavor) = .FALSE.
      if ( op%opt_nondiag >0 )  iflavor_d=0
      if ( op%opt_nondiag==0 )  iflavor_d=iflavor
      CALL GreenHyboffdiag_measHybrid(op%Greens, op%Bath%M, op%Impurity%Particles, updated,op%signvalue,iflavor_d) 

      CALL Ctqmcoffdiag_measN        (op, iflavor, updated)
      IF ( op%opt_analysis .EQ. 1 ) &
        CALL Ctqmcoffdiag_measCorrelation (op, iflavor)
      IF ( op%opt_order .GT. 0 ) &
        CALL Ctqmcoffdiag_measPerturbation(op, iflavor)
    END DO
    !CALL GreenHyboffdiag_measHybrid(op%Greens, op%Bath%M, op%Impurity%Particles, updated,op%signvalue,iflavor_d) 
    !DO iflavor = 1,flavors
    !  CALL Ctqmcoffdiag_measN        (op, iflavor, updated)
    !END DO

    IF ( MOD(isweep,modGlobalMove) .EQ. 0 ) THEN
  ! !sui!write(std_out,*) "isweep,modGlobalMove,inside",isweep,modGlobalMove
      CALL Ctqmcoffdiag_trySwap(op,swapUpdate1, swapUpdate2)
     ! !write(std_out,*) "no global move yet for non diag hybridization"
      IF ( swapUpdate1 .NE. 0 .AND. swapUpdate2 .NE. 0 ) THEN
        updated_swap(swapUpdate1) = .TRUE.
        updated_swap(swapUpdate2) = .TRUE.
      END IF
    END IF
    
    IF ( MOD(isweep,measurements) .EQ. 0 ) THEN ! default is always 
      CALL ImpurityOperator_measDE(op%Impurity,op%measDE)
      IF ( op%opt_spectra .GE. 1 .AND. MOD(isweep,measurements*op%opt_spectra) .EQ. 0 ) THEN
        op%density(1:flavors,indDensity) = op%measN(3,1:flavors)
        indDensity = indDensity+1
      END IF
    END IF

    IF ( MOD(isweep, modNoise1) .EQ. 0 ) THEN
      !modNext = isweep + modNoise2
      NRJ_new = op%measDE(1,1)
      CALL Vector_pushBack(op%measNoise(1),NRJ_new - NRJ_old1)
      NRJ_old1 = NRJ_new

      !! Try to limit accumulation error
      CALL ImpurityOperator_cleanOverlaps(op%Impurity)

      IF ( op%opt_noise .EQ. 1 ) THEN
        DO ifl1 = 1, flavors
          DO ind = 1, op%Greens%map(ifl1,ifl1)%tail
            itau = op%Greens%map(ifl1,ifl1)%listINT(ind)
            gtmp_new(itau,ifl1) = op%Greens%oper(itau,ifl1,ifl1) & 
                     +op%Greens%map(ifl1,ifl1)%listDBLE(ind)*DBLE(op%Greens%factor)
          END DO
          DO itau = 1, sp1
           CALL Vector_pushBack(op%measNoiseG(itau,ifl1,1), gtmp_new(itau,ifl1) - gtmp_old1(itau,ifl1))
           gtmp_old1(itau,ifl1) = gtmp_new(itau,ifl1)
          END DO
        END DO
      END IF
    END IF

    IF ( MOD(isweep,modNoise2) .EQ. 0 ) THEN
      NRJ_new = op%measDE(1,1)
      CALL Vector_pushBack(op%measNoise(2),NRJ_new - NRJ_old2)
      NRJ_old2 = NRJ_new
      IF ( op%opt_noise .EQ. 1 ) THEN
        DO ifl1 = 1, flavors
          DO ind = 1, op%Greens%map(ifl1,ifl1)%tail
            itau = op%Greens%map(ifl1,ifl1)%listINT(ind)
            gtmp_new(itau,ifl1) = op%Greens%oper(itau,ifl1,ifl1) & 
                  +op%Greens%map(ifl1,ifl1)%listDBLE(ind)*op%Greens%factor
          END DO
          DO itau = 1, sp1
            CALL Vector_pushBack(op%measNoiseG(itau,ifl1,2), gtmp_new(itau,ifl1) - gtmp_old2(itau,ifl1))
            gtmp_old2(itau,ifl1) = gtmp_new(itau,ifl1)
          END DO
        END DO 
      END IF

      IF ( op%rank .EQ. 0 ) THEN 
        new_percent = CEILING(DBLE(isweep)*100.d0/DBLE(itotal))
        DO ipercent = old_percent+1, new_percent 
          WRITE(op%ostream,'(A)',ADVANCE="NO") "-"
        END DO
        old_percent = new_percent
      END IF
    END IF

    IF ( op%opt_movie .EQ. 1 ) THEN
      WRITE(ilatex,'(A11,I9)') "%iteration ", isweep
      CALL ImpurityOperator_printLatex(op%Impurity,ilatex,isweep)
    END IF

  END DO

  IF ( op%rank .EQ. 0 ) THEN
    DO ipercent = old_percent+1, 100
      WRITE(op%ostream,'(A)',ADVANCE="NO") "-"
    END DO
    WRITE(op%ostream,'(A)') "|"
  END IF
 
  FREE(gtmp_new)
  FREE(gtmp_old1)
  FREE(gtmp_old2)
  FREE(updated_swap)

  IF ( op%opt_spectra .GE. 1 .AND. itotal .EQ. op%sweeps ) THEN
    IF ( endDensity .NE. indDensity-1 ) THEN
      op%density(:,endDensity) = -1.d0
    END IF
  END IF

  CALL CPU_TIME(cpu_time2)

  op%runTime = (cpu_time2 - cpu_time1)*1.05d0 ! facteur arbitraire de correction
END SUBROUTINE Ctqmcoffdiag_loop
!!***



!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_tryAddRemove
!! NAME
!!  Ctqmcoffdiag_tryAddRemove
!!
!! FUNCTION
!!  Try to add or remove a segment and an anti-segment
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_tryAddRemove(op,updated)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)             , INTENT(INOUT) :: op
!  TYPE(BathOperatoroffdiag)    , INTENT(INOUT) :: Bath 
!  TYPE(ImpurityOperator), INTENT(INOUT) :: Impurity 
  LOGICAL               , INTENT(  OUT) :: updated
!Local variables ------------------------------
  INTEGER                               :: position
  INTEGER         , DIMENSION(1:2)     :: nature ! -2 for antiseg and 1 for seg
  INTEGER                               :: i! -2 for antiseg and 1 for seg
  !INTEGER                               :: it,it1 !ii,
  DOUBLE PRECISION                      :: action
  DOUBLE PRECISION                      :: beta
  DOUBLE PRECISION                      :: time1
  DOUBLE PRECISION                      :: time2
  DOUBLE PRECISION                      :: time_avail
  DOUBLE PRECISION                      :: det_ratio,sign_det_ratio
  DOUBLE PRECISION                      :: overlap
  DOUBLE PRECISION                      :: length
  DOUBLE PRECISION                      :: signe
  DOUBLE PRECISION                      :: tail
  INTEGER                      :: tailint
  DOUBLE PRECISION                      :: signdet, signdetprev
  DOUBLE PRECISION, DIMENSION(1:2)      :: CdagC_1

  IF ( .NOT. op%set ) &
    CALL ERROR("Ctqmcoffdiag_trySegment : QMC not set                       ")

        !write(std_out,*) "      TryAddRemove start"
  nature(1) = CTQMC_SEGME
  nature(2) = CTQMC_ANTIS
  beta      = op%beta

  updated = .FALSE.
  tailint  = (op%Impurity%particles(op%Impurity%activeFlavor)%tail)
  tail  = DBLE(tailint)
  !write(std_out,*) "op%Impurity%particles(op%Impurity%activeFlavor)%tail",op%Impurity%activeFlavor,tail


  !=====================================
  ! First choose segment or antisegment
  !=====================================
  DO i = 1, 2
    signe = SIGN(1.d0,DBLE(nature(i))) 
!      -----  1: segment        signe= 1  ( CTQMC_SEGME =  1 )
!      -----  2: antisegment    signe=-1  ( CTQMC_ANTIS = -2 )
!    NB: Sign(a,b) = sign(b) * a

 !prt!if(op%prtopt==1) write(std_out,*) "       ==Starting configuration",i
 !prt!if(op%prtopt==1) write(std_out,*) "        = Segments:"
    tailint  = (op%Impurity%particles(op%Impurity%activeFlavor)%tail)
!prt!    do ii=0, op%Impurity%Particles(op%Impurity%activeFlavor)%tail
 !prt!if(op%prtopt==1)  write(std_out,*) ii, op%Impurity%Particles(op%Impurity%activeFlavor)%list(ii,1), &
!prt!&                    op%Impurity%Particles(op%Impurity%activeFlavor)%list(ii,2)
!prt!    enddo
  !sui!write(std_out,*) "        = M Matrix",op%Bath%sumtails
!prt!    do it=1,op%Bath%sumtails
    !sui!write(std_out,'(a,3x,500e10.3)') "        M start",(op%Bath%M%mat(it,it1),it1=1,op%Bath%sumtails)
!prt!    enddo
    CALL OurRng(op%seed,action)

    !==========================
    ! Add segment/antisegment
    !==========================
    IF ( action .LT. .5d0 ) THEN ! Add a segment or antisegment
     !ii  write(std_out,*) "        =try: Segment added of type",i,op%prtopt

      ! Select time1 (>0) in [0,beta]
      !==============================
      CALL OurRng(op%seed,time1)
      time1 = time1 * beta

      ! time_avail is the distance between between time1 and 
      !   - the next start of a segment for a segment addition
      !   - the next end of a segment for an antisegment addition
      ! ImpurityOperator_getAvailableTime > 0 for a segment      (signe>0) -> time_avail>0
      ! ImpurityOperator_getAvailableTime < 0 for an antisegment (signe<0) -> time_avail>0
      !====================================================================
      time_avail = ImpurityOperator_getAvailableTime(op%Impurity,time1,position) * signe
     !ii  write(std_out,*) "        =try: time_avail",time_avail,time1
      IF ( time_avail .GT. 0.d0 ) THEN

        ! Time2 is  the length of the proposed new (anti)segment
        !=======================================================
        CALL OurRng(op%seed,time2)
        IF ( time2 .EQ. 0.d0 ) CALL OurRng(op%seed,time2) ! Prevent null segment

        ! Now time2 is the time at the end of the proposed new (anti) segment
        ! time2 > time1 
        !====================================================================
        time2     = time1 + time2 * time_avail
      !sui!write(std_out,*) tailint+1,time1,time2,position
!        CALL CdagC_init(CdagC_1,time1,time2)

        ! CdagC_1 gives the stard/end times for the proposed new segment/antisegment
        ! CdagC1(C_) can  be above beta.
        !  For a      segment CdagC_1(Cdag_) = time1 < CdagC_1(C_) = time2, l=time2-time1 > 0
        !  For a anti segment CdagC_1(Cdag_) = time2 > CdagC_1(C_) = time1, l=time1-time2 < 0
        !  time2 can be above beta and thus for a      segment CdagC_1(C_   ) > beta
        !  time2 can be above beta and thus for an antisegment CdagC_1(Cdag_) > beta
        !  length > 0 for     segment
        !  length < 0 for antisegment
        !====================================================================================
        CdagC_1(Cdag_) = ((1.d0+signe)*time1+(1.d0-signe)*time2)*0.5d0
        CdagC_1(C_   ) = ((1.d0+signe)*time2+(1.d0-signe)*time1)*0.5d0
!        length    = CdagC_length(CdagC_1)
        length    = CdagC_1(C_   ) - CdagC_1(Cdag_)
        !write(std_out,*) "      try : times", CdagC_1(C_   ),CdagC_1(Cdag_)
        !write(std_out,*) "      length", length

!      -----  Computes the determinant ratio
        det_ratio = BathOperatoroffdiag_getDetAdd(op%Bath,CdagC_1,position,op%Impurity%particles) 

!      -----  Computes the overlap
        overlap   = ImpurityOperator_getNewOverlap(op%Impurity,CdagC_1)
        signdetprev  = ImpurityOperator_getsign(op%Impurity, time2, i, action, position)

        !write(std_out,*) "      overlap   ", overlap
        CALL OurRng(op%seed,time1)
        !write(std_out,*) "      Rnd", time1
        signdet=1.d0
        det_ratio=det_ratio*signdetprev
                 
        IF ( det_ratio .LT. 0.d0 ) THEN
        !sui!write(std_out,*) "         NEGATIVE DET",det_ratio,signdetprev
          det_ratio   = - det_ratio
          sign_det_ratio=-1
          op%stats(nature(i)+CTQMC_DETSI) = op%stats(nature(i)+CTQMC_DETSI) + 1.d0
       !  op%signvaluecurrent=-1.d0
        ELSE
          sign_det_ratio=1
       !  op%signvaluecurrent=+1.d0
         ! signdet=-1.d0
        !sui!write(std_out,*) "                  DET",det_ratio,signdetprev
        END IF
      !ii  write(std_out,*) "                  DET",det_ratio
       ! op%signvaluemeas=op%signvaluemeas+1.d0
        !write(std_out,*) "        .................",(time1 * (tail + 1.d0 )),beta * time_avail * det_ratio * DEXP(op%mu(op%Impurity%activeFlavor)*length + overlap)
        !write(std_out,*) "        .................",beta , time_avail , op%mu(op%Impurity%activeFlavor),op%Impurity%activeFlavor

        IF ( (time1 * (tail + 1.d0 )) &
             .LT. (beta * time_avail * det_ratio * DEXP(op%mu(op%Impurity%activeFlavor)*length + overlap) ) ) THEN
!          write(*,*) "before"
!          CALL ListCdagCoffdiag_print(op%Impurity%particles(op%Impurity%activeFlavor),6)
          CALL ImpurityOperator_add(op%Impurity,CdagC_1,position)
!          write(*,*) "after "
!          CALL ListCdagCoffdiag_print(op%Impurity%particles(op%Impurity%activeFlavor),6)
          CALL BathOperatoroffdiag_setMAdd(op%bath,op%Impurity%particles) 
          op%stats(nature(i)+CTQMC_ADDED) = op%stats(nature(i)+CTQMC_ADDED)  + 1.d0
          updated = .TRUE. .OR. updated
          tail = tail + 1.d0
          tailint = tailint + 1
!          read(*,*) time1
          !ii  write(6,*) "        Accepted addition, new conf is",time1
         !prt! do ii=0, op%Impurity%Particles(op%Impurity%activeFlavor)%tail
          !prt!if(op%prtopt==1)  write(6,*) ii, op%Impurity%Particles(op%Impurity%activeFlavor)%list(ii,1),&
!prt!&                          op%Impurity%Particles(op%Impurity%activeFlavor)%list(ii,2)
         !prt! enddo
        !sui!write(6,*) "        = M Matrix"
         !prt! do it=1,op%Bath%sumtails
          !sui!write(6,*) "        M new",(op%Bath%M%mat(it,it1),it1=1,op%Bath%sumtails)
         !prt! enddo

          IF ( sign_det_ratio .LT. 0.d0 ) op%signvalue=-op%signvalue
        !sui!write(6,*) "                  signvalue",op%signvalue
        ELSE
     !ii  write(6,*) "        Refused      addition: proba",time1
        END IF 
      ELSE
    !sui!write(6,*) "        Refused      addition: time_avail <0"
      END IF 

    !========================================
    ! Remove segment/antisegment
    !========================================
    ELSE ! Remove a segment among the segment of the flavor activeflavor
      !ii if(op%prtopt==1)  write(6,*) "        =try: Segment removed of type",i
      IF ( tail .GT. 0.d0 ) THEN
        CALL OurRng(op%seed,time1)
        position = INT(((time1 * tail) + 1.d0) * signe )
        !prt!if(op%prtopt==1)  write(6,*) "         position",position 
        time_avail = ImpurityOperator_getAvailedTime(op%Impurity,position)
        det_ratio  = BathOperatoroffdiag_getDetRemove(op%Bath,position)
        !write(6,*) "        det_ratio", det_ratio
        CdagC_1    = ImpurityOperator_getSegment(op%Impurity,position)
!        length     = CdagC_length(CdagC_1)
        length     = CdagC_1(C_) - CdagC_1(Cdag_)
        !write(6,*) "        length   ", length
        overlap    = ImpurityOperator_getNewOverlap(op%Impurity,CdagC_1)
        !write(6,*) "        overlap  ", overlap
        CALL OurRng(op%seed,time1)
        !write(6,*) "        Random   ",time1
        signdetprev = ImpurityOperator_getsign(op%Impurity, time2, i, action, position)
        det_ratio=det_ratio*signdetprev
        signdet=1.d0
        IF ( det_ratio .LT. 0.d0 ) THEN
        !sui!write(6,*) "         NEGATIVE DET",det_ratio,signdetprev
          det_ratio   = -det_ratio
          sign_det_ratio=-1
!          op%seg_sign = op%seg_sign + 1.d0
          op%stats(nature(i)+CTQMC_DETSI) = op%stats(nature(i)+CTQMC_DETSI) + 1.d0
          signdet=-1.d0
        ELSE 
          sign_det_ratio=1
        !sui!write(6,*) "                  DET",det_ratio,signdetprev
        END IF
       !ii  write(6,*) "                  DET",det_ratio
        IF ( (time1 * beta * time_avail * DEXP(op%mu(op%Impurity%activeFlavor)*length+overlap)) &
             .LT. (tail * det_ratio ) ) THEN
          CALL ImpurityOperator_remove(op%Impurity,position)
          CALL BathOperatoroffdiag_setMRemove(op%Bath,op%Impurity%particles) 
          !op%seg_removed = op%seg_removed  + 1.d0
          op%stats(nature(i)+CTQMC_REMOV) = op%stats(nature(i)+CTQMC_REMOV)  + 1.d0
          updated = .TRUE. .OR. updated
          tail = tail -1.d0
          tailint = tailint -1
          !ii  write(6,*) "        Accepted removal, new conf is:",time1
      !prt!    do ii=0, op%Impurity%Particles(op%Impurity%activeFlavor)%tail
          !prt!if(op%prtopt==1)  write(6,*) ii, op%Impurity%Particles(op%Impurity%activeFlavor)%list(ii,1),&
!prt!&                          op%Impurity%Particles(op%Impurity%activeFlavor)%list(ii,2)
      !prt!    enddo
        !sui!write(6,*) "        = M Matrix"
       !prt!   do it=1,op%Bath%sumtails
          !sui!write(6,*) "        M new",(op%Bath%M%mat(it,it1),it1=1,op%Bath%sumtails)
        !prt!  enddo
          IF ( sign_det_ratio .LT. 0.d0 ) op%signvalue=-op%signvalue
        !sui!write(6,*) "                  signvalue",op%signvalue
        ELSE
     !ii  write(6,*) "        Refused      removal",time1
        END IF
      ELSE
      !sui!write(6,*) "        Refused      removal: no segment available"
      END IF
    END IF
    !========================================
    ! End Add/Remove Antisegment
    !========================================
  END DO
END SUBROUTINE Ctqmcoffdiag_tryAddRemove
!!***

!!SUBROUTINE Ctqmcoffdiag_trySegment(op,updated)
!!!  TYPE(BathOperatoroffdiag)    , INTENT(INOUT) :: Bath 
!!!  TYPE(ImpurityOperator), INTENT(INOUT) :: Impurity 
!!  LOGICAL               , INTENT(INOUT) :: updated
!!  INTEGER                               :: position
!!  DOUBLE PRECISION                      :: action
!!  DOUBLE PRECISION                      :: beta
!!  DOUBLE PRECISION                      :: time1
!!  DOUBLE PRECISION                      :: time2
!!  DOUBLE PRECISION                      :: time_avail
!!  DOUBLE PRECISION                      :: det_ratio
!!  DOUBLE PRECISION                      :: overlap
!!  DOUBLE PRECISION                      :: length
!!  DOUBLE PRECISION                      :: tail
!!  DOUBLE PRECISION, DIMENSION(1:2)      :: CdagC_1
!!
!!  IF ( .NOT. op%set ) &
!!    CALL ERROR("Ctqmcoffdiag_trySegment : QMC not set                       ")
!!
!!  beta     =  op%beta
!!  tail     = DBLE(op%Impurity%particles(op%Impurity%activeFlavor)%tail)
!!
!!  CALL RANDOM_NUMBER(action)
!!  
!!  updated = .FALSE.
!!
!!  IF ( action .LT. .5d0 ) THEN ! Ajout de segment
!!    CALL RANDOM_NUMBER(time1)
!!    time1 = time1 * beta
!!    time_avail = ImpurityOperator_getAvailableTime(op%Impurity,time1,position)
!!    IF ( time_avail .GT. 0.d0 ) THEN
!!      CALL RANDOM_NUMBER(time2)
!!      time2     = time1 + time2 * time_avail
!!!      CALL CdagC_init(CdagC_1,time1,time2)
!!      CdagC_1(Cdag_) = time1
!!      CdagC_1(C_   ) = time2
!!!      length    = CdagC_length(CdagC_1)
!!      length    = time2 - time1 
!!      det_ratio = BathOperatoroffdiag_getDetAdd(op%Bath,CdagC_1,position,op%Impurity%particles(op%Impurity%activeFlavor))
!!      overlap   = ImpurityOperator_getNewOverlap(op%Impurity,CdagC_1)
!!      CALL RANDOM_NUMBER(time1)
!!      IF ( det_ratio .LT. 0.d0 ) THEN
!!        det_ratio   = -det_ratio
!!        op%stats(CTQMC_SEGME+CTQMC_DETSI) = op%stats(CTQMC_SEGME+CTQMC_DETSI) + SIGN(1.d0,det_ratio)
!!      END IF
!!      IF ( (time1 * (tail + 1.d0 )) &
!!           .LT. (beta * time_avail * det_ratio * DEXP(op%mu(op%Impurity%activeFlavor)*length + overlap) ) ) THEN
!!        write(*,*) position
!!        CALL ImpurityOperator_add(op%Impurity,CdagC_1,position)
!!        CALL BathOperatoroffdiag_setMAdd(op%bath,op%Impurity%particles(op%Impurity%activeFlavor))
!!        op%stats(CTQMC_SEGME+CTQMC_ADDED) = op%stats(CTQMC_SEGME+CTQMC_ADDED)  + 1.d0
!!        updated = .TRUE.
!!      END IF 
!!    END IF 
!!
!!  ELSE ! Supprimer un segment
!!    IF ( tail .GT. 0.d0 ) THEN
!!      CALL RANDOM_NUMBER(time1)
!!      position = INT(time1 * tail) + 1
!!      time_avail = ImpurityOperator_getAvailedTime(op%Impurity,position)
!!      det_ratio  = BathOperatoroffdiag_getDetRemove(op%Bath,position)
!!      CdagC_1    = ImpurityOperator_getSegment(op%Impurity,position)
!!!      length     = CdagC_length(CdagC_1)
!!      length     = CdagC_1(C_) - CdagC_1(Cdag_)
!!      overlap    = ImpurityOperator_getNewOverlap(op%Impurity,CdagC_1)
!!      CALL RANDOM_NUMBER(time1)
!!      IF ( det_ratio .LT. 0.d0 ) THEN
!!        det_ratio   = -det_ratio
!!!        op%seg_sign = op%seg_sign + 1.d0
!!        op%stats(CTQMC_SEGME+CTQMC_DETSI) = op%stats(CTQMC_SEGME+CTQMC_DETSI) + SIGN(1.d0,det_ratio)
!!      END IF
!!      IF ( (time1 * beta * time_avail * DEXP(op%mu(op%Impurity%activeFlavor)*length+overlap)) &
!!           .LT. (tail * det_ratio ) ) THEN
!!        write(*,*) position
!!        CALL ImpurityOperator_remove(op%Impurity,position)
!!        CALL BathOperatoroffdiag_setMRemove(op%Bath,op%Impurity%particles(op%Impurity%activeFlavor))
!!        !op%seg_removed = op%seg_removed  + 1.d0
!!        op%stats(CTQMC_SEGME+CTQMC_REMOV) = op%stats(CTQMC_SEGME+CTQMC_REMOV)  + 1.d0
!!        updated = .TRUE.
!!      END IF
!!    END IF
!!  END IF
!!END SUBROUTINE Ctqmcoffdiag_trySegment
!!***
!!
!!SUBROUTINE Ctqmcoffdiag_tryAntiSeg(op, updated)
!!!  TYPE(BathOperatoroffdiag)    , INTENT(INOUT) :: Bath 
!!!  TYPE(ImpurityOperator), INTENT(INOUT) :: Impurity 
!!  LOGICAL               , INTENT(INOUT) :: updated
!!  INTEGER                               :: position
!!  DOUBLE PRECISION                      :: action
!!  DOUBLE PRECISION                      :: beta
!!  DOUBLE PRECISION                      :: time1
!!  DOUBLE PRECISION                      :: time2
!!  DOUBLE PRECISION                      :: time_avail
!!  DOUBLE PRECISION                      :: det_ratio
!!  DOUBLE PRECISION                      :: overlap
!!  DOUBLE PRECISION                      :: length
!!  DOUBLE PRECISION                      :: tail
!!  DOUBLE PRECISION, DIMENSION(1:2)      :: CdagC_1
!!
!!  IF ( .NOT. op%set ) &
!!    CALL ERROR("Ctqmcoffdiag_trySegment : QMC not set                       ")
!!
!!  beta     =  op%beta
!!  tail     =  DBLE(op%Impurity%particles(op%Impurity%activeFlavor)%tail)
!!
!!  CALL RANDOM_NUMBER(action)
!!
!!  updated = .FALSE.
!!
!!  IF ( action .LT. .5d0 ) THEN ! Ajout d'un antiseg
!!    CALL RANDOM_NUMBER(time1)
!!    time1 = time1 * beta
!!    time_avail = ImpurityOperator_getAvailableTime(op%Impurity,time1,position)
!!    IF ( time_avail .LT. 0.d0 ) THEN
!!      CALL RANDOM_NUMBER(time2)
!!      time2     = time1 - time2 * time_avail
!!!      CALL CdagC_init(CdagC_1,time2,time1)
!!      CdagC_1(Cdag_) = time2
!!      CdagC_1(C_   ) = time1
!!!      length    = CdagC_length(CdagC_1) ! /!\ length is negative
!!      length = time1 - time2
!!      det_ratio = BathOperatoroffdiag_getDetAdd(op%Bath,CdagC_1,position,op%Impurity%particles(op%Impurity%activeFlavor))
!!      overlap   = ImpurityOperator_getNewOverlap(op%Impurity,CdagC_1) ! OK
!!      CALL RANDOM_NUMBER(time1)
!!      IF ( det_ratio .LT. 0.d0 ) THEN
!!        det_ratio    = -det_ratio
!!!        op%anti_sign = op%anti_sign + 1.d0
!!        op%stats(CTQMC_ANTIS+CTQMC_DETSI) = op%stats(CTQMC_ANTIS+CTQMC_DETSI) + SIGN(1.d0,det_ratio)
!!      END IF
!!      IF ( (time1 * (tail + 1.d0 )) & 
!!           .LT. (beta * ABS(time_avail) * det_ratio * DEXP(op%mu(op%Impurity%activeFlavor)*length + overlap) ) ) THEN
!!        CALL ImpurityOperator_add(op%Impurity,CdagC_1,position) 
!!        !write(*,*) position
!!        CALL BathOperatoroffdiag_setMAdd(op%bath,op%Impurity%particles(op%Impurity%activeFlavor)) 
!!        !op%anti_added = op%anti_added  + 1.d0
!!        op%stats(CTQMC_ANTIS+CTQMC_ADDED) = op%stats(CTQMC_ANTIS+CTQMC_ADDED)  + 1.d0
!!        updated = .TRUE.
!!      END IF 
!!    END IF 
!!  ELSE ! Supprimer un antiseg
!!    IF ( tail .GT. 0.d0 ) THEN
!!      CALL RANDOM_NUMBER(time1)
!!      position = -(INT(time1 * tail) + 1)
!!      time_avail = ImpurityOperator_getAvailedTime(op%Impurity,position)!OK
!!      det_ratio  = BathOperatoroffdiag_getDetRemove(op%Bath,position)!OK
!!      CdagC_1    = ImpurityOperator_getSegment(op%Impurity,position)!OK 
!!!      length     = CdagC_length(CdagC_1) ! /!\ length is negative
!!      length = CdagC_1(C_) - CdagC_1(Cdag_)
!!      overlap    = ImpurityOperator_getNewOverlap(op%Impurity,CdagC_1) !OK
!!      CALL RANDOM_NUMBER(time1)
!!      IF ( det_ratio .LT. 0.d0 ) THEN
!!        det_ratio   = -det_ratio
!!!        op%anti_sign = op%anti_sign + 1.d0
!!        op%stats(CTQMC_ANTIS+CTQMC_DETSI) = op%stats(CTQMC_ANTIS+CTQMC_DETSI) + SIGN(1.d0,det_ratio)
!!      END IF
!!      IF ( (time1 * beta * time_avail * DEXP(op%mu(op%Impurity%activeFlavor)*length+overlap)) &
!!           .LT. (tail * det_ratio ) ) THEN
!!        CALL ImpurityOperator_remove(op%Impurity,position)
!!        !write(*,*) position
!!        CALL BathOperatoroffdiag_setMRemove(op%Bath,op%Impurity%particles(op%Impurity%activeFlavor))
!!        !op%anti_removed = op%anti_removed  + 1.d0
!!        op%stats(CTQMC_ANTIS+CTQMC_REMOV) = op%stats(CTQMC_ANTIS+CTQMC_REMOV)  + 1.d0
!!        updated = .TRUE.
!!      END IF
!!    END IF
!!  END IF
!!END SUBROUTINE Ctqmcoffdiag_tryAntiSeg
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_trySwap
!! NAME
!!  Ctqmcoffdiag_trySwap
!!
!! FUNCTION
!!  try a global move (swap to flavors)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_trySwap(op,flav_i,flav_j)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)           , INTENT(INOUT) :: op
!  TYPE(BathOperatoroffdiag)    , INTENT(INOUT) :: Bath 
!  TYPE(ImpurityOperator), INTENT(INOUT) :: Impurity 
  INTEGER               , INTENT(  OUT) :: flav_i
  INTEGER               , INTENT(  OUT) :: flav_j
!Local variables ------------------------------
  INTEGER :: flavor_i
  INTEGER :: flavor_j !,ii,it,it1 !,iflavor
  DOUBLE PRECISION :: rnd
  DOUBLE PRECISION :: lengthi
  DOUBLE PRECISION :: lengthj
  DOUBLE PRECISION :: overlapic1
  DOUBLE PRECISION :: overlapjc1
  DOUBLE PRECISION :: overlapic2
  DOUBLE PRECISION :: overlapjc2
  !DOUBLE PRECISION :: detic1
  !DOUBLE PRECISION :: detjc1
  !DOUBLE PRECISION :: detic2
  !DOUBLE PRECISION :: detjc2
  DOUBLE PRECISION :: det_ratio,detnew,detold
  DOUBLE PRECISION :: local_ratio
 ! TYPE(BathOperatoroffdiag)  :: Bathnew


  !CALL RANDOM_NUMBER(rnd)
  CALL OurRng(op%seed,rnd)
  flavor_i = NINT(rnd*DBLE(op%flavors-1.d0))+1
  !CALL RANDOM_NUMBER(rnd)
  CALL OurRng(op%seed,rnd)
  flavor_j = NINT(rnd*DBLE(op%flavors-1.d0))+1
  !ii write(6,'(a,2i4)') "--------------- new swap --------------------------------",flavor_i,flavor_j
  
  flav_i = 0
  flav_j = 0
  !ii   do iflavor=1,op%flavors
  !ii     write(6,*) "BEFORE  GMOVE For flavor", iflavor,"size is",op%Impurity%particles(iflavor)%tail," and  Conf is :"
  !ii     do ii=1, op%Impurity%Particles(iflavor)%tail
  !ii       write(6,'(i4,100f12.3)') ii, op%Impurity%Particles(iflavor)%list(ii,1),&
  !ii  &                   op%Impurity%Particles(iflavor)%list(ii,2)
  !ii     enddo
  !ii   enddo
  !ii   write(6,*) "        = M Matrix"
  !ii   write(6,'(a,2x,100(i12))') "Flavor=",((iflavor,it=1,op%Impurity%particles(iflavor)%tail),iflavor=1,op%flavors)
  !ii   write(6,'(i21,100i12)') ((it,it=1,op%Impurity%particles(iflavor)%tail),iflavor=1,op%flavors)
  !ii   do it=1,op%Bath%sumtails
  !ii     write(6,'(a,100f12.3)') "        M before",(op%Bath%M%mat(it,it1),it1=1,op%Bath%sumtails)
  !ii   enddo

  ! todoba this part
  IF ( flavor_i .NE. flavor_j ) THEN
    !CALL BathOperatoroffdiag_init(Bathnew, op%flavors, op%samples, op%beta, 0, op%opt_nondiag)
    ! On tente d'intervertir i et j
    ! Configuration actuelle :

    op%modGlobalMove(2) = op%modGlobalMove(2)+1
    ! ===========================================
    ! First use M matrix to compute determinant
    ! ===========================================
    detold     = BathOperatoroffdiag_getDetF(op%Bath) ! use op%Bath%M

    ! ===========================================
    ! Second build M_update matrix to compute determinant after.
    ! ===========================================
    !CALL ListCdagCoffdiag_print(particle)
    call BathOperatoroffdiag_recomputeM(op%Bath,op%impurity%particles,flavor_i,flavor_j) ! compute op%Bath%M_update
    detnew     = BathOperatoroffdiag_getDetF(op%Bath,option=1) ! use op%Bath%M_update

    lengthi    = ImpurityOperator_measN(op%Impurity,flavor_i)
    lengthj    = ImpurityOperator_measN(op%Impurity,flavor_j)
    overlapic1 = ImpurityOperator_overlapFlavor(op%Impurity,flavor_i)
    overlapjc1 = ImpurityOperator_overlapFlavor(op%Impurity,flavor_j)
    ! lengths unchanged
    overlapic2 = ImpurityOperator_overlapSwap(op%Impurity,flavor_i,flavor_j)
    overlapjc2 = ImpurityOperator_overlapSwap(op%Impurity,flavor_j,flavor_i)

!    IF ( detic1*detjc1 .EQ. detic2*detjc2 ) THEN
!      det_ratio = 1.d0
!    ELSE IF ( detic1*detjc1 .EQ. 0.d0 ) THEN
!      det_ratio = detic2*detjc2 ! evite de diviser par 0 si pas de segment
!    ELSE

    det_ratio = detnew/detold ! because the determinant is the determinant of F
   !ii  write(6,*) "det_ratio, detold,detnew",det_ratio, detold,detnew, detold/detnew

!    END IF
    local_ratio = DEXP(-overlapic2*overlapjc2+overlapic1*overlapjc1 &
                      +(lengthj-lengthi)*(op%mu(flavor_i)-op%mu(flavor_j)))
   !ii  write(6,*) "local_ratio",local_ratio

    ! Wloc = exp(muN-Uo)
    !CALL RANDOM_NUMBER(rnd)
    CALL OurRng(op%seed,rnd)
    IF ( rnd .LT. local_ratio*det_ratio ) THEN ! swap accepted
   !ii    write(6,*) "        = M Matrix before swap"
   !ii    write(6,'(a,2x,100(i12))') "Flavor=",((iflavor,it=1,op%Impurity%particles(iflavor)%tail),iflavor=1,op%flavors)
   !ii    write(6,'(i21,100i12)') ((it,it=1,op%Impurity%particles(iflavor)%tail),iflavor=1,op%flavors)
   !ii    do it=1,op%Bath%sumtails
   !ii      write(6,'(a,100f12.3)') "        M after ",(op%Bath%M%mat(it,it1),it1=1,op%Bath%sumtails)
   !ii    enddo
   !ii    do it=1,op%Bath%sumtails
   !ii      write(6,'(a,100f12.3)') " update M after ",(op%Bath%M_update%mat(it,it1),it1=1,op%Bath%sumtails)
   !ii    enddo
   !ii    write(6,*) "Gmove accepted",rnd,local_ratio*det_ratio
      CALL ImpurityOperator_swap(op%Impurity, flavor_i,flavor_j)
      CALL BathOperatoroffdiag_swap    (op%Bath    , flavor_i,flavor_j) !  use op%Bath%M_update to built new op%Bath%M
      
      op%swap = op%swap + 1.d0
      flav_i = flavor_i
      flav_j = flavor_j
    ELSE
   !ii   write(6,*) "Gmove refused",rnd,local_ratio*det_ratio
!      CALL WARN("Swap refused")
!      WRITE(op%ostream,'(6E24.14)') local_ratio, det_ratio, detic1, detjc1, detic2, detjc2
    END IF
   ! CALL BathOperatoroffdiag_destroy(Bathnew)
  END IF
 !ii  do iflavor=1,op%flavors
 !ii    write(6,*) "AFTER   GMOVE For flavor", iflavor,"size is",op%Impurity%particles(iflavor)%tail," and  Conf is :"
 !ii    do ii=1, op%Impurity%Particles(iflavor)%tail
 !ii      write(6,'(15x,i4,100f12.3)') ii, op%Impurity%Particles(iflavor)%list(ii,1),&
 !ii &                   op%Impurity%Particles(iflavor)%list(ii,2)
 !ii    enddo
 !ii  enddo
 !ii  write(6,*) "        = M Matrix"
 !ii  write(6,'(a,2x,100(i12))') "Flavor=",((iflavor,it=1,op%Impurity%particles(iflavor)%tail),iflavor=1,op%flavors)
 !ii  write(6,'(i21,100i12)') ((it,it=1,op%Impurity%particles(iflavor)%tail),iflavor=1,op%flavors)
 !ii  do it=1,op%Bath%sumtails
 !ii    write(6,'(a,100f12.3)') "        M after ",(op%Bath%M%mat(it,it1),it1=1,op%Bath%sumtails)
 !ii  enddo
 !ii  do it=1,op%Bath%sumtails
 !ii    write(6,'(a,100f12.3)') " update M after ",(op%Bath%M_update%mat(it,it1),it1=1,op%Bath%sumtails)
 !ii  enddo

END SUBROUTINE Ctqmcoffdiag_trySwap
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_measN
!! NAME
!!  Ctqmcoffdiag_measN
!!
!! FUNCTION
!!  measures the number of electron
!!  by taking into account the value for the move before before this one
!!  with the correct weight.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_measN(op, iflavor, updated)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)             , INTENT(INOUT)     :: op
  !TYPE(ImpurityOperator), INTENT(IN   )     :: impurity
  INTEGER               , INTENT(IN   )     :: iflavor
  LOGICAL               , INTENT(IN   )     :: updated

!  IF ( .NOT. op%set ) &
!    CALL ERROR("Ctqmcoffdiag_measN : QMC not set                           ")

  
  IF ( updated .EQV. .TRUE. ) THEN
!  --- accumulate occupations with values op%measN(3,iflavor) from the last measurements with the corresponding weight
!  ---  op*measN(4,iflavor)
    op%measN(1,iflavor) = op%measN(1,iflavor) + op%measN(3,iflavor)*op%measN(4,iflavor)

!  --- Compute total number of new measurements 
    op%measN(2,iflavor) = op%measN(2,iflavor) + op%measN(4,iflavor)

!  --- Compute the occupation for this configuration (will be put in
!  --- op%measN(1,iflavor) at the next occurence of updated=.true.), with
!  --- the corresponding weight  op%measN(4,iflavor) (we do not now it yet)
    op%measN(3,iflavor) = ImpurityOperator_measN(op%impurity)

!  --- set weight: as update=true, it is a new measurement , so put it to one
    op%measN(4,iflavor) = 1.d0

  ELSE
!  --- increased the count so that at new move, we will be able to update measN(1) correctly.
    op%measN(4,iflavor) = op%measN(4,iflavor) + 1.d0
  END IF
END SUBROUTINE Ctqmcoffdiag_measN
!!***

!#ifdef CTCtqmcoffdiag_ANALYSIS
!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_measCorrelation
!! NAME
!!  Ctqmcoffdiag_measCorrelation
!!
!! FUNCTION
!!  measure all correlations in times for a flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_measCorrelation(op, iflavor)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)             , INTENT(INOUT)       :: op
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

  IF ( .NOT. op%set ) &
    CALL ERROR("Ctqmcoffdiag_measCorrelation : QMC not set                 ")
    !write(6,*) "not available"
    stop

  size = op%impurity%particles(op%impurity%activeFlavor)%tail
  beta = op%beta

  IF ( size .EQ. 0 ) RETURN
  
  inv_dt = op%inv_dt

  DO iCdag = 1, size ! first segments
    tCdag  = op%impurity%particles(op%impurity%activeFlavor)%list(iCdag,Cdag_)
    tC     = op%impurity%particles(op%impurity%activeFlavor)%list(iCdag,C_   )
    index = INT( ( (tC - tCdag)  * inv_dt ) + .5d0 ) + 1
    op%measCorrelation(index,1,iflavor) = op%measCorrelation(index,1,iflavor) + 1.d0
    MODCYCLE(iCdag+1,size,iCdagBeta)
    index = INT( ( ( &
                    op%impurity%particles(op%impurity%activeFlavor)%list(iCdagBeta,Cdag_) - tC &
                    + AINT(DBLE(iCdag)/DBLE(size))*beta &
                   )  * inv_dt ) + .5d0 ) + 1
    IF ( index .LT. 1 .OR. index .GT. op%samples+1 ) THEN
      CALL WARN("Ctqmcoffdiag_measCorrelation : bad index line 1095         ")
    ELSE
      op%measCorrelation(index,2,iflavor) = op%measCorrelation(index,2,iflavor) + 1.d0
    END IF
!    DO iC = 1, size
!      tC = impurity%particles(impurity%activeFlavor)%list(C_,iC)
!      time = tC - tCdag
!      IF ( time .LT. 0.d0 ) time = time + beta
!      index = INT( ( time * inv_dt ) + .5d0 ) + 1
!      op%measCorrelation(index,3,iflavor) = op%measCorrelation(index,3,iflavor) + 1.d0
!    END DO
    DO iC = 1, size!  op%Greens(iflavor)%index_old%tail 
!todoba        op%measCorrelation(op%Greens(iflavor)%map%listINT(iC+(iCdag-1)*size),3,iflavor) = &
!todoba        op%measCorrelation(op%Greens(iflavor)%map%listINT(iC+(iCdag-1)*size),3,iflavor) + 1.d0
    END DO
  END DO

END SUBROUTINE Ctqmcoffdiag_measCorrelation
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_measPerturbation
!! NAME
!!  Ctqmcoffdiag_measPerturbation
!!
!! FUNCTION
!!  measure perturbation order
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_measPerturbation(op, iflavor)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)             , INTENT(INOUT)     :: op
  !TYPE(ImpurityOperator), INTENT(IN   )     :: impurity
  INTEGER               , INTENT(IN   )     :: iflavor
!Local variables ------------------------------
  INTEGER                                   :: index

  IF ( .NOT. op%set ) &
    CALL ERROR("Ctqmcoffdiag_measiPerturbation : QMC not set               ")

  index = op%impurity%particles(op%impurity%activeFlavor)%tail + 1
  IF ( index .LE. op%opt_order ) &
    op%measPerturbation(index,iflavor) = op%measPerturbation(index,iflavor) + 1.d0
  IF ( index == 1 ) THEN
    IF (op%impurity%particles(iflavor)%list(0,C_) < op%impurity%particles(iflavor)%list(0,Cdag_) ) THEN
      op%meas_fullemptylines(1,iflavor) = op%meas_fullemptylines(1,iflavor) + 1.d0
    ELSE
      op%meas_fullemptylines(2,iflavor) = op%meas_fullemptylines(2,iflavor) + 1.d0
    ENDIF
  ENDIF

END SUBROUTINE Ctqmcoffdiag_measPerturbation
!!***
!#endif

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_getResult
!! NAME
!!  Ctqmcoffdiag_getResult
!!
!! FUNCTION
!!  reduce everything to get the result of the simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_getResult(op)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)  , INTENT(INOUT)                    :: op
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
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: measN_1
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
  DOUBLE PRECISION                                       :: signvaluemeassum
!  INTEGER                                       :: fin
#ifdef HAVE_MPI
  INTEGER                                       :: ierr
#endif
  INTEGER                                       :: sizeoper,nbprocs,myrank
  DOUBLE PRECISION                              :: inv_size
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: buffer 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: buffer2,buffer2s
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: fullempty
  TYPE(FFTHyb) :: FFTmrka

  IF ( .NOT. op%done ) &
    CALL ERROR("Ctqmcoffdiag_getResult : Simulation not run                ")

  flavors     =  op%flavors
  inv_flavors = 1.d0 / DBLE(flavors)


  inv_size = 1.d0 / DBLE(op%size)
  sp1 = 0
  spAll = 0

!#ifdef CTCtqmcoffdiag_CHECK
  IF ( op%opt_check .GT. 0 ) THEN
    op%errorImpurity = ImpurityOperator_getError(op%Impurity) * inv_flavors 
    op%errorBath     = BathOperatoroffdiag_getError    (op%Bath    ) * inv_flavors 
  END IF
!#endif

  MALLOC(alpha,(1,1))
  MALLOC(beta,(1,1))
  MALLOC(buffer,(1,1))
  IF ( op%opt_noise .EQ. 1) THEN
    FREEIF(alpha)
    MALLOC(alpha,(1:op%samples+1,1:flavors))
    FREEIF(beta)
    MALLOC(beta,(1:op%samples+1,1:flavors))
  END IF

  IF ( op%have_MPI .EQV. .TRUE.) THEN 
    sp1   = 0
    spAll = sp1 + flavors + 6 

!#ifdef CTCtqmcoffdiag_ANALYSIS
    IF ( op%opt_analysis .EQ. 1 ) &
      spAll = spAll + 3*sp1 
    IF ( op%opt_order .GT. 0 ) &
      spAll = spAll + op%opt_order 
    IF ( op%opt_noise .EQ. 1 ) &
      spAll = spAll + 2*(op%samples + 1)
!#endif

    FREEIF(buffer)
    MALLOC(buffer,(1:spAll,1:MAX(2,flavors)))
  END IF

!  op%seg_added    = op%seg_added    * inv_flavors 
!  op%seg_removed  = op%seg_removed  * inv_flavors
!  op%seg_sign     = op%seg_sign     * inv_flavors
!  op%anti_added   = op%anti_added   * inv_flavors
!  op%anti_removed = op%anti_removed * inv_flavors
!  op%anti_sign    = op%anti_sign    * inv_flavors
  op%stats(:) = op%stats(:) * inv_flavors

  DO iflavor = 1, flavors
    ! Accumulate last values of  N (see also ctqmc_measn)
    op%measN(1,iflavor) = op%measN(1,iflavor) + op%measN(3,iflavor)*op%measN(4,iflavor)
    op%measN(2,iflavor) = op%measN(2,iflavor) + op%measN(4,iflavor)
    ! Reduction
    op%measN(1,iflavor)  = op%measN(1,iflavor) / ( op%measN(2,iflavor) * op%beta )
    ! Correction
!#ifdef CTCtqmcoffdiag_ANALYSIS
    IF ( op%opt_order .GT. 0 ) &
      op%measPerturbation(:   ,iflavor) = op%measPerturbation(:,iflavor) &
                                    / SUM(op%measPerturbation(:,iflavor))
    IF ( op%opt_order .GT. 0 ) &
      op%meas_fullemptylines(:   ,iflavor) = op%meas_fullemptylines(:,iflavor) &
                                    / SUM(op%meas_fullemptylines(:,iflavor))
    !write(6,*) "sum fullempty",iflavor,op%meas_fullemptylines(:,iflavor)

    IF ( op%opt_analysis .EQ. 1 ) THEN
      op%measCorrelation (:,1,iflavor) = op%measCorrelation  (:,1,iflavor) &
                                    / SUM(op%measCorrelation (:,1,iflavor)) &
                                    * op%inv_dt 
      op%measCorrelation (:,2,iflavor) = op%measCorrelation  (:,2,iflavor) &
                                    / SUM(op%measCorrelation (:,2,iflavor)) &
                                    * op%inv_dt 
      op%measCorrelation (:,3,iflavor) = op%measCorrelation  (:,3,iflavor) &
                                    / SUM(op%measCorrelation (:,3,iflavor)) &
                                    * op%inv_dt 
    END IF
!#endif
    IF ( op%opt_noise .EQ. 1 ) THEN
      TabX(1) = DBLE(op%modNoise2)
      TabX(2) = DBLE(op%modNoise1)
      DO itau = 1, op%samples+1
        op%measNoiseG(itau,iflavor,2)%vec = -op%measNoiseG(itau,iflavor,2)%vec*op%inv_dt &  
                                           /(op%beta*DBLE(op%modNoise2))
        op%measNoiseG(itau,iflavor,1)%vec = -op%measNoiseG(itau,iflavor,1)%vec*op%inv_dt &  
                                           /(op%beta*DBLE(op%modNoise1))
        n2 = op%measNoiseG(itau,iflavor,2)%tail
        TabY(1) = Stat_deviation(op%measNoiseG(itau,iflavor,2)%vec(1:n2))!*SQRT(n2/(n2-1))
        n1 = op%measNoiseG(itau,iflavor,1)%tail
        TabY(2) = Stat_deviation(op%measNoiseG(itau,iflavor,1)%vec(1:n1))!*SQRT(n1/(n1-1))
        CALL Stat_powerReg(TabX,SQRT(2.d0*LOG(2.d0))*TabY,alpha(itau,iflavor),beta(itau,iflavor),r)
        ! ecart type -> 60%
        ! largeur a mi-hauteur d'une gaussienne -> sqrt(2*ln(2))*sigma
      END DO
    END IF

  END DO
!sui!write(6,*) "getresults"
  CALL GreenHyboffdiag_measHybrid(op%Greens, op%Bath%M, op%Impurity%Particles, .TRUE.,op%signvalue)
  CALL GreenHyboffdiag_getHybrid(op%Greens)
 ! write(6,*) "op%measN",op%measN(1,:)
  MALLOC(measN_1,(flavors))
  do iflavor=1,flavors
    measN_1(iflavor)=op%measN(1,iflavor)
  enddo
  CALL GreenHyboffdiag_setN(op%Greens, measN_1(:))
  FREE(measN_1)

! todoab case _nd and _d are not completely described.
  FREEIF(buffer2)
  FREEIF(buffer2s)
  sizeoper=size(op%Greens%oper,1)
  !write(6,*) "sss",size(op%Greens%oper,1),sizeoper
  !write(6,*) "sss",size(op%Greens%oper,2),flavors
  !write(6,*) "sss",size(op%Greens%oper,3),flavors
  MALLOC(buffer2,(1:sizeoper,flavors,flavors))
  MALLOC(buffer2s,(1:sizeoper,flavors,flavors))
  MALLOC(fullempty,(2,flavors))
      !sui!write(6,*) "greens1"
  IF ( op%have_MPI .EQV. .TRUE. ) THEN 
      !sui!write(6,*) "greens2"
    fullempty=0.d0
    buffer2 = op%Greens%oper
    !write(6,*) "buffer2",(op%Greens%oper(1,n1,n1),n1=1,flavors)
    buffer2s= 0.d0
    do iflavor=1,flavors
      do itau=1,sizeoper
    !sui!write(6,*) "greens",iflavor,itau,op%Greens%oper(itau,iflavor,iflavor)
      enddo
    enddo
   !write(6,*) "beforempi",op%Greens%oper(1,1,1) ,buffer2(1,1,1)
#ifdef HAVE_MPI
   CALL MPI_COMM_SIZE(op%MY_COMM,nbprocs,ierr)
   CALL MPI_COMM_RANK(op%MY_COMM,myrank,ierr)
#endif
  !write(6,*) "procs",nbprocs,myrank
  END IF
  last = sp1

  op%measDE(:,:) = op%measDE(:,:) * DBLE(op%measurements) /(DBLE(op%sweeps)*op%beta)

  n1 = op%measNoise(1)%tail
  n2 = op%measNoise(2)%tail

  ! On utilise freqs comme tableau de regroupement
  ! Gather de Noise1
  IF ( op%have_MPI .EQV. .TRUE. ) THEN
    MALLOC(counts,(1:op%size))
    MALLOC(displs,(1:op%size))
    FREEIF(freqs)
    MALLOC(freqs,(1:op%size*n1))
    freqs = 0.d0
    freqs(n1*op%rank+1:n1*(op%rank+1)) = op%measNoise(1)%vec(1:n1) 
    counts(:) = n1
    displs(:) = (/ ( iflavor*n1, iflavor=0, op%size-1 ) /)
#ifdef HAVE_MPI
    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DOUBLE_PRECISION, &
                        freqs, counts, displs, &
                        MPI_DOUBLE_PRECISION, op%MY_COMM, ierr)
#endif
    n1 = op%size*n1
    CALL Vector_setSize(op%measNoise(1),n1)
    op%measNoise(1)%vec(1:n1) = freqs(:)
    ! Gather de Noise2
    FREE(freqs)
    MALLOC(freqs,(1:op%size*n2))
    freqs = 0.d0
    freqs(n2*op%rank+1:n2*(op%rank+1)) = op%measNoise(2)%vec(1:n2) 
    counts(:) = n2
    displs(:) = (/ ( iflavor*n2, iflavor=0, op%size-1 ) /)
#ifdef HAVE_MPI
    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        freqs, counts, displs, &
                        MPI_DOUBLE_PRECISION, op%MY_COMM, ierr)
#endif
    n2 = op%size*n2
    CALL Vector_setSize(op%measNoise(2),n2)
    op%measNoise(2)%vec(1:n2) = freqs(:)
    FREE(counts)
    FREE(displs)
    FREE(freqs)
  END IF
  !n1 = op%measNoise(1)%tail
  !n2 = op%measNoise(2)%tail

  ! Transformation des paquets pour que ca fit a CTQMC_SLICE(1|2)
  IF ( n1 .GT. CTQMC_SLICE1 ) THEN
    itau = n1/CTQMC_SLICE1
    MALLOC(freqs,(1:n1/itau))
    DO debut=1, n1/itau
      freqs(debut)=SUM(op%measNoise(1)%vec((debut-1)*itau+1:itau*debut))
    END DO
    freqs(:) = freqs(:)/DBLE(itau)
    op%modNoise1 = op%modNoise1*itau
    n1 = n1/itau
    CALL Vector_setSize(op%measNoise(1),n1)
    op%measNoise(1)%vec(1:n1) = freqs(:)
    FREE(freqs)
  END IF
  IF ( n2 .GT. CTQMC_SLICE1*CTQMC_SLICE2 ) THEN
    itau = n2/(CTQMC_SLICE1*CTQMC_SLICE2)
    MALLOC(freqs,(1:n2/itau))
    DO debut=1, n2/itau
      freqs(debut)=SUM(op%measNoise(2)%vec((debut-1)*itau+1:itau*debut))
    END DO
    freqs(:) = freqs(:)/DBLE(itau)
    op%modNoise2 = op%modNoise2*itau
    n2 = n2/itau
    CALL Vector_setSize(op%measNoise(2),n2)
    op%measNoise(2)%vec(1:n2) = freqs(:)
    FREE(freqs)
  END IF
  ! On peut s'amuser avec nos valeur d'energies
  !MALLOC(TabX,(1:20))
  !MALLOC(TabY,(1:20))

  TabX(1) = DBLE(op%modNoise2)
  TabX(2) = DBLE(op%modNoise1)

  ! Il faut calculer pour chaque modulo 10 ecarts type sur les donnes acquises
  op%measNoise(1)%vec(1:n1) = op%measNoise(1)%vec(1:n1)/(op%beta*DBLE(op%modNoise1))*DBLE(op%measurements)
  op%measNoise(2)%vec(1:n2) = op%measNoise(2)%vec(1:n2)/(op%beta*DBLE(op%modNoise2))*DBLE(op%measurements)
!  CALL Vector_print(op%measNoise(1),op%rank+70)
!  CALL Vector_print(op%measNoise(2),op%rank+50)
!  DO iflavor=1,10
!    debut = (iflavor-1)*n2/10+1
!    fin   = iflavor*n2/10
!    TabY(iflavor) = Stat_deviation(op%measNoise(2)%vec(debut:fin))
!    debut = (iflavor-1)*n1/10+1
!    fin   = iflavor*n1/10
!    TabY(10+iflavor) = Stat_deviation(op%measNoise(1)%vec(debut:fin))
!  END DO
!!  TabY(1:n) = (op%measNoise(2)%vec(1:n)   &
!!              )
!!             !/(op%beta*DBLE(op%modNoise2))*DBLE(op%measurements) &
!!             !- op%measDE(1,1))
!!  TabY(op%measNoise(2)%tail+1:n+op%measNoise(2)%tail) = (op%measNoise(1)%vec(1:n)   &
!!               )
!!             ! /(op%beta*DBLE(op%modNoise1))*DBLE(op%measurements) &
!!             ! - op%measDE(1,1))
!  IF ( op%rank .EQ. 0 ) THEN
!    DO iflavor=1,20
!      write(45,*) TabX(iflavor), TabY(iflavor)
!    END DO
!  END IF
!


  TabY(1) = Stat_deviation(op%measNoise(2)%vec(1:n2))!*SQRT(n2/(n2-1))
!!  write(op%rank+10,*) TabX(2)
!!  write(op%rank+40,*) TabX(1)
!!  CALL Vector_print(op%measNoise(1),op%rank+10)
!!  CALL Vector_print(op%measNoise(2),op%rank+40)
!!  CLOSE(op%rank+10)
!!  CLOSE(op%rank+40)
  TabY(2) = Stat_deviation(op%measNoise(1)%vec(1:n1))!*SQRT(n1/(n1-1))
!!  ! Ecart carre moyen ~ ecart type mais non biaise. Serait moins precis. Aucun
  ! impact sur la pente, juste sur l'ordonnee a l'origine.

  CALL Stat_powerReg(TabX,SQRT(2.d0*LOG(2.d0))*TabY,a,b,r)
!  FREE(TabX)
!  FREE(TabY)
  ! ecart type -> 60%
  ! largeur a mi-hauteur d'une gaussienne -> sqrt(2*ln(2))*sigma

  !op%measDE(1,1) = SUM(op%measNoise(1)%vec(1:op%measNoise(1)%tail))/(DBLE(op%measNoise(1)%tail*op%modNoise1)*op%beta)
  !op%measDE(2:flavors,1:flavors) = op%measDE(2:flavors,1:flavors) /(DBLE(op%sweeps)*op%beta)
  CALL ImpurityOperator_getErrorOverlap(op%Impurity,op%measDE)
  ! Add the difference between true calculation and quick calculation of the
  ! last sweep overlap to measDE(2,2)
  !op%measDE = op%measDE * DBLE(op%measurements) 
  IF ( op%have_MPI .EQV. .TRUE. ) THEN 
    IF ( op%opt_analysis .EQ. 1 ) THEN
      buffer(last+1:last+sp1,:) = op%measCorrelation(:,1,:)
      last = last + sp1
      buffer(last+1:last+sp1,:) = op%measCorrelation(:,2,:)
      last = last + sp1
      buffer(last+1:last+sp1,:) = op%measCorrelation(:,3,:)
      last = last + sp1
    END IF
    IF ( op%opt_order .GT. 0 ) THEN
      buffer(last+1:last+op%opt_order, :) = op%measPerturbation(:,:)
      last = last + op%opt_order
    END IF
    IF ( op%opt_noise .EQ. 1 ) THEN
      buffer(last+1:last+op%samples+1,:) = alpha(:,:)
      last = last + op%samples + 1
      buffer(last+1:last+op%samples+1,:) = beta(:,:)
      last = last + op%samples + 1
    END IF
!  op%measDE(2,2) = a*EXP(b*LOG(DBLE(op%sweeps*op%size)))
    buffer(spall-(flavors+5):spAll-6,:) = op%measDE(:,:)
!    buffer(spAll  ,1) = op%seg_added   
!    buffer(spAll-1,1) = op%seg_removed 
!    buffer(spAll-2,1) = op%seg_sign    
!    buffer(spAll  ,2) = op%anti_added  
!    buffer(spAll-1,2) = op%anti_removed
!    buffer(spAll-2,2) = op%anti_sign   
    buffer(spAll  ,1) = op%stats(1)
    buffer(spAll-1,1) = op%stats(2)
    buffer(spAll-2,1) = op%stats(3)
    buffer(spAll  ,2) = op%stats(4)
    buffer(spAll-1,2) = op%stats(5)
    buffer(spAll-2,2) = op%stats(6)
    buffer(spAll-3,1) = op%swap
    buffer(spAll-3,2) = DBLE(op%modGlobalMove(2))
    buffer(spAll-4,1) = a
    buffer(spAll-4,2) = b
!#ifdef CTCtqmcoffdiag_CHECK
    buffer(spAll-5,1) = op%errorImpurity
    buffer(spAll-5,2) = op%errorBath 
    signvaluemeassum = 0
!#endif

#ifdef HAVE_MPI
   !write(6,*) "bufferbefore",buffer(1,1)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, buffer, spAll*flavors, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, op%MY_COMM, ierr)
   !write(6,*) "bufferafter",buffer(1,1)
   ! CALL MPI_ALLREDUCE(MPI_IN_PLACE, buffer2, sp1*flavors*flavors, &
   !                  MPI_DOUBLE_PRECISION, MPI_SUM, op%MY_COMM, ierr)
    CALL MPI_ALLREDUCE( buffer2, buffer2s, sizeoper*flavors*flavors, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, op%MY_COMM, ierr)
   !write(6,*) "justaftermpi",op%Greens%oper(1,1,1) ,buffer2s(1,1,1)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, op%runTime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             op%MY_COMM, ierr)
    CALL MPI_ALLREDUCE(op%Greens%signvaluemeas, signvaluemeassum , 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
             op%MY_COMM, ierr)
    IF ( op%opt_order .GT. 0 ) THEN
      CALL MPI_ALLREDUCE(op%meas_fullemptylines, fullempty, 2*flavors, MPI_DOUBLE_PRECISION, MPI_SUM, &
               op%MY_COMM, ierr)
    ENDIF
#endif

  
    buffer          = buffer * inv_size
    op%measDE(:,:)  = buffer(spall-(flavors+5):spAll-6,:)
!    op%seg_added    = buffer(spAll  ,1)
!    op%seg_removed  = buffer(spAll-1,1)
!    op%seg_sign     = buffer(spAll-2,1)
!    op%anti_added   = buffer(spAll  ,2)
!    op%anti_removed = buffer(spAll-1,2)
!    op%anti_sign    = buffer(spAll-2,2)
    op%stats(1)    = buffer(spAll  ,1)
    op%stats(2)    = buffer(spAll-1,1)
    op%stats(3)    = buffer(spAll-2,1)
    op%stats(4)    = buffer(spAll  ,2)
    op%stats(5)    = buffer(spAll-1,2)
    op%stats(6)    = buffer(spAll-2,2)
    op%swap         = buffer(spAll-3,1)
    op%modGlobalMove(2) = NINT(buffer(spAll-3,2))
    a               = buffer(spAll-4,1) 
    b               = buffer(spAll-4,2)
!!#ifdef CTCtqmcoffdiag_CHECK
    op%errorImpurity= buffer(spAll-5,1) 
    op%errorBath    = buffer(spAll-5,2)   
!#endif

   ! DO iflavor = 1, flavors
   !   op%Greens(iflavor)%oper          = buffer(1:sp1          , iflavor)
   ! END DO
    op%Greens%oper = buffer2s/float(nbprocs)
   ! write(6,*) "buffer2s",(op%Greens%oper(1,n1,n1),n1=1,flavors)
    op%Greens%signvaluemeas = signvaluemeassum/float(nbprocs)
    !sui!write(6,*) "nbprocs",nbprocs,op%Greens%signvaluemeas
    op%Greens%oper = op%Greens%oper / op%Greens%signvaluemeas
   ! write(6,*) "buffer3s",(op%Greens%oper(1,n1,n1),n1=1,flavors)
    IF ( op%opt_order .GT. 0 ) THEN
      op%meas_fullemptylines= fullempty/float(nbprocs)
    ENDIF
    do iflavor=1,flavors
      do itau=1,sizeoper
    !sui!write(6,*) "greens_av",iflavor,itau,op%Greens%oper(itau,iflavor,iflavor)
      enddo
    enddo
   !write(6,*) "aftermpi",op%Greens%oper(1,1,1) ,buffer2s(1,1,1)
    last = sp1
    IF ( op%opt_analysis .EQ. 1 ) THEN
      op%measCorrelation(:,1,:) = buffer(last+1:last+sp1,:) 
      last = last + sp1
      op%measCorrelation(:,2,:) = buffer(last+1:last+sp1,:) 
      last = last + sp1
      op%measCorrelation(:,3,:) = buffer(last+1:last+sp1,:) 
      last = last + sp1
    END IF
    IF ( op%opt_order .GT. 0 ) THEN
      op%measPerturbation(:,:) = buffer(last+1:last+op%opt_order, :)
      last = last + op%opt_order
    END IF
    IF ( op%opt_noise .EQ. 1 ) THEN
      alpha(:,:) = buffer(last+1:last+op%samples+1,:)
      last = last + op%samples + 1
      beta(:,:) = buffer(last+1:last+op%samples+1,:)
      last = last + op%samples + 1
    END IF
  END IF
  DO iflavor = 1, flavors
    ! complete DE matrix
    op%measDE(iflavor, iflavor+1:flavors) = op%measDE(iflavor+1:flavors,iflavor)
  END DO
  FREE(buffer)
  FREE(buffer2)
  FREE(buffer2s)
  FREE(fullempty)

  IF ( op%opt_spectra .GE. 1 ) THEN
    endDensity = SIZE(op%density,2)
    IF ( op%density(1,endDensity) .EQ. -1.d0 ) &
      endDensity = endDensity - 1
    CALL FFTHyb_init(FFTmrka,endDensity,DBLE(op%thermalization)/DBLE(op%measurements*op%opt_spectra))
    ! Not very Beauty 
    MALLOC(freqs,(1:FFTmrka%size/2))
    DO iflavor = 1, flavors
      ! mean value is removed to supress the continue composent 
      CALL FFTHyb_setData(FFTmrka,op%density(iflavor,1:endDensity)/op%beta+op%Greens%oper(op%samples+1,iflavor,iflavor))
      CALL FFTHyb_run(FFTmrka,1)
      CALL FFTHyb_getData(FFTmrka,endDensity,op%density(iflavor,:),freqs)
    END DO
    op%density(flavors+1,:) = -1.d0
    op%density(flavors+1,1:FFTmrka%size/2) = freqs
    CALL FFTHyb_destroy(FFTmrka)
    FREE(freqs)
  END IF

  op%a_Noise = a
  op%b_Noise = b
  IF ( op%opt_noise .EQ. 1 ) THEN
    op%abNoiseG(1,:,:) = alpha
    op%abNoiseG(2,:,:) = beta
  END IF
  FREE(alpha)
  FREE(beta)

END SUBROUTINE Ctqmcoffdiag_getResult
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_symmetrizeGreen
!! NAME
!!  Ctqmcoffdiag_symmetrizeGreen
!!
!! FUNCTION
!!  optionnaly symmetrize the green functions
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_symmetrizeGreen(op, syms)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)                     , INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN   ) :: syms
!Local variables ------------------------------
  !INTEGER :: iflavor1
  !INTEGER :: iflavor2
  !INTEGER :: flavors
  !DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: green_tmp
  !DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:  ) :: n_tmp

  ABI_UNUSED((/syms(1,1), op%swap/))

!  flavors = op%flavors
!  IF ( SIZE(syms,1) .NE. flavors .OR. SIZE(syms,2) .NE. flavors ) THEN
!    CALL WARNALL("Ctqmcoffdiag_symmetrizeGreen : wrong opt_sym -> not symmetrizing")
!    RETURN
!  END IF
! 
!  MALLOC(green_tmp,(1:op%samples+1,flavors))
!  green_tmp(:,:) = 0.d0
!  MALLOC(n_tmp,(1:flavors))
!  n_tmp(:) = 0.d0
!  DO iflavor1=1, flavors
!    DO iflavor2=1,flavors
!      green_tmp(:,iflavor1) = green_tmp(:,iflavor1) &
!                             + syms(iflavor2,iflavor1) * op%Greens(iflavor2)%oper(:)
!      n_tmp(iflavor1) = n_tmp(iflavor1) &
!                             + syms(iflavor2,iflavor1) * op%measN(1,iflavor2)
!    END DO
!  END DO
!  DO iflavor1=1, flavors
!    op%Greens(iflavor1)%oper(:) = green_tmp(:,iflavor1)
!    op%measN(1,iflavor1)          = n_tmp(iflavor1)
!  END DO
!  FREE(green_tmp)
!  FREE(n_tmp)
END SUBROUTINE Ctqmcoffdiag_symmetrizeGreen
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_getGreen
!! NAME
!!  Ctqmcoffdiag_getGreen
!!
!! FUNCTION
!!  Get the full green functions in time and/or frequency
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_getGreen(op, Gtau, Gw)

!Arguments ------------------------------------
 USE m_GreenHyboffdiag
  TYPE(Ctqmcoffdiag)          , INTENT(INOUT)    :: op
  DOUBLE PRECISION, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: Gtau
  COMPLEX(KIND=8), DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: Gw
!Local variables ------------------------------
  !INTEGER                            :: itime
  INTEGER                            :: iflavor1
  INTEGER                            :: iflavor1b !,iflavor,iflavorbis
  INTEGER                            :: iflavor2
  INTEGER                            :: iflavor3
  INTEGER                            :: flavors,tail
  INTEGER                            :: ifreq,itime
  DOUBLE PRECISION :: u1 
  DOUBLE PRECISION :: u2
  DOUBLE PRECISION :: u3
  DOUBLE PRECISION :: Un
  DOUBLE PRECISION :: UUnn,iw !omega,
  CHARACTER(LEN=4)                   :: cflavors
  CHARACTER(LEN=50)                  :: string
  TYPE(GreenHyboffdiag)                     :: F_tmp

  flavors = op%flavors
  DO iflavor1 = 1, flavors
    u1 = 0.d0
    u2 = 0.d0
    u3 = 0.d0
    DO iflavor2 = 1, flavors
      IF ( iflavor2 .EQ. iflavor1 ) CYCLE
      Un = op%Impurity%mat_U(iflavor2,iflavor1) * op%measN(1,iflavor2)
!      Un = op%Impurity%mat_U(iflavor2,iflavor1) * (op%Greens%oper(1,iflavor2,iflavor2) + 1.d0)
      !write(6,*) "forsetmoments",iflavor1,iflavor2,(op%Greens%oper(1,iflavor2,iflavor2) + 1.d0), Un
      u1 = u1 + Un 
      u2 = u2 + Un*op%Impurity%mat_U(iflavor2,iflavor1) 
      DO iflavor3 = 1, flavors
        IF ( iflavor3 .EQ. iflavor2 .OR. iflavor3 .EQ. iflavor1 ) CYCLE
        UUnn = (op%Impurity%mat_U(iflavor2,iflavor1)*op%Impurity%mat_U(iflavor3,iflavor1)) * &
&                                                    op%measDE(iflavor2,iflavor3) 
        u2 = u2 + UUnn 
      END DO
    END DO  
     ! write(6,*) "u1,u2",u1,u2

    DO iflavor1b = 1, flavors
      u3 =-(op%Impurity%mat_U(iflavor1,iflavor1b))*op%Greens%oper(1,iflavor1,iflavor1b)
      ! u3=U_{1,1b}*G_{1,1b}
      CALL GreenHyboffdiag_setMoments(op%Greens,iflavor1,iflavor1b,u1,u2,u3)
    END DO ! iflavor1b

  END DO ! iflavor1

  IF ( PRESENT( Gtau ) ) THEN
    DO iflavor1 = 1, flavors
      DO iflavor2 = 1, flavors
        Gtau(1:op%samples,iflavor1,iflavor2) = op%Greens%oper(1:op%samples,iflavor1,iflavor2)
      END DO  
    END DO ! iflavor1
  END IF
! !--------- Write Occupation matrix before Gtau
!  write(ostream,'(17x,a)') "Occupation matrix"
!  write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
!  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
!  do iflavor=1, op%flavors
!    write(ostream,'(7x,i10,a,30f10.4)') iflavor,"|",(-op%Greens%oper(op%samples,iflavor,iflavorbis),iflavorbis=1,op%flavors)
!  enddo
!  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors),ch10
! !------------------------------------------------------------------------------------------
! !--------- Write Occupation matrix Gtau
!  write(ostream,'(17x,a)') "Occupation matrix"
!  write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
!  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
!  do iflavor=1, op%flavors
!    write(ostream,'(7x,i10,a,30f10.4)') iflavor,"|",(Gtau(op%samples,iflavor,iflavorbis),iflavorbis=1,op%flavors)
!  enddo
!  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors),ch10
! !------------------------------------------------------------------------------------------

!================================================
  if(3==4) then
!================================================
    DO iflavor1 = 1, flavors
      DO iflavor1b = 1, flavors
       !call nfourier3(op%Greens%oper(1:op%samples,iflavor1,iflavor1b),Gw(1:op%samples,iflavor1,iflavor1b),iflavor1==iflavor1b,op%Greens%samples,op%Greens%samples-1,op%Greens%beta,1.d0,op%Greens%Mk(iflavor1,iflavor1b,1),op%Greens%Mk(iflavor1,iflavor1b,2),op%Greens%Mk(iflavor1,iflavor1b,3))
      END DO  
    END DO ! iflavor1
!    ============== write Gomega_nd.dat
    if(op%rank==0) then
    OPEN(UNIT=44, FILE="Gomega_nd_nfourier2.dat")
    WRITE(cflavors,'(I4)') 2*(flavors*flavors+1)
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'E15.5)'
    !write(6,*) " op%Greens%Wmax", op%Greens%Wmax
    do  iflavor1=1, flavors 
      do  iflavor1b=1, flavors 
    write(44,*) "#op%Greens%Mk(iflavor1,iflavor2,1",op%Greens%Mk(iflavor1,iflavor1b,:)
        DO ifreq = 1, op%samples
!      !write(6,string) (DBLE(ifreq)*2-1)*3.1415/op%Greens%beta, &
!      (/ ((real(Gw(ifreq,iflavor1,iflavor1b)),imag(Gw(ifreq,iflavor1,iflavor1b)), iflavor1=1, flavors),iflavor1b=1,flavors) /)
!      WRITE(44,string) (DBLE(ifreq)*2.d0-1.d0)*3.1415926/op%Greens%beta, &
        iw=aimag(Gw(ifreq,op%flavors,op%flavors+1))
         WRITE(44,string) aimag(Gw(ifreq,op%flavors,op%flavors+1)),&
         real(Gw(ifreq,iflavor1,iflavor1b)),aimag(Gw(ifreq,iflavor1,iflavor1b)),&
            ( -op%Greens%Mk(iflavor1,iflavor1b,2) )/(iw*iw) , (op%Greens%Mk(iflavor1,iflavor1b,1))/iw!-op%Greens%Mk(iflavor1,iflavor1b,3)/(iw*iw))/iw 
      !   WRITE(102,*) aimag(Gw(ifreq,op%flavors,op%flavors+1)), (op%Greens%Mk(iflavor1,iflavor1b,1))/iw,op%Greens%Mk(iflavor1,iflavor1b,1),iw
        END DO
         WRITE(44,*) 
      END DO
    END DO
    close(44)
    endif
!================================================
  endif
!================================================
       !!write(6,*) "present gw", present(gw)
  IF ( PRESENT( Gw ) ) THEN
     !!write(6,*) "size gw",SIZE(Gw,DIM=2) ,flavors+1 
    IF ( SIZE(Gw,DIM=3) .EQ. flavors+1 ) THEN
     ! CALL GreenHyboffdiag_forFourier(op%Greens, Gomega=Gw, omega=Gw(:,op%flavors,op%flavors+1))
      CALL GreenHyboffdiag_forFourier(op%Greens, Gomega=Gw, omega=Gw(:,op%flavors,op%flavors+1))
      !write(6,*) "1"
      !IF ( op%rank .EQ. 0 ) write(20,*) Gw(:,iflavor1)
    ELSE IF ( SIZE(Gw,DIM=3) .EQ. flavors ) THEN  
      CALL GreenHyboffdiag_forFourier(op%Greens,Gomega=Gw)
      !write(6,*) "2"
    ELSE
      CALL WARNALL("Ctqmcoffdiag_getGreen : Gw is not valid                    ")
      CALL GreenHyboffdiag_forFourier(op%Greens,Wmax=op%Wmax)
      !write(6,*) "3"
    END IF
  ELSE
    CALL GreenHyboffdiag_forFourier(op%Greens,Wmax=op%Wmax)
  END IF
!  ============== write Gomega_nd.dat
!================================================
!  if(3==4) then
!================================================
  if(op%rank==0.and.3==4) then
  OPEN(UNIT=44, FILE="Gomega_nd.dat")
  WRITE(cflavors,'(I4)') 2*(flavors*flavors+1)
  string = '(1x,'//TRIM(ADJUSTL(cflavors))//'E15.5)'
  !write(6,*) " op%Greens%Wmax", op%Greens%Wmax
  do  iflavor1=1, flavors 
    do  iflavor1b=1, flavors 
  write(44,*) "#op%Greens%Mk(iflavor1,iflavor2,1",op%Greens%Mk(iflavor1,iflavor1b,:)
      DO ifreq = 1, SIZE(Gw,1)    
!    !write(6,string) (DBLE(ifreq)*2-1)*3.1415/op%Greens%beta, &
!    (/ ((real(Gw(ifreq,iflavor1,iflavor1b)),imag(Gw(ifreq,iflavor1,iflavor1b)), iflavor1=1, flavors),iflavor1b=1,flavors) /)
!    WRITE(44,string) (DBLE(ifreq)*2.d0-1.d0)*3.1415926/op%Greens%beta, &
      iw=aimag(Gw(ifreq,op%flavors,op%flavors+1))
       WRITE(44,string) aimag(Gw(ifreq,op%flavors,op%flavors+1)),&
       real(Gw(ifreq,iflavor1,iflavor1b)),aimag(Gw(ifreq,iflavor1,iflavor1b)),&
          ( -op%Greens%Mk(iflavor1,iflavor1b,2) )/(iw*iw) , &
&          (op%Greens%Mk(iflavor1,iflavor1b,1)-op%Greens%Mk(iflavor1,iflavor1b,3)/(iw*iw))/iw 
      END DO
       WRITE(44,*) 
    END DO
  END DO
  endif
!================================================
!  endif
!================================================


!  ==============================
  ! --- Initialize F_tmp 
  !write(6,*) "10"

  IF ( op%have_MPI .EQV. .TRUE. ) THEN
    !CALL GreenHyboffdiag_init(F_tmp,op%samples,op%beta,op%flavors,MY_COMM=op%MY_COMM)
    CALL GreenHyboffdiag_init(F_tmp,op%samples,op%beta,flavors)
    !write(6,*) "10a"
  ELSE
    CALL GreenHyboffdiag_init(F_tmp,op%samples,op%beta,flavors)
    !write(6,*) "10b"
  END IF

  !write(6,*) "11"
!  CALL GreenHyboffdiag_setOperW(F_tmp,Gw)

  tail = op%samples 
  F_tmp%Wmax=op%samples ! backFourier only works for linear freq: calculation of A and etc..
  MALLOC(F_tmp%oper_w,(1:tail,op%flavors,op%flavors))
  F_tmp%oper_w(1:tail,1:F_tmp%nflavors,1:F_tmp%nflavors) = Gw(1:tail,1:F_tmp%nflavors,1:F_tmp%nflavors)
  !write(6,*) "example",F_tmp%oper_w(1,1,1)
  !write(6,*) "example",Gw(1,1,1)
  F_tmp%setW = .TRUE.
  !write(6,*) size(F_tmp%oper_w,1)
  !write(6,*) size(F_tmp%oper_w,2)
  !write(6,*) size(F_tmp%oper_w,3)
  !write(6,*) size(Gw,1)
  !write(6,*) size(Gw,2)
  !write(6,*) size(Gw,3)

  !write(6,*) "eee", (2.d0*DBLE(ifreq)-1.d0) * 3.1415/op%beta,real(F_tmp%oper_w(1,1,1)),imag(F_tmp%oper_w(1,1,1))
!================================================
  if(3==4) then
!================================================
    OPEN(UNIT=3337, FILE="Gomega_nd2.dat")
    do  iflavor1=1, flavors 
      do  iflavor1b=1, flavors 
        do  ifreq=1, tail
!         write(3337,*) (2.d0*DBLE(ifreq)-1.d0) * 3.1415/op%beta,real(F_tmp%oper_w(ifreq,iflavor1,iflavor1b)),&
!   &     imag(F_tmp%oper_w(ifreq,iflavor1,iflavor1b))
          write(3337,*) aimag(Gw(ifreq,op%flavors,op%flavors+1)), real(F_tmp%oper_w(ifreq,iflavor1,iflavor1b)),&
 &        aimag(F_tmp%oper_w(ifreq,iflavor1,iflavor1b))

      !    omega=(2.d0*DBLE(ifreq)-1.d0) * 3.1415/op%beta
 !        F_tmp%oper_w(ifreq,iflavor1,iflavor1b)=0.1**2/Gw(ifreq,op%flavors,op%flavors+1)
        enddo 
        write(3337,*)
      enddo 
    enddo
    close(3337)
!================================================
  endif
!================================================

  !write(6,*) "12",F_tmp%Wmax

!  CALL GreenHyboffdiag_backFourier(F_tmp,func="green")

  !write(6,*) "13"

!================================================
  if(3==4) then
!================================================
    OPEN(UNIT=48, FILE="Gtau_nd_2.dat")
!    --- Print full non diagonal Gtau in Gtau_nd.dat
    WRITE(cflavors,'(I4)') flavors*flavors+1
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'ES22.14)'
    DO itime = 1, op%samples+1
      WRITE(48,string) DBLE(itime-1)*op%beta/DBLE(op%samples), &
 &    ((F_tmp%oper(itime,iflavor1,iflavor1b), iflavor1=1, flavors),iflavor1b=1,flavors)
    END DO
!================================================
  endif
!================================================

  CALL GreenHyboffdiag_destroy(F_tmp)

  !FREE(F_tmp%oper_w)
!  ==============================
END SUBROUTINE Ctqmcoffdiag_getGreen
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_getD
!! NAME
!!  Ctqmcoffdiag_getD
!!
!! FUNCTION
!!  get double occupation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_getD(op, D)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)       , INTENT(IN ) :: op
  DOUBLE PRECISION, INTENT(OUT) :: D
!Local variables ------------------------------
  INTEGER                       :: iflavor1
  INTEGER                       :: iflavor2

  D = 0.d0

  DO iflavor1 = 1, op%flavors
    DO iflavor2 = iflavor1+1, op%flavors
      D = D + op%measDE(iflavor2,iflavor1)
    END DO
  END DO
  !IF ( op%rank .EQ. 0 ) THEN
  !  DO iflavor1 = 1, op%flavors
  !    DO iflavor2 = iflavor1+1, op%flavors
  !     write(4533,*) op%measDE(iflavor2,iflavor1)k
  !     write(4534,*) op%Impurity%mat_U(iflavor2,iflavor1)k
  !    END DO
  !  END DO

  !ENDIF

END SUBROUTINE Ctqmcoffdiag_getD
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_getE
!! NAME
!!  Ctqmcoffdiag_getE
!!
!! FUNCTION
!!  get interaction energy and noise on it
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_getE(op,E,noise)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)       , INTENT(IN ) :: op
  DOUBLE PRECISION, INTENT(OUT) :: E
  DOUBLE PRECISION, INTENT(OUT) :: Noise

  E = op%measDE(1,1)  
  Noise = op%a_Noise*(DBLE(op%sweeps)*DBLE(op%size))**op%b_Noise
END SUBROUTINE Ctqmcoffdiag_getE
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printAll
!! NAME
!!  Ctqmcoffdiag_printAll
!!
!! FUNCTION
!!  print different functions computed during the simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printAll(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op

  IF ( .NOT. op%done ) &
    CALL WARNALL("Ctqmcoffdiag_printAll : Simulation not run                 ")

!sui!write(6,*) "op%stats",op%stats
  CALL Ctqmcoffdiag_printQMC(op)

  CALL Ctqmcoffdiag_printGreen(op)

  CALL Ctqmcoffdiag_printD(op)

!  CALL Ctqmcoffdiag_printE(op)

!#ifdef CTCtqmcoffdiag_ANALYSIS
  CALL Ctqmcoffdiag_printPerturbation(op)

  CALL Ctqmcoffdiag_printCorrelation(op)
!#endif

  CALL Ctqmcoffdiag_printSpectra(op)

END SUBROUTINE Ctqmcoffdiag_printAll
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printQMC
!! NAME
!!  Ctqmcoffdiag_printQMC
!!
!! FUNCTION
!!  print ctqmc statistics
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printQMC(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  INTEGER                  :: ostream
  INTEGER                  :: iflavor,iflavorbis,iorder
  DOUBLE PRECISION         :: sweeps
  DOUBLE PRECISION         :: invSweeps
  CHARACTER(LEN=2)         :: a
  CHARACTER(LEN=15)        :: string

  !IF ( op%rank .NE. 0) RETURN
  IF ( op%rank .NE. MOD(op%size,op%size)) RETURN

  ostream   = op%ostream
  sweeps    = DBLE(op%sweeps)
  invSweeps = 1.d0/sweeps

  WRITE(ostream,'(1x,F13.0,A11,F10.2,A12,I5,A5)') sweeps*DBLE(op%size), " sweeps in ", op%runTime, &
                 " seconds on ", op%size, " CPUs"
  WRITE(ostream,'(A28,F6.2)') "Segments added        [%] : ", op%stats(4)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Segments removed      [%] : ", op%stats(5)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Segments <0 sign      [%] : ", op%stats(6)*invSweeps*100.d0
  !WRITE(ostream,'(A28,F12.2)') "Number of meas        [%] : ", op%stats(6)
  WRITE(ostream,'(A28,F6.2)') "Anti-segments added   [%] : ", op%stats(1)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Anti-segments removed [%] : ", op%stats(2)*invSweeps*100.d0
  WRITE(ostream,'(A28,F6.2)') "Anti-segments <0 sign [%] : ", op%stats(3)*invSweeps*100.d0
  !WRITE(ostream,'(A28,F12.2)') "Sum of sign       [%] : ", op%stats(3)
  WRITE(ostream,'(A28,F13.2)') "Signe value               : ", op%Greens%signvaluemeas
  IF ( op%modGlobalMove(1) .LT. op%sweeps + 1 ) THEN
    WRITE(ostream,'(A28,F6.2)') "Global Move           [%] : ", op%swap         *invSweeps*100.d0*op%modGlobalMove(1)
    WRITE(ostream,'(A28,F6.2)') "Global Move Reduced   [%] : ", op%swap         / DBLE(op%modGlobalMove(2))*100.d0
  END IF
!#ifdef CTCtqmcoffdiag_CHECK
  IF ( op%opt_check .EQ. 1 .OR. op%opt_check .EQ. 3 ) &
    WRITE(ostream,'(A28,E22.14)') "Impurity test         [%] : ", op%errorImpurity*100.d0
  IF ( op%opt_check .GE. 2 ) &
      WRITE(ostream,'(A28,E22.14)') "Bath     test         [%] : ", op%errorBath    *100.d0
!#endif
  WRITE(ostream,'(A28,ES22.14,A5,ES21.14)') "<Epot>                [U] : ", op%measDE(1,1), " +/- ",&
!#ifdef HAVE_MPI
                                                              op%a_Noise*(sweeps*DBLE(op%size))**op%b_Noise
!#else
!                                                              op%a_Noise*(sweeps)**op%b_Noise
!#endif
 !--------- Write double occupation between all pairs of orbitals --------------------------
  write(ostream,'(17x,a)') "Double occupation between pairs of orbitals"
  write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
  do iflavor=1, op%flavors
    write(ostream,'(7x,i10,a,30f10.4)') iflavor,"|",(op%measDE(iflavor,iflavorbis),iflavorbis=1,op%flavors)
  enddo
  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors),ch10
 !------------------------------------------------------------------------------------------

 !--------- Write number of segments for each orbitals
 ! write(ostream,'(a)') "Number of segments for each orbitals"
 ! write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
 ! write(ostream,'(17x,30a)') ("----------",iflavorbis=1,op%flavors)
 ! do iflavor=1, op%flavors
 !   write(ostream,'(i17,a,30f10.4)') iflavor,"|",(op%Impurity%particles(IT)%tail
 ! enddo
 ! write(ostream,'(17x,30a)') ("----------",iflavorbis=1,op%flavors)
 !------------------------------------------------------------------------------------------
 !--------- Write G(L)
  write(ostream,'(17x,a)') "G(L)"
  write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
  do iflavor=1, op%flavors
    write(ostream,'(7x,i10,a,30f10.4)') iflavor,"|",(op%Greens%oper(op%samples,iflavor,iflavorbis),iflavorbis=1,op%flavors)
  enddo
  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors),ch10
 !------------------------------------------------------------------------------------------
 !--------- Write G(1)
  write(ostream,'(17x,a)') "G(1)"
  write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
  do iflavor=1, op%flavors
    write(ostream,'(7x,i10,a,30f10.4)') iflavor,"|",(op%Greens%oper(1,iflavor,iflavorbis),iflavorbis=1,op%flavors)
  enddo
  write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors),ch10
 !------------------------------------------------------------------------------------------

  WRITE(ostream,'(A28,F8.4,A3,F7.4)') "Noise                 [U] : ", op%a_Noise, " x^", op%b_Noise
  WRITE(ostream,'(A28,E10.2)')  "Niquist puls.     [/beta] : ", ACOS(-1.d0)*op%inv_dt
  WRITE(ostream,'(A28,E22.14)') "Max Acc. Epot Error   [U] : ", op%measDE(2,2)/(op%beta*op%modNoise1*2.d0)*sweeps
  
  !WRITE(ostream,'(A28,F7.4,A3,F7.4,A4,E20.14)') "Noise            [G(tau)] : ", op%a_Noise(2), "x^", op%b_Noise(2), " -> ", &
                                                              !op%a_Noise(2)*(sweeps*DBLE(op%size))**op%b_Noise(2)
 !----- PERTURBATION ORDER------------------------------------------------------------------
  IF ( op%opt_order .GT. 0 ) THEN 
    write(ostream,*) 
    WRITE(a,'(I2)') op%flavors
    string = '(A28,'//TRIM(ADJUSTL(a))//'(1x,I3))'
    WRITE(ostream,string) "Perturbation orders       : ",(/ (MAXLOC(op%measPerturbation(:, iflavor))-1, iflavor=1, op%flavors) /)
    write(ostream,'(17x,a)') "order of Perturbation for flavors"
    write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
    write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
    write(ostream,'(12x,a,30i10)') " max ",(/ (MAXLOC(op%measPerturbation(:, iflavor))-1, iflavor=1, op%flavors) /)
    write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
    do iorder=0, op%opt_order-1
      write(ostream,'(7x,i10,a,30f10.4)') iorder,"|",(op%measPerturbation(iorder+1,iflavor),iflavor=1,op%flavors)
    enddo
  END IF
 !------------------------------------------------------------------------------------------
 !----- PERTURBATION ORDER------------------------------------------------------------------
  IF ( op%opt_order .GT. 0 ) THEN 
    write(ostream,*) 
    write(ostream,'(17x,a)') "Proportion of full and empty orbital for order 0"
    write(ostream,'(17x,30i10)') (iflavorbis,iflavorbis=1,op%flavors)
    write(ostream,'(18x,30a)') ("----------",iflavorbis=1,op%flavors)
    write(ostream,'(2x,a,30f10.4)') " full  orbital |",(op%meas_fullemptylines(1,iflavor),iflavor=1,op%flavors)
    write(ostream,'(2x,a,30f10.4)') " empty orbital |",(op%meas_fullemptylines(2,iflavor),iflavor=1,op%flavors)
  END IF
 !------------------------------------------------------------------------------------------
  !CALL FLUSH(op%ostream)
  IF ( ABS(((op%stats(4) *invSweeps*100.d0) / (op%stats(5) *invSweeps*100.d0) - 1.d0)) .GE. 0.02d0 &
   .OR. ABS(((op%stats(1)*invSweeps*100.d0) / (op%stats(2)*invSweeps*100.d0) - 1.d0)) .GE. 0.02d0 ) &
    THEN 
    CALL WARNALL("Ctqmcoffdiag_printQMC : bad statistic according to moves. Increase sweeps")
  END IF
  IF ( ABS(op%b_Noise+0.5)/0.5d0 .GE. 0.05d0 ) &
    CALL WARNALL("Ctqmcoffdiag_printQMC : bad statistic according to Noise. Increase sweeps")
!  IF ( ISNAN(op%a_Noise) .OR. ISNAN(op%a_Noise) ) &
!    CALL WARNALL("Ctqmcoffdiag_printQMC : NaN appeared. Increase sweeps    ")


END SUBROUTINE Ctqmcoffdiag_printQMC
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printGreen
!! NAME
!!  Ctqmcoffdiag_printGreen
!!
!! FUNCTION
!!  print green functions
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printGreen(op, oFileIn)

!Arguments ------------------------------------
  use m_io_tools, only : flush_unit
  TYPE(Ctqmcoffdiag)        , INTENT(IN)    :: op
  INTEGER  , OPTIONAL, INTENT(IN)    :: oFileIn
!Local variables ------------------------------
  INTEGER                            :: oFile
  INTEGER                            :: itime
  INTEGER                            :: sp1
  INTEGER                            :: iflavor,iflavorb
  INTEGER                            :: flavors, iflavor2 !,iflavor1,
  CHARACTER(LEN=4)                   :: cflavors
  CHARACTER(LEN=50)                  :: string
  DOUBLE PRECISION                   :: dt
  DOUBLE PRECISION                   :: sweeps

  !IF ( op%rank .NE. MOD(1,op%size)) RETURN
  IF ( op%rank .NE. MOD(op%size+1,op%size)) RETURN

  oFile = 40
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Gtau.dat")
  END IF
  OPEN(UNIT=43, FILE="Gtau_nd.dat")
  rewind(43)
  sp1     =  op%samples
  dt      =  op%beta / DBLE(sp1)
  sp1     =  sp1 + 1
  flavors =  op%flavors
  sweeps = DBLE(op%sweeps)*DBLE(op%size)

  IF ( op%opt_noise .EQ. 1) THEN
    WRITE(cflavors,'(I4)') (2*flavors+1)*2
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'ES22.14)'
    DO itime = 1, sp1
      WRITE(oFile,string) DBLE(itime-1)*dt, &
      (/ (op%Greens%oper(itime,iflavor,iflavor), iflavor=1, flavors) /), &
      (/ (op%abNoiseG(1,itime,iflavor)*(sweeps)**op%abNoiseG(2,itime,iflavor), iflavor=1, flavors) /)
    END DO
  ELSE
    WRITE(cflavors,'(I4)') (flavors+1)*2
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'ES22.14)'
    DO itime = 1, sp1
      WRITE(45,string) DBLE(itime-1)*dt, &
      (/ (op%Greens%oper(itime,iflavor,iflavor), iflavor=1, flavors) /)
      WRITE(oFile,string) DBLE(itime-1)*dt, &
      (/ (op%Greens%oper(itime,iflavor,iflavor), iflavor=1, flavors) /)
      WRITE(46,*) DBLE(itime-1)*dt, &
      & (/ ((op%Greens%oper(itime,iflavor,iflavorb), iflavor=1, flavors),iflavorb=1,flavors) /)
    END DO
    DO itime = 1, sp1
      WRITE(47,*) DBLE(itime-1)*dt, &
      & (/ ((op%Greens%oper(itime,iflavor,iflavorb), iflavor=1, flavors),iflavorb=1,flavors) /)
    END DO
!  --- Print full non diagonal Gtau in Gtau_nd.dat
    WRITE(cflavors,'(I4)') (flavors*flavors+1)
    write(47,*) "cflavors",cflavors
    string = '(1x,'//TRIM(ADJUSTL(cflavors))//'ES22.14)'
    write(47,*) string
    DO itime = 1, sp1
      WRITE(43,string) DBLE(itime-1)*dt, &
      & (/ ((op%Greens%oper(itime,iflavor,iflavorb), iflavorb=1, flavors),iflavor=1,flavors) /)
      WRITE(44,*) DBLE(itime-1)*dt, &
      & (op%Greens%oper(itime,iflavor,iflavor), iflavor=1, flavors)
      WRITE(44,string) DBLE(itime-1)*dt, &
      & (op%Greens%oper(itime,iflavor,iflavor), iflavor=1, flavors)
    END DO
      WRITE(43,*) 
  END IF
     DO iflavor = 1, flavors
       DO iflavor2 = 1, flavors
           write(4436,*) "#",iflavor,iflavor2
         do  itime=1,sp1
           write(4436,*) DBLE(itime-1)*dt,real(op%Greens%oper(itime,iflavor,iflavor2))
         enddo
           write(4436,*) 
       END DO
     END DO
     close(4436)

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)
  CLOSE(43)
  CLOSE(44)
  CLOSE(45)
  CLOSE(46)
  CLOSE(47)
  !call flush_unit(43)

END SUBROUTINE Ctqmcoffdiag_printGreen
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printD
!! NAME
!!  Ctqmcoffdiag_printD
!!
!! FUNCTION
!!  print individual double occupancy
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printD(op,oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)          , INTENT(IN)    :: op
  INTEGER  , OPTIONAL, INTENT(IN)    :: oFileIn
!Local variables ------------------------------
  INTEGER                            :: oFile
  INTEGER                            :: iflavor1
  INTEGER                            :: iflavor2

  !IF ( op%rank .NE. MOD(2,op%size)) RETURN
  IF ( op%rank .NE. MOD(op%size+2,op%size)) RETURN

  oFile = 41
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="D.dat")
  END IF

  DO iflavor1 = 1, op%flavors
    DO iflavor2 = iflavor1+1, op%flavors
      WRITE(oFile,'(1x,A8,I4,A1,I4,A3,ES21.14)') "Orbitals", iflavor1, "-", iflavor2, " : ", op%measDE(iflavor2,iflavor1)
    END DO
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmcoffdiag_printD
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printE
!! NAME
!!  Ctqmcoffdiag_printE
!!
!! FUNCTION
!!  print energy and noise 
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printE(op,oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)          , INTENT(IN)    :: op
  INTEGER  , OPTIONAL, INTENT(IN)    :: oFileIn
!Local variables ------------------------------
  INTEGER                            :: oFile
  DOUBLE PRECISION                   :: E
  DOUBLE PRECISION                   :: Noise

  !IF ( op%rank .NE. MOD(3,op%size)) RETURN
  IF ( op%rank .NE. MOD(op%size+3,op%size)) RETURN

  oFile = 42
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="BetaENoise.dat")
  END IF

  CALL Ctqmcoffdiag_getE(op,E,Noise)

  WRITE(oFile,'(1x,F3.2,A2,ES21.14,A2,ES21.14)') op%beta, "  ", E, "  ",  Noise

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmcoffdiag_printE
!!***

!#ifdef CTCtqmcoffdiag_ANALYSIS
!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printPerturbation
!! NAME
!!  Ctqmcoffdiag_printPerturbation
!!
!! FUNCTION
!!  print perturbation order
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printPerturbation(op, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)          , INTENT(IN)           :: op
  INTEGER  , OPTIONAL,  INTENT(IN)          :: oFileIn
!Local variables-------------------------------
  INTEGER                                   :: oFile
  INTEGER                                   :: iorder
  INTEGER                                   :: order
  INTEGER                                   :: iflavor
  INTEGER                                   :: flavors
  CHARACTER(LEN=2)                          :: a
  CHARACTER(LEN=50)                         :: string

  !IF ( op%rank .NE. MOD(4,op%size)) RETURN
  IF ( op%rank .NE. MOD(op%size+4,op%size)) RETURN
  IF ( op%opt_order .LE. 0 ) RETURN

  oFile = 43
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Perturbation.dat")
  END IF
    
  order        =  op%opt_order
  flavors      =  op%flavors

  WRITE(a,'(I2)') flavors
  string = '(I5,'//TRIM(ADJUSTL(a))//'F19.15)'
  DO iorder = 1, order
    WRITE(oFile,string) iorder-1, &
                (/ (op%measPerturbation(iorder, iflavor), iflavor=1, flavors) /)
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)
END SUBROUTINE Ctqmcoffdiag_printPerturbation
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printCorrelation
!! NAME
!!  Ctqmcoffdiag_printCorrelation
!!
!! FUNCTION
!!  print correlation fonctions
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printCorrelation(op, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)          , INTENT(IN)             :: op
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

  !IF ( op%rank .NE. MOD(5,op%size)) RETURN
  IF ( op%rank .NE. MOD(op%size+5,op%size)) RETURN
  IF ( op%opt_analysis .NE. 1 ) RETURN

  oFile = 44
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Correlation.dat")
  END IF

  sp1         =  op%samples
  dt          =  op%beta / sp1
  sp1         =  sp1 + 1
  flavors     =  op%flavors

  i = 3*flavors + 1
  WRITE(a,'(I2)') i
  WRITE(oFile,*) "# time  (/ (segement, antiseg, correl), i=1, flavor/)"
  string = '(1x,'//TRIM(ADJUSTL(a))//'F19.15)'
  DO itime = 1, sp1
    WRITE(oFile,string) DBLE(itime-1)*dt, &
                   (/ ( &
                   (/ ( op%measCorrelation(itime, i, iflavor), i=1,3) /) &
                   , iflavor=1, flavors) /)
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmcoffdiag_printCorrelation
!!***
!#endif

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_printSpectra
!! NAME
!!  Ctqmcoffdiag_printSpectra
!!
!! FUNCTION
!!  print fourier transform of time evolution of number of electrons
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_printSpectra(op, oFileIn)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag)          , INTENT(IN)             :: op
  INTEGER  , OPTIONAL, INTENT(IN)             :: oFileIn
!Local variables ------------------------------
  INTEGER                                     :: oFile
  INTEGER                                     :: flavors
  INTEGER                                     :: indDensity
  INTEGER                                     :: endDensity
  CHARACTER(LEN=4)                            :: a
  CHARACTER(LEN=16)                           :: formatSpectra

  !IF ( op%rank .NE. MOD(6,op%size)) RETURN
  IF ( op%opt_spectra .LT. 1 ) RETURN

  oFile = 45+op%rank
  a ="0000"
  WRITE(a,'(I4)') op%rank
  IF ( PRESENT(oFileIn) ) THEN
    oFile = oFileIn
  ELSE
    OPEN(UNIT=oFile, FILE="Markov_"//TRIM(ADJUSTL(a))//".dat")
  END IF

  flavors     =  op%flavors
  WRITE(a,'(I4)') flavors+1
  formatSpectra ='(1x,'//TRIM(ADJUSTL(a))//'ES22.14)'
  WRITE(oFile,*) "# freq[/hermalization] FFT"

  endDensity = SIZE(op%density,2)
  DO WHILE ( op%density(flavors+1,endDensity) .EQ. -1 )
    endDensity = endDensity -1
  END DO

  DO indDensity = 1, endDensity
    WRITE(oFile,formatSpectra) op%density(flavors+1,indDensity), op%density(1:flavors,indDensity)
  END DO

  IF ( .NOT. PRESENT(oFileIn) ) CLOSE(oFile)

END SUBROUTINE Ctqmcoffdiag_printSpectra
!!***

!!****f* ABINIT/m_Ctqmcoffdiag/Ctqmcoffdiag_destroy
!! NAME
!!  Ctqmcoffdiag_destroy
!!
!! FUNCTION
!!  destroy and deallocate all variables
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=ctqmc
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

SUBROUTINE Ctqmcoffdiag_destroy(op)

!Arguments ------------------------------------
  TYPE(Ctqmcoffdiag), INTENT(INOUT) :: op
!Local variables ------------------------------
  !INTEGER                  :: iflavor
  INTEGER                  :: flavors
  INTEGER                  :: i
  INTEGER                  :: j
  INTEGER                  :: k

  flavors = op%flavors

  CALL ImpurityOperator_destroy(op%Impurity)
  CALL BathOperatoroffdiag_destroy(op%Bath)
  CALL Vector_destroy(op%measNoise(1))
  CALL Vector_destroy(op%measNoise(2))

!sui!write(6,*) "before greenhyb_destroy in ctmqc_destroy"
  CALL GreenHyboffdiag_destroy(op%Greens)
!#ifdef CTCtqmcoffdiag_ANALYSIS
  FREEIF(op%measCorrelation)
  FREEIF(op%measPerturbation)
  FREEIF(op%meas_fullemptylines)
  FREEIF(op%measN)
  FREEIF(op%measDE)
  FREEIF(op%mu)
  FREEIF(op%hybri_limit)
  FREEIF(op%abNoiseG)
  IF ( ALLOCATED(op%measNoiseG) ) THEN
    DO i=1,2
      DO j = 1, op%flavors
        DO k= 1, op%samples+1
          CALL Vector_destroy(op%measNoiseG(k,j,i))
        END DO
      END DO
    END DO
    DT_FREE(op%measNoiseG)
  END IF
  FREEIF(op%density)
!#endif
  op%ostream        = 0
  op%istream        = 0
 
  op%sweeps         = 0
  op%thermalization = 0
  op%flavors        = 0
  op%samples        = 0
  op%beta           = 0.d0
!  op%seg_added      = 0.d0
!  op%anti_added     = 0.d0
!  op%seg_removed    = 0.d0
!  op%anti_removed   = 0.d0
!  op%seg_sign       = 0.d0
!  op%anti_sign      = 0.d0
  op%stats          = 0.d0
  op%swap           = 0.d0


  op%set  = .FALSE.
  op%done = .FALSE.
  op%init = .FALSE.
END SUBROUTINE Ctqmcoffdiag_destroy
!!***

END MODULE m_Ctqmcoffdiag
!!***

#undef CTQMC_SLICE1
#undef CTQMC_SLICE2
#undef CTQMC_SEGME 
#undef CTQMC_ANTIS 
#undef CTQMC_ADDED 
#undef CTQMC_REMOV 
#undef CTQMC_DETSI 

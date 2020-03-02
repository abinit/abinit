
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_GreenHyboffdiag
!! NAME
!!  m_GreenHyboffdiag
!! 
!! FUNCTION 
!!  Manage a green function for one orbital
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (B.Amadon, J. Denier and J. Bieder)
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
MODULE m_GreenHyboffdiag

 USE m_global
 USE m_MatrixHyb
 USE m_Vector
 USE m_VectorInt
 USE m_ListCdagC
 USE m_MapHyb
#ifdef HAVE_MPI2
 USE mpi
#endif
 
 IMPLICIT NONE


 public ::  GreenHyboffdiag_init
 public ::  GreenHyboffdiag_reset
 public ::  GreenHyboffdiag_clear
 public ::  GreenHyboffdiag_setOperW
 public ::  GreenHyboffdiag_measHybrid
 public ::  GreenHyboffdiag_getHybrid
 public ::  GreenHyboffdiag_setN
 public ::  GreenHyboffdiag_setMuD1
 public ::  GreenHyboffdiag_setMoments
 public ::  GreenHyboffdiag_backFourier
 public ::  GreenHyboffdiag_forFourier
 public ::  GreenHyboffdiag_print
 public ::  GreenHyboffdiag_destroy
 public ::  nfourier3

!!***

!!****t* m_GreenHyboffdiag/GreenHyboffdiag
!! NAME
!!  GreenHyboffdiag
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

 TYPE GreenHyboffdiag

  LOGICAL :: set = .FALSE.
   ! True if variable of type GreenHyboffdiag is initialized

  LOGICAL :: setT = .FALSE.
   ! True if variable oper contains data

  LOGICAL :: setW = .FALSE.
   ! True if variable oper_w contains data

  LOGICAL :: have_MPI = .FALSE.
   ! True if MPI is used.

  INTEGER :: setMk = 0
   ! setMk=0 is moments for Fourier transform are not computed

  INTEGER :: samples
   ! samples=imaginary time slices (dmftqmc_l+1)

  INTEGER :: measurements    
   ! number of measurements for the Green's function

  INTEGER :: factor
   ! if the move is not accepted, the statistic weight has to be
   ! increased for the current configuration.

  INTEGER :: MY_COMM
   ! MPI Communicator

  INTEGER :: size
   ! size=1 

  INTEGER :: rank
   ! rank=0 

  INTEGER :: Wmax
   ! samples-1 if frequency Green's function

  INTEGER :: iTech
   ! Precise if Frequency Green's function is computed or not

  INTEGER :: nflavors
   ! Number of flavors 

  DOUBLE PRECISION :: beta
   ! Inverse of temperature

  DOUBLE PRECISION :: inv_beta
   ! Temperature

  DOUBLE PRECISION :: delta_t
   ! 1/inv_dt

  DOUBLE PRECISION :: inv_dt
   ! (samples-1)/beta
  DOUBLE PRECISION :: signvaluemeas

  DOUBLE PRECISION :: signvalueold

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: oper
   ! oper(samples)

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: omega
   ! omega(Wmax) 

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Mk
   ! Moments for FT

  COMPLEX(KIND=8)  , ALLOCATABLE, DIMENSION(:,:,:) :: oper_w 
   ! Frequency Green's function

  COMPLEX(KIND=8)  , ALLOCATABLE, DIMENSION(:) :: oper_w_old
   ! Old frequency Green's function (not used)

  TYPE(Vector)                            :: oper_old          
   ! useless data

  TYPE(VectorInt)                         :: index_old          
   ! useless data

  TYPE(MapHyb), ALLOCATABLE, DIMENSION(:,:)  :: map
   ! value of time and Green's functions computed in GreenHyboffdiag_measHybrid
   ! These values are used to fill op%oper in the same routine.

 END TYPE GreenHyboffdiag
!!***

CONTAINS
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_init
!! NAME
!!  GreenHyboffdiag_init
!!
!! FUNCTION
!!  Initialize and allocate
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  samples=imaginary time slices
!!  beta=inverse temperature
!!  iTech=SHOULD NOT BE USED => BUGGY
!!  MY_COMM=mpi_communicator
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

SUBROUTINE GreenHyboffdiag_init(op, samples, beta,nflavors,iTech,MY_COMM)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)     , INTENT(INOUT) :: op
  INTEGER         , INTENT(IN   ) :: samples
  DOUBLE PRECISION, INTENT(IN   ) :: beta
  INTEGER , INTENT(IN   ) :: nflavors
  !INTEGER          , INTENT(IN   ) :: Wmax
  INTEGER, OPTIONAL, INTENT(IN   ) :: iTech
  INTEGER, OPTIONAL, INTENT(IN   ) :: MY_COMM
!Local variables ------------------------------
  INTEGER                         :: iflavor,iflavorbis,sp1
  DOUBLE PRECISION                :: dt
#ifdef HAVE_MPI
  INTEGER          :: ierr
#endif

  IF ( PRESENT(MY_COMM)) THEN
#ifdef HAVE_MPI
    op%have_MPI = .TRUE.
    op%MY_COMM = MY_COMM
    CALL MPI_Comm_rank(op%MY_COMM, op%rank, ierr)
    CALL MPI_Comm_size(op%MY_COMM, op%size, ierr)
#else
    CALL WARN("GreenHyboffdiag_init : MPI is not used                                    ")
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

  sp1             = samples + 1
  op%samples      = sp1
  op%measurements = 0
  op%nflavors     = nflavors
  op%beta         = beta
  op%inv_beta     = 1.d0 / beta
  op%inv_dt       = DBLE(samples) * op%inv_beta
  dt              = 1.d0 / op%inv_dt
  op%delta_t      = dt
  !op%Wmax         = Wmax
  op%Wmax         = -1
  FREEIF(op%Mk)
  MALLOC(op%Mk,(nflavors,nflavors,3))
  FREEIF(op%oper)
  MALLOC(op%oper,(sp1,nflavors,nflavors))
  ! If we want to measure in frequences
  ! let assume we first have "samples" frequences
  IF ( PRESENT(iTech) ) THEN
    op%iTech = iTech
    SELECT CASE (op%iTech)
    CASE (GREENHYB_TAU)  ! omega
      op%iTech = GREENHYB_TAU
    CASE (GREENHYB_OMEGA)  ! omega
      op%Wmax = samples
      FREEIF(op%oper_w)
      MALLOC(op%oper_w,(1:op%Wmax,nflavors,nflavors))
      FREEIF(op%oper_w_old)
      MALLOC(op%oper_w_old,(1:op%Wmax))
      op%oper_w     = CMPLX(0.d0,0.d0,8)
      op%oper_w_old = CMPLX(0.d0,0.d0,8)
      FREEIF(op%omega)
      MALLOC(op%omega,(1:op%Wmax))
      op%omega = (/ ((2.d0 * DBLE(sp1) - 1.d0)*ACOS(-1.d0)*op%inv_beta, sp1=1, op%Wmax) /)
    END SELECT
  ELSE
    op%iTech = GREENHYB_TAU
  END IF
  ! end if
  CALL Vector_init(op%oper_old,10000)
  CALL VectorInt_init(op%index_old,10000)
  DT_FREEIF(op%map)
  MALLOC(op%map,(nflavors,nflavors))
  do iflavor=1,nflavors
    do iflavorbis=1,nflavors
      CALL MapHyb_init(op%map(iflavor,iflavorbis),10000)
    enddo
  enddo

  op%oper       = 0.d0
  op%signvaluemeas = 0.d0
  op%signvalueold = 0.d0
  op%set        = .TRUE.
  op%factor     = 1
  op%setMk      = 0
  op%Mk         = 0.d0
END SUBROUTINE GreenHyboffdiag_init
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_reset
!! NAME
!!  GreenHyboffdiag_reset
!!
!! FUNCTION
!!  reset green function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
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

SUBROUTINE GreenHyboffdiag_reset(op)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)     , INTENT(INOUT) :: op

  CALL GreenHyboffdiag_clear(op)
  op%setMk        = 0
  op%Mk           = 0.d0
  op%setT         = .FALSE.
  op%setW         = .FALSE.
END SUBROUTINE GreenHyboffdiag_reset
!!***


!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_clear
!! NAME
!!  GreenHyboffdiag_clear
!!
!! FUNCTION
!!  clear green function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
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

SUBROUTINE GreenHyboffdiag_clear(op)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)     , INTENT(INOUT) :: op
  INTEGER :: iflavor,iflavorbis

  !CALL Vector_clear(op%oper_old)
  !CALL VectorInt_clear(op%index_old)
  do iflavor=1,op%nflavors
    do iflavorbis=1,op%nflavors
      CALL MapHyb_clear(op%map(iflavor,iflavorbis))
    enddo
  enddo
  op%measurements = 0
  IF ( ALLOCATED(op%oper) ) &
  op%oper         = 0.d0
  op%signvaluemeas = 0.d0
  op%signvalueold = 1.d0
  IF ( op%iTech .EQ. GREENHYB_OMEGA ) THEN
    IF ( ALLOCATED(op%oper_w) ) &
    op%oper_w       = CMPLX(0.d0,0.d0,8)
    IF ( ALLOCATED(op%oper_w_old) ) &
    op%oper_w_old   = CMPLX(0.d0,0.d0,8)
  END IF
  op%factor       = 0
END SUBROUTINE GreenHyboffdiag_clear
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_setOperW
!! NAME
!!  GreenHyboffdiag_setOperW
!!
!! FUNCTION
!!  set Green function in frequencies
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  Gomega=Input values
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

SUBROUTINE GreenHyboffdiag_setOperW(op, Gomega)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)          , INTENT(INOUT) :: op
  COMPLEX(KIND=8), DIMENSION(:,:,:), INTENT(IN   ) :: Gomega
!Loval variables ------------------------------
  INTEGER :: tail

  tail = SIZE(Gomega,1)
  IF ( .NOT. op%set ) &
    CALL ERROR("GreenHyboffdiag_setOperW : Uninitialized GreenHyboffdiag structure")
  IF ( ALLOCATED(op%oper_w) ) THEN
    IF ( SIZE(op%oper_w) .NE. tail ) THEN
      FREE(op%oper_w)
      MALLOC(op%oper_w,(1:tail,op%nflavors,op%nflavors))
    END IF
  ELSE
    MALLOC(op%oper_w,(1:tail,op%nflavors,op%nflavors))
  END IF
  op%oper_w(:,:,:) = Gomega(:,:,:)
  op%Wmax = tail
  op%setW = .TRUE.
END SUBROUTINE GreenHyboffdiag_setOperW
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_measHybrid
!! NAME
!!  GreenHyboffdiag_measHybrid
!!
!! FUNCTION
!!  Measure Green's function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  Mmatrix=M matrix for the current flavor
!!  ListCdagC_1=list of all creator and annhilator operators
!!  updated=should we accumulate or not
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

SUBROUTINE GreenHyboffdiag_measHybrid(op, Mmatrix, ListCdagC_1, updated,signvalue,activeflavor)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)    , INTENT(INOUT) :: op
  TYPE(MatrixHyb)   , INTENT(IN   ) :: Mmatrix
  TYPE(ListCdagC)   , INTENT(IN   ) :: ListCdagC_1(op%nflavors)
  DOUBLE PRECISION  , INTENT(IN   ) :: signvalue
  LOGICAL        , INTENT(IN   ) :: updated
  INTEGER, OPTIONAL  , INTENT(IN  ) :: activeflavor
!Local variables ------------------------------
  INTEGER                        :: iC
  INTEGER                        :: iCdag
  INTEGER                        :: tail
  INTEGER                        :: tailbis
  !INTEGER                        :: index
  INTEGER                        :: idx_old
  INTEGER                        :: old_size
!  INTEGER                        :: omegaSamples
!  INTEGER                        :: iomega
  INTEGER                        :: iflavor
  INTEGER                        :: iflavorbis
  INTEGER                        :: iC_m,iC_m_add
  INTEGER                        :: iCdag_m,iCdag_m_add
  INTEGER                        :: stail !,ii
!  DOUBLE PRECISION               :: pi_invBeta
  DOUBLE PRECISION               :: mbeta_two
  DOUBLE PRECISION               :: beta 
  DOUBLE PRECISION               :: beta_tc
  DOUBLE PRECISION               :: tcbeta_tc
  DOUBLE PRECISION               :: inv_dt 
  DOUBLE PRECISION               :: tC,tc_phys
  DOUBLE PRECISION               :: tCdag
  DOUBLE PRECISION               :: time
  DOUBLE PRECISION               :: signe,signe2
  DOUBLE PRECISION               :: argument
  INTEGER                        :: iflavorbegin,iflavorend,prtopt
  !DOUBLE PRECISION               :: taupi_invbeta
  !COMPLEX(KIND=8)                   :: cargument
  !COMPLEX(2*8)                   :: base_exp
  !COMPLEX(2*8)                   :: increm_exp
  !write(6,*) "measHybrid"
  prtopt=0
  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_measHybrid : green operator not set         ")
  stail=0
  do iflavor=1,op%nflavors
    stail=stail + ListCdagC_1(iflavor)%tail
  enddo
  iflavorbegin = 1
  iflavorend   = op%nflavors

  if(present(activeflavor)) then
    if(activeflavor.ne.0) then
!sui!write(6,*) "measHybrid activeflavor",activeflavor
      iflavorbegin = activeflavor
      iflavorend   = activeflavor
    endif
  endif

  IF ( stail .NE. Mmatrix%tail ) &
      CALL ERROR("GreenHyboffdiag_measHybrid : ListCdagC & M unconsistent     ")

  IF ( updated .EQV. .TRUE. ) THEN ! NEW change in the configuration
    ! FIXME SHOULD be much more faster
    

    SELECT CASE(op%iTech)
    CASE (GREENHYB_TAU)
      argument = DBLE(op%factor)
!     At the beginning old_size=0, then it increases
!     until 
!     for all values of iC, increment green%oper with the value of the
!     Green's function in listDBLE(iC) obtained from previous iteration
!     (below)
      ! ===============================================================
      ! An update has been done. So the Green's function will change
      ! It is thus the good moment to store the previous Green's
      ! function with argument, the number of times this Green's
      ! function has been constant
      ! ===============================================================
      DO iflavor=1, op%nflavors
        DO iflavorbis=1, op%nflavors
          old_size = op%map(iflavor,iflavorbis)%tail
          !write(6,*) "size listDBLE",size(op%map(iflavor,iflavorbis)%listDBLE)
        !sui!write(6,*) " measHybrid",old_size,iflavor,iflavorbis
          DO iC = 1, old_size
            if(op%map(iflavor,iflavorbis)%listINT(iC)==0) then
              !write(6,*) "listINT(iC)=",iC,op%map(iflavor,iflavorbis)%listINT(iC)
            endif
            !write(6,*) " measHybrid  iflavor,iflavorbis,iC listINT ",iflavor,iflavorbis,iC,op%map(iflavor,iflavorbis)%listINT(iC)
               !write(6,*) " measHybrid  listDBLE ",iflavor,iflavorbis,iC,op%map(iflavor,iflavorbis)%listDBLE(iC),argument
            op%oper(op%map(iflavor,iflavorbis)%listINT(iC),iflavor,iflavorbis) =                &
                           op%oper(op%map(iflavor,iflavorbis)%listINT(iC),iflavor,iflavorbis) &
                         + op%map(iflavor,iflavorbis)%listDBLE(iC) * op%signvalueold *  argument
           !if(op%map(iflavor,iflavorbis)%listINT(iC)==1.and.iflavor==iflavorbis) then
          !  if(iflavor==iflavorbis) then
          !   !sui!write(6,*) "G(0)", op%map(iflavor,iflavorbis)%listINT(iC),op%map(iflavor,iflavorbis)%listDBLE(iC) * op%signvalueold,op%oper(op%map(iflavor,iflavorbis)%listINT(iC),iflavor,iflavorbis),iflavor
          !  endif
          !  if(iflavor==1.and.iflavorbis==6.and.op%map(iflavor,iflavorbis)%listINT(iC)==1) then
          !          !prt!if(prtopt==1) write(6,*) "G16(0)", op%map(iflavor,iflavorbis)%listDBLE(iC) * op%signvalueold*argument,op%oper(op%map(iflavor,iflavorbis)%listINT(iC),iflavor,iflavorbis),op%signvalueold
          !  endif
          !  if(iflavor==6.and.iflavorbis==1.and.op%map(iflavor,iflavorbis)%listINT(iC)==1) then
          !          !prt!if(prtopt==1) write(6,*) "G61(0)", op%map(iflavor,iflavorbis)%listDBLE(iC) * op%signvalueold*argument,op%oper(op%map(iflavor,iflavorbis)%listINT(iC),iflavor,iflavorbis),op%signvalueold
          !  endif
          END DO
      ! tail**2 is the number of possible t-t'
      ! MapHyb_setSize with resize map tail*tail will thus be the new
      ! op%map%tail
      ! update size of map and map%tail
          CALL MapHyb_setSize(op%map(iflavor,iflavorbis),&
&          ListCdagC_1(iflavor)%tail*ListCdagC_1(iflavorbis)%tail)
        END DO
      END DO
      op%signvaluemeas = op%signvaluemeas + op%signvalueold * argument
      op%measurements = op%measurements + op%factor
    !sui!write(6,*) "   measurements", op%measurements
         !sui! write(6,*) "                  signvaluemeas",op%signvaluemeas,op%signvalueold*argument
         !sui! write(6,*) "                  signvaluemeas/measurements",op%signvaluemeas/op%measurements
  
      ! This is new measurement, thus op%factor should be put to one
      op%factor = 1


      ! initialized index idx_old for the doubles loops over flavors and segments.

      ! setup usefull quantities
      beta   =  op%beta
      mbeta_two = -(beta*0.5d0)
      inv_dt =  op%inv_dt

      ! WARNING time is not the time but just a temporary variable.
      ! Index Time has been calculated previously and is in mat_tau

      ! initialized index for each annihilation time of a segment for a given flavor

      ! initialized index for each creation time of a segment for another  flavor

      iC_m=0
      iC_m_add=0
      DO iflavor=1,op%nflavors
        tail=ListCdagC_1(iflavor)%tail
        !write(6,*) " measHybrid  iflavor",iflavor,tail

        iCdag_m=0
        iCdag_m_add=0
        DO iflavorbis=1,op%nflavors
          tailbis=ListCdagC_1(iflavorbis)%tail
          !write(6,*) " measHybrid  iflavorbis",iflavorbis,tailbis
          idx_old = 0

          DO iC  = 1, tail
            ! tC is the annihilation (C_) time for segment iC and flavor iflavor
            !-------------------------------------------------------------------
            tC   = ListCdagC_1(iflavor)%list(iC,C_)

            !iC_m=iC_m+1  ! for Mmatrix%mat
            ! For each flavor iflavor, iC start at \sum_{iflavor1<iflavor} tail(iflavor1)
            ! It thus explains the presence of iC_m_add (same below for iCdag_m_add)
            ! ---------------------------------------------------------------------------------
            iC_m=iC_m_add+iC
            beta_tc = beta - tC
            tcbeta_tc = tC * beta_tc

             !write(6,*) " measHybrid  iC_m",iC_m
             !write(6,*) " measHybrid  tailbis",tailbis
            DO iCdag = 1, tailbis
              !iCdag_m=iCdag_m+1
              iCdag_m=iCdag_m_add+iCdag
             !write(6,*) " measHybrid  iCdag_m",iCdag_m

              ! tCdag is the creation time for segment iCdag and flavor iflavorbis
              tCdag  = ListCdagC_1(iflavorbis)%list(iCdag,Cdag_)

!  ---        time is equivalent to  time=(tc-tcdag)*(beta-tc) and is only
!  ---        useful for signe
              time = tcbeta_tc - tCdag*beta_tc
  
              !signe = SIGN(1.d0,time)
              !time = time + (signe-1.d0)*mbeta_two
              !signe = signe * SIGN(1.d0,beta-tC)
              !signe = SIGN(1.d0,time) * SIGN(1.d0,beta-tC)
              tc_phys=tc
              if(tc>beta) tc_phys=tc-beta
              signe2=SIGN(1.d0,tc_phys-tcdag)

              if(iflavor==iflavorbis) signe = SIGN(1.d0,time)
              if(iflavor/=iflavorbis) signe = signe2
             ! signe = SIGN(1.d0,tc-tcdag)
              ! --- tc>tcdag and beta>tc signe=1  ! segment in the middle or antisegment at the edge
              !                                   ! tc-tcdag  > 0
              ! --- tc<tcdag and beta<tc signe=1  ! never
              ! --- tc>tcdag and beta<tc signe=-1 ! segment  at the edges
              !                                   ! tc'-tcdag < 0 (with tc'=tc-beta)  -> signe < 0
              ! --- tc<tcdag and beta>tc signe=-1 ! antisegment in the middle 
              !                                   ! tc-tcdag  < 0 (with tc'=tc-beta)  -> signe < 0
              ! 22/09/14:
              ! ListCdagC_1 is the list of segment, so we are dealing
              ! only with segment here. However all combination of Cdag
              ! and C are taken, this it is possible that tc<tcdag

              ! 21/10/14: Wagt is important are the true times (between
              ! 0 and beta). If tauC>tauCdag signe=+1
              !              If tauC<tauCdag signe=-1
              ! if(tc<tcdag.and.(iflavor==iflavorbis)) then
              !   write(6,*)  ListCdagC_1(iflavorbis)%tail
              !   do ii=1, ListCdagC_1(iflavorbis)%tail
              !     write(6,*)  ListCdagC_1(iflavorbis)%list(ii,1), ListCdagC_1(iflavorbis)%list(ii,2)
              !   enddo
              !   write(6,*) "tc<tcdag", tc,tcdag,beta,iflavor,iflavorbis
              !   stop
              ! endif

              if(tc-tcdag>beta) then
             !   write(6,*) " tc-tcdag > beta ", tcdag-tc,beta
              endif
              !if(tc>beta) then
              !  write(6,*) " TC>BETA"
              !  write(6,*) " iflavor,iflavorbis",iflavor,iflavorbis
              !  write(6,*) " ic,icdag          ",ic,icdag
              !  write(6,*) " signe             ",signe
              !  write(6,*) " tc,tcdag          ",tc,tcdag
              !  write(6,*) " Mmatrix%mat       ",Mmatrix%mat(iCdag_m,iC_m)
              !  write(6,*) " Mmatrix%mat_tau   ",Mmatrix%mat_tau(iCdag_m,iC_m)
              !endif
              ! Si iflavor/=iflavorbis, tc-tcdag can be negative..so in
              ! this case, on should add beta to tc-tcdag with the minus
              ! sign. NOT DONE HERE??
  
!             ----- Compute the Green's function as the value of the matrix M for times iCdag and iC.
              argument = signe*Mmatrix%mat(iCdag_m,iC_m)
  
              !index = INT( ( time * inv_dt ) + 1.5d0 )
              !IF (index .NE. Mmatrix%mat_tau(iCdag,iC)) THEN
              !  WRITE(*,*) index, Mmatrix%mat_tau(iCdag,iC)
              !!  CALL ERROR("Plantage")
              !END IF
  
              idx_old = idx_old + 1

              ! --- define the  value of listDBLE as a function of idx_old
              op%map(iflavor,iflavorbis)%listDBLE(idx_old) = argument
              !write(6,*) " measHybrid  listDBLE2 ",iflavor,iflavorbis,idx_old,argument
              !op%map%listINT(idx_old)  = index

              ! --- define the new corresponding value of listINT(idx_old) from mat_tau (integers)
              ! --- idx_old has no meaning but listINT(idx_old) has.
              op%map(iflavor,iflavorbis)%listINT(idx_old)  = Mmatrix%mat_tau(iCdag_m,iC_m)
             !write(6,*) " measHybrid  idx_old listINT ",idx_old,op%map(iflavor,iflavorbis)%listINT(idx_old)
             !write(6,*) " measHybrid  iCdag_m, iC_m, mat_tau",iCdag_m,iC_m,Mmatrix%mat_tau(iCdag_m,iC_m)
            !  if(iflavor==1.and.iflavorbis==2.and.op%map(iflavor,iflavorbis)%listINT(idx_old)==1) then
            !    !prt!if(prtopt==1) write(6,*) "---------------------------"
            !    !prt!if(prtopt==1) write(6,*) "GG12(0)", op%map(iflavor,iflavorbis)%listINT(idx_old),op%map(iflavor,iflavorbis)%listDBLE(idx_old),tcdag,tc,signe,signe2
            !    !prt!if(prtopt==1) write(6,*) "       ", tc-tcdag,tc_phys-tcdag
            !    do ii=1, tail
            !      !prt!if(prtopt==1)  write(6,*) ii, ListCdagC_1(iflavor)%list(ii,1), ListCdagC_1(iflavor)%list(ii,2)
            !    enddo
            !    do ii=1, tailbis
            !      !prt!if(prtopt==1)  write(6,*) ii, ListCdagC_1(iflavorbis)%list(ii,1), ListCdagC_1(iflavorbis)%list(ii,2)
            !    enddo
            !    !prt!if(prtopt==1) write(6,*) "---------------------------"
            !  endif
            !  if(iflavor==1.and.iflavorbis==2.and.op%map(iflavor,iflavorbis)%listINT(idx_old)==999) then
            !    !prt!if(prtopt==1) write(66,*) "---------------------------"
            !    !prt!if(prtopt==1) write(66,*) "GG12(0)", op%map(iflavor,iflavorbis)%listINT(idx_old),op%map(iflavor,iflavorbis)%listDBLE(idx_old),tcdag,tc,signe,signe2
            !    !prt!if(prtopt==1) write(66,*) "       ", tc-tcdag,tc_phys-tcdag
            !    do ii=1, tail
            !      !prt!if(prtopt==1)  write(66,*) ii, ListCdagC_1(iflavor)%list(ii,1), ListCdagC_1(iflavor)%list(ii,2)
            !    enddo
            !    do ii=1, tailbis
            !      !prt!if(prtopt==1)  write(66,*) ii, ListCdagC_1(iflavorbis)%list(ii,1), ListCdagC_1(iflavorbis)%list(ii,2)
            !    enddo
            !    !prt!if(prtopt==1) write(66,*) "---------------------------"
            !  endif
            !  !if(iflavor==2.and.iflavorbis==1.and.op%map(iflavor,iflavorbis)%listINT(idx_old)==1) then
              !   !prt!if(prtopt==1) write(6,*) "GG21(0)", op%map(iflavor,iflavorbis)%listINT(idx_old),op%map(iflavor,iflavorbis)%listDBLE(idx_old),tcdag,tc,signe
              !endif


            END DO
          END DO
        !  do ii=1,tail*tailbis
        !   !write(6,*) " measHybrid  ii,op%map(iflavor,iflavorbis)%listINT(ii)", ii,op%map(iflavor,iflavorbis)%listINT(ii)
        !  enddo
          iCdag_m_add=iCdag_m_add+tailbis
        END DO ! iflavorbis
       iC_m_add=iC_m_add+tail
      END DO ! iflavor
      op%signvalueold = signvalue
    CASE (GREENHYB_OMEGA)
    !  argument = DBLE(op%factor)
    !  DO iomega = 1, omegaSamples
    !    op%oper_w(iomega) = op%oper_w(iomega) + op%oper_w_old(iomega) * argument
    !  END DO
    !  op%measurements = op%measurements + op%factor

    !  op%factor = 1
    !  beta   =  op%beta
    !  mbeta_two = -(beta*0.5d0)
    !  pi_invBeta = ACOS(-1.d0)/beta
    !  omegaSamples = op%samples-1
    !  DO iC  = 1, tail
    !    tC   = ListCdagC_1%list(iC,C_)
    !    DO iCdag = 1, tail
    !      tCdag  = ListCdagC_1%list(iCdag,Cdag_)
    !      time = tC - tCdag

    !      signe = SIGN(1.d0,time)
    !      time = time + (signe-1.d0)*mbeta_two
    !      signe = signe * SIGN(1.d0,beta-tC)
    !      argument = signe*Mmatrix%mat(iCdag,iC)

    !      DO iomega = 1, omegaSamples
    !        !op%oper_w_old(iomega) = Mmatrix%mat_tau(iCdag,iC)*CMPLX(0.d0,argument)
    !        op%oper_w_old(iomega) = EXP(CMPLX(0.d0,op%omega(iomega)*time))*CMPLX(0.d0,argument)
    !      END DO
    !    END DO
    !  END DO
    END SELECT
  ELSE
    op%factor = op%factor + 1
  END IF
END SUBROUTINE GreenHyboffdiag_measHybrid
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_getHybrid
!! NAME
!!  GreenHyboffdiag_getHybrid
!!
!! FUNCTION
!!  reduce green function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
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

SUBROUTINE GreenHyboffdiag_getHybrid(op)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag), INTENT(INOUT) :: op

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_getHybrid : green operator not set          ")

  SELECT CASE(op%iTech)
  CASE (GREENHYB_TAU)
    op%oper = -(op%oper * op%inv_beta) / (DBLE(op%measurements) * op%delta_t)
  !sui!write(6,*) "measurements",op%measurements,op%delta_t,op%inv_beta
  !sui!write(6,*) "signevaluemeas meas",op%signvaluemeas,op%measurements
    op%signvaluemeas = op%signvaluemeas / DBLE(op%measurements)
   ! print*, "op%oper",op%oper(1,1,1)
  !sui!write(6,*) "signevaluemeas/meas",op%signvaluemeas
   ! print*, "signevaluemeas/meas",op%signvaluemeas
    op%setT = .TRUE.
  CASE (GREENHYB_OMEGA)
    op%oper_w = -(op%oper_w * op%inv_beta) / (DBLE(op%measurements) * op%delta_t)
    op%setW = .TRUE.
    CALL GreenHyboffdiag_backFourier(op)
  END SELECT

END SUBROUTINE GreenHyboffdiag_getHybrid
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_setN
!! NAME
!!  GreenHyboffdiag_setN
!!
!! FUNCTION
!!  impose number of electrons for this flavor
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  N=number of electrons
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

SUBROUTINE GreenHyboffdiag_setN(op,N)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)    , INTENT(INOUT)    :: op
  DOUBLE PRECISION  , INTENT(IN   )    :: N(op%nflavors)
  INTEGER :: iflavor,iflavor2
  !COMPLEX(KIND=8) :: tmpoper

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_setN: green op%operator not set                ")
  DO iflavor=1, op%nflavors
   ! write(6,*) "iflavor",-N(iflavor)*op%signvaluemeas ,2*op%oper(op%samples,iflavor,iflavor),(N(iflavor)-1.d0)*op%signvaluemeas
    ! the mulplication by signvaluemeas is necessary because N is
    ! exactly the number of electrons in the flavor iflavor whereas
    ! op%oper is not exact, because it still has to be divided by
    ! signvaluemeas after the MPIREDUCE
    op%oper(1,iflavor,iflavor) = (N(iflavor) - 1.d0)*op%signvaluemeas
    op%oper(op%samples,iflavor,iflavor) = - N(iflavor)*op%signvaluemeas
    !op%oper(op%samples,iflavor,iflavor) = 2*op%oper(op%samples,iflavor,iflavor)
    !op%oper(1,iflavor,iflavor) = 2*op%oper(1,iflavor,iflavor)
    DO iflavor2=1, op%nflavors
      if(iflavor/=iflavor2) then
              ! UNEXPLAINED but MANDATORY to have exact results for U=0 nspinor=4 with pawspnorb=0
              ! Correction: The fact 2 is necessary for edge points because the points are at the
              ! edges. 
              ! It is of course necessary to fulfill exact results (U=0).
    !tmpoper=(op%oper(op%samples,iflavor,iflavor2)-op%oper(1,iflavor,iflavor2))
    !op%oper(op%samples,iflavor,iflavor2) = tmpoper
    !op%oper(1,iflavor,iflavor2) = -tmpoper
    !op%oper(op%samples,iflavor,iflavor2) = 2*op%oper(op%samples,iflavor,iflavor2)
    !op%oper(1,iflavor,iflavor2) = 2*op%oper(1,iflavor,iflavor2)
        op%oper(op%samples,iflavor,iflavor2) = 2*op%oper(op%samples,iflavor,iflavor2)
        op%oper(1,iflavor,iflavor2) = 2*op%oper(1,iflavor,iflavor2)
      endif
    ENDDO
  ENDDO
END SUBROUTINE GreenHyboffdiag_setN
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_setMuD1
!! NAME
!!  GreenHyboffdiag_setMuD1
!!
!! FUNCTION
!!  Set first moments for G
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  mu=energy level (irrespectige with fermi level)
!!  d1=first moment of hybridization function ("K")
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

SUBROUTINE GreenHyboffdiag_setMuD1(op,iflavor,iflavor2,mu,d1)

!Arguments ------------------------------------
!Arguments ------------------------------------
!scalars
  DOUBLE PRECISION, INTENT(IN   ) :: mu
  DOUBLE PRECISION, INTENT(IN   ) :: d1
  INTEGER         , INTENT(IN   ) :: iflavor
  INTEGER         , INTENT(IN   ) :: iflavor2
!type
  TYPE(GreenHyboffdiag)  , INTENT(INOUT) :: op
!Local variables ------------------------------------
  DOUBLE PRECISION                :: mu2
!*********************************************************************

  ABI_UNUSED((/d1/))
  
  mu2=0
  if(iflavor==iflavor2) mu2=mu

  if(iflavor==iflavor2) then
    op%Mk(iflavor,iflavor2,3) = -d1-(mu*mu)
    !op%Mk(iflavor,iflavor2,3) = -(mu*mu)
    op%Mk(iflavor,iflavor2,2) = -mu
 !sui!write(6,*) "setmud1",iflavor,iflavor2, op%Mk(iflavor,iflavor2,2), op%Mk(iflavor,iflavor2,3)
  else
    op%Mk(iflavor,iflavor2,3) = 0.d0
    op%Mk(iflavor,iflavor2,2) = 0.d0
  endif
  op%setMk = op%setMk + 1
!write(6,*) "mom1",op%Mk(iflavor,iflavor2,:)
END SUBROUTINE GreenHyboffdiag_setMuD1
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_setMoments
!! NAME
!!  GreenHyboffdiag_setMoments
!!
!! FUNCTION
!!  Compute full moments
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Greenb
!!  u1_iflavor1=\sum_{iflavor2} U_{iflavor2,iflavor1} N_iflavor2
!!    (useful for first moment)  
!!  u2=\sum_{iflavor1,iflavor2,iflavor3} U_{iflavor1,iflavor2} N_iflavor2
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
!! CHI
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE GreenHyboffdiag_setMoments(op,iflavor1,iflavor1b,u1,u2,u3)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)  , INTENT(INOUT) :: op
  DOUBLE PRECISION, INTENT(IN   ) :: u1
  DOUBLE PRECISION, INTENT(IN   ) :: u2
  DOUBLE PRECISION, INTENT(IN   ) :: u3
  INTEGER         , INTENT(IN   ) :: iflavor1
  INTEGER         , INTENT(IN   ) :: iflavor1b
  
  if(iflavor1==iflavor1b) then
    op%Mk(iflavor1,iflavor1b,1) = -1.d0
!   c_a(3)=-d1-mu*mu-2(-mu)(\sum_{b.ne.a} Uab nb) 
    op%Mk(iflavor1,iflavor1b,3) = op%Mk(iflavor1,iflavor1b,3) - 2.d0*(op%Mk(iflavor1,iflavor1b,2)*u1)

!   c_a(2)=-mu+\sum_{b.ne.a} Uab n_b
    op%Mk(iflavor1,iflavor1b,2) = op%Mk(iflavor1,iflavor1b,2) + u1
 !sui!write(6,*) "setmiments",iflavor1,iflavor1b,u1

!   c_a(3)=c_a(3) + \sum Uab^2 nb + \sum Uba Uca <nbnc> 
!   ie c_a(3)=-d1+mu*mu-2mu*\sumb Uab nb + \sum Uab^2 nb + \sum Uba Uca <nbnc> 
    op%Mk(iflavor1,iflavor1b,3) = op%Mk(iflavor1,iflavor1b,3) - u2
  else
    op%Mk(iflavor1,iflavor1b,1) = 0.d0
    op%Mk(iflavor1,iflavor1b,2) = u3
    op%Mk(iflavor1,iflavor1b,3) = 0.d0
  endif
!write(6,*) "mom",iflavor1,iflavor1b, op%Mk(iflavor1,iflavor1b,:)

  op%setMk = op%setMk + 1

END SUBROUTINE GreenHyboffdiag_setMoments
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_backFourier
!! NAME
!!  GreenHyboffdiag_backFourier
!!
!! FUNCTION
!!  perform back fourier transform
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  dvgc=divergence parameter
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

SUBROUTINE GreenHyboffdiag_backFourier(op,dvgc,func,hybri_limit,opt_hybri_limit)

 use m_fstrings,       only : int2char4

#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(GreenHyboffdiag)            , INTENT(INOUT) :: op
  DOUBLE PRECISION, OPTIONAL, INTENT(IN   ) :: dvgc
  CHARACTER(len=5)  ,OPTIONAL, INTENT(IN) :: func
  COMPLEX(KIND=8), DIMENSION(op%nflavors,op%nflavors), OPTIONAL, INTENT(IN) :: hybri_limit
  INTEGER, OPTIONAL, INTENT(IN) :: opt_hybri_limit
!Local variables ------------------------------
  INTEGER :: itau
  INTEGER :: iomega
  INTEGER :: omegaSamples
  INTEGER :: tauSamples
  INTEGER :: tauBegin
  INTEGER :: tauEnd
  INTEGER :: delta
  INTEGER :: residu
  INTEGER :: iflavor1
  INTEGER :: iflavor2,unitnb !,unitnb1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: counts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: displs
  DOUBLE PRECISION :: A,AA ! Correction factor
  COMPLEX(KIND=8) :: B,BB ! Correction factor
  COMPLEX(KIND=8) :: C !,CC ! Correction factor
  DOUBLE PRECISION :: inv_beta
  DOUBLE PRECISION :: pi_invBeta
  DOUBLE PRECISION :: two_invBeta
  DOUBLE PRECISION :: minusDt
  DOUBLE PRECISION :: minusOmegaTau
  DOUBLE PRECISION :: omegaa
  DOUBLE PRECISION :: minusTau
  DOUBLE PRECISION :: sumTerm
  DOUBLE PRECISION :: pi
  DOUBLE PRECISION :: twoPi
  DOUBLE PRECISION :: correction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Domega
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: A_omega
  COMPLEX(KIND=8) , ALLOCATABLE, DIMENSION(:) :: C_omega
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: opertau
  CHARACTER(len=5) :: funct
  character(len=4) :: tag_proc
  character(len=30) :: tmpfil

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_backFourier : Uninitialized GreenHyboffdiag structure")
  IF ( op%setW .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_backFourier : no G(iw)")
  
  funct="hybri"
  if(present(func)) funct=func
!sui!write(6,*) funct
  inv_beta     = op%inv_beta
  two_invBeta  = 2.d0 * inv_beta
  minusDt      = - op%delta_t
  omegaSamples = op%Wmax
  tauSamples   = op%samples-1
  pi         = ACOS(-1.d0)
  twoPi        = 2.d0 * pi
  pi_invBeta = pi * inv_beta
!sui!write(6,*) "omegaSamples",omegaSamples
  MALLOC(Domega,(1:omegaSamples))
  MALLOC(A_omega,(1:omegaSamples))
  MALLOC(C_omega,(1:omegaSamples))
  IF ( op%rank .EQ. 0 ) THEN
    !DO iflavor1 = 1, op%nflavors
    !  DO iflavor2 = 1, op%nflavors
    !    write(22236,*) "#",iflavor1,iflavor2
    !    do  iomega=1,op%Wmax
    !      write(22236,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%oper_w(iomega,iflavor1,iflavor2))
    !    enddo
    !    write(22236,*) 
    !  ENDDO
    !ENDDO
  ENDIF

  op%oper = 0.d0

  DO iflavor1 = 1, op%nflavors
    DO iflavor2 = 1, op%nflavors
  ! --  compute limit of function G*(i\omega_n)
      if(funct=="hybri") then
        IF ( PRESENT(dvgc) ) THEN
          A = dvgc
        ELSE
          A = AIMAG(op%oper_w(omegaSamples,iflavor1,iflavor2))&! A = \lim_\infty Imag(G(iwn))*(wn)
            *(2.d0*DBLE(omegaSamples)-1.d0) * pi_invBeta
          AA = AIMAG(op%oper_w(omegaSamples-10,iflavor1,iflavor2))&! A = \lim_\infty Imag(G(iwn))*(wn)
            *(2.d0*DBLE(omegaSamples-10)-1.d0) * pi_invBeta
          B = op%oper_w(omegaSamples,iflavor1,iflavor2)&! A = \lim_\infty (G(iwn))*(iwn)
            *(2.d0*DBLE(omegaSamples)-1.d0) *pi_invBeta*cmplx(0.d0,1.d0,kind=8)
          BB = op%oper_w(omegaSamples-10,iflavor1,iflavor2)&! A = \lim_\infty (G(iwn))*(iwn)
            *(2.d0*DBLE(omegaSamples-10)-1.d0) *pi_invBeta*cmplx(0.d0,1.d0,kind=8)
        !sui!write(6,*) "B=",iflavor1,iflavor2,B,BB
        END IF
      else if(iflavor1==iflavor2.and.funct=="green") then
        A = -1.d0
      else if(iflavor1/=iflavor2.and.funct=="green") then
        A = 0.d0
      endif
    !sui!write(6,*) "A=",iflavor1,iflavor2,A,AA
      C=cmplx(-A,0.d0,kind=8)
      if(present(hybri_limit)) then
        if(present(opt_hybri_limit)) then
          if(opt_hybri_limit==1) C= (hybri_limit(iflavor1,iflavor2))
        !sui!write(6,*) "C=                         ",C
        endif
      endif 



  ! --  correction on G(tau=0) is thus
      correction = -C*0.5d0
    
  ! --  built frequency mesh
      Domega = (/ ((2.d0 * DBLE(iomega) - 1.d0)*pi_invbeta, iomega=1, omegaSamples) /)

  ! --  built asymptotic function
      !if(present(hybri_limit)) then
      C_omega = C / (Domega*cmplx(0.d0,1.d0,kind=8))
      !else
      !  A_omega = A / Domega
      !  C_omega=cmplx(0.d0,A_omega,kind=8)
      !endif
      !write(6,*) "AC 1",A_omega(2),C_omega(2)

      IF ( op%rank .EQ. 0 ) THEN
       ! write(236,*) "#",iflavor1,iflavor2
       ! write(237,*) "#",iflavor1,iflavor2
       ! write(2236,*) "#",iflavor1,iflavor2
       ! write(2237,*) "#",iflavor1,iflavor2
       ! write(238,*) "#",iflavor1,iflavor2
       ! do  iomega=1,op%Wmax
       !   write(2236,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%oper_w(iomega,iflavor1,iflavor2))
       !   write(2237,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(op%oper_w(iomega,iflavor1,iflavor2))
       !   write(236,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%oper_w(iomega,iflavor1,iflavor2)-C_omega(iomega))
       !   write(237,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,imag(op%oper_w(iomega,iflavor1,iflavor2)-C_omega(iomega))
       !   write(238,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,C_omega(iomega)
       ! enddo
       ! write(236,*) 
       ! write(237,*) 
       ! write(2236,*) 
       ! write(2237,*) 
       ! write(238,*) 
      ENDIF
      IF ( op%rank .EQ. 1 ) THEN
       ! write(22360,*) "#",iflavor1,iflavor2
       ! do  iomega=1,op%Wmax
       !   write(22360,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%oper_w(iomega,iflavor1,iflavor2))
       ! enddo
       ! write(22360,*) 
      ENDIF
      !open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
      !write(unitnb,*) "#",iflavor1,iflavor2
      !do  iomega=1,op%Wmax
      !  write(unitnb,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%oper_w(iomega,iflavor1,iflavor2))
      !enddo
      !write(unitnb,*) 
  ! --  built time mesh
      IF (op%have_MPI .EQV. .TRUE.) THEN
        delta = tauSamples / op%size
        residu = tauSamples - op%size*delta
        IF ( op%rank .LT. op%size - residu ) THEN
          tauBegin = 1 + op%rank*delta
          tauEnd   = (op%rank + 1)*delta
        ELSE
!          tauBegin = (op%size-residu)*delta + 1 + (op%rank-op%size+residu)*(delta+1)
          tauBegin = 1 + op%rank*(delta + 1) -op%size + residu
          tauEnd = tauBegin + delta
        END IF
        MALLOC(counts,(1:op%size))
        MALLOC(displs,(1:op%size))
        counts = (/ (delta, iTau=1, op%size-residu), &
                    (delta+1, iTau=op%size-residu+1, op%size) /)
        displs(1)=0
        DO iTau = 2, op%size
          displs(iTau) = displs(iTau-1) + counts (iTau-1)
        END DO
      ELSE
        tauBegin = 1
        tauEnd   = tauSamples
      END IF
      MALLOC(opertau,(1:tauSamples+1))
      do iomega=1,omegaSamples
       ! write(6,*) iomega, imag(op%oper_w(iomega,iflavor1,iflavor2)), A_omega(iomega) ,"#diff"
      enddo
      unitnb=70000+op%rank
      call int2char4(op%rank,tag_proc)
      tmpfil = 'counts'//tag_proc
     ! open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
     ! write(unitnb,*) "#",iflavor1,iflavor2
     ! do  itau=1,op%size
     ! write(unitnb,*)  itau,counts(itau),displs(itau)
     ! enddo
     ! write(unitnb,*) 

      unitnb=10000+op%rank
      call int2char4(op%rank,tag_proc)
      tmpfil = 'oper_w'//tag_proc
     ! open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
     ! write(unitnb,*) "#",iflavor1,iflavor2,C
     ! ! C_omega et oper_w differents Domega identique. Est ce du a des
     !! ! diago differentes   pour chaque procs dans qmc_prep_ctqmc
     !! do  iomega=1,op%Wmax
     !! write(unitnb,*)  (2.d0*DBLE(iomega)-1.d0) * pi_invBeta,real(op%oper_w(iomega,iflavor1,iflavor2)),C_omega(iomega),Domega(iomega)
     ! enddo
     ! write(unitnb,*) 

     ! unitnb=40000+op%rank
     ! unitnb1=50000+op%rank
     ! call int2char4(op%rank,tag_proc)
     ! tmpfil = 'tauend'//tag_proc
     ! open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
     ! tmpfil = 'taubegin'//tag_proc
     ! open (unit=unitnb1,file=trim(tmpfil),status='unknown',form='formatted')
     ! write(unitnb,*) "#",iflavor1,iflavor2
     ! write(unitnb1,*) "#",iflavor1,iflavor2

  ! -- compute Fourier transformation
      opertau=0.d0
      DO itau = tauBegin, tauEnd
      !DO itau = max(tauBegin-1,1), tauEnd
        minusTau = DBLE(itau -1) * minusDt
        DO iomega = 1, omegaSamples
          omegaa         = Domega(iomega)
          minusOmegaTau = MOD(omegaa*minusTau, TwoPi)
          sumTerm       = REAL(( op%oper_w(iomega,iflavor1,iflavor2) &
                          -  C_omega(iomega) ) &
                          !- CMPLX(0.d0, A_omega(iomega),8) ) &
                          * EXP( CMPLX(0.d0, minusOmegaTau, 8)))
          opertau(itau)  = opertau(itau) + sumTerm
!         Domega et minusomegatau identique MAIS oper_w different
            !write(unitnb,*) iomega,Domega(iomega),real(C_omega(iomega)),imag(C_omega(iomega))
          !if(itau==tauend) then
          !  write(unitnb,*) iomega, sumTerm,opertau(itau),minusOmegaTau,op%oper_w(iomega,iflavor1,iflavor2),Domega(iomega)
          !endif 
          !if(itau==max(tauBegin-1,1)) then
          !  write(unitnb1,*) iomega, sumTerm,opertau(itau),minusOmegaTau,op%oper_w(iomega,iflavor1,iflavor2),Domega(iomega)
          !endif 

        END DO
         ! if(itau==tauEnd) write(unitnb,*) 
         ! if(itau==max(tauBegin-1,1)) write(unitnb1,*) 
!        if(iflavor1==iflavor2) then
          opertau(itau) = correction + two_invBeta*opertau(itau)
         ! if(itau==tauend) then
         !   write(unitnb,*) "final", opertau(itau),correction
         ! endif 
         ! if(itau==max(tauBegin-1,1)) then
         !   write(unitnb1,*) "final",opertau(itau),correction
         ! endif 
             !write(66666,*) itau, opertau(itau),correction
!        else
!          opertau(itau) =              &
!             two_invBeta*opertau(itau)
!        endif
      END DO
      !unitnb=20000+op%rank
      !call int2char4(op%rank,tag_proc)
      !tmpfil = 'opertau'//tag_proc
      !open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
      !write(unitnb,*) "#",iflavor1,iflavor2,tauBegin,tauEnd
      !!do  itau=tauBegin, tauEnd
      !do  itau=1,tauSamples
      !  write(unitnb,*)    itau,opertau(itau)
      !enddo
      !write(unitnb,*) 
      !opertau(tauBegin-1)=0.d0
      !opertau(tauEnd+1)=0.d0
             !write(66666,*)
  
  ! -- Gather
      IF ( op%have_MPI .EQV. .TRUE. ) THEN
! rassembler les resultats
#ifdef HAVE_MPI
        CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                          opertau, counts, displs, &
                          MPI_DOUBLE_PRECISION, op%MY_COMM, residu)
#endif
        FREE(counts)
        FREE(displs)
      END IF
     ! unitnb=30000+op%rank
     ! call int2char4(op%rank,tag_proc)
     ! tmpfil = 'opertau_MPI_'//tag_proc
     ! open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
     ! write(unitnb,*) "#",iflavor1,iflavor2
     ! do  itau=tauBegin, tauEnd
     !   write(unitnb,*)    itau,opertau(itau)
     ! enddo
     ! write(unitnb,*) 

  ! -- Add correction for discontinuity.
!      if(iflavor1==iflavor2) then
        !G(0+)-G(0-)=G(0+)+G(beta-)=A
        opertau(tauSamples+1) = -real(C) - opertau(1)
      !sui!write(6,*) "BackFourier",opertau(tauSamples+1),opertau(1),real(C)

        op%setT = .TRUE.
!      endif
      op%oper(:,iflavor1,iflavor2)=opertau(:)
      FREE(opertau)
    END DO ! iflavor2
  END DO ! iflavor1
  ! -- End loop over flavors.

  FREE(Domega)
  FREE(A_omega)
  FREE(C_omega)
  close(236)
  close(237)

END SUBROUTINE GreenHyboffdiag_backFourier
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_forFourier
!! NAME
!!  GreenHyboffdiag_forFourier
!!
!! FUNCTION
!!  perform forward fourier transform
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  Wmax=linear maximum frequency
!!
!! OUTPUT
!!  Gomega=Results for omega frequencies
!!  omega=ask frequencies
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

SUBROUTINE GreenHyboffdiag_forFourier(op, Gomega, omega, Wmax)
!Arguments ------------------------------------

#ifdef HAVE_MPI1
include 'mpif.h'
#endif
  TYPE(GreenHyboffdiag)             , INTENT(INOUT) :: op
  COMPLEX(KIND=8), DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: Gomega  ! INOUT for MPI
  COMPLEX(KIND=8), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: omega  
  INTEGER                 , OPTIONAL, INTENT(IN   ) :: Wmax   
  INTEGER :: i
  INTEGER :: j
  INTEGER :: iflavor1
  INTEGER :: iflavor2
  INTEGER :: nflavors
  INTEGER :: L
  INTEGER :: Lspline
  INTEGER :: Nom
  INTEGER :: omegaBegin
  INTEGER :: omegaEnd
  INTEGER :: deltaw
  INTEGER :: residu
  INTEGER, ALLOCATABLE, DIMENSION(:) :: counts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: displs
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: tau
  DOUBLE PRECISION :: delta
  DOUBLE PRECISION :: deltabis
  DOUBLE PRECISION :: inv_delta
  DOUBLE PRECISION :: inv_delta2
  DOUBLE PRECISION :: omdeltabis
  DOUBLE PRECISION :: tmp
  DOUBLE PRECISION :: xpi
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  diag
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  diagL
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  lastR
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  lastC
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  XM
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  X2
  DOUBLE PRECISION :: iw
  COMPLEX(KIND=8) :: iwtau
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: Gwtmp  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: omegatmp
  nflavors=op%nflavors

!sui!write(6,*) " Fourier transformation begin"

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_forFourier : Uninitialized GreenHyboffdiag structure")
  IF ( op%setT .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_forFourier : no G(tau)")
  !write(6,*) "op%setMk=", op%setMk
  IF ( op%setMk .NE. 2*nflavors*nflavors ) &
    CALL WARNALL("GreenHyboffdiag_forFourier : green does not have moments    ")

  L  = op%samples

  xpi=acos(-1.d0)                !!! XPI=PI
  beta = op%beta
  Nom  = op%Wmax
  IF ( PRESENT(Gomega) ) THEN
    Nom = SIZE(Gomega,1)    
    !IF ( op%rank .EQ. 0 ) &
      !!write(6,*) "size Gomega", Nom
  END IF
  IF ( PRESENT(omega) ) THEN
    IF ( PRESENT(Gomega) .AND. SIZE(omega) .NE. Nom ) THEN
      CALL ERROR("GreenHyboffdiag_forFourier : sizes mismatch              ")               
    !ELSE 
      !Nom = SIZE(omega)
    END IF 
  END IF
  IF ( .NOT. PRESENT(Gomega) .AND. .NOT. PRESENT(omega) ) THEN
    IF ( PRESENT(Wmax) ) THEN
      Nom=Wmax
    ELSE
      CALL ERROR("GreenHyboffdiag_forFourier : Missing argument Wmax")
    END IF
  END IF

  !!IF ( ALLOCATED(op%oper_w) ) THEN
  !!  IF ( SIZE(op%oper_w,1) .NE. Nom ) THEN
  !!    FREE(op%oper_w)
  !!    MALLOC(op%oper_w,(1:Nom,nflavors,nflavors))
  !!  END IF
  !!ELSE
  !!  MALLOC(op%oper_w,(1:Nom,nflavors,nflavors))
  !!END IF

  !!write(6,*) "PRESENT(GOMEGA)", PRESENT(GOMEGA)
  !!write(6,*) "PRESENT(OMEGA)", PRESENT(OMEGA)
  !call flush(6)

  delta=op%delta_t
  inv_delta = op%inv_dt
  inv_delta2 = inv_delta*inv_delta
 
  MALLOC(diagL,(L-1))
  MALLOC(lastR,(L-1))
  MALLOC(diag,(L))
  MALLOC(lastC,(L-1))

!(cf Stoer) for the spline interpolation : 
! second derivatives XM solution of A*XM=B.
!A=(2.4.2.11) of Stoer&Bulirsch + 2 limit conditions
!The LU decomposition of A is known explicitly; 

  diag (1) = 4.d0 ! 1.d0 *4.d0 factor 4 added for conditionning
  diagL(1) = 0.25d0 !1.d0/4.d0
  lastR(1) = -0.5d0 ! -2.d0/4.d0
  lastC(1) = 4.d0 ! 1.d0*4.d0
  diag (2) = 4.d0
  diagL(2) = 0.25d0
  lastR(2) = -0.25d0
  lastC(2) = -1.d0

  DO i = 3, L-2
    tmp = 4.d0 - diagL(i-1)
    diagL(i) = 1.d0 / tmp
  END DO
  DO i = 3, L-2
    diag (i) = 1.d0 / diagL(i)
    lastR(i) = -(lastR(i-1)*diagL(i))
    lastC(i) = -(lastC(i-1)*diagL(i-1))
  END DO
  
  tmp = 1.d0/diag(L-2)
  diag (L-1) = 4.d0 - tmp
  lastR(L-1) = (1.d0 - lastR(L-2))/ diag(L-1)
  !diagL(L-1) = lastR(L-1)
  diagL(L-1) = 0.d0 ! for the Lq=B resolution
  !lastC(L-1) = 1.d0 - lastC(L-2)*diagL(L-1) ! equivalent to the next line
  lastC(L-1) = 1.d0 - (lastC(L-2)*lastR(L-1)) ! True value
  diag (L  ) = 2.d0! - DOT_PRODUCT( lastR , lastC )
  tmp = 0.d0
  DO i = 1, L-1
    tmp = tmp + lastR(i)*lastC(i)
  END DO
  diag (L  ) = diag (L  ) - tmp
  lastC(L-1) = lastC(L-1)-1.d0 ! 1 is removed for the u.XM=q resolution

  MALLOC(XM,(L))
  MALLOC(Gwtmp,(1:Nom))

  Lspline = L-1
  MALLOC(X2,(1:Lspline+1)) ! We impose L = Nom

  IF ( op%have_MPI .EQV. .TRUE. ) THEN  
    deltaw = Nom / op%size
    residu = Nom - op%size*deltaw
    IF ( op%rank .LT. op%size - residu ) THEN
      omegaBegin = 1 + op%rank*deltaw
      omegaEnd   = (op%rank + 1)*deltaw
    ELSE
  !    tauBegin = (op%size-residu)*deltaw + 1 + (op%rank-op%size+residu)*(deltaw+1)
      omegaBegin = 1 + op%rank*(deltaw + 1) -op%size + residu
      omegaEnd = omegaBegin + deltaw
    END IF
    MALLOC(counts,(1:op%size))
    MALLOC(displs,(1:op%size))
    counts = (/ (deltaw, i=1, op%size-residu), &
                (deltaw+1, i=op%size-residu+1, op%size) /)
    displs(1)=0
    DO i = 2, op%size
      displs(i) = displs(i-1) + counts (i-1)
    END DO
  ELSE
    omegaBegin = 1
    omegaEnd   = Nom 
  END IF

!  op%Mk(iflavor1,iflavor2,1) = 0.d0
!  DO iflavor1 = 1, nflavors
!    op%Mk(iflavor1,iflavor1,1) = -1.d0
!  ENDDO
!  op%Mk(:,:,3) = 0.d0

  MALLOC(omegatmp,(omegaBegin:omegaEnd))
  IF ( PRESENT(omega) ) THEN
    omegatmp(omegaBegin:omegaEnd) = (/ (AIMAG(omega(i)),i=omegaBegin,omegaEnd) /)
  ELSE
    omegatmp(omegaBegin:omegaEnd) = (/ ((((2.d0*DBLE(i)-1.d0)*xpi)/Beta), i=omegaBegin,omegaEnd) /)
  END IF

  DO iflavor1 = 1, nflavors
    DO iflavor2 = 1, nflavors
   ! write(6,*) "   Moments:",op%Mk(iflavor1,iflavor2,:),iflavor1,iflavor2
  
! construct the B vector from A.Xm=B
      XM(1) = 4.d0*op%Mk(iflavor1,iflavor2,3)
      XM(L) = (6.d0 * inv_delta) * ( op%Mk(iflavor1,iflavor2,2) - ( &
        (op%oper(2,iflavor1,iflavor2)-op%oper(1,iflavor1,iflavor2)) + &
        (op%oper(L,iflavor1,iflavor2)-op%oper(L-1,iflavor1,iflavor2)) ) * inv_delta )
!    built generic second derivative of oper
!sui!write(6,*)  "XM 1 L",XM(1),XM(L),op%Mk(iflavor1,iflavor2,2),op%Mk(iflavor1,iflavor2,3)
      DO i = 2, L-1
        XM(i) = (6.d0 * inv_delta2) * ( (op%oper(i+1,iflavor1,iflavor2) &
          - 2.d0 * op%oper(i,iflavor1,iflavor2)) &
          +        op%oper(i-1,iflavor1,iflavor2) )
    !sui!write(6,*) "XM",i,XM(i),op%oper(i,iflavor1,iflavor2)
      END DO

! Find second derivatives XM: Solve the system
! SOLVING Lq= XM 
!  q = XM 
      do j=1,L-1
          XM(j+1)=XM(j+1)-(diagL(j)*XM(j))
          XM(L)  =XM(L)  -(lastR(j)*XM(j))
      end do


! SOLVING U.XM=q 
!  XM = q
      do j=L-1,2,-1
       XM(j+1)  = XM(j+1) / diag(j+1)
       XM(j)= (XM(j)-(XM(L)*lastC(j)))-XM(j+1)
      end do
      XM(2)  = XM(2) / diag(2)
      XM(1) = (XM(1)-XM(L)*lastC(1)) / diag(1)



      !Construct L2 second derivative from known derivatives XM
      deltabis = beta / DBLE(Lspline)
      DO i = 1, Lspline
        tau = deltabis * DBLE(i-1)
        j = ((L-1)*(i-1))/Lspline + 1!INT(tau * inv_delta) + 1
        X2(i) = inv_delta * ( XM(j)*(DBLE(j)*delta - tau ) + XM(j+1)*(tau - DBLE(j-1)*delta) )
      END DO
      X2(Lspline+1) = XM(L)


       DO i = omegaBegin, omegaEnd
         iw = omegatmp(i)
         omdeltabis = iw*deltabis
         Gwtmp(i)=CMPLX(0.d0,0.d0,8)
         DO j=2, Lspline ! We impose  L+1 = Nom
           iwtau = CMPLX(0.d0,omdeltabis*DBLE(j-1),8)
           Gwtmp(i) = Gwtmp(i) + EXP(iwtau) * CMPLX((X2(j+1) + X2(j-1))-2.d0*X2(j),0.d0,8)
           !write(6,*) "ww",i,j,Gwtmp(i),X2(j),iwtau
         END DO
         Gwtmp(i) = Gwtmp(i)/CMPLX(((iw*iw)*(iw*iw)*deltabis),0.d0,8) &
         + CMPLX( ( ((X2(2)-X2(1))+(X2(Lspline+1)-X2(Lspline)))/((iw*iw)*deltabis) -op%Mk(iflavor1,iflavor2,2) ) &
         /(iw*iw) , (op%Mk(iflavor1,iflavor2,1)-op%Mk(iflavor1,iflavor2,3)/(iw*iw))/iw , 8) 
                   !+ CMPLX( (X2(2)-X2(1))+(X2(Lspline+1)-X2(Lspline)), 0.d0, 8 ) ) &
                   !   / (((iw*iw)*(iw*iw))*CMPLX(deltabis,0.d0,8)) &
                   !- CMPLX(op%Mk(1),0.d0,8)/iw  &
                   !+ CMPLX(op%Mk(2),0.d0,8)/(iw*iw) &
                   !- CMPLX(op%Mk(3),0.d0,8)/((iw*iw)*iw) 
         !IF ( op%rank .EQ. 0 )  write(12819,*) iw,gwtmp(i)
       END DO
       !call flush(12819)
       IF ( op%have_MPI .EQV. .TRUE. ) THEN
#ifdef HAVE_MPI
        CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                          Gwtmp  , counts, displs, &
                          MPI_DOUBLE_COMPLEX, op%MY_COMM, residu)
#endif
      END IF
      IF ( PRESENT(Gomega) ) THEN
        Gomega(:,iflavor1,iflavor2) = Gwtmp(:)
      END IF
      op%setW = .TRUE.
    ENDDO ! iflavor1
  ENDDO ! iflavor2
  !!op%oper_w=Gomega
  do iflavor1=1,nflavors
    !sui!write(6,*)  iflavor1
      do i=1,Nom
      !write(6,*) "w",i,op%oper_w(i,iflavor1,iflavor1)
      enddo
  enddo
  FREE(Gwtmp)
  FREE(diagL)
  FREE(lastR)
  FREE(diag)
  FREE(lastC)
  FREE(XM)
  FREE(omegatmp)
  FREE(X2)
  FREE(counts)
  FREE(displs)

END SUBROUTINE GreenHyboffdiag_forFourier
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_print
!! NAME
!!  GreenHyboffdiag_print
!!
!! FUNCTION
!!  print Green function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
!!  ostream=file stream
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

SUBROUTINE GreenHyboffdiag_print(op, ostream)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag), INTENT(IN) :: op
  INTEGER, OPTIONAL , INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                        :: ostream_val
  INTEGER                        :: isample
  INTEGER                        :: samples
  INTEGER                        :: iflavor1
  INTEGER                        :: iflavor2
  

  IF ( op%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyboffdiag_print : green op%operator not set              ")

  IF ( PRESENT(ostream) ) THEN
    ostream_val = ostream
  ELSE
    ostream_val = 66
    OPEN(UNIT=ostream_val,FILE="Green.dat")
  END IF

  samples =  op%samples 

  DO iflavor1=1,op%nflavors
    DO iflavor2=1,op%nflavors
    WRITE(ostream_val,'(a,i3,a,i3,a)')  "## (iflavor1,iflavor2)= (", iflavor1,",",iflavor2,")"
      DO isample = 1, samples
      WRITE(ostream_val,*) DBLE(isample-1)*op%delta_t, op%oper(isample,iflavor1,iflavor2)
      END DO
      WRITE(ostream_val,*) 
    END DO
  END DO

  IF ( .NOT. PRESENT(ostream) ) &
    CLOSE(ostream_val)
END SUBROUTINE GreenHyboffdiag_print
!!***

!!****f* ABINIT/m_GreenHyboffdiag/GreenHyboffdiag_destroy
!! NAME
!!  GreenHyboffdiag_destroy
!!
!! FUNCTION
!!  destroy green function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=Green
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

SUBROUTINE GreenHyboffdiag_destroy(op)

!Arguments ------------------------------------
  TYPE(GreenHyboffdiag), INTENT(INOUT) :: op
  INTEGER :: iflavor,iflavorbis

  op%set          = .FALSE.
  op%setT         = .FALSE.
  op%setW         = .FALSE.
  op%samples      = 0
  op%measurements = 0
  op%beta         = 0.d0
  op%inv_beta     = 0.d0
  op%inv_dt       = 0.d0
  op%delta_t      = 0.d0
  CALL VectorInt_destroy(op%index_old)
  CALL Vector_destroy(op%oper_old)
  do iflavor=1,op%nflavors
    do iflavorbis=1,op%nflavors
     !sui!write(6,*) "test",iflavor,iflavorbis
      CALL MapHyb_destroy(op%map(iflavor,iflavorbis))
    enddo
  enddo
  DT_FREEIF(op%map)
  FREEIF(op%oper)
  FREEIF(op%Mk)
  FREEIF(op%oper_w)
  FREEIF(op%oper_w_old)
  FREEIF(op%omega)
END SUBROUTINE GreenHyboffdiag_destroy

!!***
! This routine contains direct and inverse fourier transformation
! It is a modification of a routine of the GNU GPL 
! code available on http://dmft.rutgers.edu/ and
! described in the RMP 2006 paper written by
! G.Kotliar, S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti
!=======+=========+=========+=========+=========+=========+=========+=$
!       TYPE   : SUBROUTINE
!       PROGRAM: nfourier3
!       PURPOSE: fourier-transform the natural-spline interpolation
!                of function Green(tau) 
!                calculate function Green(omega)
!       I/O    :
!       VERSION: 2-16-92
!                29-Nov-95 removal of minimal bug concerning 
!                          DIMENSION of rindata
!       COMMENT: cf J. Stoer R. Bulirsch, Introduction to numerical
!                analysis (Springer, New York, 1980)
!=======+=========+=========+=========+=========+=========+=========+=$
!
      SUBROUTINE nfourier3(rindata,coutdata,lflag,Iwmax,L,Beta,AA,c1,c2,c3)

! use m_profiling
!       include 'param.dat'
!       use defs_basis
!Arguments ------------------------------------
       integer, intent(in) :: Iwmax,L
       logical, intent(in)  :: lflag
       real*8, intent(in) :: beta,AA,c1,c2,c3
       real*8, intent(in) :: rindata(L)
       complex*16, intent(out) :: coutdata(Iwmax+1)
!Local variables ------------------------------
       integer    :: i,j,k,p
       real*8     :: rincopy(L+1),a(L),b(L),c(L),d(L),u(L+1), q(L+1),XM(L+1)
       complex*16 :: cdummy,explus,ex,j_dpc
       real*8     :: one,two,zero,three,six,tau,xpi,delta,om !wn,
       complex*16     :: czero
!***********************************************

       ABI_UNUSED((/aa, c3/))
       ABI_UNUSED((/lflag/))
       czero=cmplx(0.d0,0.d0)
       zero=0.d0
       one=1.d0
       two=2.d0
       three=3.d0
       six=6.d0
       j_dpc=dcmplx(0.d0,1.d0)
       xpi = ACOS(-One)
       delta = beta/float(L)
       ! c2 devrait etre nul en symetrie particule trou.
       DO i = 1,L
          tau=beta/dble(L)*dble(i-1)
          !worksrincopy(i) = rindata(i)-c1/two!-c2/4.d0*(-Beta+2*tau)+c3/4.d0*(beta*tau-tau*tau)
          rincopy(i) = rindata(i)-c1/two-c2/4.d0*(-Beta+2.d0*tau)!+c3/4.d0*(beta*tau-tau*tau)
         ! rincopy(i) = -c2/4.d0*(-Beta+2*tau)
       !   write(99,*) i,rindata(i)
       !   write(98,*) i,tau,rincopy(i),rindata(i),rindata(i)-c1/two,-c1/two,-c2/4.d0*(-Beta+2.d0*tau)
       !   write(97,*) i,(-Beta+two*tau),c2/4.d0,2.d0*tau,-c2/4.d0*(-Beta+2.d0*tau)
       ENDDO
       !   write(99,*) 
       !   write(98,*) 
!       if(lflag) then
!         rincopy(L+1) = AA-rindata(1)
!       else 
         rincopy(L+1) = -rindata(1)
!       endif
       !DO i = 1,L+1
       !   write(999,*) i,rincopy(i)
       !ENDDO
       !write(6,*) lflag,Iwmax,L,Beta,delta
!       Three = Two+One
!       six = Two*Three
     
!c
!c     spline interpolation:  the spline is given by
!c     G(tau) = a(i) + b(i) (tau-tau_i) + c(i) ( )^2 + d(i) ( )^3
!c     The following formulas are taken directly from  Stoer and
!c     Bulirsch p. 102
!c
       q(1) = Zero
       u(1) = Zero
       DO k = 2,L
          p = q(k-1)/Two+Two
          q(k)=-One/Two/p
!     this is equation 2.4.2.10 or Bulirsch for dn. here uk=dn
          u(k)=Three/delta**2*(rincopy(k+1)+rincopy(k-1)-Two*rincopy(k))
          u(k)=(u(k)-u(k-1)/Two)/p
       ENDDO
       XM(L+1) = 0
       DO k = L,1,-1
          XM(k) = q(k)*XM(k+1)+u(k)
       ENDDO
!c
!c     The following formulas are taken directly from  Stoer and
!c     Bulirsch p. 98 second edition.
!c     a b c d are the spline coefficients.
!c     XM(j) is the second derivative at node j
!c     

       DO j = 1, L
          a(j) = rincopy(j)
          c(j) = XM(j)/Two
          b(j) = (rincopy(j+1)-rincopy(j))/delta - &
     &       (Two*XM(j)+XM(j+1))*delta/6.
          d(j) = (XM(j+1)-XM(j))/(6.*delta)
       ENDDO

!c
!c     The Spline multiplied by the exponential can now be exlicitely
!c     integrated. The following formulas were obtained using
!c     MATHEMATICA
!c
        DO i = 0,Iwmax
           om = (Two*(i)+One)*xpi/Beta
           coutdata(i+1) = czero
           DO j = 1,L
              cdummy = j_dpc*om*delta*j
              explus = exp(cdummy)
              cdummy = j_dpc*om*delta*(j-1)
              ex = exp(cdummy)
              coutdata(i+1) = coutdata(i+1) + explus*(&
     &         ( -six* d(j) )/om**4 + &
     &         ( Two*j_dpc*c(j) + six*delta*j_dpc*d(j)  )/om**3 +&
     &         ( b(j)+ Two*delta*c(j)+ three*delta**2*d(j) )/om**2 +&
     &         (- j_dpc*a(j) - delta*j_dpc*b(j) - delta**2*j_dpc*c(j) -&
     &         delta**3*j_dpc*d(j))/om)
 
              coutdata(i+1) = coutdata(i+1) + ex*(&
     &        six*d(j)/om**4 - Two*j_dpc*c(j)/om**3 &
     &        -b(j)/om**2 + j_dpc*a(j)/om)
           ENDDO
          !write(100,*) i,real(coutdata(i+1)),imag(coutdata(i+1))
        ENDDO
!        DO i = 0,Iwmax
!          wn=3.1415926/beta*dble(2*i+1)
!          write(101,*) wn,real(coutdata(i+1)),imag(coutdata(i+1)),c1/wn,-c2/(wn*wn),-c3/(wn*wn*wn)
!          coutdata(i+1)=coutdata(i+1)+cmplx(-c2/(wn*wn),c1/(wn))!-c3/(wn*wn*wn))
! !works         coutdata(i+1)=coutdata(i+1)+cmplx(0.d0,c1/(wn))
!          write(103,*) wn,real(coutdata(i+1)),imag(coutdata(i+1)),c1/wn,-c2/(wn*wn),-c3/(wn*wn*wn)
!        ENDDO
!          write(101,*)
!          write(103,*)
        end subroutine nfourier3

END MODULE m_GreenHyboffdiag
!!***

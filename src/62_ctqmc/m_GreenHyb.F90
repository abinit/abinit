
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_GreenHyb
!! NAME
!!  m_GreenHyb
!! 
!! FUNCTION 
!!  Manage a green function for one orbital
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
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
MODULE m_GreenHyb
USE m_Global
USE m_MatrixHyb
USE m_Vector
USE m_VectorInt
USE m_ListCdagC
USE m_MapHyb
#ifdef HAVE_MPI2
USE mpi
#endif
IMPLICIT NONE

!!***

PRIVATE

!!****t* m_GreenHyb/GreenHyb
!! NAME
!!  GreenHyb
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

TYPE, PUBLIC :: GreenHyb
  LOGICAL _PRIVATE :: set = .FALSE.
  LOGICAL _PRIVATE :: setT = .FALSE.
  LOGICAL _PRIVATE :: setW = .FALSE.
  LOGICAL _PRIVATE :: have_MPI = .FALSE.
  INTEGER _PRIVATE :: setMk = 0
  INTEGER _PRIVATE :: samples
  INTEGER _PRIVATE :: measurements    
  INTEGER          :: factor
  INTEGER _PRIVATE :: MY_COMM
  INTEGER _PRIVATE :: size
  INTEGER _PRIVATE :: rank
  INTEGER _PRIVATE :: Wmax
  INTEGER _PRIVATE :: iTech
  DOUBLE PRECISION _PRIVATE :: beta
  DOUBLE PRECISION _PRIVATE :: inv_beta
  DOUBLE PRECISION _PRIVATE :: delta_t
  DOUBLE PRECISION _PRIVATE :: inv_dt
  DOUBLE PRECISION, ALLOCATABLE , DIMENSION(:)            :: oper
  DOUBLE PRECISION, ALLOCATABLE , DIMENSION(:)   _PRIVATE :: omega
  DOUBLE PRECISION              , DIMENSION(1:3) _PRIVATE :: Mk
  COMPLEX(KIND=8)  , ALLOCATABLE, DIMENSION(:)   _PRIVATE :: oper_w 
  COMPLEX(KIND=8)  , ALLOCATABLE, DIMENSION(:)   _PRIVATE :: oper_w_old
  TYPE(MapHyb)          :: this
END TYPE GreenHyb
!!***

PUBLIC :: GreenHyb_init
PUBLIC :: GreenHyb_clear
PUBLIC :: GreenHyb_reset
PUBLIC :: GreenHyb_setOperW
PUBLIC :: GreenHyb_measHybrid
PUBLIC :: GreenHyb_getHybrid
PUBLIC :: GreenHyb_setN
PUBLIC :: GreenHyb_setMuD1
PUBLIC :: GreenHyb_setMoments
PUBLIC :: GreenHyb_backFourier
PUBLIC :: GreenHyb_forFourier
PUBLIC :: GreenHyb_print
PUBLIC :: GreenHyb_destroy

CONTAINS
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_init
!! NAME
!!  GreenHyb_init
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
!!  this=Green
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

SUBROUTINE GreenHyb_init(this, samples, beta,iTech,MY_COMM)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(GreenHyb)     , INTENT(INOUT) :: this
  INTEGER         , INTENT(IN   ) :: samples
  DOUBLE PRECISION, INTENT(IN   ) :: beta
  !INTEGER          , INTENT(IN   ) :: Wmax
  INTEGER, OPTIONAL, INTENT(IN   ) :: iTech
  INTEGER, OPTIONAL, INTENT(IN   ) :: MY_COMM
!Local variables ------------------------------
  INTEGER                         :: sp1
  DOUBLE PRECISION                :: dt
#ifdef HAVE_MPI
  INTEGER          :: ierr
#endif

  IF ( PRESENT(MY_COMM)) THEN
#ifdef HAVE_MPI
    this%have_MPI = .TRUE.
    this%MY_COMM = MY_COMM
    CALL MPI_Comm_rank(this%MY_COMM, this%rank, ierr)
    CALL MPI_Comm_size(this%MY_COMM, this%size, ierr)
#else
    CALL WARN("GreenHyb_init : MPI is not used                                    ")
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

  sp1             = samples + 1
  this%samples      = sp1
  this%measurements = 0
  this%beta         = beta
  this%inv_beta     = 1.d0 / beta
  this%inv_dt       = DBLE(samples) * this%inv_beta
  dt              = 1.d0 / this%inv_dt
  this%delta_t      = dt
  !this%Wmax         = Wmax
  this%Wmax         = -1
  FREEIF(this%oper)
  MALLOC(this%oper,(1:sp1))
  ! If we want to measure in frequences
  ! let assume we first have "samples" frequences
  IF ( PRESENT(iTech) ) THEN
    this%iTech = iTech
    SELECT CASE (this%iTech)
    CASE (GREENHYB_TAU)  ! omega
      this%iTech = GREENHYB_TAU
    CASE (GREENHYB_OMEGA)  ! omega
      this%Wmax = samples
      FREEIF(this%oper_w)
      MALLOC(this%oper_w,(1:this%Wmax))
      FREEIF(this%oper_w_old)
      MALLOC(this%oper_w_old,(1:this%Wmax))
      this%oper_w     = CMPLX(0.d0,0.d0,8)
      this%oper_w_old = CMPLX(0.d0,0.d0,8)
      FREEIF(this%omega)
      MALLOC(this%omega,(1:this%Wmax))
      this%omega = (/ ((2.d0 * DBLE(sp1) - 1.d0)*ACOS(-1.d0)*this%inv_beta, sp1=1, this%Wmax) /)
    END SELECT
  ELSE
    this%iTech = GREENHYB_TAU
  END IF
  ! end if
  !CALL Vector_init(this%oper_old,10000)
  !CALL VectorInt_init(this%index_old,10000)
  CALL MapHyb_init(this%this,10000)

  this%oper       = 0.d0
  this%set        = .TRUE.
  this%factor     = 1
  this%setMk      = 0
  this%Mk         = 0.d0
END SUBROUTINE GreenHyb_init
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_clear
!! NAME
!!  GreenHyb_clear
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
!!  this=Green
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

SUBROUTINE GreenHyb_clear(this)

!Arguments ------------------------------------
  TYPE(GreenHyb)     , INTENT(INOUT) :: this

  !CALL Vector_clear(this%oper_old)
  !CALL VectorInt_clear(this%index_old)
  CALL MapHyb_clear(this%this)
  this%measurements = 0
  IF ( ALLOCATED(this%oper) ) &
  this%oper         = 0.d0
  IF ( this%iTech .EQ. GREENHYB_OMEGA ) THEN
    IF ( ALLOCATED(this%oper_w) ) &
    this%oper_w       = CMPLX(0.d0,0.d0,8)
    IF ( ALLOCATED(this%oper_w_old) ) &
    this%oper_w_old   = CMPLX(0.d0,0.d0,8)
  END IF
  this%factor       = 0
END SUBROUTINE GreenHyb_clear
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_reset
!! NAME
!!  GreenHyb_reset
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
!!  this=Green
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

SUBROUTINE GreenHyb_reset(this)

!Arguments ------------------------------------
  TYPE(GreenHyb)     , INTENT(INOUT) :: this

  CALL GreenHyb_clear(this)
  this%setMk        = 0
  this%Mk           = 0.d0
  this%setT         = .FALSE.
  this%setW         = .FALSE.
END SUBROUTINE GreenHyb_reset
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_setOperW
!! NAME
!!  GreenHyb_setOperW
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
!!  this=Green
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

SUBROUTINE GreenHyb_setOperW(this, Gomega)

!Arguments ------------------------------------
  TYPE(GreenHyb)          , INTENT(INOUT) :: this
  COMPLEX(KIND=8), DIMENSION(:), INTENT(IN   ) :: Gomega
!Loval variables ------------------------------
  INTEGER :: tail

  tail = SIZE(Gomega)
  IF ( .NOT. this%set ) &
    CALL ERROR("GreenHyb_setOperW : Uninitialized GreenHyb structure")
  IF ( ALLOCATED(this%oper_w) ) THEN
    IF ( SIZE(this%oper_w) .NE. tail ) THEN
      FREE(this%oper_w)
      MALLOC(this%oper_w,(1:tail))
    END IF
  ELSE
    MALLOC(this%oper_w,(1:tail))
  END IF
  this%oper_w(:) = Gomega(:)
  this%Wmax = tail
  this%setW = .TRUE.
END SUBROUTINE GreenHyb_setOperW
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_measHybrid
!! NAME
!!  GreenHyb_measHybrid
!!
!! FUNCTION
!!  Measure Green function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=Green
!!  Mthis=M this for the current flavor
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

SUBROUTINE GreenHyb_measHybrid(this, Mthis, ListCdagC_1, updated)

!Arguments ------------------------------------
  TYPE(GreenHyb)    , INTENT(INOUT) :: this
  TYPE(MatrixHyb)   , INTENT(IN   ) :: Mthis
  TYPE(ListcdagC), INTENT(IN   ) :: ListCdagC_1
  LOGICAL        , INTENT(IN   ) :: updated
!Local variables ------------------------------
  INTEGER                        :: iC
  INTEGER                        :: iCdag
  INTEGER                        :: tail
  !INTEGER                        :: index
  INTEGER                        :: idx_old
  INTEGER                        :: old_size
  INTEGER                        :: omegaSamples
  INTEGER                        :: iomega
  DOUBLE PRECISION               :: pi_invBeta
  DOUBLE PRECISION               :: mbeta_two
  DOUBLE PRECISION               :: beta 
  DOUBLE PRECISION               :: beta_tc
  DOUBLE PRECISION               :: tcbeta_tc
  DOUBLE PRECISION               :: inv_dt 
  DOUBLE PRECISION               :: tC
  DOUBLE PRECISION               :: tCdag
  DOUBLE PRECISION               :: time
  DOUBLE PRECISION               :: signe
  DOUBLE PRECISION               :: argument
  !DOUBLE PRECISION               :: taupi_invbeta
  !COMPLEX(KIND=8)                   :: cargument
  !COMPLEX(2*8)                   :: base_exp
  !COMPLEX(2*8)                   :: increm_exp

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_measHybrid : green operator not set         ")

  tail = ListCdagC_1%tail
  IF ( tail .NE. Mthis%tail ) &
    CALL ERROR("GreenHyb_measHybrid : ListCdagC & M unconsistent     ")

  IF ( updated .EQV. .TRUE. ) THEN ! NEW change in the configuration
    ! FIXME SHOULD be much more faster
    
      old_size = this%this%tail
    SELECT CASE(this%iTech)
    CASE (GREENHYB_TAU)
      argument = DBLE(this%factor)
      DO iC = 1, old_size
        this%oper(this%this%listINT(iC)) = this%oper(this%this%listINT(iC)) + this%this%listDBLE(iC) * argument
      END DO
      this%measurements = this%measurements + this%factor
  
      CALL MapHyb_setSize(this%this,tail*tail)
      this%factor = 1
      idx_old = 0
      beta   =  this%beta
      mbeta_two = -(beta*0.5d0)
      inv_dt =  this%inv_dt
      ! WARNING time is not the time but just a temporary variable.
      ! Index Time has been calculated previously and is in mat_tau
      DO iC  = 1, tail
        tC   = ListCdagC_1%list(iC,C_)
        beta_tc = beta - tC
        tcbeta_tc = tC * beta_tc
        DO iCdag = 1, tail
          tCdag  = ListCdagC_1%list(iCdag,Cdag_)
          time = tcbeta_tc - tCdag*beta_tc
  
          !signe = SIGN(1.d0,time)
          !time = time + (signe-1.d0)*mbeta_two
    
          !signe = signe * SIGN(1.d0,beta-tC)

          !signe = SIGN(1.d0,time) * SIGN(1.d0,beta-tC)
          signe = SIGN(1.d0,time)
  
          argument = signe*Mthis%mat(iCdag,iC)
  
          !index = INT( ( time * inv_dt ) + 1.5d0 )
          !IF (index .NE. Mthis%mat_tau(iCdag,iC)) THEN
          !  WRITE(*,*) index, Mthis%mat_tau(iCdag,iC)
          !!  CALL ERROR("Plantage")
          !END IF
  
          idx_old = idx_old + 1
          this%this%listDBLE(idx_old) = argument
          !this%this%listINT(idx_old)  = index
          this%this%listINT(idx_old)  = Mthis%mat_tau(iCdag,iC)
        END DO
      END DO
    CASE (GREENHYB_OMEGA)
      omegaSamples = this%Wmax
      argument = DBLE(this%factor)
      DO iomega = 1, omegaSamples
        this%oper_w(iomega) = this%oper_w(iomega) + this%oper_w_old(iomega) * argument
      END DO
      this%measurements = this%measurements + this%factor

      this%factor = 1
      beta   =  this%beta
      mbeta_two = -(beta*0.5d0)
      pi_invBeta = ACOS(-1.d0)/beta
      DO iC  = 1, tail
        tC   = ListCdagC_1%list(iC,C_)
        DO iCdag = 1, tail
          tCdag  = ListCdagC_1%list(iCdag,Cdag_)
          time = tC - tCdag

          signe = SIGN(1.d0,time)
          time = time + (signe-1.d0)*mbeta_two
          signe = signe * SIGN(1.d0,beta-tC)
          argument = signe*Mthis%mat(iCdag,iC)

          DO iomega = 1, omegaSamples
            !this%oper_w_old(iomega) = Mthis%mat_tau(iCdag,iC)*CMPLX(0.d0,argument)
            this%oper_w_old(iomega) = EXP(CMPLX(0.d0,this%omega(iomega)*time,8))*CMPLX(0.d0,argument,8)
          END DO
        END DO
      END DO
    END SELECT
  ELSE
    this%factor = this%factor + 1
  END IF
END SUBROUTINE GreenHyb_measHybrid
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_getHybrid
!! NAME
!!  GreenHyb_getHybrid
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
!!  this=Green
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

SUBROUTINE GreenHyb_getHybrid(this)

!Arguments ------------------------------------
  TYPE(GreenHyb), INTENT(INOUT) :: this

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_getHybrid : green operator not set          ")

  SELECT CASE(this%iTech)
  CASE (GREENHYB_TAU)
    this%oper = -(this%oper * this%inv_beta) / (DBLE(this%measurements) * this%delta_t)
    this%setT = .TRUE.
  CASE (GREENHYB_OMEGA)
    this%oper_w = -(this%oper_w * this%inv_beta) / (DBLE(this%measurements) * this%delta_t)
    this%setW = .TRUE.
    CALL GreenHyb_backFourier(this)
  END SELECT

END SUBROUTINE GreenHyb_getHybrid
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_setN
!! NAME
!!  GreenHyb_setN
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
!!  this=Green
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

SUBROUTINE GreenHyb_setN(this,N)

!Arguments ------------------------------------
  TYPE(GreenHyb)     , INTENT(INOUT) :: this
  DOUBLE PRECISION, INTENT(IN   ) :: N

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_setN: green this%operator not set                ")
  this%oper(1) = N - 1.d0
  this%oper(this%samples) = - N
END SUBROUTINE GreenHyb_setN
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_setMuD1
!! NAME
!!  GreenHyb_setMuD1
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
!!  this=Green
!!  mu=energy level (irrespectige with fermi level)
!!  d1=firt moment of hybridization function
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

SUBROUTINE GreenHyb_setMuD1(this,mu,d1)

!Arguments ------------------------------------
  TYPE(GreenHyb)  , INTENT(INOUT) :: this
  DOUBLE PRECISION, INTENT(IN   ) :: mu
  DOUBLE PRECISION, INTENT(IN   ) :: d1

  this%Mk(3) = -d1-(mu*mu)
  this%Mk(2) = -mu
  this%setMk = this%setMk + 1
END SUBROUTINE GreenHyb_setMuD1
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_setMoments
!! NAME
!!  GreenHyb_setMoments
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
!!  this=Green
!!  u1=interaction energi like
!!  u2=idem order2
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

SUBROUTINE GreenHyb_setMoments(this,u1,u2)

!Arguments ------------------------------------
  TYPE(GreenHyb)  , INTENT(INOUT) :: this
  DOUBLE PRECISION, INTENT(IN   ) :: u1
  DOUBLE PRECISION, INTENT(IN   ) :: u2
  
  this%Mk(1) = -1.d0
  this%Mk(3) = this%Mk(3) - 2.d0*(this%Mk(2)*u1)
  this%Mk(2) = this%Mk(2) + u1
  this%Mk(3) = this%Mk(3) - u2

  this%setMk = this%setMk + 1

END SUBROUTINE GreenHyb_setMoments
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_backFourier
!! NAME
!!  GreenHyb_backFourier
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
!!  this=Green
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

SUBROUTINE GreenHyb_backFourier(this,dvgc)


#ifdef HAVE_MPI1
include 'mpif.h'
#endif
!Arguments ------------------------------------
  TYPE(GreenHyb)            , INTENT(INOUT) :: this
  DOUBLE PRECISION, OPTIONAL, INTENT(IN   ) :: dvgc
!Local variables ------------------------------
  INTEGER :: itau
  INTEGER :: iomega
  INTEGER :: omegaSamples
  INTEGER :: tauSamples
  INTEGER :: tauBegin
  INTEGER :: tauEnd
  INTEGER :: delta
  INTEGER :: residu
  INTEGER, ALLOCATABLE, DIMENSION(:) :: counts
  INTEGER, ALLOCATABLE, DIMENSION(:) :: displs
  DOUBLE PRECISION :: A ! Correction factor
  DOUBLE PRECISION :: inv_beta
  DOUBLE PRECISION :: pi_invBeta
  DOUBLE PRECISION :: two_invBeta
  DOUBLE PRECISION :: minusDt
  DOUBLE PRECISION :: minusOmegaTau
  DOUBLE PRECISION :: omega
  DOUBLE PRECISION :: minusTau
  DOUBLE PRECISION :: sumTerm
  DOUBLE PRECISION :: pi
  DOUBLE PRECISION :: twoPi
  DOUBLE PRECISION :: correction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Domega
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: A_omega

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_backFourier : Uninitialized GreenHyb structure")
  IF ( this%setW .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_backFourier : no G(iw)")
  
  inv_beta     = this%inv_beta
  two_invBeta  = 2.d0 * inv_beta
  minusDt      = - this%delta_t
  omegaSamples = this%Wmax
  tauSamples   = this%samples-1

  this%oper = 0.d0

  pi         = ACOS(-1.d0)
  twoPi        = 2.d0 * pi
  pi_invBeta = pi * inv_beta

  IF ( PRESENT(dvgc) ) THEN
    A = dvgc
  ELSE
    A = AIMAG(this%oper_w(omegaSamples)) &  ! A = \lim_\infty G(w)
             *(2.d0*DBLE(omegaSamples)-1.d0) * pi_invBeta
  END IF

  correction = A*0.5d0

  MALLOC(Domega,(1:omegaSamples))
  MALLOC(A_omega,(1:omegaSamples))
  Domega = (/ ((2.d0 * DBLE(iomega) - 1.d0)*pi_invbeta, iomega=1, omegaSamples) /)
  A_omega = A / Domega
  IF (this%have_MPI .EQV. .TRUE.) THEN
    delta = tauSamples / this%size
    residu = tauSamples - this%size*delta
    IF ( this%rank .LT. this%size - residu ) THEN
      tauBegin = 1 + this%rank*delta
      tauEnd   = (this%rank + 1)*delta
    ELSE
!      tauBegin = (this%size-residu)*delta + 1 + (this%rank-this%size+residu)*(delta+1)
      tauBegin = 1 + this%rank*(delta + 1) -this%size + residu
      tauEnd = tauBegin + delta
    END IF
    MALLOC(counts,(1:this%size))
    MALLOC(displs,(1:this%size))
    counts = (/ (delta, iTau=1, this%size-residu), &
                (delta+1, iTau=this%size-residu+1, this%size) /)
    displs(1)=0
    DO iTau = 2, this%size
      displs(iTau) = displs(iTau-1) + counts (iTau-1)
    END DO
  ELSE
    tauBegin = 1
    tauEnd   = tauSamples
  END IF
  DO iTau = tauBegin, tauEnd
    minusTau = DBLE(itau -1) * minusDt
    DO iomega = 1, omegaSamples
      omega         = Domega(iomega)
      minusOmegaTau = MOD(omega*minusTau, TwoPi)
      sumTerm       = REAL(( this%oper_w(iomega) - CMPLX(0.d0, A_omega(iomega),8) ) &
                      * EXP( CMPLX(0.d0, minusOmegaTau, 8)))
      this%oper(itau)    = this%oper(itau) + sumTerm
    END DO
    this%oper(itau) = correction + two_invBeta*this%oper(itau)
  END DO
  IF ( this%have_MPI .EQV. .TRUE. ) THEN
! rassembler les resultats
#ifdef HAVE_MPI
    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                      this%oper, counts, displs, &
                      MPI_DOUBLE_PRECISION, this%MY_COMM, residu)
#endif
    FREE(counts)
    FREE(displs)
  END IF
  this%oper(tauSamples+1) = A - this%oper(1) !G(0+)-G(0-)=G(0+)+G(beta-)=A
  this%setT = .TRUE.
  FREE(Domega)
  FREE(A_omega)

END SUBROUTINE GreenHyb_backFourier
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_forFourier
!! NAME
!!  GreenHyb_forFourier
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
!!  this=Green
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

SUBROUTINE GreenHyb_forFourier(this, Gomega, omega, Wmax)
!Arguments ------------------------------------

#ifdef HAVE_MPI1
include 'mpif.h'
#endif
  TYPE(GreenHyb)             , INTENT(INOUT) :: this
  COMPLEX(KIND=8), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: Gomega  ! INOUT for MPI
  COMPLEX(KIND=8), DIMENSION(:), OPTIONAL, INTENT(IN   ) :: omega  
  INTEGER                 , OPTIONAL, INTENT(IN   ) :: Wmax   
  INTEGER :: i
  INTEGER :: j
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

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_forFourier : Uninitialized GreenHyb structure")
  IF ( this%setT .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_forFourier : no G(tau)")
  IF ( this%setMk .NE. 2 ) &
    CALL WARNALL("GreenHyb_forFourier : green does not have moments    ")

  L  = this%samples

  xpi=acos(-1.d0)                !!! XPI=PI
  beta = this%beta
  Nom  = this%Wmax
  IF ( PRESENT(Gomega) ) THEN
    Nom = SIZE(Gomega)    
    !IF ( this%rank .EQ. 0 ) &
      !write(6,*) "size Gomega", Nom
  END IF
  IF ( PRESENT(omega) ) THEN
    IF ( PRESENT(Gomega) .AND. SIZE(omega) .NE. Nom ) THEN
      CALL ERROR("GreenHyb_forFourier : sizes mismatch              ")               
    !ELSE 
      !Nom = SIZE(omega)
    END IF 
  END IF
  IF ( .NOT. PRESENT(Gomega) .AND. .NOT. PRESENT(omega) ) THEN
    IF ( PRESENT(Wmax) ) THEN
      Nom=Wmax
    ELSE
      CALL ERROR("GreenHyb_forFourier : Missing argument Wmax")
    END IF
  END IF

  IF ( ALLOCATED(this%oper_w) ) THEN
    IF ( SIZE(this%oper_w) .NE. Nom ) THEN
      FREE(this%oper_w)
      MALLOC(this%oper_w,(1:Nom))
    END IF
  ELSE
    MALLOC(this%oper_w,(1:Nom))
  END IF

  !write(6,*) "PRESENT(GOMEGA)", PRESENT(GOMEGA)
  !write(6,*) "PRESENT(OMEGA)", PRESENT(OMEGA)
  !call flush(6)

  delta=this%delta_t
  inv_delta = this%inv_dt
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
  
! construct the B vector from A.Xm=B
  MALLOC(XM,(L))
  XM(1) = 4.d0*this%Mk(3)
  XM(L) = (6.d0 * inv_delta) * ( this%Mk(2) - ( &
          (this%oper(2)-this%oper(1)) + &
          (this%oper(L)-this%oper(L-1)) ) * inv_delta )
  DO i = 2, L-1
    XM(i) = (6.d0 * inv_delta2) * ( (this%oper(i+1) &
                          - 2.d0 * this%oper(i)) &
                          +        this%oper(i-1) )
  END DO

! Find second derivatives XM: Solve the system
! SOLVING Lq= XM 
!  q = XM 
  do j=1,L-1
      XM(j+1)=XM(j+1)-(diagL(j)*XM(j))
      XM(L)  =XM(L)  -(lastR(j)*XM(j))
  end do
  FREE(diagL)
  FREE(lastR)

! SOLVING U.XM=q 
!  XM = q
  do j=L-1,2,-1
   XM(j+1)  = XM(j+1) / diag(j+1)
   XM(j)= (XM(j)-(XM(L)*lastC(j)))-XM(j+1)
  end do
  XM(2)  = XM(2) / diag(2)
  XM(1) = (XM(1)-XM(L)*lastC(1)) / diag(1)

  FREE(diag)
  FREE(lastC)

  Lspline = L-1
  MALLOC(X2,(1:Lspline+1)) ! We impose L = Nom
  !Construct L2 second derivative from known derivatives XM
  deltabis = beta / DBLE(Lspline)
  DO i = 1, Lspline
    tau = deltabis * DBLE(i-1)
    j = ((L-1)*(i-1))/Lspline + 1!INT(tau * inv_delta) + 1
    X2(i) = inv_delta * ( XM(j)*(DBLE(j)*delta - tau ) + XM(j+1)*(tau - DBLE(j-1)*delta) )
  END DO
  X2(Lspline+1) = XM(L)
  FREE(XM)

  IF ( this%have_MPI .EQV. .TRUE. ) THEN  
    deltaw = Nom / this%size
    residu = Nom - this%size*deltaw
    IF ( this%rank .LT. this%size - residu ) THEN
      omegaBegin = 1 + this%rank*deltaw
      omegaEnd   = (this%rank + 1)*deltaw
    ELSE
  !    tauBegin = (this%size-residu)*deltaw + 1 + (this%rank-this%size+residu)*(deltaw+1)
      omegaBegin = 1 + this%rank*(deltaw + 1) -this%size + residu
      omegaEnd = omegaBegin + deltaw
    END IF
    MALLOC(counts,(1:this%size))
    MALLOC(displs,(1:this%size))
    counts = (/ (deltaw, i=1, this%size-residu), &
                (deltaw+1, i=this%size-residu+1, this%size) /)
    displs(1)=0
    DO i = 2, this%size
      displs(i) = displs(i-1) + counts (i-1)
    END DO
  ELSE
    omegaBegin = 1
    omegaEnd   = Nom 
  END IF

  this%Mk(1) = -1.d0
  MALLOC(omegatmp,(omegaBegin:omegaEnd))
  IF ( PRESENT(omega) ) THEN
    omegatmp(omegaBegin:omegaEnd) = (/ (AIMAG(omega(i)),i=omegaBegin,omegaEnd) /)
  ELSE
    omegatmp(omegaBegin:omegaEnd) = (/ ((((2.d0*DBLE(i)-1.d0)*xpi)/Beta), i=omegaBegin,omegaEnd) /)
  END IF
  MALLOC(Gwtmp,(1:Nom))
  DO i = omegaBegin, omegaEnd
    iw = omegatmp(i)
    omdeltabis = iw*deltabis
    Gwtmp(i)=CMPLX(0.d0,0.d0,8)
    DO j=2, Lspline ! We impose  L+1 = Nom
      iwtau = CMPLX(0.d0,omdeltabis*DBLE(j-1),8)
      Gwtmp(i) = Gwtmp(i) + EXP(iwtau) * CMPLX((X2(j+1) + X2(j-1))-2.d0*X2(j),0.d0,8)
    END DO
    Gwtmp(i) = Gwtmp(i)/CMPLX(((iw*iw)*(iw*iw)*deltabis),0.d0,8) &

              + CMPLX( &
                ( ((X2(2)-X2(1))+(X2(Lspline+1)-X2(Lspline)))/((iw*iw)*deltabis) -this%Mk(2) ) &
                /(iw*iw) &
               , &
                (this%Mk(1)-this%Mk(3)/(iw*iw))/iw &
               , 8) 
              !+ CMPLX( (X2(2)-X2(1))+(X2(Lspline+1)-X2(Lspline)), 0.d0, 8 ) ) &
              !   / (((iw*iw)*(iw*iw))*CMPLX(deltabis,0.d0,8)) &
              !- CMPLX(this%Mk(1),0.d0,8)/iw  &
              !+ CMPLX(this%Mk(2),0.d0,8)/(iw*iw) &
              !- CMPLX(this%Mk(3),0.d0,8)/((iw*iw)*iw) 
    !IF ( this%rank .EQ. 0 )  write(12819,*) iw,gwtmp(i)
  END DO
  FREE(omegatmp)
  !call flush(12819)
  FREE(X2)
  IF ( this%have_MPI .EQV. .TRUE. ) THEN
#ifdef HAVE_MPI
    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                      Gwtmp  , counts, displs, &
                      MPI_DOUBLE_COMPLEX, this%MY_COMM, residu)
#endif
    FREE(counts)
    FREE(displs)
  END IF
  IF ( PRESENT(Gomega) ) THEN
    Gomega = Gwtmp
  END IF
  this%oper_w = Gwtmp
  this%setW = .TRUE.
  FREE(Gwtmp)
END SUBROUTINE GreenHyb_forFourier
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_print
!! NAME
!!  GreenHyb_print
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
!!  this=Green
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

SUBROUTINE GreenHyb_print(this, ostream)

!Arguments ------------------------------------
  TYPE(GreenHyb), INTENT(IN) :: this
  INTEGER, OPTIONAL , INTENT(IN) :: ostream
!Local variables ------------------------------
  INTEGER                        :: ostream_val
  INTEGER                        :: i
  INTEGER                        :: samples
  

  IF ( this%set .EQV. .FALSE. ) &
    CALL ERROR("GreenHyb_print : green this%operator not set              ")

  IF ( PRESENT(ostream) ) THEN
    ostream_val = ostream
  ELSE
    ostream_val = 66
    OPEN(UNIT=ostream_val,FILE="Green.dat")
  END IF

  samples =  this%samples 

  DO i = 1, samples
    WRITE(ostream_val,*) DBLE(i-1)*this%delta_t, this%oper(i)
  END DO

  IF ( .NOT. PRESENT(ostream) ) &
    CLOSE(ostream_val)
END SUBROUTINE GreenHyb_print
!!***

!!****f* ABINIT/m_GreenHyb/GreenHyb_destroy
!! NAME
!!  GreenHyb_destroy
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
!!  this=Green
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

SUBROUTINE GreenHyb_destroy(this)

!Arguments ------------------------------------
  TYPE(GreenHyb), INTENT(INOUT) :: this

  this%set          = .FALSE.
  this%setT         = .FALSE.
  this%setW         = .FALSE.
  this%samples      = 0
  this%measurements = 0
  this%beta         = 0.d0
  this%inv_beta     = 0.d0
  this%inv_dt       = 0.d0
  this%delta_t      = 0.d0
  !CALL VectorInt_destroy(this%index_old)
  !CALL Vector_destroy(this%oper_old)
  CALL MapHyb_destroy(this%this)
  FREEIF(this%oper)
  FREEIF(this%oper_w)
  FREEIF(this%oper_w_old)
  FREEIF(this%omega)
END SUBROUTINE GreenHyb_destroy
!!***

END MODULE m_GreenHyb
!!***

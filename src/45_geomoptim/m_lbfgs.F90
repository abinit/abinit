!!****m* ABINIT/m_lbfgs
!! NAME
!!  m_lbfgs
!!
!! FUNCTION
!!  This module provides several routines for the application of a
!!  Limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS) minimization algorithm.
!!  The working routines were based on the original implementation of J. Nocera available on netlib.org
!!  They have been reshaped and translated into modern fortran here.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2020 ABINIT group (FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lbfgs

 use defs_basis
 use m_abicore

 implicit none

type,public :: lbfgs_internal
 integer              :: lbfgs_status
 integer              :: ndim
 integer              :: history_record
 integer              :: iter
 real(dp)             :: gtol
 real(dp),allocatable :: diag(:)
 real(dp),allocatable :: work(:)
 real(dp)             :: line_stp
 real(dp)             :: line_stpmin
 real(dp)             :: line_stpmax
 integer              :: line_info
 integer              :: line_infoc
 integer              :: line_nfev
 real(dp)             :: line_dginit
 real(dp)             :: line_finit
 real(dp)             :: line_stx
 real(dp)             :: line_fx
 real(dp)             :: line_dgx
 real(dp)             :: line_sty
 real(dp)             :: line_fy
 real(dp)             :: line_dgy
 real(dp)             :: line_stmin
 real(dp)             :: line_stmax
 logical              :: line_bracket
 logical              :: line_stage1
end type lbfgs_internal

type(lbfgs_internal),save,public :: lbfgs_plan

!!***

contains

!----------------------------------------------------------------------

!!****f* m_lbfgs/lbfgs_init
!! NAME
!! lbfgs_init
!!
!! FUNCTION
!!   Initialize the internal object lbfgs_internal for LBFGS minimization
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_pred_bfgs
!!
!! CHILDREN
!!
!! SOURCE

subroutine lbfgs_init(ndim,history_record,diag_guess)

integer,intent(in)   :: ndim
integer,intent(in)   :: history_record
real(dp),intent(in)  :: diag_guess(ndim)

integer :: nwork

 lbfgs_plan%lbfgs_status = 0
 lbfgs_plan%iter   = 0
 lbfgs_plan%ndim = ndim
 lbfgs_plan%history_record = history_record
 ABI_MALLOC(lbfgs_plan%diag,(ndim))

 nwork = ndim * ( 2 * history_record + 1 ) + 2 * history_record
 ABI_MALLOC(lbfgs_plan%work,(nwork))

 lbfgs_plan%gtol = 0.9
 lbfgs_plan%line_stpmin = 1.0e-20
 lbfgs_plan%line_stpmax = 1.0e+20
 lbfgs_plan%line_stp    = 1.0

 lbfgs_plan%diag(:) = diag_guess(:)

end subroutine lbfgs_init
!!***


!----------------------------------------------------------------------

!!****f* m_lbfgs/lbfgs_destroy
!! NAME
!! lbfgs_destroy
!!
!! FUNCTION
!!   Free the memory of the internal object lbfgs_internal for LBFGS minimization
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_pred_bfgs
!!
!! CHILDREN
!!
!! SOURCE

subroutine lbfgs_destroy()

 if(allocated (lbfgs_plan%work)) then
   ABI_FREE(lbfgs_plan%work)
 end if
 if(allocated (lbfgs_plan%diag)) then
   ABI_FREE(lbfgs_plan%diag)
 end if

end subroutine lbfgs_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_lbfgs/lbfgs_execute
!! NAME
!! lbfgs_execute
!!
!! FUNCTION
!!   Perform one-step of LBFGS minimization
!!   all the internal information are stored in the lbfgs_internal object
!!
!! INPUTS
!!   x: input and output position vector (atomic reduced coordinates + cell parameters)
!!   f: total energy
!!   gradf: gradient of the total energy (=negative forces)
!!
!! OUTPUT
!!
!! PARENTS
!!      pred_lbfgs
!!
!! CHILDREN
!!
!! SOURCE

function lbfgs_execute(x,f,gradf)

real(dp),intent(inout) :: x(lbfgs_plan%ndim)
real(dp),intent(in)    :: f
real(dp),intent(in)    :: gradf(lbfgs_plan%ndim)
integer                :: lbfgs_execute

 call lbfgs(lbfgs_plan%ndim, lbfgs_plan%history_record, x, f, gradf, lbfgs_plan%diag, lbfgs_plan%work, lbfgs_plan%lbfgs_status, &
       lbfgs_plan%gtol, lbfgs_plan%line_stpmin, lbfgs_plan%line_stpmax, lbfgs_plan%line_stp, lbfgs_plan%iter,                   &
       lbfgs_plan%line_info, lbfgs_plan%line_nfev,                                                   &
       lbfgs_plan%line_dginit, lbfgs_plan%line_finit,                                                &
       lbfgs_plan%line_stx,  lbfgs_plan%line_fx,  lbfgs_plan%line_dgx,                               &
       lbfgs_plan%line_sty,  lbfgs_plan%line_fy,  lbfgs_plan%line_dgy,                               &
       lbfgs_plan%line_stmin,  lbfgs_plan%line_stmax,                                                &
       lbfgs_plan%line_bracket, lbfgs_plan%line_stage1, lbfgs_plan%line_infoc)


!lbfgs_execute = lbfgs_plan%lbfgs_status
 lbfgs_execute = lbfgs_plan%line_info

end function lbfgs_execute
!!***


!----------------------------------------------------------------------

!!****f* m_lbfgs/lbfgs
!! NAME
!! lbfgs
!!
!! FUNCTION
!!   Perform the LBFGS step
!!   Fortran90 rewritting of the original subroutine by J. Nocera
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_lbfgs
!!
!! CHILDREN
!!
!! SOURCE

subroutine lbfgs(N,M,X,F,G,DIAG,W,IFLAG,      &
                 GTOL,STPMIN,STPMAX,STP,ITER, &
                 INFO, NFEV,                  &
                 LINE_DGINIT,LINE_FINIT,      &
                 LINE_STX,LINE_FX,LINE_DGX,   &
                 LINE_STY,LINE_FY,LINE_DGY,   &
                 LINE_STMIN,LINE_STMAX,       &
                 LINE_BRACKT,LINE_STAGE1,LINE_INFOC)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: LINE_INFOC
 integer,intent(inout)  :: ITER,IFLAG,INFO,NFEV
 integer,intent(in)     :: N,M
 real(dp),intent(inout) :: GTOL
 real(dp),intent(in)    :: STPMIN,STPMAX
 real(dp),intent(inout) :: STP
 real(dp),intent(in)    :: F
 real(dp),intent(inout) :: LINE_DGINIT,LINE_FINIT
 real(dp),intent(inout) :: LINE_STX,LINE_FX,LINE_DGX
 real(dp),intent(inout) :: LINE_STY,LINE_FY,LINE_DGY
 real(dp),intent(inout) :: LINE_STMIN,LINE_STMAX
 logical,intent(inout)  :: LINE_BRACKT,LINE_STAGE1
!arrays
 real(dp),intent(inout) :: X(N),DIAG(N),W(N*(2*M+1)+2*M)
 real(dp),intent(in)    :: G(N)

!Local variables-------------------------------
!scalars
 real(dp) :: FTOL,YS,YY,SQ,YR,BETA
 integer :: POINT,ISPT,IYPT,MAXFEV, &
         BOUND,NPT,CP,I,INMC,IYCN,ISCN
!***************************************************************************


!
! Initialize
!-----------

! Parameters for line search routine
 FTOL = 1.0D-4
 MAXFEV = 20

 ISPT = N + 2 * M
 IYPT = ISPT + N * M
 POINT = MAX( 0 , MOD(ITER-1,M) )
 NPT = POINT * N
 ITER  = ITER + 1
 BOUND = MIN( ITER-1 , M)


 !
 ! Entering the subroutine with a new position and gradient
 ! or entering for the first time ever
 if( IFLAG /= 1 ) then
   W(ISPT+1:ISPT+N) = -G(1:N) * DIAG(1:N)

 else

   call MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,MAXFEV,INFO,NFEV, &
               DIAG,GTOL,STPMIN,STPMAX,LINE_DGINIT,LINE_FINIT, &
               LINE_STX,LINE_FX,LINE_DGX, &
               LINE_STY,LINE_FY,LINE_DGY, &
               LINE_STMIN,LINE_STMAX, &
               LINE_BRACKT,LINE_STAGE1,LINE_INFOC)
   !
   ! Compute the new step and gradient change
   !
   NPT = POINT * N
   W(ISPT+NPT+1:ISPT+NPT+N) = STP * W(ISPT+NPT+1:ISPT+NPT+N)
   W(IYPT+NPT+1:IYPT+NPT+N) = G(1:N) - W(1:N)

   YS = DOT_PRODUCT( W(IYPT+NPT+1:IYPT+NPT+N) , W(ISPT+NPT+1:ISPT+NPT+N) )
   YY = DOT_PRODUCT( W(IYPT+NPT+1:IYPT+NPT+N) , W(IYPT+NPT+1:IYPT+NPT+N) )
   DIAG(1:N)= YS / YY

!
!  COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!  "Updating quasi-Newton matrices with limited storage",
!  Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!  ---------------------------------------------------------
!
   POINT = MODULO(ITER - 1,M)
   CP = POINT
   if (POINT == 0) CP = M

   W(N+CP) = one / YS
   W(1:N)  = -G(1:N)

   CP = POINT
   do I= 1,BOUND
     CP = CP - 1
     if (CP ==  -1) CP = M - 1
     SQ = DOT_PRODUCT(W(ISPT+CP*N+1:ISPT+CP*N+N),W(1:N))
     INMC = N + M + CP + 1
     IYCN = IYPT + CP * N
     W(INMC)= W(N+CP+1) * SQ
     W(1:N) = W(1:N) - W(INMC) * W(IYCN+1:IYCN+N)
   enddo

   W(1:N) = DIAG(1:N) * W(1:N)

   do I=1,BOUND
     YR = DOT_PRODUCT(W(IYPT+CP*N+1:IYPT+CP*N+N),W(1:N))
     BETA = W(N+CP+1) * YR
     INMC = N + M + CP + 1
     BETA = W(INMC) - BETA
     ISCN = ISPT + CP * N
     W(1:N) = W(1:N) + BETA * W(ISCN+1:ISCN+N)
     CP = CP + 1
     if (CP == M) CP = 0
   enddo

!
!  STORE THE NEW SEARCH DIRECTION
   W(ISPT+POINT*N+1:ISPT+POINT*N+N) = W(1:N)

 endif

!
! Obtain the one-dimensional minimizer of the function
! by using the line search routine mcsrch
!----------------------------------------------------
 NFEV = 0
 STP = one
 W(1:N) = G(1:N)

 INFO  = 0

 call MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,MAXFEV,INFO,NFEV, &
             DIAG,GTOL,STPMIN,STPMAX,LINE_DGINIT,LINE_FINIT, &
             LINE_STX,LINE_FX,LINE_DGX, &
             LINE_STY,LINE_FY,LINE_DGY, &
             LINE_STMIN,LINE_STMAX, &
             LINE_BRACKT,LINE_STAGE1,LINE_INFOC)

 if (INFO  ==  -1) then
   IFLAG = 1
   return
 else
   IFLAG = -1
   return
 endif

end subroutine lbfgs
!!***

!----------------------------------------------------------------------

!!****f* m_lbfgs/mcsrch
!! NAME
!! mcsrch
!!
!! FUNCTION
!!   Perform the line minimization step
!!   Fortran90 rewritting of the original subroutine by J. Nocera
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_lbfgs
!!
!! CHILDREN
!!
!! SOURCE

subroutine mcsrch(N,X,F,G,S,STP,FTOL,MAXFEV,INFO,NFEV,WA, &
                  GTOL,STPMIN,STPMAX,DGINIT,FINIT, &
                  STX,FX,DGX,STY,FY,DGY,STMIN,STMAX, &
                  BRACKT,STAGE1,INFOC)

!Arguments ------------------------------------
!scalars
 integer,intent(in)     :: N,MAXFEV
 integer,intent(inout)  :: INFO,NFEV
 integer,intent(inout)  :: INFOC
 real(dp),intent(in)    :: GTOL,STPMIN,STPMAX
 real(dp),intent(in)    :: F,FTOL
 real(dp),intent(inout) :: STP,DGINIT,FINIT
 real(dp),intent(inout) :: STX,FX,DGX
 real(dp),intent(inout) :: STY,FY,DGY
 real(dp),intent(inout) :: STMIN,STMAX
 logical,intent(inout) :: BRACKT,STAGE1
!arrays
 real(dp),intent(in)     :: G(N)
 real(dp),intent(inout)  :: X(N),S(N),WA(N)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: XTOL=1.0e-17_dp
 real(dp),parameter :: P5     = 0.50_dp
 real(dp),parameter :: P66    = 0.66_dp
 real(dp),parameter :: XTRAPF = 4.00_dp
 real(dp) :: DG,DGM,DGTEST,DGXM,DGYM, &
        FTEST1,FM,FXM,FYM,WIDTH,WIDTH1
!***************************************************************************

 DGTEST = FTOL * DGINIT
 WIDTH = STPMAX - STPMIN
 WIDTH1 = WIDTH / P5

 ! Is it a first entry (info == 0)
 ! or a second entry (info == -1)?
 if( INFO == -1 ) then

   ! Reset INFO
   INFO = 0

   NFEV = NFEV + 1
   DG = SUM( G(:) * S(:) )
   FTEST1 = FINIT + STP * DGTEST
!
!  TEST FOR CONVERGENCE.
!
   if ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) &
      .OR. INFOC  ==  0) INFO = 6
   if (STP  ==  STPMAX .AND. &
       F <= FTEST1 .AND. DG <= DGTEST) INFO = 5
   if (STP  ==  STPMIN .AND.  &
       (F > FTEST1 .OR. DG >= DGTEST)) INFO = 4
   if (NFEV >= MAXFEV) INFO = 3
   if (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX) INFO = 2
   if (F <= FTEST1 .AND. ABS(DG) <= GTOL*(-DGINIT)) INFO = 1
!
!  CHECK FOR TERMINATION.
!
   if (INFO /= 0) return
!
!  IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
!  FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
!
   if (STAGE1 .AND. F <= FTEST1 .AND. &
       DG >= MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
!
!  A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
!  WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
!  FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
!  DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
!  OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
!
   if (STAGE1 .AND. F <= FX .AND. F > FTEST1) then
!
!     DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
!
      FM = F - STP * DGTEST
      FXM = FX - STX * DGTEST
      FYM = FY - STY * DGTEST
      DGM = DG - DGTEST
      DGXM = DGX - DGTEST
      DGYM = DGY - DGTEST
!
!     CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!     AND TO COMPUTE THE NEW STEP.
!
      call mcstep(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,BRACKT,STMIN,STMAX,INFOC)
!
!     RESET THE FUNCTION AND GRADIENT VALUES FOR F.
!
      FX = FXM + STX * DGTEST
      FY = FYM + STY * DGTEST
      DGX = DGXM + DGTEST
      DGY = DGYM + DGTEST
   else
!
!     CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!     AND TO COMPUTE THE NEW STEP.
!
      call mcstep(STX,FX,DGX,STY,FY,DGY,STP,F,DG,BRACKT,STMIN,STMAX,INFOC)
   endif
!
!  FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
!  INTERVAL OF UNCERTAINTY.
!
   if (BRACKT) then
      if (ABS(STY-STX) >= P66 * WIDTH1) STP = STX + P5 * (STY - STX)
      WIDTH1 = WIDTH
      WIDTH = ABS(STY-STX)
   endif

 else

   INFOC = 1
!
!  CHECK THE INPUT PARAMETERS FOR ERRORS.
!
   if ( STP <= zero .OR. FTOL < zero .OR.  &
       GTOL < zero .OR. XTOL < zero .OR. STPMIN < zero &
       .OR. STPMAX < STPMIN ) return
!
!  COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
!  AND CHECK THAT S IS A DESCENT DIRECTION.
!
   DGINIT = DOT_PRODUCT( G , S )

   if (DGINIT > zero) then
     return
   endif
!
!  INITIALIZE LOCAL VARIABLES.
!

   BRACKT = .FALSE.
   STAGE1 = .TRUE.
   NFEV = 0
   FINIT = F
   DGTEST = FTOL * DGINIT
   WA(:) = X(:)


!
!  THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
!  FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
!  THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
!  FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
!  THE INTERVAL OF UNCERTAINTY.
!  THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
!  FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
!
   STX = zero
   FX = FINIT
   DGX = DGINIT
   STY = zero
   FY = FINIT
   DGY = DGINIT
 endif

!
!SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
!TO THE PRESENT INTERVAL OF UNCERTAINTY.
!
 if (BRACKT) then
    STMIN = MIN(STX,STY)
    STMAX = MAX(STX,STY)
 else
    STMIN = STX
    STMAX = STP + XTRAPF*(STP - STX)
 endif
!
!FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
!
 STP = MAX(STPMIN,STP)
 STP = MIN(STP,STPMAX)
!
!IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
!STP BE THE LOWEST POINT OBTAINED SO FAR.
!
 if ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) &
    .OR. NFEV >= MAXFEV-1 .OR. INFOC  ==  0 &
    .OR. (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX)) STP = STX

!
!Evaluate the function and gradient at STP
!and compute the directional derivative.
!We return to main program to obtain F and G.
!
 X(:) = WA(:) + STP * S(:)

 INFO = -1

end subroutine mcsrch
!!***

!----------------------------------------------------------------------

!!****f* m_lbfgs/mcstep
!! NAME
!! mcstep
!!
!! FUNCTION
!!   Perform the step choice in line minimization
!!   Fortran90 rewritting of the original subroutine by J. Nocera
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_lbfgs
!!
!! CHILDREN
!!
!! SOURCE

subroutine mcstep(STX,FX,DX,STY,FY,DY,STP,FP,DG,BRACKT,STPMIN,STPMAX,INFO)

!Arguments ------------------------------------
!scalars
 integer,intent(inout)  :: INFO
 real(dp),intent(in)     :: FP
 real(dp),intent(inout)  :: STX,FX,DX,STY,FY,DY,STP,DG,STPMIN,STPMAX
 logical,intent(inout) :: BRACKT

!Local variables-------------------------------
!scalars
 logical BOUND
 real(dp) GAM,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
!***************************************************************************

 INFO = 0
!
! CHECK THE INPUT PARAMETERS FOR ERRORS.
!
 IF ((BRACKT .AND. (STP <= MIN(STX,STY) .OR. &
     STP >= MAX(STX,STY))) .OR.  &
     DX*(STP-STX) >= 0.0 .OR. STPMAX < STPMIN) RETURN
!
! Determine if the derivatives have opposite sign
!
 SGND = DG * ( DX / ABS(DX) )

! FIRST CASE. A HIGHER FUNCTION VALUE.
! THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
! TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
! ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
!
 IF (FP > FX) THEN
    INFO = 1
    BOUND = .TRUE.
    THETA = 3*(FX - FP)/(STP - STX) + DX + DG
    S = MAX(ABS(THETA),ABS(DX),ABS(DG))
    GAM = S * SQRT( (THETA/S)**2 - (DX/S)*(DG/S) )
    IF (STP < STX) GAM = -GAM
    P = (GAM - DX) + THETA
    Q = ((GAM - DX) + GAM) + DG
    R = P / Q
    STPC = STX + R*(STP - STX)
    STPQ = STX + ( ( DX / ( ( FX - FP ) / ( STP - STX ) + DX ) ) / 2 ) * ( STP - STX )
    IF (ABS(STPC-STX) < ABS(STPQ-STX)) THEN
       STPF = STPC
    ELSE
      STPF = STPC + (STPQ - STPC) / 2
    END IF
    BRACKT = .TRUE.
!
! SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
! OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
! STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
! THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
!
 ELSE IF (SGND < 0.0) THEN
    INFO = 2
    BOUND = .FALSE.
    THETA = 3*(FX - FP)/(STP - STX) + DX + DG
    S = MAX(ABS(THETA),ABS(DX),ABS(DG))
    GAM = S * SQRT( (THETA/S)**2 - (DX/S)*(DG/S) )
    IF (STP > STX) GAM = -GAM
    P = (GAM - DG) + THETA
    Q = ((GAM - DG) + GAM) + DX
    R = P/Q
    STPC = STP + R*(STX - STP)
    STPQ = STP + (DG/(DG-DX))*(STX - STP)
    IF (ABS(STPC-STP) > ABS(STPQ-STP)) THEN
       STPF = STPC
    ELSE
       STPF = STPQ
    END IF
    BRACKT = .TRUE.
!
! THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
! SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
! THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
! IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
! IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
! EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
! COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
! CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
!
 ELSE IF (ABS(DG) < ABS(DX)) THEN
    INFO = 3
    BOUND = .TRUE.
    THETA = 3*(FX - FP)/(STP - STX) + DX + DG
    S = MAX(ABS(THETA),ABS(DX),ABS(DG))
!
!   THE CASE GAM = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
!   TO INFINITY IN THE DIRECTION OF THE STEP.
!
    GAM = S * SQRT( MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DG/S)) )
    IF (STP > STX) GAM = -GAM
    P = (GAM - DG) + THETA
    Q = (GAM + (DX - DG)) + GAM
    R = P/Q
    IF (R < 0.0 .AND. GAM .NE. 0.0) THEN
       STPC = STP + R*(STX - STP)
    ELSE IF (STP > STX) THEN
       STPC = STPMAX
    ELSE
       STPC = STPMIN
    END IF
    STPQ = STP + (DG/(DG-DX))*(STX - STP)
    IF (BRACKT) THEN
       IF (ABS(STP-STPC) < ABS(STP-STPQ)) THEN
          STPF = STPC
       ELSE
          STPF = STPQ
       END IF
    ELSE
       IF (ABS(STP-STPC) > ABS(STP-STPQ)) THEN
          STPF = STPC
       ELSE
          STPF = STPQ
       END IF
    END IF
!
! FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
! SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
! NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
! IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
!
 ELSE
    INFO = 4
    BOUND = .FALSE.
    IF (BRACKT) THEN
       THETA = 3*(FP - FY)/(STY - STP) + DY + DG
       S = MAX(ABS(THETA),ABS(DY),ABS(DG))
       GAM = S * SQRT( (THETA/S)**2 - (DY/S)*(DG/S) )
       IF (STP > STY) GAM = -GAM
       P = (GAM - DG) + THETA
       Q = ((GAM - DG) + GAM) + DY
       R = P/Q
       STPC = STP + R*(STY - STP)
       STPF = STPC
    ELSE IF (STP > STX) THEN
       STPF = STPMAX
    ELSE
       STPF = STPMIN
    END IF
 END IF

!
! Update the interval of uncertainty. this update does not
! depend on the new step or the case analysis above.
!
 IF (FP > FX) THEN
    STY = STP
    FY = FP
    DY = DG
 ELSE
    IF (SGND < 0.0) THEN
       STY = STX
       FY = FX
       DY = DX
    END IF
    STX = STP
    FX = FP
    DX = DG
 END IF

!
! Compute the new step and safeguard it.
!
 STPF = MIN(STPMAX,STPF)
 STPF = MAX(STPMIN,STPF)
 STP = STPF
 IF (BRACKT .AND. BOUND) THEN
    IF (STY > STX) THEN
       STP = MIN( STX + 0.66 * (STY-STX) , STP)
    ELSE
       STP = MAX( STX + 0.66 * (STY-STX) , STP)
    END IF
 END IF


end subroutine mcstep
!!***

end module m_lbfgs
!!***

MODULE BoundedLeastSquares

! Summary:
 
! BVLS solves linear least-squares problems with upper and lower bounds on the
! variables, using an active set strategy.  It is documented in the J. of
! Computational Statistics, and can be used iteratively to solve minimum
! l-1, l-2 and l-infinity fitting problems.


! Statement:
! This algorithm may be used freely for non-commercial purposes, and may be
! freely distributed for non-commercial purposes.  The authors do not warrant
! the software in any way: use it at your own risk.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-11-23  Time: 23:41:42


IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS

!=======================================================================

SUBROUTINE bvls(key, m, n, a, b, bl, bu, x, istate, loopa)
!=======================================================================

! N.B. Arguments W, ACT & ZZ have been removed.

INTEGER, INTENT(IN)        :: key
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(m,n)
REAL (dp), INTENT(IN OUT)  :: b(m)
REAL (dp), INTENT(IN)      :: bl(n)
REAL (dp), INTENT(IN)      :: bu(n)
REAL (dp), INTENT(OUT)     :: x(n)
INTEGER, INTENT(IN OUT)    :: istate(n+1)
INTEGER, INTENT(OUT)       :: loopa

!--------------------Bounded Variable Least Squares---------------------

!        Robert L. Parker and Philip B. Stark    Version 3/19/90

!  Robert L. Parker                           Philip B. Stark
!  Scripps Institution of Oceanography        Department of Statistics
!  University of California, San Diego        University of California
!  La Jolla CA 92093                          Berkeley CA 94720-3860
!  rlparker@ucsd.edu                          stark@stat.berkeley.edu

!  Copyright of this software is reserved by the authors; however, this
!  algorithm and subroutine may be used freely for non-commercial
!  purposes, and may be distributed freely for non-commercial purposes.

!  The authors do not warrant this software in any way: use it at your
!  own risk.


!  See the article ``Bounded Variable Least Squares:  An Algorithm and
!  Applications'' by P.B. Stark and R.L. Parker, in the journal:
!  Computational Statistics, vol.10(2), (1995) for further description and
!  applications to minimum l-1, l-2 and l-infinity fitting problems, as well
!  as finding bounds on linear functionals subject to bounds on variables and
!  fitting linear data within l-1, l-2 or l-infinity measures of misfit.

!  BVLS solves the problem:

!          min  || a.x - b ||     such that   bl <= x <= bu
!                            2
!    where
!               x  is an unknown n-vector
!               a  is a given m by n matrix
!               b  is a given  m-vector
!               bl is a given n-vector of lower bounds on the components of x.
!               bu is a given n-vector of upper bounds on the components of x.


!-----------------------------------------------------------------------
!    Input parameters:

!  m, n, a, b, bl, bu   see above.   Let mm=min(m,n).

!  If key = 0, the subroutine solves the problem from scratch.

!  If key > 0 the routine initializes using the user's guess about which
!   components of  x  are `active', i.e. are stricly within their bounds,
!   which are at their lower bounds, and which are at their upper bounds.
!   This information is supplied through the array istate.  istate(n+1)
!   should contain the total number of components at their bounds (the
!   `bound variables').  The absolute values of the first nbound=istate(n+1)
!   entries of  istate  are the indices of these `bound' components of  x.
!   The sign of istate(j), j=1,..., nbound, indicates whether  x(|istate(j)|)
!   is at its upper or lower bound.  istate(j) is positive if the component
!   is at its upper bound, negative if the component is at its lower bound.
!   istate(j), j=nbound+1,...,n  contain the indices of the components of  x
!   that are active (i.e. are expected to lie strictly within their bounds).
!   When key > 0, the routine initially sets the active components to the
!   averages of their upper and lower bounds:
!   x(j)=(bl(j)+bu(j))/2, for j in the active set.

!-----------------------------------------------------------------------
!    Output parameters:

!  x       the solution vector.

!  w(1)    the minimum 2-norm || a.x-b ||.

!  istate  vector indicating which components of  x  are active and which are
!          at their bounds (see the previous paragraph).
!          istate can be supplied to the routine to give it a good starting
!          guess for the solution.

!  loopA   number of iterations taken in the main loop, Loop A.

!-----------------------------------------------------------------------
!    Working  arrays:

!  w      dimension n.               act      dimension m*(mm+2).
!  zz     dimension m.               istate   dimension n+1.

!-----------------------------------------------------------------------
!  Method: active variable method along the general plan of NNLS by
!  Lawson & Hanson, "Solving Least Squares Problems", Prentice-Hall 1974.
!  See Algorithm 23.10.  Step numbers in comment statements refer to their
!  scheme.
!  For more details and further uses, see the article
!  "Bounded Variable Least Squares:  An Algorithm and Applications"
!  by Stark and Parker in 1995 Computational Statistics.

!-----------------------------------------------------------------------
!  A number of measures are taken to enhance numerical reliability:

! 1. As noted by Lawson and Hanson, roundoff errors in the computation of the
!   gradient of the misfit may cause a component on the bounds to appear to
!   want to become active, yet when the component is added to the active set,
!   it moves away from the feasible region.  In this case the component is not
!   made active, the gradient of the misfit with respect to a change in that
!   component is set to zero, and the program returns to the Kuhn-Tucker test.
!   Flag  ifrom5  is used in this test, which occurs at the end of Step 6.


! 2. When the least-squares minimizer after Step 6 is infeasible, it is used
!   in a convex interpolation with the previous solution to obtain a feasible
!   vector.  The constant in this interpolation is supposed to put at least
!   one component of  x  on a bound.  There can be difficulties:

! 2a. Sometimes, due to roundoff, no interpolated component ends up on
!   a bound.  The code in Step 11 uses the flag  jj, computed in Step 8,
!   to ensure that at least the component that determined the
!   interpolation constant  alpha  is moved to the appropriate bound.
!   This guarantees that what Lawson and Hanson call `Loop B' is finite.

! 2b. The code in Step 11 also incorporates Lawson and Hanson's feature that
!   any components remaining infeasible at this stage (which must be due to
!   roundoff) are moved to their nearer bound.


! 3. If the columns of  a  passed to qr are linearly dependent, the new
!   potentially active component is not introduced: the gradient of the
!   misfit with respect to that component is set to zero, and control
!   returns to the Kuhn-Tucker test.


! 4. When some of the columns of  a  are approximately linearly dependent,
!   we have observed cycling of active components: a component just moved
!   to a bound desires immediately to become active again; qr allows it
!   to become active and a different component is moved to its bound.
!   This component immediately wants to become active, which qr allows, and
!   the original component is moved back to its bound.  We have taken two
!   steps to avoid this problem:

! 4a. First, the column of the matrix  a  corresponding to the new
!   potentially active component is passed to qr as the last column of
!   its matrix.  This ordering tends to make a component recently moved
!   to a bound fail the test mentioned in (1), above.

! 4b. Second, we have incorporated a test that prohibits short cycles.
!   If the most recent successful change to the active set was to move
!   the component x(jj) to a bound, x(jj) is not permitted to reenter
!   the solution at this stage.  This test occurs just after checking
!   the Kuhn-Tucker conditions, and uses the flag  jj, set in Step 8.
!   The flag  jj  is reset after Step 6 if Step 6 was entered from
!   Step 5 indicating that a new component has successfully entered the
!   active set.  The test for resetting  jj  uses the flag  ifrom5,
!   which will not equal zero in case Step 6 was entered from Step 5.

!     dimension w(n), act(m,min(m,n)+2), zz(m), istate(n+1)

! Local variables

REAL (dp), PARAMETER  :: eps = 1.0E-11_dp
INTEGER               :: i, iact, ifrom5, it, j, jj, k, k1, kk, ks, mm, mm1, &
                         nact, nbound, noldb
REAL (dp)             :: act(m,m+2), alf, alpha, bad, bdiff, bnorm, bound, &
                         bsq, obj, resq, ri, sj, w(n), worst, zz(m)

!----------------------First Executable Statement-----------------------

!  Step 1.  Initialize everything--active and bound sets, initial values, etc.

!  Initialize flags, etc.
mm = MIN(m,n)
mm1 = mm + 1
jj = 0
ifrom5 = 0

!  Check consistency of given bounds  bl, bu.
bdiff = 0.0_dp
DO  j = 1, n
  bdiff = MAX(bdiff, bu(j)-bl(j))
  IF (bl(j) > bu(j)) THEN
    WRITE(*, *) ' Inconsistent bounds in BVLS. '
    STOP
  END IF
END DO
IF (bdiff == 0.0_dp) THEN
  WRITE(*, *) ' No free variables in BVLS--check input bounds.'
  STOP
END IF

!  In a fresh initialization (key = 0) bind all variables at their lower bounds.
!  If (key != 0), use the supplied  istate  vector to initialize the variables.
!  istate(n+1) contains the number of bound variables.  The absolute values of
!  the first nbound=istate(n+1) entries of  istate  are the indices of the
!  bound variables.  The sign of each entry determines whether the indicated
!  variable is at its upper (positive) or lower (negative) bound.
IF (key == 0) THEN
  nbound = n
  nact = 0
  DO  j = 1, nbound
    istate(j) = -j
  END DO
ELSE
  nbound = istate(n+1)
END IF
nact = n - nbound
IF (nact > mm) THEN
  WRITE(*, *) ' Too many active variables in BVLS starting solution!'
  STOP
END IF
DO  k = 1, nbound
  j = ABS(istate(k))
  IF (istate(k) < 0) x(j) = bl(j)
  IF (istate(k) > 0) x(j) = bu(j)
END DO

!  In a warm start (key != 0) initialize the active variables to (bl+bu)/2.
!   This is needed in case the initial qr results in active variables
!   out-of-bounds and Steps 8-11 get executed the first time through.
DO  k = nbound + 1, n
  kk = istate(k)
  x(kk) = (bu(kk)+bl(kk)) / 2
END DO

!  Compute bnorm, the norm of the data vector b, for reference.
bsq = SUM( b(1:m)**2 )
bnorm = SQRT(bsq)

!-----------------------------Main Loop---------------------------------

!  Initialization complete.  Begin major loop (Loop A).
DO  loopa = 1, 3 * n
  
!  Step 2.
!  Initialize the negative gradient vector w(*).
  obj = 0.0_dp
  w(1:n) = 0.0_dp
  
!  Compute the residual vector b-a.x , the negative gradient vector
!   w(*), and the current objective value obj = || a.x - b ||.
!   The residual vector is stored in the mm+1'st column of act(*,*).
  DO  i = 1, m
    ri = b(i) - DOT_PRODUCT( a(i,1:n), x(1:n) )
    obj = obj + ri ** 2
    w(1:n) = w(1:n) + a(i,1:n) * ri
    act(i,mm1) = ri
  END DO
  
!  Converged?  Stop if the misfit << || b ||, or if all components are
!   active (unless this is the first iteration from a warm start).
  IF (SQRT(obj) <= bnorm*eps.OR.(loopa > 1.AND.nbound == 0)) THEN
    istate(n+1) = nbound
    w(1) = SQRT(obj)
    RETURN
  END IF
  
!  Add the contribution of the active components back into the residual.
  DO  k = nbound + 1, n
    j = istate(k)
    act(1:m,mm1) = act(1:m,mm1) + a(1:m,j) * x(j)
  END DO
  
!  The first iteration in a warm start requires immediate qr.
  IF (loopa == 1 .AND. key /= 0) GO TO 150
  
!  Steps 3, 4.
!  Find the bound element that most wants to be active.
  120 worst = 0.0_dp
  it = 1
  DO  j = 1, nbound
    ks = ABS(istate(j))
    bad = w(ks) * SIGN(1,istate(j))
    IF (bad < worst) THEN
      it = j
      worst = bad
      iact = ks
    END IF
  END DO
  
!  Test whether the Kuhn-Tucker condition is met.
  IF (worst >= 0.0_dp) THEN
    istate(n+1) = nbound
    w(1) = SQRT(obj)
    RETURN
  END IF
  
!  The component  x(iact)  is the one that most wants to become active.
!   If the last successful change in the active set was to move x(iact)
!   to a bound, don't let x(iact) in now: set the derivative of the misfit
!   with respect to x(iact) to zero and return to the Kuhn-Tucker test.
  IF (iact == jj) THEN
    w(jj) = 0.0_dp
    GO TO 120
  END IF
  
!  Step 5.
!  Undo the effect of the new (potentially) active variable on the
!   residual vector.
  IF (istate(it) > 0) bound = bu(iact)
  IF (istate(it) < 0) bound = bl(iact)
  act(1:m,mm1) = act(1:m,mm1) + bound * a(1:m,iact)
  
!  Set flag ifrom5, indicating that Step 6 was entered from Step 5.
!   This forms the basis of a test for instability: the gradient calculation
!   shows that x(iact) wants to join the active set; if qr puts x(iact) beyond
!   the bound from which it came, the gradient calculation was in error and
!   the variable should not have been introduced.
  ifrom5 = istate(it)
  
!  Swap the indices (in istate) of the new active variable and the
!   rightmost bound variable; `unbind' that location by decrementing nbound.
  istate(it) = istate(nbound)
  nbound = nbound - 1
  nact = nact + 1
  istate(nbound+1) = iact
  
  IF (mm < nact) THEN
    WRITE(*, *) ' Too many free variables in BVLS!'
    STOP
  END IF
  
!  Step 6.
!  Load array  act  with the appropriate columns of  a  for qr.  For added
!   stability, reverse the column ordering so that the most recent addition to
!   the active set is in the last column.  Also copy the residual vector from
!   act(., mm1) into act(., mm1+1).
  150 DO  i = 1, m
    act(i,mm1+1) = act(i,mm1)
    DO  k = nbound + 1, n
      j = istate(k)
      act(i, nact+1-k+nbound) = a(i,j)
    END DO
  END DO
  
  CALL qr(m, nact, act, act(:,mm1+1), zz, resq)
  
!  Test for linear dependence in qr, and for an instability that moves the
!   variable just introduced away from the feasible region (rather than into
!   the region or all the way through it).
!   In either case, remove the latest vector introduced from the active set
!   and adjust the residual vector accordingly.
!   Set the gradient component (w(iact)) to zero and return to the Kuhn-Tucker
!   test.
  IF (resq < 0.0_dp .OR. (ifrom5 > 0 .AND. zz(nact) > bu(iact)) .OR.   &
      (ifrom5 < 0 .AND. zz(nact) < bl(iact))) THEN
    nbound = nbound + 1
    istate(nbound) = istate(nbound) * SIGN(1.0D0, x(iact)-bu(iact))
    nact = nact - 1
    act(1:m,mm1) = act(1:m,mm1) - x(iact) * a(1:m,iact)
    ifrom5 = 0
    w(iact) = 0.0_dp
    GO TO 120
  END IF
  
!  If Step 6 was entered from Step 5 and we are here, a new variable
!   has been successfully introduced into the active set; the last
!   variable that was fixed at a bound is again permitted to become active.
  IF (ifrom5 /= 0) jj = 0
  ifrom5 = 0
  
!   Step 7.  Check for strict feasibility of the new qr solution.
  DO  k = 1, nact
    k1 = k
    j = istate(k+nbound)
    IF (zz(nact+1-k) < bl(j) .OR. zz(nact+1-k) > bu(j)) GO TO 210
  END DO
  DO  k = 1, nact
    j = istate(k+nbound)
    x(j) = zz(nact+1-k)

  END DO
!  New iterate is feasible; back to the top.
  CYCLE
  
!  Steps 8, 9.
  210 alpha = 2.0_dp
  alf = alpha
  DO  k = k1, nact
    j = istate(k+nbound)
    IF (zz(nact+1-k) > bu(j)) alf = (bu(j)-x(j)) / ( zz(nact+1-k)-x(j))
    IF (zz(nact+1-k) < bl(j)) alf = (bl(j)-x(j)) / ( zz(nact+1-k)-x(j))
    IF (alf < alpha) THEN
      alpha = alf
      jj = j
      sj = SIGN(1.0_dp, zz(nact+1-k)-bl(j))
    END IF
  END DO
  
!  Step 10
  DO  k = 1, nact
    j = istate(k+nbound)
    x(j) = x(j) + alpha * (zz(nact+1-k)-x(j))
  END DO
  
!  Step 11.
!  Move the variable that determined alpha to the appropriate bound.
!   (jj is its index; sj is + if zz(jj)> bu(jj), - if zz(jj)<bl(jj) ).
!   If any other component of  x  is infeasible at this stage, it must
!   be due to roundoff.  Bind every infeasible component and every
!   component at a bound to the appropriate bound.  Correct the
!   residual vector for any variables moved to bounds.  Since at least
!   one variable is removed from the active set in this step, Loop B
!   (Steps 6-11) terminates after at most  nact  steps.
  noldb = nbound
  DO  k = 1, nact
    j = istate(k+noldb)
    IF (bu(j)-x(j) <= 0.0_dp .OR. (j == jj .AND. sj > 0.0_dp)) THEN

!  Move x(j) to its upper bound.
      x(j) = bu(j)
      istate(k+noldb) = istate(nbound+1)
      istate(nbound+1) = j
      nbound = nbound + 1
      act(1:m,mm1) = act(1:m,mm1) - bu(j) * a(1:m,j)
    ELSE IF (x(j)-bl(j) <= 0.0_dp .OR. (j == jj .AND. sj < 0.0_dp)) THEN

!  Move x(j) to its lower bound.
      x(j) = bl(j)
      istate(k+noldb) = istate(nbound+1)
      istate(nbound+1) = -j
      nbound = nbound + 1
      act(1:m,mm1) = act(1:m,mm1) - bl(j) * a(1:m,j)
    END IF
  END DO
  nact = n - nbound
  
!  If there are still active variables left repeat the qr; if not,
!    go back to the top.
  IF (nact > 0) GO TO 150
  
END DO

WRITE(*, *) ' BVLS fails to converge! '
STOP
END SUBROUTINE bvls

!======================================================================

SUBROUTINE qr(m, n, a, b, x, resq)
!======================================================================

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(m,n)
REAL (dp), INTENT(IN OUT)  :: b(m)
REAL (dp), INTENT(OUT)     :: x(n)
REAL (dp), INTENT(OUT)     :: resq

!$$$$ calls no other routines
!  Relies on FORTRAN77 do-loop conventions!
!  Solves over-determined least-squares problem  ax ~ b
!  where  a  is an  m by n  matrix,  b  is an m-vector.
!  resq  is the sum of squared residuals of optimal solution.  Also used
!  to signal error conditions - if -2 , system is underdetermined, if
!  -1,  system is singular.
!  Method - successive Householder rotations.  See Lawson & Hanson -
!  Solving Least Squares Problems (1974).
!  Routine will also work when m=n.
!*****   CAUTION -  a and b  are overwritten by this routine.

! Local variables

REAL (dp) :: const, total, dot, qv1, sq, u1
INTEGER   :: i, j, j1, jj

resq = -2.0
IF (m < n) RETURN
resq = -1.0

!   Loop ending on 1800 rotates  a  into upper triangular form.
DO  j = 1, n

!   Find constants for rotation and diagonal entry.
  sq = SUM( a(j:m,j) ** 2 )
  IF (sq == 0.0_dp) RETURN
  qv1 = -SIGN(SQRT(sq), a(j,j))
  u1 = a(j,j) - qv1
  a(j,j) = qv1
  j1 = j + 1

!  Rotate remaining columns of sub-matrix.
  DO  jj = j1, n
    dot = u1 * a(j,jj) + DOT_PRODUCT( a(j1:m,jj), a(j1:m,j) )
    const = dot / ABS(qv1*u1)
    a(j1:m,jj) = a(j1:m,jj) - const * a(j1:m,j)
    a(j,jj) = a(j,jj) - const * u1
  END DO

!  Rotate  b  vector.
  dot = u1 * b(j) + DOT_PRODUCT( a(j1:m,j), b(j1:m) )
  const = dot / ABS(qv1*u1)
  b(j) = b(j) - const * u1
  b(j1:m) = b(j1:m) - const * a(j1:m,j)
END DO

!  Solve triangular system by back-substitution.
DO  i = n, 1, -1
  total = b(i) - DOT_PRODUCT( a(i,i+1:n), x(i+1:n) )
  IF (a(i,i) == 0.0_dp) RETURN
  x(i) = total / a(i,i)
END DO

!  Find residual in overdetermined case.
resq = SUM( b(n+1:m) ** 2 )

RETURN
END SUBROUTINE qr
!______________________________________________________________________

END MODULE BoundedLeastSquares



!PROGRAM Test_BVLS
!USE BoundedLeastSquares

!IMPLICIT NONE
!INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
!
!INTEGER, PARAMETER  :: ncases = 100, ncols = 10
!REAL (dp)           :: a(ncases,ncols), y(ncases), bl(ncols), bu(ncols),  &
!                       beta(ncols), x(ncols), e
!INTEGER             :: case, istate(ncols+1), j, key, loopa

! Generate artificial data satisfying
!   Y = A.beta + noise

!CALL RANDOM_NUMBER(beta)
!beta = 4.0*(beta - 0.5)
!CALL RANDOM_NUMBER(a)
!DO case = 1, ncases
!  CALL RANDOM_NUMBER(e)
!  y(case) = DOT_PRODUCT( a(case, :), beta ) + 0.1*(e - 0.5)
!END DO
!key = 0
!bl = 0.0_dp
!bu = 1.0_dp
!CALL bvls(key, ncases, ncols, a, y, bl, bu, x, istate, loopa)
!
!WRITE(*, *) ' Column   Original beta   Solution'
!DO j = 1, ncols
!  WRITE(*, '(i5, 2f14.3)') j, beta(j), x(j)
!END DO
!WRITE(*, '(a, i4)') ' No. of iterations = ', loopa
!WRITE(*, *)

!STOP
!END PROGRAM Test_BVLS


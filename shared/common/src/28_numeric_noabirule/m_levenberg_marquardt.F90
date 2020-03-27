!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_levenberg_marquardt
!! NAME
!! m_levenberg_marquardt
!!
!! FUNCTION
!!  Module containing routines for performing Levenberg-Marquardt
!!  nonlinear least-squares fit. These routines are the conversions
!!  to F90 of the Minpack F77 routines by Alan Miller
!!
!!  TODO: For some reason, the routine lmder, which uses analytic
!!        evaluation of the Jacobian, did not work in testing (bug?).
!!        lmdif proved stable, but potentially uses more funtion
!!        evaluations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_levenberg_marquardt

 use defs_basis
 use m_abicore
! MINPACK routines which are used by both LMDIF & LMDER
! 25 October 2001:
!    Changed INTENT of iflag in several places to IN OUT.
!    Changed INTENT of fvec to IN OUT in user routine FCN.
!    Removed arguments diag and qtv from LMDIF & LMDER.
!    Replaced several DO loops with array operations.
! amiller @ bigpond.net.au

 implicit none

 private

 public :: lmdif1
 public :: lmdif
 public :: lm_fit_print_info

CONTAINS  !==============================================================================
!!***

!!****f* m_levenberg_marquardt/lmdif1
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine lmdif1(fcn, m, n, x, fvec, tol, info, iwa)

 ! code converted using to_f90 by alan miller
 ! date: 1999-12-11  time: 00:51:44

 ! n.b. arguments wa & lwa have been removed.
 integer, intent(in)        :: m
 integer, intent(in)        :: n
 real (dp), intent(in out)  :: x(:)
 real (dp), intent(out)     :: fvec(:)
 real (dp), intent(in)      :: tol
 integer, intent(out)       :: info
 integer, intent(out)       :: iwa(:)

 ! external fcn

 interface
   subroutine fcn(m, n, x, fvec, iflag, y, nfqre, nfqim, bsign, csign, multi_x_exp)
     implicit none
     integer, parameter         :: dp = KIND(1.0d0)
     integer, intent(in)        :: m, n
     real (dp), intent(in)      :: x(n)
     real (dp), intent(in out)  :: fvec(m)
     integer, intent(in out)    :: iflag
     integer, optional, intent(in) :: nfqre,nfqim,bsign,csign,multi_x_exp
     real (dp), optional, intent(in)      :: y(m)
   end subroutine fcn
 end interface

 !  **********

 !  subroutine lmdif1

 !  The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear
 !  functions in n variables by a modification of the Levenberg-Marquardt
 !  algorithm.  This is done by using the more general least-squares
 !  solver lmdif.  The user must provide a subroutine which calculates the
 !  functions.  The jacobian is then calculated by a forward-difference
 !  approximation.

 !  the subroutine statement is

 !    subroutine lmdif1(fcn, m, n, x, fvec, tol, info, iwa)

 !  where

 !    fcn is the name of the user-supplied subroutine which calculates
 !      the functions.  fcn must be declared in an external statement in the
 !      user calling program, and should be written as follows.

 !      subroutine fcn(m, n, x, fvec, iflag)
 !      integer m, n, iflag
 !      REAL (dp) x(n), fvec(m)
 !      ----------
 !      calculate the functions at x and return this vector in fvec.
 !      ----------
 !      return
 !      end

 !      the value of iflag should not be changed by fcn unless
 !      the user wants to terminate execution of lmdif1.
 !      In this case set iflag to a negative integer.

 !    m is a positive integer input variable set to the number of functions.

 !    n is a positive integer input variable set to the number of variables.
 !      n must not exceed m.

 !    x is an array of length n.  On input x must contain an initial estimate
 !      of the solution vector.  On output x contains the final estimate of
 !      the solution vector.

 !    fvec is an output array of length m which contains
 !      the functions evaluated at the output x.

 !    tol is a nonnegative input variable.  Termination occurs when the
 !      algorithm estimates either that the relative error in the sum of
 !      squares is at most tol or that the relative error between x and the
 !      solution is at most tol.

 !    info is an integer output variable.  If the user has terminated execution,
 !      info is set to the (negative) value of iflag.  See description of fcn.
 !      Otherwise, info is set as follows.

 !      info = 0  improper input parameters.

 !      info = 1  algorithm estimates that the relative error
 !                in the sum of squares is at most tol.

 !      info = 2  algorithm estimates that the relative error
 !                between x and the solution is at most tol.

 !      info = 3  conditions for info = 1 and info = 2 both hold.

 !      info = 4  fvec is orthogonal to the columns of the
 !                jacobian to machine precision.

 !      info = 5  number of calls to fcn has reached or exceeded 200*(n+1).

 !      info = 6  tol is too small. no further reduction in
 !                the sum of squares is possible.

 !      info = 7  tol is too small.  No further improvement in
 !                the approximate solution x is possible.

 !    iwa is an integer work array of length n.

 !    wa is a work array of length lwa.

 !    lwa is a positive integer input variable not less than m*n+5*n+m.

 !  subprograms called

 !    user-supplied ...... fcn

 !    minpack-supplied ... lmdif

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********
 integer   :: maxfev, mode, nfev, nprint
 real (dp) :: epsfcn, ftol, gtol, xtol, fjac(m,n)
 real (dp), parameter :: factor = 100._dp , zero = 0.0_dp

 info = 0

 !     check the input parameters for errors.

 IF (n <= 0 .OR. m < n .OR. tol < zero) GO TO 10

 !     call lmdif.

 maxfev = 200*(n + 1)
 ftol = tol
 xtol = tol
 gtol = zero
 epsfcn = zero
 mode = 1
 nprint = 0
 CALL lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,   &
&           mode, factor, nprint, info, nfev, fjac, iwa)
 IF (info == 8) info = 4

 10 RETURN

 !     last card of subroutine lmdif1.

 end subroutine lmdif1
!!***

!----------------------------------------------------------------------

!!****f* m_levenberg_marquardt/lmdif
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,  &
                  mode, factor, nprint, info, nfev, fjac, ipvt)

 ! N.B. Arguments LDFJAC, DIAG, QTF, WA1, WA2, WA3 & WA4 have been removed.
 integer, intent(in)        :: m
 integer, intent(in)        :: n
 real (dp), intent(in out)  :: x(:)
 real (dp), intent(out)     :: fvec(:)
 real (dp), intent(in)      :: ftol
 real (dp), intent(in)      :: xtol
 real (dp), intent(in out)  :: gtol
 integer, intent(in out)    :: maxfev
 real (dp), intent(in out)  :: epsfcn
 integer, intent(in)        :: mode
 real (dp), intent(in)      :: factor
 integer, intent(in)        :: nprint
 integer, intent(out)       :: info
 integer, intent(out)       :: nfev
 real (dp), intent(out)     :: fjac(:,:)    ! fjac(ldfjac,n)
 integer, intent(out)       :: ipvt(:)

 ! external fcn

 interface
   subroutine fcn(m, n, x, fvec, iflag, y, nfqre, nfqim, bsign, csign, multi_x_exp)
     implicit none
     integer, parameter         :: dp = KIND(1.0d0)
     integer, intent(in)        :: m, n
     real (dp), intent(in)      :: x(n)
     real (dp), intent(in out)  :: fvec(m)
     integer, intent(in out)    :: iflag
     integer, optional, intent(in) :: nfqre,nfqim,bsign,csign,multi_x_exp
     real (dp), optional, intent(in)      :: y(m)
   end subroutine fcn
 end interface

 !  **********

 !  subroutine lmdif

 !  The purpose of lmdif is to minimize the sum of the squares of m nonlinear
 !  functions in n variables by a modification of the Levenberg-Marquardt
 !  algorithm.  The user must provide a subroutine which calculates the
 !  functions.  The jacobian is then calculated by a forward-difference
 !  approximation.

 !  the subroutine statement is

 !    subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
 !                     diag, mode, factor, nprint, info, nfev, fjac,
 !                     ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)

 ! N.B. 7 of these arguments have been removed in this version.

 !  where

 !    fcn is the name of the user-supplied subroutine which calculates the
 !      functions.  fcn must be declared in an external statement in the user
 !      calling program, and should be written as follows.

 !      subroutine fcn(m, n, x, fvec, iflag)
 !      integer m, n, iflag
 !      REAL (dp) x(:), fvec(m)
 !      ----------
 !      calculate the functions at x and return this vector in fvec.
 !      ----------
 !      return
 !      end

 !      the value of iflag should not be changed by fcn unless
 !      the user wants to terminate execution of lmdif.
 !      in this case set iflag to a negative integer.

 !    m is a positive integer input variable set to the number of functions.

 !    n is a positive integer input variable set to the number of variables.
 !      n must not exceed m.

 !    x is an array of length n.  On input x must contain an initial estimate
 !      of the solution vector.  On output x contains the final estimate of the
 !      solution vector.

 !    fvec is an output array of length m which contains
 !      the functions evaluated at the output x.

 !    ftol is a nonnegative input variable.  Termination occurs when both the
 !      actual and predicted relative reductions in the sum of squares are at
 !      most ftol.  Therefore, ftol measures the relative error desired
 !      in the sum of squares.

 !    xtol is a nonnegative input variable.  Termination occurs when the
 !      relative error between two consecutive iterates is at most xtol.
 !      Therefore, xtol measures the relative error desired in the approximate
 !      solution.

 !    gtol is a nonnegative input variable.  Termination occurs when the cosine
 !      of the angle between fvec and any column of the jacobian is at most
 !      gtol in absolute value.  Therefore, gtol measures the orthogonality
 !      desired between the function vector and the columns of the jacobian.

 !    maxfev is a positive integer input variable.  Termination occurs when the
 !      number of calls to fcn is at least maxfev by the end of an iteration.

 !    epsfcn is an input variable used in determining a suitable step length
 !      for the forward-difference approximation.  This approximation assumes
 !      that the relative errors in the functions are of the order of epsfcn.
 !      If epsfcn is less than the machine precision, it is assumed that the
 !      relative errors in the functions are of the order of the machine
 !      precision.

 !    diag is an array of length n.  If mode = 1 (see below), diag is
 !      internally set.  If mode = 2, diag must contain positive entries that
 !      serve as multiplicative scale factors for the variables.

 !    mode is an integer input variable.  If mode = 1, the variables will be
 !      scaled internally.  If mode = 2, the scaling is specified by the input
 !      diag. other values of mode are equivalent to mode = 1.

 !    factor is a positive input variable used in determining the initial step
 !      bound.  This bound is set to the product of factor and the euclidean
 !      norm of diag*x if nonzero, or else to factor itself.  In most cases
 !      factor should lie in the interval (.1,100.). 100. is a generally
 !      recommended value.

 !    nprint is an integer input variable that enables controlled printing of
 !      iterates if it is positive.  In this case, fcn is called with iflag = 0
 !      at the beginning of the first iteration and every nprint iterations
 !      thereafter and immediately prior to return, with x and fvec available
 !      for printing.  If nprint is not positive, no special calls
 !      of fcn with iflag = 0 are made.

 !    info is an integer output variable.  If the user has terminated
 !      execution, info is set to the (negative) value of iflag.
 !      See description of fcn.  Otherwise, info is set as follows.

 !      info = 0  improper input parameters.

 !      info = 1  both actual and predicted relative reductions
 !                in the sum of squares are at most ftol.

 !      info = 2  relative error between two consecutive iterates <= xtol.

 !      info = 3  conditions for info = 1 and info = 2 both hold.

 !      info = 4  the cosine of the angle between fvec and any column of
 !                the Jacobian is at most gtol in absolute value.

 !      info = 5  number of calls to fcn has reached or exceeded maxfev.

 !      info = 6  ftol is too small. no further reduction in
 !                the sum of squares is possible.

 !      info = 7  xtol is too small. no further improvement in
 !                the approximate solution x is possible.

 !      info = 8  gtol is too small. fvec is orthogonal to the
 !                columns of the jacobian to machine precision.

 !    nfev is an integer output variable set to the number of calls to fcn.

 !    fjac is an output m by n array. the upper n by n submatrix
 !      of fjac contains an upper triangular matrix r with
 !      diagonal elements of nonincreasing magnitude such that

 !             t     t           t
 !            p *(jac *jac)*p = r *r,

 !      where p is a permutation matrix and jac is the final calculated
 !      Jacobian.  Column j of p is column ipvt(j) (see below) of the
 !      identity matrix. the lower trapezoidal part of fjac contains
 !      information generated during the computation of r.

 !    ldfjac is a positive integer input variable not less than m
 !      which specifies the leading dimension of the array fjac.

 !    ipvt is an integer output array of length n.  ipvt defines a permutation
 !      matrix p such that jac*p = q*r, where jac is the final calculated
 !      jacobian, q is orthogonal (not stored), and r is upper triangular
 !      with diagonal elements of nonincreasing magnitude.
 !      Column j of p is column ipvt(j) of the identity matrix.

 !    qtf is an output array of length n which contains
 !      the first n elements of the vector (q transpose)*fvec.

 !    wa1, wa2, and wa3 are work arrays of length n.

 !    wa4 is a work array of length m.

 !  subprograms called

 !    user-supplied ...... fcn

 !    minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac

 !    fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********
 integer   :: i, iflag, iter, j, l
 real (dp) :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm
 real(dp)  :: par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm
 real (dp) :: diag(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m)
 real (dp), parameter :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp
 real(dp),parameter ::   p25 = 0.25_dp, p75 = 0.75_dp, p0001 = 0.0001_dp, zero = 0.0_dp

 !     epsmch is the machine precision.

 epsmch = EPSILON(zero)

 info = 0
 iflag = 0
 nfev = 0

 !     check the input parameters for errors.

 IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
&     .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
 IF (mode /= 2) GO TO 20
 DO  j = 1, n
   IF (diag(j) <= zero) GO TO 300
 END DO

 !     evaluate the function at the starting point and calculate its norm.

 20 iflag = 1
 CALL fcn(m, n, x, fvec, iflag)
 nfev = 1
 IF (iflag < 0) GO TO 300
 fnorm = enorm(m, fvec)

 !     initialize levenberg-marquardt parameter and iteration counter.

 par = zero
 iter = 1

 !     beginning of the outer loop.

 !        calculate the jacobian matrix.

 30 iflag = 2
 CALL fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
 nfev = nfev + n
 IF (iflag < 0) GO TO 300

 !        If requested, call fcn to enable printing of iterates.

 IF (nprint <= 0) GO TO 40
 iflag = 0
 IF (MOD(iter-1,nprint) == 0) THEN
   CALL fcn(m, n, x, fvec, iflag)
 END IF
 IF (iflag < 0) GO TO 300

 !        Compute the qr factorization of the jacobian.

 40 CALL qrfac(m, n, fjac, .true., ipvt, wa1, wa2)

 !        On the first iteration and if mode is 1, scale according
 !        to the norms of the columns of the initial jacobian.

 IF (iter /= 1) GO TO 80
 IF (mode == 2) GO TO 60
 DO  j = 1, n
   diag(j) = wa2(j)
   IF (wa2(j) == zero) diag(j) = one
 END DO

 !        On the first iteration, calculate the norm of the scaled x
 !        and initialize the step bound delta.

 60 wa3(1:n) = diag(1:n)*x(1:n)
 xnorm = enorm(n, wa3)
 delta = factor*xnorm
 IF (delta == zero) delta = factor

 !        Form (q transpose)*fvec and store the first n components in qtf.

 80 wa4(1:m) = fvec(1:m)
 DO  j = 1, n
   IF (fjac(j,j) == zero) GO TO 120
   sum = DOT_PRODUCT( fjac(j:m,j), wa4(j:m) )
   temp = -sum/fjac(j,j)
   DO  i = j, m
     wa4(i) = wa4(i) + fjac(i,j)*temp
   END DO
   120 fjac(j,j) = wa1(j)
   qtf(j) = wa4(j)
 END DO

 !        compute the norm of the scaled gradient.

 gnorm = zero
 IF (fnorm == zero) GO TO 170
 DO  j = 1, n
   l = ipvt(j)
   IF (wa2(l) == zero) CYCLE
   sum = zero
   DO  i = 1, j
     sum = sum + fjac(i,j)*(qtf(i)/fnorm)
   END DO
   gnorm = MAX(gnorm, ABS(sum/wa2(l)))
 END DO

 !        test for convergence of the gradient norm.

 170 IF (gnorm <= gtol) info = 4
 IF (info /= 0) GO TO 300

 !        rescale if necessary.

 IF (mode == 2) GO TO 200
 DO  j = 1, n
   diag(j) = MAX(diag(j), wa2(j))
 END DO

 !        beginning of the inner loop.

 !           determine the Levenberg-Marquardt parameter.

 200 CALL lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

 !           store the direction p and x + p. calculate the norm of p.

 DO  j = 1, n
   wa1(j) = -wa1(j)
   wa2(j) = x(j) + wa1(j)
   wa3(j) = diag(j)*wa1(j)
 END DO
 pnorm = enorm(n, wa3)

 !           on the first iteration, adjust the initial step bound.

 IF (iter == 1) delta = MIN(delta, pnorm)

 !           evaluate the function at x + p and calculate its norm.

 iflag = 1
 CALL fcn(m, n, wa2, wa4, iflag)
 nfev = nfev + 1
 IF (iflag < 0) GO TO 300
 fnorm1 = enorm(m, wa4)

 !           compute the scaled actual reduction.

 actred = -one
 IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

 !           Compute the scaled predicted reduction and
 !           the scaled directional derivative.

 DO  j = 1, n
   wa3(j) = zero
   l = ipvt(j)
   temp = wa1(l)
   DO  i = 1, j
     wa3(i) = wa3(i) + fjac(i,j)*temp
   END DO
 END DO
 temp1 = enorm(n,wa3)/fnorm
 temp2 = (SQRT(par)*pnorm)/fnorm
 prered = temp1**2 + temp2**2/p5
 dirder = -(temp1**2 + temp2**2)

 !           compute the ratio of the actual to the predicted reduction.

 ratio = zero
 IF (prered /= zero) ratio = actred/prered

 !           update the step bound.

 IF (ratio <= p25) THEN
   IF (actred >= zero) temp = p5
   IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
   IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
   delta = temp*MIN(delta,pnorm/p1)
   par = par/temp
 ELSE
   IF (par /= zero .AND. ratio < p75) GO TO 260
   delta = pnorm/p5
   par = p5*par
 END IF

 !           test for successful iteration.

 260 IF (ratio < p0001) GO TO 290

 !           successful iteration. update x, fvec, and their norms.

 DO  j = 1, n
   x(j) = wa2(j)
   wa2(j) = diag(j)*x(j)
 END DO
 fvec(1:m) = wa4(1:m)
 xnorm = enorm(n, wa2)
 fnorm = fnorm1
 iter = iter + 1

 !           tests for convergence.

 290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
 IF (delta <= xtol*xnorm) info = 2
 IF (ABS(actred) <= ftol .AND. prered <= ftol  &
&     .AND. p5*ratio <= one .AND. info == 2) info = 3
 IF (info /= 0) GO TO 300

 !           tests for termination and stringent tolerances.

 IF (nfev >= maxfev) info = 5
 IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
&     .AND. p5*ratio <= one) info = 6
 IF (delta <= epsmch*xnorm) info = 7
 IF (gnorm <= epsmch) info = 8
 IF (info /= 0) GO TO 300

 !           end of the inner loop. repeat if iteration unsuccessful.

 IF (ratio < p0001) GO TO 200

 !        end of the outer loop.

 GO TO 30

 !     termination, either normal or user imposed.

 300 IF (iflag < 0) info = iflag
 iflag = 0
 IF (nprint > 0) THEN
   CALL fcn(m, n, x, fvec, iflag)
 END IF
 RETURN

 !     last card of subroutine lmdif.

 end subroutine lmdif
!!***

!----------------------------------------------------------------------

!!****f* m_levenberg_marquardt/lmpar
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine lmpar(n, r, ipvt, diag, qtb, delta, par, x, sdiag)

 ! Code converted using TO_F90 by Alan Miller
 ! Date: 1999-12-09  Time: 12:46:12

 ! N.B. Arguments LDR, WA1 & WA2 have been removed.
 integer, intent(in)        :: n
 real (dp), intent(in out)  :: r(:,:)
 integer, intent(in)        :: ipvt(:)
 real (dp), intent(in)      :: diag(:)
 real (dp), intent(in)      :: qtb(:)
 real (dp), intent(in)      :: delta
 real (dp), intent(out)     :: par
 real (dp), intent(out)     :: x(:)
 real (dp), intent(out)     :: sdiag(:)

 !  **********

 !  subroutine lmpar

 !  Given an m by n matrix a, an n by n nonsingular diagonal matrix d,
 !  an m-vector b, and a positive number delta, the problem is to determine a
 !  value for the parameter par such that if x solves the system

 !        a*x = b ,     sqrt(par)*d*x = 0 ,

 !  in the least squares sense, and dxnorm is the euclidean
 !  norm of d*x, then either par is zero and

 !        (dxnorm-delta) <= 0.1*delta ,

 !  or par is positive and

 !        abs(dxnorm-delta) <= 0.1*delta .

 !  This subroutine completes the solution of the problem if it is provided
 !  with the necessary information from the r factorization, with column
 !  qpivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
 !  q has orthogonal columns, and r is an upper triangular matrix with diagonal
 !  elements of nonincreasing magnitude, then lmpar expects the full upper
 !  triangle of r, the permutation matrix p, and the first n components of
 !  (q transpose)*b.
 !  On output lmpar also provides an upper triangular matrix s such that

 !         t   t                   t
 !        p *(a *a + par*d*d)*p = s *s .

 !  s is employed within lmpar and may be of separate interest.

 !  Only a few iterations are generally needed for convergence of the algorithm.
 !  If, however, the limit of 10 iterations is reached, then the output par
 !  will contain the best value obtained so far.

 !  the subroutine statement is

 !    subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, wa1,wa2)

 !  where

 !    n is a positive integer input variable set to the order of r.

 !    r is an n by n array. on input the full upper triangle
 !      must contain the full upper triangle of the matrix r.
 !      On output the full upper triangle is unaltered, and the
 !      strict lower triangle contains the strict upper triangle
 !      (transposed) of the upper triangular matrix s.

 !    ldr is a positive integer input variable not less than n
 !      which specifies the leading dimension of the array r.

 !    ipvt is an integer input array of length n which defines the
 !      permutation matrix p such that a*p = q*r. column j of p
 !      is column ipvt(j) of the identity matrix.

 !    diag is an input array of length n which must contain the
 !      diagonal elements of the matrix d.

 !    qtb is an input array of length n which must contain the first
 !      n elements of the vector (q transpose)*b.

 !    delta is a positive input variable which specifies an upper
 !      bound on the euclidean norm of d*x.

 !    par is a nonnegative variable. on input par contains an
 !      initial estimate of the levenberg-marquardt parameter.
 !      on output par contains the final estimate.

 !    x is an output array of length n which contains the least
 !      squares solution of the system a*x = b, sqrt(par)*d*x = 0,
 !      for the output par.

 !    sdiag is an output array of length n which contains the
 !      diagonal elements of the upper triangular matrix s.

 !    wa1 and wa2 are work arrays of length n.

 !  subprograms called

 !    minpack-supplied ... dpmpar,enorm,qrsolv

 !    fortran-supplied ... ABS,MAX,MIN,SQRT

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********
 integer   :: iter, j, jm1, jp1, k, l, nsing
 real (dp) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, sum, temp
 real (dp) :: wa1(n), wa2(n)
 real (dp), parameter :: p1 = 0.1_dp, p001 = 0.001_dp, zero = 0.0_dp

 !     dwarf is the smallest positive magnitude.

 dwarf = TINY(zero)

 !     compute and store in x the gauss-newton direction. if the
 !     jacobian is rank-deficient, obtain a least squares solution.

 nsing = n
 DO  j = 1, n
   wa1(j) = qtb(j)
   IF (r(j,j) == zero .AND. nsing == n) nsing = j - 1
   IF (nsing < n) wa1(j) = zero
 END DO

 DO  k = 1, nsing
   j = nsing - k + 1
   wa1(j) = wa1(j)/r(j,j)
   temp = wa1(j)
   jm1 = j - 1
   wa1(1:jm1) = wa1(1:jm1) - r(1:jm1,j)*temp
 END DO

 DO  j = 1, n
   l = ipvt(j)
   x(l) = wa1(j)
 END DO

 !     initialize the iteration counter.
 !     evaluate the function at the origin, and test
 !     for acceptance of the gauss-newton direction.

 iter = 0
 wa2(1:n) = diag(1:n)*x(1:n)
 dxnorm = enorm(n, wa2)
 fp = dxnorm - delta
 IF (fp <= p1*delta) GO TO 220

 !     if the jacobian is not rank deficient, the newton
 !     step provides a lower bound, parl, for the zero of
 !     the function.  Otherwise set this bound to zero.

 parl = zero
 IF (nsing < n) GO TO 120
 DO  j = 1, n
   l = ipvt(j)
   wa1(j) = diag(l)*(wa2(l)/dxnorm)
 END DO
 DO  j = 1, n
   sum = DOT_PRODUCT( r(1:j-1,j), wa1(1:j-1) )
   wa1(j) = (wa1(j) - sum)/r(j,j)
 END DO
 temp = enorm(n,wa1)
 parl = ((fp/delta)/temp)/temp

 !     calculate an upper bound, paru, for the zero of the function.

 120 DO  j = 1, n
   sum = DOT_PRODUCT( r(1:j,j), qtb(1:j) )
   l = ipvt(j)
   wa1(j) = sum/diag(l)
 END DO
 gnorm = enorm(n,wa1)
 paru = gnorm/delta
 IF (paru == zero) paru = dwarf/MIN(delta,p1)

 !     if the input par lies outside of the interval (parl,paru),
 !     set par to the closer endpoint.

 par = MAX(par,parl)
 par = MIN(par,paru)
 IF (par == zero) par = gnorm/dxnorm

 !     beginning of an iteration.

 150 iter = iter + 1

 !        evaluate the function at the current value of par.

 IF (par == zero) par = MAX(dwarf, p001*paru)
 temp = SQRT(par)
 wa1(1:n) = temp*diag(1:n)
 CALL qrsolv(n, r, ipvt, wa1, qtb, x, sdiag)
 wa2(1:n) = diag(1:n)*x(1:n)
 dxnorm = enorm(n, wa2)
 temp = fp
 fp = dxnorm - delta

 !        if the function is small enough, accept the current value
 !        of par. also test for the exceptional cases where parl
 !        is zero or the number of iterations has reached 10.

 IF (ABS(fp) <= p1*delta .OR. parl == zero .AND. fp <= temp  &
&     .AND. temp < zero .OR. iter == 10) GO TO 220

 !        compute the newton correction.

 DO  j = 1, n
   l = ipvt(j)
   wa1(j) = diag(l)*(wa2(l)/dxnorm)
 END DO
 DO  j = 1, n
   wa1(j) = wa1(j)/sdiag(j)
   temp = wa1(j)
   jp1 = j + 1
   wa1(jp1:n) = wa1(jp1:n) - r(jp1:n,j)*temp
 END DO
 temp = enorm(n,wa1)
 parc = ((fp/delta)/temp)/temp

 !        depending on the sign of the function, update parl or paru.

 IF (fp > zero) parl = MAX(parl,par)
 IF (fp < zero) paru = MIN(paru,par)

 !        compute an improved estimate for par.

 par = MAX(parl, par+parc)

 !        end of an iteration.

 GO TO 150

 !     termination.

 220 IF (iter == 0) par = zero
 RETURN

 !     last card of subroutine lmpar.

 end subroutine lmpar
!!***

!----------------------------------------------------------------------

!!****f* m_levenberg_marquardt/qrfac
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine qrfac(m, n, a, pivot, ipvt, rdiag, acnorm)

 ! Code converted using TO_F90 by Alan Miller
 ! Date: 1999-12-09  Time: 12:46:17

 ! N.B. Arguments LDA, LIPVT & WA have been removed.
 integer, intent(in)        :: m
 integer, intent(in)        :: n
 real (dp), intent(in out)  :: a(:,:)
 logical, intent(in)        :: pivot
 integer, intent(out)       :: ipvt(:)
 real (dp), intent(out)     :: rdiag(:)
 real (dp), intent(out)     :: acnorm(:)

 !  **********

 !  subroutine qrfac

 !  This subroutine uses Householder transformations with column pivoting
 !  (optional) to compute a qr factorization of the m by n matrix a.
 !  That is, qrfac determines an orthogonal matrix q, a permutation matrix p,
 !  and an upper trapezoidal matrix r with diagonal elements of nonincreasing
 !  magnitude, such that a*p = q*r.  The householder transformation for
 !  column k, k = 1,2,...,min(m,n), is of the form

 !                        t
 !        i - (1/u(k))*u*u

 !  where u has zeros in the first k-1 positions.  The form of this
 !  transformation and the method of pivoting first appeared in the
 !  corresponding linpack subroutine.

 !  the subroutine statement is

 !    subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)

 ! N.B. 3 of these arguments have been omitted in this version.

 !  where

 !    m is a positive integer input variable set to the number of rows of a.

 !    n is a positive integer input variable set to the number of columns of a.

 !    a is an m by n array.  On input a contains the matrix for
 !      which the qr factorization is to be computed.  On output
 !      the strict upper trapezoidal part of a contains the strict
 !      upper trapezoidal part of r, and the lower trapezoidal
 !      part of a contains a factored form of q (the non-trivial
 !      elements of the u vectors described above).

 !    lda is a positive integer input variable not less than m
 !      which specifies the leading dimension of the array a.

 !    pivot is a logical input variable.  If pivot is set true,
 !      then column pivoting is enforced.  If pivot is set false,
 !      then no column pivoting is done.

 !    ipvt is an integer output array of length lipvt.  ipvt
 !      defines the permutation matrix p such that a*p = q*r.
 !      Column j of p is column ipvt(j) of the identity matrix.
 !      If pivot is false, ipvt is not referenced.

 !    lipvt is a positive integer input variable.  If pivot is false,
 !      then lipvt may be as small as 1.  If pivot is true, then
 !      lipvt must be at least n.

 !    rdiag is an output array of length n which contains the
 !      diagonal elements of r.

 !    acnorm is an output array of length n which contains the norms of the
 !      corresponding columns of the input matrix a.
 !      If this information is not needed, then acnorm can coincide with rdiag.

 !    wa is a work array of length n.  If pivot is false, then wa
 !      can coincide with rdiag.

 !  subprograms called

 !    minpack-supplied ... dpmpar,enorm

 !    fortran-supplied ... MAX,SQRT,MIN

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********
 integer   :: i, j, jp1, k, kmax, minmn
 real (dp) :: ajnorm, epsmch, sum, temp, wa(n)
 real (dp), parameter :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

 !     epsmch is the machine precision.

 epsmch = EPSILON(zero)

 !     compute the initial column norms and initialize several arrays.

 DO  j = 1, n
   acnorm(j) = enorm(m,a(1:,j))
   rdiag(j) = acnorm(j)
   wa(j) = rdiag(j)
   IF (pivot) ipvt(j) = j
 END DO

 !     Reduce a to r with Householder transformations.

 minmn = MIN(m,n)
 DO  j = 1, minmn
   IF (.NOT.pivot) GO TO 40

 !        Bring the column of largest norm into the pivot position.

   kmax = j
   DO  k = j, n
     IF (rdiag(k) > rdiag(kmax)) kmax = k
   END DO
   IF (kmax == j) GO TO 40
   DO  i = 1, m
     temp = a(i,j)
     a(i,j) = a(i,kmax)
     a(i,kmax) = temp
   END DO
   rdiag(kmax) = rdiag(j)
   wa(kmax) = wa(j)
   k = ipvt(j)
   ipvt(j) = ipvt(kmax)
   ipvt(kmax) = k

 !     Compute the Householder transformation to reduce the
 !     j-th column of a to a multiple of the j-th unit vector.

   40 ajnorm = enorm(m-j+1, a(j:,j))
   IF (ajnorm == zero) CYCLE
   IF (a(j,j) < zero) ajnorm = -ajnorm
   a(j:m,j) = a(j:m,j)/ajnorm
   a(j,j) = a(j,j) + one

 !     Apply the transformation to the remaining columns and update the norms.

   jp1 = j + 1
   DO  k = jp1, n
     sum = DOT_PRODUCT( a(j:m,j), a(j:m,k) )
     temp = sum/a(j,j)
     a(j:m,k) = a(j:m,k) - temp*a(j:m,j)
     IF (.NOT.pivot .OR. rdiag(k) == zero) CYCLE
     temp = a(j,k)/rdiag(k)
     rdiag(k) = rdiag(k)*SQRT(MAX(zero, one-temp**2))
     IF (p05*(rdiag(k)/wa(k))**2 > epsmch) CYCLE
     rdiag(k) = enorm(m-j, a(jp1:,k))
     wa(k) = rdiag(k)
   END DO
   rdiag(j) = -ajnorm
 END DO
 RETURN

 !     last card of subroutine qrfac.

 end subroutine qrfac
!!***

!----------------------------------------------------------------------

!!****f* m_levenberg_marquardt/qrsolv
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 SUBROUTINE qrsolv(n, r, ipvt, diag, qtb, x, sdiag)

 ! N.B. Arguments LDR & WA have been removed.
 integer, intent(in)        :: n
 real (dp), intent(in out)  :: r(:,:)
 integer, intent(in)        :: ipvt(:)
 real (dp), intent(in)      :: diag(:)
 real (dp), intent(in)      :: qtb(:)
 real (dp), intent(out)     :: x(:)
 real (dp), intent(out)     :: sdiag(:)

 !  **********

 !  subroutine qrsolv

 !  Given an m by n matrix a, an n by n diagonal matrix d, and an m-vector b,
 !  the problem is to determine an x which solves the system

 !        a*x = b ,     d*x = 0 ,

 !  in the least squares sense.

 !  This subroutine completes the solution of the problem if it is provided
 !  with the necessary information from the qr factorization, with column
 !  pivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
 !  q has orthogonal columns, and r is an upper triangular matrix with diagonal
 !  elements of nonincreasing magnitude, then qrsolv expects the full upper
 !  triangle of r, the permutation matrix p, and the first n components of
 !  (q transpose)*b.  The system a*x = b, d*x = 0, is then equivalent to

 !               t       t
 !        r*z = q *b ,  p *d*p*z = 0 ,

 !  where x = p*z. if this system does not have full rank,
 !  then a least squares solution is obtained.  On output qrsolv
 !  also provides an upper triangular matrix s such that

 !         t   t               t
 !        p *(a *a + d*d)*p = s *s .

 !  s is computed within qrsolv and may be of separate interest.

 !  the subroutine statement is

 !    subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)

 ! N.B. Arguments LDR and WA have been removed in this version.

 !  where

 !    n is a positive integer input variable set to the order of r.

 !    r is an n by n array.  On input the full upper triangle must contain
 !      the full upper triangle of the matrix r.
 !      On output the full upper triangle is unaltered, and the strict lower
 !      triangle contains the strict upper triangle (transposed) of the
 !      upper triangular matrix s.

 !    ldr is a positive integer input variable not less than n
 !      which specifies the leading dimension of the array r.

 !    ipvt is an integer input array of length n which defines the
 !      permutation matrix p such that a*p = q*r.  Column j of p
 !      is column ipvt(j) of the identity matrix.

 !    diag is an input array of length n which must contain the
 !      diagonal elements of the matrix d.

 !    qtb is an input array of length n which must contain the first
 !      n elements of the vector (q transpose)*b.

 !    x is an output array of length n which contains the least
 !      squares solution of the system a*x = b, d*x = 0.

 !    sdiag is an output array of length n which contains the
 !      diagonal elements of the upper triangular matrix s.

 !    wa is a work array of length n.

 !  subprograms called

 !    fortran-supplied ... ABS,SQRT

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********
 integer   :: i, j, k, kp1, l, nsing
 real (dp) :: cos, cotan, qtbpj, sin, sum, tan, temp, wa(n)
 real (dp), parameter :: p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

 !     Copy r and (q transpose)*b to preserve input and initialize s.
 !     In particular, save the diagonal elements of r in x.

 DO  j = 1, n
   r(j:n,j) = r(j,j:n)
   x(j) = r(j,j)
   wa(j) = qtb(j)
 END DO

 !     Eliminate the diagonal matrix d using a givens rotation.

 DO  j = 1, n

 !        Prepare the row of d to be eliminated, locating the
 !        diagonal element using p from the qr factorization.

   l = ipvt(j)
   IF (diag(l) == zero) CYCLE
   sdiag(j:n) = zero
   sdiag(j) = diag(l)

 !     The transformations to eliminate the row of d modify only a single
 !     element of (q transpose)*b beyond the first n, which is initially zero.

   qtbpj = zero
   DO  k = j, n

 !        Determine a givens rotation which eliminates the
 !        appropriate element in the current row of d.

     IF (sdiag(k) == zero) CYCLE
     IF (ABS(r(k,k)) < ABS(sdiag(k))) THEN
       cotan = r(k,k)/sdiag(k)
       SIN = p5/SQRT(p25 + p25*cotan**2)
       COS = SIN*cotan
     ELSE
       TAN = sdiag(k)/r(k,k)
       COS = p5/SQRT(p25 + p25*TAN**2)
       SIN = COS*TAN
     END IF

 !        Compute the modified diagonal element of r and
 !        the modified element of ((q transpose)*b,0).

     r(k,k) = COS*r(k,k) + SIN*sdiag(k)
     temp = COS*wa(k) + SIN*qtbpj
     qtbpj = -SIN*wa(k) + COS*qtbpj
     wa(k) = temp

 !        Accumulate the tranformation in the row of s.

     kp1 = k + 1
     DO  i = kp1, n
       temp = COS*r(i,k) + SIN*sdiag(i)
       sdiag(i) = -SIN*r(i,k) + COS*sdiag(i)
       r(i,k) = temp
     END DO
   END DO

 !     Store the diagonal element of s and restore
 !     the corresponding diagonal element of r.

   sdiag(j) = r(j,j)
   r(j,j) = x(j)
 END DO

 !     Solve the triangular system for z.  If the system is singular,
 !     then obtain a least squares solution.

 nsing = n
 DO  j = 1, n
   IF (sdiag(j) == zero .AND. nsing == n) nsing = j - 1
   IF (nsing < n) wa(j) = zero
 END DO

 DO  k = 1, nsing
   j = nsing - k + 1
   sum = DOT_PRODUCT( r(j+1:nsing,j), wa(j+1:nsing) )
   wa(j) = (wa(j) - sum)/sdiag(j)
 END DO

 !     Permute the components of z back to components of x.

 DO  j = 1, n
   l = ipvt(j)
   x(l) = wa(j)
 END DO
 RETURN

 !     last card of subroutine qrsolv.

 END SUBROUTINE qrsolv
!!***

!----------------------------------------------------------------------

!!****f* m_levenberg_marquardt/enorm
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function enorm(n,x) result(fn_val)

 ! Code converted using TO_F90 by Alan Miller
 ! Date: 1999-12-09  Time: 12:45:34

!Arguments ------------------------------------
!scalars
 real (dp)              :: fn_val
 integer, intent(in)    :: n
!arrays
 real (dp), intent(in)  :: x(n)

!Local variables-------------------------------
 integer   :: i
 real (dp) :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
 real (dp), parameter :: one = 1.0_dp, zero = 0.0_dp, rdwarf = 3.834e-20_dp, rgiant = 1.304e+19_dp
! *********************************************************************

 !  function enorm

 !  given an n-vector x, this function calculates the euclidean norm of x.

 !  the euclidean norm is computed by accumulating the sum of squares in
 !  three different sums.  The sums of squares for the small and large
 !  components are scaled so that no overflows occur.  Non-destructive
 !  underflows are permitted.  Underflows and overflows do not occur in the
 !  computation of the unscaled sum of squares for the intermediate
 !  components.  The definitions of small, intermediate and large components
 !  depend on two constants, rdwarf and rgiant.  The main restrictions on
 !  these constants are that rdwarf**2 not underflow and rgiant**2 not
 !  overflow.  The constants given here are suitable for every known computer.

 !  the function statement is

 !    REAL (dp) function enorm(n,x)

 !  where

 !    n is a positive integer input variable.

 !    x is an input array of length n.

 !  subprograms called

 !    fortran-supplied ... ABS,SQRT

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********


 s1 = zero
 s2 = zero
 s3 = zero
 x1max = zero
 x3max = zero
 floatn = n
 agiant = rgiant/floatn
 DO  i = 1, n
   xabs = ABS(x(i))
   IF (xabs > rdwarf .AND. xabs < agiant) GO TO 70
   IF (xabs <= rdwarf) GO TO 30

 !              sum for large components.

   IF (xabs <= x1max) GO TO 10
   s1 = one + s1*(x1max/xabs)**2
   x1max = xabs
   GO TO 20

   10 s1 = s1 + (xabs/x1max)**2

   20 GO TO 60

 !              sum for small components.

   30 IF (xabs <= x3max) GO TO 40
   s3 = one + s3*(x3max/xabs)**2
   x3max = xabs
   GO TO 60

   40 IF (xabs /= zero) s3 = s3 + (xabs/x3max)**2

   60 CYCLE

 !           sum for intermediate components.

   70 s2 = s2 + xabs**2
 END DO

 !     calculation of norm.

 IF (s1 == zero) GO TO 100
 fn_val = x1max*SQRT(s1 + (s2/x1max)/x1max)
 GO TO 120

 100 IF (s2 == zero) GO TO 110
 IF (s2 >= x3max) fn_val = SQRT(s2*(one + (x3max/s2)*(x3max*s3)))
 IF (s2 < x3max) fn_val = SQRT(x3max*((s2/x3max) + (x3max*s3)))
 GO TO 120

 110 fn_val = x3max*SQRT(s3)

 120 RETURN

 !     last card of function enorm.

 END FUNCTION enorm
!!***

!----------------------------------------------------------------------

!!****f* m_levenberg_marquardt/fdjac2
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 SUBROUTINE fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)

 ! Code converted using TO_F90 by Alan Miller
 ! Date: 1999-12-09  Time: 12:45:44

 ! N.B. Arguments LDFJAC & WA have been removed.
 integer, intent(in)        :: m
 integer, intent(in)        :: n
 real (dp), intent(in out)  :: x(n)
 real (dp), intent(in)      :: fvec(m)
 real (dp), intent(out)     :: fjac(:,:)    ! fjac(ldfjac,n)
 integer, intent(in out)    :: iflag
 real (dp), intent(in)      :: epsfcn

 interface
   subroutine fcn(m, n, x, fvec, iflag, y, nfqre, nfqim, bsign, csign, multi_x_exp)
     implicit none
     integer, parameter         :: dp = KIND(1.0d0)
     integer, intent(in)        :: m, n
     real (dp), intent(in)      :: x(n)
     real (dp), intent(in out)  :: fvec(m)
     integer, intent(in out)    :: iflag
     integer, optional, intent(in) :: nfqre,nfqim,bsign,csign,multi_x_exp
     real (dp), optional, intent(in)      :: y(m)
   end subroutine fcn
 end interface

 !  **********

 !  subroutine fdjac2

 !  this subroutine computes a forward-difference approximation
 !  to the m by n jacobian matrix associated with a specified
 !  problem of m functions in n variables.

 !  the subroutine statement is

 !    subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)

 !  where

 !    fcn is the name of the user-supplied subroutine which calculates the
 !      functions.  fcn must be declared in an external statement in the user
 !      calling program, and should be written as follows.

 !      subroutine fcn(m,n,x,fvec,iflag)
 !      integer m,n,iflag
 !      REAL (dp) x(n),fvec(m)
 !      ----------
 !      calculate the functions at x and
 !      return this vector in fvec.
 !      ----------
 !      return
 !      end

 !      the value of iflag should not be changed by fcn unless
 !      the user wants to terminate execution of fdjac2.
 !      in this case set iflag to a negative integer.

 !    m is a positive integer input variable set to the number of functions.

 !    n is a positive integer input variable set to the number of variables.
 !      n must not exceed m.

 !    x is an input array of length n.

 !    fvec is an input array of length m which must contain the
 !      functions evaluated at x.

 !    fjac is an output m by n array which contains the
 !      approximation to the jacobian matrix evaluated at x.

 !    ldfjac is a positive integer input variable not less than m
 !      which specifies the leading dimension of the array fjac.

 !    iflag is an integer variable which can be used to terminate
 !      the execution of fdjac2.  see description of fcn.

 !    epsfcn is an input variable used in determining a suitable step length
 !      for the forward-difference approximation.  This approximation assumes
 !      that the relative errors in the functions are of the order of epsfcn.
 !      If epsfcn is less than the machine precision, it is assumed that the
 !      relative errors in the functions are of the order of the machine
 !      precision.

 !    wa is a work array of length m.

 !  subprograms called

 !    user-supplied ...... fcn

 !    minpack-supplied ... dpmpar

 !    fortran-supplied ... ABS,MAX,SQRT

 !  argonne national laboratory. minpack project. march 1980.
 !  burton s. garbow, kenneth e. hillstrom, jorge j. more

 !  **********
 integer   :: j
 real (dp) :: eps, epsmch, h, temp, wa(m)
 real (dp), parameter :: zero = 0.0_dp

 !     epsmch is the machine precision.

 epsmch = EPSILON(zero)

 eps = SQRT(MAX(epsfcn, epsmch))
 DO  j = 1, n
   temp = x(j)
   h = eps*ABS(temp)
   IF (h == zero) h = eps
   x(j) = temp + h
   CALL fcn(m, n, x, wa, iflag)
   IF (iflag < 0) EXIT
   x(j) = temp
   fjac(1:m,j) = (wa(1:m) - fvec(1:m))/h
 END DO

 RETURN

!     last card of subroutine fdjac2.

 end subroutine fdjac2
!!***

!!****f* m_levenberg_marquardt/lm_fit_print_info
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine lm_fit_print_info(info,msg)

  integer, intent(in) :: info
  character(len=500)  :: msg

  select case (info)
  case (:-1)
     write(msg,'(a,i0)') 'STOP LM fit: fcn in LM fit returned info = ', -info
  case (0)
     write(msg,'(a)') 'ERROR LM fit: improper values for input parameters in LM fit'
  case (1:3)
     write(msg,'(a)') 'LM fit converged'
  case (4)
     write(msg,'(a)') 'ERROR LM fit: residuals orthogonal to the jacobian, there may be an error in fcn.'
  case (5)
     write(msg,'(a)') 'ERROR LM fit: too many calls to fcn, either slow convergence, or an error in opt.'
  case (6:7)
     write(msg,'(a)') 'ERROR LM fit: tol was set too small'
  case default
     write(msg,'(a,i0)') 'ERROR LM fit: unknown info = ',info
  end select

  return

end subroutine lm_fit_print_info

end module m_levenberg_marquardt
!!***

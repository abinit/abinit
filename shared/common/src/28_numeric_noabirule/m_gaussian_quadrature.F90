!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gaussian_quadrature
!! NAME
!! m_gaussian_quadrature
!!
!! FUNCTION
!!  .
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (JLJ, BR, MC)
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

module m_gaussian_quadrature
!----------------------------------------------------------------------------------------------------
!
! This module contains routines to generate the weights and positions at which to evaluate
! an integrand in order to produce a Gaussian quadrature of a given type.
!
! The code here was obtained from the internet (see header to subroutine cdgqf for further details).
! The code was slightly modified to use doubles instead of floats in order to reach greater
! accuracy.
!
!----------------------------------------------------------------------------------------------------

 use m_abicore
 use m_errors

 use defs_basis, only : dp, tol20, std_out
 use m_io_tools, only : open_file

 implicit none

 private
!!***

!!***
 public :: gaussian_quadrature_gegenbauer, gaussian_quadrature_legendre
 public :: get_frequencies_and_weights_legendre
 public :: cgqf
!!***

contains


!!****f* m_hamiltonian/gaussian_quadrature_gegenbauer
!! NAME
!!  gaussian_quadrature_gegenbauer
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

 subroutine gaussian_quadrature_gegenbauer(integrand,integral,order)
 !----------------------------------------------------------------------------------------------------
 !
 ! This subroutine performs a gaussian quadrature of given order of the input function integrand
 ! and stores the result in the variable integral.
 !
 ! The integrand function is assumed to be defined on the range from 0 to infinity, and it is
 ! assumed to to behave as x^2 near x=0, and 1/x^4 as x-> infinity.
 !
 ! The Gegenbauer integration method is of the form:
 !
 ! int_{-1}^{1} dt (1-t)^2(1+t)^2 f(t) =  sum_{i} w_i f(t_i)
 !
 !
 ! In order to perform such an integral, make a change of variable x = (1-t)/(1+t); this implies
 !
 !              int_0^infty dx I(x) = int_{-1}^{1} dt g(t),     g(t) = 2/(1+t)^2 I((1-t)/(1+t))
 !
 ! It is easy to show that g(1-delta) ~ delta^2 and g(-1+delta) ~ delta^2 for delta small.
 ! Thus, we can define
 !
 !              int_{-1}^{1} dt g(t) = int_{-1}^{1} dt (1+t)^2(1-t)^2 f(t), f(t) = g(t)/(1-t^2)^2
 !
 !
 !----------------------------------------------------------------------------------------------------
  implicit none


  ! input and output variables
  integer,   intent(in)  :: order

  interface
        function integrand(x)
        integer, parameter :: dp=kind(1.0d0)
        real(dp) integrand
        real(dp),intent(in) :: x
        end function integrand
  end interface

  real(dp),  intent(out) :: integral


  real(dp) :: t(order), w(order)


  integer  :: i

  integer  :: gaussian_kind
  real(dp) :: alpha, beta, a, b
  real(dp) :: x, y, f

! *************************************************************************


  ! Set the arguments and weights
  gaussian_kind = 3
  alpha=  2.0_dp
  beta =  2.0_dp
  a    = -1.0_dp
  b    =  1.0_dp
  call cgqf ( order, gaussian_kind, alpha, beta, a, b, t, w )



  ! call the integrant the right number of times

  integral = 0.0_dp

  do i = 1, order
        x = (1.0_dp-t(i))/(1.0_dp+t(i))

        y = integrand(x)

        f = 2.0_dp/(1.0_dp+t(i))**2/(1.0_dp-t(i)**2)**2*y

        integral = integral + w(i)*f

  end do ! i

 end subroutine gaussian_quadrature_gegenbauer
!!***

!!****f* m_hamiltonian/gaussian_quadrature_legendre
!! NAME
!!  gaussian_quadrature_legendre
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

 subroutine gaussian_quadrature_legendre(integrand,integral,order)
 !----------------------------------------------------------------------------------------------------
 !
 ! This subroutine performs a gaussian quadrature of given order of the input function integrand
 ! and stores the result in the variable integral.
 !
 ! The integrand function is assumed to be defined on the range from 0 to infinity, and it is
 ! assumed to to behave as x^2 near x=0, and 1/x^4 as x-> infinity.
 !
 ! The Legendre integration method is of the form:
 !
 ! int_{-1}^{1} dt f(t) =  sum_{i} w_i f(t_i)
 !
 !
 ! In order to perform such an integral, make a change of variable x = (1-t)/(1+t); this implies
 !
 !              int_0^infty dx I(x) = int_{-1}^{1} dt f(t),     f(t) = 2/(1+t)^2 I((1-t)/(1+t))
 !
 !
 !
 !----------------------------------------------------------------------------------------------------
  implicit none


  ! input and output variables
  integer,   intent(in)  :: order

  interface
        function integrand(x)
        integer, parameter :: dp=kind(1.0d0)
        real(dp) integrand
        real(dp),intent(in) :: x
        end function integrand
  end interface

  real(dp),  intent(out) :: integral


  real(dp) :: t(order), w(order)


  integer  :: i

  integer  :: gaussian_kind
  real(dp) :: alpha, beta, a, b
  real(dp) :: x, y, f

! *************************************************************************


  ! Set the arguments and weights
  gaussian_kind = 1
  alpha=  0.0_dp
  beta =  0.0_dp
  a    = -1.0_dp
  b    =  1.0_dp
  call cgqf ( order, gaussian_kind, alpha, beta, a, b, t, w )



  ! call the integrant the right number of times

  integral = 0.0_dp

  do i = 1, order
        x = (1.0_dp-t(i))/(1.0_dp+t(i))

        y = integrand(x)

        f = 2.0_dp/(1.0_dp+t(i))**2*y

        integral = integral + w(i)*f

  end do ! i

 end subroutine gaussian_quadrature_legendre
!!***


!!****f* m_hamiltonian/get_frequencies_and_weights_legendre
!! NAME
!!  get_frequencies_and_weights_legendre
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

 subroutine get_frequencies_and_weights_legendre(order,frequencies,weights)
 !----------------------------------------------------------------------------------------------------
 !
 ! This subroutine generates the frequencies and weights to integrate using legendre Gaussian
 ! quadrature for a given order.
 !
 ! The integrand function is assumed to be defined on the range from 0 to infinity, and it is
 ! assumed to to behave as x^2 near x=0, and 1/x^4 as x-> infinity.
 !
 ! The Legendre integration method is of the form:
 !
 ! int_{-1}^{1} dt f(t) =  sum_{i} w_i f(t_i)
 !
 ! In order to perform such an integral, make a change of variable
 !
 !              omega = (1-t)/(1+t)
 !
 ! and thus
 !       int_0^infty domega I(omega) = int_{-1}^{1} dt f(t),    f(t) = 2/(1+t)^2 I((1-t)/(1+t))
 !
 ! and so the weights are defined as weights(i) = 2/(1+ti)^2*wi.
 !
 !
 !----------------------------------------------------------------------------------------------------
  implicit none
  ! input and output variables
  integer,   intent(in)   :: order
  real(dp),  intent(out)  :: frequencies(order)
  real(dp),  intent(out)  :: weights(order)


  real(dp) :: t(order), w(order)


  integer  :: i

  integer  :: gaussian_kind
  real(dp) :: alpha, beta, a, b
  real(dp) :: x,  f

! *************************************************************************


  ! Set the arguments and weights
  gaussian_kind = 1
  alpha =  0.0_dp
  beta  =  0.0_dp
  a     = -1.0_dp
  b     =  1.0_dp

  call cgqf ( order, gaussian_kind, alpha, beta, a, b, t, w )


  ! call the integrant the right number of times

  do i = 1, order

        x = (1.0_dp-t(i))/(1.0_dp+t(i))
        f = 2.0_dp/(1.0_dp+t(i))**2

        frequencies(i) = x
        weights(i)     = w(i)*f

  end do ! i

 end subroutine get_frequencies_and_weights_legendre
!!***

!!****f* m_hamiltonian/cdgqf
!! NAME
!!  cdgqf
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine cdgqf ( nt, gaussian_kind, alpha, beta, t, wts )

!*****************************************************************************80
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with a classical weight function with default values for A and B,
!    and only simple knots.
!
!    There are no moments checks and no printing is done.
!
!    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, real ( dp ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( dp ) BETA, the value of Beta, if needed.
!
!    Output, real ( dp ) T(NT), the knots.
!
!    Output, real ( dp ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( dp ) aj(nt)
  real ( dp ) alpha
  real ( dp ) beta
  real ( dp ) bj(nt)
  integer ( kind = 4 ) gaussian_kind
  real ( dp ) t(nt)
  real ( dp ) wts(nt)
  real ( dp ) zemu

! *************************************************************************

  call parchk ( gaussian_kind, 2 * nt, alpha, beta )
!
!  Get the Jacobi matrix and zero-th moment.
!
  call class_matrix ( gaussian_kind, nt, alpha, beta, aj, bj, zemu )
!
!  Compute the knots and weights.
!
  call sgqf ( nt, aj, bj, zemu, t, wts )

  return
end subroutine cdgqf
!!***

!!****f* m_hamiltonian/cgqf
!! NAME
!!  cgqf
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_efmas,m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine cgqf ( nt, gaussian_kind, alpha, beta, a, b, t, wts )

!*****************************************************************************80
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!    The user may specify the interval (A,B).
!
!    Only simple knots are produced.
!
!    Use routine EIQFS to evaluate this quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( dp ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( dp ) BETA, the value of Beta, if needed.
!
!    Input, real ( dp ) A, B, the interval endpoints, or
!    other parameters.
!
!    Output, real ( dp ) T(NT), the knots.
!
!    Output, real ( dp ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( dp ) a
  real ( dp ) alpha
  real ( dp ) b
  real ( dp ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) gaussian_kind
  integer ( kind = 4 ), allocatable :: mlt(:)
  integer ( kind = 4 ), allocatable :: ndx(:)
  real ( dp ) t(nt)
  real ( dp ) wts(nt)

! *************************************************************************

!
!  Compute the Gauss quadrature formula for default values of A and B.
!
  call cdgqf ( nt, gaussian_kind, alpha, beta, t, wts )
!
!  Prepare to scale the quadrature formula to other weight function with
!  valid A and B.
!
  ABI_ALLOCATE ( mlt, (1:nt) )

  mlt(1:nt) = 1

  ABI_ALLOCATE ( ndx, (1:nt) )

  do i = 1, nt
    ndx(i) = i
  end do

  call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, gaussian_kind, alpha, beta, a, b )

  ABI_DEALLOCATE ( mlt )
  ABI_DEALLOCATE ( ndx )

  return
end subroutine cgqf
!!***

!!****f* m_hamiltonian/ch_cap
!! NAME
!!  ch_cap
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character              c
  integer ( kind = 4 ) itemp

! *************************************************************************

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end subroutine ch_cap
!!***

!!****f* m_hamiltonian/ch_eqi
!! NAME
!!  ch_eqi
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

! *************************************************************************

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end function ch_eqi
!!***

!!****f* m_hamiltonian/ch_to_digit
!! NAME
!!  ch_to_digit
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

! *************************************************************************

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end subroutine ch_to_digit
!!***

!!****f* m_hamiltonian/class_matrix
!! NAME
!!  class_matrix
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine class_matrix ( gaussian_kind, m, alpha, beta, aj, bj, zemu )

!*****************************************************************************80
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!    This routine computes the diagonal AJ and sub-diagonal BJ
!    elements of the order M tridiagonal symmetric Jacobi matrix
!    associated with the polynomials orthogonal with respect to
!    the weight function specified by KIND.
!
!    For weight functions 1-7, M elements are defined in BJ even
!    though only M-1 are needed.  For weight function 8, BJ(M) is
!    set to zero.
!
!    The zero-th moment of the weight function is returned in ZEMU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
!
!    Input, real ( dp ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( dp ) BETA, the value of Beta, if needed.
!
!    Output, real ( dp ) AJ(M), BJ(M), the diagonal and subdiagonal
!    of the Jacobi matrix.
!
!    Output, real ( dp ) ZEMU, the zero-th moment.
!
  implicit none

  integer ( kind = 4 ) m

  real ( dp ) a2b2
  real ( dp ) ab
  real ( dp ) aba
  real ( dp ) abi
  real ( dp ) abj
  real ( dp ) abti
  real ( dp ) aj(m)
  real ( dp ) alpha
  real ( dp ) apone
  real ( dp ) beta
  real ( dp ) bj(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) gaussian_kind
  real ( dp ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( dp ) temp
  real ( dp ) temp2
  real ( dp ) zemu

! *************************************************************************

  temp = epsilon ( temp )

  call parchk ( gaussian_kind, 2 * m - 1, alpha, beta )

  temp2 = 0.5D+00

  if ( 500.0D+00 * temp < abs ( ( r8_gamma_gq ( temp2 ) )**2 - pi ) ) then
    MSG_ERROR('Gamma function does not match machine parameters.')
  end if

  if ( gaussian_kind == 1 ) then

    ab = 0.0D+00

    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( gaussian_kind == 2 ) then

    zemu = pi

    aj(1:m) = 0.0D+00

    bj(1) =  sqrt ( 0.5D+00 )
    bj(2:m) = 0.5D+00

  else if ( gaussian_kind == 3 ) then

    ab = alpha * 2.0D+00
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma_gq ( alpha + 1.0D+00 )**2 &
      / r8_gamma_gq ( ab + 2.0D+00 )

    aj(1:m) = 0.0D+00
    bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    do i = 2, m
      bj(i) = i * ( i + ab ) / ( 4.0D+00 * ( i + alpha )**2 - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( gaussian_kind == 4 ) then

    ab = alpha + beta
    abi = 2.0D+00 + ab
    zemu = 2.0D+00**( ab + 1.0D+00 ) * r8_gamma_gq ( alpha + 1.0D+00 ) &
      * r8_gamma_gq ( beta + 1.0D+00 ) / r8_gamma_gq ( abi )
    aj(1) = ( beta - alpha ) / abi
    bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
      / ( ( abi + 1.0D+00 ) * abi * abi )
    a2b2 = beta * beta - alpha * alpha

    do i = 2, m
      abi = 2.0D+00 * i + ab
      aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
      abi = abi**2
      bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
        / ( ( abi - 1.0D+00 ) * abi )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( gaussian_kind == 5 ) then

    zemu = r8_gamma_gq ( alpha + 1.0D+00 )

    do i = 1, m
      aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
      bj(i) = i * ( i + alpha )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( gaussian_kind == 6 ) then

    zemu = r8_gamma_gq ( ( alpha + 1.0D+00 ) / 2.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( gaussian_kind == 7 ) then

    ab = alpha
    zemu = 2.0D+00 / ( ab + 1.0D+00 )

    aj(1:m) = 0.0D+00

    do i = 1, m
      abi = i + ab * mod ( i, 2 )
      abj = 2 * i + ab
      bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
    end do
    bj(1:m) =  sqrt ( bj(1:m) )

  else if ( gaussian_kind == 8 ) then

    ab = alpha + beta
    zemu = r8_gamma_gq ( alpha + 1.0D+00 ) * r8_gamma_gq ( - ( ab + 1.0D+00 ) ) &
      / r8_gamma_gq ( - beta )
    apone = alpha + 1.0D+00
    aba = ab * apone
    aj(1) = - apone / ( ab + 2.0D+00 )
    bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
    do i = 2, m
      abti = ab + 2.0D+00 * i
      aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
      aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
    end do

    do i = 2, m - 1
      abti = ab + 2.0D+00 * i
      bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
        / ( abti**2 ) * ( ab + i ) / ( abti + 1.0D+00 )
    end do

    bj(m) = 0.0D+00
    bj(1:m) =  sqrt ( bj(1:m) )

  end if

  return
end subroutine class_matrix
!!***

!!****f* m_hamiltonian/imtqlx
!! NAME
!!  imtqlx
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( dp ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( dp ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( dp ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( dp ) b
  real ( dp ) c
  real ( dp ) d(n)
  real ( dp ) e(n)
  real ( dp ) f
  real ( dp ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( dp ) p
  real ( dp ) prec
  real ( dp ) r
  real ( dp ) s
  real ( dp ) z(n)

! *************************************************************************

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( std_out, '(a)' ) ' '
        write ( std_out, '(a)' ) 'IMTQLX - Fatal error!'
        write ( std_out, '(a)' ) '  Iteration limit exceeded.'
        write ( std_out, '(a,i8)' ) '  J = ', j
        write ( std_out, '(a,i8)' ) '  L = ', l
        write ( std_out, '(a,i8)' ) '  M = ', m
        write ( std_out, '(a,i8)' ) '  N = ', n
        MSG_ERROR("Aborting now")
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end subroutine imtqlx
!!***

!!****f* m_hamiltonian/parchk
!! NAME
!!  parchk
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine parchk ( gaussian_kind, m, alpha, beta )

!*****************************************************************************80
!
!! PARCHK checks parameters ALPHA and BETA for classical weight functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
!
!    Input, integer ( kind = 4 ) M, the order of the highest moment to
!    be calculated.  This value is only needed when KIND = 8.
!
!    Input, real ( dp ) ALPHA, BETA, the parameters, if required
!    by the value of KIND.
!
  implicit none

  real ( dp ) alpha
  real ( dp ) beta
  integer ( kind = 4 ) gaussian_kind
  integer ( kind = 4 ) m
  real ( dp ) tmp

! *************************************************************************

  if ( gaussian_kind <= 0 ) then
    MSG_ERROR('PARCHK - Fatal error: KIND <= 0.')
  end if
!
!  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
!
  if ( 3 <= gaussian_kind .and. alpha <= -1.0D+00 ) then
    MSG_ERROR('PARCHK - Fatal error! 3 <= KIND and ALPHA <= -1.')
  end if
!
!  Check BETA for Jacobi.
!
  if ( gaussian_kind == 4 .and. beta <= -1.0D+00 ) then
    MSG_ERROR('PARCHK - Fatal error: KIND == 4 and BETA <= -1.0.')
  end if
!
!  Check ALPHA and BETA for rational.
!
  if ( gaussian_kind == 8 ) then
    tmp = alpha + beta + m + 1.0D+00
    if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
      MSG_ERROR('PARCHK - Fatal error!  KIND == 8 but condition on ALPHA and BETA fails.')
    end if
  end if

  return
end subroutine parchk
!!***


!!****f* m_hamiltonian/r8_gamma_gq
!! NAME
!!  r8_gamma_gq
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

function r8_gamma_gq ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( dp ) X, the argument of the function.
!
!    Output, real ( dp ) R8_GAMMA, the value of the function.
!
  implicit none

  real ( dp ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( dp ), parameter :: eps = 2.22D-16
  real ( dp ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( dp ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( dp ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( dp ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( dp ) r8_gamma_gq
  real ( dp ) res
  real ( dp ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( dp ) sum
  real ( dp ) x
  real ( dp ), parameter :: xbig = 171.624D+00
  real ( dp ) xden
  real ( dp ), parameter :: xinf = 1.0D+30
  real ( dp ), parameter :: xminin = 2.23D-308
  real ( dp ) xnum
  real ( dp ) y
  real ( dp ) y1
  real ( dp ) ysq
  real ( dp ) z

! *************************************************************************

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if (abs(res) > tol20) then !if ( res /= 0.0D+00 ) then

      if ( abs(y1 - aint ( y1 * 0.5D+00 ) * 2.0D+00) > tol20 ) then ! if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma_gq = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma_gq = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, dp )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma_gq = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( abs(fact - 1.0D+00) > tol20 ) then !if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma_gq = res

  return
end function r8_gamma_gq
!!***

!!****f* m_hamiltonian/r8mat_write
!! NAME
!!  r8mat_write
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( dp ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  character(len=500) :: msg
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( dp ) table(m,n)

! *************************************************************************

!
!  Open the file.
!
  if (open_file(output_filename, msg, newunit=output_unit, status='replace') /= 0) then
    MSG_ERROR(msg)
  end if

!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

end subroutine r8mat_write
!!***


!!****f* m_hamiltonian/rule_write
!! NAME
!!  rule_write
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine rule_write ( order, x, w, r, filename )

!*****************************************************************************80
!
!! RULE_WRITE writes a quadrature rule to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( dp ) X(ORDER), the abscissas.
!
!    Input, real ( dp ) W(ORDER), the weights.
!
!    Input, real ( dp ) R(2), defines the region.
!
!    Input, character ( len = * ) FILENAME, specifies the output.
!    'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining weights,
!    abscissas, and region.
!
  implicit none

  integer ( kind = 4 ) order

  character ( len = * ) filename
  character ( len = 255 ) filename_r
  character ( len = 255 ) filename_w
  character ( len = 255 ) filename_x
  real ( dp ) r(2)
  real ( dp ) w(order)
  real ( dp ) x(order)

! *************************************************************************

  filename_w = trim ( filename ) // '_w.txt'
  filename_x = trim ( filename ) // '_x.txt'
  filename_r = trim ( filename ) // '_r.txt'

  write ( std_out, '(a)' ) ' '
  write ( std_out, '(a)' ) '  Creating quadrature files.'
  write ( std_out, '(a)' ) ' '
  write ( std_out, '(a)' ) '  "Root" file name is   "' // trim ( filename ) // '".'
  write ( std_out, '(a)' ) ' '
  write ( std_out, '(a)' ) '  Weight file will be   "' // trim ( filename_w ) // '".'
  write ( std_out, '(a)' ) '  Abscissa file will be "' // trim ( filename_x ) // '".'
  write ( std_out, '(a)' ) '  Region file will be   "' // trim ( filename_r ) // '".'

  call r8mat_write ( filename_w, 1, order, w )
  call r8mat_write ( filename_x, 1, order, x )
  call r8mat_write ( filename_r, 1, 2,     r )

  return
end subroutine rule_write
!!***


!!****f* m_hamiltonian/s_to_i4
!! NAME
!!  s_to_i4
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  character ( len = * ) s

! *************************************************************************

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end subroutine s_to_i4
!!***

!!****f* m_hamiltonian/s_to_r8
!! NAME
!!  s_to_r8
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( dp ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  real ( dp ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( dp ) rbot
  real ( dp ) rexp
  real ( dp ) rtop
  character ( len = * ) s

! *************************************************************************

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, dp )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, dp )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( std_out, '(a)' ) ' '
    write ( std_out, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( std_out, '(a)' ) '  Illegal or nonnumeric input:'
    write ( std_out, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, dp ) &
        / real ( jbot, dp ) )
    end if
  end if

  dval = real ( isgn, dp ) * rexp * rtop / rbot

  return
end subroutine s_to_r8
!!***

!!****f* m_hamiltonian/scqf
!! NAME
!!  scqf
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, gaussian_kind, alpha, beta, a, &
  b )

!*****************************************************************************80
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!    The arrays WTS and SWTS may coincide.
!
!    The arrays T and ST may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( dp ) T(NT), the original knots.
!
!    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
!
!    Input, real ( dp ) WTS(NWTS), the weights.
!
!    Input, integer ( kind = 4 ) NWTS, the number of weights.
!
!    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
!    For more details see the comments in CAWIQ.
!
!    Output, real ( dp ) SWTS(NWTS), the scaled weights.
!
!    Output, real ( dp ) ST(NT), the scaled knots.
!
!    Input, integer ( kind = 4 ) KIND, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( dp ) ALPHA, the value of Alpha, if needed.
!
!    Input, real ( dp ) BETA, the value of Beta, if needed.
!
!    Input, real ( dp ) A, B, the interval endpoints.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nwts

  real ( dp ) a
  real ( dp ) al
  real ( dp ) alpha
  real ( dp ) b
  real ( dp ) be
  real ( dp ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) gaussian_kind
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mlt(nt)
  integer ( kind = 4 ) ndx(nt)
  real ( dp ) p
  real ( dp ) shft
  real ( dp ) slp
  real ( dp ) st(nt)
  real ( dp ) swts(nwts)
  real ( dp ) t(nt)
  real ( dp ) temp
  real ( dp ) tmp
  real ( dp ) wts(nwts)

! *************************************************************************

  temp = epsilon ( temp )

  call parchk ( gaussian_kind, 1, alpha, beta )

  if ( gaussian_kind == 1 ) then

    al = 0.0D+00
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      MSG_ERROR('SCQF - Fatal error! |B - A| too small.')
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( gaussian_kind == 2 ) then

    al = -0.5D+00
    be = -0.5D+00

    if ( abs ( b - a ) <= temp ) then
      MSG_ERROR('SCQF - Fatal error! |B - A| too small.')
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( gaussian_kind == 3 ) then

    al = alpha
    be = alpha

    if ( abs ( b - a ) <= temp ) then
      MSG_ERROR('SCQF - Fatal error!  |B - A| too small.')
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( gaussian_kind == 4 ) then

    al = alpha
    be = beta

    if ( abs ( b - a ) <= temp ) then
      MSG_ERROR("SCQF - Fatal error! |B - A| too small.")
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( gaussian_kind == 5 ) then

    if ( b <= 0.0D+00 ) then
      MSG_ERROR('SCQF - Fatal error!  B <= 0')
    end if

    shft = a
    slp = 1.0D+00 / b
    al = alpha
    be = 0.0D+00

  else if ( gaussian_kind == 6 ) then

    if ( b <= 0.0D+00 ) then
      MSG_ERROR('SCQF - Fatal error! B <= 0.')
    end if

    shft = a
    slp = 1.0D+00 / sqrt ( b )
    al = alpha
    be = 0.0D+00

  else if ( gaussian_kind == 7 ) then

    al = alpha
    be = 0.0D+00

    if ( abs ( b - a ) <= temp ) then
      MSG_ERROR('|B - A| too small.')
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  else if ( gaussian_kind == 8 ) then

    if ( a + b <= 0.0D+00 ) then
      MSG_ERROR('  A + B <= 0.')
    end if

    shft = a
    slp = a + b
    al = alpha
    be = beta

  else if ( gaussian_kind == 9 ) then

    al = 0.5D+00
    be = 0.5D+00

    if ( abs ( b - a ) <= temp ) then
      MSG_ERROR('|B - A| too small.')
    end if

    shft = ( a + b ) / 2.0D+00
    slp = ( b - a ) / 2.0D+00

  end if

  p = slp**( al + be + 1.0D+00 )

  do k = 1, nt

    st(k) = shft + slp * t(k)
    l = abs ( ndx(k) )

    if ( l /= 0 ) then
      tmp = p
      do i = l, l + mlt(k) - 1
        swts(i) = wts(i) * tmp
        tmp = tmp * slp
      end do
    end if

  end do

  return
end subroutine scqf
!!***


!!****f* m_hamiltonian/sgqf
!! NAME
!!  sgqf
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gaussian_quadrature
!!
!! CHILDREN
!!      imtqlx
!!
!! SOURCE

subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*****************************************************************************80
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!    This routine computes all the knots and weights of a Gauss quadrature
!    formula with simple knots from the Jacobi matrix and the zero-th
!    moment of the weight function, using the Golub-Welsch technique.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of knots.
!
!    Input, real ( dp ) AJ(NT), the diagonal of the Jacobi matrix.
!
!    Input/output, real ( dp ) BJ(NT), the subdiagonal of the Jacobi
!    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
!
!    Input, real ( dp ) ZEMU, the zero-th moment of the weight function.
!
!    Output, real ( dp ) T(NT), the knots.
!
!    Output, real ( dp ) WTS(NT), the weights.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( dp ) aj(nt)
  real ( dp ) bj(nt)
  real ( dp ) t(nt)
  real ( dp ) wts(nt)
  real ( dp ) zemu

! *************************************************************************

!
!  Exit if the zero-th moment is not positive.
!
  if ( zemu <= 0.0D+00 ) then
    MSG_ERROR('SGQF - Fatal error! ZEMU <= 0.')
  end if
!
!  Set up vectors for IMTQLX.
!
  t(1:nt) = aj(1:nt)

  wts(1) = sqrt ( zemu )
  wts(2:nt) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt)**2

  return
end subroutine sgqf
!!***

end module m_gaussian_quadrature
!!***

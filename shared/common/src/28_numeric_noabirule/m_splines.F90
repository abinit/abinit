!!****m* ABINIT/m_splines
!! NAME
!!  m_splines
!!
!! FUNCTION
!!  This module contains routines for spline interpolation.
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (YP, BAmadon)
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

module m_splines

 use defs_basis
 use m_abicore
 use m_errors

 use m_fstrings, only : sjoin, itoa, ftoa
 !use m_time,   only : timab

 implicit none

 public :: splfit
 public :: spline
 public :: spline_bicubic
 public :: spline_c
 public :: spline_complex
 public :: spline_integrate
 public :: splint
 public :: splint_complex

 !FIXME deprecated
 public :: intrpl

! *************************************************************************

contains
!!***

!----------------------------------------------------------------------

!!****f* m_splines/splfit
!! NAME
!!  splfit
!!
!! FUNCTION
!!  Evaluate cubic spline fit to get function values on input set of ORDERED, UNFORMLY SPACED points.
!!  Optionally gives derivatives (first and second) at those points too.
!!  If point lies outside the range of arg, assign the extremal
!!  point values to these points, and zero derivative.
!!
!! INPUTS
!!  arg(numarg)=equally spaced arguments (spacing de) for data to which spline was fit.
!!  fun(numarg,2)=function values to which spline was fit and spline
!!   fit to second derivatives (from Numerical Recipes spline).
!!  ider=  see above
!!  newarg(numnew)=new values of arguments at which function is desired.
!!  numarg=number of arguments at which spline was fit.
!!  numnew=number of arguments at which function values are desired.
!!
!! OUTPUT
!!  derfun(numnew)=(optional) values of first or second derivative of function.
!!   This is only computed for ider=1 or 2; otherwise derfun not used.
!!  newfun(numnew)=values of function at newarg(numnew).
!!   This is only computed for ider=0 or 1.
!!
!! NOTES
!!        if ider=0, compute only the function (contained in fun)
!!        if ider=1, compute the function (contained in fun) and its first derivative (in derfun)
!!        if ider=2, compute only the second derivative of the function (in derfun)
!!
!! PARENTS
!!      m_mkffnl,m_mklocl,m_occ,m_pawpwij,m_psptk
!!
!! CHILDREN
!!
!! SOURCE

subroutine splfit(arg, derfun, fun, ider, newarg, newfun, numarg, numnew)

 integer, intent(in) :: ider, numarg, numnew
 real(dp), intent(in) :: arg(numarg), fun(numarg,2), newarg(numnew)
 real(dp), intent(out) :: derfun(numnew)
 real(dp), intent(inout) :: newfun(numnew)

!Local variables---------------------------------------
 integer :: i,jspl
 real(dp) :: argmin,de,d,aa,bb,cc,dd,de2_dby_six,de_dby_six
 !real(dp) :: tsec(2)

! *************************************************************************

 ! Keep track of time spent in mkffnl
 !call timab(1905, 1, tsec)

 ! argmin is smallest x value in spline fit; de is uniform spacing of spline argument
 argmin = arg(1)
 de = (arg(numarg) - argmin) / dble(numarg-1)
 de2_dby_six = de**2 / six
 de_dby_six = de / six

 if (de < tol12) then
   ABI_ERROR(sjoin('spacing should be strictly positive, while de is: ', ftoa(de)))
 endif

 jspl = -1

 ! Do one loop for no grads, other for grads
 select case (ider)
 case (0)

  ! Spline index loop for no grads:
  do i=1,numnew
    if (newarg(i) >= arg(numarg)) then
      ! function values are being requested outside range of data.',a1,'
      ! Function and slope will be set to values at upper end of data.

      newfun(i) = fun(numarg,1)

    else if (newarg(i) <= arg(1)) then
      newfun(i) = fun(1,1)

    else
      jspl = 1 + int((newarg(i) - argmin)/de)
      d = newarg(i) - arg(jspl)
      bb = d / de
      aa = one - bb
      cc = aa*(aa**2 -one) * de2_dby_six
      dd = bb*(bb**2 -one) * de2_dby_six
      newfun(i)= aa * fun(jspl,1) + bb*fun(jspl+1,1) + cc*fun(jspl,2) + dd*fun(jspl+1,2)
    end if
  enddo

 case (1)

   ! Spline index loop includes grads:
   do i=1,numnew

     if (newarg(i) >= arg(numarg)) then
       newfun(i) = fun(numarg,1)
       derfun(i) = zero

     else if (newarg(i) <= arg(1)) then
       newfun(i) = fun(1,1)
       derfun(i) = zero

     else
       ! cubic spline interpolation:
       jspl = 1 + int((newarg(i) - arg(1)) / de)
       d = newarg(i) - arg(jspl)
       bb = d / de
       aa = one - bb
       cc = aa*(aa**2 - one) * de2_dby_six
       dd = bb*(bb**2 - one) * de2_dby_six
       newfun(i) = aa*fun(jspl,1) + bb*fun(jspl+1,1) + cc*fun(jspl,2) + dd*fun(jspl+1,2)
       ! spline fit to first derivative:
       ! note correction of Numerical Recipes sign error
       derfun(i) = (fun(jspl+1,1)-fun(jspl,1)) / de +    &
          (-(3.d0*aa**2 -one) * fun(jspl,2) + (3.d0*bb**2 -one) * fun(jspl+1,2)) * de_dby_six

     end if
   enddo

 case (2)

   do i=1,numnew

     if (newarg(i) >= arg(numarg)) then
       derfun(i) = zero

     else if (newarg(i) <= arg(1)) then
       derfun(i) = zero

     else
       ! cubic spline interpolation:
       jspl = 1 + int((newarg(i) - argmin) / de)
       d = newarg(i) - arg(jspl)
       bb = d / de
       aa = one - bb
       ! second derivative of spline (piecewise linear function)
       derfun(i) = aa*fun(jspl,2) + bb*fun(jspl+1,2)

     end if
   enddo

 case default
   ABI_ERROR(sjoin("Invalid ider:", itoa(ider)))
 end select

 !call timab(1905, 2, tsec)

end subroutine splfit
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline
!! NAME
!!  spline
!!
!! FUNCTION
!!  SPLINE (originally SPLINE_CUBIC_SET) computes the second derivatives
!!  of a cubic spline.
!!
!! INPUTS
!!    Input, integer N, the number of data points; N must be at least 2.
!!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the
!!    spline will actually be linear.
!!
!!    Input, double precision T(N), the knot values, that is, the points where data
!!    is specified.  The knot values should be distinct, and increasing.
!!
!!    Input, double precision Y(N), the data values to be interpolated.
!!
!!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!!
!! OUTPUT
!!    Output, double precision YPP(N), the second derivatives of the cubic spline.
!!    Work space, double precision DIAG(N) - should be removed ...
!!
!! PARENTS
!!      m_a2ftr,m_bader,m_dens,m_entropyDMFT,m_mkrho,m_occ,m_outscfcv
!!      m_paw_atomorb,m_paw_init,m_paw_mkrho,m_paw_slater,m_predict_string
!!      m_psp1,m_psp5,m_psp6,m_psp8,m_psp9,m_psp_hgh,m_psptk,m_screening_driver
!!      m_sigc,m_special_funcs,m_spin_current,m_splines,m_upf2abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine spline( t, y, n, ybcbeg, ybcend, ypp )

!*******************************************************************************
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    07 February 1999
!    28 November 2004 XGonze : double precision
!                              make arguments similar to the Numeric Recipes routine
!                              also use algorithmics similar to the Numeric Recipes routine
!
!  Author:
!
!    John Burkardt
!    (XGonze got it from http://www.psc.edu/~burkardt/src/spline/spline.html)
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2.
!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the
!    spline will actually be linear.
!
!    Input, double precision T(N), the knot values, that is, the points where data
!    is specified.  The knot values should be distinct, and increasing.
!
!    Input, double precision Y(N), the data values to be interpolated.
!
!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!
!    Output, double precision YPP(N), the second derivatives of the cubic spline.
!
!    Work space, double precision DIAG(N) - should be removed ...
!
!
!    XG041127 : In the initial implementation, one had the control on
!     IBCBEG and IBCEND. Now, they are determined by the values
!     of YBCBEG, YBCEND. Option 2 has been disabled.
!
!    Input, integer IBCBEG, left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, integer IBCEND, right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.

  integer, intent(in) :: n
  real(dp), intent(in) :: t(n)
  real(dp), intent(in) :: y(n)
  real(dp), intent(in) :: ybcbeg
  real(dp), intent(in) :: ybcend

  real(dp), intent(out) :: ypp(n)

  integer :: ibcbeg
  integer :: ibcend
  integer :: i,k
  real(dp) :: ratio,pinv
  real(dp), allocatable :: tmp(:)
!
!  Check.
!
  if ( n <= 1 ) then
    write(std_out,* ) ' '
    write(std_out,* ) 'SPLINE_CUBIC_SET - Fatal error!'
    write(std_out,* ) '  The number of knots must be at least 2.'
    write(std_out,* ) '  The input value of N = ', n
    ABI_ERROR("Fatal error")
  end if

  ABI_MALLOC(tmp,(n))

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write(std_out,* ) ' '
      write(std_out,* ) 'SPLINE_CUBIC_SET - Fatal error!'
      write(std_out,* ) '  The knots must be strictly increasing, but'
      write(std_out,* ) '  T(',  i,') = ', t(i)
      write(std_out,* ) '  T(',i+1,') = ', t(i+1)
      ABI_ERROR("Fatal error")
    end if
  end do
!
!  XG041127
  ibcbeg=1 ; ibcend=1
  if(ybcbeg>1.0d+30)ibcbeg=0
  if(ybcend>1.0d+30)ibcend=0
!
!  Set the first and last equations.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.d0
    tmp(1) = 0.d0
  else if ( ibcbeg == 1 ) then
    ypp(1) = -0.5d0
    tmp(1) = (3.d0/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1))-ybcbeg)
  end if
  if ( ibcend == 0 ) then
    ypp(n) = 0.d0
    tmp(n) = 0.d0
  else if ( ibcend == 1 ) then
    ypp(n) = 0.5d0
    tmp(n) = (3.d0/(t(n)-t(n-1)))*(ybcend-(y(n)-y(n-1))/(t(n)-t(n-1)))
  end if

!
!  Set the intermediate equations.
!
  do i=2,n-1
   ratio=(t(i)-t(i-1))/(t(i+1)-t(i-1))
   pinv = 1.0d0/(ratio*ypp(i-1) + 2.0d0)
   ypp(i) = (ratio-1.0d0)*pinv
   tmp(i)=(6.0d0*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1)) &
&    /(t(i)-t(i-1)))/(t(i+1)-t(i-1))-ratio*tmp(i-1))*pinv
   if (abs(tmp(i))<1.d5*tiny(0.d0)) tmp(i)=0.d0   !MT20050927
  enddo

! Solve the equations
  ypp(n) = (tmp(n)-ypp(n)*tmp(n-1))/(ypp(n)*ypp(n-1)+1.0d0)
  do k=n-1,1,-1
   ypp(k)=ypp(k)*ypp(k+1)+tmp(k)
  enddo

  ABI_FREE(tmp)
end subroutine spline
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline_bicubic
!! NAME
!!  spline_bicubic
!!
!! FUNCTION
!!  Generates coefficients for bicubic spline interpolation.
!!
!! INPUTS
!!  n1 = length of first dimension
!!  n2 = length of second dimension
!!  x1 = positions on first dimension
!!  x2 = positions on second dimension
!!  y = function values on the (x1,x2) grid
!!  der1_x1 = first derivative of y wrt x1
!!  der1_x2 = first derivative of y wrt x2
!!  der2_x1x2 = second-order cross-derivative of y wrt x1x2
!!
!! OUTPUT
!!  spl_c = spline coefficients
!!
!! NOTES
!!  Adapted from Numerical Recipes and libbci.
!!
!! PARENTS
!!      m_xc_vdw
!!
!! CHILDREN
!!
!! SOURCE

subroutine spline_bicubic(n1,n2,x1,x2,y,der1_x1,der1_x2,der2_x1x2,spl_c)

  integer,intent(in)  :: n1,n2
  real(dp),intent(in) :: x1(n1),x2(n2),y(n1,n2)
  real(dp),intent(in) :: der1_x1(n1,n2),der1_x2(n1,n2),der2_x1x2(n1,n2)
  real(dp),intent(out):: spl_c(4,4,n1,n2)

  integer :: i1,i2
  real(dp) :: dx1,dx2,wt(16,16),z(16)

  data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4, &
&          8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6, &
&          2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3, &
&          2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0, &
&          -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0, &
&          -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2, &
&          -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

  ! Set coefficients for i1<n1 and i2<n2
  do i2 = 1,n2-1
    do i1 = 1,n1-1
      dx1 = x1(i1+1) - x1(i1)
      dx2 = x2(i2+1) - x2(i2)
      z(1)  = y(i1,i2)
      z(2)  = y(i1+1,i2)
      z(3)  = y(i1+1,i2+1)
      z(4)  = y(i1,i2+1)
      z(5)  = der1_x1(i1,i2) * dx1
      z(6)  = der1_x1(i1+1,i2) * dx1
      z(7)  = der1_x1(i1+1,i2+1) * dx1
      z(8)  = der1_x1(i1,i2+1) * dx1
      z(9)  = der1_x2(i1,i2) * dx2
      z(10) = der1_x2(i1+1,i2) * dx2
      z(11) = der1_x2(i1+1,i2+1) * dx2
      z(12) = der1_x2(i1,i2+1) * dx2
      z(13) = der2_x1x2(i1,i2) * dx1 * dx2
      z(14) = der2_x1x2(i1+1,i2) * dx1 * dx2
      z(15) = der2_x1x2(i1+1,i2+1) * dx1 * dx2
      z(16) = der2_x1x2(i1,i2+1) * dx1 * dx2
      z = matmul(wt,z)
      spl_c(:,:,i1,i2) = reshape(z,(/4,4/),order=(/2,1/))
    end do
  end do

! Set coefficients for i1=n1 and i2=n2 (valid only at the border)
  spl_c(:,:,n1,:) = 0
  spl_c(:,:,:,n2) = 0
  spl_c(1,1,n1,:) = y(n1,:)
  spl_c(1,1,:,n2) = y(:,n2)

end subroutine spline_bicubic
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline_c
!! NAME
!!  spline_c
!!
!! FUNCTION
!!  Computes the spline of a complex function.
!!
!! INPUTS
!!  nomega_lo   = number of point in the non regular grid (e.g.  !logarithmic)
!!  nomega_li   = number of point in the regular grid on which the  spline is computed
!!  omega_lo    = value of freq on the 1st grid
!!  omega_li    = value of freq on the 2nd grid
!!  tospline_lo = function on the 1st grid
!!
!! OUTPUT
!!  splined_lo  = spline  (on the 2nd grid)
!!
!! PARENTS
!!      m_green
!!
!! CHILDREN
!!
!! SOURCE

subroutine spline_c( nomega_lo, nomega_li, omega_lo, omega_li, splined_li, tospline_lo)

!Arguments --------------------------------------------
!scalars
 integer, intent(in) :: nomega_lo, nomega_li
 real(dp), intent(in) :: omega_lo(nomega_lo)
 real(dp), intent(in) :: omega_li(nomega_li)
 complex(dpc), intent(in) :: tospline_lo(nomega_lo)
 complex(dpc), intent(out) :: splined_li(nomega_li)

!Local variables---------------------------------------
!scalars
 complex(dpc) :: ybcbeg, ybcend
 complex(dpc), allocatable :: ysplin2_lo(:)

 ybcbeg=czero
 ybcend=czero

 ABI_MALLOC(ysplin2_lo,(nomega_lo))
 call spline_complex(omega_lo, tospline_lo, nomega_lo, ybcbeg, ybcend, ysplin2_lo)
 call splint_complex( nomega_lo, omega_lo, tospline_lo,ysplin2_lo, nomega_li, omega_li, splined_li)
 ABI_FREE(ysplin2_lo)

end subroutine spline_c
!!***

!----------------------------------------------------------------------

!!****f* m_splines/spline_complex
!! NAME
!!  spline_complex
!!
!! FUNCTION
!!  spline_complex interfaces the usual spline routine in a case of a
!!  complex function
!!
!! INPUTS
!!    Input, integer N, the number of data points; N must be at least 2.
!!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the
!!    spline will actually be linear.
!!
!!    Input, double precision T(N), the knot values, that is, the points where data
!!    is specified.  The knot values should be distinct, and increasing.
!!
!!    Input, complex Y(N), the data values to be interpolated.
!!
!!    Input, complex YBCBEG, YBCEND, the values to be used in the boundary
!!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!!
!! OUTPUT
!!    Output, complex YPP(N), the second derivatives of the cubic spline.
!!
!! PARENTS
!!      m_paw_dmft,m_splines
!!
!! CHILDREN
!!
!! SOURCE

subroutine spline_complex( t, y, n, ybcbeg, ybcend, ypp )

 integer, intent(in) :: n
 real(dp), intent(in) :: t(n)
 complex(dpc), intent(in) :: y(n)
 complex(dpc), intent(in) :: ybcbeg
 complex(dpc), intent(in) :: ybcend
 complex(dpc), intent(out) :: ypp(n)

 real(dp), allocatable :: y_r(:)
 real(dp) :: ybcbeg_r
 real(dp) :: ybcend_r
 real(dp), allocatable :: ypp_r(:)
 real(dp), allocatable :: y_i(:)
 real(dp) :: ybcbeg_i
 real(dp) :: ybcend_i
 real(dp), allocatable :: ypp_i(:)

 ABI_MALLOC(y_r,(n))
 ABI_MALLOC(ypp_r,(n))
 ABI_MALLOC(y_i,(n))
 ABI_MALLOC(ypp_i,(n))
 y_r=real(y)
 y_i=aimag(y)    !vz_d
 ybcbeg_r=real(ybcbeg)
 ybcbeg_i=aimag(ybcbeg)    !vz_d
 ybcend_r=real(ybcend)
 ybcend_i=aimag(ybcend)    !vz_d
 call spline( t, y_r, n, ybcbeg_r, ybcend_r, ypp_r )
 call spline( t, y_i, n, ybcbeg_i, ybcend_i, ypp_i )
 ypp=cmplx(ypp_r,ypp_i)
 ABI_FREE(y_r)
 ABI_FREE(ypp_r)
 ABI_FREE(y_i)
 ABI_FREE(ypp_i)

end subroutine spline_complex
!!***

!----------------------------------------------------------------------

!!****f* m_splines/splint
!! NAME
!!  splint
!!
!! FUNCTION
!!  Compute spline interpolation. There is no hypothesis
!!  about the spacing of the input grid points.
!!
!! INPUTS
!!  nspline: number of grid points of input mesh
!!  xspline(nspline): input mesh
!!  yspline(nspline): function on input mesh
!!  ysplin2(nspline): second derivative of yspline on input mesh
!!  nfit: number of points of output mesh
!!  xfit(nfit): output mesh
!!
!! OUTPUT
!!  yfit(nfit): function on output mesh
!!  [ierr]=A non-zero value is used to signal that some points in xfit exceed xspline(nspline).
!!    The input value is incremented by the number of such points.
!!
!! PARENTS
!!      m_a2ftr,m_cut3d,m_entropyDMFT,m_epjdos,m_mkrho,m_outscfcv,m_paw_atomorb
!!      m_paw_mkrho,m_paw_slater,m_predict_string,m_psp6,m_psptk
!!      m_screening_driver,m_sigc,m_special_funcs,m_spin_current,m_splines
!!      m_wvl_rho
!!
!! CHILDREN
!!
!! SOURCE

subroutine splint(nspline,xspline,yspline,ysplin2,nfit,xfit,yfit,ierr)

 integer, intent(in) :: nfit, nspline
 integer,optional,intent(out) :: ierr
 real(dp), intent(in) :: xspline(nspline)
 real(dp), intent(in) :: yspline(nspline)
 real(dp), intent(in) :: ysplin2(nspline)
 real(dp), intent(in) :: xfit(nfit)
 real(dp), intent(out) :: yfit(nfit)

!local
 integer :: left,i,k,right,my_err
 real(dp) :: delarg,invdelarg,aa,bb

!source

 my_err=0

 left = 1
 do i=1, nfit
   yfit(i)=0.d0  ! Initialize for the unlikely event that rmax exceed r(mesh)
   !
   do k=left+1, nspline
     if(xspline(k) >= xfit(i)) then
       if(xspline(k-1) <= xfit(i)) then
         right = k
         left = k-1
       else
         if (k-1.eq.1 .and. i.eq.1) then
           ABI_ERROR('xfit(1) < xspline(1)')
           !my_err=my_err+1
           !exit
         else
           ABI_ERROR('xfit not properly ordered')
         end if
       end if
       delarg= xspline(right) - xspline(left)
       invdelarg= 1.0d0/delarg
       aa= (xspline(right)-xfit(i))*invdelarg
       bb= (xfit(i)-xspline(left))*invdelarg

       yfit(i) = aa*yspline(left) + bb*yspline(right)    &
&               +( (aa*aa*aa-aa)*ysplin2(left) +         &
&                  (bb*bb*bb-bb)*ysplin2(right) ) *delarg*delarg/6.0d0
       exit
     end if
   end do ! k
   !
   if (k==nspline+1) my_err=my_err+1 ! xfit not found
 end do ! i

 if (PRESENT(ierr)) ierr=my_err

end subroutine splint
!!***

!----------------------------------------------------------------------

!!****f* m_splines/splint_complex
!! NAME
!!  splint_complex
!!
!! FUNCTION
!!  Interface to the usual splint to compute *complex* spline interpolation. There is no hypothesis
!!  about the spacing of the input grid points.
!!
!! INPUTS
!!  nspline: number of grid points of input mesh
!!  xspline(nspline): input mesh
!!  yspline(nspline): complex function on input mesh
!!  ysplin2(nspline): second derivative of yspline on input mesh
!!  nfit: number of points of output mesh
!!  xfit(nfit): output mesh
!!
!! OUTPUT
!!  yfit(nfit): complex function on output mesh
!!
!! PARENTS
!!      m_paw_dmft,m_splines
!!
!! CHILDREN
!!
!! SOURCE

subroutine splint_complex (nspline,xspline,yspline,ysplin2,nfit,xfit,yfit)

 integer, intent(in) :: nfit, nspline
 real(dp), intent(in) :: xspline(nspline)
 complex(dpc), intent(in) :: yspline(nspline)
 complex(dpc), intent(in) :: ysplin2(nspline)
 real(dp), intent(in) :: xfit(nfit)
 complex(dpc), intent(out) :: yfit(nfit)

 real(dp), allocatable :: ysplin2_r(:)
 real(dp), allocatable :: ysplin2_i(:)
 real(dp), allocatable :: yspline_r(:)
 real(dp), allocatable :: yspline_i(:)
 real(dp), allocatable :: yfit_r(:)
 real(dp), allocatable :: yfit_i(:)

 ABI_MALLOC(yspline_r,(nspline))
 ABI_MALLOC(yspline_i,(nspline))
 ABI_MALLOC(ysplin2_r,(nspline))
 ABI_MALLOC(ysplin2_i,(nspline))
 ABI_MALLOC(yfit_r,(nfit))
 ABI_MALLOC(yfit_i,(nfit))

!local

!source
 yspline_r=real(yspline)
 yspline_i=aimag(yspline)    !vz_d
 ysplin2_r=real(ysplin2)
 ysplin2_i=aimag(ysplin2)    !vz_d
 call splint (nspline,xspline,yspline_r,ysplin2_r,nfit,xfit,yfit_r)
 call splint (nspline,xspline,yspline_i,ysplin2_i,nfit,xfit,yfit_i)
 yfit=cmplx(yfit_r,yfit_i)
 ABI_FREE(yspline_r)
 ABI_FREE(yspline_i)
 ABI_FREE(ysplin2_r)
 ABI_FREE(ysplin2_i)
 ABI_FREE(yfit_r)
 ABI_FREE(yfit_i)

end subroutine splint_complex
!!***

!!****f* m_splines/spline_integrate
!! NAME
!!  spline_integrate
!!
!! FUNCTION
!!  Calculates an integral using cubic spline interpolation.
!!
!! INPUTS
!!  npts= number of grid points of input mesh
!!  dx= step of input mesh
!!  integrand= function on input mesh
!!
!! OUTPUT
!!  integral= integral of the input function
!!
!! PARENTS
!!      test_spline_integrate
!!
!! CHILDREN
!!
!! SOURCE

subroutine spline_integrate(integral,npts,dx,integrand)

 integer,intent(in) :: npts
 real(dp),intent(out) :: integral
 real(dp),intent(in) :: dx,integrand(npts)

 integer :: ix
 real(dp) :: ptmp,sf(npts),sf_der2(npts),sf_mesh(npts),utmp(npts)

 ! Prepare mesh
 forall (ix=1:npts) sf_mesh(ix) = (ix - 1) * dx

 ! Calculate second derivative of integrand (adapted from Numercial Recipes)
 sf_der2(1) = zero
 sf_der2(npts) = zero
 utmp(1) = zero

 do ix=2,npts-1
  ptmp = half * sf_der2(ix-1) + two
  sf_der2(ix) = (half - one) / ptmp
  utmp(ix) = (three * (integrand(ix+1) + integrand(ix-1) - &
&  two*integrand(ix)) / (dx**2) - half * utmp(ix-1)) / ptmp
 end do
 do ix=npts-1,1,-1
  sf_der2(ix) = sf_der2(ix) * sf_der2(ix+1) + utmp(ix)
 end do

 ! Actually calculate integral
 sf(:) = integrand(:) * dx
 integral = (sf(1) + sf(npts)) / 2.0_dp - &
&           (sf_der2(1) + sf_der2(npts)) / 24.0_dp + &
&           sum(sf(2:npts-1)) - sum(sf_der2(2:npts-1)) / 12.0_dp

end subroutine spline_integrate
!!***

!!****f* m_splines/intrpl
!! NAME
!!  intrpl
!!
!! FUNCTION
!!
!!  DOUBLE PRECISION INTERPOLATION OF A SINGLE VALUED FUNCTION.
!!  THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
!!  GIVEN  AS ORDINATES OF INPUT DATA POINTS IN AN X-Y PLANE
!!  AND FOR A GIVEN SET OF X VALUES(ABSCISSAE),THE VALUES OF
!!  A SINGLE VALUED FUNCTION Y=Y(X).
!!
!!  THE SUBROUTINE ALSO CALCULATES FIRST DERIVATIVES DV(X) AND
!!  SECOND DERIVATIVE DV2(X)
!
!!  THE INPUT PARAMETERS ARE;
!!
!!  L=NUMBER OF DATA POINTS
!!  (MUST BE TWO OR GREATER)
!!  X=ARRAY OF DIMENSION L STORING THE X VALUES
!!  OF INPUT DATA POINTS (IN ASCENDING ORDER)
!!  Y=ARRAY OF DIMENSION L STORING THE Y VALUES OF INPUT DATA POINTS
!!  N=NUMBER OF POINTS AT WHICH INTERPOLATION OF THE Y-VALUES
!!  IS REQUIRED (MUST BE 1 OR GREATER)
!!  U=ARRAY OF DIMENSION N STORING THE X VALUES
!!  OF THE DESIRED POINTS
!!
!!  THE OUTPUT PARAMETER IS V=ARRAY OF DIMENSION N WHERE THE
!!  INTERPOLATED Y VALUES ARE TO BE DISPLAYED
!!
!! INPUTS
!!  CUBIC SPLINE INTERPOLATION
!!
!! OUTPUT
!!
!! NOTES
!!   This routine is deprecated and will be replaced by the other routines of this module.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

SUBROUTINE INTRPL(L,X,Y,N,U,V,dv,dv2,ideriv)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
!
      PARAMETER (NQQ=12000)

      COMMON/QQ/ QQ(4,NQQ)
      DIMENSION X(L),Y(L),U(N),V(N),DV(NQQ),DV2(NQQ)
      EQUIVALENCE (P0,X3),(Q0,Y3),(Q1,T3)
      REAL*8 M1,M2,M3,M4,M5
      EQUIVALENCE (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),&
     & (J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
!
!     PRELIMINARY PROCESSING

      L0=L
      LM1=L0-1
      LM2=LM1-1
      LP1=L0+1
      N0=N
      IF(N0.GT.NQQ) THEN
          NQQV=NQQ
          write(std_out,2089) NQQV,N0
!          CALL EXIT
      END IF
      IF(LM2.LT.0) GO TO 90
      IF(N0.LE.0) GO TO 91
      DO 11 I=2,L0

!     IF(X(I-1)-X(I))11,95,96
      IF(X(I-1)-X(I).EQ.0.0D0) GO TO 95
      IF(X(I-1)-X(I).GT.0.0D0) GO TO 96
   11 CONTINUE
      IPV=0
!
!***  MAIN LOOP
      FINT=0.0D0
      DO 80 K=1,N0
      UK=U(K)
!
!***  ROUTINE TO LOCATE THE DESIRED POINT
       IF(UK.GE.X(L0)) GO TO 26
      IF(UK.LT.X(1)) GO TO 25
      IMN=2
      IMX=L0
   21 I=(IMN+IMX)/2
      IF(UK.GE.X(I)) GO TO 23
      IMX=I
      GO TO 24
   23 IMN=I+1
   24 IF(IMX.GT.IMN) GO TO 21
      I=IMX
      GO TO 30
   25 I=1
      GO TO 30
   26 I=LP1
      GO TO 30
!
!***  CHECK IF I=IPV
   30 IF(I.EQ.IPV) GO TO 70
      IPV=I
!
!***  ROUTINES TO PICK UP NECESSARY X AND Y VALUES AND TO
!***  ESTIMATE THEM IF NECESSARY
      J=I
      IF(J.EQ.1) J=2
      IF(J.EQ.LP1) J=L0
      X3=X(J-1)
      Y3=Y(J-1)
      X4=X(J)
      Y4=Y(J)
      A3=X4-X3
      M3=(Y4-Y3)/A3
      IF(LM2.EQ.0) GO TO 43
      IF(J.EQ.2) GO TO 41
      X2=X(J-2)
      Y2=Y(J-2)
      A2=X3-X2
      M2=(Y3-Y2)/A2
      IF(J.EQ.L0) GO TO 42
   41 X5=X(J+1)
      Y5=Y(J+1)
      A4=X5-X4
      M4=(Y5-Y4)/A4
      IF(J.EQ.2) M2=M3+M3-M4
      GO TO 45
   42 M4=M3+M3-M2
      GO TO 45
   43 M2=M3
   45 IF(J.LE.3) GO TO 46
      A1=X2-X(J-3)
      M1=(Y2-Y(J-3))/A1
      GO TO 47
   46 M1=M2+M2-M3
   47 IF(J.GE.LM1) GO TO 48
      A5=X(J+2)-X5
      M5=(Y(J+2)-Y5)/A5
      GO TO 50
   48 M5=M4+M4-M3
!
!***  NUMERICAL DIFFERENTIATION
   50 IF(I.EQ.LP1) GO TO 52
      W2=ABS(M4-M3)
      W3=ABS(M2-M1)
      SW=W2+W3
      IF(SW.NE.0.0) GO TO 51
      W2=0.5D0
      W3=0.5D0
      SW=1.0D0
   51 T3=(W2*M2+W3*M3)/SW
      IF(I.EQ.1) GO TO 54
   52 W3=ABS(M5-M4)
      W4=ABS(M3-M2)
      SW=W3+W4
      IF(SW.NE.0.0) GO TO 53
      W3=0.5D0
      W4=0.5D0
      SW=1.0D0
   53 T4=(W3*M3+W4*M4)/SW
      IF(I.NE.LP1) GO TO 60
      T3=T4
      SA=A2+A3
      T4=0.5D0*(M4+M5-A2*(A2-A3)*(M2-M3)/(SA*SA))
      X3=X4
      Y3=Y4
      A3=A2
      M3=M4
      GO TO 60
   54 T4=T3
      SA=A3+A4
      T3=0.5D0*(M1+M2-A4*(A3-A4)*(M3-M4)/(SA*SA))
      X3=X3-A4
      Y3=Y3-M2*A4
      A3=A4
      M3=M2
!
!***  COMPUTATION OF THE POLYNOMIAL
   60 Q2=(2.0D0*(M3-T3)+M3-T4)/A3
      Q3=(-M3-M3+T3+T4)/(A3*A3)
   70 DX=UK-P0
      V(K)=Q0+DX*(Q1+DX*(Q2+DX*Q3))

      IF(IDERIV.EQ.0) GO TO 80
      DV(K)=Q1+DX*(2.0D0*Q2+DX*3.0D0*Q3)
      DV2(k)=6.0D0*Q3*DX+2.d0*Q2
      QQ(1,K)=Q0
      QQ(2,K)=Q1
      QQ(3,K)=Q2
      QQ(4,K)=Q3
   80  CONTINUE
      RETURN
!
!***  ERROR EXIT
   90 write(std_out,2090)
      GO TO 99
   91 write(std_out,2091)
      GO TO 99
   95 write(std_out,2095)
      GO TO 97
   96 write(std_out,2096)
   97 write(std_out,2097)I,X(I)
   99 write(std_out,2099) L0,N0
      RETURN
!
!***  FORMAT STATEMENTS
 2089  FORMAT( 'WARNING ERROR IN INTRPL. MAX ALLOWED VALUE OF N0 IS',&
     & I3,' HERE N0 IS',I3)
 2090  FORMAT(1X/' N = 1 OR LESS.'/)
 2091  FORMAT(1X/' N = 0 OR LESS.'/)
 2095  FORMAT(1X/' IDENTICAL X VALUES.'/)
 2096  FORMAT(1X/' X VALUES OUT OF SEQUENCE.'/)
 2097  FORMAT(4X,'I =',I7,10X,6X,'X(I) =',E12.3)
 2099  FORMAT(4X,'L =',I7,10X,3X,'N =',I7/ &
     & ' ERROR DETECTED IN ROUTINE INTRPL')
!
END subroutine intrpl
!!***

end module m_splines
!!***

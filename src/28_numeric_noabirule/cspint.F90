!{\src2tex{textfont=tt}}
!!****f* ABINIT/cspint
!! NAME
!!  cspint
!!
!! FUNCTION 
!!  CSPINT estimates the integral of a tabulated function.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES 
!!
!!    The routine is given the value of a function F(X) at a set of 
!!    nodes XTAB, and estimates
!!
!!      Integral ( A <= X <= B ) F(X) DX
!!
!!    by computing the cubic natural spline S(X) that interpolates
!!    F(X) at the nodes, and then computing
!!
!!      Integral ( A <= X <= B ) S(X) DX
!!
!!    exactly.
!!
!!    Other output from the program includes the definite integral
!!    from X(1) to X(I) of S(X), and the coefficients necessary for
!!    the user to evaluate the spline S(X) at any point.
!!
!!  Modified:
!!
!!    30 October 2000
!!
!!  Reference:
!!
!!    Philip Davis and Philip Rabinowitz,
!!    Methods of Numerical Integration,
!!    Blaisdell Publishing, 1967.
!!
!!  Parameters:
!!
!!    Input, real (dp) FTAB(NTAB), contains the tabulated values of
!!    the function, FTAB(I) = F(XTAB(I)).
!!
!!    Input, real (dp) XTAB(NTAB), contains the points at which the
!!    function was evaluated.  The XTAB's must be distinct and
!!    in ascending order.
!!
!!    Input, integer NTAB, the number of entries in FTAB and
!!    XTAB.  NTAB must be at least 3.
!!
!!    Input, real (dp) A, lower limit of integration.
!!
!!    Input, real (dp) B, upper limit of integration.
!!
!!    Output, real (dp) Y(3,NTAB), will contain the coefficients
!!    of the interpolating natural spline over each subinterval.
!!
!!    For XTAB(I) <= X <= XTAB(I+1),
!!
!!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!!                   + Y(2,I)*(X-XTAB(I))**2
!!                   + Y(3,I)*(X-XTAB(I))**3
!!
!!    Output, real (dp) E(NTAB), E(I) = the definite integral from
!!    XTAB(1) to XTAB(I) of S(X).
!!
!!    Workspace, real (dp) WORK(NTAB).
!!
!!    Output, real (dp) RESULT, the estimated value of the integral.
!!
!!
!! PARENTS
!!      m_xc_vdw,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cspint ( ftab, xtab, ntab, a, b, y, e, work, result )

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cspint'
!End of the abilint section

  implicit none

  integer, intent(in) :: ntab

  real(dp), intent(in) :: a
  real(dp), intent(in) :: b
  real(dp), intent(inout) :: e(ntab)
  real(dp), intent(in) :: ftab(ntab)
  integer :: i
  integer :: j
  real(dp) :: r
  real(dp), intent(out) :: result
  real(dp) :: s
  real(dp) :: term
  real(dp) :: u
  real(dp), intent(inout) :: work(ntab)
  real(dp), intent(in) :: xtab(ntab)
  real(dp), intent(inout) :: y(3,ntab)

  if ( ntab < 3 ) then
    write(std_out,'(a)' ) ' '
    write(std_out,'(a)' ) 'CSPINT - Fatal error!'
    write(std_out,'(a,i6)' ) '  NTAB must be at least 3, but input NTAB = ',ntab
    MSG_ERROR("Aborting now")
  end if
 
  do i = 1, ntab-1
 
    if ( xtab(i+1) <= xtab(i) ) then
      write(std_out,'(a)' ) ' '
      write(std_out,'(a)' ) 'CSPINT - Fatal error!'
      write(std_out,'(a)' ) '  Nodes not in strict increasing order.'
      write(std_out,'(a,i6)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
      write(std_out,'(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
      write(std_out,'(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
      MSG_ERROR("Aborting now")
    end if
 
  end do
 
  s = zero
  do i = 1, ntab-1
    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
    y(2,i) = r - s
    s = r
  end do
 
  result = zero
  s = zero
  r = zero
  y(2,1) = zero
  y(2,ntab) = zero
 
  do i = 2, ntab-1
    y(2,i) = y(2,i) + r * y(2,i-1)
    work(i) = two * ( xtab(i-1) - xtab(i+1) ) - r * s
    s = xtab(i+1) - xtab(i)
    r = s / work(i)
  end do
 
  do j = 2, ntab-1
    i = ntab+1-j
    y(2,i) = ( ( xtab(i+1) - xtab(i) ) * y(2,i+1) - y(2,i) ) / work(i)
  end do
 
  do i = 1, ntab-1
    s = xtab(i+1) - xtab(i)
    r = y(2,i+1) - y(2,i)
    y(3,i) = r / s
    y(2,i) = three * y(2,i)
    y(1,i) = ( ftab(i+1) - ftab(i) ) / s - ( y(2,i) + r ) * s
  end do
 
  e(1) = 0.0D+00
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    term = ((( y(3,i) * quarter * s + y(2,i) * third ) * s &
      + y(1,i) * half ) * s + ftab(i) ) * s
    e(i+1) = e(i) + term
  end do
!
!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!
  r = a
  u = one
 
  do j = 1, 2
!
!  The endpoint is less than or equal to XTAB(1).
!
    if ( r <= xtab(1) ) then
      result = result-u*((r-xtab(1))*y(1,1)*half +ftab(1))*(r-xtab(1))
!
!  The endpoint is greater than or equal to XTAB(NTAB).
!
    else if ( xtab(ntab) <= r ) then

      result = result -u * ( e(ntab) + ( r - xtab(ntab) ) &
        * ( ftab(ntab) + half * ( ftab(ntab-1) &
        + ( xtab(ntab) - xtab(ntab-1) ) * y(1,ntab-1) ) &
        * ( r - xtab(ntab) )))
!
!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!
    else

      do i = 1, ntab-1
 
        if ( r <= xtab(i+1) ) then
          r = r-xtab(i)
          result = result-u*(e(i)+(((y(3,i)*quarter*r+y(2,i)*third)*r &
            +y(1,i)*half )*r+ftab(i))*r)
          go to 120
        end if
 
      end do
 
    end if
 
  120   continue
 
    u = -one
    r = b
 
  end do
 
  return

end subroutine cspint
!!***

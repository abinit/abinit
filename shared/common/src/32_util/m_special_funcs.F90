!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_special_funcs
!! NAME
!! m_special_funcs
!!
!! FUNCTION
!! This module contains routines and functions used to
!! evaluate special functions frequently needed in Abinit.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, MT, FB, XG, MVer, FJ, NH, GZ, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_special_funcs

 use defs_basis
 use m_abicore
 use m_errors
 use m_splines

 use m_fstrings,        only : sjoin, ftoa
 use m_numeric_tools,   only : arth, simpson

 implicit none

 private

 public :: clp               ! x-1, if x>1/2, x+1, if x<-1/2
 public :: factorial         ! Calculates N! returning a real.
 public :: permutations      ! Returns N!/(N-k) if N>=0 and N-k>0 else 0.
 public :: binomcoeff        ! Binominal coefficient n!/(n-k)!
 public :: laguerre          ! Laguerre Polynomial(x,n,a).
 public :: RadFnH            ! Atomic radial function(r,n,l,Z).
 public :: iradfnh           ! Norm of atomic radial function(a,b,n,l,Z).
 public :: gaussian          ! Normalized Gaussian distribution.
 public :: lorentzian        ! Approximate Dirac Delta with lorentzian
 public :: abi_derf          ! Evaluates the error function in real(dp).
 public :: abi_derfc         ! Evaluates the complementary error function in real(dp).
 public :: gamma_function    ! Computes the gamma function
 public :: besjm             ! Spherical bessel function of order nn. Handles nn=0,1,2,3,4, or 5 only.
 public :: sbf8              ! Computes set of spherical bessel functions using accurate algorithm
 public :: k_fermi           ! Fermi wave vector corresponding to the local value of the real space density rhor.
 public :: k_thfermi         ! Thomas-Fermi wave vector corresponding to the local value of the real space density rhor
 public :: levi_civita_3     ! Return Levi-Civita tensor of rank 3
!!***

!!****t* m_special_funcs/jlspline_t
!! NAME
!! jlspline_t
!!
!! FUNCTION
!!  Object used to interpolate Bessel functions
!!
!! SOURCE

 public :: fermi_dirac        ! Fermi Dirac distribution
 public :: bose_einstein      ! Bose Einstein distribution

 type,public :: jlspline_t

   integer :: nx
   ! number of points on linear mesh used in spline.

   integer :: mlang
   ! mlang= max angular momentum + 1

   real(dp) :: delta
   ! Step of linear mesh.

   real(dp) :: maxarg
   ! max arg value.

   real(dp),allocatable :: xx(:)
   ! xx(nx)
   ! coordinates of points belonging to the grid

   real(dp),allocatable :: bess_spl(:,:)
   ! bess_spl(nx,mlang)
   ! bessel functions computed on the linear mesh

   real(dp),allocatable :: bess_spl_der(:,:)
   ! bess_spl_der(nx,mlang)
   ! the second derivatives of the cubic spline.

 end type jlspline_t
!!***

 public :: jlspline_new         ! Create new object.
 public :: jlspline_free        ! Free memory.
 public :: jlspline_integral    ! Compute integral.

!!****t* m_special_funcs/gspline_t
!! NAME
!! gspline_t
!!
!! FUNCTION
!!  Object used to interpolate the gaussian approximant and its primitive with cubic spline.
!!  Particularly useful if we are computing DOSes with many k-points/bands
!!  because one can significantly decrease the number of calls to exponential functions.
!!
!! SOURCE

 type,public :: gspline_t

   integer :: nspline
    ! Number of points used in spline table.

   real(dp) :: sigma
    ! Broadening parameter.

   real(dp) :: xmin,xmax
    ! Min and max x in spline mesh. Only positive xs are stored in memory
    ! The values at -x are reconstructed by symmetry.
    ! xmin is usually zero, xmax is the point where the gaussian == tol16.
    ! g(x) is set to zero if x > xmin.

   real(dp) :: step, stepm1, step2div6
    ! Step of the linear mesh used in spline and associated coeffients.

   real(dp),allocatable :: xvals(:)
    ! xvals(nspline)
    ! The xvalues used in the spline

   real(dp),allocatable :: svals(:,:)
    ! svals(nspline,4)
    ! Internal tables with spline data.

 end type gspline_t
!!***

 public :: gspline_new       ! Creation method.
 public :: gspline_eval      ! Evaluate interpolant
 public :: gspline_free      ! Free memory.

CONTAINS  !===========================================================
!!***

!!****f* m_special_funcs/clp
!! NAME
!! clp
!!
!! FUNCTION
!! clp(x)= x-1, if x>1/2
!!         x+1, if x<-1/2
!!
!! INPUTS
!!  x= input variable
!!
!! OUTPUT
!!  clp= resulting function
!!
!! PARENTS
!!      nhatgrid
!!
!! SOURCE

pure function clp(x)

!Arguments ------------------------------------
!scalars
 real(dp) :: clp
 real(dp),intent(in) :: x

! **********************************************************************

 if(x > half) then
   clp=x-one
 elseif(x < -half) then
   clp=x+one
 else
   clp=x
 end if

end function clp
!!***

!!****f* m_special_funcs/factorial
!! NAME
!! factorial
!!
!! FUNCTION
!! Calculates N!. Returns a (dp) real.
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   factorial= n! (real)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function factorial(nn)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: factorial

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: ff

! *********************************************************************

 ff=one
 do ii=2,nn
   ff=ff*ii
 end do

 factorial=ff

end function factorial
!!***

!!****f* m_special_funcs/permutations
!! NAME
!! permutations
!!
!! FUNCTION
!! Returns N!/(N-k)!  if N>=0 and N-k>0
!!                    otherwise 0 is returned
!! Output is real
!!
!! INPUTS
!!   kk=number k to use
!!   nn=number N to use
!!
!! OUTPUT
!!   permutations= n!/(n-k)! (real)
!!
!! PARENTS
!!      green_atomic_hubbard
!!
!! CHILDREN
!!
!! SOURCE

pure function permutations(nn,kk)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kk,nn
 real(dp) :: permutations

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: pp

! *********************************************************************

 if ((nn>=0).and.((nn-kk)>=0)) then
   pp=one
   do ii=nn-kk+1,nn
     pp=pp*ii
   end do
 else
   pp=zero
 end if

 permutations=pp

end function permutations
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/binomcoeff
!! NAME
!! factorial
!!
!! FUNCTION
!! Calculates n!/( k!* (n-k)!). Returns a real (dp)
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   binomcoeff= n!/( k!* (n-k)!)  (real dp)
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!! SOURCE

elemental function binomcoeff(n,k)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: n,k
 real(dp) :: binomcoeff

! *********************************************************************

 binomcoeff=factorial(n)/(factorial(k)*factorial(n-k))

end function binomcoeff
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/laguerre
!! NAME
!! laguerre
!!
!! FUNCTION
!! Laguerre(x,n,a). Returns a (dp) real.
!!
!! INPUTS
!!   x position
!!   n order of laguerre polynomial
!!   a
!!
!! OUTPUT
!!   Laguerre(x,n,a) (dp)
!!
!! PARENTS
!!
!!
!! CHILDREN
!!   factorial
!!
!! SOURCE

function laguerre(x,n,a)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in),optional :: n,a
 real(dp)                    :: laguerre
 real(dp),intent(in)         :: x

!Local variables ---------------------------------------
!scalars
 integer :: ii, nn, aa

!arrays
 real(dp),allocatable :: ff(:)

! *********************************************************************

 if (present(n)) then
   nn=n
 else
   nn=1
 end if

 if (present(a)) then
   aa=a
 else
   aa=0
 end if
 ABI_ALLOCATE(ff,(nn+1))
 ff=0.0_dp
 ff=(/ (binomcoeff(nn+aa,nn-ii)*((-1.0_dp)*x)**ii/factorial(ii) ,ii=0,nn) /)
 laguerre=sum(ff)

 ABI_DEALLOCATE(ff)

end function laguerre
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/RadFnH
!! NAME
!! RadFnH
!!
!! FUNCTION
!!  RadFnH(r,n,l,Z) radial function of atomic wavefunction with nuclear charge Z.
!!  for quantum number n, and l.
!!  Default: Fe 3d function. Returns a (dp) real.
!!
!! INPUTS
!!   r radius
!!   n principal quantum number
!!   l quantum number
!!
!! OUTPUT
!!  RadFnH(r,n,l,Z) (dp)
!!
!! PARENTS
!!
!!
!! CHILDREN
!!  Laguerre
!!  factorial
!!
!! SOURCE


function RadFnH(r,n,l,Z)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in),optional :: n,l
 real(dp) :: RadFnH
 real(dp),intent(in) :: r
 real(dp),intent(in),optional :: Z

!Local variables ---------------------------------------
!scalars
 integer   :: nn,ll
 real(dp)  :: ff,rr,ZZ

! *********************************************************************

 if (present(n)) then
   nn=n
 else
   nn=3
 end if

 if (present(l)) then
   ll=l
 else
   ll=2
 end if

 if (present(Z)) then
   ZZ=Z
 else
   ZZ=28.0_dp
 end if

 rr=ZZ*r/nn
 ff=exp(log(ZZ*1.0_dp)*(3.0_dp/2.0_dp))*2/nn**2
 ff=ff*sqrt(factorial(nn-ll-1)/factorial(nn+ll))*(2*rr)**ll
 RadFnH=ff*exp(-1*rr)*laguerre(2*rr,nn-ll-1,2*ll+1)

end function RadFnH
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/IRadFnH
!! NAME
!! IRadFnH
!!
!! FUNCTION
!!  IRadFnH(a,b,n,l,Z): Integral of radial function of atomic wavefunction between a and b.
!!  recursive programming using simpson's rule
!!  iteration depth of m=8 corresponds to relative error of 10^(-12).
!!
!! INPUTS
!!   a lower limit for integration
!!   b upper limit for integration
!!   n principal quantum number
!!   l quantum number
!    Z nuclear charge
!!
!! OUTPUT
!!  IRadFnH(a,b,n,l,Z) (dp)
!!
!! PARENTS
!!
!!
!! CHILDREN
!!  Laguerre
!!  factorial
!!  RadFnH
!!
!! SOURCE

recursive function IRadFnH(a,b,n,l,Z,m) result(x)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in),optional  :: n,l,m
 real(dp),intent(in):: a
 real(dp),intent(in),optional :: b,Z

!Local variables ---------------------------------------
!scalars
 integer   :: nn,ll,mm
 real(dp)  :: h,bb,ZZ,x

! *********************************************************************

 if (present(n)) then
   nn=n
 else
   nn=3
 end if

 if (present(l)) then
   ll=l
 else
   ll=2
 end if

 if (present(Z)) then
   ZZ=Z
 else
   ZZ=28
 end if

 if (present(b)) then
   bb=b
 else
   bb=100.0_dp
 end if

 if (present(m)) then
   mm=m
 else
   mm=0
 end if

 h=(bb-a)/2.0_dp
 if (mm<8) then
  !h=2*h/exp(1.0_dp)
  x=IRadFnH(a,a+h,nn,ll,ZZ,mm+1)+IRadFnH(a+h,bb,nn,ll,ZZ,mm+1)
 else
  x=RadFnH(a,nn,ll,ZZ)**2*a**2+4.0_dp*RadFnH(a+h,nn,ll,ZZ)**2*(a+h)**2
  x=h/3.0_dp*(x+RadFnH(bb,nn,ll,ZZ)**2*bb**2)
 end if

end function IRadFnH
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/gaussian
!! NAME
!! gaussian
!!
!! FUNCTION
!!  Return the values of the normalized Gaussian distribution
!!    Gauss(arg,sigma) =  1/(sigma SQRT(2*pi)) e^{-arg**2/(2*sigma**2)}
!!
!! INPUTS
!!   arg=Argument of the Gaussian.
!!   sigma=Standard deviation
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function gaussian(arg, sigma)

!Arguments ---------------------------------------------
!scalars
 real(dp),intent(in) :: arg,sigma
 real(dp) :: gaussian

!Local variables ---------------------------------------
!scalars
 real(dp) :: xx

! *********************************************************************

 xx = arg / (sqrt2 * sigma)
 gaussian = exp(-xx*xx) / (sigma * sqrt(two_pi))

end function gaussian
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/lorentzian
!! NAME
!! lorentzian
!!
!! FUNCTION
!!  Approximate Dirac Delta with lorentzian
!!
!! INPUTS
!!   arg=Argument of the lorentzian.
!!   sigma=Broadening factor
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function lorentzian(arg, sigma)

!Arguments ---------------------------------------------
!scalars
 real(dp),intent(in) :: arg, sigma
 real(dp) :: lorentzian

! *********************************************************************

 lorentzian = piinv * sigma / (arg ** 2 + sigma ** 2)

end function lorentzian
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/abi_derf
!! NAME
!! abi_derf
!!
!! FUNCTION
!! Evaluates the error function in real(dp).
!! Same implementation as imsl.
!! Simple mod of derfc.F90
!!
!! INPUTS
!! yy
!!
!! OUTPUT
!! derf_yy= error function of yy
!!
!! PARENTS
!!      evdw_wannier,wvl_wfs_set
!!
!! CHILDREN
!!
!! SOURCE

elemental function abi_derf(yy) result(derf_yy)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: yy
 real(dp) :: derf_yy

!Local variables-------------------------------
 integer          ::  done,ii,isw
! coefficients for 0.0 <= yy < .477
 real(dp), parameter :: &
&  pp(5)=(/ 113.8641541510502e0_dp, 377.4852376853020e0_dp,  &
&           3209.377589138469e0_dp, .1857777061846032e0_dp,  &
&           3.161123743870566e0_dp /)
 real(dp), parameter :: &
&  qq(4)=(/ 244.0246379344442e0_dp, 1282.616526077372e0_dp,  &
&           2844.236833439171e0_dp, 23.60129095234412e0_dp/)
! coefficients for .477 <= yy <= 4.0
 real(dp), parameter :: &
&  p1(9)=(/ 8.883149794388376e0_dp, 66.11919063714163e0_dp,  &
&           298.6351381974001e0_dp, 881.9522212417691e0_dp,  &
&           1712.047612634071e0_dp, 2051.078377826071e0_dp,  &
&           1230.339354797997e0_dp, 2.153115354744038e-8_dp, &
&           .5641884969886701e0_dp /)
 real(dp), parameter :: &
&  q1(8)=(/ 117.6939508913125e0_dp, 537.1811018620099e0_dp,  &
&           1621.389574566690e0_dp, 3290.799235733460e0_dp,  &
&           4362.619090143247e0_dp, 3439.367674143722e0_dp,  &
&           1230.339354803749e0_dp, 15.74492611070983e0_dp/)
  ! coefficients for 4.0 < y,
 real(dp), parameter :: &
&  p2(6)=(/ -3.603448999498044e-01_dp, -1.257817261112292e-01_dp,   &
&           -1.608378514874228e-02_dp, -6.587491615298378e-04_dp,   &
&           -1.631538713730210e-02_dp, -3.053266349612323e-01_dp/)
 real(dp), parameter :: &
&  q2(5)=(/ 1.872952849923460e0_dp   , 5.279051029514284e-01_dp,    &
&           6.051834131244132e-02_dp , 2.335204976268692e-03_dp,    &
&           2.568520192289822e0_dp /)
 real(dp), parameter :: &
&  sqrpi=.5641895835477563e0_dp, xbig=13.3e0_dp, xlarge=6.375e0_dp, xmin=1.0e-10_dp
 real(dp) ::  res,xden,xi,xnum,xsq,xx

! ******************************************************************

 xx = yy
 isw = 1
!Here change the sign of xx, and keep track of it thanks to isw
 if (xx<0.0e0_dp) then
   isw = -1
   xx = -xx
 end if

 done=0

!Residual value, if yy < -6.375e0_dp
 res=-1.0e0_dp

!abs(yy) < .477, evaluate approximation for erfc
 if (xx<0.477e0_dp) then
!  xmin is a very small number
   if (xx<xmin) then
     res = xx*pp(3)/qq(3)
   else
     xsq = xx*xx
     xnum = pp(4)*xsq+pp(5)
     xden = xsq+qq(4)
     do ii = 1,3
       xnum = xnum*xsq+pp(ii)
       xden = xden*xsq+qq(ii)
     end do
     res = xx*xnum/xden
   end if
   if (isw==-1) res = -res
   done=1
 end if

!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
 if (xx<=4.0e0_dp .and. done==0 ) then
   xsq = xx*xx
   xnum = p1(8)*xx+p1(9)
   xden = xx+q1(8)
   do ii=1,7
     xnum = xnum*xx+p1(ii)
     xden = xden*xx+q1(ii)
   end do
   res = xnum/xden
   res = res* exp(-xsq)
   if (isw.eq.-1) then
     res = res-1.0e0_dp
   else
     res=1.0e0_dp-res
   end if
   done=1
 end if

!y > 13.3e0_dp
 if (isw > 0 .and. xx > xbig .and. done==0 ) then
   res = 1.0e0_dp
   done=1
 end if

!4.0 < yy < 13.3e0_dp  .or. -6.375e0_dp < yy < -4.0
!evaluate minimax approximation for erfc
 if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
   xsq = xx*xx
   xi = 1.0e0_dp/xsq
   xnum= p2(5)*xi+p2(6)
   xden = xi+q2(5)
   do ii = 1,4
     xnum = xnum*xi+p2(ii)
     xden = xden*xi+q2(ii)
   end do
   res = (sqrpi+xi*xnum/xden)/xx
   res = res* exp(-xsq)
   if (isw.eq.-1) then
     res = res-1.0e0_dp
   else
     res=1.0e0_dp-res
   end if
 end if

!All cases have been investigated
 derf_yy = res

end function abi_derf
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/abi_derfc
!! NAME
!! abi_derfc
!!
!! FUNCTION
!! Evaluates the complementary error function in real(dp).
!! Same implementation as imsl.
!!
!! INPUTS
!! yy
!!
!! OUTPUT
!! derfc_yy=complementary error function of yy
!!
!! PARENTS
!!      ewald,ewald2,dfpt_ewald,elt_ewald,ewald9,make_efg_ion,psp2lo,wvl_wfs_set
!!
!! CHILDREN
!!
!! SOURCE

elemental function abi_derfc(yy) result(derfc_yy)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: yy
 real(dp) :: derfc_yy

!Local variables-------------------------------
 integer          ::  done,ii,isw
! coefficients for 0.0 <= yy < .477
 real(dp), parameter :: &
&  pp(5)=(/ 113.8641541510502e0_dp, 377.4852376853020e0_dp,  &
&           3209.377589138469e0_dp, .1857777061846032e0_dp,  &
&           3.161123743870566e0_dp /)
 real(dp), parameter :: &
&  qq(4)=(/ 244.0246379344442e0_dp, 1282.616526077372e0_dp,  &
&           2844.236833439171e0_dp, 23.60129095234412e0_dp/)
! coefficients for .477 <= yy <= 4.0
 real(dp), parameter :: &
&  p1(9)=(/ 8.883149794388376e0_dp, 66.11919063714163e0_dp,  &
&           298.6351381974001e0_dp, 881.9522212417691e0_dp,  &
&           1712.047612634071e0_dp, 2051.078377826071e0_dp,  &
&           1230.339354797997e0_dp, 2.153115354744038e-8_dp, &
&           .5641884969886701e0_dp /)
 real(dp), parameter :: &
&  q1(8)=(/ 117.6939508913125e0_dp, 537.1811018620099e0_dp,  &
&           1621.389574566690e0_dp, 3290.799235733460e0_dp,  &
&           4362.619090143247e0_dp, 3439.367674143722e0_dp,  &
&           1230.339354803749e0_dp, 15.74492611070983e0_dp/)
 ! coefficients for 4.0 < y,
 real(dp), parameter :: &
&  p2(6)=(/ -3.603448999498044e-01_dp, -1.257817261112292e-01_dp,   &
&           -1.608378514874228e-02_dp, -6.587491615298378e-04_dp,   &
&           -1.631538713730210e-02_dp, -3.053266349612323e-01_dp/)
 real(dp), parameter :: &
&  q2(5)=(/ 1.872952849923460e0_dp   , 5.279051029514284e-01_dp,    &
&           6.051834131244132e-02_dp , 2.335204976268692e-03_dp,    &
&           2.568520192289822e0_dp /)
 real(dp), parameter :: &
&  sqrpi=.5641895835477563e0_dp, xbig=13.3e0_dp, xlarge=6.375e0_dp, xmin=1.0e-10_dp
 real(dp) ::  res,xden,xi,xnum,xsq,xx

!******************************************************************

 xx = yy
 isw = 1
!Here change the sign of xx, and keep track of it thanks to isw
 if (xx<0.0e0_dp) then
   isw = -1
   xx = -xx
 end if

 done=0

!Residual value, if yy < -6.375e0_dp
 res=2.0e0_dp

!abs(yy) < .477, evaluate approximation for erfc
 if (xx<0.477e0_dp) then
!  xmin is a very small number
   if (xx<xmin) then
     res = xx*pp(3)/qq(3)
   else
     xsq = xx*xx
     xnum = pp(4)*xsq+pp(5)
     xden = xsq+qq(4)
     do ii = 1,3
       xnum = xnum*xsq+pp(ii)
       xden = xden*xsq+qq(ii)
     end do
     res = xx*xnum/xden
   end if
   if (isw==-1) res = -res
   res = 1.0e0_dp-res
   done=1
 end if

!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
 if (xx<=4.0e0_dp .and. done==0 ) then
   xsq = xx*xx
   xnum = p1(8)*xx+p1(9)
   xden = xx+q1(8)
   do ii=1,7
     xnum = xnum*xx+p1(ii)
     xden = xden*xx+q1(ii)
   end do
   res = xnum/xden
   res = res* exp(-xsq)
   if (isw.eq.-1) res = 2.0e0_dp-res
   done=1
 end if

!y > 13.3e0_dp
 if (isw > 0 .and. xx > xbig .and. done==0 ) then
   res = 0.0e0_dp
   done=1
 end if

!4.0 < yy < 13.3e0_dp  .or. -6.375e0_dp < yy < -4.0
!evaluate minimax approximation for erfc
 if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
   xsq = xx*xx
   xi = 1.0e0_dp/xsq
   xnum= p2(5)*xi+p2(6)
   xden = xi+q2(5)
   do ii = 1,4
     xnum = xnum*xi+p2(ii)
     xden = xden*xi+q2(ii)
   end do
   res = (sqrpi+xi*xnum/xden)/xx
   res = res* exp(-xsq)
   if (isw.eq.-1) res = 2.0e0_dp-res
 end if

!All cases have been investigated
 derfc_yy = res

end function abi_derfc
!!***

!!****f* ABINIT/GAMMA_FUNCTION
!! NAME
!!  GAMMA_FUNCTION
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_rec,m_spin_current
!!
!! CHILDREN
!!
!! SOURCE

subroutine GAMMA_FUNCTION(X,GA)

!       ====================================================
!       Purpose: This program computes the gamma function
!                Gamma(x) using subroutine GAMMA
!       Examples:
!                   x          Gamma(x)
!                ----------------------------
!                  1/3       2.678938534708
!                  0.5       1.772453850906
!                 -0.5      -3.544907701811
!                 -1.5       2.363271801207
!                  5.0      24.000000000000
!       ====================================================
!
!  This routine was downloaded from UIUC:
!  http://jin.ece.uiuc.edu/routines/routines.html
!
!  The programs appear to accompany a book "Computation of Special
!  Functions" (1996) John Wiley and Sons, but are distributed online
!  by the authors. Exact copyright should be checked.
!
!  Authors / copyright:
!     Shanjie Zhang and Jianming Jin
!     Proposed contact is:  j-jin1@uiuc.edu
!
!  20 October 2008:
!     Incorporated into ABINIT by M. Verstraete
!
!
!
!       ==================================================
!       Purpose: Compute the gamma function Gamma(x)
!       Input :  x  --- Argument of Gamma(x)
!                       ( x is not equal to 0,-1,-2, etc )
!       Output:  GA --- Gamma(x)
!       ==================================================
!

  ! arguments

  real(dp),intent(in) :: x
  real(dp),intent(out) :: ga

  ! local variables
  integer :: k,m
  real(dp) :: m1,z,r,gr
  real(dp) :: G(26)

  ! source code:

  ! initialization of reference data
  G=(/1.0D0,0.5772156649015329D0, &
     &  -0.6558780715202538D0, -0.420026350340952D-1, &
     &   0.1665386113822915D0,-.421977345555443D-1, &
     &  -.96219715278770D-2, .72189432466630D-2, &
     &  -.11651675918591D-2, -.2152416741149D-3, &
     &   .1280502823882D-3, -.201348547807D-4, &
     &  -.12504934821D-5, .11330272320D-5, &
     &  -.2056338417D-6, .61160950D-8, &
     &   .50020075D-8, -.11812746D-8, &
     &   .1043427D-9, .77823D-11, &
     &  -.36968D-11, .51D-12, &
     &  -.206D-13, -.54D-14, .14D-14, .1D-15/)


  ! for the integer case, do explicit factorial
  if (X==int(X)) then
    if (X > 0.0D0) then
      GA=1.0D0
      M1=X-1
      do K=2,int(M1)
        GA=GA*K
      end do
    else
      GA=1.0D+300
    end if
  ! for the integer case, do explicit factorial
  else
    if (abs(X) > 1.0D0) then
      Z=abs(X)
      M=int(Z)
      R=1.0D0
      do K=1,M
        R=R*(Z-K)
      end do
      Z=Z-M
    else
      Z=X
    end if
    GR=G(26)
    do K=25,1,-1
      GR=GR*Z+G(K)
    end do
    GA=1.0D0/(GR*Z)
    if (abs(X) > 1.0D0) then
      GA=GA*R
      if (X < 0.0D0) GA=-PI/(X*GA*SIN(PI*X))
    end if
  end if
  return

end subroutine GAMMA_FUNCTION
!!***

!!****f* m_special_funcs/besjm
!! NAME
!! besjm
!!
!! FUNCTION
!! Spherical bessel function of order nn. Handles nn=0,1,2,3,4, or 5 only.
!!
!! INPUTS
!!  arg= scaling to be applied to xx(nx)
!!  nn=order of spherical bessel function (only 0 through 5 allowed)
!!  cosx(1:nx)=cosines of arg*xx(1:nx)
!!  xx(1:nx)=set of dimensionless arguments of function
!!  nx=number of arguments
!!  sinx(1:nx)=sines of arg*xx(1:nx)
!!
!! OUTPUT
!!  besjx(1:nx)=returned values
!!
!! NOTES
!! besj(nn,y)=$ j_{nn}(y) =(\frac{\pi}{2y})^{\frac{1}{2}}J(nn+\frac{1}{2},y)$
!! where J=Bessel function of the first kind.
!! besjm compute multiple values, and relies on precomputed values of sin and cos of y.
!! The argument y is arg*xx(ix), for ix from 1 to nx
!! The values of xx must be positive, and ordered by increasing order
!! At small arg, the higher orders have so much cancellation that the
!! analytic expression is very poor computationally.  In that case we
!! use a rational polynomial approximation.
!!
!! PARENTS
!!      m_mlwfovlp,m_psp1,m_special_funcs
!!
!! CHILDREN
!!
!! SOURCE

subroutine besjm(arg,besjx,cosx,nn,nx,sinx,xx)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn,nx
 real(dp),intent(in) :: arg
!arrays
 real(dp),intent(in) :: cosx(nx),sinx(nx),xx(nx)
 real(dp),intent(out) :: besjx(nx)

!Local variables-------------------------------
!scalars
 integer :: ix,switchx
!Series or rational polynomial coefficients
 real(dp),parameter :: b01=1.d0/6.d0,b02=1.d0/120.d0,b03=1.d0/5040.d0
 real(dp),parameter :: b04=1.d0/362880.d0,b11=0.8331251468724171d-1
 real(dp),parameter :: b12=0.2036961284395412d-2,b13=0.1932970379901801d-4
 real(dp),parameter :: b14=0.6526053169009489d-7,b21=0.5867824627555163d-1
 real(dp),parameter :: b22=0.1152501878595934d-2,b23=0.1011071389414764d-4
 real(dp),parameter :: b24=0.4172322111421287d-7,b25=0.6790616688656543d-10
 real(dp),parameter :: b31=0.439131885807176d-1,b32=0.6813139609887099d-3
 real(dp),parameter :: b33=0.4899103784264755d-5,b34=0.17025590795625d-7
 real(dp),parameter :: b35=0.2382642910613347d-10,b41=0.3587477991030971d-1
 real(dp),parameter :: b42=0.4833719855268907d-3,b43=0.3238388977796242d-5
 real(dp),parameter :: b44=0.1171802513125112d-7,b45=0.223261650431992d-10
 real(dp),parameter :: b46=.1800045587335951d-13,b51=0.295232406376567d-1
 real(dp),parameter :: b52=0.3359864457080573d-3,b53=0.19394750603618d-5
 real(dp),parameter :: b54=0.6143166228216219d-8,b55=0.10378501636108d-10
 real(dp),parameter :: b56=.749975122872713d-14
 real(dp),parameter :: c11=0.1668748531275829d-1,c12=0.1342812442426702d-3
 real(dp),parameter :: c13=0.6378249315355233d-6,c14=0.1573564527360138d-8
 real(dp),parameter :: c21=0.127503251530198d-1,c22=0.7911240539893565d-4
 real(dp),parameter :: c23=0.3044380758068054d-6,c24=0.7439837832363479d-9
 real(dp),parameter :: c25=0.9515065658793124d-12,c31=0.1164236697483795d-1
 real(dp),parameter :: c32=0.654858636312224d-4,c33=0.2265576367562734d-6
 real(dp),parameter :: c34=0.4929905563217352d-9,c35=0.555120465710914d-12
 real(dp),parameter :: c41=0.9579765544235745d-2,c42=0.4468999977536864d-4
 real(dp),parameter :: c43=0.1315634305905896d-6,c44=0.2615492488301639d-9
 real(dp),parameter :: c45=0.3387473312408129d-12,c46=.2280866204624012d-15
 real(dp),parameter :: c51=0.8938297823881763d-2,c52=0.3874149021633025d-4
 real(dp),parameter :: c53=0.1054692715135225d-6,c54=0.192879620987602d-9
 real(dp),parameter :: c55=0.2284469423833734d-12,c56=0.139729234332572d-15
 real(dp),parameter :: ffnth=1.d0/15.d0,o10395=1.d0/10395d0,oo105=1.d0/105.d0
 real(dp),parameter :: oo945=1.d0/945.d0
 real(dp) :: bot,rr,rsq,top
 character(len=500) :: message

! *************************************************************************

 if (nn==0) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=1.d-1) then
       rsq=rr*rr
       besjx(ix)=1.d0-rsq*(b01-rsq*(b02-rsq*(b03-rsq*b04)))
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     besjx(ix)=sinx(ix)/rr
   end do

 else if (nn==1) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=1.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b11-rsq*(b12-rsq*(b13-rsq*b14)))
       bot=1.d0+rsq*(c11+rsq*(c12+rsq*(c13+rsq*c14)))
       besjx(ix)=third*rr*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=(sinx(ix)-rr*cosx(ix))/rsq
   end do

 else if (nn==2) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=2.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b21-rsq*(b22-rsq*(b23-rsq*(b24-rsq*b25))))
       bot=1.d0+rsq*(c21+rsq*(c22+rsq*(c23+rsq*(c24+rsq*c25))))
       besjx(ix)=ffnth*rsq*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=((3.d0-rsq)*sinx(ix)-3.d0*rr*cosx(ix))/(rr*rsq)
   end do

 else if (nn==3) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=2.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b31-rsq*(b32-rsq*(b33-rsq*(b34-rsq*b35))))
       bot=1.d0+rsq*(c31+rsq*(c32+rsq*(c33+rsq*(c34+rsq*c35))))
       besjx(ix)=rr*rsq*oo105*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=( (15.d0-6.d0*rsq)*sinx(ix)&
&     + rr*(rsq-15.d0)  *cosx(ix) ) /(rsq*rsq)
   end do

 else if (nn==4) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=4.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b41-rsq*(b42-rsq*(b43-rsq*(b44-rsq*(b45-rsq*b46)))))
       bot=1.d0+rsq*(c41+rsq*(c42+rsq*(c43+rsq*(c44+rsq*(c45+rsq*c46)))))
       besjx(ix)=rsq*rsq*oo945*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=( (105.d0-rsq*(45.d0-rsq)) *sinx(ix)&
&     + rr * (10.d0*rsq-105.d0)  *cosx(ix) ) /(rsq*rsq*rr)
   end do

 else if (nn==5) then

   switchx=nx+1
   do ix=1,nx
     rr=arg*xx(ix)
     if (rr<=4.d0) then
       rsq=rr*rr
       top=1.d0-rsq*(b51-rsq*(b52-rsq*(b53-rsq*(b54-rsq*(b55-rsq*b56)))))
       bot=1.d0+rsq*(c51+rsq*(c52+rsq*(c53+rsq*(c54+rsq*(c55+rsq*c56)))))
       besjx(ix)=rsq*rsq*rr*o10395*top/bot
     else
       switchx=ix
       exit
     end if
   end do

   do ix=switchx,nx
     rr=arg*xx(ix)
     rsq=rr*rr
     besjx(ix)=( (945.d0-rsq*(420.d0-rsq*15.d0)) *sinx(ix)&
&     + rr * (945.d0-rsq*(105.d0-rsq))  *cosx(ix) ) /(rsq*rsq*rr)
   end do

 else
   write(message, '(a,i0,a)' )' besjm only defined for nn in [0,5]; input was nn=',nn,'.'
   MSG_BUG(message)
 end if

end subroutine besjm
!!***

!!****f* m_special_funcs/sbf8
!! NAME
!! sbf8
!!
!! FUNCTION
!! Computes set of spherical bessel functions using accurate algorithm
!! based on downward recursion in order and normalization sum.
!! Power series used at small arguments.
!!
!! INPUTS
!!  nm=maximum angular momentum wanted + one
!!  xx=argument of sbf
!!
!! OUTPUT
!!  sb_out(nm)=values of spherical bessel functions for l=0,nm-1
!!
!! PARENTS
!!      m_forctqmc,m_paw_overlap,m_positron,m_psptk
!!
!! CHILDREN
!!
!! SOURCE

subroutine sbf8(nm,xx,sb_out)

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: nm
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: sb_out(nm)

!Local variables-------------------------------
!scalars
 integer :: nlim,nn
 real(dp) :: fn,sn,xi,xn,xs
!arrays
 real(dp),allocatable :: sb(:)

! *************************************************************************

 if(xx<= 1.0e-36_dp) then
!  zero argument section
   sb_out(:)=zero
   sb_out(1)=one
 else if(xx<1.e-3_dp) then
!  small argument section
   xn=one
   xs=half*xx**2
   do nn=1,nm
     sb_out(nn)=xn*(one - xs*(one - xs/(4*nn+6))/(2*nn+1))
     xn=xx*xn/(2*nn+1)
   end do
 else
!  recursion method
   if(xx<one) then
     nlim=nm+int(15.0e0_dp*xx)+1
   else
     nlim=nm+int(1.36e0_dp*xx)+15
   end if
   ABI_ALLOCATE(sb,(nlim+1))
   nn=nlim
   xi=one/xx
   sb(nn+1)=zero
   sb(nn)=1.e-18_dp
   sn=dble(2*nn-1)*1.e-36_dp
   do nn=nlim-1,1,-1
     sb(nn)=dble(2*nn+1)*xi*sb(nn+1) - sb(nn+2)
   end do
   do nn=1,nlim-1
     sn=sn + dble(2*nn-1)*sb(nn)*sb(nn)
   end do
   fn=1.d0/sqrt(sn)
   sb_out(:)=fn*sb(1:nm)
   ABI_DEALLOCATE(sb)
 end if

end subroutine sbf8
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/fermi_dirac
!! NAME
!!  fermi_dirac
!!
!! FUNCTION
!!  Returns the Fermi Dirac distribution for T and energy wrt Fermi level
!!  presumes everything is in Hartree!!!! Not Kelvin for T
!!
!! INPUTS
!!  energy = electron energy level
!!  mu = chemical potential
!!  temperature = T
!!
!! PARENTS
!!
!! SOURCE

function fermi_dirac(energy, mu, temperature)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: energy, mu, temperature
 real(dp) :: fermi_dirac

!Local variables-------------------------------
!scalars
 real(dp) :: arg

! *************************************************************************

 fermi_dirac = zero
 if (temperature > tol12) then
   arg = (energy-mu)/temperature
   if(arg < -600._dp)then ! far below Ef
     fermi_dirac = one
   else if (arg < 600._dp)then ! around Ef
     fermi_dirac = one / (exp(arg)  + one)
   end if
 else  ! T is too small - just step function
   if (mu-energy > tol12) fermi_dirac = one
 end if

end function fermi_dirac
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/bose_einstein
!! NAME
!!  bose_einstein
!!
!! FUNCTION
!!  Returns the Bose Einstein distribution for T and energy
!!  presumes everything is in Hartree!!!! Not Kelvin for T
!!
!! INPUTS
!!  energy = electron energy level
!!  temperature = T
!!
!! PARENTS
!!
!! SOURCE

function bose_einstein(energy, temperature)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: energy, temperature
 real(dp) :: bose_einstein

!Local variables-------------------------------
!scalars
 real(dp) :: arg
 character(len=500) :: message

! *************************************************************************

 bose_einstein = zero
 if (temperature > tol12) then
   arg = energy/temperature
   if(arg > tol12 .and. arg < 600._dp)then
     bose_einstein = one / (exp(arg)  - one)
   else if (arg < tol12) then
     write(message,'(a)') 'No Bose Einstein for negative energies'
     MSG_WARNING(message)
   end if
 else if (arg < tol12) then
   write(message,'(a)') 'No Bose Einstein for negative or 0 T'
   MSG_WARNING(message)
 end if


end function bose_einstein
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/k_fermi
!! NAME
!!  k_fermi
!!
!! FUNCTION
!!  Returns the Fermi wave vector corresponding to the local value of the real space density rhor.
!!
!! INPUTS
!!  rhor=Local density in real space.
!!
!! PARENTS
!!
!! SOURCE

elemental function k_fermi(rhor)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rhor
 real(dp) :: k_fermi

!Local variables-------------------------------
!scalars
 real(dp),parameter :: pisq=pi**2

! *************************************************************************

 k_fermi = (three*pisq*rhor)**third

end function k_fermi
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/k_thfermi
!! NAME
!!  k_thfermi
!!
!! FUNCTION
!!  Returns the Thomas-Fermi wave vector corresponding to the local value of the real space density rhor.
!!
!! INPUTS
!!  rhor=Local density in real space.
!!
!! PARENTS
!!
!! SOURCE

elemental function k_thfermi(rhor)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rhor
 real(dp) :: k_thfermi

!Local variables-------------------------------
!scalars
 real(dp),parameter :: pisq=pi**2

! *************************************************************************

 k_thfermi = SQRT(four*k_fermi(rhor)*piinv)

end function k_thfermi
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/levi_civita_3
!! NAME
!!  levi_civita_3
!!
!! FUNCTION
 !! Return Levi-Civita tensor of rank 3
!!
!! PARENTS
!!
!! SOURCE

pure function levi_civita_3() result(ee)

!Arguments ------------------------------------
 integer :: ee(3,3,3)

! *************************************************************************

 ee = 0
 ee(1,2,3) = 1
 ee(2,3,1) = 1
 ee(3,1,2) = 1
 !
 ee(3,2,1) = -1
 ee(1,3,2) = -1
 ee(2,1,3) = -1

end function levi_civita_3
!!***

!!****f* m_special_funcs/jlspline_new
!! NAME
!! jlspline_new
!!
!! FUNCTION
!! Pre-calculate the j_v(y) for recip_ylm on regular grid
!!     NOTE: spherical Bessel function small j!
!!
!! INPUTS
!!  nx = max number of points on grid for integral
!!  delta = space between integral arguments
!!  mlang= max angular momentum
!!
!! OUTPUT
!!  bess_spl=array of integrals
!!  bess_spl_der=array of derivatives of integrals
!!  xx=coordinates of points belonging to the grid
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      besjm,spline
!!
!! SOURCE

type(jlspline_t) function jlspline_new(nx, delta, mlang) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,mlang
 real(dp),intent(in) :: delta

!Local variables -------------------------
!scalars
 integer :: ix,ll
 real(dp) :: yp1,ypn
!arrays
 real(dp),allocatable :: cosbessx(:),sinbessx(:)

! *********************************************************************

 if (nx < 2) then
   MSG_ERROR('need more than one point for the interpolation routines')
 end if

 new%nx = nx; new%mlang = mlang; new%delta = delta; new%maxarg = (nx-1) * delta
 ABI_MALLOC(new%xx, (nx))
 ABI_MALLOC(new%bess_spl, (nx, mlang))
 ABI_MALLOC(new%bess_spl_der, (nx, mlang))

 !-----------------------------------------------------------------
 !Bessel function into array
 !-----------------------------------------------------------------
 ! integration grid is nfiner times finer than the interpolation grid
 ABI_MALLOC(sinbessx, (nx))
 ABI_MALLOC(cosbessx, (nx))

 ! could be done by chain rule for cos sin (is it worth it?) but
 ! precision problems as numerical errors are propagated.
 do ix=1,nx
   new%xx(ix) = (ix-1) * delta
   sinbessx(ix) = sin(new%xx(ix))
   cosbessx(ix) = cos(new%xx(ix))
 end do

 ! fill bess_spl array
 do ll=0,mlang-1
   call besjm(one,new%bess_spl(:,ll+1),cosbessx,ll,nx,sinbessx,new%xx)

   ! call spline to get 2nd derivative (reuse in splint later)
   yp1 = zero; ypn = zero
   call spline(new%xx, new%bess_spl(:,ll+1), nx, yp1, ypn, new%bess_spl_der(:,ll+1))
 end do

!write(std_out,*) ' bess funct  0   1   2   3   4'
!do ix=1,nx
!write(std_out,*) xx(ix), (new%bess_spl(ix,ll),ll=1,mlang)
!end do

 ABI_FREE(sinbessx)
 ABI_FREE(cosbessx)

end function jlspline_new
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/jlspline_free
!! NAME
!! jlspline_free
!!
!! FUNCTION
!!  deallocate memory
!!
!! PARENTS
!!      m_cut3d,m_epjdos
!!
!! CHILDREN
!!
!! SOURCE

subroutine jlspline_free(jlspl)

!Arguments ------------------------------------
 type(jlspline_t),intent(inout) :: jlspl

! *********************************************************************

 if (allocated(jlspl%xx)) then
   ABI_FREE(jlspl%xx)
 end if
 if (allocated(jlspl%bess_spl)) then
   ABI_FREE(jlspl%bess_spl)
 end if
 if (allocated(jlspl%bess_spl_der)) then
   ABI_FREE(jlspl%bess_spl_der)
 end if

end subroutine jlspline_free
!!***

!----------------------------------------------------------------------

!!****f* m_special_funcs/jlspline_integral
!! NAME
!! jlspline_integral
!!
!! INPUTS
!!
!! OUTPUT
!!
!! FUNCTION
!!
!! PARENTS
!!
!! SOURCE

real(dp) function jlspline_integral(jlspl, il, qq, powr, nr, rcut)  result(res)

!Arguments ------------------------------------
 integer,intent(in) :: il,nr,powr
 real(dp),intent(in) :: qq, rcut
 type(jlspline_t),intent(in) :: jlspl

!Local variables ---------------------------------------
 integer :: ierr
 real(dp) :: step
!arrays
 real(dp):: xfit(nr),yfit(nr),rr(nr)
! *********************************************************************

 step = rcut / (nr - 1)
 rr = arth(zero, step, nr)
 xfit = qq * rr
 call splint(jlspl%nx, jlspl%xx, jlspl%bess_spl(:,il), jlspl%bess_spl_der(:,il), nr, xfit, yfit, ierr=ierr)

 if (ierr /= 0) then
   write(std_out,*)"qq, rcut, qq*rcut, maxarg", qq, rcut, qq*rcut, jlspl%maxarg
   write(std_out,*)"x[0], x[-1]",jlspl%xx(1),jlspl%xx(jlspl%nx)
   write(std_out,*)"minval xfit: ",minval(xfit)
   write(std_out,*)"maxval xfit: ",maxval(xfit)
   MSG_ERROR("splint returned ierr != 0")
 end if

 if (powr /= 1) yfit = yfit * (rr ** powr)
 res = simpson(step, yfit)

end function jlspline_integral
!!***

!!****f* m_special_funcs/gspline_new
!! NAME
!!  gspline_new
!!
!! FUNCTION
!!  Build object to spline the gaussian approximant and its primitive.
!!
!! INPUTS
!!  sigma=Broadening parameter.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type (gspline_t) function gspline_new(sigma) result(new)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: sigma

!Local variables ------------------------------
!scalars
 integer :: ii
 real(dp) :: ybcbeg, ybcend
! *************************************************************************

 new%nspline = 5 * 1024; new%sigma = sigma
 ABI_CHECK(sigma > zero, sjoin("invalid sigma:", ftoa(sigma)))
 new%xmin = zero
 new%xmax = sigma * sqrt(-log(sigma * sqrt(pi) * tol12)) ! gauss(xmax) = tol12
 new%step = (new%xmax - new%xmin) / (new%nspline - 1)
 new%stepm1 = one / new%step; new%step2div6 = new%step**2 / six

 ABI_MALLOC(new%xvals, (new%nspline))
 do ii=1,new%nspline
   new%xvals(ii) = new%xmin + (ii-1) * new%step
 end do
 new%xmax = new%xvals(new%nspline)

 ! Spline the gaussian approximant.
 ABI_MALLOC(new%svals, (new%nspline, 4))
 new%svals(:, 1) = gaussian(new%xvals, sigma)
 ybcbeg = - (two * new%xmin / sigma**2) * new%svals(1,1)
 ybcend = - (two * new%xmax / sigma**2) * new%svals(new%nspline,1)
 call spline(new%xvals, new%svals(:,1), new%nspline, ybcbeg, ybcend, new%svals(:,2))

 ! Spline the primitive: 1/2 [1 + erf(x/sigma)]
 new%svals(:, 3) = half * (one + abi_derf(new%xvals / new%sigma))
 call spline(new%xvals, new%svals(:,3), new%nspline, new%svals(1,1), new%svals(new%nspline, 1), new%svals(:,4))
 !do ii=1,new%nspline; write(98,*)new%xvals(ii),new%svals(ii,3),new%svals(ii,4); end do

end function gspline_new
!!***

!!****f* m_special_funcs/gspline_eval
!! NAME
!!  gspline_eval
!!
!! FUNCTION
!!  Evaluate the gaussian approximant and its primitive at (xmesh - x0)
!!
!! INPUTS
!!  self<gspline_t>=Object used to spline the gaussian approximant
!!  x0=Shift to be given to xmesh
!!  nx=Number of points in input mesh.
!!  xmesh(nx)=Frequency points (not necessarly linear).
!!
!! OUTPUT
!!  weights(nx,2)=First slice contains the gaussian approximant on xmesh.
!!   The second slice stores the primitive.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine gspline_eval(self, x0, nx, xmesh, weights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx
 real(dp),intent(in) :: x0
 type(gspline_t),intent(in) :: self
!arrays
 real(dp),intent(in) :: xmesh(nx)
 real(dp),intent(out) :: weights(nx,2)

!Local variables ------------------------------
!scalars
 integer :: ix,jspl
 real(dp) :: xx,absx,aa,bb,cc,dd
 logical :: isneg
 !real(dp) :: int_values(nx)

! *************************************************************************

 do ix=1,nx
   xx = xmesh(ix) - x0; absx = abs(xx); isneg = xx < zero
   if (absx >= self%xmax) then
     ! Region in which gauss(x) is negligible.
     weights(ix,1) = zero
     if (isneg) then
       weights(ix,2) = zero
     else
       weights(ix,2) = one
     end if
   else
     ! Spline functions at |x| and recover the value at x:
     ! g(x) = g(-x); G(-x) = 1 - G(x)
     jspl = 1 + int((absx - self%xmin) * self%stepm1); dd = absx - self%xvals(jspl)
     bb = dd * self%stepm1
     aa = one - bb
     cc = aa*(aa**2-one) * self%step2div6
     dd = bb*(bb**2-one) * self%step2div6

     weights(ix,1) = aa*self%svals(jspl,1) + bb*self%svals(jspl+1,1) + cc*self%svals(jspl,2) + dd*self%svals(jspl+1,2)
     weights(ix,2) = aa*self%svals(jspl,3) + bb*self%svals(jspl+1,3) + cc*self%svals(jspl,4) + dd*self%svals(jspl+1,4)
     if (isneg) weights(ix,2) = one - weights(ix,2)
   end if
 end do

 !call simpson_int(nx,xmesh(2) - xmesh(1),weights(:,1),int_values)
 !do ix=1,nx
 !  write(99,*)xmesh(ix), weights(ix,1), gaussian(xx, self%sigma), weights(ix,2), int_values(ix)
 !end do

end subroutine gspline_eval
!!***

!!****f* m_special_funcs/gspline_free
!! NAME
!!  gspline_free
!!
!! FUNCTION
!!  Free dynamic memory
!!
!! INPUTS
!!  self<gspline_t>=Object used to spline the gaussian approximant
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gspline_free(self)

!Arguments ------------------------------------
!scalars
 type(gspline_t),intent(inout) :: self

! *************************************************************************

 if (allocated(self%xvals)) then
   ABI_FREE(self%xvals)
 end if
 if (allocated(self%svals)) then
   ABI_FREE(self%svals)
 end if

end subroutine gspline_free
!!***

end module m_special_funcs
!!***

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
!! Copyright (C) 2008-2016 ABINIT group (MG,MT,FB,XG,FJ,NH,GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_special_funcs

 use defs_basis
 use m_profiling_abi
 use m_errors

 implicit none

 private

 public :: clp               ! x-1, if x>1/2, x+1, if x<-1/2
 public :: factorial         ! Calculates N! returning a real.
 public :: permutations      ! Returns N!/(N-k) if N>=0 and N-k>0 else 0.
 public :: binomcoeff        ! Binominal coefficient n!/(n-k)!
 public :: laguerre          ! Laguerre Polynomial(x,n,a). 
 public :: RadFnH            ! Atomic radial function(r,n,l,Z).
 public :: iradfnh           ! Norm of atomic radial function(a,b,n,l,Z).  
 public :: dirac_delta       ! Approximate Dirac delta with normal distribution.
 public :: gaussian          ! Normalized Gaussian distribution.
 public :: abi_derf          ! Evaluates the error function in real(dp).
 public :: abi_derfc         ! Evaluates the complementary error function in real(dp).
 public :: phim              ! Computes Phi_m[theta]=Sqrt[2] cos[m theta],      if m>0
                             !                       Sqrt[2] sin[Abs(m) theta], if m<0
                             !                       1                        , if m=0
 public :: k_fermi           ! Fermi wave vector corresponding to the local value of the real space density rhor.
 public :: k_thfermi         ! Thomas-Fermi wave vector corresponding to the local value of the real space density rhor

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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clp'
!End of the abilint section

 implicit none

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
!!      setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

elemental function factorial(nn)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'factorial'
!End of the abilint section

 implicit none

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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'permutations'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'binomcoeff'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'laguerre'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'RadFnH'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'IRadFnH'
!End of the abilint section

 implicit none

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

!!****f* m_special_funcs/dirac_delta
!! NAME
!! dirac_delta
!!
!! FUNCTION
!!  Approximate Dirac delta with normal distribution.
!!    delta(x,sigma) = 1/(sigma sqrt(pi)) e^{-x**2/(sigma**2)}
!!
!! INPUTS
!!   arg=Argument of the approximated Delta.
!!   sigma=Broadening factor.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function dirac_delta(arg,sigma)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dirac_delta'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 real(dp),intent(in) :: arg,sigma
 real(dp) :: dirac_delta

!Local variables ---------------------------------------
!scalars
 real(dp) :: xx

! *********************************************************************

 xx=arg/sigma
 dirac_delta = exp(-xx*xx) / (sigma*sqrt(pi))

end function dirac_delta
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

elemental function gaussian(arg,sigma)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaussian'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 real(dp),intent(in) :: arg,sigma
 real(dp) :: gaussian

!Local variables ---------------------------------------
!scalars
 real(dp) :: xx

! *********************************************************************

 xx=arg/(sqrt2*sigma)
 gaussian = exp(-xx*xx) / (sigma*sqrt(two_pi))

end function gaussian
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_derf'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_derfc'
!End of the abilint section

 implicit none

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

!!****f* ABINIT/phim
!! NAME
!! phim
!!
!! FUNCTION
!! Computes Phi_m[theta]=Sqrt[2] cos[m theta],      if m>0
!!                       Sqrt[2] sin[Abs(m) theta], if m<0
!!                       1                        , if m=0
!!
!! INPUTS
!!  costeta= cos(theta)  (theta= input angle)
!!  mm = index m
!!  sinteta= sin(theta)  (theta= input angle)
!!
!! OUTPUT
!!  phim= Phi_m(theta) (see above)
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!     setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

pure function phim(costheta,sintheta,mm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phim'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mm
 real(dp) :: phim
 real(dp),intent(in) :: costheta,sintheta

! *********************************************************************

 if (mm==0)  phim=one
 if (mm==1)  phim=sqrt2*costheta
 if (mm==-1) phim=sqrt2*sintheta
 if (mm==2)  phim=sqrt2*(costheta*costheta-sintheta*sintheta)
 if (mm==-2) phim=sqrt2*two*sintheta*costheta
 if (mm==3)  phim=sqrt2*&
& (costheta*(costheta*costheta-sintheta*sintheta)&
& -sintheta*two*sintheta*costheta)
 if (mm==-3) phim=sqrt2*&
& (sintheta*(costheta*costheta-sintheta*sintheta)&
& +costheta*two*sintheta*costheta)

 end function phim
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_fermi'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_thfermi'
!End of the abilint section

 implicit none

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

END MODULE m_special_funcs
!!***

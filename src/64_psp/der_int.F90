!{\src2tex{textfont=tt}}
!!****f* ABINIT/der_int
!! NAME
!! der_int
!!
!! FUNCTION
!! Given input function f(i) on Teter radial grid, and grid spacing
!! dr(i), compute function derivative df/dr on points from 0 to n.
!! Integrate function f(i) on grid r(i) from r(0) to r(nlast).
!! Note that array dimensions start at 0.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  f(0 to nlast)=function values on grid
!!  r(0 to nlast)=radial grid points
!!  dr(0 to nlast)=INVERSE of spacing on grid
!!  nlast=radial grid point for upper limit
!!
!! OUTPUT
!!  df(0 to n)=derivative $ \frac{df}{dr}$ on grid
!!  smf= $ \int_{r(0)}^{r(nlast)} f(r) dr $.
!!
!! PARENTS
!!      psp1lo,psp1nl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine der_int(ff,df,rr,dr,nlast,smf)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'der_int'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!nmax sets standard number of grid points ! SHOULD BE REMOVED
!scalars
 integer,parameter :: nmax=2000
 integer,intent(in) :: nlast
 real(dp),intent(out) :: smf
!no_abirules
!Note that dimension here starts at 0
 real(dp), intent(in) :: dr(0:nmax),ff(0:nmax),rr(0:nmax)
 real(dp), intent(out) :: df(0:nmax)

!Local variables-------------------------------
!scalars
 integer :: jj
 real(dp),parameter :: div12=1.d0/12.d0
 real(dp) :: hh
 character(len=500) :: message

! *************************************************************************

!Check that nlast lie within 0 to nmax
 if (nlast<0.or.nlast>nmax) then
   write(message, '(a,i12,a,i12)' )&
&   ' nlast=',nlast,' lies outside range [0,nmax] with dimension nmax=',nmax
   MSG_BUG(message)
 end if

!Compute derivatives at lower end, near r=0
 df(0)=-25.d0/12.d0*ff(0)+4.d0*ff(1)-3.d0*ff(2)+4.d0/3.d0*ff(3)&
& -1.d0/4.d0*ff(4)
 df(1)=-1.d0/4.d0*ff(0)-5.d0/6.d0*ff(1)+3.d0/2.d0*ff(2)&
& -1.d0/2.d0*ff(3)+1.d0/12.d0*ff(4)

!Run over range from just past r=0 to near r(n), using central differences
 do jj=2,nlast-2
   df(jj)=(ff(jj-2)-8.d0*(ff(jj-1)-ff(jj+1))-ff(jj+2))*div12
 end do

!Compute derivative at upper end of range
 if (nlast < 4) then 
   message = ' der_int: ff does not have enough elements. nlast is too low'
   MSG_ERROR(message)
 end if

 df(nlast-1)=-1.d0/12.d0*ff(nlast-4)&
& +1.d0/2.d0*ff(nlast-3)&
& -3.d0/2.d0*ff(nlast-2)&
& +5.d0/6.d0*ff(nlast-1)&
& +1.d0/4.d0*ff(nlast)
 df(nlast)=1.d0/4.d0*ff(nlast-4)&
& -4.d0/3.d0*ff(nlast-3)&
& +3.d0*ff(nlast-2)&
& -4.d0*ff(nlast-1)&
& +25.d0/12.d0*ff(nlast)

!Apply correct normalization over full range
 do jj=0,nlast
   df(jj)=df(jj)*dr(jj)
 end do

 smf=0.d0
 do jj=0,nlast-1
   hh=rr(jj+1)-rr(jj)
   smf=smf+hh*(6.d0*(ff(jj)+ff(jj+1))+hh*(df(jj)-df(jj+1)))
 end do
 smf=smf/12.d0

end subroutine der_int
!!***

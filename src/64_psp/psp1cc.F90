!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp1cc
!! NAME
!! psp1cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 1 and 5 of pseudopotentials.
!! WARNING : the fifth derivate is actually set to zero
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG, DCA, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fchrg=magnitude of the core charge correction
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! NOTES
!! This is a revised expression for core density (5 Nov 1992) :
!! density(r)=fchrg*gg(xx)
!! with
!! $ gg(xx)=(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4 xx^2)(1-xx^2)})^2 $
!! and
!! $ xx=\frac{r}{rchrg}=\frac{r}{xcccrc/3.0d0}=3*\frac{r}{xcccrc}=3*yy $
!!
!! Code for gg(xx), gp(xx), and gpp(xx) has been tested by numerical
!! derivatives--looks ok. gpp(x) should still be rewritten.
!! The argument of xccc1d is assumed to be normalized, and to vary
!! from yy=0 to 1 (from r=0 to r=xcccrc, or from xx=0 to 3)
!! Thus :
!!{{\ \begin{equation}
!! xccc1d(yy)=fchrg*[\frac{\sin(2*\pi*(3yy))}
!! {(6*\pi*(3yy))(1-4*(3yy)^2)(1-(3yy)^2)}]^2
!!\end{equation} }}
!!
!! WARNINGS
!! Warning : the fifth derivative is not yet delivered.
!!
!! PARENTS
!!      psp1in,psp5in
!!
!! CHILDREN
!!      gg1cc,gp1cc,gpp1cc,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp1cc(fchrg,n1xccc,xccc1d)

 use defs_basis
 use m_splines
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp1cc'
 use interfaces_64_psp, except_this_one => psp1cc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1xccc
 real(dp),intent(in) :: fchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1xccc,ider
 real(dp) :: der1,dern,factor,gg1cc_xx,gp1cc_xx,gpp1cc_xx,xx
 character(len=500) :: message
!arrays
 real(dp),allocatable :: ff(:),ff2(:),work(:),yy(:)

! *************************************************************************

 ABI_ALLOCATE(ff,(n1xccc))
 ABI_ALLOCATE(ff2,(n1xccc))
 ABI_ALLOCATE(work,(n1xccc))
 ABI_ALLOCATE(yy,(n1xccc))

 if(n1xccc > 1)then
   factor=one/dble(n1xccc-1)
   do i1xccc=1,n1xccc
     yy(i1xccc)=(i1xccc-1)*factor
   end do
 else
   write(message, '(a,i0)' )' n1xccc should larger than 1, while it is n1xccc=',n1xccc
   MSG_BUG(message)
 end if

!Initialization, to avoid some problem with some compilers
 xccc1d(1,:)=zero ; xccc1d(n1xccc,:)=zero

!Take care of each derivative separately
 do ider=0,2

   if(ider==0)then
!    Generate spline fitting for the function gg
     do i1xccc=1,n1xccc
       xx=three*yy(i1xccc)
       call gg1cc(gg1cc_xx,xx)
       ff(i1xccc)=fchrg*gg1cc_xx
     end do
!    Complete with derivatives at end points
     der1=zero
     call gp1cc(gp1cc_xx,three)
     dern=three*fchrg*gp1cc_xx
   else if(ider==1)then
!    Generate spline fitting for the function gp
     do i1xccc=1,n1xccc
       xx=three*yy(i1xccc)
       call gp1cc(gp1cc_xx,xx)
       ff(i1xccc)=three*fchrg*gp1cc_xx
     end do
!    Complete with derivatives at end points, already estimated
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   else if(ider==2)then
!    Generate spline fitting for the function gpp
!    (note : the function gpp has already been estimated, for the spline
!    fitting of the function gg, but it is replaced here by the more
!    accurate analytic derivative)
     do i1xccc=1,n1xccc
       xx=three*yy(i1xccc)
       call gpp1cc(gpp1cc_xx,xx)
       ff(i1xccc)=9.0_dp*fchrg*gpp1cc_xx
     end do
!    Complete with derivatives of end points
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   end if

!  Produce second derivative numerically, for use with splines
   call spline(yy,ff,n1xccc,der1,dern,ff2)
   xccc1d(:,ider+1)=ff(:)
   xccc1d(:,ider+3)=ff2(:)

 end do

 xccc1d(:,6)=zero

!DEBUG
!write(std_out,*)' psp1cc : output of core charge density and derivatives '
!write(std_out,*)'   yy          gg           gp  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) yy(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   yy          gpp          gg2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) yy(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   yy          gp2          gpp2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) yy(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(yy)

end subroutine psp1cc
!!***

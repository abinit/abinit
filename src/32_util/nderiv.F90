!{\src2tex{textfont=tt}}
!!****f* ABINIT/nderiv
!! NAME
!! nderiv
!!
!! FUNCTION
!! Given an input function y(x) on a regular grid,
!! compute its first or second derivative.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  hh= radial step
!!  ndim= radial mesh size
!!  yy(ndim)= input function
!!  norder= order of derivation (1 or 2)
!!
!! OUTPUT
!!  zz(ndim)= first or second derivative of y
!!
!! PARENTS
!!      upf2abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nderiv(hh,yy,zz,ndim,norder)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nderiv'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndim,norder
 real(dp),intent(in) :: hh
!arrays
 real(dp),intent(in) :: yy(ndim)
 real(dp),intent(out) :: zz(ndim)

!Local variables ---------------------------------------
!scalars
 integer :: ier,ii
 real(dp) :: aa,bb,cc,h1,y1

! *********************************************************************

!Initialization (common to 1st and 2nd derivative)
 h1=one/(12.d0*hh)
 y1=yy(ndim-4)

!FIRST DERIVATIVE
!================
 if (norder==1) then

!  Prepare differentiation loop
   bb=h1*(-25.d0*yy(1)+48.d0*yy(2)-36.d0*yy(3)+16.d0*yy(4)-3.d0*yy(5))
   cc=h1*(-3.d0*yy(1)-10.d0*yy(2)+18.d0*yy(3)-6.d0*yy(4)+yy(5))
!  Start differentiation loop
   do ii=5,ndim
     aa=bb;bb=cc
     cc=h1*(yy(ii-4)-yy(ii)+8.d0*(yy(ii-1)-yy(ii-3)))
     zz(ii-4)=aa
   end do
!  Normal exit
   ier=0
   aa=h1*(-y1+6.d0*yy(ndim-3)-18.d0*yy(ndim-2)+10.d0*yy(ndim-1)+3.d0*yy(ndim))
   zz(ndim)=h1*(3.d0*y1-16.d0*yy(ndim-3)+36.d0*yy(ndim-2) -48.d0*yy(ndim-1)+25.d0*yy(ndim))
   zz(ndim-1)=aa
   zz(ndim-2)=cc
   zz(ndim-3)=bb

!  SECOND DERIVATIVE
!  =================
 else
   h1=h1/hh
!  Prepare differentiation loop
   bb=h1*(35.d0*yy(1)-104.d0*yy(2)+114.d0*yy(3)-56.d0*yy(4)+11.d0*yy(5))
   cc=h1*(11.d0*yy(1)-20.d0*yy(2)+6.d0*yy(3)+4.d0*yy(4)-yy(5))
!  Start differentiation loop
   do ii=5,ndim
     aa=bb;bb=cc
     cc=h1*(-yy(ii-4)-yy(ii)+16.d0*(yy(ii-1)+yy(ii-3))-30.d0*yy(ii-2))
     zz(ii-4)=aa
   end do
!  Normal exit
   ier=0
   aa=h1*(-y1+4.d0*yy(ndim-3)+6.d0*yy(ndim-2)-20.d0*yy(ndim-1)+11.d0*yy(ndim))
   zz(ndim)=h1*(11.d0*y1-56.d0*yy(ndim-3)+114.d0*yy(ndim-2) -104.d0*yy(ndim-1)+35.d0*yy(ndim))
   zz(ndim-1)=aa
   zz(ndim-2)=cc
   zz(ndim-3)=bb

 end if !norder

end subroutine nderiv
!!***

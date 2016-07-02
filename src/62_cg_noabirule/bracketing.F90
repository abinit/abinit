!{\src2tex{textfont=tt}}
!!****f* ABINIT/bracketing
!! NAME
!! bracketing
!!
!! FUNCTION
!! bracket a minimun of a function f
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! dp_dum_vdp: the function of which the mimimum should be bracketted
!!
!!
!! OUTPUT
!! b= last member of the bracketing triplet a < x < b
!! fa,fx,fb= value of the function at dp_dum_vdp(v(:)+y*grad(:))
!!
!!
!! SIDE EFFECTS
!! v: the initial vector for the function (return unchanged)
!! grad: the direction on which the bracketting is to be performed (return unchanged)
!! a,x: two members of the bracketing triplet (see b)
!!
!! PARENTS
!!      linmin
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine bracketing (nv1,nv2,dp_dum_v2dp,v,grad,a,x,b,fa,fx,fb)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bracketing'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 include "dummy_functions.inc"
!scalars
 integer,intent(in) :: nv1,nv2
 real(dp),intent(inout) :: a,x
 real(dp),intent(out) :: b,fa,fb,fx
!arrays
 real(dp),intent(inout) :: grad(nv1,nv2),v(nv1,nv2)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: maglimit=10000.0_dp
 real(dp) :: c,fu,q,r,u,ulim

! *************************************************************************

 fa=dp_dum_v2dp(nv1,nv2,v(:,:)+(a*grad(:,:)))
 fx=dp_dum_v2dp(nv1,nv2,(x*grad(:,:))+v(:,:))
 if(fx > fa) then
  c=a
  a=x
  x=c
  c=fa
  fa=fx
  fx=c
 end if
 b=x+gold*(x-a)
 fb=dp_dum_v2dp(nv1,nv2,(b*grad(:,:))+v(:,:))
 do
  if (fx <= fb) return
  r=(x-a)*(fx-fb)
  q=(x-b)*(fx-fa)
  u=x-((x-b)*q-(x-a)*r)/(two*sign(max(abs(q-r),smallest_real),q-r))
  ulim=x+maglimit*(b-x)
  if((x-u)*(u-b) > zero) then
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
   if(fu < fb) then
    a=x
    fa=fx
    x=u
    fx=fu
    return
   else if (fx < fu) then
    b=u
    fb=fu
    return
   end if
   u=b+gold*(b-x)
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  else if((b-u)*(u-ulim) > zero) then
   fu=dp_dum_v2dp(nv1,nv2,u*grad(:,:)+v(:,:))
   if(fu<fb) then
    x=b
    b=u
    u=b+gold*(b-x)
    fx=fb
    fb=fu
    fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
   end if
  else if((u-ulim)*(ulim-b) >= zero) then
   u=ulim
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  else
   u=b+gold*(b-x)
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  end if
  a=x
  x=b
  b=u
  fa=fx
  fx=fb
  fb=fu
 end do

end subroutine bracketing
!!***

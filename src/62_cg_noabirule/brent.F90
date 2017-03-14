!{\src2tex{textfont=tt}}
!!****f* ABINIT/brent
!! NAME
!! brent
!!
!! FUNCTION
!! minimizes a function along a line
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! dp_dum_vdp: function  to be minimized (return a dp from a vector of dp)
!! vdp_dum_vdp: derivative of the function (return a vector of dp from a vector of dp)
!! itmax: number of iterations allowed
!! tol: tolerance on error. It depend on the precision of the numbers
!! (usualy chosen as sqrt(max precision available with your floating point reresentation))
!! ax,xx,bx: a bracketing triplet around the minimum to be find
!! OUTPUT
!! xmin: value such that dp_dum_vdp(v(:)+xmin*grad(:)) is minimum
!! brent:  dp_dum_vdp(v(:)+xmin*grad(:))
!!
!! SIDE EFFECTS
!! grad(:): direction along which the minimization is performed
!! v(:): starting and ending point of the minimization
!!
!! PARENTS
!! linmin
!!
!! CHILDREN
!! dotproduct
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function brent(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,itmax,v,grad,ax,xx,bx,tol,xmin)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'brent'
 use interfaces_62_cg_noabirule, except_this_one => brent
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 include "dummy_functions.inc"
!scalars
 integer,intent(in) :: itmax,nv1,nv2
 real(dp) :: brent
 real(dp),intent(in) :: ax,bx,tol,xx
 real(dp),intent(out) :: xmin
!arrays
 real(dp),intent(inout) :: grad(nv1,nv2),v(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: iter
 real(dp) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,vv,w
 real(dp) :: x,xm,zeps
 logical :: ok1,ok2,ok3,ok4

!************************************************************************
 zeps=epsilon(ax*real(1e-2,dp))
 a=min(ax,bx)
 b=max(ax,bx)
 vv=xx
 w=xx
 x=xx
 e=zero
 fx=dp_dum_v2dp(nv1,nv2,x*grad(:,:)+v(:,:))
 fv=fx
 fw=fx
!the function sub_dum_dp_v2dp_v2dp must do the equivalent of
!v(:,:)=v(:,:)+(grad(:,:)*x)
!but for instance renormilizing the density if brent is used on a density...
!vp(:,:) = v(:,:)
!sub_dum_dp_v2dp_v2dp(x,grad(:,:),vp(:,:)
!dx=dotproduct(v2dp_dum_v2dp(vp(:,:)),grad(:,:))
 dx=dotproduct(nv1,nv2,v2dp_dum_v2dp(nv1,nv2,v(:,:)+x*grad(:,:)),grad(:,:))
 dv=dx
 dw=dx
 do iter=1,itmax
  xm=half*(a+b)
  tol1=tol*abs(x)+zeps
  tol2=two*tol1
  if(abs(x-xm) <= (tol2-half*(b-a))) then
   exit
  end if
  if(abs(e) > tol1) then
   d1=two*(b-a)
   d2=d1
   if(dw /= dx) d1=(w-x)*dx/(dx-dw)
   if(dv /= dx) d2=(vv-x)*dx/(dx-dv)
   u1=x+d1
   u2=x+d2
   ok1=((a-u1)*(u1-b)>zero).and.(dx*d1<=zero)
   ok2=((a-u2)*(u2-b)>zero).and.(dx*d2<=zero)
   olde=e
   e=d
   if(ok1.or.ok2) then
    if(ok1.and.ok2) then
     d=merge(d1,d2,abs(d1)<abs(d2))
    else
     d=merge(d1,d2,ok1)
    end if
    if(abs(d)<=abs(half*olde)) then
     u=x+d
     if(((u-a)<tol2).or.((b-u)<tol2)) d=sign(tol1,xm-x)
    else
     e=merge(a,b,dx>=zero)-x
     d=half*e
    end if
   else
    e=merge(a,b,dx>=zero)-x
    d=half*e
   end if
  else
   e=merge(a,b,dx>=zero)-x
   d=half*e
  end if

  if(abs(d) >=tol1)then
   u=x+d
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  else
   u=x+sign(tol1,d)
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
   if(fu>fx) then
    exit
   end if
  end if
  du=dotproduct(nv1,nv2,v2dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:)),grad(:,:))
  if(fu<=fx)then
   if(u>=x)then
    a=x
   else
    b=x
   end if
   vv=w
   fv=fw
   dv=dw
   w=x
   fw=fx
   dw=dx
   x=u
   dx=du
   fx=fu
  else
   if(u<x) then
    a=u
   else
    b=u
   end if
   ok3=(w==x).or.(fu.le.fw)
   ok4=(vv==w).or.(vv==x).or.(fu.lt.fv)
   if(ok3) then
    vv=w
    fv=fw
    dv=dw
    w=u
    fw=fu
    dw=du
   else if( ok4 ) then
    vv=u
    fv=fu
    dv=du
   end if
  end if
 end do
 xmin=x
!the function sub_dum_dp_v2dp_v2dp must do the equivalent of
!v(:,:)=v(:,:)+(grad(:,:)*x)
!but for instance renormilizing the density if brent is used on a density...
 call sub_dum_dp_v2dp_v2dp(nv1,nv2,x,grad(:,:),v(:,:))
 brent=fx
end function brent
!!***

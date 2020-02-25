!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_integrals
!! NAME
!!  m_integrals
!!
!! FUNCTION
!!  Helper functions to compute integrals
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (Camilo Espejo)
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

module m_integrals

 use defs_basis
 use m_errors
 use m_abicore

 use m_numeric_tools,   only : simpson_int

 implicit none

 private
!!***

 public :: radsintr
!!***

contains
!!***

!!****f* ABINIT/radsintr
!! NAME
!!  radsintr
!!
!! FUNCTION
!! Computes the sine transformation of a radial function F(r) to V(q).
!! Computes integrals using corrected Simpson integration on a linear grid.
!!
!! INPUTS
!!  funr(mrgrid)=F(r) on radial grid
!!  mqgrid=number of grid points in q from 0 to qmax
!!  mrgrid=number of grid points in r from 0 to rmax
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rgrid(mrgrid)=r grid values (bohr)
!!
!! OUTPUT
!!  funq(mqgrid)=\int_0^inf 4\pi\frac{\sin(2\pi r)}{2\pi r}r^2F(r)dr
!!  yq1, yqn: d/dq (F(q)) at the ends of the interval
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  This routine is a modified version of /src/65_psp/psp7lo.F90
!!
!! PARENTS
!!      m_xc_vdw,test_radsintr
!!
!! CHILDREN
!!      simpson_int
!!
!! SOURCE

subroutine radsintr(funr,funq,mqgrid,mrgrid,qgrid,rgrid,yq1,yqn)

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: mqgrid,mrgrid
 real(dp), intent(out) :: yq1,yqn
!arrays
 real(dp), intent(in)  :: funr(mrgrid),qgrid(mqgrid),rgrid(mrgrid)
 real(dp), intent(out) :: funq(mqgrid)
!Local variables-------------------------------
!scalars
 integer :: iq,ir,irmax
 real(dp) :: arg,r0tor1,r1torm,rmax,rmtoin,rstep
 logical :: begin_r0
!arrays
 real(dp),allocatable :: ff(:),intg(:),rzf(:)

! *************************************************************************

!DEBUG
!write (std_out,*) ' radsintr : enter'
!ENDDEBUG

 rstep=rgrid(2)-rgrid(1);rmax=min(20._dp,rgrid(mrgrid))
 irmax=int(tol8+rmax/rstep)+1

!Particular case of a null fuction to transform
 if (maxval(abs(funr(1:irmax)))<=1.e-20_dp) then
   funq=zero;yq1=zero;yqn=zero
   return
 end if

 ABI_ALLOCATE(ff,(mrgrid))
 ABI_ALLOCATE(intg,(mrgrid))
 ABI_ALLOCATE(rzf,(mrgrid))
 ff=zero;rzf=zero

!Is mesh beginning with r=0 ?
 begin_r0=(rgrid(1)<1.e-20_dp)

!Store r.F(r)
 do ir=1,irmax
   ff(ir)=rgrid(ir)*funr(ir) !ff is part of the integrand
 end do

!===========================================
!=== Compute v(q) for q=0 separately
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) r0tor1=(funr(1)*third)*rgrid(1)**three

!Integral from r1 to rmax
 do ir=1,irmax
   if (abs(ff(ir))>1.e-20_dp) then
     rzf(ir)=rgrid(ir)*ff(ir) !sin(2\pi*q*r)/(2*\pi*q*r)-->1 for q=0
   end if                     !so the integrand is 4*\pi*r^2*F(r)
 end do
 call simpson_int(mrgrid,rstep,rzf,intg)
 r1torm=intg(mrgrid)

!Integral from rmax to infinity
!This part is neglected... might be improved.
 rmtoin=zero

!Sum of the three parts

 funq(1)=four_pi*(r0tor1+r1torm+rmtoin)

!===========================================
!=== Compute v(q) for other q''s
!===========================================

!Loop over q values
 do iq=2,mqgrid
   arg=two_pi*qgrid(iq)

!  Integral from 0 to r1 (only if r1<>0)
   r0tor1=zero
   if (.not.begin_r0) r0tor1=(funr(1)/arg) &
&   *(sin(arg*rgrid(1))/arg-rgrid(1)*cos(arg*rgrid(1)))

!  Integral from r1 to rmax
   rzf=zero
   do ir=1,irmax
     if (abs(ff(ir))>1.e-20_dp) rzf(ir)=sin(arg*rgrid(ir))*ff(ir)
   end do
   call simpson_int(mrgrid,rstep,rzf,intg)
   r1torm=intg(mrgrid)

!  Integral from rmax to infinity
!  This part is neglected... might be improved.
   rmtoin=zero

!  Store F(q)
   funq(iq) = two/qgrid(iq)*(r0tor1+r1torm+rmtoin)
 end do


!===========================================
!=== Compute derivatives of F(q)
!=== at ends of interval
!===========================================

!yq(0)=zero
 yq1=zero

!yp(qmax)=$
 arg=two_pi*qgrid(mqgrid)

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero
 if (.not.begin_r0) r0tor1=(funr(1)/arg) &
& *(sin(arg*rgrid(1))*(rgrid(1)**2-three/arg**2) &
& +three*rgrid(1)*cos(arg*rgrid(1))/arg)

!Integral from r1 to rmax
 do ir=1,irmax
   if (abs(ff(ir))>1.e-20_dp) rzf(ir)=(rgrid(ir)*cos(arg*rgrid(ir)) &
&   -sin(arg*rgrid(ir))/arg)*ff(ir)
 end do
 call simpson_int(mrgrid,rstep,rzf,intg)
 r1torm=intg(mrgrid)

!Integral from rmax to infinity
!This part is neglected... might be improved.
 rmtoin=zero

!Sum of the three parts
 yqn=(four*pi/qgrid(mqgrid))*(r0tor1+r1torm+rmtoin)

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(intg)
 ABI_DEALLOCATE(rzf)

end subroutine radsintr
!!***

end module m_integrals
!!***

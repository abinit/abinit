!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp2nl
!! NAME
!! psp2nl
!!
!! FUNCTION
!! Goedecker-Teter-Hutter nonlocal pseudopotential (from preprint of 1996).
!! Uses Gaussians for fully nonlocal form, analytic expressions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  h1p=factor defining strength of 1st projector for l=1 channel
!!  h1s=factor defining strength of 1st projector for l=0 channel
!!  h2s=factor defining strength of 2nd projector for l=0 channel
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mqgrid=number of grid points for qgrid
!!  qgrid(mqgrid)=array of |G| values
!!  rrp=core radius for p channel (bohr)
!!  rrs=core radius for s channel (bohr)
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!   and each projector
!!
!! PARENTS
!!      psp2in
!!
!! CHILDREN
!!      spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp2nl(ekb,ffspl,h1p,h1s,h2s,lnmax,mqgrid,qgrid,rrp,rrs)

 use defs_basis
 use m_splines
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2nl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lnmax,mqgrid
 real(dp),intent(in) :: h1p,h1s,h2s,rrp,rrs
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(inout) :: ekb(lnmax),ffspl(mqgrid,2,lnmax) !vz_i

!Local variables-------------------------------
!scalars
 integer :: iln,iqgrid
 real(dp) :: qmax,yp1,ypn
!arrays
 real(dp),allocatable :: work(:)

! *************************************************************************

 ABI_ALLOCATE(work,(mqgrid))

!Kleinman-Bylander energies ekb were set to zero in calling program

!Compute KB energies
 iln=0
 if (abs(h1s)>1.d-12) then
   iln=iln+1
   ekb(iln)=h1s*32.d0*rrs**3*(pi**(2.5d0)/(4.d0*pi)**2)
 end if
 if (abs(h2s)>1.d-12) then
   iln=iln+1
   ekb(iln) =h2s*(128.d0/15.d0)*rrs**3*(pi**(2.5d0)/(4.d0*pi)**2)
 end if
 if (abs(h1p)>1.d-12) then
   iln=iln+1
   ekb(iln)=h1p*(64.d0/3.d0)*rrp**5*(pi**(2.5d0)/(4.d0*pi)**2)
 end if

!Compute KB form factor
 iln=0

!l=0 first projector
 if (abs(h1s)>1.d-12) then
   iln=iln+1
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,iln)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2)
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1=0.d0
   qmax=qgrid(mqgrid)
   ypn=-4.d0*pi**2*qmax*rrs**2*exp(-0.5d0*(two_pi*qmax*rrs)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))
!  else
!  or else put first projector nonlocal correction at l=0 to 0
!  ffspl(:,:,iln)=0.0d0
 end if

!l=0 second projector
 if (abs(h2s)>1.d-12) then
   iln=iln+1
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,iln)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2) * &
&     (3.d0-(two_pi*qgrid(iqgrid)*rrs)**2)
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1=0.d0
   qmax=qgrid(mqgrid)
   ypn=4.d0*pi**2*qmax*rrs**2*exp(-0.5d0*(two_pi*qmax*rrs)**2) * &
&   (-5.d0+(two_pi*qmax*rrs)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))
!  else if(mproj>=2)then
!  or else put second projector nonlocal correction at l=0 to 0
!  ffspl(:,:,iln)=0.0d0
 end if

!l=1 first projector
 if (abs(h1p)>1.d-12) then
   iln=iln+1
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,iln)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * &
&     (two_pi*qgrid(iqgrid))
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1=two_pi
   qmax=qgrid(mqgrid)
   ypn=-two_pi*((two_pi*qmax*rrp)**2-1.d0) * exp(-0.5d0*(two_pi*qmax*rrp)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))
!  else if(mpsang>=2)then
!  or else put first projector l=1 nonlocal correction to 0
!  ffspl(:,:,iln)=0.0d0
 end if

 ABI_DEALLOCATE(work)

end subroutine psp2nl
!!***

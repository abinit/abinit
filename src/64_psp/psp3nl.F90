!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp3nl
!! NAME
!! psp3nl
!!
!! FUNCTION
!! Hartwigsen-Goedecker-Hutter nonlocal pseudopotential (from preprint of 1998).
!! Uses Gaussians for fully nonlocal form, analytic expressions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, FRD, XG, GMR, PT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  h11s=factor defining strength of 1st projector for l=0 channel
!!  h22s=factor defining strength of 2nd projector for l=0 channel
!!  h33s=factor defining strength of 3rd projector for l=0 channel
!!  h11p=factor defining strength of 1st projector for l=1 channel
!!  h22p=factor defining strength of 2nd projector for l=1 channel
!!  h33p=factor defining strength of 2nd projector for l=1 channel
!!  h11d=factor defining strength of 1st projector for l=2 channel
!!  h22d=factor defining strength of 2nd projector for l=2 channel
!!  h33d=factor defining strength of 2nd projector for l=2 channel
!!  h11f=factor defining strength of 1st projector for l=3 channel
!!  mproj=maximum number of projectors in any channel
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for qgrid
!!  qgrid(mqgrid)=array of |G| values
!!  rrd=core radius for d channel (bohr)
!!  rrf=core radius for f channel (bohr)
!!  rrp=core radius for p channel (bohr)
!!  rrs=core radius for s channel (bohr)
!!
!! OUTPUT
!!  ekb(mpsang,mproj)=Kleinman-Bylander energies
!!  ffspl(mqgrid,2,mpssang,mproj)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projectors
!!
!! PARENTS
!!      psp3in
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp3nl(ekb,ffspl,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&
&                  h33d,h11f,mproj,mpsang,mqgrid,qgrid,rrd,rrf,rrp,rrs)

 use defs_basis
 use m_splines
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp3nl'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mproj,mpsang,mqgrid
 real(dp),intent(in) :: h11d,h11f,h11p,h11s,h22d,h22p,h22s,h33d,h33p,h33s,rrd
 real(dp),intent(in) :: rrf,rrp,rrs
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(mpsang,mproj),ffspl(mqgrid,2,mpsang,mproj)

!Local variables-------------------------------
!scalars
 integer :: info,iproj,iqgrid,ldz,mu,nproj,nu
 real(dp) :: qmax
 character(len=500) :: message
 character :: jobz,uplo
!arrays
 real(dp) :: ap(2,9),rwork1(9),work1(2,9),ww(3),yp1j(3),ypnj(3)
 real(dp),allocatable :: ppspl(:,:,:,:),uu(:,:),work(:),zz(:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' psp3nl : enter '
!stop
!ENDDEBUG

 ABI_ALLOCATE(ppspl,(mqgrid,2,mpsang,mproj))
 ABI_ALLOCATE(work,(mqgrid))

 qmax=qgrid(mqgrid)
 jobz='v'
 uplo='u'
 ekb(:,:)=0.0d0

!---------------------------------------------------------------
!Treat s channel

 nproj=0
 ap(:,:)=0.0d0
!If there is at least one s-projector
 if  ( abs(h11s) >= 1.0d-8 ) then
   nproj=1 ; ldz=1 ; ap(1,1)=h11s
 end if
 nproj=1
!If there is a second projector
 if  ( abs(h22s) >= 1.0d-8 ) then
   nproj=2 ; ldz=2 ; ap(1,3)=h22s
   ap(1,2)=-0.5d0*sqrt(0.6d0)*h22s
 end if
!If there is a third projector
 if ( abs(h33s) >= 1.0d-8 ) then
   nproj=3 ; ldz=3 ; ap(1,6)=h33s
   ap(1,4)=0.5d0*sqrt(5.d0/21.d0)*h33s
   ap(1,5)=-0.5d0*sqrt(100.d0/63.d0)*h33s
 end if

 if(nproj/=0)then

   ABI_ALLOCATE(uu,(nproj,nproj))
   ABI_ALLOCATE(zz,(2,nproj,nproj))

   if (nproj > 1) then
     call ZHPEV(jobz,uplo,nproj,ap,ww,zz,ldz,work1,rwork1,info)
     uu(:,:)=zz(1,:,:)
   else
     ww(1)=h11s
     uu(1,1)=1.0d0
   end if

!  Initialization of ekb, and spline fitting
   do iproj=1,nproj
     ekb(1,iproj)=ww(iproj)*32.d0*(rrs**3)*(pi**2.5d0)/(4.d0*pi)**2
     if(iproj==1)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,1,1)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2)
       end do
       yp1j(1)=0.d0
       ypnj(1)=-(two_pi*rrs)**2*qmax*exp(-0.5d0*(two_pi*qmax*rrs)**2)
     else if(iproj==2)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,1,2)=2.0d0/sqrt(15.0d0)     &
&         *exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2) &
&         *( 3.d0-(two_pi*qgrid(iqgrid)*rrs)**2 )
       end do
       yp1j(2)=0.0d0
       ypnj(2)=2.0d0/sqrt(15.0d0)*(two_pi*rrs)**2*qmax &
&       *exp(-0.5d0*(two_pi*qmax*rrs)**2) * (-5.d0+(two_pi*qmax*rrs)**2)
     else if(iproj==3)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,1,3)=(4.0d0/3.0d0)/sqrt(105.0d0)*&
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2) * &
&         (15.0d0-10.0d0*(two_pi*qgrid(iqgrid)*rrs)**2 + &
&         (two_pi*qgrid(iqgrid)*rrs)**4)
       end do
       yp1j(3)=0.0d0
       ypnj(3)=(4.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrs)**2) * &
&       (two_pi*rrs)**2*qmax*(-35.0d0+14d0*(two_pi*qmax*rrs)**2-(two_pi*qmax*rrs)**4)
     end if
     call spline(qgrid,ppspl(:,1,1,iproj),mqgrid,&
&     yp1j(iproj),ypnj(iproj),ppspl(:,2,1,iproj))
   end do

!  Linear combination using the eigenvectors
   ffspl(:,:,1,:)=0.0d0
   do mu=1,nproj
     do nu=1,nproj
       do iqgrid=1,mqgrid
         ffspl(iqgrid,1:2,1,mu)=ffspl(iqgrid,1:2,1,mu) &
&         +uu(nu,mu)*ppspl(iqgrid,1:2,1,nu)
       end do
     end do
   end do

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(zz)

!  End condition on nproj(/=0)
 end if

!DEBUG
!write(std_out,*)' psp3nl : after s channel '
!stop
!ENDDEBUG

!--------------------------------------------------------------------
!Now treat p channel

 nproj=0
 ap(:,:)=0.0d0
!If there is at least one projector
 if  ( abs(h11p) >= 1.0d-8 ) then
   nproj=1 ; ldz=1 ; ap(1,1)=h11p
 end if
!If there is a second projector
 if  ( abs(h22p) >= 1.0d-8 ) then
   nproj=2 ; ldz=2 ; ap(1,3)=h22p
   ap(1,2)=-0.5d0*sqrt(5.d0/7.d0)*h22p
 end if
!If there is a third projector
 if ( abs(h33p) >= 1.0d-8 ) then
   nproj=3 ; ldz=3 ; ap(1,6)=h33p
   ap(1,4)= (1.d0/6.d0)*sqrt(35.d0/11.d0)*h33p
   ap(1,5)=-(1.d0/6.d0)*(14.d0/sqrt(11.d0))*h33p
 end if

 if(nproj/=0)then

   ABI_ALLOCATE(uu,(nproj,nproj))
   ABI_ALLOCATE(zz,(2,nproj,nproj))

   if (nproj > 1) then
     call ZHPEV(jobz,uplo,nproj,ap,ww,zz,ldz,work1,rwork1,info)
     uu(:,:)=zz(1,:,:)
   else
     ww(1)=h11p
     uu(1,1)=1.0d0
   end if

!  Initialization of ekb, and spline fitting
   do iproj=1,nproj
     ekb(2,iproj)=ww(iproj)*64.d0*(rrp**5)*(pi**2.5d0)/(4.d0*pi)**2
     if(iproj==1)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,2,1)=(1.0d0/sqrt(3.0d0))* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * (two_pi*qgrid(iqgrid))
       end do
       yp1j(1)=two_pi*(1.0d0/sqrt(3.0d0))
       ypnj(1)=-two_pi*((two_pi*qmax*rrp)**2-1.d0)*exp(-0.5d0*(two_pi*qmax*rrp)**2)*&
&       (1.0d0/sqrt(3.0d0))
     else if(iproj==2)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,2,2)=(2.0d0/sqrt(105.0d0))* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * &
&         (two_pi*qgrid(iqgrid))*(5.0d0-(two_pi*qgrid(iqgrid)*rrp)**2)
       end do
       yp1j(2)=(5.0d0*two_pi)*(2.0d0/sqrt(105.0d0))
       ypnj(2)=(2.0d0/sqrt(105.0d0))*two_pi*exp(-0.5d0*(two_pi*qmax*rrp)**2)* &
&       (-8*(two_pi*qmax*rrp)**2 + (two_pi*qmax*rrp)**4 + 5.0d0)
     else if(iproj==3)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,2,3)=(4.0d0/3.0d0)/sqrt(1155d0)*&
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * &
&         (two_pi*qgrid(iqgrid))*&
&         (35.0d0-14.0d0*(two_pi*qgrid(iqgrid)*rrp)**2+(two_pi*qgrid(iqgrid)*rrp)**4)
       end do
       yp1j(3)=(35.0d0*two_pi)*(4.0d0/3.0d0)/sqrt(1155.0d0)
       ypnj(3)=(4.0d0/3.0d0)/sqrt(1155.0d0)*two_pi*exp(-0.5d0*(two_pi*qmax*rrp)**2)* &
&       (35.0d0-77.0d0*(two_pi*qmax*rrp)**2+19.0d0*(two_pi*qmax*rrp)**4 - &
&       (two_pi*qmax*rrp)**6)
     end if
     call spline(qgrid,ppspl(:,1,2,iproj),mqgrid,&
&     yp1j(iproj),ypnj(iproj),ppspl(:,2,2,iproj))
   end do

!  Linear combination using the eigenvectors
   ffspl(:,:,2,:)=0.0d0
   do mu=1,nproj
     do nu=1,nproj
       do iqgrid=1,mqgrid
         ffspl(iqgrid,1:2,2,mu)=ffspl(iqgrid,1:2,2,mu) &
&         +uu(nu,mu)*ppspl(iqgrid,1:2,2,nu)
       end do
     end do
   end do

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(zz)

!  End condition on nproj(/=0)
 end if

!DEBUG
!write(std_out,*)' psp3nl : after p channel '
!stop
!ENDDEBUG

!-----------------------------------------------------------------------
!Now treat d channel.

 nproj=0
 ap(:,:)=0.0d0
!If there is at least one projector
 if  ( abs(h11d) >= 1.0d-8 ) then
   nproj=1 ; ldz=1 ; ap(1,1)=h11d
 end if
!If there is a second projector
 if  ( abs(h22d) >= 1.0d-8 ) then
   nproj=2 ; ldz=2 ; ap(1,3)=h22d
   ap(1,2)=-0.5d0*sqrt(7.d0/9.d0)*h22d
 end if
!If there is a third projector. Warning : only two projectors are allowed.
 if ( abs(h33d) >= 1.0d-8 ) then
   write(message, '(a,a,a)' )&
&   '  only two d-projectors are allowed ',ch10,&
&   '  Action : check your pseudopotential file.'
   MSG_ERROR(message)
!  nproj=3 ; ldz=3 ; ap(1,6)=h33d
!  ap(1,4)= 0.5d0*sqrt(63.d0/143.d0)*h33d
!  ap(1,5)= -0.5d0*(18.d0/sqrt(143.d0))*h33d
 end if

 if(nproj/=0)then

   ABI_ALLOCATE(uu,(nproj,nproj))
   ABI_ALLOCATE(zz,(2,nproj,nproj))

   if (nproj > 1) then
     call ZHPEV(jobz,uplo,nproj,ap,ww,zz,ldz,work1,rwork1,info)
     uu(:,:)=zz(1,:,:)
   else
     ww(1)=h11d
     uu(1,1)=1.0d0
   end if

!  Initialization of ekb, and spline fitting
   do iproj=1,nproj
     ekb(3,iproj)=ww(iproj)*128.d0*(rrd**7)*(pi**2.5d0)/(4.d0*pi)**2
     if(iproj==1)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,3,1)=(1.0d0/sqrt(15.0d0))* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrd)**2) * (two_pi*qgrid(iqgrid))**2
       end do
       yp1j(1)=0.0d0
       ypnj(1)=(1.0d0/sqrt(15.0d0))*(two_pi**2)*&
&       exp(-0.5d0*(two_pi*qmax*rrd)**2)*qmax*(2d0-(two_pi*qmax*rrd)**2)
     else if(iproj==2)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,3,2)=(2.0d0/3.0d0)/sqrt(105.0d0)* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrd)**2) * &
&         ((two_pi*qgrid(iqgrid))**2)*(7.0d0-(two_pi*qgrid(iqgrid)*rrd)**2)
       end do
       yp1j(2)=0.0d0
       ypnj(2)=(2.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrd)**2)* &
&       qmax*(two_pi**2)*( (two_pi*qmax*rrd)**4 - 11.0d0*(two_pi*qmax*rrd)**2 + 14.0d0)
     end if
     call spline(qgrid,ppspl(:,1,3,iproj),mqgrid,&
&     yp1j(iproj),ypnj(iproj),ppspl(:,2,3,iproj))
   end do

!  Linear combination using the eigenvectors
   ffspl(:,:,3,:)=0.0d0
   do mu=1,nproj
     do nu=1,nproj
       do iqgrid=1,mqgrid
         ffspl(iqgrid,1:2,3,mu)=ffspl(iqgrid,1:2,3,mu) &
&         +uu(nu,mu)*ppspl(iqgrid,1:2,3,nu)
       end do
     end do
   end do

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(zz)

!  End condition on nproj(/=0)
 end if

!DEBUG
!write(std_out,*)' psp3nl : after d channel '
!stop
!ENDDEBUG

!-----------------------------------------------------------------------
!Treat now f channel (max one projector ! - so do not use ppspl)

!l=3 first projector
 if (abs(h11f)>1.d-12) then
   ekb(4,1)=h11f*(256.0d0/105.0d0)*(rrf**9)*(pi**2.5d0)/(4.d0*pi)**2
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,4,1)=(two_pi*qgrid(iqgrid))**3* &
&     exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrf)**2)
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1j(1)=0d0
   ypnj(1)=(two_pi**3)*qmax**2*exp(-0.5d0*(two_pi*qmax*rrf)**2)*&
&   (3.0d0-(two_pi*qmax*rrf)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,4,1),mqgrid,&
&   yp1j(1),ypnj(1),ffspl(:,2,4,1))
 end if

!-----------------------------------------------------------------------

 ABI_DEALLOCATE(ppspl)
 ABI_DEALLOCATE(work)

!DEBUG
!write(std_out,*)' psp3nl : exit '
!stop
!ENDDEBUG

end subroutine psp3nl
!!***

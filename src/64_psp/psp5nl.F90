!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp5nl
!! NAME
!! psp5nl
!!
!! FUNCTION
!! Make Kleinman-Bylander form factors f_l(q) for each l from
!! 0 to lmax; Vloc is assumed local potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, FrD, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  al=grid spacing in exponent for radial grid
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!   Usually lmax=1, sometimes = 0 (e.g. for oxygen); lmax <= 2 allowed.
!!  mmax=number of radial grid points for atomic grid
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for q grid
!!  qgrid(mqgrid)=values at which form factors are returned
!!  rad(mmax)=radial grid values
!!  vloc(mmax)=local pseudopotential on radial grid
!!  vpspll(mmax,3)=nonlocal pseudopotentials for each l on radial grid
!!  wfll(mmax,3)=reference state wavefunctions on radial grid
!!                mmax and mqgrid
!!
!! OUTPUT
!!  ekb(mpsang)=Kleinman-Bylander energy,
!!             {{\\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!               {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!               \end{equation} }}
!!                for each l
!!  ffspl(mqgrid,2,mpsang)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!
!! NOTES
!! u_l(r) is reference state wavefunction (input as wf);
!! j_l(q) is a spherical Bessel function;
!! dV_l(r) = vpsp_l(r)-vloc(r) for angular momentum l;
!! f_l(q) = $ \int_0^{rmax}[j_l(2\pi q r) u_l(r) dV_l(r) r dr]/\sqrt{dvms}$
!! where dvms = $\int_0^{rmax} [(u_l(r) dV_l(r))^2 dr]$ is the mean
!! square value of the nonlocal correction for angular momentum l.
!! Xavier Gonze s E_KB = $ dvms/\int_0^{rmax}[(u_l(r))^2 dV_l(r) dr]$.
!! This is the eigenvalue of the Kleinman-Bylander operator and sets
!! the energy scale of the nonlocal psp corrections.
!!
!! PARENTS
!!      psp5in,psp6in
!!
!! CHILDREN
!!      ctrap,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp5nl(al,ekb,ffspl,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll)

 use defs_basis
 use m_profiling_abi
 use m_splines
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp5nl'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,mmax,mpsang,mqgrid
 real(dp),intent(in) :: al
!arrays
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vloc(mmax),vpspll(mmax,mpsang)
 real(dp),intent(in) :: wfll(mmax,mpsang)
 real(dp),intent(out) :: ekb(mpsang),ffspl(mqgrid,2,mpsang)

!Local variables-------------------------------
!DEBUG
!real(dp) :: norm,wf
!ENDDEBUG
!scalars
 integer,parameter :: dpsang=5
 integer :: iq,ir,lp1
 real(dp) :: arg,bessel,dvwf,qr,result,yp1,ypn,ztor1
 character(len=500) :: message
!arrays
 real(dp) :: ckb(dpsang),dvms(dpsang),eta(dpsang),renorm(dpsang)
 real(dp),allocatable :: work1(:),work2(:),work3(:),work4(:)

!*************************************************************************

!l=0,1,2 and 3 spherical Bessel functions
!The accuracy of the bes1, bes2, bes3 functions for small arguments
!may be insufficient. In the present version
!of the routines, some care is taken with the value of the argument.
!If smaller than 1.d-3, a two terms
!Taylor series expansion is prefered.
! bes0(arg)=sin(arg)/arg
! bes1(arg)=(sin(arg)-arg*cos(arg))/arg**2
! bes2(arg)=( (3.0d0-arg**2)*sin(arg)-&
!& 3.0d0*arg*cos(arg) )      /arg**3

! bes3(arg)=(15.d0*sin(arg)-15.d0*arg*cos(arg) &
!& -6.d0*arg**2*sin(arg)+arg**3*cos(arg) )/arg**4

!Zero out Kleinman-Bylander energies ekb
 ekb(:)=0.0d0

 ABI_ALLOCATE(work1,(mmax))
 ABI_ALLOCATE(work2,(mmax))
 ABI_ALLOCATE(work3,(mmax))
 ABI_ALLOCATE(work4,(mmax))

!Allow for no nonlocal correction (lmax=-1)
 if (lmax/=-1) then

!  Check that lmax is within allowed range
   if (lmax<0.or.lmax>3) then
     write(message, '(a,i12,a,a,a,a,a,a,a)' )&
&     'lmax=',lmax,' is not an allowed value.',ch10,&
&     'Allowed values are -1 for no nonlocal correction or else',ch10,&
&     '0, 1,2 or 3 for maximum l nonlocal correction.',ch10,&
&     'Action : check the input atomic psp data file for lmax.'
     MSG_ERROR(message)
   end if

!  Compute normalizing integrals eta=<dV> and mean square
!  nonlocal psp correction dvms=<dV^2>
!  "dvwf" consistently refers to dV(r)*wf(r) where dV=nonlocal correction
   do lp1=1,lmax+1

!    integral from 0 to r1
     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)
     ztor1=(wfll(1,lp1)*dvwf)*rad(1)/dble(2*(lp1-1)+3)
!    integrand for r1 to r(mmax) (incl extra factor of r)
     do ir=1,mmax
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)
       work1(ir)=rad(ir)*(wfll(ir,lp1)*dvwf)
     end do
!    do integral by corrected trapezoidal integration
     call ctrap(mmax,work1,al,result)
     eta(lp1)=ztor1+result

     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)
     ztor1=dvwf**2*rad(1)/dble(2*(lp1-1)+3)
     do ir=1,mmax
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)
       work1(ir)=rad(ir)*(dvwf**2)
     end do
     call ctrap(mmax,work1,al,result)
     dvms(lp1)=ztor1+result

!    DEBUG
!    Compute the norm of wfll
!    wf=wfll(1,lp1)
!    ztor1=wf**2*rad(1)/dble(2*(lp1-1)+3)
!    do ir=1,mmax
!    wf=wfll(ir,lp1)
!    work1(ir)=rad(ir)*(wf**2)
!    end do
!    call ctrap(mmax,work1,al,result)
!    norm=ztor1+result
!    write(std_out,*)' lp1, norm',lp1,norm
!    ENDDEBUG

!    If dvms is not 0 for any given angular momentum l,
!    compute Xavier Gonze s definition of the Kleinman-Bylander
!    energy E_KB = dvms/eta.  In this case also renormalize
!    the projection operator to u_KB(r)=$u_l(r)*dV(r)/\sqrt{dvms}$.
!    This means dvwf gets multiplied by the normalization factor
!    "renorm"=$1/\sqrt{dvms}$ as seen below.
     if (dvms(lp1)/=0.0d0) then
       ekb(lp1)=dvms(lp1)/eta(lp1)
       renorm(lp1)=1.0d0/sqrt(dvms(lp1))
!      ckb is Kleinman-Bylander "cosine" (Xavier Gonze)
       ckb(lp1)=eta(lp1)/sqrt(dvms(lp1))
     else
       ekb(lp1)=0.0d0
     end if

   end do

!  l=0 form factor if ekb(1) not 0 (lmax always at least 0)
   if (ekb(1)/=0.0d0) then

!    do q=0 separately
     lp1=1
!    0 to r1 integral
     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
     ztor1=(rad(1)*dvwf)*rad(1)/3.0d0
!    integrand
     do ir=1,mmax
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
       work1(ir)=rad(ir)*(rad(ir)*dvwf)
     end do
     call ctrap(mmax,work1,al,result)
     ffspl(1,1,1)=ztor1+result

!    do rest of q points
     do iq=2,mqgrid
       arg=two_pi*qgrid(iq)
!      0 to r1 integral
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       ztor1=(bes0_psp5(arg*rad(1))*rad(1)*dvwf)*rad(1)/3.0d0
!      integrand
       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
         work1(ir)=rad(ir)*(rad(ir)*bes0_psp5(arg*rad(ir))*dvwf)
       end do
       call ctrap(mmax,work1,al,result)
       ffspl(iq,1,1)=ztor1+result
     end do

!    Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
!    yp1=0 for l=0
     yp1=0.0d0
!    ypn=$ \int [2\pi r (-bes1(2\pi r q)) wf(r) dV(r) r dr]$
     arg=two_pi*qgrid(mqgrid)
     dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
     qr=arg*rad(1)
     if(qr<1.d-3)then
       bessel=(10.d0-qr*qr)*qr/30.0d0
     else
       bessel=bes1_psp5(qr)
     end if
!    ztor1=(-bes1(arg*rad(1))*two_pi*rad(1)*r(1)*dvwf)*rad(1)/5.0d0
     ztor1=(-bessel*two_pi*rad(1)*rad(1)*dvwf)*rad(1)/5.0d0
     do ir=1,mmax
       qr=arg*rad(ir)
       if(qr<1.d-3)then
         bessel=(10.d0-qr*qr)*qr/30.0d0
       else
         bessel=bes1_psp5(qr)
       end if
       dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!      work(ir)=rad(ir)*(-bes1(arg*rad(ir))*two_pi*rad(ir)*rad(ir)*dvwf)
       work1(ir)=rad(ir)*(-bessel*two_pi*rad(ir)*rad(ir)*dvwf)
     end do
     call ctrap(mmax,work1,al,result)
     ypn=ztor1+result

!    Fit spline to get second derivatives by spline fit
     call spline(qgrid,ffspl(1,1,1),mqgrid,yp1,ypn,ffspl(1,2,1))

   else
!    or else put nonlocal correction at l=0 to 0
     ffspl(:,:,1)=0.0d0
   end if

!  Finished if lmax=0 (highest nonlocal correction)
!  Do l=1 form factor if ekb(2) not 0 and lmax>=1
   if (lmax>0)then
     if(ekb(2)/=0.0d0) then

       lp1=2
!      do q=0 separately: f_1(q=0) vanishes !
       ffspl(1,1,2)=0.0d0

!      do rest of q points
       do iq=2,mqgrid
         arg=two_pi*qgrid(iq)
         dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
         qr=arg*rad(1)
         if(qr<1.d-3)then
           bessel=(10.d0-qr*qr)*qr/30.0d0
         else
           bessel=bes1_psp5(qr)
         end if
!        ztor1=(bes1(arg*rad(1))*rad(1)*dvwf)*rad(1)/5.0d0
         ztor1=(bessel*rad(1)*dvwf)*rad(1)/5.0d0

         do ir=1,mmax
           qr=arg*rad(ir)
           if(qr<1.d-3)then
             bessel=(10.d0-qr*qr)*qr/30.0d0
           else
             bessel=bes1_psp5(qr)
           end if
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
           work2(ir)=rad(ir)*(rad(ir)*bessel*dvwf)
         end do

         call ctrap(mmax,work2,al,result)
         ffspl(iq,1,2)=ztor1+result
       end do

!      Compute yp1,ypn for l=1
!      yp1=$\displaystyle \int [2\pi r^2 wf(r) dV(r)]/3$
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       ztor1=((two_pi*rad(1)**2)*dvwf)*rad(1)/(3.0d0*5.0d0)
       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
         work2(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf/3.0d0)
       end do
       call ctrap(mmax,work2,al,result)
       yp1=ztor1+result
!      ypn=$\int [2\pi r^2 wf(r) dV(r) (j_0(x)-(2/x)j_1(x)) dr]$
!      where x=2 Pi qgrid(mqgrid) r
       arg=two_pi*qgrid(mqgrid)
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       qr=arg*rad(1)
       if(qr<1.d-3)then
         bessel=(10.d0-3.0d0*qr*qr)/30.0d0
       else
         bessel=bes0_psp5(qr)-2.d0*bes1_psp5(qr)/qr
       end if
!      ztor1=( (two_pi*rad(1)**2)*dvwf* (bes0(arg*rad(1))-
!      2.0d0*bes1(arg*rad(1))/(arg*rad(1))) ) * rad(1)/5.0d0
       ztor1=( (two_pi*rad(1)**2)*dvwf*bessel)*  rad(1)/5.0d0

       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
         qr=arg*rad(ir)
         if(qr<1.d-3)then
           bessel=(10.d0-3.0d0*qr*qr)/30.0d0
         else
           bessel=bes0_psp5(qr)-2.d0*bes1_psp5(qr)/qr
         end if
!        work(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*
!        (bes0(arg*rad(ir))-2.d0*bes1(arg*rad(ir))/(arg*rad(ir))) )
         work2(ir)=rad(ir)*(two_pi*rad(ir)**2)*dvwf*bessel
       end do
       call ctrap(mmax,work2,al,result)
       ypn=ztor1+result

!      Fit spline for l=1 Kleinman-Bylander form factor
       call spline(qgrid,ffspl(1,1,2),mqgrid,yp1,ypn,ffspl(1,2,2))

     else
!      or else put form factor to 0 for l=1
       ffspl(:,:,2)=0.0d0
     end if
!    Endif condition of lmax>0
   end if

!  Finished if lmax=1 (highest nonlocal correction)
!  Do l=2 nonlocal form factor if eta(3) not 0 and lmax>=2
   if (lmax>1)then
     if(ekb(3)/=0.0d0) then

       lp1=3
!      do q=0 separately; f_2(q=0) vanishes
       ffspl(1,1,3)=0.0d0

!      do rest of q points
       do iq=2,mqgrid
         arg=two_pi*qgrid(iq)
         dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
         qr=arg*rad(1)
         if(qr<1.d-3)then
           bessel=qr*qr/15.0d0-qr**4/210.0d0
         else
           bessel=bes2_psp5(qr)
         end if
!        ztor1=(bes2(arg*rad(1))*rad(1)*dvwf)*rad(1)/7.0d0
         ztor1=(bessel*rad(1)*dvwf)*rad(1)/7.0d0
         do ir=1,mmax
           qr=arg*rad(ir)
           if(qr<1.d-3)then
             bessel=qr*qr/15.0d0-qr**4/210.0d0
           else
             bessel=bes2_psp5(qr)
           end if
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!          work(ir)=rad(ir)*(r(ir)*bes2(arg*rad(ir))*dvwf)
           work3(ir)=rad(ir)*(rad(ir)*bessel*dvwf)
         end do
         call ctrap(mmax,work3,al,result)
         ffspl(iq,1,3)=ztor1+result
       end do

!      Compute yp1,ypn for l=2
!      yp1=0 for l=2
       yp1=0.0d0
!      ypn=$\int [2 \pi r^2 wf(r) dV(r) (j_1(x)-(3/x)j_2(x)) dr]$
!      where x=2 Pi qgrid(mqgrid) r
       arg=two_pi*qgrid(mqgrid)
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       qr=arg*rad(1)
       if(qr<1.d-3)then
         bessel=qr*2.0d0/15.0d0-qr**3*4.0d0/210.0d0
       else
         bessel=bes1_psp5(qr)-3.0d0*bes2_psp5(qr)/qr
       end if
!      ztor1=( (two_pi*rad(1)**2)*dvwf* (bes1(arg*rad(1))-
!      3.0d0*bes2(arg*rad(1))/(arg*rad(1))) ) * rad(1)/7.0d0
       ztor1=( (two_pi*rad(1)**2)*dvwf* bessel ) * rad(1)/7.0d0
       do ir=1,mmax
         qr=arg*rad(ir)
         if(qr<1.d-3)then
           bessel=qr*2.0d0/15.0d0-qr**3*4.0d0/210.0d0
         else
           bessel=bes1_psp5(qr)-3.0d0*bes2_psp5(qr)/qr
         end if
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!        work3(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*
!        (bes1(arg*rad(ir))-3.d0*bes2(arg*rad(ir))/(arg*rad(ir))) )
         work3(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*bessel)
       end do
       call ctrap(mmax,work3,al,result)
       ypn=ztor1+result

!      Fit spline for l=2 Kleinman-Bylander form factor
       call spline(qgrid,ffspl(1,1,3),mqgrid,yp1,ypn,ffspl(1,2,3))

     else
!      or else put form factor to 0 for l=1
       ffspl(:,:,3)=0.0d0
     end if
!    Endif condition of lmax>1
   end if

!  Finished if lmax=2 (highest nonlocal correction)
!  Do l=3 nonlocal form factor if eta(4) not 0 and lmax>=3
   if (lmax>2)then
     if(ekb(4)/=0.0d0) then

       lp1=4
!      do q=0 separately; f_3(q=0) vanishes
       ffspl(1,1,4)=0.0d0

!      do rest of q points
       do iq=2,mqgrid
         arg=two_pi*qgrid(iq)
         dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
         qr=arg*rad(1)
         if(qr<1.d-3)then
           bessel=qr*qr*qr/105.0d0-qr**5/1890.0d0+qr**7/83160.0d0
         else
           bessel=bes3_psp5(qr)
         end if
!        ztor1=(bes3(arg*rad(1))*rad(1)*dvwf)*rad(1)/9.0d0
         ztor1=(bessel*rad(1)*dvwf)*rad(1)/9.0d0
         do ir=1,mmax
           qr=arg*rad(ir)
           if(qr<1.d-3)then
             bessel=qr*qr*qr/105.0d0-qr**5/1890.0d0+qr**7/83160.0d0
           else
             bessel=bes3_psp5(qr)
           end if
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!          work(ir)=rad(ir)*(rad(ir)*bes3(arg*rad(ir))*dvwf)
           work4(ir)=rad(ir)*(rad(ir)*bessel*dvwf)
         end do
         call ctrap(mmax,work4,al,result)
         ffspl(iq,1,4)=ztor1+result
       end do

!      Compute yp1,ypn for l=3
!      yp1=0 for l=3
       yp1=0.0d0
!      ypn=$\int [2\pi r^2 wf(r) dV(r) (j_2(x)-(4/x)j_3(x)) dr]$
!      where x=2 Pi qgrid(mqgrid) r
       arg=two_pi*qgrid(mqgrid)
       dvwf=(vpspll(1,lp1)-vloc(1))*wfll(1,lp1)*renorm(lp1)
       qr=arg*rad(1)
       if(qr<1.d-3)then
         bessel=3.d0*qr**2/105.0d0-5.d0*qr**4/1890.0d0+7.d0*qr**6/83160.0d0
       else
         bessel=bes2_psp5(qr)-4.0d0*bes3_psp5(qr)/qr
       end if
!      ztor1=( (two_pi*rad(1)**2)*dvwf* (bes2(arg*rad(1))-
!      3.0d0*bes3(arg*rad(1))/(arg*rad(1))) ) * rad(1)/9.0d0
       ztor1=( (two_pi*rad(1)**2)*dvwf* bessel ) * rad(1)/9.0d0
       do ir=1,mmax
         qr=arg*rad(ir)
         if(qr<1.d-3)then
           bessel=3.d0*qr**2/105.0d0-5.d0*qr**4/1890.0d0+7.d0*qr**6/83160.0d0
         else
           bessel=bes2_psp5(qr)-4.0d0*bes3_psp5(qr)/qr
         end if
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
!        work4(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*
!        (bes2(arg*rad(ir))-4.d0*bes3(arg*rad(ir))/(arg*rad(ir))) )
         work4(ir)=rad(ir)*((two_pi*rad(ir)**2)*dvwf*bessel)
       end do
       call ctrap(mmax,work4,al,result)
       ypn=ztor1+result

!      Fit spline for l=3 Kleinman-Bylander form factor
       call spline(qgrid,ffspl(1,1,4),mqgrid,yp1,ypn,ffspl(1,2,4))

     else
!      or else put form factor to 0 for l=3
       ffspl(:,:,4)=0.0d0
     end if
!    Endif condition of lmax>2
   end if

!  Endif condition lmax/=-1
 end if

!DEBUG
!write(std_out,*) 'EKB=',(ekb(iq),iq=1,3)
!write(std_out,*) 'COSKB=',(ckb(iq),iq=1,3)
!ENDDEBUG

 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)
 ABI_DEALLOCATE(work3)
 ABI_DEALLOCATE(work4)

 contains
 
   function  bes0_psp5(arg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bes0_psp5'
!End of the abilint section

   real(dp) :: bes0_psp5,arg
   bes0_psp5=sin(arg)/arg
 end function bes0_psp5

   function bes1_psp5(arg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bes1_psp5'
!End of the abilint section

   real(dp) :: bes1_psp5,arg
   bes1_psp5=(sin(arg)-arg*cos(arg))/arg**2
 end function bes1_psp5

   function bes2_psp5(arg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bes2_psp5'
!End of the abilint section

   real(dp) :: bes2_psp5,arg
   bes2_psp5=( (3.0d0-arg**2)*sin(arg)-&
&   3.0d0*arg*cos(arg) )      /arg**3
 end function bes2_psp5

   function bes3_psp5(arg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bes3_psp5'
!End of the abilint section

   real(dp) :: bes3_psp5, arg
   bes3_psp5=(15.d0*sin(arg)-15.d0*arg*cos(arg) &
&   -6.d0*arg**2*sin(arg)+arg**3*cos(arg) )/arg**4
 end function bes3_psp5

end subroutine psp5nl
!!***

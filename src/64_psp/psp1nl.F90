!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp1nl
!! NAME
!! psp1nl
!!
!! FUNCTION
!! Make Kleinman-Bylander form factors f_l(q) for each l from
!! 0 to lmax; Vloc is assumed local potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dr(mmax)=inverse of grid spacing for radial grid
!!  lloc=angular momentum of local channel (avoid doing integrals for this l)
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!  mmax=number of radial grid points for atomic grid
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for q grid
!!  qgrid(mqgrid)=values at which form factors are returned
!!  rad(mmax)=radial grid values
!!  vloc(mmax)=local pseudopotential on radial grid
!!  vpspll(mmax,lmax+1)=nonlocal pseudopotentials for each l on radial grid
!!  wfll(mmax,lmax+1)=reference state wavefunctions on radial grid
!!  wksincos(mmax,2,2)=contains sine and cosine of 2*pi*r(:)*dq and 2*pi*r(:)*q
!!    at input :  wksincos(:,1,1)=sine of 2*pi*r(:)*dq
!!                wksincos(:,2,1)=cosine of 2*pi*r(:)*dq
!!    wksincos(:,:,2) is not initialized, will be used inside the routine
!!
!! OUTPUT
!!  ekb(mpsang)=Kleinman-Bylander energy,
!!              {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!              {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!               \end{equation} }}
!!              for each l
!!  ffspl(mqgrid,2,mpsang)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!
!! NOTES
!! u_l(r) is reference state wavefunction (input as wfll);
!! j_l(q) is a spherical Bessel function;
!! dV_l(r) = vpsp_l(r)-vloc(r) for angular momentum l;
!! f_l(q) =$ \int_0^{rmax}[j_l(2\pi q r) u_l(r) dV_l(r) r dr]/\sqrt{dvms}$
!! where dvms=$\displaystyle \int_0^{rmax}[(u_l(r) dV_l(r))^2 dr]$ is the mean
!! square value of the nonlocal correction for angular momentum l.
!! E_KB = $\displaystyle \frac{dvms}{\int_0^{rmax}[(u_l(r))^2 dV_l(r) dr]}$.
!! This is the eigenvalue of the Kleinman-Bylander operator and sets
!! the energy scale of the nonlocal psp corrections.
!! Bessel functions replaced by besj, which accomodates args near 0.
!!
!! PARENTS
!!      psp1in
!!
!! CHILDREN
!!      besjm,der_int,sincos,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp1nl(dr,ekb,ffspl,lloc,lmax,mmax,mpsang,mqgrid,&
&                  qgrid,rad,vloc,vpspll,wfll,wksincos)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_splines

 use m_special_funcs,   only : besjm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp1nl'
 use interfaces_64_psp, except_this_one => psp1nl
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,mmax,mpsang,mqgrid
!arrays
 real(dp),intent(in) :: dr(mmax),qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(in) :: vpspll(mmax,mpsang),wfll(mmax,mpsang)
 real(dp),intent(inout) :: wksincos(mmax,2,2)
 real(dp),intent(out) :: ekb(mpsang),ffspl(mqgrid,2,mpsang)

!Local variables-------------------------------
!scalars
 integer,parameter :: dpsang=5
 integer :: iq,ir,irmax,lp1
 real(dp) :: dvwf,result,test,tpiq,yp1,ypn
 character(len=500) :: message
!arrays
 real(dp) :: ckb(dpsang),dvms(dpsang),eta(dpsang),renorm(dpsang)
 real(dp),allocatable :: besjx(:),work1(:),work2(:),work3(:),work4(:),work5(:)
 real(dp),allocatable :: work_spl(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' psp1nl : enter'
!stop
!ENDDEBUG

!Zero out Kleinman-Bylander energies ekb
 ekb(:)=0.0d0
!Zero out eta and other parameters too (so 0 s show up in output later)
 eta(:)=0.0d0
 dvms(:)=0.0d0
 ckb(:)=0.0d0

!Allow for no nonlocal correction (lmax=-1)
 if (lmax/=-1) then

!  Check that lmax is within allowed range
   if (lmax<0.or.lmax>3) then
     write(message, '(a,i12,a,a,a,a,a,a,a)' )&
&     'lmax=',lmax,' is not an allowed value.',ch10,&
&     'Allowed values are -1 for no nonlocal correction or else',ch10,&
&     '0, 1, 2, or 3 for maximum l nonlocal correction.',ch10,&
&     'Action: check the input atomic psp data file for lmax.'
     MSG_ERROR(message)
   end if

!  Compute normalizing integrals eta=<dV> and mean square
!  nonlocal psp correction dvms=<dV^2>
!  "dvwf" consistently refers to dV(r)*wf(r) where dV=nonlocal correction

   ABI_ALLOCATE(work1,(mmax+1))
   ABI_ALLOCATE(work2,(mmax+1))
   ABI_ALLOCATE(work_spl,(mqgrid))
   ABI_ALLOCATE(work5,(mmax))
   ABI_ALLOCATE(besjx,(mmax))

   do lp1=1,lmax+1

!    Only do the work if nonlocal correction is nonzero
     if (lp1 /= lloc+1) then

!      integrand for 0 to r(mmax)
       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)
         work1(ir)=wfll(ir,lp1)*dvwf
       end do

!      do integral
!      first need derivative of function; note use of
!      shifted indices to accomodate Mike Teter s choice of 0:mmax-1
       call der_int(work1,work2,rad,dr,mmax-1,result)
       eta(lp1)=result

!      DEBUG
!      write(std_out,*)' psp1nl : write eta(lp1)'
!      write(std_out,*)result
!      do ir=1,mmax,61
!      write(std_out,*)vpspll(ir,lp1),vloc(ir),wfll(ir,lp1)
!      end do
!      write(std_out,*)
!      do ir=1,mmax,61
!      write(std_out,*)work1(ir),rad(ir),dr(ir)
!      end do
!      ENDDEBUG

       do ir=1,mmax
         dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)
         work1(ir)=dvwf**2
       end do
       call der_int(work1,work2,rad,dr,mmax-1,result)

       dvms(lp1)=result

!      If dvms is not 0 for any given angular momentum l,
!      compute Xavier Gonze s definition of the Kleinman-Bylander
!      energy E_KB = dvms/eta.  In this case also renormalize
!      the projection operator to u_KB(r)=$u_l(r) dV(r)/\sqrt{dvms}$.
!      This means dvwf gets multiplied by the normalization factor
!      "renorm"=$1/\sqrt{dvms}$ as seen below.
!      With dvwf=dV(r)*wf(r) for wf(r)=``radial'' wf, the integrand
!      for each angular momentum l is
!      Bessel_l(2 $\pi$ q r) * wf(r) * dV(r) * r;
!      NOTE presence of extra r in integrand.

       if (dvms(lp1)/=0.0d0) then
         ekb(lp1)=dvms(lp1)/eta(lp1)
         renorm(lp1)=1.0d0/sqrt(dvms(lp1))
!        ckb is Kleinman-Bylander "cosine" (Xavier Gonze)
         ckb(lp1)=eta(lp1)/sqrt(dvms(lp1))
       else
         ekb(lp1)=0.0d0
       end if
     end if
   end do

!  Loop on angular momenta
   do lp1=1,lmax+1

!    Compute form factor if ekb(lp1) not 0
     if (ekb(lp1)/=0.0d0) then

!      do q=0 separately, non-zero if l=0
       if(lp1==1)then
         do ir=1,mmax
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
           work1(ir)=rad(ir)*dvwf
         end do
         call der_int(work1,work2,rad,dr,mmax-1,result)
         ffspl(1,1,lp1)=result
       else
!        For l non-zero, f(q=0) vanishes !
         ffspl(1,1,lp1)=0.0d0
       end if

!      Prepare loop over q values
       irmax=mmax+1
       do ir=mmax,2,-1
         test=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)*rad(ir)
         work5(ir)=test
         work1(ir)=0.0d0
!        Will ignore tail within decade of machine precision
         if ((10.0d0+abs(test))==10.0d0 .and. irmax==ir+1) then
           irmax=ir
         end if
       end do
!      Increase irmax a bit
       irmax=irmax+4
!      Ask irmax to be lower than mmax
       if(irmax>mmax-1)irmax=mmax-1

       ABI_ALLOCATE(work3,(irmax-1))
       ABI_ALLOCATE(work4,(irmax-1))

!      Loop over q values
       do iq=2,mqgrid
         tpiq=two_pi*qgrid(iq)
         call sincos(iq,irmax,mmax,wksincos,rad,tpiq)
         work3(:)=wksincos(2:irmax,2,2) !Temporary array (Intel compiler compatibility)
         work4(:)=wksincos(2:irmax,1,2) !Temporary array (Intel compiler compatibility)

!        Handle r=0 separately
         work1(1)=0.0d0
         call besjm(tpiq,besjx(2:irmax),work3,(lp1-1),irmax-1,work4,rad(2:irmax))
         do ir=2,irmax
           work1(ir)=besjx(ir)*work5(ir)
         end do
!        do integral
         call der_int(work1,work2,rad,dr,irmax,result)
         ffspl(iq,1,lp1)=result
       end do

!      Compute yp1=derivative of f(q) at q=0
       if(lp1/=2)then
!        For l/=1, yp1=0
         yp1=0.0d0
       else
!        For l=1, yp1=Int [2 Pi r^2 wf(r) dV(r)]/3
         do ir=1,irmax
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
           work1(ir)=(two_pi*rad(ir)**2)*dvwf/3.0d0
         end do
         call der_int(work1,work2,rad,dr,irmax,result)
         yp1=result
       end if

!      Compute ypn=derivative of f(q) at q=qgrid(mqgrid)
       tpiq=two_pi*qgrid(mqgrid)
!      Treat ir=1, r=0, separately
       work1(1)=0.0d0
!      Here, must distinguish l==0 from others
       if(lp1==1)then
!        l==0 : ypn=$\int [2\pi r (-bes1(2\pi r q)) wf(r) dV(r) r dr]$
!        The sine and cosine of the last point were computed in the previous loop
!        So, there is no need to call sincos. Note that the rank of besj is 1.
         call besjm(tpiq,besjx(2:irmax),work3,1,irmax-1,work4,rad(2:irmax))
         do ir=2,irmax
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
           work1(ir)=-besjx(ir)*two_pi*rad(ir)*rad(ir)*dvwf
         end do
       else
!        l==1 : ypn=$\int [2\pi r^2 wf(r) dV(r) (j_0(x)-(2/x)j_1(x)) dr]$
!        l==2 : ypn=$\int [2\pi r^2 wf(r) dV(r) (j_1(x)-(3/x)j_2(x)) dr]$
!        l==3 : ypn=$\int [2\pi r^2 wf(r) dV(r) (j_2(x)-(4/x)j_3(x)) dr]$
!        The sine and cosine of the last point were computed in the previous loop
!        Store first previously computed value with besj of order l, then use
!        besj of order l-1 (=lp1-2)
         work1(2:irmax)=besjx(2:irmax)
         call besjm(tpiq,besjx(2:irmax),work3,(lp1-2),irmax-1,work4,rad(2:irmax))
         do ir=2,irmax
           dvwf=(vpspll(ir,lp1)-vloc(ir))*wfll(ir,lp1)*renorm(lp1)
           work1(ir)=(two_pi*rad(ir)**2)*dvwf*&
&           ( besjx(ir) - ( dble(lp1)*work1(ir)/(tpiq*rad(ir)) ) )
         end do
       end if
!      work1 is ready for integration
       call der_int(work1,work2,rad,dr,irmax,result)
       ypn=result

!      Fit spline to get second derivatives by spline fit
       call spline(qgrid,ffspl(:,1,lp1),mqgrid,yp1,ypn,&
&       ffspl(:,2,lp1))

       ABI_DEALLOCATE(work3)
       ABI_DEALLOCATE(work4)

     else

!      KB energy is zero, put nonlocal correction at l=0 to 0
       ffspl(:,:,lp1)=0.0d0

     end if

!    End loop on angular momenta
   end do

   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
   ABI_DEALLOCATE(work_spl)
   ABI_DEALLOCATE(work5)
   ABI_DEALLOCATE(besjx)

!  End of lmax/=-1 condition
 end if

end subroutine psp1nl
!!***

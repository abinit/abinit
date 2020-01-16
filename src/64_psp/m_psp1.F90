!!****m* ABINIT/m_psp1
!! NAME
!!  m_psp1
!!
!! FUNCTION
!!  Initialize pspcod=1 or 4 pseudopotential (Teter format)
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, FrD, MT)
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

module m_psp1

 use defs_basis
 use m_errors
 use m_abicore
 use m_splines

 use m_special_funcs,   only : besjm
 use m_psptk,           only : psp1cc

 implicit none

 private
!!***

 public :: psp1in       ! Initialize pspcod=1 or 4 pseudopotential (Teter format)
!!***

contains
!!***

!!****f* m_psp1/psp1in
!! NAME
!! psp1in
!!
!! FUNCTION
!! Initialize pspcod=1 or 4 pseudopotential (Teter format):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! INPUTS
!!  dq= spacing of the q-grid
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pspcod=pseudopotential type
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=atomic number of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  ekb1(mpsang)= Kleinman-Bylander energy from the psp file, for iproj=1
!!  ekb2(mpsang)= Kleinman-Bylander energy from the psp file, for iproj=2
!!  epsatm=$(4\pi) \int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  epspsp(mpsang)=values of epsatm for different angular momenta, from the psp file
!!  e990(mpsang)=ecut at which 0.99 of the kinetic energy is recovered
!!  e999(mpsang)=ecut at which 0.999 of the kinetic energy is recovered
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  qchrg is the total (integrated) core charge
!!  rcpsp(mpsang)=cut-off radius for each angular momentum
!!  rms(mpsang)=root mean square of the KB psp
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!  xcccrc=XC core correction cutoff radius (bohr)
!!
!! NOTES
!! there are only minor differences in the two formats
!! 1) With pspcod=1, even for the LOCAL angular momentum, there is
!!    a block for the wfs (can be set to zero, though)
!! 2) The core charge density differs: for pspcod=1, it is a
!!    revised expression for core density of 5 Nov 1992, while
!!    for pspcod=4, it is an older expression, of 7 May 1992 .
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      psp1cc,psp1lo,psp1nl,psp4cc,spline,wrtout
!!
!! SOURCE

subroutine psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&
&                  e990,e999,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mqgrid,nproj,n1xccc,pspcod,&
&                  qchrg,qgrid,rcpsp,rms,useylm,vlspl,xcccrc,xccc1d,&
&                  zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mmax,mpsang,mqgrid,n1xccc,pspcod
 integer,intent(in) :: useylm
 real(dp),intent(in) :: dq,zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpsang)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: e990(mpsang),e999(mpsang),ekb(lnmax),ekb1(mpsang)
 real(dp),intent(out) :: ekb2(mpsang),epspsp(mpsang)
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)
 real(dp),intent(out) :: rcpsp(mpsang),rms(mpsang),vlspl(mqgrid,2)
 real(dp),intent(inout) :: xccc1d(n1xccc,6)

!Local variables-------------------------------
!scalars
 integer :: ii,iln,index,ipsang,kk,lhigh,ll,mm,nlmax
 real(dp) :: arg,dq2pi,fchrg,rchrg,xx,yp1,ypn
 character(len=500) :: message,errmsg
!arrays
 real(dp),allocatable :: drad(:),ekb_tmp(:,:),ffspl_tmp(:,:,:,:),rad(:),vloc(:)
 real(dp),allocatable :: vpspll(:,:),wfll(:,:),wksincos(:,:,:),work_space(:)
 real(dp),allocatable :: work_spl1(:),work_spl2(:)

! ***************************************************************************

!Note: Teter s grid is hard-coded at mmax=2001
!mmax was read from the pseudopotential file in the calling routine
 if (mmax/=2001) then
   write(message, '(a,i12,a,a,a,a)' )&
&   'Using Teter grid (pspcod=1 and 4) but mmax=',mmax,ch10,&
&   'mmax must be 2001 for Teter grid.',ch10,&
&   'Action: check your pseudopotential input file.'
   MSG_ERROR(message)
 end if

!File format of formatted Teter psp input (the 3 first lines
!have already been read in calling -pspatm- routine) :

!(1) title (character) line
!(2) znucl,zion,pspdat
!(3) pspcod,pspxc,lmax,lloc,mmax,r2well
!For each angular momentum :
!(4) ll,e990(ll),e999(ll),nproj(ll),rcpsp(ll)
!(5) rms(ll),ekb1(ll),ekb2(ll),epspsp(ll)
!(6) rchrg,fchrg,qchrg
!(7) ll
!(8) (vpsp(j,ll),j=0,nmax)
!Then for iproj=1 to 2
!for ll=0,lmax
!(10) ll
!(11) ((upsp(j,ll,iproj),j=0,nmax)

 do ipsang=1,lmax+1

   read (tmp_unit,*,err=10,iomsg=errmsg) ll,e990(ipsang),e999(ipsang),nproj(ipsang),rcpsp(ipsang)
   write(message, '(i5,2f8.3,i5,f12.7,t47,a)' ) &
&   ipsang-1,e990(ipsang),e999(ipsang),nproj(ipsang),rcpsp(ipsang),&
&   'l,e99.0,e99.9,nproj,rcpsp'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*,err=10,iomsg=errmsg) rms(ipsang),ekb1(ipsang),ekb2(ipsang),epspsp(ipsang)
   write(message, '(4f13.8,t55,a)' ) &
&   rms(ipsang),ekb1(ipsang),ekb2(ipsang),epspsp(ipsang),&
&   '   rms, ekb1, ekb2, epsatm'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end do

!Initialize array indlmn array giving l,m,n,lm,ln,s for i=lmn
 index=0;iln=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
   if(nproj(ipsang)>0)then
     ll=ipsang-1
     do kk=1,nproj(ipsang)
       iln=iln+1
       do mm=1,2*ll*useylm+1
         index=index+1
         indlmn(1,index)=ll
         indlmn(2,index)=mm-ll*useylm-1
         indlmn(3,index)=kk
         indlmn(4,index)=ll*ll+(1-useylm)*ll+mm
         indlmn(5,index)=iln
         indlmn(6,index)=1
       end do
     end do
   end if
 end do

 read (tmp_unit,*,err=10,iomsg=errmsg) rchrg,fchrg,qchrg
 write(message, '(3f20.14,t64,a)' ) rchrg,fchrg,qchrg,'rchrg,fchrg,qchrg'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 ! Generate core charge function and derivatives, if needed
 if(fchrg>1.0d-15)then
   if(pspcod==1)then
     call psp1cc(fchrg,n1xccc,xccc1d)
     ! The core charge function for pspcod=1 becomes zero beyond 3*rchrg only.
     ! Thus xcccrc must be set equal to 3*rchrg .
     xcccrc=3*rchrg
   else if(pspcod==4)then
     call psp4cc(fchrg,n1xccc,xccc1d)
     ! For pspcod=4, the core charge cut off exactly beyond rchrg
     xcccrc=rchrg
   end if
 else
   xcccrc=0.0d0
   xccc1d(:,:)=0.0d0
 end if

!--------------------------------------------------------------------
!Will now proceed at the reading of pots and wfs, as well as their treatment

!vpspll(:,1),...,vpspll(:,4)=nonlocal pseudopotentials
!vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
!rad(:)=radial grid r(i)
!drad(:)= inverse of d(r(i))/d(i) for radial grid
!wfll(:,1),...,wfll(:,4)=reference config. wavefunctions

 ABI_ALLOCATE(vloc,(mmax))
 ABI_ALLOCATE(vpspll,(mmax,mpsang))
 if(lmax==-1) vpspll(:,:)=zero

!(1) Read atomic pseudopotential for each l, filling up array vpspll
!Note: put each l into vpspll(:,l+1)
 do ipsang=1,lmax+1
   read (tmp_unit,*,err=10,iomsg=errmsg) ll
   read (tmp_unit,*,err=10,iomsg=errmsg) (vpspll(ii,ipsang),ii=1,mmax)
 end do

!Copy appropriate nonlocal psp for use as local one
 vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )

!(2) Create radial grid, and associated quantities
 ABI_ALLOCATE(rad,(mmax))
 ABI_ALLOCATE(drad,(mmax))
 ABI_ALLOCATE(wksincos,(mmax,2,2))

!Teter grid--need both r and dr in this case
 do ii=0,mmax-1
   xx=dble(ii)/dble(mmax-1)
   rad (ii+1)=100.d0*(xx+.01d0)**5-1.d-8
   drad(ii+1)=500.d0*(xx+.01d0)**4/dble(mmax-1)
 end do

!here compute sin(r(:)*dq) and cos(r(:)*dq)
!NOTE: also invert dr !!
 dq2pi=2.0d0*pi*dq
 do ii=1,mmax
   arg=dq2pi*rad(ii)
   drad(ii)=1.0d0/drad(ii)
   wksincos(ii,1,1)=sin(arg)
   wksincos(ii,2,1)=cos(arg)
 end do

!(3)Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform) to get q^2 V(q).
 ABI_ALLOCATE(work_space,(mqgrid))
 ABI_ALLOCATE(work_spl1,(mqgrid))
 ABI_ALLOCATE(work_spl2,(mqgrid))
 call psp1lo(drad,epsatm,mmax,mqgrid,qgrid,&
& work_spl1,rad,vloc,wksincos,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 call spline (qgrid,work_spl1,mqgrid,yp1,ypn,work_spl2)
 vlspl(:,1)=work_spl1(:)
 vlspl(:,2)=work_spl2(:)

 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl1)
 ABI_DEALLOCATE(work_spl2)

!(4)Take care of non-local part

!Zero out all Kleinman-Bylander energies to initialize
 ekb(:)=0.0d0
!write(std_out,*)' psp1in : before nonlocal corrections '
!write(std_out,*)' psp1in : lloc, lmax = ',lloc,lmax

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 else

!  Proceed to make Kleinman-Bylander form factors for each l up to lmax

!  Read wavefunctions for each l up to lmax
   ABI_ALLOCATE(wfll,(mmax,mpsang))
   do ipsang=1,lmax+1
!    For pspcod==4, wfs for the local angular momentum are not written
     if (nproj(ipsang)/=0 .or. pspcod==1) then
       read (tmp_unit,*,err=10,iomsg=errmsg) ll
       if (ipsang/=ll+1) then
         write(message, '(a,a,a,a,a,a,2i6,a,a)' )&
&         'Pseudopotential input file does not have',ch10,&
&         'angular momenta in order expected for first projection',&
&         'operator.',ch10,' Values are ',ipsang-1,ll,ch10,&
&         'Action: check your pseudopotential input file.'
         MSG_ERROR(message)
       end if
       read (tmp_unit,*,err=10,iomsg=errmsg) wfll(:,ipsang)

     else
       wfll(:,ipsang)=0.0d0
     end if

   end do
!  ----------------------------------------------------------------------
!  Compute KB form factors and fit splines

!  nlmax is highest l for which a nonlocal correction is being computed
   nlmax=lmax
   if (lloc==lmax) nlmax=lmax-1
!  write(std_out,*)' psp1in : lmax,lloc=',lmax,lloc
   ABI_ALLOCATE(ekb_tmp,(mpsang,2))
   ABI_ALLOCATE(ffspl_tmp,(mqgrid,2,nlmax+1,2))

   call psp1nl(drad,ekb_tmp(:,1),ffspl_tmp(:,:,:,1),lloc,&
&   nlmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll,wksincos)

!  Read second wavefunction for second projection operator
!  (only read cases where nproj(ll)=2) --also find highest l for which nproj(l)=2
   lhigh=-1
   do ipsang=1,min(lmax+1,mpsang)
     if (nproj(ipsang)==2) then
       lhigh=ipsang-1
       read (tmp_unit,*,err=10,iomsg=errmsg) ll
       if (ipsang/=ll+1) then
         write(message, '(a,a,a,a,a,a,2i6,a,a)' )&
&         'Pseudopotential input file does not have',ch10,&
&         'angular momenta in order expected for second projection',&
&         'operator.',ch10,' Values are ',ipsang-1,ll,ch10,&
&         'Action: check your pseudopotential input file.'
         MSG_ERROR(message)
       end if
       read (tmp_unit,*,err=10,iomsg=errmsg) wfll(:,ipsang)

     else
       wfll(:,ipsang)=0.0d0

     end if
   end do

!  Compute KB form factors and fit splines for second wf if any

   if (lhigh>-1) then
     call psp1nl(drad,ekb_tmp(:,2),ffspl_tmp(:,:,:,2),lloc,&
&     lhigh,mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll,wksincos)
   end if

!  Convert ekb and ffspl
   iln=0
   do ii=1,lmnmax
     kk=indlmn(5,ii)
     if (kk>iln) then
       iln=kk
       ekb(kk)=ekb_tmp(1+indlmn(1,ii),indlmn(3,ii))
!      write(std_out,*)' psp1in : lmnmax,ii,indlmn(1,ii)=',lmnmax,ii,indlmn(1,ii)
       ffspl(:,:,kk)=ffspl_tmp(:,:,1+indlmn(1,ii),indlmn(3,ii))
     end if
   end do

   ABI_DEALLOCATE(ekb_tmp)
   ABI_DEALLOCATE(ffspl_tmp)
   ABI_DEALLOCATE(wfll)
 end if

 ABI_DEALLOCATE(vpspll)
 ABI_DEALLOCATE(rad)
 ABI_DEALLOCATE(drad)
 ABI_DEALLOCATE(vloc)
 ABI_DEALLOCATE(wksincos)

 return

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine psp1in
!!***

!!****f* m_psp1/psp1lo
!! NAME
!! psp1lo
!!
!! FUNCTION
!! Compute sine transform to transform from v(r) to q^2 v(q)
!! using subroutines related to Teter atomic structure grid.
!!
!! INPUTS
!!  drad(mmax)=inverse of r grid spacing at each point
!!  mmax=number of radial r grid points (Teter grid)
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rad(mmax)=r grid values (bohr).
!!  vloc(mmax)=v(r) on radial grid.
!!  wksincos(mmax,2,2)=contains sine and cosine of 2*pi*r(:)*dq and 2*pi*r(:)*q
!!    at input :  wksincos(:,1,1)=sine of 2*pi*r(:)*dq
!!                wksincos(:,2,1)=cosine of 2*pi*r(:)*dq
!!    wksincos(:,:,2) is not initialized, will be used inside the routine
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm= $4\pi \int[r^2 (v(r)+Zv/r) dr]$
!!  q2vq(mqgrid)=$q^2 v(q)$
!!  =$\displaystyle -Zv/\pi+q^2 4\pi\int(\frac{\sin(2\pi q r)}{2 \pi q r})(r^2 v(r)+r Zv)dr$.
!!  yp1,ypn=derivative of q^2 v(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      psp1in
!!
!! CHILDREN
!!      der_int,sincos
!!
!! SOURCE

subroutine psp1lo(drad,epsatm,mmax,mqgrid,qgrid,q2vq,rad,&
&  vloc,wksincos,yp1,ypn,zion)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,mqgrid
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm,yp1,ypn
!arrays
 real(dp),intent(in) :: drad(mmax),qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(inout) :: wksincos(mmax,2,2)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables-------------------------------
!scalars
 integer,parameter :: mma0=2001
 integer :: iq,ir,irmax
 real(dp),parameter :: scale=10.0d0
 real(dp) :: result,test,tpiq
!arrays
 real(dp) :: wk(mma0),wk1(mma0),wk2(mma0)

! *************************************************************************

!Do q=0 separately (compute epsatm)
!Set up integrand for q=0: Int[r^2 (V(r)+Zv/r) dr]
!Treat r=0 by itself
 wk(1)=0.0d0

 do ir=2,mmax
!  (at large r do not want prefactor of r^2 and should see
!  V(r)+Zv/r go to 0 at large r)
   test=vloc(ir)+zion/rad(ir)
!  write(std_out,'(i4,3es20.10)' )ir,rad(ir),test,rad(ir)*test
!  In this routine, NO cut-off radius is imposed : the input
!  vloc MUST be in real(dp) to obtain numerically
!  accurate values. The error can be on the order of 0.001 Ha !
   if (abs(test)<1.0d-20) then
     wk(ir)=0.0d0
   else
     wk(ir)=rad(ir)*(rad(ir)*vloc(ir)+zion)
   end if
 end do
!Do integral from 0 to r(max) (disregard contrib beyond r(max)
!(need numerical derivatives to do integral)
!Use mmax-1 to convert to Teter s dimensioning starting at 0
 call der_int(wk,wk2,rad,drad,mmax-1,result)

 epsatm=4.d0*pi*(result)
!q=0 value of integral is -zion/Pi + q^2 * epsatm = -zion/Pi
 q2vq(1)=-zion/pi

!Prepare loop over q values
 irmax=mmax+1
 do ir=mmax,2,-1
   test=vloc(ir)+zion/rad(ir)
   wk1(ir)=test*rad(ir)
!  Will ignore tail within decade of machine precision
   if ((scale+abs(test))==scale .and. irmax==ir+1) then
     irmax=ir
   end if
 end do
!Increase irmax a bit : this is copied from psp1nl
 irmax=irmax+4
 if(irmax>mmax)irmax=mmax

!Loop over q values
 do iq=2,mqgrid
   tpiq=two_pi*qgrid(iq)
   call sincos(iq,irmax,mmax,wksincos,rad,tpiq)
!  set up integrand Sin(2Pi q r)(rV(r)+Zv) for integral
!$\displaystyle -Zv/\pi + q^2 4\pi \int[\frac{\sin(2\pi q r)}{2\pi q r}(r^2 v(r)+r Zv)dr]$.
!  Handle r=0 separately
   wk(1)=0.0d0
   do ir=2,irmax
     wk(ir)=wksincos(ir,1,2)*wk1(ir)
   end do
!  do integral from 0 to r(max)
   if(irmax>mmax-1)irmax=mmax-1

   call der_int(wk,wk2,rad,drad,irmax,result)
!  store q^2 v(q)
   q2vq(iq)=-zion/pi+2.d0*qgrid(iq)*result
 end do

!Compute derivatives of q^2 v(q) at ends of interval
 yp1=0.0d0
!ypn=$\displaystyle 2\int_0^\infty (\sin (2\pi qmax r)+(2\pi qmax r)\cos (2\pi qmax r)(r V(r)+Z)dr]$
!integral from r(mmax) to infinity is overkill; ignore
!set up integrand
!Handle r=0 separately
 wk(1)=0.0d0
 tpiq=two_pi*qgrid(mqgrid)
 do ir=2,mmax
   test=vloc(ir)+zion/rad(ir)
!  Ignore contributions within decade of machine precision
   if ((scale+abs(test))==scale) then
     wk(ir)=0.0d0
   else
     wk(ir)=(sin(tpiq*rad(ir))+tpiq*rad(ir)*cos(tpiq*rad(ir))) * &
&     (rad(ir)*vloc(ir)+zion)
   end if
 end do
 call der_int(wk,wk2,rad,drad,mmax-1,result)

 ypn=2.0d0*result

end subroutine psp1lo
!!***

!!****f* m_psp1/psp1nl
!! NAME
!! psp1nl
!!
!! FUNCTION
!! Make Kleinman-Bylander form factors f_l(q) for each l from
!! 0 to lmax; Vloc is assumed local potential.
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

subroutine psp1nl(dr,ekb,ffspl,lloc,lmax,mmax,mpsang,mqgrid,&
&                  qgrid,rad,vloc,vpspll,wfll,wksincos)

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

   end do !    End loop on angular momenta

   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
   ABI_DEALLOCATE(work_spl)
   ABI_DEALLOCATE(work5)
   ABI_DEALLOCATE(besjx)
 end if !  End of lmax/=-1 condition

end subroutine psp1nl
!!***

!!****f* m_psp1/der_int
!! NAME
!! der_int
!!
!! FUNCTION
!! Given input function f(i) on Teter radial grid, and grid spacing
!! dr(i), compute function derivative df/dr on points from 0 to n.
!! Integrate function f(i) on grid r(i) from r(0) to r(nlast).
!! Note that array dimensions start at 0.
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

subroutine der_int(ff,df,rr,dr,nlast,smf)

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

!!****f* m_psp1/sincos
!! NAME
!! sincos
!!
!! FUNCTION
!! Update the sine and cosine values, needed inside the
!! pseudopotential routines psp1lo and psp1nl.
!!
!! INPUTS
!!  iq  = number of current wavevector q
!!  irmax = number of values  of r on the radial grid to be computed
!!  mmax = dimension of pspwk and rad
!!  pspwk(:,1,1) and pspwk(:,2,1) : sine and cosine of 2$\pi$ dq * rad
!!  pspwk(:,1,2) and pspwk(:,2,2) : sine and cosine of 2$\pi$ previous q * rad
!!  rad(mmax) radial grid
!!  tpiq = 2 $\pi$ * current wavevector q
!!
!! OUTPUT
!!  pspwk(*,1,2) and pspwk(*,2,2) : sine and cosine of 2$\pi$ current q * rad
!!
!! NOTES
!! The speed was a special concern, so iterative computation
!! based on addition formula is possible. Interestingly,
!! this algorithm places strong constraints on accuracy,
!! so this routine is machine-dependent.
!!
!! PARENTS
!!      psp1lo,psp1nl
!!
!! CHILDREN
!!
!! SOURCE

subroutine sincos(iq,irmax,mmax,pspwk,rad,tpiq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,irmax,mmax
 real(dp),intent(in) :: tpiq
!arrays
 real(dp),intent(in) :: rad(mmax)
 real(dp),intent(inout) :: pspwk(mmax,2,2)

!Local variables-------------------------------
!scalars
 integer :: ir,nstep
 real(dp) :: prevcos,prevsin
 logical :: testmipspro


! *************************************************************************

 if(iq==2)then

!  Here set up the sin and cos at iq=2
   do ir=2,irmax
     pspwk(ir,1,2)=pspwk(ir,1,1)
     pspwk(ir,2,2)=pspwk(ir,2,1)
   end do

 else
!
!  The sensitivity of the algorithm to changes of nstep
!  has been tested : for all the machines except SGI - R10000 ,
!  either using only the hard way, or
!  using up to nstep=40 causes changes at the level
!  of 1.0d-16 in the total energy. Larger values of
!  nstep might be possible, but the associated residual
!  is already very small ! The accelerated computation of
!  sine and cosine is essential for a good speed on IBM, but,
!  fortunately, on the SGI - R10000 the normal computation is fast enough.

   testmipspro=.false.
   nstep=40
   if(iq-(iq/nstep)*nstep == 0 .or. testmipspro)then

!    Every nstep steps, uses the hard way
     do ir=2,irmax
       pspwk(ir,1,2)=sin(tpiq*rad(ir))
       pspwk(ir,2,2)=cos(tpiq*rad(ir))
     end do

   else

!    Here the fastest way, iteratively
     do ir=2,irmax
       prevsin=pspwk(ir,1,2)
       prevcos=pspwk(ir,2,2)
       pspwk(ir,1,2)=prevsin*pspwk(ir,2,1)+prevcos*pspwk(ir,1,1)
       pspwk(ir,2,2)=prevcos*pspwk(ir,2,1)-prevsin*pspwk(ir,1,1)
     end do

   end if

 end if ! iq==2

end subroutine sincos
!!***

!!****f* m_psp1/psp4cc
!! NAME
!! psp4cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 4 of pseudopotentials.
!! This is a even polynomial of 24th order for core density,
!! that is cut off exactly beyond rchrg.
!! It has been produced on 7 May 1992 by M. Teter.
!!
!! INPUTS
!!  fchrg=magnitude of the core charge correction
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! NOTES
!! The argument of xccc1d is assumed to be normalized, and to vary
!! from xx=0 to 1 (from r=0 to r=xcccrc)
!!
!! WARNINGS
!! the fifth derivative is not yet delivered.
!!
!! PARENTS
!!      psp1in
!!
!! CHILDREN
!!      spline
!!
!! SOURCE

subroutine psp4cc(fchrg,n1xccc,xccc1d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1xccc
 real(dp),intent(in) :: fchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1xccc,ider
 real(dp),parameter :: a10=-0.1156854803757563d5,a12=+0.2371534625455588d5
 real(dp),parameter :: a14=-0.3138755797827918d5,a16=+0.2582842713241039d5
 real(dp),parameter :: a18=-0.1200356429115204d5,a20=+0.2405099057118771d4
 real(dp),parameter :: a2=-0.8480751097855989d1,a4=+0.9684600878284791d2
 real(dp),parameter :: a6=-0.7490894651588015d3,a8=+0.3670890998130434d4
 real(dp) :: der1,dern,factor
 character(len=500) :: message
!arrays
 real(dp),allocatable :: ff(:),ff2(:),work(:),xx(:)
 real(dp) :: x

! *************************************************************************

 ABI_ALLOCATE(ff,(n1xccc))
 ABI_ALLOCATE(ff2,(n1xccc))
 ABI_ALLOCATE(work,(n1xccc))
 ABI_ALLOCATE(xx,(n1xccc))


 if(n1xccc > 1)then
   factor=1.0d0/dble(n1xccc-1)
   do i1xccc=1,n1xccc
     xx(i1xccc)=(i1xccc-1)*factor
   end do
 else
   write(message, '(a,i0)' )'  n1xccc should larger than 1, while it is n1xccc=',n1xccc
   MSG_BUG(message)
 end if

!Initialization, to avoid some problem with some compilers
 xccc1d(1,:)=zero ; xccc1d(n1xccc,:)=zero

!Take care of each derivative separately
 do ider=0,2

   if(ider==0)then
!    Generate spline fitting for the function gg
     do i1xccc=1,n1xccc
!      ff(i1xccc)=fchrg*gg(xx(i1xccc))
       ff(i1xccc)=fchrg*gg_psp4(xx(i1xccc))
     end do
!    Complete with derivatives at end points
     der1=0.0d0
!    dern=fchrg*gp(1.0d0)
     dern=fchrg*gp_psp4(1.0d0)
   else if(ider==1)then
!    Generate spline fitting for the function gp
     do i1xccc=1,n1xccc
!      ff(i1xccc)=fchrg*gp(xx(i1xccc))
       ff(i1xccc)=fchrg*gp_psp4(xx(i1xccc))
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
       x=xx(i1xccc)
       ff(i1xccc)=fchrg*(gpp_1_psp4(x)+gpp_2_psp4(x)+gpp_3_psp4(x))
!      ff(i1xccc)=fchrg*gpp(xx(i1xccc))
     end do
!    Complete with derivatives of end points
     der1=xccc1d(1,ider+2)
     dern=xccc1d(n1xccc,ider+2)
   end if

!  Produce second derivative numerically, for use with splines
   call spline(xx,ff,n1xccc,der1,dern,ff2)
   xccc1d(:,ider+1)=ff(:)
   xccc1d(:,ider+3)=ff2(:)
 end do

 xccc1d(:,6)=zero

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(xx)

!DEBUG
!write(std_out,*)' psp1cc : output of core charge density and derivatives '
!write(std_out,*)'   xx          gg           gp  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   xx          gpp          gg2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   xx          gp2          gpp2  '
!do i1xccc=1,n1xccc
!write(std_out,'(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 contains

   function gg_psp4(x)
!Expression of 7 May 1992
   real(dp) :: gg_psp4
   real(dp),intent(in) :: x
   gg_psp4=(1.d0+x**2*(a2 +x**2*(a4 +x**2*(a6 +x**2*(a8 + &
&   x**2*(a10+x**2*(a12+x**2*(a14+x**2*(a16+ &
&   x**2*(a18+x**2*(a20)))))))))))          *(1.0d0-x**2)**2
 end function gg_psp4

   function gp_psp4(x)
!gp(x) is the derivative of gg(x) wrt x
   real(dp) :: gp_psp4
   real(dp),intent(in) :: x
   gp_psp4=2.d0*x*((a2+x**2*(2.d0*a4+x**2*(3.d0*a6+x**2*(              &
&   4.d0*a8+x**2*(5.d0*a10+x**2*(6.d0*a12+x**2*(                     &
&   7.d0*a14+x**2*(8.d0*a16+x**2*(9.d0*a18+x**2*(10.d0*a20))))))))))*&
&   (1.d0-x**2)**2                                                &
&   -2.0d0*(1.d0+x**2*(a2 +x**2*(a4 +x**2*(a6 +x**2*(a8 +            &
&   x**2*(a10+x**2*(a12+x**2*(a14+x**2*(a16+            &
&   x**2*(a18+x**2*a20))))))))))        *(1.0d0-x**2) )
 end function gp_psp4

   function gpp_1_psp4(x)
!gpp(x) is the second derivative of gg(x) wrt x
   real(dp) :: gpp_1_psp4
   real(dp),intent(in) :: x
   gpp_1_psp4= ( 2.d0*a4+ x**2*(3.d0*2.d0*a6 +x**2*(               &
&   4.d0*3.d0*a8+ x**2*(5.d0*4.d0*a10+x**2*(               &
&   6.d0*5.d0*a12+x**2*(7.d0*6.d0*a14+x**2*(               &
&   8.d0*7.d0*a16+x**2*(9.d0*8.d0*a18+x**2*(               &
&   10.d0*9.d0*a20)                                        &
&   ))))))))*(2.d0*x*(1.d0-x**2))**2
 end function gpp_1_psp4

   function gpp_2_psp4(x)

   real(dp) :: gpp_2_psp4
   real(dp),intent(in) :: x
   gpp_2_psp4=(a2+x**2*(2.d0*a4+x**2*(3.d0*a6+x**2*(                 &
&   4.d0*a8 +x**2*(5.d0*a10+x**2*(6.d0*a12+x**2*(          &
&   7.d0*a14+x**2*(8.d0*a16+x**2*(9.d0*a18+x**2*(          &
&   10.d0*a20)                                             &
&   )))))))))*(1.d0-x**2)*2*(1.d0-9.d0*x**2)
 end function gpp_2_psp4

   function gpp_3_psp4(x)

   real(dp) :: gpp_3_psp4
   real(dp),intent(in) :: x
   gpp_3_psp4=(1.d0+x**2*(a2 +x**2*(a4 +x**2*(a6 +x**2*(a8 +         &
&   x**2*(a10+x**2*(a12+x**2*(a14+x**2*(a16+         &
&   x**2*(a18+x**2*a20                               &
&   ))))))))))*(1.0d0-3.d0*x**2)*(-4.d0)
 end function gpp_3_psp4

end subroutine psp4cc
!!***

end module m_psp1
!!***

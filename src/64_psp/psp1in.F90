!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp1in
!! NAME
!! psp1in
!!
!! FUNCTION
!! Initialize pspcod=1 or 4 pseudopotential (Teter format):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, FrD, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&
&                  e990,e999,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mqgrid,nproj,n1xccc,pspcod,&
&                  qchrg,qgrid,rcpsp,rms,useylm,vlspl,xcccrc,xccc1d,&
&                  zion,znucl)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp1in'
 use interfaces_14_hidewrite
 use interfaces_64_psp, except_this_one => psp1in
!End of the abilint section

 implicit none

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
&   'Action : check your pseudopotential input file.'
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
 write(message, '(3f20.14,t64,a)' ) rchrg,fchrg,qchrg,&
& 'rchrg,fchrg,qchrg'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Generate core charge function and derivatives, if needed
 if(fchrg>1.0d-15)then
   if(pspcod==1)then
     call psp1cc(fchrg,n1xccc,xccc1d)
!    The core charge function for pspcod=1
!    becomes zero beyond 3*rchrg only. Thus xcccrc must be set
!    equal to 3*rchrg .
     xcccrc=3*rchrg
   else if(pspcod==4)then
     call psp4cc(fchrg,n1xccc,xccc1d)
!    For pspcod=4, the core charge cut off exactly beyond rchrg
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
!  DEBUG
!  write(std_out,*) 'END OF READING PSP',ll,'OK'
!  ENDDEBUG

 end do

!Copy appropriate nonlocal psp for use as local one
 vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )

!DEBUG
!write(std_out,*) 'VLOC=',vloc(1),vloc(2),vloc(3)
!write(std_out,*) 'VLOC=',vloc(4),vloc(5),vloc(6)
!ENDDEBUG

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

!DEBUG
!write(std_out,*) 'RADIAL GRID CREATED'
!ENDDEBUG

!here compute sin(r(:)*dq) and cos(r(:)*dq)
!NOTE : also invert dr !!
 dq2pi=2.0d0*pi*dq
 do ii=1,mmax
   arg=dq2pi*rad(ii)
   drad(ii)=1.0d0/drad(ii)
   wksincos(ii,1,1)=sin(arg)
   wksincos(ii,2,1)=cos(arg)
 end do

!(3)Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).
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

!DEBUG
!write(std_out,*)' psp1in : before nonlocal corrections '
!write(std_out,*)' psp1in : lloc, lmax = ',lloc,lmax
!if(.true.)stop
!ENDDEBUG

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then

   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 else

!  Proceed to make Kleinman-Bylander form factors for
!  each l up to lmax

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
&         'Action : check your pseudopotential input file.'
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

!  DEBUG
!  write(std_out,*)' psp1in : lmax,lloc=',lmax,lloc
!  ENDDEBUG
   ABI_ALLOCATE(ekb_tmp,(mpsang,2))
   ABI_ALLOCATE(ffspl_tmp,(mqgrid,2,nlmax+1,2))

   call psp1nl(drad,ekb_tmp(:,1),ffspl_tmp(:,:,:,1),lloc,&
&   nlmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll,wksincos)

!  Read second wavefunction for second projection operator
!  (only read cases where nproj(ll)=2)
!  --also find highest l for which nproj(l)=2

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
&         'Action : check your pseudopotential input file.'
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
!      DEBUG
!      write(std_out,*)' psp1in : lmnmax,ii,indlmn(1,ii)=',lmnmax,ii,indlmn(1,ii)
!      ENDDEBUG
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

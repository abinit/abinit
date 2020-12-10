!!****m* ABINIT/m_psp5
!! NAME
!!  m_psp5
!!
!! FUNCTION
!! Initialize pspcod=5 ("Phoney pseudopotentials" with Hamman grid):
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, FrD, FJ, MT)
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

module m_psp5

 use defs_basis
 use m_splines
 use m_errors
 use m_abicore

 use m_psptk,           only : psp1cc, psp5lo, psp5nl

 implicit none

 private
!!***

 public :: psp5in
!!***

contains
!!***

!!****f* ABINIT/psp5in
!! NAME
!! psp5in
!!
!! FUNCTION
!! Initialize pspcod=5 ("Phoney pseudopotentials" with Hamman grid):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! INPUTS
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pspso= spin orbit signal
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=nuclear number of atom as specified in psp file
!!
!!  OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!             if any, spin-orbit components begin at l=mpsang+1
!!  ekb1(mpssoang)= Kleinman-Bylander energy from the psp file, for iproj=1
!!  ekb2(mpssoang)= Kleinman-Bylander energy from the psp file, for iproj=2
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r} dr]$ (hartree)
!!  epspsp(mpssoang)=values of epsatm for different angular momenta, from the psp file
!!  e990(mpssoang)=ecut at which 0.99 of the kinetic energy is recovered
!!  e999(mpssoang)=ecut at which 0.999 of the kinetic energy is recovered
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpssoang)=number of projection functions for each angular momentum
!!  qchrg is the total (integrated) core charge
!!  rcpsp(mpssoang)=cut-off radius for each angular momentum
!!  rms(mpssoang)=root mean square of the KB psp
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      psp1cc,psp5lo,psp5nl,spline,wrtout
!!
!! SOURCE

subroutine psp5in(ekb,ekb1,ekb2,epsatm,epspsp,e990,e999,ffspl,indlmn,&
&                  lloc,lmax,lmnmax,lnmax,mmax,mpsang,mpssoang,mqgrid,&
&                  nproj,n1xccc,pspso,qchrg,qgrid,rcpsp,rms,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mmax,mpsang,mpssoang,mqgrid
 integer,intent(in) :: n1xccc,pspso,useylm
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
!arrays
 integer,intent(out) :: indlmn(6,lmnmax) !vz_i
 integer,intent(inout) :: nproj(mpssoang) !vz_i
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: e990(mpssoang),e999(mpssoang),ekb(lnmax)
 real(dp),intent(out) :: ekb1(mpssoang),ekb2(mpssoang),epspsp(mpssoang)
 real(dp),intent(out) :: rcpsp(mpssoang),rms(mpssoang) !vz_i
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax) !vz_i
 real(dp),intent(out) :: vlspl(mqgrid,2) !vz_i
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ii,iln,index,ipsang,kk,lhigh,ll,mm,mproj,nn,nso,pspso0
 real(dp) :: al,fchrg,r1,rchrg,yp1,ypn
 logical :: test
 character(len=500) :: message,errmsg
!arrays
 real(dp),allocatable :: ekb_so(:),ekb_sr(:),ekb_tmp(:,:),ffspl_so(:,:,:)
 real(dp),allocatable :: ffspl_sr(:,:,:),ffspl_tmp(:,:,:,:),rad(:),vloc(:)
 real(dp),allocatable :: vpspll(:,:),vpspll_so(:,:),wfll(:,:),wfll_so(:,:)
 real(dp),allocatable :: work_space(:),work_spl(:)

! ***************************************************************************

!File format of formatted Phoney psp input (the 3 first lines
!have already been read in calling -pspatm- routine) :

!(1) title (character) line
!(2) znucl,zion,pspdat
!(3) pspcod,pspxc,lmax,lloc,mmax,r2well
!(4) r1,al,pspso
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

!Read fourth line of the file ; parameter pspso is optional (but
!this is not treated correctly by all machines - problems with SGI )
 pspso0=1
 read (tmp_unit,fmt=*,err=50,end=50) r1,al,pspso0
 50 continue

 if(pspso0/=1 .and. pspso0/=2)then
   write(message, '(3a,i0,2a)' )&
&   'Problem reading the fourth line of pseudopotential file.',ch10,&
&   'The parameter pspso should be 1 or 2, but it is pspso= ',pspso0,ch10,&
&   'Action: check your pseudopotential input file.'
   ABI_ERROR(message)
 end if

 write(message, '(2es16.6,t47,a)' ) r1,al,'r1 and al (Hamman grid)'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 if (pspso0/=1) then
   write(message,'(a)') ' Pseudopotential is in spin-orbit format '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 if (pspso/=0.and.pspso0==1) then
   write(message, '(a,a,a,a,a)' )&
&   'The treatment of spin-orbit interaction is required (pspso/=0)',ch10,&
&   'but pseudopotential file format cannot contain spin-orbit information !',ch10,&
&   'Action: check your pseudopotential input file.'
   ABI_ERROR(message)
 end if

 nso=1;if (pspso0/=1) nso=2
 do ipsang=1,(nso*lmax)+1
   read (tmp_unit,*, err=10, iomsg=errmsg) ll,e990(ipsang),e999(ipsang),nproj(ipsang),rcpsp(ipsang)
   write(message, '(i5,2f8.3,i5,f12.7,t47,a)' ) &
&   ll,e990(ipsang),e999(ipsang),nproj(ipsang),rcpsp(ipsang),'l,e99.0,e99.9,nproj,rcpsp'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   read (tmp_unit,*, err=10, iomsg=errmsg) rms(ipsang),ekb1(ipsang),ekb2(ipsang),epspsp(ipsang)
   write(message, '(4f13.8,t55,a)' ) &
&   rms(ipsang),ekb1(ipsang),ekb2(ipsang),epspsp(ipsang),'   rms, ekb1, ekb2, epsatm'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do

!If pspso/=0 and nproj/=2, forces nproj to be 2
!(for compatibility with old psp-file format)
 if (pspso/=0.and.lmax>0) then
   test=.false.
   do ipsang=1,(nso*lmax)+1
     if (ipsang>1.and.nproj(ipsang)/=2) then
       test=.true.;nproj(ipsang)=2
     end if
   end do
   if (test) then
     write(message, '(a,a,a,a,a)' )&
&     'Pseudopotential file is spin-orbit (pspso=2)',ch10,&
&     'and number of projector for l/=0 is not 2 !',ch10,&
&     'It has been forced to 2.'
     call wrtout(std_out,message,'COLL')
     ABI_WARNING(message)
   end if
 end if

!mproj=maxval(nproj(1:lmax+1))
!mjv 10/2008: I believe this is correct. Perhaps unnecessary if the normal
!projectors are always more numerous, but should be conservative anyway with
!maxval.
 mproj=maxval(nproj)
 index=0;iln=0;indlmn(:,:)=0
 do nn=1,nso
   do ipsang=1+(nn-1)*(lmax+1),nn*lmax+1
     if (nproj(ipsang)>0) then
       ll=ipsang-(nn-1)*lmax-1
       do kk=1,nproj(ipsang)
         iln=iln+1
         do mm=1,2*ll*useylm+1
           index=index+1
           indlmn(1,index)=ll
           indlmn(2,index)=mm-ll*useylm-1
           indlmn(3,index)=kk
           indlmn(4,index)=ll*ll+(1-useylm)*ll+mm
           indlmn(5,index)=iln
           indlmn(6,index)=nn
         end do
       end do
     end if
   end do
 end do

 read (tmp_unit,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
 write(message, '(3f20.14,t64,a)' ) rchrg,fchrg,qchrg,'rchrg,fchrg,qchrg'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Generate core charge function and derivatives, if needed
 xcccrc=zero
 if(n1xccc>0)then
!  Use the revised expression of 5 Nov 1992, also used for format=1.
   call psp1cc(fchrg,n1xccc,xccc1d)
   xcccrc=3*rchrg
 end if

!--------------------------------------------------------------------
!Will now proceed at the reading of pots and wfs, as well as their treatment

!vpspll(:,1),...,vpspll(:,4)=nonlocal pseudopotentials
!vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
!rad(:)=radial grid r(i)
!wfll(:,1),...,wfll(:,4)=reference config. wavefunctions
 ABI_MALLOC(vloc,(mmax))
 ABI_MALLOC(vpspll,(mmax,mpsang))

!(1) Read atomic pseudopotential for each l, filling up array vpspll
!Note: put each l into vpspll(:,l+1)

 if (pspso==0) then

!  --NON SPIN-ORBIT
   do ipsang=1,lmax+1
     read (tmp_unit,*, err=10, iomsg=errmsg) ll
     read (tmp_unit,*, err=10, iomsg=errmsg) (vpspll(ii,ipsang),ii=1,mmax)
!    write(std_out,*) 'END OF READING PSP',ll,'OK'
   end do
 else

!  --SPIN-ORBIT
   ABI_MALLOC(vpspll_so,(mmax,mpsang))
   read (tmp_unit,*, err=10, iomsg=errmsg) ll
   read (tmp_unit,*, err=10, iomsg=errmsg) (vpspll(ii,1),ii=1,mmax)
   vpspll_so(:,1)=0.0d0
   do ipsang=2,lmax+1
     read (tmp_unit,*, err=10, iomsg=errmsg) ll
     read (tmp_unit,*, err=10, iomsg=errmsg) (vpspll(ii,ipsang),ii=1,mmax)
     read (tmp_unit,*, err=10, iomsg=errmsg) ll
     read (tmp_unit,*, err=10, iomsg=errmsg) (vpspll_so(ii,ipsang),ii=1,mmax)
   end do
 end if

!Copy appropriate nonlocal psp for use as local one
 if (pspso==0) then
   vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )
 else
   if(lloc<=0) then
     vloc( 1:mmax ) = vpspll( 1:mmax , -lloc+1 )
   else
     vloc( 1:mmax ) = vpspll_so( 1:mmax , lloc+1 )
   end if
 end if
!DEBUG
!write(std_out,*) 'VLOC=',vloc(1),vloc(2),vloc(3)
!write(std_out,*) 'VLOC=',vloc(4),vloc(5),vloc(6)
!ENDDEBUG


!(2) Create radial grid, and associated quantities

!Now compute Hamman Grid
 ABI_MALLOC(rad,(mmax))
 do ii=1,mmax
   rad (ii)=r1*exp(dble(ii-1)*al)
 end do
!DEBUG
!write(std_out,*) 'HAMMAN RADIAL GRID r1 and al',r1,al
!write(std_out,*) 'rad(1)=',rad(1)
!write(std_out,*) 'rad(10)=',rad(10)
!write(std_out,*) 'rad(100)=',rad(100)
!ENDDEBUG


!(3)Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).

 call psp5lo(al,epsatm,mmax,mqgrid,qgrid,&
& vlspl(:,1),rad,vloc,yp1,ypn,zion)


!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_MALLOC(work_space,(mqgrid))
 ABI_MALLOC(work_spl,(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)

 ABI_FREE(work_space)
 ABI_FREE(work_spl)

!(4)Take care of non-local part

!DEBUG
!write(std_out,*)' psp5in : before nonlocal corrections '
!write(std_out,*)' psp5in : lloc, lmax = ',lloc,lmax
!ENDDEBUG

!Zero out all Kleinman-Bylander energies to initialize
 ekb(:)=0.0d0

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then

   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 else

!  Proceed to make Kleinman-Bylander form factors for
!  each l up to lmax

!  Read wavefunctions for each l up to lmax
   ABI_MALLOC( wfll,(mmax,mpsang))
!  -----------------------------------------------------------------

   if (pspso==0) then

!    --NON SPIN-ORBIT
     do ipsang=1,lmax+1
       if (nproj(ipsang)/=0) then
         read (tmp_unit,*, err=10, iomsg=errmsg) ll
         if (ipsang/=ll+1) then
           write(message, '(a,a,a,a,a,a,2i6,a,a)' )&
&           'Pseudopotential input file does not have',ch10,&
&           'angular momenta in order expected for first projection',&
&           'operator.',ch10,' Values are ',ipsang-1,ll,ch10,&
&           'Action: check your pseudopotential input file.'
           ABI_ERROR(message)
         end if
         read (tmp_unit,*, err=10, iomsg=errmsg) wfll(:,ipsang)
       else
         wfll(:,ipsang)=0.0d0
       end if
     end do
   else

!    --SPIN-ORBIT
     ABI_MALLOC(wfll_so,(mmax,mpsang))
     if (nproj(1)/=0) then
       read (tmp_unit,*,err=10,iomsg=errmsg) ll
       read (tmp_unit,*,err=10,iomsg=errmsg) wfll(:,1)
     else
       wfll(:,1)=0.0d0
     end if
     wfll_so(:,1)=0.0d0
     do ipsang=2,lmax+1
       if (nproj(ipsang)/=0) then
         read (tmp_unit,*,err=10,iomsg=errmsg) ll
         read (tmp_unit,*,err=10,iomsg=errmsg) wfll(:,ipsang)
         read (tmp_unit,*,err=10,iomsg=errmsg) ll
         read (tmp_unit,*,err=10,iomsg=errmsg) wfll_so(:,ipsang)
       else
         wfll(:,ipsang)=0.0d0
         wfll_so(:,ipsang)=0.0d0
       end if
     end do
   end if

!  ----------------------------------------------------------------------
!  Compute KB form factors and fit splines
   ABI_MALLOC(ekb_tmp,(mpssoang,max(nso,mproj)))
   ABI_MALLOC(ffspl_tmp,(mqgrid,2,mpssoang,max(nso,mproj)))
   ekb_tmp(:,:)=0.d0

   ABI_MALLOC(ekb_sr,(mpsang))
   ABI_MALLOC(ffspl_sr,(mqgrid,2,mpsang))
   call psp5nl(al,ekb_sr(:),ffspl_sr(:,:,:),lmax,mmax,mpsang,mqgrid,&
&   qgrid,rad,vloc,vpspll,wfll)
   ekb_tmp(1:mpsang,1)=ekb_sr(1:mpsang)
   ffspl_tmp(:,:,1:mpsang,1)=ffspl_sr(:,:,1:mpsang)

   if (pspso/=0) then
     ABI_MALLOC(ekb_so,(mpsang))
     ABI_MALLOC(ffspl_so,(mqgrid,2,mpsang))
     call psp5nl(al,ekb_so,ffspl_so,lmax,mmax,mpsang,mqgrid,&
&     qgrid,rad,vloc,vpspll_so,wfll_so)
     ekb_tmp(mpsang+1:mpssoang,1)=ekb_so(2:mpsang)
     do ipsang=2,lmax+1
       if((ekb_sr(ipsang)*ekb_so(ipsang))<0.0) then
         ABI_ERROR('BIG PROBLEM WITH THE SPIN ORBIT IN PSP5NL')
       end if
     end do

     if(lloc<0) ekb_sr(-lloc+1)=ekb_so(-lloc+1)
     if(lloc<0) ekb_tmp(-lloc+1,1)=ekb_tmp(-lloc+1+lmax,1)
     if(lloc>0) ekb_so(lloc+1)=ekb_sr(lloc+1)
     if(lloc>0) ekb_tmp(lmax+lloc+1,1)=ekb_tmp(lloc+1,1)
     do ipsang=1,mpssoang
       if(ekb_tmp(ipsang,1)>0) ekb_tmp(ipsang,1)= 1.d0
       if(ekb_tmp(ipsang,1)<0) ekb_tmp(ipsang,1)=-1.d0
     end do

!    v_ion is calculated in ffspl_tmp(:,:,1:mpsang,1) and v_so in
!    ffspl_tmp(:,:,mpsang+1:mpssoang,1) taking into account sqrt(ekb)
     do i1=1,mqgrid
       do i2=1,2
         ffspl_tmp(i1,i2,1,1)=ffspl_sr(i1,i2,1)*sqrt(abs(ekb_sr(1)))
         do ipsang=2,mpsang
           ffspl_tmp(i1,i2,ipsang,1)=((ffspl_sr(i1,i2,ipsang)*&
&           sqrt(abs(ekb_sr(ipsang)))*(ipsang-1))+&
&           (ffspl_so(i1,i2,ipsang)*&
&           sqrt(abs(ekb_so(ipsang)))*(ipsang)))&
&           /(2.d0*ipsang-1)
           ffspl_tmp(i1,i2,mpsang+ipsang-1,1)=(-ffspl_sr(i1,i2,ipsang)*&
&           sqrt(abs(ekb_sr(ipsang)))+&
&           ffspl_so(i1,i2,ipsang)*&
&           sqrt(abs(ekb_so(ipsang))))*2.d0&
&           /(2.d0*ipsang-1)
         end do
       end do
     end do
     ABI_FREE(ekb_so)
     ABI_FREE(ffspl_so)
     ABI_FREE(vpspll_so)
     ABI_FREE(wfll_so)

!    The non local contribution is written as quadratic form of the vector
!    V=(v_ion,v_so)
!    t_V (Q1+Q2 L.S) V
!    with Q1= (1      0   )    et   Q2=(0     1 )
!    (0  l(l+1)/4)            (1   -1/2)
!    The LS independent part is already diagonal. V is therefore built
!    putting v_so in the second projector of ffspl for the non spin-orbit
!    part and taking the eigenvalues of Q1 as new ekb (apart the sign)
     do ipsang=2,mpsang
       do i1=1,mqgrid
         do i2=1,2
           ffspl_tmp(i1,i2,ipsang,2)= ffspl_tmp(i1,i2,mpsang+ipsang-1,1)
         end do
       end do
       ekb_tmp(ipsang,2)=ekb_tmp(mpsang+ipsang-1,1)*ipsang*(ipsang-1)*0.25d0
     end do

!    For the spin orbit part, after diagonalisation of Q2, the eigenvectors
!    are: ((1-sqrt(17))/4  , 1) and ((1+sqrt(17))/4 ,1)
!    The passage matrix is therefore P=((1-sqrt(17))/4  (1+sqrt(17))/4)
!    (    1                 1       )
!    t_P*Q2*P=( -sqrt(17)/2   0    )
!    ( 0        sqrt(17)/2)
!    The diagonal values are the new ekb and the new ffspl are
!    P^-1 (v_ion)
!    (v_so )
     do ipsang=2,mpsang
       do i1=1,mqgrid
         do i2=1,2
           ffspl_tmp(i1,i2,mpsang+ipsang-1,1)=-2.d0/sqrt(17.d0)*&
&           (ffspl_tmp(i1,i2,ipsang,1)-&
&           ((sqrt(17.d0)+1)*0.25d0)*&
           ffspl_tmp(i1,i2,ipsang,2))
           ffspl_tmp(i1,i2,mpsang+ipsang-1,2)=2.d0/sqrt(17.d0)*&
&           (ffspl_tmp(i1,i2,ipsang,1)+&
&           ((sqrt(17.d0)-1)*0.25d0)*&
&           ffspl_tmp(i1,i2,ipsang,2))
         end do
       end do
       ekb_tmp(mpsang+ipsang-1,1)=-(sqrt(17.d0)*0.5d0)*ekb_tmp(ipsang,1)
       ekb_tmp(mpsang+ipsang-1,2)= (sqrt(17.d0)*0.5d0)*ekb_tmp(ipsang,1)
     end do

   end if

   ABI_FREE(ekb_sr)
   ABI_FREE(ffspl_sr)

!  FJ WARNING : No spin orbit if nproj>1
   if (pspso==0) then

!    Read second wavefunction for second projection operator
!    (only read cases where nproj(ll)=2)
!    --also find highest l for which nproj(l)=2
     lhigh=-1
     do ipsang=1,min(lmax+1,mpsang)
       if (nproj(ipsang)==2) then
         lhigh=ipsang-1
         read (tmp_unit,*, err=10, iomsg=errmsg) ll
         if (ipsang/=ll+1) then
           write(message, '(a,a,a,a,a,a,2i6,a,a)' )&
&           'Pseudopotential input file does not have',ch10,&
&           'angular momenta in order expected for second projection',&
&           'operator.',ch10,' Values are ',ipsang-1,ll,ch10,&
&           'Action: check your pseudopotential input file.'
           ABI_ERROR(message)
         end if
         read (tmp_unit,*, err=10, iomsg=errmsg) wfll(:,ipsang)
!        DEBUG
!        write(std_out,*) 'WF second',ipsang-1,wfll(1,ipsang),wfll(2,ipsang),wfll(3,ipsang)
!        ENDDEBUG
       else
         wfll(:,ipsang)=0.0d0
       end if

     end do

!    Compute KB form factors and fit splines for second wf if any
     if (lhigh>-1) then
       call psp5nl(al,ekb_tmp(:,2),ffspl_tmp(:,:,:,2),lmax,&
&       mmax,mpsang,mqgrid,qgrid,rad,vloc,vpspll,wfll)
     end if

   end if

!  Convert ekb and ffspl
   iln=0
   do ii=1,lmnmax
     kk=indlmn(5,ii)
     if (kk>iln) then
       iln=kk
       ll=indlmn(1,ii);nn=indlmn(3,ii)
       if (indlmn(6,ii)==1) then
         ekb(kk)=ekb_tmp(1+ll,nn)
         ffspl(:,:,kk)=ffspl_tmp(:,:,1+ll,nn)
       else
         ekb(kk)=ekb_tmp(mpsang+ll,nn)
         ffspl(:,:,kk)=ffspl_tmp(:,:,mpsang+ll,nn)
       end if
     end if
   end do

   ABI_FREE(ekb_tmp)
   ABI_FREE(ffspl_tmp)
   ABI_FREE(wfll)

!  end of if concerning lloc
 end if

 ABI_FREE(vpspll)
 ABI_FREE(rad)
 ABI_FREE(vloc)

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp5in
!!***

end module m_psp5
!!***

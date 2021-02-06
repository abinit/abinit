!!****m* ABINIT/m_mkffnl
!! NAME
!!  m_mkffnl
!!
!! FUNCTION
!! Make FFNL, nonlocal form factors, for each type of atom up to ntypat
!! and for each angular momentum.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MT, DRH)
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

module m_mkffnl

 use defs_basis
 use m_abicore
 use m_errors
 use m_splines
 use m_xmpi

 use m_time,     only : timab
 use m_kg,       only : mkkin
 use m_sort,     only : sort_dp

 implicit none

 private
!!***

 public :: mkffnl
!!***

contains
!!***

!!****f* ABINIT/mkffnl
!! NAME
!! mkffnl
!!
!! FUNCTION
!! Make FFNL, nonlocal form factors, for each type of atom up to ntypat
!! and for each angular momentum.
!! When Legendre polynomials are used in the application of the
!!   nonlocal operator, FFNLs depend on (l,n) components; in this
!!   case, form factors are real and divided by |k+G|^l;
!! When spherical harmonics are used, FFNLs depend on (l,m,n)
!!   components; in this case, form factors are multiplied by Ylm(k+G).
!!
!! INPUTS
!!  dimekb=second dimension of ekb (see ekb)
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  ekb(dimekb,ntypat*(1-usepaw))=(Real) Kleinman-Bylander energies (hartree)
!!                                ->NORM-CONSERVING PSPS ONLY
!!  ffspl(mqgrid,2,lnmax,ntypat)=form factors and spline fit to 2nd derivative
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  ider=0=>no derivative wanted; 1=>1st derivative wanted; 2=>1st and 2nd derivatives wanted
!!  idir=ONLY WHEN YLMs ARE USED:
!!       When 1st derivative has to be computed:  (see more info below)
!!       - Determine the direction(s) of the derivatives(s)
!!       - Determine the set of coordinates (reduced or cartesians)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
!!                                                     or i=lmn (if useylm=1)
!!  kg(3,npw)=integer coordinates of planewaves in basis sphere for this k point.
!!  kpg(npw,nkpg)= (k+G) components (only if useylm=1)
!!  kpt(3)=reduced coordinates of k point
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=size of q (or |G|) grid for f(q)
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere
!!  ntypat=number of types of atoms
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  pspso(ntypat)=spin-orbit characteristics for each atom type (1, 2, or 3)
!!  qgrid(mqgrid)=uniform grid of q values from 0 to qmax
!!  rmet(3,3)=real space metric (bohr**2)
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  ylm   (npw,mpsang*mpsang*useylm)=real spherical harmonics for each G and k point
!!  ylm_gr(npw,3,mpsang*mpsang*useylm)=gradients of real spherical harmonics wrt (k+G)
!! [comm]=MPI communicator. Default: xmpi_comm_self.
!!
!! OUTPUT
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=described below
!! [request]=Used in conjunction with [comm] to perform non-blocking xmpi_isum_ip. Client code must
!!  wait on request before using ffnl. If not present, blocking API is used.
!!
!! NOTES
!!  Uses spline fit ffspl provided by Numerical Recipes spline subroutine.
!!  Form factor $f_l(q)$ is defined by
!!   \begin{equation}
!!  \textrm{f}_l(q)=\frac{1}{dvrms} \int_0^\infty [j_l(2 \pi r q) u_l(r) dV(r) r dr]
!!   \end{equation}
!!   where u_l(r)=reference state wavefunction, dV(r)=nonlocal psp
!!   correction, j_l(arg)=spherical Bessel function for angular momentum l,
!!   and
!!   \begin{equation}
!!    \textrm{dvrms} =  \int_0^\infty [(u_l(r) dV(r))^2 dr])^{1/2}
!!   \end{equation}
!!   which is square root of mean square dV, i.e.
!!     $ (\langle (dV)^2 \rangle)^{1/2} $ .
!!   This routine is passed f_l(q) in spline form in the array ffspl and then
!!   constructs the values of $f_l(q)$ on the relevant (k+G) in array ffnl.
!!   The evaluation of the integrals defining ffspl was done in mkkbff.
!!
!!  Delivers the following (for each atom type t, or itypat):
!!   --------------------------
!!   Using Legendre polynomials in the application of nl operator:
!!     ffnl are real.
!!     ffnl(ig,1,(l,0,n),itypat) $= f_ln(k+G)/|k+G|^l $
!!     === if ider>=1
!!       ffnl(ig,2,(l,0,n),itypat) $=(fprime_ln(k+G)-l*f_ln(k+G)/|k+G|)/|k+G|^(l+1) $
!!     === if ider=2
!!       ffnl(ig,3,(l,0,n),itypat) $=(fprimeprime_ln(k+G)-(2l+1)*fprime_ln(k+G)/|k+G|
!!                                   +l(l+2)*f_ln(k+G)/|k+G|**2)/|k+G|^(l+2)
!!   --------------------------
!!   Using spherical harmonics in the application of nl operator:
!!     ffnl are real (we use REAL spherical harmonics).
!!     ffnl(ig,1,(l,m,n),itypat) = ffnl_1
!!                              $= f_ln(k+G) * Y_lm(k+G) $
!!     === if ider>=1
!!     --if (idir==0)
!!       ffnl(ig,1+i,(l,m,n),itypat) = dffnl_i = 3 reduced coord. of d(ffnl_1)/dK^cart
!!         $= fprime_ln(k+G).Y_lm(k+G).(k+G)^red_i/|k+G|+f_ln(k+G).(dY_lm/dK^cart)^red_i $
!!         for i=1..3
!!     --if (1<=idir<=3)
!!       ffnl(ig,2,(l,m,n),itypat)= cart. coordinate idir of d(ffnl_1)/dK^red
!!                                = Sum_(mu,nu) [ Gprim(mu,idir) Gprim(mu,nu) dffnl_nu ]
!!     --if (idir==4)
!!       ffnl(ig,1+i,(l,m,n),itypat)= 3 cart. coordinates of d(ffnl_1)/dK^red
!!                                  = Sum_(mu,nu) [ Gprim(mu,i) Gprim(mu,nu) dffnl_nu ]
!!     --if (-6<idir<-1)
!!       ffnl(ig,2,(l,m,n),itypat)=1/2 [d(ffnl)/dK^cart_mu K^cart_nu + d(ffnl)/dK^cart_nu K^cart_mu]
!!                                with d(ffnl)/dK^cart_i = Sum_nu [ Gprim(nu,i) dffnl_nu ]
!!                                for |idir|->(mu,nu) (1->11,2->22,3->33,4->32,5->31,6->21)
!!     --if (idir==-7)
!!       ffnl(ig,2:7,(l,m,n),itypat)=1/2 [d(ffnl)/dK^cart_mu K^cart_nu + d(ffnl)/dK^cart_nu K^cart_mu]
!!                                with d(ffnl)/dK^cart_i = Sum_nu [ Gprim(nu,i) dffnl_nu ]
!!                                for all (mu,nu) (6 independant terms)
!!     === if ider==2
!!     --if (idir==0)
!!       ffnl(ig,4+i,(l,m,n),itypat) = d2ffnl_mu,nu = 6 reduced coord. of d2(ffnl_1)/dK^cart.dK^cart
!!        for all i=(mu,nu) (6 independant terms)
!!     --if (idir==4)
!!       ffnl(ig,4+i,(l,m,n),itypat) = d2ffnl_i =6 cart. coordinates of d2(ffnl_1)/dK^red.dK^red
!!        for all i=(mu,nu) (6 independant terms)
!!        = Sum_(mu1,mu2,mu3,mu4) [ Gprim(mu1,mu) Gprim(mu2,nu) Gprim(mu1,mu3) Gprim(mu2,mu4) d2ffnl_mu3,mu4 ]
!!   --------------------------
!!
!!  1) l may be 0, 1, 2, or 3 in this version.
!!
!!  2) Norm-conserving psps : only FFNL for which ekb is not zero are calculated.
!!
!!  3) Each expression above approaches a constant as $|k+G| \rightarrow 0 $.
!!     In the cases where $|k+G|$ is in the denominator, there is always a
!!     factor of $(k+G)_mu$ multiplying the ffnl term where it is actually used,
!!     so that we may replace the ffnl term by any constant when $|k+G| = 0$.
!!     Below we replace 1/0 by 1/tol10, thus creating an arbitrary constant
!!     which will later be multiplied by 0.
!!
!! TODO
!!  Some parts can be rewritten with BLAS1 calls.
!!
!! PARENTS
!!      m_cgprj,m_d2frnl,m_dfpt_nstwf,m_dfpt_scfcv,m_dfptnl_pert,m_dft_energy
!!      m_fock_getghc,m_forstr,m_getgh1c,m_io_kss,m_ksdiago,m_nonlop_test
!!      m_orbmag,m_pead_nl_loop,m_vkbr,m_vtorho,m_wfd
!!
!! CHILDREN
!!      mkkin,splfit,timab
!!
!! SOURCE

subroutine mkffnl(dimekb, dimffnl, ekb, ffnl, ffspl, gmet, gprimd, ider, idir, indlmn, &
                   kg, kpg, kpt, lmnmax, lnmax, mpsang, mqgrid, nkpg, npw, ntypat, pspso, &
                   qgrid, rmet, usepaw, useylm, ylm, ylm_gr, &
                   comm, request) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimekb,dimffnl,ider,idir,lmnmax,lnmax,mpsang,mqgrid,nkpg
 integer,intent(in) :: npw,ntypat,usepaw,useylm
 integer,optional,intent(in) :: comm
 integer ABI_ASYNC, optional,intent(out):: request
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg(3,npw),pspso(ntypat)
 real(dp),intent(in) :: ekb(dimekb,ntypat*(1-usepaw))
 real(dp),intent(in) :: ffspl(mqgrid,2,lnmax,ntypat),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpg(npw,nkpg),kpt(3),qgrid(mqgrid),rmet(3,3)
 real(dp),intent(in) :: ylm(:,:),ylm_gr(:,:,:)
 real(dp),intent(out) :: ffnl(npw,dimffnl,lmnmax,ntypat)
 ! MG: Should be ABI_ASYNC due to optional non-Blocking API but NAG complains
 ! Error: m_d2frnl.F90, line 600: Array section FFNL_STR(:,:,:,:,MU) supplied for dummy FFNL (no. 4) of MKFFNL,
 ! the dummy is ASYNCHRONOUS but not assumed-shape
 ! so we declare request as ASYNCHRONOUS

!Local variables-------------------------------
!scalars
 integer :: ider_tmp,iffnl,ig,ig0,il,ilm,ilmn,iln,iln0,im,itypat,mu,mua,mub,nlmn,nu,nua,nub
 integer :: nprocs, my_rank, cnt, ierr
 real(dp),parameter :: renorm_factor=0.5d0/pi**2,tol_norm=tol10
 real(dp) :: ecut,ecutsm,effmass_free,fact,kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3,rmetab,yp1
 logical :: testnl=.false.
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 real(dp) :: rprimd(3,3),tsec(2)
 real(dp),allocatable :: dffnl_cart(:,:),dffnl_red(:,:),dffnl_tmp(:)
 real(dp),allocatable :: d2ffnl_cart(:,:),d2ffnl_red(:,:),d2ffnl_tmp(:)
 real(dp),allocatable :: kpgc(:,:),kpgn(:,:),kpgnorm(:),kpgnorm_inv(:),wk_ffnl1(:)
 real(dp),allocatable :: wk_ffnl2(:),wk_ffnl3(:),wk_ffspl(:,:)

! *************************************************************************

 ! Keep track of time spent in mkffnl
 call timab(16, 1, tsec)

 nprocs = 1; my_rank = 0
 if (present(comm)) then
   nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 end if

 ! Compatibility tests
 if (mpsang>4) then
   write(msg,'(a,i0,a,a)')&
   'Called with mpsang > 4, =',mpsang,ch10,&
   'This subroutine will not accept lmax+1 > 4.'
   ABI_BUG(msg)
 end if
 if (idir<-7.or.idir>4) then
   ABI_BUG('Called with idir<-6 or idir>4 !')
 end if
 if (useylm==0) then
   iffnl=1+ider
 else
   iffnl=1
   if (ider>=1) then
     if (idir==0) iffnl=iffnl+3
     if (idir/=0) iffnl=iffnl+1
     if (idir==4) iffnl=iffnl+2
     if (idir==-7) iffnl=iffnl+5
   end if
   if (ider==2) then
     if (idir==0) iffnl=iffnl+6
     if (idir==4) iffnl=iffnl+6
   end if
 end if
 if (iffnl/=dimffnl) then
   write(msg,'(2(a,i1),a,i2)') 'Incompatibility between ider, idir and dimffnl : ider = ',ider,&
                               ' idir = ',idir,' dimffnl = ',dimffnl
   ABI_BUG(msg)
 end if
 if (useylm==1) then
   ABI_CHECK(size(ylm,1)==npw,'BUG: wrong ylm size (1)')
   ABI_CHECK(size(ylm,2)==mpsang**2,'BUG: wrong ylm size (2)')
   if(ider>0)then
     ABI_CHECK(size(ylm_gr,1)==npw,'BUG: wrong ylm_gr size (1)')
     ABI_CHECK(size(ylm_gr,2)>=3+6*(ider/2),'BUG: wrong ylm_gr size (2)')
     ABI_CHECK(size(ylm_gr,3)==mpsang**2,'BUG: wrong ylm_gr size (3)')
   end if
 end if

 ! Get (k+G) and |k+G|
 ABI_MALLOC(kpgnorm,(npw))
 ABI_MALLOC(kpgnorm_inv,(npw))

 ig0=-1 ! index of |k+g|=0 vector

 if (useylm==1) then
   ABI_MALLOC(kpgc,(npw,3))
   if (ider>=1) then
     ABI_MALLOC(kpgn,(npw,3))
   end if
   if (nkpg<3) then
!$OMP PARALLEL DO PRIVATE(ig,kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3)
     do ig=1,npw
       kpg1=kpt(1)+dble(kg(1,ig))
       kpg2=kpt(2)+dble(kg(2,ig))
       kpg3=kpt(3)+dble(kg(3,ig))
       kpgc1=kpg1*gprimd(1,1)+kpg2*gprimd(1,2)+kpg3*gprimd(1,3)
       kpgc2=kpg1*gprimd(2,1)+kpg2*gprimd(2,2)+kpg3*gprimd(2,3)
       kpgc3=kpg1*gprimd(3,1)+kpg2*gprimd(3,2)+kpg3*gprimd(3,3)
       kpgc(ig,1)=kpgc1
       kpgc(ig,2)=kpgc2
       kpgc(ig,3)=kpgc3
       kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
       if (kpgnorm(ig)<=tol_norm) ig0=ig
       if (ider>=1) then
         kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol_norm)
         kpgn(ig,1)=kpg1*kpgnorm_inv(ig)
         kpgn(ig,2)=kpg2*kpgnorm_inv(ig)
         kpgn(ig,3)=kpg3*kpgnorm_inv(ig)
       end if
     end do
   else
!$OMP PARALLEL DO PRIVATE(ig,kpgc1,kpgc2,kpgc3)
     do ig=1,npw
       kpgc1=kpg(ig,1)*gprimd(1,1)+kpg(ig,2)*gprimd(1,2)+kpg(ig,3)*gprimd(1,3)
       kpgc2=kpg(ig,1)*gprimd(2,1)+kpg(ig,2)*gprimd(2,2)+kpg(ig,3)*gprimd(2,3)
       kpgc3=kpg(ig,1)*gprimd(3,1)+kpg(ig,2)*gprimd(3,2)+kpg(ig,3)*gprimd(3,3)
       kpgc(ig,1)=kpgc1
       kpgc(ig,2)=kpgc2
       kpgc(ig,3)=kpgc3
       kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
       if (kpgnorm(ig)<=tol_norm) ig0=ig
       if (ider>=1) then
         kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol_norm)
         kpgn(ig,1:3)=kpg(ig,1:3)*kpgnorm_inv(ig)
       end if
     end do
   end if
 else
   if (nkpg<3) then
     ecut=huge(0.0d0)*0.1d0;ecutsm=zero;effmass_free=one
     ! Note that with ecutsm=0, the right kinetic energy is computed
     call mkkin(ecut,ecutsm,effmass_free,gmet,kg,kpgnorm,kpt,npw,0,0)
!$OMP PARALLEL DO
     do ig=1,npw
       kpgnorm(ig)=sqrt(renorm_factor*kpgnorm(ig))
       kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol_norm)
       if (kpgnorm(ig)<=tol_norm) ig0=ig
     end do
   else
!$OMP PARALLEL DO PRIVATE(ig,kpgc1,kpgc2,kpgc3)
     do ig=1,npw
       kpgc1=kpg(ig,1)*gprimd(1,1)+kpg(ig,2)*gprimd(1,2)+kpg(ig,3)*gprimd(1,3)
       kpgc2=kpg(ig,1)*gprimd(2,1)+kpg(ig,2)*gprimd(2,2)+kpg(ig,3)*gprimd(2,3)
       kpgc3=kpg(ig,1)*gprimd(3,1)+kpg(ig,2)*gprimd(3,2)+kpg(ig,3)*gprimd(3,3)
       kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
       kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol_norm)
       if (kpgnorm(ig)<=tol_norm) ig0=ig
     end do
   end if
 end if

 ! Need rprimd in some cases
 if (ider>=1.and.useylm==1.and.ig0>0) then
   do mu=1,3
     do nu=1,3
       rprimd(mu,nu)=gprimd(mu,1)*rmet(1,nu)+gprimd(mu,2)*rmet(2,nu)+gprimd(mu,3)*rmet(3,nu)
     end do
   end do
 end if

 ! Allocate several temporary arrays
 ABI_MALLOC(wk_ffnl1,(npw))
 ABI_MALLOC(wk_ffnl2,(npw))
 ABI_MALLOC(wk_ffnl3,(npw))
 ABI_MALLOC(wk_ffspl,(mqgrid,2))

 if (ider>=1.and.useylm==1) then
   ABI_MALLOC(dffnl_red,(npw,3))
   if (idir/=0) then
     ABI_MALLOC(dffnl_cart,(npw,3))
   end if
   if (idir>0) then
     ABI_MALLOC(dffnl_tmp,(npw))
   end if
 end if
 if (ider>=2 .and. useylm==1) then
   ABI_MALLOC(d2ffnl_red,(npw,6))
   if (idir==4) then
     ABI_MALLOC(d2ffnl_cart,(npw,6))
     ABI_MALLOC(d2ffnl_tmp,(npw))
   end if
 end if

 ! Loop over types of atoms
 ffnl = zero; cnt = 0
 do itypat=1,ntypat

   ! Loop over (l,m,n) values
   iln0=0; nlmn=count(indlmn(3,:,itypat)>0)

   do ilmn=1,nlmn
     il=indlmn(1,ilmn,itypat)
     im=indlmn(2,ilmn,itypat)
     ilm =indlmn(4,ilmn,itypat)
     iln =indlmn(5,ilmn,itypat)
     iffnl=ilmn;if (useylm==0) iffnl=iln

     ! Special case: we enter the loop in case of spin-orbit calculation
     ! even if the psp has no spin-orbit component.
     if (indlmn(6,ilmn,itypat) ==1 .or. pspso(itypat) /=0) then

       ! Compute FFNL only if ekb>0 or paw
       if (usepaw==1) testnl=.true.
       if (usepaw==0) testnl=(abs(ekb(iln,itypat))>tol_norm)

       if (testnl) then
         cnt = cnt + 1
         if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism (optional)
         !
         ! Store form factors (from ffspl)
         ! -------------------------------
         ! MG: This part is an hotspot of in the EPH code due to the large number of k-points used
         ! To improve memory locality, I tried to:
         !      1) call a new version of splfit that operates on wk_ffspl with shape: (2,mqgrid)
         !      2) pass a sorted kpgnorm array and then rearrange the output spline
         ! but I didn't manage to make it significantly faster.
         ! For the time being, we rely on MPI-parallelsism via the optional MPI communicator.
         if (iln > iln0) then
           wk_ffspl(:,:)=ffspl(:,:,iln,itypat)
           ider_tmp = min(ider, 1)
           call splfit(qgrid,wk_ffnl2,wk_ffspl,ider_tmp,kpgnorm,wk_ffnl1,mqgrid,npw)
           if (ider == 2) then
             call splfit(qgrid,wk_ffnl3,wk_ffspl,ider,kpgnorm,wk_ffnl1,mqgrid,npw)
           end if
         end if

         ! Store FFNL and FFNL derivatives
         ! -------------------------------

         ! =========================================================================
         ! A-USE OF SPHER. HARMONICS IN APPLICATION OF NL OPERATOR:
         ! ffnl(K,l,m,n)=fnl(K).Ylm(K)
         ! --if (idir==0)
         ! ffnl_prime(K,1:3,l,m,n)=3 reduced coordinates of d(ffnl)/dK^cart
         ! =fnl_prime(K).Ylm(K).K^red_i/|K|+fnl(K).(dYlm/dK^cart)^red_i
         ! --if (0<idir<4)
         ! ffnl_prime(K,l,m,n)=cart. coordinate idir of d(ffnl)/dK^red
         ! --if (idir==4)
         ! ffnl_prime(K,l,m,n)=3 cart. coordinates of d(ffnl)/dK^red
         ! --if (-7<=idir<0) - |idir|=(mu,nu) (1->11,2->22,3->33,4->32,5->31,6->21)
         ! ffnl_prime(K,l,m,n)=1/2 [d(ffnl)/dK^cart_mu K^cart_nu + d(ffnl)/dK^cart_nu K^cart_mu]
         ! ffnl_prime_prime(K,l,m,n)=6 reduced coordinates of d2(ffnl)/dK^cart.dK^cart

         if (useylm==1) then
!$OMP PARALLEL DO
           do ig=1,npw
             ffnl(ig,1,iffnl,itypat)=ylm(ig,ilm)*wk_ffnl1(ig)
           end do

           if (ider>=1) then
!$OMP PARALLEL DO COLLAPSE(2)
             do mu=1,3
               do ig=1,npw
                 dffnl_red(ig,mu)=ylm(ig,ilm)*wk_ffnl2(ig)*kpgn(ig,mu)+ylm_gr(ig,mu,ilm)*wk_ffnl1(ig)
               end do
             end do
             ! Special cases |k+g|=0
             if (ig0>0) then
               do mu=1,3
                 dffnl_red(ig0,mu)=zero
                 if (il==1) then
                   !Retrieve 1st-deriv. of ffnl at q=zero according to spline routine
                   yp1=(wk_ffspl(2,1)-wk_ffspl(1,1))/qgrid(2)-sixth*qgrid(2)*(two*wk_ffspl(1,2)+wk_ffspl(2,2))
                   fact=yp1*sqrt(three/four_pi)
                   if (im==-1) dffnl_red(ig0,mu)=fact*rprimd(2,mu)
                   if (im== 0) dffnl_red(ig0,mu)=fact*rprimd(3,mu)
                   if (im==+1) dffnl_red(ig0,mu)=fact*rprimd(1,mu)
                 end if
               end do
             end if
             if (idir==0) then
!$OMP PARALLEL DO COLLAPSE(2)
               do mu=1,3
                 do ig=1,npw
                   ffnl(ig,1+mu,iffnl,itypat)=dffnl_red(ig,mu)
                 end do
               end do
             else
               dffnl_cart=zero
!$OMP PARALLEL DO COLLAPSE(2)
               do mu=1,3
                 do ig=1,npw
                   do nu=1,3
                     dffnl_cart(ig,mu)=dffnl_cart(ig,mu)+dffnl_red(ig,nu)*gprimd(mu,nu)
                   end do
                 end do
               end do
               if (idir>=1.and.idir<=3) then
                 dffnl_tmp=zero
!$OMP PARALLEL PRIVATE(nu,ig)
!$OMP DO
                 do ig=1,npw
                   do nu=1,3
                     dffnl_tmp(ig)=dffnl_tmp(ig) + dffnl_cart(ig,nu)*gprimd(nu,idir)
                   end do
                 end do
!$OMP END DO
!$OMP WORKSHARE
                 ffnl(:,2,iffnl,itypat)=dffnl_tmp(:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
               else if (idir==4) then
                 do mu=1,3
!$OMP PARALLEL PRIVATE(nu,ig)
!$OMP WORKSHARE
                   dffnl_tmp=zero
!$OMP END WORKSHARE
!$OMP DO
                   do ig=1,npw
                     do nu=1,3
                       dffnl_tmp(ig)=dffnl_tmp(ig) + dffnl_cart(ig,nu)*gprimd(nu,mu)
                     end do
                   end do
!$OMP END DO
!$OMP WORKSHARE
                   ffnl(:,1+mu,iffnl,itypat)=dffnl_tmp(:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
                 end do
               else if (idir/=-7) then
                 mu=abs(idir);mua=alpha(mu);mub=beta(mu)
!$OMP PARALLEL DO
                 do ig=1,npw
                   ffnl(ig,2,iffnl,itypat)=0.5d0* (dffnl_cart(ig,mua)*kpgc(ig,mub) + dffnl_cart(ig,mub)*kpgc(ig,mua))
                 end do
               else if (idir==-7) then
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(mua, mub)
                 do mu=1,6
                   do ig=1,npw
                     mua=alpha(mu);mub=beta(mu)
                     ffnl(ig,1+mu,iffnl,itypat)=0.5d0 * (dffnl_cart(ig,mua)*kpgc(ig,mub) + dffnl_cart(ig,mub)*kpgc(ig,mua))
                   end do
                 end do
               end if
             end if
           end if

           if (ider==2) then
             do mu=1,6
               mua=alpha(mu);mub=beta(mu)
               rmetab=rmet(mua,mub)
!$OMP PARALLEL DO
               do ig=1,npw
                 d2ffnl_red(ig,mu)= &
                 ylm_gr(ig,3+mu,ilm)*wk_ffnl1(ig) &
                 + (rmetab-kpgn(ig,mua)*kpgn(ig,mub))*ylm(ig,ilm)*wk_ffnl2(ig)*kpgnorm_inv(ig) &
                 + ylm(ig,ilm)*kpgn(ig,mua)*kpgn(ig,mub)*wk_ffnl3(ig) &
                 + (ylm_gr(ig,mua,ilm)*kpgn(ig,mub)+ylm_gr(ig,mub,ilm)*kpgn(ig,mua))*wk_ffnl2(ig)
               end do
               ! Special cases |k+g|=0
               if (ig0>0) then
                 d2ffnl_red(ig0,mu)=zero
                 if (il==0) then
                   d2ffnl_red(ig0,mu)=wk_ffspl(1,2)*rmetab/sqrt(four_pi)
                 end if
                 if (il==2) then
                   fact=wk_ffspl(1,2)*quarter*sqrt(15._dp/pi)
                   if (im==-2) d2ffnl_red(ig0,mu)=fact*(rprimd(1,mua)*rprimd(2,mub)+rprimd(2,mua)*rprimd(1,mub))
                   if (im==-1) d2ffnl_red(ig0,mu)=fact*(rprimd(2,mua)*rprimd(3,mub)+rprimd(3,mua)*rprimd(2,mub))
                   if (im==+1) d2ffnl_red(ig0,mu)=fact*(rprimd(1,mua)*rprimd(3,mub)+rprimd(3,mua)*rprimd(1,mub))
                   if (im==+2) d2ffnl_red(ig0,mu)=fact*(rprimd(1,mua)*rprimd(1,mub)-rprimd(2,mua)*rprimd(2,mub))
                   if (im== 0) d2ffnl_red(ig0,mu)=(fact/sqrt3)*(two*rprimd(3,mua)*rprimd(3,mub) &
                                                  -rprimd(1,mua)*rprimd(1,mub)-rprimd(2,mua)*rprimd(2,mub))
                 end if
               end if
             end do
             if (idir==0) then
!$OMP PARALLEL DO COLLAPSE(2)
               do mu=1,6
                 do ig=1,npw
                   ffnl(ig,4+mu,iffnl,itypat)=d2ffnl_red(ig,mu)
                 end do
               end do
             else if (idir==4) then
               d2ffnl_cart=zero
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(mu,mua,mub,ig,nu,nua,nub)
               do mu=1,6
                 do ig=1,npw
                   mua=alpha(mu);mub=beta(mu)
                   do nua=1,3
                     do nub=1,3
                       nu=gamma(nua,nub)
                       d2ffnl_cart(ig,mu)=d2ffnl_cart(ig,mu)+d2ffnl_red(ig,nu)*gprimd(mua,nua)*gprimd(mub,nub)
                     end do
                   end do
                 end do
               end do
               do mu=1,6
                 mua=alpha(mu);mub=beta(mu)
!$OMP PARALLEL PRIVATE(nu,nua,nub,ig)
!$OMP WORKSHARE
                 d2ffnl_tmp=zero
!$OMP END WORKSHARE
!$OMP DO
                 do ig=1,npw
                   do nua=1,3
                     do nub=1,3
                       nu=gamma(nua,nub)
                       d2ffnl_tmp(ig)=d2ffnl_tmp(ig)+d2ffnl_cart(ig,nu)*gprimd(nua,mua)*gprimd(nub,mub)
                     end do
                   end do
                 end do
!$OMP END DO
!$OMP WORKSHARE
                 ffnl(:,4+mu,iffnl,itypat)=d2ffnl_tmp(:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
               end do
             end if
           end if

           ! =========================================================================
           ! B-USE OF LEGENDRE POLYNOMIAL IN APPLICATION OF NL OPERATOR:
           ! ffnl(K,l,n)=fnl(K)/|K|^l
           ! ffnl_prime(K,l,n)=(fnl_prime(K)-l*fnl(K)/|K|)/|K|^(l+1)
           ! ffnl_prime_prime(K,l,n)=(fnl_prime_prime(K)-(2*l+1)*fnl_prime(K)/|K|
           ! +l*(l+2)*fnl(K)/|K|^2)/|K|^(l+2)
         else if (iln>iln0) then

           if (il==0) then
!$OMP PARALLEL DO
             do ig=1,npw
               ffnl(ig,1,iffnl,itypat)=wk_ffnl1(ig)
             end do
           else
!$OMP PARALLEL DO
             do ig=1,npw
               ffnl(ig,1,iffnl,itypat)=wk_ffnl1(ig)*kpgnorm_inv(ig)**il
             end do
           end if
           if (ider>=1) then
!$OMP PARALLEL DO
             do ig=1,npw
               ffnl(ig,2,iffnl,itypat)= (wk_ffnl2(ig)-dble(il)*wk_ffnl1(ig)*kpgnorm_inv(ig))*kpgnorm_inv(ig)**(il+1)
             end do
             if (ider==2) then
!$OMP PARALLEL DO
               do ig=1,npw
                 ffnl(ig,3,iffnl,itypat)= (wk_ffnl3(ig)-       &
                   dble(2*il+1)*wk_ffnl2(ig)*kpgnorm_inv(ig)+   &
                   dble(il*(il+2))*wk_ffnl1(ig)*kpgnorm_inv(ig)**2)*kpgnorm_inv(ig)**(il+2)
               end do
             end if
           end if

         end if  ! Use of Ylm or not

       else
         ! No NL part
!$OMP PARALLEL DO COLLAPSE(2)
         do mu=1,dimffnl
           do ig=1,npw
             ffnl(ig,mu,iffnl,itypat)=zero
           end do
         end do

       end if ! testnl (a nonlocal part exists)
     end if ! special case: spin orbit calc. & no spin-orbit psp

     if (iln > iln0) iln0 = iln

   end do ! loop over (l,m,n) values
 end do ! loop over atom types

 ABI_FREE(kpgnorm_inv)
 ABI_FREE(kpgnorm)
 ABI_FREE(wk_ffnl1)
 ABI_FREE(wk_ffnl2)
 ABI_FREE(wk_ffnl3)
 ABI_FREE(wk_ffspl)

 ! Optional deallocations.
 ABI_SFREE(kpgc)
 ABI_SFREE(kpgn)
 ABI_SFREE(dffnl_red)
 ABI_SFREE(d2ffnl_red)
 ABI_SFREE(dffnl_cart)
 ABI_SFREE(d2ffnl_cart)
 ABI_SFREE(dffnl_tmp)
 ABI_SFREE(d2ffnl_tmp)

 if (nprocs > 1) then
   ! Blocking/non-blocking depending on the presence of request.
   if (present(request)) then
     call xmpi_isum_ip(ffnl, comm, request, ierr)
   else
     call xmpi_sum(ffnl, comm, ierr)
   end if
 else
   if (present(request)) request = xmpi_request_null
 end if

 call timab(16, 2, tsec)

end subroutine mkffnl
!!***

end module m_mkffnl
!!***

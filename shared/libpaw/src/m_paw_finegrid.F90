!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_finegrid
!! NAME
!!  m_paw_finegrid
!!
!! FUNCTION
!!  This module contains a set of routines to compute various quantities
!!  on the fine grid around a given atom.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_finegrid

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

 use m_pawtab,      only : pawtab_type
 use m_paw_sphharm, only : initylmr
 use m_paw_numeric, only : paw_jbessel,paw_splint,paw_uniform_splfit,paw_sort_dp

 implicit none

 private

!public procedures.
 public :: pawgylm      ! g_l(r-R)*Y_lm(r-R) (and derivatives)
 public :: pawgylmg     ! Fourier transform of g_l(r-R)*Y_lm(r-R), plane-waves case
 public :: pawrfgd_fft  ! r-R, plane-waves case
 public :: pawrfgd_wvl  ! r-R, wavelets case
 public :: pawexpiqr    ! exp(i.q.(r-R))

!declarations for the whole module (were needed to replace the statement functions)
!MG: Why this? Global variables are powerful but extremely DANGEROUS
integer,private,save :: lambda
real(dp),private,save :: pi_over_rshp,sigma
real(dp),private,save,allocatable :: alpha(:,:),qq(:,:)
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_finegrid/pawgylm
!! NAME
!! pawgylm
!!
!! FUNCTION
!! Compute g_l(r-R)*Y_lm(r-R) (and derivatives) on the fine (rectangular) grid
!! around one atom (g_l=radial shape function).
!! R is the position of the atom
!!
!! INPUTS
!!  lm_size=number of lm components to be calculated
!!  nfgd= number of (fine grid) FFT points in the paw sphere around current atom
!!  optgr0= 1 if g_l(r-R)*Y_lm(r-R) are computed
!!  optgr1= 1 if first derivatives of g_l(r-R)*Y_lm(r-R) are computed
!!  optgr2= 1 if second derivatives of g_l(r-R)*Y_lm(r-R) are computed
!!  pawtab <type(pawtab_type)>=paw tabulated starting data for current atom
!!  rfgd(3,nfgd)= coordinates of r-R on the fine rect. grid around current atom
!!
!! OUTPUT
!!  if (optgr0==1)
!!    gylm(nfgd,lm_size)= g_l(r-R)*Y_lm(r-R) around current atom
!!  if (optgr1==1)
!!    gylmgr(3,nfgd,lm_size)= derivatives of g_l(r-R)*Y_lm(r-R) wrt cart. coordinates
!!  if (optgr2==1)
!!    gylmgr2(6,nfgd,lm_size)= second derivatives of g_l(r-R)*Y_lm(r-R) wrt cart. coordinates
!!
!! PARENTS
!!      m_paw_denpot,m_paw_dfpt,m_paw_nhat,m_pawdij
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawgylm(gylm,gylmgr,gylmgr2,lm_size,nfgd,optgr0,optgr1,optgr2,pawtab,rfgd)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lm_size,nfgd,optgr0,optgr1,optgr2
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(in) :: rfgd(:,:)
 real(dp),intent(out) :: gylm(nfgd,optgr0*lm_size)
 real(dp),intent(out) :: gylmgr(3,nfgd,optgr1*lm_size)
 real(dp),intent(out) :: gylmgr2(6,nfgd,optgr2*lm_size)

!Local variables ------------------------------
!scalars
 integer :: ic,ilm,izero,l_size,ll,normchoice,option,shape_type
 real(dp) :: arg
 real(dp) :: jbes1,jbes2,jbesp1,jbesp2,jbespp1,jbespp2,rcut
 real(dp) :: splfact
 logical :: compute_gr0,compute_gr1,compute_gr2
 character(len=500) :: msg
!arrays
 integer,allocatable :: isort(:)
 real(dp),parameter :: ffact(1:9)=(/1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,&
&                                   135135._dp,2027025._dp,34459425._dp/)
 real(dp),parameter :: toldev=tol3
 real(dp) :: ss(3)
 real(dp),allocatable :: cc(:,:),d2gfact(:,:),d2shpfuncnum(:,:),dgfact(:,:)
 real(dp),allocatable :: dshpfuncnum(:,:),gfact(:,:)
 real(dp),allocatable :: rnrm(:),rnrm_inv(:),rnrm_sort(:)
 real(dp),allocatable :: shpfuncnum(:,:),work(:),ylmr(:,:),ylmrgr(:,:,:)

! *************************************************************************

 if (optgr0==0.and.optgr1==0.and.optgr2==0) return
 if (nfgd==0) return

!Compatibility test
!==========================================================
 if (size(rfgd)/=3*nfgd) then
   msg='rfgd array must be allocated at rfgd(3,nfgd)!'
   MSG_BUG(msg)
 end if
 if (pawtab%lcut_size>9) then
   msg='l_size>10 forbidden!'
   MSG_BUG(msg)
 end if
 if (pawtab%shape_type==1.and.pawtab%shape_lambda<2) then
   msg='Exponent lambda of gaussian shape function must be > 1!'
   MSG_ERROR(msg)
 end if

!Initializations
!==========================================================
!Options for computation
 compute_gr0=(optgr0==1.or.optgr1==1.or.optgr2==1)
 compute_gr1=(optgr1==1.or.optgr2==1)
 compute_gr2=(optgr2==1)
 l_size=pawtab%lcut_size

!Norms of vectors around the atom
 LIBPAW_ALLOCATE(rnrm,(nfgd))
 izero=-1
 do ic=1,nfgd
   rnrm(ic)=sqrt(rfgd(1,ic)**2+rfgd(2,ic)**2+rfgd(3,ic)**2)
   if (rnrm(ic)<=tol10) izero=ic  ! Has to be consistent with initylmr !!
 end do

!Initializations
 if (optgr0==1) gylm=zero
 if (optgr1==1) gylmgr=zero
 if (optgr2==1) gylmgr2=zero

!Some definitions concerning shape function g_l(r)
 shape_type=pawtab%shape_type
 sigma=pawtab%shape_sigma;lambda=pawtab%shape_lambda
 pi_over_rshp=pi/pawtab%rshp
 rcut=tol12+pawtab%rshp
 if (shape_type==3) then
   LIBPAW_ALLOCATE(alpha,(2,l_size))
   LIBPAW_ALLOCATE(qq,(2,l_size))
   do ll=1,l_size
     alpha(1:2,ll)=pawtab%shape_alpha(1:2,ll)
     qq(1:2,ll)=pawtab%shape_q(1:2,ll)
   end do
 end if

!If needed, sort selected radii by increasing norm
 if (shape_type==-1) then
   LIBPAW_ALLOCATE(isort,(nfgd))
   LIBPAW_ALLOCATE(rnrm_sort,(nfgd))
   do ic=1,nfgd
     isort(ic)=ic
   end do
   rnrm_sort(1:nfgd)=rnrm(1:nfgd)
   call paw_sort_dp(nfgd,rnrm_sort,isort,tol16)
 end if

!If shape function is "numeric", spline it onto selected radii
 if (shape_type==-1) then
   LIBPAW_ALLOCATE(work,(nfgd))
   if (compute_gr0) then
     LIBPAW_ALLOCATE(shpfuncnum,(nfgd,l_size))
     do ll=1,l_size
       call paw_splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%shapefunc(:,ll),&
&       pawtab%dshpfunc(:,ll,2),nfgd,rnrm_sort,work)
       do ic=1,nfgd
         shpfuncnum(isort(ic),ll)=work(ic)
       end do
     end do
   end if
   if(compute_gr1) then
     LIBPAW_ALLOCATE(dshpfuncnum,(nfgd,l_size))
     do ll=1,l_size
       call paw_splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%dshpfunc(:,ll,1),&
&       pawtab%dshpfunc(:,ll,3),nfgd,rnrm_sort,work)
       do ic=1,nfgd
         dshpfuncnum(isort(ic),ll)=work(ic)
       end do
     end do
   end if
   if(compute_gr2) then
     LIBPAW_ALLOCATE(d2shpfuncnum,(nfgd,l_size))
     do ll=1,l_size
       call paw_splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%dshpfunc(:,ll,2),&
&       pawtab%dshpfunc(:,ll,4),nfgd,rnrm_sort,work)
       do ic=1,nfgd
         d2shpfuncnum(isort(ic),ll)=work(ic)
       end do
     end do
   end if
   LIBPAW_DEALLOCATE(work)
 end if

 if (shape_type==-1)  then
   LIBPAW_DEALLOCATE(isort)
   LIBPAW_DEALLOCATE(rnrm_sort)
 end if

!If needed, compute limits at r=0 of shape function and derivatives
 if (izero>0) then
   LIBPAW_ALLOCATE(cc,(3,min(l_size,3)))
   cc=zero
   if (shape_type==-1) then
     splfact=(pawtab%rad_for_spline(4)-pawtab%rad_for_spline(1))&
&     /(pawtab%rad_for_spline(3)-pawtab%rad_for_spline(2))
   end if
   do ll=1,min(l_size,3)
!    cc(2,l) is g_prime(0)
     if (optgr0==1.or.optgr1==1.or.optgr2==1) then
       if (shape_type==-1) then
         ss(1:3)=pawtab%shapefunc(2:4,ll)/pawtab%rad_for_spline(2:4)**(ll-1)
         cc(1,ll)=ss(3)+(ss(1)-ss(2))*splfact
       else if (shape_type==1.or.shape_type==2) then
         cc(1,ll)=one
       else if (shape_type==3) then
         cc(1,ll)=(alpha(1,ll)*qq(1,ll)**(ll-1) &
&         +alpha(2,ll)*qq(2,ll)**(ll-1))/ffact(ll)
       end if
       cc(1,ll)=cc(1,ll)*pawtab%gnorm(ll)
     end if
!    cc(2,l) is g_prime(0)
     if (optgr1==1.or.optgr2==1) then
       if (shape_type==-1) then
         ss(1:3)=(ss(1:3)-cc(1,ll))/pawtab%rad_for_spline(2:4)
         cc(2,ll)=ss(3)+(ss(1)-ss(2))*splfact
       else if (shape_type==1.and.lambda==1) then
         cc(2,ll)=-one/sigma
       else
         cc(2,ll)=zero
       end if
       cc(2,ll)=cc(2,ll)*pawtab%gnorm(ll)
     end if
!    cc(3,l) is g_prime_prime(0)
     if (optgr2==1) then
       if (shape_type==-1) then
         ss(1:3)=(ss(1:3)-cc(2,ll))/pawtab%rad_for_spline(2:4)
         cc(3,ll)=two*(ss(3)+(ss(1)-ss(2))*splfact)
       else if (shape_type==1) then
         if (lambda==1) cc(3,ll)=one/sigma**2
         if (lambda==2) cc(3,ll)=-two/sigma**2
         if (lambda >2) cc(3,ll)=zero
       else if (shape_type==2) then
         cc(3,ll)=-(two/three)*pi_over_rshp**2
       else if (shape_type==3) then
         cc(3,ll)=-(alpha(1,ll)*qq(1,ll)**(ll+1) &
&         +alpha(2,ll)*qq(2,ll)**(ll+1))/ffact(ll+1)
       end if
       cc(3,ll)=cc(3,ll)*pawtab%gnorm(ll)
     end if
   end do
 end if

!Y_lm(r-R) calculation
!==========================================================
 normchoice=1 ; option=max(optgr0,2*optgr1,3*optgr2)
 if(compute_gr0)  then
   LIBPAW_ALLOCATE(ylmr,(l_size**2,nfgd))
 end if
 if(compute_gr1.and.(.not.compute_gr2))  then
   LIBPAW_ALLOCATE(ylmrgr,(3,l_size**2,nfgd))
 end if
 if(compute_gr2)  then
   LIBPAW_ALLOCATE(ylmrgr,(9,l_size**2,nfgd))
 end if
 if (compute_gr0.and.(.not.compute_gr1).and.(.not.compute_gr2)) then
   call initylmr(l_size,normchoice,nfgd,rnrm,option,rfgd,ylmr)
 else
   call initylmr(l_size,normchoice,nfgd,rnrm,option,rfgd,ylmr,ylmrgr)
 end if

!gl(r) and derivatives calculation for l>=0
!==========================================================
!Compute gl(r), gl_prime(r)/r and (gl_prime_prime(r)-gl_prime(r)/r)/r**2
 if (compute_gr0)  then
   LIBPAW_BOUND2_ALLOCATE(gfact,BOUNDS(1,nfgd),BOUNDS(0,l_size-1))
   gfact(:,:)=zero
 end if
 if (compute_gr1)  then
   LIBPAW_BOUND2_ALLOCATE(dgfact,BOUNDS(1,nfgd),BOUNDS(0,l_size-1))
   dgfact(:,:)=zero
 end if
 if (compute_gr2)  then
   LIBPAW_BOUND2_ALLOCATE(d2gfact,BOUNDS(1,nfgd),BOUNDS(0,l_size-1))
   d2gfact(:,:)=zero
 end if
 if(compute_gr1) then
   LIBPAW_ALLOCATE(rnrm_inv,(nfgd))
   do ic=1,nfgd
     if (ic/=izero) rnrm_inv(ic)=one/rnrm(ic)
   end do
   if (izero>0) rnrm_inv(izero)=zero
 end if

!----- type -1 -----
 if (shape_type==-1) then
   if (compute_gr0) then
     do ll=0,l_size-1
       do ic=1,nfgd
         if (rnrm(ic)<=rcut) then
           gfact(ic,ll)=shpfuncnum(ic,ll+1)
         end if
       end do
     end do
   end if
   if (compute_gr1) then
     do ll=0,l_size-1
       do ic=1,nfgd
         if (rnrm(ic)<=rcut) then
           dgfact(ic,ll)=dshpfuncnum(ic,ll+1)*rnrm_inv(ic)
         end if
       end do
     end do
   end if
   if(compute_gr2) then
     do ll=0,l_size-1
       do ic=1,nfgd
         if (rnrm(ic)<=rcut) then
           d2gfact(ic,ll)=(d2shpfuncnum(ic,ll+1)-dgfact(ic,ll))*rnrm_inv(ic)**2
         end if
       end do
    end do
   end if

!  ----- type 1 or 2 -----
 else if (shape_type==1.or.shape_type==2) then
!  FIRST COMPUTE FACTORS FOR l=0
   if (optgr0==1.and.optgr1==0.and.optgr2==0) then
     if (shape_type==1) then
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc1_0(arg)
         else if (arg<=rcut) then
           gfact(ic,0)=shapefunc1(arg)
         end if
       end do
     else ! shape_type==2
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc2_0(arg)
         else if (arg<=rcut) then
           gfact(ic,0)=shapefunc2(arg)
         end if
       end do
     end if
   else if (optgr1==1.and.optgr2==0) then
     if (shape_type==1) then
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc1_0(arg)
           if (lambda==2) then
             dgfact(ic,0)=dshpfunc1_ovr_0_2(arg)
           else ! lambda>2
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
           end if
         else if (arg<=rcut) then
           gfact(ic,0)=shapefunc1(arg)
           dgfact(ic,0)=dshpfunc1(arg)*rnrm_inv(ic)
         end if
       end do
     else ! shape_type==2
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc2_0(arg)
           dgfact(ic,0)=dshpfunc2_ovr_0(arg)
         else if (arg<=rcut) then
           gfact(ic,0)=shapefunc2(arg)
           dgfact(ic,0)=dshpfunc2(arg)*rnrm_inv(ic)
         end if
       end do
     end if
   else if (optgr2==1) then
     if (shape_type==1) then
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc1_0(arg)
           if (lambda==2) then
             dgfact(ic,0)=dshpfunc1_ovr_0_2(arg)
             d2gfact(ic,0)=d2shpfunc1_ovr2_0_2(arg)
           else if (lambda==3) then
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
             if (ic/=izero) then
               d2gfact(ic,0)=d2shpfunc1_ovr2_0_3(arg)
             else
               d2gfact(ic,0)=zero ! Diverging case
             end if
           else if (lambda==4) then
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
             d2gfact(ic,0)=d2shpfunc1_ovr2_0_4(arg)
           else ! lambda>4
             dgfact(ic,0)=dshpfunc1_ovr_0(arg)
             d2gfact(ic,0)=d2shpfunc1_ovr2_0(arg)
           end if
         else if (arg<=rcut) then
           gfact(ic,0)=shapefunc1(arg)
           dgfact(ic,0)=dshpfunc1(arg)*rnrm_inv(ic)
           d2gfact(ic,0)=(d2shpfunc1(arg)-dgfact(ic,0))*rnrm_inv(ic)**2
         end if
       end do
     else ! shape_type==2
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<toldev) then
           gfact(ic,0)=shapefunc2_0(arg)
           dgfact(ic,0)=dshpfunc2_ovr_0(arg)
           d2gfact(ic,0)=d2shpfunc2_ovr2_0(arg)
         else if (arg<=rcut) then
           gfact(ic,0)=shapefunc2(arg)
           dgfact(ic,0)=dshpfunc2(arg)*rnrm_inv(ic)
           d2gfact(ic,0)=(d2shpfunc2(arg)-dgfact(ic,0))*rnrm_inv(ic)**2
         end if
       end do
     end if
   end if

!  THEN COMPUTE FACTORS FOR l>0 (from l=0)
   if (compute_gr0) then
     if (l_size>1) then
       do ll=1,l_size-1
         do ic=1,nfgd
           gfact(ic,ll)=pawtab%gnorm(ll+1)*gfact(ic,0)*rnrm(ic)**ll
         end do
       end do
     end if
   end if
   if (compute_gr1) then
     if (l_size>1) then
       do ic=1,nfgd
         dgfact(ic,1)=pawtab%gnorm(2)*(gfact(ic,0)*rnrm_inv(ic)+dgfact(ic,0)*rnrm(ic))
       end do
     end if
     if (l_size>2) then
       do ic=1,nfgd
         dgfact(ic,2)=pawtab%gnorm(3)*(two*gfact(ic,0)+dgfact(ic,0)*rnrm(ic)**2)
       end do
     end if
     if (l_size>3) then
       do ll=3,l_size-1
         do ic=1,nfgd
           dgfact(ic,ll)=pawtab%gnorm(ll+1) &
&           *(dble(ll)*gfact(ic,0)*rnrm(ic)**(ll-2)+dgfact(ic,0)*rnrm(ic)**ll)
         end do
       end do
     end if
   end if
   if (compute_gr2) then
     if (l_size>1) then
       do ic=1,nfgd
         d2gfact(ic,1)=pawtab%gnorm(2) &
&         *(-gfact(ic,0)*rnrm_inv(ic)**3+two*dgfact(ic,0)*rnrm_inv(ic)+d2gfact(ic,0)*rnrm(ic))
       end do
     end if
     if (l_size>2) then
       do ic=1,nfgd
         d2gfact(ic,2)=pawtab%gnorm(3)*(four*dgfact(ic,0)+d2gfact(ic,0)*rnrm(ic)**2)
       end do
     end if
     if (l_size>3) then
       do ic=1,nfgd
         d2gfact(ic,3)=pawtab%gnorm(4) &
&         *(three*gfact(ic,0)*rnrm_inv(ic)+6._dp*dgfact(ic,0)*rnrm(ic)+ d2gfact(ic,0)*rnrm(ic)**3)
       end do
     end if
     if (l_size>4) then
       do ic=1,nfgd
         d2gfact(ic,4)=pawtab%gnorm(5) &
&         *(8._dp*gfact(ic,0)+8._dp*dgfact(ic,0)*rnrm(ic)**2+d2gfact(ic,0)*rnrm(ic)**4)
       end do
     end if
     if (l_size>5) then
       do ll=5,l_size-1
         do ic=1,nfgd
           d2gfact(ic,ll)=pawtab%gnorm(ll+1) &
&           *(dble(ll*(ll-2))*gfact(ic,0)*rnrm(ic)**(ll-4) &
&            +dble(2*ll)*dgfact(ic,0)*rnrm(ic)**(ll-2) &
&            +d2gfact(ic,0)*rnrm(ic)**ll)
         end do
       end do
     end if
   end if
   if (compute_gr0) gfact(:,0)=gfact(:,0)*pawtab%gnorm(1)
   if (compute_gr1) dgfact(:,0)=dgfact(:,0)*pawtab%gnorm(1)
   if (compute_gr2) d2gfact(:,0)=d2gfact(:,0)*pawtab%gnorm(1)

!  ----- type 3 -----
 else if (shape_type==3) then
   if (optgr0==1.and.optgr1==0.and.optgr2==0) then
     do ll=0,l_size-1
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<=rcut) then
           call paw_jbessel(jbes1,jbesp1,jbespp1,ll,0,qq(1,1+ll)*arg)
           call paw_jbessel(jbes2,jbesp2,jbespp2,ll,0,qq(2,1+ll)*arg)
           gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
         end if
       end do
     end do
   else if (optgr1==1.and.optgr2==0) then
     do ll=0,l_size-1
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<=rcut) then
           call paw_jbessel(jbes1,jbesp1,jbespp1,ll,1,qq(1,1+ll)*arg)
           call paw_jbessel(jbes2,jbesp2,jbespp2,ll,1,qq(2,1+ll)*arg)
           gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
           dgfact(ic,ll)=dshpfunc3(jbesp1,jbesp2,ll)*rnrm_inv(ic)
         end if
       end do
     end do
     if (izero>0.and.l_size>=1)  dgfact(izero,0)=-(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))/three
     if (izero>0.and.l_size>=3)  dgfact(izero,2)=two/15._dp*(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))
!    Note: for l=1, dgfact is diverging - d2gfact is diverging for l<4
   else if (optgr2==1) then
     do ll=0,l_size-1
       do ic=1,nfgd
         arg=rnrm(ic)
         if (arg<=rcut) then
           call paw_jbessel(jbes1,jbesp1,jbespp1,ll,2,qq(1,1+ll)*arg)
           call paw_jbessel(jbes2,jbesp2,jbespp2,ll,2,qq(2,1+ll)*arg)
           gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
           dgfact(ic,ll)=dshpfunc3(jbesp1,jbesp2,ll)*rnrm_inv(ic)
           d2gfact(ic,ll)=(d2shpfunc3(jbespp1,jbespp2,ll)-dgfact(ic,ll))*rnrm_inv(ic)**2
         end if
       end do
     end do
     if (izero>0.and.l_size>=1)  dgfact(izero,0)=-(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))/three
     if (izero>0.and.l_size>=3)  dgfact(izero,2)=two/15._dp*(alpha(1,1)*qq(1,1)+alpha(2,1)*qq(2,1))
!    Note: for l=1, dgfact is diverging - d2gfact is diverging for l<4
   end if
 end if

!g_l(r-R)*Y_lm(r-R) calculation
!==========================================================
 if (optgr0==1) then

   do ll=0,l_size-1
     do ilm=ll**2+1,min((ll+1)**2,lm_size)
       do ic=1,nfgd
         gylm(ic,ilm)=gfact(ic,ll)*ylmr(ilm,ic)
       end do
     end do
   end do

!  Special value at r-R=0  (supposing shapefunc(r)->C.r**l when r->0)
   if (izero>0) then
     gylm(izero,1:lm_size)=zero
     if (lm_size>=1) gylm(izero,1)=ylmr(1,izero)*cc(1,1)
   end if

 end if

!d/dr{g_l(r-R)*Y_lm(r-R)} calculation
!==========================================================
 if(optgr1==1) then

   do ll=0,l_size-1
     do ilm=ll**2+1,min((ll+1)**2,lm_size)
       do ic=1,nfgd
         gylmgr(1:3,ic,ilm)=gfact(ic,ll)*ylmrgr(1:3,ilm,ic)&
&         +dgfact(ic,ll)*rfgd(1:3,ic)*ylmr(ilm,ic)
       end do
     end do
   end do

!  Special values at r-R=0  (supposing shapefunc(r)->C.r**l when r->0)
   if (izero>0) then
     gylmgr(1:3,izero,1:lm_size)=zero
     if (lm_size>=1) then
       arg=cc(2,1)/sqrt(four_pi)
       gylmgr(1:3,izero,1)=arg
     end if
     if (lm_size>=2) then
       arg=cc(1,2)*sqrt(three/four_pi)
       gylmgr(2,izero,2)=arg
       if (lm_size>=3) gylmgr(3,izero,3)=arg
       if (lm_size>=4) gylmgr(1,izero,4)=arg
     end if
   end if

 end if

!d2/dridrj{g_l(r-R)*Y_lm(r-R)} calculation
!==========================================================
 if(optgr2==1) then

   do ll=0,l_size-1
     do ilm=ll**2+1,min((ll+1)**2,lm_size)
       do ic=1,nfgd
         gylmgr2(1,ic,ilm)=gfact(ic,ll)*ylmrgr(4,ilm,ic) &
&         +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(1,ic)*ylmrgr(1,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(1,ic)*rfgd(1,ic)
         gylmgr2(2,ic,ilm)=gfact(ic,ll)*ylmrgr(5,ilm,ic) &
&         +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(2,ic)*ylmrgr(2,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(2,ic)*rfgd(2,ic)
         gylmgr2(3,ic,ilm)=gfact(ic,ll)*ylmrgr(6,ilm,ic) &
&         +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(3,ic)*ylmrgr(3,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(3,ic)
         gylmgr2(4,ic,ilm)=gfact(ic,ll)*ylmrgr(7,ilm,ic) &
&         +dgfact(ic,ll)*(rfgd(3,ic)*ylmrgr(2,ilm,ic)+rfgd(2,ic)*ylmrgr(3,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(2,ic)
         gylmgr2(5,ic,ilm)=gfact(ic,ll)*ylmrgr(8,ilm,ic) &
&         +dgfact(ic,ll)*(rfgd(3,ic)*ylmrgr(1,ilm,ic)+rfgd(1,ic)*ylmrgr(3,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(1,ic)
         gylmgr2(6,ic,ilm)=gfact(ic,ll)*ylmrgr(9,ilm,ic) &
&         +dgfact(ic,ll)*(rfgd(1,ic)*ylmrgr(2,ilm,ic)+rfgd(2,ic)*ylmrgr(1,ilm,ic)) &
&         +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(1,ic)*rfgd(2,ic)
       end do
     end do
   end do

!  Special values at r-R=0  (supposing shapefunc(r)->C.r**l when r->0)
   if (izero>0) then
     gylmgr2(1:6,izero,1:lm_size)=zero
     if (lm_size>=1) then
       arg=cc(3,1)/sqrt(four_pi)
       gylmgr2(1:3,izero,1)=arg
     end if
     if (lm_size>=2) then
       arg=cc(2,2)*sqrt(three/four_pi)
       gylmgr2(2,izero,2)=two*arg
       gylmgr2(4,izero,2)=    arg
       if (lm_size>=3) then
         gylmgr2(1,izero,3)=two*arg
         gylmgr2(3,izero,3)=two*arg
       end if
       if (lm_size>=4) then
         gylmgr2(5,izero,4)=arg
         gylmgr2(6,izero,4)=arg
       end if
     end if
     if (lm_size>=5) then
       arg=cc(1,3)*sqrt(15._dp/four_pi)
       gylmgr2(6,izero,5)=arg
       if (lm_size>=6) gylmgr2(4,izero,6)=arg
       if (lm_size>=7) then
         gylmgr2(1,izero,7)=   -arg/sqrt3
         gylmgr2(2,izero,7)=   -arg/sqrt3
         gylmgr2(3,izero,7)=two*arg/sqrt3
       end if
       if (lm_size>=8) gylmgr2(5,izero,8)=arg
       if (lm_size>=9) then
         gylmgr2(1,izero,9)= arg
         gylmgr2(2,izero,9)=-arg
       end if
     end if
   end if

 end if

!Memory deallocation
!==========================================================
 LIBPAW_DEALLOCATE(rnrm)
 if (allocated(cc)) then
   LIBPAW_DEALLOCATE(cc)
 end if
 if (compute_gr0)  then
   LIBPAW_DEALLOCATE(gfact)
 end if
 if (compute_gr1)  then
   LIBPAW_DEALLOCATE(dgfact)
 end if
 if (compute_gr2)  then
   LIBPAW_DEALLOCATE(d2gfact)
 end if
 if (compute_gr1)  then
   LIBPAW_DEALLOCATE(rnrm_inv)
 end if
 if (shape_type==3)  then
   LIBPAW_DEALLOCATE(alpha)
   LIBPAW_DEALLOCATE(qq)
 end if
 if (compute_gr0)  then
   LIBPAW_DEALLOCATE(ylmr)
 end if
 if (compute_gr1)  then
   LIBPAW_DEALLOCATE(ylmrgr)
 end if
 if (shape_type==-1) then
   if (compute_gr0)  then
     LIBPAW_DEALLOCATE(shpfuncnum)
   end if
   if (compute_gr1)  then
     LIBPAW_DEALLOCATE(dshpfuncnum)
   end if
   if (compute_gr2)  then
     LIBPAW_DEALLOCATE(d2shpfuncnum)
   end if
 end if

! -----------------------------------------------------------------
!Small functions related to analytical expression of shape function
 CONTAINS
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/shapefunc1
!  shapefunc1 is g(x) (gaussian)
   function shapefunc1(arg)

     real(dp) :: shapefunc1
     real(dp),intent(in) :: arg
     shapefunc1=exp(-(arg/sigma)**lambda)
   end function shapefunc1
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/shapefunc1_0
!  shapefunc1_0 is g(x) (gaussian) for small x
   function shapefunc1_0(arg)

     real(dp) :: shapefunc1_0
     real(dp),intent(in) :: arg
     shapefunc1_0=one-(arg/sigma)**lambda+half*(arg/sigma)**(2*lambda)-(arg/sigma)**(3*lambda)/6._dp
   end function shapefunc1_0
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/shapefunc2
!  shapefunc2 is g(x) (sinc2)
   function shapefunc2(arg)

     real(dp) :: shapefunc2
     real(dp),intent(in) :: arg
     shapefunc2=(sin(pi_over_rshp*arg)/(pi_over_rshp*arg))**2
   end function shapefunc2
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/shapefunc2_0
!  shapefunc2_0 is g(x) (sinc2) for small x
   function shapefunc2_0(arg)

     real(dp) :: shapefunc2_0
     real(dp),intent(in) :: arg
     shapefunc2_0=one-(pi_over_rshp*arg)**2/three+two*(pi_over_rshp*arg)**4/45._dp
   end function shapefunc2_0
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/shapefunc3
!  shapefunc3 is g(x) (Bessel)
   function shapefunc3(jbes1,jbes2,argl)

     integer,intent(in) :: argl
     real(dp) :: shapefunc3
     real(dp),intent(in) :: jbes1,jbes2
     shapefunc3= alpha(1,1+argl)*jbes1+alpha(2,1+argl)*jbes2
   end function shapefunc3
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/dshpfunc1
!  dshpfunc1(x) is g_prime(x) (gaussian)
   function dshpfunc1(arg)

     real(dp) :: dshpfunc1
     real(dp),intent(in) :: arg
     dshpfunc1=-lambda/sigma*(arg/sigma)**(lambda-1)*exp(-(arg/sigma)**lambda)
   end function dshpfunc1
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/dshpfunc1_ovr_0
!  dshpfunc1_ovr_0(x) is g_prime(x)/x (gaussian) for small x and lambda>2
   function dshpfunc1_ovr_0(arg)

     real(dp) :: dshpfunc1_ovr_0
     real(dp),intent(in) :: arg
     dshpfunc1_ovr_0=-lambda/sigma**2*((arg/sigma)**(lambda-2)-(arg/sigma)**(2*lambda-2))
   end function dshpfunc1_ovr_0
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/dshpfunc1_ovr_0_2
!  dshpfunc1_ovr_0_2(x) is g_prime(x)/x (gaussian) for small x and lambda=2
   function dshpfunc1_ovr_0_2(arg)

     real(dp) :: dshpfunc1_ovr_0_2
     real(dp),intent(in) :: arg
     dshpfunc1_ovr_0_2=-two/sigma**2*(one-(arg/sigma)**2+half*(arg/sigma)**4)
   end function dshpfunc1_ovr_0_2
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/dshpfunc2
!  dshpfunc2(x) is g_prime(x) (sinc2)
   function dshpfunc2(arg)

     real(dp) :: dshpfunc2
     real(dp),intent(in) :: arg
     dshpfunc2=two*pi_over_rshp*sin(pi_over_rshp*arg)/(pi_over_rshp*arg)**3&
&              *(pi_over_rshp*arg*cos(pi_over_rshp*arg)-sin(pi_over_rshp*arg))
   end function dshpfunc2
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/dshpfunc2_ovr_0
!  dshpfunc2_ovr_0(x) is g_prime(x)/x (sinc2) for small x
   function dshpfunc2_ovr_0(arg)

     real(dp) :: dshpfunc2_ovr_0
     real(dp),intent(in) :: arg
     dshpfunc2_ovr_0=-two*pi_over_rshp**2/3._dp+8._dp*pi_over_rshp**4*arg**2/45._dp
   end function dshpfunc2_ovr_0
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/dshpfunc3
!  dshpfunc3(x) is g_prime(x) (Bessel)
   function dshpfunc3(jbesp1,jbesp2,argl)

     integer :: argl
     real(dp) :: dshpfunc3
     real(dp),intent(in) :: jbesp1,jbesp2
     dshpfunc3= alpha(1,1+argl)*qq(1,1+argl)*jbesp1 &
&                              +alpha(2,1+argl)*qq(2,1+argl)*jbesp2
   end function dshpfunc3
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc1
!  d2shpfunc1(x) is g_prime_prime(x) (gaussian)
   function d2shpfunc1(arg)

     real(dp) :: d2shpfunc1
     real(dp),intent(in) :: arg
     d2shpfunc1=lambda/(sigma**2)*(lambda*(arg/sigma)**(2*lambda-2) &
&               -(lambda-1)*(arg/sigma)**(lambda-2))*exp(-(arg/sigma)**lambda)
   end function d2shpfunc1
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc1_ovr2_0
!  d2shpfunc1_ovr2_0(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda>4
   function d2shpfunc1_ovr2_0(arg)

     real(dp) :: d2shpfunc1_ovr2_0
     real(dp),intent(in) :: arg
     d2shpfunc1_ovr2_0=-lambda/(sigma**4)*((lambda-2)*(arg/sigma)**(lambda-4) &
&                                          -(lambda-1)*two*(arg/sigma)**(2*lambda-4))
   end function d2shpfunc1_ovr2_0
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc1_ovr2_0_2
!  d2shpfunc1_ovr2_0_2(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda==2
   function d2shpfunc1_ovr2_0_2(arg)

     real(dp) :: d2shpfunc1_ovr2_0_2
     real(dp),intent(in) :: arg
     d2shpfunc1_ovr2_0_2=four/(sigma**4)*(one-(arg/sigma)**2)
   end function d2shpfunc1_ovr2_0_2
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc1_ovr2_0_3
!  d2shpfunc1_ovr2_0_3(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda==3
   function d2shpfunc1_ovr2_0_3(arg)

     real(dp) :: d2shpfunc1_ovr2_0_3
     real(dp),intent(in) :: arg
     d2shpfunc1_ovr2_0_3=-three/arg/sigma**3+12._dp*arg**2/sigma**6-half*21._dp*arg**5/sigma**9
   end function d2shpfunc1_ovr2_0_3
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc1_ovr2_0_4
!  d2shpfunc1_ovr2_0_4(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (gaussian) for small x and lambda==4
   function d2shpfunc1_ovr2_0_4(arg)

     real(dp) :: d2shpfunc1_ovr2_0_4
     real(dp),intent(in) :: arg
     d2shpfunc1_ovr2_0_4=-8._dp/(sigma**4)*(one-three*(arg/sigma)**4)
   end function d2shpfunc1_ovr2_0_4
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc2
!  d2shpfunc2(x) is g_prime_prime(x) (sinc2)
   function d2shpfunc2(arg)

     real(dp) :: d2shpfunc2
     real(dp),intent(in) :: arg
     d2shpfunc2=two/(pi_over_rshp**2*arg**4)* &
&               (pi_over_rshp**2*arg**2*(cos(pi_over_rshp*arg))**2 &
&               +(three-pi_over_rshp**2*arg**2)*(sin(pi_over_rshp*arg))**2 &
&               -four*pi_over_rshp*arg*cos(pi_over_rshp*arg)*sin(pi_over_rshp*arg))
   end function d2shpfunc2
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc2_ovr2_0
!  d2shpfunc2_ovr2_0(x) is (g_prime_prime(x)-g_prime(x)/x)/x**2 (sinc2) for small x
   function d2shpfunc2_ovr2_0(arg)

     real(dp) :: d2shpfunc2_ovr2_0
     real(dp),intent(in) :: arg
     d2shpfunc2_ovr2_0=16._dp/45._dp*pi_over_rshp**4-8._dp/105._dp*pi_over_rshp**6*arg**2 &
&                      +41._dp/6300._dp*pi_over_rshp**8*arg**4
   end function d2shpfunc2_ovr2_0
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/d2shpfunc3
!  d2shpfunc3(x) is g_prime_prime(x) (Bessel)
   function d2shpfunc3(jbespp1,jbespp2,argl)

     integer,intent(in) :: argl
     real(dp) :: d2shpfunc3
     real(dp),intent(in) :: jbespp1,jbespp2
     d2shpfunc3= alpha(1,1+argl)*(qq(1,1+argl)**2)*jbespp1 &
&                                 +alpha(2,1+argl)*(qq(2,1+argl)**2)*jbespp2
   end function d2shpfunc3
! ------------------------------------------------

end subroutine pawgylm
!!***

!----------------------------------------------------------------------

!!****f* m_paw_finegrid/pawgylmg
!! NAME
!! pawgylmg
!!
!! FUNCTION
!! PAW: Compute Fourier transform of each g_l(r).Y_lm(r) function
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,npw)=integer coordinates of planewaves in basis sphere for this k point.
!!  kpg(npw,nkpg)= (k+G) components (only if useylm=1)
!!  kpt(3)=reduced coordinates of k point
!!  lmax=1+max. value of l angular momentum
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere
!!  ntypat=number of types of atoms
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ylm(npw,lmax**2)=real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  gylmg(npw,lmax**2,ntypat)=Fourier transform of each g_l(r).Y_lm(r) function
!!
!! PARENTS
!!      m_suscep_stat
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawgylmg(gprimd,gylmg,kg,kpg,kpt,lmax,nkpg,npw,ntypat,pawtab,ylm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,nkpg,npw,ntypat
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gprimd(3,3),kpg(npw,nkpg),kpt(3)
 real(dp),intent(in) :: ylm(npw,lmax**2)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

 real(dp),intent(out) :: gylmg(npw,lmax**2,ntypat)

!Local variables-------------------------------
!scalars
 integer :: ig,ilm,itypat,ll,l0,mm,mqgrid
 real(dp) :: kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3
!arrays
 real(dp),allocatable :: glg(:),qgrid(:),kpgnorm(:),shpf(:,:),work(:)

! *************************************************************************

!Get |k+G|:
 LIBPAW_ALLOCATE(kpgnorm,(npw))
 if (nkpg<3) then
   do ig=1,npw
     kpg1=kpt(1)+dble(kg(1,ig));kpg2=kpt(2)+dble(kg(2,ig));kpg3=kpt(3)+dble(kg(3,ig))
     kpgc1=kpg1*gprimd(1,1)+kpg2*gprimd(1,2)+kpg3*gprimd(1,3)
     kpgc2=kpg1*gprimd(2,1)+kpg2*gprimd(2,2)+kpg3*gprimd(2,3)
     kpgc3=kpg1*gprimd(3,1)+kpg2*gprimd(3,2)+kpg3*gprimd(3,3)
     kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
   end do
 else
   do ig=1,npw
     kpgc1=kpg(ig,1)*gprimd(1,1)+kpg(ig,2)*gprimd(1,2)+kpg(ig,3)*gprimd(1,3)
     kpgc2=kpg(ig,1)*gprimd(2,1)+kpg(ig,2)*gprimd(2,2)+kpg(ig,3)*gprimd(2,3)
     kpgc3=kpg(ig,1)*gprimd(3,1)+kpg(ig,2)*gprimd(3,2)+kpg(ig,3)*gprimd(3,3)
     kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
   end do
 end if

 LIBPAW_ALLOCATE(glg,(npw))
 LIBPAW_ALLOCATE(work,(npw))

!Loop over types of atoms
 do itypat=1,ntypat

   mqgrid=pawtab(itypat)%mqgrid_shp
   LIBPAW_ALLOCATE(qgrid,(mqgrid))
   LIBPAW_ALLOCATE(shpf,(mqgrid,2))
   qgrid(1:mqgrid)=pawtab(itypat)%qgrid_shp(1:mqgrid)

!  Loops over (l,m) values
   do ll=0,pawtab(itypat)%lcut_size-1
     l0=ll**2+ll+1

     shpf(1:mqgrid,1:2)=pawtab(itypat)%shapefncg(1:mqgrid,1:2,1+ll)
     call paw_uniform_splfit(qgrid,work,shpf,0,kpgnorm,glg,mqgrid,npw)

     do mm=-ll,ll
       ilm=l0+mm

       gylmg(1:npw,ilm,itypat)=ylm(1:npw,ilm)*glg(1:npw)

!      End loops over (l,m) values
     end do
   end do

!  End loop over atom types
   LIBPAW_DEALLOCATE(qgrid)
   LIBPAW_DEALLOCATE(shpf)
 end do

 LIBPAW_DEALLOCATE(kpgnorm)
 LIBPAW_DEALLOCATE(glg)
 LIBPAW_DEALLOCATE(work)

end subroutine pawgylmg
!!***

!----------------------------------------------------------------------

!!****f* m_paw_finegrid/pawrfgd_fft
!! NAME
!! pawrfgd_fft
!!
!! FUNCTION
!! Determine each point of the (fine) rectangular grid
!! around a given atom and compute r-R vectors.
!! R is the position of the atom.
!!
!! INPUTS
!!  [fft_distrib(n3)]= (optional) index of processes which own fft planes in 3rd dimension
!!  [fft_index(n3)]= (optional) local fft indexes for current process
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  [me_fft]= (optional) my rank in the FFT MPI communicator
!!  n1,n2,n3= sizes of the FFT grid (entire simulation cell)
!!  rcut= radius of the sphere around the atom
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol= unit cell volume
!!  xred(3)= reduced coordinates of the atom
!!
!! OUTPUT
!!  ifftsph(nfgd)= FFT index (fine grid) of the points in the sphere around current atom
!!  nfgd= number of points in the sphere around current atom
!!  rfgd(3,nfgd)= cartesian coordinates of r-R.
!!
!! PARENTS
!!      m_paw_dfpt,m_paw_nhat
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrfgd_fft(ifftsph,gmet,n1,n2,n3,nfgd,rcut,rfgd,rprimd,ucvol,xred, &
&                      fft_distrib,fft_index,me_fft) ! optional arguments

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3
 integer,intent(out) :: nfgd
 integer,optional,intent(in) :: me_fft
 real(dp),intent(in) :: rcut,ucvol
!arrays
 integer,target,optional,intent(in) :: fft_distrib(n3),fft_index(n3)
 integer,allocatable,intent(out) :: ifftsph(:)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),xred(3)
 real(dp),allocatable,intent(out) :: rfgd(:,:)

!Local variables ------------------------------
!scalars
 integer,parameter :: ishift=15
 integer :: i1,i2,i3,ifft_local,ix,iy,iz,izloc,me_fft_,n1a,n1b,n2a,n2b,n3a,n3b,ncmax
 real(dp),parameter :: delta=0.99_dp
 real(dp) :: difx,dify,difz,rr1,rr2,rr3,r2,r2cut,rx,ry,rz
 character(len=500) :: msg
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 integer,pointer :: fft_distrib_(:),fft_index_(:)
 real(dp),allocatable :: rfgd_tmp(:,:)

! *************************************************************************

!Define a "box" around the atom
 r2cut=1.0000001_dp*rcut**2
 rr1=sqrt(r2cut*gmet(1,1))
 rr2=sqrt(r2cut*gmet(2,2))
 rr3=sqrt(r2cut*gmet(3,3))
 n1a=int((xred(1)-rr1+ishift)*n1+delta)-ishift*n1
 n1b=int((xred(1)+rr1+ishift)*n1      )-ishift*n1
 n2a=int((xred(2)-rr2+ishift)*n2+delta)-ishift*n2
 n2b=int((xred(2)+rr2+ishift)*n2      )-ishift*n2
 n3a=int((xred(3)-rr3+ishift)*n3+delta)-ishift*n3
 n3b=int((xred(3)+rr3+ishift)*n3      )-ishift*n3

!Get the distrib associated with this fft_grid
 if (present(fft_distrib).and.present(fft_index).and.present(me_fft)) then
   me_fft_=me_fft ; fft_distrib_ => fft_distrib ; fft_index_ => fft_index
 else
   me_fft_=0
   LIBPAW_POINTER_ALLOCATE(fft_distrib_,(n3))
   LIBPAW_POINTER_ALLOCATE(fft_index_,(n3))
   fft_distrib_=0;fft_index_=(/(i3,i3=1,n3)/)
 end if

!Temporary allocate "large" arrays
 ncmax=1+int(1.5_dp*(n1*n2*n3)*four_pi/(three*ucvol)*rcut**3)
 LIBPAW_ALLOCATE(ifftsph_tmp,(ncmax))
 LIBPAW_ALLOCATE(rfgd_tmp,(3,ncmax))

!Set number of points to zero
 nfgd=0

!Loop over FFT points
 do i3=n3a,n3b
   iz=mod(i3+ishift*n3,n3)
   if (fft_distrib_(iz+1)==me_fft_) then
     izloc=fft_index_(iz+1) - 1
     difz=dble(i3)/dble(n3)-xred(3)
     do i2=n2a,n2b
       iy=mod(i2+ishift*n2,n2)
       dify=dble(i2)/dble(n2)-xred(2)
       do i1=n1a,n1b
         ix=mod(i1+ishift*n1,n1)
         difx=dble(i1)/dble(n1)-xred(1)

!        Compute r-R
         rx=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
         ry=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
         rz=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
         r2=rx**2+ry**2+rz**2

!        Select matching points
         if (r2 <= r2cut) then
           ifft_local=1+ix+n1*(iy+n2*izloc)
           if (ifft_local>0) then
             nfgd=nfgd+1
             if (nfgd>ncmax) then
               msg='Number of fft points around atom exceeds max. allowed!'
               MSG_BUG(msg)
             end if
             rfgd_tmp(1,nfgd)=rx
             rfgd_tmp(2,nfgd)=ry
             rfgd_tmp(3,nfgd)=rz
             ifftsph_tmp(nfgd)=ifft_local
           end if
         end if

!      End of loops
       end do
     end do
   end if
 end do

!Now fill output arrays
 if (allocated(ifftsph)) then
   LIBPAW_DEALLOCATE(ifftsph)
 end if
 if (allocated(rfgd)) then
   LIBPAW_DEALLOCATE(rfgd)
 end if
 LIBPAW_ALLOCATE(ifftsph,(nfgd))
 LIBPAW_ALLOCATE(rfgd,(3,nfgd))
 ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)
 rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)

!Release temporary memory
 LIBPAW_DEALLOCATE(ifftsph_tmp)
 LIBPAW_DEALLOCATE(rfgd_tmp)
 if (.not.present(fft_distrib).or..not.present(fft_index)) then
   LIBPAW_POINTER_DEALLOCATE(fft_distrib_)
   LIBPAW_POINTER_DEALLOCATE(fft_index_)
 end if

end subroutine pawrfgd_fft
!!***

!----------------------------------------------------------------------

!!****f* m_paw_finegrid/pawrfgd_wvl
!! NAME
!! pawrfgd_wvl
!!
!! FUNCTION
!! Determine each point of the (fine) rectangular grid
!! around a given atom and compute r-R vectors.
!! R is the position of the atom.
!!
!! INPUTS
!!  geocode= code for geometry (boundary conditions)
!!  hh(3)=fine grid spacing
!!  i3s= TO BE COMPLETED
!!  n1,n2,n3= TO BE COMPLETED
!!  n1i,n2i,n3pi= TO BE COMPLETED
!!  rcut= radius of the sphere around the atom
!!  rloc= cut-off radius for local psp?
!!  shift= TO BE COMPLETED
!!  xred(3)= cartesian coordinates of the atom
!!
!! OUTPUT
!!  ifftsph(nfgd)= FFT index (fine grid) of the points in the sphere around current atom
!!  nfgd= number of points in the sphere around current atom
!!  rfgd(3,nfgd)= cartesian coordinates of r-R.
!!
!! PARENTS
!!      m_paw_nhat
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawrfgd_wvl(geocode,hh,ifftsph,i3s,n1,n1i,n2,n2i,n3,n3pi,&
&                      nfgd,rcut,rloc,rfgd,shift,xcart)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: i3s,n1,n1i,n2,n2i,n3,n3pi,shift
 integer,intent(out) :: nfgd
 real(dp),intent(in) :: rcut,rloc
 character(1),intent(in) :: geocode
!arrays
 integer,allocatable,intent(out) :: ifftsph(:)
 real(dp),intent(in) :: hh(3),xcart(3)
 real(dp),allocatable,intent(out) :: rfgd(:,:)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,iex,iey,iez,ind,isx,isy,isz,j1,j2,j3
 integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,ncmax
 logical :: gox,goy,goz,perx,pery,perz
 real(dp) :: cutoff,r2,r2cut,rx,ry,rz,xx,yy,zz
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 real(dp),allocatable :: rfgd_tmp(:,:)

! *************************************************************************

!Data for periodicity in the three directions
 perx=(geocode/='F')
 pery=(geocode=='P')
 perz=(geocode/='F')
 call my_ext_buffers(perx,nbl1,nbr1)
 call my_ext_buffers(pery,nbl2,nbr2)
 call my_ext_buffers(perz,nbl3,nbr3)

!Define a "box" around the atom
 cutoff=10.d0*rloc
 r2cut=1.0000001_dp*rcut**2
 rx=xcart(1)
 ry=xcart(2)
 rz=xcart(3)
 isx=floor((rx-cutoff)/hh(1))
 isy=floor((ry-cutoff)/hh(2))
 isz=floor((rz-cutoff)/hh(3))
 iex=ceiling((rx+cutoff)/hh(1))
 iey=ceiling((ry+cutoff)/hh(2))
 iez=ceiling((rz+cutoff)/hh(3))

!Temporary allocate "large" arrays
!  use factor 1+int(1.1*, for safety reasons
 ncmax=1
 if (n3pi>0) ncmax=1+int((rcut/hh(1)+1.0)*(rcut/hh(2)+1.0)*(rcut/hh(3)+1.0)*four_pi/three)
 LIBPAW_ALLOCATE(ifftsph_tmp,(ncmax))
 LIBPAW_ALLOCATE(rfgd_tmp,(3,ncmax))

!Set number of points to zero
 nfgd=0

!Loop over WVL points
 do i3=isz,iez
   zz=real(i3,kind=8)*hh(3)-rz
   call my_ind_positions(perz,i3,n3,j3,goz)
   j3=j3+nbl3+1
   do i2=isy,iey
     yy=real(i2,kind=8)*hh(2)-ry
     call my_ind_positions(pery,i2,n2,j2,goy)
     do i1=isx,iex
       xx=real(i1,kind=8)*hh(1)-rx
       call my_ind_positions(perx,i1,n1,j1,gox)
       r2=xx**2+yy**2+zz**2
       if (j3>=i3s.and.j3<=i3s+n3pi-1.and.goy.and.gox) then

!        Select matching points
         if (r2<=r2cut) then
           ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
           nfgd=nfgd+1
           rfgd_tmp(:,nfgd)=[xx,yy,zz]
           ifftsph_tmp(nfgd)=shift+ind
         end if

!      End of loops
       end if
     end do
   end do
 end do

!Now fill output arrays
 if (allocated(ifftsph)) then
   LIBPAW_DEALLOCATE(ifftsph)
 end if
 if (allocated(rfgd)) then
   LIBPAW_DEALLOCATE(rfgd)
 end if
 LIBPAW_ALLOCATE(ifftsph,(nfgd))
 LIBPAW_ALLOCATE(rfgd,(3,nfgd))
 ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)
 rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)

!Release temporary memory
 LIBPAW_DEALLOCATE(ifftsph_tmp)
 LIBPAW_DEALLOCATE(rfgd_tmp)

!*********************************************************************
!Small functions related to boundary conditions
 contains
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/my_ind_positions
   subroutine my_ind_positions(periodic,i,n,j,go)

     integer,intent(in) :: i,n
     logical,intent(in) :: periodic
     integer,intent(out) :: j
     logical,intent(out) :: go
     if (periodic) then
       j=modulo(i,2*n+2) ; go=.true.
     else
       j=i ; go=(i>=-14.and.i<=2*n+16)
     end if
   end subroutine my_ind_positions
!!***
! ------------------------------------------------
!!****f* m_paw_finegrid/my_ext_buffers
   subroutine my_ext_buffers(periodic,nl,nr)

     logical, intent(in) :: periodic
     integer, intent(out) :: nl,nr
     if (periodic) then
       nl=0 ; nr=0
     else
       nl=14 ; nr=15
     end if
   end subroutine my_ext_buffers
! ------------------------------------------------

end subroutine pawrfgd_wvl
!!***

!----------------------------------------------------------------------

!!****f* m_paw_finegrid/pawexpiqr
!! NAME
!! pawexpiqr
!!
!! FUNCTION
!! Compute exp(i.q.r) for each point of the (fine) rectangular grid
!! around a given atomic site. R is the position of the atom.
!! Used for the determination of phonons at non-zero q wavevector.
!!
!! INPUTS
!!  gprimd(3,3)= dimensional primitive translations for reciprocal space
!!  nfgd= number of (fine grid) FFT points in the paw sphere around current atom
!!  qphon(3)= wavevector of the phonon
!!  rfgd(3,nfgd)= coordinates of r-R on the fine grid around current atom
!!  xred(3)= reduced atomic coordinates
!!
!! OUTPUT
!!  expiqr(2,nfgd)= exp(i.q.r) around the current atom
!!                                 Not allocated if q=0 !
!!
!! PARENTS
!!      m_paw_dfpt,m_paw_nhat,m_pawdij,m_respfn_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawexpiqr(expiqr,gprimd,nfgd,qphon,rfgd,xred)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nfgd
!arrays
 real(dp),intent(in) :: gprimd(3,3),qphon(3),xred(3)
 real(dp),intent(in) :: rfgd(:,:)
 real(dp),intent(out) :: expiqr(2,nfgd)

!Local variables ------------------------------
!scalars
 integer :: ic
 logical :: qne0
 real(dp) :: phase,phase_xred,qx,qy,qz
 character(len=500) :: msg
!arrays

! *************************************************************************

 if (size(rfgd)/=3*nfgd) then
   msg='rfgd array must be allocated!'
   MSG_BUG(msg)
 end if

 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)

!Compute q in cartesian coordinates
 if (qne0) then
   qx=gprimd(1,1)*qphon(1)+gprimd(1,2)*qphon(2)+gprimd(1,3)*qphon(3)
   qy=gprimd(2,1)*qphon(1)+gprimd(2,2)*qphon(2)+gprimd(2,3)*qphon(3)
   qz=gprimd(3,1)*qphon(1)+gprimd(3,2)*qphon(2)+gprimd(3,3)*qphon(3)
   phase_xred=two_pi*(qphon(1)*xred(1)+qphon(2)*xred(2)+qphon(3)*xred(3))
 end if

!Compute exp(i.q.r)
 if (qne0) then
   do ic=1,nfgd
     phase=two_pi*(qx*rfgd(1,ic)+qy*rfgd(2,ic)+qz*rfgd(3,ic)) + phase_xred
     expiqr(1,ic)=cos(phase)
     expiqr(2,ic)=sin(phase)
   end do
 end if

end subroutine pawexpiqr
!!***

!----------------------------------------------------------------------

END MODULE m_paw_finegrid
!!***

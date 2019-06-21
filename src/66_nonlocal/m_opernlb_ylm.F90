!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_opernlb_ylm
!! NAME
!!  m_opernlb_ylm
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (MT)
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

module m_opernlb_ylm

 use defs_basis
 use m_abicore
 use m_errors
#if defined HAVE_OPENMP
 use OMP_LIB
#endif

 implicit none

 private
!!***

 public :: opernlb_ylm
!!***

contains
!!***

!!****f* ABINIT/opernlb_ylm
!! NAME
!! opernlb_ylm
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   from projected scalars to reciprocal space.
!! * Operate with the non-local projectors and the overlap matrix,
!!   from projected scalars to reciprocal space.
!!
!! INPUTS
!!  choice=chooses possible output (see below)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_dgxdt(ndgxdt_fac) = used only when cplex = 1
!!    cplex_dgxdt(i)=1 if dgxdt(1,i,:,:) is real, 2 if it is pure imaginary
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Vnl (NL operator)
!!  dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfacrelated to Sij (overlap)
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))= reduced projected scalars related to Sij (overlap)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2) or (choice=22,signs=2)
!!                        - k point direction in the case (choice=5, 51, 52 and signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!                        - strain component (1:9) in the case (choice=33,signs=2) 
!!                        - (1:9) components to specify the atom to be moved and the second q-gradient 
!!                          direction in the case (choice=25,signs=2)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  kpg(npw,nkpg)=(k+G) components (if nkpg=3).
!!                (k+G) Cartesian components for choice=33
!!  matblk=dimension of the array ph3d
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0 or 3)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  [qdir]= optional, direction of the q-gradient (only for choice=22, choice=25 and choice=33)
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! --if (paw_opt=0)
!!    vectout(2,npwout*my_nspinor*ndat)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.  
!!      if (choice=22) <G|d2V_nonlocal/d(atm. pos)dq|vect_in> (at q=0)
!!      if (choice=25) <G|d3V_nonlocal/d(atm. pos)dqdq|vect_in> (at q=0)
!!      if (choice=33) <G|d2V_nonlocal/d(strain)dq|vect_in> (at q=0)
!! --if (paw_opt=0, 1 or 4)
!!    vect(2,npwout*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1)  <G|V_nonlocal|vect_in>
!!      if (choice=2)  <G|dV_nonlocal/d(atm. pos)|vect_in>
!!      if (choice=3)  <G|dV_nonlocal/d(strain)|vect_in>
!!      if (choice=5)  <G|dV_nonlocal/d(k)|vect_in>
!!      if (choice=51) <G|d(right)V_nonlocal/d(k)|vect_in>
!!      if (choice=52) <G|d(left)V_nonlocal/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)V_nonlocal/d(k)|vect_in>
!!      if (choice=54) <G|d[d(right)V_nonlocal/d(k)]/d(atm. pos)|vect_in>
!!      if (choice=8)  <G|d2V_nonlocal/d(k)d(k)|vect_in>
!!      if (choice=81) <G|d[d(right)V_nonlocal/d(k)]/d(k)|vect_in>
!!  if (paw_opt=2)
!!    vect(2,npwout*nspinor)=final vector in reciprocal space:
!!      if (choice=1)  <G|V_nonlocal-lamdba.(I+S)|vect_in> (note: not including <G|I|c>)
!!      if (choice=2)  <G|d[V_nonlocal-lamdba.(I+S)]/d(atm. pos)|vect_in>
!!      if (choice=3)  <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_in>
!!      if (choice=5)  <G|d[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=54) <G|d[d(right)V_nonlocal/d(k)]/d(atm. pos)|vect_in>
!!      if (choice=8)  <G|d2[V_nonlocal-lamdba.(I+S)]/d(k)d(k)|vect_in>
!!      if (choice=81) <G|d[d(right[V_nonlocal-lamdba.(I+S)]/d(k)]/d(k)|vect_in>
!! --if (paw_opt=3 or 4)
!!    svect(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1)  <G|I+S|vect_in> (note: not including <G|I|c>)
!!      if (choice=2)  <G|dS/d(atm. pos)|vect_in>
!!      if (choice=3)  <G|dS/d(strain)|vect_in>
!!      if (choice=5)  <G|dS/d(k)|vect_in>
!!      if (choice=51) <G|d(right)S/d(k)|vect_in>
!!      if (choice=52) <G|d(left)S/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)S/d(k)|vect_in>
!!      if (choice=54) <G|d[d(right)V_nonlocal/d(k)]/d(atm. pos)|vect_in>
!!      if (choice=7)  <G|sum_i[p_i><p_i]|vect_in>
!!      if (choice=8)  <G|d2S/d(k)d(k)|vect_in>
!!      if (choice=81) <G|d[d(right)S/d(k)]/d(k)|vect_in>
!!
!! NOTES
!! 1-The openMP version is different from the standard version:
!!   the standard version is more effifient on one CPU core.
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!!
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE

subroutine opernlb_ylm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
&                      d2gxdtfac,d2gxdtfac_sij,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
&                      ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nd2gxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
&                      nspinor,paw_opt,ph3d,svect,ucvol,vect,qdir)

 implicit none

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,idir,matblk,ndgxdtfac,nd2gxdtfac,nincat
 integer,intent(in) :: nkpg,nlmn,npw,nspinor,paw_opt
 integer,intent(in),optional :: qdir
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) ::  cplex_dgxdt(ndgxdtfac),cplex_d2gxdt(nd2gxdtfac),indlmn(6,nlmn),nloalg(3)
 real(dp),intent(in) :: d2gxdtfac(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
 real(dp),intent(in) :: d2gxdtfac_sij(cplex,nd2gxdtfac,nlmn,nincat*(paw_opt/3),nspinor)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: svect(:,:),vect(:,:)
!Local variables-------------------------------
!Arrays
!scalars
 integer :: fdb,fdf,ia,ialpha,iaph3d,ibeta,ic,idelta,idelgam,igamma
 integer :: ii,il,ilmn,ipw,ipwshft,ispinor,jc,nthreads,ffnl_dir1,ffnl_dir(3)
 real(dp) :: scale,two_piinv,wt
 logical :: parity
!arrays
 integer,parameter :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,parameter :: idir1(9)=(/1,1,1,2,2,2,3,3,3/),idir2(9)=(/1,2,3,1,2,3,1,2,3/)
 integer,parameter :: nalpha(9)=(/1,2,3,3,3,2,2,1,1/),nbeta(9)=(/1,2,3,2,1,1,3,3,2/)
 real(dp),allocatable :: d2gxdtfac_(:,:,:),d2gxdtfacs_(:,:,:),dgxdtfac_(:,:,:),dgxdtfacs_(:,:,:),gxfac_(:,:),gxfacs_(:,:)
! real(dp),allocatable :: kpg(:,:)
 complex(dpc),allocatable :: ztab(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Nothing to do when choice=4, 6 or 23
 if (choice==4.or.choice==6.or.choice==23) return

!DDK not compatible with istwkf > 1
 if(cplex==1.and.(any(cplex_dgxdt(:)==2).or.any(cplex_d2gxdt(:)==2)))then
   MSG_BUG("opernlb_ylm+ddk not compatible with istwfk>1")
 end if

!Inits
 wt=four_pi/sqrt(ucvol)
 nthreads=1
#if defined HAVE_OPENMP
 nthreads=OMP_GET_NUM_THREADS()
#endif

 if (paw_opt/=3) then
   ABI_ALLOCATE(gxfac_,(2,nlmn))
   gxfac_(:,:)=zero
   if (choice>1) then
     ABI_ALLOCATE(dgxdtfac_,(2,ndgxdtfac,nlmn))
     if(ndgxdtfac>0) dgxdtfac_(:,:,:)=zero
   end if
   if (choice==54.or.choice==8.or.choice==81.or.choice==33) then
     ABI_ALLOCATE(d2gxdtfac_,(2,nd2gxdtfac,nlmn))
     if(nd2gxdtfac>0) d2gxdtfac_(:,:,:)=zero
   end if
 end if
 if (paw_opt>=3) then
   ABI_ALLOCATE(gxfacs_,(2,nlmn))
   gxfacs_(:,:)=zero
   if (choice>1) then
     ABI_ALLOCATE(dgxdtfacs_,(2,ndgxdtfac,nlmn))
     if (ndgxdtfac>0) dgxdtfacs_(:,:,:)=zero
   end if
   if (choice==54.or.choice==8.or.choice==81) then
     ABI_ALLOCATE(d2gxdtfacs_,(2,nd2gxdtfac,nlmn))
     if (nd2gxdtfac>0) d2gxdtfacs_(:,:,:)=zero
   end if
 end if

#ifdef MR_DEV
if (choice==33) two_piinv=1.0_dp/two_pi
#endif

 ABI_ALLOCATE(ztab,(npw))

!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================
 if (nthreads==1) then

!  Loop on spinorial components
   do ispinor=1,nspinor
     ipwshft=(ispinor-1)*npw

!    Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(2)>0) iaph3d=ia+ia3-1
!      Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxfac_(2,ilmn)=zero
           else
             gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
             else
               gxfac_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfac_(1,ii,ilmn)= scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
                   dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic =  cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfac_(jc,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfac_(jc,ii,ilmn)= scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
         if (choice==54.or.choice==8.or.choice==81.or.choice==33) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 d2gxdtfac_(1:cplex_fac,1:nd2gxdtfac,ilmn)=scale*d2gxdtfac(1:cplex_fac,1:nd2gxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,nd2gxdtfac
                   ic = cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfac_(ic,ii,ilmn)=scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                   d2gxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,nd2gxdtfac
                   d2gxdtfac_(1,ii,ilmn)= scale*d2gxdtfac(2,ii,ilmn,ia,ispinor)
                   d2gxdtfac_(2,ii,ilmn)=-scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,nd2gxdtfac
                   ic =  cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     d2gxdtfac_(jc,ii,ilmn)=-scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     d2gxdtfac_(jc,ii,ilmn)= scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
       end if

!      Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       if (paw_opt>=3) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
             if (cplex==1) gxfacs_(2,ilmn)=zero
           else
             gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
             if (cplex==2) then
               gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
             else
               gxfacs_(1,ilmn)=zero
             end if
           end if
         end do
         if (choice>1) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfacs_(1,ii,ilmn)= scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfacs_(jc,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfacs_(jc,ii,ilmn)= scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
         if (choice==54.or.choice==8.or.choice==81) then
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 d2gxdtfacs_(1:cplex,1:nd2gxdtfac,ilmn)=scale*d2gxdtfac_sij(1:cplex,1:nd2gxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,nd2gxdtfac
                   ic = cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfacs_(ic,ii,ilmn)=scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   d2gxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,nd2gxdtfac
                   d2gxdtfacs_(1,ii,ilmn)= scale*d2gxdtfac_sij(2,ii,ilmn,ia,ispinor)
                   d2gxdtfacs_(2,ii,ilmn)=-scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,nd2gxdtfac
                   ic = cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     d2gxdtfacs_(jc,ii,ilmn)=-scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     d2gxdtfacs_(jc,ii,ilmn)= scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
         end if
       end if

!      Compute <g|Vnl|c> (or derivatives) for each plane wave:

       if (paw_opt/=3) then

         ztab(:)=czero

!        ------
         if (choice==1) then ! <g|Vnl|c>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==2) then ! derivative w.r.t. atm. pos
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir)*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

#ifdef MR_DEV
!        ------
         if (choice==22) then ! mixed derivative w.r.t. atm. pos and q vector (at q=0)
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+qdir
           if (idir==qdir) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end if
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+kpg(:,idir)*ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             ztab(:)=ztab(:)-ffnl(:,ffnl_dir1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
           ztab(:)=ztab(:)*two_pi
         end if
#endif

#ifdef MR_DEV
!        ------
         if (choice==25) then ! mixed derivative w.r.t. atm. pos and two q vectors (at q=0)
           !Use same notation as the notes for clarity
           ialpha=nalpha(idir)
           idelta=nbeta(idir)
           igamma=qdir
           idelgam=gamma(idelta,igamma)
           if (ialpha==igamma) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1+idelta,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end if
           if (ialpha==idelta) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1+igamma,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end if
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+kpg(:,ialpha)*ffnl(:,4+idelgam,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             ztab(:)=ztab(:)-ffnl(:,4+idelgam,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
           ztab(:)=ztab(:)*two_pi
         end if
#endif

!        ------
         if (choice==3) then ! derivative w.r.t. strain
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           if (idir<=3) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
&               *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp)&
&               -ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           else
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
&               -ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end if
         end if

#ifdef MR_DEV
!        ------
         if (choice==33) then ! mixed derivative w.r.t. strain and q vector (at q=0)
           !Use same notation as the notes for clarity
           ibeta=nalpha(idir)
           idelta=nbeta(idir)
           igamma=qdir
           idelgam=gamma(idelta,igamma)
           if (ibeta==igamma) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+onehalf*ffnl(:,1+idelta,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               ztab(:)=ztab(:)+half*ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp)
             end do
           end if
           if (ibeta==idelta) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+onehalf*ffnl(:,1+igamma,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               ztab(:)=ztab(:)+half*ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end if
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+kpg(:,ibeta)*ffnl(:,4+idelgam,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             ztab(:)=ztab(:)+ffnl(:,1+idelta,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)
             ztab(:)=ztab(:)+ffnl(:,1+igamma,ilmn)*cmplx(d2gxdtfac_(1,2,ilmn),d2gxdtfac_(2,2,ilmn),kind=dp)
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(d2gxdtfac_(1,3,ilmn),d2gxdtfac_(2,3,ilmn),kind=dp)
           end do
           ztab(:)=ztab(:)*two_piinv
         end if
#endif

!        ------
         if (choice==5) then ! full derivative w.r.t. k
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1        ,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)&
&             +ffnl(:,ffnl_dir1,ilmn)*cmplx(   gxfac_(1,  ilmn),   gxfac_(2  ,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
!                                                -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) + &
&             ffnl(:,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) - &
&             ffnl(:,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==54) then ! mixed derivative w.r.t. atm. pos and (right) k
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfac_(2,1,ilmn),-dgxdtfac_(1,1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir1(idir))*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==8) then ! full second order derivative w.r.t. k
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=idir1(idir); ffnl_dir(2)=idir2(idir)
           ffnl_dir(3) = gamma(ffnl_dir(1),ffnl_dir(2))
           do ilmn=1,nlmn
             ztab(:)=ztab(:) &
&             +ffnl(:,4+ffnl_dir(3),ilmn)*cmplx(    gxfac_(1,  ilmn),    gxfac_(2,  ilmn),kind=dp)&
&             +ffnl(:,1            ,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)&
&             +ffnl(:,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfac_(1,2,ilmn), dgxdtfac_(2,2,ilmn),kind=dp)&
&             +ffnl(:,1+ffnl_dir(2),ilmn)*cmplx( dgxdtfac_(1,1,ilmn), dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==81) then
           ! partial second order derivative w.r.t. k
           ! full derivative w.r.t. k1, right derivative w.r.t. k2
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=1; if(dimffnl>2) ffnl_dir(1)=idir1(idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) &
&             +ffnl(:,1            ,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)&
&             +ffnl(:,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfac_(1,1,ilmn), dgxdtfac_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)

         vect(1,1+ipwshft:npw+ipwshft)=vect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
         vect(2,1+ipwshft:npw+ipwshft)=vect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))

       end if

!      Compute <g|S|c> (or derivatives) for each plane wave:

       if (paw_opt>=3) then

         ztab(:)=czero

!        ------
         if (choice==1) then ! <g|S|c>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==2) then ! derivative w.r.t. atm. pos
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir)*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==3) then ! derivative w.r.t. strain
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           if (idir<=3) then
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)&
&               *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
&               -ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           else
             do ilmn=1,nlmn
               ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
&               -ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end if
         end if

!        ------
         if (choice==5) then ! full derivative w.r.t. k
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)&
&             +ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==51) then ! right derivative: <G|p>S<dp/dk|psi>
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==52) then ! left derivative: <G|dp/dk>S<p|psi>
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>S<dp/dk_(idir-1)|psi>
!                                                -<G|dp/dk_(idir-1)>S<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) + &
&             ffnl(:,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) - &
&             ffnl(:,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==54) then ! mixed derivative w.r.t. atm. pos and k
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(dgxdtfacs_(2,1,ilmn),-dgxdtfacs_(1,1,ilmn),kind=dp)
           end do
           ztab(:)=two_pi*kpg(:,idir1(idir))*ztab(:)
           do ilmn=1,nlmn
             ztab(:)=ztab(:)+ffnl(:,1,ilmn)*cmplx(d2gxdtfacs_(1,1,ilmn),d2gxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==8) then ! full second order derivative w.r.t. k
           ffnl_dir(1)=idir1(idir); ffnl_dir(2)=idir2(idir)
           ffnl_dir(3) = gamma(ffnl_dir(1),ffnl_dir(2))
           do ilmn=1,nlmn
             ztab(:)=ztab(:) &
&             +ffnl(:,4+ffnl_dir(3),ilmn)*cmplx(    gxfacs_(1,  ilmn),    gxfacs_(2,  ilmn),kind=dp)&
&             +ffnl(:,1            ,ilmn)*cmplx(d2gxdtfacs_(1,1,ilmn),d2gxdtfacs_(2,1,ilmn),kind=dp)&
&             +ffnl(:,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfacs_(1,2,ilmn), dgxdtfacs_(2,2,ilmn),kind=dp)&
&             +ffnl(:,1+ffnl_dir(2),ilmn)*cmplx( dgxdtfacs_(1,1,ilmn), dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if

!        ------
         if (choice==81) then
           ! partial second order derivative w.r.t. k
           ! full derivative w.r.t. k1, right derivative w.r.t. k2
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=1; if(dimffnl>2) ffnl_dir(1)=idir1(idir)
           do ilmn=1,nlmn
             ztab(:)=ztab(:) &
&             +ffnl(:,1            ,ilmn)*cmplx(d2gxdtfacs_(1,1,ilmn),d2gxdtfacs_(2,1,ilmn),kind=dp)&
&             +ffnl(:,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfacs_(1,1,ilmn), dgxdtfacs_(2,1,ilmn),kind=dp)
           end do
         end if


!        ------
         ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
         svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
         svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
       end if

!      End loop on atoms
     end do
   end do !  End loop on spinors


!  ==========================================================================
!  ========== OPENMP VERSION ================================================
!  ==========================================================================
 else

!  Loop on spinorial components
   do ispinor=1,nspinor
     ipwshft=(ispinor-1)*npw

!    Loop on atoms (blocking)
     do ia=1,nincat
       iaph3d=ia;if (nloalg(2)>0) iaph3d=ia+ia3-1

!      Scale gxfac with 4pi/sqr(omega).(-i)^l
       if (paw_opt/=3) then
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfac_(1:cplex_fac,ilmn)=scale*gxfac(1:cplex_fac,ilmn,ia,ispinor)
             if (cplex_fac==1) gxfac_(2,ilmn)=zero
           else
             gxfac_(2,ilmn)=-scale*gxfac(1,ilmn,ia,ispinor)
             if (cplex_fac==2) then
               gxfac_(1,ilmn)=scale*gxfac(2,ilmn,ia,ispinor)
             else
               gxfac_(1,ilmn)=zero
             end if
           end if
         end do
!$OMP END DO
         if (choice>1) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 dgxdtfac_(1:cplex_fac,1:ndgxdtfac,ilmn)=scale*dgxdtfac(1:cplex_fac,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfac_(2,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   dgxdtfac_(1,ii,ilmn)= scale*dgxdtfac(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic =  cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfac_(jc,ii,ilmn)=-scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfac_(jc,ii,ilmn)= scale*dgxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
!$OMP END DO
         end if
         if (choice==54.or.choice==8.or.choice==81.or.choice==33) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex_fac==2)then
                 d2gxdtfac_(1:cplex_fac,1:nd2gxdtfac,ilmn)=scale*d2gxdtfac(1:cplex_fac,1:nd2gxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,nd2gxdtfac
                   ic = cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfac_(ic,ii,ilmn)=scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                   d2gxdtfac_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex_fac==2)then
                 do ii=1,nd2gxdtfac
                   d2gxdtfac_(1,ii,ilmn)= scale*d2gxdtfac(2,ii,ilmn,ia,ispinor)
                   d2gxdtfac_(2,ii,ilmn)=-scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,nd2gxdtfac
                   ic =  cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfac_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     d2gxdtfac_(jc,ii,ilmn)=-scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                   else
                     d2gxdtfac_(jc,ii,ilmn)= scale*d2gxdtfac(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
!$OMP END DO
         end if
!$OMP END PARALLEL
       end if

!      Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
       if (paw_opt>=3) then
!$OMP PARALLEL PRIVATE(ilmn,il,parity,scale,ii,ic,jc)
!$OMP DO
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (parity) then
             gxfacs_(1:cplex,ilmn)=scale*gxfac_sij(1:cplex,ilmn,ia,ispinor)
             if (cplex==1) gxfacs_(2,ilmn)=zero
           else
             gxfacs_(2,ilmn)=-scale*gxfac_sij(1,ilmn,ia,ispinor)
             if (cplex==2) then
               gxfacs_(1,ilmn)=scale*gxfac_sij(2,ilmn,ia,ispinor)
             else
               gxfacs_(1,ilmn)=zero
             end if
           end if
         end do
!$OMP END DO
         if (choice>1) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 dgxdtfacs_(1:cplex,1:ndgxdtfac,ilmn)=scale*dgxdtfac_sij(1:cplex,1:ndgxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,ndgxdtfac
                   dgxdtfacs_(2,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   dgxdtfacs_(1,ii,ilmn)= scale*dgxdtfac_sij(2,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,ndgxdtfac
                   ic = cplex_dgxdt(ii) ; jc = 3-ic
                   dgxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     dgxdtfacs_(jc,ii,ilmn)=-scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     dgxdtfacs_(jc,ii,ilmn)= scale*dgxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
!$OMP END DO
         end if
         if (choice==54.or.choice==8.or.choice==81) then
!$OMP DO
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (parity) then
               if(cplex==2)then
                 d2gxdtfacs_(1:cplex,1:nd2gxdtfac,ilmn)=scale*d2gxdtfac_sij(1:cplex,1:nd2gxdtfac,ilmn,ia,ispinor)
               else
                 do ii=1,nd2gxdtfac
                   ic = cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfacs_(ic,ii,ilmn)=scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   d2gxdtfacs_(jc,ii,ilmn)=zero
                 end do
               end if
             else
               if(cplex==2)then
                 do ii=1,nd2gxdtfac
                   d2gxdtfacs_(1,ii,ilmn)= scale*d2gxdtfac_sij(2,ii,ilmn,ia,ispinor)
                   d2gxdtfacs_(2,ii,ilmn)=-scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                 end do
               else
                 do ii=1,nd2gxdtfac
                   ic = cplex_d2gxdt(ii) ; jc = 3-ic
                   d2gxdtfacs_(ic,ii,ilmn)=zero
                   if(ic==1)then
                     d2gxdtfacs_(jc,ii,ilmn)=-scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   else
                     d2gxdtfacs_(jc,ii,ilmn)= scale*d2gxdtfac_sij(1,ii,ilmn,ia,ispinor)
                   end if
                 end do
               end if
             end if
           end do
!$OMP END DO
         end if
!$OMP END PARALLEL
       end if

!      Compute <g|Vnl|c> (or derivatives) for each plane wave:
       if (paw_opt/=3) then
!$OMP PARALLEL PRIVATE(ipw,ilmn,fdf,fdb,ffnl_dir1)

!        ------
         if (choice==1) then ! <g|Vnl|c>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==2) then ! derivative w.r.t. atm. pos
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(2,ilmn),-gxfac_(1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir)*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

#ifdef MR_DEV
!        ------
         else if (choice==22) then ! mixed derivative w.r.t. atm. pos and q vector (at q=0)
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+qdir
!$OMP DO
           do ipw=1,npw
             if (idir==qdir) then
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end if
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+kpg(ipw,idir)*ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               ztab(ipw)=ztab(ipw)-ffnl(ipw,ffnl_dir1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
             ztab(ipw)=ztab(ipw)*two_pi
           end do
!$OMP END DO
#endif

#ifdef MR_DEV
!        ------
         else if (choice==25) then ! mixed derivative w.r.t. atm. pos and thwo q vectors (at q=0)
           !Use same notation as the notes for clarity
           ialpha=nalpha(idir)
           idelta=nbeta(idir)
           igamma=qdir
           idelgam=gamma(idelta,igamma)
!$OMP DO
           do ipw=1,npw
             if (ialpha==igamma) then
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1+idelta,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end if
             if (ialpha==idelta) then
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1+igamma,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end if
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+kpg(ipw,ialpha)*ffnl(ipw,4+idelgam,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               ztab(ipw)=ztab(ipw)-ffnl(ipw,4+idelgam,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
             ztab(ipw)=ztab(ipw)*two_pi
           end do
!$OMP END DO
#endif
!        ------
         else if (choice==3) then ! derivative w.r.t. strain
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           if (idir<=3) then
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn) &
&                 *cmplx(dgxdtfac_(1,1,ilmn)-gxfac_(1,ilmn),dgxdtfac_(2,1,ilmn)-gxfac_(2,ilmn),kind=dp) &
&                 -ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           else
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) &
&                 -ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           end if

#ifdef MR_DEV
!        ------
         else if (choice==33) then ! mixed derivative w.r.t. strain and q vector (at q=0)
           !Use same notation as the notes for clarity
           ibeta=nalpha(idir)
           idelta=nbeta(idir)
           igamma=qdir
           idelgam=gamma(idelta,igamma)
!$OMP DO
           do ipw=1,npw
             if (ibeta==igamma) then
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+onehalf*ffnl(ipw,1+idelta,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
                 ztab(ipw)=ztab(ipw)+half*ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp)
               end do
             end if
             if (ibeta==idelta) then
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+onehalf*ffnl(ipw,1+igamma,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
                 ztab(ipw)=ztab(ipw)+half*ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
               end do
             end if
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+kpg(ipw,ibeta)*ffnl(ipw,4+idelgam,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1+idelta,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1+igamma,ilmn)*cmplx(d2gxdtfac_(1,2,ilmn),d2gxdtfac_(2,2,ilmn),kind=dp)
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(d2gxdtfac_(1,3,ilmn),d2gxdtfac_(2,3,ilmn),kind=dp)
             end do
             ztab(ipw)=ztab(ipw)*two_piinv
           end do
!$OMP END DO
#endif

!        ------
         else if (choice==5) then ! full derivative w.r.t. k
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp) &
&               +ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==51) then ! right derivative: <G|p>V<dp/dk|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==52) then ! left derivative: <G|dp/dk>V<p|psi>
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfac_(1,ilmn),gxfac_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>V<dp/dk_(idir-1)|psi>
!                                                     -<G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,fdf,ilmn)*cmplx(dgxdtfac_(1,2,ilmn),dgxdtfac_(2,2,ilmn),kind=dp) &
&               -ffnl(ipw,fdb,ilmn)*cmplx(dgxdtfac_(1,1,ilmn),dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==54) then ! mixed derivative w.r.t. atm. pos and k
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfac_(2,1,ilmn),-dgxdtfac_(1,1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir1(idir))*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==8) then ! full second order derivative w.r.t. k
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=idir1(idir); ffnl_dir(2)=idir2(idir)
           ffnl_dir(3) = gamma(ffnl_dir(1),ffnl_dir(2))
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,4+ffnl_dir(3),ilmn)*cmplx(    gxfac_(1,  ilmn),    gxfac_(2,  ilmn),kind=dp)&
&               +ffnl(ipw,1            ,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)&
&               +ffnl(ipw,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfac_(1,2,ilmn), dgxdtfac_(2,2,ilmn),kind=dp)&
&               +ffnl(ipw,1+ffnl_dir(2),ilmn)*cmplx( dgxdtfac_(1,1,ilmn), dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==81) then
           ! partial second order derivative w.r.t. k
           ! full derivative w.r.t. k1, right derivative w.r.t. k2
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=1; if(dimffnl>2) ffnl_dir(1)=idir1(idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,1            ,ilmn)*cmplx(d2gxdtfac_(1,1,ilmn),d2gxdtfac_(2,1,ilmn),kind=dp)&
&               +ffnl(ipw,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfac_(1,1,ilmn), dgxdtfac_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else
!$OMP WORKSHARE
           ztab(:)=czero
!$OMP END WORKSHARE
         end if


!        ------
!$OMP DO
         do ipw=1,npw
           ztab(ipw)=ztab(ipw)*cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           vect(1,ipw+ipwshft)=vect(1,ipw+ipwshft)+real(ztab(ipw))
           vect(2,ipw+ipwshft)=vect(2,ipw+ipwshft)+aimag(ztab(ipw))
         end do
!$OMP END DO

!$OMP END PARALLEL
       end if

!      Compute <g|S|c> (or derivatives) for each plane wave:
       if (paw_opt>=3) then
!$OMP PARALLEL PRIVATE(ilmn,ipw,fdf,fdb,ffnl_dir1)

!        ------
         if (choice==1) then ! <g|S|c>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==2) then ! derivative w.r.t. atm. pos
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(gxfacs_(2,ilmn),-gxfacs_(1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir)*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==3) then ! derivative w.r.t. strain
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
           if (idir<=3) then
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn) &
&                 *cmplx(dgxdtfacs_(1,1,ilmn)-gxfacs_(1,ilmn),dgxdtfacs_(2,1,ilmn)-gxfacs_(2,ilmn),kind=dp)&
&                 -ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           else
!$OMP DO
             do ipw=1,npw
               ztab(ipw)=czero
               do ilmn=1,nlmn
                 ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) &
&                 -ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
               end do
             end do
!$OMP END DO
           end if

!        ------
         else if (choice==5) then ! full derivative w.r.t. k
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp) &
&               +ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==51) then ! right derivative: <G|p>S<dp/dk|psi>
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==52) then ! left derivative: <G|dp/dk>S<p|psi>
           ffnl_dir1=2; if(dimffnl>2) ffnl_dir1=1+idir
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,ffnl_dir1,ilmn)*cmplx(gxfacs_(1,ilmn),gxfacs_(2,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==53) then ! twist derivative: <G|dp/dk_(idir+1)>S<dp/dk_(idir-1)|psi> -
!          <G|dp/dk_(idir-1)>V<dp/dk_(idir+1)|psi>
           fdf = ffnl_dir_dat(2*idir-1)
           fdb = ffnl_dir_dat(2*idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,fdf,ilmn)*cmplx(dgxdtfacs_(1,2,ilmn),dgxdtfacs_(2,2,ilmn),kind=dp) &
&               -ffnl(ipw,fdb,ilmn)*cmplx(dgxdtfacs_(1,1,ilmn),dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==54) then ! mixed derivative w.r.t. atm. pos and k
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(dgxdtfacs_(2,1,ilmn),-dgxdtfacs_(1,1,ilmn),kind=dp)
             end do
             ztab(ipw)=two_pi*kpg(ipw,idir1(idir))*ztab(ipw)
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw)+ffnl(ipw,1,ilmn)*cmplx(d2gxdtfacs_(1,1,ilmn),d2gxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==8) then ! full second order derivative w.r.t. k
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=idir1(idir); ffnl_dir(2)=idir2(idir)
           ffnl_dir(3) = gamma(ffnl_dir(1),ffnl_dir(2))
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,4+ffnl_dir(3),ilmn)*cmplx(    gxfacs_(1,  ilmn),    gxfacs_(2,  ilmn),kind=dp)&
&               +ffnl(ipw,1            ,ilmn)*cmplx(d2gxdtfacs_(1,1,ilmn),d2gxdtfacs_(2,1,ilmn),kind=dp)&
&               +ffnl(ipw,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfacs_(1,2,ilmn), dgxdtfacs_(2,2,ilmn),kind=dp)&
&               +ffnl(ipw,1+ffnl_dir(2),ilmn)*cmplx( dgxdtfacs_(1,1,ilmn), dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else if (choice==81) then
           ! partial second order derivative w.r.t. k
           ! full derivative w.r.t. k1, right derivative w.r.t. k2
           !idir= (xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9)
           ffnl_dir(1)=1; if(dimffnl>2) ffnl_dir(1)=idir1(idir)
!$OMP DO
           do ipw=1,npw
             ztab(ipw)=czero
             do ilmn=1,nlmn
               ztab(ipw)=ztab(ipw) &
&               +ffnl(ipw,1            ,ilmn)*cmplx(d2gxdtfacs_(1,1,ilmn),d2gxdtfacs_(2,1,ilmn),kind=dp)&
&               +ffnl(ipw,1+ffnl_dir(1),ilmn)*cmplx( dgxdtfacs_(1,1,ilmn), dgxdtfacs_(2,1,ilmn),kind=dp)
             end do
           end do
!$OMP END DO

!        ------
         else
!$OMP WORKSHARE
           ztab(:)=czero
!$OMP END WORKSHARE
         end if


!        ------
!        The OMP WORKSHARE directive doesn't have a good performance with Intel Compiler
!        !$OMP WORKSHARE
!        ztab(:)=ztab(:)*cmplx(ph3d(1,:,iaph3d),-ph3d(2,:,iaph3d),kind=dp)
!        !$OMP END WORKSHARE
!        !$OMP WORKSHARE
!        svect(1,1+ipwshft:npw+ipwshft)=svect(1,1+ipwshft:npw+ipwshft)+real(ztab(:))
!        svect(2,1+ipwshft:npw+ipwshft)=svect(2,1+ipwshft:npw+ipwshft)+aimag(ztab(:))
!        !$OMP END WORKSHARE
!$OMP DO
         do ipw=1,npw
           ztab(ipw)=ztab(ipw)*cmplx(ph3d(1,ipw,iaph3d),-ph3d(2,ipw,iaph3d),kind=dp)
           svect(1,ipw+ipwshft)=svect(1,ipw+ipwshft)+real(ztab(ipw))
           svect(2,ipw+ipwshft)=svect(2,ipw+ipwshft)+aimag(ztab(ipw))
         end do
!$OMP END DO
!$OMP END PARALLEL
       end if

!      End loop on atoms
     end do
!    End loop on spinors
   end do

!  ==========================================================================
 end if

 ABI_DEALLOCATE(ztab)

 if (paw_opt/=3) then
   ABI_DEALLOCATE(gxfac_)
   if (choice>1) then
     ABI_DEALLOCATE(dgxdtfac_)
   end if
   if (choice==54.or.choice==8.or.choice==81) then
     ABI_DEALLOCATE(d2gxdtfac_)
   end if
 end if
 if (paw_opt>=3) then
   ABI_DEALLOCATE(gxfacs_)
   if (choice>1) then
     ABI_DEALLOCATE(dgxdtfacs_)
   end if
   if (choice==54.or.choice==8.or.choice==81) then
     ABI_DEALLOCATE(d2gxdtfacs_)
   end if
 end if

 DBG_EXIT("COLL")

#if !defined HAVE_OPENMP
!Fake use of unused variable
 if (.false.) write(std_out,*) ipw
#endif

end subroutine opernlb_ylm
!!***

end module m_opernlb_ylm
!!***

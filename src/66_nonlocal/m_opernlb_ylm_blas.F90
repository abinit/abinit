!!****m* ABINIT/m_opernlb_ylm_blas
!! NAME
!!  m_opernlb_ylm_blas
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (MT)
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

module m_opernlb_ylm_blas

 use defs_basis
 use m_abicore
 use m_errors
#if defined HAVE_OPENMP
 use OMP_LIB
#endif

 implicit none

 private
!!***

 public :: opernlb_ylm_blas
!!***

contains
!!***

!!****f* ABINIT/opernlb_ylm_blas
!! NAME
!! opernlb_ylm_blas
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
!!      nonlop_ylm_blas
!!
!! CHILDREN
!!
!! SOURCE

subroutine opernlb_ylm_blas(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
&                      d2gxdtfac,d2gxdtfac_sij,dgxdtfac,dgxdtfac_sij,dimffnl,ffnl,gxfac,gxfac_sij,&
&                      ia3,idir,indlmn,kpg,matblk,ndgxdtfac,nd2gxdtfac,nincat,nkpg,nlmn,nloalg,npw,&
&                      nspinor,paw_opt,ph3d,svect,ucvol,vect,qdir)

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
 integer :: ii,il,ilmn,ipw,jpw,ipwshft,ispinor,jc,nthreads,ffnl_dir1,ffnl_dir(3)
 real(dp) :: scale,two_piinv,wt
 logical :: parity
!arrays
 integer,parameter :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,parameter :: idir1(9)=(/1,1,1,2,2,2,3,3,3/),idir2(9)=(/1,2,3,1,2,3,1,2,3/)
 integer,parameter :: nalpha(9)=(/1,2,3,3,3,2,2,1,1/),nbeta(9)=(/1,2,3,2,1,1,3,3,2/)
 real(dp),allocatable :: d2gxdtfac_(:,:,:),d2gxdtfacs_(:,:,:),dgxdtfac_(:,:,:),dgxdtfacs_(:,:,:),gxfac_(:,:),gxfacs_(:,:)
 real(dp),allocatable :: scalr(:),scali(:)
! complex(dpc),allocatable :: ztab(:)
 complex(dp) :: ctmp, cil(4)

! *************************************************************************

 DBG_ENTER("COLL")

!Nothing to do when choice=4, 6 or 23
 if (choice==4.or.choice==6.or.choice==23) return

 if (abs(choice)>1) then
   MSG_ERROR('Only abs(choice)<=0 is available for now.')
 end if
 if (cplex/=2) then
   MSG_ERROR('Only cplex=2 is available for now.')
 end if
 if (cplex_fac/=2) then
   MSG_ERROR('Only cplex_fac=2 is available for now.')
 end if
! if (istwf_k/=1) then
!   MSG_ERROR('Only istwf_k=1 is available for now.')
! end if
 if (paw_opt/=4) then
   MSG_ERROR('Only paw_opt=4 is available for now.')
 end if
 !
!DDK not compatible with istwkf > 1
 if(cplex==1.and.(any(cplex_dgxdt(:)==2).or.any(cplex_d2gxdt(:)==2)))then
   MSG_BUG("opernlb_ylm_blas+ddk not compatible with istwfk>1")
 end if

!Inits
 wt=four_pi/sqrt(ucvol)
 nthreads=1
#if defined HAVE_OPENMP
 nthreads=OMP_GET_NUM_THREADS()
#endif
 if (nthreads>1) then
   MSG_ERROR('Only nthreads=1 is available for now.')
 end if

 if (paw_opt/=3) then
   ABI_ALLOCATE(gxfac_,(nlmn,2))
!   gxfac_(:,:)=zero
!   if (choice>1) then
!     ABI_ALLOCATE(dgxdtfac_,(2,ndgxdtfac,nlmn))
!     if(ndgxdtfac>0) dgxdtfac_(:,:,:)=zero
!   end if
!   if (choice==54.or.choice==8.or.choice==81.or.choice==33) then
!     ABI_ALLOCATE(d2gxdtfac_,(2,nd2gxdtfac,nlmn))
!     if(nd2gxdtfac>0) d2gxdtfac_(:,:,:)=zero
!   end if
 end if
 if (paw_opt>=3) then
   ABI_ALLOCATE(gxfacs_,(nlmn,2))
!   gxfacs_(:,:)=zero
!   if (choice>1) then
!     ABI_ALLOCATE(dgxdtfacs_,(2,ndgxdtfac,nlmn))
!     if (ndgxdtfac>0) dgxdtfacs_(:,:,:)=zero
!   end if
!   if (choice==54.or.choice==8.or.choice==81) then
!     ABI_ALLOCATE(d2gxdtfacs_,(2,nd2gxdtfac,nlmn))
!     if (nd2gxdtfac>0) d2gxdtfacs_(:,:,:)=zero
!   end if
 end if

if (choice==33) two_piinv=1.0_dp/two_pi

! ABI_ALLOCATE(ztab,(npw))
 ABI_ALLOCATE(scalr,(npw))
 ABI_ALLOCATE(scali,(npw))

! (-i)^l
 cil(1) = ( 1.0_DP, 0.0_DP) * wt
 cil(2) = ( 0.0_DP,-1.0_DP) * wt
 cil(3) = (-1.0_DP, 0.0_DP) * wt
 cil(4) = ( 0.0_DP, 1.0_DP) * wt

!==========================================================================
!========== STANDARD VERSION ==============================================
!==========================================================================

!Loop on spinorial components
 do ispinor=1,nspinor
   ipwshft=(ispinor-1)*npw

!  Loop on atoms (blocking)
   do ia=1,nincat
     iaph3d=ia;if (nloalg(2)>0) iaph3d=ia+ia3-1
!    Scale gxfac with 4pi/sqr(omega).(-i)^l
     if (paw_opt/=3) then
       do ilmn=1,nlmn
         il=mod(indlmn(1,ilmn),4)+1
         ctmp = cil(il) * cmplx( gxfac(1,ilmn,ia,ispinor), gxfac(2,ilmn,ia,ispinor), kind=DP )
         gxfac_(ilmn,1) =  real(ctmp)
         gxfac_(ilmn,2) = aimag(ctmp)
       end do
     end if

!    Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
     if (paw_opt>=3) then
       do ilmn=1,nlmn
         il=mod(indlmn(1,ilmn),4)+1
         ctmp = cil(il) * cmplx( gxfac_sij(1,ilmn,ia,ispinor), gxfac_sij(2,ilmn,ia,ispinor), kind=DP )
         gxfacs_(ilmn,1) =  real(ctmp)
         gxfacs_(ilmn,2) = aimag(ctmp)
       end do
     end if

!    Compute <g|Vnl|c> (or derivatives) for each plane wave:

     if (paw_opt/=3) then

       call DGEMV('N',npw,nlmn,1.0_DP,ffnl(:,1,:),npw,gxfac_(:,1),1,0.0_DP,scalr,1)
       call DGEMV('N',npw,nlmn,1.0_DP,ffnl(:,1,:),npw,gxfac_(:,2),1,0.0_DP,scali,1)

       do ipw=1,npw
         jpw=ipw+ipwshft
         vect(1,jpw)=vect(1,jpw)+scalr(ipw)*ph3d(1,ipw,iaph3d)+scali(ipw)*ph3d(2,ipw,iaph3d)
         vect(2,jpw)=vect(2,jpw)-scalr(ipw)*ph3d(2,ipw,iaph3d)+scali(ipw)*ph3d(1,ipw,iaph3d)
       end do

     end if

!    Compute <g|S|c> (or derivatives) for each plane wave:

     if (paw_opt>=3) then

       call DGEMV('N',npw,nlmn,1.0_DP,ffnl(:,1,:),npw,gxfacs_(:,1),1,0.0_DP,scalr,1)
       call DGEMV('N',npw,nlmn,1.0_DP,ffnl(:,1,:),npw,gxfacs_(:,2),1,0.0_DP,scali,1)

       do ipw=1,npw
         jpw=ipw+ipwshft
         svect(1,jpw)=svect(1,jpw)+scalr(ipw)*ph3d(1,ipw,iaph3d)+scali(ipw)*ph3d(2,ipw,iaph3d)
         svect(2,jpw)=svect(2,jpw)-scalr(ipw)*ph3d(2,ipw,iaph3d)+scali(ipw)*ph3d(1,ipw,iaph3d)
       end do

     end if

!    End loop on atoms
   end do
 end do !  End loop on spinors

! ABI_DEALLOCATE(ztab)
 ABI_DEALLOCATE(scalr)
 ABI_DEALLOCATE(scali)

 if (paw_opt/=3) then
   ABI_DEALLOCATE(gxfac_)
   if (choice>1) then
     ABI_DEALLOCATE(dgxdtfac_)
   end if
   if (choice==54.or.choice==8.or.choice==81.or.choice==33) then
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

end subroutine opernlb_ylm_blas
!!***

end module m_opernlb_ylm_blas
!!***

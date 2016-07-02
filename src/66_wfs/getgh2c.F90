!{\src2tex{textfont=tt}}
!!****f* ABINIT/getgh2c
!!
!! NAME
!! getgh2c
!!
!! FUNCTION
!! Compute <G|H^(2)|C> (or <G|H^(2)-Eps.S^(2)|C>) for input vector |C> expressed in reciprocal space.
!! (H^(2) is the 2nd-order pertubed Hamiltonian, S^(2) is the 2nd-order perturbed overlap operator).
!! Result is put in array gh2c.
!! If required, part of <G|K(2)+Vnonlocal^(2)|C> not depending on VHxc^(2) is also returned in gvnl2.
!! If required, <G|S^(2)|C> is returned in gs2c (S=overlap - PAW only)
!!
!! COPYRIGHT
!! Copyright (C) 2015-2016 ABINIT group (MT,JLJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cwavef(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C>
!!  ddkinpw(npw)=derivative of the (modified) kinetic energy for each plane wave at k (Hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  lambda=real use to apply H^(2)-lambda.S^(2)
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  npw,npw1=number of planewaves in basis sphere (should be identical as q=0)
!!  optlocal=0: local part of H^(2) is not computed in gh2c=<G|H^(2)|C>
!!           1: local part of H^(2) is computed in gh2c=<G|H^(2)|C>
!!  optnl=0: non-local part of H^(2) is not computed in gh2c=<G|H^(2)|C>
!!        1: non-local part of H^(2) depending on VHxc^(2) is not computed in gh2c=<G|H^(2)|C>
!!        2: non-local part of H^(2) is totally computed in gh2c=<G|H^(2)|C>
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 2nd-order Hamiltonian at k,k+q
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H^(2)|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S^(2)|C> have to be computed in gs2c in addition to gh2c
!!                       if -1, matrix elements <G|H^(2)-lambda.S^(2)|C> have to be computed in gh2c (gs2c not used)
!!  tim_getgh2c=timing code of the calling subroutine (can be set to 0 if not attributed)
!!  usevnl=1 if gvnl2=(part of <G|K^(2)+Vnl^(2)-lambda.S^(2)|C> not depending on VHxc^(2)) has to be input/output
!!
!! OUTPUT
!! gh2c(2,npw*nspinor)= <G|H^(2)|C> or  <G|H^(2)-lambda.S^(2)|C>
!!                     (only kinetic+non-local parts if optlocal=0)
!! if (usevnl==1)
!!  gvnl2(2,npw1*nspinor)=  part of <G|K^(2)+Vnl^(2)|C> not depending on VHxc^(2)              (sij_opt/=-1)
!!                       or part of <G|K^(2)+Vnl^(2)-lambda.S^(2)|C> not depending on VHxc^(2) (sij_opt==-1)
!! if (sij_opt=1)
!!  gs2c(2,npw*nspinor)=<G|S^(2)|C> (S=overlap).
!!
!! PARENTS
!!      m_rf2
!!
!! CHILDREN
!!      nonlop
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getgh2c(cwavef,cwaveprj,gh2c,gs2c,gs_hamkq,gvnl2,idir,ipert,lambda,&
&                  mpi_enreg,optlocal,optnl,rf_hamkq,sij_opt,tim_getgh2c,usevnl)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawcprj,     only : pawcprj_type
 use m_hamiltonian, only : gs_hamiltonian_type,rf_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getgh2c'
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,optlocal,optnl,sij_opt,tim_getgh2c,usevnl
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq
!arrays
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(inout),target :: gvnl2(:,:)
 real(dp),intent(out) :: gh2c(:,:),gs2c(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,idir1,idir2,idirc,ipw,ipws,ispinor
 integer :: my_nspinor,natom,nnlout=1,npw,npw1,paw_opt,signs,tim_nonlop
 logical :: has_kin,has_vnl
 real(dp) :: enlout_dum(1)
 character(len=500) :: msg
!arrays
! integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: alpha(9)=(/1,2,3,2,1,1,3,3,2/),beta(9)=(/1,2,3,3,3,2,2,1,1/)
!real(dp) :: tsec(2)
 real(dp) :: svectout_dum(1,1)
 real(dp),ABI_CONTIGUOUS pointer :: gvnl2_(:,:)
 real(dp), pointer :: ddkinpw(:),kinpw1(:)

! *********************************************************************

!Keep track of total time spent in getgh2c
!call timab(196+tim_getgh2c,1,tsec)

!======================================================================
!== Initialisations and compatibility tests
!======================================================================

 npw  =gs_hamkq%npw_k
 npw1 =gs_hamkq%npw_kp
 natom=gs_hamkq%natom

!Compatibility tests
 if(ipert/=natom+10.and.ipert/=natom+11)then
   msg='only ipert<>natom+10/natom+11 implemented!'
   MSG_BUG(msg)
 end if
 if(gs_hamkq%usepaw==0.and.gs_hamkq%usecprj==1)then
   msg='usecprj==1 not allowed for NC psps!'
   MSG_BUG(msg)
 end if
 if (mpi_enreg%paral_spinor==1) then
   msg='Not compatible with parallelization over spinorial components!'
   MSG_BUG(msg)
 end if
 if (gs_hamkq%nvloc>1) then
   msg='Not compatible with nvloc=4 (non-coll. magnetism)!'
   MSG_BUG(msg)
 end if

!Check sizes
 my_nspinor=max(1,gs_hamkq%nspinor/mpi_enreg%nproc_spinor)
 if (size(cwavef)<2*npw*my_nspinor) then
   msg='wrong size for cwavef!'
   MSG_BUG(msg)
 end if
 if (size(gh2c)<2*npw1*my_nspinor) then
   msg='wrong size for gh2c!'
   MSG_BUG(msg)
 end if
 if (usevnl/=0) then
   if (size(gvnl2)<2*npw1*my_nspinor) then
     msg='wrong size for gvnl2!'
     MSG_BUG(msg)
   end if
 end if
 if (sij_opt==1) then
   if (size(gs2c)<2*npw1*my_nspinor) then
     msg='wrong size for gs2c!'
     MSG_BUG(msg)
   end if
 end if
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj/=0) then
   if (size(cwaveprj)<gs_hamkq%natom*my_nspinor) then
     msg='wrong size for cwaveprj!'
     MSG_BUG(msg)
   end if
 end if

 tim_nonlop=8
 if (tim_getgh2c==1.and.ipert<=natom) tim_nonlop=7
 if (tim_getgh2c==2.and.ipert<=natom) tim_nonlop=5
 if (tim_getgh2c==1.and.ipert> natom) tim_nonlop=8
 if (tim_getgh2c==2.and.ipert> natom) tim_nonlop=5
 if (tim_getgh2c==3                 ) tim_nonlop=0

! if (ipert==natom+10) then
 idir1=alpha(idir);idir2=beta(idir)
! else if (ipert==natom+11) then
!   idir1=alpha(3+idir);idir2=beta(3+idir)
! end if

!======================================================================
!== Apply the 2nd-order local potential to the wavefunction
!======================================================================

!Not implemented

 if (ipert/=natom+10.and.ipert/=natom+11.and.optlocal>0) then
   msg='local part not implemented'
   MSG_BUG(msg)

 else
!  In the case of ddk operator, no local contribution (also because no self-consistency)
!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gh2c(:,ipw)=zero
   end do

 end if

!======================================================================
!== Apply the 2st-order non-local potential to the wavefunction
!======================================================================

 has_vnl=(ipert==natom+10.or.ipert==natom+11)
 
!Use of gvnl2 depends on usevnl
 if (usevnl==1) then
   gvnl2_ => gvnl2
 else
   ABI_ALLOCATE(gvnl2_,(2,npw1*my_nspinor))
 end if

!k-point perturbation
!  -------------------------------------------
 if (has_vnl.and.(optnl>0.or.sij_opt/=0)) then

   idirc=3*(idir1-1)+idir2 !xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9, (xyz,xyz)=(idir1,idir2)

   if (gs_hamkq%usepaw==1) then
!    cpopt=-1+5*gs_hamkq%usecprj -> cpopt=-1
!    usecprj should be used again
     choice=8; signs=2; cpopt=-1; paw_opt=1; if (sij_opt/=0) paw_opt=sij_opt+3
     call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_hamkq,idirc,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,gs2c,tim_nonlop,cwavef,gvnl2_)
   else
     choice=8; signs=2; cpopt=-1 ; paw_opt=0
     call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_hamkq,idirc,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,svectout_dum,tim_nonlop,cwavef,gvnl2_)
   end if

!No non-local part
!-------------------------------------------
 else

   if (optnl>=1) then
 !$OMP PARALLEL DO
     do ipw=1,npw1*my_nspinor
       gvnl2_(:,ipw)=zero
     end do
   end if
   if (sij_opt/=0) gs2c=zero

 end if

!======================================================================
!== Apply the 2nd-order kinetic operator to the wavefunction
!======================================================================

 has_kin=(ipert==natom+10.or.ipert==natom+11)

!k-point perturbation
!-------------------------------------------
 if (associated(gs_hamkq%kinpw_kp)) then
   kinpw1 => gs_hamkq%kinpw_kp
 else if (optnl>=1.or.has_kin) then
   msg='need kinpw1 allocated!'
   MSG_BUG(msg)
 end if
 if (associated(rf_hamkq%ddkinpw_k)) then
   ddkinpw => rf_hamkq%ddkinpw_k
 else if (has_kin) then
   msg='need ddkinpw allocated!'
   MSG_BUG(msg)
 end if

 if (has_kin) then
   do ispinor=1,my_nspinor
 !$OMP PARALLEL DO PRIVATE(ipw,ipws) SHARED(cwavef,ispinor,gvnl2_,ddkinpw,kinpw1,npw,my_nspinor)
     do ipw=1,npw
       ipws=ipw+npw*(ispinor-1)
       if(kinpw1(ipw)<huge(zero)*1.d-11)then
         gvnl2_(1,ipws)=gvnl2_(1,ipws)+ddkinpw(ipw)*cwavef(1,ipws)
         gvnl2_(2,ipws)=gvnl2_(2,ipws)+ddkinpw(ipw)*cwavef(2,ipws)
       else
         gvnl2_(1,ipws)=zero
         gvnl2_(2,ipws)=zero
       end if
     end do

   end do
 end if

!======================================================================
!== Sum contributions to get the application of H^(2) to the wf
!======================================================================
!Also filter the wavefunctions for large modified kinetic energy

!Add non-local+kinetic to local part
 if (optnl>=1.or.has_kin) then
   do ispinor=1,my_nspinor
     ipws=(ispinor-1)*npw1
 !$OMP PARALLEL DO PRIVATE(ipw) SHARED(gh2c,gvnl2_,kinpw1,ipws,npw1)
     do ipw=1+ipws,npw1+ipws
       if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
         gh2c(1,ipw)=gh2c(1,ipw)+gvnl2_(1,ipw)
         gh2c(2,ipw)=gh2c(2,ipw)+gvnl2_(2,ipw)
       else
         gh2c(1,ipw)=zero
         gh2c(2,ipw)=zero
       end if
     end do
   end do
 end if

 if (usevnl==1) then
   nullify(gvnl2_)
 else
   ABI_DEALLOCATE(gvnl2_)
 end if

!call timab(196+tim_getgh2c,2,tsec)

end subroutine getgh2c
!!***

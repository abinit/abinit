!!****m* ABINIT/m_opernld_ylm
!! NAME
!!  m_opernld_ylm
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MT)
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

module m_opernld_ylm

 use defs_basis
 use m_errors
 use m_abicore

 implicit none

 private
!!***

 public :: opernld_ylm
!!***

contains
!!***

!!****f* ABINIT/opernld_ylm
!! NAME
!! opernld_ylm
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   in order to get contributions to energy/forces/stress/dyn.matrix/elst tens.
!!   from projected scalars
!! * Operate with the non-local projectors and the overlap matrix Sij
!!   in order to get contributions to <c|S|c>
!!   from projected scalars
!!
!! INPUTS
!!  choice=chooses possible output
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)=2nd gradients of projected scalars
!!  dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)=gradients of projected scalars
!!  dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)=gradients of reduced projected scalars
!!                                                    related to Vnl (NL operator)
!!  dgxdtfac_sij(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)=gradients of reduced projected scalars
!!                                                        related to Sij (overlap)
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars
!!  gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  gxfac_sij(cplex,nlmn,nincat,nspinor)= reduced projected scalars related to Sij (overlap)
!!  ia3=gives the absolute number of the first atom in the subset presently treated
!!  natom=number of atoms in cell
!!  nd2gxdt=second dimension of d2gxdt
!!  ndgxdt=second dimension of dgxdt
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nnlout=dimension of enlout
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)             -> the energy
!!      if choice=2 : enlout(3*natom)       -> 1st deriv. of energy wrt atm. pos (forces)
!!      if choice=3 : enlout(6)             -> 1st deriv. of energy wrt strain (stresses)
!!      if choice=4 : enlout(6*natom)       -> 2nd deriv. of energy wrt 2 atm. pos (dyn. mat.)
!!      if choice=23: enlout(6+3*natom)     -> 1st deriv. of energy wrt atm. pos (forces) and
!!                                             1st deriv. of energy wrt strain (stresses)
!!      if choice=24: enlout(9*natom)       -> 1st deriv. of energy wrt atm. pos (forces) and
!!                                             2nd deriv. of energy wrt 2 atm. pos (dyn. mat.)
!!      if choice=5 : enlout(3)             -> 1st deriv. of energy wrt k
!!      if choice=53: enlout(3)             -> 1st deriv. (twist) of energy wrt k
!!      if choice=54: enlout(18*natom)      -> 2nd deriv. of energy wrt atm. pos and right k (Born eff. charge)
!!      if choice=55: enlout(36)            -> 2nd deriv. of energy wrt strain and right k (piezoelastic tensor)
!!      if choice=6 : enlout(36+18*natom)   -> 2nd deriv. of energy wrt 2 strains (elast. tensor) and
!!                                             2nd deriv. of energy wrt to atm. pos and strain (internal strain)
!!      if choice=8 : enlout(6)             -> 2nd deriv. of energy wrt 2 k
!!      if choice=81: enlout(18)            -> 2nd deriv. of energy wrt k and right k
!! --If (paw_opt==3)
!!      if choice=1 : enlout(1)             -> contribution to <c|S|c> (note: not including <c|c>)
!!      if choice=2 : enlout(3*natom)       -> contribution to <c|dS/d_atm.pos|c>
!!      if choice=54: enlout(18*natom)      -> 2nd deriv. of energy wrt atm. pos and right k (Born eff. charge)
!!      if choice=55: enlout(36)            -> 2nd deriv. of energy wrt strain and right k (piezoelastic tensor)
!!      if choice=8 : enlout(6)             -> 2nd deriv. of energy wrt 2 k
!!      if choice=81: enlout(18)            -> 2nd deriv. of energy wrt k and right k
!! --If (paw_opt==4)
!!      not available
!!
!! NOTES
!! Operate for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE

subroutine opernld_ylm(choice,cplex,cplex_fac,ddkk,dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,&
&                      enlk,enlout,fnlk,gx,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,&
&                      ndgxdtfac,nincat,nlmn,nnlout,nspinor,paw_opt,strnlk,&
&                      enlout_im)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,cplex_fac,ia3,natom,nd2gxdt,ndgxdt
 integer,intent(in) :: ndgxdtfac,nincat,nlmn,nnlout,nspinor,paw_opt
 real(dp),intent(inout) :: enlk
!arrays
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: gx(cplex,nlmn,nincat,nspinor),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(inout) :: ddkk(6),enlout(nnlout),fnlk(3*natom),strnlk(6)
 real(dp),intent(inout),optional :: enlout_im(nnlout)

!Local variables-------------------------------
!scalars
 integer :: ia,iashift,ilmn,iplex,ishift,ispinor,mu,mua,mua1,mua2,mub,mushift,mut,muu,nu,nushift
 real(dp) :: dummy
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,parameter :: twist_dir(6)=(/2,3,3,1,1,2/)
 real(dp) :: d2gx(cplex),enlj(6*cplex),gxfacj(cplex)
 real(dp),allocatable :: enljj(:)
 complex(dpc),allocatable :: cft(:,:), cfu(:,:)

! *************************************************************************

 ABI_CHECK(cplex_fac>=cplex,'BUG: invalid cplex_fac<cplex!')

 if (paw_opt==0.or.paw_opt==1.or.paw_opt==2) then

!  ============== Accumulate the non-local energy ===============
   if (choice==1) then
     if (present(enlout_im).and.cplex==2) then ! cplex=cplex_fac=2
       do ispinor=1,nspinor
         do ia=1,nincat
           do ilmn=1,nlmn
             enlout   (1)=enlout   (1)+gxfac(1,ilmn,ia,ispinor)*gx(1,ilmn,ia,ispinor)
             enlout   (1)=enlout   (1)+gxfac(2,ilmn,ia,ispinor)*gx(2,ilmn,ia,ispinor)
             enlout_im(1)=enlout_im(1)+gxfac(2,ilmn,ia,ispinor)*gx(1,ilmn,ia,ispinor)
             enlout_im(1)=enlout_im(1)-gxfac(1,ilmn,ia,ispinor)*gx(2,ilmn,ia,ispinor)
           end do
         end do
       end do
     else if (present(enlout_im).and.cplex_fac==2) then ! cplex=1,cplex_fac=2
       do ispinor=1,nspinor
         do ia=1,nincat
           do ilmn=1,nlmn
             enlout   (1)=enlout   (1)+gxfac(1,ilmn,ia,ispinor)*gx(1,ilmn,ia,ispinor)
             enlout_im(1)=enlout_im(1)+gxfac(2,ilmn,ia,ispinor)*gx(1,ilmn,ia,ispinor)
           end do
         end do
       end do
     else ! only the real part is needed or the imaginary part is zero
       do ispinor=1,nspinor
         do ia=1,nincat
           do ilmn=1,nlmn
             do iplex=1,cplex
               enlout(1)=enlout(1)+gxfac(iplex,ilmn,ia,ispinor)*gx(iplex,ilmn,ia,ispinor)
             end do
           end do
         end do
       end do
     end if
   end if

!  ============ Accumulate the forces contributions =============
   if (choice==2.or.choice==23.or.choice==24) then
     ishift=0;if (choice==23) ishift=6
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:3)=zero
         iashift=3*(ia+ia3-2)+ishift
         do ilmn=1,nlmn
           do mu=1,3
             dummy = zero  ! Dummy needed here to get the correct forces with intel -O3
             do iplex=1,cplex
               !enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu+ishift,ilmn,ia,ispinor)
               dummy=dummy+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu+ishift,ilmn,ia,ispinor)
             end do
             enlj(mu)=enlj(mu)+dummy
           end do
         end do
         enlout(iashift+1:iashift+3)=enlout(iashift+1:iashift+3)+two*enlj(1:3)
       end do
     end do
   end if

!  ======== Accumulate the stress tensor contributions ==========
   if (choice==3.or.choice==23) then
     enlj(1:6)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           gxfacj(1:cplex)=gxfac(1:cplex,ilmn,ia,ispinor)
           do iplex=1,cplex
             enlk=enlk+gxfacj(iplex)*gx(iplex,ilmn,ia,ispinor)
           end do
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfacj(iplex)*dgxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
         end do
       end do
     end do
     enlout(1:6)=enlout(1:6)+two*enlj(1:6)
   end if

!  ====== Accumulate the dynamical matrix contributions =========
   if (choice==4.or.choice==24) then
     ishift=0;if (choice==24) ishift=3*natom
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:6)=zero
         iashift=6*(ia+ia3-2)+ishift
         do ilmn=1,nlmn
           do mu=1,6
             mua=alpha(mu);mub=beta(mu)
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)&
&               +dgxdtfac(iplex,mub,ilmn,ia,ispinor)*dgxdt(iplex,mua,ilmn,ia,ispinor)
             end do
           end do
         end do
         enlout(iashift+1:iashift+6)=enlout(iashift+1:iashift+6)+two*enlj(1:6)
       end do
     end do
   end if

!  ======== Accumulate the contributions of derivatives of E wrt to k ==========
   if (choice==5) then
     enlj(1:3)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         if(cplex==2)then
           do ilmn=1,nlmn
             do mu=1,3
               do iplex=1,cplex
                 enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu,ilmn,ia,ispinor)
               end do
             end do
           end do
!        If cplex=1, dgxdt is pure imaginary; thus there is no contribution
         else if (cplex_fac==2) then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(mu)=enlj(mu)+gxfac(2,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         end if
       end do
     end do
     enlout(1:3)=enlout(1:3)+two*enlj(1:3)
   end if

!  ======== Accumulate the contributions of partial derivatives of E wrt to k ==========
!  Choice 51: right derivative wrt to k ; Choice 52: left derivative wrt to k
   if (choice==51.or.choice==52) then
     enlj(1:6)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         if(cplex==2)then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(2*mu-1)=enlj(2*mu-1)+gxfac(1,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor) &
&               +gxfac(2,ilmn,ia,ispinor)*dgxdt(2,mu,ilmn,ia,ispinor)
               enlj(2*mu  )=enlj(2*mu  )+gxfac(1,ilmn,ia,ispinor)*dgxdt(2,mu,ilmn,ia,ispinor) &
&               -gxfac(2,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         else if (cplex_fac==2) then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(2*mu-1)=enlj(2*mu-1)+gxfac(2,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
               enlj(2*mu  )=enlj(2*mu  )+gxfac(1,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         else if (cplex_fac==1) then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(2*mu  )=enlj(2*mu  )+gxfac(1,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         end if
       end do
     end do
     if (choice==52) then
       enlj(2)=-enlj(2);enlj(4)=-enlj(4);enlj(6)=-enlj(6)
     end if
     enlout(1:6)=enlout(1:6)+enlj(1:6)
   end if

!  ======== Accumulate the contributions of twist derivatives of E wrt to k ==========
   if (choice==53) then
     enlj(1:3)=zero
     ABI_ALLOCATE(cft,(3,nlmn))
     ABI_ALLOCATE(cfu,(3,nlmn))
!    If cplex=1, dgxdt is pure imaginary;
!    If cplex_fac=1, dgxdtfac is pure imaginary;
     do ispinor=1,nspinor
       do ia=1,nincat
         if(cplex==2)then
           cft(1:3,1:nlmn)=cmplx(dgxdt(1,1:3,1:nlmn,ia,ispinor),dgxdt(2,1:3,1:nlmn,ia,ispinor))
         else
           cft(1:3,1:nlmn)=cmplx(zero,dgxdt(1,1:3,1:nlmn,ia,ispinor))
         end if
         if(cplex_fac==2)then
           cfu(1:3,1:nlmn)=cmplx(dgxdtfac(1,1:3,1:nlmn,ia,ispinor),dgxdtfac(2,1:3,1:nlmn,ia,ispinor))
         else
           cfu(1:3,1:nlmn)=cmplx(zero,dgxdtfac(1,1:3,1:nlmn,ia,ispinor))
         end if
         do ilmn=1,nlmn
           do mu=1,3
             mut = twist_dir(2*mu-1)
             muu = twist_dir(2*mu)
             enlj(mu) = enlj(mu) + aimag(conjg(cft(mut,ilmn))*cfu(muu,ilmn))
             enlj(mu) = enlj(mu) - aimag(conjg(cft(muu,ilmn))*cfu(mut,ilmn))
           end do
         end do
       end do
     end do
     enlout(1:3)=enlout(1:3)+enlj(1:3)
     ABI_DEALLOCATE(cft)
     ABI_DEALLOCATE(cfu)
   end if

!  ====== Accumulate the effective charges contributions =========
   if (choice==54) then
     ABI_ALLOCATE(enljj,(18))
     do ispinor=1,nspinor
       do ia=1,nincat
         enljj(1:18)=zero
         iashift=18*(ia+ia3-2)
!        If cplex=1, dgxdt is real for atm. pos, pure imaginary for k;
!        If cplex_fac=1, dgxdtfac is pure imaginary for k;
         if(cplex==2.and.cplex_fac==2) then
           do ilmn=1,nlmn
             mu=1;nu=1
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 enljj(nu)=enljj(nu) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(1,3+mub,ilmn,ia,ispinor) &
&                 +dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac(2,3+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor) &
&                 +gxfac(2,ilmn,ia,ispinor)*d2gxdt(2,mu,ilmn,ia,ispinor)
                 enljj(nu+1)=enljj(nu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(2,3+mub,ilmn,ia,ispinor) &
&                 -dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac(1,3+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gxdt(2,mu,ilmn,ia,ispinor) &
&                 -gxfac(2,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor)
                 mu=mu+1;nu=nu+2
               end do
             end do
           end do
         else if(cplex==1.and.cplex_fac==2)then
           do ilmn=1,nlmn
             mu=1;nu=1
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 enljj(nu)=enljj(nu) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(1,3+mub,ilmn,ia,ispinor) &
&                 +gxfac(2,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor)
                 enljj(nu+1)=enljj(nu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(2,3+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor)
                 mu=mu+1;nu=nu+2
               end do
             end do
           end do
         else if(cplex==1.and.cplex_fac==1)then
           do ilmn=1,nlmn
             mu=1;nu=1
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 enljj(nu+1)=enljj(nu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(1,3+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor)
                 mu=mu+1;nu=nu+2
               end do
             end do
           end do
         end if
         enlout(iashift+1:iashift+18)=enlout(iashift+1:iashift+18)+enljj(1:18)
       end do
     end do
     ABI_DEALLOCATE(enljj)
   end if

!  ====== Accumulate the piezoelectric tensor contributions =========
   if (choice==55) then
     ABI_ALLOCATE(enljj,(36))
     do ispinor=1,nspinor
       do ia=1,nincat
         enljj(1:36)=zero;enlj(:)=zero
!        If cplex=1, dgxdt is real for strain, pure imaginary for k;
!        If cplex_fac=1, dgxdtfac is pure imaginary for k;
         if(cplex==2.and.cplex_fac==2) then
           do ilmn=1,nlmn
!            First compute 2nd-derivative contribution
             mu=1
             do mua=1,6 ! strain (lambda,nu)
               mua1=alpha(mua) ! (nu)
               mua2=beta(mua)  ! (lambda)
               do mub=1,3 ! k (mu)
                 muu=3*(gamma(mua1,mub)-1)+mua2
                 mut=3*(gamma(mua2,mub)-1)+mua1
                 d2gx(1:cplex)=half*(d2gxdt(1:cplex,muu,ilmn,ia,ispinor) &
&                 +d2gxdt(1:cplex,mut,ilmn,ia,ispinor))
                 enljj(mu)=enljj(mu) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(1,6+mub,ilmn,ia,ispinor) &
&                 +dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac(2,6+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gx(1)+gxfac(2,ilmn,ia,ispinor)*d2gx(2)
                 enljj(mu+1)=enljj(mu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(2,6+mub,ilmn,ia,ispinor) &
&                 -dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac(1,6+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gx(2)-gxfac(2,ilmn,ia,ispinor)*d2gx(1)
                 mu=mu+2
               end do
             end do
!            Then store 1st-derivative contribution
             mu=1
             do nu=1,3
               enlj(mu  )=enlj(mu  )+gxfac(1,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor) &
&               +gxfac(2,ilmn,ia,ispinor)*dgxdt(2,6+nu,ilmn,ia,ispinor)
               enlj(mu+1)=enlj(mu+1)+gxfac(1,ilmn,ia,ispinor)*dgxdt(2,6+nu,ilmn,ia,ispinor) &
&               -gxfac(2,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor)
               mu=mu+2
             end do
           end do
         else if(cplex==1.and.cplex_fac==2)then
           do ilmn=1,nlmn
!            First compute 2nd-derivative contribution
             mu=1
             do mua=1,6 ! strain (lambda,nu)
               mua1=alpha(mua) ! (nu)
               mua2=beta(mua)  ! (lambda)
               do mub=1,3 ! k (mu)
                 muu=3*(gamma(mua1,mub)-1)+mua2
                 mut=3*(gamma(mua2,mub)-1)+mua1
                 d2gx(1)=half*(d2gxdt(1,muu,ilmn,ia,ispinor)+d2gxdt(1,mut,ilmn,ia,ispinor))
                 enljj(mu)=enljj(mu) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(1,6+mub,ilmn,ia,ispinor) &
&                 +gxfac(2,ilmn,ia,ispinor)*d2gx(1)
                 enljj(mu+1)=enljj(mu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(2,6+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gx(1)
                 mu=mu+2
               end do
             end do
!            Then store 1st-derivative contribution
             mu=1
             do nu=1,3
               enlj(mu  )=enlj(mu  )+gxfac(2,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor)
               enlj(mu+1)=enlj(mu+1)+gxfac(1,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor)
               mu=mu+2
             end do
           end do
         else if(cplex==1.and.cplex_fac==1)then
           do ilmn=1,nlmn
             mu=1
             do mua=1,6 ! strain (lambda,nu)
               mua1=alpha(mua) ! (nu)
               mua2=beta(mua)  ! (lambda)
               do mub=1,3 ! k (mu)
                 muu=3*(gamma(mua1,mub)-1)+mua2
                 mut=3*(gamma(mua2,mub)-1)+mua1
                 d2gx(1)=half*(d2gxdt(1,muu,ilmn,ia,ispinor)+d2gxdt(1,mut,ilmn,ia,ispinor))
                 enljj(mu+1)=enljj(mu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac(1,6+mub,ilmn,ia,ispinor) &
&                 +gxfac(1,ilmn,ia,ispinor)*d2gx(1)
                 mu=mu+2
               end do
             end do
!            Then store 1st-derivative contribution
             mu=1
             do nu=1,3
               enlj(mu+1)=enlj(mu+1)+gxfac(1,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor)
               mu=mu+2
             end do
           end do
         end if
         enlout(1:36)=enlout(1:36)+enljj(1:36)
         ddkk(1:6)=ddkk(1:6)+enlj(1:6)
       end do
     end do
     ABI_DEALLOCATE(enljj)
   end if

!  ======= Accumulate the elastic tensor contributions ==========
   if (choice==6) then
     do ispinor=1,nspinor
       do ia=1,nincat
         iashift=3*(ia+ia3-2)
         do ilmn=1,nlmn
           do iplex=1,cplex
             enlk=enlk+gxfac(iplex,ilmn,ia,ispinor)*gx(iplex,ilmn,ia,ispinor)
           end do
           enlj(1:3)=zero
           do mu=1,3
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,6+mu,ilmn,ia,ispinor)
             end do
           end do
           fnlk(iashift+1:iashift+3)=fnlk(iashift+1:iashift+3)+two*enlj(1:3)
           enlj(1:6)=zero
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
           strnlk(1:6)=strnlk(1:6)+two*enlj(1:6)
           do mub=1,6
             mushift=6*(mub-1);nushift=(3*natom+6)*(mub-1)
             do mua=1,6
               mu=mushift+mua;nu=nushift+mua
               do iplex=1,cplex
                 enlout(nu)=enlout(nu)+two* &
&                 (gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)&
&                 +dgxdtfac(iplex,mua,ilmn,ia,ispinor)*dgxdt(iplex,mub,ilmn,ia,ispinor))
               end do
             end do
             mushift=36+3*(mub-1);nushift=6+iashift+(3*natom+6)*(mub-1)
             do mua=1,3
               mu=mushift+mua;nu=nushift+mua
               do iplex=1,cplex
                 enlout(nu)=enlout(nu)+two* &
&                 (gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)&
&                 +dgxdtfac(iplex,mub,ilmn,ia,ispinor)*dgxdt(iplex,6+mua,ilmn,ia,ispinor))
               end do
             end do
           end do
         end do
       end do
     end do
   end if

!  ======== Accumulate the contributions of 2nd-derivatives of E wrt to k ==========
   if (choice==8) then
     ABI_ALLOCATE(cft,(3,nlmn))
     ABI_ALLOCATE(cfu,(3,nlmn))
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:6)=zero
         do ilmn=1,nlmn
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
         end do
!        If cplex=1, dgxdt is pure imaginary;
!        If cplex_fac=1, dgxdtfac is pure imaginary;
         if(cplex==2)then
           cft(1:3,1:nlmn)=cmplx(dgxdt(1,1:3,1:nlmn,ia,ispinor),dgxdt(2,1:3,1:nlmn,ia,ispinor))
         else
           cft(1:3,1:nlmn)=cmplx(zero,dgxdt(1,1:3,1:nlmn,ia,ispinor))
         end if
         if(cplex_fac==2)then
           cfu(1:3,1:nlmn)=cmplx(dgxdtfac(1,1:3,1:nlmn,ia,ispinor),dgxdtfac(2,1:3,1:nlmn,ia,ispinor))
         else
           cfu(1:3,1:nlmn)=cmplx(zero,dgxdtfac(1,1:3,1:nlmn,ia,ispinor))
         end if
         do ilmn=1,nlmn
           do mu=1,6
             mua=alpha(mu);mub=beta(mu)
             enlj(mu)=enlj(mu)+real(conjg(cfu(mub,ilmn))*cft(mua,ilmn))
           end do
         end do
         enlout(1:6)=enlout(1:6)+two*enlj(1:6)
       end do
     end do
     ABI_DEALLOCATE(cft)
     ABI_DEALLOCATE(cfu)
   end if

!  ======== Accumulate the contributions of partial 2nd-derivatives of E wrt to k ==========
!  Full derivative wrt to k1, right derivative wrt to k2
   if (choice==81) then
     ABI_ALLOCATE(cft,(3,nlmn))
     ABI_ALLOCATE(cfu,(6,nlmn))
     ABI_ALLOCATE(enljj,(18))
     do ispinor=1,nspinor
       do ia=1,nincat
         enljj(1:18)=zero
         if(cplex_fac==2)then !If cplex_fac=1, gxfac is pure real
           cft(1,1:nlmn)=cmplx(gxfac(1,1:nlmn,ia,ispinor),gxfac(2,1:nlmn,ia,ispinor))
         else
           cft(1,1:nlmn)=cmplx(gxfac(1,1:nlmn,ia,ispinor),zero)
         end if
         if(cplex==2)then !If cplex=1, d2gxdt is pure real
           cfu(1:6,1:nlmn)=cmplx(d2gxdt(1,1:6,1:nlmn,ia,ispinor),d2gxdt(2,1:6,1:nlmn,ia,ispinor))
         else
           cfu(1:6,1:nlmn)=cmplx(d2gxdt(1,1:6,1:nlmn,ia,ispinor),zero)
         end if
         do ilmn=1,nlmn
           do mu=1,3
             do nu=1,3
               muu=3*(mu-1)+nu ; mut=gamma(mu,nu)
               enljj(2*muu-1)=enljj(2*muu-1)+ real(conjg(cft(1,ilmn))*cfu(mut,ilmn))
               enljj(2*muu  )=enljj(2*muu  )+aimag(conjg(cft(1,ilmn))*cfu(mut,ilmn))
             end do
           end do
         end do
         if(cplex==2)then !If cplex=1, dgxdt is pure imaginary
           cft(1:3,1:nlmn)=cmplx(dgxdt(1,1:3,1:nlmn,ia,ispinor),dgxdt(2,1:3,1:nlmn,ia,ispinor))
         else
           cft(1:3,1:nlmn)=cmplx(zero,dgxdt(1,1:3,1:nlmn,ia,ispinor))
         end if
         if(cplex_fac==2)then !If cplex_fac=1, dgxdtfac is pure imaginary
           cfu(1:3,1:nlmn)=cmplx(dgxdtfac(1,1:3,1:nlmn,ia,ispinor),dgxdtfac(2,1:3,1:nlmn,ia,ispinor))
         else
           cfu(1:3,1:nlmn)=cmplx(zero,dgxdtfac(1,1:3,1:nlmn,ia,ispinor))
         end if
         do ilmn=1,nlmn
           do mu=1,3
             do nu=1,3
               muu=3*(mu-1)+nu
               enljj(2*muu-1)=enljj(2*muu-1)+ real(conjg(cft(mu,ilmn))*cfu(nu,ilmn))
               enljj(2*muu  )=enljj(2*muu  )+aimag(conjg(cft(mu,ilmn))*cfu(nu,ilmn))
             end do
           end do
         end do
         enlout(1:18)=enlout(1:18)+enljj(1:18)
       end do
     end do
     ABI_DEALLOCATE(cft)
     ABI_DEALLOCATE(cfu)
     ABI_DEALLOCATE(enljj)
   end if

 end if

 if (paw_opt==3) then

!  ============== Accumulate contribution to <c|S|c> ===============
   if (choice==1) then
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           do iplex=1,cplex
             enlout(1)=enlout(1)+gxfac_sij(iplex,ilmn,ia,ispinor)*gx(iplex,ilmn,ia,ispinor)
           end do
         end do
       end do
     end do
   end if

!  ============== Accumulate contribution to <c|dS/d_atm_pos|c> ===============
   if (choice==2.or.choice==23) then
     ishift=0;if (choice==23) ishift=6
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:3)=zero
         iashift=3*(ia+ia3-2)
         do ilmn=1,nlmn
           do mu=1,3
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac_sij(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu+ishift,ilmn,ia,ispinor)
             end do
           end do
         end do
         enlout(iashift+1:iashift+3)=enlout(iashift+1:iashift+3)+two*enlj(1:3)
       end do
     end do
   end if

!  ============== Accumulate contribution to <c|dS/d_strain|c> ===============
   if (choice==3.or.choice==23) then
     enlj(1:6)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           gxfacj(1:cplex)=gxfac_sij(1:cplex,ilmn,ia,ispinor)
           do iplex=1,cplex
             enlk=enlk+gxfacj(iplex)*gx(iplex,ilmn,ia,ispinor)
           end do
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfacj(iplex)*dgxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
         end do
       end do
     end do
     enlout(1:6)=enlout(1:6)+two*enlj(1:6)
   end if

!  ======== Accumulate the contributions of derivatives of <c|S|c> wrt to k ==========
   if (choice==5) then
     enlj(1:3)=zero
!    If cplex=1, gxfac is real and dgxdt is pure imaginary; thus there is no contribution
     if(cplex==2)then
       do ispinor=1,nspinor
         do ia=1,nincat
           do ilmn=1,nlmn
             do mu=1,3
               do iplex=1,cplex
                 enlj(mu)=enlj(mu)+gxfac_sij(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu,ilmn,ia,ispinor)
               end do
             end do
           end do
         end do
       end do
     end if
     enlout(1:3)=enlout(1:3)+two*enlj(1:3)
   end if

!  ====== Accumulate the contributions of left or right derivatives of <c|S|c> wrt to k ==========
!  Choice 51: right derivative wrt to k ; Choice 52: left derivative wrt to k
   if (choice==51.or.choice==52) then
     enlj(1:6)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         if(cplex==2)then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(2*mu-1)=enlj(2*mu-1)+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor) &
&               +gxfac_sij(2,ilmn,ia,ispinor)*dgxdt(2,mu,ilmn,ia,ispinor)
               enlj(2*mu  )=enlj(2*mu  )+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(2,mu,ilmn,ia,ispinor) &
&               -gxfac_sij(2,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         else if (cplex_fac==2) then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(2*mu-1)=enlj(2*mu-1)+gxfac_sij(2,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
               enlj(2*mu  )=enlj(2*mu  )+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         else if (cplex_fac==1) then
           do ilmn=1,nlmn
             do mu=1,3
               enlj(2*mu  )=enlj(2*mu  )+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(1,mu,ilmn,ia,ispinor)
             end do
           end do
         end if
       end do
     end do
     if (choice==52) then
       enlj(2)=-enlj(2);enlj(4)=-enlj(4);enlj(6)=-enlj(6)
     end if
     enlout(1:6)=enlout(1:6)+enlj(1:6)
   end if

!  ====== Accumulate contribution to <c|d2S/d_atm_pos d_left_k|c> =========
   if (choice==54) then
     ABI_ALLOCATE(enljj,(18))
     do ispinor=1,nspinor
       do ia=1,nincat
         enljj(1:18)=zero
         iashift=18*(ia+ia3-2)
         if(cplex==2) then
           do ilmn=1,nlmn
             mu=1;nu=1
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 enljj(nu)=enljj(nu) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac_sij(1,3+mub,ilmn,ia,ispinor) &
&                 +dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac_sij(2,3+mub,ilmn,ia,ispinor) &
&                 +gxfac_sij(1,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor) &
&                 +gxfac_sij(2,ilmn,ia,ispinor)*d2gxdt(2,mu,ilmn,ia,ispinor)

                 enljj(nu+1)=enljj(nu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac_sij(2,3+mub,ilmn,ia,ispinor) &
&                 -dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac_sij(1,3+mub,ilmn,ia,ispinor) &
&                 +gxfac_sij(1,ilmn,ia,ispinor)*d2gxdt(2,mu,ilmn,ia,ispinor) &
&                 -gxfac_sij(2,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor)

                 mu=mu+1;nu=nu+2
               end do
             end do
           end do
!        If cplex=1, dgxdt, d2gxdt and dgxdtfac_sij are real for atm. pos, pure imaginary for k
         else
           do ilmn=1,nlmn
             mu=1;nu=1
             do mua=1,3 ! atm. pos
               do mub=1,3 ! k
                 enljj(nu+1)=enljj(nu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac_sij(1,3+mub,ilmn,ia,ispinor) &
&                 +gxfac_sij(1,ilmn,ia,ispinor)*d2gxdt(1,mu,ilmn,ia,ispinor)
                 mu=mu+1;nu=nu+2
               end do
             end do
           end do
         end if
         enlout(iashift+1:iashift+18)=enlout(iashift+1:iashift+18)+enljj(1:18)
       end do
     end do
     ABI_DEALLOCATE(enljj)
   end if

!  ====== Accumulate contribution to <c|d2S/d_dstrain d_right_k|c> =========
   if (choice==55) then
     ABI_ALLOCATE(enljj,(36))
     do ispinor=1,nspinor
       do ia=1,nincat
         enljj(1:36)=zero;enlj(:)=zero
!        If cplex=1, dgxdt is real for strain, pure imaginary for k;
!        If cplex_fac=1, dgxdtfac is pure imaginary for k;
         if(cplex==2.and.cplex_fac==2) then
           do ilmn=1,nlmn
!            First compute 2nd-derivative contribution
             mu=1
             do mua=1,6 ! strain (lambda,nu)
               mua1=alpha(mua) ! (nu)
               mua2=beta(mua)  ! (lambda)
               do mub=1,3 ! k (mu)
                 muu=3*(gamma(mua1,mub)-1)+mua2
                 mut=3*(gamma(mua2,mub)-1)+mua1
                 d2gx(1:cplex)=half*(d2gxdt(1:cplex,muu,ilmn,ia,ispinor) &
&                 +d2gxdt(1:cplex,mut,ilmn,ia,ispinor))
                 enljj(mu)=enljj(mu) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac_sij(1,6+mub,ilmn,ia,ispinor) &
&                 +dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac_sij(2,6+mub,ilmn,ia,ispinor) &
&                 +gxfac_sij(1,ilmn,ia,ispinor)*d2gx(1)+gxfac_sij(2,ilmn,ia,ispinor)*d2gx(2)
                 enljj(mu+1)=enljj(mu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac_sij(2,6+mub,ilmn,ia,ispinor) &
&                 -dgxdt(2,mua,ilmn,ia,ispinor)*dgxdtfac_sij(1,6+mub,ilmn,ia,ispinor) &
&                 +gxfac_sij(1,ilmn,ia,ispinor)*d2gx(2)-gxfac_sij(2,ilmn,ia,ispinor)*d2gx(1)
                 mu=mu+2
               end do
             end do
!            Then store 1st-derivative contribution
             mu=1
             do nu=1,3
               enlj(mu  )=enlj(mu  )+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor) &
&               +gxfac_sij(2,ilmn,ia,ispinor)*dgxdt(2,6+nu,ilmn,ia,ispinor)
               enlj(mu+1)=enlj(mu+1)+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(2,6+nu,ilmn,ia,ispinor) &
&               -gxfac_sij(2,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor)
               mu=mu+2
             end do
           end do
!        If cplex=1, dgxdt, d2gxdt and dgxdtfac_sij are real for atm. pos, pure imaginary for k
         else
           do ilmn=1,nlmn
             mu=1
             do mua=1,6 ! strain (lambda,nu)
               mua1=alpha(mua) ! (nu)
               mua2=beta(mua)  ! (lambda)
               do mub=1,3 ! k (mu)
                 muu=3*(gamma(mua1,mub)-1)+mua2
                 mut=3*(gamma(mua2,mub)-1)+mua1
                 d2gx(1)=half*(d2gxdt(1,muu,ilmn,ia,ispinor)+d2gxdt(1,mut,ilmn,ia,ispinor))
                 enljj(mu+1)=enljj(mu+1) &
&                 +dgxdt(1,mua,ilmn,ia,ispinor)*dgxdtfac_sij(1,6+mub,ilmn,ia,ispinor) &
&                 +gxfac_sij(1,ilmn,ia,ispinor)*d2gx(1)
                 mu=mu+2
               end do
             end do
!            Then store 1st-derivative contribution
             mu=1
             do nu=1,3
               enlj(mu+1)=enlj(mu+1)+gxfac_sij(1,ilmn,ia,ispinor)*dgxdt(1,6+nu,ilmn,ia,ispinor)
               mu=mu+2
             end do
           end do
         end if
         enlout(1:36)=enlout(1:36)+enljj(1:36)
         ddkk(1:6)=ddkk(1:6)+enlj(1:6)
       end do
     end do
     ABI_DEALLOCATE(enljj)
   end if

!  ======  Accumulate contribution to <c|d2S/d_k d_k|c> =========
   if (choice==8) then
     ABI_ALLOCATE(cft,(3,nlmn))
     ABI_ALLOCATE(cfu,(3,nlmn))
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:6)=zero
         do ilmn=1,nlmn
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac_sij(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
         end do
!        If cplex=1, dgxdt is pure imaginary, dgxdtfac_sij is pure imaginary;
         if(cplex==2)then
           cft(1:3,1:nlmn)=cmplx(dgxdt(1,1:3,1:nlmn,ia,ispinor),dgxdt(2,1:3,1:nlmn,ia,ispinor))
           cfu(1:3,1:nlmn)=cmplx(dgxdtfac_sij(1,1:3,1:nlmn,ia,ispinor),dgxdtfac_sij(2,1:3,1:nlmn,ia,ispinor))
         else
           cft(1:3,1:nlmn)=cmplx(zero,dgxdt(1,1:3,1:nlmn,ia,ispinor))
           cfu(1:3,1:nlmn)=cmplx(zero,dgxdtfac_sij(1,1:3,1:nlmn,ia,ispinor))
         end if
         do ilmn=1,nlmn
           do mu=1,6
             mua=alpha(mu);mub=beta(mu)
             enlj(mu)=enlj(mu)+real(conjg(cfu(mub,ilmn))*cft(mua,ilmn))
           end do
         end do
         enlout(1:6)=enlout(1:6)+two*enlj(1:6)
       end do
     end do
     ABI_DEALLOCATE(cft)
     ABI_DEALLOCATE(cfu)
   end if

!  ======  Accumulate contribution to <c|d/d_k[d(right)S/d_k]|c> =========
!  Full derivative wrt to k1, right derivative wrt to k2
   if (choice==81) then
     ABI_ALLOCATE(cft,(3,nlmn))
     ABI_ALLOCATE(cfu,(6,nlmn))
     ABI_ALLOCATE(enljj,(18))
     do ispinor=1,nspinor
       do ia=1,nincat
         enljj(1:18)=zero
         if(cplex_fac==2)then !If cplex_fac=1, gxfac is pure real
           cft(1,1:nlmn)=cmplx(gxfac_sij(1,1:nlmn,ia,ispinor),gxfac_sij(2,1:nlmn,ia,ispinor))
         else
           cft(1,1:nlmn)=cmplx(gxfac_sij(1,1:nlmn,ia,ispinor),zero)
         end if
         if(cplex==2)then !If cplex=1, d2gxdt is pure real
           cfu(1:6,1:nlmn)=cmplx(d2gxdt(1,1:6,1:nlmn,ia,ispinor),d2gxdt(2,1:6,1:nlmn,ia,ispinor))
         else
           cfu(1:6,1:nlmn)=cmplx(d2gxdt(1,1:6,1:nlmn,ia,ispinor),zero)
         end if
         do ilmn=1,nlmn
           do mu=1,3
             do nu=1,3
               muu=3*(mu-1)+nu ; mut=gamma(mu,nu)
               enljj(2*muu-1)=enljj(2*muu-1)+ real(conjg(cft(1,ilmn))*cfu(mut,ilmn))
               enljj(2*muu  )=enljj(2*muu  )+aimag(conjg(cft(1,ilmn))*cfu(mut,ilmn))
             end do
           end do
         end do
         if(cplex==2)then !If cplex=1, dgxdt is pure imaginary
           cft(1:3,1:nlmn)=cmplx(dgxdt(1,1:3,1:nlmn,ia,ispinor),dgxdt(2,1:3,1:nlmn,ia,ispinor))
         else
           cft(1:3,1:nlmn)=cmplx(zero,dgxdt(1,1:3,1:nlmn,ia,ispinor))
         end if
         if(cplex_fac==2)then !If cplex_fac=1, dgxdtfac is pure imaginary
           cfu(1:3,1:nlmn)=cmplx(dgxdtfac_sij(1,1:3,1:nlmn,ia,ispinor),dgxdtfac_sij(2,1:3,1:nlmn,ia,ispinor))
         else
           cfu(1:3,1:nlmn)=cmplx(zero,dgxdtfac_sij(1,1:3,1:nlmn,ia,ispinor))
         end if
         do ilmn=1,nlmn
           do mu=1,3
             do nu=1,3
               muu=3*(mu-1)+nu
               enljj(2*muu-1)=enljj(2*muu-1)+ real(conjg(cft(mu,ilmn))*cfu(nu,ilmn))
               enljj(2*muu  )=enljj(2*muu  )+aimag(conjg(cft(mu,ilmn))*cfu(nu,ilmn))
             end do
           end do
         end do
         enlout(1:18)=enlout(1:18)+enljj(1:18)
       end do
     end do
     ABI_DEALLOCATE(cft)
     ABI_DEALLOCATE(cfu)
     ABI_DEALLOCATE(enljj)
   end if

 end if

end subroutine opernld_ylm
!!***

end module m_opernld_ylm
!!***

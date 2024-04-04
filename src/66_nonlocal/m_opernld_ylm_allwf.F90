!!****m* ABINIT/m_opernld_ylm_allwf
!! NAME
!!  m_opernld_ylm_allwf
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opernld_ylm_allwf

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_geometry, only : strconv
 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: opernld_ylm_allwf
!!***

contains
!!***

!!****f* ABINIT/opernld_ylm_allwf
!! NAME
!! opernld_ylm_allwf
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
!!      if choice=53: enlout(3)             -> 1st deriv. (twist) of energy wrt k
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
!! SOURCE

subroutine opernld_ylm_allwf(choice,cplex,cplex_fac,ddkk,&
&                      dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,&
&                      enlk,enlout,gx,gxfac,gxfac_sij,ndat,nd2gxdt,ndgxdt,&
&                      ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
&                      nattyp,gpu_option)

 ! Arguments ------------------------------------
 ! scalars
 integer,intent(in) :: choice,paw_opt,ntypat,ndgxdtfac,nd2gxdt,ndgxdt
 integer,intent(in) :: cplex,cplex_fac,ndat,nnlout,nspinor,nprojs,lmnmax,gpu_option

 ! arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),nattyp(ntypat)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdtfac(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex_fac,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: gx(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex_fac,nprojs,ndat*nspinor)
 real(dp),intent(inout) :: enlout(nnlout*ndat)
 real(dp),intent(inout) :: enlk(ndat),ddkk(6,ndat)

 ! locals
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer :: force_shift, shift, nattyp_i
 integer :: itypat, ilmn, ia, idat, igrad, ii, nlmn, iend, ibeg, iatm, iashift
 integer :: mua, mub, nu, mu, mua1, mua2, muu, mut
 real(dp) :: esum,esumi
 real(dp) :: d2gx(cplex)

 ABI_UNUSED((/gpu_option/))

 enlout=zero
 if(paw_opt < 3) then
   if(choice==1.or.choice==3.or.choice==23) then
     shift=0
     iatm=0
     esum=zero
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       ibeg = shift+1
       iend = shift+nattyp(itypat)*nlmn
       nattyp_i = nattyp(itypat)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& MAP(to:enlk,gxfac,gx) &
       !$OMP& PRIVATE(idat,esum) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat*nspinor
         esum=zero
         !$OMP PARALLEL DO REDUCTION(+:esum) PRIVATE(ia,ilmn,ii)
         do ia=1,nattyp_i
           do ilmn=1,nlmn
             do ii=1,cplex
               esum=esum +gxfac(ii,shift+(ia-1)*nlmn+ilmn,idat) &
               &         *gx   (ii,shift+(ia-1)*nlmn+ilmn,idat)
             end do
           end do
         end do
         enlk(idat) = enlk(idat) + esum
       end do

       shift = shift + nattyp(itypat)*nlmn
       iatm = iatm+nattyp(itypat)
     end do
     if (choice==1) then
       !$OMP TARGET PARALLEL DO MAP(to:enlout,enlk) PRIVATE(idat) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
         enlout(idat)=enlk(idat)
       end do
     end if
   end if ! choice=1/3/23

!  ======== Accumulate the stress tensor contributions ==========
   if (choice==3.or.choice==23) then
     shift=0
     iatm=0
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       ibeg = shift+1
       iend = shift+nattyp(itypat)*nlmn
       nattyp_i = nattyp(itypat)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& MAP(to:enlout,gxfac,dgxdt) &
       !$OMP& PRIVATE(idat,igrad,esum) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat*nspinor
         do igrad=1,6
           esum=zero
           !$OMP PARALLEL DO REDUCTION(+:esum) PRIVATE(ia,ilmn,ii)
           do ia=1,nattyp_i
             !Following loops are a [D][Z]DOT
             do ilmn=1,nlmn
               do ii=1,cplex
                 esum=esum +gxfac(ii,shift+(ia-1)*nlmn+ilmn,idat) &
                 &         *dgxdt(ii,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + igrad,idat)
               end do
             end do
           end do
           enlout((idat-1)*nnlout+igrad) = enlout((idat-1)*nnlout+igrad) + two*esum
         end do
       end do

       shift = shift + nattyp(itypat)*nlmn
       iatm = iatm+nattyp(itypat)
     end do
   end if

!  ============ Accumulate the forces contributions =============
   if (choice==2.or.choice==23) then
     shift=0
     iatm=0
     force_shift=0; if(choice==23) force_shift=6
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       nattyp_i = nattyp(itypat)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& MAP(to:enlout,gxfac,dgxdt) &
       !$OMP& PRIVATE(idat,igrad,ia,esum) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat*nspinor
         do ia=1,nattyp_i
           do igrad=1,3
             !Following loops are a [D][Z]DOT
             esum=zero
             !$OMP PARALLEL DO REDUCTION(+:esum) PRIVATE(ilmn,ii)
             do ilmn=1,nlmn
               do ii=1,cplex
                 esum=esum +gxfac(ii,shift+(ia-1)*nlmn+ilmn,idat) &
                 &         *dgxdt(ii,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + (igrad+force_shift) ,idat)
               end do
             end do
             enlout((idat-1)*nnlout + force_shift + (iatm+ia-1)*3 + igrad)= &
             &             enlout((idat-1)*nnlout + force_shift + (iatm+ia-1)*3 + igrad) + two*esum
           end do
         end do
       end do

       shift = shift + nattyp(itypat)*nlmn
       iatm = iatm+nattyp(itypat)
     end do
   end if
 end if ! paw_opt < 3

 if(paw_opt==3) then

!  ====== Accumulate contribution to <c|d2S/d_atm_pos d_left_k|c> =========
   if (choice==54) then
     shift=0
     iatm=0
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       nattyp_i = nattyp(itypat)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& MAP(to:enlout,gxfac_sij,dgxdt,dgxdtfac_sij,d2gxdt) &
       !$OMP& PRIVATE(idat,igrad,ia,esum,esumi,mu,nu,mua,mub,iashift) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat*nspinor
       do ia=1,nattyp_i
         iashift=18*(ia+iatm-1)
         if(cplex==2) then
           do mua=1,3 ! atm. pos
             do mub=1,3 ! k
               mu=(mua-1)*3+mub
               nu=(mua-1)*6+mub*2-1
               esum=zero; esumi=zero
               !$OMP PARALLEL DO REDUCTION(+:esum,esumi) PRIVATE(ilmn)
               do ilmn=1,nlmn
                 esum=esum &
&                 +dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)&
                    *dgxdtfac_sij(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 3+mub,idat) &
&                 +dgxdt(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)&
&                   *dgxdtfac_sij(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 3+mub,idat) &
&                 +gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*d2gxdt(1,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + mu,idat) &
&                 +gxfac_sij(2,shift+(ia-1)*nlmn+ilmn,idat)*d2gxdt(2,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + mu,idat)

                 esumi=esumi &
&                 +dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)&
                    *dgxdtfac_sij(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 3+mub,idat) &
&                 -dgxdt(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)&
                    *dgxdtfac_sij(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 3+mub,idat) &
&                 +gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*d2gxdt(2,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + mu,idat) &
&                 -gxfac_sij(2,shift+(ia-1)*nlmn+ilmn,idat)*d2gxdt(1,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + mu,idat)
               end do
               enlout(nnlout*(idat-1)+iashift+nu)=enlout(nnlout*(idat-1)+iashift+nu)+esum
               enlout(nnlout*(idat-1)+iashift+nu+1)=enlout(nnlout*(idat-1)+iashift+nu+1)+esumi
             end do
           end do
!        If cplex=1, dgxdt, d2gxdt and dgxdtfac_sij are real for atm. pos, pure imaginary for k
         else
           do mua=1,3 ! atm. pos
             do mub=1,3 ! k
               mu=(mua-1)*3+mub
               nu=(mua-1)*6+mub*2-1
               esumi=zero
               !$OMP PARALLEL DO REDUCTION(+:esumi) PRIVATE(ilmn)
               do ilmn=1,nlmn
                 esumi=esumi &
&                 +dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)&
&                   *dgxdtfac_sij(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 3+mub,idat) &
&                 +gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*d2gxdt(1,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + mu,idat)
               end do
               enlout(nnlout*(idat-1)+iashift+nu+1)=enlout(nnlout*(idat-1)+iashift+nu+1)+esumi
             end do
           end do
         end if
       end do ! ia
       end do ! idat

       shift = shift + nattyp(itypat)*nlmn
       iatm = iatm+nattyp(itypat)
     end do ! itypat
   end if

!  ====== Accumulate contribution to <c|d2S/d_dstrain d_right_k|c> =========
   if (choice==55) then
     shift=0
     iatm=0
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       nattyp_i = nattyp(itypat)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& MAP(to:enlout,gxfac_sij,dgxdt,dgxdtfac_sij,d2gxdt) &
       !$OMP& PRIVATE(idat,igrad,ia,esum,esumi,mu,nu,mua,mub,mua1,mua2,muu,mut) &
       !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
       do idat=1,ndat
       do ia=1,nattyp_i
!        If cplex=1, dgxdt is real for strain, pure imaginary for k;
!        If cplex_fac=1, dgxdtfac is pure imaginary for k;
         if(cplex==2.and.cplex_fac==2) then
!          First compute 2nd-derivative contribution
           do mua=1,6 ! strain (lambda,nu)
             do mub=1,3 ! k (mu)
               mu=(mua-1)*6+mub*2-1
               esum=zero; esumi=zero
               mua1=alpha(mua) ! (nu)
               mua2=beta(mua)  ! (lambda)
               muu=3*(gamma(mua1,mub)-1)+mua2
               mut=3*(gamma(mua2,mub)-1)+mua1
               !$OMP PARALLEL DO REDUCTION(+:esum,esumi) PRIVATE(ilmn,d2gx)
               do ilmn=1,nlmn
                 d2gx(1:cplex)=half*(d2gxdt(1:cplex,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + muu,idat) &
&                 +d2gxdt(1:cplex,nd2gxdt*shift + (ia-1)*nlmn*nd2gxdt + (ilmn-1)*nd2gxdt + mut,idat))
                 esum=esum &
&                 +dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)*dgxdtfac_sij(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+mub,idat) &
&                 +dgxdt(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)*dgxdtfac_sij(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+mub,idat) &
&                 +gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*d2gx(1)+gxfac_sij(2,shift+(ia-1)*nlmn+ilmn,idat)*d2gx(2)
                 esumi=esumi &
&                 +dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)*dgxdtfac_sij(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+mub,idat) &
&                 -dgxdt(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + mua,idat)*dgxdtfac_sij(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+mub,idat) &
&                 +gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*d2gx(2)-gxfac_sij(2,shift+(ia-1)*nlmn+ilmn,idat)*d2gx(1)
               end do
               enlout(nnlout*(idat-1)+mu)   = enlout(nnlout*(idat-1)+mu)+esum
               enlout(nnlout*(idat-1)+mu+1) = enlout(nnlout*(idat-1)+mu+1)+esumi
             end do
           end do
!          Then store 1st-derivative contribution
           do nu=1,3
             mu=nu*2-1
             esum=zero; esumi=zero
             !$OMP PARALLEL DO REDUCTION(+:esum,esumi) PRIVATE(ilmn)
             do ilmn=1,nlmn
               esum=esum+gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+nu,idat) &
&               +gxfac_sij(2,shift+(ia-1)*nlmn+ilmn,idat)*dgxdt(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+nu,idat)
               esumi=esumi+gxfac_sij(1,shift+(ia-1)*nlmn+ilmn,idat)*dgxdt(2,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+nu,idat) &
&               -gxfac_sij(2,shift+(ia-1)*nlmn+ilmn,idat)*dgxdt(1,ndgxdt*shift + (ia-1)*nlmn*ndgxdt + (ilmn-1)*ndgxdt + 6+nu,idat)
             end do
             ddkk(mu,idat)  = ddkk(mu,idat)+esum
             ddkk(mu+1,idat)= ddkk(mu+1,idat)+esumi
           end do
!        If cplex=1, dgxdt, d2gxdt and dgxdtfac_sij are real for atm. pos, pure imaginary for k
         else
          ! ABI_BUG("Not implemented")
          ! do ilmn=1,nlmn
          !   mu=1
          !   do mua=1,6 ! strain (lambda,nu)
          !     mua1=alpha(mua) ! (nu)
          !     mua2=beta(mua)  ! (lambda)
          !     do mub=1,3 ! k (mu)
          !       muu=3*(gamma(mua1,mub)-1)+mua2
          !       mut=3*(gamma(mua2,mub)-1)+mua1
          !       d2gx(1)=half*(d2gxdt(1,muu,ilmn,ia,idat)+d2gxdt(1,mut,ilmn,ia,idat))
          !       enljj(mu+1)=enljj(mu+1) &
          !        +dgxdt(1,mua,ilmn,ia,idat)*dgxdtfac_sij(1,6+mub,ilmn,ia,idat) &
          !        +gxfac_sij(1,ilmn,ia,idat)*d2gx(1)
          !       mu=mu+2
          !     end do
          !   end do
!         !   Then store 1st-derivative contribution
          !   mu=1
          !   do nu=1,3
          !     enlj(mu+1)=enlj(mu+1)+gxfac_sij(1,ilmn,ia,idat)*dgxdt(1,6+nu,ilmn,ia,idat)
          !     mu=mu+2
          !   end do
          ! end do
         end if
       end do ! ia
       end do ! idat

       shift = shift + nattyp(itypat)*nlmn
       iatm = iatm+nattyp(itypat)
     end do ! itypat
   end if

 end if ! paw_opt == 3

end subroutine opernld_ylm_allwf
!!***

end module m_opernld_ylm_allwf
!!***

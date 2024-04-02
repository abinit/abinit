!!****m* ABINIT/m_opernld_ylm_allwf_cpu
!! NAME
!!  m_opernld_ylm
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

module m_opernld_ylm_allwf_cpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_geometry, only : strconv
 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: opernld_ylm_allwf_cpu
!!***

contains
!!***

!!****f* ABINIT/opernld_ylm_allwf_cpu
!! NAME
!! opernld_ylm_allwf_cpu
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

subroutine opernld_ylm_allwf_cpu(choice,cplex,cplex_fac,&
&                      dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,&
&                      enlk,enlout,gx,gxfac,gxfac_sij,ndat,nd2gxdt,ndgxdt,&
&                      ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
&                      nattyp)

 ! Arguments ------------------------------------
 ! scalars
 integer,intent(in) :: choice,paw_opt,ntypat,ndgxdtfac,nd2gxdt,ndgxdt
 integer,intent(in) :: cplex,cplex_fac,ndat,nnlout,nspinor,nprojs,lmnmax

 ! arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),nattyp(ntypat)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdtfac(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex_fac,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: gx(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex_fac,nprojs,ndat*nspinor)
 real(dp),intent(inout) :: enlout(nnlout*ndat)
 real(dp),intent(inout) :: enlk(ndat)

 ! locals
 integer :: enlout_shift, force_shift, atom_force_shift, grad_shift, proj_shift, nnlout_test
 integer :: itypat, ilmn, ia, idat, ierr, igrad, ii, nlmn
 real(dp) :: esum
 real(dp) :: work(6)
 integer :: i1,i2,i3

 enlout=zero
 enlk=zero
 if(choice==1.or.choice==3.or.choice==23) then
   do idat=1,ndat*nspinor
     proj_shift=0
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       do ia=1,nattyp(itypat)
         !Following loops are a [D][Z]DOT
         esum=zero
         do ilmn=1,nlmn
           do ii=1,cplex
             esum=esum +gxfac(ii,proj_shift+ilmn,idat) &
             &         *gx   (ii,proj_shift+ilmn,idat)
           end do
         end do
         proj_shift=proj_shift+nlmn
         enlk(idat) = enlk(idat) + esum
       end do
     end do
   end do
   if (choice==1) enlout(1:ndat)=enlk(1:ndat)
 end if ! choice=1/3/23
 if(choice==2.or.choice==3.or.choice==23) then
   do idat=1,ndat*nspinor
     proj_shift=0 ; grad_shift=0
     enlout_shift=(idat-1)*nnlout
     atom_force_shift=merge(6,0,choice==23)
     force_shift=merge(6,0,choice==23)
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       do ia=1,nattyp(itypat)
         if (choice==3.or.choice==23) then
           do igrad=1,6
             !Following loops are a [D][Z]DOT
             esum=zero
             do ilmn=1,nlmn
               do ii=1,cplex
                 esum=esum +gxfac(ii,proj_shift+ilmn,idat) &
                 &         *dgxdt(ii,grad_shift+(ilmn-1)*ndgxdt+igrad,idat)
               end do
             end do
             enlout(enlout_shift+igrad)=enlout(enlout_shift+igrad) + two*esum
           end do
         end if
         if (choice==2.or.choice==23) then
           do igrad=1,3
             !Following loops are a [D][Z]DOT
             esum=zero
             do ilmn=1,nlmn
               do ii=1,cplex
                 esum=esum +gxfac(ii,proj_shift+ilmn,idat) &
                   &       *dgxdt(ii,grad_shift+(ilmn-1)*ndgxdt+force_shift+igrad,idat)
               end do
             end do
             enlout(enlout_shift+atom_force_shift+igrad)= &
             &    enlout(enlout_shift+atom_force_shift+igrad) + two*esum
           end do
           atom_force_shift=atom_force_shift+3
         end if
         grad_shift=grad_shift+nlmn*ndgxdt
         proj_shift=proj_shift+nlmn
       end do
     end do
   end do
 end if ! choice=2, 3 or 23


end subroutine opernld_ylm_allwf_cpu

end module m_opernld_ylm_allwf_cpu

!!****m* ABINIT/m_opernld_ylm_allwf_ompgpu
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

module m_opernld_ylm_allwf_ompgpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_geometry, only : strconv
 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: opernld_ylm_allwf_ompgpu
!!***

contains
!!***

!!****f* ABINIT/opernld_ylm_allwf_ompgpu
!! NAME
!! opernld_ylm_allwf_ompgpu
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

subroutine opernld_ylm_allwf_ompgpu(choice,cplex,cplex_fac,&
&                      dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,&
&                      enlout,gx,gxfac,gxfac_sij,ndat,nd2gxdt,ndgxdt,&
&                      ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
&                      gprimd,nattyp,mpi_enreg)

 ! Arguments ------------------------------------
 ! scalars
 integer,intent(in) :: choice,paw_opt,ntypat,ndgxdtfac,nd2gxdt,ndgxdt
 integer,intent(in) :: cplex,cplex_fac,ndat,nnlout,nspinor,nprojs,lmnmax
 real(dp),intent(inout) :: enlout(nnlout*ndat)
 type(MPI_type),intent(in) :: mpi_enreg

 ! arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),nattyp(ntypat)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdtfac(cplex,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: dgxdtfac_sij(cplex_fac,ndgxdt*nprojs,ndat*nspinor)
 real(dp),intent(in) :: gx(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac(cplex,nprojs,ndat*nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex_fac,nprojs,ndat*nspinor)

 ! locals
 integer :: enlout_shift, force_shift, grad_shift, proj_shift, nnlout_test, shift, nattyp_i
 integer :: itypat, ilmn, ia, idat, ierr, igrad, ii, nlmn, iend, ibeg, iatm
 real(dp) :: esum
 real(dp) :: work(6)
 real(dp), allocatable :: enlk(:)

 enlout=zero
 !$OMP TARGET ENTER DATA MAP(to:nattyp)
 if(choice==1.or.choice==3.or.choice==23) then
   ABI_MALLOC(enlk,(ndat))
   enlk=zero
   !$OMP TARGET ENTER DATA MAP(to:enlk)
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
     !$OMP& PRIVATE(idat,esum)
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
     !$OMP TARGET PARALLEL DO MAP(to:enlout,enlk) PRIVATE(idat)
     do idat=1,ndat
       enlout(idat)=enlk(idat)
     end do
   end if
 end if ! choice=1/3/23
 if(choice==2.or.choice==3.or.choice==23) then
   grad_shift=merge(9,6,choice==23)
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
       !$OMP& PRIVATE(idat,igrad,esum)
       do idat=1,ndat*nspinor
         do igrad=1,6
           esum=zero
           !$OMP PARALLEL DO REDUCTION(+:esum) PRIVATE(ia,ilmn,ii)
           do ia=1,nattyp_i
             !Following loops are a [D][Z]DOT
             do ilmn=1,nlmn
               do ii=1,cplex
                 esum=esum +gxfac(ii,shift+(ia-1)*nlmn+ilmn,idat) &
                 &         *dgxdt(ii,grad_shift*shift + (ia-1)*nlmn*grad_shift + (igrad-1)*nlmn +ilmn,idat)
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

   if (choice==2.or.choice==23) then
     shift=0
     iatm=0
     force_shift=merge(6,0,choice==23)
     do itypat=1, ntypat
       nlmn=count(indlmn(3,:,itypat)>0)
       nattyp_i = nattyp(itypat)

       !$OMP TARGET TEAMS DISTRIBUTE &
       !$OMP& MAP(to:enlout,gxfac,dgxdt) &
       !$OMP& PRIVATE(idat,igrad,ia,esum)
       do idat=1,ndat*nspinor
         do ia=1,nattyp_i
           do igrad=1,3
             !Following loops are a [D][Z]DOT
             esum=zero
             !$OMP PARALLEL DO REDUCTION(+:esum) PRIVATE(ilmn,ii)
             do ilmn=1,nlmn
               do ii=1,cplex
                 esum=esum +gxfac(ii,shift+(ia-1)*nlmn+ilmn,idat) &
                 &         *dgxdt(ii,grad_shift*shift+(ia-1)*nlmn*grad_shift+(igrad-1+force_shift)*nlmn +ilmn,idat)
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
 end if ! choice=2, 3 or 23

 ! Reduction in case of parallelism
 if (mpi_enreg%paral_spinor==1) then
   if (size(enlout)>0) then
     !$OMP TARGET UPDATE FROM(enlout)
     call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
   end if
   if (choice==3.or.choice==23) then
     !$OMP TARGET UPDATE FROM(enlk)
     call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
   end if
 end if

 ! Derivatives wrt strain
 !  - Convert from reduced to cartesian coordinates
 !  - Substract volume contribution
 if ((choice==3.or.choice==23).and.paw_opt<=3) then
   !$OMP TARGET UPDATE FROM(enlout,enlk)
   do idat=1,ndat
     enlout_shift=(idat-1)*nnlout
     call strconv(enlout(enlout_shift+1:enlout_shift+6),gprimd,work)
     enlout(enlout_shift+1:enlout_shift+3)=(work(1:3)-enlk(idat))
     enlout(enlout_shift+4:enlout_shift+6)= work(4:6)
   end do
 end if

 if (allocated(enlk)) then
   !$OMP TARGET EXIT DATA MAP(delete:enlk,nattyp)
   ABI_FREE(enlk)
 end if
 !$OMP TARGET EXIT DATA MAP(delete:nattyp)

end subroutine opernld_ylm_allwf_ompgpu

end module m_opernld_ylm_allwf_ompgpu

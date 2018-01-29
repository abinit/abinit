!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_sygra
!!
!! NAME
!! dfpt_sygra
!!
!! FUNCTION
!! Symmetrize derivatives of energy with respect to coordinates,
!! as appearing in phonon calculations ...
!! Unsymmetrized gradients are input as deunsy; symmetrized grads are then
!! placed in desym.
!! If nsym=1 simply copy deunsy into desym (only symmetry is identity).
!! The index of the initial perturbation is needed, in case there is a change
!! of atom position (moved in another cell) due to the symmetry operation.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (DCA,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  deunsy(2,3,natom)=unsymmetrized gradients wrt dimensionless tn (hartree)
!!  note : there is a real and a imaginary part ...
!!  indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!   label which gets rotated into given atom by given symmetry
!!   (first three elements are related primitive translation--
!!   see symatm where this is computed)
!!  nsym=number of symmetry operators in group
!!  ipert=index of the initial perturbation
!!  qpt(3)= wavevector of the phonon, in reduced coordinates
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!    reciprocal space primitive translations--see comments below
!!
!! OUTPUT
!! desym(2,3,natom)=symmetrized gradients wrt dimensionless tn (hartree)
!!
!! NOTES
!! Written by X. Gonze starting from sygrad, written by D.C. Allan :
!!    introduction of the q vector for phonon symmetrization
!! This routine should once be merged with sygrad ...
!!
!! PARENTS
!!      dfpt_nstdy,dfpt_nstpaw
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_sygra(natom,desym,deunsy,indsym,ipert,nsym,qpt,symrec)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_sygra'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ipert,natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym)
 real(dp),intent(in) :: deunsy(2,3,natom),qpt(3)
 real(dp),intent(out) :: desym(2,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,ind,isym,mu
 real(dp) :: arg,im,re,sumi,sumr

! *********************************************************************

!DEBUG
!write(std_out,*)' dfpt_sygra : enter '
!write(std_out,*)' dfpt_sygra : qpt(:)',qpt(:)
!do ia=1,natom
!do mu=1,3
!write(std_out,*)' dfpt_sygra : deunsy(:2,mu,ia)',deunsy(:2,mu,ia)
!enddo
!enddo
!ENDDEBUG

 if (nsym==1) then

!  Only symmetry is identity so simply copy
   desym(:,:,:)=deunsy(:,:,:)

 else

!  Actually conduct symmetrization
!  DEBUG
!  write(std_out,*)' dfpt_sygra : desym(:2,:3,:natom),qpt(:)',desym(:2,:3,:natom),qpt(:)
!  ENDDEBUG
   do ia=1,natom
!    DEBUG
!    write(std_out,*)' dfpt_sygra : ia=',ia
!    ENDDEBUG
     do mu=1,3
       sumr=zero
       sumi=zero
!      DEBUG
!      write(std_out,*)' dfpt_sygra : mu=',mu
!      ENDDEBUG
       do isym=1,nsym
         ind=indsym(4,isym,ia)
!        Must shift the atoms back to the unit cell.
!        arg=two_pi*( qpt(1)*indsym(1,isym,ia)&
!        &         +qpt(2)*indsym(2,isym,ia)&
!        &         +qpt(3)*indsym(3,isym,ia) )
!        Selection of non-zero q point, to avoid ipert being outside the 1 ... natom range
         if(qpt(1)**2+qpt(2)**2+qpt(3)**2 > tol16)then
           arg=two_pi*( qpt(1)*(indsym(1,isym,ia)-indsym(1,isym,ipert))&
&           +qpt(2)* (indsym(2,isym,ia)-indsym(2,isym,ipert))&
&           +qpt(3)* (indsym(3,isym,ia)-indsym(3,isym,ipert)))
         else
           arg=zero
         end if

         re=dble(symrec(mu,1,isym))*deunsy(1,1,ind)+&
&         dble(symrec(mu,2,isym))*deunsy(1,2,ind)+&
&         dble(symrec(mu,3,isym))*deunsy(1,3,ind)
         im=dble(symrec(mu,1,isym))*deunsy(2,1,ind)+&
&         dble(symrec(mu,2,isym))*deunsy(2,2,ind)+&
&         dble(symrec(mu,3,isym))*deunsy(2,3,ind)
         sumr=sumr+re*cos(arg)-im*sin(arg)
         sumi=sumi+re*sin(arg)+im*cos(arg)
!        sumr=sumr+re
!        sumi=sumi+im
!        DEBUG
!        write(std_out,*)' dfpt_sygra : isym,indsym(4,isym,ia),arg,re,im,sumr,sumi',&
!        &      isym,indsym(4,isym,ia),arg,re,im,sumr,sumi
!        ENDDEBUG
       end do
       desym(1,mu,ia)=sumr/dble(nsym)
       desym(2,mu,ia)=sumi/dble(nsym)
!      DEBUG
!      write(std_out,*)' dfpt_sygra : desym(:,mu,ia)',desym(:,mu,ia)
!      ENDDEBUG
     end do
   end do
 end if

end subroutine dfpt_sygra
!!***

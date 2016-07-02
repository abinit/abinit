!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_sydy
!!
!! NAME
!! dfpt_sydy
!!
!! FUNCTION
!! Symmetrize dynamical matrix (eventually diagonal wrt to the atoms)
!! Unsymmetrized dynamical matrix   is  input as dyfrow;
!! symmetrized dynamical matrix is then  placed in sdyfro.
!! If nsym=1 simply copy dyfrow   into sdyfro.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=1 if dynamical matrix is real, 2 if it is complex
!!  dyfrow(3,3,natom,1+(natom-1)*nondiag)=unsymmetrized dynamical matrix
!!  indsym(4,msym*natom)=indirect indexing array: for each
!!   isym,iatom, fourth element is label of atom into which iatom is sent by
!!   INVERSE of symmetry operation isym; first three elements are the primitive
!!   translations which must be subtracted after the transformation to get back
!!   to the original unit cell.
!!  natom=number of atoms in cell.
!!  nondiag=0 if dynamical matrix is     diagonal with respect to atoms
!           1 if dynamical matrix is non diagonal with respect to atoms
!!  nsym=number of symmetry operators in group.
!!  qphon(3)=wavevector of the phonon
!!  symq(4,2,nsym)=1 if symmetry preserves present qpoint. From littlegroup_q
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on real
!!    space primitive translations (integers).
!!
!! OUTPUT
!!  sdyfro(3,3,natom,1+(natom-1)*nondiag)=symmetrized dynamical matrix
!!
!! NOTES
!! Symmetrization of gradients with respect to reduced
!! coordinates tn is conducted according to the expression
!! $[d(e)/d(t(n,a))]_{symmetrized} = (1/Nsym)*Sum(S)*symrec(n,m,S)*
!!              [d(e)/d(t(m,b))]_{unsymmetrized}$
!! where $t(m,b)= (symrel^{-1})(m,n)*(t(n,a)-tnons(n))$ and tnons
!! is a possible nonsymmorphic translation.  The label "b" here
!! refers to the atom which gets rotated into "a" under symmetry "S".
!! symrel is the symmetry matrix in real space, which is the inverse
!! transpose of symrec.  symrec is the symmetry matrix in reciprocal
!! space.  $sym_{cartesian} = R * symrel * R^{-1} = G * symrec * G^{-1}$
!! where the columns of R and G are the dimensional primitive translations
!! in real and reciprocal space respectively.
!! Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      dfpt_dyfro
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_sydy(cplex,dyfrow,indsym,natom,nondiag,nsym,qphon,sdyfro,symq,symrec)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_sydy'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,natom,nondiag,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symq(4,2,nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: dyfrow(cplex,3,3,natom,1+(natom-1)*nondiag),qphon(3)
 real(dp),intent(out) :: sdyfro(cplex,3,3,natom,1+(natom-1)*nondiag)

!Local variables -------------------------
!scalars
 integer :: ia,indi,indj,isym,ja,kappa,mu,natom_nondiag,nsym_used,nu
 logical :: qeq0
 real(dp) :: arg,div,phasei,phaser
!arrays
 real(dp) :: work(cplex,3,3)

! *********************************************************************

 if (nsym==1) then

!  Only symmetry is identity so simply copy
   sdyfro(:,:,:,:,:)=dyfrow(:,:,:,:,:)

 else

!  Actually carry out symmetrization
   sdyfro(:,:,:,:,:)=zero
   qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-14)
!  === Diagonal dyn. matrix OR q=0
   if (nondiag==0.or.qeq0) then
     natom_nondiag=1;if (nondiag==1) natom_nondiag=natom
     do ja=1,natom_nondiag
       do ia=1,natom
         do isym=1,nsym
           indi=indsym(4,isym,ia)
           indj=1;if (nondiag==1) indj=indsym(4,isym,ja)
           work(:,:,:)=zero
           do mu=1,3
             do nu=1,3
               do kappa=1,3
                 work(:,mu,kappa)=work(:,mu,kappa)+symrec(mu,nu,isym)*dyfrow(:,nu,kappa,indi,indj)
               end do
             end do
           end do
           do mu=1,3
             do nu=1,3
               do kappa=1,3
                 sdyfro(:,kappa,mu,ia,ja)=sdyfro(:,kappa,mu,ia,ja)+symrec(mu,nu,isym)*work(:,kappa,nu)
               end do
             end do
           end do
         end do
       end do
     end do
     div=one/dble(nsym)
     sdyfro(:,:,:,:,:)=div*sdyfro(:,:,:,:,:)
!    === Non diagonal dyn. matrix AND q<>0
   else
     do ja=1,natom
       do ia=1,natom
         nsym_used=0
         do isym=1,nsym
           if (symq(4,1,isym)==1) then
             arg=two_pi*(qphon(1)*(indsym(1,isym,ia)-indsym(1,isym,ja)) &
&             +qphon(2)*(indsym(2,isym,ia)-indsym(2,isym,ja)) &
&             +qphon(3)*(indsym(3,isym,ia)-indsym(3,isym,ja)))
             phaser=cos(arg);phasei=sin(arg)
             nsym_used=nsym_used+1
             indi=indsym(4,isym,ia)
             indj=indsym(4,isym,ja)
             work(:,:,:)=zero
             do mu=1,3
               do nu=1,3
                 do kappa=1,3
                   work(:,mu,kappa)=work(:,mu,kappa)+symrec(mu,nu,isym)*dyfrow(:,nu,kappa,indi,indj)
                 end do
               end do
             end do
             do mu=1,3
               do nu=1,3
                 do kappa=1,3
                   sdyfro(1,kappa,mu,ia,ja)=sdyfro(1,kappa,mu,ia,ja) &
&                   +symrec(mu,nu,isym)*(work(1,kappa,nu)*phaser-work(2,kappa,nu)*phasei)
                 end do
               end do
             end do
             if (cplex==2) then
               do mu=1,3
                 do nu=1,3
                   do kappa=1,3
                     sdyfro(2,kappa,mu,ia,ja)=sdyfro(2,kappa,mu,ia,ja) &
&                     +symrec(mu,nu,isym)*(work(1,kappa,nu)*phasei+work(2,kappa,nu)*phaser)
                   end do
                 end do
               end do
             end if
           end if
         end do
         div=one/dble(nsym_used)
         sdyfro(:,:,:,ia,ja)=div*sdyfro(:,:,:,ia,ja)
       end do
     end do
   end if

 end if

end subroutine dfpt_sydy
!!***

!    CODE TO BE EVENTUALLY REUSED
!    Sym preserves direction and atom
!    if (symq(1,1,isym)==0.and.symq(2,1,isym)==0.and.symq(3,1,isym)==0.and.symq(4,1,isym)==1)then
!      if (ipert==indsym(4,isym,ipert)) then
!        tok=1
!        do idir1=1,3
!          if ((idir1==idir.and.symrec(idir,idir1,isym)/=1).or.&
! &            (idir1/=idir.and.symrec(idir,idir1,isym)/=0)) tok=0
!        end do
!      end if
!    end if
!    div=one/dble(count(symq(4,1,:)==1))

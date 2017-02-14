!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_dfptnl_accrhoij
!!
!! NAME
!! paw_dfptnl_accrhoij
!!
!! FUNCTION
!! Accumulate the PAW quantities rhoij (augmentation occupancies)
!! or their 1-st order change or their gradient vs r
!! Add the contribution of a given k-point and band
!! Remember: for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (sorted-->random), inverse of atindx.
!!  cplex: if 1, WFs (or 1st-order WFs) are REAL, if 2, COMPLEX
!!  cwaveprj(natom,nspinor)= LEFT wave function at given n,k
!!                         projected with non-local projectors: cwaveprj=<p_i|Cnk>
!!  cwaveprj1(natom,nspinor)= RIGHT wave function at n,k,q
!!                          projected with non-local projectors: cwaveprj1=<p_i|C1nk,q>
!!                          * USED for RF  : C1nk is the first-order wave function
!!                          * USED for DMFT: C1nk is the RIGHT wave function
!!                          * NOT USED in usual GS case; can be set to cwaveprj in that case
!!  ipert=index of perturbation (RF only, i.e. option=2)
!!  isppol=index of current spin component
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  nspinor=number of spinorial components (on current proc)
!!  occ_k=occupation number for current band n,k
!!  option: choice of calculation:
!!          1: update rhoij (Ground-State)
!!          2: update 1-st order rhoij (Response Function) according to ipert
!!          3: update gradients of rhoij with respect to r,strain of both
!!  usetimerev=.TRUE. if time-reversal symmetry is used (WF(-k)=Conjg[WF(k)])
!!  wtk_k=weight assigned to current k-point
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= GS: paw rhoij occupancies and related data
!!                                         RF: 1-st order paw rhoij occupancies and related data
!!  On output, has been updated with the contribution of current n,k
!!    === option=1:
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)}
!!    === option=2:
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[cprj_nk(ii)].cprj1_nk,q(jj)
!!                                 conjugate[cprj_nk(jj)].cprj1_nk,q(ii)}
!!          + Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dlambda].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dlambda)}
!!    === option=3:
!!        pawrhoij(:)%grhoij(lmn2_size,mu,nspden)=   (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dr_mu].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dr_mu)}
!!
!! PARENTS
!!      d2frnl,dfpt_accrho,energy,pawmkrhoij,posdoppler,wfd_pawrhoij
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine paw_dfptnl_accrhoij(atindx,cplex,cwaveprj0_pert1,cwaveprj0_pert2,&
&                       cwaveprj1_pert12,cwaveprj1_pert21,ipert1,ipert2,isppol,my_natom,natom,&
&                       nspinor,occ_k,pawrhoij,usetimerev,wtk_k,occ_k_2, &
&                       comm_atom,mpi_atmtab ) ! optional (parallelism)


 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_pawrhoij,   only : pawrhoij_type
 use m_pawcprj,    only : pawcprj_type
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_dfptnl_accrhoij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ipert1,ipert2,isppol,my_natom,natom,nspinor
 integer,optional,intent(in) :: comm_atom
 logical,intent(in) :: usetimerev
 real(dp),intent(in) :: occ_k,wtk_k
 real(dp),optional,intent(in) :: occ_k_2
!arrays
 integer,intent(in) :: atindx(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawcprj_type),intent(in) :: cwaveprj0_pert1(natom,nspinor),cwaveprj0_pert2(natom,nspinor)
 type(pawcprj_type),intent(in) :: cwaveprj1_pert12(natom,nspinor),cwaveprj1_pert21(natom,nspinor)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iatm,iatom,iatom1,ilmn,iplex,j0lmn,jlmn,klmn,klmn_im,klmn_re
 integer :: mu,my_comm_atom,ncpgr,nspden_rhoij
 logical :: compute_impart,compute_impart_cplex,substract_diagonal
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: ro11_im,ro11_re,ro12_im,ro12_re,ro21_im,ro21_re,ro22_im,ro22_re,weight,weight_2
 character(len=500) :: message
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: cpi0(2,nspinor),d1cpi0(2,nspinor),d2cpi0(2,nspinor)
 real(dp) :: cpj0(2,nspinor),d1cpj0(2,nspinor),d2cpj0(2,nspinor)
 real(dp) :: cpi1(2,nspinor),d2cpi1(2,nspinor)
 real(dp) :: cpj1(2,nspinor),d2cpj1(2,nspinor)
 real(dp) :: cpi2(2,nspinor),d1cpi2(2,nspinor)
 real(dp) :: cpj2(2,nspinor),d1cpj2(2,nspinor)

! ***********************************************************************

 DBG_ENTER("COLL")

 if (my_natom==0) return

 ncpgr=1
 if (ipert1<0.or.ipert1>natom+2.or.ipert2<0.or.ipert2>natom+2) then
   message = 'Wrong inputs. Necessary conditions on ipert1 or ipert2 : 0 <= ipert <= natom+2'
   MSG_BUG(message)
 end if

! ncpgr=0
! if (option==2.and.(ipert<=natom.or.ipert==natom+3.or.ipert==natom+4)) ncpgr=1
! if (option==3) ncpgr=cwaveprj(1,1)%ncpgr
!!Tests
! if(option==2.and.(ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11)) then
!   message = ' not relevant for ipert=natom+1 or ipert=natom+10 or ipert=natom+11 !'
!   MSG_BUG(message)
! end if
! if(option==2.and.cwaveprj(1,1)%ncpgr<ncpgr) then
!   message = ' Error on cwaveprj1 factors derivatives !'
!   MSG_BUG(message)
! end if
! if(option==3.and.cwaveprj(1,1)%ncpgr/=ncpgr) then
!   message = ' Error on cwaveprj factors derivatives !'
!   MSG_BUG(message)
! end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

 weight=wtk_k*occ_k
 if (pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1.and.nspinor==1) weight=half*weight

!!  ==================================================================
!!  === Accumulate (n,k) contribution to partial 2nd-order rhoij   ===
!!  ==================================================================

!compute_impart=(pawrhoij(1)%cplex==2)
 compute_impart_cplex=((pawrhoij(1)%cplex==2).and.(cplex==2))
!substract_diagonal=(ipert==natom+3)

!Accumulate :   < Psi^(pert1) | p_i^(0) > < p_j^(0) | Psi^(pert2) >
!             + < Psi^(pert2) | p_i^(0) > < p_j^(0) | Psi^(pert1) >
 if (nspinor==1) then
   do iatom=1,my_natom
     iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
     iatm=atindx(iatom1)
     cplex_rhoij=pawrhoij(iatom)%cplex
     do jlmn=1,pawrhoij(iatom)%lmn_size
       j0lmn=jlmn*(jlmn-1)/2
       cpj1(1:2,1)=cwaveprj1_pert12(iatm,1)%cp(1:2,jlmn)   ! < p_j^(0) | Psi^(pert1) >
       cpj2(1:2,1)=cwaveprj1_pert21(iatm,1)%cp(1:2,jlmn)   ! < p_j^(0) | Psi^(pert2) >
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         klmn_re=cplex_rhoij*(klmn-1)+1
         cpi1(1:2,1)=cwaveprj1_pert12(iatm,1)%cp(1:2,ilmn) ! < p_i^(0) | Psi^(pert1) >
         cpi2(1:2,1)=cwaveprj1_pert21(iatm,1)%cp(1:2,ilmn) ! < p_i^(0) | Psi^(pert2) >
         ro11_re=zero
         do iplex=1,cplex
           ro11_re=ro11_re+cpi1(iplex,1)*cpj2(iplex,1)+cpj1(iplex,1)*cpi2(iplex,1)
         end do
         pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
         if (compute_impart_cplex) then
           klmn_im=klmn_re+1
           ro11_im=        cpi1(1,1)*cpj2(2,1)-cpi1(2,1)*cpj2(1,1)
           ro11_im=ro11_im+cpj1(1,1)*cpi2(2,1)-cpj1(2,1)*cpi2(1,1)
           pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
         end if
       end do
     end do
   end do
 else ! nspinor=2
   MSG_BUG("paw_dfptnl_accrhoij is not implemented for nspinor=2")
 end if

!Accumulate :   < Psi^(pert1) | p_i^(pert2) > < p_j^(0)     | Psi^(0)     >
!             + < Psi^(pert1) | p_i^(0)     > < p_j^(pert2) | Psi^(0)     >
!             + < Psi^(0)     | p_i^(pert2) > < p_j^(0)     | Psi^(pert1) >
!             + < Psi^(0)     | p_i^(0)     > < p_j^(pert2) | Psi^(pert1) >
 if (ipert2>0.and.ipert2<=natom) then
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)  =cwaveprj0_pert2 (iatm,1)% cp(1:2,  jlmn) ! < p_j^(0)     | Psi^(0)     >
         d2cpj0(1:2,1)=cwaveprj0_pert2 (iatm,1)%dcp(1:2,1,jlmn) ! < p_j^(pert2) | Psi^(0)     >
         cpj1(1:2,1)  =cwaveprj1_pert12(iatm,1)% cp(1:2,  jlmn) ! < p_j^(0)     | Psi^(pert1) >
         d2cpj1(1:2,1)=cwaveprj1_pert12(iatm,1)%dcp(1:2,1,jlmn) ! < p_j^(pert2) | Psi^(pert1) >
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)  =cwaveprj0_pert2 (iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(0)     >
           d2cpi0(1:2,1)=cwaveprj0_pert2 (iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert2) | Psi^(0)     >
           cpi1(1:2,1)  =cwaveprj1_pert12(iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(pert1) >
           d2cpi1(1:2,1)=cwaveprj1_pert12(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert2) | Psi^(pert1) >
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+d2cpi1(iplex,1)*  cpj0(iplex,1)
             ro11_re=ro11_re+  cpi1(iplex,1)*d2cpj0(iplex,1)
             ro11_re=ro11_re+d2cpi0(iplex,1)*  cpj1(iplex,1)
             ro11_re=ro11_re+  cpi0(iplex,1)*d2cpj1(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=        d2cpi1(1,1)*  cpj0(2,1)-d2cpi1(2,1)*  cpj0(1,1)
             ro11_im=ro11_im+  cpi1(1,1)*d2cpj0(2,1)-  cpi1(2,1)*d2cpj0(1,1)
             ro11_im=ro11_im+d2cpi0(1,1)*  cpj1(2,1)-d2cpi0(2,1)*  cpj1(1,1)
             ro11_im=ro11_im+  cpi0(1,1)*d2cpj1(2,1)-  cpi0(2,1)*d2cpj1(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     MSG_BUG("paw_dfptnl_accrhoij is not implemented for nspinor=2")
   end if
 end if

!Accumulate : < Psi^(0)     | p_i^(pert1) > < p_j^(0)     | Psi^(pert2) >
!           + < Psi^(pert2) | p_i^(pert1) > < p_j^(0)     | Psi^(0)     >
!           + < Psi^(0)     | p_i^(pert1) > < p_j^(0)     | Psi^(pert2) >
!           + < Psi^(pert2) | p_i^(pert1) > < p_j^(0)     | Psi^(0)     >
 if (ipert1>0.and.ipert1<=natom) then
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)  =cwaveprj0_pert1 (iatm,1)% cp(1:2,  jlmn) ! < p_j^(0)     | Psi^(0)     >
         d1cpj0(1:2,1)=cwaveprj0_pert1 (iatm,1)%dcp(1:2,1,jlmn) ! < p_j^(pert1) | Psi^(0)     >
         cpj2(1:2,1)  =cwaveprj1_pert21(iatm,1)% cp(1:2,  jlmn) ! < p_j^(0)     | Psi^(pert2) >
         d1cpj2(1:2,1)=cwaveprj1_pert21(iatm,1)%dcp(1:2,1,jlmn) ! < p_j^(pert1) | Psi^(pert2) >
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)  =cwaveprj0_pert2 (iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(0)     >
           d1cpi0(1:2,1)=cwaveprj0_pert2 (iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(0)     >
           cpi2(1:2,1)  =cwaveprj1_pert12(iatm,1)% cp(1:2,  ilmn) ! < p_i^(0)     | Psi^(pert2) >
           d1cpi2(1:2,1)=cwaveprj1_pert12(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(pert2) >
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+d1cpi2(iplex,1)*  cpj0(iplex,1)
             ro11_re=ro11_re+  cpi2(iplex,1)*d1cpj0(iplex,1)
             ro11_re=ro11_re+d1cpi0(iplex,1)*  cpj2(iplex,1)
             ro11_re=ro11_re+  cpi0(iplex,1)*d1cpj2(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=        d1cpi2(1,1)*  cpj0(2,1)-d1cpi2(2,1)*  cpj0(1,1)
             ro11_im=ro11_im+  cpi2(1,1)*d1cpj0(2,1)-  cpi2(2,1)*d1cpj0(1,1)
             ro11_im=ro11_im+d1cpi0(1,1)*  cpj2(2,1)-d1cpi0(2,1)*  cpj2(1,1)
             ro11_im=ro11_im+  cpi0(1,1)*d1cpj2(2,1)-  cpi0(2,1)*d1cpj2(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     MSG_BUG("paw_dfptnl_accrhoij is not implemented for nspinor=2")
   end if
 end if
!  End

!Accumulate :   < Psi^(0) | p_i^(pert1) > < p_j^(pert2) | Psi^(0) >
!             + < Psi^(0) | p_i^(pert2) > < p_j^(pert1) | Psi^(0) >
 if (ipert1>0.and.ipert1<=natom.and.ipert2>0.and.ipert2<=natom) then
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         d1cpj0(1:2,1)=cwaveprj0_pert1(iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert1) | Psi^(0) >
         d2cpj0(1:2,1)=cwaveprj0_pert2(iatm,1)%dcp(1:2,1,jlmn)   ! < p_j^(pert2) | Psi^(0) >
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           d1cpi0(1:2,1)=cwaveprj0_pert1(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(0) >
           d2cpi0(1:2,1)=cwaveprj0_pert2(iatm,1)%dcp(1:2,1,ilmn) ! < p_i^(pert1) | Psi^(0) >
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+d1cpi0(iplex,1)*d2cpj0(iplex,1)+d2cpi0(iplex,1)*d1cpj0(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=        d1cpi0(1,1)*d2cpj0(2,1)-d1cpi0(2,1)*d2cpj0(1,1)
             ro11_im=ro11_im+d2cpi0(1,1)*d1cpj0(2,1)-d2cpi0(2,1)*d1cpj0(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     MSG_BUG("paw_dfptnl_accrhoij is not implemented for nspinor=2")
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine paw_dfptnl_accrhoij
!!***

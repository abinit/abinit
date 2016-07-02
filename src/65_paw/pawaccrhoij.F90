!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawaccrhoij
!!
!! NAME
!! pawaccrhoij
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

 subroutine pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj1,ipert,isppol,my_natom,natom,&
&                       nspinor,occ_k,option,pawrhoij,usetimerev,wtk_k,occ_k_2, &
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
#define ABI_FUNC 'pawaccrhoij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ipert,isppol,my_natom,natom,nspinor,option
 integer,optional,intent(in) :: comm_atom
 logical,intent(in) :: usetimerev
 real(dp),intent(in) :: occ_k,wtk_k
 real(dp),optional,intent(in) :: occ_k_2
!arrays
 integer,intent(in) :: atindx(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawcprj_type),intent(in) :: cwaveprj(natom,nspinor),cwaveprj1(natom,nspinor)
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
 real(dp) :: cpi0(2,nspinor),cpi1(2,nspinor),cpj0(2,nspinor),cpj1(2,nspinor)
 real(dp) :: dcpi0(2,nspinor,9),dcpj0(2,nspinor,9)

! ***********************************************************************

 DBG_ENTER("COLL")

 if (my_natom==0) return

 ncpgr=0
 if (option==2.and.(ipert<=natom.or.ipert==natom+3.or.ipert==natom+4)) ncpgr=1
 if (option==3) ncpgr=cwaveprj(1,1)%ncpgr
!Tests
 if(option==2.and.(ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11)) then
   message = ' not relevant for ipert=natom+1 or ipert=natom+10 or ipert=natom+11 !'
   MSG_BUG(message)
 end if
 if(option==2.and.cwaveprj(1,1)%ncpgr<ncpgr) then
   message = ' Error on cwaveprj1 factors derivatives !'
   MSG_BUG(message)
 end if
 if(option==3.and.cwaveprj(1,1)%ncpgr/=ncpgr) then
   message = ' Error on cwaveprj factors derivatives !'
   MSG_BUG(message)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

 weight=wtk_k*occ_k
 weight_2=zero
 if(present(occ_k_2).and.nspinor==2) weight_2=wtk_k*occ_k_2
 if (pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1.and.nspinor==1) weight=half*weight
 if (pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1.and.nspinor==1.and.present(occ_k_2)) weight_2=half*weight_2
 if (option==1) then

!  ==================================================================
!  === OPTION 1: Accumulate (n,k) contribution to rhoij =============
!  ==================================================================
   compute_impart=((.not.usetimerev).and.(pawrhoij(1)%cplex==2))
   compute_impart_cplex=((compute_impart).and.(cplex==2))
   if (nspinor==1) then
     cplex_rhoij=pawrhoij(1)%cplex
     if (cplex_rhoij==1) then
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             cpi0(1:cplex,1)=cwaveprj1(iatm,1)%cp(1:cplex,ilmn)
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+cpi0(iplex,1)*cpj0(iplex,1)
             end do
             pawrhoij(iatom)%rhoij_(klmn,isppol)=pawrhoij(iatom)%rhoij_(klmn,isppol)+weight*ro11_re
           end do
         end do
       end do
     else
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             klmn_re=cplex_rhoij*(klmn-1)+1
             cpi0(1:cplex,1)=cwaveprj1(iatm,1)%cp(1:cplex,ilmn)
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+cpi0(iplex,1)*cpj0(iplex,1)
             end do
             pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
             if (compute_impart_cplex) then
               klmn_im=klmn_re+1
               ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
               pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
             end if
           end do
         end do
       end do
     end if
   else ! nspinor=2
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       nspden_rhoij=pawrhoij(iatom)%nspden
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
         cpj0(1:cplex,2)=cwaveprj(iatm,2)%cp(1:cplex,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:cplex,1)=cwaveprj1(iatm,1)%cp(1:cplex,ilmn)
           cpi0(1:cplex,2)=cwaveprj1(iatm,2)%cp(1:cplex,ilmn)
           ro11_re=zero;ro22_re=zero
           ro12_re=zero;ro21_re=zero
           ro12_im=zero;ro21_im=zero
           do iplex=1,cplex
             ro11_re=ro11_re+cpi0(iplex,1)*cpj0(iplex,1)
             ro22_re=ro22_re+cpi0(iplex,2)*cpj0(iplex,2)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
           if (nspden_rhoij>1) then
             do iplex=1,cplex
               ro12_re=ro12_re+cpi0(iplex,2)*cpj0(iplex,1)
               ro21_re=ro21_re+cpi0(iplex,1)*cpj0(iplex,2)
             end do
             pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
             pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
             if (cplex==2) then
               ro12_im=cpi0(1,2)*cpj0(2,1)-cpi0(2,2)*cpj0(1,1)
               ro21_im=cpi0(1,1)*cpj0(2,2)-cpi0(2,1)*cpj0(1,2)
               pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
             end if
           end if
           if (compute_impart) then
             klmn_im=klmn_re+1
             if (nspden_rhoij>1) pawrhoij(iatom)%rhoij_(klmn_im,3)=pawrhoij(iatom)%rhoij_(klmn_im,3)+weight*(ro12_re-ro21_re)
             if (cplex==2) then
               ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
               ro22_im=cpi0(1,2)*cpj0(2,2)-cpi0(2,2)*cpj0(1,2)
               pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
               if (nspden_rhoij>1) then
                 pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_im,2)+weight*(ro12_im+ro21_im)
               end if
               if (present(occ_k_2).and.nspinor==2) then
                 pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight_2*(-ro11_im-ro22_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight_2*( ro11_re+ro22_re)
                 pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight_2*(-ro21_im-ro12_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_im,2)+weight_2*( ro21_re+ro12_re)
                 pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight_2*(-ro12_re+ro21_re)
                 pawrhoij(iatom)%rhoij_(klmn_im,3)=pawrhoij(iatom)%rhoij_(klmn_im,3)+weight_2*(-ro12_im+ro21_im)
                 pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight_2*(-ro11_im+ro22_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight_2*( ro11_re-ro22_re)
               end if
             end if
           end if
         end do
       end do
     end do
   end if

 else if (option==2) then

!  ==================================================================
!  === OPTION 2: Accumulate (n,k) contribution to 1st-order rhoij ===
!  ==================================================================

   compute_impart=(pawrhoij(1)%cplex==2)
   compute_impart_cplex=((pawrhoij(1)%cplex==2).and.(cplex==2))
   substract_diagonal=(ipert==natom+3)

!  Accumulate (n,k) contribution to rhoij1
!  due to derivative of wave-function
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,jlmn)
         cpj1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,ilmn)
           cpi1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,ilmn)
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+cpi0(iplex,1)*cpj1(iplex,1)+cpj0(iplex,1)*cpi1(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             ro11_im=cpi0(1,1)*cpj1(2,1)-cpi0(2,1)*cpj1(1,1)+cpj0(1,1)*cpi1(2,1)-cpj0(2,1)*cpi1(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       nspden_rhoij=pawrhoij(iatom)%nspden
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,jlmn)
         cpj0(1:2,2)=cwaveprj (iatm,2)%cp(1:2,jlmn)
         cpj1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,jlmn)
         cpj1(1:2,2)=cwaveprj1(iatm,2)%cp(1:2,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,ilmn)
           cpi0(1:2,2)=cwaveprj (iatm,2)%cp(1:2,ilmn)
           cpi1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,ilmn)
           cpi1(1:2,2)=cwaveprj1(iatm,2)%cp(1:2,ilmn)
           ro11_re=zero;ro22_re=zero
           ro12_re=zero;ro21_re=zero
           ro12_im=zero;ro21_im=zero
           do iplex=1,cplex
             ro11_re=cpj0(iplex,1)*cpi1(iplex,1)+cpi0(iplex,1)*cpj1(iplex,1)
             ro22_re=cpj0(iplex,2)*cpi1(iplex,2)+cpi0(iplex,2)*cpj1(iplex,2)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
           if (nspden_rhoij>1) then
             do iplex=1,cplex
               ro12_re=cpj0(iplex,1)*cpi1(iplex,2)+cpi0(iplex,2)*cpj1(iplex,1)
               ro21_re=cpj0(iplex,2)*cpi1(iplex,1)+cpi0(iplex,1)*cpj1(iplex,2)
             end do
             pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
             pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
             if (cplex==2) then
               ro12_im=cpj0(1,1)*cpi1(2,2)-cpi1(1,1)*cpj0(2,2)+cpi0(1,2)*cpj1(2,1)-cpj1(1,2)*cpi0(2,1)
               ro21_im=cpj0(1,2)*cpi1(2,1)-cpi1(1,2)*cpj0(2,1)+cpi0(1,1)*cpj1(2,2)-cpj1(1,1)*cpi0(2,2)
               pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
             end if
           end if
           if (compute_impart) then
             klmn_im=klmn_re+1
             if (nspden_rhoij>1) pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro12_re-ro21_re)
             if (cplex==2) then
               ro11_im=cpj0(1,1)*cpi1(2,1)-cpi1(1,1)*cpj0(2,1)+cpi0(1,1)*cpj1(2,1)-cpj1(1,1)*cpi0(2,1)
               ro22_im=cpj0(1,2)*cpi1(2,2)-cpi1(1,2)*cpj0(2,2)+cpi0(1,2)*cpj1(2,2)-cpj1(1,2)*cpi0(2,2)
               pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
               if (nspden_rhoij>1) then
                 pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
                 pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_im+ro21_im)
               end if
             end if
           end if
         end do
       end do
     end do
   end if

!  Accumulate (n,k) contribution to rhoij1
!  due to derivative of projectors
   if (ipert/=natom+2) then
     if (nspinor==1) then
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         if (ipert<=natom.and.iatom/=ipert) cycle
         cplex_rhoij=pawrhoij(iatom)%cplex
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,jlmn)
           dcpj0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             klmn_re=cplex_rhoij*(klmn-1)+1
             cpi0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,ilmn)
             dcpi0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,ilmn)
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dcpi0(iplex,1,1)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,1)
             end do
             if (substract_diagonal) then
               do iplex=1,cplex
                 ro11_re=ro11_re-cpi0(iplex,1)*cpj0(iplex,1)
               end do
             end if
             pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
!            This imaginary part does not have to be computed
!            It is cancelled because rho_ij+rho_ji is stored in rho_ij
!            if (compute_impart_cplex) then
!            klmn_im=klmn_re+1
!            ro11_im=dcpi0(1,1,1)*cpj0(2,1)-dcpi0(2,1,1)*cpj0(1,1)+cpi0(1,1)*dcpj0(2,1,1)-cpi0(2,1)*dcpj0(1,1,1)
!            if (substract_diagonal) then
!            ro11_im=ro11_im-cpi0(1,1)*cpj0(2,1)+cpi0(2,1)*cpj0(1,1)
!            end if
!            pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
!            end if
           end do
         end do
       end do
     else ! nspinor=2
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         if (ipert<=natom.and.iatom/=ipert) cycle
         cplex_rhoij=pawrhoij(iatom)%cplex
         nspden_rhoij=pawrhoij(iatom)%nspden
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,jlmn)
           dcpj0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,jlmn)
           cpj0 (1:2,2)  =cwaveprj(iatm,2)%cp (1:2  ,jlmn)
           dcpj0(1:2,2,1)=cwaveprj(iatm,2)%dcp(1:2,1,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             klmn_re=cplex_rhoij*(klmn-1)+1
             cpi0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,ilmn)
             dcpi0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,ilmn)
             cpi0 (1:2,2)  =cwaveprj(iatm,2)%cp (1:2  ,ilmn)
             dcpi0(1:2,2,1)=cwaveprj(iatm,2)%dcp(1:2,1,ilmn)
             ro11_re=zero;ro22_re=zero
             ro12_re=zero;ro21_re=zero
             ro12_im=zero;ro21_im=zero
             do iplex=1,cplex
               ro11_re=dcpi0(iplex,1,1)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,1)
               ro22_re=dcpi0(iplex,2,1)*cpj0(iplex,2)+cpi0(iplex,2)*dcpj0(iplex,2,1)
             end do
             if (substract_diagonal) then
               do iplex=1,cplex
                 ro11_re=ro11_re-cpi0(iplex,1)*cpj0(iplex,1)
                 ro22_re=ro22_re-cpi0(iplex,2)*cpj0(iplex,2)
               end do
             end if
             pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
             if (nspden_rhoij>1) then
               do iplex=1,cplex
                 ro12_re=dcpi0(iplex,2,1)*cpj0(iplex,1)+cpi0(iplex,2)*dcpj0(iplex,1,1)
                 ro21_re=dcpi0(iplex,1,1)*cpj0(iplex,2)+cpi0(iplex,1)*dcpj0(iplex,2,1)
               end do
               if (substract_diagonal) then
                 do iplex=1,cplex
                   ro12_re=ro12_re-cpi0(iplex,2)*cpj0(iplex,1)
                   ro21_re=ro21_re-cpi0(iplex,1)*cpj0(iplex,2)
                 end do
               end if
               pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
               pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
               if (cplex==2) then
                 ro12_im=dcpi0(1,2,1)*cpj0(2,1)-dcpi0(2,2,1)*cpj0(1,1)+cpi0(1,2)*dcpj0(2,1,1)-cpi0(2,2)*dcpj0(1,1,1)
                 ro21_im=dcpi0(1,1,1)*cpj0(2,2)-dcpi0(2,1,1)*cpj0(1,2)+cpi0(1,1)*dcpj0(2,2,1)-cpi0(2,1)*dcpj0(1,2,1)
                 if (substract_diagonal) then
                   ro12_im=ro12_im-cpi0(1,2)*cpj0(2,1)+cpi0(2,2)*cpj0(1,1)
                   ro21_im=ro21_im-cpi0(1,1)*cpj0(2,2)+cpi0(2,1)*cpj0(1,2)
                 end if
                 pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
               end if
             end if
             if (compute_impart) then
               klmn_im=klmn_re+1
               if (nspden_rhoij>1) pawrhoij(iatom)%rhoij_(klmn_im,3)=pawrhoij(iatom)%rhoij_(klmn_im,3)+weight*(ro12_re-ro21_re)
               if (cplex==2) then
                 ro11_im=dcpi0(1,1,1)*cpj0(2,1)-dcpi0(2,1,1)*cpj0(1,1)+cpi0(1,1)*dcpj0(2,1,1)-cpi0(2,1)*dcpj0(1,1,1)
                 ro22_im=dcpi0(1,2,1)*cpj0(2,2)-dcpi0(2,2,1)*cpj0(1,2)+cpi0(1,2)*dcpj0(2,2,1)-cpi0(2,2)*dcpj0(1,2,1)
                 if (substract_diagonal) then
                   ro11_im=ro11_im-cpi0(1,1)*cpj0(2,1)+cpi0(2,1)*cpj0(1,1)
                   ro22_im=ro22_im-cpi0(1,2)*cpj0(2,2)+cpi0(2,2)*cpj0(1,2)
                 end if
                 pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
                 if (nspden_rhoij>1) then
                   pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
                   pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_im,2)+weight*(ro12_im+ro21_im)
                 end if
               end if
             end if
           end do
         end do
       end do
     end if
   end if

 else if (option==3) then

!  ==================================================================
!  === OPTION 3: Accumulate (n,k) contribution to drhoij/dr =========
!  ==================================================================

   compute_impart=((.not.usetimerev).and.(pawrhoij(1)%cplex==2))
   compute_impart_cplex=((compute_impart).and.(cplex==2))
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:cplex,1)         =cwaveprj(iatm,1)%cp (1:cplex,jlmn)
         dcpj0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:cplex,1)         =cwaveprj(iatm,1)%cp (1:cplex,ilmn)
           dcpi0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,ilmn)
           do mu=1,ncpgr
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dcpi0(iplex,1,mu)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,mu)
             end do
             pawrhoij(iatom)%grhoij(mu,klmn_re,isppol)=pawrhoij(iatom)%grhoij(mu,klmn_re,isppol)+weight*ro11_re
           end do
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             do mu=1,ncpgr
               ro11_im=dcpi0(1,1,mu)*cpj0(2,1)+cpi0(1,1)*dcpj0(2,1,mu)-dcpi0(2,1,mu)*cpj0(1,1)-cpi0(2,1)*dcpj0(1,1,mu)
               pawrhoij(iatom)%grhoij(mu,klmn_im,isppol)=pawrhoij(iatom)%grhoij(mu,klmn_im,isppol)+weight*ro11_im
             end do
           end if
         end do
       end do
     end do
   else ! nspinor=2
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex
       nspden_rhoij=pawrhoij(iatom)%nspden
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,jlmn)
         cpj0(1:cplex,2)     =cwaveprj(iatm,2)%cp (1:cplex,jlmn)
         dcpj0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,jlmn)
         dcpj0(1:cplex,2,1:ncpgr)=cwaveprj(iatm,2)%dcp(1:cplex,1:ncpgr,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           klmn_im=klmn_re+1
           cpi0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,ilmn)
           cpi0(1:cplex,2)     =cwaveprj(iatm,2)%cp (1:cplex,ilmn)
           dcpi0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,ilmn)
           dcpi0(1:cplex,2,1:ncpgr)=cwaveprj(iatm,2)%dcp(1:cplex,1:ncpgr,ilmn)
           do mu=1,ncpgr
             ro11_re=zero;ro22_re=zero
             ro12_re=zero;ro21_re=zero
             ro12_im=zero;ro21_im=zero
             do iplex=1,cplex
               ro11_re=dcpi0(iplex,1,mu)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,mu)
               ro22_re=dcpi0(iplex,2,mu)*cpj0(iplex,2)+cpi0(iplex,2)*dcpj0(iplex,2,mu)
             end do
             pawrhoij(iatom)%grhoij(mu,klmn_re,1)=pawrhoij(iatom)%grhoij(mu,klmn_re,1)+weight*(ro11_re+ro22_re)
             if (nspden_rhoij>1) then
               do iplex=1,cplex
                 ro12_re=ro12_re+dcpi0(iplex,2,mu)*cpj0(iplex,1)+cpi0(iplex,2)*dcpj0(iplex,1,mu)
                 ro21_re=ro21_re+dcpi0(iplex,1,mu)*cpj0(iplex,2)+cpi0(iplex,1)*dcpj0(iplex,2,mu)
               end do
               pawrhoij(iatom)%grhoij(mu,klmn_re,4)=pawrhoij(iatom)%grhoij(mu,klmn_re,4)+weight*(ro11_re-ro22_re)
               pawrhoij(iatom)%grhoij(mu,klmn_re,2)=pawrhoij(iatom)%grhoij(mu,klmn_re,2)+weight*(ro12_re+ro21_re)
               if (cplex==2) then
                 ro12_im=dcpi0(1,2,mu)*cpj0(2,1)+cpi0(1,2)*dcpj0(2,1,mu)-dcpi0(2,2,mu)*cpj0(1,1)-cpi0(2,2)*dcpj0(1,1,mu)
                 ro21_im=dcpi0(1,1,mu)*cpj0(2,2)+cpi0(1,1)*dcpj0(2,2,mu)-dcpi0(2,1,mu)*cpj0(1,2)-cpi0(2,1)*dcpj0(1,2,mu)
                 pawrhoij(iatom)%grhoij(mu,klmn_re,3)=pawrhoij(iatom)%grhoij(mu,klmn_re,3)+weight*(ro21_im-ro12_im)
               end if
             end if
             if (compute_impart) then
               if (nspden_rhoij>1) then
                 pawrhoij(iatom)%grhoij(mu,klmn_im,3)=pawrhoij(iatom)%grhoij(mu,klmn_im,3)+weight*(ro12_re-ro21_re)
               end if
               if (cplex==2) then
                 ro11_im=dcpi0(1,1,mu)*cpj0(2,1)+cpi0(1,1)*dcpj0(2,1,mu)-dcpi0(2,1,mu)*cpj0(1,1)-cpi0(2,1)*dcpj0(1,1,mu)
                 ro22_im=dcpi0(1,2,mu)*cpj0(2,2)+cpi0(1,2)*dcpj0(2,2,mu)-dcpi0(2,2,mu)*cpj0(1,2)-cpi0(2,2)*dcpj0(1,2,mu)
                 pawrhoij(iatom)%grhoij(mu,klmn_im,1)=pawrhoij(iatom)%grhoij(mu,klmn_im,1)+weight*(ro11_im+ro22_im)
                 if (nspden_rhoij>1) then
                   pawrhoij(iatom)%grhoij(mu,klmn_im,4)=pawrhoij(iatom)%grhoij(mu,klmn_im,4)+weight*(ro11_im-ro22_im)
                   pawrhoij(iatom)%grhoij(mu,klmn_im,2)=pawrhoij(iatom)%grhoij(mu,klmn_im,2)+weight*(ro12_im+ro21_im)
                 end if
               end if
             end if
           end do
         end do
       end do
     end do
   end if

!  End
 end if ! option

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawaccrhoij
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawnhatfr
!!
!! NAME
!! pawnhatfr
!!
!! FUNCTION
!! PAW: Compute frozen part of charge compensation density nhat
!!      nhatfr(r)=Sum_ij,lm[rhoij_ij.q_ij^l.(g_l(r).Y_lm(r))^(1)]
!!      Depends on q wave vector but not on first-order wave-function.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  ider=0: computes frozen part of compensation density
!!       1: computes frozen part of compensation density and cartesian gradients
!!  idir=direction of atomic displacement (in case of phonons perturb.)
!!  ipert=nindex of perturbation
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= Ground-State paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translations for real space
!!
!! OUTPUT
!!  pawfgrtab(iatom)%nhatfr(nfgd,nspden)
!!                  frozen part of charge compensation density (inside PAW spheres)
!!                  =Sum_ij,lm[rhoij_ij.q_ij^l.(g_l(r).Y_lm(r))^(1)]
!!  === If ider==1
!!  pawfgrtab(iatom)%nhatfrgr(3,nfgd,nspden)
!!                  gradients of frozen part of charge compensation density (inside PAW spheres)
!!                  =Sum_ij,lm[rhoij_ij.q_ij^l . d/dr((g_l(r).Y_lm(r))^(1))]
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_scfcv,pawmknhat
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawgylm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawnhatfr(ider,idir,ipert,my_natom,natom,nspden,ntypat,&
&                    pawang,pawfgrtab,pawrhoij,pawtab,rprimd, &
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_pawrhoij,     only : pawrhoij_type
 use m_paw_finegrid, only : pawgylm
 use m_paral_atom,   only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawnhatfr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ider,idir,ipert,my_natom,natom,nspden,ntypat
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: rprimd(3,3)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,iatom_tot,ic,ils,ilslm,irhoij,isel,ispden,istr,itypat,jrhoij
 integer :: klm,klmn,lm_size,lmn2_size,lm0,lmax,lmin,mua,mub,mm,mu,my_comm_atom,nfgd,nu,optgr0,optgr1,optgr2
 logical :: my_atmtab_allocated,my_pert,paral_atom
 real(dp) :: contrib,ro
!arrays
 integer,parameter :: voigt(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: nhatfr_tmp(:,:),nhatfrgr_tmp(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Only relevant for atomic displacement and strain perturbation
 if (ipert>natom.and.ipert/=natom+3.and.ipert/=natom+4) return

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Compatibility tests
 if (my_natom>0) then
   if ((pawfgrtab(1)%gylm_allocated==0.or.pawfgrtab(1)%gylmgr_allocated==0).and. &
&   pawfgrtab(1)%rfgd_allocated==0) then
     MSG_BUG('  pawfgrtab()%rfgd array must be allocated  !')
   end if
 end if

 my_pert = (ipert<=natom).or.ipert==natom+3.or.ipert==natom+4

!Get correct index of strain pertubation
 if (ipert==natom+3) istr = idir
 if (ipert==natom+4) istr = idir + 3
 
!Loops over  atoms
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

!  Eventually allocate frozen nhat points
   if (my_pert) then
     if (pawfgrtab(iatom)%nhatfr_allocated==0) then
       if (allocated(pawfgrtab(iatom)%nhatfr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%nhatfr,(pawfgrtab(iatom)%nfgd,nspden))
       pawfgrtab(iatom)%nhatfr_allocated=1
     end if
     if (ider==1.and.pawfgrtab(iatom)%nhatfrgr_allocated==0) then
       if (allocated(pawfgrtab(iatom)%nhatfrgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfrgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%nhatfrgr,(3,pawfgrtab(iatom)%nfgd,nspden))
       pawfgrtab(iatom)%nhatfrgr_allocated=1
     end if
   end if

!  Select if frozen part of nhat exists for the current perturbation
   if ((.not.my_pert).or.(pawfgrtab(iatom)%nhatfr_allocated==0)) cycle

!  Some atom-dependent quantities
   itypat=pawfgrtab(iatom)%itypat
   lm_size=pawfgrtab(iatom)%l_size**2
   lmn2_size=pawtab(itypat)%lmn2_size

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   nfgd=pawfgrtab(iatom)%nfgd
   if ((pawfgrtab(iatom)%gylmgr_allocated==0).or. &
&   (pawfgrtab(iatom)%gylmgr2_allocated==0.and.ider==1)) then
     optgr0=0;optgr1=0;optgr2=0
     if(ipert==natom+3.or.ipert==natom+4)then
       if (pawfgrtab(iatom)%gylm_allocated==0) then
         if (allocated(pawfgrtab(iatom)%gylm))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
         pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
       end if
     end if
     if (pawfgrtab(iatom)%gylmgr_allocated==0) then
       if (allocated(pawfgrtab(iatom)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,lm_size))
       pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
     end if
     if (ider==1.and.pawfgrtab(iatom)%gylmgr2_allocated==0) then
       if (allocated(pawfgrtab(iatom)%gylmgr2))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(6,nfgd,lm_size))
       pawfgrtab(iatom)%gylmgr2_allocated=2;optgr2=1
     end if
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),pawfgrtab(iatom)%rfgd)
     end if
   end if


!  ============ Phonons ====================================
   if (ipert<=natom) then

!    Loop over spin components
     do ispden=1,nspden

       ABI_ALLOCATE(nhatfr_tmp,(3,nfgd))
       nhatfr_tmp=zero
       if (ider==1) then
         ABI_ALLOCATE(nhatfrgr_tmp,(3,nfgd,3))
         nhatfrgr_tmp=zero
       end if
       
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         klm =pawtab(itypat)%indklmn(1,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         if (nspden/=2) then
           ro=pawrhoij(iatom)%rhoijp(jrhoij,ispden)
         else
           if (ispden==1) then
             ro=pawrhoij(iatom)%rhoijp(jrhoij,1)+pawrhoij(iatom)%rhoijp(jrhoij,2)
           else if (ispden==2) then
             ro=pawrhoij(iatom)%rhoijp(jrhoij,1)
           end if
         end if
         ro=pawtab(itypat)%dltij(klmn)*ro
         do ils=lmin,lmax,2
           lm0=ils**2+ils+1
           do mm=-ils,ils
             ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
             if (isel>0) then
               do ic=1,nfgd
                 do mu=1,3
                   contrib=-ro*pawtab(itypat)%qijl(ilslm,klmn)&
&                   *pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                   nhatfr_tmp(mu,ic)=nhatfr_tmp(mu,ic)+contrib
                 end do
               end do
               if (ider==1) then
                 do ic=1,nfgd
                   do nu=1,3
                     do mu=1,3
                       contrib=-ro*pawtab(itypat)%qijl(ilslm,klmn) &
&                       *pawfgrtab(iatom)%gylmgr2(voigt(mu,nu),ic,ilslm)
                       nhatfrgr_tmp(mu,ic,nu)=nhatfrgr_tmp(mu,ic,nu)+contrib
                     end do
                   end do
                 end do
               end if
             end if
           end do
         end do
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do

!      Convert from cartesian to reduced coordinates
       do ic=1,nfgd
         pawfgrtab(iatom)%nhatfr(ic,ispden)= &
&         rprimd(1,idir)*nhatfr_tmp(1,ic) &
&         +rprimd(2,idir)*nhatfr_tmp(2,ic) &
&         +rprimd(3,idir)*nhatfr_tmp(3,ic)
       end do
       if (ider==1) then
         do nu=1,3
           do ic=1,nfgd
             pawfgrtab(iatom)%nhatfrgr(nu,ic,ispden)= &
&             rprimd(1,idir)*nhatfrgr_tmp(1,ic,nu) &
&             +rprimd(2,idir)*nhatfrgr_tmp(2,ic,nu) &
&             +rprimd(3,idir)*nhatfrgr_tmp(3,ic,nu)
           end do
         end do
       end if
       ABI_DEALLOCATE(nhatfr_tmp)
       if (ider==1) then
         ABI_DEALLOCATE(nhatfrgr_tmp)
       end if
!      End loop over spin components
     end do ! ispden


!  ============ Elastic tensor ===============================
   else if (ipert==natom+3.or.ipert==natom+4) then
!    Loop over spin components
     pawfgrtab(iatom)%nhatfr(:,:) = zero
     do ispden=1,nspden
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         klm =pawtab(itypat)%indklmn(1,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         if (nspden/=2) then
           ro=pawrhoij(iatom)%rhoijp(jrhoij,ispden)
         else
           if (ispden==1) then
             ro=pawrhoij(iatom)%rhoijp(jrhoij,1)+pawrhoij(iatom)%rhoijp(jrhoij,2)
           else if (ispden==2) then
             ro=pawrhoij(iatom)%rhoijp(jrhoij,1)
           end if
         end if
         ro=pawtab(itypat)%dltij(klmn)*ro
         do ils=lmin,lmax,2
           lm0=ils**2+ils+1
           do mm=-ils,ils
             ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
             if (isel>0) then
!              Sum{[Q_ij_q^LM^(1)]}
               do ic=1,nfgd
                 mua=alpha(istr);mub=beta(istr)
                 pawfgrtab(iatom)%nhatfr(ic,ispden) = pawfgrtab(iatom)%nhatfr(ic,ispden)+&
&                 ro*pawtab(itypat)%qijl(ilslm,klmn)*half*(&
&                 pawfgrtab(iatom)%gylmgr(mua,ic,ilslm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&                 +pawfgrtab(iatom)%gylmgr(mub,ic,ilslm)*pawfgrtab(iatom)%rfgd(mua,ic))
               end do
!              Add volume contribution
               if(istr<=3)then
                 do ic=1,nfgd
                   pawfgrtab(iatom)%nhatfr(ic,ispden) = pawfgrtab(iatom)%nhatfr(ic,ispden)+&
&                   ro*pawtab(itypat)%qijl(ilslm,klmn)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
               if (ider==1) then
                 MSG_ERROR("nhatgr not implemented for strain perturbationxs")
!                 do ic=1,nfgd
!                   do nu=1,6
!                     do mu=1,6
!                       contrib=-ro*pawtab(itypat)%qijl(ilslm,klmn) &
!&                       *pawfgrtab(iatom)%gylmgr2(voigt(mu,nu),ic,ilslm)
!                       nhatfrgr_tmp(mu,ic,nu)=nhatfrgr_tmp(mu,ic,nu)+contrib
!                     end do
!                   end do
!                 end do
               end if
             end if
           end do
         end do
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do
     end do ! ispden
   end if

!  Eventually free temporary space for g_l(r).Y_lm(r) gradients and exp(-i.q.r)
   if (pawfgrtab(iatom)%gylmgr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
     pawfgrtab(iatom)%gylmgr_allocated=0
   end if
   if (pawfgrtab(iatom)%gylmgr2_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(0,0,0))
     pawfgrtab(iatom)%gylmgr2_allocated=0
   end if

!  End loop on atoms
 end do


!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawnhatfr
!!***

!!****m* ABINIT/m_paw_nhat
!! NAME
!!  m_paw_nhat
!!
!! FUNCTION
!!  This module contains several routines related to the PAW compensation
!!    charge density (i.e. n^hat(r)).
!!
!! COPYRIGHT
!! Copyright (C) 2018-2020 ABINIT group (FJ, MT, MG, TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_nhat

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes,  only : MPI_type
 use m_time,         only : timab
 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_pawrhoij,     only : pawrhoij_type
 use m_pawcprj,      only : pawcprj_type
 use m_paw_finegrid, only : pawgylm,pawrfgd_fft,pawrfgd_wvl,pawexpiqr
 use m_paral_atom,   only : get_my_atmtab, free_my_atmtab
 use m_distribfft,   only : distribfft_type,init_distribfft_seq,destroy_distribfft
 use m_geometry,     only : xred2xcart
 use m_cgtools,      only : mean_fftr
 use m_mpinfo,       only : set_mpi_enreg_fft,unset_mpi_enreg_fft,initmpi_seq
 use m_fft,          only : zerosym, fourwf, fourdp
 use m_paw_lmn,      only : klmn2ijlmn

 implicit none

 private

!public procedures.
 public :: pawmknhat        ! Compute compensation charge density on the real space (fine) grid
 public :: pawmknhat_psipsi ! Compute compensation charge density associated to the product of two WF
 public :: pawnhatfr        ! Compute frozen part of 1st-order compensation charge density nhat^(1) (DFPT)
 public :: pawsushat        ! Compute contrib. to the product of two WF from compensation charge density
 public :: nhatgrid         ! Determine points of the (fine) grid that are located around atoms - PW version
 public :: wvl_nhatgrid     ! Determine points of the (fine) grid that are located around atoms - WVL version

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nhat/pawmknhat
!! NAME
!! pawmknhat
!!
!! FUNCTION
!! PAW only:
!! Compute compensation charge density (and derivatives) on the fine FFT grid
!! Can also compute first-order compensation charge density (RF calculations)
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  ider= 0: nhat(r) is computed
!!        1: cartesian derivatives of nhat(r) are computed
!!        2: nhat(r) and derivatives are computed
!!  idir=direction of atomic displacement (in case of atomic displ. perturb.)
!!  ipert=index of perturbation; must be 0 for ground-state calculations
!!  izero=if 1, unbalanced components of nhat(g) have to be set to zero
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhatgrdim= -PAW only- 0 if pawgrnhat array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         (1st-order occupancies if ipert>0)
!!  pawrhoij0(my_natom) <type(pawrhoij_type)>= GS paw rhoij occupancies and related data (used only if ipert>0)
!!                                          set equat to pawrhoij for GS calculations
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon (RF only)
!!  rprimd(3,3)=dimensional primitive translations for real space
!!  ucvol=volume of the unit cell
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    compch_fft=compensation charge inside spheres computed over fine fft grid
!!    pawnhat(nfft,ispden)=nhat on fine rectangular grid
!!  === if ider=1 or 2
!!    pawgrnhat(nfft,ispden,3)=derivatives of nhat on fine rectangular grid (and derivatives)
!!
!! PARENTS
!!      bethe_salpeter,dfpt_scfcv,energy,nres2vres,odamix,paw_qpscgw,pawmkrho
!!      respfn,scfcv,screening,setup_positron,sigma
!!
!! CHILDREN
!!      destroy_distribfft,fourdp,free_my_atmtab,get_my_atmtab
!!      init_distribfft_seq,initmpi_seq,mean_fftr,pawexpiqr,pawgylm,pawnhatfr
!!      set_mpi_enreg_fft,timab,unset_mpi_enreg_fft,xmpi_sum,zerosym
!!
!! SOURCE

subroutine pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,&
&          my_natom,natom,nfft,ngfft,nhatgrdim,nspden,ntypat,pawang,pawfgrtab,&
&          pawgrnhat,pawnhat,pawrhoij,pawrhoij0,pawtab,qphon,rprimd,ucvol,usewvl,xred,&
&          mpi_atmtab,comm_atom,comm_fft,mpi_comm_wvl,me_g0,paral_kgb,distribfft) ! optional arguments

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ider,idir,ipert,izero,my_natom,natom,nfft
 integer,intent(in)  :: usewvl
 integer,intent(in) :: nhatgrdim,nspden,ntypat
 integer,optional,intent(in) :: me_g0,comm_atom,comm_fft,mpi_comm_wvl,paral_kgb
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: compch_fft
 type(distribfft_type),optional,intent(in),target :: distribfft
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: pawgrnhat(cplex*nfft,nspden,3*nhatgrdim)
 real(dp),intent(inout) :: pawnhat(cplex*nfft,nspden) !vz_i
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom),pawrhoij0(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iatom,iatom_tot,ic,ierr,ii,ils,ilslm,iq0,irhoij,ispden,itypat
 integer :: jc,jrhoij,kc,klm,klmn,lmax,lmin,lm_size,mfgd,mm,mpi_comm_sphgrid
 integer :: my_comm_atom,my_comm_fft,nfgd,nfftot,option,optgr0,optgr1,optgr2,paral_kgb_fft
 logical :: compute_grad,compute_nhat,my_atmtab_allocated,need_frozen,paral_atom,qeq0
 logical :: compute_phonons,has_phase
 type(distribfft_type),pointer :: my_distribfft
 type(mpi_type) :: mpi_enreg_fft
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: ro(cplex),ro_ql(cplex),tmp_compch_fft(nspden),tsec(2)
 real(dp),allocatable :: pawgrnhat_atm(:,:),pawnhat_atm(:),work(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 compute_nhat=(ider==0.or.ider==2)
 compute_grad=(ider==1.or.ider==2)
 compute_phonons=(ipert>0.and.ipert<=natom)

!Compatibility tests
 qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
 if (present(comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
     MSG_BUG('Need paral_kgb and me_g0 with comm_fft !')
   end if
 end if
 if(ider>0.and.nhatgrdim==0) then
   MSG_BUG('Gradients of nhat required but not allocated!')
 end if
 if (my_natom>0) then
   if(nspden>1.and.nspden/=pawrhoij(1)%nspden) then
     MSG_BUG('Wrong values for nspden and pawrhoij%nspden!')
   end if
   if(nspden>1.and.nspden/=pawfgrtab(1)%nspden) then
     MSG_BUG('Wrong values for nspden and pawfgrtab%nspden!')
   end if
   if(pawrhoij(1)%qphase<cplex) then
     MSG_BUG('Must have pawrhoij()%qphase >= cplex!')
   end if
   if (compute_phonons.and.(.not.qeq0)) then
     if (pawfgrtab(1)%rfgd_allocated==0) then
       MSG_BUG('pawfgrtab()%rfgd array must be allocated!')
     end if
     if (compute_grad.and.(.not.compute_nhat)) then
       MSG_BUG('When q<>0, nhat gradients need nhat!')
     end if
   end if
 end if

!nhat1 does not have to be computed for ddk or d2dk
 if (ipert==natom+1.or.ipert==natom+10) return

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

!Initialisations
 if ((.not.compute_nhat).and.(.not.compute_grad)) return
 mfgd=zero;if (my_natom>0) mfgd=maxval(pawfgrtab(1:my_natom)%nfgd)
 if (compute_nhat) then
   ABI_ALLOCATE(pawnhat_atm,(cplex*mfgd))
   pawnhat=zero
 end if
 if (compute_grad) then
   ABI_ALLOCATE(pawgrnhat_atm,(cplex*mfgd,3))
   pawgrnhat=zero
 end if

!mpi communicators for spherical grid:
 mpi_comm_sphgrid=xmpi_comm_self !no communicators passed
 if(present(comm_fft) .and. usewvl==0) mpi_comm_sphgrid=comm_fft
 if(present(mpi_comm_wvl) .and. usewvl==1) mpi_comm_sphgrid=mpi_comm_wvl

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat=pawrhoij(iatom)%itypat
   lm_size=pawfgrtab(iatom)%l_size**2
   need_frozen=((compute_nhat).and.(ipert==iatom_tot.or.ipert==natom+3.or.ipert==natom+4))
   nfgd=pawfgrtab(iatom)%nfgd
   cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
   iq0=cplex_rhoij*pawrhoij(iatom)%lmn2_size

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&   ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)).or.&
&   ((compute_grad.and.need_frozen).and.(pawfgrtab(iatom)%gylmgr2_allocated==0))) then
     optgr0=0;optgr1=0;optgr2=0
     if ((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
     end if
     if ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
     end if
     if ((compute_grad.and.need_frozen).and.(pawfgrtab(iatom)%gylmgr2_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylmgr2))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(6,nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr2_allocated=2;optgr2=1
     end if
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),pawfgrtab(iatom)%rfgd)
     end if
   end if


!  Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
   if (compute_phonons.and.(.not.qeq0).and.pawfgrtab(iatom)%expiqr_allocated==0) then
     if (allocated(pawfgrtab(iatom)%expiqr))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
     call pawexpiqr(pawfgrtab(iatom)%expiqr,gprimd,nfgd,qphon,&
&     pawfgrtab(iatom)%rfgd,xred(:,iatom_tot))
     pawfgrtab(iatom)%expiqr_allocated=2
   end if
   has_phase=(compute_phonons.and.pawfgrtab(iatom)%expiqr_allocated/=0)

!  Eventually compute frozen part of nhat for the current atom (if not already done)
   if ((need_frozen).and.((pawfgrtab(iatom)%nhatfr_allocated==0).or.&
&   (compute_grad.and.pawfgrtab(iatom)%nhatfrgr_allocated==0))) then
     if (allocated(pawfgrtab(iatom)%nhatfr))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%nhatfr,(nfgd,pawfgrtab(iatom)%nspden))
     option=0;pawfgrtab(iatom)%nhatfr_allocated=2
     if (compute_grad) then
       option=1
       if (allocated(pawfgrtab(iatom)%nhatfrgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfrgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%nhatfrgr,(3,nfgd,pawfgrtab(iatom)%nspden))
       pawfgrtab(iatom)%nhatfrgr_allocated=2
     end if
     call pawnhatfr(option,idir,ipert,1,natom,nspden,ntypat,pawang,pawfgrtab(iatom),&
&                   pawrhoij0(iatom),pawtab,rprimd)
   end if

!  ------------------------------------------------------------------------
!  ----- Loop over density components
!  ------------------------------------------------------------------------

   do ispden=1,nspden

     if (compute_nhat) pawnhat_atm(1:cplex*nfgd)=zero
     if (compute_grad) pawgrnhat_atm(1:cplex*nfgd,1:3)=zero

!    ------------------------------------------------------------------------
!    ----- Loop over ij channels (basis components)
!    ------------------------------------------------------------------------
     jrhoij=1
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)
       lmax=pawtab(itypat)%indklmn(4,klmn)

!      Retrieve rhoij
       if (pawrhoij(iatom)%nspden/=2) then
         ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)
         if (cplex==2) ro(2)=pawrhoij(iatom)%rhoijp(iq0+jrhoij,ispden)
       else
         if (ispden==1) then
           ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,1)+pawrhoij(iatom)%rhoijp(jrhoij,2)
           if (cplex==2) ro(2)=pawrhoij(iatom)%rhoijp(iq0+jrhoij,1)+pawrhoij(iatom)%rhoijp(iq0+jrhoij,2)
         else if (ispden==2) then
           ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,1)
           if (cplex==2) ro(2)=pawrhoij(iatom)%rhoijp(iq0+jrhoij,1)
         end if
       end if
       ro(1:cplex)=pawtab(itypat)%dltij(klmn)*ro(1:cplex)

       if (compute_nhat) then
         if (cplex==1) then
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   pawnhat_atm(ic)=pawnhat_atm(ic)+ro_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
             end do
           end do
         else
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1:2)=ro(1:2)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   jc=2*ic-1
                   pawnhat_atm(jc:jc+1)=pawnhat_atm(jc:jc+1)+ro_ql(1:2)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
             end do
           end do
         end if
       end if

       if (compute_grad) then
         if (cplex==1) then
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   pawgrnhat_atm(ic,1)=pawgrnhat_atm(ic,1)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                   pawgrnhat_atm(ic,2)=pawgrnhat_atm(ic,2)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                   pawgrnhat_atm(ic,3)=pawgrnhat_atm(ic,3)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
                 end do
               end if
             end do
           end do
         else
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1:2)=ro(1:2)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   jc=2*ic-1
                   pawgrnhat_atm(jc:jc+1,1)=pawgrnhat_atm(jc:jc+1,1) &
&                   +ro_ql(1:2)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                   pawgrnhat_atm(jc:jc+1,2)=pawgrnhat_atm(jc:jc+1,2) &
&                   +ro_ql(1:2)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                   pawgrnhat_atm(jc:jc+1,3)=pawgrnhat_atm(jc:jc+1,3) &
&                   +ro_ql(1:2)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
                 end do
               end if
             end do
           end do
         end if
       end if

!      ------------------------------------------------------------------------
!      ----- End loop over ij channels
!      ------------------------------------------------------------------------
       jrhoij=jrhoij+cplex_rhoij
     end do

!    If RF calculation, add frozen part of 1st-order compensation density
     if (need_frozen) then
       if (cplex==1) then
         do ic=1,nfgd
           pawnhat_atm(ic)=pawnhat_atm(ic)+pawfgrtab(iatom)%nhatfr(ic,ispden)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1
           pawnhat_atm(jc)=pawnhat_atm(jc)+pawfgrtab(iatom)%nhatfr(ic,ispden)
         end do
       end if
       if (compute_grad) then
         if (cplex==1) then
           do ic=1,nfgd
             pawgrnhat_atm(ic,1)=pawgrnhat_atm(ic,1)+pawfgrtab(iatom)%nhatfrgr(1,ic,ispden)
             pawgrnhat_atm(ic,2)=pawgrnhat_atm(ic,2)+pawfgrtab(iatom)%nhatfrgr(2,ic,ispden)
             pawgrnhat_atm(ic,3)=pawgrnhat_atm(ic,3)+pawfgrtab(iatom)%nhatfrgr(3,ic,ispden)
           end do
         else
           do ic=1,nfgd
             jc=2*ic-1
             pawgrnhat_atm(jc,1)=pawgrnhat_atm(jc,1)+pawfgrtab(iatom)%nhatfrgr(1,ic,ispden)
             pawgrnhat_atm(jc,2)=pawgrnhat_atm(jc,2)+pawfgrtab(iatom)%nhatfrgr(2,ic,ispden)
             pawgrnhat_atm(jc,3)=pawgrnhat_atm(jc,3)+pawfgrtab(iatom)%nhatfrgr(3,ic,ispden)
           end do
         end if
       end if
     end if

!    If needed, multiply eventually by exp(-i.q.r) phase
     if (has_phase) then
       if (cplex==1) then
         do ic=1,nfgd
           pawnhat_atm(ic)=pawnhat_atm(ic)*pawfgrtab(iatom)%expiqr(1,ic)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
           ro(1:2)=pawnhat_atm(jc:jc+1)
           pawnhat_atm(jc  )=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           pawnhat_atm(jc+1)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
       if (compute_grad) then
         if (cplex==1) then
           do ic=1,nfgd
             pawgrnhat_atm(ic,1:3)=pawgrnhat_atm(ic,1:3)*pawfgrtab(iatom)%expiqr(1,ic)
           end do
         else
           do ic=1,nfgd
             jc=2*ic-1
!            dn^hat(r)/dr_i * exp(-i.q.r)
             ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
             ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
             do ii=1,3
               ro(1:2)=pawgrnhat_atm(jc:jc+1,ii)
               pawgrnhat_atm(jc  ,ii)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
               pawgrnhat_atm(jc+1,ii)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
             end do
!            -i.q_i * [n^hat(r).exp(-i.q.r)]
             ro(1:2)=pawnhat_atm(jc:jc+1)
             do ii=1,3
               pawgrnhat_atm(jc  ,ii)=pawgrnhat_atm(kc  ,ii)+qphon(ii)*ro(2)
               pawgrnhat_atm(jc+1,ii)=pawgrnhat_atm(kc+1,ii)-qphon(ii)*ro(1)
             end do
           end do
         end if
       end if
     end if

!    Add the contribution of the atom to the compensation charge
     if (compute_nhat) then
       if (cplex==1) then
         do ic=1,nfgd
           kc=pawfgrtab(iatom)%ifftsph(ic)
           pawnhat(kc,ispden)=pawnhat(kc,ispden)+pawnhat_atm(ic)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1;kc=2*pawfgrtab(iatom)%ifftsph(ic)-1
           pawnhat(kc:kc+1,ispden)=pawnhat(kc:kc+1,ispden)+pawnhat_atm(jc:jc+1)
         end do
       end if
     end if
     if (compute_grad) then
       if (cplex==1) then
         do ic=1,nfgd
           kc=pawfgrtab(iatom)%ifftsph(ic)
           pawgrnhat(kc,ispden,1:3)=pawgrnhat(kc,ispden,1:3)+pawgrnhat_atm(ic,1:3)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1;kc=2*pawfgrtab(iatom)%ifftsph(ic)-1
           do ii=1,3
             pawgrnhat(kc:kc+1,ispden,ii)=pawgrnhat(kc:kc+1,ispden,ii)+pawgrnhat_atm(jc:jc+1,ii)
           end do
         end do
       end if
     end if
!    ------------------------------------------------------------------------
!    ----- End loop over density components
!    ------------------------------------------------------------------------
   end do

   if (pawfgrtab(iatom)%gylm_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(0,0))
     pawfgrtab(iatom)%gylm_allocated=0
   end if
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
   if (pawfgrtab(iatom)%nhatfr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     ABI_ALLOCATE(pawfgrtab(iatom)%nhatfr,(0,0))
     pawfgrtab(iatom)%nhatfr_allocated=0
   end if
   if (pawfgrtab(iatom)%nhatfrgr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfrgr)
     ABI_ALLOCATE(pawfgrtab(iatom)%nhatfrgr,(0,0,0))
     pawfgrtab(iatom)%nhatfrgr_allocated=0
   end if
   if (pawfgrtab(iatom)%expiqr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(0,0))
     pawfgrtab(iatom)%expiqr_allocated=0
   end if

!  ------------------------------------------------------------------------
!  ----- End loop over atoms
!  ------------------------------------------------------------------------
 end do

!----- Free some memory
 if (compute_nhat) then
   ABI_DEALLOCATE(pawnhat_atm)
 end if
 if (compute_grad) then
   ABI_DEALLOCATE(pawgrnhat_atm)
 end if

!----- Reduction in case of parallelism
 if (paral_atom) then
   call timab(48,1,tsec)
   if (compute_nhat) then
     call xmpi_sum(pawnhat,my_comm_atom,ierr)
   end if
   if (compute_grad) then
     call xmpi_sum(pawgrnhat,my_comm_atom,ierr)
   end if
   call timab(48,2,tsec)
 end if

!----- Avoid unbalanced g-components numerical errors
 if (izero==1.and.compute_nhat.and.usewvl==0) then
!  Create fake mpi_enreg to wrap fourdp
   if (present(distribfft)) then
     my_distribfft => distribfft
   else
     ABI_DATATYPE_ALLOCATE(my_distribfft,)
     call init_distribfft_seq(my_distribfft,'f',ngfft(2),ngfft(3),'fourdp')
   end if
   call initmpi_seq(mpi_enreg_fft)
   ABI_DATATYPE_DEALLOCATE(mpi_enreg_fft%distribfft)
   if (present(comm_fft)) then
     call set_mpi_enreg_fft(mpi_enreg_fft,comm_fft,my_distribfft,me_g0,paral_kgb)
     my_comm_fft=comm_fft;paral_kgb_fft=paral_kgb
   else
     my_comm_fft=xmpi_comm_self;paral_kgb_fft=0;
     mpi_enreg_fft%distribfft => my_distribfft
   end if
!  do FFT
   ABI_ALLOCATE(work,(2,nfft))
   do ispden=1,min(2,nspden)
     call fourdp(cplex,work,pawnhat(:,ispden),-1,mpi_enreg_fft,nfft,1,ngfft,0)
     call zerosym(work,2,ngfft(1),ngfft(2),ngfft(3),comm_fft=my_comm_fft,distribfft=my_distribfft)
     call fourdp(cplex,work,pawnhat(:,ispden),+1,mpi_enreg_fft,nfft,1,ngfft,0)
   end do
   ABI_DEALLOCATE(work)
!  Destroy fake mpi_enreg
   call unset_mpi_enreg_fft(mpi_enreg_fft)
   if (.not.present(distribfft)) then
     call destroy_distribfft(my_distribfft)
     ABI_DATATYPE_DEALLOCATE(my_distribfft)
   end if
 end if

!----- Computation of compensation charge over real space grid
 if (compute_nhat.and.ipert==0) then
   nfftot=PRODUCT(ngfft(1:3))
   call mean_fftr(pawnhat,tmp_compch_fft,nfft,nfftot,1,mpi_comm_sphgrid)
   compch_fft = tmp_compch_fft(1)
   compch_fft=compch_fft*ucvol
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawmknhat
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nhat/pawmknhat_psipsi
!! NAME
!! pawmknhat_psipsi
!!
!! FUNCTION
!! PAW only:
!! Compute on the fine FFT grid the compensation charge density (and derivatives) associated
!! to the product of two wavefunctions n_{12}(r) = \Psi_1* \Psi_2. Based on pawmknhat.
!!
!! INPUTS
!!  cprj1(natom,nspinor), cprj2(natom,nspinor) <type(pawcprj_type)>=
!!   projected input wave functions <Proj_i|Cnk> with all NL projectors corresponding to
!!   the \Psi_1 and \Psi_2, respectively.
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  ider= 0: nhat(r) is computed
!!        1: cartesian derivatives of nhat(r) are computed
!!        2: nhat(r) and derivatives are computed
!!        3: nhat(r) and gradients of nhat  wrt atomic coordinates are computed
!!        Note: ider>0 not compatible with ipert>0
!!  izero=if 1, unbalanced components of nhat(g) have to be set to zero
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat12_grdim= 0 if grnhat12 array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    nhat12(2,nfft,nspinor**2)=nhat on fine rectangular grid*exp(iqr)
!!  === if ider=1 or 2
!!    grnhat12(nfft,nspinor**2,3)=gradient of (nhat*exp(iqr)) on fine rectangular grid (derivative versus r)
!!  === if ider=3
!!    grnhat_12(nfft,nspinor**2,3,natom*(ider/3))=derivatives of nhat on fine rectangular grid versus R*exp(iqr)
!!
!! PARENTS
!!      calc_sigx_me,fock_getghc,prep_calc_ucrpa
!!
!! CHILDREN
!!      destroy_distribfft,fourdp,free_my_atmtab,get_my_atmtab
!!      init_distribfft_seq,initmpi_seq,pawexpiqr,pawgylm,set_mpi_enreg_fft
!!      timab,unset_mpi_enreg_fft,xmpi_sum,zerosym
!!
!! SOURCE

subroutine pawmknhat_psipsi(cprj1,cprj2,ider,izero,my_natom,natom,nfft,ngfft,nhat12_grdim,&
&          nspinor,ntypat,pawang,pawfgrtab,grnhat12,nhat12,pawtab, &
&          gprimd,grnhat_12,qphon,xred,atindx,mpi_atmtab,comm_atom,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ider,izero,my_natom,natom,nfft,nhat12_grdim,ntypat,nspinor
 integer,optional,intent(in) :: me_g0,comm_fft,paral_kgb
 integer,optional,intent(in) :: comm_atom
 type(distribfft_type),optional,intent(in),target :: distribfft
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 integer,optional,intent(in) ::atindx(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),optional, intent(in) ::gprimd(3,3),qphon(3),xred(3,natom)
 real(dp),intent(out) :: grnhat12(2,nfft,nspinor**2,3*nhat12_grdim)
 real(dp),optional,intent(out) :: grnhat_12(2,nfft,nspinor**2,3,natom*(ider/3))
 real(dp),intent(out) :: nhat12(2,nfft,nspinor**2)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(pawcprj_type),intent(in) :: cprj1(natom,nspinor),cprj2(natom,nspinor)

!Local variables ---------------------------------------
!scalars
 integer :: iatm,iatom,iatom_tot,ic,ierr,ils,ilslm,isp1,isp2,isploop,itypat,jc,klm,klmn
 integer :: lmax,lmin,lm_size,mm,my_comm_atom,my_comm_fft,optgr0,optgr1,paral_kgb_fft
 integer :: cplex,ilmn,jlmn,lmn_size,lmn2_size
 real(dp) :: re_p,im_p
 logical :: compute_grad,compute_grad1,compute_nhat,my_atmtab_allocated,paral_atom,qeq0,compute_phonon,order
 type(distribfft_type),pointer :: my_distribfft
 type(mpi_type) :: mpi_enreg_fft
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 integer,pointer :: my_atmtab(:)
 real(dp) :: rdum(1),cpf(2),cpf_ql(2),tsec(2),ro(2),ro_ql(2),nhat12_atm(2,nfft,nspinor**2)
 real(dp),allocatable :: work(:,:), qijl(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if (present(comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
     MSG_BUG('Need paral_kgb and me_g0 with comm_fft!')
   end if
   if (present(paral_kgb)) then
     if (paral_kgb/=0) then
       MSG_BUG('paral_kgb/=0 not coded!')
     end if
   end if
 end if
 if (ider>0.and.nhat12_grdim==0) then
!   MSG_BUG('Gradients of nhat required but not allocated !')
 end if
 if (nspinor==2) then
   MSG_BUG('nspinor==2 not coded!')
 end if

 compute_phonon=.false.;qeq0=.false.
 if (present(gprimd).and.present(qphon).and.present(xred)) compute_phonon=.true.
 if (compute_phonon) qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
 if (present(atindx)) order=.true.
!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Initialisations
 compute_nhat=(ider==0.or.ider==2.or.ider==3)
 compute_grad=(ider==1.or.ider==2)
 compute_grad1=(ider==3)
 if ((.not.compute_nhat).and.(.not.compute_grad)) return

 if (compute_nhat) nhat12=zero
 if (compute_grad) grnhat12=zero
 if (compute_grad1) grnhat_12=zero

 if (compute_grad) then
!   MSG_BUG('compute_grad not tested!')
 end if

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   iatm=iatom_tot
   if (order) iatm=atindx(iatom_tot)
   itypat    = pawfgrtab(iatom)%itypat
   lm_size   = pawfgrtab(iatom)%l_size**2
   lmn_size  = pawtab(itypat)%lmn_size
   lmn2_size = pawtab(itypat)%lmn2_size
   ABI_ALLOCATE(qijl,(lm_size,lmn2_size))
   qijl=zero
   qijl=pawtab(itypat)%qijl
   if (compute_nhat) nhat12_atm=zero

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&   (((compute_grad).or.(compute_grad1)).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then
     optgr0=0; optgr1=0
     if ((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
     end if
     if (((compute_grad).or.(compute_grad1)).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
     end if
     if (optgr0+optgr1>0) then
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,rdum,&
&       lm_size,pawfgrtab(iatom)%nfgd,optgr0,optgr1,0,pawtab(itypat),&
&       pawfgrtab(iatom)%rfgd)
     end if

   end if
   if (compute_phonon.and.(.not.qeq0).and.(pawfgrtab(iatom)%expiqr_allocated==0)) then
     if (allocated(pawfgrtab(iatom)%expiqr))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,pawfgrtab(iatom)%nfgd))
     call pawexpiqr(pawfgrtab(iatom)%expiqr,gprimd,pawfgrtab(iatom)%nfgd,qphon,&
&     pawfgrtab(iatom)%rfgd,xred(:,iatom_tot))
     pawfgrtab(iatom)%expiqr_allocated=2
   end if

   do isploop=1,nspinor**2    ! Loop over density components of the compensation charge.
!    TODO Here we might take advantage of symmetry relations between the four components if nspinor==2
     isp1=spinor_idxs(1,isploop)
     isp2=spinor_idxs(2,isploop)

     do klmn=1,lmn2_size  ! Loop over ij channels of this atom type.
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)  ! abs(il-jl)
       lmax=pawtab(itypat)%indklmn(4,klmn)  ! il+jl
       ilmn=pawtab(itypat)%indklmn(7,klmn)
       jlmn=pawtab(itypat)%indklmn(8,klmn)
!      call klmn2ijlmn(klmn,lmn_size,ilmn,jlmn)  ! This mapping should be stored in pawtab_type

!      Retrieve the factor due to the PAW projections.
       re_p =  cprj1(iatm,isp1)%cp(1,ilmn) * cprj2(iatm,isp2)%cp(1,jlmn) &
&       +cprj1(iatm,isp1)%cp(2,ilmn) * cprj2(iatm,isp2)%cp(2,jlmn) &
&       +cprj1(iatm,isp1)%cp(1,jlmn) * cprj2(iatm,isp2)%cp(1,ilmn) &
&       +cprj1(iatm,isp1)%cp(2,jlmn) * cprj2(iatm,isp2)%cp(2,ilmn)

       im_p =  cprj1(iatm,isp1)%cp(1,ilmn) * cprj2(iatm,isp2)%cp(2,jlmn) &
&       -cprj1(iatm,isp1)%cp(2,ilmn) * cprj2(iatm,isp2)%cp(1,jlmn) &
&       +cprj1(iatm,isp1)%cp(1,jlmn) * cprj2(iatm,isp2)%cp(2,ilmn) &
&       -cprj1(iatm,isp1)%cp(2,jlmn) * cprj2(iatm,isp2)%cp(1,ilmn)

       cpf(1)=re_p*pawtab(itypat)%dltij(klmn)*half
       cpf(2)=im_p*pawtab(itypat)%dltij(klmn)*half

       if (compute_nhat) then
         do ils=lmin,lmax,2   ! Sum over (L,M)
           do mm=-ils,ils
             ilslm=ils*ils+ils+mm+1
             if (pawang%gntselect(ilslm,klm)>0) then
               cpf_ql(1)=cpf(1)*qijl(ilslm,klmn)
               cpf_ql(2)=cpf(2)*qijl(ilslm,klmn)
               do ic=1,pawfgrtab(iatom)%nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 nhat12_atm(1,jc,isploop)=nhat12_atm(1,jc,isploop)+cpf_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 nhat12_atm(2,jc,isploop)=nhat12_atm(2,jc,isploop)+cpf_ql(2)*pawfgrtab(iatom)%gylm(ic,ilslm)
               end do
             end if
           end do
         end do
       end if ! compute_nhat

       if (compute_grad) then
         do ils=lmin,lmax,2  ! Sum over (L,M)
           do mm=-ils,ils
             ilslm=ils*ils+ils+mm+1
             if (pawang%gntselect(ilslm,klm)>0) then
               cpf_ql(1)=cpf(1)*qijl(ilslm,klmn)
               cpf_ql(2)=cpf(2)*qijl(ilslm,klmn)
               do ic=1,pawfgrtab(iatom)%nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 grnhat12(1,jc,isploop,1)=grnhat12(1,jc,isploop,1)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 grnhat12(1,jc,isploop,2)=grnhat12(1,jc,isploop,2)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 grnhat12(1,jc,isploop,3)=grnhat12(1,jc,isploop,3)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)

                 grnhat12(2,jc,isploop,1)=grnhat12(2,jc,isploop,1)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 grnhat12(2,jc,isploop,2)=grnhat12(2,jc,isploop,2)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 grnhat12(2,jc,isploop,3)=grnhat12(2,jc,isploop,3)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
               end do
             end if
           end do
         end do
       end if ! compute_grad
       if (compute_grad1) then
         do ils=lmin,lmax,2  ! Sum over (L,M)
           do mm=-ils,ils
             ilslm=ils*ils+ils+mm+1
             if (pawang%gntselect(ilslm,klm)>0) then
               cpf_ql(1)=cpf(1)*qijl(ilslm,klmn)
               cpf_ql(2)=cpf(2)*qijl(ilslm,klmn)
               do ic=1,pawfgrtab(iatom)%nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 grnhat_12(1,jc,isploop,1,iatom)=grnhat_12(1,jc,isploop,1,iatom)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 grnhat_12(1,jc,isploop,2,iatom)=grnhat_12(1,jc,isploop,2,iatom)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 grnhat_12(1,jc,isploop,3,iatom)=grnhat_12(1,jc,isploop,3,iatom)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)

                 grnhat_12(2,jc,isploop,1,iatom)=grnhat_12(2,jc,isploop,1,iatom)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 grnhat_12(2,jc,isploop,2,iatom)=grnhat_12(2,jc,isploop,2,iatom)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 grnhat_12(2,jc,isploop,3,iatom)=grnhat_12(2,jc,isploop,3,iatom)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
               end do
             end if
           end do
         end do
       end if ! compute_grad1
     end do  ! klmn (ij channels)
!    If needed, multiply eventually by exp(-i.q.r) phase
     if (compute_nhat) then
       if(compute_phonon.and.(.not.qeq0).and.pawfgrtab(iatom)%expiqr_allocated/=0) then
         do ic=1,pawfgrtab(iatom)%nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)= pawfgrtab(iatom)%expiqr(2,ic)
           ro(1:2)=nhat12_atm(1:2,jc,isploop)
           nhat12_atm(1,jc,isploop)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           nhat12_atm(2,jc,isploop)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
     end if

     if (compute_grad) then
       if(compute_phonon.and.(.not.qeq0).and.pawfgrtab(iatom)%expiqr_allocated/=0) then
         do ic=1,pawfgrtab(iatom)%nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)= pawfgrtab(iatom)%expiqr(2,ic)
           ro(1)=grnhat12(1,jc,isploop,1)-qphon(1)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat12(2,jc,isploop,1)+qphon(1)*nhat12_atm(1,jc,isploop)
           grnhat12(1,jc,isploop,1)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat12(2,jc,isploop,1)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat12(1,jc,isploop,2)-qphon(2)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat12(2,jc,isploop,2)+qphon(2)*nhat12_atm(1,jc,isploop)
           grnhat12(1,jc,isploop,2)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat12(2,jc,isploop,2)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat12(1,jc,isploop,3)-qphon(3)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat12(2,jc,isploop,3)+qphon(3)*nhat12_atm(1,jc,isploop)
           grnhat12(1,jc,isploop,3)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat12(2,jc,isploop,3)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
     end if
     if (compute_grad1) then
       if(compute_phonon.and.(.not.qeq0).and.pawfgrtab(iatom)%expiqr_allocated/=0) then
         do ic=1,pawfgrtab(iatom)%nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)= pawfgrtab(iatom)%expiqr(2,ic)
           ro(1)=grnhat_12(1,jc,isploop,1,iatom)
           ro(2)=grnhat_12(2,jc,isploop,1,iatom)
           grnhat_12(1,jc,isploop,1,iatom)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat_12(2,jc,isploop,1,iatom)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat_12(1,jc,isploop,2,iatom)
           ro(2)=grnhat_12(2,jc,isploop,2,iatom)
           grnhat_12(1,jc,isploop,2,iatom)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat_12(2,jc,isploop,2,iatom)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat_12(1,jc,isploop,3,iatom)
           ro(2)=grnhat_12(2,jc,isploop,3,iatom)
           grnhat_12(1,jc,isploop,3,iatom)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat_12(2,jc,isploop,3,iatom)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
     end if
   end do ! isploop (density components of the compensation charge)
! accumlate nhat12 for all the atoms
   if (compute_nhat) nhat12=nhat12+nhat12_atm

   if (pawfgrtab(iatom)%gylm_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(0,0))
     pawfgrtab(iatom)%gylm_allocated=0
   end if
   if (pawfgrtab(iatom)%gylmgr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
     pawfgrtab(iatom)%gylmgr_allocated=0
   end if
   ABI_DEALLOCATE(qijl)
   if (pawfgrtab(iatom)%expiqr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(0,0))
     pawfgrtab(iatom)%expiqr_allocated=0
   end if

 end do ! iatom

 if (compute_grad1) grnhat_12=-grnhat_12

!----- Reduction in case of parallelism -----!
 if (paral_atom)then
   call timab(48,1,tsec)
   if (compute_nhat) then
     call xmpi_sum(nhat12,my_comm_atom,ierr)
   end if
   if (compute_grad) then
     call xmpi_sum(grnhat12,my_comm_atom,ierr)
   end if
   if (compute_grad1) then
     call xmpi_sum(grnhat_12,my_comm_atom,ierr)
   end if
   call timab(48,2,tsec)
 end if

!----- Avoid unbalanced g-components numerical errors -----!

 if (izero==1.and.compute_nhat) then
!  Create fake mpi_enreg to wrap fourdp
   if (present(distribfft)) then
     my_distribfft => distribfft
   else
     ABI_DATATYPE_ALLOCATE(my_distribfft,)
     call init_distribfft_seq(my_distribfft,'f',ngfft(2),ngfft(3),'fourdp')
   end if
   call initmpi_seq(mpi_enreg_fft)
   ABI_DATATYPE_DEALLOCATE(mpi_enreg_fft%distribfft)
   if (present(comm_fft)) then
     call set_mpi_enreg_fft(mpi_enreg_fft,comm_fft,my_distribfft,me_g0,paral_kgb)
     my_comm_fft=comm_fft;paral_kgb_fft=paral_kgb
   else
     my_comm_fft=xmpi_comm_self;paral_kgb_fft=0;
     mpi_enreg_fft%distribfft => my_distribfft
   end if
!  Do FFT
   ABI_ALLOCATE(work,(2,nfft))
   cplex=2
   do isp1=1,MIN(2,nspinor**2)
     call fourdp(cplex,work,nhat12(:,:,isp1),-1,mpi_enreg_fft,nfft,1,ngfft,0)
     call zerosym(work,cplex,ngfft(1),ngfft(2),ngfft(3),comm_fft=my_comm_fft,distribfft=my_distribfft)
     call fourdp(cplex,work,nhat12(:,:,isp1),+1,mpi_enreg_fft,nfft,1,ngfft,0)
   end do
   ABI_DEALLOCATE(work)
!  Destroy fake mpi_enreg
   call unset_mpi_enreg_fft(mpi_enreg_fft)
   if (.not.present(distribfft)) then
     call destroy_distribfft(my_distribfft)
     ABI_DATATYPE_DEALLOCATE(my_distribfft)
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawmknhat_psipsi
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nhat/pawnhatfr
!!
!! NAME
!! pawnhatfr
!!
!! FUNCTION
!! PAW: Compute frozen part of 1st-order compensation charge density nhat^(1)
!!      nhatfr(r)=Sum_ij,lm[rhoij_ij.q_ij^l.(g_l(r).Y_lm(r))^(1)]
!!      Depends on q wave vector but not on first-order wave-function.
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

subroutine pawnhatfr(ider,idir,ipert,my_natom,natom,nspden,ntypat,&
&                    pawang,pawfgrtab,pawrhoij,pawtab,rprimd, &
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

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
     MSG_BUG('pawnhatfr: pawfgrtab()%rfgd array must be allocated!')
   end if
   if (pawrhoij(1)%qphase/=1) then
     MSG_BUG('pawnhatfr: not supposed to be called with qphase=2!')
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
         jrhoij=jrhoij+pawrhoij(iatom)%cplex_rhoij
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
         jrhoij=jrhoij+pawrhoij(iatom)%cplex_rhoij
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

!----------------------------------------------------------------------

!!****f* m_paw_nhat/pawsushat
!! NAME
!! pawsushat
!!
!! FUNCTION
!! PAW only, for susceptibility matrix:
!! Compute contribution to the product of two wavefunctions (exchange charge density)
!! from hat (compensation charge) density (in reciprocal space and eventually in real space):
!!    sushat_{ij,R}(g)=Sum_{L}[Q^L_ijR(g)]
!!
!! INPUTS
!!  atindx(natom)=index table for atoms, inverse of atindx
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!                                 WARNING: cprj(iatom,:) ARE SORTED BY ATOM TYPE !!!
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for small FFT sphere.
!!  gylmg_diel(npwdiel,lmax_diel**2,ntypat)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!  iband1,iband2= indices of the bands concerned with
!!  ispinor1,ispinor2= indices of spinorial components concerned with
!!  istwf_k=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FT components
!!  natom=number of atoms in cell
!!  nband=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  nfftdiel=number of FFT grid points for the small (diel) grid
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  optreal=0 if WF product has to be output in reciprocal space
!!          1 if WF product has to be output in real space
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!
!! SIDE EFFECTS
!!  === if optreal=0
!!  wfprod(2,npwdiel)=PAW contrib. to product of two wavefunctions (iband1,iband2):
!!                    is added (in reciprocal space)
!!  === if optreal=1
!!  wfraug(2,ndiel4,ndiel5,ndiel6)=PAW contrib. to product of two wavefunctions (iband1,iband2)
!!                                 is added (in real space)
!!
!! PARENTS
!!      susk,suskmm
!!
!! CHILDREN
!!      destroy_distribfft,fourwf,free_my_atmtab,get_my_atmtab
!!      init_distribfft_seq,initmpi_seq,set_mpi_enreg_fft,unset_mpi_enreg_fft
!!      xmpi_sum
!!
!! SOURCE

subroutine pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband1,iband2,ispinor1,ispinor2,istwf_k,kg_diel,&
&                    lmax_diel,mgfftdiel,natom,nband,ndiel4,ndiel5,ndiel6,&
&                    ngfftdiel,npwdiel,nspinor,ntypat,optreal,&
&                    pawang,pawtab,ph3d_diel,typat,wfprod,wfraug, &
&                    mpi_atmtab,comm_atom,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments (parallelism)

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iband1,iband2,ispinor1,ispinor2,istwf_k,lmax_diel,mgfftdiel
 integer,intent(in) :: natom,nband,ndiel4,ndiel5,ndiel6,npwdiel,nspinor
 integer,intent(in) :: ntypat,optreal
 integer,optional,intent(in) :: me_g0,comm_atom,comm_fft,paral_kgb
 type(distribfft_type),optional,intent(in),target :: distribfft
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: kg_diel(3,npwdiel),ngfftdiel(18),typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom)
 real(dp),intent(inout) :: wfprod(2,npwdiel*(1-optreal))
 real(dp),intent(inout) :: wfraug(2,ndiel4,ndiel5,ndiel6*optreal)
 type(pawcprj_type),intent(in) :: cprj_k(natom,nspinor*nband)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex,iatm,iatom,iatom_tot,ibsp1,ibsp2,ierr,il,ilmn,ils,ilslm,ipw
 integer :: itypat,j0lmn,jlmn,klm,klmn,lmax,lmin,mm,my_comm_atom,my_comm_fft,my_natom,tim_fourwf
 real(dp) :: phil1,phil2,sgn,weight_dum,wf1,wf2
 logical :: my_atmtab_allocated,parity,paral_atom
 type(distribfft_type),pointer :: my_distribfft
 type(mpi_type) :: mpi_enreg_fft
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: ro(2),ro_ql(2)
 real(dp),allocatable :: dummy(:,:),wfprod_paw(:,:),wfraug_paw(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (present(comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
     MSG_BUG('Need paral_kgb and me_g0 with comm_fft !')
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom)
 my_natom=natom;if (paral_atom) my_natom=size(my_atmtab)

 cplex=1;if (istwf_k>1) cplex=2
 ABI_ALLOCATE(wfprod_paw,(2,npwdiel))
 wfprod_paw(:,:)=zero
 ibsp1=(iband1-1)*nspinor+ispinor1
 ibsp2=(iband2-1)*nspinor+ispinor2

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   iatm=atindx(iatom_tot)
   itypat=typat(iatom_tot)

!  ------------------------------------------------------------------------
!  ----- Loop over ij channels (basis components)
!  ------------------------------------------------------------------------
   do jlmn=1,pawtab(itypat)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)
       lmax=pawtab(itypat)%indklmn(4,klmn)

       ro(1)=cprj_k(iatm,ibsp1)%cp(1,ilmn)*cprj_k(iatm,ibsp2)%cp(1,jlmn)
       if (cplex==2) then
         ro(1)=ro(1)+cprj_k(iatm,ibsp1)%cp(2,ilmn)*cprj_k(iatm,ibsp2)%cp(2,jlmn)
         ro(2)=cprj_k(iatm,ibsp1)%cp(2,ilmn)*cprj_k(iatm,ibsp2)%cp(1,jlmn) &
&         -cprj_k(iatm,ibsp1)%cp(1,ilmn)*cprj_k(iatm,ibsp2)%cp(2,jlmn)
       end if
       ro(1:cplex)=ro(1:cplex)*pawtab(itypat)%dltij(klmn)

       do ils=lmin,lmax,2
         il=mod(ils,4);parity=(mod(il,2)==0)
         sgn=one;if (il>1) sgn=-one

         do mm=-ils,ils
           ilslm=ils*ils+ils+mm+1
           if (pawang%gntselect(ilslm,klm)>0) then

             ro_ql(1:cplex)=pawtab(itypat)%qijl(ilslm,klmn)*ro(1:cplex)

!            Compute: Sum_{ijR} [ cpi* cpj qij^l (-i)^l g_l(g) S_lm(g) ]

             if (cplex==1) then
               if (parity) then
                 do ipw=1,npwdiel
                   phil1= sgn*ph3d_diel(1,ipw,iatm)     ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(2,ipw,iatm)
                   wf1= phil1*ro_ql(1)                  ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               else
                 do ipw=1,npwdiel
                   phil1=-sgn*ph3d_diel(2,ipw,iatm)  ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(1,ipw,iatm)
                   wf1= phil1*ro_ql(1)               ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               end if

             else

               if (parity) then
                 do ipw=1,npwdiel
                   phil1= sgn*ph3d_diel(1,ipw,iatm)     ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(2,ipw,iatm)
                   wf1=phil1*ro_ql(1)+phil2*ro_ql(2)    ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=phil1*ro_ql(2)-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               else
                 do ipw=1,npwdiel
                   phil1=-sgn*ph3d_diel(2,ipw,iatm)     ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(1,ipw,iatm)
                   wf1=phil1*ro_ql(1)+phil2*ro_ql(2)    ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=phil1*ro_ql(2)-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               end if

             end if
           end if
         end do
       end do

!      ----- End loop over ij channels
     end do
   end do

!  ----- End loop over atoms
 end do

!Reduction in case of parallelism over atoms
 if (paral_atom) then
   call xmpi_sum(wfprod_paw,my_comm_atom,ierr)
 end if

 if (optreal==0) then

!  === Output in reciprocal space
   wfprod(:,:)=wfprod(:,:)+wfprod_paw(:,:)

 else
!  === Output in reciprocal space
   tim_fourwf=17;weight_dum=0
!  Create fake mpi_enreg to wrap fourdp
   if (present(distribfft)) then
     my_distribfft => distribfft
   else
     ABI_DATATYPE_ALLOCATE(my_distribfft,)
     call init_distribfft_seq(my_distribfft,'c',ngfftdiel(2),ngfftdiel(3),'fourwf')
   end if
   call initmpi_seq(mpi_enreg_fft)
   ABI_DATATYPE_DEALLOCATE(mpi_enreg_fft%distribfft)
   if (present(comm_fft)) then
     call set_mpi_enreg_fft(mpi_enreg_fft,comm_fft,my_distribfft,me_g0,paral_kgb)
     my_comm_fft=comm_fft
     mpi_enreg_fft%paral_kgb = paral_kgb
   else
     my_comm_fft=xmpi_comm_self
     mpi_enreg_fft%paral_kgb = 0
     mpi_enreg_fft%distribfft => my_distribfft
   end if
!  do FFT
   ABI_ALLOCATE(wfraug_paw,(2,ndiel4,ndiel5,ndiel6))
   call fourwf(1,dummy,wfprod_paw,dummy,wfraug_paw,gbound_diel,gbound_diel,&
&   istwf_k,kg_diel,kg_diel,mgfftdiel,mpi_enreg_fft,1,ngfftdiel,1,npwdiel,&
&   ndiel4,ndiel5,ndiel6,0,tim_fourwf,weight_dum,weight_dum)
   wfraug(:,:,:,:)=wfraug(:,:,:,:)+wfraug_paw(:,:,:,:)
   ABI_DEALLOCATE(wfraug_paw)
   call unset_mpi_enreg_fft(mpi_enreg_fft)
   if (.not.present(distribfft)) then
     call destroy_distribfft(my_distribfft)
     ABI_DATATYPE_DEALLOCATE(my_distribfft)
   end if
 end if

 ABI_DEALLOCATE(wfprod_paw)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawsushat
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nhat/nhatgrid
!! NAME
!! nhatgrid
!!
!! FUNCTION
!! Determine parts of the rectangular (fine) grid that are contained
!! inside spheres around atoms (used to compute n_hat density).
!! If corresponding option is selected, compute also g_l(r)*Y_lm(r)
!! (and derivatives) on this grid (g_l=radial shape function).
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  distribfft<type(distribfft_type)>=--optional-- contains all the informations related
!!                                    to the FFT parallelism and plane sharing
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in unit cell
!!  optcut= option for the cut-off radius of spheres:
!!          if optcut=0, cut-off radius=pawtab%rshp=cut-off radius of compensation charge
!!          if optcut=1, cut-off radius=pawtab%rpaw=radius of PAW augmentation regions
!!  optgr0= 1 if g_l(r)*Y_lm(r) are computed
!!  optgr1= 1 if first derivatives of g_l(r)*Y_lm(r) are computed
!!  optgr2= 1 if second derivatives of g_l(r)*Y_lm(r) are computed
!!  optrad= 1 if vectors (r-r_atom) on the fine grid around atoms have to be stored
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type (integer) for each atom
!!  typord=1 if the output is ordered by type of atoms, 0 otherwise
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  pawfgrtab(natom)%ifftsph(nfgd)=FFT index (fine grid) of a points in paw spheres around each atom
!!  pawfgrtab(natom)%nfgd= number of (fine grid) FFT points in paw spheres around atoms
!!  if (optgr0==1)
!!    pawfgrtab(natom)%gylm(nfgd,l_size**2)= g_l(r)*Y_lm(r) around each atom
!!  if (optgr1==1)
!!    pawfgrtab(natom)%gylmgr(3,nfgd,l_size**2)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    pawfgrtab(natom)%gylmgr2(6,nfgd,l_size**2)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optrad==1)
!!    pawfgrtab(natom)%rfgd(3,nfgd)= coordinates of r-r_atom around each atom
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,classify_bands,denfgr,exc_plot,m_wfd
!!      pawmkaewf,respfn,scfcv,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawgylm,pawrfgd_fft,timab
!!
!! SOURCE

subroutine nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfft,ntypat,&
& optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,typat,ucvol,xred, &
& mpi_atmtab,comm_atom,comm_fft,distribfft,typord) ! optional arguments (parallelism)

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat,optcut,optgr0,optgr1,optgr2,optrad
 integer,optional,intent(in) :: comm_atom,comm_fft,typord
 real(dp),intent(in) :: ucvol
 type(distribfft_type),optional,target,intent(in)  :: distribfft
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 integer,intent(in),target :: atindx1(natom),nattyp(ntypat)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),xred(3,natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ------------------------------
!scalars
 integer :: i3,iat,iatm,iatom,iatom_,iatom_tot,itypat,lm_size,me_fft,my_comm_atom,n1,n2,n3,nfgd
 logical :: grid_found,my_atmtab_allocated,paral_atom
 real(dp) :: rcut
 character(len=500) :: msg
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 integer,pointer :: my_atindx1(:),my_atmtab(:),my_nattyp(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: rfgd_tmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(559,1,tsec)
 if (my_natom==0) return

!Set up parallelism over FFT
 me_fft=0
 if (present(comm_fft)) then
   me_fft=xmpi_comm_rank(comm_fft)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)
 if (paral_atom) then
   ABI_ALLOCATE(my_atindx1,(natom))
   ABI_ALLOCATE(my_nattyp,(ntypat))
   my_atindx1(:)=0;my_nattyp(:)=0
   iat=1
   do itypat=1,ntypat
     if (my_natom>0) then
       do iatom=1,my_natom
         if(typat(my_atmtab(iatom))==itypat)then
           my_nattyp(itypat)=my_nattyp(itypat)+1
           my_atindx1(iat)=iatom
           iat=iat+1
         end if
       end do
     end if
   end do
 else
   my_atindx1 => atindx1
   my_nattyp => nattyp
 end if

!Get the distrib associated with this fft_grid
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 if (present(distribfft)) then
   grid_found=.false.
   if (n2 == distribfft%n2_coarse) then
     if (n3== size(distribfft%tab_fftdp3_distrib)) then
       fftn3_distrib => distribfft%tab_fftdp3_distrib
       ffti3_local => distribfft%tab_fftdp3_local
       grid_found=.true.
     end if
   end if
   if (n2 == distribfft%n2_fine) then
     if (n3 == size(distribfft%tab_fftdp3dg_distrib)) then
       fftn3_distrib => distribfft%tab_fftdp3dg_distrib
       ffti3_local => distribfft%tab_fftdp3dg_local
       grid_found = .true.
     end if
   end if
   if (.not.(grid_found)) then
     msg='Unable to find an allocated distrib for this fft grid!'
     MSG_BUG(msg)
   end if
 else
   ABI_ALLOCATE(fftn3_distrib,(n3))
   ABI_ALLOCATE(ffti3_local,(n3))
   fftn3_distrib=0;ffti3_local=(/(i3,i3=1,n3)/)
 end if

!Loop over types of atom
!-------------------------------------------
 iatm=0
 do itypat=1,ntypat

   if (optcut==1) then
     rcut=pawtab(itypat)%rpaw
   else
     rcut=pawtab(itypat)%rshp
   end if

!  Loop over atoms
!  -------------------------------------------
   do iat=1,my_nattyp(itypat)
     iatm=iatm+1;iatom=my_atindx1(iatm)
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     iatom_=iatom;if(present(typord)) iatom_=merge(iatm,iatom,typord==1)
     lm_size=pawfgrtab(iatom_)%l_size**2

!    ------------------------------------------------------------------
!    A-Determine FFT points and r-R vectors around the atom
!    ------------------------------------------------------------------

     call pawrfgd_fft(ifftsph_tmp,gmet,n1,n2,n3,nfgd,rcut,rfgd_tmp,rprimd,ucvol,&
&     xred(:,iatom_tot),fft_distrib=fftn3_distrib,fft_index=ffti3_local,me_fft=me_fft)

!    Allocate arrays defining sphere (and related data) around current atom
     if (allocated(pawfgrtab(iatom_)%ifftsph)) then
       ABI_DEALLOCATE(pawfgrtab(iatom_)%ifftsph)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom_)%ifftsph,(nfgd))
     pawfgrtab(iatom_)%nfgd=nfgd
     pawfgrtab(iatom_)%ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)

     if (optrad==1) then
       if (allocated(pawfgrtab(iatom_)%rfgd))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%rfgd)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%rfgd,(3,nfgd))
       pawfgrtab(iatom_)%rfgd_allocated=1
       pawfgrtab(iatom_)%rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)
     end if

     if (optgr0==1) then
       if (allocated(pawfgrtab(iatom_)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%gylm,(nfgd,lm_size))
       pawfgrtab(iatom_)%gylm_allocated=1
     end if

     if (optgr1==1) then
       if (allocated(pawfgrtab(iatom_)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%gylmgr,(3,nfgd,lm_size))
       pawfgrtab(iatom_)%gylmgr_allocated=1
     end if

     if (optgr2==1) then
       if (allocated(pawfgrtab(iatom_)%gylmgr2))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%gylmgr2)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%gylmgr2,(6,nfgd,lm_size))
       pawfgrtab(iatom_)%gylmgr2_allocated=1
     end if

!    ------------------------------------------------------------------
!    B-Calculate g_l(r-R)*Y_lm(r-R) for each r around the atom R
!    ------------------------------------------------------------------
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom_)%gylm,pawfgrtab(iatom_)%gylmgr,pawfgrtab(iatom_)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),rfgd_tmp(:,1:nfgd))
     end if

!    End loops over types/atoms
!    -------------------------------------------
     ABI_DEALLOCATE(ifftsph_tmp)
     ABI_DEALLOCATE(rfgd_tmp)
   end do
 end do

!Destroy atom tables used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 if (paral_atom) then
   ABI_DEALLOCATE(my_atindx1)
   ABI_DEALLOCATE(my_nattyp)
 end if

 if (.not.present(distribfft)) then
   ABI_DEALLOCATE(fftn3_distrib)
   ABI_DEALLOCATE(ffti3_local)
 end if

 call timab(559,2,tsec)

 DBG_EXIT("COLL")

end subroutine nhatgrid
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nhat/wvl_nhatgrid
!! NAME
!! wvl_nhatgrid
!!
!! FUNCTION
!! Determine parts of the rectangular (fine) grid that are contained
!! inside spheres around atoms (used to compute n_hat density).
!! If corresponding option is selected, compute also g_l(r)*Y_lm(r)
!! (and derivatives) on this grid (g_l=radial shape function).
!!
!! INPUTS
!!
!! OUTPUT
!!  pawfgrtab(natom)%ifftsph(nfgd)=FFT index (fine grid) of a points in paw spheres around each atom
!!  pawfgrtab(natom)%nfgd= number of (fine grid) FFT points in paw spheres around atoms
!!  if (optgr0==1)
!!    pawfgrtab(natom)%gylm(nfgd,l_size**2)= g_l(r)*Y_lm(r) around each atom
!!  if (optgr1==1)
!!    pawfgrtab(natom)%gylmgr(3,nfgd,l_size**2)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    pawfgrtab(natom)%gylmgr2(6,nfgd,l_size**2)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optrad==1)
!!    pawfgrtab(natom)%rfgd(3,nfgd)= coordinates of r-r_atom around each atom
!!
!! NOTES
!!   PENDING: ADD PARALELLISM OVER ATOMS:
!!   COPY NHATGRID
!!
!! PARENTS
!!      afterscfloop,scfcv
!!
!! CHILDREN
!!      pawgylm,pawrfgd_wvl,timab,xred2xcart
!!
!! SOURCE

subroutine wvl_nhatgrid(atindx1,geocode,h,i3s,natom,natom_tot,&
& nattyp,ntypat,n1,n1i,n2,n2i,n3,n3pi,optcut,optgr0,optgr1,optgr2,optrad,&
& pawfgrtab,pawtab,psppar,rprimd,shift,xred)

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: i3s,natom,natom_tot,ntypat,optcut,optgr0,optgr1,optgr2,optrad
 integer,intent(in) :: n1,n2,n3,n1i,n2i,n3pi,shift
 real(dp),intent(in) :: h(3)
 character(1),intent(in) :: geocode
!integer,intent(in),optional :: mpi_comm_wvl
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 real(dp),intent(in) :: psppar(0:4,0:6,ntypat),rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ------------------------------
!scalars
!buffer to be added at the end of the last dimension of an array to control bounds_check
 integer :: iat,iatm,iatom,iatom_tot,itypat,lm_size,nfgd
 real(dp) :: rloc,rshp,xcart(3,natom)
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 real(dp) :: hh(3) !fine grid spacing for wavelets
 real(dp) :: tsec(2)
 real(dp),allocatable :: rfgd_tmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

#if !defined HAVE_BIGDFT
 BIGDFT_NOTENABLED_ERROR()
#endif

 call timab(559,1,tsec)

!Set up parallelism for wvl
!for debug: use me_wvl=xmpi_comm_rank(MPI_COMM_WORLD)
!if (present(mpi_comm_wvl)) then
!me_wvl=xmpi_comm_rank(mpi_comm_wvl)
!nproc_wvl=xmpi_comm_size(mpi_comm_wvl)
!else
!me_wvl=0;nproc_wvl=1
!end if
!Pending: parallelism over atoms: see nhatgrid

 if (natom_tot<natom) then   ! This test has to be remove when natom_tot is used
   MSG_BUG(' natom_tot<natom !')
 end if

!Fine grid
 hh(:)=0.5d0*h(:)

!Compute xcart from xred
 call xred2xcart(natom,rprimd,xcart,xred)

!Loop over types of atom
 iatm=0
 do itypat=1,ntypat

   rloc=psppar(0,0,itypat)
   if (optcut==1) then
     rshp=pawtab(itypat)%rpaw
   else
     rshp=pawtab(itypat)%rshp
   end if

!  Loop over atoms
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)
     iatom_tot=iatom; !if (paral_atom) iatom_tot=my_atmtab(iatom)
     lm_size=pawfgrtab(iatom)%l_size**2

!    Determine FFT points and r-R vectors around the atom
     call pawrfgd_wvl(geocode,hh,ifftsph_tmp,i3s,n1,n1i,n2,n2i,n3,n3pi,nfgd,rshp,rloc,&
&     rfgd_tmp,shift,xcart(:,iatom_tot))

!    Allocate arrays defining sphere (and related data) around current atom
     if (allocated(pawfgrtab(iatom)%ifftsph)) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%ifftsph)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%ifftsph,(nfgd))
     pawfgrtab(iatom)%nfgd=nfgd
     pawfgrtab(iatom)%ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)

     if (optrad==1) then
       if (allocated(pawfgrtab(iatom)%rfgd)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%rfgd)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%rfgd,(3,nfgd))
       pawfgrtab(iatom)%rfgd_allocated=1
       pawfgrtab(iatom)%rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)
     end if

     if (optgr0==1) then
       if (allocated(pawfgrtab(iatom)%gylm)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
       pawfgrtab(iatom)%gylm_allocated=1
     end if

     if (optgr1==1) then
       if (allocated(pawfgrtab(iatom)%gylmgr)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,lm_size))
       pawfgrtab(iatom)%gylmgr_allocated=1
     end if

     if (optgr2==1) then
       if (allocated(pawfgrtab(iatom)%gylmgr2)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(6,nfgd,lm_size))
       pawfgrtab(iatom)%gylmgr2_allocated=1
     end if

!    Calculate g_l(r-R)*Y_lm(r-R) for each r around the atom R
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),rfgd_tmp(:,1:nfgd))
     end if

!    End loops over types/atoms
     ABI_DEALLOCATE(ifftsph_tmp)
     ABI_DEALLOCATE(rfgd_tmp)
   end do
 end do

 call timab(559,2,tsec)

 DBG_EXIT("COLL")

end subroutine wvl_nhatgrid
!!***

!----------------------------------------------------------------------

END MODULE m_paw_nhat
!!***


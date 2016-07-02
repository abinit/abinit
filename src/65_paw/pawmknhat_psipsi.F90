!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmknhat_psipsi
!! NAME
!! pawmknhat_psipsi
!!
!! FUNCTION
!! PAW only:
!! Compute on the fine FFT grid the compensation charge density (and derivatives) associated
!! to the product of two wavefunctions n_{12}(r) = \Psi_1* \Psi_2. Based on pawmknhat.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MG, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat12_grdim= 0 if grnhat12 array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    nhat12(2,nfft,nspinor**2)=nhat on fine rectangular grid
!!  === if ider=1 or 2
!!    grnhat12(nfft,nspinor**2,3)=derivatives of nhat on fine rectangular grid (and derivatives)
!!  === if ider=3
!!    grnhat12(nfft,nspinor**2,3)=derivatives of nhat on fine rectangular grid (and derivatives)
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmknhat_psipsi(cprj1,cprj2,ider,izero,my_natom,natom,nfft,ngfft,nhat12_grdim,&
&          nspinor,ntypat,pawang,pawfgrtab,grnhat12,nhat12,pawtab, &
&          gprimd,grnhat_12,qphon,xred,mpi_atmtab,comm_atom,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use defs_abitypes,    only : mpi_type
 use m_mpinfo,         only : set_mpi_enreg_fft,unset_mpi_enreg_fft
 use m_lmn_indices,    only : klmn2ijlmn
 use m_pawang,         only : pawang_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgrtab,      only : pawfgrtab_type
 use m_pawcprj,        only : pawcprj_type
 use m_paw_finegrid,   only : pawgylm,pawexpiqr
 use m_paral_atom,     only : get_my_atmtab, free_my_atmtab
 use m_distribfft,     only : distribfft_type, init_distribfft_seq, destroy_distribfft

! use m_lmn_indices,  only : klmn2ijlmn

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmknhat_psipsi'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

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
 integer :: iatom,iatom_tot,ic,ierr,ils,ilslm,isp1,isp2,isploop,itypat,jc,klm,klmn
 integer :: lmax,lmin,lm_size,mm,my_comm_atom,my_comm_fft,optgr0,optgr1,paral_kgb_fft
 integer :: cplex,ilmn,jlmn,lmn_size,lmn2_size
 real(dp) :: re_p,im_p
 logical :: compute_grad,compute_grad1,compute_nhat,my_atmtab_allocated,paral_atom,qeq0,compute_phonon
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
     MSG_BUG('Need paral_kgb and me_g0 with comm_fft !')
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
!       call klmn2ijlmn(klmn,lmn_size,ilmn,jlmn)  ! This mapping should be stored in pawtab_type


!      Retrieve the factor due to the PAW projections.
       re_p =  cprj1(iatom_tot,isp1)%cp(1,ilmn) * cprj2(iatom_tot,isp2)%cp(1,jlmn) &
&       +cprj1(iatom_tot,isp1)%cp(2,ilmn) * cprj2(iatom_tot,isp2)%cp(2,jlmn) &
&       +cprj1(iatom_tot,isp1)%cp(1,jlmn) * cprj2(iatom_tot,isp2)%cp(1,ilmn) &
&       +cprj1(iatom_tot,isp1)%cp(2,jlmn) * cprj2(iatom_tot,isp2)%cp(2,ilmn)

       im_p =  cprj1(iatom_tot,isp1)%cp(1,ilmn) * cprj2(iatom_tot,isp2)%cp(2,jlmn) &
&       -cprj1(iatom_tot,isp1)%cp(2,ilmn) * cprj2(iatom_tot,isp2)%cp(1,jlmn) &
&       +cprj1(iatom_tot,isp1)%cp(1,jlmn) * cprj2(iatom_tot,isp2)%cp(2,ilmn) &
&       -cprj1(iatom_tot,isp1)%cp(2,jlmn) * cprj2(iatom_tot,isp2)%cp(1,ilmn)

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
               cpf_ql(1)=cpf(1)*pawtab(itypat)%qijl(ilslm,klmn)
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
               cpf_ql(1)=cpf(1)*pawtab(itypat)%qijl(ilslm,klmn)
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
       if(compute_phonon.and.pawfgrtab(iatom)%expiqr_allocated/=0) then
         do ic=1,pawfgrtab(iatom)%nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
           ro(1:2)=nhat12_atm(1:2,jc,isploop)
           nhat12_atm(1,jc,isploop)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           nhat12_atm(2,jc,isploop)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
     end if
     if (compute_grad) then
       if(compute_phonon.and.pawfgrtab(iatom)%expiqr_allocated/=0) then
         do ic=1,pawfgrtab(iatom)%nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
           ro(1)=grnhat12(1,jc,isploop,1)+qphon(1)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat12(2,jc,isploop,1)-qphon(1)*nhat12_atm(1,jc,isploop)
           grnhat12(1,jc,isploop,1)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat12(2,jc,isploop,1)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat12(1,jc,isploop,2)+qphon(2)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat12(2,jc,isploop,2)-qphon(2)*nhat12_atm(1,jc,isploop)
           grnhat12(1,jc,isploop,2)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat12(2,jc,isploop,2)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat12(1,jc,isploop,3)+qphon(3)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat12(2,jc,isploop,3)-qphon(3)*nhat12_atm(1,jc,isploop)
           grnhat12(1,jc,isploop,3)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat12(2,jc,isploop,3)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
     end if
     if (compute_grad1) then
       if(compute_phonon.and.pawfgrtab(iatom)%expiqr_allocated/=0) then
         do ic=1,pawfgrtab(iatom)%nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
           ro(1)=grnhat_12(1,jc,isploop,1,iatom)+qphon(1)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat_12(2,jc,isploop,1,iatom)-qphon(1)*nhat12_atm(1,jc,isploop)
           grnhat_12(1,jc,isploop,1,iatom)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat_12(2,jc,isploop,1,iatom)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat_12(1,jc,isploop,2,iatom)+qphon(2)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat_12(2,jc,isploop,2,iatom)-qphon(2)*nhat12_atm(1,jc,isploop)
           grnhat_12(1,jc,isploop,2,iatom)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           grnhat_12(2,jc,isploop,2,iatom)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
           ro(1)=grnhat_12(1,jc,isploop,3,iatom)+qphon(3)*nhat12_atm(2,jc,isploop)
           ro(2)=grnhat_12(2,jc,isploop,3,iatom)-qphon(3)*nhat12_atm(1,jc,isploop)
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
     call fourdp(cplex,work,nhat12(:,:,isp1),-1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
     call zerosym(work,cplex,ngfft(1),ngfft(2),ngfft(3),comm_fft=my_comm_fft,distribfft=my_distribfft)
     call fourdp(cplex,work,nhat12(:,:,isp1),+1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
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

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_grad_berry
!! NAME
!! make_grad_berry
!!
!! FUNCTION
!! compute gradient contribution from berry phase in finite
!! electric field case
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=input wavefunctions
!!  cgq(2,mcgq) = wavefunctions at neighboring k points
!!  cprj_k(natom,nband_k*usepaw)=cprj at this k point
!!  dimlmn(natom)=lmn_size for each atom in input order
!!  dimlmn_srt(natom)=lmn_size for each atom sorted by type
!!  direc(2,npw*nspinor)=gradient vector
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  iband=index of band currently being treated
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point currently being treated
!!  isppol=spin polarization currently treated
!!  natom=number of atoms in cell.
!!  mband =maximum number of bands
!!  mpw=maximum dimensioned size of npw
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array
!!  mkgq = second dimension of pwnsfacq
!!  nkpt=number of k points
!!  mpi_enreg=information about MPI parallelization
!!  npw=number of planewaves in basis sphere at given k.
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=number of spin polarizations
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point (electric field, MPI //)
!!
!! OUTPUT
!! grad_berry(2,npw*nspinor) :: contribution to gradient in finite electric field case
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = variables related to Berry phase calculations (see initberry.f)
!!
!! NOTES
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawcprj_symkn,smatrix,smatrix_k_paw
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_grad_berry(cg,cgq,cprj_k,detovc,dimlmn,dimlmn_srt,direc,dtefield,grad_berry,&
&                          gs_hamk,iband,icg,ikpt,isppol,mband,mcg,mcgq,mkgq,mpi_enreg,mpw,natom,nkpt,&
&                          npw,npwarr,nspinor,nsppol,pwind,pwind_alloc,pwnsfac,pwnsfacq)

 use defs_abitypes
 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_efield

 use m_pawcprj,     only : pawcprj_type, pawcprj_get, pawcprj_alloc, pawcprj_free, pawcprj_copy, pawcprj_symkn
 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_grad_berry'
 use interfaces_32_util
 use interfaces_65_paw
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!scalars
 integer,intent(in) :: iband,icg,ikpt,isppol,mband,mcg,mcgq
 integer,intent(in) :: mkgq,mpw,natom,nkpt,npw,nspinor,nsppol,pwind_alloc
 type(gs_hamiltonian_type),intent(in) :: gs_hamk
 type(efield_type),intent(inout) :: dtefield
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 integer,intent(in) :: dimlmn(natom),dimlmn_srt(natom)
 integer,intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cg(2,mcg),cgq(2,mcgq)
 real(dp),intent(inout) :: direc(2,npw*nspinor)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp),intent(out) :: detovc(2,2,3),grad_berry(2,npw*nspinor)
 type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%mband_occ*gs_hamk%usepaw*dtefield%nspinor)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,ddkflag,dimenlc1,dimenlr1,iatom,icg1,icp2,idum1
 integer :: idir,ifor,ikgf,ikptf,ikpt2,ikpt2f,ipw,i_paw_band,ispinor,itrs,itypat,job
 integer :: klmn,mcg1_k,mcg_q,nbo,npw_k2,nspinortot,paw_opt,shiftbd,signs
 real(dp) :: fac
 character(len=500) :: message
!arrays
 integer :: pwind_k(npw),sflag_k(dtefield%mband_occ)
 real(dp) :: cg1_k(2,npw*nspinor),dtm_k(2),pwnsfac_k(4,mpw)
 real(dp) :: smat_k(2,dtefield%mband_occ,dtefield%mband_occ)
 real(dp) :: smat_inv(2,dtefield%mband_occ,dtefield%mband_occ),svectout_dum(2,0)
 real(dp) :: dummy_enlout(0)
 real(dp),allocatable :: cgq_k(:,:),enl_rij(:,:,:),grad_berry_ev(:,:)
 real(dp),allocatable :: qijbkk(:,:,:),smat_k_paw(:,:,:)
! type(pawcprj_type) :: cprj_dum(1,1) ! was used in on-site dipole, now suppressed
! 15 June 2012 J Zwanziger
 type(pawcprj_type),allocatable :: cprj_kb(:,:),cprj_band_srt(:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)


! *********************************************************************

!DBG_ENTER("COLL")

 nbo = dtefield%mband_occ

!allocations

!Electric field: compute the gradient of the Berry phase part of the energy functional.
!See PRL 89, 117602 (2002), grad_berry(:,:) is the second term of Eq. (4)
 grad_berry(:,:) = zero
 job = 11 ; shiftbd = 1
 mcg_q = mpw*mband*nspinor
 mcg1_k = npw*nspinor

 if (gs_hamk%usepaw /= 0) then
   dimenlr1 = gs_hamk%lmnmax*(gs_hamk%lmnmax+1)/2
   dimenlc1 = 2*dimenlr1
   ABI_ALLOCATE(qijbkk,(dimenlc1,natom,nspinor**2))
   ABI_ALLOCATE(enl_rij,(nspinor*dimenlr1,natom,nspinor**2))
   ABI_ALLOCATE(smat_k_paw,(2,nbo,nbo))
   ABI_ALLOCATE(grad_berry_ev,(2,npw*nspinor))
   enl_rij = zero
   qijbkk = zero
   smat_k_paw = zero
   ABI_DATATYPE_ALLOCATE(cprj_kb,(natom,nbo*nspinor))
   call pawcprj_alloc(cprj_kb,0,dimlmn)
   ABI_DATATYPE_ALLOCATE(cprj_band_srt,(natom,nspinor))
   call pawcprj_alloc(cprj_band_srt,0,dimlmn_srt)
   if (nkpt /= dtefield%fnkpt) then
     ABI_DATATYPE_ALLOCATE(cprj_fkn,(natom,nbo*nspinor))
     ABI_DATATYPE_ALLOCATE(cprj_ikn,(natom,nbo*nspinor))
     call pawcprj_alloc(cprj_fkn,0,dimlmn)
     call pawcprj_alloc(cprj_ikn,0,dimlmn)
   else
     ABI_DATATYPE_ALLOCATE(cprj_fkn,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_ikn,(0,0))
   end if
 else
   ABI_ALLOCATE(qijbkk,(0,0,0))
   ABI_ALLOCATE(enl_rij,(0,0,0))
   ABI_ALLOCATE(smat_k_paw,(0,0,0))
   ABI_ALLOCATE(grad_berry_ev,(0,0))
   ABI_DATATYPE_ALLOCATE(cprj_kb,(0,0))
   ABI_DATATYPE_ALLOCATE(cprj_band_srt,(0,0))
   ABI_DATATYPE_ALLOCATE(cprj_fkn,(0,0))
   ABI_DATATYPE_ALLOCATE(cprj_ikn,(0,0))
 end if

 ikptf = dtefield%i2fbz(ikpt)
 ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind

 do idir = 1, 3
!  skip idir values for which efield_dot(idir)=0
   if (abs(dtefield%efield_dot(idir)) < tol12) cycle
!  Implicitly, we use the gradient multiplied by the number of k points in the FBZ
   fac = dtefield%efield_dot(idir)*dble(dtefield%fnkpt)/&
&   (dble(dtefield%nstr(idir))*four_pi)
   do ifor = 1, 2
!    Handle dtefield%i2fbz properly and ask whether t.r.s. is used
     ikpt2f = dtefield%ikpt_dk(ikptf,ifor,idir)
     if (dtefield%indkk_f2ibz(ikpt2f,6) == 1) then
       itrs = 10
     else
       itrs = 0
     end if
     ikpt2 = dtefield%indkk_f2ibz(ikpt2f,1)
     npw_k2 = npwarr(ikpt2)
     ABI_ALLOCATE(cgq_k,(2,nbo*nspinor*npw_k2))
     pwind_k(1:npw) = pwind(ikgf+1:ikgf+npw,ifor,idir)
     pwnsfac_k(1:2,1:npw) = pwnsfac(1:2,ikgf+1:ikgf+npw)
     sflag_k(:) = dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir)
     smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir)
     if (mpi_enreg%nproc_cell > 1) then
       icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
       cgq_k(:,1:nbo*nspinor*npw_k2) = &
&       cgq(:,icg1+1:icg1+nbo*nspinor*npw_k2)
       idum1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
     else
       icg1 = dtefield%cgindex(ikpt2,isppol)
       cgq_k(:,1:nbo*nspinor*npw_k2) = &
&       cg(:,icg1+1:icg1+nbo*nspinor*npw_k2)
       idum1 = dtefield%fkgindex(ikpt2f)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
     end if
     if (gs_hamk%usepaw == 1) then
       icp2=nbo*(ikpt2-1)*nspinor
       call pawcprj_get(gs_hamk%atindx1,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,isppol,&
&       nbo,dtefield%fnkpt,natom,nbo,nbo,nspinor,nsppol,0,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
         call pawcprj_copy(cprj_kb,cprj_ikn)
         call pawcprj_symkn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&         dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
&         dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
&         dtefield%lmax,dtefield%lmnmax,mband,natom,nbo,nspinor,&
&         dtefield%nsym,gs_hamk%ntypat,gs_hamk%typat,dtefield%zarot)
         call pawcprj_copy(cprj_fkn,cprj_kb)
       end if
       call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,mband,natom,smat_k_paw,gs_hamk%typat)
     end if

     icg1 = 0 ; ddkflag = 1
     call smatrix(cg,cgq_k,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,&
&     job,iband,mcg,mcg_q,mcg1_k,iband,mpw,nbo,dtefield%nband_occ(isppol),&
&     npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&     shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
     ABI_DEALLOCATE(cgq_k)
     detovc(:,ifor,idir) = dtm_k(:) !store the determinant of the overlap
     if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
       write(message,'(3a,i5,a,i3,a,a,a)') &
&       '  (electric field)',ch10,&
&       '  For k-point #',ikpt,' and band # ',iband,',',ch10,&
&       '  the determinant of the overlap matrix is found to be 0. Fixing...'
!      REC try this:
       write(std_out,*)message,dtm_k(1:2)
       if(abs(dtm_k(1))<=1d-12)dtm_k(1)=1d-12
       if(abs(dtm_k(2))<=1d-12)dtm_k(2)=1d-12
       write(std_out,*)' Changing to:',dtm_k(1:2)
!      REC       MSG_BUG(message)
     end if

     if (gs_hamk%usepaw == 1) then
!      this loop applies discretized derivative of projectors
!      note that qijb_kk is sorted by input atom order, but nonlop wants it sorted by type
       do iatom = 1, natom
         itypat = gs_hamk%typat(gs_hamk%atindx1(iatom))
         do klmn = 1, dtefield%lmn2_size(itypat)
!          note: D_ij-like terms have 4 spinor components: 11, 22, 12, and 21. Here the qijb is diagonal
!          in spin space so only the first two are nonzero and they are equal
           do ispinor = 1, nspinor
             qijbkk(2*klmn-1,iatom,ispinor) = dtefield%qijb_kk(1,klmn,gs_hamk%atindx1(iatom),idir)
             qijbkk(2*klmn,  iatom,ispinor) = dtefield%qijb_kk(2,klmn,gs_hamk%atindx1(iatom),idir)
             if (ifor > 1) qijbkk(2*klmn,iatom,ispinor) = -qijbkk(2*klmn,iatom,ispinor)
           end do
         end do ! end loop over lmn2_size
       end do ! end loop over natom

       choice = 1
       signs = 2
       paw_opt = 1
       cpopt = 2 ! use cprj_kb in memory
       nspinortot=min(2,nspinor*(1+mpi_enreg%paral_spinor))
       do i_paw_band = 1, nbo

         call pawcprj_get(gs_hamk%atindx,cprj_band_srt,cprj_kb,natom,i_paw_band,0,ikpt,1,&
&         isppol,nbo,1,natom,1,nbo,nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

         ! Pass dummy_enlout to avoid aliasing (enl, enlout)
         call nonlop(choice,cpopt,cprj_band_srt,dummy_enlout,gs_hamk,idir,(/zero/),mpi_enreg,1,0,&
&         paw_opt,signs,svectout_dum,0,direc,grad_berry_ev,enl=qijbkk)

!        Add i*fac*smat_inv(i_paw_band,iband)*grad_berry_ev to the gradient
         do ipw = 1, npw*nspinor

           grad_berry(1,ipw) = grad_berry(1,ipw) - &
&           fac*(smat_inv(2,i_paw_band,iband)*grad_berry_ev(1,ipw) + &
&           smat_inv(1,i_paw_band,iband)*grad_berry_ev(2,ipw))

           grad_berry(2,ipw) = grad_berry(2,ipw) + &
&           fac*(smat_inv(1,i_paw_band,iband)*grad_berry_ev(1,ipw) - &
&           smat_inv(2,i_paw_band,iband)*grad_berry_ev(2,ipw))

         end do
       end do
     end if ! end if PAW

!    Add i*fac*cg1_k to the gradient
     do ipw = 1, npw*nspinor
       grad_berry(1,ipw) = grad_berry(1,ipw) - fac*cg1_k(2,ipw)
       grad_berry(2,ipw) = grad_berry(2,ipw) + fac*cg1_k(1,ipw)
     end do
     fac = -1._dp*fac
     dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir) = sflag_k(:)
     dtefield%sflag(iband,ikpt+(isppol-1)*nkpt,ifor,idir) = 0
     dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir) = smat_k(:,:,:)
   end do  ! ifor

!  if (gs_hamk%usepaw == 1) then
!  !    call nonlop to apply on-site dipole <EV> part to direc
!  !    note that rij is sorted by input atom order, but nonlop wants it sorted by type
!  do iatom = 1, natom
!  itypat = gs_hamk%typat(gs_hamk%atindx1(iatom))
!  do klmn = 1, dtefield%lmn2_size(itypat)
!  !        note: D_ij-like terms have 4 spinor components: 11, 22, 12, and 21. Here the enl_rij is diagonal
!  !        in spin space so only the first two are nonzero and they are equal
!  do ispinor = 1, nspinor
!  if (nspinor == 1) then
!  enl_rij(klmn,iatom,ispinor) = dtefield%rij(klmn,itypat,idir)
!  else
!  enl_rij(2*klmn-1,iatom,ispinor) = dtefield%rij(klmn,itypat,idir)
!  end if
!  end do
!  end do ! end loop over lmn2_size
!  end do ! end loop over natom
!  cpopt = -1 ! compute cprj inside nonlop because we do not have them for direc
!  call nonlop(choice,cpopt,cprj_dum,dummy_enlout,gs_hamk,idir,zero,mpi_enreg,1,0,&
!  &           paw_opt,signs,svectout_dum,0,direc,grad_berry_ev,enl=enl_rij)
!  grad_berry(:,:) = grad_berry(:,:) - dtefield%efield_dot(idir)*grad_berry_ev(:,:)/two_pi
!  end if

 end do    ! idir

!deallocations
 if(gs_hamk%usepaw /= 0) then
   call pawcprj_free(cprj_kb)
   call pawcprj_free(cprj_band_srt)
   if (nkpt /= dtefield%fnkpt) then
     call pawcprj_free(cprj_fkn)
     call pawcprj_free(cprj_ikn)
   end if
 end if
 ABI_DEALLOCATE(grad_berry_ev)
 ABI_DEALLOCATE(qijbkk)
 ABI_DEALLOCATE(enl_rij)
 ABI_DEALLOCATE(smat_k_paw)
 ABI_DATATYPE_DEALLOCATE(cprj_kb)
 ABI_DATATYPE_DEALLOCATE(cprj_band_srt)
 ABI_DATATYPE_DEALLOCATE(cprj_fkn)
 ABI_DATATYPE_DEALLOCATE(cprj_ikn)

!DBG_EXIT("COLL")

end subroutine make_grad_berry
!!***

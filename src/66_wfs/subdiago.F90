!{\src2tex{textfont=tt}}
!!****f* ABINIT/subdiago
!! NAME
!! subdiago
!!
!! FUNCTION
!! This routine diagonalizes the Hamiltonian in the eigenfunction subspace
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  istwf_k=input parameter that describes the storage of wfs
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  subham(nband_k*(nband_k+1))=Hamiltonian expressed in the WFs subspace
!!  subovl(nband_k*(nband_k+1)*use_subovl)=overlap matrix expressed in the WFs subspace
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  me_g0=1 if this processors has G=0, 0 otherwise.
!!
!! OUTPUT
!!  eig_k(nband_k)=array for holding eigenvalues (hartree)
!!  evec(2*nband_k,nband_k)=array for holding eigenvectors
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=wavefunctions
!!  gsc(2,mgsc)=<g|S|c> matrix elements (S=overlap)
!!
!! PARENTS
!!      rayleigh_ritz,vtowfk
!!
!! CHILDREN
!!      abi_xcopy,abi_xgemm,abi_xhpev,abi_xhpgv,cg_normev,hermit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine subdiago(cg,eig_k,evec,gsc,icg,igsc,istwf_k,&
&                   mcg,mgsc,nband_k,npw_k,nspinor,paral_kgb,&
&                   subham,subovl,use_subovl,usepaw,me_g0)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_cgtools
 use m_linalg_interfaces
 use m_abi_linalg
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'subdiago'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: icg,igsc,istwf_k,mcg,mgsc,nband_k,npw_k,me_g0
 integer,intent(in) :: nspinor,paral_kgb,use_subovl,usepaw
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1)),subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp),intent(out) :: eig_k(nband_k),evec(2*nband_k,nband_k)
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc)

!Local variables-------------------------------
 integer :: iband,ii,ierr,rvectsize,vectsize,use_slk
 character(len=500) :: message
 ! real(dp) :: tsec(2)
 real(dp),allocatable :: evec_tmp(:,:),subovl_tmp(:),subham_tmp(:)
 real(dp),allocatable :: work(:,:)
 real(dp),allocatable :: blockvectora(:,:),blockvectorb(:,:),blockvectorc(:,:)

! *********************************************************************

 if (paral_kgb<0) then
   MSG_BUG('paral_kgb should be positive ')
 end if

 ! 1 if Scalapack version is used.
 use_slk = paral_kgb 

 rvectsize=npw_k*nspinor
 vectsize=2*rvectsize;if (me_g0==1) vectsize=vectsize-1

!Impose Hermiticity on diagonal elements of subham (and subovl, if needed)
! MG FIXME: In these two calls we are aliasing the args
 call hermit(subham,subham,ierr,nband_k)
 if (use_subovl==1) then
   call hermit(subovl,subovl,ierr,nband_k)
 end if

!Diagonalize the Hamitonian matrix
 if(istwf_k==2) then
   ABI_ALLOCATE(evec_tmp,(nband_k,nband_k))
   ABI_ALLOCATE(subham_tmp,(nband_k*(nband_k+1)/2))
   subham_tmp=subham(1:nband_k*(nband_k+1):2)
   evec_tmp=zero
   if (use_subovl==1) then
     ABI_ALLOCATE(subovl_tmp,(nband_k*(nband_k+1)/2))
     subovl_tmp=subovl(1:nband_k*(nband_k+1):2)
!    TO DO: Not sure this one has been fully tested
     call abi_xhpgv(1,'V','U',nband_k, subham_tmp,subovl_tmp, eig_k,evec_tmp,istwf_k=istwf_k,use_slk=use_slk)
     ABI_DEALLOCATE(subovl_tmp)
   else
     call abi_xhpev('V','U',nband_k,subham_tmp,eig_k,evec_tmp,istwf_k=istwf_k,use_slk=use_slk)
   end if
   evec(:,:)=zero;evec(1:2*nband_k:2,:) =evec_tmp
   ABI_DEALLOCATE(evec_tmp)
   ABI_DEALLOCATE(subham_tmp)
 else
   if (use_subovl==1) then
     call abi_xhpgv(1,'V','U',nband_k,subham,subovl,eig_k,evec,istwf_k=istwf_k,use_slk=use_slk)
   else
     call abi_xhpev('V','U',nband_k,subham,eig_k,evec,istwf_k=istwf_k,use_slk=use_slk)
   end if
 end if

!Normalize each eigenvector and set phase:
!The problem with minus/plus signs might be present also if .not. use_subovl
!if(use_subovl == 0) then
 call cg_normev(evec,nband_k,nband_k)
!end if

 if(istwf_k==2)then
   do iband=1,nband_k
     do ii=1,nband_k
       if(abs(evec(2*ii,iband))>1.0d-10)then
         write(message,'(3a,2i0,2es16.6,a,a)')ch10,&
&         ' subdiago: For istwf_k=2, observed the following element of evec :',ch10,&
&         iband,ii,evec(2*ii-1,iband),evec(2*ii,iband),ch10,'  with a non-negligible imaginary part.'
         MSG_BUG(message)
       end if
     end do
   end do
 end if
 
!=====================================================
!Carry out rotation of bands C(G,n) according to evecs
! ZGEMM if istwfk==1, DGEMM if istwfk==2
!=====================================================
 if (istwf_k==2) then

   ABI_STAT_ALLOCATE(blockvectora,(vectsize,nband_k), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in blockvectora")
   ABI_STAT_ALLOCATE(blockvectorb,(nband_k,nband_k), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in blockvectorb")
   ABI_STAT_ALLOCATE(blockvectorc,(vectsize,nband_k), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in blockvectorc")

   do iband=1,nband_k
     if (me_g0 == 1) then       
       call abi_xcopy(1,cg(1,cgindex_subd(iband)),1,blockvectora(1,iband),1)
       call abi_xcopy(rvectsize-1,cg(1,cgindex_subd(iband)+1),2,blockvectora(2,iband),1)
       call abi_xcopy(rvectsize-1,cg(2,cgindex_subd(iband)+1),2,blockvectora(rvectsize+1,iband),1)
     else
       call abi_xcopy(rvectsize,cg(1,cgindex_subd(iband)),2,blockvectora(1,iband),1)
       call abi_xcopy(rvectsize,cg(2,cgindex_subd(iband)),2,blockvectora(rvectsize+1,iband),1)
     end if
     call abi_xcopy(nband_k,evec(2*iband-1,1),2*nband_k,blockvectorb(iband,1),nband_k)
   end do

   call abi_xgemm('N','N',vectsize,nband_k,nband_k,&
&   cone,blockvectora,vectsize,blockvectorb,nband_k,czero,blockvectorc,vectsize)

   do iband=1,nband_k
     if (me_g0 == 1) then
       call abi_xcopy(1,blockvectorc(1,iband),1,cg(1,cgindex_subd(iband)),1)
       call abi_xcopy(rvectsize-1,blockvectorc(2,iband),1,cg(1,cgindex_subd(iband)+1),2)
       call abi_xcopy(rvectsize-1,blockvectorc(rvectsize+1,iband),1,cg(2,cgindex_subd(iband)+1),2)
     else
       call abi_xcopy(rvectsize,blockvectorc(1,iband),1,cg(1,cgindex_subd(iband)),2)
       call abi_xcopy(rvectsize,blockvectorc(rvectsize+1,iband),1,cg(2,cgindex_subd(iband)),2)
     end if
   end do

!  If paw, musb also rotate S.C(G,n):
   if (usepaw==1) then

     do iband=1,nband_k
       if (me_g0 == 1) then
         call abi_xcopy(1,gsc(1,gscindex_subd(iband)),1,blockvectora(1,iband),1)
         call abi_xcopy(rvectsize-1,gsc(1,gscindex_subd(iband)+1),2,blockvectora(2,iband),1)
         call abi_xcopy(rvectsize-1,gsc(2,gscindex_subd(iband)+1),2,blockvectora(rvectsize+1,iband),1)
       else
         call abi_xcopy(rvectsize  ,gsc(1,gscindex_subd(iband)),2,blockvectora(1,iband),1)
         call abi_xcopy(rvectsize  ,gsc(2,gscindex_subd(iband)),2,blockvectora(rvectsize+1,iband),1)
       end if
       call abi_xcopy(nband_k,evec(2*iband-1,1),2*nband_k,blockvectorb(iband,1),nband_k)
     end do

     call abi_xgemm('N','N',vectsize,nband_k,nband_k,&
&     cone,blockvectora,vectsize,blockvectorb,nband_k,czero,blockvectorc,vectsize)

     do iband=1,nband_k
       if (me_g0 == 1) then
         call abi_xcopy(1,blockvectorc(1,iband),1,gsc(1,gscindex_subd(iband)),1)
         call abi_xcopy(rvectsize-1,blockvectorc(2,iband),1,gsc(1,gscindex_subd(iband)+1),2)
         call abi_xcopy(rvectsize-1,blockvectorc(rvectsize+1,iband),1,gsc(2,gscindex_subd(iband)+1),2)
       else
         call abi_xcopy(rvectsize,blockvectorc(1,iband),1,gsc(1,gscindex_subd(iband)),2)
         call abi_xcopy(rvectsize,blockvectorc(rvectsize+1,iband),1,gsc(2,gscindex_subd(iband)),2)
       end if
     end do

   end if

   ABI_DEALLOCATE(blockvectora)
   ABI_DEALLOCATE(blockvectorb)
   ABI_DEALLOCATE(blockvectorc)

 else

   ABI_STAT_ALLOCATE(work,(2,npw_k*nspinor*nband_k), ierr)
   ABI_CHECK(ierr==0, "out-of-memory in work")

!  MG: Do not remove this initialization.
!  telast_06 stops in fxphase on inca_debug and little_buda (very very strange, due to atlas?)
   work=zero

   call abi_xgemm('N','N',npw_k*nspinor,nband_k,nband_k,cone, &
&   cg(:,icg+1:npw_k*nspinor*nband_k+icg),npw_k*nspinor, &
&   evec,nband_k,czero,work,npw_k*nspinor,x_cplx=2)

   call abi_xcopy(npw_k*nspinor*nband_k,work(1,1),1,cg(1,1+icg),1,x_cplx=2)
   
!  If paw, must also rotate S.C(G,n):
   if (usepaw==1) then
     call abi_xgemm('N','N',npw_k*nspinor,nband_k,nband_k,cone, &
&     gsc(:,1+igsc:npw_k*nspinor*nband_k+igsc),npw_k*nspinor, &
&     evec,nband_k,czero,work,npw_k*nspinor,x_cplx=2)
     call abi_xcopy(npw_k*nspinor*nband_k, work(1,1),1,gsc(1,1+igsc),1,x_cplx=2)   
   end if

   ABI_DEALLOCATE(work)
 end if

 contains

   function cgindex_subd(iband)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgindex_subd'
!End of the abilint section

   integer :: iband,cgindex_subd
   cgindex_subd=npw_k*nspinor*(iband-1)+icg+1
 end function cgindex_subd
   function gscindex_subd(iband)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gscindex_subd'
!End of the abilint section

   integer :: iband,gscindex_subd
   gscindex_subd=npw_k*nspinor*(iband-1)+igsc+1
 end function gscindex_subd

end subroutine subdiago
!!***

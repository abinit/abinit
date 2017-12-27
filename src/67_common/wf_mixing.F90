!{\src2tex{textfont=tt}}
!!****f* ABINIT/wf_mixing
!!
!! NAME
!! wf_mixing
!!
!! FUNCTION
!! Mixing of wavefunctions in the outer loop of a double loop SCF approach.
!! Different algorithms are implemented, depending on the value of wfmixalg.
!!
!! COPYRIGHT
!! Copyright (C) 2017 ABINIT group (XG,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  istep=number of call the routine (usually the outer loop in the SCF double loop)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of cprj array 
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(dtset%ntypat)=number of atoms of each type in cell.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wf_mixing(atindx1,cg,cprj,dtset,istep,mcg,mcprj,mpi_enreg,&
& nattyp,npwarr,pawtab,scf_history)

 use defs_basis
 use defs_abitypes
 use m_scf_history
 use m_xmpi
 use m_profiling_abi
 use m_errors
 use m_cgtools

 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_get, pawcprj_lincom, &
&                      pawcprj_free, pawcprj_axpby, pawcprj_put, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wf_mixing'
 use interfaces_32_util
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: hermitian
 integer :: ibd,ibg,iblockbd,iblockbd1,icg,icgb
 integer :: ierr,ii,ikpt,inc,indh,inplace
 integer :: isize,isppol,istwf_k,kk,me_distrb,my_nspinor
 integer :: nband_k,nbdmix,npw_k,ntypat,ortalgo,spaceComm_band,usepaw,wfmixalg
 !character(len=500) :: message
!arrays
 real(dp) :: alpha,beta
 integer,allocatable :: bufsize(:),bufsize_wf(:),bufdisp(:),bufdisp_wf(:)
 integer,allocatable :: ipiv(:),dimcprj(:),npw_block(:),npw_disp(:)
 real(dp),allocatable :: al(:,:),cwavef(:,:),cwavefh(:,:)
 real(dp),allocatable :: dum(:,:)
 real(dp),allocatable :: dmn(:,:,:),dmn_debug(:,:,:),mmn(:,:,:)
 real(dp),allocatable :: smn(:,:,:),smn_(:,:,:)
 real(dp),allocatable :: work(:,:),work1(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kh(:,:),cprj_k3(:,:)
!DEBUG
!real(dp),allocatable :: cg_ref(:,:)
!type(pawcprj_type),allocatable :: cprj_ref(:,:)
!ENDDEBUG

! *************************************************************************

!DEBUG
!write(std_out,*)' wf_mixing : enter '
!write(std_out,*)' istep,scf_history%alpha=',istep,scf_history%alpha
!write(std_out,*)' cg(1:2,1:2)=',cg(1:2,1:2)
!write(std_out,*)' scf_history%cg(1:2,1:2,1)=',scf_history%cg(1:2,1:2,1)
!ABI_ALLOCATE(cg_ref,(2,mcg))
!cg_ref(:,:)=cg(:,:)
!ABI_DATATYPE_ALLOCATE(cprj_ref,(dtset%natom,mcprj))
!cprj_ref(:,:)=cprj(:,:)
!ENDDEBUG

 if (istep==0) return

 ntypat=dtset%ntypat
 usepaw=dtset%usepaw
 wfmixalg=scf_history%wfmixalg

!Useful array
 ABI_ALLOCATE(dimcprj,(dtset%natom))
 if (usepaw==1) then
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

!Index for the wavefunction stored in scf_history
 indh=1

!DEBUG
! if(istep==1)then
!   scf_history%cg(:,:,1)=cg(:,:)
!   if(usepaw==1) then
!     scf_history%cprj(:,:,1)=cprj(:,:)
!   end if
! else
!   cg(:,:)=scf_history%cg(:,:,1)
!   if(usepaw==1) then
!     cprj(:,:)=scf_history%cprj(:,:,1)
!   end if
! endif
! return
!ENDDEBUG

!First step
 if (istep==1 .or. (wfmixalg==2 .and. abs(scf_history%alpha-one)<tol8) ) then
   scf_history%cg(:,:,1)=cg(:,:)
   if(usepaw==1) then
     scf_history%cprj(:,:,1)=cprj(:,:)
   end if

 else
!From 2nd step

!  Init parallelism
   me_distrb=mpi_enreg%me_kpt
   spaceComm_band=xmpi_comm_self

   icg=0
   ibg=0

!  LOOP OVER SPINS
   my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       nbdmix=dtset%nbandhf
!      nbdmix=nband_k   ! This gives back the previous behavior, if needed.
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)

!      Allocate arrays for a wave-function (or a block of WFs)
       ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
       ABI_ALLOCATE(cwavefh,(2,npw_k*my_nspinor))

!      Space biorthogonalization

!      Loop over bands 
       icgb=icg

       ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,my_nspinor*nbdmix))
       ABI_DATATYPE_ALLOCATE(cprj_kh,(dtset%natom,my_nspinor*nbdmix))
       if(usepaw==1) then
         call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nbdmix,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_alloc(cprj_kh,scf_history%cprj(1,1,indh)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_kh,scf_history%cprj(:,:,indh),dtset%natom,1,ibg,ikpt,1,isppol,&
&         dtset%mband,dtset%mkmem,dtset%natom,nbdmix,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       ABI_ALLOCATE(smn,(2,nbdmix,nbdmix))

!DEBUG
!      write(std_out,*)' Compute the S matrix, whose matrix elements are scalar products.'
!ENDDEBUG
       hermitian=0
!      Note that cprj_k might have been used instead of cprj , in which case the first ibg should have been 0
!      Note that cprj_kh might have been used instead of scf_history%cprj(:,:,indh)a , in which case the second ibg should have been 0
       call dotprod_set_cgcprj(atindx1,cg,scf_history%cg(:,:,indh),cprj,scf_history%cprj(:,:,indh),dimcprj,hermitian,&
&        ibg,ibg,icg,icg,ikpt,isppol,istwf_k,dtset%mband,mcg,mcg,mcprj,mcprj,dtset%mkmem,&
&        mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn,usepaw)


!      Invert S matrix, that is NOT hermitian. 
!      Calculate M=S^-1
       ABI_ALLOCATE(mmn,(2,nbdmix,nbdmix))
       mmn=zero
       do kk=1,nbdmix
         mmn(1,kk,kk)=one
       end do

       ABI_ALLOCATE(smn_,(2,nbdmix,nbdmix))
       ABI_ALLOCATE(ipiv,(nbdmix))
!      For debugging purposes, the smn matrix is preserved (later, M * S = 1 is checked ...)
!      The smn_ arrays stores a copy of the smn array, and will be destroyed by the following inverse call
       smn_=smn
       call zgesv(nbdmix,nbdmix,smn_,nbdmix,ipiv,mmn,nbdmix,ierr)
       ABI_DEALLOCATE(ipiv)
       ABI_DEALLOCATE(smn_)
!DEBUG
       if(ierr/=0)then
         write(std_out,*)' wf_mixing : the call to cgesv general inversion routine returned an error code ierr=',ierr
         stop
       endif
!ENDDEBUG

       inplace=1
       call lincom_cgcprj(mmn,scf_history%cg(:,:,indh),cprj_kh,dimcprj,&
&         icg,inplace,mcg,my_nspinor*nband_k,dtset%natom,nbdmix,nbdmix,npw_k,my_nspinor,usepaw)

       ABI_DEALLOCATE(mmn)

!      Wavefunction extrapolation, simple mixing case
       ibd=0  
       inc=npw_k*my_nspinor
       do iblockbd=1,nbdmix
!        scf_history%alpha contains dtset%wfmix in the simple mixing case.
         cg(:,icg+1+ibd:icg+inc+ibd)=scf_history%cg(:,1+icg+ibd:icg+ibd+inc,indh)&
&          +scf_history%alpha*(cg(:,icg+1+ibd:ibd+icg+inc)-scf_history%cg(:,1+icg+ibd:icg+ibd+inc,indh))
         if(usepaw==1) then
           alpha=one-scf_history%alpha
           beta=scf_history%alpha
           call pawcprj_axpby(alpha,beta,cprj_kh(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
         endif
         ibd=ibd+inc
       end do ! end loop on iblockbd

!      Back to usual orthonormalization 
       call cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf_k,mcg,my_nspinor*nband_k,dtset%mkmem,&
&        mpi_enreg,dtset%natom,nattyp,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,usepaw)

!      Store the newly extrapolated wavefunctions, orthonormalized, in scf_history
       ibd=0
       inc=npw_k*my_nspinor
       do iblockbd=1,nbdmix
         scf_history%cg(:,1+icg+ibd:icg+ibd+inc,indh)=cg(:,icg+1+ibd:ibd+icg+inc)
         if(usepaw==1) then
           call pawcprj_put(atindx1,cprj_k,scf_history%cprj(:,:,indh),dtset%natom,1,ibg,ikpt,1,isppol,&
&           dtset%mband,dtset%mkmem,dtset%natom,nbdmix,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
         end if
         ibd=ibd+inc
       end do ! end loop on iblockbd

       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(cwavefh)
       ABI_DEALLOCATE(smn)
       if(usepaw==1) then
         call pawcprj_free(cprj_k)
         call pawcprj_free(cprj_kh)
       end if
       ABI_DATATYPE_DEALLOCATE(cprj_k)
       ABI_DATATYPE_DEALLOCATE(cprj_kh)

       ibg=ibg+my_nspinor*nband_k
       icg=icg+my_nspinor*nband_k*npw_k

!      End big k point loop
     end do
!    End loop over spins
   end do

 end if ! istep>=2

!DEBUG
! write(std_out,*)' wf_mixing : exit '
! write(std_out,*)' cg(1:2,1:2)=',cg(1:2,1:2)
! write(std_out,*)' scf_history%cg(1:2,1:2,1)=',scf_history%cg(1:2,1:2,1)
! ABI_DEALLOCATE(cg_ref)
! ABI_DATATYPE_DEALLOCATE(cprj_ref)
!ENDDEBUG

 ABI_DEALLOCATE(dimcprj)

end subroutine wf_mixing
!!***

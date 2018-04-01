!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_bandfft_tabs
!! NAME
!! prep_bandfft_tabs
!!
!! FUNCTION
!! This routine transpose various tabs needed in bandfft parallelization
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FB,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ikpt=index of the k-point
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!
!! SIDE EFFECTS
!!  bandfft_kpt tabs (defined in m_bandfft_kpt module)
!!
!! PARENTS
!!      energy,fock2ACE,forstrnps,vtorho
!!
!! CHILDREN
!!      bandfft_kpt_init2,bandfft_kpt_set_ikpt,mkkpg,timab,xmpi_allgatherv
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine prep_bandfft_tabs(gs_hamk,ikpt,mkmem,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_time,        only : timab
 use m_kg,          only : mkkpg
 use m_bandfft_kpt, only : bandfft_kpt,bandfft_kpt_set_ikpt,bandfft_kpt_init2
 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_bandfft_tabs'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 integer,intent(in) :: ikpt,mkmem
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: dimffnl,ierr,ikpt_this_proc,ipw,lmnmax,matblk,ndatarecv,nkpg,npw_k,ntypat,spaceComm
 logical :: tabs_allocated
 real(dp) :: tsec(2)
 character(len=500)   :: message
 integer, allocatable :: recvcounts(:),rdispls(:)
 integer, allocatable :: recvcountsloc(:),rdisplsloc(:)
 real(dp),allocatable :: ffnl_gather(:,:,:,:),ffnl_little(:,:,:,:),ffnl_little_gather(:,:,:,:)
 real(dp),allocatable :: kinpw_gather(:),kpg_k_gather(:,:)
 real(dp),allocatable :: ph3d_gather(:,:,:),ph3d_little(:,:,:),ph3d_little_gather(:,:,:)

! *********************************************************************

 call timab(575,1,tsec)

 ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
 tabs_allocated=((bandfft_kpt(ikpt_this_proc)%flag1_is_allocated==1).and.&
& (ikpt_this_proc <= mkmem).and.(ikpt_this_proc/=0))

 if (.not.tabs_allocated) then
   message = ' the bandfft tabs are not allocated !'
   MSG_BUG(message)
 end if

 ntypat=gs_hamk%ntypat
 lmnmax=gs_hamk%lmnmax
 matblk=gs_hamk%matblk
 npw_k=gs_hamk%npw_k
 dimffnl=size(gs_hamk%ffnl_k,2)
 nkpg=size(gs_hamk%kpg_k,2)

 ABI_ALLOCATE(rdispls,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(recvcounts,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(rdisplsloc,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(recvcountsloc,(mpi_enreg%nproc_band))

 spaceComm    =mpi_enreg%comm_band
 ndatarecv    =bandfft_kpt(ikpt_this_proc)%ndatarecv
 rdispls(:)   =bandfft_kpt(ikpt_this_proc)%rdispls(:)
 recvcounts(:)=bandfft_kpt(ikpt_this_proc)%recvcounts(:)

!---- Process FFNL
 if (associated(gs_hamk%ffnl_k)) then
   ABI_ALLOCATE(ffnl_gather,(ndatarecv,dimffnl,lmnmax,ntypat))
   ABI_ALLOCATE(ffnl_little,(dimffnl,lmnmax,ntypat,npw_k))
   ABI_ALLOCATE(ffnl_little_gather,(dimffnl,lmnmax,ntypat,ndatarecv))
   do ipw=1,npw_k
     ffnl_little(:,:,:,ipw)=gs_hamk%ffnl_k(ipw,:,:,:)
   end do
   recvcountsloc(:)=recvcounts(:)*dimffnl*lmnmax*ntypat
   rdisplsloc(:)=rdispls(:)*dimffnl*lmnmax*ntypat
   call xmpi_allgatherv(ffnl_little,npw_k*dimffnl*lmnmax*ntypat,ffnl_little_gather,&
&   recvcountsloc(:),rdisplsloc,spaceComm,ierr)
   do ipw=1,ndatarecv
     ffnl_gather(ipw,:,:,:)=ffnl_little_gather(:,:,:,ipw)
   end do
   ABI_DEALLOCATE(ffnl_little)
   ABI_DEALLOCATE(ffnl_little_gather)
 else
   ABI_ALLOCATE(ffnl_gather,(0,0,0,0))
 end if

!---- Process PH3D
 if (associated(gs_hamk%ph3d_k)) then
   ABI_ALLOCATE(ph3d_gather,(2,ndatarecv,matblk))
   ABI_ALLOCATE(ph3d_little,(2,matblk,npw_k))
   ABI_ALLOCATE(ph3d_little_gather,(2,matblk,ndatarecv))
   recvcountsloc(:)=recvcounts(:)*2*matblk
   rdisplsloc(:)=rdispls(:)*2*matblk
   do ipw=1,npw_k
     ph3d_little(:,:,ipw)=gs_hamk%ph3d_k(:,ipw,:)
   end do
   call xmpi_allgatherv(ph3d_little,npw_k*2*matblk,ph3d_little_gather,&
&   recvcountsloc(:),rdisplsloc,spaceComm,ierr)
   do ipw=1,ndatarecv
     ph3d_gather(:,ipw,:)=ph3d_little_gather(:,:,ipw)
   end do
   ABI_DEALLOCATE(ph3d_little_gather)
   ABI_DEALLOCATE(ph3d_little)
 else
   ABI_ALLOCATE(ph3d_gather,(0,0,0))
 end if

!---- Process KPG_K
 if (associated(gs_hamk%kpg_k)) then
   ABI_ALLOCATE(kpg_k_gather,(ndatarecv,nkpg))
   if (nkpg>0) then
     call mkkpg(bandfft_kpt(ikpt_this_proc)%kg_k_gather,kpg_k_gather,gs_hamk%kpt_k,nkpg,ndatarecv)
   end if
 else
   ABI_ALLOCATE(kpg_k_gather,(0,0))
 end if

!---- Process KINPW
 if (associated(gs_hamk%kinpw_k)) then
   ABI_ALLOCATE(kinpw_gather,(ndatarecv))
   recvcountsloc(:)=recvcounts(:)
   rdisplsloc(:)=rdispls(:)
   call xmpi_allgatherv(gs_hamk%kinpw_k,npw_k,kinpw_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)
 else
   ABI_ALLOCATE(kinpw_gather,(0))
 end if

 ABI_DEALLOCATE(recvcounts)
 ABI_DEALLOCATE(rdispls)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)

 call bandfft_kpt_init2(bandfft_kpt,dimffnl,ffnl_gather,ikpt_this_proc,kinpw_gather,kpg_k_gather,&
& lmnmax,matblk,mkmem,ndatarecv,nkpg,ntypat,ph3d_gather)

 ABI_DEALLOCATE(ffnl_gather)
 ABI_DEALLOCATE(ph3d_gather)
 ABI_DEALLOCATE(kinpw_gather)
 ABI_DEALLOCATE(kpg_k_gather)

!---- Store current kpt index
 call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

 call timab(575,2,tsec)

end subroutine prep_bandfft_tabs
!!***

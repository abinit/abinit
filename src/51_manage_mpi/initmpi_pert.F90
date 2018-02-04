!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_pert
!! NAME
!!  initmpi_pert
!!
!! FUNCTION
!!  Creates group for Parallelization over Perturbations.
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2018 ABINIT group (FJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      mpi_setup
!!
!! CHILDREN
!!      get_npert_rbz,xmpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine initmpi_pert(dtset,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_pert'
 use interfaces_51_manage_mpi, except_this_one => initmpi_pert
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
!scalars
 integer:: iprocmin,irank,npert,nproc_per_cell,nrank,numproc
 integer,allocatable :: ranks(:)
 character(len=500) :: msg
!arrays
 integer,pointer :: nkpt_rbz(:)
 real(dp),pointer :: nband_rbz(:,:)

! ***********************************************************************

 if (mpi_enreg%me_pert<0) then
   msg='Error in MPI distribution! Change your proc(s) distribution or use autoparal>0.'
   MSG_ERROR(msg)
 end if

 call get_npert_rbz(dtset,nband_rbz,nkpt_rbz,npert)

 if (dtset%nppert>=1) then
   if (mpi_enreg%comm_cell/=mpi_enreg%comm_world) then
     call xmpi_comm_free(mpi_enreg%comm_cell)
   end if
   mpi_enreg%comm_cell=mpi_enreg%comm_world

   ! These values will be properly set in set_pert_comm
   mpi_enreg%me_cell=mpi_enreg%me
   mpi_enreg%nproc_cell=mpi_enreg%nproc

   if (mpi_enreg%me>=0) then
     nproc_per_cell=mpi_enreg%nproc/dtset%nppert
     ABI_ALLOCATE(ranks,(dtset%nppert))
     iprocmin=mod(mpi_enreg%me,nproc_per_cell)
     ranks=(/((iprocmin+(irank-1)*nproc_per_cell),irank=1,dtset%nppert)/)
     mpi_enreg%comm_pert=xmpi_subcomm(mpi_enreg%comm_world,dtset%nppert,ranks)
     ABI_DEALLOCATE(ranks)
     mpi_enreg%me_pert=xmpi_comm_rank(mpi_enreg%comm_pert)
     mpi_enreg%nproc_pert=dtset%nppert
     if (iprocmin==0.and.mpi_enreg%me_pert==0.and.mpi_enreg%me/=0) then
       MSG_BUG('Error on me_pert!')
     end if
!    Define mpi_enreg%distrb_pert
     ABI_ALLOCATE(mpi_enreg%distrb_pert,(npert))
     nrank=0
     do irank=1,npert
       nrank=nrank+1
       mpi_enreg%distrb_pert(irank)=mod(nrank,dtset%nppert)-1
       if (mpi_enreg%distrb_pert(irank)==-1) mpi_enreg%distrb_pert(irank)=dtset%nppert-1
     end do
     ! Make sure that subrank 0 is working on the last perturbation
     ! Swap the ranks if necessary
     numproc=mpi_enreg%distrb_pert(npert)
     if(numproc/=0) then
       do irank=1,npert
         if (mpi_enreg%distrb_pert(irank)==numproc) mpi_enreg%distrb_pert(irank)=-2
         if (mpi_enreg%distrb_pert(irank)==0) mpi_enreg%distrb_pert(irank)=-3
       end do
       do irank=1,npert
         if (mpi_enreg%distrb_pert(irank)==-2) mpi_enreg%distrb_pert(irank)=0
         if (mpi_enreg%distrb_pert(irank)==-3) mpi_enreg%distrb_pert(irank)=numproc
       end do
     end if
!    Communicator over one cell
     ABI_ALLOCATE(ranks,(nproc_per_cell))
     iprocmin=(mpi_enreg%me/nproc_per_cell)*nproc_per_cell
     ranks=(/((iprocmin+irank-1),irank=1,nproc_per_cell)/)
     mpi_enreg%comm_cell_pert=xmpi_subcomm(mpi_enreg%comm_world,nproc_per_cell,ranks)
     ABI_DEALLOCATE(ranks)
   end if

 else  !nppert<=1
   mpi_enreg%nproc_pert=1
   mpi_enreg%comm_pert=xmpi_comm_self
   mpi_enreg%me_pert=0
   ABI_ALLOCATE(mpi_enreg%distrb_pert,(npert))
   mpi_enreg%distrb_pert(:)=0
 end if

 ABI_DEALLOCATE(nband_rbz)
 ABI_DEALLOCATE(nkpt_rbz)

end subroutine initmpi_pert
!!***

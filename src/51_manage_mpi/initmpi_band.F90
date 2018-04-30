!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_band
!! NAME
!!  initmpi_band
!!
!! FUNCTION
!!  Initializes the mpi informations for band parallelism (paralbd=1).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_enreg= informations about MPI parallelization
!!  nband(nkpt*nsppol)= number of bands per k point, for each spin
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for polarized
!!
!! OUTPUT
!!  mpi_enreg=information about MPI parallelization
!!  mpi_enreg%comm_band=communicator of BAND set
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_band(mpi_enreg,nband,nkpt,nsppol)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_band'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: ii,ikpt,iproc_min,iproc_max,irank,isppol
 integer :: me,nband_k,nproc,nbsteps,nrank,nstates,spacecomm
 character(len=500) :: msg
!arrays
 integer,allocatable :: ranks(:)

! ***********************************************************************

 mpi_enreg%comm_band=xmpi_comm_self

 if (mpi_enreg%paralbd==1.and.xmpi_paral==1) then

!  Comm_kpt is supposed to treat spins, k-points and bands
   spacecomm=mpi_enreg%comm_kpt
   nproc=mpi_enreg%nproc_kpt
   me=mpi_enreg%me_kpt

   nstates=sum(nband(1:nkpt*nsppol))
   nbsteps=nstates/nproc
   if (mod(nstates,nproc)/=0) nbsteps=nbsteps+1

   if (nbsteps<maxval(nband(1:nkpt*nsppol))) then

     nrank=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         ii=ikpt+(isppol-1)*nkpt
         nband_k=nband(ii)
         if (nbsteps<nband_k) then
           iproc_min=minval(mpi_enreg%proc_distrb(ikpt,:,isppol))
           iproc_max=maxval(mpi_enreg%proc_distrb(ikpt,:,isppol))
           if ((me>=iproc_min).and.(me<=iproc_max)) then
             nrank=iproc_max-iproc_min+1
             if (.not.allocated(ranks)) then
               ABI_ALLOCATE(ranks,(nrank))
               if (nrank>0) ranks=(/((iproc_min+irank-1),irank=1,nrank)/)
             else if (nrank/=size(ranks)) then
               msg='Number of bands per proc should be the same for all k-points!'
               MSG_BUG(msg)
             end if
           end if
         end if
       end do
     end do
     if (.not.allocated(ranks)) then
       ABI_ALLOCATE(ranks,(0))
     end if

     mpi_enreg%comm_band=xmpi_subcomm(spacecomm,nrank,ranks)

     ABI_DEALLOCATE(ranks)
   end if
 end if

end subroutine initmpi_band
!!***

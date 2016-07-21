!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_atom
!! NAME
!!  initmpi_atom
!!
!! FUNCTION
!!  Initializes the mpi informations for parallelism over atoms (PAW).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg= informations about MPI parallelization
!!
!! OUTPUT
!!  mpi_enreg= informations about MPI parallelization
!!    comm_atom                 =communicator over atoms
!!    nproc_atom                =size of the communicator over atoms
!!    my_natom                  =number of atoms treated by current proc
!!    my_atmtab(mpi_enreg%natom)=indexes of the atoms treated by current processor
!!
!! PARENTS
!!      m_paral_pert,mpi_setup
!!
!! CHILDREN
!!      get_my_atmtab,get_my_natom
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_atom(dtset,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_paral_atom, only : get_my_natom, get_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_atom'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays

!Local variables-------------------------------
!scalars
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
 integer :: iatom

! ***********************************************************************

 DBG_ENTER("COLL")

 mpi_enreg%nproc_atom=1
 mpi_enreg%comm_atom=xmpi_comm_self
 mpi_enreg%my_natom=dtset%natom
 if (associated(mpi_enreg%my_atmtab))then
   ABI_DEALLOCATE(mpi_enreg%my_atmtab)
 end if
 nullify(mpi_enreg%my_atmtab)

 if (xmpi_paral==0) then 
   mpi_enreg%nproc_atom=0
   ABI_ALLOCATE(mpi_enreg%my_atmtab,(0))
   return
 end if


!Check compatibility
 if (dtset%paral_atom>0) then
   msg=''
   if (dtset%usepaw==0)  msg= 'Parallelisation over atoms not compatible with usepaw=0 !'
   if (dtset%usedmft==1) msg=' Parallelisation over atoms not compatible with usedmft=1 !'
   if (dtset%usewvl==1)  msg= 'Parallelisation over atoms not compatible with usewvl=1 !'
   if (dtset%prtden>1.and.dtset%paral_kgb<=0) &
&   msg= 'Parallelisation over atoms not compatible with prtden>1 (PAW AE densities) !'
   if (dtset%optdriver/=RUNL_GSTATE.and.dtset%optdriver/=RUNL_RESPFN) &
&   msg=' Parallelisation over atoms only compatible with GS or RF !'
   if (dtset%macro_uj/=0)msg=' Parallelisation over atoms not compatible with macro_uj!=0 !'
   if (msg/='') then
     MSG_ERROR(msg)
   end if
 end if

 if (mpi_enreg%comm_atom==xmpi_comm_null) then
   mpi_enreg%nproc_atom=0;mpi_enreg%my_natom=0
   ABI_ALLOCATE(mpi_enreg%my_atmtab,(0))
   return
 end if

 if (dtset%paral_atom>0) then

!  Build correct atom communicator
   if (dtset%optdriver==RUNL_GSTATE.and.dtset%paral_kgb==1) then
     mpi_enreg%comm_atom=mpi_enreg%comm_kptband
   else
     mpi_enreg%comm_atom=mpi_enreg%comm_cell
   end if

!  Get number of processors sharing the atomic data distribution
   mpi_enreg%nproc_atom=xmpi_comm_size(mpi_enreg%comm_atom)

!  Get local number of atoms
   call get_my_natom(mpi_enreg%comm_atom,mpi_enreg%my_natom,dtset%natom)
   paral_atom=(mpi_enreg%my_natom/=dtset%natom)

!  Build atom table
   if (mpi_enreg%my_natom>0.and.paral_atom) then
     my_atmtab_allocated=.false.
     call get_my_atmtab(mpi_enreg%comm_atom,mpi_enreg%my_atmtab,my_atmtab_allocated, &
&     paral_atom,dtset%natom)
   else if (.not.paral_atom) then
     ABI_ALLOCATE(mpi_enreg%my_atmtab,(dtset%natom))
     mpi_enreg%my_atmtab(1:dtset%natom)=(/(iatom, iatom=1,dtset%natom)/)
   else if (mpi_enreg%my_natom==0) then
     ABI_ALLOCATE(mpi_enreg%my_atmtab,(0))
   end if

 end if

 DBG_EXIT("COLL")

end subroutine initmpi_atom
!!***

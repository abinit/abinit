!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_supercell
!! NAME
!!  initmpi_supercell
!!
!! FUNCTION
!!  Initializes the mpi informations for parallelism over supercell.
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
!!  mpi_enreg= informations about MPI parallelization
!!
!! OUTPUT
!! This is for the parallelisation over the supercell
!!
!!  mpi_enreg= informations about MPI parallelization
!!   comm_supercell =  Communicator over all processors treating the same cell
!!   me_supercell =  Index of my processor in the comm. over one cell
!!   my_ncell     =  Number of cell treated by current proc
!!   my_cells(:)  = Number of the cells in the supercell treat by this CPU
!!   my_index_cells(:,:) = indexes of the cells in the supercell treat by this CPU
!!
!! PARENTS
!!     
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_supercell(cell_number,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_supercell'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: cell_number(3)
!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,icell,ii
 integer :: ncell_alone,ncell
 character(len=500) :: msg
!array
! ***********************************************************************

 DBG_ENTER("COLL")

!set the number of cell
 ncell = product(cell_number(:))

! Do some checks
  if (any(cell_number <= 0).or.ncell<=0) then
    write(msg,'(a,a)')' No supercell found for setting'
    MSG_ERROR(msg)
  end if

!Determine the number of cell for each CPU
 ncell_alone = mod(ncell,mpi_enreg%nproc)
 mpi_enreg%my_ncell = aint(real(ncell,sp)/mpi_enreg%nproc)  
 if(mpi_enreg%me >= (mpi_enreg%nproc-ncell_alone)) then
   mpi_enreg%my_ncell = mpi_enreg%my_ncell  + 1
 end if

!Allocation of array
 ABI_ALLOCATE(mpi_enreg%my_cells,(mpi_enreg%my_ncell))
 ABI_ALLOCATE(mpi_enreg%my_index_cells,(mpi_enreg%my_ncell,3))
 mpi_enreg%my_cells = zero
 mpi_enreg%my_index_cells = zero

 do icell = 1,mpi_enreg%my_ncell
   if(mpi_enreg%me >= (mpi_enreg%nproc-ncell_alone))then
     mpi_enreg%my_cells(icell)=(aint(real(ncell,sp)/mpi_enreg%nproc))*(mpi_enreg%me)+&
&                              (mpi_enreg%me - (mpi_enreg%nproc-ncell_alone)) + icell
   else
     mpi_enreg%my_cells(icell)=(mpi_enreg%my_ncell)*(mpi_enreg%me)  + icell
   end if
 end do
 
 icell = 0
 ii = 0
 do i1 = 1,cell_number(1)
   do i2 = 1,cell_number(2)
     do i3 = 1,cell_number(3)
       ii = ii +1
       if(any(mpi_enreg%my_cells==ii))then
         icell=icell+1
         mpi_enreg%my_index_cells(icell,1) = i1; 
         mpi_enreg%my_index_cells(icell,2) = i2; 
         mpi_enreg%my_index_cells(icell,3) = i3;
       end if
     end do
   end do
 end do

 DBG_EXIT("COLL")

end subroutine initmpi_supercell
!!***
